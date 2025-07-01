from typing import List
import numpy as np
from copy import deepcopy
import json

from qiskit_algorithms.optimizers import CG, SLSQP, L_BFGS_B, COBYLA, SPSA
from qiskit_nature.second_q.operators import FermionicOp

from copy import deepcopy
from Ternary_Tree.ucc.abstractucc import Molecule
from Ternary_Tree.qiskit_interface.circuit_provider import *
from my_utils import Timer
import multiprocessing as mp
from Ternary_Tree import logger
import sys
# from pyscf import lib

import logging

# logger.setLevel(logging.INFO)

# logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
logger.addHandler(handler)
logger.setLevel(logging.INFO)
# logger.info("hello")
# logger.warning("warn hello")
# log = lib.logger.Logger(sys.stdout, 4)
# log.info('info level')
# log.verbose = 3
# log.info('info level')
# log.note('note level')

class H2_H2(Molecule):
    def __init__(self, R: float):
        super().__init__(geometry=f'H 0 0 0; H 0 0 1.23; H {R} 0 0; H {R} 0 1.23', 
                         num_electrons=(2,2), 
                         active_orbitals=[0,1,2,3], 
                         basis='sto-3g')
        
    def set_distance(self, R: float):
        self.geometry = f'H 0 0 0; H 0 0 1.23; H {R} 0 0; H {R} 0 1.23'
        
H2_4 = Molecule(geometry='H 0 0 0; H 0 0 0.7349', num_electrons=(1,1), active_orbitals=[0,1], basis='sto-3g')
H2_8 = Molecule(geometry='H 0 0 0; H 0 0 0.7349', num_electrons=(1,1), active_orbitals=[0,1,2,3], basis='6-31g')
LiH_8 = Molecule(geometry='H 0 0 0; Li 0 0 1.5459', num_electrons=(2,2), active_orbitals=[0,1,2,5], basis='sto-3g')
DH4 = H2_H2(1.23)
        

n = 100
def learning_rate(n=n, c=0.5):
    k = 0
    while k < n:
        k += 1
        yield c/k**0.3

def perturabation(n=n, a=0.01):
    k = 0
    while k < n:
        k += 1
        yield a/k**0.3
    
optimizers = [(L_BFGS_B(maxiter=80, ftol=0.0000001), 'L_BFGS_B')]


class vqeData:
    def __init__(self, 
                 file_name_to_write: str,
                  molecule: Molecule,
                  optimizer: any,
                  reps: int=1,
                  noise_type: str="",
                  probs: np.ndarray=1 - np.flip(np.geomspace(0.00002, (0.001), 15)),
                  device: str="CPU",
                  ):
        self.file_name = file_name_to_write 
        self.molecule = molecule 
        self.optimizer = optimizer 
        self.reps = reps
        self.noise_type = noise_type 
        self.probs = probs 
        self.circ_prov = CircuitProvider(reps=self.reps, molecule=self.molecule)
        self.ref_value = numpy_energy(self.circ_prov.fermionic_op, self.circ_prov.uccgsd)
        self.data = []
        self.device = device
    
    def update_prov(self):
        self.circ_prov.update_molecule(self.molecule)
        self.ref_value = numpy_energy(self.circ_prov.fermionic_op, self.circ_prov.uccgsd)

def eval_additional_observable(circuit: QuantumCircuit, parameters, estimator,  mapper ):
    num_spin_orbitals = circuit.num_qubits
    def num_observable(num_spin_orbitals, mapper):
        ferm_observable_N = FermionicOp(
            {f"+_{i} -_{i}": 1 for i in range(num_spin_orbitals)},
            num_spin_orbitals=num_spin_orbitals
        )
        return mapper.map(ferm_observable_N)
    circuit = circuit.assign_parameters(parameters)
    n_obs = num_observable(num_spin_orbitals, mapper)
    res = estimator.run([circuit]*2, observables=[n_obs, n_obs @ n_obs]).result()
    return list(res.values) + [res.values[1] - res.values[0]**2]


@Timer.attach_timer("thread_timer")
def to_thread(namet, vqe_data: vqeData, r:float=0, is_rust=False):
    init_point = None
    probs = vqe_data.probs
    data = []
    name, circ, op_mapper = vqe_data.circ_prov.get_circ(namet)
    mapper=JordanWignerMapper()
    if namet == Circuits.bk():
        mapper = BravyiKitaevMapper()
    if is_rust:
        name = name + "_lex"
    for index, prob in enumerate(probs):
        circs = CircSim(circ, op_mapper[0], prob, vqe_data.noise_type, init_point)
        # energy, parameters = circs.run_adapt_vqe(vqe_data.optimizer[0], 
                                                #  vqe_data.device, 
                                                #  reps=1, 
                                                #  is_rust=is_rust, 
                                                #  cp=vqe_data.circ_prov, 
                                                #  mapper=mapper)

        energy, parameters, est = circs.run_qiskit_vqe(vqe_data.optimizer[0], vqe_data.device, reps=1)
        additional_res = eval_additional_observable(circs.circ, parameters, est, op_mapper[1])
        logger.info(f"{energy:.5f} hf: {name=}, {noise=}, {index=}")
        data.append({
            "name": name, 
            "ref_ener": vqe_data.ref_value, 
            "energy": energy, 
            "param": parameters, 
            "addition_res:": additional_res,
            "optimizer": vqe_data.optimizer[1],
            "gate_count": circs.circ.count_ops(),
            "dist": r,
            "prob": prob,
            "noise": vqe_data.noise_type
        })
    return data
        
def run_vqe(name: str, vqe_data: vqeData, data: dict, r:float=0, is_rust=False):
    result = to_thread(name, vqe_data, r, is_rust)
    if is_rust:
        name = name + "_lex"
    if name + vqe_data.noise_type not in data:
        data.setdefault(name + vqe_data.noise_type, result)
    else:
        data[name + vqe_data.noise_type] = data[name + vqe_data.noise_type] + result
    logger.info(f"thread_timer: {Timer.timers["thread_timer"]:.4f} sec")
    return result

if __name__ == "__main__":
    mult = [0.01, 0.1, 0.5,   1, 2, 3, 5, 10]
    for noise in ["D", "X", "Y", "Z"]:
    # for noise in ["", "D", "X","Y","Z"]:
    # for noise in ["X", ]:
    
        vqe_data=vqeData(
                "data/LiH_8",
                LiH_8,
                optimizers[0],
                reps=1,
                # probs=mult,
                probs=1 - np.flip(np.geomspace(0.00001, (0.006), 10)),
                noise_type=noise,
                device="CPU",
            )

        # circ_names = Circuits.get_circs_names()[:]
        # circ_names = Circuits.get_circs_names()[:2] + Circuits.get_circs_names()[4:] +  [Circuits.swap_2xn_alt()]
        circ_names = Circuits.get_circs_names()[:-2]
        # namet = Circuits.swap_2xn_alt()
        # ## circ_names = circ_names[:2] + circ_names[-3:] 
        # ## circ_names = circ_names[2:3]
        # name, circ, op = vqe_data.circ_prov.get_circ(namet)
        # # name, circ, op = vqe_data.circ_prov.get_circ(circ_names[-2])
        # sn = to_excitaions(circ, circ.excitation_pos, circ.excitation_pos)
        # # print(sn)
        # # print(circ.decompose(reps=2))
        # # print(transpile(sn, basis_gates=["rzz", "cz", "rx", "rz"], optimization_level=3))
        # circs = CircSim(sn, op, 0.99999, vqe_data.noise_type)
        # energy, parameters = circs.run_adapt_vqe(vqe_data.optimizer[0], 
        #                                         vqe_data.device, 
        #                                         reps=1, 
        #                                         )
        # circs.run_qiskit_vqe(vqe_data.optimizer[0], vqe_data.device,)
        
        manager = mp.Manager()
        data = manager.dict()
        # for r in R:
        r = 1.23
        # vqe_data.molecule.set_distance(r)
        # vqe_data.update_prov()
        procs = []
        for name in circ_names:
            # if (name == "bk_lex" or name == "jw_lex"):
                # procs.append(mp.Process(target=run_vqe, args=(name[0:2], vqe_data, data, r, True) ))
            # else:
            procs.append(mp.Process(target=run_vqe, args=(name, vqe_data, data, r) ))
        for proc in procs:
            proc.start()
        for proc in procs:
            proc.join()
        for proc in procs:
            proc.close()
            
        for name in circ_names:
            with open(get_file_name(vqe_data.file_name, vqe_data.noise_type, name), 'w') as file:
                json.dump(data[name + vqe_data.noise_type], file, indent=4)
