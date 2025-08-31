from typing import List
import numpy as np
from copy import deepcopy
import json

from qiskit.quantum_info import DensityMatrix
from qiskit_algorithms.optimizers import CG, SLSQP, L_BFGS_B, COBYLA, SPSA
from qiskit_nature.second_q.operators import FermionicOp

from copy import deepcopy
from Ternary_Tree.ucc.abstractucc import Molecule
from Ternary_Tree.qiskit_interface.circuit_provider import *
from my_utils import Timer
import multiprocessing as mp

import logging


logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
logger.addHandler(handler)
logger.setLevel(logging.INFO)


class H2_H2(Molecule):
    def __init__(self, R: float):
        super().__init__(geometry=f'H 0 0 0; H 0 0 1.23; H {R} 0 0; H {R} 0 1.23', 
                         num_electrons=(2,2), 
                         active_orbitals=[0,1,2,3], 
                         basis='sto-3g')
        
    def set_distance(self, R: float):
        self.geometry = f'H 0 0 0; H 0 0 1.23; H {R} 0 0; H {R} 0 1.23'
        
DH4 = H2_H2(1.23)
H2_4 = Molecule(geometry='H 0 0 0; H 0 0 0.7349', num_electrons=(1,1), active_orbitals=[0,1], basis='sto-3g')
H2_8 = Molecule(geometry='H 0 0 0; H 0 0 0.7349', num_electrons=(1,1), active_orbitals=[0,1,2,3], basis='6-31g')
LiH_8 = Molecule(geometry='H 0 0 0; Li 0 0 1.5459', num_electrons=(2,2), active_orbitals=[0,1,2,5], basis='sto-3g')
LiH_10 = Molecule(geometry='H 0 0 0; Li 0 0 1.5459', num_electrons=(2,2), active_orbitals=[0,1,2,3,5], basis='sto-3g')
        

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
    
optimizers = [(L_BFGS_B(maxiter=150, ftol=0.00000001), 'L_BFGS_B')]


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

def eval_additional_observable(circuit: QuantumCircuit, parameters, estimator: Estimator,  sim: AerSimulator, mapper ):
    num_spin_orbitals = circuit.num_qubits
    def num_observable(num_spin_orbitals, mapper):
        ferm_observable_N = FermionicOp(
            {f"+_{i} -_{i}": 1 for i in range(num_spin_orbitals)},
            num_spin_orbitals=num_spin_orbitals
        )
        return mapper.map(ferm_observable_N)
    
    circuit = circuit.assign_parameters(parameters)
    circuit.save_density_matrix() 
    n_obs = num_observable(num_spin_orbitals, mapper)
    res = estimator.run([circuit]*2, observables=[n_obs, n_obs @ n_obs]).result()
    result = sim.run(circuit).result()
    rho = result.data(0)['density_matrix']
    purity = (rho.data @ rho.data).trace().real
    return list(res.values) + [res.values[1] - res.values[0]**2] + [purity]


@Timer.attach_timer("thread_timer")
def to_thread(namet, vqe_data: vqeData, r:float=0):
    init_point = None
    probs = vqe_data.probs
    data = []
    name, circ, op_mapper = vqe_data.circ_prov.get_circ(namet)
    for index, prob in enumerate(probs):
        circs = CircSim(circ, op_mapper[0], prob, vqe_data.noise_type, init_point)
        energy, parameters, est, sim, cb = circs.run_qiskit_vqe(vqe_data.optimizer[0], vqe_data.device, reps=1)
        additional_res = eval_additional_observable(circs.circ, parameters, est, sim, op_mapper[1])
        logger.info(f"{energy:.5f} hf: {name=}, {noise=}, {index=}")
        data.append({
            "name": name, 
            "ref_ener": vqe_data.ref_value, 
            "energy": energy, 
            "energy_array": cb.energy_array,
            "param": parameters, 
            "addition_res:": additional_res,
            "optimizer": vqe_data.optimizer[1],
            "gate_count": circs.circ.count_ops(),
            "dist": r,
            "prob": prob,
            "noise": vqe_data.noise_type
        })
    return data
        
def run_vqe(name: str, vqe_data: vqeData, data: dict, r:float=0):
    result = to_thread(name, vqe_data, r)
    if name + vqe_data.noise_type not in data:
        data.setdefault(name + vqe_data.noise_type, result)
    else:
        data[name + vqe_data.noise_type] = data[name + vqe_data.noise_type] + result
    logger.info(f"thread_timer: {Timer.timers["thread_timer"]:.4f} sec")
    return result

if __name__ == "__main__":
    mult = [0.000005, 0.00001, 0.0005, 0.001, 1][0:-1]
    # for noise in ["D", "X", "Y", "Z"]:
    for noise in ["sc", "ion", "D", "X","Y","Z"][1:2]:
    # for noise in ["D"]:
    # for noise in [""]:
        if noise in {"sc", "ion"}:
            probs = mult[:4]
        else:
            probs = 1 - np.flip(np.geomspace(0.000001, (0.0002), 5))
        vqe_data=vqeData(
                "data_10/LiH_10",
                LiH_10,
                optimizers[0],
                reps=1,
                probs=probs,
                noise_type=noise,
                device="GPU",
            )

        circ_names = Circuits.get_circs_names()[0:2] + Circuits.get_circs_names()[6:7]
        # + Circuits.get_circs_names()[6:7]
        print(circ_names) 
        # manager = mp.Manager()
        # data = manager.dict()
        data = {}
        r = 1.23
        # vqe_data.molecule.set_distance(r)
        # vqe_data.update_prov()
        procs = []
        for name in circ_names:
            run_vqe(name, vqe_data, data, r)
            # procs.append(mp.Process(target=run_vqe, args=(name, vqe_data, data, r) ))
        # for proc in procs:
        #     proc.start()
        # for proc in procs:
        #     proc.join()
        # for proc in procs:
        #     proc.close()
            
            # for name in circ_names:
            with open(get_file_name(vqe_data.file_name, vqe_data.noise_type, name), 'w') as file:
                json.dump(data[name + vqe_data.noise_type], file, indent=4)
