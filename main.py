from typing import List
import numpy as np
from copy import deepcopy
import logging 
from pyscf import gto, scf, cc
from pyscf.gto import Mole
# from qiskit_algorithms.optimizers import SPSA ,CG, SLSQP, L_BFGS_B, COBYLA
from qiskit_algorithms.optimizers import CG, SLSQP, L_BFGS_B, COBYLA, SPSA, GradientDescent

from scipy.optimize import BFGS
from Ternary_Tree.optimizer.soap import SOAP
from ExcitationSolve.excitationsolve.excitation_solve_qiskit import ExcitationSolveQiskit
import json
from copy import deepcopy
from Ternary_Tree.ucc.abstractucc import Molecule
from utils import *
from my_utils import Timer
import multiprocessing as mp

import sys
from pyscf import lib


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
        
optimizers = [(SLSQP(maxiter=50, ftol=0), 'SLSQP')]
# optimizers = [(GradientDescent(maxiter=50, ftol=0), 'GD')]
# optimizers = [(COBYLA(maxiter=100, tol=0), 'COBYLA')]

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
    
optimizers = [(SPSA(maxiter=n, blocking=True, allowed_increase=0.00001), 'SPSA')]
optimizers = [(SOAP(maxfev=n), 'SOAP')]

optimizers = [(L_BFGS_B(maxiter=450, ftol=0.0000001), 'L_BFGS_B')]
#optimizers = [(ExcitationSolveQiskit(maxiter=100), 'ES')]
# optimizers = [(BFGS(maxiter=50, ftol=0.000001), 'L_BFGS_B')]


R = np.linspace(1,1.6,12)
class vqeData:
    def __init__(self, 
                 file_name_to_read: str,
                 file_name_to_write: str,
                  molecule: Molecule,
                  optimizer: any,
                  reps: int=1,
                  noise_type: str="",
                  probs: np.ndarray=np.flip(np.geomspace(0.00002, (0.001), 1)),
                  device: str="CPU",
                  ):
        self.file_name_to_write = file_name_to_write
        self.file_name_to_ideal = file_name_to_read 
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

@Timer.attach_timer("thread_timer")
def to_thread(namet, vqe_data: vqeData, r:float=0, index: int=0, init_point=None, reps=1):
    # probs = vqe_data.probs
    # if vqe_data.noise_type != "":
    try:
        ref_en = -4.29884
        with open(get_file_name(vqe_data.file_name_to_ideal, "", namet), "r") as file:
            js = json.load(file)[index]
            init_point = js["param"]
            ref_en = js["ref_ener"]
            ref_en = -4.29884
            # init_point = None
            # r = js["dist"]
            if isinstance(init_point, float):
                init_point = np.array([init_point])
        if vqe_data.noise_type == "":
            probs = [0]
        else:
            probs = vqe_data.probs
            
        # pass
    except FileNotFoundError:
        pass
        print("RUNNING NOISY SIM WITHOUT INIT POINT")
    data = []
    name, circ, op = vqe_data.circ_prov.get_circ(namet)
    # print(circ)
    for prob in probs:
        print(prob)
        print(vqe_data.noise_type)
        circs = CircSim(circ, op, prob, vqe_data.noise_type, init_point)
        energy, parameters = circs.run_qiskit_vqe(vqe_data.optimizer[0], vqe_data.device, reps=reps, random_init=True)
        data.append({
            "name": name, 
            "ref_ener": vqe_data.ref_value, 
            "energy": energy, 
            "param": parameters, 
            "optimizer": vqe_data.optimizer[1],
            "gate_count": circ.count_ops(),
            "dist": r,
            "prob": prob,
            "noise": vqe_data.noise_type
        })
    return data
        
def run_vqe(name: str, vqe_data: vqeData, data: dict, r:float=0, index:int=0):
    init_point = None
    reps = 1
    try:
        init_point = data[name + ""][-1]["param"]
        # reps = 1
    except KeyError:
        print("[thyz]")
        pass
    # noise = vqe_data.noise_type 
    # vqe_data.noise_type = ""
    # result = to_thread(name, vqe_data, r, index, init_point, reps)
    # init_point = result[0]["param"]
    # vqe_data.noise_type = noise
    result = to_thread(name, vqe_data, r, index, init_point, reps)
    if name + vqe_data.noise_type not in data:
        data.setdefault(name + vqe_data.noise_type, result)
    else:
        data[name + vqe_data.noise_type] = data[name + vqe_data.noise_type] + result
    print("thread_timer" + ": ", Timer.timers["thread_timer"])
    return result

def run_ccsd(mol: Molecule):
    import pyscf.cc as cc
    mol: Mole = mol.to_pyscf_mol()
    mf = scf.RHF(mol).run()
    # Note that the line following these comments could be replaced by
    # mycc = cc.CCSD(mf)
    # mycc.kernel()
    mycc = cc.CCSD(mf).run()
    print('NUC', mf.energy_nuc())
    print('CCSD total energy', mycc.e_tot - mf.energy_nuc())
    et = mycc.ccsd_t()
    print('CCSD(T) total energy', mycc.e_tot + et)
    # with open(get_file_name(vqe_data.file_name_to_write, vqe_data.noise_type, name), 'w') as file:
        # pass
        # json.dump(mycc, file, indent=4)
    return mycc.e_tot - mf.energy_nuc()

if __name__ == "__main__":
    # logging.basicConfig(level=logging.DEBUG)  # Enable logging

    circ_names = Circuits.get_circs_names()[-3:]
    manager = mp.Manager()
    data = manager.dict()
    # for noise in ["D", "X", "Z"]:
    procs = []
    vqe_data=vqeData(
                # "datah2h2_bfgs/DH4",
                # "data04/H2_4",
                # "datah2_4/H2_4",
                # "data04/H2_4",
                # "datah2h2copy/DH4",
                "datah2h2ExcSolProb/DH4",
                "datah2h2noise_level/DH4",
                DH4,
                optimizers[0],
                reps=2,
                probs=1 - np.flip(np.geomspace((0.001), 0.000005, 15)),
                noise_type="",
                device="CPU",
            )
    # for noise in ["D", "X", "Y", "Z", "sc"]:
    R = [1.23]
    procs = []
    for noise in ["X", "Z"]:
        vqe_data.noise_type = noise        
        for index, r in enumerate(R):
            # vqe_data.molecule.set_distance(r)
            # r = 0.7349
            # index = 0
            # vqe_data.update_prov()
            # print("hello")
            for name in circ_names:
            # vqe_data2.molecule.set_distance(r)
            # vqe_data2.update_prov()
                # procs.append(mp.Process(target=run_vqe, args=(name, vqe_data, data)))
                procs.append(mp.Process(target=run_vqe, args=(name, vqe_data, data, r, index)))
                
    for proc in procs:
        proc.start()
    for proc in procs:
        proc.join()
    for proc in procs:
        proc.close()
            
    # for noise in ["D", "X", "Y", "Z"]:
        for name in circ_names:
            with open(get_file_name(vqe_data.file_name_to_write, noise, name), 'w') as file:
                json.dump(data[name + noise], file, indent=4)
                
    