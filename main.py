from typing import List
import numpy as np

from qiskit_algorithms.optimizers import SPSA ,CG, SLSQP, L_BFGS_B, COBYLA

import json
from copy import deepcopy
from Ternary_Tree.optimizer.soap import SOAP
from Ternary_Tree.ucc.abstractucc import Molecule
from utils import *
from my_utils import Timer
import threading
from multiprocessing import Process
import multiprocessing as mp
from functools import partial
from multiprocessing import Process

H2_4 = Molecule(geometry='H 0 0 0; H 0 0 0.7349', num_electrons=(1,1), active_orbitals=[0,1], basis='sto-3g')
H2_8 = Molecule(geometry='H 0 0 0; H 0 0 0.7349', num_electrons=(1,1), active_orbitals=[0,1,2,3], basis='6-31g')
LiH_8 = Molecule(geometry='H 0 0 0; Li 0 0 1.5459', num_electrons=(2,2), active_orbitals=[0,1,2,5], basis='sto-3g')


# optimizers = [(SLSQP(maxiter=200, ftol=0), 'SLSQP')]
optimizers = [(L_BFGS_B(maxiter=1000, ftol=0), 'L_BFGS_B')]

class vqeData:
    def __init__(self, 
                 file_name_to_read: str,
                  molecule: Molecule,
                  optimizer: any,
                  reps: int=1,
                  noise_type: str="",
                  probs: np.ndarray=np.flip(np.geomspace(0.00002, (0.001), 1)),
                  device: str="CPU",
                  ):
        self.file_name_to_ideal = file_name_to_read 
        self.molecule = molecule 
        self.optimizer = optimizer 
        self.circ_prov = CircuitProvider(reps=reps, molecule=molecule)
        self.noise_type = noise_type 
        self.probs = probs 
        self.ref_value = numpy_energy(self.circ_prov.fermionic_op, self.circ_prov.uccgsd)
        self.data = []
        self.device = device
        

@Timer.attach_timer("thread_timer")
def to_thread(namet, vqe_data: vqeData):
    # vqe_data=deepcopy(vqe_data), 
    name, circ, op = vqe_data.circ_prov.get_circ(namet)
    init_point = None
    probs = [0]
    if vqe_data.noise_type != "":
        probs = vqe_data.probs
        try:
            with open(vqe_data.file_name_to_ideal + "_" + name +  ".json", "r") as file:
                init_point = json.load(file)[0]["param"]
        except FileNotFoundError:
            print("RUNNING NOISY SIM WITHOUT INIT POINT")
    data = []
    for prob in probs:
        circs = CircSim(circ, op, prob, vqe_data.noise_type, init_point)
        energy, parameters = circs.run_qiskit_vqe(vqe_data.optimizer[0], vqe_data.device)
        data.append({
            "name": name, 
            "ref_ener": vqe_data.ref_value, 
            "energy": energy, 
            "param": parameters[-1], 
            "optimizer": vqe_data.optimizer[1],
            "prob": prob,
            "noise": vqe_data.noise_type
        })
    return data
        
def run_vqe(name: str, vqe_data: vqeData):
    data = []
    result = to_thread(name, vqe_data)
    data.extend(result)
    with open(vqe_data.file_name_to_ideal + f"_{vqe_data.noise_type}" + name +".json", 'w') as file:
        json.dump(data, file, indent=4)
    print("thread_timer" + ": ", Timer.timers["thread_timer"])
    return result

if __name__ == "__main__":
    
    vqe_data=vqeData(
            "data03/H2_8",
            H2_8,
            optimizers[0],
            reps=1,
            probs=1 - np.flip(np.geomspace(0.00002, (0.002), 10)),
            noise_type="X",
            device="CPU",
        )
    vqe_data2=vqeData(
            "data03/H2_8",
            H2_8,
            optimizers[0],
            reps=1,
            probs=1 - np.flip(np.geomspace(0.00002, (0.002), 10)),
            noise_type="Z",
            device="GPU",
        )
    circ_names = circ_order()[0:2]
    # circ_names = circ_order()[1:2] + circ_order()[3:4]
    procs= []
    for name in circ_names:
        procs.append(mp.Process(target=run_vqe, args=(name, vqe_data,)))
    for name in circ_names:
        procs.append(mp.Process(target=run_vqe, args=(name, vqe_data2,)))
    for proc in procs:
        proc.start()
    for proc in procs:
        proc.join()
        
    
