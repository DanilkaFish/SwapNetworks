from typing import List

from qiskit_algorithms.optimizers import SPSA ,CG, SLSQP, L_BFGS_B, COBYLA

import json

from Ternary_Tree.optimizer.soap import SOAP
from Ternary_Tree.UCC.AbstractUCC import Molecule
from utils import *


H2_4 = Molecule(geometry='H 0 0 0; H 0 0 0.7349', num_electrons=(1,1), active_orbitals=[0,1], basis='sto-3g')
H2_8 = Molecule(geometry='H 0 0 0; H 0 0 0.7349', num_electrons=(1,1), active_orbitals=[0,1,2,3], basis='6-31g')
LiH_8 = Molecule(geometry='H 0 0 0; Li 0 0 1.5459', num_electrons=(2,2), active_orbitals=[0,1,2,5], basis='sto-3g')


optimizers = [(SLSQP(maxiter=100, ftol=0), 'SLSQP')]


def run_noisy_vqe(file_name_to_write: str,
                  file_name_to_read: str,
                  molecule: Molecule,
                  optimizer: any,
                  circ_names: List[str]=circ_order(),
                  noise_type: str="D",
                  probs: np.ndarray=np.flip(np.geomspace(0.00002, (0.001), 7)),
                  reps: int=1
                  ):
    data = []
    with open(file_name_to_write, 'w') as file:
        with open(file_name_to_read, 'r') as rf:
            datajson = json.load(rf)
            circ_prov = CircuitProvider(reps=reps, molecule=molecule)
            ref_value = numpy_energy(circ_prov.fermionic_op, circ_prov.ucc)
            for name, circ, op in circ_prov(circ_names):
                for prob in probs:
                    init_point = datajson[index]["param"]
                    print(init_point)
                    circs = CircSim(circ, op, True, 1 - prob, noise_type, q1_noise=False, q2_noise=True, init_point=init_point)
                    energy, parameters = circs.run_qiskit_vqe(optimizer[0])
                    data.append({"name": name, 
                                 "ref_ener": ref_value, 
                                 "energy": energy, 
                                 "param" : parameters[-1], 
                                 "optimizer": optimizer[1], 
                                 "prob": prob, 
                                 "noise gates": circ.count_ops()["cx"]})
                index += 1
            json.dump(data, file, indent=4)

# def write_data():
#     noises = ["D", "X", "Y", "Z"]
#     id_en = [None]*n
#     for noise in noises:
#         with open(f'data/NewSimPar{noise}8.json', 'w') as file:
#             data = []
#             with open('data/SimParId8.json', 'r') as rf:
#                 index = 0
#                 datajson = json.load(rf)
#                 for mol in molecules[1:2]:
#                     for opt in optimizers[0:1]:
#                         for rep in reps:
#                             circ_prov = CircuitProvider(reps=rep, active_orbitals=mol[2], num_electrons=mol[1], geometry=mol[0], basis=mol[3])
#                             ref_value = numpy_energy(circ_prov.fermionic_op, circ_prov.ucc)
#                             for name, circ, op in circ_prov.__iter__(circ_order()):
#                                 for prob in np.flip(np.geomspace(0.00002, (0.001), 7)):
#                                     init_point = datajson[index]["param"]
#                                     print(init_point)
#                                     circs = CircSim(circ, op, True, 1 - prob, noise, q1_noise=False, q2_noise=True, init_point=init_point)
#                                     energy, _, values, parameters = circs.run_qiskit_vqe(opt[0])
#                                     data.append({"name": name, "ref_ener": ref_value, "energy": energy, "param" : parameters[-1], "optimizer": opt[1], "prob": prob, "noise gates": circ.count_ops()["cx"]})
#                                 index += 1
#             json.dump(data, file, indent=4)

# def pure_write_data():
#     noises = ["Id"]
#     import json
#     id_en = [None]*n
#     for noise in noises:
#         with open(f'data/SimPar{noise}8.json', 'w') as file:
#             data = []
#             index = 0
#             for mol in molecules[1:2]:
#                 for opt in optimizers[0:1]:
#                     for rep in reps:
#                         iter_num = [None]*n
#                         # circs, ops, ref_value, init_point, enhf, circ_prov = get_circs_and_ops(get_cp=True, reps=rep, active_orbitals=mol[2], num_electrons=mol[1], geometry=mol[0], basis=mol[3])
#                         circ_prov = CircuitProvider(reps=rep, active_orbitals=mol[2], num_electrons=mol[1], geometry=mol[0], basis=mol[3])
#                         ref_value = numpy_energy(circ_prov.fermionic_op, circ_prov.ucc)
#                         # circs.append(circ_prov.get_yordan_dynamic())
#                         # ops.append(ops[0])
#                         # circs.append(circ_prov.get_zyx_yordan_dynamic())
#                         # ops.append(ops[1])
#                         for name, circ, op in circ_prov.__iter__(circ_order()):
#                             circs = CircSim(circ, op, False)
#                             # energy, _, values, parameters = [0,0,[0],[0]]
#                             energy, _, values, parameters = circs.run_qiskit_vqe(opt[0])
#                             data.append({"name": name, "ref_ener": ref_value, "energy": energy, "param" : parameters[-1], "optimizer": opt[1]})
#                             index += 1
#             json.dump(data, file, indent=4)

# pure_write_data()
# write_data()

