from qiskit_algorithms.optimizers import SPSA ,CG, SLSQP, L_BFGS_B ,COBYLA
from qiskit.quantum_info import Statevector
from qiskit import QuantumCircuit
from qiskit.quantum_info.operators import SparsePauliOp
from qiskit.circuit.library.standard_gates import IGate, XGate, ZGate, YGate

from Ternary_Tree.optimizer.soap import SOAP

from utils import *

def f(*args, **kwargs):
    return None
# SPSA ,CG, SLSQP, L_BFGS_B ,COBYLA = [f for _ in range(5)]
# SOAP = f

# a0=0.07
# af=0.07
# # af = 0    
# b0 = 0.00001
# bf = 0.00001
# bf = 0
# # n = 12


a0=0.1
af=0.004
# af = 0
b0 = 0.012
bf = 0.0004
# bf = 0
m = 1
def a(i=0):
    def _a(i=0):
        while True:
            i += 1
            _a0 = a0/(0.2*m+1)
            _af = af/(0.2*m+1)
            yield (_a0 - _af)/i**(0.602) + _af
    return _a()
def b():
    def _b(i=0):
        while True:
            i += 1
            _bf = bf/(2*m+1)
            _b0 = b0/(2*m+1)
            yield (_b0 - _bf)/i**(0.101) + _bf
    return _b()

counts = [0]
def call(*args, **kwargs):
    counts[0] += 1
# choice of optimizer
optimizer = SLSQP(maxiter=5000)
optimizer = CG(maxiter=50, eps=0.0001, gtol=0.000000001)
# optimizer = L_BFGS_B(maxiter=300, eps=0.00001,ftol=0.000)

# optimizer = SPSA(maxiter=800, learning_rate=a, perturbation=b,  callback=call)
# optimizer = SLSQP(maxiter=50, ftol=0)
# optimizer = COBYLA(maxiter=1000, tol=0)
# optimizer = SOAP(500)  

# circuit_init
number_of_reps = 1
# circ_prov = CircuitProvider(reps=number_of_reps, active_orbitals=[0,1,2,3])




# en_hf = hf(thebest_qc)
id_ener = []
noise_ener = []
k = 0
n = 7+3

id_ener = [[] for i in range(n)]
noise_ener = [[] for i in range(n)]

_counts = [0 for i in range(n)]
_values = [0 for i in range(n)]


basises = ['sto-3g', '6-31g']
# basises = ['sto-3g']
geometry='H 0 0 0; Li 0 0 0.7349'
active_orbitals=[0, 1]
num_electrons=(1, 1)
# SPSA(maxiter=3000, learning_rate=a, perturbation=b, callback=call, perturbation_dims=10),
optimizers = [(SLSQP(maxiter=100, ftol=0), 'SLSQP')]
            #   (SLSQP(maxiter=800, ftol=0), 'SLSQP'), 
                # ]
                # (COBYLA(maxiter=1000, tol=0), 'COBYLA'),
                # (SPSA(maxiter=500, learning_rate=a, perturbation=b,  callback=call), 'SPSA'),
                # (SOAP(4000), "SOAP")]

# optimizers = [(SPSA(maxiter=500, learning_rate=a, perturbation=b,  callback=call), 'SPSA')]
# optimizers = [(SOAP(4000), "SOAP")]
# optimizers = [(L_BFGS_B(maxiter=1000, eps=0.0000001,ftol=0.0000000001), "LBFGSB")]
# optimizers = [(SLSQP(maxiter=500, ftol=0), 'SLSQP')]


molecules = [['H 0 0 0; H 0 0 0.7349', (1,1), [0,1], basises[0]], ['H 0 0 0; H 0 0 0.7349', (1,1), [0,1,2,3], basises[1]], 
                ['H 0 0 0; Li 0 0 1.5459', (2,2), [0,1,2,5], basises[0]]]
molecules = [('H 0 0 0; H 0 0 0.7349', (1,1), [0,1], basises[1]), ('H 0 0 0; H 0 0 0.7349', (1,1), [0,1,2,3], basises[1])]
# molecules = [('H 0 0 0; Li 0 0 1.5459', (2,2), [0,1,2,5], basises[0])]
reps = [1, 2, 3]
reps = [1]

def write_data():
    noises = ["D", "X", "Y", "Z"]
    import json
    id_en = [None]*n
    for noise in noises:
        with open(f'data/NewSimPar{noise}8.json', 'w') as file:
            data = []
            with open('data/SimParId8.json', 'r') as rf:
                index = 0
                datajson = json.load(rf)
                for mol in molecules[1:2]:
                    for opt in optimizers[0:1]:
                        for rep in reps:
                            iter_num = [None]*n
                            # circs, ops, ref_value, init_point, enhf, circ_prov = get_circs_and_ops(get_cp=True, reps=rep, active_orbitals=mol[2], num_electrons=mol[1], geometry=mol[0], basis=mol[3])
                            circ_prov = CircuitProvider(reps=rep, active_orbitals=mol[2], num_electrons=mol[1], geometry=mol[0], basis=mol[3])
                            ref_value = numpy_energy(circ_prov.fermionic_op, circ_prov.ucc)
                            # circs.append(circ_prov.get_yordan_dynamic())
                            # ops.append(ops[0])
                            # circs.append(circ_prov.get_zyx_yordan_dynamic())
                            # ops.append(ops[1])
                            for name, circ, op in circ_prov.__iter__(circ_order()):
                                for prob in np.flip(np.geomspace(0.00002, (0.001), 7)):
                                    init_point = datajson[index]["param"]
                                    print(init_point)
                                    circs = CircSim(circ, op, True, 1 - prob, noise, q1_noise=False, q2_noise=True, init_point=init_point)
                                
                                    energy, _, values, parameters = circs.run_qiskit_vqe(opt[0])
                                    # data.append({"mol": mol, "act_orb": mol[2], "ref_ener": ref_value, "energy": id_en[index], "iter_num": iter_num, "param" : parameters[-1]})
                                    data.append({"name": name, "ref_ener": ref_value, "energy": energy, "param" : parameters[-1], "optimizer": opt[1], "prob": prob, "noise gates": circ.count_ops()["cx"]})
                                # circs = CircSim(circ, op, False)
                                # # energy, _, values, parameters = [0,0,[0],[0]]
                                # energy, _, values, parameters = circs.run_qiskit_vqe(opt[0])
                                # data.append({"name": name, "ref_ener": ref_value, "energy": energy, "param" : parameters[-1], "optimizer": opt[1]})
                                index += 1
            json.dump(data, file, indent=4)

def pure_write_data():
    noises = ["Id"]
    import json
    id_en = [None]*n
    for noise in noises:
        with open(f'data/SimPar{noise}8.json', 'w') as file:
            data = []
            index = 0
            for mol in molecules[1:2]:
                for opt in optimizers[0:1]:
                    for rep in reps:
                        iter_num = [None]*n
                        # circs, ops, ref_value, init_point, enhf, circ_prov = get_circs_and_ops(get_cp=True, reps=rep, active_orbitals=mol[2], num_electrons=mol[1], geometry=mol[0], basis=mol[3])
                        circ_prov = CircuitProvider(reps=rep, active_orbitals=mol[2], num_electrons=mol[1], geometry=mol[0], basis=mol[3])
                        ref_value = numpy_energy(circ_prov.fermionic_op, circ_prov.ucc)
                        # circs.append(circ_prov.get_yordan_dynamic())
                        # ops.append(ops[0])
                        # circs.append(circ_prov.get_zyx_yordan_dynamic())
                        # ops.append(ops[1])
                        for name, circ, op in circ_prov.__iter__(circ_order()):
                            circs = CircSim(circ, op, False)
                            # energy, _, values, parameters = [0,0,[0],[0]]
                            energy, _, values, parameters = circs.run_qiskit_vqe(opt[0])
                            data.append({"name": name, "ref_ener": ref_value, "energy": energy, "param" : parameters[-1], "optimizer": opt[1]})
                            index += 1
            json.dump(data, file, indent=4)

# pure_write_data()
write_data()

# Noise_ops = [(None, "D"), (ZGate(), "Z"), (XGate(), "X"), (YGate(),"Y")]
Noise_ops = [(None, "D")]

# def test():
    # import numpy as np
    # from qiskit import QuantumCircuit, transpile
    # import json

    # f = open('./numbers.json')
    # data_ideal = json.load(f)
    # # for d in data_ideal:
    # #     print(d)
    # # print(data_ideal)
    # print(len(data_ideal))
    # X = np.linspace(0, 1, 500)
    # E = []
    # toto = []
    # l = 0
    # EE = {} 

    # for mol in molecules[:]:
    #     for opt in optimizers[:]:
    #         for rep in reps:
    #             iter_num = [None]*n
    #             circs, ops, ref_value, init_point, en_hf, circ_prov = get_circs_and_ops(get_cp=True, reps=rep, active_orbitals=mol[2], num_electrons=mol[1], geometry=mol[0], basis=mol[3])
    #             circs.append(circ_prov.get_yordan_dynamic())
    #             ops.append(ops[0])
    #             # a,b,c,d = ideal_energy(circ, op, ref_value, opt[0], init_point,  en_hf)
    #             # print(a,c[-1],d[-1])
    #             # print(circ)
    #             k = 0
    #             for circ, op in zip(circs[0:1], ops[0:1]):
    #                 k += 1
    #                 id_par = data_ideal[l]['param']
    #                 par = circ.parameters
    #                 for noise_op, name in Noise_ops:
    #                     E = []
    #                     ansatz = circ.assign_parameters({el: id_par[i] for i, el in enumerate(par)})
    #                     print(ansatz)
    #                     print(op)
    #                     # energy, _, values, paramters = ideal_energy(circ, op, ref_value, opt[0], init_point,  en_hf)
    #                     for prob in [0.99, 0.999,0.9999, 0.99999, 0.999999, 0.9999999]:
    #                     # for prob in [0.99999, 0.999999, 0.9999999]:
    #                         est = get_device_noise_estimator(8, noise_op=noise_op, prob=prob)

    #                         # x = 0.11174828549146505
    #                         # print({el: id_par[i] for i, el in enumerate(par)})
    #                         energy = est.run(ansatz, op).result()
    #                         E.append(energy.values[0])

    #                         # energy, _, values, _ = VQE_energy_with_noise(circ, op, ref_value, opt[0], init_point,  en_hf, est)
    #                         # E.append(energy)
    #                     EE[name + "_" + str(k)] = E
    #                 l += 1

    #             for key in EE:
    #                 print(key, ": ", EE[key])

    # return EE

# EE = test()
# print(EE)
# import pylab    
# pylab.xlabel("error")
# pylab.ylabel("Energy")
# # pylab.title("Convergence with no noise")
# # pylab.ylim(bottom=-1.90, top=-1.30)
# # pylab.rcParams["figure.figsize"] = (12, 4)

# def print_last_data(counts, values, name="hello"):
#     pylab.plot(counts, values, label=name)
    
# # print_last_data(X, E, name=str(i))
# # print(min(E))
# # pylab.legend()
# # pylab.savefig("toto2")


# a = [[-1.0031434741388943, -1.6905585840512296, -1.8393852284472576, -1.855560165924037, -1.8571912570353855],
#     [-1.444400023964278, -1.803600922873114, -1.8517715150181235, -1.856809923221062, -1.8573163472207639],
#     [-1.003119841516724, -1.690559509101338, -1.839385656780147, -1.855560541931011, -1.8571916277828961]]
# probs = [1 - prob for prob in  [0.9, 0.99,0.999,0.9999, 0.99999]]
# labels = ["12 gates", "38 gates", "fixed parms"]
# for i in range(3):
#     pylab.plot(probs[1:], a[i][1:], label=labels[i])
# pylab.xscale('log')
# pylab.legend()
# pylab.savefig("toto3")



