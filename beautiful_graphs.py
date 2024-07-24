
from Ternary_Tree.CircAlgorithm.UpUCCSDG import UpUCCSDG
from qiskit_algorithms import VQE
from functools import partial
from scipy.optimize import minimize
from qiskit.primitives import Estimator
from qiskit.quantum_info import Statevector
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit import QuantumCircuit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



from Ternary_Tree.OpenFermionSQ.Ternary_Tree_mapper import TernaryTreeMapper


basis="6-311+g"
names = 'H 0 0 0; Li 0 0 '
active_orbitals = [0, 1, 2, 3]

# active_orbitals = None

def init(dist=0.714):
    global ucc, fermionic_op
    ucc = UpUCCSDG(geometry=names + str(dist), basis=basis, active_orbitals=active_orbitals, num_electrons=(2, 1))
    fermionic_op = ucc.mol.hamiltonian.second_q_op()

def jw_mapping():
    mapper = JordanWignerMapper()
    qubit_jw_op = mapper.map(fermionic_op)
    return qubit_jw_op

def tt_mapping():
    print(ucc.tt)
    fermionic_op = ucc.mol.hamiltonian.second_q_op()
    mapper = TernaryTreeMapper(ucc.tt)
    qubit_tt_op = mapper.map(fermionic_op)
    return qubit_tt_op

def test():
    from numpy.random import random
    from qiskit import transpile
    qc = ucc.get_parametrized_circuit()
    # print(ucc.get_alpha_excitations())
    # print(ucc.tt)
    qc = QuantumCircuit.from_instructions(qc)
    # print(qc.draw(output='latex_source'))
    # print(qc)
    par = qc.parameters
    qc = qc.assign_parameters({el: random(1)[0] for el in par})
    circ = transpile(qc.decompose(reps=5))
    # print(transpile(qc.decompose(reps=2)))
    num = 0
    for inst in circ.data:
        if len(inst.qubits) == 2:
            num += 1
    
    print(f"num = {num}")
    print(f"depth = ", circ.depth() )
    print(len(ucc.maj_exc_par))
    return num, circ.depth(), len(ucc.maj_exc_par)

def plot_smth(datas):
    n_qubits = np.array([4,6,8,10,12,14,16,18,20,22,24,26])
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for index, data in enumerate(datas):
        cnot = data[0, :]
        depth = data[1, :]
        excitat = data[2, :]
        if index == 0:
            ax.plot(n_qubits*2, cnot/excitat, "--d",color="blue", markersize=6,label="JW")
        elif index == 1:
            ax.plot(n_qubits*2, cnot/excitat, "--^", color="blue", markersize=6,label="BK")
        elif index == 2: 
            ax.plot(n_qubits*2, cnot/excitat, "--v", color="blue", markersize=6,label="JW OPT")
        elif index == 3: 
            ax.plot(n_qubits*2, cnot/excitat, "--d", color="red", markersize=6,label="JW lexic")
        elif index == 4: 
            ax.plot(n_qubits*2, cnot/excitat, "--^", color="red", markersize=6,label="BK lexic")
        elif index == 5: 
            ax.plot(n_qubits*2, cnot/excitat, "--v", color="red", markersize=6,label="JW OPT lexic")
        elif index == 6: 
            ax.plot(n_qubits*2, cnot/excitat, "--*", color="green", markersize=6,label="Dynamic")
    ax.set_xlabel('Количество кубитов')
    ax.set_ylabel('Среднее количество CNOT гейтов')
    ax.set_xlim(4, 56)
    # ax.set_yscale('log')
    # ax.set_ylim(2, 4)
    ax.set_xticks(n_qubits*2)
    # ax.set_yticks([2.0, 2.5,3.0,3.5,4.0])
    ax.grid()
    # ax.plot(R[k:f],E_HF[k:f], label="Hartree Fock (classic initial energy)")
    ax.legend()
    # ax.set_title(label="Mоделирование молекулы H_2")
    plt.savefig("Wieght.png", dpi=400)

def plot_depth(datas):
    n_qubits = np.array([4,6,8,10,12,14,16,18,20,22,24,26])
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    for index, data in enumerate(datas):
        cnot = data[0, :]
        depth = data[1, :]
        excitat = data[2, :]
        if index == 0:
            ax.plot(n_qubits*2, depth, "--d",color="blue", markersize=6,label="JW")
        elif index == 1:
            ax.plot(n_qubits*2, depth, "--^", color="blue", markersize=6,label="BK")
        elif index == 2: 
            ax.plot(n_qubits*2, depth, "--v", color="blue", markersize=6,label="JW OPT")
        elif index == 3: 
            ax.plot(n_qubits*2, depth, "--d", color="red", markersize=6,label="JW lexic")
        elif index == 4: 
            ax.plot(n_qubits*2, depth, "--^", color="red", markersize=6,label="BK lexic")
        elif index == 5: 
            ax.plot(n_qubits*2, depth, "--v", color="red", markersize=6,label="JW OPT lexic")
        elif index == 6: 
            ax.plot(n_qubits*2, depth, "--*", color="green", markersize=6,label="Dynamic")
    ax.set_xlabel('Количество кубитов')
    ax.set_ylabel('Глубина 1 слоя k-UpUCCSD анзаца')
    ax.set_yscale('log')
    ax.set_xlim(4, 56)
    # ax.set_ylim(2, 4)
    ax.set_xticks(n_qubits*2)
    # ax.set_yticks([2.0, 2.5,3.0,3.5,4.0])
    ax.grid()
    # ax.plot(R[k:f],E_HF[k:f], label="Hartree Fock (classic initial energy)")
    ax.legend()
        # ax.set_title(label="Mоделирование молекулы H_2")
    plt.savefig("Depth.png", dpi=400)


import numpy as np

n = 10
displ = 0
data = np.empty([3,n])
for n in range(1, n+1):
    # active_orbitals = [i for i in range(displ, n*2 + displ)]
    # init(dist=1)
    # data[:, n - 1] = np.array(test())
    pass




data = np.array([[28, 162, 348, 704, 1022, 1446, 1932, 2672, 3250, 3966],
                [35,  99, 175, 241, 315, 377, 429, 511,  593,  653],
                [8,  64, 152, 272, 424, 608, 824, 1072, 1352, 1664]])
# plot_smth(data)
# plot_depth(data)

N = [4,6,8,10,12,14]
num_temrs = [64, 152, 272, 424, 608, 824]

# data_jw = np.array([[204, 503, 948, 1525, 2236, 3084, 4068],
#                     [328, 633, 1516, 2424, 3568, 4918, 6484],
#                     [72, 180, 336, 540, 792, 1092, 1440]])

# data_bk = np.array([[222, 772, 1484, 2995, 4631, 6530, 8388],
#                     [325, 1053, 1950, 3795, 5719, 8043, 10065],
#                     [72, 180, 336, 540, 792, 1092, 1440]])

# data_thebest = np.array([[146, 358, 656, 1048, 1526, 2090, 2744],
#                          [82, 151, 284, 268, 301, 359, 384],
#                          [72, 180, 336, 540, 792, 1092, 1440]])

# plot_depth([data_jw, data_bk, data_thebest])
# plot_smth([data_jw, data_bk, data_thebest])

cx_lex_jw = [187, 481, 903, 1349, 2131, 2937, 3871]
depth_lex_jw = [296, 728, 1328, 2096, 3032, 4136, 5408]

cx_lex_bk = [214, 750, 1398, 2799, 4306, 6192, 7831]
depth_lex_bk = [313, 1017, 1823, 3459, 5245, 7476, 9205]



data_bk_lex = np.array([[222, 772, 1484, 2995, 4631, 6530, 8388, 12441, 16202, 20472, 24450, 29804],
                        [325, 1053, 1950, 3792, 5719, 8043, 10065, 14830, 19125, 24094, 27563, 34601],
                        [72, 180, 336, 540, 792, 1092, 1440, 1836, 2280, 2772, 3312, 3900]])
data_bk = np.array([[340, 1396, 2622, 5259, 8384, 11878, 15614, 22813, 30081, 38219, 46761, 56091],
                    [455, 1728, 2966, 6639, 10269, 14605, 17068, 27272, 36068, 45796, 55014, 66460],
                    [72, 180, 336, 540, 792, 1092, 1440, 1836, 2280, 2772, 3312, 3900]])

data_jw_opt_lex = np.array([[204, 508, 948, 1524, 2236, 3084, 4068, 5188, 6444, 7836, 9364, 11028],
                        [328, 814, 1516, 2434, 3568, 4918, 6484, 8266, 10264, 12478, 14908, 17554],
                        [72, 180, 336, 540, 792, 1092, 1440, 1836, 2280, 2772, 3312, 3900]])
data_jw_opt = np.array([[325, 977, 2117, 3886, 6399, 9797, 14195, 19734, 26529, 34721, 44425, 55782],
                    [391, 1025, 1951, 3455, 5439, 8237, 11737, 16259, 21745, 28449, 36362, 45665],
                    [72, 180, 336, 540, 792, 1092, 1440, 1836, 2280, 2772, 3312, 3900]])

data_jw_lex = np.array([[202, 569, 1196, 2147, 3486, 5277, 7584, 10471, 14002, 18241, 23252, 29099],
                            [322, 889, 1812, 3155, 4982, 7357, 10344, 14007, 18410, 23617, 29692, 36699],
                            [72, 180, 336, 540, 792, 1092, 1440, 1836, 2280, 2772, 3312, 3900]])
data_jw = np.array([[351, 1223, 2510, 4850, 8326, 13162, 19582, 27810, 38070, 50586, 65582, 83282],
                    [481, 1473, 3235, 5929, 9761, 14914, 21591, 29986, 40294, 52686, 67372, 84539],
                    [72, 180, 336, 540, 792, 1092, 1440, 1836, 2280, 2772, 3312, 3900]])
data_thebest = np.array([[146, 358, 656, 1048, 1526, 2090, 2744, 3492, 4326, 5246, 6258, 7352],
                         [82, 151, 284, 268, 301, 359, 384, 501, 516, 583, 599, 682],
                         [72, 180, 336, 540, 792, 1092, 1440, 1836, 2280, 2772, 3312, 3900]])

plot_depth([data_jw, data_bk, data_jw_opt, data_jw_lex, data_bk_lex, data_jw_opt_lex, data_thebest])
plot_smth([data_jw, data_bk, data_jw_opt, data_jw_lex, data_bk_lex, data_jw_opt_lex, data_thebest])



# bk_lex_depth:  [325, 1053, 1950, 3792, 5719, 8043, 10065, 14830, 19125, 24094, 27563, 34601]
# bk_lex_cx:  [222, 772, 1484, 2995, 4631, 6530, 8388, 12441, 16202, 20472, 24450, 29804]
# bk_depth:  [455, 1728, 2966, 6639, 10269, 14605, 17068, 27272, 36068, 45796, 55014, 66460]
# bk_depth:  [340, 1396, 2622, 5259, 8384, 11878, 15614, 22813, 30081, 38219, 46761, 56091]
# num_maj:  [72, 180, 336, 540, 792, 1092, 1440, 1836, 2280, 2772, 3312, 3900]
# 
# # jw_lex_depth:  [328, 814, 1516, 2434, 3568, 4918, 6484, 8266, 10264, 12478, 14908, 17554]
# jw_lex_cx:  [204, 508, 948, 1524, 2236, 3084, 4068, 5188, 6444, 7836, 9364, 11028]
# jw_depth:  [391, 1025, 1951, 3455, 5439, 8237, 11737, 16259, 21745, 28449, 36362, 45665]
# jw_cx:  [325, 977, 2117, 3886, 6399, 9797, 14195, 19734, 26529, 34721, 44425, 55782]
# num_maj:  [72, 180, 336, 540, 792, 1092, 1440, 1836, 2280, 2772, 3312, 3900]

# jw_lex_depth:  [322, 889, 1812, 3155, 4982, 7357, 10344, 14007, 18410, 23617, 29692, 36699]
# jw_lex_cx:  [202, 569, 1196, 2147, 3486, 5277, 7584, 10471, 14002, 18241, 23252, 29099]
# jw_depth:  [481, 1473, 3235, 5929, 9761, 14914, 21591, 29986, 40294, 52686, 67372, 84539]
# jw_cx:  [351, 1223, 2510, 4850, 8326, 13162, 19582, 27810, 38070, 50586, 65582, 83282]
# num_maj:  [72, 180, 336, 540, 792, 1092, 1440, 1836, 2280, 2772, 3312, 3900]