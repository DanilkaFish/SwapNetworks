from Ternary_Tree.ucc.upgccsd import UpGCCSD
from qiskit.quantum_info import Statevector
from qiskit import transpile
from Ternary_Tree.ucc.abstractucc import Molecule
from utils import CircuitProvider, Circuits
import numpy as np


basis_gates = ["u3", "cx"]
# basis_gates1 = [ "h","rz","cx"]
# basis_gates1 = basis_gates
# basis_gates2 = ["u3", 'cx']

basis = '6-311G'
# geometry='Be 0 0 0; Be 0 0 0.739'
geometry='H 0 0 0; H 0 0 0.739'
active_orbitals=[0, 1]
num_electrons=(1, 1)
H2_4 = Molecule(geometry='H 0 0 0; H 0 0 0.7349', num_electrons=(1,1), active_orbitals=[0,1], basis='sto-3g')

def base_test():
    ucc = UpGCCSD(molecule=H2_4)
    cirq, mtoq = ucc.swap2xn(1)
    # cirq = ucc.swap_gen(1, LadExcNames.YORDAN())
    print(cirq)
    params = cirq.parameters
    params = {par: 1 for par in params}
    cirq.assign_parameters(params, inplace=True)
    # print(cirq.decompose())
    state = Statevector.from_label("0000")
    res = state.evolve(cirq)
    probs = res.probabilities_dict()
    for key in probs:
        if abs(probs[key]) > 0.000001:
            print(key, ": ", probs[key])
    print(probs)

def res(angle):
    active_orbitals = [j for j in range(2)]
    ucc = UpGCCSD(geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons)
    cirq = ucc.swap_gen()
    state = Statevector.from_label("0000")
    par = cirq.parameters[0]
    # cirq.assign_parameters({par:angle})
    cirq1 = cirq.bind_parameters([0])
    print(state.evolve(cirq1).probabilities_dict())
    cirq = cirq.bind_parameters([angle])
    state = state.evolve(cirq)
    print(cirq.decompose())
    print(state.probabilities_dict())
    
def scaleability():
    mol = Molecule(geometry='Li 0 0 0; Li 0 0 1.5459', num_electrons=(2,2), active_orbitals=[0,1,2,5], basis='6-311g')
    circ_names = Circuits.get_circs_names()[2:4]
    double = {circ: [] for circ in circ_names}
    single = {circ: [] for circ in circ_names}
    depth = {circ: [] for circ in circ_names}
    num = {circ: [] for circ in circ_names}
    pauli_num = {circ: [] for circ in circ_names}
    for i in range(4,28,2):
        mol.active_orbitals = list(range(i))
        try:
            circ_prov = CircuitProvider(reps=1, molecule=mol)
            for name in circ_names:
                circ = circ_prov.get_circ(name)[1]
                circ = transpile(circ, basis_gates=basis_gates, optimization_level=3)
                # print(cirq)
                # num.append(i*4)
                # print(circ.count_ops())
                double[name].append(circ.count_ops()["cx"])
                single[name].append(circ.count_ops()["u3"])
                depth[name].append(circ.depth())
                pauli_num[name].append(i*(i-1)/2*10)
                # pauli_num.append(2*i*(2*i-1)/2*10)
                # print(name, ":", circ.count_ops())
                # print(name, ":", circ.depth())
        except IndexError:
            pass
    for name in circ_names:
        print(name, " double:", double[name])
        print(name, " single:", single[name])
        print(name, "depth:", depth[name])

    print(Circuits.get_circs_names()[0], "pauli num:", pauli_num[name])
    # print(num)
    # print(double)
    # print(single)
    # print(depth)
    # print(pauli_num)


# def print_circuit():
#     mol = ('Be 0 0 0; Be 0 0 0.7349', (1, 1), [0,1,2,3,4,5], "6-311G")
#     # molecules = [('H 0 0 0; Li 0 0 1.5459', (2,2), [0,1,2,5], basises[0])]
#     reps = [1, 2, 3]
#     rep = 1
#     cp = CircuitProvider(reps=rep, active_orbitals=mol[2], num_electrons=mol[1], geometry=mol[0], basis=mol[3])
#     cp.get_circ_with_mapping(TernaryTreeMapper(cp.ucc.get_jw_opt()), lexic=True, init=False)
scaleability()
# base_test()
# print_circuit()



# jw  double: [413, 1340, 3087, 5910, 10065, 15808, 23395, 33082, 45125, 59780, 77303, 97950]
# jw  single: [370, 963, 1832, 2977, 4398, 6095, 8068, 10317, 12842, 15643, 18720, 22073]
# jw depth: [527, 1550, 3346, 6133, 10103, 15478, 22457, 31244, 42070, 55093, 70602, 88710]
# bk  double: [394, 1528, 2881, 6064, 9546, 13594, 17771, 25778, 34489, 43765, 53287, 64024]
# bk  single: [298, 876, 1583, 2700, 4046, 5623, 7222, 9428, 11922, 14560, 17523, 20663]
# bk depth: [497, 1840, 3257, 7118, 10723, 15592, 18397, 29528, 38484, 48184, 57396, 70769]
# jw_lex  double: [149, 392, 819, 1339, 2004, 2994, 3996, 5630, 6903, 9276, 11351, 14902]
# jw_lex  single: [181, 416, 888, 1403, 2071, 3021, 4063, 5523, 6820, 8864, 10597, 13873]
# jw_lex depth: [204, 477, 882, 1350, 1980, 2811, 3750, 5097, 5885, 7876, 9318, 12194]
# bk_lex  double: [126, 407, 804, 1406, 2269, 3109, 4413, 6111, 7957, 10523, 12695, 19096]
# bk_lex  single: [148, 428, 784, 1458, 2110, 3079, 4194, 6062, 7234, 9845, 11714, 17083]
# bk_lex depth: [164, 498, 866, 1439, 2223, 2911, 4158, 5512, 6684, 8643, 10680, 15599]
# swap 2xn  double: [100, 246, 456, 730, 1068, 1470, 1936, 2466, 3060, 3718, 4440, 5226]
# swap 2xn  single: [174, 420, 772, 1230, 1794, 2464, 3240, 4122, 5110, 6204, 7404, 8710]
# swap 2xn depth: [67, 99, 131, 163, 195, 227, 259, 291, 323, 355, 387, 419]
# swap gen short  double: [142, 342, 626, 994, 1446, 1982, 2602, 3306, 4094, 4966, 5922, 6962]
# swap gen short  single: [240, 570, 1040, 1650, 2400, 3290, 4320, 5490, 6800, 8250, 9840, 11570]
# swap gen short depth: [99, 152, 205, 258, 311, 364, 417, 470, 523, 576, 629, 682]
# swap gen yor  double: [136, 327, 598, 949, 1380, 1891, 2482, 3153, 3904, 4735, 5646, 6637]
# swap gen yor  single: [260, 618, 1128, 1790, 2604, 3570, 4688, 5958, 7380, 8954, 10680, 12558]
# swap gen yor depth: [126, 192, 258, 324, 390, 456, 522, 588, 654, 720, 786, 852]
# jw pauli num: [60.0, 150.0, 280.0, 450.0, 660.0, 910.0, 1200.0, 1530.0, 1900.0, 2310.0, 2760.0, 3250.0]



# jw_lex  double: [127, 351, 737, 1237, 1894, 2865, 3898, 5451, 6605, 8906, 10882, 14075]
# jw_lex  single: [147, 349, 784, 1286, 1926, 2857, 3908, 5306, 6501, 8475, 10125, 13134]
# jw_lex depth: [177, 437, 802, 1267, 1871, 2694, 3639, 4929, 5641, 7571, 8939, 11586]
# bk_lex  double: [112, 372, 748, 1301, 2168, 2936, 4156, 5915, 7547, 10021, 12209, 18485]
# bk_lex  single: [122, 376, 710, 1331, 1977, 2868, 3936, 5833, 6824, 9392, 11243, 16497]
# bk_lex depth: [148, 458, 817, 1348, 2118, 2791, 3968, 5355, 6373, 8249, 10319, 15205]
# jw pauli num: [60.0, 150.0, 280.0, 450.0, 660.0, 910.0, 1200.0, 1530.0, 1900.0, 2310.0, 2760.0, 3250.0]