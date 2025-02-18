from Ternary_Tree.ucc.upgccsd import UpGCCSD, LadExcNames
from qiskit.quantum_info import Statevector
from qiskit import transpile
from Ternary_Tree.ucc.abstractucc import Molecule
import numpy as np

from utils import CircuitProvider

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
    double = []
    single = []
    depth = []
    num = []
    pauli_num = []
    for i in range(1,2):
        active_orbitals = [j for j in range(2*i)]
        ucc = UpGCCSD(geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons)
        cirq = ucc.swap_gen()
        cirq = transpile(cirq, basis_gates=basis_gates, optimization_level=3)
        # print(cirq)
        num.append(i*4)
        double.append(cirq.count_ops()["cx"])
        single.append(cirq.count_ops()["u3"])
        depth.append(cirq.depth())
        pauli_num.append(2*i*(2*i-1)/2*10)
        print(cirq.count_ops())
        print(cirq.depth())

    print(num)
    print(double)
    print(single)
    print(depth)
    print(pauli_num)


# def print_circuit():
#     mol = ('Be 0 0 0; Be 0 0 0.7349', (1, 1), [0,1,2,3,4,5], "6-311G")
#     # molecules = [('H 0 0 0; Li 0 0 1.5459', (2,2), [0,1,2,5], basises[0])]
#     reps = [1, 2, 3]
#     rep = 1
#     cp = CircuitProvider(reps=rep, active_orbitals=mol[2], num_electrons=mol[1], geometry=mol[0], basis=mol[3])
#     cp.get_circ_with_mapping(TernaryTreeMapper(cp.ucc.get_jw_opt()), lexic=True, init=False)
# scaleability()
base_test()
# print_circuit()
# 2xn lattice
# num = np.array([4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52])
# double = np.array([14, 100, 246, 456, 730, 1068, 1470, 1936, 2466, 3060, 3718, 4440, 5226])
# single = np.array([26, 178, 429, 788, 1255, 1830, 2513, 3304, 4203, 5210, 6325, 7548, 8879])
# depth = np.array([15, 67, 99, 131, 163, 195, 227, 259, 291, 323, 355, 387, 419])
# pauli_num = np.array([10.0, 60.0, 150.0, 280.0, 450.0, 660.0, 910.0, 1200.0, 1530.0, 1900.0, 2310.0, 2760.0, 3250.0])

# yordan
# num = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52]
# double = [15, 136, 327, 598, 949, 1380, 1891, 2482, 3153, 3904, 4735, 5646, 6637]
# single = [30, 262, 620, 1130, 1792, 2606, 3572, 4690, 5960, 7382, 8956, 10682, 12560]
# depth = [25, 126, 192, 258, 324, 390, 456, 522, 588, 654, 720, 786, 852]
# pauli_num = [10.0, 60.0, 150.0, 280.0, 450.0, 660.0, 910.0, 1200.0, 1530.0, 1900.0, 2310.0, 2760.0, 3250.0]

# short
# num = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52]
# double = [16, 142, 342, 626, 994, 1446, 1982, 2602, 3306, 4094, 4966, 5922, 6962]
# single = [25, 247, 582, 1059, 1678, 2439, 3342, 4387, 5574, 6903, 8374, 9987, 11742]
# depth = [18, 99, 152, 205, 258, 311, 364, 417, 470, 523, 576, 629, 682]
# pauli_num = [10.0, 60.0, 150.0, 280.0, 450.0, 660.0, 910.0, 1200.0, 1530.0, 1900.0, 2310.0, 2760.0, 3250.0]

# for i in range(0, 13):
#     print("(", num[i],",", double[i]/pauli_num[i], ")", end=" ")


# res(0.78)