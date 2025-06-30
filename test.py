from Ternary_Tree.ucc.upgccsd import UpGCCSD, UCCSD
from qiskit.quantum_info import Statevector
from qiskit import transpile
from Ternary_Tree.ucc.abstractucc import Molecule
from Ternary_Tree.qiskit_interface.circuit_provider import CircuitProvider, Circuits
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
H2_8 = Molecule(geometry='H 0 0 0; H 0 0 0.7349', num_electrons=(1,1), active_orbitals=[0,1,2,3], basis='6-31g')

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
    circ_names = Circuits.get_circs_names()[0:4]
    double = {circ: [] for circ in circ_names}
    single = {circ: [] for circ in circ_names}
    depth = {circ: [] for circ in circ_names}
    num = {circ: [] for circ in circ_names}
    pauli_num = {circ: [] for circ in circ_names}
    for i in range(4,14,2):
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

if __name__ == "__main__":
    # scaleability()
    ucc = UCCSD(molecule=H2_8)
    circ, mapping = ucc.cyclic_algorithm()
    # pars = circ.parameters[:]
    # qc = circ.to_excitaions(pars)
    print(circ)
    # na = 3
    # cycl_double = []
    # cycl_depth = []
    # for i in range(4,16,2):
    #     nb = i - na
    #     s = na*nb*2
    #     sg = i*(i - 1)
    #     d = 2*na*(na-1)*nb*(nb-1)/4 + na*na*nb*nb
    #     dg = 2*i*(i - 1)*i*(i-1)/4 + i*(i-1)*i*(i-1)
    #     # print(s)
    #     # print(d)
    #     print(s*2 + d*8)
    #     swap = (4/6*i*i*i + i*i  + 1/3*i - 14)*4
    #     cycl_double.append(swap + 12*d)
    #     cycl_depth.append(swap/4*7 + 10*d)

    # print("cycl_depth",cycl_depth)
    # print("cycl_double",cycl_double)