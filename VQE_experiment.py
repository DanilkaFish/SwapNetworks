from Ternary_Tree.CircAlgorithm.UpUCCSDG import UpUCCSDG
from qiskit_algorithms import VQE
from functools import partial
from scipy.optimize import minimize
from qiskit.primitives import Estimator
from qiskit.quantum_info import Statevector
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit import QuantumCircuit
from qiskit_nature.second_q.circuit.library.ansatzes import UCCSD





from Ternary_Tree.OpenFermionSQ.Ternary_Tree_mapper import TernaryTreeMapper
basis="6-311g"
basis="6-31G"
# basis='STO-3G'
names = 'H 0 0 0; H 0 0 '
active_orbitals = [0,1]
# basis="6-311g"
# names = 'Li 0 0 0; Li 0 0 '
# active_orbitals = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
# active_orbitals = None

def init(dist=0.714):
    global ucc, fermionic_op
    ucc = UpUCCSDG(geometry=names + str(dist), basis=basis, active_orbitals=active_orbitals, num_electrons=(1,1))
    fermionic_op = ucc.mol.hamiltonian.second_q_op()

def jw_mapping():
    mapper = JordanWignerMapper()
    qubit_jw_op = mapper.map(fermionic_op)
    return qubit_jw_op

def tt_mapping():
    # print(ucc.tt)
    fermionic_op = ucc.mol.hamiltonian.second_q_op()
    mapper = TernaryTreeMapper(ucc.tt)
    qubit_tt_op = mapper.map(fermionic_op)
    print(qubit_tt_op)
    return qubit_tt_op

def vqe():
    qc = ucc.get_parametrized_circuit(1)
    ham = tt_mapping()
    estimator = Estimator()
    optimizer = partial(minimize, method="L-BFGS-B")
    res = VQE(estimator, qc, optimizer)
    energy = res.compute_minimum_eigenvalue(ham)
    return energy.eigenvalue + ucc.mol.nuclear_repulsion_energy

def jw_full_vqe():
    qubit_mapper = JordanWignerMapper()
    ham = jw_mapping()
    qc_init = QuantumCircuit(8)
    qc_init.x(0)
    qc_init.x(4)

    qc = UCCSD(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), qubit_mapper=qubit_mapper)
    qc = qc_init.compose(qc)
    ham = jw_mapping()
    estimator = Estimator()
    optimizer = partial(minimize, method="L-BFGS-B")
    res = VQE(estimator, qc, optimizer)
    # state = [0]*2**qc.num_qubits
    # state[0] = 1
    # state = Statevector(state)
    # state = state.evolve(qc)
    # energy = state.expectation_value(ham)
    energy = res.compute_minimum_eigenvalue(ham)
    return energy.eigenvalue + ucc.mol.nuclear_repulsion_energy
    # return energy + ucc.mol.nuclear_repulsion_energy

    
def evolution():
    qc = ucc.get_parametrized_circuit()
    print(qc)
    ham = tt_mapping()
    qc = qc[:]
    qc = QuantumCircuit.from_instructions(qc)
    par = qc.parameters
    qc = qc.assign_parameters({el: 0 for el in par})
    state = [0]*2**qc.num_qubits
    state[0] = 1
    state = Statevector(state)
    state = state.evolve(qc)
    energy = state.expectation_value(ham)
    print("nulpar = ", energy + ucc.mol.nuclear_repulsion_energy)
    return state
n = 12
bk_lex_depth = [0]*n
bk_depth = [0]*n
bk_lex_cx = [0]*n
bk_cx = [0]*n
num_maj = [0]*n
def test(k=0):
    from numpy.random import random
    from qiskit import transpile
    qc = ucc.get_parametrized_circuit()
    ham = tt_mapping()
    qc = QuantumCircuit.from_instructions(qc)
    # print(qc.draw(output='latex_source'))
    # print(qc.decompose(reps=0))

    par = qc.parameters
    qc = qc.assign_parameters({el: random() for el in par})
    qc = qc.decompose(reps=4)
    circ = transpile(qc,  optimization_level=3)
    # print(circ)
    # print(circ)
    # print(transpile(qc.decompose(reps=2)))
    # num = 0
    # for inst in circ.data:
    #     if len(inst.qubits) == 2:
    #         num += 1
    # print(f"num = {num}")
    bk_lex_depth[k-2] = circ.depth() 
    bk_depth[k-2] = circ.decompose(reps=3).count_ops()
    print(f"depth = ", circ.depth() )
    print(circ.decompose(reps=3).count_ops())
    # state = [0]*2**qc.num_qubits
    # state[0] = 1
    # state = Statevector(state)
    # state = state.evolve(qc)
    # energy = state.expectation_value(ham)
    # print("nulpar = ", energy + ucc.mol.nuclear_repulsion_energy)

import pyscf
def energy_classic(bond_length, at_name, nmodes = None, basis ="6-31G"):

    mol = pyscf.M(
        atom = names + str(bond_length),
        basis = basis)
    mf = mol.HF().run()
    mycc = mf.CCSD().run()
    myci = mf.CISD().run()
    return mycc.e_corr, mf.e_tot, myci.e_tot

map_dict = {"JW": JordanWignerMapper, "TT": TernaryTreeMapper}


from numpy import ndarray, linspace

N = 1
GE4 = [0]*N
GE4_JW = [0]*N
R = linspace(0.4, 2.5, N)
R[0] = 0.714
E_CISD = [0]*N
E_CCSD = [0]*N

E_HF = [0]*N  
k = -1

for k in range(N):
# for k in range(2,n + 2):
    # print(n)
    # active_orbitals = [i for i in range(k*2)]
    init(dist=R[0])
    # GE4_JW[n] = jw_full_vqe()

    # GE4[k] = vqe()

    # print(GE4[n])
    # init(dist=R[n]) 
    evolution()
    # test(k)12
    E_CCSD[k], E_HF[k], E_CISD[k] = energy_classic(R[k], "H", basis=basis)
    # print("E_HF = ", E_HF[n])
print(f"pair_UCC: {GE4[0]}")
print(f"UCC: {GE4_JW[0]}")
print(f"HF: {E_HF[0]}")

print(f"CCSD: {E_CCSD[0]}")
print(f"CISD: {E_CISD[0]}")

# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# fig = plt.figure()
# ax = fig.add_subplot(1,1,1)
# f = N
# k = 0
# # ax.plot(R[k:f],GE4[k:f],label="TernaryTree in UCCSD, n_qubits = 4")
# ax.set_xlabel('Межъядерное расстояние (А)')
# ax.set_ylabel('Энергия (эВ)')
# ax.plot(R[k:f],GE4[k:f], label="TernaryTree in UpUCCSD, n_qubits = 8")
# # ax.plot(R[k:f],GE8[k:f],label="TernaryTree in UCCSD, n_qubits = 8")
# ax.plot(R[k:f],E_HF[k:f], label="HarteeFork")
# ax.grid()
# # ax.plot(R[k:f],E_HF[k:f], label="Hartree Fock (classic initial energy)")
# ax.plot(R[k:f],E_CISD[k:f], label="Classic CCSD")
# ax.legend()
# ax.set_title(label="Mоделирование молекулы H_2")
# plt.savefig("figure.png")




