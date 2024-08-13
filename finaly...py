from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.circuit.library.ansatzes import UCC
from qiskit_nature.second_q.algorithms.initial_points import HFInitialPoint
from qiskit.circuit import CircuitInstruction, instruction, Gate
from qiskit.quantum_info import PauliList, Pauli
from qiskit.quantum_info.operators.symplectic import SparsePauliOp
from qiskit import QuantumCircuit
from qiskit_algorithms.optimizers import SLSQP
from qiskit.primitives import Estimator
from qiskit.circuit.library import HamiltonianGate, PauliEvolutionGate
from qiskit_algorithms import TrotterQRTE, TimeEvolutionProblem
from qiskit.primitives import Estimator
from qiskit import transpile
from Paulihedral_v2.Paulihedral_new.benchmark.mypauli import *
from Paulihedral_v2.Paulihedral_new.parallel_bl import depth_oriented_scheduling, gate_count_oriented_scheduling
from Paulihedral_v2.Paulihedral_new.tools import *
from qiskit import transpile
import time
from Paulihedral_v2.Paulihedral_new import synthesis_FT

from qiskit_nature.second_q.circuit.library.ansatzes.utils import generate_fermionic_excitations
from Ternary_Tree.CircAlgorithm.UpUCCSDG import UpUCCSDG
from Ternary_Tree.OpenFermionSQ.Ternary_Tree_mapper import TernaryTreeMapper
from qiskit_nature.second_q.mappers import BravyiKitaevMapper

basis="6-311g"
names = 'Li 0 0 0; Be 0 0 '
n = 15
bk_lex_depth = [0]*n
bk_depth = [0]*n
bk_lex_cx = [0]*n
bk_cx = [0]*n
num_maj = [0]*n

for k in range(8,n + 2):
    # if k >13: 
    active_orbitals = [i for i in range(k*2)]
    # else:
    #     active_orbitals = [i for i in range(6)]
    n = len(active_orbitals)
    # active_orbitals = None
    def initialize_geometry(distance: float = 0.714) -> None:
        """
        Initialize the UpUCCSDG object with the given geometry and basis.
        """
        global ucc, fermionic_hamiltonian
        ucc = UpUCCSDG(
            geometry=f"{names} {distance}",
            basis=basis,
            active_orbitals=active_orbitals,
            multiplicity=0,
            num_electrons=(2, 2)
        )
        fermionic_hamiltonian = ucc.mol.hamiltonian.second_q_op()
    try:
        initialize_geometry()
        al = ucc.get_alpha_excitations()
        be = ucc.get_beta_excitations()
        do = ucc.get_double_excitations()




        # print(al)
        # print(be)
        # print(do)

        toto = ucc.get_parametrized_circuit()
        # print(toto)
        # qubit_mapper = TernaryTreeMapper(ucc.tt)
        # qubit_mapper = BravyiKitaevMapper()
        # qubit_mapper = JordanWignerMapper()
        # new_ucc = UCC(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), excitations=my_gen, qubit_mapper=qubit_mapper)
        # parr = []
        # from numpy.random import shuffle
        # for gate in new_ucc.decompose(reps=2):
        #     parr.append([pauliString(gate.operation.name[-1-n*2:-1], 1.0)])
        # shuffle(parr)
        # ctime = time.time
        # nq = len(parr[0][0])
        # t0 = ctime()
        # a1 = gate_count_oriented_scheduling(parr)
        # qc = synthesis_FT.block_opt_FT(a1)
        from numpy.random import random
        par = toto.parameters
        qc = toto.assign_parameters({el: random() for el in par})
        qc2 = transpile(qc, basis_gates=['u3', 'cx'], optimization_level=3)
        bk_lex_cx[k-2] = qc2.decompose(reps=3).count_ops()["cx"]
        bk_lex_depth[k-2] = qc2.depth()
        # print_qc(qc2)


        # from numpy.random import random
        # par = new_ucc.parameters
        # qc = new_ucc.assign_parameters({el: random() for el in par})
        # circ = transpile(qc,  optimization_level=3, basis_gates=['id', 'rz', 'rx','cx','h','s'])
        # circ = circ.decompose(reps=3)

        # bk_cx[k-2] = circ.decompose(reps=3).count_ops()["cx"]
        # bk_depth[k-2] = circ.depth()
    except Exception:
        print("totot")
print("jw_lex_depth: ", bk_lex_depth)
print("jw_lex_cx: ", bk_lex_cx)
print("jw_depth: ", bk_depth)
print("jw_cx: ", bk_cx)
print("num_maj: ", num_maj)



