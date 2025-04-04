from __future__ import annotations
from typing import Tuple 


from qiskit_nature.second_q.mappers import JordanWignerMapper, BravyiKitaevMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCC

from qiskit import transpile, QuantumCircuit
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit.circuit.parametervector import ParameterVector
from qiskit.circuit.library.standard_gates import IGate, XGate, ZGate, YGate

from qiskit_aer.noise import (NoiseModel, QuantumError, ReadoutError,
    pauli_error, depolarizing_error, thermal_relaxation_error)
# from qiskit_aer.primitives import EstimatorV2 as Estimator
from qiskit_aer.primitives import Estimator
from qiskit_algorithms import NumPyMinimumEigensolver


from numpy.random import shuffle

from Ternary_Tree.ucc.abstractucc import Molecule
from Ternary_Tree.ucc.upgccsd import UpGCCSD, LadExcImpl
from Ternary_Tree.utils import MajoranaContainer, MajoranaMapper
from Ternary_Tree.qiskit_interface import VQEV1 as VQE
# from Ternary_Tree.qiskit_interface import VQEV2 as VQE

from Paulihedral_v2.Paulihedral_new.benchmark.mypauli import *
from Paulihedral_v2.Paulihedral_new.parallel_bl import depth_oriented_scheduling, gate_count_oriented_scheduling
from Paulihedral_v2.Paulihedral_new.tools import *
from Paulihedral_v2.Paulihedral_new import synthesis_FT


noise_dict_qiskit = {"I": IGate(),"X":XGate(), "Y": YGate(), "Z": ZGate()}


class SwapCircNames:
    SWAP2XN = ("swap2xn", LadExcImpl.CNOT12())
    SWAPGENYORDAN = ("swap_gen", LadExcImpl.YORDAN())
    SWAPGENSHORT = ("swap_gen", LadExcImpl.SHORT())
        
def circ_order():
    return ["jw", "bk", "swap 2xn", "swap gen yor", "swap gen short"]
    # return ["swap 2xn", "swap gen yor", "swap gen short"]
    # return ["swap_sh", "swap_sh_inv", "jw", "jw_lex", "bk", "bk_lex", "jw_opt", "jw_opt_lexic", "swap_yo", "swap_yo_inv", "swap gen"]

def eq_alpha_beta(qc):
    n = qc.num_qubits
    k = n//2
    theta = ParameterVector("θ", qc.num_parameters - k*(k-1)//2)
    params = [theta[i - k*(k-1)//2] for i in range(k*(k-1), 3*k*(k-1)//2)] +\
            [theta[i] for i in range(k*(k-1)//2)] + [theta[i] for i in range(k*(k-1)//2)] 
    qc.assign_parameters(params, inplace=True)

def print_params(ansatz):
    print(ansatz.count_ops())
    print(ansatz.depth())
    




class CircuitProvider:
    def __init__(self,
                 reps=1,
                 molecule: Molecule=Molecule(),
                 basis_gates=["u", "cx"]
                 ):
        self.reps = reps
        self.num_electrons = molecule.num_electrons
        self.basis_gates = basis_gates

        self.uccgsd = UpGCCSD(molecule=molecule)
        self.al = self.uccgsd.get_alpha_excitations()
        self.be = self.uccgsd.get_beta_excitations()
        self.do = self.uccgsd.get_double_excitations()
        self.fermionic_op = self.uccgsd.mol.hamiltonian.second_q_op()

    def to_qiskit_excitations(self, **kwargs):
        ls1 = []
        ls2 = []
        for el in self.al:
            ls1.append(((el[1],), (el[0],)))
        for el in self.be:
            ls1.append(((el[1],), (el[0],)))
        for el in self.do:
            if el[0] < el[2]:
                ls2.append(((el[3],el[2]), (el[1],el[0])))
        return ls2 + ls1

    def get_swap_circuit(self, name: Tuple[str,str]) -> Tuple[QuantumCircuit, SparsePauliOp]:
        method, double = name
        cirq, mtoq = getattr(self.uccgsd, method)(self.reps, double)
        ansatz = transpile(cirq, basis_gates=self.basis_gates, optimization_level=3)
        print_params(ansatz)
        # eq_alpha_beta(qc)
        return ansatz, ucc_ham(self.fermionic_op, mtoq)
    
    def get_circ_via_mapping(self, qubit_mapper,  init=True):
        num_electrons = (self.uccgsd.n_alpha, self.uccgsd.n_beta)
        initial_state = None if not init else HartreeFock(self.uccgsd.n_spatial, 
                                                        num_electrons,
                                                        qubit_mapper)
        qc = UCC(
            self.uccgsd.n_spatial, 
            num_electrons,
            excitations=self.to_qiskit_excitations, 
            qubit_mapper=qubit_mapper, 
            initial_state=initial_state
            )  
        
        for i in range(self.reps-1):
            new_ucc = self.qiskit_ucc(qubit_mapper, init=False)
            theta = ParameterVector("θ" + str(i), new_ucc.num_parameters)
            new_ucc.assign_parameters(theta, inplace=True)
            qc = qc.compose(new_ucc)
        ansatz= transpile(qc, basis_gates=self.basis_gates, optimization_level=3).decompose(reps=3)
        print_params(ansatz)
        return ansatz


    #TODO
    # def get_jw_opt_ansatz(self):
    #     qubit_mapper=TernaryTreeMapper(self.jw_opt_map)
    #     qc = QuantumCircuit(self.dynamic_circ.num_qubits)
    #     num_electrons = self.num_electrons[0]    
    #     for i in range(self.dynamic_circ.num_qubits):
    #         if (self.jw_opt_map[i][0].num <= 2*num_electrons) or (self.uccgsd.n_spatial*2 + 1<= 
    #                 self.jw_opt_map[i][0].num <= 2*self.uccgsd.n_spatial + 2*num_electrons):
    #             qc.x(i)
    #         else:
    #             qc.id(i)
    #     ansatz = qc.compose(self.get_circ_via_mapping(qubit_mapper, lexic=False, init=False))
    #     ansatz = transpile(ansatz, basis_gates=basis_gates, optimization_level=3)
    #     return ansatz.decompose(reps=1)


    #TODO
    # def get_jw_opt_lexic_ansatz(self):
    #     qubit_mapper = TernaryTreeMapper(self.jw_opt_map)
    #     qc = QuantumCircuit(self.dynamic_circ.num_qubits)
    #     num_electrons = self.num_electrons[0]    
    #     for i in range(self.dynamic_circ.num_qubits):
    #         if (self.jw_opt_map[i][0].num <= 2*num_electrons) or (self.uccgsd.n_spatial*2 + 1<= 
    #                 self.jw_opt_map[i][0].num <= 2*self.uccgsd.n_spatial + 2*num_electrons):
    #             qc.x(i)
    #         else:
    #             qc.id(i)
    #     qc = qc.compose(self.get_circ_via_mapping(qubit_mapper, lexic=True, init=False))
    #     ansatz = transpile(qc, basis_gates=basis_gates, optimization_level=3)
    #     return ansatz

    def get_circ_op(self, name):
        match name:
            case "jw":
                return (self.get_circ_via_mapping(JordanWignerMapper()), jw_ham(self.fermionic_op))
            case "bk":
                return (self.get_circ_via_mapping(BravyiKitaevMapper()), bk_ham(self.fermionic_op))
            case "swap 2xn":
                return self.get_swap_circuit(SwapCircNames.SWAP2XN)
            case "swap gen yor":
                return self.get_swap_circuit(SwapCircNames.SWAPGENYORDAN)
            case "swap gen short":
                return self.get_swap_circuit(SwapCircNames.SWAPGENSHORT)
            # case "jw_lex":
                # return (self.get_jw_lexic(), jw_ham(self.fermionic_op))
            # case "bk_lex":
            #     return (self.get_bk_lexic(), bk_ham(self.fermionic_op))
            # case "jw_opt":
            #     return (self.get_jw_opt_ansatz(), ucc_ham(self.fermionic_op, self.jw_opt_map))
            # case "jw_opt_lexic":
            #     return (self.get_jw_opt_lexic_ansatz(), ucc_ham(self.fermionic_op, self.jw_opt_map))
            

    def get_circ(self, name) -> Tuple[str, QuantumCircuit, SparsePauliOp]:
        return name, *self.get_circ_op(name)

    def __iter__(self, name_list=circ_order()):
        for name in name_list:
            yield name, *self.get_circ_op(name)


class CircSim:
    def __init__(self, circ, op, noise_par=0.999, noise_type="D",  init_point=None):
        self.circ = circ
        self.op = op
        
        self.hf = hf(circ, op)
        if init_point is None:
            self.init_point = [0 for _ in circ.parameters]
        else:
            self.init_point = init_point
        self.noise_par = noise_par
        self.noise_type = noise_type

    def run_qiskit_vqe(self, optimizer, device="CPU"):
        counts, values, params = [], [], []
        def store_intermediate_result(eval_count, parameters, mean, std):
            counts.append(eval_count)
            values.append(mean)
            params.append(parameters)
        if self.noise_type == "":
            est = Estimator(
                run_options={"seed": 170, "shots": None, },
                approximation=True,
                backend_options={"device": device},
            )
            vqe = VQE(est, self.circ, optimizer=optimizer, callback=store_intermediate_result, initial_point=self.init_point)
            result = vqe.compute_minimum_eigenvalue(operator=self.op)
            print(f"VQE on Aer qasm simulatfrom qiskit.quantum_info import Statevector (no noise): {result.eigenvalue.real:.5f}")
            return result.eigenvalue.real, params
        else:
            noise_est = get_qiskit_device_noise_estimator(self.circ.num_qubits, 
                                                          noise_op=self.noise_type, 
                                                          prob=self.noise_par)
            
            vqe = VQE(noise_est, self.circ, optimizer=optimizer, callback=store_intermediate_result, initial_point=self.init_point)
            result = vqe.compute_minimum_eigenvalue(operator=self.op)
            print(f"VQE on Aer qasm simulator (with noise): {result.eigenvalue.real:.5f}")
            return result.eigenvalue.real, params
        
    def run_qulacs_sim(parameters):
        pass

def jw_ham(fermionic_op):
    mapper = JordanWignerMapper()
    qubit_jw_op = mapper.map(fermionic_op)
    return qubit_jw_op


def bk_ham(fermionic_op):
    mapper = BravyiKitaevMapper()
    qubit_bk_op = mapper.map(fermionic_op)
    return qubit_bk_op


def ucc_ham(fermionic_op, mtoq: MajoranaContainer):
    mapper = MajoranaMapper(mtoq)
    return mapper.map(fermionic_op)

    

def get_qiskit_device_noise_estimator(n_qubits, noise_op,  prob):
    # coupling_map = [(x,y) for x in range(n_qubits) for y in range(n_qubits) if x != y]
    noise_model = NoiseModel(basis_gates=["u", "cx"])
    if noise_op=="D":
        error1 = depolarizing_error(1 - prob, 1)
        error2 = depolarizing_error(1 - prob, 2)
    else:
        error1 = QuantumError([(noise_dict_qiskit["I"], prob), (noise_dict_qiskit[noise_op], 1 - prob)])
        error2 = error1.tensor(error1)
    noise_model.add_all_qubit_quantum_error(error2, ['cx'])

    noisy_estimator = Estimator(
                run_options={"seed": 170, "shots": None, },
                approximation=True,
                backend_options={
                        "noise_model": noise_model,
                        "basis_gates": ["u","cx"],
                        "method": "density_matrix",
                        "device": "GPU",
                        },
            )
    return noisy_estimator


def numpy_energy(fermionic_op, ucc):
    numpy_solver = NumPyMinimumEigensolver()
    result = numpy_solver.compute_minimum_eigenvalue(operator=bk_ham(fermionic_op))
    ref_value = result.eigenvalue.real
    print(f"Reference value: {ref_value + ucc.mol.nuclear_repulsion_energy:.5f}")
    print(f"Reference value: {ref_value:.5f}")

    return ref_value 

def hf(ansatz, op):
    par = ansatz.parameters
    qc = ansatz.assign_parameters({el: 0 for el in par})
    state = [0]*2**qc.num_qubits
    state[0] = 1
    state = Statevector(state)
    state = state.evolve(qc)
    energy = state.expectation_value(op)
    print("nulpar = ", energy)
    return energy








# else:
#             qc = self.qiskit_ucc(qubit_mapper, init=init)
#             if init:
#                 qc = HartreeFock(self.uccgsd.n_spatial, (self.uccgsd.n_alpha, self.uccgsd.n_beta), qubit_mapper)
#                 qc._build()
#             else:
#                 qc = QuantumCircuit(self.uccgsd.n_qubits)
#             for i in range(self.reps):
#                 new_ucc = self.qiskit_ucc(qubit_mapper)
#                 theta = ParameterVector("θ" + str(i), new_ucc.num_parameters)
#                 new_ucc.assign_parameters(theta, inplace=True)
#                 parr = []
#                 cir = new_ucc.decompose(reps=1)
#                 # print(cir[0].decompose())
#                 names = []
#                 for gate in new_ucc.decompose(reps=1):
#                     for pauli in gate.operation.operator.paulis:
#                         parr.append([pauliString(pauli.__str__(), 1.0)])
#                         names.append(pauli.__str__())
#                 shuffle(parr)
#                 a1 = gate_count_oriented_scheduling(parr)