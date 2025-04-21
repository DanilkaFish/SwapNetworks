from __future__ import annotations
from typing import Tuple 
from numpy.random import shuffle


from qiskit_nature.second_q.mappers import JordanWignerMapper, BravyiKitaevMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCC
from qiskit_nature.second_q.mappers.fermionic_mapper import FermionicMapper

from qiskit import transpile, QuantumCircuit
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit.circuit.parametervector import ParameterVector
from qiskit.circuit.library.standard_gates import IGate, XGate, ZGate, YGate
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.synthesis.evolution import synth_pauli_network_rustiq
from qiskit.providers.fake_provider import Fake20QV1, Fake5QV1, GenericBackendV2

from qiskit_aer.noise import (NoiseModel, QuantumError, kraus_error,
    pauli_error, depolarizing_error, thermal_relaxation_error)
from qiskit_aer import AerSimulator
# from qiskit_aer.primitives import EstimatorV2 as Estimator
from qiskit_aer.primitives import Estimator
from qiskit.primitives import BackendEstimator
from qiskit.transpiler import generate_preset_pass_manager
from qiskit_algorithms import NumPyMinimumEigensolver

from Ternary_Tree.ucc.abstractucc import Molecule
from Ternary_Tree.ucc.upgccsd import UpGCCSD, LadExcImpl
from Ternary_Tree.utils import MajoranaContainer, MajoranaMapper
from Ternary_Tree.qiskit_interface import VQEV1 as VQE
# from Ternary_Tree.qiskit_interface.vqe_qiskit_v2 import MyEstimator as Estimator
# from qiskit_aer.library.save_instructions.save_density_matrix import *
# QuantumCircuit.save_density_matrix = save_density_matrix

noise_dict_qiskit = {"I": IGate(),"X": XGate(), "Y": YGate(), "Z": ZGate()}
import numpy as np
tensors_dict = {"X": np.array([[0,1], [1,0]]),"Y":  np.array([[0,-1j], [1j,0]]), "Z":  np.array([[1,0], [0,-1]]),
                "I":  np.array([[1,0], [0,1]])}

trasnpile_backend = Fake20QV1()
trasnpile_backend = Fake5QV1()
# trasnpile_backend = GenericBackendV2(4)
# Create a pass manager for circuit transpilation
# pass_manager = generate_preset_pass_manager(optimization_level=3, backend=trasnpile_backend)
pass_manager = None
# Set the pre-initialization stage of the pass manager with passes suggested by ffsim

# pass_manager.pre_init = ffsim.qiskit.PRE_INIT
def get_file_name(name, noise, method):
    return name + f"_{noise}" + method +".json"

class SwapCircNames:
    SWAP2XN = ("swap2xn", LadExcImpl.CNOT12())
    SWAPGENYORDAN = ("swap_gen", LadExcImpl.YORDAN())
    SWAPGENSHORT = ("swap_gen", LadExcImpl.SHORT())
        

class Circuits:
    @staticmethod
    def jw():
        return "jw"
    
    @staticmethod
    def bk():
        return "bk"
        
    @staticmethod
    def jw_lex():
        return "jw_lex"
        
    @staticmethod
    def bk_lex():
        return "bk_lex"

    @staticmethod
    def swap_2xn():
        return "swap 2xn"

    @staticmethod
    def swap_gen_yor():
        return "swap gen yor"    

    @staticmethod
    def swap_gen_short():
        return "swap gen short"
    
    
    @staticmethod
    def get_circs_names():
        circs = []
        circs.append(Circuits.jw())
        circs.append(Circuits.bk())
        circs.append(Circuits.jw_lex())
        circs.append(Circuits.bk_lex())
        circs.append(Circuits.swap_2xn())
        circs.append(Circuits.swap_gen_short())
        circs.append(Circuits.swap_gen_yor())
        return circs

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
        return "sd"
        for el in self.al:
            ls1.append(((el[1],), (el[0],)))
        for el in self.be:
            ls1.append(((el[1],), (el[0],)))
        for el in self.do:
            if el[0] < el[2]:
                ls2.append(((el[3],el[2]), (el[1],el[0])))
        return ls2
        # return ls2 + ls1

    def get_swap_circuit(self, name: Tuple[str,str]) -> Tuple[QuantumCircuit, SparsePauliOp]:
        method, double = name
        cirq, mtoq = getattr(self.uccgsd, method)(self.reps, double)
        if pass_manager is not None:
            ansatz = pass_manager.run(cirq)
        else:
            ansatz = transpile(cirq,  basis_gates=self.basis_gates, optimization_level=3)
        print_params(ansatz)
        # eq_alpha_beta(qc)
        return ansatz, ucc_ham(self.fermionic_op, mtoq)
    
    def get_ucc(self, qubit_mapper,  init=False):
        num_electrons = (self.uccgsd.n_alpha, self.uccgsd.n_beta)
        initial_state = None if not init else HartreeFock(self.uccgsd.n_spatial, 
                                                        num_electrons,
                                                        qubit_mapper)
        return UCC(
            self.uccgsd.n_spatial, 
            num_electrons,
            excitations=self.to_qiskit_excitations(), 
            qubit_mapper=qubit_mapper, 
            initial_state=initial_state
            ) 
        
    def get_circ_via_mapping(self, qubit_mapper,  init=True):
        qc = self.get_ucc(qubit_mapper, True)
        for i in range(self.reps-1):
            new_ucc = self.get_ucc(qubit_mapper, init=False)
            theta = ParameterVector("θ" + str(i), new_ucc.num_parameters)
            new_ucc.assign_parameters(theta, inplace=True)
            qc = qc.compose(new_ucc)
        if pass_manager is not None:
            ansatz = pass_manager.run(qc)
        else:
            ansatz = transpile(qc,  basis_gates=self.basis_gates, optimization_level=3).decompose(reps=3)
        # ansatz.save_density_matrix()
        # qc = QuantumCircuit(1)
        # qc = QuantumCircuit.from_instructions([*ansatz])
        # qc.save_density_matrix()
        # print(qc)
        print_params(ansatz)
        return ansatz
        
    # TODO reps
    def get_rust_circ(self, qubit_mapper: FermionicMapper):
        circ = HartreeFock(self.uccgsd.n_spatial, (self.uccgsd.n_alpha, self.uccgsd.n_beta), qubit_mapper)
        circ._build()
        
        qc = self.get_ucc(qubit_mapper, init=None)
        nq = qc.num_qubits
        rq = list(range(nq))
        parr = []
        for gate in qc.decompose(reps=1):
            i = 0
            for pauli, coef in zip(gate.operation.operator.paulis, gate.operation.operator.coeffs):
                parr.append((str(pauli), list(reversed(rq)), coef*gate.params[0]))
        print(len(parr))
        circ.compose(synth_pauli_network_rustiq(nq, 
                                                    parr, 
                                                    optimize_count=True, 
                                                    preserve_order=True, 
                                                    upto_phase=True, 
                                                    upto_clifford=False, 
                                                    resynth_clifford_method=0
                                                ),
                                                inplace=True)
        if pass_manager is not None:
            ansatz = pass_manager.run(qc)
        else:
            ansatz = transpile(circ.decompose(reps=3), basis_gates=self.basis_gates, optimization_level=3).decompose(reps=3)
        print_params(ansatz)
        return circ

    def get_circ_op(self, name):
        if name == Circuits.jw():
            return (self.get_circ_via_mapping(JordanWignerMapper()), jw_ham(self.fermionic_op))
        elif name == Circuits.bk():
            return (self.get_circ_via_mapping(BravyiKitaevMapper()), bk_ham(self.fermionic_op))
        elif name == Circuits.swap_2xn():
            return self.get_swap_circuit(SwapCircNames.SWAP2XN)
        elif name == Circuits.swap_gen_yor():
            return self.get_swap_circuit(SwapCircNames.SWAPGENYORDAN)
        elif name == Circuits.swap_gen_short():
            return self.get_swap_circuit(SwapCircNames.SWAPGENSHORT)
        elif name == Circuits.jw_lex():
            return (self.get_rust_circ(JordanWignerMapper()), jw_ham(self.fermionic_op))
        elif name == Circuits.bk_lex():
            return (self.get_rust_circ(BravyiKitaevMapper()), bk_ham(self.fermionic_op))
            
    def get_circ(self, name) -> Tuple[str, QuantumCircuit, SparsePauliOp]:
        return name, *self.get_circ_op(name)

    def __iter__(self, name_list=Circuits.get_circs_names()):
        for name in name_list:
            yield name, *self.get_circ_op(name)


class CircSim:
    def __init__(self, circ, op, noise_par=0.999, noise_type="D",  init_point=None):
        self.circ = circ
        self.op = op
        
        # self.hf = hf(circ, op)
        self.hf = 0
        if init_point is None:
            self.init_point = [0 for _ in circ.parameters]
        else:
            self.init_point = init_point
        self.noise_par = noise_par
        self.noise_type = noise_type

    def run_qiskit_vqe(self, optimizer, device="CPU"):
        params = []

        if self.noise_type == "":
            est = Estimator(
                run_options={"seed": 170, "shots": None, },
                approximation=True,
                backend_options={"device": device},
            )
        else:
            est = get_qiskit_device_noise_estimator(
                                                    noise_op=self.noise_type, 
                                                    prob=self.noise_par,
                                                    device=device
                                                    )
            
        vqe = VQE(est, self.circ, optimizer=optimizer, initial_point=self.init_point)
        result = vqe.compute_minimum_eigenvalue(operator=self.op)
        print(f"VQE on Aer qasm simulator (with noise): {result.eigenvalue.real:.5f}")
        return result.eigenvalue.real, list(result.optimal_parameters.values())
        
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

def get_qiskit_device_noise_estimator(noise_op,  prob, device):
    noise_model = NoiseModel(basis_gates=["u", "cx"])
    if noise_op=="D":
        error1 = depolarizing_error(1 - prob, 1)
        error2 = depolarizing_error(1 - prob, 2)
    else:
        error1 = QuantumError([(noise_dict_qiskit["I"], prob), (noise_dict_qiskit[noise_op], 1 - prob)])
        # print(error1)
        error2 = error1.tensor(error1)
        op = tensors_dict[noise_op]
        # op = tensors_dict["I"]
        I = tensors_dict["I"]
        error2 = kraus_error([np.sqrt(prob) * np.kron(I, I), np.sqrt(1-prob) * np.kron(op, op)])
    noise_model.add_all_qubit_quantum_error(error2, ['cx'])
    if trasnpile_backend is not None:
        noise_model = NoiseModel.from_backend(trasnpile_backend)
    noisy_estimator = Estimator(
                run_options={"seed": 170, "shots": None, },
                approximation=True,
                backend_options={
                        "noise_model": noise_model,
                        "basis_gates": ["u","cx"],
                        "method": "density_matrix",
                        "device": device,
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


if __name__ == "__main__":
    nq = 4
    rq = list(range(4))
    parr = []
    theta = ParameterVector("θ" + str(0), 1)
    for key, coef in {"YYIZ": 1,"XXIZ": 1,"YYZI": 1,"XXZI": 1,"IZYY": -1,"IZXX": -1,"ZIYY": -1,"ZIXX": -1}.items():
        # circ.append(PauliEvolutionGate(pauli, coef*gate.params[0]), rq)
        parr.append((str(key), rq, theta[0]))
    circ = synth_pauli_network_rustiq(nq, parr, optimize_count=False, preserve_order=False, upto_phase=True, resynth_clifford_method=0)
    print(circ)
    # print(circ)
    # ansatz = transpile(circ.decompose(reps=3), trasnpile_backend, basis_gates=["u","cx"], optimization_level=3).decompose(reps=3)
    # print_params(ansatz)