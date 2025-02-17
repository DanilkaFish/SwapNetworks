from __future__ import annotations
from typing import Tuple 


from qiskit_nature.second_q.mappers import JordanWignerMapper, BravyiKitaevMapper
from qiskit_nature.second_q.circuit.library.ansatzes import UCC, UCCSD
from qiskit_nature.second_q.circuit.library import HartreeFock
from qiskit.providers.fake_provider import GenericBackendV2
from qiskit.primitives import Estimator

from qiskit import transpile, QuantumCircuit
from qiskit.circuit.parametervector import ParameterVector
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel
from qiskit_aer.primitives import EstimatorV2 as Estimator
# from qiskit_algorithms.utils import algorithm_globals
from qiskit_algorithms import VQE, NumPyMinimumEigensolver

from numpy.random import shuffle
# from Ternary_Tree.UCC.UpGCCSD import UpUCCSDG
from Ternary_Tree.UCC.AbstractUCC import Molecule
from Ternary_Tree.UCC.UpGCCSDopt import UpUCCSDG, LadExcNames
from Paulihedral_v2.Paulihedral_new.benchmark.mypauli import *
from Paulihedral_v2.Paulihedral_new.parallel_bl import depth_oriented_scheduling, gate_count_oriented_scheduling
from Paulihedral_v2.Paulihedral_new.tools import *
from Paulihedral_v2.Paulihedral_new import synthesis_FT
from qiskit.quantum_info import Statevector


from qiskit_aer.noise import (NoiseModel, QuantumError, ReadoutError,
    pauli_error, depolarizing_error, thermal_relaxation_error)

from qiskit.circuit.library.standard_gates import IGate, XGate, ZGate, YGate

seed = 172
basis_gates = ["u", "cx"]

basis = '6-31G'
geometry='H 0 0 0; H 0 0 0.739'
active_orbitals=[0, 1]
num_electrons=(1, 1)

def ucc_parts(al,be,do, **kwargs):
    ls1 = []
    ls2 = []
    for el in al:
        ls1.append(((el[1],), (el[0],)))
    for el in be:
        ls1.append(((el[1],), (el[0],)))
    for el in do:
        if el[0] < el[2]:
            ls2.append(((el[3],el[2]), (el[1],el[0])))
    return ls2 + ls1

class SwapCircNames:
    SWAP2XN = ("swap2xn", LadExcNames.CNOT12())
    SWAPGENYORDAN = ("swap_gen", LadExcNames.YORDAN())
    SWAPGENSHORT = ("swap_gen", LadExcNames.SHORT())
        
def circ_order():
    return ["jw", "bk",   "swap 2xn", "swap gen yor", "swap gen short"]
    # return ["swap_sh", "swap_sh_inv", "jw", "jw_lex", "bk", "bk_lex", "jw_opt", "jw_opt_lexic", "swap_yo", "swap_yo_inv", "swap gen"]

def get_ucc_ops(self, **kwargs):
    return ucc_parts(self.al, self.be, self.do, **kwargs)


class CircuitProvider:
    def __init__(self,
                 reps=1,
                 molecule: Molecule=Molecule(),
                 basis_gates=basis_gates
                 ):
        self.reps = reps
        self.num_electrons = num_electrons
        self.basis_gates = basis_gates

        self.ucc = UpUCCSDG(molecule=molecule)
        self.al = self.ucc.get_alpha_excitations()
        self.be = self.ucc.get_beta_excitations()
        self.do = self.ucc.get_double_excitations()
        self.fermionic_op = self.ucc.mol.hamiltonian.second_q_op()
        # self.swap2xn = self.ucc_opt.get_parametrized_circuit(reps)
        # self.ucc_opt_tt = self.ucc_opt.mtoq
        # self.swapgen_yordan = self.ucc_opt.get_paramaetrized_circuit_generalized(reps, type=1)
        # self.swapgen_short = self.ucc_opt.get_paramaetrized_circuit_generalized(reps, type=0)
        # self.ucc_gen_tt = self.ucc_opt.mtoq
        # self.ucc.swap2xn(reps)

        # self.dynamic_map = self.ucc.tt
        # self.jw_opt_map = self.ucc.get_jw_opt()
        # self.jw_zyx_map = None
        # self.noise_type=noise_type
        # self.noise_probs=noise_probs
        # self.kwargs = {"num_spatial_orbitals"}

    def get_ucc_ops(self, **kwargs):
        return ucc_parts(self.al, self.be, self.do, **kwargs)

    def swap_circuit(self, name: Tuple[str,str]) -> Tuple[QuantumCircuit, None]:
        method, double = name
        cirq, mtoq = getattr(self.ucc, method)(self.reps, double)
        ansatz = transpile(cirq, basis_gates=self.basis_gates, optimization_level=3)
        print(ansatz.count_ops())
        print(ansatz.depth())
        return ansatz, mtoq
    
    def qiskit_ucc(self, qubit_mapper, init=False):
        if init:
            qc = UCC(self.ucc.n_spatial, 
                 (self.ucc.n_alpha, self.ucc.n_beta), 
                 excitations=self.get_ucc_ops, 
                 qubit_mapper=qubit_mapper, 
                 initial_state=HartreeFock(self.ucc.n_spatial, 
                                           (self.ucc.n_alpha, self.ucc.n_beta), 
                                           qubit_mapper))            
        else:
            qc = UCC(self.ucc.n_spatial, 
                        (self.ucc.n_alpha, self.ucc.n_beta), 
                        excitations=self.get_ucc_ops, 
                        qubit_mapper=qubit_mapper, 
                        )

        # n = qc.num_qubits
        # k = n//2
        # theta = ParameterVector("θ", qc.num_parameters - k*(k-1)//2)
        # params = [theta[i - k*(k-1)//2] for i in range(k*(k-1), 3*k*(k-1)//2)] +\
        #      [theta[i] for i in range(k*(k-1)//2)] + [theta[i] for i in range(k*(k-1)//2)] 
        # qc.assign_parameters(params, inplace=True)
        return qc

    def get_circ_with_mapping(self, qubit_mapper, lexic=False, init=True):
        if not lexic:
            qc = self.qiskit_ucc(qubit_mapper, init=init)
            
            for i in range(self.reps-1):
                new_ucc = self.qiskit_ucc(qubit_mapper, init=False)
                theta = ParameterVector("θ" + str(i), new_ucc.num_parameters)
                new_ucc.assign_parameters(theta, inplace=True)
                qc = qc.compose(new_ucc)
        else:
            qc = self.qiskit_ucc(qubit_mapper,init=init)
            if init:
                qc = HartreeFock(self.ucc.n_spatial, (self.ucc.n_alpha, self.ucc.n_beta), qubit_mapper)
                qc._build()
            else:
                qc = QuantumCircuit(self.ucc.n_qubits)
            for i in range(self.reps):
                new_ucc = self.qiskit_ucc(qubit_mapper)
                theta = ParameterVector("θ" + str(i), new_ucc.num_parameters)
                new_ucc.assign_parameters(theta, inplace=True)
                parr = []
                print(new_ucc)
                cir = new_ucc.decompose(reps=1)
                print(cir[0].operation.operator.paulis)
                # print(cir[0].decompose())
                names = []
                for gate in new_ucc.decompose(reps=1):
                    for pauli in gate.operation.operator.paulis:
                        parr.append([pauliString(pauli.__str__(), 1.0)])
                        names.append(pauli.__str__())
                print(names)
                shuffle(parr)
                a1 = gate_count_oriented_scheduling(parr)
                print(a1)
        qc = transpile(qc, basis_gates=self.basis_gates, optimization_level=3)
        ansatz = qc.decompose(reps=3)
        print(ansatz.count_ops())
        print(ansatz.depth())
        return ansatz


    def get_jw(self):
        qubit_mapper=JordanWignerMapper()
        circ = self.get_circ_with_mapping(qubit_mapper)
        return circ
    
    def get_bk(self):
        qubit_mapper=BravyiKitaevMapper()
        return self.get_circ_with_mapping(qubit_mapper)

    def get_jw_lexic(self):
        qubit_mapper=JordanWignerMapper()
        return self.get_circ_with_mapping(qubit_mapper, lexic=True)

    def get_bk_lexic(self):
        qubit_mapper=BravyiKitaevMapper()
        return self.get_circ_with_mapping(qubit_mapper, lexic=True)

    #TODO
    # def get_jw_opt_ansatz(self):
    #     qubit_mapper=TernaryTreeMapper(self.jw_opt_map)
    #     qc = QuantumCircuit(self.dynamic_circ.num_qubits)
    #     num_electrons = self.num_electrons[0]    
    #     for i in range(self.dynamic_circ.num_qubits):
    #         if (self.jw_opt_map[i][0].num <= 2*num_electrons) or (self.ucc.n_spatial*2 + 1<= 
    #                 self.jw_opt_map[i][0].num <= 2*self.ucc.n_spatial + 2*num_electrons):
    #             qc.x(i)
    #         else:
    #             qc.id(i)
    #     ansatz = qc.compose(self.get_circ_with_mapping(qubit_mapper, lexic=False, init=False))
    #     ansatz = transpile(ansatz, basis_gates=basis_gates, optimization_level=3)
    #     return ansatz.decompose(reps=1)


    #TODO
    # def get_jw_opt_lexic_ansatz(self):
    #     qubit_mapper = TernaryTreeMapper(self.jw_opt_map)
    #     qc = QuantumCircuit(self.dynamic_circ.num_qubits)
    #     num_electrons = self.num_electrons[0]    
    #     for i in range(self.dynamic_circ.num_qubits):
    #         if (self.jw_opt_map[i][0].num <= 2*num_electrons) or (self.ucc.n_spatial*2 + 1<= 
    #                 self.jw_opt_map[i][0].num <= 2*self.ucc.n_spatial + 2*num_electrons):
    #             qc.x(i)
    #         else:
    #             qc.id(i)
    #     qc = qc.compose(self.get_circ_with_mapping(qubit_mapper, lexic=True, init=False))
    #     ansatz = transpile(qc, basis_gates=basis_gates, optimization_level=3)
    #     return ansatz

    def get_circ_op(self, name):
        match name:
            case "jw":
                return (self.get_jw(), jw_ham(self.fermionic_op))
            case "jw_lex":
                return (self.get_jw_lexic(), jw_ham(self.fermionic_op))
            case "bk":
                return (self.get_bk(), bk_ham(self.fermionic_op))
            case "bk_lex":
                return (self.get_bk_lexic(), bk_ham(self.fermionic_op))
            # case "jw_opt":
            #     return (self.get_jw_opt_ansatz(), ucc_ham(self.fermionic_op, self.jw_opt_map))
            # case "jw_opt_lexic":
            #     return (self.get_jw_opt_lexic_ansatz(), ucc_ham(self.fermionic_op, self.jw_opt_map))
            case "swap 2xn":
                return self.swap_circuit(SwapCircNames.SWAP2XN)
            case "swap gen yor":
                return self.swap_circuit(SwapCircNames.SWAP2XN)
            case "swap gen short":
                return self.swap_circuit(SwapCircNames.SWAP2XN)

    def __iter__(self, name_list=circ_order()):
        for name in name_list:
            yield name, *self.get_circ_op(name)


class CircSim:
    def __init__(self, circ, op, is_noise=False, noise_par=0.999, noise_type="D", q1_noise=False, q2_noise=True, init_point=None):
        self.circ = circ
        self.op = op
        
        self.hf = hf(circ, op)
        if init_point is None:
            self.init_point = [0 for _ in circ.parameters]
        else:
            self.init_point = init_point
        self.is_noise = is_noise
        self.noise_par = noise_par
        self.noise_type = noise_type
        self.q1_noise = q1_noise
        self.q2_noise = q2_noise

    def run_qiskit_vqe(self, optimizer):
        counts, values, params = [], [], []
        def store_intermediate_result(eval_count, parameters, mean, std):
            counts.append(eval_count)
            values.append(mean)
            params.append(parameters)
        if not self.is_noise:
            est = Estimator()
            est.options.run_options={"shots": None}
            # print(self.circ)
            vqe = VQE(est, self.circ, optimizer=optimizer, callback=store_intermediate_result, initial_point=self.init_point)
            result = vqe.compute_minimum_eigenvalue(operator=self.op)
            print(f"VQE on Aer qasm simulatfrom qiskit.quantum_info import Statevector (no noise): {result.eigenvalue.real:.5f}")
            print(f"Delta from hf is {(result.eigenvalue.real - self.hf):.5f}")
            return result.eigenvalue.real, counts, values, params
        else:
            noise_est = get_qiskit_device_noise_estimator(self.circ.num_qubits, 
                                                          noise_op=self.noise_type, 
                                                          q1=self.q1_noise, 
                                                          q2=self.q2_noise,
                                                          prob=self.noise_par)
            
            vqe = VQE(noise_est, self.circ, optimizer=optimizer, callback=store_intermediate_result, initial_point=self.init_point)
            result = vqe.compute_minimum_eigenvalue(operator=self.op)
            print(f"VQE on Aer qasm simulator (with noise): {result.eigenvalue.real:.5f}")
            print(f"Delta from hf is {(result.eigenvalue.real - self.hf):.5f}")
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


# def ucc_ham(fermionic_op, tree_map):
#     mapper = TernaryTreeMapper(tree_map)
#     qubit_tt_op = mapper.map(fermionic_op)
#     return qubit_tt_op

    
noise_dict_qiskit = {"I": IGate(),"X":XGate(), "Y": YGate(), "Z": ZGate()}

def get_qiskit_device_noise_estimator(n_qubits, noise_op, q1, q2, prob):
    # coupling_map = [(x,y) for x in range(n_qubits) for y in range(n_qubits) if x != y]
    noise_model = NoiseModel(basis_gates=basis_gates)
    if noise_op=="D":
        error1 = depolarizing_error(1 - prob, 1)
        error2 = depolarizing_error(1 - prob, 2)
    else:
        error1 = QuantumError([(noise_dict_qiskit["I"], prob), (noise_dict_qiskit[noise_op], 1 - prob)])
        error2 = error1.tensor(error1)
    if q1:
        for i in range(n_qubits):
            noise_model.add_quantum_error(error1, ['u'], [i])
    if q2:
        # noise_model.add_quantum_error(error1, ['u'], [i])
        noise_model.add_all_qubit_quantum_error(error2, ['cx'])


    # noisy_estimator = AerSimulator(
    #         noise_model=noise_model,
    #         basis_gates=["u","cx"],
    #         # backend_options = {
    #         method="density_matrix",
    #         shots=None
    #         # "coupling_map": coupling_map,
    #         # "noise_model": noise_model,
    #         # }
    # )
    noisy_estimator = Estimator(options={"backend_options":{
            "noise_model": noise_model,
            "basis_gates": ["u","cx"],
            "method": "density_matrix",
            }
            }
            )
    # noisy_estimator.options.backend_options = {
    #         # "method": "density_matrix",
    #         # "coupling_map": coupling_map,
    #         "noise_model": noise_model,
    # }
    
    noisy_estimator.options.run_options={"shots": None,}
        # transpile_options={"seed_transpiler": seed},
    # noisy_estimator.options.
    print(noisy_estimator.options)
    return noisy_estimator


def numpy_energy(fermionic_op, ucc):
    numpy_solver = NumPyMinimumEigensolver()
    result = numpy_solver.compute_minimum_eigenvalue(operator=bk_ham(fermionic_op))
    ref_value = result.eigenvalue.real
    print(f"Reference value: {ref_value + ucc.mol.nuclear_repulsion_energy:.5f}")
    print(f"Reference value: {ref_value:.5f}")

    return ref_value 


def ideal_energy(ansatz, op, ref_value, optimizer, init_point, en_hf=None):
    counts, values, params = [], [], []
    def store_intermediate_result(eval_count, parameters, mean, std):
        counts.append(eval_count)
        values.append(mean)
        params.append(parameters)
        # print(parameters)
    est = Estimator(approximation=True, 
                       run_options={"shots": None})
    est.options.run_options={"shots": None}
    vqe = VQE(est, ansatz, optimizer=optimizer, callback=store_intermediate_result, initial_point=init_point)
    result = vqe.compute_minimum_eigenvalue(operator=op)
    print(f"VQE on Aer qasm simulatfrom qiskit.quantum_info import Statevector (no noise): {result.eigenvalue.real:.5f}")
    print(f"Delta init from reference energy value is {(values[0] - ref_value):.5f}")
    print(f"Delta from hf is {(result.eigenvalue.real - en_hf):.5f}")
    return result.eigenvalue.real, counts, values, params


# def VQE_energy_with_noise(ansatz, op, ref_value, optimizer, init_point, en_hf, noise_est=None, prob=0.999):
#     # noiseless_estimator = Estimator()
#     counts, values, params = [], [], []
#     def store_intermediate_result(eval_count, parameters, mean, std):
#         counts.append(eval_count)
#         values.append(mean)
#         params.append(parameters)
#     if noise_est is None:
#         noise_est = get_qiskit_device_noise_estimator(ansatz.num_qubits, prob=prob)
#     vqe = VQE(noise_est, ansatz, optimizer=optimizer, callback=store_intermediate_result, initial_point=init_point)
#     result = vqe.compute_minimum_eigenvalue(operator=op)

#     print(f"VQE on Aer qasm simulator (with noise): {result.eigenvalue.real:.5f}")
#     # print(f"Delta from reference energy value is {(result.eigenvalue.real - ref_value):.5f}")
#     # print(f"Delta init from reference energy value is {(values[0] - ref_value):.5f}")
#     # print(f"Delta from hf is {(result.eigenvalue.real - en_hf):.5f}")

#     return result.eigenvalue.real, counts, values, params


# def get_circs_and_ops(get_cp=False, reps=1,**kwargs):
#     circ_prov = CircuitProvider(reps=reps, **kwargs)
#     thebest_qc = circ_prov.get_dynamic()
#     thebest_qc_inversed = circ_prov.get_zyx_dynamic()
    
#     # init_point[0] =  0.07
#     jw_qc = circ_prov.get_jw()
#     from numpy.random import uniform
#     init_point = [0 for _ in jw_qc.parameters]
#     # init_point = [uniform(-0.1, 0.1, 1) for _ in jw_qc.parameters]
#     # init_point = None
    # jw_lexic_qc = circ_prov.get_jw_lexic()
    # bk_qc = circ_prov.get_bk()
    # bk_lexic_qc = circ_prov.get_bk_lexic()
    # # jw_qc = None
    # # jw_lexic_qc = None
    # # bk_qc = None
    # # bk_lexic_qc = None
    # jw_opt_qc = circ_prov.get_jw_opt_ansatz()
    # jw_opt_lexic_qc = circ_prov.get_jw_opt_lexic_ansatz()
    # # jw_opt_qc = None
    # # jw_opt_lexic_qc = None
    # print(circ_prov.jw_zyx_map)
    # en_hf = hf(thebest_qc, ucc_ham(circ_prov.fermionic_op, circ_prov.dynamic_map))
    # circs = [thebest_qc, thebest_qc_inversed, jw_qc, jw_lexic_qc, bk_qc, bk_lexic_qc, jw_opt_qc, jw_opt_lexic_qc]
    # ops = [ucc_ham(circ_prov.fermionic_op, circ_prov.dynamic_map),ucc_ham(circ_prov.fermionic_op, circ_prov.jw_zyx_map), jw_ham(circ_prov.fermionic_op), jw_ham(circ_prov.fermionic_op), bk_ham(circ_prov.fermionic_op), 
    #     bk_ham(circ_prov.fermionic_op), ucc_ham(circ_prov.fermionic_op, circ_prov.jw_opt_map), ucc_ham(circ_prov.fermionic_op, circ_prov.jw_opt_map)]
    # ops = del_identity(ops)
    # ref_value = numpy_energy(circ_prov.fermionic_op, circ_prov.ucc)
    # if get_cp:
    #     return circs, ops, ref_value, init_point, en_hf, circ_prov
    # else:
    #     return circs, ops, ref_value, init_point, en_hf

# circs, ops, ref_value, init_point,en_hf = get_circs_and_ops()

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