from qiskit_nature.second_q.mappers import JordanWignerMapper, BravyiKitaevMapper
from qiskit_nature.second_q.circuit.library.ansatzes import UCC, UCCSD
from qiskit_nature.second_q.circuit.library import HartreeFock
from qiskit.providers.fake_provider import GenericBackendV2
from qiskit.primitives import Estimator

from qiskit import transpile, QuantumCircuit
from qiskit.circuit.parametervector import ParameterVector
from qiskit_aer.noise import NoiseModel
from qiskit_aer.primitives import Estimator as AerEstimator
# from qiskit_algorithms.utils import algorithm_globals
from qiskit_algorithms import VQE, NumPyMinimumEigensolver

from numpy.random import shuffle
from copy import copy
from Ternary_Tree.OpenFermionSQ.Ternary_Tree_mapper import TernaryTreeMapper
from Ternary_Tree.CircAlgorithm.UpUCCSDG import UpUCCSDG
from Paulihedral_v2.Paulihedral_new.benchmark.mypauli import *
from Paulihedral_v2.Paulihedral_new.parallel_bl import depth_oriented_scheduling, gate_count_oriented_scheduling
from Paulihedral_v2.Paulihedral_new.tools import *
from Paulihedral_v2.Paulihedral_new import synthesis_FT
from qiskit.quantum_info import Statevector


from qiskit_aer.noise import (NoiseModel, QuantumError, ReadoutError,
    pauli_error, depolarizing_error, thermal_relaxation_error)

from qiskit.circuit.library.standard_gates import IGate, XGate, ZGate, YGate

seed = 172
# algorithm_globals.random_seed = seed
basis_gates = ["u3", "cx"]
# basis_gates1 = [ "h","rz","cx"]
# basis_gates1 = basis_gates
# basis_gates2 = ["u3", 'cx']

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

def circ_order():
    return ["swap_sh", "swap_sh_inv", "jw", "jw_lex", "bk", "bk_lex", "jw_opt", "jw_opt_lexic", "swap_yo", "swap_yo_inv"]

def get_ucc_ops(self, **kwargs):
    return ucc_parts(self.al, self.be, self.do, **kwargs)


class CircuitProvider:
    def __init__(self,
                 reps=1,
                 geometry=geometry, 
                 basis=basis, 
                 active_orbitals=active_orbitals, 
                 num_electrons=num_electrons,
                 basis_gates=basis_gates,
                 noise_type='D',
                 noise_probs=[0.999]):
        self.reps = reps
        self.ucc = UpUCCSDG(geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons)
        self.num_electrons = num_electrons
        self.fermionic_op = self.ucc.mol.hamiltonian.second_q_op()
        self.dynamic_circ = self.ucc.get_parametrized_circuit(reps)
        self.ucc.get_parametrized_circuit(reps)

        self.dynamic_map = self.ucc.tt
        self.jw_opt_map = self.ucc.get_jw_opt()
        self.jw_zyx_map = None
        self.basis_gates = basis_gates
        self.al = self.ucc.get_alpha_excitations()
        self.be = self.ucc.get_beta_excitations()
        self.do = self.ucc.get_double_excitations()
        self.noise_type=noise_type
        self.noise_probs=noise_probs
        # self.kwargs = {"num_spatial_orbitals"}

    def get_ucc_ops(self, **kwargs):
        return ucc_parts(self.al, self.be, self.do, **kwargs)


    def get_dynamic(self):
        # ansatz = transpile(self.dynamic_circ, basis_gates=basis_gates1, optimization_level=3)
        ansatz = transpile(self.dynamic_circ, basis_gates=self.basis_gates, optimization_level=3)
        print(ansatz.count_ops())
        print(ansatz.depth())
        return ansatz

    def get_yordan_dynamic(self):
        ansatz = transpile(self.ucc.get_parametrized_circuit(self.reps, type=1), basis_gates=self.basis_gates, optimization_level=3)
        print(ansatz.count_ops())
        print(ansatz.depth())
        return ansatz

    def get_zyx_dynamic(self):
        ansatz = transpile(self.ucc.get_parametrized_circuit(self.reps, type=2), basis_gates=self.basis_gates, optimization_level=3)
        self.jw_zyx_map = self.ucc.tt
        print(ansatz.count_ops())
        print(ansatz.depth())
        return ansatz
    
    def get_zyx_yordan_dynamic(self):
        ansatz = transpile(self.ucc.get_parametrized_circuit(self.reps, type=3), basis_gates=self.basis_gates, optimization_level=3)
        self.jw_zyx_map = self.ucc.tt
        print(ansatz.count_ops())
        print(ansatz.depth())
        return ansatz
    
    def get_ucc(self, qubit_mapper, init=False):
        if init:
            qc = UCC(self.ucc.num_spatial_orbitals, 
                 (self.ucc.num_alpha, self.ucc.num_beta), 
                 excitations=self.get_ucc_ops, 
                 qubit_mapper=qubit_mapper, 
                 initial_state=HartreeFock(self.ucc.num_spatial_orbitals, 
                                           (self.ucc.num_alpha, self.ucc.num_beta), 
                                           qubit_mapper))            
        else:
            qc = UCC(self.ucc.num_spatial_orbitals, 
                        (self.ucc.num_alpha, self.ucc.num_beta), 
                        excitations=self.get_ucc_ops, 
                        qubit_mapper=qubit_mapper, 
                        )

        n = qc.num_qubits
        k = n//2
        theta = ParameterVector("θ", qc.num_parameters - k*(k-1)//2)
        params = [theta[i - k*(k-1)//2] for i in range(k*(k-1), 3*k*(k-1)//2)] +\
             [theta[i] for i in range(k*(k-1)//2)] + [theta[i] for i in range(k*(k-1)//2)] 
        qc.assign_parameters(params, inplace=True)
        return qc

    def get_circ_with_mapping(self, qubit_mapper, lexic=False, init=True):
        if not lexic:
            qc = self.get_ucc(qubit_mapper, init=init)
            
            for i in range(self.reps-1):
                new_ucc = self.get_ucc(qubit_mapper, init=False)
                theta = ParameterVector("θ" + str(i), new_ucc.num_parameters)
                new_ucc.assign_parameters(theta, inplace=True)
                qc = qc.compose(new_ucc)
        else:
            qc = self.get_ucc(qubit_mapper,init=init)
            if init:
                qc = HartreeFock(self.ucc.num_spatial_orbitals, (self.ucc.num_alpha, self.ucc.num_beta), qubit_mapper)
                qc._build()
            else:
                qc = QuantumCircuit(self.ucc.n_qubits)
            for i in range(self.reps):
                new_ucc = self.get_ucc(qubit_mapper)
                theta = ParameterVector("θ" + str(i), new_ucc.num_parameters)
                new_ucc.assign_parameters(theta, inplace=True)
                parr = []
                for gate in new_ucc.decompose(reps=2):
                    parr.append([pauliString(gate.operation.name[-1-self.ucc.n_qubits:-1], 1.0)])
                    shuffle(parr)
                    nq = len(parr[0][0])
                    # length = nq//2 # `length' is a hyperparameter, and can be adjusted for best performance
                    a1 = gate_count_oriented_scheduling(parr)
                for pauli in a1:
                    for gate in new_ucc.decompose(reps=2):
                        if gate.operation.name[-1-self.ucc.n_qubits:-1] == str(pauli[0][0]):
                            qc = qc.compose(gate[0])
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

    def get_jw_opt_ansatz(self):
        qubit_mapper=TernaryTreeMapper(self.jw_opt_map)
        qc = QuantumCircuit(self.dynamic_circ.num_qubits)
        num_electrons = self.num_electrons[0]    
        for i in range(self.dynamic_circ.num_qubits):
            if (self.jw_opt_map[i][0].num <= 2*num_electrons) or (self.ucc.num_spatial_orbitals*2 + 1<= 
                    self.jw_opt_map[i][0].num <= 2*self.ucc.num_spatial_orbitals + 2*num_electrons):
                # qc.h(i)
                # qc.s(i)
                # qc.z(i)
                qc.x(i)
            else:
                qc.id(i)
                # qc.h(i)
                # qc.s(i)
        ansatz = qc.compose(self.get_circ_with_mapping(qubit_mapper, lexic=False, init=False))
        ansatz = transpile(ansatz, basis_gates=basis_gates, optimization_level=3)
        return ansatz.decompose(reps=1)


    def get_jw_opt_lexic_ansatz(self):
        qubit_mapper = TernaryTreeMapper(self.jw_opt_map)
        qc = QuantumCircuit(self.dynamic_circ.num_qubits)
        num_electrons = self.num_electrons[0]    
        for i in range(self.dynamic_circ.num_qubits):
            if (self.jw_opt_map[i][0].num <= 2*num_electrons) or (self.ucc.num_spatial_orbitals*2 + 1<= 
                    self.jw_opt_map[i][0].num <= 2*self.ucc.num_spatial_orbitals + 2*num_electrons):
                qc.x(i)
            else:
                qc.id(i)
        qc = qc.compose(self.get_circ_with_mapping(qubit_mapper, lexic=True, init=False))
        ansatz = transpile(qc, basis_gates=basis_gates, optimization_level=3)
        return ansatz

    def get_circ_op(self, name):
        match name:
            case "swap_sh":
                return (self.get_dynamic(), ucc_ham(self.fermionic_op, self.dynamic_map))
            case "swap_sh_inv":
                return (self.get_zyx_dynamic(), ucc_ham(self.fermionic_op, self.jw_zyx_map))
            case "jw":
                return (self.get_jw(), jw_ham(self.fermionic_op))
            case "jw_lex":
                return (self.get_jw_lexic(), jw_ham(self.fermionic_op))
            case "bk":
                return (self.get_bk(), bk_ham(self.fermionic_op))
            case "bk_lex":
                return (self.get_bk_lexic(), bk_ham(self.fermionic_op))
            case "jw_opt":
                return (self.get_jw_opt_ansatz(), ucc_ham(self.fermionic_op, self.jw_opt_map))
            case "jw_opt_lexic":
                return (self.get_jw_opt_lexic_ansatz(), ucc_ham(self.fermionic_op, self.jw_opt_map))
            case "swap_yo":
                return (self.get_yordan_dynamic(), ucc_ham(self.fermionic_op, self.dynamic_map))
            case "swap_yo_inv":
                return (self.get_zyx_yordan_dynamic(), ucc_ham(self.fermionic_op, self.jw_zyx_map))
    
    def __iter__(self, name_list=circ_order()):
        for name in name_list:
            yield name, *self.get_circ_op(name)

class CircSim:
    def __init__(self, circ, op, is_noise=False, noise_par=0.999, noise_type="D", q1_noise=False, q2_noise=True):
        self.circ = circ
        self.op = op
        self.hf = hf(circ, op)
        self.init_point = [0 for _ in circ.parameters]
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
            # print(parameters)
        # optimizer = SLSQP(maxiter=200, ftol=0)
        if not self.is_noise:
            est = AerEstimator(approximation=True, 
                            run_options={"shots": None})
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
            return result.eigenvalue.real, counts, values, params
        
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


def ucc_ham(fermionic_op, tree_map):
    mapper = TernaryTreeMapper(tree_map)
    qubit_tt_op = mapper.map(fermionic_op)
    return qubit_tt_op

    
noise_dict_qiskit = {"X":XGate(), "Y": YGate(), "Z": ZGate()}

def get_qiskit_device_noise_estimator(n_qubits, noise_op, q1, q2, prob):
    # coupling_map = [(x,y) for x in range(n_qubits) for y in range(n_qubits) if x != y]
    noise_model = NoiseModel(basis_gates=basis_gates)
    if noise_op=="D":
        error1 = depolarizing_error(1 - prob, 1)
        error2 = depolarizing_error(1 - prob, 2)
    else:
        error1 = QuantumError([(noise_dict_qiskit("I"), prob), (noise_dict_qiskit(noise_op), 1 - prob)])
        error2 = QuantumError([(noise_dict_qiskit("I"), prob), (noise_dict_qiskit(noise_op), 1 - prob)])
    if q1:
        for i in range(n_qubits):
            noise_model.add_quantum_error(error1, ['u3'], [i])
    if q2:
        noise_model.add_all_qubit_quantum_error(error2, ['cx'])


    noisy_estimator = AerEstimator(
        backend_options={
            "method": "density_matrix",
            # "coupling_map": coupling_map,
            "noise_model": noise_model,
        },
        run_options={"shots": None},
        # transpile_options={"seed_transpiler": seed},
        approximation=True
    )
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
    est = AerEstimator(approximation=True, 
                       run_options={"shots": None})

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