from qiskit_nature.second_q.mappers import JordanWignerMapper, BravyiKitaevMapper
from qiskit_nature.second_q.circuit.library.ansatzes import UCC, UCCSD
from qiskit_nature.second_q.circuit.library import HartreeFock
from qiskit.providers.fake_provider import GenericBackendV2
from qiskit.primitives import Estimator

from qiskit import transpile, QuantumCircuit
from qiskit.circuit.parametervector import ParameterVector
from qiskit_aer.noise import NoiseModel
from qiskit_aer.primitives import Estimator as AerEstimator
from qiskit_algorithms.utils import algorithm_globals
from qiskit_algorithms import VQE, NumPyMinimumEigensolver

from numpy.random import shuffle
from copy import copy
from Ternary_Tree.OpenFermionSQ.Ternary_Tree_mapper import TernaryTreeMapper
from Ternary_Tree.CircAlgorithm.UpUCCSDG import UpUCCSDG
from Paulihedral_v2.Paulihedral_new.benchmark.mypauli import *
from Paulihedral_v2.Paulihedral_new.parallel_bl import depth_oriented_scheduling, gate_count_oriented_scheduling
from Paulihedral_v2.Paulihedral_new.tools import *
from Paulihedral_v2.Paulihedral_new import synthesis_FT


from qiskit_aer.noise import (NoiseModel, QuantumError, ReadoutError,
    pauli_error, depolarizing_error, thermal_relaxation_error)

from qiskit.circuit.library.standard_gates import IGate, XGate, ZGate

seed = 172
algorithm_globals.random_seed = seed
basis_gates = ["u3", "cx"]
basis = '6-31G'
geometry='H 0 0 0; H 0 0 0.739'
active_orbitals=[0, 1]
num_electrons=(1, 1)


def my_generation(al,be,do, **kwargs):
    ls1 = []
    ls2 = []
    for el in al:
        ls1.append(((el[1],), (el[0],)))
    for el in be:
        ls1.append(((el[1],), (el[0],)))
    for el in do:
        if el[0] < el[2]:
            ls2.append(((el[3],el[2]), (el[1],el[0])))
    return ls2
    return ls2 + ls1

def my_gen(self, **kwargs):
    return my_generation(self.al, self.be, self.do, **kwargs)

class CircuitProvider:
    def __init__(self,
                 reps=1,
                 geometry=geometry, 
                 basis=basis, 
                 active_orbitals=active_orbitals, 
                 num_electrons=num_electrons,
                 basis_gates=basis_gates):
        self.reps = reps
        self.ucc = UpUCCSDG(geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons)
        self.num_electrons = num_electrons
        self.fermionic_op = self.ucc.mol.hamiltonian.second_q_op()
        self.dynamic_circ = self.ucc.get_parametrized_circuit(reps)
        self.dynamic_map = self.ucc.tt
        self.jw_opt_map = self.ucc.get_jw_opt()
        self.basis_gates = basis_gates
        self.al = self.ucc.get_alpha_excitations()
        self.be = self.ucc.get_beta_excitations()
        self.do = self.ucc.get_double_excitations()
        # self.kwargs = {"num_spatial_orbitals"}

    def my_gen(self, **kwargs):
        return my_generation(self.al, self.be, self.do, **kwargs)


    def get_dynamic(self):
        ansatz = transpile(self.dynamic_circ, basis_gates=self.basis_gates, optimization_level=3)
        print(ansatz.count_ops())
        return ansatz


    def get_ucc(self, qubit_mapper, init=False):
        # if init:
        #     return UCCSD(self.ucc.num_spatial_orbitals, 
        #          (self.ucc.num_alpha, self.ucc.num_beta), 
        #          qubit_mapper=qubit_mapper, 
        #          initial_state=HartreeFock(self.ucc.num_spatial_orbitals, 
        #                                    (self.ucc.num_alpha, self.ucc.num_beta), 
        #                                    qubit_mapper))
        # else:
        #     return UCCSD(self.ucc.num_spatial_orbitals, 
        #                 (self.ucc.num_alpha, self.ucc.num_beta), 
        #                 qubit_mapper=qubit_mapper, 
        #                 )
        if init:
            qc = UCC(self.ucc.num_spatial_orbitals, 
                 (self.ucc.num_alpha, self.ucc.num_beta), 
                 excitations=self.my_gen, 
                 qubit_mapper=qubit_mapper, 
                 initial_state=HartreeFock(self.ucc.num_spatial_orbitals, 
                                           (self.ucc.num_alpha, self.ucc.num_beta), 
                                           qubit_mapper))            
        else:
            qc = UCC(self.ucc.num_spatial_orbitals, 
                        (self.ucc.num_alpha, self.ucc.num_beta), 
                        excitations=self.my_gen, 
                        qubit_mapper=qubit_mapper, 
                        )

        n = qc.num_qubits
        k = n//2
        # theta = ParameterVector("θ", qc.num_parameters - k*(k-1)//2)
        # params = [theta[i - k*(k-1)//2] for i in range(k*(k-1), 3*k*(k-1)//2)] +\
        #      [theta[i] for i in range(k*(k-1)//2)] + [theta[i] for i in range(k*(k-1)//2)] 
        # qc.assign_parameters(params, inplace=True)
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
        ansatz = transpile(qc, basis_gates=basis_gates, optimization_level=3)
        # ansatz = qc
        print(ansatz.count_ops())
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
                qc.x(i)
            else:
                qc.id(i)
        qc = qc.compose(self.get_circ_with_mapping(qubit_mapper, lexic=False, init=False))
        return qc

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
        return qc

def jw_ham(fermionic_op):
    mapper = JordanWignerMapper()
    qubit_jw_op = mapper.map(fermionic_op)
    return qubit_jw_op


def bk_ham(fermionic_op):
    mapper = BravyiKitaevMapper()
    qubit_jw_op = mapper.map(fermionic_op)
    return qubit_jw_op


def ucc_ham(fermionic_op, tree_map):
    mapper = TernaryTreeMapper(tree_map)
    qubit_tt_op = mapper.map(fermionic_op)
    return qubit_tt_op

def create_noise(n_qubits, coupling_map, prob=0.9999):
    noise_model = NoiseModel(basis_gates=basis_gates)

    error1 = depolarizing_error(0.0001, 1)
    noise_ops = [(IGate(), prob), (XGate(), 1 - prob)]
    error1 = QuantumError(noise_ops)
    # for i in range(n_qubits):
    #     noise_model.add_quantum_error(error1, ['u3'], [i])
    # error2 = depolarizing_error(0.0001, 2)
    # for x,y in coupling_map:
    #     noise_model.add_quantum_error(error2, ["cx"], [x,y])
    # noise_model.add_all_qubit_quantum_error(error2, ['cx'])
    noise_model.add_all_qubit_quantum_error(error1, ['u3'])

    # Print noise model info
    return noise_model

def get_device_noise_estimator(n_qubits, prob=0.9999):
    coupling_map = [(x,y) for x in range(n_qubits) for y in range(n_qubits) if x != y]

    # device = GenericBackendV2(num_qubits=n_qubits, seed=50)

    # noise_model = NoiseModel.from_backend(device)
    noise_model = create_noise(n_qubits, coupling_map, prob)
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
    print(noise_model)
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


def VQE_energy_with_noise(ansatz, op, ref_value, optimizer, init_point, en_hf):
    # noiseless_estimator = Estimator()
    counts, values, params = [], [], []
    def store_intermediate_result(eval_count, parameters, mean, std):
        counts.append(eval_count)
        values.append(mean)
        params.append(parameters)
    vqe = VQE(get_device_noise_estimator(ansatz.num_qubits), ansatz, optimizer=optimizer, callback=store_intermediate_result, initial_point=init_point)
    result = vqe.compute_minimum_eigenvalue(operator=op)

    print(f"VQE on Aer qasm simulator (with noise): {result.eigenvalue.real:.5f}")
    print(f"Delta from reference energy value is {(result.eigenvalue.real - ref_value):.5f}")
    print(f"Delta init from reference energy value is {(values[0] - ref_value):.5f}")
    print(f"Delta from hf is {(result.eigenvalue.real - en_hf):.5f}")

    return result.eigenvalue.real, counts, values, params
