from __future__ import annotations
from typing import Tuple 
import logging
from numpy.random import shuffle
import numpy as np
import pandas as pd
from qiskit_nature.second_q.mappers import JordanWignerMapper, BravyiKitaevMapper
from qiskit_nature.second_q.operators import FermionicOp 
from qiskit_nature.second_q.circuit.library import HartreeFock, UCC
from qiskit_nature.second_q.mappers.fermionic_mapper import FermionicMapper

from qiskit import transpile, QuantumCircuit
from qiskit.transpiler import CouplingMap
from qiskit.quantum_info import SparsePauliOp, Statevector, Operator
from qiskit.circuit import Parameter, ParameterVector, Delay
# from qiskit.circuit.parametervector import 
from qiskit.circuit.library.standard_gates import IGate, XGate, ZGate, YGate, RZZGate ,CZGate
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.synthesis.evolution import synth_pauli_network_rustiq
from qiskit.transpiler.passes.synthesis.hls_plugins import PauliEvolutionSynthesisRustiq
from qiskit.providers.fake_provider import Fake20QV1, Fake5QV1, GenericBackendV2
from qiskit_aer.noise import (NoiseModel, QuantumError, kraus_error, RelaxationNoisePass,
    pauli_error, depolarizing_error, thermal_relaxation_error)
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import Estimator
from qiskit_aer.primitives.sampler import Sampler
from qiskit.primitives import BackendEstimator
from qiskit.transpiler import generate_preset_pass_manager
from qiskit_algorithms import NumPyMinimumEigensolver

from ..ucc.abstractucc import Molecule
from ..ucc.upgccsd import UpGCCSD, LadExcImpl
# from ..utils import MajoranaContainer, MajoranaMapper
from .mapper import MajoranaMapper, MajoranaContainer
from ..qiskit_interface import VQEV1 as VQE
from .qiskit_circ import QiskitCirc, to_excitaions

# from logger import logger

logger = logging.getLogger(__name__)


noise_dict_qiskit = {"I": IGate(),"X": XGate(), "Y": YGate(), "Z": ZGate()}
tensors_dict = {"X": np.array([[0,1], [1,0]]),"Y":  np.array([[0,-1j], [1j,0]]), "Z":  np.array([[1,0], [0,-1]]),
                "I":  np.array([[1,0], [0,1]])}

trasnpile_backend = Fake20QV1()
trasnpile_backend = Fake5QV1()
pass_manager = None

# pass_manager.pre_init = ffsim.qiskit.PRE_INIT
def get_file_name(name, noise, method):
    return name + f"_{noise}" + method +".json"

class SwapCircNames:
    SWAP2XN = ("swap2xn12", LadExcImpl.CNOT12xyz())
    # SWAP2XN_ALT = ("swap2xn_alt", LadExcImpl.CNOT12zyx())
    SWAPGENYORDAN = ("swap_gen", LadExcImpl.YORDAN())
    SWAPGENSHORT = ("swap_gen", LadExcImpl.SHORT())
    SWAP2XNYORDAN = ("swap2xnferm", LadExcImpl.YORDAN())
    SWAP2XNSHORT = ("swap2xnferm", LadExcImpl.SHORT())
    SWAPGENMAJYORDAN = ("swapgenmaj", LadExcImpl.YORDAN())
    SWAPGENMAJSHORT = ("swapgenmaj", LadExcImpl.SHORT())

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
    def swap_2xn_alt():
        return "swap 2xn alt"

    @staticmethod
    def swap_gen_yor():
        return "swap gen yor"    

    @staticmethod
    def swap_gen_short():
        return "swap gen short"
    
    @staticmethod
    def swap_2xn_short():
        return "swap 2xn short"
    
    @staticmethod
    def swap_2xn_12cnot():
        return "swap 2xn short"
    
    @staticmethod
    def swap_2xn_yor():
        return "swap 2xn yor"

    @staticmethod
    def swap_gen_maj_yor():
        return "swap gen maj yor"    

    @staticmethod
    def swap_gen_maj_short():
        return "swap gen maj short"
    
    @staticmethod
    def get_circs_names():
        circs = []
        circs.append(Circuits.jw())
        circs.append(Circuits.bk())
        circs.append(Circuits.jw_lex())
        circs.append(Circuits.bk_lex())
        circs.append(Circuits.swap_2xn())
        circs.append(Circuits.swap_gen_yor())
        circs.append(Circuits.swap_gen_short())
        circs.append(Circuits.swap_2xn_yor())
        circs.append(Circuits.swap_2xn_short())
        circs.append(Circuits.swap_gen_maj_yor())
        circs.append(Circuits.swap_gen_maj_short())
        return circs


def eq_alpha_beta(qc, reps=1):
    n = qc.num_qubits
    k = n//2
    # params = []
    # for i in range(reps):
    #     theta = ParameterVector("θ" + str(i), qc.num_parameters - k*(k-1)//2)
    #     params = params + [theta[i - k*(k-1)//2] for i in range(k*(k-1), 3*k*(k-1)//2)] +\
    #          [theta[i] for i in range(k*(k-1)//2)] + [theta[i] for i in range(k*(k-1)//2)] 
    # qc.assign_parameters(params, inplace=True)

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

    def update_molecule(self, molecule):
        self.uccgsd = UpGCCSD(molecule=molecule)
        self.fermionic_op = self.uccgsd.mol.hamiltonian.second_q_op()
        
    def to_qiskit_excitations(self, **kwargs):
        ls1 = []
        ls2 = []
        # return "sd"
        for el in self.al:
            ls1.append(((el[1],), (el[0],)))
        for el in self.be:
            ls1.append(((el[1],), (el[0],)))
        for el in self.do:
            if el[0] < el[2]:
                ls2.append(((el[3],el[2]), (el[1],el[0])))
        # ls = [((3,2), ())]
        # return ls2
        # if len(self.do) == 6:
        #     return [((1, 5), (0, 4)), ((1,), (0,)), ((5,), (4,)), 
        #             ((3, 7), (2, 6)), ((3,), (2,)), ((7,), (6,)),
        #             ((3, 7), (0, 4)), ((3,), (0,)), ((7,), (4,)),
        #             ((3, 7), (1, 5)), ((3,), (1,)), ((7,), (5,)),
        #             ((2, 6), (0, 4)), ((2,), (0,)), ((6,), (4,)),
        #             ((2, 6), (1, 5)), ((2,), (1,)), ((6,), (5,))]
        return ls2 + ls1  

    def get_swap_circuit(self, name: Tuple[str,str]) -> Tuple[QuantumCircuit, SparsePauliOp]:
        method, double = name
        cirq, mtoq = getattr(self.uccgsd, method)(self.reps, double)
        return cirq, ucc_ham(self.fermionic_op, mtoq)
    
    def get_ucc(self, qubit_mapper,  init=False):
        num_electrons = (self.uccgsd.n_alpha, self.uccgsd.n_beta)
        initial_state = None if not init else HartreeFock(self.uccgsd.n_spatial, 
                                                        num_electrons,
                                                        qubit_mapper)
        return UCC(
            self.uccgsd.n_spatial, 
            num_electrons,
            excitations=self.to_qiskit_excitations, 
            qubit_mapper=qubit_mapper, 
            initial_state=initial_state,
            reps=self.reps
            ) 
        
    def get_circ_via_mapping(self, qubit_mapper,  init=True):
        qc = self.get_ucc(qubit_mapper, init)
        eq_alpha_beta(qc, reps=self.reps)
        qc = qc.decompose(reps=1)
        pars_pos = {}
        for index, instr in enumerate(qc):
            if len(instr.params) == 1:
                if instr.params[0] in pars_pos:
                    pars_pos[instr.params[0]].append(index)
                else:
                    pars_pos[instr.params[0]] = [index]
        qc.excitation_pos = pars_pos
        return qc

    def optimize_circ(self, circ: QuantumCircuit, qubit_mapper):
        nq = circ.num_qubits
        init = []
        pauli = []
        for instr in circ:
            if len(instr.params) == 0:
                init.append(instr)
            else:
                pauli.append(instr)
        circ = QuantumCircuit.from_instructions(init, qubits=circ.qubits)
        rq = list(range(nq))
        parr = []
        for gate in pauli:
            for pauli, coef in zip(gate.operation.operator.paulis, gate.operation.operator.coeffs):
                parr.append((str(pauli), list(reversed(rq)), coef*gate.params[0]))
        
        circ.compose(synth_pauli_network_rustiq(nq, 
                                                    parr, 
                                                    optimize_count=False, 
                                                    preserve_order=False, 
                                                    upto_phase=True, 
                                                    upto_clifford=False, 
                                                    resynth_clifford_method=1
                                                ),
                                                inplace=True)
        return circ
    

    def get_rust_circ(self, qubit_mapper: FermionicMapper):
        logger.info("runnint rust sythesis")
        circ = HartreeFock(self.uccgsd.n_spatial, (self.uccgsd.n_alpha, self.uccgsd.n_beta), qubit_mapper)
        circ._build()
        
        qc = self.get_ucc(qubit_mapper, init=None)
        nq = qc.num_qubits
        rq = list(range(nq))
        parr = []
        for gate in qc.decompose(reps=1):
            
            for pauli, coef in zip(gate.operation.operator.paulis, gate.operation.operator.coeffs):
                parr.append((str(pauli), list(reversed(rq)), coef*gate.params[0]))
        circ.compose(synth_pauli_network_rustiq(nq, 
                                                    parr, 
                                                    optimize_count=True,
                                                    preserve_order=True, 
                                                    upto_phase=True, 
                                                    upto_clifford=False, 
                                                    resynth_clifford_method=1
                                                ),
                                                inplace=True)

        
        if pass_manager is not None:
            ansatz = pass_manager.run(qc)
        else:
            ansatz = transpile(circ.decompose(reps=3), basis_gates=self.basis_gates, optimization_level=3).decompose(reps=3)
        return circ

    def get_circ_op(self, name):
        if name == Circuits.jw():
            return (self.get_circ_via_mapping(JordanWignerMapper()), jw_ham(self.fermionic_op))
        elif name == Circuits.bk():
            return (self.get_circ_via_mapping(BravyiKitaevMapper()), bk_ham(self.fermionic_op))
        
        elif name == Circuits.swap_gen_yor():
            return self.get_swap_circuit(SwapCircNames.SWAPGENYORDAN)
        elif name == Circuits.swap_gen_short():
            return self.get_swap_circuit(SwapCircNames.SWAPGENSHORT)
        elif name == Circuits.swap_2xn_yor():
            return self.get_swap_circuit(SwapCircNames.SWAP2XNYORDAN)
        elif name == Circuits.swap_2xn_short():
            return self.get_swap_circuit(SwapCircNames.SWAP2XNSHORT)
        elif name == Circuits.jw_lex():
            return (self.get_rust_circ(JordanWignerMapper()), jw_ham(self.fermionic_op))
        elif name == Circuits.bk_lex():
            return (self.get_rust_circ(BravyiKitaevMapper()), bk_ham(self.fermionic_op))
        elif name == Circuits.swap_gen_maj_yor():
            return self.get_swap_circuit(SwapCircNames.SWAPGENMAJYORDAN)
        elif name == Circuits.swap_gen_maj_short():
            return self.get_swap_circuit(SwapCircNames.SWAPGENMAJSHORT)
        elif name == Circuits.swap_2xn():
            return self.get_swap_circuit(SwapCircNames.SWAP2XN)
        # elif name == Circuits.swap_2xn_alt():
        #     return self.get_swap_circuit(SwapCircNames.SWAP2XN_ALT)
        
    def get_circ(self, name) -> Tuple[str, QuantumCircuit, SparsePauliOp]:
        circ, op = self.get_circ_op(name)
        # logger.info(f"\n{circ.decompose()}")

        return name, circ, op

    def __iter__(self, name_list=Circuits.get_circs_names()):
        for name in name_list:
            yield name, *self.get_circ_op(name)

class Callback:
    def __init__(self):
        self._energy_array = []

    def __call__(self, step: int, params: np.ndarray, energy: float, metadata: dict):
        self._energy_array.append(energy)

    @property
    def energy_array(self):
        return self._energy_array
    

class CircSim:
    def __init__(self, circ: QiskitCirc, op: SparsePauliOp, noise_par=0.9999, noise_type="D",  init_point=None, s_basis=["u3"], d_basis=["cx"]):
        self.circ = circ
        self.op = op
        # self.s_basis = []
        self.hf = hf(circ, op)
        if init_point is None:
            self.init_point = [0 for _ in circ.parameters]
        else:
            self.init_point = init_point
        self.noise_par = noise_par
        logger.info("starting transpilation")
        instr_dur = []

        if noise_type in {"sc", "ion"}:
            instr_dur = []
            if noise_type == "ion":
                td = 650000
                ts = 63000
                cp = None
            else:
                td = 68
                ts = 0
                cp = coupling_map_2xn(self.circ.num_qubits//2)
            for i in range(circ.num_qubits):
                instr_dur.append(("rx", [i], ts))
                instr_dur.append(("rz", [i], ts))
                for j in range( circ.num_qubits):
                    if i != j:
                        instr_dur.append(("rzz", [i, j], td))
                        instr_dur.append(("cz", [i, j], td))
            if noise_type == "ion":
                mapped = transpile(self.circ.decompose(reps=3),
                                   basis_gates=["cz", "rzz", "rx", "rz"],
                                   optimization_level=3)
                serialized = serialize_two_qubit_gates(mapped)
                self.circ = transpile(serialized,
                                      basis_gates=["cz", "rzz", "rx", "rz"],
                                    optimization_level=1, 
                                    coupling_map=cp,
                                    instruction_durations=instr_dur,
                                    layout_method='trivial',  # or 'dense' if you prefer
                                    routing_method='basic',
                                    scheduling_method="asap")
                logger.info(f"{self.circ}")
            else:
                self.circ = transpile(self.circ.decompose(reps=3), 
                                    basis_gates=["cz", "rzz", "rx", "rz"],
                                    optimization_level=3, 
                                    coupling_map=cp,
                                    instruction_durations=instr_dur,
                                    layout_method='trivial',  # or 'dense' if you prefer
                                    routing_method='basic',
                                    scheduling_method="asap"
                                    )
        else:
            self.circ = transpile(self.circ.decompose(reps=4), 
                              basis_gates=[*s_basis, *d_basis], 
                              optimization_level=2, 
                                coupling_map=coupling_map_2xn(self.circ.num_qubits//2),
                                # instruction_durations=instr_dur,
                                layout_method='trivial',  # or 'dense' if you prefer
                                routing_method='basic',
                                seed_transpiler=42
                              )
            # self.circ.excitation_pos = ep
            # self.circ.delay(100, unit="dt")
            
        logger.info(f"circ number of operators = {self.circ.count_ops()}")
        if self.circ.layout is not None:
            layout = self.circ.layout.final_index_layout()
            self.op = self.op.apply_layout(layout)
        self.noise_type = noise_type

    
    def run_qiskit_vqe(self, optimizer, device="CPU", reps=1):
        
        if self.noise_type == "":
            est = Estimator(
                run_options={"seed": 170, "shots": None, },
                approximation=True,
                backend_options={"device": device},
            )
        elif self.noise_type == "sc":
            est = get_noise_estiamtor_from_csv(self.noise_par, device, self.circ.num_qubits)
            sim = get_noise_estiamtor_from_csv(self.noise_par, device, self.circ.num_qubits, sim=True)
        elif self.noise_type == "ion":
            est = get_ion_noise_estimator(self.noise_par, device, self.circ.num_qubits)
            sim = get_ion_noise_estimator(self.noise_par, device, self.circ.num_qubits, sim=True)
        else:
            est = get_qiskit_device_noise_estimator(
                                                    noise_op=self.noise_type, 
                                                    prob=self.noise_par,
                                                    device=device
                                                    )
            sim = get_qiskit_device_noise_estimator(
                                                    noise_op=self.noise_type, 
                                                    prob=self.noise_par,
                                                    device=device,
                                                    sim=True
                                                    )
        cb = Callback()
        vqe = VQE(est, self.circ, optimizer=optimizer, initial_point=self.init_point, callback=cb)
        result = vqe.compute_minimum_eigenvalue(operator=self.op)
        for i in range(reps-1):
            vqe.initial_point = np.random.rand(len(self.init_point)) - 0.5
            
            _res = vqe.compute_minimum_eigenvalue(operator=self.op)
            if result.eigenvalue.real > _res.eigenvalue.real:
                result = _res
                logger.info(f"{vqe.initial_point=}")
        # print(f"VQE on Aer qasm simulator (with noise): {result.eigenvalue.real:.5f}")
        return result.eigenvalue.real, list(result.optimal_parameters.values()), est, sim, cb
        
    def run_adapt_vqe(self, optimizer, device="CPU", reps=1, is_rust=False, cp: CircuitProvider=None, mapper=None):
        par_used = {par: np.random.random() - 0.5 for par in self.circ.parameters if par not in self.circ.excitation_pos}
        pars_pull = set(self.circ.excitation_pos)
        # init_point = [0]
        if self.noise_type == "":
            est = Estimator(
                run_options={"seed": 170, "shots": None, },
                approximation=True,
                backend_options={"device": device},
            )
        elif self.noise_type == "sc":
            est = get_noise_estiamtor_from_csv(self.noise_par, device)
        elif self.noise_type == "ion":
            est = get_ion_noise_estimator(self.noise_par, device)
        else:
            est = get_qiskit_device_noise_estimator(
                                                    noise_op=self.noise_type, 
                                                    prob=self.noise_par,
                                                    device=device
                                                    )
        def get_param(en=0) -> None:
            for parameter in pars_pull:
                pars = list(par_used.keys())
                pars.append(parameter)
                qc = to_excitaions(self.circ, pars, self.circ.excitation_pos)
                if is_rust:
                    qc = cp.optimize_circ(qc, qubit_mapper=mapper)

                qc = qc.assign_parameters(par_used, inplace=False)
                # if self.noise_type == "sc":
                    # qc = transpile_to_sc(qc)
                # print(qc.num_parameters)
                vqe = VQE(est, qc, optimizer=optimizer, initial_point=[0])
                res = vqe.compute_minimum_eigenvalue(operator=self.op)
                if res.eigenvalue < en:
                    en = res.eigenvalue
                    value = list(res.optimal_parameters.values())[0]
                    par = parameter
            par_used[par] = value
            pars_pull.discard(par)
        en = 0
        _res = 0
        counter = 0
        while pars_pull:
            if len(par_used) == 0:
                get_param(0)
            qc = to_excitaions(self.circ, par_used, self.circ.excitation_pos)
            # qc = self.circ.to_excitaions(par_used)
            init_point = [par_used[par] for par in qc.parameters]
            # print(init_point)
            # qc = transpile(qc,  basis_gates=["cx", "u"], optimization_level=3).decompose(reps=3)
            if is_rust:
                qc = cp.optimize_circ(qc, qubit_mapper=mapper)
            if self.noise_type == "sc":
                qc = transpile_to_sc(qc)
            # print(qc)
            vqe = VQE(est, qc, optimizer=optimizer, initial_point=init_point)
            result = vqe.compute_minimum_eigenvalue(operator=self.op)
            if result.eigenvalue < en:
                en = result.eigenvalue
                _res = result
                counter = 0
            else:
                if counter == 4:
                    break
                counter += 1
            # for index, par in enumerate(qc.parameters):
            par_used = result.optimal_parameters
            # print(qc.count_ops())
            # print(en)
            get_param(0)
        return _res.eigenvalue.real, list(_res.optimal_parameters.values())

def coupling_map_2xn(n):
    cm = [(i, i + 1) for i in range(n-1)]
    cm = cm + [(i, i + 1) for i in range(n, 2*n-1)]
    cm = cm + [(i, 2*n - 1 - i) for i in range(0, n)]
    cm = cm + [(j,i) for (i,j) in cm]
    cm = CouplingMap(cm)
    return None
    return cm

def jw_ham(fermionic_op):
    mapper = JordanWignerMapper()
    qubit_jw_op = mapper.map(fermionic_op)
    return qubit_jw_op, mapper

def bk_ham(fermionic_op):
    mapper = BravyiKitaevMapper()
    qubit_bk_op = mapper.map(fermionic_op)
    return qubit_bk_op, mapper

def ucc_ham(fermionic_op, mtoq: MajoranaContainer):
    mapper = MajoranaMapper(mtoq)
    return mapper.map(fermionic_op), mapper

def get_qiskit_device_noise_estimator(noise_op,  prob, device, sim=False, s_basic=["u"], d_basis=["cx"], shots=None) ->Estimator:
    basis_gates = s_basic + d_basis
    noise_model = NoiseModel(basis_gates=basis_gates)
    if noise_op=="D":
        error1 = depolarizing_error(1 - prob, 1)
        error2 = depolarizing_error(1 - prob, 2)
    else:
        error1 = QuantumError([(noise_dict_qiskit["I"], prob), (noise_dict_qiskit[noise_op], 1 - prob)])
        error2 = error1.tensor(error1)
        op = tensors_dict[noise_op]
        I = tensors_dict["I"]
        error2 = kraus_error([np.sqrt(prob) * np.kron(I, I), np.sqrt(1-prob) * np.kron(op, op)])
    noise_model.add_all_qubit_quantum_error(error2, d_basis)

    if shots is not None:
        noisy_estimator = Estimator(
                run_options={"seed": 170, "shots": shots, },
                backend_options={
                        "noise_model": noise_model,
                        "basis_gates": basis_gates,
                        "method": "density_matrix",
                        "device": device,
                        },
            )
    else:
        noisy_estimator = Estimator(
                run_options={"seed": 170, "shots": None, },
                approximation=True,
                backend_options={
                        "noise_model": noise_model,
                        "basis_gates": basis_gates,
                        "method": "density_matrix",
                        "device": device,
                        },
            )

    if sim:
        return AerSimulator(
                        noise_model=noise_model,
                        basis_gates=basis_gates,
                        method="density_matrix",
                        device=device,
                        )
    else:
        return noisy_estimator

def numpy_energy(fermionic_op, ucc):
    numpy_solver = NumPyMinimumEigensolver()
    result = numpy_solver.compute_minimum_eigenvalue(operator=bk_ham(fermionic_op)[0])
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
    return energy

def mean(df, name="CZ error"):
    # Извлекаем и разбираем колонку "CZ error"
    cz_column = df[name].dropna()

    # Список для хранения всех значений ошибок CZ
    cz_errors = []

    # Парсинг строки формата "X:Y;Z:W" и извлечение только значений Y, W, ...
    for row in cz_column:
        pairs = row.split(';')
        for pair in pairs:
            if ':' in pair:
                _, value = pair.split(':')
                cz_errors.append(float(value))

    # Рассчитываем среднее значение ошибок CZ
    average_cz_error = sum(cz_errors) / len(cz_errors)

    return average_cz_error

def serialize_two_qubit_gates(circ: QuantumCircuit) -> QuantumCircuit:
    qc = QuantumCircuit(*circ.qregs, *circ.cregs, name=circ.name)
    for inst, qargs, cargs in circ.data:
        if inst.num_qubits == 2 and inst.name != "barrier":
            qc.barrier(*qc.qubits)              # не даём ничему “перепрыгнуть” вперёд
            qc.append(inst, qargs, cargs)       # сам 2-кубитный гейт
            qc.barrier(*qc.qubits)              # и не даём ничему встать рядом
        else:
            qc.append(inst, qargs, cargs)
    return qc




def get_noise_estiamtor_from_csv(mult, device, nq, sim=False):
    file_name = "./Ternary_Tree/qiskit_interface/ibm_kingston_calibrations_2025-07-02T15_30_16Z.csv"
    basis_gates = ["cz", "rzz", "rx", "rz"]
    noise_model = NoiseModel(basis_gates=basis_gates)
    df = pd.read_csv(file_name)
    T1 = df["T1 (us)"]
    T2 = df["T2 (us)"]
    T1 = T1.mean()/mult
    T2 = T2.mean()/mult
    logger.info(f"{T1=}")
    logger.info(f"{T2=}")
    U = df["Pauli-X error"].mean()*mult
    logger.info(f"{U=}")
    # U = 0.0003*mult
    CX = mean(df)*mult
    # U = mean(df)*mult
    logger.info(f"{CX=}")
    error1 = depolarizing_error(4./3*(1 - (1-U)**2), 1)
    error2 = depolarizing_error(4./3*(1 - (1-CX)**2), 2)
    t1s = [T1 for prop in range(nq)]
    t2s = [T2 for prop in range(nq)]
    delay_pass = RelaxationNoisePass(
        t1s=[np.inf if x is None else x*1000 for x in t1s],
        t2s=[np.inf if x is None else x*1000 for x in t2s],
        dt=1,
        op_types=[Delay, CZGate, RZZGate],
    )
    noise_model._custom_noise_passes.append(delay_pass)
    noise_model.add_all_qubit_quantum_error(error2, ["cz", "rzz"])
    noise_model.add_all_qubit_quantum_error(error1, ["rz", "rx"])
    noisy_estimator = Estimator(
                run_options={"seed": 170, "shots": None, },
                approximation=True,
                backend_options={
                        "noise_model": noise_model,
                        "basis_gates": basis_gates,
                        "method": "density_matrix",
                        "device": device,
                        },
            )
    if sim:
        return AerSimulator(
                        noise_model=noise_model,
                        basis_gates=basis_gates,
                        method="density_matrix",
                        device=device,
                        )
    else:
        return noisy_estimator

def get_ion_noise_estimator(mult, device, nq, sim=False):
    basis_gates = ["cz", "rzz", "rx", "rz"]
    noise_model = NoiseModel(basis_gates=basis_gates)
    T1 = 188/mult*1000000
    T2 = 0.95/mult*1000000
    logger.info(f"{T1=}")
    logger.info(f"{T2=}")
    CX = 0.0062*mult
    U = 0.0002*mult
    logger.info(f"real {CX=}")
    error1 = depolarizing_error(4./3*(1 - (1-U)**2), 1)
    error2 = depolarizing_error(16./15*(1 - (1-CX)**2), 2)
    t1s = [T1 for prop in range(nq)]
    t2s = [T2 for prop in range(nq)]
    delay_pass = RelaxationNoisePass(
        t1s=[np.inf if x is None else x*1000 for x in t1s],
        t2s=[np.inf if x is None else x*1000 for x in t2s],
        dt=1,
        op_types=[Delay, CZGate, RZZGate],
    )
    noise_model._custom_noise_passes.append(delay_pass)
    noise_model.add_all_qubit_quantum_error(error2, ["cz", "rzz"])
    noise_model.add_all_qubit_quantum_error(error1, ["rz", "rx"])
    noisy_estimator = Estimator(
                run_options={"seed": 170, "shots": None, },
                approximation=True,
                backend_options={
                        "noise_model": noise_model,
                        "basis_gates": basis_gates,
                        "method": "density_matrix",
                        "device": device,
                        },
            )
    if sim:
        return AerSimulator(
                        noise_model=noise_model,
                        basis_gates=basis_gates,
                        method="density_matrix",
                        device=device,
                        )
    else:
        return noisy_estimator

def transpile_to_sc(circ):
    instr_dur = []
    td = 68
    ts = 0
    for i in range(circ.num_qubits):
        instr_dur.append(("rx", [i], ts))
        instr_dur.append(("rz", [i], ts))
        for j in range( circ.num_qubits):
            if i != j:
                instr_dur.append(("rzz", [i, j], td))
                instr_dur.append(("cz", [i, j], td))
    qc = transpile(circ, 
                    basis_gates=["cz", "rzz", "rx", "rz"], 
                    optimization_level=3, 
                    instruction_durations=instr_dur,
                #   dt=1,
                    scheduling_method="asap"
                    )
    return qc


if __name__ == "__main__":
    nq = 4
    rq = list(range(4))
    parr = []
    theta = ParameterVector("θ" + str(0), 1)
    for key, coef in {"YYIZ": 1,"XXIZ": 1,"YYZI": 1,"XXZI": 1,"IZYY": -1,"IZXX": -1,"ZIYY": -1,"ZIXX": -1}.items():
        parr.append((str(key), rq, theta[0]))
    circ = synth_pauli_network_rustiq(nq, parr, optimize_count=False, preserve_order=False, upto_phase=True, resynth_clifford_method=0)
    # ansatz = transpile(circ.decompose(reps=3), trasnpile_backend, basis_gates=["u","cx"], optimization_level=3).decompose(reps=3)