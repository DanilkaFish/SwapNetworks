from __future__ import annotations
import logging
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
import sys

from qiskit_nature.second_q.mappers import JordanWignerMapper, BravyiKitaevMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCC
from qiskit_nature.second_q.mappers.fermionic_mapper import FermionicMapper

from qiskit import transpile, QuantumCircuit
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit.circuit import Delay
from qiskit_algorithms.gradients import FiniteDiffEstimatorGradient
from qiskit.circuit.library.standard_gates import (
    CZGate,
    IGate,
    RXGate,
    RZGate,
    RZZGate,
    XGate,
    YGate,
    ZGate,
)
from qiskit.synthesis.evolution import synth_pauli_network_rustiq
from qiskit_aer.noise import (
    NoiseModel,
    QuantumError,
    RelaxationNoisePass,
    depolarizing_error,
    kraus_error,
    thermal_relaxation_error,
)
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import Estimator
from qiskit_algorithms import NumPyMinimumEigensolver

from ..ucc.abstractucc import Molecule
from ..ucc.upgccsd import UpGCCSD, LadExcImpl
from .mapper import MajoranaMapper, MajoranaContainer
from .vqe_qiskit_v2 import VQEV1 as VQE
from .qiskit_circ import QiskitCirc, to_excitaions

logger = logging.getLogger(__name__)
CALIBRATION_DIR = Path(__file__).resolve().parent


NOISE_GATES = {"I": IGate(), "X": XGate(), "Y": YGate(), "Z": ZGate()}
SC_SINGLE_QUBIT_DURATION = 0
SC_TWO_QUBIT_DURATION = 68
ION_SINGLE_QUBIT_DURATION = 63000
ION_TWO_QUBIT_DURATION = 650000
SCHEDULED_NOISE_TYPES = {"sc_scheduled", "ion_scheduled"}
HARDWARE_NOISE_TYPES = {"sc", "ion", *SCHEDULED_NOISE_TYPES}
SCHEDULED_RELAXATION_OP_TYPES = [Delay, RXGate, RZGate, CZGate, RZZGate]
PAULI_MATRICES = {
    "X": np.array([[0, 1], [1, 0]]),
    "Y": np.array([[0, -1j], [1j, 0]]),
    "Z": np.array([[1, 0], [0, -1]]),
    "I": np.array([[1, 0], [0, 1]]),
}

def get_file_name(name, noise, method):
    return name + f"_{noise}" + method +".json"

class SwapCircNames:
    SWAP2XN = ("swap2xn12", LadExcImpl.CNOT12xyz())
    SWAPGENYORDAN = ("swap_gen", LadExcImpl.YORDAN())
    SWAPGENYORDAN2XN = ("swap_gen", LadExcImpl.YORDAN2XN())
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
    def swap_gen_yor_2xn():
        return "swap gen yor 2xn"  
    
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
        return [
            Circuits.jw(),
            Circuits.bk(),
            Circuits.jw_lex(),
            Circuits.bk_lex(),
            Circuits.swap_gen_yor(),
            Circuits.swap_gen_yor_2xn(),
            Circuits.swap_2xn(),
            Circuits.swap_gen_short(),
            Circuits.swap_2xn_yor(),
            Circuits.swap_2xn_short(),
            Circuits.swap_gen_maj_yor(),
            Circuits.swap_gen_maj_short(),
        ]
    

class CircuitProvider:
    def __init__(self,
                 reps=1,
                 molecule: Molecule=None,
                 basis_gates=None
                 ):
        if molecule is None:
            molecule = Molecule()
        self.reps = reps
        self.num_electrons = molecule.num_electrons
        self.basis_gates = ["u", "cx"] if basis_gates is None else basis_gates

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
        logger.info("running rust synthesis")
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
        return circ

    def get_circ_op(self, name):
        if name == Circuits.jw():
            return (self.get_circ_via_mapping(JordanWignerMapper()), jw_ham(self.fermionic_op))
        elif name == Circuits.bk():
            return (self.get_circ_via_mapping(BravyiKitaevMapper()), bk_ham(self.fermionic_op))
        
        elif name == Circuits.swap_gen_yor():
            return self.get_swap_circuit(SwapCircNames.SWAPGENYORDAN)
        elif name == Circuits.swap_gen_yor_2xn():
            return self.get_swap_circuit(SwapCircNames.SWAPGENYORDAN2XN)
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
        raise ValueError(f"Unknown circuit name: {name}")
        
    def get_circ(self, name) -> Tuple[str, QuantumCircuit, SparsePauliOp]:
        circ, op = self.get_circ_op(name)
        return name, circ, op

    def __iter__(self, name_list=None):
        if name_list is None:
            name_list = Circuits.get_circs_names()
        for name in name_list:
            yield name, *self.get_circ_op(name)

class Callback:
    def __init__(self):
        self._energy_array = []

    def __call__(self, step: int, params: np.ndarray, energy: float, metadata: dict):
        sys.stdout.write(f"\rProgress: {step}, Energy: {energy}%")
        sys.stdout.flush()
        self._energy_array.append(energy)

    @property
    def energy_array(self):
        return self._energy_array
    

def hardware_instruction_durations(num_qubits, single_qubit_duration, two_qubit_duration):
    instr_dur = []
    for i in range(num_qubits):
        instr_dur.append(("rx", [i], single_qubit_duration))
        instr_dur.append(("rz", [i], single_qubit_duration))
        for j in range(num_qubits):
            if i != j:
                instr_dur.append(("rzz", [i, j], two_qubit_duration))
                instr_dur.append(("cz", [i, j], two_qubit_duration))
    return instr_dur


class CircSim:
    def __init__(self, circ: QiskitCirc, op: SparsePauliOp, noise_par=0.9999, noise_type="D", init_point=None, s_basis=None, d_basis=None):
        s_basis = ["u3"] if s_basis is None else s_basis
        d_basis = ["cx"] if d_basis is None else d_basis
        source_parameters = list(circ.parameters)
        self.circ = circ
        self.op = op
        self.hf = hf(circ, op)
        self.noise_par = noise_par
        logger.info("starting transpilation")
        if noise_type in HARDWARE_NOISE_TYPES:
            if noise_type in {"ion", "ion_scheduled"}:
                td = ION_TWO_QUBIT_DURATION
                ts = ION_SINGLE_QUBIT_DURATION
                cp = None
            else:
                td = SC_TWO_QUBIT_DURATION
                ts = SC_SINGLE_QUBIT_DURATION
                cp = coupling_map_2xn(self.circ.num_qubits//2)
            instr_dur = hardware_instruction_durations(circ.num_qubits, ts, td)
            hardware_basis_gates = ["cz", "rzz", "rx", "rz"]
            decomposed = self.circ.decompose(reps=3)
            if noise_type == "ion":
                mapped = transpile(decomposed,
                                   basis_gates=hardware_basis_gates,
                                   optimization_level=3)
                serialized = serialize_two_qubit_gates(mapped)
                self.circ = transpile(serialized,
                                      basis_gates=hardware_basis_gates,
                                    optimization_level=1, 
                                    coupling_map=cp,
                                    instruction_durations=instr_dur,
                                    layout_method='trivial',
                                    routing_method='basic',
                                    scheduling_method="asap")
            else:
                mapped = transpile(decomposed,
                                    basis_gates=hardware_basis_gates,
                                    optimization_level=3,
                                    coupling_map=cp,
                                    layout_method='trivial',
                                    routing_method='basic',
                                    )
                self.circ = transpile(mapped,
                                    basis_gates=hardware_basis_gates,
                                    optimization_level=0,
                                    coupling_map=cp,
                                    instruction_durations=instr_dur,
                                    layout_method='trivial',
                                    routing_method='basic',
                                    scheduling_method="asap"
                                    )
        else:
            self.circ = transpile(self.circ.decompose(reps=4), 
                              basis_gates=[*s_basis, *d_basis], 
                              optimization_level=2, 
                                coupling_map=coupling_map_2xn(self.circ.num_qubits//2),
                                layout_method='trivial',
                                routing_method='basic',
                                seed_transpiler=42
                              )
            
        logger.info(f"circ number of operators = {self.circ.count_ops()}")
        if self.circ.layout is not None:
            layout = self.circ.layout.final_index_layout()
            self.op = self.op.apply_layout(layout)
        self.noise_type = noise_type
        self.init_point = self._initial_point(init_point, source_parameters)

    def _initial_point(self, init_point, source_parameters):
        target_parameters = list(self.circ.parameters)
        if init_point is None:
            return [0 for _ in target_parameters]
        if isinstance(init_point, dict):
            init_by_name = {getattr(key, "name", key): value for key, value in init_point.items()}
            return [
                init_point.get(parameter, init_by_name.get(parameter.name, 0))
                for parameter in target_parameters
            ]

        values = list(init_point)
        if len(values) != len(target_parameters):
            raise ValueError(
                f"Initial point has {len(values)} values, but circuit has {len(target_parameters)} parameters"
            )
        if source_parameters == target_parameters:
            return values

        values_by_parameter = dict(zip(source_parameters, values))
        values_by_name = {parameter.name: value for parameter, value in values_by_parameter.items()}
        return [
            values_by_parameter.get(parameter, values_by_name.get(parameter.name, 0))
            for parameter in target_parameters
        ]

    def parameter_map(self, values):
        return dict(zip(self.circ.parameters, values))
    
    def run_qiskit_vqe(self, optimizer, device="CPU", reps=1):
        
        if self.noise_type == "":
            est = Estimator(
                run_options={"seed": 170, "shots": None, },
                approximation=True,
                backend_options={"device": device},
            )
            sim = AerSimulator(method="density_matrix", device=device)
        elif self.noise_type == "sc":
            est = get_noise_estimator_from_csv(self.noise_par, device)
            sim = get_noise_estimator_from_csv(self.noise_par, device, sim=True)
        elif self.noise_type == "sc_scheduled":
            est = get_noise_estimator_from_csv(
                self.noise_par,
                device,
                scheduled=True,
                num_qubits=self.circ.num_qubits,
            )
            sim = get_noise_estimator_from_csv(
                self.noise_par,
                device,
                sim=True,
                scheduled=True,
                num_qubits=self.circ.num_qubits,
            )
        elif self.noise_type == "ion":
            est = get_ion_noise_estimator(self.noise_par, device, self.circ.num_qubits)
            sim = get_ion_noise_estimator(self.noise_par, device, self.circ.num_qubits, sim=True)
        elif self.noise_type == "ion_scheduled":
            est = get_ion_noise_estimator(
                self.noise_par,
                device,
                self.circ.num_qubits,
                scheduled=True,
            )
            sim = get_ion_noise_estimator(
                self.noise_par,
                device,
                self.circ.num_qubits,
                sim=True,
                scheduled=True,
            )
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
        grad = FiniteDiffEstimatorGradient(est, method="forward", epsilon=1e-7)
        vqe = VQE(est, self.circ, optimizer=optimizer, gradient=grad, initial_point=self.init_point, callback=cb)
        result = vqe.compute_minimum_eigenvalue(operator=self.op)
        for i in range(reps-1):
            vqe.initial_point = np.random.rand(len(self.init_point)) - 0.5
            
            _res = vqe.compute_minimum_eigenvalue(operator=self.op)
            if result.eigenvalue.real > _res.eigenvalue.real:
                result = _res
                logger.info(f"{vqe.initial_point=}")
        result_by_name = {parameter.name: value for parameter, value in result.optimal_parameters.items()}
        parameters = []
        for parameter in self.circ.parameters:
            if parameter in result.optimal_parameters:
                parameters.append(result.optimal_parameters[parameter])
            else:
                parameters.append(result_by_name[parameter.name])
        return result.eigenvalue.real, parameters, est, sim, cb
        
    def run_adapt_vqe(self, optimizer, device="CPU", reps=1, is_rust=False, cp: CircuitProvider=None, mapper=None):
        par_used = {par: np.random.random() - 0.5 for par in self.circ.parameters if par not in self.circ.excitation_pos}
        pars_pull = set(self.circ.excitation_pos)
        if self.noise_type == "":
            est = Estimator(
                run_options={"seed": 170, "shots": None, },
                approximation=True,
                backend_options={"device": device},
            )
        elif self.noise_type == "sc":
            est = get_noise_estimator_from_csv(self.noise_par, device)
        elif self.noise_type == "sc_scheduled":
            est = get_noise_estimator_from_csv(
                self.noise_par,
                device,
                scheduled=True,
                num_qubits=self.circ.num_qubits,
            )
        elif self.noise_type == "ion":
            est = get_ion_noise_estimator(self.noise_par, device, self.circ.num_qubits)
        elif self.noise_type == "ion_scheduled":
            est = get_ion_noise_estimator(
                self.noise_par,
                device,
                self.circ.num_qubits,
                scheduled=True,
            )
        else:
            est = get_qiskit_device_noise_estimator(
                                                    noise_op=self.noise_type, 
                                                    prob=self.noise_par,
                                                    device=device
                                                    )
        def get_param(en=0) -> None:
            best_parameter = None
            best_value = None
            for parameter in pars_pull:
                pars = list(par_used.keys())
                pars.append(parameter)
                qc = to_excitaions(self.circ, pars, self.circ.excitation_pos)
                if is_rust:
                    qc = cp.optimize_circ(qc, qubit_mapper=mapper)

                qc = qc.assign_parameters(par_used, inplace=False)
                vqe = VQE(est, qc, optimizer=optimizer, initial_point=[0])
                res = vqe.compute_minimum_eigenvalue(operator=self.op)
                if res.eigenvalue < en:
                    en = res.eigenvalue
                    best_value = list(res.optimal_parameters.values())[0]
                    best_parameter = parameter
            if best_parameter is not None:
                par_used[best_parameter] = best_value
                pars_pull.discard(best_parameter)
        en = 0
        _res = 0
        counter = 0
        while pars_pull:
            if len(par_used) == 0:
                get_param(0)
            qc = to_excitaions(self.circ, par_used, self.circ.excitation_pos)
            init_point = [par_used[par] for par in qc.parameters]
            if is_rust:
                qc = cp.optimize_circ(qc, qubit_mapper=mapper)
            if self.noise_type in {"sc", "sc_scheduled"}:
                qc = transpile_to_sc(qc)
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
            par_used = result.optimal_parameters
            get_param(0)
        return _res.eigenvalue.real, list(_res.optimal_parameters.values())

def coupling_map_2xn(n):
    return None

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

def get_qiskit_device_noise_estimator(noise_op, prob, device, sim=False, s_basic=None, d_basis=None, shots=None) -> Estimator:
    s_basic = ["u"] if s_basic is None else s_basic
    d_basis = ["cx"] if d_basis is None else d_basis
    basis_gates = s_basic + d_basis
    noise_model = NoiseModel(basis_gates=basis_gates)
    if noise_op=="D":
        error1 = depolarizing_error(1 - prob, 1)
        error2 = depolarizing_error(1 - prob, 2)
    else:
        error1 = QuantumError([(NOISE_GATES["I"], prob), (NOISE_GATES[noise_op], 1 - prob)])
        error2 = error1.tensor(error1)
        op = PAULI_MATRICES[noise_op]
        I = PAULI_MATRICES["I"]
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

def mean_gate_error(df, name="CZ error"):
    cz_column = df[name].dropna()
    cz_errors = []
    for row in cz_column:
        pairs = row.split(';')
        for pair in pairs:
            if ':' in pair:
                _, value = pair.split(':')
                cz_errors.append(float(value))
    average_cz_error = sum(cz_errors) / len(cz_errors)
    return average_cz_error

def serialize_two_qubit_gates(circ: QuantumCircuit) -> QuantumCircuit:
    qc = QuantumCircuit(*circ.qregs, *circ.cregs, name=circ.name)
    for inst, qargs, cargs in circ.data:
        if inst.num_qubits == 2 and inst.name != "barrier":
            qc.barrier(*qc.qubits)              
            qc.append(inst, qargs, cargs)
            qc.barrier(*qc.qubits)              
        else:
            qc.append(inst, qargs, cargs)
    return qc


def add_scheduled_relaxation(noise_model, t1_us, t2_us, num_qubits):
    if num_qubits is None:
        raise ValueError("num_qubits is required for scheduled relaxation noise")

    noise_model._custom_noise_passes.append(
        RelaxationNoisePass(
            t1s=[t1_us * 1e-6 for _ in range(num_qubits)],
            t2s=[t2_us * 1e-6 for _ in range(num_qubits)],
            dt=1e-9,
            op_types=SCHEDULED_RELAXATION_OP_TYPES,
        )
    )


def get_noise_estimator_from_csv(mult, device, sim=False, scheduled=False, num_qubits=None):
    file_name = CALIBRATION_DIR / "ibm_kingston_calibrations_2025-07-02T15_30_16Z.csv"
    basis_gates = ["cz", "rzz", "rx", "rz"]
    noise_model = NoiseModel(basis_gates=basis_gates)
    df = pd.read_csv(file_name)
    T1 = df["T1 (us)"].mean() / mult
    T2 = df["T2 (us)"].mean() / mult
    logger.info(f"{T1=}")
    logger.info(f"{T2=}")
    U = df["Pauli-X error"].mean() * mult
    logger.info(f"{U=}")
    CX = mean_gate_error(df) * mult
    logger.info(f"{CX=}")

    error1 = depolarizing_error(4./2 * U, 1)
    error2 = depolarizing_error(4./3 * CX, 2)
    if scheduled:
        add_scheduled_relaxation(noise_model, T1, T2, num_qubits)
    else:
        relax1q = thermal_relaxation_error(T1 * 1000, T2 * 1000, SC_SINGLE_QUBIT_DURATION)
        relax2q = thermal_relaxation_error(T1 * 1000, T2 * 1000, SC_TWO_QUBIT_DURATION)
        relax2q_both = relax2q.expand(relax2q)
        error1 = error1.compose(relax1q)
        error2 = error2.compose(relax2q_both)

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


def get_ion_noise_estimator(mult, device, nq, sim=False, scheduled=False):
    basis_gates = ["cz", "rzz", "rx", "rz"]
    noise_model = NoiseModel(basis_gates=basis_gates)
    T1 = 188/mult*1000000
    T2 = 0.95/mult*1000000
    logger.info(f"{T1=}")
    logger.info(f"{T2=}")
    CX = 0.0062*mult
    U = 0.0002*mult
    logger.info(f"real {CX=}")
    error1 = depolarizing_error(4./2*U, 1)
    error2 = depolarizing_error(4/3*CX, 2)
    if scheduled:
        add_scheduled_relaxation(noise_model, T1, T2, nq)
    else:
        relax1q = thermal_relaxation_error(
            T1 * 1000,
            T2 * 1000,
            ION_SINGLE_QUBIT_DURATION,
        )
        relax2q = thermal_relaxation_error(
            T1 * 1000,
            T2 * 1000,
            ION_TWO_QUBIT_DURATION,
        )
        relax2q_both = relax2q.expand(relax2q)
        error1 = error1.compose(relax1q)
        error2 = error2.compose(relax2q_both)
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
    instr_dur = hardware_instruction_durations(
        circ.num_qubits,
        SC_SINGLE_QUBIT_DURATION,
        SC_TWO_QUBIT_DURATION,
    )
    hardware_basis_gates = ["cz", "rzz", "rx", "rz"]
    mapped = transpile(circ,
                    basis_gates=hardware_basis_gates,
                    optimization_level=3,
                    )
    qc = transpile(mapped,
                    basis_gates=hardware_basis_gates,
                    optimization_level=0,
                    instruction_durations=instr_dur,
                    scheduling_method="asap"
                    )
    return qc
