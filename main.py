import json
import logging
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import Estimator
from qiskit_algorithms.optimizers import L_BFGS_B
from qiskit_nature.second_q.operators import FermionicOp

from Ternary_Tree.qiskit_interface.circuit_provider import (
    CircSim,
    CircuitProvider,
    Circuits,
    get_file_name,
    numpy_energy,
)
from Ternary_Tree.ucc.abstractucc import Molecule
from utils import Timer


logger = logging.getLogger(__name__)
if not logger.handlers:
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(levelname)s %(message)s"))
    logger.addHandler(handler)
logger.setLevel(logging.INFO)


class H2H2(Molecule):
    def __init__(self, distance: float):
        super().__init__(
            geometry=f"H 0 0 0; H 0 0 1.23; H {distance} 0 0; H {distance} 0 1.23",
            num_electrons=(2, 2),
            active_orbitals=[0, 1, 2, 3],
            basis="sto-3g",
        )

    def set_distance(self, distance: float):
        self.geometry = f"H 0 0 0; H 0 0 1.23; H {distance} 0 0; H {distance} 0 1.23"


DH4 = H2H2(1.23)
H2_4 = Molecule(
    geometry="H 0 0 0; H 0 0 0.7349",
    num_electrons=(1, 1),
    active_orbitals=[0, 1],
    basis="sto-3g",
)
H2_8 = Molecule(
    geometry="H 0 0 0; H 0 0 0.7349",
    num_electrons=(1, 1),
    active_orbitals=[0, 1, 2, 3],
    basis="6-31g",
)
LiH_8 = Molecule(
    geometry="H 0 0 0; Li 0 0 1.5459",
    num_electrons=(2, 2),
    active_orbitals=[0, 1, 2, 5],
    basis="sto-3g",
)
LiH_10 = Molecule(
    geometry="H 0 0 0; Li 0 0 1.5459",
    num_electrons=(2, 2),
    active_orbitals=[0, 1, 2, 3, 5],
    basis="sto-3g",
)
LiH_12 = Molecule(
    geometry="H 0 0 0; Li 0 0 1.5459",
    num_electrons=(2, 2),
    active_orbitals=[0, 1, 2, 3, 4, 5],
    basis="sto-3g",
)
LiH_14 = Molecule(
    geometry="H 0 0 0; Li 0 0 1.5459",
    num_electrons=(2, 2),
    active_orbitals=list(range(7)),
    basis="6-31g",
)


OPTIMIZER = (L_BFGS_B(maxiter=150, ftol=0.00000001), "L_BFGS_B")
EXPERIMENT_MOLECULE = LiH_14
AVAILABLE_NOISES = (
    "sc",
    "sc_scheduled",
    "sc_scheduled_1q2q",
    "sc_scheduled_2q",
    "ion",
    "ion_scheduled",
    "ion_scheduled_1q2q",
    "ion_scheduled_2q",
    "D",
    "X",
    "Y",
    "Z",
)
EXPERIMENT_NOISES = ("D", "X", "Y", "Z")
EXPERIMENT_CIRCUITS = (Circuits.jw_lex(), Circuits.bk_lex(), Circuits.swap_2xn(),Circuits.swap_2xn_yor(),Circuits.swap_gen_yor())[2:]
OUTPUT_PREFIX = "data_revised/LIH14"
DEVICE = "GPU"
REPS = 1
DISTANCE = 1.23
NOISY_INITIAL_POINT_MODE = "noiseless"  # "zero", "noiseless", or "previous-noise"
NOISELESS_VQE_RESTARTS = 1
NOISY_VQE_RESTARTS = 1
CACHE_NOISELESS_INITIAL_POINTS = True

SC_ION_MULTIPLIERS = np.array([0.000005, 0.00001, 0.0005, 0.001])
PAULI_PROBS = 1 - np.flip(np.geomspace(0.0000001, 0.00002, 5))
INITIAL_POINT_MODES = {"zero", "noiseless", "previous-noise"}
NOISELESS_INITIAL_POINT_CACHE = {}


def learning_rate(n=100, c=0.5):
    for k in range(1, n + 1):
        yield c / k**0.3


def perturbation(n=100, a=0.01):
    for k in range(1, n + 1):
        yield a / k**0.3


class VQEData:
    def __init__(
        self,
        file_name_to_write: str,
        molecule: Molecule,
        optimizer,
        reps: int = 1,
        noise_type: str = "",
        probs=None,
        device: str = "GPU",
    ):
        self.file_name = file_name_to_write
        self.molecule = molecule
        self.optimizer = optimizer
        self.reps = reps
        self.noise_type = noise_type
        self.probs = probs if probs is not None else 1 - np.flip(np.geomspace(0.00002, 0.001, 15))
        self.circ_prov = CircuitProvider(reps=self.reps, molecule=self.molecule)
        self.ref_value = numpy_energy(self.circ_prov.fermionic_op, self.circ_prov.uccgsd)
        self.device = device

    def update_provider(self):
        self.circ_prov.update_molecule(self.molecule)
        self.ref_value = numpy_energy(self.circ_prov.fermionic_op, self.circ_prov.uccgsd)


def build_probabilities(noise: str):
    if noise.startswith(("sc", "ion")):
        return SC_ION_MULTIPLIERS
    return PAULI_PROBS


def validate_initial_point_mode():
    if NOISY_INITIAL_POINT_MODE not in INITIAL_POINT_MODES:
        raise ValueError(
            f"Unknown NOISY_INITIAL_POINT_MODE={NOISY_INITIAL_POINT_MODE!r}; "
            f"expected one of {sorted(INITIAL_POINT_MODES)}"
        )


def number_observable(num_spin_orbitals: int, mapper):
    ferm_observable = FermionicOp(
        {f"+_{i} -_{i}": 1 for i in range(num_spin_orbitals)},
        num_spin_orbitals=num_spin_orbitals,
    )
    return mapper.map(ferm_observable)


def eval_additional_observable(
    circuit: QuantumCircuit,
    parameters,
    estimator: Estimator,
    simulator: AerSimulator,
    mapper,
):
    circuit = circuit.assign_parameters(parameters)
    circuit.save_density_matrix()

    n_obs = number_observable(circuit.num_qubits, mapper)
    estimator_result = estimator.run([circuit] * 2, observables=[n_obs, n_obs @ n_obs]).result()
    simulator_result = simulator.run(circuit).result()
    rho = simulator_result.data(0)["density_matrix"]
    purity = (rho.data @ rho.data).trace().real

    n_mean = estimator_result.values[0]
    n2_mean = estimator_result.values[1]
    return [n_mean, n2_mean, n2_mean - n_mean**2, purity]


def build_noiseless_initial_point(circuit_name: str, circuit: QuantumCircuit, operator, vqe_data: VQEData):
    cache_key = (
        vqe_data.file_name,
        circuit_name,
        vqe_data.reps,
        vqe_data.optimizer[1],
        vqe_data.device,
    )
    if CACHE_NOISELESS_INITIAL_POINTS and cache_key in NOISELESS_INITIAL_POINT_CACHE:
        return NOISELESS_INITIAL_POINT_CACHE[cache_key]

    logger.info("Finding noiseless initial point for %s", circuit_name)
    simulator = CircSim(circuit, operator, noise_type="")
    energy, parameters, _, _, callback = simulator.run_qiskit_vqe(
        vqe_data.optimizer[0],
        vqe_data.device,
        reps=NOISELESS_VQE_RESTARTS,
    )
    logger.info("%.5f noiseless: name=%s", energy, circuit_name)
    initial_point = simulator.parameter_map(parameters)
    result = {
        "energy": energy,
        "energy_array": callback.energy_array,
        "param": parameters,
    }
    if CACHE_NOISELESS_INITIAL_POINTS:
        NOISELESS_INITIAL_POINT_CACHE[cache_key] = initial_point, result
    return initial_point, result


@Timer.attach_timer("thread_timer")
def evaluate_circuit(circuit_name: str, vqe_data: VQEData, distance: float = 0):
    validate_initial_point_mode()
    data = []
    name, circuit, op_mapper = vqe_data.circ_prov.get_circ(circuit_name)
    operator, mapper = op_mapper
    initial_point = None
    noiseless_result = None

    if NOISY_INITIAL_POINT_MODE in {"noiseless", "previous-noise"}:
        initial_point, noiseless_result = build_noiseless_initial_point(
            name,
            circuit,
            operator,
            vqe_data,
        )

    for index, probability in enumerate(vqe_data.probs):
        simulator = CircSim(
            circuit,
            operator,
            probability,
            vqe_data.noise_type,
            init_point=initial_point,
        )
        energy, parameters, estimator, aer_simulator, callback = simulator.run_qiskit_vqe(
            vqe_data.optimizer[0],
            vqe_data.device,
            reps=NOISY_VQE_RESTARTS,
        )
        additional_result = eval_additional_observable(
            simulator.circ,
            parameters,
            estimator,
            aer_simulator,
            mapper,
        )

        logger.info(
            "%.5f hf: name=%s, noise=%s, index=%s",
            energy,
            name,
            vqe_data.noise_type,
            index,
        )
        data.append(
            {
                "name": name,
                "ref_ener": vqe_data.ref_value,
                "energy": energy,
                "energy_array": callback.energy_array,
                "param": parameters,
                "initial_point_mode": NOISY_INITIAL_POINT_MODE,
                "noiseless_result": noiseless_result,
                "addition_res:": additional_result,
                "optimizer": vqe_data.optimizer[1],
                "gate_count": simulator.circ.count_ops(),
                "dist": distance,
                "prob": probability,
                "noise": vqe_data.noise_type,
            }
        )
        if NOISY_INITIAL_POINT_MODE == "previous-noise":
            initial_point = simulator.parameter_map(parameters)
    return data


def json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def save_results(file_name: str, results):
    path = Path(file_name)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as file:
        json.dump(results, file, indent=4, default=json_default)


def run_vqe(circuit_name: str, vqe_data: VQEData, distance: float = 0):
    result = evaluate_circuit(circuit_name, vqe_data, distance)
    file_name = get_file_name(vqe_data.file_name, vqe_data.noise_type, circuit_name)
    save_results(file_name, result)
    logger.info("thread_timer: %.4f sec", Timer.timers["thread_timer"])
    return result


def build_experiment(noise: str):
    return VQEData(
        OUTPUT_PREFIX,
        EXPERIMENT_MOLECULE,
        OPTIMIZER,
        reps=REPS,
        probs=build_probabilities(noise),
        noise_type=noise,
        device=DEVICE,
    )


def run_configured_experiment():
    for noise in EXPERIMENT_NOISES:
        vqe_data = build_experiment(noise)
        logger.info("Running circuits for noise=%s: %s", noise, ", ".join(EXPERIMENT_CIRCUITS))
        for circuit_name in EXPERIMENT_CIRCUITS:
            run_vqe(circuit_name, vqe_data, DISTANCE)


if __name__ == "__main__":
    run_configured_experiment()
