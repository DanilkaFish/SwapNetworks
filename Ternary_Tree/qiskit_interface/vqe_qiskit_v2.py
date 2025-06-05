
"""The variational quantum eigensolver algorithm."""

from __future__ import annotations

import logging
from time import time
from collections.abc import Callable
from collections import defaultdict
from typing import Any
from warnings import warn

import numpy as np


from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info.operators.base_operator import BaseOperator

from qiskit_algorithms import VQE

from qiskit.circuit import QuantumCircuit
from qiskit.primitives import BaseEstimator
from qiskit.quantum_info.operators.base_operator import BaseOperator

# from qiskit_aer.primitives import Estimator, _pauli_expval_with_variance
from qiskit_aer.primitives.estimator import Estimator, _pauli_expval_with_variance, _ExperimentManager

from qiskit_algorithms.gradients import BaseEstimatorGradient
from qiskit.quantum_info import Pauli, PauliList, SparsePauliOp
from qiskit.primitives import BaseEstimatorV1, EstimatorResult
from qiskit.result.models import ExperimentResult

from qiskit_algorithms.exceptions import AlgorithmError
from qiskit_algorithms.list_or_dict import ListOrDict
from qiskit_algorithms.optimizers import Optimizer, Minimizer, OptimizerResult
from qiskit_algorithms.variational_algorithm import VariationalAlgorithm, VariationalResult
from qiskit_algorithms.minimum_eigensolvers import MinimumEigensolver, MinimumEigensolverResult, VQEResult
from qiskit_algorithms.observables_evaluator import estimate_observables
from qiskit_algorithms.utils import validate_initial_point, validate_bounds
from qiskit_algorithms.utils.set_batching import _set_default_batchsize
from qiskit.quantum_info import DensityMatrix
from qiskit import QuantumRegister

from my_utils import Timer

logger = logging.getLogger(__name__)


class VQEV2(VQE):
 
    def _get_evaluate_energy(
        self,
        ansatz: QuantumCircuit,
        operator: BaseOperator,
    ) -> Callable[[np.ndarray], np.ndarray | float]:

        num_parameters = ansatz.num_parameters

        # avoid creating an instance variable to remain stateless regarding results
        eval_count = 0

        @Timer.attach_timer("evaluate_v2_timer")
        def evaluate_energy(parameters: np.ndarray) -> np.ndarray | float:
            nonlocal eval_count

            # handle broadcasting: ensure parameters is of shape [array, array, ...]
            parameters = np.reshape(parameters, (-1, num_parameters)).tolist()
            batch_size = len(parameters)
 
            try:
                job = self.estimator.run(tuple(zip(batch_size * [ansatz], batch_size * [operator], parameters)))
                estimator_result = job.result()
            except Exception as exc:
                raise KeyError("The primitive job to evaluate the gradient failed!") from exc
            values = [item.data.evs for item in estimator_result]
            if self.callback is not None:
                metadata = estimator_result.metadata
                for params, value, meta in zip(parameters, values, metadata):
                    eval_count += 1
                    self.callback(eval_count, params, value, meta)

            energy = values[0] if len(values) == 1 else values

            return energy

        return evaluate_energy

class VQEV1(VQE):
 
    def _get_evaluate_energy(
        self,
        ansatz: QuantumCircuit,
        operator: BaseOperator,
    ) -> Callable[[np.ndarray], np.ndarray | float]:
        """Returns a function handle to evaluate the energy at given parameters for the ansatz.
        This is the objective function to be passed to the optimizer that is used for evaluation.

        Args:
            ansatz: The ansatz preparing the quantum state.
            operator: The operator whose energy to evaluate.

        Returns:
            A callable that computes and returns the energy of the hamiltonian of each parameter.

        Raises:
            AlgorithmError: If the primitive job to evaluate the energy fails.
        """
        num_parameters = ansatz.num_parameters

        # avoid creating an instance variable to remain stateless regarding results
        eval_count = 0

        @Timer.attach_timer("evaluate_v1_timer")
        def evaluate_energy(parameters: np.ndarray) -> np.ndarray | float:
            nonlocal eval_count

            # handle broadcasting: ensure parameters is of shape [array, array, ...]
            parameters = np.reshape(parameters, (-1, num_parameters)).tolist()
            batch_size = len(parameters)

            try:
                job = self.estimator.run(batch_size * [ansatz], batch_size * [operator], parameters)
                estimator_result = job.result()
            except Exception as exc:
                raise KeyError("The primitive job to evaluate the energy failed!") from exc

            values = estimator_result.values

            if self.callback is not None:
                metadata = estimator_result.metadata
                for params, value, meta in zip(parameters, values, metadata):
                    eval_count += 1
                    self.callback(eval_count, params, value, meta)

            energy = values[0] if len(values) == 1 else values

            return energy

        return evaluate_energy
    
class MyEstimator(Estimator):
    def __init__(
        self,
        *,
        backend_options: dict | None = None,
        transpile_options: dict | None = None,
        run_options: dict | None = None,
        approximation: bool = False,
        num_electrons=2,
        num_qubits=4,
        skip_transpilation: bool = False,
        abelian_grouping: bool = True,
    ):
        self.num_electrons = num_electrons
        self.num_qubits = num_qubits
        self.mask = self._gen_mask()
        super().__init__(backend_options=backend_options, 
                         transpile_options=transpile_options,
                         run_options=run_options,
                         approximation=approximation,
                         skip_transpilation=skip_transpilation,
                         abelian_grouping=abelian_grouping)
    def stabilize_via_z_sum(self, dm: DensityMatrix):
        dm._data.flat[self.mask] = 0
        dm._data = dm.data/dm.data.trace()
                         
    def _gen_mask(self):
        size = 1 << self.num_qubits
        ne = self.num_electrons
        indices = np.array([(i,j) for i in range(size) for j in range(size) if i.bit_count() != ne or j.bit_count() != ne])
        flat_indices = indices[:, 0] * size + indices[:, 1]
        return flat_indices

    def _compute_with_approximation(
        self, circuits, observables, parameter_values, run_options, seed
    ):
        # Key for cache
        key = (tuple(circuits), tuple(observables), self._approximation)
        shots = run_options.pop("shots", None)

        # Create expectation value experiments.
        if key in self._cache:  # Use a cache
            parameter_binds = defaultdict(dict)
            for i, j, value in zip(circuits, observables, parameter_values):
                self._validate_parameter_length(value, i)
                for k, v in zip(self._parameters[i], value):
                    if k in parameter_binds[(i, j)]:
                        parameter_binds[(i, j)][k].append(v)
                    else:
                        parameter_binds[(i, j)][k] = [v]
            experiment_manager = self._cache[key]
            experiment_manager.parameter_binds = list(parameter_binds.values())
        else:
            self._transpile_circuits(circuits)
            experiment_manager = _ExperimentManager()
            for i, j, value in zip(circuits, observables, parameter_values):
                self._validate_parameter_length(value, i)
                if (i, j) in experiment_manager.keys:
                    key_index = experiment_manager.keys.index((i, j))
                    circuit = experiment_manager.experiment_circuits[key_index]
                else:
                    circuit = (
                        self._circuits[i].copy()
                        if self._skip_transpilation
                        else self._transpiled_circuits[i].copy()
                    )
                    observable = self._observables[j]
                    if shots is None:
                        circuit.save_expectation_value(observable, self._layouts[i])
                    else:
                        for term_ind, pauli in enumerate(observable.paulis):
                            circuit.save_expectation_value(
                                pauli, self._layouts[i], label=str(term_ind)
                            )
                experiment_manager.append(
                    key=(i, j),
                    parameter_bind=dict(zip(self._parameters[i], value)),
                    experiment_circuit=circuit,
                )

            self._cache[key] = experiment_manager
        result = self._backend.run(
            experiment_manager.experiment_circuits,
            parameter_binds=experiment_manager.parameter_binds,
            **run_options,
        ).result()
        dm: DensityMatrix = result.data()['density_matrix']
        self.stabilize_via_z_sum(dm)
        
        obs: SparsePauliOp = self._observables[0]
        expectation_values = [dm.expectation_value(obs)]
        if shots is not None:
            raise ValueError("only \"shot is None\" is possible")
            expectation_values = [
                result.data(i)["expectation_value"] for i in experiment_manager.experiment_indices
            ]
        metadata = [
            {"simulator_metadata": result.results[i].metadata}
            for i in experiment_manager.experiment_indices
        ]
        

        return EstimatorResult(np.real_if_close(expectation_values), metadata)

class _MyPostProcessing:
    def __init__(
        self,
        result_indices: list[int],
        paulis: list[PauliList],
        coeffs: list[list[float]],
    ):
        self._result_indices = result_indices
        self._paulis = paulis
        self._coeffs = [np.array(c) for c in coeffs]

    def run(self, results: list[ExperimentResult]) -> tuple[float, dict]:
        """Coumpute expectation value.

        Args:
            results: list of ExperimentResult.

        Returns:
            tuple of an expectation value and metadata.
        """
        combined_expval = 0.0
        combined_var = 0.0
        simulator_metadata = []
        for c_i, paulis, coeffs in zip(self._result_indices, self._paulis, self._coeffs):
            if c_i is None:
                # Observable is identity
                expvals, variances = np.array([1], dtype=complex), np.array([0], dtype=complex)
                shots = 0
            else:
                result = results[c_i]
                count = result.data.counts
                shots = sum(count.values())
                # added for compatibility with qiskit 1.4 (metadata as attribute)
                # and qiskit 2.0 (header as dict)
                try:
                    basis = result.header.metadata["basis"]
                except AttributeError:
                    basis = result.header["metadata"]["basis"]
                indices = np.where(basis.z | basis.x)[0]
                measured_paulis = PauliList.from_symplectic(
                    paulis.z[:, indices], paulis.x[:, indices], 0
                )
                expvals, variances = _pauli_expval_with_variance(count, measured_paulis)
                simulator_metadata.append(result.metadata)
            combined_expval += np.dot(expvals, coeffs)
            combined_var += np.dot(variances, coeffs**2)
        metadata = {
            "shots": shots,
            "variance": np.real_if_close(combined_var).item(),
            "simulator_metadata": simulator_metadata,
        }
        return combined_expval, metadata


def _paulis2basis(paulis: PauliList) -> Pauli:
    return Pauli(
        (
            np.logical_or.reduce(paulis.z),  # pylint:disable=no-member
            np.logical_or.reduce(paulis.x),  # pylint:disable=no-member
        )
    )
