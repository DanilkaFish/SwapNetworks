
"""The variational quantum eigensolver algorithm."""

from __future__ import annotations

import logging
from time import time
from collections.abc import Callable
from typing import Any

import numpy as np

from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info.operators.base_operator import BaseOperator

from qiskit_algorithms import VQE
# from qiskit import BaseEstimator

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
    
  