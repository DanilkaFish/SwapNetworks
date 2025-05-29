from __future__ import annotations

from collections import deque
from collections.abc import Iterator
from typing import Callable, Any, SupportsFloat
import logging
import warnings
from time import time

import numpy as np

from scipy.optimize import OptimizeResult
from qiskit_algorithms.utils import algorithm_globals

from qiskit_algorithms.optimizers.optimizer import Optimizer, OptimizerSupportLevel, OptimizerResult, POINT

# number of function evaluations, parameters, loss, stepsize, accepted
CALLBACK = Callable[[int, np.ndarray, float, SupportsFloat, bool], None]
TERMINATIONCHECKER = Callable[[int, np.ndarray, float, SupportsFloat, bool], bool]

class SOAP(Optimizer):
    def __init__(
        self,
        maxfev: int = 100,
        callback: CALLBACK | None = None,
        termination_checker: TERMINATIONCHECKER | None = None,
    ) -> None:
        super().__init__()

        # general optimizer arguments
        self.maxfev = maxfev
        self.callback = callback
        self.termination_checker = termination_checker

        # runtime arguments
        self._nfev: int | None = None  # the number of function evaluations
        self._smoothed_hessian: np.ndarray | None = None  # smoothed average of the Hessians

    def minimize(
        self,
        fun: Callable[[POINT], float],
        x0: POINT,
        jac: Callable[[POINT], POINT] | None = None,
        bounds: list[tuple[float, float]] | None = None,
    ) -> OptimizerResult:
        nfev = 0
        nit = 0
        def _fun(_x):
            nonlocal nfev
            nfev += 1
            return fun(_x,)

        trajectory = [x0.copy()]
        vec_list = []
        # direction order
        metric = np.abs(x0)
        for i in np.argsort(metric)[::-1]:
            vec = np.zeros_like(x0)
            vec[i] = 1
            vec_list.append(vec)
        vec_list_copy = vec_list.copy()

        e_list = [_fun(trajectory[-1])]
        offset_list = []
        diff_list = []

        scale = 0.05

        while nfev < self.maxfev:
            if len(vec_list) != 0:
                vec = vec_list[0]
                vec_list = vec_list[1:]
            else:
                vec_list = vec_list_copy.copy()
                # continue
                if len(trajectory) < len(vec_list_copy):
                    continue
                p0 = trajectory[-1 - len(vec_list_copy)]
                f0 = e_list[-1 - len(vec_list_copy)]
                pn = trajectory[-1]
                fn = e_list[-1]
                fe = _fun(2 * pn - p0)
                if fe > f0:
                    continue
                average_direction = pn - p0
                if np.allclose(average_direction, 0):
                    continue
                average_direction /= np.linalg.norm(average_direction)
                replace_idx = np.argmax(np.abs(diff_list[-len(vec_list_copy) :]))
                df = np.abs(diff_list[-len(vec_list_copy) :][replace_idx])
                if 2 * (f0 - 2 * fn + fe) * (f0 - fn - df) ** 2 > (f0 - fe) ** 2 * df:
                    continue
                del vec_list[replace_idx]
                vec_list = [average_direction] + vec_list
                vec_list_copy = vec_list.copy()
                continue

            vec_normed = vec / np.linalg.norm(vec)
            x = [-scale, 0, scale]
            es = [None, e_list[-1], None]
            for j in [0, -1]:
                es[j] = _fun(trajectory[-1] + x[j] * vec_normed)
            if np.argmin(es) == 0:
                x = [-scale * 4, -scale, 0, scale]
                es = [None, es[0], es[1], es[2]]
                for j in [0]:
                    es[j] = _fun(trajectory[-1] + x[j] * vec_normed)
            elif np.argmin(es) == 2:
                x = [-scale, 0, scale, scale * 4]
                es = [es[0], es[1], es[2], None]
                for j in [-1]:
                    es[j] = _fun(trajectory[-1] + x[j] * vec_normed)
            else:
                assert np.argmin(es) == 1
            a, b, c = np.polyfit(x, es, 2)
            if np.argmin(es) not in [0, 3]:
                offset = b / 2 / a
                e = c - b**2 / 4 / a
            else:
                # print(a, b)
                offset = -x[np.argmin(es)]
                e = np.argmin(es)
            offset_list.append(offset)
            trajectory.append(trajectory[-1] - offset * vec_normed)
            if len(es) == 3:
                e_list.append(e)
                # nit -= 1
            else:
                e_list.append(_fun(trajectory[-1]))
            diff_list.append(e_list[-1] - e_list[-2])

            if self.callback is not None:
                self.callback(np.copy(x0))

            nit += 1
        print(e_list[-5:])
        print(_fun(trajectory[-1]))
        return OptimizeResult(fun=e_list[-1], x=trajectory[-1], nit=nit, nfev=nfev, success=True)

    def get_support_level(self):
        """Get the support level dictionary."""
        return {
            "gradient": OptimizerSupportLevel.ignored,
            "bounds": OptimizerSupportLevel.ignored,
            "initial_point": OptimizerSupportLevel.required,
        }
 