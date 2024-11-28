from __future__ import annotations

from qiskit_algorithms.optimizers.optimizer import Optimizer, OptimizerSupportLevel, OptimizerResult, POINT
import numpy as np

from scipy.optimize import OptimizeResult


class SOAP(Optimizer):
    def __init__(
        self,
        maxiter: int = 100,
    ) -> None:
        r"""
        Args:
        """
        super().__init__()
        self.maxiter = maxiter

    def minimize(self, fun, x0, args=(), callback=None, **kwargs):
        """
        Scipy Optimizer interface for sequantial optimization with
        approximate parabola (SOAP)

        Parameters
        ----------
        fun : callable ``f(x, *args)``
            Function to be optimized.
        x0 : ndarray, shape (n,)
            Initial guess. Array of real elements of size (n,),
            where 'n' is the number of independent variables.
        args : tuple, optional
            Extra arguments passed to the objective function.
        maxfev : int
            Maximum number of function evaluations to perform.
            Default: 2000.
        callback : callable, optional
            Called after each iteration.

        Returns
        -------
        res : OptimizeResult
            The optimization result represented as a SciPy ``OptimizeResult`` object.
            Important attributes are: ``x`` the solution array.
        """

        nfev = 0
        nit = 0

        def _fun(_x):
            nonlocal nfev
            nfev += 1
            return fun(_x, *args)

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

        scale = 0.5

        while nfev < self.maxiter:
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
            # if len(es) == 3:
            #     # print(e)
            #     e_list.append(e)
            # else:
            #     print("here")
                # e_list.append(_fun(trajectory[-1]))
            e_list.append(_fun(trajectory[-1]))
            diff_list.append(e_list[-1] - e_list[-2])

            if callback is not None:
                callback(np.copy(trajectory[-1]))

            nit += 1

        result = OptimizerResult()
        result.x = trajectory[-1]
        result.fun = e_list[-1]
        result.fun = _fun(trajectory[-1])
        result.nfev = nfev
        result.nit = nit
        return result
        # return OptimizeResult(fun=e_list[-1], x=trajectory[-1], nit=nit, nfev=nfev, success=True)

    def get_support_level(self):
        """Get the support level dictionary."""
        return {
            "gradient": OptimizerSupportLevel.ignored,
            "bounds": OptimizerSupportLevel.ignored,
            "initial_point": OptimizerSupportLevel.required,
        }
    # def minimize(
    #     self,
    #     fun: Callable[[POINT], float],
    #     x0: POINT,
    #     jac: Callable[[POINT], POINT] | None = None,
    #     bounds: list[tuple[float, float]] | None = None,
    # ) -> OptimizerResult:
    #     # ensure learning rate and perturbation are correctly set: either none or both
    #     # this happens only here because for the calibration the loss function is required
    #     x0 = np.asarray(x0)
    #     if self.learning_rate is None and self.perturbation is None:
    #         get_eta, get_eps = self.calibrate(fun, x0, max_evals_grouped=self._max_evals_grouped)
    #     else:
    #         get_eta, get_eps = _validate_pert_and_learningrate(
    #             self.perturbation, self.learning_rate
    #         )
    #     eta, eps = get_eta(), get_eps()

    #     lse_solver = self.lse_solver
    #     if self.lse_solver is None:
    #         lse_solver = np.linalg.solve

    #     # prepare some initials
    #     x = np.asarray(x0)
    #     if self.initial_hessian is None:
    #         self._smoothed_hessian = np.identity(x.size)
    #     else:
    #         self._smoothed_hessian = self.initial_hessian

    #     self._nfev = 0

    #     # if blocking is enabled we need to keep track of the function values
    #     if self.blocking:
    #         fx = fun(x)  # pylint: disable=invalid-name

    #         self._nfev += 1
    #         if self.allowed_increase is None:
    #             self.allowed_increase = 2 * self.estimate_stddev(
    #                 fun, x, max_evals_grouped=self._max_evals_grouped
    #             )

    #     logger.info("SPSA: Starting optimization.")
    #     start = time()

    #     # keep track of the last few steps to return their average
    #     last_steps = deque([x])

    #     # use a local variable and while loop to keep track of the number of iterations
    #     # if the termination checker terminates early
    #     k = 0
    #     while k < self.maxiter:
    #         k += 1
    #         iteration_start = time()
    #         # compute update
    #         fx_estimate, update = self._compute_update(fun, x, k, next(eps), lse_solver)

    #         # trust region
    #         if self.trust_region:
    #             norm = np.linalg.norm(update)
    #             if norm > 1:  # stop from dividing by 0
    #                 update = update / norm

    #         # compute next parameter value
    #         update = update * next(eta)
    #         x_next = x - update
    #         fx_next = None

    #         # blocking
    #         if self.blocking:
    #             self._nfev += 1
    #             fx_next = fun(x_next)

    #             if fx + self.allowed_increase <= fx_next:  # accept only if loss improved
    #                 if self.callback is not None:
    #                     self.callback(
    #                         self._nfev,  # number of function evals
    #                         x_next,  # next parameters
    #                         fx_next,  # loss at next parameters
    #                         np.linalg.norm(update),  # size of the update step
    #                         False,
    #                     )  # not accepted

    #                 logger.info(
    #                     "Iteration %s/%s rejected in %s.",
    #                     k,
    #                     self.maxiter + 1,
    #                     time() - iteration_start,
    #                 )
    #                 continue
    #             fx = fx_next  # pylint: disable=invalid-name

    #         logger.info(
    #             "Iteration %s/%s done in %s.", k, self.maxiter + 1, time() - iteration_start
    #         )

    #         if self.callback is not None:
    #             # if we didn't evaluate the function yet, do it now
    #             if not self.blocking:
    #                 self._nfev += 1
    #                 fx_next = fun(x_next)

    #             self.callback(
    #                 self._nfev,  # number of function evals
    #                 x_next,  # next parameters
    #                 fx_next,  # loss at next parameters
    #                 np.linalg.norm(update),  # size of the update step
    #                 True,
    #             )  # accepted

    #         # update parameters
    #         x = x_next

    #         # update the list of the last ``last_avg`` parameters
    #         if self.last_avg > 1:
    #             last_steps.append(x_next)
    #             if len(last_steps) > self.last_avg:
    #                 last_steps.popleft()

    #         if self.termination_checker is not None:
    #             fx_check = fx_estimate if fx_next is None else fx_next
    #             if self.termination_checker(
    #                 self._nfev, x_next, fx_check, np.linalg.norm(update), True
    #             ):
    #                 logger.info("terminated optimization at {k}/{self.maxiter} iterations")
    #                 break

    #     logger.info("SPSA: Finished in %s", time() - start)

    #     if self.last_avg > 1:
    #         x = np.mean(np.asarray(last_steps), axis=0)

    #     result = OptimizerResult()
    #     result.x = x
    #     result.fun = fun(x)
    #     result.nfev = self._nfev
    #     result.nit = k

    #     return result