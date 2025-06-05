from __future__ import annotations
from functools import lru_cache
from typing import Tuple, List
import copy

from qiskit.quantum_info.operators import SparsePauliOp, Pauli
from qiskit_nature import QiskitNatureError
from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.mappers.fermionic_mapper import FermionicMapper
from qiskit_nature.second_q.operators import SparseLabelOp

import numpy as np

from .pauli import MajoranaContainer

def pauli_table(mc: MajoranaContainer):
    """
    This function is used by MajoranaMapper to obtain pauli_table
    """
    pauli_table: List[Tuple[Pauli, complex]] = []
    n = len(mc) // 2
    for i in range(0, len(mc), 2):
        pauli1 = mc[i]
        pauli2 = mc[i + 1]
        pauli_table.append(((Pauli((pauli1.bsf[n:], pauli1.bsf[:n], pauli1.pow)),
                             Pauli((pauli2.bsf[n:], pauli2.bsf[:n], pauli2.pow))), 
                             (0.5, 0.5j)))
    return pauli_table

class MajoranaMapper(FermionicMapper):  # pylint: disable=missing-class-docstring

    def __init__(self, pauli_container: MajoranaContainer):
        """The Teranry Tree fermion-to-qubit mapping. Standart mapping is alpha_beta_tree"""
        self.pauli_container = pauli_container
        super().__init__()


    def pauli_table(self, nmodes):
        pt = copy.copy(pauli_table(self.pauli_container))
        return pt


    @lru_cache(maxsize=32)
    def sparse_pauli_operators(self, nmodes: int) -> tuple[list[SparsePauliOp], list[SparsePauliOp]]:
        """Generates the cached :class:`.SparsePauliOp` terms.

        This uses :meth:`.QubitMapper.pauli_table` to construct a list of operators used to
        translate the second-quantization symbols into qubit operators.

        Args:
            nmodes: the number of modes for which to generate the operators.

        Returns:
            Two lists stored in a tuple, consisting of the creation and annihilation  operators,
            applied on the individual modes.
        """
        times_creation_op = []
        times_annihilation_op = []

        for paulis, coef in self.pauli_table(nmodes):
            real_part = SparsePauliOp(paulis[0], coeffs=coef[0])
            imag_part = SparsePauliOp(paulis[1], coeffs=coef[1])

            # The creation operator is given by 0.5*(X - 1j*Y)
            creation_op = real_part - imag_part
            times_creation_op.append(creation_op)

            # The annihilation operator is given by 0.5*(X + 1j*Y)
            annihilation_op = real_part + imag_part
            times_annihilation_op.append(annihilation_op)
        return (times_creation_op, times_annihilation_op)

    #     @classmethod
    def mode_based_mapping(cls, second_q_op: SparseLabelOp, nmodes: int) -> SparsePauliOp:
        """Utility method to map a `SparseLabelOp` to a `SparsePauliOp` using a pauli table.

        Args:
            second_q_op: the `SparseLabelOp` to be mapped.
            nmodes: the number of modes for which to generate the operators.

        Returns:
            The `SparsePauliOp` corresponding to the problem-Hamiltonian in the qubit space.

        Raises:
            QiskitNatureError: If number length of pauli table does not match the number
                of operator modes, or if the operator has unexpected label content
        """
        times_creation_op, times_annihilation_op = cls.sparse_pauli_operators(nmodes)

        # make sure ret_op_list is not empty by including a zero op
        ret_op_list = [SparsePauliOp("I" * nmodes, coeffs=[0])]
        if not isinstance(second_q_op, list):
            n_second_q_op = [second_q_op]
        else:
            n_second_q_op = second_q_op
        ans = []
        for ops in n_second_q_op:
            ret_op_list = []
            for terms, coeff in ops.terms():
                # 1. Initialize an operator list with the identity scaled by the `coeff`
                ret_op = SparsePauliOp("I" * nmodes, coeffs=np.array([coeff]))

                # Go through the label and replace the fermion operators by their qubit-equivalent, then
                # save the respective Pauli string in the pauli_str list.
                for term in terms:
                    char = term[0]
                    if char == "":
                        break
                    position = int(term[1])
                    if char == "+":
                        ret_op = ret_op.compose(times_creation_op[position], front=True)
                    elif char == "-":
                        ret_op = ret_op.compose(times_annihilation_op[position], front=True)
                    # catch any disallowed labels
                    else:
                        raise QiskitNatureError(
                            f"FermionicOp label included '{char}'. Allowed characters: I, N, E, +, -"
                        )
                ret_op_list.append(ret_op)
            ans.append(SparsePauliOp.sum(ret_op_list).simplify())
        if not isinstance(second_q_op, list):
            return ans[0]
        else:
            return ans


    def map(self, second_q_op: FermionicOp,
        *,
        register_length: int | None = None,) -> SparsePauliOp:
        if isinstance(second_q_op,list):
            register_length=second_q_op[0].register_length
        else:
            register_length=second_q_op.register_length
        return self.mode_based_mapping(second_q_op, nmodes=register_length)
    