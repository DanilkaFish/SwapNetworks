import numpy as np
from qiskit.quantum_info.operators import Pauli
import sys


from ..Tree.TernaryTree import TernaryTree


def _pauli_table_TT(tt: TernaryTree):
    """
    This function is used by InitstateTTInfo for vacuum preparation.
    """


    n_qubits = tt.n_qubits
    branches = tt.branches()

    def majorana_to_str(num):
        pauli_list = ["I"] * n_qubits
        maj, sign = tt.get_majorana(num + 1)
        for el in maj:
            pauli_list[el] = maj[el]
        return ''.join(pauli_list), sign

    pauli_table_str = []
    coef = []
    for i in range(n_qubits):
        str1, sign1 = majorana_to_str(2 * i)
        str2, sign2 = majorana_to_str(2 * i + 1)
        pauli_table_str.append(
            (str1, str2))
        coef.append([0.5 * sign1, 0.5j * sign2])
    return pauli_table_str, coef


def pauli_table_TT(**kwargs):
    """
    This function is used by TernaryTreeMapper to obtain pauli_table
    """
    pauli_table = []
    pauli_table_str, coefs = _pauli_table_TT(**kwargs)
    for i in range(len(pauli_table_str)):
        pauli_table.append(((Pauli(pauli_table_str[i][0][::-1]), Pauli(pauli_table_str[i][1][::-1])), coefs[i]))
    
    return pauli_table


def _pauli_table_JW(nmodes, num_list=None):
    pt = []
    for i in range(nmodes):
        a_z = np.asarray([1] * i + [0] + [0] * (nmodes - i - 1), dtype=bool)
        a_x = np.asarray([0] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
        b_z = np.asarray([1] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
        b_x = np.asarray([0] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
        pt.append(["Z" * i + "X" + "I" * (nmodes - i - 1), "Z" * i + "Y" + "I" * (nmodes - i - 1)])
    return [pt, [1, -1j] * len(pt)]


def pauli_table_JW(nmodes, **kwargs):
    pauli_table = []
    for i in range(nmodes):
        a_z = np.asarray([1] * i + [0] + [0] * (nmodes - i - 1), dtype=bool)
        a_x = np.asarray([0] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
        b_z = np.asarray([1] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
        b_x = np.asarray([0] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
        # c_z = np.asarray([0] * i + [1] + [0] * (nmodes - i - 1), dtype=bool)
        # c_x = np.asarray([0] * nmodes, dtype=bool)
        pauli_table.append((Pauli((a_z, a_x)), Pauli((b_z, b_x))))
    return pauli_table


def pauli_table_BK(nmodes, **kwargs):
    pauli_table = []

    def parity_set(j, n):
        """
        Computes the parity set of the j-th orbital in n modes.

        Args:
            j (int) : the orbital index
            n (int) : the total number of modes

        Returns:
            numpy.ndarray: Array of mode indices
        """
        indices = np.array([])
        if n % 2 != 0:
            return indices

        if j < n / 2:
            indices = np.append(indices, parity_set(j, n / 2))
        else:
            indices = np.append(
                indices, np.append(parity_set(j - n / 2, n / 2) + n / 2, n / 2 - 1)
            )
        return indices

    def update_set(j, n):
        """
        Computes the update set of the j-th orbital in n modes.

        Args:
            j (int) : the orbital index
            n (int) : the total number of modes

        Returns:
            numpy.ndarray: Array of mode indices
        """
        indices = np.array([])
        if n % 2 != 0:
            return indices
        if j < n / 2:
            indices = np.append(indices, np.append(n - 1, update_set(j, n / 2)))
        else:
            indices = np.append(indices, update_set(j - n / 2, n / 2) + n / 2)
        return indices

    def flip_set(j, n):
        """
        Computes the flip set of the j-th orbital in n modes.

        Args:
            j (int) : the orbital index
            n (int) : the total number of modes

        Returns:
            numpy.ndarray: Array of mode indices
        """
        indices = np.array([])
        if n % 2 != 0:
            return indices
        if j < n / 2:
            indices = np.append(indices, flip_set(j, n / 2))
        elif n / 2 <= j < n - 1:
            indices = np.append(indices, flip_set(j - n / 2, n / 2) + n / 2)
        else:
            indices = np.append(
                np.append(indices, flip_set(j - n / 2, n / 2) + n / 2), n / 2 - 1
            )
        return indices

    pauli_table = []
    # FIND BINARY SUPERSET SIZE
    bin_sup = 1
    while nmodes > np.power(2, bin_sup):
        bin_sup += 1
    # DEFINE INDEX SETS FOR EVERY FERMIONIC MODE
    update_sets = []
    update_pauli = []

    parity_sets = []
    parity_pauli = []

    flip_sets = []

    remainder_sets = []
    remainder_pauli = []
    for j in range(nmodes):

        update_sets.append(update_set(j, np.power(2, bin_sup)))
        update_sets[j] = update_sets[j][update_sets[j] < nmodes]

        parity_sets.append(parity_set(j, np.power(2, bin_sup)))
        parity_sets[j] = parity_sets[j][parity_sets[j] < nmodes]

        flip_sets.append(flip_set(j, np.power(2, bin_sup)))
        flip_sets[j] = flip_sets[j][flip_sets[j] < nmodes]

        remainder_sets.append(np.setdiff1d(parity_sets[j], flip_sets[j]))

        update_pauli.append(Pauli((np.zeros(nmodes, dtype=bool), np.zeros(nmodes, dtype=bool))))
        parity_pauli.append(Pauli((np.zeros(nmodes, dtype=bool), np.zeros(nmodes, dtype=bool))))
        remainder_pauli.append(
            Pauli((np.zeros(nmodes, dtype=bool), np.zeros(nmodes, dtype=bool)))
        )
        for k in range(nmodes):
            if np.in1d(k, update_sets[j]):
                update_pauli[j].x[k] = True
            if np.in1d(k, parity_sets[j]):
                parity_pauli[j].z[k] = True
            if np.in1d(k, remainder_sets[j]):
                remainder_pauli[j].z[k] = True

        x_j = Pauli((np.zeros(nmodes, dtype=bool), np.zeros(nmodes, dtype=bool)))
        x_j.x[j] = True
        y_j = Pauli((np.zeros(nmodes, dtype=bool), np.zeros(nmodes, dtype=bool)))
        y_j.z[j] = True
        y_j.x[j] = True
        pauli_table.append(
            (
                parity_pauli[j] & x_j & update_pauli[j],
                remainder_pauli[j] & y_j & update_pauli[j],
            )
        )

    # PauliList has the phase information unlike deprecated PauliTable.
    # Here, phase is unnecessary, so the following removes phase.
    for pauli1, pauli2 in pauli_table:
        pauli1.phase = 0
        pauli2.phase = 0
    return pauli_table
