from .pauli_tables import pauli_table_JW, pauli_table_BK, pauli_table_TT
import copy
import numpy as np
from qiskit.quantum_info.operators import Pauli

map_dict = {"JW": pauli_table_JW, "BK": pauli_table_BK, "TT": pauli_table_TT}


def get_all_branches_from_pauli_tables(nmodes, mapper_name="JW"):
    pauli_table = map_dict[mapper_name](nmodes)
    list_tables = enumerate_pauli_table(pauli_table)
    for i in range(len(list_tables)):
        list_tables[i] = puali_table_to_branches(list_tables[i])
    return list_tables


def enumerate_pauli_table(pauli_table):
    s = []
    k1 = []
    k2 = []
    lim = min(len(pauli_table) // 2, 6)
    pauli_table_half1 = pauli_table[:len(pauli_table) // 2]
    pauli_table_half2 = pauli_table[len(pauli_table) // 2:]

    def numerate(pauli_table1, pauli_table2, k1, k2):
        flag = True
        for i, branch in enumerate(pauli_table1):
            pauli_table10 = copy.copy(pauli_table1)
            pauli_table10.pop(i)
            pauli_table20 = copy.copy(pauli_table2)
            pauli_table20.pop(i)
            numerate(pauli_table10, pauli_table20, k1 + [pauli_table1[i]], k2 + [pauli_table2[i]])
            flag = False

        if flag:
            for i in range(lim, len(pauli_table) // 2):
                k1 += [pauli_table_half1[i]]
                k2 += [pauli_table_half2[i]]
            #             s.append(k)
            s.append(k1 + k2)

    numerate(pauli_table_half1[:lim], pauli_table_half2[:lim], k1, k2)
    return s


def puali_table_to_branches(pauli_table):
    n_qubits = len(pauli_table)
    branches = []
    for a in pauli_table:
        g1 = []
        g2 = []
        for n in range(n_qubits):
            if a[0][n].to_label() != 'I':
                g1.append([n, a[0][n].to_label()])
                g2.append([n, a[1][n].to_label()])
        branches.append(g1)
        branches.append(g2)
    return branches
