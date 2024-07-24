# from .BaseTree import BaseTernaryTree, gate_number_to_name
from .TreeErrors import TreeStructureError

from .Nodes import LostNum, QubitNum, NodeContacts, BranchNum, ROOT
from .bk import bk_majorana_operators
gate_name_to_number = {'X': 0, 'Y': 1, 'Z': 2}
gate_number_to_name = {0: 'X', 1: 'Y', 2: 'Z'}


def _build_alphabeta_tree(
        tt,
        nmodes: int
    ):
    if nmodes > 0:

        t = [(0,)]
        if nmodes % 2 == 0:
            n_qubits = nmodes // 2
            t1 = []
            level = 0
            qubits = list(range(1, n_qubits + 1))
            while 3 ** level + (3 ** level - 1) // 2 < n_qubits:
                t1 = t1 + [tuple(qubits[(3 ** level - 1) // 2:3 ** level + (3 ** level - 1) // 2])]
                level += 1

            t1 = t1 + [tuple(
                qubits[(3 ** level - 1) // 2:n_qubits] + [None] * (3 ** level + (3 ** level - 1) // 2 - n_qubits))]

            t2 = []
            level = 0
            qubits = list(range(n_qubits + 1, 2 * n_qubits + 1))

            while 3 ** level + (3 ** level - 1) // 2 < n_qubits:
                t2 = t2 + [tuple(qubits[(3 ** level - 1) // 2: 3 ** level + (3 ** level - 1) // 2])]
                level += 1

            t2 = t2 + [tuple(qubits[(3 ** level - 1) // 2:n_qubits])]
            t = t + [t1[0] + (None,) + t2[0]]
            for i in range(1, len(t1)):
                t = t + [t1[i] + t2[i]]

            tt.nodes = t
            tt.delete_node(n_qubits)
            tt.update_branchnum()
            node = tt.root
            while not tt[node][2].is_last:
                node = tt[node][2]
            tt[node][2] = LostNum()
            # tt[n_qubits][2] = LostNum()
            # tt[2*n_qubits][2] = LostNum()
            tt.update_branchnum(renumerate=False)
            tt.check_branch_numeration()
        else:
            raise TreeStructureError("There is should be even number of modes")
    else:
        tt._nodes = {}


def _build_full_tt(self, n_qubits):
    if n_qubits > 0:
        self.root = QubitNum(0)
        t = []
        level = 0
        qubits = list(range(0, n_qubits))
        while 3 ** level + (3 ** level - 1) // 2 < n_qubits:
            t = t + [tuple(qubits[(3 ** level - 1) // 2:3 ** level + (3 ** level - 1) // 2])]
            level += 1

        t = t + [tuple(qubits[(3 ** level - 1) // 2:n_qubits])]
        self.nodes = t
        self.update_branchnum(renumerate = True)
        self[QubitNum(n_qubits - 1)][2] = LostNum()
    else:
        self._nodes = {}


def _build_jw_tree(self, n_qubits):
    if n_qubits > 0:
        self.root = QubitNum(0)
        self._nodes[self.root] = NodeContacts(parent=ROOT,
                                              childs=[BranchNum(1), BranchNum(2), 1])
        for i in range(1, n_qubits):
            self.nodes[QubitNum(i)] = NodeContacts(parent=i-1,
                                                   childs=[BranchNum(2*i+1), BranchNum(2*i+2), i+1])
        self._nodes[QubitNum(n_qubits - 1)][2] = LostNum()
        self.update_branchnum()
        self[QubitNum(n_qubits - 1)][2] = LostNum()
        self.check_branch_numeration()
    else:
        self._nodes = {}


def _build_bk_tree(self, n_qubits):
    # TODO
    def down(_branches, node, edge = None):
        nonlocal pos

        if len(_branches) <= 1:
            if len(_branches) == 1:
                self[node][edge] = BranchNum(pos)
                pos += 1
        else:
            qubits_set = {sigma[0] for sigma in _branches[0]}
            for branch in _branches[1:]:
                qubits_set = qubits_set.intersection({sigma[0] for sigma in branch})

            new_node = QubitNum(list(qubits_set)[0])

            self[new_node] = NodeContacts(node)
            if node == ROOT:
                self.root = new_node
            else:
                self[node][edge] = new_node
            for i in range(3):
                new_branches = []
                for branch in _branches:
                    _branch = []
                    flag = False
                    for sigma in branch:
                        if sigma == (new_node.num, gate_number_to_name[i]):
                            flag = True
                        else:
                            _branch = _branch + [sigma]
                    if flag:
                        new_branches = new_branches + [_branch]
                down(new_branches, new_node, i)

    self._nodes = {}

    if n_qubits > 0:
        branches = bk_majorana_operators(n_qubits)
        pos = 1
        self._nodes = {}
        down(branches, ROOT)
    self.check_branch_numeration()
