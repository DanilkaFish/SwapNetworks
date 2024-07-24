from typing import Dict

from .BaseTree import BaseTernaryTree, TreeStructureError
from .Nodes import QubitNum, BranchNum, LostNum, NodeNum

gate_name = {0: 'X', 1: 'Y', 2: 'Z'}
gate_index = {'X': 0, 'Y': 1, 'Z': 2}
dict_prod = {"II": "I", "XX": "I", "YY": "I", "ZZ": "I",
             "XY": 'Z', "YX": 'Z', "XZ": 'Y', "ZX": 'Y', "YZ": 'X', "ZY": 'X',
             "IX": "X", "XI": "X", "YI": "Y", "IY": "Y", "IZ": "Z", "ZI": "Z"}

dict_prod_coef = {"II" : ["I",1] , "XX" :  ["I",1] , "YY" :  ["I",1] ,"ZZ" :  ["I",1] ,
             "XY" :  ["Z", 1j] ,"YX" : ["Z", -1j],"XZ" : ["Y",-1j], "ZX" : ["Y",1j],"YZ" : ["X",1j],"ZY" : ["X",-1j],
            "IX" : ["X",1], "XI" : ["X",1], "YI" : ["Y",1],"IY" : ["Y",1],"IZ" : ["Z",1],"ZI" : ["Z",1]}


def sign_prod(pauli1, pauli2, coef1=1):
    coefpr1 = 1j
    coefpr2 = 1j

    for qubit in pauli1:
        if qubit in pauli2:
            coefpr1 = coefpr1 * (dict_prod_coef[pauli1[qubit] + pauli2[qubit]][1])
            coefpr2 = coefpr2 * (dict_prod_coef[pauli2[qubit] + pauli1[qubit]][1])
#    if abs((coefpr2 - coefpr1).real) < 0.1 and abs((coefpr2 - coefpr1).imag) < 0.1:
    if coefpr1 == coefpr2:
        return coef1
    elif coefpr1.real > 0:
        return -coef1
    else:
        return coef1


class TernaryTree(BaseTernaryTree):

    def tree_up_operator(self, base, index) -> Dict[int, str]:
        base = QubitNum(base)
        maj = {}
        while not self.is_root(base):
            maj[base.num] = gate_name[index]
            child = base
            base = self.parent(base)
            index = self[base].child_index(child)
        maj[base.num] = gate_name[index]
        return maj

    def find_branch_node(self, maj_num: int|NodeNum) -> tuple[QubitNum, int, int]:
        """
        return:
            (parent qubit node, edge index, +-1)
        """
        base = None
        for node in self.nodes:
            for index, child in enumerate(self[node]):
                if child == BranchNum(maj_num):
                    return (node, index, child.sign)
        raise ValueError("Tree doesn't have gamma_" + str(maj_num))

    def get_majorana(
            self,
            maj_num: int
    ) -> tuple[Dict[int, str], int]:
        """
        return: (maj, sign)
        """
        base, index, sign = self.find_branch_node(maj_num)
        maj = self.tree_up_operator(base,index)
        return (maj, sign)

    def branch_transposition(
            self,
            node1: QubitNum | int,
            edge1: int,
            node2: QubitNum | int,
            edge2: int
    ) -> Dict[int, str]:
        """
        Transpose subtrees or branches  with nodes numeration's preserving.
        node1 = QubitNum -- node's number in tree
        edge1 = 0|1|2 -- edge which connects node1 with its parent
        node2 = QubitNum -- node's number in tree
        edge2 = 0|1|2 -- edge which connects node2 with its parent
        Return: pauli operator implementing this transformation.
        """
        def unsigned_prod(gamma1, gamma2):
            gamma = {}
            for node in gamma1:
                if node in gamma2:
                    gamma[node] = dict_prod[gamma1[node]+ gamma2[node]]
                else:
                    gamma[node] = gamma1[node]
            for node in gamma2:
                if node not in gamma1:
                    gamma[node] = gamma2[node]
            return gamma

        if node1 == node2 and edge1 == edge2:
            return ""
        _edge1 = edge1
        _edge2 = edge2
        _node1 = node1
        _node2 = node2

        gamma1 = self.tree_up_operator(node1, edge1)
        gamma2 = self.tree_up_operator(node2, edge2)
        if ((gamma2.get(node1) == gate_name[edge1]) or
            (gamma1.get(node2) == gate_name[edge2])):
            raise KeyError("Impossible transposition")

        gamma = unsigned_prod(gamma1, gamma2)

        for node in self.nodes:
            for child in self[node].childs:
                if isinstance(child,BranchNum):
                    child.sign = sign_prod(gamma, *self.get_majorana(child))

        # Tree transposition
        aux_node = self[_node1][_edge1]
        self[_node1][_edge1] = self[_node2][_edge2]
        if not self[_node2][_edge2].is_last:
            self[self[_node2][_edge2]].parent = _node1
        self[_node2][_edge2] = aux_node
        if not aux_node.is_last:
            self[aux_node].parent = _node2
        return gamma

    def to0vac(self) ->list[Dict[int, str]]:
        """
        Use algorithm described in the graduate work
        return: list of used pauli operators
        """
        def check_xy():
            """
            find unoccupied branch for transposition with Z...Z branch
            """
            for node in self.nodes:
                for index, child in enumerate(self[node]):
                    if not child:
                        return [node, index]
            raise TreeStructureError("Here is now lost branch")

        def descent_before_not_z(_node, _edge):
            _edge = (_edge + 1) % 2
            if not self[_node][_edge].is_last:
                _node = self[_node][_edge]
                _edge = 2
                while not self[_node][_edge].is_last:
                    _node = self[_node][_edge]
            return _node, _edge

        def rise_before_not_z(_node, _edge):
            while _edge == 2:
                # Подъем до не Z edge
                n = _node
                _node = self.parent(_node)
                _edge = self[_node].childs.index(n)
            return _node, _edge

        s = []
        # find Z...Z branch to eliminate it
        child = self[self.root][2]
        while not child.is_last:
            parent = child
            child = self[child][2]
        if child:
            first, index = check_xy()
            s.append(self.branch_transposition(first, index, parent, 2))

        nnodes = self.n_qubits

        for i in range(nnodes):
            branches = self.branches()
            try:
                num = {branches[j][0]: j for j in range(2 * nnodes)}
            except IndexError:
                raise IndexError("You should use 2n branches for numeration, where n = ",
                                 self.n_qubits, "but you use", len(branches))
            except:
                raise Exception("What's the hell??")

            node1, edge1 = branches[num[2 * i + 1]][1][-1]
            node2, edge2 = branches[num[2 * i + 2]][1][-1]
            edge1 = gate_index[edge1]
            edge2 = gate_index[edge2]

            node1, edge1 = rise_before_not_z(node1, edge1)
            node1, edge1 = descent_before_not_z(node1, edge1)

            r = self.branch_transposition(node1, edge1, node2, edge2)
            if len(r) > 0:
                s.append(r)
        return s


def prod_pauli_strings(pauli1: list[tuple[int, str]],
                       pauli2: list[tuple[int, str]]) -> tuple[list[tuple[int, str]], int]:
    if not isinstance(pauli1, dict):
        pauli1 = {gate[0]: gate[1] for gate in pauli1}
    if not isinstance(pauli2, dict):
        pauli2 = {gate[0]: gate[1] for gate in pauli2}
    pauli = {}
    coef = 1
    for qubit in set(list(pauli1.keys()) + list(pauli2.keys())):
        pauli[qubit], sign = dict_prod_coef[pauli1.get(qubit, "I") + pauli2.get(qubit, "I")]
        coef = coef * sign
    return list(pauli.items()), coef


def pauli_weight(tt: TernaryTree, indexes: list[int]) -> float:
    branches = tt.branches()
    d = {}
    for i in indexes:
        d[i] = d.get(i, 0) + 1

    indexes = [key for key in d if d[key] % 2 == 1]
    if len(indexes) > 0:
        branches = [branch[1] for branch in branches if branch[0] in indexes]
        prod = branches[0]
        for pauli in branches[1:]:
            prod, _ = prod_pauli_strings(prod, pauli)
        # print(indexes,":  ", len([p for p in prod if p[1] != "I"]))
        return len([p for p in prod if p[1] != "I"])
    return 0.
