from __future__ import annotations

import warnings

from .Nodes import BranchNum, QubitNum, LostNum, NodeContacts, ROOT
from .tree_builder import _build_alphabeta_tree, _build_full_tt, _build_jw_tree, _build_bk_tree, _build_jw_zyx_tree
from .TreeErrors import TreeStructureError

# self.gate_number_to_name = {2: 'X', 0: 'Y', 1: 'Z'}
# self.gate_name_to_number = {'X': 2, 'Y': 0, 'Z': 1}
# self.gate_name_to_number = {'X': 1, 'Y': 2, 'Z': 0}
# self.gate_number_to_name = {1: 'X', 2: 'Y', 0: 'Z'}

dict_prod = {"II": "I", "XX": "I", "YY": "I", "ZZ": "I",
             "XY": 'Z', "YX": 'Z', "XZ": 'Y', "ZX": 'Y', "YZ": 'X', "ZY": 'X',
             "IX": "X", "XI": "X", "YI": "Y", "IY": "Y", "IZ": "Z", "ZI": "Z"}


class BaseTernaryTree:
    """ This is the class for representation ternary tree and its manipulations"""

    def __init__(
            self,
            n_qubits: int = 0,
            nodes: [[int]] = None,
    ):
        """
        Arg:
            n_qubits: number of qubits (number of tree nodes)
            nodes: tree's structure.
            Example:
                [[1], -- first level (root node)
                [2,3,4], -- second level
                [5,_,_,6,_,_,7,8,9]] -- third level
                       1
                     / | \
                    2  3  4
                   /  /  /|\
                  5  6  7 8 9
        Fields:
            self.root: QubitNum -- root of the tree
            self.nodes: {QubitNum : NodeContacts} -- dict of the tree's nodes
        """
        self.gate_name_to_number = {'X': 0, 'Y': 1, 'Z': 2}
        self.gate_number_to_name = {0: 'X', 1: 'Y', 2: 'Z'}
        self.root = None
        self._nodes = {}
        if nodes is not None:
            self.nodes = nodes
        else:
            self.build_jw_tree(n_qubits)

    def build_full_tt(self, n_qubits):
        _build_full_tt(self, n_qubits)

    def build_jw_zyx_tree(self, n_qubits):
        _build_jw_zyx_tree(self, n_qubits)
        
    def build_jw_tree(self, n_qubits):
        _build_jw_tree(self, n_qubits)

    def build_bk_tree(self, n_qubits):
        _build_bk_tree(self, n_qubits)

    def build_alphabeta_tree(self, nmodes):
        _build_alphabeta_tree(self, nmodes)

    @property
    def nodes(self):
        return self._nodes

    @nodes.setter
    def nodes(self, nodes: [[int]]):
        self._nodes = {}
        last_level = [None]
        max_nodes = 1
        for index, level in enumerate(nodes):
            minim = min(len(level), max_nodes)
            number_nodes = 0
            for i in range(minim):
                if isinstance(level[i], int):
                    number_nodes += 1
                    self.add_node(node=level[i],
                                  its_parent=last_level[i // 3],
                                  pos=i % 3)
            last_level = [i for i in level if not (i is None)]
            max_nodes = number_nodes * 3

    @property
    def n_qubits(self) -> int:
        return len(self.nodes)

    def parent(self, node: QubitNum | int) -> QubitNum:
        return self[node].parent

    def is_root(self, node: QubitNum | int) -> bool:
        return (self[node].parent == ROOT)

    def height(self) -> int:
        """
        Return tree's height
        """
        branches = self.branches()
        height = max([len(branch[1]) for branch in branches])
        return height

    def parent_contacts(self, node) -> NodeContacts:
        return self[self[node].parent]

    def parent_to_node_edge(self, node) -> int:
        return self.parent_contacts(node).childs.index(node)

    def __getitem__(self, key: QubitNum | int) -> NodeContacts:
        if QubitNum(key) in self.nodes:
            return self.nodes[QubitNum(key)]
        else:
            raise TreeStructureError("Ref to non existing node")

    def __setitem__(self, key: QubitNum | int, val: NodeContacts):
        self.nodes[QubitNum(key)] = val
        return val

    def __iter__(self):
        return self.nodes.__iter__()

    def add_node(self, node, its_parent, pos=0):
        if self.n_qubits == 0:
            self.root = QubitNum(node)
            self[self.root] = NodeContacts(ROOT)
        else:
            node = QubitNum(node)
            its_parent = QubitNum(its_parent)
            if node not in self.nodes and its_parent in self.nodes:
                if isinstance(self[its_parent][pos], (BranchNum, LostNum)):
                    self[node] = NodeContacts(its_parent)
                    self[its_parent][pos] = node
                else:
                    raise KeyError("Parent Node " + str(its_parent) +
                                   " already has child on " + str(pos))
            else:
                raise KeyError("Self.nodes already has node " + str(node) +
                               " or haven't its_parent " + str(its_parent))

    def delete_node(self, nodes: [Union[int | QubitNum]] | int | QubitNum = None):
        """
        Remove nodes and its childs from parent child
        """
        if isinstance(nodes, (QubitNum, int)):
            self._delete_node(nodes)
        else:
            for node in nodes:
                self._delete_node(node)

    def down_generator(self,
                       node: QubitNum | int = None,
                       do_after: callable = lambda *args, **kwargs: None,
                       **kwargs):
        """
        Generator to pass the tree with the branhces.
        node: Node to start come through
        do_after: action to do after recursive exit
        """
        if node is None:
            node = self.root
        if not node.is_last:
            for index, child in enumerate(self[node]):
                yield index, child, node
                yield from self.down_generator(child,
                                               do_after=do_after,
                                               **kwargs)
        do_after(node, **kwargs)

    def _delete_node(self, node: int | QubitNum = 0):
        """
        Remove node and its childs from nodes
        """
        node = QubitNum(node)

        def do_after(_node, **kwargs):
            if _node in self.nodes:
                self._nodes.pop(_node)

        if node in self.nodes:
            if self[node].parent == ROOT:
                self._nodes = {}
            else:
                index = self.parent_to_node_edge(node)
                self.parent_contacts(node)[index] = LostNum()
                for obj in self.down_generator(node, do_after=do_after):
                    pass

    def update_branchnum(self, enum_list: list = None, renumerate=False):
        if renumerate:
            num = 1
            for parent in self.nodes:
                for i, child in enumerate(self[parent]):
                    if child.is_last:
                        self[parent][i] = BranchNum(num)
                        num += 1
            return True

        if enum_list is None:
            enum_list = range(1, self.n_qubits * 2 + 1)

        pos = 0
        try:
            for i, child, parent in self.down_generator(self.root):
                if isinstance(child, BranchNum):
                    self[parent][i] = BranchNum(enum_list[pos])
                    pos += 1
        except IndexError:
            raise TreeStructureError("Numeration should be from 1 to " + str(len(self.branches())))

        self.check_branch_numeration()

    def check_branch_numeration(self):
        bnum = []
        for node in self:
            for i, child in enumerate(self[node]):
                if isinstance(child, BranchNum):
                    bnum.append(child.num)
        if set(bnum) != set(list(range(1, len(bnum) + 1))):
            raise TreeStructureError("Numeration of the branches should be from 1 to "
                                     + str(len(bnum)) + " without repetitions")

    def _branches(self) -> [{BranchNum: (QubitNum, str)}]:
        branches = {}
        branch = []

        def do_before(parent, child, index, branch, branches):

            if isinstance(child, BranchNum):
                branches[child] = branch + [(parent, self.gate_number_to_name[index])]
            branch.append((parent, self.gate_number_to_name[index]))

        def do_after(node, branch, **kwargs):
            if branch:
                branch.pop(-1)

        for index, child, parent in self.down_generator(self.root,
                                                        do_after=do_after,
                                                        do_before=do_before,
                                                        branch=branch):
            do_before(parent, child, index, branch, branches)
        return branches

    # TODO
    def branches(self) -> [tuple[int, list[tuple[int, str]]]]:
        """
        Convert nodes to list of branches. Used by pauli tables
        """
        branches = []
        _branches = self._branches()
        for branchnode in _branches:
            branches.append([branchnode.num, [(gate[0].num, gate[1]) for gate in _branches[branchnode]]])
        return branches

    def __str__(self):
        def to_str(gate):
            return str(gate[0].num) + gate[1]

        L = self.height()
        s = self._branches()
        branch_nodes = s.keys()
        s = s.values()
        pr = ''
        k = 4
        for l in range(0, L):
            for branch in s:
                if len(branch) > l:
                    pr += to_str(branch[l]) + ' ' * (k - len(branch[l]))
                else:
                    pr += " " * k
            pr = pr + '\n'
        for num in branch_nodes:
            pr += str(num.sign * num.num) + ' ' * (k - len(str(num.sign * num.num)))

        return pr + '\n ---------------------------------------------------------------'
