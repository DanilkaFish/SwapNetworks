import unittest

# from BaseTree import BaseTernaryTree, QubitNum, BranchNum, NodeContacts, LostNum
# from TernaryTree import TernaryTree, pauli_weight, prod_pauli_strings
from .tree import tree

# from TreeVizualization import draw
from openfermion.ops import FermionOperator
from openfermion.utils.operator_utils import count_qubits


class Test_tree(unittest.TestCase):
    def setUp(self):
        self.term = FermionOperator(((3, 1), (1, 0)))
        self.n_qubits = count_qubits(self.term)
        # my_term = FermionOperator('3^ 1')

    def test_tree(self):
        print(tree(self.term, self.n_qubits))



if __name__ == "__main__":
    unittest.main()