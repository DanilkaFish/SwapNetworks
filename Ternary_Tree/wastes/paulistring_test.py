import unittest
from Ternary_Tree.UCC.paulistring import PauliStringSign, prod_pauli_strings

class TestPauliStringSign(unittest.TestCase):

    def setUp(self):
        self.pss1 = PauliStringSign({0: 'X', 1: 'Y'})
        self.pss2 = PauliStringSign({0: 'Y', 1: 'X'})
        self.pss3 = PauliStringSign({1: 'Y', 2: 'Y'})

    def test_init(self):
        # Test initialization with no arguments
        pss = PauliStringSign()
        self.assertEqual(pss.ps, {})
        self.assertEqual(pss.sign, 1)

        # Test initialization with ps argument
        pss = PauliStringSign({0: 'X', 1: 'Y'})
        self.assertEqual(pss.ps, {0: 'X', 1: 'Y'})
        self.assertEqual(pss.sign, 1)

        # Test initialization with sign argument
        pss = PauliStringSign(sign=-1)
        self.assertEqual(pss.ps, {})
        self.assertEqual(pss.sign, -1)
        
        # Test initialization with sign argument
        pss = PauliStringSign([(0,"X"), (1,"Y"), (2,"Z")], -1)
        self.assertEqual(pss.ps, {0:"X", 1:"Y", 2:"Z"})
        self.assertEqual(pss.sign, -1)

    def test_setter(self):
        # Test setting ps attribute
        pss = PauliStringSign()
        pss.ps = {0: 'X', 1: 'Y'}
        self.assertEqual(pss.ps, {0: 'X', 1: 'Y'})

        # Test setting sign attribute
        pss = PauliStringSign()
        pss.sign = -1
        self.assertEqual(pss.sign, -1)

    def test_prod_same_qubits(self):
        pss = prod_pauli_strings(self.pss1, self.pss2)
        self.assertEqual(pss.ps, {0: 'Z', 1: 'Z'})
        self.assertEqual(pss.sign, 1)

    def test_prod_diff_qubits(self):
        pss = prod_pauli_strings(self.pss2, self.pss3)
        self.assertEqual(pss.ps, {0: 'Y', 1: 'Z', 2: 'Y'})
        self.assertEqual(pss.sign, 1j)

if __name__ == '__main__':
    unittest.main()