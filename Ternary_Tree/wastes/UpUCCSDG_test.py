import unittest
from unittest.mock import MagicMock
from Ternary_Tree.UCC.UpGCCSD import UpUCCSDG

class TestUpUCCSDG(unittest.TestCase):
    def setUp(self):
        self.upuccsdg = UpUCCSDG(active_orbitals=[0,1,2,3,4,5])
        # self.upuccsdg.n_qubits = 4
        # self.upuccsdg.num_alpha = 2
        # self.upuccsdg.num_beta = 2
        # self.upuccsdg.tt = MagicMock()
        # self.upuccsdg.get_excitations = MagicMock(return_value={})

    def test_get_parametrized_circuit_default(self):
        cirq = self.upuccsdg.get_parametrized_circuit()
        print(cirq)
        
    # def test_get_parametrized_circuit_custom_layers(self):
    #     cirq = self.upuccsdg.get_parametrized_circuit(number_of_layers=2)
        # self.assertIsInstance(cirq, MyCirq)

if __name__ == '__main__':
    unittest.main()
