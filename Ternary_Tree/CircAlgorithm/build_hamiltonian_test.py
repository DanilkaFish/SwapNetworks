import unittest
from UpUCCSDG import UpUCCSDG


class TestAbstractAnsatz(unittest.TestCase):

    def setUp(self):
        self.ucc = UpUCCSDG()

    def test_initialization(self):
        self.ucc
        pass


if __name__ == "__main__":
    unittest.main()



