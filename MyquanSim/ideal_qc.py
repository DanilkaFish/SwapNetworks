import numpy as np
from .instructions import Instruction, CZgate, XgateNoisyC, XgateNoisyD, Xgate, U3, CXgate, U3_without_noise
from itertools import product 
from numpy.typing import ArrayLike 
    

class QuantumCircuit:
    def __init__(self, num_qubits: int):
        self.num_qubits = num_qubits
        self.instructions = []

    def run(self, psi_init, counts: int=3):
        psi_init = psi_init.reshape([2]*self.num_qubits)
        psi_array = np.full((counts, *psi_init.shape), psi_init)

        for instr in self.instructions:
            psi_array = instr.evolve_states(psi_array)
        return psi_array
    
    def cz(self, ctrl_tar: list[int, int]):
        self.add_instr(CZgate(ctrl_tar))
    
    def cx(self, ctrl_tar: list[int,int]):
        self.add_instr(CXgate(ctrl_tar))

    def rx(self, qubits: ArrayLike, theta: float):
        self.add_instr(Xgate(qubits, theta))

    def u3(self, qubits, params):
        self.add_instr(U3(qubits, *params))

    def u3_witout_noise(self, qubits, params):
        self.add_instr(U3_without_noise(qubits, *params))

    def rxd(self, qubits: ArrayLike, theta=0., sigma=0.0001):
        self.add_instr(XgateNoisyD(qubits, theta, sigma))

    def rxc(self, qubits: ArrayLike, theta=0., sigma=0.0001):
        self.add_instr(XgateNoisyC(qubits, theta, sigma))

    def add_instr(self, instr: Instruction):
        self.instructions.append(instr)

    def __str__(self):
        return str(self.instructions)

    def x(self, qubits: ArrayLike):
        self.add_instr(Xgate(qubits, np.pi))
        
def get_probs(psi_array: np.ndarray, ops: np.ndarray) -> tuple[np.ndarray, list[tuple]]:
    """
    Args:
        psi_array: array of evolved states
        ops: array of 2-dim projectors. len(ops) corresponds to number of qubits
    Return:
        probs: array of probabilities with shape=(counts, 2**n_qubits)
        res: array of measurement results
    """
    id = np.eye(2)
    counts = psi_array.shape[0]
    len_psi = np.prod(psi_array.shape[1:])
    probs = np.zeros((len_psi, counts))

    psi_array = psi_array.reshape([counts, len_psi])
    ops = [(op, id - op) for op in ops]
    res = list(product([0, 1], repeat=int(np.log2(len_psi))))

    for index, el in enumerate(product(*ops)):
        meas = np.array(el[0])
        for op in el[1:]:
            meas = np.kron(meas, op)
        probs[index] = np.abs(np.einsum("ki,ij,kj->k", psi_array.conj(), meas, psi_array))
    
    return probs.transpose((1, 0)), res

            
    