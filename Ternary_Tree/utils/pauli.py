from __future__ import annotations
from typing import List, Optional, Tuple, Union
from copy import copy

from numpy.typing import ArrayLike
import numpy as np


class PauliBin:
    I = (0,0)
    X = (1,0)
    Z = (0,1)
    Y = (1,1)
    stob = {"I": (0,0), "X": (1,0), "Y": (1,1), "Z": (0,1)}
    btos = {(0,0): "I", (1,0): "X", (1,1): "Y", (0,1): "Z"}


class Pauli:
    """
        P = i^{pow - k_1 * k_2} X^{k1} Z^{k2},
        bsf = (k1,k2) 
    """
    def __init__(self, bsf: np.ndarray[bool], pow: int=0):
        self.bsf = bsf
        self.pow = pow

    @classmethod
    def from_str(cls, letters: str, qubits: Optional[List[int]]=None, n_qubits: Optional[int]=None, pow=0):
        if qubits is None:
            qubits = list(range(len(letters)))
        if n_qubits is None:
            n_qubits = len(letters)
        bsf = (letters, qubits, n_qubits)
        return cls(bsf, pow)
    
    @property
    def bsf(self) -> np.ndarray[bool]:
        return self._bsf

    @property
    def bsf_pairs(self):
        n = len(self.bsf) // 2
        return zip(self.bsf[:n], self.bsf[n:])

    @bsf.setter
    def bsf(self, obj: Union[Tuple[str, List[int], Optional[int]], np.ndarray]):
        if isinstance(obj, np.ndarray):
            self._bsf: np.ndarray[bool] = copy(obj)
            return 
        
        letters, qubits, n_qubits = obj
        self._bsf: np.array[bool] = np.zeros(2*n_qubits, dtype=bool)
        for let, num in zip(letters, qubits):
            self._bsf[num], self._bsf[num + n_qubits] = PauliBin.stob[let]

    def __len__(self):
        return len(self.bsf) // 2
    
    def cx(self, ctrl: int, trg: int):
        n = len(self.bsf) // 2
        left = Pauli(np.zeros(2*n, dtype=bool))
        right = Pauli(np.zeros(2*n, dtype=bool))
        if self.bsf[ctrl] == 1:
            left.bsf[trg] = 1
        if self.bsf[trg + n] == 1:
            right.bsf[ctrl + n] = 1
        fin = left * self * right
        self.bsf = fin.bsf
        self.pow = fin.pow

    def __mul__(self, pauli: Pauli):
        new_bsf = self.bsf ^ pauli.bsf
        n = len(self.bsf) // 2
        k1, k2 = self.bsf[:n], self.bsf[n:]
        p1, p2 = pauli.bsf[:n], pauli.bsf[n:]
        pow = (self.pow  + pauli.pow + sum(2*k2*p1 + k1*k2 + p1*p2  - new_bsf[:n] * new_bsf[n:] )) % 4
        return Pauli(new_bsf, pow)
    
    def get_label_qubs(self):
        label, qubits = "", []
        for index, let in enumerate(self._to_str()):
            if let != "I":
                label += let
                qubits.append(index)
        return label, qubits
    
    def get_label_carr(self, qubits):
        label = ""
        pauli = self._to_str()
        for qub in qubits:
            label += pauli[qub]
        return label

    def _to_str(self):
        return ''.join([PauliBin.btos[bin] for bin in self.bsf_pairs])

    def __str__(self):
        return f"e^(i\\pi {self.pow})" + "*" +  self._to_str()
    
    def __repr__(self):
        return f"e^(i\\pi {self.pow})" + "*" +  self._to_str()
        

class PauliContainer:
    def __init__(self, lst: ArrayLike[Pauli]):
        self._paulis: ArrayLike[Pauli] = list(lst)
        self.n = len(self._paulis[0])

    def __len__(self):
        return len(self._paulis)
    
    def __getitem__(self, i: int) -> Pauli:
        return self._paulis[i]
    
    def transform(self, pauli: Pauli):
        for index, pl in enumerate(self._paulis):
            new_pl1 = pl*pauli
            new_pl2 = pauli*pl
            if ((new_pl1.pow - new_pl2.pow) % 4 == 2):
                self._paulis[index] = new_pl1

    def cx_transform(self, ctrl, trg):
        for pauli in self._paulis:
            pauli.cx(ctrl, trg)

    def single_transform(self, qub, name, pow: 0 | 2=0):
        letters = "I"* self.n
        letters[qub] = name
        pauli = Pauli.from_str(letters, pow=pow)
        self.transform(pauli)

    def __str__(self):
        s: str = ""
        for index, pauli in enumerate(self._paulis):
            s += str(pauli) + '\n'
        return s
        

class MajoranaContainer(PauliContainer):
    def __init__(self, pauli_list: ArrayLike[Pauli], qubs):
        self.qubs = qubs
        super().__init__(pauli_list)
        
    @classmethod
    def jw(cls, n_qubits: int, encoding="xyz"):
        qubs = list(range(2*n_qubits))
        encoding = encoding.upper()
        X = Pauli.from_str(encoding[0]).bsf
        Y = Pauli.from_str(encoding[1]).bsf
        Z = Pauli.from_str(encoding[2]).bsf
        bsfx = np.zeros(n_qubits*2, dtype=bool)
        bsfy = np.zeros(n_qubits*2, dtype=bool)
        pauli_list = []
        for i in range(n_qubits):
            bsfx[i], bsfx[i + n_qubits] = X
            bsfy[i], bsfy[i + n_qubits] = Y
            pauli_list.append(Pauli(bsfx))
            pauli_list.append(Pauli(bsfy))
            bsfx[i], bsfx[i + n_qubits] = Z
            bsfy[i],bsfy[i + n_qubits] = Z
        return cls(pauli_list, qubs)
    
    @classmethod
    def jw_alt(cls, n_qubits: int):
        qubs = list(range(2*n_qubits))
        bsfx = np.zeros(n_qubits*2, dtype=bool)
        bsfy = np.zeros(n_qubits*2, dtype=bool)
        pauli_list = []
        for i in range(n_qubits):
            bsfx[i], bsfx[i + n_qubits] = 0, 1
            bsfy[i], bsfy[i + n_qubits] = 1, 1
            pauli_list.append(Pauli(bsfx))
            pauli_list.append(Pauli(bsfy))
            bsfx[i], bsfx[i + n_qubits] = 1, 0
            bsfy[i],bsfy[i + n_qubits] = 1, 0
        return cls(pauli_list, qubs)
    
    
    def transpose(self, qub1: int, i1: int, qub2: int,  i2: int) -> Pauli:
        pauli: Pauli = self[self.qubs[2*qub1 + i1]] * self[self.qubs[2*qub2 + i2]]
        self.transform(pauli)
        self.qubs[2*qub1 + i1], self.qubs[2*qub2 + i2] = self.qubs[2*qub2 + i2], self.qubs[2*qub1 + i1] 
        return pauli
    
    def renumerate(self, list_enum):
        r"""
            pauli[ list_enum[i] ] = pauli[ i ] -- wrong
        """
        _paulis = [None]*len(list_enum)
        _qubs = [None]*len(list_enum)
        for index, num in enumerate(list_enum):
            _paulis[self.qubs[num]] = self._paulis[self.qubs[index]]
            _qubs[index] = self.qubs[num]
        
        self._paulis = _paulis
        self.qubs = _qubs
            
    def get_by_qubs(self, qub: int) -> Tuple[int,int]:
        return (self.qubs[2*qub], self.qubs[2*qub + 1])
    
    def __str__(self):
        s: str = ""
        for index, pauli in enumerate(self._paulis):
            s += f"{self.qubs[index]} : " + f"{index} : " + str(pauli) + '\n'
        return s
    
    def __repr__(self):
        return str(self)
    
if __name__ == "__main__":
    pl1 = Pauli.from_str("XY")
    pl2 = Pauli.from_str("XZ")
    pl1.cx(0,1)
    print(pl1)
    pl2.cx(0,1)
    print(pl2)
    pl2.cx(0,1)
    print(pl2)
    pl1.cx(1,0)
    print(pl1)
    # print(pl1, pl2, pl1*pl2)
    # print(pl1, pl2, pl2*pl1)
    