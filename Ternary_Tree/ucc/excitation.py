from __future__ import annotations
from typing import List, Tuple, Dict, Generator
from abc import ABC, abstractmethod
import itertools

from numpy.typing import ArrayLike


class Excitation(ABC):
    
    @abstractmethod
    def maj_range(self) -> List[MajExcitation]:
        pass

    def __getitem__(self, i: int):
        return self._op[i]
    
    def __len__(self):
        return len(self.op)
    
    def __eq__(self, exc: MajExcitation):
        return exc.op == self.op
    
    def __repr__(self):
        return str(self.op)
    
    def __hash__(self):
        return self.op.__hash__()
    
class MajExcitation(Excitation):
    @property
    def op(self):
        return self._op
    
    @op.setter
    def op(self, op: Tuple[int,...]):
        self._op = _arranging(op)
        self.sign = _parity(op)
    
    def maj_range(self):
        return [self.op]
    
class LadExcitation(Excitation):
    @property
    def op(self):
        return self._op
    
    @op.setter
    def op(self, op: Tuple[int,...]):
        self._op = op
        self.sign = 1
    
    @abstractmethod
    def maj_range(self):
        pass

class SingleMajExcitation(MajExcitation):
    
    def __init__(self, i: int, j: int, sign: int=1):
        self.op = (i,j)
        self.sign = self.sign * sign
    
    
class DoubleMajExcitation(MajExcitation):
    def __init__(self, i: int, j: int, k: int, l: int, sign: int=1):
        self.op = (i, j, k, l)
        self.sign = self.sign * sign

    
class SingleLadExcitation(LadExcitation):
    
    def __init__(self, i:int, j: int, sign: int=1):
        self.op = (i,j)
        self.sign = self.sign * sign
        
    def maj_range(self) -> Generator[Tuple[int, int], int]:
        """
        a_i^+ a_j - a_j^+ a_i
        """
        yield SingleMajExcitation(2*self.op[0], 2*self.op[1], 1)
        yield SingleMajExcitation(2*self.op[0] + 1, 2*self.op[1] + 1, 1)


class DoubleLadExcitation(LadExcitation):
    def __init__(self, i: int, j: int, k: int, l: int, sign: int=1):
        self.op = (i, j, k, l)
        self.sign = self.sign * sign
    
    def maj_range(self):
        dcoef = [[1, 1j], [1, 1j], [1, -1j], [1, -1j]]
        for elems in itertools.product([0, 1], repeat=4):
            coef = 1j
            for i in range(4):
                coef *= (dcoef[i][elems[i] - 1])
            if coef.imag == 0:
                j = [2 * self.op[i] + elems[i] for i in range(4)]
                if any(self.op.count(x) % 2 for x in set(self.op)):
                    yield DoubleMajExcitation(*j, 1j*coef.real)


def _parity(t: tuple):
    if len(set(t)) == 1:
        return 0
    par = 1
    for index, el1 in enumerate(t):
        for el2 in t[index + 1:]:
            if el1 > el2:
                par = par * (-1)
    return par

def _arranging(t: tuple):
    new_t = tuple(sorted(t))
    if len(new_t) == 4:
        new_t = new_t[:2] if new_t[2] == new_t[3] else new_t
    new_t = new_t[2:] if new_t[1] == new_t[0] else new_t
    return new_t
