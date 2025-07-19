from __future__ import annotations
from typing import List, Tuple, Dict, Generator
from abc import ABC, abstractmethod
import itertools

from numpy.typing import ArrayLike


class MajExcitation:
    def __init__(self, op: Tuple[int,...], sign: int=1):
        self.op = op
        self.sign = self.sign * sign
        
    def maj_range(self) -> List[MajExcitation]:
        return [self.op]

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
    
    @property
    def op(self):
        return self._op
    
    @op.setter
    def op(self, op: Tuple[int,...]):
        self._op = tuple(sorted(op))
        self.sign = _parity(op)

    
class LadExcitation(ABC):
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
    
    @property
    def op(self):
        return self._op
    
    @op.setter
    def op(self, op: Tuple[int,...]):
        self._op = op
        self.sign = 1

    
class SingleLadExcitation(LadExcitation):
    """
        a_i^+ a_j - a_j^+ a_ = (i, j)
    """
    def __init__(self, init: Tuple[int], finit: Tuple[int], sign: int=1):
        self.op = (*finit, *init)
        self.sign = self.sign * sign
        
    def maj_range(self) -> Generator[Tuple[int, int], int]:
        """
        a_i^+ a_j - a_j^+ a_i
        """
        yield MajExcitation((2*self.op[0], 2*self.op[1]), 1)
        yield MajExcitation((2*self.op[0] + 1, 2*self.op[1] + 1), 1)


class DoubleLadExcitation(LadExcitation):
    """
        a_i^+ a_j^+ a_k a_l - a_i a_j a_k^+ a_l^+ = (i, j, k, l)
    """
    def __init__(self, init: Tuple[int, int], finit: Tuple[int, int], sign: int=1):
        self.op = (*init, *finit )
        self.sign = self.sign * sign
    
    def maj_range(self):
        dcoef = [[1, 1j], [1, 1j], [1, -1j], [1, -1j]]
        for elems in itertools.product([0, 1], repeat=4):
            coef = 1
            for i in range(4):
                coef *= (dcoef[i][elems[i] - 1])
            if coef.real == 0:
                j = [2 * self.op[i] + elems[i] for i in range(4)]
                if any(self.op.count(x) % 2 for x in set(self.op)):
                    yield MajExcitation(j, coef)
                    

def _parity(t: tuple):
    if len(set(t)) == 1:
        return 0
    par = 1
    for index, el1 in enumerate(t):
        for el2 in t[index + 1:]:
            if el1 > el2:
                par = par * (-1)
    return par
