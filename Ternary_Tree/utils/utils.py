from __future__ import annotations
from typing import Dict
from copy import deepcopy
from numpy.typing import ArrayLike

from .excitation import LadExcitation, MajExcitation, SingleLadExcitation
from time import time


class Parameter:
    """
        var = coef * name
    """
    def __init__(self, name: str='t', coef=1):
        self.name = name
        self.coef = coef

    def __repr__(self):
        return self.name
        
    def __str__(self):
        return self.name
        
    def __mul__(self, val):
        return Parameter(self.name, self.coef*val)
        
    def __rmul__(self, val):
        return Parameter(self.name, self.coef*val)
        
def static_vars(**kwargs):
    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func
    return decorate

def lad2maj(ladder_excitations: ArrayLike[LadExcitation], 
                name: str="t_"
                ) -> Dict[MajExcitation, Parameter]:
    ladder_exc_par: Dict[LadExcitation, Parameter] = {}
        
    for op in ladder_excitations:
        ladder_exc_par[op] = Parameter(name + ','.join([str(i) for i in op]))
        
    maj_exc_par: Dict[MajExcitation, Parameter] = {}
    for op in ladder_exc_par:
        for majop in op.maj_range():
            if majop in maj_exc_par:
                raise KeyError("twice simplification")
            maj_exc_par[majop] = ladder_exc_par[op] 
            maj_exc_par[majop].coef = majop.sign.real
    return maj_exc_par

def alpha2beta(maj_alpha_par_exc: Dict[MajExcitation, Parameter], 
                    n_qubits: int
                    ) -> Dict[MajExcitation, Parameter]:
    toto = deepcopy(maj_alpha_par_exc)
    for key in maj_alpha_par_exc:
        if isinstance(key, MajExcitation):
            toto[MajExcitation((key[0] + n_qubits, key[1] + n_qubits), sign=key.sign)] = toto[key]
        else:
            toto[SingleLadExcitation((key[1] + n_qubits//2,),(key[0] + n_qubits//2,), sign=key.sign)] = toto[key]
    return toto      

def lad2lad(ladder_excitations: ArrayLike[LadExcitation],
            name: str="t_"
            ) -> Dict[LadExcitation, Parameter]:
    ladder_exc_par = {}
    for op in ladder_excitations:
        ladder_exc_par[op] = Parameter(name + ','.join([str(i) for i in op]))
    return ladder_exc_par