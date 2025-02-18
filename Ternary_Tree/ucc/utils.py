from __future__ import annotations
from typing import Dict
from copy import deepcopy

from qiskit.circuit import Parameter
from numpy.typing import ArrayLike

from .excitation import LadExcitation, MajExcitation

def lad2maj(ladder_exciations: ArrayLike[LadExcitation], 
                name: str="t_"
                ) -> Dict[MajExcitation, Parameter]:
    ladder_exc_par: Dict[LadExcitation, Parameter] = {}
    for op in ladder_exciations:
        ladder_exc_par[op] = Parameter(name + ','.join([str(i) for i in op]))
        
    maj_exc_par: Dict[MajExcitation, Parameter] = {}
    for op in ladder_exc_par:
        for majop in op.maj_range():
            maj_exc_par[majop] = maj_exc_par.get(majop, 0)  + majop.sign.real * ladder_exc_par[op] 
    return maj_exc_par

def alpha2beta(maj_alpha_par_exc: Dict[MajExcitation, Parameter], 
                    n_qubits: int
                    ) -> Dict[MajExcitation, Parameter]:
    toto = deepcopy(maj_alpha_par_exc)
    for key in maj_alpha_par_exc:
        toto[(key[0] + n_qubits, key[1] + n_qubits)] = toto[key]
    return toto      

def lad2lad(ladder_excitations: ArrayLike[LadExcitation],
            name: str="t_"
            ) -> Dict[LadExcitation, Parameter]:
    ladder_exc_par = {}
    for op in ladder_excitations:
            ladder_exc_par[op] = Parameter(name + ','.join([str(i) for i in op]))
    return ladder_exc_par