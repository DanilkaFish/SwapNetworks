from typing import List
import numpy as np
from copy import deepcopy
from functools import reduce

from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
from qiskit.circuit import Parameter
from qiskit_aer.primitives import Estimator
from qiskit_aer.primitives.sampler import Sampler
# from qiskit_algorithms.optimizers import SPSA ,CG, SLSQP, L_BFGS_B, COBYLA
from qiskit_algorithms.optimizers import CG, SLSQP, L_BFGS_B, COBYLA, SPSA
# from qiskit.opflow.primitive_ops import PauliSumOp
from qiskit_nature.second_q.operators import FermionicOp
# from qiskit_nature.converters.second_quantization.qubit_converter import QubitConverter

import json
from copy import deepcopy
from Ternary_Tree.ucc.abstractucc import Molecule
from Ternary_Tree.qiskit_interface.circuit_provider import *
from Ternary_Tree.utils.circ_wrapper import ExcitationImpl, CircWrapper
from Ternary_Tree.utils.pauli import MajoranaContainer
from Ternary_Tree.utils.excitation import DoubleLadExcitation
from Ternary_Tree.utils.utils import lad2maj
from Ternary_Tree.utils.mapper import MajoranaMapper
from Ternary_Tree.qiskit_interface.circuit_provider import get_qiskit_device_noise_estimator
from my_utils import Timer
import multiprocessing as mp

import sys
from pyscf import lib

H2_4 = Molecule(geometry='H 0 0 0; H 0 0 0.7349', num_electrons=(1,1), active_orbitals=[0,1], basis='sto-3g')

def double_circ(mtoq):
    circ = QuantumCircuit(4)
    exc = DoubleLadExcitation((1, 0), (3, 2))
    maj_exc = lad2maj([exc])
    list_signs = ExcitationImpl.get_pauli_double_ex_yordan()
    for maj in maj_exc:
        pauli = reduce(lambda x,y: x * y, [mtoq[q] for q in maj.op])
        label = pauli.get_label_carr([0,1,2,3])
        list_signs[label] = (list_signs[label] * maj.sign * 1j**pauli.pow).imag
    ExcitationImpl.double_ex_short(circ, 0,1,2,3, Parameter("t"), list_signs)
    # ExcitationImpl.double_ex_yordan(circ, 0,1,2,3, Parameter("t"), list_signs)
    return circ

def num_observable(num_spin_orbitals):
    ferm_observable_N = FermionicOp(
        {f"+_{i} -_{i}": 1 for i in range(4)},
        num_spin_orbitals=4
    )
    return MajoranaMapper(mtoq).map(ferm_observable_N)

def compose(mtoq):
    init_circ = QuantumCircuit(4)
    ExcitationImpl.jw_init_state((0,2), init_circ, mtoq, encoding="xyz")
    # ExcitationImpl.jw_init_state((0,2), init_circ, mtoq, encoding="xyz")
    double = double_circ(mtoq)
    init_circ.compose(double, inplace=True, wrap=True)
    init_circ = transpile(init_circ, basis_gates=["u", "cx"])
    init_circ.assign_parameters([0.6 for i in init_circ.parameters], inplace=True)
    return init_circ

def disp(vs):
    return vs[1] - vs[0]**2

if __name__ == "__main__":
    init_state = Statevector.from_label("0000")
    up = UpGCCSD(H2_4)
    init_circ, mtoq = up.swap2xn(1, LadExcImpl.CNOT12xyz())
    # init_circ, mtoq = up.swap_gen(1, LadExcImpl.SHORT())
    init_circ.assign_parameters([-0.1117466680568747 for i in init_circ.parameters], inplace=True)
    mapper = MajoranaMapper(mtoq)
    # mtoq = MajoranaContainer.jw(4)
    n_obs = num_observable(4)
    n_obs2 = n_obs @ n_obs
    # mapper.map(fermionic_op)
    # # print(n_obs * n_obs)
    est = get_qiskit_device_noise_estimator("Z", 0.999, device="CPU")
    # init_circ = compose(mtoq)
    
    # print(init_circ)
    # fin_state = init_state.evolve(init_circ)
    # print(fin_state.probabilities())
    # sampler = Sampler(est)
    res = est.run([init_circ]*3, observables=[n_obs, n_obs2,  mapper.map(up.mol.hamiltonian.second_q_op())]).result()
    print(res.values)
    print(disp(res.values))