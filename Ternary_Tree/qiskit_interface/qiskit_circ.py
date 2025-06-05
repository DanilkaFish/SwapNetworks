from abc import abstractmethod
from numpy import pi
import logging

from qiskit.circuit import Parameter, QuantumCircuit, ParameterExpression, Instruction
from qiskit.circuit.library import PauliEvolutionGate, IGate
from qiskit.quantum_info import Pauli as qiskit_pauli
from typing import List, Dict
from ..utils import Parameter as MyParameter
from ..utils import static_vars, Pauli, CircWrapper
from ..utils.circ_wrapper import *

logger = logging.getLogger(__name__)


@static_vars(d_par={})
def to_qiskit_parameter(par: MyParameter) -> ParameterExpression:
    return to_qiskit_parameter.d_par.setdefault(par.name, par.coef * Parameter(par.name))


class QiskitCirc(CircWrapper, QuantumCircuit):
    def __init__(self, qc, name=None):
        self.paulis = []
        self.swaps = []
        self.excitation_pos = dict()
        super().__init__(qc, name=name)

    def rotation(self, 
              pauli: Pauli,
              par: MyParameter,
              ) -> None:
        par = to_qiskit_parameter(par)
        par = (1j**pauli.pow).imag * par
        label, qubits = pauli.get_label_qubs()
        if len(qubits) <=2 and isinstance(par, (ParameterExpression, Parameter)):
            getattr(self,"r" + label[0].lower() * len(qubits))(par, *qubits)
            # if list(par.parameters)[0] not in self.excitation_pos:
            #     self.excitation_pos[list(par.parameters)[0]] = [len(self) - 1]
            # else:
            #     self.excitation_pos[list(par.parameters)[0]].append(len(self) - 1)

    def any_rot(self, 
              pauli: Pauli,
              par: MyParameter,
              ) -> None:
        par = to_qiskit_parameter(par)
        label, qubits = pauli.get_label_qubs()
        op = qiskit_pauli(label)
        gate = PauliEvolutionGate(op, par)
        gate.name = label
        self.append(gate, qubits)
        self.paulis.append((pauli, par))
        
    def mswap(self, pauli: Pauli) -> None:
        par = (1j**pauli.pow).imag * pi/4
        label, qubits = pauli.get_label_qubs()
        op = qiskit_pauli(label)
        gate = PauliEvolutionGate(op, par)
        gate.name = label
        self.append(gate, qubits)
        self.paulis.append((pauli, par))

    def fswap(self, qubits):
        q0, q1 = qubits
        circ = QuantumCircuit(2)
        circ.rx(pi/2, q0)
        circ.ry(pi/2, q1)
        circ.cx(1,0)
        circ.rz(-pi/2, q0)
        circ.ry(pi/2, q1)
        circ.cx(1,0)
        circ.rx(pi/2, q0)
        circ.ry(pi/2, q1)
        circ.z(1)
        circ.name = "fswap"
        self.compose(circ, qubits=qubits, inplace=True, wrap=True)

    def single_excitation(self, method_name, qubits, par: MyParameter, list_signs):
        parameter =  to_qiskit_parameter(par)
        q0, q1 = 0,1
        circ = QuantumCircuit(2)
        getattr(CircWrapper, method_name)(circ, q0, q1, parameter, list_signs)
        self.compose(circ, qubits=qubits, inplace=True, wrap=True)
        # if list(circ.parameters)[0] not in self.excitation_pos:
        #     self.excitation_pos[list(circ.parameters)[0]] = [len(self) - 1]
        # else:
        #     self.excitation_pos[list(circ.parameters)[0]].append(len(self) - 1)

    def double_excitation(self, method_name, qubits, par: MyParameter, list_signs):
        parameter = to_qiskit_parameter(par)
        q0, q1, q2, q3 = 0,1,2,3
        circ = QuantumCircuit(4)
        getattr(CircWrapper, method_name)(circ, q0, q1, q2, q3, parameter, list_signs)
        circ.name = par.name
        self.compose(circ, qubits=qubits, inplace=True, wrap=True)
        self.excitation_pos[circ.parameters[0]] = [len(self) - 1]


def to_excitaions(qc, parameters, excitation_pos) -> QuantumCircuit:
    circ = qc.copy()
    for parameter in excitation_pos:
        if parameter not in parameters:
            for item in excitation_pos[parameter]:
                for qub in circ[item].qubits[0:1]:
                    circ.data[item] = (IGate(), [qub], [])
    return circ