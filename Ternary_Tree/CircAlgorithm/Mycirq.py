from typing import Dict

from qiskit.circuit import Parameter, QuantumCircuit, ParameterExpression, QuantumRegister
from numpy import pi
from qiskit.circuit.library import HamiltonianGate, TwoLocal, n_local
from qiskit.quantum_info import Operator

QuantumCircuit()
class MyCirq(QuantumCircuit):

    def __init__(self, qc, name=None):
        self.paulis = []
        self.swaps = []
        super().__init__(qc, name=name)

    def pauli(self, 
              pl: Dict[int, str],
              par: float | Parameter = pi/4
              ):
        
        if not isinstance(pl, dict):
            pl = {gate[0]: gate[1] for gate in pl}
        pauli = {}
        qubits = []
        for el in pl:
            if pl[el] != "I":
                pauli[el] = pl[el]
                qubits.append(el)
        qubits = list(reversed(qubits))
        if (len(qubits) == 2) and (isinstance(par, Parameter) or isinstance(par, ParameterExpression)):
            self.rzz(par, qubits[0], qubits[1])
            # gate.name = '0'

        if (len(qubits) == 1)  and (isinstance(par, Parameter) or isinstance(par, ParameterExpression)):
            self.rz(par, qubits[0])

        if isinstance(par, float):
            label = ''.join(pauli.values())
            op = Operator.from_label(label)
            qubits = list(reversed([i for i in pauli.keys()]))
            gate = HamiltonianGate(op, par)

            # gate.name = 'pi/4'
            gate.name = "%s" %label
            self.append(gate, qubits)
            self.paulis.append((pauli, par))



    def pauli_str(self, 
                  pl: Dict[int, str]
                  ):
        if not isinstance(pl, dict):
            pl = {gate[0]: gate[1] for gate in pl}
        pauli = {}
        for el in pl:
            if pl[el] != "I":
                pauli[el] = pl[el]

        for qubit, gate in pauli.items():
            if gate == "X":
                self.x(qubit)
            if gate == "Y":
                self.y(qubit)
            if gate == "Z":
                self.z(qubit)
