from typing import Dict

from qiskit.circuit import Parameter, QuantumCircuit, ParameterExpression, QuantumRegister
from numpy import pi
from qiskit.circuit.library import HamiltonianGate, TwoLocal, n_local
from qiskit.quantum_info import Operator
from .paulistring import PauliStringSign

class MyCirq(QuantumCircuit):

    def __init__(self, qc, name=None):
        self.paulis = []
        self.swaps = []
        super().__init__(qc, name=name)

    def rotation(self, 
              pss: PauliStringSign,
              par: Parameter,
              ):
        par = pss.sign*par*(-1j)
        pauli = list(pss.ps.values())
        els = sorted(list(pss.ps.items()), key=lambda x: x[0]) 
        pauli = [el[1] for el in els]
        qubits = [el[0] for el in els]
        if len(qubits) <=2 and isinstance(par, (ParameterExpression, Parameter)):
            getattr(self,"r" + pauli[0].lower() * len(qubits))(par, *qubits)


    def maj_swap(self, pss: PauliStringSign):
        par = pss.sign*pi/4
        pauli = list(pss.ps.values())
        els = sorted(list(pss.ps.items()), key=lambda x: x[0]) 
        pauli = [el[1] for el in els]
        qubits = list(reversed([el[0] for el in els]))
        # qubits = [el[0] for el in els]
        label = ''.join(pauli)
        op = Operator.from_label(label)
        gate = HamiltonianGate(op, par)
        gate.name = label
        self.append(gate, qubits)
        self.paulis.append((pauli, par))

    def fswap(self, qubits):
        q0, q1 = qubits
        circ = QuantumCircuit(2)
        circ.rx(pi/2, 0)
        circ.ry(pi/2, 1)
        circ.cx(1,0)
        circ.rz(-pi/2, 0)
        circ.ry(pi/2, 1)
        circ.cx(1,0)
        circ.rx(pi/2, 0)
        circ.ry(pi/2, 1)
        circ.z(1)
        circ.name = "fswap"
        self.compose(circ, qubits=qubits, inplace=True, wrap=True)

    def double_ex_opt(self, qubits, parameter, list_signs):
        q0, q1, q2, q3 = 1,0,2,3
        circ = QuantumCircuit(4)
        circ.cx(q1, q0)
        circ.cx(q2, q3)
        circ.ry(pi/2, q2)
        circ.cx(q2, q1)
        circ.ry(parameter*list_signs["XXXY"]*0.25, q2)
        circ.ry(-parameter*list_signs["YXXX"]*0.25, q1)
        circ.cx(q0, q1)
        circ.cx(q3, q2)
        circ.ry(parameter*list_signs["YXYY"]*0.25, q1)
        circ.ry(-parameter*list_signs["YYXY"]*0.25, q2)
        circ.cx(q0, q2)
        circ.cx(q3, q1)
        circ.ry(parameter*list_signs["XYYY"]*0.25, q1)
        circ.ry(-parameter*list_signs["YYYX"]*0.25, q2)
        circ.cx(q0, q1)
        circ.cx(q3, q2)
        circ.ry(-parameter*list_signs["XYXX"]*0.25, q1)
        circ.ry(parameter*list_signs["XXYX"]*0.25, q2)
        circ.cx(q0, q2)
        circ.cx(q3, q1)
        circ.cx(q2, q1)
        circ.ry(-pi/2, q2)
        circ.cx(q1, q0)
        circ.cx(q2, q3)
        circ.name = parameter.name

        self.compose(circ, qubits=qubits, inplace=True, wrap=True)
    
    def double_ex_zyx_opt(self, qubits, parameter, list_signs):
        q0, q1, q2, q3 = 1,0,2,3
        circ = QuantumCircuit(4)
        circ.cx(q0, q1)
        circ.cx(q3, q2)
        circ.ry(pi/2, q2)
        circ.cx(q1, q2)
        circ.ry(parameter*list_signs["ZZZY"]*0.25, q2)
        circ.ry(-parameter*list_signs["YZZZ"]*0.25, q1)
        circ.cx(q1, q0)
        circ.cx(q2, q3)
        circ.ry(parameter*list_signs["YZYY"]*0.25, q1)
        circ.ry(-parameter*list_signs["YYZY"]*0.25, q2)
        circ.cx(q2, q0)
        circ.cx(q1, q3)
        circ.ry(parameter*list_signs["ZYYY"]*0.25, q1)
        circ.ry(-parameter*list_signs["YYYZ"]*0.25, q2)
        circ.cx(q1, q0)
        circ.cx(q2, q3)
        circ.ry(-parameter*list_signs["ZYZZ"]*0.25, q1)
        circ.ry(parameter*list_signs["ZZYZ"]*0.25, q2)
        circ.cx(q2, q0)
        circ.cx(q1, q3)
        circ.cx(q1, q2)
        circ.ry(-pi/2, q2)
        circ.cx(q0, q1)
        circ.cx(q3, q2)
        circ.name = parameter.name

        self.compose(circ, qubits=qubits, inplace=True, wrap=True)
        
    def double_ex_yordan(self, qubits, parameter, list_signs):
        circ = QuantumCircuit(4)
        q0, q1, q2, q3 = 0,1,2,3
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.x(q1), circ.x(q3)
        circ.cx(q0, q2)
        circ.h(q1), circ.ry(parameter*0.25, q0)
        circ.cx(q0, q1)
        circ.h(q3), circ.ry(-parameter*0.25, q0)
        circ.cx(q0, q3)
        circ.ry(parameter*0.25, q0)
        circ.cx(q0, q1)
        circ.h(q2), circ.ry(-parameter*0.25, q0)
        circ.cx(q0, q2)
        circ.ry(parameter*0.25, q0)
        circ.cx(q0, q1)
        circ.ry(-parameter*0.25, q0)
        circ.cx(q0, q3)
        circ.h(q3), circ.ry(parameter*0.25, q0)
        circ.cx(q0, q1)
        circ.h(q1), circ.ry(-parameter*0.25, q0)
        circ.rz(pi/2, q2)
        circ.cx(q0, q2)
        # circ.h(q2)
        # circ.cx(q0, q2)
        circ.rz(pi/2, q2), circ.rz(-pi/2, q0)
        circ.ry(pi/2, q2)
        circ.x(q1), circ.x(q3)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.name = parameter.name
        self.compose(circ, qubits=qubits, inplace=True, wrap=True)

    def double_ex_zyx_yordan(self, qubits, parameter, list_signs):
        circ = QuantumCircuit(4)
        q0, q1, q2, q3 = 0,1,2,3
        circ.cx(q1, q0), circ.cx(q3, q2)
        circ.x(q1), circ.x(q3)
        circ.cx(q2, q0)
        circ.h(q1), circ.ry(parameter*0.25*list_signs["YZZZ"], q0)
        circ.cx(q1, q0)
        circ.h(q3), circ.ry(-parameter*0.25*list_signs["ZYZZ"], q0)
        circ.cx(q3, q0)
        circ.ry(parameter*0.25*list_signs["ZYYY"], q0)
        circ.cx(q1, q0)
        circ.h(q2), circ.ry(-parameter*0.25*list_signs["YZYY"], q0)
        circ.cx(q2, q0)
        circ.ry(parameter*0.25*list_signs["ZZZY"], q0)
        circ.cx(q1, q0)
        circ.ry(-parameter*0.25*list_signs["YYZY"], q0)
        circ.cx(q3, q0)
        circ.h(q3), circ.ry(parameter*0.25*list_signs["YYYZ"], q0)
        circ.cx(q1, q0)
        circ.h(q1), circ.ry(-parameter*0.25*list_signs["ZZYZ"], q0)
        circ.rz(pi/2, q2)
        circ.cx(q2, q0)
        circ.h(q2)
        circ.cx(q2, q0)
        # circ.rz(pi/2, q2), circ.rz(-pi/2, q0)
        # circ.ry(pi/2, q2)
        circ.x(q1), circ.x(q3)
        circ.cx(q1, q0), circ.cx(q3, q2)
        circ.name = parameter.name
        self.compose(circ, qubits=qubits, inplace=True, wrap=True)
        
        
    def fermionic_2swap(self, qubits):
        """
        Implements the 2-fermionic swap gate.

        Args:
            qubits (list[int]): The indices of the qubits involved in the swap.
        """
        # Define the qubit indices
        q0, q1, q2, q3 = 0,1,2,3

        # Define a quantum circuit
        circ = QuantumCircuit(4)

        # Apply the necessary gates
        circ.ry(pi/2, q1)
        circ.cx(q2, q1)
        circ.ry(pi/2, q2)
        circ.rx(pi/2, q0)
        circ.rx(pi/2, q3)
        circ.cx(q2, q0)
        circ.cx(q1, q3)
        circ.rz(-pi/2, q0)
        circ.ry(pi/2, q1)
        circ.ry(pi/2, q2)
        circ.rz(-pi/2, q3)
        circ.cx(q2, q0)
        circ.cx(q1, q3)
        circ.rx(pi/2, q0)
        circ.ry(-pi/2, q2)
        circ.rx(-pi/2, q3)
        circ.cx(q2, q1)
        circ.ry(pi/2, q1)
        circ.y(q2)
        circ.x(q3)
        circ.name = "2fswap"

        # Combine the circuit with the current circuit
        self.compose(circ, qubits=qubits, inplace=True, wrap=True)

