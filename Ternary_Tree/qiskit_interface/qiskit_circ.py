from abc import abstractmethod
from numpy import pi
from qiskit.circuit import Parameter, QuantumCircuit, ParameterExpression
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import Pauli as qiskit_pauli
from typing import List, Dict
from ..utils import Parameter as MyParameter
from ..utils import static_vars, Pauli, CircWrapper

  
@static_vars(d_par={})
def to_qiskit_parameter(par: MyParameter) -> ParameterExpression:
    return to_qiskit_parameter.d_par.setdefault(par.name, par.coef * Parameter(par.name))


class QiskitCirc(CircWrapper, QuantumCircuit):
    def __init__(self, qc, name=None):
        self.paulis = []
        self.swaps = []
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


    @staticmethod
    def get_pauli_single_ex():
        return {"XY": -1, "YX": 1}

    def single_ex(self, qubits, par: MyParameter, list_signs):
        parameter =  to_qiskit_parameter(par)
        q0, q1 = 0,1
        circ = QuantumCircuit(2)
        circ.ry(pi/2, q0)
        circ.cx(q0,q1)
        circ.ry(parameter*list_signs["XY"], q0)
        circ.ry(parameter*list_signs["YX"], q1)
        circ.cx(q0,q1)
        circ.ry(-pi/2, q0)
        circ.name = par.name
        self.compose(circ, qubits=qubits, inplace=True, wrap=True)

    @staticmethod
    def get_pauli_double_ex_short():
        return  {"XXXY": 1, "XXYX": 1, "XYXX": -1, "YXXX": -1, "YYYX": -1, "YYXY": -1, "YXYY": 1, "XYYY": 1}
        
    def double_ex_short(self, 
                        qubits: List[int], 
                        par: MyParameter, 
                        list_signs: Dict[str, int]):
        parameter = to_qiskit_parameter(par)
        q0, q1, q2, q3 = 0,1,2,3
        circ = QuantumCircuit(4)
        circ.cx(q1, q0)
        circ.cx(q2, q3)
        circ.ry(pi/2, q2)
        circ.cx(q2, q1)
        circ.ry(parameter*list_signs["XXXY"]*0.25, q2)
        circ.ry(parameter*list_signs["YXXX"]*0.25, q1)
        circ.cx(q0, q1)
        circ.cx(q3, q2)
        circ.ry(parameter*list_signs["YXYY"]*0.25, q1)
        circ.ry(parameter*list_signs["YYXY"]*0.25, q2)
        circ.cx(q0, q2)
        circ.cx(q3, q1)
        circ.ry(parameter*list_signs["XYYY"]*0.25, q1)
        circ.ry(parameter*list_signs["YYYX"]*0.25, q2)
        circ.cx(q0, q1)
        circ.cx(q3, q2)
        circ.ry(parameter*list_signs["XYXX"]*0.25, q1)
        circ.ry(parameter*list_signs["XXYX"]*0.25, q2)
        circ.cx(q0, q2)
        circ.cx(q3, q1)
        circ.cx(q2, q1)
        circ.ry(-pi/2, q2)
        circ.cx(q1, q0)
        circ.cx(q2, q3)
        circ.name = par.name

        self.compose(circ, qubits=qubits, inplace=True, wrap=True)
    
    def double_ex_zyx_opt(self, qubits, par, list_signs):
        parameter =  to_qiskit_parameter(par)
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
        circ.name = par.name

        self.compose(circ, qubits=qubits, inplace=True, wrap=True)

    @staticmethod
    def get_pauli_double_ex_yordan():
        return {"XXXY": -1, "XXYX": 1, "XYXX": -1, "YXXX": 1, "YYYX": 1, "YYXY": -1, "YXYY": 1, "XYYY": -1}
        
    def double_ex_yordan(self, qubits, par, list_signs):
        parameter =  to_qiskit_parameter(par)
        circ = QuantumCircuit(4)
        q0, q1, q2, q3 = 3,2,1,0
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.x(q1), circ.x(q3)
        circ.cx(q0, q2)
        circ.h(q1), circ.ry(list_signs["YXXX"]*parameter*0.25, q0)
        circ.cx(q0, q1)
        circ.h(q3), circ.ry(list_signs["XYXX"]*parameter*0.25, q0)
        circ.cx(q0, q3)
        circ.ry(list_signs["XYYY"]*parameter*0.25, q0)
        circ.cx(q0, q1)
        circ.h(q2), circ.ry(list_signs["YXYY"]*parameter*0.25, q0)
        circ.cx(q0, q2)
        circ.ry(list_signs["XXXY"]*parameter*0.25, q0)
        circ.cx(q0, q1)
        circ.ry(list_signs["YYXY"]*parameter*0.25, q0)
        circ.cx(q0, q3)
        circ.h(q3), circ.ry(list_signs["YYYX"]*parameter*0.25, q0)
        circ.cx(q0, q1)
        circ.h(q1), circ.ry(list_signs["XXYX"]*parameter*0.25, q0)
        circ.rz(pi/2, q2)
        circ.cx(q0, q2)
        # circ.h(q2)
        # circ.cx(q0, q2)
        circ.rz(pi/2, q2), circ.rz(-pi/2, q0)
        circ.ry(pi/2, q2)
        circ.x(q1), circ.x(q3)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.name = par.name
        self.compose(circ, qubits=qubits, inplace=True, wrap=True)

    def double_ex_zyx_yordan(self, qubits, par, list_signs):
        parameter =  to_qiskit_parameter(par)
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
        circ.name = par.name
        self.compose(circ, qubits=qubits, inplace=True, wrap=True)
    
    @staticmethod
    def get_pauli_double_ex_12cnot():
        return {"YYIZ": 1,"XXIZ": 1,"YYZI": 1,"XXZI": 1,"IZYY": -1,"IZXX": -1,"ZIYY": -1,"ZIXX": -1}

    # TODO
    def double_ex_12cnot(self, qubits, par, list_signs):
        
        parameter = to_qiskit_parameter(par)
        circ = QuantumCircuit(4)
        q0, q1, q2, q3 = 0,1,2,3
        circ.ry(pi/2, q2)
        circ.rz(pi/2, q3)
        circ.rz(pi/2, q0)
        circ.ry(pi/2, q0)
        circ.cx(q1, q2)
        circ.cx(q0, q1), circ.cx(q2, q3)
        # list_signs = {"YYIZ": 1,"XXIZ": 1,"YYZI": 1,"XXZI": 1,"IZYY": -1,"IZXX": -1,"ZIYY": -1,"ZIXX": -1}

        circ.ry(parameter*list_signs["XXZI"]*0.25, q0)
        circ.ry(parameter*list_signs["YYZI"]*0.25, q1)
        circ.ry(parameter*list_signs["IZYY"]*0.25, q2)
        circ.ry(parameter*list_signs["IZXX"]*0.25, q3)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.cx(q1, q2), circ.cx(q3, q0)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.ry(parameter*list_signs["XXIZ"]*0.25, q0)
        circ.ry(parameter*list_signs["YYIZ"]*0.25, q1)
        circ.ry(parameter*list_signs["ZIYY"]*0.25, q2)
        circ.ry(parameter*list_signs["ZIXX"]*0.25, q3)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.cx(q3, q0)
        circ.ry(-pi/2, q2)
        circ.ry(-pi/2, q0)
        circ.rz(-pi/2, q3)
        circ.rz(-pi/2, q0)
        circ.name = par.name
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

        self.compose(circ, qubits=qubits, inplace=True, wrap=True)