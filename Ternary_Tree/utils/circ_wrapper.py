from __future__ import annotations
 
from typing import Union, List, Dict
from abc import abstractmethod
from numpy import pi
from qiskit.circuit import Parameter, QuantumCircuit, ParameterExpression
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import Pauli as qiskit_pauli

from .pauli import Pauli
from .parameter import Parameter as MyParameter
from .utils import static_vars


class CircWrapper:
    @abstractmethod
    def mswap(self, pauli: Pauli):
        pass
    
    @abstractmethod
    def fswap(self, pauli: Pauli):
        pass

    @abstractmethod
    def rx(self, par: float | MyParameter):
        pass
    
    @abstractmethod
    def ry(self, par: float | MyParameter):
        pass

    @abstractmethod
    def rz(self, par: float | MyParameter):
        pass
    
    @abstractmethod
    def cz(self, par: float | MyParameter):
        pass
    
    @abstractmethod
    def cx(self, par: float | MyParameter):
        pass

    @abstractmethod
    def rotation(self, 
              pauli: Pauli,
              par: MyParameter,
              ) -> None:
        pass
    
    @abstractmethod
    def fswap(self, qubits):
        pass

    @abstractmethod
    def single_excitation(self, qubits, parameter: MyParameter, list_signs):
        pass

    @abstractmethod
    def get_pauli_double_ex_12cnot_alt():
        pass

    @abstractmethod
    def double_excitation(self, qubits, par, list_signs):
        pass

  
@static_vars(d_par={})
def to_qiskit_parameter(par: MyParameter) -> ParameterExpression:
    return to_qiskit_parameter.d_par.setdefault(par.name, par.coef * Parameter(par.name))


# def mswap(self, pauli: Pauli) -> None:
#     par = (1j**pauli.pow).imag * pi/4
#     label, qubits = pauli.get_label_qubs()
#     op = qiskit_pauli(label)
#     gate = PauliEvolutionGate(op, par)
#     gate.name = label
#     self.append(gate, qubits)
#     self.paulis.append((pauli, par))

class ExcitationImpl:
    
    @staticmethod
    def get_pauli_single_ex():
        return {"XY": 1, "YX": 1}

    @staticmethod
    def get_pauli_double_ex_short():
        return  {"XXXY": 1, "XXYX": 1, "XYXX": 1, "YXXX": 1, "YYYX": 1, "YYXY": 1, "YXYY": 1, "XYYY": 1}
    
    @staticmethod
    def get_pauli_double_ex_yordan():
        return {"XXXY": 1, "XXYX": 1, "XYXX": 1, "YXXX": 1, "YYYX": 1, "YYXY": 1, "YXYY": 1, "XYYY": 1}
            
    @staticmethod
    def get_pauli_double_ex_12cnot():
        return {"YYIZ": 1,"XXIZ": 1,"YYZI": 1,"XXZI": 1,"IZYY": 1,"IZXX": 1,"ZIYY": 1,"ZIXX": 1}
        
    @staticmethod
    def fswap(circ, q0, q1):
        circ.rx(pi/2, q0)
        circ.ry(pi/2, q1)
        circ.cx(1,0)
        circ.rz(-pi/2, q0)
        circ.ry(pi/2, q1)
        circ.cx(1,0)
        circ.rx(pi/2, q0)
        circ.ry(pi/2, q1)
        circ.z(q1)

    @staticmethod
    def single_ex(circ: CircWrapper, q0, q1, parameter, list_signs):
        circ.ry(pi/2, q0)
        circ.cx(q0, q1)
        circ.ry(-parameter * list_signs["XY"], q0)
        circ.ry(parameter * list_signs["YX"], q1)
        circ.cx(q0, q1)
        circ.ry(-pi/2, q0)

    @staticmethod
    def single_ex_alt(circ: CircWrapper, q0, q1, parameter, list_signs):
        circ.ry(pi/2, q0)
        circ.cz(q0, q1)
        circ.ry(-parameter*list_signs["ZY"], q0)
        circ.ry(parameter*list_signs["YZ"], q1)
        circ.cz(q0, q1)
        circ.ry(-pi/2, q0)

    @staticmethod
    def double_ex_short(circ,
                        q0, q1, q2, q3,
                        parameter,
                        list_signs):
        circ.cx(q1, q0)
        circ.cx(q2, q3)
        circ.ry(pi/2, q2)
        circ.cx(q2, q1)
        circ.ry( parameter*list_signs["XXXY"]*0.25, q2)
        circ.ry(-parameter*list_signs["YXXX"]*0.25, q1)
        circ.cx(q0, q1)
        circ.cx(q3, q2)
        circ.ry( parameter*list_signs["YXYY"]*0.25, q1)
        circ.ry(-parameter*list_signs["YYXY"]*0.25, q2)
        circ.cx(q0, q2)
        circ.cx(q3, q1)
        circ.ry( parameter*list_signs["XYYY"]*0.25, q1)
        circ.ry(-parameter*list_signs["YYYX"]*0.25, q2)
        circ.cx(q0, q1)
        circ.cx(q3, q2)
        circ.ry(-parameter*list_signs["XYXX"]*0.25, q1)
        circ.ry( parameter*list_signs["XXYX"]*0.25, q2)
        circ.cx(q0, q2)
        circ.cx(q3, q1)
        circ.cx(q2, q1)
        circ.ry(-pi/2, q2)
        circ.cx(q1, q0)
        circ.cx(q2, q3)

    @staticmethod
    def double_ex_alt_opt(circ,
                        q0, q1, q2, q3,
                        parameter,
                        list_signs):
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

    @staticmethod
    def double_ex_yordan(circ,
                        q0, q1, q2, q3,
                        parameter,
                        list_signs):
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.x(q1), circ.x(q3)
        circ.cx(q0, q2)
        circ.h(q1), circ.ry(list_signs["YXXX"]*parameter*0.25, q0)
        circ.cx(q0, q1)
        circ.h(q3), circ.ry(-list_signs["XYXX"]*parameter*0.25, q0)
        circ.cx(q0, q3)
        circ.ry(-list_signs["XYYY"]*parameter*0.25, q0)
        circ.cx(q0, q1)
        circ.h(q2), circ.ry(list_signs["YXYY"]*parameter*0.25, q0)
        circ.cx(q0, q2)
        circ.ry(-list_signs["XXXY"]*parameter*0.25, q0)
        circ.cx(q0, q1)
        circ.ry(-list_signs["YYXY"]*parameter*0.25, q0)
        circ.cx(q0, q3)
        circ.h(q3), circ.ry(list_signs["YYYX"]*parameter*0.25, q0)
        circ.cx(q0, q1)
        circ.h(q1), circ.ry(list_signs["XXYX"]*parameter*0.25, q0)
        circ.rz(pi/2, q2)
        circ.cx(q0, q2)
        circ.rz(pi/2, q2), circ.rz(-pi/2, q0)
        circ.ry(pi/2, q2)
        circ.x(q1), circ.x(q3)
        circ.cx(q0, q1), circ.cx(q2, q3)

    @staticmethod
    def double_ex_12cnot(circ,
                        q0, q1, q2, q3,
                        parameter,
                        list_signs):
        circ.ry(pi/2, q2)
        circ.rz(pi/2, q3)
        circ.rz(pi/2, q0)
        circ.ry(pi/2, q0)
        circ.cx(q1, q2)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.ry( parameter*list_signs["XXZI"]*0.25, q0)
        circ.ry( parameter*list_signs["YYZI"]*0.25, q1)
        circ.ry(-parameter*list_signs["IZYY"]*0.25, q2)
        circ.ry(-parameter*list_signs["IZXX"]*0.25, q3)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.cx(q1, q2), circ.cx(q3, q0)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.ry( parameter*list_signs["XXIZ"]*0.25, q0)
        circ.ry( parameter*list_signs["YYIZ"]*0.25, q1)
        circ.ry(-parameter*list_signs["ZIYY"]*0.25, q2)
        circ.ry(-parameter*list_signs["ZIXX"]*0.25, q3)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.cx(q3, q0)
        circ.ry(-pi/2, q2)
        circ.ry(-pi/2, q0)
        circ.rz(-pi/2, q3)
        circ.rz(-pi/2, q0)

    @staticmethod
    def double_ex_12cnot_alt(circ,
                            q0, q1, q2, q3,
                            parameter,
                            list_signs):
        circ.rz(pi/2, q0), circ.ry(pi/2, q1), circ.ry(pi/2, q2), circ.rz(pi/2, q3)
        circ.cz(q1, q2)
        circ.rx(pi/2, q1), circ.rx(pi/2, q2)
        circ.cz(q0, q1), circ.cz(q2, q3)
        circ.rx(-parameter*list_signs["ZZXI"]*0.25, q1)
        circ.rx( parameter*list_signs["YYXI"]*0.25, q0)
        circ.rx( parameter*list_signs["IXYY"]*0.25, q3)
        circ.rx(-parameter*list_signs["IXZZ"]*0.25, q2)
        circ.cz(q0, q1), circ.cz(q2, q3)
        circ.rx(pi/2, q0), circ.rx(-pi/2, q1), circ.rx(-pi/2, q2), circ.rx(pi/2, q3) 
        circ.cz(q1, q2), circ.cz(q3, q0)
        circ.rx(pi/2, q0), circ.rx(pi/2, q1), circ.rx(pi/2, q2), circ.rx(pi/2, q3) 
        circ.cz(q0, q1), circ.cz(q2, q3)
        circ.rx(-parameter*list_signs["ZZIX"]*0.25, q1)
        circ.rx(-parameter*list_signs["YYIX"]*0.25, q0)
        circ.rx(-parameter*list_signs["XIYY"]*0.25, q3)
        circ.rx(-parameter*list_signs["XIZZ"]*0.25, q2)
        circ.cz(q0, q1), circ.cz(q2, q3)
        circ.rx(-pi/2, q0), circ.rx(-pi/2, q3)
        circ.cz(q3, q0)
        circ.rx(-pi/2, q0), circ.rx(-pi/2, q1), circ.rx(-pi/2, q2), circ.rx(-pi/2, q3) 
        circ.rz(-pi/2, q0), circ.ry(-pi/2, q1), circ.ry(-pi/2, q2), circ.rz(-pi/2, q3)