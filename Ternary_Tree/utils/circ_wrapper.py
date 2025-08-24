from __future__ import annotations
 
from typing import Union, List, Dict
import logging
from abc import abstractmethod


from numpy import pi
import numpy as np
import scipy


from qiskit.circuit import Gate
from qiskit.circuit import QuantumCircuit

def iswap(theta, phi):
    # exp(-i θ/2 (cosφ XX + sinφ YY))
    xx = np.kron([[0,1],[1,0]], [[0,1],[1,0]])
    yy = np.kron([[0,-1j],[1j,0]], [[0,-1j],[1j,0]])
    h = np.cos(phi)*xx + np.sin(phi)*yy
    return scipy.linalg.expm(-1j * theta/2 * h)

class MS2Gate(Gate):
    def __init__(self, theta, phi):
        super().__init__("ms2", 2, [theta, phi])

    def _define(self):
        return None

    def to_matrix(self):
        return iswap(*self.params)
    

class LadExcImpl:
    @staticmethod
    def YORDAN():
        return "double_ex_yordan"    

    @staticmethod
    def CNOT12():
        return "double_ex_12cnot" 

    @staticmethod
    def SINGLE():
        return "single_ex"
    
    @staticmethod
    def SHORT():
        return "double_ex_short"
     

    
from .utils import static_vars
from .utils import Parameter as MyParameter
from .pauli import Pauli, MajoranaContainer


logger = logging.getLogger(__name__)

class CircWrapper:
    @property
    def num_qubits(self):
        pass

    @abstractmethod
    def mswap(self, pauli: Pauli):
        pass
    
    @abstractmethod
    def fswap(self, pauli: Pauli):
        pass

    @abstractmethod
    def rx(self, par: float | MyParameter, qub):
        pass
    
    @abstractmethod
    def ry(self, par: float | MyParameter, qub):
        pass

    @abstractmethod
    def rz(self, par: float | MyParameter, qub):
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
    def single_excitation(self, qubits, parameter, list_signs):
        pass

    @abstractmethod
    def get_pauli_double_ex_12cnot_alt():
        pass

    @abstractmethod
    def double_excitation(self, qubits, par, list_signs):
        pass

  


def cx_via_bridge(circ, a, b, c):
    """Реализует CX(a -> c) через промежуточный b на линейной связи a-b-c."""
    # Шаблон 4 CX (оптимальнее, чем 2 SWAP = 6 CX)
    circ.cx(a, b)
    circ.cx(b, c)
    circ.cx(a, b)
    circ.cx(b, c)

class ExcitationImpl:

    @staticmethod
    def jw_init_state(electrons, circ: CircWrapper, mtoq: MajoranaContainer, encoding="xyz"):
        n = circ.num_qubits
        if encoding[-1] == "z":
            prep_method = lambda qub, angle: circ.rx(pi + angle, qub)
        elif encoding[-1] == "y":
            prep_method = lambda qub, angle: circ.ry(pi/2 + angle, qub)
        elif encoding[-1] == "x":
            prep_method = lambda qub, angle: circ.rz(pi/2 + angle, qub)
        if encoding[0:2] in {"xy", "yz", "zx"}:
            angle = 0
        else:
            angle = pi

        for i in range(n):
            if mtoq.get_by_qubs(i)[0]/2 in electrons:
                prep_method(i, angle)
            else:
                prep_method(i, angle + pi)

    @staticmethod
    def majorana_init_state(electrons: List[int], circ: CircWrapper, mtoq: MajoranaContainer, encoding: str="xyz"):
        n = circ.num_qubits
        if encoding[-1] == "z":
            prep_method = lambda qub, angle: circ.rx(pi + angle, qub)
        elif encoding[-1] == "y":
            prep_method = lambda qub, angle: circ.ry(pi/2 + angle, qub)
        elif encoding[-1] == "x":
            prep_method = lambda qub, angle: circ.rz(pi/2 + angle, qub)
        if encoding[0:2] in {"xy", "yz", "zx"}:
            angle = 0
        else:
            angle = pi

        for i in range(n):
            if mtoq.get_by_qubs(i)[0]/2 in electrons:
                prep_method(i, angle)
            else:
                prep_method(i, angle + pi)
        for i in range(2):
            for j in range(i*n//2, i*n//2 + n//2 - 1, 2):
                circ.mswap(mtoq.transpose(j, 1, j + 1, 0))
        logger.info(f"preparing majorana init state for {encoding=} and {electrons=}...")
        
    @staticmethod
    def get_pauli_single_ex():
        return {"XY": 1, "YX": 1}
            
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
    def get_pauli_double_ex_short():
        return {'XXXY': 1, 'XXYX': 1, 'XYXX': -1, 'YXXX': -1, 'YYYX': -1, 'YYXY': -1, 'YXYY': 1, 'XYYY': 1}
    
    @staticmethod
    def double_ex_short(circ,
                        q0, q1, q2, q3,
                        parameter,
                        list_signs):
        """
            Implementation of exp(-parameter/8("XXXY + XXYX + XYXX + YXXX + YYYX + YYXY + YXYY + XYYY))
        """
        circ.cx(q1, q0)
        circ.cx(q2, q3)
        circ.ry(pi/2, q2)
        circ.cx(q2, q1)
        circ.ry(parameter*list_signs["XYXX"]*0.25, q1)
        circ.ry( parameter*list_signs["XXYX"]*0.25, q2)
    
        circ.cx(q0, q1)
        circ.cx(q3, q2)
        circ.ry( parameter*list_signs["XXXY"]*0.25, q2)
        circ.ry(parameter*list_signs["YXXX"]*0.25, q1)
        
        circ.cx(q0, q2)
        circ.cx(q3, q1)
        circ.ry(parameter*list_signs["YYXY"]*0.25, q2)
        circ.ry( parameter*list_signs["YXYY"]*0.25, q1)
        
        circ.cx(q0, q1)
        circ.cx(q3, q2)
        circ.ry( parameter*list_signs["XYYY"]*0.25, q1)
        circ.ry(parameter*list_signs["YYYX"]*0.25, q2)
        circ.cx(q0, q2)
        circ.cx(q3, q1)
        circ.cx(q2, q1)
        circ.ry(-pi/2, q2)
        circ.cx(q1, q0)
        circ.cx(q2, q3)


    @staticmethod
    def get_pauli_double_linear_short():
        return {'XXXY': 1, 'XXYX': 1, 'XYXX': -1, 'YXXX': -1, 'YYYX': 1, 'YYXY': 1, 'YXYY': 1, 'XYYY': 1}
    
    @staticmethod
    def double_linear_short(circ,
                            q0, q1, q2, q3,
                            parameter,
                            list_signs):
        circ.cx(q1, q0)
        circ.cx(q2, q3)
        circ.ry(pi/2, q2)
        circ.cx(q2, q1)
        circ.ry(parameter*list_signs["XYXX"]*0.25, q1)
        circ.ry( parameter*list_signs["XXYX"]*0.25, q2)
        circ.cx(q0, q1)
        circ.cx(q3, q2)
        circ.ry( parameter*list_signs["XXXY"]*0.25, q2)
        circ.ry(parameter*list_signs["YXXX"]*0.25, q1)
        circ.cx(q2, q1)
        circ.ry(pi/2, q1)
        circ.cx(q0, q1)
        circ.cx(q3, q2)
        circ.ry(pi/2, q2)
        circ.cx(q1, q2)
        circ.ry(parameter*list_signs["YYYX"]*0.25, q2)
        circ.ry( parameter*list_signs["XYYY"]*0.25, q1)
        circ.cx(q0, q1)
        circ.cx(q3, q2)
        circ.ry( parameter*list_signs["YYXY"]*0.25, q2)
        circ.ry(parameter*list_signs["YXYY"]*0.25, q1)
    
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
    def get_pauli_double_ex_yordan():
        return {"XXXY": -1, "XXYX": 1, "XYXX": -1, "YXXX": 1, "YYYX": 1, "YYXY": -1, "YXYY": 1, "XYYY": -1}
    
    @staticmethod
    def double_ex_yordan(circ,
                        q0, q1, q2, q3,
                        parameter,
                        list_signs):
        
        q0,q1 = q1,q0
        nlist = {}
        olist = ExcitationImpl.get_pauli_double_ex_yordan()
        for key, val in list_signs.items():
            nkey = key[1::-1] + key[2:]
            nlist[nkey] = val * olist[key] * olist[nkey]
        list_signs = nlist
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.x(q1), circ.x(q3)
        circ.cx(q0, q2)
        circ.h(q1), circ.ry(list_signs["YXXX"]*parameter*0.25, q0)
        circ.cx(q0, q1)
        circ.h(q3), circ.ry(list_signs["XYXX"]*parameter*0.25, q0)
        # cx_via_bridge(circ, q0, q2, q3)
        circ.cx(q0, q3)
        circ.ry(list_signs["XYYY"]*parameter*0.25, q0)
        circ.cx(q0, q1)
        circ.h(q2), circ.ry(list_signs["YXYY"]*parameter*0.25, q0)
        circ.cx(q0, q2)
        circ.ry(list_signs["XXXY"]*parameter*0.25, q0)
        circ.cx(q0, q1)
        circ.ry(list_signs["YYXY"]*parameter*0.25, q0)
        # cx_via_bridge(circ, q0, q2, q3)
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
    def get_pauli_double_ex_12cnot():
        return { "YYIZ" : 1, "XXIZ" : 1, "YYZI" : 1, "XXZI" : 1,"IZYY" : -1,"IZXX" : -1, "ZIYY" : -1,"ZIXX" : -1}
    
    @staticmethod
    def double_ex_12cnot(circ: QuantumCircuit,
                        q0, q1, q2, q3,
                        parameter,
                        list_signs):
        q2,q3 = q3,q2
        olist = ExcitationImpl.get_pauli_double_ex_12cnot()
        nlist = {}
        for key, val in list_signs.items():
            nkey = key[:2] + key[-1:-3:-1]
            
            nlist[nkey] = val * olist[key] * olist[nkey]
            nlist[nkey] = val
        list_signs = nlist
        circ.ry(pi/2, q2)
        circ.rz(pi/2, q3)
        circ.rz(pi/2, q0)
        circ.ry(pi/2, q0)
        circ.cx(q1, q2)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.ry( parameter*list_signs["XXZI"]*0.25, q0)
        circ.ry( parameter*list_signs["YYZI"]*0.25, q1)
        circ.ry( parameter*list_signs["IZYY"]*0.25, q2)
        circ.ry( parameter*list_signs["IZXX"]*0.25, q3)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.cx(q1, q2), circ.cx(q3, q0)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.ry( parameter*list_signs["XXIZ"]*0.25, q0)
        circ.ry( parameter*list_signs["YYIZ"]*0.25, q1)
        circ.ry( parameter*list_signs["ZIYY"]*0.25, q2)
        circ.ry( parameter*list_signs["ZIXX"]*0.25, q3)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.cx(q3, q0)
        circ.ry(-pi/2, q2)
        circ.ry(-pi/2, q0)
        circ.rz(-pi/2, q3)
        circ.rz(-pi/2, q0)

    @staticmethod
    def double_ex_12cnot_ion(circ: QuantumCircuit,
                        q0, q1, q2, q3,
                        parameter,
                        list_signs):
        from qiskit.circuit.library import MSGate
        circ.cz(q1, q2)
        if (list_signs["YYZI"] == -1):
            pass
        if (list_signs["XXZI"] == -1):
            pass
        ms = MSGate(num_qubits=2)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.ry( parameter*list_signs["XXZI"]*0.25, q0)
        circ.ry( parameter*list_signs["YYZI"]*0.25, q1)
        circ.ry( parameter*list_signs["IZYY"]*0.25, q2)
        circ.ry( parameter*list_signs["IZXX"]*0.25, q3)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.cx(q1, q2), circ.cx(q3, q0)
        circ.cx(q0, q1), circ.cx(q2, q3)
        circ.ry( parameter*list_signs["XXIZ"]*0.25, q0)
        circ.ry( parameter*list_signs["YYIZ"]*0.25, q1)
        circ.ry( parameter*list_signs["ZIYY"]*0.25, q2)
        circ.ry( parameter*list_signs["ZIXX"]*0.25, q3)
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
        

