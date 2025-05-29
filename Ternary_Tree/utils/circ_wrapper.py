
 
from typing import Union
from abc import abstractmethod
from numpy import pi



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
     

class Parameter:
    """
        var = coef * name
    """
    def __init__(self, name: str='t'):
        self.name = name
        self.coef = 1

    def __repr__(self):
        return self.name
        
    def __str__(self):
        return self.name
    
MyParameter = Parameter


from .pauli import Pauli

class CircWrapper:
    @abstractmethod
    def mswap(self, pauli: Pauli):
        pass
    
    @abstractmethod
    def fswap(self, pauli: Pauli):
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
    
    @staticmethod
    def get_pauli_single_ex():
        return {"XY": -1, "YX": 1}

    @staticmethod
    def get_pauli_double_ex_short():
        return  {"XXXY": 1, "XXYX": 1, "XYXX": -1, "YXXX": -1, "YYYX": -1, "YYXY": -1, "YXYY": 1, "XYYY": 1}
    
    @staticmethod
    def get_pauli_double_ex_yordan():
        return {"XXXY": -1, "XXYX": 1, "XYXX": -1, "YXXX": 1, "YYYX": 1, "YYXY": -1, "YXYY": 1, "XYYY": -1}
            
    @staticmethod
    def get_pauli_double_ex_12cnot():
        return {"YYIZ": 1,"XXIZ": 1,"YYZI": 1,"XXZI": 1,"IZYY": -1,"IZXX": -1,"ZIYY": -1,"ZIXX": -1}

    @abstractmethod
    def single_ex(self, qubits, parameter: MyParameter, list_signs):
        pass
        
    @abstractmethod
    def double_ex_short(self, qubits, par, list_signs):
        pass
    
    @abstractmethod
    def double_ex_zyx_opt(self, qubits, par, list_signs):
        pass

    @abstractmethod
    def double_ex_yordan(self, qubits, par, list_signs):
        pass

    @abstractmethod
    def double_ex_zyx_yordan(self, qubits, par, list_signs):
        pass

    @abstractmethod
    def double_ex_12cnot(self, qubits, par, list_signs):
        pass

    @abstractmethod
    def fermionic_2swap(self, qubits):
        pass
  


