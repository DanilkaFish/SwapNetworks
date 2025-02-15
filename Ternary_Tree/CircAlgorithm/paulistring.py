from __future__ import annotations


dict_prod_coef = {"II" : ["I",1] , "XX" :  ["I",1] , "YY" :  ["I",1] ,"ZZ" :  ["I",1] ,
             "XY" :  ["Z", 1j] ,"YX" : ["Z", -1j],"XZ" : ["Y",-1j], "ZX" : ["Y",1j],"YZ" : ["X",1j],"ZY" : ["X",-1j],
            "IX" : ["X",1], "XI" : ["X",1], "YI" : ["Y",1],"IY" : ["Y",1],"IZ" : ["Z",1],"ZI" : ["Z",1]}

class PauliStringSign:
    def __init__(self, ps=None, sign=1):
        self.ps = ps
        self.sign = sign

    @property
    def ps(self):
        return self._ps

    @ps.setter
    def ps(self, ps):
        if ps is None:
            self._ps = {}
        elif isinstance(ps, dict):
            self._ps = {gate[0]: gate[1] for gate in ps.items() if gate[1] != "I"}
        else:
            self._ps = {gate[0]: gate[1] for gate in ps if gate[1] != "I"}
        
    def __getitem__(self, key: int):
        return self.ps.get(key, "I")
    
    def __setitem__(self, key: int, value: (str, int) | str):
        if value[0] != "I":
            if isinstance(value, str):
                self.ps[key] = value
            else:
                self.ps[key], sign_value = value
                self.sign *= sign_value

    @property
    def qubits(self):
        return list(self.ps.keys())

    @property    
    def pauli_str(self):
        return "".join(list(self.ps.values()))
    
    
def prod_pauli_strings(pss1: PauliStringSign,
                       pss2: PauliStringSign) -> PauliStringSign:

    pss = PauliStringSign(sign=pss1.sign*pss2.sign)
    for qubit in set(pss1.qubits + pss2.qubits):
        pss[qubit] = dict_prod_coef[pss1[qubit] + pss2[qubit]]
    return pss