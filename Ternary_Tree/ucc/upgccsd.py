from __future__ import annotations
from functools import reduce
from typing import Tuple, List, Dict, Optional


from .abstractucc import AbstractUCC
from ..qiskit_interface import QiskitCirc 
from ..utils import CircWrapper, Pauli, Parameter
from ..utils import lad2lad, lad2maj, MajoranaContainer, MajExcitation 
from ..utils import SingleLadExcitation, DoubleLadExcitation,  LadExcitation

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
     
    
class UpGCCSD(AbstractUCC):
    def __init__(self,  **kwargs):
        super().__init__(**kwargs)

    def get_alpha_excitations(self) -> List[SingleLadExcitation]:
        """
        Method to get spin preserving single excitation operators via ladder
        operators for alpha orbitals
        """
        alpha = range(self.n_spatial)
        return [SingleLadExcitation((i,), (j,)) for i in alpha for j in alpha if i > j]
 
    def get_beta_excitations(self) -> List[SingleLadExcitation]:
        """
        Method to get spin preserving single excitation operators via ladder
        operators for beta orbitals
        """
        beta = range(self.n_spatial, 2*self.n_spatial)
        return [SingleLadExcitation((i,), (j,)) for i in beta for j in beta if i > j]

    def get_double_excitations(self) -> List[DoubleLadExcitation]:
        """
        Method to get spin preserving pair double excitation operators via ladder
        operators
        """
        alpha = range(self.n_spatial)
        beta_index_offset = self.n_spatial
        ab = [DoubleLadExcitation((a + beta_index_offset, a), (b + beta_index_offset, b))
              for a in alpha for b in alpha if a < b]
        return ab

    def swap2xn(self, 
                number_of_layers: int=1,
                method_name: str=LadExcImpl.CNOT12(),
                circ: Optional[CircWrapper]=None
                ) -> Tuple[CircWrapper, MajoranaContainer]:
                
        # circuit generation
        n = self.n_qubits
        if circ is None:
            circ = QiskitCirc(n)
        mtoq = MajoranaContainer.jw(n)
        single_ladder_exc = self.get_alpha_excitations() + self.get_beta_excitations()
        double_ladder_exc = self.get_double_excitations()
        

        # supporting functions 
        def init_state():
            num_electrons = self.n_alpha  
            for i in range(self.n_qubits):
                if (mtoq.get_by_qubs(i)[0] < 2*num_electrons) or (self.n_spatial * 2 <= 
                        mtoq.get_by_qubs(i)[0]  < 2*self.n_spatial + 2*num_electrons):
                    circ.x(i)
                else:
                    circ.id(i)
            for i in range(2):
                for j in range(i*n//2, i*n//2 + n//2 - 1, 2):
                    mswap(j, 1, j + 1, 0)
                    
        def single_exc_layer():
            func = lambda i: MajExcitation((mtoq.get_by_qubs(i)[0], mtoq.get_by_qubs(i)[1]))
            _ = [append_maj_S(func(i)) for i in range(n)]

        def mswap_layer(parity):
            for i in range(2):
                for j in range(i * n//2, i * n//2 + n//2 - 1, 2):
                    mswap(j, 1, j + 1, 0)
            for i in range(2):
                for j in range(i * n//2 + 1, i * n//2 + n//2 - 2, 2):
                    mswap(j, 1, j + 1, 0)   
            
            for i in range(parity % 2, n//2 - parity % 2):
                mswap(i, 1, i, 0)
            for i in range(parity % 2 + n//2, n - parity % 2):
                mswap(i, 1, i, 0)

        def double_exc_layer(parity):
            for i in range(parity % 2, n//2 - 1 - parity % 2, 2):
                append_lad_D([i, i + 1, i + n//2, i + n//2 + 1])

        def mswap(qub1, disp1, qub2, disp2):
            return _mswap(qub1, disp1, qub2, disp2, circ, mtoq)
        
        def append_maj_S(ls: MajExcitation):
            return _append_maj_S(ls, sp_maj_exc, circ, mtoq)

        def append_lad_D(qubits: List[int]):
            return _append_lad(qubits,  circ, mtoq, method_name, dp_lad_exc)    
        

        init_state()
        for k in range(number_of_layers):
            sp_maj_exc: Dict[MajExcitation, Parameter] = lad2maj(single_ladder_exc, name="t" + str(k)+ "_")
            dp_lad_exc: Dict[LadExcitation, Parameter] = lad2lad(double_ladder_exc, name="t" + str(k)+ "_")
            # if n == 4:
            #     single_exc_layer()
            #     double_exc_layer(0)
            #     for i in range(2):
            #         for j in range(i*n//2, i*n//2 + n//2 - 1, 2):
            #             mswap(j, 1, j + 1, 0)
            # else:
            for i in range(n//2):
                single_exc_layer()
                double_exc_layer(i)
                mswap_layer(i+1)
            print(dp_lad_exc)
            # print(sp_lad_exc)
        return circ, mtoq
    
    def swap_gen(self, 
                number_of_layers: int=1, 
                method_name: str=LadExcImpl.YORDAN(),
                circ: Optional[CircWrapper]=None
                ) -> Tuple[CircWrapper, MajoranaContainer]:
        n = self.n_qubits
        if circ is None:
            circ = QiskitCirc(n)
        mtoq = MajoranaContainer.jw(n)
        list_enum = []
        for j in range(n//4):
            for i in range(2*j, 2*j + 2):
                list_enum.extend([2*i, 2*i + 1])
            for i in range(2*j, 2*j + 2):
                list_enum.extend([2*i + n, 2 * i + n + 1])
        if (n%4 != 0):
            list_enum.extend([n - 2, n -1, 2*n-2, 2*n - 1])

        # mtoq.qubs = list_enum
        mtoq.renumerate(list_enum)
        single_ladder_exc = self.get_alpha_excitations() + self.get_beta_excitations()
        double_ladder_exc = self.get_double_excitations()

        def append_lad_S(i):
            _append_lad([i, i + 1], circ, mtoq, LadExcImpl.SINGLE(), sp_lad_exc )
        
        def append_lad_D(i):
            return _append_lad([i, i + 1, i + 2, i + 3], circ, mtoq, method_name, dp_lad_exc )

        def init_state():
            num_electrons = self.n_alpha  
            for i in range(self.n_qubits):
                if (mtoq.get_by_qubs(i)[0] < 2*num_electrons) or (self.n_spatial * 2 <= 
                        mtoq.get_by_qubs(i)[0]  < 2*self.n_spatial + 2*num_electrons):
                    circ.x(i)
                else:
                   circ.id(i)
                    
        def single_exc_layer(parity):
            for i in range(0, n-1, 2):
                append_lad_S(i)
            for i in range(1, n-1, 2):
                append_lad_S(i)
                

        def fswap(fq, circ=circ):
            # list_trans = [(fq, 0, fq + 1, 0), (fq, 1, fq + 1, 1)]
            _mswap(fq, 0, fq + 1, 0, circ, mtoq)
            _mswap(fq, 1, fq + 1, 1, circ, mtoq)
            # circ.fswap([fq, fq + 1])

        def fswap_layer():
            for i in range(0, n-1, 2):
                fswap(i)
            for i in range(1, n-1, 2):
                fswap(i)

        def double_exc_layer(parity):
            if (parity % 2 == 0):
                for i in range(0, n-3,4):
                    append_lad_D(i)
            else:
                for i in range(2, n-3,4):
                    append_lad_D(i)
                
            for i in range(0, n-3,4):
                append_lad_D(i)

        init_state()
        for k in range(number_of_layers):
            sp_lad_exc: Dict[LadExcitation, Parameter] = lad2lad(single_ladder_exc, name="t" + str(k)+ "_")
            dp_lad_exc: Dict[LadExcitation, Parameter] = lad2lad(double_ladder_exc, name="t" + str(k)+ "_")
            if n == 4:
                single_exc_layer(0)
                double_exc_layer(0)
            else:
                for i in range(n//2):
                    single_exc_layer(i)
                    double_exc_layer(i)
                    fswap_layer()
            print(dp_lad_exc)
            print(sp_lad_exc)
        return circ, mtoq
    
    def get_jw_opt(self) -> MajoranaContainer:
        n = self.n_qubits
        mtoq = MajoranaContainer.jw(n)
            
        def fswap(fq):
            list_trans = [(fq, 0, fq + 1, 0), (fq, 1, fq + 1, 1)]
            _ = [mtoq.transpose(*el) for el in list_trans] 

        def double_exc_tree():
            _ = [fswap(n//2 - l + 2*i - 1) for l in range(n // 2 - 1) for i in range(l + 1)]
        double_exc_tree()
        return mtoq
    
    
def _mswap(qub1, disp1, qub2, disp2, circ: CircWrapper, mtoq: MajoranaContainer):
    circ.mswap(mtoq.transpose(qub1, disp1, qub2, disp2))
    
def _append_maj_S(ls: MajExcitation,
                 sp_maj_exc: Dict[MajExcitation, Parameter],
                 circ: CircWrapper, 
                 mtoq: MajoranaContainer) -> bool:
    try:
        paulis = [mtoq[i] for i in ls]
        circ.rotation(reduce(lambda x,y: x * y, paulis), sp_maj_exc.pop(ls))
        return True
    except KeyError:
        return False
    
def _append_lad(qubits: Tuple[int,...],
                circ: CircWrapper, 
                mtoq: MajoranaContainer, 
                method_name: str,
                lad_exc: Dict[LadExcitation, Parameter]) -> bool:
    maj_set = []
    list_signs = getattr(CircWrapper, "get_pauli_" + method_name)()
    for i in qubits:
        maj_set = maj_set + list(mtoq.get_by_qubs(i))
    maj_set = list(sorted(maj_set))
    try:
        if len(qubits) == 4:
            exc = DoubleLadExcitation((maj_set[4]//2, maj_set[0]//2), (maj_set[6]//2, maj_set[2]//2))
        else:
            exc = SingleLadExcitation((maj_set[2]//2,), (maj_set[0]//2,))
        maj_exc: List[MajExcitation] = lad2maj([exc])
        for maj in maj_exc:
            pauli = reduce(lambda x,y: x * y, [mtoq[q] for q in maj.op])
            label = pauli.get_label_carr(qubits)
            list_signs[label] = (list_signs[label] * maj.sign * 1j**pauli.pow ).imag
            if len(pauli.get_label_qubs()[0]) > len(label):
                raise KeyError
        getattr(circ, method_name)(qubits, lad_exc.pop(exc), list_signs)
    except KeyError:
        return False
