from __future__ import annotations
from functools import reduce
from typing import Tuple, List, Dict, Optional
import logging

from numpy import pi
logger = logging.getLogger(__name__)
# from .logger import logger, x


from .abstractucc import AbstractUCC, Molecule
from ..qiskit_interface import QiskitCirc
from ..utils.circ_wrapper import ExcitationImpl
from ..utils import CircWrapper, Pauli, Parameter, alpha2beta
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
    def CNOT12xyz():
        return "double_ex_12cnot" 
    
    @staticmethod
    def SINGLE():
        return "single_ex"
    
    @staticmethod
    def SHORT():
        return "double_ex_short"
     
    
class UpGCCSD(AbstractUCC):

    def __init__(self, molecule: Molecule=Molecule()):
        super().__init__(molecule)

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
        bio = self.n_spatial
        ab = [DoubleLadExcitation((a + bio, a), (b + bio, b))
              for a in alpha for b in alpha if a < b]
        return ab

    def swapgenmaj(self, 
            number_of_layers: int=1, 
            method_name: str=LadExcImpl.YORDAN(),
            circ: Optional[CircWrapper]=None
            ) -> Tuple[CircWrapper, MajoranaContainer]:
 
        logger.info("Generalized Majorana swap circuit generation...")
        n = self.n_qubits
        if circ is None:
            circ = QiskitCirc(n)
        mtoq = MajoranaContainer.jw(n)
        if ((n == 8 or n == 10) and self.n_alpha != 1):
            list_enum = [i for i in range(n*2)]
            for i in range(0, 2*n, n):
                for j in range(2):
                    list_enum[2 + i + j], list_enum[4 + i + j] = list_enum[4 + i + j], list_enum[2 + i + j]
            mtoq.renumerate(list_enum)

        logger.info(f"{mtoq=}")
        list_enum = []
        for j in range(n//4):
            for i in range(2*j, 2*j + 2):
                list_enum.extend([2*i, 2*i + 1])
            for i in range(2*j, 2*j + 2):
                list_enum.extend([2*i + n, 2 * i + n + 1])
        if (n%4 != 0):
            list_enum.extend([n - 2, n -1, 2*n-2, 2*n - 1])
        mtoq.renumerate(list_enum)

        single_ladder_exc = self.get_alpha_excitations() + self.get_beta_excitations()
        double_ladder_exc = self.get_double_excitations()
        
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
                   
        def append_maj_S(ls: MajExcitation):
            return _append_maj_S(ls, sp_maj_exc, circ, mtoq)
               
        def single_exc_layer():
                    func = lambda i: MajExcitation((mtoq.get_by_qubs(i)[0], mtoq.get_by_qubs(i)[1]))
                    _ = [append_maj_S(func(i)) for i in range(n)]

        def mswap(qub1, disp1, qub2, disp2):
            return _mswap(qub1, disp1, qub2, disp2, circ, mtoq)
        
        def mswap_single_exc_layer(parity):
            parity = parity % 2
            for i in range(parity, n - 1, 2):
                mswap(i,0,i + 1, 1)
            single_exc_layer()
            for i in range(parity, n - 1, 2):
                mswap(i,1,i + 1, 0)
                mswap(i,0,i,1)
                mswap(i + 1,0,i + 1,1)
            

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
            sp_maj_exc = lad2maj(self.get_alpha_excitations() + self.get_beta_excitations(), name="t" + str(k)+ "_")
            dp_lad_exc: Dict[LadExcitation, Parameter] = lad2lad(double_ladder_exc, name="t" + str(k)+ "_")
            if n == 0:
                double_exc_layer(0)
            else:
                for i in range(n//2):
                    double_exc_layer(i)
                    mswap_single_exc_layer(0)
                    mswap_single_exc_layer(1)

        logger.info(f"{dp_lad_exc=}")
        logger.info(f"{sp_maj_exc=}")
        return circ, mtoq
    
    
    def swap2xnferm(self, 
                number_of_layers: int=1, 
                method_name: str=LadExcImpl.YORDAN(),
                circ: Optional[CircWrapper]=None
                ) -> Tuple[CircWrapper, MajoranaContainer]:
 
        logger.info("Fermion swap2xn circuit generation...")
        n = self.n_qubits
        circ = circ
        if circ is None:
            circ = QiskitCirc(n)
        mtoq = MajoranaContainer.jw(n)
        if ((n == 8 or n == 10) and self.n_alpha != 1):
            list_enum = [i for i in range(n*2)]
            for i in range(0, 2*n, n):
                for j in range(2):
                    list_enum[2 + i + j], list_enum[4 + i + j] = list_enum[4 + i + j], list_enum[2 + i + j]
                    # list_enum[3 + i + j], list_enum[5 + i + j] = list_enum[5 + i + j], list_enum[3 + i + j]
            mtoq.renumerate(list_enum)
        logger.info(f"{mtoq=}")
        double_ladder_exc = self.get_double_excitations()

        def init_state():
            num_electrons = self.n_alpha  
            for i in range(self.n_qubits):
                if (mtoq.get_by_qubs(i)[0] < 2*num_electrons) or (self.n_spatial * 2 <= 
                        mtoq.get_by_qubs(i)[0]  < 2*self.n_spatial + 2*num_electrons):
                    circ.x(i)
                else:
                   circ.id(i)
                   
        def single_exc_layer():
            func = lambda i: MajExcitation((mtoq.get_by_qubs(i)[0], mtoq.get_by_qubs(i)[1]))
            _ = [append_maj_S(func(i)) for i in range(n)]

        def mswap_single_exc_layer(parity):
            parity = parity % 2
            for spin in range(2):
                for i in range(n//2*spin + parity, n//2 + n//2*spin - 1, 2):
                    mswap(i,0,i + 1, 1)
            single_exc_layer()
            for spin in range(2):
                for i in range(n//2*spin + parity, n//2 + n//2*spin - 1, 2):
                    mswap(i,1,i + 1, 0)
                    mswap(i,0,i,1)
                    mswap(i + 1,0,i + 1,1)

        def double_exc_layer(parity):
            parity = parity % 2
                # for i in range(0, n-3,4):
            for i in range(parity, n//2 - 1, 2):
                qubits = [i, i + 1, i + n // 2, i + 1 + n // 2 ]
                append_lad_D(qubits)

        def mswap(qub1, disp1, qub2, disp2):
            return _mswap(qub1, disp1, qub2, disp2, circ, mtoq)
        
        def append_maj_S(ls: MajExcitation):
            return _append_maj_S(ls, sp_maj_exc, circ, mtoq)

        def append_lad_D(qubits: List[int]):
            return _append_lad(qubits,  circ, mtoq, method_name, dp_lad_exc)   
        
        init_state()
        for k in range(number_of_layers):
            sp_maj_exc = lad2maj(self.get_alpha_excitations() + self.get_beta_excitations(), name="t" + str(k)+ "_")
            dp_lad_exc: Dict[LadExcitation, Parameter] = lad2lad(double_ladder_exc, name="t" + str(k)+ "_")
            if n == 0:
                double_exc_layer(0)
            else:
                for i in range(n//2):
                    double_exc_layer(i)
                    mswap_single_exc_layer(i)
        logger.info(f"{dp_lad_exc=}")
        logger.info(f"{sp_maj_exc=}")
        return circ, mtoq
    
    def swap2xn12(self, 
                number_of_layers: int=1,
                method_name: str=LadExcImpl.CNOT12xyz(),
                encoding: str="xyz", 
                circ: Optional[QiskitCirc]=None
                ) -> Tuple[CircWrapper, MajoranaContainer]:
        logger.info("Fermion swap2xn circuit generation...")
        n = self.n_qubits
        circ = circ
        if circ is None:
            circ = QiskitCirc(n)
        mtoq = MajoranaContainer.jw(n)
        list_enum = [i for i in range(n*2)]
        if ((n == 8 or n == 10) and self.n_alpha != 1):
            list_enum = [i for i in range(n*2)]
            for i in range(0, 2*n, n):
                for j in range(2):
                    list_enum[2 + i + j], list_enum[4 + i + j] = list_enum[4 + i + j], list_enum[2 + i + j]
                    # list_enum[3 + i + j], list_enum[5 + i + j] = list_enum[5 + i + j], list_enum[3 + i + j]
            # mtoq.renumerate(list_enum)
        for i in range(n//2, 3*n//4):
            for j in range(2):
                list_enum[2*i + j], list_enum[2*(n - i + n//2 - 1) + j] = list_enum[2*(n - i + n//2 - 1) + j], list_enum[2*i + j], 
        mtoq.renumerate(list_enum)
        logger.info(f"{list_enum=}")
        logger.info(f"{mtoq=}")

        double_ladder_exc = self.get_double_excitations()

        def init_state():
            num_electrons = self.n_alpha  
            for i in range(self.n_qubits):
                if (mtoq.get_by_qubs(i)[0] < 2*num_electrons) or (self.n_spatial * 2 <= 
                        mtoq.get_by_qubs(i)[0]  < 2*self.n_spatial + 2*num_electrons):
                    circ.x(i)
                else:
                   circ.id(i)
                   
        def single_exc_layer():
            func = lambda i: MajExcitation((mtoq.get_by_qubs(i)[0], mtoq.get_by_qubs(i)[1]))
            _ = [append_maj_S(func(i)) for i in range(n)]

        def mswap_single_exc_layer(parity):
            parity = parity % 2
            for spin in range(2):
                for i in range(n//2*spin + parity, n//2 + n//2*spin - 1, 2):
                    mswap(i,0,i + 1, 1)
            single_exc_layer()
            for spin in range(2):
                for i in range(n//2*spin + parity, n//2 + n//2*spin - 1, 2):
                    mswap(i,1,i + 1, 0)
                    mswap(i,0,i,1)
                    mswap(i + 1,0,i + 1,1)

        def double_exc_layer(parity):
            parity = parity % 2
            for i in range(parity, n//2 - 1, 2):
                qubits = [i, i + 1, n - 1 - (i), n - 2 - i  ]
                # raise KeyError()
                append_lad_D(qubits)

        def mswap_layer(parity, layer):
            parity = parity % 2
            if layer == 0:
                for i in range(parity, n//2  - 1, 2):
                    mswap(i,0,i + 1, 1)
                for i in range(n - 1 - parity, n//2 , -2):
                    mswap(i, 0, i - 1, 1)
            else:
                for i in range(parity, n//2 - 1, 2):
                    mswap(i,1,i + 1, 0)
                    mswap(i,0,i,1)
                    mswap(i + 1,0,i + 1,1)
                    
                for i in range(n - 1 - parity, n//2 , -2):
                    mswap(i,1,i - 1, 0)
                    mswap(i,0,i,1)
                    mswap(i - 1,1,i - 1,0)

        def mswap(qub1, disp1, qub2, disp2):
            return _mswap(qub1, disp1, qub2, disp2, circ, mtoq)
        
        def append_maj_S(ls: MajExcitation):
            return _append_maj_S(ls, sp_maj_exc, circ, mtoq)

        def append_lad_D(qubits: List[int]):
            return _append_lad(qubits,  circ, mtoq, method_name, dp_lad_exc)   
        
        init_state()
        for k in range(number_of_layers):
            sp_maj_exc = lad2maj(self.get_alpha_excitations() + self.get_beta_excitations(), name="t" + str(k)+ "_")
            dp_lad_exc: Dict[LadExcitation, Parameter] = lad2lad(double_ladder_exc, name="t" + str(k)+ "_")
            # if n == 0:
            #     double_exc_layer(0)
            # else:
            for i in range(n//2):
                mswap_layer(i, 0)
                double_exc_layer(i)
                single_exc_layer()
                mswap_layer(i, 1)
        logger.info(f"{dp_lad_exc=}")
        logger.info(f"{sp_maj_exc=}")
        return circ, mtoq

    def swap2xn(self, 
                number_of_layers: int=1,
                method_name: str=LadExcImpl.CNOT12xyz(),
                encoding: str="xyz", 
                circ: Optional[QiskitCirc]=None
                ) -> Tuple[CircWrapper, MajoranaContainer]:
        logger.info("Swap2xn circuit generation...")
        encoding = encoding.lower()
        # method_name = getattr(LadExcImpl, method_name)()
        n = self.n_qubits
        circ = circ
        if circ is None:
            circ = QiskitCirc(n)
        mtoq = MajoranaContainer.jw(n, encoding=encoding)
        single_ladder_exc = self.get_alpha_excitations() + self.get_beta_excitations()
        double_ladder_exc = self.get_double_excitations()

        def init_state():
            num_electrons = self.n_alpha  
            nels = [i for i in range(num_electrons)] + [i + n//2 for i in range(num_electrons)]
            ExcitationImpl.majorana_init_state(nels, circ, mtoq, encoding)
                    
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
            sp_maj_exc = lad2maj(self.get_alpha_excitations() + self.get_beta_excitations(), name="t" + str(k)+ "_")
            dp_lad_exc: Dict[LadExcitation, Parameter] = lad2lad(double_ladder_exc, name="t" + str(k)+ "_")
            if n == 0:
                single_exc_layer()
                double_exc_layer(0)
            else:
                for i in range(n//2):
                    double_exc_layer(i)
                    single_exc_layer()
                    mswap_layer(i+1)
        logger.info(f"{dp_lad_exc=}")
        logger.info(f"{sp_maj_exc=}")
        
        return circ, mtoq
    
    def swap_gen(self, 
                number_of_layers: int=1, 
                method_name: str=LadExcImpl.YORDAN(),
                circ: Optional[CircWrapper]=None
                ) -> Tuple[CircWrapper, MajoranaContainer]:
        logger.info("Fermionic swap circuit generation...")
        n = self.n_qubits
        if circ is None:
            circ = QiskitCirc(n)
        mtoq = MajoranaContainer.jw(n)
        if ((n == 8 or n == 10) and self.n_alpha != 1):
            list_enum = [i for i in range(n*2)]
            for i in range(0, 2*n, n):
                for j in range(2):
                    list_enum[2 + i + j], list_enum[4 + i + j] = list_enum[4 + i + j], list_enum[2 + i + j]
            mtoq.renumerate(list_enum)
        list_enum = []
        for j in range(n//4):
            for i in range(2*j, 2*j + 2):
                list_enum.extend([2*i, 2*i + 1])
            for i in range(2*j, 2*j + 2):
                list_enum.extend([2*i + n, 2 * i + n + 1])
        if (n%4 != 0):
            list_enum.extend([n - 2, n -1, 2*n-2, 2*n - 1])
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
                    
        def single_exc_layer():
            for i in range(0, n-1, 2):
                append_lad_S(i)
            for i in range(1, n-1, 2):
                append_lad_S(i)

        def fswap(fq, circ=circ):
            _mswap(fq, 0, fq + 1, 0, circ, mtoq)
            _mswap(fq, 1, fq + 1, 1, circ, mtoq)

        def fswap_layer():
            for i in range(1, n-1, 2):
                fswap(i)
            for i in range(0, n-1, 2):
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
            sp_lad_exc = lad2lad(self.get_alpha_excitations() + self.get_beta_excitations(), name="t" + str(k)+ "_")
            dp_lad_exc: Dict[LadExcitation, Parameter] = lad2lad(double_ladder_exc, name="t" + str(k)+ "_")
            if n == 4:
                # single_exc_layer()
                double_exc_layer(0)
                single_exc_layer()
            else:
                for i in range(n//2):
                    double_exc_layer(i)
                    single_exc_layer()
                    fswap_layer()
        # logger.info(f"generalized:\n{circ}")
        logger.info(f"{dp_lad_exc=}")
        logger.info(f"{sp_lad_exc=}")
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
    
    
class UCCSD(AbstractUCC):
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
        bio = self.n_spatial
        aa = []
        bb = []
        for a in alpha:
            for b in alpha[a + 1:]:
                for c in alpha[b + 1:]:
                    for d in alpha[c + 1:]:
                        aa.append(DoubleLadExcitation((a, b),(c, d) ))
                        # aa.append(DoubleLadExcitation((a, c),(b, d) ))
                        # aa.append(DoubleLadExcitation((a, d),(b, c) ))
                        bb.append(DoubleLadExcitation((a + bio, b + bio),(c + bio, d + bio) ))
                        # bb.append(DoubleLadExcitation((a + bio, c + bio),(b + bio, d + bio) ))
                        # bb.append(DoubleLadExcitation((a + bio, d + bio),(b + bio, c + bio) ))
        ab = []
        for a in alpha:
            for b in alpha[a + 1:]:
                for c in alpha:
                    for d in alpha[c + 1:]:
                        ab.append(DoubleLadExcitation((a, c + bio), (b, d + bio) ))
        return aa + ab + bb


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
    for i in qubits:
        maj_set = maj_set + list(mtoq.get_by_qubs(i))
    maj_set = list(sorted(maj_set))
    if len(qubits) == 4:
        excs = []
        excs.append(DoubleLadExcitation((maj_set[4]//2, maj_set[0]//2), (maj_set[6]//2, maj_set[2]//2)))
        # excs.append(DoubleLadExcitation((maj_set[2]//2, maj_set[0]//2), (maj_set[6]//2, maj_set[4]//2)))
        # excs.append(DoubleLadExcitation((maj_set[6]//2, maj_set[0]//2), (maj_set[2]//2, maj_set[4]//2)))
        method = circ.double_excitation
    else:
        excs= [SingleLadExcitation((maj_set[2]//2,), (maj_set[0]//2,))]
        method = circ.single_excitation
    for exc in excs:
        list_signs = getattr(ExcitationImpl, "get_pauli_" + method_name)()
        try:
            maj_exc: List[MajExcitation] = lad2maj([exc])
            for maj in maj_exc:
                pauli = reduce(lambda x,y: x * y, [mtoq[q] for q in maj.op])
                label = pauli.get_label_carr(qubits)
                list_signs[label] = (list_signs[label] * maj.sign * 1j**pauli.pow).imag
                if len(pauli.get_label_qubs()[0]) > len(label):
                    raise KeyError
            method(method_name, qubits, lad_exc.pop(exc), list_signs)
        except KeyError:
            pass
