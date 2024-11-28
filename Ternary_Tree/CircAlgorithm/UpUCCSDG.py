from __future__ import annotations

import itertools
from . import TernaryTree
from .AbstractUCC import AbstractUCC
from .Mycirq import MyCirq
import copy
from qiskit import QuantumRegister
from itertools import chain
from functools import reduce
from .paulistring import prod_pauli_strings, PauliStringSign


dict_prods = {(1, 1, 1, 2) : 1j, (1,1,2,1): 1j, (1,2,1,1): -1j, (2,1,1,1): -1j, (1,2,2,2): 1j, (2,1,2,2): 1j, (2,2,1,2): -1j,  (2,2,2,1): -1j}


class UpUCCSDG(AbstractUCC):
    def __init__(self,  **kwargs):
        super().__init__(**kwargs)

    def get_alpha_excitations(self) -> list[tuple[int, int]]:
        """
        Method to get spin preserving single excitation operators via ladder
        operators for alpha orbitals
        """
        alpha = range(self.num_spatial_orbitals)
        return [(i, j) for i in alpha for j in alpha if i > j]
 
    def get_beta_excitations(self) -> list[tuple[int, int]]:
        """
        Method to get spin preserving single excitation operators via ladder
        operators for beta orbitals
        """
        beta = range(self.num_spatial_orbitals, 2*self.num_spatial_orbitals)
        return [(i, j) for i in beta for j in beta if i > j]


    def get_double_excitations(self):
        """
        Method to get spin preserving pair double excitation operators via ladder
        operators
        """
        alpha = range(self.num_spatial_orbitals)
        beta_index_offset = self.num_spatial_orbitals
        ab = [(a + beta_index_offset, a, b + beta_index_offset, b)
              for a in alpha for b in alpha if a < b]
        return ab

    def get_parametrized_circuit(self, number_of_layers=1, type=0):
        single_ladder_alpha_exc = self.get_alpha_excitations()
        single_ladder_beta_exc = self.get_beta_excitations()
        single_ladder_exc = single_ladder_alpha_exc + single_ladder_beta_exc

        double_ladder_exc = self.get_double_excitations()
        n = self.n_qubits
        qr = QuantumRegister(n)
        cirq = MyCirq(qr)
            
        self.tt = self.get_jw_opt(type)
        
        def append_maj_single_exc(cirq, tt, ls: tuple) -> bool:
            try:
                psss = [PauliStringSign(*tt.get_majorana(i)) for i in ls]
                cirq.rotation(reduce(prod_pauli_strings, psss), sp_maj_exc.pop(ls))
                return True
            except KeyError:
                return False
        
        def append_lad_double_exc(cirq: MyCirq, tt, i) -> bool:
            qubits = [i, i + 1, i + 2, i + 3]
            maj_nums = [self.tt[qubit][0].num for qubit in qubits] + [self.tt[qubit][1].num for qubit in qubits]
            init_orb = min(maj_nums) // 2
            fin_orb = (max(maj_nums) - n - 1) // 2 
            list_signs = {"XXXY": 1, "XXYX": 1, "XYXX": 1, "YXXX": 1, "YYYX": 1, "YYXY": 1, "YXYY": 1, "XYYY": 1}
            if type in {2}:
                list_signs = {"ZZZY": -1, "ZZYZ": -1, "ZYZZ": 1, "YZZZ": 1, "YYYZ": -1, "YYZY": -1, "YZYY": 1, "ZYYY": 1}
            if type in {3}:
                list_signs = {"ZZZY": -1, "ZZYZ": 1, "ZYZZ": 1, "YZZZ": -1, "YYYZ": 1, "YYZY": -1, "YZYY": -1, "ZYYY": 1}
            try:
                for key in dict_prods:
                    psss = [PauliStringSign(*tt.get_majorana(j1 + j2)) for j1,j2 in zip(key, (2*init_orb, n + 2*init_orb, 2*fin_orb, n + 2*fin_orb))]
                    pss = reduce(prod_pauli_strings, psss)
                    if set(qubits) == set(pss.qubits):
                        pauli = sorted(pss.ps.items(), key=lambda x: x[0])
                        pauli = "".join([el[1] for el in pauli])
                        list_signs[pauli] = (list_signs[pauli]*dict_prods[key]*pss.sign*(-1j)).real
                    else:
                        return False
                if type==1:
                    qubits = [i, i + 1, i + 2, i + 3]
                    cirq.double_ex_yordan(qubits, dp_lad_exc.pop((init_orb + n//2, init_orb, fin_orb + n//2, fin_orb)), list_signs)
                if type==0:
                    cirq.double_ex_opt(qubits, dp_lad_exc.pop((init_orb + n//2, init_orb, fin_orb + n//2, fin_orb)), list_signs)
                if type==2:
                    cirq.double_ex_zyx_opt(qubits, dp_lad_exc.pop((init_orb + n//2, init_orb, fin_orb + n//2, fin_orb)), list_signs)
                if type==3:
                    print(list_signs)
                    cirq.double_ex_zyx_yordan(qubits, dp_lad_exc.pop((init_orb + n//2, init_orb, fin_orb + n//2, fin_orb)), list_signs)
            except KeyError:
                print("HERERERE")
                return False


        def maj_swap(node1, edge1, node2, edge2, cirq=cirq, unsigned=False):
            cirq.maj_swap(PauliStringSign(self.tt.branch_transposition(node1, edge1, node2, edge2, unsigned)))


        def fswap(fq, cirq=cirq):
            list_trans = [(fq, 0, fq + 1, 0), (fq, 1, fq + 1, 1)]
            _ = [self.tt.branch_transposition(*el, unsigned=True) for el in list_trans]  
            cirq.fswap([fq, fq + 1])
            # _ = [cirq.maj_swap(PauliStringSign(self.tt.branch_transposition(*el))) for el in list_trans]  

        
        def f2swap(fq: int, cirq: MyCirq = cirq):
            list_trans = [(fq, 0, fq + 2, 0), (fq, 1, fq + 2, 1), (fq+1,0,fq+3,0),(fq+1,1,fq+3,1)]
            _ = [self.tt.branch_transposition(*el, unsigned=True) for el in list_trans]  
            cirq.fermionic_2swap([fq, fq + 1, fq + 2, fq + 3])

        
        def init_state():
            # _ = [cirq.id(i) for i in range(n)]
            # _ = [cirq.x(i) if ((i < self.num_alpha) or (0 <= (i - n//2)< self.num_beta)) 
            #                 else cirq.id(i) for i in range(n)]

            num_electrons = self.num_alpha  
            for i in range(self.n_qubits):
                if (self.tt[i][0].num <= 2*num_electrons) or (self.num_spatial_orbitals * 2 + 1<= 
                        self.tt[i][0].num <= 2*self.num_spatial_orbitals + 2*num_electrons):
                    cirq.x(i)
                    # TODO
                else:
                    cirq.id(i)
            if type in {2,3}:
                for i in range(self.n_qubits):
                    cirq.h(i)
                    cirq.z(i)
        def single_prep(cirq):
            res, disp = 0, 0
            for i in range(0, n//2):
                maj_swap(2*i, 0, 2*i + 1, 1)
            for k in range(1, n//4):
                for i in range(0, n//4 - k):
                    fswap(2*i + k, cirq)
                    fswap(2*i + k + n//2, cirq)


        def single_exc():
            def _maj():
                func = lambda i: tuple(sorted([self.tt[i][0].num, self.tt[i][1].num]))
                _ = [append_maj_single_exc(cirq, self.tt, func(i)) for i in range(n)]

            def ferm_trans(i, j):
                _ = [(_maj(), maj_swap(i,n,j,m)) for n,m in ((0,1),(1,1),(1,0))]
            
            def last_ferm_trans(i, j):
                _ = [(_maj(), maj_swap(i,n,j,m)) for n,m in ((0,1),(1,1))]
                _ = _maj()
            # cirq.barrier()
            _maj()
            # change
            for block in range(4):
                _ = [ferm_trans(j + block*n//4, j + block*n//4 + 1) for i in range(n//4)
                        for j in range(i % 2, n//4 - 1, 2)]
                # _ = [last_ferm_trans(j + block*n//4, j + block*n//4 + 1) for j in range((n//4-1) % 2, n//4 - 1, 2)]
        
        def single_anti_prep():
            res, disp = 1, n//4 - 1
            for k in range(n//4 - 1, 0, -1):
                for i in range(0, n//4 - k):
                    fswap(2*i + k, cirq)
                    fswap(2*i + k + n//2, cirq)
    
            for i in range(0, n//2):
                maj_swap(2*i, 0, 2*i + 1, 1)

        def double_exc_tree():
            # def tr(i, l):
                # maj_swap(n//2 - l + 2*i - 1, int(l>=2*i), n//2 - l + 2*i, int(l>2*i))

            _ = [fswap(n//2 - l + 2*i - 1) for l in range(n // 2 - 1) for i in range(l + 1)]
            # _ = [tr(i, l) for l in range(n // 2 - 2, -1 , -1) for i in range(l + 1)]
            # cirq.barrier()
            # _ = [maj_swap for l in range(n // 2) for i in range(l + 1)]
        def double_anti_exc_tree():
            _ = [fswap(n//2 - l + 2*i - 1) for l in range(n // 2 - 2, -1 , -1) for i in range(l + 1)]

        def double_exc_layer():
            # f2swap(fq)
            # cirq.double_ex_opt()
            def layer(init: int):
                _ = [append_lad_double_exc(cirq, self.tt, i) for i in range(init, n-3, 4)]
            def swap_layer(init: int):
                _ = [f2swap(i, cirq) for i in range(init, n-3, 4)]
            
            layer(0)
            layer(2)
            for _ in range(n//4 - 1):
                swap_layer(0)
                swap_layer(2)
                layer(0)
                layer(2)
        

        k = 0
        init_state()
        # sp_maj_exc = copy.deepcopy(self.to_par_maj_excitations(single_ladder_exc, name="t" + str(k) + "_"))
        dp_lad_exc = copy.deepcopy(self.to_par_ladder_exciations(double_ladder_exc, name="t" + str(k) + "_"))
        sp_maj_exc = self.to_par_maj_exitations_comp(
            copy.deepcopy(self.to_par_maj_excitations(single_ladder_alpha_exc, name="t" + str(k) + "_")),
            n
        )
        double_exc_layer()
        double_anti_exc_tree()
        single_prep(cirq)
        # cirq.barrier()
        single_exc()
        for k in range(1, number_of_layers):
            dp_lad_exc = copy.deepcopy(self.to_par_ladder_exciations(double_ladder_exc, name="t" + str(k) + "_"))
            sp_maj_exc = self.to_par_maj_exitations_comp(
                copy.deepcopy(self.to_par_maj_excitations(single_ladder_alpha_exc, name="t" + str(k) + "_")),
                n
            )
            single_anti_prep()
            double_exc_tree()
            double_exc_layer()
            double_anti_exc_tree()
            single_prep(cirq)
            single_exc()

        # print(sp_maj_exc, dp_lad_exc)
        return cirq

    def get_jw_opt(self, type=0):
        n = self.n_qubits
        tt = TernaryTree(n)
        if type in {2,3}:
            tt.gate_name_to_number = {'Z': 0, "Y": 1, 'X': 2}
            tt.gate_number_to_name = {0: 'Z', 1: 'Y', 2: 'X'}
            
        def fswap(fq):
            list_trans = [(fq, 0, fq + 1, 0), (fq, 1, fq + 1, 1)]
            _ = [tt.branch_transposition(*el, unsigned=True) for el in list_trans] 

        def double_exc_tree():
            _ = [fswap(n//2 - l + 2*i - 1) for l in range(n // 2 - 1) for i in range(l + 1)]
        double_exc_tree()
        return tt