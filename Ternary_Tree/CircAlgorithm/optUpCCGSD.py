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
dict_prods_single = {(1,1) : 1, (2, 2) : 1}

class UpUCCSDG_opt(AbstractUCC):
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

    def get_parametrized_circuit(self, number_of_layers=1):
        # sp_maj_exc = self.to_par_maj_exitations_comp(
        #     copy.deepcopy(self.to_par_maj_excitations(single_ladder_alpha_exc, name="t" + str(k) + "_")),
        #     n
        # )
        
        single_ladder_alpha_exc = self.get_alpha_excitations()
        single_ladder_beta_exc = self.get_beta_excitations()
        single_ladder_exc = single_ladder_alpha_exc + single_ladder_beta_exc

        double_ladder_exc = self.get_double_excitations()
        n = self.n_qubits
        qr = QuantumRegister(n)
        cirq = MyCirq(qr)
            
        self.tt = TernaryTree(n)

        k = 0
        dp_lad_exc = copy.deepcopy(self.to_par_ladder_exciations(double_ladder_exc, name="t" + str(k) + "_"))
        sp_maj_exc = self.to_par_maj_exitations_comp(
            copy.deepcopy(self.to_par_maj_excitations(single_ladder_alpha_exc, name="t" + str(k) + "_")),
            n
        )

        def maj_swap(node1, edge1, node2, edge2, cirq=cirq, unsigned=False):
            cirq.maj_swap(PauliStringSign(self.tt.branch_transposition(node1, edge1, node2, edge2, unsigned)))

        def append_maj_single_exc(cirq: MyCirq, tt: TernaryTree, ls: tuple) -> bool:
            try:
                psss = [PauliStringSign(*tt.get_majorana(i)) for i in ls]
                cirq.rotation(reduce(prod_pauli_strings, psss), sp_maj_exc.pop(ls))
                return True
            except KeyError:
                return False
            
        def append_lad_double_exc(cirq: MyCirq, tt, a1, a2, b1, b2) -> bool:
            qubits = [a1, a2, b1, b2]
            
            maj_nums = [self.tt[qubit][0].num for qubit in qubits] + [self.tt[qubit][1].num for qubit in qubits]
            init_orb = min(maj_nums) // 2
            fin_orb = (max(maj_nums) - n - 1) // 2 
            list_signs = {"YYIZ": 1,"XXIZ": 1,"YYZI": 1,"XXZI": 1,"IZYY": -1,"IZXX": -1,"ZIYY": -1,"ZIXX": -1}
            try:
                for key in dict_prods:
                    psss = [PauliStringSign(*tt.get_majorana(j1 + j2)) for j1,j2 in zip(key, (2*init_orb, n + 2*init_orb, 2*fin_orb, n + 2*fin_orb))]
                    pss = reduce(prod_pauli_strings, psss)
                    _pauli = sorted(pss.ps.items(), key=lambda x: x[0])
                    pauli = ""
                    i = 0 
                    if len(_pauli) == 3:
                        for el in _pauli:
                            if (el[0] != qubits[i]):
                                pauli += "I"
                                i += 1
                            pauli += el[1]
                            i += 1
                        if len(pauli) == 3:
                            pauli += "I"
                    list_signs[pauli] = (list_signs[pauli]*dict_prods[key]*pss.sign*(-1j)).real
                cirq.double_ex_12cnot(qubits, dp_lad_exc.pop((init_orb + n//2, init_orb, fin_orb + n//2, fin_orb)), list_signs)
            except KeyError:
                print("HERERERE")
                return False

        def init_state():
            num_electrons = self.num_alpha  
            for i in range(self.n_qubits):
                if (self.tt[i][0].num <= 2*num_electrons) or (self.num_spatial_orbitals * 2 + 1<= 
                        self.tt[i][0].num <= 2*self.num_spatial_orbitals + 2*num_electrons):
                    cirq.x(i)
                    # TODO
                else:
                    cirq.id(i)
            for i in range(2):
                for j in range(i*n//2, i*n//2 + n//2 - 1, 2):
                    maj_swap(j, 1, j + 1, 0)
                    # maj_swap(j+1, 0, j+1, 1)
            # cirq.barrier()
                    
        def single_exc_layer():
            func = lambda i: tuple(sorted([self.tt[i][0].num, self.tt[i][1].num]))
            _ = [append_maj_single_exc(cirq, self.tt, func(i)) for i in range(n)]

        def mswap_layer(parity):
            for i in range(2):
                for j in range(i*n//2, i*n//2 + n//2 - 1, 2):
                    maj_swap(j, 1, j + 1, 0)
            for i in range(2):
                for j in range(i*n//2 + 1, i*n//2 + n//2 - 2, 2):
                    maj_swap(j, 1, j + 1, 0)   
            
            for i in range(parity%2,n//2 - parity%2):
                maj_swap(i, 1, i, 0)
            for i in range(parity%2 + n//2,n - parity%2):
                maj_swap(i, 1, i, 0)

            # cirq.barrier()
        def double_exc_layer(parity):
            for i in range(parity % 2, n//2 - 1 - parity % 2, 2):
                append_lad_double_exc(cirq, self.tt, i, i + 1, i + n//2, i + n//2 + 1)

        init_state()
        if n == 4:
            single_exc_layer()
            double_exc_layer(0)
            # for i in range(2):
            #     for j in range(i*n//2, i*n//2 + n//2 - 1, 2):
            #         maj_swap(j, 1, j + 1, 0)
        else:
            for i in range(n//2):
                single_exc_layer()
                double_exc_layer(i)
                mswap_layer(i+1)
        print(sp_maj_exc)
        print(dp_lad_exc)
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
    
    def get_paramaetrized_circuit_generalized(self, reps=1, type=0):
        single_ladder_alpha_exc = self.get_alpha_excitations()
        single_ladder_beta_exc = self.get_beta_excitations()
        single_ladder_exc = single_ladder_alpha_exc + single_ladder_beta_exc
        double_ladder_exc = self.get_double_excitations()
        n = self.n_qubits
        qr = QuantumRegister(n)
        cirq = MyCirq(qr)
            
        self.tt = TernaryTree(n)
        numeration = []
        # numeration_alpha = []
        for i in range(1, n + 1, 2):
            numeration = numeration + [i, i + 1, i + n, i + n + 1]
        for i in range(2, 2*n, 8):
            numeration[i:i+2], numeration[i + 2:i + 4] = numeration[i + 2:i + 4], numeration[i:i+2]
        # print(numeration)
        self.tt.update_branchnum(numeration)
        # print(self.tt)
        k = 0
        dp_lad_exc = copy.deepcopy(self.to_par_ladder_exciations(double_ladder_exc, name="t" + str(k) + "_"))

        sp_lad_exc = copy.deepcopy(self.to_par_ladder_exciations(single_ladder_alpha_exc))
        print(sp_lad_exc)
        sp_lad_exc = copy.deepcopy(self.to_sab_together_ladder_excitations(single_ladder_alpha_exc))
        print(sp_lad_exc)
        sp_maj_exc = self.to_par_maj_exitations_comp(
            copy.deepcopy(self.to_par_maj_excitations(single_ladder_alpha_exc, name="t" + str(k) + "_")),
            n
        )

        def append_single(cirq: MyCirq, tt: TernaryTree, i) -> bool:
            try:
                qubits = [i, i + 1]
                maj_nums = [self.tt[qubit][0].num for qubit in qubits] + [self.tt[qubit][1].num for qubit in qubits]
                list_signs = {"XY": 1, "YX": 1}
                init_orb = min(maj_nums) // 2
                fin_orb = (max(maj_nums) - 1) // 2 
                # psss = [PauliStringSign(*tt.get_majorana(i)) for i in ls]
                # cirq.rotation(reduce(prod_pauli_strings, psss), sp_maj_exc.pop(ls))
                for key in dict_prods_single:
                    psss = [PauliStringSign(*tt.get_majorana(j1 + j2)) for j1,j2 in zip(key, (2*init_orb, 2*fin_orb))]
                    pss = reduce(prod_pauli_strings, psss)
                    if set(qubits) == set(pss.qubits):
                        pauli = sorted(pss.ps.items(), key=lambda x: x[0])
                        pauli = "".join([el[1] for el in pauli])
                        list_signs[pauli] = (list_signs[pauli]*dict_prods_single[key]*pss.sign*(-1j)).real
                    else:
                        return False
                    cirq.single_ex(qubits, sp_lad_exc.pop((fin_orb, init_orb)), list_signs)
                # pss = reduce(prod_pauli_strings, psss)
                # cirq.pauli(pss.pauli_str, pss.qubits)
                return True
            except KeyError:
                return False
            
            
        def append_lad_double_exc(cirq: MyCirq, tt, i) -> bool:
            qubits = [i, i + 1, i + 2, i + 3]
            maj_nums = [self.tt[qubit][0].num for qubit in qubits] + [self.tt[qubit][1].num for qubit in qubits]
            init_orb = min(maj_nums) // 2
            fin_orb = (max(maj_nums) - n - 1) // 2 
            list_signs = {"XXXY": 1, "XXYX": 1, "XYXX": 1, "YXXX": 1, "YYYX": 1, "YYXY": 1, "YXYY": 1, "XYYY": 1}
            # list_signs = {"XXXY": -1, "XXYX": -1, "XYXX": 1, "YXXX": 1, "YYYX": 1, "YYXY": 1, "YXYY": -1, "XYYY": -1}
            if type == 2:
                list_signs = {"ZZZY": -1, "ZZYZ": -1, "ZYZZ": 1, "YZZZ": 1, "YYYZ": -1, "YYZY": -1, "YZYY": 1, "ZYYY": 1}
            if type == 3:
                list_signs = {"ZZZY": -1, "ZZYZ": 1, "ZYZZ": 1, "YZZZ": -1, "YYYZ": 1, "YYZY": -1, "YZYY": -1, "ZYYY": 1}
            try:
                if (fin_orb <= 0) or (init_orb >=n//2):
                    return False
                for key in dict_prods:
                    psss = [PauliStringSign(*tt.get_majorana(j1 + j2)) for j1,j2 in zip(key, (2*init_orb, n + 2*init_orb, 2*fin_orb, n + 2*fin_orb))]
                    pss = reduce(prod_pauli_strings, psss)
                    if set(qubits) == set(pss.qubits):
                        pauli = sorted(pss.ps.items(), key=lambda x: x[0])
                        pauli = pss.ps.items()
                        pauli = "".join([el[1] for el in pauli])
                        list_signs[pauli] = (list_signs[pauli]*dict_prods[key]*pss.sign*(-1j)).real
                    else:
                        return False
                print(list_signs)
                if type==1:
                    cirq.double_ex_yordan(qubits, dp_lad_exc.pop((init_orb + n//2, init_orb, fin_orb + n//2, fin_orb)), list_signs)
                if type==0:
                    cirq.double_ex_opt(qubits, dp_lad_exc.pop((init_orb + n//2, init_orb, fin_orb + n//2, fin_orb)), list_signs)
                if type==2:
                    cirq.double_ex_zyx_opt(qubits, dp_lad_exc.pop((init_orb + n//2, init_orb, fin_orb + n//2, fin_orb)), list_signs)
                if type==3:
                    cirq.double_ex_zyx_yordan(qubits, dp_lad_exc.pop((init_orb + n//2, init_orb, fin_orb + n//2, fin_orb)), list_signs)
            except KeyError:
                # print("HERERERE")
                return False

        def init_state():
            num_electrons = self.num_alpha  
            for i in range(self.n_qubits):
                if (self.tt[i][0].num <= 2*num_electrons) or (self.num_spatial_orbitals * 2 + 1<= 
                        self.tt[i][0].num <= 2*self.num_spatial_orbitals + 2*num_electrons):
                    cirq.x(i)
                    # TODO
                else:
                    cirq.id(i)
                    
        def single_exc_layer(parity):
            for i in range(0, n-1, 2):
                append_single(cirq, self.tt, i)
            for i in range(1, n-1, 2):
                append_single(cirq, self.tt, i)

        def fswap(fq, cirq=cirq):
            list_trans = [(fq, 0, fq + 1, 0), (fq, 1, fq + 1, 1)]
            _ = [self.tt.branch_transposition(*el, unsigned=True) for el in list_trans]  
            cirq.fswap([fq, fq + 1])

        def fswap_layer():
            # cirq.barrier()
            for i in range(0, n-1, 2):
                fswap(i)
            for i in range(1, n-1, 2):
                fswap(i)

            # cirq.barrier()

        def double_exc_layer(parity):
            if (parity % 2 == 0):
                for i in range(0, n-3,4):
                    append_lad_double_exc(cirq, self.tt, i)
            else:
                for i in range(2, n-3,4):
                    append_lad_double_exc(cirq, self.tt, i)
                
            for i in range(0, n-3,4):
                append_lad_double_exc(cirq, self.tt, i)

        init_state()
        if n == 4:
            pass
            single_exc_layer(0)
            double_exc_layer(0)
            # for i in range(2):
            #     for j in range(i*n//2, i*n//2 + n//2 - 1, 2):
            #         maj_swap(j, 1, j + 1, 0)
        else:
            for i in range(n//2):
                single_exc_layer(i)
                double_exc_layer(i)
                fswap_layer()
        print(sp_lad_exc)
        print(dp_lad_exc)
        print(self.tt)
        return cirq