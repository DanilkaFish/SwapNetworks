from __future__ import annotations

import itertools
from . import TernaryTree, prod_pauli_strings
from .AbstractUCC import AbstractUCC
from .Mycirq import MyCirq
import copy
from qiskit import QuantumRegister

class UpUCCSDG(AbstractUCC):
    def __init__(self,  **kwargs):
        super().__init__(**kwargs)
        self.tt = TernaryTree(self.num_spin_orbitals)

    def get_alpha_excitations(self) -> list[tuple[int, int]]:
        """
        Method to get spin preserving single excitation operators via ladder
        operators for alpha orbitals
        """
        # alpha_occ = range(self.num_alpha)
        # alpha_unocc = range(self.num_alpha, self.num_spatial_orbitals)
        # return list(itertools.product(alpha_unocc, alpha_occ))
        alpha = range(self.num_spatial_orbitals)
        return [op for op in itertools.product(alpha, alpha) if op[0]>op[1]]
        # return list(itertools.product(alpha, alpha))

    def get_beta_excitations(self) -> list[tuple[int, int]]:
        """
        Method to get spin preserving single excitation operators via ladder
        operators for beta orbitals
        """
        # beta_index_offset = self.num_spatial_orbitals
        # beta_occ = range(beta_index_offset, beta_index_offset + self.num_beta)
        # beta_unocc = range(beta_index_offset + self.num_beta, self.num_spin_orbitals)
        
        # return list(itertools.product(beta_unocc, beta_occ))
        beta = range(self.num_spatial_orbitals, 2*self.num_spatial_orbitals)
        
        return [op for op in itertools.product(beta, beta) if op[0]>op[1]]

    def get_double_excitations(self):
        """
        Method to get spin preserving pair double excitation operators via ladder
        operators
        """
        alpha = range(self.num_spatial_orbitals)
        beta_index_offset = self.num_spatial_orbitals
        # ab = list(itertools.product(alpha, alpha))
        ab = [(x,y) for x in alpha for y in alpha if x >= y]
        ab = [(el[0] + beta_index_offset, el[0], el[1] + beta_index_offset, el[1]) for el in ab]
        return ab              
    

    def get_parametrized_circuit(self, number_of_layers=1):
        n = self.n_qubits
        qr = QuantumRegister(n)
        cirq = MyCirq(qr)
        self.tt = TernaryTree(n) # JW tree
        maj_exc_par = None
        numeration = [None]*n*2

        for i in range(n//4):
            numeration[4*i: 4*i + 4] = [2*i + 1, 2*i + 2, 2*i + n//2 + 1, 2*i + n//2 + 2]
        for i in range(n//2, n//2 + n//4):
            numeration[4*i - n: 4*i - n + 4] = [2*i + 1, 2*i + 2, 2*i + n//2 + 1, 2*i + n//2 + 2]

        def append_maj_exc(cirq, tt, ls: tuple) -> bool:
            if ls in maj_exc_par:
                coef = 1
                maj = [tt.get_majorana(i) for i in ls]
                for el in maj:
                    coef = coef * el[1]
                maj = tuple([el[0] for el in maj])
                prod = maj[0]

                for pauli in maj[1:]:
                    prod, sign = prod_pauli_strings(prod, pauli)
                    coef = coef * sign

                prod = ((gate[0], gate[1]) for gate in prod)
                cirq.pauli(prod, coef * maj_exc_par[ls])
                maj_exc_par.pop(ls)
                return True
            return False
        
        def transposition(node1, edge1, node2, edge2, cirq=cirq):
            cirq.pauli(self.tt.branch_transposition(node1, edge1, node2, edge2))
        
        def init_state():
            for i in range(n):
                cirq.id(i)
            for i in range(n):
                if ((i < self.num_alpha) or (0 <= (i - n//2)< self.num_beta)):
                    cirq.x(i)
                else:
                    cirq.id(i)

        def single_prep():
            res = 0
            disp = 0
            for k in range(n//2,0,-2):
                for i in range(0, k//2):
                    transposition(2*i + disp, res, 2*i + 1 + disp, (res+1)%2)
                    transposition(2*i + n//2 + disp, res, 2*i + n//2 + 1  + disp, (res+1)%2)
                for i in range(0, (k-1)//2):
                    transposition(2*i + disp + 1, res, 2*i + 2 + disp, (res+1)%2)
                    transposition(2*i + disp + n//2 + 1, res, 2*i + n//2 + 2 + disp, (res+1)%2)
                res = (res + 1)%2
                disp += 1

        def single_exc():
            def _maj():
                for i in range(n):
                    ls = sorted([self.tt[i][0].num, self.tt[i][1].num])
                    append_maj_exc(cirq, self.tt, tuple(ls))
            def ferm_trans(i, j):
                _maj(), transposition(i,0,j,1)
                _maj(), transposition(i,1,j,1)
                _maj(), transposition(i,1,j,0)
            _maj()
            for block in range(4):
                for i in range(n//4):
                    j = i % 2
                    while j + 1 < n//4:
                        ferm_trans(j + block*n//4, j + block*n//4 + 1)
                        j = j + 2
        
        def single_anti_prep():
            res = 1
            disp = n//4 - 1
            for k in range(2,n//2, 2):
                for i in range(0, (k-1)//2):
                    transposition(2*i + disp + 1, res, 2*i + 2 + disp, (res+1)%2)
                    transposition(2*i + disp + n//2 + 1, res, 2*i + n//2 + 2 + disp, (res+1)%2)
                for i in range(0, k//2):
                    transposition(2*i + disp, res, 2*i + 1 + disp, (res+1)%2)
                    transposition(2*i + n//2 + disp, res, 2*i + n//2 + 1  + disp, (res+1)%2)
                res = (res + 1)%2
                disp -= 1
            for i in range(0, (n//2-1)//2):
                transposition(2*i + disp + 1, res, 2*i + 2 + disp, (res+1)%2)
                transposition(2*i + disp + n//2 + 1, res, 2*i + n//2 + 2 + disp, (res+1)%2)

        def double_exc_tree():
            def tr(i, l):
                if l == 2*i:
                    transposition(n//2 - l + 2*i - 1, 1, n//2 - l + 2*i, 0)
                elif 2*i < l:
                    transposition(n//2 - l + 2*i - 1, 1, n//2 - l + 2*i, 1)
                else:
                    transposition(n//2 - l + 2*i - 1, 0, n//2 - l + 2*i, 0)

            for l in range(n // 2):
                for i in range(l + 1):
                    tr(i, l)
            for l in range(n // 2 - 2, -1 , -1):
                for i in range(l + 1):
                    tr(i, l)
            # cirq.barrier()

        def double_exc_layer():
            def inv(i):
                pauli = self.tt.branch_transposition(2*i, 1, 2*i + 1, 1)
                cirq.pauli(pauli)  
            def append_maj(first, second, ancilla_cirq):
                ls = sorted([self.tt[first][0].num, self.tt[first][1].num,
                                self.tt[second][0].num, self.tt[second][1].num])
                return int(append_maj_exc(ancilla_cirq, self.tt, tuple(ls)))

            k = 1
            parities = [0 for _ in range(n//2)]
            while k < n//2:
                i = 0
                while i < n//2:
                    for _ in range(k):
                        if i < n//2 and parities[i] == 0:
                            inv(i)
                            parities[i] = 1
                        i += 1
                    for _ in range(k):
                        if i < n//2 and parities[i] == 1:
                            inv(i)
                            parities[i] = 0 
                        i += 1

                zero_pos = [0 for _ in range(n//2)]
                one_pos = [0 for _ in range(n//2)]
                for _ in range(2):
                    i = 0
                    j = 0
                    for index, par in enumerate(parities):
                        if par == 0:
                            zero_pos[i] = index
                            i += 1
                        else:
                            one_pos[j] = index
                            j += 1  
                    for dist in range(n//4):
                        ancilla_cirq = MyCirq(n)
                        counts = 0
                        for r in range(n//2):
                            first = 2*one_pos[r % j]
                            second = 2*zero_pos[(r + dist) % (i)]
                            counts += append_maj(first, second, ancilla_cirq)
                            counts += append_maj(first + 1, second + 1, ancilla_cirq)
                            counts += append_maj(first, second + 1, ancilla_cirq)
                            counts += append_maj(first + 1, second, ancilla_cirq)

                        ancilla_cirq.name = str(counts) + " 2-gate"
                        if counts > 0:
                            cirq.compose(ancilla_cirq, inplace=True, wrap=True)
                    for i in range(n//2):
                        inv(i)
                        parities[i] = (parities[i] + 1) % 2

                    i += 1
                k = k * 2
            for i in range(n//2):
                if parities[i] != 0:
                    inv(i)
        init_state()
        single_prep()
        # print("init_state", self.tt)
        for k in range(number_of_layers):
            maj_exc_par = copy.deepcopy(self.get_excitations(name="t" + str(k) + "_"))
            single_exc()
            single_anti_prep()
            # print("before_double", self.tt)

            double_exc_tree()
            double_exc_layer()
            print(maj_exc_par)
            if k < number_of_layers-1:
                double_exc_tree()
                single_anti_prep()

            # print("after_double", self.tt)
            # print("------------------", maj_exc_par)
            # single_anti_prep()
        # print(cirq)
        return cirq
