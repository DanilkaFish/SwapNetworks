
def get_parametrized_circuit_w(self):
    n = self.num_spin_orbitals
    N = n//2
    cirq = MyCirq(n)
    for i in range(n):
        cirq.id(i)
    for i in range(self.num_alpha):
        cirq.x(i + 1)
    for i in range(self.num_alpha,N - self.num_alpha + 1):
        cirq.id(i)
    for i in range(self.num_alpha):
        cirq.x(i + N + 1)
    for i in range(self.num_alpha,N - self.num_alpha + 1):
        cirq.id(i + N)
    # cirq.s(3)

    qubs = [i for i in range(n-1,-1, -1)]
    self.tt = TernaryTree(self.num_spin_orbitals)

    def append_maj_exc(maj_exc_par, tt, ls: tuple):
        coef = 1
        maj = [tt.get_majorana(i) for i in ls]
        for el in maj:
            coef = coef * el[1]
        maj = [el[0] for el in maj]

        maj = tuple(maj)
        if ls in maj_exc_par:
            prod = maj[0]
            for pauli in maj[1:]:
                prod, sign = prod_pauli_strings(prod, pauli)
                coef = coef * sign
            prod = ((qubs.index(gate[0]), gate[1]) for gate in prod)
            cirq.pauli(prod, qubs, coef * maj_exc_par[ls])
            # maj_exc_par.pop(ls)
        # single layer
    for i in range(1):
        for maj in self.maj_exc_par:
            if len(maj) == 2:
                append_maj_exc(self.maj_exc_par, self.tt, maj)  
        for maj in self.maj_exc_par:
            if len(maj) == 4:
                append_maj_exc(self.maj_exc_par, self.tt, maj) 
        self.maj_exc_par = self.get_excitations(name=str(i))
    cirq.s(3)
    
    return cirq  

    
def get_real_parametrized_cirquit(self):
    n = self.num_spin_orbitals
    cirq = MyCirq(n)
    self.tt = TernaryTree(n) # JW tree
    numeration = [None]*n*2
    for i in range(n//4):
        numeration[4*i: 4*i + 4] = [2*i + 1, 2*i + 2, 2*i + n//2 + 1, 2*i + n//2 + 2]
    for i in range(n//2, n//2 + n//4):
        numeration[4*i - n: 4*i - n + 4] = [2*i + 1, 2*i + 2, 2*i + n//2 + 1, 2*i + n//2 + 2]

    self.tt.update_branchnum(numeration)
    maj_exc_par = copy.deepcopy(self.maj_exc_par)

    def append_maj_exc(cirq, tt, ls: tuple):
        coef = 1
        maj = [tt.get_majorana(i) for i in ls]
        for el in maj:
            coef = coef * el[1]
        maj = [el[0] for el in maj]
        maj = tuple(maj)
        if ls in maj_exc_par:
            prod = maj[0]
            for pauli in maj[1:]:
                prod, sign = prod_pauli_strings(prod, pauli)
                coef = coef * sign
            prod = ((gate[0], gate[1]) for gate in prod)
            cirq.pauli(prod, coef * maj_exc_par[ls])
            maj_exc_par.pop(ls)
            return True
        return False
    
    def init_state():
        for i in range(n):
            if (i % 2 == 0) and ((i // 2 < self.num_alpha) or (n // 2<= i//2 < n // 2 + self.num_beta)):
                cirq.x(i)
            else:
                cirq.id(i)
        # for i in range(self.num_alpha):
        #     cirq.x(2 * i)
        # for i in range(self.num_spatial_orbitals - self.num_alpha):
        #     cirq.id(2 * i + 1)
        # for i in range(self.num_beta):
        #     cirq.x(2 * i + self.num_spatial_orbitals)
        # for i in range(self.num_spatial_orbitals - self.num_beta ):
        #     cirq.id(2 * i + 1 + self.num_spatial_orbitals)

    def single_exc_tree():
        for i in range(n // 4):
            pauli = self.tt.branch_transposition(2 * i, 1, 2 * i + 1, 0)
            cirq.pauli_str(pauli)
            pauli = self.tt.branch_transposition(n//2 + 2*i, 1, n//2 + 2*i + 1, 0)
            cirq.pauli_str(pauli)

    def single_exc_layer():
        for _ in range(n//2):
            for i in range(n):
                ls = sorted([self.tt[i][0].num, self.tt[i][1].num])
                append_maj_exc(cirq, self.tt, tuple(ls))

            for i in range(n//4):
                pauli = self.tt.branch_transposition(2 * i, 0, 2 * i + 1, 0)
                cirq.pauli_str(pauli)
                pauli = self.tt.branch_transposition(n//2 + 2*i, 0, n//2 + 2*i + 1, 0)
                cirq.pauli_str(pauli)

            for i in range((n - 2)//4 ):
                pauli = self.tt.branch_transposition(2 * i + 1, 0, 2 * i + 2, 0)
                cirq.pauli_str(pauli)
                pauli = self.tt.branch_transposition(n//2 + 2*i + 1, 0, n//2 + 2*i + 2, 0)
                cirq.pauli_str(pauli)

    def double_exc_tree():
        for l in range(n // 2):
            for i in range(l + 1):
                if l == 2*i:
                    pauli = self.tt.branch_transposition(n//2 - l + 2*i - 1, 1, n//2 - l + 2*i, 0)
                    cirq.pauli_str(pauli)
                elif 2*i < l:
                    pauli = self.tt.branch_transposition(n//2 - l + 2*i - 1, 1, n//2 - l + 2*i, 1)
                    cirq.pauli_str(pauli)
                else:
                    pauli = self.tt.branch_transposition(n//2 - l + 2*i - 1, 0, n//2 - l + 2*i, 0)
                    cirq.pauli_str(pauli)

        for l in range(n // 2 - 2, -1 , -1):
            for i in range(l + 1):
                if l == 2*i:
                    pauli = self.tt.branch_transposition(n//2 - l + 2*i - 1, 1, n//2 - l + 2*i, 0)
                    cirq.pauli_str(pauli)
                elif 2*i < l:
                    pauli = self.tt.branch_transposition(n//2 - l + 2*i - 1, 1, n//2 - l + 2*i, 1)
                    cirq.pauli_str(pauli)
                else:
                    pauli = self.tt.branch_transposition(n//2 - l + 2*i - 1, 0, n//2 - l + 2*i, 0)
                    cirq.pauli_str(pauli)
    
    def double_exc_layer():
        def inv(i):
            pauli = self.tt.branch_transposition(2*i, 1, 2*i + 1, 1)
            cirq.pauli(pauli)  

        k = 1
        parities = [0 for _ in range(n//2)]
        while k <= n // 4 :
            i = 0
            while i < n // 2:
                for _ in range(k):
                    if i < n//2 and parities[i] != 1:
                        inv(i)
                        parities[i] = 1
                    i += 1
                for _ in range(k):
                    if i < n//2 and parities[i] != 0:
                        inv(i)
                        parities[i] = 0 
                    i += 1

            print(parities)
            zero_pos = [0 for _ in range(n//4)]
            one_pos = [0 for _ in range(n//4)]
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
                    for r in range(n//4):
                        first = 2*one_pos[r]
                        second = 2*zero_pos[(r + dist) % (n//4)]
                        ls = sorted([self.tt[first][0].num, self.tt[first][1].num,
                                        self.tt[second][0].num, self.tt[second][1].num])
                        flag = append_maj_exc(ancilla_cirq, self.tt, tuple(ls))
                        if flag:
                            counts += 1
                        ls = sorted([self.tt[first + 1][0].num, self.tt[first + 1][1].num,
                                        self.tt[second][0].num, self.tt[second][1].num])
                        flag = append_maj_exc(ancilla_cirq, self.tt, tuple(ls))
                        if flag:
                            counts += 1
                        ls = sorted([self.tt[first][0].num, self.tt[first][1].num,
                                        self.tt[second + 1][0].num, self.tt[second + 1][1].num])
                        flag = append_maj_exc(ancilla_cirq, self.tt, tuple(ls))
                        if flag:
                            counts += 1
                        ls = sorted([self.tt[first+1][0].num, self.tt[first+1][1].num,
                                        self.tt[second+1][0].num, self.tt[second+1][1].num])
                        flag = append_maj_exc(ancilla_cirq, self.tt, tuple(ls))
                        if flag:
                            counts += 1

                        print(first, second)
                    ancilla_cirq.name = str(counts) + " 2-gate"
                    cirq.compose(ancilla_cirq, inplace=True, wrap=True)
                for i in range(n//2):
                    inv(i)
                    parities[i] = (parities[i] + 1) % 2

            k = k * 2

        
    init_state()
    single_exc_tree()
    single_exc_layer()
    double_exc_tree()
    double_exc_layer()
    print(maj_exc_par)
    return cirq

def get_parametrized_circuit(self):
    n = self.num_spin_orbitals
    cirq = MyCirq(n)

    qubs = [i for i in range(n - 1, -1, -1)]
    self.tt = TernaryTree(self.num_spin_orbitals)
    
    def append_maj_exc(maj_exc_par, tt, ls: tuple):
        coef = 1
        maj = [tt.get_majorana(i) for i in ls]
        for el in maj:
            coef = coef * el[1]
        maj = [el[0] for el in maj]

        maj = tuple(maj)
        if ls in maj_exc_par:
            prod = maj[0]
            for pauli in maj[1:]:
                prod, sign = prod_pauli_strings(prod, pauli)
                coef = coef * sign
            prod = ((qubs.index(gate[0]), gate[1]) for gate in prod)
            cirq.pauli(prod, qubs, coef * maj_exc_par[ls])
            maj_exc_par.pop(ls)

    def single_prep1():
        for i in range(n // 4):
            pauli = self.tt.branch_transposition(2 * i, 0, 2 * i + 1, 0)
            cirq.pauli(pauli, qubs)
            pauli = self.tt.branch_transposition(n // 2 + 2 * i, 0, n // 2 + 2 * i + 1, 0)
            cirq.pauli(pauli, qubs)

    def single_prep2():
        for i in range((n + 2) // 4 - 1):
            pauli = self.tt.branch_transposition(2 * i + 1, 0, 2 * i + 2, 0)
            cirq.pauli(pauli, qubs)
            pauli = self.tt.branch_transposition(n // 2 + 2 * i + 1, 0, n // 2 + 2 * i + 2, 0)
            cirq.pauli(pauli, qubs)

    def single_layer(maj_exc_par, tt):
        for i in range(n):
            ls = sorted([tt[qubs[i]][0].num, tt[qubs[i]][1].num])
            append_maj_exc(maj_exc_par, tt, tuple(ls))

    def double_preparation():
        num = []
        for el in [[2 * i + 1, 2 * i + 2, 2 * i + n + 1, 2 * i + n + 2] for i in range(n // 2)]:
            num = num + el
        # TODO
        for l in range(n // 2 - 1):
            for i in range(l + 1):
                pauli = self.tt.branch_transposition(n // 2 - l + 2 * i - 1, 0, n // 2 - l + 2 * i, 0)
                cirq.pauli(pauli, qubs)
                pauli = self.tt.branch_transposition(n // 2 - l + 2 * i - 1, 1, n // 2 - l + 2 * i, 1)
                cirq.pauli(pauli, qubs)

        for i in range(self.tt.n_qubits // 2):
            pauli = self.tt.branch_transposition(2 * i, 1, 2 * i + 1, 0)
            cirq.pauli(pauli, qubs)

    def double_antipreparation():
        num = []
        for el in [[2 * i + 1, 2 * i + 2, 2 * i + n + 1, 2 * i + n + 2] for i in range(n // 2)]:
            num = num + el

        for i in reversed(range(self.tt.n_qubits // 2)):
            pauli = self.tt.branch_transposition(2 * i, 1, 2 * i + 1, 0)
            cirq.pauli(pauli, qubs)

        # TODO
        for l in reversed(range(n // 2 - 1)):
            for i in reversed(range(l + 1)):
                pauli = self.tt.branch_transposition(n // 2 - l + 2 * i - 1, 0, n // 2 - l + 2 * i, 0)
                cirq.pauli(pauli, qubs=qubs)
                pauli = self.tt.branch_transposition(n // 2 - l + 2 * i - 1, 1, n // 2 - l + 2 * i, 1)
                cirq.pauli(pauli,qubs=qubs)

    def swap():
        for i in range(n // 2):
            qubs[2 * i], qubs[2 * i + 1] = qubs[2 * i + 1], qubs[2 * i]
            cirq.swap(2 * i, 2 * i + 1)
        for i in range(n // 2 - 1):
            qubs[2 * i + 1], qubs[2 * i + 2] = qubs[2 * i + 2], qubs[2 * i + 1]
            cirq.swap(2 * i + 1, 2 * i + 2)

    def double_layer(maj_exc_par, tt):
        for i in range(n // 2):
            ls = sorted([tt[qubs[2 * i + 1]][1].num, tt[qubs[2 * i + 1]][0].num,
                            tt[qubs[2 * i]][1].num, tt[qubs[2 * i]][0].num])
            append_maj_exc(maj_exc_par, tt, tuple(ls))

        for i in range(n // 2 - 1):
            ls = sorted([tt[qubs[2 * i + 2]][0].num, tt[qubs[2 * i + 2]][1].num,
                            tt[qubs[2 * i + 1]][0].num, tt[qubs[2 * i + 1]][1].num])
            append_maj_exc(maj_exc_par, tt, tuple(ls))

    N = self.mol.num_spatial_orbitals
    for l in range(3):
        for par in self.maj_exc_par:
            # print(par," : ", self.maj_exc_par)
            pass
        maj_exc_par = copy.deepcopy(self.maj_exc_par)
        for i in range(n):
            cirq.id(i)
        for i in range(self.num_alpha):
            cirq.x(qubs[i])
        for i in range(self.num_alpha,N - self.num_alpha + 1):
            cirq.id(qubs[i])
        for i in range(self.num_alpha):
            cirq.x(qubs[i + N])
        for i in range(self.num_alpha,N - self.num_alpha + 1):
            cirq.id(qubs[i + N])
        # qubs = list(reversed(qubs))
        for _ in range(n // 2):
            single_layer(maj_exc_par, self.tt)
            single_prep1()
            single_prep2()

        double_preparation()

        for _ in range(n // 2):
            double_layer(maj_exc_par, self.tt)
            swap()

        for i in range(n // 2):
            pauli = self.tt.branch_transposition(qubs[2 * i], 1, qubs[2 * i + 1], 1)
            cirq.pauli(pauli, qubs)

        for _ in range(n // 2):
            double_layer(maj_exc_par, self.tt)
            swap()
        double_antipreparation()
    return cirq


    def get_parametrized_circuit(self):
            n = self.n_qubits
            qr = QuantumRegister(n)
            cirq = MyCirq(qr)
            self.tt = TernaryTree(n) # JW tree
            maj_exc_par = copy.deepcopy(self.maj_exc_par)
            numeration = [None]*n*2

            for i in range(n//4):
                numeration[4*i: 4*i + 4] = [2*i + 1, 2*i + 2, 2*i + n//2 + 1, 2*i + n//2 + 2]
            for i in range(n//2, n//2 + n//4):
                numeration[4*i - n: 4*i - n + 4] = [2*i + 1, 2*i + 2, 2*i + n//2 + 1, 2*i + n//2 + 2]

            self.tt.update_branchnum(numeration)

            def append_maj_exc(cirq, tt, ls: tuple) -> bool:
                coef = 1
                maj = [tt.get_majorana(i) for i in ls]
                for el in maj:
                    coef = coef * el[1]
                maj = tuple([el[0] for el in maj])

                if ls in maj_exc_par:
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
                    if (i % 2 == 0) and ((i // 2 < self.num_alpha) or (0 <= (i - n//2)//2 < self.num_beta)):
                        cirq.x(i)
                    else:
                        cirq.id(i)

            def single_exc_tree():
                for i in range(n // 4):
                    transposition(2 * i, 1, 2 * i + 1, 0)
                    transposition(n//2 + 2 * i, 1, n//2 + 2 * i + 1, 0)
                
            def single_exc_layer():
                def _maj():
                    for i in range(n):
                        ls = sorted([self.tt[i][0].num, self.tt[i][1].num])
                        append_maj_exc(cirq, self.tt, tuple(ls))

                def move2_forward(pair, layer=0):
                    def tr(i):
                        transposition(pair[i], layer, pair[i] + 1, layer)
                        pair[i] += 1
                    tr(1), tr(0), tr(1), tr(0)

                def move2_back(pair, layer=0):
                    def tr(i):
                        transposition(pair[i] - 1, layer, pair[i], layer)
                        pair[i] -= 1
                    tr(0), tr(1), tr(0), tr(1)

                def forward(init=0, layer=0, num=self.num_alpha):
                    _maj()
                    nodes = [[init, init + 1]]
                    for _ in range(n//4):
                        if len(nodes) < num:
                            if nodes[0][0] == 4 + init:
                                nodes.insert(0, [init, init + 1])
                        if nodes[-1][-1] == n // 2 - 1 + init:
                            nodes.pop(-1)
                            num -= 1
                        for pair in nodes:
                            move2_forward(pair, layer)
                        _maj()

                def back(init=0, layer=0, num=self.num_alpha):
                    nodes = [[n//2 - 2 + init, n//2 - 1 + init]]
                    for _ in range(n//4):
                        if len(nodes) < num:
                            if nodes[-1][0] == n//2 - 6 + init:
                                nodes.append([n//2 - 2 + init, n//2 - 1 + init])
                        if nodes[0][0] == init:
                            nodes.pop(0)
                            num -= 1
                        for pair in nodes:
                            move2_back(pair, layer)
                        _maj()

                forward(0), forward(n//2, num=self.num_beta)
                for i in range(n//2 - self.num_alpha*2):
                    transposition(i, 0, i, 1)
                for i in range(n//2, n - self.num_beta*2):
                    transposition(i, 0, i, 1)

                back(0), back(n//2, num=self.num_beta)
                forward(0, 1, self.num_alpha)
                forward(n//2, 1, self.num_beta)
                for i in range(self.num_alpha*2, n//2):
                    transposition(i, 0, i, 1)
                for i in range(n//2 + self.num_beta*2, n):
                    transposition(i, 0, i, 1)

            def double_exc_tree():
                def tr(i, l):
                    if l == 2*i:
                        transposition(n//2 - l + 2*i - 1, 1, n//2 - l + 2*i, 0)
                    elif 2*i < l:
                        transposition(n//2 - l + 2*i - 1, 1, n//2 - l + 2*i, 1)
                    else:
                        transposition(n//2 - l + 2*i - 1, 0, n//2 - l + 2*i, 0)

                cirq.barrier()
                for l in range(n // 2):
                    for i in range(l + 1):
                        tr(i, l)
                for l in range(n // 2 - 2, -1 , -1):
                    for i in range(l + 1):
                        tr(i, l)
                cirq.barrier()

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
                        print("(i, j) = ", (i,j))
                        print("parities = ", parities)
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
            single_exc_tree()
            single_exc_layer()
            double_exc_tree()
            double_exc_layer()
            print(maj_exc_par)

            return cirq


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