from itertools import combinations


dict_prod_coef = {"II" : ["I",1] , "XX" :  ["I",1] , "YY" :  ["I",1] ,"ZZ" :  ["I",1] ,
             "XY" :  ["Z", 1j] ,"YX" : ["Z", -1j],"XZ" : ["Y",-1j], "ZX" : ["Y",1j],"YZ" : ["X",1j],"ZY" : ["X",-1j],
            "IX" : ["X",1], "XI" : ["X",1], "YI" : ["Y",1],"IY" : ["Y",1],"IZ" : ["Z",1],"ZI" : ["Z",1]}

dict_cx_trans = {"II": "II", "IX": "IX", "ZI": "ZI", "ZX": "ZX", "XI": "XX", "XX": "XI", "IZ": "ZZ", "ZZ": "IZ", "YI": "YX", "YX": "YI", "XZ": "YY", "YY": "XZ", 
                 "XY": "YZ", "YZ": "XY", "ZY": "IY", "IY":"ZY"}
def prod(pauli1, pauli2):
    ops = ""
    coeff = 1
    for el1, el2 in zip(pauli1, pauli2):
        op, coef = dict_prod_coef[el1 + el2]
        coeff *= coef
        ops += op
    return ops, coeff

def weight(pauli: str):
    i = 0
    for s in pauli:
        if s !="I":
            i += 1
    return i

class clifford_generator:
    def __init__(self,
                 n=1,
                 gen_dict=None):
        
        if gen_dict is None:
            self.gen_dict = self.get_jw_majorana(n)
        else:
            self.gen_dict = gen_dict
        self.prods = {i: (i,) for i in range(2*n)}

    def get_jw_majorana(self, n):
        gen_dict = {}
        for i in range(n):
            pauliX, pauliY = "", ""
            for _ in range(i):
                pauliX += "Z"
                pauliY += "Z"
            pauliX, pauliY = pauliX + "X", pauliY + "Y"

            for _ in range(i + 1, n):
                pauliX += "I"
                pauliY += "I"

            gen_dict[2*i],  gen_dict[2*i + 1] = pauliX, pauliY
        return gen_dict

        
    def get_full_clifford_basis(self):
        prods = {}
        maj_tuple = tuple([i for i in range(len(self))])
        for i in maj_tuple:
            for els in combinations(maj_tuple, i + 1):
                op = "I" * (len(self) // 2)
                for el in els:
                    op = prod(op, self.gen_dict[el])[0]
                prods[els] = op
        return prods
    
    def transform_maj(self, pauli: str, num=None):
        c = 1
        if pauli != "cx":
            for key in self.gen_dict:
                ops, coef = prod(self.gen_dict[key], pauli)
                c = coef*c
                if coef.real == 0:
                    self.gen_dict[key] = ops
        else:
            for key in self.gen_dict:
                pauli = self.gen_dict[key]
                tr = dict_cx_trans[pauli[num[0]] + pauli[num[1]]]
                pauli = list(pauli)
                pauli[num[0]] = tr[0]
                pauli[num[1]] = tr[1]
                self.gen_dict[key] = "".join(pauli)

    def __len__(self):
        return len(self.gen_dict)

    def __str__(self):
        s = ''
        # for i in range(len(self) // 2):
        #     s += f"({2*i}, {2*i + 1}): {self.gen_dict[2*i]}, {self.gen_dict[2*i + 1]}\n"
        for key in self.gen_dict:
            s += f"({key}): {self.gen_dict[key]}\n"
        return s
    
    def __iter__(self):
        return self.gen_dict.__iter__()


def get_pauli(char='I', pos=0, numq=4):
    pauli = ["I" for _ in range(numq)]
    pauli[pos] = char
    # print("".join(pauli))
    return "".join(pauli)
 
def check_anticom(pauli1, pauli2):
    # print(f"{prod(pauli1, pauli2)[1]} : {prod(pauli2, pauli1)[1]}")
    if prod(pauli1, pauli2)[1] != prod(pauli2, pauli1)[1]:
        return True
    return False


def num_basis(n):
    cg = clifford_generator(n)
    pauli = cg.get_full_clifford_basis()
    counts = 0
    for item in combinations(pauli, 2*n):
        flag = True
        for el in combinations(item, 2):
            if not check_anticom(pauli[el[0]], pauli[el[1]]):
                flag = False
        if flag:
            counts += 1
            print(item)
    print(counts)
def num_anticom(n):
    cg = clifford_generator(n)
    pauli = cg.get_full_clifford_basis()
    pauli_0  = (0,1,2,3)
    counts = 0
    counts2 = 0
    for pauli_1 in pauli:
        for pauli_4 in pauli:
            counts2 = 0
            for pauli_3 in pauli:
                counts = 0
                for pauli_2 in pauli:
                    flag = True
                    item = [pauli_0, pauli_1, pauli_2, pauli_3, pauli_4]
                    for el in combinations(item, 2):
                        if not check_anticom(pauli[el[0]], pauli[el[1]]):
                            flag = False
                    if flag:
                        counts += 1
                if counts > 0 :
                    counts2 += 1
                    print(counts)
            if counts2 > 0 :
                print(counts2)
