from itertools import combinations


dict_prod_coef = {"II" : ["I",1] , "XX" :  ["I",1] , "YY" :  ["I",1] ,"ZZ" :  ["I",1] ,
             "XY" :  ["Z",1j] ,"YX" : ["Z",-1j],"XZ" : ["Y",-1j], "ZX" : ["Y",1j],"YZ" : ["X",1j],"ZY" : ["X",-1j],
            "IX" : ["X",1], "XI" : ["X",1], "YI" : ["Y",1],"IY" : ["Y",1],"IZ" : ["Z",1],"ZI" : ["Z",1]}

dict_cx_trans = {"II": ["II", 1], "IX": ["IX", 1],"ZI": ["ZI", 1],"ZX": ["ZX", 1],"XI": ["XX", 1],"XX": ["XI", 1],
                 "IZ": ["ZZ", 1],"ZZ": ["IZ", 1],"YI": ["YX", 1],"YX": ["YI", 1],"XZ": ["YY", -1],"YY": ["XZ", -1],
                 "XY": ["YZ", 1],"YZ": ["XY", 1],"ZY": ["IY", 1],"IY": ["ZY", 1]}

dict_cz_trans = {"II": ["II", 1],"IZ": ["IZ", 1],"ZI": ["ZI", 1],"ZZ": ["ZZ", 1],"XX": ["YY", 1],"YY": ["XX", 1],"XY": ["YX", 1],"YX": ["XY", 1],"XZ": ["XI", 1],"XI": ["XZ", 1],"ZX": ["IX", 1],"IX": ["ZX", 1],
                 "YI": ["YZ", 1],"YZ": ["YI", 1],"ZY": ["IY", 1],"IY": ["ZY", 1]}
dict_h_trans = {"I": ["I", 1], "X": ["Z", 1], "Y": ["Y", -1], "Z": ["X", 1],}
dict_x_trans = {"I": ["I", 1], "X": ["X", 1], "Y": ["Y", -1], "Z": ["Z", -1],}
                 
def prod(pauli1, pauli2):
    ops = ""
    coeff = 1j
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

    @property
    def gen_dict(self):
        return self._gen_dict
    
    @gen_dict.setter
    def gen_dict(self, smth):
        self._gen_dict = {}
        for key in smth:
            if isinstance(smth[key], list):
                self._gen_dict[key] = smth[key]
            else:
                self._gen_dict[key] = [smth[key], 1]


    def get_jw_majorana(self, n):
        gen_dict = {}
        for i in range(n):
            pauliX, pauliY = "",""
            for _ in range(i):
                pauliX += "Z"
                pauliY += "Z"
            pauliX, pauliY = pauliX + "X",pauliY + "Y"

            for _ in range(i + 1, n):
                pauliX += "I"
                pauliY += "I"

            gen_dict[2*i],  gen_dict[2*i + 1] = [pauliX, 1], [pauliY, 1]
        return gen_dict

        
    # def get_full_clifford_basis(self):
    #     prods = {}
    #     maj_tuple = tuple([i for i in range(len(self))])
    #     for i in maj_tuple:
    #         for els in combinations(maj_tuple, i + 1):
    #             op = "I" * (len(self) // 2)
    #             for el in els:
    #                 op = prod(op, self.gen_dict[el])[0]
    #             prods[els] = op
    #     return prods
    

    def transform_maj(self, pauli: str, num=None, coef=1):
        c = coef
        if pauli == "cx":
            # num = list(reversed(num))
            for key in self.gen_dict:
                pauli = self.gen_dict[key][0]
                tr, coef = dict_cx_trans[pauli[num[0]] + pauli[num[1]]]
                pauli = list(pauli)
                pauli[num[0]] = tr[0]
                pauli[num[1]] = tr[1]
                self.gen_dict[key] = ["".join(pauli), c*coef*self.gen_dict[key][1]]
        elif pauli == "cz":
            for key in self.gen_dict:
                pauli = self.gen_dict[key][0]
                tr, coef = dict_cz_trans[pauli[num[0]] + pauli[num[1]]]
                pauli = list(pauli)
                pauli[num[0]] = tr[0]
                pauli[num[1]] = tr[1]
                self.gen_dict[key] = ["".join(pauli), c*coef*self.gen_dict[key][1]] 
        elif pauli == "h":
            for key in self.gen_dict:
                pauli = self.gen_dict[key][0]
                tr, coef = dict_h_trans[pauli[num[0]]]
                pauli = list(pauli)
                pauli[num[0]] = tr[0]
                self.gen_dict[key] = ["".join(pauli), c*coef*self.gen_dict[key][1]] 
        elif pauli == "x":
            for key in self.gen_dict:
                pauli = self.gen_dict[key][0]
                tr, coef = dict_x_trans[pauli[num[0]]]
                pauli = list(pauli)
                pauli[num[0]] = tr[0]
                self.gen_dict[key] = ["".join(pauli), c*coef*self.gen_dict[key][1]] 
        else:
            for key in self.gen_dict:
                ops, coef = prod(self.gen_dict[key][0], pauli)
                if coef.imag == 0:
                    self.gen_dict[key] = [ops, c*int(coef.real)*self.gen_dict[key][1]]


    def __len__(self):
        return len(self.gen_dict)

    def __str__(self):
        s = ''
        # for i in range(len(self) // 2):
        #     s += f"({2*i}, {2*i + 1}): {self.gen_dict[2*i]}, {self.gen_dict[2*i + 1]}\n"
        for key in self.gen_dict:
            if self.gen_dict[key][1] > 0:
                s += f"({key}): +{self.gen_dict[key][0]}\n"
            else:
                s += f"({key}): -{self.gen_dict[key][0]}\n"

        return s
    
    def __iter__(self):
        return self.gen_dict.__iter__()


def get_pauli(char='I', pos=0, numq=6):
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

from tqdm import tqdm


if __name__ == "__main__":
    n = 4
    # gen_dict = {0 : "YYZI",1: "YYIZ",2: "XXZI", 3: "XXIZ", 4: "ZIXX", 5: "IZXX", 6: "ZIYY",7: "IZYY"}
    gen_dict = {0 : "XXXY",1: "XXYX",2: "XYXX", 3: "YXXX", 4: "YYYX", 5: "YYXY", 6: "YXYY",7: "XYYY"}
    # gen_dict = {0 : "ZIZI",1: "ZIIZ",2: "IZZI", 3: "IZIZ", 4: "XXZI", 5: "XXIZ", 6: "YYZI", 7: "YYIZ"}
    # gen_dict = {0 : "XIII",1: "YIII",2: "ZXII", 3: "ZYII", 4: "ZZXI", 5: "ZZYI", 6: "ZZZX", 7: "ZZZY", 8: "ZZZZ"}
    # gen_dict = {0 : "YZXI",1: "XZYI",2: "IYZX", 3: "IXZY"}
    # gen_dict = {0 : "YZXI",1: "XZYI",2: "IYZX", 3: "IXZY"}
    # gen_dict = {0:"XX", 1: "YY"}
    
    # gen_dict = {0 : "YZZXII",1: "XZZYII",2: "IYZZXI", 3: "IXZZYI", 4: "IIXZZY", 5: "IIYZZX"}

    cg = clifford_generator(n, gen_dict)
    print(cg)
    # cg.transform_maj("z", (0,))
    # cg.transform_maj("z", (1,))
    # cg.transform_maj("cx", (0,1))
    # cg.transform_maj('h', (0,))
    # cg.transform_maj('h', (1,))
    # cg.transform_maj('h', (2,))
    # cg.transform_maj('h', (3,))
    # ---------------------------------------------
    cg.transform_maj("cx", (2,3))
    cg.transform_maj("cx", (0,1))
    cg.transform_maj('x', (3,))
    cg.transform_maj('h', (3,))
    cg.transform_maj('x', (1,))
    cg.transform_maj('h', (1,))
    cg.transform_maj("cx", (0,2))
    print(cg)
    cg.transform_maj("cx", (0,1))
    print(cg)
    cg.transform_maj("cx", (0,3))
    print(cg)
    cg.transform_maj("cx", (0,1))
    print(cg)
    cg.transform_maj('h', (2,))
    cg.transform_maj("cx", (0,2))
    print(cg)
    cg.transform_maj("cx", (0,1))
    print(cg)
    cg.transform_maj("cx", (0,3))
    print(cg)
    cg.transform_maj('h', (3,))
    cg.transform_maj("cx", (0,1))
    print(cg)
    cg.transform_maj("cx", (0,2))
    cg.transform_maj('h', (2,))

    cg.transform_maj("cx", (0,2))
    print(cg)
    cg.transform_maj('x', (3,))
    cg.transform_maj('h', (1,))
    cg.transform_maj('x', (1,))
    cg.transform_maj("cx", (2,3))
    cg.transform_maj("cx", (0,1))
    print(cg)
    # # ---------------------------------------------
    # # cg.transform_maj('h', (0,))
    # # cg.transform_maj('h', (1,))
    # # cg.transform_maj('h', (2,))
    # # cg.transform_maj('h', (3,))

    # cg.transform_maj(get_pauli("Y", 1))
    # # cg.transform_maj(get_pauli("Y", 2))
    # # cg.transform_maj(get_pauli("Y", 3))
    # cg.transform_maj("cx", (2, 1))
    # cg.transform_maj("cx", (0, 1))
    # cg.transform_maj("cx", (3, 2))
    # cg.transform_maj(get_pauli("X", 1))
    # cg.transform_maj(get_pauli("X", 2))
    # cg.transform_maj(get_pauli("Z", 0))
    # cg.transform_maj(get_pauli("Z", 3))


    # print(cg)
    # cg.transform_maj("cx", (1,3))
    # cg.transform_maj("cx", (0,2))

    # print(cg)
    # cg.transform_maj(get_pauli("Z", 1))
    # cg.transform_maj(get_pauli("Z", 2))
    # cg.transform_maj("cx", (1,2))

    # print(cg)
    # cg.transform_maj(get_pauli("X", 0))
    # # cg.transform_maj(get_pauli("Z", 1))
    # cg.transform_maj(get_pauli("Y", 2))    
    # cg.transform_maj(get_pauli("X", 3))
    # print(cg)
    # cg.transform_maj("cx", (1,3))
    # cg.transform_maj("cx", (2,0))
    # print(cg)
    # cg.transform_maj(get_pauli("Z", 0), coef=-1)
    # cg.transform_maj(get_pauli("Y", 1), coef=1)
    # cg.transform_maj(get_pauli("Y", 2), coef=1)
    # cg.transform_maj(get_pauli("Z", 3), coef=-1)

    # cg.transform_maj("cx", (1,3))
    # cg.transform_maj("cx", (2,0))
    # cg.transform_maj(get_pauli("X", 0), coef=1)
    # cg.transform_maj(get_pauli("Y", 2), coef=-1)
    # cg.transform_maj(get_pauli("X", 3), coef=-1)
    # cg.transform_maj("cx", (2,1))
    # cg.transform_maj(get_pauli("Y", 1), coef=1)
    # print(cg)


    # cg.transform_maj(get_pauli("Y", 1))
    # # # cg.transform_maj(get_pauli("Y", 2))
    # # # cg.transform_maj(get_pauli("Y", 3))
    # cg.transform_maj("cx", (2,1))
    # # print(cg)
    # cg.transform_maj(get_pauli("X", 0))
    # # cg.transform_maj(get_pauli("Z", 1))
    # cg.transform_maj(get_pauli("Y", 2))    
    # cg.transform_maj(get_pauli("X", 3))
    # print(cg)
    # cg.transform_maj("cx", (1,3))
    # cg.transform_maj("cx", (2,0))
    # print(cg)
    # cg.transform_maj(get_pauli("Z", 0), coef=-1)
    # cg.transform_maj(get_pauli("Y", 1), coef=1)
    # cg.transform_maj(get_pauli("Y", 2), coef=1)
    # cg.transform_maj(get_pauli("Z", 3), coef=-1)

    # cg.transform_maj("cx", (1,3))
    # cg.transform_maj("cx", (2,0))
    # cg.transform_maj(get_pauli("X", 0), coef=1)
    # cg.transform_maj(get_pauli("Y", 2), coef=-1)
    # cg.transform_maj(get_pauli("X", 3), coef=-1)
    # cg.transform_maj("cx", (2,1))
    # cg.transform_maj(get_pauli("Y", 1), coef=1)
    # print(cg)

    # cg.transform_maj(get_pauli("Y", 3))

    # gen_dict = {0 : "XY",1: "YX"}
    # --------------------------------
    # cg = clifford_generator(n)
    # print(cg)
    # cg.transform_maj(get_pauli("X", 0))    
    # cg.transform_maj(get_pauli("Y", 1))
    # cg.transform_maj("cx", (1 ,0))
    # print(cg)
    # cg.transform_maj(get_pauli("Y", 1))    
    # cg.transform_maj(get_pauli("Z", 0), coef=-1)    
    # cg.transform_maj("cx", (1,0))
    # cg.transform_maj(get_pauli("Y", 1), coef=1)
    # cg.transform_maj(get_pauli("X", 0), coef=1)    
    # ------------------------------
    # print(cg)
    # cg.transform_maj(get_pauli("X", 0), coef=-1)    
    # cg.transform_maj(get_pauli("X", 1), coef=-1)
    # cg.transform_maj("cz", (1,0))
    # cg.transform_maj(get_pauli("X", 1), coef=-1)    
    # cg.transform_maj(get_pauli("X", 0), coef=-1)    
    # cg.transform_maj("cz", (1,0))
    # cg.transform_maj(get_pauli("X", 1), coef=-1)    
    # cg.transform_maj(get_pauli("X", 0), coef=-1)    
    # cg.transform_maj(get_pauli("Z", 1), coef=-1)
    # cg.transform_maj(get_pauli("Z", 0), coef=-1)    
        
    
    
    # cg.transform_maj(get_pauli("X", 0), coef=-1)    
    # cg.transform_maj(get_pauli("Y", 1), coef=-1)
    # cg.transform_maj("cx", (1 ,0))
    # print(cg)
    # cg.transform_maj(get_pauli("Y", 1), coef=-1)    
    # cg.transform_maj(get_pauli("Z", 0), coef=1)    
    # cg.transform_maj("cx", (1,0))
    # cg.transform_maj(get_pauli("Y", 1), coef=-1)
    # cg.transform_maj(get_pauli("X", 0), coef=-1)   
    # print(cg)











    
    # cg.transform_maj(get_pauli("X",0))
    # cg.transform_maj(get_pauli("X",3))
    # # print(cg)
    # cg.transform_maj("cz", (0,1))
    # cg.transform_maj("cz", (2,3))
    # print(cg)
    # cg.transform_maj(get_pauli("X",1))
    # cg.transform_maj(get_pauli("X",2))
    # # print(cg)
    # cg.transform_maj("cz", (1,2))
    # # cg.transform_maj("cz", (2,3))
    # print(cg)
    # cg.transform_maj(get_pauli("X",0))
    # cg.transform_maj(get_pauli("X",3))
    # # print(cg)
    # cg.transform_maj("cz", (0,1))
    # cg.transform_maj("cz", (2,3))
    # print(cg)
    # cg.transform_maj("cz", (0,2))
    # cg.transform_maj("cz", (1,3))
    # print(cg)
    # # cg.transform_maj("cz", (0,1))
    # # cg.transform_maj("cz", (2,3))
    # # print(cg)
    # # cg.transform_maj("cz", (1,2))
    # # cg.transform_maj(get_pauli("X",0))
    # # cg.transform_maj(get_pauli("X",3))
    # # cg.transform_maj(get_pauli("X",1))
    # # cg.transform_maj(get_pauli("X",2))
    # # print(cg)
    # # cg.transform_maj("cz", (0,1))
    # # cg.transform_maj("cz", (2,3))
    # # print(cg)

