from itertools import combinations
import sys


CX_TRANSFORM = {
    "II": ("II", 1),
    "IX": ("IX", 1),
    "IY": ("ZY", 1),
    "IZ": ("ZZ", 1),
    "XI": ("XX", 1),
    "YI": ("YX", 1),
    "ZI": ("ZI", 1),
    "XX": ("XI", 1),
    "ZZ": ("IZ", 1),
    "ZY": ("IY", 1),
    "YX": ("YI", 1),
    "XY": ("YZ", 1),
    "YZ": ("XY", 1),
    "XZ": ("YY", -1),
    "YY": ("XZ", -1),
    "ZX": ("ZX", 1),
}

CZ_TRANSFORM = {
    "II": ("II", 1),
    "IZ": ("IZ", 1),
    "ZI": ("ZI", 1),
    "ZZ": ("ZZ", 1),
    "XX": ("YY", -1),
    "YY": ("XX", -1),
    "XY": ("YX", 1),
    "YX": ("XY", 1),
    "XZ": ("XI", 1),
    "XI": ("XZ", 1),
    "ZX": ("IX", 1),
    "IX": ("ZX", 1),
    "YI": ("YZ", 1),
    "YZ": ("YI", 1),
    "ZY": ("IY", 1),
    "IY": ("ZY", 1),
}

H_TRANSFORM = {
    "I": ("I", 1),
    "X": ("Z", 1),
    "Y": ("Y", -1),
    "Z": ("X", 1),
}

X_TRANSFORM = {
    "I": ("I", 1),
    "X": ("X", 1),
    "Y": ("Y", -1),
    "Z": ("Z", -1),
}

PRODUCT = {
    "II": ("I", 1),
    "XX": ("I", 1),
    "YY": ("I", 1),
    "ZZ": ("I", 1),
    "XY": ("Z", 1j),
    "YX": ("Z", -1j),
    "XZ": ("Y", -1j),
    "ZX": ("Y", 1j),
    "YZ": ("X", 1j),
    "ZY": ("X", -1j),
    "IX": ("X", 1),
    "XI": ("X", 1),
    "YI": ("Y", 1),
    "IY": ("Y", 1),
    "IZ": ("Z", 1),
    "ZI": ("Z", 1),
}

DEFAULT_GENERATORS = {
    0: ("YYIZ", -1),
    1: ("XXIZ", -1),
    2: ("YYZI", 1),
    3: ("XXZI", 1),
    4: ("IZYY", 1),
    5: ("IZXX", 1),
    6: ("ZIYY", -1),
    7: ("ZIXX", -1),
}


def prod(pauli1: str, pauli2: str):
    ops = []
    coeff = 1
    for op1, op2 in zip(pauli1, pauli2):
        op, coef = PRODUCT[op1 + op2]
        ops.append(op)
        coeff *= coef
    return "".join(ops), coeff


def weight(pauli: str):
    return sum(op != "I" for op in pauli)


def jw_majoranas(n: int):
    gen_dict = {}
    for index in range(n):
        prefix = "Z" * index
        suffix = "I" * (n - index - 1)
        gen_dict[2 * index] = [prefix + "X" + suffix, 1]
        gen_dict[2 * index + 1] = [prefix + "Y" + suffix, 1]
    return gen_dict


class CliffordGenerator:
    def __init__(self, n=1, gen_dict=None):
        self.gen_dict = jw_majoranas(n) if gen_dict is None else gen_dict

    @property
    def gen_dict(self):
        return self._gen_dict

    @gen_dict.setter
    def gen_dict(self, values):
        self._gen_dict = {}
        for key, value in values.items():
            if isinstance(value, (list, tuple)):
                self._gen_dict[key] = [value[0], value[1]]
            else:
                self._gen_dict[key] = [value, 1]

    def get_jw_majorana(self, n):
        return jw_majoranas(n)

    def copy(self):
        return CliffordGenerator(gen_dict={key: value[:] for key, value in self.gen_dict.items()})

    def transform_maj(self, pauli: str, num=None, coef=1):
        if pauli == "cx":
            self._transform_two_qubit(num, coef, CX_TRANSFORM)
        elif pauli == "cz":
            self._transform_two_qubit(num, coef, CZ_TRANSFORM)
        elif pauli == "x":
            self._transform_one_qubit(num, coef, X_TRANSFORM)
        elif pauli == "h":
            self._transform_one_qubit(num, coef, H_TRANSFORM)
        else:
            self._multiply_by_pauli(pauli, coef)

    def _transform_two_qubit(self, qubits, coef, transform):
        q0, q1 = qubits
        for key, value in self.gen_dict.items():
            pauli, sign = value
            transformed, transform_sign = transform[pauli[q0] + pauli[q1]]
            pauli = list(pauli)
            pauli[q0] = transformed[0]
            pauli[q1] = transformed[1]
            self.gen_dict[key] = ["".join(pauli), coef * transform_sign * sign]

    def _transform_one_qubit(self, qubits, coef, transform):
        qubit = qubits[0]
        for key, value in self.gen_dict.items():
            pauli, sign = value
            transformed, transform_sign = transform[pauli[qubit]]
            pauli = list(pauli)
            pauli[qubit] = transformed
            self.gen_dict[key] = ["".join(pauli), coef * transform_sign * sign]

    def _multiply_by_pauli(self, pauli, coef):
        for key, value in self.gen_dict.items():
            ops, product_coef = prod(value[0], pauli)
            if product_coef.real == 0:
                self.gen_dict[key] = [ops, coef * int(-product_coef.imag) * value[1]]

    def __len__(self):
        return len(self.gen_dict)

    def __iter__(self):
        return iter(self.gen_dict)

    def __str__(self):
        rows = []
        for key, value in self.gen_dict.items():
            sign = "+" if value[1] > 0 else "-"
            rows.append(f"({key}): {sign}{value[0]}")
        return "\n".join(rows) + "\n"


def get_pauli(char="I", pos=0, numq=4):
    pauli = ["I" for _ in range(numq)]
    pauli[pos] = char
    return "".join(pauli)


def check_anticom(pauli1, pauli2):
    return prod(pauli1, pauli2)[1] != prod(pauli2, pauli1)[1]


def full_clifford_basis(generator: CliffordGenerator):
    basis = {}
    maj_tuple = tuple(range(len(generator)))
    for size in range(1, len(maj_tuple) + 1):
        for elements in combinations(maj_tuple, size):
            op = "I" * (len(generator) // 2)
            for element in elements:
                op = prod(op, generator.gen_dict[element][0])[0]
            basis[elements] = op
    return basis


def count_anticommuting_bases(n):
    generator = CliffordGenerator(n)
    paulis = full_clifford_basis(generator)
    bases = []
    for item in combinations(paulis, 2 * n):
        if all(check_anticom(paulis[left], paulis[right]) for left, right in combinations(item, 2)):
            bases.append(item)
    return bases


def num_basis(n):
    bases = count_anticommuting_bases(n)
    for basis in bases:
        print(basis)
    print(len(bases))


def num_anticom(n):
    generator = CliffordGenerator(n)
    paulis = full_clifford_basis(generator)
    pauli_0 = (0, 1, 2, 3)
    for pauli_1 in paulis:
        for pauli_4 in paulis:
            counts2 = 0
            for pauli_3 in paulis:
                counts = 0
                for pauli_2 in paulis:
                    item = [pauli_0, pauli_1, pauli_2, pauli_3, pauli_4]
                    if all(check_anticom(paulis[left], paulis[right]) for left, right in combinations(item, 2)):
                        counts += 1
                if counts > 0:
                    counts2 += 1
                    print(counts)
            if counts2 > 0:
                print(counts2)


def default_generator():
    return CliffordGenerator(4, DEFAULT_GENERATORS)


def apply_pauli(generator, char, pos, coef=1, numq=4):
    generator.transform_maj(get_pauli(char, pos, numq), coef=coef)


def decompose_short(generator=None):
    generator = default_generator() if generator is None else generator
    snapshots = [str(generator)]
    generator.transform_maj("cz", (1, 2))
    snapshots.append(str(generator))
    generator.transform_maj("cz", (1, 2))
    generator.transform_maj("cz", (0, 3))
    snapshots.append(str(generator))
    return generator, snapshots


def decompose_linear_short(generator=None):
    generator = default_generator() if generator is None else generator
    generator.transform_maj("cx", (1, 0))
    generator.transform_maj("cx", (2, 3))
    apply_pauli(generator, "Y", 2)
    generator.transform_maj("cx", (2, 1))
    generator.transform_maj("cx", (0, 1))
    generator.transform_maj("cx", (3, 2))
    return generator


def decompose_zyx_12cnot(generator=None):
    generator = default_generator() if generator is None else generator
    apply_pauli(generator, "Z", 0)
    apply_pauli(generator, "Z", 3)
    apply_pauli(generator, "Y", 1)
    apply_pauli(generator, "Y", 2)
    generator.transform_maj("cz", (1, 2))
    apply_pauli(generator, "X", 1)
    apply_pauli(generator, "X", 2)
    generator.transform_maj("cz", (0, 1))
    generator.transform_maj("cz", (2, 3))
    generator.transform_maj("cz", (2, 3))
    generator.transform_maj("cz", (0, 1))
    apply_pauli(generator, "X", 0)
    apply_pauli(generator, "X", 1, coef=-1)
    apply_pauli(generator, "X", 2, coef=-1)
    apply_pauli(generator, "X", 3)
    generator.transform_maj("cz", (1, 2))
    generator.transform_maj("cz", (0, 3))
    apply_pauli(generator, "X", 0)
    apply_pauli(generator, "X", 1)
    apply_pauli(generator, "X", 2)
    apply_pauli(generator, "X", 3)
    generator.transform_maj("cz", (2, 3))
    generator.transform_maj("cz", (0, 1))
    generator.transform_maj("cz", (0, 1))
    generator.transform_maj("cz", (2, 3))
    apply_pauli(generator, "X", 3, coef=-1)
    apply_pauli(generator, "X", 0, coef=-1)
    generator.transform_maj("cz", (0, 3))
    apply_pauli(generator, "X", 0, coef=-1)
    apply_pauli(generator, "X", 3, coef=-1)
    apply_pauli(generator, "X", 1, coef=-1)
    apply_pauli(generator, "X", 2, coef=-1)
    apply_pauli(generator, "Z", 0, coef=-1)
    apply_pauli(generator, "Z", 3, coef=-1)
    apply_pauli(generator, "Y", 1, coef=-1)
    apply_pauli(generator, "Y", 2, coef=-1)
    return generator


def decompose_12cnot(generator=None):
    generator = default_generator() if generator is None else generator
    apply_pauli(generator, "Z", 0)
    apply_pauli(generator, "Z", 3)
    apply_pauli(generator, "Y", 0)
    apply_pauli(generator, "Y", 2)
    generator.transform_maj("cx", (1, 2))
    generator.transform_maj("cx", (0, 1))
    generator.transform_maj("cx", (2, 3))
    generator.transform_maj("cx", (2, 3))
    generator.transform_maj("cx", (0, 1))
    generator.transform_maj("cx", (1, 2))
    generator.transform_maj("cx", (3, 0))
    generator.transform_maj("cx", (2, 3))
    generator.transform_maj("cx", (0, 1))
    return generator


def decompose_yordan(generator=None):
    generator = default_generator() if generator is None else generator
    generator.transform_maj("cx", (0, 1))
    generator.transform_maj("cx", (2, 3))
    generator.transform_maj("x", (1,))
    generator.transform_maj("h", (1,))
    generator.transform_maj("x", (3,))
    generator.transform_maj("h", (3,))
    generator.transform_maj("cx", (0, 2))
    generator.transform_maj("h", (2,))
    generator.transform_maj("cx", (0, 1))
    generator.transform_maj("cx", (0, 3))
    generator.transform_maj("cx", (0, 1))
    generator.transform_maj("cx", (0, 2))
    generator.transform_maj("cx", (0, 1))
    generator.transform_maj("cx", (0, 3))
    generator.transform_maj("cx", (0, 1))
    generator.transform_maj("h", (1,))
    generator.transform_maj("x", (1,))
    generator.transform_maj("h", (3,))
    generator.transform_maj("x", (3,))
    apply_pauli(generator, "Z", 2)
    generator.transform_maj("cx", (0, 2))
    apply_pauli(generator, "Z", 0, coef=-1)
    apply_pauli(generator, "Z", 2)
    apply_pauli(generator, "Y", 2)
    generator.transform_maj("cx", (0, 1))
    generator.transform_maj("cx", (2, 3))
    return generator


DECOMPOSITIONS = {
    "short": decompose_short,
    "linear-short": decompose_linear_short,
    "zyx-12cnot": decompose_zyx_12cnot,
    "12cnot": decompose_12cnot,
    "yordan": decompose_yordan,
}


def run_decomposition(name="short"):
    if name not in DECOMPOSITIONS:
        choices = ", ".join(sorted(DECOMPOSITIONS))
        raise ValueError(f"Unknown decomposition {name!r}. Choices: {choices}")

    result = DECOMPOSITIONS[name]()
    if name == "short":
        _, snapshots = result
        for snapshot in snapshots:
            print(snapshot)
    else:
        print(result)


if __name__ == "__main__":
    run_decomposition(sys.argv[1] if len(sys.argv) > 1 else "short")
