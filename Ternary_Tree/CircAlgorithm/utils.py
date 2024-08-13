import itertools
from .paulistring import PauliStringSign

def single_exc(op) -> tuple[tuple[int, int], int]:
    if op[1] != op[0]:
        if len(set(op)) != 1:
            yield tuple([2 * op[0] + 1, 2 * op[1] + 1]), 1
            yield tuple([2 * op[0] + 2, 2 * op[1] + 2]), 1


def double_exc(op) -> tuple[tuple[int, int, int, int], int]:
    dcoef = [[1, 1j], [1, 1j], [1, -1j], [1, -1j]]
    for elems in itertools.product([1, 2], repeat=4):
        coef = 1j
        for i in range(4):
            coef *= (dcoef[i][elems[i] - 1])
        if coef.imag == 0:
            j = [2 * op[i] + elems[i] for i in range(4)]
            if any(op.count(x) % 2 for x in set(op)):
                yield tuple(j), coef.real

def _parity(t: tuple):
    if len(set(t)) == 1:
        return 0
    par = 1
    for index, el1 in enumerate(t):
        for el2 in t[index + 1:]:
            if el1 > el2:
                par = par * (-1)
    return par

def _arranging(t: tuple):
    new_t = tuple(sorted(t))
    if len(new_t) == 4:
        new_t = new_t[:2] if new_t[2] == new_t[3] else new_t
    new_t = new_t[2:] if new_t[1] == new_t[0] else new_t
    return new_t

def simplification(maj_exc):
    maj = {}
    for key in maj_exc:
        new_key, coef = _arranging(key), _parity(key)
        maj[new_key] = maj.get(new_key, 0) + coef * maj_exc[key]
    return maj


def _get_double_excitations(self):
    """
    Method to get spin preserving double excitation operators via ladder
    operators
    """
    alpha_occ = range(self.num_alpha)
    alpha_unocc = range(self.num_alpha, self.num_spatial_orbitals)

    beta_index_offset = self.num_spatial_orbitals
    beta_occ = range(beta_index_offset, beta_index_offset + self.num_beta)
    beta_unocc = range(beta_index_offset + self.num_beta, self.num_spin_orbitals)

    aa = list(itertools.product(alpha_unocc, alpha_unocc, alpha_occ, alpha_occ))
    bb = list(itertools.product(beta_unocc, beta_unocc, beta_occ, beta_occ))
    ab = list(itertools.product(beta_unocc, alpha_unocc, beta_occ, alpha_occ))

    return aa + bb + ab


