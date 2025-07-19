from math import pi
import numpy as np
from scipy.linalg import expm
import itertools

# Базовые матрицы
I = np.array([[1,0],[0,1]], dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)

# Однокубитная RY и CNOT
def RY(theta):
    return np.cos(theta/2)*I - 1j*np.sin(theta/2)*Y

def CNOT(control, target, num_qubits=4):
    # Строим полный 2^n×2^n оператор
    P0 = np.array([[1,0],[0,0]], dtype=complex)
    P1 = np.array([[0,0],[0,1]], dtype=complex)
    ops = []
    for i in range(num_qubits):
        if i == control:
            ops.append(P0)
        else:
            ops.append(I)
    term0 = ops[0]
    for op in ops[1:]:
        term0 = np.kron(term0, op)
    ops = []
    for i in range(num_qubits):
        if i == control:
            ops.append(P1)
        elif i == target:
            ops.append(X)
        else:
            ops.append(I)
    term1 = ops[0]
    for op in ops[1:]:
        term1 = np.kron(term1, op)
    return term0 + term1

# Перемножение цепочки gate_list
def apply_gates(gate_list):
    U = np.eye(2**4, dtype=complex)
    for G in gate_list:
        U = G @ U
    return U

# Собираем цепочку double_ex_short в матричном виде
def U_double_ex_short(theta, signs):
    gates = []
    # Helper to add RY on qubit j
    def ry_q(theta, j):
        mat = 1
        for k in range(4):
            mat = np.kron(mat, RY(theta)) if k == j else np.kron(mat, I)
        return mat

    # Список операций по коду
    ops = [
        ('cx', 1,0), ('cx', 2,3), ('ry', pi/2, 2), ('cx', 2,1),
        ('ry', theta*signs["XXXY"]*0.25, 2),
        ('ry', theta*signs["YXXX"]*0.25, 1),
        ('cx', 0,1), ('cx', 3,2),
        ('ry', theta*signs["YXYY"]*0.25, 1),
        ('ry', theta*signs["YYXY"]*0.25, 2),
        ('cx', 0,2), ('cx', 3,1),
        ('ry', theta*signs["XYYY"]*0.25, 1),
        ('ry', theta*signs["YYYX"]*0.25, 2),
        ('cx', 0,1), ('cx', 3,2),
        ('ry', theta*signs["XYXX"]*0.25, 1),
        ('ry', theta*signs["XXYX"]*0.25, 2),
        ('cx', 0,2), ('cx', 3,1),
        ('cx', 2,1), ('ry', -pi/2, 2),
        ('cx', 1,0), ('cx', 2,3),
    ]
    for op in ops:
        if op[0] == 'cx':
            gates.append(CNOT(op[1], op[2]))
        elif op[0] == 'ry':
            gates.append(ry_q(op[1], op[2]))
    return apply_gates(gates)

# Гамильтониан
pauli_keys = ["XXXY","XXYX","XYXX","YXXX","YYYX","YYXY","YXYY","XYYY"]
def H_matrix(signs):
    H = np.zeros((16,16), dtype=complex)
    for key, s in signs.items():
        mats = {'X':X,'Y':Y,'Z':Z}
        op = mats[key[0]]
        for ch in key[1:]:
            op = np.kron(op, mats[ch])
        H += s*op
    return H

# Бруткфорс поиска
theta = 0.223
for combo in itertools.product([-1,1], repeat=8):
    signs = dict(zip(pauli_keys, combo))
    U_circ = U_double_ex_short(theta, signs)
    U_exp = expm(-1j*theta/8 * H_matrix(signs))
    # глобальная фаза
    phi = np.vdot(U_circ.flatten(), U_exp.flatten())
    phi /= abs(phi)
    delta = np.linalg.norm(U_circ - phi*U_exp)
    if delta < 1e-6:
        print("Найден подходящий набор знаков:")
        print(signs)
        break
else:
    print("Подходящий набор не найден.")