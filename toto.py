from math import pi
import itertools
import numpy as np
from scipy.linalg import expm
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Operator, Pauli
from qiskit.quantum_info.operators import PauliList

# Словарь знаков из get_pauli_double_ex_short()
# signs = {
#     "XXXY":  1, "XXYX":  1, "XYXX": -1, "YXXX": -1,
#     "YYYX": -1, "YYXY": -1, "YXYY":  1, "XYYY":  1
# }
short = {'XXXY': -1.0, 'XXYX': 1.0, 'XYXX': 1.0, 'YXXX': 1.0, 'YYYX': 1.0, 'YYXY': -1.0, 'YXYY': -1.0, 'XYYY': -1.0}

yordan = {'XXXY': 1.0, 'XXYX': -1.0, 'XYXX': -1.0, 'YXXX': -1.0, 'YYYX': -1.0, 'YYXY': 1.0, 'YXYY': 1.0, 'XYYY': 1.0}

def yordan_signs():
    return {"XXXY": -1, "XXYX": 1, "XYXX": -1, "YXXX": 1, "YYYX": 1, "YYXY": -1, "YXYY": 1, "XYYY": -1}

def double_ex_yordan(circ,
                    q0, q1, q2, q3,
                    parameter,
                    list_signs):
    circ.cx(q0, q1), circ.cx(q2, q3)
    circ.x(q1), circ.x(q3)
    circ.cx(q0, q2)
    circ.h(q1), circ.ry(list_signs["YXXX"]*parameter*0.25, q0)
    circ.cx(q0, q1)
    circ.h(q3), circ.ry(list_signs["XYXX"]*parameter*0.25, q0)
    circ.cx(q0, q3)
    circ.ry(list_signs["XYYY"]*parameter*0.25, q0)
    circ.cx(q0, q1)
    circ.h(q2), circ.ry(list_signs["YXYY"]*parameter*0.25, q0)
    circ.cx(q0, q2)
    circ.ry(list_signs["XXXY"]*parameter*0.25, q0)
    circ.cx(q0, q1)
    circ.ry(list_signs["YYXY"]*parameter*0.25, q0)
    circ.cx(q0, q3)
    circ.h(q3), circ.ry(list_signs["YYYX"]*parameter*0.25, q0)
    circ.cx(q0, q1)
    circ.h(q1), circ.ry(list_signs["XXYX"]*parameter*0.25, q0)
    circ.rz(pi/2, q2)
    circ.cx(q0, q2)
    circ.rz(pi/2, q2), circ.rz(-pi/2, q0)
    circ.ry(pi/2, q2)
    circ.x(q1), circ.x(q3)
    circ.cx(q0, q1), circ.cx(q2, q3)

from math import pi
import numpy as np
from scipy.linalg import expm
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator, Pauli
from qiskit.quantum_info.operators import PauliList

signs = {'XXXY': 1, 'XXYX': 1, 'XYXX': -1, 'YXXX': -1, 'YYYX': -1, 'YYXY': -1, 'YXYY': 1, 'XYYY': 1}

def double_ex_short(circ, q0, q1, q2, q3, parameter, list_signs):        
    # list_signs = {''.join(reversed(l)): list_signs[l] for l in list_signs}
    # print(list_signs)
    circ.cx(q1, q0)
    circ.cx(q2, q3)
    circ.ry(pi/2, q2)
    circ.cx(q2, q1)
    circ.ry(parameter*list_signs["XYXX"]*0.25, q1)
    circ.ry( parameter*list_signs["XXYX"]*0.25, q2)
 
    circ.cx(q0, q1)
    circ.cx(q3, q2)
    circ.ry( parameter*list_signs["XXXY"]*0.25, q2)
    circ.ry(parameter*list_signs["YXXX"]*0.25, q1)
    
    circ.cx(q0, q2)
    circ.cx(q3, q1)
    circ.ry(parameter*list_signs["YYXY"]*0.25, q2)
    circ.ry( parameter*list_signs["YXYY"]*0.25, q1)
    
    circ.cx(q0, q1)
    circ.cx(q3, q2)
    circ.ry( parameter*list_signs["XYYY"]*0.25, q1)
    circ.ry(parameter*list_signs["YYYX"]*0.25, q2)
    circ.cx(q0, q2)
    circ.cx(q3, q1)
    circ.cx(q2, q1)
    circ.ry(-pi/2, q2)
    circ.cx(q1, q0)
    circ.cx(q2, q3)

def get_clif_circuit():
    circ = QuantumCircuit(4)
    circ.cx(3,2)
    circ.cx(0,1)
    circ.cx(1,2)
    circ.ry(pi/2, 1)
    circ.cx(2,3)
    circ.cx(1,0)
    # circ.rz(pi/2, 0)
    # circ.rz(pi/2, 1)
    # circ.rz(pi/2, 2)
    # circ.rz(pi/2, 3)
    # circ.cx(2,3)
    # circ.cx(1,0)
    # circ.ry(pi/2, 1)
    # circ.cx(1,2)
    # circ.cx(3,2)
    # circ.cx(0,1)
    return circ 


def build_short_circuit(theta):
    qc = QuantumCircuit(4)
    qc.cx
    # double_ex_short(qc, 3,2,1,0, theta, signs)
    double_ex_short(qc, 0,1,2,3, theta, signs)
    
    return qc
def build_yordan_circuit(theta):
    qc = QuantumCircuit(4)
    double_ex_yordan(qc, 0,1,2,3, theta, yordan_signs())
    return qc


def hamiltonian_matrix(signs_dict):
    # Построим H = sum s_alpha * P_alpha
    # each P_alpha — tensor of 4 Паули
    pauli_ops = {
        'I': np.array([[1,0],[0,1]], dtype=complex),
        'X': np.array([[0,1],[1,0]], dtype=complex),
        'Y': np.array([[0,-1j],[1j,0]], dtype=complex),
        'Z': np.array([[1,0],[0,-1]], dtype=complex),
    }
    H = np.zeros((16,16), dtype=complex)
    for pauli_str, s in signs_dict.items():
        # tensor product of single-qubit matrices
        op = pauli_ops[pauli_str[0]]
        for j in pauli_str[1:]:
            op = np.kron( pauli_ops[j], op)
        H += op
    return H

def unitary_from_exponential(H, theta):
    # exp(-i * theta/8 * H)
    return expm(-1j * theta/8 * H)

def compare_unitaries(U1, U2, tol=1e-8):
    # Найдём фазовый множитель ϕ = <U1,U2> / dim
    dim = U1.shape[0]
    phi = np.vdot(U1.flatten(), U2.flatten())/np.abs(np.vdot(U1.flatten(), U2.flatten()))
    # Проверим норму разности после учёта фазы
    delta = np.linalg.norm(U1 - phi.conj() * U2)
    return delta, phi

if __name__ == "__main__":
    theta = 0.2  # произвольный параметр
    # 1) Схема → унитарная матрица
    # for combo in itertools.product([-1,1], repeat=8):
    #     signs = {key: combo[i] for i, key in enumerate(signs.keys())}
    #     qc = build_short_circuit(theta)
    #     qcy = build_yordan_circuit(theta)
    #     U_circ = Operator(qc).data
    #     U_circy = Operator(qcy).data
    #     # 2) Гамильтониан → expm
    #     H = hamiltonian_matrix(signs)
    #     U_exp = unitary_from_exponential(H, theta)

    #     # 3) Сравнение
    #     delta, phi = compare_unitaries(U_exp,  U_circ)
    #     # print(f"||U_circ - ϕ·U_exp|| = {delta:.3e}")
    #     # print(f"Глобальная фаза ϕ = {phi}")
    #     if delta < 1e-6:
    #         print(signs)
    #         print(qc)
    #         print("✔ Схема совпадает с экспонентой Гамильтониана (с учётом глобальной фазы).")
    #         # break
    #     else:
    #         pass
    #         # print("✘ Есть расхождение!")
    circ = get_clif_circuit()
    print(circ)
    print(transpile(circ,optimization_level=3).decompose())
    # tra