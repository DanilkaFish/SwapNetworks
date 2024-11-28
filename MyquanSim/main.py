import numpy as np
from ideal_qc import QuantumCircuit, get_probs
import matplotlib.pyplot as plt
from itertools import product

def psi_Z_init(num_qubits):
    psi_0 = np.array([1, 0])
    psi_init = np.copy(psi_0)
    for _ in range(num_qubits - 1):
        psi_init = np.tensordot(psi_init, psi_0, axes=0)
    return psi_init

def psi_X_init(num_qubits):
    psi_0 = np.array([1, -1])/np.sqrt(2)
    psi_init = np.copy(psi_0)
    for _ in range(num_qubits - 1):
        psi_init = np.tensordot(psi_init, psi_0, axes=0)
    return psi_init

def get_cumulative_mean(x: np.ndarray) -> np.ndarray:
    """Get cumulative mean.

    :param x: Tensor, axis0 - n repeats, axis1 n experiments.
    :return: Cumulative mean.
    """
    cm_list = np.cumsum(x, axis=0)
    n = np.arange(1, x.shape[0] + 1)
    return np.einsum("r...,r->r...", cm_list, 1 / n)

def get_qc_p_th(qc: QuantumCircuit, psi_init, meas):
    instructions = qc.instructions
    thetas_sigmas = []
    for instr in instructions:
        if instr.name in {"rxc", "rxd"}:
            for _ in instr.qubits:
                thetas_sigmas.append((instr.theta0, instr.sigma2))

    probs = np.zeros(2**qc.num_qubits)

    for el in product([0, 1], repeat=len(thetas_sigmas)):
        qci = QuantumCircuit(qc.num_qubits)
        l = 0
        pgates = 1
        for instr in instructions:
            if instr.name in {"rxc", "rxd"}:
                for qubit in instr.qubits:
                    qci.rx([qubit], thetas_sigmas[l][0] + el[l]*np.pi)
                    pgates *= (1 + (-1)**el[l] * np.exp( - thetas_sigmas[l][1]/2)) / 2
                    l += 1
            else:
                qci.add_instr(instr)
            
        psi = qci.run(psi_init, counts=1)
        probs += get_probs(psi, meas)[0][0] * pgates
    return probs


if __name__ == "__main__":
    N = 5
    K = 2
    counts = 100
    n_exp = 10
    theta = np.pi/5
    sigma = np.pi/5

    qubs = [i for i in range(N)]
    qcd = QuantumCircuit(N)
    qcc = QuantumCircuit(N)
    for _ in range(K):
        qcd.rxd(qubs, theta, sigma)
        qcc.rxc(qubs, theta, sigma)
        for i in range(N - 2):
            qcd.cz(i, i + 1)
            qcc.cz(i, i + 1)
    measZ = np.array([[1, 0],[0, 0]])
    
    psi_init = psi_Z_init(N)
    ccumprobs = np.empty([n_exp, counts, 2**N])
    dcumprobs = np.empty([n_exp, counts, 2**N])

    for i in range(n_exp):
        psi_array = qcd.run(psi_init, counts=counts)
        probs, res = get_probs(psi_array, np.full((N, *measZ.shape), measZ))
        ccumprobs[i] = get_cumulative_mean(probs)

        psi_array = qcc.run(psi_init, counts=counts)
        probs, res = get_probs(psi_array, np.full((N, *measZ.shape), measZ))
        dcumprobs[i] = get_cumulative_mean(probs)

    p_theory = get_qc_p_th(qcc, psi_init, np.full((N, *measZ.shape), measZ))

    f, axs = plt.subplots(1, 2)
    index = 0
    border = 0.05
    axs[0].plot(ccumprobs.transpose([2, 1, 0])[index])
    axs[0].hlines(p_theory[index], 0, counts, color="black", linestyle="dashed")
    axs[0].set_ylim(p_theory[index] - border, p_theory[index] + border)

    axs[1].plot(dcumprobs.transpose([2, 1, 0])[index])
    axs[1].hlines(p_theory[index], 0, counts, color="black", linestyle="dashed")
    axs[1].set_ylim(p_theory[index] - border, p_theory[index] + border)

    plt.savefig("try.png", dpi=400)