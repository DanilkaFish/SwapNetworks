from qiskit_algorithms.optimizers import SPSA ,CG, SLSQP, L_BFGS_B ,COBYLA
from qiskit import QuantumCircuit
from qiskit.quantum_info.operators import SparsePauliOp
from qiskit.circuit.library.standard_gates import IGate, XGate, ZGate, YGate

from Ternary_Tree.optimizer.soap import SOAP

from utils import *
from MyquanSim.noise_qc import NoisyCircuit

# molecules = [('H 0 0 0; H 0 0 0.7349', (1,1), [0,1], '6-31g'), ('H 0 0 0; H 0 0 0.7349', (1,1), [0,1,2,3], '6-31g')]
# molecules = [('H 0 0 0; H 0 0 0.7349', (1,1), [0,1,2,3 ], '6-31g')]
molecules = [('H 0 0 0; H 0 0 0.7349', (1,1), [0,1 ], '6-31g')]
reps = [1]

n = 8
Noise_ops = ["D", "Z", "X", "Y"]
# Noise_ops = ["X"]
# Noise_ops = ['P','A']
def qiskit_to_mycirc(qiskit_circ, ls) -> NoisyCircuit:
    n_qubits = qiskit_circ.num_qubits
    my_circ = NoisyCircuit(n_qubits)
    counter_cx = 0
    counter_s = 0

    for inst in qiskit_circ:
        if inst.operation.name in {'u3', 'u'}:
            # if counter_s in ls:
                # params = [-i for i in inst.operation.params]

            my_circ.u3([n_qubits - 1 - inst.qubits[0].index], inst.operation.params)
            # else:
                # my_circ.u3_witout_noise([n_qubits - 1 - inst.qubits[0].index], inst.operation.params)

            counter_s += 1
        else:
            my_circ.cx([n_qubits - 1 - inst.qubits[0].index, n_qubits - 1 - inst.qubits[1].index])
            counter_cx += 1
    return my_circ

def qiskit_to_circ_add_one_gate(qiskit_circ, ls, i=0) -> NoisyCircuit:
    n_qubits = qiskit_circ.num_qubits
    my_circ = NoisyCircuit(n_qubits)
    counter_cx = 0
    counter_s = 0

    for inst in qiskit_circ:
        if inst.operation.name in {'u3', 'u'}:
            if counter_s == i:
                my_circ.xflip([n_qubits - 1 - inst.qubits[0].index], 1)
            counter_s += 1
            my_circ.u3([n_qubits - 1 - inst.qubits[0].index], inst.operation.params)
        else:
            my_circ.cx([n_qubits - 1 - inst.qubits[0].index, n_qubits - 1 - inst.qubits[1].index])
            counter_cx += 1
    return my_circ

def qiskit_op_to_matrix(qiskit_op):
    sparse = qiskit_op.to_matrix()
    return np.array(sparse)


def trace(rho, op, op_kron=False):
    size = len(rho.shape)//2
    rho_axes = np.concatenate(((np.arange(size)*2 + 1, np.arange(size)*2)))
    if not op_kron:
        op_axes = np.concatenate((np.arange(size)*2, np.arange(size)*2 + 1,))
    else:
        op_axes = np.arange(2*size)
    return float(np.tensordot(rho, op.reshape([2]*size*2), axes=(rho_axes, op_axes)).real)


def test():
    import numpy as np
    from qiskit import QuantumCircuit, transpile
    import json

    f = open('./data/DoubleParam.json')
    data_ideal = json.load(f)
    E = []
    toto = []
    l = 0
    EE = {} 
    ls = list(range(0,12))
    ls = []
    probs = np.flip(np.geomspace(0.00005, (0.05), 6))
    # probs = [1-i for i in [0.99, 0.999,0.9999, 0.99999, 0.999999, 0.9999999]]
    # probs = [1-i for i in [0.99]]
    
    
    
    for mol in molecules[:]:
        for rep in reps:
            iter_num = [None]*n
            circs, ops, ref_value, init_point, en_hf, circ_prov = get_circs_and_ops(get_cp=True, reps=rep, active_orbitals=mol[2], num_electrons=mol[1], geometry=mol[0], basis=mol[3])
            circs.append(circ_prov.get_yordan_dynamic())
            ops.append(ops[0])
            circs.append(circ_prov.get_zyx_yordan_dynamic())
            ops.append(ops[1])
            k = 0
            print(len(circs))
            for circ, op in zip(circs[:], ops[:]):
                
                k += 1
                id_par = data_ideal[l]['param']

                par = circ.parameters
                cir = circ.assign_parameters({el: id_par[i] for i, el in enumerate(par)})
                # ansatz = qiskit_to_circ_add_one_gate(cir, ls,i)
                ansatz = qiskit_to_mycirc(cir, ls)
                
                rho = np.zeros([2**ansatz.num_qubits, 2**ansatz.num_qubits])
                rho[0][0] = 1
                for name in Noise_ops:
                    E = []
                    for prob in probs[:]:
                        # for i in range(circ.count_ops()['u3']):
                        
                        cir = ansatz.add_noise(prob, name,  single=False, double=True)
                        # cir = ansatz.add_damping(gate_time=prob,time_dec=10, name=name)
                        # print(cir)
                        # print("---------------")
                        # ls = [0,0,0,0]
                        # for inst in cir.instructions:
                        #     if inst.name[0:4] == 'damp':
                        #         for qubit in inst.qubits:
                        #             ls[qubit] += float(inst.name[5:])
                        #     # print(inst)
                        # print(ls)
                        rho_evolved = cir.matrix_run(rho)
                        E.append(trace(rho_evolved, qiskit_op_to_matrix(op), op_kron=True))
                        # E.append(trace(rho_evolved, np.eye(2**ansatz.num_qubits), op_kron=True))

                    EE[name + "_" + str(k)] = E
                l += 1

            for key in EE:
                print(key, ": ", EE[key])
    return EE

EE = test()