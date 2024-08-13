from qiskit_algorithms.optimizers import SPSA ,CG, SLSQP, L_BFGS_B ,COBYLA
from qiskit.quantum_info import Statevector
from qiskit import QuantumCircuit

from Ternary_Tree.optimizer.soap import SOAP

from utils import *

# a0=0.07
# af=0.07
# # af = 0    
# b0 = 0.00001
# bf = 0.00001
# bf = 0
# # n = 12


a0=0.1
af=0.004
# af = 0
b0 = 0.012
bf = 0.0004
# bf = 0
m = 1
def a(i=0):
    def _a(i=0):
        while True:
            i += 1
            _a0 = a0/(0.2*m+1)
            _af = af/(0.2*m+1)
            yield (_a0 - _af)/i**(0.602) + _af
    return _a()
def b():
    def _b(i=0):
        while True:
            i += 1
            _bf = bf/(2*m+1)
            _b0 = b0/(2*m+1)
            yield (_b0 - _bf)/i**(0.101) + _bf
    return _b()

counts = [0]
def call(*args, **kwargs):
    counts[0] += 1
# choice of optimizer
optimizer = SLSQP(maxiter=5000)
optimizer = CG(maxiter=500, eps=0.0001, gtol=0.000000001)
optimizer = L_BFGS_B(maxiter=300, eps=0.00001,ftol=0.000)

optimizer = SPSA(maxiter=800, learning_rate=a, perturbation=b,  callback=call)
# optimizer = SLSQP(maxiter=50, ftol=0)
# optimizer = COBYLA(maxiter=1000, tol=0)
# optimizer = SOAP(500)  

# circuit_init
number_of_reps = 1
# circ_prov = CircuitProvider(reps=number_of_reps, active_orbitals=[0,1,2,3])


def hf(ansatz, op):
    par = ansatz.parameters
    qc = ansatz.assign_parameters({el: 0 for el in par})
    state = [0]*2**qc.num_qubits
    state[0] = 1
    state = Statevector(state)
    state = state.evolve(qc)
    energy = state.expectation_value(op)
    print("nulpar = ", energy)
    return energy

def get_circs_and_ops(**kwargs):
    circ_prov = CircuitProvider(reps=number_of_reps, **kwargs)
    thebest_qc = circ_prov.get_dynamic()
    # init_point[0] =  0.07
    jw_qc = circ_prov.get_jw()
    from numpy.random import uniform
    init_point = [0 for _ in jw_qc.parameters]
    # init_point = [uniform(-0.1, 0.1, 1) for _ in jw_qc.parameters]
    # init_point = None
    jw_lexic_qc = circ_prov.get_jw_lexic()
    bk_qc = circ_prov.get_bk()
    bk_lexic_qc = circ_prov.get_bk_lexic()
    # jw_qc = None
    # jw_lexic_qc = None
    # bk_qc = None
    # bk_lexic_qc = None
    jw_opt_qc = circ_prov.get_jw_opt_ansatz()
    jw_opt_lexic_qc = circ_prov.get_jw_opt_lexic_ansatz()
    # jw_opt_qc = None
    # jw_opt_lexic_qc = None
    en_hf = hf(thebest_qc, ucc_ham(circ_prov.fermionic_op, circ_prov.dynamic_map))
    circs = [thebest_qc, jw_qc, jw_lexic_qc, bk_qc, bk_lexic_qc, jw_opt_qc, jw_opt_lexic_qc]
    ops = [ucc_ham(circ_prov.fermionic_op, circ_prov.dynamic_map), jw_ham(circ_prov.fermionic_op), jw_ham(circ_prov.fermionic_op), bk_ham(circ_prov.fermionic_op), 
        bk_ham(circ_prov.fermionic_op), ucc_ham(circ_prov.fermionic_op, circ_prov.jw_opt_map), ucc_ham(circ_prov.fermionic_op, circ_prov.jw_opt_map)]
    ref_value = numpy_energy(circ_prov.fermionic_op, circ_prov.ucc)
    return circs, ops, ref_value, init_point, en_hf

# circs, ops, ref_value, init_point,en_hf = get_circs_and_ops()

def hf(ansatz, op):
    par = ansatz.parameters
    qc = ansatz.assign_parameters({el: 0 for el in par})
    state = [0]*2**qc.num_qubits
    state[0] = 1
    state = Statevector(state)
    state = state.evolve(qc)
    energy = state.expectation_value(op)
    print("nulpar = ", energy)
    return energy

# en_hf = hf(thebest_qc)
id_ener = []
noise_ener = []
k = 0
n = 7

id_ener = [[] for i in range(n)]
noise_ener = [[] for i in range(n)]

_counts = [0 for i in range(n)]
_values = [0 for i in range(n)]


basises = ['sto-3g', '6-31g']
# basises = ['sto-3g']
geometry='H 0 0 0; Li 0 0 0.7349'
active_orbitals=[0, 1]
num_electrons=(1, 1)
# SPSA(maxiter=3000, learning_rate=a, perturbation=b, callback=call, perturbation_dims=10),
optimizers = [(SLSQP(maxiter=200, ftol=0), 'SLSQP'), 
                (COBYLA(maxiter=500, tol=0), 'COBYLA'),
                (SPSA(maxiter=500, learning_rate=a, perturbation=b,  callback=call), 'SPSA'),
                (SOAP(4000), "SOAP")]

# optimizers = [(SPSA(maxiter=500, learning_rate=a, perturbation=b,  callback=call), 'SPSA')]
optimizers = [(SOAP(1000), "SOAP")]

molecules = [['H 0 0 0; H 0 0 0.7349', (1,1), [0,1], basises[0]], ['H 0 0 0; H 0 0 0.7349', (1,1), [0,1,2,3], basises[1]], 
                ['H 0 0 0; Li 0 0 1.5459', (2,2), [0,1,2,5], basises[0]]]
molecules = [('H 0 0 0; H 0 0 0.7349', (1,1), [0,1], basises[0])] 

def write_data():
    data = []
    import json 
    id_en = [None]*n
    with open('numbers.json', 'w') as file:
        for mol in molecules[:]:
            # for basis in basises[:]:
            for opt in optimizers[:]:
                # for orb in active_orbitals[:]:
                iter_num = [None]*n
                # print("herheehr")
                circs, ops, ref_value, init_point, en_hf = get_circs_and_ops(active_orbitals=mol[2], num_electrons=mol[1], geometry=mol[0], basis=mol[3])
                for index, obj in enumerate(zip(circs[:], ops[:])):
                    circ, op = obj
                    print(circ.count_ops())
                    print(init_point)
                    if index == 0: 
                        id_en[index], _, values, _ = VQE_energy_with_noise(circ, op, ref_value, opt[0], init_point[:1],  en_hf)
                    else: 
                        id_en[index], _, values, _ = VQE_energy_with_noise(circ, op, ref_value, opt[0], init_point,  en_hf)

                    # id_en[index], _, values, _ = ideal_energy(circ, op, ref_value, opt[0], init_point, en_hf)
                    for index2, val in enumerate(values):
                        if abs(val - ref_value) < 0.00016:
                            iter_num[index] = index2 
                            break
                data.append({"mol": mol, "opt": opt[1], "act_orb": mol[2], "hf": en_hf.real, "ref_ener": ref_value, "energy": id_en[:], "iter_num": iter_num})
        json.dump(data, file, indent=4)

X = np.linspace(0, 1, 500)
E = []
toto = []
for mol in molecules[:]:
    for opt in optimizers[:]:
        iter_num = [None]*n
        circs, ops, ref_value, init_point, en_hf = get_circs_and_ops(active_orbitals=mol[2], num_electrons=mol[1], geometry=mol[0], basis=mol[3])
        par = circs[0].parameters
        circ = circs[0]
        op = ops[0]
        a,b,c,d = ideal_energy(circ, op, ref_value, opt[0], init_point,  en_hf)
        print(a,c[-1],d[-1])

        for prob in [0.9, 0.99,0.999,0.9999, 0.99999, 0.999999]:
            est = get_device_noise_estimator(4,prob)
            x = 0.11174828549146505
            ansatz = circs[0].assign_parameters({el: x for el in par})
            energy = est.run(ansatz, ops[0]).result()
            E.append(energy.values[0])
            # E = []
            # for x in X:
            #     ansatz = circs[0].assign_parameters({el: x for el in par})
            #     energy = est.run(ansatz, ops[0]).result()
            #     # print("nulpar = ", energy)
                # E.append(energy.values[0])

print(E)
import pylab    
pylab.xlabel("error")
pylab.ylabel("Energy")
# pylab.title("Convergence with no noise")
# pylab.ylim(bottom=-1.90, top=-1.30)
# pylab.rcParams["figure.figsize"] = (12, 4)

def print_last_data(counts, values, name="hello"):
    pylab.plot(counts, values, label=name)
    
# print_last_data(X, E, name=str(i))
# print(min(E))
# pylab.legend()
# pylab.savefig("toto2")


a = [[-1.0031434741388943, -1.6905585840512296, -1.8393852284472576, -1.855560165924037, -1.8571912570353855],
    [-1.444400023964278, -1.803600922873114, -1.8517715150181235, -1.856809923221062, -1.8573163472207639],
    [-1.003119841516724, -1.690559509101338, -1.839385656780147, -1.855560541931011, -1.8571916277828961]]
probs = [1 - prob for prob in  [0.9, 0.99,0.999,0.9999, 0.99999]]
labels = ["12 gates", "38 gates", "fixed parms"]
for i in range(3):
    pylab.plot(probs[1:], a[i][1:], label=labels[i])
pylab.xscale('log')
pylab.legend()
pylab.savefig("toto3")



