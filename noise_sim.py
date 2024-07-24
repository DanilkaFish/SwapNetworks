from qiskit_algorithms.optimizers import SPSA ,CG, SLSQP, L_BFGS_B
from qiskit.quantum_info import Statevector
from qiskit import QuantumCircuit


from utils import *

a0=0.05
af=0.0005
# af = 0
b0 = 0.1
bf = 0.0002
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
# optimizer = SPSA(maxiter=3000, learning_rate=a, perturbation=b, callback=call, perturbation_dims=10)
optimizer = SLSQP(maxiter=200, ftol=0)

# circuit_init
number_of_reps = 2
circ_prov = CircuitProvider(reps=number_of_reps)
thebest_qc = circ_prov.get_dynamic()
init_point = [0 for _ in thebest_qc.parameters]
# init_point[0] =  0.07
jw_qc = circ_prov.get_jw()
jw_lexic_qc = circ_prov.get_jw_lexic()
bk_qc = circ_prov.get_bk()
bk_lexic_qc = circ_prov.get_bk_lexic()
# jw_opt_qc = get_jw_opt_ucc_ansatz(reps=number_of_reps)
# jw_opt_lexic_qc = get_jw_opt_lexic_ucc_ansatz(reps=number_of_reps)
jw_opt_qc = None
jw_opt_lexic_qc = None

circs = [thebest_qc, jw_qc, jw_lexic_qc, bk_qc, bk_lexic_qc, jw_opt_qc, jw_opt_lexic_qc]
ops = [ucc_ham(circ_prov.fermionic_op, circ_prov.jw_opt_map), jw_ham(circ_prov.fermionic_op), jw_ham(circ_prov.fermionic_op), bk_ham(circ_prov.fermionic_op), 
       bk_ham(circ_prov.fermionic_op), ucc_ham(circ_prov.fermionic_op, circ_prov.jw_opt_map), ucc_ham(circ_prov.fermionic_op, circ_prov.jw_opt_map)]
ref_value = numpy_energy(circ_prov.fermionic_op, circ_prov.ucc)

def hf(ansatz, op=ucc_ham(circ_prov.fermionic_op, circ_prov.jw_opt_map)):
    par = ansatz.parameters
    qc = ansatz.assign_parameters({el: 0 for el in par})
    state = [0]*2**qc.num_qubits
    state[0] = 1
    state = Statevector(state)
    state = state.evolve(qc)
    energy = state.expectation_value(op)
    print("nulpar = ", energy)
    return energy

en_hf = hf(thebest_qc)
id_ener = []
noise_ener = []
k = 0
n = 5
id_ener = [[] for i in range(n)]
noise_ener = [[] for i in range(n)]

circ = circs[3]
op = ops[3]


a0=0.07
af=0.07
# af = 0    
b0 = 0.00001
bf = 0.00001
# bf = 0
# n = 12

_counts = [0 for i in range(n)]
_values = [0 for i in range(n)]
# for k in range(n):
    # m = (k + 4) * 1
    # m = 1
for circ, op in zip(circs[:-2], ops[:-2]):
    print(k)
    counts[0] = 0
    # optimizer = SPSA(maxiter=3000, learning_rate=a, perturbation=b, callback=call, perturbation_dims=21)

    id_en, _counts[k], _values[k], params = ideal_energy(circ, op, ref_value, optimizer, init_point, en_hf)
    print("counts", counts[0])
    print(params[-1])
    print(params[3000])

    print(len(params[0]))
    id_ener[k].append(id_en)
    counts[0] = 0
    k+= 1
    # id_en, _, _, _ = VQE_energy_with_noise(circ, op, ref_value, optimizer, init_point, en_hf) 
    # print("counts", counts[0])
    # noise_ener[k].append(id_en)

print(id_ener)
print(noise_ener)

import pylab    
pylab.xlabel("Eval count")
pylab.ylabel("Energy")
pylab.title("Convergence with no noise")
pylab.ylim(bottom=-1.872, top=-1.830)
pylab.rcParams["figure.figsize"] = (12, 4)

def print_last_data(counts, values, name="hello"):
    pylab.plot(counts, values, label=name)

for i in range(n):
    print_last_data(_counts[i], _values[i], name=str(i))
pylab.legend()
pylab.savefig("toto2")
