from Ternary_Tree.ucc.upgccsd import UpGCCSD, UCCSD
from qiskit.quantum_info import Statevector
from qiskit import transpile
from Ternary_Tree.ucc.abstractucc import Molecule
from Ternary_Tree.qiskit_interface.circuit_provider import CircuitProvider, Circuits
import numpy as np
import logging

# logger.setLevel(logging.INFO)

# logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
logger.addHandler(handler)
logger.setLevel(logging.INFO)

basis_gates = ["u3", "cx"]
# basis_gates1 = [ "h","rz","cx"]
# basis_gates1 = basis_gates
# basis_gates2 = ["u3", 'cx']

basis = '6-311G'
# geometry='Be 0 0 0; Be 0 0 0.739'
geometry='H 0 0 0; H 0 0 0.739'
active_orbitals=[0, 1]
num_electrons=(1, 1)
H2_4 = Molecule(geometry='H 0 0 0; H 0 0 0.7349', num_electrons=(1,1), active_orbitals=[0,1], basis='sto-3g')
H2_8 = Molecule(geometry='H 0 0 0; H 0 0 0.7349', num_electrons=(1,1), active_orbitals=[0,1,2,3], basis='6-31g')

def base_test():
    ucc = UpGCCSD(molecule=H2_4)
    cirq, mtoq = ucc.swap2xn(1)
    # cirq = ucc.swap_gen(1, LadExcNames.YORDAN())
    print(cirq)
    params = cirq.parameters
    params = {par: 1 for par in params}
    cirq.assign_parameters(params, inplace=True)
    # print(cirq.decompose())
    state = Statevector.from_label("0000")
    res = state.evolve(cirq)
    probs = res.probabilities_dict()
    for key in probs:
        if abs(probs[key]) > 0.000001:
            print(key, ": ", probs[key])
    print(probs)

def res(angle):
    active_orbitals = [j for j in range(2)]
    ucc = UpGCCSD(geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons)
    cirq = ucc.swap_gen()
    state = Statevector.from_label("0000")
    par = cirq.parameters[0]
    # cirq.assign_parameters({par:angle})
    cirq1 = cirq.bind_parameters([0])
    print(state.evolve(cirq1).probabilities_dict())
    cirq = cirq.bind_parameters([angle])
    state = state.evolve(cirq)
    print(cirq.decompose())
    print(state.probabilities_dict())
    
def scaleability():
    mol = Molecule(geometry='Li 0 0 0; Li 0 0 1.5459', num_electrons=(2,2), active_orbitals=[0,1,2,5], basis='6-311g')
    circ_names = Circuits.get_circs_names()[:2]
    double = {circ: [] for circ in circ_names}
    single = {circ: [] for circ in circ_names}
    depth = {circ: [] for circ in circ_names}
    num = {circ: [] for circ in circ_names}
    pauli_num = {circ: [] for circ in circ_names}
    for i in range(4,14,2):
        mol = Molecule(geometry='Li 0 0 0; Li 0 0 1.5459', num_electrons=(2,2), active_orbitals=list(range(i)), basis='6-311g')
        try:
            logger.info(f"{i=}")
            circ_prov = CircuitProvider(reps=1, molecule=mol)
            for name in circ_names:
                circ = circ_prov.get_circ(name)[1]
                circ = transpile(circ, basis_gates=basis_gates, optimization_level=3)
                # print(cirq)
                # num.append(i*4)
                # print(circ.count_ops())
                double[name].append(circ.count_ops()["cx"])
                single[name].append(circ.count_ops()["u3"])
                depth[name].append(circ.depth())
                pauli_num[name].append(i*(i-1)/2*10)
                # pauli_num.append(2*i*(2*i-1)/2*10)
                # print(name, ":", circ.count_ops())
                # print(name, ":", circ.depth())
        except IndexError:
            pass
    for name in circ_names:
        print(name, " double:", double[name])
        print(name, " single:", single[name])
        print(name, "depth:", depth[name])

    print(Circuits.get_circs_names()[0], "pauli num:", pauli_num[name])
    # print(num)
    # print(double)
    # print(single)
    # print(depth)
    # print(pauli_num)

if __name__ == "__main__":
    # scaleability()
    jw  = np.array([413, 1340, 3087, 5910, 10065])
    bk  = np.array([394, 1528, 2881, 6064, 9546])
    p =   np.array([60.0, 150.0, 280.0, 450.0, 660.0])
    q = np.array([8,12,16,20,24])
    print(" ".join([ f"({el[0]}, {el[1]})" for el in zip(q, jw/p)]))
    print(" ".join([ f"({el[0]}, {el[1]})" for el in zip(q, bk/p)]))
    # ucc = UCCSD(molecule=H2_8)
    # circ, mapping = ucc.cyclic_algorithm()
    # pars = circ.parameters[:]
    # qc = circ.to_excitaions(pars)
    # print(circ)
    # na = 3
    # cycl_double = []
    # cycl_depth = []
    # for i in range(4,16,2):
    #     nb = i - na
    #     s = na*nb*2
    #     sg = i*(i - 1)
    #     d = 2*na*(na-1)*nb*(nb-1)/4 + na*na*nb*nb
    #     dg = 2*i*(i - 1)*i*(i-1)/4 + i*(i-1)*i*(i-1)
    #     # print(s)
    #     # print(d)
    #     print(s*2 + d*8)
    #     swap = (4/6*i*i*i + i*i  + 1/3*i - 14)*4
    #     cycl_double.append(swap + 12*d)
    #     cycl_depth.append(swap/4*7 + 10*d)
    # dy = [109, 163, 217, 271, 325, 379, 433, 487, 541, 595, 649, 703]
    # ds = [85, 127, 169, 211, 253, 295, 337, 379, 421, 463, 505, 547]
    # gy = [108, 270, 504, 810, 1188, 1638, 2160, 2754, 3420, 4158, 4968, 5850]
    # gs = [102, 255, 476, 765, 1122, 1547, 2040, 2601, 3230, 3927, 4692, 5525] 
    # pauli = [60.0, 150.0, 280.0, 450.0, 660.0, 910.0, 1200.0, 1530.0, 1900.0, 2310.0, 2760.0, 3250.0]
    # depth_yor = ""
    # depth_short = ""
    # gate_short = ""
    # gate_yor = ""
    # for index, i in enumerate(range(4, 12, 2)):
    #     depth_yor += f"({2*i}, {dy[index]}) "
    #     depth_short += f"({2*i}, {ds[index]}) "
    #     gate_yor += f"({2*i}, {gy[index]/pauli[index]:.5f}) "
    #     gate_short += f"({2*i}, {gs[index]/pauli[index]:.5f}) "
    # print(f"{depth_yor=}")
    # print(f"{depth_short=}")
    # print(f"{gate_yor=}")
    # print(f"{gate_short=}")
        
    # print("cycl_depth",cycl_depth)
    # print("cycl_double",cycl_double)