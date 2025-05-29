import logging
import matplotlib.pyplot as plt
import numpy as np
from qiskit_algorithms.minimum_eigensolvers import VQE
from qiskit_nature.second_q.mappers import JordanWignerMapper
import qiskit_algorithms.optimizers as optim
from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
from qiskit.primitives import Estimator  # , StatevectorEstimator
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_algorithms import VQE
import qiskit_algorithms.optimizers as optim
from qiskit.primitives import Estimator
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from pyscf import gto, scf, fci
import pennylane as qml
from ExcitationSolve.excitationsolve.excitation_solve_qiskit import ExcitationSolveQiskit

if __name__ == "__main__":
    molname = "H2"
    bondlength = 0.742

    basis = "STO-3G"

    # atom = "; ".join([f"{a} {', '.join([str(x) for x in p.tolist()])}" for a, p in zip(symbols, geometry)])
    symbols = ["H", "H"]
    geometry = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]])
    molecule = qml.qchem.Molecule(symbols, geometry, charge=0, mult=1, basis_name='sto-3g')
    H, qubits = qml.qchem.molecular_hamiltonian(molecule, method='pyscf')
    # data = data[0]
    # molecule = data.molecule
    # H = data.hamiltonian
    electrons = molecule.n_electrons
    charge = molecule.charge
    qubits = molecule.n_orbitals * 2
    fci_energy = -1.87
    unit = DistanceUnit.BOHR  # Angstrom or Bohr

    mol_pyscf = gto.M(atom='H 0 0 0; H 0 0 0.7349', basis=basis, charge=charge, unit=unit.value)
    electrons = mol_pyscf.nelectron
    rhf = scf.RHF(mol_pyscf)
    hf_energy = rhf.kernel()

    fci_calc = fci.FCI(mol_pyscf, rhf.mo_coeff)
    fci_energy, ci_vector = fci_calc.kernel()
    print(f"FCI Energy: {fci_energy} Ha")

    driver = PySCFDriver(
        atom='H 0 0 0; H 0 0 0.7349',
        basis=basis,
        charge=charge,
        # unit=unit,
    )

    problem = driver.run()

    mapper = JordanWignerMapper()
    energy = problem.hamiltonian
    offset = np.sum(list(energy.constants.values()))
    fermionic_op = energy.second_q_op()
    pauli_sum_op = mapper.map(fermionic_op)
    print(pauli_sum_op)
    logging.info("pauli_sum_op.num_qubits: %f", pauli_sum_op.num_qubits)

    initial_state = HartreeFock(
        num_spatial_orbitals=problem.num_spatial_orbitals,
        num_particles=problem.num_particles,
        qubit_mapper=mapper,
    )

    ansatz = UCCSD(
        num_spatial_orbitals=problem.num_spatial_orbitals,
        num_particles=problem.num_particles,
        qubit_mapper=mapper,
        reps=1,
        initial_state=initial_state,
        generalized=True,
        preserve_spin=True,
        include_imaginary=False,
    )
    # Fix UCCSD operator order: first doubles, then singles
    doubles = [op for (op, excitation) in zip(ansatz.operators, ansatz.excitation_list) if len(excitation[0]) == 2]
    singles = [op for (op, excitation) in zip(ansatz.operators, ansatz.excitation_list) if len(excitation[0]) == 1]
    ansatz.operators = doubles + singles
    ansatz._invalidate()
    # print(ansatz.decompose(reps=4))
    # print(ansatz.decompose(reps=2))
    logging.info("ansatz.num_qubits: %i", ansatz.num_qubits)

    optimizers = [
        optim.L_BFGS_B(maxiter=1e6, ftol=1e-6),
        optim.COBYLA(maxiter=1e6, tol=1e-6),
        ExcitationSolveQiskit(maxiter=10, hf_energy=0, save_parameters=True),
    ]
    estimator = Estimator()
    # estimator = StatevectorEstimator()

    logging.info("HF-result: %s", estimator.run(initial_state, pauli_sum_op).result().values)

    counts_dict = {}
    values_dict = {}

    for optimizer in optimizers[-1:]:
        print(f"Using {type(optimizer).__name__}...")

        counts = []
        values = []

        def store_intermediate_result(eval_count, parameters, mean, std):
            print(
                f"Optimizer evaluation #{eval_count}, mean: {mean}, Diff. to ref.:{np.abs(mean + offset - fci_energy)}",
                end="\r",
            )
            counts.append(eval_count)
            values.append(mean + offset)

        solver = VQE(
            estimator,
            ansatz,
            optimizer,
            callback=store_intermediate_result,
            initial_point=np.zeros((ansatz.num_parameters,)),
        )
        logging.info("VQE defined")
        logging.info("Solving VQE...")
        vqe_result = solver.compute_minimum_eigenvalue(pauli_sum_op)
        logging.info("Solved")
        elec_struc_result = problem.interpret(vqe_result)

        print(elec_struc_result)

        if isinstance(optimizer, ExcitationSolveQiskit):
            counts = optimizer.nfevs
            values = optimizer.energies + offset

        counts_dict[type(optimizer).__name__] = np.array(counts)
        values_dict[type(optimizer).__name__] = np.array(values)

    plt.figure()
    for key, vals in values_dict.items():
        counts = counts_dict[key]

        plt.plot(counts, np.abs(vals - fci_energy), marker=".", linestyle="-", label=key)
    plt.yscale("log")
    plt.grid()
    plt.legend()
