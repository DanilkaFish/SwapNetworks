

# def get_thebest_ucc_ansatz(reps=1, 
#                            geometry=geometry, 
#                            basis=basis, 
#                            active_orbitals=active_orbitals, 
#                            num_electrons=num_electrons):
#     ucc = UpUCCSDG(geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons)
#     fermionic_op = ucc.mol.hamiltonian.second_q_op()
#     ansatz = ucc.get_parametrized_circuit(reps)
#     tree_map = copy(ucc.tt)
#     ansatz = transpile(ansatz, basis_gates=basis_gates, optimization_level=3)
#     print(ansatz.count_ops())
#     # print(ansatz.decompose(reps=4))
#     return ucc, ansatz, fermionic_op, tree_map


def get_jw_ucc_ansatz(reps=1, geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons):
    ucc = UpUCCSDG(geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons)
    al = ucc.get_alpha_excitations()
    be = ucc.get_beta_excitations()
    do = ucc.get_double_excitations()
    def my_gen(**kwargs):
        return my_generation(al,be, do, **kwargs)
    qubit_mapper=JordanWignerMapper()
    qc = UCC(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), excitations=my_gen, qubit_mapper=qubit_mapper, 
                  initial_state=HartreeFock(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), qubit_mapper))
    for i in range(reps-1):
        new_ucc = UCC(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), excitations=my_gen, qubit_mapper=qubit_mapper)
        theta = ParameterVector("θ" + str(i), new_ucc.num_parameters)
        new_ucc.assign_parameters(theta, inplace=True)
        qc = qc.compose(new_ucc)

    ansatz = transpile(qc, basis_gates=basis_gates, optimization_level=3)
    print(ansatz.count_ops())
    return ansatz


def get_jw_lexic_ucc_ansatz(reps=1, geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons):
    ucc = UpUCCSDG(geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons)
    al = ucc.get_alpha_excitations()
    be = ucc.get_beta_excitations()
    do = ucc.get_double_excitations()
    def my_gen(**kwargs):
        return my_generation(al,be, do, **kwargs)
    qubit_mapper=JordanWignerMapper()
    qc = HartreeFock(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), qubit_mapper)
    qc._build()
    
    for i in range(reps):
        new_ucc = UCC(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), excitations=my_gen, qubit_mapper=qubit_mapper)
        theta = ParameterVector("θ" + str(i), new_ucc.num_parameters)
        new_ucc.assign_parameters(theta, inplace=True)
        parr = []
        for gate in new_ucc.decompose(reps=2):
            parr.append([pauliString(gate.operation.name[-1-ucc.n_qubits:-1], 1.0)])
            shuffle(parr)
            nq = len(parr[0][0])
            # length = nq//2 # `length' is a hyperparameter, and can be adjusted for best performance
            a1 = gate_count_oriented_scheduling(parr)
        for pauli in a1:
            for gate in new_ucc.decompose(reps=2):
                if gate.operation.name[-1-ucc.n_qubits:-1] == str(pauli[0][0]):
                    qc = qc.compose(gate[0])
    ansatz = transpile(qc, basis_gates=basis_gates, optimization_level=3)
    print(ansatz.count_ops())
    return ansatz


def get_bk_ucc_ansatz(reps=1, geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons):
    ucc = UpUCCSDG(geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons)
    al = ucc.get_alpha_excitations()
    be = ucc.get_beta_excitations()
    do = ucc.get_double_excitations()
    def my_gen(**kwargs):
        return my_generation(al,be, do, **kwargs)
    qubit_mapper=BravyiKitaevMapper()
    qc = UCC(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), excitations=my_gen, qubit_mapper=qubit_mapper, 
                  initial_state=HartreeFock(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), qubit_mapper))
    for i in range(reps-1):
        new_ucc = UCC(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), excitations=my_gen, qubit_mapper=qubit_mapper)
        theta = ParameterVector("θ" + str(i), new_ucc.num_parameters)
        new_ucc.assign_parameters(theta, inplace=True)
        qc = qc.compose(new_ucc)
    ansatz = transpile(qc, basis_gates=basis_gates, optimization_level=3)
    print(ansatz.count_ops())
    return ansatz

def get_bk_lexic_ucc_ansatz(reps=1, geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons):
    ucc = UpUCCSDG(geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons)
    al = ucc.get_alpha_excitations()
    be = ucc.get_beta_excitations()
    do = ucc.get_double_excitations()
    def my_gen(**kwargs):
        return my_generation(al,be, do, **kwargs)
    qubit_mapper=BravyiKitaevMapper()
    
    qc = HartreeFock(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), qubit_mapper)
    qc._build()
    for i in range(reps):
        new_ucc = UCC(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), excitations=my_gen, qubit_mapper=qubit_mapper)
        theta = ParameterVector("θ" + str(i), new_ucc.num_parameters)
        new_ucc.assign_parameters(theta, inplace=True)
        parr = []
        for gate in new_ucc.decompose(reps=2):
            parr.append([pauliString(gate.operation.name[-1-ucc.n_qubits:-1], 1.0)])
            shuffle(parr)
            nq = len(parr[0][0])
            # length = nq//2 # `length' is a hyperparameter, and can be adjusted for best performance
            a1 = gate_count_oriented_scheduling(parr)
        for pauli in a1:
            for gate in new_ucc.decompose(reps=2):
                if gate.operation.name[-1-ucc.n_qubits:-1] == str(pauli[0][0]):
                    qc = qc.compose(gate[0])
    ansatz = transpile(qc, basis_gates=basis_gates, optimization_level=3)
    print(ansatz.count_ops())
    return ansatz


def get_jw_opt_ucc_ansatz(reps=1, geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons):
    ucc = UpUCCSDG(geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons)
    qc = ucc.get_parametrized_circuit(number_of_layers=1)
    al = ucc.get_alpha_excitations()
    be = ucc.get_beta_excitations()
    do = ucc.get_double_excitations()
    def my_gen(**kwargs):
        return my_generation(al,be, do, **kwargs)
    qubit_mapper=TernaryTreeMapper(copy(ucc.tt))
    qc = UCC(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), excitations=my_gen, qubit_mapper=qubit_mapper, 
                  initial_state=HartreeFock(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), qubit_mapper))
    for i in range(reps-1):
        new_ucc = UCC(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), excitations=my_gen, qubit_mapper=qubit_mapper)
        theta = ParameterVector("θ" + str(i), new_ucc.num_parameters)
        new_ucc.assign_parameters(theta, inplace=True)
        qc = qc.compose(new_ucc)
    ansatz = transpile(qc, basis_gates=basis_gates, optimization_level=3)
    return ansatz

def get_jw_opt_lexic_ucc_ansatz(reps=1, geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons):
    ucc = UpUCCSDG(geometry=geometry, basis=basis, active_orbitals=active_orbitals, num_electrons=num_electrons)
    qc = ucc.get_parametrized_circuit(number_of_layers=1)
    al = ucc.get_alpha_excitations()
    be = ucc.get_beta_excitations()
    do = ucc.get_double_excitations()
    def my_gen(**kwargs):
        return my_generation(al,be, do, **kwargs)
    qubit_mapper=TernaryTreeMapper(ucc.tt)
    qc = HartreeFock(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), qubit_mapper)
    qc._build()
    
    for i in range(reps):
        new_ucc = UCC(ucc.num_spatial_orbitals, (ucc.num_alpha, ucc.num_beta), excitations=my_gen, qubit_mapper=qubit_mapper)
        theta = ParameterVector("θ" + str(i), new_ucc.num_parameters)
        new_ucc.assign_parameters(theta, inplace=True)
        parr = []
        for gate in new_ucc.decompose(reps=2):
            parr.append([pauliString(gate.operation.name[-1-ucc.n_qubits:-1], 1.0)])
            shuffle(parr)
            nq = len(parr[0][0])
            # length = nq//2 # `length' is a hyperparameter, and can be adjusted for best performance
            a1 = gate_count_oriented_scheduling(parr)
        for pauli in a1:
            for gate in new_ucc.decompose(reps=2):
                if gate.operation.name[-1-ucc.n_qubits:-1] == str(pauli[0][0]):
                    qc = qc.compose(gate[0])
    ansatz = transpile(qc, basis_gates=basis_gates, optimization_level=3)
    return ansatz