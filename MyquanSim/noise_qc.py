from .ideal_qc import QuantumCircuit
from MyquanSim.noise_instruction import Depolarization,Xflip,Yflip,Zflip, AmpDamp, PhaseDamp
import numpy as np
class NoisyCircuit(QuantumCircuit):
    
    def matrix_run(self, rho_init):
        rho_init = rho_init.reshape([2,2]*self.num_qubits)
        rho = rho_init
        for instr in self.instructions:
            rho = instr.evolve_matrix(rho)
        return rho

    def add_noise(self, prob, name, single=True, double=False):
        qc = NoisyCircuit(self.num_qubits)
        instr_new = []
        dic = {"D": Depolarization, "X": Xflip, "Y": Yflip, "Z": Zflip}
        for instr in self.instructions:
            if instr.name == "U3" and single:
                instr_new.append(dic[name](instr.qubits, prob))
            if instr.name in {"CX", "CZ"} and double:
                for qubit in instr.qubits:
                    instr_new.append(dic[name]([qubit], prob))
            instr_new.append(instr)
        qc.instructions = instr_new
        return qc

    def add_damping(self, time_dec=25, gate_time=0.1,name='A'):
        qc = NoisyCircuit(self.num_qubits)
        instr_new = []
        free_layer = [0 for i in range(self.num_qubits)]
        dic = {"A": AmpDamp, "P": PhaseDamp}
        # pseudo_free_layer = [0 for i in range(self.num_qubits)]
        damp_time = [0 for i in range(self.num_qubits)]
        def update_layer(qubits):
            max_layer = max([free_layer[i] for i in qubits]) 
            for i in qubits:
                damp_time[i] += (max_layer - free_layer[i] )*gate_time 
                free_layer[i] = max_layer +1
            return max_layer + 1

        for instr in self.instructions:
            
            max_layer = update_layer(instr.qubits)
            for i in instr.qubits:
                instr_new.append(dic[name]([i], np.exp(-damp_time[i]/time_dec)))
                damp_time[i] = gate_time
            instr_new.append(instr)
        qc.instructions = instr_new
        return qc

    def dep(self, qubits, prob):
        self.add_instr(Depolarization(qubits, prob))

    def xflip(self, qubits, prob):
        self.add_instr(Xflip(qubits, prob))

    def yflip(self, qubits, prob):
        self.add_instr(Yflip(qubits, prob))

    def zflip(self, qubits, prob):
        self.add_instr(Zflip(qubits, prob))