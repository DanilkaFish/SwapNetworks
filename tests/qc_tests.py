import unittest
from MyquanSim.ideal_qc import QuantumCircuit
from MyquanSim.noise_qc import NoisyCircuit
from MyquanSim.instructions import *
from MyquanSim.noise_instruction import *
import numpy as np

def prod3(psi1, psi2, psi3):
    psi = np.tensordot(psi1, psi2, axes=0)
    psi = np.tensordot(psi, psi3, axes=0)
    return psi

def prod5(psi1, psi2, psi3, psi4, psi5):
    psi = np.tensordot(psi1, psi2, axes=0)
    psi = np.tensordot(psi, psi3, axes=0)
    psi = np.tensordot(psi, psi4, axes=0)
    psi = np.tensordot(psi, psi5, axes=0)
    return psi

def trace(rho, op, op_kron=False):
    size = len(rho.shape)//2
    rho_axes = np.concatenate(((np.arange(size)*2 + 1, np.arange(size)*2)))
    if not op_kron:
        op_axes = np.concatenate((np.arange(size)*2, np.arange(size)*2 + 1,))
    else:
        op_axes = np.arange(2*size)
    return np.tensordot(rho, op, axes=(rho_axes, op_axes))

class test_Instruction(unittest.TestCase):
    def setUp(self):
        self.rho_Z_0 = np.tensordot(np.array([1, 0]), np.array([1, 0]), axes=0)
        a = np.array([0, 1])
        self.rho_Z_1 = np.tensordot(a, a.conj(), axes=0)
        a = np.array([1, 1j])/np.sqrt(2)
        self.rho_Y_0 = np.tensordot(a, a.conj(), axes=0)
        a = np.array([1, -1j])/np.sqrt(2)
        self.rho_Y_1 = np.tensordot(a, a.conj(), axes=0)
        a = np.array([1, 1])/np.sqrt(2)
        self.rho_X_0 = np.tensordot(a, a.conj(), axes=0)
        a = np.array([1, -1])/np.sqrt(2)
        self.rho_X_1 = np.tensordot(a, a.conj(), axes=0)
    
    def test_one_xgate(self):
        X = Xgate([0], np.pi)
        self.assertEqual(np.round(np.trace(self.rho_X_0 @ X.evolve_matrix(self.rho_X_0)), 5), 1)
        self.assertEqual(np.round(np.trace(self.rho_X_1 @ X.evolve_matrix(self.rho_X_1)), 5), 1)
        self.assertEqual(np.round(np.trace(self.rho_Z_0 @ X.evolve_matrix(self.rho_Z_0)), 5), 0)
        self.assertEqual(np.round(np.trace(self.rho_Z_1 @ X.evolve_matrix(self.rho_Z_1)), 5), 0)
        self.assertEqual(np.round(np.trace(self.rho_Y_0 @ X.evolve_matrix(self.rho_Y_0)), 5), 0)
        self.assertEqual(np.round(np.trace(self.rho_Y_1 @ X.evolve_matrix(self.rho_Y_1)), 5), 0)
        
        self.assertEqual(np.round(trace(self.rho_X_0, X.evolve_matrix(self.rho_X_0)), 5), 1)
        self.assertEqual(np.round(trace(self.rho_X_1, X.evolve_matrix(self.rho_X_1)), 5), 1)
        self.assertEqual(np.round(trace(self.rho_Z_0, X.evolve_matrix(self.rho_Z_0)), 5), 0)
        self.assertEqual(np.round(trace(self.rho_Z_1, X.evolve_matrix(self.rho_Z_1)), 5), 0)
        self.assertEqual(np.round(trace(self.rho_Y_0, X.evolve_matrix(self.rho_Y_0)), 5), 0)
        self.assertEqual(np.round(trace(self.rho_Y_1, X.evolve_matrix(self.rho_Y_1)), 5), 0)
        self.assertEqual(np.round(trace(self.rho_Y_0, X.evolve_matrix(self.rho_Y_1)), 5), 1)
        

    def test_multi_xgate(self):
        two_rho_Z = np.tensordot(self.rho_Z_0, self.rho_Z_0, axes=0)
        two_rho_X = np.tensordot(self.rho_X_0, self.rho_X_0, axes=0)
        three_rho_Z = np.tensordot(two_rho_Z, self.rho_Z_0, axes=0)
        three_rho_X = np.tensordot(two_rho_X, self.rho_X_0, axes=0)
        gate0 = Xgate([0], np.pi)
        gate1 = Xgate([1], np.pi)
        gate2 = Xgate([2], np.pi)

        correctZ = gate0.ops() @ self.rho_Z_0 @ gate0.ops().conj().T  
        evolvedZ2_0 = np.tensordot(correctZ, self.rho_Z_0,  axes=0)
        evolvedZ2_1 = np.tensordot(self.rho_Z_0, correctZ, axes=0)
        evolvedZ2_2 = np.tensordot(self.rho_Z_0, np.tensordot(self.rho_Z_0, correctZ, axes=0), axes=0)

        self.assertEqual(np.round(trace(evolvedZ2_1, gate1.evolve_matrix(two_rho_Z)), 5), 1)
        self.assertEqual(np.round(trace(evolvedZ2_0, gate0.evolve_matrix(two_rho_Z)), 5), 1)
        self.assertEqual(np.round(trace(evolvedZ2_2, gate2.evolve_matrix(three_rho_Z)), 5), 1)

        correctX = gate0.ops() @ self.rho_X_0 @ gate0.ops().conj().T  
        evolvedX2_0 = np.tensordot(correctX, self.rho_X_0,  axes=0)
        evolvedX2_1 = np.tensordot(self.rho_X_0, correctX, axes=0)
        evolvedX2_2 = np.tensordot(self.rho_X_0, np.tensordot(self.rho_X_0, correctX, axes=0), axes=0)

        self.assertEqual(np.round(trace(evolvedX2_1, gate1.evolve_matrix(two_rho_X)), 5), 1)
        self.assertEqual(np.round(trace(evolvedX2_0, gate0.evolve_matrix(two_rho_X)), 5), 1)
        self.assertEqual(np.round(trace(evolvedX2_2, gate2.evolve_matrix(three_rho_X)), 5), 1)
    
    def test_cx(self):
        gate = CXgate([0,1])
        Z0Z0 = np.tensordot(self.rho_Z_0, self.rho_Z_0, axes=0)
        Z0Z1 = np.tensordot(self.rho_Z_0, self.rho_Z_1, axes=0)
        Z1Z0 = np.tensordot(self.rho_Z_1, self.rho_Z_0, axes=0)
        Z1Z1 = np.tensordot(self.rho_Z_1, self.rho_Z_1, axes=0)

        # print(Z0Z0)
        # print(gate.evolve_matrix(Z0Z0))

        self.assertEqual(np.allclose(Z0Z0, gate.evolve_matrix(Z0Z0)), 1)
        self.assertEqual(np.allclose(Z1Z1, gate.evolve_matrix(Z1Z0)), 1)
        self.assertEqual(np.allclose(Z0Z1, gate.evolve_matrix(Z0Z1)), 1)
        self.assertEqual(np.allclose(Z1Z0, gate.evolve_matrix(Z1Z1)), 1)
    
    def text_xcx_2qubit(self):
        qc = NoisyCircuit(2)
        qc.x([0])
        qc.cx([0,1])
        Z0Z0 = np.tensordot(self.rho_Z_0, self.rho_Z_0, axes=0)
        Z0Z1 = np.tensordot(self.rho_Z_0, self.rho_Z_1, axes=0)
        Z1Z0 = np.tensordot(self.rho_Z_1, self.rho_Z_0, axes=0)
        Z1Z1 = np.tensordot(self.rho_Z_1, self.rho_Z_1, axes=0)
        self.assertEqual(np.allclose(Z1Z1, qc.matrix_run(Z0Z0)), 1)
        qc.x([1])
        self.assertEqual(np.allclose(Z1Z0, qc.matrix_run(Z0Z0)), 1)
        self.assertEqual(np.allclose(Z0Z0, qc.matrix_run(Z1Z1)), 1)
        self.assertEqual(np.allclose(Z1Z1Z0, qc.matrix_run(Z0Z0Z0)), 1)     
        self.assertEqual(np.allclose(Z1Z0Z0, qc.matrix_run(Z0Z0Z0)), 1)
        self.assertEqual(np.allclose(Z0Z0Z0, qc.matrix_run(Z1Z1Z0)), 1)


    def test_xcx_3qubit(self):
        qc = NoisyCircuit(3)
        qc.x([1])
        qc.cx([0,2])

        Z0Z0Z0 = np.tensordot(np.tensordot(self.rho_Z_0, self.rho_Z_0, axes=0), self.rho_Z_0, axes=0) 
        Z0Z1Z0 = np.tensordot(np.tensordot(self.rho_Z_0, self.rho_Z_1, axes=0), self.rho_Z_0, axes=0)
        Z1Z0Z0 = np.tensordot(np.tensordot(self.rho_Z_1, self.rho_Z_0, axes=0), self.rho_Z_0, axes=0)
        Z1Z1Z0 = np.tensordot(np.tensordot(self.rho_Z_1, self.rho_Z_1, axes=0), self.rho_Z_0, axes=0)
        Z0Z0Z1 = np.tensordot(np.tensordot(self.rho_Z_0, self.rho_Z_0, axes=0), self.rho_Z_1, axes=0) 
        Z0Z1Z1 = np.tensordot(np.tensordot(self.rho_Z_0, self.rho_Z_1, axes=0), self.rho_Z_1, axes=0)
        Z1Z0Z1 = np.tensordot(np.tensordot(self.rho_Z_1, self.rho_Z_0, axes=0), self.rho_Z_1, axes=0)
        Z1Z1Z1 = np.tensordot(np.tensordot(self.rho_Z_1, self.rho_Z_1, axes=0), self.rho_Z_1, axes=0)

        self.assertEqual(np.allclose(Z0Z1Z0, qc.matrix_run(Z0Z0Z0)), 1)     
        self.assertEqual(np.allclose(Z1Z1Z1, qc.matrix_run(Z1Z0Z0)), 1)
        self.assertEqual(np.allclose(Z1Z0Z0, qc.matrix_run(Z1Z1Z1)), 1)

    def test_one_noisy_gate(self):
        D = Depolarization([0], 1)
        mixed = np.array([[0.5, 0],[0, 0.5]])
        ev = D.evolve_matrices([self.rho_X_0])[0]
        ev = ev/ev[0][0]*abs(ev[0][0])
        print(ev)
        print(mixed)
        self.assertEqual(np.allclose(mixed, ev), 1)

    def test_u3(self):
        qc1 = NoisyCircuit(2)
        qc1.x([0])
        qc1.x([1])
        qc2 = NoisyCircuit(2)
        params = [np.pi, 0, np.pi]
        qc2.u3([0], [np.pi, 0, np.pi])
        qc2.u3([1], [np.pi, 0, np.pi])
        rho = np.zeros([2**qc1.num_qubits, 2**qc1.num_qubits])
        rho[0][0] = 1
        rho1 = qc1.matrix_run(rho)
        rho2 = qc2.matrix_run(rho)
        print(rho1)
        print(rho2)


class rxd_gate_almost_ideal(unittest.TestCase):
    def setUp(self):
        self.psi_Z_0 = np.array([1, 0])
        self.psi_Z_1 = np.array([0, 1])
        self.psi_Y_0 = np.array([1, 1j])/np.sqrt(2)
        self.psi_Y_1 = np.array([1, -1j])/np.sqrt(2)
        self.psi_X_0 = np.array([1, 1])/np.sqrt(2)
        self.psi_X_1 = np.array([1, -1])/np.sqrt(2)
        
    def test_one_qubit_Z(self):
        self.qc = QuantumCircuit(1)
        self.qc.rxd([0], np.pi, 0.000001)
        self.assertEqual(np.round(abs(self.qc.run(self.psi_Z_0)[0].conj() @ self.psi_Z_1), 5), 1)
        self.assertEqual(np.round(abs(self.qc.run(self.psi_Z_1)[0].conj() @ self.psi_Z_0), 5), 1)
        self.assertEqual(np.round(abs(self.qc.run(self.psi_X_0)[0].conj() @ self.psi_X_0), 5), 1)
        self.assertEqual(np.round(abs(self.qc.run(self.psi_X_1)[0].conj() @ self.psi_X_1), 5), 1)
        self.assertEqual(np.round(abs(self.qc.run(self.psi_Y_1)[0].conj() @ self.psi_Y_0), 5), 1)
        self.assertEqual(np.round(abs(self.qc.run(self.psi_Y_0)[0].conj() @ self.psi_Y_1), 5), 1)

    def test_three_qubit_cirq(self):

        self.qc = QuantumCircuit(3)
        self.qc.rxd([0, 1, 2], np.pi, 0.000001)
        test_psi = prod3(self.psi_Z_0, self.psi_Y_0, self.psi_X_0)
        xxx_test_psi = prod3(self.psi_Z_1, self.psi_Y_1, self.psi_X_0).reshape(8)
        evolved_psi = self.qc.run(test_psi, 1)[0].reshape(8)
        self.assertEqual(np.round(abs(evolved_psi.conj() @ xxx_test_psi), 5), 1)

        self.qc = QuantumCircuit(3)
        self.qc.rxd([0, 2], np.pi, 0.000001)
        test_psi = prod3(self.psi_Z_0, self.psi_Y_0, self.psi_Y_0)
        xxx_test_psi = prod3(self.psi_Z_1, self.psi_Y_0, self.psi_Y_1).reshape(8)
        evolved_psi = self.qc.run(test_psi, 1)[0].reshape(8)
        self.assertEqual(np.round(abs(evolved_psi.conj() @ xxx_test_psi), 5), 1)

    def two_layer(self):
        self.qc = QuantumCircuit(5)
        self.qc.rxd([0, 1, 2], np.pi, 0.000001)
        self.qc.rxd([0, 1, 2], np.pi, 0.000001)

        test_psi = prod5(self.psi_Z_1, self.psi_Y_0, self.psi_Y_1, self.psi_Z_0, self.psi_X_0)
        xxx_test_psi = test_psi.reshape(32)
        evolved_psi = self.qc.run(test_psi, 1)[0].reshape(32)
        self.assertEqual(np.round(abs(evolved_psi.conj() @ xxx_test_psi), 5), 1)


class test_CZ_gate(unittest.TestCase):
    def setUp(self):
        self.psi_Z_0 = np.array([1, 0])
        self.psi_Z_1 = np.array([0, 1])
        self.psi_Y_0 = np.array([1, 1j])/np.sqrt(2)
        self.psi_Y_1 = np.array([1, -1j])/np.sqrt(2)
        self.psi_X_0 = np.array([1, 1])/np.sqrt(2)
        self.psi_X_1 = np.array([1, -1])/np.sqrt(2)
        
    def test_CZ(self):
        self.qc = QuantumCircuit(3)
        self.qc.cz(0, 1)
        
        test_psi = prod3(self.psi_Z_1, self.psi_X_1, self.psi_Y_0)
        xxx_test_psi = prod3(self.psi_Z_1, self.psi_X_0, self.psi_Y_0).reshape(8)
        evolved_psi = self.qc.run(test_psi, 1)[0].reshape(8)        
        self.assertEqual(np.round(abs(evolved_psi.conj() @ xxx_test_psi), 5), 1)

        test_psi = prod3(self.psi_Z_0, self.psi_X_1, self.psi_Y_0)
        xxx_test_psi = prod3(self.psi_Z_0, self.psi_X_1, self.psi_Y_0).reshape(8)
        evolved_psi = self.qc.run(test_psi, 1)[0].reshape(8)
        self.assertEqual(np.round(abs(evolved_psi.conj() @ xxx_test_psi), 5), 1)
       
        self.qc.cz(0, 2)
        test_psi = prod3(self.psi_Z_1, self.psi_X_1, self.psi_Y_0)
        xxx_test_psi = prod3(self.psi_Z_1, self.psi_X_0, self.psi_Y_1).reshape(8)
        evolved_psi = self.qc.run(test_psi, 1)[0].reshape(8)
        self.assertEqual(np.round(abs(evolved_psi.conj() @ xxx_test_psi), 5), 1)

class test_Noise_circuit(unittest.TestCase):
    def setUp(self):
        self.rho_Z_0 = np.tensordot(np.array([1, 0]), np.array([1, 0]), axes=0)
        a = np.array([0, 1])
        self.rho_Z_0 = np.tensordot(a, a, axes=0)
        a = np.array([1, 1j])/np.sqrt(2)
        self.phi_Y_0 = np.tensordot(a, a, axes=0)
        a = np.array([1, -1j])/np.sqrt(2)
        self.phi_Y_1 = np.tensordot(a, a, axes=0)
        a = np.array([1, 1])/np.sqrt(2)
        self.phi_X_0 = np.tensordot(a, a, axes=0)
        a = np.array([1, -1])/np.sqrt(2)
        self.phi_X_1 = np.tensordot(a, a, axes=0)

    # def test_matrix_evolution(self):
    #     self.qc = NoisyCircuit(2)
    #     self.qc.cz(0,1)
    #     evolved_rho_Z_0 = self.qc.run(self.rho_Z_0, 1)[0]
if __name__ == "__main__":
    unittest.main()