from abc import ABC, abstractmethod
import numpy as np
from numpy.typing import ArrayLike 

from .instructions import Instruction

X = np.array([[0,1],[1,0]])
Y = np.array([[0,-1j],[1j,0]])
Z = np.array([[1,0],[0,-1]])
I = np.eye(2)


class NoisyInstruction(Instruction, ABC):
    def __init__(self, 
                 qubits: ArrayLike,
                 prob:float):
        super().__init__(qubits)
        self.prob = prob
        
    def evolve_matrices(self, rho_array):
        kraus, probs = self.get_kraus_and_probs()
        a = 0 
        ax_down_rho = (self.axes_rho - 1)*2 - self.axes_rho[0] + 2
        for k, p in zip(kraus, probs):
            a = np.moveaxis(np.tensordot(np.moveaxis(np.tensordot(rho_array, k.T, axes=(ax_down_rho, self.axes_ops - 1)),
                                                                        source=-1*(self.axes_ops), 
                                                                        destination=reversed(ax_down_rho)),
                                        k.conj().T, 
                                        axes=(ax_down_rho + 1, self.axes_ops - 1)),
                            source=-1*(self.axes_ops), 
                            destination=reversed(ax_down_rho + 1))*p + a
        return a
            
    @abstractmethod
    def get_kraus_and_probs():
        pass

    def ops(self, counts):
        pass


class Depolarization(NoisyInstruction):
    def get_kraus_and_probs(self):
        kraus = [I,X,Y,Z]
        probs = [1-3/4*self.prob, self.prob/4, self.prob/4, self.prob/4]
        return kraus, probs

class Xflip(NoisyInstruction):
    def get_kraus_and_probs(self):
        kraus = [I,X]
        probs = [1-self.prob, self.prob]
        return kraus, probs


class Zflip(NoisyInstruction):
    def get_kraus_and_probs(self):
        kraus = [I,Z]
        probs = [1-self.prob, self.prob]
        return kraus, probs


class Yflip(NoisyInstruction):
    def get_kraus_and_probs(self):
        kraus = [I,Y]
        probs = [1-self.prob, self.prob]
        return kraus, probs


class AmpDamp(NoisyInstruction):
    def get_kraus_and_probs(self):
        x = self.prob
        kraus = [np.array([[1,0],[0,x]]), np.array([[0, np.sqrt(1 - x**2)],[0,0]])]
        probs = [1,1]
        return kraus, probs

    @property
    def name(self):
        return 'damp ' + str(1 - self.prob**2) 


class PhaseDamp(NoisyInstruction):
    def get_kraus_and_probs(self):
        x = self.prob
        kraus = [I,Z]
        probs = [x**2,1-x**2]
        return kraus, probs

    @property
    def name(self):
        return 'damp ' + str((1 - self.prob)/2) 