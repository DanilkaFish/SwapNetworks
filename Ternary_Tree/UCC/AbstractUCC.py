from __future__ import annotations

from qiskit.circuit import Parameter, QuantumCircuit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
from .utils import lad2maj
from abc import ABC, abstractmethod
from typing import List, Tuple, Dict, Union, Optional
from .excitation import LadExcitation, DoubleLadExcitation, SingleLadExcitation


class Molecule:
    def __init__(self,
                 geometry='H 0 0 0; Li 0 0 0.7414',
                 basis: str="6-31G",
                 multiplicity: int=1,
                 charge: int=0,
                 active_orbitals: Optional[List[Tuple[int,int]]]=None,
                 num_electrons: Optional[Union[int, Tuple[int, int]]]=None
                 ):
        self.geometry: str = geometry
        self.basis: str = basis
        self.multiplicity: int = multiplicity
        self.charge: int = charge
        self.active_orbitals: Optional[List[Tuple[int,int]]] = active_orbitals
        self.num_electrons: Optional[Union[int, Tuple[int, int]]] = num_electrons

    def get():
        return 1
    
class AbstractUCC(ABC):
    def __init__(self, molecule: Molecule=Molecule()):
        driver = PySCFDriver(atom=molecule.geometry,
                              basis=molecule.basis,
                              spin=(molecule.multiplicity - 1)//2,
                              charge=molecule.charge)

        self.mol = driver.run()
        if molecule.active_orbitals is not None:
            if molecule.num_electrons is None:
                molecule.um_electrons = self.n_alpha + self.n_beta
            transformer = ActiveSpaceTransformer(num_electrons=molecule.num_electrons, 
                                                 num_spatial_orbitals=len(molecule.active_orbitals), 
                                                 active_orbitals=molecule.active_orbitals)
            self.mol = transformer.transform(self.mol)

        self.fermionic_op = self.mol.hamiltonian.second_q_op()
        self.n_qubits = int(self.mol.num_spin_orbitals)
        self.maj_exc_par = self.get_excitations()

    @property
    def n_alpha(self) -> int:
        return self.mol.num_alpha

    @property
    def n_beta(self) -> int:
        return self.mol.num_beta

    @property
    def n_spatial(self) -> int:
        return self.n_qubits//2

    @property
    def n_spin(self) -> int:
        return self.n_qubits
    
    def get_excitations(self, name: str='t_') -> Dict[LadExcitation, Parameter]:
        """
        Method to get parametrized UpCCGSD excitations via majorana operators
        """
        ladder_exc = self.get_alpha_excitations()
        ladder_exc += self.get_beta_excitations()
        ladder_exc += self.get_double_excitations()
        return lad2maj(ladder_exc)

    @abstractmethod
    def get_alpha_excitations(self) -> list[SingleLadExcitation]:
        pass

    @abstractmethod
    def get_beta_excitations(self) -> list[SingleLadExcitation]:
        pass

    @abstractmethod
    def get_double_excitations(self) -> list[DoubleLadExcitation]:
        pass

    @abstractmethod
    def swap2xn(self) -> QuantumCircuit:
        pass

