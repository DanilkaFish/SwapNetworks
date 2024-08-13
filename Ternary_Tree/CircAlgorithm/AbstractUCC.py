from __future__ import annotations

from qiskit.circuit import Parameter, QuantumCircuit
from openfermionpyscf import run_pyscf
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
from .utils import simplification, single_exc, double_exc
from abc import ABC, abstractmethod
from copy import deepcopy

class AbstractUCC(ABC):
    def __init__(self,
                 geometry: str= 'H 0 0 0; Li 0 0 0.7414',
                 basis: str = "6-31G",
                 multiplicity: int = 1,
                 charge: int = 0,
                 active_orbitals: list[list[int,int]] = None,
                 num_electrons: int | tuple[int, int] = None
                 ):
        driver = PySCFDriver(atom=geometry,
                              basis=basis,
                              spin= (multiplicity - 1)//2,
                              charge=charge)
        self.mol = driver.run()
        if active_orbitals is not None:
            num_spatial = len(active_orbitals)
            if num_electrons is None:
                num_electrons = self.mol.num_alpha + self.mol.num_beta
            transformer = ActiveSpaceTransformer(num_electrons=num_electrons, 
                                                num_spatial_orbitals=num_spatial, 
                                                active_orbitals=active_orbitals)
            self.mol = transformer.transform(self.mol)

        self.fermionic_op = self.mol.hamiltonian.second_q_op()
        self.n_qubits = int(self.mol.num_spin_orbitals)
        self.maj_exc_par = self.get_excitations()


    @property
    def num_alpha(self):
        return self.mol.num_alpha

    @property
    def num_beta(self):
        return self.mol.num_beta

    @property
    def num_spatial_orbitals(self):
        return self.n_qubits//2

    @property
    def num_spin_orbitals(self):
        return self.n_qubits

    def to_par_maj_excitations(self, ladder_exciations, name="t_" ):
        ladder_exc_par = {}
        for op in ladder_exciations:
            ladder_exc_par[op] = Parameter(name + ','.join([str(i) for i in op]))
            
        maj_exc_par = {}
        for op in ladder_exc_par:
            if len(op) == 2:
                for new_op, sign in single_exc(op):
                    maj_exc_par[new_op] = sign * ladder_exc_par[op] + maj_exc_par.get(new_op, 0)
            if len(op) == 4:
                for new_op, sign in double_exc(op):
                    maj_exc_par[new_op] = sign * ladder_exc_par[op] + maj_exc_par.get(new_op, 0)

        maj_exc_par = simplification(maj_exc_par)
        return maj_exc_par

    def to_par_maj_exitations_comp(self, maj_alpha_par_exc, n, name="t_", ):
        toto = deepcopy(maj_alpha_par_exc)
        for key in maj_alpha_par_exc:
            toto[(key[0] + n, key[1] + n)] = toto[key]
        return toto      

    def to_par_ladder_exciations(self, ladder_exciations, name="t_" ):
        ladder_exc_par = {}
        for op in ladder_exciations:
            ladder_exc_par[op] = Parameter(name + ','.join([str(i) for i in op]))
        return ladder_exc_par

    
    def get_excitations(self, name='t_') -> {(int, int) | (int, int)}:
        """
        Method to get parametrized UpCCGSD excitations via majorana operators
        """
        ladder_exc_par = {}

        ladder_exc = self.get_alpha_excitations()
        ladder_exc += self.get_beta_excitations()
        ladder_exc += self.get_double_excitations()
        for op in ladder_exc:
            ladder_exc_par[op] = Parameter(name + ','.join([str(i) for i in op]))
            
        maj_exc_par = {}
        for op in ladder_exc_par:
            if len(op) == 2:
                for new_op, sign in single_exc(op):
                    maj_exc_par[new_op] = sign * ladder_exc_par[op] + maj_exc_par.get(new_op, 0)
            if len(op) == 4:
                for new_op, sign in double_exc(op):
                    maj_exc_par[new_op] = sign * ladder_exc_par[op] + maj_exc_par.get(new_op, 0)

        maj_exc_par = simplification(maj_exc_par)
        return maj_exc_par

    @abstractmethod
    def get_alpha_excitations(self) -> list[tuple[int, int]]:
        pass

    @abstractmethod
    def get_beta_excitations(self) -> list[tuple[int, int]]:
        pass

    @abstractmethod
    def get_double_excitations(self) -> list[tuple[int, int, int, int]]:
        pass

    @abstractmethod
    def get_parametrized_circuit(self) -> QuantumCircuit:
        pass

