#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

"""
Module for interaction with the results.tag file.
"""

from __future__ import annotations
import dftbplus_ptools.tagreader as reader
import numpy as np

class Output:
    """Class for reading dftb+ outputs"""

    def __init__(self, filename: str = 'results.tag')  -> None:
        """Initialises the Output class

        Args:
            filename (str): name of file
        """
        self._filename = filename
        self._resultstag = self._read_resultstag()

    def _read_resultstag(self) -> dict | None:
        """Function for accessing results.tag file."""
        try:
            return reader.results_access(filename=self._filename)
        except OSError:
            return None


    def get_cell_volume(self) -> np.ndarray | None:
        """Help function for reading 'unit cell volume' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['cell_volume']
            except KeyError:
                return None


    def get_end_coords(self) -> np.ndarray | None:
        """Help function for reading 'final geometry' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['end_coords']
            except KeyError:
                return None


    def get_exc_energies_sqr(self) -> np.ndarray | None:
        """Help function for reading 'excitation energies in Casida formalism'
        from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['exc_energies_sqr']
            except KeyError:
                return None


    def get_exc_forces(self) -> np.ndarray | None:
        """Help function for reading 'excited state force contributions'
        from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['exc_forces']
            except KeyError:
                return None


    def get_exc_oscillator(self) -> np.ndarray | None:
        """Help function for reading 'oscillator strength for excitations'
        from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['exc_oscillator']
            except KeyError:
                return None


    def get_exc_transdip(self) -> np.ndarray | None:
        """Help function for reading 'Transition dipole moments for
        excitations' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['exc_transdip']
            except KeyError:
                return None


    def get_coupling_vectors(self) -> np.ndarray | None:
        """Help function for reading 'nonadiabatic coupling vector, H' from
        results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['coupling_vectors']
            except KeyError:
                return None


    def get_forces(self) -> np.ndarray | None:
        """Help function for reading 'ground state total forces' from
        results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['forces']
            except KeyError:
                return None


    def get_forces_ext_charges(self) -> np.ndarray | None:
        """Help function for reading 'forces on any external charges present'
        from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['forces_ext_charges']
            except KeyError:
                return None


    def get_fermi_level(self) -> np.ndarray | None:
        """Help function for reading 'Fermi level(s)' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['fermi_level']
            except KeyError:
                return None


    def get_number_of_electrons(self) -> np.ndarray | None:
        """Help function for reading 'number of electrons' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['number_of_electrons']
            except KeyError:
                return None


    def get_eigenvalues(self) -> np.ndarray | None:
        """Help function for reading 'eigenvalues/single particle states' from
        results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['eigenvalues']
            except KeyError:
                return None


    def get_filling(self) -> np.ndarray | None:
        """Help function for reading 'filling of the eigenstates' from
        results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['filling']
            except KeyError:
                return None


    def get_gibbs_energy(self) -> np.ndarray | None:
        """Help function for reading 'Gibbs free energy for finite pressure
        periodic systems' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['gibbs_energy']
            except KeyError:
                return None


    def get_gross_atomic_charges(self) -> np.ndarray | None:
        """Help function for reading 'Gross atomic charges' from
        results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['gross_atomic_charges']
            except KeyError:
                return None


    def get_cm5_atomic_charges(self) -> np.ndarray | None:
        """Help function for reading 'Charge model 5 corrected atomic gross
        charges' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['cm5_atomic_charges']
            except KeyError:
                return None


    def get_gross_atomic_spins(self) -> np.ndarray | None:
        """Help function for reading 'Gross atomic spin polarizations' from
        results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['gross_atomic_spins']
            except KeyError:
                return None


    def get_hessian_numerical(self) -> np.ndarray | None:
        """Help function for reading 'numerically calculated second
        derivatives matrix' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['hessian_numerical']
            except KeyError:
                return None


    def get_final_energy(self) -> np.ndarray | None:
        """Help function for reading 'final energy components after real-time
        propagation' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['final_energy']
            except KeyError:
                return None


    def get_final_dipole_moment(self) -> np.ndarray | None:
        """Help function for reading 'final dipole moment vector after
        real-time propagation' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['final_dipole_moment']
            except KeyError:
                return None


    def get_final_td_charges(self) -> np.ndarray | None:
        """Help function for reading 'final negative gross atomic Mulliken
        charges after real-time propagation' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['final_td_charges']
            except KeyError:
                return None


    def get_final_ehrenfest_forc(self) -> np.ndarray | None:
        """Help function for reading 'final forces components after real-time
        (Ehrenfest) propagation' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['final_ehrenfest_forc']
            except KeyError:
                return None


    def get_final_ehrenfest_geom(self) -> np.ndarray | None:
        """Help function for reading 'final geometry after real-time
        (Ehrenfest) propagation' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['final_ehrenfest_geom']
            except KeyError:
                return None


    def get_final_ehrenfest_velo(self) -> np.ndarray | None:
        """Help function for reading 'final velocities after real-time
        (Ehrenfest) propagation' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['final_ehrenfest_velo']
            except KeyError:
                return None


    def get_final_td_proj_occ(self) -> np.ndarray | None:
        """Help function for reading 'final molecular orbitals occupations
        after real-time (Ehrenfest) propagation' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['final_td_proj_occ']
            except KeyError:
                return None


    def get_sum_bond_pops(self) -> np.ndarray | None:
        """Help function for reading 'Sum of bond populaion values (should be
        number of electrons)' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['sum_bond_pops']
            except KeyError:
                return None


    def get_mermin_energy(self) -> np.ndarray | None:
        """Help function for reading 'total energy including electron TS
        contribution' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['mermin_energy']
            except KeyError:
                return None


    def get_orbital_charges(self) -> np.ndarray | None:
        """Help function for reading 'Mulliken charges' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['orbital_charges']
            except KeyError:
                return None


    def get_pm_localisation(self) -> np.ndarray | None:
        """Help function for reading 'Pipek-Mezey localisation score of single
        particle levels' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['pm_localisation']
            except KeyError:
                return None


    def get_stress(self) -> np.ndarray | None:
        """Help function for reading 'total stress tensor for periodic
        geometries' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['stress']
            except KeyError:
                return None


    def get_total_tunneling(self) -> np.ndarray | None:
        """Help function for reading 'total tunneling vector' from
        results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['total_tunneling']
            except KeyError:
                return None


    def get_total_localdos(self) -> np.ndarray | None:
        """Help function for reading 'total projected DOS vector' from
        results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['total_localdos']
            except KeyError:
                return None


    def get_local_currents(self) -> np.ndarray | None:
        """Help function for reading 'total bond currents' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['local_currents']
            except KeyError:
                return None


    def get_total_energy(self) -> np.ndarray | None:
        """Help function for reading 'total internal energy' from
        results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['total_energy']
            except KeyError:
                return None


    def get_extrapolated0_energy(self) -> np.ndarray | None:
        """Help function for reading 'total internal energy extrapolated to
        0 K' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['extrapolated0_energy']
            except KeyError:
                return None


    def get_forcerelated_energy(self) -> np.ndarray | None:
        """Help function for reading 'Energy, which if differentiated gives -
        force' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['forcerelated_energy']
            except KeyError:
                return None


    def get_internal_efield(self) -> np.ndarray | None:
        """Help function for reading 'Internal electric field' from
        results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['internal_efield']
            except KeyError:
                return None


    def get_external_efield(self) -> np.ndarray | None:
        """Help function for reading 'External electric field' from
        results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['external_efield']
            except KeyError:
                return None


    def get_staticPolResponse(self) -> np.ndarray | None:
        """Help function for reading 'Static electric polarizability from
        linear response/perturbation' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['staticPolResponse']
            except KeyError:
                return None


    def get_staticChargeReponse(self) -> np.ndarray | None:
        """Help function for reading 'Static gross charge (Mulliken) response
        from linear response/perturbation' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['staticChargeReponse']
            except KeyError:
                return None


    def get_dEidEfield(self) -> np.ndarray | None:
        """Help function for reading 'Derivatives of ground state single
        particle eigenvalues wrt. k' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['dEidEfield']
            except KeyError:
                return None


    def get_neFermi(self) -> np.ndarray | None:
        """Help function for reading 'Number of electrons at the Fermi energy'
        from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['neFermi']
            except KeyError:
                return None


    def get_dEfdE(self) -> np.ndarray | None:
        """Help function for reading 'Derivative of the Fermi energy with
        respect to electric field' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['dEfdE']
            except KeyError:
                return None


    def get_dEidVons(self) -> np.ndarray | None:
        """Help function for reading 'Derivatives of ground state single
        particle eigenvalues wrt. onsite potentials' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['dEidVons']
            except KeyError:
                return None


    def get_dEidV(self) -> np.ndarray | None:
        """Help function for reading 'Derivatives of ground state single
        particle eigenvalues wrt. potential at an atom' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['dEidV']
            except KeyError:
                return None


    def get_dqdV(self) -> np.ndarray | None:
        """Help function for reading 'Static gross charge (Mulliken) response
        with respect to potential at an atom' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['dqdV']
            except KeyError:
                return None


    def get_dqnetdV(self) -> np.ndarray | None:
        """Help function for reading 'Static net charge (onsite) response with
        respect to potential at an atom' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['dqnetdV']
            except KeyError:
                return None


    def get_2e_add_rem_energies(self) -> np.ndarray | None:
        """Help function for reading 'two-electron addition/removal energies
        in ppRPA formalism' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['2e_add-rem_energies']
            except KeyError:
                return None


    def get_atomic_masses(self) -> np.ndarray | None:
        """Help function for reading 'atomic masses' from results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['atomic_masses']
            except KeyError:
                return None


    def get_atomic_dipole_moments(self) -> np.ndarray | None:
        """Help function for reading 'Atomic dipole moments' from
        results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['atomic_dipole_moments']
            except KeyError:
                return None


    def get_dipole_moments(self) -> np.ndarray | None:
        """Help function for reading 'Total dipole moment' from
        results.tag"""
        if self._resultstag is not None:
            try:
                return self._resultstag['dipole_moments']
            except KeyError:
                return None
