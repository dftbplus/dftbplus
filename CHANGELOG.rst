**********
Change Log
**********

Notable project changes since release 1.3.1 (2017-02-22).


Unreleased
==========


Fixed
-----

- Onsite and +U potentials in real time-propagation, which was broken
  in October 2019 by commit 11abba39b


22.1 (2022-05-25)
=================

Added
-----

- Real time electronic dynamics for xTB Hamiltonian

- Real time electronic dynamics for range separated DFTB

- Support for MPI-parallel GPU accelerated calculations via ELPA/ELSI library

- (Optionally) rescale externally applied fields and dipole moments
  when implicit solvents are used

- Enable lattice constraints in new geometry optimization driver

- Dynamic polarizability and response kernel at finite frequencies

- API call for CM5 charges

- Numerical Hessian calculation can be split over multiple runs


Changed
-------

- PLUMED simulations may deliver due to an incompatible change in version 2.8.0
  of the external PLUMED library slightly different results as before. See also
  the `change log of PLUMED 2.8
  <https://www.plumed.org/doc-v2.8/user-doc/html/_c_h_a_n_g_e_s-2-8.html>`_.

- Allow electric fields in periodic systems even when interactions
  cross the sawtooth in the field

- Allow printing of dipole moments, even in cases where the absolute
  value is ill-defined (charged systems or periodic cases), but its
  derivative may be meaningful.

- Use the DFTB+ xyz writer for the modes program, removing the
  XMakemol output option.

- Re-enable q=0 (sawtooth) electric fields for periodic/helical structures


Fixed
-----

- incorrect atomic mass unit for xTB calculations

- electronic temperature read for Green's function solver

- MPI code for spin polarised metallic perturbation at q=0 for spin
  polarized molecules with processor groups


21.2 (2021-12-13)
=================

Added
-----

- On-site potentials added

- Support for extended tight binding (xTB) Hamiltonian via tblite library

- DFTBPLUS_PARAM_DIR for searching Slater-Koster parameter files, solvation
  parameter files, and xTB parameter files

- Atomic potential responses (enables atom resolved response kernel evaluation
  and condensed Fukui functions)

- Internal changes for response evaluation for DFTB ground state hamiltonians
  (except self-consistent dispersion) with molecular, periodic and helical
  boundary conditions.

- Stratmann solver for excited state, including range separated calculations

- Rational function geometry optimization driver

- ChIMES force field corrections of the repulsive potentials implemented

- New geometry optimization drivers with coupled cartesian and lattice parameter
  optimization


Changed
-------

- Source tree reorganised to match the `Fortran package manager
  <https://fpm.fortran-lang.org/>`_ preferred structure.

- Updated parser version to 10.

- Replace backend to implement DFT-D3 dispersion correction.
  Use `s-dftd3 <https://github.com/awvwgk/simple-dftd3>`_ instead of
  `dftd3-lib <https://github.com/dftbplus/dftd3-lib>`_.
  Option ``WITH_DFTD3`` is removed and replaced with ``WITH_SDFTD3``.


Fixed
-----

- CM5 correction added with incorrect sign to charge populations

- External fields disabled for XLBOMD

- self-consistent DFT-D4 uses populations instead of partial charges
  in potential shift, energy expression and derivatives

- Number of electrons for Fixed / spin-common Fermi energies and transport in
  results.tag

- D3(BJ)-ATM calculator was not being passed the exponent for ATM zero damping
  calculations

- LBFGS implementation fixed in new geometry optimization driver


21.1 (2021-05-12)
=================

Added
-----

- Conductor like screening model (COSMO) implicit solvation model for SCC
  calculations

- Printout of cavity information as a cosmo file

- Extended syntax for selecting atoms in HSD input

- Static coupled perturbed response for homogeneous electric fields (evaluating
  molecular electric polarisability)


Changed
-------

- DFT-D4 can now be evaluated self-consistently within the SCC procedure

- Self-consistent DFT-D4 with REKS

- Upgraded to libMBD 0.12.1 (TS-forces are calculated analytically)


Fixed
-----

- Fix bug in binary eigenvector output in non-MPI builds (only eigenvectors
  belonging to the first k-point and spin channel were stored)

- Fix transpose of lattice vectors on return from iPI (thanks to Bingqing Cheng
  and Edgar Engel)


20.2.1 (2020-12-07)
===================

Fixed
-----

- Lattice derivatives are now correctly written into detailed.out

- Upgraded to libNEGF version 1.0.1 fixing usage of uninitialized variables

- Removed '-heap-arrays' option from ifort compiler options to work around Intel
  compiler bug causing steadily increasing memory consumption during long runs


20.2 (2020-11-17)
=================

Added
-----

- Many body and Tkatchenko-Scheffler dispersion

- Delta DFTB for lowest singlet excitated state

- Electron transport for system with colinear spin polarisation

- Phonon transport calculations with new code

- Linear response gradients for spin polarisation

- FIRE geometry optimizer

- Simple D3-dispersion implementation (can be used without needing the external
  D3-library)


Changed
-------

- MPI parallelisation for UFF, Slater-Kirkwood and DFT-D4 dispersion

- OMP parallelisation for UFF and Slater-Kirkwood dispersion

- Option to take quasi-Newton steps in lBFGS (set as default)

- CMake cache variable names in accordance with CMake devel documentation


Fixed
-----

- Stress tensor is now calculated with Slater-Kirkwood dispersion

- Cube format closer to the files expected by several external tools


20.1 (2020-07-22)
=================

Added
-----

- REKS (spin-Restricted Ensemble Kohn-Sham) calculations for ground and
  low-lying exited states

- Support for meta-dynamics in MD via the Plumed library

- Option to set mass of atoms in the modes code input file (syntax matches
  existing DFTB+ feature)

- Use of processor groups with transport calculations, enabling better
  parallelism for systems that need k-points

- Reading of input coordinates in XYZ format

- Reading of input coordinates in the VASP POSCAR format

- The DFT-D4 dispersion model

- Helical geometries supported for non-SCC calculations

- Generalised Born (GB) and Analytical Linearised Poisson-Boltzmann (ALPB)
  implicit solvation models for SCC calculations

- Non-polar solvent accessible surface area solvation model

- Particle-particle random-phase approximation available for suitable excitation
  calculations

- Range separated excited state calculations for spin free singlet systems

- New algorithm for the ground state range-separated hamiltonian

- Real time electronic and coupled electron-ion Ehrenfest dynamics


Changed
-------

- New build system using CMake (the old makefile system has been retired)

- Input in GEN format now strictly follows the description in the manual

- Versioned format for transport contact shift files (backward compatible), also
  enables the Fermi energy to be read directly from the contact file.

- Removed residual XML input (leaving detailed.xml export, depreciating the
  undocumented <<! tag in HSD)

- Output of energies clarified (total energy when electron entropy is not
  available, Mermin free energy when it is and force related energy when the
  energy associated with Helmann-Feynman forces is available)

- API extended for MPI parallel calculations and interfaces added to obtain API
  version and DFTB+ release.

- Poisson solver available without libNEGF enabled compilation

- Parser input can now be set according to the code release version (20.1)


Fixed
-----

- Correct update of block Mulliken population for onsite correction with
  range-separation hybrid DFTB.

- MD temperature profiles that do not start with an initial constant temperature

- Free energy for PEXSI calculations

- ELSI calculations for spin-orbit and onsite corrected corrections


19.1 (2019-07-01)
=================

Added
-----

- Non-equilibrium Green's function transport.

- Use of the ELSI library.

- Ability to perform ground state MD with excitation energies.

- Caching for transition charges in excited state.

- DFTB+ can be compiled as a library and accessed via high level API (version
  number is in the file api/mm/API_VERSION below the main directory).

- Onsite corrected hamiltonian for ground state energies.

- Range-separated hybrid DFTB.

- GPU acceleration using the MAGMA library for eigensolution. WARNING: this is
  currently an experimental feature, so should be used with care.

- Labelling of atomic orbital choices in output.

- Halogen X correction.


Changed
-------

- Updated parser version to 7.


Fixed
-----

- Orbital-resolved projected eigenstates (shell-resolved were correct)

- Corrected Orbital to Shell naming conventions


18.2 (2018-08-19)
=================

Added
-----

- Option for removing translational and rotational degrees of freedom in modes.

- H5 correction for hydrogen bonds.


Changed
-------

- Updated parser version to 6.

- Syntax for H5 and DampedHX corrections for hydrogen bonds unified.


Fixed
-----

- Compilation when socket interface disabled.

- Stress tensor evaluation for 3rd order DFTB.

- Tollerance keyword typo.

- Corrected erroneous Lennard-Jones-dispersion for periodic cases (broken since
  release 1.3)

- Forces/stresses for dual spin orbit.


18.1 (2018-03-02)
=================

Added
-----

- MPI-parallelism.

- Various user settings for MPI-parallelism.

- Improved thread-parallelism.

- LBGFS geometry driver.

- Evaluation of electrostatic potentials at specified points in space.

- Blurred external charges for periodic systems.

- Option to read/write restart charges as ASCII text.

- Timer for collecting timings and printing them at program end.

- Tolerance of Ewald summation can be set in user input.

- Shutdown possibility when using socket driver.

- Header for code prints release / git commit version information.

- Warning when downloading license incompatible external components.

- Tool straingen for distorting gen-files.


Changed
-------

- Using allocatables instead of pointers where possible.

- Change to use the Fypp-preprocessor.

- Excited state (non-force) properties for multiple excitations.

- Broyden-mixer does not use file I/O.

- Source code documentation is Ford-compatible.

- Various refactorings to improve on modularity and code clarity.


Fixed
-----

- Keyword Atoms in modes_in.hsd consider only the first specified entry.

- Excited window selection in Cassida time-dependent calculation.

- Formatting of eigenvalues and fillings in detailed.out and band.out

- iPI socket interface with cluster geometries fixed (protocol contains
  redundant lattice information in these cases).


17.1 (2017-06-16)
=================

Added
-----

- Add dptools toolkit.


Changed
-------

- Convert to LGPL 3 license.

- Restructure source tree.

- Streamline autotest suite and build system.


Fixed
-----

- Skip irrelevant tests that give false positives for particular compilation
  modes.

- Make geometry writing in gen and xyz files consistent.
