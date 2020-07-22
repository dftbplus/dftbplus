**********
Change Log
**********

Notable project changes since release 1.3.1 (2017-02-22).


Unreleased
==========

Added
-----


Changed
-------


Fixed
-----



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
