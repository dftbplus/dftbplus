**********
Change Log
**********

Notable project changes since release 1.3.1 (2017-02-22).


Unreleased
=========

Added
-----

- MPI-parallelism.

- Improved thread-parallelism.

- Tool straingen for distorting gen-files.

- Shutdown possibility when using socket driver.

- Evaluation of electrostatic potentials at specified points in space.

- Blurred external charges for periodic systems.

- Option to read/write restart charges as ASCII text.

- Header for code prints release / git commit version information.

Changed
-------

- Using allocatables instead of pointers whenever possible.

- Change to use the Fypp-preprocessor.

- Excited state (non-force) properties for multiple excitations.

Fixed
-----

- Keyword Atoms in modes_in.hsd consider only the first specified entry.

- Excited window selection in Cassida time-dependent calculation.

- iPI interface with cluster geometries fixed (protocol contains redundant
  lattice information in these cases).

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
