#
# Toolchain file for
#
# GNU compiler
#
# Notes:
#
#  * Settings here should work out of the box on Ubuntu (tested on 18.4). Other build environments
#    may need some fine tuning.
#
#  * CMake format: Command line options (e.g. compiler flags) space separated, other kind
#    of lists semicolon separated.
#
#  * Variables containing library search paths are empty by default. The CMAKE_PREFIX_PATH
#    environment variable should be set up correctly, so that CMake can find those libraries
#    automatically. If that is not the case, override those variables to add search paths
#    manually
#


#
# Fortran compiler settings
#
set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
  CACHE STRING "Build type independent Fortran compiler flags")

set(Fortran_FLAGS_RELEASE "-O2 -funroll-all-loops"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_RELWITHDEBINFO "-g ${Fortran_FLAGS_RELEASE}"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_DEBUG "-g -Wall -std=f2008ts -pedantic -fbounds-check"
  CACHE STRING "Fortran compiler flags for Debug build")

# Use intrinsic Fortran 2008 erf/erfc functions
set(INTERNAL_ERFC CACHE BOOL 0)

set(FYPP_FLAGS "" CACHE STRING "Fypp preprocessor flags")


#
# C compiler settings
#
set(C_FLAGS "${CMAKE_C_FLAGS}"
  CACHE STRING "Build type independent C compiler flags")

set(C_FLAGS_RELEASE "-O2 -funroll-all-loops"
  CACHE STRING  "C compiler flags for Release build")

set(C_FLAGS_RELWITDEBINFO "-g ${C_FLAGS_RELEASE}"
  CACHE STRING  "C compiler flags for RelWithDebInfo build")

set(C_FLAGS_DEBUG "-g -Wall -pedantic -fbounds-check"
  CACHE STRING "C compiler flags for Debug build")


#
# External libraries
#

# NOTE: CMake searches for the external libraries in the standard system library paths and in the
# paths defined in the CMAKE_PREFIX_PATH environment variable. Make sure, CMAKE_PREFIX_PATH contains
# the path to your libraries in case those are not installed in the standard locations (e.g. when
# using custom installed libraries or environment modules in HPC centers).

# LAPACK and BLAS
#set(LAPACK_LIBRARY "openblas" CACHE STRING "LAPACK and BLAS libraries to link")

# ARPACK -- only needed when built with ARPACK support
#set(ARPACK_LIBRARY "arpack" CACHE STRING "Arpack libraries")

# ScaLAPACK -- only needed for MPI-parallel build
#set(SCALAPACK_LIBRARY "scalapack-openmpi" CACHE STRING "Scalapack libraries to link")

# PLUMED -- only needed when compiled with PLUMED support
#set(PLUMED_LIBRARY "plumed;plumedKernel" CACHE STRING "Libraries to link for PLUMED support")

# MAGMA -- only needed when compiled with GPU support
#set(MAGMA_LIBRARY "magma" CACHE STRING "Magma library")
#set(MAGMA_INCLUDE_DIRECTORY "" CACHE STRING "Directories to scan for MAGMA include files")
