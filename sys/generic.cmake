#
# Toolchain file for
#
# Generic build environment
#
# This is a generic template which probably will not work on your system out of the box. You should
# either modify it or override the variables via command line options to make it to
# work. Alternatively, have a look at the specialized toolchain files in this folder as they may
# give you a better starting point for your build environment.
#
# Notes:
#
#  * CMake format: Command line options (e.g. compiler flags) space separated, other kind
#    of lists semicolon separated.
#
#  * Variables containing library search paths are empty by default. The CMAKE_PREFIX_PATH
#    environment variable should be set up correctly, so that CMake can find those libraries
#    automatically. If that is not the case, override those variables to add search paths
#    manually


#
# Fortran compiler settings
#
set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
  CACHE STRING "Build type independent Fortran compiler flags")

set(Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_RELWITHDEBINFO "${Fortran_FLAGS_RELWITHDEBINFO}"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
  CACHE STRING "Fortran compiler flags for Debug build")

# Use intrinsic Fortran 2008 erf/erfc functions
set(INTERNAL_ERFC CACHE BOOL 0)

set(FYPP_FLAGS ""
  CACHE STRING "Flags for the preprocessor")


#
# C compiler settings
#
set(C_FLAGS "${CMAKE_C_FLAGS}"
  CACHE STRING "Build type independent C compiler flags")

set(C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}"
  CACHE STRING  "C compiler flags for Release build")

set(C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
  CACHE STRING "C compiler flags for Debug build")


#
# External libraries
#

# NOTE: Libraries with CMake export files (e.g. ELSI and if the HYBRID_CONFIG_METHODS variable
# contains the "Find" method also libNEGF, libMBD, ScalapackFx and MpiFx) are included by searching
# for the CMake export file in the paths defined in the CMAKE_PREFIX_PATH **environment**
# variable. Make sure your CMAKE_PREFIX_PATH variable is set up accordingly.

# LAPACK and BLAS
#set(LAPACK_LIBRARY "lapack;blas" CACHE STRING "LAPACK and BLAS libraries to link")
#set(LAPACK_LIBRARY_DIR "" CACHE STRING
#  "Directories where LAPACK and BLAS libraries can be found")

# ARPACK -- only needed when built with ARPACK support
#set(ARPACK_LIBRARY "arpack" CACHE STRING "Arpack libraries")
#set(ARPACK_LIBRARY_DIR "" CACHE STRING "Directories where Arpack library can be found")

# ScaLAPACK -- only needed for MPI-parallel build
#set(SCALAPACK_LIBRARY "scalapack" CACHE STRING "Scalapack libraries to link")
#set(SCALAPACK_LIBRARY_DIR "" CACHE STRING "Directories where Scalapack libraries can be found")

# NOTE: The libraries below provide Pkg-Conf export files.  If your PKG_CONFIG_PATH environment
# variable has been set up correctly (containing the paths to these libraries), no adjustment should
# be necessary below.

# PLUMED -- only needed when compiled with PLUMED support
#set(PLUMED_LIBRARY "plumed;plumedKernel" CACHE STRING "Libraries to link for PLUMED support")
#set(PLUMED_LIBRARY_DIR "" CACHE STRING "Directories to scan for PLUMED libraries")

# MAGMA -- only needed when compiled with GPU support
#set(MAGMA_LIBRARY "magma" CACHE STRING "Magma library")
#set(MAGMA_LIBRARY_DIR "" CACHE STRING "Directories to scan for MAGMA library")
#set(MAGMA_INCLUDE_DIRECTORY "" CACHE STRING "Directories to scan for MAGMA include files")
