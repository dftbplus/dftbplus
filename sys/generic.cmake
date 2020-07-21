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

set(Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
  CACHE STRING "Fortran compiler flags for Debug build")

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

# LAPACK and BLAS
set(LAPACK_LIBRARIES "lapack;blas" CACHE STRING "LAPACK and BLAS libraries to link")
set(LAPACK_LIBRARY_DIRS "" CACHE STRING
  "Directories where LAPACK and BLAS libraries can be found")

# ARPACK -- only needed when built with ARPACK support
set(ARPACK_LIBRARIES "arpack" CACHE STRING "Arpack libraries")
set(ARPACK_LIBRARY_DIRS "" CACHE STRING "Directories where Arpack library can be found")

# ScaLAPACK -- only needed for MPI-parallel build
set(SCALAPACK_LIBRARIES "scalapack" CACHE STRING "Scalapack libraries to link")
set(SCALAPACK_LIBRARY_DIRS "" CACHE STRING "Directories where Scalapack libraries can be found")

# ELSI -- only needed when compiled with ELSI support
set(ELSI_ROOT "" CACHE STRING "Root directory of the ELSI installation")

set(ELSI_EXTERNAL_LIBRARIES "" CACHE STRING
  "Any EXTERNAL libraries ELSI needs apart of its own libraries (and scalapack)")
set(ELSI_EXTERNAL_LIBRARY_DIRS "" CACHE STRING
  "Directories where ELSI external libraries can be found")

# PEXSI -- only needed when ELSI was compiled with PEXSI support
# Note: PEXSI usually needs explicit linking of the standard C++ library. Make sure to
#     provide the library path to that C++ standard library, which was used to compile PEXSI.
set(PEXSI_EXTERNAL_LIBRARIES "stdc++" CACHE STRING
  "Any EXTERNAL libraries PEXSI needs apart of its own libraries")
set(PEXSI_EXTERNAL_LIBRARY_DIRS "" CACHE STRING  "Directories with PEXSI external libraries")

# PLUMED -- only needed when compiled with PLUMED support
set(PLUMED_LIBRARIES "plumed;plumedKernel" CACHE STRING "Libraries to link for PLUMED support")
set(PLUMED_LIBRARY_DIRS "" CACHE STRING "Directories to scan for PLUMED libraries")

# MAGMA -- only needed when compiled with GPU support
set(MAGMA_LIBRARIES "magma" CACHE STRING "Magma library")
set(MAGMA_LIBRARY_DIRS "" CACHE STRING "Directories to scan for MAGMA library")
set(MAGMA_INCLUDE_DIRS "" CACHE STRING "Directories to scan for MAGMA include files")

# Any other library needed to be linked or considered as include
set(OTHER_LIBRARIES "" CACHE STRING "Other libraries to link")
set(OTHER_LIBRARY_DIRS "" CACHE STRING "Directories where the other libraries can be found")
set(OTHER_INCLUDE_DIRS "" CACHE STRING "Other include directories to consider")
