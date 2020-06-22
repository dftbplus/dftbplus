#
# Toolchain file example for
#
# NAG compiler
#
# Note the CMake format: Command line options (e.g. compiler flags) space separated, other kind
# of lists semicolon separated.
#

#
# Fortran compiler settings
#
if(WITH_MPI)
  set(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "Fortran compiler")
else()
  set(CMAKE_Fortran_COMPILER "nagfor" CACHE STRING "Fortran compiler")
endif()

set(CMAKE_Fortran_FLAGS "-ieee=full" CACHE STRING "General Fortran flags")

set(CMAKE_Fortran_FLAGS_RELEASE "-O2" CACHE STRING
  "Specific Fortran flags for Release (production) mode")

set(FYPP_FLAGS "-DINTERNAL_ERFC -DEXP_TRAP" CACHE STRING "Flags for the preprocessor")


#
# C compiler settings
#
set(CMAKE_C_COMPILER "gcc" CACHE STRING "C compiler")

set(CMAKE_C_FLAGS "" CACHE STRING "General C flags")

set(CMAKE_C_FLAGS_RELEASE "-O2 -funroll-all-loops -fall-intrinsics" CACHE STRING
  "Specific C flags for Release mode")


#
# External libraries
#

# LAPACK and BLAS
set(LAPACK_LIBRARIES "lapack;blas" CACHE STRING "LAPACK and BLAS libraries to link")
set(LAPACK_LIBRARY_DIRS "" CACHE STRING
  "Directories where LAPACK and BLAS libraries can be found")

# ARPACK -- only needed when built with ARPACK support
set(ARPACK_LIBRARIES "arpack" CACHE STRING "Arpack library")
set(ARPACK_LIBRARY_DIRS "" CACHE STRING "Directories where Arpack library can be found")

# ScaLAPACK -- only needed for MPI-parallel build
set(SCALAPACK_LIBRARIES "scalapack" CACHE STRING "Scalapack libraries to link")
set(SCALAPACK_LIBRARY_DIRS "" CACHE STRING
  "Directories where Scalapack libraries can be found")

# ELSI -- only needed when compiled with ELSI support
set(ELSI_ROOT "/opt/elsi" CACHE STRING "Root directory of the ELSI installation")

set(ELSI_EXTERNAL_LIBRARIES "" CACHE STRING
  "Any EXTERNAL libraries ELSI needs apart of its own libraries (and scalapack)")
set(ELSI_EXTERNAL_LIBRARY_DIRS "" CACHE STRING
  "Directories where ELSI external libraries can be found")

# PEXSI -- only needed when ELSI was compiled with PEXSI support
# Note: PEXSI usually needs explicit linking of the standard C++ library. Make sure to
#     provide the library path to that C++ standard library, which was used to compile PEXSI.
set(PEXSI_EXTERNAL_LIBRARIES "stdc++" CACHE STRING
  "Any EXTERNAL libraries PEXSI needs apart of its own libraries")
set(PEXSI_EXTERNAL_LIBRARY_DIRS "/usr/lib/gcc/x86_64-linux-gnu/7" CACHE STRING
  "Directories with PEXSI external libraries")

# PLUMED -- only needed when compiled with PLUMED support
set(PLUMED_LIBRARIES "plumed;plumedKernel" CACHE STRING "Libraries to link for PLUMED support")
set(PLUMED_LIBRARY_DIRS "" CACHE STRING "Directories to scan for PLUMED libraries")

# MAGMA -- only needed when compiled with GPU support
set(MAGMA_LIBRARIES "magma" CACHE STRING "Magma library")
set(MAGMA_LIBRARY_DIRS "" CACHE STRING "Directories to scan for MAGMA library")
set(MAGMA_INCLUDE_DIRS "/opt/magma/include" CACHE STRING
  "Directories to scan for MAGMA include files")

# Any other library needed to be linked or considered as include
set(OTHER_LIBRARIES "" CACHE STRING "Other libraries to link")
set(OTHER_LIBRARY_DIRS "" CACHE STRING "Directories where the other libraries can be found")
set(OTHER_INCLUDE_DIRS "" CACHE STRING "Other include directories to consider")


#
# Debug settings (for developers)
#
set(CMAKE_Fortran_FLAGS_DEBUG "-f2008 -g -nan -C=all" CACHE STRING
  "Specific Fortran flags for Debug mode")

set(CMAKE_C_FLAGS_DEBUG "-g -Wall -pedantic -fall-intrinsics -fbounds-check" CACHE STRING
  "Specific C flags for Debug mode")

