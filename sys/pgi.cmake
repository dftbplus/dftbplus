#
# Toolchain file example for
#
# PGI compiler
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
  set(CMAKE_Fortran_COMPILER "pgfortran" CACHE STRING "Fortran compiler")
endif()

# Make sure '-Mallocatable=03' is among the options
set(CMAKE_Fortran_FLAGS "-Mallocatable=03" CACHE STRING "General Fortran flags")

set(CMAKE_Fortran_FLAGS_RELEASE "-O2" CACHE STRING
  "Specific Fortran flags for Release (production) mode")

set(FYPP_FLAGS "-DEMULATE_F08_MATH" CACHE STRING "Flags for the preprocessor")


#
# C compiler settings
#
set(CMAKE_C_COMPILER "pgcc" CACHE STRING "C compiler")

set(CMAKE_C_FLAGS "" CACHE STRING "General C flags")

set(CMAKE_C_FLAGS_RELEASE "-O2" CACHE STRING "Specific C flags for Release mode")

#
# External libraries
#
set(PGI_LIBDIR "/opt/pgi/linux86-64/18.10/lib" CACHE STRING
  "Directory containing PGI libraries")

# LAPACK and BLAS
set(LAPACK_LIBRARIES "lapack blas" CACHE STRING "LAPACK and BLAS libraries to link")
set(LAPACK_LIBRARY_DIRS "${PGI_LIBDIR}" CACHE STRING
  "Directories where LAPACK and BLAS libraries can be found")

# ARPACK -- only needed when built with ARPACK support
set(ARPACK_LIBRARIES "arpack" CACHE STRING "Arpack library")
set(ARPACK_LIBRARY_DIRS "" CACHE STRING "Directories where Arpack library can be found")

# ScaLAPACK -- only needed for MPI-parallel build
set(SCALAPACK_LIBRARIES "scalapack" CACHE STRING "Scalapack libraries to link")
set(SCALAPACK_LIBRARY_DIRS "${PGI_LIBDIR}/scalapack/scalapack-2.0.2/openmpi-2.1.2/lib"
  CACHE STRING "Directories where Scalapack libraries can be found")

#
# NOTE: PGI can not compile the ELSI library, therefore, no ELSI settings are
# given here. If you manage to set up an ELSI library with the PGI-compiler and link it agains
# DFTB+, let us know, so that we can include appropriate settings here.
#

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


#
# Debug settings
#
set(CMAKE_Fortran_FLAGS_DEBUG "-g -C -Mchkptr -traceback" CACHE STRING
  "Specific Fortran flags for Debug mode")

set(CMAKE_C_FLAGS_DEBUG "-g" CACHE STRING "Specific C flags for Debug mode")
