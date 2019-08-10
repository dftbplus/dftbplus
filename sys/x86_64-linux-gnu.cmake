#
# Fortran compiler settings
#
if(WITH_MPI)
  set(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "Fortran compiler")
else()
  set(CMAKE_Fortran_COMPILER "gfortran" CACHE STRING "Fortran compiler")
endif()

set(CMAKE_Fortran_FLAGS "" CACHE STRING "General Fortran flags")

set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -mavx -funroll-all-loops" CACHE STRING
  "Specific Fortran flags for Release (production) mode")

set(FYPP_FLAGS "-n --line-marker-format=gfortran5" CACHE STRING "Flags for the preprocessor")


#
# C compiler settings
#
set(CMAKE_C_COMPILER "gcc" CACHE STRING "C compiler")

set(CMAKE_C_FLAGS "" CACHE STRING "General C flags")

set(CMAKE_C_FLAGS_RELEASE "-O2 -funroll-all-loops" CACHE STRING "Specific C flags for Release mode")


#
# External libraries
#

# LAPACK and BLAS
set(LAPACK_LIBRARIES "lapack blas" CACHE STRING "LAPACK and BLAS libraries to link")
set(LAPACK_LIBRARY_DIRS "" CACHE STRING
  "Directories where LAPACK and BLAS libraries can be found")

# ARPACK -- only needed when built with ARPACK support
set(ARPACK_LIBRARIES "arpack" CACHE STRING "Arpack libraries")
set(ARPACK_LIBRARY_DIRS "" CACHE STRING "Directories where Arpack library can be found")

# ScaLAPACK -- only needed for MPI-parallel build
set(SCALAPACK_LIBRARIES "scalapack" CACHE STRING "Scalapack libraries to link")
set(SCALAPACK_LIBRARY_DIRS "" CACHE STRING
  "Directories where Scalapack libraries can be found")

# ELSI -- only needed when compiled with ELSI support
set(ELSI_ROOT "/opt/elsi" CACHE STRING "Root directory of the ELSI installation")
set(ELSI_LIBRARIES "elsi OMM MatrixSwitch elpa NTPoly fortjson" CACHE STRING
  "ELSI libraries to link")
set(ELSI_LIBRARY_DIRS "${ELSI_ROOT}/lib" CACHE STRING "Directories with ELSI libraries")
set(ELSI_INCLUDE_DIRS "${ELSI_ROOT}/include" CACHE STRING "ELSI include directories")

# PEXSI -- only needed when compiled with PEXSI support
# Note: PEXSI usually needs explicit linking of the standard C++ library. Make sure to
#     provide the library path to that C++ standard library, which was used when PEXSI was
#     compiled.
set(PEXSI_LIBRARIES
  "pexsi superlu_dist ptscotchparmetis ptscotch ptscotcherr scotchmetis scotch scotcherr stdc++"
  CACHE STRING "PEXSI libraries to link")
set(PEXSI_LIBRARY_DIRS "${ELSI_ROOT}/lib /usr/lib/gcc/x86_64-linux-gnu/7" CACHE STRING
  "Directories with PEXSI libraries")
set(PEXSI_INCLUDE_DIRS "${ELSI_ROOT}/include" CACHE STRING "PEXSI include directories")

# Any other library needed to be linked or considered as include
set(OTHER_LIBRARIES "" CACHE STRING "Other libraries to link")
set(OTHER_LIBRARY_DIRS "" CACHE STRING "Directories where the other libraries can be found")
set(OTHER_INCLUDE_DIRS "" CACHE STRING "Other include directories to consider")


#
# Test settings
#

# Number of OMP-threads used by one test process
set(TEST_OMP_THREADS "1" CACHE STRING "Number of threads to use for each test")

# Number of MPI-processed launched for each test
set(TEST_MPI_PROCS "1" CACHE STRING "Number of mpi processes to use for each test")

# Command line used to launch the test code
if(WITH_MPI)
  set(TEST_RUNNER "env OMP_NUM_THREADS=${TEST_OMP_THREADS} mpiexec -n ${TEST_MPI_PROCS}")
else()
  set(TEST_RUNNER "env OMP_NUM_THREADS=${TEST_OMP_THREADS}")
endif()


#
# Debug settings (for developers)
#
set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -std=f2008 -pedantic -fbounds-check" CACHE STRING
  "Specific Fortran flags for Debug mode")

set(CMAKE_C_FLAGS_DEBUG "-g -Wall -pedantic -fall-intrinsics -fbounds-check" CACHE STRING
  "Specific C flags for Debug mode")
