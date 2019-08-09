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


#
# C compiler settings
#
set(CMAKE_C_COMPILER "gcc" CACHE STRING "C compiler")

set(CMAKE_C_FLAGS "" CACHE STRING "General C flags")

set(CMAKE_C_FLAGS_RELEASE "-O2 -funroll-all-loops" CACHE STRING "Specific C flags for Release mode")


#
# External libraries
#
if(WITH_MPI)

  set(SCALAPACK_LIBRARIES "scalapack" CACHE STRING "Scalapack libraries to link")
  set(SCALAPACK_LIBRARY_DIRS "" CACHE STRING
    "Directories where Scalapack libraries can be found")

  set(LAPACK_LIBRARIES "lapack blas" CACHE STRING "LAPACK and BLAS libraries to link")
  set(LAPACK_LIBRARY_DIRS "" CACHE STRING
    "Directories where LAPACK and BLAS libraries can be found")

  set(OTHER_LIBRARIES "" CACHE STRING "Other libraries to link")
  set(OTHER_LIBRARY_DIRS "" CACHE STRING "Directories where the other libraries can be found")

else(WITH_MPI)

  set(LAPACK_LIBRARIES "lapack blas" CACHE STRING "LAPACK and BLAS libraries to link")
  set(LAPACK_LIBRARY_DIRS "" CACHE STRING
    "Directories where LAPACK and BLAS libraries can be found")

  set(OTHER_LIBRARIES "" CACHE STRING "Other libraries to link")
  set(OTHER_LIBRARY_DIRS "" CACHE STRING "Directories where the other libraries can be found")

endif(WITH_MPI)


if(WITH_ARPACK)
  set(ARPACK_LIBRARIES "arpack" CACHE STRING "Arpack libraries")
  set(ARPACK_LIBRARY_DIRS "" CACHE STRING "Directories where Arpack library can be found")
endif(WITH_ARPACK)

set(TEST_OMP_THREADS "1" CACHE STRING "Number of threads to use for each test")

set(TEST_MPI_PROCS "1" CACHE STRING "Number of mpi processes to use for each test")


#
# Debug settings
#
set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -std=f2008 -pedantic -fbounds-check" CACHE STRING
  "Specific Fortran flags for Debug mode")

set(CMAKE_C_FLAGS_DEBUG "-g -Wall -pedantic -fall-intrinsics -fbounds-check" CACHE STRING
  "Specific C flags for Debug mode")
