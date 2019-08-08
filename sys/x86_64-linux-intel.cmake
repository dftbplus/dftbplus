#
# Fortran compiler settings
#
if(WITH_MPI)
  set(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "Fortran compiler")
else()
  set(CMAKE_Fortran_COMPILER "ifort" CACHE STRING "Fortran compiler")
endif()

set(CMAKE_Fortran_FLAGS "-standard-semantics -heap-arrays 10" CACHE LIST "General Fortran flags")

set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -ip" CACHE LIST
  "Specific Fortran flags for Release (production) mode")


#
# C compiler settings
#
set(CMAKE_C_COMPILER "icc" CACHE STRING "C compiler")

set(CMAKE_C_FLAGS "" CACHE LIST "General C flags")

set(CMAKE_C_FLAGS_RELEASE "-O2 -ip" CACHE LIST "Specific C flags for Release mode")


#
# External libraries
#
if(WITH_MPI)

  set(SCALAPACK_LIBRARIES "mkl_scalapack_lp64;mkl_blacs_intelmpi_lp64" CACHE STRING
    "Scalapack libraries to link")
  set(SCALAPACK_LIBRARY_PATHS "$ENV{MKLROOT}/lib/intel64" CACHE STRING
    "Directories where Scalapack libraries can be found")

  set(LAPACK_LIBRARIES "mkl_intel_lp64;mkl_intel_thread;mkl_core" CACHE STRING
    "LAPACK and BLAS libraries to link")
  set(LAPACK_LIBRARY_PATHS "$ENV{MKLROOT}/lib/intel64" CACHE
    STRING "Directories where LAPACK and BLAS libraries can be found")

  set(OTHER_LIBRARIES "" CACHE STRING "Other libraries to link")
  set(OTHER_LIBRARY_PATHS "" CACHE STRING "Directories where the other libraries can be found")

else(WITH_MPI)

  set(LAPACK_LIBRARIES "mkl_intel_lp64;mkl_intel_thread;mkl_core" CACHE STRING
    "LAPACK and BLAS libraries to link")
  set(LAPACK_LIBRARY_PATHS "$ENV{MKLROOT}/lib/intel64" CACHE
    STRING "Directories where LAPACK and BLAS libraries can be found")

  set(OTHER_LIBRARIES "" CACHE STRING "Other libraries to link")
  set(OTHER_LIBRARY_PATHS "" CACHE STRING "Directories where the other libraries can be found")

endif(WITH_MPI)


if(WITH_ARPACK)
  set(ARPACK_LIBRARIES "arpack" CACHE STRING "Arpack library")
  set(ARPACK_LIBRARY_PATHS "" CACHE STRING "Directories where Arpack library can be found")
endif(WITH_ARPACK)

set(TEST_OMP_THREADS "1" CACHE STRING "Number of threads to use for each test")

set(TEST_MPI_PROCS "1" CACHE STRING "Number of mpi processes to use for each test")


#
# Debug settings
#
set(CMAKE_Fortran_FLAGS_DEBUG "-g -warn all -stand f08 -check -diag-error-limit 1 -traceback"
  CACHE LIST "Specific Fortran flags for Debug mode")

set(CMAKE_C_FLAGS_DEBUG "-g -warn all -check")
