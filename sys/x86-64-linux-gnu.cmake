set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/volatile/local" CACHE STRING "Install prefix")

set(CMAKE_Fortran_COMPILER "gfortran" CACHE STRING "Fortran compiler")
set(CMAKE_C_COMPILER "gcc" CACHE STRING "C compiler")

set(CMAKE_Fortran_FLAGS "-O2 -mavx" CACHE LIST "Fortran flags")
set(CMAKE_C_FLAGS "-O2" CACHE STRING "C flags")





#set(CMAKE_Fortran_COMPILER gfortran)
#
#set(CMAKE_C_COMPILER gcc)
#
#set(CMAKE_Fortran_FLAGS_RELEASE -O2 -funroll-all-loops -fopenmp)
#set(CMAKE_C_FLAGS_RELEASE -O2 -funroll-all-loops -fall-intrinsics)
#
#set(CMAKE_Fortran_FLAGS_DEBUG -fopenmp -g -Wall -std=f2008 -pedantic -fbounds-check
#  -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-function -Wno-uninitialized)
#set(CMAKE_C_FLAGS_DEBUG -g -Wall -pedantic -fbounds-check)
#
#set(LAPACK_LIBRARIES -llapack -lblas)
#
#set(TEST_OMP_THREADS 1)
#set(TEST_MPI_PROCS 1)
#set(TESTRUNNER env OMP_NUM_THREADS=${TEST_OMP_THREADS})
#
#set(ARPACK_LIBRARIES arpack)
#set(ARPACK_LIBRARY_DIRS "")
