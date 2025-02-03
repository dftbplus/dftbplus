#
# Toolchain file for
#
# Intel compiler, MKL library
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


#
# Fortran compiler settings
#
if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "IntelLLVM")
  set(Fortran_FLAGS_RELEASE "-O2"
    CACHE STRING "Fortran compiler flags for Release build")
else()
  set(Fortran_FLAGS_RELEASE "-O2 -ip"
    CACHE STRING "Fortran compiler flags for Release build")
endif()

set(Fortran_FLAGS_RELWITHDEBINFO "-g ${Fortran_FLAGS_RELEASE}"
  CACHE STRING "Fortran compiler flags for Release build")

# Note: uninit only works reliably if all linked libraries were compiled using this flag
set(Fortran_FLAGS_DEBUG "-g -O0 -warn all -stand f08 -check all,nouninit -diag-error-limit 1 -traceback"
  CACHE STRING "Fortran compiler flags for Debug build")

# Use intrinsic Fortran 2008 erf/erfc functions
set(INTERNAL_ERFC CACHE BOOL 0)

set(FYPP_FLAGS "" CACHE STRING "Flags for the preprocessor")


#
# C compiler settings
#
set(C_FLAGS "${CMAKE_C_FLAGS}"
  CACHE STRING "Build type independent C compiler flags")

if("${CMAKE_C_COMPILER_ID}" MATCHES "IntelLLVM")
  set(C_FLAGS_RELEASE "-O2"
    CACHE STRING  "C compiler flags for Release build")
else()
  set(C_FLAGS_RELEASE "-O2 -ip"
    CACHE STRING  "C compiler flags for Release build")
endif()

set(C_FLAGS_DEBUG "-g -Wall"
  CACHE STRING "C compiler flags for Debug build")


#
# External libraries
#

# NOTE: Libraries with CMake export files (e.g. ELSI and if the HYBRID_CONFIG_METHODS variable
# contains the "Find" method also libNEGF, libMBD, ScalapackFx and MpiFx) are included by searching
# for the export file in the paths defined in the CMAKE_PREFIX_PATH **environment** variable. Make
# sure your CMAKE_PREFIX_PATH variable is set up accordingly.

# LAPACK and BLAS
# (if the BLAS library contains the LAPACK functions, set LAPACK_LIBRARY to "NONE")

# if(WITH_OMP)
#   set(BLAS_LIBRARY "mkl_intel_lp64;mkl_intel_thread;mkl_core" CACHE STRING "BLAS library to link")
# else()
#   set(BLAS_LIBRARY "mkl_intel_lp64;mkl_sequential;mkl_core" CACHE STRING "BLAS libraries to link")
# endif()
# set(BLAS_LIBRARY_DIR "$ENV{MKLROOT}/lib/intel64" CACHE STRING
#     "Directories where BLAS libraries can be found")

# Automatic CMake LAPACK/BLAS finder settings. If they don't work out, you may try to use the
# manual settings above instead.
if ("${BLAS_LIBRARY}" STREQUAL "")
  # In order to load the library via dlopen() in Python, special MKL library must be linked.
  if(BUILD_SHARED_LIBS AND ENABLE_DYNAMIC_LOADING)
    if(WITH_MPI)
      message(FATAL_ERROR
        "Don't know how to link MKL as shared library with MPI and dlopen (Python API) support. "
        "Set BLAS_LIBRARY, LAPACK_LIBRARY and SCALAPACK_LIBRARY manually, if you know how.")
    endif()
    if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.17")
      set(BLA_VENDOR Intel10_64_dyn)
    else()
      set(BLAS_LIBRARY "mkl_rt")
      set(BLAS_LIBRARY_DIR "$ENV{MKLROOT}/lib/intel64" CACHE STRING
        "Directories where BLAS libraries can be found")
    endif()
  elseif(WITH_OMP)
    set(BLA_VENDOR Intel10_64lp)
  else()
    set(BLA_VENDOR Intel10_64lp_seq)
  endif()
endif()

#set(LAPACK_LIBRARY_DIR "$ENV{MKLROOT}/lib/intel64" CACHE STRING
#    "Directories where LAPACK libraries can be found")
set(LAPACK_LIBRARY "NONE")

# ARPACK -- only needed when built with ARPACK support
#set(ARPACK_LIBRARY "arpack" CACHE STRING "ARPACK library (with path if necessary)")

# PARPACK -- only needed when built with ARPACK and MPI support
#set(PARPACK_LIBRARY "parpack" CACHE STRING "PARPACK library (with path if necessary)")

# ScaLAPACK -- only needed for MPI-parallel build
set(SCALAPACK_LIBRARY "mkl_scalapack_lp64;mkl_blacs_intelmpi_lp64" CACHE STRING
  "Scalapack libraries to link")
set(SCALAPACK_LIBRARY_DIR "$ENV{MKLROOT}/lib/intel64" CACHE STRING
  "Directories where Scalapack libraries can be found")

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
