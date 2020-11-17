#
# Toolchain file for
#
# GNU compiler with MKL-library using static linking (distributed binary package)
#
# FOR EXPERTS ONLY! This is a highly customized toolchain file to enforce static build with the
# MKL-library within an ideal build environment. It may fail or even misbehave on your system for
# various reasons.
#

# Switching on all components compatible with static linking and binary distribution.
set(WITH_OMP TRUE CACHE BOOL "Whether OpenMP thread parallisation should be enabled" FORCE)
set(WITH_TRANSPORT TRUE CACHE BOOL "Whether transport via libNEGF should be included." FORCE)
set(WITH_API TRUE CACHE BOOL "Whether API should be built" FORCE)
set(WITH_ARPACK TRUE CACHE BOOL
  "Whether the ARPACK library should be included (needed for TD-DFTB)" FORCE)
set(WITH_PLUMED TRUE CACHE BOOL
  "Whether metadynamics via the PLUMED2 library should be allowed for" FORCE)
set(WITH_MBD TRUE CACHE BOOL "Whether MBD library should be included" FORCE)
set(WITH_SOCKETS TRUE CACHE BOOL "Whether socket communication should be allowed for" FORCE)

# Switching on all components not compatible with static linking or binary distribution.
set(WITH_MPI FALSE CACHE BOOL "Whether DFTB+ should support MPI-parallelism" FORCE)
set(WITH_ELSI FALSE CACHE BOOL "Whether DFTB+ with MPI-parallelism should use the ELSI libraries"
  FORCE)
set(WITH_GPU FALSE CACHE BOOL
  "Whether DFTB+ should support GPU-acceleration via the MAGMA-library" FORCE)
# Not compatible due to license issues (GPL vs LGPL)
set(WITH_DFTD3 FALSE CACHE BOOL "Whether the DFTD3 library should be included" FORCE)

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

# Make sure everything is statically linked
set(BUILD_SHARED_LIBS FALSE CACHE BOOL "Whether the libraries built should be shared" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "-static" CACHE STRING "Linker flags" FORCE)

# Prevent OpenMP detection by providing the compiler flags directly
# (Otherwise CMake links the shared OpenMP libraries)
add_library(OpenMP::OpenMP_Fortran INTERFACE IMPORTED)
target_compile_options(OpenMP::OpenMP_Fortran INTERFACE "-fopenmp")
target_link_options(OpenMP::OpenMP_Fortran INTERFACE "-fopenmp")


#
# External libraries
#

# NOTE: Libraries with CMake export files (e.g. ELSI and if the HYBRID_CONFIG_METHODS variable
# contains the "Find" method also libNEGF, libMBD, ScalapackFx and MpiFx) are included by searching
# for the export file in the paths defined in the CMAKE_PREFIX_PATH **environment** variable. Make
# sure your CMAKE_PREFIX_PATH variable is set up accordingly.

# LAPACK and BLAS
set(LAPACK_LIBRARY "-Wl,--start-group -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group"
  CACHE STRING "LAPACK libraries")
#set(LAPACK_LIBRARY_DIR "" CACHE STRING
#  "Directories where LAPACK and BLAS libraries can be found")

# ARPACK -- only needed when built with ARPACK support
set(ARPACK_LIBRARY "libarpack.a" CACHE STRING "Arpack libraries")
#set(ARPACK_LIBRARY_DIR "" CACHE STRING "Directories where Arpack library can be found")

# ScaLAPACK -- only needed for MPI-parallel build
#set(SCALAPACK_LIBRARY "scalapack-openmpi" CACHE STRING "Scalapack libraries to link")
#set(SCALAPACK_LIBRARY_DIR "" CACHE STRING "Directories where Scalapack libraries can be found")

# NOTE: The libraries below provide Pkg-Conf export files.  If your PKG_CONFIG_PATH environment
# variable has been set up correctly (containing the paths to these libraries), no adjustment should
# be necessary below.

# PLUMED -- only needed when compiled with PLUMED support
set(PLUMED_LIBRARY "libplumed.a;libplumedWrapper.a;libstdc++.a;libz.a"
      CACHE STRING "Plumed libraries")
#set(PLUMED_LIBRARY_DIR "" CACHE STRING "Directories to scan for PLUMED libraries")

# MAGMA -- only needed when compiled with GPU support
#set(MAGMA_LIBRARY "magma" CACHE STRING "Magma library")
#set(MAGMA_LIBRARY_DIR "" CACHE STRING "Directories to scan for MAGMA library")
#set(MAGMA_INCLUDE_DIRECTORY "" CACHE STRING "Directories to scan for MAGMA include files")
