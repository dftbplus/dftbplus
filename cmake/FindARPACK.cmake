# see https://github.com/dune-project/dune-istl/blob/master/cmake/modules/FindARPACK.cmake
# for the possibility of a proper FindARPACK module.

cmake_push_check_state(RESET)

if(ARPACK_LIBRARIES)
  message("-- Using user settings for ARPACK: '${ARPACK_LIBRARIES}'")
else()
  message(FATAL_ERROR "You must specify ARPACK_LIBRARIES when building with ARPACK")
endif()

cmake_pop_check_state()
