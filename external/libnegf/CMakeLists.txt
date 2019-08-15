set(includes "-I${CMAKE_CURRENT_BINARY_DIR}/ext_system")
if(WITH_MPI)
  set(mpi "Y")
  set(includes "${includes} -I${CMAKE_BINARY_DIR}/external/mpifx/origin/src/include")
else()
  set(mpi "N")
endif()

string(REGEX REPLACE ";" " " fypp_flags "${FYPP_FLAGS}")

add_custom_command(
  OUTPUT libnegf-${ARCH}.a libsparskit.a libsyscalls.a
  COMMAND make -C ${CMAKE_CURRENT_BINARY_DIR} -f ${CMAKE_CURRENT_SOURCE_DIR}/origin/Makefile.lib MPI=${mpi} FPP=${FYPP} FPPOPT=${fypp_flags} FXX=${CMAKE_Fortran_COMPILER} FXXOPT=${CMAKE_Fortran_FLAGS} SRCDIR=${CMAKE_CURRENT_SOURCE_DIR}/origin OBJDIR=${CMAKE_CURRENT_BINARY_DIR} INCLUDES=${includes} ARCH=${ARCH}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  VERBATIM)

add_custom_target(build_libnegf DEPENDS libnegf-${ARCH}.a libsparskit.a libsyscalls.a)

add_library(libnegf INTERFACE)

add_dependencies(libnegf build_libnegf)

set(moddirs ${CMAKE_CURRENT_BINARY_DIR}/libnegf ${CMAKE_CURRENT_BINARY_DIR}/ext_system)
set(libfiles
  ${CMAKE_CURRENT_BINARY_DIR}/libnegf-${ARCH}.a
  ${CMAKE_CURRENT_BINARY_DIR}/libsparskit.a
  ${CMAKE_CURRENT_BINARY_DIR}/libsyscalls.a)

target_include_directories(libnegf INTERFACE ${moddirs})
target_link_libraries(libnegf INTERFACE ${libfiles})

if(WITH_MPI)
  add_dependencies(build_libnegf mpifx)
  target_link_libraries(libnegf INTERFACE mpifx)
endif()

foreach(libfile IN LISTS libfiles)
  install(FILES ${libfile} DESTINATION ${INSTALL_LIB_DIR})
endforeach()

foreach(moddir IN LISTS moddirs)
  install(DIRECTORY "${moddir}/" DESTINATION "${INSTALL_MOD_DIR}" FILES_MATCHING PATTERN "*.mod")
endforeach()
