set(includes "-I${CMAKE_CURRENT_BINARY_DIR}/ext_system")
if(WITH_MPI)
  set(mpi "Y")
  set(includes "${includes} -I${CMAKE_BINARY_DIR}/external/mpifx/src/include")
else()
  set(mpi "N")
endif()

string(REGEX REPLACE ";" " " fypp_flags "${FYPP_FLAGS}")

if(WITH_MPI)
  set(dependencies mpifx)
else()
  set(dependencies)
endif()

include(ExternalProject)
ExternalProject_Add(
  LibNegfProject
  DEPENDS ${dependencies}
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/origin
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make -C ${CMAKE_CURRENT_BINARY_DIR} -f ${CMAKE_CURRENT_SOURCE_DIR}/origin/Makefile.lib MPI=${mpi} FPP=${FYPP} FPPOPT=${fypp_flags} FXX=${CMAKE_Fortran_COMPILER} FXXOPT=${CMAKE_Fortran_FLAGS} SRCDIR=${CMAKE_CURRENT_SOURCE_DIR}/origin OBJDIR=${CMAKE_CURRENT_BINARY_DIR} INCLUDES=${includes} ARCH=${ARCH}
  INSTALL_COMMAND ""
  TEST_COMMAND "")

add_library(libnegf INTERFACE)
add_dependencies(libnegf LibNegfProject)

target_include_directories(libnegf INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/libnegf>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/ext_system>
  $<INSTALL_INTERFACE:${INSTALL_MOD_DIR}>)

target_link_libraries(libnegf INTERFACE 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/libnegf-${ARCH}.a>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/libsparskit.a>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/libsyscalls.a>
  $<INSTALL_INTERFACE:${INSTALL_LIB_DIR}/libnegf-${ARCH}.a>
  $<INSTALL_INTERFACE:${INSTALL_LIB_DIR}/libsparskit.a>
  $<INSTALL_INTERFACE:${INSTALL_LIB_DIR}/libsyscalls.a>)

if(WITH_MPI)
  target_link_libraries(libnegf INTERFACE mpifx)
endif()


set(libfiles libnegf-${ARCH}.a libsparskit.a libsyscalls.a)
dftbp_prepend_path(libfiles_build ${CMAKE_CURRENT_BINARY_DIR} ${libfiles})
foreach(libfile IN LISTS libfiles_build)
  install(FILES ${libfile} DESTINATION ${INSTALL_LIB_DIR})
endforeach()

set(moddirs ${CMAKE_CURRENT_BINARY_DIR}/libnegf ${CMAKE_CURRENT_BINARY_DIR}/ext_system)
foreach(moddir IN LISTS moddirs)
  install(DIRECTORY "${moddir}/" DESTINATION "${INSTALL_MOD_DIR}" FILES_MATCHING PATTERN "*.mod")
endforeach()

install(TARGETS libnegf EXPORT dftbplus-targets)
