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
  BUILD_COMMAND make -C ${CMAKE_CURRENT_BINARY_DIR} -f ${CMAKE_CURRENT_SOURCE_DIR}/origin/Makefile.lib MPI=${mpi} FPP=${FYPP} FPPOPT=${fypp_flags} FXX=${CMAKE_Fortran_COMPILER} FXXOPT=${CMAKE_Fortran_FLAGS} SRCDIR=${CMAKE_CURRENT_SOURCE_DIR}/origin OBJDIR=${CMAKE_CURRENT_BINARY_DIR} INCLUDES=${includes} ARCH=dftbplus
  INSTALL_COMMAND ""
  TEST_COMMAND "")

add_library(libnegf INTERFACE)
add_dependencies(libnegf LibNegfProject)

target_include_directories(libnegf INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/libnegf>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/ext_system>
  $<INSTALL_INTERFACE:${INSTALL_MOD_DIR}>)

target_link_libraries(libnegf INTERFACE 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/libnegf-dftbplus.a>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/libsparskit.a>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/libsyscalls.a>
  $<INSTALL_INTERFACE:${INSTALL_LIB_DIR}/libnegf-dftbplus.a>
  $<INSTALL_INTERFACE:${INSTALL_LIB_DIR}/libsparskit.a>
  $<INSTALL_INTERFACE:${INSTALL_LIB_DIR}/libsyscalls.a>)

if(WITH_MPI)
  target_link_libraries(libnegf INTERFACE mpifx)
endif()


set(libfiles libnegf-dftbplus.a libsparskit.a libsyscalls.a)
dftbp_add_prefix(${CMAKE_CURRENT_BINARY_DIR}/ "${libfiles}" libfiles_build)
foreach(libfile IN LISTS libfiles_build)
  install(FILES ${libfile} DESTINATION ${INSTALL_LIB_DIR})
endforeach()

set(moddirs ${CMAKE_CURRENT_BINARY_DIR}/libnegf ${CMAKE_CURRENT_BINARY_DIR}/ext_system)
foreach(moddir IN LISTS moddirs)
  install(DIRECTORY "${moddir}/" DESTINATION "${INSTALL_MOD_DIR}" FILES_MATCHING PATTERN "*.mod")
endforeach()

install(TARGETS libnegf EXPORT dftbplus-targets)

set(exported_compiled_libraries negf-dftbplus sparskit syscalls)
# Make sure, mpifx also appears after libnegf in the export file when compiled with MPI support.
# (it will appear twice, but that does hopefully less harm, than only once at the wrong location)
if(WITH_MPI)
  set(exported_compiled_libraries ${exported_compiled_libraries} mpifx)
endif()

set(EXPORTED_COMPILED_LIBRARIES ${EXPORTED_COMPILED_LIBRARIES} ${exported_compiled_libraries}
  PARENT_SCOPE)

