cmake_minimum_required(VERSION 3.5)

set(DFTBPLUS_VERSION "19.2")

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(dftbpUtils)

include(${CMAKE_SOURCE_DIR}/config.cmake)

dftbp_ensure_config_consistency()

project(dftbplus VERSION ${DFTBPLUS_VERSION} LANGUAGES Fortran C)

set(FYPP "${PROJECT_SOURCE_DIR}/external/fypp/bin/fypp" CACHE FILEPATH "Fypp preprocessor")
dftbp_convert_to_list(FYPP_FLAGS)
dftbp_add_fypp_defines(FYPP_FLAGS)

set(PYTHON_INTERPRETER "python" CACHE STRING
  "Python interpreter to use for installing and test python components")

dftbp_get_release_name(RELEASE)


if(WITH_MPI)
  find_package(MPI REQUIRED)
endif()

if(WITH_OMP)
  find_package(OpenMP REQUIRED)
  string(APPEND CMAKE_C_FLAGS " ${OpenMP_C_FLAGS}")
  string(APPEND CMAKE_Fortran_FLAGS " ${OpenMP_Fortran_FLAGS}")
  string(APPEND CMAKE_EXE_LINKER_FLAGS " ${OpenMP_EXE_LINKER_FLAGS}")
endif()


#
# Check external libraries
#
dftbp_convert_to_list(SCALAPACK_LIBRARIES SCALAPACK_LIBRARY_DIRS LAPACK_LIBRARIES
  LAPACK_LIBRARY_DIRS ARPACK_LIBRARIES ARPACK_LIBRARY_DIRS ELSI_EXTERNAL_LIBRARIES
  ELSI_EXTERNAL_LIBRARY_DIRS PEXSI_EXTERNAL_LIBRARIES PEXSI_EXTERNAL_LIBRARY_DIRS
  OTHER_LIBRARIES OTHER_LIBRARY_DIRS)

set(EXPORTED_EXTERNAL_LIBRARIES)
set(EXPORTED_EXTERNAL_LIBRARY_DIRS)
set(EXPORTED_EXTERNAL_PACKAGES)
set(EXPORTED_EXPORTED_COMPILED_LIBRARIES)

if(WITH_MPI)
  dftbp_create_library_targets("${SCALAPACK_LIBRARIES}" "${SCALAPACK_LIBRARY_DIRS}")
  list(APPEND EXPORTED_EXTERNAL_LIBRARIES ${SCALAPACK_LIBRARIES})
  list(APPEND EXPORTED_EXTERNAL_LIBRARY_DIRS ${SCALAPACK_LIBRARY_DIRS})
endif()

if(WITH_ELSI)
  find_package(elsi REQUIRED PATHS ${ELSI_ROOT})
  dftbp_create_library_targets("${ELSI_EXTERNAL_LIBRARIES}" "${ELSI_EXTERNAL_LIBRARY_DIRS}")
  list(APPEND EXPORTED_EXTERNAL_PACKAGES elsi)
  if(TARGET elsi::pexsi)
    set(ELSI_WITH_PEXSI TRUE)
  else()
    set(ELSI_WITH_PEXSI FALSE)
  endif()
  if(ELSI_WITH_PEXSI)
    dftbp_create_library_targets("${PEXSI_EXTERNAL_LIBRARIES}" "${PEXSI_EXTERNAL_LIBRARY_DIRS}")
  endif()
endif()

dftbp_create_library_targets("${LAPACK_LIBRARIES}" "${LAPACK_LIBRARY_DIRS}")
list(APPEND EXPORTED_EXTERNAL_LIBRARIES ${LAPACK_LIBRARIES})
list(APPEND EXPORTED_EXTERNAL_LIBRARY_DIRS ${LAPACK_LIBRARY_DIRS})

if(WITH_ARPACK)
  dftbp_create_library_targets("${ARPACK_LIBRARIES}" "${ARPACK_LIBRARY_DIRS}")
  list(APPEND EXPORTED_EXTERNAL_LIBRARIES ${ARPACK_LIBRARIES})
  list(APPEND EXPORTED_EXTERNAL_LIBRARY_DIRS ${ARPACK_LIBRARY_DIRS})
endif()

dftbp_create_library_targets("${OTHER_LIBRARIES}" "${OTHER_LIBRARY_DIRS}")
list(APPEND EXPORTED_EXTERNAL_LIBRARIES ${OTHER_LIBRARIES})
list(APPEND EXPORTED_EXTERNAL_LIBRARY_DIRS ${OTHER_LIBRARY_DIRS})


#
# Add optional external components
#
if(BUILD_API)
  set(exclude)
else()
  set(exclude EXCLUDE_FROM_ALL)
endif()

add_subdirectory(external/xmlf90 EXCLUDE_FROM_ALL)

if(WITH_SOCKETS)
  add_subdirectory(external/fsockets EXCLUDE_FROM_ALL)
endif()

if(WITH_DFTD3)
  add_subdirectory(external/dftd3 ${exclude})
endif()

if(WITH_MPI)
  add_subdirectory(external/mpifx ${exclude})
  add_subdirectory(external/scalapackfx ${exclude})
endif()

if(WITH_TRANSPORT)
  add_subdirectory(external/mudpack ${exclude})
  add_subdirectory(external/libnegf ${exclude})
endif()


#
# Add internal components
#
add_subdirectory(prog)
add_subdirectory(tools/dptools)

enable_testing()
add_subdirectory(test)


#
# Export package info
# 
if(BUILD_API)
  install(EXPORT dftbplus-targets FILE DftbPlusTargets.cmake NAMESPACE DftbPlus::
    DESTINATION ${INSTALL_LIB_DIR}/cmake/DftbPlus)
  
  configure_file(utils/export/DftbPlusConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/DftbPlusConfig.cmake @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/DftbPlusConfig.cmake DESTINATION
    ${INSTALL_LIB_DIR}/cmake/DftbPlus)

  set(PKGCONFIG_LIBS -ldftbplus)
  foreach(complib IN LISTS EXPORTED_COMPILED_LIBRARIES)
    set(PKGCONFIG_LIBS "${PKGCONFIG_LIBS} -l${complib}")
  endforeach()
  set(PKGCONFIG_LIBS_PRIVATE)
  foreach(extlibdir IN LISTS EXPORTED_EXTERNAL_LIBRARY_DIRS)
    set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -L${extlibdir}")
  endforeach()
  foreach(extlib IN LISTS EXPORTED_EXTERNAL_LIBRARIES)
    set(PKGCONFIG_LIBS_PRIVATE "${PKGCONFIG_LIBS_PRIVATE} -l${extlib}")
  endforeach()
  string(REPLACE ";" " " PKGCONFIG_REQUIRES "${EXPORTED_EXTERNAL_PACKAGES}")
  configure_file(utils/export/dftbplus.pc.in ${CMAKE_CURRENT_BINARY_DIR}/dftbplus.pc @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/dftbplus.pc DESTINATION ${INSTALL_LIB_DIR}/pkgconfig)
  
  #set(API_VERSION "0.2")
  #include(CMakePackageConfigHelpers)
  #write_basic_package_version_file(
  #  ${INSTALL_CMAKE_DIR}/DftbPlusConfigVersion.cmake
  #  VERSION ${API_VERSION}
  #  COMPATIBILITY SameMajorVersion)
endif()
