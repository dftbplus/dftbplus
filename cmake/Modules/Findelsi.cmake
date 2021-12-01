find_package(PkgConfig QUIET)
pkg_check_modules(ELSI QUIET elsi)
if(ELSI_FOUND)
  message(STATUS "Found 'elsi' via pkg-config")

  add_library(elsi::elsi INTERFACE IMPORTED)
  target_link_libraries(
    elsi::elsi
    INTERFACE
    "${ELSI_LINK_LIBRARIES}"
  )
  target_include_directories(
    elsi::elsi
    INTERFACE
    "${ELSI_INCLUDE_DIRS}"
  )

  # Set this for ELSI detection in libmbd
  set(ELSI_LIBRARIES "${ELSI_LINK_LIBRARIES}")

  if("${ELSI_LINK_LIBRARIES}" MATCHES "pexsi")
    add_library(elsi::pexsi INTERFACE IMPORTED)
  endif()

  # DFTB+ checks for the lowercase variable name
  set(elsi_VERSION "${ELSI_VERSION}")
endif()
