set(installed_scripts
  calc_timeprop_maxpoldir
  calc_timeprop_spectrum)

install(FILES ${installed_scripts} DESTINATION "${CMAKE_INSTALL_BINDIR}"
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
