add_executable(buildwire buildwire.f90)
install(TARGETS buildwire DESTINATION ${INSTALL_BIN_DIR})

add_executable(flux flux.f90)
install(TARGETS flux DESTINATION ${INSTALL_BIN_DIR})

add_executable(makecube makecube.f90)
install(TARGETS makecube DESTINATION ${INSTALL_BIN_DIR})
