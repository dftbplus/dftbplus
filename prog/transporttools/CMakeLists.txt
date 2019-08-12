add_executable("buildwire" buildwire.f90)
dftbp_install_target(buildwire)

add_executable("flux" flux.f90)
dftbp_install_target(flux)

add_executable("makecube" makecube.f90)
dftbp_install_target(makecube)
