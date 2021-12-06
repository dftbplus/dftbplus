# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!
# Note: ChIMES currently can only be invoked as a submodule!

set(WITH_FORTRAN08_API True)

add_subdirectory(${CHIMESCALC_SOURCE_DIR} ${CHIMESCALC_BINARY_DIR})

add_library(ChimesCalc::ChimesCalc_Fortran INTERFACE IMPORTED GLOBAL)
target_link_libraries(ChimesCalc::ChimesCalc_Fortran INTERFACE ChimesCalc_Fortran)
