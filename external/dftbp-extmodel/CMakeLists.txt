# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

add_subdirectory(${EXTERNALMODEL_SOURCE_DIR} ${EXTERNALMODEL_BINARY_DIR})

add_library(ExternalModel::ExternalModel INTERFACE IMPORTED GLOBAL)
target_link_libraries(ExternalModel::ExternalModel INTERFACE ExternalModel)
