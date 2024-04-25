add_library(plugin_Fe2 SHARED plugin.c)
target_include_directories(plugin_Fe2 PRIVATE ${CMAKE_SOURCE_DIR}/src/dftbp/plugins)
