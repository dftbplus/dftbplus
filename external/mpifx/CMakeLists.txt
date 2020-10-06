set(BUILD_EXPORTED_TARGETS_ONLY TRUE)
add_subdirectory(origin)
add_library(MpiFx::MpiFx INTERFACE IMPORTED GLOBAL)
target_link_libraries(MpiFx::MpiFx INTERFACE MpiFx)
