set(BUILD_TESTING FALSE)
add_subdirectory(origin)
add_library(Negf::Negf INTERFACE IMPORTED GLOBAL)
target_link_libraries(Negf::Negf INTERFACE Negf)
