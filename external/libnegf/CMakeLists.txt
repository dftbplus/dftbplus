macro(before_source_config)
  set(BUILD_TESTING FALSE)
endmacro()

macro(after_source_config)
  add_library(Negf::Negf INTERFACE IMPORTED GLOBAL)
  target_link_libraries(Negf::Negf INTERFACE Negf)
endmacro()

dftbp_config_hybrid_source_dependency(Negf NEGF)
