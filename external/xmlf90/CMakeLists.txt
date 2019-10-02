set(sources
  buffer.f90
  charset.f90
  converters.f90
  debug.f90
  dictionary.f90
  dom_attribute.f90
  dom_debug.f90
  dom_document.f90
  dom_element.f90
  dom_error.f90
  dom_namednodemap.f90
  dom_node.f90
  dom_nodelist.f90
  dom_parse.f90
  dom_types.f90
  dom_utils.f90
  elstack.f90
  entities.f90
  flib_dom.f90
  flib_sax.f90
  flib_wxml.f90
  fsm.f90
  io.f90
  reader.f90
  strings.f90
  wxml_buffer.f90
  wxml_core.f90
  wxml_dictionary.f90
  wxml_elstack.f90
  wxml_text.f90
  xml_error.f90
  xml_parser.f90)

add_library(xmlf90_objlib OBJECT ${sources})

set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(xmlf90_objlib PROPERTIES Fortran_MODULE_DIRECTORY ${includedir})
set_target_properties(xmlf90_objlib PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})

target_include_directories(xmlf90_objlib PUBLIC ${includedir})
