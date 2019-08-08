set(sources
  flib_dom.f90
  flib_sax.f90
  flib_wxml.f90
  m_buffer.f90
  m_charset.f90
  m_converters.f90
  m_debug.f90
  m_dictionary.f90
  m_dom_attribute.f90
  m_dom_debug.f90
  m_dom_document.f90
  m_dom_element.f90
  m_dom_error.f90
  m_dom_namednodemap.f90
  m_dom_node.f90
  m_dom_nodelist.f90
  m_dom_parse.f90
  m_dom_types.f90
  m_dom_utils.f90
  m_elstack.f90
  m_entities.f90
  m_fsm.f90
  m_io.f90
  m_reader.f90
  m_strings.f90
  m_wxml_buffer.f90
  m_wxml_core.f90
  m_wxml_dictionary.f90
  m_wxml_elstack.f90
  m_wxml_text.f90
  m_xml_error.f90
  m_xml_parser.f90)

add_library(xmlf90_objlib OBJECT ${sources})
add_library(xmlf90 $<TARGET_OBJECTS:xmlf90_objlib>)

set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(xmlf90_objlib PROPERTIES Fortran_MODULE_DIRECTORY ${includedir})

target_include_directories(xmlf90_objlib INTERFACE ${includedir})
target_include_directories(xmlf90 INTERFACE
  $<TARGET_PROPERTY:xmlf90_objlib,INTERFACE_INCLUDE_DIRECTORIES>)
