set(SOURCES
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

add_library(xmlf90 ${SOURCES})

set(INCLUDEDIR ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(xmlf90 PROPERTIES Fortran_MODULE_DIRECTORY ${INCLUDEDIR})

target_include_directories(xmlf90 INTERFACE ${INCLUDEDIR})
