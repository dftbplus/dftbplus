!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Exporting the functionality we use from the library XMLF90.
module dftbp_extlibs_xmlf90
  use xmlf90_flib_dom, only : append, append_to_string, appendChild, assignment(=), char,&
      & createDocumentNode, createElement, createTextNode, destroyNode, destroyNodeList,&
      & ELEMENT_NODE, fnode, fnodeList, getAttribute, getAttributeNode, getFirstChild, getItem1,&
      & getLastChild, getLength, getNextSibling, getNodeName, getNodeType, getNodeValue,&
      & getParentNode, getPreviousSibling, item, len, normalize, operator(==), parsefile,&
      & prepend_to_string, removeAttribute, removeChild, replaceChild, resize_string, setAttribute,&
      & setTagName, string, TEXT_NODE, textNodeName, trim
  use xmlf90_flib_wxml, only : xml_AddPCData, xml_ADDXMLDeclaration, xml_Close, xml_EndElement,&
      & xml_NewElement, xml_OpenFile, xmlf_t
  implicit none
  public

end module dftbp_extlibs_xmlf90
