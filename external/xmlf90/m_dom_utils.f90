module m_dom_utils

  use m_dom_types
  use m_dom_element
  use m_dom_document
  use m_dom_node
  use m_dom_namednodemap
  use m_dom_debug
  use m_strings

  use flib_wxml

implicit none
  
  public :: dumpTree
  public :: xmlize

  private

CONTAINS

  subroutine dumpTree(startNode)

    type(fnode), pointer :: startNode   

    character(len=50) :: indent = " "
    integer           :: indent_level
    type(string)      :: s

    indent_level = 0
    call dump2(startNode)

  contains

    recursive subroutine dump2(input)
      type(fnode), pointer :: input
      type(fnode), pointer :: temp     
      temp => input
      do while(associated(temp))
        call getNodeName(temp, s)
         write(*,'(3a,i3)') indent(1:indent_level), &
                        char(s), " of type ", &
                        getNodeType(temp)
         if (hasChildNodes(temp)) then
            indent_level = indent_level + 3
            call dump2(getFirstChild(temp))
            indent_level = indent_level - 3
         endif
         temp => getNextSibling(temp)
      enddo

    end subroutine dump2

  end subroutine dumpTree
!----------------------------------------------------------------

  subroutine xmlize(startNode,fname)

    type(fnode), pointer :: startNode   
    character(len=*), intent(in) :: fname

    type(xmlf_t)  :: xf
    type(string)  :: s, sv, sn       ! to avoid memory leaks

    call xml_OpenFile(fname,xf)
    call dump_xml(startNode)
    call xml_Close(xf)

  contains

    recursive subroutine dump_xml(input)
      type(fnode), pointer         :: input
!
!     Just this node and its descendants, no siblings. 
!     Of course, the document node has only children...
!
      type(fnode), pointer         :: node, attr
      type(fnamedNodeMap), pointer :: attr_map
      integer  ::  i

      node => input
      do
         if (.not. associated(node)) exit
         select case (getNodeType(node))

          case (DOCUMENT_NODE)

             call xml_AddXMLDeclaration(xf)
             if (hasChildNodes(node)) call dump_xml(getFirstChild(node))

          case (ELEMENT_NODE)

             call getNodeName(node, s)
             call xml_NewElement(xf,char(s))
             attr_map => getAttributes(node)
             do i = 0, getLength(attr_map) - 1
                attr => item(attr_map,i)
                call getNodeName(attr, sn)
                call getNodeValue(attr, sv)
                call xml_AddAttribute(xf, char(sn),char(sv))
             enddo
             if (hasChildNodes(node)) call dump_xml(getFirstChild(node))
             call getNodeName(node, s)
             call xml_EndElement(xf,char(s))

          case (TEXT_NODE)
             
             call getNodeValue(node,s)
             call xml_AddPcdata(xf,char(s))

          case (CDATA_SECTION_NODE)
             
             call getNodeValue(node, s)
             call xml_AddCdataSection(xf,char(s))

          case (COMMENT_NODE)
             
             call getNodeValue(node, s)
             call xml_AddComment(xf,char(s))

        end select
        if (associated(node,StartNode)) exit  ! In case we request the 
                                              ! dumping of a single element,
                                              ! do not do siblings
        node => getNextSibling(node)
     enddo

    end subroutine dump_xml

  end subroutine xmlize

end module m_dom_utils
