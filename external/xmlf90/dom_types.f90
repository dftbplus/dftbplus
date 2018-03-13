Module xmlf90_dom_types

  use xmlf90_strings

  implicit none

  private

  !-------------------------------------------------------   
  ! A GENERIC NODE
  !-------------------------------------------------------   
  type, public :: fnode
     type(string)         :: nodeName
     type(string)         :: nodeValue
!!!     character(len=200)    :: nodeName  = ""
!!!     character(len=200)    :: nodeValue = ""
     integer              :: nc              = 0 
     integer              :: nodeType        = 0
     type(fnode), pointer :: parentNode      => null()
     type(fnodeList), pointer :: childNodes  => null()  ! New
     type(fnode), pointer :: firstChild      => null()
     type(fnode), pointer :: lastChild       => null()
     type(fnode), pointer :: previousSibling => null()
     type(fnode), pointer :: nextSibling     => null()
     type(fnode), pointer :: ownerDocument   => null()
     type(fnamedNodeMap), pointer :: attributes => null()
  end type fnode

  !-----------------------------------------------------------
  !  ONE WAY TO IMPLEMENT A NAMEDNODEMAP  (dictionary)
  !-----------------------------------------------------------

  ! Linked list of name/node pairs, with overall length variable

  type, public :: fnamedNode
     type(string)                   :: name
!!!     character(len=100)            :: name
     type(fnode), pointer          :: node => null()
     type(fnamedNode), pointer     :: next => null()
  end type fnamedNode

  type, public :: fnamedNodeMap
     integer :: length = 0
     type(fnamedNode), pointer  :: head => null()
     type(fnamedNode), pointer  :: tail => null()
  end type fnamedNodeMap

  !-----------------------------------------------------------
  !  ONE WAY TO IMPLEMENT A NODELIST 
  !-----------------------------------------------------------

  type, public :: flistNode
     type(fnode), pointer          :: node => null()
     type(flistNode), pointer      :: next => null()
  end type flistNode

  type, public :: fnodeList
     integer                      :: length = 0
     type(flistNode), pointer     :: head => null()
     type(flistNode), pointer     :: tail => null()
     type(flistNode), pointer     :: pCache => null()
     integer                      :: iCache = 0
  end type fnodeList

!========================================================================
  integer, save, private          :: allocated_nodes = 0
!========================================================================

  !-------------------------------------------------------   
  ! NODETYPES
  !-------------------------------------------------------   
  integer, parameter, public :: ELEMENT_NODE                = 1
  integer, parameter, public :: ATTRIBUTE_NODE              = 2
  integer, parameter, public :: TEXT_NODE                   = 3
  integer, parameter, public :: CDATA_SECTION_NODE          = 4
  integer, parameter, public :: ENTITY_REFERENCE_NODE       = 5
  integer, parameter, public :: ENTITY_NODE                 = 6
  integer, parameter, public :: PROCESSING_INSTRUCTION_NODE = 7
  integer, parameter, public :: COMMENT_NODE                = 8
  integer, parameter, public :: DOCUMENT_NODE               = 9
  integer, parameter, public :: DOCUMENT_TYPE_NODE          = 10
  integer, parameter, public :: DOCUMENT_FRAGMENT_NODE      = 11
  integer, parameter, public :: NOTATION_NODE               = 12

  public :: node_class
  public :: createNode
  public :: destroyNode
  public :: destroyNamedNodeMap
  public :: destroyNodeList
  public :: getNumberofAllocatedNodes

CONTAINS

  function getNumberofAllocatedNodes() result(n)
    integer   :: n

    n = allocated_nodes
  end function getNumberofAllocatedNodes

!--------------------------------------------------------------
  function createNode() result(node)
    type(fnode), pointer  :: node

    allocate(node)
    allocated_nodes = allocated_nodes + 1

  end function createNode
!--------------------------------------------------------------

  function node_class(nodetype) result(class)
    integer, intent(in) :: nodetype
    character(len=10)  ::     class

    select case(nodetype)
    case(ELEMENT_NODE)
       class = "element"
    case(ATTRIBUTE_NODE)
       class = "attribute"
    case(TEXT_NODE)
       class = "text"
    case(COMMENT_NODE)
       class = "comment"
    case(DOCUMENT_NODE)
       class = "document"
    end select
  end function node_class

  subroutine destroyNamedNodeMap(nodemap, destroyNodes)
    type(fnamedNodeMap), pointer :: nodemap
    logical, intent(in), optional :: destroyNodes

    type(fnamednode), pointer  :: nnp, nnp2
    type(fnode), pointer       :: ghost
    logical :: tDestroyNodes

    if (present(destroyNodes)) then
      tDestroyNodes = destroyNodes
    else
      tDestroyNodes = .false.
    end if
    
    if (.not. associated(nodemap)) return
    nnp => nodemap%head
    do while (associated(nnp))
       ghost => nnp%node
       nnp2 => nnp%next
       if (tDestroyNodes) then
         call destroyNode(ghost)
       end if
       deallocate(nnp)
       nnp => nnp2
     enddo
     deallocate(nodemap)
     
  end subroutine destroyNamedNodeMap

  subroutine destroyNodeList(nodelist)
    type(fnodeList), pointer :: nodelist

    type(flistnode), pointer   :: p, ghost
    
    if (.not. associated(nodelist)) return
    p => nodelist%head
    do while (associated(p))
       ghost => p
       p => p%next
       deallocate(ghost)
     enddo
     deallocate(nodelist)
    
  end subroutine destroyNodeList

  recursive subroutine destroyNode(node)
    type(fnode), pointer  :: node
    
    type(fnode), pointer  :: np, ghost
    
    np => node
    do while (associated(np))
      if (associated(np%firstChild)) then
          call destroyNode(np%firstChild)
        endif
        if (associated(np%attributes)) then
          call destroyNamedNodeMap(np%attributes, destroyNodes=.true.)
        end if
       if (associated(np%previousSibling)) then
         np%previousSibling%nextSibling => np%nextSibling
       end if
       if (associated(np%nextSibling)) then
         np%nextSibling%previousSibling => np%previousSibling
       end if
       if (associated(np%parentNode)) then
         if (associated(np%parentNode%firstChild,np)) then
            np%parentNode%firstChild => null()
          end if
          if (associated(np%parentNode%lastChild,np)) then
            np%parentNode%lastChild => null()
          end if
        endif
       if (associated(np,node)) then    
          deallocate(np)
          allocated_nodes = allocated_nodes - 1
          EXIT                           ! do not destroy siblings
        else
          ghost => np
          np => np%nextSibling
          deallocate(ghost)
          allocated_nodes = allocated_nodes - 1
       endif
    enddo
    node => null()     ! superfluous ?

  end subroutine destroyNode

end module xmlf90_dom_types

