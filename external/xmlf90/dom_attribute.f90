module xmlf90_dom_attribute

use xmlf90_dom_types
use xmlf90_dom_node
use xmlf90_strings

implicit none

private
  !-------------------------------------------------------  
  ! METHODS FOR ATTRIBUTE NODES
  !------------------------------------------------------- 

  public :: getName
  public :: getValue
  public :: setValue

  interface getValue
    module procedure attrib_getValue
  end interface

  interface setValue
    module procedure attrib_setValue
  end interface

  interface getName
    module procedure attrib_getName
  end interface

CONTAINS

  subroutine attrib_getName(attribute, name)
    type(fnode), intent(in) :: attribute
    type(string), intent(inout) :: name

    if (attribute % nodeType == ATTRIBUTE_NODE) then
       name = attribute%nodeName
    else
      name = ''
    endif

  end subroutine attrib_getName
  
  !-----------------------------------------------------------

  subroutine attrib_getValue(attribute, value)
    type(fnode), intent(in) :: attribute
    type(string), intent(inout):: value

    if (attribute % nodeType == ATTRIBUTE_NODE) then
       value = attribute%nodeValue
    else
       value = ''
    endif
    
  end subroutine attrib_getValue

  !-----------------------------------------------------------

  subroutine attrib_setValue(attribute, value)

    character(len=*), intent(in) :: value
    type(fnode), pointer  :: attribute

    if (attribute % nodeType == ATTRIBUTE_NODE) then
       call setNodeValue(attribute,value)
    endif

  end subroutine attrib_setValue

  !-----------------------------------------------------------


!!! NB Is this a good idea?
!!! NB pure functions have no side effects

  pure function attr_name_len(attribute)
    type(fnode), intent(in) :: attribute
    integer :: attr_name_len
    if (attribute % nodeType == ATTRIBUTE_NODE) then
       attr_name_len = len_trim(attribute % nodeName)
    else
       attr_name_len = 0
    end if
  end function attr_name_len
  
  pure function attr_val_len(attribute)   
    type(fnode), intent(in) :: attribute
    integer :: attr_val_len
    if (attribute % nodeType == ATTRIBUTE_NODE) then
       attr_val_len = len_trim(attribute % nodeValue)
    else
       attr_val_len = 0
    end if
  end function attr_val_len


end module xmlf90_dom_attribute
