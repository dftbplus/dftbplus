module xmlf90_dom_error

implicit none
  
  !-------------------------------------------------------
  ! EXCEPTION CODES
  !-------------------------------------------------------
  integer, parameter, public :: INDEX_SIZE_ERR              = 1
  integer, parameter, public :: DOMSTRING_SIZE_ERR          = 2
  integer, parameter, public :: HIERARCHY_REQUEST_ERR       = 3
  integer, parameter, public :: WRONG_DOCUMENT_ERR          = 4
  integer, parameter, public :: INVALID_CHARACTER_ERR       = 5
  integer, parameter, public :: NO_DATA_ALLOWED_ERR         = 6
  integer, parameter, public :: NO_MODIFICATION_ALLOWED_ERR = 7
  integer, parameter, public :: NOT_FOUND_ERR               = 8
  integer, parameter, public :: NOT_SUPPORTED_ERR           = 9
  integer, parameter, public :: INUSE_ATTRIBUTE_ERR         = 10
  integer, parameter, public :: INVALID_STATE_ERR           = 11
  integer, parameter, public :: SYNTAX_ERR                  = 12
  integer, parameter, public :: INVALID_MODIFICATION_ERR    = 13
  integer, parameter, public :: NAMESPACE_ERR               = 14
  integer, parameter, public :: INVALID_ACCESS_ERR          = 15
  integer, parameter, public :: VALIDATION_ERR              = 16
  integer, parameter, public :: TYPE_MISMATCH_ERR           = 17

CONTAINS

  subroutine dom_error(name,code,msg)
    character(len=*), intent(in) :: name, msg
    integer, intent(in)          :: code

    print *, "***ERROR***"
    print *, "Routine ", trim(name), ":", trim(msg)
    print *, 1.0 / sin(3.141592654)
  end subroutine dom_error

end module xmlf90_dom_error
