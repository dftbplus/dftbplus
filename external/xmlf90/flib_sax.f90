module xmlf90_flib_sax

!
! Stub module to gather all the functionality needed by the user
! 
! In future m_dictionary and m_converters could be exported by
! other parts of a more general fortran library.
!
! m_xml_error is necessary in order to use a custom error handler.
!
use xmlf90_dictionary
use xmlf90_xml_parser
use xmlf90_converters
use xmlf90_xml_error

implicit none

public

end module xmlf90_flib_sax

