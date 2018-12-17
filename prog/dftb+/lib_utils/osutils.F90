!> Contains utilities interacting with the operating system
module dftbp_osutils
  use iso_c_binding, only : c_char, c_size_t, c_int
  implicit none
  private

  public :: getHostName

  interface
    function c_gethostname(name, maxNameLength) result(status) bind(C, name='gethostname')
      import :: c_char, c_int, c_size_t
      character(kind=c_char), intent(out) :: name(*)
      integer(c_size_t), value, intent(in) :: maxNameLength
      integer(c_int) :: status
    end function c_gethostname

    function c_strnlen(buffer, maxLength) result(length) bind(C, name='strnlen')
      import :: c_char, c_size_t
      character(kind=c_char), intent(in) :: buffer(*)
      integer(c_size_t), value, intent(in) :: maxLength
      integer(c_size_t) :: length
    end function c_strnlen
  end interface


contains

  !> Returns the host name.
  function getHostName() result(hostName)

    !> Name of the host, or the empty string, if could not be determined.
    character(:), allocatable :: hostName

    integer(c_size_t), parameter :: maxHostNameLength = 256
    character(len=maxHostNameLength, kind=c_char) :: buffer
    integer(c_int) :: status

    status = c_gethostname(buffer, maxHostNameLength)
    if (status == 0) then
      hostName = buffer(1:c_strnlen(buffer, maxHostNameLength))
    else
      hostName = ""
    end if

  end function getHostName

end module dftbp_osutils
