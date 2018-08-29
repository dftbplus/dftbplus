#:include "common.fypp"

!> Imports the functionality of libMBD.
module mbd_module
#:if WITH_MBD
  use mbd_api
#:endif
  implicit none

#:if not WITH_MBD
  ! Dummy empty type
  type :: mbd_input
  end type mbd_input

  ! Dummy empty type
  type :: mbd_calculation
  end type mbd_calculation
#:endif

end module mbd_module
