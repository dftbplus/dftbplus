#:include "common.fypp"

!> Imports the functionality of libMBD.
module mbd_module
#:if WITH_MBD
  use mbd, TMbdInit => mbd_input_t, TMbd => mbd_calc_t
#:endif
  implicit none
  private

  public:: TMbdInit, TMbd

#:if not WITH_MBD
  ! Dummy empty type
  type :: TMbdInit
  end type TMbdInit

  ! Dummy empty type
  type :: TMbd
  end type TMbd
#:endif

end module mbd_module
