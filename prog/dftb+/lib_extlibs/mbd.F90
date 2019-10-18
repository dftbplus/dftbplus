#:include "common.fypp"

!> Imports the functionality of libMBD.
module dftbp_mbd
#:if WITH_MBD
  use libmbd
#:endif
  implicit none
  private

#:if WITH_MBD
  public:: TMbdInit, TMbd
  public:: TMbd_init, TMbd_destruct
#:else
  public :: TMbd, TMbdInit
#:endif

#:if not WITH_MBD
  ! Dummy empty type
  type :: TMbdInit
  end type TMbdInit

  ! Dummy empty type
  type :: TMbd
  end type TMbd
#:endif

end module dftbp_mbd
