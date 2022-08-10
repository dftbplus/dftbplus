!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Provides a dummy external hamiltonian if no external library is linked
module dftbp_externalham
  use dftbp_common_status, only : TStatus
  implicit none

  private
  public :: TExternalHamiltonian, hamProvides

  type TExternalHamiltonian

  end type TExternalHamiltonian

contains

  subroutine hamprovides(status)

    !> Status of operation
    type(TStatus), intent(out) :: status

    @:RAISE_ERROR(status, -1, "Dummy hamiltonian present, non-functioning calculation")

  end subroutine hamprovides

end module dftbp_externalham
