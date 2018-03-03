!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Exporting the functionality we use from the library dftd3.
module dftd3_module
  use dftd3_api
  implicit none
  private

  public :: dftd3_input, dftd3_calc
  public :: dftd3_init, dftd3_set_params, dftd3_set_functional
  public :: dftd3_dispersion, dftd3_pbc_dispersion
  public :: get_atomic_number
  public :: withDftD3

#:if WITH_DFTD3


  !> Whether code was built with DFTD3 support
  logical, parameter :: withDftD3 = .true.

#:else


  !> Whether code was built with DFTD3 support
  logical, parameter :: withDftD3 = .false.

#:endif

end module dftd3_module
