!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains geometrical boundary condition information on the calculation
module dftbp_boundarycond
  use dftbp_angmomentum, only : rotateZ
  use dftbp_quaternions, only : rotate3
  use dftbp_constants, only : pi
  use dftbp_accuracy, only : dp
  use dftbp_commontypes, only : TOrbitals
  implicit none

  private
  public :: zAxis

  !> z direction vector for rotation
  real(dp), parameter :: zAxis(3) = [0.0_dp,0.0_dp,1.0_dp]

end module dftbp_boundarycond
