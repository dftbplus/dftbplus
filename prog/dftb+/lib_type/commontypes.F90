!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains some widely used types (at the moment only TOrbitals)
module commontypes
  use orbitals
  use parallelks
  implicit none
  private

  public :: TOrbitals
  public :: TParallelKS, TParallelKS_init

end module commontypes
