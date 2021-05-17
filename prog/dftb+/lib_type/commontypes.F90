!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains some widely used types (at the moment only TOrbitals)
module dftbp_commontypes
  use dftbp_orbitals, only : TOrbitals
  use dftbp_parallelks, only : TParallelKS, TParallelKS_init
  implicit none
  
  private
  public :: TOrbitals
  public :: TParallelKS, TParallelKS_init

end module dftbp_commontypes
