!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains some widely used types
module dftbp_type_commontypes
  use dftbp_type_orbitals, only : TOrbitals
  use dftbp_type_parallelks, only : TParallelKS, TParallelKS_init
  implicit none

  private
  public :: TOrbitals
  public :: TParallelKS, TParallelKS_init

end module dftbp_type_commontypes
