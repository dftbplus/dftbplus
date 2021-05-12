!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains some widely used types (at the moment only TOrbitals)
module dftbp_type_commontypes
  use dftbp_type_orbitals
  use dftbp_type_parallelks
  implicit none
  private

  public :: TOrbitals
  public :: TParallelKS, TParallelKS_init

end module dftbp_type_commontypes
