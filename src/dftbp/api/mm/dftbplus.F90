!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> DFTB+ library
module dftbplus
  use dftbp_capi  ! does not export anything but needed for bind(C) routines
  use dftbp_hsdapi, only : dumpHsd, fnode, setChild, setChildValue
  use dftbp_mmapi, only : convertAtomTypesToSpecies, getDftbPlusApi, getDftbPlusBuild,&
      & getMaxAngFromSlakoFile, TDftbPlus, TDftbPlus_destruct, TDftbPlus_init, TDftbPlusAtomList,&
      & TDftbPlusInput, TDftbPlusInput_destruct, TQDepExtPotGen
  implicit none

end module dftbplus
