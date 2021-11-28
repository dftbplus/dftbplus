!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'pointerlist.fypp'

!> Implements interface for the repulsive (force-field like) potential
module dftbp_dftb_repulsive_repulsivelist
  use dftbp_dftb_repulsive_repulsive, only : TRepulsive
  implicit none

  @:declare_pointer_list(NAME=TRepulsiveList, TYPE=class(TRepulsive))

contains

  @:implement_pointer_list(NAME=TRepulsiveList, TYPE=class(TRepulsive))

end module dftbp_dftb_repulsive_repulsivelist
