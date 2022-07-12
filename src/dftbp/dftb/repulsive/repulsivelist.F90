!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'allocatablelist.fypp'

!> Implements interface for the repulsive (force-field like) potential
module dftbp_dftb_repulsive_repulsivelist
  use dftbp_dftb_repulsive_repulsive, only : TRepulsive
  implicit none

  @:declare_allocatable_list(NAME=TRepulsiveList, TYPE=class(TRepulsive))

contains

  @:implement_allocatable_list(NAME=TRepulsiveList, TYPE=class(TRepulsive))

end module dftbp_dftb_repulsive_repulsivelist
