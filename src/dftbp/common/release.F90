!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Release of code constants and types
module dftbp_common_release
  implicit none

  public

  !> Release name of the code
  character(len=*), parameter :: releaseName = '${RELEASE}$'

  !> Year of release
  integer, parameter :: releaseYear = 2024

  !> Mapping between input version and parser version
  type :: TVersionMap
    !> Named version of parser input
    character(10) :: inputVersion
    !> Corresponding numerical version of parser input
    integer :: parserVersion
  end type TVersionMap

end module dftbp_common_release
