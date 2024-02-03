!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a two body repulsive potential interface
module dftbp_dftb_repulsive_pairrepulsive
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TPairRepulsive, TPairRepulsiveItem


  !> Abstract type for the interface of pairwise repulsives
  type, abstract :: TPairRepulsive
  contains
  procedure(TPairRepulsive_getCutoff), deferred :: getCutoff
  procedure(TPairRepulsive_getValue), deferred :: getValue
  end type TPairRepulsive


  !> Wrapper for a TPairRepulsive class item
  type :: TPairRepulsiveItem
    class(TPairRepulsive), allocatable :: item
  end type TPairRepulsiveItem


  abstract interface

    !> Returns the real-space cutoff of the two-body repulsive
    function TPairRepulsive_getCutoff(this) result(cutoff)
      import :: TPairRepulsive, dp
      implicit none

      !> Instance
      class(TPairRepulsive), intent(in) :: this

      !> Real space cutoff
      real(dp) :: cutoff

    end function TPairRepulsive_getCutoff


    !> Returns energy of the two-body repulsive for a given distance
    subroutine TPairRepulsive_getValue(this, rr, energy, dEnergy, d2Energy)
      import :: TPairRepulsive, dp
      implicit none

      !> Instance
      class(TPairRepulsive), intent(in) :: this

      !> Distance between interacting atoms
      real(dp), intent(in) :: rr

      !> Energy contribution
      real(dp), optional, intent(out) :: energy

      !> First derivative Energy contribution
      real(dp), optional, intent(out) :: dEnergy

      !> Second derivative energy contribution
      real(dp), optional, intent(out) :: d2Energy

    end subroutine TPairRepulsive_getValue

  end interface


end module dftbp_dftb_repulsive_pairrepulsive
