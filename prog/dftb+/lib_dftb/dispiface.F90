!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Common interface for all dispersion modules.
!!
module dispiface
  use accuracy, only : dp
  use periodic, only : TNeighborList
  implicit none

  !> Interface for classes providing dispersion.
  !!
  type, abstract :: DispersionIface
  contains
    procedure(updateCoordsIface), deferred :: updateCoords
    procedure(updateLatVecsIface), deferred :: updateLatVecs
    procedure(getRCutoffIface), deferred :: getRCutoff
    procedure(getEnergiesIface), deferred :: getEnergies
    procedure(addGradientsIface), deferred :: addGradients
    procedure(getStressIface), deferred :: getStress
  end type DispersionIface


  abstract interface
    subroutine updateCoordsIface(this, neigh, img2CentCell, coords, species0)
      import :: DispersionIface, TNeighborList, dp
      class(DispersionIface), intent(inout) :: this
      type(TNeighborList), intent(in) :: neigh
      integer, intent(in) :: img2CentCell(:)
      real(dp), intent(in) :: coords(:,:)
      integer, intent(in) ::  species0(:)
    end subroutine updateCoordsIface

    subroutine updateLatVecsIface(this, latVecs)
      import :: DispersionIface, dp
      class(DispersionIface), intent(inout) :: this
      real(dp), intent(in) :: latVecs(:,:)
    end subroutine updateLatVecsIface

    subroutine getEnergiesIface(this, energies)
      import :: DispersionIface, dp
      class(DispersionIface), intent(inout) :: this
      real(dp), intent(out) :: energies(:)
    end subroutine getEnergiesIface

    subroutine addGradientsIface(this, gradients)
      import :: DispersionIface, dp
      class(DispersionIface), intent(inout) :: this
      real(dp), intent(inout) :: gradients(:,:)
    end subroutine addGradientsIface

    subroutine getStressIface(this, stress)
      import :: DispersionIface, dp
      class(DispersionIface), intent(inout) :: this
      real(dp), intent(out) :: stress(:,:)
    end subroutine getStressIface

    function getRCutoffIface(this) result(cutoff)
      import :: DispersionIface, dp
      class(DispersionIface), intent(inout) :: this
      real(dp) :: cutoff
    end function getRCutoffIface
  end interface

end module dispiface
