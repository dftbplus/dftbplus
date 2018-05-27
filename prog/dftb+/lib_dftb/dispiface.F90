!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Common interface for all dispersion modules.
module dispiface
  use accuracy, only : dp
  use periodic, only : TNeighbourList
  implicit none
  private
  
  public :: DispersionIface

  !> Interface for classes providing dispersion.
  type, abstract :: DispersionIface
  contains

    !> update internal copy of coordinates
    procedure(updateCoordsIface), deferred :: updateCoords

    !> update internal copy of lattice vectors
    procedure(updateLatVecsIface), deferred :: updateLatVecs

    !> get real space cutoff
    procedure(getRCutoffIface), deferred :: getRCutoff

    !> get energy contributions
    procedure(getEnergiesIface), deferred :: getEnergies

    !> get force contributions
    procedure(addGradientsIface), deferred :: addGradients

    !> get stress tensor contributions
    procedure(getStressIface), deferred :: getStress
  end type DispersionIface


  abstract interface

    !> Update internal stored coordinate
    subroutine updateCoordsIface(this, neigh, img2CentCell, coords, species0)
      import :: DispersionIface, TNeighbourList, dp

      !> data structure
      class(DispersionIface), intent(inout) :: this

      !> list of neighbours to atoms
      type(TNeighbourList), intent(in) :: neigh

      !> image to central cell atom index
      integer, intent(in) :: img2CentCell(:)

      !> atomic coordinates
      real(dp), intent(in) :: coords(:,:)

      !> central cell chemical species
      integer, intent(in) :: species0(:)
    end subroutine updateCoordsIface


    !> update internal copy of lattice vectors
    subroutine updateLatVecsIface(this, latVecs)
      import :: DispersionIface, dp

      !> data structure
      class(DispersionIface), intent(inout) :: this

      !> lattice vectors
      real(dp), intent(in) :: latVecs(:,:)
    end subroutine updateLatVecsIface


    !> get energy contributions
    subroutine getEnergiesIface(this, energies)
      import :: DispersionIface, dp

      !> data structure
      class(DispersionIface), intent(inout) :: this

      !> energy contributions for each atom
      real(dp), intent(out) :: energies(:)
    end subroutine getEnergiesIface


    !> get force contributions
    subroutine addGradientsIface(this, gradients)
      import :: DispersionIface, dp

      !> data structure
      class(DispersionIface), intent(inout) :: this

      !> gradient contributions for each atom
      real(dp), intent(inout) :: gradients(:,:)
    end subroutine addGradientsIface


    !> get stress tensor contributions
    subroutine getStressIface(this, stress)
      import :: DispersionIface, dp

      !> data structure
      class(DispersionIface), intent(inout) :: this

      !> Stress tensor contributions
      real(dp), intent(out) :: stress(:,:)
    end subroutine getStressIface


    !> Distance cut off for dispersion interactions
    function getRCutoffIface(this) result(cutoff)
      import :: DispersionIface, dp

      !> data structure
      class(DispersionIface), intent(inout) :: this

      !> resulting cutoff
      real(dp) :: cutoff
    end function getRCutoffIface

  end interface

end module dispiface
