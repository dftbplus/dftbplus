!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Abstract interface for classical, pure geometry dependent correction schemes
module dftbp_clscorr
  use dftbp_accuracy, only : dp
  use dftbp_periodic, only : TNeighbourList
  implicit none
  private

  public :: TClassicalCorrection


  !> Interface to a classical, pure geometry dependent correction scheme
  type, abstract :: TClassicalCorrection
  contains

    !> Returns the distance over which the classical correction decays
    procedure(getRCutOff), deferred :: getRCutOff

    !> Get energy contributions from classical correction
    procedure(getEnergies), deferred :: getEnergies

    !> Gradient contribution from the classical term
    procedure(addGradients), deferred :: addGradients

    !> The stress tensor contribution from the classical term
    procedure(getStress), deferred :: getStress

  end type TClassicalCorrection


  abstract interface
    !> Returns the distance over which the classical correction decays
    function getRCutOff(this)
      import :: TClassicalCorrection, dp

      !> Instance of the correction
      class(TClassicalCorrection), intent(inout) :: this

      !> Returned distance
      real(dp) :: getRCutOff

    end function getRCutOff


    !> Get energy contributions from classical correction
    subroutine getEnergies(this, atomE, coords, species, neigh, img2CentCell)
      import :: TClassicalCorrection, TNeighbourList, dp

      !> Instance of the correction
      class(TClassicalCorrection), intent(in) :: this

      !> Resulting  energy contributions
      real(dp), intent(out) :: atomE(:)

      !> Current coordinates
      real(dp), intent(in) :: coords(:,:)

      !> Species of the atoms
      integer, intent(in) :: species(:)

      !> Neighbour list
      type(TNeighbourList), intent(in) :: neigh

      !> Updated mapping to central cell.
      integer, intent(in) :: img2CentCell(:)

    end subroutine getEnergies


    !> Gradient contribution from the classical term
    subroutine addGradients(this, derivs, coords, species, neigh, img2CentCell)
      import :: TClassicalCorrection, TNeighbourList, dp

      !> Instance of the correction
      class(TClassicalCorrection), intent(in) :: this

      !> Derivatives to add contribution to to
      real(dp), intent(inout) :: derivs(:,:)

      !> Current coordinates
      real(dp), intent(in) :: coords(:,:)

      !> Chemical species of atoms
      integer, intent(in) :: species(:)

      !> Neighbour list
      type(TNeighbourList), intent(in) :: neigh

      !> Updated mapping to central cell.
      integer, intent(in) :: img2CentCell(:)

    end subroutine addGradients


    !> The stress tensor contribution from the classical term
    subroutine getStress(this, st, coords, neigh, species, img2CentCell, cellVol)
      import :: TClassicalCorrection, TNeighbourList, dp

      !> Instance of the correction
      class(TClassicalCorrection), intent(in) :: this

      !> Stress tensor
      real(dp), intent(out) :: st(:, :)

      !> Coordinates (x,y,z, all atoms including possible images)
      real(dp), intent(in) :: coords(:, :)

      !> Neighbour list
      type(TNeighbourList), intent(in) :: neigh

      !> Species of atoms in the central cell.
      integer, intent(in) :: species(:)

      !> Indexing array for periodic image atoms
      integer, intent(in) :: img2CentCell(:)

      !> Cell volume
      real(dp), intent(in) :: cellVol

    end subroutine getStress
  end interface


end module dftbp_clscorr
