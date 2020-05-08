!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines to calculate repulsive pair contributions to energy and forces
module dftbp_repulsive
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_repcont
  implicit none
  private

  public :: getERep, getERepDeriv, getERepDeriv2nd


contains


  !> Subroutine for repulsive energy contributions for each atom
  subroutine getERep(reslt, coords, nNeighbourRep, iNeighbours, species, repCont, img2CentCell)

    !> Energy for each atom.
    real(dp), intent(out) :: reslt(:)

    !> coordinates (x,y,z, all atoms including possible images)
    real(dp), intent(in) :: coords(:,:)

    !> Number of neighbours for atoms in the central cell
    integer, intent(in) :: nNeighbourRep(:)

    !> Index of neighbours for a given atom.
    integer, intent(in) :: iNeighbours(0:,:)

    !> Species of atoms in the central cell.
    integer, intent(in) :: species(:)

    !> Container for repulsive potentials.
    type(TRepCont), intent(in) :: repCont

    !> Index of each atom in the central cell which the atom is mapped on.
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iNeigh, iAt2, iAt2f
    real(dp) :: vect(3), dist, intermed

  @:ASSERT(size(reslt) == size(nNeighbourRep))

    reslt(:) = 0.0_dp
    do iAt1 = 1, size(nNeighbourRep)
      do iNeigh = 1, nNeighbourRep(iAt1)
        iAt2 = iNeighbours(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        vect(:) = coords(:,iAt1) - coords(:,iAt2)
        dist = sqrt(sum(vect**2))
        call getEnergy(repCont, intermed, dist, species(iAt1), species(iAt2))
        reslt(iAt1) = reslt(iAt1) + 0.5_dp * intermed
        if (iAt2f /= iAt1) then
          reslt(iAt2f) = reslt(iAt2f) + 0.5_dp * intermed
        end if
      end do
    end do

  end subroutine getERep


  !> Subroutine for force contributions of the repulsives.
  subroutine getERepDeriv(reslt, coords, nNeighbourRep, iNeighbours, species, repCont, img2CentCell)

    !> Energy for each atom.
    real(dp), intent(out) :: reslt(:,:)

    !> coordinates (x,y,z, all atoms including possible images)
    real(dp), intent(in) :: coords(:,:)

    !> Number of neighbours for atoms in the central cell
    integer, intent(in) :: nNeighbourRep(:)

    !> Index of neighbours for a given atom.
    integer, intent(in) :: iNeighbours(0:,:)

    !> Species of atoms in the central cell.
    integer, intent(in) :: species(:)

    !> Container for repulsive potentials.
    type(TRepCont), intent(in) :: repCont

    !> Index of each atom in the central cell which the atom is mapped on.
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iNeigh, iAt2, iAt2f
    real(dp) :: vect(3), intermed(3)

  @:ASSERT(size(reslt,dim=1) == 3)

    reslt(:,:) = 0.0_dp
    do iAt1 = 1, size(nNeighbourRep)
      lpNeigh: do iNeigh = 1, nNeighbourRep(iAt1)
        iAt2 = iNeighbours(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        if (iAt2f == iAt1) then
          cycle lpNeigh
        end if
        vect(:) = coords(:,iAt1) - coords(:,iAt2)
        call getEnergyDeriv(repCont, intermed, vect, species(iAt1), species(iAt2))
        reslt(:,iAt1) = reslt(:,iAt1) + intermed(:)
        reslt(:,iAt2f) = reslt(:,iAt2f) - intermed(:)
      end do lpNeigh
    end do

  end subroutine getERepDeriv


  !> Subroutine for second derivative contributions from the repulsives.
  subroutine getERepDeriv2nd(reslt,coords, nNeighbors, iNeighbors, species, repCont, img2CentCell)

    !> Second derivatives for each atom.
    real(dp), intent(out) :: reslt(:,:,:)

    !> Coordinates (x,y,z, all atoms including possible images)
    real(dp), intent(in) :: coords(:,:)

    !> Number of neighbors for atoms in the central cell
    integer, intent(in) :: nNeighbors(:)

    !> Index of neighbors for a given atom.
    integer, intent(in) :: iNeighbors(0:,:)

    !> Species of atoms in the central cell.
    integer, intent(in) :: species(:)

    !> Container for repulsive potentials.
    type(TRepCont), intent(in) :: repCont

    !> Index of each atom in the central cell which the atom is mapped on to.
    integer, intent(in) :: img2CentCell(:)

    integer  :: iAt1, iNeigh, iAt2, iAt2f
    real(dp) :: vect(3), junk(3), intermed(3,3)

  @:ASSERT(size(reslt,dim=1) == 3)
  @:ASSERT(size(reslt,dim=2) == 3)
  @:ASSERT(size(reslt,dim=3) == size(nNeighbors))

    reslt = 0.0_dp
    do iAt1 = 1, size(nNeighbors)
      lpNeigh: do iNeigh = 1, nNeighbors(iAt1)
        iAt2 = iNeighbors(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)

        if (iAt2f == iAt1) then
          ! q=0 case
          cycle lpNeigh
        end if

        vect = coords(:,iAt1) - coords(:,iAt2)

        call getEnergyDeriv(repCont, junk, vect, species(iAt1), species(iAt2),intermed)

        reslt(:,:,iAt1) = reslt(:,:,iAt1) + intermed
        reslt(:,:,iAt2f) = reslt(:,:,iAt2f) + intermed

      end do lpNeigh
    end do

  end subroutine getERepDeriv2nd


end module dftbp_repulsive
