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
  use dftbp_constants, only : pi
  use dftbp_boundarycond, only : zAxis
  use dftbp_repcont
  use dftbp_quaternions, only : rotate3
  implicit none
  private

  public :: getERep, getERepDeriv

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
  subroutine getERepDeriv(reslt, coords, nNeighbourRep, iNeighbours, species, repCont,&
      & img2CentCell, tHelical)

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

    !> Optional signalling of helical operations
    logical, intent(in), optional :: tHelical

    integer :: iAt1, iNeigh, iAt2, iAt2f
    real(dp) :: vect(3), intermed(3), theta
    logical :: tHelix

    tHelix = .false.
    if (present(tHelical)) then
      tHelix = tHelical
    end if

    reslt(:,:) = 0.0_dp

    if (tHelix) then
      do iAt1 = 1, size(nNeighbourRep)
        lpNeighH: do iNeigh = 1, nNeighbourRep(iAt1)
          iAt2 = iNeighbours(iNeigh,iAt1)
          iAt2f = img2CentCell(iAt2)
          vect(:) = coords(:,iAt1) - coords(:,iAt2)
          call getEnergyDeriv(repCont, intermed, vect, species(iAt1), species(iAt2))
          reslt(:,iAt1) = reslt(:,iAt1) + intermed(:)
          theta = - atan2(coords(2,iAt2),coords(1,iAt2)) &
              & + atan2(coords(2,iAt2f),coords(1,iAt2f))
          theta = mod(theta,2.0_dp*pi)
          call rotate3(intermed, theta, zAxis)
          if (iAt2f /= iAt1) then
            reslt(:,iAt2f) = reslt(:,iAt2f) - intermed(:)
          end if
        end do lpNeighH
      end do

    else

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

    end if

  end subroutine getERepDeriv

end module dftbp_repulsive
