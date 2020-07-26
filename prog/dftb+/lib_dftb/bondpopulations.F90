!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Evaluates bond populations from sparse matrices. To do: add evaluation of bond orders, perhaps
!> DOI: 10.1039/c7ra07400j can be adapted for DFTB.
module dftbp_bondpops
  use dftbp_accuracy, only : dp

  private
  public :: addPairWiseBondInfo

contains

  !> Calculates properties per atom pair. If given the overlap as sparseMat, returns electron count
  !> on atoms and in bonds (summing to total number of electrons). H0 returns the non-SCC bond
  !> energy.
  !> Note: In periodic systems, the bond contributions from image atoms are included in the
  !> central-cell bond.
  subroutine addPairWiseBondInfo(info, rhoPrim, sparseMat, iSquare, iNeighbour, nNeighbourSK,&
      & img2CentCell, iSparseStart)

    !> pairwise contribution
    real(dp), allocatable, intent(inout) :: info(:,:)

    !> sparse density matrix (only spin-unpolarized)
    real(dp), intent(in) :: rhoPrim(:)

    !> non-scc hamiltonian or overlap matrix in sparse format
    real(dp), intent(in) :: sparseMat(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> Atomic neighbour data
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(0:,:)

    integer :: iAt1, iAt2, iAt2f, nOrb1, nOrb2, iOrig, iNeigh, iOrb1, iOrb2
    integer :: nAtom

    nAtom = size(iSquare) - 1
    !$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(RUNTIME) REDUCTION(+:info)&
    !$OMP& PRIVATE(iOrb1, nOrb1, iNeigh, iOrig, iAt2, iAt2f, iOrb2, nOrb2, tmp)
    do iAt1 = 1, nAtom
      iOrb1 = iSquare(iAt1)
      nOrb1 = iSquare(iAt1+1) - iOrb1
      do iNeigh = 0, nNeighbourSK(iAt1)
        iOrig = iSparseStart(iNeigh,iAt1) + 1
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iOrb2 = iSquare(iAt2f)
        nOrb2 = iSquare(iAt2f+1) - iOrb2
        info(iAt2f,iAt1) = info(iAt2f,iAt1)&
            & + sum(rhoPrim(iOrig:iOrig+nOrb1*nOrb2-1) * sparseMat(iOrig:iOrig+nOrb1*nOrb2-1))
      end do
    end do
    !$OMP  END PARALLEL DO

    ! fill other triangle
    do iAt1 = 1, nAtom - 1
      info(iAt1, iAt1+1:) = info(iAt1+1:, iAt1)
    end do

  end subroutine addPairWiseBondInfo

end module dftbp_bondpops
