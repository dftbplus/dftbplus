!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for linear response derivative calculations using perturbation methods with respect to
!> atom locations
module dftbp_perturbxderivs
  use dftbp_accuracy
  use dftbp_constants
  use dftbp_globalenv
  use dftbp_message
  use dftbp_commontypes
  use dftbp_potentials
  use dftbp_scc
  use dftbp_orbitalequiv
  use dftbp_populations
  use dftbp_spin
  use dftbp_thirdorder_module, only : ThirdOrder
  use dftbp_dftbplusu
  use dftbp_onsitecorrection
  use dftbp_mainio
  use dftbp_shift
  use dftbp_mixer
  use dftbp_finitethelper
  use dftbp_scalapackfx
  use dftbp_environment
  use dftbp_periodic
  use dftbp_densedescr
  use dftbp_sparse2dense
  use dftbp_nonscc, only : NonSccDiff
  use dftbp_rotateDegenerateOrbs, only : TDegeneracyTransform
  use dftbp_taggedoutput
  use dftbp_wrappedintr
  use dftbp_taggedoutput, only : TTaggedWriter
  use dftbp_slakocont
#:if WITH_MPI
  use dftbp_mpifx
#:endif
#:if WITH_SCALAPACK
  use dftbp_scalafxext
#:else
  use dftbp_blasroutines
#:endif

  implicit none

  private
  public :: dPsidx

contains

  !> Calculate derivatives with respect to atomic positions
  subroutine dPsidx(env, parallelKS, eigvals, eigVecsReal, eigVecsCplx, ham, over, skHamCont,&
      & skOverCont, nonSccDeriv, orb, nAtom, species, neighbourList, nNeighbourSK, denseDesc,&
      & iSparseStart, img2CentCell, coords, kPoint, kWeight, cellVec, iCellVec, latVec,&
      & taggedWriter, tWriteAutoTest, autoTestTagFile, tWriteTaggedOut, taggedResultsFile,&
      & tWriteDetailedOut, fdDetailedOut, eigValsDeltaq, eigVecsCplxDeltaq)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> ground state eigenvectors in real case
    real(dp), intent(in), allocatable :: eigVecsReal(:,:,:)

    !> ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Container for SK Hamiltonian integrals
    type(OSlakoCont), intent(in) :: skHamCont

    !> Container for SK overlap integrals
    type(OSlakoCont), intent(in) :: skOverCont

    !> method for calculating derivatives of S and H0
    type(NonSccDiff), intent(in) :: nonSccDeriv

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> chemical species
    integer, intent(in) :: species(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Lattice vectors for the unit cell
    real(dp), intent(in) :: latVec(:,:)

    !> Tagged writer object
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> should regression test data be written
    logical, intent(in) :: tWriteAutoTest

    !> File name for regression data
    character(*), intent(in) :: autoTestTagFile

    !> should machine readable output data be written
    logical, intent(in) :: tWriteTaggedOut

    !> File name for machine readable results data
    character(*), intent(in) :: taggedResultsFile

    !> should detailed.out be written to
    logical, intent(in) :: tWriteDetailedOut

    !> File id for detailed.out
    integer, intent(in) :: fdDetailedOut

    !> Eigenvalue at k+q for each level, kpoint and spin channel
    real(dp), intent(in), allocatable :: eigValsDeltaq(:,:,:)

    !> Ground state complex eigenvectors at k+q
    complex(dp), intent(in), allocatable :: eigvecsCplxDeltaq(:,:,:)

    integer :: ii, iAt
    real(dp) :: sPrime(size(over),3)

  #:if WITH_SCALAPACK

    call error("MPI parallel not added yet")

  #:else

    if (allocated(eigVecsReal)) then

      do iAt = 1, nAtom

        call nonSccDeriv%getFirstDeriv(Sprime, skOverCont, coords, species, iAt, orb, nNeighbourSK,&
            & neighbourList%iNeighbour, iSparseStart, img2centcell)

        write(stdOut, *) 'Atom ', iAt
        do ii = 1, 3
          write(stdOut, *)sPrime(:,ii)
          write(stdOut, *)
        end do

      end do

      !subroutine getFirstDeriv(this, deriv, skCont, coords, species, iAt1, orb, nNeighbourSK,&
      !    & iNeighbours, iPair, img2centcell)



    else
      
      if (allocated(eigvecsCplxDeltaq)) then

        call error("q/=0 Not added yet")

      else

        call error("q=0 not added yet")

      end if

    end if

  #:endif

  end subroutine dPsidx

end module dftbp_perturbxderivs
