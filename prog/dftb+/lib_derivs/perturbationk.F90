!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for linear response derivative calculations using perturbation methods
module dftbp_perturbkderivs
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
  use dftbp_rotateDegenerateOrbs, only : TDegeneracyTransform
  use dftbp_taggedoutput
  use dftbp_wrappedintr
  use dftbp_taggedoutput, only : TTaggedWriter
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
  public :: dPsidK

contains

  subroutine dPsidK(env, parallelKS, eigvals, eigVecsCplx, ham, over, orb, nAtom, species,&
      & neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell, coord, kPoint, kWeight,&
      & cellVec, iCellVec, latVec, taggedWriter, tWriteAutoTest, autoTestTagFile, tWriteTaggedOut,&
      & taggedResultsFile, tWriteDetailedOut, fdDetailedOut)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

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
    real(dp), intent(in) :: coord(:,:)

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

    integer :: nIndepHam, nSpin, nOrbs, nKpts, iAt1, iNeigh, iAt2, iAt2f, ii, jj, kk, iCell, iSpin
    integer :: iDir, iK, iOrb, iKS, iS, iCart
    real(dp), allocatable :: dHam(:,:), dOver(:), dEi(:,:,:,:)
    real(dp) :: vec(3)
    complex(dp), allocatable :: dPsi(:,:,:,:), dHamSqr(:,:), dOverSqr(:,:), workLocal(:,:)
    complex(dp), allocatable :: hamSqr(:,:,:), overSqr(:,:,:), eCi(:,:), work2Local(:,:)
    complex(dp), allocatable :: work3Local(:,:), eigvecsTransformed(:,:)
    logical :: tTransformed
    integer :: fdResults

    !complex(dp), allocatable :: U(:,:)
    !type(wrappedCmplx2), allocatable :: UBlock(:)
    !integer, allocatable :: blockRange(:,:)

    type(TDegeneracyTransform) :: transform

  #:if WITH_SCALAPACK
    ! need distributed matrix descriptors
    integer :: desc(DLEN_), nn

    type(blocklist) :: blocks
    integer :: iLoc, iGlob, jGlob, blockSize


  #:endif

    if (.not.allocated(eigVecsCplx)) then
      call error("Missing complex eigenvectors")
    end if

    nOrbs = size(eigVals,dim=1)
    nSpin = size(eigVals,dim=3)
    nKpts = size(eigVals,dim=2)

    allocate(dEi(nOrbs, nKpts, nSpin, 3))
    dEi(:,:,:,:) = 0.0_dp

  #:if WITH_SCALAPACK


    nn = denseDesc%fullSize
    call scalafx_getdescriptor(env%blacs%orbitalGrid, nn, nn, env%blacs%rowBlockSize,&
        & env%blacs%columnBlockSize, desc)

    allocate(dHamSqr(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2)))
    allocate(dOverSqr(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2)))
    allocate(eCi(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2)))
    allocate(workLocal(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2)))
    allocate(work2local(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2)))

    call transform%init()

    do iKS = 1, parallelKS%nLocalKS


      iK = parallelKS%localKS(1, iKS)
      iS = parallelKS%localKS(2, iKS)


      call blocks%init(env%blacs%orbitalGrid, desc, "c")

      do iCart = 1, 3

        eigvecsTransformed = eigVecsCplx(:,:,iKS)

        call unpackHSdkBlacs(env%blacs, ham(:,iS), kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, dHamSqr,&
            & iCart)

        call unpackHSdkBlacs(env%blacs, over, kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, dOverSqr,&
            & iCart)

        ! e_i |c_i>
        do ii = 1, size(blocks)
          call blocks%getblock(ii, iGlob, iLoc, blockSize)
          do jj = 0, min(blockSize - 1, nOrbs - iGlob)
            eCi(:, iLoc + jj) = eigVals(iGlob + jj, iK, iS) * eigVecsCplx(:, iLoc + jj, iKS)
          end do
        end do

        ! H' <c|
        call pblasfx_phemm(dHamSqr, denseDesc%blacsOrbSqr, eigVecsCplx(:,:,iKS),&
            & denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr)

        ! H' - e S' <c|
        call pblasfx_phemm(dOverSqr, denseDesc%blacsOrbSqr, eCi,&
            & denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr, alpha=(-1.0_dp,0.0_dp),&
            & beta=(1.0_dp,0.0_dp))

        ! |c> H' - e S' <c|
        call pblasfx_pgemm(workLocal, denseDesc%blacsOrbSqr, eigVecsCplx(:,:,iKS),&
            & denseDesc%blacsOrbSqr, work2local, denseDesc%blacsOrbSqr, transa="C")

        call transform%generateUnitary(env, work2local, eigvals(:,iK,iS),&
            & eigVecsTransformed, denseDesc, tTransformed)

        ! now have states orthogonalised agains the operator in degenerate cases, |c~>

        !if (tTransformed) then
          ! re-form |c> H' - e S' <c| with the transformed vectors

          ! e_i |c~_i>
          do ii = 1, size(blocks)
            call blocks%getblock(ii, iGlob, iLoc, blockSize)
            do jj = 0, min(blockSize - 1, nOrbs - iGlob)
              eCi(:, iLoc + jj) = eigVals(iGlob + jj, iK, iS) * eigVecsTransformed(:, iLoc + jj)
            end do
          end do

          ! H' <c~|
          call pblasfx_phemm(dHamSqr, denseDesc%blacsOrbSqr, eigVecsTransformed,&
              & denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr)

          ! H' - e S' <c~|
          call pblasfx_phemm(dOverSqr, denseDesc%blacsOrbSqr, eCi,&
              & denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr, alpha=(-1.0_dp,0.0_dp),&
              & beta=(1.0_dp,0.0_dp))

          ! |c~> H' - e S' <c~|
          call pblasfx_pgemm(workLocal, denseDesc%blacsOrbSqr, eigVecsTransformed,&
              & denseDesc%blacsOrbSqr, work2local, denseDesc%blacsOrbSqr, transa="C")
        !end if

        if (allocated(dEi)) then
          ! derivative of eigenvalues stored in diagonal of matrix work2Local, from <c|h'|c>
          do jj = 1, size(work2Local,dim=2)
            jGlob = scalafx_indxl2g(jj, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
                & env%blacs%orbitalGrid%ncol)
            do ii = 1, size(work2Local,dim=1)
              iGlob = scalafx_indxl2g(ii, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_),&
                  & env%blacs%orbitalGrid%nrow)
              if (iGlob == jGlob) then
                write(*,*)iGlob,real(work2Local(ii,jj),dp) * Hartree__eV
                dEi(iGlob, iK, iS, iCart) = real(work2Local(ii,jj),dp)
              end if
            end do
          end do
        end if
        write(*,*)



      end do

    end do

    if (allocated(dEi)) then
      call mpifx_allreduceip(env%mpi%globalComm, dEi, MPI_SUM)
    end if

  #:else

    allocate(dHamSqr(nOrbs, nOrbs))
    allocate(dOverSqr(nOrbs, nOrbs))
    allocate(eCi(nOrbs, nOrbs))

    allocate(workLocal(nOrbs, nOrbs))
    allocate(work2Local(nOrbs, nOrbs))
    allocate(work3Local(nOrbs, nOrbs))

    call transform%init()

    do iKS = 1, parallelKS%nLocalKS

      iK = parallelKS%localKS(1, iKS)
      iS = parallelKS%localKS(2, iKS)

      ! e_i <c_i|
      do iOrb = 1, nOrbs
        eCi(:,iOrb) = eigvals(iOrb,iK,iS) * eigvecsCplx(:,iOrb,iKS)
      end do

      do iCart = 1, 3

        dOverSqr(:,:) = cmplx(0,0,dp)
        dHamSqr(:,:) = cmplx(0,0,dp)
        call unpackHSdk(dOverSqr, over, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
            & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell, iCart)
        call unpackHSdk(dHamSqr, ham(:,iS), kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
            & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell, iCart)

        ! form |c> H' - e S' <c|
        call hemm(workLocal, 'l', dHamSqr, eigVecsCplx(:,:,iKS))
        call hemm(workLocal, 'l', dOverSqr, eCi, alpha=(-1.0_dp,0.0_dp), beta=(1.0_dp,0.0_dp))
        workLocal(:,:) = matmul(transpose(conjg(eigVecsCplx(:,:,iKS))), workLocal)

        call transform%generateUnitary(workLocal, eigvals(:,iK,iS))
        call transform%degenerateTransform(workLocal)

        ! diagonal elements of workLocal are now derivatives of eigenvalues if needed
        if (allocated(dEi)) then
          do ii = 1, nOrbs
            dEi(ii, iK, iS, iCart) = real(workLocal(ii,ii), dp)
          end do
        end if

        work2Local(:,:) = eigVecsCplx(:,:,iKS)
        call transform%applyUnitary(work2Local)

        do ii = 1, nOrbs
          workLocal(ii,ii) = 0.0_dp
          do jj = ii + 1, nOrbs
            if (.not.transform%degenerate(ii,jj)) then
              workLocal(jj,ii) = workLocal(jj,ii) / (eigvals(ii,iK,iS) - eigvals(jj,iK,iS))
              workLocal(ii, jj) = -workLocal(jj,ii)
            end if
          end do
        end do

        work3Local(:,:) = matmul(work2Local, workLocal)

        write(stdOut,*)'Derivatives of ci'
        do ii = 1, nOrbs
          write(stdOut,"(16E12.2)")work3Local(:,ii)
        end do
        write(stdOut,*)

        work2Local(:,1) = sum(conjg(eigVecsCplx(:,:,iKS)) * work3Local, dim = 1)
        write(stdOut,*)'product'
        do ii = 1, nOrbs
          write(stdOut,*)ii, work2Local(ii,1)
        end do
        write(stdOut,*)

      end do

    end do

    call transform%destroy()

  #:endif

    if (allocated(dEi) .and. env%tGlobalMaster) then
      if (tWriteAutoTest) then
        open(newunit=fdResults, file=trim(autoTestTagFile), position="append")
        call taggedWriter%write(fdResults, tagLabels%dEigenDk, dEi)
        close(fdResults)
      end if
      if (tWriteTaggedOut) then
        open(newunit=fdResults, file=trim(taggedResultsFile), action="write", status="old",&
            & position="append")
        call taggedWriter%write(fdResults, tagLabels%dEigenDk, dEi)
        close(fdResults)
      end if
      if (tWriteDetailedOut) then
        write(fdDetailedOut,*)'Linear response derivatives or eigenvalues (eV) wrt to k(x,y,z)'
        do iK = 1, nKpts
          do iS = 1, nSpin
            write(fdDetailedOut,"(1X,A,I2,A,I0)")'Spin ',iS,' kpoint ', iK
            do iOrb = 1, nOrbs
              write(fdDetailedOut,"(I6,3F10.3)")iOrb, Hartree__eV * dEi(iOrb, iK, iS, :)
            end do
          end do
        end do
        write(fdDetailedOut,*)
      end if
    end if

  end subroutine dPsidK

end module dftbp_perturbkderivs
