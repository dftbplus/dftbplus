!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module to wrap around the process of converting from a Hamiltonian and overlap in sparse form
!> into eigenvectors
module eigenvects
  use assert
  use accuracy
  use eigensolver
  use message
  use solvertypes
#:if WITH_SCALAPACK
  use scalapackfx
#:if WITH_ELSI
  use ELSI
#:endif
#:endif
  use parallelks
  use solvers
  implicit none
  private

  public :: diagDenseMtx
#:if WITH_SCALAPACK
  public :: diagDenseMtxBlacs
#:endif


  interface diagDenseMtx
    module procedure diagDenseRealMtx
    module procedure diagDenseComplexMtx
  end interface diagDenseMtx

#:if WITH_SCALAPACK

  interface diagDenseMtxBlacs
    module procedure diagDenseRealMtxBlacs
    module procedure diagDenseCplxMtxBlacs
  end interface diagDenseMtxBlacs

#:endif

contains

  !> Diagonalizes a sparse represented Hamiltonian and overlap to give the eigenvectors and values,
  !> as well as often the Cholesky factorized overlap matrix (due to a side effect of lapack)
  subroutine diagDenseRealMtx(electronicSolver, jobz, HSqrReal, SSqrReal, eigen)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> type of eigen-problem, either 'V'/'v' with vectors or 'N'/'n' eigenvalues only
    character, intent(in) :: jobz

    !> Large square matrix for the resulting eigenvectors
    real(dp), intent(inout) :: HSqrReal(:,:)

    !> Large square matrix for the overlap workspace, often overwritten with the Cholesky factorized
    !> form.
    real(dp), intent(inout) :: SSqrReal(:,:)

    !> Eigen values.
    real(dp), intent(out) :: eigen(:)


    @:ASSERT(size(HSqrReal, dim=1) == size(HSqrReal, dim=2))
    @:ASSERT(all(shape(HSqrReal) == shape(SSqrReal)))
    @:ASSERT(size(HSqrReal, dim=1) == size(eigen))
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

    select case(electronicSolver%iSolver)
    case(electronicSolverTypes%QR)
      call hegv(HSqrReal,SSqrReal,eigen,'L',jobz)
    case(electronicSolverTypes%divideandconquer)
      call hegvd(HSqrReal,SSqrReal,eigen,'L',jobz)
    case(electronicSolverTypes%relativelyrobust)
      call gvr(HSqrReal,SSqrReal,eigen,'L',jobz)
    case default
      call error('Unknown eigensolver')
    end select

  end subroutine diagDenseRealMtx


  !> Diagonalizes a sparse represented Hamiltonian and overlap with k-points to give the
  !> eigenvectors and values, as well as often the Cholesky factorized overlap matrix (due to a side
  !> effect of lapack)
  subroutine diagDenseComplexMtx(electronicSolver, jobz, HSqrCplx, SSqrCplx, eigen)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> type of eigen-problem, either 'V'/'v' vectors or 'N'/'n' eigenvalues only
    character, intent(in) :: jobz

    !> Large square matrix for the resulting eigenvectors
    complex(dp), intent(inout) :: HSqrCplx(:,:)

    !> Large square matrix for the overlap workspace, overwritten with the Cholesky factorized form.
    complex(dp), intent(inout) :: SSqrCplx(:,:)

    !> The eigenvalues of the matrices
    real(dp), intent(out) :: eigen(:)


    @:ASSERT(size(HSqrCplx, dim=1) == size(HSqrCplx, dim=2))
    @:ASSERT(all(shape(HSqrCplx) == shape(SSqrCplx)))
    @:ASSERT(size(HSqrCplx, dim=1) == size(eigen))
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

    select case(electronicSolver%iSolver)
    case(electronicSolverTypes%QR)
      call hegv(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(electronicSolverTypes%divideandconquer)
      call hegvd(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(electronicSolverTypes%relativelyrobust)
      call gvr(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case default
      call error('Unknown eigensolver')
    end select

  end subroutine diagDenseComplexMtx


#:if WITH_SCALAPACK

  !> Diagonalizes a sparse represented Hamiltonian and overlap to give the eigenValsvectors and
  !> values, as well as often the Cholesky factorized overlap matrix (due to a side effect of
  !> lapack)
  subroutine diagDenseRealMtxBlacs(electronicSolver, jobz, desc, HSqr, SSqr, eigenVals, eigenVecs)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> type of eigenVals-problem, either 'V'/'v' with vectors or 'N'/'n' eigenValsvalues only
    character, intent(in) :: jobz

    !> Dense descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Large square matrix with hamiltonian
    real(dp), intent(inout) :: HSqr(:,:)

    !> Large square matrix for the overlap workspace, often overwritten with the Cholesky factorized
    !> form.
    real(dp), intent(inout) :: SSqr(:,:)

    !> Eigenvalues.
    real(dp), intent(out) :: eigenVals(:)

    !> Eigenvectors
    real(dp), intent(out) :: eigenVecs(:,:)

    real(dp), save, allocatable :: SSqrTmp(:,:)

    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

    if (electronicSolver%tCholeskiiDecomposed(1)) then
      ! Matrix was already decomposed, so use stored copy (possibly 2 spin channels, but same
      ! overlap matrix)
    else
      ! if this is the first call, allocate space
      if (.not.allocated(ssqrTmp)) then
        allocate(ssqrTmp(size(ssqr,dim=1), size(ssqr,dim=2)))
      end if
      SSqrTmp = SSqr
    end if

    select case(electronicSolver%iSolver)
    case(electronicSolverTypes%QR)
      call scalafx_psygv(HSqr, desc, SSqrTmp, desc, eigenVals, eigenVecs, desc, uplo="L",&
          & jobz=jobz, skipchol=electronicSolver%tCholeskiiDecomposed(1))
    case(electronicSolverTypes%divideandconquer)
      call scalafx_psygvd(HSqr, desc, SSqrTmp, desc, eigenVals, eigenVecs, desc, uplo="L",&
          & jobz="V", allocfix=.true., skipchol=electronicSolver%tCholeskiiDecomposed(1))
    case(electronicSolverTypes%relativelyrobust)
      call scalafx_psygvr(HSqr, desc, SSqrTmp, desc, eigenVals, eigenVecs, desc, uplo="L",&
          & jobz="V", skipchol=electronicSolver%tCholeskiiDecomposed(1))
    case(electronicSolverTypes%elpa)
    #:if WITH_ELSI
      ! ELPA solver, returns eigenstates
      ! note, this only factorises overlap on first call - no skipchol equivalent
      call elsi_ev_real(electronicSolver%elsiHandle, HSqr, SSqrTmp, eigenVals, eigenVecs)
    #:else
      call error('Requires the ELSI library')
    #:endif
    case default
      call error('Unknown eigensolver')
    end select

    ! Choleskii factors of S now definitely avalable
    electronicSolver%tCholeskiiDecomposed = .true.

  end subroutine diagDenseRealMtxBlacs



  !> Diagonalizes a sparse represented Hamiltonian and overlap to give the eigenValsvectors and
  !> values, as well as often the Cholesky factorized overlap matrix (due to a side effect of
  !> lapack)
  subroutine diagDenseCplxMtxBlacs(electronicSolver, parallelKS, iKS, jobz, desc, HSqr, SSqr,&
      & eigenVals, eigenVecs)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> K-points and spins in use
    type(TParallelKS), intent(in) :: parallelKS

    !> Index for current k-point/spin
    integer, intent(in) :: iKS

    !> type of eigenVals-problem, either 'V'/'v' with vectors or 'N'/'n' eigenValsvalues only
    character, intent(in) :: jobz

    !> Dense descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Large square matrix for the resulting eigenValsvectors
    complex(dp), intent(inout) :: HSqr(:,:)

    !> Large square matrix for the overlap workspace, often overwritten with the Cholesky factorized
    !> form.
    complex(dp), intent(inout) :: SSqr(:,:)

    !> Eigenvalues.
    real(dp), intent(out) :: eigenVals(:)

    !> Eigenvectors
    complex(dp), intent(out) :: eigenVecs(:,:)

    complex(dp), save, allocatable :: SSqrTmp(:,:,:)

    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

    ! if this is the first call, allocate space
    if (.not.allocated(ssqrTmp)) then
      allocate(ssqrTmp(size(ssqr,dim=1), size(ssqr,dim=2), parallelKS%nLocalKS))
    end if

    if (electronicSolver%tCholeskiiDecomposed(iKS)) then
      ! Matrix was already decomposed, so use stored copy (possibly 2 spin channels, but can not
      ! easily use the same overlap as may not be in same group or ordered in any particular
      ! pattern)
      SSqr = SSqrTmp(:,:,iKS)
    end if

    select case(electronicSolver%iSolver)
    case(electronicSolverTypes%QR)
      call scalafx_phegv(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz, &
          & skipchol=electronicSolver%tCholeskiiDecomposed(iKS))
    case(electronicSolverTypes%divideandconquer)
      call scalafx_phegvd(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz,&
          & allocfix=.true., skipchol=electronicSolver%tCholeskiiDecomposed(iKS))
    case(electronicSolverTypes%relativelyrobust)
      call scalafx_phegvr(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz,&
          & skipchol=electronicSolver%tCholeskiiDecomposed(iKS))
    case(electronicSolverTypes%elpa)
  #:if WITH_ELSI
      ! ELPA solver, returns eigenstates
      ! note, this only factorises overlap on first call - no skipchol equivalent
      call elsi_ev_complex(electronicSolver%elsiHandle, HSqr, SSqr, eigenVals, eigenVecs)
  #:else
      call error('Requires the ELSI library')
  #:endif
    case default
      call error('Unknown eigensolver')
    end select

    if (.not.electronicSolver%tCholeskiiDecomposed(iKS)) then
      ! store the factorised overlap matrix for this iKS value
      SSqrTmp(:,:,iKS) = SSqr
    end if

    ! Choleskii factors of S now avalable
    electronicSolver%tCholeskiiDecomposed(iKS) = .true.

  end subroutine diagDenseCplxMtxBlacs

#:endif

end module eigenvects
