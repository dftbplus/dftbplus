!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains some miscellaneous eigenvector phase related routines
module dftbp_math_phase
  use dftbp_common_accuracy, only : dp, elecTolMax
  use dftbp_common_environment, only : TEnvironment
  use dftbp_math_degeneracy, only : TDegeneracyFind
  use dftbp_math_blasroutines, only : symm
  use dftbp_math_matrixops, only : orthonormalizeVectors
#:if WITH_SCALAPACK
  use dftbp_common_blacsenv, only : TBlacsEnv
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip
  use dftbp_extlibs_scalapackfx, only : blacsgrid, blacsfx_gemr2d, scalafx_indxl2g,&
      & scalafx_getlocalshape, DLEN_, CSRC_, NB_
  use dftbp_type_densedescr, only : TDenseDescr
#:endif

  implicit none

  private
  public :: phaseLock

  !> Enforce (an arbitrary) phase convention on eigenvectors, including for denerate states
  interface phaseLock
    module procedure phaseLock_real
  #:if WITH_SCALAPACK
    module procedure phaseLock_real_BLACS
  #:endif
    module procedure phaseLockGeneral_real
  end interface phaseLock

contains


  !> Apply a phase convention to a set of eigenvectors, including defined unitary choice for
  !! degenerate sets of vectors.
  subroutine phaseLock_real(env, evec, eval, tolerance)

    !> Computatonal environment settings
    type(TEnvironment), intent(in) :: env

    !> Eigenvectors
    real(dp), intent(inout) :: evec(:,:)

    !> Eigenvalues
    real(dp), intent(in) :: eval(:)

    !> Optional tolerance
    real(dp), intent(in), optional :: tolerance

    integer :: ii, jj, iDegen, nDegen,  maxDegen, nDegenVectors
    integer, allocatable :: degenerate(:,:) ! (start & end, n degenerate sub-spaces)
    logical :: isDegenerate
    real(dp) :: tol
    real(dp), parameter :: tolMax = 100.0_dp * epsilon(1.0_dp)
    real(dp), allocatable :: vTrial(:, :)
    type(TDegeneracyFind) :: DegeneracyFind

    tol = elecTolMax
    if (present(tolerance)) tol = tolerance

    ! set phase of each row of the matrix
    call simplePhase(eVec, tol)

    call degeneracyFind%init(tol)
    call degeneracyFind%degeneracyTest(eval, isDegenerate)

    ! Harmless to use the above simplePhase transformation on degenerate subspaces, even though that
    ! still leaves them with an undetermined unitary choice, so wastes a small amount of time.
    if (.not. isDegenerate) then

      return

    endif

    degenerate = degeneracyFind%degenerateRanges()
    maxDegen = maxval(degenerate(2,:) - degenerate(1,:)) + 1

    allocate(vTrial(size(evec, dim=1), maxDegen), source=0.0_dp)

    lpDegen : do iDegen = 1, size(degenerate, dim=2)

      nDegen = degenerate(2,iDegen) - degenerate(1,iDegen) + 1

      if (nDegen < 2) cycle lpDegen

      nDegenVectors = 0

      lpStartIndex : do ii = 1, size(evec, dim=1)

        if (sum(abs(evec(ii,degenerate(1,iDegen):degenerate(2,iDegen)))) <= tolMax) then

          ! Zero projection (L1) of +ve ii-th cartesian unit vector onto the subspace, not a
          ! suitable candidate vector.
          cycle lpStartIndex

        else

          vTrial(:,nDegenVectors+1) = 0.0_dp
          do jj = degenerate(1,iDegen), degenerate(2,iDegen)
            ! project the defined unit vector onto the orthonormal vectors of the subspace
            vTrial(:,nDegenVectors+1) = vTrial(:,nDegenVectors+1) + evec(:, jj) * evec(ii, jj)
          end do
          vTrial(:,nDegenVectors+1) = vTrial(:,nDegenVectors+1)&
              & / sqrt(sum(vTrial(:,nDegenVectors+1)**2))

          if (nDegenVectors == 0) then
            ! first vector in the subspace
            nDegenVectors = 1

          else

            call orthonormalizeVectors(env, nDegenVectors+1, nDegenVectors+1,&
                & vTrial(:,:nDegenVectors+1), nullTol = 1000.0_dp * epsilon(0.0_dp))

            if (sum(vTrial(:,nDegenVectors+1)**2) > 0.0_dp) then
              ! Orthonormal state is inside the subspace
              nDegenVectors = nDegenVectors + 1
            end if

          end if

        end if

        if (nDegenVectors == nDegen) exit lpStartIndex ! Enough vectors have been found

      end do lpStartIndex

      eVec(:, degenerate(1,iDegen):degenerate(2,iDegen)) = vTrial(:,:nDegenVectors)

      ! ensure the vectors spanning the subspace are orthogonal to all of the other eigenvectors
      call orthonormalizeVectors(env, degenerate(1,iDegen), degenerate(2,iDegen), eVec,&
          & orthogonaliseFrom = degenerate(2,iDegen) + 1)

    end do lpDegen

  end subroutine phaseLock_real


  !> Apply a phase convention to a set of generalised eigenvectors, including defined unitary
  !! choice for degenerate sets of vectors.
  subroutine phaseLockGeneral_real(env, evec, overlap, eval, tolerance)

    !> Computatonal environment settings
    type(TEnvironment), intent(in) :: env

    !> Eigenvectors
    real(dp), intent(inout) :: evec(:,:)

    !> Overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    !> Eigenvalues
    real(dp), intent(in) :: eval(:)

    !> Optional tolerance
    real(dp), intent(in), optional :: tolerance

    integer :: ii, jj, iDegen, nDegen,  maxDegen, nDegenVectors
    integer, allocatable :: degenerate(:,:) ! (start & end, n degenerate sub-spaces)
    logical :: isDegenerate
    real(dp) :: tol
    real(dp), parameter :: tolMax = 100.0_dp * epsilon(1.0_dp)
    real(dp), allocatable :: vTrial(:, :), sVTrial(:,:), sEvec(:,:)
    type(TDegeneracyFind) :: DegeneracyFind

    tol = elecTolMax
    if (present(tolerance)) tol = tolerance

    ! set phase of each row of the matrix
    call simplePhase(eVec, tol)

    call degeneracyFind%init(tol)
    call degeneracyFind%degeneracyTest(eval, isDegenerate)

    ! Harmless to use the above simplePhase transformation on degenerate subspaces, even though that
    ! still leaves them with an undetermined unitary choice, so wastes a small amount of time.
    if (.not. isDegenerate) then

      return

    endif

    degenerate = degeneracyFind%degenerateRanges()
    maxDegen = maxval(degenerate(2,:) - degenerate(1,:)) + 1

    allocate(vTrial(size(evec, dim=1), maxDegen), source=0.0_dp)
    allocate(sVTrial(size(evec, dim=1), maxDegen), source=0.0_dp)

    lpDegen : do iDegen = 1, size(degenerate, dim=2)

      nDegen = degenerate(2,iDegen) - degenerate(1,iDegen) + 1

      if (nDegen < 2) cycle lpDegen

      allocate(sEvec(size(eVec, dim=1), size(eVec, dim=2)))
      call symm(sEvec, 'L', overlap, eVec)

      nDegenVectors = 0

      lpStartIndex : do ii = 1, size(evec, dim=1)

        vTrial(:,nDegenVectors+1) = 0.0_dp
        do jj = degenerate(1,iDegen), degenerate(2,iDegen)
          ! project the defined unit vector onto the orthonormal vectors of the subspace
          vTrial(:,nDegenVectors+1) = vTrial(:,nDegenVectors+1) + evec(:, jj) * evec(ii, jj)
        end do

        vTrial(:,nDegenVectors+1) = vTrial(:,nDegenVectors+1)&
            & / sqrt(sum(vTrial(:,nDegenVectors+1)**2))

        if (nDegenVectors == 0) then
          ! first vector in the subspace
          nDegenVectors = 1

        else

          call orthonormalizeVectors(env, nDegenVectors+1, nDegenVectors+1,&
              & vTrial(:,:nDegenVectors+1), nullTol = 1000.0_dp * epsilon(0.0_dp))

          if (sum(vTrial(:,nDegenVectors+1)**2) > 0.0_dp) then
            ! Orthonormal state is inside the subspace
            nDegenVectors = nDegenVectors + 1
          end if

        end if


        if (nDegenVectors == nDegen) exit lpStartIndex ! Enough vectors have been found

      end do lpStartIndex

      eVec(:, degenerate(1,iDegen):degenerate(2,iDegen)) = vTrial(:,:nDegenVectors)

      ! ensure the vectors spanning the subspace are orthogonal to all of the other eigenvectors
      call orthonormalizeVectors(env, degenerate(1,iDegen), degenerate(2,iDegen), eVec,&
          & orthogonaliseFrom = degenerate(2,iDegen) + 1)

    end do lpDegen

  end subroutine phaseLockGeneral_real

#:if WITH_SCALAPACK

  !> Apply a phase convention to a set of eigenvectors, including defined unitary choice for
  !! degenerate sets of vectors.  Note : need to also include generalised eigenvectors (with an
  !! overlap matrix).
  subroutine phaseLock_real_BLACS(env, evec, eval, tolerance, denseDesc)

    !> Computatonal environment settings
    type(TEnvironment), intent(in) :: env

    !> Eigenvectors
    real(dp), intent(inout) :: evec(:,:)

    !> Eigenvalues
    real(dp), intent(in) :: eval(:)

    !> Optional tolerance
    real(dp), intent(in), optional :: tolerance

    !> Descriptor of the dense matrix
    type(TDenseDescr), intent(in) :: denseDesc

    real(dp) :: tmpSum
    integer :: jGlob, nLocalRows, nLocalCols
    real(dp), allocatable :: evecComm(:)
    ! Work matrix for entire rows of the matrix on a given rank
    real(dp), allocatable :: workRowMatrix(:,:)
    ! Work matrix for entire columns of the matrix on a given rank
    real(dp), allocatable :: workColMatrix(:,:)

    integer :: ii, jj, iDegen, nDegen,  maxDegen, nDegenVectors
    integer, allocatable :: degenerate(:,:) ! (start & end, n degenerate sub-spaces)
    logical :: isDegenerate
    real(dp) :: tol
    real(dp), parameter :: tolMax = 100.0_dp * epsilon(1.0_dp)
    real(dp), allocatable :: vTrial(:, :)
    type(TDegeneracyFind) :: DegeneracyFind

    tol = elecTolMax
    if (present(tolerance)) tol = tolerance

    call scalafx_getlocalshape(env%blacs%rowOrbitalGrid, denseDesc%blacsRowSqr, nLocalRows,&
        & nLocalCols)
    allocate(workRowMatrix(nLocalRows, nLocalCols), source=0.0_dp)
    ! redistribute global matrix to row form
    call blacsfx_gemr2d(denseDesc%nOrb, denseDesc%nOrb, eVec, 1, 1, denseDesc%blacsOrbSqr,&
        & workRowMatrix, 1, 1, denseDesc%blacsRowSqr, env%blacs%orbitalGrid%ctxt)

    ! set phase of each row of the matrix
    call simplePhase(workRowMatrix, tol)

    call degeneracyFind%init(tol)
    call degeneracyFind%degeneracyTest(eval, isDegenerate)

    ! Harmless to use the above simplePhase transformation on degenerate subspaces, even though that
    ! still leaves them with an undetermined unitary choice, so wastes a small amount of time.
    if (.not. isDegenerate) then

      ! redistribute processed row matrix back to block cyclic form. No degenerate cases processed,
      ! so matrices still in row storage
      call blacsfx_gemr2d(denseDesc%nOrb, denseDesc%nOrb, workRowMatrix, 1, 1,&
          & denseDesc%blacsRowSqr, eVec, 1, 1, denseDesc%blacsOrbSqr, env%blacs%orbitalGrid%ctxt)

      return

    endif

    degenerate = degeneracyFind%degenerateRanges()
    maxDegen = maxval(degenerate(2,:) - degenerate(1,:)) + 1

    ! convert row stored matrix to column storage
    call scalafx_getlocalshape(env%blacs%colOrbitalGrid, denseDesc%blacsColSqr, nLocalRows,&
        & nLocalCols)
    allocate(workColMatrix(nLocalRows, nLocalCols), source=0.0_dp)
    call blacsfx_gemr2d(denseDesc%nOrb, denseDesc%nOrb, workRowMatrix, 1, 1, denseDesc%blacsRowSqr,&
        & workColMatrix, 1, 1, denseDesc%blacsColSqr, env%blacs%orbitalGrid%ctxt)

    allocate(vTrial(size(workColMatrix, dim=1), maxDegen), source=0.0_dp)

    lpMpiDegen : do iDegen = 1, size(degenerate, dim=2)

      nDegen = degenerate(2,iDegen) - degenerate(1,iDegen) + 1

      if (nDegen < 2) cycle lpMpiDegen

      nDegenVectors = 0

      lpStartIndexMpi : do ii = 1, denseDesc%nOrb

        allocate(evecComm(degenerate(1,iDegen):degenerate(2,iDegen)), source=0.0_dp)
        tmpSum = 0.0_dp

        do jj = 1, size(workRowMatrix, dim=2)
          jGlob = scalafx_indxl2g(jj, denseDesc%blacsRowSqr(NB_), env%blacs%rowOrbitalGrid%mycol,&
              & denseDesc%blacsRowSqr(CSRC_), env%blacs%rowOrbitalGrid%ncol)
          if (jGlob >= degenerate(1,iDegen) .and. jGlob <= degenerate(2,iDegen)) then
            tmpSum = tmpSum + abs(workRowMatrix(ii, jj))
            evecComm(jGlob) = workRowMatrix(ii, jj)
          end if
        end do
        call mpifx_allreduceip(env%mpi%groupComm, evecComm, MPI_SUM)

        if (tmpSum <= tolMax) then
          ! Zero projection (L1) of +ve ii-th cartesian unit vector onto the subspace, not a
          ! suitable candidate vector.
          cycle lpStartIndexMpi

        else

          call mpifx_allreduceip(env%mpi%groupComm, evecComm, MPI_SUM)
          vTrial(:,nDegenVectors+1) = 0.0_dp
          do jj = degenerate(1,iDegen), degenerate(2,iDegen)
            ! project the defined unit vector onto the orthonormal vectors of the subspace
            vTrial(:,nDegenVectors+1) = vTrial(:,nDegenVectors+1)&
                & + workColMatrix(:, jj) * evecComm(jj)
          end do

          ! Normalize the resulting vector.
          tmpSum = dot_product(vTrial(:,nDegenVectors+1), vTrial(:,nDegenVectors+1))
          call mpifx_allreduceip(env%mpi%groupComm, tmpSum, MPI_SUM)
          vTrial(:,nDegenVectors+1) = vTrial(:,nDegenVectors+1) / sqrt(tmpSum)

          if (nDegenVectors == 0) then
            ! first vector in the subspace
            nDegenVectors = 1
          else

            call orthonormalizeVectors(env, nDegenVectors+1, nDegenVectors+1,&
                & vTrial(:,:nDegenVectors+1), nullTol = 1000.0_dp * epsilon(0.0_dp))

            tmpSum = dot_product(vTrial(:,nDegenVectors+1), vTrial(:,nDegenVectors+1))
            call mpifx_allreduceip(env%mpi%groupComm, tmpSum, MPI_SUM)
            if (tmpSum > 0.0_dp) then
              ! Orthonormal state is inside the subspace
              nDegenVectors = nDegenVectors + 1
            end if

          end if

        end if

        if (nDegenVectors == nDegen) exit lpStartIndexMpi ! Enough vectors have been found

      end do lpStartIndexMpi

      workColMatrix(:, degenerate(1,iDegen):degenerate(2,iDegen)) = vTrial(:,:nDegenVectors)

      ! ensure the vectors spanning the subspace are orthogonal to all of the other eigenvectors
      call orthonormalizeVectors(env, degenerate(1,iDegen), degenerate(2,iDegen), workColMatrix,&
          & orthogonaliseFrom = degenerate(2,iDegen) + 1)

    end do lpMpiDegen

    ! redistribute processed column matrix back to block cyclic form
    call blacsfx_gemr2d(denseDesc%nOrb, denseDesc%nOrb, workColMatrix, 1, 1,&
        & denseDesc%blacsColSqr, eVec, 1, 1, denseDesc%blacsOrbSqr, env%blacs%orbitalGrid%ctxt)


  end subroutine phaseLock_real_BLACS

#:endif

  !> Simple phase convention for eigenvectors to make first non-negligible element positive for each
  !! row.
  subroutine simplePhase(matrix, tol)

    !> Matrix to process
    real(dp), intent(inout) :: matrix(:,:)

    !> Tolerance for testing non-zero elements
    real(dp), intent(in) :: tol

    integer :: iLev, nLev, iOrb

    nLev = size(matrix, dim=2)

    do iLev = 1, nLev
      lpPhase : do iOrb = 1, size(matrix, dim=1)
        if (abs(matrix(iOrb, iLev)) >= tol) then ! test needs work
          matrix(:, iLev) = sign(1.0_dp, matrix(iOrb, iLev)) * matrix(:, iLev)
          exit lpPhase
        end if
      end do lpPhase
    end do

  end subroutine simplePhase

end module dftbp_math_phase
