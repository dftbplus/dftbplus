!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Simple matrix operations for which LAPACK does not have a direct call
module dftbp_math_matrixops
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_schedule, only : assembleChunks
#:if WITH_SCALAPACK
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip
  use dftbp_extlibs_scalapackfx, only : DLEN_, CSRC_, RSRC_, MB_, NB_, pblasfx_ptranc,&
      & pblasfx_ptran, pblasfx_psymv, scalafx_indxl2g
#:else
  use dftbp_math_blasroutines, only : symm, hemm
#:endif

  implicit none

  private
  public :: adjointLowerTriangle, orthonormalizeVectors, iterOrthonorm
  public :: generalOrthonormalizeVectors
#:if WITH_SCALAPACK
  public :: adjointLowerTriangle_BLACS
#:endif

  !> Copy lower triangle into the upper triangle of a square matrix, obeying hermitian symmetry if
  !! appropriate
  interface adjointLowerTriangle
    module procedure symmetrizeSquareMatrix
    module procedure hermitianSquareMatrix
  end interface adjointLowerTriangle

  !> Orthogonalise a specfied set of rows in a matrix
  interface orthonormalizeVectors
    module procedure realGramSchmidt
    module procedure cmplxGramSchmidt
  end interface orthonormalizeVectors

  interface iterOrthonorm
    module procedure iGS_real
    module procedure iGS_cmplx
  end interface iterOrthonorm

  interface generalOrthonormalizeVectors
    module procedure realGeneralGramSchmidt
    module procedure cmplxGeneralGramSchmidt
  end interface generalOrthonormalizeVectors

#:if WITH_SCALAPACK

  !> Copy lower triangle of distributed matrix into the upper triangle, obeying hermitian symmetry
  !! if appropriate
  interface adjointLowerTriangle_BLACS
    module procedure symmetrize_BLACS
    module procedure hermitian_BLACS
  end interface adjointLowerTriangle_BLACS

#:endif

contains


  !> Copy lower triangle to upper for a square matrix.
  subroutine symmetrizeSquareMatrix(matrix)

    !> matrix to symmetrize
    real(dp), intent(inout) :: matrix(:,:)
    integer :: ii, matSize

    matSize = size(matrix, dim = 1)
    do ii = 1, matSize -1
      matrix(ii, ii+1:) = matrix(ii+1:, ii)
    end do

  end subroutine symmetrizeSquareMatrix


  !> Copy lower triangle to upper for a square matrix with Hermitian symmetry
  subroutine hermitianSquareMatrix(matrix)

    !> matrix to symmetrize
    complex(dp), intent(inout) :: matrix(:,:)
    integer :: ii, matSize

    matSize = size(matrix, dim = 1)
    do ii = 1, matSize -1
      matrix(ii, ii+1:) = conjg(matrix(ii+1:, ii))
    end do

  end subroutine hermitianSquareMatrix


#:if WITH_SCALAPACK

  !> Copy upper triangle into lower triangle of distributed matrix
  subroutine symmetrize_BLACS(desc, myCol, myRow, nCol, nRow, matrix)

    !> BLACS matrix descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Column of the current process in the BLACS grid
    integer, intent(in) :: myCol

    !> Row of the current process in the BLACS grid
    integer, intent(in) :: myRow

    !> Number. of process columns in the grid
    integer, intent(in) :: nCol

    !> Number. of process rows in the grid
    integer, intent(in) :: nRow

    !> Matrix to symmetrize
    real(dp), intent(inout) :: matrix(:,:)

    real(dp), allocatable :: work(:,:)
    integer :: ii, jj, iGlob, jGlob

    allocate(work(size(matrix, dim=1), size(matrix, dim=2)), source = 0.0_dp)
    call pblasfx_ptran(matrix, desc, work, desc)
    do jj = 1, size(matrix,dim=2)
      jGlob = scalafx_indxl2g(jj, desc(NB_), mycol, desc(CSRC_), ncol)
      do ii = 1, size(matrix,dim=1)
        iGlob = scalafx_indxl2g(ii, desc(MB_), myrow, desc(RSRC_), nrow)
        if (iGlob < jGlob) then
          ! copy transposed work matrix into the lower triangle of the main matrix
          matrix(ii, jj) = work(ii, jj)
        end if
      end do
    end do

  end subroutine symmetrize_BLACS


  !> Copy upper triangle into lower triangle of distributed hermitian matrix
  subroutine hermitian_BLACS(desc, myCol, myRow, nCol, nRow, matrix)

    !> BLACS matrix descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Column of the current process in the BLACS grid
    integer, intent(in) :: myCol

    !> Row of the current process in the BLACS grid
    integer, intent(in) :: myRow

    !> Number. of process columns in the grid
    integer, intent(in) :: nCol

    !> Number. of process rows in the grid
    integer, intent(in) :: nRow

    !> Matrix to hermitian symmetrize
    complex(dp), intent(inout) :: matrix(:,:)

    complex(dp), allocatable :: work(:,:)
    integer :: ii, jj, iGlob, jGlob

    allocate(work(size(matrix, dim=1), size(matrix, dim=2)), source = (0.0_dp, 0.0_dp))
    call pblasfx_ptranc(matrix, desc, work, desc)
    do jj = 1, size(matrix,dim=2)
      jGlob = scalafx_indxl2g(jj, desc(NB_), mycol, desc(CSRC_), ncol)
      do ii = 1, size(matrix,dim=1)
        iGlob = scalafx_indxl2g(ii, desc(MB_), myrow, desc(RSRC_), nrow)
        if (iGlob < jGlob) then
          ! copy transposed work matrix into the lower triangle of the main matrix
          matrix(ii, jj) = work(ii, jj)
        elseif (iGlob == jGlob) then
          ! diagonal should be real for hermitian matrix
          matrix(ii, jj) = real(matrix(ii, jj),dp)
        end if
      end do
    end do

  end subroutine hermitian_BLACS

#:endif


  !> Perform modified Gram-Schmidt orthonormalization of the vectors in the columns of
  !! vecs(:,start:end), while also keeping them orthogonal to vecs(:,:start-1) (which are assumed to
  !! already be orthogonal)
  subroutine realGramSchmidt(env, iStart, iEnd, vecs, nullTol, orthogonaliseFrom)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Starting place in vectors to work from
    integer, intent(in) :: iStart

    !> Ending place in vectors
    integer, intent(in) :: iEnd

    !> Vectors, iStart:iEnd to be orthogonalized against 1:iEnd vectors
    real(dp), intent(inout) :: vecs(:,:)

    !> Optional tolerance to drop orthogonalised vectors (before normalisation)
    real(dp), intent(in), optional :: nullTol

    !> If we need to orthogonalize against fixed vectors after iEnd
    integer, intent(in), optional :: orthogonaliseFrom

    integer :: ii, jj
    real(dp) :: dummyReal, nullTol_

    @:ASSERT(iStart <= iEnd)

    if (present(nullTol)) then
      nullTol_ = nullTol
    else
      ! Used in a test against (supposedly) positive definite value
      nullTol_ = -1.0_dp
    end if

    ! Obviously, not optimal in terms of communication, can be optimized if necessary. Assumes vecs
    ! are column block distributed (not block cyclic column and row) if parallel.
    do ii = iStart, iEnd

      do jj = 1, ii - 1
        dummyReal = dot_product(vecs(:,jj), vecs(:,ii))
        call assembleChunks(env, dummyReal)
        vecs(:,ii) = vecs(:,ii) - dummyReal * vecs(:,jj)
      end do

      if (present(orthogonaliseFrom)) then
        @:ASSERT(orthogonaliseFrom > iEnd)
        do jj = orthogonaliseFrom, size(vecs, dim=2)
          dummyReal = dot_product(vecs(:,jj), vecs(:,ii))
          call assembleChunks(env, dummyReal)
          vecs(:,ii) = vecs(:,ii) - dummyReal * vecs(:,jj)
        end do
      end if

      ! Normalize the resulting vector.
      dummyReal = dot_product(vecs(:,ii), vecs(:,ii))
      call assembleChunks(env, dummyReal)

      if (dummyReal > nullTol_) then
        vecs(:,ii) = vecs(:,ii) / sqrt(dummyReal)
      else
        ! In the null space
        vecs(:,ii) = 0.0_dp
      end if

    end do

  end subroutine realGramSchmidt


  !> Perform modified Gram-Schmidt orthonormalization of the vectors in the columns of
  !! vecs(:,start:end), while also keeping them orthogonal to vecs(:,:start-1) (which are assumed to
  !! already be orthogonal)
  subroutine cmplxGramSchmidt(env, iStart, iEnd, vecs, nullTol, orthogonaliseFrom)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Starting place in vectors to work from
    integer, intent(in) :: iStart

    !> Ending place in vectors
    integer, intent(in) :: iEnd

    !> Vectors, iStart:iEnd to be orthogonalized against 1:iEnd vectors
    complex(dp), intent(inout) :: vecs(:,:)

    !> Optional tolerance to drop orthogonalised vectors (before normalisation)
    real(dp), intent(in), optional :: nullTol

    !> If we need to orthogonalize against fixed vectors after iEnd
    integer, intent(in), optional :: orthogonaliseFrom

    integer :: ii, jj
    real(dp) :: dummyReal, nullTol_
    complex(dp) :: dummyCmplx

    if (present(nullTol)) then
      nullTol_ = nullTol
    else
      ! Used in a test against (supposedly) positive definite value
      nullTol_ = -1.0_dp
    end if

    ! Obviously, not optimal in terms of communication, can be optimized if necessary. Assumes vecs
    ! are column distributed (not column and row distributed) if parallel.
    do ii = iStart, iEnd

      do jj = 1, ii - 1
        dummyCmplx = dot_product(vecs(:,jj), vecs(:,ii))
        call assembleChunks(env, dummyCmplx)
        vecs(:,ii) = vecs(:,ii) - dummyCmplx * vecs(:,jj)
      end do

      if (present(orthogonaliseFrom)) then
        @:ASSERT(orthogonaliseFrom > iEnd)
        do jj = orthogonaliseFrom, size(vecs, dim=2)
          dummyCmplx = dot_product(vecs(:,jj), vecs(:,ii))
          call assembleChunks(env, dummyCmplx)
          vecs(:,ii) = vecs(:,ii) - dummyCmplx * vecs(:,jj)
        end do
      end if

      ! Normalize the resulting vector.
      dummyReal = real(dot_product(vecs(:,ii), vecs(:,ii)),dp)
      call assembleChunks(env, dummyReal)

      if (dummyReal > nullTol_) then
        vecs(:,ii) = vecs(:,ii) / sqrt(dummyReal)
      else
        ! null space
        vecs(:,ii) = cmplx(0,0,dp)
      end if

    end do

  end subroutine cmplxGramSchmidt


  !> Perform modified Gram-Schmidt orthonormalization of the vectors in the columns of
  !! vecs(:,start:end), while also keeping them orthogonal to vecs(:,:start-1) (which are assumed to
  !! already be orthogonal)
  subroutine realGeneralGramSchmidt(env, iStart, iEnd, vecs, sVecs, nullTol, orthogonaliseFrom)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Starting place in vectors to work from
    integer, intent(in) :: iStart

    !> Ending place in vectors
    integer, intent(in) :: iEnd

    !> Vectors, iStart:iEnd to be orthogonalized against 1:iEnd vectors
    real(dp), intent(inout) :: vecs(:,:)

    !> Overlap times vectors, iStart:iEnd to be orthogonalized against 1:iEnd vectors
    real(dp), intent(inout) :: sVecs(:,:)

    !> Optional tolerance to drop orthogonalised vectors (before normalisation)
    real(dp), intent(in), optional :: nullTol

    !> If we need to orthogonalize against fixed vectors after iEnd
    integer, intent(in), optional :: orthogonaliseFrom

    integer :: ii, jj
    real(dp) :: dummyReal, nullTol_

    @:ASSERT(iStart <= iEnd)

    if (present(nullTol)) then
      nullTol_ = nullTol
    else
      ! Used in a test against (supposedly) positive definite value
      nullTol_ = -1.0_dp
    end if

    ! Obviously, not optimal in terms of communication, can be optimized if necessary. Assumes vecs
    ! are column block distributed (not block cyclic column and row) if parallel.
    do ii = iStart, iEnd

      do jj = 1, ii - 1
        dummyReal = dot_product(vecs(:,jj), sVecs(:,ii))
        call assembleChunks(env, dummyReal)
        vecs(:,ii) = vecs(:,ii) - dummyReal * vecs(:,jj)
        sVecs(:,ii) = sVecs(:,ii) - dummyReal * sVecs(:, jj)
      end do

      if (present(orthogonaliseFrom)) then
        @:ASSERT(orthogonaliseFrom > iEnd)
        do jj = orthogonaliseFrom, size(vecs, dim=2)
          dummyReal = dot_product(vecs(:,jj), sVecs(:,ii))
          call assembleChunks(env, dummyReal)
          vecs(:,ii) = vecs(:,ii) - dummyReal * vecs(:,jj)
          sVecs(:,ii) = sVecs(:,ii) - dummyReal * sVecs(:,jj)
        end do
      end if

      ! Normalize the resulting vector.
      dummyReal = dot_product(vecs(:,ii), sVecs(:,ii))
      call assembleChunks(env, dummyReal)

      if (dummyReal > nullTol_) then
        vecs(:,ii) = vecs(:,ii) / sqrt(dummyReal)
        sVecs(:,ii) = sVecs(:,ii) / sqrt(dummyReal)
      else
        ! In the null space
        vecs(:,ii) = 0.0_dp
        sVecs(:,ii) = 0.0_dp
      end if

    end do

  end subroutine realGeneralGramSchmidt


  !> Perform modified Gram-Schmidt orthonormalization of the vectors in the columns of
  !! vecs(:,start:end), while also keeping them orthogonal to vecs(:,:start-1) (which are assumed to
  !! already be orthogonal)
  subroutine cmplxGeneralGramSchmidt(env, iStart, iEnd, vecs, sVecs, nullTol, orthogonaliseFrom)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Starting place in vectors to work from
    integer, intent(in) :: iStart

    !> Ending place in vectors
    integer, intent(in) :: iEnd

    !> Vectors, iStart:iEnd to be orthogonalized against 1:iEnd vectors
    complex(dp), intent(inout) :: vecs(:,:)

    !> Overlap times vectors, iStart:iEnd to be orthogonalized against 1:iEnd vectors
    complex(dp), intent(inout) :: sVecs(:,:)

    !> Optional tolerance to drop orthogonalised vectors (before normalisation)
    real(dp), intent(in), optional :: nullTol

    !> If we need to orthogonalize against fixed vectors after iEnd
    integer, intent(in), optional :: orthogonaliseFrom

    integer :: ii, jj
    real(dp) :: dummyReal, nullTol_
    complex(dp) :: dummyCmplx

    if (present(nullTol)) then
      nullTol_ = nullTol
    else
      ! Used in a test against (supposedly) positive definite value
      nullTol_ = -1.0_dp
    end if

    ! Obviously, not optimal in terms of communication, can be optimized if necessary. Assumes vecs
    ! are column distributed (not column and row distributed) if parallel.
    do ii = iStart, iEnd

      do jj = 1, ii - 1
        dummyCmplx = dot_product(vecs(:,jj), sVecs(:,ii))
        call assembleChunks(env, dummyCmplx)
        vecs(:,ii) = vecs(:,ii) - dummyCmplx * vecs(:,jj)
        sVecs(:,ii) = sVecs(:,ii) - dummyCmplx * sVecs(:,jj)
      end do

      if (present(orthogonaliseFrom)) then
        @:ASSERT(orthogonaliseFrom > iEnd)
        do jj = orthogonaliseFrom, size(vecs, dim=2)
          dummyCmplx = dot_product(vecs(:,jj), sVecs(:,ii))
          call assembleChunks(env, dummyCmplx)
          vecs(:,ii) = vecs(:,ii) - dummyCmplx * vecs(:,jj)
          sVecs(:,ii) = sVecs(:,ii) - dummyCmplx * sVecs(:,jj)
        end do
      end if

      ! Normalize the resulting vector.
      dummyReal = real(dot_product(vecs(:,ii), sVecs(:,ii)),dp)
      call assembleChunks(env, dummyReal)

      if (dummyReal > nullTol_) then
        vecs(:,ii) = vecs(:,ii) / sqrt(dummyReal)
        sVecs(:,ii) = sVecs(:,ii) / sqrt(dummyReal)
      else
        ! null space
        vecs(:,ii) = cmplx(0,0,dp)
        sVecs(:,ii) = cmplx(0,0,dp)
      end if

    end do

  end subroutine cmplxGeneralGramSchmidt


  !> Performs (unmodified) iterative Gram-Schmidt orthonormalization
  !! Lingen, F.J. (2000), Efficient Gram–Schmidt orthonormalisation on parallel
  !! computers. Commun. Numer. Meth. Engng., 16: 57-66.
  !! https://doi.org/10.1002/(SICI)1099-0887(200001)16:1%3C57::AID-CNM320%3E3.0.CO;2-I
  subroutine iGS_real(env, iStart, iEnd, vecs)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Starting place in vectors to work from
    integer, intent(in) :: iStart

    !> Ending place in vectors
    integer, intent(in) :: iEnd

    !> Vectors to be orthogonalized against 1:iEnd vectors
    real(dp), intent(inout) :: vecs(:,:)

    integer :: ii, iNew, iOld, m, n, kk
    real(dp), allocatable :: q(:,:), p(:)
    real(dp), parameter :: alpha = 0.5_dp ! value from paper
    real(dp) :: norm(0:1)

  #:if WITH_SCALAPACK

  #:else

    m = size(vecs, dim=1)
    n = iEnd
    allocate(p(n))
    allocate(q(m, 0:1))

    if (iStart == 1) then
      ! nothing before this to orthogonalise against, so just try to normalise
      norm(0) = sqrt(dot_product(vecs(:, iStart), vecs(:, iStart)))
      if (norm(0) < epsilon(1.0_dp)) then
        vecs(:, iStart) = 0.0_dp
      else
        vecs(:, iStart) = vecs(:, iStart) / norm(0)
      end if
    end if

    do ii = max(iStart, 2), iEnd

      iOld = 0
      iNew = 1
      q(:, iOld) = vecs(:, ii)
      norm(iOld) = sqrt(dot_product(q(:, iOld), q(:, iOld)))
      kk = 0
      lpIter: do while (.true.)
        kk = kk + 1
        p(:ii-1) = matmul(transpose(vecs(:, :ii-1)), q(:, iOld))
        q(:,iNew) = q(:, iOld) - matmul(vecs(:, :ii-1), p(:ii-1))

        norm(iNew) = sqrt(dot_product(q(:, iNew), q(:, iNew)))
        if (norm(iNew) < epsilon(1.0_dp)) then
          q(:,iNew) = 0.0_dp
          norm(iNew) = 1.0_dp
          exit lpIter
        end if
        if (norm(iNew) > alpha * norm(iOld) ) then
          ! L2 norm test passed
          exit lpIter
        end if

        iNew = mod(iNew + 1, 2)
        iOld = mod(iOld + 1, 2)

      end do lpIter

      vecs(:,ii) = q(:,iNew) / norm(iNew)

    end do

  #:endif

  end subroutine iGS_real


  !> Performs (unmodified) iterative Gram-Schmidt orthonormalization
  !! Lingen, F.J. (2000), Efficient Gram–Schmidt orthonormalisation on parallel
  !! computers. Commun. Numer. Meth. Engng., 16: 57-66.
  !! https://doi.org/10.1002/(SICI)1099-0887(200001)16:1%3C57::AID-CNM320%3E3.0.CO;2-I
  subroutine iGS_cmplx(env, iStart, iEnd, vecs)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Starting place in vectors to work from
    integer, intent(in) :: iStart

    !> Ending place in vectors
    integer, intent(in) :: iEnd

    !> Vectors to be orthogonalized against 1:iEnd vectors
    complex(dp), intent(inout) :: vecs(:,:)

    integer :: ii, iNew, iOld, m, n, kk
    complex(dp), allocatable :: q(:,:), p(:)
    real(dp), parameter :: alpha = 0.5_dp ! value from paper
    real(dp) :: norm(0:1)

  #:if WITH_SCALAPACK

  #:else

    m = size(vecs, dim=1)
    n = iEnd
    allocate(p(n))
    allocate(q(m, 0:1))

    if (iStart == 1) then
      ! nothing before this to orthogonalise against, so just try to normalise
      norm(0) = sqrt(real(dot_product(vecs(:, iStart), vecs(:, iStart)), dp))
      if (norm(0) < epsilon(1.0_dp)) then
        vecs(:, iStart) = cmplx(0,0,dp)
      else
        vecs(:, iStart) = vecs(:, iStart) / norm(0)
      end if
    end if

    do ii = max(iStart, 2), iEnd

      iOld = 0
      iNew = 1
      q(:, iOld) = vecs(:, ii)
      norm(iOld) = sqrt(real(dot_product(q(:, iOld), q(:, iOld)), dp))
      kk = 0
      lpIter: do while (.true.)
        kk = kk + 1
        p(:ii-1) = matmul(transpose(conjg(vecs(:, :ii-1))), q(:, iOld))
        q(:,iNew) = q(:, iOld) - matmul(vecs(:, :ii-1), p(:ii-1))

        norm(iNew) = sqrt(real(dot_product(q(:, iNew), q(:, iNew)), dp))
        if (norm(iNew) > alpha * norm(iOld) ) then
          ! L2 norm test passed
          exit lpIter
        end if

        iNew = mod(iNew + 1, 2)
        iOld = mod(iOld + 1, 2)

      end do lpIter

      vecs(:,ii) = q(:,iNew) / norm(iNew)

    end do

  #:endif

  end subroutine iGS_cmplx

end module dftbp_math_matrixops
