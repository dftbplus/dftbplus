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
#:if WITH_SCALAPACK
  use dftbp_extlibs_scalapackfx, only : DLEN_, CSRC_, RSRC_, MB_, NB_, pblasfx_ptranc,&
      & pblasfx_ptran, scalafx_indxl2g
#:endif

  implicit none

  private
  public :: adjointLowerTriangle
#:if WITH_SCALAPACK
  public :: adjointLowerTriangle_BLACS

  !> Copy lower triangle of distributed matrix into the upper triangle, obeying hermitian symmetry
  !! if appropriate
  interface adjointLowerTriangle_BLACS
    module procedure symmetrize_BLACS
    module procedure hermitian_BLACS
  end interface adjointLowerTriangle_BLACS

#:endif

  !> Copy lower triangle into the upper triangle of a square matrix, obeying hermitian symmetry if
  !! appropriate
  interface adjointLowerTriangle
    module procedure symmetrizeSquareMatrix
    module procedure hermitianSquareMatrix
  end interface adjointLowerTriangle


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
          matrix(ii, jj) = work(ii, jj)
        elseif (iGlob == jGlob) then
          matrix(ii, jj) = real(matrix(ii, jj),dp)
        end if
      end do
    end do

  end subroutine hermitian_BLACS

#:endif

end module dftbp_math_matrixops
