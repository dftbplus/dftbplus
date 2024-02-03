!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"
#:include "error.fypp"

!> Common mathematical operations built out of multiple scalapack calls
module dftbp_math_scalafxext
  use dftbp_common_accuracy, only : lc, dp
  use dftbp_io_message, only : error
#:if WITH_SCALAPACK
  use dftbp_common_status, only : TStatus
  use dftbp_extlibs_mpifx, only : mpifx_allreduceip, mpifx_comm, MPI_SUM
  use dftbp_extlibs_scalapackfx, only : MB_, NB_, RSRC_, CSRC_, DLEN_, blacsfx_gsum,&
      & scalafx_ppotrf, scalafx_ppotri, blacsgrid, scalafx_indxl2g, scalafx_getlocalshape
#:endif
  implicit none

  private
  public :: psymmatinv, phermatinv
#:if WITH_SCALAPACK
  public :: distrib2replicated
#:endif

contains

#:if WITH_SCALAPACK

  !> Inversion of a symmetric matrix
  subroutine psymmatinv(desc, aa, status, uplo)

    !> Matrix descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Matrix to invert on entry, inverted matrix on exit
    real(dp), intent(inout) :: aa(:,:)

    !> Status flag
    type(TStatus), intent(out) :: status

    !> Whether upper or lower triangle is specified in the matrix ("U" or "L", default: "L")
    character, intent(in), optional :: uplo

    integer :: info

    call scalafx_ppotrf(aa, desc, uplo=uplo, info=info)
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(status, -1, "('scalafx_ppotrf failed in psymmatinv (info: ',I0,')')",&
          & info)
    end if
    call scalafx_ppotri(aa, desc, uplo=uplo, info=info)
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(status, -1, "('scalafx_ppotri failed in psymmatinv (info: ',I0,')')",&
          & info)
    end if

  end subroutine psymmatinv


  !> Inversion of a hermitian matrix
  subroutine phermatinv(desc, aa, status, uplo)

    !> Matrix descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Matrix to invert on entry, inverted matrix on exit
    complex(dp), intent(inout) :: aa(:,:)

    !> Status flag
    type(TStatus), intent(out) :: status

    !> Whether upper or lower triangle is specified in the matrix ("U" or "L", default: "L")
    character, intent(in), optional :: uplo

    integer :: info

    call scalafx_ppotrf(aa, desc, uplo=uplo, info=info)
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(status, -1, "('scalafx_ppotrf failed in phermatinv (info: ',I0,')')",&
          & info)
    end if
    call scalafx_ppotri(aa, desc, uplo=uplo, info=info)
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(status, -1, "('scalafx_ppotri failed in phermatinv (info: ',I0,')')",&
          & info)
    end if

  end subroutine phermatinv


  !> Collect distributed BLACS array into duplicated arrays on each processor
  subroutine distrib2replicated(grid, desc, locArrayPart, glbDuplicatedArray)

    !> Group grid for the matrix
    type(blacsgrid) :: grid

    !> Dense matrix descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Local part of BLACS distributed array
    real(dp), intent(in) :: locArrayPart(:,:)

    !> Globally duplicated entire array
    real(dp), intent(out) :: glbDuplicatedArray(:,:)

    integer :: iLoc, jLoc, iGlb, jGlb, nLocalRows, nLocalCols, ierr

    call scalafx_getlocalshape(grid, desc, nLocalRows, nLocalCols)
    glbDuplicatedArray(:,:) = 0.0_dp
    do iLoc = 1, nLocalRows
      iGlb = scalafx_indxl2g(iLoc, desc(MB_), grid%myrow, desc(RSRC_), grid%nrow)
      do jLoc = 1, nLocalCols
        jGlb = scalafx_indxl2g(jLoc, desc(NB_), grid%mycol, desc(CSRC_), grid%ncol)
        glbDuplicatedArray(iGlb,jGlb) = locArrayPart(iLoc,jLoc)
      end do
    end do
    call blacsfx_gsum(grid, glbDuplicatedArray, rdest=-1, cdest=-1)

  end subroutine distrib2replicated

#:else

  !> Stub routine
  subroutine psymmatinv()
  end subroutine psymmatinv

  !> Stub routine
  subroutine phermatinv()
  end subroutine phermatinv

  !> Stub routine
  subroutine distrib2replicated()
  end subroutine distrib2replicated

#:endif

end module dftbp_math_scalafxext
