!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
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
  use dftbp_extlibs_scalapackfx, only : DLEN_, scalafx_ppotrf, scalafx_ppotri
#:endif
  implicit none

  private
  public :: psymmatinv, phermatinv

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

#:else

  subroutine psymmatinv()
  end subroutine psymmatinv

  subroutine phermatinv()
  end subroutine phermatinv

#:endif

end module dftbp_math_scalafxext
