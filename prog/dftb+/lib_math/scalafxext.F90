!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

!> Common mathematical operations built out of multiple scalapack calls
module dftbp_math_scalafxext
  use dftbp_common_accuracy, only : lc, dp
  use dftbp_io_message, only : error
#:if WITH_SCALAPACK
  use dftbp_extlibs_scalapackfx, only : DLEN_, scalafx_ppotrf, scalafx_ppotri
#:endif
  implicit none
  
  private
  public :: psymmatinv, phermatinv

contains

#:if WITH_SCALAPACK

  !> Inversion of a symmetric matrix
  subroutine psymmatinv(desc, aa, uplo, info)

    !> Matrix descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Matrix to invert on entry, inverted matrix on exit
    real(dp), intent(inout) :: aa(:,:)

    !> Whether upper or lower triangle is specified in the matrix ("U" or "L", default: "L")
    character, intent(in), optional :: uplo

    !> Info flag: if non-zero, inversion failed
    integer, intent(out), optional :: info

    integer :: info0
    character(lc) :: buffer

    call scalafx_ppotrf(aa, desc, uplo=uplo, info=info0)
    if (info0 /= 0) then
      if (present(info)) then
        info = info0
        return
      else
        write(buffer, "(A,I0,A)") "scalafx_ppotrf failed in psymmatinv (info: ", info0, ")"
        call error(buffer)
      end if
    end if

    call scalafx_ppotri(aa, desc, uplo=uplo, info=info)
    if (info0 /= 0) then
      if (present(info)) then
        info = info0
        return
      else
        write(buffer, "(A,I0,A)") "scalafx_ppotri failed in psymmatinv (info: ", info0, ")"
        call error(buffer)
      end if
    end if

  end subroutine psymmatinv

#:else

  subroutine psymmatinv()
  end subroutine psymmatinv

#:endif


#:if WITH_SCALAPACK

  !> Inversion of a hermitian matrix
  subroutine phermatinv(desc, aa, uplo, info)

    !> Matrix descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Matrix to invert on entry, inverted matrix on exit
    complex(dp), intent(inout) :: aa(:,:)

    !> Whether upper or lower triangle is specified in the matrix ("U" or "L", default: "L")
    character, intent(in), optional :: uplo

    !> Info flag: if non-zero, inversion failed
    integer, intent(out), optional :: info

    integer :: info0
    character(lc) :: buffer

    call scalafx_ppotrf(aa, desc, uplo=uplo, info=info0)
    if (info0 /= 0) then
      if (present(info)) then
        info = info0
        return
      else
        write(buffer, "(A,I0,A)") "scalafx_ppotrf failed in psymmatinv (info: ", info0, ")"
        call error(buffer)
      end if
    end if

    call scalafx_ppotri(aa, desc, uplo=uplo, info=info)
    if (info0 /= 0) then
      if (present(info)) then
        info = info0
        return
      else
        write(buffer, "(A,I0,A)") "scalafx_ppotri failed in psymmatinv (info: ", info0, ")"
        call error(buffer)
      end if
    end if

  end subroutine phermatinv

#:else

  subroutine phermatinv()
  end subroutine phermatinv

#:endif


end module dftbp_math_scalafxext
