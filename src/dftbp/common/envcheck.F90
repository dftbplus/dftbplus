!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Contains functionality to check environment settings.
module dftbp_common_envcheck

  use, intrinsic :: iso_c_binding, only : c_char, c_int
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_globalenv, only : stdOut
  use dftbp_io_message, only : warning
  implicit none

  private
  public :: checkStackSize


  interface
    ! Interface for checking stacksize through C-routine.
    subroutine get_stacksize(cStack, iErr) bind(C, name="get_stacksize_c")

      import :: c_char, c_int
      implicit none

      !> Current stacksize
      integer(kind=c_int), intent(out) :: cStack

      !> Error status
      integer(kind=c_int), intent(out) :: iErr

    end subroutine get_stacksize
  end interface


contains

  !> Checks stacksize settings for optimal user experience/performance.
  subroutine checkStackSize(env)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !! Current stacksize
    integer :: cStack

    !! Error status
    integer :: iErr

    call get_stacksize(cStack, iErr)

    if (iErr /= 0) then
      write(stdOut, "(A,':',T30,A,I0,A)") "Current stacksize", "N/A (error code: ", iErr, ")"
    else
      if (cStack == -1 .or. cStack == 0) then
        write(stdOut, "(A,':',T30,A)") "Current stacksize", "unlimited"
      else
        write(stdOut, "(A,':',T30,I0,A)") "Current stacksize",&
            & nint(real(cStack, dp) / 1024.0_dp**2), " [Mb] (recommended: unlimited)"
        call warning("Current stacksize not set to unlimited or hard limit, which might cause"&
            & // new_line("A") // "   random crashes (e.g. segmentation faults). It is advised to&
            & unlimit the" // new_line("A") // "   stacksize by issuing 'ulimit -s unlimited'&
            & (Linux) or setting it to the " // new_line("A") // "   hard limit by 'ulimit -s hard'&
            & (Mac) in advance.")
      end if
    end if

  end subroutine checkStackSize

end module dftbp_common_envcheck
