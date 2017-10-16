!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Interface with plumed routines.
module plumed

  public

  interface


    !> creats a plumed object for use
    subroutine plumed_f_gcreate() bind(c)

    end subroutine plumed_f_gcreate


    !> requests plumed perform some action
    subroutine plumed_f_gcmd(key, val) bind(c)
      
      !> key string determines what plumed will be doing
      character(len=*), intent(in) :: key

      !> the actual value being passed to plumed to be read or acted upon
      unspecified_type, intent(inout) :: val(*)

    end subroutine plumed_f_gcmd


    !> finalises plumed
    subroutine plumed_f_gfinalize() bind(c)

    end subroutine plumed_f_gfinalize 


  end interface

end module plumed
