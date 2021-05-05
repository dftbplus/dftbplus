!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Implements various wrapped data types for use in creating ragged multi-dimensional arrays.
module dftbp_wrappedintr
  use dftbp_accuracy
  implicit none
  private

#:set FLAVOURS = [('logical','Logical',''), ('integer', 'Int', ''), ('real', 'Real', '(dp)'),&
  & ('complex', 'Cmplx', '(dp)') ]

#:for _, SUFFIX, _ in FLAVOURS
#:for DIM in [('1'), ('2')]

  public :: Twrapped${SUFFIX}$${DIM}$

#:endfor
#:endfor

#:for TYPE, NAME, PREC in FLAVOURS
#:for DIM, ARRAY in [('1',':'), ('2', ':,:')]

  !> ${DIM}$ dimensional ${TYPE}$
  type :: Twrapped${NAME}$${DIM}$
    ${TYPE}$${PREC}$, allocatable :: data(${ARRAY}$)
  end type Twrapped${NAME}$${DIM}$

#:endfor
#:endfor

end module dftbp_wrappedintr
