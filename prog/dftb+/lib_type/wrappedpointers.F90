!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Contains simple derived types wrapping pointers to scalars and arrays.
module wrappedpointers
  use accuracy
  implicit none
  private

  public :: PIntR0, PIntR1, PIntR2
  public :: PRealR0, PRealR1, PRealR2

  !!* Wrapper for rank 0 integer.
  type PIntR0
    integer, pointer :: pt
  end type PIntR0

  !!* Wrapper for rank 1 integer.
  type PIntR1
    integer, pointer :: pt(:)
  end type PIntR1
  
  !!* Wrapper for rank 2 integer.
  type PIntR2
    integer, pointer :: pt(:,:)
  end type PIntR2

  !!* Wrapper for rank 0 real.
  type PRealR0
    real(dp), pointer :: pt
  end type PRealR0
  
  !!* Wrapper for rank 1 real.
  type PRealR1
    real(dp), pointer :: pt(:)
  end type PRealR1

  !!* Wrapper for rank 2 real.
  type PRealR2
    real(dp), pointer :: pt(:,:)
  end type PRealR2

end module wrappedpointers
