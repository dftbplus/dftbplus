!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Exporting mpifx functionality if compiled with mpi support, otherwise empty.
module dftbp_extlibs_mpifx
#:if WITH_MPI
  use libmpifx_module
#:endif
  implicit none
  public

end module dftbp_extlibs_mpifx
