!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Proxy module to ddCOSMO implementation
module dftbp_extlibs_ddcosmo
  use ddcosmo_core, only : adjrhs, calcv, ddupdate, fdoga, fdoka, fdokb, hsnorm, intrhs, prtsph,&
      & TDomainDecomposition, TDomainDecomposition_init, TDomainDecompositionInput, wghpot
  use ddcosmo_solver, only : hnorm, jacobi_diis, ldm1x, lstarx, lx
  implicit none

  public

end module dftbp_extlibs_ddcosmo
