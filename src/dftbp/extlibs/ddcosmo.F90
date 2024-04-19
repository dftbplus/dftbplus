!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Proxy module to ddCOSMO implementation
module dftbp_extlibs_ddcosmo
  use ddcosmo_core, only : TDomainDecomposition, TDomainDecompositionInput, &
      & TDomainDecomposition_init, hsnorm, calcv, intrhs, prtsph, adjrhs, wghpot, &
      & ddupdate, fdoka, fdokb, fdoga
  use ddcosmo_solver, only : jacobi_diis, lx, lstarx, ldm1x, hnorm
  implicit none

  public

end module dftbp_extlibs_ddcosmo
