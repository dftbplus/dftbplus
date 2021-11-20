!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Global names for output files for the main program.
module dftbp_dftbplus_outputfiles
  implicit none

  public

  !> Tagged output files (machine readable)
  character(*), parameter :: autotestTag = "autotest.tag"

  !> Detailed user output
  character(*), parameter :: userOut = "detailed.out"

  !> File for band structure and filling information
  character(*), parameter :: bandOut = "band.out"

  !> File for derivatives of band structure
  character(*), parameter :: derivEBandOut = "dE_band.out"

  !> File for derivatives of band structure
  character(*), parameter :: derivVBandOut = "dV_band.out"

  !> File accumulating data during an MD run
  character(*), parameter :: mdOut = "md.out"

  !> Machine readable tagged output
  character(*), parameter :: resultsTag = "results.tag"

  !> Second derivative of the energy with respect to atomic positions
  character(*), parameter :: hessianOut = "hessian.out"

  !> file name prefix for charge data
  character(*), parameter :: fCharges = "charges"

  !> file to stop code during geometry driver
  character(*), parameter :: fStopDriver = "stop_driver"

  !> file to stop code during scc cycle
  character(*), parameter :: fStopSCC = "stop_scc"

  !> file name for shift data
  character(*), parameter :: fShifts = "shifts.dat"

end module dftbp_dftbplus_outputfiles
