!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Labels for various numerically enumerated solvers
module solvertypes
  implicit none

  ! Charge mixers

  !> Simple charge mixer
  integer, parameter :: mixerSimple = 1
  !> Anderson charge mixer
  integer, parameter :: mixerAnderson = 2
  !> Broyden charge mixer
  integer, parameter :: mixerBroyden = 3
  !> DIIS charge mixer
  integer, parameter :: mixerDIIS = 4


  ! geometry optimisers

  !> Null case
  integer, parameter :: optNull = 0
  !> Steepest descent optimiser
  integer, parameter :: optSD = 1
  !> Conjugate gradient optimiser
  integer, parameter :: optCG = 2
  !> gDIIS optimiser
  integer, parameter :: optDIIS = 3
  !> LBFGS optimiser
  integer, parameter :: optLBFGS = 4


  ! force evaluation method

  !> Conventional forces
  integer, parameter :: forceOrig = 0
  !> convergence corrected at 0 temperature
  integer, parameter :: forceDynT0 = 1
  !> convergence corrected at finite temperature
  integer, parameter :: forceDynT = 2


  ! eigen- or alternative solvers

  !> QR solver
  integer, parameter :: solverQR = 1
  !> Divide and conquer eigensolver
  integer, parameter :: solverDAC = 2
  !> Relatively robust representation eigensolver
  integer, parameter :: solverRR = 3
  !> Green's function solver
  integer, parameter :: solverGF = 4
  !> Transport only transmission solution
  integer, parameter :: solverOnlyTransport = 5


  ! electrostatic solution method

  !> Softened coulombic with gamma
  integer, parameter :: gammaf = 0
  !> Poisson equation solver
  integer, parameter :: poisson = 1

end module solvertypes
