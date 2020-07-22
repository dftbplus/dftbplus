!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> REKS and SI-SA-REKS formulation in DFTB as developed by Lee et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet appart from Gamma point.
!> * Orbital potentials or spin-orbit or external E-field does not work yet.
!> * Only for closed shell system.
!> * Onsite corrections are not included in this version
module dftbp_reks

  use dftbp_rekscommon
  use dftbp_reksen
  use dftbp_reksfon
  use dftbp_reksinterface
  use dftbp_reksio
  use dftbp_reksvar

  implicit none

  private

  !> In REKS method, there is a symmetry for the microstates due to the restricted scheme.
  !> For the reduce of memory allocation, I make two representation used in REKS.
  !>
  !> Notation : 1u means up spin channel of 1st microstate
  !>
  !> 1. my_ud : It has only up spin channel for the microstates
  !>            1 = 1u, 2 = 2u, 3 = 3u, 4 = 4u, 5 = 5u, 6 = 6u.
  !> 2. my_qm : It has sum or difference between the microstates which have symmetry between them.
  !>            1 = 1u + 1d (= 1u + 1u), 2 = 2u + 2d (= 2u + 2u), 3 = 3u + 3d (= 3u + 4u),
  !>            4 = 3u - 3d (= 3u - 4u), 5 = 5u + 5d (= 5u + 6u), 6 = 5u - 5d (= 5u - 6u).

  !> dftbp_rekscommon modules used in main.F90
  public :: checkGammaPoint
  public :: qm2udL, ud2qmL
  public :: qmExpandL!, udExpandL

  !> dftbp_reksen modules used in main.F90, mainio.F90
  public :: constructMicrostates, calcWeights
  public :: activeOrbSwap, getFilling, calcSaReksEnergy
  public :: getFockandDiag, guessNewEigvecs
  public :: setReksTargetEnergy

  !> dftbp_reksfon module used in main.F90
  public :: optimizeFons

  !> dftbp_reksinterface modules used in main.F90
  public :: getStateInteraction, getReksEnProperties
  public :: getReksGradients, getReksGradProperties
  public :: getReksStress

  !> dftbp_reksio modules used in main.F90
  public :: printReksMicrostates, printSaReksEnergy, printReksSAInfo

  !> dftbp_reksvar module used in main.F90, mainio.F90, inputdata.F90, initprogram.F90, parser.F90
  public :: TReksInp, TReksCalc, REKS_init, reksTypes

end module dftbp_reks
