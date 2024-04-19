!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> REKS and SI-SA-REKS formulation in DFTB as developed by Lee et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet apart from Gamma point.
!> * Orbital potentials or spin-orbit or external E-field does not work yet.
!> * Only for closed shell system.
!> * Onsite corrections are not included in this version
module dftbp_reks_reks
  use dftbp_reks_rekscommon, only : checkGammaPoint, qm2udL, qmExpandL, ud2qmL
  use dftbp_reks_reksen, only : activeOrbSwap, calcSaReksEnergy, calcWeights, constructMicrostates,&
      & getFilling, getFockandDiag, guessNewEigvecs, setReksTargetEnergy
  use dftbp_reks_reksfon, only : optimizeFons
  use dftbp_reks_reksinterface, only : getReksEnProperties, getReksGradProperties,&
      & getReksGradients, getReksStress, getStateInteraction
  use dftbp_reks_reksio, only : printReksMicrostates, printReksSAInfo, printSaReksEnergy
  use dftbp_reks_reksvar, only : REKS_init, reksTypes, TReksCalc, TReksInp

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

  !> dftbp_reks_rekscommon modules used in main.F90
  public :: checkGammaPoint
  public :: qm2udL, ud2qmL
  public :: qmExpandL!, udExpandL

  !> dftbp_reks_reksen modules used in main.F90, mainio.F90
  public :: constructMicrostates, calcWeights
  public :: activeOrbSwap, getFilling, calcSaReksEnergy
  public :: getFockandDiag, guessNewEigvecs
  public :: setReksTargetEnergy

  !> dftbp_reks_reksfon module used in main.F90
  public :: optimizeFons

  !> dftbp_reks_reksinterface modules used in main.F90
  public :: getStateInteraction, getReksEnProperties
  public :: getReksGradients, getReksGradProperties
  public :: getReksStress

  !> dftbp_reks_reksio modules used in main.F90
  public :: printReksMicrostates, printSaReksEnergy, printReksSAInfo

  !> dftbp_reks_reksvar module used in main.F90, mainio.F90, inputdata.F90, initprogram.F90, parser.F90
  public :: TReksInp, TReksCalc, REKS_init, reksTypes

end module dftbp_reks_reks
