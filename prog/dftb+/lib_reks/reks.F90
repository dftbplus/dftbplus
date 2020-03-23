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
  use dftbp_reksproperty
  use dftbp_reksvar
!  import dftbp_rekscommon
!  import dftbp_reksen
!  import dftbp_reksfon
!  import dftbp_reksinterface
!  import dftbp_reksproperty
!  import dftbp_reksvar

  implicit none

  private

  !> dftbp_rekscommon modules used in main.F90
  public :: checkGammaPoint
  public :: qm2udL, ud2qmL
  public :: qmExpandL!, udExpandL

  !> dftbp_reksen modules used in main.F90
  public :: constructMicrostates, calcWeights
  public :: activeOrbSwap, getFilling, calcSaReksEnergy
  public :: getFockandDiag, guessNewEigvecs
  public :: adjustEigenval, solveSecularEqn

  !> dftbp_reksfon module used in main.F90
  public :: optimizeFONs

  !> dftbp_reksinterface modules used in main.F90
  public :: getReksGradients, getReksGradProperties
  public :: getReksStress

  !> dftbp_rekspreoperty modules used in main.F90
  public :: getUnrelaxedDMandTDP
  public :: getDipoleIntegral, getDipoleMomentMatrix, getReksOsc

  !> dftbp_reksvar module used in main.F90
  public :: TReksCalc

end module dftbp_reks
