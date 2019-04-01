!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines to add additiona to repulsive pair contributions
module halogenX
  use assert
  use accuracy, only : dp
  use vdwdata
  use message, only : warning
  implicit none
  private

  public :: TRepAdditionsInp
  public :: TRepAdditions, TRepAdditions_init


  !> Type for input variables to set repulsive pairwise additions to the DFTB model
  type :: TRepAdditionsInp

    !> Correction from http://dx.doi.org/10.1021/ct5009137
    logical :: tHalogen

  end type TRepAdditionsInp


  !> Type for repulsive pairwise additions
  type :: TRepAdditions

    private

    !> Correction from http://dx.doi.org/10.1021/ct5009137
    logical :: tHalogen
    real(dp) :: c(3)
    real(dp), allocatable :: dij(:,:)

  contains

    procedure :: getCutOff
    procedure :: getEnergy
    procedure :: getForces
    procedure :: getStress

  end type TRepAdditions

  !> Switching radius as multipliers of vdw radii
  real(dp), parameter :: R0 = 0.7_dp
  real(dp), parameter :: R1 = 0.8_dp

  !> Parameters from table 3 of 10.1021/ct5009137 in AA
  real(dp), parameter :: dab(6) = [1.237_dp, 1.099_dp, 1.313_dp, 1.526_dp, 1.349_dp, 1.521_dp]

  !> Atomic numbers of pairs of corrected atoms
  integer, parameter :: pairs(6,2) = reshape([16,16,16,15,15,15,17,35,53,17,35,53],[6,2],&
      & order=[1,2])

  !> C constants from table 3 of 10.1021/ct5009137 in kcal/mol, some weird power of inverse distance
  !> and dimensionless respectively
  real(dp), parameter :: c(3) = [7.761_dp, 0.050_dp, 4.518_dp]

contains


  subroutine TRepAdditions_init(this, inp, species, speciesNames)

    !>
    integer, intent(in) :: species(:)

    !> Names of the atom types
    character(*), intent(in) :: speciesNames(:)

  end subroutine TRepAdditions_init


  !> Returns the distance over which the halogen correction decays
  pure function getCutOff(speciesNames)

    !> Names of the atom types
    character(*), intent(in) :: speciesNames(:)

    integer :: nSpecies, iSp1m iSp2
    real(dp) :: dMax, radius
    character(mc) :: spName1, spName2

    real(dp), parameter :: minInteraction = 1.0E-20_dp

    nSpecies = size(speciesNames)

    do iSp1 = 1, size(speciesNames)
      spName1 = speciesNames(iSp1)
      call getVdwData(speciesNames(iSp1), radius, found=tFoundRadius)
      do iSp2 = 1, size(speciesNames)

        spName2 = speciesNames(iSp2)


        
      end do
    end do

    ! add on the extra distance over which
    dmax = dmax + (-log(2.0_dp * minInteraction / c1) / c2)**(1.0_dp/c3)

  end function getCutOff


  !> DFTB3-X term (Eqn. 5 of 10.1021/ct5009137)
  pure function fx(R, c, dab)

    real(dp), intent(in) :: R
    real(dp), intent(in) :: c(3)
    real(dp), intent(in) :: dab

    fx = 0.5_dp * c(1) * exp(-c(2) * (R - dab)**c(3))

  end function fx


  !> Switch function (Eqn. 8 of 10.1021/ct5009137)
  pure function halogenSigma(R) result(out)
    real(dp), intent(in) :: R
    real(dp) :: out

    real(dp) :: x

    if (R < R0) then
      out = 1.0_dp
    elseif (R > R1) then
      out = 0.0_dp
    else
      x = (R - R0) / (R1 - R0)
      out = -20.0_dp*x**7 +70.0_dp*x**6 -84.0_dp*x**5 +35.0_dp*x**4
    end if

  end function halogenSigma


  !> Derivative of switch function (Eqn. 12 of 10.1021/ct5009137)
  pure function halogendSigma(R) result(out)
    real(dp), intent(in) :: R
    real(dp) :: out

    real(dp) :: x

    if (R > R0 .and. R < R1) then
      x = (R - R0) / (R1 - R0)
      out = -140.0_dp*x**6 +420.0_dp*x**5 -420.0_dp*x**4 +140.0_dp*x**3
    end if

  end function halogendSigma


end module halogenX
