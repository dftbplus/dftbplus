!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines to add addition to repulsive pair contributions involving halogens
module halogenX
  use assert
  use accuracy, only : dp, mc
  use vdwdata
  use constants, only : AA__Bohr
  use message
  implicit none
  private

  public :: THalogenX, THalogenX_init

  !> Type for repulsive pairwise additions
  type :: THalogenX

    private

  contains

    procedure :: getCutOff
    !procedure :: getEnergy
    !procedure :: getForces
    !procedure :: getStress

  end type THalogenX

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

  !> Initialise structure
  subroutine THalogenX_init(this, species, speciesNames)

    type(THalogenX), intent(out) :: this

    !> Species for each atom
    integer, intent(in) :: species(:)

    !> Names of the atom types
    character(*), intent(in) :: speciesNames(:)

    logical :: tHalogen
    integer :: iSp1, iSp2
    character(mc) :: spName1, spName2

    tHalogen = .false.
    do iSp1 = 1, maxval(species)
      spName1 = speciesNames(iSp1)
      if (.not. any(spName1 == ["N","O"])) then
        cycle
      end if
      do iSp2 = 1, maxval(species)
        spName2 = speciesNames(iSp2)
        if (.not. any(spName2 == ["Cl","Br"]) .and. spName2 /= "I") then
          cycle
        end if
        write(*,*)iSp1,iSp2,trim(spName1),trim(spName2)
        tHalogen = .true.
      end do
    end do

    if (.not. tHalogen) then
      call error("No suitable O-X or N-X halogen combinations")
    end if

  end subroutine THalogenX_init


  !> Returns the distance over which the halogen correction decays
  function getCutOff(this, speciesNames)

    !> instance of the correction
    class(THalogenX), intent(in) :: this

    !> Names of the atom types
    character(*), intent(in) :: speciesNames(:)

    real(dp) :: getCutOff

    real(dp) :: dMax

    ! energy in kcal/mol for a pair at the cutoff distance
    real(dp), parameter :: minInteraction = 1.0E-14_dp

    if ( .not. (any(speciesNames == ["Cl","Br"]) .or. any(speciesNames == "I")) .and.&
        &.not. any(speciesNames == ["N","O"]) ) then
      call error("No relevant combinations of atomic elements present for halogen correction")
    end if

    ! add on the extra distance over which the interaction decays
    dMax = maxval(dab) + (-log(2.0_dp * minInteraction / c(1)) / c(2))**(1.0_dp/c(3))

    getCutOff = dMax * AA__Bohr

  end function getCutOff


  !> DFTB3-X term (Eqn. 5 of 10.1021/ct5009137)
  pure function fx(R, c, dab)

    real(dp) :: fx

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
