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
  use periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use message
  implicit none
  private

  public :: THalogenX, THalogenX_init

  !> Type for repulsive pairwise additions
  type :: THalogenX

    private

    integer, allocatable :: relevantSpecies(:,:)
    real(dp) :: maxDab = 0.0_dp
    integer :: nAtom
    real(dp) :: cutOff = 0.0_dp

  contains

    procedure :: getRCutOff
    procedure :: getEnergies
    !procedure :: addGradients
    !procedure :: getStress

  end type THalogenX

  ! energy in kcal/mol for pair truncation
  real(dp), parameter :: minInteraction = 1.0E-14_dp

  !> Switching radii as multipliers of vdw radii
  real(dp), parameter :: R0 = 0.7_dp
  real(dp), parameter :: R1 = 0.8_dp

  !> Parameters from table 3 of 10.1021/ct5009137 in AA
  real(dp), parameter :: dab(6) = [1.237_dp, 1.099_dp, 1.313_dp, 1.526_dp, 1.349_dp, 1.521_dp]

  !> c constants from table 3 of 10.1021/ct5009137 in kcal/mol, some weird power of inverse distance
  !> and dimensionless respectively
  real(dp), parameter :: c(3) = [7.761_dp, 0.050_dp, 4.518_dp]

contains

  !> Initialise structure
  subroutine THalogenX_init(this, species0, speciesNames)

    type(THalogenX), intent(out) :: this

    !> Species for each atom in central cell
    integer, intent(in) :: species0(:)

    !> Names of the atom types
    character(*), intent(in) :: speciesNames(:)

    logical :: tHalogen
    integer :: iSp1, iSp2, ii, jj, nSpecies
    character(mc) :: spName1, spName2

    nSpecies = maxval(species0)
    this%nAtom = size(species0)
    tHalogen = .false.
    allocate(this%relevantSpecies(nSpecies, nSpecies))
    this%relevantSpecies(:,:) = 0
    this%maxDab = 0.0_dp
    do iSp1 = 1, nSpecies
      spName1 = speciesNames(iSp1)
      if (.not. any(spName1 == ["N","O"])) then
        cycle
      end if
      do iSp2 = 1, nSpecies
        spName2 = speciesNames(iSp2)
        if (.not. any(spName2 == ["Cl","Br"]) .and. spName2 /= "I") then
          cycle
        end if

        select case(trim(spName1)//trim(spName2))
        case('OCl')
          this%relevantSpecies(iSp1,iSp2) = 1
        case('OBr')
          this%relevantSpecies(iSp1,iSp2) = 2
        case('OI')
          this%relevantSpecies(iSp1,iSp2) = 3
        case('NCl')
          this%relevantSpecies(iSp1,iSp2) = 4
        case('NBr')
          this%relevantSpecies(iSp1,iSp2) = 5
        case('NI')
          this%relevantSpecies(iSp1,iSp2) = 6
        end select
        this%maxDab = max(this%maxDab, dab(this%relevantSpecies(iSp1,iSp2)))
        tHalogen = .true.
      end do
    end do

    if (.not. tHalogen) then
      call error("No suitable O-X or N-X halogen combinations for this correction")
    end if

  end subroutine THalogenX_init


  !> Returns the distance over which the halogen correction decays
  function getRCutOff(this)

    !> instance of the correction
    class(THalogenX), intent(inout) :: this

    real(dp) :: getRCutOff

    real(dp) :: cutoff

    ! Distance over which the interaction decays to minInteraction
    this%cutoff = this%maxDab + (-log(2.0_dp * minInteraction / c(1)) / c(2))**(1.0_dp/c(3))
    this%cutoff = this%cutoff * AA__Bohr
    getRCutOff = this%cutoff

  end function getRCutOff


  subroutine getEnergies(this, coords, neigh, img2CentCell)

    !> instance of the correction
    class(THalogenX), intent(in) :: this

    !> Current coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neigh

    !> Updated mapping to central cell.
    integer, intent(in) :: img2CentCell(:)

    integer, allocatable :: nNeigh(:)
    integer :: iAt1, iNeig, iAt2, iAt2f

    allocate(nNeigh(this%nAtom))
    call getNrOfNeighboursForAll(nNeigh, neigh, this%cutoff)



  end subroutine getEnergies

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
