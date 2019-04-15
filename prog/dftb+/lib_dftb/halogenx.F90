!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines to add addition to repulsive pair contributions involving halogens
module dftbp_halogenX
  use dftbp_assert
  use dftbp_accuracy, only : dp, mc
  use dftbp_vdwdata
  use dftbp_constants, only : AA__Bohr, Bohr__AA, kcal_mol__Hartree
  use dftbp_periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use dftbp_message
  implicit none
  private

  public :: THalogenX, THalogenX_init

  !> Type for repulsive pairwise additions
  type :: THalogenX

    private

    integer, allocatable :: relevantSpecies(:,:)
    real(dp), allocatable :: radii(:)
    real(dp) :: maxDab = 0.0_dp
    integer :: nAtom
    real(dp) :: cutOff = 0.0_dp

  contains

    procedure :: getRCutOff
    procedure :: getEnergies
    procedure :: addGradients
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
    allocate(this%relevantSpecies(nSpecies, nSpecies))
    this%relevantSpecies(:,:) = 0
    allocate(this%radii(nSpecies))
    this%radii(:) = 0.0_dp

    tHalogen = .false.
    this%maxDab = 0.0_dp
    do iSp1 = 1, nSpecies
      spName1 = speciesNames(iSp1)
      if (any(spName1 == ["N","O","I"]) .or. any(spName1 == ["Cl","Br"])) then
        call getVdwData(spName1, this%radii(iSp1))
      end if
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
        this%relevantSpecies(iSp2,iSp1) = this%relevantSpecies(iSp1,iSp2)
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


  subroutine getEnergies(this, atomE, coords, species, neigh, img2CentCell)

    !> instance of the correction
    class(THalogenX), intent(in) :: this

    real(dp), intent(out) :: atomE(:)

    !> Current coordinates
    real(dp), intent(in) :: coords(:,:)

    integer, intent(in) :: species(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neigh

    !> Updated mapping to central cell.
    integer, intent(in) :: img2CentCell(:)

    integer, allocatable :: nNeigh(:)
    integer :: iAt1, iNeigh, iAt2, iAt2f, iSp1, iSp2
    real(dp) :: r, rvdw, eTmp

    allocate(nNeigh(this%nAtom))
    call getNrOfNeighboursForAll(nNeigh, neigh, this%cutoff)

    atomE(:) = 0.0_dp

    do iAt1 = 1, this%nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeigh(iAt1)
        iAt2 = neigh%iNeighbour(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        if (this%relevantSpecies(iSp1,iSp2) /= 0) then
          r = sqrt(neigh%neighDist2(iNeigh, iAt1))
          rvdw = this%radii(iSp1) + this%radii(iSp2)
          eTmp = halogenSigma(r,rvdw) * fx(r, dab(this%relevantSpecies(iSp1,iSp2)))
          eTmp = eTmp&
              & + (1.0_dp-halogenSigma(r,rvdw)) * fx(R0*rvdw, dab(this%relevantSpecies(iSp1,iSp2)))
          eTmp = eTmp * 0.5_dp * kcal_mol__Hartree
          atomE(iAt1) = atomE(iAt1) + eTmp
          atomE(iAt2f) = atomE(iAt2f) + eTmp
        end if
      end do
    end do

  end subroutine getEnergies


  subroutine addGradients(this, derivs, coords, species, neigh, img2CentCell)

    !> instance of the correction
    class(THalogenX), intent(in) :: this

    !> Derivatives to add contribution to to
    real(dp), intent(inout) :: derivs(:,:)

    !> Current coordinates
    real(dp), intent(in) :: coords(:,:)

    integer, intent(in) :: species(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neigh

    !> Updated mapping to central cell.
    integer, intent(in) :: img2CentCell(:)

    integer, allocatable :: nNeigh(:)
    integer :: iAt1, iNeigh, iAt2, iAt2f, iSp1, iSp2
    real(dp) :: r, rvdw, eTmp, fTmp(3)

    allocate(nNeigh(this%nAtom))
    call getNrOfNeighboursForAll(nNeigh, neigh, this%cutoff)

    do iAt1 = 1, this%nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeigh(iAt1)
        iAt2 = neigh%iNeighbour(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        if (iAt2f == iAt1) then
          cycle
        end if
        iSp2 = species(iAt2f)
        if (this%relevantSpecies(iSp1,iSp2) /= 0) then
          r = sqrt(neigh%neighDist2(iNeigh, iAt1))
          rvdw = this%radii(iSp1) + this%radii(iSp2)

          eTmp = halogendSigma(r,rvdw) * fx(r, dab(this%relevantSpecies(iSp1,iSp2)))&
              & +halogenSigma(r,rvdw) * dfx(r, dab(this%relevantSpecies(iSp1,iSp2)))
          eTmp = eTmp&
              & +(1.0_dp-halogenSigma(r,rvdw)) * dfx(R0*rvdw, dab(this%relevantSpecies(iSp1,iSp2)))&
              & -halogendSigma(r,rvdw) * fx(R0*rvdw, dab(this%relevantSpecies(iSp1,iSp2)))
          eTmp = eTmp * 0.5_dp * kcal_mol__Hartree * Bohr__AA

          fTmp(:) = ( coords(:,iAt2f) - coords(:,iAt1) ) * eTmp / r

          derivs(:,iAt1) = derivs(:,iAt1) + fTmp
          derivs(:,iAt2f) = derivs(:,iAt2f) - fTmp

        end if
      end do
    end do

  end subroutine addGradients

  !> DFTB3-X term (Eqn. 5 of 10.1021/ct5009137)
  pure function fx(R, dab)

    !> result in kcal/mol
    real(dp) :: fx

    !> distance in a.u.
    real(dp), intent(in) :: R

    !> distance cut-off in AA
    real(dp), intent(in) :: dab

    fx = 0.5_dp * c(1) * exp(-c(2) * (R * Bohr__AA - dab)**c(3))

  end function fx


  !> Derivative of DFTB3-X term wrt. R
  pure function dfx(R, dab)

    !> result in kcal/mol AA
    real(dp) :: dfx

    !> distance in a.u.
    real(dp), intent(in) :: R

    !> distance cut-off in AA
    real(dp), intent(in) :: dab

    dfx = fx(R, dab)
    dfx = dfx * c(2) * c(3) * (R * Bohr__AA - dab)**(c(3)-1.0_dp)

  end function dfx


  !> Switch function (Eqn. 8 of 10.1021/ct5009137)
  pure function halogenSigma(R, rvdw) result(out)

    !> Distance in a.u.
    real(dp), intent(in) :: R

    !> Van der Waals radii sum in a.u.
    real(dp), intent(in) :: rvdw

    !> result fraction
    real(dp) :: out

    real(dp) :: x

    if (R < R0*rvdw) then
      out = 1.0_dp
    elseif (R > R1*rvdw) then
      out = 0.0_dp
    else
      x = (R - R0) / (R1 - R0)
      out = -20.0_dp*x**7 +70.0_dp*x**6 -84.0_dp*x**5 +35.0_dp*x**4
    end if

  end function halogenSigma


  !> Derivative of switch function (Eqn. 12 of 10.1021/ct5009137)
  pure function halogendSigma(R, rvdw) result(out)

    !> Distance in a.u.
    real(dp), intent(in) :: R

    !> Van der Waals radii sum in a.u.
    real(dp), intent(in) :: rvdw

    !> result
    real(dp) :: out

    real(dp) :: x

    if (R > R0*rvdw .and. R < R1*rvdw) then
      x = (R - R0) / (R1 - R0)
      out = -140.0_dp*x**6 +420.0_dp*x**5 -420.0_dp*x**4 +140.0_dp*x**3
    end if

  end function halogendSigma

end module dftbp_halogenX
