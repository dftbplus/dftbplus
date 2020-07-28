!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines to add addition to repulsive pair contributions involving halogens
!> from doi: 10.1021/ct5009137
module dftbp_halogenx
  use dftbp_assert
  use dftbp_accuracy, only : dp, mc
  use dftbp_vdwdata
  use dftbp_constants, only : AA__Bohr, Bohr__AA, kcal_mol__Hartree
  use dftbp_periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use dftbp_message
  implicit none
  private

  public :: THalogenX, THalogenX_init
  public :: halogenXSpecies1, halogenXSpecies2


  !> Type for repulsive pairwise additions
  type :: THalogenX

    private

    integer, allocatable :: interactionType(:,:)
    real(dp), allocatable :: radii(:)
    real(dp) :: maxDab = 0.0_dp
    integer :: nAtom
    real(dp) :: cutOff = 0.0_dp

  contains

    procedure :: getRCutOff
    procedure :: getEnergies
    procedure :: addGradients
    procedure :: getStress

  end type THalogenX

  !> Possible types for the first species in the interaction
  character(*), parameter :: halogenXSpecies1(2) = [character(1) ::&
      & 'O', 'N']

  !> Possible types for the second species in the interaction
  character(*), parameter :: halogenXSpecies2(3) = [character(2) ::&
      & 'Cl', 'Br', 'I']

  !> energy in kcal/mol for pair truncation
  real(dp), parameter :: minInteraction = 1.0E-14_dp

  !> Switching radii as multipliers of vdw radii
  real(dp), parameter :: R0 = 0.7_dp
  real(dp), parameter :: R1 = 0.8_dp

  !> Parameters from table 3 of doi: 10.1021/ct5009137 in AA
  !> O-Cl, O-Br, O-I, N-Cl, N-Br, N-I
  real(dp), parameter :: dab(size(halogenXSpecies1) * size(halogenXSpecies2)) =&
      & [1.237_dp, 1.099_dp, 1.313_dp, 1.526_dp, 1.349_dp, 1.521_dp]

  !> c constants from table 3 of doi: 10.1021/ct5009137 in kcal/mol, some weird power of inverse
  !> distance and dimensionless respectively
  real(dp), parameter :: c(3) = [7.761_dp, 0.050_dp, 4.518_dp]

contains


  !> Initialise structure
  subroutine THalogenX_init(this, species0, speciesNames)

    !> Instance
    type(THalogenX), intent(out) :: this

    !> Species for each atom in central cell
    integer, intent(in) :: species0(:)

    !> Names of the atom types
    character(*), intent(in) :: speciesNames(:)

    integer :: iSp1, iSp2, nSpecies
    character(mc) :: spName1, spName2

    nSpecies = maxval(species0)
    this%nAtom = size(species0)
    allocate(this%interactionType(nSpecies, nSpecies))
    allocate(this%radii(nSpecies))

    this%radii(:) = 0.0_dp
    do iSp1 = 1, nSpecies
      spName1 = speciesNames(iSp1)
      if (any(spName1 == halogenXSpecies1) .or. any(spName1 == halogenXSpecies2)) then
        call getVdwData(spName1, this%radii(iSp1))
      end if
    end do

    this%maxDab = 0.0_dp
    this%interactionType(:,:) = 0
    do iSp1 = 1, nSpecies
      spName1 = speciesNames(iSp1)
      if (.not. any(spName1 == halogenXSpecies1)) then
        cycle
      end if
      do iSp2 = 1, nSpecies
        spName2 = speciesNames(iSp2)
        if (.not. any(spName2 == halogenXSpecies2)) then
          cycle
        end if
        select case(trim(spName1) // trim(spName2))
        case('OCl')
          this%interactionType(iSp1,iSp2) = 1
        case('OBr')
          this%interactionType(iSp1,iSp2) = 2
        case('OI')
          this%interactionType(iSp1,iSp2) = 3
        case('NCl')
          this%interactionType(iSp1,iSp2) = 4
        case('NBr')
          this%interactionType(iSp1,iSp2) = 5
        case('NI')
          this%interactionType(iSp1,iSp2) = 6
        end select
        this%interactionType(iSp2, iSp1) = this%interactionType(iSp1, iSp2)
        this%maxDab = max(this%maxDab, dab(this%interactionType(iSp1, iSp2)))
      end do
    end do

    if (all(this%interactionType == 0)) then
      call error("No suitable (O-X or N-X) halogen combinations present for this correction")
    end if

    ! Distance over which the interaction decays to minInteraction
    this%cutoff = this%maxDab + (-log(2.0_dp * minInteraction / c(1)) / c(2))**(1.0_dp/c(3))
    this%cutoff = this%cutoff * AA__Bohr ! in a.u.

  end subroutine THalogenX_init


  !> Returns the distance over which the halogen correction decays
  function getRCutOff(this)

    !> instance of the correction
    class(THalogenX), intent(inout) :: this

    !> Returned distance
    real(dp) :: getRCutOff

    getRCutOff = this%cutoff

  end function getRCutOff


  !> Get energy contributions from halogen-X correction
  subroutine getEnergies(this, atomE, coords, species, neigh, img2CentCell)

    !> instance of the correction
    class(THalogenX), intent(in) :: this

    !> Resulting  energy contributions
    real(dp), intent(out) :: atomE(:)

    !> Current coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neigh

    !> Updated mapping to central cell.
    integer, intent(in) :: img2CentCell(:)

    ! Local variables
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
        if (this%interactionType(iSp1,iSp2) /= 0) then

          ! values in AA
          r = sqrt(neigh%neighDist2(iNeigh, iAt1)) * Bohr__AA
          rvdw = (this%radii(iSp1) + this%radii(iSp2)) * Bohr__AA

          eTmp = halogenSigma(r,rvdw) * fx(r, dab(this%interactionType(iSp1,iSp2)))
          eTmp = eTmp&
              & + (1.0_dp-halogenSigma(r,rvdw)) * fx(R0*rvdw, dab(this%interactionType(iSp1,iSp2)))

          ! convert back to a.u.
          eTmp = eTmp * kcal_mol__Hartree

          atomE(iAt1) = atomE(iAt1) + eTmp
          atomE(iAt2f) = atomE(iAt2f) + eTmp

        end if
      end do
    end do

  end subroutine getEnergies


  !> Gradient contribution from the halogen-X term
  subroutine addGradients(this, derivs, coords, species, neigh, img2CentCell)

    !> instance of the correction
    class(THalogenX), intent(in) :: this

    !> Derivatives to add contribution to to
    real(dp), intent(inout) :: derivs(:,:)

    !> Current coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Chemical species of atoms
    integer, intent(in) :: species(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neigh

    !> Updated mapping to central cell.
    integer, intent(in) :: img2CentCell(:)

    ! Local variables
    integer, allocatable :: nNeigh(:)
    integer :: iAt1, iNeigh, iAt2, iAt2f, iSp1, iSp2
    real(dp) :: rvdw, fTmp(3)

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
        if (this%interactionType(iSp1,iSp2) /= 0) then

          rvdw = (this%radii(iSp1) + this%radii(iSp2))

          call getEnergyDeriv(fTmp, coords(:,iAt2) - coords(:,iAt1), rvdw,&
              & this%interactionType(iSp1,iSp2))

          derivs(:,iAt1) = derivs(:,iAt1) + fTmp
          derivs(:,iAt2f) = derivs(:,iAt2f) - fTmp

        end if
      end do
    end do

  end subroutine addGradients


  !> Pairwise derivative terms
  subroutine getEnergyDeriv(intermed, vec, rvdw, iDab)

    !> Pair force component
    real(dp), intent(out) :: intermed(3)

    !> vector between atoms
    real(dp), intent(in) :: vec(3)

    !> Van der Waals radii
    real(dp), intent(in) :: rvdw

    !> Chemical combination for parameters
    integer, intent(in) :: iDab

    real(dp) :: rvdwAA, r, rTmp, fTmp

    rvdwAA = rvdw * Bohr__AA
    rTmp = sqrt(sum(vec**2))
    r = rTmp * Bohr__AA

    fTmp = halogendSigma(r,rvdwAA) * fx(r, dab(iDab)) + halogenSigma(r,rvdwAA) * dfx(r, dab(iDab)) &
        & -halogendSigma(r,rvdwAA) * fx(R0*rvdwAA, dab(iDab))

    intermed(:) = 2.0_dp * kcal_mol__Hartree * Bohr__AA * fTmp * vec / rTmp

  end subroutine getEnergyDeriv


  !> The stress tensor contribution from the halogen-X term
  subroutine getStress(this, st, coords, neigh, species, img2CentCell, cellVol)

    !> instance of the correction
    class(THalogenX), intent(in) :: this

    !> stress tensor
    real(dp), intent(out) :: st(:,:)

    !> coordinates (x,y,z, all atoms including possible images)
    real(dp), intent(in) :: coords(:,:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neigh

    !> Species of atoms in the central cell.
    integer, intent(in) :: species(:)

    !> indexing array for periodic image atoms
    integer, intent(in) :: img2CentCell(:)

    !> cell volume.
    real(dp), intent(in) :: cellVol

    integer :: iAt1, iNeigh, iAt2, iAt2f, ii, iSp1, iSp2
    real(dp) :: vect(3), intermed(3), prefac, rvdw
    integer, allocatable :: nNeigh(:)

    st(:,:) = 0.0_dp

    allocate(nNeigh(this%nAtom))
    call getNrOfNeighboursForAll(nNeigh, neigh, this%cutoff)

    do iAt1 = 1, this%nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeigh(iAt1)
        iAt2 = neigh%iNeighbour(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        if (this%interactionType(iSp1,iSp2) /= 0) then
          vect(:) = coords(:,iAt2) - coords(:,iAt1)
          rvdw = (this%radii(iSp1) + this%radii(iSp2))
          call getEnergyDeriv(intermed, vect, rvdw, this%interactionType(iSp1,iSp2))

          if (iAt1 == iAt2f) then
            prefac = 0.5_dp
          else
            prefac = 1.0_dp
          end if
          do ii = 1, 3
            st(:, ii) = st(:, ii) + prefac * intermed * vect(ii)
          end do
        end if
      end do
    end do

    st = st / cellVol

  end subroutine getStress


  !> DFTB3-X term (Eqn. 5 of doi: 10.1021/ct5009137)
  pure function fx(R, dab)

    !> result in kcal/mol
    real(dp) :: fx

    !> distance in AA
    real(dp), intent(in) :: R

    !> distance cut-off in AA
    real(dp), intent(in) :: dab

    fx = 0.5_dp * c(1) * exp(-c(2) * ((R - dab)**c(3)))

  end function fx


  !> Derivative of DFTB3-X term wrt. R
  pure function dfx(R, dab)

    !> result in kcal/mol AA
    real(dp) :: dfx

    !> distance in AA
    real(dp), intent(in) :: R

    !> distance cut-off in AA
    real(dp), intent(in) :: dab

    dfx = fx(R, dab)
    dfx = dfx * c(2) * c(3) * (R - dab)**(c(3)-1.0_dp)

  end function dfx


  !> Switch function (Eqn. 8 of doi: 10.1021/ct5009137)
  pure function halogenSigma(R, rvdw) result(out)

    !> Distance in a.u.
    real(dp), intent(in) :: R

    !> Van der Waals radii sum in a.u.
    real(dp), intent(in) :: rvdw

    !> result fraction
    real(dp) :: out

    real(dp) :: x

    if (R < R0*rvdw) then
      out = 0.0_dp
    elseif (R > R1*rvdw) then
      out = 1.0_dp
    else
      x = (R - R0*rvdw) / ((R1 - R0)*rvdw)
      out = -20.0_dp*x**7 +70.0_dp*x**6 -84.0_dp*x**5 +35.0_dp*x**4
    end if

  end function halogenSigma


  !> Derivative of switch function (Eqn. 12 of doi: 10.1021/ct5009137)
  pure function halogendSigma(R, rvdw) result(out)

    !> Distance in a.u.
    real(dp), intent(in) :: R

    !> Van der Waals radii sum in a.u.
    real(dp), intent(in) :: rvdw

    !> result
    real(dp) :: out

    real(dp) :: x

    if (R > R0*rvdw .and. R < R1*rvdw) then
      x = (R - R0*rvdw) / ((R1 - R0)*rvdw)
      out = -140.0_dp*x**6 +420.0_dp*x**5 -420.0_dp*x**4 +140.0_dp*x**3
      out = out / ((R1 - R0)*rvdw)
    else
      out = 0.0_dp
    end if

  end function halogendSigma

end module dftbp_halogenx
