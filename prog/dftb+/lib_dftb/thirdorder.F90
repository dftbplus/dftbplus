!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines implementing the full 3rd order DFTB.
module dftbp_thirdorder
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_commontypes, only : TOrbitals
  use dftbp_shortgamma, only : expGammaCutoff
  use dftbp_periodic, only : TNeighbourList, getNrOfNeighbours
  use dftbp_charges
  implicit none
  private

  public :: TThirdOrderInp, TThirdOrder, ThirdOrder_init


  !> Input for the 3rd order module
  type TThirdOrderInp

    !> Orbital information
    type(TOrbitals), pointer :: orb

    !> Hubbard U values. Shape: [nShell, nSpecies]
    real(dp), allocatable :: hubbUs(:,:)

    !> Hubbard U derivatives. Shape: [nShell, nSpecies]
    real(dp), allocatable :: hubbUDerivs(:,:)

    !> Whether interaction should be damped when given atom is involved.
    logical, allocatable :: damped(:)

    !> Exponention of the damping
    real(dp) :: dampExp

    !> Whether third order should be considered shell resolved. If not, only the first shell of each
    !> atom in hubbUs and hubbUDerivs is used.
    logical :: shellResolved

  end type TThirdOrderInp


  !> Internal status of third order.
  type TThirdOrder
    integer :: nSpecies, nAtoms, mShells, mShellsReal
    logical :: shellResolved
    integer, allocatable :: nShells(:)
    real(dp), allocatable :: UU(:,:)
    real(dp), allocatable :: dUdQ(:,:)
    real(dp), allocatable :: shift1(:,:), shift2(:,:), shift3(:)
    real(dp), allocatable :: chargesPerAtom(:)
    real(dp), allocatable :: chargesPerShell(:,:)
    real(dp), allocatable :: cutoffs(:,:)
    integer, allocatable :: nNeigh(:,:), nNeighMax(:)
    real(dp), allocatable :: gamma3ab(:,:,:,:), gamma3ba(:,:,:,:)
    logical, allocatable :: damped(:)
    real(dp) :: dampExp
    real(dp) :: maxCutoff
  contains
    procedure :: getCutoff
    procedure :: updateCoords
    procedure :: updateCharges
    procedure :: getShifts
    procedure :: getEnergyPerAtom
    procedure :: getEnergyPerAtomXlbomd
    procedure :: addGradientDc
    procedure :: addGradientDcXlbomd
    procedure :: addStressDc
  end type TThirdOrder

contains


  !> Initializes instance.
  subroutine ThirdOrder_init(this, inp)

    !> Instance.
    type(TThirdOrder), intent(out) :: this

    !> Input data.
    type(TThirdOrderInp), intent(in) :: inp

    this%nAtoms = size(inp%orb%nOrbAtom)
    this%mShellsReal = inp%orb%mShell
    this%nSpecies = size(inp%hubbUs, dim=2)
    this%shellResolved = inp%shellResolved
    if (this%shellResolved) then
      this%mShells = this%mShellsReal
      this%nShells = inp%orb%nShell
      this%UU = inp%hubbUs
      this%dUdQ = inp%hubbUDerivs
    else
      this%mShells = 1
      allocate(this%nShells(this%nSpecies))
      this%nShells(:) = 1
      this%UU = inp%hubbUs(1:1, :)
      this%dUdQ = inp%hubbUDerivs(1:1, :)
    end if

    allocate(this%cutoffs(this%nSpecies, this%nSpecies))
    call calcCutoffs(this%UU, this%nShells, this%cutoffs)
    this%maxCutoff = maxval(this%cutoffs)

    allocate(this%nNeigh(this%nSpecies, this%nAtoms))
    allocate(this%nNeighMax(this%nAtoms))
    allocate(this%chargesPerAtom(this%nAtoms))
    allocate(this%chargesPerShell(this%mShells, this%nAtoms))
    allocate(this%shift1(this%mShells, this%nAtoms))
    allocate(this%shift2(this%mShells, this%nAtoms))
    allocate(this%shift3(this%nAtoms))
    allocate(this%gamma3ab(this%mShells, this%mShells, 0, this%nAtoms))
    allocate(this%gamma3ba(this%mShells, this%mShells, 0, this%nAtoms))
    allocate(this%damped(this%nSpecies))
    this%damped = inp%damped
    this%dampExp = inp%dampExp

  end subroutine ThirdOrder_init


  !> Returns real space cutoff.
  function getCutoff(this) result(cutoff)

    !> Instance
    class(TThirdOrder), intent(inout) :: this

    !> Cutoff
    real(dp) :: cutoff

    cutoff = this%maxCutoff

  end function getCutoff


  !> Updates data structures if there are changed coordinates for the instance.
  subroutine updateCoords(this, neighList, species)

    !> Instance.
    class(TThirdOrder), intent(inout) :: this

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighList

    !> Species for all atoms, shape: [nAllAtom].
    integer, intent(in) :: species(:)

    integer :: iNeigh, iAt1, iAt2, iSp1, iSp2, iSh1, iSh2
    logical :: damping
    real(dp) :: rr

    this%nNeigh(:,:) = 0
    do iAt1 = 1, this%nAtoms
      iSp1 = species(iAt1)
      do iSp2 = 1, this%nSpecies
        this%nNeigh(iSp2, iAt1) = getNrOfNeighbours(neighList, this%cutoffs(iSp2, iSp1), iAt1)
      end do
    end do
    this%nNeighMax = maxval(this%nNeigh, dim=1)

    if (size(this%gamma3ab, dim=3) < maxval(this%nNeighMax) + 1) then
      deallocate(this%gamma3ab)
      deallocate(this%gamma3ba)
      allocate(this%gamma3ab(this%mShells, this%mShells, 0:maxval(this%nNeighMax), this%nAtoms))
      allocate(this%gamma3ba(this%mShells, this%mShells, 0:maxval(this%nNeighMax), this%nAtoms))
    end if
    this%gamma3ab(:,:,:,:) = 0.0_dp
    this%gamma3ba(:,:,:,:) = 0.0_dp
    do iAt1 = 1, this%nAtoms
      iSp1 = species(iAt1)
      do iNeigh = 0, this%nNeighMax(iAt1)
        iAt2 = neighList%iNeighbour(iNeigh, iAt1)
        iSp2 = species(iAt2)
        if (iNeigh <= this%nNeigh(iSp2, iAt1)) then
          rr = sqrt(neighList%neighDist2(iNeigh, iAt1))
          damping = this%damped(iSp1) .or. this%damped(iSp2)
          do iSh1 = 1, this%nShells(iSp1)
            do iSh2 = 1, this%nShells(iSp2)
              this%gamma3ab(iSh2, iSh1, iNeigh, iAt1) = gamma3(this%UU(iSh1, iSp1),&
                  & this%UU(iSh2, iSp2), this%dUdQ(iSh1, iSp1), rr, damping, this%dampExp)
              this%gamma3ba(iSh2, iSh1, iNeigh, iAt1) = gamma3(this%UU(iSh2, iSp2),&
                  & this%UU(iSh1, iSp1), this%dUdQ(iSh2, iSp2), rr, damping, this%dampExp)
            end do
          end do
        end if
      end do
    end do

  end subroutine updateCoords


  !> Updates with changed charges for the instance.
  subroutine updateCharges(this, species, neighList, qq, q0, img2CentCell, orb)

    !> Instance
    class(TThirdOrder), intent(inout) :: this

    !> Species, shape: [nAtom]
    integer, intent(in) :: species(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighList

    !> Orbital charges.
    real(dp), intent(in) :: qq(:,:,:)

    !> Reference orbital charges.
    real(dp), intent(in) :: q0(:,:,:)

    !> Mapping on atoms in central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    real(dp), allocatable :: chargesPerShell(:,:)
    integer :: iAt1, iAt2f, iSp1, iSp2, iSh1, iSh2, iNeigh

    @:ASSERT(size(species) == this%nAtoms)
    @:ASSERT(size(qq, dim=2) == this%nAtoms)
    @:ASSERT(size(q0, dim=2) == this%nAtoms)

    if (this%shellResolved) then
      call getSummedCharges(species, orb, qq, q0, dQAtom=this%chargesPerAtom,&
          & dQShell=this%chargesPerShell)
    else
      ! First (only) component of this%chargesPerShell contains atomic charge
      allocate(chargesPerShell(this%mShellsReal, this%nAtoms))
      call getSummedCharges(species, orb, qq, q0, dQAtom=this%chargesPerAtom,&
          & dQShell=chargesPerShell)
      this%chargesPerShell(1,:) = sum(chargesPerShell, dim=1)
    end if

    this%shift1(:,:) = 0.0_dp
    this%shift2(:,:) = 0.0_dp
    this%shift3(:) = 0.0_dp

    do iAt1 = 1, this%nAtoms
      iSp1 = species(iAt1)
      do iNeigh = 0, this%nNeighMax(iAt1)
        iAt2f = img2CentCell(neighList%iNeighbour(iNeigh, iAt1))
        iSp2 = species(iAt2f)
        do iSh1 = 1, this%nShells(iSp1)
          do iSh2 = 1, this%nShells(iSp2)
            this%shift1(iSh1, iAt1) = this%shift1(iSh1, iAt1)&
                & + this%gamma3ab(iSh2, iSh1, iNeigh, iAt1) * this%chargesPerShell(iSh2, iAt2f)&
                & * this%chargesPerAtom(iAt1)
            this%shift2(iSh1, iAt1) = this%shift2(iSh1, iAt1)&
                & + this%gamma3ba(iSh2, iSh1, iNeigh, iAt1) * this%chargesPerShell(iSh2, iAt2f)&
                & * this%chargesPerAtom(iAt2f)
            this%shift3(iAt1) = this%shift3(iAt1)&
                & + this%gamma3ab(iSh2, iSh1, iNeigh, iAt1) * this%chargesPerShell(iSh2, iAt2f)&
                & * this%chargesPerShell(iSh1, iAt1)
            if (iAt2f /= iAt1) then
              this%shift1(iSh2, iAt2f) = this%shift1(iSh2, iAt2f)&
                  & + this%gamma3ba(iSh2, iSh1, iNeigh, iAt1) * this%chargesPerShell(iSh1, iAt1)&
                  & * this%chargesPerAtom(iAt2f)
              this%shift2(iSh2, iAt2f) = this%shift2(iSh2, iAt2f)&
                  & + this%gamma3ab(iSh2, iSh1, iNeigh, iAt1) * this%chargesPerShell(iSh1, iAt1)&
                  & * this%chargesPerAtom(iAt1)
              this%shift3(iAt2f) = this%shift3(iAt2f)&
                  & + this%gamma3ba(iSh2, iSh1, iNeigh, iAt1) * this%chargesPerShell(iSh1, iAt1)&
                  & * this%chargesPerShell(iSh2, iAt2f)
            end if
          end do
        end do
      end do
    end do
    this%shift1(:,:) = 1.0_dp / 3.0_dp * this%shift1
    this%shift2(:,:) = 1.0_dp / 3.0_dp * this%shift2
    this%shift3(:) = 1.0_dp / 3.0_dp * this%shift3

  end subroutine updateCharges


  !> Returns shifts per atom.
  subroutine getShifts(this, shiftPerAtom, shiftPerShell)

    !> Instance.
    class(TThirdOrder), intent(inout) :: this

    !> Shift per atom.
    real(dp), intent(out) :: shiftPerAtom(:)

    !> Shift per shell.
    real(dp), intent(out) :: shiftPerShell(:,:)

    @:ASSERT(size(shiftPerAtom) == this%nAtoms)
    @:ASSERT(size(shiftPerShell, dim=1) == this%mShellsReal)

    if (this%shellResolved) then
      shiftPerAtom(:) = this%shift3
      shiftPerShell(:,:) = this%shift1 + this%shift2
    else
      shiftPerAtom(:) = this%shift1(1,:) + this%shift2(1,:) + this%shift3
      shiftPerShell(:,:) = 0.0_dp
    end if

  end subroutine getShifts


  !> Returns energy per atom.
  subroutine getEnergyPerAtom(this, energyPerAtom)
    class(TThirdOrder), intent(inout) :: this
    real(dp), intent(out) :: energyPerAtom(:)

    @:ASSERT(size(energyPerAtom) == this%nAtoms)

    energyPerAtom(:) = (1.0_dp / 3.0_dp) * (&
        & sum((this%shift1 + this%shift2) * this%chargesPerShell, dim=1)&
        & + this%shift3 * this%chargesPerAtom)

  end subroutine getEnergyPerAtom


  !> Returns the energy per atom for linearized 3rd order Hamiltonian.
  !> Note: When using the linear XLBOMD form, charges should not be updated via updateCharges after
  !> the diagonalization, so that the shift vectors remain the ones built with the input
  !> charges. However, since for calculating energy/forces, the output charges are needed, they must
  !> be passed explicitely here.
  subroutine getEnergyPerAtomXlbomd(this, qOut, q0, species, orb, energyPerAtom)

    !> Instance.
    class(TThirdOrder), intent(inout) :: this

    !> Output populations determined after the diagonalization.
    real(dp), intent(in) :: qOut(:,:,:)

    !> Reference populations.
    real(dp), intent(in) :: q0(:,:,:)

    !> Species of each atom.
    integer, intent(in) :: species(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Energy per atom for linearized case.
    real(dp), intent(out) :: energyPerAtom(:)

    real(dp), allocatable :: qOutAtom(:), qOutShell(:,:), qOutShellTmp(:,:)

    allocate(qOutAtom(this%nAtoms))
    allocate(qOutShell(this%mShells, this%nAtoms))
    if (this%shellResolved) then
      call getSummedCharges(species, orb, qOut, q0, dQAtom=qOutAtom, dQShell=qOutShell)
    else
      ! First (only) component of qOutShell contains atomic charge
      allocate(qOutShellTmp(this%mShellsReal, this%nAtoms))
      call getSummedCharges(species, orb, qOut, q0, dQAtom=qOutAtom, dQShell=qOutShellTmp)
      qOutShell(1,:) = sum(qOutShellTmp, dim=1)
    end if
    energyPerAtom(:) = sum(this%shift1 * qOutShell, dim=1)&
        & + sum(this%shift2 * (qOutShell - this%chargesPerShell), dim=1)&
        & + this%shift3 * (qOutAtom - this%chargesPerAtom)

  end subroutine getEnergyPerAtomXlbomd


  !> Add gradient component resulting from the derivative of the potential.
  subroutine addGradientDc(this, neighList, species, coords, img2CentCell, derivs)

    !> Instance.
    class(TThirdOrder), intent(inout) :: this

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighList

    !> Specie for each atom.
    integer, intent(in) :: species(:)

    !> Coordinate of each atom.
    real(dp), intent(in) :: coords(:,:)

    !> Mapping of atoms to cetnral cell.
    integer, intent(in) :: img2CentCell(:)

    !> Gradient on exit.
    real(dp), intent(inout) :: derivs(:,:)

    integer :: iAt1, iAt2, iAt2f, iSp1, iSp2, iSh1, iSh2, iNeigh
    real(dp) :: rab, tmp, tmp3(3)
    logical :: damping

    do iAt1 = 1, this%nAtoms
      iSp1 = species(iAt1)
      do iNeigh = 1, this%nNeighMax(iAt1)
        iAt2 = neighList%iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        if (iAt1 == iAt2f .or. iNeigh > this%nNeigh(iSp2, iAt1)) then
          cycle
        end if
        rab = sqrt(neighList%neighDist2(iNeigh, iAt1))
        damping = this%damped(iSp1) .or. this%damped(iSp2)
        tmp = 0.0_dp
        do iSh1 = 1, this%nShells(iSp1)
          do iSh2 = 1, this%nShells(iSp2)
            tmp = tmp + this%chargesPerShell(iSh1, iAt1) * this%chargesPerShell(iSh2, iAt2f)&
                & * (this%chargesPerAtom(iAt1)&
                & * gamma3pR(this%UU(iSh1, iSp1), this%UU(iSh2, iSp2),&
                & this%dUdQ(iSh1, iSp1), rab, damping, this%dampExp)&
                & + this%chargesPerAtom(iAt2f)&
                & * gamma3pR(this%UU(iSh2, iSp2), this%UU(iSh1, iSp1),&
                & this%dUdQ(iSh2, iSp2), rab, damping, this%dampExp))
          end do
        end do
        tmp3(:) = tmp / (3.0_dp * rab) * (coords(:, iAt1) - coords(:, iAt2))
        derivs(:, iAt1) = derivs(:, iAt1) + tmp3
        derivs(:, iAt2f) = derivs(:, iAt2f) - tmp3
      end do
    end do

  end subroutine addGradientDc


  !> Add gradient component resulting from the derivative of the potential for the linearized
  !> (XLBOMD) case.
  subroutine addGradientDcXlbomd(this, neighList, species, coords, img2CentCell, qOut, q0, orb,&
      & derivs)

    !> Instance.
    class(TThirdOrder), intent(inout) :: this

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighList

    !> Specie for each atom.
    integer, intent(in) :: species(:)

    !> Coordinate of each atom.
    real(dp), intent(in) :: coords(:,:)

    !> Mapping of atoms to cetnral cell.
    integer, intent(in) :: img2CentCell(:)

    !> Output populations determined after the diagonalization.
    real(dp), intent(in) :: qOut(:,:,:)

    !> Reference populations.
    real(dp), intent(in) :: q0(:,:,:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Modified gradient on exit.
    real(dp), intent(inout) :: derivs(:,:)

    real(dp), allocatable :: qOutAtom(:), qOutShell(:,:), qOutShellTmp(:,:)
    real(dp), allocatable :: qDiffShell(:,:), qDiffAtom(:)
    integer :: iAt1, iAt2, iAt2f, iSp1, iSp2, iSh1, iSh2, iNeigh
    real(dp) :: gammaDeriv1, gammaderiv2, rab, tmp
    real(dp) :: tmp3(3)
    logical :: damping

    allocate(qOutAtom(this%nAtoms))
    allocate(qOutShell(this%mShells, this%nAtoms))
    allocate(qDiffAtom(this%nAtoms))
    allocate(qDiffShell(this%mShells, this%nAtoms))
    if (this%shellResolved) then
      call getSummedCharges(species, orb, qOut, q0, dQAtom=qOutAtom, dQShell=qOutShell)
    else
      ! First (only) component of qOutShell contains atomic charge
      allocate(qOutShellTmp(this%mShellsReal, this%nAtoms))
      call getSummedCharges(species, orb, qOut, q0, dQAtom=qOutAtom, dQShell=qOutShellTmp)
      qOutShell(1,:) = sum(qOutShellTmp, dim=1)
    end if

    qDiffAtom(:) = qOutAtom - this%chargesPerAtom
    qDiffShell(:,:) = qOutShell - this%chargesPerShell
    do iAt1 = 1, this%nAtoms
      iSp1 = species(iAt1)
      do iNeigh = 1, this%nNeighMax(iAt1)
        iAt2 = neighList%iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        if (iAt1 == iAt2f .or. iNeigh > this%nNeigh(iSp2, iAt1)) then
          cycle
        end if
        rab = sqrt(neighList%neighDist2(iNeigh, iAt1))
        damping = this%damped(iSp1) .or. this%damped(iSp2)
        tmp = 0.0_dp
        do iSh1 = 1, this%nShells(iSp1)
          do iSh2 = 1, this%nShells(iSp2)
            gammaDeriv1 =&
                & gamma3pR(this%UU(iSh1, iSp1), this%UU(iSh2, iSp2),&
                & this%dUdQ(iSh1, iSp1), rab, damping, this%dampExp)
            gammaDeriv2 =&
                & gamma3pR(this%UU(iSh2, iSp2), this%UU(iSh1, iSp1),&
                & this%dUdQ(iSh2, iSp2), rab, damping, this%dampExp)
            tmp = tmp + gammaDeriv1&
                & * (this%chargesPerAtom(iAt1) * this%chargesPerShell(iSh2, iAt2f)&
                & * qOutShell(iSh1, iAt1)&
                & + this%chargesPerAtom(iAt1) * this%chargesPerShell(iSh1, iAt1)&
                & * qDiffShell(iSh2, iAt2f)&
                & + this%chargesPerShell(iSh1, iAt1) * this%chargesPerShell(iSh2, iAt2f)&
                & * qDiffAtom(iAt1))
            tmp = tmp + gammaDeriv2&
                & * (this%chargesPerAtom(iAt2f) * this%chargesPerShell(iSh1, iAt1)&
                & * qOutShell(iSh2, iAt2f)&
                & + this%chargesPerAtom(iAt2f) * this%chargesPerShell(iSh2, iAt2f)&
                & * qDiffShell(iSh1, iAt1)&
                & + this%chargesPerShell(iSh2, iAt2f) * this%chargesPerShell(iSh1, iAt1)&
                & * qDiffAtom(iAt2f))
          end do
        end do
        tmp3 = tmp / (3.0_dp * rab) * (coords(:, iAt1) - coords(:, iAt2))
        derivs(:, iAt1) = derivs(:, iAt1) + tmp3
        derivs(:, iAt2f) = derivs(:, iAt2f) - tmp3
      end do
    end do

  end subroutine addGradientDcXlbomd


  !> Add stress component resulting from the derivative of the potential.
  subroutine addStressDc(this, neighList, species, coords, img2CentCell, cellVol, st)

    !> Instance.
    class(TThirdOrder), intent(inout) :: this

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighList

    !> Specie for each atom.
    integer, intent(in) :: species(:)

    !> Coordinate of each atom.
    real(dp), intent(in) :: coords(:,:)

    !> Mapping of atoms to cetnral cell.
    integer, intent(in) :: img2CentCell(:)

    !> cell volume.
    real(dp), intent(in) :: cellVol

    !> Gradient on exit.
    real(dp), intent(inout) :: st(:,:)

    integer :: iAt1, iAt2, iAt2f, iSp1, iSp2, iSh1, iSh2, iNeigh, ii
    real(dp) :: rab, tmp, tmp3(3), stTmp(3,3), prefac, vect(3)
    logical :: damping

    stTmp(:,:) = 0.0_dp
    do iAt1 = 1, this%nAtoms
      iSp1 = species(iAt1)
      do iNeigh = 1, this%nNeighMax(iAt1)
        iAt2 = neighList%iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        if (iNeigh > this%nNeigh(iSp2, iAt1)) then
          cycle
        end if
        vect(:) = coords(:,iAt1) - coords(:,iAt2)
        if (iAt1 == iAt2f) then
          prefac = 0.5_dp
        else
          prefac = 1.0_dp
        end if
        rab = sqrt(neighList%neighDist2(iNeigh, iAt1))
        damping = this%damped(iSp1) .or. this%damped(iSp2)
        tmp = 0.0_dp
        do iSh1 = 1, this%nShells(iSp1)
          do iSh2 = 1, this%nShells(iSp2)
            tmp = tmp + this%chargesPerShell(iSh1, iAt1) * this%chargesPerShell(iSh2, iAt2f)&
                & * (this%chargesPerAtom(iAt1)&
                & * gamma3pR(this%UU(iSh1, iSp1), this%UU(iSh2, iSp2),&
                & this%dUdQ(iSh1, iSp1), rab, damping, this%dampExp)&
                & + this%chargesPerAtom(iAt2f)&
                & * gamma3pR(this%UU(iSh2, iSp2), this%UU(iSh1, iSp1),&
                & this%dUdQ(iSh2, iSp2), rab, damping, this%dampExp))
          end do
        end do
        tmp3(:) = tmp / (3.0_dp * rab) * (coords(:, iAt1) - coords(:, iAt2))
        do ii = 1, 3
          stTmp(:, ii) = stTmp(:, ii) + prefac * tmp3(:) * vect(ii)
        end do
      end do
    end do

    st = st - stTmp / cellVol

  end subroutine addStressDc


! Private routines


  !> calculate short range cut off distance
  subroutine calcCutoffs(hubbUs, nShells, cutoffs)

    !> Hubard U values
    real(dp), intent(in) :: hubbUs(:,:)

    !> Shells on atoms
    integer, intent(in) :: nShells(:)

    !> resulting cutoff distances
    real(dp), intent(out) :: cutoffs(:,:)

    integer :: nSpecies
    real(dp) :: cutoff
    real(dp), allocatable :: minUs(:)
    integer :: iSp1, iSp2

    nSpecies = size(cutoffs, dim=1)
    allocate(minUs(nSpecies))
    do iSp1 = 1, nSpecies
      minUs(iSp1) = minval(hubbUs(1:nShells(iSp1), iSp1))
    end do
    do iSp1 = 1, nSpecies
      do iSp2 = iSp1, nSpecies
        cutoff = expGammaCutoff(minUs(iSp2), minUs(iSp1))
        cutoffs(iSp2, iSp1) = cutoff
        cutoffs(iSp1, iSp2) = cutoff
      end do
    end do

  end subroutine calcCutoffs


  !> Gamma_AB = dgamma_AB/dUa * (dUa/dQa)
  function gamma3(Ua, Ub, dUa, rab, damping, xi) result(res)
    real(dp), intent(in) :: Ua
    real(dp), intent(in) :: Ub
    real(dp), intent(in) :: dUa
    real(dp), intent(in) :: rab
    logical, intent(in) :: damping
    real(dp), intent(in) :: xi
    real(dp) :: res

    res = gamma2pU(Ua, Ub, rab, damping, xi) * dUa

  end function gamma3


  !> dGamma_AB/dr
  function gamma3pR(Ua, Ub, dUa, rab, damping, xi) result(res)
    real(dp), intent(in) :: Ua, Ub, dUa, rab
    logical, intent(in) :: damping
    real(dp), intent(in) :: xi
    real(dp) :: res

    res = gamma2pUpR(Ua, Ub, rab, damping, xi) * dUa

  end function gamma3pR


  !> dgamma_AB/dUa
  !> Sign convention: routine delivers dgamma_AB/dUa with the right sign.
  !> Energy contributions must be therefore summed with *positive* sign.
  function gamma2pU(Ua, Ub, rab, damping, xi) result(res)
    real(dp), intent(in) :: Ua, Ub, rab
    logical, intent(in) :: damping
    real(dp), intent(in) :: xi
    real(dp) :: res

    real(dp) :: tauA, tauB, tau, uu

    tauA = 3.2_dp * Ua      ! 16/5 * Ua
    tauB = 3.2_dp * Ub

    if (rab < tolSameDist) then
      if (abs(Ua - Ub) < minHubDiff) then
        ! Limiting case for dG/dU with Ua=Ub and Rab = 0
        res = 0.5_dp
      else
        res = dGdUr0(tauA, tauB)
      end if
    else if (abs(Ua - Ub) < minHubDiff) then
      tau = 0.5_dp * (tauA + tauB)
      res = -3.2_dp * shortpT_2(tau, rab)
      if (damping) then
        uu = 0.5_dp * (Ua + Ub)
        res = res * hh(uu, uu, rab, xi) - short_2(tau, rab) * hpU(uu, uu, rab, xi)
      end if
    else
      res = -3.2_dp * shortpT_1(tauA, tauB, rab)
      if (damping) then
        res = res * hh(Ua, Ub, rab, xi) - short_1(tauA, tauB, rab) * hpU(Ua, Ub, rab, xi)
      end if
    end if

  end function gamma2pU


  !> d^2gamma_AB/dUa*dr
  !> Sign convention: routine delivers d^2gamma_AB/dUa*dr with the right sign.
  !> Gradient contributions must be therefore summed with *positive* sign.
  function gamma2pUpR(Ua, Ub, rab, damping, xi) result(res)
    real(dp), intent(in) :: Ua, Ub, rab
    logical, intent(in) :: damping
    real(dp), intent(in) :: xi
    real(dp) :: res

    real(dp) :: tauA, tauB, tau, uu

    tauA = 3.2_dp * Ua
    tauB = 3.2_dp * Ub

    if (rab < tolSameDist) then
      res = 0.0_dp
    else if (abs(Ua - Ub) < minHubDiff) then
      tau = 0.5_dp * (tauA + tauB)
      res = -3.2_dp * shortpTpR_2(tau, rab)
      if (damping) then
        uu = 0.5_dp * (Ua + Ub)
        res = res * hh(uu, uu, rab, xi)&
            & - 3.2_dp * shortpT_2(tau, rab) * hpR(uu, uu, rab, xi)&
            & - shortpR_2(tau, rab) * hpU(uu, uu, rab, xi)&
            & - short_2(tau, rab) * hpUpR(uu, uu, rab, xi)
      end if
    else
      res = -3.2_dp *shortpTpR_1(tauA, tauB, rab)
      if (damping) then
        res = res * hh(Ua, Ub, rab, xi)&
            & - 3.2_dp * shortpT_1(tauA, tauB, rab) * hpR(Ua, Ub, rab, xi)&
            & - shortpR_1(tauA, tauB, rab) * hpU(Ua, Ub, rab, xi)&
            & - short_1(tauA, tauB, rab) * hpUpR(Ua, Ub, rab, xi)
      end if
    end if

  end function gamma2pUpR


  !> \frac{d\gamma}{dU_{l_a}} for r = 0
  !> Eq S7 in Gaus et al. (2015) JCTC 11:4205-4219, DOI: 10.1021/acs.jctc.5b00600
  function dGdUr0(tauA, tauB) result(res)
    real(dp), intent(in) :: tauA, tauB
    real(dp) :: res

    real(dp) :: invTauSum

    invTauSum = 1.0_dp / (tauA + tauB)

    res = 1.6_dp * invTauSum * (tauB&
      & + invTauSum * (-tauA * tauB&
      & + invTauSum * (2.0_dp * tauA * tauB**2&
      & + invTauSum * (-3.0_dp * tauA**2 * tauB**2))))

  end function dGdUr0


  !> S1(tauA,tauB,r): Short range SCC when tauA <> tauB and r <> 0
  function short_1(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = exp(-tauA * rab) * ff(tauA, tauB, rab) + exp(-tauB * rab) * ff(tauB, tauA, rab)

  end function short_1


  !> S2(tau,r), short range SCC when tauA = tauB = tau and r <> 0.
  function short_2(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = exp(-tau * rab) * gg(tau, rab)

  end function short_2


  !> dS1(tauA,tauB,r)/dtauA
  function shortpT_1(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = exp(-tauA * rab) * (fpT1(tauA, tauB, rab) - rab * ff(tauA, tauB, rab))&
        & + exp(-tauB * rab) * fpT2(tauB, tauA, rab)

  end function shortpT_1


  !> dS2(tauA,tauB,r)/dtauA
  function shortpT_2(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = exp(-tau * rab) * (gpT(tau, rab) - rab * gg(tau, rab))

  end function shortpT_2


  !> dS1(tauA,tauB,r)/dr
  function shortpR_1(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = exp(-tauA * rab) * (fpR(tauA, tauB, rab) - tauA *  ff(tauA, tauB, rab))&
        & + exp(-tauB * rab) * (fpR(tauB, tauA, rab) - tauB * ff(tauB, tauA, rab))

  end function shortpR_1


  !> dS2(tauA,tauB,r)/dr
  function shortpR_2(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = exp(-tau * rab) * (gpR(tau, rab) - tau * gg(tau, rab))

  end function shortpR_2


  !> d^2S1(tauA,tauB,r)/dtauA*dr
  function shortpTpR_1(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = exp(-tauA * rab) * (ff(tauA, tauB, rab) * (tauA * rab - 1.0_dp)&
        & - tauA * fpT1(tauA, tauB, rab) + fpT1pR(tauA, tauB, rab) - rab * fpR(tauA, tauB, rab))&
        & + exp(-tauB * rab) * (fpT2pR(tauB, tauA, rab) - tauB * fpT2(tauB, tauA, rab))

  end function shortpTpR_1


  !> d^2S2(tau,r)/dtau*dr
  function shortpTpR_2(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = exp(-tau * rab) * ((tau * rab - 1.0_dp) * gg(tau, rab) - tau * gpT(tau, rab)&
        & + gpTpR(tau, rab) - rab * gpR(tau, rab))

  end function shortpTpR_2


  !> f(tauA,tauB,r)
  function ff(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = 0.5_dp * tauA * tauB**4 / (tauA**2 - tauB**2)**2&
        & - (tauB**6 - 3.0_dp * tauA**2 * tauB**4) / ((tauA**2 - tauB**2)**3 * rab)

  end function ff


  !> df(tauA,tauB,r)/dtauA
  function fpT1(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = -0.5_dp * (tauB**6 + 3.0_dp * tauA**2 * tauB**4) / (tauA**2 - tauB**2)**3&
        & - 12.0_dp * tauA**3 * tauB**4 / ((tauA**2 - tauB**2)**4 * rab)

  end function fpT1


  !> df(tauB,tauA,rab)/dtauA
  function fpT2(tauB, tauA , rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = 2.0_dp * tauB**3 * tauA**3 / (tauB**2 - tauA**2)**3&
        & + 12.0_dp * tauB**4 * tauA**3 / ((tauB**2 - tauA**2)**4 * rab)

  end function fpT2


  !> df(tauA, tauB,r)/dr
  function fpR(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = (tauB**6 - 3.0_dp * tauB**4 * tauA**2) / (rab**2 * (tauA**2 - tauB**2)**3)

  end function fpR


  !> d^2f(tauA,tauB,r)/dtauA*dr
  function fpT1pR(tauA, tauB, rab) result(res)
    real(dp), intent(in) :: tauA, tauB, rab
    real(dp) :: res

    res = 12.0_dp * tauA**3 * tauB**4 / (rab**2 * (tauA**2 - tauB**2)**4)

  end function fpT1pR


  !> d^2f(tauB,tauA,r)/dtauA*dr
  function fpT2pR(tauB, tauA, rab) result(res)
    real(dp), intent(in) :: tauB, tauA, rab
    real(dp) :: res

    res = -12.0_dp * tauA**3 * tauB**4 / (rab**2 * (tauA**2 - tauB**2)**4)

  end function fpT2pR


  !> g(tau,r)
  function gg(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = 1.0_dp / (48.0_dp * rab) * (48.0_dp + 33.0_dp * tau * rab + 9.0_dp * tau**2 * rab**2&
        & + tau**3 * rab**3)

  end function gg


  !> dg(tau,rab)/dtau
  function gpT(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = 1.0_dp / 48.0_dp * (33.0_dp + 18.0_dp * tau * rab + 3.0_dp * tau**2 * rab**2)

  end function gpT


  !> dg(tau,r)/dr
  function gpR(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = (-48.0_dp + 9.0_dp * (tau * rab)**2 + 2.0_dp * (tau * rab)**3) / (48.0_dp * rab**2)

  end function gpR


  !> d^2g(tau,r)/dtau*dr
  function gpTpR(tau, rab) result(res)
    real(dp), intent(in) :: tau, rab
    real(dp) :: res

    res = (3.0_dp * tau + tau**2 * rab) / 8.0_dp

  end function gpTpR


  !> Damping: h(Ua,Ub)
  function hh(Ua, Ub, rab, xi) result(res)
    real(dp), intent(in) :: Ua, Ub, rab, xi
    real(dp) :: res

    res = exp(-(0.5_dp * (Ua + Ub))**xi * rab**2)

  end function hh


  !> dh(Ua,Ub)/dUa
  function hpU(Ua, Ub, rab, xi) result(res)
    real(dp), intent(in) :: Ua, Ub, rab, xi
    real(dp) :: res

    res = -0.5_dp * xi * rab**2 * (0.5_dp * (Ua + Ub))**(xi - 1.0_dp) * hh(Ua, Ub, rab, xi)

  end function hpU


  !> dh(Ua,Ub)/dr
  function hpR(Ua, Ub, rab, xi) result(res)
    real(dp), intent(in) :: Ua, Ub, rab, xi
    real(dp) :: res

    res = -2.0_dp * rab * (0.5_dp * (Ua + Ub))**xi * hh(Ua, Ub, rab, xi)

  end function hpR


  !> dh(Ua,Ub)/dUa*dr
  function hpUpR(Ua, Ub, rab, xi) result(res)
    real(dp), intent(in) :: Ua, Ub, rab, xi
    real(dp) :: res

    res = xi * rab * (0.5_dp * (Ua + Ub))**(xi - 1.0_dp)&
        & * (rab**2 * (0.5_dp * (Ua + Ub))**xi - 1.0_dp) * hh(Ua, Ub, rab, xi)

  end function hpUpR

end module dftbp_thirdorder
