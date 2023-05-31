!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> update the SCC hamiltonian
module dftbp_dftb_hamiltonian
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_environment, only : TEnvironment
  use dftbp_dftb_dftbplusu, only : TDftbU
  use dftbp_dftb_dispersions, only : TDispersionIface
  use dftbp_common_environment, only : TEnvironment
  use dftbp_dftb_extfields, only : TEField
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_potentials, only : TPotentials
  use dftbp_dftb_scc, only : TScc
  use dftbp_dftb_shift, only : addShift, totalShift, addOnSiteShift, addAtomicMultipoleShift
  use dftbp_dftb_spin, only : getSpinShift
  use dftbp_dftb_spinorbit, only : getDualSpinOrbitShift
  use dftbp_dftb_thirdorder, only : TThirdOrder
  use dftbp_extlibs_tblite, only : TTBLite
  use dftbp_io_message, only : error
  use dftbp_solvation_solvation, only : TSolvation
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_integral, only : TIntegral
  use dftbp_type_multipole, only : TMultipole
  implicit none

  private
  public :: resetExternalPotentials, getSccHamiltonian, mergeExternalPotentials
  public :: resetInternalPotentials, addChargePotentials
  public :: addBlockChargePotentials, TRefExtPot

  !> Container for external potentials
  type :: TRefExtPot
    real(dp), allocatable :: atomPot(:,:)
    real(dp), allocatable :: shellPot(:,:,:)
    real(dp), allocatable :: potGrad(:,:)
  end type TRefExtPot

contains


  !> Sets the external potential components to zero
  subroutine resetExternalPotentials(refExtPot, potential)

    !> Reference external potential (usually set via API)
    type(TRefExtPot), intent(in) :: refExtPot

    !> Potential contributions
    type(TPotentials), intent(inout) :: potential

    if (allocated(refExtPot%atomPot)) then
      potential%extAtom(:,:) = refExtPot%atomPot
    else
      potential%extAtom(:,:) = 0.0_dp
    end if
    if (allocated(refExtPot%shellPot)) then
      potential%extShell(:,:,:) = refExtPot%shellPot
    else
      potential%extShell(:,:,:) = 0.0_dp
    end if
    potential%extBlock(:,:,:,:) = 0.0_dp
    if (allocated(refExtPot%potGrad)) then
      potential%extGrad(:,:) = refExtPot%potGrad
    else
      potential%extGrad(:,:) = 0.0_dp
    end if
    if (allocated(potential%extOnSiteAtom)) then
      potential%extOnSiteAtom(:,:) = 0.0_dp
    end if
    if (allocated(potential%extDipoleAtom)) then
      potential%extDipoleAtom(:, :) = 0.0_dp
    end if

  end subroutine resetExternalPotentials


  !> Merges atomic and shell resolved external potentials into blocked one
  subroutine mergeExternalPotentials(orb, species, potential)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species for atoms
    integer, intent(in) :: species(:)

    !> Potential energy contributions
    type(TPotentials), intent(inout) :: potential

    call totalShift(potential%extShell, potential%extAtom, orb, species)
    call totalShift(potential%extBlock, potential%extShell, orb, species)

  end subroutine mergeExternalPotentials


  !> Reset internal potential related quantities
  subroutine resetInternalPotentials(tDualSpinOrbit, xi, orb, species, potential)

    !> Is dual spin orbit being used (block potentials)
    logical, intent(in) :: tDualSpinOrbit

    !> Spin orbit constants if required
    real(dp), allocatable, intent(in) :: xi(:,:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> chemical species
    integer, intent(in) :: species(:)

    !> potentials in the system
    type(TPotentials), intent(inout) :: potential

    @:ASSERT(.not. tDualSpinOrbit .or. allocated(xi))

    potential%intAtom(:,:) = 0.0_dp
    potential%intShell(:,:,:) = 0.0_dp
    potential%intBlock(:,:,:,:) = 0.0_dp
    potential%orbitalBlock(:,:,:,:) = 0.0_dp
    potential%iOrbitalBlock(:,:,:,:) = 0.0_dp
    if (tDualSpinOrbit) then
      call getDualSpinOrbitShift(potential%iOrbitalBlock, xi, orb, species)
    end if
    if (allocated(potential%intOnSiteAtom)) then
      potential%intOnSiteAtom(:,:) = 0.0_dp
    end if
    if (allocated(potential%dipoleAtom)) then
      potential%dipoleAtom(:,:) = 0.0_dp
    end if
    if (allocated(potential%quadrupoleAtom)) then
      potential%quadrupoleAtom(:,:) = 0.0_dp
    end if

  end subroutine resetInternalPotentials


  !> Add potentials coming from electrostatics (in various boundary conditions and models), possibly
  !> spin, and where relevant dispersion
  subroutine addChargePotentials(env, sccCalc, tblite, updateScc, qInput, q0, chargePerShell,&
      & orb, multipole, species, neighbourList, img2CentCell, spinW, solvation, thirdOrd,&
      & dispersion, potential)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> Library interface handler
    type(TTBLite), intent(inout), allocatable :: tblite

    !> Whether the charges in the scc calculator should be updated before obtaining the potential
    logical, intent(in) :: updateScc

    !> Input atomic populations
    real(dp), intent(in) :: qInput(:,:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> charges per atomic shell
    real(dp), intent(in) :: chargePerShell(:,:,:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Multipole information
    type(TMultipole), intent(in) :: multipole

    !> species of all atoms
    integer, target, intent(in) :: species(:)

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> Solvation mode
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Potentials acting
    type(TPotentials), intent(inout) :: potential

    !> Dispersion interactions object
    class(TDispersionIface), allocatable, intent(inout) :: dispersion

    ! local variables
    real(dp), allocatable :: atomPot(:,:)
    real(dp), allocatable :: shellPot(:,:,:)
    real(dp), allocatable :: dipPot(:,:), quadPot(:,:)
    integer, pointer :: pSpecies0(:)
    integer :: nAtom, nSpin

    nAtom = size(qInput, dim=2)
    nSpin = size(qInput, dim=3)
    pSpecies0 => species(1:nAtom)

    allocate(atomPot(nAtom, nSpin), source=0.0_dp)
    allocate(shellPot(orb%mShell, nAtom, nSpin), source=0.0_dp)
    if (allocated(potential%dipoleAtom)) then
      allocate(dipPot, mold=potential%dipoleAtom)
      dipPot(:, :) = 0.0_dp
    end if
    if (allocated(potential%quadrupoleAtom)) then
      allocate(quadPot, mold=potential%quadrupoleAtom)
      quadPot(:, :) = 0.0_dp
    end if

    if (allocated(sccCalc)) then
      if (updateScc) then
        call sccCalc%updateCharges(env, qInput, orb, species, q0)
      end if
      call sccCalc%updateShifts(env, orb, species, neighbourList%iNeighbour, img2CentCell)
      call sccCalc%getShiftPerAtom(atomPot(:,1))
      call sccCalc%getShiftPerL(shellPot(:,:,1))

      if (allocated(potential%coulombShell)) then
        ! need to retain the just electrostatic contributions to the potential for a contact
        ! calculation or similar
        potential%coulombShell(:,:,:) = shellPot(:,:,1:1)
        call totalShift(potential%coulombShell, atomPot(:,1:1), orb, species)
      end if
    end if

    if (allocated(dispersion)) then
      call dispersion%updateCharges(env, pSpecies0, neighbourList, qInput, q0, img2CentCell, orb)
      call dispersion%addPotential(atomPot(:,1))
    end if

    potential%intAtom(:,1) = potential%intAtom(:,1) + atomPot(:,1)
    potential%intShell(:,:,1) = potential%intShell(:,:,1) + shellPot(:,:,1)

    if (allocated(tblite)) then
      call tblite%updateCharges(env, species, neighbourList, qInput, q0, &
          & multipole%dipoleAtom, multipole%quadrupoleAtom, img2CentCell, orb)
      call tblite%getShifts(atomPot(:,1), shellPot(:,:,1), dipPot, quadPot)
      potential%intAtom(:,1) = potential%intAtom(:,1) + atomPot(:,1)
      potential%intShell(:,:,1) = potential%intShell(:,:,1) + shellPot(:,:,1)
      if (allocated(potential%dipoleAtom)) then
        potential%dipoleAtom(:,:) = potential%dipoleAtom + dipPot
      end if
      if (allocated(potential%quadrupoleAtom)) then
        potential%quadrupoleAtom(:,:) = potential%quadrupoleAtom + quadPot
      end if
    end if

    if (allocated(thirdOrd)) then
      call thirdOrd%updateCharges(pSpecies0, neighbourList, qInput, q0, img2CentCell, orb)
      call thirdOrd%getShifts(atomPot(:,1), shellPot(:,:,1))
      potential%intAtom(:,1) = potential%intAtom(:,1) + atomPot(:,1)
      potential%intShell(:,:,1) = potential%intShell(:,:,1) + shellPot(:,:,1)
    end if

    if (allocated(solvation)) then
      call solvation%updateCharges(env, pSpecies0, neighbourList, qInput, q0, img2CentCell, orb)
      call solvation%getShifts(atomPot(:,1), shellPot(:,:,1))
      potential%intAtom(:,1) = potential%intAtom(:,1) + atomPot(:,1)
      potential%intShell(:,:,1) = potential%intShell(:,:,1) + shellPot(:,:,1)
    end if

    if (nSpin /= 1 .and. allocated(spinW)) then
      shellPot(:,:,:) = 0.0_dp
      call getSpinShift(shellPot(:,:,2:), chargePerShell(:,:,2:), species, orb, spinW)
      potential%intShell(:,:,2:) = potential%intShell(:,:,2:) + shellPot(:,:,2:)
    end if

    call totalShift(potential%intShell, potential%intAtom, orb, species)
    call totalShift(potential%intBlock, potential%intShell, orb, species)

  end subroutine addChargePotentials


  !> Add potentials coming from on-site block of the dual density matrix.
  subroutine addBlockChargePotentials(qBlockIn, qiBlockIn, dftbU, tImHam, species, orb, potential)

    !> block input charges
    real(dp), allocatable, intent(in) :: qBlockIn(:,:,:,:)

    !> imaginary part
    real(dp), allocatable, intent(in) :: qiBlockIn(:,:,:,:)

    !> is this a +U calculation
    type(TDftbU), intent(in), allocatable :: dftbU

    !> does the hamiltonian have an imaginary part in real space?
    logical, intent(in) :: tImHam

    !> chemical species of all atoms
    integer, intent(in) :: species(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> potentials acting in system
    type(TPotentials), intent(inout) :: potential

    if (allocated(dftbU)) then
      if (tImHam) then
        call dftbU%getDftbUShift(potential%orbitalBlock, potential%iorbitalBlock, qBlockIn,&
            & qiBlockIn, species,orb)
      else
        call dftbU%getDftbUShift(potential%orbitalBlock, qBlockIn, species, orb)
      end if
      potential%intBlock = potential%intBlock + potential%orbitalBlock
    end if

  end subroutine addBlockChargePotentials


  !> Returns the Hamiltonian for the given scc iteration
  subroutine getSccHamiltonian(env, H0, ints, nNeighbourSK, neighbourList, species, orb, iSparseStart,&
      & img2CentCell, potential, isREKS, ham, iHam)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> non-SCC hamiltonian (sparse)
    real(dp), intent(in) :: H0(:)

    !> Integral container
    type(TIntegral), intent(in) :: ints

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> list of atomic neighbours
    type(TNeighbourList), intent(in) :: neighbourList

    !> species of atoms
    integer, intent(in) :: species(:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> image atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> potential acting on sustem
    type(TPotentials), intent(in) :: potential

    !> Is this DFTB/SSR formalism
    logical, intent(in) :: isREKS

    !> resulting hamitonian (sparse)
    real(dp), intent(inout) :: ham(:,:)

    !> imaginary part of hamiltonian (if required, signalled by being allocated)
    real(dp), allocatable, intent(inout) :: iHam(:,:)

    integer :: nAtom
    real(dp), allocatable :: dipoleAtom(:, :)

    nAtom = size(orb%nOrbAtom)

    if (.not. isREKS) then
      ham(:,:) = 0.0_dp
    end if

    call addShift(env, ham, ints%overlap, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
        & iSparseStart, nAtom, img2CentCell, potential%intBlock, .not. isREKS)

    if (.not. isREKS) then
      ham(:,1) = ham(:,1) + h0
    end if

    if (allocated(potential%intOnSiteAtom)) then
      call addOnSiteShift(ham, ints%overlap, species, orb, iSparseStart, nAtom,&
          & potential%intOnSiteAtom)
    end if
    if (allocated(potential%extOnSiteAtom)) then
      call addOnSiteShift(ham, ints%overlap, species, orb, iSparseStart, nAtom,&
          & potential%extOnSiteAtom)
    end if

    if (allocated(potential%dipoleAtom)) then
      dipoleAtom = potential%dipoleAtom
      if (allocated(potential%extDipoleAtom)) then
        dipoleAtom(:, :) = dipoleAtom + potential%extDipoleAtom
      end if
      call addAtomicMultipoleShift(ham, ints%dipoleBra, ints%dipoleKet, nNeighbourSK, &
          & neighbourList%iNeighbour, species, orb, iSparseStart, nAtom, img2CentCell, &
          & dipoleAtom)
    end if

    if (allocated(potential%quadrupoleAtom)) then
      call addAtomicMultipoleShift(ham, ints%quadrupoleBra, ints%quadrupoleKet, nNeighbourSK, &
          & neighbourList%iNeighbour, species, orb, iSparseStart, nAtom, img2CentCell, &
          & potential%quadrupoleAtom)
    end if

    if (allocated(iHam)) then
      iHam(:,:) = 0.0_dp
      call addShift(env, iHam, ints%overlap, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
          & iSparseStart, nAtom, img2CentCell, potential%iorbitalBlock, .true.)
    end if

  end subroutine getSccHamiltonian


end module dftbp_dftb_hamiltonian
