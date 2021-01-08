!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Evaluate energies
module dftbp_getenergies
  use dftbp_accuracy, only : dp, lc
  use dftbp_assert
  use dftbp_energytypes, only : TEnergies
  use dftbp_populations
  use dftbp_commontypes, only : TOrbitals
  use dftbp_periodic, only : TNeighbourList
  use dftbp_potentials, only : TPotentials
  use dftbp_shift, only : add_shift, total_shift
  use dftbp_spin, only : getSpinShift
  use dftbp_spinorbit, only : getDualSpinOrbitShift, getDualSpinOrbitEnergy
  use dftbp_dftbplusu, only : TDftbU
  use dftbp_message, only : error
  use dftbp_thirdorder, only : TThirdOrder
  use dftbp_environment, only : TEnvironment
  use dftbp_scc, only : TScc
  use dftbp_rangeseparated, only : TRangeSepFunc
  use dftbp_qdepextpotproxy, only : TQDepExtPotProxy
  use dftbp_onsitecorrection
  use dftbp_dispiface
#:if WITH_MBD
  use dftbp_dispmbd, only: TDispMbd
#:endif
  use dftbp_solvation, only : TSolvation
  use dftbp_repcont
  use dftbp_repulsive
  use dftbp_reks, only : TReksCalc
  use dftbp_determinants, only : TDftbDeterminants, determinants
  implicit none

  private
  public :: calcEnergies, calcRepulsiveEnergy, calcDispersionEnergy, sumEnergies

contains


  !> Calculates various energy contribution that can potentially update for the same geometry
  subroutine calcEnergies(sccCalc, qOrb, q0, chargePerShell, species, tExtField, isXlbomd, dftbU,&
      & tDualSpinOrbit, rhoPrim, H0, orb, neighbourList, nNeighbourSK, img2CentCell, iSparseStart,&
      & cellVol, extPressure, TS, potential, energy, thirdOrd, solvation, rangeSep, reks,&
      & qDepExtPot, qBlock, qiBlock, xi, iAtInCentralRegion, tFixEf, Ef, onSiteElements)

    !> SCC module internal variables
    type(TScc), allocatable, intent(in) :: sccCalc

    !> Electrons in each atomic orbital
    real(dp), intent(in) :: qOrb(:,:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> electrons in each atomi shell
    real(dp), intent(in) :: chargePerShell(:,:,:)

    !> chemical species
    integer, intent(in) :: species(:)

    !> is an external electric field present
    logical, intent(in) :: tExtField

    !> Is the extended Lagrangian being used for MD
    logical, intent(in) :: isXlbomd

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> Is dual spin orbit being used
    logical, intent(in) :: tDualSpinOrbit

    !> density matrix in sparse storage
    real(dp), intent(in) :: rhoPRim(:,:)

    !> non-self-consistent hamiltonian
    real(dp), intent(in) :: H0(:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> neighbour list
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours within cut-off for each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> image to real atom mapping
    integer, intent(in) :: img2CentCell(:)

    !> index for sparse large matrices
    integer, intent(in) :: iSparseStart(:,:)

    !> unit cell volume
    real(dp), intent(in) :: cellVol

    !> external pressure
    real(dp), intent(in) :: extPressure

    !> electron entropy contribution
    real(dp), intent(in) :: TS(:)

    !> potentials acting
    type(TPotentials), intent(in) :: potential

    !> energy contributions
    type(TEnergies), intent(inout) :: energy

    !> 3rd order settings
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> Solvation model
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> Data from rangeseparated calculations
    type(TRangeSepFunc), intent(inout), allocatable :: rangeSep

    !> data type for REKS
    type(TReksCalc), allocatable, intent(inout) :: reks

    !> Proxy for querying Q-dependant external potentials
    type(TQDepExtPotProxy), intent(inout), allocatable :: qDepExtPot

    !> block (dual) atomic populations
    real(dp), intent(in), allocatable :: qBlock(:,:,:,:)

    !> Imaginary part of block atomic populations
    real(dp), intent(in), allocatable :: qiBlock(:,:,:,:)

    !> Spin orbit constants
    real(dp), intent(in), allocatable :: xi(:,:)

    !> Atoms over which to sum the total energies
    integer, intent(in) :: iAtInCentralRegion(:)

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> If tFixEf is .true. contains reservoir chemical potential, otherwise the Fermi levels found
    !> from the given number of electrons
    real(dp), intent(inout) :: Ef(:)

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    integer :: nSpin
    real(dp) :: nEl(2)

    nSpin = size(qOrb, dim=3)

    ! Tr[H0 * Rho] can be done with the same algorithm as Mulliken-analysis
    energy%atomNonSCC(:) = 0.0_dp
    call mulliken(energy%atomNonSCC, rhoPrim(:,1), H0, orb, neighbourList%iNeighbour, nNeighbourSK,&
        & img2CentCell, iSparseStart)
    energy%EnonSCC = sum(energy%atomNonSCC(iAtInCentralRegion))

    energy%atomExt(:) = 0.0_dp
    if (tExtField) then
      energy%atomExt(:) = energy%atomExt&
          & + sum(qOrb(:,:,1) - q0(:,:,1), dim=1) * potential%extAtom(:,1)
    end if
    if (allocated(qDepExtPot)) then
      call qDepExtPot%addEnergy(energy%atomExt)
    end if
    energy%Eext = sum(energy%atomExt)

    if (allocated(sccCalc)) then
      if (isXlbomd) then
        call sccCalc%getEnergyPerAtomXlbomd(species, orb, qOrb, q0, energy%atomSCC)
      else
        call sccCalc%getEnergyPerAtom(energy%atomSCC)
      end if
      energy%Escc = sum(energy%atomSCC(iAtInCentralRegion))

      if (nSpin > 1) then
        energy%atomSpin(:) = 0.5_dp * sum(sum(potential%intShell(:,:,2:nSpin)&
            & * chargePerShell(:,:,2:nSpin), dim=1), dim=2)
        energy%Espin = sum(energy%atomSpin(iAtInCentralRegion))
      end if
    end if

    if (allocated(thirdOrd)) then
      if (isXlbomd) then
        call thirdOrd%getEnergyPerAtomXlbomd(qOrb, q0, species, orb, energy%atom3rd)
      else
        call thirdOrd%getEnergyPerAtom(energy%atom3rd)
      end if
      energy%e3rd = sum(energy%atom3rd(iAtInCentralRegion))
    end if

    if (allocated(solvation)) then
      call solvation%getEnergies(energy%atomSolv)
      energy%eSolv = sum(energy%atomSolv(iAtInCentralRegion))
    end if

    if (allocated(onSiteElements)) then
      call getEons(energy%atomOnSite, qBlock, qiBlock, q0, onSiteElements, species, orb)
      energy%eOnSite = sum(energy%atomOnSite)
    end if

    if (allocated(dftbU)) then
      if (allocated(qiBlock)) then
        call dftbU%getEnergy(energy%atomDftbu, qBlock, species, orb, qiBlock)
      else
        call dftbU%getEnergy(energy%atomDftbu, qBlock, species, orb)
      end if
      energy%Edftbu = sum(energy%atomDftbu(iAtInCentralRegion))
    end if

    if (tDualSpinOrbit) then
      energy%atomLS(:) = 0.0_dp
      call getDualSpinOrbitEnergy(energy%atomLS, qiBlock, xi, orb, species)
      energy%ELS = sum(energy%atomLS(iAtInCentralRegion))
    end if

    ! Add exchange contribution for range separated calculations
    if (allocated(rangeSep) .and. .not. allocated(reks)) then
      energy%Efock = 0.0_dp
      call rangeSep%addLREnergy(energy%Efock)
    end if

    ! Free energy contribution if attached to an electron reservoir
    if (tFixEf) then
      if (nSpin == 2) then
        nEl(:) = sum(sum(qOrb(:,iAtInCentralRegion,:),dim=1),dim=1)
        nEl(1) = 0.5_dp * ( nEl(1) + nEl(2) )
        nEl(2) = nEl(1) - nEl(2)
        energy%NEf = sum(nEl(:2) * Ef(:2))
      else
        nEl = 0.0_dp
        nEl(1) = sum(qOrb(:,iAtInCentralRegion,1))
        energy%NEf = nEl(1) * Ef(1)
      end if
    else
      energy%NEf = 0.0_dp
    end if

    energy%pV = cellVol * extPressure

    ! Electronic entropy term
    energy%TS = TS

  end subroutine calcEnergies


  !> Calculates repulsive energy for current geometry
  subroutine calcRepulsiveEnergy(coord, species, img2CentCell, nNeighbourRep, neighbourList,&
      & pRepCont, Eatom, Etotal, iAtInCentralRegion)

    !> All atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> All atoms chemical species
    integer, intent(in) :: species(:)

    !> Image atom indices to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> Number of neighbours for each atom within the repulsive distance
    integer, intent(in) :: nNeighbourRep(:)

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Repulsive interaction data
    type(TRepCont), intent(in) :: pRepCont

    !> Energy for each atom
    real(dp), intent(out) :: Eatom(:)

    !> Total energy
    real(dp), intent(out) :: Etotal

    !> atoms in the central cell (or device region if transport)
    integer, intent(in) :: iAtInCentralRegion(:)

    call getERep(Eatom, coord, nNeighbourRep, neighbourList%iNeighbour, species, pRepCont,&
        & img2CentCell)
    Etotal = sum(Eatom(iAtInCentralRegion))

  end subroutine calcRepulsiveEnergy


  !> Calculates dispersion energy for current geometry.
  subroutine calcDispersionEnergy(dispersion, Eatom, Etotal, iAtInCentralRegion)

    !> dispersion interactions
    class(TDispersionIface), intent(inout) :: dispersion

    !> energy per atom
    real(dp), intent(out) :: Eatom(:)

    !> total energy
    real(dp), intent(out) :: Etotal

    !> atoms in the central cell (or device region if transport)
    integer, intent(in) :: iAtInCentralRegion(:)

    call dispersion%getEnergies(Eatom)
  #:if WITH_MBD
    select type (dispersion)
    type is (TDispMbd)
      call dispersion%checkError()
    end select
  #:endif
    Etotal = sum(Eatom(iAtInCentralRegion))

  end subroutine calcDispersionEnergy


  !> Sums together components of final energies
  subroutine sumEnergies(energy)

    !> energy contributions
    type(TEnergies), intent(inout) :: energy

    energy%Eelec = energy%EnonSCC + energy%ESCC + energy%Espin + energy%ELS + energy%Edftbu&
        & + energy%Eext + energy%e3rd + energy%eOnSite + energy%ESolv + energy%Efock

    energy%atomElec(:) = energy%atomNonSCC + energy%atomSCC + energy%atomSpin + energy%atomDftbu&
        & + energy%atomLS + energy%atomExt + energy%atom3rd + energy%atomOnSite &
        & + energy%atomSolv
    energy%atomTotal(:) = energy%atomElec + energy%atomRep + energy%atomDisp + energy%atomHalogenX
    energy%Etotal = energy%Eelec + energy%Erep + energy%eDisp + energy%eHalogenX
    energy%EMermin = energy%Etotal - sum(energy%TS)
    ! energy extrapolated to 0 K
    energy%Ezero = energy%Etotal - 0.5_dp * sum(energy%TS)
    energy%EGibbs = energy%EMermin + energy%pV

    ! Free energy of system, with contribution if attached to an electron reservoir
    ! negative sign due to electron charge
    energy%EForceRelated = energy%EGibbs  - energy%NEf

  end subroutine sumEnergies

end module dftbp_getenergies
