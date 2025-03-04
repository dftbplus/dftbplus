!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Evaluate energies
module dftbp_dftb_getenergies
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_densitymatrix, only : TDensityMatrix
  use dftbp_dftb_determinants, only : TDftbDeterminants, determinants
  use dftbp_dftb_dftbplusu, only : TDftbU
  use dftbp_dftb_dispiface, only : TDispersionIface
  use dftbp_dftb_energytypes, only : TEnergies
  use dftbp_dftb_hybridxc, only : THybridXcFunc
  use dftbp_dftb_onsitecorrection, only : getEons
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_populations, only : mulliken
  use dftbp_dftb_potentials, only : TPotentials
  use dftbp_dftb_rangeseponscorr, only : TRangeSepOnsCorrFunc
  use dftbp_dftb_repulsive_repulsive, only : TRepulsive
  use dftbp_dftb_scc, only : TScc
  use dftbp_dftb_spinorbit, only : getDualSpinOrbitShift, getDualSpinOrbitEnergy
  use dftbp_dftb_thirdorder, only : TThirdOrder
  use dftbp_dftbplus_qdepextpotproxy, only : TQDepExtPotProxy
  use dftbp_extlibs_tblite, only : TTBLite
  use dftbp_io_message, only : error
  use dftbp_reks_reks, only : TReksCalc
  use dftbp_solvation_solvation, only : TSolvation
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_multipole, only : TMultipole
#:if WITH_MBD
  use dftbp_dftb_dispmbd, only: TDispMbd
#:endif
  implicit none

  private
  public :: calcEnergies, calcDispersionEnergy, sumEnergies

contains


  !> Calculates various energy contribution that can potentially update for the same geometry
  subroutine calcEnergies(env, sccCalc, tblite, qOrb, q0, chargePerShell, multipole, species,&
      & isExtField, isXlbomd, dftbU, tDualSpinOrbit, rhoPrim, H0, orb, neighbourList,&
      & nNeighbourSK, img2CentCell, iSparseStart, cellVol, extPressure, TS, potential,&
      & energy, thirdOrd, solvation, hybridXc, rsOnsCorr, reks, qDepExtPot, qBlock,&
      & qiBlock, xi, iAtInCentralRegion, tFixEf, Ef, tRealHS, onSiteElements, errStatus,&
      & qNetAtom, vOnSiteAtomInt, vOnSiteAtomExt, densityMatrix, kWeights, localKS)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> SCC module internal variables
    type(TScc), allocatable, intent(in) :: sccCalc

    !> Library interface handler
    type(TTBLite), allocatable, intent(inout) :: tblite

    !> Electrons in each atomic orbital
    real(dp), intent(in) :: qOrb(:,:,:)

    !> Reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> Electrons in each atomi shell
    real(dp), intent(in) :: chargePerShell(:,:,:)

    !> Multipole moments
    type(TMultipole), intent(in) :: multipole

    !> Chemical species
    integer, intent(in) :: species(:)

    !> Is an external field present
    logical, intent(in) :: isExtField

    !> Is the extended Lagrangian being used for MD
    logical, intent(in) :: isXlbomd

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> Is dual spin orbit being used
    logical, intent(in) :: tDualSpinOrbit

    !> Density matrix in sparse storage
    real(dp), intent(in) :: rhoPRim(:,:)

    !> Non-self-consistent hamiltonian
    real(dp), intent(in) :: H0(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Neighbour list
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours within cut-off for each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> Image to real atom mapping
    integer, intent(in) :: img2CentCell(:)

    !> Index for sparse large matrices
    integer, intent(in) :: iSparseStart(:,:)

    !> Unit cell volume
    real(dp), intent(in) :: cellVol

    !> External pressure
    real(dp), intent(in) :: extPressure

    !> Electron entropy contribution
    real(dp), intent(in) :: TS(:)

    !> Potentials acting
    type(TPotentials), intent(in) :: potential

    !> Energy contributions
    type(TEnergies), intent(inout) :: energy

    !> 3rd order settings
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> Solvation model
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> Data from hybrid xc-functional calculations
    class(THybridXcFunc), intent(inout), allocatable :: hybridXc

    !> Onsite correction data with range-separated functional
    type(TRangeSepOnsCorrFunc), allocatable, intent(inout) :: rsOnsCorr

    !> Data type for REKS
    type(TReksCalc), allocatable, intent(inout) :: reks

    !> Proxy for querying Q-dependant external potentials
    type(TQDepExtPotProxy), intent(inout), allocatable :: qDepExtPot

    !> Block (dual) atomic populations
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

    !> True, if overlap and Hamiltonian are real-valued
    logical, intent(in) :: tRealHS

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !> Net atom populations
    real(dp), intent(in), optional :: qNetAtom(:)

    !> On-site only (internal) potential
    real(dp), intent(in), optional :: vOnSiteAtomInt(:,:)

    !> On-site only (external) potential
    real(dp), intent(in), optional :: vOnSiteAtomExt(:,:)

    !> Holds real and complex delta density matrices and pointers
    type(TDensityMatrix), intent(in), optional :: densityMatrix

    !> The k-point weights
    real(dp), intent(in), optional :: kWeights(:)

    !> The (K, S) tuples of the local processor group (localKS(1:2,iKS))
    !> Usage: iK = localKS(1, iKS); iS = localKS(2, iKS)
    integer, intent(in), optional :: localKS(:,:)

    integer :: nSpin
    real(dp) :: nEl(2)

    nSpin = size(qOrb, dim=3)

    ! Tr[H0 * Rho] can be done with the same algorithm as Mulliken-analysis
    energy%atomNonSCC(:) = 0.0_dp
    call mulliken(env, energy%atomNonSCC, rhoPrim(:,1), H0, orb, neighbourList%iNeighbour,&
        & nNeighbourSK, img2CentCell, iSparseStart)
    energy%EnonSCC = sum(energy%atomNonSCC(iAtInCentralRegion))

    energy%atomExt(:) = 0.0_dp
    if (isExtField) then
      energy%atomExt(:) = energy%atomExt&
          & + sum(qOrb(:,:,1) - q0(:,:,1), dim=1) * potential%extAtom(:,1)
      if (allocated(potential%extDipoleAtom) .and. allocated(multipole%dipoleAtom)) then
        energy%atomExt(:) = energy%atomExt &
            & + sum(potential%extDipoleAtom * multipole%dipoleAtom(:, :, 1), 1)
      end if
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
    end if

    if (allocated(sccCalc) .or. allocated(tblite)) then
      if (nSpin > 1) then
        energy%atomSpin(:) = 0.5_dp * sum(sum(potential%intShell(:,:,2:nSpin)&
          & * chargePerShell(:,:,2:nSpin), dim=1), dim=2)
        energy%Espin = sum(energy%atomSpin(iAtInCentralRegion))
      end if
    end if

    if (allocated(tblite)) then
      call tblite%getEnergies(energy%atomSCC)
      energy%Escc = sum(energy%atomSCC(iAtInCentralRegion))
    end if

    if (present(qNetAtom)) then
      if (present(vOnSiteAtomExt)) then
        energy%atomExt = energy%atomExt + (qNetAtom - sum(q0(:,:,1),dim=1)) * vOnSiteAtomExt(:,1)
        energy%Eext = energy%Eext + sum((qNetAtom - sum(q0(:,:,1),dim=1)) * vOnSiteAtomExt(:,1))
      end if
      if (present(vOnSiteAtomInt)) then
        energy%atomScc = energy%atomScc + (qNetAtom - sum(q0(:,:,1),dim=1)) * vOnSiteAtomInt(:,1)
        energy%EScc = energy%EScc + sum((qNetAtom - sum(q0(:,:,1),dim=1)) * vOnSiteAtomInt(:,1))
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
    if (allocated(hybridXc) .and. .not. allocated(reks)) then
      if (tRealHS) then
        call hybridXc%getHybridEnergy_real(env, energy%Efock)
      else
        if ((.not. present(densityMatrix)) .or. (.not. present(kWeights))) then
          @:RAISE_ERROR(errStatus, -1, "Missing expected array(s) for hybrid xc-functional&
              & calculation.")
        end if
        call hybridXc%getHybridEnergy_kpts(env, localKS, densityMatrix%iKiSToiGlobalKS, kWeights,&
            & densityMatrix%deltaRhoOutCplx, energy%Efock)
      end if
    end if

    ! Add long-range onsite contribution from range separated calculations
    if (allocated(rsOnsCorr)) then
      energy%EfockOnSite = 0.0_dp
      call rsOnsCorr%addLrOcEnergy(energy%EfockOnSite)
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


  !> Calculates dispersion energy for current geometry.
  subroutine calcDispersionEnergy(dispersion, Eatom, Etotal, iAtInCentralRegion)

    !> Dispersion interactions
    class(TDispersionIface), intent(inout) :: dispersion

    !> Energy per atom
    real(dp), intent(out) :: Eatom(:)

    !> Total energy
    real(dp), intent(out) :: Etotal

    !> Atoms in the central cell (or device region if transport)
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

    !> Energy contributions
    type(TEnergies), intent(inout) :: energy

    energy%Eelec = energy%EnonSCC + energy%ESCC + energy%Espin + energy%ELS + energy%Edftbu&
        & + energy%Eext + energy%e3rd + energy%eOnSite + energy%ESolv + energy%Efock&
        & + energy%EfockOnSite

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

end module dftbp_dftb_getenergies
