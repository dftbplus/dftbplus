!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module to wrap around the different energy components in the DFTB total energy expression
module dftbp_energytypes
  use dftbp_assert
  use dftbp_accuracy
  implicit none
  private

  public :: TEnergies, TEnergies_init


  !> Data type to store components of the energy as named variables instead of
  !> in the old arrays - makes extending energy expression easier.
  type TEnergies

    !> repulsive energy
    real(dp) :: Erep = 0.0_dp

    !> Non-SCC energy
    real(dp) :: EnonSCC = 0.0_dp

    !> SCC energy
    real(dp) :: ESCC = 0.0_dp

    !> spin energy
    real(dp) :: Espin = 0.0_dp

    !> range-separation energy
    real(dp) :: Efock = 0.0_dp

    !> spin orbit energy
    real(dp) :: ELS = 0.0_dp

    !> DFTB+U energy
    real(dp) :: Edftbu = 0.0_dp

    !> energy in external field
    real(dp) :: Eext = 0.0_dp

    !> total electronic energy
    real(dp) :: Eelec = 0.0_dp

    !> Dispersion energy
    real(dp) :: eDisp = 0.0_dp

    !> Number of electrons times Fermi energy for cases of the system being connected to reservoir
    real(dp) :: NEf = 0.0_dp

    !> Presure times volume for periodic systems
    real(dp) :: pV = 0.0_dp

    !> Electronic entropy times temperature
    real(dp), allocatable :: TS(:)

    !> band structure energy
    real(dp), allocatable :: EBand(:)

    !> zero temperature estimated band energy
    real(dp), allocatable :: E0(:)

    !> Onsite correction energy
    real(dp) :: eOnSite = 0.0_dp

    !> Halogen-X correction energy
    real(dp) :: eHalogenX = 0.0_dp

    !> Total 3rd order
    real(dp) :: e3rd = 0.0_dp

    !> Solvation free energy
    real(dp) :: ESolv = 0.0_dp

    !> Excitation energy
    real(dp) :: Eexcited = 0.0_dp

    !> total energy (Erep+Etotal)
    real(dp) :: Etotal = 0.0_dp

    !> Total Mermin energy (note that this may be evaluated even when the TS term cannot be set, so
    !> contains the same as Etotal in those cases)
    real(dp) :: EMermin = 0.0_dp

    !> Zero temperature extrapolated energy
    real(dp) :: Ezero = 0.0_dp

    !> Gibbs free energy
    real(dp) :: EGibbs = 0.0_dp

    !> Kinetic energy
    real(dp) :: EKin = 0.0_dp

    !> Total Mermin energy including kinetic energy
    real(dp) :: EMerminKin = 0.0_dp

    !> Gibbs free energy including kinetic energy
    real(dp) :: EGibbsKin = 0.0_dp

    !> Energy or free energy which is related to the forces via the Helmann-Feynman theorem. This is
    !> used for example in geometry optimisation or energetic comparisions.
    real(dp) :: EForceRelated = 0.0_dp

    !> atom resolved repulsive
    real(dp), allocatable :: atomRep(:)

    !> atom resolved non-SCC
    real(dp), allocatable :: atomNonSCC(:)

    !> atom resolved SCC
    real(dp), allocatable :: atomSCC(:)

    !> atom resolved spin
    real(dp), allocatable :: atomSpin(:)

    !> atom resolved spin orbit
    real(dp), allocatable :: atomLS(:)

    !> atom resolved DFTB+U
    real(dp), allocatable :: atomDftbu(:)

    !> atom resolved external field
    real(dp), allocatable :: atomExt(:)

    !> atom resolved electronic total
    real(dp), allocatable :: atomElec(:)

    !> atom resolved dispersion
    real(dp), allocatable :: atomDisp(:)

    !> atom onsite correction energies
    real(dp), allocatable :: atomOnSite(:)

    !> atom halogen bond correction energies
    real(dp), allocatable :: atomHalogenX(:)

    !> atom resolved 3rd order
    real(dp), allocatable :: atom3rd(:)

    !> atom resolved solvation free energy
    real(dp), allocatable :: atomSolv(:)

    !> atom resolved total
    real(dp), allocatable :: atomTotal(:)

    !> data structure initialised
    logical :: tInitialised = .false.

  end type TEnergies


contains


  !> Allocates storage for the energy components
  subroutine TEnergies_init(this, nAtom, nSpin)

    !> data structure to allocate
    type(TEnergies), intent(out) :: this

    !> number of atoms needed for atom resolved arrays
    integer, intent(in) :: nAtom

    !> Number of independent spins
    integer, intent(in) :: nSpin

    allocate(this%TS(nSpin))
    allocate(this%E0(nSpin))
    allocate(this%EBand(nSpin))

    this%TS(:) = 0.0_dp
    this%E0(:) = 0.0_dp
    this%EBand(:) = 0.0_dp

    allocate(this%atomRep(nAtom))
    allocate(this%atomNonSCC(nAtom))
    allocate(this%atomSCC(nAtom))
    allocate(this%atomSpin(nAtom))
    allocate(this%atomLS(nAtom))
    allocate(this%atomDftbu(nAtom))
    allocate(this%atomExt(nAtom))
    allocate(this%atomElec(nAtom))
    allocate(this%atomDisp(nAtom))
    allocate(this%atomOnSite(nAtom))
    allocate(this%atomHalogenX(nAtom))
    allocate(this%atom3rd(nAtom))
    allocate(this%atomSolv(nAtom))
    allocate(this%atomTotal(nAtom))
    this%atomRep(:) = 0.0_dp
    this%atomNonSCC(:) = 0.0_dp
    this%atomSCC(:) = 0.0_dp
    this%atomSpin(:) = 0.0_dp
    this%atomLS(:) = 0.0_dp
    this%atomDftbu(:) = 0.0_dp
    this%atomExt(:) = 0.0_dp
    this%atomElec(:) = 0.0_dp
    this%atomDisp(:) = 0.0_dp
    this%atomOnSite(:) = 0.0_dp
    this%atomHalogenX(:) = 0.0_dp
    this%atom3rd(:) = 0.0_dp
    this%atomSolv(:) = 0.0_dp
    this%atomTotal(:) = 0.0_dp

    this%Erep = 0.0_dp
    this%EnonSCC = 0.0_dp
    this%ESCC = 0.0_dp
    this%Espin = 0.0_dp
    this%Efock = 0.0_dp
    this%ELS = 0.0_dp
    this%Edftbu = 0.0_dp
    this%Eext = 0.0_dp
    this%Eelec = 0.0_dp
    this%EDisp = 0.0_dp
    this%EOnSite = 0.0_dp
    this%EHalogenX = 0.0_dp
    this%E3rd = 0.0_dp
    this%ESolv = 0.0_dp
    this%Etotal = 0.0_dp
    this%EMermin = 0.0_dp
    this%EGibbs = 0.0_dp
    this%EKin = 0.0_dp
    this%EMerminKin = 0.0_dp
    this%EGibbsKin = 0.0_dp

  end subroutine TEnergies_init

end module dftbp_energytypes
