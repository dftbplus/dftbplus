!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module to wrap around the different energy components in the DFTB total energy expression
module energies
  use assert
  use accuracy
  implicit none
  private

  public :: TEnergies, init


  !> Data type to store components of the energy as named variables instead of
  !> in the old arrays - makes extending energy expression easier.
  type TEnergies

    !> repulsive energy
    real(dp) :: Erep    = 0.0_dp

    !> Non-SCC energy
    real(dp) :: EnonSCC = 0.0_dp

    !> SCC energy
    real(dp) :: ESCC    = 0.0_dp

    !> spin energy
    real(dp) :: Espin   = 0.0_dp

    !> spin orbit energy
    real(dp) :: ELS     = 0.0_dp

    !> DFTB+U energy
    real(dp) :: Edftbu  = 0.0_dp

    !> energy in external field
    real(dp) :: Eext    = 0.0_dp

    !> total electronic energy
    real(dp) :: Eelec   = 0.0_dp

    !> Dispersion energy
    real(dp) :: eDisp   = 0.0_dp

    !> Total 3rd order
    real(dp) :: e3rd    = 0.0_dp

    !> Excitation energy
    real(dp) :: Eexcited = 0.0_dp

    !> total energy (Erep+Etotal)
    real(dp) :: Etotal  = 0.0_dp

    !> Total Mermin energy
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

    !> atom resolved 3rd order
    real(dp), allocatable :: atom3rd(:)

    !> atom resolved total
    real(dp), allocatable :: atomTotal(:)

    !> data structure initialised
    logical :: tInitialised = .false.

  end type TEnergies


  !> initialise the data type for storing energies
  interface init
    module procedure Energies_init
  end interface init

contains


  !> Allocates storage for the energy components
  subroutine Energies_init(self, nAtom)

    !> data structure to allocate
    type(TEnergies), intent(out) :: self

    !> number of atoms needed for atom resolved arrays
    integer, intent(in) :: nAtom

    allocate(self%atomRep(nAtom))
    allocate(self%atomNonSCC(nAtom))
    allocate(self%atomSCC(nAtom))
    allocate(self%atomSpin(nAtom))
    allocate(self%atomLS(nAtom))
    allocate(self%atomDftbu(nAtom))
    allocate(self%atomExt(nAtom))
    allocate(self%atomElec(nAtom))
    allocate(self%atomDisp(nAtom))
    allocate(self%atom3rd(nAtom))
    allocate(self%atomTotal(nAtom))
    self%atomRep(:) = 0.0_dp
    self%atomNonSCC(:) = 0.0_dp
    self%atomSCC(:) = 0.0_dp
    self%atomSpin(:) = 0.0_dp
    self%atomLS(:) = 0.0_dp
    self%atomDftbu(:) = 0.0_dp
    self%atomExt(:) = 0.0_dp
    self%atomElec(:) = 0.0_dp
    self%atomDisp(:) = 0.0_dp
    self%atom3rd(:) = 0.0_dp
    self%atomTotal(:) = 0.0_dp

    self%Erep = 0.0_dp
    self%EnonSCC = 0.0_dp
    self%ESCC = 0.0_dp
    self%Espin = 0.0_dp
    self%ELS = 0.0_dp
    self%Edftbu = 0.0_dp
    self%Eext = 0.0_dp
    self%Eelec = 0.0_dp
    self%EDisp = 0.0_dp
    self%E3rd = 0.0_dp
    self%Etotal = 0.0_dp
    self%EMermin = 0.0_dp
    self%EGibbs = 0.0_dp
    self%EKin = 0.0_dp
    self%EMerminKin = 0.0_dp
    self%EGibbsKin = 0.0_dp

  end subroutine Energies_init

end module energies
