!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Module to wrap around the different energy components in the DFTB
!!* total energy expression
module energies
  use assert
  use accuracy
  implicit none

  public :: TEnergies, init

  private

  !!* data type to store components of the energy as named variables instead of
  !!* in the old arrays - makes extending energy expression easier.
  type TEnergies
    real(dp) :: Erep    = 0.0_dp !* repulsive energy
    real(dp) :: EnonSCC = 0.0_dp !* Non-SCC energy
    real(dp) :: ESCC    = 0.0_dp !* SCC energy
    real(dp) :: Espin   = 0.0_dp !* spin energy
    real(dp) :: ELS     = 0.0_dp !* spin orbit energy
    real(dp) :: Edftbu  = 0.0_dp !* DFTB+U energy
    real(dp) :: Eext    = 0.0_dp !* energy in external field
!    real(dp) :: EextB   = 0.0_dp !* energy in external B field
    real(dp) :: Eelec   = 0.0_dp !* total electronic energy
    real(dp) :: eDisp   = 0.0_dp !* Dispersion energy
    real(dp) :: e3rd    = 0.0_dp !* Total 3rd order
    real(dp) :: Eexcited = 0.0_dp !* Excitation energy
    real(dp) :: Etotal  = 0.0_dp !* total energy (Erep+Etotal)

    real(dp), allocatable :: atomRep(:)      !* atom resolved repulsive
    real(dp), allocatable :: atomNonSCC(:)   !* atom resolved non-SCC
    real(dp), allocatable :: atomSCC(:)      !* atom resolved SCC
    real(dp), allocatable :: atomSpin(:)     !* atom resolved spin
    real(dp), allocatable :: atomLS(:)       !* atom resolved spin orbit
    real(dp), allocatable :: atomDftbu(:)    !* atom resolved DFTB+U
    real(dp), allocatable :: atomExt(:)      !* atom resolved external field
!    real(dp), allocatable :: atomExtB(:)     !* atom resolved external B field
    real(dp), allocatable :: atomElec(:)     !* atom resolved electronic total
    real(dp), allocatable :: atomDisp(:)     !* atom resolved dispersion
    real(dp), allocatable :: atom3rd(:)      !* atom resolved 3rd order
    real(dp), allocatable :: atomTotal(:)    !* atom resolved total
    logical :: tInitialised = .false.
  end type TEnergies

  interface init
    module procedure Energies_init
  end interface init


contains

  !!* Allocates storage for the energy components
  !!* @param self data structure to allocate
  !!* @param nAtom number of atoms needed for atom resolved arrays
  subroutine Energies_init(self,nAtom)
    type(TEnergies), intent(out) :: self
    integer, intent(in) :: nAtom

    @:ASSERT(.not. self%tInitialised)

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
!    self%atomExtB(:) = 0.0_dp
    self%atomElec(:) = 0.0_dp
    self%atomDisp(:) = 0.0_dp
    self%atom3rd(:) = 0.0_dp
    self%atomTotal(:) = 0.0_dp
    self%tInitialised = .true.

    self%Erep = 0.0_dp
    self%EnonSCC = 0.0_dp
    self%ESCC = 0.0_dp
    self%Espin = 0.0_dp
    self%ELS = 0.0_dp
    self%Edftbu = 0.0_dp
    self%Eext = 0.0_dp
!    self%EextB = 0.0_dp
    self%Eelec = 0.0_dp
    self%eDisp = 0.0_dp
    self%e3rd = 0.0_dp
    self%Etotal = 0.0_dp

  end subroutine Energies_init


end module energies
