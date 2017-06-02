!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Module to wrap around the different energy components in the DFTB
!!* total energy expression
module energies
#include "allocate.h"
#include "assert.h"
  use accuracy
  implicit none

  public :: TEnergies, create, destroy

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

    real(dp), pointer :: atomRep(:)      !* atom resolved repulsive
    real(dp), pointer :: atomNonSCC(:)   !* atom resolved non-SCC
    real(dp), pointer :: atomSCC(:)      !* atom resolved SCC
    real(dp), pointer :: atomSpin(:)     !* atom resolved spin
    real(dp), pointer :: atomLS(:)       !* atom resolved spin orbit
    real(dp), pointer :: atomDftbu(:)    !* atom resolved DFTB+U
    real(dp), pointer :: atomExt(:)      !* atom resolved external field
!    real(dp), pointer :: atomExtB(:)     !* atom resolved external B field
    real(dp), pointer :: atomElec(:)     !* atom resolved electronic total
    real(dp), pointer :: atomDisp(:)     !* atom resolved dispersion
    real(dp), pointer :: atom3rd(:)      !* atom resolved 3rd order
    real(dp), pointer :: atomTotal(:)    !* atom resolved total
    logical :: tInitialised = .false.
  end type TEnergies

  interface create
    module procedure Energies_create
  end interface

  interface destroy
    module procedure Energies_destroy
  end interface


contains

  !!* Allocates storage for the energy components
  !!* @param self data structure to allocate
  !!* @param nAtom number of atoms needed for atom resolved arrays
  subroutine Energies_create(self,nAtom)
    type(TEnergies), intent(out) :: self
    integer, intent(in) :: nAtom

    ASSERT(.not. self%tInitialised)

    INITALLOCATE_PARR(self%atomRep,(nAtom))
    INITALLOCATE_PARR(self%atomNonSCC,(nAtom))
    INITALLOCATE_PARR(self%atomSCC,(nAtom))
    INITALLOCATE_PARR(self%atomSpin,(nAtom))
    INITALLOCATE_PARR(self%atomLS,(nAtom))
    INITALLOCATE_PARR(self%atomDftbu,(nAtom))
    INITALLOCATE_PARR(self%atomExt,(nAtom))
!    INITALLOCATE_PARR(self%atomExtB,(nAtom))
    INITALLOCATE_PARR(self%atomElec,(nAtom))
    INITALLOCATE_PARR(self%atomDisp, (nAtom))
    INITALLOCATE_PARR(self%atom3rd, (nAtom))
    INITALLOCATE_PARR(self%atomTotal,(nAtom))
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

  end subroutine Energies_create



  !!* De-allocates storage for the energy components
  subroutine Energies_destroy(self)
    type(TEnergies), intent(inout) :: self

    ASSERT(self%tInitialised)

    self%tInitialised = .false.
    DEALLOCATE_PARR(self%atomRep)
    DEALLOCATE_PARR(self%atomNonSCC)
    DEALLOCATE_PARR(self%atomSCC)
    DEALLOCATE_PARR(self%atomSpin)
    DEALLOCATE_PARR(self%atomDftbu)
    DEALLOCATE_PARR(self%atomExt)
!    DEALLOCATE_PARR(self%atomExtB)
    DEALLOCATE_PARR(self%atomElec)
    DEALLOCATE_PARR(self%atomDisp)
    DEALLOCATE_PARR(self%atomTotal)
    DEALLOCATE_PARR(self%atom3rd)

  end subroutine Energies_destroy


end module energies
