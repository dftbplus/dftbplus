!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Data type to handle multipole moments
module dftbp_type_multipole
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TMultipole, TMultipole_init


  !> Container to store multipole moments
  type :: TMultipole

    !> Cumulative atomic dipole moments
    real(dp), allocatable :: dipoleAtom(:, :, :)

    !> Cumulative atomic quadrupole moments
    real(dp), allocatable :: quadrupoleAtom(:, :, :)

  end type TMultipole


contains


  !> Create new container for multipole moments
  subroutine TMultipole_init(this, nAtom, nDipole, nQuadrupole, nSpin)

    !> Instance of the multipole moments
    type(TMultipole), intent(out) :: this

    !> Number of atoms in the central cell
    integer, intent(in) :: nAtom

    !> Number of dipole moment components
    integer, intent(in) :: nDipole

    !> Number of quadrupole moment components
    integer, intent(in) :: nQuadrupole

    !> Number of spin channels
    integer, intent(in) :: nSpin

    if (nDipole > 0) then
      allocate(this%dipoleAtom(nDipole, nAtom, nSpin), source=0.0_dp)
    end if

    if (nQuadrupole > 0) then
      allocate(this%quadrupoleAtom(nQuadrupole, nAtom, nSpin), source=0.0_dp)
    end if

  end subroutine TMultipole_init


end module dftbp_type_multipole
