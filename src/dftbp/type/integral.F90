!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Data types to handle overlap related integrals
module dftbp_type_integral
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TIntegral, TIntegral_init


  !> Container to store overlap related integrals
  type :: TIntegral

    !> Overlap integrals in atomic block sparse form
    real(dp), allocatable :: overlap(:)

    !> Dipole integrals with operator on ket function
    real(dp), allocatable :: dipoleKet(:, :)

    !> Dipole integrals with operator on bra function
    real(dp), allocatable :: dipoleBra(:, :)

    !> Quadrupole integrals with operator on ket function
    real(dp), allocatable :: quadrupoleKet(:, :)

    !> Quadrupole integrals with operator on bra function
    real(dp), allocatable :: quadrupoleBra(:, :)

    !> Real Hamiltonian integrals in atomic block sparse form
    real(dp), allocatable :: hamiltonian(:, :)

    !> Imaginary Hamiltonian integrals in atomic block sparse form
    real(dp), allocatable :: iHamiltonian(:, :)

  end type TIntegral


contains


  !> Initializier for integral container
  subroutine TIntegral_init(this, nSpin, tReHam, tImHam, nDipole, nQuadrupole)

    !> Instance of the integral container
    type(TIntegral), intent(out) :: this

    !> Number of spins channels in the system
    integer, intent(in) :: nSpin

    !> Allocate space for real Hamiltonian
    logical, intent(in) :: tReHam

    !> Allocate space for imaginary Hamiltonian
    logical, intent(in) :: tImHam

    !> Number of dipole moment components included
    integer, intent(in) :: nDipole

    !> Number of quadrupole moment components included
    integer, intent(in) :: nQuadrupole

    if (tReHam) then
      allocate(this%hamiltonian(0, nSpin))
    end if
    if (tImHam) then
      allocate(this%iHamiltonian(0, nSpin))
    end if
    allocate(this%overlap(0))
    if (nDipole > 0) then
      allocate(this%dipoleKet(nDipole, 0))
      allocate(this%dipoleBra(nDipole, 0))
    end if
    if (nQuadrupole > 0) then
      allocate(this%quadrupoleKet(nQuadrupole, 0))
      allocate(this%quadrupoleBra(nQuadrupole, 0))
    end if

  end subroutine TIntegral_init


end module dftbp_type_integral
