!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines for the additional local potential, which should enforce charge constraints.
!>
!> Note: this also has the same functional form as onsite 3rd order SCC contributions
module dftbp_dftb_chargeconstr
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TChrgConstr, TChrgConstr_init


  !> Constraint object
  type TChrgConstr
    private

    ! Instance initialised?
    logical :: tInit_ = .false.

    ! Number of atoms
    integer :: nAtom_

    ! Exponent of potential
    integer :: kappa_

    ! Target charges
    real(dp), allocatable :: refCharges_(:)

    ! Prefactor for constraint
    real(dp), allocatable :: prefactors_(:)

    ! Potential from constraint
    real(dp), allocatable :: shift_(:)

  contains

    procedure :: buildShift
    procedure :: addShiftPerAtom
    procedure :: addEnergyPerAtom

  end type TChrgConstr


contains


  !> Initializes
  subroutine TChrgConstr_init(this, inp, kappa)

    !> Instance of a constraint
    type(TChrgConstr), intent(out) :: this

    !> Array containing reference charges and prefactors (nAtom, 2)
    real(dp), intent(in) :: inp(:,:)

    !> exponent of the local potential to add
    integer, intent(in) :: kappa

    @:ASSERT(.not. this%tInit_)
    @:ASSERT(size(inp, dim=1) > 0)
    @:ASSERT(size(inp, dim=2) == 2)

    this%nAtom_ = size(inp, dim=1)
    allocate(this%refCharges_(this%nAtom_))
    allocate(this%prefactors_(this%nAtom_))
    allocate(this%shift_(this%nAtom_))
    this%refCharges_ = inp(:,1)
    this%prefactors_ = inp(:,2)
    this%kappa_ = kappa
    this%tInit_ = .true.

  end subroutine TChrgConstr_init


  !> build the shift (potential)
  subroutine buildShift(this, chargesPerAtom)

    !> Instance of a constraint
    class(TChrgConstr), intent(inout) :: this

    !> Atomic charges
    real(dp), intent(in) :: chargesPerAtom(:)

    @:ASSERT(this%tInit_)
    @:ASSERT(size(chargesPerAtom) == size(this%shift_))

    this%shift_(:) = real(this%kappa_, dp) * this%prefactors_&
        & * (chargesPerAtom - this%refCharges_)**(this%kappa_ - 1)

  end subroutine buildShift


  !> Add the shift onto a supplied vector
  subroutine addShiftPerAtom(this, shiftPerAtom)

    !> Instance of a constraint
    class(TChrgConstr), intent(in) :: this

    !> Shift to append onto
    real(dp), intent(inout) :: shiftPerAtom(:)

    @:ASSERT(this%tInit_)
    @:ASSERT(size(shiftPerAtom) == this%nAtom_)

    shiftPerAtom(:) = shiftPerAtom + this%shift_

  end subroutine addShiftPerAtom


  !> Energy associated with constrain violation
  subroutine addEnergyPerAtom(this, energyPerAtom, chargesPerAtom)

    !> Instance of a constraint
    class(TChrgConstr), intent(in) :: this

    !> Energy per atom from constraint
    real(dp), intent(inout) :: energyPerAtom(:)

    !> Instance of a constraint
    real(dp), intent(in) :: chargesPerAtom(:)

    @:ASSERT(this%tInit_)
    @:ASSERT(size(energyPerAtom) == this%nAtom_)
    @:ASSERT(size(energyPerAtom) == this%nAtom_)

    energyPerAtom(:) = energyPerAtom&
        & + this%shift_ * (chargesPerAtom - this%refCharges_) / real(this%kappa_, dp)

  end subroutine addEnergyPerAtom

end module dftbp_dftb_chargeconstr
