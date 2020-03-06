!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines for the additional local potential, which should enforce charge constraints.
!>
!> Note: this also has the same functional form as 3rd order SCC contributions
module dftbp_chargeconstr
  use dftbp_assert
  use dftbp_accuracy
  implicit none
  private

  public :: TChrgConstr, init
  public :: buildShift, addShiftPerAtom, addEnergyPerAtom


  !> Constraint object
  type TChrgConstr
    private

    !> Instance initialised?
    logical :: tInit = .false.

    !> Number of atoms
    integer :: nAtom

    !> Exponent of potential
    integer :: kappa

    !> Target charges
    real(dp), allocatable :: refCharges(:)

    !> Prefactor for constraint
    real(dp), allocatable :: prefactors(:)

    !> Potential from constraint
    real(dp), allocatable :: shift(:)
  end type TChrgConstr


  !> Initialise
  interface init
    module procedure ChrgConstr_init
  end interface


  !> build the shift
  interface buildShift
    module procedure ChrgConstr_buildShift
  end interface


  !> add the shift in
  interface addShiftPerAtom
    module procedure ChrgConstr_addShiftPerAtom
  end interface


  !> energy contributions
  interface addEnergyPerAtom
    module procedure ChrgConstr_addEnergyPerAtom
  end interface

contains


  !> Initializes
  subroutine ChrgConstr_init(sf, inp, kappa)

    !> Instance of a constraint
    type(TChrgConstr), intent(inout) :: sf

    !> Array contining reference charges and prefactors (nAtom, 2)
    real(dp), intent(in) :: inp(:,:)

    !> exponent of the local potential to add
    integer, intent(in) :: kappa

    @:ASSERT(.not. sf%tInit)
    @:ASSERT(size(inp, dim=1) > 0)
    @:ASSERT(size(inp, dim=2) == 2)

    sf%nAtom = size(inp, dim=1)
    allocate(sf%refCharges(sf%nAtom))
    allocate(sf%prefactors(sf%nAtom))
    allocate(sf%shift(sf%nAtom))
    sf%refCharges = inp(:,1)
    sf%prefactors = inp(:,2)
    sf%kappa = kappa
    sf%tInit = .true.

  end subroutine ChrgConstr_init


  !> build the shift (potential)
  subroutine ChrgConstr_buildShift(sf, chargesPerAtom)

    !> Instance of a constraint
    type(TChrgConstr), intent(inout) :: sf

    !> Atomic charges
    real(dp), intent(in) :: chargesPerAtom(:)

    @:ASSERT(sf%tInit)
    @:ASSERT(size(chargesPerAtom) == size(sf%shift))

    sf%shift = real(sf%kappa, dp) * sf%prefactors * (chargesPerAtom - sf%refCharges)**(sf%kappa - 1)

  end subroutine ChrgConstr_buildShift


  !> Add the shift onto a supplied vector
  subroutine ChrgConstr_addShiftPerAtom(sf, shiftPerAtom)

    !> Instance of a constraint
    type(TChrgConstr), intent(in) :: sf

    !> Shift to append onto
    real(dp), intent(inout) :: shiftPerAtom(:)

    @:ASSERT(sf%tInit)
    @:ASSERT(size(shiftPerAtom) == sf%nAtom)

    shiftPerAtom = shiftPerAtom + sf%shift

  end subroutine ChrgConstr_addShiftPerAtom


  !> Energy associated with constrain violation
  subroutine ChrgConstr_addEnergyPerAtom(sf, energyPerAtom, chargesPerAtom)

    !> Instance of a constraint
    type(TChrgConstr), intent(in) :: sf

    !> Energy per atom from constraint
    real(dp), intent(inout) :: energyPerAtom(:)

    !> Instance of a constraint
    real(dp), intent(in) :: chargesPerAtom(:)

    @:ASSERT(sf%tInit)
    @:ASSERT(size(energyPerAtom) == sf%nAtom)
    @:ASSERT(size(energyPerAtom) == sf%nAtom)

    energyPerAtom = energyPerAtom + sf%shift * (chargesPerAtom - sf%refCharges) / real(sf%kappa, dp)

  end subroutine ChrgConstr_addEnergyPerAtom

end module dftbp_chargeconstr
