!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Simplified C-interface with callbacks for population dependant external potential generators.
module dftbp_qdepextpotgenc
  use, intrinsic :: iso_c_binding
  use dftbp_accuracy, only : dp
  use dftbp_qdepextpotgen, only : TQDepExtPotGen
  implicit none
  private

  public :: getExtPotIfaceC, getExtPotGradIfaceC
  public :: TQDepExtPotGenC, TQDepExtPotGenC_init

  !> Interface to the routine which calculates the external potential due to charges
  abstract interface

    !> Interface to set up external potential
    subroutine getExtPotIfaceC(refPtr, dQAtom, extPotAtom) bind(C)
      import :: c_double, c_ptr

      !> Reference pointer
      type(c_ptr), value, intent(in) :: refPtr

      !> Net number of electrons on each atom (note: positive number means electron excess)
      real(c_double), intent(in)  :: dQAtom(*)

      !> Potential on each atom (note: positive number means electron repulsion)
      real(c_double), intent(out) :: extPotAtom(*)

    end subroutine getExtPotIfaceC


    !> Interface to set up gradient of external potential
    subroutine getExtPotGradIfaceC(refPtr, dQAtom, extPotAtomGrad) bind(C)
      import :: c_double, c_ptr

      !> Reference pointer
      type(c_ptr), value, intent(in) :: refPtr

      !> Net number of electrons on each atom (note: positive number means electron excess)
      real(c_double), intent(in) :: dQAtom(*)

      !> Gradient of the potential on each atom (note: positive number means electron repulsion)
      real(c_double), intent(out) :: extPotAtomGrad(3, *)

    end subroutine getExtPotGradIfaceC

  end interface


  !> builds on the charge dependent external interface type from dftbp_qdepextpotgen for the C API
  type, extends(TQDepExtPotGen) :: TQDepExtPotGenC
    private
    type(c_ptr) :: refPtr
    procedure(getExtPotIfaceC), nopass, pointer :: getExtPot
    procedure(getExtPotGradIfaceC), nopass, pointer :: getExtPotGrad
  contains
    procedure :: getExternalPot => TDepExtPotGenC_getExternalPot
    procedure :: getExternalPotGrad => TQDepExtPotGenC_getExternalPotGrad
  end type TQDepExtPotGenC


contains


  !> Initialise an external charge dependent external potential within this type
  subroutine TQDepExtPotGenC_init(this, refPtr, extPotFunc, extPotGradFunc)

    !> Instance
    type(TQDepExtPotGenC), intent(out) :: this

    !> pointer to the C routine for the external potential
    type(c_ptr), intent(in) :: refPtr

    !> function for the potential evaluation
    procedure(getExtPotIfaceC), pointer, intent(in) :: extPotFunc

    !> function for the gradient of the potential
    procedure(getExtPotGradIfaceC), pointer, intent(in) :: extPotGradFunc

    this%getExtPot => extPotFunc
    this%getExtPotGrad => extPotGradFunc
    this%refPtr = refPtr

  end subroutine TQDepExtPotGenC_init



  !> extra routine a charge dependent external potential
  subroutine TDepExtPotGenC_getExternalPot(this, chargePerAtom, chargePerShell, extPotAtom,&
      & extPotShell)

    !> Instance.
    class(TQDepExtPotGenC), intent(inout) :: this

    !> Net number of electrons on the atom with respect to its reference configuration, the
    !> neutral atom.  Shape: [nAtom].
    real(dp), intent(in) :: chargePerAtom(:)

    !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
    real(dp), intent(in) :: chargePerShell(:,:)

    !> Atom dependent external potential contribution. Shape: [nAtom]
    real(dp), intent(out) :: extPotAtom(:)

    !> Shell-resolved external potential contribution. Shape: [mShell, nAtom].
    real(dp), intent(out) :: extPotShell(:,:)

    call this%getExtPot(this%refPtr, chargePerAtom, extPotAtom)
    ! currently only atom resolved, no non-local l-dependent part
    extPotShell(:,:) = 0.0_dp

  end subroutine TDepExtPotGenC_getExternalPot



  !> extra routine for interfacing gradients from a charge dependent external potential
  subroutine TQDepExtPotGenC_getExternalPotGrad(this, chargePerAtom, chargePerShell, extPotGrad)

    !> Class instance.
    class(TQDepExtPotGenC), intent(inout) :: this

    !> Net number of electrons on the atom with respect to its reference configuration, the
    !> neutral atom.  Shape: [nAtom].
    real(dp), intent(in) :: chargePerAtom(:)

    !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
    real(dp), intent(in) :: chargePerShell(:,:)

    !> Gradient of the potential at each atom. Shape: [3, nAtom].
    real(dp), intent(out) :: extPotGrad(:,:)

    call this%getExtPotGrad(this%refPtr, chargePerAtom, extPotGrad)

  end subroutine TQDepExtPotGenC_getExternalPotGrad


end module dftbp_qdepextpotgenc
