!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains the interface for an external population dependent potential.
module dftbp_qdepextpotgen
  use dftbp_accuracy, only : dp
  implicit none
  private

  public :: TQDepExtPotGen, TQDepExtPotGenWrapper

  !> Base class for generating external population dependent potentials.
  !>
  !> It should be extended whenever DFTB+ needs to be coupled with an external potential provider.
  !> Additional attributes and methods can be added in order to ease the conversion between the
  !> external tool and DFTB+.
  !>
  type, abstract :: TQDepExtPotGen
  contains
    !> get the external potential
    procedure(getExtPotIface), deferred :: getExternalPot
    !> get the gradient of the potential wrt DFTB atom positions
    procedure(getExtPotGradIface), deferred :: getExternalPotGrad
  end type TQDepExtPotGen


  !> Wrapper around TQDepExtPotGen to allow for building arrays.
  type :: TQDepExtPotGenWrapper

    !> TQDepExtPotGen instance to wrap.
    class(TQDepExtPotGen), allocatable :: instance

  end type TQDepExtPotGenWrapper


  abstract interface

    !> Called, when DFTB+ needs the value of the external potential at the position of the atoms.
    !>
    !> The routine is called in every SCC iteration so that the external potential can be updated
    !> according to the atom charges. The actual implementation should return both potential
    !> contributions (or zero them out, if not needed). The atom and shell-resolved contribution
    !> will be then added in DFTB+.
    !>
    !> Note: External potential is defined as the external potential the electrons feel, so in case
    !> of an electrostatic potential you would have to invert its sign.
    !>
    subroutine getExtPotIface(this, chargePerAtom, chargePerShell, extPotAtom, extPotShell)
      import :: TQDepExtPotGen, dp

      !> Instance.
      class(TQDepExtPotGen), intent(inout) :: this

      !> Net number of electrons on the atom with respect to its reference configuration, the
      !> neutral atom.  Shape: [nAtom].
      real(dp), intent(in) :: chargePerAtom(:)

      !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
      real(dp), intent(in) :: chargePerShell(:,:)

      !> Atom dependent external potential contribution. Shape: [nAtom]
      real(dp), intent(out) :: extPotAtom(:)

      !> Shell-resolved external potential contribution. Shape: [mShell, nAtom].
      real(dp), intent(out) :: extPotShell(:,:)

    end subroutine getExtPotIface


    !> Called when DFTB needs the gradient of the external potential at the position of the atoms.
    !>
    !>
    !> The routine is only called once after finishing the SCC iteration, provied the calculator has
    !> been set up to calculate forces.
    !> Note: External potential is defined as the external potential the electrons feel, so in case
    !> of an electrostatic potential you would have to invert its sign.
    !>
    subroutine getExtPotGradIface(this, chargePerAtom, chargePerShell, extPotGrad)
      import :: TQDepExtPotGen, dp

      !> Class instance.
      class(TQDepExtPotGen), intent(inout) :: this

      !> Net number of electrons on the atom with respect to its reference configuration, the
      !> neutral atom.  Shape: [nAtom].
      real(dp), intent(in) :: chargePerAtom(:)

      !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
      real(dp), intent(in) :: chargePerShell(:,:)

      !> Gradient of the potential at each atom. Shape: [3, nAtom].
      real(dp), intent(out) :: extPotGrad(:,:)

    end subroutine getExtPotGradIface

  end interface


end module dftbp_qdepextpotgen
