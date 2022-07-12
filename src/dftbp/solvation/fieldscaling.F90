!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Scale electrostatic field properties in presence of a solvent. Note, expressions assume spherical
!> solvent cavity contains the system.
module dftbp_solvation_fieldscaling
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : warning
  use dftbp_solvation_solvation, only : TSolvation
  implicit none

  private
  public :: TScaleExtEField, init_TScaleExtEField

  type TScaleExtEField

    !> Should dipoles / fields be re-scaled for external dielectric?
    logical :: isRescaled = .false.

    !> External media dielectric constant
    real(dp) :: eps_r

    !> Is inf supported by arithmetic model
    logical :: is_finite

  contains

    procedure, private :: scaledExtEField_0
    procedure, private :: scaledExtEField_1
    procedure, private :: scaledSoluteDipole_1
    procedure, private :: scaledSoluteDipole_2
    procedure, private :: scaledSoluteDipole_3
    generic :: scaledExtEField => scaledExtEField_0, scaledExtEField_1
    generic :: scaledSoluteDipole => scaledSoluteDipole_1, scaledSoluteDipole_2,&
        & scaledSoluteDipole_3

  end type TScaleExtEField

contains


  !> Initialise type
  subroutine init_TScaleExtEField(this, solvent, isRescaled)

    !> Instance
    type(TScaleExtEField), intent(out) :: this

    !> Solvent model
    class(TSolvation), allocatable, intent(in) :: solvent

    !> Is scaling of field applied
    logical, intent(in) :: isRescaled

    this%isRescaled = isRescaled
    this%eps_r = 1.0_dp
    if (allocated(solvent)) then
      this%isRescaled = this%isRescaled
      if (this%isRescaled .and. .not.solvent%isEFieldModified()) then
        call warning("Field rescaling requested for solvent model with no electrostatic effects,&
            & turning scaling off")
        this%isRescaled= .false.
      end if
    else
      this%isRescaled = .false.
    end if
    if (this%isRescaled) then
      this%eps_r = solvent%getEpsilon_r()
    end if

    this%is_finite = .false.
    if (ieee_is_finite(this%eps_r)) then
      if (this%eps_r < huge(this%eps_r)) then
        this%is_finite = .true.
      end if
    end if

  end subroutine init_TScaleExtEField


  !> Scale incident uniform electric field from outside the cavity
  pure function scaledExtEField_0(this, E0) result(E)

    !> Instance
    class(TScaleExtEField), intent(in) :: this

    !> Electric field amplitude
    real(dp), intent(in) :: E0

    !> Resulting electric field inside cavity
    real(dp) :: E

    E = E0
    if (this%isRescaled) then
      if (this%is_finite) then
        E = E * 3.0_dp * this%eps_r / (2.0_dp * this%eps_r + 1.0_dp)
      else
        E = E * 3.0_dp / 2.0_dp
      end if
    end if

  end function scaledExtEField_0


  !> Scale incident uniform electric field from outside the cavity
  pure function scaledExtEField_1(this, E0) result(E)

    !> Instance
    class(TScaleExtEField), intent(in) :: this

    !> Electric field amplitude
    real(dp), intent(in) :: E0(:)

    !> Resulting electric field inside cavity
    real(dp) :: E(size(E0))

    real(dp) :: Er

    E(:) = E0
    if (this%isRescaled) then
      if (this%is_finite) then
        E(:) = E * 3.0_dp * this%eps_r / (2.0_dp * this%eps_r + 1.0_dp)
      else
        E(:) = E * 3.0_dp / 2.0_dp
      end if
    end if

  end function scaledExtEField_1


#:for LABEL, DIMS in [('1', ':'), ('2', ':,:'), ('3', ':,:,:')]

  !> Scale ${LABEL}$ dim array of dipole moment(s) of solute. Note, these expressions are for a
  !> point dipole inside a spherical cavity.
  pure function scaledSoluteDipole_${LABEL}$(this, mu0) result (mu)

    !> Instance
    class(TScaleExtEField), intent(in) :: this

    !> Dipole moment of solute charges
    real(dp), intent(in) :: mu0(${DIMS}$)

    real(dp), allocatable :: mu(${DIMS}$)

    mu = mu0
    if (this%isRescaled) then
      if (this%is_finite) then
        mu(${DIMS}$) = 3.0_dp * this%eps_r * mu / (2.0_dp * this%eps_r + 1.0_dp)
      else
        mu(${DIMS}$) = 3.0_dp * mu / 2.0_dp
      end if
    end if

  end function scaledSoluteDipole_${LABEL}$

#:endfor

end module dftbp_solvation_fieldscaling
