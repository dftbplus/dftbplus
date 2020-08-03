!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> MBD/TS dispersion model.
module dftbp_dispmbd
  use dftbp_accuracy, only: dp, mc
  use dftbp_assert
  use dftbp_dispiface, only: TDispersionIface
  use dftbp_environment, only: TEnvironment
  use dftbp_periodic, only: TNeighbourList
  use dftbp_simplealgebra, only: determinant33
  use mbd, only: TDispMbdInp => mbd_input_t, mbd_calc_t

  implicit none

  private
  public :: TDispMbdInp, TDispMbd

  type, extends(TDispersionIface) :: TDispMbd
    private
    type(mbd_calc_t), allocatable :: calculator
    integer :: nAtom
    real(dp) :: cellVol

  contains

    procedure :: init
    procedure :: updateCoords
    procedure :: updateLatVecs
    procedure :: getEnergies
    procedure :: addGradients
    procedure :: getStress
    procedure :: getRCutoff
    procedure :: updateOnsiteCharges
  end type TDispMbd

contains

  subroutine init(this, inp)
    class(TDispMbd), intent(out) :: this
    type(TDispMbdInp), intent(in) :: inp

    @:ASSERT(.not. allocated(this%calculator))
    allocate(this%calculator)
    call this%calculator%init(inp)
    this%nAtom = size(inp%coords, 2)
    if (allocated(inp%lattice_vectors)) then
      this%cellVol = abs(determinant33(inp%lattice_vectors))
    end if
  end subroutine init

  subroutine updateCoords(this, env, neigh, img2CentCell, coords, species0)
    class(TDispMbd), intent(inout) :: this
    type(TEnvironment), intent(in) :: env
    type(TNeighbourList), intent(in) :: neigh
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species0(:)

    @:ASSERT(allocated(this%calculator))
    call this%calculator%update_coords(coords(:,:size(species0)))
  end subroutine updateCoords

  subroutine updateLatVecs(this, latVecs)
    class(TDispMbd), intent(inout) :: this
    real(dp), intent(in) :: latVecs(:,:)

    @:ASSERT(allocated(this%calculator))
    call this%calculator%update_lattice_vectors(latVecs)
  end subroutine updateLatVecs

  subroutine getEnergies(this, energies)
    class(TDispMbd), intent(inout) :: this
    real(dp), intent(out) :: energies(:)

    real(dp) :: energy

    @:ASSERT(allocated(this%calculator))
    call this%calculator%evaluate_vdw_method(energy)
    energies(:) = energy / this%nAtom
  end subroutine getEnergies

  subroutine addGradients(this, gradients)
    class(TDispMbd), intent(inout) :: this
    real(dp), intent(inout) :: gradients(:,:)

    real(dp), allocatable :: gradients_mbd(:, :)

    @:ASSERT(allocated(this%calculator))
    allocate (gradients_mbd(3, size(gradients, 2)))
    call this%calculator%get_gradients(gradients_mbd)
    gradients = gradients + gradients_mbd
  end subroutine addGradients

  subroutine getStress(this, stress)
    class(TDispMbd), intent(inout) :: this
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(allocated(this%calculator))
    call this%calculator%get_lattice_stress(stress)
    stress = -stress/this%cellVol
  end subroutine getStress

  real(dp) function getRCutoff(this) result(cutoff)
    class(TDispMbd), intent(inout) :: this

    cutoff = 0.0_dp
  end function getRCutoff

  subroutine updateOnsiteCharges(this, qOnsite, orb, referenceN0, speciesName, species0)
    use dftbp_commontypes, only : TOrbitals
    use mbd_vdw_param, only: species_index

    class(TDispMbd), intent(inout) :: this
    real(dp), intent(in) :: qOnsite(:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: referenceN0(:,:)
    character(mc), intent(in) :: speciesName(:)
    integer, intent(in) :: species0(:)

    real(dp), allocatable :: cpa(:), free_charges(:)
    integer :: nAtom, i_atom, i_spec

    nAtom = size(qOnsite)
    allocate(free_charges(nAtom))
    do i_atom = 1, nAtom
        i_spec = species0(i_atom)
        free_charges(i_atom) = sum(referenceN0(1:orb%nShell(i_spec), i_spec))
    end do
    cpa = 1d0 + (qOnsite-free_charges)/dble(species_index(speciesName(species0)))
    call this%calculator%update_vdw_params_from_ratios(cpa)
  end subroutine updateOnsiteCharges
end module dftbp_dispmbd
