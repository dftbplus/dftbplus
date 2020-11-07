!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> MBD/TS dispersion model.
module dftbp_dispmbd
  use dftbp_accuracy, only: dp, mc, lc
  use dftbp_assert
  use dftbp_commontypes, only : TOrbitals
  use dftbp_constants, only: symbolToNumber
  use dftbp_dispiface, only: TDispersionIface
  use dftbp_environment, only: TEnvironment
  use dftbp_globalenv, only: stdOut
  use dftbp_periodic, only: TNeighbourList
  use dftbp_simplealgebra, only: determinant33
  use dftbp_typegeometry, only: TGeometry
  use mbd, only: TDispMbdInp => mbd_input_t, mbd_calc_t

  implicit none

  private
  public :: TDispMbdInp, TDispMbd, TDispMbd_init

  type, extends(TDispersionIface) :: TDispMbd
    private

    !> calculator to evaluate dispersion
    type(mbd_calc_t), allocatable :: calculator

    !> number of atoms
    integer :: nAtom

    !> unit cell volume
    real(dp) :: cellVol

    !> is the correction self-consistent or post-hoc
    logical :: isPostHoc

    !> energies for atoms
    real(dp), allocatable :: energies(:)

    !> gradient contribution
    real(dp), allocatable :: gradients(:,:)

    !> stress tensor component
    real(dp) :: stress(3,3)

    !> atomic numbers
    integer, allocatable :: izp(:)

    !> are the coordinates current?
    logical :: chargesUpdated

    !> is the coordinates current?
    logical :: energyUpdated

    !> are the gradients current?
    logical :: gradientsUpdated

    !> is the stress current?
    logical :: stressUpdated

    !> caught error code from Libmbd, zero if no error
    integer :: errCode

    !> caught error origin from Libmbd
    character(len=mc) :: errOrigin

    !> caught error message from Libmbd
    character(len=lc) :: errMessage

  contains

    procedure :: updateCoords
    procedure :: updateLatVecs
    procedure :: getEnergies
    procedure :: addGradients
    procedure :: getStress
    procedure :: getRCutoff
    procedure :: updateOnsiteCharges
    procedure :: energyAvailable
    procedure :: checkError

  end type TDispMbd

contains

  !> Initialize instance of MBD calculation
  subroutine TDispMbd_init(this, inp, geom, isPostHoc)

    !> Instance
    type(TDispMbd), intent(out) :: this

    !> MBD input structure
    type(TDispMbdInp), intent(inout) :: inp

    !> geometry of the system
    type(TGeometry), intent(in) :: geom

    !> Should MBD be evaluated after the SCC cycle updates (T) or during (F)? In principle F case
    !> returns a potential contribution (feature to be added). Assumes true if absent
    logical, intent(in), optional :: isPostHoc

    integer :: iSp

    @:ASSERT(.not. allocated(this%calculator))
    allocate(this%calculator)
    inp%printer => mbdPrinter
    call this%calculator%init(inp)
    call this%calculator%get_exception(this%errCode, this%errOrigin, this%errMessage)
    this%nAtom = size(inp%coords, 2)
    this%chargesUpdated = .false.
    this%energyUpdated = .false.
    this%gradientsUpdated = .false.
    this%stressUpdated = .false.
    this%errCode = 0
    allocate(this%energies(this%nAtom))
    allocate(this%gradients(3,this%nAtom))
    if (present(isPostHoc)) then
      this%isPostHoc = isPostHoc
    else
      this%isPostHoc = .true.
    end if
    if (allocated(inp%lattice_vectors)) then
      this%cellVol = abs(determinant33(inp%lattice_vectors))
    end if
    allocate (this%izp(size(geom%speciesNames)))
    do iSp = 1, size(geom%speciesNames)
      this%izp(iSp) = symbolToNumber(geom%speciesNames(iSp))
    end do
    if (any(this%izp == 0)) then
      this%errCode = -1
      this%errMessage = 'Only standard elements are supported'
    end if

  end subroutine TDispMbd_init


  !> Update atomic coordinates
  subroutine updateCoords(this, env, neigh, img2CentCell, coords, species0, stat)

    !> Instance
    class(TDispMbd), intent(inout) :: this

    !> Computational environment
    type(TEnvironment), intent(in) :: env

    !> Neighbour list
    type(TNeighbourList), intent(in) :: neigh

    !> Mapping into atoms of central cell
    integer, intent(in) :: img2CentCell(:)

    !> Atomic coordinates (including any periodic images)
    real(dp), intent(in) :: coords(:,:)

    !> Species of atoms in the central cell
    integer, intent(in) :: species0(:)

    !> Status of operation
    integer, intent(out), optional :: stat

    @:ASSERT(allocated(this%calculator))
    call this%calculator%update_coords(coords(:,:size(species0)))

    ! mark properties as requiring re-evaluation
    this%chargesUpdated = .false.
    this%energyUpdated = .false.
    this%gradientsUpdated = .false.
    this%stressUpdated = .false.

    call this%checkError(stat)

  end subroutine updateCoords


  !> Update the lattice vectors for periodic geometries
  subroutine updateLatVecs(this, latVecs)

    !> Instance
    class(TDispMbd), intent(inout) :: this

    !> Lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    @:ASSERT(allocated(this%calculator))
    call this%calculator%update_lattice_vectors(latVecs)

    ! mark properties as requiring re-evaluation
    this%chargesUpdated = .false.
    this%energyUpdated = .false.
    this%gradientsUpdated = .false.
    this%stressUpdated = .false.

  end subroutine updateLatVecs


  !> Return MBD atomic energies (using cached values if possible)
  subroutine getEnergies(this, energies)

    !> Instance
    class(TDispMbd), intent(inout) :: this

    !> Resulting dispersion energies
    real(dp), intent(out) :: energies(:)

    real(dp) :: energy

    @:ASSERT(allocated(this%calculator))

    if (.not.this%chargesUpdated) then
      ! cannot evaluate, so returning 0. This could happen if used post-hoc with an energy return
      ! requested inside SCC loop
      energies(:) = 0.0_dp
      this%energyUpdated = .false.
      return
    end if

    if (this%energyUpdated) then
      energies(:) = this%energies
    else
      call this%calculator%evaluate_vdw_method(energy)
      call this%calculator%get_exception(this%errCode, this%errOrigin, this%errMessage)
      energies(:) = energy / this%nAtom ! replace if MBD library gives atom resolved energies
      this%energies(:) = energies
      this%energyUpdated = .true.
    end if

  end subroutine getEnergies


  !> Adds MBD forces to a gradient (using cached values if possible)
  subroutine addGradients(this, gradients)

    !> Instance
    class(TDispMbd), intent(inout) :: this

    !> Gradients to be modified
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(allocated(this%calculator))

    if (.not.this%chargesUpdated) then
      return
    end if

    if (.not. this%gradientsUpdated) then
      call this%calculator%get_gradients(this%gradients)
      this%gradientsUpdated = .true.
    end if
    gradients(:,:) = gradients + this%gradients

  end subroutine addGradients


  !> Return the MBD stress (using cached values if possible)
  subroutine getStress(this, stress)

    !> Instance
    class(TDispMbd), intent(inout) :: this

    !> MBD stress contribution
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(allocated(this%calculator))

    if (.not.this%chargesUpdated) then
      stress(:,:) = 0.0_dp
      return
    end if

    if (this%stressUpdated) then
      stress(:,:) = this%stress
    else
      call this%calculator%get_lattice_stress(stress)
      stress(:,:) = -stress/this%cellVol
      this%stress(:,:) = stress
      this%stressUpdated = .true.
    end if

  end subroutine getStress


  !> Spatial cutoff for evaluation MBD as needed by DFTB+ neighbour lists
  real(dp) function getRCutoff(this) result(cutoff)

    !> Instance
    class(TDispMbd), intent(inout) :: this

    cutoff = 0.0_dp

  end function getRCutoff


  !> Update charges in the MBD model
  subroutine updateOnsiteCharges(this, qNetAtom, orb, referenceN0, species0, tCanUseCharges)

    !> Instance
    class(TDispMbd), intent(inout) :: this

    !> Net charges
    real(dp), intent(in), allocatable :: qNetAtom(:)

    !> Atomic orbital data
    type(TOrbitals), intent(in) :: orb

    !> Reference neutral atom data
    real(dp), intent(in) :: referenceN0(:,:)

    !> Chemical species of the atoms
    integer, intent(in) :: species0(:)

    !> Are these charges from a converged SCC calculation/are suitable to evaluate MBD from
    !> (i.e. non-converged but can be used)
    logical, intent(in) :: tCanUseCharges

    real(dp), allocatable :: cpa(:), free_charges(:)
    integer :: nAtom, i_atom, i_spec

    @:ASSERT(allocated(qNetAtom))

    if (tCanUseCharges .or. .not.this%isPostHoc) then
      ! update charges as they are either converged/suitable or this correction is being used
      ! self-consistently

      nAtom = size(qNetAtom)
      allocate(free_charges(nAtom))
      do i_atom = 1, nAtom
        i_spec = species0(i_atom)
        free_charges(i_atom) = sum(referenceN0(1:orb%nShell(i_spec), i_spec))
      end do
      cpa = 1.0_dp + (qNetAtom-free_charges)/this%izp(species0)
      call this%calculator%update_vdw_params_from_ratios(cpa)

      ! dependent properties will need re-evaluation
      this%energyUpdated = .false.
      this%gradientsUpdated = .false.
      this%stressUpdated = .false.

      ! charges have been updated though, so are available for property evaluations
      this%chargesUpdated = .true.
    else
      this%chargesUpdated = .false.
    end if

  end subroutine updateOnsiteCharges


  !> Is the dispersion energy available for use in the main code after calling getEnergies
  function energyAvailable(this)

    !> data structure
    class(TDispMbd), intent(in) :: this

    !> result (dummy for most dispersion models)
    logical :: energyAvailable

    energyAvailable = this%energyUpdated

  end function energyAvailable


  !> Raises error if it was previously caught
  subroutine checkError(this, err)

    !> data structure
    class(TDispMbd), intent(in) :: this

    !> Error code return, 0 if no problems
    integer, intent(out), optional :: err

    if (this%errCode /= 0) then
      @:ERROR_HANDLING(err, this%errCode, this%errMessage)
    end if
  end subroutine checkError


  !> Printer procedure passed to Libmbd
  subroutine mbdPrinter(str)

    !> message
    character(len=*), intent(in) :: str

    write(stdOut, "(A,A)") '* Libmbd: ', str

  end subroutine mbdPrinter


end module dftbp_dispmbd
