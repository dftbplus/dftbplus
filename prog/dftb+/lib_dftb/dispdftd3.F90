!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> DFT-D3 dispersion model.
!!
module dispdftd3_module
  use assert
  use accuracy
  use dispiface
  use dftd3_module
  implicit none
  private

  public :: DispDftD3Inp, DispDftD3, DispDftD3_init


  !> Input structure for the DFT-D3 initialization.
  type :: DispDftD3Inp
    real(dp) :: s6, s8, a1, a2, sr6, sr8, alpha6
    logical :: tBeckeJohnson
    integer :: version
    real(dp) :: cutoff, cutoffCN
    logical :: threebody, numgrad
  end type DispDftD3Inp


  !> Internal state of the DFT-D3
  type, extends(DispersionIface) :: DispDftD3
    private
    type(dftd3_calc), allocatable :: calculator
    integer :: nAtom
    real(dp) :: dispE
    real(dp), allocatable :: gradients(:,:)
    real(dp) :: stress(3, 3), latVecs(3, 3)
    integer, allocatable :: izp(:)
    logical :: tPeriodic
    logical :: tCoordsUpdated = .false.
  contains
    procedure :: updateCoords
    procedure :: updateLatVecs
    procedure :: getEnergies
    procedure :: addGradients
    procedure :: getStress
    procedure :: getRCutoff
  end type DispDftD3


contains

  !> Inits a DispDftD3 instance.
  !!
  !! \param this  Initialised instance at return.
  !! \param inp  Specific input parameters for DFT-D3 Grimme.
  !! \param nAtom  Nr. of atoms in the system.
  !! \param species0  Species of every atom in the unit cell.
  !! \param speciesNames  Names of species.
  !! \param latVecs  Lattice vectors, if the system is periodic.
  !!
  subroutine DispDftD3_init(this, inp, nAtom, species0, speciesNames, latVecs)
    type(DispDftD3), intent(out) :: this
    type(DispDftD3Inp), intent(in) :: inp
    integer, intent(in) :: nAtom
    integer, intent(in) :: species0(:)
    character(*), intent(in) :: speciesNames(:)
    real(dp), intent(in), optional :: latVecs(:,:)

    type(dftd3_input) :: d3inp
    integer :: iAt

    @:ASSERT(.not. allocated(this%calculator))

    this%tPeriodic = present(latVecs)
    if (this%tPeriodic) then
      this%latVecs(:,:) = latVecs
    end if
    this%nAtom = nAtom

    d3inp%threebody = inp%threebody
    d3inp%numgrad = inp%numgrad
    d3inp%cutoff = inp%cutoff
    d3inp%cutoff_cn = inp%cutoffCN

    allocate(this%calculator)
    call dftd3_init(this%calculator, d3inp)
    if (inp%tBeckeJohnson) then
      call dftd3_set_params(this%calculator, [inp%s6, inp%a1, inp%s8, &
          & inp%a2, 0.0_dp], 4)
    else
      call dftd3_set_params(this%calculator, [inp%s6, inp%sr6, inp%s8, &
          & inp%sr8, inp%alpha6], 3)
    end if

    allocate(this%izp(nAtom))
    do iAt = 1, this%nAtom
      this%izp(iAt) =  get_atomic_number(speciesNames(species0(iAt)))
    end do
    allocate(this%gradients(3, nAtom))

  end subroutine DispDftD3_init


  !> Notifies the objects about changed coordinates.
  !!
  !! \param neigh  Updated neighbor list.
  !! \param img2CentCell  Updated mapping to central cell.
  !! \param coord  Updated coordinates.
  !! \param species0  Species of the atoms in the unit cell.
  !!
  subroutine updateCoords(this, neigh, img2CentCell, coords, species0)
    class(DispDftD3), intent(inout) :: this
    type(TNeighborList), intent(in) :: neigh
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) ::  species0(:)

    @:ASSERT(allocated(this%calculator))

    if (this%tPeriodic) then
      ! dftd3 calculates the periodic images by itself -> only coords in
      ! central cell must be passed.
      call dftd3_pbc_dispersion(this%calculator, coords(:, 1:this%nAtom), &
          & this%izp, this%latVecs, this%dispE, this%gradients, this%stress)
    else
      call dftd3_dispersion(this%calculator, coords, this%izp, this%dispE, &
          & this%gradients)
    end if
    this%tCoordsUpdated = .true.

  end subroutine updateCoords


  !> Notifies the object about updated lattice vectors.
  !!
  !! \param latVecs  New lattice vectors
  !!
  subroutine updateLatVecs(this, latVecs)
    class(DispDftD3), intent(inout) :: this
    real(dp), intent(in) :: latVecs(:,:)

    @:ASSERT(all(shape(latvecs) == shape(this%latvecs)))

    this%latVecs(:,:) = latVecs
    this%tCoordsUpdated = .false.

  end subroutine updateLatVecs


  !> Returns the atomic resolved energies due to the dispersion.
  !!
  !! \param energies  Contains the atomic energy contributions on exit.
  !!
  subroutine getEnergies(this, energies)
    class(DispDftD3), intent(inout) :: this
    real(dp), intent(out) :: energies(:)

    @:ASSERT(allocated(this%calculator))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(size(energies) == this%nAtom)

    ! DftD3 only delivers total energy, so we distribute it evenly over all
    ! atoms.
    energies(:) = this%dispE / real(this%nAtom, dp)

  end subroutine getEnergies


  !> Adds the atomic gradients to the provided vector.
  !!
  !! \param gradients  The vector to increase by the gradients.
  !!
  subroutine addGradients(this, gradients)
    class(DispDftD3), intent(inout) :: this
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(allocated(this%calculator))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    gradients(:,:) = gradients + this%gradients

  end subroutine addGradients


  !> Returns the stress tensor.
  !!
  !! \param stress tensor from the dispersion
  !!
  subroutine getStress(this, stress)
    class(DispDftD3), intent(inout) :: this
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(allocated(this%calculator))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))

    stress(:,:) = this%stress

  end subroutine getStress


  !> Estimates the real space cutoff of the dispersion interaction.
  !!
  !! \return Cutoff
  !!
  function getRCutoff(this) result(cutoff)
    class(DispDftD3), intent(inout) :: this
    real(dp) :: cutoff

    ! Since dftd3-routine uses its own real space summation routines, we do not
    ! need any real space neighbors from neighbour list -> return 0 cutoff
    cutoff = 0.0_dp

  end function getRCutoff


end module dispdftd3_module
