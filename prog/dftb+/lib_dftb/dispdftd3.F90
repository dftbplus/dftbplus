!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> DFT-D3 dispersion model.
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
    ! D3H5 - additional H-H repulsion
    logical :: hhrepulsion
    ! D3H5 end
  end type DispDftD3Inp


  !> Internal state of the DFT-D3 dispersion
  type, extends(DispersionIface) :: DispDftD3
    private

    !> calculator to evaluate dispersion
    type(dftd3_calc), allocatable :: calculator

    !> number of atoms
    integer :: nAtom

    !> energy
    real(dp) :: dispE

    !> force contributions
    real(dp), allocatable :: gradients(:,:)

    !> lattice vectors if periodic
    real(dp) :: latVecs(3, 3)

    !> stress tensor
    real(dp) :: stress(3, 3)

    !> atomic nuber
    integer, allocatable :: izp(:)

    !> is this periodic
    logical :: tPeriodic

    !> are  the coordinates current?
    logical :: tCoordsUpdated = .false.

    ! D3H5 - additional H-H repulsion
    logical :: tHHRepulsion
    ! D3H5 end
  contains

    !> update internal store of coordinates
    procedure :: updateCoords

    !> update internal store of lattice vectors
    procedure :: updateLatVecs

    !> return energy contribution
    procedure :: getEnergies

    !> return force contribution
    procedure :: addGradients

    !> return stress tensor contribution
    procedure :: getStress

    !> cutoff distance in real space for dispersion
    procedure :: getRCutoff

    !> add D3H5 H-H repulsion to the results
    procedure :: addHHRepulsion
  end type DispDftD3

contains


  !> Inits a DispDftD3 instance.
  subroutine DispDftD3_init(this, inp, nAtom, species0, speciesNames, latVecs)

    !> Initialised instance at return.
    type(DispDftD3), intent(out) :: this

    !> Specific input parameters for DFT-D3 Grimme.
    type(DispDftD3Inp), intent(in) :: inp

    !> Nr. of atoms in the system.
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell.
    integer, intent(in) :: species0(:)

    !> Names of species.
    character(*), intent(in) :: speciesNames(:)

    !> Lattice vectors, if the system is periodic.
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
    ! D3H5 - additional H-H repulsion
    this%tHHRepulsion = inp%hhrepulsion
    ! D3H5 end

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
  subroutine updateCoords(this, neigh, img2CentCell, coords, species0)

    !> Instance of stress data
    class(DispDftD3), intent(inout) :: this

    !> Updated neighbor list.
    type(TNeighborList), intent(in) :: neigh

    !> Updated mapping to central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Updated coordinates.
    real(dp), intent(in) :: coords(:,:)

    !> Species of the atoms in the unit cell.
    integer, intent(in) :: species0(:)

    @:ASSERT(allocated(this%calculator))

    if (this%tPeriodic) then
      ! dftd3 calculates the periodic images by itself -> only coords in central cell must be
      ! passed.
      call dftd3_pbc_dispersion(this%calculator, coords(:, 1:this%nAtom), &
          & this%izp, this%latVecs, this%dispE, this%gradients, this%stress)
    else
      call dftd3_dispersion(this%calculator, coords, this%izp, this%dispE, &
          & this%gradients)
      ! D3H5 - additional H-H repulsion
      if (this%tHHRepulsion) then
        call this%addHHRepulsion(coords)
      end if
      ! D3H5 end
    end if
    this%tCoordsUpdated = .true.

  end subroutine updateCoords


  !> Notifies the object about updated lattice vectors.
  subroutine updateLatVecs(this, latVecs)

    !> Instance of stress data
    class(DispDftD3), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    @:ASSERT(all(shape(latvecs) == shape(this%latvecs)))

    this%latVecs(:,:) = latVecs
    this%tCoordsUpdated = .false.

  end subroutine updateLatVecs


  !> Returns the atomic resolved energies due to the dispersion.
  subroutine getEnergies(this, energies)

    !> Instance of stress data
    class(DispDftD3), intent(inout) :: this

    !> Contains the atomic energy contributions on exit.
    real(dp), intent(out) :: energies(:)

    @:ASSERT(allocated(this%calculator))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(size(energies) == this%nAtom)

    ! DftD3 only delivers total energy, so we distribute it evenly over all atoms.
    energies(:) = this%dispE / real(this%nAtom, dp)

  end subroutine getEnergies


  !> Adds the atomic gradients to the provided vector.
  subroutine addGradients(this, gradients)

    !> Instance of stress data
    class(DispDftD3), intent(inout) :: this

    !> The vector to increase by the gradients.
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(allocated(this%calculator))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    gradients(:,:) = gradients + this%gradients

  end subroutine addGradients


  !> Returns the stress tensor.
  subroutine getStress(this, stress)

    !> Instance of stress data
    class(DispDftD3), intent(inout) :: this

    !> stress tensor from the dispersion
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(allocated(this%calculator))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))

    stress(:,:) = this%stress

  end subroutine getStress


  !> Estimates the real space cutoff of the dispersion interaction.
  function getRCutoff(this) result(cutoff)

    !> Instance of stress data
    class(DispDftD3), intent(inout) :: this

    !> Resulting cutoff
    real(dp) :: cutoff

    ! Since dftd3-routine uses its own real space summation routines, we do not need any real space
    ! neighbors from neighbour list -> return 0 cutoff
    cutoff = 0.0_dp

  end function getRCutoff

  ! D3H5 - additional H-H repulsion
  !> Add the H-H repulsion for D3H5
  subroutine addHHRepulsion(this, coords)
    !> Instance of D3
    class(DispDftD3), intent(inout) :: this
    !> Current coordinates
    real(dp), intent(in) :: coords(:,:)

    integer :: i, j, c
    real(dp) :: r, repE, dEdR

    ! Parameters
    real(dp) :: k = 0.30_dp * 0.001593601371772_dp
    real(dp) :: e = 14.31_dp
    real(dp) :: r0 = 2.35_dp

    @:ASSERT(all(shape(coords) == [3, this%nAtom]))

    repE = 0.0_dp
    do i = 1, this%nAtom
      if (this%izp(i) == 1) then
        do j = i + 1, this%nAtom
          if (this%izp(j) == 1) then
            ! Distance in Angstrom
            r = sqrt((coords(1,i) - coords(1,j))**2 + (coords(2,i) - coords(2,j))**2 &
                & + (coords(3,i) - coords(3,j))**2) * 0.5291772083_dp
            ! Repulsion energy in a.u.
            repE = repE + k * (1.0_dp - 1.0_dp/(1.0_dp + exp(-e*(r/r0 - 1.0_dp))))
            ! Derivative
            dEdR = -k * (1.0 / (1.0 + exp(-e*(r/r0-1.0)))**2 * e/r0 * exp(-e*(r/r0-1.0)))
          end if
        end do
      end if
    end do 
    ! Add the repulsion energy to the result
    this%dispE = this%dispE + repE
  end subroutine addHHRepulsion
  ! D3H5 end

end module dispdftd3_module
