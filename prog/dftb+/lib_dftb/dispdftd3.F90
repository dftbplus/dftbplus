!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> DFT-D3 dispersion model.
module dftbp_dispdftd3
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_dispiface
  use dftbp_dftd3
  use dftbp_environment, only : TEnvironment
  use dftbp_periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use dftbp_simplealgebra, only : determinant33
  use dftbp_constants
  implicit none
  private

  public :: TDispDftD3Inp, TDispDftD3, DispDftD3_init
  public :: hhRepCutOff

  ! cut-off distance in Bohr for the H-H repulsion term dropping below 1E-10
  real(dp), parameter :: hhRepCutOff = 10.0_dp

  !> Input structure for the DFT-D3 initialization.
  type :: TDispDftD3Inp
    real(dp) :: s6, s8, a1, a2, sr6, sr8, alpha6
    logical :: tBeckeJohnson
    integer :: version
    real(dp) :: cutoff, cutoffCN
    logical :: threebody, numgrad

    !> D3H5 - additional H-H repulsion
    logical :: hhrepulsion

  end type TDispDftD3Inp


  !> Internal state of the DFT-D3 dispersion
  type, extends(TDispersionIface) :: TDispDftD3
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

    !> are the coordinates current?
    logical :: tCoordsUpdated = .false.

    !> D3H5 - additional H-H repulsion
    logical :: tHHRepulsion

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

  end type TDispDftD3

contains


  !> Inits a DispDftD3 instance.
  subroutine DispDftD3_init(this, inp, nAtom, species0, speciesNames, latVecs)

    !> Initialised instance at return.
    type(TDispDftD3), intent(out) :: this

    !> Specific input parameters for DFT-D3 Grimme.
    type(TDispDftD3Inp), intent(in) :: inp

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

    this%tHHRepulsion = inp%hhrepulsion

    allocate(this%calculator)
    call dftd3_init(this%calculator, d3inp)
    if (inp%tBeckeJohnson) then
      call dftd3_set_params(this%calculator, [inp%s6, inp%a1, inp%s8, inp%a2, 0.0_dp], 4)
    else
      call dftd3_set_params(this%calculator, [inp%s6, inp%sr6, inp%s8, inp%sr8, inp%alpha6], 3)
    end if

    allocate(this%izp(nAtom))
    do iAt = 1, this%nAtom
      this%izp(iAt) =  get_atomic_number(speciesNames(species0(iAt)))
    end do
    allocate(this%gradients(3, nAtom))

  end subroutine DispDftD3_init


  !> Notifies the objects about changed coordinates.
  subroutine updateCoords(this, env, neigh, img2CentCell, coords, species0)

    !> Instance of DFTD3 data
    class(TDispDftD3), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Updated neighbour list.
    type(TNeighbourList), intent(in) :: neigh

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
      call dftd3_pbc_dispersion(this%calculator, coords(:, 1:this%nAtom), this%izp, this%latVecs,&
          & this%dispE, this%gradients, this%stress)
    else
      call dftd3_dispersion(this%calculator, coords, this%izp, this%dispE, this%gradients)
    end if

    if (this%tHHRepulsion) then
      ! D3H5 - additional H-H repulsion contributions to energy, forces and stress (if periodic)
      call this%addHHRepulsion(coords, neigh, img2CentCell)
    end if

    this%tCoordsUpdated = .true.

  end subroutine updateCoords


  !> Notifies the object about updated lattice vectors.
  subroutine updateLatVecs(this, latVecs)

    !> Instance of DFTD3 data
    class(TDispDftD3), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    @:ASSERT(all(shape(latvecs) == shape(this%latvecs)))

    this%latVecs(:,:) = latVecs
    this%tCoordsUpdated = .false.

  end subroutine updateLatVecs


  !> Returns the atomic resolved energies due to the dispersion.
  subroutine getEnergies(this, energies)

    !> Instance of DFTD3 data
    class(TDispDftD3), intent(inout) :: this

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

    !> Instance of DFTD3 data
    class(TDispDftD3), intent(inout) :: this

    !> The vector to increase by the gradients.
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(allocated(this%calculator))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    gradients(:,:) = gradients + this%gradients

  end subroutine addGradients


  !> Returns the stress tensor.
  subroutine getStress(this, stress)

    !> Instance of DFTD3 data
    class(TDispDftD3), intent(inout) :: this

    !> stress tensor from the dispersion
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(allocated(this%calculator))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))

    stress(:,:) = this%stress

  end subroutine getStress


  !> Estimates the real space cutoff of the dispersion interaction.
  function getRCutoff(this) result(cutoff)

    !> Instance of DFTD3 data
    class(TDispDftD3), intent(inout) :: this

    !> Resulting cutoff
    real(dp) :: cutoff

    ! Since dftd3-routine uses its own real space summation routines, we do not need any real space
    ! neighbours from neighbour list, unless H-H repulsion is active.
    if (this%tHHRepulsion) then
      cutoff = hhRepCutOff
    else
      cutoff = 0.0_dp
    end if

  end function getRCutoff


  !> Add the additional H-H repulsion for D3H5
  subroutine addHHRepulsion(this, coords, neigh, img2CentCell)

    !> Instance of DFTD3 data
    class(TDispDftD3), intent(inout) :: this

    !> Current coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neigh

    !> Updated mapping to central cell.
    integer, intent(in) :: img2CentCell(:)

    ! Parameters (as published in dx.doi.org/10.1021/ct200751e )
    real(dp), parameter :: kk = 0.30_dp * kcal_mol__Hartree  ! s_HH
    real(dp), parameter :: ee = 14.31_dp  ! e_HH
    real(dp), parameter :: r0 = 2.35_dp * AA__Bohr  ! r_0,HH

    integer :: iAt1, iNeigh, iAt2, iAt2f, ii
    real(dp) :: rr, repE, dEdR, dCdR(3), cellVol, stressTmp(3,3), vect(3), prefac
    integer, allocatable :: nNeigh(:)

    if (this%tPeriodic) then
      stressTmp(:,:) = 0.0_dp
    end if

    allocate(nNeigh(this%nAtom))
    call getNrOfNeighboursForAll(nNeigh, neigh, HHRepCutOff)

    repE = 0.0_dp
    do iAt1 = 1, this%nAtom
      if (this%izp(iAt1) /= 1) then
        cycle
      end if
      do iNeigh = 1, nNeigh(iAt1)
        iAt2 = neigh%iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        if (this%izp(iAt2f) /= 1) then
          cycle
        end if
        if (iAt1 == iAt2f) then
          prefac = 0.5_dp
        else
          prefac = 1.0_dp
        end if
        vect(:) = coords(:,iAt1) - coords(:,iAt2)
        rr = sqrt(sum(vect**2))
        repE = repE + kk * (1.0_dp - 1.0_dp / (1.0_dp + exp(-ee * (rr / r0 - 1.0_dp))))
        dEdR = -kk * ee * exp(-ee * (rr / r0 -1.0_dp))&
            & / (r0 * (1.0_dp + exp(-ee * (rr / r0 - 1.0_dp)))**2)
        dCdR(:) = (coords(:,iAt1) - coords(:,iAt2)) / rr
        this%gradients(:,iAt1) = this%gradients(:,iAt1) + dEdR * dCdR
        this%gradients(:,iAt2f) = this%gradients(:,iAt2f) - dEdR * dCdR

        if (this%tPeriodic) then
          do ii = 1, 3
            stressTmp(:, ii) = stressTmp(:, ii) + prefac * vect(ii) * dEdR * dCdR
          end do
        end if
      end do
    end do

    this%dispE = this%dispE + repE

    if (this%tPeriodic) then
      cellVol = abs(determinant33(this%latVecs))
      this%stress = this%stress - stressTmp / cellVol
    end if

  end subroutine addHHRepulsion

end module dftbp_dispdftd3
