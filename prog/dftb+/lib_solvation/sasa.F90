!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Non-polar solvent accessible surface area (SASA) contributions
module dftbp_sasa
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_blasroutines, only : gemv
  use dftbp_constants, only : pi
  use dftbp_commontypes, only : TOrbitals
  use dftbp_environment, only : TEnvironment
  use dftbp_lebedev, only : getAngGrid, gridSize
  use dftbp_message, only : error
  use dftbp_periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use dftbp_schedule, only : distributeRangeInChunks, assembleChunks
  use dftbp_simplealgebra, only : determinant33
  use dftbp_solvation, only : TSolvation
  implicit none
  private

  public :: TSASACont, TSASAInput, TSASACont_init
  public :: writeSASAContInfo


  !> Input parameters to initialize the solvent accessible surface area model
  type :: TSASAInput

    !> Real space cutoff
    real(dp) :: sOffset = 0.0_dp

    !> Grid for numerical integration of atomic surfaces
    integer :: gridSize

    !> Probe radius for solvent molecules in atomic surfaces calculations
    real(dp) :: probeRad

    !> Smoothing dielectric function parameters
    real(dp) :: smoothingPar

    !> Real space cut-offs for surface area contribution
    real(dp) :: tolerance

    !> Van-der-Waals radii
    real(dp), allocatable :: vdwRad(:)

    !> Surface tension for each species
    real(dp), allocatable :: surfaceTension(:)

  end type TSASAInput


  !> Data for the solvent accessible surface area model
  type, extends(TSolvation) :: TSASACont

    !> number of atoms
    integer :: nAtom = 0

    !> solvation free energy
    real(dp), allocatable :: energies(:)

    !> lattice vectors if periodic
    real(dp) :: latVecs(3, 3) = 0.0_dp

    !> Volume of the unit cell
    real(dp) :: volume = 0.0_dp

    !> stress tensor
    real(dp) :: stress(3, 3) = 0.0_dp

    !> is this periodic
    logical :: tPeriodic

    !> are the coordinates current?
    logical :: tCoordsUpdated = .false.

    !> are the charges current?
    logical :: tChargesUpdated = .false.

    !> Real space cutoff
    real(dp) :: sCutoff = 0.0_dp

    !> Real space cut-offs for surface area contribution
    real(dp) :: tolerance

    !> Smoothing dielectric function parameters
    real(dp) :: smoothingPar(3)

    !> Van-der-Waals radii + probe radius
    real(dp), allocatable :: probeRad(:)

    !> Surface tension for each species
    real(dp), allocatable :: surfaceTension(:)

    !> Thresholds for smooth numerical integration
    real(dp), allocatable :: thresholds(:, :)

    !> Radial weight
    real(dp), allocatable :: radWeight(:)

    !> Angular grid for surface integration
    real(dp), allocatable :: angGrid(:, :)

    !> Weights of grid points for surface integration
    real(dp), allocatable :: angWeight(:)

    !> Solvent accessible surface area
    real(dp), allocatable :: sasa(:)

    !> Derivative of solvent accessible surface area w.r.t. coordinates
    real(dp), allocatable :: dsdr(:, :, :)

    !> Derivative of solvent accessible surface area w.r.t. strain deformations
    real(dp), allocatable :: dsdL(:, :, :)

  contains

    !> update internal copy of coordinates
    procedure :: updateCoords

    !> update internal copy of lattice vectors
    procedure :: updateLatVecs

    !> get real space cutoff
    procedure :: getRCutoff

    !> get energy contributions
    procedure :: getEnergies

    !> get force contributions
    procedure :: addGradients

    !> get stress tensor contributions
    procedure :: getStress

    !> Updates with changed charges for the instance
    procedure :: updateCharges

    !> Returns shifts per atom
    procedure :: getShifts

  end type TSASACont


contains


  !> Initialize solvent accessible surface area model from input data
  subroutine TSASACont_init(this, input, nAtom, species0, speciesNames, latVecs)

    !> Initialised instance at return
    type(TSASACont), intent(out) :: this

    !> Specific input parameters for solvent accessible surface area model
    type(TSASAInput), intent(in) :: input

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    integer :: iAt1, iSp1, nSpecies, stat

    nSpecies = size(speciesNames)

    this%tPeriodic = present(latVecs)
    if (this%tPeriodic) then
      call this%updateLatVecs(LatVecs)
    end if
    this%nAtom = nAtom

    allocate(this%energies(nAtom))
    allocate(this%sasa(nAtom))
    allocate(this%dsdr(3, nAtom, nAtom))
    allocate(this%dsdL(3, 3, nAtom))

    allocate(this%angGrid(3, gridSize(input%gridSize)))
    allocate(this%angWeight(gridSize(input%gridSize)))
    call getAngGrid(input%gridSize, this%angGrid, this%angWeight, stat)
    if (stat /= 0) then
      call error("Could not initialize angular grid for SASA model")
    end if

    allocate(this%probeRad(nSpecies))
    allocate(this%surfaceTension(nAtom))
    allocate(this%thresholds(2, nSpecies))
    allocate(this%radWeight(nSpecies))
    this%probeRad(:) = input%probeRad + input%vdwRad
    do iAt1 = 1, nAtom
      iSp1 = species0(iAt1)
      this%surfaceTension(iAt1) = input%surfaceTension(iSp1) * 4.0e-5_dp * pi
    end do
    call getIntegrationParam(nSpecies, input%smoothingPar, this%smoothingPar, &
        & this%probeRad, this%thresholds, this%radWeight)

    this%sCutoff = 2.0_dp * (maxval(this%probeRad) + input%smoothingPar) + input%sOffset
    this%tolerance = input%tolerance

    this%tCoordsUpdated = .false.
    this%tChargesUpdated = .false.

  end subroutine TSASACont_init


  !> Get parameters for smooth numerical integration
  subroutine getIntegrationParam(nSpecies, w, ah, rad, thr, wrp)

    !> Number of species
    integer, intent(in) :: nSpecies

    !> Smoothing parameter
    real(dp), intent(in) :: w

    !> Smoothing parameters
    real(dp), intent(out) :: ah(3)

    !> Probe radii
    real(dp), intent(in) :: rad(:)

    !> Numerical thresholds
    real(dp), intent(out) :: thr(:, :)

    !> Radial weight
    real(dp), intent(out) :: wrp(:)

    integer :: iSp1
    real(dp) :: rm, rp

    ah(1) = 0.5_dp
    ah(2) = 3.0_dp/(4.0_dp*w)
    ah(3) = -1.0_dp/(4.0_dp*w*w*w)

    do iSp1 = 1, nSpecies
      rm = rad(iSp1) - w
      rp = rad(iSp1) + w
      thr(1, iSp1) = rm**2
      thr(2, iSp1) = rp**2
      wrp(iSp1) = (0.25_dp/w + 3.0_dp*ah(3)*(0.2_dp*rp*rp-0.5_dp*rp*rad(iSp1) &
          & + rad(iSp1)*rad(iSp1)/3.0_dp))*rp*rp*rp - (0.25_dp/w &
          & + 3.0_dp*ah(3)*(0.2_dp*rm*rm - 0.5_dp*rm*rad(iSp1) &
          & + rad(iSp1)*rad(iSp1)/3.0_dp))*rm*rm*rm
    end do

  end subroutine getIntegrationParam


  !> Print the solvation model used
  subroutine writeSASAContInfo(unit, solvation)

    !> Formatted unit for IO
    integer, intent(in) :: unit

    !> Solvation model
    type(TSASACont), intent(in) :: solvation

    write(unit, '(a, ":", t30, es14.6, 1x, a, t50, i14, 1x, a)') "Grid points", &
        & solvation%nAtom*real(size(solvation%angWeight, dim=1), dp), "total", &
        & size(solvation%angWeight, dim=1), "per atom"

  end subroutine writeSASAContInfo


  !> Update internal stored coordinates
  subroutine updateCoords(this, env, neighList, img2CentCell, coords, species0)

    !> Data structure
    class(TSASACont), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> List of neighbours to atoms
    type(TNeighbourList), intent(in) :: neighList

    !> Image to central cell atom index
    integer, intent(in) :: img2CentCell(:)

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    integer, allocatable :: nNeigh(:)

    allocate(nNeigh(this%nAtom))
    call getNrOfNeighboursForAll(nNeigh, neighList, this%sCutoff)

    call getSASA(env, nNeigh, neighList%iNeighbour, img2CentCell, species0, &
        & coords, this%tolerance, this%smoothingPar, this%probeRad, &
        & this%thresholds, this%radWeight, this%angGrid, this%angWeight, &
        & this%sasa, this%dsdr)
    this%energies(:) = this%sasa * this%surfaceTension

    this%tCoordsUpdated = .true.
    this%tChargesUpdated = .false.

  end subroutine updateCoords


  !> Update internal copy of lattice vectors
  subroutine updateLatVecs(this, latVecs)

    !> Data structure
    class(TSASACont), intent(inout) :: this

    !> Lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    @:ASSERT(this%tPeriodic)
    @:ASSERT(all(shape(latvecs) == shape(this%latvecs)))

    this%volume = abs(determinant33(latVecs))
    this%latVecs(:,:) = latVecs

    this%tCoordsUpdated = .false.
    this%tChargesUpdated = .false.

  end subroutine updateLatVecs


  !> Get energy contributions
  subroutine getEnergies(this, energies)

    !> data structure
    class(TSASACont), intent(inout) :: this

    !> energy contributions for each atom
    real(dp), intent(out) :: energies(:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(size(energies) == this%nAtom)

    energies(:) = this%energies

  end subroutine getEnergies


  !> Get force contributions
  subroutine addGradients(this, env, neighList, species, coords, img2CentCell, gradients)

    !> Data structure
    class(TSASACont), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Neighbour list
    type(TNeighbourList), intent(in) :: neighList

    !> Specie for each atom
    integer, intent(in) :: species(:)

    !> Coordinate of each atom
    real(dp), intent(in) :: coords(:,:)

    !> Mapping of atoms to cetnral cell
    integer, intent(in) :: img2CentCell(:)

    !> Gradient contributions for each atom
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    call gemv(gradients, this%dsdr, this%surfaceTension, beta=1.0_dp)

  end subroutine addGradients


  !> get stress tensor contributions
  subroutine getStress(this, stress)

    !> data structure
    class(TSASACont), intent(inout) :: this

    !> Stress tensor contributions
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))
    @:ASSERT(this%tPeriodic)
    @:ASSERT(this%volume > 0.0_dp)

    stress(:,:) = this%stress / this%volume

  end subroutine getStress


  !> Distance cut off for solvent accessible surface area calculations
  function getRCutoff(this) result(cutoff)

    !> data structure
    class(TSASACont), intent(inout) :: this

    !> resulting cutoff
    real(dp) :: cutoff

    cutoff = this%sCutoff

  end function getRCutoff


  !> Updates with changed charges for the instance.
  subroutine updateCharges(this, env, species, neighList, qq, q0, img2CentCell, orb)

    !> Data structure
    class(TSASACont), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Species, shape: [nAtom]
    integer, intent(in) :: species(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighList

    !> Orbital charges.
    real(dp), intent(in) :: qq(:,:,:)

    !> Reference orbital charges.
    real(dp), intent(in) :: q0(:,:,:)

    !> Mapping on atoms in central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    @:ASSERT(this%tCoordsUpdated)

    this%tChargesUpdated = .true.

  end subroutine updateCharges


  !> Returns shifts per atom
  subroutine getShifts(this, shiftPerAtom, shiftPerShell)

    !> Data structure
    class(TSASACont), intent(inout) :: this

    !> Shift per atom
    real(dp), intent(out) :: shiftPerAtom(:)

    !> Shift per shell
    real(dp), intent(out) :: shiftPerShell(:,:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(size(shiftPerAtom) == this%nAtom)
    @:ASSERT(size(shiftPerShell, dim=2) == this%nAtom)

    shiftPerAtom(:) = 0.0_dp
    shiftPerShell(:,:) = 0.0_dp

  end subroutine getShifts


  !> Calculate solvent accessible surface area for every atom
  subroutine getSASA(env, nNeighbour, iNeighbour, img2CentCell, &
      & species, coords, tolerance, smoothingPar, probeRad, thresholds, &
      & radWeight, angGrid, angWeight, sasa, dsdr)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist
    integer, intent(in) :: iNeighbour(0:, :)

    !> Mapping into the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Species, shape: [nAtom]
    integer, intent(in) :: species(:)

    !> Current atomic positions
    real(dp), intent(in) :: coords(:,:)

    !> Real space cut-offs for surface area contribution
    real(dp) :: tolerance

    !> Smoothing dielectric function parameters
    real(dp) :: smoothingPar(3)

    !> Van-der-Waals radii + probe radius
    real(dp), allocatable :: probeRad(:)

    !> Thresholds for smooth numerical integration
    real(dp), allocatable :: thresholds(:, :)

    !> Radial weight
    real(dp), allocatable :: radWeight(:)

    !> Angular grid for surface integration
    real(dp), allocatable :: angGrid(:, :)

    !> Weights of grid points for surface integration
    real(dp), allocatable :: angWeight(:)

    !> Solvent accessible surface area
    real(dp), allocatable :: sasa(:)

    !> Derivative of solvent accessible surface area w.r.t. coordinates
    real(dp), allocatable :: dsdr(:, :, :)

    integer :: iAt1, iSp1, iAt2, iAt2f, iSp2, iNeigh, ip
    integer :: iAtFirst, iAtLast, nAtom, mNeighbour, nEval
    real(dp) :: vec(3), dist2, dist
    real(dp) :: uj, ah3uj2
    real(dp) :: sasaij, dsasaij
    real(dp) :: rsas, sasai, point(3), sasap, wsa, dGr(3)
    real(dp), allocatable :: grds(:,:), derivs(:,:)
    integer, allocatable :: grdi(:)

    nAtom = size(nNeighbour)
    mNeighbour = maxval(nNeighbour)

    allocate(derivs(3,nAtom))
    allocate(grds(3,mNeighbour))
    allocate(grdi(mNeighbour))

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)
    sasa(:) = 0.0_dp
    dsdr(:, :, :) = 0.0_dp

    !$omp parallel do default(none) schedule(runtime) reduction(+:sasa, dsdr) &
    !$omp shared(iAtFirst, iAtLast, species, coords, probeRad, angGrid) &
    !$omp shared(nNeighbour, iNeighbour, img2CentCell, thresholds, smoothingPar) &
    !$omp private(iAt1, iSp1, rsas, derivs, sasai, ip, point, nEval, grds, grdi) &
    !$omp private(sasap, iNeigh, iAt2, iAt2f, iSp2, vec, dist2, dist, uj, ah3uj2) &
    !$omp private(dsasaij, sasaij, wsa, dGr) shared(radWeight, angWeight, tolerance)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)

      rsas = probeRad(iSp1)

      ! initialize storage
      derivs(:, :) = 0.0_dp
      sasai = 0.0_dp

      ! loop over grid points
      do ip = 1, size(angGrid, dim=2)
        ! grid point position
        point(:) = coords(:,iAt1) + rsas * angGrid(:, ip)

        ! atomic surface function at the grid point
        nEval = 0
        grds(:, :) = 0.0_dp
        grdi(:) = 0
        sasap = 1.0_dp
        do iNeigh = 1, nNeighbour(iAt1)
          iAt2 = iNeighbour(iNeigh, iAt1)
          iAt2f = img2CentCell(iAt2)
          iSp2 = species(iAt2f)
          ! compute the distance to the atom
          vec(:) = point(:) - coords(:, iAt2)
          dist2 = vec(1)**2 + vec(2)**2 + vec(3)**2
          ! if within the outer cut-off compute
          if (dist2 < thresholds(2,iSp2)) then
            if (dist2 < thresholds(1,iSp2)) then
              sasap = 0.0_dp
              exit
            else
              dist = sqrt(dist2)
              uj = dist - probeRad(iSp2)
              ah3uj2 = smoothingPar(3)*uj*uj
              dsasaij = smoothingPar(2)+3.0_dp*ah3uj2
              sasaij =  smoothingPar(1)+(smoothingPar(2)+ah3uj2)*uj

              ! accumulate the molecular surface
              sasap = sasap*sasaij
              ! compute the gradient wrt the neighbor
              dsasaij = dsasaij/(sasaij*dist)
              nEval = nEval+1
              grdi(nEval) = iAt2f
              grds(:, nEval) = dsasaij*vec(:)
            end if
          end if
        end do

        if (sasap > tolerance) then
          ! numerical quadrature weight
          wsa = angWeight(ip)*radWeight(iSp1)*sasap
          ! accumulate the surface area
          sasai = sasai + wsa
          ! accumulate the surface gradient
          do iNeigh = 1, nEval
            iAt2f = grdi(iNeigh)
            dGr(:) = wsa * grds(:, iNeigh)
            derivs(:, iAt1) = derivs(:, iAt1) + dGr(:)
            derivs(:, iAt2f) = derivs(:, iAt2f) - dGr(:)
          end do
        end if
      end do

      sasa(iAt1) = sasai
      dsdr(:,:,iAt1) = derivs

    end do
    !$omp end parallel do

    call assembleChunks(env, sasa)
    call assembleChunks(env, dsdr)

  end subroutine getSASA


end module dftbp_sasa
