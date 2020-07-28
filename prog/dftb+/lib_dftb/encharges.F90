!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implementation of the electronegativity equilibration charge model used
!> for the charge scaling in DFT-D4.
!>
!> This implementation is general enough to be used outside of DFT-D4.
module dftbp_encharges
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_constants, only : pi
  use dftbp_errorfunction, only : erfwrap
  use dftbp_coulomb, only : ewaldReal, ewaldReciprocal, derivStressEwaldRec, &
      & getMaxGEwald, getOptimalAlphaEwald
  use dftbp_coordnumber, only : TCNCont, TCNInput, init
  use dftbp_blasroutines, only : hemv, gemv, gemm
  use dftbp_lapackroutines, only : symmatinv
  use dftbp_periodic, only : TNeighbourList, getNrOfNeighboursForAll, getLatticePoints
  use dftbp_simplealgebra, only : determinant33, invert33
  implicit none
  private

  public :: TEeqCont, TEeqInput, init
  public :: getEEQcharges

  real(dp), parameter :: sqrtpi = sqrt(pi)
  real(dp), parameter :: sqrt2pi = sqrt(2.0_dp/pi)


  !> EEQ parametrisation data
  type :: TEeqParam

    !> Electronegativities for EEQ model
    real(dp), allocatable :: chi(:)

    !> Chemical hardnesses for EEQ model
    real(dp), allocatable :: gam(:)

    !> Charge widths for EEQ model
    real(dp), allocatable :: rad(:)

    !> CN scaling for EEQ model
    real(dp), allocatable :: kcn(:)

  end type


  !> Input for the EEQ model
  type, extends(TEeqParam) :: TEeqInput

    !> Net charge of the system.
    real(dp) :: nrChrg

    !> Parameter for Ewald summation.
    real(dp) :: parEwald

    !> Ewald tolerance
    real(dp) :: tolEwald

    !> Cutoff for real-space summation under PBCs
    real(dp) :: cutoff

    !> Input for coordination number
    type(TCNInput) :: cnInput

  end type TEeqInput


  !> Container containing all calculation data for the EEQ model
  type :: TEeqCont

    !> number of atoms
    integer :: nAtom

    !> EEQ parametrisation
    type(TEeqParam) :: param

    !> Net charge of the system.
    real(dp) :: nrChrg

    !> Cutoff for real-space summation under PBCs
    real(dp) :: cutoff

    !> lattice vectors if periodic
    real(dp) :: latVecs(3, 3)

    !> Volume of the unit cell
    real(dp) :: vol

    !> evaluate Ewald parameter
    logical :: tAutoEwald

    !> Ewald tolerance
    real(dp) :: tolEwald

    !> Parameter for Ewald summation.
    real(dp) :: parEwald

    !> is this periodic
    logical :: tPeriodic

    !> Contains the points included in the reciprocal sum.
    !> The set should not include the origin or inversion related points.
    real(dp), allocatable :: recPoint(:, :)

    !> Coordination number container
    type(TCNCont) :: cnCont

    !> are the coordinates current?
    logical :: tCoordsUpdated

    !> electrostatic energy
    real(dp), allocatable :: energies(:)

    !> force contributions
    real(dp), allocatable :: gradients(:, :)

    !> stress tensor
    real(dp), allocatable :: stress(:, :)

    !> partial charges
    real(dp), allocatable :: charges(:)

    !> derivative of partial charges w.r.t. coordinates
    real(dp), allocatable :: dqdr(:, :, :)

    !> derivative of partial charges w.r.t. strain deformations
    real(dp), allocatable :: dqdL(:, :, :)

  contains

    !> update internal store of coordinates
    procedure :: updateCoords

    !> update internal store of lattice vectors
    procedure :: updateLatVecs

    !> return energy contribution
    procedure :: getEnergies

    !> return energy contribution
    procedure :: getCharges

    !> return force contribution
    generic :: addGradients => addGradientsEnergy, addGradientsCharges

    !> return force contribution for electrostatic energy
    procedure, private :: addGradientsEnergy

    !> return force contribution for charge derivatives
    procedure, private :: addGradientsCharges

    !> return stress tensor contribution
    generic :: addStress => addStressEnergy, addStressCharges

    !> return stress tensor contribution for electrostatic energy
    procedure, private :: addStressEnergy

    !> return stress tensor contribution for charge derivatives
    procedure, private :: addStressCharges

    !> cutoff distance in real space for EEQ
    procedure :: getRCutoff

  end type TEeqCont


  !> initialize container from input
  interface init
    module procedure :: initialize
  end interface init


contains


  subroutine initialize(this, input, withEnergies, withCharges, nAtom, latVecs)

    !> Instance of EEQ container
    type(TEeqCont), intent(out) :: this

    !> Input for the EEQ model
    type(TEeqInput), intent(in) :: input

    !> Enable electrostatic energy calculation + derivatives
    logical, intent(in) :: withEnergies

    !> Enable partial charge calculation + derivatives
    logical, intent(in) :: withCharges

    !> Nr. of atoms in the system.
    integer, intent(in) :: nAtom

    !> Lattice vectors, if the system is periodic.
    real(dp), intent(in), optional :: latVecs(:, :)

    this%tPeriodic = present(latVecs)

    if (this%tPeriodic) then
      call init(this%cnCont, input%cnInput, nAtom, latVecs)
    else
      call init(this%cnCont, input%cnInput, nAtom)
    end if

    this%param = input%TEeqParam
    this%cutoff = input%cutoff
    this%nrChrg = input%nrChrg

    if (this%tPeriodic) then
      this%tAutoEwald = input%parEwald <= 0.0_dp
      this%tolEwald = input%tolEwald
      this%parEwald = input%parEwald

      call this%updateLatVecs(latVecs)
    end if

    this%nAtom = nAtom

    if (withEnergies) then
      allocate(this%energies(nAtom))
      allocate(this%gradients(3, nAtom))
      allocate(this%stress(3, 3))
    end if

    if (withCharges) then
      allocate(this%charges(nAtom))
      allocate(this%dqdr(3, nAtom, nAtom))
      allocate(this%dqdL(3, 3, nAtom))
    end if

    this%tCoordsUpdated = .false.

  end subroutine initialize


  !> Notifies the objects about changed coordinates.
  subroutine updateCoords(this, neigh, img2CentCell, coords, species)

    !> Instance of EEQ container
    class(TEeqCont), intent(inout) :: this

    !> Updated neighbour list.
    type(TNeighbourList), intent(in) :: neigh

    !> Updated mapping to central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Updated coordinates.
    real(dp), intent(in) :: coords(:,:)

    !> Species of the atoms in the unit cell.
    integer, intent(in) :: species(:)

    integer, allocatable :: nNeigh(:)

    call this%cnCont%updateCoords(neigh, img2CentCell, coords, species)

    allocate(nNeigh(this%nAtom))
    call getNrOfNeighboursForAll(nNeigh, neigh, this%cutoff)

    call getEEQCharges(this%nAtom, coords, species, this%nrChrg, nNeigh, &
        & neigh%iNeighbour, neigh%neighDist2, img2CentCell, this%recPoint, this%parEwald, &
        & this%vol, this%param%chi, this%param%kcn, this%param%gam, this%param%rad, &
        & this%cnCont%cn, this%cnCont%dcndr, this%cnCont%dcndL, &
        & this%energies, this%gradients, this%stress, &
        & this%charges, this%dqdr, this%dqdL)

    this%tCoordsUpdated = .true.

  end subroutine updateCoords


  !> Notifies the object about updated lattice vectors.
  subroutine updateLatVecs(this, latVecs)

    !> Instance of EEQ container
    class(TEeqCont), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:, :)

    real(dp) :: recVecs(3, 3), maxGEwald

    @:ASSERT(this%tPeriodic)
    @:ASSERT(all(shape(latVecs) == shape(this%latVecs)))

    this%latVecs(:, :) = latVecs
    this%vol = determinant33(latVecs)
    call invert33(recVecs, latVecs, this%vol)
    this%vol = abs(this%vol)
    recVecs(:, :) = 2.0_dp * pi * transpose(recVecs)

    if (this%tAutoEwald) then
      this%parEwald = getOptimalAlphaEwald(latVecs, recVecs, this%vol, this%tolEwald)
    end if
    maxGEwald = getMaxGEwald(this%parEwald, this%vol, this%tolEwald)
    call getLatticePoints(this%recPoint, recVecs, latVecs/(2.0_dp*pi), maxGEwald,&
        & onlyInside=.true., reduceByInversion=.true., withoutOrigin=.true.)
    this%recPoint(:, :) = matmul(recVecs, this%recPoint)

    call this%cnCont%updateLatVecs(LatVecs)

    this%tCoordsUpdated = .false.

  end subroutine updateLatVecs


  !> Returns the atomic resolved energies due to the EEQ.
  subroutine getEnergies(this, energies)

    !> Instance of EEQ container
    class(TEeqCont), intent(inout) :: this

    !> Contains the atomic energy contributions on exit.
    real(dp), intent(out) :: energies(:)

    @:ASSERT(allocated(this%energies))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(size(energies) == this%nAtom)

    energies(:) = this%energies

  end subroutine getEnergies


  !> Returns the atomic resolved energies due to the EEQ.
  subroutine getCharges(this, charges)

    !> Instance of EEQ container
    class(TEeqCont), intent(inout) :: this

    !> Contains the atomic partial charges on exit.
    real(dp), intent(out) :: charges(:)

    @:ASSERT(allocated(this%charges))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(size(charges) == this%nAtom)

    charges(:) = this%charges

  end subroutine getCharges


  !> Adds the atomic gradients to the provided vector.
  subroutine addGradientsEnergy(this, gradients)

    !> Instance of EEQ container
    class(TEeqCont), intent(inout) :: this

    !> The vector to increase by the gradients.
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(allocated(this%gradients))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    gradients(:,:) = gradients + this%gradients

  end subroutine addGradientsEnergy


  !> Adds the atomic gradients to the provided vector.
  subroutine addGradientsCharges(this, gradients, dEdq, beta)

    !> Instance of EEQ container
    class(TEeqCont), intent(inout) :: this

    !> The vector to increase by the gradients.
    real(dp), intent(inout) :: gradients(:,:)

    !> Derivative of the energy expression w.r.t. the partial charges
    real(dp), intent(in) :: dEdq(:)

    !> Optional scaling factor
    real(dp), intent(in), optional :: beta
    real(dp) :: beta0

    if (present(beta)) then
      beta0 = beta
    else
      beta0 = 1.0_dp
    end if

    @:ASSERT(allocated(this%dqdr))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))
    @:ASSERT(size(dEdq) == this%nAtom)

    call gemv(gradients, this%dqdr, dEdq, beta=beta0)

  end subroutine addGradientsCharges


  !> Returns the stress tensor.
  subroutine addStressEnergy(this, stress)

    !> Instance of EEQ container
    class(TEeqCont), intent(inout) :: this

    !> stress tensor from the EEQ
    real(dp), intent(inout) :: stress(:,:)

    @:ASSERT(allocated(this%stress))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))

    stress(:,:) = stress(:,:) + this%stress

  end subroutine addStressEnergy


  !> Returns the stress tensor.
  subroutine addStressCharges(this, stress, dEdq, beta)

    !> Instance of EEQ container
    class(TEeqCont), intent(inout) :: this

    !> stress tensor from the EEQ
    real(dp), intent(inout) :: stress(:,:)

    !> Derivative of the energy expression w.r.t. the partial charges
    real(dp), intent(in) :: dEdq(:)

    !> Optional scaling factor
    real(dp), intent(in), optional :: beta
    real(dp) :: beta0

    if (present(beta)) then
      beta0 = beta
    else
      beta0 = 1.0_dp
    end if

    @:ASSERT(allocated(this%dqdL))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))
    @:ASSERT(size(dEdq) == this%nAtom)

    call gemv(stress, this%dqdL, dEdq, beta=beta0)

  end subroutine addStressCharges


  !> Estimates the real space cutoff of the EEQ interaction.
  function getRCutoff(this) result(cutoff)

    !> Instance of EEQ container
    class(TEeqCont), intent(inout) :: this

    !> Resulting cutoff
    real(dp) :: cutoff

    cutoff = max(this%cutoff, this%cnCont%getRCutoff())

  end function getRCutoff

  !> Generates full interaction matrix for Gaussian charge distributions.
  subroutine getCoulombMatrixCluster(nAtom, coords, species, gam, rad, aMat)

    !> number of atoms
    integer, intent(in) :: nAtom

    !> List of atomic coordinates.
    real(dp), intent(in) :: coords(:, :)

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> element-specific chemical hardnesses
    real(dp), intent(in) :: gam(:)

    !> element-specific charge widths / atomic radii
    real(dp), intent(in) :: rad(:)

    !> Interaction Matrix for each atom pair.
    real(dp), intent(out) :: aMat(:, :)

    integer :: iAt1, iSp1, iAt2f, iSp2
    real(dp) :: dist, vec(3), eta12

    aMat(:, :) = 0.0_dp

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP SHARED(nAtom, species, gam, rad, coords, aMat) &
    !$OMP PRIVATE(iSp1, iAt2f, iSp2, vec, dist, eta12)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      aMat(iAt1, iAt1) = gam(iSp1) + sqrt2pi / rad(iSp1)
      do iAt2f = iAt1 + 1, nAtom
        iSp2 = species(iAt2f)
        vec(:) = coords(:, iAt1) - coords(:, iAt2f)
        dist = sqrt(sum(vec**2))
        eta12 = 1.0_dp / sqrt(rad(iSp1)**2 + rad(iSp2)**2)
        aMat(iAt2f, iAt1) = erfwrap(eta12 * dist) / dist
        aMat(iAt1, iAt2f) = aMat(iAt2f, iAt1)
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine getCoulombMatrixCluster


  !> Derivatives of interaction matrix for Gaussian charge distributions.
  subroutine getCoulombDerivsCluster(nAtom, coords, species, rad, qVec, aMatdr, aTrace)

    !> Number of atoms.
    integer, intent(in) :: nAtom

    real(dp), intent(in) :: coords(:, :)

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> element-specific charge widths / atomic radii
    real(dp), intent(in) :: rad(:)

    !> List of charges on each atom.
    real(dp), intent(in) :: qVec(:)

    !> Contains the derivative on exit.
    real(dp), intent(out) :: aMatdr(:, :, :)

    !> Contains the `trace' derivative on exit.
    real(dp), intent(out) :: aTrace(:, :)

    integer :: iAt1, iAt2f, iSp1, iSp2
    real(dp) :: dist, vec(3), aTmp(3), arg, eta12

    aMatdr(:,:,:) = 0.0_dp
    aTrace(:,:) = 0.0_dp

    !$OMP PARALLEL DEFAULT(NONE) REDUCTION(+:aTrace, aMatdr) &
    !$OMP SHARED(nAtom, species, rad, coords, qVec) &
    !$OMP PRIVATE(iAt2f, iSp1, iSp2, vec, dist, aTmp, arg, eta12)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      do iAt2f = 1, iAt1 - 1
        iSp2 = species(iAt2f)
        vec(:) = coords(:, iAt1) - coords(:, iAt2f)
        dist = sqrt(sum(vec**2))
        eta12 = 1.0_dp / sqrt(rad(iSp1)**2 + rad(iSp2)**2)
        arg = dist * eta12
        aTmp = vec * (2.0_dp * eta12 / sqrtpi * exp(-arg*arg) - erfwrap(arg)/dist) / dist**2
        aTrace(:, iAt1) = aTrace(:, iAt1) + aTmp * qVec(iAt2f)
        aTrace(:, iAt2f) = aTrace(:, iAt2f) - aTmp * qVec(iAt1)
        aMatdr(:, iAt1, iAt2f) = aMatdr(:, iAt1, iAt2f) + aTmp * qVec(iAt1)
        aMatdr(:, iAt2f, iAt1) = aMatdr(:, iAt2f, iAt1) - aTmp * qVec(iAt2f)
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine getCoulombDerivsCluster


  !> Generates full interaction matrix for Gaussian charge distributions.
  subroutine getCoulombMatrixPeriodic(nAtom, coords, species, nNeighbour, iNeighbour, neighDist2,&
      & img2CentCell, recPoint, gam, rad, alpha, volume, aMat)

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coords(:, :)

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist.
    integer, intent(in) :: iNeighbour(0:, :)

    !> Square distances of the neighbours.
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Mapping into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Contains the points included in the reciprocal sum.
    !  The set should not include the origin or inversion related points.
    real(dp), intent(in) :: recPoint(:, :)

    !> element-specific chemical hardnesses
    real(dp), intent(in) :: gam(:)

    !> element-specific charge widths / atomic radii
    real(dp), intent(in) :: rad(:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> Matrix of 1/R values for each atom pair.
    real(dp), intent(out) :: aMat(:, :)

    aMat(:, :) = 0.0_dp

    ! Real space part of the Ewald sum.
    call addRealSpaceContribs(nAtom, species, nNeighbour, iNeighbour, neighDist2, img2CentCell,&
        & gam, rad, alpha, aMat)

    ! Reciprocal space part of the Ewald sum.
    call addEwaldContribs(nAtom, coords, recPoint, alpha, volume, aMat)

  end subroutine getCoulombMatrixPeriodic


  !> Derivatives of interaction matrix for Gaussian charge distributions.
  subroutine getCoulombDerivsPeriodic(nAtom, coords, species, nNeighbour, iNeighbour, neighDist2,&
      & img2CentCell, recPoint, alpha, volume, rad, qVec, aMatdr, aMatdL, aTrace)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coords(:, :)

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist.
    integer, intent(in) :: iNeighbour(0:, :)

    !> Square distances of the neighbours.
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Mapping into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Contains the points included in the reciprocal sum.
    !  The set should not include the origin or inversion related points.
    real(dp), intent(in) :: recPoint(:, :)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> element-specific charge widths / atomic radii
    real(dp), intent(in) :: rad(:)

    !> List of charges on each atom
    real(dp), intent(in) :: qVec(:)

    !> Contains the derivative on exit.
    real(dp), intent(out) :: aMatdr(:, :, :)

    !> Contains the strain derivative on exit.
    real(dp), intent(out) :: aMatdL(:, :, :)

    !> Contains the `trace' derivative on exit.
    real(dp), intent(out) :: aTrace(:, :)

    @:ASSERT(volume > 0.0_dp)

    aMatdr(:, :, :) = 0.0_dp
    aMatdL(:, :, :) = 0.0_dp
    aTrace(:, :) = 0.0_dp

    ! d(1/R)/dr real space
    call addRealSpaceDerivs(nAtom, coords, species, nNeighbour, iNeighbour, neighDist2,&
        & img2CentCell, alpha, rad, qVec, aMatdr, aMatdL, aTrace)

    ! d(1/R)/dr reciprocal space
    call addEwaldDerivs(nAtom, coords, recPoint, alpha, volume, qVec, aMatdr, aMatdL, aTrace)

  end subroutine getCoulombDerivsPeriodic


  !> Reciprocal space contributions to interaction matrix.
  subroutine addEwaldContribs(nAtom, coords, recPoint, alpha, volume, aMat)

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coords(:, :)

    !> Contains the points included in the reciprocal sum.
    !  The set should not include the origin or inversion related points.
    real(dp), intent(in) :: recPoint(:, :)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> Interaction Matrix for each atom pair.
    real(dp), intent(inout) :: aMat(:, :)

    integer :: iAt1, iAt2f
    real(dp) :: vec(3), rTerm

    @:ASSERT(volume > 0.0_dp)

    ! Reciprocal space part of the Ewald sum.
    ! Workaround:nagfor 7.0 with combined DO and PARALLEL
    !$OMP PARALLEL DO DEFAULT(NONE) REDUCTION(+:aMat)&
    !$OMP& SHARED(nAtom, alpha, volume, coords, recPoint)&
    !$OMP& PRIVATE(iAt2f, vec, rTerm)&
    !$OMP& SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      aMat(iAt1, iAt1) = aMat(iAt1, iAt1) - alpha / sqrt(pi) + pi / (volume * alpha**2)
      do iAt2f = iAt1, nAtom
        vec(:) = coords(:, iAt1)-coords(:, iAt2f)
        rTerm = ewaldReciprocal(vec, recPoint, alpha, volume) - pi / (volume * alpha**2)
        aMat(iAt2f, iAt1) = aMat(iAt2f, iAt1) + rTerm
        aMat(iAt1, iAt2f) = aMat(iAt1, iAt2f) + rTerm
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine addEwaldContribs


  !> Reciprocal space contributions to derivative of interaction matrix.
  subroutine addEwaldDerivs(nAtom, coords, recPoint, alpha, volume, qVec, aMatdr, aMatdL, aTrace)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coords(:, :)

    !> Contains the points included in the reciprocal sum.
    !  The set should not include the origin or inversion related points.
    real(dp), intent(in) :: recPoint(:, :)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> List of charges on each atom
    real(dp), intent(in) :: qVec(:)

    !> Contains the derivative on exit.
    real(dp), intent(inout) :: aMatdr(:, :, :)

    !> Contains the strain derivative on exit.
    real(dp), intent(inout) :: aMatdL(:, :, :)

    !> Contains the `trace' derivative on exit.
    real(dp), intent(inout) :: aTrace(:, :)

    integer :: iAt1, iAt2f
    real(dp) :: vec(3), aTmp(3), sigma(3, 3)

    !$OMP PARALLEL DEFAULT(NONE) REDUCTION(+:aTrace, aMatdr, aMatdL) &
    !$OMP SHARED(nAtom, coords, recPoint, alpha, volume, qVec) &
    !$OMP PRIVATE(iAt2f, vec, aTmp, sigma)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      do iAt2f = iAt1, nAtom
        vec(:) = coords(:, iAt1) - coords(:, iAt2f)
        call derivStressEwaldRec(vec, recPoint, alpha, volume, aTmp, sigma)
        aTrace(:, iAt1) = aTrace(:, iAt1) + aTmp * qVec(iAt2f)
        aTrace(:, iAt2f) = aTrace(:, iAt2f) - aTmp * qVec(iAt1)
        aMatdr(:, iAt1, iAt2f) = aMatdr(:, iAt1, iAt2f) + aTmp * qVec(iAt1)
        aMatdr(:, iAt2f, iAt1) = aMatdr(:, iAt2f, iAt1) - aTmp * qVec(iAt2f)
        aMatdL(:, :, iAt1) = aMatdL(:, :, iAt1) + sigma * qVec(iAt2f)
        if (iAt1 /= iAt2f) then
          aMatdL(:, :, iAt2f) = aMatdL(:, :, iAt2f) + sigma * qVec(iAt1)
        end if
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine addEwaldDerivs


  !> Real space contributions to interaction matrix.
  subroutine addRealSpaceContribs(nAtom, species, nNeighbour, iNeighbour, neighDist2, img2CentCell,&
      & gam, rad, alpha, aMat)

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist.
    integer, intent(in) :: iNeighbour(0:, :)

    !> Square distances of the neighbours.
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Mapping into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> element-specific chemical hardnesses
    real(dp), intent(in) :: gam(:)

    !> element-specific charge widths / atomic radii
    real(dp), intent(in) :: rad(:)

    !> Interaction Matrix for each atom pair.
    real(dp), intent(inout) :: aMat(:, :)

    real(dp) :: dist, eta12, rTerm
    integer :: iAt1, iAt2, iAt2f, iSp1, iSp2, iNeigh

    ! Workaround:nagfor 7.0 with combined DO and PARALLEL
    !$OMP PARALLEL DO DEFAULT(NONE) REDUCTION(+:aMat)&
    !$OMP& SHARED(nAtom, species, gam, rad, nNeighbour, iNeighbour, img2CentCell, neighDist2)&
    !$OMP& SHARED(alpha)&
    !$OMP& PRIVATE(iNeigh, iAt2, iAt2f, iSp1, iSp2, dist, rTerm, eta12)&
    !$OMP& SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      aMat(iAt1, iAt1) = aMat(iAt1, iAt1) + gam(iSp1) + sqrt2pi/rad(iSp1)
      do iNeigh = 1, nNeighbour(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        dist = sqrt(neighDist2(iNeigh, iAt1))
        eta12 = 1.0_dp/sqrt(rad(iSp1)**2 + rad(iSp2)**2)
        rTerm = (erfwrap(eta12 * dist) - erfwrap(alpha * dist)) / dist
        aMat(iAt2f, iAt1) = aMat(iAt2f, iAt1) + rTerm
        aMat(iAt1, iAt2f) = aMat(iAt1, iAt2f) + rTerm
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine addRealSpaceContribs


  !> Real space contributions to derivative of interaction matrix.
  subroutine addRealSpaceDerivs(nAtom, coords, species, nNeighbour, iNeighbour, neighDist2,&
      & img2CentCell, alpha, rad, qVec, aMatdr, aMatdL, aTrace)

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coords(:, :)

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist.
    integer, intent(in) :: iNeighbour(0:, :)

    !> Square distances of the neighbours.
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Mapping into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> element-specific charge widths / atomic radii
    real(dp), intent(in) :: rad(:)

    !> List of charges on each atom.
    real(dp), intent(in) :: qVec(:)

    !> Contains the derivative on exit.
    real(dp), intent(inout) :: aMatdr(:, :, :)

    !> Contains the strain derivative on exit.
    real(dp), intent(inout) :: aMatdL(:, :, :)

    !> Contains the `trace' derivative on exit.
    real(dp), intent(inout) :: aTrace(:, :)

    real(dp) :: dist, vec(3), aTmp(3), arg, eta12, ewl, sigma(3, 3)
    integer :: iAt1, iAt2, iAt2f, iSp1, iSp2, iNeigh

    !$OMP PARALLEL DEFAULT(NONE) REDUCTION(+:aTrace, aMatdr, aMatdL) &
    !$OMP SHARED(nAtom, species, nNeighbour, iNeighbour, coords, img2CentCell) &
    !$OMP SHARED(neighDist2, rad, qVec, alpha) &
    !$OMP PRIVATE(iAt2, iAt2f, iSp1, iSp2, vec, dist, aTmp, arg, ewl, eta12, sigma)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeighbour(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        vec(:) = coords(:, iAt1) - coords(:, iAt2)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        dist = sqrt(neighDist2(iNeigh, iAt1))
        eta12 = 1.0_dp/sqrt(rad(iSp1)**2 + rad(iSp2)**2)
        arg = dist * eta12
        ewl = dist * alpha
        aTmp = ((2*eta12 / sqrtpi * exp(-arg *arg) - erfwrap(arg) / dist) &
            & - (2*alpha / sqrtpi * exp(-ewl *ewl) - erfwrap(ewl) / dist)) * vec / dist**2
        aTrace(:, iAt1) = aTrace(:, iAt1) + aTmp * qVec(iAt2f)
        aTrace(:, iAt2f) = aTrace(:, iAt2f) - aTmp * qVec(iAt1)
        aMatdr(:, iAt1, iAt2f) = aMatdr(:, iAt1, iAt2f) + aTmp * qVec(iAt1)
        aMatdr(:, iAt2f, iAt1) = aMatdr(:, iAt2f, iAt1) - aTmp * qVec(iAt2f)
        sigma(:,:) = spread(aTmp, 1, 3) * spread(vec, 2, 3)
        aMatdL(:, :, iAt1) = aMatdL(:, :, iAt1) + sigma * qVec(iAt2f)
        if (iAt1 /= iAt2f) then
          aMatdL(:, :, iAt2f) = aMatdL(:, :, iAt2f) + sigma * qVec(iAt1)
        end if
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine addRealSpaceDerivs


  !> Electronegativity equilibration charge model
  subroutine getEEQCharges(nAtom, coords, species, charge, nNeighbour, iNeighbour, neighDist2,&
      & img2CentCell, recPoint, alpha, volume, chi, kcn, gam, rad, cn, dcndr, dcndL, energies,&
      & gradients, stress, qAtom, dqdr, dqdL)

    !> number of atoms
    integer, intent(in) :: nAtom

    !> List of atomic coordinates.
    real(dp), intent(in) :: coords(:, :)

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> Total molecular charge.
    real(dp), intent(in) :: charge

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist.
    integer, intent(in) :: iNeighbour(0:, :)

    !> Square distances of the neighbours.
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Mapping into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Contains the points included in the reciprocal sum.
    !  The set should not include the origin or inversion related points.
    real(dp), allocatable, intent(in) :: recPoint(:, :)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> element-specific electronegativity
    real(dp), intent(in) :: chi(:)

    !> element-specific chemical hardnesses
    real(dp), intent(in) :: gam(:)

    !> element-specific CN scaling constant
    real(dp), intent(in) :: kcn(:)

    !> element-specific charge widths / atomic radii
    real(dp), intent(in) :: rad(:)

    !> Error function coordination number.
    real(dp), intent(in) :: cn(:)

    !> Derivative of the CN with respect to the Cartesian coordinates.
    real(dp), intent(in) :: dcndr(:, :, :)

    !> Derivative of the CN with respect to strain deformations.
    real(dp), intent(in) :: dcndL(:, :, :)

    !> Updated energy vector at return
    real(dp), intent(inout), optional :: energies(:)

    !> Updated gradient vector at return
    real(dp), intent(inout), optional :: gradients(:, :)

    !> Updated stress tensor at return
    real(dp), intent(inout), optional :: stress(:, :)

    !> Atomic partial charges.
    real(dp), intent(out), optional :: qAtom(:)

    !> Derivative of partial charges w.r.t. cartesian coordinates.
    real(dp), intent(out), optional :: dqdr(:, :, :)

    !> Derivative of partial charges w.r.t. strain deformations.
    real(dp), intent(out), optional :: dqdL(:, :, :)

    real(dp), parameter :: small = 1.0e-14_dp

    !> Matrix of 1/R values for each atom pair.
    real(dp), allocatable :: aMat(:, :)
    real(dp), allocatable :: xVec(:)
    real(dp), allocatable :: xFac(:)
    real(dp), allocatable :: qVec(:)

    real(dp), allocatable :: aInv(:, :)

    real(dp), allocatable :: aMatdr(:, :, :)
    real(dp), allocatable :: aMatdL(:, :, :)
    real(dp), allocatable :: aTrace(:, :)

    logical :: tPeriodic
    integer :: nDim
    integer :: iAt1, iSp1, ii
    real(dp) :: tmp

    tPeriodic = allocated(recPoint)

    nDim = nAtom + 1
    allocate(xVec(nDim), xFac(nAtom), aMatdr(3, nAtom, nDim), aMat(nDim, nDim), &
        & aInv(nDim, nDim), qVec(nDim), aMatdL(3, 3, nDim), aTrace(3, nAtom))

    ! Step 1: contruct RHS for linear system
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      tmp = kcn(iSp1) / (sqrt(cn(iAt1)) + small)
      xVec(iAt1) = -chi(iSp1) + tmp*cn(iAt1)
      xFac(iAt1) = -0.5_dp*tmp
    end do
    xVec(nDim) = charge

    ! Step 2: construct interaction matrix
    if (tPeriodic) then
      call getCoulombMatrixPeriodic(nAtom, coords, species, nNeighbour, iNeighbour, neighDist2,&
          & img2CentCell, recPoint, gam, rad, alpha, volume, aMat)
    else
      call getCoulombMatrixCluster(nAtom, coords, species, gam, rad, aMat)
    end if

    aMat(nDim, 1:nAtom) = 1.0_dp
    aMat(1:nAtom, nDim) = 1.0_dp
    aMat(nDim, nDim) = 0.0_dp

    ! Step 3: invert linear system
    aInv(:, :) = aMat
    call symmatinv(aInv)
    qVec(:) = 0.0_dp
    call hemv(qVec, aInv, xVec)

    ! Step 4: return atom resolved energies if requested, xVec is scratch
    if (present(energies)) then
      call hemv(xVec, aMat, qVec, alpha=0.5_dp, beta=-1.0_dp)
      energies(:) = energies(:) + xVec(:nAtom) * qVec(:nAtom)
    end if
    deallocate(xVec) ! free xVec to avoid using it later

    ! Step 5: get derivative of interaction matrix
    if (tPeriodic) then
      call getCoulombDerivsPeriodic(nAtom, coords, species, nNeighbour, iNeighbour, neighDist2,&
          & img2CentCell, recPoint, alpha, volume, rad, qVec, aMatdr, aMatdL, aTrace)
    else
      call getCoulombDerivsCluster(nAtom, coords, species, rad, qVec, aMatdr, aTrace)
    end if
    do iAt1 = 1, nAtom
      aMatdr(:, :, iAt1) = aMatdr(:, :, iAt1) + dcndr(:, :, iAt1) * xFac(iAt1)
      aMatdL(:, :, iAt1) = aMatdL(:, :, iAt1) + dcndL(:, :, iAt1) * xFac(iAt1)
    end do

    if (present(gradients)) then
      call gemv(gradients, aMatdr, qVec, beta=1.0_dp)
    end if

    if (present(stress)) then
      call gemv(stress, aMatdL, qVec, beta=1.0_dp)
    end if

    do iAt1 = 1, nAtom
      aMatdr(:, iAt1, iAt1) = aMatdr(:, iAt1, iAt1) + aTrace(:, iAt1)
    end do

    if (present(dqdr) .or. present(dqdL)) then
      ! Symmetrise inverted matrix
      do ii = 1, nDim
        aInv(ii, ii+1:nDim) = aInv(ii+1:nDim, ii)
      end do
    end if

    if (present(dqdr)) then
      dqdr(:, :, :) = 0.0_dp
      call gemm(dqdr, aMatdr, aInv, alpha=-1.0_dp)
    end if

    if (present(dqdL)) then
      dqdL(:, :, :) = 0.0_dp
      call gemm(dqdL, aMatdL, aInv, alpha=-1.0_dp)
    end if

    if (present(qAtom)) then
      qAtom(:) = qVec(:nAtom)
    end if

  end subroutine getEEQCharges


end module dftbp_encharges
