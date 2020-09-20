!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> A simple reimplementation of DFT-D3
module dftbp_simpledftd3
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_blasroutines, only : gemv
  use dftbp_constants, only : pi, symbolToNumber
  use dftbp_coordnumber, only : TCNCont, TCNInput, init
  use dftbp_dftd3param
  use dftbp_dftd4param, only : getSqrtZr4r2
  use dftbp_environment, only : TEnvironment
  use dftbp_dispiface, only : TDispersionIface
  use dftbp_periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use dftbp_schedule, only : distributeRangeInChunks, assembleChunks
  use dftbp_simplealgebra, only : determinant33, invert33
  implicit none
  private

  public :: TSimpleDftD3, TSimpleDftD3Input, init


  !> Method parameters for DFT-D3 calculation.
  type :: TDftD3Param

    !> Scaling parameter for dipole-dipole coefficients.
    real(dp) :: s6 = 1.0_dp

    !> Scaling parameter for dipole-quadrupole coefficients.
    real(dp) :: s8

    !> Scaling parameter for quadrupole-quadrupole coefficients.
    real(dp) :: s10 = 0.0_dp

    !> Scaling parameter for <r4>/<r2> expectation value based critical radii.
    real(dp) :: a1

    !> Constant offset of critical radii.
    real(dp) :: a2

    !> Exponent of for the zero-damping function used for non-addititive triple dipole
    !> contributions.
    real(dp) :: alpha = 14.0_dp

    !> Gaussian weighting factor for interpolation of dispersion coefficients.
    real(dp) :: weightingFactor = 4.0_dp

  end type TDftD3Param


  !> Input for the DFT-D3 model
  type, extends(TDftD3Param) :: TSimpleDftD3Input

    !> Cutoff radius for dispersion interactions.
    real(dp) :: cutoffInter = 64.0_dp

    !> Coordination number specific input
    type(TCNInput) :: cnInput

  end type TSimpleDftD3Input


  !> Internal state of the DFT-D3 dispersion.
  type, extends(TDispersionIface) :: TSimpleDftD3
    private

    !> reference systems for dispersion
    type(TDftD3Ref) :: ref

    !> parameters for evaluating the dispersion energy
    type(TDftD3Param) :: param

    !> number of atoms
    integer :: nAtom

    !> energy
    real(dp), allocatable :: energies(:)

    !> force contributions
    real(dp), allocatable :: gradients(:,:)

    !> lattice vectors if periodic
    real(dp) :: latVecs(3, 3)

    !> Volume of the unit cell
    real(dp) :: vol

    !> Strain derivatives
    real(dp) :: sigma(3, 3)

    !> is this periodic
    logical :: tPeriodic

    !> are the coordinates current?
    logical :: tCoordsUpdated

    !> Coordination number
    type(TCNCont) :: cnCont

    !> Cutoff radius for dispersion interactions.
    real(dp) :: cutoffInter = 64.0_dp

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

  end type TSimpleDftD3


  !> Initialize DFT-D3 calculator from input
  interface init
    module procedure :: TSimpleDftD3_init
  end interface init


contains


  !> Initialize SimpleDftD3 instance.
  subroutine TSimpleDftD3_init(this, input, nAtom, species0, speciesNames, latVecs)

    !> Initialized instance of D4 dispersion model.
    type(TSimpleDftD3), intent(out) :: this

    !> Specific input parameters for damping function.
    type(TSimpleDftD3Input), intent(in) :: input

    !> Nr. of atoms in the system.
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell.
    integer, intent(in) :: species0(:)

    !> Names of species.
    character(*), intent(in) :: speciesNames(:)

    !> Lattice vectors, if the system is periodic.
    real(dp), intent(in), optional :: latVecs(:, :)

    real(dp) :: recVecs(3, 3), maxGEwald
    integer :: iAt1

    this%tPeriodic = present(latVecs)

    this%param = input%TDftD3Param
    this%cutoffInter = input%cutoffInter

    call init(this%ref, speciesNames)

    if (this%tPeriodic) then
      call init(this%cnCont, input%cnInput, nAtom, latVecs)
    else
      call init(this%cnCont, input%cnInput, nAtom)
    end if

    this%tCoordsUpdated = .false.

    if (this%tPeriodic) then
      call this%updateLatVecs(latVecs)
    end if

    this%nAtom = nAtom

    allocate(this%energies(nAtom))
    allocate(this%gradients(3, nAtom))

  end subroutine TSimpleDftD3_init


  !> Notifies the objects about changed coordinates.
  subroutine updateCoords(this, env, neigh, img2CentCell, coords, species0, stat)

    !> Instance of DFTD4 data
    class(TSimpleDftD3), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Updated neighbour list.
    type(TNeighbourList), intent(in) :: neigh

    !> Updated mapping to central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Updated coordinates.
    real(dp), intent(in) :: coords(:,:)

    !> Status of operation
    integer, intent(out), optional :: stat

    ! Species of the atoms in the unit cell.
    integer, intent(in) :: species0(:)

    ! Nr. of neighbours for each atom
    integer, allocatable :: nNeighbour(:)
    real(dp) :: sigma(3, 3)

    @:ASSERT(allocated(this%energies))
    @:ASSERT(allocated(this%gradients))

    if (present(stat)) then
      stat = 0
    end if

    allocate(nNeighbour(this%nAtom))

    this%energies(:) = 0.0_dp
    this%gradients(:, :) = 0.0_dp
    this%sigma(:, :) = 0.0_dp

    call this%cnCont%updateCoords(neigh, img2CentCell, coords, species0)

    call getNrOfNeighboursForAll(nNeighbour, neigh, this%cutoffInter)
    call dispersionGradient(env, this%ref, this%param, nNeighbour, coords, &
        & species0, neigh, img2CentCell, this%cnCont%cn, this%cnCont%dcndr, &
        & this%cnCont%dcndL, this%energies, this%gradients, this%sigma)

    this%tCoordsUpdated = .true.

  end subroutine updateCoords


  !> Notifies the object about updated lattice vectors.
  subroutine updateLatVecs(this, latVecs)

    !> Instance of DFTD4 data
    class(TSimpleDftD3), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:, :)

    real(dp) :: recVecs(3, 3), maxGEwald

    @:ASSERT(this%tPeriodic)
    @:ASSERT(all(shape(latvecs) == shape(this%latvecs)))

    this%latVecs = latVecs
    this%vol = abs(determinant33(latVecs))

    call this%cnCont%updateLatVecs(LatVecs)

    this%tCoordsUpdated = .false.

  end subroutine updateLatVecs


  !> Returns the atomic resolved energies due to the dispersion.
  subroutine getEnergies(this, energies)

    !> Instance of DFTD4 data
    class(TSimpleDftD3), intent(inout) :: this

    !> Contains the atomic energy contributions on exit.
    real(dp), intent(out) :: energies(:)

    @:ASSERT(allocated(this%energies))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(size(energies) == this%nAtom)

    energies(:) = this%energies

  end subroutine getEnergies


  !> Adds the atomic gradients to the provided vector.
  subroutine addGradients(this, gradients)

    !> Instance of DFTD4 data
    class(TSimpleDftD3), intent(inout) :: this

    !> The vector to increase by the gradients.
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(allocated(this%gradients))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    gradients(:,:) = gradients + this%gradients

  end subroutine addGradients


  !> Returns the stress tensor.
  subroutine getStress(this, stress)

    !> Instance of DFTD4 data
    class(TSimpleDftD3), intent(inout) :: this

    !> stress tensor from the dispersion
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(this%tPeriodic)
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))

    stress(:,:) = this%sigma / this%vol

  end subroutine getStress


  !> Estimates the real space cutoff of the dispersion interaction.
  function getRCutoff(this) result(cutoff)

    !> Instance of DFTD4 data
    class(TSimpleDftD3), intent(inout) :: this

    !> Resulting cutoff
    real(dp) :: cutoff

    cutoff = max(this%cutoffInter, this%cnCont%getRCutoff())

  end function getRCutoff


  !> Calculate the weights of the reference systems.
  !> It also saves the derivative w.r.t. coordination number for later use.
  subroutine weightReferences(env, calc, param, nAtom, species, cn, gwVec, gwdcn)

    !> Computation environment
    type(TEnvironment), intent(in) :: env

    !> DFT-D dispersion model.
    type(TDftD3Ref), intent(in) :: calc

    !> DFT-D dispersion model.
    type(TDftD3Param), intent(in) :: param

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> Coordination number of every atom.
    real(dp), intent(in) :: cn(:)

    !> weighting function for the atomic reference systems
    real(dp), intent(out) :: gwVec(:, :)

    !> derivative of the weighting function w.r.t. the coordination number
    real(dp), intent(out) :: gwdcn(:, :)

    integer :: iAtFirst, iAtLast, iAt1, iSp1, iRef1
    real(dp) :: norm, dnorm, wf, gw, expw, expd, gwk, dgwk

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    gwVec(:,:) = 0.0_dp
    gwdcn(:,:) = 0.0_dp

    !$omp parallel do schedule(runtime) default(none) &
    !$omp shared(gwVec, gwdcn, calc, param, iAtFirst, iAtLast, species, cn) &
    !$omp private(iAt1, iSp1, norm, dnorm, iRef1, wf, gw, expw, expd, gwk, dgwk)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      norm = 0.0_dp
      dnorm = 0.0_dp
      do iRef1 = 1, calc%numberOfReferences(iSp1)
        wf = param%weightingFactor
        gw = weightCN(wf, cn(iAt1), calc%referenceCN(iRef1, iSp1))
        norm = norm + gw
        dnorm = dnorm + 2*wf*(calc%referenceCN(iRef1, iSp1) - cn(iAt1)) * gw
      end do
      norm = 1.0_dp / norm
      do iRef1 = 1, calc%numberOfReferences(iSp1)
        wf = param%weightingFactor
        gw = weightCN(wf, cn(iAt1), calc%referenceCN(iRef1, iSp1))
        expw = gw
        expd = 2.0_dp * wf * (calc%referenceCN(iRef1, iSp1) - cn(iAt1)) * gw

        gwk = expw * norm
        if (ieee_is_nan(gwk)) then
          if (maxval(calc%referenceCN(:calc%numberOfReferences(iSp1), iSp1))&
              & == calc%referenceCN(iRef1, iSp1)) then
            gwk = 1.0_dp
          else
            gwk = 0.0_dp
          end if
        end if
        gwVec(iRef1, iAt1) = gwk

        dgwk = expd * norm - expw * dnorm * norm**2
        if (ieee_is_nan(dgwk)) then
          dgwk = 0.0_dp
        end if
        gwdcn(iRef1, iAt1) = dgwk
      end do
    end do
    !$omp end parallel do

    call assembleChunks(env, gwVec)
    call assembleChunks(env, gwdcn)

  end subroutine weightReferences


  !> Actual implementation of the dispersion energy, gradients and stress tensor.
  subroutine dispersionGradient(env, calc, param, nNeighbour, coords, species, &
      & neigh, img2CentCell, cn, dcndr, dcndL, energies, gradients, sigma)

    !> Computation environment
    type(TEnvironment), intent(in) :: env

    !> DFT-D dispersion model.
    type(TDftD3Ref), intent(in) :: calc

    !> DFT-D dispersion model.
    type(TDftD3Param), intent(in) :: param

    !> Nr. of neighbours for every atom
    integer, intent(in) :: nNeighbour(:)

    !> Coordinates of the atoms (including images)
    real(dp), intent(in) :: coords(:, :)

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> Updated neighbour list.
    type(TNeighbourList), intent(in) :: neigh

    !> Mapping into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Coordination number of every atom.
    real(dp), intent(in) :: cn(:)

    !> Derivative of coordination number w.r.t. the cartesian coordinates.
    real(dp), intent(in) :: dcndr(:, :, :)

    !> Derivative of coordination number w.r.t. the strain deformations.
    real(dp), intent(in) :: dcndL(:, :, :)

    !> Updated energy vector at return
    real(dp), intent(inout) :: energies(:)

    !> Updated gradient vector at return
    real(dp), intent(inout) :: gradients(:, :)

    !> Updated sigma tensor at return
    real(dp), intent(inout) :: sigma(:, :)

    integer :: nRef, nAtom, iAtFirst, iAtLast
    integer :: iAt1, iSp1, iNeigh, iAt2, iSp2, iAt2f
    real(dp) :: dc6, dc6dcn1, dc6dcn2
    real(dp) :: vec(3), grad(3), dEr, dGr, dSr(3, 3)
    real(dp) :: rc, r1, r2, r4, r5, r6, r8, r10, rc1, rc2, rc6, rc8, rc10
    real(dp) :: f6, df6, f8, df8, f10, df10

    real(dp), allocatable :: gwVec(:, :), gwdcn(:, :), dEdcn(:)
    real(dp), allocatable :: c6(:, :), dc6dcn(:, :)
    real(dp), allocatable :: localDeriv(:,:), localSigma(:, :), localEnergies(:)

    nAtom = size(nNeighbour)
    nRef = maxval(calc%numberOfReferences(species))
    allocate(gwVec(nRef, nAtom), gwdcn(nRef, nAtom), c6(nAtom, nAtom), dc6dcn(nAtom, nAtom),&
        & dEdcn(nAtom))

    dEdcn(:) = 0.0_dp

    call weightReferences(env, calc, param, nAtom, species, cn, gwVec, gwdcn)

    call getAtomicC6(env, calc, nAtom, species, gwVec, gwdcn, c6, dc6dcn)

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    allocate(localEnergies(nAtom), localDeriv(3, nAtom), localSigma(3, 3))
    localEnergies(:) = 0.0_dp
    localDeriv(:,:) = 0.0_dp
    localSigma(:,:) = 0.0_dp

    !$omp parallel do default(none) schedule(runtime) &
    !$omp reduction(+:localEnergies, localDeriv, localSigma, dEdcn) &
    !$omp shared(iAtFirst, iAtLast, species, nNeighbour, neigh, coords) &
    !$omp shared(img2CentCell, c6, dc6dcn, calc, param) &
    !$omp private(iAt1, iSp1, iNeigh, iAt2, vec, iAt2f, iSp2, r2, r1, r4, r5) &
    !$omp private(r6, r8, r10, rc1, rc2, rc6, rc8, rc10, dc6, dc6dcn1, dc6dcn2) &
    !$omp private(rc, dEr, dGr, grad, dSr, f6, f8, f10, df6, df8, df10)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeighbour(iAt1)
        iAt2 = neigh%iNeighbour(iNeigh, iAt1)
        vec(:) = coords(:, iAt1) - coords(:, iAt2)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        r2 = neigh%neighDist2(iNeigh, iAt1)
        r1 = sqrt(r2)
        r4 = r2 * r2
        r5 = r4 * r1
        r6 = r2 * r4
        r8 = r2 * r6
        r10 = r2 * r8

        rc = 3.0_dp * calc%sqrtZr4r2(iSp1) * calc%sqrtZr4r2(iSp2)
        rc1 = param%a1 * sqrt(rc) + param%a2
        rc2 = rc1 * rc1
        rc6 = rc2 * rc2 * rc2
        rc8 = rc2 * rc6
        rc10 = rc2 * rc8

        dc6 = c6(iAt1, iAt2f)
        dc6dcn1 = dc6dcn(iAt1, iAt2f)
        dc6dcn2 = dc6dcn(iAt2f, iAt1)

        f6 = 1.0_dp / (r6 + rc6)
        f8 = 1.0_dp / (r8 + rc8)
        f10 = 1.0_dp / (r10 + rc10)

        df6 = -6.0_dp * r5 * f6 * f6
        df8 = -8.0_dp * r2 * r5 * f8 * f8
        df10 = -10.0_dp * r4 * r5 * f10 * f10

        dEr = param%s6 * f6 + param%s8 * f8 * rc + param%s10 * rc * rc * 49.0_dp / 40.0_dp * f10
        dGr = param%s6 * df6 + param%s8 * df8 * rc + param%s10 * rc * rc * 49.0_dp / 40.0_dp * df10

        grad(:) = -dGr*dc6 * vec / r1
        dSr(:,:) = spread(grad, 1, 3) * spread(vec, 2, 3)

        localEnergies(iAt1) = localEnergies(iAt1) - dEr*dc6/2
        dEdcn(iAt1) = dEdcn(iAt1) - dc6dcn1 * dEr
        if (iAt1 /= iAt2f) then
          localSigma(:,:) = localSigma - dSr
          localEnergies(iAt2f) = localEnergies(iAt2f) - dEr*dc6/2
          localDeriv(:, iAt1) = localDeriv(:, iAt1) + grad
          localDeriv(:, iAt2f) = localDeriv(:, iAt2f) - grad
          dEdcn(iAt2f) = dEdcn(iAt2f) - dc6dcn2 * dEr
        else
          localSigma(:,:) = localSigma - 0.5_dp * dSr
        end if

      end do
    end do
    !$omp end parallel do

    call assembleChunks(env, localEnergies)
    call assembleChunks(env, localDeriv)
    call assembleChunks(env, localSigma)
    call assembleChunks(env, dEdcn)

    energies(:) = energies + localEnergies
    gradients(:,:) = gradients + localDeriv
    sigma(:,:) = sigma + localSigma

    ! handle CN contributions to the gradient by matrix-vector operation
    call gemv(gradients, dcndr, dEdcn, beta=1.0_dp)

    ! handle CN contributions to the sigma tensor
    call gemv(sigma, dcndL, dEdcn, beta=1.0_dp)

  end subroutine dispersionGradient


  !> calculate atomic dispersion coefficients and their derivatives w.r.t.
  !> coordination number and partial charge.
  subroutine getAtomicC6(env, calc, nAtom, species, gwVec, gwdcn, c6, dc6dcn)

    !> Computation environment
    type(TEnvironment), intent(in) :: env

    !> DFT-D dispersion model.
    type(TDftD3Ref), intent(in) :: calc

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> weighting function for the atomic reference systems
    real(dp), intent(in) :: gwVec(:, :)

    !> derivative of the weighting function w.r.t. the coordination number
    real(dp), intent(in) :: gwdcn(:, :)

    !> C6 coefficients for all atom pairs.
    real(dp), intent(out) :: c6(:, :)

    !> derivative of the C6 w.r.t. the coordination number
    real(dp), intent(out) :: dc6dcn(:, :)

    integer :: iAtFirst, iAtLast, iAt1, iAt2, iSp1, iSp2, iRef1, iRef2
    real(dp) :: refc6, dc6, dc6dcn1, dc6dcn2

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    c6(:,:) = 0.0_dp
    dc6dcn(:,:) = 0.0_dp

    !$omp parallel do default(none) schedule(runtime) shared(c6, dc6dcn) &
    !$omp shared(calc, gwVec, gwdcn, iAtFirst, iAtLast, species) &
    !$omp private(iAt1, iAt2, iSp1, iSp2, iRef1, iRef2, refc6, dc6, dc6dcn1, dc6dcn2)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      do iAt2 = 1, iAt1
        iSp2 = species(iAt2)
        dc6 = 0.0_dp
        dc6dcn1 = 0.0_dp
        dc6dcn2 = 0.0_dp
        do iRef1 = 1, calc%numberOfReferences(iSp1)
          do iRef2 = 1, calc%numberOfReferences(iSp2)
            refc6 = calc%referenceC6(iRef1, iRef2, iSp1, iSp2)
            dc6 = dc6 + gwVec(iRef1, iAt1) * gwVec(iRef2, iAt2) * refc6
            dc6dcn1 = dc6dcn1 + gwdcn(iRef1, iAt1) * gwVec(iRef2, iAt2) * refc6
            dc6dcn2 = dc6dcn2 + gwVec(iRef1, iAt1) * gwdcn(iRef2, iAt2) * refc6
          end do
        end do
        c6(iAt1, iAt2) = dc6
        c6(iAt2, iAt1) = dc6
        dc6dcn(iAt1, iAt2) = dc6dcn1
        dc6dcn(iAt2, iAt1) = dc6dcn2
      end do
    end do
    !$omp end parallel do

    call assembleChunks(env, c6)
    call assembleChunks(env, dc6dcn)

  end subroutine getAtomicC6


  !> Gaussian weight based on coordination number difference
  elemental function weightCN(wf, cn, cnref) result(cngw)

    !> weighting factor / width of the gaussian function
    real(dp), intent(in) :: wf

    !> current coordination number
    real(dp), intent(in) :: cn

    !> reference coordination number
    real(dp), intent(in) :: cnref

    !> CN-gaussian-weight
    real(dp) :: cngw

    cngw = exp(-wf * (cn - cnref)**2)

  end function weightCN


end module dftbp_simpledftd3
