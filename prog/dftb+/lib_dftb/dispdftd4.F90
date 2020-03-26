!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> implementation of the D4 dispersion model
module dftbp_dispdftd4
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_environment, only : TEnvironment
  use dftbp_blasroutines, only : gemv
  use dftbp_constants, only : pi, symbolToNumber
  use dftbp_coordnumber, only : TCNCont, init
  use dftbp_dftd4param, only : TDftD4Calculator, TDispDftD4Inp, initializeCalculator
  use dftbp_dispiface, only : TDispersionIface
  use dftbp_encharges, only : TEeqCont, init
  use dftbp_periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use dftbp_simplealgebra, only : determinant33
  implicit none
  private

  public :: TDispDftD4, TDispDftD4Inp, init


  !> Internal state of the DFT-D4 dispersion.
  type, extends(TDispersionIface) :: TDispDftD4
    private

    !> calculator to evaluate dispersion
    type(TDftD4Calculator), allocatable :: calculator

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

    !> stress tensor
    real(dp) :: stress(3, 3)

    !> is this periodic
    logical :: tPeriodic

    !> are the coordinates current?
    logical :: tCoordsUpdated

    !> EEQ model
    type(TEeqCont) :: eeqCont

    !> Coordination number
    type(TCNCont) :: cnCont

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

  end type TDispDftD4


  interface init
    module procedure :: DispDftD4_init
  end interface init


contains


  !> Initialize DispDftD4 instance.
  subroutine DispDftD4_init(this, inp, nAtom, species0, speciesNames, latVecs)

    !> Initialized instance of D4 dispersion model.
    type(TDispDftD4), intent(out) :: this

    !> Specific input parameters for damping function.
    type(TDispDftD4Inp), intent(in) :: inp

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

    if (this%tPeriodic) then
      call init(this%eeqCont, inp%eeqInput, .false., .true., nAtom, latVecs)
      call init(this%cnCont, inp%cnInput, nAtom, latVecs)
    else
      call init(this%eeqCont, inp%eeqInput, .false., .true., nAtom)
      call init(this%cnCont, inp%cnInput, nAtom)
    end if

    this%tCoordsUpdated = .false.

    if (this%tPeriodic) then
      call this%updateLatVecs(latVecs)
    end if

    this%nAtom = nAtom

    allocate(this%energies(nAtom))
    allocate(this%gradients(3, nAtom))

    allocate(this%calculator)
    call initializeCalculator(this%calculator, inp, this%nAtom, speciesNames)

  end subroutine DispDftD4_init


  !> Notifies the objects about changed coordinates.
  subroutine updateCoords(this, env, neigh, img2CentCell, coords, species0)

    !> Instance of DFTD4 data
    class(TDispDftD4), intent(inout) :: this

    !> Updated neighbour list.
    type(TNeighbourList), intent(in) :: neigh

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Updated mapping to central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Updated coordinates.
    real(dp), intent(in) :: coords(:,:)

    !> Species of the atoms in the unit cell.
    integer, intent(in) :: species0(:)

    @:ASSERT(allocated(this%calculator))

    if (this%tPeriodic) then
      call dispersionEnergy(this%calculator, this%nAtom, coords, species0, neigh, img2CentCell,&
          & this%eeqCont, this%cnCont, this%energies, this%gradients, &
          & stress=this%stress, volume=this%vol, parEwald=this%eeqCont%parEwald)
    else
      call dispersionEnergy(this%calculator, this%nAtom, coords, species0, neigh, img2CentCell,&
          & this%eeqCont, this%cnCont, this%energies, this%gradients)
    end if

    this%tCoordsUpdated = .true.

  end subroutine updateCoords


  !> Notifies the object about updated lattice vectors.
  subroutine updateLatVecs(this, latVecs)

    !> Instance of DFTD4 data
    class(TDispDftD4), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:, :)

    real(dp) :: recVecs(3, 3), maxGEwald

    @:ASSERT(this%tPeriodic)
    @:ASSERT(all(shape(latvecs) == shape(this%latvecs)))

    this%latVecs = latVecs
    this%vol = abs(determinant33(latVecs))

    call this%eeqCont%updateLatVecs(latVecs)
    call this%cnCont%updateLatVecs(LatVecs)

    this%tCoordsUpdated = .false.

  end subroutine updateLatVecs


  !> Returns the atomic resolved energies due to the dispersion.
  subroutine getEnergies(this, energies)

    !> Instance of DFTD4 data
    class(TDispDftD4), intent(inout) :: this

    !> Contains the atomic energy contributions on exit.
    real(dp), intent(out) :: energies(:)

    @:ASSERT(allocated(this%calculator))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(size(energies) == this%nAtom)

    energies(:) = this%energies

  end subroutine getEnergies


  !> Adds the atomic gradients to the provided vector.
  subroutine addGradients(this, gradients)

    !> Instance of DFTD4 data
    class(TDispDftD4), intent(inout) :: this

    !> The vector to increase by the gradients.
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(allocated(this%calculator))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    gradients(:,:) = gradients + this%gradients

  end subroutine addGradients


  !> Returns the stress tensor.
  subroutine getStress(this, stress)

    !> Instance of DFTD4 data
    class(TDispDftD4), intent(inout) :: this

    !> stress tensor from the dispersion
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(allocated(this%calculator))
    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))

    stress(:,:) = this%stress

  end subroutine getStress


  !> Estimates the real space cutoff of the dispersion interaction.
  function getRCutoff(this) result(cutoff)

    !> Instance of DFTD4 data
    class(TDispDftD4), intent(inout) :: this

    !> Resulting cutoff
    real(dp) :: cutoff

    cutoff = max(this%calculator%cutoffInter, this%cnCont%getRCutoff(),&
        & this%eeqCont%getRCutoff(), this%calculator%cutoffThree)

  end function getRCutoff


  !> Calculate the weights of the reference systems and scale them according
  !> to their partial charge difference. It also saves the derivatives w.r.t.
  !> coordination number and partial charge for later use.
  !>
  !> This subroutine produces two sets of weighting and scaling vectors, one
  !> for the input partial charges and one for the case of no charge fluctuations,
  !> meaning q=0. The normal zeta-Vector is used to scale the pair-wise dispersion
  !> terms while the so called zero-Vector is used for the non-additive terms in
  !> the ATM dispersion energy.
  subroutine weightReferences(calc, nAtom, species, cn, q, zetaVec, zeroVec, zetadq, zetadcn,&
      & zerodcn)

    !> DFT-D dispersion model.
    type(TDftD4Calculator), intent(in) :: calc

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> Coordination number of every atom.
    real(dp), intent(in) :: cn(:)

    !> Partial charge of every atom.
    real(dp), intent(in) :: q(:)

    !> weighting and scaling function for the atomic reference systems
    real(dp), intent(out) :: zetaVec(:, :)

    !> weighting and scaling function for the atomic reference systems for q=0
    real(dp), intent(out) :: zeroVec(:, :)

    !> derivative of the weight'n'scale function w.r.t. the partial charges
    real(dp), intent(out) :: zetadq(:, :)

    !> derivative of the weight'n'scale function w.r.t. the coordination number
    real(dp), intent(out) :: zetadcn(:, :)

    !> derivative of the weight'n'scale function w.r.t. the CN for q=0
    real(dp), intent(out) :: zerodcn(:, :)

    integer :: iAt1, iSp1, iRef1, iCount1
    real(dp) :: eta1, zEff1, qRef1
    real(dp) :: norm, dnorm, wf, gw, expw, expd, gwk, dgwk

    zetaVec(:,:) = 0.0_dp
    zetadq(:,:) = 0.0_dp
    zetadcn(:,:) = 0.0_dp
    zeroVec(:,:) = 0.0_dp
    zerodcn(:,:) = 0.0_dp

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP SHARED(zetaVec, zetadq, zetadcn, zeroVec, zerodcn, calc, nAtom, species, q, cn) &
    !$OMP PRIVATE(iAt1, iSp1, eta1, zEff1, norm, dnorm, iRef1, iCount1, wf, gw) &
    !$OMP PRIVATE(expw, expd, gwk, dgwk, qRef1)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      eta1 = calc%gc * calc%chemicalHardness(iSp1)
      zEff1 = calc%effectiveNuclearCharge(iSp1)
      norm = 0.0_dp
      dnorm = 0.0_dp
      do iRef1 = 1, calc%numberOfReferences(iSp1)
        do iCount1 = 1, calc%countNumber(iRef1, iSp1)
          wf = iCount1 * calc%wf
          gw = weightCN(wf, cn(iAt1), calc%referenceCN(iRef1, iSp1))
          norm = norm + gw
          dnorm = dnorm + 2*wf*(calc%referenceCN(iRef1, iSp1) - cn(iAt1)) * gw
        end do
      end do
      norm = 1.0_dp / norm
      do iRef1 = 1, calc%numberOfReferences(iSp1)
        expw = 0.0_dp
        expd = 0.0_dp
        do iCount1 = 1, calc%countNumber(iref1, iSp1)
          wf = iCount1 * calc%wf
          gw = weightCN(wf, cn(iAt1), calc%referenceCN(iRef1, iSp1))
          expw = expw + gw
          expd = expd + 2.0_dp * wf * (calc%referenceCN(iRef1, iSp1) - cn(iAt1)) * gw
        end do

        gwk = expw * norm
        if (ieee_is_nan(gwk)) then
          if (maxval(calc%referenceCN(:calc%numberOfReferences(iSp1), iSp1))&
              & == calc%referenceCN(iRef1, iSp1)) then
            gwk = 1.0_dp
          else
            gwk = 0.0_dp
          end if
        end if
        qRef1 = calc%referenceCharge(iRef1, iSp1) + zEff1
        zetaVec(iRef1, iAt1) = zetaScale(calc%ga, eta1, qRef1, q(iAt1) + zEff1) * gwk
        zeroVec(iRef1, iAt1) = zetaScale(calc%ga, eta1, qRef1, zEff1) * gwk

        zetadq(iRef1, iAt1) = dzetaScale(calc%ga, eta1, qRef1, q(iAt1)+zEff1) * gwk

        dgwk = expd * norm - expw * dnorm * norm**2
        if (ieee_is_nan(dgwk)) then
          dgwk = 0.0_dp
        end if
        zetadcn(iRef1, iAt1) = zetaScale(calc%ga, eta1, qRef1, q(iAt1) + zEff1) * dgwk
        zerodcn(iRef1, iAt1) = zetaScale(calc%ga, eta1, qRef1, zEff1) * dgwk
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine weightReferences


  !> Actual implementation of the dispersion energy, gradients and stress tensor.
  !>
  !> This subroutine is somewhat agnostic to the dispersion model, meaning that
  !> it can be used for either D4 or D3(BJ), if the inputs like coordination
  !> number, charges and scaling parameters are chosen accordingly.
  subroutine dispersionGradient(calc, nAtom, coords, species, neigh, img2CentCell, cn, dcndr,&
      & dcndL, q, dqdr, dqdL, energies, gradients, stress)

    !> DFT-D dispersion model.
    type(TDftD4Calculator), intent(in) :: calc

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

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

    !> Partial charges of every atom.
    real(dp), intent(in) :: q(:)

    !> Derivative of partial charges w.r.t. the cartesian coordinates.
    real(dp), intent(in) :: dqdr(:, :, :)

    !> Derivative of partial charges w.r.t. the strain deformations.
    real(dp), intent(in) :: dqdL(:, :, :)

    !> Updated energy vector at return
    real(dp), intent(inout) :: energies(:)

    !> Updated gradient vector at return
    real(dp), intent(inout) :: gradients(:, :)

    !> Updated stress tensor at return
    real(dp), intent(inout) :: stress(:, :)

    integer :: nRef
    integer :: iAt1, iSp1, iNeigh, iAt2, iSp2, iAt2f
    real(dp) :: dc6, dc6dcn1, dc6dcn2, dc6dq1, dc6dq2
    real(dp) :: vec(3), grad(3), dEr, dGr, sigma(3, 3)
    real(dp) :: rc, r1, r2, r4, r5, r6, r8, r10, rc1, rc2, rc6, rc8, rc10
    real(dp) :: f6, df6, f8, df8, f10, df10

    !> Nr. of neighbours for each atom
    integer, allocatable :: nNeighbour(:)

    real(dp), allocatable :: zetaVec(:, :), zeroVec(:, :)
    real(dp), allocatable :: zetadq(:, :), zerodq(:, :)
    real(dp), allocatable :: zetadcn(:, :), zerodcn(:, :)
    real(dp), allocatable :: dEdq(:), dEdcn(:)
    real(dp), allocatable :: c6(:, :), dc6dq(:, :), dc6dcn(:, :)

    nRef = maxval(calc%numberOfReferences(species))
    allocate(nNeighbour(nAtom))
    allocate(zetaVec(nRef, nAtom), zeroVec(nRef, nAtom), zetadq(nRef, nAtom),&
        & zetadcn(nRef, nAtom), zerodcn(nRef, nAtom), zerodq(nRef, nAtom),&
        & c6(nAtom, nAtom), dc6dq(nAtom, nAtom), dc6dcn(nAtom, nAtom),&
        & dEdq(nAtom), dEdcn(nAtom))

    dEdq(:) = 0.0_dp
    dEdcn(:) = 0.0_dp

    call weightReferences(calc, nAtom, species, cn, q, zetaVec, zeroVec, zetadq,&
        & zetadcn, zerodcn)

    call getAtomicC6(calc, nAtom, species, zetaVec, zetadq, zetadcn,&
        & c6, dc6dcn, dc6dq)

    call getNrOfNeighboursForAll(nNeighbour, neigh, calc%cutoffInter)
    !$OMP PARALLEL DEFAULT(NONE) REDUCTION(+:energies, gradients, stress, dEdq, dEdcn) &
    !$OMP SHARED(nAtom, species, nNeighbour, neigh, coords, img2CentCell, c6, dc6dq, dc6dcn, calc) &
    !$OMP PRIVATE(iAt1, iSp1, iNeigh, iAt2, vec, iAt2f, iSp2, r2, r1, r4, r5, r6, r8, r10) &
    !$OMP PRIVATE(rc, rc1, rc2, rc6, rc8, rc10, dc6, dc6dcn1, dc6dcn2, dc6dq1, dc6dq2) &
    !$OMP PRIVATE(dEr, dGr, grad, sigma, f6, f8, f10, df6, df8, df10)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
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
        rc1 = calc%a1 * sqrt(rc) + calc%a2
        rc2 = rc1 * rc1
        rc6 = rc2 * rc2 * rc2
        rc8 = rc2 * rc6
        rc10 = rc2 * rc8

        dc6 = c6(iAt1, iAt2f)
        dc6dcn1 = dc6dcn(iAt1, iAt2f)
        dc6dcn2 = dc6dcn(iAt2f, iAt1)
        dc6dq1 = dc6dq(iAt1, iAt2f)
        dc6dq2 = dc6dq(iAt2f, iAt1)

        f6 = 1.0_dp / (r6 + rc6)
        f8 = 1.0_dp / (r8 + rc8)
        f10 = 1.0_dp / (r10 + rc10)

        df6 = -6.0_dp * r5 * f6 * f6
        df8 = -8.0_dp * r2 * r5 * f8 * f8
        df10 = -10.0_dp * r4 * r5 * f10 * f10

        dEr = calc%s6 * f6 + calc%s8 * f8 * rc + calc%s10 * rc * rc * 49.0_dp / 40.0_dp * f10
        dGr = calc%s6 * df6 + calc%s8 * df8 * rc + calc%s10 * rc * rc * 49.0_dp / 40.0_dp * df10

        grad(:) = -dGr*dc6 * vec / r1
        sigma(:,:) = spread(grad, 1, 3) * spread(vec, 2, 3)

        energies(iAt1) = energies(iAt1) - dEr*dc6/2
        dEdcn(iAt1) = dEdcn(iAt1) - dc6dcn1 * dEr
        dEdq(iAt1) = dEdq(iAt1) - dc6dq1 * dEr
        if (iAt1 /= iAt2f) then
          stress(:,:) = stress - sigma
          energies(iAt2f) = energies(iAt2f) - dEr*dc6/2
          gradients(:, iAt1) = gradients(:, iAt1) + grad
          gradients(:, iAt2f) = gradients(:, iAt2f) - grad
          dEdcn(iAt2f) = dEdcn(iAt2f) - dc6dcn2 * dEr
          dEdq(iAt2f) = dEdq(iAt2f) - dc6dq2 * dEr
        else
          stress(:,:) = stress - 0.5_dp * sigma
        end if

      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    if (calc%s9 > 0.0_dp) then
      zerodq(:, :) = 0.0_dp  ! really make sure there is no q `dependency' from alloc
      call getNrOfNeighboursForAll(nNeighbour, neigh, calc%cutoffThree)
      call threeBodyDispersionGradient(calc, nAtom, coords, species, nNeighbour, neigh%iNeighbour,&
          & neigh%neighDist2, img2CentCell, zeroVec, zerodq, zerodcn, dEdq, dEdcn, energies,&
          & gradients, stress)
    end if

    ! handle CN and charge contributions to the gradient by matrix-vector operation
    call gemv(gradients, dcndr, dEdcn, beta=1.0_dp)
    call gemv(gradients, dqdr, dEdq, beta=1.0_dp)

    ! handle CN and charge contributions to the stress tensor
    call gemv(stress, dcndL, dEdcn, beta=1.0_dp)
    call gemv(stress, dqdL, dEdq, beta=1.0_dp)

  end subroutine dispersionGradient


  !> calculate atomic dispersion coefficients and their derivatives w.r.t.
  !> coordination number and partial charge.
  subroutine getAtomicC6(calc, nAtom, species, zetaVec, zetadq, zetadcn, c6, dc6dcn, dc6dq)

    !> DFT-D dispersion model.
    type(TDftD4Calculator), intent(in) :: calc

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> weighting and scaling function for the atomic reference systems
    real(dp), intent(in) :: zetaVec(:, :)

    !> derivative of the weight'n'scale function w.r.t. the partial charges
    real(dp), intent(in) :: zetadq(:, :)

    !> derivative of the weight'n'scale function w.r.t. the coordination number
    real(dp), intent(in) :: zetadcn(:, :)

    !> C6 coefficients for all atom pairs.
    real(dp), intent(out) :: c6(:, :)

    !> derivative of the C6 w.r.t. the partial charges
    real(dp), intent(out) :: dc6dq(:, :)

    !> derivative of the C6 w.r.t. the coordination number
    real(dp), intent(out) :: dc6dcn(:, :)

    integer :: iAt1, iAt2, iSp1, iSp2, iRef1, iRef2
    real(dp) :: refc6, dc6, dc6dcn1, dc6dcn2, dc6dq1, dc6dq2

    c6(:,:) = 0.0_dp
    dc6dcn(:,:) = 0.0_dp
    dc6dq(:,:) = 0.0_dp

    !$OMP PARALLEL DEFAULT(NONE) SHARED(c6, dc6dcn, dc6dq) &
    !$OMP SHARED(calc, zetaVec, zetadq, zetadcn, nAtom, species) &
    !$OMP PRIVATE(iAt1, iAt2, iSp1, iSp2, iRef1, iRef2, refc6) &
    !$OMP PRIVATE(dc6, dc6dcn1, dc6dcn2, dc6dq1, dc6dq2)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      do iAt2 = 1, iAt1
        iSp2 = species(iAt2)
        dc6 = 0.0_dp
        dc6dcn1 = 0.0_dp
        dc6dcn2 = 0.0_dp
        dc6dq1 = 0.0_dp
        dc6dq2 = 0.0_dp
        do iRef1 = 1, calc%numberOfReferences(iSp1)
          do iRef2 = 1, calc%numberOfReferences(iSp2)
            refc6 = calc%referenceC6(iRef1, iRef2, iSp1, iSp2)
            dc6 = dc6 + zetaVec(iRef1, iAt1) * zetaVec(iRef2, iAt2) * refc6
            dc6dcn1 = dc6dcn1 + zetadcn(iRef1, iAt1) * zetaVec(iRef2, iAt2) * refc6
            dc6dcn2 = dc6dcn2 + zetaVec(iRef1, iAt1) * zetadcn(iRef2, iAt2) * refc6
            dc6dq1 = dc6dq1 + zetadq(iRef1, iAt1) * zetaVec(iRef2, iAt2) * refc6
            dc6dq2 = dc6dq2 + zetaVec(iRef1, iAt1) * zetadq(iRef2, iAt2) * refc6
          end do
        end do
        c6(iAt1, iAt2) = dc6
        c6(iAt2, iAt1) = dc6
        dc6dcn(iAt1, iAt2) = dc6dcn1
        dc6dcn(iAt2, iAt1) = dc6dcn2
        dc6dq(iAt1, iAt2) = dc6dq1
        dc6dq(iAt2, iAt1) = dc6dq2
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine getAtomicC6


  !> Energies, derivatives and strain derivatives of the Axilrod-Teller-Muto
  !> non-additive triple-dipole contribution.
  subroutine threeBodyDispersionGradient(calc, nAtom, coords, species, nNeighbour, iNeighbour,&
      & neighDist2, img2CentCell, zetaVec, zetadq, zetadcn, dEdq, dEdcn, energies, gradients,&
      & stress)

    !> DFT-D dispersion model.
    type(TDftD4Calculator), intent(in) :: calc

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Coordinates of the atoms (including images)
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

    !> weighting and scaling function for the atomic reference systems
    real(dp), intent(in) :: zetaVec(:, :)

    !> derivative of the weight'n'scale function w.r.t. the partial charges
    real(dp), intent(in) :: zetadq(:, :)

    !> derivative of the weight'n'scale function w.r.t. the coordination number
    real(dp), intent(in) :: zetadcn(:, :)

    !> derivative of the energy w.r.t. the partial charges
    real(dp), intent(inout) :: dEdq(:)

    !> derivative of the energy w.r.t. the coordination number
    real(dp), intent(inout) :: dEdcn(:)

    !> Updated energy vector at return (misses CN and q contributions)
    real(dp), intent(inout) :: energies(:)

    !> Updated gradient vector at return (misses CN and q contributions)
    real(dp), intent(inout) :: gradients(:, :)

    !> Updated stress tensor at return (misses CN and q contributions)
    real(dp), intent(inout) :: stress(:, :)

    integer :: iAt1, iAt2, iAt3, iAt2f, iAt3f, iSp1, iSp2, iSp3
    integer :: iNeigh2, iNeigh3
    real(dp) :: vec12(3), vec13(3), vec23(3), dist12, dist13, dist23
    real(dp) :: c9, c6_12, c6_13, c6_23, rc12, rc13, rc23, rc
    real(dp) :: r1, r2, r3, r5, rr, fdmp, dfdmp, ang, dang
    real(dp) :: dEr, dG12(3), dG13(3), dG23(3), sigma(3, 3)
    real(dp) :: dc9dcn1, dc9dcn2, dc9dcn3, dc9dq1, dc9dq2, dc9dq3
    real(dp), allocatable :: c6(:, :), dc6dq(:, :), dc6dcn(:, :)

    allocate(c6(nAtom, nAtom), dc6dq(nAtom, nAtom), dc6dcn(nAtom, nAtom))
    call getAtomicC6(calc, nAtom, species, zetaVec, zetadq, zetadcn, c6, dc6dcn, dc6dq)

    !$OMP PARALLEL DEFAULT(NONE) REDUCTION(+:energies, gradients, stress, dEdq, dEdcn) &
    !$OMP SHARED(nAtom, species, nNeighbour, iNeighbour, coords, img2CentCell, neighDist2) &
    !$OMP SHARED(c6, dc6dcn, dc6dq, calc) &
    !$OMP PRIVATE(iAt1, iSp1, iNeigh2, iNeigh3, iAt2, iAt2f, iAt3, iAt3f, iSp2, iSp3) &
    !$OMP PRIVATE(vec12, vec13, vec23, dist12, dist13, dist23, c6_12, c6_13, c6_23, c9, rr) &
    !$OMP PRIVATE(rc12, rc13, rc23, rc, r2, r1, r3, r5, fdmp, ang, dfdmp, dang, dG12, dG13, dG23) &
    !$OMP PRIVATE(dEr, sigma, dc9dcn1, dc9dcn2, dc9dcn3, dc9dq1, dc9dq2, dc9dq3)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      do iNeigh2 = 1, nNeighbour(iAt1)
        iAt2 = iNeighbour(iNeigh2, iAt1)
        vec12(:) = coords(:, iAt2) - coords(:, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        dist12 = neighDist2(iNeigh2, iAt1)
        do iNeigh3 = 1, iNeigh2 - 1
          iAt3 = iNeighbour(iNeigh3, iAt1)
          vec13(:) = coords(:, iAt3) - coords(:, iAt1)
          iAt3f = img2CentCell(iAt3)
          iSp3 = species(iAt3f)
          dist13 = neighDist2(iNeigh3, iAt1)
          vec23(:) = coords(:, iAt3) - coords(:, iAt2)
          dist23 = sum(vec23**2)

          c6_12 = c6(iAt2f, iAt1)
          c6_13 = c6(iAt3f, iAt1)
          c6_23 = c6(iAt3f, iAt2f)
          c9 = -sqrt(c6_12 * c6_13 * c6_23)

          rc12 = calc%a1 * sqrt(3.0_dp * calc%sqrtZr4r2(iSp1) * calc%sqrtZr4r2(iSp2)) + calc%a2
          rc13 = calc%a1 * sqrt(3.0_dp * calc%sqrtZr4r2(iSp1) * calc%sqrtZr4r2(iSp3)) + calc%a2
          rc23 = calc%a1 * sqrt(3.0_dp * calc%sqrtZr4r2(iSp2) * calc%sqrtZr4r2(iSp3)) + calc%a2
          rc = rc12 * rc13 * rc23

          r2 = dist12 * dist13 * dist23
          r1 = sqrt(r2)
          r3 = r2 * r1
          r5 = r3 * r2

          fdmp = 1.0_dp / (1.0_dp + 6.0_dp * (rc / r1)**(calc%alpha / 3.0_dp))
          ang = 0.375_dp*(dist12 + dist23 - dist13)*(dist12 - dist23 + dist13)&
              & *(-dist12 + dist23 + dist13) / r5 + 1.0_dp / r3

          rr = ang*fdmp

          dfdmp = -2.0_dp * calc%alpha * (rc / r1)**(calc%alpha / 3.0_dp) * fdmp**2

          ! d/dr12
          dang = -0.375_dp * (dist12**3 + dist12**2 * (dist23 + dist13)&
              & + dist12 * (3.0_dp * dist23**2 + 2.0_dp * dist23*dist13 + 3.0_dp * dist13**2)&
              & - 5.0_dp * (dist23 - dist13)**2 * (dist23 + dist13)) / r5
          dG12(:) = calc%s9 * c9 * (-dang*fdmp + ang*dfdmp) / dist12 * vec12

          ! d/dr13
          dang = -0.375_dp * (dist13**3 + dist13**2 * (dist23 + dist12)&
              & + dist13 * (3.0_dp * dist23**2 + 2.0_dp * dist23 * dist12 + 3.0_dp * dist12**2)&
              & - 5.0_dp * (dist23 - dist12)**2 * (dist23 + dist12)) / r5
          dG13(:) = calc%s9 * c9 * (-dang * fdmp + ang * dfdmp) / dist13 * vec13

          ! d/dr23
          dang = -0.375_dp * (dist23**3 + dist23**2*(dist13 + dist12)&
              & + dist23 * (3.0_dp * dist13**2 + 2.0_dp * dist13 * dist12 + 3.0_dp * dist12**2)&
              & - 5.0_dp * (dist13 - dist12)**2 * (dist13 + dist12)) / r5
          dG23(:) = calc%s9 * c9 * (-dang * fdmp + ang * dfdmp) / dist23 * vec23

          dEr = calc%s9 * rr * c9 * tripleScale(iAt1, iAt2f, iAt3f)
          energies(iAt1) = energies(iAt1) - dEr / 3.0_dp
          energies(iAt2f) = energies(iAt2f) - dEr / 3.0_dp
          energies(iAt3f) = energies(iAt3f) - dEr / 3.0_dp

          gradients(:, iAt1) = gradients(:, iAt1) - dG12 - dG13
          gradients(:, iAt2f) = gradients(:, iAt2f) + dG12 - dG23
          gradients(:, iAt3f) = gradients(:, iAt3f) + dG13 + dG23

          sigma(:,:) = spread(dG12, 1, 3) * spread(vec12, 2, 3)&
              & + spread(dG13, 1, 3) * spread(vec13, 2, 3)&
              & + spread(dG23, 1, 3) * spread(vec23, 2, 3)

          stress(:,:) = stress - sigma

          dc9dcn1 = (dc6dcn(iAt1, iAt2f) / c6_12 + dc6dcn(iAt1, iAt3f) / c6_13) * 0.5_dp
          dc9dcn2 = (dc6dcn(iAt2f, iAt1) / c6_12 + dc6dcn(iAt2f, iAt3f) / c6_23) * 0.5_dp
          dc9dcn3 = (dc6dcn(iAt3f, iAt1) / c6_13 + dc6dcn(iAt3f, iAt2f) / c6_23) * 0.5_dp
          dEdcn(iAt1) = dEdcn(iAt1) - dc9dcn1 * dEr
          dEdcn(iAt2f) = dEdcn(iAt2f) - dc9dcn2 * dEr
          dEdcn(iAt3f) = dEdcn(iAt3f) - dc9dcn3 * dEr

          dc9dq1 = (dc6dq(iAt1, iAt2f) / c6_12 + dc6dq(iAt1, iAt3f) / c6_13) * 0.5_dp
          dc9dq2 = (dc6dq(iAt2f, iAt1) /c6_12 + dc6dq(iAt2f, iAt3f) / c6_23) * 0.5_dp
          dc9dq3 = (dc6dq(iAt3f, iAt1) /c6_13 + dc6dq(iAt3f, iAt2f) /c6_23) * 0.5_dp
          dEdq(iAt1) = dEdq(iAt1) - dc9dq1*dEr
          dEdq(iAt2f) = dEdq(iAt2f) - dc9dq2*dEr
          dEdq(iAt3f) = dEdq(iAt3f) - dc9dq3*dEr
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine threeBodyDispersionGradient


  !> Driver for the calculation of DFT-D4 dispersion related properties.
  subroutine dispersionEnergy(calculator, nAtom, coords, species, neigh, img2CentCell, &
      & eeqCont, cnCont, energies, gradients, stress, volume, parEwald)

    !> DFT-D dispersion model.
    type(TDftD4Calculator), intent(in) :: calculator

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Coordinates of the atoms (including images)
    real(dp), intent(in) :: coords(:, :)

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> Updated neighbour list.
    type(TNeighbourList), intent(in) :: neigh

    !> Mapping into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> EEQ model
    type(TEeqCont), intent(inout) :: eeqCont

    !> Coordination number
    type(TCNCont), intent(inout) :: cnCont

    !> Parameter for Ewald summation.
    real(dp), intent(in), optional :: parEwald

    !> Updated energy vector at return
    real(dp), intent(out) :: energies(:)

    !> Updated gradient vector at return
    real(dp), intent(out) :: gradients(:, :)

    !> Upgraded stress
    real(dp), intent(out), optional :: stress(:, :)

    !> Volume, if system is periodic
    real(dp), intent(in), optional :: volume

    real(dp) :: sigma(3, 3)
    real(dp) :: vol, parEwald0

    if (present(volume)) then
      vol = volume
    else
      vol = 0.0_dp
    end if

    if (present(parEwald)) then
       parEwald0 = parEwald
    else
       parEwald0 = 0.0_dp
    end if

    energies(:) = 0.0_dp
    gradients(:, :) = 0.0_dp
    sigma(:, :) = 0.0_dp

    call eeqCont%updateCoords(neigh, img2CentCell, coords, species)

    call cnCont%updateCoords(neigh, img2CentCell, coords, species)

    call dispersionGradient(calculator, nAtom, coords, species, neigh, img2CentCell, &
        & cnCont%cn, cnCont%dcndr, cnCont%dcndL, eeqCont%charges, eeqCont%dqdr, eeqCont%dqdL, &
        & energies, gradients, sigma)

    if (present(stress)) then
      stress(:, :) = sigma / volume
    end if

  end subroutine dispersionEnergy


  !> charge scaling function
  elemental function zetaScale(beta1, gamma1, zref, z1) result(zeta)

    !> current effective nuclear charge
    real(dp), intent(in) :: z1

    !> reference effective nuclear charge
    real(dp), intent(in) :: zref

    !> maximum scaling
    real(dp), intent(in) :: beta1

    !> steepness of the scaling function.
    real(dp), intent(in) :: gamma1

    !> Charge scaling value
    real(dp) :: zeta

    if (z1 < 0.0_dp) then
      zeta = exp(beta1)
    else
      zeta = exp(beta1 * (1.0_dp - exp(gamma1 * (1.0_dp - zref/z1))))
    end if

  end function zetaScale


  !> derivative of charge scaling function w.r.t. charge
  elemental function dzetaScale(beta1, gamma1, zref, z1) result(dzeta)

    !> current effective nuclear charge
    real(dp), intent(in) :: z1

    !> reference effective nuclear charge
    real(dp), intent(in) :: zref

    !> maximum scaling
    real(dp), intent(in) :: beta1

    !> steepness of the scaling function.
    real(dp), intent(in) :: gamma1

    !> Derivative of the charge scaling
    real(dp) :: dzeta

    if (z1 < 0.0_dp) then
      dzeta = 0.0_dp
    else
      dzeta = - beta1 * gamma1 * exp(gamma1 * (1.0_dp - zref/z1))&
          & * zetaScale(beta1, gamma1, zref, z1) * zref / (z1**2)
    end if

  end function dzetaScale


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


  !> Logic exercise to distribute a triple energy to atomwise energies.
  !>
  !> The logic problem is outsourced to this routine to determine the factor
  !> based on the three indices so code doesn't look too overloaded.
  !> There is of course the problem that E /= 3 * E/3.0_dp for the ii'i" case,
  !> but numerics won't hit that hard and there are not to many of such triples.
  elemental function tripleScale(ii, jj, kk) result(triple)

    !> Atom indices
    integer, intent(in) :: ii, jj, kk

    !> Fraction of energy
    real(dp) :: triple

    if (ii == jj) then
      if (ii == kk) then
        ! ii'i" -> 1/6
        triple = 1.0_dp/6.0_dp
      else
        ! ii'j -> 1/2
        triple = 0.5_dp
      end if
    else
      if (ii /= kk .and. jj /= kk) then
        ! ijk -> 1 (full)
        triple = 1.0_dp
      else
        ! ijj' and iji' -> 1/2
        triple = 0.5_dp
      end if
    end if

  end function tripleScale


end module dftbp_dispdftd4
