!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Implementation of the D4 dispersion model
module dftbp_dftb_dispdftd4
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi, symbolToNumber
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_schedule, only : distributeRangeInChunks, assembleChunks
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_charges, only : getSummedCharges
  use dftbp_dftb_coordnumber, only : TCNCont, init_ => init
  use dftbp_dftb_dftd4param, only : TDftD4Calc, TDispDftD4Inp, TDftD4Ref, &
      & TDftD4Calculator_init, TDftD4Ref_init
  use dftbp_dftb_dispiface, only : TDispersionIface
  use dftbp_dftb_encharges, only : TEeqCont, init_ => init
  use dftbp_dftb_periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use dftbp_math_blasroutines, only : gemv
  use dftbp_math_simplealgebra, only : determinant33
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private

  public :: TDispDftD4, TDispDftD4Inp, TDispDftD4_init, init, writeDftD4Info


  !> Dispersion data cache for self-consistent evaluation of DFT-D4
  type :: TScD4

    !> Number of references
    integer :: nRef

    !> Atomic charges
    real(dp), allocatable :: charges(:)

    !> Potential shift
    real(dp), allocatable :: shift(:)

    !> Self-consistent part of the dispersion energy
    real(dp), allocatable :: energies(:)

    !> Dispersion matrix
    real(dp), allocatable :: dispMat(:, :, :, :)

    !> Scratch storage
    real(dp), allocatable :: scratch(:, :)

    !> Zeta vector
    real(dp), allocatable :: zetaVec(:, :)

    !> Derivative of zeta vector w.r.t. atomic charges
    real(dp), allocatable :: zetadq(:, :)

  end type TScD4


  !> Internal state of the DFT-D4 dispersion
  type, extends(TDispersionIface) :: TDispDftD4
    private

    !> Calculator to evaluate dispersion
    type(TDftD4Calc) :: calc

    !> Reference systems for dispersion model
    type(TDftD4Ref) :: ref

    !> Self-consistent dispersion cache
    type(TScD4), allocatable :: sc

    !> Number of atoms
    integer :: nAtom

    !> Energy
    real(dp), allocatable :: energies(:)

    !> Force contributions
    real(dp), allocatable :: gradients(:,:)

    !> Lattice vectors if periodic
    real(dp) :: latVecs(3, 3)

    !> Volume of the unit cell
    real(dp) :: vol

    !> Stress tensor
    real(dp) :: stress(3, 3)

    !> Is this periodic
    logical :: tPeriodic

    !> Are the coordinates current?
    logical :: tCoordsUpdated

    !> Are the coordinates current?
    logical :: tChargesUpdated

    !> EEQ model
    type(TEeqCont), allocatable :: eeqCont

    !> Coordination number
    type(TCNCont) :: cnCont

  contains

    !> Update internal store of coordinates
    procedure :: updateCoords

    !> Update internal store of lattice vectors
    procedure :: updateLatVecs

    !> Return energy contribution
    procedure :: getEnergies

    !> Return force contribution
    procedure :: addGradients

    !> Return stress tensor contribution
    procedure :: getStress

    !> Cutoff distance in real space for dispersion
    procedure :: getRCutoff

    !> Updates with changed charges for the instance
    procedure :: updateCharges

    !> Returns the potential shift
    procedure :: addPotential

  end type TDispDftD4


  interface init
    module procedure :: TDispDftD4_init
  end interface init


contains


  !> Initialize DispDftD4 instance
  subroutine TDispDftD4_init(this, inp, nAtom, speciesNames, latVecs)

    !> Initialized instance of D4 dispersion model
    type(TDispDftD4), intent(out) :: this

    !> Specific input parameters for damping function
    type(TDispDftD4Inp), intent(in) :: inp

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Names of species
    character(*), intent(in) :: speciesNames(:)

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:, :)

    this%tPeriodic = present(latVecs)

    if (allocated(inp%eeqInput)) then
      allocate(this%eeqCont)
      call init_(this%eeqCont, inp%eeqInput, .false., .true., nAtom, latVecs)
    end if
    call init_(this%cnCont, inp%cnInput, nAtom, latVecs)

    this%tCoordsUpdated = .false.
    this%tChargesUpdated = .false.

    if (this%tPeriodic) then
      call this%updateLatVecs(latVecs)
    end if

    this%nAtom = nAtom

    allocate(this%energies(nAtom))
    allocate(this%gradients(3, nAtom))

    call TDftD4Calculator_init(this%calc, inp)
    call TDftD4Ref_init(this%ref, inp%izp, inp%chargeScale, inp%chargeSteepness)

    ! Calculation must be self-consistent if no charge model is available
    if (inp%selfConsistent) then
      allocate(this%sc)
      call TScD4_init(this%sc, this%ref, nAtom)
    end if

  end subroutine TDispDftD4_init


  !> Initialize dispersion data cache for self-consistent evaluation
  subroutine TScD4_init(this, ref, nAtom)

    !> Instance of the self-consistent dispersion data
    type(TScD4), intent(out) :: this

    !> DFT-D4 reference systems
    type(TDftD4Ref), intent(in) :: ref

    !> Number of atoms
    integer, intent(in) :: nAtom

    integer :: nRef

    nRef = maxval(ref%nRef)
    this%nRef = nRef

    allocate(this%charges(nAtom))
    allocate(this%shift(nAtom))
    allocate(this%energies(nAtom))
    allocate(this%dispmat(nRef, nAtom, nRef, nAtom))
    allocate(this%scratch(nRef, nAtom))
    allocate(this%zetaVec(nRef, nAtom))
    allocate(this%zetadq(nRef, nAtom))

  end subroutine TScD4_init


  !> Write information about DFT-D4 dispersion model
  subroutine writeDftD4Info(unit, this)
    !> Formatted unit for output
    integer, intent(in) :: unit

    !> Instance of the type to describe
    class(TDispDftD4), intent(in) :: this

    if (allocated(this%sc)) then
      write(unit, "(A)") "Using self-consistent DFT-D4 dispersion corrections"
    else
      write(unit, "(A)") "Using DFT-D4 dispersion corrections"
    end if
  end subroutine writeDftD4Info


  !> Notifies the objects about changed coordinates
  subroutine updateCoords(this, env, neigh, img2CentCell, coords, species0, stat)

    !> Instance of DFTD4 data
    class(TDispDftD4), intent(inout) :: this

    !> Updated neighbour list
    type(TNeighbourList), intent(in) :: neigh

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Updated mapping to central cell
    integer, intent(in) :: img2CentCell(:)

    !> Updated coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Species of the atoms in the unit cell
    integer, intent(in) :: species0(:)

    !> Status of operation
    type(TStatus), intent(out) :: stat

    integer, allocatable :: nNeighbour(:)

    if (this%tPeriodic) then
      call evalDispersion(this%calc, this%ref, env, this%nAtom, species0, coords, neigh,&
          & img2CentCell, this%eeqCont, this%cnCont, this%energies, this%gradients, stat,&
          & stress=this%stress, volume=this%vol)
    else
      call evalDispersion(this%calc, this%ref, env, this%nAtom, species0, coords, neigh,&
          & img2CentCell, this%eeqCont, this%cnCont, this%energies, this%gradients, stat)
    end if
    @:PROPAGATE_ERROR(stat)

    if (allocated(this%sc)) then
      allocate(nNeighbour(this%nAtom))
      call getNrOfNeighboursForAll(nNeighbour, neigh, this%calc%cutoffInter)
      call weightDispMat(this%calc, this%ref, env, this%nAtom, nNeighbour, neigh, &
          & img2CentCell, species0, this%cnCont%cn, this%sc%scratch, this%sc%dispMat)
    end if

    this%tCoordsUpdated = .true.
    this%tChargesUpdated = allocated(this%eeqCont)

  end subroutine updateCoords


  !> Notifies the object about updated lattice vectors
  subroutine updateLatVecs(this, latVecs)

    !> Instance of DFTD4 data
    class(TDispDftD4), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:, :)

    @:ASSERT(this%tPeriodic)
    @:ASSERT(all(shape(latvecs) == shape(this%latvecs)))

    this%latVecs = latVecs
    this%vol = abs(determinant33(latVecs))

    if (allocated(this%eeqCont)) then
      call this%eeqCont%updateLatVecs(latVecs)
    end if
    call this%cnCont%updateLatVecs(LatVecs)

    this%tCoordsUpdated = .false.
    this%tChargesUpdated = .false.

  end subroutine updateLatVecs


  !> Updates with changed charges for the instance.
  subroutine updateCharges(this, env, species, neigh, qq, q0, img2CentCell, orb)

    !> Data structure
    class(TDispDftD4), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Species, shape: [nAtom]
    integer, intent(in) :: species(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neigh

    !> Orbital populations
    real(dp), intent(in) :: qq(:,:,:)

    !> Reference orbital populations
    real(dp), intent(in) :: q0(:,:,:)

    !> Mapping on atoms in central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    @:ASSERT(this%tCoordsUpdated)

    if (allocated(this%sc)) then
      ! Obtain the atomic charges from the input populations
      call getSummedCharges(species, orb, qq, q0, dQAtom=this%sc%charges)
      ! We get populations, but we want partial charges
      this%sc%charges = -this%sc%charges

      ! Obtain the charge scaling and its derivative w.r.t. the atomic charges
      call scaleReferences(this%calc, this%ref, env, this%nAtom, species, &
          & this%sc%charges, this%sc%zetaVec, this%sc%zetadq)

      ! Contract the dispersion matrix with the charge scaling
      call gemv(this%sc%scratch, this%sc%dispMat, this%sc%zetaVec)

      ! Contraction with the charge derivative yields the potential
      this%sc%zetadq(:, :) = this%sc%scratch * this%sc%zetadq
      this%sc%shift(:) = sum(this%sc%zetadq, 1)

      ! Contraction with the charge scaling again yields the dispersion energy
      this%sc%zetaVec(:, :) = 0.5_dp * this%sc%scratch * this%sc%zetaVec
      this%sc%energies(:) = sum(this%sc%zetaVec, 1)

      ! Mark the charges as updated
      this%tChargesUpdated = .true.
    end if

  end subroutine updateCharges


  !> Returns the potential shift
  subroutine addPotential(this, vDisp)

    !> Instance of DFTD4 data
    class(TDispDftD4), intent(in) :: this

    !> Atomistic potential
    real(dp), intent(inout) :: vDisp(:)

    if (allocated(this%sc)) then
      ! Potential shift is w.r.t. potentials, but we calculate w.r.t. charges (switch sign)
      vDisp(:) = vDisp - this%sc%shift
    end if

  end subroutine addPotential


  !> Returns the atomic resolved energies due to the dispersion
  subroutine getEnergies(this, energies)

    !> Instance of DFTD4 data
    class(TDispDftD4), intent(inout) :: this

    !> Contains the atomic energy contributions on exit
    real(dp), intent(out) :: energies(:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(size(energies) == this%nAtom)

    if (allocated(this%sc)) then
      energies(:) = this%energies + this%sc%energies
    else
      energies(:) = this%energies
    end if


  end subroutine getEnergies


  !> Adds the atomic gradients to the provided vector
  subroutine addGradients(this, env, neigh, img2CentCell, coords, species0, &
      & gradients, stat)

    !> Instance of DFTD4 data
    class(TDispDftD4), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> List of neighbours to atoms
    type(TNeighbourList), intent(in) :: neigh

    !> Image to central cell atom index
    integer, intent(in) :: img2CentCell(:)

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> The vector to increase by the gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> Status of operation
    integer, intent(out), optional :: stat

    real(dp), allocatable :: scEnergies(:), scGradients(:, :), scStress(:, :)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    if (allocated(this%sc)) then
      allocate(scEnergies(this%nAtom), scGradients(3, this%nAtom))
      if (this%tPeriodic) allocate(scStress(3, 3))
      call addScDispGradient(this%calc, this%ref, env, this%nAtom, species0, coords, neigh, &
          & img2CentCell, this%sc%charges, this%cnCont, scEnergies, scGradients, scStress, &
          & this%vol, stat)
      this%gradients(:, :) = this%gradients + scGradients
      if (allocated(scStress)) then
        this%stress(:, :) = this%stress + scStress
      end if
    else
      if (present(stat)) stat = 0
    end if
    gradients(:,:) = gradients + this%gradients

  end subroutine addGradients


  !> Returns the stress tensor
  subroutine getStress(this, stress)

    !> Instance of DFTD4 data
    class(TDispDftD4), intent(inout) :: this

    !> Stress tensor from the dispersion
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))

    stress(:,:) = this%stress

  end subroutine getStress


  !> Estimates the real space cutoff of the dispersion interaction
  function getRCutoff(this) result(cutoff)

    !> Instance of DFTD4 data
    class(TDispDftD4), intent(inout) :: this

    !> Resulting cutoff
    real(dp) :: cutoff

    cutoff = max(this%calc%cutoffInter, this%cnCont%getRCutoff(),&
        & this%calc%cutoffThree)
    if (allocated(this%eeqCont)) then
      cutoff = max(cutoff, this%eeqCont%getRCutoff())
    end if

  end function getRCutoff


  !> Calculate the weights of the reference systems and scale them according
  !> to their partial charge difference. It also saves the derivatives w.r.t.
  !> coordination number and partial charge for later use.
  subroutine weightReferences(calc, ref, env, nAtom, species, cn, q, zetaVec, zetadq, zetadcn)

    !> DFT-D dispersion model
    type(TDftD4Calc), intent(in) :: calc

    !> DFT-D dispersion model
    type(TDftD4Ref), intent(in) :: ref

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Species of every atom
    integer, intent(in) :: species(:)

    !> Coordination number of every atom
    real(dp), intent(in) :: cn(:)

    !> Partial charge of every atom
    real(dp), intent(in) :: q(:)

    !> Weighting and scaling function for the atomic reference systems
    real(dp), intent(out) :: zetaVec(:, :)

    !> Derivative of the weight'n'scale function w.r.t. the partial charges
    real(dp), intent(out) :: zetadq(:, :)

    !> Derivative of the weight'n'scale function w.r.t. the coordination number
    real(dp), intent(out) :: zetadcn(:, :)

    integer :: iAt1, iSp1, iRef1, iCount1, iAtFirst, iAtLast
    real(dp) :: eta1, zEff1, qRef1
    real(dp) :: norm, dnorm, wf, gw, expw, expd, gwk, dgwk

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    zetaVec(:,:) = 0.0_dp
    zetadq(:,:) = 0.0_dp
    zetadcn(:,:) = 0.0_dp

    !$omp parallel do default(none) schedule(runtime) shared(iAtFirst, iAtLast) &
    !$omp shared(zetaVec, zetadq, zetadcn, calc, ref, species, q, cn) &
    !$omp private(iAt1, iSp1, eta1, zEff1, norm, dnorm, iRef1, iCount1, wf, gw) &
    !$omp private(expw, expd, gwk, dgwk, qRef1)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      eta1 = calc%eta(iSp1)
      zEff1 = calc%zEff(iSp1)
      norm = 0.0_dp
      dnorm = 0.0_dp
      do iRef1 = 1, ref%nRef(iSp1)
        do iCount1 = 1, ref%countNumber(iRef1, iSp1)
          wf = iCount1 * calc%wf
          gw = weightCN(wf, cn(iAt1), ref%cn(iRef1, iSp1))
          norm = norm + gw
          dnorm = dnorm + 2*wf*(ref%cn(iRef1, iSp1) - cn(iAt1)) * gw
        end do
      end do
      norm = 1.0_dp / norm
      do iRef1 = 1, ref%nRef(iSp1)
        expw = 0.0_dp
        expd = 0.0_dp
        do iCount1 = 1, ref%countNumber(iref1, iSp1)
          wf = iCount1 * calc%wf
          gw = weightCN(wf, cn(iAt1), ref%cn(iRef1, iSp1))
          expw = expw + gw
          expd = expd + 2.0_dp * wf * (ref%cn(iRef1, iSp1) - cn(iAt1)) * gw
        end do

        gwk = expw * norm
        if (ieee_is_nan(gwk)) then
          if (maxval(ref%cn(:ref%nRef(iSp1), iSp1))&
              & == ref%cn(iRef1, iSp1)) then
            gwk = 1.0_dp
          else
            gwk = 0.0_dp
          end if
        end if
        qRef1 = ref%charge(iRef1, iSp1) + zEff1
        zetaVec(iRef1, iAt1) = zetaScale(calc%ga, eta1, qRef1, q(iAt1) + zEff1) * gwk
        zetadq(iRef1, iAt1) = dzetaScale(calc%ga, eta1, qRef1, q(iAt1)+zEff1) * gwk

        dgwk = expd * norm - expw * dnorm * norm**2
        if (ieee_is_nan(dgwk)) then
          dgwk = 0.0_dp
        end if
        zetadcn(iRef1, iAt1) = zetaScale(calc%ga, eta1, qRef1, q(iAt1) + zEff1) * dgwk
      end do
    end do

    call assembleChunks(env, zetaVec)
    call assembleChunks(env, zetadq)
    call assembleChunks(env, zetadcn)

  end subroutine weightReferences


  !> Actual implementation of the dispersion energy, gradients and sigma tensor.
  !>
  !> This subroutine is somewhat agnostic to the dispersion model, meaning that
  !> it can be used for either D4 or D3(BJ), if the inputs like coordination
  !> number, charges and scaling parameters are chosen accordingly.
  subroutine dispersionGradient(calc, iAtFirst, iAtLast, nNeighbour, neigh, &
      & species, coords, img2CentCell, c6, dc6dq, dc6dcn, dEdq, dEdcn, &
      & energies, gradients, sigma)

    !> DFT-D dispersion model
    type(TDftD4Calc), intent(in) :: calc

    !> Number of atoms
    integer, intent(in) :: iAtFirst

    !> Number of atoms
    integer, intent(in) :: iAtLast

    !> Number of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Updated neighbour list
    type(TNeighbourList), intent(in) :: neigh

    !> Mapping into the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Species of every atom
    integer, intent(in) :: species(:)

    !> Coordinates of the atoms (including images)
    real(dp), intent(in) :: coords(:, :)

    !> C6 coefficients for all atom pairs
    real(dp), intent(in) :: c6(:, :)

    !> Derivative of the C6 w.r.t. the partial charges
    real(dp), intent(in) :: dc6dq(:, :)

    !> Derivative of the C6 w.r.t. the coordination number
    real(dp), intent(in) :: dc6dcn(:, :)

    !> Derivative of the energy w.r.t. the partial charges
    real(dp), intent(inout) :: dEdq(:)

    !> Derivative of the energy w.r.t. the coordination number
    real(dp), intent(inout) :: dEdcn(:)

    !> Updated energy vector at return
    real(dp), intent(inout) :: energies(:)

    !> Updated gradient vector at return
    real(dp), intent(inout) :: gradients(:, :)

    !> Updated sigma tensor at return
    real(dp), intent(inout) :: sigma(:, :)

    integer :: iAt1, iSp1, iNeigh, iAt2, iSp2, iAt2f
    real(dp) :: dc6, dc6dcn1, dc6dcn2, dc6dq1, dc6dq2
    real(dp) :: vec(3), grad(3), dEr, dGr, dSr(3, 3)
    real(dp) :: rc, r1, r2, r4, r5, r6, r8, r10, rc1, rc2, rc6, rc8, rc10
    real(dp) :: f6, df6, f8, df8, f10, df10

    !$omp parallel do default(none) schedule(runtime) &
    !$omp reduction(+:energies, gradients, sigma, dEdq, dEdcn) &
    !$omp shared(iAtFirst, iAtLast, species, nNeighbour, neigh, coords, img2CentCell, c6, dc6dq) &
    !$omp shared(dc6dcn, calc) &
    !$omp private(iAt1, iSp1, iNeigh, iAt2, vec, iAt2f, iSp2, r2, r1, r4, r5, r6, r8, r10) &
    !$omp private(rc, rc1, rc2, rc6, rc8, rc10, dc6, dc6dcn1, dc6dcn2, dc6dq1, dc6dq2) &
    !$omp private(dEr, dGr, grad, dSr, f6, f8, f10, df6, df8, df10)
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
        dSr(:,:) = spread(grad, 1, 3) * spread(vec, 2, 3)

        energies(iAt1) = energies(iAt1) - dEr*dc6/2
        dEdcn(iAt1) = dEdcn(iAt1) - dc6dcn1 * dEr
        dEdq(iAt1) = dEdq(iAt1) - dc6dq1 * dEr
        if (iAt1 /= iAt2f) then
          sigma(:,:) = sigma - dSr
          energies(iAt2f) = energies(iAt2f) - dEr*dc6/2
          gradients(:, iAt1) = gradients(:, iAt1) + grad
          gradients(:, iAt2f) = gradients(:, iAt2f) - grad
          dEdcn(iAt2f) = dEdcn(iAt2f) - dc6dcn2 * dEr
          dEdq(iAt2f) = dEdq(iAt2f) - dc6dq2 * dEr
        else
          sigma(:,:) = sigma - 0.5_dp * dSr
        end if

      end do
    end do

  end subroutine dispersionGradient


  !> Calculate atomic dispersion coefficients and their derivatives w.r.t.
  !> coordination number and partial charge.
  subroutine getAtomicC6(calc, ref, env, nAtom, species, zetaVec, zetadq, zetadcn, c6, dc6dcn, dc6dq)

    !> DFT-D dispersion model
    type(TDftD4Calc), intent(in) :: calc

    !> DFT-D dispersion model
    type(TDftD4Ref), intent(in) :: ref

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Species of every atom
    integer, intent(in) :: species(:)

    !> Weighting and scaling function for the atomic reference systems
    real(dp), intent(in) :: zetaVec(:, :)

    !> Derivative of the weight'n'scale function w.r.t. the partial charges
    real(dp), intent(in) :: zetadq(:, :)

    !> Derivative of the weight'n'scale function w.r.t. the coordination number
    real(dp), intent(in) :: zetadcn(:, :)

    !> C6 coefficients for all atom pairs
    real(dp), intent(out) :: c6(:, :)

    !> Derivative of the C6 w.r.t. the partial charges
    real(dp), intent(out) :: dc6dq(:, :)

    !> Derivative of the C6 w.r.t. the coordination number
    real(dp), intent(out) :: dc6dcn(:, :)

    integer :: iAtFirst, iAtLast, iAt1, iAt2, iSp1, iSp2, iRef1, iRef2
    real(dp) :: refc6, dc6, dc6dcn1, dc6dcn2, dc6dq1, dc6dq2

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    c6(:,:) = 0.0_dp
    dc6dcn(:,:) = 0.0_dp
    dc6dq(:,:) = 0.0_dp

    !$omp parallel do default(none) schedule(runtime) shared(c6, dc6dcn, dc6dq) &
    !$omp shared(calc, ref, zetaVec, zetadq, zetadcn, iAtFirst, iAtLast, species) &
    !$omp private(iAt1, iAt2, iSp1, iSp2, iRef1, iRef2, refc6) &
    !$omp private(dc6, dc6dcn1, dc6dcn2, dc6dq1, dc6dq2)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      do iAt2 = 1, iAt1
        iSp2 = species(iAt2)
        dc6 = 0.0_dp
        dc6dcn1 = 0.0_dp
        dc6dcn2 = 0.0_dp
        dc6dq1 = 0.0_dp
        dc6dq2 = 0.0_dp
        do iRef1 = 1, ref%nRef(iSp1)
          do iRef2 = 1, ref%nRef(iSp2)
            refc6 = ref%c6(iRef1, iRef2, iSp1, iSp2)
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

    call assembleChunks(env, c6)
    call assembleChunks(env, dc6dcn)
    call assembleChunks(env, dc6dq)

  end subroutine getAtomicC6


  !> Energies, derivatives and strain derivatives of the Axilrod-Teller-Muto
  !> non-additive triple-dipole contribution.
  subroutine threeBodyDispersionGradient(calc, iAtFirst, iAtLast, nNeighbour, &
      & neigh, img2CentCell, species, coords, c6, dc6dq, dc6dcn, dEdq, dEdcn, &
      & energies, gradients, sigma)

    !> DFT-D dispersion model
    type(TDftD4Calc), intent(in) :: calc

    !> Number of atoms
    integer, intent(in) :: iAtFirst

    !> Number of atoms
    integer, intent(in) :: iAtLast

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Updated neighbour list
    type(TNeighbourList), intent(in) :: neigh

    !> Mapping into the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Coordinates of the atoms (including images)
    real(dp), intent(in) :: coords(:, :)

    !> Species of every atom
    integer, intent(in) :: species(:)

    !> C6 coefficients for all atom pairs
    real(dp), intent(in) :: c6(:, :)

    !> Derivative of the C6 w.r.t. the partial charges
    real(dp), intent(in) :: dc6dq(:, :)

    !> Derivative of the C6 w.r.t. the coordination number
    real(dp), intent(in) :: dc6dcn(:, :)

    !> Derivative of the energy w.r.t. the partial charges
    real(dp), intent(inout) :: dEdq(:)

    !> Derivative of the energy w.r.t. the coordination number
    real(dp), intent(inout) :: dEdcn(:)

    !> Updated energy vector at return (misses CN and q contributions)
    real(dp), intent(inout) :: energies(:)

    !> Updated gradient vector at return (misses CN and q contributions)
    real(dp), intent(inout) :: gradients(:, :)

    !> Updated sigma tensor at return (misses CN and q contributions)
    real(dp), intent(inout) :: sigma(:, :)

    integer :: iAt1, iAt2, iAt3, iAt2f, iAt3f, iSp1, iSp2, iSp3
    integer :: iNeigh2, iNeigh3
    real(dp) :: vec12(3), vec13(3), vec23(3), dist12, dist13, dist23
    real(dp) :: c9, c6_12, c6_13, c6_23, rc12, rc13, rc23, rc
    real(dp) :: r1, r2, r3, r5, rr, fdmp, dfdmp, ang, dang
    real(dp) :: dEr, dG12(3), dG13(3), dG23(3), dSr(3, 3)
    real(dp) :: dc9dcn1, dc9dcn2, dc9dcn3, dc9dq1, dc9dq2, dc9dq3

    !$omp parallel do default(none) schedule(runtime) &
    !$omp reduction(+:energies, gradients, sigma, dEdq, dEdcn) &
    !$omp shared(iAtFirst, iAtLast, species, nNeighbour, coords) &
    !$omp shared(img2CentCell, neigh, c6, dc6dcn, dc6dq, calc) &
    !$omp private(iAt1, iSp1, iNeigh2, iNeigh3, iAt2, iAt2f, iAt3, iAt3f, iSp2, iSp3) &
    !$omp private(vec12, vec13, vec23, dist12, dist13, dist23, c6_12, c6_13, c6_23, c9, rr) &
    !$omp private(rc12, rc13, rc23, rc, r2, r1, r3, r5, fdmp, ang, dfdmp, dang, dG12, dG13, dG23) &
    !$omp private(dEr, dSr, dc9dcn1, dc9dcn2, dc9dcn3, dc9dq1, dc9dq2, dc9dq3)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      do iNeigh2 = 1, nNeighbour(iAt1)
        iAt2 = neigh%iNeighbour(iNeigh2, iAt1)
        vec12(:) = coords(:, iAt2) - coords(:, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        dist12 = neigh%neighDist2(iNeigh2, iAt1)
        do iNeigh3 = 1, iNeigh2 - 1
          iAt3 = neigh%iNeighbour(iNeigh3, iAt1)
          vec13(:) = coords(:, iAt3) - coords(:, iAt1)
          iAt3f = img2CentCell(iAt3)
          iSp3 = species(iAt3f)
          dist13 = neigh%neighDist2(iNeigh3, iAt1)
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

          dSr(:,:) = spread(dG12, 1, 3) * spread(vec12, 2, 3)&
              & + spread(dG13, 1, 3) * spread(vec13, 2, 3)&
              & + spread(dG23, 1, 3) * spread(vec23, 2, 3)

          sigma(:,:) = sigma - dSr

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

  end subroutine threeBodyDispersionGradient


  !> Driver for the calculation of DFT-D4 dispersion related properties.
  subroutine evalDispersion(calc, ref, env, nAtom, species, coords, neigh, img2CentCell, &
      & eeqCont, cnCont, energies, gradients, stat, stress, volume)

    !> DFT-D dispersion model
    type(TDftD4Calc), intent(in) :: calc

    !> DFT-D dispersion model
    type(TDftD4Ref), intent(in) :: ref

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Species of every atom
    integer, intent(in) :: species(:)

    !> Coordinates of the atoms (including images)
    real(dp), intent(in) :: coords(:, :)

    !> Updated neighbour list
    type(TNeighbourList), intent(in) :: neigh

    !> Mapping into the central cell
    integer, intent(in) :: img2CentCell(:)

    !> EEQ model
    type(TEeqCont), intent(inout), optional :: eeqCont

    !> Coordination number
    type(TCNCont), intent(inout) :: cnCont

    !> Updated energy vector at return
    real(dp), intent(out) :: energies(:)

    !> Updated gradient vector at return
    real(dp), intent(out) :: gradients(:, :)

    !> Status of operation
    type(TStatus), intent(out) :: stat

    !> Upgraded stress
    real(dp), intent(out), optional :: stress(:, :)

    !> Volume, if system is periodic
    real(dp), intent(in), optional :: volume

    integer :: iAtFirst, iAtLast, nRef
    real(dp) :: sigma(3, 3)
    real(dp) :: vol

    !> Nr. of neighbours for each atom
    integer, allocatable :: nNeighbour(:)

    real(dp), allocatable :: zetaVec(:, :), zero(:)
    real(dp), allocatable :: zetadq(:, :)
    real(dp), allocatable :: zetadcn(:, :)
    real(dp), allocatable :: dEdq(:), dEdcn(:)
    real(dp), allocatable :: c6(:, :), dc6dq(:, :), dc6dcn(:, :)
    real(dp), allocatable :: localEnergies(:), localDeriv(:, :), localSigma(:, :)

    if (present(volume)) then
      vol = volume
    else
      vol = 0.0_dp
    end if

    energies(:) = 0.0_dp
    gradients(:, :) = 0.0_dp
    sigma(:, :) = 0.0_dp

    if (present(eeqCont)) then
      call eeqCont%updateCoords(neigh, img2CentCell, coords, species, stat)
      @:PROPAGATE_ERROR(stat)
    end if

    call cnCont%updateCoords(neigh, img2CentCell, coords, species)

    nRef = maxval(ref%nRef)
    allocate(nNeighbour(nAtom))
    allocate(zetaVec(nRef, nAtom), zetadq(nRef, nAtom), zetadcn(nRef, nAtom),&
        & c6(nAtom, nAtom), dc6dq(nAtom, nAtom), dc6dcn(nAtom, nAtom),&
        & dEdq(nAtom), dEdcn(nAtom), zero(nAtom))

    dEdq(:) = 0.0_dp
    dEdcn(:) = 0.0_dp
    zero(:) = 0.0_dp

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)
    allocate(localEnergies(nAtom), localDeriv(3, nAtom), localSigma(3, 3))
    localEnergies(:) = 0.0_dp
    localDeriv(:,:) = 0.0_dp
    localSigma(:,:) = 0.0_dp

    if (present(eeqCont)) then
      call getNrOfNeighboursForAll(nNeighbour, neigh, calc%cutoffInter)

      call weightReferences(calc, ref, env, nAtom, species, cnCont%cn, &
          & eeqCont%charges, zetaVec, zetadq, zetadcn)

      call getAtomicC6(calc, ref, env, nAtom, species, zetaVec, zetadq, zetadcn,&
          & c6, dc6dcn, dc6dq)

      call dispersionGradient(calc, iAtFirst, iAtLast, nNeighbour, neigh, &
          & species, coords, img2CentCell, c6, dc6dq, dc6dcn, dEdq, dEdcn, &
          & localEnergies, localDeriv, localSigma)
    end if

    if (calc%s9 > 0.0_dp) then
      call getNrOfNeighboursForAll(nNeighbour, neigh, calc%cutoffThree)

      call weightReferences(calc, ref, env, nAtom, species, cnCont%cn, zero, &
          & zetaVec, zetadq, zetadcn)
      zetadq(:, :) = 0.0_dp  ! really make sure there is no q `dependency'

      call getAtomicC6(calc, ref, env, nAtom, species, zetaVec, zetadq, zetadcn, &
          & c6, dc6dcn, dc6dq)

      call threeBodyDispersionGradient(calc, iAtFirst, iAtLast, nNeighbour, &
          & neigh, img2CentCell, species, coords, c6, dc6dq, dc6dcn, &
          & dEdq, dEdcn, localEnergies, localDeriv, localSigma)
    end if

    call assembleChunks(env, localEnergies)
    call assembleChunks(env, localDeriv)
    call assembleChunks(env, localSigma)

    energies(:) = localEnergies
    gradients(:,:) = localDeriv
    sigma(:,:) = localSigma

    call assembleChunks(env, dEdcn)
    call assembleChunks(env, dEdq)

    ! handle CN and charge contributions to the gradient by matrix-vector operation
    call cnCont%addGradients(gradients, dEdcn)
    if (present(eeqCont)) then
      call eeqCont%addGradients(gradients, dEdq)
    end if

    ! handle CN and charge contributions to the sigma tensor
    call cnCont%addStress(sigma, dEdcn)
    if (present(eeqCont)) then
      call eeqCont%addStress(sigma, dEdq)
    end if

    if (present(stress)) then
      stress(:, :) = sigma / volume
    end if

  end subroutine evalDispersion


  !> Charge scaling function
  elemental function zetaScale(beta1, gamma1, zref, z1) result(zeta)

    !> Current effective nuclear charge
    real(dp), intent(in) :: z1

    !> Reference effective nuclear charge
    real(dp), intent(in) :: zref

    !> Maximum scaling
    real(dp), intent(in) :: beta1

    !> Steepness of the scaling function
    real(dp), intent(in) :: gamma1

    !> Charge scaling value
    real(dp) :: zeta

    if (z1 < 0.0_dp) then
      zeta = exp(beta1)
    else
      zeta = exp(beta1 * (1.0_dp - exp(gamma1 * (1.0_dp - zref/z1))))
    end if

  end function zetaScale


  !> Derivative of charge scaling function w.r.t. charge
  elemental function dzetaScale(beta1, gamma1, zref, z1) result(dzeta)

    !> Current effective nuclear charge
    real(dp), intent(in) :: z1

    !> Reference effective nuclear charge
    real(dp), intent(in) :: zref

    !> Maximum scaling
    real(dp), intent(in) :: beta1

    !> Steepness of the scaling function
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

    !> Weighting factor / width of the gaussian function
    real(dp), intent(in) :: wf

    !> Current coordination number
    real(dp), intent(in) :: cn

    !> Reference coordination number
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


  !> Calculate the weights of the reference systems and scale them according
  !> to their partial charge difference. It also saves the derivatives w.r.t.
  !> coordination number and partial charge for later use.
  subroutine scaleReferences(calc, ref, env, nAtom, species, q, zetaVec, zetadq)

    !> DFT-D dispersion model
    type(TDftD4Calc), intent(in) :: calc

    !> DFT-D dispersion model
    type(TDftD4Ref), intent(in) :: ref

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Species of every atom
    integer, intent(in) :: species(:)

    !> Partial charge of every atom
    real(dp), intent(in) :: q(:)

    !> Scaling function for the atomic reference systems
    real(dp), intent(out) :: zetaVec(:, :)

    !> Derivative of the scaling function w.r.t. the partial charges
    real(dp), intent(out) :: zetadq(:, :)

    integer :: iAt1, iSp1, iRef1, iAtFirst, iAtLast
    real(dp) :: eta1, zEff1, qRef1

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    zetaVec(:,:) = 0.0_dp
    zetadq(:,:) = 0.0_dp

    !$omp parallel do default(none) schedule(runtime) shared(iAtFirst, iAtLast) &
    !$omp shared(zetaVec, zetadq, calc, ref, species, q) &
    !$omp private(iAt1, iSp1, eta1, zEff1, iRef1, qRef1)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      eta1 = calc%eta(iSp1)
      zEff1 = calc%zEff(iSp1)
      do iRef1 = 1, ref%nRef(iSp1)
        qRef1 = ref%charge(iRef1, iSp1) + zEff1
        zetaVec(iRef1, iAt1) = zetaScale(calc%ga, eta1, qRef1, q(iAt1) + zEff1)
        zetadq(iRef1, iAt1) = dzetaScale(calc%ga, eta1, qRef1, q(iAt1)+zEff1)
      end do
    end do

    call assembleChunks(env, zetaVec)
    call assembleChunks(env, zetadq)

  end subroutine scaleReferences


  !> Calculate the weights of the reference systems and scale them according
  !> to their partial charge difference. It also saves the derivatives w.r.t.
  !> coordination number and partial charge for later use.
  subroutine weightDispMat(calc, ref, env, nAtom, nNeighbour, neigh, img2CentCell, &
      & species, cn, gwVec, dispMat)

    !> DFT-D dispersion model
    type(TDftD4Calc), intent(in) :: calc

    !> DFT-D dispersion model
    type(TDftD4Ref), intent(in) :: ref

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Number of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Updated neighbour list
    type(TNeighbourList), intent(in) :: neigh

    !> Mapping into the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Species of every atom
    integer, intent(in) :: species(:)

    !> Coordination number of every atom
    real(dp), intent(in) :: cn(:)

    !> Scratch storage
    real(dp), intent(out) :: gwVec(:, :)

    !> Dispersion matrix
    real(dp), intent(out) :: dispMat(:, :, :, :)

    integer :: iAt1, iSp1, iRef1, nRef1, iCount1, iAtFirst, iAtLast
    integer :: iNeigh, iAt2, iSp2, iAt2f, iRef2, nRef2
    real(dp) :: eta1, zEff1, qRef1
    real(dp) :: norm, dnorm, wf, gw, expw, expd, gwk
    real(dp) :: dEr, rc, r2, r4, r6, r8, r10, rc1, rc2, rc6, rc8, rc10
    real(dp) :: f6, f8, f10, dd

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    gwVec(:,:) = 0.0_dp

    !$omp parallel do default(none) schedule(runtime) shared(iAtFirst, iAtLast) &
    !$omp shared(gwVec, calc, ref, species, cn) &
    !$omp private(iAt1, iSp1, norm, dnorm, iRef1, iCount1, wf, gw) &
    !$omp private(expw, expd, gwk)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      norm = 0.0_dp
      dnorm = 0.0_dp
      do iRef1 = 1, ref%nRef(iSp1)
        do iCount1 = 1, ref%countNumber(iRef1, iSp1)
          wf = iCount1 * calc%wf
          gw = weightCN(wf, cn(iAt1), ref%cn(iRef1, iSp1))
          norm = norm + gw
          dnorm = dnorm + 2*wf*(ref%cn(iRef1, iSp1) - cn(iAt1)) * gw
        end do
      end do
      norm = 1.0_dp / norm
      do iRef1 = 1, ref%nRef(iSp1)
        expw = 0.0_dp
        expd = 0.0_dp
        do iCount1 = 1, ref%countNumber(iref1, iSp1)
          wf = iCount1 * calc%wf
          gw = weightCN(wf, cn(iAt1), ref%cn(iRef1, iSp1))
          expw = expw + gw
          expd = expd + 2.0_dp * wf * (ref%cn(iRef1, iSp1) - cn(iAt1)) * gw
        end do

        gwk = expw * norm
        if (ieee_is_nan(gwk)) then
          if (maxval(ref%cn(:ref%nRef(iSp1), iSp1))&
              & == ref%cn(iRef1, iSp1)) then
            gwk = 1.0_dp
          else
            gwk = 0.0_dp
          end if
        end if
        gwVec(iRef1, iAt1) = gwk
      end do
    end do

    call assembleChunks(env, gwVec)

    dispMat(:, :, :, :) = 0.0_dp

    !$omp parallel do default(none) schedule(runtime) shared(dispMat) &
    !$omp shared(iAtFirst, iAtLast, calc, ref, species, nNeighbour, neigh) &
    !$omp shared(img2CentCell, gwVec) private(iAt1, iSp1, iNeigh, iAt2) &
    !$omp private(iRef1, iRef2, nRef1, nRef2, iAt2f, iSp2, r2, r4, r6, r8, r10) &
    !$omp private(rc, rc1, rc2, rc6, rc8, rc10, dEr, f6, f8, f10, dd)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeighbour(iAt1)
        iAt2 = neigh%iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        r2 = neigh%neighDist2(iNeigh, iAt1)
        r4 = r2 * r2
        r6 = r2 * r4
        r8 = r2 * r6
        r10 = r2 * r8

        rc = 3.0_dp * calc%sqrtZr4r2(iSp1) * calc%sqrtZr4r2(iSp2)
        rc1 = calc%a1 * sqrt(rc) + calc%a2
        rc2 = rc1 * rc1
        rc6 = rc2 * rc2 * rc2
        rc8 = rc2 * rc6
        rc10 = rc2 * rc8

        f6 = 1.0_dp / (r6 + rc6)
        f8 = 1.0_dp / (r8 + rc8)
        f10 = 1.0_dp / (r10 + rc10)

        dEr = calc%s6 * f6 + calc%s8 * f8 * rc + calc%s10 * rc * rc * 49.0_dp / 40.0_dp * f10
        nRef1 = ref%nRef(iSp1)
        nRef2 = ref%nRef(iSp2)
        do iRef1 = 1, nRef1
          do iRef2 = 1, nRef2
            dd = - dEr * gwVec(iRef1, iAt1) * gwVec(iRef2, iAt2f) &
                & * ref%c6(iRef2, iRef1, iSp2, iSp1)
            !$omp atomic
            dispMat(iRef2, iAt2f, iRef1, iAt1) = dispMat(iRef2, iAt2f, iRef1, iAt1) + dd
          end do
        end do

        if (iAt1 /= iAt2) then
          do iRef2 = 1, nRef2
            do iRef1 = 1, nRef1
              dd = - dEr * gwVec(iRef2, iAt2f) * gwVec(iRef1, iAt1) &
                  & * ref%c6(iRef1, iRef2, iSp1, iSp2)
              !$omp atomic
              dispMat(iRef1, iAt1, iRef2, iAt2f) = dispMat(iRef1, iAt1, iRef2, iAt2f) + dd
            end do
          end do
        end if

      end do
    end do

    call assembleChunks(env, dispMat)

  end subroutine weightDispMat


  subroutine addScDispGradient(calc, ref, env, nAtom, species, coords, neigh, img2CentCell, &
      & charges, cnCont, energies, gradients, stress, volume, stat)

    !> DFT-D dispersion model
    type(TDftD4Calc), intent(in) :: calc

    !> DFT-D dispersion model
    type(TDftD4Ref), intent(in) :: ref

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Species of every atom
    integer, intent(in) :: species(:)

    !> Coordinates of the atoms (including images)
    real(dp), intent(in) :: coords(:, :)

    !> Updated neighbour list
    type(TNeighbourList), intent(in) :: neigh

    !> Mapping into the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Atomic partial charges
    real(dp), intent(in) :: charges(:)

    !> Coordination number
    type(TCNCont), intent(inout) :: cnCont

    !> Updated energy vector at return
    real(dp), intent(out) :: energies(:)

    !> Updated gradient vector at return
    real(dp), intent(out) :: gradients(:, :)

    !> Upgraded stress
    real(dp), intent(out), optional :: stress(:, :)

    !> Volume, if system is periodic
    real(dp), intent(in), optional :: volume

    !> Status of operation
    integer, intent(out), optional :: stat

    integer :: iAtFirst, iAtLast, nRef
    real(dp) :: sigma(3, 3)
    real(dp) :: vol

    !> Nr. of neighbours for each atom
    integer, allocatable :: nNeighbour(:)

    real(dp), allocatable :: zetaVec(:, :), zero(:)
    real(dp), allocatable :: zetadq(:, :)
    real(dp), allocatable :: zetadcn(:, :)
    real(dp), allocatable :: dEdq(:), dEdcn(:)
    real(dp), allocatable :: c6(:, :), dc6dq(:, :), dc6dcn(:, :)
    real(dp), allocatable :: localEnergies(:), localDeriv(:, :), localSigma(:, :)

    if (present(volume)) then
      vol = volume
    else
      vol = 0.0_dp
    end if

    energies(:) = 0.0_dp
    gradients(:, :) = 0.0_dp
    sigma(:, :) = 0.0_dp

    nRef = maxval(ref%nRef)
    allocate(nNeighbour(nAtom))
    allocate(zetaVec(nRef, nAtom), zetadq(nRef, nAtom), zetadcn(nRef, nAtom),&
        & c6(nAtom, nAtom), dc6dq(nAtom, nAtom), dc6dcn(nAtom, nAtom),&
        & dEdq(nAtom), dEdcn(nAtom), zero(nAtom))

    dEdq(:) = 0.0_dp
    dEdcn(:) = 0.0_dp
    zero(:) = 0.0_dp

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)
    allocate(localEnergies(nAtom), localDeriv(3, nAtom), localSigma(3, 3))
    localEnergies(:) = 0.0_dp
    localDeriv(:,:) = 0.0_dp
    localSigma(:,:) = 0.0_dp

    call getNrOfNeighboursForAll(nNeighbour, neigh, calc%cutoffInter)

    call weightReferences(calc, ref, env, nAtom, species, cnCont%cn, &
        & charges, zetaVec, zetadq, zetadcn)

    call getAtomicC6(calc, ref, env, nAtom, species, zetaVec, zetadq, zetadcn,&
        & c6, dc6dcn, dc6dq)

    call dispersionGradient(calc, iAtFirst, iAtLast, nNeighbour, neigh, &
        & species, coords, img2CentCell, c6, dc6dq, dc6dcn, dEdq, dEdcn, &
        & localEnergies, localDeriv, localSigma)

    call assembleChunks(env, localEnergies)
    call assembleChunks(env, localDeriv)
    call assembleChunks(env, localSigma)

    energies(:) = localEnergies
    gradients(:,:) = localDeriv
    sigma(:,:) = localSigma

    call assembleChunks(env, dEdcn)

    ! handle CN and charge contributions to the gradient by matrix-vector operation
    call cnCont%addGradients(gradients, dEdcn)

    ! handle CN and charge contributions to the sigma tensor
    call cnCont%addStress(sigma, dEdcn)

    if (present(stress)) then
      stress(:, :) = sigma / volume
    end if

    if (present(stat)) stat = 0

  end subroutine addScDispGradient


end module dftbp_dftb_dispdftd4
