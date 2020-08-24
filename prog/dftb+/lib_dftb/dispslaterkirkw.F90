!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Dispersion a la Slater-Kirkwood as implemented by M. Elstner in old DFTB.
!>
!> The expression as found in the old DFTB had been reimplemented. The periodic
!> case had been completely rewritten using a correct Ewald summation (instead
!> of the pure real space summation, which converges very poorly).
!>
!> See also: Elstner et al., J. Chem. Phys., 114, 5149 (2001)
!>
!> Note: The expression for C6(iAt1,iAt2) is not the same as in the reference paper by M. Elstner,
!> but the one found in the old code (implemented by him as well). Furthermore, the expression for
!> the dispersion energy in the paper (eq. 9) misses a factor of 1/2.
!>
!> Todo: The generation of the reciprocal lattice vectors should not be done localy, but somewhere
!> outside, since the Coulomb module does the same.
module dftbp_dispslaterkirkw
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_constants, only : pi
  use dftbp_dispiface
  use dftbp_dispcommon
  use dftbp_environment, only : TEnvironment
  use dftbp_lapackroutines, only : matinv
  use dftbp_message
  use dftbp_periodic, only: TNeighbourList, getNrOfNeighboursForAll, getLatticePoints
  use dftbp_schedule, only : distributeRangeInChunks, assembleChunks
  use dftbp_simplealgebra, only : determinant33
  implicit none
  private

  public :: TDispSlaKirkInp, TDispSlaKirk, DispSlaKirk_init


  !> Contains the initialisation data for the Slater-Kirkwood module.
  type :: TDispSlaKirkInp

    !> Atomic polarisabilities (nAtom)
    real(dp), allocatable :: polar(:)

    !> Van der Waals radii (nAtom)
    real(dp), allocatable :: rWaals(:)

    !> Effective charges (nAtom)
    real(dp), allocatable :: charges(:)

  end type TDispSlaKirkInp


  !> Data for the Slater-Kirkwood type dispersion.
  type, extends(TDispersionIface) :: TDispSlaKirk
    private

    !> Atomic polarisabilities (nAtom)
    real(dp), allocatable :: c6(:,:)

    !> Van der Waals radii (nAtom)
    real(dp), allocatable :: rVdW2(:,:)

    !> Nr. of atoms (without images)
    integer :: nAtom

    !> Energies
    real(dp), allocatable :: energies(:)

    !> Gradients (3, nAtom)
    real(dp), allocatable :: gradients(:,:)

    !> stress tensor components
    real(dp) :: stress(3,3)

    !> If system is periodic
    logical :: tPeriodic

    !> Real space cutoff
    real(dp) :: rCutoff

    !> Reciprocal space cutoff
    real(dp) :: gCutoff

    !> Cutoff, where damping function = 1
    real(dp) :: dampCutoff

    !> Periodic summation parameter
    real(dp) :: eta

    !> Volume of the unit cell
    real(dp) :: vol

    !> Largest pair distance
    real(dp) :: maxR

    !> Temporary dirty solution
    real(dp), allocatable :: gLatPoint(:,:)

    !> If first coordinate update done
    logical :: coordsUpdated = .false.

  contains

    !> update internal copy of coordinates
    procedure :: updateCoords

    !> update internal copy of lattice vectors
    procedure :: updateLatVecs

    !> energy contribution
    procedure :: getEnergies

    !> force contributions
    procedure :: addGradients

    !> stress tensor contribution
    procedure :: getStress

    !> real space cutoff
    procedure :: getRCutoff

  end type TDispSlaKirk


  !> Some magic constants for the damping function (see paper)
  integer, parameter :: nn_ = 7         ! N
  integer, parameter :: mm_ = 4         ! M
  real(dp), parameter :: dd_ = 3.0_dp   ! d


contains


  !> Initializes a SlaterKirkwood instance.
  subroutine DispSlaKirk_init(this, inp, latVecs)

    !> Initialized instance on exit.
    type(TDispSlaKirk), intent(out) :: this

    !> Input parameters for Slater-Kirkwood.
    type(TDispSlaKirkInp), intent(in) :: inp

    !> Lattice vectors, if the system is periodic.
    real(dp), intent(in), optional :: latVecs(:,:)

    real(dp) :: recVecs(3, 3), invRecVecs(3, 3)
    real(dp) :: tol, rTmp, c6sum
    integer :: iAt1, iAt2

    @:ASSERT(size(inp%polar) > 0)
    @:ASSERT(size(inp%polar) == size(inp%rWaals))
    @:ASSERT(size(inp%polar) == size(inp%charges))
    @:ASSERT(all(inp%polar >= 0.0_dp))
    @:ASSERT(all(inp%rWaals >= 0.0_dp))
    #:block DEBUG_CODE
    if (present(latVecs)) then
      @:ASSERT(all(shape(latVecs) == [3, 3]))
    end if
    #:endblock DEBUG_CODE

    this%nAtom = size(inp%polar)
    allocate(this%c6(this%nAtom, this%nAtom))
    allocate(this%rVdW2(this%nAtom, this%nAtom))
    this%rCutoff = 0.0_dp
    this%c6 = 0.0_dp
    this%rVdW2 = 0.0_dp
    this%maxR = 0.0_dp
    tol = epsilon(1.0_dp)
    do iAt1 = 1, this%nAtom
      if (inp%polar(iAt1) < tol .or. inp%rWaals(iAt1) < tol) then
        cycle
      end if
      do iAt2 = 1, iAt1
        this%c6(iAt2, iAt1) = 1.5_dp * inp%polar(iAt1) * inp%polar(iAt2) / (sqrt(inp%polar(iAt1)&
          & / inp%charges(iAt1)) + sqrt(inp%polar(iAt2)/ inp%charges(iAt2)))
        rTmp = (inp%rWaals(iAt1)**3 + inp%rWaals(iAt2)**3) / (inp%rWaals(iAt1)**2&
          & + inp%rWaals(iAt2)**2)
        this%rVdW2(iAt2, iAt1) = dd_ / rTmp**nn_
        this%maxR = max(this%maxR, rTmp)
        if (iAt1 /= iAt2) then
          this%c6(iAt1, iAt2) = this%c6(iAt2, iAt1)
          this%rVdW2(iAt1, iAt2) = this%rVdW2(iAt2, iAt1)
        end if
      end do
    end do
    this%rCutoff = (maxval(this%c6) / tolDispersion)**(1.0_dp / 6.0_dp)

    this%tPeriodic = present(latVecs)
    if (this%tPeriodic) then
      this%vol = abs(determinant33(latVecs))
      invRecVecs(:,:) = latVecs / (2.0_dp * pi)
      recVecs(:,:) = transpose(invRecVecs)
      call matinv(recVecs)

      ! Scaling down optimal eta (as suggested in the literature) is purely empirical, it reduces
      ! the real space summation, and seems to yield shorter execution times. (It does not influence
      ! the result.)
      this%eta =  getOptimalEta(latVecs, this%vol) / sqrt(2.0_dp)
      c6sum = sum(abs(this%c6))
      this%rCutoff = getMaxRDispersion(this%eta, c6sum, this%vol, tolDispersion)
      ! Cutoff, beyond which dispersion is purely 1/r^6 without damping
      this%dampCutoff = getDampCutoff_(this%maxR, tolDispDamp)
      this%rCutoff = max(this%rCutoff, this%dampCutoff)
      this%gCutoff = getMaxGDispersion(this%eta, c6sum, tolDispersion)
      call getLatticePoints(this%gLatPoint, recVecs, invRecVecs, this%gCutoff, onlyInside=.true.,&
        & reduceByInversion=.true., withoutOrigin=.true.)
      this%gLatPoint(:,:) = matmul(recVecs, this%gLatPoint)
    end if

    allocate(this%energies(this%nAtom))
    allocate(this%gradients(3, this%nAtom))
    this%coordsUpdated = .false.

  end subroutine DispSlaKirk_init


  !> Notifies the objects about changed coordinates.
  subroutine updateCoords(this, env, neigh, img2CentCell, coords, species0, stat)

    !> The data object for dispersion
    class(TDispSlaKirk), intent(inout) :: this

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

    !> Status of operation
    integer, intent(out), optional :: stat


    ! Neighbours for real space summation
    integer, allocatable :: nNeighReal(:)

    ! Nr. of neighbours with damping
    integer, allocatable :: nNeighDamp(:)

    if (present(stat)) then
      stat = 0
    end if

    allocate(nNeighReal(this%nAtom))
    call getNrOfNeighboursForAll(nNeighReal, neigh, this%rCutoff)
    this%energies(:) = 0.0_dp
    this%gradients(:,:) = 0.0_dp
    this%stress(:,:) = 0.0_dp
    if (this%tPeriodic) then
      ! Make Ewald summation for a pure 1/r^6 interaction
      call addDispEGr_per_atom(env, this%nAtom, coords, nNeighReal, neigh%iNeighbour,&
          & neigh%neighDist2, img2CentCell, this%c6, this%eta, this%vol, this%gLatPoint,&
          & this%energies, this%gradients, this%stress)
      ! Correct those terms, where damping is important
      allocate(nNeighDamp(this%nAtom))
      call getNrOfNeighboursForAll(nNeighDamp, neigh, this%dampCutoff)
      call addDispEnergyAndGrad_cluster(env, this%nAtom, coords, nNeighDamp, neigh%iNeighbour,&
        & neigh%neighDist2, img2CentCell, this%c6, this%rVdW2, this%energies, this%gradients,&
        & dampCorrection=-1.0_dp, stress=this%stress, vol=this%vol)
    else
      call addDispEnergyAndGrad_cluster(env, this%nAtom, coords, nNeighReal, neigh%iNeighbour,&
        & neigh%neighDist2, img2CentCell, this%c6, this%rVdW2, this%energies, this%gradients)
    end if
    this%coordsUpdated = .true.

  end subroutine updateCoords


  !> Notifies the object about updated lattice vectors.
  subroutine updateLatVecs(this, latVecs)

    !> The data object for dispersion
    class(TDispSlaKirk), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    real(dp) :: recVecs(3, 3), invRecVecs(3, 3)
    real(dp) :: c6sum

    this%vol = abs(determinant33(latVecs))
    invRecVecs(:,:) = latVecs / (2.0_dp * pi)
    recVecs(:,:) = transpose(invRecVecs)
    call matinv(recVecs)
    this%eta =  getOptimalEta(latVecs, this%vol) / sqrt(2.0_dp)
    c6sum = sum(abs(this%c6))
    this%rCutoff = getMaxRDispersion(this%eta, c6sum, this%vol, tolDispersion)
    ! Cutoff, beyond which dispersion is purely 1/r^6 without damping
    this%dampCutoff = getDampCutoff_(this%maxR, tolDispDamp)
    this%rCutoff = max(this%rCutoff, this%dampCutoff)
    this%gCutoff = getMaxGDispersion(this%eta, c6sum, tolDispersion)
    call getLatticePoints(this%gLatPoint, recVecs, invRecVecs, this%gCutoff, onlyInside=.true.,&
      & reduceByInversion=.true., withoutOrigin=.true.)
    this%gLatPoint(:,:) = matmul(recVecs, this%gLatPoint)
    this%coordsUpdated = .false.

  end subroutine updateLatVecs


  !> Returns the atomic resolved energies due to the dispersion.
  subroutine getEnergies(this, energies)

    !> The data object for dispersion
    class(TDispSlaKirk), intent(inout) :: this

    !> Contains the atomic energy contributions on exit.
    real(dp), intent(out) :: energies(:)

    @:ASSERT(this%coordsUpdated)
    @:ASSERT(size(energies) == this%nAtom)

    energies(:) = this%energies

  end subroutine getEnergies


  !> Adds the atomic gradients to the provided vector.
  subroutine addGradients(this, gradients)

    !> The data object for dispersion
    class(TDispSlaKirk), intent(inout) :: this

    !> The vector to increase by the gradients.
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(this%coordsUpdated)
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    gradients(:,:) = gradients + this%gradients

  end subroutine addGradients


  !> Returns the stress tensor.
  subroutine getStress(this, stress)

    !> The data object for dispersion
    class(TDispSlaKirk), intent(inout) :: this

    !> tensor from the dispersion
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(this%coordsUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))

    stress(:,:) = this%stress

  end subroutine getStress


  !> Estimates the real space cutoff of the dispersion interaction.
  function getRCutoff(this) result(cutoff)

    !> The data object for dispersion
    class(TDispSlaKirk), intent(inout) :: this

    !> Cutoff for the interaction
    real(dp) :: cutoff

    cutoff = this%rCutoff

  end function getRCutoff


  !> Adds the energy per atom and the gradients for the cluster case
  subroutine addDispEnergyAndGrad_cluster(env, nAtom, coords, nNeighbourSK, iNeighbour, neighDist2,&
      & img2CentCell, c6, rVdW2, energies, gradients, dampCorrection, stress, vol)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Coordinates of the atoms (including images)
    real(dp), intent(in) :: coords(:,:)

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> Neighbourlist.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Square distances of the neighbours.
    real(dp), intent(in) :: neighDist2(0:,:)

    !> Mapping into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Van der Waals coefficients (nAtom, nAtom)
    real(dp), intent(in) :: c6(:,:)

    !> Scaled inverse van der Waals radii (nAtom, nAtom)
    real(dp), intent(in) :: rVdW2(:,:)

    !> Updated energy vector at return
    real(dp), intent(inout) :: energies(:)

    !> Updated gradient vector at return
    real(dp), intent(inout) :: gradients(:,:)

    !> Adds the provided value to the damping function (use -1.0 to sum up damped 1/r^6 terms and
    !> subtract pure 1/r^6 ones, in order to correct periodic Ewald sum for the short range damped
    !> terms.)
    real(dp), intent(in), optional :: dampCorrection

    !> Stress tensor
    real(dp), intent(inout), optional :: stress(:,:)

    !> Volume of the unit cell
    real(dp), intent(in), optional :: vol

    integer :: iAtFirst, iAtLast, iAt1, iNeigh, iAt2, iAt2f, ii
    real(dp) :: dist2, dist, h0, h1, h2, rTmp
    real(dp) :: vec(3), gr(3)
    real(dp) :: corr
    real(dp), allocatable :: localDeriv(:,:), localSigma(:, :), localEnergies(:)

  #:block DEBUG_CODE
    if (present(stress)) then
      @:ASSERT(all(shape(stress) == [3, 3]))
    end if
  #:endblock DEBUG_CODE
    @:ASSERT(present(stress) .eqv. present(vol))

    if (present(dampCorrection)) then
      corr = dampCorrection
    else
      corr = 0.0_dp
    end if

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    allocate(localEnergies(nAtom), localDeriv(3, nAtom), localSigma(3, 3))
    localEnergies(:) = 0.0_dp
    localDeriv(:,:) = 0.0_dp
    localSigma(:,:) = 0.0_dp

    ! Cluster case => explicit sum of the contributions NOTE: the cluster summation also (ab)used in
    ! the periodic case, neighbours may go over the cell boundary -> img2CentCell needed for folding
    ! back.
    !$omp parallel do default(none) schedule(runtime) &
    !$omp reduction(+:localEnergies, localDeriv, localSigma) &
    !$omp shared(iAtFirst, iAtLast, nNeighbourSK, iNeighbour, img2CentCell, c6) &
    !$omp shared(neighDist2, rVdW2, coords, corr) &
    !$omp private(iAt1, iNeigh, iAt2, iAt2f, dist2, dist, h0, h1, h2, rTmp, vec, gr, ii)
    do iAt1 = iAtFirst, iAtLast
      do iNeigh = 1, nNeighbourSK(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        if (c6(iAt2f, iAt1) == 0.0_dp) then
          cycle
        end if
        dist2 = neighDist2(iNeigh, iAt1)
        dist = sqrt(dist2)
        h0 = rVdW2(iAt2f, iAt1)
        h1 = exp(-1.0_dp * h0 * dist**nn_)
        h2 = 1.0_dp - h1
        ! Energy
        rTmp = -0.5_dp * c6(iAt2f, iAt1) * (h2**mm_ + corr) / dist**6
        localEnergies(iAt1) = localEnergies(iAt1) + rTmp
        if (iAt1 /= iAt2f) then
          localEnergies(iAt2f) = localEnergies(iAt2f) + rTmp
        end if
        ! Gradients
        vec(:) = (coords(:,iAt1) - coords(:,iAt2))
        gr(:) = -c6(iAt2f, iAt1) * vec * (mm_*h2**(mm_-1)*h1*h0*nn_*dist**(nn_-8)&
          & - 6.0_dp * (h2**mm_ + corr) * dist**(-8))
        localDeriv(:,iAt1) = localDeriv(:,iAt1) + gr
        localDeriv(:,iAt2f) = localDeriv(:,iAt2f) - gr
        if (iAt1 /= iAt2f) then
          do ii = 1, 3
            localSigma(:,ii) = localSigma(:,ii) - gr * vec(ii)
          end do
        else
          do ii = 1, 3
            localSigma(:,ii) = localSigma(:,ii) - 0.5_dp * gr * vec(ii)
          end do
        end if
      end do
    end do
    !$omp end parallel do

    call assembleChunks(env, localEnergies)
    call assembleChunks(env, localDeriv)
    call assembleChunks(env, localSigma)

    energies(:) = energies + localEnergies
    gradients(:,:) = gradients + localDeriv
    if (present(stress) .and. present(vol)) then
      stress(:,:) = stress + localSigma / vol
    end if

  end subroutine addDispEnergyAndGrad_cluster


  !> Returns the distance, beyond that the damping function equals approx. 1.
  function getDampCutoff_(r0, tol) result(xx)

    !> Length scaling parameter
    real(dp), intent(in) :: r0

    !> Tolerance value.
    real(dp), intent(in) :: tol

    !> cutoff
    real(dp) :: xx

    ! solve: 1 - tol < (1-exp(-d*(r/r0)^N))^M for r and hope that the logarithm is not blowing up
    ! your computer.
    xx = r0 * (-1.0_dp/dd_ * log(1.0_dp&
      & - (1.0_dp - tol)**(1.0_dp/real(mm_,dp))))**(1.0_dp/real(nn_, dp))

  end function getDampCutoff_


end module dftbp_dispslaterkirkw
