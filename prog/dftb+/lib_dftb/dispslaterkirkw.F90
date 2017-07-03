!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Dispersion a la Slater-Kirkwood as implemented by M. Elstner in old DFTB.
!!
!! The expression as found in the old DFTB had been reimplemented. The periodic
!! case had been completely rewritten using a correct Ewald summation (instead
!! of the pure real space summation, which converges very poorly).
!!
!! See also: Elstner et al., J. Chem. Phys., 114, 5149 (2001)
!!
!! \note The expression for C6(iAt1,iAt2) is not the same as in the reference
!!     paper by M. Elstner, but the one found in the old code (implemented by
!!     him as well). Furthermore, the expression for the dispersion energy in
!!     the paper (eq. 9) misses a factor of 1/2.
!!
!! \todo The generation of the reciprocal lattice vectors should not be done
!!     localy, but somewhere outside, since the Coulomb module does the same.
!!
module dispslaterkirkw
  use assert
  use accuracy
  use simplealgebra, only : determinant33
  use lapackroutines, only : matinv
  use periodic, only: TNeighborList, getNrOfNeighborsForAll, getLatticePoints
  use constants, only : pi
  use dispiface
  use dispcommon
  use message
  implicit none
  private

  public :: DispSlaKirkInp, DispSlaKirk, DispSlaKirk_init


  !> Contains the initialisation data for the Slater-Kirkwood module.
  !!
  type :: DispSlaKirkInp
    !> Atomic polarisabilities (nAtom)
    real(dp), allocatable :: polar(:)

    !> Van der Waals radii (nAtom)
    real(dp), allocatable :: rWaals(:)

    !> Effective charges (nAtom)
    real(dp), allocatable :: charges(:)
  end type DispSlaKirkInp


  !> Data for the Slater-Kirkwood type dispersion.
  !!
  type, extends(DispersionIface) :: DispSlaKirk
    private
    real(dp), allocatable :: c6(:,:)  ! Atomic polarisabilities (nAtom)
    real(dp), allocatable :: rVdW2(:,:)  ! Van der Waals radii (nAtom)
    integer :: nAtom  ! Nr. of atoms (without images)
    real(dp), allocatable :: energies(:)  ! Energies
    real(dp), allocatable :: gradients(:,:) ! Gradients (3, nAtom)
    real(dp) :: stress(3,3)  ! stress tensor components
    logical :: tPeriodic  ! If system is periodic
    real(dp) :: rCutoff  ! Real space cutoff
    real(dp) :: gCutoff  ! Reciprocal space cutoff
    real(dp) :: dampCutoff  ! Cutoff, where damping function = 1
    real(dp) :: eta  ! Periodic summation parameter
    real(dp) :: vol  ! Volume of the unit cell
    real(dp) :: maxR
    real(dp), allocatable :: gLatPoint(:,:)  ! Temporary dirty solution
    logical :: coordsUpdated = .false.  ! If first coordinate update done
  contains
    procedure :: updateCoords
    procedure :: updateLatVecs
    procedure :: getEnergies
    procedure :: addGradients
    procedure :: getStress
    procedure :: getRCutoff
  end type DispSlaKirk


  !! Some magic constants for the damping function (see paper)
  integer, parameter :: nn_ = 7         ! N
  integer, parameter :: mm_ = 4         ! M
  real(dp), parameter :: dd_ = 3.0_dp   ! d


contains

  !> Initializes a SlaterKirkwood instance.
  !!
  !! \param this  Initialized instance on exit.
  !! \param inp Input parameters for Slater-Kirkwood.
  !! \param latVecs Lattice vectors, if the system is periodic.
  !!
  subroutine DispSlaKirk_init(this, inp, latVecs)
    type(DispSlaKirk), intent(out) :: this
    type(DispSlaKirkInp), intent(in) :: inp
    real(dp), intent(in), optional :: latVecs(:,:)

    real(dp) :: recVecs(3, 3), invRecVecs(3, 3)
    real(dp) :: tol, rTmp, c6sum
    integer :: iAt1, iAt2

    @:ASSERT(size(inp%polar) > 0)
    @:ASSERT(size(inp%polar) == size(inp%rWaals))
    @:ASSERT(size(inp%polar) == size(inp%charges))
    @:ASSERT(all(inp%polar >= 0.0_dp))
    @:ASSERT(all(inp%rWaals >= 0.0_dp))
  #:call ASSERT_CODE
    if (present(latVecs)) then
      @:ASSERT(all(shape(latVecs) == [3, 3]))
    end if
  #:endcall ASSERT_CODE

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
        this%c6(iAt2, iAt1) = 1.5_dp * inp%polar(iAt1) * inp%polar(iAt2)&
            & / (sqrt(inp%polar(iAt1)/ inp%charges(iAt1)) &
            & + sqrt(inp%polar(iAt2)/ inp%charges(iAt2)))
        rTmp = (inp%rWaals(iAt1)**3 + inp%rWaals(iAt2)**3) &
            &/ (inp%rWaals(iAt1)**2 + inp%rWaals(iAt2)**2)
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

      !! Scaling down optimal eta (as suggested in the literature) is purely
      !! empirical, it reduces the real space summation, and seems to yield
      !! shorter execution times. (It does not influence the result.)
      this%eta =  getOptimalEta(latVecs, this%vol) / sqrt(2.0_dp)
      c6sum = sum(abs(this%c6))
      this%rCutoff = getMaxRDispersion(this%eta, c6sum, this%vol, &
          & tolDispersion)
      !! Cutoff, beyond which dispersion is purely 1/r^6 without damping
      this%dampCutoff = getDampCutoff_(this%maxR, tolDispDamp)
      this%rCutoff = max(this%rCutoff, this%dampCutoff)
      this%gCutoff = getMaxGDispersion(this%eta, c6sum, tolDispersion)
      call getLatticePoints(this%gLatPoint, recVecs, invRecVecs, &
          & this%gCutoff, onlyInside=.true., reduceByInversion=.true., &
          & withoutOrigin=.true.)
      this%gLatPoint(:,:) = matmul(recVecs, this%gLatPoint)
    end if

    allocate(this%energies(this%nAtom))
    allocate(this%gradients(3, this%nAtom))
    this%coordsUpdated = .false.

  end subroutine DispSlaKirk_init


  !> Notifies the objects about changed coordinates.
  !!
  !! \param neigh  Updated neighbor list.
  !! \param img2CentCell  Updated mapping to central cell.
  !! \param coord  Updated coordinates.
  !! \param species0  Species of the atoms in the unit cell.
  !!
  subroutine updateCoords(this, neigh, img2CentCell, coords, species0)
    class(DispSlaKirk), intent(inout) :: this
    type(TNeighborList), intent(in) :: neigh
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species0(:)

    integer, allocatable :: nNeighReal(:) ! Neighbors for real space summation
    integer, allocatable :: nNeighDamp(:) ! Nr. of neighbors with damping

    allocate(nNeighReal(this%nAtom))
    call getNrOfNeighborsForAll(nNeighReal, neigh, this%rCutoff)
    this%energies(:) = 0.0_dp
    this%gradients(:,:) = 0.0_dp
    this%stress(:,:) = 0.0_dp
    if (this%tPeriodic) then
      !! Make Ewald summation for a pure 1/r^6 interaction
      call addDispEGr_per_atom(this%nAtom, coords, nNeighReal, &
          & neigh%iNeighbor, neigh%neighDist2, img2CentCell, this%c6, &
          & this%eta, this%vol, this%gLatPoint, this%energies, this%gradients, &
          & this%stress)
      !! Correct those terms, where damping is important
      allocate(nNeighDamp(this%nAtom))
      call getNrOfNeighborsForAll(nNeighDamp, neigh, this%dampCutoff)
      call addDispEnergyAndGrad_cluster(this%nAtom, coords, nNeighDamp, &
          & neigh%iNeighbor, neigh%neighDist2, img2CentCell, this%c6, &
          & this%rVdW2, this%energies, this%gradients, dampCorrection=-1.0_dp)
    else
      call addDispEnergyAndGrad_cluster(this%nAtom, coords, nNeighReal, &
          & neigh%iNeighbor, neigh%neighDist2, img2CentCell, this%c6, &
          & this%rVdW2, this%energies, this%gradients)
    end if
    this%coordsUpdated = .true.

  end subroutine updateCoords


  !> Notifies the object about updated lattice vectors.
  !!
  !! \param latVecs  New lattice vectors
  !!
  subroutine updateLatVecs(this, latVecs)
    class(DispSlaKirk), intent(inout) :: this
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
    !! Cutoff, beyond which dispersion is purely 1/r^6 without damping
    this%dampCutoff = getDampCutoff_(this%maxR, tolDispDamp)
    this%rCutoff = max(this%rCutoff, this%dampCutoff)
    this%gCutoff = getMaxGDispersion(this%eta, c6sum, tolDispersion)
    call getLatticePoints(this%gLatPoint, recVecs, invRecVecs, &
        & this%gCutoff, onlyInside=.true., reduceByInversion=.true., &
        & withoutOrigin=.true.)
    this%gLatPoint(:,:) = matmul(recVecs, this%gLatPoint)
    this%coordsUpdated = .false.

  end subroutine updateLatVecs


  !> Returns the atomic resolved energies due to the dispersion.
  !!
  !! \param energies  Contains the atomic energy contributions on exit.
  !!
  subroutine getEnergies(this, energies)
    class(DispSlaKirk), intent(inout) :: this
    real(dp), intent(out) :: energies(:)

    @:ASSERT(this%coordsUpdated)
    @:ASSERT(size(energies) == this%nAtom)

    energies(:) = this%energies(:)

  end subroutine getEnergies


  !> Adds the atomic gradients to the provided vector.
  !!
  !! \param gradients  The vector to increase by the gradients.
  !!
  subroutine addGradients(this, gradients)
    class(DispSlaKirk), intent(inout) :: this
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(this%coordsUpdated)
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    gradients(:,:) = gradients(:,:) + this%gradients(:,:)

  end subroutine addGradients


  !> Returns the stress tensor.
  !!
  !! \note The stress tensor is not calculated for this dispersion model
  !!     so the program is stopped, if this method is called.
  !!
  !! \param stress tensor from the dispersion
  !!
  subroutine getStress(this, stress)
    class(DispSlaKirk), intent(inout) :: this
    real(dp), intent(out) :: stress(:,:)

    call error("Internal error: DispSlaKirk%getStress() was called")

  end subroutine getStress


  !> Estimates the real space cutoff of the dispersion interaction.
  !!
  !! \return Cutoff
  !!
  function getRCutoff(this) result(cutoff)
    class(DispSlaKirk), intent(inout) :: this
    real(dp) :: cutoff

    cutoff = this%rCutoff

  end function getRCutoff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Adds the energy per atom and the gradients for the cluster case
  !!
  !! \param nAtom Nr. of atoms (without periodic images)
  !! \param coords Coordinates of the atoms (including images)
  !! \param nNeighbors Nr. of neighbors for each atom
  !! \param iNeighbor Neighborlist.
  !! \param neighDist2 Square distances of the neighbours.
  !! \param img2CentCell Mapping into the central cell.
  !! \param c6 Van der Waals coefficients (nAtom, nAtom)
  !! \param rVdW2 Scaled inverse van der Waals radii (nAtom, nAtom)
  !! \param energies Updated energy vector at return
  !! \param gradients Updated gradient vector at return
  !! \param dampCorrection Adds the provided value to the damping function
  !!   (use -1.0 to sum up damped 1/r^6 terms and subtract pure 1/r^6 ones, in
  !!   order to correct periodic Ewald sum for the short range damped terms.)
  !!
  subroutine addDispEnergyAndGrad_cluster(nAtom, coords, nNeighbors, &
      & iNeighbor, neighDist2, img2CentCell, c6, rVdW2, energies, gradients, &
      & dampCorrection)
    integer, intent(in) :: nAtom
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nNeighbors(:)
    integer, intent(in) :: iNeighbor(0:,:)
    real(dp), intent(in) :: neighDist2(0:,:)
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: c6(:,:)
    real(dp), intent(in) :: rVdW2(:,:)
    real(dp), intent(inout) :: energies(:)
    real(dp), intent(inout) :: gradients(:,:)
    real(dp), intent(in), optional :: dampCorrection

    integer :: iAt1, iNeigh, iAt2, iAt2f
    real(dp) :: dist2, dist, h0, h1, h2, rTmp
    real(dp) :: diff(3), gr(3)
    real(dp) :: corr

    if (present(dampCorrection)) then
      corr = dampCorrection
    else
      corr = 0.0_dp
    end if

    !! Cluster case => explicit sum of the contributions
    !! NOTE: the cluster summation also (ab)used in the periodic case, neighbors
    !! may go over the cell boundary -> img2CentCell needed for folding back.
    do iAt1 = 1, nAtom
      do iNeigh = 1, nNeighbors(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        if (c6(iAt2f, iAt1) == 0.0_dp) then
          cycle
        end if
        dist2 = neighDist2(iNeigh, iAt1)
        if (dist2 > minNeighDist2) then
          dist = sqrt(dist2)
          h0 = rVdW2(iAt2f, iAt1)
          h1 = exp(-1.0_dp * h0 * dist**nn_)
          h2 = 1.0_dp - h1
          !! Energy
          rTmp = -0.5_dp * c6(iAt2f, iAt1) * (h2**mm_ + corr) / dist**6
          energies(iAt1) = energies(iAt1) + rTmp
          if (iAt1 /= iAt2f) then
            energies(iAt2f) = energies(iAt2f) + rTmp
          end if
          !! Gradients
          diff(:) = (coords(:,iAt1) - coords(:,iAt2))
          gr(:) = -c6(iAt2f, iAt1) * diff(:) &
              &* (mm_*h2**(mm_-1)*h1*h0*nn_*dist**(nn_-8) &
              &- 6.0_dp * (h2**mm_ + corr) * dist**(-8))
          gradients(:,iAt1) = gradients(:,iAt1) + gr(:)
          gradients(:,iAt2f) = gradients(:,iAt2f) - gr(:)
        end if
      end do
    end do

  end subroutine addDispEnergyAndGrad_cluster


  !! Returns the distance, beyond that the damping function equals approx. 1.
  !!
  !! \param r0 Length scaling parameter
  !! \param tol  Tolerance value.
  !! \return cutoff
  !!
  function getDampCutoff_(r0, tol) result(xx)
    real(dp), intent(in) :: r0, tol
    real(dp) :: xx

    !! solve: 1 - tol < (1-exp(-d*(r/r0)^N))^M  for r and hope that the
    !! logarithm is not blowing up your computer.
    xx = r0 * (-1.0_dp/dd_ * log(1.0_dp &
        &- (1.0_dp - tol)**(1.0_dp/real(mm_,dp))))**(1.0_dp/real(nn_, dp))

  end function getDampCutoff_


end module dispslaterkirkw
