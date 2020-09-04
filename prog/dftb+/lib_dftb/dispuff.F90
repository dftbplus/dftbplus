!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Dispersion a la UFF, similar to Thomas Heine's approach in the deMon code.
!>
!> See L. Zheckov et al., JCTC 1, 841-847 (2005)
!>
!> Note: Periodic case could be inaccurate, if two atoms are very close to each other.
!>
!> To Do: Take the reciprocal lattice vectors from outside.
!>
module dftbp_dispuff
  use dftbp_accuracy
  use dftbp_assert
  use dftbp_constants, only: pi
  use dftbp_dispiface
  use dftbp_dispcommon
  use dftbp_environment, only : TEnvironment
  use dftbp_lapackroutines, only : matinv
  use dftbp_periodic, only: TNeighbourList, getNrOfNeighboursForAll, getLatticePoints
  use dftbp_schedule, only : distributeRangeInChunks, assembleChunks
  use dftbp_simplealgebra, only : determinant33
  implicit none
  private

  public :: TDispUffInp, TDispUff, DispUff_init


  !> Input structure for the van der Waals initialization.
  type :: TDispUffInp

    !> potential depths (sized as nSpecies)
    real(dp), allocatable :: energies(:)

    !> van der Waals radii (sized as nSpecies)
    real(dp), allocatable :: distances(:)

  end type TDispUffInp


  !> Internal state of the van der Waals dispersion module.
  type, extends(TDispersionIface) :: TDispUff
    private

    !> Nr. of atoms, species
    integer :: nAtom, nSpecies

    !> Prefactors for r^{-6}
    real(dp), allocatable :: c6(:,:)

    !> Prefactors for r^{-12}
    real(dp), allocatable :: c12(:,:)

    !> Prefactors for polynomial
    real(dp), allocatable :: cPoly(:,:,:)

    !> Switching radius
    real(dp), allocatable :: r0(:,:)

    !> Real space cutoff
    real(dp) :: rCutoff

    !> Periodic system?
    logical :: tPeriodic

    !> Volume of the unit cell
    real(dp) :: vol

    !> Ewald summation parameter
    real(dp) :: eta

    !> Ewald cutoff radii
    real(dp) :: ewaldRCut, ewaldGCut

    !> Sum of the c6 coeffs.
    real(dp) :: c6sum

    !> Reciprocal lattice vectors
    real(dp), allocatable :: gLatPoints(:,:)

    !> energies for atoms
    real(dp), allocatable :: energies(:)

    !> gradient contribution
    real(dp), allocatable :: gradients(:,:)

    !> stress tensor component
    real(dp) :: stress(3,3) = 0.0_dp

    !> If first coordinate update
    logical :: coordsUpdated = .false.

  contains

    !> Notifies the objects about changed coordinates.
    procedure :: updateCoords

    !> Notifies the object about updated lattice vectors.
    procedure :: updateLatVecs

    !> Returns the atomic resolved energies due to the dispersion.
    procedure :: getEnergies

    !> Adds the atomic gradients to the provided vector.
    procedure :: addGradients

    !> Returns the stress tensor.
    procedure :: getStress

    !> Estimates the real space cutoff of the dispersion interaction.
    procedure :: getRCutoff

  end type TDispUff


contains


  !> Inits a DispUff instance.
  subroutine DispUff_init(this, inp, nAtom, species0, latVecs)

    !> data structure to initialise
    type(TDispUff), intent(out) :: this

    !> Specific input parameters for Slater-Kirkwood.
    type(TDispUffInp), intent(in) :: inp

    !> Nr. of atoms in the system.
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell.
    integer, intent(in), optional :: species0(:)

    !> Lattice vectors, if system is periodic.
    real(dp), intent(in), optional :: latVecs(:,:)

    integer :: iSp1, iSp2, iAt1
    real(dp), allocatable :: dij(:,:), rij(:,:)
    real(dp) :: preU0, preU5, preU10, c6sum

    @:ASSERT(size(inp%energies) > 0)
    @:ASSERT(size(inp%distances) == size(inp%energies))
    @:ASSERT(all(inp%energies >= 0.0_dp))
    @:ASSERT(all(inp%distances >= 0.0_dp))
    @:ASSERT(present(latVecs) .eqv. present(species0))
  #:block DEBUG_CODE
    if (present(latVecs)) then
      @:ASSERT(all(shape(latVecs) == [3, 3]))
    end if
  #:endblock DEBUG_CODE

    this%nSpecies = size(inp%energies)
    this%nAtom = nAtom
    allocate(this%c6(this%nSpecies, this%nSpecies))
    allocate(this%c12(this%nSpecies, this%nSpecies))
    allocate(this%cPoly(3, this%nSpecies, this%nSpecies))
    allocate(this%r0(this%nSpecies, this%nSpecies))

    allocate(dij(this%nSpecies, this%nSpecies))
    allocate(rij(this%nSpecies, this%nSpecies))
    do iSp1 = 1, this%nSpecies
      do iSp2 = 1, this%nSpecies
        dij(iSp1,iSp2) = sqrt(inp%energies(iSp1) * inp%energies(iSp2))
        rij(iSp1,iSp2) = sqrt(inp%distances(iSp1) * inp%distances(iSp2))
      end do
    end do

    this%c6 = 2.0_dp * dij * rij**6
    this%c12 = dij * rij**12

    preU0 = 396.0_dp / 25.0_dp
    preU5 = 2.0_dp**(5.0_dp/6.0_dp) * 672.0_dp / 25.0_dp
    preU10 = -(2.0_dp**(2.0_dp/3.0_dp)) * 552.0_dp / 25.0_dp
    this%cPoly(1,:,:) = preU0 * dij
    this%cPoly(2,:,:) = preU5 * dij / rij**5
    this%cPoly(3,:,:) = preU10 * dij / rij**10
    this%r0 = 2.0_dp**(-1.0_dp / 6.0_dp) * rij

    this%tPeriodic = present(latVecs)
    if (this%tPeriodic) then
      ! Cutoff for the direct summation of r^(-12) terms. To be sure, it is delivering the required
      ! accuracy, dispTol is strengthened by two orders of magnitude more.
      this%rCutoff = (maxval(this%c12) / (tolDispersion * 1.0e-2_dp))**(1.0_dp / 12.0_dp)
      ! Summing with loop to avoid creation of (nAtom, nAtom) tmp array.
      c6sum = 0.0_dp
      do iAt1 = 1, nAtom
        c6sum = c6sum + sum(abs(this%c6(species0,species0(iAt1))))
      end do
      this%c6sum = c6sum
      call this%updateLatVecs(latVecs)
    else
      ! Cutoff for the direct real space summation of r^(-6) terms.
      this%rCutoff = (maxval(this%c6) / tolDispersion)**(1.0_dp / 6.0_dp)
      this%ewaldRCut = 0.0_dp
    end if

    allocate(this%energies(this%nAtom))
    allocate(this%gradients(3, this%nAtom))

  end subroutine DispUff_init


  !> Notifies the objects about changed coordinates.
  subroutine updateCoords(this, env, neigh, img2CentCell, coords, species0, stat)

    !> Instance of dispersion to update
    class(TDispUff), intent(inout) :: this

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

    !> Status of operation
    integer, intent(out), optional :: stat

    integer, allocatable :: nNeigh(:)

    if (present(stat)) then
      stat = 0
    end if

    allocate(nNeigh(this%nAtom))
    call getNrOfNeighboursForAll(nNeigh, neigh, this%rCutoff)
    if (this%tPeriodic) then
      call getDispEnergyAndGrad_cluster(env, this%nAtom, coords, species0, nNeigh,&
          & neigh%iNeighbour, neigh%neighDist2, img2CentCell, this%c6, this%c12, this%cPoly,&
          & this%r0, this%energies, this%gradients, removeR6=.true., stress=this%stress,&
          & vol=this%vol)
      call getNrOfNeighboursForAll(nNeigh, neigh, this%ewaldRCut)
      call addDispEGr_per_species(env, this%nAtom, coords, species0, nNeigh, neigh%iNeighbour,&
          & neigh%neighDist2, img2CentCell, this%c6, this%eta, this%vol, this%gLatPoints,&
          & this%energies, this%gradients, this%stress)
    else
      call getDispEnergyAndGrad_cluster(env, this%nAtom, coords, species0, nNeigh,&
          & neigh%iNeighbour, neigh%neighDist2, img2CentCell, this%c6, this%c12, this%cPoly,&
          & this%r0, this%energies, this%gradients)
    end if

    this%coordsUpdated = .true.

  end subroutine updateCoords


  !> Notifies the object about updated lattice vectors.
  subroutine updateLatVecs(this, latVecs)

    !> Instance to update
    class(TDispUff), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    real(dp) :: recVecs(3, 3), invRecVecs(3, 3)

    @:ASSERT(this%tPeriodic)
    @:ASSERT(all(shape(latVecs) == [3, 3]))

    this%vol = abs(determinant33(latVecs))
    invRecVecs(:,:) = latVecs / (2.0_dp * pi)
    recVecs(:,:) = transpose(invRecVecs)
    call matinv(recVecs)
    this%eta =  getOptimalEta(latVecs, this%vol) / sqrt(2.0_dp)
    this%ewaldRCut = getMaxRDispersion(this%eta, this%c6sum, this%vol, tolDispersion)
    this%ewaldGCut = getMaxGDispersion(this%eta, this%c6sum, tolDispersion)
    call getLatticePoints(this%gLatPoints, recVecs, invRecVecs, this%ewaldGCut, onlyInside=.true.,&
        & reduceByInversion=.true., withoutOrigin=.true.)
    this%gLatPoints = matmul(recVecs, this%gLatPoints)
    this%coordsUpdated = .false.

  end subroutine updateLatVecs


  !> Returns the atomic resolved energies due to the dispersion.
  subroutine getEnergies(this, energies)

    !> Instance of dispersion
    class(TDispUff), intent(inout) :: this

    !> Contains the atomic energy contributions on exit.
    real(dp), intent(out) :: energies(:)

    @:ASSERT(this%coordsUpdated)
    @:ASSERT(size(energies) == this%nAtom)

    energies(:) = this%energies

  end subroutine getEnergies


  !> Adds the atomic gradients to the provided vector.
  subroutine addGradients(this, gradients)

    !> Instance of dispersion
    class(TDispUff), intent(inout) :: this

    !> The vector to increase by the gradients.
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(this%coordsUpdated)
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    gradients(:,:) = gradients + this%gradients

  end subroutine addGradients


  !> Returns the stress tensor.
  subroutine getStress(this, stress)

    !> Instance of dispersion
    class(TDispUff), intent(inout) :: this

    !> tensor from the dispersion
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(this%coordsUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))

    stress = this%stress

  end subroutine getStress


  !> Estimates the real space cutoff of the dispersion interaction.
  function getRCutoff(this) result(cutoff)

    !> Instance of dispersion
    class(TDispUff), intent(inout) :: this

    !> Cutoff distance
    real(dp) :: cutoff

    cutoff = max(this%rCutoff, this%ewaldRCut)

  end function getRCutoff


  !> Returns the energy per atom and the gradients for the cluster case
  subroutine getDispEnergyAndGrad_cluster(env, nAtom, coords, species, nNeighbourSK, iNeighbour,&
      & neighDist2, img2CentCell, c6, c12, cPoly, r0, energies, gradients, removeR6, stress, vol)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Coordinates of the atoms (including images)
    real(dp), intent(in) :: coords(:,:)

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> Neighbourlist.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Square distances of the neighbours.
    real(dp), intent(in) :: neighDist2(0:,:)

    !> Mapping into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Prefactors for the r^{-6} potential
    real(dp), intent(in) :: c6(:,:)

    !> Prefactors for the r^{-12} potential
    real(dp), intent(in) :: c12(:,:)

    !> Prefactors for the polynomial part.
    real(dp), intent(in) :: cPoly(:,:,:)

    !> Distances where polynomial repr. should change to LJ.
    real(dp), intent(in) :: r0(:,:)

    !> Updated energy vector at return
    real(dp), intent(out) :: energies(:)

    !> Updated gradient vector at return
    real(dp), intent(out) :: gradients(:,:)

    !> If yes, the 1/r^6 term is substracted from every interaction.
    logical , intent(in), optional :: removeR6

    !> Stress tensor
    real(dp), intent(out), optional :: stress(:,:)

    !> Volume of the unit cell
    real(dp), intent(in), optional :: vol

    integer :: iAtFirst, iAtLast
    integer :: iAt1, iAt2, iAt2f, iSp1, iSp2, iNeigh, ii
    real(dp) :: rr, r2, r5, r6, r10, r12, k1, k2, dE, dGr, u0, u1, u2, f6
    real(dp) :: gr(3), vec(3)
    real(dp), allocatable :: localDeriv(:,:), localSigma(:, :), localEnergies(:)

  #:block DEBUG_CODE
    if (present(stress)) then
      @:ASSERT(all(shape(stress) == [3, 3]))
    end if
  #:endblock DEBUG_CODE
    @:ASSERT(present(stress) .eqv. present(vol))

    ! Cluster case => explicit sum of the contributions
    if (present(removeR6)) then
      if (removeR6) then
        f6 = 0.0_dp
      else
        f6 = 1.0_dp
      end if
    else
      f6 = 1.0_dp
    end if

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    allocate(localEnergies(nAtom), localDeriv(3, nAtom), localSigma(3, 3))
    localEnergies(:) = 0.0_dp
    localDeriv(:,:) = 0.0_dp
    localSigma(:,:) = 0.0_dp

    !$omp parallel do default(none) schedule(runtime) &
    !$omp reduction(+:localEnergies, localDeriv, localSigma) &
    !$omp shared(iAtFirst, iAtLast, species, nNeighbourSK, iNeighbour, coords) &
    !$omp shared(img2CentCell, neighDist2, c6, r0, c12, cPoly, f6) &
    !$omp private(iAt1, iSp1, iNeigh, iAt2, vec, iAt2f, iSp2, r2, rr, k1, r6) &
    !$omp private(r12, k2, dE, dGr, r10, r5, u0, u1, u2, gr, ii)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeighbourSK(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        vec(:) = coords(:,iAt1)-coords(:,iAt2)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        r2 = neighDist2(iNeigh, iAt1)
        rr = sqrt(r2)
        k1 = c6(iSp2, iSp1)
        r6 = r2 * r2 * r2
        if (rr > r0(iSp2, iSp1)) then
          ! Two atoms far enough: Lennard-Jones potential
          r12 = r6 * r6
          k2 = c12(iSp2, iSp1)
          dE = 0.5_dp * (-(k1 * f6) / r6 + k2 / r12)
          dGr = (6.0_dp * k1 * f6 / r6 - 12.0_dp * k2 / r12) / rr
        else
          ! Two atoms close: polynomial potential
          r10 = r2**5
          r5 = sqrt(r10)
          u0 = cPoly(1, iSp2, iSp1)
          u1 = cPoly(2, iSp2, iSp1)
          u2 = cPoly(3, iSp2, iSp1)
          dE = 0.5_dp * ((u0 - u1 * r5 - u2 * r10) + (1.0_dp - f6) * k1 / r6)
          dGr = (-5.0_dp * u1 * r5 - 10.0_dp * u2 * r10 - 6.0_dp * k1 * (1.0_dp - f6) / r6) / rr
        end if
        localEnergies(iAt1) = localEnergies(iAt1) + dE
        if (iAt1 /= iAt2f) then
          localEnergies(iAt2f) = localEnergies(iAt2f) + dE
        end if
        gr(:) =  dGr * (coords(:,iAt1) - coords(:,iAt2)) / rr
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

    energies(:) = localEnergies
    gradients(:,:) = localDeriv
    if (present(stress) .and. present(vol)) then
      stress(:,:) = localSigma / vol
    end if

  end subroutine getDispEnergyAndGrad_cluster


end module dftbp_dispuff
