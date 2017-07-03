!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Dispersion a la Uff, similar to Thomas Heines approach in the deMon code.
!!
!! \ref L. Zheckov et al., JCTC 1, 841-847 (2005)
!!
!! \note Periodic case could be inaccurate, if two atoms are very close to
!!   each other.
!!
!! \todo Take the reciprocal lattice vectors from outside.
!!
module dispuff_module
  use assert
  use accuracy
  use simplealgebra, only : determinant33
  use lapackroutines, only : matinv
  use periodic, only: TNeighborList, getNrOfNeighborsForAll, getLatticePoints
  use constants, only: pi
  use dispiface
  use dispcommon
  implicit none
  private

  public :: DispUffInp, DispUff, DispUff_init


  !! Input structure for the van der Waals initialization.
  type :: DispUffInp
    !> potential depths (nSpecies)
    real(dp), allocatable :: energies(:)

    !> van der Waals radii (nSpecies)
    real(dp), allocatable :: distances(:)
  end type DispUffInp


  !> Internal state of the van der Waals dispersion module.
  type, extends(DispersionIface) :: DispUff
    private
    integer :: nAtom, nSpecies  ! Nr. of atoms, species
    real(dp), allocatable :: c6(:,:)  ! Prefactors for r^-6
    real(dp), allocatable :: c12(:,:)  ! Prefactors for r^-12
    real(dp), allocatable :: cPoly(:,:,:)  ! Prefactors for polynomial
    real(dp), allocatable :: r0(:,:)  ! Switching radius
    real(dp) :: rCutoff  ! Real space cutoff
    logical :: tPeriodic  ! Periodic system?
    real(dp) :: vol  ! Volume of the unit cell
    real(dp) :: eta  ! Ewald summation parameter
    real(dp) :: ewaldRCut, ewaldGCut  ! Ewald cutoff radii
    real(dp) :: c6sum  ! Sum of the c6 coeffs.
    real(dp), allocatable :: gLatPoints(:,:) ! Reciprocal lattice vectors
    real(dp), allocatable :: energies(:)
    real(dp), allocatable :: gradients(:,:)
    real(dp) :: stress(3,3) = 0.0_dp  ! stress tensor component
    logical :: coordsUpdated = .false.  ! If first coordinate update
  contains
    procedure :: updateCoords
    procedure :: updateLatVecs
    procedure :: getEnergies
    procedure :: addGradients
    procedure :: getStress
    procedure :: getRCutoff
  end type DispUff



contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Public routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Inits a DispUff instance.
  !!
  !! \param inp Specific input parameters for Slater-Kirkwood.
  !! \param nAtom Nr. of atoms in the system.
  !! \param species0 Species of every atom in the unit cell.
  !! \param latVecs Lattice vectors, if system is periodic.
  !!
  subroutine DispUff_init(this, inp, nAtom, species0, latVecs)
    type(DispUff), intent(out) :: this
    type(DispUffInp), intent(in) :: inp
    integer, intent(in) :: nAtom
    integer, intent(in), optional :: species0(:)
    real(dp), intent(in), optional :: latVecs(:,:)

    integer :: iSp1, iSp2, iAt1
    real(dp), allocatable :: dij(:,:), rij(:,:)
    real(dp) :: preU0, preU5, preU10, c6sum

    @:ASSERT(size(inp%energies) > 0)
    @:ASSERT(size(inp%distances) == size(inp%energies))
    @:ASSERT(all(inp%energies >= 0.0_dp))
    @:ASSERT(all(inp%distances >= 0.0_dp))
    @:ASSERT(present(latVecs) .eqv. present(species0))
  #:call ASSERT_CODE
    if (present(latVecs)) then
      @:ASSERT(all(shape(latVecs) == [3, 3]))
    end if
  #:endcall ASSERT_CODE

    this%nSpecies = size(inp%energies)
    this%nAtom = nAtom
    allocate(this%c6(this%nSpecies, this%nSpecies))
    allocate(this%c12(this%nSpecies, this%nSpecies))
    allocate(this%cPoly(3, this%nSpecies, this%nSpecies))
    allocate(this%r0(this%nSpecies, this%nSpecies))

    allocate(dij(this%nSpecies, this%nSpecies))
    allocate(rij(this%nSpecies, this%nSpecies))
    forall (iSp1=1:this%nSpecies, iSp2=1:this%nSpecies)
      dij(iSp1,iSp2) = sqrt(inp%energies(iSp1) * inp%energies(iSp2))
      rij(iSp1,iSp2) = sqrt(inp%distances(iSp1) * inp%distances(iSp2))
    end forall

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
      ! Cutoff for the direct summation of r^(-12) terms. To be sure, it is
      ! delivering the required accuracy, dispTol is strengthened by two orders
      ! of magnitude more.
      this%rCutoff = (maxval(this%c12) &
          & / (tolDispersion * 1.0e-2_dp))**(1.0_dp / 12.0_dp)
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
  !!
  !! \param neigh  Updated neighbor list.
  !! \param img2CentCell  Updated mapping to central cell.
  !! \param coord  Updated coordinates.
  !! \param species0  Species of the atoms in the unit cell.
  !!
  subroutine updateCoords(this, neigh, img2CentCell, coords, species0)
    class(DispUff), intent(inout) :: this
    type(TNeighborList), intent(in) :: neigh
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species0(:)

    integer, allocatable :: nNeigh(:)

    allocate(nNeigh(this%nAtom))
    call getNrOfNeighborsForAll(nNeigh, neigh, this%rCutoff)
    if (this%tPeriodic) then
      call getDispEnergyAndGrad_cluster(this%nAtom, coords, species0, nNeigh, &
          & neigh%iNeighbor, neigh%neighDist2, img2CentCell, this%c6, &
          & this%c12, this%cPoly, this%r0, this%energies, this%gradients, &
          & removeR6=.true., stress=this%stress, vol=this%vol)
      call getNrOfNeighborsForAll(nNeigh, neigh, this%ewaldRCut)
      call addDispEGr_per_species(this%nAtom, coords, species0, nNeigh, &
          &neigh%iNeighbor, neigh%neighDist2, img2CentCell, &
          &this%c6, this%eta, this%vol, this%gLatPoints, this%energies, &
          & this%gradients, this%stress)
    else
      call getDispEnergyAndGrad_cluster(this%nAtom, coords, species0, nNeigh, &
          & neigh%iNeighbor, neigh%neighDist2, img2CentCell, this%c6, &
          & this%c12, this%cPoly, this%r0, this%energies, this%gradients)
    end if

    this%coordsUpdated = .true.

  end subroutine updateCoords


  !> Notifies the object about updated lattice vectors.
  !!
  !! \param latVecs  New lattice vectors
  !!
  subroutine updateLatVecs(this, latVecs)
    class(DispUff), intent(inout) :: this
    real(dp), intent(in) :: latVecs(:,:)

    real(dp) :: recVecs(3, 3), invRecVecs(3, 3)

    @:ASSERT(this%tPeriodic)
    @:ASSERT(all(shape(latVecs) == [3, 3]))

    this%vol = abs(determinant33(latVecs))
    invRecVecs(:,:) = latVecs / (2.0_dp * pi)
    recVecs(:,:) = transpose(invRecVecs)
    call matinv(invRecVecs)
    this%eta =  getOptimalEta(latVecs, this%vol) / sqrt(2.0_dp)
    this%ewaldRCut = getMaxRDispersion(this%eta, this%c6sum, this%vol, &
        & tolDispersion)
    this%ewaldGCut = getMaxGDispersion(this%eta, this%c6sum, tolDispersion)
    call getLatticePoints(this%gLatPoints, recVecs, invRecVecs, &
        &this%ewaldGCut, onlyInside=.true., reduceByInversion=.true., &
        &withoutOrigin=.true.)
    this%gLatPoints = matmul(recVecs, this%gLatPoints)
    this%coordsUpdated = .false.

  end subroutine updateLatVecs


  !> Returns the atomic resolved energies due to the dispersion.
  !!
  !! \param energies  Contains the atomic energy contributions on exit.
  !!
  subroutine getEnergies(this, energies)
    class(DispUff), intent(inout) :: this
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
    class(DispUff), intent(inout) :: this
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(this%coordsUpdated)
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    gradients(:,:) = gradients(:,:) + this%gradients(:,:)

  end subroutine addGradients


  !> Returns the stress tensor.
  !!
  !! \param stress tensor from the dispersion
  !!
  subroutine getStress(this, stress)
    class(DispUff), intent(inout) :: this
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(this%coordsUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))

    stress = this%stress

  end subroutine getStress


  !> Estimates the real space cutoff of the dispersion interaction.
  !!
  !! \return Cutoff
  !!
  function getRCutoff(this) result(cutoff)
    class(DispUff), intent(inout) :: this
    real(dp) :: cutoff

    cutoff = max(this%rCutoff, this%ewaldRCut)

  end function getRCutoff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Returns the energy per atom and the gradients for the cluster case
  !!
  !! \param nAtom Nr. of atoms (without periodic images)
  !! \param coords Coordinates of the atoms (including images)
  !! \param species Species of every atom.
  !! \param nNeighbors Nr. of neighbors for each atom
  !! \param iNeighbor Neighborlist.
  !! \param neighDist2 Square distances of the neighbours.
  !! \param img2CentCell Mapping into the central cell.
  !! \param c6 Prefactors for the r^-6 potential
  !! \param c12 Prefactors for the r^-12 potential
  !! \param cPoly Prefactors for the polynomial part.
  !! \param r0 Distances where polynomial repr. should change to LJ.
  !! \param energies Updated energy vector at return
  !! \param gradients Updated gradient vector at return
  !! \param removeR6 If yes, the 1/r^6 term is substracted from every
  !!   interaction.
  !!
  subroutine getDispEnergyAndGrad_cluster(nAtom, coords, species, nNeighbors, &
      &iNeighbor, neighDist2, img2CentCell, c6, c12, cPoly, r0, energies, &
      &gradients, removeR6, stress, vol)
    integer, intent(in) :: nAtom
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: nNeighbors(:), iNeighbor(0:,:)
    real(dp), intent(in) :: neighDist2(0:,:)
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: c6(:,:), c12(:,:), cPoly(:,:,:), r0(:,:)
    real(dp), intent(out) :: energies(:), gradients(:,:)
    logical , intent(in), optional :: removeR6
    real(dp), intent(out), optional :: stress(:,:)
    real(dp), intent(in), optional :: vol

    integer :: iAt1, iAt2, iAt2f, iSp1, iSp2, iNeigh, ii
    real(dp) :: rr, r2, r5, r6, r10, r12, k1, k2, dE, dGr, u0, u1, u2, f6
    real(dp) :: gr(3), vec(3)

  #:call ASSERT_CODE
    if (present(stress)) then
      @:ASSERT(all(shape(stress) == [3, 3]))
    end if
  #:endcall ASSERT_CODE

    !! Cluster case => explicit sum of the contributions
    if (present(removeR6)) then
      if (removeR6) then
        f6 = 0.0_dp
      else
        f6 = 1.0_dp
      end if
    else
      f6 = 1.0_dp
    end if
    energies = 0.0_dp
    gradients = 0.0_dp
    if (present(stress)) then
      stress = 0.0_dp
    end if
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeighbors(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        vec = coords(:,iAt1)-coords(:,iAt2)
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
        elseif (rr > minNeighDist) then
          ! Two atoms close: polynomial potential
          r10 = r2**5
          r5 = sqrt(r10)
          u0 = cPoly(1, iSp2, iSp1)
          u1 = cPoly(2, iSp2, iSp1)
          u2 = cPoly(3, iSp2, iSp1)
          dE = 0.5_dp * ((u0 - u1 * r5 - u2 * r10) + (1.0_dp - f6) * k1 / r6)
          dGr = (-5.0_dp * u1 * r5 - 10.0_dp * u2 * r10 &
              &- 6.0_dp * k1 * (1.0_dp - f6) / r6) / rr
        else
          ! Two atoms at the same position -> forget it
          dE = 0.0_dp
          dGr = 0.0_dp
        end if
        energies(iAt1) = energies(iAt1) + dE
        if (iAt1 /= iAt2f) then
          energies(iAt2f) = energies(iAt2f) + dE
        end if
        gr(:) =  dGr * (coords(:,iAt1) - coords(:,iAt2)) / rr
        gradients(:,iAt1) = gradients(:,iAt1) + gr
        gradients(:,iAt2f) = gradients(:,iAt2f) - gr
        if (present(stress)) then
          if (iAt1 /= iAt2f) then
            do ii = 1, 3
              stress(:,ii) = stress(:,ii) - gr * vec(ii) / vol
            end do
          else
            do ii = 1, 3
              stress(:,ii) = stress(:,ii) - 0.5_dp * gr * vec(ii) / vol
            end do
          end if
        end if
      end do
    end do

  end subroutine getDispEnergyAndGrad_cluster

end module dispuff_module
