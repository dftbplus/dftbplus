!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains code to calculate the H0 Hamiltonian and overlap matrix and their
!> derivatives.
module dftbp_nonscc
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_sk
  use dftbp_slakocont
  use dftbp_commontypes
  use dftbp_message
  use dftbp_schedule
  use dftbp_environment
  implicit none
  private

  public :: buildH0, buildS
  public :: TNonSccDiff, NonSccDiff_init
  public :: diffTypes

  ! Workaround: pgfortran 16.5
  ! Keep default value as a constant, as using the expression directly
  ! for the default initialisation triggers ICE.
  real(dp), parameter :: DELTA_X_DIFF_DEFAULT = epsilon(1.0_dp)**0.25_dp

  !> Contains settings for the derivation of the non-scc contribution.
  type :: TNonSccDiff
    private
    ! default of a reasonable choice for round off and a second order formula
    real(dp) :: deltaXDiff = DELTA_X_DIFF_DEFAULT
    integer :: diffType
  contains

    !> evaluate first derivative
    procedure :: getFirstDeriv

    !> evaluate second derivative
    procedure :: getSecondDeriv

  end type TNonSccDiff

  !> Namespace for possible differentiation methods
  type :: TDiffTypesEnum
    integer :: finiteDiff
    integer :: richardson
  end type TDiffTypesEnum

  !> Actual values for diffTypes.
  type(TDiffTypesEnum), parameter :: diffTypes = TDiffTypesEnum(1, 2)

contains

  !> Driver for making the non-SCC Hamiltonian in the primitive sparse format.
  subroutine buildH0(env, ham, skHamCont, selfegy, coords, nNeighbourSK, iNeighbours, species,&
      & iPair, orb)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Returns the non-self-consistent Hamiltonian
    real(dp), intent(out) :: ham(:)

    !> Container for the SlaKo Hamiltonian integrals
    type(TSlakoCont), intent(in) :: skHamCont

    !> On-site energies for each species
    real(dp), intent(in) :: selfegy(:,:)

    !> List of all coordinates, including possible periodic images of atoms
    real(dp), intent(in) :: coords(:,:)

    !> Number of surrounding neighbours for each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> List of surrounding neighbours for each atom
    integer, intent(in) :: iNeighbours(0:,:)

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Shift vector, where the interaction between two atoms
    integer, intent(in) :: iPair(0:,:)

    !> Information about the orbitals in the system
    type(TOrbitals), intent(in) :: orb

    integer :: nAtom, iAt1, iSp1, ind, iOrb1, iAtFirst, iAtLast

    nAtom = size(nNeighbourSK)
    ham(:) = 0.0_dp

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    ! Put the on-site energies into the Hamiltonian,
    ! and <lm|l'm'> = delta_l,l' * delta_m',m' for the overlap  !'
    !$OMP PARALLEL DO PRIVATE(iAt1, iSp1, ind, iOrb1) DEFAULT(SHARED) SCHEDULE(RUNTIME)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      ind = iPair(0, iAt1) + 1
      do iOrb1 = 1, orb%nOrbAtom(iAt1)
        ham(ind) = selfegy(orb%iShellOrb(iOrb1, iSp1), iSp1)
        ind = ind + orb%nOrbAtom(iAt1) + 1
      end do
    end do
    !$OMP  END PARALLEL DO

    call buildDiatomicBlocks(iAtFirst, iAtLast, skHamCont, coords, nNeighbourSK, iNeighbours,&
        & species, iPair, orb, ham)

    call assembleChunks(env, ham)

  end subroutine buildH0

  !> Driver for making the overlap matrix in the primitive sparse format.
  subroutine buildS(env, over, skOverCont, coords, nNeighbourSK, iNeighbours, species, iPair, orb)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Returns the overlap
    real(dp), intent(out) :: over(:)

    !> Container for the SlaKo overlap integrals
    type(TSlakoCont), intent(in) :: skOverCont

    !> List of all coordinates, including possible periodic images of atoms
    real(dp), intent(in) :: coords(:,:)

    !> Number of surrounding neighbours for each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> List of surrounding neighbours for each atom
    integer, intent(in) :: iNeighbours(0:,:)

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Shift vector, where the interaction between two atoms starts in the sparse format.
    integer, intent(in) :: iPair(0:,:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    integer :: nAtom, iAt1, iSp1, ind, iOrb1, iAtFirst, iAtLast

    nAtom = size(nNeighbourSK)
    over(:) = 0.0_dp

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    ! Put 1.0 for the diagonal elements of the overlap.
    !$OMP PARALLEL DO PRIVATE(iAt1, iSp1, ind, iOrb1) DEFAULT(SHARED) SCHEDULE(RUNTIME)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      ind = iPair(0,iAt1) + 1
      do iOrb1 = 1, orb%nOrbAtom(iAt1)
        over(ind) = 1.0_dp
        ind = ind + orb%nOrbAtom(iAt1) + 1
      end do
    end do
    !$OMP  END PARALLEL DO

    call buildDiatomicBlocks(iAtFirst, iAtLast, skOverCont, coords, nNeighbourSK, iNeighbours,&
        & species, iPair, orb, over)

    call assembleChunks(env, over)

  end subroutine buildS

  !> Initializes a differentiator for the non-scc contributions.
  !> Note: Second derivative can not be calculated currently via Richardson
  !> interpolation, so the finite difference method is invoked instead.
  subroutine NonSccDiff_init(this, diffType, deltaXDiff)

    !> Initialised instance on exit.
    type(TNonSccDiff), intent(out) :: this

    !> Type of the differentiator: diffTypes%finiteDiff or diffTypes%richardson
    integer, intent(in) :: diffType

    !> Displacement for finite difference differentiation.
    real(dp), intent(in), optional :: deltaXDiff

    if (all([diffTypes%finiteDiff, diffTypes%richardson] /= diffType)) then
      call error("Invalid differentiator type in NonSccDiff_init")
    end if
    this%diffType = diffType
    if (present(deltaXDiff)) then
      this%deltaXDiff = deltaXDiff
    end if

  end subroutine NonSccDiff_init

  !> Calculates the first derivative of H0 or S.
  subroutine getFirstDeriv(this, deriv, skCont,coords, species, atomI, atomJ, orb)

    !> Instance
    class(TNonSccDiff), intent(in) :: this

    !> Derivative of H or S diatomic block, with respect to x,y,z (last index).
    real(dp), intent(out) :: deriv(:,:,:)

    !> Container for the SK integrals
    type(TSlakoCont), intent(in) :: skCont

    !> List of all coordinates, including possible periodic images of atoms
    real(dp), intent(in) :: coords(:,:)

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> The first atom in the diatomic block
    integer, intent(in) :: atomI

    !> The second atom in the diatomic block
    integer, intent(in) :: atomJ

    !> Orbital informations
    type(TOrbitals), intent(in) :: orb

    select case (this%diffType)
    case (diffTypes%finiteDiff)
      call getFirstDerivFiniteDiff(deriv, skCont, coords, species, atomI, atomJ, orb,&
          & this%deltaXDiff)
    case (diffTypes%richardson)
      call getFirstDerivRichardson(deriv, skCont, coords, species, atomI, atomJ, orb)
    end select

  end subroutine getFirstDeriv

  !> Calculates the numerical second derivative of a diatomic block of H0 or S.
  subroutine getSecondDeriv(this, deriv, skCont, coords, species, atomI, atomJ, orb)

    !> Instance.
    class(TNonSccDiff), intent(in) :: this

    !> Derivative of the diatomic block, with respect to x,y,z (last 2 indices)
    real(dp), intent(out) :: deriv(:,:,:,:)

    !> Container for SK Hamiltonian integrals
    type(TSlakoCont), intent(in) :: skCont

    !> List of all coordinates, including possible  periodic images of atoms.
    real(dp), intent(in) :: coords(:,:)

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> First atom in the diatomic block
    integer, intent(in) :: atomI

    !> Second atom in the diatomic block
    integer, intent(in) :: atomJ

    !> Orbital informations
    type(TOrbitals), intent(in) :: orb

    ! Note, second derivatives are not Richardson interpolated yet, so finite difference code is
    ! invoked instead.
    select case (this%diffType)
    case (diffTypes%finiteDiff)
      call getSecondDerivFiniteDiff(deriv, skCont, coords, species, atomI, atomJ, orb,&
          & this%deltaXDiff)
    case (diffTypes%richardson)
      call error("Richardson second derivative extrapolation not implemented")
    end select

  end subroutine getSecondDeriv


  !> Helper routine to calculate the diatomic blocks for the routines buildH0 and buildS.
  subroutine buildDiatomicBlocks(firstAtom, lastAtom, skCont, coords, nNeighbourSK, iNeighbours,&
      & species, iPair, orb, out)
    integer, intent(in) :: firstAtom
    integer, intent(in) :: lastAtom
    type(TSlakoCont), intent(in) :: skCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nNeighbourSK(:)
    integer, intent(in) :: iNeighbours(0:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: iPair(0:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(inout) :: out(:)

    real(dp) :: tmp(orb%mOrb,orb%mOrb)
    real(dp) :: vect(3), dist
    real(dp) :: interSK(getMIntegrals(skCont))
    integer :: nOrb1, nOrb2
    integer :: iAt1, iAt2, iSp1, iSp2, iNeigh1, ind

    ! Do the diatomic blocks for each of the atoms with its nNeighbourSK
    !$OMP PARALLEL DO PRIVATE(iAt1,iSp1,nOrb1,iNeigh1,iAt2,iSp2,nOrb2,ind,vect,dist,tmp,interSK) &
    !$OMP& DEFAULT(SHARED) SCHEDULE(RUNTIME)
    do iAt1 = firstAtom, lastAtom
      iSp1 = species(iAt1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      do iNeigh1 = 1, nNeighbourSK(iAt1)
        iAt2 = iNeighbours(iNeigh1, iAt1)
        iSp2 = species(iAt2)
        nOrb2 = orb%nOrbSpecies(iSp2)
        ind = iPair(iNeigh1,iAt1)
        vect(:) = coords(:,iAt2) - coords(:,iAt1)
        dist = sqrt(sum(vect**2))
        vect(:) = vect / dist
        call getSKIntegrals(skCont, interSK, dist, iSp1, iSp2)
        call rotateH0(tmp, interSK, vect(1), vect(2), vect(3), iSp1, iSp2, orb)
        out(ind + 1 : ind + nOrb2 * nOrb1) = reshape(tmp(1:nOrb2, 1:nOrb1), [nOrb2 * nOrb1])
      end do
    end do
    !$OMP  END PARALLEL DO

  end subroutine buildDiatomicBlocks


  !> Calculates the numerical derivative of a diatomic block H0 or S by finite differences.
  subroutine getFirstDerivFiniteDiff(deriv, skCont, coords, species, atomI, atomJ, orb, deltaXDiff)
    real(dp), intent(out) :: deriv(:,:,:)
    type(TSlakoCont), intent(in) :: skCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: atomI, atomJ
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: deltaXDiff

    ! interpolated H integs.
    real(dp) :: interSk(getMIntegrals(skCont))

    real(dp) :: vect(3), dist
    integer :: ii, jj
    integer :: sp1, sp2
    real(dp) :: tmp(size(deriv,dim=1), size(deriv,dim=2),2,3)

    @:ASSERT(size(deriv, dim=3) == 3)

    deriv(:,:,:) = 0.0_dp

    sp1 = species(atomI)
    sp2 = species(atomJ)

    do jj = 1, 3
      do ii = 1, 2
        vect(:) = coords(:,atomJ) - coords(:,atomI)
        vect(jj) = vect(jj) - real(2 * ii - 3, dp) * deltaXDiff
        dist = sqrt(sum(vect**2))
        vect(:) = vect / dist
        call getSKIntegrals(skCont, interSk, dist, sp1, sp2)
        call rotateH0(tmp(:,:,ii,jj), interSk, vect(1), vect(2), vect(3), sp1, sp2, orb)
      end do
    end do

    do ii = 1, 3
      deriv(:,:,ii) = (tmp(:,:,2,ii) - tmp(:,:,1,ii)) / (2.0_dp * deltaXDiff)
    end do

  end subroutine getFirstDerivFiniteDiff

  !> Calculates the numerical derivative of a diatomic block H0 or S by Richardsons method.
  subroutine getFirstDerivRichardson(deriv, skCont, coords, species, atomI, atomJ, orb)
    real(dp), intent(out) :: deriv(:,:,:)
    type(TSlakoCont), intent(in) :: skCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: atomI, atomJ
    type(TOrbitals), intent(in) :: orb

    integer, parameter :: maxrows = 20
    integer, parameter :: maxcolumns = 3
    real(dp), parameter :: tol = 1.0E-11_dp ! epsilon ?

    real(dp) :: interSk(getMIntegrals(skCont))  ! interpolated S integs.
    real(dp) :: vect(3), dist
    integer :: sp1, sp2
    real(dp) :: tmp(size(deriv, dim=1), size(deriv, dim=2),2)
    real(dp) :: diff
    real(dp) :: dd(size(deriv, dim=1), size(deriv, dim=2), 0:maxcolumns, 0:maxrows)
    logical :: tConverged(size(deriv, dim=1), size(deriv, dim=2))
    real(dp) :: hh
    integer :: ii, iCart, kk, ll, mm, nk

    @:ASSERT(size(deriv, dim=3) == 3)

    deriv(:,:,:) = 0.0_dp

    sp1 = species(atomI)
    sp2 = species(atomJ)
    do iCart = 1, 3

      ! initial finite difference, chosen as an exactly representable value in floating point
      hh = 0.125_dp
      dd(:,:,:,:) = 0.0_dp

      do kk = 1, 2
        vect(:) = coords(:,atomJ) - coords(:,atomI)
        vect(iCart) = vect(iCart) - real(2 * kk - 3, dp) * hh
        dist = sqrt(sum(vect(:)**2))
        vect(:) = vect / dist
        call getSKIntegrals(skCont, interSk, dist, sp1, sp2)
        call rotateH0(tmp(:,:,kk), interSk, vect(1), vect(2), vect(3), sp1, sp2, orb)
      end do

      ! row number
      ii = 0
      dd(:,:,0,ii) = (tmp(:,:,2) - tmp(:,:,1)) / (2.0_dp * hh)
      tConverged(:,:) = .false.

      ! want value of ii to be consistent however the loop exits, so using a while instead of for
      ! loop
      do while (ii < maxrows .and. .not. all(tConverged))
        ! using powers of 2, so should be exactly representable
        hh = hh / 2.0_dp
        if (hh <= epsilon(1.0_dp)**0.25_dp) then
          ! rounding error on order of epsilon/h^2 and central difference formula truncation error
          ! of O(h^2), so want to keep above epsilon/h^4. If the function gets large, then
          ! max(1,|f|) epsilon^(1/4) might be better
          exit
        end if
        ii = ii + 1
        do kk = 1, 2
          vect(:) = coords(:,atomJ) - coords(:,atomI)
          vect(iCart) = vect(iCart) - real(2*kk-3,dp) * hh
          dist = sqrt(sum(vect**2))
          vect(:) = vect / dist
          call getSKIntegrals(skCont, interSk, dist, sp1, sp2)
          call rotateH0(tmp(:,:,kk), interSk, vect(1), vect(2), vect(3), sp1, sp2, orb)
        end do
        where (.not.tConverged)
          dd(:,:,0,ii) = (tmp(:,:,2) - tmp(:,:,1)) / (2.0_dp * hh)
        end where
        mm = min(ii, maxcolumns)
        do kk = 1, mm
          ! as central difference derivative formula
          nk = 2 * kk
          where (.not.tConverged)
            dd(:,:,kk,ii) = ((2**(nk)) * dd(:,:,kk-1,ii) - dd(:,:,kk-1,ii-1)) / real(2**nk - 1, dp)
          end where
        end do
        do ll = 1, size(dd, dim=2)
          do kk = 1, size(dd, dim=1)
            if (tConverged(kk, ll)) then
              cycle
            end if
            diff = abs(dd(kk,ll,mm-1,ii) - dd(kk,ll,mm,ii))
            if (diff <= tol * dd(kk,ll,mm,ii)) then
              ! copy converged value over all entries of the table for this element
              dd(kk,ll,:,:) = dd(kk,ll,mm,ii)
              ! and mark it as converged
              tConverged(kk,ll) = .true.
            end if
          end do
        end do
      end do
      deriv(:,:,iCart) = dd(:,:,mm,ii)
    end do

  end subroutine getFirstDerivRichardson

  !> Contains code to calculate the numerical second derivative of a diatomic block of the H0
  !> Hamiltonian and overlap.
  subroutine getSecondDerivFiniteDiff(deriv, skCont, coords, species, atomI, atomJ, orb, deltaXDiff)
    real(dp), intent(out) :: deriv(:,:,:,:)
    type(TSlakoCont), intent(in) :: skCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: atomI, atomJ
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: deltaXDiff

    real(dp) :: interSk(getMIntegrals(skCont))   ! interpolated integs.
    real(dp) :: vect(3), dist
    integer :: ii, jj, kk, ll
    integer :: sp1, sp2
    real(dp) :: tmp(size(deriv, dim=1), size(deriv, dim=2))

    ! second derivative weights for d2 F /dx2
    real(dp), parameter :: stencil(-1:1) = [1.0_dp, -2.0_dp, 1.0_dp]

    @:ASSERT(size(deriv, dim=3) == 3)
    @:ASSERT(size(deriv, dim=4) == 3)

    deriv(:,:,:,:) = 0.0_dp

    sp1 = species(atomI)
    sp2 = species(atomJ)

    do ii = 1, 3
      do jj = 1, 3

        if (ii == jj) then
          ! treat diagonal case separately
          cycle
        end if

        do kk = -1, 1, 2
          do ll = -1, 1, 2

            vect(:) = coords(:,atomJ) - coords(:,atomI)
            vect(ii) = vect(ii) + real(kk, dp) * deltaXDiff
            vect(jj) = vect(jj) + real(ll, dp) * deltaXDiff

            dist = sqrt(sum(vect**2))
            vect(:) = vect / dist

            call getSKIntegrals(skCont, interSk, dist, sp1, sp2)
            call rotateH0(tmp,interSk,vect(1),vect(2),vect(3),sp1,sp2,orb)
            deriv(:,:,jj,ii) = deriv(:,:,jj,ii) + real(kk * ll, dp) * tmp
          end do
        end do
      end do

      do jj = -1, 1

        vect(:) = coords(:,atomJ) - coords(:,atomI)
        vect(ii) = vect(ii) + real(jj, dp) * deltaXDiff

        dist = sqrt(sum(vect**2))
        vect(:) = vect / dist

        call getSKIntegrals(skCont, interSk, dist, sp1, sp2)
        call rotateH0(tmp,interSk,vect(1),vect(2),vect(3),sp1,sp2,orb)
        deriv(:,:,ii,ii) = deriv(:,:,ii,ii) + stencil(jj) * tmp

      end do

    end do

    deriv = deriv / deltaXDiff**2

    ! Off diagonal should be scaled by 1/4
    do ii = 1, 3
      do jj = ii + 1, 3
        deriv(:,:,ii,jj) = 0.25_dp * deriv(:,:,ii,jj)
        deriv(:,:,jj,ii) = 0.25_dp * deriv(:,:,jj,ii)
      end do
    end do

  end subroutine getSecondDerivFiniteDiff

end module dftbp_nonscc
