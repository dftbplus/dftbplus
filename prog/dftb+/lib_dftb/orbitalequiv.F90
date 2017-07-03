!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains routines to manipulate orbital equivalency relations.
!!* @desc An orbital equivalency relation is a mapping, which maps the orbitals
!!* of the atoms onto a one dimensional vector, where equivalent orbitals are
!!* mapped on the same element in the vector. Two orbitals are equivalent, if
!!* charge can be transferred between the orbitals, without changing the
!!* resulting Hamiltonian or the total energy. The mapping is an
!!* (mmAng, nAtom, nSpin) shaped integer array, where the integer for
!!* (iOrb, iAtom, iSpin) specifies the position in the 1D vector for orbital
!!* iOrb on atom iAtom for spin iSpin. Values must be positive integers and
!!* continuous. Zeros in the mapping vector stand for non-existent orbitals.
module orbitalequiv
  use assert
  use accuracy
  use commontypes
  implicit none
  private

  public :: OrbitalEquiv_merge, OrbitalEquiv_reduce, OrbitalEquiv_expand

contains

  !!* Merges two equivalency arrays by finding the intersection in the
  !!* equivalences.
  !!* @param equiv1 First equivalency array.
  !!* @param equiv2 Second equivalency array
  !!* @param orb Contains information about the atomic orbitals in the system.
  !!* @param equivNew New equivalency array on return.
  !!* @note Equivalences marked below 0 are mapped to 0 in the final
  !!* equivalence matrix
  subroutine OrbitalEquiv_merge(equiv1, equiv2, orb, equivNew)
    integer, intent(in) :: equiv1(:,:,:)
    integer, intent(in) :: equiv2(:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(out) :: equivNew(:,:,:)

    integer :: nAtom, nSpin
    integer :: newInd, iS, iAt, iOrb
    logical, allocatable :: mask(:,:,:), tmpMask(:,:,:)

    nAtom = size(equiv1, dim=2)
    nSpin = size(equiv1, dim=3)

    @:ASSERT(all(shape(equiv1) == shape(equiv2)))
    @:ASSERT(all(shape(equiv1) == shape(equivNew)))
    @:ASSERT(size(equiv1, dim=1) >= orb%mOrb)
    @:ASSERT(all((equiv1 /= 0) .eqv. (equiv2 /=0)))

    allocate(mask(size(equiv1, dim=1), size(equiv1, dim=2), size(equiv1, dim=3)))
    allocate(tmpMask(size(equiv1, dim=1), size(equiv1, dim=2), size(equiv1, dim=3)))

    mask(:,:,:) = (equiv1 > 0 .and. equiv2 > 0) ! True for the elements to be
    !  processed
    equivNew(:,:,:) = 0
    newInd = 1                          ! Position in the reduced vector
    do iS = 1, nSpin
      do iAt = 1, nAtom
        do iOrb = 1, orb%nOrbAtom(iAt)
          !! If element had been already processed, skip
          if (.not. mask(iOrb, iAt, iS)) then
            cycle
          end if
          !! Get mask for orbitals equivalent to (iOrb, iAt, iS) in both arrays
          tmpMask = (equiv1 == equiv1(iOrb, iAt, iS) &
              &.and. equiv2 == equiv2(iOrb, iAt, iS))
          !! If equivalent orbitals exist, map them to reduced vector
          if (any(tmpMask)) then
            where (tmpMask)
              equivNew = newInd
            end where
            newInd = newInd + 1
            mask = mask .and. .not. tmpMask
          end if
        end do
      end do
    end do

  end subroutine OrbitalEquiv_merge

  !!* Reduces passed orbital property by summing up on equivalent orbitals.
  !!* @param input The vector containing the values for all orbitals
  !!* @param equiv Vector containing equivalency relations between the orbitals.
  !!* @param orb Contains information about the atomic orbitals in the system
  !!* @param output Contains the equivalent orbital sums at exit as vector.
  !!* @note equivalences of 0 or below are not reduced
  subroutine OrbitalEquiv_reduce(input, equiv, orb, output)
    real(dp), intent(in) :: input(:,:,:)
    integer, intent(in) :: equiv(:,:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(out) :: output(:)

    integer :: nAtom, nSpin
    integer :: iS, iOrb, iAt

    nAtom = size(input, dim=2)
    nSpin = size(input, dim=3)

    @:ASSERT(size(input, dim=1) == orb%mOrb)
    @:ASSERT(all(shape(equiv) == (/ orb%mOrb, nAtom, nSpin /)))
    @:ASSERT(size(output) == maxval(equiv))

    output(:) = 0.0_dp
    do iS = 1, nSpin
      do iAt = 1, nAtom
        do iOrb = 1, orb%nOrbAtom(iAt)
          if (equiv(iOrb, iAt, iS) > 0) then
            output(equiv(iOrb, iAt, iS)) = output(equiv(iOrb, iAt, iS)) &
                & + input(iOrb, iAt, iS)
          end if
        end do
      end do
    end do

  end subroutine OrbitalEquiv_reduce



  !!* Expands a reduced vector by putting every element of it into the first
  !!* corresponding orbital and putting zero for all other equivalent orbitals.
  !!* @param input Reduced vector.
  !!* @param equiv Equivalency relation.
  !!* @param orb Contains information about the atomic orbitals in the system
  !!* @param output Expanded array.
  !!* @note equivalences of 0 or below are not expanded
  subroutine OrbitalEquiv_expand(input, equiv, orb, output)
    real(dp), intent(in) :: input(:)
    integer, intent(in) :: equiv(:,:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(out) :: output(:,:,:)

    integer :: nSpin, nAtom
    integer ::iS, iAt, iOrb
    logical, allocatable :: mask(:)

    nSpin = size(output, dim=3)
    nAtom = size(output, dim=2)

    @:ASSERT(all(shape(equiv) == shape(output)))
    @:ASSERT(maxval(equiv) == size(input))

    allocate(mask(0:size(input)))

    mask(:) = .true.
    output(:,:,:) = 0.0_dp
    do iS = 1, nSpin
      do iAt = 1, nAtom
        do iOrb = 1, orb%nOrbAtom(iAt)
          if (mask(equiv(iOrb, iAt, iS))) then
            if (equiv(iOrb, iAt, iS) > 0) then
              output(iOrb, iAt, iS) = input(equiv(iOrb, iAt, iS))
              mask(equiv(iOrb, iAt, iS)) = .false.
            end if
          end if
        end do
      end do
    end do

  end subroutine OrbitalEquiv_expand


end module orbitalequiv
