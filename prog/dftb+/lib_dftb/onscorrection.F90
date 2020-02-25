!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module dftbp_onsitecorrection
  use dftbp_accuracy
  use dftbp_assert
  use dftbp_commontypes
  use dftbp_message
  use dftbp_nonscc, only : TNonSccDiff
  use dftbp_slakocont
  implicit none
  private

  public :: addOnsShift, getEOns, getOnsME
  public :: ons_getOrbitalEquiv, ons_blockIndx, onsBlock_reduce, onsBlock_expand

contains

  !> Add the block shift due to onsite matrix element contributions
  subroutine addOnsShift(potential, iPotential, qBlock, qiBlock, q0, onsMEs, species, orb)

    !> resulting onsite matrix elements
    real(dp), intent(inout) :: potential(:,:,:,:)

    !> resulting onsite matrix elements (imaginary part)
    real(dp), allocatable, intent(inout) :: iPotential(:,:,:,:)

    !> Block charges
    real(dp), intent(in) :: qBlock(:,:,:,:)

    !> Block charges (imaginary part)
    real(dp), intent(in), allocatable :: qiBlock(:,:,:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in) :: onsMEs(:,:,:,:)

    !> species of each atom
    integer, intent(in) :: species(:)

    !> Information about the orbitals in the system
    type(TOrbitals), intent(in) :: orb

    integer :: iAt, nAt, iSp, iSpin, nSpin, ud, iSh, jSh, iOrb, nOrb
    real(dp), allocatable :: tmpME(:,:,:), tmpBlock(:,:)
    real(dp) :: qSumL, degeneracy

    ! factors for q and the 3 possible m channels
    real(dp), parameter :: factor(4) = [1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp]

    nAt = size(potential, dim=3)
    nSpin = size(potential, dim=4)
    allocate(tmpME(orb%mOrb,orb%mOrb,2))
    allocate(tmpBlock(orb%mOrb,orb%mOrb))

    do iAt = 1, nAt
      iSp = species(iAt)
      nOrb = orb%nOrbAtom(iAt)

      call getOnsME(orb, iSp, onsMEs, nOrb, tmpME)

      do iSpin = 1, nSpin
        tmpBlock(:,:) = 0.0_dp
        ! extract the relevant charge parts
        tmpBlock(:nOrb, :nOrb) = qBlock(:nOrb, :nOrb, iAt, iSpin)
        ! diagonal adjusted by reference charges
        do iOrb = 1, nOrb
          tmpBlock(iOrb, iOrb) = tmpBlock(iOrb, iOrb) - q0(iOrb, iAt, iSpin)
        end do

        ! (lambda_ss \pm lambda_st) Delta P^\pm, note that ss' is already zero
        tmpBlock(:nOrb,:nOrb) = tmpBlock(:nOrb,:nOrb) *&
            & (tmpME(:nOrb,:nOrb,1) + factor(iSpin)*tmpME(:nOrb,:nOrb,2))

        ! rotational invariance corection for diagonal part
        do iSh = 1, orb%nShell(iSp)
          degeneracy = real(2*orb%angShell(iSh, iSp) + 1, dp)
          qSumL = 0.0_dp
          do iOrb = orb%posShell(iSh, iSp), orb%posShell(iSh + 1, iSp) - 1
            qSumL = qSumL + tmpBlock(iOrb,iOrb)
          end do
          qSumL = qSumL / degeneracy
          do iOrb = orb%posShell(iSh, iSp), orb%posShell(iSh + 1, iSp) - 1
            tmpBlock(iOrb,iOrb) = tmpBlock(iOrb,iOrb) - qSumL
          end do
        end do

        potential(:nOrb,:nOrb,iAt,iSpin) = potential(:nOrb,:nOrb,iAt,iSpin)&
            & + tmpBlock(:nOrb,:nOrb)

      end do

      if (allocated(qiBlock)) then
        do iSpin = 1, nSpin
          tmpBlock(:,:) = 0.0_dp
          ! extract the relevant charge parts
          tmpBlock(:nOrb, :nOrb) = qiBlock(:nOrb, :nOrb, iAt, iSpin)
          ! diagonal should be zero, as skew symmetric, but just in case
          do iOrb = 1, nOrb
            tmpBlock(iOrb, iOrb) = 0.0_dp
          end do

          ! (lambda_ss \pm lambda_st) Delta P^\pm
          tmpBlock(:nOrb,:nOrb) = tmpBlock(:nOrb,:nOrb) *&
              & (tmpME(:nOrb,:nOrb,1) + factor(iSpin)*tmpME(:nOrb,:nOrb,2))

          iPotential(:nOrb,:nOrb,iAt,iSpin) = iPotential(:nOrb,:nOrb,iAt,iSpin)&
              & + tmpBlock(:nOrb,:nOrb)

        end do
      end if

    end do

  end subroutine addOnsShift


  !> Evaluate an atomic block of on-site correction matrix elements
  subroutine getOnsME(orb, iSp, onsMEs, nOrb, blockME)

    !> Information about the orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> species of atom
    integer, intent(in) :: iSp

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in) :: onsMEs(:,:,:,:)

    !> orbital information
    integer, intent(in) :: nOrb

    !> resulting block of elements
    real(dp), intent(out) :: blockME(:,:,:)

    integer :: ud, iSh, jSh

    blockME(:,:,:) = 0.0_dp
    do ud = 1,2
      ! loop running over same spin and different spin

      do iSh = 1, orb%nShell(iSp)
        do jSh = 1, orb%nShell(iSp)

          if (iSh == jSh .and. orb%angShell(jSh, iSp) == 0) then
            ! ss' on same shell
            cycle
          end if

          blockME(orb%posShell(jSh, iSp) : orb%posShell(jSh + 1, iSp) - 1,&
              & orb%posShell(iSh, iSp) : orb%posShell(iSh + 1, iSp) - 1, ud) =&
              & onsMEs(jSh, iSh, ud, iSp)

        end do
      end do

      ! make symmetric just in case
      blockME(:nOrb, :nOrb, ud) = 0.5_dp * (blockME(:nOrb, :nOrb, ud)&
          & + transpose(blockME(:nOrb, :nOrb, ud)))
    end do

  end subroutine getOnsME

  !> get the onsite energy correction
  subroutine getEons(Eons, qBlock, qiBlock, q0, onsMEs, species, orb)

    !> Onsite energy correction
    real(dp), intent(out) :: Eons(:)

    !> Block charges
    real(dp), intent(in) :: qBlock(:,:,:,:)

    !> Block charges imaginary part
    real(dp), intent(in), allocatable :: qiBlock(:,:,:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> onsite matrix elements for shells (diagonals are xx' elements and ss' are ignored)
    real(dp), intent(in) :: onsMEs(:,:,:,:)

    !> species of each atom
    integer, intent(in) :: species(:)

    !> Information about the orbitals in the system
    type(TOrbitals), intent(in) :: orb

    real(dp), allocatable :: shift(:,:,:,:), iShift(:,:,:,:)

    allocate(shift(orb%mOrb, orb%mOrb, size(qBlock, dim=3), size(qBlock, dim=4)))
    shift(:,:,:,:) = 0.0_dp
    if (allocated(qiBlock)) then
      allocate(iShift(orb%mOrb, orb%mOrb, size(qBlock, dim=3), size(qBlock, dim=4)))
      iShift(:,:,:,:) = 0.0_dp
    end if
    call addOnsShift(shift, iShift, qBlock, qiBlock, q0, onsMEs, species, orb)
    Eons(:) = 0.5_dp*sum(sum(sum(shift(:,:,:,:)*qBlock(:,:,:,:),dim=1),dim=1),dim=2)
    if (allocated(qiBlock)) then
      Eons(:) = Eons(:) + 0.5_dp*sum(sum(sum(iShift(:,:,:,:)*qiBlock(:,:,:,:),dim=1),dim=1),dim=2)
    end if

  end subroutine getEons


  !> Returns the equivalence between the orbitals in the onsite correction
  subroutine ons_getOrbitalEquiv(equiv, orb, species)

    !> The equivalence vector on return
    integer, intent(out) :: equiv(:,:,:)

    !> Information about the orbitals and their angular momenta
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer, intent(in) :: species(:)

    integer :: nAtom, iCount, iSpin, nSpin
    integer :: iAt, iSp, ii, jj, kk, iStart, iEnd

    nAtom = size(equiv, dim=2)
    nSpin = size(equiv, dim=3)

    @:ASSERT(size(equiv, dim=1) == orb%mOrb)

    equiv(:,:,:) = 0

    iCount = 0
    ! set blocks to be full of unique orbitals
    do iSpin = 1, nSpin
      do iAt = 1, nAtom
        iSp = species(iAt)
        do ii = 1, orb%nOrbSpecies(iSp)
          iCount = iCount + 1
          equiv(ii, iAt, iSpin) = iCount
        end do
      end do
    end do

  end subroutine ons_getOrbitalEquiv


  !> Returns the index for packing the relevant parts of atomic blocks into a 1D array
  subroutine ons_blockIndx(iEqBlock, iEqBlockLS, count, orb)

    !> The mapping array on return
    integer, intent(out) :: iEqBlock(:,:,:,:)

    !> The equivalence vector for imagninary parts on return
    integer, intent(inout), allocatable :: iEqBlockLS(:,:,:,:)

    !> Number of prior entries in 1D array holding regular charges
    integer, intent(in) :: count

    !> Information about the orbitals and their angular momenta
    type(TOrbitals), intent(in) :: orb

    integer :: nAtom, nSpin, iCount
    integer :: iAt, iSp, iSpecies
    integer :: iStart1, iEnd1, iStart2, iEnd2
    integer :: ii, jj, kk, ll, ik

    nAtom = size(iEqBlock, dim=3)
    nSpin = size(iEqBlock, dim=4)

    @:ASSERT(size(iEqBlock, dim=1) == orb%mOrb)
    @:ASSERT(size(iEqBlock, dim=2) == orb%mOrb)

    iEqBlock(:,:,:,:) = 0

    iCount = count
    do iSp = 1, nSpin
      do iAt = 1, nAtom
        do ii = 1, orb%nOrbAtom(iAt)
          do jj = ii+1, orb%nOrbAtom(iAt)
            iCount = iCount + 1
            iEqBlock(jj,ii,iAt,iSp) = iCount
          end do
        end do
      end do
    end do

    if (allocated(iEqBlockLS)) then
      iEqBlockLS(:,:,:,:) = 0
      do iSp = 1, nSpin
        do iAt = 1, nAtom
          do ii = 1, orb%nOrbAtom(iAt)
            do jj = ii+1, orb%nOrbAtom(iAt)
              iCount = iCount + 1
              iEqBlockLS(jj,ii,iAt,iSp) = iCount
            end do
          end do
        end do
      end do
    end if

  end subroutine ons_blockIndx


  !> Adds blocks onto end of a 1D vector
  subroutine onsBlock_reduce(input, equiv, orb, output, skew)

    !> unpacked data
    real(dp), intent(in) :: input(:,:,:,:)

    !> equivalences
    integer, intent(in) :: equiv(:,:,:,:)

    !> Information about the orbitals and their angular momenta
    type(TOrbitals), intent(in) :: orb

    !> 1D array with appended data
    real(dp), intent(inout) :: output(:)

    !> is skew symmetry required
    logical, optional, intent(in) :: skew

    integer :: nAtom, nSpin
    integer :: iS, iOrb1, iOrb2, iAt
    logical :: iSkew

    nAtom = size(input, dim=3)
    nSpin = size(input, dim=4)
    @:ASSERT(size(input, dim=1) == orb%mOrb)
    @:ASSERT(size(input, dim=2) == orb%mOrb)
    @:ASSERT(all(shape(equiv) == (/ orb%mOrb, orb%mOrb, nAtom, nSpin /)))

    if (present(skew)) then
      iSkew = skew
    else
      iSkew = .false.
    end if

    do iS = 1, nSpin
      do iAt = 1, nAtom
        do iOrb1 = 1, orb%nOrbAtom(iAt)
          do iOrb2 = 1, orb%nOrbAtom(iAt)
            !if (iOrb1 == iOrb2) then
            !  cycle
            !end if
            if (equiv(iOrb1, iOrb2, iAt, iS) > 0) then
              if (iSkew) then
                output(equiv(iOrb1, iOrb2, iAt, iS)) = &
                    & 0.5_dp*( input(iOrb1, iOrb2, iAt, iS) &
                    &  - input(iOrb2, iOrb1, iAt, iS) )
              else
                output(equiv(iOrb1, iOrb2, iAt, iS)) = &
                    & 0.5_dp*( input(iOrb1, iOrb2, iAt, iS) &
                    &  + input(iOrb2, iOrb1, iAt, iS) )
              end if
            end if
          end do
        end do
      end do
    end do

  end subroutine onsBlock_reduce


  !> Extract orbital blocks from the end of a 1D vector
  subroutine onsblock_expand(input, blockEquiv, orb, output, orbEquiv, skew)

    !> 1D array of packed data
    real(dp), intent(in) :: input(:)

    !> equivalences for blocks on atomic sites
    integer, intent(in) :: blockEquiv(:,:,:,:)

    !> Information about the orbitals and their angular momenta
    type(TOrbitals), intent(in) :: orb

    !> unpacked data
    real(dp), intent(out) :: output(:,:,:,:)

    !> equivalences for atoms
    integer, intent(in),optional :: orbEquiv(:,:,:)

    !> is skew symmetry required
    logical, optional, intent(in) :: skew

    integer :: nAtom, nSpin
    integer :: iAt, iSp
    integer :: iStart1, iEnd1, iStart2, iEnd2
    integer :: ii, jj, kk, ll, ik
    logical :: iSkew

    nAtom = size(output, dim=3)
    nSpin = size(output, dim=4)

    if (present(skew)) then
      iSkew = skew
    else
      iSkew = .false.
    end if

    @:ASSERT(size(output, dim=1) == orb%mOrb)
    @:ASSERT(size(output, dim=2) == orb%mOrb)
  #:call ASSERT_CODE
    if (present(orbEquiv)) then
      @:ASSERT(all(shape(orbEquiv) == (/ orb%mOrb, nAtom, nSpin /)))
    end if
  #:endcall ASSERT_CODE
    @:ASSERT(all(shape(blockEquiv) == shape(output)))

    output = 0.0_dp

    do iSp = 1, nSpin
      do iAt = 1, nAtom
        do ii = 1, orb%nOrbAtom(iAt)
          if (present(orbEquiv) .and. .not. iSkew) then
            output(ii,ii,iAt,iSp) = input(orbEquiv(ii,iAt,iSp))
          end if
          do jj = ii + 1, orb%nOrbAtom(iAt)
            output(jj,ii,iAt,iSp) = input(blockEquiv(jj,ii,iAt,iSp))
            if (iSkew) then
              output(ii,jj,iAt,iSp) = -input(blockEquiv(jj,ii,iAt,iSp))
            else
              output(ii,jj,iAt,iSp) = input(blockEquiv(jj,ii,iAt,iSp))
            end if
          end do
        end do
      end do
    end do

  end subroutine onsblock_expand

end module dftbp_onsitecorrection
