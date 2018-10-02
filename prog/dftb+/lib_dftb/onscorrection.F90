!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module onsitecorrection
  use accuracy
  use assert
  use commontypes
  use message
  use nonscc, only : NonSccDiff
  use slakocont
  implicit none
  private

  public :: addOnsShift, getEOns
  public :: ons_getOrbitalEquiv, ons_blockIndx, onsBlock_reduce, onsBlock_expand

contains

  subroutine addOnsShift(potential, qBlock, orb, ons_en, species, q0)
    real(dp), intent(inout) :: potential(:,:,:,:)
    real(dp), intent(in) :: qBlock(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: ons_en(:,:)
    integer, intent(in) :: species(:)
    real(dp), intent(in) :: q0(:,:,:)

    integer :: nOrb, nAtom, nSpin
    integer :: iAt, iSp, iSpin, sigma, ud
    integer :: mu
    real(dp) :: ons_conv(orb%mOrb,orb%mOrb,2)
    real(dp), allocatable :: tmpBlock(:,:), tmpQ(:,:,:)
    real(dp) :: factor

    nAtom = size(potential, dim=3)
    nSpin = size(potential, dim=4)
    allocate(tmpBlock(orb%mOrb,orb%mOrb))
    allocate(tmpQ(orb%mOrb, nAtom, nSpin))

    tmpBlock(:,:) = 0.0_dp
    tmpQ(:,:,:) = 0.0_dp
    do iSpin = 1,nSpin
      if (iSpin == 1) then
        factor = 1.0_dp
      else
        factor = -1.0_dp
      end if
      do iAt = 1, nAtom
        nOrb = orb%nOrbAtom(iAt)
        iSp  = species(iAt)
        call getOnsME(orb, iSp, ons_en, ons_conv)
        tmpBlock(:nOrb,:nOrb) = qBlock(:nOrb,:nOrb,iAt,iSpin)
        do mu = 1, nOrb
          tmpQ(mu, iAt, iSpin) = tmpBlock(mu,mu)
          tmpBlock(mu,mu) = 0.0_dp
        end do
        potential(:nOrb,:nOrb,iAt,iSpin) = potential(:nOrb,:nOrb,iAt,iSpin)&
            & + tmpBlock(:nOrb,:nOrb)*(ons_conv(:nOrb,:nOrb,1)+factor*ons_conv(:nOrb,:nOrb,2))
      end do
    end do

    ! rotational invariant correction
    call addRIShift(potential, tmpQ, q0, orb, ons_en, species)

  end subroutine addOnsShift


  ! Note, d orbitals were missing in original code, not added back yet:
  subroutine addRIshift(potential, q, q0, orb, ons_en, species)
    real(dp), intent(inout)        :: potential(:,:,:,:)
    real(dp), intent(in)           :: q(:,:,:), q0(:,:,:)
    type(TOrbitals), intent(in)    :: orb
    real(dp), intent(in)           :: ons_en(:,:)
    integer, intent(in)            :: species(:)

    integer :: nOrb, nAtom, nSpin
    integer :: iAt, iSpin, iSp, iOrb
    real(dp) :: qDiff( orb%mOrb, size(potential, dim=3), size(potential, dim=4))
    real(dp) :: factor
    real(dp), parameter :: onethird = 1.0_dp/3.0_dp

    nAtom = size(potential, dim=3)
    nSpin = size(potential, dim=4)

    do iOrb = 1, orb%mOrb
      qDiff(iOrb,:,:) = sum(q0(2:4,:,:),dim=1) - sum(q(2:4,:,:),dim=1)
    end do
    qDiff(2:4,:,:) = qDiff(2:4,:,:) + 3.0_dp * ( q(2:4,:,:) - q0(2:4,:,:) )

    do iSpin = 1,nSpin
      if (iSpin == 1) then
        factor = 1.0_dp
      else
        factor = -1.0_dp
      end if
      do iAt = 1, nAtom
        nOrb = orb%nOrbAtom(iAt)
        iSp  = species(iAt)
        if (nOrb > 1) then
          do iOrb = 2, 4
            potential(iOrb,iOrb,iAt,iSpin) = potential(iOrb,iOrb,iAt,iSpin)&
                & + onethird * qDiff(iOrb,iAt,iSpin) * ( ons_en(iSp,3) + factor*ons_en(iSp,4) )
          end do
        end if
      end do
    end do

  end subroutine addRIshift


  subroutine getEons(Eons,qBlock,q,q0,orb,ons_en,species)
    real(dp), intent(out)       :: Eons(:)
    real(dp), intent(in)        :: q(:,:,:), q0(:,:,:)
    real(dp), intent(in)        :: qBlock(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in)        :: ons_en(:,:)
    integer, intent(in)         :: species(:)

    integer :: nOrb, nAtom, nSpin
    integer :: iAt, iSp, iSpin
    integer :: mu, nu
    real(dp) :: ons_conv(orb%mOrb,orb%mOrb,2)
    real(dp) :: factor
    real(dp) :: Eri(size(Eons))

    real(dp), allocatable :: shift(:,:,:,:)

    nAtom = size(qBlock, dim=3)
    nSpin = size(qBlock, dim=4)

    Eons(:) = 0.0_dp

    allocate(shift(orb%mOrb, orb%mOrb, nAtom, nSpin))

    shift(:,:,:,:) = 0.0_dp
    call addOnsShift(shift, qBlock, orb, ons_en, species, q0)

    Eons(:) = 0.5_dp*sum(sum(sum(shift(:,:,:,:)*qBlock(:,:,:,:),dim=1),dim=1),dim=2)

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
  subroutine ons_blockIndx(iEqBlock, count, orb)

    !> The mapping array on return
    integer, intent(out) :: iEqBlock(:,:,:,:)

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


  ! Internal routines

    subroutine getOnsME(orb, iSp, ons_en, ons_conv)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iSp
    real(dp), intent(in) :: ons_en(:,:)
    real(dp), intent(out) :: ons_conv(:,:,:)

    integer :: iOrb1, iOrb2, nOrb, ud

    ons_conv(:,:,:) = 0.0_dp

    nOrb = orb%nOrbSpecies(iSp)
    if (nOrb > 9) then
      call error("Onsite correction does not work for atoms containing 'f' orbitals yet.")
    end if

    do iOrb1 =1, nOrb-1
      do iOrb2 = iOrb1+1, nOrb
        do ud = 1,2

          if (iOrb1 == 1) then
            ons_conv(iOrb1,iOrb2,ud) = ons_en(iSp,ud)
            if (iOrb2 > iOrb1 + 3) then
              ons_conv(iOrb1,iOrb2,ud)= ons_en(iSp,4+ud)
            endif
            !__________d orbital_______________
          elseif (iOrb1 > 4) then
            ons_conv(iOrb1,iOrb2,ud) = ons_en(iSp,8+ud)
            !__________________________________
          else
            ons_conv(iOrb1,iOrb2,ud) = ons_en(iSp,2+ud)
            !for d orbitals:
            if (iOrb2 > 4) then
              ons_conv(iOrb1,iOrb2,ud) = ons_en(iSp,6+ud)
            end if
          end if

          ons_conv(iOrb2,iOrb1,ud) = ons_conv(iOrb1,iOrb2,ud)

        end do
      end do
    end do

  end subroutine getOnsME

end module onsitecorrection
