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
  public

contains

  subroutine addOnsShift(potential, qBlock, orb, ons_en, species)
    real(dp), intent(inout) :: potential(:,:,:,:)
    real(dp), intent(in) :: qBlock(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: ons_en(:,:)
    integer, intent(in) :: species(:)

    integer :: nOrb, nAtom, nSpin
    integer :: iAt, iSp, iSpin, sigma, ud
    integer :: mu
    real(dp) :: ons_conv(orb%mOrb,orb%mOrb,2)
    real(dp) :: factor

    nAtom = size(potential, dim=3)
    nSpin = size(potential, dim=4)

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
        do mu = 1, nOrb-1
          potential(mu+1:nOrb,mu,iAt,iSpin) = potential(mu+1:nOrb,mu,iAt,iSpin)&
              & + qBlock(mu+1:nOrb,mu,iAt,iSpin)&
              & * ( ons_conv(mu+1:nOrb,mu,1) + factor*ons_conv(mu+1:nOrb,mu,2) )
        end do
        do mu = 1, nOrb
          potential(mu,mu+1:nOrb,iAt,iSpin) = potential(mu+1:nOrb,mu,iAt,iSpin)
        end do
      end do
    end do

  end subroutine addOnsShift


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
    call addOnsShift(shift, qBlock, orb, ons_en, species)

    Eons(:) = 0.5_dp*sum(sum(sum(shift(:,:,:,:)*qBlock(:,:,:,:),dim=1),dim=1),dim=2)

    call getEri(Eri,q,q0,orb,ons_en,species)
    Eons(:) = Eons + Eri

  end subroutine getEons


  subroutine getEri(Eri, q, q0, orb, ons_en, species)
    real(dp), intent(out)          :: Eri(:)
    real(dp), intent(in)           :: q(:,:,:), q0(:,:,:)
    type(TOrbitals), intent(in)    :: orb
    real(dp), intent(in)           :: ons_en(:,:)
    integer, intent(in)            :: species(:)

    integer  :: nOrb, nAtom, nSpin
    integer  :: iAt, iSp, iSpin
    integer  :: mu, nu
    real(dp) :: factor, fact, onethird, onefifth
    real(dp) :: qDiff(orb%mOrb,size(q, dim=2),size(q, dim=3))

    nAtom = size(q, dim=2)
    nSpin = size(q, dim=3)

    onethird = 1.0_dp/3.0_dp
    onefifth = 1.0_dp/5.0_dp

    Eri = 0.0_dp
    qDiff = q - q0

    do iSpin = 1,nSpin
      factor = -1.0_dp
      if (iSpin == 1) factor = 1.0_dp
      do iAt = 1, nAtom
        nOrb = orb%nOrbAtom(iAt)
        iSp  = species(iAt)
        if (nOrb > 1) then
          do mu = 2, 4
            do nu = mu, 4
              fact = -1.0_dp
              if (mu == nu) fact = 1.0_dp
              Eri(iAt) = Eri(iAt) + fact*onethird &
                  & *qDiff(mu,iAt,iSpin)*qDiff(nu,iAt,iSpin)&
                  & *( ons_en(iSp,3) + factor*ons_en(iSp,4) )
            end do
          end do
        end if
        !d-orbital (will be optimized)
        if (nOrb > 4) then
          do mu = 5, nOrb
            do nu = mu, nOrb
              fact = -1.0_dp
              if (mu == nu) fact = 2.0_dp
              Eri(iAt) = Eri(iAt) + fact*onefifth &
                  & *qDiff(mu,iAt,iSpin)*qDiff(nu,iAt,iSpin)&
                  & *( ons_en(iSp,9) + factor*ons_en(iSp,10) )
            end do
          end do
        end if
      end do
    end do

  end subroutine getEri


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


  subroutine Ons_blockIndx(iEqBlockOns, count, orb)
    integer, intent(out)        :: iEqBlockOns(:,:,:,:)
    integer, intent(in)         :: count
    type(TOrbitals), intent(in) :: orb

    integer :: nAtom, nSpin, nOrb, iCount
    integer :: iAt, iSp, iOrb1, iOrb2

    nAtom = size(iEqBlockOns, dim=3)
    nSpin = size(iEqBlockOns, dim=4)
    @:ASSERT(nSpin == 1 .or. nSpin == 2)

    @:ASSERT(size(iEqBlockOns, dim=1) == orb%mOrb)
    @:ASSERT(size(iEqBlockOns, dim=2) == orb%mOrb)

    iEqBlockOns(:,:,:,:) = 0

    iCount = count
    do iSp = 1, nSpin
      do iAt = 1, nAtom
        nOrb = orb%nOrbAtom(iAt)
        do iOrb1 = 1, nOrb - 1
          do iOrb2 = iOrb1 + 1, nOrb
            iCount = iCount + 1
            iEqBlockOns(iOrb1, iOrb2, iAt, iSp) = iCount
          end do
        end do
      end do
    end do

  end subroutine Ons_blockIndx


  subroutine Onsblock_reduce(input, equiv, orb, output)
    real(dp), intent(in) :: input(:,:,:,:)
    integer, intent(in) :: equiv(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(inout) :: output(:)

    integer :: nAtom, nSpin, nOrb
    integer :: iS, iOrb1, iOrb2, iAt
    integer :: pos

    nAtom = size(input, dim=3)
    nSpin = size(input, dim=4)
    @:ASSERT(size(input, dim=1) == orb%mOrb)
    @:ASSERT(size(input, dim=2) == orb%mOrb)
    @:ASSERT(all(shape(equiv) == (/ orb%mOrb, orb%mOrb, nAtom, nSpin /)))
    @:ASSERT(nSpin == 1 .or. nSpin == 2)

    do iS = 1, nSpin
      do iAt = 1, nAtom
        nOrb = orb%nOrbAtom(iAt)
        do iOrb1 = 1, nOrb
          do iOrb2 = 1, nOrb
            pos = equiv(iOrb1, iOrb2, iAt, iS)
            if (pos > 0) then
              output(pos) = 0.5_dp*( input(iOrb1, iOrb2, iAt, iS) &
                                &  + input(iOrb2, iOrb1, iAt, iS) )
            end if
          end do
        end do
      end do
    end do

  end subroutine OnsBlock_reduce


  subroutine Onsblock_expand(input, equiv, orb, output)
    real(dp), intent(in) :: input(:)
    integer, intent(in) :: equiv(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(inout) :: output(:,:,:,:)

    integer :: nAtom, nSpin, nOrb
    integer :: iS, iOrb1, iOrb2, iAt
    logical :: mask(size(input))
    integer :: pos

    nAtom = size(output, dim=3)
    nSpin = size(output, dim=4)
    @:ASSERT(size(output, dim=1) == orb%mOrb)
    @:ASSERT(size(output, dim=2) == orb%mOrb)
    @:ASSERT(all(shape(equiv) == (/ orb%mOrb, orb%mOrb, nAtom, nSpin /)))
    @:ASSERT(nSpin == 1 .or. nSpin == 2)

    mask(:) = .true.
    do iS = 1, nSpin
      do iAt = 1, nAtom
        nOrb = orb%nOrbAtom(iAt)
        do iOrb1 = 1, nOrb
          do iOrb2 = 1, nOrb
            pos = equiv(iOrb1, iOrb2, iAt, iS)
            if (pos > 0) then
              if (mask(pos)) then
                output(iOrb1, iOrb2, iAt, iS) = input(pos)
                output(iOrb2, iOrb1, iAt, iS) = input(pos)
                mask(pos) = .false.
              end if
            end if
          end do
        end do
      end do
    end do

  end subroutine Onsblock_expand

end module onsitecorrection
