!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains code to perform the sk rotations of matrix elements from the parameterization
!> orientation along <0,0,1> to the one needed in the actual calculation.
!>
!> To do: Transformations to give the derivatives with respect to ll, mm and nn.
!> Base on "Compact expression for the angular dependence of tight-binding hamiltonian matrix
!> elements", A. V. Podolskiy and P. Vogl, Phys. Rev.  B 69, 233101 (2004).
!>
!> Caveat: Only angular momenta up to f are currently allowed
module sk
  use assert
  use accuracy
  use commontypes
  implicit none

  private

  public :: rotateH0


  !> Maximal angular momentum, for which rotations are present
  integer, parameter :: mAngRot_ = 3

contains


  !> Driver for making the non-SCC hhamiltonian or overlap matrices for a given diatomic block
  !> Caveat: Only angular momenta up to f are currently allowed
  subroutine rotateH0(hh, skIntegs, ll, mm, nn, iSp1, iSp2, orb)

    !> the rectangular matrix containing the resulting diatomic matrix elements
    real(dp), intent(out) :: hh(:,:)

    !> Slater-Koster table for dimer of species i-j
    real(dp), intent(in), target :: skIntegs(:)

    !> directional cosine ll
    real(dp), intent(in) :: ll

    !> directional cosine mm
    real(dp), intent(in) :: mm

    !> directional cosine nn
    real(dp), intent(in) :: nn

    !> Chemical species of atom i
    integer, intent(in) :: iSp1

    !> chemical species of atom j
    integer, intent(in) :: iSp2

    !> Information about the orbitals of chemical species in the system.
    type(TOrbitals), intent(in) :: orb

    integer :: iCol, iRow, ind, iSh1, iSh2
    integer :: ang1, ang2, nOrb1, nOrb2
    real(dp), pointer :: pSK(:)
    real(dp) :: tmpH(2*mAngRot_+1,2*mAngRot_+1)

    @:ASSERT(maxval(orb%angShell) <=  mAngRot_)
    @:ASSERT(all(shape(hh) >= (/ orb%nOrbSpecies(iSp1), orb%nOrbSpecies(iSp2) /)))

    hh(:,:) = 0.0_dp
    ind = 1
    iCol = 1
    do iSh1 = 1, orb%nShell(iSp1)
      ang1 = orb%angShell(iSh1, iSp1)
      nOrb1 = 2 * ang1 + 1
      iRow = 1
      do iSh2 = 1, orb%nShell(iSp2)
        ang2 = orb%angShell(iSh2, iSp2)
        nOrb2 = 2 * ang2 + 1
        @:ASSERT(size(skIntegs) >= ind + min(ang1,ang2))
        pSK => skIntegs(ind:ind+min(ang1,ang2))
        select case (ang1)
        case (0)
          select case (ang2)
          case (0)
            call ss(tmpH,pSK)
          case (1)
            call sp(tmpH,ll,mm,nn,pSK)
          case (2)
            call sd(tmpH,ll,mm,nn,pSK)
          case (3)
            call sf(tmpH,ll,mm,nn,pSK)
          end select
        case (1)
          select case (ang2)
          case (0)
            call sp(tmpH,ll,mm,nn,pSK)
          case (1)
            call pp(tmpH,ll,mm,nn,pSK)
          case (2)
            call pd(tmpH,ll,mm,nn,pSK)
          case (3)
            call pf(tmpH,ll,mm,nn,pSK)
          end select
        case (2)
          select case (ang2)
          case(0)
            call sd(tmpH,ll,mm,nn,pSK)
          case(1)
            call pd(tmpH,ll,mm,nn,pSK)
          case (2)
            call dd(tmpH,ll,mm,nn,pSK)
          case (3)
            call df(tmpH,ll,mm,nn,pSK)
          end select
        case (3)
          select case (ang2)
          case(0)
            call sf(tmpH,ll,mm,nn,pSK)
          case(1)
            call pf(tmpH,ll,mm,nn,pSK)
          case(2)
            call df(tmpH,ll,mm,nn,pSK)
          case(3)
            call ff(tmpH,ll,mm,nn,pSK)
          end select
        end select

        if (ang1 <= ang2) then
          hh(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1) = tmpH(1:nOrb2,1:nOrb1)
        else
          hh(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1) = (-1.0_dp)**(ang1+ang2) &
              &* transpose(tmpH(1:nOrb1,1:nOrb2))
        end if
        ind = ind + min(ang1,ang2) + 1
        iRow = iRow + nOrb2
      end do
      iCol = iCol + nOrb1
    end do

  end subroutine rotateH0


  !> rotation routine for interaction of an s orbital with an s orbital
  subroutine ss(hh, sk)

    !> dimeric block to put the results in to
    real(dp), intent(inout) :: hh(:,:)

    !> Slater-Koster table for dimer element of the Slater-Koster table
    real(dp), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 1)
    @:ASSERT(all(shape(hh) >= (/ 1, 1 /)))

    hh(1,1) = sk(1)

  end subroutine ss


  !> rotation routine for interaction of an s orbital with a p orbital
  subroutine sp(hh, ll, mm, nn, sk)

    !> dimeric block to put the results in to
    real(dp), intent(inout) :: hh(:,:)

    !> directional cosine ll
    real(dp), intent(in) :: ll

    !> directional cosine mm
    real(dp), intent(in) :: mm

    !> directional cosine nn
    real(dp), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    real(dp), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 1)
    @:ASSERT(all(shape(hh) >= (/ 3, 1 /)))

    hh(1,1) = mm*sk(1)
    hh(2,1) = nn*sk(1)
    hh(3,1) = ll*sk(1)

  end subroutine sp


  !> rotation routine for interaction of an s orbital with a d orbital
  subroutine sd(hh, ll, mm, nn, sk)

    !> dimeric block to put the results in to
    real(dp), intent(inout) :: hh(:,:)

    !> directional cosine ll
    real(dp), intent(in) :: ll

    !> directional cosine mm
    real(dp), intent(in) :: mm

    !> directional cosine nn
    real(dp), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    real(dp), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 1)
    @:ASSERT(all(shape(hh) >= (/ 5, 1 /)))

    hh(1,1) = ll*mm*sqrt(3.0_dp)*sk(1)
    hh(2,1) = mm*sqrt(3.0_dp)*nn*sk(1)
    hh(3,1) = (3.0_dp/ 2.0_dp*nn**2- 1.0_dp/ 2.0_dp)*sk(1)
    hh(4,1) = ll*sqrt(3.0_dp)*nn*sk(1)
    hh(5,1) = (2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(3.0_dp)*sk(1)/2.0_dp

  end subroutine sd


  !> rotation routine for interaction of an s orbital with an f orbital
  subroutine sf(hh, ll, mm, nn, sk)

    !> dimeric block to put the results in to
    real(dp), intent(inout) :: hh(:,:)

    !> directional cosine ll
    real(dp), intent(in) :: ll

    !> directional cosine mm
    real(dp), intent(in) :: mm

    !> directional cosine nn
    real(dp), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    real(dp), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 1)
    @:ASSERT(all(shape(hh) >= (/ 1, 7 /)))

    hh(1,1) = sqrt(2.0_dp)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(5.0_dp)&
        &*sk(1)/ 4.0_dp
    hh(2,1) = ll*mm*sqrt(15.0_dp)*nn*sk(1)
    hh(3,1) = sqrt(2.0_dp)*mm*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)/ 4.0_dp
    hh(4,1) = (nn*(5.0_dp*nn**2-3.0_dp)*sk(1))/ 2.0_dp
    hh(5,1) = sqrt(2.0_dp)*ll*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)/4.0_dp
    hh(6,1) = (2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(15.0_dp)*nn*sk(1)/ 2.0_dp
    hh(7,1) = sqrt(2.0_dp)*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp

  end subroutine sf


  !> rotation routine for interaction of a p orbital with a p orbital
  subroutine pp(hh, ll, mm, nn, sk)

    !> dimeric block to put the results in to
    real(dp), intent(inout) :: hh(:,:)

    !> directional cosine ll
    real(dp), intent(in) :: ll

    !> directional cosine mm
    real(dp), intent(in) :: mm

    !> directional cosine nn
    real(dp), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    real(dp), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 2)
    @:ASSERT(all(shape(hh) >= (/ 3, 3 /)))

    hh(1,1) = (1.0_dp-nn**2-ll**2)*sk(1)+(nn**2+ll**2)*sk(2)
    hh(2,1) = nn*mm*sk(1)-nn*mm*sk(2)
    hh(3,1) = ll*mm*sk(1)-ll*mm*sk(2)
    hh(1,2) = hh(2,1)
    hh(2,2) = nn**2*sk(1)+(1.0_dp-nn**2)*sk(2)
    hh(3,2) = nn*ll*sk(1)-nn*ll*sk(2)
    hh(1,3) = hh(3,1)
    hh(2,3) = hh(3,2)
    hh(3,3) = ll**2*sk(1)+(1.0_dp-ll**2)*sk(2)

  end subroutine pp


  !> rotation routine for interaction of a p orbital with a d orbital
  subroutine pd(hh, ll, mm, nn, sk)

    !> dimeric block to put the results in to
    real(dp), intent(inout) :: hh(:,:)

    !> directional cosine ll
    real(dp), intent(in) :: ll

    !> directional cosine mm
    real(dp), intent(in) :: mm

    !> directional cosine nn
    real(dp), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    real(dp), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 2)
    @:ASSERT(all(shape(hh) >= (/ 3, 5 /)))

    hh(1,1) = -(-1.0_dp+nn**2+ll**2)*ll*sqrt(3.0_dp)&
        &*sk(1)+((2.0_dp*nn**2+2.0_dp*ll**2-1.0_dp)*ll*sk(2))
    hh(2,1) = -(-1.0_dp+nn**2+ll**2)*sqrt(3.0_dp)*nn*&
        &sk(1)+((2.0_dp*nn**2+2.0_dp*ll**2-1.0_dp)*nn*sk(2))
    hh(3,1) = mm*(3.0_dp*nn**2-1.0_dp)*sk(1)/2.0_dp&
        &-sqrt(3.0_dp)*(nn**2)*mm*sk(2)
    hh(4,1) = mm*ll*sqrt(3.0_dp)*nn*sk(1)-2.0_dp*ll*mm*nn*sk(2)
    hh(5,1) = mm*(2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(3.0_dp)*sk(1)/2.0_dp&
        &-(nn**2+2.0_dp*ll**2)*mm*sk(2)
    hh(1,2) = ll*mm*nn*sqrt(3.0_dp)*sk(1)-2.0_dp*nn*ll*mm*sk(2)
    hh(2,2) = mm*(nn**2)*sqrt(3.0_dp)*sk(1)&
        &-(2.0_dp*nn**2-1.0_dp)*mm*sk(2)
    hh(3,2) = (nn*(3.0_dp*nn**2-1.0_dp)*sk(1))/2.0_dp&
        &-nn*sqrt(3.0_dp)*(-1.0_dp+nn**2)*sk(2)
    hh(4,2) = ll*nn**2*sqrt(3.0_dp)*sk(1)-(2.0_dp*nn**2-1.0_dp)*ll*&
        &sk(2)
    hh(5,2) = (2.0_dp*ll**2-1.0_dp+nn**2)*nn*sqrt(3.0_dp)*sk(1)/2.0_dp&
        &-(nn*(2.0_dp*ll**2-1.0_dp+nn**2)*sk(2))
    hh(1,3) = (ll**2)*mm*sqrt(3.0_dp)*sk(1)&
        &-(2.0_dp*ll**2-1.0_dp)*mm*sk(2)
    hh(2,3) = ll*mm*sqrt(3.0_dp)*nn*sk(1)-2.0_dp*mm*ll*nn*sk(2)
    hh(3,3) = (ll*(3.0_dp*nn**2-1.0_dp)*sk(1))/2.0_dp&
        &-sqrt(3.0_dp)*(nn**2)*ll*sk(2)
    hh(4,3) = ll**2*sqrt(3.0_dp)*nn*sk(1)-(2.0_dp*ll**2-1.0_dp)*nn&
        &*sk(2)
    hh(5,3) = ll*(2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(3.0_dp)*sk(1)/2.0_dp&
        &-((nn**2-2.0_dp+2.0_dp*ll**2)*ll*sk(2))

  end subroutine pd


  !> rotation routine for interaction of a p orbital with an f orbital
  subroutine pf(hh, ll, mm, nn, sk)

    !> dimeric block to put the results in to
    real(dp), intent(inout) :: hh(:,:)

    !> directional cosine ll
    real(dp), intent(in) :: ll

    !> directional cosine mm
    real(dp), intent(in) :: mm

    !> directional cosine nn
    real(dp), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    real(dp), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 2)
    @:ASSERT(all(shape(hh) >= (/ 3, 7 /)))

    hh(1,1) = -(-1.0_dp+nn**2+ll**2)*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(2.0_dp)&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp+sqrt(15.0_dp)&
        &*(nn**4-nn**2+5.0_dp*nn**2*ll**2-3.0_dp*ll**2+4.0_dp*ll**4)&
        &*sk(2)/4.0_dp
    hh(2,1) = -(-1.0_dp+nn**2+ll**2)*ll*sqrt(15.0_dp)*nn*sk(1)&
        &+(3.0_dp*nn**2+3.0_dp*ll**2-2.0_dp)*ll*nn*sqrt(10.0_dp)*sk(2)&
        &/2.0_dp
    hh(3,1) = -(-1.0_dp+nn**2+ll**2)*sqrt(2.0_dp)*sqrt(3.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)/4.0_dp+(15.0_dp/4.0_dp&
        &*(nn**4)+ 15.0_dp/ 4.0_dp*(ll**2)*(nn**2)- 11.0_dp&
        &/4.0_dp*(nn**2)-(ll**2)/ 4.0_dp)*sk(2)
    hh(4,1) = mm*nn*(5.0_dp*nn**2-3.0_dp)*sk(1)/2.0_dp&
        &-(5.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*sqrt(2.0_dp)*nn*mm&
        &*sk(2)/4.0_dp
    hh(5,1) = mm*ll*sqrt(2.0_dp)*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)*sk(1)&
        &/4.0_dp-(15.0_dp*nn**2-1.0_dp)*ll*mm*sk(2)/4.0_dp
    hh(6,1) = mm*(2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(15.0_dp)*nn*&
        &sk(1)/2.0_dp-(3.0_dp*nn**2+6.0_dp*ll**2-1.0_dp)*nn*mm&
        &*sqrt(10.0_dp)*sk(2)/4.0_dp
    hh(7,1) = mm*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*sqrt(2.0_dp)&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp-ll*mm*sqrt(15.0_dp)&
        &*(3.0_dp*nn**2+4.0_dp*ll**2-1.0_dp)*sk(2)/4.0_dp
    hh(1,2) = sqrt(2.0_dp)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*nn&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp-mm*(4.0_dp*ll**2-1.0_dp+nn**2)&
        &*sqrt(15.0_dp)*nn*sk(2)/4.0_dp
    hh(2,2) = ll*mm*(nn**2)*sqrt(15.0_dp)*sk(1)-(3.0_dp*nn**2-1.0_dp)*ll*mm&
        &*sqrt(10.0_dp)*sk(2)/2.0_dp
    hh(3,2) = sqrt(2.0_dp)*mm*nn*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)*sk(1)&
        &/4.0_dp-(15.0_dp*nn**2-11.0_dp)*nn*mm*sk(2)/4.0_dp
    hh(4,2) = (nn**2*(5.0_dp*nn**2-3.0_dp)*sk(1))/2.0_dp&
        &-(5.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*sqrt(2.0_dp)*(-1.0_dp+nn**2)&
        &*sk(2)/4.0_dp
    hh(5,2) = sqrt(2.0_dp)*ll*nn*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)*sk(1)&
        &/4.0_dp-(15.0_dp*nn**2-11.0_dp)*nn*ll*sk(2)/4.0_dp
    hh(6,2) = (2.0_dp*ll**2-1.0_dp+nn**2)*(nn**2)*sqrt(15.0_dp)*sk(1)&
        &/2.0_dp-(3.0_dp*nn**2-1.0_dp)*(2.0_dp*ll**2-1.0_dp+nn**2)&
        &*sqrt(10.0_dp)*sk(2)/4.0_dp
    hh(7,2) = sqrt(2.0_dp)*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*nn&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp-ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)&
        &*sqrt(15.0_dp)*nn*sk(2)/4.0_dp
    hh(1,3) = ll*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(2.0_dp)*sqrt(5.0_dp)&
        &*sk(1)/4.0_dp-ll*mm*sqrt(15.0_dp)*(nn**2+4.0_dp*ll**2-3.0_dp)&
        &*sk(2)/4.0_dp
    hh(2,3) = (ll**2)*mm*sqrt(15.0_dp)*nn*sk(1)-(3.0_dp*ll**2-1.0_dp)*nn*mm&
        &*sqrt(10.0_dp)*sk(2)/2.0_dp
    hh(3,3) = ll*mm*sqrt(2.0_dp)*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)/4.0_dp-(15.0_dp*nn**2-1.0_dp)*ll*mm*sk(2)/4.0_dp
    hh(4,3) = (ll*nn*(5.0_dp*nn**2-3.0_dp)*sk(1))/2.0_dp&
        &-(5.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*sqrt(2.0_dp)*nn*ll*sk(2)&
        &/4.0_dp
    hh(5,3) = ll**2*sqrt(2.0_dp)*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)/4.0_dp+(-15.0_dp/4.0_dp*ll**2*(nn**2)+5.0_dp&
        &/4.0_dp*(nn**2)- 1.0_dp/4.0_dp+ll**2/4.0_dp)*sk(2)
    hh(6,3) = ll*(2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(15.0_dp)*nn*sk(1)&
        &/2.0_dp-(3.0_dp*nn**2+6.0_dp*ll**2-5.0_dp)*nn*ll*sqrt(10.0_dp)&
        &*sk(2)/4.0_dp
    hh(7,3) = (ll**2)*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*sqrt(2.0_dp)&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp-sqrt(15.0_dp)*(3.0_dp*ll**2*&
        &nn**2-nn**2-5.0_dp*ll**2+4.0_dp*ll**4+1.0_dp)*sk(2)/4.0_dp

  end subroutine pf


  !> rotation routine for interaction of a d orbital with a d orbital
  subroutine dd(hh, ll, mm, nn, sk)

    !> dimeric block to put the results in to
    real(dp), intent(inout) :: hh(:,:)

    !> directional cosine ll
    real(dp), intent(in) :: ll

    !> directional cosine mm
    real(dp), intent(in) :: mm

    !> directional cosine nn
    real(dp), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    real(dp), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 3)
    @:ASSERT(all(shape(hh) >= (/ 5, 5 /)))

    hh(1,1) = -3.0_dp*ll**2*(-1.0_dp+nn**2+ll**2)*sk(1)&
        &+(4.0_dp*ll**2*nn**2-nn**2+4.0_dp*ll**4-4.0_dp*ll**2+1.0_dp)*sk(2)&
        &+(-ll**2*nn**2+nn**2+ll**2-ll**4)*sk(3)
    hh(2,1) = -3.0_dp*ll*(-1.0_dp+nn**2+ll**2)*nn*sk(1)&
        &+(4.0_dp*nn**2+4.0_dp*ll**2-3.0_dp)*nn*ll*sk(2)&
        &-ll*(nn**2+ll**2)*nn*sk(3)
    hh(3,1) = ll*mm*sqrt(3.0_dp)*(3.0_dp*nn**2-1.0_dp)*sk(1)/2.0_dp&
        &-2.0_dp*sqrt(3.0_dp)*mm*ll*(nn**2)*sk(2)&
        &+ll*mm*(nn**2+1.0_dp)*sqrt(3.0_dp)*sk(3)/2.0_dp
    hh(4,1) = 3.0_dp*(ll**2)*mm*nn*sk(1)-(4.0_dp*ll**2-1.0_dp)*nn*mm&
        &*sk(2)+mm*(-1.0_dp+ll**2)*nn*sk(3)
    hh(5,1) = 3.0_dp/2.0_dp*mm*ll*(2.0_dp*ll**2-1.0_dp+nn**2)*sk(1)&
        &-2.0_dp*mm*ll*(2.0_dp*ll**2-1.0_dp+nn**2)*sk(2)&
        &+mm*ll*(2.0_dp*ll**2-1.0_dp+nn**2)*sk(3)/2.0_dp
    hh(1,2) = hh(2,1)
    hh(2,2) = -3.0_dp*(-1.0_dp+nn**2+ll**2)*nn**2*sk(1)&
        &+(4.0_dp*nn**4-4.0_dp*nn**2+4.0_dp*ll**2*nn**2+1.0_dp-ll**2)*sk(2)&
        &-(-1.0_dp+nn)*(nn**3+nn**2+ll**2*nn+ll**2)*sk(3)
    hh(3,2) = mm*sqrt(3.0_dp)*nn*(3.0_dp*nn**2-1.0_dp)*sk(1)/2.0_dp&
        &-nn*sqrt(3.0_dp)*(2.0_dp*nn**2-1.0_dp)*mm*sk(2)&
        &+(-1.0_dp+nn**2)*sqrt(3.0_dp)*nn*mm*sk(3)/2.0_dp
    hh(4,2) = 3.0_dp*mm*ll*(nn**2)*sk(1)-(4.0_dp*nn**2-1.0_dp)*mm*ll&
        &*sk(2)+ll*mm*(-1.0_dp+nn**2)*sk(3)
    hh(5,2) = 3.0_dp/ 2.0_dp*mm*(2.0_dp*ll**2-1.0_dp+nn**2)*nn*sk(1)&
        &-(2.0_dp*nn**2-1.0_dp+4.0_dp*ll**2)*nn*mm*sk(2)&
        &+mm*(nn**2+2.0_dp*ll**2+1.0_dp)*nn*sk(3)/2.0_dp
    hh(1,3) = hh(3,1)
    hh(2,3) = hh(3,2)
    hh(3,3) = ((3.0_dp*nn**2-1.0_dp)**2*sk(1))/4.0_dp&
        &-(3.0_dp*(-1.0_dp+nn**2)*nn**2*sk(2))&
        &+3.0_dp/4.0_dp*((-1.0_dp+nn**2)**2)*sk(3)
    hh(4,3) = ll*(3.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*nn*sk(1)/2.0_dp&
        &-(2.0_dp*nn**2-1.0_dp)*ll*nn*sqrt(3.0_dp)*sk(2)&
        &+nn*ll*sqrt(3.0_dp)*(-1.0_dp+nn**2)*sk(3)/2.0_dp
    hh(5,3) = (2.0_dp*ll**2-1.0_dp+nn**2)*(3.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)&
        &*sk(1)/4.0_dp-(2.0_dp*ll**2-1.0_dp+nn**2)*(nn**2)&
        &*sqrt(3.0_dp)*sk(2)+sqrt(3.0_dp)*(2.0_dp*ll**2-1.0_dp+nn**2)&
        &*(nn**2+1.0_dp)*sk(3)/4.0_dp
    hh(1,4) = hh(4,1)
    hh(2,4) = hh(4,2)
    hh(3,4) = hh(4,3)
    hh(4,4) = 3.0_dp*ll**2*nn**2*sk(1)+(-4.0_dp*ll**2*nn**2+nn**2+ll**2)&
        &*sk(2)+(-1.0_dp+nn)*(-nn+ll**2*nn-1.0_dp+ll**2)*sk(3)
    hh(5,4) = 3.0_dp/2.0_dp*ll*(2.0_dp*ll**2-1.0_dp+nn**2)*nn*sk(1)&
        &-((2.0_dp*nn**2-3.0_dp+4.0_dp*ll**2)*nn*ll*sk(2))&
        &+(ll*(nn**2-3.0_dp+2.0_dp*ll**2)*nn*sk(3))/2.0_dp
    hh(1,5) = hh(5,1)
    hh(2,5) = hh(5,2)
    hh(3,5) = hh(5,3)
    hh(4,5) = hh(5,4)
    hh(5,5) = 3.0_dp/4.0_dp*((2.0_dp*ll**2-1.0_dp+nn**2)**2)*sk(1)&
        &+((-nn**4+nn**2-4.0_dp*ll**2*nn**2-4.0_dp*ll**4&
        &+4.0_dp*ll**2)*sk(2))+((nn**4)/4.0_dp&
        &+(ll**2*nn**2)+(nn**2)/2.0_dp+1.0_dp/4.0_dp-(ll**2)+(ll**4))*sk(3)

  end subroutine dd


  !> rotation routine for interaction of a d orbital with an f orbital
  subroutine df(hh, ll, mm, nn, sk)

    !> dimeric block to put the results in to
    real(dp), intent(inout) :: hh(:,:)

    !> directional cosine ll
    real(dp), intent(in) :: ll

    !> directional cosine mm
    real(dp), intent(in) :: mm

    !> directional cosine nn
    real(dp), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    real(dp), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 3)
    @:ASSERT(all(shape(hh) >= (/ 5, 7 /)))

    hh(1,1) = -ll*(-1.0_dp+nn**2+ll**2)*(4.0_dp*ll**2-1.0_dp+nn**2)&
        &*sqrt(6.0_dp)*sqrt(5.0_dp)*sk(1)/4.0_dp&
        &+sqrt(15.0_dp)*ll*(2.0_dp*nn**4-5.0_dp*nn**2+10.0_dp*ll**2*nn**2&
        &+3.0_dp-10.0_dp*ll**2+8.0_dp*ll**4)*sk(2)/4.0_dp&
        &-sqrt(6.0_dp)*ll*(nn**4+5.0_dp*ll**2*nn**2-4.0_dp*nn**2+1.0_dp&
        &+4.0_dp*ll**4-5.0_dp*ll**2)*sk(3)/4.0_dp
    hh(2,1) = -3.0_dp*(ll**2)*(-1.0_dp+nn**2+ll**2)*sqrt(5.0_dp)*nn*sk(1)&
        &+(6.0_dp*ll**2*nn**2-nn**2+1.0_dp+6.0_dp*ll**4-6.0_dp*ll**2)&
        &*sqrt(10.0_dp)*nn*sk(2)/2.0_dp-(nn*(3.0_dp*ll**2*nn**2-2.0_dp*nn**2&
        &+1.0_dp-3.0_dp*ll**2+3.0_dp*ll**4)*sk(3))
    hh(3,1) = -3.0_dp/4.0_dp*ll*(-1.0_dp+nn**2+ll**2)*sqrt(2.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)+((30.0_dp*nn**4+30.0_dp*ll**2*nn**2&
        &-27.0_dp*nn**2-2.0_dp*ll**2+1.0_dp)*ll*sk(2))/4.0_dp&
        &-ll*sqrt(10.0_dp)*(3.0_dp*nn**4+3.0_dp*ll**2*nn**2+ll**2-1.0_dp)&
        &*sk(3)/4.0_dp
    hh(4,1) = ll*mm*sqrt(3.0_dp)*nn*(5.0_dp*nn**2-3.0_dp)*sk(1)/2.0_dp&
        &-(5.0_dp*nn**2-1.0_dp)*nn*ll*mm*sqrt(6.0_dp)*sk(2)/2.0_dp&
        &+ll*mm*(nn**2+1.0_dp)*sqrt(15.0_dp)*nn*sk(3)/2.0_dp
    hh(5,1) = 3.0_dp/4.0_dp*(ll**2)*mm*sqrt(2.0_dp)*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)-(30.0_dp*ll**2*nn**2-5.0_dp*nn**2-2.0_dp*ll**2+1.0_dp)*mm&
        &*sk(2)/4.0_dp+mm*sqrt(10.0_dp)*(3.0_dp*ll**2*nn**2-2.0_dp*nn**2&
        &+ll**2)*sk(3)/4.0_dp
    hh(6,1) = 3.0_dp/2.0_dp*ll*mm*(2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(5.0_dp)*nn&
        &*sk(1)-3.0_dp/2.0_dp*nn*sqrt(10.0_dp)*mm*ll&
        &*(2.0_dp*ll**2-1.0_dp+nn**2)*sk(2)+3.0_dp/2.0_dp*ll*mm&
        &*(2.0_dp*ll**2-1.0_dp+nn**2)*nn*sk(3)
    hh(7,1) = (ll**2)*mm*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*sqrt(6.0_dp)&
        &*sqrt(5.0_dp)*sk(1)/ 4.0_dp-sqrt(15.0_dp)*mm*&
        &(6.0_dp*ll**2*nn**2-nn**2+1.0_dp+8.0_dp*ll**4-6.0_dp*ll**2)&
        &*sk(2)/4.0_dp+sqrt(6.0_dp)*mm*(3.0_dp*ll**2*nn**2-2.0_dp*nn**2&
        &+4.0_dp*ll**4-3.0_dp*ll**2)*sk(3)/4.0_dp
    hh(1,2) = -(-1.0_dp+nn**2+ll**2)*(4.0_dp*ll**2-1.0_dp+nn**2)&
        &*sqrt(6.0_dp)*nn*sqrt(5.0_dp)*sk(1)/4.0_dp&
        &+sqrt(15.0_dp)*nn*(2.0_dp*nn**4-3.0_dp*nn**2+10.0_dp*ll**2*nn**2&
        &+1.0_dp+8.0_dp*ll**4-8.0_dp*ll**2)*sk(2)/4.0_dp&
        &-sqrt(6.0_dp)*nn*(nn**4+5.0_dp*ll**2*nn**2-1.0_dp+4.0_dp*ll**4-ll**2)&
        &*sk(3)/4.0_dp
    hh(2,2) = -3.0_dp*(-1.0_dp+nn**2+ll**2)*ll*(nn**2)*sqrt(5.0_dp)*sk(1)&
        &+(6.0_dp*nn**4-6.0_dp*nn**2+6.0_dp*ll**2*nn**2+1.0_dp-ll**2)*ll&
        &*sqrt(10.0_dp)*sk(2)/2.0_dp-(ll*(3.0_dp*nn**4+3.0_dp*ll**2*nn**2&
        &-3.0_dp*nn**2+1.0_dp-2.0_dp*ll**2)*sk(3))
    hh(3,2) = -3.0_dp/4.0_dp*(-1.0_dp+nn**2+ll**2)*sqrt(2.0_dp)*nn&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)+((30.0_dp*nn**4-37.0_dp*nn**2&
        &+30.0_dp*ll**2*nn**2+11.0_dp-12.0_dp*ll**2)*nn*sk(2))&
        &/4.0_dp-(-1.0_dp+nn)*sqrt(10.0_dp)*(3.0_dp*nn**3+3.0_dp*nn**2-nn&
        &+3.0_dp*ll**2*nn-1.0_dp+3.0_dp*ll**2)*nn*sk(3)/4.0_dp
    hh(4,2) = mm*sqrt(3.0_dp)*(nn**2)*(5.0_dp*nn**2-3.0_dp)*sk(1)&
        &/2.0_dp-(5.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*sqrt(2.0_dp)*&
        &(2.0_dp*nn**2-1.0_dp)*mm*sk(2)/4.0_dp+(-1.0_dp+nn**2)&
        &*sqrt(15.0_dp)*(nn**2)*mm*sk(3)/2.0_dp
    hh(5,2) = 3.0_dp/4.0_dp*mm*ll*sqrt(2.0_dp)*nn*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)-3.0_dp/2.0_dp*(5.0_dp*nn**2-2.0_dp)*mm*ll*nn&
        &*sk(2)+3.0_dp/4.0_dp*nn*sqrt(10.0_dp)*ll*mm*(-1.0_dp+nn**2)&
        &*sk(3)
    hh(6,2) = 3.0_dp/2.0_dp*mm*(2.0_dp*ll**2-1.0_dp+nn**2)*(nn**2)&
        &*sqrt(5.0_dp)*sk(1)-(6.0_dp*nn**4+12.0_dp*ll**2*nn**2-5.0_dp*nn**2&
        &+1.0_dp-2.0_dp*ll**2)*mm*sqrt(10.0_dp)*sk(2)/4.0_dp&
        &+mm*(3.0_dp*nn**4-nn**2+6.0_dp*ll**2*nn**2-4.0_dp*ll**2)*sk(3)&
        &/2.0_dp
    hh(7,2) = mm*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*sqrt(6.0_dp)*nn&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp-mm*ll*sqrt(15.0_dp)*nn&
        &*(3.0_dp*nn**2-2.0_dp+4.0_dp*ll**2)*sk(2)/2.0_dp&
        &+sqrt(6.0_dp)*mm*ll*nn*(3.0_dp*nn**2+1.0_dp+4.0_dp*ll**2)&
        &*sk(3)/4.0_dp
    hh(1,3) = sqrt(2.0_dp)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)&
        &*(3.0_dp*nn**2-1.0_dp)*sqrt(5.0_dp)*sk(1)/8.0_dp&
        &-3.0_dp/4.0_dp*(nn**2)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(5.0_dp)&
        &*sk(2)+3.0_dp/8.0_dp*sqrt(2.0_dp)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)&
        &*(nn**2+1.0_dp)*sk(3)
    hh(2,3) = ll*mm*(3.0_dp*nn**2-1.0_dp)*sqrt(15.0_dp)*nn*sk(1)/2.0_dp&
        &-(3.0_dp*nn**2-1.0_dp)*ll*nn*mm*sqrt(30.0_dp)*sk(2)/2.0_dp&
        &+sqrt(3.0_dp)*ll*mm*nn*(3.0_dp*nn**2-1.0_dp)*sk(3)/2.0_dp
    hh(3,3) = sqrt(2.0_dp)*mm*(3.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)/8.0_dp-(15.0_dp*nn**2-11.0_dp)&
        &*mm*(nn**2)*sqrt(3.0_dp)*sk(2)/4.0_dp&
        &+(3.0_dp*nn**3+3.0_dp*nn**2-nn-1.0_dp)*(-1.0_dp+nn)*sqrt(5.0_dp)&
        &*sqrt(2.0_dp)*mm*sqrt(3.0_dp)*sk(3)/8.0_dp
    hh(4,3) = ((3.0_dp*nn**2-1.0_dp)*nn*(5.0_dp*nn**2-3.0_dp)*sk(1))&
        &/4.0_dp-3.0_dp/4.0_dp*(5.0_dp*nn**2-1.0_dp)*(-1.0_dp+nn**2)*nn&
        &*sqrt(2.0_dp)*sk(2)+3.0_dp/4.0_dp*((-1.0_dp+nn**2)**2)&
        &*sqrt(5.0_dp)*nn*sk(3)
    hh(5,3) = sqrt(2.0_dp)*ll*(3.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*&
        &(5.0_dp*nn**2-1.0_dp)*sk(1)/8.0_dp-(15.0_dp*nn**2-11.0_dp)&
        &*ll*(nn**2)*sqrt(3.0_dp)*sk(2)/4.0_dp+(3.0_dp*nn**3+3.0_dp*nn**2&
        &-nn-1.0_dp)*(-1.0_dp+nn)*sqrt(5.0_dp)*sqrt(2.0_dp)*ll*sqrt(3.0_dp)&
        &*sk(3)/8.0_dp
    hh(6,3) = (2.0_dp*ll**2-1.0_dp+nn**2)*(3.0_dp*nn**2-1.0_dp)*sqrt(15.0_dp)&
        &*nn*sk(1)/4.0_dp-(3.0_dp*nn**2-1.0_dp)*(2.0_dp*ll**2-1.0_dp+nn**2)&
        &*nn*sqrt(30.0_dp)*sk(2)/4.0_dp+sqrt(3.0_dp)&
        &*(2.0_dp*ll**2-1.0_dp+nn**2)*nn*(3.0_dp*nn**2-1.0_dp)*sk(3)/&
        &4.0_dp
    hh(7,3) = sqrt(2.0_dp)*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)&
        &*(3.0_dp*nn**2-1.0_dp)*sqrt(5.0_dp)*sk(1)/8.0_dp&
        &-3.0_dp/4.0_dp*(nn**2)*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)&
        &*sqrt(5.0_dp)*sk(2)+ 3.0_dp/8.0_dp*sqrt(2.0_dp)*ll&
        &*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*(nn**2+1.0_dp)*sk(3)
    hh(1,4) = ll*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(6.0_dp)*nn&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp-ll*mm*sqrt(15.0_dp)*nn&
        &*(nn**2-2.0_dp+4.0_dp*ll**2)*sk(2)/2.0_dp&
        &+sqrt(6.0_dp)*mm*ll*nn*(nn**2+4.0_dp*ll**2-5.0_dp)*sk(3)/4.0_dp
    hh(2,4) = 3.0_dp*(ll**2)*mm*(nn**2)*sqrt(5.0_dp)*sk(1)&
        &-(6.0_dp*ll**2.0_dp*nn**2-nn**2-ll**2)*mm*sqrt(10.0_dp)*sk(2)&
        &/2.0_dp+mm*(-2.0_dp*nn**2+3.0_dp*ll**2.0_dp*nn**2+1.0_dp-2.0_dp*ll**2)&
        &*sk(3)
    hh(3,4) = 3.0_dp/4.0_dp*ll*mm*sqrt(2.0_dp)*nn*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)-3.0_dp/2.0_dp*(5.0_dp*nn**2-2.0_dp)*ll*mm*nn*&
        &sk(2)+3.0_dp/4.0_dp*nn*sqrt(10.0_dp)*mm*ll*(-1.0_dp+nn**2)&
        &*sk(3)
    hh(4,4) = ll*sqrt(3.0_dp)*(nn**2)*(5.0_dp*nn**2-3.0_dp)*sk(1)/2.0_dp&
        &-(5.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*sqrt(2.0_dp)*(2.0_dp*nn**2-1.0_dp)&
        &*ll*sk(2)/4.0_dp+(-1.0_dp+nn**2)*sqrt(15.0_dp)*(nn**2)*ll&
        &*sk(3)/2.0_dp
    hh(5,4) = 3.0_dp/ 4.0_dp*ll**2.0_dp*sqrt(2.0_dp)*nn*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)-(30.0_dp*ll**2.0_dp*(nn**2)-(5.0_dp*nn**2)&
        &-12.0_dp*ll**2+1.0_dp)*nn*sk(2)/4.0_dp+(-1.0_dp+nn)&
        &*sqrt(10.0_dp)*(-(2.0_dp*nn)+3.0_dp*ll**2.0_dp*nn-2.0_dp+3.0_dp*ll**2)&
        &*nn*sk(3)/4.0_dp
    hh(6,4) = 3.0_dp/2.0_dp*ll*(2.0_dp*ll**2-1.0_dp+nn**2)*(nn**2)&
        &*sqrt(5.0_dp)&
        &*sk(1)-(6.0_dp*nn**4+12.0_dp*ll**2.0_dp*nn**2-9.0_dp*nn**2&
        &+1.0_dp-2.0_dp*ll**2)*ll*sqrt(10.0_dp)*sk(2)/4.0_dp&
        &+(ll*(3.0_dp*nn**4-9.0_dp*nn**2+6.0_dp*ll**2.0_dp*nn**2+4.0_dp&
        &-4.0_dp*ll**2)*sk(3))/2.0_dp
    hh(7,4) = (ll**2)*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*sqrt(6.0_dp)*nn&
        &*sqrt(5.0_dp)*sk(1)/4.0_dp-sqrt(15.0_dp)*nn&
        &*(6.0_dp*ll**2.0_dp*nn**2-nn**2+1.0_dp+8.0_dp*ll**4-8.0_dp*ll**2)&
        &*sk(2)/4.0_dp+sqrt(6.0_dp)*nn*(3.0_dp*ll**2.0_dp*nn**2&
        &-2.0_dp*nn**2-7.0_dp*ll**2+2.0_dp+4.0_dp*ll**4)*sk(3)/4.0_dp
    hh(1,5) = (2.0_dp*ll**2-1.0_dp+nn**2)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)&
        &*sqrt(6.0_dp)*sqrt(5.0_dp)*sk(1)/8.0_dp-sqrt(15.0_dp)*mm&
        &*(nn**4-nn**2+6.0_dp*ll**2.0_dp*nn**2+8.0_dp*ll**4-6.0_dp*ll**2)&
        &*sk(2)/4.0_dp+sqrt(6.0_dp)*mm*(nn**4+6.0_dp*ll**2.0_dp*nn**2&
        &+2.0_dp*nn**2+1.0_dp+8.0_dp*ll**4-6.0_dp*ll**2)*sk(3)/8.0_dp
    hh(2,5) = 3.0_dp/2.0_dp*(2.0_dp*ll**2-1.0_dp+nn**2)*ll*mm*sqrt(5.0_dp)*nn&
        &*sk(1)-3.0_dp/2.0_dp*nn*sqrt(10.0_dp)*(2.0_dp*ll**2-1.0_dp&
        &+nn**2)*ll*mm*sk(2)+ 3.0_dp/2.0_dp*(2.0_dp*ll**2-1.0_dp+nn**2)*ll*mm&
        &*nn*sk(3)
    hh(3,5) = 3.0_dp/8.0_dp*(2.0_dp*ll**2-1.0_dp+nn**2)*mm*sqrt(2.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)-(15.0_dp*nn**4+30.0_dp*ll**2*nn**2&
        &-11.0_dp*nn**2-2.0_dp*ll**2)*mm*sk(2)/4.0_dp&
        &+mm*sqrt(10.0_dp)*(3.0_dp*nn**4+2.0_dp*nn**2+6.0_dp*ll**2*nn**2-1.0_dp&
        &+2.0_dp*ll**2)*sk(3)/8.0_dp
    hh(4,5) = (2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(3.0_dp)*nn*(5.0_dp*nn**2-3.0_dp)&
        &*sk(1)/4.0_dp-(5.0_dp*nn**2-1.0_dp)*nn*(2.0_dp*ll**2-1.0_dp&
        &+nn**2)*sqrt(6.0_dp)*sk(2)/4.0_dp+(2.0_dp*ll**2-1.0_dp+nn**2)&
        &*(nn**2+1.0_dp)*sqrt(15.0_dp)*nn*sk(3)/4.0_dp
    hh(5,5) = 3.0_dp/8.0_dp*(2.0_dp*ll**2-1.0_dp+nn**2)*ll*sqrt(2.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)-((15.0_dp*nn**4+30.0_dp*ll**2*nn**2&
        &-21.0_dp*nn**2+2.0_dp-2.0_dp*ll**2)*ll*sk(2))/4.0_dp&
        &+ll*sqrt(10.0_dp)*(3.0_dp*nn**4+6.0_dp*ll**2*nn**2-6.0_dp*nn**2&
        &+2.0_dp*ll**2-1.0_dp)*sk(3)/8.0_dp
    hh(6,5) = 3.0_dp/4.0_dp*((2.0_dp*ll**2-1.0_dp+nn**2)**2)*sqrt(5.0_dp)*nn*&
        &sk(1)-(3.0_dp*nn**4+12.0_dp*ll**2*nn**2-4.0_dp*nn**2+12.0_dp*ll**4&
        &+1.0_dp-12.0_dp*ll**2)*sqrt(10.0_dp)*nn*sk(2)/4.0_dp&
        &+(nn*(3.0_dp*nn**4+12.0_dp*ll**2*nn**2+2.0_dp*nn**2-12.0_dp*ll**2&
        &-1.0_dp+12.0_dp*ll**4)*sk(3))/4.0_dp
    hh(7,5) = (2.0_dp*ll**2-1.0_dp+nn**2)*ll*(4.0_dp*ll**2-3.0_dp&
        &+3.0_dp*nn**2)*sqrt(6.0_dp)*sqrt(5.0_dp)*sk(1)/8.0_dp&
        &-sqrt(15.0_dp)*ll*(3.0_dp*nn**4+10.0_dp*ll**2*nn**2-5.0_dp*nn**2&
        &-10.0_dp*ll**2+8.0_dp*ll**4+2.0_dp)*sk(2)/4.0_dp+sqrt(6.0_dp)*ll&
        &*(3.0_dp*nn**4+10.0_dp*ll**2*nn**2-2.0_dp*nn**2+3.0_dp+8.0_dp*ll**4&
        &-10.0_dp*ll**2)*sk(3)/8.0_dp

  end subroutine df


  !> rotation routine for interaction of an f orbital with an f orbital
  subroutine ff(hh, ll, mm, nn, sk)

    !> dimeric block to put the results in to
    real(dp), intent(inout) :: hh(:,:)

    !> directional cosine ll
    real(dp), intent(in) :: ll

    !> directional cosine mm
    real(dp), intent(in) :: mm

    !> directional cosine nn
    real(dp), intent(in) :: nn

    !> Slater-Koster table for dimer element of the Slater-Koster table
    real(dp), intent(in) :: sk(:)

    @:ASSERT(size(sk) == 4)
    @:ASSERT(all(shape(hh) >= (/ 7, 7 /)))

    hh(1,1) = - 5.0_dp/ 8.0_dp*(-1.0_dp+nn**2+ll**2)*((4.0_dp*ll**2-1.0_dp&
        &+nn**2)**2)*sk(1)+(15.0_dp/16.0_dp*(nn**6)- 15.0_dp&
        &/8.0_dp*(nn**4)+135.0_dp/ 16.0_dp*(nn**4)*(ll**2)-135.0_dp/8.0_dp&
        &*(ll**2)*(nn**2)+15.0_dp/16.0_dp*(nn**2)+45.0_dp/2.0_dp*(nn**2)&
        &*(ll**4)+135.0_dp/16.0_dp*(ll**2)-45.0_dp/2.0_dp*(ll**4)&
        &+(15.0_dp*ll**6))*sk(2)+(- 3.0_dp/ 8.0_dp*(nn**6)&
        &- 27.0_dp/8.0_dp*(nn**4)&
        &*(ll**2)-3.0_dp/8.0_dp*(nn**4)+3.0_dp/8.0_dp*(nn**2)-(9.0_dp*nn**2&
        &*ll**4)+27.0_dp/4.0_dp*(ll**2)*(nn**2)+3.0_dp/8.0_dp-(6.0_dp*ll**6)&
        &-27.0_dp/8.0_dp*(ll**2)+(9.0_dp*ll**4))*sk(3)&
        &+((nn**6)/16.0_dp+3.0_dp/8.0_dp*(nn**4)+9.0_dp/16.0_dp*(nn**4)*(ll**2)&
        &-9.0_dp/8.0_dp*(ll**2)*(nn**2)+9.0_dp/16.0_dp*(nn**2)&
        &+3.0_dp/2.0_dp*(nn**2)*(ll**4)+9.0_dp/16.0_dp*(ll**2)-3.0_dp/2.0_dp&
        &*(ll**4)+(ll**6))*sk(4)
    hh(2,1) = - 5.0_dp/4.0_dp*(-1.0_dp+nn**2+ll**2)*(4.0_dp*ll**2-1.0_dp&
        &+nn**2)*ll*sqrt(6.0_dp)*nn*sk(1)+5.0_dp/8.0_dp*sqrt(6.0_dp)*nn*ll&
        &*(3.0_dp*nn**4+15.0_dp*ll**2*nn**2-7.0_dp*nn**2+4.0_dp-15.0_dp*ll**2&
        &+12.0_dp*ll**4)*sk(2)-sqrt(6.0_dp)*nn*ll*(3.0_dp*nn**4&
        &+15.0_dp*ll**2*nn**2-10.0_dp*nn**2+5.0_dp-15.0_dp*ll**2+12.0_dp*ll**4)&
        &*sk(3)/4.0_dp+ll*sqrt(6.0_dp)*(nn**4+5.0_dp*ll**2*nn**2&
        &-5.0_dp*nn**2+4.0_dp*ll**4-5.0_dp*ll**2)*nn*sk(4)/8.0_dp
    hh(3,1) = -(-1.0_dp+nn**2+ll**2)*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(5.0_dp)&
        &*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)*sk(1)/8.0_dp&
        &+sqrt(15.0_dp)*(15.0_dp*nn**6-26.0_dp*nn**4+75.0_dp*ll**2*nn**4&
        &-70.0_dp*ll**2*nn**2+11.0_dp*nn**2+60.0_dp*ll**4*nn**2-4.0_dp*ll**4&
        &+3.0_dp*ll**2)*sk(2)/16.0_dp-sqrt(15.0_dp)*(3.0_dp*nn**6-nn**4&
        &+15.0_dp*ll**2*nn**4-3.0_dp*nn**2-2.0_dp*ll**2*nn**2+12.0_dp&
        &*ll**4*nn**2+1.0_dp+4.0_dp*ll**4-5.0_dp*ll**2)*sk(3)/8.0_dp&
        &+sqrt(15.0_dp)*(nn**6+2.0_dp*nn**4+5.0_dp*ll**2*nn**4-3.0_dp*nn**2&
        &+6.0_dp*ll**2*nn**2+4.0_dp*ll**4*nn**2-3.0_dp*ll**2+4.0_dp*ll**4)&
        &*sk(4)/16.0_dp
    hh(4,1) = mm*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(2.0_dp)&
        &*sqrt(5.0_dp)*nn*(5.0_dp*nn**2-3.0_dp)*sk(1)/8.0_dp-3.0_dp&
        &/16.0_dp*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*sqrt(5.0_dp)*nn*sqrt(2.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(2)+3.0_dp/8.0_dp*sqrt(10.0_dp)*mm&
        &*(4.0_dp*ll**2-1.0_dp+nn**2)*(nn**2+1.0_dp)*nn*sk(3)&
        &-sqrt(5.0_dp)*sqrt(2.0_dp)*(nn**2+3.0_dp)*nn*(4.0_dp*ll**2-1.0_dp&
        &+nn**2)*mm*sk(4)/16.0_dp
    hh(5,1) = mm*(4.0_dp*ll**2-1.0_dp+nn**2)*ll*sqrt(5.0_dp)*sqrt(3.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)/8.0_dp-mm*sqrt(15.0_dp)*ll&
        &*(15.0_dp*nn**4+60.0_dp*ll**2*nn**2-26.0_dp*nn**2+3.0_dp-4.0_dp*ll**2)&
        &*sk(2)/ 16.0_dp+mm*sqrt(15.0_dp)*ll*(3.0_dp*nn**4&
        &+12.0_dp*ll**2*nn**2-10.0_dp*nn**2+4.0_dp*ll**2-1.0_dp)*sk(3)&
        &/8.0_dp-mm*sqrt(15.0_dp)*ll*(nn**4+4.0_dp*ll**2*nn**2-6.0_dp*nn**2&
        &+4.0_dp*ll**2-3.0_dp)*sk(4)/16.0_dp
    hh(6,1) = 5.0_dp/8.0_dp*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*(2.0_dp*ll**2&
        &-1.0_dp+nn**2)*sqrt(6.0_dp)*nn*sk(1)-5.0_dp/16.0_dp*sqrt(6.0_dp)*mm*nn&
        &*(3.0_dp*nn**4+18.0_dp*ll**2*nn**2-4.0_dp*nn**2+24.0_dp*ll**4+1.0_dp&
        &-18.0_dp*ll**2)*sk(2)+sqrt(6.0_dp)*mm*nn*(3.0_dp*nn**4&
        &+18.0_dp*ll**2*nn**2+2.0_dp*nn**2+24.0_dp*ll**4-18.0_dp*ll**2-1.0_dp)&
        &*sk(3)/8.0_dp-mm*sqrt(6.0_dp)*(nn**4+6.0_dp*ll**2*nn**2&
        &+4.0_dp*nn**2+3.0_dp-6.0_dp*ll**2+8.0_dp*ll**4)*nn*sk(4)/16.0_dp
    hh(7,1) = 5.0_dp/8.0_dp*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*mm&
        &*(4.0_dp*ll**2-1.0_dp+nn**2)*sk(1)-15.0_dp/16.0_dp*ll&
        &*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)&
        &*sk(2)+3.0_dp/8.0_dp*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*mm&
        &*(4.0_dp*ll**2-1.0_dp+nn**2)*sk(3)-ll*(4.0_dp*ll**2-3.0_dp&
        &+3.0_dp*nn**2)*mm*(4.0_dp*ll**2-1.0_dp+nn**2)*sk(4)/16.0_dp
    hh(1,2) = hh(2,1)
    hh(2,2) = -(15.0_dp*ll**2*(-1.0_dp+nn**2+ll**2)*nn**2*sk(1))&
        &+(45.0_dp/2.0_dp*(ll**2)*(nn**4)-5.0_dp/2.0_dp*(nn**4)&
        &-(25.0_dp*ll**2*nn**2)+45.0_dp/2.0_dp*(ll**4)*(nn**2)+5.0_dp/2.0_dp&
        &*(nn**2)+5.0_dp/2.0_dp*(ll**2)-5.0_dp/2.0_dp*(ll**4))*sk(2)&
        &+((-9.0_dp*ll**2*nn**4+4.0_dp*nn**4+13.0_dp*ll**2*nn**2-9.0_dp*ll**4&
        &*nn**2-4.0_dp*nn**2-4.0_dp*ll**2+4.0_dp*ll**4+1.0_dp)*sk(3))&
        &+(3.0_dp/2.0_dp*(ll**2)*(nn**4)-3.0_dp/2.0_dp*(nn**4)+3.0_dp/2.0_dp&
        &*(nn**2)+3.0_dp/2.0_dp*(ll**4)*(nn**2)-(3.0_dp*ll**2*nn**2)&
        &+3.0_dp/2.0_dp*(ll**2)-3.0_dp/2.0_dp*(ll**4))*sk(4)
    hh(3,2) = - 3.0_dp/4.0_dp*ll*(-1.0_dp+nn**2+ll**2)*sqrt(10.0_dp)*nn&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)+(45.0_dp*nn**4-53.0_dp*nn**2&
        &+45.0_dp*ll**2*nn**2+12.0_dp-13.0_dp*ll**2)*nn*ll*sqrt(10.0_dp)&
        &*sk(2)/8.0_dp-sqrt(10.0_dp)*ll*(9.0_dp*nn**4+9.0_dp*ll**2*nn**2&
        &-10.0_dp*nn**2+3.0_dp-5.0_dp*ll**2)*nn*sk(3)/4.0_dp&
        &+3.0_dp/8.0_dp*ll*sqrt(2.0_dp)*(-1.0_dp+nn)*sqrt(5.0_dp)&
        &*(nn**3+nn**2+ll**2*nn+ll**2)*nn*sk(4)
    hh(4,2) = ll*mm*sqrt(15.0_dp)*(nn**2)*(5.0_dp*nn**2-3.0_dp)*sk(1)&
        &/2.0_dp-(5.0_dp*nn**2-1.0_dp)*(3.0_dp*nn**2-1.0_dp)*ll*mm&
        &*sqrt(15.0_dp)*sk(2)/4.0_dp+ll*mm*(nn**2)*(3.0_dp*nn**2-1.0_dp)&
        &*sqrt(15.0_dp)*sk(3)/2.0_dp-ll*mm*sqrt(3.0_dp)*(nn**3+nn**2+nn+1.0_dp)&
        &*sqrt(5.0_dp)*(-1.0_dp+nn)*sk(4)/4.0_dp
    hh(5,2) = 3.0_dp/4.0_dp*(ll**2)*mm*sqrt(10.0_dp)*nn*(5.0_dp*nn**2-1.0_dp)&
        &*sk(1)-(45.0_dp*ll**2*nn**2-5.0_dp*nn**2-13.0_dp*ll**2+1.0_dp)*nn&
        &*mm*sqrt(10.0_dp)*sk(2)/8.0_dp+sqrt(10.0_dp)*mm*(-4.0_dp*nn**2&
        &+9.0_dp*ll**2*nn**2+2.0_dp-5.0_dp*ll**2)*nn*sk(3)/4.0_dp&
        &-3.0_dp/8.0_dp*mm*sqrt(2.0_dp)*(-1.0_dp+nn)*sqrt(5.0_dp)&
        &*(-nn+ll**2*nn-1.0_dp+ll**2)*nn*sk(4)
    hh(6,2) = 15.0_dp/2.0_dp*ll*mm*(2.0_dp*ll**2-1.0_dp+nn**2)*(nn**2)&
        &*sk(1)-5.0_dp/4.0_dp*(9.0_dp*nn**2-1.0_dp)*(2.0_dp*ll**2&
        &-1.0_dp+nn**2)*ll*mm*sk(2)+(2.0_dp*ll**2-1.0_dp+nn**2)*ll*mm&
        &*(9.0_dp*nn**2-4.0_dp)*sk(3)/2.0_dp-3.0_dp/4.0_dp&
        &*(-1.0_dp+nn**2)*(2.0_dp*ll**2-1.0_dp+nn**2)*ll*mm*sk(4)
    hh(7,2) = 5.0_dp/4.0_dp*(ll**2)*mm*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)&
        &*sqrt(6.0_dp)*nn*sk(1)-5.0_dp/8.0_dp*sqrt(6.0_dp)*mm*nn&
        &*(9.0_dp*nn**2*ll**2-nn**2+1.0_dp+12.0_dp*ll**4-9.0_dp*ll**2)&
        &*sk(2)+sqrt(6.0_dp)*mm*nn*(9.0_dp*nn**2*ll**2-4.0_dp*nn**2&
        &+12.0_dp*ll**4-9.0_dp*ll**2+2.0_dp)*sk(3)/4.0_dp&
        &-mm*sqrt(6.0_dp)*(3.0_dp*nn**2*ll**2-3.0_dp*nn**2-1.0_dp-3.0_dp*ll**2&
        &+4.0_dp*ll**4)*nn*sk(4)/8.0_dp
    hh(1,3) = hh(3,1)
    hh(2,3) = hh(3,2)
    hh(3,3) = -3.0_dp/8.0_dp*(-1.0_dp+nn**2+ll**2)*((5.0_dp*nn**2-1.0_dp)**2)&
        &*sk(1)+(225.0_dp/16.0_dp*(nn**6)-165.0_dp/8.0_dp*(nn**4)&
        &+225.0_dp/16.0_dp*(ll**2)*(nn**4)+121.0_dp/16.0_dp*(nn**2)-65.0_dp&
        &/8.0_dp*(nn**2)*(ll**2)+(ll**2)/16.0_dp)*sk(2)-5.0_dp/8.0_dp&
        &*(-1.0_dp+nn)*(9.0_dp*nn**5+9.0_dp*nn**4-6.0_dp*nn**3+9.0_dp*ll**2&
        &*nn**3-6.0_dp*nn**2+9.0_dp*nn**2*ll**2+nn-ll**2*nn+1.0_dp-ll**2)*sk(3)&
        &+15.0_dp/16.0_dp*(-nn**2+nn**4+nn**2*ll**2-ll**2)*(-1.0_dp+nn**2)&
        &*sk(4)
    hh(4,3) = mm*sqrt(2.0_dp)*sqrt(3.0_dp)*(5.0_dp*nn**2-1.0_dp)*nn&
        &*(5.0_dp*nn**2-3.0_dp)*sk(1)/8.0_dp-(5.0_dp*nn**2-1.0_dp)&
        &*sqrt(3.0_dp)*sqrt(2.0_dp)*(15.0_dp*nn**2-11.0_dp)*nn*mm&
        &*sk(2)/ 16.0_dp+5.0_dp/8.0_dp*(-1.0_dp+nn**2)*nn*sqrt(3.0_dp)&
        &*(3.0_dp*nn**2-1.0_dp)*sqrt(2.0_dp)*mm*sk(3)-5.0_dp/16.0_dp&
        &*(-2.0_dp*nn**2+1.0_dp+nn**4)*sqrt(2.0_dp)*nn*sqrt(3.0_dp)*mm&
        &*sk(4)
    hh(5,3) = 3.0_dp/8.0_dp*mm*ll*((5.0_dp*nn**2-1.0_dp)**2)*sk(1)&
        &-(225.0_dp*nn**4-130.0_dp*nn**2+1.0_dp)*mm*ll*sk(2)/16.0_dp&
        &+5.0_dp/8.0_dp*(9.0_dp*nn**3+9.0_dp*nn**2-nn-1.0_dp)*ll*(-1.0_dp+nn)&
        &*mm*sk(3)- 15.0_dp/16.0_dp*((-1.0_dp+nn**2)**2)*mm*ll*sk(4)
    hh(6,3) = 3.0_dp/8.0_dp*mm*(2.0_dp*ll**2-1.0_dp+nn**2)*(5.0_dp*nn**2&
        &-1.0_dp)*sqrt(10.0_dp)*nn*sk(1)-(45.0_dp*nn**4+90.0_dp*nn**2*ll**2&
        &-48.0_dp*nn**2+11.0_dp-26.0_dp*ll**2)*nn*mm*sqrt(10.0_dp)&
        &*sk(2)/16.0_dp+sqrt(10.0_dp)*mm*(9.0_dp*nn**4-6.0_dp*nn**2&
        &+18.0_dp*nn**2*ll**2-10.0_dp*ll**2+1.0_dp)*nn*sk(3)/8.0_dp&
        &-3.0_dp/16.0_dp*mm*(-1.0_dp+nn)*sqrt(5.0_dp)*sqrt(2.0_dp)*(nn**3+nn**2&
        &+2.0_dp*nn*ll**2+nn+2.0_dp*ll**2+1.0_dp)*nn*sk(4)
    hh(7,3) = mm*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*sqrt(3.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sqrt(5.0_dp)*sk(1)/8.0_dp-mm*ll*sqrt(15.0_dp)&
        &*(45.0_dp*nn**4+60.0_dp*nn**2*ll**2-38.0_dp*nn**2+1.0_dp-4.0_dp*ll**2)&
        &*sk(2)/16.0_dp+mm*ll*sqrt(15.0_dp)*(9.0_dp*nn**4+2.0_dp*nn**2&
        &+12.0_dp*nn**2*ll**2+4.0_dp*ll**2-3.0_dp)*sk(3)/8.0_dp-mm&
        &*sqrt(15.0_dp)*ll*(3.0_dp*nn**4+4.0_dp*nn**2*ll**2+6.0_dp*nn**2-1.0_dp&
        &+4.0_dp*ll**2)*sk(4)/ 16.0_dp
    hh(1,4) = hh(4,1)
    hh(2,4) = hh(4,2)
    hh(3,4) = hh(4,3)
    hh(4,4) = (nn**2*(5.0_dp*nn**2-3.0_dp)**2*sk(1))/4.0_dp-3.0_dp&
        &/8.0_dp*(-1.0_dp+nn**2)*((5.0_dp*nn**2-1.0_dp)**2)*sk(2)&
        &+15.0_dp/4.0_dp*(nn**2)*((-1.0_dp+nn**2)**2)*sk(3)-5.0_dp&
        &/8.0_dp*(-1.0_dp+nn**2)*(-2.0_dp*nn**2+1.0_dp+nn**4)*sk(4)
    hh(5,4) = sqrt(2.0_dp)*ll*nn*(5.0_dp*nn**2-3.0_dp)*sqrt(3.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sk(1)/8.0_dp-(15.0_dp*nn**2-11.0_dp)&
        &*nn*ll*(5.0_dp*nn**2-1.0_dp)*sqrt(3.0_dp)*sqrt(2.0_dp)&
        &*sk(2)/16.0_dp+5.0_dp/8.0_dp*(3.0_dp*nn**3&
        &+3.0_dp*nn**2-nn-1.0_dp)*(-1.0_dp+nn)*sqrt(2.0_dp)*ll*nn*sqrt(3.0_dp)&
        &*sk(3)-5.0_dp/16.0_dp*nn*sqrt(3.0_dp)*ll*((-1.0_dp+nn**2)**2)&
        &*sqrt(2.0_dp)*sk(4)
    hh(6,4) = (2.0_dp*ll**2-1.0_dp+nn**2)*(nn**2)*(5.0_dp*nn**2-3.0_dp)&
        &*sqrt(15.0_dp)*sk(1)/4.0_dp-(3.0_dp*nn**2-1.0_dp)&
        &*(2.0_dp*ll**2-1.0_dp+nn**2)*(5.0_dp*nn**2-1.0_dp)*sqrt(15.0_dp)&
        &*sk(2)/8.0_dp+sqrt(15.0_dp)*(nn**2)*(2.0_dp*ll**2-1.0_dp+nn**2)&
        &*(3.0_dp*nn**2-1.0_dp)*sk(3)/4.0_dp-sqrt(5.0_dp)&
        &*(2.0_dp*ll**2-1.0_dp+nn**2)*sqrt(3.0_dp)*(-1.0_dp+nn**4)*sk(4)&
        &/8.0_dp
    hh(7,4) = sqrt(2.0_dp)*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*nn&
        &*(5.0_dp*nn**2-3.0_dp)*sqrt(5.0_dp)*sk(1)/8.0_dp-3.0_dp&
        &/16.0_dp*sqrt(2.0_dp)*(5.0_dp*nn**2-1.0_dp)*ll*(4.0_dp*ll**2-3.0_dp&
        &+3.0_dp*nn**2)*sqrt(5.0_dp)*nn*sk(2)+3.0_dp/8.0_dp&
        &*sqrt(10.0_dp)*nn*ll*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*(nn**2+1.0_dp)&
        &*sk(3)-sqrt(2.0_dp)*sqrt(5.0_dp)*ll*(4.0_dp*ll**2-3.0_dp&
        &+3.0_dp*nn**2)*nn*(nn**2+3.0_dp)*sk(4)/16.0_dp
    hh(1,5) = hh(5,1)
    hh(2,5) = hh(5,2)
    hh(3,5) = hh(5,3)
    hh(4,5) = hh(5,4)
    hh(5,5) = 3.0_dp/8.0_dp*(ll**2)*((5.0_dp*nn**2-1.0_dp)**2)*sk(1)&
        &+(-225.0_dp/ 16.0_dp*(ll**2)*(nn**4)+25.0_dp/16.0_dp*(nn**4)&
        &+65.0_dp/8.0_dp*(ll**2)*(nn**2)-5.0_dp/8.0_dp*(nn**2)+1.0_dp/16.0_dp&
        &-(ll**2)/16.0_dp)*sk(2)+5.0_dp/8.0_dp*(-1.0_dp+nn)&
        &*(9.0_dp*ll**2*nn**3-4.0_dp*nn**3+9.0_dp*ll**2*nn**2&
        &-4.0_dp*nn**2-ll**2*nn&
        &-ll**2)*sk(3)- 15.0_dp/16.0_dp*(ll**2*nn**2+1.0_dp-nn**2-ll**2)&
        &*(-1.0_dp+nn**2)*sk(4)
    hh(6,5) = 3.0_dp/ 8.0_dp*ll*(2.0_dp*ll**2-1.0_dp+nn**2)*&
        &(5.0_dp*nn**2-1.0_dp)*sqrt(10.0_dp)*nn*sk(1)-(45.0_dp*nn**4&
        &+90.0_dp*ll**2*nn**2-68.0_dp*nn**2+15.0_dp-26.0_dp*ll**2)*nn*ll&
        &*sqrt(10.0_dp)*sk(2)/ 16.0_dp+sqrt(10.0_dp)*ll&
        &*(9.0_dp*nn**4-22.0_dp*nn**2+18.0_dp*ll**2*nn**2+9.0_dp-10.0_dp*ll**2)&
        &*nn*sk(3)/ 8.0_dp-3.0_dp/ 16.0_dp*ll*(-1.0_dp+nn)*sqrt(5.0_dp)&
        &*sqrt(2.0_dp)*(nn**3+nn**2+2.0_dp*ll**2*nn-3.0_dp*nn+2.0_dp*ll**2&
        &-3.0_dp)*nn*sk(4)
    hh(7,5) = (ll**2)*(4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)*sqrt(3.0_dp)&
        &*(5.0_dp*nn**2-1.0_dp)*sqrt(5.0_dp)*sk(1)/8.0_dp&
        &-sqrt(15.0_dp)*(-5.0_dp*nn**4+45.0_dp*ll**2*nn**4-58.0_dp*ll**2*nn**2&
        &+6.0_dp*nn**2+60.0_dp*ll**4*nn**2-1.0_dp+5.0_dp*ll**2-4.0_dp*ll**4)&
        &*sk(2)/16.0_dp+sqrt(15.0_dp)*(-4.0_dp*nn**4+9.0_dp*ll**2*nn**4&
        &-14.0_dp*ll**2*nn**2+12.0_dp*ll**4*nn**2+4.0_dp*nn**2-3.0_dp*ll**2&
        &+4.0_dp*ll**4)*sk(3)/8.0_dp-sqrt(15.0_dp)&
        &*(3.0_dp*ll**2.0_dp*nn**4-3.0_dp*nn**4+2*nn**2+4.0_dp*ll**4*nn**2&
        &-6.0_dp*ll**2*nn**2+1.0_dp-5.0_dp*ll**2+4.0_dp*ll**4)*sk(4)&
        &/16.0_dp
    hh(1,6) = hh(6,1)
    hh(2,6) = hh(6,2)
    hh(3,6) = hh(6,3)
    hh(4,6) = hh(6,4)
    hh(5,6) = hh(6,5)
    hh(6,6) = 15.0_dp/4.0_dp*((2.0_dp*ll**2-1.0_dp+nn**2)**2)*(nn**2)&
        &*sk(1)+(- 45.0_dp/8.0_dp*(nn**6)-45.0_dp/2.0_dp*(ll**2)*(nn**4)&
        &+75.0_dp/8.0_dp*(nn**4)+(25.0_dp*ll**2*nn**2)-35.0_dp/8.0_dp*(nn**2)&
        &-45.0_dp/2.0_dp*(ll**4)*(nn**2)-5.0_dp/2.0_dp*(ll**2)+5.0_dp/8.0_dp&
        &+5.0_dp/2.0_dp*(ll**4))*sk(2)+(9.0_dp/4.0_dp*(nn**6)&
        &+(9.0_dp*ll**2*nn**4)-3.0_dp/2.0_dp*(nn**4)-(13.0_dp*ll**2*nn**2)&
        &+(nn**2)/4.0_dp+(9.0_dp*ll**4*nn**2)-(4.0_dp*ll**4)+(4.0_dp*ll**2))&
        &*sk(3)+(- 3.0_dp/8.0_dp*(nn**6)-3.0_dp/2.0_dp*(ll**2)*(nn**4)&
        &-3.0_dp/8.0_dp*(nn**4)+ 3.0_dp/8.0_dp*(nn**2)+(3.0_dp*ll**2*nn**2)&
        &-3.0_dp/2.0_dp*(ll**4)*(nn**2)+3.0_dp/2.0_dp*(ll**4)-3.0_dp/2.0_dp&
        &*(ll**2)+3.0_dp/8.0_dp)*sk(4)
    hh(7,6) = 5.0_dp/8.0_dp*(2.0_dp*ll**2-1.0_dp+nn**2)*ll*(4.0_dp*ll**2&
        &-3.0_dp+3.0_dp*nn**2)*sqrt(6.0_dp)*nn*sk(1)-5.0_dp/ 16.0_dp&
        &*sqrt(6.0_dp)*ll*nn*(9.0_dp*nn**4+30.0_dp*ll**2*nn**2-16.0_dp*nn**2&
        &-30.0_dp*ll**2&
        &+24.0_dp*ll**4+7.0_dp)*sk(2)+sqrt(6.0_dp)*ll*nn*(9.0_dp*nn**4&
        &+30.0_dp*ll**2*nn**2-10.0_dp*nn**2+24.0_dp*ll**4-30.0_dp*ll**2+5.0_dp)&
        &*sk(3)/8.0_dp-sqrt(6.0_dp)*ll*(3.0_dp*nn**4+10.0_dp*ll**2*nn**2&
        &+8.0_dp*ll**4+5.0_dp-10.0_dp*ll**2)*nn*sk(4)/16.0_dp
    hh(1,7) = hh(7,1)
    hh(2,7) = hh(7,2)
    hh(3,7) = hh(7,3)
    hh(4,7) = hh(7,4)
    hh(5,7) = hh(7,5)
    hh(6,7) = hh(7,6)
    hh(7,7) = 5.0_dp/8.0_dp*(ll**2)*((4.0_dp*ll**2-3.0_dp+3.0_dp*nn**2)**2)&
        &*sk(1)+(-135.0_dp/16.0_dp*(ll**2)*(nn**4)+15.0_dp/16.0_dp&
        &*(nn**4)-45.0_dp/2.0_dp*(ll**4)*(nn**2)-15.0_dp/8.0_dp*(nn**2)&
        &+135.0_dp/8.0_dp*(ll**2)*(nn**2)+45.0_dp/2.0_dp*(ll**4)&
        &+15.0_dp/16.0_dp-135.0_dp/16.0_dp*(ll**2)-(15.0_dp*ll**6))*sk(2)&
        &+(27.0_dp/8.0_dp*(ll**2)*(nn**4)-3.0_dp/2.0_dp*(nn**4)+3.0_dp/2.0_dp&
        &*(nn**2)+(9.0_dp*ll**4*nn**2)-27.0_dp/4.0_dp*(ll**2)*(nn**2)+27.0_dp&
        &/8.0_dp*(ll**2)-(9.0_dp*ll**4)+(6.0_dp*ll**6))*sk(3)&
        &+(9.0_dp/16.0_dp*(nn**4)-9.0_dp/ 16.0_dp*(ll**2)*(nn**4)-3.0_dp/2.0_dp&
        &*(ll**4)*(nn**2)+3.0_dp/8.0_dp*(nn**2)+9.0_dp/8.0_dp*(ll**2)*(nn**2)&
        &+3.0_dp/2.0_dp*(ll**4)+1.0_dp/16.0_dp-9.0_dp/16.0_dp*(ll**2)-(ll**6))&
        &*sk(4)

  end subroutine ff

end module sk
