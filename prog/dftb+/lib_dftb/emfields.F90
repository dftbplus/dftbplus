!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Module for external electromagnetic fields - currently scalar magnetic field
module dftbp_emfields
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_constants
  use dftbp_globalenv, only : stdOut
  use dftbp_angmomentum, only : Loperators
  use dftbp_simplealgebra, only : cross3
  use dftbp_commontypes, only : TOrbitals
  implicit none

  private
  public :: shiftB

  !!* The shift from the scalar potential part of the B field
  interface shiftB
    module procedure shiftB_
  end interface

contains

!  !!* @note Using symmetric gauge, A = -.5 r x B
!  subroutine MagField(H0,iH0,coords,img2centCell,orb,iPair,species, &
!      & BField,nAtom,nNeighbour,iNeighbour)
!    real(dp), intent(inout)     :: H0(:)
!    real(dp), intent(inout)     :: iH0(:)
!    real(dp), intent(in)        :: coords(:,:)
!    integer,  intent(in)        :: img2CentCell(:)
!    type(TOrbitals), intent(in) :: orb
!    integer,  intent(in)        :: iPair(0:,:)
!    integer,  intent(in)        :: species(:)
!    real(dp), intent(in)        :: BField(3)
!    integer, intent(in)         :: nAtom
!    integer, intent(in)         :: nNeighbour(:)
!    integer, intent(in)         :: iNeighbour(0:,:)
!
!    real(dp) :: tmpA(3,2), phase
!    integer :: iAt1, iAt2, iAt2f, iOrig, iSp1, iSp2, nOrb1, nOrb2
!    integer :: iNeigh
!
!    ASSERT(size(nNeighbour)==nAtom)
!    ASSERT(size(iNeighbour,dim=2)==nAtom)
!    ASSERT(size(species)>=maxval(iNeighbour))
!    ASSERT(size(species)<=size(img2CentCell))
!    ASSERT(size(iPair,dim=1)>=(maxval(nNeighbour)+1))
!    ASSERT(size(iPair,dim=2)==nAtom)
!
!    do iAt1 = 1, nAtom
!      iSp1 = species(iAt1)
!      nOrb1 = orb%nOrbSpecies(iSp1)
!            do iNeigh = 1, nNeighbour(iAt1)
!        iAt2 = iNeighbour(iNeigh, iAt1)
!        iAt2f = img2CentCell(iAt2)
!        iSp2 = species(iAt2f)
!        nOrb2 = orb%nOrbSpecies(iSp2)
!        iOrig = iPair(iNeigh, iAt1)
!        call cross3(tmpA(:,2),coords(:,iAt2)-coords(:,iAt1), &
!            & coords(:,iAt2)+coords(:,iAt1))
!        phase = -.025_dp * alpha_fs &
!            & *dot_product(tmpA(:,2),BField)
!        write(stdout, *)'phase :',phase
!        iH0(iOrig+1:iOrig+nOrb2*nOrb1) = sin(phase) * &
!            & H0(iOrig+1:iOrig+nOrb2*nOrb1)
!        H0(iOrig+1:iOrig+nOrb2*nOrb1) = cos(phase) * &
!            & H0(iOrig+1:iOrig+nOrb2*nOrb1)
!      end do
!    end do
!
!  end subroutine MagField


  !!* Constructs shift potential for scalar potential part of megnetic field
  !!* @param shift block shift from the potential
  !!* @param iShift imaginary block shift from the potential
  !!* @param BFieldStrength magnetic field strength - atomi CGS units
  !!* @param BfieldVector magnetic field direction
  !!* @param orb Angular momentum information about the orbitals
  !!* @param species list of the species for each atom
  subroutine shiftB_(shift, iShift, BFieldStrength, BfieldVector, orb, species)
    real(dp), intent(inout)     :: shift(:,:,:,:)
    real(dp), intent(inout)     :: iShift(:,:,:,:)
    real(dp), intent(in)        :: BFieldStrength
    real(dp), intent(in)        :: BfieldVector(3)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: species(:)

    integer :: iAt, nAtom, iSpin, nSpin, iSp, iSh, iOrb, nSpecies
    integer :: ii, jj, kk, ll, mm, iStart, iEnd
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    real(dp), allocatable :: SpeciesL(:,:,:,:)

    nAtom = size(shift,dim=3)
    nSpin = size(shift,dim=4)
    nSpecies = maxval(species(1:nAtom))

    ! spin Zeeman part
    select case(nSpin)
    case(2) ! z aligned electron spins
      do iAt = 1, nAtom
        iSp = species(iAt)
        do iSh = 1, orb%nShell(iSp)
          do iOrb = orb%posShell(iSh,iSp),orb%posShell(iSh+1,iSp)-1
            shift(iOrb,iOrb,iAt,2) = shift(iOrb,iOrb,iAt,2) - &
                & gfac * mu_B * BFieldStrength * BfieldVector(2)
          end do
        end do
      end do
    case(4)
      do iSpin = 2, nSpin
        do iAt = 1, nAtom
          iSp = species(iAt)
          do iSh = 1, orb%nShell(iSp)
            do iOrb = orb%posShell(iSh,iSp),orb%posShell(iSh+1,iSp)-1
              shift(iOrb,iOrb,iAt,iSpin) = shift(iOrb,iOrb,iAt,iSpin) -&
                  & gfac * mu_B * BFieldStrength * BfieldVector(iSpin-1)
            end do
          end do
        end do
      end do
    end select

    ! Orbital Zeeman part

    allocate(SpeciesL(orb%mOrb,orb%mOrb,3,nSpecies))
    SpeciesL = 0.0_dp
    allocate(Lz(orb%mOrb,orb%mOrb))
    allocate(Lplus(orb%mOrb,orb%mOrb))
    do ii = 1, nSpecies
      do jj = 1, orb%nShell(ii)
        Lz = 0.0_dp
        Lplus = 0.0_dp
        kk = orb%angShell(jj,ii)
        call loperators(Lplus(1:2*kk+1,1:2*kk+1),Lz(1:2*kk+1,1:2*kk+1),kk)
        speciesL(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1, &
            & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,1,ii) &
            & = aimag(Lplus(1:2*kk+1,1:2*kk+1))
        speciesL(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1, &
            & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,2,ii) &
            & = -real(Lplus(1:2*kk+1,1:2*kk+1))
        speciesL(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1, &
            & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,3,ii) &
            & = aimag(Lz(1:2*kk+1,1:2*kk+1))
      end do
    end do

    do ii = 1, nAtom
      iSp = species(ii)
      mm = orb%nOrbSpecies(iSp)
      do jj = 1, orb%nShell(iSp)
        iStart = orb%posShell(jj,iSp)
        iEnd = orb%posShell(jj+1,iSp)-1
        do ll = 1, 3
          iShift(iStart:iEnd,iStart:iEnd,ii,1) = &
              & iShift(iStart:iEnd,iStart:iEnd,ii,1) - &
              & mu_B * BFieldStrength*BfieldVector(ll) * &
              &  SpeciesL(iStart:iEnd,iStart:iEnd,ll,iSp)
        end do
      end do
    end do

  end subroutine shiftB_

end module dftbp_emfields
