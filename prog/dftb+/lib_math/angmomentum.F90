!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Angular momentum related routines
module angmomentum
  use assert
  use accuracy, only : dp
  use qm, only : unitary
  use commontypes, only : TOrbitals

  implicit none

  private
  public :: Loperators, getL

  ! construct Lz and L+ in the tesseral spherical hamonics basis for a given
  ! value of l
  interface Loperators
    module procedure operators
  end interface

  ! Calculate atomic angular momentum for each shell
  interface getL
    module procedure onsite
    module procedure dual
  end interface

contains

  !!* Returns L+ and Lz in the tesseral spherical Harmonics basis
  !!* used in DFTB+
  !!* @param Lplus L+ operator
  !!* @param Lz Lz operator
  !!* @param l value of the orbital momentum to construct these matrices
  subroutine operators(Lplus,Lz,l)
    complex(dp),intent(out) :: Lplus(0:,0:)
    complex(dp),intent(out) :: Lz(0:,0:)
    integer, intent(in)     :: l

    integer :: m ! magnetic quantum number
    complex(dp), parameter :: i = (0.0_dp,1.0_dp)
    complex(dp), allocatable :: u(:,:)

    @:ASSERT(l >= 0)
    @:ASSERT(all(shape(Lplus)==shape(Lz)))
    @:ASSERT(size(Lplus,dim=1)==2*l+1)
    @:ASSERT(size(Lplus,dim=2)==2*l+1)

    ! Lz in usual spherical harmonic basis
    Lz = 0.0_dp
    do m = -l, l
      Lz(l+m,l+m) = real(m,dp)
    end do

    ! L+ in usual spherical harmonic basis
    Lplus = 0.0_dp
    do m = -l, l-1
      Lplus(l+m+1,l+m) = sqrt(real(l*(l+1)-m*(m+1),dp))
    end do

    allocate(u(0:2*l,0:2*l))

    ! unitary transformation from $Y_{lm}$ to $\overline{Y}_{lm}$
    u(:,:) = 0.0_dp
    do m = 1, l
      u(l+m,l+m) = sqrt(0.5_dp) * real(mod(m+1,2)-mod(m,2),dp)
      u(l+m,l-m) = sqrt(0.5_dp) * 1.0_dp
      u(l-m,l+m) = -sqrt(0.5_dp) * i * real(mod(m,2)-mod(m+1,2),dp)
      u(l-m,l-m) = -sqrt(0.5_dp) * i
    end do
    u(l,l) = 1.0_dp

    call unitary(Lz,u)
    call unitary(Lplus,u)

  end subroutine operators

  !!* Calculates the on-site orbital angular momentum
  !!* @param Lshell resulting orbital angular momentum
  !!* @param iAtomStart Offset array in the square matrix.
  !!* @param orb Information about the orbitals in the system.
  !!* @param species Species of the atoms
  subroutine onsite(Lshell, rho, iAtomStart, orb, species)
    real(dp), intent(out)       :: Lshell(:,:,:)
    complex(dp), intent(in)     :: rho(:,:)
    integer,  intent(in)        :: iAtomStart(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: species(:)

    integer :: nAtom, nSpecies, nOrb, iSp
    integer :: ii, jj, kk, ll, iStart, iEnd
    complex(dp), allocatable :: SpeciesL(:,:,:,:)
    complex(dp), allocatable :: L(:,:,:)
    complex(dp), allocatable :: Lplus(:,:)
    complex(dp), allocatable :: tmpBlock(:,:)

    complex(dp), parameter :: i = (0.0_dp,1.0_dp)

    nAtom = size(Lshell,dim=3)
    nSpecies = maxval(species(1:nAtom))
    nOrb = size(rho,dim=1)

    @:ASSERT(size(rho, dim=1) == size(rho, dim=2))
    @:ASSERT(size(iAtomStart) == nAtom+1)
    @:ASSERT(mod(nOrb,2)==0)
    nOrb = nOrb / 2

    allocate(SpeciesL(orb%mOrb,orb%mOrb,3,nSpecies))
    SpeciesL = 0.0_dp
    allocate(L(orb%mOrb,orb%mOrb,3))
    allocate(Lplus(orb%mOrb,orb%mOrb))
    allocate(tmpBlock(orb%mOrb,orb%mOrb))
    do iSp = 1, nSpecies
      do jj = 1, orb%nShell(iSp)
        L = 0.0_dp
        Lplus = 0.0_dp
        tmpBlock = 0.0_dp
        kk = orb%angShell(jj,iSp)
        iStart = orb%posShell(jj,iSp)
        iEnd = orb%posShell(jj+1,iSp)-1
        call loperators(Lplus(1:2*kk+1,1:2*kk+1),tmpBlock(1:2*kk+1,1:2*kk+1),kk)
        L(iStart:iEnd,iStart:iEnd,3) = tmpBlock(1:2*kk+1,1:2*kk+1)
        tmpBlock(1:2*kk+1,1:2*kk+1) = transpose(conjg(Lplus(1:2*kk+1,1:2*kk+1)))
        L(iStart:iEnd,iStart:iEnd,1) = 0.5_dp * &
            & (Lplus(1:2*kk+1,1:2*kk+1) + tmpBlock(1:2*kk+1,1:2*kk+1))
        L(iStart:iEnd,iStart:iEnd,2) = 0.5_dp * i * &
            & (tmpBlock(1:2*kk+1,1:2*kk+1) - Lplus(1:2*kk+1,1:2*kk+1))
        SpeciesL(iStart:iEnd,iStart:iEnd,1:3,iSp) = &
            & L(iStart:iEnd,iStart:iEnd,1:3)
      end do
    end do

    Lshell = 0.0_dp

    do ii = 1, nAtom
      iSp = species(ii)
      jj = orb%nOrbSpecies(iSp)

      ! I block
      tmpBlock = 0.0_dp
      tmpBlock(1:jj,1:jj) = 0.5_dp * ( rho(iAtomStart(ii):iAtomStart(ii+1)-1, &
          & iAtomStart(ii):iAtomStart(ii+1)-1) &
          & + rho(nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1, &
          & nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1) )
      do ll = 1, orb%nOrbSpecies(iSp)
        tmpBlock(ll,ll+1:) = conjg(tmpBlock(ll+1:,ll)) ! Hermitize
      end do
      do ll = 1, orb%nShell(iSp)
        iStart = orb%posShell(ll,iSp)
        iEnd = orb%posShell(ll+1,iSp)-1
        do kk = 1, 3
          Lshell(kk,ll,ii) = Lshell(kk,ll,ii) + &
              & real(sum(SpeciesL(iStart:iEnd,iStart:iEnd,kk,iSp)* &
              & transpose(tmpBlock(iStart:iEnd,iStart:iEnd))), dp)
        end do
      end do

    end do

  end subroutine onsite

  !!* Calculates the on-site orbital angular momentum for dual populations
  !!* @param Lshell resulting orbital angular momentum
  !!* @param qBlockSkew Antisymmetric Mulliken block populations for imaginary
  !!* coefficients of Pauli matrics
  !!* @param orb Information about the orbitals in the system.
  !!* @param species Species of the atoms
  subroutine dual(Lshell, qBlockSkew, orb, species)
    real(dp), intent(out)       :: Lshell(:,:,:)
    real(dp), intent(in)        :: qBlockSkew(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: species(:)

    integer :: nAtom, nSpecies, iSp
    integer :: ii, jj, kk, ll, mm, iStart, iEnd
    real(dp), allocatable :: SpeciesL(:,:,:,:)
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    real(dp), allocatable :: tmpBlock(:,:)

    complex(dp), parameter :: i = (0.0_dp,1.0_dp)

    nAtom = size(LShell,dim=3)
    nSpecies = maxval(species(1:nAtom))

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

    allocate(tmpBlock(orb%mOrb,orb%mOrb))

    Lshell = 0.0_dp
    do ii = 1, nAtom
      iSp = species(ii)
      mm = orb%nOrbSpecies(iSp)
      tmpBlock(:,:) = 0.0_dp
      tmpBlock(1:mm,1:mm) = qBlockSkew(1:mm,1:mm,ii,1)! identity part
      do jj = 1, orb%nShell(iSp)
        iStart = orb%posShell(jj,iSp)
        iEnd = orb%posShell(jj+1,iSp)-1
        do ll = 1, 3
          Lshell(ll,jj,ii) = Lshell(ll,jj,ii) &
              & - sum(SpeciesL(iStart:iEnd,iStart:iEnd,ll,iSp) &
              &  * transpose(tmpBlock(iStart:iEnd,iStart:iEnd)))
        end do
      end do
    end do

  end subroutine dual

end module angmomentum
