!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Routines for spin orbit coupling
module spinorbit
  use assert
  use accuracy, only : dp
  use angmomentum, only : Loperators
  use commontypes, only : TOrbitals

  implicit none

  private
  public :: getEnergySpinOrbit, shiftLS

  !!* Interfaces for spin orbit energies in either the onsite (just local part
  !!* of the density matrix) or dual (Mulliken projected density matrix)
  interface  getEnergySpinOrbit
    module procedure onsite
    module procedure dual
  end interface

contains

  !!* Calculates the spin orbit energy for on-site L.S coupling
  !!* @param Eatom returned energy for each atom
  !!* @param rho Density matrix in Packed format
  !!* @param iAtomStart Offset array in the square matrix.
  !!* @param xi spin orbit constants for each shell of each species
  !!* @param orb Information about the orbitals in the system.
  !!* @param species Species of the atoms
  subroutine onsite(Eatom, rho, iAtomStart, xi, orb, species)
    real(dp), intent(out)       :: Eatom(:)
    complex(dp), intent(in)     :: rho(:,:)
    integer,  intent(in)        :: iAtomStart(:)
    real(dp), intent(in)        :: xi(:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: species(:)

    integer :: nAtom, nSpecies, nOrb
    integer :: ii, jj, kk, ll
    complex(dp), allocatable :: SpeciesZ(:,:,:)
    complex(dp), allocatable :: SpeciesPlus(:,:,:)
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    complex(dp), allocatable :: tmpBlock(:,:)

    nAtom = size(Eatom,dim=1)
    nSpecies = maxval(species(1:nAtom))
    nOrb = size(rho,dim=1)

    @:ASSERT(size(rho, dim=1) == size(rho, dim=2))
    @:ASSERT(size(iAtomStart) == nAtom+1)
    @:ASSERT(size(xi,dim=2) == nSpecies)
    @:ASSERT(size(xi,dim=1) == orb%mShell)
    @:ASSERT(mod(nOrb,2)==0)
    nOrb = nOrb / 2
    @:ASSERT(iAtomStart(nAtom+1)==nOrb+1)

    allocate(SpeciesZ(orb%mOrb,orb%mOrb,nSpecies))
    SpeciesZ = 0.0_dp
    allocate(SpeciesPlus(orb%mOrb,orb%mOrb,nSpecies))
    SpeciesPlus = 0.0_dp
    allocate(Lz(orb%mOrb,orb%mOrb))
    allocate(Lplus(orb%mOrb,orb%mOrb))
    allocate(tmpBlock(orb%mOrb,orb%mOrb))
    do ii = 1, nSpecies
      do jj = 1, orb%nShell(ii)
        Lz = 0.0_dp
        Lplus = 0.0_dp
        kk = orb%angShell(jj,ii)
        call loperators(Lplus(1:2*kk+1,1:2*kk+1),Lz(1:2*kk+1,1:2*kk+1),kk)
        SpeciesZ(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1, &
            & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,ii) &
            & = 0.5_dp*xi(jj,ii)*Lz(1:2*kk+1,1:2*kk+1)
        SpeciesPlus(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1, &
            & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,ii) &
            & = 0.5_dp*xi(jj,ii)*Lplus(1:2*kk+1,1:2*kk+1)
      end do
    end do

    Eatom = 0.0_dp

    do ii = 1, nAtom

      jj = species(ii)
      kk = orb%nOrbSpecies(jj)
      ! uu block
      tmpBlock(1:kk,1:kk) = 0.5_dp*rho(iAtomStart(ii):iAtomStart(ii+1)-1, &
          & iAtomStart(ii):iAtomStart(ii+1)-1)
      do ll = 1, orb%nOrbSpecies(jj)
        tmpBlock(ll,ll+1:)=conjg(tmpBlock(ll+1:,ll)) ! Hermitize
      end do
      Eatom(ii) = Eatom(ii)&
          & + real(sum(SpeciesZ(1:kk,1:kk,jj) * conjg(tmpBlock(1:kk,1:kk))))

      ! dd block
      tmpBlock(1:kk,1:kk) = &
          & 0.5_dp*rho(nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1, &
          & nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1)
      do ll = 1, orb%nOrbSpecies(jj)
        tmpBlock(ll,ll+1:)=conjg(tmpBlock(ll+1:,ll)) ! Hermitize
      end do
      Eatom(ii) = Eatom(ii)&
          & - real(sum(SpeciesZ(1:kk,1:kk,jj) * conjg(tmpBlock(1:kk,1:kk))))

      ! ud block
      tmpBlock(1:kk,1:kk) = & ! two ud/du blocks so omit 0.5 factor
          & rho(nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1, &
          & iAtomStart(ii):iAtomStart(ii+1)-1)
      Eatom(ii) = Eatom(ii)&
          & + real(sum(SpeciesPlus(1:kk,1:kk,jj) * conjg(tmpBlock(1:kk,1:kk))),&
          & dp)

    end do

  end subroutine onsite

  !!* Calculates the spin orbit energy and angular momentum for dual L.S
  !!* coupling
  !!* @param Eatom returned energy for each atom
  !!* @param qBlockSkew Antisymmetric Mulliken block populations for imaginary
  !!* coefficients of Pauli matrics
  !!* @param xi spin orbit constants for each shell of each species
  !!* @param orb Information about the orbitals in the system.
  !!* @param species Species of the atoms
  subroutine dual(Eatom, qBlockSkew, xi, orb, species)
    real(dp), intent(out)       :: Eatom(:)
    real(dp), intent(in)        :: qBlockSkew(:,:,:,:)
    real(dp), intent(in)        :: xi(:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: species(:)

    real(dp) :: total
    integer :: nAtom, nSpecies, iSp
    integer :: ii, jj, kk
    real(dp), allocatable :: SpeciesL(:,:,:,:)
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    complex(dp), allocatable :: tmpBlock(:,:)

    complex(dp), parameter :: i = (0.0_dp,1.0_dp)

    nAtom = size(Eatom,dim=1)
    nSpecies = maxval(species(1:nAtom))
    @:ASSERT(size(xi,dim=2) == nSpecies)
    @:ASSERT(size(xi,dim=1) == orb%mShell)

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
            & = 0.5_dp*xi(jj,ii)*aimag(Lplus(1:2*kk+1,1:2*kk+1))
        speciesL(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1, &
            & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,2,ii) &
            & = -0.5_dp*xi(jj,ii)*real(Lplus(1:2*kk+1,1:2*kk+1))
        speciesL(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1, &
            & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,3,ii) &
            & = 0.5_dp*xi(jj,ii)*aimag(Lz(1:2*kk+1,1:2*kk+1))
      end do
    end do

    allocate(tmpBlock(orb%mOrb,orb%mOrb))

    Eatom = 0.0_dp

    do ii = 1, nAtom

      iSp = species(ii)
      jj = orb%nOrbSpecies(iSp)

      ! Lz.Sz
      tmpBlock(:,:) = 0.0_dp
      tmpBlock(1:jj,1:jj) = qBlockSkew(1:jj,1:jj,ii,4)

      total = 0.0_dp
      Eatom(ii) = Eatom(ii)&
          & - real(sum(transpose(tmpBlock) * SpeciesL(:,:,3,iSp)))

      ! (Lx.Sx + Ly.Sy).
      tmpBlock(:,:) = 0.0_dp
      tmpBlock(1:jj,1:jj) = (qBlockSkew(1:jj,1:jj,ii,3) &
          & -i * qBlockSkew(1:jj,1:jj,ii,2))
      total = 0.0_dp
      Eatom(ii) = Eatom(ii)&
          & - real(sum(transpose(tmpBlock)&
          & * (i * SpeciesL(:,:,1,iSp) + SpeciesL(:,:,2,iSp) )))

    end do

  end subroutine dual

  !!* Constructs shift potential for spin-orbit
  !!* @param shift block shift from the potential
  !!* @param xi spin orbit constants for each shell of each species
  !!* @param orb Information about the orbitals in the system.
  !!* @param species Species of the atoms
  subroutine shiftLS(shift, xi, orb, species)
    real(dp), intent(inout)       :: shift(:,:,:,:)
    real(dp), intent(in)        :: xi(:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: species(:)

    integer :: nAtom, nSpecies
    integer :: ii, jj, kk, iSpin
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    real(dp), allocatable    :: tmpShift(:,:,:,:)

    @:ASSERT(size(shift,dim=1)==orb%mOrb)
    @:ASSERT(size(shift,dim=2)==orb%mOrb)
    nAtom = size(shift,dim=3)
    @:ASSERT(size(shift,dim=4)==4)
    nSpecies = maxval(species(1:nAtom))
    @:ASSERT(size(species)>=nAtom)
    @:ASSERT(size(xi,dim=2) == nSpecies)
    @:ASSERT(size(xi,dim=1) == orb%mShell)

    allocate(tmpShift(orb%mOrb,orb%mOrb,nSpecies,4))
    allocate(Lz(orb%mOrb,orb%mOrb))
    allocate(Lplus(orb%mOrb,orb%mOrb))

    tmpShift(:,:,:,:) = 0.0_dp

    do ii = 1, nSpecies
      do jj = 1, orb%nShell(ii)
        Lz = 0.0_dp
        Lplus = 0.0_dp
        kk = orb%angShell(jj,ii)
        call loperators(Lplus(1:2*kk+1,1:2*kk+1),Lz(1:2*kk+1,1:2*kk+1),kk)
        tmpShift(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1, &
            & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,ii,4) &
            & = 0.5_dp*xi(jj,ii)*aimag(Lz(1:2*kk+1,1:2*kk+1))
        tmpShift(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1, &
            & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,ii,2) &
            & = 0.5_dp*xi(jj,ii)*aimag(Lplus(1:2*kk+1,1:2*kk+1))
        tmpShift(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1, &
            & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,ii,3) &
            & = -0.5_dp*xi(jj,ii)*real(Lplus(1:2*kk+1,1:2*kk+1))
      end do
    end do

    shift(:,:,:,:) = 0.0_dp
    do iSpin = 2, 4
      do ii = 1, nAtom
        shift(:,:,ii,iSpin) = tmpShift(:,:,species(ii),iSpin)
      end do
    end do

  end subroutine shiftLS

end module spinorbit
