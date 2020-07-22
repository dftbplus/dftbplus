!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines for spin orbit coupling
module dftbp_spinorbit
#:if WITH_SCALAPACK
  use dftbp_scalapackfx
#:endif
  use dftbp_environment
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_constants, only : imag
  use dftbp_angmomentum, only : getLOperators
  use dftbp_commontypes, only : TOrbitals
  use dftbp_densedescr
  implicit none
  private

  public :: getOnsiteSpinOrbitEnergy, addOnsiteSpinOrbitHam
  public :: getDualSpinOrbitEnergy, getDualSpinOrbitShift


contains


  !> Calculates the spin orbit energy for on-site L.S coupling
  subroutine getOnsiteSpinOrbitEnergy(env, Eatom, rho, denseDesc, xi, orb, species)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> returned energy for each atom
    real(dp), intent(out) :: Eatom(:)

    !> Density matrix in dense format
    complex(dp), intent(in) :: rho(:,:)

    !> Offset array in the square matrix.
    type(TDenseDescr), intent(in) :: denseDesc

    !> spin orbit constants for each shell of each species
    real(dp), intent(in) :: xi(:,:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms
    integer, intent(in) :: species(:)

    integer :: nAtom, nSpecies, nOrb, nOrbSp
    integer :: iSp, iAt, iOrb, iOrbStart, iOrbEnd
    complex(dp), allocatable :: speciesZ(:,:,:)
    complex(dp), allocatable :: speciesPlus(:,:,:)
    complex(dp), allocatable :: tmpBlock(:,:)

    nAtom = size(Eatom, dim=1)
    nSpecies = maxval(species(1:nAtom))
    nOrb = orb%nOrb

    @:ASSERT(size(denseDesc%iAtomStart) == nAtom + 1)
    @:ASSERT(size(xi, dim=2) == nSpecies)
    @:ASSERT(size(xi, dim=1) == orb%mShell)
    @:ASSERT(denseDesc%iAtomStart(nAtom + 1) == nOrb + 1)

    allocate(speciesZ(orb%mOrb, orb%mOrb, nSpecies))
    allocate(speciesPlus(orb%mOrb, orb%mOrb, nSpecies))
    do iSp = 1, nSpecies
      call getLSOperatorsForSpecies(orb, xi, iSp, speciesZ(:,:,iSp), speciesPlus(:,:, iSp))
    end do

    Eatom(:) = 0.0_dp

    allocate(tmpBlock(orb%mOrb, orb%mOrb))
    do iAt = 1, nAtom
      iSp = species(iAt)
      nOrbSp = orb%nOrbSpecies(iSp)
      iOrbStart = denseDesc%iAtomStart(iAt)
      iOrbEnd = denseDesc%iAtomStart(iAt + 1) - 1

      ! uu block
    #:if WITH_SCALAPACK
      call scalafx_cpg2l(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, iOrbStart, iOrbStart, rho,&
          & tmpBlock(1:nOrbSp, 1:nOrbSp))
    #:else
      tmpBlock(1:nOrbSp, 1:nOrbSp) = rho(iOrbStart:iOrbEnd, iOrbStart:iOrbEnd)
    #:endif
      tmpBlock(:,:) = 0.5_dp * tmpBlock
      do iOrb = 1, orb%nOrbSpecies(iSp)
        tmpBlock(iOrb, iOrb + 1 :) = conjg(tmpBlock(iOrb + 1 :, iOrb))
      end do
      Eatom(iAt) = Eatom(iAt) + real(sum(speciesZ(1:nOrbSp, 1:nOrbSp, iSp)&
          & * conjg(tmpBlock(1:nOrbSp, 1:nOrbSp))))

      ! dd block
    #:if WITH_SCALAPACK
      call scalafx_cpg2l(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, nOrb + iOrbStart,&
          & nOrb + iOrbStart, rho, tmpBlock(1:nOrbSp, 1:nOrbSp))
    #:else
      tmpBlock(1:nOrbSp, 1:nOrbSp) = rho(nOrb + iOrbStart : nOrb + iOrbEnd,&
          & nOrb + iOrbStart : nOrb + iOrbEnd)
    #:endif
      tmpBlock(:,:) = 0.5_dp * tmpBlock
      do iOrb = 1, orb%nOrbSpecies(iSp)
        tmpBlock(iOrb, iOrb + 1 :) = conjg(tmpBlock(iOrb + 1 :, iOrb))
      end do
      Eatom(iAt) = Eatom(iAt)&
          & - real(sum(speciesZ(1:nOrbSp, 1:nOrbSp, iSp) * conjg(tmpBlock(1:nOrbSp, 1:nOrbSp))))

      ! ud block
      ! two ud/du blocks so omit 0.5 factor
    #:if WITH_SCALAPACK
      call scalafx_cpg2l(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, nOrb + iOrbStart, iOrbStart,&
          & rho, tmpBlock(1:nOrbSp, 1:nOrbSp))
    #:else
      tmpBlock(1:nOrbSp, 1:nOrbSp) = &
          & rho(nOrb + iOrbStart : nOrb + iOrbEnd, iOrbStart : iOrbEnd)
    #:endif
      Eatom(iAt) = Eatom(iAt)&
          & + real(sum(speciesPlus(1:nOrbSp , 1:nOrbSp, iSp)&
          & * conjg(tmpBlock(1:nOrbSp, 1:nOrbSp))), dp)
    end do

  end subroutine getOnsiteSpinOrbitEnergy


  !> Adds spin-orbit contribution to dense Hamiltonian (for non-dual spin-orbit model).
  subroutine addOnsiteSpinOrbitHam(env, xi, species, orb, denseDesc, HSqrCplx)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Spin orbit constants for each species
    real(dp), intent(in) :: xi(:,:)

    !> chemical species
    integer, intent(in) :: species(:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> index array for atomic blocks in dense matrices
    type(TDenseDescr), intent(in) :: denseDesc

    !> Dense hamiltonian matrix (2 component)
    complex(dp), intent(inout) :: HSqrCplx(:,:)

    complex(dp), allocatable :: speciesZ(:,:,:), speciesPlus(:,:,:)
    integer :: nAtom, nSpecies, nOrb, nOrbSp, iOrbStart, iOrbEnd
    integer :: iSp, iAt

    nAtom = size(orb%nOrbAtom)
    nSpecies = maxval(species(1:nAtom))
    nOrb = orb%nOrb
    allocate(speciesZ(orb%mOrb, orb%mOrb, nSpecies))
    allocate(speciesPlus(orb%mOrb, orb%mOrb, nSpecies))

    do iSp = 1, nSpecies
      call getLSOperatorsForSpecies(orb, xi, iSp, speciesZ(:,:,iSp), speciesPlus(:,:,iSp))
    end do

    do iAt = 1, nAtom
      iSp = species(iAt)
      nOrbSp = orb%nOrbSpecies(iSp)
      iOrbStart = denseDesc%iAtomStart(iAt)
    #:if WITH_SCALAPACK
      call scalafx_addl2g(env%blacs%orbitalGrid, speciesZ(1:nOrbSp, 1:nOrbSp, iSp),&
          & denseDesc%blacsOrbSqr, iOrbStart, iOrbStart, HSqrCplx)
      call scalafx_addl2g(env%blacs%orbitalGrid, -speciesZ(1:nOrbSp, 1:nOrbSp, iSp),&
          & denseDesc%blacsOrbSqr, nOrb + iOrbStart, nOrb + iOrbStart, HSqrCplx)
      call scalafx_addl2g(env%blacs%orbitalGrid, speciesPlus(1:nOrbSp, 1:nOrbSp, iSp),&
          & denseDesc%blacsOrbSqr, nOrb + iOrbStart, iOrbStart, HSqrCplx)
      ! other triangle
      call scalafx_addl2g(env%blacs%orbitalGrid,&
          & transpose(conjg(speciesPlus(1:nOrbSp, 1:nOrbSp, iSp))),&
          & denseDesc%blacsOrbSqr, iOrbStart, nOrb + iOrbStart, HSqrCplx)
    #:else
      iOrbEnd = denseDesc%iAtomStart(iAt + 1) - 1
      HSqrCplx(iOrbStart:iOrbEnd, iOrbStart:iOrbEnd) = &
          & HSqrCplx(iOrbStart:iOrbEnd, iOrbStart:iOrbEnd) + speciesZ(1:nOrbSp, 1:nOrbSp, iSp)
      HSqrCplx(nOrb + iOrbStart : nOrb + iOrbEnd, nOrb + iOrbStart : nOrb + iOrbEnd) = &
          & HSqrCplx(nOrb + iOrbStart : nOrb + iOrbEnd, nOrb + iOrbStart : nOrb + iOrbEnd) &
          & - speciesZ(1:nOrbSp, 1:nOrbSp, iSp)
      HSqrCplx(nOrb + iOrbStart : nOrb + iOrbEnd, iOrbStart:iOrbEnd) = &
          & HSqrCplx(nOrb + iOrbStart : nOrb + iOrbEnd, iOrbStart:iOrbEnd)&
          & + speciesPlus(1:nOrbSp, 1:nOrbSp, iSp)
    #:endif
    end do

  end subroutine addOnsiteSpinOrbitHam


  !> Calculates the spin orbit energy and angular momentum for dual L.S coupling
  subroutine getDualSpinOrbitEnergy(Eatom, qBlockSkew, xi, orb, species)

    !> returned energy for each atom
    real(dp), intent(out) :: Eatom(:)

    !> Antisymmetric Mulliken block populations for imaginary coefficients of Pauli matrics
    real(dp), intent(in) :: qBlockSkew(:,:,:,:)

    !> spin orbit constants for each shell of each species
    real(dp), intent(in) :: xi(:,:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms
    integer, intent(in) :: species(:)

    integer :: nAtom, nSpecies
    integer :: iAt, iSp, nOrbSp
    complex(dp) :: speciesZ(orb%mOrb, orb%mOrb), speciesPlus(orb%mOrb, orb%mOrb)
    real(dp), allocatable :: speciesL(:,:,:,:)
    complex(dp), allocatable :: tmpBlock(:,:)

    nAtom = size(Eatom,dim=1)
    nSpecies = maxval(species(1:nAtom))
    @:ASSERT(size(xi,dim=2) == nSpecies)
    @:ASSERT(size(xi,dim=1) == orb%mShell)

    allocate(speciesL(orb%mOrb, orb%mOrb, 3, nSpecies))
    do iSp = 1, nSpecies
      call getLSOperatorsForSpecies(orb, xi, iSp, speciesZ, speciesPlus)
      speciesL(:, :, 1, iSp) = aimag(speciesPlus)
      speciesL(:, :, 2, iSp) = -real(speciesPlus)
      speciesL(:, :, 3, iSp) = aimag(speciesZ)
    end do

    Eatom(:) = 0.0_dp
    allocate(tmpBlock(orb%mOrb,orb%mOrb))
    do iAt = 1, nAtom
      iSp = species(iAt)
      nOrbSp = orb%nOrbSpecies(iSp)

      ! Lz.Sz
      tmpBlock(:,:) = 0.0_dp
      tmpBlock(1:nOrbSp, 1:nOrbSp) = qBlockSkew(1:nOrbSp, 1:nOrbSp, iAt, 4)
      Eatom(iAt) = Eatom(iAt) - real(sum(transpose(tmpBlock) * speciesL(:, :, 3, iSp)))

      ! (Lx.Sx + Ly.Sy).
      tmpBlock(:,:) = 0.0_dp
      tmpBlock(1:nOrbSp, 1:nOrbSp) =&
          & qBlockSkew(1:nOrbSp, 1:nOrbSp, iAt, 3) - imag * qBlockSkew(1:nOrbSp, 1:nOrbSp, iAt, 2)
      Eatom(iAt) = Eatom(iAt)&
          & - real(sum(transpose(tmpBlock) * (imag * speciesL(:,:,1,iSp) + speciesL(:,:,2,iSp))))
    end do

  end subroutine getDualSpinOrbitEnergy


  !> Constructs shift potential for spin-orbit
  subroutine getDualSpinOrbitShift(shift, xi, orb, species)

    !> block shift from the potential
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> spin orbit constants for each shell of each species
    real(dp), intent(in) :: xi(:,:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms
    integer, intent(in) :: species(:)

    integer :: nAtom, nSpecies
    integer :: iSp, iSpin
    real(dp), allocatable :: tmpShift(:,:,:,:)
    complex(dp) :: speciesZ(orb%mOrb, orb%mOrb), speciesPlus(orb%mOrb, orb%mOrb)

    @:ASSERT(size(shift,dim=1)==orb%mOrb)
    @:ASSERT(size(shift,dim=2)==orb%mOrb)
    nAtom = size(shift,dim=3)
    @:ASSERT(size(shift,dim=4)==4)
    nSpecies = maxval(species(1:nAtom))
    @:ASSERT(size(species)>=nAtom)
    @:ASSERT(size(xi,dim=2) == nSpecies)
    @:ASSERT(size(xi,dim=1) == orb%mShell)

    allocate(tmpShift(orb%mOrb, orb%mOrb, nSpecies, 4))

    tmpShift(:,:,:,:) = 0.0_dp
    do iSp = 1, nSpecies
      call getLSOperatorsForSpecies(orb, xi, iSp, speciesZ, speciesPlus)
      tmpShift(:, :, iSp, 4) = aimag(speciesZ)
      tmpShift(:, :, iSp, 2) = aimag(speciesPlus)
      tmpShift(:, :, iSp, 3) = -real(speciesPlus)
    end do

    shift(:,:,:,:) = 0.0_dp
    do iSpin = 2, 4
      do iSp = 1, nAtom
        shift(:, :, iSp, iSpin) = tmpShift(:, :, species(iSp), iSpin)
      end do
    end do

  end subroutine getDualSpinOrbitShift


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !> Returns 1/2 xi * Lz and 1/2 xi * Lplus for a given species.
  subroutine getLSOperatorsForSpecies(orb, xi, iSpecies, speciesZ, speciesPlus)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Spin coupling constants
    real(dp), intent(in) :: xi(:,:)

    !> Species to get the operators for
    integer, intent(in) :: iSpecies

    !> Species specific Lz operator
    complex(dp), intent(out) :: speciesZ(:,:)

    !> Species specific L+ operator
    complex(dp), intent(out) :: speciesPlus(:,:)

    complex(dp) :: Lz(orb%mOrb, orb%mOrb), Lplus(orb%mOrb, orb%mOrb)
    integer :: iShell, ll, nOrbShell, iOrbStart, iOrbEnd

    @:ASSERT(all(shape(speciesZ) == [orb%mOrb, orb%mOrb]))
    @:ASSERT(all(shape(speciesZ) == shape(speciesPlus)))

    speciesZ(:,:) = 0.0_dp
    speciesPlus(:,:) = 0.0_dp
    do iShell = 1, orb%nShell(iSpecies)
      ll = orb%angShell(iShell, iSpecies)
      if (ll == 0) then
        cycle
      end if
      nOrbShell = 2 * ll + 1
      iOrbStart = orb%posShell(iShell, iSpecies)
      iOrbEnd = orb%posShell(iShell + 1, iSpecies) - 1
      call getLOperators(ll, Lplus(1:nOrbShell, 1:nOrbShell), Lz(1:nOrbShell, 1:nOrbShell))
      speciesZ(iOrbStart:iOrbEnd, iOrbStart:iOrbEnd) =&
          & 0.5_dp * xi(iShell, iSpecies) * Lz(1:nOrbShell, 1:nOrbShell)
      speciesPlus(iOrbStart:iOrbEnd, iOrbStart:iOrbEnd) =&
          & 0.5_dp * xi(iShell, iSpecies) * Lplus(1:nOrbShell, 1:nOrbShell)
    end do

  end subroutine getLSOperatorsForSpecies


end module dftbp_spinorbit
