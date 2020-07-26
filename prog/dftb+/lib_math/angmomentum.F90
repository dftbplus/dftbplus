!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Angular momentum related routines
module dftbp_angmomentum
#:if WITH_SCALAPACK
  use dftbp_scalapackfx
#:endif
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_constants, only : imag
  use dftbp_qm
  use dftbp_commontypes, only : TOrbitals
  use dftbp_environment
  use dftbp_densedescr
  implicit none
  private

  public :: getLOperators, getLOperatorsForSpecies, getLOnsite, getLDual, rotateZ

  !> Construct matrix for rotation of orbitals around the z axis in the tesseral spherical hamonics
  !> basis
  interface rotateZ
    module procedure zrot_onel
    module procedure zrot_manyl
  end interface rotateZ


contains


  !> Returns L^+ and L_z in the tesseral spherical Harmonics basis used in DFTB+
  subroutine getLOperators(ll, Lplus, Lz)

    !> value of the orbital momentum to construct these matrices
    integer, intent(in) :: ll

    !> L^+ operator
    complex(dp),intent(out) :: Lplus(0:, 0:)

    !> L_z operator
    complex(dp),intent(out) :: Lz(0:, 0:)

    ! magnetic quantum number
    integer :: mm

    complex(dp), allocatable :: uu(:,:)

    @:ASSERT(ll >= 0)
    @:ASSERT(all(shape(Lplus) == shape(Lz)))
    @:ASSERT(size(Lplus, dim=1) == 2 * ll + 1)
    @:ASSERT(size(Lplus, dim=2) == 2 * ll + 1)

    ! L_z in usual spherical harmonic basis
    Lz(:,:) = 0.0_dp
    do mm = -ll, ll
      Lz(ll + mm, ll + mm) = real(mm, dp)
    end do

    ! L^+ in usual spherical harmonic basis
    Lplus(:,:) = 0.0_dp
    do mm = -ll, ll - 1
      Lplus(ll + mm + 1, ll + mm) = sqrt(real(ll * (ll + 1) - mm * (mm + 1), dp))
    end do

    allocate(uu(0 : 2 * ll, 0 : 2 * ll))

    ! unitary transformation from Y_{lm} to \overline{Y}_{lm}
    uu(:,:) = 0.0_dp
    do mm = 1, ll
      uu(ll + mm, ll  + mm) = sqrt(0.5_dp) * real(mod(mm + 1, 2)-mod(mm, 2), dp)
      uu(ll + mm, ll - mm) = sqrt(0.5_dp)
      uu(ll - mm, ll + mm) = -sqrt(0.5_dp) * imag * real(mod(mm, 2) - mod(mm + 1, 2), dp)
      uu(ll - mm, ll - mm) = -sqrt(0.5_dp) * imag
    end do
    uu(ll, ll) = 1.0_dp

    call makeSimilarityTrans(Lz, uu)
    call makeSimilarityTrans(Lplus, uu)

  end subroutine getLOperators


  !> Returns Lz and Lplus for a given species.
  subroutine getLOperatorsForSpecies(orb, iSpecies, speciesZ, speciesPlus)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species to get the operators for
    integer, intent(in) :: iSpecies

    !> Species specific Lz operator
    complex(dp), intent(out) :: speciesZ(:,:)

    !> Species specific L+ operator
    complex(dp), intent(out) :: speciesPlus(:,:)

    integer :: iShell, ll, iOrbStart, iOrbEnd

    @:ASSERT(all(shape(speciesZ) == shape(speciesPlus)))

    speciesZ(:,:) = 0.0_dp
    speciesPlus(:,:) = 0.0_dp
    do iShell = 1, orb%nShell(iSpecies)
      ll = orb%angShell(iShell, iSpecies)
      iOrbStart = orb%posShell(iShell, iSpecies)
      iOrbEnd = orb%posShell(iShell + 1, iSpecies) - 1
      call getLOperators(ll, speciesPlus(iOrbStart:iOrbEnd, iOrbStart:iOrbEnd),&
          & speciesZ(iOrbStart:iOrbEnd, iOrbStart:iOrbEnd))
    end do

  end subroutine getLOperatorsForSpecies


  !> Calculates the on-site orbital angular momentum
  subroutine getLOnsite(env, Lshell, rho, denseDesc, orb, species)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> resulting orbital angular momentum (cartesian component, atomic shell, atom)
    real(dp), intent(out) :: Lshell(:,:,:)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms
    integer, intent(in) :: species(:)

    integer :: nAtom, nSpecies, nOrb, nOrbSp
    integer :: iSp, iAt, iShell, iOrb, iOrbStart, iOrbEnd, kk
    complex(dp), allocatable :: speciesL(:,:,:,:)
    complex(dp) :: tmpBlock(orb%mOrb, orb%mOrb)

    nAtom = size(Lshell, dim=3)
    nSpecies = maxval(species(1:nAtom))
    nOrb = size(rho, dim=1)

    @:ASSERT(size(denseDesc%iAtomStart) == nAtom + 1)
    @:ASSERT(mod(nOrb,2) == 0)
    nOrb = nOrb / 2

    allocate(speciesL(orb%mOrb, orb%mOrb, 3, nSpecies))
    speciesL(:,:,:,:) = 0.0_dp
    do iSp = 1, nSpecies
      call localGetLOperatorsForSpecies(orb, iSp, speciesL(:,:,:,iSp))
    end do

    Lshell(:,:,:) = 0.0_dp

    do iAt = 1, nAtom
      iSp = species(iAt)
      nOrbSp = orb%nOrbSpecies(iSp)
      iOrbStart = denseDesc%iAtomStart(iAt)
      iOrbEnd = denseDesc%iAtomStart(iAt + 1) - 1

      ! I block
      tmpBlock(:,:) = 0.0_dp
    #:if WITH_SCALAPACK
      call scalafx_addg2l(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, iOrbStart, iOrbStart, rho,&
          & tmpBlock(1:nOrbSp, 1:nOrbSp))
      call scalafx_addg2l(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, nOrb + iOrbStart,&
          & nOrb + iOrbStart, rho, tmpBlock(1:nOrbSp, 1:nOrbSp))
    #:else
      tmpBlock(1:nOrbSp, 1:nOrbSp) = rho(iOrbStart:iOrbEnd, iOrbStart:iOrbEnd) &
          & + rho(nOrb + iOrbStart : nOrb + iOrbEnd, nOrb + iOrbStart : nOrb + iOrbEnd)
    #:endif
      tmpBlock(:,:) = 0.5_dp * tmpBlock
      ! Hermitize
      do iOrb = 1, nOrbSp
        tmpBlock(iOrb, iOrb + 1 :) = conjg(tmpBlock(iOrb + 1 :, iOrb))
      end do
      do iShell = 1, orb%nShell(iSp)
        iOrbStart = orb%posShell(iShell, iSp)
        iOrbEnd = orb%posShell(iShell + 1, iSp) - 1
        do kk = 1, 3
          Lshell(kk, iShell, iAt) = Lshell(kk, iShell, iAt) - &
              & real(sum(speciesL(iOrbStart:iOrbEnd, iOrbStart:iOrbEnd, kk, iSp)&
              & * transpose(conjg(tmpBlock(iOrbStart:iOrbEnd, iOrbStart:iOrbEnd)))), dp)
        end do
      end do

    end do

  end subroutine getLOnsite


  !> Calculates the on-site orbital angular momentum for dual populations
  subroutine getLDual(Lshell, qBlockSkew, orb, species)

    !> resulting orbital angular momentum (cartesian component, atomic shell, atom)
    real(dp), intent(out) :: Lshell(:,:,:)

    !> Antisymmetric Mulliken block populations for imaginary coefficients of Pauli matrics
    real(dp), intent(in) :: qBlockSkew(:,:,:,:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms
    integer, intent(in) :: species(:)

    integer :: nAtom, nSpecies, nOrbSp
    integer :: iAt, iSp, iOrb, iOrbStart, iOrbEnd, kk
    complex(dp), allocatable :: L(:,:,:)
    real(dp), allocatable :: speciesL(:,:,:,:)
    real(dp), allocatable :: tmpBlock(:,:)

    complex(dp), parameter :: i = (0.0_dp,1.0_dp)

    nAtom = size(LShell, dim=3)
    nSpecies = maxval(species(1:nAtom))

    allocate(speciesL(orb%mOrb, orb%mOrb, 3, nSpecies))
    allocate(L(orb%mOrb, orb%mOrb, 3))
    do iSp = 1, nSpecies
      call localGetLOperatorsForSpecies(orb, iSp, L)
      speciesL(:, :, :, iSp) = aimag(L)
    end do

    allocate(tmpBlock(orb%mOrb, orb%mOrb))

    Lshell(:,:,:) = 0.0_dp
    do iAt = 1, nAtom
      iSp = species(iAt)
      nOrbSp = orb%nOrbSpecies(iSp)
      tmpBlock(:,:) = 0.0_dp
      ! Identity part
      tmpBlock(1:nOrbSp, 1:nOrbSp) = qBlockSkew(1:nOrbSp, 1:nOrbSp, iAt, 1)
      do iOrb = 1, orb%nShell(iSp)
        iOrbStart = orb%posShell(iOrb, iSp)
        iOrbEnd = orb%posShell(iOrb + 1, iSp) - 1
        do kk = 1, 3
          Lshell(kk, iOrb, iAt) = Lshell(kk, iOrb, iAt)&
              & - sum(speciesL(iOrbStart:iOrbEnd, iOrbStart:iOrbEnd, kk, iSp)&
              & * transpose(tmpBlock(iOrbStart:iOrbEnd, iOrbStart:iOrbEnd)))
        end do
      end do
    end do

  end subroutine getLDual


  !> Returns L_{x,y,z} in the tesseral spherical Harmonics basis used in DFTB+ for a given value of
  !> l
  subroutine localLOperators(ll, L)

    !> value of the orbital momentum to construct these matrices
    integer, intent(in) :: ll

    !> L_{x,y,z} operator
    complex(dp),intent(out) :: L(0:, 0:, :)

    ! L^+ operator
    complex(dp) :: Lplus(0:2*ll, 0:2*ll)

    ! L^- operator
    complex(dp) :: Lminus(0:2*ll, 0:2*ll)

    ! unitary matrix from spherical harmonics to tesseral harmonics
    complex(dp) :: uu(0:2*ll, 0:2*ll)

    ! magnetic quantum number
    integer :: mm

    integer :: iCart

    @:ASSERT(ll >= 0)
    @:ASSERT(all(shape(L) == [2 * ll + 1, 2 * ll + 1, 3]))

    L(:,:,:) = cmplx(0,0,dp)

    ! L_z in usual spherical harmonic basis
    do mm = -ll, ll
      L(ll + mm, ll + mm, 3) = real(mm, dp)
    end do

    ! L^+ in usual spherical harmonic basis
    Lplus(:,:) = cmplx(0,0,dp)
    do mm = -ll, ll - 1
      Lplus(ll + mm + 1, ll + mm) = sqrt(real(ll * (ll + 1) - mm * (mm + 1), dp))
    end do

    ! L^- in usual spherical harmonic basis
    Lminus(:,:) = cmplx(0,0,dp)
    do mm = -ll + 1, ll
      Lminus(ll + mm - 1, ll + mm) = sqrt(real(ll * (ll + 1) - mm * (mm - 1), dp))
    end do

    ! L_x
    L(:,:,1) = (Lplus + Lminus) / 2.0_dp
    ! L_y
    L(:,:,2) = (Lplus - Lminus) / (2.0_dp * imag)

    ! unitary transformation from Y_{lm} to \overline{Y}_{lm}
    uu(:,:) = cmplx(0,0,dp)
    do mm = 1, ll
      uu(ll + mm, ll  + mm) = sqrt(0.5_dp) * real(mod(mm + 1, 2)-mod(mm, 2), dp)
      uu(ll + mm, ll - mm) = sqrt(0.5_dp)
      uu(ll - mm, ll + mm) = -sqrt(0.5_dp) * imag * real(mod(mm, 2) - mod(mm + 1, 2), dp)
      uu(ll - mm, ll - mm) = -sqrt(0.5_dp) * imag
    end do
    ! m = 0 case
    uu(ll, ll) = 1.0_dp

    ! convert to tesseral form
    do iCart = 1, 3
      call makeSimilarityTrans(L(:,:,iCart), uu)
    end do

  end subroutine localLOperators


  !> Returns L_{x,y,z} in the tesseral spherical Harmonics basis for a give species.
  subroutine localGetLOperatorsForSpecies(orb, iSpecies, speciesL)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species to get the operators for
    integer, intent(in) :: iSpecies

    !> Species specific L operator
    complex(dp), intent(out) :: speciesL(:,:,:)

    integer :: iShell, ll, iOrbStart, iOrbEnd

    speciesL(:,:, :) = 0.0_dp
    do iShell = 1, orb%nShell(iSpecies)
      ll = orb%angShell(iShell, iSpecies)
      iOrbStart = orb%posShell(iShell, iSpecies)
      iOrbEnd = orb%posShell(iShell + 1, iSpecies) - 1
      call localLOperators(ll, speciesL(iOrbStart:iOrbEnd, iOrbStart:iOrbEnd, :))
    end do

  end subroutine localGetLOperatorsForSpecies


  !> Constructs a matrix to rotate tesseral spherical harmonic orbitals of angular momentum l around
  !> the z axis by phi radians
  pure subroutine zrot_onel(zmat,l, phi)

    !> resulting real unitary transformation matrix
    real(dp), intent(out) :: zmat(:,:)

    !> l value of angular momentum
    integer, intent(in)   :: l

    !> angle of rotation in radians
    real(dp), intent(in)  :: phi

    integer  :: m ! magnetic quantum number

    zmat(:,:) = 0.0_dp
    zmat(l+1,l+1) = 1.0_dp ! l_z = 0

    do m = 1, l
      zmat(m+l+1,m+l+1) = cos(m*phi)
      zmat(-m+l+1,-m+l+1) = cos(m*phi)
      zmat(m+l+1,-m+l+1) = -sin(m*phi)
      zmat(-m+l+1,m+l+1) = sin(m*phi)
    end do

  end subroutine zrot_onel


  !> Constructs a matrix to rotate tesseral spherical harmonic orbitals of angular momentum l around
  !> the z axis by phi radians
  pure subroutine zrot_manyl(zmat,l, phi)

    !> resulting real unitary transformation matrix
    real(dp), intent(out) :: zmat(:,:)

    !> l value of angular momentum
    integer, intent(in)   :: l(:)

    !> angle of rotation in radians
    real(dp), intent(in)  :: phi

    integer :: il, iStart, iEnd

    zmat(:,:) = 0.0_dp

    iStart = 1
    do il = 1, size(l)
      iEnd = iStart + 2*l(il)
      call zrot_onel(zmat(iStart:iEnd,iStart:iEnd),l(il), phi)
      iStart = iEnd + 1
    end do

  end subroutine zrot_manyl

end module dftbp_angmomentum
