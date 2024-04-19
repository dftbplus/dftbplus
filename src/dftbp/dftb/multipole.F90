!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Routines implementing the (One-center approximation) multipole expansion for the 2nd order DFTB.
module dftbp_dftb_multipole
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_globalenv, only : stdOut
  use dftbp_dftb_nonscc, only : TNonSccDiff
  use dftbp_dftb_shortgammafuncs, only : expGammaPrime, expGammaDoublePrime, expGammaTriplePrime,&
      & expGammaQuadruplePrime, expGammaQuintuplePrime
  use dftbp_dftb_slakocont, only : TSlakoCont
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_math_matrixops, only : adjointLowerTriangle
  use dftbp_dftb_periodic, only : TNeighbourList, getNrOfNeighbours
  implicit none
  private

  public :: TDftbMultiPoleInp, TDftbMultiPole, dftbMultiPole_init

  !> Input for the MultiPole module
  type TDftbMultiPoleInp

    !> Orbital information
    integer :: nOrb, nSpin, nSpecies
    type(TOrbitals), pointer :: orb

    !> Hubbard U values for atoms
    real(dp), allocatable :: hubbu(:)

    !> Species of atoms
    integer, allocatable :: species(:)

    !> One-center atomic intergral data
    real(dp), allocatable :: atomicDIntgrlScaling(:)
    real(dp), allocatable :: atomicQIntgrlScaling(:)
    real(dp), allocatable :: atomicOnsiteScaling(:)
    real(dp), allocatable :: atomicSXPxIntgrl(:)
    real(dp), allocatable :: atomicPxXDxxyyIntgrl(:)
    real(dp), allocatable :: atomicPxXDzzIntgrl(:)
    real(dp), allocatable :: atomicPyYDxxyyIntgrl(:)
    real(dp), allocatable :: atomicPzZDzzIntgrl(:)
    real(dp), allocatable :: atomicSXXSIntgrl(:)
    real(dp), allocatable :: atomicPxXXPxIntgrl(:)
    real(dp), allocatable :: atomicPyXXPyIntgrl(:)
    real(dp), allocatable :: atomicSXXDxxyyIntgrl(:)
    real(dp), allocatable :: atomicSXXDzzIntgrl(:)
    real(dp), allocatable :: atomicSYYDxxyyIntgrl(:)
    real(dp), allocatable :: atomicSZZDzzIntgrl(:)
    real(dp), allocatable :: atomicDxyXXDxyIntgrl(:)
    real(dp), allocatable :: atomicDyzXXDyzIntgrl(:)
    real(dp), allocatable :: atomicDxxyyXXDzzIntgrl(:)
    real(dp), allocatable :: atomicDzzXXDzzIntgrl(:)
    real(dp), allocatable :: atomicDxxyyYYDzzIntgrl(:)
    real(dp), allocatable :: atomicDzzZZDzzIntgrl(:)
    real(dp), allocatable :: atomicDxzXZDzzIntgrl(:)
    real(dp), allocatable :: atomicDyzYZDxxyyIntgrl(:)

  end type TDftbMultiPoleInp


  !> Internal status of TMultiPole.
  type TDftbMultiPole
    integer :: nAtoms, nSpecies, mShells, mShellsReal, nOrb, nSpin
    integer, allocatable :: nShells(:)

    !> Hubbard U values for atoms
    real(dp), allocatable :: hubbu(:)

    !> Species of atoms
    integer, allocatable :: species(:), nOrbSpecies(:)
    integer, allocatable :: onDQOCharges(:,:)

    !> coordinates of the atom
    real(dp), allocatable :: coords(:,:)

    real(dp), allocatable :: atomicDIntgrl(:,:,:,:)
    real(dp), allocatable :: atomicQIntgrl(:,:,:,:,:)

    real(dp), allocatable :: atomicOnsiteScaling(:)

    !> Mulliken charge per atom
    real(dp), allocatable :: deltaMAtom(:)
    !> Dipole charge per atom
    real(dp), allocatable :: deltaDAtom(:,:)
    !> Quadrupole charge per atom
    real(dp), allocatable :: deltaQAtom(:,:,:)

    !> evaluated for the E, Atom1, Atom2 at each geometry step
    real(dp), allocatable :: f10AB(:,:,:)
    real(dp), allocatable :: f20AB(:,:,:,:)
    real(dp), allocatable :: f30AB(:,:,:,:,:)
    real(dp), allocatable :: f40AB(:,:,:,:,:,:)
    !> evaluated for the gradient, Atom1, Atom2 at each geometry step
    real(dp), allocatable :: f50AB(:,:,:,:,:,:,:)

    !> add for the H and the E
    real(dp), allocatable :: pot10x1Atom(:)
    real(dp), allocatable :: pot20x2Atom(:)
    real(dp), allocatable :: pot10x0Atom(:,:)
    real(dp), allocatable :: pot11x1Atom(:,:)
    real(dp), allocatable :: pot21x2Atom(:,:)
    real(dp), allocatable :: pot20x0Atom(:,:,:)
    real(dp), allocatable :: pot21x1Atom(:,:,:)
    real(dp), allocatable :: pot22x2Atom(:,:,:)
    !> add for the gradient
    real(dp), allocatable :: pot30x0Atom(:,:,:,:)
    real(dp), allocatable :: pot31x1Atom(:,:,:,:)
    real(dp), allocatable :: pot32x2Atom(:,:,:,:)

    !> total energy
    real(dp) :: energyTT, energyMD, energyDD, energyMQ, energyDQ, energyQQ

  contains
    procedure :: updateCoords
    procedure :: updateDQPotentials
    procedure :: updateDeltaDQAtom
    procedure :: addMultiPoleHamiltonian
    procedure :: addMultiPoleEnergy
    procedure :: addMultiPoleGradients
    procedure :: addAtomicDipoleMoment
    procedure :: addAtomicQuadrupoleMoment
  end type TDftbMultiPole

contains


  !> Initializes instance.
  subroutine dftbMultiPole_init(this, inp)

    !> Instance.
    type(TDftbMultiPole), intent(out) :: this

    !> Input data.
    type(TDftbMultiPoleInp), intent(in) :: inp

    integer, parameter :: xyzLen = 3
    integer, parameter :: icx = 1, icy = 2, icz = 3
    integer, parameter :: ios = 1, iopy = 2, iopz = 3, iopx = 4
    integer, parameter :: iodxy = 5, iodyz = 6, iodzz = 7, iodxz = 8, iodxxyy = 9
    real(dp), parameter :: tolZero = 1.0e-15_dp
    integer :: nAtoms, maxNOrb, iAt1, iAt2, nOrb1, nOrb2, ii, jj, iSp1, iSp2
    integer :: mu, nu, mm, nn
    real(dp) tmpIntgrl, tmpTrace

    this%nAtoms = size(inp%orb%nOrbAtom)
    this%nSpecies = inp%nSpecies
    allocate(this%nOrbSpecies(this%nSpecies))
    this%nOrbSpecies(:) = inp%orb%nOrbSpecies(:)
    this%mShells = 1
    this%mShellsReal = inp%orb%mShell
    this%nOrb = inp%nOrb
    maxNOrb = maxval(inp%orb%nOrbSpecies(:))
    this%nSpin = inp%nSpin
    allocate(this%nShells(this%nSpecies))
    this%nShells(:) = 1

    allocate(this%hubbu(size(inp%hubbu)))
    this%hubbu(:) = inp%hubbu

    allocate(this%species(this%nAtoms))
    this%species(:) = inp%species

    allocate(this%onDQOCharges(3, this%nSpecies))
    this%onDQOCharges(1,:) = 1
    this%onDQOCharges(2,:) = 1
    this%onDQOCharges(3,:) = 0

    allocate(this%coords(xyzLen, this%nAtoms))
    this%coords(:,:) = 0.0_dp

    allocate(this%atomicDIntgrl(xyzLen, maxNOrb, maxNOrb, this%nSpecies))
    allocate(this%atomicQIntgrl(xyzLen, xyzLen, maxNOrb, maxNOrb, this%nSpecies))
    this%atomicDIntgrl(:,:,:,:) = 0.0_dp
    this%atomicQIntgrl(:,:,:,:,:) = 0.0_dp

    allocate(this%atomicOnsiteScaling(this%nSpecies))
    this%atomicOnsiteScaling(:) = 1.0_dp


    do iSp1 = 1, this%nSpecies
      nOrb1 = this%nOrbSpecies(iSp1)

      !Dipole
      if ( nOrb1 >= 4 ) then
        !assign <S|X|Px>
        tmpIntgrl = inp%atomicSXPxIntgrl(iSp1)
        this%atomicDIntgrl(icx,ios,iopx,iSp1) = tmpIntgrl
        this%atomicDIntgrl(icy,ios,iopy,iSp1) = tmpIntgrl
        this%atomicDIntgrl(icz,ios,iopz,iSp1) = tmpIntgrl
      end if

      if ( nOrb1 >= 9 ) then
        !assign <Px|X|Dxx-yy>
        tmpIntgrl = inp%atomicPxXDxxyyIntgrl(iSp1)
        this%atomicDIntgrl(icx,iopx,iodxxyy,iSp1) = tmpIntgrl
        this%atomicDIntgrl(icx,iopy,iodxy,iSp1) = tmpIntgrl
        this%atomicDIntgrl(icx,iopz,iodxz,iSp1) = tmpIntgrl
        this%atomicDIntgrl(icy,iopx,iodxy,iSp1) = tmpIntgrl
        this%atomicDIntgrl(icy,iopz,iodyz,iSp1) = tmpIntgrl
        this%atomicDIntgrl(icz,iopx,iodxz,iSp1) = tmpIntgrl
        this%atomicDIntgrl(icz,iopy,iodyz,iSp1) = tmpIntgrl

        !assign <Px|X|Dzz>
        tmpIntgrl = inp%atomicPxXDzzIntgrl(iSp1)
        this%atomicDIntgrl(icx,iopx,iodzz,iSp1) = tmpIntgrl
        this%atomicDIntgrl(icy,iopy,iodzz,iSp1) = tmpIntgrl

        !assign <Py|Y|Dxx-yy>
        tmpIntgrl = inp%atomicPyYDxxyyIntgrl(iSp1)
        this%atomicDIntgrl(icy,iopy,iodxxyy,iSp1) = tmpIntgrl

        !assign <Pz|Z|Dzz>
        tmpIntgrl = inp%atomicPzZDzzIntgrl(iSp1)
        this%atomicDIntgrl(icz,iopz,iodzz,iSp1) = tmpIntgrl
      end if

      !Quadrupole
      if ( nOrb1 >= 1 ) then
        !assign <S|XX|S>
        tmpIntgrl = inp%atomicSXXSIntgrl(iSp1)
        this%atomicQIntgrl(icx,icx,ios,ios,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icy,icy,ios,ios,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icz,icz,ios,ios,iSp1) = tmpIntgrl
      end if

      if ( nOrb1 >= 4 ) then
        !assign <Px|XX|Px>
        tmpIntgrl = inp%atomicPxXXPxIntgrl(iSp1)
        this%atomicQIntgrl(icx,icx,iopx,iopx,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icy,icy,iopy,iopy,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icz,icz,iopz,iopz,iSp1) = tmpIntgrl

        !assign <Py|XX|Py>
        tmpIntgrl = inp%atomicPyXXPyIntgrl(iSp1)
        this%atomicQIntgrl(icx,icx,iopy,iopy,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icx,icx,iopz,iopz,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icy,icy,iopx,iopx,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icy,icy,iopz,iopz,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icz,icz,iopx,iopx,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icz,icz,iopy,iopy,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icx,icy,iopx,iopy,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icx,icz,iopx,iopz,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icy,icz,iopy,iopz,iSp1) = tmpIntgrl
      end if

      if ( nOrb1 >= 9 ) then
        !assign <S|XX|Dxx-yy>
        tmpIntgrl = inp%atomicSXXDxxyyIntgrl(iSp1)
        this%atomicQIntgrl(icx,icx,ios,iodxxyy,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icx,icy,ios,iodxy,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icx,icz,ios,iodxz,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icy,icz,ios,iodyz,iSp1) = tmpIntgrl

        !assign <S|XX|Dzz>
        tmpIntgrl = inp%atomicSXXDzzIntgrl(iSp1)
        this%atomicQIntgrl(icx,icx,ios,iodzz,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icy,icy,ios,iodzz,iSp1) = tmpIntgrl

        !assign <S|YY|Dxx-yy>
        tmpIntgrl = inp%atomicSYYDxxyyIntgrl(iSp1)
        this%atomicQIntgrl(icy,icy,ios,iodxxyy,iSp1) = tmpIntgrl

        !assign <S|ZZ|Dzz>
        tmpIntgrl = inp%atomicSZZDzzIntgrl(iSp1)
        this%atomicQIntgrl(icz,icz,ios,iodzz,iSp1) = tmpIntgrl

        !assign <Dxy|XX|Dxy>
        tmpIntgrl = inp%atomicDxyXXDxyIntgrl(iSp1)
        this%atomicQIntgrl(icx,icx,iodxy,iodxy,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icx,icx,iodxz,iodxz,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icx,icx,iodxxyy,iodxxyy,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icy,icy,iodxy,iodxy,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icy,icy,iodyz,iodyz,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icy,icy,iodxxyy,iodxxyy,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icz,icz,iodxz,iodxz,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icz,icz,iodyz,iodyz,iSp1) = tmpIntgrl

        !assign <Dyz|XX|Dyz>
        tmpIntgrl = inp%atomicDyzXXDyzIntgrl(iSp1)
        this%atomicQIntgrl(icx,icx,iodyz,iodyz,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icy,icy,iodxz,iodxz,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icz,icz,iodxy,iodxy,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icz,icz,iodxxyy,iodxxyy,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icx,icy,iodxz,iodyz,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icx,icz,iodxy,iodyz,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icx,icz,iodxz,iodxxyy,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icy,icz,iodxy,iodxz,iSp1) = tmpIntgrl

        !assign <Dxx-yy|XX|Dzz>
        tmpIntgrl = inp%atomicDxxyyXXDzzIntgrl(iSp1)
        this%atomicQIntgrl(icx,icx,iodxxyy,iodzz,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icx,icy,iodxy,iodzz,iSp1) = tmpIntgrl

        !assign <Dzz|XX|Dzz>
        tmpIntgrl = inp%atomicDzzXXDzzIntgrl(iSp1)
        this%atomicQIntgrl(icx,icx,iodzz,iodzz,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icy,icy,iodzz,iodzz,iSp1) = tmpIntgrl

        !assign <Dxx-yy|YY|Dzz>
        tmpIntgrl = inp%atomicDxxyyYYDzzIntgrl(iSp1)
        this%atomicQIntgrl(icy,icy,iodxxyy,iodzz,iSp1) = tmpIntgrl

        !assign <Dzz|ZZ|Dzz>
        tmpIntgrl = inp%atomicDzzZZDzzIntgrl(iSp1)
        this%atomicQIntgrl(icz,icz,iodzz,iodzz,iSp1) = tmpIntgrl

        !assign <Dxz|XZ|Dzz>
        tmpIntgrl = inp%atomicDxzXZDzzIntgrl(iSp1)
        this%atomicQIntgrl(icx,icz,iodxz,iodzz,iSp1) = tmpIntgrl
        this%atomicQIntgrl(icy,icz,iodyz,iodzz,iSp1) = tmpIntgrl

        !assign <Dyz|YZ|Dxx-yy>
        tmpIntgrl = inp%atomicDyzYZDxxyyIntgrl(iSp1)
        this%atomicQIntgrl(icy,icz,iodyz,iodxxyy,iSp1) = tmpIntgrl
      end if
    end do

    do iSp1 = 1, this%nSpecies
      nOrb1 = this%nOrbSpecies(iSp1)
      do mm = 1, nOrb1
        do nn = 1, nOrb1
          do ii = 1, xyzLen
            if (abs(this%atomicDIntgrl(ii,mm,nn,iSp1)) > tolZero ) then
              this%atomicDIntgrl(ii,nn,mm,iSp1) = this%atomicDIntgrl(ii,mm,nn,iSp1)
            end if
            do jj = 1, xyzLen
              if (abs(this%atomicQIntgrl(ii,jj,mm,nn,iSp1)) > tolZero ) then
                this%atomicQIntgrl(ii,jj,nn,mm,iSp1) = this%atomicQIntgrl(ii,jj,mm,nn,iSp1)
                this%atomicQIntgrl(jj,ii,mm,nn,iSp1) = this%atomicQIntgrl(ii,jj,mm,nn,iSp1)
                this%atomicQIntgrl(jj,ii,nn,mm,iSp1) = this%atomicQIntgrl(ii,jj,mm,nn,iSp1)
              end if
            end do
          end do
        end do
      end do
    end do

   !> scaling atomic Intgrls
    do iSp1 = 1, this%nSpecies
      this%atomicDIntgrl(:,:,:,iSp1) = inp%atomicDIntgrlScaling(iSp1)&
          & * this%atomicDIntgrl(:,:,:,iSp1)
      this%atomicQIntgrl(:,:,:,:,iSp1) = inp%atomicQIntgrlScaling(iSp1)&
          & * this%atomicQIntgrl(:,:,:,:,iSp1)
    end do

   !> Remove Trace of atomicQIntgrl
    do iSp1 = 1, this%nSpecies
      nOrb1 = this%nOrbSpecies(iSp1)
      do mm = 1, nOrb1
        do nn = 1, nOrb1
          tmpTrace = (this%atomicQIntgrl(1,1,mm,nn,iSp1) + this%atomicQIntgrl(2,2,mm,nn,iSp1)&
              & + this%atomicQIntgrl(3,3,mm,nn,iSp1)) / 3.0_dp
          do ii = 1, xyzLen
            this%atomicQIntgrl(ii,ii,mm,nn,iSp1) = this%atomicQIntgrl(ii,ii,mm,nn,iSp1) - tmpTrace
          end do
        end do
      end do
    end do

    do iSp1 = 1, this%nSpecies
      if (maxval(abs(this%atomicDIntgrl(:,:,:,iSp1))) < tolZero ) then
        this%onDQOCharges(1,iSp1) = 0
      end if
      if (maxval(abs(this%atomicQIntgrl(:,:,:,:,iSp1))) < tolZero ) then
        this%onDQOCharges(2,iSp1) = 0
      end if
    end do

    this%atomicOnsiteScaling(:) = inp%atomicOnsiteScaling(:)


    allocate(this%deltaMAtom(this%nAtoms))
    allocate(this%deltaDAtom(xyzLen, this%nAtoms))
    allocate(this%deltaQAtom(xyzLen, xyzLen, this%nAtoms))
    allocate(this%f10AB(xyzLen, this%nAtoms, this%nAtoms))
    allocate(this%f20AB(xyzLen, xyzLen, this%nAtoms, this%nAtoms))
    allocate(this%f30AB(xyzLen, xyzLen, xyzLen, this%nAtoms, this%nAtoms))
    allocate(this%f40AB(xyzLen, xyzLen, xyzLen, xyzLen, this%nAtoms, this%nAtoms))
    allocate(this%f50AB(xyzLen, xyzLen, xyzLen, xyzLen, xyzLen, this%nAtoms, this%nAtoms))

    allocate(this%pot10x1Atom(this%nAtoms))
    allocate(this%pot20x2Atom(this%nAtoms))
    allocate(this%pot10x0Atom(xyzLen, this%nAtoms))
    allocate(this%pot11x1Atom(xyzLen, this%nAtoms))
    allocate(this%pot21x2Atom(xyzLen, this%nAtoms))
    allocate(this%pot20x0Atom(xyzLen, xyzLen, this%nAtoms))
    allocate(this%pot21x1Atom(xyzLen, xyzLen, this%nAtoms))
    allocate(this%pot22x2Atom(xyzLen, xyzLen, this%nAtoms))
    allocate(this%pot30x0Atom(xyzLen, xyzLen, xyzLen, this%nAtoms))
    allocate(this%pot31x1Atom(xyzLen, xyzLen, xyzLen, this%nAtoms))
    allocate(this%pot32x2Atom(xyzLen, xyzLen, xyzLen, this%nAtoms))

    this%deltaMAtom(:) = 0.0_dp
    this%deltaDAtom(:,:) = 0.0_dp
    this%deltaQAtom(:,:,:) = 0.0_dp
    this%f10AB(:,:,:) = 0.0_dp
    this%f20AB(:,:,:,:) = 0.0_dp
    this%f30AB(:,:,:,:,:) = 0.0_dp
    this%f40AB(:,:,:,:,:,:) = 0.0_dp
    this%f50AB(:,:,:,:,:,:,:) = 0.0_dp

    this%pot10x1Atom(:) = 0.0_dp
    this%pot20x2Atom(:) = 0.0_dp
    this%pot10x0Atom(:,:) = 0.0_dp
    this%pot11x1Atom(:,:) = 0.0_dp
    this%pot21x2Atom(:,:) = 0.0_dp
    this%pot20x0Atom(:,:,:) = 0.0_dp
    this%pot21x1Atom(:,:,:) = 0.0_dp
    this%pot22x2Atom(:,:,:) = 0.0_dp
    this%pot30x0Atom(:,:,:,:) = 0.0_dp
    this%pot31x1Atom(:,:,:,:) = 0.0_dp
    this%pot32x2Atom(:,:,:,:) = 0.0_dp

  end subroutine dftbMultiPole_init


  !> Updates data structures if there are changed coordinates for the instance.
  subroutine updateCoords(this, coords)

    !> class instance
    class(TDftbMultiPole), intent(inout) :: this

    !> list of atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    integer, parameter :: xyzLen = 3

    integer :: nAtoms, iAt1, iAt2, iSp1, iSp2, nOrb1, nOrb2, ii, jj, ll, mm, nn, kk
    real(dp) :: rab, u1, u2, coeffTerm1, coeffTerm2, coeffTerm3
    real(dp) :: gammaPrime, gammaDoublePrime, gammaTriplePrime, gammaQuadruplePrime,&
        & gammaQuintuplePrime

    real(dp) :: vRabx3(xyzLen), workVx3(xyzLen)
    real(dp) :: mI3x3(xyzLen, xyzLen), mRRab3x3(xyzLen, xyzLen), workM3x3(xyzLen, xyzLen),&
        & workM3x3x3(xyzLen, xyzLen, xyzLen)
    real(dp) :: mRoIab3x3x3(xyzLen, xyzLen, xyzLen), mIotRab3x3x3(xyzLen, xyzLen, xyzLen),&
        & mIoRab3x3x3(xyzLen, xyzLen, xyzLen)
    real(dp) :: mRIotoRab3x3x3(xyzLen, xyzLen, xyzLen), mRRRab3x3x3(xyzLen, xyzLen, xyzLen)

    real(dp) :: mI3x3otI3x3(xyzLen, xyzLen, xyzLen, xyzLen), mI3x3ohI3x3(xyzLen, xyzLen, xyzLen,&
        & xyzLen)
    real(dp) :: mI3x3oI3x3(xyzLen, xyzLen, xyzLen, xyzLen), mI3x3otohoI3x3(xyzLen, xyzLen, xyzLen,&
        & xyzLen)
    real(dp) :: workM3x3x3x3(xyzLen, xyzLen, xyzLen, xyzLen), f40Term1M3x3x3x3(xyzLen, xyzLen,&
        & xyzLen, xyzLen)
    real(dp) :: f40Term3M3x3x3x3(xyzLen, xyzLen, xyzLen, xyzLen), f40Term2M3x3x3x3(xyzLen, xyzLen,&
        & xyzLen, xyzLen)
    real(dp) :: workM3x3x3x3x3(xyzLen, xyzLen, xyzLen, xyzLen, xyzLen),&
        & f50Term1M3x3x3x3x3(xyzLen, xyzLen, xyzLen, xyzLen, xyzLen)
    real(dp) :: f50Term3M3x3x3x3x3(xyzLen, xyzLen, xyzLen, xyzLen, xyzLen),&
        & f50Term2M3x3x3x3x3(xyzLen, xyzLen, xyzLen, xyzLen, xyzLen)

    this%coords(:,:) = coords
    this%f10AB(:,:,:) = 0.0_dp
    this%f20AB(:,:,:,:) = 0.0_dp
    this%f30AB(:,:,:,:,:) = 0.0_dp
    this%f40AB(:,:,:,:,:,:) = 0.0_dp
    this%f50AB(:,:,:,:,:,:,:) = 0.0_dp

    mI3x3(:,:) = 0.0_dp
    do ii = 1, xyzLen
      mI3x3(ii, ii) = 1.0_dp
    end do
    call makeA3x3oB3x3(mI3x3oI3x3, mI3x3, mI3x3)
    call makeA3x3ohB3x3(mI3x3ohI3x3, mI3x3, mI3x3)
    call makeA3x3otB3x3(mI3x3otI3x3, mI3x3, mI3x3)
    mI3x3otohoI3x3(:,:,:,:) = mI3x3otI3x3(:,:,:,:) + mI3x3ohI3x3(:,:,:,:) + mI3x3oI3x3(:,:,:,:)

    nAtoms = this%nAtoms

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(iAt1, iSp1, u1, iAt2, iSp2, u2, rab, vRabx3, workVx3) &
    !$OMP PRIVATE(gammaPrime, gammaDoublePrime, gammaTriplePrime, gammaQuadruplePrime) &
    !$OMP PRIVATE(gammaQuintuplePrime, coeffTerm1, coeffTerm2, coeffTerm3) &
    !$OMP PRIVATE(workM3x3, workM3x3x3, mRRab3x3) &
    !$OMP PRIVATE(mRoIab3x3x3, mIotRab3x3x3, mIoRab3x3x3, mRIotoRab3x3x3, mRRRab3x3x3) &
    !$OMP PRIVATE(workM3x3x3x3, f40Term1M3x3x3x3, f40Term2M3x3x3x3, f40Term3M3x3x3x3) &
    !$OMP PRIVATE(workM3x3x3x3x3, f50Term1M3x3x3x3x3, f50Term2M3x3x3x3x3, f50Term3M3x3x3x3x3)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtoms
      iSp1 = this%species(iAt1)
      u1 = this%hubbu(iSp1)
      do iAt2 = 1, iAt1 - 1
        iSp2 = this%species(iAt2)
        u2 = this%hubbu(iSp2)
        vRabx3(:) = this%coords(:,iAt1) - this%coords(:,iAt2)
        rab = norm2(vRabx3(:))

        !> f10AB
        gammaPrime = -1.0_dp / rab**2 - expGammaPrime(rab, u1, u2)
        coeffTerm1 = gammaPrime / rab

        workVx3(:) = coeffTerm1 * vRabx3(:)
        this%f10AB(:,iAt2,iAt1) = workVx3(:)
        this%f10AB(:,iAt1,iAt2) = -workVx3(:)

        !> f20AB
        gammaDoublePrime = 2.0_dp / rab**3 - expGammaDoublePrime(rab, u1, u2)
        coeffTerm1 = gammaPrime / rab
        coeffTerm2 = gammaDoublePrime / rab**2 - gammaPrime / rab**3
        coeffTerm1 = coeffTerm1 / 2.0_dp
        coeffTerm2 = coeffTerm2 / 2.0_dp

        call makeA3oB3(mRRab3x3, vRabx3, vRabx3)
        workM3x3(:,:) = coeffTerm1 * mI3x3(:,:) &
                     &+ coeffTerm2 * mRRab3x3(:,:)
        this%f20AB(:,:,iAt2,iAt1) = workM3x3(:,:)
        this%f20AB(:,:,iAt1,iAt2) = workM3x3(:,:)

        !> f30AB
        gammaTriplePrime = -6.0_dp / rab**4 - expGammaTriplePrime(rab, u1, u2)
        coeffTerm1 = gammaDoublePrime / rab**2 - gammaPrime / rab**3
        coeffTerm2 = gammaTriplePrime / rab**3 - 3.0_dp * gammaDoublePrime / rab**4&
            & + 3.0_dp * gammaPrime / rab**5
        coeffTerm1 = coeffTerm1 / 6.0_dp
        coeffTerm2 = coeffTerm2 / 6.0_dp

        do ll = 1, xyzLen
          do jj = 1, xyzLen
            do ii = 1, xyzLen
              mRoIab3x3x3(ii,jj,ll) = vRabx3(ii) * mI3x3(jj,ll)
              mIotRab3x3x3(ii,jj,ll) = mI3x3(ii,ll) * vRabx3(jj)
              mIoRab3x3x3(ii,jj,ll) = mI3x3(ii,jj) * vRabx3(ll)
              mRRRab3x3x3(ii,jj,ll) = mRRab3x3(ii,jj) * vRabx3(ll)
            end do
          end do
        end do
        mRIotoRab3x3x3(:,:,:) = mRoIab3x3x3(:,:,:) + mIotRab3x3x3(:,:,:) + mIoRab3x3x3(:,:,:)
        workM3x3x3(:,:,:) = coeffTerm1 * mRIotoRab3x3x3(:,:,:) &
                         &+ coeffTerm2 * mRRRab3x3x3(:,:,:)
        this%f30AB(:,:,:,iAt2,iAt1) = workM3x3x3(:,:,:)
        this%f30AB(:,:,:,iAt1,iAt2) = -workM3x3x3(:,:,:)

        !> for f40AB and f50AB
        if ((this%onDQOCharges(1, iSp1) == 0 .and. this%onDQOCharges(2, iSp1) == 0 ) &
           & .or. (this%onDQOCharges(1, iSp2) == 0 .and. this%onDQOCharges(2, iSp2) == 0)) then
          cycle
        end if

        !> f40AB
        gammaQuadruplePrime = 24.0_dp / rab**5 - expGammaQuadruplePrime(rab, u1, u2)
        coeffTerm1 = gammaDoublePrime / rab**2 - gammaPrime / rab**3
        coeffTerm2 = gammaTriplePrime / rab**3 &
                  &- 3.0_dp * gammaDoublePrime / rab**4 + 3.0_dp * gammaPrime / rab**5
        coeffTerm3 = gammaQuadruplePrime / rab**4 - 6.0_dp * gammaTriplePrime / rab**5 &
                  &+ 15.0_dp * gammaDoublePrime / rab**6 - 15.0_dp * gammaPrime / rab**7
        coeffTerm1 = coeffTerm1 / 24.0_dp
        coeffTerm2 = coeffTerm2 / 24.0_dp
        coeffTerm3 = coeffTerm3 / 24.0_dp


        f40Term1M3x3x3x3(:,:,:,:) = mI3x3otohoI3x3(:,:,:,:)
        f40Term2M3x3x3x3(:,:,:,:) = 0.0_dp
        f40Term3M3x3x3x3(:,:,:,:) = 0.0_dp
        do mm = 1, xyzLen
          do ll = 1, xyzLen
            do jj = 1, xyzLen
              do ii = 1, xyzLen
                f40Term2M3x3x3x3(ii,jj,ll,mm) = f40Term2M3x3x3x3(ii,jj,ll,mm) &
                                                &+ mRRab3x3(ii,jj) * mI3x3(ll,mm) &
                                                &+ mRRab3x3(ii,ll) * mI3x3(jj,mm) &
                                                &+ mRRab3x3(ii,mm) * mI3x3(jj,ll) &
                                                &+ mI3x3(ii,jj) * mRRab3x3(ll,mm) &
                                                &+ mI3x3(ii,ll) * mRRab3x3(jj,mm) &
                                                &+ mI3x3(ii,mm) * mRRab3x3(jj,ll)
                f40Term3M3x3x3x3(ii,jj,ll,mm) = f40Term3M3x3x3x3(ii,jj,ll,mm) &
                                                &+ mRRab3x3(ii,jj) * mRRab3x3(ll,mm)
              end do
            end do
          end do
        end do
        workM3x3x3x3(:,:,:,:) = coeffTerm1 * f40Term1M3x3x3x3(:,:,:,:) &
                             &+ coeffTerm2 * f40Term2M3x3x3x3(:,:,:,:) &
                             &+ coeffTerm3 * f40Term3M3x3x3x3(:,:,:,:)
        this%f40AB(:,:,:,:,iAt2,iAt1) = workM3x3x3x3(:,:,:,:)
        this%f40AB(:,:,:,:,iAt1,iAt2) = workM3x3x3x3(:,:,:,:)

        !> f50AB
        gammaQuintuplePrime = -120.0_dp / rab**6 - expGammaQuintuplePrime(rab, u1, u2)
        coeffTerm1 = gammaTriplePrime / rab**3 &
                  &- 3.0_dp * gammaDoublePrime / rab**4 + 3.0_dp * gammaPrime / rab**5
        coeffTerm2 = gammaQuadruplePrime / rab**4 - 6.0_dp * gammaTriplePrime / rab**5 &
                  &+ 15.0_dp * gammaDoublePrime / rab**6 - 15.0_dp * gammaPrime / rab**7
        coeffTerm3 = gammaQuintuplePrime / rab**5 &
                   &- 10.0_dp * gammaQuadruplePrime / rab**6 + 45.0_dp * gammaTriplePrime / rab**7 &
                   &- 105.0_dp * gammaDoublePrime / rab**8 + 105.0_dp * gammaPrime / rab**9
        coeffTerm1 = coeffTerm1 / 120.0_dp
        coeffTerm2 = coeffTerm2 / 120.0_dp
        coeffTerm3 = coeffTerm3 / 120.0_dp

        f50Term1M3x3x3x3x3(:,:,:,:,:) = 0.0_dp
        f50Term2M3x3x3x3x3(:,:,:,:,:) = 0.0_dp
        f50Term3M3x3x3x3x3(:,:,:,:,:) = 0.0_dp
        do nn = 1, xyzLen
          do mm = 1, xyzLen
            do ll = 1, xyzLen
              do jj = 1, xyzLen
                do ii = 1, xyzLen
                  f50Term1M3x3x3x3x3(ii,jj,ll,mm,nn) = f50Term1M3x3x3x3x3(ii,jj,ll,mm,nn) &
                                                  &+ mRoIab3x3x3(ii,jj,ll) * mI3x3(mm,nn) &
                                                  &+ mRoIab3x3x3(ii,jj,mm) * mI3x3(ll,nn) &
                                                  &+ mRoIab3x3x3(ii,jj,nn) * mI3x3(ll,mm) &
                                                  &+ mI3x3(ii,jj) * mRIotoRab3x3x3(ll,mm,nn) &
                                                  &+ mI3x3(ii,ll) * mRIotoRab3x3x3(jj,mm,nn) &
                                                  &+ mI3x3(ii,mm) * mRIotoRab3x3x3(jj,ll,nn) &
                                                  &+ mI3x3(ii,nn) * mRIotoRab3x3x3(jj,ll,mm)
                  f50Term2M3x3x3x3x3(ii,jj,ll,mm,nn) = f50Term2M3x3x3x3x3(ii,jj,ll,mm,nn) &
                                                  &+ mRRRab3x3x3(ii,jj,ll) * mI3x3(mm,nn) &
                                                  &+ mRoIab3x3x3(ii,jj,nn) * mRRab3x3(ll,mm) &
                                                  &+ mI3x3(ii,jj) * mRRRab3x3x3(ll,mm,nn) &
                                                  &+ mI3x3(ii,ll) * mRRRab3x3x3(jj,mm,nn) &
                                                  &+ mI3x3(ii,mm) * mRRRab3x3x3(jj,ll,nn) &
                                                  &+ mI3x3(ii,nn) * mRRRab3x3x3(jj,ll,mm) &
                                                  &+ mRRab3x3(ii,jj) * mIoRab3x3x3(ll,mm,nn) &
                                                  &+ mRRab3x3(ii,ll) * mIoRab3x3x3(jj,mm,nn) &
                                                  &+ mRRab3x3(ii,mm) * mIoRab3x3x3(jj,ll,nn) &
                                                  &+ mRRab3x3(ii,jj) * mIotRab3x3x3(ll,mm,nn)
                  f50Term3M3x3x3x3x3(ii,jj,ll,mm,nn) = f50Term3M3x3x3x3x3(ii,jj,ll,mm,nn) &
                                                  &+ mRRRab3x3x3(ii,jj,ll) * mRRab3x3(mm,nn)
                end do
              end do
            end do
          end do
        end do
        workM3x3x3x3x3(:,:,:,:,:) = coeffTerm1 * f50Term1M3x3x3x3x3(:,:,:,:,:) &
                                 &+ coeffTerm2 * f50Term2M3x3x3x3x3(:,:,:,:,:) &
                                 &+ coeffTerm3 * f50Term3M3x3x3x3x3(:,:,:,:,:)

        this%f50AB(:,:,:,:,:,iAt2,iAt1) = workM3x3x3x3x3(:,:,:,:,:)
        this%f50AB(:,:,:,:,:,iAt1,iAt2) = -workM3x3x3x3x3(:,:,:,:,:)
      end do

      coeffTerm1 = -this%atomicOnsiteScaling(iSp1) * (3.2_dp * u1)**3 / 96.0_dp
      this%f20AB(:,:,iAt1,iAt1) = coeffTerm1 * mI3x3(:,:)

      coeffTerm1 = this%atomicOnsiteScaling(iSp1) * (3.2_dp * u1)**5 / 5760.0_dp
      this%f40AB(:,:,:,:,iAt1,iAt1) = coeffTerm1 * mI3x3otohoI3x3(:,:,:,:)

    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine updateCoords


  subroutine updateDeltaDQAtom(this, deltaRho, iSquare, orb, overlap)

    !> class instance
    class(TDftbMultiPole), intent(inout) :: this

    !> Square (unpacked) density matrix
    real(dp), dimension(:,:), target, intent(in) :: deltaRho

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtoms)
    integer, dimension(:), intent(in) :: iSquare

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    real(dp), allocatable :: tmpOvr(:,:), tmpDRho(:,:)
    integer, parameter :: xyzLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    integer :: matrixSize, nAtoms, maxNOrb
    integer :: iAt1, iAt2, iSp1, iSp2, nOrb1, nOrb2
    integer :: ii, jj, mu, nu, mm, nn, iStartOrb, iEndOrb

    real(dp) :: tmpVx3(xyzLen), tmpM3x3(xyzLen, xyzLen)
    real(dp) :: dPSSdP, tmpTrace

    matrixSize = size(overlap, dim = 1)
    nAtoms = this%nAtoms
    maxNOrb = maxval(orb%nOrbSpecies(:))

    allocate(tmpOvr(matrixSize, matrixSize))
    allocate(tmpDRho(matrixSize, matrixSize))

    tmpOvr(:,:) = overlap
    call adjointLowerTriangle(tmpOvr)

    tmpDRho(:,:) = deltaRho
    call symmetrizeSquareMatrix(tmpDRho)

    this%deltaDAtom(:,:) = 0.0_dp
    this%deltaQAtom(:,:,:) = 0.0_dp

    tmpVx3(:) = 0.0_dp
    tmpM3x3(:,:) = 0.0_dp

    do iAt1 = 1, nAtoms
      iSp1 = this%species(iAt1)
      if ((this%onDQOCharges(1, iSp1) == 0) .and. (this%onDQOCharges(2, iSp1) == 0)) then
        cycle
      end if

      iStartOrb  = iSquare(iAt1)
      iEndOrb = iSquare(iAt1 + 1) - 1
      tmpVx3(:) = 0.0_dp
      tmpM3x3(:,:) = 0.0_dp
      do mu = iStartOrb, iEndOrb
        mm = mu - iStartOrb + 1
        do nu = iStartOrb, iEndOrb
          nn = nu - iStartOrb + 1
          dPSSdP = sum(tmpDRho(mu,:) * tmpOvr(:,nu)) + sum(tmpOvr(nu,:) * tmpDRho(:,mu))
          tmpVx3(:) = tmpVx3(:) + this%atomicDIntgrl(:,mm,nn,iSp1) * dPSSdP
          tmpM3x3(:,:) = tmpM3x3(:,:) + this%atomicQIntgrl(:,:,mm,nn,iSp1) * dPSSdP
        end do
      end do
      this%deltaDAtom(:,iAt1) = this%deltaDAtom(:,iAt1) + 0.5_dp * tmpVx3(:)
      this%deltaQAtom(:,:,iAt1) = this%deltaQAtom(:,:,iAt1) + 0.5_dp * tmpM3x3(:,:)
    end do

    deallocate(tmpOvr, tmpDRho)

  end subroutine updateDeltaDQAtom


  subroutine updateDQPotentials(this, deltaMAtom)

    !> class instance
    class(TDftbMultiPole), intent(inout) :: this

    !> Delta charge per atom
    real(dp), intent(in) :: deltaMAtom(:)

    integer, parameter :: xyzLen = 3
    integer :: nAtoms, iAt1, iAt2, ii, jj, iSp1, iSp2
    real(dp) :: vRabx3(xyzLen), tmpVx3(xyzLen)

    this%deltaMAtom(:) = deltaMAtom
    nAtoms = this%nAtoms

    !> Calculate pot10x0, and pot20x0
    this%pot10x0Atom(:,:) = 0.0_dp
    this%pot20x0Atom(:,:,:) = 0.0_dp

    !$OMP PARALLEL DO PRIVATE(iAt1, iAt2) DEFAULT(SHARED) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtoms
      do iAt2 = 1, nAtoms
        this%pot10x0Atom(:,iAt1) = this%pot10x0Atom(:,iAt1)&
            & + this%f10AB(:,iAt2,iAt1) * this%deltaMAtom(iAt2)
        this%pot20x0Atom(:,:,iAt1) = this%pot20x0Atom(:,:,iAt1)&
            & + this%f20AB(:,:,iAt2,iAt1) * this%deltaMAtom(iAt2)
      end do
    end do
    !$OMP END PARALLEL DO

    !> Calculate pot10x1, pot20x2, pot11x1, pot21x2, and pot21x1
    this%pot10x1Atom(:) = 0.0_dp
    this%pot20x2Atom(:) = 0.0_dp
    this%pot11x1Atom(:,:) = 0.0_dp
    this%pot21x2Atom(:,:) = 0.0_dp
    this%pot21x1Atom(:,:,:) = 0.0_dp

    !$OMP PARALLEL DO PRIVATE(iAt1, iAt2, iSp2, ii, jj) DEFAULT(SHARED) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtoms
      do iAt2 = 1, nAtoms
        iSp2 = this%species(iAt2)
        if ((this%onDQOCharges(1, iSp2) == 0) .and. (this%onDQOCharges(2, iSp2) == 0)) then
          cycle
        end if
        this%pot10x1Atom(iAt1) = this%pot10x1Atom(iAt1)&
            & + sum(this%f10AB(:,iAt2,iAt1) * this%deltaDAtom(:,iAt2))
        this%pot20x2Atom(iAt1) = this%pot20x2Atom(iAt1)&
            & + sum(this%f20AB(:,:,iAt2,iAt1) * this%deltaQAtom(:,:,iAt2))
        do ii = 1, xyzLen
          this%pot11x1Atom(ii,iAt1) = this%pot11x1Atom(ii,iAt1)&
              & - 2.0_dp * sum(this%f20AB(:,ii,iAt2,iAt1) * this%deltaDAtom(:,iAt2))
          this%pot21x2Atom(ii,iAt1) = this%pot21x2Atom(ii,iAt1)&
              & - 3.0_dp * sum(this%f30AB(:,:,ii,iAt2,iAt1) * this%deltaQAtom(:,:,iAt2))
          do jj = 1, xyzLen
            this%pot21x1Atom(jj,ii,iAt1) = this%pot21x1Atom(jj,ii,iAt1)&
                & - 3.0_dp * sum(this%f30AB(:,jj,ii,iAt2,iAt1) * this%deltaDAtom(:,iAt2))
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    !> Calculate pot21x2
    this%pot22x2Atom(:,:,:) = 0.0_dp

    !$OMP PARALLEL DO PRIVATE(iAt1, iSp1, iAt2, iSp2, ii, jj) DEFAULT(SHARED) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtoms
      iSp1 = this%species(iAt1)
      if ((this%onDQOCharges(1, iSp1) == 0) .and. (this%onDQOCharges(2, iSp1) == 0)) then
        cycle
      end if
      do iAt2 = 1, nAtoms
        iSp2 = this%species(iAt2)
        if ((this%onDQOCharges(1, iSp2) == 0) .and. (this%onDQOCharges(2, iSp2) == 0)) then
          cycle
        end if
        do ii = 1, xyzLen
          do jj = 1, xyzLen
            this%pot22x2Atom(jj,ii,iAt1) = this%pot22x2Atom(jj,ii,iAt1)&
                & + 6.0_dp * sum(this%f40AB(:,:,jj,ii,iAt2,iAt1) * this%deltaQAtom(:,:,iAt2))
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine updateDQPotentials


  subroutine addMultiPoleHamiltonian(this, env, over, iNeighbour, nNeighbourSK, iSquare, iPair,&
      & orb, hamiltonian, overlap)

    !> class instance
    class(TDftbMultiPole), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Sparse (packed) overlap matrix.
    real(dp), dimension(:), intent(in) :: over

    !> Neighbour indices.
    integer, dimension(0:,:), intent(in) :: iNeighbour

    !> Nr. of neighbours for each atom.
    integer, dimension(:), intent(in) :: nNeighbourSK

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtoms)
    integer, dimension(:), intent(in) :: iSquare

    !> iPair Position of each (neighbour, atom) pair in the sparse matrix.
    !> Shape: (0:maxNeighbour,nAtoms)
    integer, dimension(0:,:), intent(in) :: iPair

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), dimension(:,:), intent(inout), target :: hamiltonian

    !> square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    integer, parameter :: xyzLen = 3, iStart = 1, iEnd = 2, iNOrb = 3
    real(dp) :: tmpadS(xyzLen), tmpSad(xyzLen)
    real(dp) :: tmpaQS(xyzLen, xyzLen), tmpSaQ(xyzLen, xyzLen)
    real(dp) :: tmpCoAtMu1(xyzLen, xyzLen), tmpCoAtNu1(xyzLen, xyzLen)
    real(dp) :: tmpCoAtMu2(xyzLen, xyzLen), tmpCoAtNu2(xyzLen, xyzLen)
    real(dp), allocatable :: tmpOvr(:,:), deltaHam(:,:), tmp

    integer :: matrixSize, nAtoms, maxNOrb
    integer :: iAtMu, iAtNu, iNeigh, iSpMu, iSpNu, nOrbMu, nOrbNu
    integer :: mu, nu, kappa, mm, nn, kk
    integer :: jj, descMuiStart, descMuiEnd, descNuiStart, descNuiEnd
    integer :: iSp1, iSp2, nOrb1, nOrb2
    real(dp), parameter :: tolZero = 1.0e-15_dp

    matrixSize = size(hamiltonian, dim = 1)
    nAtoms = this%nAtoms
    maxNOrb = maxval(orb%nOrbSpecies(:))

    allocate(tmpOvr(matrixSize, matrixSize))
    allocate(deltaHam(matrixSize, matrixSize))
    tmpOvr(:,:) = overlap
    call adjointLowerTriangle(tmpOvr)
    deltaHam(:,:) = 0.0_dp

    loopAtMu: do iAtMu = 1, nAtoms
      iSpMu = this%species(iAtMu)
      nOrbMu = orb%nOrbSpecies(iSpMu)
     !descMu = getDescriptor(iAtMu, iSquare)
      descMuiStart  = iSquare(iAtMu)
      descMuiEnd = iSquare(iAtMu + 1) - 1

     !loopAtNu: do iNeigh = 1, nNeighbourSK(iAtMu)
     !  iAtNu = iNeighbour(iNeigh, iAtMu)
     !  if (iAtNu > iAtMu) then
     !    cycle
     !  end if
      loopAtNu: do iAtNu = 1, iAtMu
        iSpNu = this%species(iAtNu)
        nOrbNu = orb%nOrbSpecies(iSpNu)
       !descNu = getDescriptor(iAtNu, iSquare)
        descNuiStart  = iSquare(iAtNu)
        descNuiEnd = iSquare(iAtNu + 1) - 1

        tmp=maxval(abs(tmpOvr(descMuiStart:descMuiEnd,descNuiStart:descNuiEnd)))
        if (tmp < tolZero) then
          cycle
        end if

        loopmu: do mu = descMuiStart, descMuiEnd
          mm = mu - descMuiStart + 1
          loopnu: do nu = descNuiStart, descNuiEnd
            nn = nu - descNuiStart + 1
            tmpadS(:) = 0.0_dp
            tmpSad(:) = 0.0_dp
            tmpaQS(:,:) = 0.0_dp
            tmpSaQ(:,:) = 0.0_dp
            do kk = 1, nOrbMu
              kappa = descMuiStart + kk - 1
              tmpadS(:) = tmpadS(:) + this%atomicDIntgrl(:,mm,kk,iSpMu) * tmpOvr(kappa,nu)
              tmpaQS(:,:) = tmpaQS(:,:) + this%atomicQIntgrl(:,:,mm,kk,iSpMu) * tmpOvr(kappa,nu)
            end do
            do kk = 1, nOrbNu
              kappa = descNuiStart + kk - 1
              tmpSad(:) = tmpSad(:) + tmpOvr(mu,kappa) * this%atomicDIntgrl(:,kk,nn,iSpNu)
              tmpSaQ(:,:) = tmpSaQ(:,:) + tmpOvr(mu,kappa) * this%atomicQIntgrl(:,:,kk,nn,iSpNu)
            end do

            !>add Monopole-Dipole contribution
            deltaHam(mu,nu) = deltaHam(mu,nu) - (this%pot10x1Atom(iAtMu)&
                & + this%pot10x1Atom(iAtNu)) * tmpOvr(mu,nu)
            deltaHam(mu,nu) = deltaHam(mu,nu) + sum(this%pot10x0Atom(:,iAtMu) * tmpadS(:)) &
                                            & + sum(this%pot10x0Atom(:,iAtNu) * tmpSad(:))

            !>add Dipole-Dipole contribution
            deltaHam(mu,nu) = deltaHam(mu,nu) + sum(this%pot11x1Atom(:,iAtMu) * tmpadS(:)) &
                                            & + sum(this%pot11x1Atom(:,iAtNu) * tmpSad(:))

            !>add Monopole-Quadrupole contribution
            deltaHam(mu,nu) = deltaHam(mu,nu) + (this%pot20x2Atom(iAtMu)&
                & + this%pot20x2Atom(iAtNu)) * tmpOvr(mu,nu)
            deltaHam(mu,nu) = deltaHam(mu,nu) + sum(this%pot20x0Atom(:,:,iAtMu) * tmpaQS(:,:)) &
                                            & + sum(this%pot20x0Atom(:,:,iAtNu) * tmpSaQ(:,:))

            !>add Dipole-Quadrupole contribution
            deltaHam(mu,nu) = deltaHam(mu,nu) - sum(this%pot21x2Atom(:,iAtMu) * tmpadS(:)) &
                                            & - sum(this%pot21x2Atom(:,iAtNu) * tmpSad(:))
            deltaHam(mu,nu) = deltaHam(mu,nu) + sum(this%pot21x1Atom(:,:,iAtMu) * tmpaQS(:,:)) &
                                            & + sum(this%pot21x1Atom(:,:,iAtNu) * tmpSaQ(:,:))

            !>add Quadrupole-Quadrupole contribution
            deltaHam(mu,nu) = deltaHam(mu,nu) + sum(this%pot22x2Atom(:,:,iAtMu) * tmpaQS(:,:)) &
                                            & + sum(this%pot22x2Atom(:,:,iAtNu) * tmpSaQ(:,:))

          end do loopnu
        end do loopmu
      end do loopAtNu
    end do loopAtMu

    hamiltonian(:,:) = hamiltonian + 0.5_dp * deltaHam

    deallocate(tmpOvr, deltaHam)


  end subroutine addMultiPoleHamiltonian


  !> Add the MultiPole-Energy contribution to the total energy
  subroutine addMultiPoleEnergy(this, energyPerAtom, energyMD, energyDD, energyMQ, energyDQ,&
      & energyQQ, energyTT)

    !> RangeSep class instance
    class(TDftbMultiPole), intent(inout) :: this

    real(dp), intent(out) :: energyPerAtom(:)

    !> total energy
    real(dp), intent(inout) :: energyMD, energyDD, energyMQ, energyDQ, energyQQ, energyTT

    call getEnergyPerAtom(this, energyPerAtom)

    energyMD = this%energyMD
    energyDD = this%energyDD
    energyMQ = this%energyMQ
    energyDQ = this%energyDQ
    energyQQ = this%energyQQ
    energyTT = this%energyTT

  end subroutine addMultiPoleEnergy


  !> Returns energy per atom.
  subroutine getEnergyPerAtom(this, energyPerAtom)

    class(TDftbMultiPole), intent(inout) :: this

    real(dp), intent(out) :: energyPerAtom(:)

    real(dp), allocatable :: EnergyMDAtom(:), EnergyDDAtom(:), EnergyMQAtom(:)
    real(dp), allocatable :: EnergyDQAtom(:), EnergyQQAtom(:)
    integer, parameter :: xyzLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    integer :: nAtoms
    integer :: iAt1
    @:ASSERT(size(energyPerAtom) == this%nAtoms)

    nAtoms = this%nAtoms
    allocate(EnergyMDAtom(nAtoms))
    allocate(EnergyDDAtom(nAtoms))
    allocate(EnergyMQAtom(nAtoms))
    allocate(EnergyDQAtom(nAtoms))
    allocate(EnergyQQAtom(nAtoms))
    EnergyMDAtom(:) = 0.0_dp
    EnergyDDAtom(:) = 0.0_dp
    EnergyMQAtom(:) = 0.0_dp
    EnergyDQAtom(:) = 0.0_dp
    EnergyQQAtom(:) = 0.0_dp

    do iAt1 = 1, nAtoms
      !>Monopole-Dipole contribution
      EnergyMDAtom(iAt1) = sum(this%pot10x0Atom(:,iAt1) * this%deltaDAtom(:,iAt1))
      !>Dipole-Dipole contribution
      EnergyDDAtom(iAt1) = 0.5_dp * sum(this%pot11x1Atom(:,iAt1) * this%deltaDAtom(:,iAt1))
      !>Monopole-Quadrupole contribution
      EnergyMQAtom(iAt1) = sum(this%pot20x0Atom(:,:,iAt1) * this%deltaQAtom(:,:,iAt1))
      !>Dipole-Quadrupole contribution
      EnergyDQAtom(iAt1) = sum(this%pot21x1Atom(:,:,iAt1) * this%deltaQAtom(:,:,iAt1))
      !>Quadrupole-Quadrupole contribution
      EnergyQQAtom(iAt1) = 0.5_dp * sum(this%pot22x2Atom(:,:,iAt1) * this%deltaQAtom(:,:,iAt1))
    end do

    energyPerAtom(:) = EnergyMDAtom + EnergyDDAtom + EnergyMQAtom + EnergyDQAtom + EnergyQQAtom
    this%energyMD = sum(EnergyMDAtom)
    this%energyDD = sum(EnergyDDAtom)
    this%energyMQ = sum(EnergyMQAtom)
    this%energyDQ = sum(EnergyDQAtom)
    this%energyQQ = sum(EnergyQQAtom)
    this%energyTT = this%energyMD + this%energyDD + this%energyMQ  + this%energyDQ  + this%energyQQ

    deallocate(EnergyMDAtom)
    deallocate(EnergyDDAtom)
    deallocate(EnergyMQAtom)
    deallocate(EnergyDQAtom)
    deallocate(EnergyQQAtom)


  end subroutine getEnergyPerAtom

  subroutine addMultiPoleGradients(this, gradients, derivator, rhoSqr, skHamCont, skOverCont,&
      & coords, species, orb, iSquare, ovrlapMat, iNeighbour, nNeighbourSK)

    !> class instance
    class(TDftbMultiPole), intent(inout) :: this

    !> energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> density matrix
    real(dp), intent(in) :: rhoSqr(:,:,:)

    !> sparse hamiltonian (non-scc)
    type(TSlakoCont), intent(in) :: skHamCont

    !> sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> chemical species of atoms
    integer, intent(in) :: species(:)

    !> orbital information for system
    type(TOrbitals), intent(in) :: orb

    !> index for dense arrays
    integer, intent(in) :: iSquare(:)

    !> overlap matrix
    real(dp), intent(in) :: ovrlapMat(:,:)

    !> neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> number of atoms neighbouring each site where the overlap is non-zero
    integer, intent(in) :: nNeighbourSK(:)

    !> differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    real(dp), allocatable :: tmpOvr(:,:), tmpRho(:,:), derivs(:,:)

    integer, parameter :: xyzLen = 3, iStart = 1, iEnd = 2, iNOrb = 3
    integer :: matrixSize, nAtoms, maxNOrb
    integer :: iAt1, iAt2, iSp1, iSp2
    integer :: ii, jj, ll, mu, nu, kappa, mm, nn, kk
    real(dp) :: sPrimeTmp12(orb%mOrb,orb%mOrb,3)
    real(dp) :: tmpDerivMuNu, tmpDeriv(xyzLen), tmpVx3(xyzLen), tmpM3x3(xyzLen, xyzLen)

    real(dp):: tmp
    real(dp), parameter :: tolZero = 1.0e-15_dp

    matrixSize = size(rhoSqr, dim=1)
    nAtoms = this%nAtoms
    maxNOrb = maxval(orb%nOrbSpecies(:))

    allocate(tmpOvr(matrixSize, matrixSize))
    allocate(tmpRho(matrixSize, matrixSize))
    allocate(derivs(xyzLen, size(gradients, dim = 2)))
    derivs = 0.0_dp
    tmpDeriv = 0.0_dp
    tmpOvr(:,:) = ovrlapMat
    call adjointLowerTriangle(tmpOvr)
    tmpRho(:,:) = rhoSqr(:,:,1)
    call symmetrizeSquareMatrix(tmpRho)

    !> Calculate pot30x0, pot31x1, and pot32x2
    this%pot30x0Atom(:,:,:,:) = 0.0_dp
    this%pot31x1Atom(:,:,:,:) = 0.0_dp
    this%pot32x2Atom(:,:,:,:) = 0.0_dp

    !$OMP PARALLEL DO PRIVATE(iAt1, iSp1, iAt2, iSp2, ii, jj) DEFAULT(SHARED) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtoms
      iSp1 = this%species(iAt1)
      do iAt2 = 1, nAtoms
        iSp2 = this%species(iAt2)
        this%pot30x0Atom(:,:,:,iAt1) = this%pot30x0Atom(:,:,:,iAt1)&
            & + this%f30AB(:,:,:,iAt2,iAt1) * this%deltaMAtom(iAt2)
        if ((this%onDQOCharges(1, iSp1) == 0 .and. this%onDQOCharges(2, iSp1) == 0 ) &
           & .or. (this%onDQOCharges(1, iSp2) == 0 .and. this%onDQOCharges(2, iSp2) == 0)) then
          cycle
        end if
        do ii = 1, xyzLen
          do jj = 1, xyzLen
            do ll = 1, xyzLen
              this%pot31x1Atom(ll,jj,ii,iAt1) = this%pot31x1Atom(ll,jj,ii,iAt1)&
                  & - 4.0_dp * sum(this%f40AB(:,ll,jj,ii,iAt2,iAt1) * this%deltaDAtom(:,iAt2))
              this%pot32x2Atom(ll,jj,ii,iAt1) = this%pot32x2Atom(ll,jj,ii,iAt1)&
                  & + 10.0_dp * sum(this%f50AB(:,:,ll,jj,ii,iAt2,iAt1) * this%deltaQAtom(:,:,iAt2))
            end do
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO


    derivs(:,:) = 0.0_dp

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(iAt1, iSp1, iAt2, iSp2) &
    !$OMP PRIVATE(ii, mu, mm, nu, nn, kappa, kk ) &
    !$OMP PRIVATE(tmp, tmpDerivMuNu) &
    !$OMP PRIVATE(sPrimeTmp12, tmpDeriv, tmpVx3, tmpM3x3)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, this%nAtoms
      iSp1 = species(iAt1)

      !> M-D contribution
      derivs(:,iAt1) = derivs(:,iAt1) + this%pot11x1Atom(:,iAt1) * this%deltaMAtom(iAt1)
      !> M-Q contribution
      derivs(:,iAt1) = derivs(:,iAt1) - this%pot21x2Atom(:,iAt1) * this%deltaMAtom(iAt1)
      do ii = 1, xyzLen
        !> M-D contribution
        derivs(ii,iAt1) = derivs(ii,iAt1)&
            & + 2.0_dp * sum(this%pot20x0Atom(:,ii,iAt1) * this%deltaDAtom(:,iAt1))
        !> D-D contribution
        derivs(ii,iAt1) = derivs(ii,iAt1)&
            & + 2.0_dp * sum(this%pot21x1Atom(:,ii,iAt1) * this%deltaDAtom(:,iAt1))
        !> D-Q contribution
        derivs(ii,iAt1) = derivs(ii,iAt1)&
            & + 2.0_dp * sum(this%pot22x2Atom(:,ii,iAt1) * this%deltaDAtom(:,iAt1))
        !> M-Q contribution
        derivs(ii,iAt1) = derivs(ii,iAt1)&
            & + 3.0_dp * sum(this%pot30x0Atom(:,:,ii,iAt1) * this%deltaQAtom(:,:,iAt1))
        !> D-Q contribution
        derivs(ii,iAt1) = derivs(ii,iAt1)&
            & + 3.0_dp * sum(this%pot31x1Atom(:,:,ii,iAt1) * this%deltaQAtom(:,:,iAt1))
        !> Q-Q contribution
        derivs(ii,iAt1) = derivs(ii,iAt1)&
            & + 3.0_dp * sum(this%pot32x2Atom(:,:,ii,iAt1) * this%deltaQAtom(:,:,iAt1))
      end do

      tmpDeriv = 0.0_dp
      do iAt2 = 1, this%nAtoms
        iSp2 = species(iAt2)
        if (iAt2 == iAt1) then
          cycle
        end if

        tmp=maxval(abs(tmpOvr(iSquare(iAt1):iSquare(iAt1 + 1) - 1,&
            & iSquare(iAt2):iSquare(iAt2 + 1) - 1)))
        if (tmp < tolZero) then
          cycle
        end if

        sPrimeTmp12 = 0.0_dp
        call derivator%getFirstDeriv(sPrimeTmp12, skOverCont, coords, species, iAt1, iAt2, orb)

        tmp=maxval(abs(sPrimeTmp12(:,:,:)))
        if (tmp < tolZero) then
          cycle
        end if

        do mu = iSquare(iAt1), iSquare(iAt1 + 1) - 1
          mm = mu - iSquare(iAt1) + 1
          do nu = iSquare(iAt2), iSquare(iAt2 + 1) - 1
            nn = nu - iSquare(iAt2) + 1

            tmpDerivMuNu = 0.0_dp
            !> M-D contribution
            tmpDerivMuNu = tmpDerivMuNu - (this%pot10x1Atom(iAt1)&
                & + this%pot10x1Atom(iAt2)) * tmpRho(nu,mu)
            !> M-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + (this%pot20x2Atom(iAt1)&
                & + this%pot20x2Atom(iAt2)) * tmpRho(nu,mu)

            tmpVx3 = 0.0_dp
            tmpM3x3 = 0.0_dp
            do kappa = iSquare(iAt1), iSquare(iAt1 + 1) - 1
              kk = kappa - iSquare(iAt1) + 1
              tmpVx3(:) = tmpVx3(:) + tmpRho(nu,kappa) * this%atomicDIntgrl(:,kk,mm,iSp1)
              tmpM3x3(:,:) = tmpM3x3(:,:) + tmpRho(nu,kappa) * this%atomicQIntgrl(:,:,kk,mm,iSp1)
            end do
            !> M-D contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot10x0Atom(:,iAt1) * tmpVx3(:))
            !> D-D contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot11x1Atom(:,iAt1) * tmpVx3(:))
            !> D-Q contribution
            tmpDerivMuNu = tmpDerivMuNu - sum(this%pot21x2Atom(:,iAt1) * tmpVx3(:))
            !> M-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot20x0Atom(:,:,iAt1) * tmpM3x3(:,:))
            !> D-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot21x1Atom(:,:,iAt1) * tmpM3x3(:,:))
            !> Q-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot22x2Atom(:,:,iAt1) * tmpM3x3(:,:))

            tmpVx3 = 0.0_dp
            tmpM3x3 = 0.0_dp
            do kappa = iSquare(iAt2), iSquare(iAt2 + 1) - 1
              kk = kappa - iSquare(iAt2) + 1
              tmpVx3(:) = tmpVx3(:) + tmpRho(kappa,mu) * this%atomicDIntgrl(:,nn,kk,iSp2)
              tmpM3x3(:,:) = tmpM3x3(:,:) + tmpRho(kappa,mu) * this%atomicQIntgrl(:,:,nn,kk,iSp2)
            end do
            !> M-D contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot10x0Atom(:,iAt2) * tmpVx3(:))
            !> D-D contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot11x1Atom(:,iAt2) * tmpVx3(:))
            !> D-Q contribution
            tmpDerivMuNu = tmpDerivMuNu - sum(this%pot21x2Atom(:,iAt2) * tmpVx3(:))
            !> M-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot20x0Atom(:,:,iAt2) * tmpM3x3(:,:))
            !> D-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot21x1Atom(:,:,iAt2) * tmpM3x3(:,:))
            !> Q-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot22x2Atom(:,:,iAt2) * tmpM3x3(:,:))

            tmpDeriv(:) = tmpDeriv(:) + tmpDerivMuNu * sPrimeTmp12(nn,mm,:)

          end do
        end do
      end do
      derivs(:,iAt1) = derivs(:,iAt1) + tmpDeriv(:)

    end do
    !$OMP END DO
    !$OMP END PARALLEL

    gradients(:,:) = gradients(:,:) + derivs(:,:)

    deallocate(tmpOvr, tmpRho, derivs)


  end subroutine addMultiPoleGradients


  subroutine addAtomicDipoleMoment(this, dipoleMoment)

    !> class instance
    class(TDftbMultiPole), intent(inout) :: this

    !> energy gradients
    real(dp), intent(inout) :: dipoleMoment(:)

    integer :: nAtoms
    integer :: iAt1

    nAtoms = this%nAtoms
    do iAt1 = 1, nAtoms
      dipoleMoment(:) = dipoleMoment(:) - this%deltaDAtom(:,iAt1)
    end do

  end subroutine addAtomicDipoleMoment


  subroutine addAtomicQuadrupoleMoment(this, quadrupoleMoment)

    !> class instance
    class(TDftbMultiPole), intent(inout) :: this

    !> energy gradients
    real(dp), intent(inout) :: quadrupoleMoment(:,:)

    integer :: nAtoms
    integer :: iAt1

    nAtoms = this%nAtoms
    do iAt1 = 1, nAtoms
      quadrupoleMoment(:,:) = quadrupoleMoment(:,:) - this%deltaQAtom(:,:,iAt1)
    end do

  end subroutine addAtomicQuadrupoleMoment


  subroutine makeA3oB3(M3x3, A3, B3)

    real(dp), intent(out) :: M3x3(:,:)
    real(dp), intent(in) :: A3(:), B3(:)
    integer :: ii, jj, dimSize

    dimSize=3
    @:ASSERT(size(A3) == dimSize)
    @:ASSERT(size(B3) == dimSize)
    @:ASSERT(size(M3x3, dim=1) == dimSize)
    @:ASSERT(size(M3x3, dim=2) == dimSize)
    do ii = 1, dimSize
      do jj = 1, dimSize
        M3x3(ii, jj) = A3(ii) * B3(jj)
      end do
    end do

  end subroutine makeA3oB3


  subroutine makeA3oB3oC3(M3x3x3, A3, B3, C3)

    real(dp), intent(out) :: M3x3x3(:,:,:)
    real(dp), intent(in) :: A3(:), B3(:), C3(:)
    integer :: ii, jj, ll, dimSize

    dimSize=3
    @:ASSERT(size(A3) == dimSize)
    @:ASSERT(size(B3) == dimSize)
    @:ASSERT(size(C3) == dimSize)
    @:ASSERT(size(M3x3x3, dim=1) == dimSize)
    @:ASSERT(size(M3x3x3, dim=2) == dimSize)
    @:ASSERT(size(M3x3x3, dim=3) == dimSize)
    do ii = 1, dimSize
      do jj = 1, dimSize
        do ll = 1, dimSize
          M3x3x3(ii,jj,ll) = A3(ii) * B3(jj) * C3(ll)
        end do
      end do
    end do

  end subroutine makeA3oB3oC3


  subroutine makeA3oB3x3(M3x3x3, A3, B3x3)

    real(dp), intent(out) :: M3x3x3(:,:,:)
    real(dp), intent(in) :: A3(:), B3x3(:,:)
    integer :: ii, jj, ll, dimSize

    dimSize=3
    @:ASSERT(size(A3) == dimSize)
    @:ASSERT(size(B3x3, dim=1) == dimSize)
    @:ASSERT(size(B3x3, dim=2) == dimSize)
    @:ASSERT(size(M3x3x3, dim=1) == dimSize)
    @:ASSERT(size(M3x3x3, dim=2) == dimSize)
    @:ASSERT(size(M3x3x3, dim=3) == dimSize)
    do ii = 1, dimSize
      do jj = 1, dimSize
        do ll = 1, dimSize
          M3x3x3(ii,jj,ll) = A3(ii) * B3x3(jj,ll)
        end do
      end do
    end do

  end subroutine makeA3oB3x3


  subroutine makeA3x3oB3(M3x3x3, A3x3, B3)

    real(dp), intent(out) :: M3x3x3(:,:,:)
    real(dp), intent(in) :: A3x3(:,:), B3(:)
    integer :: ii, jj, ll, dimSize

    dimSize=3
    @:ASSERT(size(B3) == dimSize)
    @:ASSERT(size(A3x3, dim=1) == dimSize)
    @:ASSERT(size(A3x3, dim=2) == dimSize)
    @:ASSERT(size(M3x3x3, dim=1) == dimSize)
    @:ASSERT(size(M3x3x3, dim=2) == dimSize)
    @:ASSERT(size(M3x3x3, dim=3) == dimSize)
    do ii = 1, dimSize
      do jj = 1, dimSize
        do ll = 1, dimSize
          M3x3x3(ii,jj,ll) = A3x3(ii,jj) * B3(ll)
        end do
      end do
    end do

  end subroutine makeA3x3oB3


  subroutine makeA3x3otB3(M3x3x3, A3x3, B3)

    real(dp), intent(out) :: M3x3x3(:,:,:)
    real(dp), intent(in) :: A3x3(:,:), B3(:)
    integer :: ii, jj, ll, dimSize

    dimSize=3
    @:ASSERT(size(B3) == dimSize)
    @:ASSERT(size(A3x3, dim=1) == dimSize)
    @:ASSERT(size(A3x3, dim=2) == dimSize)
    @:ASSERT(size(M3x3x3, dim=1) == dimSize)
    @:ASSERT(size(M3x3x3, dim=2) == dimSize)
    @:ASSERT(size(M3x3x3, dim=3) == dimSize)
    do ii = 1, dimSize
      do jj = 1, dimSize
        do ll = 1, dimSize
          M3x3x3(ii,jj,ll) = A3x3(ii,ll) * B3(jj)
        end do
      end do
    end do

  end subroutine makeA3x3otB3


  subroutine makeA3x3oB3x3(M3x3x3x3, A3x3, B3x3)

    real(dp), intent(out) :: M3x3x3x3(:,:,:,:)
    real(dp), intent(in) :: A3x3(:,:), B3x3(:,:)
    integer :: ii, jj, ll, mm, dimSize

    dimSize=3
    @:ASSERT(size(B3x3, dim=1) == dimSize)
    @:ASSERT(size(B3x3, dim=2) == dimSize)
    @:ASSERT(size(A3x3, dim=1) == dimSize)
    @:ASSERT(size(A3x3, dim=2) == dimSize)
    @:ASSERT(size(M3x3x3x3, dim=1) == dimSize)
    @:ASSERT(size(M3x3x3x3, dim=2) == dimSize)
    @:ASSERT(size(M3x3x3x3, dim=3) == dimSize)
    @:ASSERT(size(M3x3x3x3, dim=4) == dimSize)
    do ii = 1, dimSize
      do jj = 1, dimSize
        do ll = 1, dimSize
          do mm = 1, dimSize
            M3x3x3x3(ii,jj,ll,mm) = A3x3(ii,jj) * B3x3(ll,mm)
          end do
        end do
      end do
    end do

  end subroutine makeA3x3oB3x3


  subroutine makeA3x3ohB3x3(M3x3x3x3, A3x3, B3x3)

    real(dp), intent(out) :: M3x3x3x3(:,:,:,:)
    real(dp), intent(in) :: A3x3(:,:), B3x3(:,:)
    integer :: ii, jj, ll, mm, dimSize

    dimSize=3
    @:ASSERT(size(B3x3, dim=1) == dimSize)
    @:ASSERT(size(B3x3, dim=2) == dimSize)
    @:ASSERT(size(A3x3, dim=1) == dimSize)
    @:ASSERT(size(A3x3, dim=2) == dimSize)
    @:ASSERT(size(M3x3x3x3, dim=1) == dimSize)
    @:ASSERT(size(M3x3x3x3, dim=2) == dimSize)
    @:ASSERT(size(M3x3x3x3, dim=3) == dimSize)
    @:ASSERT(size(M3x3x3x3, dim=4) == dimSize)
    do ii = 1, dimSize
      do jj = 1, dimSize
        do ll = 1, dimSize
          do mm = 1, dimSize
            M3x3x3x3(ii,jj,ll,mm) = A3x3(ii,ll) * B3x3(jj,mm)
          end do
        end do
      end do
    end do

  end subroutine makeA3x3ohB3x3


  subroutine makeA3x3otB3x3(M3x3x3x3, A3x3, B3x3)

    real(dp), intent(out) :: M3x3x3x3(:,:,:,:)
    real(dp), intent(in) :: A3x3(:,:), B3x3(:,:)
    integer :: ii, jj, ll, mm, dimSize

    dimSize=3
    @:ASSERT(size(B3x3, dim=1) == dimSize)
    @:ASSERT(size(B3x3, dim=2) == dimSize)
    @:ASSERT(size(A3x3, dim=1) == dimSize)
    @:ASSERT(size(A3x3, dim=2) == dimSize)
    @:ASSERT(size(M3x3x3x3, dim=1) == dimSize)
    @:ASSERT(size(M3x3x3x3, dim=2) == dimSize)
    @:ASSERT(size(M3x3x3x3, dim=3) == dimSize)
    @:ASSERT(size(M3x3x3x3, dim=4) == dimSize)
    do ii = 1, dimSize
      do jj = 1, dimSize
        do ll = 1, dimSize
          do mm = 1, dimSize
            M3x3x3x3(ii,jj,ll,mm) = A3x3(ii,mm) * B3x3(jj,ll)
          end do
        end do
      end do
    end do

  end subroutine makeA3x3otB3x3

  !> copy lower triangle to upper for a square matrix
  subroutine symmetrizeSquareMatrix(matrix)

    !> matrix to symmetrize
    real(dp), intent(inout) :: matrix(:,:)
    integer :: ii, matSize

    @:ASSERT(size(matrix, dim=1) == size(matrix, dim=2))
    matSize = size(matrix, dim = 1)

   !!$OMP PARALLEL DO PRIVATE(ii) DEFAULT(SHARED) SCHEDULE(RUNTIME)
    do ii = 1, matSize - 1
      matrix(ii,ii+1:matSize) = matrix(ii + 1 : matSize, ii)
    end do
   !!$OMP END PARALLEL DO

  end subroutine symmetrizeSquareMatrix


  !> location of relevant atomic block indices in a dense matrix
  pure function getDescriptor(iAt, iSquare) result(desc)

    !> relevant atom
    integer, intent(in) :: iAt

    !> indexing array for start of atom orbitals
    integer, intent(in) :: iSquare(:)

    !> resulting location ranges
    integer :: desc(3)

    desc(:) = [iSquare(iAt), iSquare(iAt + 1) - 1, iSquare(iAt + 1) - iSquare(iAt)]

  end function getDescriptor

end module dftbp_dftb_multipole
