!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Routines implementing the (One-center approximation) multipole expansion for the 2nd order DFTB.
module dftbp_dftb_multiexpan
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_globalenv, only : stdOut
  use dftbp_dftb_nonscc, only : TNonSccDiff
  use dftbp_dftb_shortgammafuncs, only : expGammaPrime, expGammaDoublePrime, expGammaTriplePrime,&
      & expGammaQuadruplePrime, expGammaQuintuplePrime
  use dftbp_dftb_slakocont, only : TSlakoCont
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_multipole, only : TMultipole
  use dftbp_math_matrixops, only : adjointLowerTriangle
  use dftbp_dftb_periodic, only : TNeighbourList, getNrOfNeighbours
  implicit none
  private

  public :: TDftbMultiExpanInp, TDftbMultiExpan, dftbMultiExpan_init

  !> Input for the MultiExpan module
  type TDftbMultiExpanInp

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

  end type TDftbMultiExpanInp

  !> Internal status of TMultiExpan.
  type TDftbMultiExpan
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
    procedure :: pullDeltaDQAtom
    procedure :: pushDeltaDQAtom
    procedure :: updateDeltaDQAtom
    procedure :: updateDQPotentials
    procedure :: addMultiExpanHamiltonian
    procedure :: addMultiExpanEnergy
    procedure :: addMultiExpanGradients
    procedure :: addAtomicDipoleMoment
    procedure :: addAtomicQuadrupoleMoment
    procedure :: getMultiExpanInfo
    procedure :: getOrbitalEquiv
  end type TDftbMultiExpan

  !> Number of dipole components used in tblite library (x, y, z)
  integer, parameter :: dimDipole = 3

  !> Number of quadrupole components used in tblite library (xx, xy, yy, xz, yz, zz)
  integer, parameter :: dimQuadrupole = 6

contains

  !> Get information on required multipolar contributions
  subroutine getMultiExpanInfo(this, nDipole, nQuadrupole)

    !> Data structure
    class(TDftbMultiExpan), intent(in) :: this

    !> Number of dipole moment components
    integer, intent(out) :: nDipole

    !> Number of quadrupole moment components
    integer, intent(out) :: nQuadrupole

    integer, parameter :: dimDipole = 3
    integer, parameter :: dimQuadrupole = 6

    nDipole = dimDipole
    nQuadrupole = dimQuadrupole
  end subroutine getMultiExpanInfo


  !> Returns the equivalence to get the correct mixing of charge dependent contributions
  subroutine getOrbitalEquiv(this, equivDip, equivQuad)

    !> Data structure
    class(TDftbMultiExpan), intent(inout) :: this

    !> The equivalence vector for cumulative atomic dipole populations
    integer, intent(out) :: equivDip(:,:)

    !> The equivalence vector for cumulative atomic quadrupole populations
    integer, intent(out) :: equivQuad(:,:)

    integer :: nAtom, iCount, iSpin, nSpin, iAt, iSp, ii

    nAtom = size(equivDip, dim=2)
    nSpin = 1

    equivDip(:,:) = 0
    equivQuad(:,:) = 0

    ! equivDip
    iCount = 0
    do iAt = 1, nAtom
      do ii = 1, dimDipole
        iCount = iCount + 1
        equivDip(ii, iAt) = iCount
      end do
    end do

    ! equivQuad
    iCount = 0
    do iAt = 1, nAtom
      do ii = 1, dimQuadrupole
        iCount = iCount + 1
        equivQuad(ii, iAt) = iCount
      end do
    end do

  end subroutine getOrbitalEquiv


  !> Initializes instance.
  subroutine dftbMultiExpan_init(this, inp)

    !> Instance.
    type(TDftbMultiExpan), intent(out) :: this

    !> Input data.
    type(TDftbMultiExpanInp), intent(in) :: inp

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

  end subroutine dftbMultiExpan_init


  !> Updates data structures if there are changed coordinates for the instance.
  subroutine updateCoords(this, coords)

    !> class instance
    class(TDftbMultiExpan), intent(inout) :: this

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

  subroutine updateDeltaDQAtom(this, over, rho, orb, iNeighbour, nNeighbourSK, img2CentCell, iPair, q0)

    !> class instance
    class(TDftbMultiExpan), intent(inout) :: this

    !> Overlap matrix in packed format
    real(dp), intent(in) :: over(:)

    !> Density matrix in Packed format
    real(dp), intent(in) :: rho(:)

    !> Information about the orbitals
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours of each real atom (central cell)
    integer, intent(in) :: iNeighbour(0:,:)

    !> List of neighbours for each atom, starting at 0 for itself
    integer, intent(in) :: nNeighbourSK(:)

    !> Indexing array to convert images of atoms back into their number in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:,:)

    real(dp), allocatable :: mqPerOrbital(:,:)
    real(dp), allocatable :: dqPerOrbital(:,:,:)
    real(dp), allocatable :: qqPerOrbital(:,:,:,:)
    real(dp), allocatable :: mq(:)
    real(dp), allocatable :: dq(:,:)
    real(dp), allocatable :: qq(:,:,:)

    real(dp), allocatable :: tmpOvr(:,:), tmpDRho(:,:)
    integer, parameter :: xyzLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    integer :: matrixSize, nAtoms, maxNOrb
    integer :: iAt1, iAt2, iSp1, iSp2
    integer :: ii, jj, mu, nu, mm, nn, kk, iStartOrb, iEndOrb

    real(dp) :: tmpVx3(xyzLen), tmpM3x3(xyzLen, xyzLen)
    real(dp) :: dPSSdP, tmpTrace

    integer :: iOrig
    integer :: iIter, iNeigh
    integer :: nAtom, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: sqrTmp(orb%mOrb,orb%mOrb)
    real(dp) :: sqrTmpRho(orb%mOrb,orb%mOrb)
    real(dp) :: sqrTmpOver(orb%mOrb,orb%mOrb)
    real(dp) :: mulTmp(orb%mOrb**2)
    integer, allocatable :: iterIndices(:)

    ! nAtom = size(orb%nOrbAtom)
    nAtom = this%nAtoms
    @:ASSERT(size(over) == size(rho))

    ! matrixSize = size(overlap, dim = 1)
    ! maxNOrb = maxval(orb%nOrbSpecies(:))
    ! maxNOrb = orb%mOrb

    allocate(mqPerOrbital(orb%mOrb, nAtom))
    allocate(mq(nAtom))
    mqPerOrbital(:,:) = 0.0_dp
    mq(:) = 0.0_dp

    this%deltaDAtom(:,:) = 0.0_dp
    this%deltaQAtom(:,:,:) = 0.0_dp

    do iAtom1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAtom1)
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        iOrig = iPair(iNeigh, iAtom1) + 1

        sqrTmp(:,:) = 0.0_dp
        mulTmp(:) = 0.0_dp
        mulTmp(1:nOrb1*nOrb2) = over(iOrig:iOrig+nOrb1*nOrb2-1) * rho(iOrig:iOrig+nOrb1*nOrb2-1)
        sqrTmp(1:nOrb2,1:nOrb1) = reshape(mulTmp(1:nOrb1*nOrb2), (/nOrb2,nOrb1/))
        mqPerOrbital(1:nOrb1,iAtom1) = mqPerOrbital(1:nOrb1,iAtom1) + sum(sqrTmp(1:nOrb2,1:nOrb1), dim=1)
        ! Add contribution to the other triangle sum, using the symmetry, but only when off diagonal
        if (iAtom1 /= iAtom2f) then
          mqPerOrbital(1:nOrb2,iAtom2f) = mqPerOrbital(1:nOrb2,iAtom2f) + sum(sqrTmp(1:nOrb2,1:nOrb1), dim=2)
        end if
      end do

      iSp1 = this%species(iAtom1)
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        iSp2 = this%species(iAtom2f)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        iOrig = iPair(iNeigh, iAtom1) + 1

        sqrTmpOver(:,:) = 0.0_dp
        sqrTmpRho(:,:) = 0.0_dp
        sqrTmpOver(1:nOrb2,1:nOrb1) = reshape(over(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
        sqrTmpRho(1:nOrb2,1:nOrb1) = reshape(rho(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))

        tmpVx3(:) = 0.0_dp
        tmpM3x3(:,:) = 0.0_dp
        do mu = 1, nOrb1
          do nu =  1, nOrb1
            do kk =  1, nOrb2
              tmpVx3(:) =  tmpVx3(:) + sqrTmpRho(kk,mu) * sqrTmpOver(kk,nu) * this%atomicDIntgrl(:,nu,mu,iSp1)
              tmpM3x3(:,:) = tmpM3x3(:,:) + sqrTmpRho(kk,mu) * sqrTmpOver(kk,nu) * this%atomicQIntgrl(:,:,nu,mu,iSp1)
            end do
          end do
        end do
        this%deltaDAtom(:,iAtom1) = this%deltaDAtom(:,iAtom1) + tmpVx3(:)
        this%deltaQAtom(:,:,iAtom1) = this%deltaQAtom(:,:,iAtom1) + tmpM3x3(:,:)

        ! Add contribution to the other triangle sum, using the symmetry, but only when off diagonal
        if (iAtom1 /= iAtom2f) then
          tmpVx3(:) = 0.0_dp
          tmpM3x3(:,:) = 0.0_dp
          do mu = 1, nOrb2
            do nu =  1, nOrb2
              do kk =  1, nOrb1
                tmpVx3(:) =  tmpVx3(:) + sqrTmpRho(mu,kk) * sqrTmpOver(nu,kk) * this%atomicDIntgrl(:,nu,mu,iSp2)
                tmpM3x3(:,:) = tmpM3x3(:,:) + sqrTmpRho(mu,kk) * sqrTmpOver(nu,kk) * this%atomicQIntgrl(:,:,nu,mu,iSp2)
              end do
            end do
          end do
          this%deltaDAtom(:,iAtom2f) = this%deltaDAtom(:,iAtom2f) + tmpVx3(:)
          this%deltaQAtom(:,:,iAtom2f) = this%deltaQAtom(:,:,iAtom2f) + tmpM3x3(:,:)
        end if
      end do

    end do
    mq(:) = mq + sum(mqPerOrbital, dim=1)

    ! do iAtom1 = 1, nAtom
    !   nOrb1 = orb%nOrbAtom(iAtom1)
    !   iSp1 = this%species(iAtom1)
    !   tmpVx3(:) = 0.0_dp
    !   tmpM3x3(:,:) = 0.0_dp
    !   do mu = 1, nOrb1
    !     tmpVx3(:) = tmpVx3(:) + this%atomicDIntgrl(:,mu,mu,iSp1) * q0(mu,iAtom1,1)
    !     tmpM3x3(:,:) = tmpM3x3(:,:) + this%atomicQIntgrl(:,:,mu,mu,iSp1) * q0(mu,iAtom1,1)
    !   end do
    !   this%deltaDAtom(:,iAtom1) = this%deltaDAtom(:,iAtom1) - tmpVx3(:)
    !   this%deltaQAtom(:,:,iAtom1) = this%deltaQAtom(:,:,iAtom1) - tmpM3x3(:,:)
    ! end do

  end subroutine updateDeltaDQAtom

  subroutine pullDeltaDQAtom(this, multiExpanData)

    !> class instance
    class(TDftbMultiExpan), intent(inout) :: this

    !> Multipole moments to pull
    type(TMultipole), intent(inout) :: multiExpanData

    integer, parameter :: xyzLen = 3
    integer :: nAtoms
    integer :: iAt1, iSp1, iSpin

    real(dp) :: tmpVx3(xyzLen), tmpM3x3(xyzLen, xyzLen)

    iSpin = 1
    nAtoms = this%nAtoms
    do iAt1 = 1, nAtoms
      iSp1 = this%species(iAt1)
      if ((this%onDQOCharges(1, iSp1) == 0) .and. (this%onDQOCharges(2, iSp1) == 0)) then
        cycle
      end if

      this%deltaDAtom(:,iAt1) = multiExpanData%dipoleAtom(:,iAt1,iSpin)
      ! this%deltaDAtom(1,iAt1) = multiExpanData%dipoleAtom(1,iAt1,iSpin)
      ! this%deltaDAtom(2,iAt1) = multiExpanData%dipoleAtom(2,iAt1,iSpin)
      ! this%deltaDAtom(3,iAt1) = multiExpanData%dipoleAtom(3,iAt1,iSpin)
      tmpM3x3(:,:) = 0.0_dp
      tmpM3x3(1,1) = multiExpanData%quadrupoleAtom(1,iAt1,iSpin)
      tmpM3x3(1,2) = multiExpanData%quadrupoleAtom(2,iAt1,iSpin)
      ! tmpM3x3(2,1) = tmpM3x3(1,2)
      tmpM3x3(2,2) = multiExpanData%quadrupoleAtom(3,iAt1,iSpin)
      tmpM3x3(1,3) = multiExpanData%quadrupoleAtom(4,iAt1,iSpin)
      ! tmpM3x3(3,1) = tmpM3x3(1,3)
      tmpM3x3(2,3) = multiExpanData%quadrupoleAtom(5,iAt1,iSpin)
      ! tmpM3x3(3,2) = tmpM3x3(2,3)
      tmpM3x3(3,3) = multiExpanData%quadrupoleAtom(6,iAt1,iSpin)
      call symmetrizeSquareMatrix(tmpM3x3)
      this%deltaQAtom(:,:,iAt1) = tmpM3x3(:,:)
    end do

  end subroutine pullDeltaDQAtom

  subroutine pushDeltaDQAtom(this, multiExpanData)

    !> class instance
    class(TDftbMultiExpan), intent(inout) :: this

    !> Multipole moments push
    type(TMultipole), intent(inout) :: multiExpanData

    integer, parameter :: xyzLen = 3
    integer :: nAtoms
    integer :: iAt1, iSp1, iSpin

    real(dp) :: tmpVx3(xyzLen), tmpM3x3(xyzLen, xyzLen)

    iSpin = 1
    nAtoms = this%nAtoms
    do iAt1 = 1, nAtoms
      iSp1 = this%species(iAt1)
      if ((this%onDQOCharges(1, iSp1) == 0) .and. (this%onDQOCharges(2, iSp1) == 0)) then
        cycle
      end if

      multiExpanData%dipoleAtom(:,iAt1,iSpin) = this%deltaDAtom(:,iAt1)
      ! multiExpanData%dipoleAtom(1,iAt1,iSpin) = this%deltaDAtom(1,iAt1)
      ! multiExpanData%dipoleAtom(2,iAt1,iSpin) = this%deltaDAtom(2,iAt1)
      ! multiExpanData%dipoleAtom(3,iAt1,iSpin) = this%deltaDAtom(3,iAt1)
      tmpM3x3(:,:) = this%deltaQAtom(:,:,iAt1)
      multiExpanData%quadrupoleAtom(1,iAt1,iSpin) = tmpM3x3(1,1) 
      multiExpanData%quadrupoleAtom(2,iAt1,iSpin) = tmpM3x3(1,2) 
      multiExpanData%quadrupoleAtom(3,iAt1,iSpin) = tmpM3x3(2,2) 
      multiExpanData%quadrupoleAtom(4,iAt1,iSpin) = tmpM3x3(1,3) 
      multiExpanData%quadrupoleAtom(5,iAt1,iSpin) = tmpM3x3(2,3) 
      multiExpanData%quadrupoleAtom(6,iAt1,iSpin) = tmpM3x3(3,3) 
    end do

  end subroutine pushDeltaDQAtom

  subroutine updateDQPotentials(this, deltaMAtom)

    !> class instance
    class(TDftbMultiExpan), intent(inout) :: this

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


  subroutine addMultiExpanHamiltonian(this, ham, over, nNeighbour, iNeighbour, &
      & species, orb, iPair, nAtom, img2CentCell)

    !> class instance
    class(TDftbMultiExpan), intent(inout) :: this

    !> The resulting Hamiltonian contribution.
    real(dp), intent(inout) :: ham(:,:)

    !> The overlap matrix
    real(dp), intent(in) :: over(:)

    !> Number of neighbours surrounding each atom.
    integer, intent(in) :: nNeighbour(:)

    !> List of neighbours for each atom.
    integer, intent(in) :: iNeighbour(0:,:)

    !> List of the species of each atom.
    integer, intent(in) :: species(:)

    !> Contains Information about the atomic orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> Indexing array for the Hamiltonian.
    integer, intent(in) :: iPair(0:,:)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Index mapping atoms onto the central cell atoms.
    integer, intent(in) :: img2CentCell(:)

    integer, parameter :: xyzLen = 3
    integer :: iAt1, iAt2, img, ind, nBlk, iBlk, iSp1, iSp2, iOrb1, iOrb2, iNeigh
    integer :: iOrig
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    integer :: mu, nu, kappa, mm, nn, kk
    real(dp) :: sqrTmp(orb%mOrb,orb%mOrb)
    real(dp) :: sqrTmpRho(orb%mOrb,orb%mOrb)
    real(dp) :: sqrTmpOver(orb%mOrb,orb%mOrb)
    real(dp) :: sqrTmpOverT(orb%mOrb,orb%mOrb)
    real(dp) :: sqrTmpHam(orb%mOrb,orb%mOrb)
    real(dp) :: mulTmp(orb%mOrb**2)
    real(dp) :: tmpVx3(xyzLen), tmpM3x3(xyzLen, xyzLen)
    real(dp) :: tmpadS(xyzLen), tmpSad(xyzLen)
    real(dp) :: tmpaQS(xyzLen, xyzLen), tmpSaQ(xyzLen, xyzLen)
    real(dp) :: tmpCoAtMu1(xyzLen, xyzLen), tmpCoAtNu1(xyzLen, xyzLen)
    real(dp) :: tmpCoAtMu2(xyzLen, xyzLen), tmpCoAtNu2(xyzLen, xyzLen)

    @:ASSERT(size(nNeighbour)==nAtom)
    @:ASSERT(size(iNeighbour,dim=2)==nAtom)
    @:ASSERT(size(species)>=maxval(iNeighbour))
    @:ASSERT(size(species)<=size(img2CentCell))
    @:ASSERT(size(iPair,dim=1)>=(maxval(nNeighbour)+1))
    @:ASSERT(size(iPair,dim=2)==nAtom)

    do iAtom1 = 1, nAtom
      iSp1 = species(iAtom1)
      nOrb1 = orb%nOrbAtom(iAtom1)
      do iNeigh = 0, nNeighbour(iAtom1)
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        iSp2 = this%species(iAtom2f)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        iOrig = iPair(iNeigh, iAtom1) + 1

        tmpVx3(:) = 0.0_dp
        tmpM3x3(:,:) = 0.0_dp
        sqrTmpOver(:,:) = 0.0_dp
        sqrTmpOverT(:,:) = 0.0_dp
        sqrTmpOver(1:nOrb2,1:nOrb1) = reshape(over(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
        sqrTmpOverT(1:nOrb1,1:nOrb2) = transpose(sqrTmpOver(1:nOrb2,1:nOrb1)) 

        sqrTmpHam(:,:) = 0.0_dp
        do mu = 1, nOrb1
          do nu =  1, nOrb2
            tmpadS(:) = 0.0_dp
            tmpSad(:) = 0.0_dp
            tmpaQS(:,:) = 0.0_dp
            tmpSaQ(:,:) = 0.0_dp
            do kk = 1, nOrb1
              tmpadS(:) = tmpadS(:) + this%atomicDIntgrl(:,kk,mu,iSp1) * sqrTmpOverT(kk,nu)
              tmpaQS(:,:) = tmpaQS(:,:) + this%atomicQIntgrl(:,:,kk,mu,iSp1) * sqrTmpOverT(kk,nu)
            end do
            do kk = 1, nOrb2
              tmpSad(:) = tmpSad(:) + this%atomicDIntgrl(:,kk,nu,iSp2) * sqrTmpOver(kk,mu)
              tmpSaQ(:,:) = tmpSaQ(:,:) + this%atomicQIntgrl(:,:,kk,nu,iSp2) * sqrTmpOver(kk,mu)
            end do

            ! add Monopole-Dipole contribution
            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) - (this%pot10x1Atom(iAtom1)&
                & + this%pot10x1Atom(iAtom2f)) * sqrTmpOver(nu,mu)

            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) + sum(this%pot10x0Atom(:,iAtom1) * tmpadS(:)) &
                                            & + sum(this%pot10x0Atom(:,iAtom2f) * tmpSad(:))

            !>add Dipole-Dipole contribution
            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) + sum(this%pot11x1Atom(:,iAtom1) * tmpadS(:)) &
                                            & + sum(this%pot11x1Atom(:,iAtom2f) * tmpSad(:))
 
            !>add Monopole-Quadrupole contribution
            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) + (this%pot20x2Atom(iAtom1)&
                & + this%pot20x2Atom(iAtom2f)) * sqrTmpOver(nu,mu)
            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) + sum(this%pot20x0Atom(:,:,iAtom1) * tmpaQS(:,:)) &
                                            & + sum(this%pot20x0Atom(:,:,iAtom2f) * tmpSaQ(:,:))

            !>add Dipole-Quadrupole contribution
            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) - sum(this%pot21x2Atom(:,iAtom1) * tmpadS(:)) &
                                            & - sum(this%pot21x2Atom(:,iAtom2f) * tmpSad(:))
            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) + sum(this%pot21x1Atom(:,:,iAtom1) * tmpaQS(:,:)) &
                                            & + sum(this%pot21x1Atom(:,:,iAtom2f) * tmpSaQ(:,:))

            !>add Quadrupole-Quadrupole contribution
            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) + sum(this%pot22x2Atom(:,:,iAtom1) * tmpaQS(:,:)) &
                                            & + sum(this%pot22x2Atom(:,:,iAtom2f) * tmpSaQ(:,:))
          end do
        end do

        ham(iOrig:iOrig+nOrb1*nOrb2-1, 1) = ham(iOrig:iOrig+nOrb1*nOrb2-1, 1)&
           & + 0.5_dp * reshape(sqrTmpHam(1:nOrb2,1:nOrb1), (/nOrb2*nOrb1/))
      end do
    end do

  end subroutine addMultiExpanHamiltonian

  !> Add the MultiExpan-Energy contribution to the total energy
  subroutine addMultiExpanEnergy(this, energyPerAtom, energyMD, energyDD, energyMQ, energyDQ,&
      & energyQQ, energyTT)

    !> RangeSep class instance
    class(TDftbMultiExpan), intent(inout) :: this

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

  end subroutine addMultiExpanEnergy


  !> Returns energy per atom.
  subroutine getEnergyPerAtom(this, energyPerAtom)

    class(TDftbMultiExpan), intent(inout) :: this

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

  subroutine addMultiExpanGradients(this, derivs, derivator, skOverCont, rho, species,&
      & iNeighbour, nNeighbourSK, img2CentCell, iPair, orb, coords)

    !> class instance
    class(TDftbMultiExpan), intent(inout) :: this

    !> Gradient on exit.
    real(dp), intent(inout) :: derivs(:,:)

    !> Differentiatior for the non-scc components
    class(TNonSccDiff), intent(in) :: derivator

    !> Container for SK overlap integrals
    type(TSlakoCont) :: skOverCont

    !> Density matrix in Packed format
    real(dp), intent(in) :: rho(:)

    !> Specie for each atom.
    integer, intent(in) :: species(:)

    !> Neighbour list for atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours of each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> Mapping of atoms to cetnral cell.
    integer, intent(in) :: img2CentCell(:)

    !> Indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Information about the shells and orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Coordinate of each atom.
    real(dp), intent(in) :: coords(:,:)

    integer, parameter :: xyzLen = 3
    integer :: iAt1, iAt2, iAt2f
    integer :: iOrig, iSpin, nSpin, nAtom, iNeigh, iSp1, iSp2
    integer :: nOrb1, nOrb2
    integer :: ii, jj, ll, mu, nu, kappa, mm, nn, kk

    real(dp) :: sqrDMTmp(orb%mOrb,orb%mOrb), sqrDMTmpT(orb%mOrb,orb%mOrb)
    real(dp) :: sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: derivTmp(3)
    real(dp) :: rab, vRabx3(3), tmp, tmp3(3)

    real(dp) :: sPrimeTmp12(orb%mOrb,orb%mOrb,3)
    real(dp) :: tmpDerivMuNu, tmpDeriv(xyzLen), tmpVx3(xyzLen), tmpM3x3(xyzLen, xyzLen)


    ! real(dp) :: sqrTmp(orb%mOrb,orb%mOrb)
    ! real(dp) :: sqrTmpRho(orb%mOrb,orb%mOrb)
    ! real(dp) :: sqrTmpOver(orb%mOrb,orb%mOrb)
    ! real(dp) :: sqrTmpOverT(orb%mOrb,orb%mOrb)
    ! real(dp) :: sqrTmpHam(orb%mOrb,orb%mOrb)
    ! real(dp) :: mulTmp(orb%mOrb**2)
    real(dp) :: tmpPad(xyzLen), tmpadP(xyzLen)
    real(dp) :: tmpPaQ(xyzLen, xyzLen), tmpaQP(xyzLen, xyzLen)

    !> Calculate pot30x0, pot31x1, and pot32x2
    this%pot30x0Atom(:,:,:,:) = 0.0_dp
    this%pot31x1Atom(:,:,:,:) = 0.0_dp
    this%pot32x2Atom(:,:,:,:) = 0.0_dp

    !$OMP PARALLEL DO PRIVATE(iAt1, iSp1, iAt2, iSp2, ii, jj, ll) DEFAULT(SHARED) SCHEDULE(RUNTIME)
    do iAt1 = 1, this%nAtoms
      iSp1 = species(iAt1)
      do iAt2 = 1, this%nAtoms
        iSp2 = species(iAt2)
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

    do iAt1 = 1, this%nAtoms
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

      iSp1 = species(iAt1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      do iNeigh = 1, nNeighbourSK(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        if (iAt1 == iAt2f) then
          cycle
        end if
        nOrb2 = orb%nOrbSpecies(iSp2)
        iOrig = iPair(iNeigh,iAt1) + 1
        sqrDMTmp(1:nOrb2,1:nOrb1) = reshape(rho(iOrig:iOrig+nOrb1*nOrb2-1), [nOrb2,nOrb1])
        sqrDMTmpT(1:nOrb1,1:nOrb2) = transpose(sqrDMTmp(1:nOrb2,1:nOrb1)) 
        call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, species, iAt1, iAt2, orb)

        tmpDeriv(:) = 0.0_dp
        derivTmp(:) = 0.0_dp
        ! note factor of 2 for implicit summation over lower triangle of density matrix:
        ! do ii = 1, 3
        !   derivTmp(ii) = 2.0_dp * (&
        !       & - sum(sqrEDMTmp(1:nOrb2,1:nOrb1)*sPrimeTmp(1:nOrb2,1:nOrb1,ii)))
        ! end do

        ! iSpin = 1
        ! do ii = 1, 3
        !   shiftSprime(1:nOrb2,1:nOrb1) = 0.5_dp * (&
        !       & matmul(sPrimeTmp(1:nOrb2,1:nOrb1,ii), shift(1:nOrb1,1:nOrb1,iAt1,iSpin) )&
        !       & + matmul(shift(1:nOrb2,1:nOrb2,iAt2f,iSpin), sPrimeTmp(1:nOrb2,1:nOrb1,ii)))
        !   ! again factor of 2 from lower triangle, cf published force expressions for SCC:
        !   derivTmp(ii) = derivTmp(ii) + 2.0_dp * ( sum(shiftSprime(1:nOrb2,1:nOrb1) *&
        !       & reshape(rho(iOrig:iOrig+nOrb1*nOrb2-1), [nOrb2,nOrb1]) ) )
        ! end do

        ! vRabx3(:) = coords(:,iAt1) - coords(:,iAt2)
        ! rab = norm2(vRabx3(:))
        derivTmp(:) = 0.0_dp
        do mu = 1, nOrb1
          do nu =  1, nOrb2

            tmpPad(:) = 0.0_dp
            tmpadP(:) = 0.0_dp
            tmpPaQ(:,:) = 0.0_dp
            tmpaQP(:,:) = 0.0_dp
            do kk = 1, nOrb1
              tmpPad(:) = tmpPad(:) + this%atomicDIntgrl(:,kk,mu,iSp1) * sqrDMTmpT(kk,nu)
              tmpPaQ(:,:) = tmpPaQ(:,:) + this%atomicQIntgrl(:,:,kk,mu,iSp1) * sqrDMTmpT(kk,nu)
            end do
            do kk = 1, nOrb2
              tmpadP(:) = tmpadP(:) + this%atomicDIntgrl(:,kk,nu,iSp2) * sqrDMTmp(kk,mu)
              tmpaQP(:,:) = tmpaQP(:,:) + this%atomicQIntgrl(:,:,kk,nu,iSp2) * sqrDMTmp(kk,mu)
            end do

            tmpDerivMuNu = 0.0_dp
            !> M-D contribution
            tmpDerivMuNu = tmpDerivMuNu - (this%pot10x1Atom(iAt1) + this%pot10x1Atom(iAt2))&
                & * sqrDMTmp(nu,mu)
            !> M-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + (this%pot20x2Atom(iAt1) + this%pot20x2Atom(iAt2))&
                & * sqrDMTmp(nu,mu)

            !> M-D contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot10x0Atom(:,iAt1) * tmpPad(:))
            !> D-D contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot11x1Atom(:,iAt1) * tmpPad(:))
            !> D-Q contribution
            tmpDerivMuNu = tmpDerivMuNu - sum(this%pot21x2Atom(:,iAt1) * tmpPad(:))
            !> M-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot20x0Atom(:,:,iAt1) * tmpPaQ(:,:))
            !> D-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot21x1Atom(:,:,iAt1) * tmpPaQ(:,:))
            !> Q-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot22x2Atom(:,:,iAt1) * tmpPaQ(:,:))
            
            !> M-D contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot10x0Atom(:,iAt2) * tmpadP(:))
            !> D-D contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot11x1Atom(:,iAt2) * tmpadP(:))
            !> D-Q contribution
            tmpDerivMuNu = tmpDerivMuNu - sum(this%pot21x2Atom(:,iAt2) * tmpadP(:))
            !> M-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot20x0Atom(:,:,iAt2) * tmpaQP(:,:))
            !> D-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot21x1Atom(:,:,iAt2) * tmpaQP(:,:))
            !> Q-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot22x2Atom(:,:,iAt2) * tmpaQP(:,:))

            ! tmpDeriv(:) = tmpDeriv(:) + tmpDerivMuNu * sPrimeTmp12(nn,mm,:)
            derivTmp(:) = derivTmp(:) + tmpDerivMuNu * sPrimeTmp(nu,mu,:)

          end do
        end do
        ! forces from atom 1 on atom 2f and 2f onto 1
        derivs(:,iAt1) = derivs(:,iAt1) + derivTmp
        derivs(:,iAt2f) = derivs(:,iAt2f) - derivTmp
      end do
    end do
  end subroutine addMultiExpanGradients


  subroutine addAtomicDipoleMoment(this, dipoleMoment)

    !> class instance
    class(TDftbMultiExpan), intent(inout) :: this

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
    class(TDftbMultiExpan), intent(inout) :: this

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

end module dftbp_dftb_multiexpan
