!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Routines implementing the (One-center approximation) multipole expansion for the 2nd order DFTB.
module dftbp_dftb_mdftb
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_globalenv, only : stdOut
  use dftbp_dftb_nonscc, only : TNonSccDiff
  use dftbp_dftb_periodic, only : TNeighbourList, getNrOfNeighbours
  use dftbp_dftb_shortgammafuncs, only : expGammaPrime, expGammaDoublePrime, expGammaTriplePrime,&
      & expGammaQuadruplePrime, expGammaQuintuplePrime
  use dftbp_dftb_slakocont, only : TSlakoCont
  use dftbp_io_message, only : error
  use dftbp_math_matrixops, only : adjointLowerTriangle
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_multipole, only : TMultipole
  implicit none
  private

  public :: TMdftbAtomicIntegrals, TMdftbInp, TMdftb, TMdftb_init

  !> Store one-center atomic integrals required by mdftb
  type TMdftbAtomicIntegrals
    real(dp), allocatable :: DScaling(:)
    real(dp), allocatable :: QScaling(:)
    real(dp), allocatable :: SXPx(:)
    real(dp), allocatable :: PxXDxxyy(:)
    real(dp), allocatable :: PxXDzz(:)
    real(dp), allocatable :: PyYDxxyy(:)
    real(dp), allocatable :: PzZDzz(:)
    real(dp), allocatable :: SXXS(:)
    real(dp), allocatable :: PxXXPx(:)
    real(dp), allocatable :: PyXXPy(:)
    real(dp), allocatable :: SXXDxxyy(:)
    real(dp), allocatable :: SXXDzz(:)
    real(dp), allocatable :: SYYDxxyy(:)
    real(dp), allocatable :: SZZDzz(:)
    real(dp), allocatable :: DxyXXDxy(:)
    real(dp), allocatable :: DyzXXDyz(:)
    real(dp), allocatable :: DxxyyXXDzz(:)
    real(dp), allocatable :: DzzXXDzz(:)
    real(dp), allocatable :: DxxyyYYDzz(:)
    real(dp), allocatable :: DzzZZDzz(:)
    real(dp), allocatable :: DxzXZDzz(:)
    real(dp), allocatable :: DyzYZDxxyy(:)
  end type TMdftbAtomicIntegrals


  !> Input for the MultiExpan module
  type TMdftbInp

    !> Orbital information
    integer :: nOrb = 0, nSpin = 0, nSpecies = 0
    type(TOrbitals), pointer :: orb => null()

    !> Hubbard U values for atoms
    real(dp), allocatable :: hubbu(:)

    !> Species of atoms
    integer, allocatable :: species(:)

    !> One-center atomic integral data
    type(TMdftbAtomicIntegrals), allocatable :: mdftbAtomicIntegrals

  end type TMdftbInp


  !> Internal management for the TMultiExpan.
  type TMdftb
    integer :: nAtoms = 0, nSpecies = 0, mOrb = 0, nOrb = 0, nSpin = 0

    !> Hubbard U values for atoms
    real(dp), allocatable :: hubbu(:)

    !> Species of atoms
    integer, allocatable :: species(:), nOrbSpecies(:)
    !> Whether a species has dipole or quadrupole on-site charges (2, nSpecies)
    logical, allocatable :: hasOnsiteCharges(:,:)

    !> Atomic dipole integrals
    real(dp), allocatable :: atomicDIntgrl(:,:,:,:)
    !> Atomic quadrupole integrals
    real(dp), allocatable :: atomicQIntgrl(:,:,:,:,:)

    !> Mulliken charge per atom
    real(dp), allocatable :: deltaMAtom(:)
    !> Dipole charge per atom
    real(dp), allocatable :: deltaDAtom(:,:)
    !> Quadrupole charge per atom
    real(dp), allocatable :: deltaQAtom(:,:,:)

    !> Evaluated for the E, Atom1, Atom2 at each geometry step
    real(dp), allocatable :: f10AB(:,:,:)
    real(dp), allocatable :: f20AB(:,:,:,:)
    real(dp), allocatable :: f30AB(:,:,:,:,:)
    real(dp), allocatable :: f40AB(:,:,:,:,:,:)

    !> Evaluated for the gradient, Atom1, Atom2 at each geometry step
    real(dp), allocatable :: f50AB(:,:,:,:,:,:,:)

    !> Add for the H and the E
    real(dp), allocatable :: pot10x1Atom(:)
    real(dp), allocatable :: pot20x2Atom(:)
    real(dp), allocatable :: pot10x0Atom(:,:)
    real(dp), allocatable :: pot11x1Atom(:,:)
    real(dp), allocatable :: pot21x2Atom(:,:)
    real(dp), allocatable :: pot20x0Atom(:,:,:)
    real(dp), allocatable :: pot21x1Atom(:,:,:)
    real(dp), allocatable :: pot22x2Atom(:,:,:)

    !> Add for the gradient
    real(dp), allocatable :: pot30x0Atom(:,:,:,:)
    real(dp), allocatable :: pot31x1Atom(:,:,:,:)
    real(dp), allocatable :: pot32x2Atom(:,:,:,:)

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
  end type TMdftb


  !> Number of dipole components (x, y, z)
  integer, parameter :: dimDipole = 3

  !> Number of quadrupole components (xx, xy, yy, xz, yz, zz)
  integer, parameter :: dimQuadrupole = 6


contains

  !> Get information on required multipolar contributions
  subroutine getMultiExpanInfo(this, nDipole, nQuadrupole)

    !> Data structure
    class(TMdftb), intent(in) :: this

    !> Number of dipole moment components
    integer, intent(out) :: nDipole

    !> Number of quadrupole moment components
    integer, intent(out) :: nQuadrupole

    nDipole = dimDipole
    nQuadrupole = dimQuadrupole

  end subroutine getMultiExpanInfo


  !> Returns the equivalence to get the correct mixing of charge dependent contributions
  subroutine getOrbitalEquiv(this, equivDip, equivQuad)

    !> Data structure
    class(TMdftb), intent(inout) :: this

    !> The equivalence vector for cumulative atomic dipole populations
    integer, intent(out) :: equivDip(:,:)

    !> The equivalence vector for cumulative atomic quadrupole populations
    integer, intent(out) :: equivQuad(:,:)

    integer :: nAtom, ii

    nAtom = size(equivDip, dim=2)

    equivDip(:,:) = 0
    equivQuad(:,:) = 0
    equivDip(:,:) = reshape([(ii, ii = 1, dimDipole * nAtom)], [dimDipole, nAtom])
    equivQuad(:,:) = reshape([(ii, ii = 1, dimQuadrupole * nAtom)], [dimQuadrupole, nAtom])

  end subroutine getOrbitalEquiv


  !> Initializes instance.
  subroutine TMdftb_init(this, inp)

    !> Instance.
    type(TMdftb), intent(out) :: this

    !> Input data.
    type(TMdftbInp), intent(in) :: inp

    integer, parameter :: icx = 1, icy = 2, icz = 3
    integer, parameter :: ios = 0
    integer, parameter :: iopy = 0, iopz = 1, iopx = 2
    integer, parameter :: iodxy = 0, iodyz = 1, iodzz = 2, iodxz = 3, iodxxyy = 4
    real(dp), parameter :: minIntgrl = 1.0e-9_dp
    integer :: nAtoms, mOrb, nSpecies, nOrb1, iSp1, ii, jj, mm, nn
    integer :: ang1, ang2, iSh1, iSh2, iOrbAng1, iOrbAng2
    real(dp) tmpIntgrl, tmpAvgTrace

    this%nAtoms = size(inp%orb%nOrbAtom)
    this%nSpecies = inp%nSpecies
    this%nOrb = inp%nOrb
    this%mOrb = inp%orb%mOrb
    this%nSpin = inp%nSpin
    this%hubbu = inp%hubbu
    this%species = inp%species
    this%nOrbSpecies = inp%orb%nOrbSpecies

    nAtoms = this%nAtoms
    mOrb = this%mOrb
    nSpecies = this%nSpecies

    allocate(this%atomicDIntgrl(3, mOrb, mOrb, nSpecies), source=0.0_dp)
    allocate(this%atomicQIntgrl(3, 3, mOrb, mOrb, nSpecies), source=0.0_dp)
    
    if (maxval(inp%orb%angShell) >= 3) then
        call error("DFTB multipole expansion currently unsupported for chemical elements&
            & having angular moments higher than 2")
    end if

    do iSp1 = 1, nSpecies
      do iSh1 = 1, inp%orb%nShell(iSp1)
        ang1 = inp%orb%angShell(iSh1, iSp1)
        iOrbAng1 = inp%orb%posShell(iSh1, iSp1)
        do iSh2 = 1, inp%orb%nShell(iSp1)
          ang2 = inp%orb%angShell(iSh2, iSp1)
          iOrbAng2 = inp%orb%posShell(iSh2, iSp1)

          if (ang1 == 0 .and. ang2 == 0) then
            ! Assign Quadrupole <S|XX|S>
            tmpIntgrl = inp%mdftbAtomicIntegrals%SXXS(iSp1)
            this%atomicQIntgrl(icx, icx, iOrbAng1+ios, iOrbAng2+ios, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icy, icy, iOrbAng1+ios, iOrbAng2+ios, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icz, icz, iOrbAng1+ios, iOrbAng2+ios, iSp1) = tmpIntgrl

          else if (ang1 == 0 .and. ang2 == 1) then
            ! Assign Dipole <S|X|Px>
            tmpIntgrl = inp%mdftbAtomicIntegrals%SXPx(iSp1)
            this%atomicDIntgrl(icx, iOrbAng1+ios, iOrbAng2+iopx, iSp1) = tmpIntgrl
            this%atomicDIntgrl(icy, iOrbAng1+ios, iOrbAng2+iopy, iSp1) = tmpIntgrl
            this%atomicDIntgrl(icz, iOrbAng1+ios, iOrbAng2+iopz, iSp1) = tmpIntgrl

          else if (ang1 == 0 .and. ang2 == 2) then
            ! Assign Quadrupole <S|XX|Dxx-yy>
            tmpIntgrl = inp%mdftbAtomicIntegrals%SXXDxxyy(iSp1)
            this%atomicQIntgrl(icx, icx, iOrbAng1+ios, iOrbAng2+iodxxyy, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icx, icy, iOrbAng1+ios, iOrbAng2+iodxy, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icx, icz, iOrbAng1+ios, iOrbAng2+iodxz, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icy, icz, iOrbAng1+ios, iOrbAng2+iodyz, iSp1) = tmpIntgrl

            ! Assign Quadrupole <S|XX|Dzz>
            tmpIntgrl = inp%mdftbAtomicIntegrals%SXXDzz(iSp1)
            this%atomicQIntgrl(icx, icx, iOrbAng1+ios, iOrbAng2+iodzz, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icy, icy, iOrbAng1+ios, iOrbAng2+iodzz, iSp1) = tmpIntgrl

            ! Assign Quadrupole <S|YY|Dxx-yy>
            tmpIntgrl = inp%mdftbAtomicIntegrals%SYYDxxyy(iSp1)
            this%atomicQIntgrl(icy, icy, iOrbAng1+ios, iOrbAng2+iodxxyy, iSp1) = tmpIntgrl

            ! Assign Quadrupole <S|ZZ|Dzz>
            tmpIntgrl = inp%mdftbAtomicIntegrals%SZZDzz(iSp1)
            this%atomicQIntgrl(icz, icz, iOrbAng1+ios, iOrbAng2+iodzz, iSp1) = tmpIntgrl

          else if (ang1 == 1 .and. ang2 == 1) then
            ! Assign Quadrupole <Px|XX|Px>
            tmpIntgrl = inp%mdftbAtomicIntegrals%PxXXPx(iSp1)
            this%atomicQIntgrl(icx, icx, iOrbAng1+iopx, iOrbAng2+iopx, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icy, icy, iOrbAng1+iopy, iOrbAng2+iopy, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icz, icz, iOrbAng1+iopz, iOrbAng2+iopz, iSp1) = tmpIntgrl

            ! Assign Quadrupole <Py|XX|Py>
            tmpIntgrl = inp%mdftbAtomicIntegrals%PyXXPy(iSp1)
            this%atomicQIntgrl(icx, icx, iOrbAng1+iopy, iOrbAng2+iopy, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icx, icx, iOrbAng1+iopz, iOrbAng2+iopz, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icy, icy, iOrbAng1+iopx, iOrbAng2+iopx, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icy, icy, iOrbAng1+iopz, iOrbAng2+iopz, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icz, icz, iOrbAng1+iopx, iOrbAng2+iopx, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icz, icz, iOrbAng1+iopy, iOrbAng2+iopy, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icx, icy, iOrbAng1+iopx, iOrbAng2+iopy, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icx, icz, iOrbAng1+iopx, iOrbAng2+iopz, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icy, icz, iOrbAng1+iopy, iOrbAng2+iopz, iSp1) = tmpIntgrl

          else if (ang1 == 1 .and. ang2 == 2) then
            ! Assign Dipole <Px|X|Dxx-yy>
            tmpIntgrl = inp%mdftbAtomicIntegrals%PxXDxxyy(iSp1)
            this%atomicDIntgrl(icx, iOrbAng1+iopx, iOrbAng2+iodxxyy, iSp1) = tmpIntgrl
            this%atomicDIntgrl(icx, iOrbAng1+iopy, iOrbAng2+iodxy, iSp1) = tmpIntgrl
            this%atomicDIntgrl(icx, iOrbAng1+iopz, iOrbAng2+iodxz, iSp1) = tmpIntgrl
            this%atomicDIntgrl(icy, iOrbAng1+iopx, iOrbAng2+iodxy, iSp1) = tmpIntgrl
            this%atomicDIntgrl(icy, iOrbAng1+iopz, iOrbAng2+iodyz, iSp1) = tmpIntgrl
            this%atomicDIntgrl(icz, iOrbAng1+iopx, iOrbAng2+iodxz, iSp1) = tmpIntgrl
            this%atomicDIntgrl(icz, iOrbAng1+iopy, iOrbAng2+iodyz, iSp1) = tmpIntgrl

            ! Assign Dipole <Px|X|Dzz>
            tmpIntgrl = inp%mdftbAtomicIntegrals%PxXDzz(iSp1)
            this%atomicDIntgrl(icx, iOrbAng1+iopx, iOrbAng2+iodzz, iSp1) = tmpIntgrl
            this%atomicDIntgrl(icy, iOrbAng1+iopy, iOrbAng2+iodzz, iSp1) = tmpIntgrl

            ! Assign Dipole <Py|Y|Dxx-yy>
            tmpIntgrl = inp%mdftbAtomicIntegrals%PyYDxxyy(iSp1)
            this%atomicDIntgrl(icy, iOrbAng1+iopy, iOrbAng2+iodxxyy, iSp1) = tmpIntgrl

            ! Assign Dipole <Pz|Z|Dzz>
            tmpIntgrl = inp%mdftbAtomicIntegrals%PzZDzz(iSp1)
            this%atomicDIntgrl(icz, iOrbAng1+iopz, iOrbAng2+iodzz, iSp1) = tmpIntgrl

          else if (ang1 == 2 .and. ang2 == 2) then
            ! Assign Quadrupole <Dxy|XX|Dxy>
            tmpIntgrl = inp%mdftbAtomicIntegrals%DxyXXDxy(iSp1)
            this%atomicQIntgrl(icx, icx, iOrbAng1+iodxy, iOrbAng2+iodxy, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icx, icx, iOrbAng1+iodxz, iOrbAng2+iodxz, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icx, icx, iOrbAng1+iodxxyy, iOrbAng2+iodxxyy, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icy, icy, iOrbAng1+iodxy, iOrbAng2+iodxy, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icy, icy, iOrbAng1+iodyz, iOrbAng2+iodyz, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icy, icy, iOrbAng1+iodxxyy, iOrbAng2+iodxxyy, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icz, icz, iOrbAng1+iodxz, iOrbAng2+iodxz, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icz, icz, iOrbAng1+iodyz, iOrbAng2+iodyz, iSp1) = tmpIntgrl

            ! Assign Quadrupole <Dyz|XX|Dyz>
            tmpIntgrl = inp%mdftbAtomicIntegrals%DyzXXDyz(iSp1)
            this%atomicQIntgrl(icx, icx, iOrbAng1+iodyz, iOrbAng2+iodyz, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icy, icy, iOrbAng1+iodxz, iOrbAng2+iodxz, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icz, icz, iOrbAng1+iodxy, iOrbAng2+iodxy, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icz, icz, iOrbAng1+iodxxyy, iOrbAng2+iodxxyy, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icx, icy, iOrbAng1+iodxz, iOrbAng2+iodyz, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icx, icz, iOrbAng1+iodxy, iOrbAng2+iodyz, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icx, icz, iOrbAng1+iodxz, iOrbAng2+iodxxyy, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icy, icz, iOrbAng1+iodxy, iOrbAng2+iodxz, iSp1) = tmpIntgrl

            ! Assign Quadrupole <Dxx-yy|XX|Dzz>
            tmpIntgrl = inp%mdftbAtomicIntegrals%DxxyyXXDzz(iSp1)
            this%atomicQIntgrl(icx, icx, iOrbAng1+iodxxyy, iOrbAng2+iodzz, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icx, icy, iOrbAng1+iodxy, iOrbAng2+iodzz, iSp1) = tmpIntgrl

            ! Assign Quadrupole <Dzz|XX|Dzz>
            tmpIntgrl = inp%mdftbAtomicIntegrals%DzzXXDzz(iSp1)
            this%atomicQIntgrl(icx, icx, iOrbAng1+iodzz, iOrbAng2+iodzz, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icy, icy, iOrbAng1+iodzz, iOrbAng2+iodzz, iSp1) = tmpIntgrl

            ! Assign Quadrupole <Dxx-yy|YY|Dzz>
            tmpIntgrl = inp%mdftbAtomicIntegrals%DxxyyYYDzz(iSp1)
            this%atomicQIntgrl(icy, icy, iOrbAng1+iodxxyy, iOrbAng2+iodzz, iSp1) = tmpIntgrl

            ! Assign Quadrupole <Dzz|ZZ|Dzz>
            tmpIntgrl = inp%mdftbAtomicIntegrals%DzzZZDzz(iSp1)
            this%atomicQIntgrl(icz, icz, iOrbAng1+iodzz, iOrbAng2+iodzz, iSp1) = tmpIntgrl

            ! Assign Quadrupole <Dxz|XZ|Dzz>
            tmpIntgrl = inp%mdftbAtomicIntegrals%DxzXZDzz(iSp1)
            this%atomicQIntgrl(icx, icz, iOrbAng1+iodxz, iOrbAng2+iodzz, iSp1) = tmpIntgrl
            this%atomicQIntgrl(icy, icz, iOrbAng1+iodyz, iOrbAng2+iodzz, iSp1) = tmpIntgrl

            ! Assign Quadrupole <Dyz|YZ|Dxx-yy>
            tmpIntgrl = inp%mdftbAtomicIntegrals%DyzYZDxxyy(iSp1)
            this%atomicQIntgrl(icy, icz, iOrbAng1+iodyz, iOrbAng2+iodxxyy, iSp1) = tmpIntgrl
          end if
        end do
      end do
    end do

    do iSp1 = 1, nSpecies
      nOrb1 = this%nOrbSpecies(iSp1)
      do mm = 1, nOrb1
        do nn = 1, nOrb1
          do ii = 1, 3
            if (abs(this%atomicDIntgrl(ii,mm,nn,iSp1)) >= minIntgrl) then
              this%atomicDIntgrl(ii,nn,mm,iSp1) = this%atomicDIntgrl(ii,mm,nn,iSp1)
            else
              this%atomicDIntgrl(ii,mm,nn,iSp1) = 0.0_dp
            end if
            do jj = 1, 3
              if (abs(this%atomicQIntgrl(ii,jj,mm,nn,iSp1)) >= minIntgrl) then
                this%atomicQIntgrl(ii,jj,nn,mm,iSp1) = this%atomicQIntgrl(ii,jj,mm,nn,iSp1)
                this%atomicQIntgrl(jj,ii,mm,nn,iSp1) = this%atomicQIntgrl(ii,jj,mm,nn,iSp1)
                this%atomicQIntgrl(jj,ii,nn,mm,iSp1) = this%atomicQIntgrl(ii,jj,mm,nn,iSp1)
              else
                this%atomicQIntgrl(ii,jj,mm,nn,iSp1) = 0.0_dp
              end if
            end do
          end do
        end do
      end do
    end do

   ! Scale atomic integrals
    do iSp1 = 1, nSpecies
      this%atomicDIntgrl(:,:,:,iSp1) = inp%mdftbAtomicIntegrals%DScaling(iSp1)&
          & * this%atomicDIntgrl(:,:,:,iSp1)
      this%atomicQIntgrl(:,:,:,:,iSp1) = inp%mdftbAtomicIntegrals%QScaling(iSp1)&
          & * this%atomicQIntgrl(:,:,:,:,iSp1)
    end do

   ! Remove Trace of atomic quadrupole integrals
    do iSp1 = 1, nSpecies
      nOrb1 = this%nOrbSpecies(iSp1)
      do mm = 1, nOrb1
        do nn = 1, nOrb1
          tmpAvgTrace = (this%atomicQIntgrl(1,1,mm,nn,iSp1) + this%atomicQIntgrl(2,2,mm,nn,iSp1)&
              & + this%atomicQIntgrl(3,3,mm,nn,iSp1)) / 3.0_dp
          do ii = 1, 3
            this%atomicQIntgrl(ii,ii,mm,nn,iSp1) = this%atomicQIntgrl(ii,ii,mm,nn,iSp1)&
              & - tmpAvgTrace
          end do
        end do
      end do
    end do

    allocate(this%hasOnsiteCharges(2, nSpecies), source=.true.)
    do iSp1 = 1, nSpecies
      this%hasOnsiteCharges(1, iSp1) = (maxval(abs(this%atomicDIntgrl(:,:,:,iSp1))) >= minIntgrl)
      this%hasOnsiteCharges(2, iSp1) = (maxval(abs(this%atomicQIntgrl(:,:,:,:,iSp1))) >= minIntgrl)
    end do

    allocate(this%deltaMAtom(nAtoms), source=0.0_dp)
    allocate(this%deltaDAtom(3, nAtoms), source=0.0_dp)
    allocate(this%deltaQAtom(3, 3, nAtoms), source=0.0_dp)
    allocate(this%f10AB(3, nAtoms, nAtoms), source=0.0_dp)
    allocate(this%f20AB(3, 3, nAtoms, nAtoms), source=0.0_dp)
    allocate(this%f30AB(3, 3, 3, nAtoms, nAtoms), source=0.0_dp)
    allocate(this%f40AB(3, 3, 3, 3, nAtoms, nAtoms), source=0.0_dp)
    allocate(this%f50AB(3, 3, 3, 3, 3, nAtoms, nAtoms), source=0.0_dp)

    allocate(this%pot10x1Atom(nAtoms), source=0.0_dp)
    allocate(this%pot20x2Atom(nAtoms), source=0.0_dp)
    allocate(this%pot10x0Atom(3, nAtoms), source=0.0_dp)
    allocate(this%pot11x1Atom(3, nAtoms), source=0.0_dp)
    allocate(this%pot21x2Atom(3, nAtoms), source=0.0_dp)
    allocate(this%pot20x0Atom(3, 3, nAtoms), source=0.0_dp)
    allocate(this%pot21x1Atom(3, 3, nAtoms), source=0.0_dp)
    allocate(this%pot22x2Atom(3, 3, nAtoms), source=0.0_dp)
    allocate(this%pot30x0Atom(3, 3, 3, nAtoms), source=0.0_dp)
    allocate(this%pot31x1Atom(3, 3, 3, nAtoms), source=0.0_dp)
    allocate(this%pot32x2Atom(3, 3, 3, nAtoms), source=0.0_dp)

  end subroutine TMdftb_init


  !> Updates data structures if there are changed coordinates for the instance.
  subroutine updateCoords(this, coords)

    !> Class instance
    class(TMdftb), intent(inout) :: this

    !> List of atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    integer :: nAtoms, iAt1, iAt2, iSp1, iSp2, ii, jj, ll, mm, nn
    real(dp) :: rab, u1, u2, coeffTerm1, coeffTerm2, coeffTerm3
    real(dp) :: gammaPrime, gammaDoublePrime, gammaTriplePrime, gammaQuadruplePrime,&
        & gammaQuintuplePrime

    real(dp) :: vRabx3(3), workVx3(3)
    real(dp) :: mI3x3(3,3), mRRab3x3(3,3), workM3x3(3,3), workM3x3x3(3,3,3)
    real(dp) :: mRoIab3x3x3(3,3,3), mIotRab3x3x3(3,3,3), mIoRab3x3x3(3,3,3)
    real(dp) :: mRIotoRab3x3x3(3,3,3), mRRRab3x3x3(3,3,3)

    real(dp) :: mI3x3otI3x3(3,3,3,3), mI3x3ohI3x3(3,3,3,3)
    real(dp) :: mI3x3oI3x3(3,3,3,3), mI3x3otohoI3x3(3,3,3,3)
    real(dp) :: workM3x3x3x3(3,3,3,3), f40Term1M3x3x3x3(3,3,3,3)
    real(dp) :: f40Term3M3x3x3x3(3,3,3,3), f40Term2M3x3x3x3(3,3,3,3)
    real(dp) :: workM3x3x3x3x3(3,3,3,3,3), f50Term1M3x3x3x3x3(3,3,3,3,3)
    real(dp) :: f50Term3M3x3x3x3x3(3,3,3,3,3), f50Term2M3x3x3x3x3(3,3,3,3,3)

    this%f10AB(:,:,:) = 0.0_dp
    this%f20AB(:,:,:,:) = 0.0_dp
    this%f30AB(:,:,:,:,:) = 0.0_dp
    this%f40AB(:,:,:,:,:,:) = 0.0_dp
    this%f50AB(:,:,:,:,:,:,:) = 0.0_dp

    mI3x3(:,:) = 0.0_dp
    do ii = 1, 3
      mI3x3(ii, ii) = 1.0_dp
    end do
    call outerProductA3x3OB3x3(mI3x3oI3x3, mI3x3, mI3x3)
    call outerProductA3x3OhB3x3(mI3x3ohI3x3, mI3x3, mI3x3)
    call outerProductA3x3OtB3x3(mI3x3otI3x3, mI3x3, mI3x3)
    mI3x3otohoI3x3(:,:,:,:) = mI3x3otI3x3 + mI3x3ohI3x3 + mI3x3oI3x3

    nAtoms = this%nAtoms

    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP PRIVATE(ii, jj, mm, ll, iAt1, iSp1, u1, iAt2, iSp2, u2, rab, vRabx3, workVx3) &
    !$OMP PRIVATE(gammaPrime, gammaDoublePrime, gammaTriplePrime, gammaQuadruplePrime) &
    !$OMP PRIVATE(gammaQuintuplePrime, coeffTerm1, coeffTerm2, coeffTerm3) &
    !$OMP PRIVATE(workM3x3, workM3x3x3, mRRab3x3) &
    !$OMP PRIVATE(mRoIab3x3x3, mIotRab3x3x3, mIoRab3x3x3, mRIotoRab3x3x3, mRRRab3x3x3) &
    !$OMP PRIVATE(workM3x3x3x3, f40Term1M3x3x3x3, f40Term2M3x3x3x3, f40Term3M3x3x3x3) &
    !$OMP PRIVATE(workM3x3x3x3x3, f50Term1M3x3x3x3x3, f50Term2M3x3x3x3x3, f50Term3M3x3x3x3x3) &
    !$OMP SHARED(this, coords, nAtoms, mi3x3, mi3x3otohoi3x3) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtoms
      iSp1 = this%species(iAt1)
      u1 = this%hubbu(iSp1)
      do iAt2 = 1, iAt1 - 1
        iSp2 = this%species(iAt2)
        u2 = this%hubbu(iSp2)
        vRabx3(:) = coords(:,iAt1) - coords(:,iAt2)
        rab = norm2(vRabx3)

        ! f10AB
        gammaPrime = -1.0_dp / rab**2 - expGammaPrime(rab, u1, u2)
        coeffTerm1 = gammaPrime / rab

        workVx3(:) = coeffTerm1 * vRabx3
        this%f10AB(:,iAt2, iAt1) = workVx3
        this%f10AB(:,iAt1, iAt2) = -workVx3

        ! f20AB
        gammaDoublePrime = 2.0_dp / rab**3 - expGammaDoublePrime(rab, u1, u2)
        coeffTerm1 = gammaPrime / rab
        coeffTerm2 = gammaDoublePrime / rab**2 - gammaPrime / rab**3
        coeffTerm1 = coeffTerm1 / 2.0_dp
        coeffTerm2 = coeffTerm2 / 2.0_dp

        call outerProductA3OB3(mRRab3x3, vRabx3, vRabx3)
        workM3x3(:,:) = coeffTerm1 * mI3x3 + coeffTerm2 * mRRab3x3
        this%f20AB(:,:,iAt2, iAt1) = workM3x3
        this%f20AB(:,:,iAt1, iAt2) = workM3x3

        ! f30AB
        gammaTriplePrime = -6.0_dp / rab**4 - expGammaTriplePrime(rab, u1, u2)
        coeffTerm1 = gammaDoublePrime / rab**2 - gammaPrime / rab**3
        coeffTerm2 = gammaTriplePrime / rab**3 - 3.0_dp * gammaDoublePrime / rab**4&
            & + 3.0_dp * gammaPrime / rab**5
        coeffTerm1 = coeffTerm1 / 6.0_dp
        coeffTerm2 = coeffTerm2 / 6.0_dp

        do ll = 1, 3
          do jj = 1, 3
            do ii = 1, 3
              mRoIab3x3x3(ii, jj, ll) = vRabx3(ii) * mI3x3(jj, ll)
              mIotRab3x3x3(ii, jj, ll) = mI3x3(ii, ll) * vRabx3(jj)
              mIoRab3x3x3(ii, jj, ll) = mI3x3(ii, jj) * vRabx3(ll)
              mRRRab3x3x3(ii, jj, ll) = mRRab3x3(ii, jj) * vRabx3(ll)
            end do
          end do
        end do
        mRIotoRab3x3x3(:,:,:) = mRoIab3x3x3 + mIotRab3x3x3 + mIoRab3x3x3
        workM3x3x3(:,:,:) = coeffTerm1 * mRIotoRab3x3x3 + coeffTerm2 * mRRRab3x3x3
        this%f30AB(:,:,:,iAt2, iAt1) = workM3x3x3
        this%f30AB(:,:,:,iAt1, iAt2) = -workM3x3x3

        ! f40AB and f50AB
        if (all(.not. this%hasOnsiteCharges(:, iSp1)) &
            & .or. all(.not. this%hasOnsiteCharges(:, iSp2))) then
          cycle
        end if

        ! f40AB
        gammaQuadruplePrime = 24.0_dp / rab**5 - expGammaQuadruplePrime(rab, u1, u2)
        coeffTerm1 = gammaDoublePrime / rab**2 - gammaPrime / rab**3
        coeffTerm2 = gammaTriplePrime / rab**3 &
            & - 3.0_dp * gammaDoublePrime / rab**4 + 3.0_dp * gammaPrime / rab**5
        coeffTerm3 = gammaQuadruplePrime / rab**4 - 6.0_dp * gammaTriplePrime / rab**5 &
            & + 15.0_dp * gammaDoublePrime / rab**6 - 15.0_dp * gammaPrime / rab**7
        coeffTerm1 = coeffTerm1 / 24.0_dp
        coeffTerm2 = coeffTerm2 / 24.0_dp
        coeffTerm3 = coeffTerm3 / 24.0_dp


        f40Term1M3x3x3x3(:,:,:,:) = mI3x3otohoI3x3
        f40Term2M3x3x3x3(:,:,:,:) = 0.0_dp
        f40Term3M3x3x3x3(:,:,:,:) = 0.0_dp
        do mm = 1, 3
          do ll = 1, 3
            do jj = 1, 3
              do ii = 1, 3
                f40Term2M3x3x3x3(ii, jj, ll, mm) = f40Term2M3x3x3x3(ii, jj, ll, mm) &
                    & + mRRab3x3(ii, jj) * mI3x3(ll, mm) + mRRab3x3(ii, ll) * mI3x3(jj, mm) &
                    & + mRRab3x3(ii, mm) * mI3x3(jj, ll) + mI3x3(ii, jj) * mRRab3x3(ll, mm) &
                    & + mI3x3(ii, ll) * mRRab3x3(jj, mm) + mI3x3(ii, mm) * mRRab3x3(jj, ll)
                f40Term3M3x3x3x3(ii, jj, ll, mm) = f40Term3M3x3x3x3(ii, jj, ll, mm) &
                    & + mRRab3x3(ii, jj) * mRRab3x3(ll, mm)
              end do
            end do
          end do
        end do
        workM3x3x3x3(:,:,:,:) = coeffTerm1 * f40Term1M3x3x3x3 + coeffTerm2 * f40Term2M3x3x3x3 &
            & + coeffTerm3 * f40Term3M3x3x3x3
        this%f40AB(:,:,:,:,iAt2, iAt1) = workM3x3x3x3
        this%f40AB(:,:,:,:,iAt1, iAt2) = workM3x3x3x3

        ! f50AB
        gammaQuintuplePrime = -120.0_dp / rab**6 - expGammaQuintuplePrime(rab, u1, u2)
        coeffTerm1 = gammaTriplePrime / rab**3 &
            & - 3.0_dp * gammaDoublePrime / rab**4 + 3.0_dp * gammaPrime / rab**5
        coeffTerm2 = gammaQuadruplePrime / rab**4 - 6.0_dp * gammaTriplePrime / rab**5 &
            & + 15.0_dp * gammaDoublePrime / rab**6 - 15.0_dp * gammaPrime / rab**7
        coeffTerm3 = gammaQuintuplePrime / rab**5 &
            & - 10.0_dp * gammaQuadruplePrime / rab**6 + 45.0_dp * gammaTriplePrime / rab**7 &
            & - 105.0_dp * gammaDoublePrime / rab**8 + 105.0_dp * gammaPrime / rab**9
        coeffTerm1 = coeffTerm1 / 120.0_dp
        coeffTerm2 = coeffTerm2 / 120.0_dp
        coeffTerm3 = coeffTerm3 / 120.0_dp

        f50Term1M3x3x3x3x3(:,:,:,:,:) = 0.0_dp
        f50Term2M3x3x3x3x3(:,:,:,:,:) = 0.0_dp
        f50Term3M3x3x3x3x3(:,:,:,:,:) = 0.0_dp
        do nn = 1, 3
          do mm = 1, 3
            do ll = 1, 3
              do jj = 1, 3
                do ii = 1, 3
                  f50Term1M3x3x3x3x3(ii, jj, ll, mm, nn) = f50Term1M3x3x3x3x3(ii, jj, ll, mm, nn) &
                      & + mRoIab3x3x3(ii, jj, ll) * mI3x3(mm, nn) &
                      & + mRoIab3x3x3(ii, jj, mm) * mI3x3(ll, nn) &
                      & + mRoIab3x3x3(ii, jj, nn) * mI3x3(ll, mm) &
                      & + mI3x3(ii, jj) * mRIotoRab3x3x3(ll, mm, nn) &
                      & + mI3x3(ii, ll) * mRIotoRab3x3x3(jj, mm, nn) &
                      & + mI3x3(ii, mm) * mRIotoRab3x3x3(jj, ll, nn) &
                      & + mI3x3(ii, nn) * mRIotoRab3x3x3(jj, ll, mm)
                  f50Term2M3x3x3x3x3(ii, jj, ll, mm, nn) = f50Term2M3x3x3x3x3(ii, jj, ll, mm, nn) &
                      & + mRRRab3x3x3(ii, jj, ll) * mI3x3(mm, nn) &
                      & + mRoIab3x3x3(ii, jj, nn) * mRRab3x3(ll, mm) &
                      & + mI3x3(ii, jj) * mRRRab3x3x3(ll, mm, nn) &
                      & + mI3x3(ii, ll) * mRRRab3x3x3(jj, mm, nn) &
                      & + mI3x3(ii, mm) * mRRRab3x3x3(jj, ll, nn) &
                      & + mI3x3(ii, nn) * mRRRab3x3x3(jj, ll, mm) &
                      & + mRRab3x3(ii, jj) * mIoRab3x3x3(ll, mm, nn) &
                      & + mRRab3x3(ii, ll) * mIoRab3x3x3(jj, mm, nn) &
                      & + mRRab3x3(ii, mm) * mIoRab3x3x3(jj, ll, nn) &
                      & + mRRab3x3(ii, jj) * mIotRab3x3x3(ll, mm, nn)
                  f50Term3M3x3x3x3x3(ii, jj, ll, mm, nn) = f50Term3M3x3x3x3x3(ii, jj, ll, mm, nn) &
                      & + mRRRab3x3x3(ii, jj, ll) * mRRab3x3(mm, nn)
                end do
              end do
            end do
          end do
        end do
        workM3x3x3x3x3(:,:,:,:,:) = coeffTerm1 * f50Term1M3x3x3x3x3 &
            & + coeffTerm2 * f50Term2M3x3x3x3x3 + coeffTerm3 * f50Term3M3x3x3x3x3

        this%f50AB(:,:,:,:,:,iAt2, iAt1) = workM3x3x3x3x3
        this%f50AB(:,:,:,:,:,iAt1, iAt2) = - workM3x3x3x3x3
      end do

      coeffTerm1 = - (3.2_dp * u1)**3 / 96.0_dp
      this%f20AB(:,:,iAt1, iAt1) = coeffTerm1 * mI3x3

      coeffTerm1 = (3.2_dp * u1)**5 / 5760.0_dp
      this%f40AB(:,:,:,:,iAt1, iAt1) = coeffTerm1 * mI3x3otohoI3x3

    end do
    !$OMP END PARALLEL DO

  end subroutine updateCoords


  !> Updates the deltaDAtom and deltaQAtom arrays for the instance.
  subroutine updateDeltaDQAtom(this, over, rho, orb, iNeighbour, nNeighbourSK, img2CentCell, iPair)

    !> Class instance
    class(TMdftb), intent(inout) :: this

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

    real(dp), allocatable :: tmpOvr(:,:), tmpDRho(:,:)
    integer, parameter :: iStart = 1, iEnd = 2, iNOrb = 3

    integer :: nAtom, iAtom1, iAtom2, iAtom2f, iSp1, iSp2, nOrb1, nOrb2, iOrig, iNeigh
    integer :: ii, jj, mu, nu, kk
    real(dp) :: mulTmp(orb%mOrb**2), sqrTmp(orb%mOrb, orb%mOrb)
    real(dp) :: sqrTmpOver(orb%mOrb, orb%mOrb), sqrTmpRho(orb%mOrb, orb%mOrb)
    real(dp) :: tmpVx3(3), tmpM3x3(3,3)

    nAtom = this%nAtoms
    @:ASSERT(size(over) == size(rho))

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
        sqrTmp(1:nOrb2, 1:nOrb1) = reshape(mulTmp(1:nOrb1*nOrb2), [nOrb2, nOrb1])
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
        sqrTmpOver(1:nOrb2, 1:nOrb1) = reshape(over(iOrig:iOrig+nOrb1*nOrb2-1), [nOrb2, nOrb1])
        sqrTmpRho(1:nOrb2, 1:nOrb1) = reshape(rho(iOrig:iOrig+nOrb1*nOrb2-1), [nOrb2, nOrb1])

        tmpVx3(:) = 0.0_dp
        tmpM3x3(:,:) = 0.0_dp
        do mu = 1, nOrb1
          do nu =  1, nOrb1
            do kk =  1, nOrb2
              tmpVx3(:) =  tmpVx3 + sqrTmpRho(kk, mu) * sqrTmpOver(kk, nu)&
                  & * this%atomicDIntgrl(:,nu, mu, iSp1)
              tmpM3x3(:,:) = tmpM3x3 + sqrTmpRho(kk, mu) * sqrTmpOver(kk, nu)&
                  & * this%atomicQIntgrl(:,:,nu, mu, iSp1)
            end do
          end do
        end do
        this%deltaDAtom(:,iAtom1) = this%deltaDAtom(:,iAtom1) + tmpVx3
        this%deltaQAtom(:,:,iAtom1) = this%deltaQAtom(:,:,iAtom1) + tmpM3x3

        ! Add contribution to the other triangle sum, using the symmetry, but only when off diagonal
        if (iAtom1 /= iAtom2f) then
          tmpVx3(:) = 0.0_dp
          tmpM3x3(:,:) = 0.0_dp
          do mu = 1, nOrb2
            do nu =  1, nOrb2
              do kk =  1, nOrb1
                tmpVx3(:) =  tmpVx3 + sqrTmpRho(mu, kk) * sqrTmpOver(nu, kk)&
                    & * this%atomicDIntgrl(:,nu, mu, iSp2)
                tmpM3x3(:,:) = tmpM3x3 + sqrTmpRho(mu, kk) * sqrTmpOver(nu, kk)&
                    & * this%atomicQIntgrl(:,:,nu, mu, iSp2)
              end do
            end do
          end do
          this%deltaDAtom(:,iAtom2f) = this%deltaDAtom(:,iAtom2f) + tmpVx3
          this%deltaQAtom(:,:,iAtom2f) = this%deltaQAtom(:,:,iAtom2f) + tmpM3x3
        end if
      end do

    end do

  end subroutine updateDeltaDQAtom


  !> Receives the deltaDAtom and deltaQAtom to update mdftb.
  subroutine pullDeltaDQAtom(this, multiExpanData)

    !> class instance
    class(TMdftb), intent(inout) :: this

    !> Multipole moments to pull
    type(TMultipole), intent(in) :: multiExpanData

    integer :: nAtoms
    integer :: iAt1, iSp1, iSpin
    real(dp) :: tmpVx3(3), tmpM3x3(3,3)

    iSpin = 1
    nAtoms = this%nAtoms
    this%deltaDAtom(:,:) = multiExpanData%dipoleAtom(:,:,iSpin)
    do iAt1 = 1, nAtoms
      ! Quadrupole components used (xx, xy, yy, xz, yz, zz)
      tmpM3x3(:,:) = 0.0_dp
      tmpM3x3(1,1) = multiExpanData%quadrupoleAtom(1, iAt1, iSpin)
      tmpM3x3(2,1) = multiExpanData%quadrupoleAtom(2, iAt1, iSpin)
      tmpM3x3(2,2) = multiExpanData%quadrupoleAtom(3, iAt1, iSpin)
      tmpM3x3(3,1) = multiExpanData%quadrupoleAtom(4, iAt1, iSpin)
      tmpM3x3(3,2) = multiExpanData%quadrupoleAtom(5, iAt1, iSpin)
      tmpM3x3(3,3) = multiExpanData%quadrupoleAtom(6, iAt1, iSpin)
      call adjointLowerTriangle(tmpM3x3)
      this%deltaQAtom(:,:,iAt1) = tmpM3x3
    end do

  end subroutine pullDeltaDQAtom


  !> Broadcasts the deltaDAtom and deltaQAtom from mdftb.
  subroutine pushDeltaDQAtom(this, multiExpanData)

    !> class instance
    class(TMdftb), intent(inout) :: this

    !> Multipole moments push
    type(TMultipole), intent(inout) :: multiExpanData

    integer :: nAtoms
    integer :: iAt1, iSp1, iSpin

    real(dp) :: tmpVx3(3), tmpM3x3(3,3)

    iSpin = 1
    nAtoms = this%nAtoms
    multiExpanData%dipoleAtom(:,:,iSpin) = this%deltaDAtom
    do iAt1 = 1, nAtoms
      ! Quadrupole components used (xx, xy, yy, xz, yz, zz)
      tmpM3x3(:,:) = this%deltaQAtom(:,:,iAt1)
      multiExpanData%quadrupoleAtom(1, iAt1, iSpin) = tmpM3x3(1,1) 
      multiExpanData%quadrupoleAtom(2, iAt1, iSpin) = tmpM3x3(2,1) 
      multiExpanData%quadrupoleAtom(3, iAt1, iSpin) = tmpM3x3(2,2) 
      multiExpanData%quadrupoleAtom(4, iAt1, iSpin) = tmpM3x3(3,1) 
      multiExpanData%quadrupoleAtom(5, iAt1, iSpin) = tmpM3x3(3,2) 
      multiExpanData%quadrupoleAtom(6, iAt1, iSpin) = tmpM3x3(3,3) 
    end do

  end subroutine pushDeltaDQAtom


  !> Updates the mdftb potentials.
  subroutine updateDQPotentials(this, deltaMAtom)

    !> class instance
    class(TMdftb), intent(inout) :: this

    !> Delta charge per atom
    real(dp), intent(in) :: deltaMAtom(:)

    integer :: nAtoms, iAt1, iAt2, ii, jj, iSp1, iSp2

    this%deltaMAtom(:) = deltaMAtom
    nAtoms = this%nAtoms

    ! Calculate pot10x0, and pot20x0
    this%pot10x0Atom(:,:) = 0.0_dp
    this%pot20x0Atom(:,:,:) = 0.0_dp

    !$OMP PARALLEL DO PRIVATE(iAt1, iAt2) DEFAULT(SHARED) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtoms
      do iAt2 = 1, nAtoms
        this%pot10x0Atom(:,iAt1) = this%pot10x0Atom(:,iAt1)&
            & + this%f10AB(:,iAt2, iAt1) * this%deltaMAtom(iAt2)
        this%pot20x0Atom(:,:,iAt1) = this%pot20x0Atom(:,:,iAt1)&
            & + this%f20AB(:,:,iAt2, iAt1) * this%deltaMAtom(iAt2)
      end do
    end do
    !$OMP END PARALLEL DO

    ! Calculate pot10x1, pot20x2, pot11x1, pot21x2, and pot21x1
    this%pot10x1Atom(:) = 0.0_dp
    this%pot20x2Atom(:) = 0.0_dp
    this%pot11x1Atom(:,:) = 0.0_dp
    this%pot21x2Atom(:,:) = 0.0_dp
    this%pot21x1Atom(:,:,:) = 0.0_dp

    !$OMP PARALLEL DO PRIVATE(iAt1, iAt2, iSp2, ii, jj) DEFAULT(SHARED) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtoms
      do iAt2 = 1, nAtoms
        iSp2 = this%species(iAt2)
        if (all(.not. this%hasOnsiteCharges(:, iSp2))) then
          cycle
        end if
        this%pot10x1Atom(iAt1) = this%pot10x1Atom(iAt1)&
            & + sum(this%f10AB(:,iAt2, iAt1) * this%deltaDAtom(:,iAt2))
        this%pot20x2Atom(iAt1) = this%pot20x2Atom(iAt1)&
            & + sum(this%f20AB(:,:,iAt2, iAt1) * this%deltaQAtom(:,:,iAt2))
        do ii = 1, 3
          this%pot11x1Atom(ii, iAt1) = this%pot11x1Atom(ii, iAt1)&
              & - 2.0_dp * sum(this%f20AB(:,ii, iAt2, iAt1) * this%deltaDAtom(:,iAt2))
          this%pot21x2Atom(ii, iAt1) = this%pot21x2Atom(ii, iAt1)&
              & - 3.0_dp * sum(this%f30AB(:,:,ii, iAt2, iAt1) * this%deltaQAtom(:,:,iAt2))
          do jj = 1, 3
            this%pot21x1Atom(jj, ii, iAt1) = this%pot21x1Atom(jj, ii, iAt1)&
                & - 3.0_dp * sum(this%f30AB(:,jj, ii, iAt2, iAt1) * this%deltaDAtom(:,iAt2))
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    ! Calculate pot21x2
    this%pot22x2Atom(:,:,:) = 0.0_dp

    !$OMP PARALLEL DO PRIVATE(iAt1, iSp1, iAt2, iSp2, ii, jj) DEFAULT(SHARED) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtoms
      iSp1 = this%species(iAt1)
      if (all(.not. this%hasOnsiteCharges(:, iSp1))) then
        cycle
      end if
      do iAt2 = 1, nAtoms
        iSp2 = this%species(iAt2)
        if (all(.not. this%hasOnsiteCharges(:, iSp2))) then
          cycle
        end if
        do ii = 1, 3
          do jj = 1, 3
            this%pot22x2Atom(jj, ii, iAt1) = this%pot22x2Atom(jj, ii, iAt1)&
                & + 6.0_dp * sum(this%f40AB(:,:,jj, ii, iAt2, iAt1) * this%deltaQAtom(:,:,iAt2))
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine updateDQPotentials


  !> Adds the multipole expansion contribution to the Hamiltonian.
  subroutine addMultiExpanHamiltonian(this, ham, over, nNeighbour, iNeighbour, species, orb,&
      & iPair, nAtom, img2CentCell)

    !> class instance
    class(TMdftb), intent(inout) :: this

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

    integer :: mu, nu, kappa, mm, nn, kk
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: iSp1, iSp2, iNeigh, iOrig
    integer :: nOrb1, nOrb2
    real(dp) :: sqrTmpOver(orb%mOrb, orb%mOrb)
    real(dp) :: sqrTmpOverT(orb%mOrb, orb%mOrb)
    real(dp) :: sqrTmpHam(orb%mOrb, orb%mOrb)
    real(dp) :: tmpVx3(3), tmpM3x3(3,3)
    real(dp) :: tmpadS(3), tmpSad(3)
    real(dp) :: tmpaQS(3,3), tmpSaQ(3,3)

    @:ASSERT(size(nNeighbour)==nAtom)
    @:ASSERT(size(iNeighbour, dim=2)==nAtom)
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
        sqrTmpOver(1:nOrb2, 1:nOrb1) = reshape(over(iOrig:iOrig+nOrb1*nOrb2-1), [nOrb2, nOrb1])
        sqrTmpOverT(1:nOrb1, 1:nOrb2) = transpose(sqrTmpOver(1:nOrb2, 1:nOrb1)) 

        sqrTmpHam(:,:) = 0.0_dp
        do mu = 1, nOrb1
          do nu =  1, nOrb2
            tmpadS(:) = 0.0_dp
            tmpSad(:) = 0.0_dp
            tmpaQS(:,:) = 0.0_dp
            tmpSaQ(:,:) = 0.0_dp
            do kk = 1, nOrb1
              tmpadS(:) = tmpadS + this%atomicDIntgrl(:,kk,mu,iSp1) * sqrTmpOverT(kk,nu)
              tmpaQS(:,:) = tmpaQS + this%atomicQIntgrl(:,:,kk,mu,iSp1) * sqrTmpOverT(kk,nu)
            end do
            do kk = 1, nOrb2
              tmpSad(:) = tmpSad + this%atomicDIntgrl(:,kk,nu,iSp2) * sqrTmpOver(kk,mu)
              tmpSaQ(:,:) = tmpSaQ + this%atomicQIntgrl(:,:,kk,nu,iSp2) * sqrTmpOver(kk,mu)
            end do

            ! Add Monopole-Dipole contribution
            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) - (this%pot10x1Atom(iAtom1)&
                & + this%pot10x1Atom(iAtom2f)) * sqrTmpOver(nu,mu)
            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) + sum(this%pot10x0Atom(:,iAtom1) * tmpadS) &
                & + sum(this%pot10x0Atom(:,iAtom2f) * tmpSad)

            ! Add Dipole-Dipole contribution
            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) + sum(this%pot11x1Atom(:,iAtom1) * tmpadS) &
                & + sum(this%pot11x1Atom(:,iAtom2f) * tmpSad)
 
            ! Add Monopole-Quadrupole contribution
            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) + (this%pot20x2Atom(iAtom1)&
                & + this%pot20x2Atom(iAtom2f)) * sqrTmpOver(nu,mu)
            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) + sum(this%pot20x0Atom(:,:,iAtom1) * tmpaQS) &
                & + sum(this%pot20x0Atom(:,:,iAtom2f) * tmpSaQ)

            ! Add Dipole-Quadrupole contribution
            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) - sum(this%pot21x2Atom(:,iAtom1) * tmpadS) &
                & - sum(this%pot21x2Atom(:,iAtom2f) * tmpSad)
            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) + sum(this%pot21x1Atom(:,:,iAtom1) * tmpaQS) &
                & + sum(this%pot21x1Atom(:,:,iAtom2f) * tmpSaQ)

            ! Add Quadrupole-Quadrupole contribution
            sqrTmpHam(nu,mu) = sqrTmpHam(nu,mu) + sum(this%pot22x2Atom(:,:,iAtom1) * tmpaQS) &
                & + sum(this%pot22x2Atom(:,:,iAtom2f) * tmpSaQ)
          end do
        end do

        ham(iOrig:iOrig+nOrb1*nOrb2-1, 1) = ham(iOrig:iOrig+nOrb1*nOrb2-1, 1)&
           & + 0.5_dp * reshape(sqrTmpHam(1:nOrb2, 1:nOrb1), [nOrb2*nOrb1])
      end do
    end do

  end subroutine addMultiExpanHamiltonian


  !> Add the MultiExpan-Energy contribution to the total energy.
  subroutine addMultiExpanEnergy(this, energyPerAtomTT, energyMD, energyDD, energyMQ, energyDQ,&
      & energyQQ, energyTT)

    !> Mdftb class instance
    class(TMdftb), intent(in) :: this

    !> Energy per atom
    real(dp), intent(out) :: energyPerAtomTT(:)

    !> Energies for the different contributions
    real(dp), intent(out) :: energyMD, energyDD, energyMQ, energyDQ, energyQQ, energyTT

    real(dp), allocatable :: EnergyMDAtom(:), EnergyDDAtom(:), EnergyMQAtom(:)
    real(dp), allocatable :: EnergyDQAtom(:), EnergyQQAtom(:)

    integer :: nAtoms
    integer :: iAt1

    nAtoms = this%nAtoms
    @:ASSERT(size(energyPerAtomTT) == nAtoms)

    allocate(EnergyMDAtom(nAtoms), source=0.0_dp)
    allocate(EnergyDDAtom(nAtoms), source=0.0_dp)
    allocate(EnergyMQAtom(nAtoms), source=0.0_dp)
    allocate(EnergyDQAtom(nAtoms), source=0.0_dp)
    allocate(EnergyQQAtom(nAtoms), source=0.0_dp)

    do iAt1 = 1, nAtoms
      ! Monopole-Dipole contribution
      EnergyMDAtom(iAt1) = sum(this%pot10x0Atom(:,iAt1) * this%deltaDAtom(:,iAt1))
      ! Dipole-Dipole contribution
      EnergyDDAtom(iAt1) = 0.5_dp * sum(this%pot11x1Atom(:,iAt1) * this%deltaDAtom(:,iAt1))
      ! Monopole-Quadrupole contribution
      EnergyMQAtom(iAt1) = sum(this%pot20x0Atom(:,:,iAt1) * this%deltaQAtom(:,:,iAt1))
      ! Dipole-Quadrupole contribution
      EnergyDQAtom(iAt1) = sum(this%pot21x1Atom(:,:,iAt1) * this%deltaQAtom(:,:,iAt1))
      ! Quadrupole-Quadrupole contribution
      EnergyQQAtom(iAt1) = 0.5_dp * sum(this%pot22x2Atom(:,:,iAt1) * this%deltaQAtom(:,:,iAt1))
    end do

    energyPerAtomTT(:) = EnergyMDAtom + EnergyDDAtom + EnergyMQAtom + EnergyDQAtom + EnergyQQAtom
    energyMD = sum(EnergyMDAtom)
    energyDD = sum(EnergyDDAtom)
    energyMQ = sum(EnergyMQAtom)
    energyDQ = sum(EnergyDQAtom)
    energyQQ = sum(EnergyQQAtom)
    energyTT = energyMD + energyDD + energyMQ  + energyDQ  + energyQQ

    deallocate(EnergyMDAtom, EnergyDDAtom, EnergyMQAtom, EnergyDQAtom, EnergyQQAtom)

  end subroutine addMultiExpanEnergy


  !> Add the mdftb contribution to the gradients.
  subroutine addMultiExpanGradients(this, derivs, derivator, skOverCont, rho, species,&
      & iNeighbour, nNeighbourSK, img2CentCell, iPair, orb, coords)

    !> Class instance
    class(TMdftb), intent(inout) :: this

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

    integer :: ii, jj, ll, mu, nu, kappa, mm, nn, kk
    integer :: nOrb1, nOrb2
    integer :: iAt1, iAt2, iAt2f
    integer :: iOrig, iSpin, nSpin, nAtom, iNeigh, iSp1, iSp2

    real(dp) :: sqrDMTmp(orb%mOrb,orb%mOrb), sqrDMTmpT(orb%mOrb,orb%mOrb)
    real(dp) :: sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: derivTmp(3)
    real(dp) :: tmpDerivMuNu

    real(dp) :: tmpPad(3), tmpadP(3)
    real(dp) :: tmpPaQ(3,3), tmpaQP(3,3)

    ! Calculate pot30x0, pot31x1, and pot32x2
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
        if (all(.not. this%hasOnsiteCharges(:, iSp1)) &
            & .or. all(.not. this%hasOnsiteCharges(:, iSp2))) then
          cycle
        end if
        do ii = 1, 3
          do jj = 1, 3
            do ll = 1, 3
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
      ! M-D contribution
      derivs(:,iAt1) = derivs(:,iAt1) + this%pot11x1Atom(:,iAt1) * this%deltaMAtom(iAt1)
      ! M-Q contribution
      derivs(:,iAt1) = derivs(:,iAt1) - this%pot21x2Atom(:,iAt1) * this%deltaMAtom(iAt1)
      do ii = 1, 3
        ! D-M contribution
        derivs(ii,iAt1) = derivs(ii,iAt1)&
            & + 2.0_dp * sum(this%pot20x0Atom(:,ii,iAt1) * this%deltaDAtom(:,iAt1))
        ! D-D contribution
        derivs(ii,iAt1) = derivs(ii,iAt1)&
            & + 2.0_dp * sum(this%pot21x1Atom(:,ii,iAt1) * this%deltaDAtom(:,iAt1))
        ! D-Q contribution
        derivs(ii,iAt1) = derivs(ii,iAt1)&
            & + 2.0_dp * sum(this%pot22x2Atom(:,ii,iAt1) * this%deltaDAtom(:,iAt1))
        ! Q-M contribution
        derivs(ii,iAt1) = derivs(ii,iAt1)&
            & + 3.0_dp * sum(this%pot30x0Atom(:,:,ii,iAt1) * this%deltaQAtom(:,:,iAt1))
        ! Q-D contribution
        derivs(ii,iAt1) = derivs(ii,iAt1)&
            & + 3.0_dp * sum(this%pot31x1Atom(:,:,ii,iAt1) * this%deltaQAtom(:,:,iAt1))
        ! Q-Q contribution
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

        derivTmp(:) = 0.0_dp
        do mu = 1, nOrb1
          do nu =  1, nOrb2

            tmpPad(:) = 0.0_dp
            tmpadP(:) = 0.0_dp
            tmpPaQ(:,:) = 0.0_dp
            tmpaQP(:,:) = 0.0_dp
            do kk = 1, nOrb1
              tmpPad(:) = tmpPad + this%atomicDIntgrl(:,kk,mu,iSp1) * sqrDMTmpT(kk,nu)
              tmpPaQ(:,:) = tmpPaQ + this%atomicQIntgrl(:,:,kk,mu,iSp1) * sqrDMTmpT(kk,nu)
            end do
            do kk = 1, nOrb2
              tmpadP(:) = tmpadP + this%atomicDIntgrl(:,kk,nu,iSp2) * sqrDMTmp(kk,mu)
              tmpaQP(:,:) = tmpaQP + this%atomicQIntgrl(:,:,kk,nu,iSp2) * sqrDMTmp(kk,mu)
            end do

            tmpDerivMuNu = 0.0_dp
            ! M-D contribution
            tmpDerivMuNu = tmpDerivMuNu - (this%pot10x1Atom(iAt1) + this%pot10x1Atom(iAt2))&
                & * sqrDMTmp(nu,mu)
            ! M-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + (this%pot20x2Atom(iAt1) + this%pot20x2Atom(iAt2))&
                & * sqrDMTmp(nu,mu)

            ! D-M contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot10x0Atom(:,iAt1) * tmpPad)
            ! D-D contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot11x1Atom(:,iAt1) * tmpPad)
            ! D-Q contribution
            tmpDerivMuNu = tmpDerivMuNu - sum(this%pot21x2Atom(:,iAt1) * tmpPad)
            ! M-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot20x0Atom(:,:,iAt1) * tmpPaQ)
            ! D-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot21x1Atom(:,:,iAt1) * tmpPaQ)
            ! Q-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot22x2Atom(:,:,iAt1) * tmpPaQ)
            
            ! M-D contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot10x0Atom(:,iAt2) * tmpadP)
            ! D-D contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot11x1Atom(:,iAt2) * tmpadP)
            ! D-Q contribution
            tmpDerivMuNu = tmpDerivMuNu - sum(this%pot21x2Atom(:,iAt2) * tmpadP)
            ! M-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot20x0Atom(:,:,iAt2) * tmpaQP)
            ! D-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot21x1Atom(:,:,iAt2) * tmpaQP)
            ! Q-Q contribution
            tmpDerivMuNu = tmpDerivMuNu + sum(this%pot22x2Atom(:,:,iAt2) * tmpaQP)

            derivTmp(:) = derivTmp + tmpDerivMuNu * sPrimeTmp(nu,mu,:)

          end do
        end do
        ! forces from atom 1 on atom 2f and 2f onto 1
        derivs(:,iAt1) = derivs(:,iAt1) + derivTmp
        derivs(:,iAt2f) = derivs(:,iAt2f) - derivTmp
      end do
    end do

  end subroutine addMultiExpanGradients


  !> Adds the atomic dipole moments to the total dipole moments.
  subroutine addAtomicDipoleMoment(this, dipoleMoment)

    !> Class instance
    class(TMdftb), intent(inout) :: this

    !> The total dipole moments to be updated.
    real(dp), intent(inout) :: dipoleMoment(:)

    dipoleMoment(:) = dipoleMoment - sum(this%deltaDAtom, dim=2)

  end subroutine addAtomicDipoleMoment


  !> Adds the atomic quadrupole moments to the total quadrupole moments.
  subroutine addAtomicQuadrupoleMoment(this, quadrupoleMoment)

    !> Class instance
    class(TMdftb), intent(inout) :: this

    !> The total quadrupole moments to be updated.
    real(dp), intent(inout) :: quadrupoleMoment(:,:)

    quadrupoleMoment(:,:) = quadrupoleMoment - sum(this%deltaQAtom, dim=3)

  end subroutine addAtomicQuadrupoleMoment


  !> Computes the outer product of two 3-element vectors: M(i,j) = A(i) * B(j)
  subroutine outerProductA3OB3(M3x3, A3, B3)

    !> Output 3×3 matrix containing the outer product result.
    real(dp), intent(out) :: M3x3(:,:)

    !> Input 3-element vector.
    real(dp), intent(in)  :: A3(:)

    !> Input 3-element vector.
    real(dp), intent(in)  :: B3(:)

    @:ASSERT(size(A3) == 3)
    @:ASSERT(size(B3) == 3)
    @:ASSERT(size(M3x3, dim=1) == 3)
    @:ASSERT(size(M3x3, dim=2) == 3)
  
    M3x3(:,:) = spread(A3, 2, 3) * spread(B3, 1, 3)
  
  end subroutine outerProductA3OB3


  !> Computes the outer product of two 3×3 matrices: M(i,j,l,m) = A(i,j) * B(l,m)
  subroutine outerProductA3x3OB3x3(M3x3x3x3, A3x3, B3x3)

    !> Output 4-dimensional tensor (3×3×3×3) containing the outer product result.
    real(dp), intent(out) :: M3x3x3x3(:,:,:,:)

    !> Input 3×3 matrix A3x3.
    real(dp), intent(in)  :: A3x3(:,:)

    !> Input 3×3 matrix B3x3.
    real(dp), intent(in)  :: B3x3(:,:)

    @:ASSERT(size(B3x3, dim=1) == 3)
    @:ASSERT(size(B3x3, dim=2) == 3)
    @:ASSERT(size(A3x3, dim=1) == 3)
    @:ASSERT(size(A3x3, dim=2) == 3)
    @:ASSERT(size(M3x3x3x3, dim=1) == 3)
    @:ASSERT(size(M3x3x3x3, dim=2) == 3)
    @:ASSERT(size(M3x3x3x3, dim=3) == 3)
    @:ASSERT(size(M3x3x3x3, dim=4) == 3)
  
    M3x3x3x3 = spread(spread(A3x3, 3, 3), 4, 3) * spread(spread(B3x3, 1, 3), 1, 3)
  
  end subroutine outerProductA3x3OB3x3
  
  
  !> Computes the outer product of two 3×3 matrices: M(i,j,l,m) = A3x3(i,l) * B3x3(j,m)
  subroutine outerProductA3x3OhB3x3(M3x3x3x3, A3x3, B3x3)

    !> Output 4-dimensional tensor (3×3×3×3) containing the outer product result.
    real(dp), intent(out) :: M3x3x3x3(:,:,:,:)
  
    !> Input 3×3 matrix A3x3.
    real(dp), intent(in)  :: A3x3(:,:)
  
    !> Input 3×3 matrix B3x3.
    real(dp), intent(in)  :: B3x3(:,:)
  
    integer :: ii, jj, ll, mm

    @:ASSERT(size(B3x3, dim=1) == 3)
    @:ASSERT(size(B3x3, dim=2) == 3)
    @:ASSERT(size(A3x3, dim=1) == 3)
    @:ASSERT(size(A3x3, dim=2) == 3)
    @:ASSERT(size(M3x3x3x3, dim=1) == 3)
    @:ASSERT(size(M3x3x3x3, dim=2) == 3)
    @:ASSERT(size(M3x3x3x3, dim=3) == 3)
    @:ASSERT(size(M3x3x3x3, dim=4) == 3)
    do mm = 1, 3
      do ll = 1, 3
        do jj = 1, 3
          do ii = 1, 3
            M3x3x3x3(ii,jj,ll,mm) = A3x3(ii,ll) * B3x3(jj,mm)
          end do
        end do
      end do
    end do

  end subroutine outerProductA3x3OhB3x3


  !> Computes the outer product of two 3×3 matrices: M(i,j,l,m) = A3x3(i,m) * B3x3(j,l)
  subroutine outerProductA3x3OtB3x3(M3x3x3x3, A3x3, B3x3)

    !> Output 4-dimensional tensor (3×3×3×3) containing the outer product result.
    real(dp), intent(out) :: M3x3x3x3(:,:,:,:)

    !> Input 3×3 matrix A3x3.
    real(dp), intent(in)  :: A3x3(:,:)

    !> Input 3×3 matrix B3x3.
    real(dp), intent(in)  :: B3x3(:,:)

    integer :: ii, jj, ll, mm

    @:ASSERT(size(B3x3, dim=1) == 3)
    @:ASSERT(size(B3x3, dim=2) == 3)
    @:ASSERT(size(A3x3, dim=1) == 3)
    @:ASSERT(size(A3x3, dim=2) == 3)
    @:ASSERT(size(M3x3x3x3, dim=1) == 3)
    @:ASSERT(size(M3x3x3x3, dim=2) == 3)
    @:ASSERT(size(M3x3x3x3, dim=3) == 3)
    @:ASSERT(size(M3x3x3x3, dim=4) == 3)
    do mm = 1, 3
      do ll = 1, 3
        do jj = 1, 3
          do ii = 1, 3
            M3x3x3x3(ii,jj,ll,mm) = A3x3(ii,mm) * B3x3(jj,ll)
          end do
        end do
      end do
    end do

  end subroutine outerProductA3x3OtB3x3


end module dftbp_dftb_mdftb
