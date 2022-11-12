!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines to add contribution of range-separated hybrid functional to onsite correction
!> from doi: 10.1021/acs.jctc.2c00037
module dftbp_dftb_rangeseponscorr
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment, globalTimers
  use dftbp_dftb_sparse2dense, only : symmetrizeHS
  use dftbp_math_blasroutines, only : gemm, gemv
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: TRangeSepOnsCorrFunc, RangeSepOnsCorrFunc_init


  !> Onsite correction with range-separated hybrid functional module
  type :: TRangeSepOnsCorrFunc
    private

    !> Is this spin restricted (F) or unrestricted (T)
    logical :: tSpin

  contains

    procedure :: addLrOcHamiltonian

  end type TRangeSepOnsCorrFunc


contains


  !> Intitialize the onsite correction with range-separated hybrid functional module
  subroutine RangeSepOnsCorrFunc_init(this, tSpin)

    !> class instance
    type(TRangeSepOnsCorrFunc), intent(out) :: this

    !> Is this spin restricted (F) or unrestricted (T)
    logical, intent(in) :: tSpin

    call initialize(this, tSpin)

  contains

    !> initialise data structures
    subroutine initialize(this, tSpin)

      !> class instance
      class(TRangeSepOnsCorrFunc), intent(out) :: this

      !> Is this spin restricted (F) or unrestricted (T)
      logical, intent(in) :: tSpin

      this%tSpin = tSpin

    end subroutine initialize

  end subroutine RangeSepOnsCorrFunc_init


  !> Update Hamiltonian with long-range onsite contribution using matrix-matrix multiplications
  !>
  !> The routine provides a matrix-matrix multiplication based implementation of
  !> Eq. 11 in https://doi.org/10.1021/acs.jctc.2c00037
  !>
  subroutine addLrOcHamiltonian(this, env, orb, iSquare, species, onSiteElements,&
      & overlap, densSqr, HH)

    !> class instance
    class(TRangeSepOnsCorrFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, dimension(:), intent(in) :: iSquare

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> Correction to energy from on-site matrix elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    !> Square (unpacked) overlap matrix.
    real(dp), intent(in) :: overlap(:,:)

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: densSqr(:,:)

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout) :: HH(:,:)

    real(dp), allocatable :: Smat(:,:)
    real(dp), allocatable :: Dmat(:,:)
    real(dp), allocatable :: Omat0(:,:)
    real(dp), allocatable :: OmatRI(:,:)
    real(dp), allocatable :: HlrOC(:,:)

    integer :: nOrb

    nOrb = size(overlap,dim=1)

    allocate(Smat(nOrb,nOrb))
    allocate(Dmat(nOrb,nOrb))
    allocate(Omat0(nOrb,nOrb))
    allocate(OmatRI(nOrb,nOrb))
    allocate(HlrOC(nOrb,nOrb))

    call env%globalTimer%startTimer(globalTimers%rangeSepOnsCorrH)
    call allocateAndInit(this, orb, iSquare, species, onSiteElements, overlap, densSqr,&
        & Smat, Dmat, Omat0, OmatRI)
    call evaluateHamiltonian(this, Smat, Dmat, Omat0, OmatRI, HlrOC)
    HH(:,:) = HH + HlrOC
    call env%globalTimer%stopTimer(globalTimers%rangeSepOnsCorrH)

  contains

    !> Set up storage and get orbital-by-orbital onsite constant matrix
    subroutine allocateAndInit(this, orb, iSquare, species, onSiteElements, overlap,&
        & densSqr, Smat, Dmat, Omat0, OmatRI)

      !> class instance
      type(TRangeSepOnsCorrFunc), intent(inout) :: this

      !> Atomic orbital information
      type(TOrbitals), intent(in) :: orb

      !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
      integer, dimension(:), intent(in) :: iSquare

      !> species of all atoms in the system
      integer, intent(in) :: species(:)

      !> Correction to energy from on-site matrix elements
      real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

      !> Square (unpacked) overlap matrix.
      real(dp), intent(in) :: overlap(:,:)

      !> Square (unpacked) density matrix
      real(dp), intent(in) :: densSqr(:,:)

      !> Symmetrized square overlap matrix
      real(dp), intent(out) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(out) :: Dmat(:,:)

      !> Symmetrized square onsite constant matrix except diagonal elements
      real(dp), intent(out) :: Omat0(:,:)

      !> Symmetrized square onsite constant matrix for RI loss correction
      real(dp), intent(out) :: OmatRI(:,:)

      real(dp) :: fac
      integer :: nAtom, iAtom, ii, jj, iOrb, jOrb, iSp1, iSp2

      nAtom = size(iSquare,dim=1) - 1

      ! Symmetrize overlap and density matrices
      Smat(:,:) = overlap
      call symmetrizeHS(Smat)
      Dmat(:,:) = densSqr
      call symmetrizeHS(Dmat)

      ! Set onsite constant matrices
      Omat0(:,:) = 0.0_dp
      OmatRI(:,:) = 0.0_dp
      do iAtom = 1, nAtom

        ii = iSquare(iAtom)
        jj = iSquare(iAtom) + orb%nOrbAtom(iAtom) - 1

        do iOrb = ii, jj

          if (iOrb - ii + 1 <= 1) then
            iSp1 = 1
          else if (iOrb - ii + 1 > 1 .and. iOrb - ii + 1 <= 4) then
            iSp1 = 2
          end if

          do jOrb = ii, jj

            if (jOrb - ii + 1 <= 1) then
              iSp2 = 1
            else if (jOrb - ii + 1 > 1 .and. jOrb - ii + 1 <= 4) then
              iSp2 = 2
            end if

            ! OC contribution except diagonal elements
            if (iOrb /= jOrb) then
              Omat0(iOrb,jOrb) = onSiteElements(iSp1,iSp2,3,species(iAtom))
            end if

            ! RI loss correction for p orbitals
            if (iSp1 == 2 .and. iSp2 == 2) then
              if (iOrb == jOrb) then
                fac = 2.0_dp
              else
                fac = -1.0_dp
              end if
              OmatRI(iOrb,jOrb) = fac * onSiteElements(iSp1,iSp2,3,species(iAtom))
            end if

          end do

        end do

      end do

    end subroutine allocateAndInit


    !> Evaluate the hamiltonian using GEMM operations
    subroutine evaluateHamiltonian(this, Smat, Dmat, Omat0, OmatRI, HlrOC)

      !> class instance
      type(TRangeSepOnsCorrFunc), intent(inout) :: this

      !> Symmetrized square overlap matrix
      real(dp), intent(in) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(in) :: Dmat(:,:)

      !> Symmetrized square onsite constant matrix except diagonal elements
      real(dp), intent(in) :: Omat0(:,:)

      !> Symmetrized square onsite constant matrix for RI loss correction
      real(dp), intent(in) :: OmatRI(:,:)

      !> Symmetrized long-range onsite-corrected (lrOC) Hamiltonian matrix
      real(dp), intent(out) :: HlrOC(:,:)

      real(dp), allocatable :: PS(:,:)
      real(dp), allocatable :: tmpMat(:,:)
      real(dp), allocatable :: Hmat(:,:)
      real(dp), allocatable :: tmpVec(:)
      real(dp), allocatable :: Hvec(:)

      integer :: iOrb, jOrb

      allocate(PS(nOrb,nOrb))
      allocate(tmpMat(nOrb,nOrb))
      allocate(Hmat(nOrb,nOrb))
      allocate(tmpVec(nOrb))
      allocate(Hvec(nOrb))

      call gemm(PS, Dmat, Smat)

      HlrOC(:,:) = 0.0_dp

      ! OC contribution: 1st term
      tmpMat(:,:) = Dmat * Smat
      call gemm(Hmat, Omat0, tmpMat)
      Hvec(:) = sum(Hmat,dim=2)
      do iOrb = 1, nOrb
        do jOrb = 1, nOrb
          HlrOC(iOrb,jOrb) = HlrOC(iOrb,jOrb) + Smat(iOrb,jOrb) * &
              & (Hvec(iOrb) + Hvec(jOrb))
        end do
      end do

      ! OC contribution: 2nd term
      call gemm(Hmat, Smat, PS)
      HlrOC(:,:) = HlrOC + Omat0 * Hmat

      Hmat(:,:) = Dmat * Omat0
      call gemm(tmpMat, Hmat, Smat)
      call gemm(Hmat, Smat, tmpMat)
      HlrOC(:,:) = HlrOC + Hmat

      ! OC contribution: 3rd term
      tmpMat(:,:) = PS * Omat0
      call gemm(Hmat, tmpMat, Smat)
      HlrOC(:,:) = HlrOC + Hmat

      ! OC contribution: 4th term
      tmpMat(:,:) = transpose(PS) * Omat0
      call gemm(Hmat, Smat, tmpMat)
      HlrOC(:,:) = HlrOC + Hmat

      ! OC contribution: 5th term
      call gemm(tmpMat, Smat, PS)
      tmpVec(:) = 0.0_dp
      do iOrb = 1, nOrb
        tmpVec(iOrb) = tmpMat(iOrb,iOrb)
      end do
      call gemv(Hvec, Omat0, tmpVec)
      do iOrb = 1, nOrb
       HlrOC(iOrb,iOrb) = HlrOC(iOrb,iOrb) + Hvec(iOrb)
      end do

      ! OC contribution: 6th term
      tmpVec(:) = 0.0_dp
      do iOrb = 1, nOrb
        tmpVec(iOrb) = Dmat(iOrb,iOrb)
      end do
      call gemv(Hvec, Omat0, tmpVec)
      Hmat(:,:) = 0.0_dp
      do iOrb = 1, nOrb
        Hmat(iOrb,iOrb) = Hvec(iOrb)
      end do
      call gemm(tmpMat, Hmat, Smat)
      call gemm(Hmat, Smat, tmpMat)
      HlrOC(:,:) = HlrOC + Hmat

      ! RI loss correction term
      tmpMat(:,:) = transpose(PS) * OmatRI
      call gemm(Hmat, tmpMat, Smat)
      HlrOC(:,:) = HlrOC + HMat * 2.0_dp / 3.0_dp

      call gemm(tmpMat, Smat, PS)
      Hmat(:,:) = tmpMat * OmatRI
      HlrOC(:,:) = HlrOC + HMat * 2.0_dp / 3.0_dp

      Hmat(:,:) = Dmat * OmatRI
      call gemm(tmpMat, Hmat, Smat)
      call gemm(Hmat, Smat, tmpMat)
      HlrOC(:,:) = HlrOC + HMat * 2.0_dp / 3.0_dp

      tmpMat(:,:) = PS * OmatRI
      call gemm(Hmat, Smat, tmpMat)
      HlrOC(:,:) = HlrOC + HMat * 2.0_dp / 3.0_dp

      if (this%tSpin) then
        HlrOC(:,:) = -0.25_dp * HlrOC
      else
        HlrOC(:,:) = -0.125_dp * HlrOC
      end if

    end subroutine evaluateHamiltonian

  end subroutine addLrOcHamiltonian

end module dftbp_dftb_rangeseponscorr
