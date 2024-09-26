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
  use dftbp_dftb_nonscc, only : TNonSccDiff
  use dftbp_dftb_hybridxc, only : hybridXcFunc, hybridXcAlgo, hybridXcGammaTypes
  use dftbp_dftb_slakocont, only : TSlakoCont
  use dftbp_math_matrixops, only : adjointLowerTriangle
  use dftbp_io_message, only : error
  use dftbp_math_blasroutines, only : gemm, gemv
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: TRangeSepOnsCorrFunc, RangeSepOnsCorrFunc_init


  !> Onsite correction with range-separated hybrid functional module
  type :: TRangeSepOnsCorrFunc
    private

    !> Symmetrized square onsite constant matrix except diagonal elements
    real(dp), allocatable :: Omat0(:,:)

    !> Symmetrized square onsite constant matrix for RI loss correction
    real(dp), allocatable :: OmatRI(:,:)

    !> total onsite correction energy from range-separated functional
    real(dp) :: lrOcEnergy

    !> Is this spin restricted (F) or unrestricted (T)
    logical :: tSpin

    !> Algorithm for RSH-Hamiltonian construction (and, if applicable, force calculation)
    integer :: hybridXcAlg

    !> Hybrid xc-functional type, as extracted from SK-file(s)
    integer :: hybridXcType

    !> Gamma function type (mostly for periodic cases; 'full' for non-periodic systems)
    integer :: gammaType

  contains

    procedure :: addLrOcHamiltonian
    procedure :: addLrOcEnergy
    procedure :: addLrOcGradients

  end type TRangeSepOnsCorrFunc


contains


  !> Intitialize the onsite correction with range-separated hybrid functional module
  subroutine RangeSepOnsCorrFunc_init(this, orb, iSquare, species, onSiteElements, tSpin,&
      hybridXcAlg, hybridXcType, gammaType)

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

    !> Is this spin restricted (F) or unrestricted (T)
    logical, intent(in) :: tSpin

    !> lr-Hamiltonian construction algorithm
    integer, intent(in) :: hybridXcAlg

    !> Hybrid xc-functional type, as extracted from SK-file(s)
    integer, intent(in) :: hybridXcType

    !> Gamma function type (mostly for periodic cases)
    integer, intent(in) :: gammaType

    call checkRequirements(this, orb, hybridXcAlg, hybridXcType, gammaType)
    call initAndAllocate(this, orb, iSquare, species, onSiteElements, tSpin)

  contains

    !> initialise data structures
    subroutine initAndAllocate(this, orb, iSquare, species, onSiteElements, tSpin)

      !> class instance
      class(TRangeSepOnsCorrFunc), intent(inout) :: this

      !> Atomic orbital information
      type(TOrbitals), intent(in) :: orb

      !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
      integer, dimension(:), intent(in) :: iSquare

      !> species of all atoms in the system
      integer, intent(in) :: species(:)

      !> Correction to energy from on-site matrix elements
      real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

      !> Is this spin restricted (F) or unrestricted (T)
      logical, intent(in) :: tSpin

      real(dp) :: fac
      integer :: nAtom, iAtom, nOrb, ii, jj, iOrb, jOrb, iSh1, iSh2

      nOrb = orb%nOrb
      nAtom = size(iSquare,dim=1) - 1

      this%lrOcEnergy = 0.0_dp
      this%tSpin = tSpin

      ! Set onsite constant matrices
      allocate(this%Omat0(nOrb,nOrb))
      allocate(this%OmatRI(nOrb,nOrb))
      this%Omat0(:,:) = 0.0_dp
      this%OmatRI(:,:) = 0.0_dp

      do iAtom = 1, nAtom

        ii = iSquare(iAtom)
        jj = iSquare(iAtom) + orb%nOrbAtom(iAtom) - 1

        do iOrb = ii, jj

          ! Find l-shell for iOrb
          iSh1 = orb%iShellOrb(iOrb-ii+1,species(iAtom))

          do jOrb = ii, jj

            ! Find l-shell for jOrb
            iSh2 = orb%iShellOrb(jOrb-ii+1,species(iAtom))

            ! OC contribution
            this%Omat0(iOrb,jOrb) = onSiteElements(iSh1,iSh2,3,species(iAtom))

            ! RI loss correction
            if (iSh1 == iSh2 .and. iSh1 > 1) then
              this%OmatRI(iOrb,jOrb) = onSiteElements(iSh1,iSh2,3,species(iAtom)) &
                  & / (2.0_dp * iSh1 - 1.0_dp) / 2.0_dp
            end if

          end do

        end do

      end do

    end subroutine initAndAllocate


    !> Test for option consistency
    subroutine checkRequirements(this, orb, hybridXcAlg, hybridXcType, gammaType)

      !> class instance
      class(TRangeSepOnsCorrFunc), intent(inout) :: this

      !> Atomic orbital information
      type(TOrbitals), intent(in) :: orb

      !> lr-Hamiltonian construction algorithm
      integer, intent(in) :: hybridXcAlg

      !> Hybrid xc-functional type, as extracted from SK-file(s)
      integer, intent(in) :: hybridXcType

      !> Gamma function type (mostly for periodic cases)
      integer, intent(in) :: gammaType

      this%hybridXcAlg = hybridXcAlg
      this%hybridXcType = hybridXcType
      this%gammaType = gammaType

      ! Check for current restrictions
      if (orb%mOrb > 4) then
        call error("Onsite correction with range separated functional only works with&
            & the system consisting of up to p orbitals")
      end if

      if (.not. this%hybridXcType == hybridXcFunc%lc) then
        call error("Long-range corrected functional only works with onsite corrections")
      end if

      if (.not. this%hybridXcAlg == hybridXcAlgo%matrixBased) then
        call error("Onsite correction with range separated functional only works with&
            & matrix based algorithm")
      end if

      if (.not. this%gammaType == hybridXcGammaTypes%full) then
        call error("Full, unaltered gamma function only works with onsite corrections")
      end if

    end subroutine checkRequirements

  end subroutine RangeSepOnsCorrFunc_init


  !> Interface routine.
  subroutine addLrOcHamiltonian(this, env, overlap, densSqr, HH)

    !> class instance
    class(TRangeSepOnsCorrFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Square (unpacked) overlap matrix.
    real(dp), intent(in) :: overlap(:,:)

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: densSqr(:,:)

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout) :: HH(:,:)

    call env%globalTimer%startTimer(globalTimers%rangeSepOnsCorrH)
    select case(this%hybridXcAlg)
    case (hybridXcAlgo%thresholdBased)
      ! not supported at the moment
    case (hybridXcAlgo%neighbourBased)
      ! not supported at the moment
    case (hybridXcAlgo%matrixBased)
      call addLrOcHamiltonianMatrix(this, overlap, densSqr, HH)
    end select
    call env%globalTimer%stopTimer(globalTimers%rangeSepOnsCorrH)

  end subroutine addLrOcHamiltonian


  !> Update Hamiltonian with long-range onsite contribution using matrix-matrix multiplications
  !>
  !> The routine provides a matrix-matrix multiplication based implementation of
  !> Eq. 11 in https://doi.org/10.1021/acs.jctc.2c00037
  !>
  subroutine addLrOcHamiltonianMatrix(this, overlap, densSqr, HH)

    !> class instance
    class(TRangeSepOnsCorrFunc), intent(inout) :: this

    !> Square (unpacked) overlap matrix.
    real(dp), intent(in) :: overlap(:,:)

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: densSqr(:,:)

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout) :: HH(:,:)

    real(dp), allocatable :: Smat(:,:)
    real(dp), allocatable :: Dmat(:,:)
    real(dp), allocatable :: HlrOC(:,:)

    call allocateAndInit(this, overlap, densSqr, Smat, Dmat, HlrOC)
    call evaluateHamiltonian(this, Smat, Dmat, HlrOC)
    HH(:,:) = HH + HlrOC
    this%lrOcEnergy = this%lrOcEnergy + 0.5_dp * sum(Dmat * HlrOC)

  contains

    !> allocate and initialise some necessary arrays
    subroutine allocateAndInit(this, overlap, densSqr, Smat, Dmat, HlrOC)

      !> class instance
      type(TRangeSepOnsCorrFunc), intent(inout) :: this

      !> Square (unpacked) overlap matrix.
      real(dp), intent(in) :: overlap(:,:)

      !> Square (unpacked) density matrix
      real(dp), intent(in) :: densSqr(:,:)

      !> Symmetrized square overlap matrix
      real(dp), allocatable, intent(inout) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), allocatable, intent(inout) :: Dmat(:,:)

      !> Symmetrized long-range onsite-corrected (lrOC) Hamiltonian matrix
      real(dp), allocatable, intent(inout) :: HlrOC(:,:)

      integer :: nOrb

      nOrb = size(overlap,dim=1)

      allocate(Smat(nOrb,nOrb))
      allocate(Dmat(nOrb,nOrb))
      allocate(HlrOC(nOrb,nOrb))

      ! Symmetrize overlap and density matrices
      Smat(:,:) = overlap
      call adjointLowerTriangle(Smat)
      Dmat(:,:) = densSqr
      call adjointLowerTriangle(Dmat)

    end subroutine allocateAndInit


    !> Evaluate the hamiltonian using GEMM operations
    subroutine evaluateHamiltonian(this, Smat, Dmat, HlrOC)

      !> class instance
      type(TRangeSepOnsCorrFunc), intent(inout) :: this

      !> Symmetrized square overlap matrix
      real(dp), intent(in) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(in) :: Dmat(:,:)

      !> Symmetrized long-range onsite-corrected (lrOC) Hamiltonian matrix
      real(dp), intent(out) :: HlrOC(:,:)

      real(dp), allocatable :: PS(:,:)
      real(dp), allocatable :: tmpMat(:,:)
      real(dp), allocatable :: Hmat(:,:)
      real(dp), allocatable :: tmpVec(:)
      real(dp), allocatable :: Hvec(:)

      real(dp) :: fac, facRI
      integer :: nOrb, iOrb, jOrb

      nOrb = size(Smat,dim=1)

      allocate(PS(nOrb,nOrb))
      allocate(tmpMat(nOrb,nOrb))
      allocate(Hmat(nOrb,nOrb))
      allocate(tmpVec(nOrb))
      allocate(Hvec(nOrb))

      if (this%tSpin) then
        fac = -0.25_dp
        facRI = 1.0_dp
      else
        fac = -0.125_dp
        facRI = 0.5_dp
      end if

      call gemm(PS, Dmat, Smat)

      HlrOC(:,:) = 0.0_dp

      ! OC contribution: external overlap term
      tmpMat(:,:) = Dmat * Smat
      tmpVec(:) = sum(tmpMat,dim=2)
      call gemv(Hvec, this%Omat0, tmpVec)
      do iOrb = 1, nOrb
        do jOrb = 1, nOrb
          HlrOC(iOrb,jOrb) = HlrOC(iOrb,jOrb) + fac * Smat(iOrb,jOrb) * &
              & (Hvec(iOrb) + Hvec(jOrb))
        end do
      end do

      ! OC & RI contributions: S(P*O)S + SPS*O
      tmpMat(:,:) = fac * Dmat * this%Omat0 + facRI * Dmat * this%OmatRI
      call gemm(Hmat, tmpMat, Smat)
      call gemm(HlrOC, Smat, Hmat, beta=1.0_dp)

      call gemm(tmpMat, Smat, PS)
      HlrOC(:,:) = HlrOC + tmpMat * (fac * this%Omat0 + facRI * this%OmatRI)

      ! OC & RI contributions: S(PS*O) and its transpose part
      tmpMat(:,:) = fac * PS * this%Omat0 + facRI * transpose(PS) * this%OmatRI
      call gemm(HlrOC, tmpMat, Smat, beta=1.0_dp)

      call gemm(HlrOC, Smat, tmpMat, transB="T", beta=1.0_dp)

      ! OC contribution: diagonal element for Hamiltonian
      tmpMat(:,:) = PS * Smat
      tmpVec(:) = sum(tmpMat,dim=1)
      call gemv(Hvec, this%Omat0, tmpVec)
      do iOrb = 1, nOrb
        HlrOC(iOrb,iOrb) = HlrOC(iOrb,iOrb) + fac * Hvec(iOrb)
      end do

      ! OC contribution: diagonal element for density
      tmpVec(:) = 0.0_dp
      do iOrb = 1, nOrb
        tmpVec(iOrb) = Dmat(iOrb,iOrb)
      end do
      call gemv(Hvec, this%Omat0, tmpVec)
      Hmat(:,:) = 0.0_dp
      do iOrb = 1, nOrb
        Hmat(iOrb,iOrb) = Hvec(iOrb)
      end do
      call gemm(tmpMat, Hmat, Smat)
      call gemm(HlrOC, Smat, tmpMat, alpha=fac, beta=1.0_dp)

    end subroutine evaluateHamiltonian

  end subroutine addLrOcHamiltonianMatrix


  !> Add the onsite contribution originating from range-seprated functional to the total energy
  subroutine addLrOcEnergy(this, energy)

    !> class instance
    class(TRangeSepOnsCorrFunc), intent(inout) :: this

    !> total energy
    real(dp), intent(inout) :: energy

    energy = energy + this%lrOcEnergy
    ! hack for spin unrestricted calculation
    this%lrOcEnergy = 0.0_dp

  end subroutine addLrOcEnergy


  !> Interface routine.
  subroutine addLrOcGradients(this, gradients, derivator, skOverCont, coords, nNeighbourSK,&
      & iNeighbour, iSquare, species, orb, densSqr, overlap)

    !> class instance
    class(TRangeSepOnsCorrFunc), intent(inout) :: this

    !> energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> number of atoms neighbouring each site where the overlap is non-zero
    integer, intent(in) :: nNeighbourSK(:)

    !> neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, dimension(:), intent(in) :: iSquare

    !> species of all atoms
    integer, intent(in) :: species(:)

    !> orbital information
    type(TOrbitals), intent(in) :: orb

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: densSqr(:,:,:)

    !> square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    select case(this%hybridXcAlg)
    case (hybridXcAlgo%thresholdBased)
      ! not supported at the moment
    case (hybridXcAlgo%neighbourBased)
      ! not supported at the moment
    case (hybridXcAlgo%matrixBased)
      call addLrOcGradientsMatrix(this, gradients, derivator, skOverCont, coords, nNeighbourSK,&
          & iNeighbour, iSquare, species, orb, densSqr, overlap)
    end select

  end subroutine addLrOcGradients


  !> Update gradients with long-range onsite contribution
  !>
  !> The routine provides a matrix-matrix multiplication based implementation of
  !> Eq. 29 in https://doi.org/10.1021/acs.jctc.2c00037
  !>
  subroutine addLrOcGradientsMatrix(this, gradients, derivator, skOverCont, coords, nNeighbourSK,&
      & iNeighbour, iSquare, species, orb, densSqr, overlap)

    !> class instance
    class(TRangeSepOnsCorrFunc), intent(inout) :: this

    !> energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> number of atoms neighbouring each site where the overlap is non-zero
    integer, intent(in) :: nNeighbourSK(:)

    !> neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, dimension(:), intent(in) :: iSquare

    !> species of all atoms
    integer, intent(in) :: species(:)

    !> orbital information
    type(TOrbitals), intent(in) :: orb

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: densSqr(:,:,:)

    !> square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    real(dp), allocatable :: Smat(:,:)
    real(dp), allocatable :: Dmat(:,:,:)

    real(dp), allocatable :: PS(:,:,:)
    real(dp), allocatable :: tmpMat(:,:)
    real(dp), allocatable :: Hmat(:,:)
    real(dp), allocatable :: tmpVec(:)
    real(dp), allocatable :: Hvec(:)

    real(dp), allocatable :: shiftSqr(:,:,:)
    real(dp), allocatable :: sPrimeTmp(:,:,:)
    real(dp), allocatable :: shift(:,:)
    real(dp), allocatable :: tmpForce(:)

    real(dp) :: fac, facRI
    integer :: nOrb, nSpin, nAtom
    integer :: iOrb, mu, nu, iSpin, iAtA, iAtB, iNeighA, ii
    integer :: iAtMu1, iAtMu2, iAtNu1, iAtNu2

    nOrb = size(overlap,dim=1)
    nSpin = size(densSqr,dim=3)
    nAtom = size(iSquare,dim=1) - 1

    allocate(PS(nOrb,nOrb,nSpin))
    allocate(tmpMat(nOrb,nOrb))
    allocate(Hmat(nOrb,nOrb))
    allocate(tmpVec(nOrb))
    allocate(Hvec(nOrb))

    allocate(shiftSqr(nOrb,nOrb,nSpin))
    allocate(sPrimeTmp(orb%mOrb,orb%mOrb,3))
    allocate(shift(orb%mOrb,orb%mOrb))
    allocate(tmpForce(3))

    ! Initialize several matrices for calculation of gradients
    call allocateAndInit(this, overlap, densSqr, Smat, Dmat)

    if (this%tSpin) then
      fac = -0.25_dp
      facRI = 1.0_dp
    else
      fac = -0.125_dp
      facRI = 0.5_dp
    end if

    do iSpin = 1, nSpin
      call gemm(PS(:,:,iSpin), Dmat(:,:,iSpin), Smat)
    end do

    ! shift values
    shiftSqr(:,:,:) = 0.0_dp

    do iSpin = 1, nSpin

      ! OC contribution: diagonal element for PS
      tmpVec(:) = 0.0_dp
      do iOrb = 1, nOrb
        tmpVec(iOrb) = PS(iOrb,iOrb,iSpin)
      end do
      call gemv(Hvec, this%Omat0, tmpVec, alpha=fac)
      do mu = 1, nOrb
        do nu = 1, nOrb
          shiftSqr(mu,nu,iSpin) = shiftSqr(mu,nu,iSpin)&
              & + Dmat(mu,nu,iSpin) * (Hvec(mu) + Hvec(nu))
        end do
      end do

      ! OC contribution: diagonal element for density
      tmpVec(:) = 0.0_dp
      do iOrb = 1, nOrb
        tmpVec(iOrb) = Dmat(iOrb,iOrb,iSpin)
      end do
      call gemv(Hvec, this%Omat0, tmpVec, alpha=fac)
      do mu = 1, nOrb
        do nu = 1, nOrb
          shiftSqr(mu,nu,iSpin) = shiftSqr(mu,nu,iSpin)&
              & + (PS(mu,nu,iSpin) * Hvec(nu) + PS(nu,mu,iSpin) * Hvec(mu))
        end do
      end do

      ! OC & RI contributions: PS(P*O) + (P*O)SP
      tmpMat(:,:) = fac * Dmat(:,:,iSpin) * this%Omat0 + &
          & facRI * Dmat(:,:,iSpin) * this%OmatRI
      call gemm(Hmat, PS(:,:,iSpin), tmpMat)
      shiftSqr(:,:,iSpin) = shiftSqr(:,:,iSpin) + (transpose(Hmat) + Hmat)

      ! OC contribution: (SP*O)P + P(PS*O)
      ! RI contribution: (PS*O)P + P(SP*O)
      tmpMat(:,:) = fac * PS(:,:,iSpin) * this%Omat0 + &
          & facRI * transpose(PS(:,:,iSpin)) * this%OmatRI
      call gemm(Hmat, Dmat(:,:,iSpin), tmpMat)
      shiftSqr(:,:,iSpin) = shiftSqr(:,:,iSpin) + (transpose(Hmat) + Hmat)

    end do

    ! Compute gradients originating from onsite contribution with range separated hybrid functional

    ! sum A
    loopA: do iAtA = 1, nAtom
      ! sum B
      loopB: do iNeighA = 1, nNeighbourSK(iAtA)
        iAtB = iNeighbour(iNeighA, iAtA)

        ! evaluate the ovr_prime
        sPrimeTmp(:,:,:) = 0.0_dp
        if ( iAtA /= iAtB ) then

          call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, species, iAtA, iAtB, orb)

          ! orbital index for atom A and B
          iAtMu1 = iSquare(iAtA)
          iAtMu2 = iSquare(iAtA+1) - 1
          iAtNu1 = iSquare(iAtB)
          iAtNu2 = iSquare(iAtB+1) - 1

          ! shift values for gradient
          shift(:,:) = 0.0_dp
          do iSpin = 1, nSpin
            do mu = iAtMu1, iAtMu2
              do nu = iAtNu1, iAtNu2
                shift(nu-iAtNu1+1,mu-iAtMu1+1) = shift(nu-iAtNu1+1,mu-iAtMu1+1)&
                    & + shiftSqr(mu,nu,iSpin) + shiftSqr(nu,mu,iSpin)
              end do
            end do
          end do

          tmpForce(:) = 0.0_dp
          do ii = 1, 3
            tmpForce(ii) = sum(shift(:,:)*sPrimeTmp(:,:,ii))
          end do

          ! forces from atom A on atom B and B onto A
          gradients(:,iAtA) = gradients(:,iAtA) + tmpForce
          gradients(:,iAtB) = gradients(:,iAtB) - tmpForce

        end if

      end do loopB
    end do loopA

  contains

    !> allocate and initialise some necessary arrays
    subroutine allocateAndInit(this, overlap, densSqr, Smat, Dmat)

      !> class instance
      type(TRangeSepOnsCorrFunc), intent(inout) :: this

      !> Square (unpacked) overlap matrix.
      real(dp), intent(in) :: overlap(:,:)

      !> Square (unpacked) density matrix
      real(dp), intent(in) :: densSqr(:,:,:)

      !> Symmetrized square overlap matrix
      real(dp), allocatable, intent(inout) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), allocatable, intent(inout) :: Dmat(:,:,:)

      integer :: iSpin

      allocate(Smat(nOrb,nOrb))
      allocate(Dmat(nOrb,nOrb,nSpin))

      ! Symmetrize overlap and density matrices
      Smat(:,:) = overlap
      call adjointLowerTriangle(Smat)
      Dmat(:,:,:) = densSqr
      do iSpin = 1, nSpin
        call adjointLowerTriangle(Dmat(:,:,iSpin))
      end do

    end subroutine allocateAndInit

  end subroutine addLrOcGradientsMatrix

end module dftbp_dftb_rangeseponscorr
