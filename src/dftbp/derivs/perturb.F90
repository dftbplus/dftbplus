!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Module for linear response derivative calculations using perturbation methods at q=0
module dftbp_derivs_perturb
  use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
  use dftbp_common_accuracy, only : dp, lc, mc
  use dftbp_common_constants, only : AA__Bohr, Bohr__AA, Hartree__eV, quaternionName
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_file, only : closeFile, openFile, TFileDescr
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_status, only : TStatus
  use dftbp_derivs_fermihelper, only : deltamn
  use dftbp_derivs_linearresponse, only : dRhoCmplx, dRhoFermiChangeCmplx, dRhoFermiChangePauli,&
      & dRhoFermiChangeReal, dRhoPauli, dRhoReal
  use dftbp_derivs_rotatedegen, only : TRotateDegen, TRotateDegen_init
  use dftbp_dftb_blockpothelper, only : appendBlockReduced
  use dftbp_dftb_boundarycond, only : boundaryCondsEnum, TBoundaryConds
  use dftbp_dftb_coulomb, only : calcInvRPrimeAsymm
  use dftbp_dftb_dftbplusu, only : plusUFunctionals, TDftbU
  use dftbp_dftb_hybridxc, only : THybridXcFunc
  use dftbp_dftb_nonscc, only : TNonSccDiff
  use dftbp_dftb_onsitecorrection, only : addOnsShift, onsblock_expand
  use dftbp_dftb_orbitalequiv, only : OrbitalEquiv_expand, OrbitalEquiv_reduce
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_populations, only : denseMullikenReal, getChargePerShell, getOnsitePopulation,&
      & mulliken
  use dftbp_dftb_potentials, only : TPotentials, TPotentials_init
  use dftbp_dftb_scc, only : TScc
  use dftbp_dftb_shift, only : addShift, totalShift
  use dftbp_dftb_slakocont, only : TSlakoCont
  use dftbp_dftb_spin, only : getSpinShift, qm2ud, ud2qm
  use dftbp_dftb_thirdorder, only : TThirdOrder
  use dftbp_dftbplus_mainio, only : writeDerivBandOut, permitivityPrint
  use dftbp_dftbplus_outputfiles, only : derivVBandOut
  use dftbp_io_commonformats, only : format2U
  use dftbp_io_message, only : warning
  use dftbp_io_taggedoutput, only : tagLabels, TTaggedWriter
  use dftbp_mixer_mixer, only : TMixerReal
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_parallelks, only : TParallelKS
#:if WITH_SCALAPACK
  use dftbp_dftb_populations, only : denseMullikenRealBlacs
  use dftbp_extlibs_mpifx, only : mpifx_allreduceip, mpifx_bcast, MPI_SUM
  use dftbp_extlibs_scalapackfx, only : DLEN_, blocklist, scalafx_getdescriptor, size
#:endif
  use dftbp_dftb_sparse2dense, only : unpackHS
  implicit none

  private
  public :: TPerturbInp, TResponse, TResponse_init, responseSolverTypes


  !> Input type for perturbation calculations
  type :: TPerturbInp

    !> Is this is a static electric field perturbation calculation
    logical :: isEPerturb = .false.

    !> Is the response kernel (and frontier eigenvalue derivatives) calculated by perturbation
    logical :: isRespKernelPert = .false.

    !> Are derivatives with respect to atomic positions required
    logical :: isAtomCoordPerturb = .false.

    !> List of atoms with which to differentiate (optional, if not allocated, all of them)
    integer, allocatable :: indWrtAtoms(:)

    !> Are this derivatives with respect to external charges required?
    logical :: isExtChargeDeriv = .false.

    !> List of external charges with which to differentiate (optional, if not allocated, all of
    !! them)
    integer, allocatable :: indWrtCharges(:)

    !> Tolerance, as a multiple of eps, for idenfifying need for degenerate perturbation theory
    real(dp) :: tolDegenDFTBPT = 128.0_dp

    !> Frequencies for perturbation (0 being static case)
    real(dp), allocatable :: dynEFreq(:)

    !> Frequency dependent perturbation eta
    real(dp), allocatable :: etaFreq

    !> Is the response kernel evaluated at the RPA level, or (if SCC) self-consistent
    logical :: isRespKernelRPA = .false.

    !> Frequencies for perturbation (0 being static case)
    real(dp), allocatable :: dynKernelFreq(:)

    !> Self-consistency tolerance for perturbation (if SCC)
    real(dp) :: perturbSccTol = 0.0_dp

    !> Maximum iterations for perturbation theory
    integer :: maxPerturbIter = 1

    !> Require converged perturbation (if true, terminate on failure, otherwise return NaN for
    !! non-converged)
    logical :: isPerturbConvRequired = .true.

  end type TPerturbInp


  !> Namespace for possible perturbation solver methods
  type :: TPerturbSolverTypesEnum

    !> Coupled perturbed equations solved by sum over unperturbed spectra
    integer :: spectralSum = 1

  end type TPerturbSolverTypesEnum


  !> Actual values for perturbSolverTypes.
  type(TPerturbSolverTypesEnum), parameter :: responseSolverTypes = TPerturbSolverTypesEnum()


  !> Response solver
  type TResponse

    !> Solver choice for response evaluation
    integer :: iSolver

    !> Tolerance for degeneracy between eigenvalues
    real(dp) :: tolDegen

    !> Small complex value for frequency dependent cases
    complex(dp), allocatable :: eta

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical :: isEfFixed

  contains

    !> Perturbation with respect to electric field
    procedure :: wrtEField

    !> Perturbation with respect to potentials at atomic sites
    procedure :: wrtVAtom

    !> Derivative with respect to external charges
    procedure :: dxExtCharges

    !> Derivative with respect to atomic displacements
    procedure :: dxAtom

  end type TResponse

  !> Direction labels
  character(len=1), parameter :: direction(3) = ['x','y','z']

contains

  !> Initialise response solver(s)
  subroutine TResponse_init(this, iSolver, isEfFixed, tolDegen, eta)

    !> Type to initialise
    type(TResponse), intent(out) :: this

    !> Choice of method to solve equations
    integer, intent(in) :: iSolver

    !> Is a fixed Fermi energy used
    logical, intent(in) :: isEfFixed

    !> Tolerance for using degenerate perturbation of eigen-states
    real(dp), intent(in) :: tolDegen

    !> Small constant for imaginary part of finite frequency calculations
    real(dp), intent(in), optional :: eta

    this%iSolver = iSolver

    this%tolDegen = tolDegen

    this%isEfFixed = isEfFixed

    if (present(eta)) then
      allocate(this%eta)
      this%eta = cmplx(0.0_dp, eta, dp)
    end if

  end subroutine TResponse_init


  !> Perturbation at q=0 with respect to a static electric field
  subroutine wrtEField(this, env, parallelKS, filling, eigvals, eigVecsReal, eigVecsCplx, ham,&
      & over, boundaryCond, orb, nAtom, species, neighbourList, nNeighbourSK, denseDesc,&
      & iSparseStart, img2CentCell, coord, coord0, sccCalc, maxSccIter, sccTol, isSccConvRequired,&
      & nMixElements, nIneqMixElements, iEqOrbitals, tempElec, Ef, spinW, thirdOrd, dftbU,&
      & iEqBlockDftbu, onsMEs, iEqBlockOnSite, hybridXc, nNeighbourCam, chrgMixerReal, kPoint,&
      & kWeight, iCellVec, cellVec, polarisability, dEi, dqOut, neFermi, dEfdE, errStatus, omega)

    !> Instance
    class(TResponse), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> The k-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Filling
    real(dp), intent(in) :: filling(:,:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in), allocatable :: eigVecsReal(:,:,:)

    !> Ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> Sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Boundary conditions on the calculation
    type(TBoundaryConds), intent(in) :: boundaryCond

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> Chemical species
    integer, intent(in) :: species(:)

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Atomic coordinates in central cell
    real(dp), intent(in) :: coord0(:,:)

    !> SCC module internal variables
    type(TScc), intent(inout), allocatable :: sccCalc

    !> Maximal number of SCC iterations
    integer, intent(in) :: maxSccIter

    !> Tolerance for SCC convergence
    real(dp), intent(in) :: sccTol

    !> Use converged derivatives of charges
    logical, intent(in) :: isSccConvRequired

    !> Nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
    !> blocks (if a DFTB+U or onsite correction calculation)
    integer, intent(in) :: nMixElements

    !> Nr. of inequivalent charges
    integer, intent(in) :: nIneqMixElements

    !> Equivalence relations between orbitals
    integer, intent(in), allocatable :: iEqOrbitals(:,:,:)

    !> Onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> Third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> Equivalence mapping for dual charge blocks
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> Data for hybrid functional calculation
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Number of neighbours for each of the atoms for the exchange contributions of CAM functionals
    integer, intent(in), allocatable :: nNeighbourCam(:)

    !> Charge mixing object
    class(TMixerReal), intent(inout), allocatable :: chrgMixerReal

    !> The k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Electric polarisability
    real(dp), allocatable, intent(inout) :: polarisability(:,:,:)

    !> Derivatives of eigenvalues, if required
    real(dp), allocatable, intent(inout) :: dEi(:,:,:,:)

    !> Derivatives of Mulliken charges, if required
    real(dp), allocatable, intent(inout) :: dqOut(:,:,:,:)

    !> Number of electrons at the Fermi energy (if metallic)
    real(dp), allocatable, intent(inout) :: neFermi(:)

    !> Derivative of the Fermi energy (if metallic)
    real(dp), allocatable, intent(out) :: dEfdE(:,:)

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    !> Driving frequencies (including potentially 0 for static)
    real(dp), intent(in) :: omega(:)

    integer :: iAt, iCart
    integer :: nSpin, nKpts, nOrbs, nIndepHam

    ! maximum allowed number of electrons in a single particle state
    real(dp) :: maxFill

    integer, allocatable :: nFilled(:,:), nEmpty(:,:)
    real(dp), allocatable :: dEiTmp(:,:,:), dEfdETmp(:)

    integer :: iOmega

    ! matrices for derivatives of terms in hamiltonian and outputs
    real(dp), allocatable :: dHam(:,:), idHam(:,:), dRho(:,:), idRho(:,:)
    real(dp) :: dqIn(orb%mOrb,nAtom, size(ham, dim=2))

    real(dp), allocatable :: dqBlockIn(:,:,:,:), SSqrReal(:,:)
    real(dp), allocatable :: dqBlockOut(:,:,:,:)

    ! derivative of potentials
    type(TPotentials) :: dPotential

    logical :: tSccCalc
    logical, allocatable :: tMetallic(:,:)

    real(dp), allocatable :: dPsiReal(:,:,:)
    complex(dp), allocatable :: dPsiCmplx(:,:,:,:)

    ! variables used for hybrid functional contributions, note this stays in the up/down
    ! representation throughout if spin polarised
    real(dp), pointer :: dRhoOutSqr(:,:,:), dRhoInSqr(:,:,:)
    real(dp), allocatable, target :: dRhoOut(:), dRhoIn(:)

    !> For transformation in the  case of degeneracies
    type(TRotateDegen), allocatable :: degenTransform(:)

    write(stdOut,*)
    write(stdOut,*)'Perturbation calculation with respect to applied electric field'
    write(stdOut,*)

    call init_perturbation(parallelKS, this%tolDegen, nOrbs, nKpts, nSpin, nIndepHam, maxFill,&
        & filling, ham, nFilled, nEmpty, dHam, dRho, idHam, idRho, degenTransform, hybridXc,&
        & sSqrReal, over, neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell,&
        & dRhoOut, dRhoIn, dRhoInSqr, dRhoOutSqr, dPotential, orb, nAtom, tMetallic, neFermi,&
        & eigvals, tempElec, Ef, kWeight)

    tSccCalc = allocated(sccCalc)

    if (allocated(dftbU) .or. allocated(onsMEs)) then
      allocate(dqBlockIn(orb%mOrb,orb%mOrb,nAtom,nSpin))
      allocate(dqBlockOut(orb%mOrb,orb%mOrb,nAtom,nSpin))
    end if

    call TPotentials_init(dPotential,orb,nAtom,nSpin,0,0)

    ! If derivatives of eigenvalues are needed
    if (allocated(dEi)) then
      dEi(:,:,:,:) = 0.0_dp
      allocate(dEiTmp(nOrbs, nKPts, nIndepHam))
    end if

    if (any(tMetallic)) then
      allocate(dEfdE(nIndepHam, 3))
      allocate(dEfdETmp(nIndepHam))
    end if

    ! if derivatives of valence wavefunctions needed. Note these will have an arbitrary set of
    ! global phases
    ! if (allocated(eigVecsReal)) then
    !   allocate(dPsiReal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2), nIndepHam, 3))
    ! else
    !   allocate(dPsiCmplx(size(eigvecsCplx,dim=1), size(eigvecsCplx,dim=2), nKpts, nIndepHam, 3))
    ! end if

    dqOut(:,:,:,:) = 0.0_dp

    ! Electric field polarisation direction
    ! note: could MPI parallelise over this in principle
    lpCart: do iCart = 1, 3

      if (boundaryCond%iBoundaryCondition == boundaryCondsEnum%cluster) then
        write(stdOut,*)"Polarisabilty for field along ", trim(quaternionName(iCart+1))
      end if

      ! set outside loop, as in time dependent case if adjacent frequencies are similar this should
      ! converge a bit faster
      dqIn(:,:,:) = 0.0_dp
      if (allocated(dftbU) .or. allocated(onsMEs)) then
        dqBlockIn(:,:,:,:) = 0.0_dp
        dqBlockOut(:,:,:,:) = 0.0_dp
      end if

      dPotential%extAtom(:,:) = 0.0_dp
      dPotential%extShell(:,:,:) = 0.0_dp
      dPotential%extBlock(:,:,:,:) = 0.0_dp

      ! derivative wrt to electric field as a perturbation
      do iAt = 1, nAtom
        dPotential%extAtom(iAt,1) = coord0(iCart,iAt)
      end do
      call totalShift(dPotential%extShell, dPotential%extAtom, orb, species)
      call totalShift(dPotential%extBlock, dPotential%extShell, orb, species)

      do iOmega = 1, size(omega)

        if (allocated(dEfdETmp)) then
          dEfdETmp(:) = 0.0_dp
        end if

        dEiTmp(:,:,:) = 0.0_dp
        if (allocated(polarisability)) then
          call response(env, parallelKS, dPotential, nAtom, orb, species, neighbourList,&
              & nNeighbourSK, img2CentCell, iSparseStart, denseDesc, over, iEqOrbitals, sccCalc,&
              & sccTol, isSccConvRequired, maxSccIter, chrgMixerReal, nMixElements,&
              & nIneqMixElements, dqIn, dqOut(:,:,:,iCart), hybridXc, nNeighbourCam, sSqrReal,&
              & dRhoInSqr, dRhoOutSqr, dRhoIn, dRhoOut, nSpin, maxFill, spinW, thirdOrd, dftbU,&
              & iEqBlockDftbu, onsMEs, iEqBlockOnSite, dqBlockIn, dqBlockOut, eigVals,&
              & degenTransform, dEiTmp, dEfdETmp, filling, Ef, this%isEfFixed, dHam, idHam, dRho,&
              & idRho, tempElec, tMetallic, neFermi, nFilled, nEmpty, kPoint, kWeight, cellVec,&
              & iCellVec, eigVecsReal, eigVecsCplx, dPsiReal, dPsiCmplx, coord, errStatus,&
              & omega(iOmega), dDipole=polarisability(:, iCart, iOmega), eta=this%eta)
        else
          call response(env, parallelKS, dPotential, nAtom, orb, species, neighbourList,&
              & nNeighbourSK, img2CentCell, iSparseStart, denseDesc, over, iEqOrbitals, sccCalc,&
              & sccTol, isSccConvRequired, maxSccIter, chrgMixerReal, nMixElements,&
              & nIneqMixElements, dqIn, dqOut(:,:,:,iCart), hybridXc, nNeighbourCam, sSqrReal,&
              & dRhoInSqr, dRhoOutSqr, dRhoIn, dRhoOut, nSpin, maxFill, spinW, thirdOrd, dftbU,&
              & iEqBlockDftbu, onsMEs, iEqBlockOnSite, dqBlockIn, dqBlockOut, eigVals,&
              & degenTransform, dEiTmp, dEfdETmp, filling, Ef, this%isEfFixed, dHam, idHam, dRho,&
              & idRho, tempElec, tMetallic, neFermi, nFilled, nEmpty, kPoint, kWeight, cellVec,&
              & iCellVec, eigVecsReal, eigVecsCplx, dPsiReal, dPsiCmplx, coord, errStatus,&
              & omega(iOmega), eta=this%eta)
        end if
        @:PROPAGATE_ERROR(errStatus)

        if (allocated(dEfdE)) then
          dEfdE(:,iCart) = dEfdETmp
        end if
        if (allocated(dEiTmp)) then
          dEi(:,:,:,iCart) = dEiTmp
        end if

        if (any(tMetallic)) then
          write(stdOut,*)
          write(stdOut,"(A,2E20.12)")'d E_f / d E_'//trim(quaternionName(iCart+1))//':',&
              & dEfdE(:,iCart)
          write(stdOut,*)
        end if

      end do

    end do lpCart

  #:if WITH_SCALAPACK
    if (allocated(dEi)) then
      call mpifx_allreduceip(env%mpi%globalComm, dEi, MPI_SUM)
    end if
  #:endif

    if (boundaryCond%iBoundaryCondition == boundaryCondsEnum%cluster) then
      call permitivityPrint(stdOut, polarisability, omega)
    end if

  end subroutine wrtEField


  !> Response with respect to a potential at atomic sites
  subroutine wrtVAtom(this, env, parallelKS, isAutotestWritten, autotestTagFile,&
      & isTagResultsWritten, resultsTagFile, taggedWriter, isBandWritten, fdDetailedOut, filling,&
      & eigvals, eigVecsReal, eigVecsCplx, ham, over, orb, nAtom, species, neighbourList,&
      & nNeighbourSK, denseDesc, iSparseStart, img2CentCell, isRespKernelRPA, sccCalc, maxSccIter,&
      & sccTol, isSccConvRequired, nMixElements, nIneqMixElements, iEqOrbitals, tempElec, Ef,&
      & spinW, thirdOrd, dftbU, iEqBlockDftbu, onsMEs, iEqBlockOnSite, hybridXc, nNeighbourCam,&
      & chrgMixerReal, kPoint, kWeight, iCellVec, cellVec, nEFermi, errStatus, omega,&
      & isHelical, coord)

    !> Instance
    class(TResponse), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> The k-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Should regression test data be written
    logical, intent(in) :: isAutotestWritten

    !> Name of output file
    character(*), intent(in) :: autotestTagFile

    !> Produce machine readable results
    logical, intent(in) :: isTagResultsWritten

    !> Name of output file
    character(*), intent(in) :: resultsTagFile

    !> Tagged writer object
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Should eigenvalue (band) data derivatives be written to disc
    logical, intent(in) :: isBandWritten

    !> File descriptor for the human readable output
    type(TFileDescr), intent(in) :: fdDetailedOut

    !> Filling
    real(dp), intent(in) :: filling(:,:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in), allocatable :: eigVecsReal(:,:,:)

    !> Ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> Sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> Chemical species
    integer, intent(in) :: species(:)

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> SCC module internal variables
    type(TScc), intent(inout), allocatable :: sccCalc

    !> Should the kernel be evaluated at the RPA level (non-SCC) or self-consistent
    logical, intent(in) :: isRespKernelRPA

    !> Maximal number of SCC iterations
    integer, intent(in) :: maxSccIter

    !> Tolerance for SCC convergence
    real(dp), intent(in) :: sccTol

    !> Use converged derivatives of charges
    logical, intent(in) :: isSccConvRequired

    !> Nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
    !> blocks (if a DFTB+U or onsite correction calculation)
    integer, intent(in) :: nMixElements

    !> Nr. of inequivalent charges
    integer, intent(in) :: nIneqMixElements

    !> Equivalence relations between orbitals
    integer, intent(in), allocatable :: iEqOrbitals(:,:,:)

    !> Onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> Third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> Equivalence mapping for dual charge blocks
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> Data for hybrid functional calculation
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Number of neighbours for each of the atoms for the exchange contributions of CAM functionals
    integer, intent(in), allocatable :: nNeighbourCam(:)

    !> Charge mixing object
    class(TMixerReal), intent(inout), allocatable :: chrgMixerReal

    !> The k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Number of electrons at the Fermi energy (if metallic)
    real(dp), allocatable, intent(inout) :: neFermi(:)

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    !> Driving frequencies (including potentially 0 for static)
    real(dp), intent(in) :: omega(:)

    !> Is the geometry helical
    logical, intent(in), optional :: isHelical

    !> Coordinates of all atoms including images
    real(dp), intent(in), optional :: coord(:,:)

    integer :: iS, iK, iAt, jAt
    integer :: nSpin, nKpts, nOrbs, nIndepHam

    ! maximum allowed number of electrons in a single particle state
    real(dp) :: maxFill

    integer, allocatable :: nFilled(:,:), nEmpty(:,:)

    real(dp), allocatable :: dEi(:,:,:), dEiTmp(:,:,:,:,:)

    ! matrices for derivatives of terms in hamiltonian and outputs
    real(dp), allocatable :: dHam(:,:), idHam(:,:), dRho(:,:), idRho(:,:)
    real(dp) :: dqIn(orb%mOrb,nAtom,size(ham, dim=2))

    real(dp), allocatable :: dqBlockIn(:,:,:,:), SSqrReal(:,:)
    real(dp), allocatable :: dqBlockOut(:,:,:,:)

    ! derivative of potentials
    type(TPotentials) :: dPotential

    logical :: tSccCalc
    logical, allocatable :: tMetallic(:,:)

    real(dp), allocatable :: dPsiReal(:,:,:)
    complex(dp), allocatable :: dPsiCmplx(:,:,:,:)

    real(dp), allocatable :: dqNetAtom(:), dqNetAtomTmp(:,:,:)

    ! used for hybrid functional contributions, note this stays in the up/down representation
    ! throughout if spin polarised
    real(dp), pointer :: dRhoOutSqr(:,:,:), dRhoInSqr(:,:,:)
    real(dp), allocatable, target :: dRhoOut(:), dRhoIn(:)

    !> For transformation in the  case of degeneracies
    type(TRotateDegen), allocatable :: degenTransform(:)

    real(dp), allocatable :: dEf(:), dqOut(:,:,:), dqOutTmp(:,:,:,:,:)

    integer :: nIter, iOmega
    character(mc) :: atLabel

    logical :: isSccRequired

    type(TFileDescr) :: fd

    if (isRespKernelRPA) then
      nIter = 1
      isSccRequired = .false.
    else
      nIter = maxSccIter
      isSccRequired = isSccConvRequired
    end if

    call init_perturbation(parallelKS, this%tolDegen, nOrbs, nKpts, nSpin, nIndepHam, maxFill,&
        & filling, ham, nFilled, nEmpty, dHam, dRho, idHam, idRho, degenTransform, hybridXc,&
        & sSqrReal, over, neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell,&
        & dRhoOut, dRhoIn, dRhoInSqr, dRhoOutSqr, dPotential, orb, nAtom, tMetallic, neFermi,&
        & eigvals, tempElec, Ef, kWeight)

    write(stdOut,*)
    write(stdOut,*)'Perturbation calculation of atomic polarisability kernel'
    write(stdOut,*)

    allocate(dqOut(orb%mOrb, nAtom, nSpin))
    allocate(dqNetAtom(nAtom))
    allocate(dEi(nOrbs, nKpts, nSpin))
    if (isAutotestWritten.or.isTagResultsWritten) then
      allocate(dqOutTmp(orb%mOrb, nAtom, nSpin, nAtom, size(omega)))
      allocate(dqNetAtomTmp(nAtom, nAtom, size(omega)))
      dqOutTmp(:,:,:,:,:) = 0.0_dp
      dqNetAtomTmp(:,:,:) = 0.0_dp
    end if
    if (isTagResultsWritten) then
      allocate(dEiTmp(nOrbs, nKpts, nSpin, nAtom, size(omega)))
      dEiTmp(:,:,:,:,:) = 0.0_dp
    end if
    if (any(tMetallic)) then
      allocate(dEf(nIndepHam))
    end if

    tSccCalc = allocated(sccCalc)

    if (allocated(dftbU) .and. allocated(onsMEs)) then
      @:RAISE_ERROR(errStatus, -1, "Onsite corrected and DFTB+U terms currently not compatible for&
          & perturbation")
    end if
    if (allocated(dftbU) .or. allocated(onsMEs)) then
      allocate(dqBlockIn(orb%mOrb,orb%mOrb,nAtom,nSpin))
      allocate(dqBlockOut(orb%mOrb,orb%mOrb,nAtom,nSpin))
    end if

    ! if derivatives of valence wavefunctions needed. Note these will have an arbitrary set of
    ! global phases
    ! if (allocated(eigVecsReal)) then
    !   allocate(dPsiReal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2), nIndepHam, 3))
    ! else
    !   allocate(dPsiCmplx(size(eigvecsCplx,dim=1), size(eigvecsCplx,dim=2), nKpts, nIndepHam, 3))
    ! end if

    lpAtom: do iAt = 1, nAtom

      write(stdOut,*)'Derivative with respect to potential at atom ', iAt

      dqOut(:,:,:) = 0.0_dp
      dqIn(:,:,:) = 0.0_dp
      if (allocated(dftbU) .or. allocated(onsMEs)) then
        dqBlockIn(:,:,:,:) = 0.0_dp
        dqBlockOut(:,:,:,:) = 0.0_dp
      end if
      dPotential%extAtom(:,:) = 0.0_dp
      dPotential%extAtom(iAt,1) = 1.0_dp

      dPotential%extShell(:,:,:) = 0.0_dp
      dPotential%extBlock(:,:,:,:) = 0.0_dp
      call totalShift(dPotential%extShell, dPotential%extAtom, orb, species)
      call totalShift(dPotential%extBlock, dPotential%extShell, orb, species)

      lpOmega: do iOmega = 1, size(omega)

        dEi(:,:,:) = 0.0_dp
        call response(env, parallelKS, dPotential, nAtom, orb, species, neighbourList,&
            & nNeighbourSK, img2CentCell, iSparseStart, denseDesc, over, iEqOrbitals, sccCalc,&
            & sccTol, isSccRequired, nIter, chrgMixerReal, nMixElements, nIneqMixElements, dqIn,&
            & dqOut, hybridXc, nNeighbourCam, sSqrReal, dRhoInSqr, dRhoOutSqr, dRhoIn, dRhoOut,&
            & nSpin, maxFill, spinW, thirdOrd, dftbU, iEqBlockDftbu, onsMEs, iEqBlockOnSite,&
            & dqBlockIn, dqBlockOut, eigVals, degenTransform, dEi, dEf, filling, Ef,&
            & this%isEfFixed, dHam, idHam, dRho, idRho, tempElec, tMetallic, neFermi, nFilled,&
            & nEmpty, kPoint, kWeight, cellVec, iCellVec, eigVecsReal, eigVecsCplx, dPsiReal,&
            & dPsiCmplx, coord, errStatus, omega(iOmega), isHelical=isHelical, eta=this%eta)
          @:PROPAGATE_ERROR(errStatus)

      #:if WITH_SCALAPACK
        ! Add up and distribute eigenvalue derivatives from each processor
        call mpifx_allreduceip(env%mpi%globalComm, dEi, MPI_SUM)
      #:endif

        if (isTagResultsWritten) then
          dEiTmp(:,:,:,iAt,iOmega) = dEi
        end if

        write(stdOut,*)'Frontier orbital derivatives'
        do iS = 1, nIndepHam
          do iK = 1, nKpts
            write(stdOut,*)dEi(nFilled(iS, iK), iK, iS), dEi(nEmpty(iS, iK), iK, iS)
          end do
        end do

        call getOnsitePopulation(dRho(:,1), orb, iSparseStart, dqNetAtom)
        write(stdOut,*)'Derivatives of Mulliken and on-site (net) populations'
        do jAt = 1, nAtom
          write(stdOut,*)jAt, sum(dqOut(:,jAt,1)), dqNetAtom(jAt)
        end do

        if (isAutotestWritten.or.isTagResultsWritten) then
          dqOutTmp(:,:,:,iAt,iOmega) = dqOut
          dqNetAtomTmp(:,iAt,iOmega) = dqNetAtom
        end if

        if (isBandWritten) then
          write(atLabel,"(A,I0)")'ATOM ',iAt
          if (iAt == 1) then
            call writeDerivBandOut(derivVBandOut, dEi, kWeight, preLabel=atLabel)
          else
            call writeDerivBandOut(derivVBandOut, dEi, kWeight, isFileAppended=.true.,&
                & preLabel=atLabel)
          end if
        end if

        if (fdDetailedOut%isConnected()) then
          if (abs(omega(iOmega)) > epsilon(0.0_dp)) then
            write(fdDetailedOut%unit, format2U)"Response at omega = ", omega(iOmega), ' H ',&
                & omega(iOmega) * Hartree__eV, ' eV'
          else
            write(fdDetailedOut%unit, "(A)")"Static response:"
          end if
          write(fdDetailedOut%unit, "(A,I0)")'Derivatives wrt. a potential at atom ', iAt
          write(fdDetailedOut%unit, "(1X,A)")'Frontier orbital energy derivatives (a.u.)'
          write(fdDetailedOut%unit, "(1X,A,T14,A,T28,A)")"Spin Kpt","Last filled","First empty"
          do iS = 1, nIndepHam
            do iK = 1, nKpts
              write(fdDetailedOut%unit, "(1X,I2,I4,2F14.6)")iS, iK, dEi(nFilled(iS, iK), iK, iS),&
                  & dEi(nEmpty(iS, iK), iK, iS)
            end do
          end do
          write(fdDetailedOut%unit, *) 'Atomic population derivatives (a.u.)'
          write(fdDetailedOut%unit, "(1X,A,T10,A,T22,A)")"Atom","Mulliken","On-site"
          do jAt = 1, nAtom
            write(fdDetailedOut%unit, "(I5, 2F12.6)")jAt, sum(dqOut(:,jAt,1)), dqNetAtom(jAt)
          end do
          write(fdDetailedOut%unit, *)
        end if

      end do lpOmega

    end do lpAtom

    if (isAutotestWritten) then
      call openFile(fd, autotestTagFile, mode="a")
      call taggedWriter%write(fd%unit, tagLabels%dqdV, dqOutTmp)
      call taggedWriter%write(fd%unit, tagLabels%dqnetdV, dqNetAtomTmp)
      call closeFile(fd)
    end if
    if (isTagResultsWritten) then
      call openFile(fd, resultsTagFile, mode="a")
      call taggedWriter%write(fd%unit, tagLabels%dEigenDV, dEiTmp)
      call taggedWriter%write(fd%unit, tagLabels%dqdV, dqOutTmp)
      call taggedWriter%write(fd%unit, tagLabels%dqnetdV, dqNetAtomTmp)
      call closeFile(fd)
    end if

  end subroutine wrtVAtom


  !> Response with respect to external charge positions (at q=0)
  subroutine dxExtCharges(this, env, parallelKS, isAutotestWritten, autotestTagFile,&
      & isTagResultsWritten, resultsTagFile, taggedWriter, isBandWritten, tWriteDetailedOut,&
      & fdDetailedOut, filling, eigvals, eigVecsReal, eigVecsCplx, ham, over, orb, nAtom, species,&
      & neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell, sccCalc, maxSccIter,&
      & sccTol, nMixElements, nIneqMixElements, iEqOrbitals, tempElec, Ef, spinW, thirdOrd, dftbU,&
      & iEqBlockDftbu, onsMEs, iEqBlockOnSite, hybridXc, nNeighbourCam, chrgMixerReal, isPeriodic,&
      & coord, kPoint, kWeight, iCellVec, cellVec, nEFermi, wrtWhichCharges, nCombinedCharges,&
      & wrtCombinedCharges, combinedJacobian, dqdxExt, errStatus, omega, isHelical)

    !> Instance
    class(TResponse), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> The k-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Should regression test data be written
    logical, intent(in) :: isAutotestWritten

    !> Name of output file
    character(*), intent(in) :: autotestTagFile

    !> Produce machine readable results
    logical, intent(in) :: isTagResultsWritten

    !> Name of output file
    character(*), intent(in) :: resultsTagFile

    !> Tagged writer object
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Should eigenvalue (band) data derivatives be written to disc
    logical, intent(in) :: isBandWritten

    !> Should detailed.out be written to
    logical, intent(in) :: tWriteDetailedOut

    !> File descriptor for the human readable output
    type(TFileDescr), intent(in) :: fdDetailedOut

    !> Filling
    real(dp), intent(in) :: filling(:,:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in), allocatable :: eigVecsReal(:,:,:)

    !> Ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> Sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> Chemical species
    integer, intent(in) :: species(:)

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> SCC module internal variables
    type(TScc), intent(inout), allocatable :: sccCalc

    !> Maximal number of SCC iterations
    integer, intent(in) :: maxSccIter

    !> Tolerance for SCC convergence
    real(dp), intent(in) :: sccTol

    !> Nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
    !> blocks (if a DFTB+U or onsite correction calculation)
    integer, intent(in) :: nMixElements

    !> Nr. of inequivalent charges
    integer, intent(in) :: nIneqMixElements

    !> Equivalence relations between orbitals
    integer, intent(in), allocatable :: iEqOrbitals(:,:,:)

    !> Onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> Third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> Equivalence mapping for dual charge blocks
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> Data for hybrid functional calculation
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Number of neighbours for each of the atoms for the exchange contributions of CAM functionals
    integer, intent(in), allocatable :: nNeighbourCam(:)

    !> Charge mixing object
    class(TMixerReal), intent(inout), allocatable :: chrgMixerReal

    !> Is this a periodic geometry
    logical, intent(in) :: isPeriodic

    !> Coordinates of all atoms including images
    real(dp), intent(in) :: coord(:,:)

    !> The k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Number of electrons at the Fermi energy (if metallic)
    real(dp), allocatable, intent(inout) :: neFermi(:)

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    !> Driving frequencies (including potentially 0 for static)
    real(dp), intent(in) :: omega(:)

    !> List of which specific external charge positions (MM atoms) the DFTB charges are to be
    !! differentiated (note, incompatible with wrtCombinedCharges)
    integer, intent(in), allocatable :: wrtWhichCharges(:)

    !> Number of combined charges in each group
    integer, intent(in), allocatable :: nCombinedCharges(:)

    !> List of groups of external charge positions (MM atoms) the DFTB charges are to be
    !! differentiated (note, incompatible with wrtWhichCharges)
    integer, intent(in), allocatable :: wrtCombinedCharges(:,:)

    !> Weights for external charge positions [nCharges, iCart, nDerivs]
    real(dp), intent(in), allocatable :: combinedJacobian(:,:,:)

    !> Output, charge derivatives
    real(dp), allocatable, intent(inout), target :: dqdxExt(:,:,:)

    !> Is the geometry helical
    logical, intent(in), optional :: isHelical

    integer :: iS, iK, iExtChrgWRT, nDerivs, iCart, iAt, jAt, jCart
    integer :: nSpin, nKpts, nOrbs, nIndepHam

    ! maximum allowed number of electrons in a single particle state
    real(dp) :: maxFill

    integer, allocatable :: nFilled(:,:), nEmpty(:,:)

    real(dp), allocatable :: dEi(:,:,:), dEiTmp(:,:,:,:,:)

    ! matrices for derivatives of terms in hamiltonian and outputs
    real(dp), allocatable :: dHam(:,:), idHam(:,:), dRho(:,:), idRho(:,:)
    real(dp) :: dqIn(orb%mOrb,nAtom,size(ham, dim=2))

    real(dp), allocatable :: dqBlockIn(:,:,:,:), SSqrReal(:,:)
    real(dp), allocatable :: dqBlockOut(:,:,:,:)

    ! derivative of potentials
    type(TPotentials) :: dPotential

    logical, allocatable :: tMetallic(:,:)

    real(dp), allocatable :: dPsiReal(:,:,:)
    complex(dp), allocatable :: dPsiCmplx(:,:,:,:)

    real(dp), allocatable :: dqNetAtom(:)

    ! used for hybrid functional contributions, note this stays in the up/down representation
    ! throughout if spin polarised
    real(dp), pointer :: dRhoOutSqr(:,:,:), dRhoInSqr(:,:,:)
    real(dp), allocatable, target :: dRhoOut(:), dRhoIn(:)

    !> For transformation in the  case of degeneracies
    type(TRotateDegen), allocatable :: degenTransform(:)

    real(dp), allocatable :: dEf(:), dqOut(:,:,:)

    integer :: nIter, iOmega
    character(mc) :: atLabel

    logical :: isSccRequired

    integer :: fdResults
    real(dp), allocatable, target :: dqdxExtWork(:,:,:)
    real(dp), pointer :: pdqdxExt(:,:,:) => null()

    integer :: nExtCharge, iDeriv, jCharge
    real(dp), allocatable :: extCoord(:,:), extCharge(:), dgammaQMMM(:,:), blurWidths(:)
    real(dp), allocatable :: dgammaTmp(:,:)

    @:ASSERT(.not.(allocated(wrtWhichCharges) .and. allocated(wrtCombinedCharges)))
    @:ASSERT(allocated(wrtCombinedCharges) .eqv. allocated(nCombinedCharges))
    @:ASSERT(allocated(wrtCombinedCharges) .eqv. allocated(combinedJacobian))

    ! obtain the external charge coordinates and their values
    call sccCalc%getExternalCharges(nExtCharge, extCoord, extCharge, blurWidths=blurWidths)
    if (nExtCharge < 1) then
      write (stdOut,*) "No external charges, nothing to do in dxExtCharges."
      return
    end if

    nDerivs = countDerivs(wrtWhichCharges, nCombinedCharges, wrtCombinedCharges, nExtCharge)

    if (nDerivs < 1) then
      write (stdOut,*) "No requested derivatives wrt to external charges, nothing to do in&
          & dxExtCharges."
      return
    end if

    if (isAutotestWritten .or. isTagResultsWritten .or. tWriteDetailedOut) then
      if (allocated(dqdxExt)) then
        @:ASSERT(all(shape(dqdxExt) >= [nAtom, 3, nDerivs]))
        pdqdxExt => dqdxExt
      else
        allocate(dqdxExtWork(nAtom, 3, nDerivs))
        pdqdxExt => dqdxExtWork
      end if
    else
      if (allocated(dqdxExt)) then
        @:ASSERT(all(shape(dqdxExt) >= [nAtom, 3, nDerivs]))
        pdqdxExt => dqdxExt
      end if
    end if

    if (.not. allocated(sccCalc)) then
      @:RAISE_ERROR(errStatus, -1, "SCC currently required for external charge derivatives")
    end if

    write(stdOut,*)
    write(stdOut,*)'Perturbation calculation of derivative with respect to external charges'
    write(stdOut,*)

    nIter = maxSccIter
    isSccRequired = .true.

    call init_perturbation(parallelKS, this%tolDegen, nOrbs, nKpts, nSpin, nIndepHam, maxFill,&
        & filling, ham, nFilled, nEmpty, dHam, dRho, idHam, idRho, degenTransform, hybridXc,&
        & sSqrReal, over, neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell,&
        & dRhoOut, dRhoIn, dRhoInSqr, dRhoOutSqr, dPotential, orb, nAtom, tMetallic, neFermi,&
        & eigvals, tempElec, Ef, kWeight)

    allocate(dqOut(orb%mOrb, nAtom, nSpin))
    allocate(dqNetAtom(nAtom))
    allocate(dEi(nOrbs, nKpts, nSpin))
    if (any(tMetallic)) then
      allocate(dEf(nIndepHam))
    end if

    if (allocated(dftbU) .and. allocated(onsMEs)) then
      @:RAISE_ERROR(errStatus, -1, "Onsite corrected and DFTB+U terms currently not compatible for&
          & perturbation")
    end if
    if (allocated(dftbU) .or. allocated(onsMEs)) then
      allocate(dqBlockIn(orb%mOrb,orb%mOrb,nAtom,nSpin))
      allocate(dqBlockOut(orb%mOrb,orb%mOrb,nAtom,nSpin))
    end if

    ! derivatives of QM--MM 1/r w.r.t. coordinates of MM atoms
    allocate(dgammaQMMM(3, nAtom))

    lpExtChrg: do iDeriv = 1, nDerivs

      if (allocated(wrtWhichCharges)) then
        iExtChrgWrt = wrtWhichCharges(iDeriv)
        write(stdOut,"(A,I0)")'Derivative with respect to external charge ', iExtChrgWrt
      else if (allocated(wrtCombinedCharges)) then
        write(stdOut,"(A,I0,A)")'Derivative with respect to external charge group ',iDeriv, ':'
        write(stdOut,*)wrtCombinedCharges(:nCombinedCharges(iDeriv), iDeriv)
      else
        iExtChrgWRT = iDeriv
        write(stdOut,"(A,I0)")'Derivative with respect to external charge ', iExtChrgWRT
      end if

      if (allocated(wrtCombinedCharges)) then
        dgammaQMMM(:,:) = 0.0_dp
        allocate(dgammaTmp, source=dgammaQMMM)
        do jCharge = 1, nCombinedCharges(iDeriv)
          if (isPeriodic) then
            call calcInvRPrimeAsymm(nAtom, coord, nExtCharge, extCoord, extCharge,&
                & sccCalc%coulomb%rLatPoints_, sccCalc%coulomb%gLatPoints_, sccCalc%coulomb%alpha,&
                & sccCalc%coulomb%volume_, wrtCombinedCharges(jCharge, iDeriv), dgammaTmp,&
                & blurWidths=blurWidths)
          else
            call calcInvRPrimeAsymm(nAtom, coord, nExtCharge, extCoord, extCharge,&
                & wrtCombinedCharges(jCharge, iDeriv), dgammaTmp, blurWidths=blurWidths)
          end if
          do jCart = 1, 3
            dgammaQMMM(jCart,:) = dgammaQMMM(jCart,:)&
                & + combinedJacobian(jCharge, jCart, iDeriv) * dgammaTmp(jCart,:)
          end do
        end do
        deallocate(dgammaTmp)
      else
        if (isPeriodic) then
          call calcInvRPrimeAsymm(nAtom, coord, nExtCharge, extCoord, extCharge,&
              & sccCalc%coulomb%rLatPoints_, sccCalc%coulomb%gLatPoints_, sccCalc%coulomb%alpha,&
              & sccCalc%coulomb%volume_, iExtChrgWRT, dgammaQMMM, blurWidths=blurWidths)
        else
          call calcInvRPrimeAsymm(nAtom, coord, nExtCharge, extCoord, extCharge, iExtChrgWRT,&
              & dgammaQMMM, blurWidths=blurWidths)
        end if
      end if

      dqOut(:,:,:) = 0.0_dp
      dqIn(:,:,:) = 0.0_dp
      if (allocated(dftbU) .or. allocated(onsMEs)) then
        dqBlockIn(:,:,:,:) = 0.0_dp
        dqBlockOut(:,:,:,:) = 0.0_dp
      end if

      ! perturbation direction
      lpCart: do iCart = 1, 3

        write(stdOut,"(A,A,A,I0)")'Calculating derivative for displacement along ', &
            & trim(direction(iCart)),' for charge ', iExtChrgWRT

        dPotential%extAtom(:,:) = 0.0_dp
        dPotential%extAtom(:,1) = dgammaQMMM(iCart, :)
        dPotential%extShell(:,:,:) = 0.0_dp
        dPotential%extBlock(:,:,:,:) = 0.0_dp
        call totalShift(dPotential%extShell, dPotential%extAtom, orb, species)
        call totalShift(dPotential%extBlock, dPotential%extShell, orb, species)

        lpOmega: do iOmega = 1, size(omega)

          dEi(:,:,:) = 0.0_dp
          call response(env, parallelKS, dPotential, nAtom, orb, species, neighbourList,&
              & nNeighbourSK, img2CentCell, iSparseStart, denseDesc, over, iEqOrbitals, sccCalc,&
              & sccTol, isSccRequired, nIter, chrgMixerReal, nMixElements, nIneqMixElements, dqIn,&
              & dqOut, hybridXc, nNeighbourCam, sSqrReal, dRhoInSqr, dRhoOutSqr, dRhoIn, dRhoOut,&
              & nSpin, maxFill, spinW, thirdOrd, dftbU, iEqBlockDftbu, onsMEs, iEqBlockOnSite,&
              & dqBlockIn, dqBlockOut, eigVals, degenTransform, dEi, dEf, filling, Ef,&
              & this%isEfFixed, dHam, idHam, dRho, idRho, tempElec, tMetallic, neFermi, nFilled,&
              & nEmpty, kPoint, kWeight, cellVec, iCellVec, eigVecsReal, eigVecsCplx, dPsiReal,&
              & dPsiCmplx, coord, errStatus, omega(iOmega), isHelical=isHelical, eta=this%eta)
          @:PROPAGATE_ERROR(errStatus)

      #:if WITH_SCALAPACK
          ! Add up and distribute eigenvalue derivatives from each processor
          call mpifx_allreduceip(env%mpi%globalComm, dEi, MPI_SUM)
      #:endif

          call getOnsitePopulation(dRho(:,1), orb, iSparseStart, dqNetAtom)
          write(stdOut,*)'Derivatives of Mulliken and on-site (net) populations'
          do jAt = 1, nAtom
            write(stdOut,*)jAt, -sum(dqOut(:,jAt,1)), -dqNetAtom(jAt)
          end do

          if (associated(pdqdxExt)) then
            pdqdxExt(:, iCart, iDeriv) = -sum(dqOut(:, :, 1), dim=1)
          end if

        end do lpOmega

      end do lpCart

    end do lpExtChrg

    if (isAutotestWritten) then
      open(newunit=fdResults, file=autoTestTagFile, position="append")
      call taggedWriter%write(fdResults, tagLabels%dqdxExt, -pdqdxExt)
      close(fdResults)
    end if
    if (isTagResultsWritten) then
      open(newunit=fdResults, file=resultsTagFile, position="append")
      call taggedWriter%write(fdResults, tagLabels%dqdxExt, -pdqdxExt)
      close(fdResults)
    end if

    if (tWriteDetailedOut) then
      @:ASSERT(fdDetailedOut%isConnected())
      write(fdDetailedOut%unit,"(A)")'Derivatives of atomic charges w.r.t. coordinates of external&
          & charges'
      do iDeriv = 1, nDerivs
        if (allocated(wrtWhichCharges)) then
          iExtChrgWrt = wrtWhichCharges(iDeriv)
        else if (allocated(wrtCombinedCharges)) then
        else
          iExtChrgWRT = iDeriv
        end if
        write(fdDetailedOut%unit,"(A,I0,T10,A,T26,A,T42,A)")"Chg#", iExtChrgWrt, "x", "y", "z"
        do jAt = 1, nAtom
          write(fdDetailedOut%unit, "(A,I0,T4,3F16.8)")'At',jAt, -pdqdxExt(jAt, :, iDeriv)&
              & * AA__Bohr
        end do
        write(fdDetailedOut%unit, *)
      end do
    end if

  end subroutine dxExtCharges


  !> Static (frequency independent) perturbation at q=0 with respect to atomic postitions
  subroutine dxAtom(this, env, parallelKS, filling, eigvals, eigVecsReal, eigVecsCplx, rhoPrim,&
      & potential, qOrb, q0, ham, over, skHamCont, skOverCont, nonSccDeriv, orb, nAtom, species,&
      & speciesnames, neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell, coord,&
      & sccCalc, maxSccIter, sccTol, nMixElements, nIneqMixElements, iEqOrbitals, tempElec, Ef,&
      & tFixEf, spinW, thirdOrd, DftbU, iEqBlockDftbu, onsMEs, iEqBlockOnSite, hybridXc,&
      & nNeighbourLC, chrgMixerReal, isBandWritten, taggedWriter, isAutotestWritten,&
      & autoTestTagFile, isTagResultsWritten, taggedResultsFile, tWriteDetailedOut, fdDetailedOut,&
      & kPoint, kWeight, iCellVec, cellVec, isPeriodic, tPrintMulliken, wrtWhichCharges,&
      & nCombinedCharges, wrtCombinedCharges, combinedJacobian, errStatus, isHelical, dqdx)

    !> Instance
    class(TResponse), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Fillings of unperturbed system
    real(dp), intent(in) :: filling(:,:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in), allocatable :: eigVecsReal(:,:,:)

    !> Ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Unperturbed density matrix in sparse format
    real(dp), intent(in) :: rhoPrim(:,:)

    !> Unperturbed potentials
    type(TPotentials), intent(in) :: potential

    !> Electrons in each atomic orbital
    real(dp), intent(in) :: qOrb(:,:,:)

    !> Reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> Sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Container for SK Hamiltonian integrals
    type(TSlakoCont), intent(in) :: skHamCont

    !> Container for SK overlap integrals
    type(TSlakoCont), intent(in) :: skOverCont

    !> Method for calculating derivatives of S and H0 with respect to atom positions
    type(TNonSccDiff), intent(in) :: nonSccDeriv

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> Chemical species
    integer, intent(in) :: species(:)

    !> Label for each atomic chemical species
    character(mc), intent(in) :: speciesnames(:)

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> SCC module internal variables
    type(TScc), intent(inout), allocatable :: sccCalc

    !> Maximal number of SCC iterations
    integer, intent(in) :: maxSccIter

    !> Tolerance for SCC convergence
    real(dp), intent(in) :: sccTol

    !> Nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
    !! blocks (if tDFTBU or onsite corrections)
    integer, intent(in) :: nMixElements

    !> Nr. of inequivalent charges
    integer, intent(in) :: nIneqMixElements

    !> Equivalence relations between orbitals
    integer, allocatable, intent(in) :: iEqOrbitals(:,:,:)

    !> Onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> Spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> Third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> Equivalence mapping for dual charge blocks
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> Data for hybrid functionals
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !! functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Charge mixing object
    class(TMixerReal), allocatable, intent(inout) :: chrgMixerReal

    !> Should eigenvalue (band) data derivatives be written to disc
    logical, intent(in) :: isBandWritten

    !> Tagged writer object
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Should regression test data be written
    logical, intent(in) :: isAutotestWritten

    !> File name for regression data
    character(*), intent(in) :: autoTestTagFile

    !> Should machine readable output data be written
    logical, intent(in) :: isTagResultsWritten

    !> File name for machine readable results data
    character(*), intent(in) :: taggedResultsFile

    !> Should detailed.out be written to
    logical, intent(in) :: tWriteDetailedOut

    !> File id for detailed.out
    integer, intent(in) :: fdDetailedOut

    !> K-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Is this a periodic geometry
    logical, intent(in) :: isPeriodic

    !> Should Mulliken populations be generated/output
    logical, intent(in) :: tPrintMulliken

    !> List of which specific atom positions the DFTB charges are to be differentiated (note,
    !! incompatible with wrtCombinedCharges)
    integer, intent(in), allocatable :: wrtWhichCharges(:)

    !> Number of combined charges in each group
    integer, intent(in), allocatable :: nCombinedCharges(:)

    !> List of groups of positions the DFTB charges are to be differentiated (note, incompatible
    !! with wrtWhichCharges)
    integer, intent(in), allocatable :: wrtCombinedCharges(:,:)

    !> Weights for external charge positions [nCharges, iCart, nDerivs]
    real(dp), intent(in), allocatable :: combinedJacobian(:,:,:)

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    !> Is the geometry helical
    logical, intent(in), optional :: isHelical

    real(dp), optional, intent(inout) :: dqdx(:,:,:)

    integer :: iS, iK, iKS, iNeigh, iCart, iSCC, iLev, iSh, iSp, jAt, jAtf, iOrb, jCart
    integer :: iAt, nSpin, nKpts, nOrbs, nIndepHam, nDerivs, iDeriv, jCharge

    ! maximum allowed number of electrons in a single particle state
    real(dp) :: maxFill

    integer, allocatable :: nFilled(:,:), nEmpty(:,:)

    integer :: ii, jj, iGlob, jGlob
    integer :: iSccIter
    logical :: tStopSCC

    ! matrices for derivatives of terms in hamiltonian and outputs
    real(dp), allocatable :: dHam(:,:), idHam(:,:), dOver(:,:), dH0(:,:)
    real(dp), allocatable :: dOverTmp(:,:), dH0Tmp(:,:), dOverWork(:)

    ! overlap derivative terms in potential, omega dS + d(delta q gammma) S
    real(dp), allocatable :: sOmega(:,:,:)

    real(dp), allocatable :: drho(:,:)
    real(dp) :: drhoExtra(size(over),size(ham, dim=2))
    real(dp), allocatable :: idRho(:,:), idRhoExtra(:,:)
    real(dp) :: dqIn(orb%mOrb,nAtom,size(ham, dim=2))
    real(dp), allocatable :: dqOut(:,:,:,:,:) !(orb%mOrb, nAtom, size(ham, dim=2), 3, nAtom)
    real(dp) :: dqUpDown(orb%mOrb, size(ham, dim=2))
    real(dp) :: dqInpRed(nMixElements), dqOutRed(nMixElements)
    real(dp) :: dqDiffRed(nMixElements), sccErrorQ
    real(dp) :: dqPerShell(orb%mShell,nAtom,size(ham, dim=2))

    ! eigenvalue weighted vectors
    real(dp), allocatable :: eCiReal(:, :, :)
    complex(dp), allocatable :: eCiCplx(:, :, :)

    real(dp), allocatable :: vAt(:,:), vdgamma(:,:,:)
    real(dp), allocatable :: dGamma3(:,:), qAtom(:)

    real(dp), allocatable :: dqBlockIn(:,:,:,:), SSqrReal(:,:)
    real(dp), allocatable :: dqBlockOut(:,:,:,:)
    real(dp), allocatable :: dummy(:,:,:,:)

    ! derivative of potentials
    type(TPotentials) :: dPotential

    real(dp), allocatable :: shellPot(:,:,:), atomPot(:,:)

    logical :: tSccCalc, tConverged
    logical, allocatable :: tMetallic(:,:)

    real(dp), allocatable :: dEi(:,:,:,:,:), dEiTmp(:,:,:)
    real(dp), allocatable :: dPsiReal(:,:,:)
    complex(dp), allocatable :: dPsiCmplx(:,:,:,:)

    integer :: fdResults

    ! used for hybrid functional contributions, note this stays in the up/down representation
    ! throughout if spin polarised
    real(dp), pointer :: dRhoOutSqr(:,:,:), dRhoInSqr(:,:,:)
    real(dp), allocatable, target :: dRhoOut(:), dRhoIn(:)

    ! non-variational part of charge derivative
    real(dp), allocatable :: dqNonVariational(:,:,:), dqNonVariationalBlock(:,:,:,:)

    real(dp) :: dDipole(3)

    ! derivative of shift due to external charges w.r.t. coordinates of an atom
    real(dp), allocatable :: extShiftDerivative(:,:,:)

    !> For transformation in the  case of degeneracies
    type(TRotateDegen), allocatable :: degenTransform(:)

    !> Driving frequencies (including potentially 0 for static)
    real(dp) :: omega = 0.0_dp

    ! Number of electrons at the Fermi energy (if metallic)
    real(dp), allocatable :: neFermi(:)

    real(dp), allocatable :: bornCharges(:,:,:)

    character(lc) :: tmpStr
    logical :: areChargeDerivsNeeded

  #:if WITH_SCALAPACK
    ! need distributed matrix descriptors
    integer :: desc(DLEN_), nn

    type(blocklist) :: blocks
    integer :: blockSize, iLoc

    nn = denseDesc%fullSize
    call scalafx_getdescriptor(env%blacs%orbitalGrid, nn, nn, env%blacs%rowBlockSize,&
        & env%blacs%columnBlockSize, desc)
    call blocks%init(env%blacs%orbitalGrid, desc, "c")

  #:endif

    @:ASSERT(.not.(allocated(wrtWhichCharges) .and. allocated(wrtCombinedCharges)))
    @:ASSERT(allocated(wrtCombinedCharges) .eqv. allocated(nCombinedCharges))
    @:ASSERT(allocated(wrtCombinedCharges) .eqv. allocated(combinedJacobian))

    nDerivs =  countDerivs(wrtWhichCharges, nCombinedCharges, wrtCombinedCharges, nAtom)

    if (nDerivs < 1) then
      write (stdOut,*) "No requested derivatives wrt to external charges, nothing to do in&
          & dxAtom."
      return
    end if

    areChargeDerivsNeeded = present(dqdx)
    if (areChargeDerivsNeeded) then
      @:ASSERT(all(shape(dqdx) >= [nAtom, 3, nDerivs]))
    end if

    call init_perturbation(parallelKS, this%tolDegen, nOrbs, nKpts, nSpin, nIndepHam, maxFill,&
        & filling, ham, nFilled, nEmpty, dHam, dRho, idHam, idRho, degenTransform, hybridXc,&
        & sSqrReal, over, neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell,&
        & dRhoOut, dRhoIn, dRhoInSqr, dRhoOutSqr, dPotential, orb, nAtom, tMetallic, neFermi,&
        & eigvals, tempElec, Ef, kWeight)

    if (any(tMetallic)) then
      @:RAISE_ERROR(errStatus, -1, "Atom derivative perturbations are not currently supported for&
          & metallic systems")
    end if

    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iS = parallelKS%localKS(2, iKS)
      if (degenTransform(iKS)%isAnyDegenerate(eigvals(:, iK, iS))) then
        @:RAISE_ERROR(errStatus, -1, "Atom derivative perturbations are not currently supported for&
            & high-symmetry/degenerate systems")
      end if
    end do

    allocate(dEiTmp(nOrbs, 1, nSpin))

    allocate(dqNonVariational(orb%mOrb,nAtom,nSpin))

  #:if WITH_SCALAPACK

    if (allocated(eigVecsReal)) then

      allocate(eCiReal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2), size(eigVecsReal,dim=3)))

      ! e_i |c_i>
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iS = parallelKS%localKS(2, iKS)
        do ii = 1, size(blocks)
          call blocks%getblock(ii, iGlob, iLoc, blockSize)
          iK = parallelKS%localKS(1, iKS)
          iS = parallelKS%localKS(2, iKS)
          do jj = 0, min(blockSize - 1, nOrbs - iGlob)
            eCiReal(:, iLoc + jj, iKS) = eigVals(iGlob + jj, iK, iS)&
                & * eigVecsReal(:, iLoc + jj, iKS)
          end do
        end do
      end do

    else

      allocate(eCiCplx(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2), size(eigVecsCplx,dim=3)))

      ! e_i |c_i>
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iS = parallelKS%localKS(2, iKS)
        do ii = 1, size(blocks)
          call blocks%getblock(ii, iGlob, iLoc, blockSize)
          do jj = 0, min(blockSize - 1, nOrbs - iGlob)
            eCiCplx(:, iLoc + jj, iKS) = eigVals(iGlob + jj, iK, iS)&
                & * eigVecsCplx(:, iLoc + jj, iKS)
          end do
        end do
      end do

    end if

  #:else

    if (allocated(eigVecsReal)) then
      allocate(eCiReal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2), size(eigVecsReal,dim=3)))
      do iKS = 1, size(eigVecsReal,dim=3)
        do iOrb = 1, size(eigVecsReal,dim=2)
          eCiReal(:,iOrb, iKS) = eigVecsReal(:,iOrb, iKS) * eigVals(iOrb, 1, iKS)
        end do
      end do
    else
      allocate(eCiCplx(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2), size(eigVecsCplx,dim=3)))
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iS = parallelKS%localKS(2, iKS)
        do iOrb = 1, size(eigVecsCplx,dim=2)
          eCiCplx(:,iOrb, iKS) = eigVecsCplx(:,iOrb, iKS) * eigVals(iOrb, iK, iS)
        end do
      end do
    end if

  #:endif

    allocate(dOver(size(ham,dim=1),3))
    allocate(dH0(size(ham,dim=1),3))

    tSccCalc = allocated(sccCalc)

    if (tSccCalc) then
      if (allocated(thirdOrd)) then
        allocate(sOmega(size(ham,dim=1),nSpin,3))
        allocate(dGamma3(nAtom, nAtom))
        allocate(qAtom(nAtom))
        allocate(extShiftDerivative(3, nAtom, nSpin))
      else
        allocate(sOmega(size(ham,dim=1),nSpin,2))
      end if
      ! These are both spin-free
      allocate(vAt(nAtom, 1))
      allocate(vdgamma(orb%mShell, nAtom, 1))
    end if

    if (allocated(hybridXc)) then
    #:if WITH_SCALAPACK
      @:RAISE_ERROR(errStatus, -1, "Range separation not supported for MPI at the moment")
    #:endif
      allocate(SSqrReal(nOrbs, nOrbs))
      SSqrReal(:,:) = 0.0_dp
      call unpackHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell)
      allocate(dRhoOut(nOrbs * nOrbs * nSpin))
      dRhoOutSqr(1:nOrbs, 1:nOrbs, 1:nSpin) => dRhoOut(:nOrbs*nOrbs*nSpin)
      allocate(dRhoIn(nOrbs * nOrbs * nSpin))
      dRhoInSqr(1:nOrbs, 1:nOrbs, 1:nSpin) => dRhoIn(:nOrbs*nOrbs*nSpin)
    else
      dRhoInSqr => null()
      dRhoOutSqr => null()
    end if

    !if (tDFTBU .or. allocated(onsMEs)) then
    !  allocate(dqBlockIn(orb%mOrb,orb%mOrb,nAtom,nSpin))
    !  allocate(dqBlockOut(orb%mOrb,orb%mOrb,nAtom,nSpin))
    !  allocate(dqNonVariationalBlock(orb%mOrb,orb%mOrb,nAtom,nSpin))
    !end if

    if (allocated(spinW) .or. allocated(thirdOrd)) then
      allocate(shellPot(orb%mShell, nAtom, nSpin))
    end if
    if (allocated(thirdOrd)) then
      allocate(atomPot(nAtom, nSpin))
    end if

    if (allocated(wrtWhichCharges)) then
      allocate(dqOut(orb%mOrb, nAtom, size(ham, dim=2), 3, size(wrtWhichCharges)), source=0.0_dp)
      allocate(dEi(nOrbs, 1, nSpin, 3, size(wrtWhichCharges)), source=0.0_dp)
    else
      allocate(dqOut(orb%mOrb, nAtom, size(ham, dim=2), 3, nAtom), source=0.0_dp)
      allocate(dEi(nOrbs, 1, nSpin, 3, nAtom), source=0.0_dp)
    end if

    ! Displacement to differentiate wrt
    lpDerivs: do iDeriv = 1, nDerivs

      if (allocated(wrtWhichCharges)) then
        iAt = wrtWhichCharges(iDeriv)
        write(stdOut,"(A,I0)")'Derivative with respect to atom ', iAt
      else if (allocated(wrtCombinedCharges)) then
        write(stdOut,"(A,I0,A)")'Derivative with respect to atom group ',iDeriv, ':'
        write(stdOut,*)wrtCombinedCharges(:nCombinedCharges(iDeriv), iDeriv)
      else
        iAt = iDeriv
        write(stdOut,"(A,I0)")'Derivative with respect to atom ', iAt
      end if

      if (allocated(wrtCombinedCharges)) then
        dOver(:,:) = 0.0_dp
        dH0(:,:) = 0.0_dp
        allocate(dOverTmp, source=dOver)
        allocate(dH0Tmp, source=dH0)
        do jCharge = 1, nCombinedCharges(iDeriv)
          call nonSccDeriv%getFirstDeriv(dOverTmp, env, skOverCont, coord, species,&
              & wrtCombinedCharges(jCharge, iDeriv), orb, nNeighbourSK, neighbourList%iNeighbour,&
              & iSparseStart)
          call nonSccDeriv%getFirstDeriv(dH0Tmp, env, skHamCont, coord, species,&
              & wrtCombinedCharges(jCharge, iDeriv), orb, nNeighbourSK, neighbourList%iNeighbour,&
              & iSparseStart)
          do jCart = 1, 3
            dOver(:,jCart) = dOver(:,jCart) + combinedJacobian(jCharge, jCart, iDeriv)&
                & * dOver(:,jCart)
            dH0(:,jCart) = dH0(:,jCart) + combinedJacobian(jCharge, jCart, iDeriv) * dH0(:,jCart)
          end do
        end do
        deallocate(dOverTmp)
        deallocate(dH0Tmp)
      else
        call nonSccDeriv%getFirstDeriv(dOver, env, skOverCont, coord, species, iAt, orb,&
            & nNeighbourSK, neighbourList%iNeighbour, iSparseStart)
        call nonSccDeriv%getFirstDeriv(dH0, env, skHamCont, coord, species, iAt, orb, nNeighbourSK,&
            & neighbourList%iNeighbour, iSparseStart)
      end if

      ! perturbation direction
      lpCart: do iCart = 1, 3

        write(stdOut,"(A,A,A,I0)")'Calculating derivative for displacement along ', &
            & trim(direction(iCart)),' for atom ', iAt

        if (tSccCalc) then

          sOmega(:,:,:) = 0.0_dp
          ! First part, omega dS
          call addShift(env, sOmega(:,:,1), dOver(:,iCart), nNeighbourSK, neighbourList%iNeighbour,&
              & species, orb, iSparseStart, nAtom, img2CentCell, potential%intBlock,&
              & isInputZero=.true.)

          vdgamma(:,:,:) = 0.0_dp
          vAt(:,:) = 0.0_dp
          call sccCalc%updateCoords(env, coord, coord, species, neighbourList)
          call sccCalc%updateCharges(env, qOrb, orb, species, q0)
          call sccCalc%addPotentialDeriv(env, vAt(:,1), vdgamma, species, neighbourList%iNeighbour,&
              & img2CentCell, coord, orb, iCart, iAt)

          if (allocated(thirdOrd)) then

            ! derivatives of 3rd order Gamma matrices w.r.t. atom coordinates
            call thirdOrd%getGamma3Deriv(species, coord, iAt, iCart, dGamma3)

            ! Non-variational contribution
            ! 1/3 sum_C ( dGamma_CA/dx Dq_C^2 + 2 dGamma_AC/dx Dq_A Dq_C )
            qAtom(:) = sum(qOrb(:,:,1) - q0(:,:,1), dim=1)
            do jAt = 1, nAtom
              vAt(:, 1) = vAt(:, 1) + (dGamma3(jAt, :) * qAtom(jAt)**2.0_dp &
                  & + 2.0_dp * dGamma3(:, jAt) * qAtom(jAt) * qAtom(:))
            end do

          end if

          call totalShift(vdgamma, vAt, orb, species)

        end if

        ! non-variational part of the charge change due to basis derivatives
        if (tPrintMulliken) then
          dqNonVariational(:,:,:) = 0.0_dp

          do iS = 1, nSpin
            call mulliken(env, dqNonVariational(:,:,iS), dOver(:,iCart), rhoPrim(:,iS), orb,&
                & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
          end do
          !if (tDFTBU .or. allocated(onsMEs)) then
          !  dqNonVariationalBlock(:,:,:,:) = 0.0_dp
          !  do iS = 1, nSpin
          !    call mulliken(env, dqNonVariationalBlock(:,:,:,iS), dOver(:,iCart), rhoPrim(:,iS),&
          !        & orb, neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
          !  end do
          !end if
        end if

        dqIn(:,:,:) = 0.0_dp
        !if (tDFTBU .or. allocated(onsMEs)) then
        !  dqBlockIn(:,:,:,:) = 0.0_dp
        !  dqBlockOut(:,:,:,:) = 0.0_dp
        !end if

        if (tSccCalc) then
          call chrgMixerReal%reset(nMixElements)
          dqInpRed(:) = 0.0_dp
          dqPerShell(:,:,:) = 0.0_dp
          if (allocated(hybridXc)) then
            dRhoIn(:) = 0.0_dp
            dRhoOut(:) = 0.0_dp
          end if

          write(stdOut,"(1X,A,T12,A)")'SCC Iter','Error'
        end if

        iSccIter = 1
        tStopSCC = .false.
        lpSCC: do while (iSCCiter <= maxSccIter)

          dPotential%intAtom(:,:) = 0.0_dp
          dPotential%intShell(:,:,:) = 0.0_dp
          dPotential%intBlock(:,:,:,:) = 0.0_dp

          !if (tDFTBU .or. allocated(onsMEs)) then
          !  dPotential%orbitalBlock(:,:,:,:) = 0.0_dp
          !end if
          if (tSccCalc .and. iSCCiter>1) then
            call sccCalc%updateCharges(env, dqIn + dqNonVariational, orb, species)
            call sccCalc%updateShifts(env, orb, species, neighbourList%iNeighbour, img2CentCell)
            ! Note, should ommit external charge and constraint shifts
            call sccCalc%getShiftPerAtom(dPotential%intAtom(:,1), isOnlyInternalShifts=.true.)
            call sccCalc%getShiftPerL(dPotential%intShell(:,:,1))

            if (allocated(spinW)) then
              call getChargePerShell(dqIn + dqNonVariational, orb, species, dqPerShell)
              call getSpinShift(shellPot, dqPerShell, species, orb, spinW)
              dPotential%intShell(:,:,:) = dPotential%intShell + shellPot
            end if

            if (allocated(thirdOrd)) then
              atomPot(:,:) = 0.0_dp
              shellPot(:,:,:) = 0.0_dp
              call thirdOrd%getdShiftdQ(atomPot(:,1), shellPot(:,:,1), species, neighbourList,&
                  & dqIn + dqNonVariational, img2CentCell, orb)
              dPotential%intAtom(:,1) = dPotential%intAtom(:,1) + atomPot(:,1)
              dPotential%intShell(:,:,1) = dPotential%intShell(:,:,1) + shellPot(:,:,1)
            end if

            !if (tDFTBU) then
            ! ! note the derivatives of both FLL and pSIC are the same (pSIC, i.e. case 2 in module)
            ! call getDftbUShift(dPotential%orbitalBlock, dqBlockIn+dqNonVariationalBlock, species,&
            !     & orb, 2, UJ, nUJ, niUJ, iUJ)
            !end if
            !if (allocated(onsMEs)) then
            !  ! onsite corrections
            !  call addOnsShift(dPotential%orbitalBlock, dPotential%iOrbitalBlock,&
            !      & dqBlockIn + dqNonVariationalBlock, dummy, onsMEs, species, orb)
            !end if

          end if

          call totalShift(dPotential%intShell,dPotential%intAtom, orb, species)
          call totalShift(dPotential%intBlock,dPotential%intShell, orb, species)
          if (tSccCalc) then
            call totalShift(dPotential%intBlock(:,:,:,1:1), vdgamma, orb, species)
          end if
          !if (tDFTBU .or. allocated(onsMEs)) then
          !  dPotential%intBlock(:,:,:,:) = dPotential%intBlock + dPotential%orbitalBlock
          !end if
          dPotential%intBlock(:,:,:,:) = dPotential%intBlock + dPotential%extBlock

          if (tSccCalc) then
            ! add the (Delta q) * d gamma / dx term
            ! and add gamma * d (Delta q) / dx
            call addShift(env, sOmega(:,:,2), over, nNeighbourSK, neighbourList%iNeighbour,&
                & species, orb, iSparseStart, nAtom, img2CentCell, dpotential%intBlock,&
                & isInputZero=.true.)
          end if

          dHam(:,:) = 0.0_dp
          dHam(:,1) = dH0(:,iCart)
          if (tSccCalc) then
            dHam(:,:) = dHam + sum(sOmega, dim=3)
          end if

          if (nSpin > 1) then
            dHam(:,:) = 2.0_dp * dHam(:,:)
            if (allocated(idHam)) then
              idHam(:,:) = 2.0_dp * idHam(:,:)
            end if
          end if

          call qm2ud(dHam)
          if (allocated(idHam)) then
            call qm2ud(idHam)
          end if

          dRho(:,:) = 0.0_dp
          if (allocated(idRho)) then
            idRho(:,:) = 0.0_dp
          end if

          ! evaluate derivative of density matrix
          if (allocated(eigVecsReal)) then
            drho(:,:) = 0.0_dp
            do iKS = 1, parallelKS%nLocalKS

              iS = parallelKS%localKS(2, iKS)

              if (allocated(dRhoOut)) then
                ! replace with case that will get updated in dRhoReal
                dRhoOutSqr(:,:,iS) = dRhoInSqr(:,:,iS)
              end if

              dOverWork = dOver(:,iCart)
              dEiTmp(:,:,:) = 0.0_dp
              call dRhoReal(env, dHam, dOverWork, neighbourList, nNeighbourSK,&
                  & iSparseStart, img2CentCell, denseDesc, iKS, parallelKS, nFilled(:,1:1),&
                  & nEmpty(:,1:1), eigVecsReal, eCiReal, eigVals, filling, Ef, tempElec, orb,&
                  & drho(:,iS), dRhoOutSqr, hybridXc, over, nNeighbourLC, degenTransform(iKS),&
                  & species, dEiTmp, dPsiReal, coord, errStatus, omega, isHelical, maxFill=maxFill)
              @:PROPAGATE_ERROR(errStatus)
              dEi(:, :, :, iCart, iDeriv) = dEiTmp
            end do

          elseif (nSpin > 2) then

            @:RAISE_ERROR(errStatus, -1, "Shouldn't be here in dxAtom routine - Pauli missing")

            do iKS = 1, parallelKS%nLocalKS

              iK = parallelKS%localKS(1, iKS)

!              call dRhoPauli(env, dHam, idHam, dOver(:,igotgotCart), neighbourList, nNeighbourSK,&
!                  & iSparseStart, img2CentCell, denseDesc, parallelKS, nFilled(:, iK),&
!                  & nEmpty(:, iK), eigvecsCplx, eigVals, Ef, tempElec, orb, dRho, idRho, kPoint,&
!                  & kWeight, iCellVec, cellVec, iKS, iCart,&
!                #:if WITH_SCALAPACK
!                  & desc,&
!                #:endif
!                  & dEi, dPsiCmplx, iAt)

            end do

          else

            @:RAISE_ERROR(errStatus, -1, "Shouldn't be here in dxAtom routine")

            drho(:,:) = 0.0_dp
            do iKS = 1, parallelKS%nLocalKS

              iK = parallelKS%localKS(1, iKS)
              iS = parallelKS%localKS(2, iKS)

            end do

          end if

        #:if WITH_SCALAPACK
          ! Add up and distribute density matrix contributions from each group
          call mpifx_allreduceip(env%mpi%globalComm, dRho, MPI_SUM)
        #:endif

          call ud2qm(dRho)

          if (allocated(idRho)) then
            idRho(:,:) = maxFill * drho
            call ud2qm(idRho)
          end if

          dqOut(:,:,:,iCart, iDeriv) = 0.0_dp
          do iS = 1, nSpin
            call mulliken(env, dqOut(:,:, iS, iCart, iDeriv), over, dRho(:,iS), orb, &
                & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
            !if (allocated(dftbU) .or. allocated(onsMEs)) then
            !  dqBlockOut(:,:,:,iS) = 0.0_dp
            !  call mulliken(env, dqBlockOut(:,:,:,iS), over, dRho(:,iS), orb,&
            !      & neighbourList%iNeighbour,&
            !      & nNeighbourSK, img2CentCell, iSparseStart)
            !end if
          end do

          if (tSccCalc) then

            if (allocated(hybridXc)) then
              dqDiffRed(:) = dRhoOut - dRhoIn
            else
              dqOutRed(:) = 0.0_dp
              call OrbitalEquiv_reduce(dqOut(:,:,:,iCart, iDeriv), iEqOrbitals, orb,&
                  & dqOutRed(:nIneqMixElements))
              !if (allocated(dftbU)) then
              !  call AppendBlockReduced(dqBlockOut, iEqBlockDFTBU, orb, dqOutRed)
              !end if
              !if (allocated(onsMEs)) then
              !  call AppendBlockReduced(dqBlockOut, iEqBlockOnsite, orb, dqOutRed)
              !end if
              dqDiffRed(:) = dqOutRed - dqInpRed
            end if
            sccErrorQ = maxval(abs(dqDiffRed))

            if (maxSccIter > 1) then
              write(stdOut,"(1X,I0,T10,E20.12)")iSccIter, sccErrorQ
            end if
            tConverged = (sccErrorQ < sccTol)

            if ((.not. tConverged) .and. iSCCiter /= maxSccIter) then
              if (iSccIter == 1) then
                if (allocated(hybridXc)) then
                  dRhoIn(:) = dRhoOut
                #:if WITH_SCALAPACK
                  call denseMullikenRealBlacs(env, parallelKS, denseDesc, dRhoInSqr, SSqrReal, dqIn)
                #:else
                  call denseMullikenReal(dRhoInSqr, SSqrReal, denseDesc%iAtomStart, dqIn)
                #:endif
                else
                  dqIn(:,:,:) = dqOut(:,:,:,iCart, iDeriv)
                  dqInpRed(:) = dqOutRed
                  !if (allocated(dftbU) .or. allocated(onsMEs)) then
                  !  dqBlockIn(:,:,:,:) = dqBlockOut
                  !end if
                end if

              else

                if (allocated(hybridXc)) then
                  call chrgMixerReal%mix(dRhoIn, dqDiffRed)
                #:if WITH_SCALAPACK
                  call denseMullikenRealBlacs(env, parallelKS, denseDesc, dRhoInSqr, SSqrReal, dqIn)
                #:else
                  call denseMullikenReal(dRhoInSqr, SSqrReal, denseDesc%iAtomStart, dqIn)
                #:endif
                else
                  call chrgMixerReal%mix(dqInpRed, dqDiffRed)
                #:if WITH_SCALAPACK
                  ! Synchronise charges in order to avoid mixers that store a history drifting apart
                  call mpifx_bcast(env%mpi%globalComm, dqInpRed)
                #:endif

                  call OrbitalEquiv_expand(dqInpRed(:nIneqMixElements), iEqOrbitals, orb, dqIn)

                  !if (allocated(dftbU) .or. allocated(onsMEs)) then
                    !dqBlockIn(:,:,:,:) = 0.0_dp
                    !if (allocated(dftbU)) then
                      !call dftbU%expandBlock(dqInpRed, iEqBlockDFTBU, orb, dqBlockIn, species,&
                      !    & orbEquiv=iEqOrbitals)
                    !else
                    !  call Onsblock_expand(dqInpRed, iEqBlockOnSite, orb, dqBlockIn,&
                    !      & orbEquiv=iEqOrbitals)
                    !end if
                  !end if
                end if

              end if

              if (allocated(hybridXc)) then
                call ud2qm(dqIn)
              end if

            end if
          else
            tConverged = .true.
          end if

          if (tConverged) then
            exit lpSCC
          end if

          if (allocated(spinW)) then
            dqPerShell = 0.0_dp
            do jAt = 1, nAtom
              iSp = species(jAt)
              do iSh = 1, orb%nShell(iSp)
                dqPerShell(iSh,jAt,:nSpin) = dqPerShell(iSh,jAt,:nSpin) +&
                    & sum(dqIn(orb%posShell(iSh,iSp): orb%posShell(iSh+1,iSp)-1,jAt,:nSpin),dim=1)
              end do
            end do

          end if

          iSccIter = iSccIter +1

        end do lpSCC

        dqOut(:, :, :, iCart, iDeriv) = dqOut(:, :, :, iCart, iDeriv) + dqNonVariational

      end do lpCart

    end do lpDerivs

  #:if WITH_SCALAPACK
    call mpifx_allreduceip(env%mpi%globalComm, dEi, MPI_SUM)
  #:endif

!    write(stdOut, "(A)")'dEi/d (eV / AA)'
!    write(stdOut,"(T16,A,T32,A,T48,A)")direction
!    do iS = 1, nSpin
!      if (allocated(wrtWhichCharges)) then
!        do iDeriv = 1, size(wrtWhichCharges)
!          iAt = wrtWhichCharges(iDeriv)
!          write(stdOut, "(1X,A,I0)")'dEi/d At.', iAt
!          do iOrb = 1, size(dEi,dim=1)
!            write(stdOut, "(3F16.8)") dEi(iOrb, 1, iS, :, iDeriv) * AA__Bohr
!          end do
!        end do
!      else
!        do iAt = 1, nAtom
!          write(stdOut, "(1X,A,I0)")'dEi/d At.', iAt
!          do iOrb = 1, size(dEi,dim=1)
!            write(stdOut, "(3F16.8)") dEi(iOrb, 1, iS, :, iAt) * AA__Bohr
!          end do
!        end do
!      end if
!
!    end do

    if (isBandWritten) then
      if (allocated(wrtWhichCharges)) then
        do iDeriv = 1, size(wrtWhichCharges)
          iAt = wrtWhichCharges(iDeriv)
          do iCart = 1, 3
            write(tmpStr, "(A,1X,I0,1X,A,1X)")"Atom", iAt, "d"//direction(iCart)
            if (iAt == wrtWhichCharges(1) .and. iCart == 1) then
              call writeDerivBandOut("dBand_dx.out", dEi(:, :, :, iCart, iDeriv), [1.0_dp],.false.,&
                  & trim(tmpStr))
            else
              call writeDerivBandOut("dBand_dx.out", dEi(:, :, :, iCart, iDeriv), [1.0_dp], .true.,&
                  & trim(tmpStr))
            end if
          end do
        end do
      else
        do iAt = 1, nAtom
          do iCart = 1, 3
            write(tmpStr, "(A,1X,I0,1X,A,1X)")"Atom", iAt, "d"//direction(iCart)
            if (iAt == 1 .and. iCart == 1) then
              call writeDerivBandOut("dBand_dx.out", dEi(:, :, :, iCart, iAt), [1.0_dp],.false.,&
                  & trim(tmpStr))
            else
              call writeDerivBandOut("dBand_dx.out", dEi(:, :, :, iCart, iAt), [1.0_dp], .true.,&
                  & trim(tmpStr))
            end if
          end do
        end do
      end if
    end if

    if (tPrintMulliken) then

      write(stdOut, *)
      write(stdOut, "(A)")'Derivatives of atomic Mulliken charges with atom positions'
      if (allocated(wrtWhichCharges)) then
        do iDeriv = 1, size(wrtWhichCharges)
          iAt = wrtWhichCharges(iDeriv)
          write(stdOut,"(1X,A,I0,T10,A,T26,A,T42,A)")"At",iAt,"x","y","z"
          do iS = 1, nSpin
            do jAt = 1, nAtom
              write(stdOut, "(1X,I0,T4,3F16.8)")jAt, -sum(dqOut(:,jAt,iS,:,iDeriv),dim=1)&
                  & * AA__Bohr
            end do
            write(stdOut, *)
          end do
        end do
      else
        do iAt = 1, nAtom
          write(stdOut,"(A,I0,T10,A,T26,A,T42,A)")"At",iAt,"x","y","z"
          do iS = 1, nSpin
            do jAt = 1, nAtom
              write(stdOut, "(1X,I0,T4,3F16.8)")jAt, -sum(dqOut(:,jAt,iS,:,iAt),dim=1)&
                  & * AA__Bohr
            end do
            write(stdOut, *)
          end do
        end do
      end if

      if (tWriteDetailedOut) then
        write(fdDetailedOut, "(A)")'Derivatives of atomic Mulliken charges with atom positions'
        if (allocated(wrtWhichCharges)) then
          do iDeriv = 1, size(wrtWhichCharges)
            iAt = wrtWhichCharges(iDeriv)
            write(fdDetailedOut,"(1X,A,I0,T10,A,T26,A,T42,A)")"At",iAt,"x","y","z"
            do iS = 1, nSpin
              do jAt = 1, nAtom
                write(fdDetailedOut, "(1X,I0,T4,3F16.8)")jAt, -sum(dqOut(:,jAt,iS,:,iDeriv),dim=1)&
                    & * AA__Bohr
              end do
              write(fdDetailedOut, *)
            end do
          end do
        else
          do iAt = 1, nAtom
            write(fdDetailedOut,"(A,I0,T10,A,T26,A,T42,A)")"At",iAt,"x","y","z"
            do iS = 1, nSpin
              do jAt = 1, nAtom
                write(fdDetailedOut, "(1X,I0,T4,3F16.8)")jAt, -sum(dqOut(:,jAt,iS,:,iAt),dim=1)&
                    & * AA__Bohr
              end do
              write(fdDetailedOut, *)
            end do
          end do
        end if
      end if

      if (present(dqdx)) then
        dqdx(:,:,:) = 0.0_dp
        if (allocated(wrtWhichCharges)) then
          do iDeriv = 1, size(wrtWhichCharges)
            iAt = wrtWhichCharges(iDeriv)
            do iCart = 1, 3
              do iS = 1, nSpin
                do jAt = 1, nAtom
                  dqdx(jAt, iCart, iDeriv) = dqdx(jAt, iCart, iDeriv)&
                      & + sum(dqOut(:, jAt, iS, iCart, iDeriv))
                end do
              end do
            end do
          end do
        else
          do iAt = 1, nAtom
            do iCart = 1, 3
              do iS = 1, nSpin
                do jAt = 1, nAtom
                  dqdx(jAt, iCart, iAt) = dqdx(jAt, iCart, iAt)&
                      & + sum(dqOut(:, jAt, iS, iCart, iAt))
                end do
              end do
            end do
          end do
        end if
      end if

      if (allocated(wrtWhichCharges)) then
        allocate(bornCharges(3, 3, size(wrtWhichCharges)), source=0.0_dp)
        do iDeriv = 1, size(wrtWhichCharges)
          iAt = wrtWhichCharges(iDeriv)
          do iCart = 1, 3
            dDipole(:) = 0.0_dp
            do jCart = 1, 3
              dDipole(jCart) = sum(sum(dqOut(:,:nAtom,1,iCart,iDeriv), dim=1) * coord(jCart,:nAtom))
            end do
            dDipole(iCart) = dDipole(iCart) + sum(qOrb(:,iAt,1) - q0(:,iAt,1))
            bornCharges(:, iCart, iDeriv) = dDipole
          end do
        end do
      else
        allocate(bornCharges(3, 3, nAtom), source=0.0_dp)
        do iAt = 1, nAtom
          do iCart = 1, 3
            dDipole(:) = 0.0_dp
            do jCart = 1, 3
              dDipole(jCart) = sum(sum(dqOut(:,:nAtom,1,iCart,iAt), dim=1) * coord(jCart,:nAtom))
            end do
            dDipole(iCart) = dDipole(iCart) + sum(qOrb(:,iAt,1) - q0(:,iAt,1))
            bornCharges(:, iCart, iAt) = dDipole
          end do
        end do
      end if

      if (tWriteDetailedOut) call writeBorn(fdDetailedOut, bornCharges, wrtWhichCharges)
      call writeBorn(stdOut, bornCharges, wrtWhichCharges)

      if (isAutotestWritten) then
        open(newunit=fdResults, file=autoTestTagFile, position="append")
        call taggedWriter%write(fdResults, tagLabels%dqdx, dqOut)
        close(fdResults)
      end if
      if (isTagResultsWritten) then
        open(newunit=fdResults, file=taggedResultsFile, position="append")
        call taggedWriter%write(fdResults, tagLabels%dqdx, dqOut)
        call taggedWriter%write(fdResults, tagLabels%borncharges, bornCharges)
        close(fdResults)
      end if

    end if

  end subroutine dxAtom


  ! Internal routines

  !> Evaluates response, given the external perturbation
  subroutine response(env, parallelKS, dPotential, nAtom, orb, species, neighbourList,&
      & nNeighbourSK, img2CentCell, iSparseStart, denseDesc, over, iEqOrbitals, sccCalc, sccTol,&
      & isSccConvRequired, maxSccIter, chrgMixerReal, nMixElements, nIneqMixElements, dqIn, dqOut,&
      & hybridXc, nNeighbourCam, sSqrReal, dRhoInSqr, dRhoOutSqr, dRhoIn, dRhoOut, nSpin, maxFill,&
      & spinW, thirdOrd, dftbU, iEqBlockDftbu, onsMEs, iEqBlockOnSite, dqBlockIn, dqBlockOut,&
      & eigVals, degenTransform, dEi, dEf, filling, Ef, isEfFixed, dHam, idHam, dRho, idRho,&
      & tempElec, tMetallic, neFermi, nFilled, nEmpty, kPoint, kWeight, cellVec, iCellVec,&
      & eigVecsReal, eigVecsCplx, dPsiReal, dPsiCmplx, coord, errStatus, omega, dDipole, isHelical,&
      & eta)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> The k-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    ! derivative of potentials
    type(TPotentials), intent(inout) :: dPotential

    !> Charge mixing object
    class(TMixerReal), intent(inout), allocatable :: chrgMixerReal

    !> Nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
    !> blocks (if a DFTB+U or onsite correction calculation)
    integer, intent(in) :: nMixElements

    !> Nr. of inequivalent charges
    integer, intent(in) :: nIneqMixElements

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> SCC module internal variables
    type(TScc), intent(inout), allocatable :: sccCalc

    !> Derivative of sparse Hamiltonian
    real(dp), intent(inout) :: dHam(:,:)

    !> Derivative of imaginary part of sparse Hamiltonian
    real(dp), intent(inout), allocatable :: idHam(:,:)

    !> Derivative of sparse density matrix
    real(dp), intent(inout) :: dRho(:,:)

    !> Derivative of imaginary part of sparse density matrix
    real(dp), intent(inout), allocatable :: idRho(:,:)

    !> Maximal number of SCC iterations
    integer, intent(in) :: maxSccIter

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Number of neighbours for each of the atoms for the exchange contributions of CAM functionals
    integer, intent(in), allocatable :: nNeighbourCam(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Chemical species
    integer, intent(in) :: species(:)

    !> Spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> Third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Is there a finite density of states at the Fermi energy
    logical, intent(in) :: tMetallic(:,:)

    !> Number of electrons at the Fermi energy (if metallic)
    real(dp), allocatable, intent(in) :: neFermi(:)

    !> Tolerance for SCC convergence
    real(dp), intent(in) :: sccTol

    !> Use converged derivatives of charges
    logical, intent(in) :: isSccConvRequired

    !> Maximum allowed number of electrons in a single particle state
    real(dp), intent(in) :: maxFill

    !> Sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Equivalence relations between orbitals
    integer, intent(in), allocatable :: iEqOrbitals(:,:,:)

    !> For transformation in the  case of degeneracies
    type(TRotateDegen), allocatable, intent(inout) :: degenTransform(:)

    !> Derivative of density matrix
    real(dp), target, allocatable, intent(inout) :: dRhoIn(:)

    !> Derivative of density matrix
    real(dp), target, allocatable, intent(inout) :: dRhoOut(:)

    !> Delta density matrix response for hybrid xc-functional calculations
    real(dp), pointer :: dRhoOutSqr(:,:,:)

    !> Delta density matrix input for hybrid xc-functional calculations
    real(dp), pointer :: dRhoInSqr(:,:,:)

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> Equivalence mapping for dual charge blocks
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> Levels with at least partial filling
    integer, intent(in) :: nFilled(:,:)

    !> Levels that are at least partially empty
    integer, intent(in) :: nEmpty(:,:)

    !> Derivatives of Mulliken charges, if required
    real(dp), intent(inout) :: dqOut(:,:,:)

    !> Derivative of eigenvalues, if requested
    real(dp), intent(inout), allocatable :: dEi(:,:,:)

    !> Derivatives of Mulliken charges, if required
    real(dp), intent(inout) :: dqIn(:,:,:)

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Filling of levels
    real(dp), intent(in) :: filling(:,:,:)

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Is a fixed Fermi energy used
    logical, intent(in) :: isEfFixed

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in), allocatable :: eigVecsReal(:,:,:)

    !> Ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Derivative of Fermi energy
    real(dp), intent(inout), allocatable :: dEf(:)

    !> Derivative of block charges (input)
    real(dp), allocatable, intent(inout) :: dqBlockIn(:,:,:,:)

    !> Derivative of block charges (output)
    real(dp), allocatable, intent(inout) :: dqBlockOut(:,:,:,:)

    !> Onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Data for hybrid functional calculation
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Square matrix for overlap (if needed)
    real(dp), allocatable, intent(inout) :: sSqrReal(:,:)

    !> The k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Derivative of single particle wavefunctions (real case), if needed
    real(dp), allocatable, intent(inout) :: dPsiReal(:,:,:)

    !> Derivative of single particle wavefunctions (complex case), if needed
    complex(dp), allocatable, intent(inout) :: dPsiCmplx(:,:,:,:)

    !> Coordinates of all atoms including images
    real(dp), intent(in) :: coord(:,:)

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    !> Finite frequency, if relevant, zero otherwise
    real(dp), intent(in) :: omega

    !> Derivative of dipole
    real(dp), intent(out), optional :: dDipole(:)

    !> Is the geometry helical
    logical, intent(in), optional :: isHelical

    !> Small complex value for frequency dependent
    complex(dp), intent(in), optional :: eta

    logical :: tSccCalc, tConverged
    integer :: iSccIter
    real(dp), allocatable :: shellPot(:,:,:), atomPot(:,:)
    real(dp), allocatable :: dummy(:,:,:,:)

    real(dp) :: dqInpRed(nMixElements), dqOutRed(nMixElements), sccErrorQ
    real(dp) :: dqPerShell(orb%mShell, nAtom, nSpin)

    integer :: iAt, iKS, iK, iS, iSh, iSp
    real(dp), allocatable :: dRhoExtra(:,:), idRhoExtra(:,:)
    real(dp) :: dqDiffRed(nMixElements)

    real(dp), allocatable :: dqInBackup(:,:,:), dqBlockInBackup(:,:,:,:), dRhoInBackup(:)
    real(dp) :: quietNan

    ! temporary dumy variables for perturbation routine, as overlap not changing
    real(dp), allocatable :: dOver(:), eCiReal(:,:,:)

    tSccCalc = allocated(sccCalc)

    @:ASSERT(abs(omega) <= epsilon(0.0_dp) .or. present(eta))

    if (tSccCalc.and..not.isSccConvRequired) then
      ! store back-ups of the input charges/density matrix in case the iteration fails
      dqInBackup = dqIn
      if (allocated(dqBlockIn)) then
        dqBlockInBackup = dqBlockIn
      end if
      if (allocated(dRhoIn)) then
        dRhoInBackup = dRhoIn
      end if
    end if

    if (allocated(spinW) .or. allocated(thirdOrd)) then
      allocate(shellPot(orb%mShell, nAtom, nSpin))
    end if
    if (allocated(thirdOrd)) then
      allocate(atomPot(nAtom, nSpin))
    end if

    if (any(tMetallic)) then
      allocate(dRhoExtra(size(over),nSpin))
      if (allocated(idHam)) then
        allocate(idRhoExtra(size(over),nSpin))
      end if
    end if

    if (tSccCalc) then
      dqInpRed(:) = 0.0_dp
      dqPerShell(:,:,:) = 0.0_dp
      if (allocated(hybridXc)) then
        dRhoIn(:) = 0.0_dp
        dRhoOut(:) = 0.0_dp
      end if
    end if

    if (present(dDipole)) then
      dDipole(:) = 0.0_dp
    end if

    if (tSccCalc) then
      call chrgMixerReal%reset(size(dqInpRed))
    end if

    if (abs(omega) > epsilon(0.0_dp)) then
      write(stdOut, "(1X,A)")"Frequency dependant response calculation"
      write(stdOut, format2U)"  omega driving frequency", omega, ' H ', omega * Hartree__eV, ' eV'
    else
      write(stdOut, "(1X,A)")"Static response calculation"
    end if


    if (tSccCalc .and. maxSccIter > 1) then
      write(stdOut,"(1X,A,T12,A)")'SCC Iter','Error'
    end if

    lpSCC: do iSccIter = 1, maxSccIter

      dPotential%intAtom(:,:) = 0.0_dp
      dPotential%intShell(:,:,:) = 0.0_dp
      dPotential%intBlock(:,:,:,:) = 0.0_dp

      if (allocated(dftbU) .or. allocated(onsMEs)) then
        dPotential%orbitalBlock(:,:,:,:) = 0.0_dp
      end if

      if (tSccCalc .and. iSCCiter>1) then
        call sccCalc%updateCharges(env, dqIn, orb, species)
        call sccCalc%updateShifts(env, orb, species, neighbourList%iNeighbour, img2CentCell)

        ! Note, should ommit external charge and constraint shifts
        call sccCalc%getShiftPerAtom(dPotential%intAtom(:,1), isOnlyInternalShifts=.true.)
        call sccCalc%getShiftPerL(dPotential%intShell(:,:,1))

        if (allocated(spinW)) then
          call getChargePerShell(dqIn, orb, species, dqPerShell)
          shellPot(:,:,:) = 0.0_dp
          call getSpinShift(shellPot(:,:,2:), dqPerShell(:,:,2:), species, orb, spinW)
          dPotential%intShell(:,:,2:) = dPotential%intShell(:,:,2:) + shellPot(:,:,2:)
        end if

        if (allocated(thirdOrd)) then
          atomPot(:,:) = 0.0_dp
          shellPot(:,:,:) = 0.0_dp
          call thirdOrd%getdShiftdQ(atomPot(:,1), shellPot(:,:,1), species, neighbourList, dqIn,&
              & img2CentCell, orb)
          dPotential%intAtom(:,1) = dPotential%intAtom(:,1) + atomPot(:,1)
          dPotential%intShell(:,:,1) = dPotential%intShell(:,:,1) + shellPot(:,:,1)
        end if

        if (allocated(dftbU)) then
          ! note the derivatives of both FLL and pSIC potentials are the same (i.e. pSIC)
          call dftbU%getDftbUShift(dPotential%orbitalBlock, dqBlockIn, species, orb,&
              & plusUFunctionals%pSIC)
        end if
        if (allocated(onsMEs)) then
          ! onsite corrections
          call addOnsShift(dPotential%orbitalBlock, dPotential%iOrbitalBlock, dqBlockIn, dummy,&
              & onsMEs, species, orb)
        end if

      end if

      call totalShift(dPotential%intShell,dPotential%intAtom, orb, species)
      call totalShift(dPotential%intBlock,dPotential%intShell, orb, species)
      if (allocated(dftbU) .or. allocated(onsMEs)) then
        dPotential%intBlock(:,:,:,:) = dPotential%intBlock + dPotential%orbitalBlock
      end if
      dPotential%intBlock(:,:,:,:) = dPotential%intBlock + dPotential%extBlock

      dHam(:,:) = 0.0_dp
      call addShift(env, dHam, over, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
          & iSparseStart, nAtom, img2CentCell, dPotential%intBlock, isInputZero=.true.)

      if (nSpin > 1) then
        dHam(:,:) = 2.0_dp * dHam
        if (allocated(idHam)) then
          idHam(:,:) = 2.0_dp * idHam
        end if
      end if
      call qm2ud(dHam)
      if (allocated(idHam)) then
        call qm2ud(idHam)
      end if

      dRho(:,:) = 0.0_dp
      if (allocated(idRho)) then
        idRho(:,:) = 0.0_dp
      end if

      ! evaluate derivative of density matrix
      if (allocated(eigVecsReal)) then

        do iKS = 1, parallelKS%nLocalKS

          iS = parallelKS%localKS(2, iKS)

          if (allocated(dRhoOut)) then
            ! replace with case that will get updated in dRhoReal
            dRhoOutSqr(:,:,iS) = dRhoInSqr(:,:,iS)
          end if

          call dRhoReal(env, dHam, dOver, neighbourList, nNeighbourSK, iSparseStart, img2CentCell,&
              & denseDesc, iKS, parallelKS, nFilled, nEmpty, eigVecsReal, eCiReal, eigVals,&
              & filling, Ef, tempElec, orb, dRho(:,iS), dRhoOutSqr, hybridXc, over, nNeighbourCam,&
              & degenTransform(iKS), species, dEi, dPsiReal, coord, errStatus, omega, isHelical,&
              & eta=eta)
          if (errStatus%hasError()) then
            exit
          end if
        end do
        @:PROPAGATE_ERROR(errStatus)

      elseif (nSpin > 2) then

        do iKS = 1, parallelKS%nLocalKS

          iK = parallelKS%localKS(1, iKS)

          call dRhoPauli(env, dHam, idHam, neighbourList, nNeighbourSK, iSparseStart,&
              & img2CentCell, denseDesc, parallelKS, nFilled, nEmpty, eigvecsCplx, eigVals, Ef,&
              & tempElec, orb, dRho, idRho, kPoint, kWeight, iCellVec, cellVec, iKS,&
              & degenTransform(iKS), species, coord, dEi, dPsiCmplx, errStatus, omega, isHelical,&
              & eta=eta)
          if (errStatus%hasError()) then
            exit
          end if
        end do
        @:PROPAGATE_ERROR(errStatus)

        ! adjustment from Pauli to charge/spin
        dRho(:,:) = 2.0_dp * dRho
        if (allocated(idRho)) then
          idRho(:,:) = 2.0_dp * idRho
        end if

      else

        do iKS = 1, parallelKS%nLocalKS

          iK = parallelKS%localKS(1, iKS)

          call dRhoCmplx(env, dHam, neighbourList, nNeighbourSK, iSparseStart, img2CentCell,&
              & denseDesc, parallelKS, nFilled, nEmpty, eigvecsCplx, eigVals, Ef, tempElec, orb,&
              & dRho, kPoint, kWeight, iCellVec, cellVec, iKS, degenTransform(iKS), species,&
              & coord, dEi, dPsiCmplx, errStatus, omega, isHelical, eta=eta)
          if (errStatus%hasError()) then
            exit
          end if
        end do
        @:PROPAGATE_ERROR(errStatus)

      end if

    #:if WITH_SCALAPACK
      ! Add up and distribute density matrix contributions from each group
      call mpifx_allreduceip(env%mpi%globalComm, dRho, MPI_SUM)
    #:endif

      if (allocated(dEf)) then
        dEf(:) = 0.0_dp
      end if

      if (any(tMetallic) .and. .not. isEfFixed) then
        ! correct for Fermi level shift for q=0 fields
        dRhoExtra(:,:) = 0.0_dp
        if (allocated(idRhoExtra)) then
          idRhoExtra(:,:) = 0.0_dp
        end if
        do iKS = 1, parallelKS%nLocalKS
          iK = parallelKS%localKS(1, iKS)
          iS = parallelKS%localKS(2, iKS)

          if (.not.tMetallic(iS,iK)) then
            cycle
          end if

          dqOut(:,:,iS) = 0.0_dp
          call mulliken(env, dqOut(:,:,iS), over, dRho(:,iS), orb, &
              & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)

          dEf(iS) = -sum(dqOut(:,:, iS)) / neFermi(iS)

          if (abs(dEf(iS)) > 10.0_dp*epsilon(1.0_dp)) then
            ! Fermi level changes, so need to correct for the change in the number of charges

            if (allocated(eigVecsReal)) then

              ! real case, no k-points
              call dRhoFermiChangeReal(dRhoExtra(:, iS), env, maxFill, parallelKS, iKS,&
                  & neighbourList, nNeighbourSK, img2CentCell, iSparseStart, dEf, Ef, nFilled,&
                  & nEmpty, eigVecsReal, orb, denseDesc, tempElec, eigVals, dRhoOutSqr, species,&
                  & coord, isHelical)

            elseif (nSpin > 2) then

              ! two component wavefunction cases
              call dRhoFermiChangePauli(dRhoExtra, idRhoExtra, env, parallelKS, iKS,&
                  & kPoint, kWeight, iCellVec, cellVec, neighbourList, nNEighbourSK,&
                  & img2CentCell, iSparseStart, dEf, Ef, nFilled, nEmpty, eigVecsCplx, orb,&
                  & denseDesc, tempElec, eigVals, species, coord, errStatus, isHelical)
              @:PROPAGATE_ERROR(errStatus)

            else

              ! Complex case with k-points
              call dRhoFermiChangeCmplx(dRhoExtra, env, maxFill, parallelKS, iKS, kPoint, kWeight,&
                  & iCellVec, cellVec, neighbourList, nNEighbourSK, img2CentCell, iSparseStart,&
                  & dEf, Ef, nFilled, nEmpty, eigVecsCplx, orb, denseDesc, tempElec, eigVals,&
                  & species, coord, isHelical)

            end if

          end if

        end do

        if (nSpin > 2) then
          ! adjustment from Pauli to charge/spin
          dRhoExtra(:,:) = 2.0_dp * dRhoExtra
          if (allocated(idRhoExtra)) then
            idRhoExtra(:,:) = 2.0_dp * idRhoExtra
          end if
        end if

      #:if WITH_SCALAPACK
        ! Add up and distribute density matrix contribution from each group
        call mpifx_allreduceip(env%mpi%globalComm, dRhoExtra, MPI_SUM)
      #:endif
        dRho(:,:) = dRho + dRhoExtra

      end if

      dRho(:,:) = maxFill * drho
      if (allocated(dRhoOut)) then
        dRhoOut(:) = maxFill * dRhoOut
      end if
      call ud2qm(dRho)

      if (allocated(idRho)) then
        idRho(:,:) = maxFill * drho
        call ud2qm(idRho)
      end if

      dqOut(:,:,:) = 0.0_dp
      do iS = 1, nSpin
        call mulliken(env, dqOut(:,:,iS), over, dRho(:,iS), orb, &
            & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
        if (allocated(dftbU) .or. allocated(onsMEs)) then
          dqBlockOut(:,:,:,iS) = 0.0_dp
          call mulliken(env, dqBlockOut(:,:,:,iS), over, dRho(:,iS), orb, neighbourList%iNeighbour,&
              & nNeighbourSK, img2CentCell, iSparseStart)
        end if
      end do

      if (tSccCalc) then

        if (allocated(hybridXc)) then
          dqDiffRed(:) = dRhoOut - dRhoIn
        else
          dqOutRed(:) = 0.0_dp
          call OrbitalEquiv_reduce(dqOut, iEqOrbitals, orb, dqOutRed(:nIneqMixElements))
          if (allocated(dftbU)) then
            call AppendBlockReduced(dqBlockOut, iEqBlockDFTBU, orb, dqOutRed)
          end if
          if (allocated(onsMEs)) then
            call AppendBlockReduced(dqBlockOut, iEqBlockOnsite, orb, dqOutRed)
          end if
          dqDiffRed(:) = dqOutRed - dqInpRed
        end if
        sccErrorQ = maxval(abs(dqDiffRed))

        if (maxSccIter > 1) then
          write(stdOut,"(1X,I0,T10,E20.12)")iSCCIter, sccErrorQ
        end if
        tConverged = (sccErrorQ < sccTol)

        if ((.not. tConverged) .and. iSCCiter /= maxSccIter) then
          if (iSCCIter == 1) then
            if (allocated(hybridXc)) then
              dRhoIn(:) = dRhoOut
            #:if WITH_SCALAPACK
              call denseMullikenRealBlacs(env, parallelKS, denseDesc, dRhoInSqr, SSqrReal, dqIn)
            #:else
              call denseMullikenReal(dRhoInSqr, SSqrReal, denseDesc%iAtomStart, dqIn)
            #:endif
            else
              dqIn(:,:,:) = dqOut
              dqInpRed(:) = dqOutRed
              if (allocated(dftbU) .or. allocated(onsMEs)) then
                dqBlockIn(:,:,:,:) = dqBlockOut
              end if
            end if

          else

            if (allocated(hybridXc)) then
              call chrgMixerReal%mix(dRhoIn, dqDiffRed)
            #:if WITH_SCALAPACK
              call denseMullikenRealBlacs(env, parallelKS, denseDesc, dRhoInSqr, SSqrReal, dqIn)
            #:else
              call denseMullikenReal(dRhoInSqr, SSqrReal, denseDesc%iAtomStart, dqIn)
            #:endif
            else
              call chrgMixerReal%mix(dqInpRed, dqDiffRed)
            #:if WITH_SCALAPACK
              ! Synchronise charges in order to avoid mixers that store a history drifting apart
              call mpifx_bcast(env%mpi%globalComm, dqInpRed)
            #:endif

              call OrbitalEquiv_expand(dqInpRed(:nIneqMixElements), iEqOrbitals, orb, dqIn)

              if (allocated(dftbU) .or. allocated(onsMEs)) then
                dqBlockIn(:,:,:,:) = 0.0_dp
                if (allocated(dftbU)) then
                  call dftbU%expandBlock(dqInpRed, iEqBlockDFTBU, orb, dqBlockIn, species,&
                      & orbEquiv=iEqOrbitals)
                else
                  call Onsblock_expand(dqInpRed, iEqBlockOnSite, orb, dqBlockIn,&
                      & orbEquiv=iEqOrbitals)
                end if
              end if
            end if

          end if

          if (allocated(hybridXc)) then
            call ud2qm(dqIn)
          end if

        end if
      else
        tConverged = .true.
      end if

      if (tConverged) then
        exit lpSCC
      end if

      if (allocated(spinW)) then
        dqPerShell = 0.0_dp
        do iAt = 1, nAtom
          iSp = species(iAt)
          do iSh = 1, orb%nShell(iSp)
            dqPerShell(iSh,iAt,:nSpin) = dqPerShell(iSh,iAt,:nSpin) +&
                & sum(dqIn(orb%posShell(iSh,iSp): orb%posShell(iSh+1,iSp)-1,iAt,:nSpin),dim=1)
          end do
        end do

      end if

    end do lpSCC

    if (tSccCalc .and. .not.tConverged .and. maxSccIter > 1) then
      if (isSccConvRequired) then
        @:RAISE_ERROR(errStatus, -1, "SCC in perturbation is NOT converged, maximal SCC&
            & iterations exceeded")
      else
        quietNan = ieee_value(1.0_dp, ieee_quiet_nan)
        dRho(:,:) = quietNan
        dqOut(:,:,:) = quietNan
        if (allocated(idRho)) idRho(:,:) = quietNan
        if (allocated(dEi)) dEi(:,:,:) = quietNan
        if (allocated(dEf)) dEf(:) = quietNan
        if (allocated(dqBlockOut)) dqBlockOut(:,:,:,:) = quietNan
        if (allocated(dPsiReal)) dPsiReal(:,:,:) = quietNan
        if (allocated(dPsiCmplx)) dPsiCmplx(:,:,:,:) = quietNan
        ! restore back-ups of the input charges/density matrix in case the iteration fails
        dqIn(:,:,:) = dqInBackup
        if (allocated(dqBlockIn)) then
          dqBlockInBackup(:,:,:,:) = dqBlockIn
        end if
        if (allocated(dRhoIn)) then
          dRhoInBackup(:) = dRhoIn
        end if
        call warning("SCC in perturbation is NOT converged, maximal SCC iterations exceeded")
      end if
    end if

    if (present(dDipole)) then
      ! Note, this is origin dependent for periodic geometries:
      dDipole(:) = -matmul(coord(:,:nAtom),sum(dqOut(:,:nAtom,1),dim=1))
    end if

  end subroutine response


  !> Initialise variables for perturbation
  subroutine init_perturbation(parallelKS, tolDegen, nOrbs, nKpts, nSpin, nIndepHam, maxFill,&
      & filling, ham, nFilled, nEmpty, dHam, dRho, idHam, idRho, degenTransform, hybridXc,&
      & sSqrReal, over, neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell,&
      & dRhoOut, dRhoIn, dRhoInSqr, dRhoOutSqr, dPotential, orb, nAtom, tMetallic, neFermi,&
      & eigvals, tempElec, Ef, kWeight)

    !> The k-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Tolerance for degeneracy between eigenvalues
    real(dp), intent(in) :: tolDegen

    !> Number of orbitals
    integer, intent(out) :: nOrbs

    !> Number of k-points
    integer, intent(out) :: nKpts

    !> Number of spin channels
    integer, intent(out) :: nSpin

    !> Number of separate hamiltonians
    integer, intent(out) :: nIndepHam

    !> Maximum occupation of single particle states
    real(dp), intent(out) :: maxFill

    !> Filling of un-perturbed system
    real(dp), intent(in) :: filling(:,:,:)

    !> Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> Number of (partly) filled states in each [nIndepHam,kpt]
    integer, intent(out), allocatable :: nFilled(:,:)

    !> First (partly) empty state in each [nIndepHam,kpt]
    integer, intent(out), allocatable :: nEmpty(:,:)

    !> Derivative of hamiltonian
    real(dp), intent(out), allocatable :: dHam(:,:)

    !> Derivative of sparse density matrix
    real(dp), intent(out), allocatable :: dRho(:,:)

    !> Imaginary part of derivative of hamiltonian
    real(dp), intent(out), allocatable :: idHam(:,:)

    !> Imaginary part of derivative of sparse density matrix
    real(dp), intent(out), allocatable :: idRho(:,:)

    !> For orbital transformations in the  case of degeneracies
    type(TRotateDegen), intent(out), allocatable :: degenTransform(:)

    !> Data for hybrid functional calculation
    class(THybridXcFunc), allocatable, intent(in) :: hybridXc

    !> Square matrix for overlap (if needed in hybrid functional calculation)
    real(dp), allocatable, intent(out) :: sSqrReal(:,:)

    !> Sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Derivative of density matrix
    real(dp), target, allocatable, intent(out) :: dRhoIn(:)

    !> Derivative of density matrix
    real(dp), target, allocatable, intent(out) :: dRhoOut(:)

    !> Delta density matrix for hybrid xc-functional calculations
    real(dp), pointer :: dRhoInSqr(:,:,:)

    !> Delta density matrix for hybrid xc-functional calculations
    real(dp), pointer :: dRhoOutSqr(:,:,:)

    ! derivative of potentials
    type(TPotentials), intent(out) :: dPotential

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> Is there a finite density of states at the Fermi energy
    logical, allocatable, intent(out) :: tMetallic(:,:)

    !> Number of electrons at the Fermi energy (if metallic)
    real(dp), allocatable, intent(inout) :: neFermi(:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    integer :: iS, iK, iLev, ii

    nOrbs = size(filling,dim=1)
    nKpts = size(filling,dim=2)
    nSpin = size(ham,dim=2)

    select case(nSpin)
    case(1,4)
      nIndepHam = 1
    case(2)
      nIndepHam = 2
    end select
    select case(nSpin)
    case(1)
      maxFill = 2.0_dp
    case(2,4)
      maxFill = 1.0_dp
    end select

    allocate(nFilled(nIndepHam, nKpts))
    allocate(nEmpty(nIndepHam, nKpts))
    nFilled(:,:) = -1
    do iS = 1, nIndepHam
      do iK = 1, nKPts
        do iLev = 1, nOrbs
          if ( filling(iLev, iK, iS) < epsilon(1.0) ) then
            ! assumes Fermi filling, so above this is empty
            nFilled(iS, iK) = iLev - 1
            exit
          end if
        end do
        ! check if channel is fully filled
        if (nFilled(iS, iK) < 0) then
          nFilled(iS, iK) = nOrbs
        end if
      end do
    end do
    nEmpty(:,:) = -1
    do iS = 1, nIndepHam
      do iK = 1, nKpts
        do iLev = 1, nOrbs
          if ( abs( filling(iLev, iK, iS) - maxFill ) > epsilon(1.0)) then
            ! assumes Fermi filling, so this is filled
            nEmpty(iS, iK) = iLev
            exit
          end if
        end do
        !> Check is channel is empty
        if (nEmpty(iS, iK) < 0) then
          nEmpty(iS, iK) = 1
        end if
      end do
    end do

    allocate(dHam(size(ham,dim=1),nSpin))
    allocate(dRho(size(ham,dim=1),nSpin))
    if (nSpin == 4) then
      allocate(idHam(size(ham,dim=1),nSpin))
      allocate(idRho(size(ham,dim=1),nSpin))
      idHam(:,:) = 0.0_dp
    end if

    allocate(degenTransform(parallelKS%nLocalKS))
    do ii = 1, size(degenTransform)
      call TRotateDegen_init(degenTransform(ii), tolerance=tolDegen)
    end do

    if (allocated(hybridXc)) then
      allocate(sSqrReal(nOrbs, nOrbs))
      sSqrReal(:,:) = 0.0_dp
    #:if not WITH_SCALAPACK
      call unpackHS(sSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell)
    #:endif
      allocate(dRhoOut(nOrbs * nOrbs * nSpin))
      dRhoOutSqr(1:nOrbs, 1:nOrbs, 1:nSpin) => dRhoOut(:nOrbs*nOrbs*nSpin)
      allocate(dRhoIn(nOrbs * nOrbs * nSpin))
      dRhoInSqr(1:nOrbs, 1:nOrbs, 1:nSpin) => dRhoIn(:nOrbs*nOrbs*nSpin)
    else
      dRhoInSqr => null()
      dRhoOutSqr => null()
    end if

    call TPotentials_init(dPotential, orb, nAtom, nSpin, 0,0 )

    allocate(tMetallic(nIndepHam, nKpts))
    tMetallic(:,:) = .not.(nFilled == nEmpty -1)
    if (any(tMetallic)) then
      write(stdOut,*)'Metallic system'
      ! Density of electrons at the Fermi energy, required to correct later for shift in Fermi level
      ! at q=0 in metals
      if (allocated(neFermi)) then
        deallocate(neFermi)
      end if
      allocate(neFermi(nIndepHam))
      do iS = 1, nIndepHam
        neFermi(iS) = 0.0_dp
        do iK = 1, nKpts
          do ii = nEmpty(iS, iK), nFilled(iS, iK)
            neFermi(iS) = neFermi(iS) + kWeight(iK) * deltamn(Ef(iS), eigvals(ii,iK,iS), tempElec)
          end do
        end do
      end do
      neFermi(:) = maxFill * neFermi
      write(stdOut,*)'Density of states at the Fermi energy Nf (a.u.):', neFermi
    else
      write(stdOut,*)'Non-metallic system'
    end if

  end subroutine init_perturbation


  !> Print the Born effective charges
  subroutine writeBorn(fd, bornCharges, wrtWhichCharges)

    !> File id for data
    integer, intent(in) :: fd

    !> Born effective charges
    real(dp), intent(in) :: bornCharges(:,:,:)

    !> List of which specific atom positions the DFTB charges are to be differentiated
    integer, intent(in), allocatable :: wrtWhichCharges(:)

    integer :: ii, iAt, iCart, nAtom

    nAtom = size(bornCharges, dim=3)

    write(fd, "(A)")'Born effective charges (a.u.)'
    ! i.e. derivative of dipole moment wrt to atom positions, i.e. (d/dx) (q-q0) x, or equivalently
    ! derivative of forces wrt to a homogeneous electric field
    if (allocated(wrtWhichCharges)) then
      do ii = 1, size(wrtWhichCharges)
        iAt = wrtWhichCharges(ii)
        write(fd, "(1X,A,1X,I0)")'Atom', iAt
        do iCart = 1, 3
          write(fd, "(3F16.8)")bornCharges(:, iCart, ii)
        end do
        write(fd, *)
      end do
    else
      do iAt = 1, nAtom
        write(fd, "(1X,A,1X,I0)")'Atom', iAt
        do iCart = 1, 3
          write(fd, "(3F16.8)")bornCharges(:, iCart, iAt)
        end do
        write(fd, *)
      end do
    end if

  end subroutine writeBorn


  !> Count number of perturbation derivatives to perform
  function countDerivs(wrtWhichCharges, nCombinedCharges, wrtCombinedCharges, maxDerivs)&
      & result(nDerivs)

    !> List of which specific external charge positions (MM atoms) the DFTB charges are to be
    !! differentiated (note, incompatible with wrtCombinedCharges)
    integer, intent(in), allocatable :: wrtWhichCharges(:)

    !> Number of combined charges in each group
    integer, intent(in), allocatable :: nCombinedCharges(:)

    !> List of groups of external charge positions (MM atoms) the DFTB charges are to be
    !! differentiated (note, incompatible with wrtWhichCharges)
    integer, intent(in), allocatable :: wrtCombinedCharges(:,:)

    !> Maximum number of possible derivatives
    integer, intent(in) :: maxDerivs

    integer :: nDerivs

    nDerivs = 0
    if (allocated(wrtWhichCharges)) then
      nDerivs = size(wrtWhichCharges)
    else if (allocated(wrtCombinedCharges)) then
      nDerivs = size(wrtCombinedCharges, dim=2)
      @:ASSERT(size(nCombinedCharges) == nDerivs)
      @:ASSERT(maxval(nCombinedCharges) <= size(wrtCombinedCharges, dim=1))
    else
      nDerivs = maxDerivs
    end if
    @:ASSERT(nDerivs <= maxDerivs)

  end function countDerivs

end module dftbp_derivs_perturb
