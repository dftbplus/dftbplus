!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Module for linear response derivative calculations using perturbation methods at q=0 for fixed
!> structures
module dftbp_derivs_perturb
  use dftbp_common_accuracy, only : dp, mc
  use dftbp_common_constants, only : Hartree__eV, quaternionName
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_file, only : TFile
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_status, only : TStatus
  use dftbp_derivs_fermihelper, only : theta, deltamn, invDiff
  use dftbp_derivs_linearresponse, only : dRhoReal, dRhoFermiChangeReal,&
      & dRhoCmplx, dRhoFermiChangeCmplx, dRhoPauli, dRhoFermiChangePauli
  use dftbp_derivs_rotatedegen, only : TRotateDegen, TRotateDegen_init
  use dftbp_dftb_blockpothelper, only : appendBlockReduced
  use dftbp_dftb_dftbplusu, only : TDftbU, TDftbU_init, plusUFunctionals
  use dftbp_dftbplus_mainio, only : writeDerivBandOut
  use dftbp_dftbplus_outputfiles, only : derivVBandOut
  use dftbp_dftb_onsitecorrection, only : addOnsShift, onsblock_expand
  use dftbp_dftb_orbitalequiv, only : OrbitalEquiv_reduce, OrbitalEquiv_expand
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_populations, only : mulliken, denseMulliken, getChargePerShell, getOnsitePopulation
  use dftbp_dftb_potentials, only : TPotentials, TPotentials_init
  use dftbp_dftb_rangeseparated, only : TRangeSepFunc
  use dftbp_dftb_scc, only : TScc
  use dftbp_dftb_shift, only : addShift, totalShift
  use dftbp_dftb_spin, only : getSpinShift, ud2qm, qm2ud
  use dftbp_dftb_thirdorder, only : TThirdOrder,  TThirdOrderInp, ThirdOrder_init
  use dftbp_io_commonformats, only : format2U
  use dftbp_io_message, only : warning
  use dftbp_io_taggedoutput, only : TTaggedWriter, tagLabels
  use dftbp_mixer_mixer, only : TMixer, mix, reset
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_parallelks, only : TParallelKS, TParallelKS_init
#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : mpifx_allreduceip, MPI_SUM
#:endif
#:if WITH_SCALAPACK
  use dftbp_extlibs_scalapackfx, only : DLEN_, scalafx_getdescriptor
#:else
  use dftbp_dftb_sparse2dense, only : unpackHS
#:endif
  implicit none

  private
  public :: TResponse, TResponse_init, responseSolverTypes


  !> Namespace for possible perturbation solver methods
  type :: TPerturbSolverTypesEnum

    ! coupled perturbed solved with sum over unperturbed spectra
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

  end type TResponse


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
      & over, orb, nAtom, species, neighbourList, nNeighbourSK, denseDesc, iSparseStart,&
      & img2CentCell, coord, sccCalc, maxSccIter, sccTol, isSccConvRequired, nMixElements,&
      & nIneqMixElements, iEqOrbitals, tempElec, Ef, spinW, thirdOrd, dftbU, iEqBlockDftbu,&
      & onsMEs, iEqBlockOnSite, rangeSep, nNeighbourLC, pChrgMixer, kPoint, kWeight, iCellVec,&
      & cellVec, polarisability, dEi, dqOut, neFermi, dEfdE, errStatus, omega)

    !> Instance
    class(TResponse), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Filling
    real(dp), intent(in) :: filling(:,:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in), allocatable :: eigVecsReal(:,:,:)

    !> ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> chemical species
    integer, intent(in) :: species(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> SCC module internal variables
    type(TScc), intent(inout), allocatable :: sccCalc

    !> maximal number of SCC iterations
    integer, intent(in) :: maxSccIter

    !> Tolerance for SCC convergence
    real(dp), intent(in) :: sccTol

    !> Use converged derivatives of charges
    logical, intent(in) :: isSccConvRequired

    !> nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
    !> blocks (if a DFTB+U or onsite correction calculation)
    integer, intent(in) :: nMixElements

    !> nr. of inequivalent charges
    integer, intent(in) :: nIneqMixElements

    !> Equivalence relations between orbitals
    integer, intent(in), allocatable :: iEqOrbitals(:,:,:)

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> Third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> equivalence mapping for dual charge blocks
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Charge mixing object
    type(TMixer), intent(inout), allocatable :: pChrgMixer

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Electric polarisability
    real(dp), intent(out) :: polarisability(:, :, :)

    !> Derivatives of eigenvalues, if required
    real(dp), allocatable, intent(inout) :: dEi(:,:,:,:)

    !> Derivatives of Mulliken charges, if required
    real(dp), allocatable, intent(inout) :: dqOut(:,:,:,:)

    !> Number of electrons at the Fermi energy (if metallic)
    real(dp), allocatable, intent(inout) :: neFermi(:)

    !> Derivative of the Fermi energy (if metallic)
    real(dp), allocatable, intent(inout) :: dEfdE(:,:)

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
    real(dp) :: dqIn(orb%mOrb,nAtom,size(ham, dim=2))

    real(dp), allocatable :: dqBlockIn(:,:,:,:), SSqrReal(:,:)
    real(dp), allocatable :: dqBlockOut(:,:,:,:)

    ! derivative of potentials
    type(TPotentials) :: dPotential

    logical :: tSccCalc
    logical, allocatable :: tMetallic(:,:)

    real(dp), allocatable :: dPsiReal(:,:,:)
    complex(dp), allocatable :: dPsiCmplx(:,:,:,:)

    ! variables used for range separated contributions, note this stays in the up/down
    ! representation throughout if spin polarised
    real(dp), pointer :: dRhoOutSqr(:,:,:), dRhoInSqr(:,:,:)
    real(dp), allocatable, target :: dRhoOut(:), dRhoIn(:)

    !> For transformation in the  case of degeneracies
    type(TRotateDegen), allocatable :: transform(:)

    write(stdOut,*)
    write(stdOut,*)'Perturbation calculation of electric polarisability'
    write(stdOut,*)

    call init_perturbation(parallelKS, this%tolDegen, nOrbs, nKpts, nSpin, nIndepHam, maxFill,&
        & filling, ham, nFilled, nEmpty, dHam, dRho, idHam, idRho, transform, rangesep, sSqrReal,&
        & over, neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell, dRhoOut,&
        & dRhoIn, dRhoInSqr, dRhoOutSqr, dPotential, orb, nAtom, tMetallic, neFermi, eigvals,&
        & tempElec, Ef, kWeight)

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

    ! polarisation direction
    ! note: could MPI parallelise over this in principle
    lpCart: do iCart = 1, 3

      write(stdOut,*)"Polarisabilty for field along ", trim(quaternionName(iCart+1))

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
        dPotential%extAtom(iAt,1) = coord(iCart,iAt)
      end do
      call totalShift(dPotential%extShell, dPotential%extAtom, orb, species)
      call totalShift(dPotential%extBlock, dPotential%extShell, orb, species)

      do iOmega = 1, size(omega)

        if (allocated(dEfdETmp)) then
          dEfdETmp(:) = 0.0_dp
        end if

        dEiTmp(:,:,:) = 0.0_dp
        call response(env, parallelKS, dPotential, nAtom, orb, species, neighbourList,&
            & nNeighbourSK, img2CentCell, iSparseStart, denseDesc, over, iEqOrbitals, sccCalc,&
            & sccTol, isSccConvRequired, maxSccIter, pChrgMixer, nMixElements, nIneqMixElements,&
            & dqIn, dqOut(:,:,:,iCart), rangeSep, nNeighbourLC, sSqrReal, dRhoInSqr, dRhoOutSqr,&
            & dRhoIn, dRhoOut, nSpin, maxFill, spinW, thirdOrd, dftbU, iEqBlockDftbu, onsMEs,&
            & iEqBlockOnSite, dqBlockIn, dqBlockOut, eigVals, transform, dEiTmp, dEfdETmp, Ef,&
            & this%isEfFixed, dHam, idHam, dRho, idRho, tempElec, tMetallic, neFermi, nFilled,&
            & nEmpty, kPoint, kWeight, cellVec, iCellVec, eigVecsReal, eigVecsCplx, dPsiReal,&
            & dPsiCmplx, coord, errStatus, omega(iOmega), dDipole=polarisability(:, iCart, iOmega),&
            & eta=this%eta)
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

    write(stdOut,*)
    write(stdOut,*)'Polarisability (a.u.)'
    do iOmega = 1, size(omega)
      write(stdOut,*)
      if (abs(omega(iOmega)) > epsilon(0.0_dp)) then
        write(stdOut, format2U)"Polarisability at omega = ", omega(iOmega), ' H ',&
            & omega(iOmega) * Hartree__eV, ' eV'
      else
        write(stdOut, "(A)")"Static polarisability:"
      end if
      do iCart = 1, 3
        write(stdOut,"(3E20.12)")polarisability(:, iCart, iOmega)
      end do
    end do
    write(stdOut,*)

  end subroutine wrtEField


  !> Response with respect to a potential at atomic sites
  subroutine wrtVAtom(this, env, parallelKS, isAutotestWritten, autotestTagFile,&
      & isTagResultsWritten, resultsTagFile, taggedWriter, isBandWritten, fdDetailedOut,&
      & filling, eigvals, eigVecsReal, eigVecsCplx, ham, over, orb, nAtom, species, neighbourList,&
      & nNeighbourSK, denseDesc, iSparseStart, img2CentCell, isRespKernelRPA, sccCalc, maxSccIter,&
      & sccTol, isSccConvRequired, nMixElements, nIneqMixElements, iEqOrbitals, tempElec, Ef,&
      & spinW, thirdOrd, dftbU, iEqBlockDftbu, onsMEs, iEqBlockOnSite, rangeSep, nNeighbourLC,&
      & pChrgMixer, kPoint, kWeight, iCellVec, cellVec, nEFermi, errStatus, omega, isHelical, coord)

    !> Instance
    class(TResponse), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> should regression test data be written
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
    type(TFile), allocatable, intent(in) :: fdDetailedOut

    !> Filling
    real(dp), intent(in) :: filling(:,:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in), allocatable :: eigVecsReal(:,:,:)

    !> ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> chemical species
    integer, intent(in) :: species(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> SCC module internal variables
    type(TScc), intent(inout), allocatable :: sccCalc

    !> Should the kernel be evaluated at the RPA level (non-SCC) or self-consistent
    logical, intent(in) :: isRespKernelRPA

    !> maximal number of SCC iterations
    integer, intent(in) :: maxSccIter

    !> Tolerance for SCC convergence
    real(dp), intent(in) :: sccTol

    !> Use converged derivatives of charges
    logical, intent(in) :: isSccConvRequired

    !> nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
    !> blocks (if a DFTB+U or onsite correction calculation)
    integer, intent(in) :: nMixElements

    !> nr. of inequivalent charges
    integer, intent(in) :: nIneqMixElements

    !> Equivalence relations between orbitals
    integer, intent(in), allocatable :: iEqOrbitals(:,:,:)

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> Third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> equivalence mapping for dual charge blocks
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Charge mixing object
    type(TMixer), intent(inout), allocatable :: pChrgMixer

    !> k-points
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

    ! used for range separated contributions, note this stays in the up/down representation
    ! throughout if spin polarised
    real(dp), pointer :: dRhoOutSqr(:,:,:), dRhoInSqr(:,:,:)
    real(dp), allocatable, target :: dRhoOut(:), dRhoIn(:)

    !> For transformation in the  case of degeneracies
    type(TRotateDegen), allocatable :: transform(:)

    real(dp), allocatable :: dEf(:), dqOut(:,:,:), dqOutTmp(:,:,:,:,:)

    integer :: nIter, fd, iOmega
    character(mc) :: atLabel

    logical :: isSccRequired

    if (isRespKernelRPA) then
      nIter = 1
      isSccRequired = .false.
    else
      nIter = maxSccIter
      isSccRequired = isSccConvRequired
    end if

    call init_perturbation(parallelKS, this%tolDegen, nOrbs, nKpts, nSpin, nIndepHam, maxFill,&
        & filling, ham, nFilled, nEmpty, dHam, dRho, idHam, idRho, transform, rangesep, sSqrReal,&
        & over, neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell, dRhoOut,&
        & dRhoIn, dRhoInSqr, dRhoOutSqr, dPotential, orb, nAtom, tMetallic, neFermi, eigvals,&
        & tempElec, Ef, kWeight)

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
            & sccTol, isSccRequired, nIter, pChrgMixer, nMixElements, nIneqMixElements, dqIn,&
            & dqOut, rangeSep, nNeighbourLC, sSqrReal, dRhoInSqr, dRhoOutSqr, dRhoIn, dRhoOut,&
            & nSpin, maxFill, spinW, thirdOrd, dftbU, iEqBlockDftbu, onsMEs, iEqBlockOnSite,&
            & dqBlockIn, dqBlockOut, eigVals, transform, dEi, dEf, Ef, this%isEfFixed, dHam, idHam,&
            & dRho, idRho, tempElec, tMetallic, neFermi, nFilled, nEmpty, kPoint, kWeight, cellVec,&
            & iCellVec, eigVecsReal, eigVecsCplx, dPsiReal, dPsiCmplx, coord, errStatus,&
            & omega(iOmega), isHelical=isHelical, eta=this%eta)
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

        if (allocated(fdDetailedOut)) then
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
      open(newunit=fd, file=autotestTagFile, action="write", status="old", position="append")
      call taggedWriter%write(fd, tagLabels%dqdV, dqOutTmp)
      call taggedWriter%write(fd, tagLabels%dqnetdV, dqNetAtomTmp)
      close(fd)
    end if
    if (isTagResultsWritten) then
      open(newunit=fd, file=resultsTagFile, action="write", status="old", position="append")
      call taggedWriter%write(fd, tagLabels%dEigenDV, dEiTmp)
      call taggedWriter%write(fd, tagLabels%dqdV, dqOutTmp)
      call taggedWriter%write(fd, tagLabels%dqnetdV, dqNetAtomTmp)
      close(fd)
    end if

  end subroutine wrtVAtom


  ! Internal routines

  !> Evaluates response, given the external perturbation
  subroutine response(env, parallelKS, dPotential, nAtom, orb, species, neighbourList,&
      & nNeighbourSK, img2CentCell, iSparseStart, denseDesc, over, iEqOrbitals, sccCalc, sccTol,&
      & isSccConvRequired, maxSccIter, pChrgMixer, nMixElements, nIneqMixElements, dqIn, dqOut,&
      & rangeSep, nNeighbourLC, sSqrReal, dRhoInSqr, dRhoOutSqr, dRhoIn, dRhoOut, nSpin, maxFill,&
      & spinW, thirdOrd, dftbU, iEqBlockDftbu, onsMEs, iEqBlockOnSite, dqBlockIn, dqBlockOut,&
      & eigVals, transform, dEi, dEf, Ef, isEfFixed, dHam, idHam,  dRho, idRho, tempElec,&
      & tMetallic, neFermi, nFilled, nEmpty, kPoint, kWeight, cellVec, iCellVec, eigVecsReal,&
      & eigVecsCplx, dPsiReal, dPsiCmplx, coord, errStatus, omega, dDipole, isHelical, eta)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    ! derivative of potentials
    type(TPotentials), intent(inout) :: dPotential

    !> Charge mixing object
    type(TMixer), intent(inout), allocatable :: pChrgMixer

    !> nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
    !> blocks (if a DFTB+U or onsite correction calculation)
    integer, intent(in) :: nMixElements

    !> nr. of inequivalent charges
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

    !> maximal number of SCC iterations
    integer, intent(in) :: maxSccIter

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> chemical species
    integer, intent(in) :: species(:)

    !> spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> Third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Is there a finite density of states at the Fermi energy
    logical, intent(in) :: tMetallic(:,:)

    !> Number of electrons at the Fermi energy (if metallic)
    real(dp), allocatable, intent(inout) :: neFermi(:)

    !> Tolerance for SCC convergence
    real(dp), intent(in) :: sccTol

    !> Use converged derivatives of charges
    logical, intent(in) :: isSccConvRequired

    !> Maximum allowed number of electrons in a single particle state
    real(dp), intent(in) :: maxFill

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Equivalence relations between orbitals
    integer, intent(in), allocatable :: iEqOrbitals(:,:,:)

    !> For transformation in the  case of degeneracies
    type(TRotateDegen), allocatable, intent(inout) :: transform(:)

    !> Derivative of density matrix
    real(dp), target, allocatable, intent(inout) :: dRhoIn(:)

    !> Derivative of density matrix
    real(dp), target, allocatable, intent(inout) :: dRhoOut(:)

    !> delta density matrix for rangeseparated calculations
    real(dp), pointer :: dRhoOutSqr(:,:,:), dRhoInSqr(:,:,:)

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> equivalence mapping for dual charge blocks
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

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Is a fixed Fermi energy used
    logical, intent(in) :: isEfFixed

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in), allocatable :: eigVecsReal(:,:,:)

    !> ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Derivative of Fermi energy
    real(dp), intent(inout), allocatable :: dEf(:)

    !> Derivative of block charges (input)
    real(dp), allocatable, intent(inout) :: dqBlockIn(:,:,:,:)

    !> Derivative of block charges (output)
    real(dp), allocatable, intent(inout) :: dqBlockOut(:,:,:,:)

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Square matrix for overlap (if needed)
    real(dp), allocatable, intent(inout) :: sSqrReal(:,:)

    !> k-points
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

  #:if WITH_SCALAPACK
    ! need distributed matrix descriptors
    integer :: desc(DLEN_), nn

    nn = denseDesc%fullSize
    call scalafx_getdescriptor(env%blacs%orbitalGrid, nn, nn, env%blacs%rowBlockSize,&
        & env%blacs%columnBlockSize, desc)
  #:endif

    tSccCalc = allocated(sccCalc)

    @:ASSERT(abs(omega) <= epsilon(0.0_dp) .or. present(eta))

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
      if (allocated(rangeSep)) then
        dRhoIn(:) = 0.0_dp
        dRhoOut(:) = 0.0_dp
      end if
    end if

    if (present(dDipole)) then
      dDipole(:) = 0.0_dp
    end if

    if (tSccCalc) then
      call reset(pChrgMixer, size(dqInpRed))
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
        call sccCalc%getShiftPerAtom(dPotential%intAtom(:,1))
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
      call addShift(dHam, over, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
          & iSparseStart, nAtom, img2CentCell, dPotential%intBlock)

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

          call dRhoReal(env, dHam, neighbourList, nNeighbourSK, iSparseStart, img2CentCell,&
              & denseDesc, iKS, parallelKS, nFilled, nEmpty, eigVecsReal, eigVals, Ef, tempElec,&
              & orb, drho(:,iS), dRhoOutSqr, rangeSep, over, nNeighbourLC, transform(iKS),&
              & species,&
            #:if WITH_SCALAPACK
              & desc,&
            #:endif
              & dEi, dPsiReal, coord, errStatus, omega, isHelical, eta=eta)
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
              & transform(iKS), species, coord,&
            #:if WITH_SCALAPACK
              & desc,&
            #:endif
              & dEi, dPsiCmplx, errStatus, omega, isHelical, eta=eta)
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

          call dRhoCmplx(env, dHam, neighbourList, nNeighbourSK, iSparseStart,&
              & img2CentCell, denseDesc, parallelKS, nFilled, nEmpty, eigvecsCplx, eigVals, Ef,&
              & tempElec, orb, dRho, kPoint, kWeight, iCellVec, cellVec, iKS, transform(iKS),&
              & species, coord,&
            #:if WITH_SCALAPACK
              & desc,&
            #:endif
              & dEi, dPsiCmplx, errStatus, omega, isHelical, eta=eta)
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
          call mulliken(dqOut(:,:,iS), over, drho(:,iS), orb, &
              & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)

          dEf(iS) = -sum(dqOut(:, :, iS)) / neFermi(iS)

          if (abs(dEf(iS)) > 10.0_dp*epsilon(1.0_dp)) then
            ! Fermi level changes, so need to correct for the change in the number of charges

            if (allocated(eigVecsReal)) then

              ! real case, no k-points
              call dRhoFermiChangeReal(dRhoExtra(:, iS), env, maxFill, parallelKS, iKS,&
                  & neighbourList, nNeighbourSK, img2CentCell, iSparseStart, dEf, Ef, nFilled,&
                  & nEmpty, eigVecsReal, orb, denseDesc, tempElec, eigVals, dRhoOutSqr, species,&
                  & coord,&
                #:if WITH_SCALAPACK
                  & desc,&
                #:endif
                  & isHelical)

            elseif (nSpin > 2) then

              ! two component wavefunction cases
              call dRhoFermiChangePauli(dRhoExtra, idRhoExtra, env, parallelKS, iKS,&
                  & kPoint, kWeight, iCellVec, cellVec, neighbourList, nNEighbourSK,&
                  & img2CentCell, iSparseStart, dEf, Ef, nFilled, nEmpty, eigVecsCplx, orb,&
                  & denseDesc, tempElec, eigVals, species, coord,&
                #:if WITH_SCALAPACK
                  & desc,&
                #:endif
                  & errStatus, isHelical)
              @:PROPAGATE_ERROR(errStatus)

            else

              ! Complex case with k-points
              call dRhoFermiChangeCmplx(dRhoExtra, env, maxFill, parallelKS, iKS, kPoint, kWeight,&
                  & iCellVec, cellVec, neighbourList, nNEighbourSK, img2CentCell, iSparseStart,&
                  & dEf, Ef, nFilled, nEmpty, eigVecsCplx, orb, denseDesc, tempElec, eigVals,&
                  & species, coord,&
                #:if WITH_SCALAPACK
                  & desc,&
                #:endif
                  & isHelical)

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
        call mulliken(dqOut(:,:,iS), over, drho(:,iS), orb, &
            & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
        if (allocated(dftbU) .or. allocated(onsMEs)) then
          dqBlockOut(:,:,:,iS) = 0.0_dp
          call mulliken(dqBlockOut(:,:,:,iS), over, drho(:,iS), orb, neighbourList%iNeighbour,&
              & nNeighbourSK, img2CentCell, iSparseStart)
        end if
      end do

      if (tSccCalc) then

        if (allocated(rangeSep)) then
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
            if (allocated(rangeSep)) then
              dRhoIn(:) = dRhoOut
              call denseMulliken(dRhoInSqr, SSqrReal, denseDesc%iAtomStart, dqIn)
            else
              dqIn(:,:,:) = dqOut
              dqInpRed(:) = dqOutRed
              if (allocated(dftbU) .or. allocated(onsMEs)) then
                dqBlockIn(:,:,:,:) = dqBlockOut
              end if
            end if

          else

            if (allocated(rangeSep)) then
              call mix(pChrgMixer, dRhoIn, dqDiffRed)
              call denseMulliken(dRhoInSqr, SSqrReal, denseDesc%iAtomStart, dqIn)
            else
              call mix(pChrgMixer, dqInpRed, dqDiffRed)
            #:if WITH_MPI
              ! Synchronise charges in order to avoid mixers that store a history drifting apart
              call mpifx_allreduceip(env%mpi%globalComm, dqInpRed, MPI_SUM)
              dqInpRed(:) = dqInpRed / env%mpi%globalComm%size
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

          if (allocated(rangeSep)) then
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
        call warning("SCC in perturbation is NOT converged, maximal SCC iterations exceeded")
      end if
    end if

    if (present(dDipole)) then
      dDipole(:) = -matmul(coord(:,:nAtom),sum(dqOut(:,:nAtom,1),dim=1))
    end if

  end subroutine response


  !> Initialise variables for perturbation
  subroutine init_perturbation(parallelKS, tolDegen, nOrbs, nKpts, nSpin, nIndepHam, maxFill,&
      & filling, ham, nFilled, nEmpty, dHam, dRho, idHam, idRho, transform, rangesep, sSqrReal,&
      & over, neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell, dRhoOut, dRhoIn,&
      & dRhoInSqr, dRhoOutSqr, dPotential, orb, nAtom, tMetallic, neFermi, eigvals, tempElec, Ef,&
      & kWeight)

    !> K-points and spins to process
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
    type(TRotateDegen), intent(out), allocatable :: transform(:)

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(in) :: rangeSep

    !> Square matrix for overlap (if needed in range separated calculation)
    real(dp), allocatable, intent(out) :: sSqrReal(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Derivative of density matrix
    real(dp), target, allocatable, intent(out) :: dRhoIn(:)

    !> Derivative of density matrix
    real(dp), target, allocatable, intent(out) :: dRhoOut(:)

    !> delta density matrix for rangeseparated calculations
    real(dp), pointer :: dRhoInSqr(:,:,:)

    !> delta density matrix for rangeseparated calculations
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

    allocate(transform(parallelKS%nLocalKS))
    do ii = 1, size(transform)
      call TRotateDegen_init(transform(ii), tolDegen)
    end do

    if (allocated(rangeSep)) then
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

end module dftbp_derivs_perturb
