!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Contains the routines for initialising modes.
module modes_initmodes
  use dftbp_common_accuracy, only : dp
  use dftbp_type_typegeometryhsd, only : TGeometry
  use modes_inputdata, only : TInputData
  use dftbp_common_status, only : TStatus
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_file, only : TFileDescr, setDefaultBinaryAccess, openFile, closeFile,&
      & TOpenOptions
  use dftbp_common_globalenv, only : stdOut
#:if WITH_MAGMA
  use dftbp_common_gpuenv, only : TGpuEnv, TGpuEnv_init
#:endif
  use dftbp_io_message, only : error, warning
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_math_bisect, only : bisection
#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : mpifx_send, mpifx_bcast
#:endif
#:if WITH_SCALAPACK
  use dftbp_dftbplus_initprogram, only : getDenseDescBlacs
  use dftbp_extlibs_scalapackfx, only : CSRC_, RSRC_, MB_, NB_, scalafx_getlocalshape,&
      & scalafx_indxl2g, linecomm
#:endif
  implicit none

  private
  public :: TModesMain
  public :: setEigvecGauge
#:if WITH_MPI
  public :: setEigvecGaugeBlacs
#:endif


  type :: TModesMain

    !> Geometry
    type(TGeometry) :: geo

    !> Atomic masses to build dynamical matrix
    real(dp), allocatable :: atomicMasses(:)

    !> Dynamical matrix
    real(dp), allocatable :: dynMatrix(:,:)

    !> Eigenvalues of the dynamical matrix
    real(dp), allocatable :: eigen(:)

    !> Eigenvectors of the dynamical matrix
    real(dp), allocatable :: eigenModesScaled(:,:)

    !> Displacement vectors for every atom in every mode. Shape: [3, nAtom, nDerivs]
    real(dp), allocatable :: displ(:,:,:)

    !> Born charges matrix
    real(dp), allocatable :: bornMatrix(:)

    !> Derivatives of Born charges matrix with respect to electric field, i.e. polarizability
    !! derivatives with respect to atom locations
    real(dp), allocatable :: bornDerivsMatrix(:)

    !> Produce plots of modes
    logical :: tPlotModes

    !> Modes to produce xyz file for
    integer, allocatable :: modesToPlot(:)

    !> Number of modes being plotted
    integer, allocatable :: nModesToPlot

    !> If animating, number of cycles to show in an animation
    integer :: nCycles

    !> Steps in an animation cycle
    integer :: nSteps

    !> Number of atoms which should be moved
    integer :: nMovedAtom

    !> List of atoms in dynamical matrix
    integer, allocatable :: iMovedAtoms(:)

    !> Number of derivatives
    integer :: nDerivs

    !> Produce eigenvectors of modes, either for plotting or for property changes along mode
    !! directions
    logical :: tEigenVectors

    !> Animate mode  or as vectors
    logical :: tAnimateModes

    !> Remove translation modes
    logical :: tRemoveTranslate

    !> Remove rotation modes
    logical :: tRemoveRotate

    !> Eigensolver choice
    integer :: iSolver

    !> Dense matrix descriptor for H and S
    type(TDenseDescr), allocatable :: denseDesc

  contains

    procedure :: initProgramVariables
    procedure :: getDenseDescCommon
    procedure :: allocateDenseMatrices

  end type TModesMain

  !> Root node name of the input tree
  character(len=*), parameter :: rootTag = "modes"

  !> Input file name
  character(len=*), parameter :: hsdInput = "modes_in.hsd"

  !> Parsed output name
  character(len=*), parameter :: hsdParsedInput = "modes_pin.hsd"

  !> Version of the input document
  integer, parameter :: parserVersion = 3


contains

  !> Initializes the variables in the module based on the parsed input.
  subroutine initProgramVariables(this, input, env)

    !> Instance
    class(TModesMain), intent(inout) :: this

    !> Parsed input data
    type(TInputData), intent(inout) :: input

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

  #:if WITH_SCALAPACK
    !! Error status
    type(TStatus) :: errStatus
  #:endif

    @:ASSERT(input%tInitialized)
    @:ASSERT(allocated(input%ctrl%atomicMasses))
    if (allocated(input%ctrl%hessianFile)) then
      @:ASSERT(.not. allocated(input%ctrl%hessian))
    else
      @:ASSERT(allocated(input%ctrl%hessian))
    end if
    @:ASSERT(allocated(input%ctrl%iMovedAtoms))

    write(stdOut, "(/, A)") "Starting initialization..."
    write(stdOut, "(A80)") repeat("-", 80)

    ! set the same access for readwrite as for write (we do not open any files in readwrite mode)
    call setDefaultBinaryAccess(input%ctrl%binaryAccessTypes(1), input%ctrl%binaryAccessTypes(2),&
        & input%ctrl%binaryAccessTypes(2))

    this%geo = input%geo

    call move_alloc(input%ctrl%atomicMasses, this%atomicMasses)
    call move_alloc(input%ctrl%bornMatrix, this%bornMatrix)
    call move_alloc(input%ctrl%bornDerivsMatrix, this%bornDerivsMatrix)

    if (allocated(input%ctrl%modesToPlot)) then
      this%modesToPlot = input%ctrl%modesToPlot
      this%nModesToPlot = size(this%modesToPlot)
    end if

    this%nCycles = input%ctrl%nCycles
    this%nSteps = input%ctrl%nSteps

    call move_alloc(input%ctrl%iMovedAtoms, this%iMovedAtoms)
    this%nMovedAtom = size(this%iMovedAtoms)
    this%nDerivs = 3 * this%nMovedAtom

    this%tPlotModes = input%ctrl%tPlotModes
    this%tAnimateModes = input%ctrl%tAnimateModes
    this%tRemoveTranslate = input%ctrl%tRemoveTranslate
    this%tRemoveRotate = input%ctrl%tRemoveRotate

    this%iSolver = input%ctrl%iSolver

    this%tEigenVectors = this%tPlotModes .or. allocated(this%bornMatrix)&
        & .or. allocated(this%bornDerivsMatrix)

  #:if WITH_MPI
    call env%initMpi(input%ctrl%parallelOpts%nGroup)
  #:endif

  #:if WITH_SCALAPACK
    associate (blacsOpts => input%ctrl%parallelOpts%blacsOpts)
      call env%initBlacs(blacsOpts%blockSize, blacsOpts%blockSize, this%nDerivs, this%nMovedAtom,&
          & errStatus)
    if (errStatus%hasError()) then
      if (errStatus%code == -1) then
        call warning("Insufficient atoms for this number of MPI processors.")
      end if
      call error(errStatus%message)
    end if
    end associate
  #:endif

    call this%getDenseDescCommon()

    allocate(this%displ(3, this%geo%nAtom, this%nDerivs), source=0.0_dp)
    allocate(this%eigen(this%denseDesc%fullSize))

  #:if WITH_SCALAPACK
    associate (blacsOpts => input%ctrl%parallelOpts%blacsOpts)
      call getDenseDescBlacs(env, blacsOpts%blockSize, blacsOpts%blockSize, this%denseDesc, .false.)
    end associate
  #:endif

    call this%allocateDenseMatrices(env)

  #:if WITH_SCALAPACK
    if (allocated(input%ctrl%hessianFile)) then
      @:ASSERT(.not. allocated(input%ctrl%hessian))
      allocate(input%ctrl%hessian, mold=this%dynMatrix)
      call readHessianDirectBlacs(env, input%ctrl%hessianFile, this%denseDesc, input%ctrl%hessian)
      call dynMatFromHessianBlacs(env, this%denseDesc, this%atomicMasses, input%ctrl%hessian,&
          & this%dynMatrix)
    end if
  #:else
    @:ASSERT(allocated(this%dynMatrix))
    if (allocated(input%ctrl%hessianFile)) then
      @:ASSERT(.not. allocated(input%ctrl%hessian))
      allocate(input%ctrl%hessian(this%nDerivs, this%nDerivs))
      call readHessianDirect(input%ctrl%hessianFile, input%ctrl%hessian)
      call dynMatFromHessian(this%atomicMasses, input%ctrl%hessian, this%dynMatrix)
    else
      @:ASSERT(allocated(input%ctrl%hessian))
      call dynMatFromHessian(this%atomicMasses, input%ctrl%hessian, this%dynMatrix)
    end if
  #:endif

  end subroutine initProgramVariables


  !> Set up storage for dense matrices, either on a single processor, or as BLACS matrices.
  subroutine allocateDenseMatrices(this, env)

    !> Instance
    class(TModesMain), intent(inout) :: this

    !> Computing environment
    type(TEnvironment), intent(in) :: env

    integer :: nLocalCols, nLocalRows

  #:if WITH_SCALAPACK
    call scalafx_getlocalshape(env%blacs%orbitalGrid, this%denseDesc%blacsOrbSqr, nLocalRows,&
        & nLocalCols)
    allocate(this%eigenModesScaled(nLocalRows, nLocalCols))
  #:else
    nLocalRows = this%denseDesc%fullSize
    nLocalCols = this%denseDesc%fullSize
    if (this%tPlotModes) allocate(this%eigenModesScaled(nLocalRows, nLocalCols))
  #:endif

    allocate(this%dynMatrix(nLocalRows, nLocalCols))

  end subroutine allocateDenseMatrices


  !> Generate description of the total large square matrices.
  subroutine getDenseDescCommon(this)

    !> Instance
    class(TModesMain), intent(inout) :: this

    print *, '###########'
    if (allocated(this%denseDesc)) deallocate(this%denseDesc)
    allocate(this%denseDesc)
    allocate(this%denseDesc%iAtomStart(this%nMovedAtom + 1))
    call buildSquaredAtomIndex(this%denseDesc%iAtomStart)

    this%denseDesc%t2Component = .false.
    this%denseDesc%nOrb = this%nDerivs
    this%denseDesc%fullSize = this%nDerivs

  end subroutine getDenseDescCommon


  !> Builds an atom offset array for the dynamical matrix.
  subroutine buildSquaredAtomIndex(iAtomStart)

    !> Offset array for each atom on exit
    integer, intent(out) :: iAtomStart(:)

    integer :: ind, iAt1
    integer :: nAtom

    nAtom = size(iAtomStart) - 1

    ind = 1
    do iAt1 = 1, nAtom
      iAtomStart(iAt1) = ind
      ind = ind + 3
    end do
    iAtomStart(nAtom+1) = ind

  end subroutine buildSquaredAtomIndex


#:if WITH_SCALAPACK
  !> Reads Hessian directly from file (bypasses the slow HSD parser).
  subroutine readHessianDirectBlacs(env, fname, denseDesc, hessian)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> File name of Hessian
    character(len=*), intent(in) :: fname

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Distributed Hessian matrix
    real(dp), intent(out) :: hessian(:,:)

    real(dp), allocatable :: localHessCol(:)
    type(linecomm) :: distributor
    type(TFileDescr) :: fd
    integer :: nDerivs, iErr, iCol

    hessian(:,:) = 0.0_dp

    nDerivs = denseDesc%fullSize
    allocate(localHessCol(nDerivs))

    if (env%mpi%tGlobalLead) then
      call openFile(fd, fname, options=TOpenOptions(form='formatted', action='read'), iostat=iErr)
    end if

    call distributor%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")

    do iCol = 1, nDerivs
      if (env%mpi%tGlobalLead) then
        read(fd%unit, *, iostat=iErr) localHessCol
        if (iErr /= 0) then
          call error("Error during direct reading '" // fname // "'.")
        end if
        call distributor%setline_lead(env%blacs%orbitalGrid, iCol, localHessCol, hessian)
      else
        call distributor%setline_follow(env%blacs%orbitalGrid, iCol, hessian)
      end if
    end do

    if (env%mpi%tGlobalLead) call closeFile(fd)

  end subroutine readHessianDirectBlacs


  !> Initializes the dynamical matrix from the Hessian input.
  subroutine dynMatFromHessianBlacs(env, denseDesc, atomicMasses, hessian, dynMatrix)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Atomic masses to build dynamical matrix
    real(dp), intent(in) :: atomicMasses(:)

    !> Hessian matrix
    real(dp), intent(in) :: hessian(:,:)

    !> Dynamical matrix
    real(dp), intent(out) :: dynMatrix(:,:)

    !> Auxiliary variables
    integer :: iCount, jCount, iDeriv1, iDeriv2, iAt1, iAt2

    @:ASSERT(all(shape(hessian) == shape(dynMatrix)))

    dynMatrix(:,:) = hessian

    ! mass weight the Hessian matrix to get the dynamical matrix
    ! H_{ij} = \frac{\partial^2 \Phi}{\partial u_i \partial u_j}
    ! D_{ij} = \frac{H_{ij}}{\sqrt{m_i m_j}}
    !        = \frac{\partial^2 \Phi}{\partial w_i \partial w_j}
    ! where w_i = \sqrt{m_i} u_i
    do iCount = 1, size(dynMatrix, dim=2)
      iDeriv1 = scalafx_indxl2g(iCount, denseDesc%blacsOrbSqr(NB_), env%blacs%orbitalGrid%mycol,&
          & denseDesc%blacsOrbSqr(CSRC_), env%blacs%orbitalGrid%ncol)
      call bisection(iAt1, denseDesc%iAtomStart, iDeriv1)
      do jCount = 1, size(dynMatrix, dim=1)
        iDeriv2 = scalafx_indxl2g(jCount, denseDesc%blacsOrbSqr(MB_), env%blacs%orbitalGrid%myrow,&
            & denseDesc%blacsOrbSqr(RSRC_), env%blacs%orbitalGrid%nrow)
        call bisection(iAt2, denseDesc%iAtomStart, iDeriv2)
        dynMatrix(jCount, iCount) = dynMatrix(jCount, iCount)&
            & / (sqrt(atomicMasses(iAt1)) * sqrt(atomicMasses(iAt2)))
      end do
    end do

  end subroutine dynMatFromHessianBlacs


  !> Returns gauge-corrected eigenvectors, such that the fist non-zero coefficient of each mode is
  !! positive.
  subroutine setEigvecGaugeBlacs(env, denseDesc, eigvec)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Gauge corrected eigenvectors on exit. Shape: [iCoeff, iMode]
    real(dp), intent(inout) :: eigvec(:,:)

    !! Type for communicating a row or a column of a distributed matrix
    type(linecomm) :: collector

    !! Temporary storage for a single line of the collected, dense, square matrix
    real(dp), allocatable :: localLine(:)

    !! Auxiliary variables
    integer :: iCoeff, iDerivs, nDerivs, iCount, iGlobCol

    integer, allocatable :: prefac(:)

    nDerivs = denseDesc%fullSize
    allocate(localLine(nDerivs))
    allocate(prefac(nDerivs))

    call collector%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")

    if (env%mpi%tGlobalLead) then
      ! runs over eigenvectors
      do iDerivs = 1, nDerivs
        call collector%getline_lead(env%blacs%orbitalGrid, iDerivs, eigvec, localLine)
        ! runs over eigenvector coefficients
        lpCoeff: do iCoeff = 1, size(localLine)
          if (abs(localLine(iCoeff)) > 1.0e-10_dp) then
            prefac(iDerivs) = nint(sign(1.0_dp, localLine(iCoeff)))
            exit lpCoeff
          end if
        end do lpCoeff
      end do
    else
      do iDerivs = 1, nDerivs
        if (env%mpi%tGroupLead) then
          call collector%getline_lead(env%blacs%orbitalGrid, iDerivs, eigvec, localLine)
          call mpifx_send(env%mpi%interGroupComm, localLine, env%mpi%interGroupComm%leadrank)
        else
          call collector%getline_follow(env%blacs%orbitalGrid, iDerivs, eigvec)
        end if
      end do
    end if

    call mpifx_bcast(env%mpi%globalComm, prefac)

    ! set gauge
    do iCount = 1, size(eigvec, dim=2)
      iGlobCol = scalafx_indxl2g(iCount, denseDesc%blacsOrbSqr(NB_), env%blacs%orbitalGrid%mycol,&
          & denseDesc%blacsOrbSqr(CSRC_), env%blacs%orbitalGrid%ncol)
      eigvec(:, iCount) = eigvec(:, iCount) * real(prefac(iGlobCol), dp)
    end do

  end subroutine setEigvecGaugeBlacs

#:endif

  !> Reads Hessian directly from file (bypasses the slow HSD parser).
  subroutine readHessianDirect(fname, hessian)

    !> File name of Hessian
    character(len=*), intent(in) :: fname

    !> Hessian matrix
    real(dp), intent(out) :: hessian(:,:)

    !! File descriptor
    type(TFileDescr) :: fd

    !! Error status
    integer :: iErr

    hessian(:,:) = 0.0_dp

    call openFile(fd, fname, options=TOpenOptions(form='formatted', action='read'), iostat=iErr)
    if (iErr /= 0) then
      call error("Could not open file '" // fname // "' for direct reading." )
    end if
    read(fd%unit, *, iostat=iErr) hessian
    if (iErr /= 0) then
      call error("Error during direct reading '" // fname // "'.")
    end if
    call closeFile(fd)

  end subroutine readHessianDirect


  !> Initializes the dynamical matrix from the Hessian input.
  subroutine dynMatFromHessian(atomicMasses, hessian, dynMatrix)

    !> Atomic masses to build dynamical matrix
    real(dp), intent(in) :: atomicMasses(:)

    !> Hessian matrix
    real(dp), intent(in) :: hessian(:,:)

    !> Dynamical matrix
    real(dp), intent(out) :: dynMatrix(:,:)

    !> Auxiliary variables
    integer :: ii, jj, kk, ll, iCount, jCount

    @:ASSERT(size(hessian, dim=1) == size(hessian, dim=2))
    @:ASSERT(3 * size(atomicMasses) == size(hessian, dim=1))

    dynMatrix(:,:) = hessian

    ! mass weight the Hessian matrix to get the dynamical matrix
    ! H_{ij} = \frac{\partial^2 \Phi}{\partial u_i \partial u_j}
    ! D_{ij} = \frac{H_{ij}}{\sqrt{m_i m_j}}
    !        = \frac{\partial^2 \Phi}{\partial w_i \partial w_j}
    ! where w_i = \sqrt{m_i} u_i
    iCount = 0
    do ii = 1, size(atomicMasses)
      do kk = 1, 3
        iCount = iCount + 1
        jCount = 0
        do jj = 1, size(atomicMasses)
          do ll = 1, 3
            jCount = jCount + 1
            dynMatrix(jCount, iCount) = dynMatrix(jCount, iCount)&
                & / (sqrt(atomicMasses(ii)) * sqrt(atomicMasses(jj)))
          end do
        end do
      end do
    end do

  end subroutine dynMatFromHessian


  !> Returns gauge-corrected eigenvectors, such that the first non-zero coefficient of each mode is
  !! positive.
  subroutine setEigvecGauge(eigvec)

    !> Gauge corrected eigenvectors on exit. Shape: [iCoeff, iMode]
    real(dp), intent(inout) :: eigvec(:,:)

    !! Auxiliary variables
    integer :: iMode, iCoeff

    do iMode = 1, size(eigvec, dim=2)
      lpCoeff: do iCoeff = 1, size(eigvec, dim=1)
        if (abs(eigvec(iCoeff, iMode)) > 1e2_dp * epsilon(1.0_dp)) then
          if (sign(1.0_dp, eigvec(iCoeff, iMode)) < 0.0_dp) then
            eigvec(:, iMode) = -eigvec(:, iMode)
          end if
          exit lpCoeff
        end if
      end do lpCoeff
    end do

  end subroutine setEigvecGauge

end module modes_initmodes
