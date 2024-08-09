!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the routines for initialising modes.
module modes_initmodes
  use dftbp_common_accuracy, only : dp
  use dftbp_type_typegeometryhsd, only : TGeometry
  use modes_inputdata, only : TInputData
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_file, only : setDefaultBinaryAccess
  implicit none

  private
  public :: TModesMain
  public :: setEigvecGauge

  type :: TModesMain

    !> Geometry
    type(TGeometry) :: geo

    !> Atomic masses to build dynamical matrix
    real(dp), allocatable :: atomicMasses(:)

    !> Dynamical matrix
    real(dp), allocatable :: dynMatrix(:,:)

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

    !> Number of atoms which should be moved.
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

  contains

    procedure :: initProgramVariables

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

    @:ASSERT(input%tInitialized)
    @:ASSERT(allocated(input%ctrl%atomicMasses))
    @:ASSERT(allocated(input%ctrl%hessian))
    @:ASSERT(allocated(input%ctrl%iMovedAtoms))

    write(stdOut, "(/, A)") "Starting initialization..."
    write(stdOut, "(A80)") repeat("-", 80)

    ! set the same access for readwrite as for write (we do not open any files in readwrite mode)
    call setDefaultBinaryAccess(input%ctrl%binaryAccessTypes(1), input%ctrl%binaryAccessTypes(2),&
        & input%ctrl%binaryAccessTypes(2))

    this%geo = input%geo

    call move_alloc(input%ctrl%atomicMasses, this%atomicMasses)

    call dynMatFromHessian(this%atomicMasses, input%ctrl%hessian, this%dynMatrix)

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

    this%tEigenVectors = this%tPlotModes .or. allocated(this%bornMatrix)&
        & .or. allocated(this%bornDerivsMatrix)

  end subroutine initProgramVariables


  !> Initializes the dynamical matrix from the Hessian input.
  subroutine dynMatFromHessian(atomicMasses, hessian, dynMatrix)

    !> Atomic masses to build dynamical matrix
    real(dp), intent(in) :: atomicMasses(:)

    !> Hessian matrix
    real(dp), intent(in) :: hessian(:,:)

    !> Dynamical matrix
    real(dp), intent(out), allocatable :: dynMatrix(:,:)

    !> Auxiliary variables
    integer :: ii, jj, kk, ll, iCount, jCount

    @:ASSERT(size(hessian, dim=1) == size(hessian, dim=2))
    @:ASSERT(3 * size(atomicMasses) == size(hessian, dim=1))

    allocate(dynMatrix, source=hessian)

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
