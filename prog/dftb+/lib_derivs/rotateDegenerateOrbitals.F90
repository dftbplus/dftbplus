!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module containing routines to make linear combinations of orbitals for degenerate perturbation
!> from a hermitian/symmetric matrix
module dftbp_rotateDegenerateOrbs
  use dftbp_accuracy, only : dp
  use dftbp_eigensolver, only : heev
  use dftbp_qm
  use dftbp_message
  use dftbp_wrappedintr
  implicit none
  private

  public :: TDegeneracyTransform

#:set FLAVOURS = [('cmplx', 'complex', 'Cmplx'), ('real', 'real', 'Real')]

  !> Data type for resolving degenerate orbitals in perturbation expressions
  type :: TDegeneracyTransform
    private
#:for NAME, TYPE, LABEL in FLAVOURS

    !> Dense unitary matrix to transform orbitals
    ${TYPE}$(dp), allocatable :: ${LABEL}$U(:,:)

    !> Individual sub-blocks of unitary matrices to transform orbitals, if much smaller than U
    type(wrapped${LABEL}$2), allocatable :: ${LABEL}$UBlock(:)

#:endfor

    !> Block ranges as needed in UBlock, deallocated again by time of return if dense U case
    integer, allocatable :: blockRange(:,:)

    !> To which group of states any particular one belongs
    integer, allocatable :: degenerateGroup(:)

    !> numerical tolerance for deciding degeneracy
    real(dp) :: tolerance

    !> Sub-range of states if needed
    integer :: eiRange(2)

    !> Denominator for fraction of states in a degenerate group before the using dense transforms
    integer :: minDegenerateFraction

    !> Minimum number of states in a degenerate group before the using dense transforms
    integer :: minDegenerateStates

    !> Should the order of states (derivatives) be reversed compared to the eigensolver return
    logical :: tReverseOrder = .false.

    ! redundant variables (could test from allocation status and size), but makes code simpler to
    ! read in a few places:

    !> Number of groups of degenerate orbitals
    integer :: nGrp = -1

    !> Stores data structure for a real unitary transform
    logical :: tReal = .false.

    !> Stores data structure for a complex unitary transform
    logical :: tComplx = .false.

    !> Is any storage allocated
    logical :: tAllocateStorage = .false.

    !> Order of matrix
    integer :: nOrb = -1

    logical :: tFullUMatrix = .false.

  contains

    !> Initialises instance and set some optional parameters
    procedure :: init

    !> Generate unitary matrix for degenerate perturbation
    procedure :: generateCmplxUnitary
    procedure :: generateRealUnitary
    generic :: generateUnitary => generateCmplxUnitary, generateRealUnitary

    !> Transform a matrix to be suitable for degenerate perturbation
    procedure :: degenerateCmplxTransform
    procedure :: degenerateRealTransform
    generic :: degenerateTransform => degenerateCmplxTransform, degenerateRealTransform

    !> Applys a unitary transformation to a matrix
    procedure :: applyCmplxUnitary
    procedure :: applyRealUnitary
    generic :: applyUnitary => applyCmplxUnitary, applyRealUnitary

    !> Are a pair of states in the same degenerate group
    procedure :: degenerate

    !> Release memory and cleans up
    procedure :: destroy

  end type TDegeneracyTransform


#:if WITH_SCALAPACK
#:else
  !> Fraction of total matrix at which to use full U instead of blocked form
  integer, parameter :: maxBlockFraction = 4
#:endif

contains

  !> Initialise the structure
  subroutine init(self, smallestBlock, smallestFraction, tolerance, eiRange)

    !> Instance
    class(TDegeneracyTransform), intent(inout) :: self

    !> Smallest number of degenerate orbitals at which the dense algorithm should be used
    integer, intent(in), optional :: smallestBlock

    !> Smallest fraction of the matrix at which the dense algorithm should be used
    integer, intent(in), optional :: smallestFraction

    !> Tolerance for degeneracy testing
    real(dp), intent(in), optional :: tolerance

    !> Sub-range of states if needed, for example for metallic finite temperature in parallel gauge
    integer, intent(in), optional :: eiRange(2)

    if (present(smallestBlock)) then
      self%minDegenerateStates = smallestBlock
    else
      self%minDegenerateStates = 500
    end if

    if (present(smallestFraction)) then
      self%minDegenerateFraction = smallestFraction
    else
      self%minDegenerateFraction = 4 ! a quarter of the matrix
    end if

    if (present(tolerance)) then
      self%tolerance = tolerance
    else
      ! a few times eps, just in case of minor symmetry breaking
      self%tolerance = 128.0_dp*epsilon(0.0_dp)
    end if

    if (present(eiRange)) then
      self%eiRange(:) = eiRange
    else
      self%eiRange(:) = -1
    end if

  end subroutine init


#:if WITH_SCALAPACK

  ! to do

#:else

#:for NAME, TYPE, LABEL in FLAVOURS

  !> Set up unitary transformation of matrix for degenerate states, producing combinations that are
  !> orthogonal under the action of the matrix. This is the ${TYPE}$ case.
  subroutine generate${LABEL}$Unitary(self, matrixToProcess, ei)

    !> Instance
    class(TDegeneracyTransform), intent(inout) :: self

    !> Matrix elements between (potentially degenerate) orbitals
    ${TYPE}$(dp), intent(in) :: matrixToProcess(:,:)

    !> Eigenvalues of local block of degenerate matrix
    real(dp), intent(in) :: ei(:)

    integer :: ii, nGrp, iGrp, maxRange, nInBlock, iS, iE
    integer, allocatable :: degeneracies(:)
    ${TYPE}$(dp), allocatable :: subBlock(:,:)
    real(dp), allocatable :: eigenvals(:)
    logical :: tFullU
    integer :: eiRange(2)

  #:if TYPE == 'real'
    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: zero = 0.0_dp
  #:else
    real(dp), parameter :: one = (1.0_dp, 0.0_dp)
    real(dp), parameter :: zero = (0.0_dp, 0.0_dp)
  #:endif

    self%nOrb = size(ei)
    if (all(self%eiRange == [-1,-1] )) then
      eiRange(:) = [1, self%nOrb]
    else
      eiRange(:) = self%eiRange
    end if

    if (allocated(self%blockRange)) then
      if (any(shape(self%blockRange) /= [2, self%nOrb])) then
        deallocate(self%blockRange)
      end if
    end if
    if (.not.allocated(self%blockRange)) then
      allocate(self%blockRange(2, self%nOrb))
    end if
    self%blockRange(:,:) = 0
    if (allocated(self%degenerateGroup)) then
      if (size(self%degenerateGroup) /= self%nOrb) then
        deallocate(self%degenerateGroup)
      end if
    end if
    if (.not.allocated(self%degenerateGroup)) then
      allocate(self%degenerateGroup(self%nOrb))
    end if

    call degeneracyRanges(self%blockRange, self%nGrp, Ei, self%tolerance, eiRange,&
        & self%degenerateGroup)

    write(*,*)'groups',self%degenerateGroup

    maxRange = maxval(self%blockRange(2,:self%nGrp) - self%blockRange(1,:self%nGrp)) + 1
    if (maxRange == 1) then
      ! no transformations required
      ! also nGrp == nOrb
      return
    end if

    ! decide if the full matrix or sub-blocks are to be used
    self%tFullUMatrix = maxRange < self%nOrb / self%minDegenerateFraction&
        & .or. maxRange > self%minDegenerateStates

    ! memory set-up if needed and set to unit matrix for transformation
    if (self%tFullUMatrix) then

      ! make whole matrix as U

      if (allocated(self%${LABEL}$U)) then
        if (any(shape(self%${LABEL}$U) /= [self%nOrb, self%nOrb])) then
          deallocate(self%${LABEL}$U)
        end if
      end if
      if (.not.allocated(self%${LABEL}$U)) then
        allocate(self%${LABEL}$U(self%nOrb, self%nOrb))
      end if
      self%${LABEL}$U(:,:) = zero
      do ii = 1, self%nOrb
        self%${LABEL}$U(ii,ii) = one
      end do

    else

      ! just make individual blocks, as much smaller and faster to use

      if (allocated(self%${LABEL}$UBlock)) then
        do iGrp = 1, size(self%${LABEL}$UBlock)
          deallocate(self%${LABEL}$UBlock(iGrp)%data)
        end do
        deallocate(self%${LABEL}$UBlock)
      end if
      allocate(self%${LABEL}$UBlock(self%nGrp))
      do iGrp = 1, self%nGrp
        nInBlock = self%blockRange(2,iGrp) - self%blockRange(1,iGrp) + 1
        allocate(self%${LABEL}$UBlock(iGrp)%data(nInBlock,nInBlock))
        if (nInBlock > 1) then
          self%${LABEL}$UBlock(iGrp)%data = zero
        else
          self%${LABEL}$UBlock(iGrp)%data = one
        end if
      end do
    end if

    allocate(subBlock(maxRange,maxRange))
    allocate(eigenvals(maxRange))

    ! Form the unitary elements
    do iGrp = 1, self%nGrp
      subBlock(:,:) = zero
      eigenvals(:) = 0.0_dp
      iS = self%blockRange(1,iGrp)
      iE = self%blockRange(2,iGrp)
      nInBlock = iE - iS + 1
      if (nInBlock == 1) then
        cycle
      end if
      subBlock(:nInBlock, :nInBlock) = matrixToProcess(iS:iE, iS:iE)
      call heev(subBlock(:nInBlock, :nInBlock), eigenvals, 'L', 'V')
      if (self%tFullUMatrix) then
        if (self%tReverseOrder) then
          self%${LABEL}$U(iS:iE, iS:iE) = subBlock(:nInBlock, nInBlock:1:-1)
        else
          self%${LABEL}$U(iS:iE, iS:iE) = subBlock(:nInBlock, :nInBlock)
        end if
      else
        if (self%tReverseOrder) then
          self%${LABEL}$UBlock(iGrp)%data = subBlock(:nInBlock, nInBlock:1:-1)
        else
          self%${LABEL}$UBlock(iGrp)%data = subBlock(:nInBlock, :nInBlock)
        end if
      end if
    end do

  end subroutine generate${LABEL}$Unitary


  !> ${TYPE}$ case of orthogonalising a hermitian/symmetric matrix against degenerate perturbation
  !> operations by applying a (stored) unitary transform
  subroutine degenerate${LABEL}$Transform(self, matrixToProcess)

    !> Instance
    class(TDegeneracyTransform), intent(in) :: self

    !> Matrix elements between (potentially degenerate) orbitals
    ${TYPE}$(dp), intent(inout) :: matrixToProcess(:,:)

    integer :: iGrp, nGrp, iS, iE, ii, jj

    if (self%tFullUMatrix) then

      call makeSimiliarityTrans(matrixToProcess, self%${LABEL}$U, 'R')

    else if (allocated(self%${LABEL}$UBlock)) then

      do iGrp = 1, self%nGrp

        iS = self%blockRange(1,iGrp)
        iE = self%blockRange(2,iGrp)
        if (iS == iE) then
          cycle
        end if

        matrixToProcess(:,iS:iE) = matmul(matrixToProcess(:,iS:iE), self%${LABEL}$UBlock(iGrp)%data)

      end do

      do iGrp = 1, self%nGrp

        iS = self%blockRange(1,iGrp)
        iE = self%blockRange(2,iGrp)
        if (iS == iE) then
          matrixToProcess(iS,iS) = cmplx(real(matrixToProcess(iS,iS),dp), 0.0_dp, dp)
          cycle
        end if

      #:if TYPE == 'real'
        matrixToProcess(iS:iE,:) = matmul(transpose(self%RealUBlock(iGrp)%data),&
            & matrixToProcess(iS:iE,:))
      #:else
        matrixToProcess(iS:iE,:) = matmul(transpose(conjg(self%CmplxUBlock(iGrp)%data)),&
            & matrixToProcess(iS:iE,:))
      #:endif

        do ii = iS, iE
          matrixToProcess(ii,ii) = cmplx(real(matrixToProcess(ii,ii),dp), 0.0_dp, dp)
          do jj = ii + 1, iE
            matrixToProcess(jj,ii) = 0.0_dp
            matrixToProcess(ii,jj) = 0.0_dp
          end do
        end do

      end do

    end if

    ! clean up to exact symmetry
  #:if TYPE == 'real'
    matrixToProcess(:,:) = 0.5_dp * (matrixToProcess + transpose(matrixToProcess))
  #:else
    matrixToProcess(:,:) = 0.5_dp * (matrixToProcess + transpose(conjg(matrixToProcess)))
  #:endif

  end subroutine degenerate${LABEL}$Transform


  !> ${TYPE}$ case of transforming a unitary matrix against degenerate perturbation operations
  subroutine apply${LABEL}$Unitary(self, matrixToProcess)

    !> Instance
    class(TDegeneracyTransform), intent(in) :: self

    !> Matrix elements
    ${TYPE}$(dp), intent(inout) :: matrixToProcess(:,:)

    integer :: iGrp, nGrp, iS, iE

    if (self%tFullUMatrix) then

      #:if TYPE == 'real'
        matrixToProcess(:,:) = matmul(matrixToProcess, self%RealU)
      #:else
        matrixToProcess(:,:) = matmul(matrixToProcess, self%CmplxU)
      #:endif

    else if (allocated(self%${LABEL}$UBlock)) then

      do iGrp = 1, self%nGrp

        iS = self%blockRange(1,iGrp)
        iE = self%blockRange(2,iGrp)
        if (iS == iE) then
          cycle
        end if

      #:if TYPE == 'real'
        matrixToProcess(:, iS:iE) = matmul(matrixToProcess(:, iS:iE), self%RealUBlock(iGrp)%data)
      #:else
        matrixToProcess(:, iS:iE) = matmul(matrixToProcess(:, iS:iE), self%CmplxUBlock(iGrp)%data)
      #:endif

      end do

    end if

  end subroutine apply${LABEL}$Unitary

#:endfor

  subroutine destroy(self)

    !> Instance
    class(TDegeneracyTransform), intent(inout) :: self

    integer :: iGrp

    if (allocated(self%blockRange)) then
      deallocate(self%blockRange)
    end if
    if (allocated(self%degenerateGroup)) then
      deallocate(self%degenerateGroup)
    end if

    if (self%tFullUMatrix) then

    #:for _, _, LABEL in FLAVOURS
      if (allocated(self%${LABEL}$U)) then
        deallocate(self%${LABEL}$U)
      end if
    #:endfor

    else

    #:for _, _, LABEL in FLAVOURS
      if (allocated(self%${LABEL}$UBlock)) then
        do iGrp = 1, size(self%${LABEL}$UBlock)
          deallocate(self%${LABEL}$UBlock(iGrp)%data)
        end do
        deallocate(self%${LABEL}$UBlock)
      end if
   #:endfor

    end if

  end subroutine destroy

#:endif

  !> Returns whether states are in the same degenerate group
  pure function degenerate(self, ii, jj)

    !> Instance
    class(TDegeneracyTransform), intent(in) :: self

    !> First state
    integer, intent(in) :: ii

    !> second state
    integer, intent(in) :: jj

    !> Resulting test
    logical :: degenerate

    degenerate = (self%degenerateGroup(ii) == self%degenerateGroup(jj))

  end function degenerate


  !> Find which groups of eigenvales are degenerate to within a tolerance
  !> Note, similar process is used in Casida excited state calculations, so should spin off as its
  !> own module at some point
  subroutine degeneracyRanges(blockRange, nGrp, Ei, tol, eiRange, grpMembership)

    !> Index array for lower and upper states in degenerate group
    integer, intent(out) :: blockRange(:,:)

    !> Number of degenerate groups
    integer, intent(out) :: nGrp

    !> Eigenvalues for degeneracy testing
    real(dp), intent(in) :: Ei(:)

    !> Tolerance for degeneracy testing
    real(dp), intent(in), optional :: tol

    !> sub range of eigenvalues to process
    integer, intent(in), optional :: eiRange(2)

    integer, intent(out), optional :: grpMembership(:)

    integer :: ii, jj, nOrb, iS, iE
    real(dp) :: localTol

    if (present(tol)) then
      localTol = tol
    else
      ! a few times eps, just in case of minor symmetry breaking
      localTol = 128.0_dp * epsilon(0.0_dp)
    end if
    nOrb = size(ei)
    blockRange(:,:) = 0
    nGrp = 0

    if (present(eiRange)) then
      ! set states before group as not of interest
      iS = eiRange(1)
      iE = eiRange(2)
      if (iS == iE) then
        call error("Degeneracy range is itself degenerate, should not be here")
      end if
      do ii = 1, iS - 1
        blockRange(:, ii) = ii
      end do
    else
      iS = 1
      iE = nOrb
    end if
    nGrp = nGrp + iS - 1

    do ii = 1, nGrp
      grpMembership(ii) = ii
    end do

    ii = iS
    do while (ii <= iE)
      nGrp = nGrp + 1
      blockRange(1, nGrp) = ii
      grpMembership(ii) = nGrp
      do jj = ii + 1, iE
        ! assume sorted:
        if ( abs(ei(jj) - ei(jj-1)) > localTol) then
          exit
        end if
        grpMembership(jj) = nGrp
      end do
      ii = jj
      blockRange(2, nGrp) = jj - 1
    end do

    if (present(eiRange)) then
      ! set states after group as not of interest
      do ii = iE + 1, nOrb
        nGrp = nGrp + 1
        blockRange(:, nGrp) = ii
        grpMembership(ii) = nGrp
      end do
    end if

  end subroutine degeneracyRanges

end module dftbp_rotateDegenerateOrbs
