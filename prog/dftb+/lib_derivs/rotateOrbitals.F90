!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module containing routines to make linear combinations of orbitals for degenerate perturbation
!> from a hermitian/symmetric matrix
module dftbp_degeneratePerturb
  use dftbp_accuracy, only : dp
  use dftbp_eigensolver, only : heev
  use dftbp_qm
  use dftbp_message
  use dftbp_wrappedintr
  implicit none
  private

  public :: degenerateUnitary, degenerateTransform

#:set FLAVOURS = [('cmplx', 'complex', 'Cmplx'), ('real', 'real', 'Real')]

#:if WITH_SCALAPACK

  ! to do, but will match interfaces for serial

#:else

  !> Generate unitary matrix for degenerate perturbation
  interface degenerateUnitary
  #:for NAME, _, _ in FLAVOURS
    module procedure ${NAME}$DegenerateUnitary
  #:endfor
  end interface degenerateUnitary

  !> Transform a matrix to be suitable for degenerate perturbation
  interface degenerateTransform
  #:for NAME, _ in FLAVOURS
    module procedure ${NAME}$Transform
  #:endfor
  end interface degenerateTransform

#:endif

#:if WITH_SCALAPACK
#:else
  !> Fraction of total matrix at which to use full U instead of blocked form
  integer, parameter :: maxBlockFraction = 4
#:endif

contains

#:if WITH_SCALAPACK

  ! to do

#:else

#:for NAME, TYPE, LABEL in FLAVOURS
  !> ${NAME}$ case of orthogonalising against degenerate perturbation operation cases
  subroutine ${NAME}$Transform(matrixToProcess, ei, U, UBlock, blockRange, tol, eiRange)

    !> Matrix elements between (potentially degenerate) orbitals
    ${TYPE}$(dp), intent(inout) :: matrixToProcess(:,:)

    !> Eigenvalues of states
    real(dp), intent(in) :: ei(:)

    !> Unitary matrix to transform orbitals
    ${TYPE}$(dp), allocatable, intent(inout) :: U(:,:)

    !> Individual sub-blocks of unitary matrices to transform orbitals, if much smaller than U
    type(wrapped${LABEL}$2), allocatable, intent(inout) :: UBlock(:)

    !> Block ranges as needed in UBlock, deallocated again by time of return if dense U case
    integer, allocatable, intent(inout) :: blockRange(:,:)

    !> sub range of eigenvalues to process, for example relevant for parallel gauge in metallic
    !> systems where only partially filled states are dangerous wrt to degeneration
    integer, intent(in), optional :: eiRange(2)

    !> Optional tolerance for comparisions between states
    real(dp), intent(in), optional :: tol

    integer :: iGrp, nGrp, iS, iE

    call degenerateUnitary(U, UBlock, blockRange, matrixToProcess, ei, tol, eiRange)

    if (allocated(U)) then
      call makeSimiliarityTrans(matrixToProcess, U, 'R')
    else if (allocated(UBlock)) then
      nGrp = size(UBlock)
      do iGrp = 1, nGrp
        iS = blockRange(1,iGrp)
        iE = blockRange(2,iGrp)
        if (iS == iE) then
          cycle
        end if
        matrixToProcess(:,iS:iE) = matmul(matrixToProcess(:,iS:iE), UBlock(iGrp)%data)
      end do
      do iGrp = 1, nGrp
        iS = blockRange(1,iGrp)
        iE = blockRange(2,iGrp)
        if (iS == iE) then
          cycle
        end if
      #:if TYPE == 'real'
        matrixToProcess(iS:iE,:) = matmul(transpose(UBlock(iGrp)%data), matrixToProcess(iS:iE,:))
      #:else
        matrixToProcess(iS:iE,:) = matmul(transpose(conjg(UBlock(iGrp)%data)),&
            & matrixToProcess(iS:iE,:))
      #:endif
      end do
    end if

  end subroutine ${NAME}$Transform


  !> Transform matrix for degenerate states into combinations that are orthogonal under the action
  !> of the matrix
  subroutine ${NAME}$DegenerateUnitary(U, UBlock, blockRange, matrixToProcess, ei, tol, eiRange)

    !> Unitary matrix to transform orbitals
    ${TYPE}$(dp), intent(inout), allocatable :: U(:,:)

    !> Individual sub-blocks of unitary matrices to transform orbitals
    type(wrapped${LABEL}$2), intent(inout), allocatable :: UBlock(:)

    !> Block ranges as needed in UBlock, deallocated again by time of return if dense U case
    integer, intent(inout), allocatable :: blockRange(:,:)

    !> Matrix elements between (potentially degenerate) orbitals
    ${TYPE}$(dp), intent(in) :: matrixToProcess(:,:)

    !> Eigenvalues of local block of degenerate matrix
    real(dp), intent(in) :: ei(:)

    !> Tolerance for comparisions between states
    real(dp), intent(in), optional :: tol

    !> sub range of eigenvalues to process, for example relevant for parallel gauge in metallic
    !> systems where only partially filled states are dangerous wrt to degeneration
    integer, intent(in), optional :: eiRange(2)

    integer :: ii, nOrb, nGrp, iGrp, maxRange, nInBlock, iS, iE
    integer, allocatable :: degeneracies(:)
    ${TYPE}$(dp), allocatable :: subBlock(:,:)
    real(dp), allocatable :: eigenvals(:)
    logical :: tFullU

  #:if TYPE == 'real'
    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: zero = 0.0_dp
  #:else
    real(dp), parameter :: one = (1.0_dp, 0.0_dp)
    real(dp), parameter :: zero = (0.0_dp, 0.0_dp)
  #:endif

    nOrb = size(ei)

    if (allocated(blockRange)) then
      if (any(shape(blockRange) /= [2,nOrb])) then
        deallocate(blockRange)
      end if
    end if
    if (.not.allocated(blockRange)) then
      allocate(blockRange(2,nOrb))
    end if
    blockRange(:,:) = 0

    call degeneracyRanges(blockRange, nGrp, Ei, tol, eiRange)

    maxRange = maxval(blockRange(2,:) - blockRange(1,:)) + 1

    if (maxRange == 1) then
      ! no transformations required
      return
    end if

    ! decide if the full matrix or sub-blocks to be used
    tFullU = maxRange < nOrb / maxBlockFraction

    if (tFullU) then
      ! make whole matrix as U
      tFullU = .true.
      if (allocated(U)) then
        if (any(shape(U) /= [nOrb,nOrb])) then
          deallocate(U)
        end if
      end if
      if (.not.allocated(U)) then
        allocate(U(nOrb, nOrb))
      end if
      U(:,:) = zero
      do ii = 1, nOrb
        U(ii,ii) = one
      end do
    else
      ! just make individual blocks, as much smaller and faster to use
      tFullU = .false.
      if (allocated(UBlock)) then
        do iGrp = 1, size(UBlock)
          deallocate(UBlock(iGrp)%data)
        end do
        deallocate(UBlock)
      end if
      allocate(UBlock(nGrp))
      do iGrp = 1, nGrp
        nInBlock = blockRange(2,iGrp) - blockRange(1,iGrp) + 1
        allocate(UBlock(iGrp)%data(nInBlock,nInBlock))
        if (nInBlock > 1) then
          UBlock(iGrp)%data = zero
        else
          UBlock(iGrp)%data = one
        end if
      end do
    end if

    allocate(subBlock(maxRange,maxRange))
    allocate(eigenvals(maxRange))

    do iGrp = 1, nGrp
      subBlock(:,:) = zero
      eigenvals(:) = 0.0_dp
      iS = blockRange(1,iGrp)
      iE = blockRange(2,iGrp)
      nInBlock = iE - iS + 1
      if (nInBlock == 1) then
        cycle
      end if
      subBlock(:nInBlock, :nInBlock) = matrixToProcess(iS:iE, iS:iE)
      call heev(subBlock(:nInBlock, :nInBlock), eigenvals, 'L', 'V')
      if (tFullU) then
        U(iS:iE, iS:iE) = subBlock(:nInBlock, :nInBlock)
      else
        UBlock(iGrp)%data = subBlock(:nInBlock, :nInBlock)
      end if
    end do

  end subroutine ${NAME}$DegenerateUnitary
#:endfor

#:endif

  !> Find which groups of eigenvales are degenerate to within a tolerance
  subroutine degeneracyRanges(blockRange, nGrp, Ei, tol, eiRange)

    !> Index array for lower and upper states in degenerate group
    integer, intent(out) :: blockRange(:,:)

    !> Number of degenerate groups
    integer, intent(out) :: nGrp

    !> Eigenvalues for degeneracy testing
    real(dp), intent(in) :: Ei(:)

    !> Tolerance of testing
    real(dp), intent(in), optional :: tol

    !> sub range of eigenvalues to process
    integer, intent(in), optional :: eiRange(2)

    integer :: ii, jj, nOrb, iS, iE
    real(dp) :: localTol

    if (present(tol)) then
      localTol = tol
    else
      localTol = epsilon(0.0_dp)
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

    ii = iS
    do while (ii <= iE)
      nGrp = nGrp + 1
      blockRange(1, nGrp) = ii
      do jj = ii + 1, iE
        ! assume sorted:
        if ( abs(ei(jj) - ei(jj-1)) > localTol) then
          exit
        end if
      end do
      ii = jj
      blockRange(2, nGrp) = jj -1
    end do

    if (present(eiRange)) then
      ! set states after group as not of interest
      do ii = iE + 1, nOrb
        nGrp = nGrp + 1
        blockRange(:, nGrp) = ii
      end do
    end if

  end subroutine degeneracyRanges

end module dftbp_degeneratePerturb
