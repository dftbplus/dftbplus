!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains an DIIS mixer
!! The DIIS mixing is done by building a weighted combination over the previous input charges to
!! minimise the residue of the error.
!! Only a specified number of previous charge vectors are considered.
!! The modification based on from Kovalenko et al. (J. Comput. Chem., 20: 928-936 1999) and Patrick
!! Briddon to add a contribution from the gradient vector as well is also used.
!! In order to use the mixer you have to create and reset it.
module dftbp_mixer_diismixer
  use dftbp_common_accuracy, only : dp
  use dftbp_math_lapackroutines, only : gesv
  implicit none

#:set FLAVOURS = [('cmplx', 'complex', 'Cmplx'), ('real', 'real', 'Real')]

  private
#:for NAME, TYPE, LABEL in FLAVOURS
  public :: TDiisMixer${LABEL}$, TDiisMixer${LABEL}$_init, TDiisMixer${LABEL}$_reset
  public :: TDiisMixer${LABEL}$_mix
#:endfor


#:for NAME, TYPE, LABEL in FLAVOURS
  !> Contains the necessary data for an DIIS mixer.
  type TDiisMixer${LABEL}$
    private

    !> Initial mixing parameter
    real(dp) :: initMixParam

    !> Max. nr. of stored prev. vectors
    integer :: mPrevVector

    !> Nr. of stored previous vectors
    integer :: iPrevVector

    !> Nr. of elements in the vectors
    integer :: nElem

    !> Index for the storage
    integer :: indx

    !> Stored previous input charges
    ${TYPE}$(dp), allocatable :: prevQInput(:,:)

    !> Stored prev. charge differences
    ${TYPE}$(dp), allocatable :: prevQDiff(:,:)

    !> True if DIIS used from iteration 2 as well as mixing
    logical :: tFromStart

    !> force modification for gDIIS?
    logical :: tAddIntrpGradient

    !> Alpha factor to add in new information
    real(dp) :: alpha

    !> Holds DIIS mixed gradients from older iterations for downhill direction
    ${TYPE}$(dp), allocatable :: deltaR(:)

  end type TDiisMixer${LABEL}$
#:endfor


contains

#:for NAME, TYPE, LABEL in FLAVOURS
  !> Creates a DIIS mixer instance.
  subroutine TDiisMixer${LABEL}$_init(this, nGeneration, initMixParam, tFromStart, alpha)

    !> Pointer to an initialized DIIS mixer on exit
    type(TDiisMixer${LABEL}$), intent(out) :: this

    !> Nr. of generations (including actual) to consider
    integer, intent(in) :: nGeneration

    !> Mixing parameter for the first nGeneration cycles
    real(dp), intent(in) :: initMixParam

    !> True if using DIIS from iteration 2 as well as mixing
    logical, intent(in), optional :: tFromStart

    !> If present, fraction of extrapolated downhill direction to include in DIIS space
    real(dp), intent(in), optional :: alpha

    @:ASSERT(nGeneration >= 2)

    this%nElem = 0
    this%mPrevVector = nGeneration

    allocate(this%prevQInput(this%nElem, this%mPrevVector))
    allocate(this%prevQDiff(this%nElem, this%mPrevVector))

    this%initMixParam = initMixParam

    if (present(tFromStart)) then
      this%tFromStart = tFromStart
    else
      this%tFromStart = .false.
    end if

    if (present(alpha)) then
      this%tAddIntrpGradient = .true.
      this%alpha = alpha
      allocate(this%deltaR(this%nElem))
    else
      this%tAddIntrpGradient = .false.
      this%alpha = 0.0_dp
      allocate(this%deltaR(0))
    end if

    this%deltaR(:) = 0.0_dp

  end subroutine TDiisMixer${LABEL}$_init


  !> Makes the mixer ready for a new SCC cycle.
  subroutine TDiisMixer${LABEL}$_reset(this, nElem)

    !> DIIS mixer instance
    type(TDiisMixer${LABEL}$), intent(inout) :: this

    !> Nr. of elements in the vectors to mix
    integer, intent(in) :: nElem

    @:ASSERT(nElem > 0)

    if (nElem /= this%nElem) then
      this%nElem = nElem
      deallocate(this%prevQInput)
      deallocate(this%prevQDiff)
      allocate(this%prevQInput(this%nElem, this%mPrevVector))
      allocate(this%prevQDiff(this%nElem, this%mPrevVector))
      if (this%tAddIntrpGradient) then
        deallocate(this%deltaR)
        allocate(this%deltaR(this%nElem))
        this%deltaR(:) = 0.0_dp
      end if
    end if
    this%iPrevVector = 0
    this%indx = 0

  end subroutine TDiisMixer${LABEL}$_reset


  !> Mixes charges according to the DIIS method.
  !!
  !! Warning: The complex-valued DIIS mixer requires flattened hermitian matrices as input.
  !!   You are free to permute the individual elements of the flattened arrays as long as the same
  !!   permutation is applied to qInpResult and qDiff.
  !!   The restriction arises from the assumption that the dot-products of density matrices are
  !!   real-valued (imaginary parts add up to zero due to the hermitian property) and the linear
  !!   system of equations remains real-valued.
  subroutine TDiisMixer${LABEL}$_mix(this, qInpResult, qDiff)

    !> Pointer to the diis mixer
    type(TDiisMixer${LABEL}$), intent(inout) :: this

    !> Input charges on entry, mixed charges on exit.
    ${TYPE}$(dp), intent(inout) :: qInpResult(:)

    !> Charge difference vector between output and input charges
    ${TYPE}$(dp), intent(in) :: qDiff(:)

    real(dp), allocatable :: aa(:,:), bb(:,:)
    integer :: ii, jj

    @:ASSERT(size(qInpResult) == this%nElem)
    @:ASSERT(size(qDiff) == this%nElem)

    if (this%iPrevVector < this%mPrevVector) then
      this%iPrevVector = this%iPrevVector + 1
    end if

    call storeVectors_${NAME}$(this%prevQInput, this%prevQDiff, this%indx, qInpResult, qDiff,&
        & this%mPrevVector)

    if (this%tFromStart .or. this%iPrevVector == this%mPrevVector) then

      if (this%tAddIntrpGradient) then
        ! old DIIS estimate for downhill direction points towards current downhill direction as well
        ! as the actual vector, based on P. Briddon comments
        if (abs(dot_product(this%deltaR, qDiff)) > 0.0_dp) then
          ! mix in larger amounts of the gradient in future
          this%alpha = 1.5_dp * this%alpha
        else
          ! points the other way, mix in less
          this%alpha = 0.5 * this%alpha
        end if
      end if

      allocate(aa(this%iPrevVector + 1, this%iPrevVector + 1))
      allocate(bb(this%iPrevVector + 1, 1))

      aa(:,:) = 0.0_dp
      bb(:,:) = 0.0_dp

      ! (due to the hermitian property of our density matrices, the dot-product below is real)
      do ii = 1, this%iPrevVector
        do jj = 1, this%iPrevVector
          aa(ii, jj) = dot_product(this%prevQDiff(:, ii), this%prevQDiff(:, jj))
        end do
      end do
      aa(this%iPrevVector + 1, 1:this%iPrevVector) = -1.0_dp
      aa(1:this%iPrevVector, this%iPrevVector + 1) = -1.0_dp

      bb(this%iPrevVector + 1, 1) = -1.0_dp

      ! Solve DIIS system of linear equations
      call gesv(aa, bb)

      qInpResult(:) = 0.0_dp
      do ii = 1, this%iPrevVector
        qInpResult(:) = qInpResult + bb(ii, 1) * (this%prevQInput(:, ii) + this%prevQDiff(:, ii))
      end do

      if (this%tAddIntrpGradient) then
        ! add a fraction down the DIIS estimated gradient onto the new solution
        this%deltaR(:) = 0.0_dp
        do ii = 1, this%iPrevVector
          this%deltaR(:) = this%deltaR + bb(ii, 1) * this%prevQDiff(:, ii)
        end do
        qInpResult(:) = qInpResult - this%alpha * this%deltaR
      end if

    end if

    if (this%iPrevVector < this%mPrevVector) then
      ! First few iterations return simple mixed vector
      qInpResult(:) = qInpResult + this%initMixParam * qDiff
    end if

  end subroutine TDiisMixer${LABEL}$_mix


  !> Stores a vector pair in a limited storage.
  !! If the stack is full, oldest vector pair is overwritten.
  subroutine storeVectors_${NAME}$(prevQInp, prevQDiff, indx, qInput, qDiff, mPrevVector)

    !> Contains previous vectors of the first type
    ${TYPE}$(dp), intent(inout) :: prevQInp(:,:)

    !> Contains previous vectors of the second type
    ${TYPE}$(dp), intent(inout) :: prevQDiff(:,:)

    !> Indexing of data
    integer, intent(inout) :: indx

    !> New first vector
    ${TYPE}$(dp), intent(in) :: qInput(:)

    !> New second vector
    ${TYPE}$(dp), intent(in) :: qDiff(:)

    !> Size of the stacks.
    integer, intent(in) :: mPrevVector

    indx = mod(indx, mPrevVector) + 1
    prevQInp(:, indx) = qInput
    prevQDiff(:, indx) = qDiff

  end subroutine storeVectors_${NAME}$
#:endfor

end module dftbp_mixer_diismixer
