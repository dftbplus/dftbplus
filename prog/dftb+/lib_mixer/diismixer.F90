!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains an DIIS mixer
!> The DIIS mixing is done by building a weighted combination over the previous input charges to
!> minimise the residue of the error.  Only a specified number of previous charge vectors are
!> considered.
!> The modification based on from Kovalenko et al. (J. Comput. Chem., 20: 928-936 1999) and Patrick
!> Briddon to add a contribution from the gradient vector as well is also used
!> In order to use the mixer you have to create and reset it.
module dftbp_diismixer
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_lapackroutines, only : gesv
  implicit none

  private


  !> Contains the necessary data for an DIIS mixer
  type Tdiismixer
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
    real(dp), allocatable :: prevQInput(:,:)

    !> Stored prev. charge differences
    real(dp), allocatable :: prevQDiff(:,:)

    !> True if DIIS used from iteration 2 as well as mixing
    logical :: tFromStart

    !> force modification for gDIIS?
    logical :: tAddIntrpGradient

    !> Alpha factor to add in new information
    real(dp) :: alpha

    !> Holds DIIS mixed gradients from older iterations for downhill direction
    real(dp), allocatable :: deltaR(:)
  end type Tdiismixer


  !> Creates an DIISMixer instance
  interface init
    module procedure DIISMixer_init
  end interface init


  !> Resets the mixer
  interface reset
    module procedure DIISMixer_reset
  end interface reset


  !> Does the mixing
  interface mix
    module procedure DIISMixer_mix
  end interface mix

  public :: Tdiismixer
  public :: init, reset, mix

contains


  !> Creates an DIIS mixer instance.
  subroutine DIISMixer_init(this, nGeneration, initMixParam,tFromStart,alpha)

    !> Pointer to an initialized DIIS mixer on exit
    type(Tdiismixer), intent(out) :: this

    !> Nr. of generations (including actual) to consider
    integer, intent(in) :: nGeneration

    !> Mixing parameter for the first nGeneration cycles
    real(dp), intent(in) :: initMixParam

    !> True if using DIIS from iteration 2 as well as mixing
    logical, intent(in), optional :: tFromStart

    !> if present, fraction of extrapolated downhill direction to include in DIIS space
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

    this%deltaR = 0.0_dp

  end subroutine DIISMixer_init


  !> Makes the mixer ready for a new SCC cycle
  subroutine DIISMixer_reset(this, nElem)

    !> DIIS mixer instance
    type(Tdiismixer), intent(inout) :: this

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
        this%deltaR = 0.0_dp
      end if
    end if
    this%iPrevVector = 0
    this%indx = 0

  end subroutine DIISMixer_reset


  !> Mixes charges according to the DIIS method
  subroutine DIISMixer_mix(this, qInpResult, qDiff)

    !> Pointer to the diis mixer
    type(Tdiismixer), intent(inout) :: this

    !> Input charges on entry, mixed charges on exit.
    real(dp), intent(inout) :: qInpResult(:)

    !> Charge difference vector between output and input charges
    real(dp), intent(in) :: qDiff(:)

    real(dp), allocatable :: aa(:,:), bb(:,:)
    integer :: ii, jj

    @:ASSERT(size(qInpResult) == this%nElem)
    @:ASSERT(size(qDiff) == this%nElem)

    if (this%iPrevVector < this%mPrevVector) then
      this%iPrevVector = this%iPrevVector + 1
    end if

    call storeVectors(this%prevQInput, this%prevQDiff, this%indx, &
        &qInpResult, qDiff, this%mPrevVector)

    if (this%tFromStart .or. this%iPrevVector == this%mPrevVector) then

      if (this%tAddIntrpGradient) then
        ! old DIIS estimate for downhill direction points towards current downhill direction as well
        ! as the actual vector, based on P. Briddon comments
        if (dot_product(this%deltaR(:),qDiff) > 0.0_dp) then
          ! mix in larger amounts of the gradient in future
          this%alpha = 1.5_dp*this%alpha
        else
          ! points the other way, mix in less
          this%alpha = 0.5*this%alpha
        end if
      end if

      allocate(aa(this%iPrevVector+1, this%iPrevVector+1))
      allocate(bb(this%iPrevVector+1, 1))

      aa(:,:) = 0.0_dp
      bb(:,:) = 0.0_dp

      do ii = 1, this%iPrevVector
        do jj = 1, this%iPrevVector
          aa(ii, jj) = dot_product( this%prevQDiff(:, ii), &
              & this%prevQDiff(:, jj) )
        end do
      end do
      aa(this%iPrevVector+1, 1:this%iPrevVector) = -1.0_dp
      aa(1:this%iPrevVector, this%iPrevVector+1) = -1.0_dp

      bb(this%iPrevVector+1,1) = -1.0_dp

      ! Solve DIIS system of linear equations
      call gesv(aa, bb)

      qInpResult(:) = 0.0_dp
      do ii = 1, this%iPrevVector
        qInpResult(:) = qInpResult(:) + bb(ii,1) * ( &
            & this%prevQInput(:,ii) + this%prevQDiff(:,ii) )
      end do

      if (this%tAddIntrpGradient) then
        ! add a fraction down the DIIS estimated gradient onto the new solution
        this%deltaR = 0.0_dp
        do ii = 1, this%iPrevVector
          this%deltaR(:) = this%deltaR(:) + bb(ii,1) * this%prevQDiff(:,ii)
        end do
        qInpResult(:) = qInpResult(:) - this%alpha * this%deltaR(:)
      end if

    end if

    if (this%iPrevVector < this%mPrevVector) then
      ! First few iterations return simple mixed vector
      qInpResult(:) = qInpResult(:) + this%initMixParam * qDiff(:)
    end if

  end subroutine DIISMixer_mix


  !> Stores a vector pair in a limited storage. If the stack is full, oldest vector pair is
  !> overwritten.
  subroutine storeVectors(prevQInp, prevQDiff, indx, qInput, qDiff, &
      &mPrevVector)

    !> Contains previous vectors of the first type
    real(dp), intent(inout) :: prevQInp(:,:)

    !> Contains previous vectors of the second type
    real(dp), intent(inout) :: prevQDiff(:,:)

    !> Indexing of data
    integer, intent(inout) :: indx

    !> New first vector
    real(dp), intent(in) :: qInput(:)

    !> New second vector
    real(dp), intent(in) :: qDiff(:)

    !> Size of the stacks.
    integer, intent(in) :: mPrevVector

    indx = mod(indx, mPrevVector) + 1
    prevQInp(:,indx) = qInput(:)
    prevQDiff(:,indx) = qDiff(:)

  end subroutine storeVectors

end module dftbp_diismixer
