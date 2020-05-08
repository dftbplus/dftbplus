!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains an Anderson mixer
!>
!>   The Anderson mixing is done by building a special average over the
!>   previous input charges and over the previous charge differences
!>   separately and then linear mixing the two averaged vectors with a given
!>   mixing parameter. Only a specified amount of previous charges are
!>   considered.
!> In order to use the mixer you have to create and reset it.
module dftbp_andersonmixer
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_lapackroutines, only : gesv
  implicit none

  private


  !> Contains the necessary data for an Anderson mixer
  !>
  !> For efficiency reasons this derived type also contains the data for the limited memory storage,
  !> which stores a given number of recent vector pairs. The storage should be accessed as an array
  !> with the help of the indx(:) array. Indx(1) gives the index for the most recent stored vector
  !> pairs. (LIFO)
  type TAndersonMixer
    private

    !> General mixing parameter
    real(dp) :: mixParam

    !> Initial mixing parameter
    real(dp) :: initMixParam

    !> Convergence dependent mixing parameters
    real(dp), allocatable :: convMixParam(:,:)

    !> Symmetry breaking parameter
    real(dp) :: omega02

    !> Should symmetry be broken?
    logical :: tBreakSym

    !> Nr. of convergence dependent mixing parameters
    integer :: nConvMixParam

    !> Max. nr. of stored prev. vectors
    integer :: mPrevVector

    !> Nr. of stored previous vectors
    integer :: nPrevVector

    !> Nr. of elements in the vectors
    integer :: nElem

    !> Index array for the storage
    integer, allocatable :: indx(:)

    !> Stored previous input charges
    real(dp), allocatable :: prevQInput(:,:)

    !> Stored prev. charge differences
    real(dp), allocatable :: prevQDiff(:,:)
  end type TAndersonMixer


  !> Creates an AndersonMixer instance
  interface init
    module procedure AndersonMixer_init
  end interface init


  !> Resets the mixer
  interface reset
    module procedure AndersonMixer_reset
  end interface reset


  !> Does the mixing
  interface mix
    module procedure AndersonMixer_mix
  end interface mix

  public :: TAndersonMixer
  public :: init, reset, mix

contains


  !> Creates an Andersom mixer instance.
  subroutine AndersonMixer_init(self, nGeneration, mixParam, initMixParam, convMixParam, omega0)

    !> Initialized Anderson mixer on exit
    type(TAndersonMixer), intent(out) :: self

    !> Nr. of generations (including actual) to consider
    integer, intent(in) :: nGeneration

    !> Mixing parameter for the general case
    real(dp), intent(in) :: mixParam

    !> Mixing parameter for the first nGeneration-1 cycles
    real(dp), intent(in) :: initMixParam

    !> Convergence dependent mixing parameters. Given as 2 by n array of tolerance and mixing
    !> factors pairs. The tolerances (Euclidean norm of the charge diff. vector) must follow each

    !> other in decreasing order. Mixing parameters given here eventually override mixParam or
    !> initMixParam.
    real(dp), intent(in), optional :: convMixParam(:,:)

    !> Symmetry breaking parameter. Diagonal elements of the system of linear equations are
    !> multiplied by (1.0+omega0**2).
    real(dp), intent(in), optional :: omega0

    @:ASSERT(nGeneration >= 2)

    self%nElem = 0
    self%mPrevVector = nGeneration - 1

    allocate(self%prevQInput(self%nElem, self%mPrevVector))
    allocate(self%prevQDiff(self%nElem, self%mPrevVector))
    allocate(self%indx(self%mPrevVector))

    self%mixParam = mixParam
    self%initMixParam = initMixParam
    if (present(convMixParam)) then
      @:ASSERT(size(convMixParam, dim=1) == 2)
      self%nConvMixParam = size(convMixParam, dim=2)
      allocate(self%convMixParam(2, self%nConvMixParam))
      self%convMixParam(:,:) = convMixParam(:,:)
    else
      self%nConvMixParam = 0
    end if
    if (present(omega0)) then
      self%omega02 = omega0**2
      self%tBreakSym = .true.
    else
      self%omega02 = 0.0_dp
      self%tBreakSym = .false.
    end if

  end subroutine AndersonMixer_init


  !> Makes the mixer ready for a new SCC cycle
  subroutine AndersonMixer_reset(self, nElem)

    !> Anderson mixer instance
    type(TAndersonMixer), intent(inout) :: self

    !> Nr. of elements in the vectors to mix
    integer, intent(in) :: nElem

    integer :: ii

    @:ASSERT(nElem > 0)

    if (nElem /= self%nElem) then
      self%nElem = nElem
      deallocate(self%prevQInput)
      deallocate(self%prevQDiff)
      allocate(self%prevQInput(self%nElem, self%mPrevVector))
      allocate(self%prevQDiff(self%nElem, self%mPrevVector))
    end if
    self%nPrevVector = -1
    ! Create index array for accessing elements in the LIFO way
    do ii = 1, self%mPrevVector
      self%indx(ii) = self%mPrevVector + 1 - ii
    end do

  end subroutine AndersonMixer_reset


  !> Mixes charges according to the Anderson method
  subroutine AndersonMixer_mix(self, qInpResult, qDiff)

    !> Anderson mixer
    type(TAndersonMixer), intent(inout) :: self

    !> Input charges on entry, mixed charges on exit.
    real(dp), intent(inout) :: qInpResult(:)

    !> Charge difference
    real(dp), intent(in) :: qDiff(:)

    real(dp), allocatable :: qInpMiddle(:), qDiffMiddle(:)
    real(dp) :: mixParam
    real(dp) :: rTmp
    integer :: ii

    @:ASSERT(size(qInpResult) == self%nElem)
    @:ASSERT(size(qDiff) == self%nElem)

    if (self%nPrevVector < self%mPrevVector) then
      self%nPrevVector = self%nPrevVector + 1
      mixParam = self%initMixParam
    else
      mixParam = self%mixParam
    end if

    ! Determine mixing parameter
    rTmp = sqrt(sum(qDiff**2))
    do ii = self%nConvMixParam, 1, -1
      if (rTmp < self%convMixParam(1, ii)) then
        mixParam = self%convMixParam(2, ii)
        exit
      end if
    end do

    ! First iteration: store vectors and return simple mixed vector
    if (self%nPrevVector == 0) then
      call storeVectors(self%prevQInput, self%prevQDiff, self%indx, &
          &qInpResult, qDiff, self%mPrevVector)
      qInpResult(:) = qInpResult(:) + self%initMixParam * qDiff(:)
      return
    end if

    allocate(qInpMiddle(self%nElem))
    allocate(qDiffMiddle(self%nElem))

    ! Calculate average input charges and average charge differences
    call calcAndersonAverages(qInpMiddle, qDiffMiddle, qInpResult, &
        &qDiff, self%prevQInput, self%prevQDiff, self%nElem, self%nPrevVector, &
        &self%indx, self%tBreakSym, self%omega02)

    ! Store vectors before overwriting qInpResult
    call storeVectors(self%prevQInput, self%prevQDiff, self%indx, &
        &qInpResult, qDiff, self%mPrevVector)

    ! Mix averaged input charge and average charge difference
    qInpResult(:) = qInpMiddle(:) + mixParam * qDiffMiddle(:)

  end subroutine AndersonMixer_mix


  !> Calculates averages input charges and average charge differences according to the Anderson
  !> method.
  !>
  !> Note: The symmetry breaking is not exactly the same as in the paper of Eyert, because here it
  !> is applied to the diagonal of the "original" matrix built from the Fs and not of the "modified"
  !> matrix built from the DFs.
  subroutine calcAndersonAverages(qInpMiddle, qDiffMiddle, qInput, qDiff, prevQInp, prevQDiff, &
      & nElem, nPrevVector, indx, tBreakSym, omega02)

    !> Contains average input charge on exit
    real(dp), intent(out) :: qInpMiddle(:)

    !> Contains averages charge difference on exit
    real(dp), intent(out) :: qDiffMiddle(:)

    !> Input charge in the last iteration
    real(dp), intent(in) :: qInput(:)

    !> Charge difference in the last iteration
    real(dp), intent(in) :: qDiff(:)

    !> Input charges of the previous iterations
    real(dp), intent(in) :: prevQInp(:,:)

    !> Charge differences of the previous iterations
    real(dp), intent(in) :: prevQDiff(:,:)

    !> Nr. of elements in the charge vectors
    integer, intent(in) :: nElem

    !> Nr. of previous iterations stored
    integer, intent(in) :: nPrevVector

    !> Index array describing the reverse storage order
    integer, intent(in) :: indx(:)

    !> If symmetry of linear equation system should be broken
    logical, intent(in) :: tBreakSym

    !> Symmetry breaking constant
    real(dp), intent(in) :: omega02

    real(dp), allocatable :: aa(:,:), bb(:,:), tmp1(:)
    real(dp) :: tmp2
    integer :: ii, jj

    @:ASSERT(size(qInpMiddle) == nElem)
    @:ASSERT(size(qDiffMiddle) == nElem)
    @:ASSERT(size(qInput) == nElem)
    @:ASSERT(size(qDiff) == nElem)
    @:ASSERT(size(prevQInp, dim=1) == nElem)
    @:ASSERT(size(prevQInp, dim=2) >= nPrevVector)
    @:ASSERT(size(prevQDiff, dim=1) == nElem)
    @:ASSERT(size(prevQDiff, dim=2) >= nPrevVector)
    @:ASSERT(size(indx) >= nPrevVector)

    allocate(aa(nPrevVector, nPrevVector))
    allocate(bb(nPrevVector, 1))
    allocate(tmp1(nElem))

    ! Build the system of linear equations
    ! a(i,j) = <F(m)|F(m)-F(m-i)> - <F(m-j)|F(m)-F(m-i)>  (F ~ qDiff)
    ! b(i)   = <F(m)|F(m)-F(m-i)>                         (m ~ current iter.)
    ! Index array serves reverse indexing: indx(1) means most recent vector
    do ii = 1, nPrevVector
      tmp1(:) = qDiff(:) - prevQDiff(:, indx(ii))
      tmp2 = dot_product(qDiff, tmp1)
      bb(ii, 1) = tmp2
      do jj = 1, nPrevVector
        aa(ii, jj) = tmp2 - dot_product(prevQDiff(:, indx(jj)), tmp1)
      end do
    end do

    ! Prevent equations from beeing linearly dependent if desired
    if (tBreakSym) then
      tmp2 = (1.0_dp + omega02)
      do ii = 1, nPrevVector
        aa(ii, ii) =  tmp2 * aa(ii, ii)
      end do
    end if

    ! Solve system of linear equations
    call gesv(aa, bb)

    ! Build averages with calculated coefficients
    qDiffMiddle(:) = 0.0_dp
    do ii = 1, nPrevVector
      qDiffMiddle(:) = qDiffMiddle(:) + bb(ii,1) * prevQDiff(:,indx(ii))
    end do
    qDiffMiddle(:) = qDiffMiddle(:) + (1.0_dp - sum(bb(:,1))) * qDiff(:)

    qInpMiddle(:) = 0.0_dp
    do ii = 1, nPrevVector
      qInpMiddle(:) = qInpMiddle(:) + bb(ii,1) * prevQInp(:,indx(ii))
    end do
    qInpMiddle(:) = qInpMiddle(:) + (1.0_dp - sum(bb(:,1))) * qInput(:)

  end subroutine calcAndersonAverages


  !> Stores a vector pair in a limited storage. If the stack is full, the oldest vector pair is
  !> overwritten.
  subroutine storeVectors(prevQInp, prevQDiff, indx, qInput, qDiff, mPrevVector)

    !> Contains previous vectors of the first type
    real(dp), intent(inout) :: prevQInp(:,:)

    !> Contains previous vectors of the second type
    real(dp), intent(inout) :: prevQDiff(:,:)

    !> Indexing array to the stacks
    integer, intent(inout) :: indx(:)

    !> New first vector
    real(dp), intent(in) :: qInput(:)

    !> New second vector
    real(dp), intent(in) :: qDiff(:)

    !> Size of the stacks.
    integer, intent(in) :: mPrevVector

    integer :: tmp

    tmp = indx(mPrevVector)
    indx(2:mPrevVector) = indx(1:mPrevVector-1)
    indx(1) = tmp
    prevQInp(:,indx(1)) = qInput(:)
    prevQDiff(:,indx(1)) = qDiff(:)

  end subroutine storeVectors

end module dftbp_andersonmixer
