!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains an Anderson mixer
!!* @description
!!*   The Anderson mixing is done by building a special average over the
!!*   previous input charges and over the previous charge differences
!!*   separately and then linear mixing the two averaged vectors with a given
!!*   mixing parameter. Only a specified amount of previous charges are
!!*   considered.
!!* @note In order to use the mixer you have to create and reset it.
module andersonmixer
  use assert
  use accuracy
  use lapackroutines, only : gesv
  implicit none

  private

  !!* Contains the necessary data for an Anderson mixer
  !!* @note
  !!*   For efficiency reasons this derived type contains also the data
  !!*   for the limited storage, which stores a given number of recent vector
  !!*   pairs. The storage should be accessed as an array with the help of
  !!*   the indx(:) array. Indx(1) gives the index for the most recent stored
  !!*   vector pairs. (LIFO)
  type OAndersonMixer
    private
    real(dp) :: mixParam                    !!* General mixing parameter
    real(dp) :: initMixParam                !!* Initial mixing parameter
    real(dp), allocatable :: convMixParam(:,:)  !!* Convergence dependent mixing
                                            !!* parameters
    real(dp) :: omega02                     !!* Symmetry breaking parameter
    logical :: tBreakSym                    !!* Should symmetry be broken?
    integer :: nConvMixParam                !!* Nr. of convergence dependent
                                            !!* mixing parameters
    integer :: mPrevVector                  !!* Max. nr. of stored prev. vectors
    integer :: nPrevVector                  !!* Nr. of stored previous vectors
    integer :: nElem                        !!* Nr. of elements in the vectors
    integer, allocatable :: indx(:)             !!* Index array for the storage
    real(dp), allocatable :: prevQInput(:,:)    !!* Stored previous input charges
    real(dp), allocatable :: prevQDiff(:,:)     !!* Stored prev. charge differences
  end type OAndersonMixer


  !!* Creates an AndersonMixer instance
  interface init
    module procedure AndersonMixer_init
  end interface

  !!* Resets the mixer
  interface reset
    module procedure AndersonMixer_reset
  end interface

  !!* Does the mixing
  interface mix
    module procedure AndersonMixer_mix
  end interface


  public :: OAndersonMixer
  public :: init, reset, mix

contains

  !!* Creates an Andersom mixer instance.
  !!* @param self         Initialized Anderson mixer on exit
  !!* @param nGeneration  Nr. of generations (including actual) to consider
  !!* @param mixParam     Mixing parameter for the general case
  !!* @param initMixParam Mixing parameter for the first nGeneration-1 cycles
  !!* @param convMixParam Convergence dependent mixing parameters. Given as
  !!*   2 by n array of tolerance and mixing factors pairs. The tolerances
  !!*   (Euclidean norm of the charge diff. vector) must follow each other
  !!*   in decreasing order. Mixing parameters given here eventually override
  !!*   mixParam or initMixParam.
  !!* @param omega0       Symmetry breaking parameter. Diagonal elements of the
  !!*   system of linear equations are multiplied by (1.0+omega0**2).
  subroutine AndersonMixer_init(self, nGeneration, mixParam, initMixParam, &
      &convMixParam, omega0)
    type(OAndersonMixer), intent(out) :: self
    integer, intent(in) :: nGeneration
    real(dp), intent(in) :: mixParam
    real(dp), intent(in) :: initMixParam
    real(dp), intent(in), optional :: convMixParam(:,:)
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



  !!* Makes the mixer ready for a new SCC cycle
  !!* @param self  Anderson mixer instance
  !!* @param nElem Nr. of elements in the vectors to mix
  subroutine AndersonMixer_reset(self, nElem)
    type(OAndersonMixer), intent(inout) :: self
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
    !! Create index array for accessing elements in the LIFO way
    do ii = 1, self%mPrevVector
      self%indx(ii) = self%mPrevVector + 1 - ii
    end do

  end subroutine AndersonMixer_reset


  !!* Mixes charges according to the Anderson method
  !!* @param self       Anderson mixer
  !!* @param qInpResult Input charges on entry, mixed charges on exit.
  !!* @param qDiff      Charge difference
  subroutine AndersonMixer_mix(self, qInpResult, qDiff)
    type(OAndersonMixer), intent(inout) :: self
    real(dp), intent(inout) :: qInpResult(:)
    real(dp), intent(in)    :: qDiff(:)

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

    !! Determine mixing parameter
    rTmp = sqrt(sum(qDiff**2))
    do ii = self%nConvMixParam, 1, -1
      if (rTmp < self%convMixParam(1, ii)) then
        mixParam = self%convMixParam(2, ii)
        exit
      end if
    end do

    !! First iteration: store vectors and return simple mixed vector
    if (self%nPrevVector == 0) then
      call storeVectors(self%prevQInput, self%prevQDiff, self%indx, &
        &qInpResult, qDiff, self%mPrevVector)
      qInpResult(:) = qInpResult(:) + self%initMixParam * qDiff(:)
      return
    end if

    allocate(qInpMiddle(self%nElem))
    allocate(qDiffMiddle(self%nElem))

    !! Calculate average input charges and average charge differences
    call calcAndersonAverages(qInpMiddle, qDiffMiddle, qInpResult, &
        &qDiff, self%prevQInput, self%prevQDiff, self%nElem, self%nPrevVector, &
        &self%indx, self%tBreakSym, self%omega02)

    !! Store vectors before overwriting qInpResult
    call storeVectors(self%prevQInput, self%prevQDiff, self%indx, &
        &qInpResult, qDiff, self%mPrevVector)

    !! Mix averaged input charge and average charge difference
    qInpResult(:) = qInpMiddle(:) + mixParam * qDiffMiddle(:)

  end subroutine AndersonMixer_mix



  !!* Calculates averages input charges and average charge differences according
  !!* to the Anderson method.
  !!* @param qInpMiddle  Contains average input charge on exit
  !!* @param qDiffMiddle Contains averages charge difference on exit
  !!* @param qInput      Input charge in the last iteration
  !!* @param qDiff       Charge difference in the last iteration
  !!* @param prevQInp    Input charges of the previous iterations
  !!* @param prevQDiff   Charge differences of the previous iterations
  !!* @param nElem       Nr. of elements in the charge vectors
  !!* @param nPrevVector Nr. of previous iterations stored
  !!* @param indx        Index array describing the reverse storage order
  !!* @param tBreakSym   If symmetry of linear equation system should be broken
  !!* @param omega02     Symmetry breaking constant
  !!* @note The symmetry breaking is not exactly the same as in the paper
  !!*   of Eyert, because here it is applied to the diagonal of the "original"
  !!*   matrix built from the Fs and not of the "modified" matrix built from
  !!*   the DFs.
  subroutine calcAndersonAverages(qInpMiddle, qDiffMiddle, qInput, &
      &qDiff, prevQInp, prevQDiff, nElem, nPrevVector, indx, tBreakSym, omega02)
    real(dp), intent(out) :: qInpMiddle(:)
    real(dp), intent(out) :: qDiffMiddle(:)
    real(dp), intent(in)  :: qInput(:)
    real(dp), intent(in)  :: qDiff(:)
    real(dp), intent(in)  :: prevQInp(:,:)
    real(dp), intent(in)  :: prevQDiff(:,:)
    integer, intent(in) :: nElem
    integer, intent(in) :: nPrevVector
    integer, intent(in) :: indx(:)
    logical, intent(in) :: tBreakSym
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

    !! Build the system of linear equations
    !! a(i,j) = <F(m)|F(m)-F(m-i)> - <F(m-j)|F(m)-F(m-i)>  (F ~ qDiff)
    !! b(i)   = <F(m)|F(m)-F(m-i)>                         (m ~ current iter.)
    !! Index array serves reverse indexing: indx(1) means most recent vector
    do ii = 1, nPrevVector
      tmp1(:) = qDiff(:) - prevQDiff(:, indx(ii))
      tmp2 = dot_product(qDiff, tmp1)
      bb(ii, 1) = tmp2
      do jj = 1, nPrevVector
        aa(ii, jj) = tmp2 - dot_product(prevQDiff(:, indx(jj)), tmp1)
      end do
    end do

    !! Prevent equations from beeing linearly dependent if desired
    if (tBreakSym) then
      tmp2 = (1.0_dp + omega02)
      do ii = 1, nPrevVector
        aa(ii, ii) =  tmp2 * aa(ii, ii)
      end do
    end if

    !! Solve system of linear equations
    call gesv(aa, bb)

    !! Build averages with calculated coefficients
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



  !!* Stores a vector pair in a limited storage. If the stack is full, oldest
  !!* vector pair is overwritten.
  !!* @param prevQInp    Contains previous vectors of the first type
  !!* @param prevQDiff   Contains previous vectors of the second type
  !!* @param indx        Indexing array to the stacks
  !!* @param qInput      New first vector
  !!* @param qDiff       New second vector
  !!* @param mPrevVector Size of the stacks.
  subroutine storeVectors(prevQInp, prevQDiff, indx, qInput, qDiff, &
      &mPrevVector)
    real(dp), intent(inout) :: prevQInp(:,:)
    real(dp), intent(inout) :: prevQDiff(:,:)
    integer, intent(inout) :: indx(:)
    real(dp), intent(in) :: qInput(:)
    real(dp), intent(in) :: qDiff(:)
    integer, intent(in) :: mPrevVector

    integer :: tmp

    tmp = indx(mPrevVector)
    indx(2:mPrevVector) = indx(1:mPrevVector-1)
    indx(1) = tmp
    prevQInp(:,indx(1)) = qInput(:)
    prevQDiff(:,indx(1)) = qDiff(:)

  end subroutine storeVectors


end module andersonmixer
