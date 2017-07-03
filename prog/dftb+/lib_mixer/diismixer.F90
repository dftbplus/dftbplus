!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains an DIIS mixer
!!* @description
!!*   The DIIS mixing is done by building a weighted combination over the
!!*   previous input charges to minimise the residue of the error.  Only
!!*   a specified number of previous charge vectors are considered.
!!* The modification based on from Kovalenko et al. and Patrick Briddon to add
!!* a contribution from the gradient vector as well is also used
!!* @note In order to use the mixer you have to create and reset it.
!!* @ref Kovalenko et al. J. Comput. Chem., 20: 928–936 1999
module diismixer
  use assert
  use accuracy
  use lapackroutines, only : gesv
  implicit none

  private

  !!* Contains the necessary data for an DIIS mixer
  type ODIISMixer
    private
    real(dp) :: initMixParam                !!* Initial mixing parameter
    integer :: mPrevVector                  !!* Max. nr. of stored prev. vectors
    integer :: iPrevVector                  !!* Nr. of stored previous vectors
    integer :: nElem                        !!* Nr. of elements in the vectors
    integer :: indx                         !!* Index for the storage
    real(dp), allocatable :: prevQInput(:,:)    !!* Stored previous input charges
    real(dp), allocatable :: prevQDiff(:,:)     !!* Stored prev. charge differences
    logical :: tFromStart                   !!* True if DIIS used from iteration
                                            !!* 2 as well as mixing
    logical  :: tAddIntrpGradient           !!* force modification for gDIIS?
    real(dp) :: alpha                       !!* Alpha factor to add in new
                                            !!* information
    real(dp), allocatable :: deltaR(:)          !!* Holds DIIS mixed gradients from
                                            !!* older iterations for downhill
                                            !!* direction
  end type ODIISMixer


  !!* Creates an DIISMixer instance
  interface init
    module procedure DIISMixer_init
  end interface

  !!* Resets the mixer
  interface reset
    module procedure DIISMixer_reset
  end interface

  !!* Does the mixing
  interface mix
    module procedure DIISMixer_mix
  end interface


  public :: ODIISMixer
  public :: init, reset, mix

contains

  !!* Creates an DIIS mixer instance.
  !!* @param self         Pointer to an initialized DIIS mixer on exit
  !!* @param nGeneration  Nr. of generations (including actual) to consider
  !!* @param initMixParam Mixing parameter for the first nGeneration cycles
  !!* @param tFromStart   True if using DIIS from iteration 2 as well as mixing
  subroutine DIISMixer_init(self, nGeneration, initMixParam,tFromStart,alpha)
    type(ODIISMixer), intent(out) :: self
    integer, intent(in)            :: nGeneration
    real(dp), intent(in)           :: initMixParam
    logical, intent(in), optional  :: tFromStart
    real(dp), intent(in), optional :: alpha

    @:ASSERT(nGeneration >= 2)

    self%nElem = 0
    self%mPrevVector = nGeneration

    allocate(self%prevQInput(self%nElem, self%mPrevVector))
    allocate(self%prevQDiff(self%nElem, self%mPrevVector))

    self%initMixParam = initMixParam

    if (present(tFromStart)) then
      self%tFromStart = tFromStart
    else
      self%tFromStart = .false.
    end if

    if (present(alpha)) then
      self%tAddIntrpGradient = .true.
      self%alpha = alpha
      allocate(self%deltaR(self%nElem))
    else
      self%tAddIntrpGradient = .false.
      self%alpha = 0.0_dp
      allocate(self%deltaR(0))
    end if

    self%deltaR = 0.0_dp

  end subroutine DIISMixer_init


  !!* Makes the mixer ready for a new SCC cycle
  !!* @param self  DIIS mixer instance
  !!* @param nElem Nr. of elements in the vectors to mix
  subroutine DIISMixer_reset(self, nElem)
    type(ODIISMixer), intent(inout) :: self
    integer, intent(in) :: nElem


    @:ASSERT(nElem > 0)

    if (nElem /= self%nElem) then
      self%nElem = nElem
      deallocate(self%prevQInput)
      deallocate(self%prevQDiff)
      allocate(self%prevQInput(self%nElem, self%mPrevVector))
      allocate(self%prevQDiff(self%nElem, self%mPrevVector))
      if (self%tAddIntrpGradient) then
        deallocate(self%deltaR)
        allocate(self%deltaR(self%nElem))
        self%deltaR = 0.0_dp
      end if
    end if
    self%iPrevVector = 0
    self%indx = 0

  end subroutine DIISMixer_reset


  !!* Mixes charges according to the DIIS method
  !!* @param self       Pointer to the diis mixer
  !!* @param qInpResult Input charges on entry, mixed charges on exit.
  !!* @param qDiff      Charge difference
  subroutine DIISMixer_mix(self, qInpResult, qDiff)
    type(ODIISMixer), intent(inout) :: self
    real(dp), intent(inout) :: qInpResult(:)
    real(dp), intent(in)    :: qDiff(:)

    real(dp), allocatable :: aa(:,:), bb(:,:)
    integer :: ii, jj

    @:ASSERT(size(qInpResult) == self%nElem)
    @:ASSERT(size(qDiff) == self%nElem)

    if (self%iPrevVector < self%mPrevVector) then
      self%iPrevVector = self%iPrevVector + 1
    end if

    call storeVectors(self%prevQInput, self%prevQDiff, self%indx, &
        &qInpResult, qDiff, self%mPrevVector)

    if (self%tFromStart .or. self%iPrevVector == self%mPrevVector) then

      if (self%tAddIntrpGradient) then
        if (dot_product(self%deltaR(:),qDiff) > 0.0_dp) then ! old DIIS estimate
          ! for downhill direction points towards current downhill direction as
          ! well as the actual vector, based on P. Briddon comments
          self%alpha = 1.5_dp*self%alpha ! mix in larger amounts of the gradient
          ! in future
        else
          self%alpha = 0.5*self%alpha ! points the other way, mix in less
        end if
        ! write(*,*)'Alpha ',self%alpha
      end if

      allocate(aa(self%iPrevVector+1, self%iPrevVector+1))
      allocate(bb(self%iPrevVector+1, 1))

      aa(:,:) = 0.0_dp
      bb(:,:) = 0.0_dp

      do ii = 1, self%iPrevVector
        do jj = 1, self%iPrevVector
          aa(ii, jj) = dot_product( self%prevQDiff(:, ii), &
              & self%prevQDiff(:, jj) )
        end do
      end do
      aa(self%iPrevVector+1, 1:self%iPrevVector) = -1.0_dp
      aa(1:self%iPrevVector, self%iPrevVector+1) = -1.0_dp

      bb(self%iPrevVector+1,1) = -1.0_dp

      !! Solve system of linear equations
      call gesv(aa, bb)

      qInpResult(:) = 0.0_dp
      do ii = 1, self%iPrevVector
        qInpResult(:) = qInpResult(:) + bb(ii,1) * ( &
            & self%prevQInput(:,ii) + self%prevQDiff(:,ii) )
      end do

      if (self%tAddIntrpGradient) then ! add a fraction down the DIIS estimated
                                       ! gradient onto the new solution
        self%deltaR = 0.0_dp
        do ii = 1, self%iPrevVector
          self%deltaR(:) = self%deltaR(:) + bb(ii,1) * self%prevQDiff(:,ii)
        end do
        qInpResult(:) = qInpResult(:) - self%alpha * self%deltaR(:)
      end if

    end if

    if (self%iPrevVector < self%mPrevVector) then
      !! First few iterations return simple mixed vector
      qInpResult(:) = qInpResult(:) + self%initMixParam * qDiff(:)
    end if

  end subroutine DIISMixer_mix

  !!* Stores a vector pair in a limited storage. If the stack is full, oldest
  !!* vector pair is overwritten.
  !!* @param prevQInp    Contains previous vectors of the first type
  !!* @param prevQDiff   Contains previous vectors of the second type
  !!* @param indx        Indexing of data
  !!* @param qInput      New first vector
  !!* @param qDiff       New second vector
  !!* @param mPrevVector Size of the stacks.
  subroutine storeVectors(prevQInp, prevQDiff, indx, qInput, qDiff, &
      &mPrevVector)
    real(dp), intent(inout) :: prevQInp(:,:)
    real(dp), intent(inout) :: prevQDiff(:,:)
    integer, intent(inout) :: indx
    real(dp), intent(in) :: qInput(:)
    real(dp), intent(in) :: qDiff(:)
    integer, intent(in) :: mPrevVector

    indx = mod(indx, mPrevVector) + 1
    prevQInp(:,indx) = qInput(:)
    prevQDiff(:,indx) = qDiff(:)

  end subroutine storeVectors

end module diismixer
