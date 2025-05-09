!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set FLAVOURS = [('cmplx', 'complex', 'Cmplx'), ('real', 'real', 'Real')]

!> Contains a modified Broyden mixer.
!! The modified Broyden mixer implemented here is practically the same as the one in the old DFTB
!! code. A detailed description of the method can be found in Johnson's paper.
!! See: D.D. Johnson, PRB 38, 12807 (1988)
!! DOI: 10.1103/PhysRevB.38.12807
!! In order to use the mixer you have to create and reset it.
module dftbp_mixer_broydenmixer
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error
  use dftbp_math_lapackroutines, only : getrf, getrs
  use dftbp_mixer_mixer, only : TMixerCmplx, TMixerReal
  implicit none


  private

  public :: TBroydenMixerInp
  public :: TBroydenMixerReal, TBroydenMixerReal_init
  public :: TBroydenMixerCmplx, TBroydenMixerCmplx_init

  type :: TBroydenMixerInp
    !> Maximum nr. of iterations (max. nr. of vectors to store)
    ! This is calculated in Initprogram and needs to be updated before constructing the mixer
    integer :: maxSccIter = -1

    !> Mixing parameter
    real(dp) :: mixParam = 0.0_dp

    !> Weight for the Jacobi matrix differences
    real(dp) :: omega0 = 0.01_dp

    !> Minimum allowed weight
    real(dp) :: minWeight = 1.0_dp

    !> Maximum allowed weight
    real(dp) :: maxWeight = 1.0e5_dp

    !> Numerator of the weight
    real(dp) :: weightFac = 1.0e-2_dp
  end type TBroydenMixerInp

#:for NAME, TYPE, LABEL in FLAVOURS
  !> Contains the necessary data for a Broyden mixer.
  type, extends (TMixer${LABEL}$) :: TBroydenMixer${LABEL}$
    private

    !> Actual iteration
    integer :: iIter

    !> Nr. of maximal iterations
    integer :: mIter

    !> Nr. of elements in the vectors
    integer :: nElem

    !> Jacobi matrix differences
    real(dp) :: omega0

    !> Mixing parameter
    real(dp) :: alpha

    !> Minimal weight
    real(dp) :: minWeight

    !> Maximal weight
    real(dp) :: maxWeight

    !> Weighting factor (numerator)
    real(dp) :: weightFac

    !> Weights for prev. iterations
    real(dp), allocatable :: ww(:)

    !> Charge difference in last iteration
    ${TYPE}$(dp), allocatable :: qDiffLast(:)

    !> Input charge in last iteration
    ${TYPE}$(dp), allocatable :: qInpLast(:)

    !> Storage for the "a" matrix
    real(dp), allocatable :: aa(:,:)

    !> DF vectors
    ${TYPE}$(dp), allocatable :: dF(:,:)

    !> uu vectors
    ${TYPE}$(dp), allocatable :: uu(:,:)

    contains
      procedure :: reset => TBroydenMixer${LABEL}$_reset
      procedure :: mix1D => TBroydenMixer${LABEL}$_mix
  end type TBroydenMixer${LABEL}$
#:endfor


contains

#:for NAME, TYPE, LABEL in FLAVOURS

  !> Initializes a Broyden mixer instance.
  !! The weight associated with an iteration is calculated as weigthFac/ww where ww is the Euclidean
  !! norm of the charge difference vector. If the calculated weigth is outside of the
  !! [minWeight, maxWeight] region it is replaced with the appropriate boundary value.
  subroutine TBroydenMixer${LABEL}$_init(this, mixerInp)

    !> An initialized Broyden mixer on exit
    type(TBroydenMixer${LABEL}$), intent(out) :: this

    !> Broyden mixer input data struct
    type(TBroydenMixerInp), intent(in) :: mixerInp

    @:ASSERT(mixerInp%maxSccIter > 0)
    @:ASSERT(mixerInp%mixParam > 0.0_dp)
    @:ASSERT(mixerInp%omega0 > 0.0_dp)

    this%nElem = 0
    this%mIter = mixerInp%maxSccIter
    this%alpha = mixerInp%mixParam
    this%omega0 = mixerInp%omega0
    this%minWeight = mixerInp%minWeight
    this%maxWeight = mixerInp%maxWeight
    this%weightFac = mixerInp%weightFac
    allocate(this%ww(mixerInp%maxSccIter-1))
    allocate(this%qInpLast(this%nElem))
    allocate(this%qDiffLast(this%nElem))
    allocate(this%aa(mixerInp%maxSccIter-1, mixerInp%maxSccIter-1))
    allocate(this%dF(this%nElem, mixerInp%maxSccIter - 1))
    allocate(this%uu(this%nElem, mixerInp%maxSccIter - 1))

  end subroutine TBroydenMixer${LABEL}$_init


  !> Makes the mixer ready for a new SCC cycle.
  subroutine TBroydenMixer${LABEL}$_reset(this, nElem)

    !> Broyden mixer instance
    class(TBroydenMixer${LABEL}$), intent(inout) :: this

    !> Length of the vectors to mix
    integer, intent(in) :: nElem

    @:ASSERT(nElem > 0)

    if (nElem /= this%nElem) then
      this%nElem = nElem
      deallocate(this%qInpLast)
      deallocate(this%qDiffLast)
      allocate(this%qInpLast(this%nElem))
      allocate(this%qDiffLast(this%nElem))
      deallocate(this%dF)
      allocate(this%dF(this%nElem, this%mIter - 1))
      deallocate(this%uu)
      allocate(this%uu(this%nElem, this%mIter - 1))
    end if
    this%iIter = 0
    this%ww(:) = 0.0_dp
    this%aa(:,:) = 0.0_dp

  end subroutine TBroydenMixer${LABEL}$_reset


  !> Mixes charges according to the modified Broyden method.
  !!
  !! Warning: The complex-valued Broyden mixer requires flattened hermitian matrices as input.
  !!   You are free to permute the individual elements of the flattened arrays as long as the same
  !!   permutation is applied to qInpResult and qDiff.
  !!   The restriction arises from the assumption that the dot-products of density matrices are
  !!   real-valued (imaginary parts add up to zero due to the hermitian property) and the linear
  !!   system of equations remains real-valued.
  subroutine TBroydenMixer${LABEL}$_mix(this, qInpResult, qDiff)

    !> The Broyden mixer
    class(TBroydenMixer${LABEL}$), intent(inout) :: this

    !> Input charges on entry, mixed charges on exit
    ${TYPE}$(dp), intent(inout) :: qInpResult(:)

    !> Charge difference between output and input charges
    ${TYPE}$(dp), intent(in) :: qDiff(:)

    this%iIter = this%iIter + 1
    if (this%iIter > this%mIter) then
      call error("Broyden mixer: Maximal nr. of steps exceeded")
    end if

    call modifiedBroydenMixing${LABEL}$(qInpResult, this%qInpLast, this%qDiffLast, this%aa,&
        & this%ww, this%iIter, qDiff, this%alpha, this%omega0, this%minWeight, this%maxWeight,&
        & this%weightFac, this%nElem, this%dF, this%uu)

  end subroutine TBroydenMixer${LABEL}$_mix


  !> Does the real work for the Broyden mixer.
  subroutine modifiedBroydenMixing${LABEL}$(qInpResult, qInpLast, qDiffLast, aa, ww, nn, qDiff,&
      & alpha, omega0, minWeight, maxWeight, weightFac, nElem, dF, uu)

    !> Current input charge on entry, mixed charged on exit
    ${TYPE}$(dp), intent(inout) :: qInpResult(:)

    !> Input charge vector of the previous iterations
    ${TYPE}$(dp), intent(inout) :: qInpLast(:)

    !> Charge difference of the previous iteration
    ${TYPE}$(dp), intent(inout) :: qDiffLast(:)

    !> The matrix a (needed for the mixing)
    real(dp), intent(inout) :: aa(:,:)

    !> Weighting factors of the iterations
    real(dp), intent(inout) :: ww(:)

    !> Current iteration number
    integer, intent(in) :: nn

    !> Charge difference of the current iteration
    ${TYPE}$(dp), intent(in) :: qDiff(:)

    !> Mixing parameter
    real(dp), intent(in) :: alpha

    !> Weight for the Jacobi matrix differences
    real(dp), intent(in) :: omega0

    !> Minimal weight allowed
    real(dp), intent(in) :: minWeight

    !> Maximal weight allowed
    real(dp), intent(in) :: maxWeight

    !> Numerator of the weight
    real(dp), intent(in) :: weightFac

    !> Nr. of elements in the vectors
    integer, intent(in) :: nElem

    !> Prev. DFs
    ${TYPE}$(dp), intent(inout) :: dF(:,:)

    !> Prev. U vectors
    ${TYPE}$(dp), intent(inout) :: uu(:,:)

    real(dp), allocatable :: beta(:,:), cc(:)

    ! Current DF or U-vector
    ${TYPE}$(dp), allocatable :: dF_uu(:)

    real(dp) :: invNorm
    integer :: ii, nn_1
    integer, allocatable :: ipiv(:)

    nn_1 = nn - 1

    @:ASSERT(nn > 0)
    @:ASSERT(size(qInpResult) == nElem)
    @:ASSERT(size(qInpLast) == nElem)
    @:ASSERT(size(qDiffLast) == nElem)
    @:ASSERT(size(qDiff) == nElem)
    @:ASSERT(all(shape(aa) >= [nn_1, nn_1]))
    @:ASSERT(size(ww) >= nn_1)

    ! First iteration: simple mix and storage of qInp and qDiff
    if (nn == 1) then
      qInpLast(:) = qInpResult
      qDiffLast(:) = qDiff
      qInpResult(:) = qInpResult + alpha * qDiff
      return
    end if

    allocate(beta(nn_1, nn_1))
    allocate(cc(nn_1))
    allocate(dF_uu(nElem))
    allocate(ipiv(nn_1))

    ! Create weight factor omega for current iteration
    ww(nn_1) = sqrt(dot_product(qDiff, qDiff))
    if (ww(nn_1) > weightFac / maxWeight) then
      ww(nn_1) = weightFac / ww(nn_1)
    else
      ww(nn_1) = maxWeight
    end if
    if (ww(nn_1) < minWeight) then
      ww(nn_1) = minWeight
    end if

    ! Build |DF(m-1)> and (m is the current iteration number)
    dF_uu(:) = qDiff - qDiffLast
    invNorm = sqrt(dot_product(dF_uu, dF_uu))
    invNorm = max(invNorm, epsilon(1.0_dp))
    invNorm = 1.0_dp / invNorm
    dF_uu(:) = invNorm * dF_uu

    ! Build a, beta, c, and gamma
    ! (due to the hermitian property of our density matrices, the dot-products below are real)
    do ii = 1, nn - 2
      aa(ii, nn_1) = dot_product(dF(:,ii), dF_uu)
      aa(nn_1, ii) = aa(ii, nn_1)
      cc(ii) = ww(ii) * dot_product(dF(:,ii), qDiff)
    end do
    aa(nn_1, nn_1) = 1.0_dp
    cc(nn_1) = ww(nn_1) * dot_product(dF_uu, qDiff)

    do ii = 1, nn_1
      beta(:nn_1, ii) = ww(:nn_1) * ww(ii) * aa(:nn_1,ii)
      beta(ii, ii) = beta(ii, ii) + omega0**2
    end do

    ! LU decomposition
    call getrf(beta, ipiv)
    ! Solve system of linear equations by using the LU decomposition
    call getrs(beta, ipiv, cc, trans='t')

    ! Store |dF(m-1)>
    dF(:, nn_1) = dF_uu

    ! Create |u(m-1)>
    dF_uu(:) = alpha * dF_uu + invNorm * (qInpResult - qInpLast)

    ! Save charge vectors before overwriting
    qInpLast(:) = qInpResult
    qDiffLast(:) = qDiff

    ! Build new vector
    qInpResult(:) = qInpResult + alpha * qDiff
    do ii = 1, nn-2
      qInpResult(:) = qInpResult - ww(ii) * cc(ii) * uu(:,ii)
    end do
    qInpResult(:) = qInpResult - ww(nn_1) * cc(nn_1) * dF_uu

    ! Save |u(m-1)>
    uu(:, nn_1) = dF_uu

  end subroutine modifiedBroydenMixing${LABEL}$
#:endfor

end module dftbp_mixer_broydenmixer
