!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains a modified Broyden mixer.
!> The modified Broyden mixer implemented here is practicaly the same as the one in the old DFTB
!> code. A detailed description of the method can be found in Johnson's paper.
!> see D.D. Johnson, PRB 38, 12807 (1988)
!> In order to use the mixer you have to create and reset it.
module dftbp_broydenmixer
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_message
  use dftbp_blasroutines, only : ger
  use dftbp_lapackroutines, only : matinv
  implicit none

  private


  !> Contains the necessary data for a Broyden mixer.
  type TBroydenMixer
    private

    !> Actual iteration
    integer :: iIter

    !> Nr. of maximal iterations
    integer :: mIter

    !> Nr. of element in the vectors
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

    !> Charge difference in last iter.
    real(dp), allocatable :: qDiffLast(:)

    !> Input charge in last iteration
    real(dp), allocatable :: qInpLast(:)

    !> Storage for the "a" matrix
    real(dp), allocatable :: aa(:,:)

    !> DF vectors
    real(dp), allocatable :: dF (:,:)

    !> uu vectors
    real(dp), allocatable :: uu(:,:)
  end type TBroydenMixer


  !> Creates Broyden mixer
  interface init
    module procedure BroydenMixer_init
  end interface init


  !> Resets Broyden mixer
  interface reset
    module procedure BroydenMixer_reset
  end interface reset


  !> Does the charge mixing
  interface mix
    module procedure BroydenMixer_mix
  end interface mix


  !> Returns inverse Jacobian
  interface getInverseJacobian
    module procedure BroydenMixer_getInverseJacobian
  end interface getInverseJacobian

  public :: TBroydenMixer
  public :: init, reset, mix, getInverseJacobian

contains


  !> Creates a Broyden mixer instance.
  !> The weight associated with an iteration is calculated as weigthFac/ww where ww is the Euclidian
  !> norm of the charge difference vector. If the calculated weigth is outside of the
  !> [minWeight,maxWeight] region it is replaced with the appropriate boundary value.
  subroutine BroydenMixer_init(self, mIter, mixParam, omega0, minWeight, &
      &maxWeight, weightFac)

    !> an initialized Broyden mixer on exit
    type(TBroydenMixer), intent(out) :: self

    !> Maximum nr. of iterations (max. nr. of vectors to store)
    integer, intent(in) :: mIter

    !> Mixing parameter
    real(dp), intent(in) :: mixParam

    !> Weight for the Jacobi matrix differences
    real(dp), intent(in) :: omega0

    !> Minimal weight allowed
    real(dp), intent(in) :: minWeight

    !> Maximal weight allowed
    real(dp), intent(in) :: maxWeight

    !> Numerator of the weight
    real(dp), intent(in) :: weightFac

    @:ASSERT(mIter > 0)
    @:ASSERT(mixParam > 0.0_dp)
    @:ASSERT(omega0 > 0.0_dp)

    self%nElem = 0
    self%mIter = mIter
    self%alpha = mixParam
    self%omega0 = omega0
    self%minWeight = minWeight
    self%maxWeight = maxWeight
    self%weightFac = weightFac
    allocate(self%ww(mIter-1))
    allocate(self%qInpLast(self%nElem))
    allocate(self%qDiffLast(self%nElem))
    allocate(self%aa(mIter-1, mIter-1))
    allocate(self%dF(self%nElem, mIter - 1))
    allocate(self%uu(self%nElem, mIter - 1))

  end subroutine BroydenMixer_init


  !> Makes the mixer ready for a new SCC cycle
  subroutine BroydenMixer_reset(self, nElem)

    !> Broyden mixer instance
    type(TBroydenMixer), intent(inout) :: self

    !> Length of the vectors to mix
    integer, intent(in) :: nElem

    @:ASSERT(nElem > 0)

    if (nElem /= self%nElem) then
      self%nElem = nElem
      deallocate(self%qInpLast)
      deallocate(self%qDiffLast)
      allocate(self%qInpLast(self%nElem))
      allocate(self%qDiffLast(self%nElem))
      deallocate(self%dF)
      allocate(self%dF(self%nElem, self%mIter - 1))
      deallocate(self%uu)
      allocate(self%uu(self%nElem, self%mIter - 1))
    end if
    self%iIter = 0
    self%ww(:) = 0.0_dp
    self%aa(:,:) = 0.0_dp

  end subroutine BroydenMixer_reset


  !> Mixes charges according to the modified Broyden method
  subroutine BroydenMixer_mix(self, qInpResult, qDiff)

    !> The Broyden mixer
    type(TBroydenMixer), intent(inout) :: self

    !> Input charges on entry, mixed charges on exit.
    real(dp), intent(inout) :: qInpResult(:)

    !> Charge difference between output and input charges
    real(dp), intent(in) :: qDiff(:)

    self%iIter = self%iIter + 1
    if (self%iIter > self%mIter) then
      call error("Broyden mixer: Maximal nr. of steps exceeded")
    end if

    call modifiedBroydenMixing(qInpResult, self%qInpLast, self%qDiffLast, &
        &self%aa, self%ww, self%iIter, qDiff, self%alpha, self%omega0, &
        &self%minWeight, self%maxWeight, self%weightFac, self%nElem, &
        &self%dF, self%uu)

  end subroutine BroydenMixer_mix


  !> Does the real work for the Broyden mixer
  subroutine modifiedBroydenMixing(qInpResult, qInpLast, qDiffLast, aa, ww, &
      &nn, qDiff, alpha, omega0, minWeight, maxWeight, weightFac, nElem, &
      &dF, uu)

    !> Current input charge on entry, mixed charged on exit
    real(dp), intent(inout) :: qInpResult(:)

    !> Input charge vector of the previous iterations
    real(dp), intent(inout) :: qInpLast(:)

    !> Charge difference of the previous iteration
    real(dp), intent(inout) :: qDiffLast(:)

    !> The matrix a (needed for the mixing).
    real(dp), intent(inout) :: aa(:,:)

    !> Weighting factors of the iterations.
    real(dp), intent(inout) :: ww(:)

    !> Current iteration number
    integer, intent(in) :: nn

    !> Charge difference of the current iteration.
    real(dp), intent(in) :: qDiff(:)

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

    !> Prev. DFs.
    real(dp), intent(inout) :: dF(:,:)

    !> Prev. U vectors
    real(dp), intent(inout) :: uu(:,:)

    real(dp), allocatable :: beta(:,:), cc(:,:), gamma(:,:)

    ! Current DF or U-vector
    real(dp), allocatable :: dF_uu(:)

    real(dp) :: invNorm
    integer :: nn_1
    integer :: ii

    nn_1 = nn - 1

    @:ASSERT(nn > 0)
    @:ASSERT(size(qInpResult) == nElem)
    @:ASSERT(size(qInpLast) == nElem)
    @:ASSERT(size(qDiffLast) == nElem)
    @:ASSERT(size(qDiff) == nElem)
    @:ASSERT(all(shape(aa) >= (/ nn_1, nn_1 /)))
    @:ASSERT(size(ww) >= nn_1)

    ! First iteration: simple mix and storage of qInp and qDiff
    if (nn == 1) then
      qInpLast(:) = qInpResult(:)
      qDiffLast(:) = qDiff(:)
      qInpResult(:) = qInpResult(:) + alpha * qDiff(:)
      return
    end if

    allocate(beta(nn_1, nn_1))
    allocate(cc(1, nn_1))
    allocate(gamma(1, nn_1))
    allocate(dF_uu(nElem))

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

    ! Build |DF(m-1)> and  (m is the current iteration number)
    dF_uu(:) = qDiff(:) - qDiffLast(:)
    invNorm = sqrt(dot_product(dF_uu, dF_uu))
    !invNorm = max(invNorm, 1e-12_dp)
    invNorm = max(invNorm, epsilon(1.0_dp))
    invNorm = 1.0_dp / invNorm
    dF_uu(:) = invNorm * dF_uu(:)

    ! Build a, beta, c, and gamma
    do ii = 1, nn - 2
      aa(ii, nn_1) = dot_product(dF(:,ii), dF_uu)
      aa(nn_1, ii) = aa(ii, nn_1)
      cc(1, ii) = ww(ii) * dot_product(dF(:,ii), qDiff)
    end do
    aa(nn_1, nn_1) = 1.0_dp
    cc(1, nn_1) = ww(nn_1) * dot_product(dF_uu, qDiff)

    do ii = 1, nn_1
      beta(:nn-1, ii) = ww(:nn-1) * ww(ii) * aa(:nn-1,ii)
      beta(ii, ii) = beta(ii, ii) + omega0**2
    end do
    call matinv(beta)

    gamma = matmul(cc, beta)

    ! Store |dF(m-1)>
    dF(:, nn_1) = dF_uu

    ! Create |u(m-1)>
    dF_uu(:) = alpha * dF_uu(:) + invNorm * (qInpResult(:) - qInpLast(:))

    ! Save charge vectors before overwriting
    qInpLast(:) = qInpResult(:)
    qDiffLast(:) = qDiff(:)

    ! Build new vector
    qInpResult(:) = qInpResult + alpha * qDiff(:)
    do ii = 1, nn-2
      qInpResult(:) = qInpResult - ww(ii) * gamma(1,ii) * uu(:,ii)
    end do
    qInpResult(:) = qInpResult - ww(nn_1) * gamma(1,nn_1) * dF_uu

    ! Save |u(m-1)>
    uu(:, nn_1) = dF_uu

  end subroutine modifiedBroydenMixing


  !> return inverse of the Jacobian for the mixing process
  subroutine BroydenMixer_getInverseJacobian(self, invJac)

    !> Broyden mixer
    type(TBroydenMixer), intent(inout) :: self

    !> Inverse of the Jacobian
    real(dp), intent(out) :: invJac(:,:)

    integer :: ii, jj, kk, mm, nn
    real(dp), allocatable :: beta(:,:), zeta(:)

    @:ASSERT(all(shape(invJac) == [ self%nElem, self%nElem ]))

    mm = self%iIter - 1
    allocate(beta(mm, mm))
    allocate(zeta(self%nElem))

    ! Calculating G according to eq (14) in Johnsons paper.
    ! NOTE: The equation in the paper is incorrect, as instead of beta_ij
    ! one has to use (beta_ij * omega(i) * omega(j)) to be consistent
    ! with the mixing as given in eq. (15) and used in this mixer.
    do ii = 1, mm
      beta(:mm, ii) = self%ww(:mm) * self%ww(ii) * self%aa(:mm, ii)
      beta(ii, ii) = beta(ii, ii) + self%omega0**2
    end do
    call matinv(beta)
    do ii = 1, mm
      do jj = 1, mm
        beta(ii, jj) = beta(ii, jj) * self%ww(ii) * self%ww(jj)
      end do
    end do

    ! J^{-1}(m) = -G(m) = -G1 ...
    invJac(:,:) = 0.0_dp
    do ii = 1, self%nElem
      invJac(ii, ii) = -self%alpha
    end do

    ! ... + sum_{k=1}^m |Z_k> <dF^(k)| with |Z_k> = sum_{n}^m beta_{kn} |u(n)>
    do kk = 1, mm
      zeta(:) = 0.0_dp
      do nn = 1, mm
        zeta(:) = zeta + beta(kk, nn) * self%uu(:, nn)
      end do
      call ger(invJac, 1.0_dp, zeta, self%dF(:,kk))
    end do

    ! Normalize inverse Jacobian
    do ii = 1, self%nElem
      invJac(:,ii) = (invJac(:,ii) / sum(invJac(:,ii))) * (-self%alpha)
    end do

  end subroutine BroydenMixer_getInverseJacobian

end module dftbp_broydenmixer
