!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("coulomb")
  use dftbp_common_accuracy, only : dp, rsp
  use dftbp_common_environment, only : TEnvironment, TEnvironment_init
  use dftbp_common_globalenv, only : initGlobalEnv, destructGlobalEnv
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_coulomb, only : sumInvR, addInvRPrime, addInvRPrimePrimeClusterAsymm
  implicit none

#:contains

  #:block TEST_FIXTURE("firstDerivative")

    real(dp) :: r0(3,1), r0d(3,1), r1(3,1), q0(1), q1(1), pot(1), e0(3,1), e1(3,1), eNum(3)
    real(dp) :: eGrad(3,3,1), eGradNum(3,3,1)
    real(dp), parameter :: delta = sqrt(epsilon(0.0d0))
    integer :: ii, jj
    type(TEnvironment) :: env

    call initGlobalEnv()
    call TEnvironment_init(env)

  #:contains

    #:block TEST("potential")
      r0(:,:) = 0.0_dp
      r1(:,:) = reshape(real([0,0,1],dp), [3,1])
      q1(:) = 1.0_dp
      call sumInvR(env, 1, 1, r0, r1, q1, pot)
      @:ASSERT(all(abs(pot - 1.0_dp) <= epsilon(1.0_dp)))
    #:endblock

    #:block TEST("potentialTranslate")
      r0(:,:) = 1.0_dp
      r1(:,:) = reshape(real([1,1,2],dp), [3,1])
      q1(:) = 1.0_dp
      call sumInvR(env, 1, 1, r0, r1, q1, pot)
      @:ASSERT(all(abs(pot - 1.0_dp) <= epsilon(1.0_dp)))
    #:endblock

    #:block TEST("field")
      r0(:,:) = 0.0_dp
      r1(:,:) = reshape(real([0,0,1],dp), [3,1])
      q0(:) = -1.0_dp
      q1(:) = 1.0_dp
      eNum(:) = 0.0_dp
      do ii = 1, 3
        do jj = -1, 1, 2
          r0d(:,:) = r0
          r0d(ii,1) = r0d(ii,1) + jj * delta
          call sumInvR(env, 1, 1, r0d, r1, q1, pot)
          eNum(ii) = eNum(ii) + jj * pot(1)
        end do
      end do
      eNum(:) = eNum / (2.0_dp * delta)
      call addInvRPrime(env, 1, 1, r0, r1, q0, q1, e0, e1, .false.)
      @:ASSERT(all(abs(e1(:,1) - eNum) <= epsilon(1.0_dp)))
    #:endblock

    #:block TEST("fieldGradient")
      r0(:,:) = reshape(real([3,7,1],dp), [3,1])
      r1(:,:) = reshape(real([1,2,3],dp), [3,1])
      q0(:) = -1.0_dp
      q1(:) = 1.0_dp
      eGradNum(:,:,:) = 0.0_dp
      do ii = 1, 3
        do jj = -1, 1, 2
          r0d(:,:) = r0
          r0d(ii,1) = r0d(ii,1) + jj * delta
          e0(:,:) = 0.0_dp
          e1(:,:) = 0.0_dp
          call addInvRPrime(env, 1, 1, r0d, r1, q0, q1, e0, e1, .false.)
          eGradNum(:,ii,1) = eGradNum(:,ii,1) + jj * e0(:,1)
        end do
      end do
      eGradNum(:,:,:) = eGradNum / (2.0_dp * delta)

      call addInvRPrimePrimeClusterAsymm(env, 1, 1, r0, r1, q0, q1, eGrad)
      @:ASSERT(all(abs(eGradNum - eGrad) <= epsilon(1.0_rsp)))
    #:endblock

  #:endblock TEST_FIXTURE

#:endblock TEST_SUITE

@:TEST_DRIVER()
