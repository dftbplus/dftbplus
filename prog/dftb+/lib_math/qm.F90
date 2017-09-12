!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> contains some miscellaneous quantum mechanics related bits and pieces.
module qm
  use assert
  use accuracy, only : dp, rsp

  implicit none

  private
  public :: unitary

  !> perform a unitary transformation of a matrix X' = U X U^T*
  interface unitary
    module procedure U_cmplx
    module procedure U_real
  end interface


  !> Test if a matrix is unitary (U U^T* = 1)
  interface isunitary
    module procedure isunitary_real
    module procedure isunitary_cmplx
  end interface isunitary

contains


  !> unitary transformation of a matrix X' = U X U^T*
  subroutine U_cmplx(xx, uu)

    !> matrix in original basis, U X U^T* on return.
    complex(dp), intent(inout) :: xx(:,:)

    !> unitary matrix
    complex(dp), intent(in) :: uu(:,:)

    complex(dp) :: work(size(xx,dim=1),size(xx,dim=2))

    @:ASSERT(all(shape(xx) == shape(uu)))
    @:ASSERT(size(xx, dim=1) == size(xx, dim=2))

    @:ASSERT(isunitary(uu,epsilon(1.0_rsp)))

    work = matmul(xx, transpose(conjg(uu)))
    xx = matmul(uu, work)

  end subroutine U_cmplx


  !> unitary transformation of a matrix X' = U X U^T*
  subroutine U_real(xx, uu)

    !> matrix in original basis, U X U^T on return.
    real(dp), intent(inout) :: xx(:,:)

    !> unitary matrix
    real(dp), intent(in) :: uu(:,:)

    real(dp) :: work(size(xx,dim=1),size(xx,dim=2))

    @:ASSERT(all(shape(xx) == shape(uu)))
    @:ASSERT(size(xx, dim=1) == size(xx, dim=2))

    @:ASSERT(isunitary(uu,epsilon(1.0_rsp)))

    work = matmul(xx, transpose(uu))
    xx = matmul(uu, work)

  end subroutine U_real


  !> tests if a matrix is unitary to single precision
  function isunitary_real(U,tol) result (unitary)

    !> matrix to test
    real(dp), intent(in) :: U(:,:)

    !> tollerance for test
    real(dp), intent(in) :: tol

    !> test result
    logical :: unitary

    integer :: ii
    real(dp) :: work(size(U,dim=1),size(U,dim=2))

    @:ASSERT(size(U,dim=1) == size(U,dim=1))
    @:ASSERT(size(U,dim=1)==size(U,dim=1))

    work = matmul(U,transpose(U))

    unitary = .true.
    do ii = 1, size(work,dim=1)
      if ( abs(work(ii,ii) -1.0_dp) > tol) then
        unitary = .false.
        return
      end if
      work(ii,ii) = 0.0_dp
    end do
    if (any(abs(work) > tol )) then
      unitary = .false.
      return
    end if

  end function isunitary_real


  !> tests if a matrix is unitary to single precision
  function isunitary_cmplx(U,tol) result (unitary)

    !> matrix to test
    complex(dp), intent(in) :: U(:,:)

    !> tollerance for test
    real(dp), intent(in) :: tol

    !> test result
    logical :: unitary

    integer :: ii
    complex(dp) :: work(size(U,dim=1),size(U,dim=2))

    @:ASSERT(tol >= 0.0_dp)
    @:ASSERT(size(U,dim=1)==size(U,dim=1))

    work = matmul(U,transpose(conjg(U)))

    unitary = .true.
    do ii = 1, size(work,dim=1)
      if ( abs(work(ii,ii) -1.0_dp) > tol) then
        unitary = .false.
        return
      end if
      work(ii,ii) = 0.0_dp
    end do
    if (any(abs(work) > tol )) then
      unitary = .false.
      return
    end if

  end function isunitary_cmplx

end module qm
