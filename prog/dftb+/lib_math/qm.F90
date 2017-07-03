!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* contains some miscellaneous QM related bits and pieces.
module qm
  use assert
  use accuracy, only : dp

  implicit none

  public

  !!* constructs a commutator
  interface commutator
    module procedure C_
  end interface

  ! perform a unitary transformation on a matrix $X^\prime = U X U^\dag$
  interface unitary
    module procedure U_cmplx
    module procedure U_real
  end interface

contains

  !!* constructs a commutator for given matrices C = [A,B]
  !!* @param C result of commutator
  !!* @param first matrix
  !!* @param second matrix
  subroutine C_(C,A,B)
    complex(dp), intent(out) :: C(:,:)
    complex(dp), intent(in)  :: A(:,:)
    complex(dp), intent(in)  :: B(:,:)

    @:ASSERT(all(shape(C)==shape(A)))
    @:ASSERT(all(shape(C)==shape(B)))
    @:ASSERT(size(C,dim=1)==size(C,dim=2))

    C = matmul(A,B) - matmul(B,A)

  end subroutine C_

  !> unitary transformation on a matrix $X^\prime = U X U^\dag$
  !! \param X matrix in original basis, U X U* on return.
  !! \param U unitary matrix
  !! \todo test that U is actually unitary in an assert_env block
  subroutine U_cmplx(xx, uu)
    complex(dp), intent(inout) :: xx(:,:)
    complex(dp), intent(in) :: uu(:,:)

    complex(dp) :: work(size(xx,dim=1),size(xx,dim=2))

    @:ASSERT(all(shape(xx) == shape(uu)))
    @:ASSERT(size(xx, dim=1) == size(xx, dim=2))

    work = matmul(xx, transpose(conjg(uu)))
    xx = matmul(uu, work)

  end subroutine U_cmplx

  !> unitary transformation on a matrix $X^\prime = U X U^T$
  !! \param X matrix in original basis, U X U^T on return.
  !! \param U unitary matrix
  !! \todo test that U is actually unitary in an assert_env block
  subroutine U_real(xx, uu)
    real(dp), intent(inout) :: xx(:,:)
    real(dp), intent(in) :: uu(:,:)

    real(dp) :: work(size(xx,dim=1),size(xx,dim=2))

    @:ASSERT(all(shape(xx) == shape(uu)))
    @:ASSERT(size(xx, dim=1) == size(xx, dim=2))

    work = matmul(xx, transpose(uu))
    xx = matmul(uu, work)

  end subroutine U_real

end module qm
