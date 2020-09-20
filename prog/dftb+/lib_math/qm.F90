!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> contains some miscellaneous quantum mechanics related bits and pieces.
module dftbp_qm
  use dftbp_assert
  use dftbp_accuracy, only : dp, rsp

  implicit none

  private
  public :: makeSimilarityTrans

  !> perform a similarity (or unitary) transformation of a matrix X' = U X U^T*
  interface makeSimilarityTrans
    module procedure U_cmplx
    module procedure U_real
  end interface makeSimilarityTrans

contains


  !> unitary transformation of a matrix X' = U X U^T*
  subroutine U_cmplx(xx, uu)

    !> matrix in original basis, U X U^T* on return.
    complex(dp), intent(inout) :: xx(:,:)

    !> unitary matrix
    complex(dp), intent(in) :: uu(:,:)

    complex(dp) :: work(size(xx,dim=1), size(xx,dim=2))

    @:ASSERT(all(shape(xx) == shape(uu)))
    @:ASSERT(size(xx, dim=1) == size(xx, dim=2))

    work(:,:) = matmul(xx, transpose(conjg(uu)))
    xx(:,:) = matmul(uu, work)

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

    work = matmul(xx, transpose(uu))
    xx = matmul(uu, work)

  end subroutine U_real

end module dftbp_qm
