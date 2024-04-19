!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> contains some miscellaneous quantum mechanics related bits and pieces.
module dftbp_math_qm
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error

  implicit none

  private
  public :: makeSimilarityTrans

  !> perform a similarity (or unitary) transformation of a matrix
  interface makeSimilarityTrans
    module procedure U_cmplx
    module procedure U_real
  end interface makeSimilarityTrans

contains


  !> unitary transformation of a matrix X' = U X U^T or X' = U^T X U
  subroutine U_real(xx, uu, side)

    !> matrix in original basis, U X U^T* on return.
    real(dp), intent(inout) :: xx(:,:)

    !> unitary matrix
    real(dp), intent(in) :: uu(:,:)

    !> which transform order to use, i.e. to which side the original unitary is applied
    character(1), intent(in), optional :: side

    real(dp) :: work(size(xx,dim=1), size(xx,dim=2))
    character(1) :: iSide

    if (present(side)) then
      iSide(:) = side
    else
      iSide(:) = 'L'
    end if

    @:ASSERT(all(shape(xx) == shape(uu)))
    @:ASSERT(size(xx, dim=1) == size(xx, dim=2))

    ! should blasify:
    select case(iSide)
    case ('L', 'l')
      work(:,:) = matmul(xx, transpose(uu))
      xx(:,:) = matmul(uu, work)
    case ('R', 'r')
      work(:,:) = matmul(xx, uu)
      xx(:,:) = matmul(transpose(uu), work)
    case default
      call error("Unknown unitary transform request")
    end select

  end subroutine U_real


  !> unitary transformation of a matrix X' = U X U^T* or X' = U^T* X U
  subroutine U_cmplx(xx, uu, side)

    !> matrix in original basis, U X U^T* on return.
    complex(dp), intent(inout) :: xx(:,:)

    !> unitary matrix
    complex(dp), intent(in) :: uu(:,:)

    !> which transform order to use, i.e. which to side the original unitary is applied
    character(1), intent(in), optional :: side

    complex(dp) :: work(size(xx,dim=1), size(xx,dim=2))
    character(1) :: iSide

    if (present(side)) then
      iSide(:) = side
    else
      iSide(:) = 'L'
    end if

    @:ASSERT(all(shape(xx) == shape(uu)))
    @:ASSERT(size(xx, dim=1) == size(xx, dim=2))

    ! should blasify:
    select case(iSide)
    case ('L', 'l')
      work(:,:) = matmul(xx, transpose(conjg(uu)))
      xx(:,:) = matmul(uu, work)
    case ('R', 'r')
      work(:,:) = matmul(xx, uu)
      xx(:,:) = matmul(transpose(conjg(uu)), work)
    case default
      call error("Unknown unitary transform request")
    end select

  end subroutine U_cmplx

end module dftbp_math_qm
