!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Higher level abstraction for wrapped MPI memory windows
module dftbp_type_sharedmemory
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : MPIFX_SIZE_T, mpifx_win
#:endif

  implicit none (type, external)

  private
  public :: TArray2D, TArray2D_init

  type :: TArray2D

    integer, pointer :: data(:,:) => null()
  #:if WITH_MPI
    integer, pointer, private :: dataMemory_(:) => null()
    type(mpifx_win), private :: dataWindow_
  #:endif

  contains

    final :: TArray2D_final

  end type TArray2D


contains

  subroutine TArray2D_init(this, env, memoryDim)

    type(TArray2D), intent(out) :: this

    type(TEnvironment), intent(in) :: env

    integer, intent(in) :: memoryDim(:)

    integer :: dataLength
  #:if WITH_MPI
    integer(MPIFX_SIZE_T) :: longDataLength
  #:endif

    dataLength = product(memoryDim)

  #:if WITH_MPI
    longDataLength = int(dataLength, kind=MPIFX_SIZE_T)
    call this%dataWindow_%allocate_shared(env%mpi%nodeComm, longDataLength, this%dataMemory_)

    this%data(1:memoryDim(1), 1:memoryDim(2)) => this%dataMemory_(1:dataLength)

  #:else

    allocate(this%data(memoryDim(1), memoryDim(2)))

  #:endif

  end subroutine TArray2D_init


  subroutine TArray2D_final(this)

    type(TArray2D), intent(inout) :: this

  #:if WITH_MPI

    call this%dataWindow_%free()
    this%dataMemory_ => null()

  #:else

    deallocate(this%data)

  #:endif

    this%data => null()

  end subroutine TArray2D_final

end module dftbp_type_sharedmemory
