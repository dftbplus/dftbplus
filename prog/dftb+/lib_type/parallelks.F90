!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

!> Provides data structure for parallelising over k-points and spin.
module parallelks
  use environment
  implicit none
  private

  public :: TParallelKS, TParallelKS_init


  !> Contains information about which k-points and spins must be processed by which processor group
  type :: TParallelKS

    !> K-point and spin-channels to be processed by each processor group (groupKS(1:2,iKS,iGroup)).
    !> Note: third index (group index) starts from 0
    integer, allocatable :: groupKS(:,:,:)

    !> Number of (K, S) tuples to process for each group.
    !> Note: array index (group index) starts from 0
    integer, allocatable :: nGroupKS(:)

    !> Maximal number of KS-indices per processor group
    integer :: maxGroupKS

    !> The (K, S) tuples of the local processor group
    integer, allocatable :: localKS(:,:)

    !> Number of local (K, S) tuples to process
    integer :: nLocalKS

  end type TParallelKS

contains
  
  !> Returns the (k-point, spin) tuples to be processed by current processor grid (if parallel) or
  !> put everything in one group if serial.
  subroutine TParallelKS_init(this, env, nKpoint, nSpin)

    !> Initialised instance on exti.
    type(TParallelKS), intent(out) :: this

    !> Environenment settings
    type(TEnvironment), intent(in) :: env

    !> Number of k-points in calculation.
    integer, intent(in) :: nKpoint

    !> Number of spin channels in calculation
    integer, intent(in) :: nSpin

    integer :: nGroup, myGroup, iGroup
    integer :: maxGroupKS, nKS, res
    integer :: iS, iK

    nGroup = env%nGroup
    myGroup = env%myGroup

    nKS = nKpoint * nSpin
    maxGroupKS = nKS / nGroup
    res = nKS - maxGroupKS * nGroup
    if (res > 0) then
      maxGroupKS = maxGroupKS + 1
    end if

    allocate(this%nGroupKS(0 : nGroup - 1))
    this%nGroupKS(:) = 0
    allocate(this%groupKS(2, maxGroupKS, 0 : nGroup - 1))
    this%groupKS(:,:,:) = 0
    do iS = 1, nSpin
      do iK = 1, nKpoint
        iGroup = mod((iS - 1) * nKpoint + iK - 1, nGroup)
        this%nGroupKS(iGroup) = this%nGroupKS(iGroup) + 1
        this%groupKS(:, this%nGroupKS(iGroup), iGroup) = [iK, iS]
      end do
    end do
    this%maxGroupKS = maxGroupKS
    this%nLocalKS = this%nGroupKS(myGroup)
    this%localKS = this%groupKS(:, 1:this%nLocalKS, myGroup)

  end subroutine TParallelKS_init


end module parallelks
