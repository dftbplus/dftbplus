!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains MPI environmental settings for libNEGF.
module dftbp_negfmpienv
  use dftbp_mpienv
  use dftbp_mpifx
  use dftbp_message
  use dftbp_assert
  use libnegf, only : negf_mpi_init, negf_cart_init
  implicit none
  private

  public :: TNegfEnv, TNegfEnv_init


  ! contains libNEGF mpi communicators and possibly other settings
  type :: TNegfEnv

    type(TMpiEnv) :: negfMpiEnv

    logical :: initialized = .false.

  end type TNegfEnv


contains


  !> Initializes NEGF mpi 
  subroutine TNegfEnv_init(this, mpiEnv, nGroups)

    !> Initialized instance at exit.
    type(TNegfEnv), intent(out) :: this

    !> Initialised MPI environment
    type(TMpiEnv), intent(in) :: mpiEnv

    !> number of processors handling k-points and spin
    integer, intent(in) :: nGroups

    associate(myEnv=>this%negfMpiEnv)

    ! Invokes cartesian creation within libNEGF
    ! Caveat: we use the globalComm to store the cartesian grid communicator
    call negf_cart_init(mpiEnv%globalComm, nGroups, myEnv%globalComm, myEnv%groupComm, &
        &  myEnv%intergroupComm)

    !Now initialize TMpiEnv container (reimplementation of TMpiEnv_init)
    myEnv%groupSize = myEnv%groupComm%size
    myEnv%nGroup = myEnv%interGroupComm%size
    myEnv%myGroup = myEnv%interGroupComm%rank
print*,"Group Size:",myEnv%groupSize
print*,"nGroup:",myEnv%nGroup
print*,"myGlobalRank:", myEnv%globalComm%rank
print*,"myGroupRank:", myEnv%groupComm%rank
print*,"myGroup:", myEnv%myGroup
print*,"Check:", mpiEnv%globalComm%lead, myEnv%globalComm%lead

    if (mpiEnv%globalComm%lead .neqv. myEnv%globalComm%lead) then
      call error("Internal error: Global lead process mismatch")
    end if
    myEnv%tGlobalLead = myEnv%globalComm%lead
    myEnv%tGroupLead = myEnv%groupComm%lead
    allocate(myEnv%groupMembers(myEnv%groupSize))
    call mpifx_allgather(myEnv%groupComm, mpiEnv%globalComm%rank, myEnv%groupMembers)

    if (myEnv%tGlobalLead .and. .not. myEnv%tGroupLead) then
      call error("Internal error: Global lead process is not a group leading process")
    end if

    end associate
    
    this%initialized = .true.

  end subroutine TNegfEnv_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns rows and columns for a 2D processor grid closest possible to a square one.
  subroutine getSquareGridParams(nProc, nRow, nCol)

    !> Total number of processors
    integer, intent(in) :: nProc

    !> Number of rows in the 2D processor grid
    integer, intent(out) :: nRow

    !> Number of columns in the 2D processor grid
    integer, intent(out) :: nCol

    do nRow = int(sqrt(real(nProc))), 1, -1
      if (mod(nProc, nRow) == 0) then
        exit
      end if
    end do
    nCol = nProc / nRow

  end subroutine getSquareGridParams


end module dftbp_negfmpienv
