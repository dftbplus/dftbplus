!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"
#:include "error.fypp"


!> Contains MPI related environment settings
module dftbp_common_mpienv
  use dftbp_common_accuracy, only : lc
  use dftbp_extlibs_mpifx, only : mpifx_comm, mpifx_allgather
  use dftbp_io_message, only : error
#:if WITH_TRANSPORT
  use dftbp_extlibs_negf, only : negf_cart_init
#:endif
  implicit none
  
  private
  public :: TMpiEnv, TMpiEnv_init, TMpiEnv_final


  !> Contains MPI related environment settings
  type :: TMpiEnv

    !> Global MPI communicator
    type(mpifx_comm) :: globalComm

    !> Communicator to access processes within current group
    type(mpifx_comm) :: groupComm

    !> Communicator to access equivalent processes in other groups
    type(mpifx_comm) :: interGroupComm

    !> Size of the process groups
    integer :: groupSize

    !> Number of processor groups
    integer :: nGroup

    !> Group index of the current process (starts with 0)
    integer :: myGroup

    !> Rank of the processes in the given group (with respect of globalComm)
    integer, allocatable :: groupMembersGlobal(:)

    !> Rank of the processes in the given group (with respect of MPI_COMM_WORLD)
    integer, allocatable :: groupMembersWorld(:)

    !> Whether current process is the global lead
    logical :: tGlobalLead

    !> Whether current process is the group lead
    logical :: tGroupLead

  contains

    procedure :: mpiSerialEnv

  end type TMpiEnv


contains

  !> Initializes MPI environment.
  ! ---------------------------------------------------------------
  ! Initializes global communicator and group communicators
  ! Example:
  ! globalSize = 10
  ! nGroup = 2
  ! groupSize = 5
  !                        rank
  ! globalComm:      0 1 2 3 4 5 6 7 8 9
  ! groupComm:       0 1 2 3 4 0 1 2 3 4
  ! interGroupComm:  0 0 0 0 0 1 1 1 1 1
  ! ---------------------------------------------------------------
  ! SCALAPACK
  ! Different groups handle different kpoints/spin (iKS)
  ! All procs within a group know eigenval(:,iKS)
  ! These are distributed to all other nodes using interGroupComm
  ! eigenvec(:,:,iKS) are used to build the density matrix, DM(:,:,iKS)
  ! DM(:,:,iKS) contains kWeight(iK) and occupation(iKS)
  ! total DM(:,:) is obtained by mpiallreduce with MPI_SUM
  ! ---------------------------------------------------------------
  ! LIBNEGF
  ! Different groups handle different kpoints/spin (iKS)
  ! All procs within a group know densMat(:,:,iKS)
  ! DM(:,:,iKS) contains kWeight(iK) and occupation(iKS)
  ! total DM(:,:) is obtained by mpiallreduce with MPI_SUM
  ! ---------------------------------------------------------------
  subroutine TMpiEnv_init(this, globalMpiComm, nGroup)

    !> Initialised instance on exit
    type(TMpiEnv), intent(out) :: this

    !> The global MPI communicator (assumed to be MPI_COMM_WORLD if not specified)
    type(mpifx_comm), optional, intent(in) :: globalMpiComm

    !> Number of process groups to create
    integer, intent(in), optional :: nGroup

    if (present(globalMpiComm)) then
      this%globalComm = globalMpiComm
    else
      call this%globalComm%init()
    end if
    if (present(nGroup)) then
      this%nGroup = nGroup
    else
      this%nGroup = 1
    end if

    #:if WITH_TRANSPORT
      call setup_subgrids_negf(this)
    #:else
      call setup_subgrids_common(this)
    #:endif

    this%tGlobalLead = this%globalComm%lead
    this%tGroupLead = this%groupComm%lead

    if (this%tGlobalLead .and. .not. this%tGroupLead) then
      call error("Internal error: Global lead process is not a group leading process")
    end if

  end subroutine TMpiEnv_init


  !> Finalises the communicators in the structure supplied here
  subroutine TMpiEnv_final(this)

    !>  Initialised instance.
    type(TMpiEnv), intent(inout) :: this

    call this%interGroupComm%free()
    call this%groupComm%free()

  end subroutine TMpiEnv_final


  !> Routine to check this is a single processor instance, stopping otherwise (useful to call in
  !> purely serial codes to avid multiple copies being invoked with mpirun)
  subroutine mpiSerialEnv(this, iErr)

    !> Instance
    class(TMpiEnv), intent(in) :: this

    !> Optional error flag
    integer, intent(out), optional :: iErr

    if (this%globalComm%size > 1) then

      @:ERROR_HANDLING(iErr, -1, 'This is serial code, but invoked on multiple processors')

    end if

  end subroutine mpiSerialEnv


  !> Sets up subgrids and group communicators used with common (non-NEGF) solvers.
  subroutine setup_subgrids_common(this)

    !> Environment instance
    type(TMpiEnv), intent(inout) :: this

    integer :: myRank, myGroup
    character(lc) :: tmpStr
    type(mpifx_comm) :: mpiCommWorld

    this%groupSize = this%globalComm%size / this%nGroup
    if (this%nGroup * this%groupSize /= this%globalComm%size) then
      write(tmpStr, "(A,I0,A,I0,A)") "Number of groups (", this%nGroup,&
          & ") not compatible with number of processes (", this%globalComm%size, ")"
      call error(tmpStr)
    end if

    this%myGroup = this%globalComm%rank / this%groupSize
    myRank = mod(this%globalComm%rank, this%groupSize)
    call this%globalComm%split(this%myGroup, myRank, this%groupComm)
    allocate(this%groupMembersGlobal(this%groupSize))
    call mpifx_allgather(this%groupComm, this%globalComm%rank, this%groupMembersGlobal)

    ! Make a wrapper around MPI_COMM_WORLD and get group member ids within that descriptor
    call mpiCommWorld%init()
    allocate(this%groupMembersWorld(this%groupSize))
    call mpifx_allgather(this%groupComm, mpiCommWorld%rank, this%groupMembersWorld)

    myGroup = myRank
    myRank = this%myGroup
    call this%globalComm%split(myGroup, myRank, this%interGroupComm)

  end subroutine setup_subgrids_common


#:if WITH_TRANSPORT

  !> Sets up subgrids and group communicators used with NEGF solver.
  !>
  !> Note: it overrides the global communicator with the Cartesian global handler
  !> created by libNEGF as libNEGF uses that. At a later stage, the Cartesian
  !> communicator should be hidden in the libNEGF-DFTB+ interface instead.
  !>
  subroutine setup_subgrids_negf(this)

    !> Environment instance
    type(TMpiEnv), intent(inout) :: this

    type(mpifx_comm) :: cartComm, mpiCommWorld

    call negf_cart_init(this%globalComm, this%nGroup, cartComm, this%groupComm, this%interGroupComm)
    if (this%globalComm%lead .neqv. cartComm%lead) then
      call error("Internal error: Lead process mismatch between Cartesian and global communicator")
    end if
    this%globalComm = cartComm

    this%groupSize = this%groupComm%size
    this%myGroup = this%interGroupComm%rank

    allocate(this%groupMembersGlobal(this%groupSize))
    call mpifx_allgather(this%groupComm, this%globalComm%rank, this%groupMembersGlobal)

    ! Make a wrapper around MPI_COMM_WORLD and get group member ids within that descriptor
    call mpiCommWorld%init()
    allocate(this%groupMembersWorld(this%groupSize))
    call mpifx_allgather(this%groupComm, mpiCommWorld%rank, this%groupMembersWorld)

  end subroutine setup_subgrids_negf

#:endif

end module dftbp_common_mpienv
