!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains MPI related environment settings
module mpienv
  use accuracy, only : lc
  use mpifx
  use message
  implicit none
  private

  public :: TMpiEnv, TMpiEnv_init


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

    !> Rank within my group
    integer :: myGroupRank

    !> Global rank of the processes in the given group
    integer, allocatable :: groupMembers(:)

    !> Whether current process is the global master
    logical :: tGlobalMaster

    !> Whether current process is the group master
    logical :: tGroupMaster

    !> Number of replicas in the global comm world
    integer :: nReplicas

    !> Communicator to access processes within a replica
    type(mpifx_comm) :: intraReplicaComm

    !> Communicator between equivalent processors in replicas
    type(mpifx_comm) :: interReplicaComm

    !> Size of the group containing this replica
    integer :: replicaCommSize

    !> Replica group index the current process belongs to (starts with 0)
    integer :: myReplica

    !> Rank within my replica
    integer :: myReplicaRank

    !> Global rank of the processes in the given replica group
    integer, allocatable :: replicaMembers(:)

    !> Whether current process is the replica master
    logical :: tReplicaMaster

  end type TMpiEnv


contains

  !> Initializes MPI environment.
  subroutine TMpiEnv_init(this, nGroup, nReplicas)

    !> Initialised instance on exit
    type(TMpiEnv), intent(out) :: this

    !> Number of process groups to create
    integer, intent(in) :: nGroup

    !> Number of structure replicas
    integer, intent(in) :: nReplicas

    character(lc) :: tmpStr

    call this%globalComm%init()

    ! number of replicas (structures) in the system
    this%nReplicas = nReplicas

    ! number of processors in a replica group
    this%replicaCommSize = this%globalComm%size / this%nReplicas

    if (this%nReplicas * this%replicaCommSize /= this%globalComm%size) then
      write(tmpStr, "(A,I0,A,I0,A)") "Number of replicas (", this%nReplicas,&
          & ") not compatible with number of processes (", this%globalComm%size, ")"
      call error(tmpStr)
    end if

    ! the replica group this process belongs to
    this%myReplica = this%globalComm%rank / this%replicaCommSize

    ! rank within the replica group
    this%myReplicaRank = mod(this%globalComm%rank, this%replicaCommSize)

    ! communicator within this replica
    call this%globalComm%split(this%myReplica, this%myReplicaRank, this%intraReplicaComm)

    ! array of global process ids within this replica
    allocate(this%replicaMembers(this%replicaCommSize))
    call mpifx_allgather(this%intraReplicaComm, this%globalComm%rank, this%replicaMembers)

    ! communicator to equivalent rank processors in different replicas
    call this%globalComm%split(this%myReplicaRank, this%myReplica, this%interReplicaComm)

    ! number of groups of processors within a replica
    this%nGroup = nGroup

    ! number of processors in a group
    this%groupSize = this%intraReplicaComm%size / this%nGroup

    if (this%nGroup * this%groupSize /= this%intraReplicaComm%size) then
      if (this%nReplicas > 1) then
        write(tmpStr, "(A,I0,A,I0,A)") "Number of groups (", this%nGroup,&
            & ") not compatible with number of processes / replica (", this%intraReplicaComm%size,&
            & ")"
      else
        write(tmpStr, "(A,I0,A,I0,A)") "Number of groups (", this%nGroup,&
            & ") not compatible with number of processes (", this%intraReplicaComm%size, ")"
      end if
      call error(tmpStr)
    end if

    ! the group this process belongs to within the replica
    this%myGroup = this%intraReplicaComm%rank / this%groupSize

    ! rank within the group
    this%myGroupRank = mod(this%intraReplicaComm%rank, this%groupSize)

    ! communicator within the group
    call this%intraReplicaComm%split(this%myGroup, this%myGroupRank, this%groupComm)

    ! array of global process ids within this group
    allocate(this%groupMembers(this%groupSize))
    call mpifx_allgather(this%groupComm, this%globalComm%rank, this%groupMembers)
    ! equivalent for ids within the replica
    ! call mpifx_allgather(this%groupComm, this%intraReplicaComm%rank, this%groupMembers)

    ! communicator to equivalent rank processors in different groups within the replica
    call this%intraReplicaComm%split(this%myGroupRank, this%myGroup, this%interGroupComm)

    this%tGlobalMaster = this%globalComm%master
    this%tReplicaMaster = this%intraReplicaComm%master
    this%tGroupMaster = this%groupComm%master

    if (this%tGlobalMaster .and. .not. this%tGroupMaster) then
      call error("Internal error: Global master process is not a group master process")
    end if
    if (this%tGlobalMaster .and. .not. this%tReplicaMaster) then
      call error("Internal error: Global master process is not a replica master process")
    end if
    if (this%tReplicaMaster .and. .not. this%tGroupMaster) then
      call error("Internal error: Replica master process is not a group master process")
    end if

  end subroutine TMpiEnv_init


end module mpienv
