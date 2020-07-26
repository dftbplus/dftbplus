!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "error.fypp"

!> Contains MPI related environment settings
module dftbp_mpienv
  use dftbp_accuracy, only : lc
  use dftbp_mpifx
  use dftbp_message
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

    !> Whether current process is the global lead
    logical :: tGlobalLead

    !> Whether current process is the group lead
    logical :: tGroupLead

    !> Number of geometry replicas in the global comm world
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

    !> Whether current process is the replica lead
    logical :: tReplicaLead

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
  ! nReplicas = 1
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
  subroutine TMpiEnv_init(this, nGroup, nReplicas)

    !> Initialised instance on exit
    type(TMpiEnv), intent(out) :: this

    !> Number of process groups to create
    integer, intent(in), optional :: nGroup

    !> Number of structure replicas
    integer, intent(in), optional :: nReplicas

    character(lc) :: tmpStr

    call this%globalComm%init()

    ! number of replicas (structures) in the system
    if (present(nReplicas)) then
      this%nReplicas = nReplicas
    else
      this%nReplicas = 1
    end if
    ! number of groups within a structure (spin, k-points)
    if (present(nGroup)) then
      this%nGroup = nGroup
    else
      this%nGroup = 1
    end if

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

    this%tGlobalLead = this%globalComm%lead
    this%tReplicaLead = this%intraReplicaComm%lead
    this%tGroupLead = this%groupComm%lead

    if (this%tGlobalLead .and. .not. this%tGroupLead) then
      call error("Internal error: Global lead process is not a group leading process")
    end if
    if (this%tGlobalLead .and. .not. this%tReplicaLead) then
      call error("Internal error: Global lead process is not a replica lead process")
    end if
    if (this%tReplicaLead .and. .not. this%tGroupLead) then
      call error("Internal error: Replica lead process is not a group lead process")
    end if

  end subroutine TMpiEnv_init


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

end module dftbp_mpienv
