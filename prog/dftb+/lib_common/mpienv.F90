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
  use dftbp_mpifx, only : mpifx_comm, mpifx_allgather
  use dftbp_message, only : error
  #:if WITH_TRANSPORT
    use negf_int, only : negf_setup_mpi_communicator
  #:endif
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

    !> Global rank of the processes in the given group
    integer, allocatable :: groupMembers(:)

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
  subroutine TMpiEnv_init(this, nGroup)

    !> Initialised instance on exit
    type(TMpiEnv), intent(out) :: this

    !> Number of process groups to create
    integer, intent(in), optional :: nGroup

    call this%globalComm%init()
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

    this%groupSize = this%globalComm%size / this%nGroup
    if (this%nGroup * this%groupSize /= this%globalComm%size) then
      write(tmpStr, "(A,I0,A,I0,A)") "Number of groups (", this%nGroup,&
          & ") not compatible with number of processes (", this%globalComm%size, ")"
      call error(tmpStr)
    end if

    this%myGroup = this%globalComm%rank / this%groupSize
    myRank = mod(this%globalComm%rank, this%groupSize)
    call this%globalComm%split(this%myGroup, myRank, this%groupComm)
    allocate(this%groupMembers(this%groupSize))
    call mpifx_allgather(this%groupComm, this%globalComm%rank, this%groupMembers)

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

    type(mpifx_comm) :: cartComm

    call negf_cart_init(this%globalComm, this%nGroup, cartComm, this%groupComm, this%interGropComm)
    if (this%globalComm%lead .neqv. cartComm%lead) then
      call error("Internal error: Lead process mismatch between Cartesian and global communicator")
    end if
    this%globalComm = cartComm

    this%groupSize = this%groupComm%size
    this%myGroup = this%interGropComm%rank

    allocate(this%groupMembers(this%groupSize))
    call mpifx_allgather(this%groupComm, this%globalComm%rank, this%groupMembers)

  end subroutine setup_subgrids_negf

#:endif

end module dftbp_mpienv
