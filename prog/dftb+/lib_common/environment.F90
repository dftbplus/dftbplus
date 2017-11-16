!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains computer environment settings
module environment
#:if WITH_MPI
  use mpienv
#:endif
#:if WITH_SCALAPACK
  use blacsenv
#:endif
  implicit none
  private

  public :: TEnvironment


  !> Contains environment settings.
  type :: TEnvironment

    !> Whether this process is the master?
    logical, public :: tGlobalMaster = .true.

    !> Nr. of groups in the system
    integer, public :: nGroup = 1

    !> Id of current group (starts with 0)
    integer, public :: myGroup = 0

  #:if WITH_MPI
    !> Global mpi settings
    type(TMpiEnv), public :: mpi
  #:endif
  #:if WITH_SCALAPACK
    !> Global scalapack settings
    type(TBlacsEnv), public :: blacs
  #:endif

  contains
  #:if WITH_MPI
    procedure :: initMpi
  #:endif
  #:if WITH_SCALAPACK
    procedure :: initBlacs
  #:endif

  end type TEnvironment


contains


#:if WITH_MPI

  !> Initializes MPI environment.
  subroutine initMpi(this, nGroup)

    !> Instance
    class(TEnvironment), intent(inout) :: this

    !> Number of process groups to create
    integer, intent(in) :: nGroup

    ! MPI settings
    call TMpiEnv_init(this%mpi, nGroup)
    this%tGlobalMaster = this%mpi%tGlobalMaster
    this%nGroup = this%mpi%nGroup
    this%myGroup = this%mpi%myGroup

  end subroutine initMpi

#:endif


#:if WITH_SCALAPACK

  !> Initializes BLACS environment
  subroutine initBlacs(this, rowBlock, colBlock, nOrb, nAtom)

    !> Instance
    class(TEnvironment), intent(inout) :: this

    !> Row block size
    integer, intent(in) :: rowBlock

    !> Column block size
    integer, intent(in) :: colBlock

    !> Nr. of orbitals
    integer, intent(in) :: nOrb

    !> Nr. of atoms
    integer, intent(in) :: nAtom

    call TBlacsEnv_init(this%blacs, this%mpi, rowBlock, colBlock, nOrb, nAtom)

  end subroutine initBlacs

#:endif


end module environment
