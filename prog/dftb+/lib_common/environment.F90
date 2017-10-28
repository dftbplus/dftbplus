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
  public :: withMpi, withScalapack


  !> Contains environment settings.
  type :: TEnvironment
    private

    !> Whether this process is the master?
    logical, public :: tMaster = .true.

    !> Whether this process is supposed to do I/O
    logical, public :: tIoProc = .true.

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


  !> Whether code was compiled with MPI support
  logical, parameter :: withMpi = ${FORTRAN_LOGICAL(WITH_MPI)}$

  !> Whether code was compiled with Scalapack
  logical, parameter :: withScalapack = ${FORTRAN_LOGICAL(WITH_SCALAPACK)}$


contains


#:if WITH_MPI

  !> Initializes MPI environment.
  subroutine initMpi(this)
    class(TEnvironment), intent(inout) :: this

    ! MPI settings
    call TMpiEnv_init(this%mpi)
    this%tMaster = this%mpi%all%master
    this%tIoProc = this%tMaster

  end subroutine initMpi

#:endif


#:if WITH_SCALAPACK

  !> Initializes BLACS environment
  subroutine initBlacs(this, rowBlock, colBlock, nGroup, nOrb, nAtom)

    !> Instance
    class(TEnvironment), intent(inout) :: this

    !> Row block size
    integer, intent(in) :: rowBlock

    !> Column block size
    integer, intent(in) :: colBlock

    !> Nr. of processor groups.
    integer, intent(in) :: nGroup

    !> Nr. of orbitals
    integer, intent(in) :: nOrb

    !> Nr. of atoms
    integer, intent(in) :: nAtom

    call TBlacsEnv_init(this%blacs, rowBlock, colBlock, nGroup, nOrb, nAtom)

  end subroutine initBlacs

#:endif


end module environment
