!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Routines for the restart of the time propagation of the density matrix/atoms
module dftbp_timedep_dynamicsrestart
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_file, only : TFileDescr, TOpenOptions, openFile, closeFile
  use dftbp_common_status, only : TStatus
  #:if WITH_SCALAPACK
  use dftbp_extlibs_scalapackfx, only : linecomm, DLEN_, M_, N_, scalafx_getdescriptor
  use dftbp_type_densedescr, only: TDenseDescr
  use dftbp_type_commontypes, only : TParallelKS
  #:endif
  implicit none

  !> Version number for restart format, please increment if you change the file format (and consider
  !! adding backward compatibility functionality)
  integer, parameter :: iDumpFormat = 1

  private
  public :: writeRestartFile, readRestartFile
  #:if WITH_SCALAPACK
  public :: writeRestartFileBlacs, readRestartFileBlacs
  #:endif
  
contains

  !> Write to a restart file.
  subroutine writeRestartFile(rho, rhoOld, coord, veloc, time, dt, fileName, isAsciiFile, errStatus)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:,:)

    !> Density matrix at previous time step
    complex(dp), intent(in) :: rhoOld(:,:,:)

    !> Atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Atomic velocities
    real(dp), intent(in) :: veloc(:,:)

    !> Simulation time (in atomic units)
    real(dp), intent(in) :: time

    !> Time step being used (in atomic units)
    real(dp), intent(in) :: dt

    !> Name of the dump file
    character(len=*), intent(in) :: fileName

    !> Should restart data be written as ascii (cross platform, but potentially lower
    !! reproducibility) or binary files
    logical, intent(in) :: isAsciiFile

    !> Operation status
    type(TStatus), intent(out) :: errStatus

    type(TFileDescr) :: fd
    integer :: ii, jj, kk, iErr
    character(len=120) :: error_string

    if (isAsciiFile) then
      call openFile(fd, trim(fileName) // '.dat', mode="w", iostat=iErr)
    else
      ! Set to stream explicitely, as it was written as stream from the beginning
      call openFile(fd, trim(fileName) // '.bin',&
          & options=TOpenOptions(form='unformatted', access='stream', action='write'), iostat=iErr)
    end if

    if (iErr /= 0) then
      if (isAsciiFile) then
        write(error_string, "(A,A,A)") "Failure to open external restart file ",trim(fileName),&
            & ".dat for writing"
      else
        write(error_string, "(A,A,A)") "Failure to open external restart file ",trim(fileName),&
            & ".bin for writing"
      end if
      @:RAISE_ERROR(errStatus, -1, error_string)
    end if

    if (isAsciiFile) then

      write(fd%unit, *)iDumpFormat
      write(fd%unit, *)size(rho, dim=1), size(rho, dim=3), size(coord, dim=2), time, dt
      do ii = 1, size(rho, dim=3)
        do jj = 1, size(rho, dim=2)
          do kk = 1, size(rho, dim=1)
            write(fd%unit, *)rho(kk,jj,ii)
          end do
        end do
      end do
      do ii = 1, size(rhoOld, dim=3)
        do jj = 1, size(rhoOld, dim=2)
          do kk = 1, size(rhoOld, dim=1)
            write(fd%unit, *)rhoOld(kk,jj,ii)
          end do
        end do
      end do
      do ii = 1, size(coord, dim=2)
        write(fd%unit, *)coord(:,ii)
      end do
      do ii = 1, size(veloc, dim=2)
        write(fd%unit, *)veloc(:,ii)
      end do

    else

      write(fd%unit) iDumpFormat
      write(fd%unit) size(rho, dim=1), size(rho, dim=3), size(coord, dim=2), time, dt
      write(fd%unit) rho, rhoOld, coord, veloc

    end if

    call closeFile(fd)

  end subroutine writeRestartFile


  !> Read a restart file containing density matrix, overlap, coordinates and time step
  subroutine readRestartFile(rho, rhoOld, coord, veloc, time, dt, fileName, isAsciiFile, errStatus)

    !> Density Matrix
    complex(dp), intent(out) :: rho(:,:,:)

    !> Previous density Matrix
    complex(dp), intent(out) :: rhoOld(:,:,:)

    !> Atomic coordinates
    real(dp), intent(out) :: coord(:,:)

    !> Previous simulation elapsed time until restart file writing
    real(dp), intent(out) :: time

    !> Time step being currently used (in atomic units) for checking compatibility
    real(dp), intent(in) :: dt

    !> Name of the file to open
    character(*), intent(in) :: fileName

    !> Atomic velocities
    real(dp), intent(out) :: veloc(:,:)

    !> Should restart data be read as ascii (cross platform, but potentially lower reproducibility)
    !! or binary files
    logical, intent(in) :: isAsciiFile

    !> Operation status
    type(TStatus), intent(out) :: errStatus

    type(TFileDescr) :: fd
    integer :: ii, jj, kk, nOrb, nSpin, nAtom, version, iErr
    real(dp) :: deltaT
    logical :: isExisting
    character(len=120) :: error_string

    if (isAsciiFile) then
      inquire(file=trim(fileName)//'.dat', exist=isExisting)
      if (.not. isExisting) then
        error_string = "TD restart file " // trim(fileName)//'.dat' // " is missing"
        @:RAISE_ERROR(errStatus, -1, error_string)
      end if
    else
      inquire(file=trim(fileName)//'.bin', exist=isExisting)
      if (.not. isExisting) then
        error_string = "TD restart file " // trim(fileName)//'.bin' // " is missing"
        @:RAISE_ERROR(errStatus, -1, error_string)
      end if
    end if

    if (isAsciiFile) then
      call openFile(fd, trim(fileName)//'.dat', mode="r", iostat=iErr)
    else
      ! Set to stream explicitely, as it was written as stream from the beginning
      call openFile(fd, file=trim(fileName)//'.bin',&
          & options=TOpenOptions(form='unformatted', access='stream', action='read',&
          & position="rewind"), iostat=iErr)
    end if

    if (iErr /= 0) then
      if (isAsciiFile) then
        write(error_string, "(A,A,A)") "Failure to open external tddump file",trim(fileName), ".dat"
      else
        write(error_string, "(A,A,A)") "Failure to open external tddump file",trim(fileName), ".bin"
      end if
      @:RAISE_ERROR(errStatus, -1, error_string)
    end if

    if (isAsciiFile) then
      read(fd%unit, *)version
      if (version /= iDumpFormat) then
        @:RAISE_ERROR(errStatus, -1, "Unknown TD restart format")
      end if
      read(fd%unit, *) nOrb, nSpin, nAtom, time, deltaT
      if (nOrb /= size(rho, dim=1)) then
        write(error_string, "(A,I0,A,I0)")"Incorrect number of orbitals, ",nOrb,&
            & " in tddump file, should be ",size(rho, dim=1)
        @:RAISE_ERROR(errStatus, -1, error_string)
      end if
      if (nSpin /= size(rho, dim=3)) then
        write(error_string, "(A,I1,A,I1)")"Incorrect number of spin channels, ",nSpin,&
            & " in tddump file, should be ",size(rho, dim=3)
        @:RAISE_ERROR(errStatus, -1, error_string)
      end if
      if (nAtom /= size(coord, dim=2)) then
        write(error_string, "(A,I0,A,I0)")"Incorrect number of atoms, ",nAtom,&
            & " in tddump file, should be ", size(coord, dim=2)
        @:RAISE_ERROR(errStatus, -1, error_string)
      end if
      if (abs(deltaT - dt) > epsilon(0.0_dp)) then
        write(error_string, "(A,E14.8,A,E14.8)")"Restart file generated for time step",&
            & deltaT, " instead of current timestep of", dt
      end if
      do ii = 1, size(rho, dim=3)
        do jj = 1, size(rho, dim=2)
          do kk = 1, size(rho, dim=1)
            read(fd%unit, *)rho(kk,jj,ii)
          end do
        end do
      end do
      do ii = 1, size(rhoOld, dim=3)
        do jj = 1, size(rhoOld, dim=2)
          do kk = 1, size(rhoOld, dim=1)
            read(fd%unit, *)rhoOld(kk,jj,ii)
          end do
        end do
      end do
      do ii = 1, size(coord, dim=2)
        read(fd%unit, *)coord(:,ii)
      end do
      do ii = 1, size(veloc, dim=2)
        read(fd%unit, *)veloc(:,ii)
      end do
    else
      read(fd%unit)version
      if (version /= iDumpFormat) then
        @:RAISE_ERROR(errStatus, -1, "Unknown TD restart format")
      end if
      read(fd%unit) nOrb, nSpin, nAtom, time, deltaT
      if (nOrb /= size(rho, dim=1)) then
        write(error_string, "(A,I0,A,I0)")"Incorrect number of orbitals, ",nOrb,&
            & " in tddump file, should be ",size(rho, dim=1)
        @:RAISE_ERROR(errStatus, -1, error_string)
      end if
      if (nSpin /= size(rho, dim=3)) then
        write(error_string, "(A,I1,A,I1)")"Incorrect number of spin channels, ",nSpin,&
            & " in tddump file, should be ",size(rho, dim=3)
        @:RAISE_ERROR(errStatus, -1, error_string)
      end if
      if (nAtom /= size(coord, dim=2)) then
        write(error_string, "(A,I0,A,I0)")"Incorrect number of atoms, ",nAtom,&
            & " in tddump file, should be ", size(coord, dim=2)
        @:RAISE_ERROR(errStatus, -1, error_string)
      end if
      if (abs(deltaT - dt) > epsilon(0.0_dp)) then
        write(error_string, "(A,E14.8,A,E14.8)")"Restart file generated for time step",&
            & deltaT, " instead of current timestep of", dt
      end if
      read(fd%unit) rho, rhoOld, coord, veloc
    end if
    call closeFile(fd)

  end subroutine readRestartFile


  #:if WITH_SCALAPACK
  !> Write to a restart file the DM in distributed format
  subroutine writeRestartFileBlacs(rho, rhoOld, coord, veloc, time, dt, fileName, env, &
      & denseDesc, parallelKS, errStatus)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:,:)

    !> Density matrix at previous time step
    complex(dp), intent(in) :: rhoOld(:,:,:)

    !> Atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Atomic velocities
    real(dp), intent(in) :: veloc(:,:)

    !> Simulation time (in atomic units)
    real(dp), intent(in) :: time

    !> Time step being used (in atomic units)
    real(dp), intent(in) :: dt

    !> Name of the dump file
    character(len=*), intent(in) :: fileName

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Container for parallel distribution pattern of k-points and spins
    type(TParallelKS), allocatable :: parallelKS

    !> Operation status
    type(TStatus), intent(out) :: errStatus

    type(TFileDescr) :: fd
    integer :: iErr, nOrb, nSpin, icol, iKS, irow, icoll
    type(linecomm) :: collector
    complex(dp), allocatable :: localRhoCol(:)

    nOrb = denseDesc%fullSize
    nSpin = parallelKS%nLocalKS
    if (env%mpi%tGlobalLead) then
      call openFile(fd, trim(fileName) // '.bin',&
          & options=TOpenOptions(form='unformatted', access='stream', action='write'))
      write(fd%unit) iDumpFormat
      write(fd%unit) nOrb, nSpin, size(coord, dim=2), time, dt
      write(fd%unit) coord, veloc
    end if

    nOrb = denseDesc%fullSize
    allocate(localRhoCol(nOrb))

    call collector%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")
    do iKS = 1, parallelKS%nLocalKS
      do icol = 1, nOrb
        if (env%mpi%tGlobalLead) then
          call collector%getline_lead(env%blacs%orbitalGrid, icol, rho(:,:,iKS), localRhoCol)
          write(fd%unit)localRhoCol(:)
        else
          call collector%getline_follow(env%blacs%orbitalGrid, icol, rho(:,:,iKS))
        end if
      end do
    end do

    do iKS = 1, parallelKS%nLocalKS
      do icol = 1, nOrb
        if (env%mpi%tGlobalLead) then
          call collector%getline_lead(env%blacs%orbitalGrid, icol, rhoOld(:,:,iKS), localRhoCol)
          write(fd%unit)localRhoCol(:)
        else
          call collector%getline_follow(env%blacs%orbitalGrid, icol, rhoOld(:,:,iKS))
        end if
      end do
    end do

    if (env%mpi%tGlobalLead) then
      call closeFile(fd)
    end if

  end subroutine writeRestartFileBlacs


  !> Write to a restart file the DM in Blacs format
  subroutine readRestartFileBlacs(rho, rhoOld, coord, veloc, time, dt, fileName, env, &
      & denseDesc, parallelKS, errStatus)

    !> Density Matrix
    complex(dp), intent(out) :: rho(:,:,:)

    !> Previous density Matrix
    complex(dp), intent(out) :: rhoOld(:,:,:)

    !> Atomic coordinates
    real(dp), intent(out) :: coord(:,:)

    !> Previous simulation elapsed time until restart file writing
    real(dp), intent(out) :: time

    !> Time step being currently used (in atomic units) for checking compatibility
    real(dp), intent(in) :: dt

    !> Name of the file to open
    character(*), intent(in) :: fileName

    !> Atomic velocities
    real(dp), intent(out) :: veloc(:,:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Container for parallel distribution pattern of k-points and spins
    type(TParallelKS), allocatable :: parallelKS

    !> Operation status
    type(TStatus), intent(out) :: errStatus

    type(TFileDescr) :: fd
    integer :: iErr, version, nOrb, nSpin, nAtom, icol, iKS, irow, icoll
    character(len=120) :: error_string
    real(dp) :: deltaT
    type(linecomm) :: distributor
    complex(dp), allocatable :: localRhoCol(:)

    if (env%mpi%tGlobalLead) then
      call openFile(fd, file=trim(fileName)//'.bin',&
          & options=TOpenOptions(form='unformatted', access='stream', action='read',&
          & position="rewind"))

      read(fd%unit) version
      read(fd%unit) nOrb, nSpin, nAtom, time, deltaT
      read(fd%unit) coord, veloc
    end if

    nOrb = denseDesc%fullSize
    allocate(localRhoCol(nOrb))

    call distributor%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")
    do iKS = 1, parallelKS%nLocalKS
      do icol = 1, nOrb
        if (env%mpi%tGlobalLead) then
          read(fd%unit) localRhoCol(:)
          call distributor%setline_lead(env%blacs%orbitalGrid, icol, localRhoCol(:), rho(:,:,iKS))
        else
          call distributor%setline_follow(env%blacs%orbitalGrid, icol, rho(:,:,iKS))
        end if
      end do
    end do

    do iKS = 1, parallelKS%nLocalKS
      do icol = 1, nOrb
        if (env%mpi%tGlobalLead) then
          read(fd%unit) localRhoCol(:)
          call distributor%setline_lead(env%blacs%orbitalGrid, icol, localRhoCol(:), rhoOld(:,:,iKS))
        else
          call distributor%setline_follow(env%blacs%orbitalGrid, icol, rhoOld(:,:,iKS))
        end if
      end do
    end do

    if (env%mpi%tGlobalLead) then
      call closeFile(fd)
    end if

  end subroutine readRestartFileBlacs
  #:endif

end module dftbp_timedep_dynamicsrestart
