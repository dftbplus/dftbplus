!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> routines for the restart of the time propagation of the density matrix/atoms
module dftbp_timedep_dynamicsrestart
  use dftbp_common_accuracy, only : dp
  use dftbp_common_file, only : TFileDescr, TOpenOptions, openFile, closeFile
  use dftbp_common_status, only : TStatus
  implicit none

  !> version number for restart format, please increment if you change the file format (and consider
  !> adding backward compatibility functionality)
  integer, parameter :: iDumpFormat = 1

  private
  public :: writeRestartFile, readRestartFile

contains

  !> Write to a restart file
  subroutine writeRestartFile(rho, rhoOld, coord, veloc, time, dt, fileName, isAsciiFile, errStatus)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:,:)

    !> Density matrix at previous time step
    complex(dp), intent(in) :: rhoOld(:,:,:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> atomic velocities
    real(dp), intent(in) :: veloc(:,:)

    !> simulation time (in atomic units)
    real(dp), intent(in) :: time

    !> time step being used (in atomic units)
    real(dp), intent(in) :: dt

    !> name of the dump file
    character(len=*), intent(in) :: fileName

    !> Should restart data be written as ascii (cross platform, but potentially lower
    !> reproducibility) or binary files
    logical, intent(in) :: isAsciiFile

    !> operation status
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


  !> read a restart file containing density matrix, overlap, coordinates and time step
  subroutine readRestartFile(rho, rhoOld, coord, veloc, time, dt, fileName, isAsciiFile, errStatus)

    !> Density Matrix
    complex(dp), intent(out) :: rho(:,:,:)

    !> Previous density Matrix
    complex(dp), intent(out) :: rhoOld(:,:,:)

    !> atomic coordinates
    real(dp), intent(out) :: coord(:,:)

    !> Previous simulation elapsed time until restart file writing
    real(dp), intent(out) :: time

    !> time step being currently used (in atomic units) for checking compatibility
    real(dp), intent(in) :: dt

    !> Name of the file to open
    character(*), intent(in) :: fileName

    !> atomic velocities
    real(dp), intent(out) :: veloc(:,:)

    !> Should restart data be read as ascii (cross platform, but potentially lower reproducibility)
    !> or binary files
    logical, intent(in) :: isAsciiFile

    !> operation status
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

end module dftbp_timedep_dynamicsrestart
