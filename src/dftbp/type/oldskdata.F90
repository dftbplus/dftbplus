!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains type for representing the data stored in the old SK-file format and subroutines to read
!> that data from file.
module dftbp_type_oldskdata
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_constants, only : amu__au
  use dftbp_common_file, only : TFileDescr, openFile, closeFile
  use dftbp_dftb_hybridxc, only : THybridXcSKTag, hybridXcFunc
  use dftbp_dftb_repulsive_polyrep, only : TPolyRepInp
  use dftbp_dftb_repulsive_splinerep, only : TSplineRepInp
  use dftbp_io_message, only : error
  use dftbp_io_charmanip, only : tolower
  implicit none

  private
  public :: TOldSKData, readFromFile, readSplineRep, inquireHybridXcTag


  !> Represents the Slater-Koster data in an SK file.
  type TOldSKData

    !> Grid separation
    real(dp) :: dist

    !> Nr. of grid points
    integer :: nGrid

    !> Atomic eigenvalues.
    real(dp) :: skSelf(4)

    !> Hubbard Us
    real(dp) :: skHubbU(4)

    !> Occupations
    real(dp) :: skOcc(4)

    !> Mass of the atom
    real(dp) :: mass

    !> Table for H
    real(dp), allocatable :: skHam(:,:)

    !> Table for S
    real(dp), allocatable :: skOver(:,:)
  end type TOldSKData


  !> Reads the data from an SK-file.
  interface readFromFile
    module procedure OldSKData_readFromFile
  end interface readFromFile


  !> Reads the repulsive from an open file
  interface readSplineRep
    module procedure OldSKData_readSplineRep
  end interface readSplineRep


  !> nr. of interactions in the SK-file
  integer, parameter :: nSKInter = 20

  !> nr. of ints. in old SK-file
  integer, parameter :: nSKInterOld = 10


  !> Mapping between old (spd) and new (spdf) interactions in the SK-table
  integer, parameter :: iSKInterOld(nSKInterOld) &
      & = (/8, 9, 10, 13, 14, 15, 16, 18, 19, 20/)

contains


  !> Reads the data from an SK-file.
  subroutine OldSKData_readFromFile(skData, fileName, homo, iSp1, iSp2, splineRepIn, polyRepIn,&
      & hybridXcSK)

    !> Contains the content of the SK-file on exit
    type(TOldSKData), intent(out) :: skData

    !> Name of the file to read the data from
    character(len=*), intent(in) :: fileName

    !> Is it a homonuclear SK-file?
    logical, intent(in) :: homo

    !> Index of 1st interacting species (for error messages only)
    integer, intent(in), optional :: iSp1

    !> Index of 2nd interacting species (for error messages only)
    integer, intent(in), optional :: iSp2

    !> Repulsive spline part of the SK-file.
    type(TSplineRepInp), intent(out), optional :: splineRepIn

    !> Repulsive polynomial part of the SK-file.
    type(TPolyRepInp), intent(out), optional :: polyRepIn

    !> Reads hybrid xc-functional parameter(s) from SK file
    type(THybridXcSKTag), intent(inout), optional :: hybridXcSK

    type(TFileDescr) :: file
    character(lc) :: chDummy

    !> extended format for f orbitals
    logical :: tExtended

    integer :: nShell
    integer :: ii, iGrid
    real(dp) :: rDummy
    real(dp) :: coeffs(2:9), polyCutoff
    integer :: iostat

    @:ASSERT(present(splineRepIn) .eqv. present(iSp1))
    @:ASSERT(present(iSp1) .eqv. present(iSp2))

    call openFile(file, fileName, mode="r", ioStat=iostat)
    call checkIoError(iostat, fileName, "Unable to open file")
    rewind(file%unit)

    read (file%unit, '(A1)', iostat=iostat) chDummy
    call checkIoError(iostat, fileName, "Unable to read 1st line")
    if (chDummy == '@') then
      tExtended = .true.
      nShell = 4
    else
      tExtended = .false.
      nShell = 3
      rewind(file%unit)
    end if

    read (file%unit,*, iostat=iostat) skData%dist, skData%nGrid
    call checkIoError(iostat, fileName, "Unable to read 1st data line")
    skData%nGrid = skData%nGrid - 1
    if (homo) then
      skData%skSelf(nShell+1:) = 0.0_dp;
      skData%skHubbU(nShell+1:) = 0.0_dp;
      skData%skOcc(nShell+1:) = 0.0_dp;
      read (file%unit,*, iostat=iostat) (skData%skSelf(ii), ii = nShell, 1, -1), &
          &rDummy, &
          &(skData%skHubbU(ii), ii = nShell, 1, -1),&
          &(skData%skOcc(ii), ii = nShell, 1, -1)
      call checkIoError(iostat, fileName, "Unable to read 2nd data line")
      read (file%unit,*, iostat=iostat) skData%mass, (coeffs(ii), ii = 2, 9), &
          &polyCutoff, rDummy, (rDummy, ii = 12, 20)
      call checkIoError(iostat, fileName, "Unable to read 3rd data line")
      skData%mass = skData%mass * amu__au  ! convert to atomic units
    else
      read (file%unit,*, iostat=iostat) rDummy, (coeffs(ii), ii = 2, 9),&
          & polyCutoff, (rDummy, ii = 11, 20)
      call checkIoError(iostat, fileName, "Unable to read 1st data line")
    end if

    if (present(polyRepIn)) then
      polyRepIn%polyCoeffs(:) = coeffs(:)
      polyRepIn%cutoff = polyCutoff
    end if

    allocate(skData%skHam(skData%nGrid, nSKInter))
    allocate(skData%skOver(skData%nGrid, nSKInter))
    skData%skHam(:,:) = 0.0_dp
    skData%skOver(:,:) = 0.0_dp
    do iGrid = 1, skData%nGrid
      if (tExtended) then
        read (file%unit, *, iostat=iostat)&
            & (skData%skHam(iGrid, ii), ii = 1, nSKInter),&
            & (skData%skOver(iGrid, ii), ii = 1, nSKInter)
        call checkIoError(iostat, fileName, "Reading error for integrals")
      else
        read(file%unit,*, iostat=iostat)&
            & (skData%skHam(iGrid, iSKInterOld(ii)), ii = 1, nSKInterOld),&
            & (skData%skOver(iGrid, iSKInterOld(ii)), ii = 1, nSKInterOld)
        call checkIoError(iostat, fileName, "Reading error for integrals")
      end if
    end do

    if (present(splineRepIn)) then
      call readSplineRep(file%unit, fileName, splineRepIn, iSp1, iSp2)
    end if

    ! Read hybrid xc-functional parameter(s)
    if (present(hybridXcSK)) then
      call readHybridXcParams(file%unit, fileName, hybridXcSK)
    end if

    call closeFile(file)

  end subroutine OldSKData_readFromFile


  !> Reads hybrid xc-functional parameter(s) from an open file.
  subroutine readHybridXcParams(fp, fname, hybridXcSK)

    !> File identifier
    integer, intent(in) :: fp

    !> File name
    character(len=*), intent(in) :: fname

    !> Hybrid xc-functional parameter(s)
    type(THybridXcSKTag), intent(inout) :: hybridXcSK

    !! Error status
    integer :: iErr

    !! Temporary character storage
    character(lc) :: strDummy

    !! True, if hybrid xc-functional extra tag was found in SK-file
    logical :: isHybridXcTag

    !! Extra tag in SK-files, defining the type of hybrid xc-functional
    integer :: hybridXcTag

    call inquireHybridXcTag(fname, hybridXcTag, fp=fp)
    rewind(fp)

    ! Seek hybrid xc-functional section in SK file
    do
      read(fp, '(A)', iostat=iErr) strDummy
      if (iErr /= 0) then
        isHybridXcTag = .false.
        exit
      elseif (strDummy == "RangeSep" .or. strDummy == "GlobalHybrid") then
        isHybridXcTag = .true.
        exit
      end if
    end do

    if (.not. isHybridXcTag) then
      write(strDummy, "(A,A,A)") "Hybrid xc-functional calculation requested, but SK-file '",&
          & trim(fname), "' is not a suitable parametrization."
      call error(strDummy)
    end if

    if (hybridXcTag == hybridXcFunc%hyb) then
      hybridXcSK%omega = 0.0_dp
      hybridXcSK%camBeta = 0.0_dp
      read(fp, *, iostat=iErr) strDummy, hybridXcSK%camAlpha
      call checkioerror(iErr, fname, "Error in reading hybrid xc-functional parameter(s)")
    elseif (hybridXcTag == hybridXcFunc%lc) then
      hybridXcSK%camAlpha = 0.0_dp
      hybridXcSK%camBeta = 1.0_dp
      read(fp, *, iostat=iErr) strDummy, hybridXcSK%omega
      call checkioerror(iErr, fname, "Error in reading hybrid xc-functional parameter(s)")
    elseif (hybridXcTag == hybridXcFunc%cam) then
      read(fp, *, iostat=iErr) strDummy, hybridXcSK%omega, hybridXcSK%camAlpha, hybridXcSK%camBeta
      call checkioerror(iErr, fname, "Error in reading hybrid xc-functional parameter(s)")
    end if

    if (hybridXcTag /= hybridXcFunc%hyb .and. hybridXcSK%omega < 0.0_dp) then
      write(strDummy, "(A)") "Range-separation parameter is negative."
      call error(strDummy)
    end if

  end subroutine readHybridXcParams


  !> Reads the repulsive from an open file.
  subroutine OldSKData_readsplinerep(fp, fname, splinerepin, isp1, isp2)
    !! File identifier.
    integer, intent(in) :: fp
    !! Name of the file (for printing help messages).
    character(*), intent(in) :: fname
    !! Input structure for the spline repulsives
    type(TSplineRepInp), intent(inout) :: splinerepin
    !! Index of the first species in the repulsive (for messages)
    integer, intent(in), optional :: isp1
    !! Index of the second species in the repulsive (for messsages)
    integer, intent(in), optional :: isp2

    integer :: iostat
    integer :: nint, ii, jj
    character(lc) :: chdummy
    logical :: hasspline
    real(dp), allocatable :: xend(:)

    rewind(fp)

    ! Look for spline
    do
      read(fp, '(A)', iostat=iostat) chdummy
      if (iostat /= 0) then
        hasspline = .false.
        exit
      elseif (chdummy == "Spline") then
        hasspline = .true.
        exit
      end if
    end do

    if (.not. hasspline) then
      write(chdummy, "(A,A,A)") "No spline repulsive found in file '",&
          & trim(fname), "'"
      call error(chdummy)
    end if

    read(fp, *, iostat=iostat) nint, splinerepin%cutoff
    call checkioerror(iostat, fname, "Error in reading nint and cutoff")
    read(fp, *, iostat=iostat) (splinerepin%expcoeffs(ii), ii = 1, 3)
    call checkioerror(iostat, fname, "Error in reading exponential coeffs")
    allocate(splinerepin%xstart(nint))
    allocate(splinerepin%spcoeffs(4, nint - 1))
    allocate(xend(nint))

    do jj = 1, nint - 1
      read(fp, *, iostat=iostat) splinerepin%xstart(jj), xend(jj),&
          & (splinerepin%spcoeffs(ii,jj), ii = 1, 4)
      call checkioerror(iostat, fname, "Error in reading spline coeffs")
    end do
    read(fp, *, iostat=iostat) splinerepin%xstart(nint), xend(nint),&
        & (splinerepin%spLastCoeffs(ii), ii = 1, 6)
    call checkioerror(iostat, fname, "Error in reading last spline coeffs")
    splinerepin%cutoff = xend(nint)
    ! Check on consistenty
    do jj = 2, nint
      if (abs(xend(jj-1) - splinerepin%xstart(jj)) > 1e-8_dp) then
        if (present(isp1) .and. present(isp2)) then
          write(chdummy, "(A,I2,A,I2,A)") "Repulsive not continuous for species&
              & pair ", isp1, "-", isp2, "."
        else
          write(chdummy, "(A)") "Repulsive not continuous."
        end if
        call error(chdummy)
      end if
    end do

  end subroutine OldSKData_readsplinerep


  !> Inquires hybrid xc-functional extra tag of SK-file.
  subroutine inquireHybridXcTag(fname, hybridXcTag, fp)

    !> File name
    character(len=*), intent(in) :: fname

    !> Hybrid xc-functional extra tag, if allocated
    integer, intent(out) :: hybridXcTag

    !> File identifier, if file is already open
    integer, intent(in), optional :: fp

    !! File descriptor
    type(TFileDescr) :: file

    !! Error status
    integer :: iErr

    !! Temporary character storage
    character(lc) :: strDummy

    !> File identifier
    integer :: fd

    !! True, if hybrid xc-functional extra tag was found in SK-file
    logical :: isHybridXcTag

    if (present(fp)) then
      fd = fp
    else
      call openFile(file, fname, mode="r", ioStat=iErr)
      fd = file%unit
      call checkIoError(iErr, fname, "Unable to open file")
    end if

    rewind(fd)

    ! Seek hybrid xc-functional extra tag in SK-file
    do
      read(fd, '(A)', iostat=iErr) strDummy
      if (iErr /= 0) then
        isHybridXcTag = .false.
        exit
      elseif (strDummy == "RangeSep" .or. strDummy == "GlobalHybrid") then
        isHybridXcTag = .true.
        exit
      end if
    end do

    if (isHybridXcTag) then
      read(fd, *, iostat=iErr) strDummy
      call checkIoError(iErr, fname, "Error in reading hybrid xc-functional extra tag and method.")

      select case(tolower(trim(strDummy)))
      case ("hf")
        hybridXcTag = hybridXcFunc%hyb
      case ("lc")
        hybridXcTag = hybridXcFunc%lc
      case ("cam")
        hybridXcTag = hybridXcFunc%cam
      case default
        write(strDummy, "(A,A,A)") "Unknown hybrid xc-functional method in SK-file '",&
            & trim(fname), "'"
        call error(strDummy)
      end select
    else
      hybridXcTag = hybridXcFunc%none
    end if

    if (.not. present(fp)) call closeFile(file)

  end subroutine inquireHybridXcTag


  !> Checks for IO errors and prints message.
  subroutine checkIOError(iErr, fname, msg)

    !> Flag of the IO operation.
    integer, intent(in) :: iErr

    !> Name of the file.
    character(*), intent(in) :: fname

    !> Message to print if IO operation flag is non-zero.
    character(*), intent(in) :: msg

    if (iErr /= 0) then
      call error("IO error in file '" // trim(fname) // "': " // trim(msg))
    end if

  end subroutine checkIOError

end module dftbp_type_oldskdata
