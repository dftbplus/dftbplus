!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains type for representing the data stored in the old SK-file format and subroutines to read
!! that data from file.
module dftbp_type_oldskdata
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_constants, only : amu__au
  use dftbp_common_file, only : closeFile, openFile, TFileDescr
  use dftbp_dftb_hybridxc, only : hybridXcFunc, THybridXcSKTag
  use dftbp_dftb_repulsive_polyrep, only : TPolyRepInp
  use dftbp_dftb_repulsive_splinerep, only : TSplineRepInp
  use dftbp_io_charmanip, only : tolower
  use dftbp_io_message, only : error
  implicit none

  private
  public :: TOldSKData, readFromFile, readSplineRep, parseHybridXcTag


  !> Represents the Slater-Koster data in an SK file.
  type TOldSKData

    !> Grid separation
    real(dp) :: dist

    !> Nr. of grid points
    integer :: nGrid

    !> Atomic eigenvalues
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
    module procedure TOldSKData_readFromFile
  end interface readFromFile


  !> Reads the repulsive from an open file
  interface readSplineRep
    module procedure TOldSKData_readSplineRep
  end interface readSplineRep


  !> Nr. of interactions in the SK-file
  integer, parameter :: nSKInter = 20

  !> Nr. of ints. in old SK-file
  integer, parameter :: nSKInterOld = 10


  !> Mapping between old (spd) and new (spdf) interactions in the SK-table
  integer, parameter :: iSKInterOld(nSKInterOld) = [8, 9, 10, 13, 14, 15, 16, 18, 19, 20]


contains

  !> Reads the data from an SK-file.
  subroutine TOldSKData_readFromFile(skData, fileName, homo, iSp1, iSp2, splineRepInp, polyRepInp,&
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
    type(TSplineRepInp), intent(out), optional :: splineRepInp

    !> Repulsive polynomial part of the SK-file.
    type(TPolyRepInp), intent(out), optional :: polyRepInp

    !> Reads hybrid xc-functional parameter(s) from SK file
    type(THybridXcSKTag), intent(inout), optional :: hybridXcSK

    !! Type of hybrid xc-functional found in the SK-file
    integer :: hybridXcType

    type(TFileDescr) :: file
    character(lc) :: chDummy

    !! Extended format for f orbitals
    logical :: tExtended

    integer :: nShell
    integer :: ii, iGrid
    real(dp) :: rDummy
    real(dp) :: coeffs(2:9), polyCutoff
    integer :: iostat

    @:ASSERT(present(splineRepInp) .eqv. present(iSp1))
    @:ASSERT(present(iSp1) .eqv. present(iSp2))

    call openFile(file, fileName, mode="r", ioStat=iostat)
    call checkIoError(iostat, fileName, "Unable to open file")
    rewind(file%unit)

    read(file%unit, "(A1)", iostat=iostat) chDummy
    call checkIoError(iostat, fileName, "Unable to read 1st line")
    if (chDummy == "@") then
      tExtended = .true.
      nShell = 4
    else
      tExtended = .false.
      nShell = 3
      rewind(file%unit)
    end if

    read(file%unit, *, iostat=iostat) skData%dist, skData%nGrid
    call checkIoError(iostat, fileName, "Unable to read 1st data line")
    skData%nGrid = skData%nGrid - 1
    if (homo) then
      skData%skSelf(nShell+1:) = 0.0_dp
      skData%skHubbU(nShell+1:) = 0.0_dp
      skData%skOcc(nShell+1:) = 0.0_dp
      read(file%unit, *, iostat=iostat) (skData%skSelf(ii), ii = nShell, 1, -1), rDummy,&
          & (skData%skHubbU(ii), ii = nShell, 1, -1), (skData%skOcc(ii), ii = nShell, 1, -1)
      call checkIoError(iostat, fileName, "Unable to read 2nd data line")
      read(file%unit, *, iostat=iostat) skData%mass, (coeffs(ii), ii = 2, 9), polyCutoff, rDummy,&
          & (rDummy, ii = 12, 20)
      call checkIoError(iostat, fileName, "Unable to read 3rd data line")
      ! convert to atomic units
      skData%mass = skData%mass * amu__au
    else
      read(file%unit, *, iostat=iostat) rDummy, (coeffs(ii), ii = 2, 9), polyCutoff,&
          & (rDummy, ii = 11, 20)
      call checkIoError(iostat, fileName, "Unable to read 1st data line")
    end if

    if (present(polyRepInp)) then
      polyRepInp%polyCoeffs(:) = coeffs
      polyRepInp%cutoff = polyCutoff
    end if

    allocate(skData%skHam(skData%nGrid, nSKInter), source=0.0_dp)
    allocate(skData%skOver(skData%nGrid, nSKInter), source=0.0_dp)
    do iGrid = 1, skData%nGrid
      if (tExtended) then
        read(file%unit, *, iostat=iostat)&
            & (skData%skHam(iGrid, ii), ii = 1, nSKInter),&
            & (skData%skOver(iGrid, ii), ii = 1, nSKInter)
        call checkIoError(iostat, fileName, "Reading error for integrals")
      else
        read(file%unit, *, iostat=iostat)&
            & (skData%skHam(iGrid, iSKInterOld(ii)), ii = 1, nSKInterOld),&
            & (skData%skOver(iGrid, iSKInterOld(ii)), ii = 1, nSKInterOld)
        call checkIoError(iostat, fileName, "Reading error for integrals")
      end if
    end do

    if (present(splineRepInp)) then
      call readSplineRep(file%unit, fileName, splineRepInp, iSp1, iSp2)
    end if

    ! Read hybrid xc-functional parameter(s)
    if (present(hybridXcSK)) then
      call parseHybridXcTag(fileName, fp=file%unit, hybridXcType=hybridXcType,&
          & hybridXcSK=hybridXcSK)
      if (hybridXcType == hybridXcFunc%none) then
        write(chDummy, "(A,A,A)") "Hybrid xc-functional calculation requested, but SK-file '",&
            & trim(fileName), "' is not a suitable parametrization."
        call error(chDummy)
      end if
    end if

    call closeFile(file)

  end subroutine TOldSKData_readFromFile


  !> Reads the repulsive from an open file.
  subroutine TOldSKData_readsplinerep(fp, fname, splineRepInp, iSp1, iSp2)

    !> File identifier
    integer, intent(in) :: fp

    !> Name of the file (for printing help messages)
    character(*), intent(in) :: fname

    !> Input structure for the spline repulsives
    type(TSplineRepInp), intent(inout) :: splineRepInp

    !> Index of the first species in the repulsive (for messages)
    integer, intent(in), optional :: iSp1

    !> Index of the second species in the repulsive (for messsages)
    integer, intent(in), optional :: iSp2

    !! Error status
    integer :: iostat

    integer :: nint, ii, jj
    character(lc) :: chdummy
    logical :: hasspline
    real(dp), allocatable :: xend(:)

    rewind(fp)

    ! Look for spline
    do
      read(fp, "(A)", iostat=iostat) chdummy
      if (iostat /= 0) then
        hasspline = .false.
        exit
      elseif (chdummy == "Spline") then
        hasspline = .true.
        exit
      end if
    end do

    if (.not. hasspline) then
      write(chdummy, "(A,A,A)") "No spline repulsive found in file '", trim(fname), "'"
      call error(chdummy)
    end if

    read(fp, *, iostat=iostat) nint, splineRepInp%cutoff
    call checkioerror(iostat, fname, "Error in reading nint and cutoff")
    read(fp, *, iostat=iostat) (splineRepInp%expcoeffs(ii), ii = 1, 3)
    call checkioerror(iostat, fname, "Error in reading exponential coeffs")
    allocate(splineRepInp%xstart(nint))
    allocate(splineRepInp%spcoeffs(4, nint - 1))
    allocate(xend(nint))

    do jj = 1, nint - 1
      read(fp, *, iostat=iostat) splineRepInp%xstart(jj), xend(jj),&
          & (splineRepInp%spcoeffs(ii,jj), ii = 1, 4)
      call checkioerror(iostat, fname, "Error in reading spline coeffs")
    end do
    read(fp, *, iostat=iostat) splineRepInp%xstart(nint), xend(nint),&
        & (splineRepInp%spLastCoeffs(ii), ii = 1, 6)
    call checkioerror(iostat, fname, "Error in reading last spline coeffs")
    splineRepInp%cutoff = xend(nint)
    ! Check on consistenty
    do jj = 2, nint
      if (abs(xend(jj-1) - splineRepInp%xstart(jj)) > 1e-8_dp) then
        if (present(iSp1) .and. present(iSp2)) then
          write(chdummy, "(A,I2,A,I2,A)") "Repulsive not continuous for species pair ",&
              & iSp1, "-", iSp2, "."
        else
          write(chdummy, "(A)") "Repulsive not continuous."
        end if
        call error(chdummy)
      end if
    end do

  end subroutine TOldSKData_readsplinerep


  !> Reads hybrid xc-functional parameter(s) from an SK-file, by opening the file or accessing an
  !! already open file using its identifier (if provided).
  subroutine parseHybridXcTag(fname, fp, hybridXcTag, hybridXcType, hybridXcSK)

    !> File name
    character(len=*), intent(in) :: fname

    !> File identifier, if file is already open
    integer, intent(in), optional :: fp

    !> Hybrid xc-functional extra tag found in the SK-file
    integer, intent(out), optional :: hybridXcTag

    !> Type of hybrid xc-functional found in the SK-file
    integer, intent(out), optional :: hybridXcType

    !> Hybrid xc-functional parameter(s)
    type(THybridXcSKTag), intent(inout), optional :: hybridXcSK

    !! Error status
    integer :: iErr

    !! Temporary character storage
    character(lc) :: strDummy, tag

    !! True, if hybrid xc-functional extra tag was found in SK-file
    logical :: isHybridXcTag

    !! Hybrid xc-functional parameter(s)
    type(THybridXcSKTag) :: hybridXcSK_

    !! Hybrid xc-functional extra tag found in the SK-file
    integer :: hybridXcTag_

    !! Type of hybrid xc-functional found in the SK-file
    integer :: hybridXcType_

    !! File descriptor
    type(TFileDescr) :: file

    !! File identifier
    integer :: fd

    if (present(fp)) then
      fd = fp
    else
      call openFile(file, fname, mode="r", ioStat=iErr)
      fd = file%unit
      call checkIoError(iErr, fname, "Unable to open file")
    end if

    rewind(fd)

    ! Seek hybrid xc-functional section in SK-file
    do
      read(fd, "(A)", iostat=iErr) strDummy
      if (iErr /= 0) then
        isHybridXcTag = .false.
        exit
      elseif (strDummy == "RangeSep") then
        isHybridXcTag = .true.
        exit
      end if
    end do

    if (isHybridXcTag) then
      read(fd, "(A)", iostat=iErr) strDummy
      call checkIoError(iErr, fname, "Error in reading hybrid xc-functional extra tag and method.")
      read(strDummy, *, iostat=iErr) tag
      call checkIoError(iErr, fname, "Error in reading hybrid xc-functional extra tag and method.")

      select case(tolower(trim(tag)))
      case ("lc")
        hybridXcTag_ = hybridXcFunc%lc
      case ("cam")
        hybridXcTag_ = hybridXcFunc%cam
      case default
        write(strDummy, "(A,A,A)") "Unknown hybrid xc-functional method in SK-file '",&
            & trim(fname), "'"
        call error(strDummy)
      end select
    else
      hybridXcTag_ = hybridXcFunc%none
    end if

    select case(hybridXcTag_)
    case (hybridXcFunc%lc)
      hybridXcSK_%camAlpha = 0.0_dp
      hybridXcSK_%camBeta = 1.0_dp
      read(strDummy, *, iostat=iErr) tag, hybridXcSK_%omega
      call checkIoError(iErr, fname, "Error in reading hybrid xc-functional parameters.")
      hybridXcType_ = hybridXcFunc%lc
    case (hybridXcFunc%cam)
      read(strDummy, *, iostat=iErr) tag, hybridXcSK_%omega, hybridXcSK_%camAlpha,&
          & hybridXcSK_%camBeta
      call checkIoError(iErr, fname, "Error in reading hybrid xc-functional parameters.")
      if (abs(hybridXcSK_%camAlpha) < epsilon(1.0_dp)&
          & .and. abs(hybridXcSK_%camBeta) < epsilon(1.0_dp)) then
        hybridXcType_ = hybridXcFunc%none
      elseif (abs(hybridXcSK_%camAlpha) < epsilon(1.0_dp)&
          & .and. abs(hybridXcSK_%camBeta - 1.0_dp) < epsilon(1.0_dp)) then
        hybridXcType_ = hybridXcFunc%lc
      elseif (abs(hybridXcSK_%camBeta) < epsilon(1.0_dp)) then
        hybridXcType_ = hybridXcFunc%hyb
      else
        hybridXcType_ = hybridXcFunc%cam
      end if
    case default
      hybridXcType_ = hybridXcFunc%none
      hybridXcSK_%omega = 0.0_dp
      hybridXcSK_%camAlpha = 0.0_dp
      hybridXcSK_%camBeta = 0.0_dp
    end select

    if (hybridXcSK_%omega < 0.0_dp) then
      write(strDummy, "(A)") "Range-separation parameter is negative."
      call error(strDummy)
    end if

    if (.not. present(fp)) call closeFile(file)

    ! hand over requested information
    if (present(hybridXcTag)) hybridXcTag = hybridXcTag_
    if (present(hybridXcType)) hybridXcType = hybridXcType_
    if (present(hybridXcSK)) hybridXcSK = hybridXcSK_

  end subroutine parseHybridXcTag


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
