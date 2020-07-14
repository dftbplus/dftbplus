!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Read GBSA parametrisation data from file
module dftbp_gbsafile
  use dftbp_accuracy, only : dp
  use dftbp_born, only : TGBInput
  use dftbp_charmanip, only : newline, whiteSpaces
  use dftbp_constants, only : lc, amu__au, kg__au, AA__Bohr, kcal_mol__Hartree, &
    & symbolToNumber
  use dftbp_hsdutils, only : detailedError, detailedWarning
  use dftbp_message, only : error, warning
  use dftbp_solventdata, only : TSolventData
  use dftbp_tokenreader, only : getNextToken, TOKEN_OK
  use dftbp_xmlf90, only : fnode
  implicit none
  private

  public :: readParamGBSA


  !> Read GBSA parametrisation data
  interface readParamGBSA
    module procedure :: readParamGBSAFile
  end interface readParamGBSA


contains


  !> Read GBSA parameters from file
  subroutine readParamGBSAFile(file, input, solvent, speciesNames, node)

    !> Name of the parametrisation file
    character(len=*), intent(in) :: file

    !> Contains the input for the dispersion module on exit
    type(TGBInput), intent(out) :: input

    !> Experimental data for solvent
    type(TSolventData), intent(out) :: solvent

    !> Symbols of all species
    character(len=*), intent(in) :: speciesNames(:)

    !> Node for error handling
    type(fnode), pointer, optional :: node

    integer, parameter :: nParam = 8
    integer, parameter :: nElem = 94
    integer :: ii, lineno, unit, iStart, iErr, iSp, iZp, nSpecies
    real(dp) :: param(nParam)
    real(dp) :: descreening(nElem), surfaceTension(nElem), hBondPar(nElem)
    character(len=lc) :: line, errorStr

    open(file=file, newunit=unit, status='old', iostat=iErr, iomsg=errorStr)
    if (iErr /= 0) then
      if (present(node)) then
        call detailedError(node, "Could not open '"//trim(file)//"': "//trim(errorStr))
      else
        call error("Could not open '"//trim(file)//"': "//trim(errorStr))
      end if
    end if

    lineno = 0

    do ii = 1, 8
      call nextLine(unit, line, lineno, file, node=node)
      iStart = 1
      call getNextToken(trim(line), param(ii), iStart, iErr)
      if (iErr /= TOKEN_OK) then
        iStart = min(verify(line(iStart:), whiteSpaces) + iStart - 1, len_trim(line))
        write(errorStr, '(a, "(", i0, "):", 1x, 6a)') &
          & trim(file), lineno, "Could not read real", newline, trim(line), &
          & newline, repeat('-', max(iStart-1, 0)), '^'
        if (present(node)) then
          call detailedError(node, trim(errorStr))
        else
          call error(trim(errorStr))
        end if
      end if
      if (iStart < len_trim(line)) then
        write(errorStr, '(a, "(", i0, "):", 1x, 6a)') &
          & trim(file), lineno, "Trailing content", newline, trim(line), newline, &
          & repeat('-', max(iStart-1, 0)), '^'
        if (present(node)) then
          call detailedWarning(node, trim(errorStr))
        else
          call warning(trim(errorStr))
        end if
      end if
    end do

    do ii = 1, nElem
      call nextLine(unit, line, lineno, file, node=node)
      iStart = 1
      call getNextToken(trim(line), surfaceTension(ii), iStart, iErr)
      if (iErr /= TOKEN_OK) then
        iStart = min(verify(line(iStart:), whiteSpaces) + iStart - 1, len_trim(line))
        write(errorStr, '(a, "(", i0, "):", 1x, 6a)') trim(file), lineno, &
          & "Could not read surface tension", newline, trim(line), newline, &
          & repeat('-', max(iStart-1, 0)), '^'
        if (present(node)) then
          call detailedError(node, trim(errorStr))
        else
          call error(trim(errorStr))
        end if
      end if
      call getNextToken(trim(line), descreening(ii), iStart, iErr)
      if (iErr /= TOKEN_OK) then
        iStart = min(verify(line(iStart:), whiteSpaces) + iStart - 1, len_trim(line))
        write(errorStr, '(a, "(", i0, "):", 1x, 6a)') trim(file), lineno, &
          & "Could not read descreening", newline, trim(line), newline, &
          & repeat('-', max(iStart-1, 0)), '^'
        if (present(node)) then
          call detailedError(node, trim(errorStr))
        else
          call error(trim(errorStr))
        end if
      end if
      call getNextToken(trim(line), hBondPar(ii), iStart, iErr)
      if (iErr /= TOKEN_OK) then
        iStart = min(verify(line(iStart:), whiteSpaces) + iStart - 1, len_trim(line))
        write(errorStr, '(a, "(", i0, "):", 1x, 6a)') trim(file), lineno, &
          & "Could not read hydrogen bond strength", newline, trim(line), newline, &
          & repeat('-', max(iStart-1, 0)), '^'
        if (present(node)) then
          call detailedError(node, trim(errorStr))
        else
          call error(trim(errorStr))
        end if
      end if
      if (iStart < len_trim(line)) then
        write(errorStr, '(a, "(", i0, "):", 1x, 6a)') &
          & trim(file), lineno, "Trailing content", newline, trim(line), newline, &
          & repeat('-', max(iStart-1, 0)), '^'
        if (present(node)) then
          call detailedWarning(node, trim(errorStr))
        else
          call warning(trim(errorStr))
        end if
      end if
    end do

    allocate(input%sasaInput)

    solvent%dielectricConstant = param(1)
    solvent%molecularMass = param(2) * amu__au
    solvent%density = param(3) * 1.0e+3_dp*kg__au/(1.0e10_dp*AA__Bohr)**3
    input%bornScale = param(4)
    input%sasaInput%probeRad = param(5) * AA__Bohr
    input%freeEnergyShift = param(6) * kcal_mol__Hartree
    input%bornOffset = param(7) * AA__Bohr * 0.1_dp

    nSpecies = size(speciesNames, dim=1)
    allocate(input%descreening(nSpecies))
    allocate(input%sasaInput%surfaceTension(nSpecies))
    allocate(input%hBondPar(nSpecies))

    do iSp = 1, nSpecies
      iZp = symbolToNumber(speciesNames(iSp))
      if (iZp > 0 .and. iZp <= nElem) then
        input%descreening(iSp) = descreening(iZp)
        input%sasaInput%surfaceTension(iSp) = surfaceTension(iZp)
        input%hBondPar(iSp) = -hBondPar(iZp)**2 * kcal_mol__Hartree
      else
        write(errorStr, '(3a)') trim(file), " contains no parameters for species ", &
            & trim(speciesNames(iSp))
        if (present(node)) then
          call detailedWarning(node, trim(errorStr))
        else
          call warning(trim(errorStr))
        end if
        input%descreening(iSp) = 1.0_dp
        input%sasaInput%surfaceTension(iSp) = 0.0_dp
        input%hBondPar(iSp) = 0.0_dp
      end if
    end do

    if (all(input%hBondPar > -1.0e-14_dp)) then
      deallocate(input%hBondPar)
    end if

    close(unit)

  end subroutine readParamGBSAFile


  !> Read a whole line from a formatted IO unit
  subroutine nextLine(unit, line, lineno, file, iostat, node)

    !> IO-unit bound to the parametrisation file
    integer, intent(in) :: unit

    !> Line buffer
    character(len=*), intent(inout) :: line

    !> Current line number, will be incremented after successful read
    integer, intent(inout) :: lineno

    !> Name of the parametrisation file
    character(len=*), intent(in), optional :: file

    !> Node for error handling
    type(fnode), pointer, optional :: node

    !> Error code
    integer, intent(out), optional :: iostat

    integer :: length, iErr
    character(len=lc) :: errorStr, iomsg

    lineno = lineno + 1

    read(unit, '(a)', advance='no', iostat=iErr, iomsg=iomsg, size=length) line
    if (length >= len(line)) then
      if (present(file)) then
        write(errorStr, '(a, "(", i0, "):", 1x, a)') &
            & trim(file), lineno, "too many characters, line truncated"
      else
        write(errorStr, '(a, 1x, i0, ":", 1x, a)') &
            & "line", lineno, "too many characters, line truncated"
      end if
      if (present(node)) then
        call detailedWarning(node, trim(errorStr))
      else
        call warning(trim(errorStr))
      end if
      ! drop the rest of the line
      if (.not.is_iostat_eor(iErr)) then
        read(unit, '(a)', iostat=iErr, iomsg=iomsg)
      end if
    end if

    ! end-of-record is an implementation detail
    if (is_iostat_eor(iErr)) then
      iErr = 0
    end if

    if (present(iostat)) then
      iostat = iErr
    else
      if (iErr /= 0) then
        if (is_iostat_end(iErr)) then
          ! uncaught EOF is an error
          if (present(file)) then
            write(errorStr, '(a, "(", i0, "):", 1x, a)') &
              & trim(file), lineno, "encountered end-of-file while reading"
          else
            write(errorStr, '(a, 1x, i0, ":", 1x, a)') &
              & "line", lineno, "encountered end-of-file while reading"
          end if
          if (present(node)) then
            call detailedError(node, trim(errorStr))
          else
            call error(trim(errorStr))
          end if
        else
          ! unknown error, lets hope for something useful in iomsg
          if (present(node)) then
            call detailedError(node, trim(iomsg))
          else
            call error(trim(iomsg))
          end if
        end if
      end if
    end if

  end subroutine nextLine


end module dftbp_gbsafile
