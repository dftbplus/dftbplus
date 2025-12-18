!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Read GBSA parametrisation data from file
module dftbp_solvation_gbsafile
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_constants, only : AA__Bohr, amu__au, kcal_mol__Hartree, kg__au, symbolToNumber
  use dftbp_common_file, only : closeFile, openFile, TFileDescr
  use dftbp_extlibs_xmlf90, only : fnode
  use dftbp_io_charmanip, only : newline, whiteSpaces
  use dftbp_io_hsdutils, only : detailedError, detailedWarning
  use dftbp_io_linereader, only : TLineReader
  use dftbp_io_message, only : error, warning
  use dftbp_io_tokenreader, only : getNextToken, TOKEN_OK
  use dftbp_solvation_born, only : TGBInput
  use dftbp_solvation_solventdata, only : TSolventData
  implicit none

  private
  public :: readParamGBSA


contains


  !> Read GBSA parameters from file
  subroutine readParamGBSA(file, input, solvent, speciesNames, node)

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

    type(TFileDescr) :: fd
    type(TLineReader) :: lineReader
    integer, parameter :: nParam = 8
    integer, parameter :: nElem = 94
    integer :: ii, lineno, iStart, iErr, iSp, iZp, nSpecies
    real(dp) :: param(nParam)
    real(dp) :: descreening(nElem), surfaceTension(nElem), hBondPar(nElem)
    character(len=lc) :: errorStr
    character(:), allocatable :: line, ioMsg

    call openFile(fd, file, mode="r", ioStat=iErr, ioMsg=ioMsg)
    if (iErr /= 0) then
      if (present(node)) then
        call detailedError(node, "Could not open '"//trim(file)//"': "//trim(ioMsg))
      else
        call error("Could not open '"//trim(file)//"': "//trim(ioMsg))
      end if
    end if
    lineReader = TLineReader(fd%unit)

    lineno = 0

    do ii = 1, 8
      call nextLine(lineReader, line, lineno, file, node=node)
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
      call nextLine(lineReader, line, lineno, file, node=node)
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

    call closeFile(fd)

  end subroutine readParamGBSA


  !> Read a whole line from a formatted IO unit
  subroutine nextLine(lineReader, line, lineno, file, node)

    !> Line reader able to provide the next line from an open file
    type(TLineReader), intent(inout) :: lineReader

    !> Line buffer
    character(:), allocatable, intent(out) :: line

    !> Current line number, will be incremented after successful read
    integer, intent(inout) :: lineno

    !> Name of the parametrisation file
    character(len=*), intent(in) :: file

    !> Node (needed for generating error messages)
    type(fnode), pointer, intent(in) :: node

    integer :: iErr
    character(len=lc) :: errorStr, iomsg

    lineno = lineno + 1
    call lineReader%readLine(line, iErr)

    if (iErr /= 0) then
      write(errorStr, '(3a, i0)') "While reading file ", trim(file),&
          & "an error was encountered on line ", lineno
      call detailedError(node, trim(errorStr))
    end if

  end subroutine nextLine

end module dftbp_solvation_gbsafile
