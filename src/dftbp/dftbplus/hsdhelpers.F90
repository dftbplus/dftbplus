!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"
#:include "error.fypp"

!> HSD-parsing related helper routines.
module dftbp_dftbplus_hsdhelpers
  use dftbp_common_globalenv, only : stdOut, tIoProc
  use dftbp_common_status, only : TStatus
  use dftbp_dftbplus_inputdata, only : TInputData
  use dftbp_dftbplus_parser, only : TParserFlags, rootTag, readHsdFile, parseHsdTree
  use dftbp_extlibs_xmlf90, only : fnode, destroyNode
  use dftbp_io_hsdparser, only : dumpHSD
  use dftbp_io_hsdutils, only : getChild
  use dftbp_io_hsdutils2, only : warnUnprocessedNodes
  use dftbp_io_message, only : error
  implicit none

  private
  public :: parseHsdInput, doPostParseJobs


  !> Name of the DFTB+ input file
  character(*), parameter :: hsdFileName = "dftb_in.hsd"

  !> Name of the DFTB+ processed input file
  character(*), parameter :: hsdProcFileName = "dftb_pin.hsd"

contains

  !> Parses input file and returns initialised input structure
  subroutine parseHsdInput(input, errStatus)

    !> Input data parsed from the input file
    type(TInputData), intent(out) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: hsdTree
    type(TParserFlags) :: parserFlags

    call ensureInputFilePresence(errStatus)
    @:PROPAGATE_ERROR(errStatus)
    write(stdout, "(A)") "Reading input file '" // hsdFileName // "'"
    call readHsdFile(hsdFileName, hsdTree, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call parseHsdTree(hsdTree, input, parserFlags, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call doPostParseJobs(hsdTree, parserFlags, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call destroyNode(hsdTree)

  end subroutine parseHsdInput


  !> Checks whether input file is present and stops if not
  subroutine ensureInputFilePresence(errStatus)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    logical :: tExist

    inquire(file=hsdFileName, exist=tExist)
    if (.not. tExist) then
      @:RAISE_ERROR(errStatus, -1, "No input file '" // hsdFileName // "' not found.")
    end if

  end subroutine ensureInputFilePresence


  !> Execute parser related tasks (warning, processed input dumping) needed after parsing
  subroutine doPostParseJobs(hsdTree, parserFlags, errStatus)

    !> Tree representation of the HSD input
    type(fnode), pointer, intent(in) :: hsdTree

    !> Parser specific settings in the output
    type(TParserFlags), intent(in) :: parserFlags

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: root

    call getChild(hsdTree, rootTag, root, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! Issue warning about unprocessed nodes
    call warnUnprocessedNodes(root, errStatus, tIgnoreUnprocessed=parserFlags%tIgnoreUnprocessed)
    @:PROPAGATE_ERROR(errStatus)

    ! Dump processed tree in HSD and XML format
    if (tIoProc .and. parserFlags%tWriteHSD) then
      call dumpHSD(hsdTree, hsdProcFileName, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      write(stdout, '(/,/,A)') "Processed input in HSD format written to '" // hsdProcFileName&
          & // "'"
    end if

    ! Stop, if only parsing is required
    if (parserFlags%tStop) then
      @:RAISE_ERROR(errStatus, -1, "Keyword 'StopAfterParsing' is set to Yes. Stopping.")
    end if

  end subroutine doPostParseJobs


end module dftbp_dftbplus_hsdhelpers
