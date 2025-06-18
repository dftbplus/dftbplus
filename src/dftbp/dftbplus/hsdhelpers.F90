!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

!> HSD-parsing related helper routines.
module dftbp_dftbplus_hsdhelpers
  use dftbp_common_globalenv, only : stdOut, tIoProc
  use dftbp_common_exception, only : TException
  use dftbp_dftbplus_inputdata, only : TInputData
  use dftbp_dftbplus_parser, only : parseHsdTree, readHsdFile, rootTag, TParserFlags
  use dftbp_extlibs_xmlf90, only : destroyNode, fnode
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
  subroutine parseHsdInput(exc, input)

    !> Exception
    type(TException), allocatable, intent(out) :: exc

    !> Input data parsed from the input file
    type(TInputData), intent(out) :: input

    type(fnode), pointer :: hsdTree
    type(TParserFlags) :: parserFlags

    call ensureInputFilePresence()
    write(stdout, "(A)") "Reading input file '" // hsdFileName // "'"
    call readHsdFile(hsdFileName, hsdTree)
    call parseHsdTree(hsdTree, input, parserFlags)
    call doPostParseJobs(hsdTree, parserFlags)
    call destroyNode(hsdTree)

  end subroutine parseHsdInput


  !> Checks whether input file is present and stops if not
  subroutine ensureInputFilePresence()

    logical :: tExist

    inquire(file=hsdFileName, exist=tExist)
    if (.not. tExist) then
      call error("No input file '" // hsdFileName // "' not found.")
    end if

  end subroutine ensureInputFilePresence


  !> Execute parser related tasks (warning, processed input dumping) needed after parsing
  subroutine doPostParseJobs(hsdTree, parserFlags)

    !> Tree representation of the HSD input
    type(fnode), pointer, intent(in) :: hsdTree

    !> Parser specific settings in the output
    type(TParserFlags), intent(in) :: parserFlags

    type(fnode), pointer :: root

    call getChild(hsdTree, rootTag, root)

    ! Issue warning about unprocessed nodes
    call warnUnprocessedNodes(root, parserFlags%tIgnoreUnprocessed)

    ! Dump processed tree in HSD and XML format
    if (tIoProc .and. parserFlags%tWriteHSD) then
      call dumpHSD(hsdTree, hsdProcFileName)
      write(stdout, '(/,/,A)') "Processed input in HSD format written to '" // hsdProcFileName&
          & // "'"
    end if

    ! Stop, if only parsing is required
    if (parserFlags%tStop) then
      call error("Keyword 'StopAfterParsing' is set to Yes. Stopping.")
    end if

  end subroutine doPostParseJobs


end module dftbp_dftbplus_hsdhelpers
