!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains the routines for initialising modes.
module InitProgram
  use assert
  use HSDParser, only : parseHSD, dumpHSD, dumpHSDAsXML
  use XMLUtils
  use HSDUtils
  use HSDUtils2
  use flib_dom
  use linkedlist
  use CharManip
  use Accuracy
  use Constants
  use TypeGeometryHSD
  use Message
  use FileId
  use UnitConversion
  use OldSKData
  implicit none

  private
  save

  character(len=*), parameter :: version =  "0.01"
  character(len=*), parameter :: rootTag = "modes"
  character(len=*), parameter :: hsdInput = "modes_in.hsd"
  character(len=*), parameter :: hsdParsedInput = "modes_pin.hsd"
  character(len=*), parameter :: xmlInput = "modes_in.xml"
  character(len=*), parameter :: xmlParsedInput = "modes_pin.xml"
  integer, parameter :: parserVersion = 3

  public :: initProgramVariables, destructProgramVariables


  !! Variables from detailed.xml
  integer, public :: identity          ! Identity of the run
  type(TGeometry), public :: geo       ! Geometry

  !! Variables from the Option block
  logical, public :: tVerbose          ! If program should be verbose

  real(dp), allocatable, public :: atomicMasses(:)
  real(dp), allocatable, public :: dynMatrix(:,:)

  logical, public :: tPlotModes
  logical, public :: tAnimateModes
  logical, public :: tXmakeMol
  integer, allocatable, public :: modesToPlot(:)
  integer, public :: nModesToPlot
  integer, public :: nCycles
  integer, public, parameter :: nSteps = 10
  integer, public :: nMovedAtom  ! Number of atoms which should be moved.
  integer, allocatable, public :: iMovedAtoms(:)

  !! Locally created variables

contains

  !!* Initialise program variables
  subroutine initProgramVariables()

    type(TOldSKData) :: skData
    type(fnode), pointer :: input, root, node, tmp
    type(fnode), pointer :: value, child, child2
    type(listRealR1) :: realBuffer
    type(string) :: buffer, buffer2
    type(listString) :: lStr
    integer :: inputVersion
    integer :: ii, iSp1, iAt
    logical :: tHSD
    real(dp), allocatable :: speciesMass(:)
    type(listCharLc), allocatable :: skFiles(:)
    character(lc) :: prefix, suffix, separator, elem1, strTmp, filename
    logical :: tLower, tExist
    integer :: nDerivs
    logical :: tWriteXML, tWriteHSD ! XML or HSD output?

    !! Write header
    write (*, "(A)") repeat("=", 80)
    write (*, "(A)") "     MODES  " // version
    write (*, "(A,/)") repeat("=", 80)

    !! Read in input file as HSD or XML.
    call readHSDOrXML(hsdInput, xmlInput, rootTag, input, tHSD)
    if (tHSD) then
      write (*, "(A)") "Interpreting input file '" // hsdInput // "'"
    else
      write (*, "(A)") "Interpreting input file '" // xmlInput //  "'"
    end if
    write (*, "(A)") repeat("-", 80)
    call getChild(input, rootTag, root)

    !! Check if input version is the one, which we can handle
    call getChildValue(root, "InputVersion", inputVersion, parserVersion)
    if (inputVersion /= parserVersion) then
      call error("Version of input (" // i2c(inputVersion) // ") and parser (" &
          &// i2c(parserVersion) // ") do not match")
    end if


    call getChild(root, "Geometry", tmp)
    call readGeometry(tmp, geo)

    call getChildValue(root, "Atoms", buffer2, "1:-1", child=child, multiple=.true.)
    call convAtomRangeToInt(char(buffer2), geo%speciesNames, geo%species, &
        &child, iMovedAtoms)
    nMovedAtom = size(iMovedAtoms)
    nDerivs = 3 * nMovedAtom

    call getChild(root, "DisplayModes",child=node,requested=.false.)
    if (associated(node)) then
      tPlotModes = .true.
      call getChildValue(node, "PlotModes", buffer2, "1:-1", child=child, &
          &multiple=.true.)
      call convRangeToInt(char(buffer2), child, modesToPlot, 3 * nMovedAtom)
      nModesToPlot = size(modesToPlot)
      call getChildValue(node, "Animate", tAnimateModes, .true.)
      call getChildValue(node, "XMakeMol", tXmakeMol, .true.)
    else
      nModesToPlot = 0
      tPlotModes = .false.
      tAnimateModes = .false.
      tXmakeMol = .false.
    end if

    if (tAnimateModes.and.tXmakeMol) then
      nCycles = 1
    else
      nCycles = 3
    end if

    !! Slater-Koster files
    allocate(skFiles(geo%nSpecies))
    do iSp1 = 1, geo%nSpecies
        call init(skFiles(iSp1))
    end do

    call getChildValue(root, "SlaterKosterFiles", value, child=child)
    call getNodeName(value, buffer)
    select case(char(buffer))
    case ("type2filenames")
      call getChildValue(value, "Prefix", buffer2, "")
      prefix = unquote(char(buffer2))
      call getChildValue(value, "Suffix", buffer2, "")
      suffix = unquote(char(buffer2))
      call getChildValue(value, "Separator", buffer2, "")
      separator = unquote(char(buffer2))
      call getChildValue(value, "LowerCaseTypeName", tLower, .false.)
      do iSp1 = 1, geo%nSpecies
        if (tLower) then
          elem1 = tolower(geo%speciesNames(iSp1))
        else
          elem1 = geo%speciesNames(iSp1)
        end if
        strTmp = trim(prefix) // trim(elem1) // trim(separator) &
            &// trim(elem1) // trim(suffix)
        call append(skFiles(iSp1), strTmp)
        inquire(file=strTmp, exist=tExist)
        if (.not. tExist) then
          call detailedError(value, "SK file with generated name '" &
              &// trim(strTmp) // "' does not exist.")
        end if
      end do
    case default
      call setUnprocessed(value)
      do iSp1 = 1, geo%nSpecies
        strTmp = trim(geo%speciesNames(iSp1)) // "-" &
            &// trim(geo%speciesNames(iSp1))
        call init(lStr)
        call getChildValue(child, trim(strTmp), lStr, child=child2)
        ! We can't handle selected shells here (also not needed I guess)
        if (len(lStr) /= 1) then
          call detailedError(child2, "Incorrect number of Slater-Koster &
              &files")
        end if
        do ii = 1, len(lStr)
          call get(lStr, strTmp, ii)
          inquire(file=strTmp, exist=tExist)
          if (.not. tExist) then
            call detailedError(child2, "SK file '" // trim(strTmp) &
                &// "' does not exist'")
          end if
          call append(skFiles(iSp1), strTmp)
        end do
        call destruct(lStr)
      end do
    end select

    allocate(speciesMass(geo%nSpecies))
    do iSp1 = 1, geo%nSpecies
      call get(skFiles(iSp1), fileName, 1)
      call readFromFile(skData, fileName, .true.)
      speciesMass(iSp1) = skData%mass
      call destruct(skFiles(iSp1))
    end do

    allocate(dynMatrix(nDerivs,nDerivs))
    call getChildValue(root, "Hessian", value, "", child=child, &
        & allowEmptyValue=.true.)
    call getNodeName2(value, buffer)
    if (char(buffer) == "") then
      call error("No derivative matrix supplied!")
    else
      call init(realBuffer)
      call getChildValue(child, "", nDerivs, realBuffer)
      if (len(realBuffer)/=nDerivs) then
        call detailedError(root,"wrong number of derivatives supplied:" &
            & // i2c(len(realBuffer)) // " supplied, " &
            & // i2c(nDerivs) // " required.")
      end if
      call asArray(realBuffer, dynMatrix)
      call destruct(realBuffer)
    end if

    call getChildValue(root, "WriteHSDInput", tWriteHSD, .true.)
    call getChildValue(root, "WriteXMLInput", tWriteXML, .false.)

    !! Issue warning about unprocessed nodes
    call warnUnprocessedNodes(root,.true.)

    !! Finish parsing, dump parsed and processed input
    if (tWriteHSD) then
      call dumpHSD(input, hsdParsedInput)

      write (*, "(A)") "Processed input written as HSD to '" // hsdParsedInput &
          &//"'"
    end if
    if (tWriteXML) then
      call dumpHSDAsXML(input, xmlParsedInput)
      write (*, "(A)") "Processed input written as XML to '" // xmlParsedInput &
          &//"'"
    end if
    write (*, "(A)") repeat("-", 80)
    write (*,*)
    call destroyNode(input)


    allocate(atomicMasses(nMovedAtom))
    do iAt = 1, nMovedAtom
      atomicMasses(iAt) = speciesMass(geo%species(iMovedAtoms(iAt)))
    end do

  end subroutine initProgramVariables

  !!* Destroy the program variables created in initProgramVariables
  subroutine destructProgramVariables()

    write (*, "(/,A)") repeat("=", 80)

  end subroutine destructProgramVariables

  !!* Read in the geometry stored as xml in internal or gen format.
  !!* @param geonode Node containing the geometry
  !!* @param geo     Contains the geometry information on exit
  subroutine readGeometry(geonode, geo)
    type(fnode), pointer :: geonode
    type(TGeometry), intent(out) :: geo

    type(fnode), pointer :: child
    type(string) :: buffer

    call getChildValue(geonode, "", child)
    call getNodeName(child, buffer)
    select case (char(buffer))
    case ("genformat")
      call readTGeometryGen(child, geo)
      call removeChildNodes(geonode)
      call writeTGeometryHSD(geonode, geo)
    case default
      call readTGeometryHSD(geonode, geo)
    end select

  end subroutine readGeometry

end module InitProgram
