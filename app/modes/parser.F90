!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Fills the derived type with the input parameters from an HSD or an XML file.
module modes_parser
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_atomicmass, only : getAtomicMass
  use dftbp_common_filesystem, only : findFile, joinPathsPrettyErr, getParamSearchPaths
  use dftbp_common_globalenv, only : stdOut, withMpi, withScalapack
  use dftbp_common_release, only : releaseYear
  use dftbp_common_unitconversion, only : massUnits
  use dftbp_extlibs_xmlf90, only : fnode, fNodeList, string, char, getLength, getItem1,&
      & getNodeName, destroyNode, destroyNodeList, textNodeName
  use dftbp_io_charmanip, only : i2c, newline, tolower, unquote, newline
  use dftbp_io_formatout, only : printDftbHeader
  use dftbp_io_hsdparser, only : parseHSD, dumpHSD
  use dftbp_io_hsdutils, only : getChild, getChildValue, getChildren, getSelectedAtomIndices,&
      & getSelectedIndices, detailedError, detailedWarning
  use dftbp_io_hsdutils2, only : convertUnitHsd, setUnprocessed, warnUnprocessedNodes, getNodeName2
  use dftbp_io_message, only : error
  use dftbp_io_xmlutils, only : removeChildNodes
  use dftbp_type_linkedlist, only : TListCharLc, TListReal, TListRealR1, TListString, init,&
      & destruct, append, get, len, asArray
  use dftbp_type_oldskdata, only : TOldSkData, readFromFile
  use dftbp_type_typegeometryhsd, only : TGeometry, readTGeometryGen, readTGeometryXyz,&
      & readTGeometryHsd, readTGeometryVasp, writeTGeometryHsd
  use dftbp_common_status, only : TStatus
  use modes_inputdata, only : TInputData, TControl, TBlacsOpts, solverTypes
  use dftbp_io_hsdparser, only : parseHsd
  use dftbp_extlibs_xmlf90, only : fnode
  use dftbp_dftbplus_input_fileaccess, only : readBinaryAccessTypes
  implicit none

  private
  public :: parserVersion, rootTag
  public :: TParserFlags
  public :: readHsdFile, parseHsdTree


  !> Tag at the head of the input document tree
  character(len=*), parameter :: rootTag = "modesinput"

  !> Version of the current parser
  integer, parameter :: parserVersion = 3


  !> Container type for parser related flags.
  type TParserFlags

    !> Stop after parsing?
    logical :: tStop

    !> Continue despite unprocessed nodes
    logical :: tIgnoreUnprocessed

    !> HSD output?
    logical :: tWriteHSD

  end type TParserFlags


contains

  !> Reads the HSD input from a file.
  subroutine readHsdFile(hsdFile, hsdTree)

    !> Name of the input file
    character(*), intent(in) :: hsdFile

    !> Data tree representation of the input
    type(fnode), pointer :: hsdTree

    call parseHSD(rootTag, hsdFile, hsdTree)

  end subroutine readHsdFile


  !> Parse input from an HSD/XML file.
  subroutine parseHsdTree(hsdTree, input, parserFlags)

    !> Tree representation of the input
    type(fnode), pointer :: hsdTree

    !> Returns initialised input variables on exit
    type(TInputData), intent(out) :: input

    !> Special block containings parser related settings
    type(TParserFlags), intent(out) :: parserFlags

    type(TOldSKData) :: skData
    type(fnode), pointer :: root, node, tmp
    type(fnode), pointer :: value, child, child2, dummy
    type(TListRealR1) :: realBufferList
    type(TListReal) :: realBuffer
    type(string) :: buffer, buffer2
    type(TListString) :: lStr
    integer :: inputVersion
    integer :: ii, iSp1, iAt, nMovedAtom, iErr, nDerivs
    real(dp), allocatable :: speciesMass(:), replacementMasses(:)
    type(TListCharLc), allocatable :: skFiles(:)
    character(lc) :: prefix, suffix, separator, elem1, strTmp, str2Tmp, filename
    logical :: tLower, tDumpPHSD
    type(string), allocatable :: searchPath(:)
    character(len=:), allocatable :: strOut, strJoin, hessianFile

    call getChild(hsdTree, rootTag, root)

    ! Check if input version is the one, which we can handle
    call getChildValue(root, "InputVersion", inputVersion, parserVersion)
    if (inputVersion /= parserVersion) then
      call error("Version of input (" // i2c(inputVersion) // ") and parser ("&
          & // i2c(parserVersion) // ") do not match")
    end if

    call getChildValue(root, "ParserOptions", dummy, "", child=child, list=.true.,&
        & allowEmptyValue=.true., dummyValue=.true.)
    call readParserOptions(child, parserFlags)

    call getChild(root, "Geometry", tmp)
    call readGeometry(tmp, input%geo)

    call getChildValue(root, "RemoveTranslation", input%ctrl%tRemoveTranslate, .false.)
    call getChildValue(root, "RemoveRotation", input%ctrl%tRemoveRotate, .false.)

    call getChildValue(root, "Atoms", buffer2, "1:-1", child=child, multiple=.true.)
    call getSelectedAtomIndices(child, char(buffer2), input%geo%speciesNames, input%geo%species,&
        & input%ctrl%iMovedAtoms)
    nMovedAtom = size(input%ctrl%iMovedAtoms)
    nDerivs = 3 * nMovedAtom

    input%ctrl%tPlotModes = .false.
    input%ctrl%tAnimateModes = .false.
    call getChild(root, "DisplayModes", child=node, requested=.false.)
    if (associated(node)) then
      input%ctrl%tPlotModes = .true.
      call getChildValue(node, "PlotModes", buffer2, "1:-1", child=child, multiple=.true.)
      call getSelectedIndices(child, char(buffer2), [1, nDerivs], input%ctrl%modesToPlot)
      call getChildValue(node, "Animate", input%ctrl%tAnimateModes, .true.)
    end if

    ! Eigensolver
    call getChildValue(root, "EigenSolver", buffer2, "qr")
    select case(tolower(char(buffer2)))
    case ("qr")
      input%ctrl%iSolver = solverTypes%qr
    case ("divideandconquer")
      input%ctrl%iSolver = solverTypes%divideAndConquer
    case ("relativelyrobust")
      input%ctrl%iSolver = solverTypes%relativelyRobust
    case ("magma")
    #:if WITH_MAGMA
      call TGpuEnv_init(gpu)
    #:else
      call error("Magma-solver selected, but program was compiled without MAGMA")
    #:endif
      input%ctrl%iSolver = solverTypes%magmaEvd
    case default
      call detailedError(root, "Unknown eigensolver " // char(buffer2))
    end select

    ! Slater-Koster files
    call getParamSearchPaths(searchPath)
    strJoin = joinPathsPrettyErr(searchPath)
    allocate(speciesMass(input%geo%nSpecies), source=0.0_dp)
    do iSp1 = 1, input%geo%nSpecies
      speciesMass(iSp1) = getAtomicMass(input%geo%speciesNames(iSp1))
    end do

    call getChildValue(root, "SlaterKosterFiles", value, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true.)
    if (associated(value)) then
      allocate(skFiles(input%geo%nSpecies))
      do iSp1 = 1, input%geo%nSpecies
        call init(skFiles(iSp1))
      end do
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
        do iSp1 = 1, input%geo%nSpecies
          if (tLower) then
            elem1 = tolower(input%geo%speciesNames(iSp1))
          else
            elem1 = input%geo%speciesNames(iSp1)
          end if
          strTmp = trim(prefix) // trim(elem1) // trim(separator) // trim(elem1) // trim(suffix)
          call findFile(searchPath, strTmp, strOut)
          if (.not. allocated(strOut)) then
            call detailedError(value, "SK file with generated name '" // trim(strTmp)&
                & // "' not found." // newline // "   (search path(s): " // strJoin // ").")
          end if
          strTmp = strOut
          call append(skFiles(iSp1), strTmp)
        end do
      case default
        call setUnprocessed(value)
        call getChildValue(child, "Prefix", buffer2, "")
        prefix = unquote(char(buffer2))
        do iSp1 = 1, input%geo%nSpecies
          strTmp = trim(input%geo%speciesNames(iSp1)) // "-" // trim(input%geo%speciesNames(iSp1))
          call init(lStr)
          call getChildValue(child, trim(strTmp), lStr, child=child2)
          ! We can't handle selected shells here (also not needed)
          if (len(lStr) /= 1) then
            call detailedError(child2, "Incorrect number of Slater-Koster files")
          end if
          do ii = 1, len(lStr)
            call get(lStr, str2Tmp, ii)
            strTmp = trim(prefix) // str2Tmp
            call findFile(searchPath, strTmp, strOut)
            if (.not. allocated(strOut)) then
              call detailedError(child2, "SK file '" // trim(strTmp) // "' not found." // newline&
                  & // "   (search path(s): " // strJoin // ").")
            end if
            strTmp = strOut
            call append(skFiles(iSp1), strTmp)
          end do
          call destruct(lStr)
        end do
      end select
      do iSp1 = 1, input%geo%nSpecies
        call get(skFiles(iSp1), fileName, 1)
        call readFromFile(skData, fileName, .true.)
        speciesMass(iSp1) = skData%mass
        call destruct(skFiles(iSp1))
      end do
    end if

    call getInputMasses(root, input%geo, replacementMasses)
    allocate(input%ctrl%atomicMasses(nMovedAtom))
    do iAt = 1, nMovedAtom
      input%ctrl%atomicMasses(iAt) = speciesMass(input%geo%species(input%ctrl%iMovedAtoms(iAt)))
    end do
    if (allocated(replacementMasses)) then
      do iAt = 1, nMovedAtom
        if (replacementMasses(input%ctrl%iMovedAtoms(iAt)) >= 0.0_dp) then
          input%ctrl%atomicMasses(iAt) = replacementMasses(input%ctrl%iMovedAtoms(iAt))
        end if
      end do
    end if

    call getChildValue(root, "Hessian", value, "", child=child, allowEmptyValue=.true.)
    call getNodeName2(value, buffer)
    select case (char(buffer))
    case ("directread")
      call getChildValue(value, "File", buffer2, child=child2)
      input%ctrl%hessianFile = trim(unquote(char(buffer2)))
    case (textNodeName)
      if (withMpi) then
        call detailedError(root, "For MPI enabled builds only the 'DirectRead' method is&
            & supported.")
      end if
      call getNodeName2(value, buffer)
      call init(realBufferList)
      call getChildValue(child, "", nDerivs, realBufferList)
      if (len(realBufferList) /= nDerivs) then
        call detailedError(root, "wrong number of derivatives supplied:"&
            & // i2c(len(realBufferList)) // " supplied, " // i2c(nDerivs)&
            & // " required.")
      end if
      allocate(input%ctrl%hessian(nDerivs, nDerivs))
      call asArray(realBufferList, input%ctrl%hessian)
      call destruct(realBufferList)
    case default
      call detailedError(child, "Invalid Hessian scheme.")
    end select

    call getChild(root, "BornCharges", child, requested=.false.)
    call getNodeName2(child, buffer)
    if (char(buffer) /= "") then
      call init(realBuffer)
      call getChildValue(child, "", realBuffer)
      if (len(realBuffer) /= 3 * nDerivs) then
        call detailedError(root,"wrong number of Born charges supplied:"&
            & // i2c(len(realBuffer)) // " supplied, " // i2c(3 * nDerivs) // " required.")
      end if
      allocate(input%ctrl%bornMatrix(len(realBuffer)))
      call asArray(realBuffer, input%ctrl%bornMatrix)
      call destruct(realBuffer)
    end if

    call getChild(root, "BornDerivs", child, requested=.false.)
    call getNodeName2(child, buffer)
    if (char(buffer) /= "") then
      call init(realBuffer)
      call getChildValue(child, "", realBuffer)
      if (len(realBuffer) /= 9 * nDerivs) then
        call detailedError(root, "wrong number of Born charge derivatives supplied:"&
            & // i2c(len(realBuffer)) // " supplied, " // i2c(9 * nDerivs) // " required.")
      end if
      allocate(input%ctrl%bornDerivsMatrix(len(realBuffer)))
      call asArray(realBuffer, input%ctrl%bornDerivsMatrix)
      call destruct(realBuffer)
    end if

    call getChildValue(root, "Options", dummy, "", child=child, list=.true.,&
        & allowEmptyValue=.true., dummyValue=.true.)
    call readOptions(child, input%ctrl)

    ! read parallel calculation settings
    call readParallel(root, input%ctrl)

    ! input data strucutre has been initialised
    input%tInitialized = .true.

  end subroutine parseHsdTree


  !> Read in parser options (options not passed to the main code).
  subroutine readParserOptions(node, flags)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Contains parser flags on exit.
    type(TParserFlags), intent(out) :: flags

    call getChildValue(node, "WriteHSDInput", flags%tWriteHSD, .true.)
    if (.not. flags%tWriteHSD) then
      call detailedWarning(node, "WriteHSDInput turned off. You are not guaranteed" // newline // &
          &" to able to obtain the same results with a later version of the code!" // newline // &
          & "(the modes_pin.hsd file DOES guarantee this)")
    end if
    call getChildValue(node, "StopAfterParsing", flags%tStop, .false.)

    call getChildValue(node, "IgnoreUnprocessedNodes", flags%tIgnoreUnprocessed, .true.)

  end subroutine readParserOptions


  !> Reads the Parallel block.
  subroutine readParallel(root, ctrl)

    !> Root node eventually containing the current block
    type(fnode), intent(in), pointer :: root

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    type(fnode), pointer :: node

    call getChild(root, "Parallel", child=node, requested=.false., emptyIfMissing=withMpi)
    if (associated(node)) then
      if (.not. withMpi) then
        call detailedWarning(node, "Settings will be read but ignored (compiled without MPI&
            & support)")
      end if
      allocate(ctrl%parallelOpts)
      ! call getChildValue(node, "Groups", ctrl%parallelOpts%nGroup, 1, child=pTmp)
      ctrl%parallelOpts%nGroup = 1
      ! call getChildValue(node, "UseOmpThreads", ctrl%parallelOpts%tOmpThreads, .not. withMpi)
      ctrl%parallelOpts%tOmpThreads = .false.
      call readBlacs(node, ctrl%parallelOpts%blacsOpts)
    end if

  end subroutine readParallel


  !> Reads the Blacs block.
  subroutine readBlacs(root, blacsOpts)

    !> Root node eventually containing the current block
    type(fnode), pointer, intent(in) :: root

    !> Blacs settings
    type(TBlacsOpts), intent(inout) :: blacsOpts

    type(fnode), pointer :: node

    call getChild(root, "Blacs", child=node, requested=.false., emptyIfMissing=withScalapack)
    if (associated(node)) then
      if (.not. withScalapack) then
        call detailedWarning(node, "Settings will be read but ignored (compiled without SCALAPACK&
            & support)")
      end if
      call getChildValue(node, "BlockSize", blacsOpts%blockSize, 32)
    end if

  end subroutine readBlacs


  !> Reads the Option block.
  subroutine readOptions(node, ctrl)

    !> Node to parse
    type(fnode), pointer :: node

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    call readBinaryAccessTypes(node, ctrl%binaryAccessTypes)

  end subroutine readOptions


  !> Read in the geometry stored as gen format.
  subroutine readGeometry(geonode, geo)

    !> Node containing the geometry
    type(fnode), pointer :: geonode

    !> Contains the geometry information on exit
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
    case ("xyzformat")
      call readTGeometryXyz(child, geo)
      call removeChildNodes(geonode)
      call writeTGeometryHSD(geonode, geo)
    case ("vaspformat")
      call readTGeometryVasp(child, geo)
      call removeChildNodes(geonode)
      call writeTGeometryHSD(geonode, geo)
    case default
      call readTGeometryHSD(geonode, geo)
    end select

  end subroutine readGeometry


  !> Reads atomic masses from input file.
  subroutine getInputMasses(node, geo, masses)

    !> Relevant node of input data
    type(fnode), pointer :: node

    !> Geometry object, which contains atomic species information
    type(TGeometry), intent(in) :: geo

    !> Masses to be returned
    real(dp), allocatable, intent(out) :: masses(:)

    type(fnode), pointer :: child, child2, child3, val
    type(fnodeList), pointer :: children
    integer, allocatable :: pTmpI1(:)
    type(string) :: buffer, modifier
    real(dp) :: rTmp
    integer :: ii, jj, iAt

    call getChildValue(node, "Masses", val, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)

    ! Read individual atom specifications
    call getChildren(child, "Mass", children)
    if (getLength(children) == 0) then
      return
    end if

    allocate(masses(geo%nAtom))
    masses(:) = -1.0_dp
    do ii = 1, getLength(children)
      call getItem1(children, ii, child2)
      call getChildValue(child2, "Atoms", buffer, child=child3, multiple=.true.)
      call getSelectedAtomIndices(child3, char(buffer), geo%speciesNames, geo%species,&
          & pTmpI1)
      call getChildValue(child2, "MassPerAtom", rTmp, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), massUnits, child, rTmp)
      do jj = 1, size(pTmpI1)
        iAt = pTmpI1(jj)
        if (masses(iAt) >= 0.0_dp) then
          call detailedWarning(child3, "Previous setting for the mass  of atom" // i2c(iAt)&
              & // " overwritten")
        end if
        masses(iAt) = rTmp
      end do
      deallocate(pTmpI1)
    end do
    call destroyNodeList(children)

  end subroutine getInputMasses

end module modes_parser
