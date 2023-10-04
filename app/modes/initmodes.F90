!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Contains the routines for initialising modes.
module modes_initmodes
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_atomicmass, only : getAtomicMass
  use dftbp_common_filesystem, only : findFile, getParamSearchPath
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_release, only : releaseYear
  use dftbp_common_status, only : TStatus
  use dftbp_common_unitconversion, only : massUnits
  use dftbp_extlibs_xmlf90, only : fnode, fNodeList, string, char, getLength, getItem1,&
      & getNodeName, destroyNode, destroyNodeList
  use dftbp_io_charmanip, only : i2c, tolower, unquote
  use dftbp_io_formatout, only : printDftbHeader
  use dftbp_io_hsdparser, only : parseHSD, dumpHSD
  use dftbp_io_hsdutils, only : getChild, getChildValue, getChildren, getSelectedAtomIndices,&
      & getSelectedIndices, detailedError, detailedWarning
  use dftbp_io_hsdutils2, only : convertUnitHsd, setUnprocessed, warnUnprocessedNodes, getNodeName2
  use dftbp_io_xmlutils, only : removeChildNodes
  use dftbp_type_linkedlist, only : TListCharLc, TListReal, TListRealR1, TListString, init,&
      & destruct, append, get, len, asArray
  use dftbp_type_oldskdata, only : TOldSkData, readFromFile
  use dftbp_type_typegeometryhsd, only : TGeometry, readTGeometryGen, readTGeometryXyz,&
      & readTGeometryHsd, readTGeometryVasp, writeTGeometryHsd
  implicit none

  private
  public :: initProgramVariables
  public :: geo, atomicMasses, dynMatrix, bornMatrix, bornDerivsMatrix, modesToPlot, nModesToPlot
  public :: nCycles, nSteps, nMovedAtom, iMovedAtoms, nDerivs
  public :: tVerbose, tPlotModes, tEigenVectors, tAnimateModes, tRemoveTranslate, tRemoveRotate


  !> Program version
  character(len=*), parameter :: version = "0.03"

  !> Root node name of the input tree
  character(len=*), parameter :: rootTag = "modes"

  !> Input file name
  character(len=*), parameter :: hsdInput = "modes_in.hsd"

  !> Parsed output name
  character(len=*), parameter :: hsdParsedInput = "modes_pin.hsd"

  !> Version of the input document
  integer, parameter :: parserVersion = 3

  !> Geometry
  type(TGeometry) :: geo

  !> If program should be verbose
  logical :: tVerbose


  !> Atomic masses to build dynamical matrix
  real(dp), allocatable :: atomicMasses(:)

  !> Dynamical matrix
  real(dp), allocatable :: dynMatrix(:,:)

  !> Born charges matrix
  real(dp), allocatable :: bornMatrix(:)

  !> Derivatives of Born charges matrix with respect to electric field, i.e. polarizability
  !> derivatives with respect to atom locations
  real(dp), allocatable :: bornDerivsMatrix(:)

  !> Produce plots of modes
  logical :: tPlotModes

  !> Produce eigenvectors of modes, either for plotting or for property changes along mode
  !> directions
  logical :: tEigenVectors

  !> Animate mode  or as vectors
  logical :: tAnimateModes

  !> Remove translation modes
  logical :: tRemoveTranslate

  !> Remove rotation modes
  logical :: tRemoveRotate

  !> Modes to produce xyz file for
  integer, allocatable :: modesToPlot(:)

  !> Number of modes being plotted
  integer :: nModesToPlot

  !> If animating, number of cycles to show in an animation
  integer :: nCycles

  !> Steps in an animation cycle
  integer, parameter :: nSteps = 10

  !> Number of atoms which should be moved.
  integer :: nMovedAtom

  !> List of atoms in dynamical matrix
  integer, allocatable :: iMovedAtoms(:)

  !> Number of derivatives
  integer :: nDerivs

contains


  !> Initialise program variables.
  subroutine initProgramVariables(errStatus)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(TOldSKData) :: skData
    type(fnode), pointer :: root, node, tmp, hsdTree
    type(fnode), pointer :: value, child, child2
    type(TListRealR1) :: realBufferList
    type(TListReal) :: realBuffer
    type(string) :: buffer, buffer2
    type(TListString) :: lStr
    integer :: inputVersion
    integer :: ii, iSp1, iAt
    real(dp), allocatable :: speciesMass(:), replacementMasses(:)
    type(TListCharLc), allocatable :: skFiles(:)
    character(lc) :: prefix, suffix, separator, elem1, strTmp, filename
    logical :: tLower, tExist, tDumpPHSD
    logical :: tWriteHSD
    type(string), allocatable :: searchPath(:)
    character(len=:), allocatable :: strOut

    !! Write header
    call printDftbHeader('(MODES '// version //')', releaseYear)

    !! Read in input file as HSD
    call parseHSD(rootTag, hsdInput, hsdTree, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChild(hsdTree, rootTag, root, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    write(stdout, "(A)") "Interpreting input file '" // hsdInput // "'"
    write(stdout, "(A)") repeat("-", 80)

    !! Check if input version is the one, which we can handle
    call getChildValue(root, "InputVersion", inputVersion, errStatus, parserVersion)
    @:PROPAGATE_ERROR(errStatus)
    if (inputVersion /= parserVersion) then
      @:RAISE_ERROR(errStatus, -1, "Version of input (" // i2c(inputVersion) // ") and parser ("&
          & // i2c(parserVersion) // ") do not match")
    end if

    call getChild(root, "Geometry", tmp, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call readGeometry(tmp, geo, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(root, "RemoveTranslation", tRemoveTranslate, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(root, "RemoveRotation", tRemoveRotate, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(root, "Atoms", buffer2, errStatus, "1:-1", child=child, multiple=.true.)
    @:PROPAGATE_ERROR(errStatus)
    call getSelectedAtomIndices(child, char(buffer2), geo%speciesNames, geo%species, iMovedAtoms,&
        & errStatus)
    @:PROPAGATE_ERROR(errStatus)
    nMovedAtom = size(iMovedAtoms)
    nDerivs = 3 * nMovedAtom

    call getChild(root, "DisplayModes", node, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(node)) then
      tPlotModes = .true.
      call getChildValue(node, "PlotModes", buffer2, errStatus, "1:-1", child=child,&
          & multiple=.true.)
      @:PROPAGATE_ERROR(errStatus)
      call getSelectedIndices(child, char(buffer2), [1, 3 * nMovedAtom], modesToPlot, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      nModesToPlot = size(modesToPlot)
      call getChildValue(node, "Animate", tAnimateModes, errStatus, .true.)
      @:PROPAGATE_ERROR(errStatus)
    else
      nModesToPlot = 0
      tPlotModes = .false.
      tAnimateModes = .false.
    end if

    ! oscillation cycles in an animation
    nCycles = 3

    ! Slater-Koster files
    call getParamSearchPath(searchPath)
    allocate(speciesMass(geo%nSpecies))
    speciesMass(:) = 0.0_dp
    do iSp1 = 1, geo%nSpecies
      speciesMass(iSp1) = getAtomicMass(geo%speciesNames(iSp1))
    end do

    call getChildValue(root, "SlaterKosterFiles", value, errStatus, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(value)) then
      allocate(skFiles(geo%nSpecies))
      do iSp1 = 1, geo%nSpecies
        call init(skFiles(iSp1))
      end do
      call getNodeName(value, buffer)
      select case(char(buffer))
      case ("type2filenames")
        call getChildValue(value, "Prefix", buffer2, errStatus, "")
        @:PROPAGATE_ERROR(errStatus)
        prefix = unquote(char(buffer2))
        call getChildValue(value, "Suffix", buffer2, errStatus, "")
        @:PROPAGATE_ERROR(errStatus)
        suffix = unquote(char(buffer2))
        call getChildValue(value, "Separator", buffer2, errStatus, "")
        @:PROPAGATE_ERROR(errStatus)
        separator = unquote(char(buffer2))
        call getChildValue(value, "LowerCaseTypeName", tLower, errStatus, .false.)
        @:PROPAGATE_ERROR(errStatus)
        do iSp1 = 1, geo%nSpecies
          if (tLower) then
            elem1 = tolower(geo%speciesNames(iSp1))
          else
            elem1 = geo%speciesNames(iSp1)
          end if
          strTmp = trim(prefix) // trim(elem1) // trim(separator) // trim(elem1) // trim(suffix)
          call findFile(searchPath, strTmp, strOut)
          if (allocated(strOut)) strTmp = strOut
          call append(skFiles(iSp1), strTmp)
          inquire(file=strTmp, exist=tExist)
          if (.not. tExist) then
            call detailedError(value, "SK file with generated name '" // trim(strTmp)&
                & // "' does not exist.", errStatus)
            @:PROPAGATE_ERROR(errStatus)
          end if
        end do
      case default
        call setUnprocessed(value)
        do iSp1 = 1, geo%nSpecies
          strTmp = trim(geo%speciesNames(iSp1)) // "-" // trim(geo%speciesNames(iSp1))
          call init(lStr)
          call getChildValue(child, trim(strTmp), lStr, errStatus, child=child2)
          @:PROPAGATE_ERROR(errStatus)
          ! We can't handle selected shells here (also not needed)
          if (len(lStr) /= 1) then
            call detailedError(child2, "Incorrect number of Slater-Koster files", errStatus)
            @:PROPAGATE_ERROR(errStatus)
          end if
          do ii = 1, len(lStr)
            call get(lStr, strTmp, ii)
            inquire(file=strTmp, exist=tExist)
            if (.not. tExist) then
              call detailedError(child2, "SK file '" // trim(strTmp) // "' does not exist'",&
                  & errStatus)
              @:PROPAGATE_ERROR(errStatus)
            end if
            call append(skFiles(iSp1), strTmp)
          end do
          call destruct(lStr)
        end do
      end select
      do iSp1 = 1, geo%nSpecies
        call get(skFiles(iSp1), fileName, 1)
        call readFromFile(skData, fileName, .true.)
        speciesMass(iSp1) = skData%mass
        call destruct(skFiles(iSp1))
      end do
    end if

    call getInputMasses(root, geo, replacementMasses, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    allocate(atomicMasses(nMovedAtom))
    do iAt = 1, nMovedAtom
      atomicMasses(iAt) = speciesMass(geo%species(iMovedAtoms(iAt)))
    end do
    if (allocated(replacementMasses)) then
      do iAt = 1, nMovedAtom
        if (replacementMasses(iMovedAtoms(iAt)) >= 0.0_dp) then
          atomicMasses(iAt) = replacementMasses(iMovedAtoms(iAt))
        end if
      end do
    end if

    allocate(dynMatrix(nDerivs, nDerivs))

    tDumpPHSD = .true.

    call getChildValue(root, "Hessian", value, errStatus, "", child=child, allowEmptyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName2(value, buffer)
    if (char(buffer) == "") then
      @:RAISE_ERROR(errStatus, -1, "No derivative matrix supplied!")
    else
      call init(realBufferList)
      call getChildValue(child, "", nDerivs, realBufferList, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (len(realBufferList)/=nDerivs) then
        call detailedError(root,"wrong number of derivatives supplied:"&
            & // i2c(len(realBufferList)) // " supplied, "&
            & // i2c(nDerivs) // " required.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call asArray(realBufferList, dynMatrix)
      call destruct(realBufferList)
      tDumpPHSD = .false.
    end if

    call getChild(root, "BornCharges", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName2(child, buffer)
    if (char(buffer) /= "") then
      call init(realBuffer)
      call getChildValue(child, "", realBuffer, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (len(realBuffer) /= 3 * nDerivs) then
        call detailedError(root, "wrong number of Born charges supplied:"&
            & // i2c(len(realBuffer)) // " supplied, " // i2c(3*nDerivs) // " required.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      allocate(bornMatrix(len(realBuffer)))
      call asArray(realBuffer, bornMatrix)
      call destruct(realBuffer)
      tDumpPHSD = .false.
    end if

    call getChild(root, "BornDerivs", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName2(child, buffer)
    if (char(buffer) /= "") then
      call init(realBuffer)
      call getChildValue(child, "", realBuffer, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (len(realBuffer) /=9 * nDerivs) then
        call detailedError(root, "wrong number of Born charge derivatives supplied:"&
            & // i2c(len(realBuffer)) // " supplied, " // i2c(9 * nDerivs) // " required.",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      allocate(bornDerivsMatrix(len(realBuffer)))
      call asArray(realBuffer, bornDerivsMatrix)
      call destruct(realBuffer)
      tDumpPHSD = .false.
    end if

    call getChildValue(root, "WriteHSDInput", tWriteHSD, errStatus, tDumpPHSD)
    @:PROPAGATE_ERROR(errStatus)

    tEigenVectors = tPlotModes .or. allocated(bornMatrix) .or. allocated(bornDerivsMatrix)

    !! Issue warning about unprocessed nodes
    call warnUnprocessedNodes(root, errStatus, tIgnoreUnprocessed=.true.)
    @:PROPAGATE_ERROR(errStatus)

    !! Finish parsing, dump parsed and processed input
    if (tWriteHSD) then
      call dumpHSD(hsdTree, hsdParsedInput, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      write(stdout, "(A)") "Processed input written as HSD to '" // hsdParsedInput // "'"
    end if
    write(stdout, "(A)") repeat("-", 80)
    write(stdout, *)
    call destroyNode(hsdTree)

  end subroutine initProgramVariables


  !> Read in the geometry stored as gen format.
  subroutine readGeometry(geonode, geo, errStatus)

    !> Node containing the geometry
    type(fnode), pointer :: geonode

    !> Contains the geometry information on exit
    type(TGeometry), intent(out) :: geo

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child
    type(string) :: buffer

    call getChildValue(geonode, "", child, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(child, buffer)
    select case (char(buffer))
    case ("genformat")
      call readTGeometryGen(child, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call removeChildNodes(geonode)
      call writeTGeometryHSD(geonode, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("xyzformat")
      call readTGeometryXyz(child, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call removeChildNodes(geonode)
      call writeTGeometryHSD(geonode, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("vaspformat")
      call readTGeometryVasp(child, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call removeChildNodes(geonode)
      call writeTGeometryHSD(geonode, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case default
      call readTGeometryHSD(geonode, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine readGeometry


  !> Reads atomic masses from input file.
  subroutine getInputMasses(node, geo, masses, errStatus)

    !> Relevant node of input data
    type(fnode), pointer :: node

    !> Geometry object, which contains atomic species information
    type(TGeometry), intent(in) :: geo

    !> Masses to be returned
    real(dp), allocatable, intent(out) :: masses(:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child, child2, child3, val
    type(fnodeList), pointer :: children
    integer, allocatable :: pTmpI1(:)
    type(string) :: buffer, modifier
    real(dp) :: rTmp
    integer :: ii, jj, iAt

    call getChildValue(node, "Masses", val, errStatus, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)
    @:PROPAGATE_ERROR(errStatus)

    ! Read individual atom specifications
    call getChildren(child, "Mass", children)
    if (getLength(children) == 0) then
      return
    end if

    allocate(masses(geo%nAtom))
    masses(:) = -1.0_dp
    do ii = 1, getLength(children)
      call getItem1(children, ii, child2)
      call getChildValue(child2, "Atoms", buffer, errStatus, child=child3, multiple=.true.)
      @:PROPAGATE_ERROR(errStatus)
      call getSelectedAtomIndices(child3, char(buffer), geo%speciesNames, geo%species, pTmpI1,&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(child2, "MassPerAtom", rTmp, errStatus, modifier=modifier, child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), massUnits, child, rTmp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
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

end module modes_initmodes
