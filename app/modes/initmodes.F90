!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the routines for initialising modes.
module modes_initmodes
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_atomicmass, only : getAtomicMass
  use dftbp_common_file, only : TFileDescr, openFile, closeFile
  use dftbp_common_filesystem, only : findFile, joinPathsPrettyErr, getParamSearchPaths
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_release, only : releaseYear
  use dftbp_common_unitconversion, only : massUnits
  use dftbp_extlibs_xmlf90, only : fnode, fNodeList, string, char, getLength, getItem1,&
      & getNodeName, destroyNode, destroyNodeList, textNodeName
  use dftbp_io_charmanip, only : i2c, newline, tolower, unquote
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
  implicit none

  private
  public :: initProgramVariables
  public :: geo, atomicMasses, dynMatrix, bornMatrix, bornDerivsMatrix, modesToPlot, nModesToPlot
  public :: nCycles, nSteps, nMovedAtom, iMovedAtoms, nDerivs
  public :: tVerbose, tPlotModes, tEigenVectors, tAnimateModes, tRemoveTranslate, tRemoveRotate
  public :: setEigvecGauge


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
  subroutine initProgramVariables()

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
    character(lc) :: prefix, suffix, separator, elem1, strTmp, str2Tmp, filename
    logical :: tLower, tDumpPHSD
    logical :: tWriteHSD
    type(string), allocatable :: searchPath(:)
    character(len=:), allocatable :: strOut, strJoin, hessianFile
    type(TFileDescr) :: file
    integer :: iErr

    !! Write header
    call printDftbHeader('(MODES '// version //')', releaseYear)

    !! Read in input file as HSD
    call parseHSD(rootTag, hsdInput, hsdTree)
    call getChild(hsdTree, rootTag, root)

    write(stdout, "(A)") "Interpreting input file '" // hsdInput // "'"
    write(stdout, "(A)") repeat("-", 80)

    !! Check if input version is the one, which we can handle
    call getChildValue(root, "InputVersion", inputVersion, parserVersion)
    if (inputVersion /= parserVersion) then
      call error("Version of input (" // i2c(inputVersion) // ") and parser ("&
          & // i2c(parserVersion) // ") do not match")
    end if

    call getChild(root, "Geometry", tmp)
    call readGeometry(tmp, geo)

    call getChildValue(root, "RemoveTranslation", tRemoveTranslate, .false.)
    call getChildValue(root, "RemoveRotation", tRemoveRotate, .false.)

    call getChildValue(root, "Atoms", buffer2, "1:-1", child=child, multiple=.true.)
    call getSelectedAtomIndices(child, char(buffer2), geo%speciesNames, geo%species, iMovedAtoms)
    nMovedAtom = size(iMovedAtoms)
    nDerivs = 3 * nMovedAtom

    tPlotModes = .false.
    nModesToPlot = 0
    tAnimateModes = .false.
    call getChild(root, "DisplayModes",child=node,requested=.false.)
    if (associated(node)) then
      tPlotModes = .true.
      call getChildValue(node, "PlotModes", buffer2, "1:-1", child=child, multiple=.true.)
      call getSelectedIndices(child, char(buffer2), [1, 3 * nMovedAtom], modesToPlot)
      nModesToPlot = size(modesToPlot)
      call getChildValue(node, "Animate", tAnimateModes, .true.)
    end if

    ! oscillation cycles in an animation
    nCycles = 3

    ! Slater-Koster files
    call getParamSearchPaths(searchPath)
    strJoin = joinPathsPrettyErr(searchPath)
    allocate(speciesMass(geo%nSpecies))
    speciesMass(:) = 0.0_dp
    do iSp1 = 1, geo%nSpecies
      speciesMass(iSp1) = getAtomicMass(geo%speciesNames(iSp1))
    end do

    call getChildValue(root, "SlaterKosterFiles", value, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true.)
    if (associated(value)) then
      allocate(skFiles(geo%nSpecies))
      do iSp1 = 1, geo%nSpecies
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
        do iSp1 = 1, geo%nSpecies
          if (tLower) then
            elem1 = tolower(geo%speciesNames(iSp1))
          else
            elem1 = geo%speciesNames(iSp1)
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

        do iSp1 = 1, geo%nSpecies
          strTmp = trim(geo%speciesNames(iSp1)) // "-" // trim(geo%speciesNames(iSp1))
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
      do iSp1 = 1, geo%nSpecies
        call get(skFiles(iSp1), fileName, 1)
        call readFromFile(skData, fileName, .true.)
        speciesMass(iSp1) = skData%mass
        call destruct(skFiles(iSp1))
      end do
    end if

    call getInputMasses(root, geo, replacementMasses)
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

    call getChildValue(root, "Hessian", value, "", child=child, allowEmptyValue=.true.)
    call getNodeName2(value, buffer)
    select case (char(buffer))
    case ("directread")
      call getChildValue(value, "File", buffer2, child=child2)
      hessianFile = trim(unquote(char(buffer2)))
      call openFile(file, hessianFile, mode="r", iostat=iErr)
      if (iErr /= 0) then
        call detailedError(child2, "Could not open file '" // hessianFile&
            & // "' for direct reading." )
      end if
      read(file%unit, *, iostat=iErr) dynMatrix
      if (iErr /= 0) then
        call detailedError(child2, "Error during direct reading '" // hessianFile // "'.")
      end if
      call closeFile(file)
    case (textNodeName)
      call getNodeName2(value, buffer)
      call init(realBufferList)
      call getChildValue(child, "", nDerivs, realBufferList)
      if (len(realBufferList) /= nDerivs) then
        call detailedError(root,"wrong number of derivatives supplied:"&
            & // i2c(len(realBufferList)) // " supplied, " // i2c(nDerivs) // " required.")
      end if
      call asArray(realBufferList, dynMatrix)
      call destruct(realBufferList)
      tDumpPHSD = .false.
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
            & // i2c(len(realBuffer)) // " supplied, " // i2c(3*nDerivs) // " required.")
      end if
      allocate(bornMatrix(len(realBuffer)))
      call asArray(realBuffer, bornMatrix)
      call destruct(realBuffer)
      tDumpPHSD = .false.
    end if

    call getChild(root, "BornDerivs", child, requested=.false.)
    call getNodeName2(child, buffer)
    if (char(buffer) /= "") then
      call init(realBuffer)
      call getChildValue(child, "", realBuffer)
      if (len(realBuffer) /= 9 * nDerivs) then
        call detailedError(root,"wrong number of Born charge derivatives supplied:"&
            & // i2c(len(realBuffer)) // " supplied, " // i2c(9 * nDerivs) // " required.")
      end if
      allocate(bornDerivsMatrix(len(realBuffer)))
      call asArray(realBuffer, bornDerivsMatrix)
      call destruct(realBuffer)
      tDumpPHSD = .false.
    end if

    call getChildValue(root, "WriteHSDInput", tWriteHSD, tDumpPHSD)

    tEigenVectors = tPlotModes .or. allocated(bornMatrix) .or. allocated(bornDerivsMatrix)

    !! Issue warning about unprocessed nodes
    call warnUnprocessedNodes(root, .true.)

    !! Finish parsing, dump parsed and processed input
    if (tWriteHSD) then
      call dumpHSD(hsdTree, hsdParsedInput)
      write(stdout, "(A)") "Processed input written as HSD to '" // hsdParsedInput // "'"
    end if
    write(stdout, "(A)") repeat("-", 80)
    write(stdout, *)
    call destroyNode(hsdTree)

  end subroutine initProgramVariables


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
      call getSelectedAtomIndices(child3, char(buffer), geo%speciesNames, geo%species, pTmpI1)
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


  !> Returns gauge-corrected eigenvectors, such that the first non-zero coefficient of each mode is
  !! positive.
  subroutine setEigvecGauge(eigvec)

    !> Gauge corrected eigenvectors on exit. Shape: [iCoeff, iMode]
    real(dp), intent(inout) :: eigvec(:,:)

    !! Auxiliary variables
    integer :: iMode, iCoeff

    do iMode = 1, size(eigvec, dim=2)
      lpCoeff: do iCoeff = 1, size(eigvec, dim=1)
        if (abs(eigvec(iCoeff, iMode)) > 1e2_dp * epsilon(1.0_dp)) then
          if (sign(1.0_dp, eigvec(iCoeff, iMode)) < 0.0_dp) then
            eigvec(:, iMode) = -eigvec(:, iMode)
          end if
          exit lpCoeff
        end if
      end do lpCoeff
    end do

  end subroutine setEigvecGauge

end module modes_initmodes
