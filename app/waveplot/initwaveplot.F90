!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the routines for initialising waveplot.
module waveplot_initwaveplot

  use dftbp_common_accuracy, only : dp
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_unitconversion, only : lengthUnits
  use dftbp_extlibs_xmlf90, only : fnode, fNodeList, string, char, getLength, getItem1,&
      & getNodeName,destroyNode
  use dftbp_io_charmanip, only : i2c, unquote
  use dftbp_io_fileid, only : getFileId
  use dftbp_io_hsdparser, only : parseHSD, dumpHSD
  use dftbp_io_hsdutils, only : getChildValue, setChildValue, getChild, setChild, getChildren,&
      & getSelectedIndices, detailedError, detailedWarning
  use dftbp_io_hsdutils2, only : getModifierIndex, readHSDAsXML, warnUnprocessedNodes
  use dftbp_io_message, only : warning, error
  use dftbp_io_xmlutils, only : removeChildNodes
  use dftbp_math_simplealgebra, only : invert33
  use dftbp_type_linkedlist, only : TListIntR1, TListReal, init, destruct, len, append, asArray
  use dftbp_type_typegeometryhsd, only : TGeometry, readTGeometryGen, readTGeometryHSD,&
      & readTGeometryVasp, readTGeometryXyz, writeTGeometryHSD

  ! use waveplot_gridcache, only : TGridCache, init
  use waveplot_molorb, only : TMolecularOrbital, init, TSpeciesBasis
  use waveplot_slater, only : init


  implicit none

  private
  save

  public :: TProgramVariables
  public :: TProgramVariables_init


  !> Data type containing variables from detailed.xml
  type TXml

    !> Geometry
    type(TGeometry) :: geo

    !> Identity of the run
    integer :: identity

    !> Nr. of orbitals per state
    integer :: nOrb

    !> If eigvecs/hamiltonian is real
    logical :: tRealHam

    !> Occupations
    real(dp), allocatable :: occupations(:,:,:)

  end type TXml


  !> Data type containing variables from eigenvec.bin
  type TEig

    !> Nr. of states
    integer :: nState

    !> Real eigenvectors
    real(dp), allocatable :: eigvecsReal(:,:)

    !> Complex eigenvectors
    complex(dp), allocatable :: eigvecsCplx(:,:)

  end type TEig


  !> Data type containing variables from the Option block
  type TOption

    !> Number of (total) grid points along 3 directions
    integer :: nTotPoints(3)

    !> Number of (species) grid points along 3 directions
    integer :: nSpPoints(3)

    !> Repeat box along 3 directions
    integer :: repeatBox(3)

    !> Levels to plot
    integer, allocatable :: plottedLevels(:)

    !> K-points to plot
    integer, allocatable :: plottedKPoints(:)

    !> Spins to plot
    integer, allocatable :: plottedSpins(:)

    !> If box should filled with folded atoms
    logical :: tFillBox

    !> If coords should be folded to unit cell
    logical :: tFoldCoords

    !> If program should be verbose
    logical :: tVerbose

    !> If total charge should be plotted
    logical :: tPlotTotChrg

    !> If total charge should be calculated
    logical :: tCalcTotChrg

    !> If total spin pol. to be plotted
    logical :: tPlotTotSpin

    !> If total charge difference to be plotted
    logical :: tPlotTotDiff

    !> If atomic densities to be plotted
    logical :: tPlotAtomDens

    !> If atomic densities to be calculated
    logical :: tCalcAtomDens

    !> If charge for orbitals to be plotted
    logical :: tPlotChrg

    !> If charge difference for orbs. to be plotted
    logical :: tPlotChrgDiff

    !> If real part of the wfcs to plot.
    logical :: tPlotReal

    !> If imaginary part of the wfcs to plot
    logical :: tPlotImag

    !> Box vectors for the plotted region
    real(dp) :: boxVecs(3,3)

    !> Origin of the box
    real(dp) :: origin(3)

    !> Origin of the (total) grid in the box
    real(dp) :: totGridOrig(3)

    !> List of levels to plot, whereby insignificant occupations were filtered out
    integer, allocatable :: levelIndex(:,:)

    ! !> Origin of the (species) grids in the box
    ! real(dp), allocatable :: spGridOrigs(:,:)

  end type TOption


  !> Data type containing variables from the Basis block
  type TBasis

    !> definition of the wfcs
    type(TSpeciesBasis), allocatable :: basis(:)

  end type TBasis


  !> Data type containing variables from the AtomicNumbers block
  type TAtomicNumber

    !> species-atomic nr. corresp.
    integer, allocatable :: atomicNumbers(:)

  end type TAtomicNumber


  !> Data type containing locally created variables
  type TInternal

    !> Molecular orbital
    type(TMolecularOrbital), allocatable :: molOrb

    !> pointer to the orbital
    type(TMolecularOrbital), pointer :: pMolOrb

    !> (total) grid vectors
    real(dp) :: totGridVec(3,3)

    !> (species) grids vectors
    real(dp), allocatable :: speciesGridsVecs(:,:,:)

    !> Volume of the grid
    real(dp) :: gridVol

    !> Center coordinate of total grid
    real(dp) :: totGridCenter(3)

    !> m-resolved orbital occupations
    real(dp), allocatable :: orbitalOcc(:,:)

    !> Origins of species grids
    real(dp), allocatable :: speciesGridsOrigs(:,:)

    !> Index mappin orbital --> atom
    integer, allocatable :: orbitalToAtom(:)

    !> Index mappin orbital --> species
    integer, allocatable :: orbitalToSpecies(:)

  end type TInternal


  !> Data type containing program variables
  type TProgramVariables

    !> Data of detailed.xml
    type(TXml) :: xml

    !> Data of eigenvec.bin
    type(TEig) :: eig

    !> Data of Option block
    type(TOption) :: option

    !> Data of Basis block
    type(TBasis) :: basis

    !> Data of AtomicNumber block
    type(TAtomicNumber) :: atomicNumber

    !> Locally created data
    type(TInternal) :: internal

  end type TProgramVariables


  !> program version
  character(len=*), parameter :: version =  "0.3"

  !> root node name of the input tree
  character(len=*), parameter :: rootTag = "waveplot"

  !> input file name
  character(len=*), parameter :: hsdInput = "waveplot_in.hsd"

  !> parsed output name
  character(len=*), parameter :: hsdParsedInput = "waveplot_pin.hsd"

  !> version of the input document
  integer, parameter :: parserVersion = 3


  interface readEigenvecs
    module procedure readRealEigenvecs
    module procedure readCplxEigenvecs
  end interface readEigenvecs


contains


  !> Initialise program variables
  subroutine TProgramVariables_init(this)

    !> Container of program variables
    type(TProgramVariables), intent(inout), target :: this

    type(fnode), pointer :: root, tmp, detailed, hsdTree

    !> File with binary eigenvectors
    character(len=1024) :: eigVecBin

    !> String buffer instance
    type(string) :: strBuffer

    !> Nr. of K-Points
    integer :: nKPoint

    !> Nr. of spins
    integer :: nSpin

    integer :: inputVersion

    !> If grid should shifted by a half cell
    logical :: tShiftGrid

    logical :: tGroundState

    !! Write header
    write(stdout, "(A)") repeat("=", 80)
    write(stdout, "(A)") "     WAVEPLOT  " // version
    write(stdout, "(A,/)") repeat("=", 80)

    !! Read in input file as HSD
    call parseHSD(rootTag, hsdInput, hsdTree)
    call getChild(hsdTree, rootTag, root)

    write(stdout, "(A)") "Interpreting input file '" // hsdInput // "'"

    !! Check if input version is the one, which we can handle
    call getChildValue(root, "InputVersion", inputVersion, parserVersion)
    if (inputVersion /= parserVersion) then
      call error("Version of input (" // i2c(inputVersion) // ") and parser (" &
          &// i2c(parserVersion) // ") do not match")
    end if

    call getChildValue(root, "GroundState", tGroundState, .true.)

    !! Read data from detailed.xml
    call getChildValue(root, "DetailedXML", strBuffer)
    call readHSDAsXML(unquote(char(strBuffer)), tmp)
    call getChild(tmp, "detailedout", detailed)
    call readDetailed(this, detailed, tGroundState, nKPoint, nSpin, this%eig%nState)
    call destroyNode(tmp)

    !! Read basis
    call getChild(root, "Basis", tmp)
    call readBasis(this, tmp,  this%xml%geo%speciesNames)
    call getChildValue(root, "EigenvecBin", strBuffer)
    eigVecBin = unquote(char(strBuffer))
    call checkEigenvecs(eigVecBin, this%xml%identity)

    !! Read options
    call getChild(root, "Options", tmp)
    call readOptions(this, tmp, this%eig%nState, nKPoint, nSpin, tShiftGrid)

    !! Read eigenvectors
    if (this%xml%tRealHam) then
      allocate(this%eig%eigvecsReal(this%xml%nOrb, this%eig%nState))
      call readEigenvecs(eigVecBin, this%eig%eigvecsReal)
    else
      allocate(this%eig%eigvecsCplx(this%xml%nOrb, this%eig%nState))
      call readEigenvecs(eigVecBin, this%eig%eigvecsCplx)
    end if

    !! Issue warning about unprocessed nodes
    call warnUnprocessedNodes(root, .true.)

    !! Finish parsing, dump parsed and processed input
    call dumpHSD(hsdTree, hsdParsedInput)

    write(stdout, "(A)") "Processed input written as HSD to '" // hsdParsedInput //"'"
    write(stdout, "(A)") repeat("-", 80)
    write(stdout, *)

    call destroyNode(hsdTree)

    write(stdout, "(A)") "Doing initialisation"

    !! Initialize necessary objects
    allocate(this%internal%molOrb)
    this%internal%pMolOrb => this%internal%molOrb
    call init(this%internal%molOrb, this%xml%geo, this%basis%basis)

    !! Determine useful index mappings
    call getIndexMappings(this)

    !! Repeat boxes if necessary
    call getRepeatedBox(this)

  end subroutine TProgramVariables_init


  !> Interpret the information stored in detailed.xml
  subroutine readDetailed(this, detailed, tGroundState, nKPoint, nSpin, nState)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> Pointer to the node, containing the info
    type(fnode), intent(in), pointer :: detailed

    !> look for ground state occupations (True) or excited (False)
    logical, intent(in) :: tGroundState

    !> Nr. of K-Points
    integer, intent(out) :: nKPoint

    !> Nr. of spins
    integer, intent(out) :: nSpin

    !> Nr. of states
    integer, intent(out) :: nState

    type(fnode), pointer :: tmp, occ, spin

    !> K-Points & weights
    real(dp), allocatable :: kPointsWeights(:,:)

    integer :: iSpin, iK

    call getChildValue(detailed, "Identity", this%xml%identity)
    call getChild(detailed, "Geometry", tmp)
    call readGeometry(this, tmp)

    call getChildValue(detailed, "Real", this%xml%tRealHam)
    call getChildValue(detailed, "NrOfKPoints", nKPoint)
    call getChildValue(detailed, "NrOfSpins", nSpin)
    call getChildValue(detailed, "NrOfStates", nState)
    call getChildValue(detailed, "NrOfOrbitals", this%xml%nOrb)

    allocate(kPointsWeights(4, nKPoint))
    allocate(this%xml%occupations(nState, nKPoint, nSpin))

    call getChildValue(detailed, "KPointsAndWeights", kPointsWeights)

    if (tGroundState) then

      call getChild(detailed, "Occupations", occ)
      do iSpin = 1, nSpin
        call getChild(occ, "spin" // i2c(iSpin), spin)
        do iK = 1, nKPoint
          call getChildValue(spin, "k" // i2c(iK), this%xml%occupations(:, iK, iSpin))
        end do
      end do
      do iK = 1, nKPoint
        this%xml%occupations(:, iK, :) = this%xml%occupations(:, iK, :) *&
            & kPointsWeights(4, iK)
      end do

    else

      call getChild(detailed, "ExcitedOccupations", occ)
      do iSpin = 1, nSpin
        call getChild(occ, "spin" // i2c(iSpin), spin)
        do iK = 1, nKPoint
          call getChildValue(spin, "k" // i2c(iK), this%xml%occupations(:, iK, iSpin))
        end do
      end do
      do iK = 1, nKPoint
        this%xml%occupations(:, iK, :) = this%xml%occupations(:, iK, :) *&
            & kPointsWeights(4, iK)
      end do

    end if

  end subroutine readDetailed


  !> Read in the geometry stored as xml in internal or gen format.
  subroutine readGeometry(this, geonode)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> Node containing the geometry
    type(fnode), pointer :: geonode

    type(fnode), pointer :: child
    type(string) :: buffer

    call getChildValue(geonode, "", child)
    call getNodeName(child, buffer)

    select case (char(buffer))

    case ("genformat")
      call readTGeometryGen(child, this%xml%geo)
      call removeChildNodes(geonode)
      call writeTGeometryHSD(geonode, this%xml%geo)

    case ("xyzformat")
      call readTGeometryXyz(child, this%xml%geo)
      call removeChildNodes(geonode)
      call writeTGeometryHSD(geonode, this%xml%geo)

    case ("vaspformat")
      call readTGeometryVasp(child, this%xml%geo)
      call removeChildNodes(geonode)
      call writeTGeometryHSD(geonode, this%xml%geo)

    case default
      call readTGeometryHSD(geonode, this%xml%geo)

    end select

  end subroutine readGeometry


  !> Interpret the options.
  subroutine readOptions(this, node, nLevel, nKPoint, nSpin, tShiftGrid)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> Node containig the information
    type(fnode), pointer :: node

    !> Nr. of states in the dftb+ calc.
    integer, intent(in) :: nLevel

    !> Nr. of K-points
    integer, intent(in) :: nKPoint

    !> Nr. of spins
    integer, intent(in) :: nSpin

    !> If grid should shifted by a half cell
    logical, intent(out) :: tShiftGrid

    !> List of levels to plot, whereby insignificant occupations were filtered out
    integer, allocatable :: levelIndex(:,:)

    type(fnode), pointer :: subnode, field, value
    type(string) :: buffer, modifier
    type(TListIntR1) :: indexBuffer
    integer :: curId, curVec(3)
    integer :: ind, ii, iLevel, iKPoint, iSpin, iAtom, iSpecies
    logical :: tFound
    real(dp) :: tmpvec(3), minvals(3), maxvals(3)
    real(dp), allocatable :: mcutoffs(:)
    real(dp) :: minEdge

    character(len=63) :: warnId(3) = (/ &
        & "The external files you are providing differ from those provided",&
        & "when this input file was generated. The results you obtain with",&
        & "the current files could be, therefore, different.              "&
        & /)

    !! Warning, if processed input is read in, but eigenvectors are different
    call getChildValue(node, "Identity", curId, this%xml%identity)
    if (curId /= this%xml%identity) then
      call warning(warnId)
    end if

    call getChildValue(node, "TotalChargeDensity", this%option%tPlotTotChrg, .false.)

    if (nSpin == 2) then
      call getChildValue(node, "TotalSpinPolarisation", this%option%tPlotTotSpin, .false.)
    else
      this%option%tPlotTotSpin = .false.
    end if

    call getChildValue(node, "TotalChargeDifference", this%option%tPlotTotDiff, .false.,&
        &child=field)
    call getChildValue(node, "TotalAtomicDensity", this%option%tPlotAtomDens, .false.)
    call getChildValue(node, "ChargeDensity", this%option%tPlotChrg, .false.)
    call getChildValue(node, "ChargeDifference", this%option%tPlotChrgDiff, .false.)

    this%option%tCalcTotChrg = this%option%tPlotTotChrg .or. this%option%tPlotTotSpin .or.&
        & this%option%tPlotTotDiff .or. this%option%tPlotChrgDiff
    this%option%tCalcAtomDens = this%option%tPlotTotDiff .or. this%option%tPlotChrgDiff .or.&
        & this%option%tPlotAtomDens

    call getChildValue(node, "RealComponent", this%option%tPlotReal, .false.)
    call getChildValue(node, "ImagComponent", this%option%tPlotImag, .false., child=field)

    if (this%option%tPlotImag .and. this%xml%tRealHam) then
      call detailedWarning(field, "Wave functions are real, no imaginary part will be plotted")
      this%option%tPlotImag = .false.
    end if

    call getChildValue(node, "PlottedLevels", buffer, child=field, multiple=.true.)
    call getSelectedIndices(node, char(buffer), [1, nLevel], this%option%plottedLevels)

    if (this%xml%geo%tPeriodic) then
      call getChildValue(node, "PlottedKPoints", buffer, child=field, multiple=.true.)
      call getSelectedIndices(node, char(buffer), [1, nKPoint], this%option%plottedKPoints)
    else
      allocate(this%option%plottedKPoints(1))
      this%option%plottedKPoints(1) = 1
    end if

    call getChildValue(node, "PlottedSpins", buffer, child=field, multiple=.true.)
    call getSelectedIndices(node, char(buffer), [1, nSpin], this%option%plottedSpins)

    !! Create the list of the levels, which must be calculated explicitely
    call init(indexBuffer)

    do iSpin = 1, nSpin
      do iKPoint = 1, nKPoint
        do iLevel = 1, nLevel
          tFound = any(this%option%plottedLevels == iLevel) &
              &.and. any(this%option%plottedKPoints == iKPoint) &
              &.and. any(this%option%plottedSpins == iSpin)
          if ((.not. tFound) .and. this%option%tCalcTotChrg) then
            tFound = this%xml%occupations(iLevel, iKPoint, iSpin) > 1e-08_dp
          end if
          if (tFound) then
            call append(indexBuffer, [iLevel, iKPoint, iSpin])
          end if
        end do
      end do
    end do

    if (len(indexBuffer) == 0) then
      call error("No levels specified for plotting")
    end if

    allocate(levelIndex(3, len(indexBuffer)))
    call asArray(indexBuffer, levelIndex)
    call destruct(indexBuffer)

    allocate(this%option%levelIndex(3, size(levelIndex, dim=2)))

    ! Make sure, level entries are correctly sorted in the list
    ind = 1
    do iSpin = 1, nSpin
      do iKPoint = 1, nKPoint
        do iLevel = 1, nLevel
          curVec = [iLevel, iKPoint, iSpin]
          tFound = .false.
          lpLevelIndex: do ii = 1, size(levelIndex, dim=2)
            tFound = all(levelIndex(:, ii) == curVec)
            if (tFound) then
              exit lpLevelIndex
            end if
          end do lpLevelIndex
          if (tFound) then
            this%option%levelIndex(:, ind) = curVec(:)
            ind = ind + 1
          end if
        end do
      end do
    end do

    call getChildValue(node, "SpGridPoints", this%option%nSpPoints, child=field)
    if (any(this%option%nSpPoints <= 0)) then
      call detailedError(field, "Specified numbers must be greater than zero")
    end if

    !! Plotted region: if last (and hopefully only) childnode is not an allowed method -> assume
    !! explicit setting, parse the node "PlottedRegion" for the appropriate children.
    call getChildValue(node, "PlottedRegion", value, child=subnode)
    call getNodeName(value, buffer)

    allocate(this%internal%speciesGridsOrigs(3, this%xml%geo%nSpecies))

    select case (char(buffer))

    case ("unitcell")
      !! Unit cell for the periodic case, smallest possible cuboid for cluster
      if (this%xml%geo%tPeriodic) then
        this%option%origin(:) = [0.0_dp, 0.0_dp, 0.0_dp]
        this%option%boxVecs(:,:) = this%xml%geo%latVecs
      else
        call getChildValue(value, "MinEdgeLength", minEdge, child=field, default=1.0_dp)
        if (minEdge < 0.0_dp) then
          call detailedError(field, "Minimal edge length must be positive")
        end if
        this%option%origin = minval(this%xml%geo%coords, dim=2)
        tmpvec = maxval(this%xml%geo%coords, dim=2) - this%option%origin
        do ii = 1, 3
          if (tmpvec(ii) < minEdge) then
            this%option%origin(ii) = this%option%origin(ii) - 0.5_dp * (minEdge - tmpvec(ii))
            tmpvec(ii) = minEdge
          end if
        end do
        this%option%boxVecs(:,:) = 0.0_dp
        do ii = 1, 3
          this%option%boxVecs(ii, ii) = tmpvec(ii)
        end do
      end if

    case ("optimalcuboid")
      !! Determine optimal cuboid, so that no basis function leaks out
      call getChildValue(value, "MinEdgeLength", minEdge, child=field, default=1.0_dp)
      if (minEdge < 0.0_dp) then
        call detailedError(field, "Minimal edge length must be positive")
      end if
      allocate(mcutoffs(this%xml%geo%nSpecies))
      do iSpecies = 1, this%xml%geo%nSpecies
        mcutoffs(iSpecies) = maxval(this%basis%basis(iSpecies)%cutoffs)
        this%internal%speciesGridsOrigs(:, iSpecies) = - mcutoffs(iSpecies)
      end do
      minvals = this%xml%geo%coords(:, 1)
      maxvals = this%xml%geo%coords(:, 1)
      do iAtom = 1, this%xml%geo%nAtom
        iSpecies = this%xml%geo%species(iAtom)
        maxvals(:) = max(maxvals, this%xml%geo%coords(:, iAtom) + mcutoffs(iSpecies))
        minvals(:) = min(minvals, this%xml%geo%coords(:, iAtom) - mcutoffs(iSpecies))
      end do
      this%option%origin(:) = minvals(:)
      tmpvec(:) = maxvals(:) - minvals(:)
      do ii = 1, 3
        if (tmpvec(ii) < minEdge) then
          this%option%origin(ii) = this%option%origin(ii) - 0.5_dp * (minEdge - tmpvec(ii))
          tmpvec(ii) = minEdge
        end if
      end do
      this%option%boxVecs(:,:) = 0.0_dp
      allocate(this%internal%speciesGridsVecs(3, 3, this%xml%geo%nSpecies))
      this%internal%speciesGridsVecs(:,:,:) = 0.0_dp
      do iSpecies = 1, this%xml%geo%nSpecies
        do ii = 1, 3
          this%internal%speciesGridsVecs(ii, ii, iSpecies) = 2.0_dp * mcutoffs(iSpecies)&
              & / (real(this%option%nSpPoints(ii), dp) - 1.0_dp)
        end do
      end do
      do ii = 1, 3
        this%option%boxVecs(ii, ii) = tmpvec(ii)
      end do

    case ("origin", "box")
      !! Those nodes are part of an explicit specification
      call getChildValue(subnode, "Box", this%option%boxVecs, modifier=modifier,&
          &child=field)
      if (abs(determinant(this%option%boxVecs)) < 1e-08_dp) then
        call detailedError(field, "Vectors are linearly dependent")
      end if
      if (len(modifier) > 0) then
        ind = getModifierIndex(char(modifier), lengthUnits, field)
        this%option%boxVecs(:,:) = this%option%boxVecs(:,:) * lengthUnits(ind)%convertValue
      end if
      call getChildValue(subnode, "Origin", this%option%origin, modifier=modifier, child=field)
      if (len(modifier) > 0) then
        ind = getModifierIndex(char(modifier), lengthUnits, field)
        this%option%origin(:) = this%option%origin * lengthUnits(ind)%convertValue
      end if

    case default
      !! Object with unknown name passed
      call detailedError(value, "Invalid element name")
    end select

    !! Replace existing PlottedRegion definition
    call setChild(node, "PlottedRegion", field, replace=.true.)
    call setChildValue(field, "Origin", this%option%origin, .true.)
    call setChildValue(field, "Box", this%option%boxVecs, .true.)

    call getChildValue(node, "TotGridPoints", this%option%nTotPoints, child=field)

    if (any(this%option%nTotPoints <= 0)) then
      call detailedError(field, "Specified numbers must be greater than zero")
    end if

    call getChildValue(node, "ShiftGrid", tShiftGrid, default=.true.)

    if (this%xml%geo%tPeriodic) then
      call getChildValue(node, "FoldAtomsToUnitCell", this%option%tFoldCoords, default=.false.)
      call getChildValue(node, "FillBoxWithAtoms", this%option%tFillBox, default=.false.)
      this%option%tFoldCoords = this%option%tFoldCoords .or. this%option%tFillBox
    else
      this%option%tFillBox = .false.
      this%option%tFoldCoords = .false.
    end if

    call getChildValue(node, "RepeatBox", this%option%repeatBox, default=[1, 1, 1], child=field)

    if (.not. all(this%option%repeatBox > 0)) then
      call detailedError(field, "Indexes must be greater than zero")
    end if

    call getChildValue(node, "Verbose", this%option%tVerbose, .false.)

    !! Create grid vectors, shift them if necessary
    do ii = 1, 3
      this%internal%totGridVec(:, ii) = this%option%boxVecs(:, ii)&
          & / real(this%option%nTotPoints(ii), dp)
    end do

    if (char(buffer) == 'origin' .or. char(buffer) == 'box' .or. char(buffer) == 'unitcell') then
      call calculateSpeciesGridsVecs(this)
    end if

    if (tShiftGrid) then
      this%option%totGridOrig(:) = this%option%origin(:) +&
          & 0.5_dp * sum(this%internal%totGridVec, dim=2)
      ! do iSpecies = 1, this%xml%geo%nSpecies
      !   this%internal%speciesGridsOrigs(:, iSpecies) =&
      !       & this%internal%speciesGridsOrigs(:, iSpecies) + 0.5_dp *&
      !       & sum(this%internal%totGridVec, dim=2)
      ! end do
    else
      this%option%totGridOrig(:) = this%option%origin(:)
    end if
    this%internal%gridVol = determinant(this%internal%totGridVec)

    !! Calculate the center coordinate of the total grid
    this%internal%totGridCenter(:) = 0.0_dp
    do ii = 1, 3
      this%internal%totGridCenter(:) = this%internal%totGridCenter + 0.5_dp *&
          & real(this%option%nTotPoints(ii), dp) * this%internal%totGridVec(:, ii)
    end do

  end subroutine readOptions


  !> Calculates box vectors of species grids for given total grid vectors
  subroutine calculateSpeciesGridsVecs(this)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> cartesian total grid coordinates
    real(dp), allocatable :: subcubeCartCoords(:,:)

    !> integer total grid coordinates
    integer :: gridcoords(3,8)

    !> Inverse of total grid basis
    real(dp) :: totGridInvBasis(3,3)

    !> Box vectors for current species
    real(dp) :: speciesBoxVecs(3,3)

    !> Maximum cutoff radius for a species
    real(dp) :: maxSpeciesCutoff

    !> Integer total grid ranges that contain the cutoff subcube
    integer :: maxRanges(3)

    !> Auxiliary variables
    integer :: ii, iSpecies

    allocate(this%internal%speciesGridsVecs(3, 3, this%xml%geo%nSpecies))

    call invert33(totGridInvBasis, this%internal%totGridVec)

    do iSpecies = 1, this%xml%geo%nSpecies

      maxSpeciesCutoff = maxval(this%basis%basis(iSpecies)%cutoffs)
      subcubeCartCoords = reshape([&
          0.0_dp, 0.0_dp, 0.0_dp,&
          0.0_dp, 0.0_dp, 1.0_dp,&
          0.0_dp, 1.0_dp, 0.0_dp,&
          0.0_dp, 1.0_dp, 1.0_dp,&
          1.0_dp, 0.0_dp, 0.0_dp,&
          1.0_dp, 0.0_dp, 1.0_dp,&
          1.0_dp, 1.0_dp, 0.0_dp,&
          1.0_dp, 1.0_dp, 1.0_dp],&
          & [3, 8]) * 2.0_dp * maxSpeciesCutoff

      do ii = 1, size(subcubeCartCoords, dim=2)
        gridcoords(:, ii) = floor(matmul(subcubeCartCoords(:, ii), totGridInvBasis))
      end do

      maxRanges = maxval(gridcoords, dim=2) - minval(gridcoords, dim=2)

      do ii = 1, 3
        speciesBoxVecs(:, ii) = this%internal%totGridVec(:, ii) * maxRanges(ii)
      end do

      this%internal%speciesGridsOrigs(:, iSpecies) = [0.0_dp, 0.0_dp, 0.0_dp]
      do ii = 1, 3
        this%internal%speciesGridsVecs(:, ii, iSpecies) = speciesBoxVecs(:, ii)&
            & / real(this%option%nSpPoints(ii), dp)
        this%internal%speciesGridsOrigs(:, iSpecies) = this%internal%speciesGridsOrigs(:, iSpecies)&
            & - 0.5_dp * real(this%option%nSpPoints(ii), dp) *&
            & this%internal%speciesGridsVecs(:, ii, iSpecies)
      end do

    end do

  end subroutine calculateSpeciesGridsVecs


#:for DTYPE, NAME in [('complex', 'Cplx'), ('real', 'Real')]

  !> Read external eigenvector file (eigenvec.bin)
  subroutine read${NAME}$Eigenvecs(filename, eigenvecs, jobId)

    !> File to read eigenvectors from
    character(len=*), intent(in) :: filename

    !> Resulting eigenvectors read from file
    ${DTYPE}$(dp), intent(out) :: eigenvecs(:,:)

    !> ID of the calculation which produced the file
    integer, intent(out), optional :: jobId

    !> Number of orbitals per state
    integer :: nOrb

    integer :: fd, iOrb, dummy
    logical :: exst

    nOrb = size(eigenvecs, dim=1)

    inquire(file=filename, exist=exst)

    if (exst) then

      open(newunit=fd, file=filename, action="read", form="unformatted", position='rewind')
      read(fd) dummy

      if (present(jobId)) then
        jobId = dummy
      end if

      do iOrb = 1, nOrb
        read(fd) eigenvecs(:, iOrb)
      end do

      close(fd)

    else
      call error('no ' // filename // ' file!')
    end if

  end subroutine read${NAME}$Eigenvecs

#:endfor


  !> Read in the basis related informations
  subroutine readBasis(this, node, speciesNames)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> Node containing the basis definition
    type(fnode), pointer :: node

    !> Names of the species for which basis should be read in
    character(len=*), intent(in) :: speciesNames(:)

    character(len=len(speciesNames)) :: speciesName
    type(fnode), pointer :: speciesNode
    integer :: nSpecies
    integer :: ii

    nSpecies = size(speciesNames)

    @:ASSERT(nSpecies > 0)

    allocate(this%basis%basis(nSpecies))
    allocate(this%atomicNumber%atomicNumbers(nSpecies))

    do ii = 1, nSpecies
      speciesName = speciesNames(ii)
      call getChild(node, speciesName, speciesNode)
      call readSpeciesBasis(speciesNode, this%basis%basis(ii))
      this%atomicNumber%atomicNumbers(ii) = this%basis%basis(ii)%atomicNumber
    end do

  end subroutine readBasis


  !> Read in basis function for a species.
  subroutine readSpeciesBasis(node, spBasis)

    !> Node containing the basis definition for a species
    type(fnode), pointer :: node

    !> Contains the basis on return
    type(TSpeciesBasis), intent(out) :: spBasis

    type(fnode), pointer :: tmpNode, child
    type(fnodeList), pointer :: children
    type(TListReal) :: bufferRwf
    real(dp), allocatable :: rwf(:)
    integer :: ii

    call getChildValue(node, "AtomicNumber", spBasis%atomicNumber)
    call getChildren(node, "Orbital", children)
    spBasis%nOrb = getLength(children)

    if (spBasis%nOrb < 1) then
      call detailedError(node, "Missing orbital definitions")
    end if

    allocate(spBasis%angMoms(spBasis%nOrb))
    allocate(spBasis%occupations(spBasis%nOrb))
    allocate(spBasis%rwfs(spBasis%nOrb))
    allocate(spBasis%cutoffs(spBasis%nOrb))

    do ii = 1, spBasis%nOrb

      call getItem1(children, ii, tmpNode)
      call getChildValue(tmpNode, "AngularMomentum", spBasis%angMoms(ii))
      call getChildValue(tmpNode, "Occupation", spBasis%occupations(ii))
      call getChildValue(tmpNode, "Cutoff", spBasis%cutoffs(ii))

      call init(bufferRwf)
      call getChildValue(tmpNode, "RadialWaveFunc", bufferRwf, child=child)

      if (len(bufferRwf) == 0) then
        call detailedError(child, "Missing radial wavefunction")
      end if

      allocate(rwf(len(bufferRwf)))

      call asArray(bufferRwf, rwf)
      call destruct(bufferRwf)

      call init(spBasis%rwfs(ii), transpose(reshape(rwf, [3, size(rwf) / 3])))

      deallocate(rwf)

    end do

  end subroutine readSpeciesBasis


  !> Determines convenient index mappings
  subroutine getIndexMappings(this)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> Auxiliary arrays
    integer, allocatable :: tmp1(:), tmp2(:), tmp3(:,:)

    !> Auxiliary variables
    integer :: iAtom, iOrb, iSpecies, iAng, iL, iM, iOrbPerSpecies, mAng, ind

    allocate(this%internal%orbitalOcc(this%xml%nOrb, 1))
    allocate(this%internal%orbitalToAtom(this%xml%nOrb))
    allocate(this%internal%orbitalToSpecies(this%xml%nOrb))
    allocate(tmp1(this%xml%geo%nSpecies + 1))
    allocate(tmp2(this%xml%nOrb))
    allocate(tmp3(2, this%xml%nOrb))

    iOrb = 1

    do iAtom = 1, this%xml%geo%nAtom
      ind = 1
      iSpecies = this%xml%geo%species(iAtom)
      do iAng = 1, size(this%basis%basis(iSpecies)%angMoms)
        mAng = 2 * this%basis%basis(iSpecies)%angMoms(iAng) + 1
        this%internal%orbitalOcc(iOrb:iOrb + mAng - 1, 1) =&
            & this%basis%basis(iSpecies)%occupations(iAng) / real(mAng, dp)
        do iM = 1, mAng
          tmp3(1, iOrb) = iSpecies
          tmp3(2, iOrb) = ind
          this%internal%orbitalToAtom(iOrb) = iAtom
          iOrb = iOrb + 1
          ind = ind + 1
        end do
      end do
    end do

    ind = 1

    do iAtom = 1, this%xml%geo%nAtom
      iSpecies = this%xml%geo%species(iAtom)
      do iAng = 1, size(this%basis%basis(iSpecies)%angMoms)
        mAng = 2 * this%basis%basis(iSpecies)%angMoms(iAng) + 1
        this%internal%orbitalOcc(ind:ind + mAng - 1, 1) =&
            & this%basis%basis(iSpecies)%occupations(iAng) / real(mAng, dp)
        ind = ind + mAng
      end do
    end do

    ind = 1

    do iSpecies = 1, this%xml%geo%nSpecies
      tmp1(iSpecies) = ind
      mAng = maxval(this%basis%basis(iSpecies)%angMoms)
      do iAng = this%internal%molorb%iStos(iSpecies), this%internal%molorb%iStos(iSpecies + 1) - 1
        iL = this%internal%molorb%angMoms(iAng)
        ind = ind + 2 * iL + 1
      end do
    end do

    tmp1(this%xml%geo%nSpecies + 1) = ind


    do iOrb = 1, this%xml%nOrb
      iSpecies = tmp3(1, iOrb)
      iOrbPerSpecies = tmp3(2, iOrb)
      this%internal%orbitalToSpecies(iOrb) = tmp1(iSpecies) + iOrbPerSpecies - 1
    end do

    deallocate(tmp1)
    deallocate(tmp2)
    deallocate(tmp3)

  end subroutine getIndexMappings


  !> Repeat box if necessary
  subroutine getRepeatedBox(this)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> Temporary storage
    real(dp), allocatable :: tmpCoords(:,:)
    integer, allocatable :: tmpSpecies(:)

    !> Amount to shift
    real(dp) :: shift(3)

    !> Total number of boxes
    integer :: nBox

    !> Auxiliary variables
    integer :: i1, i2, i3, iAtom, ind

    nBox = product(this%option%repeatBox)

    if (nBox > 1) then

      ! If tFillBox is off, coordinates must be repeated here.
      ! Otherwise the part for filling with atoms will do that.
      if (.not. this%option%tFillBox) then

        allocate(tmpCoords(3, size(this%xml%geo%coords)))
        allocate(tmpSpecies(size(this%xml%geo%species)))

        ! temporarily save unrepeated coordinates and species
        tmpCoords(:,:) = this%xml%geo%coords(:,:)
        tmpSpecies(:) = this%xml%geo%species(:)

        deallocate(this%xml%geo%coords)
        deallocate(this%xml%geo%species)
        allocate(this%xml%geo%coords(3, nBox * this%xml%geo%nAtom))
        allocate(this%xml%geo%species(nBox * this%xml%geo%nAtom))

        ind = 0

        do i1 = 0, this%option%repeatBox(1) - 1
          do i2 = 0, this%option%repeatBox(2) - 1
            do i3 = 0, this%option%repeatBox(3) - 1
              shift(:) = matmul(this%option%boxVecs, real([i1, i2, i3], dp))
              do iAtom = 1, this%xml%geo%nAtom
                this%xml%geo%coords(:, ind + iAtom) = tmpCoords(:, iAtom) + shift(:)
              end do
              this%xml%geo%species(ind + 1:ind + this%xml%geo%nAtom) = tmpSpecies(:)
              ind = ind + this%xml%geo%nAtom
            end do
          end do
        end do
        this%xml%geo%nAtom = nBox * this%xml%geo%nAtom
      end if
      do i1 = 1, 3
        this%option%boxVecs(:, i1) = this%option%boxVecs(:, i1) *&
            & real(this%option%repeatBox(i1), dp)
      end do

    end if

  end subroutine getRepeatedBox


  !> Checks, if the eigenvector file has the right identity number.
  subroutine checkEigenvecs(fileName, identity)

    !> File to check
    character(len=*), intent(in) :: fileName

    !> Identity number.
    integer, intent(in) :: identity

    integer :: fd, id, iostat

    fd = getFileId()

    open(fd, file=fileName, action="read", position="rewind", form="unformatted", iostat=iostat)

    if (iostat /= 0) then
      call error("Can't open file '" // trim(fileName) // "'.")
    end if

    read (fd) id

    if (id /= identity) then
      call error("Ids for eigenvectors ("// i2c(id) //") and xml-input ("// i2c(identity) // &
          & ") don't match.")
    end if

    close(fd)

  end subroutine checkEigenvecs


  !> Determinant of a 3x3 matrix (Only temporary!)
  function  determinant(matrix)

    !> The matrix to calculate the determinant from.
    real(dp), intent(in) :: matrix(:, :)

    !> Determinant of the matrix.
    real(dp) :: determinant

    real(dp) :: tmp

    @:ASSERT(all(shape(matrix) == (/3, 3/)))

    tmp = matrix(1, 1) &
        &* (matrix(2, 2) * matrix(3, 3) - matrix(3, 2) * matrix(2, 3))
    tmp = tmp - matrix(1, 2) &
        &* (matrix(2, 1) * matrix(3, 3) - matrix(3, 1) * matrix(2, 3))
    tmp = tmp + matrix(1, 3) &
        &* (matrix(2, 1) * matrix(3, 2) - matrix(3, 1) * matrix(2, 2))

    determinant = abs(tmp)

  end function determinant

end module waveplot_initwaveplot
