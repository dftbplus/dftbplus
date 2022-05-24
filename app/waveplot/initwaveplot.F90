!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the routines for initialising waveplot.
module waveplot_initwaveplot
  use dftbp_common_accuracy, only : dp
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_status, only : TStatus
  use dftbp_common_unitconversion, only : lengthUnits
  use dftbp_dftb_boundarycond, only : boundaryConditions, TBoundaryConditions,&
      & TBoundaryConditions_init
  use dftbp_extlibs_xmlf90, only : fnode, fNodeList, string, char, getLength, getItem1,&
      & getNodeName,destroyNode
  use dftbp_io_charmanip, only : i2c, unquote
  use dftbp_io_hsdparser, only : parseHSD, dumpHSD
  use dftbp_io_hsdutils, only : getChildValue, setChildValue, getChild, setChild, getChildren,&
      & getSelectedIndices, detailedError, detailedWarning
  use dftbp_io_hsdutils2, only : convertUnitHsd, readHSDAsXML, warnUnprocessedNodes
  use dftbp_io_message, only : warning, error
  use dftbp_io_xmlutils, only : removeChildNodes
  use dftbp_type_linkedlist, only : TListIntR1, TListReal, init, destruct, len, append, asArray
  use dftbp_type_typegeometryhsd, only : TGeometry, readTGeometryGen, readTGeometryHSD,&
      & readTGeometryVasp, readTGeometryXyz, writeTGeometryHSD
  use waveplot_gridcache, only : TGridCache, init
  use waveplot_molorb, only : TMolecularOrbital, TSpeciesBasis
  implicit none

  private
  public :: TProgramVariables, TProgramVariables_init


  !> Data type containing variables from detailed.xml
  type TInput

    !> Geometry instance
    type(TGeometry) :: geo

    !> Identity of the run
    integer :: identity

    !> Nr. of orbitals per state
    integer :: nOrb

    !> True, if eigenvectors/hamiltonian is real-valued
    logical :: tRealHam

    !> Occupations
    real(dp), allocatable :: occupations(:,:,:)

  end type TInput


  !> Data type containing variables from the Option block
  type TOption

    !> Nr of grid points along 3 directions
    integer :: nPoints(3)

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

    !> Origin of the grid in the box
    real(dp) :: gridOrigin(3)

    !> List of levels to plot, whereby insignificant occupations were filtered out
    integer, allocatable :: levelIndex(:,:)

  end type TOption


  !> Data type containing variables from the Basis block
  type TBasis

    !> Definition of the wfcs
    type(TSpeciesBasis), allocatable :: basis(:)

    !> Resolution of the radial wfcs
    real(dp) :: basisResolution

  end type TBasis


  !> Data type containing variables from eigenvec.bin
  type TEig

    !> Nr. of states
    integer :: nState

    !> Real eigenvectors
    real(dp), allocatable :: eigvecsReal(:,:)

    !> Complex eigenvectors
    complex(dp), allocatable :: eigvecsCplx(:,:)

  end type TEig


  !> Data type containing variables from the AtomicNumbers block
  type TAtomicNumber

    !> Species-atomic nr. corresp.
    integer, allocatable :: atomicNumbers(:)

  end type TAtomicNumber


  !> Data type containing locally created variables
  type TInternal

    !> Molecular orbital
    type(TMolecularOrbital), allocatable :: molOrb

    !> pointer to the orbital
    type(TMolecularOrbital), pointer :: pMolOrb

    !> Grid cache
    type(TGridCache) :: grid

    !> Grid vectors
    real(dp) :: gridVec(3,3)

    !> Volume of the grid
    real(dp) :: gridVol

    !> List of levels to plot
    integer, allocatable :: levelIndex(:,:)

  end type TInternal


  !> Data type containing program variables
  type TProgramVariables

    !> Data of detailed.xml
    type(TInput) :: input

    !> Data of eigenvec.bin
    type(TEig) :: eig

    !> Data of Option block
    type(TOption) :: opt

    !> Boundary condition
    type(TBoundaryConditions) :: boundaryCond

    !> Data of Basis block
    type(TBasis) :: basis

    !> Data of AtomicNumber block
    type(TAtomicNumber) :: aNr

    !> Locally created data
    type(TInternal) :: loc

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

contains


  !> Initialise program variables
  subroutine TProgramVariables_init(this)

    !> Container of program variables
    type(TProgramVariables), intent(out), target :: this

    !> Pointers to input nodes
    type(fnode), pointer :: root, tmp, detailed, hsdTree

    !> String buffer instance
    type(string) :: strBuffer

    !> Id of input/parser version
    integer :: inputVersion

    !> Nr. of cached grids
    integer :: nCached

    !> Nr. of K-points
    integer :: nKPoint

    !> Nr. of spins
    integer :: nSpin

    !> Wether to look for ground state occupations (True) or excited (False)
    logical :: tGroundState

    !> If grid should shifted by a half cell
    logical :: tShiftGrid

    !> K-points and weights
    real(dp), allocatable :: kPointsWeights(:,:)

    !> File with binary eigenvectors
    character(len=1024) :: eigVecBin

    !> Auxiliary variable
    integer :: ii

    !> Operation status, if an error needs to be returned
    type(TStatus) :: errStatus

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
    call readDetailed(this, detailed, tGroundState, kPointsWeights)
    call destroyNode(tmp)

    nKPoint = size(kPointsWeights, dim=2)
    nSpin = size(this%input%occupations, dim=3)
    this%eig%nState = size(this%input%occupations, dim=1)

    !! Read basis
    call getChild(root, "Basis", tmp)
    call readBasis(this, tmp, this%input%geo%speciesNames)
    call getChildValue(root, "EigenvecBin", strBuffer)
    eigVecBin = unquote(char(strBuffer))
    call checkEigenvecs(eigVecBin, this%input%identity)

    !! Read options
    call getChild(root, "Options", tmp)
    call readOptions(this, tmp, this%eig%nState, nKPoint, nSpin, nCached, tShiftGrid)

    !! Issue warning about unprocessed nodes
    call warnUnprocessedNodes(root, .true.)

    !! Finish parsing, dump parsed and processed input
    call dumpHSD(hsdTree, hsdParsedInput)
    write(stdout, "(A)") "Processed input written as HSD to '" // hsdParsedInput &
        &//"'"
    write(stdout, "(A,/)") repeat("-", 80)
    call destroyNode(hsdTree)

    if (this%input%geo%tPeriodic) then
      call TBoundaryConditions_init(this%boundaryCond, boundaryConditions%pbc3d, errStatus)
    else if (this%input%geo%tHelical) then
      call TBoundaryConditions_init(this%boundaryCond, boundaryConditions%helical, errStatus)
    else
      call TBoundaryConditions_init(this%boundaryCond, boundaryConditions%cluster, errStatus)
    end if
    if (errStatus%hasError()) then
      call error(errStatus%message)
    end if

    !! Create grid vectors, shift them if necessary
    do ii = 1, 3
      this%loc%gridVec(:, ii) = this%opt%boxVecs(:, ii) / real(this%opt%nPoints(ii), dp)
    end do
    if (tShiftGrid) then
      this%opt%gridOrigin(:) = this%opt%origin(:) + 0.5_dp * sum(this%loc%gridVec, dim=2)
    else
      this%opt%gridOrigin(:) = this%opt%origin(:)
    end if
    this%loc%gridVol = determinant(this%loc%gridVec)

    write(stdout, "(A)") "Doing initialisation"

    !! Initialize necessary (molecular orbital, grid) objects
    allocate(this%loc%molOrb)
    this%loc%pMolOrb => this%loc%molOrb
    call init(this%loc%molOrb, this%input%geo, this%boundaryCond, this%basis%basis)

    call init(this%loc%grid, levelIndex=this%loc%levelIndex, nOrb=this%input%nOrb,&
        & nAllLevel=this%eig%nState, nAllKPoint=nKPoint, nAllSpin=nSpin, nCached=nCached,&
        & nPoints=this%opt%nPoints, tVerbose=this%opt%tVerbose, eigVecBin=eigVecBin,&
        & gridVec=this%loc%gridVec, origin=this%opt%gridOrigin,&
        & kPointCoords=kPointsWeights(1:3, :), tReal=this%input%tRealHam, molorb=this%loc%pMolOrb)

  end subroutine TProgramVariables_init


  !> Interpret the information stored in detailed.xml
  subroutine readDetailed(this, detailed, tGroundState, kPointsWeights)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> Pointer to the node, containing the information
    type(fnode), pointer :: detailed

    !> Wether to look for ground state occupations (True) or excited (False)
    logical, intent(in) :: tGroundState

    !> K-points and weights
    real(dp), intent(out), allocatable :: kPointsWeights(:,:)

    !> Pointers to input nodes
    type(fnode), pointer :: tmp, occ, spin

    !> Nr. of K-points
    integer :: nKPoint

    !> Nr. of spins
    integer :: nSpin

    !> Nr. of states
    integer :: nState

    !> Auxiliary variables
    integer :: iSpin, iKpoint

    call getChildValue(detailed, "Identity", this%input%identity)
    call getChild(detailed, "Geometry", tmp)
    call readGeometry(this%input%geo, tmp)

    call getChildValue(detailed, "Real", this%input%tRealHam)
    call getChildValue(detailed, "NrOfKPoints", nKPoint)
    call getChildValue(detailed, "NrOfSpins", nSpin)
    call getChildValue(detailed, "NrOfStates", nState)
    call getChildValue(detailed, "NrOfOrbitals", this%input%nOrb)

    allocate(kPointsWeights(4, nKPoint))
    allocate(this%input%occupations(nState, nKPoint, nSpin))

    call getChildValue(detailed, "KPointsAndWeights", kPointsWeights)

    if (tGroundState) then
      call getChild(detailed, "Occupations", occ)
      do iSpin = 1, nSpin
        call getChild(occ, "spin" // i2c(iSpin), spin)
        do iKpoint = 1, nKPoint
          call getChildValue(spin, "k" // i2c(iKpoint), this%input%occupations(:, iKpoint, iSpin))
        end do
      end do
      do iKpoint = 1, nKPoint
        this%input%occupations(:, iKpoint, :) = this%input%occupations(:, iKpoint, :)&
            & * kPointsWeights(4, iKpoint)
      end do
    else
      call getChild(detailed, "ExcitedOccupations", occ)
      do iSpin = 1, nSpin
        call getChild(occ, "spin" // i2c(iSpin), spin)
        do iKpoint = 1, nKPoint
          call getChildValue(spin, "k" // i2c(iKpoint), this%input%occupations(:, iKpoint, iSpin))
        end do
      end do
      do iKpoint = 1, nKPoint
        this%input%occupations(:, iKpoint, :) = this%input%occupations(:, iKpoint, :)&
            & * kPointsWeights(4, iKpoint)
      end do
    end if

  end subroutine readDetailed


  !> Read in the geometry stored as xml in internal or gen format.
  subroutine readGeometry(geo, geonode)

    !> Geometry instance
    type(TGeometry), intent(out) :: geo

    !> Node containing the geometry
    type(fnode), pointer :: geonode

    !> Pointers to input nodes
    type(fnode), pointer :: child

    !> String buffer instance
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


  !> Interpret the options.
  subroutine readOptions(this, node, nLevel, nKPoint, nSpin, nCached, tShiftGrid)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> Node containing the information
    type(fnode), pointer :: node

    !> Nr. of states in the calculation
    integer, intent(in) :: nLevel

    !> Nr. of K-points
    integer, intent(in) :: nKPoint

    !> Nr. of spins
    integer, intent(in) :: nSpin

    !> Nr. of cached grids
    integer, intent(out) :: nCached

    !> If grid should be shifted by half a cell
    logical, intent(out) :: tShiftGrid

    !> Pointer to the nodes, containing the information
    type(fnode), pointer :: subnode, field, value

    !> String buffer instances
    type(string) :: buffer, modifier

    !> Onedimensional integer-valued index list
    type(TListIntR1) :: indexBuffer

    !> Id of calculation at hand
    integer :: curId

    !> If current level is found be calculated explicitely
    logical :: tFound

    !> Warning issued, if the detailed.xml id does not match the eigenvector id
    character(len=63) :: warnId(3) = [&
        & "The external files you are providing differ from those provided", &
        & "when this input file was generated. The results you obtain with", &
        & "the current files could therefore be different.                "]

    !> Auxiliary variables
    integer :: ii, iLevel, iKPoint, iSpin, iAtom, iSpecies
    real(dp) :: tmpvec(3), minvals(3), maxvals(3)
    real(dp), allocatable :: mcutoffs(:)
    real(dp) :: minEdge

    !! Warning, if processed input is read in, but eigenvectors are different
    call getChildValue(node, "Identity", curId, this%input%identity)
    if (curId /= this%input%identity) then
      call warning(warnId)
    end if

    call getChildValue(node, "TotalChargeDensity", this%opt%tPlotTotChrg, .false.)

    if (nSpin == 2) then
      call getChildValue(node, "TotalSpinPolarisation", this%opt%tPlotTotSpin, .false.)
    else
      this%opt%tPlotTotSpin = .false.
    end if

    call getChildValue(node, "TotalChargeDifference", this%opt%tPlotTotDiff, .false., child=field)
    call getChildValue(node, "TotalAtomicDensity", this%opt%tPlotAtomDens, .false.)
    call getChildValue(node, "ChargeDensity", this%opt%tPlotChrg, .false.)
    call getChildValue(node, "ChargeDifference", this%opt%tPlotChrgDiff, .false.)

    this%opt%tCalcTotChrg = this%opt%tPlotTotChrg .or. this%opt%tPlotTotSpin&
        & .or. this%opt%tPlotTotDiff .or. this%opt%tPlotChrgDiff
    this%opt%tCalcAtomDens = this%opt%tPlotTotDiff .or. this%opt%tPlotChrgDiff&
        & .or. this%opt%tPlotAtomDens

    call getChildValue(node, "RealComponent", this%opt%tPlotReal, .false.)
    call getChildValue(node, "ImagComponent", this%opt%tPlotImag, .false., child=field)

    if (this%opt%tPlotImag .and. this%input%tRealHam) then
      call detailedWarning(field, "Wave functions are real, no imaginary part will be plotted")
      this%opt%tPlotImag = .false.
    end if

    call getChildValue(node, "PlottedLevels", buffer, child=field, multiple=.true.)
    call getSelectedIndices(node, char(buffer), [1, nLevel], this%opt%plottedLevels)

    if (this%input%geo%tPeriodic) then
      call getChildValue(node, "PlottedKPoints", buffer, child=field, multiple=.true.)
      call getSelectedIndices(node, char(buffer), [1, nKPoint], this%opt%plottedKPoints)
    else
      allocate(this%opt%plottedKPoints(1))
      this%opt%plottedKPoints(1) = 1
    end if

    call getChildValue(node, "PlottedSpins", buffer, child=field, multiple=.true.)
    call getSelectedIndices(node, char(buffer), [1, nSpin], this%opt%plottedSpins)

    !! Create the list of the levels, which must be calculated explicitely
    call init(indexBuffer)
    do iSpin = 1, nSpin
      do iKPoint = 1, nKPoint
        do iLevel = 1, nLevel
          tFound = any(this%opt%plottedLevels == iLevel)&
              & .and. any(this%opt%plottedKPoints == iKPoint)&
              & .and. any(this%opt%plottedSpins == iSpin)
          if ((.not. tFound) .and. this%opt%tCalcTotChrg) then
            tFound = this%input%occupations(iLevel, iKPoint, iSpin) > 1e-08_dp
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

    allocate(this%loc%levelIndex(3, len(indexBuffer)))
    call asArray(indexBuffer, this%loc%levelIndex)
    call destruct(indexBuffer)

    call getChildValue(node, "NrOfCachedGrids", nCached, 1, child=field)

    if (nCached < 1 .and. nCached /= -1) then
      call detailedError(field, "Value must be -1 or greater than zero.")
    end if

    if (nCached == -1) then
      nCached = size(this%loc%levelIndex, dim=2)
    end if

    !! Plotted region: if last (and hopefully only) childnode is not an allowed method -> assume
    !! explicit setting, parse the node "PlottedRegion" for the appropriate children.
    call getChildValue(node, "PlottedRegion", value, child=subnode)
    call getNodeName(value, buffer)

    select case (char(buffer))
    case ("unitcell")
      !! Unit cell for the periodic case, smallest possible cuboid for cluster
      if (this%input%geo%tPeriodic) then
        this%opt%origin(:) = [0.0_dp, 0.0_dp, 0.0_dp]
        this%opt%boxVecs(:,:) = this%input%geo%latVecs(:,:)
      else
        call getChildValue(value, "MinEdgeLength", minEdge, child=field, default=1.0_dp)
        if (minEdge < 0.0_dp) then
          call detailedError(field, "Minimal edge length must be positive")
        end if
        this%opt%origin = minval(this%input%geo%coords, dim=2)
        tmpvec = maxval(this%input%geo%coords, dim=2) - this%opt%origin
        do ii = 1, 3
          if (tmpvec(ii) < minEdge) then
            this%opt%origin(ii) = this%opt%origin(ii) - 0.5_dp * (minEdge - tmpvec(ii))
            tmpvec(ii) = minEdge
          end if
        end do
        this%opt%boxVecs(:,:) = 0.0_dp
        do ii = 1, 3
          this%opt%boxVecs(ii, ii) = tmpvec(ii)
        end do
      end if

    case ("optimalcuboid")
      !! Determine optimal cuboid, so that no basis function leaks out
      call getChildValue(value, "MinEdgeLength", minEdge, child=field, default=1.0_dp)
      if (minEdge < 0.0_dp) then
        call detailedError(field, "Minimal edge length must be positive")
      end if
      allocate(mcutoffs(this%input%geo%nSpecies))
      do iSpecies = 1 , this%input%geo%nSpecies
        mcutoffs(iSpecies) = maxval(this%basis%basis(iSpecies)%cutoffs)
      end do
      minvals = this%input%geo%coords(:,1)
      maxvals = this%input%geo%coords(:,1)
      do iAtom = 1, this%input%geo%nAtom
        iSpecies = this%input%geo%species(iAtom)
        maxvals(:) = max(maxvals, this%input%geo%coords(:, iAtom) + mcutoffs(iSpecies))
        minvals(:) = min(minvals, this%input%geo%coords(:, iAtom) - mcutoffs(iSpecies))
      end do
      this%opt%origin(:) = minvals(:)
      tmpvec(:) = maxvals(:) - minvals(:)
      do ii = 1, 3
        if (tmpvec(ii) < minEdge) then
          this%opt%origin(ii) = this%opt%origin(ii) - 0.5_dp * (minEdge - tmpvec(ii))
          tmpvec(ii) = minEdge
        end if
      end do
      this%opt%boxVecs(:,:) = 0.0_dp
      do ii = 1, 3
        this%opt%boxVecs(ii, ii) = tmpvec(ii)
      end do

    case ("origin","box")
      !! Those nodes are part of an explicit specification -> explitic specif
      call getChildValue(subnode, "Box", this%opt%boxVecs, modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), lengthUnits, field, this%opt%boxVecs)
      if (abs(determinant(this%opt%boxVecs)) < 1e-08_dp) then
        call detailedError(field, "Vectors are linearly dependent")
      end if
      call getChildValue(subnode, "Origin", this%opt%origin, modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), lengthUnits, field, this%opt%origin)

    case default
      !! Object with unknown name passed
      call detailedError(value, "Invalid element name")
    end select

    !! Replace existing PlottedRegion definition
    call setChild(node, "PlottedRegion", field, replace=.true.)
    call setChildValue(field, "Origin", this%opt%origin, replace=.true.)
    call setChildValue(field, "Box", this%opt%boxVecs, replace=.true.)

    call getChildValue(node, "NrOfPoints", this%opt%nPoints, child=field)

    if (any(this%opt%nPoints <= 0)) then
      call detailedError(field, "Specified numbers must be greater than zero")
    end if

    call getChildValue(node, "ShiftGrid", tShiftGrid, default=.true.)

    if (this%input%geo%tPeriodic) then
      call getChildValue(node, "FoldAtomsToUnitCell", this%opt%tFoldCoords, default=.false.)
      call getChildValue(node, "FillBoxWithAtoms", this%opt%tFillBox, default=.false.)
      this%opt%tFoldCoords = this%opt%tFoldCoords .or. this%opt%tFillBox
    else
      this%opt%tFillBox = .false.
      this%opt%tFoldCoords = .false.
    end if

    call getChildValue(node, "RepeatBox", this%opt%repeatBox, default=[1, 1, 1], child=field)

    if (.not. all(this%opt%repeatBox > 0)) then
      call detailedError(field, "Indexes must be greater than zero")
    end if

    call getChildValue(node, "Verbose", this%opt%tVerbose, default=.false.)

  end subroutine readOptions


  !> Read in the basis related informations
  subroutine readBasis(this, node, speciesNames)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> Node containing the basis definition
    type(fnode), pointer :: node

    !> Names of the species for which the basis should be read in
    character(len=*), intent(in) :: speciesNames(:)

    !> Name of current species
    character(len=len(speciesNames)) :: speciesName

    !> Input node instance, containing the information
    type(fnode), pointer :: speciesNode

    !> Total number of species in the system
    integer :: nSpecies

    !> Auxiliary variable
    integer :: ii

    nSpecies = size(speciesNames)

    @:ASSERT(nSpecies > 0)

    call getChildValue(node, "Resolution", this%basis%basisResolution)

    allocate(this%basis%basis(nSpecies))
    allocate(this%aNr%atomicNumbers(nSpecies))

    do ii = 1, nSpecies
      speciesName = speciesNames(ii)
      call getChild(node, speciesName, speciesNode)
      call readSpeciesBasis(speciesNode, this%basis%basisResolution, this%basis%basis(ii))
      this%aNr%atomicNumbers(ii) = this%basis%basis(ii)%atomicNumber
    end do

  end subroutine readBasis


  !> Read in basis function for a species.
  subroutine readSpeciesBasis(node, basisResolution, spBasis)

    !> Node containing the basis definition for a species
    type(fnode), pointer :: node

    !> Grid distance for discretising the basis functions
    real(dp), intent(in) :: basisResolution

    !> Contains the basis on return
    type(TSpeciesBasis), intent(out) :: spBasis

    !> Input node instances, containing the information
    type(fnode), pointer :: tmpNode, child

    !> Node list instance
    type(fnodeList), pointer :: children

    !> Real-valued buffer lists
    type(TListReal) :: bufferExps, bufferCoeffs

    !> Basis coefficients and exponents
    real(dp), allocatable :: coeffs(:), exps(:)

    !> Auxiliary variable
    integer :: ii

    call getChildValue(node, "AtomicNumber", spBasis%atomicNumber)
    call getChildren(node, "Orbital", children)
    spBasis%nOrb = getLength(children)

    if (spBasis%nOrb < 1) then
      call detailedError(node, "Missing orbital definitions")
    end if

    allocate(spBasis%angMoms(spBasis%nOrb))
    allocate(spBasis%occupations(spBasis%nOrb))
    allocate(spBasis%stos(spBasis%nOrb))
    allocate(spBasis%cutoffs(spBasis%nOrb))

    do ii = 1, spBasis%nOrb
      call getItem1(children, ii, tmpNode)
      call getChildValue(tmpNode, "AngularMomentum", spBasis%angMoms(ii))
      call getChildValue(tmpNode, "Occupation", spBasis%occupations(ii))
      call getChildValue(tmpNode, "Cutoff", spBasis%cutoffs(ii))
      call init(bufferExps)

      call getChildValue(tmpNode, "Exponents", bufferExps, child=child)
      if (len(bufferExps) == 0) then
        call detailedError(child, "Missing exponents")
      end if
      call init(bufferCoeffs)
      call getChildValue(tmpNode, "Coefficients", bufferCoeffs, child=child)
      if (len(bufferCoeffs) == 0) then
        call detailedError(child, "Missing coefficients")
      end if
      if (mod(len(bufferCoeffs), len(bufferExps)) /= 0) then
        call detailedError(child, "Number of coefficients incompatible with number of exponents")
      end if
      allocate(exps(len(bufferExps)))
      call asArray(bufferExps, exps)
      call destruct(bufferExps)
      allocate(coeffs(len(bufferCoeffs)))
      call asArray(bufferCoeffs, coeffs)
      call destruct(bufferCoeffs)
      call init(spBasis%stos(ii), reshape(coeffs, [size(coeffs) / size(exps), size(exps)]),&
          & exps, ii - 1, basisResolution, spBasis%cutoffs(ii))
      deallocate(exps)
      deallocate(coeffs)
    end do

  end subroutine readSpeciesBasis


  !> Checks, if the eigenvector file has the right identity number.
  subroutine checkEigenvecs(fileName, identity)

    !> File to check
    character(len=*), intent(in) :: fileName

    !> Identity number
    integer, intent(in) :: identity

    integer :: fd, id, iostat

    open(newunit=fd, file=fileName, action="read", position="rewind", form="unformatted", iostat=iostat)

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


  !> Determinant of a 3x3 matrix (only temporary!)
  function determinant(matrix)

    !> Matrix to calculate the determinant for
    real(dp), intent(in) :: matrix(:,:)

    !> Calculated determinant
    real(dp) :: determinant

    !> Auxiliary variable
    real(dp) :: tmp

    @:ASSERT(all(shape(matrix) == [3, 3]))

    tmp = matrix(1, 1)&
        & * (matrix(2, 2) * matrix(3, 3) - matrix(3, 2) * matrix(2, 3))
    tmp = tmp - matrix(1, 2)&
        & * (matrix(2, 1) * matrix(3, 3) - matrix(3, 1) * matrix(2, 3))
    tmp = tmp + matrix(1, 3)&
        & * (matrix(2, 1) * matrix(3, 2) - matrix(3, 1) * matrix(2, 2))

    determinant = abs(tmp)

  end function determinant

end module waveplot_initwaveplot
