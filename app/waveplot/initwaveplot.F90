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
      & getNodeName, destroyNode
  use dftbp_io_charmanip, only : i2c, unquote
  use dftbp_io_hsdparser, only : parseHSD, dumpHSD
  use dftbp_io_hsdutils, only : getChildValue, setChildValue, getChild, setChild, getChildren,&
      & getSelectedIndices, detailedError, detailedWarning
  use dftbp_io_hsdutils2, only : convertUnitHsd, readHSDAsXML, warnUnprocessedNodes
  use dftbp_io_message, only : warning, error
  use dftbp_io_xmlutils, only : removeChildNodes
  use dftbp_math_simplealgebra, only : invert33
  use dftbp_type_linkedlist, only : TListIntR1, TListReal, init, destruct, len, append, asArray
  use dftbp_type_typegeometryhsd, only : TGeometry, readTGeometryGen, readTGeometryHSD,&
      & readTGeometryVasp, readTGeometryXyz, writeTGeometryHSD
  use dftbp_io_charmanip, only : tolower

  use waveplot_molorb, only : TMolecularOrbital, init, TSpeciesBasis
  use waveplot_grids, only : modifyEigenvecs, TTabulationTypesEnum, TGridInterpolationTypesEnum
  use waveplot_slater, only : initSto
  use omp_lib

  implicit none

  private
  save

  public :: TProgramVariables
  public :: TProgramVariables_init


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

  !> Container for enumerated grid interpolation types
  type(TGridInterpolationTypesEnum), parameter :: gridInterpolTypes = TGridInterpolationTypesEnum()

  !> Container for enumerated radial wf tabulation types
  type(TTabulationTypesEnum), parameter :: rwTabulationTypes = TTabulationTypesEnum()


  !> Data type containing variables from detailed.xml
  type TInput

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

    !> k-points and weights
    real(dp), allocatable :: kPointsandWeight(:,:)

  end type TInput


  !> Data type containing variables from eigenvec.bin
  type TEig

    !> Nr. of states
    integer :: nState

    !> Real eigenvectors
    real(dp), allocatable :: eigvecsReal(:,:,:,:)

    !> Complex eigenvectors
    complex(dp), allocatable :: eigvecsCplx(:,:,:,:)

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

    !> If total spin polarization to be plotted
    logical :: tPlotTotSpin

    !> If total charge difference to be plotted
    logical :: tPlotTotDiff

    !> If atomic densities to be plotted
    logical :: tPlotAtomDens

    !> If atomic densities to be calculated
    logical :: tCalcAtomDens

    !> If charge for orbitals to be plotted
    logical :: tPlotChrg

    !> If charge difference for orbitals to be plotted
    logical :: tPlotChrgDiff

    !> If real part of the wfcs to plot
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

    !> Type for radial wavefunction tabulation.
    integer :: rwTabulationType

    !> Interpolation type for grid interpolation.
    integer :: gridInterType

    !> NUmber of parallel regions which should be used in the mapping from subgrid to the total grid
    integer :: parallelRegionNum

  end type TOption


  !> Data type containing variables from the Basis block
  type TBasis

    !> Definition of the wfcs
    type(TSpeciesBasis), allocatable :: basis(:)

    !> Resolution of the radial wfcs
    real(dp) :: basisResolution

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

    !> Index mapping: orbital --> atom
    integer, allocatable :: orbitalToAtom(:)

    !> Index mapping: orbital --> species
    integer, allocatable :: orbitalToSpecies(:)

    !> Index mapping: orbital --> angular momentum
    integer, allocatable :: orbitalToAngMoms(:)

    !> Index mapping: orbital --> magnetic quantum number
    integer, allocatable :: orbitalToM(:)

    !> Index mapping: orbital --> slater type orbital
    integer, allocatable :: orbitalToStos(:)

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

    !> If grid should shifted by a half of a grid vector
    logical :: tShiftGrid

    logical :: tGroundState

    !> Array to temporary store the eigenvectors from the eigenvec.bin file
    real(dp), allocatable :: tmparray(:,:)
    complex(dp), allocatable :: tmparrayCplx(:,:)

    !> True, if the radial WFs should be calculated explicitly, false, if they should be calculated
    !> via interpolation
    logical :: rwExplicit
    type(string) :: buffer
    character(len=:), allocatable :: charTabulation

    !> Operation status, if an error needs to be returned
    type(TStatus) :: errStatus

    ! Write header
    write(stdout, "(A)") repeat("=", 80)
    write(stdout, "(A)") "     WAVEPLOT  " // version
    write(stdout, "(A,/)") repeat("=", 80)

    ! Read in input file as HSD
    call parseHSD(rootTag, hsdInput, hsdTree)
    call getChild(hsdTree, rootTag, root)

    write(stdout, "(A)") "Interpreting input file '" // hsdInput // "'"

    ! Check if input version is the one, which we can handle
    call getChildValue(root, "InputVersion", inputVersion, parserVersion)
    if (inputVersion /= parserVersion) then
      call error("Version of input (" // i2c(inputVersion) // ") and parser (" &
          &// i2c(parserVersion) // ") do not match")
    end if

    call getChildValue(root, "GroundState", tGroundState, .true.)

    ! Read data from detailed.xml
    call getChildValue(root, "DetailedXML", strBuffer)
    call readHSDAsXML(unquote(char(strBuffer)), tmp)
    call getChild(tmp, "detailedout", detailed)
    call readDetailed(this, detailed, tGroundState, nKPoint, nSpin, this%eig%nState)
    call destroyNode(tmp)

    ! Read the type of tabulation of the radial WFs
    call getChild(root, "Options", tmp)
    call getChildValue(tmp, "RadialWFTabulation", Buffer, "explicit")

    charTabulation = tolower(unquote(char(Buffer)))

    if (.not. (charTabulation == 'polynomial' .or. charTabulation == 'trivial' &
        & .or. charTabulation == 'linear' .or. charTabulation == 'explicit')) &
        & call detailedError(tmp, 'Wrong wavefunction interpolation type specified')

    select case (charTabulation)

      case ("trivial")
        this%opt%rwTabulationType = rwTabulationTypes%trivial

      case ("linear")
        this%opt%rwTabulationType = rwTabulationTypes%linear

      case ("polynomial")
        this%opt%rwTabulationType = rwTabulationTypes%polynomial

      case ("explicit")
        this%opt%rwTabulationType = rwTabulationTypes%explicit

      case default
        this%opt%rwTabulationType = rwTabulationTypes%explicit

    end select

    if (this%opt%rwTabulationType == rwTabulationTypes%explicit) then
      rwExplicit = .true.
    else
      rwExplicit = .false.
    end if

    ! Read basis
    call getChild(root, "Basis", tmp)
    call readBasis(this, tmp,  this%input%geo%speciesNames, rwExplicit)
    call getChildValue(root, "EigenvecBin", strBuffer)
    eigVecBin = unquote(char(strBuffer))
    call checkEigenvecs(eigVecBin, this%input%identity)

    ! Read options
    call getChild(root, "Options", tmp)
    call readOptions(this, tmp, this%eig%nState, nKPoint, nSpin, tShiftGrid)

    ! Read eigenvectors
    if (this%input%tRealHam) then
      allocate(tmparray(this%input%nOrb, this%eig%nState * nKPoint * nSpin))
      allocate(this%eig%eigvecsReal(this%input%nOrb, this%eig%nState, 1, nSpin))
      call readEigenvecs(eigVecBin, tmparray)
      call modifyEigenvecs(tmparray, this%eig%eigvecsReal, 1, nSpin)
    else
      allocate(tmparrayCplx(this%input%nOrb, this%eig%nState * nKPoint * nSpin))
      allocate(this%eig%eigvecsCplx(this%input%nOrb, this%eig%nState, nKPoint, nSpin))
      call readEigenvecs(eigVecBin, tmparrayCplx)
      call modifyEigenvecs(tmparrayCplx, this%eig%eigvecsCplx, nKPoint, nSpin)
    end if

    ! Issue warning about unprocessed nodes
    call warnUnprocessedNodes(root, .true.)

    ! Finish parsing, dump parsed and processed input
    call dumpHSD(hsdTree, hsdParsedInput)

    write(stdout, "(A)") "Processed input written as HSD to '" // hsdParsedInput //"'"
    write(stdout, "(A)") repeat("-", 80)
    write(stdout, *)

    call destroyNode(hsdTree)

    write(stdout, "(A)") "Doing initialisation"
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

    ! Initialize necessary (molecular orbital, grid) objects
    allocate(this%loc%molOrb)
    this%loc%pMolOrb => this%loc%molOrb
    call init(this%loc%molOrb, this%input%geo, this%boundaryCond, this%basis%basis)

    ! Determine useful index mappings
    call getIndexMappings(this)

    ! Repeat boxes if necessary
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

    !> K-Points & weights
    real(dp), allocatable :: kPointsWeights(:,:)

    type(fnode), pointer :: tmp, occ, spin

    integer :: iSpin, iK

    call getChildValue(detailed, "Identity", this%input%identity)
    call getChild(detailed, "Geometry", tmp)
    call readGeometry(this, tmp)

    call getChildValue(detailed, "Real", this%input%tRealHam)
    call getChildValue(detailed, "NrOfKPoints", nKPoint)
    call getChildValue(detailed, "NrOfSpins", nSpin)
    call getChildValue(detailed, "NrOfStates", nState)
    call getChildValue(detailed, "NrOfOrbitals", this%input%nOrb)

    allocate(kPointsWeights(4, nKPoint))
    allocate(this%input%occupations(nState, nKPoint, nSpin))

    call getChildValue(detailed, "KPointsAndWeights", kPointsWeights)

    allocate(this%input%kPointsandWeight(4, nKPoint))
    this%input%kPointsandWeight = kPointsWeights

    if (tGroundState) then

      call getChild(detailed, "Occupations", occ)
      do iSpin = 1, nSpin
        call getChild(occ, "spin" // i2c(iSpin), spin)
        do iK = 1, nKPoint
          call getChildValue(spin, "k" // i2c(iK), this%input%occupations(:, iK, iSpin))
        end do
      end do
      do iK = 1, nKPoint
        this%input%occupations(:,iK,:) = this%input%occupations(:,iK,:) * kPointsWeights(4,iK)
      end do

    else

      call getChild(detailed, "ExcitedOccupations", occ)
      do iSpin = 1, nSpin
        call getChild(occ, "spin" // i2c(iSpin), spin)
        do iK = 1, nKPoint
          call getChildValue(spin, "k" // i2c(iK), this%input%occupations(:, iK, iSpin))
        end do
      end do
      do iK = 1, nKPoint
        this%input%occupations(:,iK,:) = this%input%occupations(:,iK,:) * kPointsWeights(4,iK)
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
      call readTGeometryGen(child, this%input%geo)
      call removeChildNodes(geonode)
      call writeTGeometryHSD(geonode, this%input%geo)

    case ("xyzformat")
      call readTGeometryXyz(child, this%input%geo)
      call removeChildNodes(geonode)
      call writeTGeometryHSD(geonode, this%input%geo)

    case ("vaspformat")
      call readTGeometryVasp(child, this%input%geo)
      call removeChildNodes(geonode)
      call writeTGeometryHSD(geonode, this%input%geo)

    case default
      call readTGeometryHSD(geonode, this%input%geo)

    end select

  end subroutine readGeometry


  !> Interpret the options.
  subroutine readOptions(this, node, nLevel, nKPoint, nSpin, tShiftGrid)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> Node containig the information
    type(fnode), pointer :: node

    !> Nr. of states in the dftb+ calculation
    integer, intent(in) :: nLevel

    !> Nr. of K-points
    integer, intent(in) :: nKPoint

    !> Nr. of spins
    integer, intent(in) :: nSpin

    !> If grid should shifted by a half of a grid vector
    logical, intent(out) :: tShiftGrid

    !> List of levels to plot, whereby insignificant occupations were filtered out
    integer, allocatable :: levelIndex(:,:)

    type(fnode), pointer :: subnode, field, value
    type(string) :: buffer, modifier
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
    integer :: ind, ii, iLevel, iKPoint, iSpin, iAtom, iSpecies
    integer :: curVec(3)
    real(dp) :: tmpvec(3), minvals(3), maxvals(3)
    real(dp), allocatable :: mcutoffs(:)
    real(dp) :: minEdge
    integer :: numThreadsDefault
    character(len=:), allocatable :: charInterpol

    ! Warning, if processed input is read in, but eigenvectors are different
    call getChildValue(node, "Identity", curId, this%input%identity)
    if (curId /= this%input%identity) then
      call warning(warnId)
    end if

    call getChildValue(node, "TotalChargeDensity", this%opt%tPlotTotChrg, .false.)

    call getChildValue(node, "GridInterpolation", Buffer, "linear")
    charInterpol = tolower(unquote(char(Buffer)))

    if (.not. (charInterpol == 'linear' .or. charInterpol == 'trivial' .or. &
        & charInterpol == 'explicit')) &
        & call detailedError(node, 'Wrong grid interpolation type specified')

    select case (charInterpol)

      case ('trivial')
        this%opt%gridInterType = gridInterpolTypes%trivial

      case ('linear')
        this%opt%gridInterType = gridInterpolTypes%linear

      case ('explicit')
        this%opt%gridInterType = gridInterpolTypes%explicit

      case default
        this%opt%gridInterType = gridInterpolTypes%linear

    end select

    if (nSpin == 2) then
      call getChildValue(node, "TotalSpinPolarisation", this%opt%tPlotTotSpin, .false.)
    else
      this%opt%tPlotTotSpin = .false.
    end if

    call getChildValue(node, "TotalChargeDifference", this%opt%tPlotTotDiff, .false.,&
        &child=field)
    call getChildValue(node, "TotalAtomicDensity", this%opt%tPlotAtomDens, .false.)
    call getChildValue(node, "ChargeDensity", this%opt%tPlotChrg, .false.)
    call getChildValue(node, "ChargeDifference", this%opt%tPlotChrgDiff, .false.)

    this%opt%tCalcTotChrg = this%opt%tPlotTotChrg .or. this%opt%tPlotTotSpin .or.&
        & this%opt%tPlotTotDiff
    this%opt%tCalcAtomDens = this%opt%tPlotTotDiff .or. this%opt%tPlotChrgDiff .or.&
        & this%opt%tPlotAtomDens

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

    ! Create the list of the levels, which must be calculated explicitely
    call init(indexBuffer)

    do iSpin = 1, nSpin
      do iKPoint = 1, nKPoint
        do iLevel = 1, nLevel
          tFound = any(this%opt%plottedLevels == iLevel) &
              &.and. any(this%opt%plottedKPoints == iKPoint) &
              &.and. any(this%opt%plottedSpins == iSpin)
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

    allocate(levelIndex(3, len(indexBuffer)))
    call asArray(indexBuffer, levelIndex)
    call destruct(indexBuffer)

    allocate(this%opt%levelIndex(3, size(levelIndex, dim=2)))
    this%opt%levelIndex = levelIndex

    call getChildValue(node, "SpGridPoints", this%opt%nSpPoints, child=field)
    if (any(this%opt%nSpPoints <= 0)) then
      call detailedError(field, "Specified numbers must be greater than zero")
    end if

    numThreadsDefault = omp_get_max_threads()

    call getChildValue(node, "ParallelRegionNum", this%opt%parallelRegionNum,&
        & default=numThreadsDefault, child=field)
    if (this%opt%parallelRegionNum <= 0) then
      call detailedError(field, "Specified numbers must be greater than zero")
    end if

    ! Plotted region: if last (and hopefully only) childnode is not an allowed method -> assume
    ! explicit setting, parse the node "PlottedRegion" for the appropriate children.
    call getChildValue(node, "PlottedRegion", value, child=subnode)
    call getNodeName(value, buffer)

    allocate(this%loc%speciesGridsOrigs(3, this%input%geo%nSpecies))

    select case (char(buffer))

    case ("unitcell")
      ! Unit cell for the periodic case, smallest possible cuboid for cluster
      if (this%input%geo%tPeriodic) then
        this%opt%origin(:) = [0.0_dp, 0.0_dp, 0.0_dp]
        this%opt%boxVecs(:,:) = this%input%geo%latVecs
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
          this%opt%boxVecs(ii,ii) = tmpvec(ii)
        end do
      end if

    case ("optimalcuboid")
      ! Determine optimal cuboid, so that no basis function leaks out
      call getChildValue(value, "MinEdgeLength", minEdge, child=field, default=1.0_dp)
      if (minEdge < 0.0_dp) then
        call detailedError(field, "Minimal edge length must be positive")
      end if
      allocate(mcutoffs(this%input%geo%nSpecies))
      do iSpecies = 1, this%input%geo%nSpecies
        mcutoffs(iSpecies) = maxval(this%basis%basis(iSpecies)%cutoffs)
        this%loc%speciesGridsOrigs(:, iSpecies) = - mcutoffs(iSpecies)
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
      allocate(this%loc%speciesGridsVecs(3, 3, this%input%geo%nSpecies))
      this%loc%speciesGridsVecs(:,:,:) = 0.0_dp
      do iSpecies = 1, this%input%geo%nSpecies
        do ii = 1, 3
          this%loc%speciesGridsVecs(ii, ii, iSpecies) = 2.0_dp * mcutoffs(iSpecies)&
              & / (real(this%opt%nSpPoints(ii), dp))
        end do
      end do
      do ii = 1, 3
        this%opt%boxVecs(ii,ii) = tmpvec(ii)
      end do

    case ("origin","box")
      ! Those nodes are part of an explicit specification -> explitic specif
      call getChildValue(subnode, "Box", this%opt%boxVecs, modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), lengthUnits, field, this%opt%boxVecs)
      if (abs(determinant(this%opt%boxVecs)) < 1e-08_dp) then
        call detailedError(field, "Vectors are linearly dependent")
      end if
      call getChildValue(subnode, "Origin", this%opt%origin, modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), lengthUnits, field, this%opt%origin)

    case default
      ! Object with unknown name passed
      call detailedError(value, "Invalid element name")

    end select

    ! Replace existing PlottedRegion definition
    call setChild(node, "PlottedRegion", field, replace=.true.)
    call setChildValue(field, "Origin", this%opt%origin, .true.)
    call setChildValue(field, "Box", this%opt%boxVecs, .true.)

    call getChildValue(node, "TotGridPoints", this%opt%nTotPoints, child=field)

    if (any(this%opt%nTotPoints <= 0)) then
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

    call getChildValue(node, "Verbose", this%opt%tVerbose, .false.)

    ! Create grid vectors, shift them if necessary
    do ii = 1, 3
      this%loc%totGridVec(:, ii) = this%opt%boxVecs(:,ii)&
          & / real(this%opt%nTotPoints(ii), dp)
    end do

    if (char(buffer) == 'origin' .or. char(buffer) == 'box' .or. char(buffer) == 'unitcell') then
      call calculateSpeciesGridsVecs(this)
    end if

    ! If desired, shift the origin by half of a grid point
    if (tShiftGrid) then
      this%opt%totGridOrig(:) = this%opt%origin(:) +&
          & 0.5_dp * sum(this%loc%totGridVec, dim=2)
      do iSpecies = 1, this%input%geo%nSpecies
        this%loc%speciesGridsOrigs(:, iSpecies) =&
            & this%loc%speciesGridsOrigs(:, iSpecies) + 0.5_dp *&
            & sum(this%loc%totGridVec, dim=2)
      end do
    else
      this%opt%totGridOrig(:) = this%opt%origin(:)
    end if
    this%loc%gridVol = determinant(this%loc%totGridVec)

    ! Calculate the center coordinate of the total grid
    this%loc%totGridCenter(:) = 0.0_dp
    do ii = 1, 3
      this%loc%totGridCenter(:) = this%loc%totGridCenter + 0.5_dp *&
          & real(this%opt%nTotPoints(ii), dp) * this%loc%totGridVec(:,ii)
    end do

  end subroutine readOptions


  !> Calculates box vectors of species grids for given total grid vectors
  subroutine calculateSpeciesGridsVecs(this)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> cartesian total grid coordinates
    real(dp), allocatable :: subcubeCartCoords(:,:)

    !> real total grid coordinates
    real(dp) :: gridrealCoords(3,8)

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

    allocate(this%loc%speciesGridsVecs(3, 3, this%input%geo%nSpecies))

    call invert33(totGridInvBasis, this%loc%totGridVec)

    do iSpecies = 1, this%input%geo%nSpecies

      maxSpeciesCutoff = maxval(this%basis%basis(iSpecies)%cutoffs)
      subcubeCartCoords = reshape([&
          & 0.0_dp, 0.0_dp, 0.0_dp,&
          & 0.0_dp, 0.0_dp, 1.0_dp,&
          & 0.0_dp, 1.0_dp, 0.0_dp,&
          & 0.0_dp, 1.0_dp, 1.0_dp,&
          & 1.0_dp, 0.0_dp, 0.0_dp,&
          & 1.0_dp, 0.0_dp, 1.0_dp,&
          & 1.0_dp, 1.0_dp, 0.0_dp,&
          & 1.0_dp, 1.0_dp, 1.0_dp],&
          & [3, 8]) * 2.0_dp * maxSpeciesCutoff

      do ii = 1, size(subcubeCartCoords, dim=2)
        gridrealCoords(:,ii) = matmul(totGridInvBasis, subcubeCartCoords(:, ii))
      end do

      maxRanges = ceiling(maxval(gridrealCoords, dim=2)) - floor(minval(gridrealCoords, dim=2))

      do ii = 1, 3
        speciesBoxVecs(:,ii) = this%loc%totGridVec(:,ii) * maxRanges(ii)
      end do

      this%loc%speciesGridsOrigs(:, iSpecies) = [0.0_dp, 0.0_dp, 0.0_dp]
      do ii = 1, 3
        this%loc%speciesGridsVecs(:, ii, iSpecies) = speciesBoxVecs(:, ii)&
            & / real(this%opt%nSpPoints(ii), dp)
        this%loc%speciesGridsOrigs(:, iSpecies) = this%loc%speciesGridsOrigs(:, iSpecies)&
            & - 0.5_dp * real(this%opt%nSpPoints(ii), dp) *&
            & this%loc%speciesGridsVecs(:, ii, iSpecies)
      end do

    end do

  end subroutine calculateSpeciesGridsVecs


  !> Read in the basis related informations.
  subroutine readBasis(this, node, speciesNames, rwExplicit)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> Node containing the basis definition
    type(fnode), pointer :: node

    !> Names of the species for which basis should be read in
    character(len=*), intent(in) :: speciesNames(:)

    !> True, if the radial WFs should be calculated explicitly, false, if they should be calculated
    !> via interpolation
    logical, intent(in) :: rwExplicit

    character(len=len(speciesNames)) :: speciesName
    type(fnode), pointer :: speciesNode
    integer :: nSpecies
    integer :: ii

    nSpecies = size(speciesNames)

    @:ASSERT(nSpecies > 0)

    if (rwExplicit .eqv. .true.) then
      call getChildValue(node, "Resolution", this%basis%basisResolution)
    end if

    allocate(this%basis%basis(nSpecies))
    allocate(this%aNr%atomicNumbers(nSpecies))

    do ii = 1, nSpecies
      speciesName = speciesNames(ii)
      call getChild(node, speciesName, speciesNode)
      call readSpeciesBasis(speciesNode, this%basis%basis(ii), &
          & rwExplicit, this%basis%basisResolution)
      this%aNr%atomicNumbers(ii) = this%basis%basis(ii)%atomicNumber
    end do

  end subroutine readBasis


  !> Read in basis function for a species.
  subroutine readSpeciesBasis(node, spBasis, rwExplicit, basisResolution)

    !> Node containing the basis definition for a species
    type(fnode), pointer :: node

    !> Contains the basis on return
    type(TSpeciesBasis), intent(out) :: spBasis

    !> Routine which should be used to tabulate the radial WFs.
    logical, intent(in) :: rwExplicit

    !> Grid distance for discretising the basis functions
    real(dp), intent(in), optional :: basisResolution

    type(fnode), pointer :: tmpNode, child
    type(fnodeList), pointer :: children

    !> Real-valued buffer lists
    type(TListReal) :: bufferExps, bufferCoeffs

    !> Basis coefficients and exponents
    real(dp), allocatable :: coeffs(:), exps(:)

    type(TListReal) :: bufferRwf
    type(string) :: buffer
    character(len=10) :: charBuffer
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
    allocate(spBasis%stos(spBasis%nOrb))
    allocate(spBasis%cutoffs(spBasis%nOrb))

    do ii = 1, spBasis%nOrb

      call getItem1(children, ii, tmpNode)
      call getChildValue(tmpNode, "AngularMomentum", spBasis%angMoms(ii))
      call getChildValue(tmpNode, "Occupation", spBasis%occupations(ii))
      call getChildValue(tmpNode, "Cutoff", spBasis%cutoffs(ii))

      if (.not. (rwExplicit)) then
        call init(bufferRwf)
        call getChildValue(tmpNode, "RadialWaveFunc", bufferRwf, child=child)

        if (len(bufferRwf) == 0) then
          call detailedError(child, "Missing radial wavefunction")
        end if

        allocate(rwf(len(bufferRwf)))

        call asArray(bufferRwf, rwf)
        call destruct(bufferRwf)

        call initSto(spBasis%stos(ii), reshape([0.0_dp, 0.0_dp], [1, 2]),&
            & [0.0_dp, 0.0_dp], ii - 1, spBasis%cutoffs(ii), &
            & rwf=transpose(reshape(rwf, [3, size(rwf) / 3])))

        deallocate(rwf)
      end if

      if (rwExplicit) then
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
        exps(:) = 0.0_dp
        call asArray(bufferExps, exps)
        call destruct(bufferExps)

        allocate(coeffs(len(bufferCoeffs)))
        coeffs = 0.0_dp
        call asArray(bufferCoeffs, coeffs)
        call destruct(bufferCoeffs)

        call initSto(spBasis%stos(ii), reshape(coeffs, [size(coeffs) / size(exps), size(exps)]),&
            & exps, ii - 1, spBasis%cutoffs(ii), basisResolution)

        deallocate(exps)
        deallocate(coeffs)
      end if

    end do

  end subroutine readSpeciesBasis


  !> Determines convenient index mappings
  subroutine getIndexMappings(this)

    !> Container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> Auxiliary arrays
    integer, allocatable :: tmp1(:), tmp2(:), tmp3(:,:)

    !> Auxiliary variables
    integer :: iAtom, iOrb, iSpecies, iAng, iL, iM, iOrbPerSpecies, mAng, ind, magnQN, &
        & angularMoment

    allocate(this%loc%orbitalOcc(this%input%nOrb, 1))
    allocate(this%loc%orbitalToAtom(this%input%nOrb))
    allocate(this%loc%orbitalToSpecies(this%input%nOrb))
    allocate(tmp1(this%input%geo%nSpecies + 1))
    allocate(tmp2(this%input%nOrb))
    allocate(tmp3(2, this%input%nOrb))
    allocate(this%loc%orbitalToAngMoms(this%input%nOrb))
    allocate(this%loc%orbitalToM(this%input%nOrb))
    allocate(this%loc%orbitalToStos(this%input%nOrb))

    iOrb = 1

    do iAtom = 1, this%input%geo%nAtom
      ind = 1
      iSpecies = this%input%geo%species(iAtom)
      do iAng = 1, size(this%basis%basis(iSpecies)%angMoms)
        mAng = 2 * this%basis%basis(iSpecies)%angMoms(iAng) + 1
        this%loc%orbitalOcc(iOrb:iOrb + mAng - 1, 1) =&
            & this%basis%basis(iSpecies)%occupations(iAng) / real(mAng, dp)
        do iM = 1, mAng
          tmp3(1, iOrb) = iSpecies
          tmp3(2, iOrb) = ind
          this%loc%orbitalToAtom(iOrb) = iAtom
          iOrb = iOrb + 1
          ind = ind + 1
        end do
      end do
    end do

    ind = 1

    do iAtom = 1, this%input%geo%nAtom
      iSpecies = this%input%geo%species(iAtom)
      do iAng = 1, size(this%basis%basis(iSpecies)%angMoms)
        mAng = 2 * this%basis%basis(iSpecies)%angMoms(iAng) + 1
        this%loc%orbitalOcc(ind:ind + mAng - 1, 1) =&
            & this%basis%basis(iSpecies)%occupations(iAng) / real(mAng, dp)
        ind = ind + mAng
      end do
    end do

    ind = 1

    do iSpecies = 1, this%input%geo%nSpecies
      tmp1(iSpecies) = ind
      mAng = maxval(this%basis%basis(iSpecies)%angMoms)
      do iAng = this%loc%molorb%iStos(iSpecies), this%loc%molorb%iStos(iSpecies + 1) - 1
        iL = this%loc%molorb%angMoms(iAng)
        ind = ind + 2 * iL + 1
      end do
    end do

    tmp1(this%input%geo%nSpecies + 1) = ind

    do iOrb = 1, this%input%nOrb
      iSpecies = tmp3(1, iOrb)
      iOrbPerSpecies = tmp3(2, iOrb)
      this%loc%orbitalToSpecies(iOrb) = tmp1(iSpecies) + iOrbPerSpecies - 1
    end do

    iOrb = 1
    do iAtom = 1, this%input%geo%nAtom
      iSpecies = this%input%geo%species(iAtom)
      do angularMoment = 1, size(this%basis%basis(iSpecies)%angMoms)
        iL = this%basis%basis(iSpecies)%angMoms(angularMoment)
        do magnQN = - iL, iL
          this%loc%orbitalToM(iOrb) = magnQN
          this%loc%orbitalToAngMoms(iOrb) = iL
          this%loc%orbitalToStos(iOrb) = iSpecies
          iOrb = iOrb + 1
        end do
      end do
    end do

    iOrb = 1
    do iAtom = 1, this%input%geo%nAtom
      iSpecies = this%input%geo%species(iAtom)
      do iAng = this%loc%molorb%iStos(iSpecies), this%loc%molorb%iStos(iSpecies + 1) - 1
        iL = this%loc%molorb%angMoms(iAng)
        do magnQN = - iL, iL
          this%loc%orbitalToStos(iOrb) = iAng
          iOrb = iOrb + 1
        end do
        ind = ind + 1
      end do
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

    nBox = product(this%opt%repeatBox)

    if (nBox > 1) then

      ! If tFillBox is off, coordinates must be repeated here.
      ! Otherwise the part for filling with atoms will do that.
      if (.not. this%opt%tFillBox) then

        allocate(tmpCoords(3, size(this%input%geo%coords)))
        allocate(tmpSpecies(size(this%input%geo%species)))

        ! temporarily save unrepeated coordinates and species
        tmpCoords(:,:) = this%input%geo%coords(:,:)
        tmpSpecies(:) = this%input%geo%species(:)

        deallocate(this%input%geo%coords)
        deallocate(this%input%geo%species)
        allocate(this%input%geo%coords(3, nBox * this%input%geo%nAtom))
        allocate(this%input%geo%species(nBox * this%input%geo%nAtom))

        ind = 0

        do i1 = 0, this%opt%repeatBox(1) - 1
          do i2 = 0, this%opt%repeatBox(2) - 1
            do i3 = 0, this%opt%repeatBox(3) - 1
              shift(:) = matmul(this%opt%boxVecs, real([i1, i2, i3], dp))
              do iAtom = 1, this%input%geo%nAtom
                this%input%geo%coords(:, ind + iAtom) = tmpCoords(:, iAtom) + shift(:)
              end do
              this%input%geo%species(ind + 1:ind + this%input%geo%nAtom) = tmpSpecies(:)
              ind = ind + this%input%geo%nAtom
            end do
          end do
        end do
        this%input%geo%nAtom = nBox * this%input%geo%nAtom
      end if
      do i1 = 1, 3
        this%opt%boxVecs(:,i1) = this%opt%boxVecs(:,i1) * real(this%opt%repeatBox(i1), dp)
      end do

    end if

  end subroutine getRepeatedBox


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

    nOrb = size(eigenvecs, dim=2)

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


  !> Checks, if the eigenvector file has the right identity number.
  subroutine checkEigenvecs(fileName, identity)

    !> File to check
    character(len=*), intent(in) :: fileName

    !> Identity number.
    integer, intent(in) :: identity

    integer :: fd, id, iostat

    open(newunit=fd, file=fileName, action="read", position="rewind", form="unformatted",&
        & iostat=iostat)

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
    real(dp), intent(in) :: matrix(:,:)

    !> Determinant of the matrix.
    real(dp) :: determinant

    real(dp) :: tmp

    @:ASSERT(all(shape(matrix) == (/3, 3/)))

    tmp = matrix(1, 1) * (matrix(2, 2) * matrix(3, 3) - matrix(3, 2) * matrix(2, 3))
    tmp = tmp - matrix(1, 2) * (matrix(2, 1) * matrix(3, 3) - matrix(3, 1) * matrix(2, 3))
    tmp = tmp + matrix(1, 3) * (matrix(2, 1) * matrix(3, 2) - matrix(3, 1) * matrix(2, 2))

    determinant = abs(tmp)

  end function determinant

end module waveplot_initwaveplot
