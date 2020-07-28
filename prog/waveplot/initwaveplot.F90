!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the routines for initialising waveplot.
module dftbp_initwaveplot
  use dftbp_assert
  use dftbp_globalenv, only : stdOut
  use dftbp_hsdparser, only : parseHSD, dumpHSD
  use dftbp_xmlutils
  use dftbp_hsdutils
  use dftbp_hsdutils2
  use xmlf90_flib_dom
  use dftbp_linkedlist
  use dftbp_charmanip
  use dftbp_accuracy
  use dftbp_constants
  use dftbp_typegeometryhsd
  use dftbp_message
  use dftbp_fileid
  use dftbp_molecularorbital
  use dftbp_gridcache
  use dftbp_unitconversion
  use dftbp_slater
  implicit none

  private
  save


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

  public :: initProgramVariables

  !! Variables from detailed.xml


  !> Identity of the run
  integer, public :: identity

  !> Geometry
  type(TGeometry), public :: geo

  !> If eigvecs/hamiltonian is real
  logical, public :: tRealHam

  !> Nr. of K-Points
  integer :: nKPoint

  !> Nr. of spins
  integer :: nSpin

  !> Nr. of states
  integer :: nState

  !> Nr. of orbitals per state
  integer, public :: nOrb

  !> K-Points & weights
  real(dp), allocatable :: kPointsWeights(:, :)

  !> Occupations
  real(dp), allocatable, public :: occupations(:,:,:)

  !! Variables from the Option block


  !> If total charge should be plotted
  logical, public :: tPlotTotChrg

  !> If total charge should be calculated
  logical, public :: tCalcTotChrg

  !> If total spin pol. to be plotted
  logical, public :: tPlotTotSpin

  !> If total charge difference to be pl.
  logical, public :: tPlotTotDiff

  !> If atomic densities to be plotted
  logical, public :: tPlotAtomDens

  !> If atomic densities to be calculated
  logical, public :: tCalcAtomDens

  !> If charge for orbitals to be plotted
  logical, public :: tPlotChrg

  !> If chrg difference for orbs. to be pl.
  logical, public :: tPlotChrgDiff

  !> If real part of the wfcs to plot.
  logical, public :: tPlotReal

  !> If imaginary part of the wfcs to plot
  logical, public :: tPlotImag

  !> Levels to plot
  integer, allocatable, public :: plottedLevels(:)

  !> K-points to plot
  integer, allocatable, public :: plottedKPoints(:)

  !> Spins to plot
  integer, allocatable, public :: plottedSpins(:)

  !> Nr. of cached grids
  integer :: nCached

  !> Box vectors for the plotted region
  real(dp), public :: boxVecs(3, 3)

  !> Origin of the box
  real(dp), public :: origin(3)

  !> Origin of the grid in the box
  real(dp), public :: gridOrigin(3)

  !> Nr of grid points along 3 directions
  integer, public :: nPoints(3)

  !> If grid should shifted by a half cell
  logical :: tShiftGrid

  !> If box should filled with folded atoms
  logical, public :: tFillBox

  !> Repeat box along 3 directions
  integer, public :: repeatBox(3)

  !> If coords should be folded to unit cell
  logical, public :: tFoldCoords

  !> If program should be verbose
  logical, public :: tVerbose

  !! Variables from the Basis block


  !> Resolution of the radial wfcs
  real(dp) :: basisResolution

  !> definition of the wfcs
  type(TSpeciesBasis), allocatable, public :: basis(:)

  !! Variables from the EigvecBin block


  !> File with binary eigenvectors
  character(len=1024) :: eigVecBin

  !! Variables from AtomicNumbers block


  !> species-atomic nr. corresp.
  integer, allocatable, public :: atomicNumbers(:)

  !! Locally created variables


  !> Molecular orbital
  type(TMolecularOrbital), allocatable, target, public :: molOrb

  !> pointer to the orbital
  type(TMolecularOrbital), pointer :: pMolOrb

  !> Grid cache
  type(TGridCache), public :: grid

  !> grid vectors
  real(dp), public :: gridVec(3, 3)

  !> List of levels to plot
  integer, allocatable :: levelIndex(:,:)

  !> Volume of the grid
  real(dp), public :: gridVol

contains


  !> Initialise program variables
  subroutine initProgramVariables()

    type(fnode), pointer :: root, tmp, detailed, hsdTree
    type(string) :: strBuffer
    integer :: inputVersion
    integer :: ii
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
    call readDetailed(detailed, tGroundState)
    call destroyNode(tmp)

    !! Read basis
    call getChild(root, "Basis", tmp)
    call readBasis(tmp, geo%speciesNames)
    call getChildValue(root, "EigenvecBin", strBuffer)
    eigVecBin = unquote(char(strBuffer))
    call checkEigenvecs(eigVecBin, identity)

    !! Read options
    call getChild(root, "Options", tmp)
    call readOptions(tmp, identity, nState, nKPoint, nSpin, occupations, &
        &tRealHam, geo%tPeriodic, basis)

    !! Issue warning about unprocessed nodes
    call warnUnprocessedNodes(root, .true.)

    !! Finish parsing, dump parsed and processed input
    call dumpHSD(hsdTree, hsdParsedInput)
    write(stdout, "(A)") "Processed input written as HSD to '" // hsdParsedInput &
        &//"'"
    write(stdout, "(A)") repeat("-", 80)
    write(stdout, *)
    call destroyNode(hsdTree)

    !! Create grid vectors, shift them if necessary
    do ii = 1, 3
      gridVec(:,ii) = boxVecs(:,ii) / real(nPoints(ii), dp)
    end do
    if (tShiftGrid) then
      gridOrigin(:) = origin(:) + 0.5 * sum(gridVec, dim=2)
    else
      gridOrigin(:) = origin(:)
    end if
    gridVol = determinant(gridVec)

    write(stdout, "(A)") "Doing initialisation"

    !! Initialize necessary objects
    allocate(molOrb)
    pMolOrb => molOrb
    call init(molOrb, geo, basis)

    call init(grid, levelIndex=levelIndex, &
        &nOrb=nOrb, nAllLevel=nState, nAllKPoint=nKPoint, nAllSpin=nSpin, &
        &nCached=nCached, nPoints=nPoints, tVerbose=tVerbose, &
        &eigVecBin=eigVecBin, gridVec=gridVec, origin=gridOrigin, &
        &kPointCoords=kPointsWeights(1:3,:), tReal=tRealHam,molorb=pMolOrb)

  end subroutine initProgramVariables


  !> Interpret the information stored in detailed.xml
  subroutine readDetailed(detailed, tGroundState)

    !> Pointer to the node, containing the info
    type(fnode), pointer :: detailed

    !> look for ground state occupations (T) or excited (F)
    logical, intent(in) :: tGroundState

    type(fnode), pointer :: tmp, occ, spin
    integer :: iSpin, iK

    call getChildValue(detailed, "Identity", identity)
    call getChild(detailed, "Geometry", tmp)
    call readGeometry(tmp, geo)
    call getChildValue(detailed, "Real", tRealHam)
    call getChildValue(detailed, "NrOfKPoints", nKPoint)
    call getChildValue(detailed, "NrOfSpins", nSpin)
    call getChildValue(detailed, "NrOfStates", nState)
    call getChildValue(detailed, "NrOfOrbitals", nOrb)
    allocate(kPointsWeights(4, nKPoint))
    call getChildValue(detailed, "KPointsAndWeights", kPointsWeights)
    allocate(occupations(nState, nKPoint, nSpin))

    if (tGroundState) then
      call getChild(detailed, "Occupations", occ)
      do iSpin = 1, nSpin
        call getChild(occ, "spin" // i2c(iSpin), spin)
        do iK = 1, nKPoint
          call getChildValue(spin, "k" // i2c(iK), occupations(:, iK, iSpin))
        end do
      end do
      do iK = 1, nKPoint
        occupations(:, iK, :) = occupations(:, iK, :) * kPointsWeights(4, iK)
      end do
    else
      call getChild(detailed, "ExcitedOccupations", occ)
      do iSpin = 1, nSpin
        call getChild(occ, "spin" // i2c(iSpin), spin)
        do iK = 1, nKPoint
          call getChildValue(spin, "k" // i2c(iK), occupations(:, iK, iSpin))
        end do
      end do
      do iK = 1, nKPoint
        occupations(:, iK, :) = occupations(:, iK, :) * kPointsWeights(4, iK)
      end do
    end if

  end subroutine readDetailed


  !> Read in the geometry stored as xml in internal or gen format.
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


  !> Interpret the options.
  subroutine readOptions(node, identity, nLevel, nKPoint, nSpin, occupations, tRealHam, tPeriodic, &
      & basis)

    !> Node containig the information
    type(fnode), pointer :: node

    !> Identity nr. of the dftb+ calculation to after-process
    integer, intent(in) :: identity

    !> Nr. of states in the dftb+ calc.
    integer, intent(in) :: nLevel

    !> Nr. of K-points
    integer, intent(in) :: nKPoint

    !> Nr. of spins
    integer, intent(in) :: nSpin

    !> Occupations in the dftb+ calculation
    real(dp), intent(in) :: occupations(:,:,:)

    !> If Hamiltonian was real or not
    logical, intent(in) :: tRealHam

    !> If the system was periodic
    logical, intent(in) :: tPeriodic

    !> Basis functions (needed for the cutoffs)
    type(TSpeciesBasis), intent(in) :: basis(:)

    type(fnode), pointer :: subnode, field, value
    type(string) :: buffer, modifier
    type(TListIntR1) :: indexBuffer
    integer :: curId
    integer :: ind, ii, iLevel, iKPoint, iSpin, iAtom, iSpecies
    logical :: tFound
    real(dp) :: tmpvec(3), minvals(3), maxvals(3)
    real(dp), allocatable :: mcutoffs(:)
    real(dp) :: minEdge
    character(len=63) :: warnId(3) = (/ &
        &"The external files you are providing differ from those provided", &
        &"when this input file was generated. The results you obtain with", &
        &"the current files could be, therefore, different.              " &
        &/)

    !! Warning, if processed input is read in, but eigenvectors are different
    call getChildValue(node, "Identity", curId, identity)
    if (curId /= identity) then
      call warning(warnId)
    end if

    call getChildValue(node, "TotalChargeDensity", tPlotTotChrg, .false.)
    if (nSpin == 2) then
      call getChildValue(node, "TotalSpinPolarisation", tPlotTotSpin, .false.)
    else
      tPlotTotSpin = .false.
    end if
    call getChildValue(node, "TotalChargeDifference", tPlotTotDiff, .false., &
        &child=field)
    call getChildValue(node, "TotalAtomicDensity", tPlotAtomDens, .false.)
    call getChildValue(node, "ChargeDensity", tPlotChrg, .false.)
    call getChildValue(node, "ChargeDifference", tPlotChrgDiff, .false.)
    tCalcTotChrg = tPlotTotChrg .or. tPlotTotSpin .or. tPlotTotDiff&
        &.or. tPlotChrgDiff
    tCalcAtomDens = tPlotTotDiff .or. tPlotChrgDiff .or. tPlotAtomDens
    call getChildValue(node, "RealComponent", tPlotReal, .false.)
    call getChildValue(node, "ImagComponent", tPlotImag, .false., child=field)
    if (tPlotImag .and. tRealHam) then
      call detailedWarning(field, &
          &"Wave functions are real, no imaginary part will be plotted")
      tPlotImag = .false.
    end if
    call getChildValue(node, "PlottedLevels", buffer, child=field, &
        &multiple=.true.)
    call convRangeToInt(char(buffer), node, plottedLevels, nLevel)
    if (tPeriodic) then
      call getChildValue(node, "PlottedKPoints", buffer, child=field, &
          &multiple=.true.)
      call convRangeToInt(char(buffer), node, plottedKPoints, nKPoint)
    else
      allocate(plottedKPoints(1))
      plottedKPoints(1) = 1
    end if
    call getChildValue(node, "PlottedSpins", buffer, child=field, &
        &multiple=.true.)
    call convRangeToInt(char(buffer), node, plottedSpins, nSpin)

    !! Create the list of the levels, which must be calculated explicitely
    call init(indexBuffer)
    do iSpin = 1, nSpin
      do iKPoint = 1, nKPoint
        do iLevel = 1, nLevel
          tFound = any(plottedLevels == iLevel) &
              &.and. any(plottedKPoints == iKPoint) &
              &.and. any(plottedSpins == iSpin)
          if ((.not. tFound) .and. tCalcTotChrg) then
            tFound = occupations(iLevel, iKPoint, iSpin) > 1e-8_dp
          end if
          if (tFound) then
            call append(indexBuffer, (/iLevel, iKPoint, iSpin /))
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

    call getChildValue(node, "NrOfCachedGrids", nCached, 1, child=field)
    if (nCached < 1 .and. nCached /= -1) then
      call detailedError(field, "Value must be -1 or greater than zero.")
    end if
    if (nCached == -1) then
      nCached = size(levelIndex, dim=2)
    end if

    !! Plotted region: if last (and hopefully only) childnode is not an allowed method -> assume
    !! explicit setting, parse the node "PlottedRegion" for the appropriate children.
    call getChildValue(node, "PlottedRegion", value, child=subnode)
    call getNodeName(value, buffer)
    select case (char(buffer))
    case ("unitcell")
      !! Unit cell for the periodic case, smallest possible cuboid for cluster
      if (tPeriodic) then
        origin(:) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
        boxVecs(:,:) = geo%latVecs(:,:)
      else
        call getChildValue(value, "MinEdgeLength", minEdge, child=field, &
            &default=1.0_dp)
        if (minEdge < 0.0_dp) then
          call detailedError(field, "Minimal edge length must be positive")
        end if
        origin = minval(geo%coords, dim=2)
        tmpvec = maxval(geo%coords, dim=2) - origin
        do ii = 1, 3
          if (tmpvec(ii) < minEdge) then
            origin(ii) = origin(ii) - 0.5_dp * (minEdge - tmpvec(ii))
            tmpvec(ii) = minEdge
          end if
        end do
        boxVecs(:,:) = 0.0_dp
        do ii = 1, 3
          boxVecs(ii, ii) = tmpvec(ii)
        end do
      end if

    case ("optimalcuboid")
      !! Determine optimal cuboid, so that no basis function leaks out
      call getChildValue(value, "MinEdgeLength", minEdge, child=field, &
          &default=1.0_dp)
      if (minEdge < 0.0_dp) then
        call detailedError(field, "Minimal edge length must be positive")
      end if
      allocate(mcutoffs(geo%nSpecies))
      do iSpecies = 1 , geo%nSpecies
        mcutoffs(iSpecies) = maxval(basis(iSpecies)%cutoffs)
      end do
      minvals = geo%coords(:,1)
      maxvals = geo%coords(:,1)
      do iAtom = 1, geo%nAtom
        iSpecies = geo%species(iAtom)
        maxvals(:) = max(maxvals, geo%coords(:, iAtom) + mcutoffs(iSpecies))
        minvals(:) = min(minvals, geo%coords(:, iAtom) - mcutoffs(iSpecies))
      end do
      origin(:) = minvals(:)
      tmpvec(:) = maxvals(:) - minvals(:)
      do ii = 1, 3
        if (tmpvec(ii) < minEdge) then
          origin(ii) = origin(ii) - 0.5_dp * (minEdge - tmpvec(ii))
          tmpvec(ii) = minEdge
        end if
      end do
      boxVecs(:,:) = 0.0_dp
      do ii = 1, 3
        boxVecs(ii, ii) = tmpvec(ii)
      end do

    case ("origin","box")
      !! Those nodes are part of an explicit specification -> explitic specif
      call getChildValue(subnode, "Box", boxVecs, modifier=modifier, &
          &child=field)
      if (abs(determinant(boxVecs)) < 1e-8_dp) then
        call detailedError(field, "Vectors are linearly dependent")
      end if
      if (len(modifier) > 0) then
        ind = getModifierIndex(char(modifier), lengthUnits, field)
        boxVecs(:,:) = boxVecs(:,:) * lengthUnits(ind)%convertValue
      end if
      call getChildValue(subnode, "Origin", origin, modifier=modifier, &
          &child=field)
      if (len(modifier) > 0) then
        ind = getModifierIndex(char(modifier), lengthUnits, field)
        origin(:) = origin(:) * lengthUnits(ind)%convertValue
      end if

    case default
      !! Object with unknown name passed
      call detailedError(value, "Invalid element name")
    end select

    !! Replace existing PlottedRegion definition
    call setChild(node, "PlottedRegion", field, replace=.true.)
    call setChildValue(field, "Origin", origin, .true.)
    call setChildValue(field, "Box", boxVecs, .true.)

    call getChildValue(node, "NrOfPoints", nPoints, child=field)
    if (any(nPoints <= 0)) then
      call detailedError(field, "Specified numbers must be greater than zero")
    end if
    call getChildValue(node, "ShiftGrid", tShiftGrid, default=.true.)
    if (tPeriodic) then
      call getChildValue(node, "FoldAtomsToUnitCell", tFoldCoords, &
          &default=.false.)
      call getChildValue(node, "FillBoxWithAtoms", tFillBox, default=.false.)
      tFoldCoords = tFoldCoords .or. tFillBox
    else
      tFillBox = .false.
      tFoldCoords = .false.
    end if
    call getChildValue(node, "RepeatBox", repeatBox, default=(/ 1, 1, 1 /), &
        &child=field)
    if (.not. all(repeatBox > 0)) then
      call detailedError(field, "Indexes must be greater than zero")
    end if
    call getChildValue(node, "Verbose", tVerbose, .false.)

  end subroutine readOptions


  !> Read in the basis related informations
  subroutine readBasis(node, speciesNames)

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

    call getChildValue(node, "Resolution", basisResolution)
    allocate(basis(nSpecies))
    allocate(atomicNumbers(nSpecies))
    do ii = 1, nSpecies
      speciesName = speciesNames(ii)
      call getChild(node, speciesName, speciesNode)
      call readSpeciesBasis(speciesNode, basisResolution, basis(ii))
      atomicNumbers(ii) = basis(ii)%atomicNumber
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

    type(fnode), pointer :: tmpNode, child
    type(fnodeList), pointer :: children
    type(TListReal) :: bufferExps, bufferCoeffs
    real(dp), allocatable :: coeffs(:), exps(:)
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
        call detailedError(child, "Number of coefficients incompatible with &
            &number of exponents")
      end if
      allocate(exps(len(bufferExps)))
      call asArray(bufferExps, exps)
      call destruct(bufferExps)
      allocate(coeffs(len(bufferCoeffs)))
      call asArray(bufferCoeffs, coeffs)
      call destruct(bufferCoeffs)
      call init(spBasis%stos(ii), &
          &reshape(coeffs, (/ size(coeffs)/size(exps), size(exps) /)), &
          &exps, ii - 1, basisResolution, spBasis%cutoffs(ii))
      deallocate(exps)
      deallocate(coeffs)
    end do

  end subroutine readSpeciesBasis


  !> Checks, if the eigenvector file has the right identity number.
  subroutine checkEigenvecs(fileName, identity)

    !> File to check
    character(len=*), intent(in) :: fileName

    !> Identity number.
    integer, intent(in) :: identity

    integer :: fd, id, iostat

    fd = getFileId()
    open(fd, file=fileName, action="read", position="rewind", &
        &form="unformatted", iostat=iostat)
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

end module dftbp_initwaveplot
