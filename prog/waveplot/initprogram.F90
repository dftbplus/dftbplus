!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains the routines for initialising waveplot.
module InitProgram
  use assert
  use HSDParser, only : parseHSD, dumpHSD, dumpHSDAsXML
  use XMLUtils
  use HSDUtils
  use HSDUtils2
  use flib_dom
  use LinkedList
  use CharManip
  use Accuracy
  use Constants
  use TypeGeometryHSD
  use Message
  use FileId
  use MolecularOrbital
  use GridCache
  use UnitConversion
  use Slater
  implicit none

  private
  save

  character(len=*), parameter :: version =  "0.3"
  character(len=*), parameter :: rootTag = "waveplot"
  character(len=*), parameter :: hsdInput = "waveplot_in.hsd"
  character(len=*), parameter :: hsdParsedInput = "waveplot_pin.hsd"
  character(len=*), parameter :: xmlInput = "waveplot_in.xml"
  character(len=*), parameter :: xmlParsedInput = "waveplot_pin.xml"
  integer, parameter :: parserVersion = 3

  public :: initProgramVariables


  !! Variables from detailed.xml
  integer, public :: identity          ! Identity of the run
  type(TGeometry), public :: geo       ! Geometry
  logical, public :: tRealHam          ! If eigvecs/hamiltonian is real
  integer :: nKPoint                   ! Nr. of K-Points
  integer :: nSpin                     ! Nr. of spins
  integer :: nState                    ! Nr. of states
  integer, public :: nOrb              ! Nr. of orbitals per state
  real(dp), allocatable :: kPointsWeights(:, :)         ! K-Points & weights
  real(dp), allocatable, public :: occupations(:,:,:)   ! Occupations

  !! Variables from the Option block
  logical, public :: tPlotTotChrg      ! If total charge should be plotted
  logical, public :: tCalcTotChrg      ! If total charge should be calculated
  logical, public :: tPlotTotSpin      ! If total spin pol. to be plotted
  logical, public :: tPlotTotDiff      ! If total charge difference to be pl.
  logical, public :: tPlotAtomDens     ! If atomic densities to be plotted
  logical, public :: tCalcAtomDens     ! If atomic densities to be calculated
  logical, public :: tPlotChrg         ! If charge for orbitals to be plotted
  logical, public :: tPlotChrgDiff     ! If chrg difference for orbs. to be pl.
  logical, public :: tPlotReal         ! If real part of the wfcs to plot.
  logical, public :: tPlotImag         ! If imaginary part of the wfcs to plot
  integer, allocatable, public :: plottedLevels(:)   ! Levels to plot
  integer, allocatable, public :: plottedKPoints(:)  ! K-points to plot
  integer, allocatable, public :: plottedSpins(:)    ! Spins to plot
  integer  :: nCached                  ! Nr. of cached grids
  real(dp), public :: boxVecs(3, 3)    ! Box vectors for the plotted region
  real(dp), public :: origin(3)        ! Origin of the box
  real(dp), public :: gridOrigin(3)    ! Origin of the grid in the box
  integer, public :: nPoints(3)        ! Nr of grid points along 3 directions
  logical :: tShiftGrid                ! If grid should shifted by a half cell
  logical, public :: tFillBox          ! If box should filled with folded atoms
  integer, public :: repeatBox(3)      ! Repeat box along 3 directions
  logical, public :: tFoldCoords       ! If coords should be folded to unit cell
  logical, public :: tVerbose          ! If program should be verbose

  !! Variables from the Basis block
  real(dp) :: basisResolution                    ! Resolution of the radial wfcs
  type(TSpeciesBasis), allocatable, public :: basis(:)   ! definition of the wfcs

  !! Variables from the EigvecBin block
  character(len=1024) :: eigVecBin     ! File with binary eigenvectors

  !! Variables from AtomicNumbers block
  integer, allocatable, public :: atomicNumbers(:)  ! species-atomic nr. corresp.

  !! Locally created variables
  type(OMolecularOrbital), allocatable, target, public :: molOrb    ! Molecular orbital
  type(OMolecularOrbital), pointer :: pMolOrb
  type(OGridCache), public :: grid                     ! Grid cache
  real(dp), public :: gridVec(3, 3)                    ! grid vectors
  integer, allocatable :: levelIndex(:,:)              ! List of levels to plot
  real(dp), public :: gridVol                          ! Volume of the grid

contains

  !!* Initialise program variables
  subroutine initProgramVariables()

    type(fnode), pointer :: input, root, tmp, detailed
    type(string) :: strBuffer
    integer :: inputVersion
    integer :: ii
    logical :: tHSD, tGroundState

    !! Write header
    write (*, "(A)") repeat("=", 80)
    write (*, "(A)") "     WAVEPLOT  " // version
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
    call dumpHSD(input, hsdParsedInput)
    write (*, "(A)") "Processed input written as HSD to '" // hsdParsedInput &
        &//"'"
    call dumpHSDAsXML(input, xmlParsedInput)
    write (*, "(A)") "Processed input written as XML to '" // xmlParsedInput &
        &//"'"
    write (*, "(A)") repeat("-", 80)
    write (*,*)
    call destroyNode(input)

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

    write (*, "(A)") "Doing initialisation"

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


  !!* Interpret the information stored in detailed.xml
  !!* @param detailed Pointer to the node, containing the info
  !!* @para tGroundState look for ground state occupations (T) or excited (F)
  subroutine readDetailed(detailed, tGroundState)
    type(fnode), pointer :: detailed
    logical, intent(in)  :: tGroundState

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



  !!* Interpret the options.
  !!* @param node        Node containig the information
  !!* @param identity    Identity nr. of the dftb+ calculation to after-process
  !!* @param nLevel      Nr. of states in the dftb+ calc.
  !!* @param nKPoint     Nr. of K-points
  !!* @param nSpin       Nr. of spins
  !!* @param occupations Occupations in the dftb+ calculation
  !!* @param tRealHam    If Hamiltonian was real or not
  !!* @param tPeriodic   If the system was periodic
  !!* @param basis       Basis functions (needed for the cutoffs)
  subroutine readOptions(node, identity, nLevel, nKPoint, nSpin, occupations, &
      &tRealHam, tPeriodic, basis)
    type(fnode), pointer :: node
    integer, intent(in) :: identity
    integer, intent(in) :: nLevel
    integer, intent(in) :: nKPoint
    integer, intent(in) :: nSpin
    real(dp), intent(in) :: occupations(:,:,:)
    logical, intent(in) :: tRealHam
    logical, intent(in) :: tPeriodic
    type(TSpeciesBasis), intent(in) :: basis(:)

    type(fnode), pointer :: subnode, field, value
    type(string) :: buffer, modifier
    type(ListIntR1) :: indexBuffer
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

    !! Plotted region: if last (and hopefully only) childnode is not an allowed
    !! method -> assume explicit setting, parse the node "PlottedRegion" for
    !! the appropriate children.
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
        boxVecs(:,:) = boxVecs(:,:) * lengthUnits(ind)%value
      end if
      call getChildValue(subnode, "Origin", origin, modifier=modifier, &
          &child=field)
      if (len(modifier) > 0) then
        ind = getModifierIndex(char(modifier), lengthUnits, field)
        origin(:) = origin(:) * lengthUnits(ind)%value
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



  !!* Read in the basis related informations
  !!* @param node        Node containing the basis definition
  !!* @param speciesNames Names of the species for which basis should be read in
  subroutine readBasis(node, speciesNames)
    type(fnode), pointer :: node
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



  !!* Read in basis function for a species.
  !!* @param node            Node containing the basis definition for a species
  !!* @param basisResolution Grid distance for discretising the basis functions
  !!* @param spBasis         Contains the basis on return
  subroutine readSpeciesBasis(node, basisResolution, spBasis)
    type(fnode), pointer :: node
    real(dp), intent(in) :: basisResolution
    type(TSpeciesBasis), intent(out) :: spBasis

    type(fnode), pointer :: tmpNode, child
    type(fnodeList), pointer :: children
    type(listReal) :: bufferExps, bufferCoeffs
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



  !!* Checks, if the eigenvector file has the right identity number.
  !!* @param fileName File to check
  !!* @param identity Identity number.
  subroutine checkEigenvecs(fileName, identity)
    character(len=*), intent(in) :: fileName
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
      call error("Ids for eigenvectors ("// i2c(id) //") and xml-input ("&
          &// i2c(identity) // ") don't match.")
    end if
    close(fd)

  end subroutine checkEigenvecs



  !!* Determinant of a 3x3 matrix (Only temporary!)
  !!* @param matrix The matrix to calculate the determinant from.
  !!* @return       Determinant of the matrix.
  real(dp) function  determinant(matrix)
    real(dp), intent(in) :: matrix(:, :)

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



end module InitProgram
