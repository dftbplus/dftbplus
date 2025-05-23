!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Routines to read/write a TGeometry type in HSD and XML format.
module dftbp_type_typegeometryhsd
  use dftbp_common_accuracy, only : dp, lc, mc
  use dftbp_common_atomicmass, only : getAtomicSymbol
  use dftbp_common_constants, only : AA__Bohr, avogadConst, Bohr__AA, pi
  use dftbp_common_globalenv, only : stdout
  use dftbp_common_unitconversion, only : angularUnits, lengthUnits
  use dftbp_extlibs_xmlf90, only : char, fnode, getNodeType, getNodeValue,&
      & flib_normalize => normalize, string, TEXT_NODE, xmlf_t
  use dftbp_io_charmanip, only : i2c, tolower
  use dftbp_io_hsdutils, only : checkError, detailedError, detailedWarning, getChildValue,&
      & getFirstTextChild, setChildValue, writeChildValue
  use dftbp_io_hsdutils2, only : convertUnitHsd, splitModifier
  use dftbp_io_message, only : error
  use dftbp_io_tokenreader, only : getNextToken, TOKEN_ERROR, TOKEN_OK
  use dftbp_math_simplealgebra, only : determinant33, invert33
  use dftbp_type_linkedlist, only : append, asArray, destruct, find, init, len, TListIntR1,&
      & TListRealR1, TListString
  use dftbp_type_typegeometry, only : normalize, reduce, setlattice, TGeometry
  implicit none

  private
  !> Types/subroutines from TypeGeometry
  public :: TGeometry, normalize
  !> Locally defined subroutines
  public :: writeTGeometryHSD, readTGeometryHSD, readTGeometryGen
  public :: readTGeometryXyz, readTGeometryVasp, readTGeometryLammps
  !> makes public subroutines from typegeometry
  public :: reduce, setlattice

  !> Writes the content of a geometry object to a dom tree or to an xml-writer
  interface writeTGeometryHSD
    module procedure writeTGeometryHSD_dom
    module procedure writeTGeometryHSD_xmlf
  end interface

contains


  !> Write the geometry in HSD format to a specified node
  subroutine writeTGeometryHSD_dom(node, geo)

    !> Node in the HSD-tree which should contain the geometry
    type(fnode), pointer :: node

    !> The geometry
    type(TGeometry), intent(in) :: geo

    call setChildValue(node, "TypeNames", geo%speciesNames, .false.)
    call setChildValue(node, "TypesAndCoordinates", &
        &reshape(geo%species, (/ 1, size(geo%species) /)), geo%coords, .false.)
    call setChildValue(node, "Periodic", geo%tPeriodic, .false.)
    if (geo%tPeriodic .or. geo%tHelical) then
      call setChildValue(node, "LatticeVectors", geo%latVecs, .false.)
      call setChildValue(node, "CoordinateOrigin", geo%origin, .false.)
    end if
    call setChildValue(node, "Helical", geo%tHelical, .false.)

  end subroutine writeTGeometryHSD_dom


  !> Write the geometry in HSD format to an xml writer
  subroutine writeTGeometryHSD_xmlf(xf, geo)

    !> Node in the HSD-tree which should contain the geometry
    type(xmlf_t), intent(inout) :: xf

    !> The geometry
    type(TGeometry), intent(in) :: geo

    call writeChildValue(xf, "TypeNames", geo%speciesNames)
    if (geo%tPeriodic .or. geo%tHelical) then
      call writeChildValue(xf, "TypesAndCoordinates", &
          &reshape(geo%species, (/ 1, size(geo%species) /)), geo%coords&
          & + spread(geo%origin, 2, size(geo%species)))
    else
      call writeChildValue(xf, "TypesAndCoordinates", &
          &reshape(geo%species, (/ 1, size(geo%species) /)), geo%coords)
    end if
    call writeChildValue(xf, "Periodic", geo%tPeriodic)
    call writeChildValue(xf, "Helical", geo%tHelical)
    if (geo%tPeriodic .or. geo%tHelical) then
      call writeChildValue(xf, "CoordinateOrigin", geo%origin)
      call writeChildValue(xf, "LatticeVectors", geo%latVecs)
    end if

  end subroutine writeTGeometryHSD_xmlf


  !> Read the geometry from a node in a HSD tree.
  subroutine readTGeometryHSD(node, geo)

    !> Node in the HSD tree containing the geomery
    type(fnode), pointer :: node

    !> Contains the geometry on exit
    type(TGeometry), intent(out) :: geo

    type(string) :: modifier, modifs(2)
    type(TListString) :: stringBuffer
    type(TListRealR1) :: realBuffer
    type(TListIntR1) :: intBuffer
    type(fnode), pointer :: child, typesAndCoords
    integer, allocatable :: tmpInt(:,:)
    real(dp) :: latvec(9), det, helVec(3)

    call getChildValue(node, "Periodic", geo%tPeriodic, default=.false.)
    call getChildValue(node, "Helical", geo%tHelical, default=.false.)
    if (geo%tPeriodic .and. geo%tHelical) then
      call error("Periodic and helical boundary conditions mutually exclusive.")
    end if
    call init(stringBuffer)
    call getChildValue(node, "TypeNames", stringBuffer)
    geo%nSpecies = len(stringBuffer)
    if (geo%nSpecies == 0) then
      call detailedError(node, "Missing species names.")
    end if
    allocate(geo%speciesNames(geo%nSpecies))
    call asArray(stringBuffer, geo%speciesNames)
    call destruct(stringBuffer)
    call init(intBuffer)
    call init(realBuffer)
    call getChildValue(node, "TypesAndCoordinates", 1, intBuffer, 3, &
        &realBuffer, modifier=modifier, child=typesAndCoords)
    geo%nAtom = len(intBuffer)
    if (geo%nAtom == 0) then
      call detailedError(typesAndCoords, "Missing coordinates")
    end if
    allocate(geo%species(geo%nAtom))
    allocate(geo%coords(3, geo%nAtom))
    allocate(tmpInt(1, geo%nAtom))
    call asArray(intBuffer, tmpInt)
    call destruct(intBuffer)
    geo%species(:) = tmpInt(1,:)
    deallocate(tmpInt)
    !! Check validity of species
    if (any(geo%species < 1 .or. geo%species > geo%nSpecies)) then
      call detailedError(typesAndCoords, "Type index must be between 1 and " &
          &// i2c(geo%nSpecies) // ".")
    end if
    call asArray(realBuffer, geo%coords)
    call destruct(realBuffer)
    geo%tFracCoord = .false.
    if (len(modifier) > 0) then
      select case(tolower(char(modifier)))
      case ("relative")
        if (.not. geo%tPeriodic) then
          call detailedError(typesAndCoords, "Relative coordinates are only &
              &allowed for periodic systems")
        end if
        geo%tFracCoord = .true.
      case default
        call convertUnitHsd(char(modifier), lengthUnits, typesAndCoords, geo%coords)
        call setChildValue(typesAndCoords, "", reshape(geo%species, [1, size(geo%species)]),&
            & geo%coords, replace=.true.)
      end select
    end if
    if (geo%tPeriodic) then
      allocate(geo%origin(3))
      if (geo%tFracCoord) then
        call getChildValue(node, "CoordinateOrigin", geo%origin, [0.0_dp,0.0_dp,0.0_dp])
      else
        call getChildValue(node, "CoordinateOrigin", geo%origin, [0.0_dp,0.0_dp,0.0_dp],&
            & modifier=modifier, child=child)
        call convertUnitHsd(char(modifier), lengthUnits, child, geo%origin, replace=.true.)
      end if
      geo%coords(:,:) = geo%coords - spread(geo%origin, 2, geo%nAtom)
      allocate(geo%latVecs(3,3))
      call getChildValue(node, "LatticeVectors", latvec, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), lengthUnits, child, latvec, replace=.true.)
      geo%latVecs(:,:) = reshape(latvec, [3, 3])
      if (geo%tFracCoord) then
        geo%coords(:,:) = matmul(geo%latVecs, geo%coords)
        geo%origin(:) = matmul(geo%latVecs, geo%origin)
      end if
      allocate(geo%recVecs2p(3, 3))
      det = determinant33(geo%latVecs)
      if (abs(det) < 1e-12_dp) then
        call detailedError(child, "Dependent lattice vectors")
      end if
      call invert33(geo%recVecs2p, geo%latVecs, det)
      geo%recVecs2p(:,:) = reshape(geo%recVecs2p, (/3, 3/), order=(/2, 1/))
    end if

    if (geo%tHelical) then
      allocate(geo%origin(3))
      call getChildValue(node, "CoordinateOrigin", geo%origin, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), lengthUnits, child, geo%origin, replace=.true.)
      geo%coords(:,:) = geo%coords - spread(geo%origin, 2, geo%nAtom)
      allocate(geo%latVecs(3, 1))
      call getChildValue(node, "LatticeVectors", helVec, modifier=modifier, child=child)
      if (len(modifier) > 0) then
        call splitModifier(char(modifier), child, modifs)
        call convertUnitHsd(char(modifs(1)), lengthUnits, child, helVec(1))
        call convertUnitHsd(char(modifs(2)), angularUnits, child, helVec(2))
      end if
      geo%latVecs(:3,1) = helVec
      if (geo%latVecs(3,1) < 1) then
        call error("Helical structure rotation order non-positive")
      end if
      allocate(geo%recVecs2p(1, 1))
      geo%recVecs2p = 2.0_dp * pi / geo%latVecs(1,1)
    end if

    call normalize(geo)

  end subroutine readTGeometryHSD


  !> Reads the geometry from a node in a HSD tree in GEN format
  subroutine readTGeometryGen(node, geo)

    !> Node containing the geometry in Gen format
    type(fnode), pointer :: node

    !> Contains the geometry on exit
    type(TGeometry), intent(out) :: geo

    type(string) :: text

    call getFirstTextChild(node, text)
    call readTGeometryGen_help(node, geo, char(text))

  end subroutine readTGeometryGen


  !> Helping routine for reading geometry from a HSD tree in GEN format
  subroutine readTGeometryGen_help(node, geo, text)

    !> Node to parse (only needed to produce proper error messages)
    type(fnode), pointer :: node

    !> Contains the geometry on exit
    type(TGeometry), intent(out) :: geo

    !> Text content of the node
    character(len=*), intent(in) :: text

    type(string) :: txt
    integer :: iStart, iErr, iEnd
    integer :: ii, iTmp, iSp
    real(dp) :: coords(3)
    type(TListString) :: speciesNames
    character(lc) :: errorStr

    ! Read first line of the gen file: Number of atoms, boundary conditions
    iStart = 1
    iEnd = nextLine(text, iStart)
    call getNextToken(text(:iEnd), geo%nAtom, iStart, iErr)
    call checkError(node, iErr, "Bad number of atoms in the first line of geometry")
    call getNextToken(text(:iEnd), txt, iStart, iErr)
    select case (char(txt))
    case("S","s")
      geo%tPeriodic = .true.
      geo%tFracCoord = .false.
      geo%tHelical = .false.
    case("F","f")
      geo%tPeriodic = .true.
      geo%tFracCoord = .true.
      geo%tHelical = .false.
    case("C", "c")
      geo%tPeriodic = .false.
      geo%tFracCoord = .false.
      geo%tHelical = .false.
    case("H", "h")
      geo%tHelical = .true.
      geo%tPeriodic = .false.
      geo%tFracCoord = .false.
    case default
      call detailedError(node, "Unknown boundary condition type '" &
          &// char(txt) // "'")
    end select
    if (iStart < iEnd) then
      call detailedError(node, "Found trailing characters in first line")
    end if
    iStart = iEnd + 1

    ! Reading the 2nd line of a gen file.
    iEnd = nextLine(text, iStart)
    call init(speciesNames)
    iErr = TOKEN_OK
    iSp = 0
    do while (iErr == TOKEN_OK)
      call getNextToken(text(:iEnd), txt, iStart, iErr)
      if (iErr == TOKEN_OK) then
        if (find(speciesNames, char(txt)) > 0) then
          call detailedError(node, "Species name '"//char(txt)//"' is not unique, check species'//&
              &//' names in second line")
        end if
        call append(speciesNames, char(txt))
      end if
    end do
    geo%nSpecies = len(speciesNames)
    if (geo%nSpecies == 0) then
      call detailedError(node, "No species are given in second line of geometry.")
    end if
    allocate(geo%speciesNames(geo%nSpecies))
    call asArray(speciesNames, geo%speciesNames)
    call destruct(speciesNames)
    iStart = iEnd + 1

    ! Read in sequential and species indices.
    allocate(geo%species(geo%nAtom))
    allocate(geo%coords(3, geo%nAtom))
    do ii = 1, geo%nAtom
      ! save atom number as string for error printout
      write(errorStr, '(i0)') ii
      iEnd = nextLine(text, iStart)
      call getNextToken(text(:iEnd), iTmp, iStart, iErr)
      call checkError(node, iErr, "Bad sequential number for atom "//trim(errorStr))
      call getNextToken(text(:iEnd), geo%species(ii), iStart, iErr)
      call checkError(node, iErr, "Bad species number for atom "//trim(errorStr))
      call getNextToken(text(:iEnd), coords, iStart, iErr)
      call checkError(node, iErr, "Bad coordinates for atom "//trim(errorStr))
      geo%coords(:, ii) = coords(:)
      if (iStart < iEnd) then
        call detailedError(node, "Found trailing characters for atom "//trim(errorStr))
      end if
      iStart = iEnd + 1
    end do
    if (geo%nSpecies /= maxval(geo%species) .or. minval(geo%species) /= 1) then
      call detailedError(node, &
          &"Nr. of species and nr. of specified elements do not match.")
    end if

    ! Read in origin an lattice vectors, if the structure is periodic
    if (geo%tPeriodic) then
      iEnd = nextLine(text, iStart)
      allocate(geo%origin(3))
      allocate(geo%latVecs(3, 3))
      call getNextToken(text(:iEnd), geo%origin, iStart, iErr)
      call checkError(node, iErr, "Invalid origin given in geometry.")
      iStart = iEnd + 1
      do ii = 1, 3
        iEnd = nextLine(text, iStart)
        call getNextToken(text(:iEnd), geo%latVecs(:, ii), iStart, iErr)
        call checkError(node, iErr, "Invalid lattice vectors in geometry.")
        iStart = iEnd + 1
      end do
    end if

    if (geo%tHelical) then
      allocate(geo%origin(3))
      iEnd = nextLine(text, iStart)
      call getNextToken(text(:iEnd), geo%origin, iStart, iErr)
      call checkError(node, iErr, 'Invalid specified helical boundary conditions: origin.')
      geo%origin(:) = geo%origin * AA__Bohr
      allocate(geo%latVecs(3, 1))
      iEnd = nextLine(text, iStart)
      call getNextToken(text(:iEnd), geo%latVecs(:, 1), iStart, iErr)
      call checkError(node, iErr, 'Invalid specified helical boundary conditions: "translation,&
          & twist angle, rotation order" should be supplied).')
      geo%latVecs(1,1) = geo%latVecs(1,1) * AA__Bohr
      geo%latVecs(2,1) = geo%latVecs(2,1) * pi / 180.0_dp
      if (geo%latVecs(3,1) < 1) then
        call error("Helical structure rotation order non-positive")
      end if
      allocate(geo%recVecs2p(1, 1))
      geo%recVecs2p = 1.0_dp / (geo%latVecs(1,1) * 2.0_dp * pi)
    end if

    ! Check if any data remains in the geometry - should be nothing left now
    if (iStart <= len(text)) then
      call detailedError(node, "Superfluous data found. Check if specified number of atoms matches&
          & the number of actually entered positions.")
    end if

    ! tests that are relevant to periodic geometries only
    if (geo%tPeriodic) then
      call setupPeriodicGeometry(node, geo)
    end if

    ! convert coords to correct internal units
    if (geo%tFracCoord) then
      geo%coords(:,:) = matmul(geo%latVecs, geo%coords)
      geo%origin(:) = matmul(geo%latVecs, geo%origin)
    else
      geo%coords = geo%coords * AA__Bohr
    end if

    if (geo%tHelical .or. geo%tPeriodic) then
      geo%coords(:,:) = geo%coords - spread(geo%origin, 2, geo%nAtom)
    end if

    call normalize(geo)

  end subroutine readTGeometryGen_help


  !> Reads the geometry from a node in a HSD tree in XYZ format
  subroutine readTGeometryXyz(node, geo)

    !> Node containing the geometry in XYZ format
    type(fnode), pointer :: node

    !> Contains the geometry on exit
    type(TGeometry), intent(out) :: geo

    type(string) :: text

    call getFirstTextChild(node, text)
    call readTGeometryXyz_help(node, geo, char(text))

  end subroutine readTGeometryXyz


  !> Helping routine for reading geometry from a HSD tree in XYZ format
  subroutine readTGeometryXyz_help(node, geo, text)

    !> Node to parse (only needed to produce proper error messages)
    type(fnode), pointer :: node

    !> Contains the geometry on exit
    type(TGeometry), intent(out) :: geo

    !> Text content of the node
    character(len=*), intent(in) :: text

    type(string) :: txt
    integer :: iStart, iOldStart, iErr, iEnd
    integer :: ii, iSp
    real(dp) :: coords(3)
    type(TListString) :: speciesNames
    character(lc) :: errorStr

    ! Read first line of the xyz file: Number of atoms
    iStart = 1
    iEnd = nextLine(text, iStart)
    call getNextToken(text(:iEnd), geo%nAtom, iStart, iErr)
    call checkError(node, iErr, "Bad number of atoms.")

    if (iStart < iEnd) then
      call checkError(node, TOKEN_ERROR, "Additional field(s) found on first line of xyz geometry")
    end if

    ! advance to next line
    iStart = iEnd + 1

    ! The parser can strip empty comment lines, so we have to try to reconstruct
    ! the original xyz file first...
    iOldStart = iStart
    iEnd = nextLine(text, iStart)

    ! check `second' line
    call getNextToken(text(:iEnd), txt, iStart, iErr)
    if (iErr == TOKEN_OK) then
      call getNextToken(text(:iEnd), coords, iStart, iErr)
    end if
    if (iErr == TOKEN_OK) then
      iStart = iOldStart  ! second line was empty or HSD commented and therefore stripped
    else
      iStart = iEnd + 1  ! second line is an actual XYZ comment line, drop it
    end if

    ! Read in sequential and species indices.
    call init(speciesNames)
    allocate(geo%species(geo%nAtom))
    allocate(geo%coords(3, geo%nAtom))
    iSp = 0
    do ii = 1, geo%nAtom
      iEnd = nextLine(text, iStart)
      call getNextToken(text(:iEnd), txt, iStart, iErr)
      write(errorStr,"(A,1X,I0)")"Bad species name for atom", ii
      call checkError(node, iErr, trim(errorStr))
      iSp = find(speciesNames, char(txt))
      if (iSp == 0) then
        call append(speciesNames, char(txt))
        iSp = len(speciesNames)
      end if
      geo%species(ii) = iSp
      call getNextToken(text(:iEnd), coords, iStart, iErr)
      write(errorStr,"(A,1X,I0)")"Bad coordinates for atom", ii
      call checkError(node, iErr, trim(errorStr))
      geo%coords(:, ii) = coords(:)
      iStart = iEnd + 1
    end do

    geo%nSpecies = len(speciesNames)
    allocate(geo%speciesNames(geo%nSpecies))
    call asArray(speciesNames, geo%speciesNames)
    call destruct(speciesNames)

    if (geo%nSpecies /= maxval(geo%species) .or. minval(geo%species) /= 1) then
      call detailedError(node, "Nr. of species and nr. of specified elements do not match.")
    end if

    if (iStart <= len(text)) then
      call detailedError(node, "Superfluous data found. Check if specified number of atoms matches&
          & the number of actually entered positions.")
    end if

    ! convert coords to correct internal units
    geo%coords = geo%coords * AA__Bohr

    ! original xyz files are always molecular boundary conditions
    geo%tPeriodic = .false.
    geo%tFracCoord = .false.
    geo%tHelical = .false.

    call normalize(geo)

  end subroutine readTGeometryXyz_help


  !> Reads the geometry from a node in a HSD tree in VASP POSCAR/CONTCAR formats
  subroutine readTGeometryVasp(node, geo)

    !> Node containing the geometry in Gen format
    type(fnode), pointer :: node

    !> Contains the geometry on exit
    type(TGeometry), intent(out) :: geo

    type(string) :: text

    call getFirstTextChild(node, text)
    call readTGeometryVasp_help(node, geo, char(text))

  end subroutine readTGeometryVasp


  !> Helping routine for reading geometry from a HSD tree in VASP format
  subroutine readTGeometryVasp_help(node, geo, text)

    !> Node to parse (only needed to produce proper error messages)
    type(fnode), pointer :: node

    !> Contains the geometry on exit
    type(TGeometry), intent(out) :: geo

    !> Text content of the node
    character(len=*), intent(in) :: text

    type(string) :: txt
    character(mc), allocatable :: vaspNames(:)
    integer :: iStart, iOldStart, iErr, iEnd
    integer :: ii, iSp, iTmp
    real(dp) :: coords(3), latVec(3), rScale
    integer, allocatable :: vaspSp(:)
    integer, allocatable :: countSp(:)
    type(TListString) :: speciesNames
    logical :: hasComment
    character(lc) :: errorStr


    ! Read `first' line of the POSCAR/CONTCAR file
    ! This is actually a comment line, but contains by user convention the atomic symbols
    ! In case it does not contain the atomic symbols, is empty or an HSD comment,
    ! we have to ignore it, therefore we attempt to read it as `second' line first
    hasComment = .false.
    iStart = 1
    iEnd = nextLine(text, iStart)
    call getNextToken(text(:iEnd), rScale, iStart, iErr)
    ! seems like we found the `second' line already, to be sure we check for another token
    if (iErr == TOKEN_OK) then
      hasComment = iStart <= iEnd
    else
      hasComment = .true.
    end if

    if (hasComment) then
      iStart = 1
      call init(speciesNames)
      iErr = TOKEN_OK
      do while(iErr == TOKEN_OK)
        call getNextToken(text(:iEnd), txt, iStart, iErr)
        if (iErr == TOKEN_OK) then
          call append(speciesNames, char(txt))
        end if
      end do
      iStart = iEnd + 1

      ! try to read the real `second' line now
      iEnd = nextLine(text, iStart)
      call getNextToken(text(:iEnd), rScale, iStart, iErr)
      call checkError(node, iErr, "Bad scaling factor in line 2 of VASP geometry")
      iStart = iEnd + 1
    end if

    if (rScale <= 0.0_dp) then
      call detailedError(node, "Scaling factor for VASP geometry must be positive and non-zero")
    end if

    ! we expect the lattice information now
    allocate(geo%origin(3))
    allocate(geo%latVecs(3, 3))
    geo%origin(:) = 0.0_dp
    do ii = 1, 3
      iEnd = nextLine(text, iStart)
      call getNextToken(text, latVec, iStart, iErr)
      call checkError(node, iErr, "Bad lattice vectors, please check lines 3-5 of the geometry")
      geo%latVecs(:, ii) = latVec(:) * rScale
      iStart = iEnd + 1
    end do

    ! Here are the number of each species listed, or since Vasp >=5.1 element symbols
    iEnd = nextLine(text, iStart)
    iOldStart = iStart
    call getNextToken(text(:iEnd), iTmp, iOldStart, iErr)
    ! Seems to be element symbols, so we prefer them over the ones from the comment line
    if (iErr /= TOKEN_OK) then
      if (hasComment) then
        call destruct(speciesNames)
      end if
      call init(speciesNames)
      iErr = TOKEN_OK
      do while(iErr == TOKEN_OK)
        call getNextToken(text(:iEnd), txt, iStart, iErr)
        if (iErr == TOKEN_OK) then
          call append(speciesNames, char(txt))
        end if
      end do
      iStart = iEnd + 1
      iEnd = nextLine(text, iStart)
    end if

    ! Now we have to deal with repeating the same species in unsorted POSCARs
    allocate(vaspNames(len(speciesNames)))
    allocate(vaspSp(len(speciesNames)))
    ! convert to array of strings, otherwise we have no reliable way to access the elements
    call asArray(speciesNames, vaspNames)
    ! reset speciesNames to be populated with the unique species present
    call destruct(speciesNames)
    call init(speciesNames)
    iSp = 0
    do ii = 1, size(vaspNames, dim=1)
      iSp = find(speciesNames, trim(vaspNames(ii)))
      if (iSp == 0) then
        call append(speciesNames, trim(vaspNames(ii)))
        iSp = len(speciesNames)
      end if
      vaspSp(ii) = iSp
    end do

    geo%nSpecies = len(speciesNames)
    if (geo%nSpecies == 0) then
      call detailedError(node, "Number of species in the VASP geometry equals zero.")
    end if
    allocate(geo%speciesNames(geo%nSpecies))
    call asArray(speciesNames, geo%speciesNames)
    call destruct(speciesNames)
    allocate(countSp(size(vaspSp, dim=1)))
    do iSp = 1, size(vaspSp, dim=1)
      call getNextToken(text(:iEnd), countSp(iSp), iStart, iErr)
      call checkError(node, iErr, "Could not read number of species in geometry")
    end do
    geo%nAtom = sum(countSp)
    allocate(geo%species(geo%nAtom))
    allocate(geo%coords(3, geo%nAtom))
    ii = 0
    do iSp = 1, size(vaspSp, dim=1)
      geo%species(ii+1:ii+countSp(iSp)) = vaspSp(iSp)
      ii = ii + countSp(iSp)
    end do
    iStart = iEnd + 1

    ! Search for the selective dynamics keyword here
    iEnd = nextLine(text, iStart)
    call getNextToken(text(:iEnd), txt, iStart, iErr)
    if (scan(char(txt), "sS") == 1) then
      iStart = iEnd + 1
      iEnd = nextLine(text, iStart)
      call getNextToken(text(:iEnd), txt, iStart, iErr)
    end if
    if (scan(char(txt), "cCkK") == 1) then
      geo%tFracCoord = .false.
    else if (scan(char(txt), "dD") == 1) then
      geo%tFracCoord = .true.
    else
      call detailedError(node, "Unknown coordinate format '"// char(txt) // "'")
    end if
    iStart = iEnd + 1

    do ii = 1, geo%nAtom
      iEnd = nextLine(text, iStart)
      call getNextToken(text, coords, iStart, iErr)
      write(errorStr,"(A,1X,I0)")"Bad coordinates for atom", ii
      call checkError(node, iErr, trim(errorStr))
      geo%coords(:, ii) = coords(:)
      iStart = iEnd + 1
    end do

    call setupPeriodicGeometry(node, geo)

    ! convert coords to correct internal units
    if (geo%tFracCoord) then
      geo%coords = matmul(geo%latVecs, geo%coords)
    else
      geo%coords = geo%coords * AA__Bohr * rScale
    end if

    geo%tPeriodic = .true.
    geo%tHelical = .false.

    call normalize(geo)

  end subroutine readTGeometryVasp_help


  !> Reads the geometry in a HSD tree in LAMMPS data file format
  subroutine readTGeometryLammps(node, geo)

    !> Node containing the geometry in Gen format
    type(fnode), pointer :: node

    !> Contains the geometry on exit
    type(TGeometry), intent(out) :: geo

    type(fnode), pointer :: child
    type(string) :: text1, text2

    call getChildValue(node, "CommandFile", child)
    if (.not. associated(child) .or. getNodeType(child) /= TEXT_NODE) then
      call detailedError(node, "Missing CommandFile for LammpsFormat")
    end if
    call getNodeValue(child, text1)

    call getChildValue(node, "DataFile", child)
    if (.not. associated(child) .or. getNodeType(child) /= TEXT_NODE) then
      call detailedError(node, "Missing DataFile for LammpsFormat")
    end if
    call getNodeValue(child, text2)

    call readTGeometryLammps_help(node, geo, char(text1), char(text2))

  end subroutine readTGeometryLammps


  !> Helping routine for reading geometry from a HSD tree in LAMMPS format
  subroutine readTGeometryLammps_help(node, geo, commandInput, dataInput)

    !> Node to parse (only needed to produce proper error messages)
    type(fnode), pointer :: node

    !> Contains the geometry on exit
    type(TGeometry), intent(out) :: geo

    !> Text contents of the command file
    character(len=*), intent(in) :: commandInput

    !> Text contents of the data file
    character(len=*), intent(in) :: dataInput

    integer :: iStart, iEnd, iErr, intValue, skipItemsForCoords, i, j
    real(dp) :: realValue(3), toAngstrom, toAtomicMassUnit
    type(string) :: text, command
    logical :: haveMasses, haveAtoms, skipMoleculeId, readIntValue, readRealValue(3)

    haveMasses = .false.
    haveAtoms = .false.
    geo%nAtom = 0
    geo%nSpecies = 0
    allocate(geo%latVecs(3, 3))
    allocate(geo%origin(3))
    geo%origin(:) = -0.5_dp
    geo%latVecs(:,:) = 0.0_dp
    do i = 1, 3
      geo%latVecs(i,i) = 1.0_dp
    end do

    toAngstrom = 0.0_dp
    toAtomicMassUnit = 0.0_dp
    skipMoleculeId = .false.
    skipItemsForCoords = 0

    ! Read command file first

    iStart = 1
    iEnd = len(commandInput)
    do
      call getNextToken(commandInput, command, iStart, iErr)
      if (iErr /= TOKEN_OK) then
        exit
      end if

      select case(char(command))
      case("pair_style")
        call getNextToken(commandInput, text, iStart, iErr)
        call checkError(node, iErr, "Error reading pair_style in command file")
        if (char(text) /= "dftbplus") then
          call detailedError(node, "pair_style must be dftbplus")
        end if
        call jumpToEndOfLine(commandInput, iStart)
      case("atom_style") ! default is atomic
        call getNextToken(commandInput, text, iStart, iErr)
        call checkError(node, iErr, "Error reading atom_style in command file")
        select case(char(text))
        case("angle", "bond", "molecular")
          skipMoleculeId = .true.
        case("full")
          skipMoleculeId = .true.
          skipItemsForCoords = 1
        case("line", "tri")
          skipMoleculeId = .true.
          skipItemsForCoords = 2
        case("charge", "dielectric", "dipole", "dpd", "mdpd")
          skipItemsForCoords = 1
        case("edpd", "ellipsoid")
          skipItemsForCoords = 2
        case("electron", "sph", "template")
          skipItemsForCoords = 3
        case("wavepacket")
          skipItemsForCoords = 6
        case("body", "mesont", "peri", "smd", "sphere", "bpm/sphere")
          call detailedError(node, "Unsupported atom_style " // char(text))
        end select
        call jumpToEndOfLine(commandInput, iStart)
      case("boundary") ! default is p p p
        do i = 1, 3
          call getNextToken(commandInput, text, iStart, iErr)
          call checkError(node, iErr, "Error reading boundary in command file")
          if (char(text) /= "p") then
            call detailedError(node, "Only periodic boundary conditions supported")
          end if
        end do
        call jumpToEndOfLine(commandInput, iStart)
      case("units") ! default is lj
        call getNextToken(commandInput, text, iStart, iErr)
        call checkError(node, iErr, "Error reading units in command file")
        select case(char(text))
        case("si")
          toAngstrom = 1.0e10_dp ! from meters
          toAtomicMassUnit = 1.0e3_dp * avogadConst ! from kilograms
        case("cgs")
          toAngstrom = 1.0e8_dp ! from centimeters
          toAtomicMassUnit = avogadConst ! from grams
        case("electron")
          toAngstrom = Bohr__AA ! from Bohr radii
          toAtomicMassUnit = 1.0_dp
        case("micro")
          toAngstrom = 1.0e4_dp ! from micrometers
          toAtomicMassUnit = 1.0e-12_dp * avogadConst ! from picograms
        case("nano")
          toAngstrom = 10.0_dp ! from nanometers
          toAtomicMassUnit = 1.0e-18_dp * avogadConst ! from attograms
        case("lj")
          call detailedError(node, "Unit system lj is not supported")
        case default
          toAngstrom = 1.0_dp
          toAtomicMassUnit = 1.0_dp
        end select
        call jumpToEndOfLine(commandInput, iStart)
      end select
    end do

    ! Now read the data file defining the geometry and the atoms

    iStart = 1
    iEnd = len(dataInput)
    do
      ! Try to read real values at the beginning of the line
      readRealValue(:) = .false.
      do i = 1, size(realValue)
        call getNextToken(dataInput, realValue(i), iStart, iErr)
        readRealValue(i) = (iErr == TOKEN_OK)
        if (.not. readRealValue(i)) then
          exit
        end if
      end do

      ! Try to convert the first real value to integer
      if (readRealValue(1)) then
        intValue = int(realValue(1))
        if (real(intValue, kind=dp) == realValue(1)) then
          readIntValue = .true.
        end if
      end if

      ! Read the actual command now
      call getNextToken(dataInput, command, iStart, iErr)
      if (iErr /= TOKEN_OK) then
        exit
      end if

      select case(char(command))
      case("atoms")
        if (.not. readIntValue) then
          call detailedError(node, "Invalid number of atoms")
        end if
        geo%nAtom = intValue
        call jumpToEndOfLine(dataInput, iStart)
      case("atom")
        call getNextToken(dataInput, text, iStart, iErr)
        if (iErr == TOKEN_OK .and. char(text) == "types") then
          if (.not. readIntValue) then
            call detailedError(node, "Invalid number of species")
          end if
          geo%nSpecies = intValue
          call jumpToEndOfLine(dataInput, iStart)
        end if
      case("xlo")
        call getNextToken(dataInput, text, iStart, iErr)
        if (iErr == TOKEN_OK .and. char(text) == "xhi") then
          if (.not. readRealValue(1) .or. .not. readRealValue(2)) then
            call detailedError(node, "Invalid values for xlo and/or xhi")
          end if
          geo%origin(1) = realValue(1)
          geo%latVecs(1,1) = realValue(2) - realValue(1)
          call jumpToEndOfLine(dataInput, iStart)
        end if
      case("ylo")
        call getNextToken(dataInput, text, iStart, iErr)
        if (iErr == TOKEN_OK .and. char(text) == "yhi") then
          if (.not. readRealValue(1) .or. .not. readRealValue(2)) then
            call detailedError(node, "Invalid values for ylo and/or yhi")
          end if
          geo%origin(2) = realValue(1)
          geo%latVecs(2,2) = realValue(2) - realValue(1)
          call jumpToEndOfLine(dataInput, iStart)
        end if
      case("zlo")
        call getNextToken(dataInput, text, iStart, iErr)
        if (iErr == TOKEN_OK .and. char(text) == "zhi") then
          if (.not. readRealValue(1) .or. .not. readRealValue(2)) then
            call detailedError(node, "Invalid values for zlo and/or zhi")
          end if
          geo%origin(3) = realValue(1)
          geo%latVecs(3,3) = realValue(2) - realValue(1)
          call jumpToEndOfLine(dataInput, iStart)
        end if
      case("xy")
        call getNextToken(dataInput, text, iStart, iErr)
        if (iErr == TOKEN_OK .and. char(text) == "xz") then
          call getNextToken(dataInput, text, iStart, iErr)
          if (iErr == TOKEN_OK .and. char(text) == "yz") then
            if (any(.not. readRealValue)) then
              call detailedError(node, "Invalid values for xy/xz/yz")
            end if
            geo%latVecs(1,2) = realValue(1) ! xy
            geo%latVecs(1,3) = realValue(2) ! xz
            geo%latVecs(2,3) = realValue(3) ! yz
          end if
        end if
      case("Masses")
        if (geo%nSpecies == 0) then
          call detailedError(node, "Missing number of species (atom types) in data file header")
        end if
        allocate(geo%speciesNames(geo%nSpecies))
        call jumpToEndOfLine(dataInput, iStart)

        do i = 1, geo%nSpecies
          call getNextToken(dataInput, intValue, iStart, iErr)
          call checkError(node, iErr, "Invalid value for species index")
          if (intValue /= i) then
            call detailedError(node, "Unexpected species index")
          end if
          call getNextToken(dataInput, realValue(1), iStart, iErr)
          call checkError(node, iErr, "Invalid value for species mass")
          geo%speciesNames(i) = trim(getAtomicSymbol(realValue(1) * toAtomicMassUnit))
          call jumpToEndOfLine(dataInput, iStart)
        end do
        haveMasses = .true.
      case("Atoms")
        if (geo%nAtom == 0) then
          call detailedError(node, "Missing number of atoms in data file header")
        end if
        allocate(geo%species(geo%nAtom))
        allocate(geo%coords(3, geo%nAtom))
        call jumpToEndOfLine(dataInput, iStart)

        do i = 1, geo%nAtom
          call getNextToken(dataInput, intValue, iStart, iErr)
          call checkError(node, iErr, "Invalid value for atom index")
          if (intValue /= i) then
            call detailedError(node, "Unexpected atom index")
          end if
          if (skipMoleculeId) then
            call getNextToken(dataInput, text, iStart, iErr)
          end if
          call getNextToken(dataInput, intValue, iStart, iErr)
          call checkError(node, iErr, "Invalid value for atom species")
          geo%species(i) = intValue
          if (geo%species(i) > geo%nSpecies) then
            call detailedError(node, "Unexpected species index")
          end if
          iEnd = nextLine(dataInput, iStart)
          do j = 1, skipItemsForCoords
            call getNextToken(dataInput(:iEnd), text, iStart, iErr)
          end do
          call getNextToken(dataInput(:iEnd), realValue(1:3), iStart, iErr)
          if (iErr == TOKEN_OK) then
            geo%coords(:,i) = realValue(:)
          else
            geo%coords(:,i) = 0.0_dp
          end if
          call jumpToEndOfLine(dataInput, iStart)
        end do
        haveAtoms = .true.
      end select
    end do

    if (.not. haveMasses) then
      call detailedError(node, "Missing Masses section in data file")
    end if
    if (.not. haveAtoms) then
      call detailedError(node, "Missing Atoms section in data file")
    end if
    if (toAngstrom == 0.0_dp) then
      call detailedError(node, "Require explicit units definition in command file")
    end if

    geo%origin(:) = toAngstrom * geo%origin(:)
    geo%latVecs(:,:) = toAngstrom * geo%latVecs(:,:)
    geo%coords(:,:) = toAngstrom * geo%coords(:,:)

    geo%tPeriodic = .true.
    geo%tFracCoord = .false.
    geo%tHelical = .false.

    call setupPeriodicGeometry(node, geo)
    geo%coords = geo%coords * AA__Bohr

    write(stdout, "(A,I8,A,I3,A)") "Read values from LAMMPS input file: ",&
        & geo%nAtom, " atoms, ", geo%nSpecies, " species"

    call normalize(geo)

  end subroutine readTGeometryLammps_help


  !> Common checks for periodic input and generation of associated information
  subroutine setupPeriodicGeometry(node, geo)

    !> Node to parse (only needed to produce proper error messages)
    type(fnode), pointer :: node

    !> Contains the geometry on exit
    type(TGeometry), intent(inout) :: geo

    real(dp) :: det

    geo%origin = geo%origin * AA__Bohr
    geo%latVecs = geo%latVecs * AA__Bohr
    if (geo%tFracCoord) then
      if (any(abs(geo%coords) > 1.0_dp)) then
        call detailedWarning(node, "Fractional coordinates with absolute value greater than one.")
      end if
    end if
    allocate(geo%recVecs2p(3, 3))
    det = determinant33(geo%latVecs)
    if (abs(det) < 1e-12_dp) then
      call detailedError(node, "Dependent lattice vectors")
    end if
    call invert33(geo%recVecs2p, geo%latVecs, det)

  end subroutine setupPeriodicGeometry


  !> Return index for next line ending within text after position iStart
  pure function nextLine(text, iStart) result(iEnd)

    !> Text content of the node
    character(len=*), intent(in) :: text

    !> Start of the text content consider
    integer, intent(in) :: iStart

    !> End of the line, as delimited by a new-line character
    integer :: iEnd

    iEnd = index(text(iStart:), new_line(text)) + iStart - 1
    if (iEnd < iStart) iEnd = len(text)

  end function nextLine


  !> Advance iStart to the end of the line, but only if currently not at the beginning of the line
  subroutine jumpToEndOfLine(text, iStart)

    !> Text content
    character(len=*), intent(in) :: text

    !> Start of the text of interest; end of the line on exit
    integer, intent(inout) :: iStart

    if (iStart > 1 .and. iStart < len(text)) then
      if (text(iStart-1:iStart-1) /= new_line(text)) then
        iStart = nextLine(text, iStart)
      end if
    end if

  end subroutine jumpToEndOfLine


end module dftbp_type_typegeometryhsd
