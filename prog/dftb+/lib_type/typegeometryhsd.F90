!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Routines to read/write a TGeometry type in HSD and XML format.
module typegeometryhsd
  use typegeometry
  use hsdutils
  use hsdutils2
  use tokenreader
  use unitconversion
  use linkedlist
  use charmanip
  use simplealgebra, only : invert33, determinant33
  use xmlf90, flib_normalize => normalize
  implicit none
  private


  !> Writes the content of a geometry object to a dom tree or to an xml-writer
  interface writeTGeometryHSD
    module procedure writeTGeometryHSD_dom
    module procedure writeTGeometryHSD_xmlf
  end interface


  !> Types/subroutines from TypeGeometry
  public :: TGeometry, normalize


  !> Locally defined subroutines
  public :: writeTGeometryHSD, readTGeometryHSD, readTGeometryGen

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
    if (geo%tPeriodic) then
      call setChildValue(node, "LatticeVectors", geo%latVecs, .false.)
    end if

  end subroutine writeTGeometryHSD_dom


  !> Write the geometry in HSD format to an xml writer
  subroutine writeTGeometryHSD_xmlf(xf, geo)

    !> Node in the HSD-tree which should contain the geometry
    type(xmlf_t), intent(inout) :: xf

    !> The geometry
    type(TGeometry), intent(in) :: geo

    call writeChildValue(xf, "TypeNames", geo%speciesNames)
    call writeChildValue(xf, "TypesAndCoordinates", &
        &reshape(geo%species, (/ 1, size(geo%species) /)), geo%coords)
    call writeChildValue(xf, "Periodic", geo%tPeriodic)
    if (geo%tPeriodic) then
      call writeChildValue(xf, "LatticeVectors", geo%latVecs)
    end if

  end subroutine writeTGeometryHSD_xmlf


  !> Read the geometry from a node in a HSD tree.
  subroutine readTGeometryHSD(node, geo)

    !> Node in the HSD tree containing the geomery
    type(fnode), pointer :: node

    !> Contains the geometry on exit
    type(TGeometry), intent(out) :: geo

    type(string) :: modifier
    integer :: ind
    type(listString) :: stringBuffer
    type(listRealR1) :: realBuffer
    type(listIntR1) :: intBuffer
    type(fnode), pointer :: child, typesAndCoords
    integer, allocatable :: tmpInt(:,:)
    real(dp) :: latvec(9), det

    call getChildValue(node, "Periodic", geo%tPeriodic, default=.false.)
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
        ind = getModifierIndex(char(modifier), lengthUnits, typesAndCoords)
        geo%coords(:,:) = geo%coords(:,:) * lengthUnits(ind)%convertValue
        call setChildValue(typesAndCoords, "", &
            &reshape(geo%species, (/ 1, size(geo%species) /)), geo%coords, &
            &replace=.true.)
      end select
    end if
    if (geo%tPeriodic) then
      allocate(geo%latVecs(3,3))
      call getChildValue(node, "LatticeVectors", latvec, modifier=modifier, &
          &child=child)
      geo%latVecs(:,:) = reshape(latvec, (/3, 3/))
      if (len(modifier) > 0) then
        ind = getModifierIndex(char(modifier), lengthUnits, child)
        geo%latVecs = geo%latVecs(:,:) * lengthUnits(ind)%convertValue
        call setChildValue(child, "", geo%latVecs, .true.)
      end if
      if (geo%tFracCoord) then
        geo%coords = matmul(geo%latVecs, geo%coords)
      end if
      allocate(geo%recVecs2p(3, 3))
      det = determinant33(geo%latVecs)
      if (abs(det) < 1e-12_dp) then
        call detailedError(child, "Dependent lattice vectors")
      end if
      call invert33(geo%recVecs2p, geo%latVecs, det)
      geo%recVecs2p(:,:) = reshape(geo%recVecs2p, (/3, 3/), order=(/2, 1/))
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
    integer :: iStart, iOldStart, iErr
    integer :: ii, iTmp
    real(dp) :: coords(3), rTmp, det
    type(listString) :: speciesNames

    ! Read first line of the gen file: Number of atoms, boundary conditions
    iStart = 1
    call getNextToken(text, geo%nAtom, iStart, iErr)
    call checkError(node, iErr, "Bad number of atoms.")
    call getNextToken(text, txt, iStart, iErr)
    select case (char(txt))
    case("S","s")
      geo%tPeriodic = .true.
      geo%tFracCoord = .false.
    case("F","f")
      geo%tPeriodic = .true.
      geo%tFracCoord = .true.
    case("C", "c")
      geo%tPeriodic = .false.
      geo%tFracCoord = .false.
    case default
      call detailedError(node, "Unknown boundary condition type '" &
          &// char(txt) // "'")
    end select

    ! Reading the 2nd line of a gen file.
    ! Since we cannot rely on line breaks, we try to read in integers. If that fails, we had a
    ! species name instead so it will be read as string.
    call init(speciesNames)
    iErr = TOKEN_ERROR
    iOldStart = iStart
    call getNextToken(text, iTmp, iStart, iErr)
    do while (iErr /= TOKEN_OK)
      call getNextToken(text, txt, iStart, iErr)
      if (iErr /= TOKEN_OK) then
        call detailedError(node, "Unexpected end of geometry information.")
      end if
      call append(speciesNames, char(txt))
      iOldStart = iStart
      call getNextToken(text, iTmp, iStart, iErr)
    end do
    geo%nSpecies = len(speciesNames)
    if (geo%nSpecies == 0) then
      call detailedError(node, "Number of species equals zero.")
    end if
    allocate(geo%speciesNames(geo%nSpecies))
    call asArray(speciesNames, geo%speciesNames)
    call destruct(speciesNames)

    ! Read in sequential and species indices.
    allocate(geo%species(geo%nAtom))
    allocate(geo%coords(3, geo%nAtom))
    iStart = iOldStart
    do ii = 1, geo%nAtom
      call getNextToken(text, iTmp, iStart, iErr)
      call checkError(node, iErr, "Bad sequentual number for an atom.")
      call getNextToken(text, geo%species(ii), iStart, iErr)
      call checkError(node, iErr, "Bad species number for an atom.")
      call getNextToken(text, coords, iStart, iErr)
      call checkError(node, iErr, "Bad coordinates for an atom.")
      geo%coords(:, ii) = coords(:)
    end do
    if (geo%nSpecies /= maxval(geo%species) .or. minval(geo%species) /= 1) then
      call detailedError(node, &
          &"Nr. of species and nr. of specified elements do not match.")
    end if

    ! Read in origin an lattice vectors, if the structure is periodic
    if (geo%tPeriodic) then
      allocate(geo%origin(3))
      allocate(geo%latVecs(3, 3))
      call getNextToken(text, geo%origin, iStart, iErr)
      call checkError(node, iErr, "Invalid origin.")
      do ii = 1, 3
        call getNextToken(text, geo%latVecs(:, ii), iStart, iErr)
        call checkError(node, iErr, "Invalid lattice vectors.")
      end do
    end if

    ! Check if any data remains in the geometry - should be nothing left now
    call getNextToken(text, rTmp, iStart, iErr)
    if (iErr /= TOKEN_EOS) then
      call detailedError(node, "Superfluous data found. Check if specified number of atoms matches&
          & the number of actually entered positions.")
    end if

    ! tests that are relevant to periodic geometries only
    if (geo%tPeriodic) then
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
    end if

    ! convert coords to correct internal units
    if (geo%tFracCoord) then
      geo%coords = matmul(geo%latVecs, geo%coords)
    else
      geo%coords = geo%coords * AA__Bohr
    end if

    call normalize(geo)

  end subroutine readTGeometryGen_help

end module typegeometryhsd
