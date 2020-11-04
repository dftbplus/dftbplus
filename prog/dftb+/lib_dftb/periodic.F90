!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines for the periodic boundary conditions and neighbour data
module dftbp_periodic
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_constants, only : pi
  use dftbp_message
  use dftbp_sorting
  use dftbp_bisect
  use dftbp_linkedlist
  use dftbp_simplealgebra, only : determinant33, invert33
  use dftbp_commontypes
  use dftbp_memman
  use dftbp_latpointiter
  use dftbp_quaternions, only : rotate3
  use dftbp_boundarycond, only : zAxis
  implicit none

  private

  public :: getCellTranslations, getLatticePoints, foldCoordToUnitCell
  public :: reallocateHS, buildSquaredAtomIndex
  public :: TNeighbourList, TNeighbourlist_init
  public :: updateNeighbourList, updateNeighbourListAndSpecies
  public :: getNrOfNeighbours, getNrOfNeighboursForAll
  public :: getSuperSampling
  public :: frac2cart, cart2frac, cyl2cart, cart2cyl
  public :: getSparseDescriptor


  !> resize sparse arrays
  interface reallocateHS
    module procedure reallocateHS_1
    module procedure reallocateHS_2
    module procedure reallocateHS_Single
  end interface reallocateHS


  !> convert fractional coordinates to cartesian
  interface frac2cart
    module procedure fractionalCartesian
  end interface frac2cart


  !> cartesian to fractional coordinates in periodic geometry
  interface cart2frac
    module procedure cartesianFractional
  end interface cart2frac

  !> routine to convert from cylindrical to cartesian coordinate systems
  interface cyl2cart
    module procedure cyl2cart_vec
    module procedure cyl2cart_array
  end interface cyl2cart

  !> routine to convert from cartesian to cylindrical coordinate systems
  interface cart2cyl
    module procedure cart2cyl_vec
  end interface cart2cyl


  !> Contains essential data for the neighbourlist
  type TNeighbourList

    !> index of neighbour atoms
    integer, allocatable :: iNeighbour(:,:)

    !> nr. of neighbours
    integer, allocatable :: nNeighbour(:)

    !> temporary array for neighbour distances
    real(dp), allocatable :: neighDist2(:,:)

    !> cutoff it was generated for
    real(dp) :: cutoff

    !> initialised data
    logical :: initialized = .false.
  end type TNeighbourList

contains


  !> Initializes a neighbourlist instance.
  subroutine TNeighbourlist_init(neighbourList, nAtom, nInitNeighbour)

    !> Neighbourlist data.
    type(TNeighbourList), intent(out) :: neighbourList

    !> Nr. of atoms in the system.
    integer, intent(in) :: nAtom

    !> Expected nr. of neighbours per atom.
    integer, intent(in) :: nInitNeighbour

    @:ASSERT(.not. neighbourList%initialized)
    @:ASSERT(nAtom > 0)
    @:ASSERT(nInitNeighbour > 0)

    allocate(neighbourList%nNeighbour(nAtom))
    allocate(neighbourList%iNeighbour(0:nInitNeighbour, nAtom))
    allocate(neighbourList%neighDist2(0:nInitNeighbour, nAtom))

    neighbourList%cutoff = -1.0_dp
    neighbourList%initialized = .true.

  end subroutine TNeighbourlist_init


  !> Calculates the translation vectors for cells, which could contain atoms interacting with any of
  !> the atoms in the central cell.
  !> This subroutine uses a simple guess to get the necessary translation vectors. This results in a
  !> set of vectors wich could for very asymmetric cells a large amount bigger than the real
  !> necessary one.
  subroutine getCellTranslations(cellVec, rCellVec, latVec, recVec2p, cutoff)

    !> Returns cell translation vectors in relative coordinates.
    real(dp), allocatable, intent(out) :: cellVec(:, :)

    !> Returns cell translation vectors in absolute units.
    real(dp), allocatable, intent(out) :: rCellVec(:,:)

    !> Lattice vectors
    real(dp), intent(in) :: latVec(:,:)

    !> Reciprocal lattice vectors in 2*pi units.
    real(dp), intent(in) :: recVec2p(:,:)

    !> Global cutoff for the diatomic interactions
    real(dp), intent(in) :: cutoff

    integer :: ii

    if (all(shape(latVec) == (/3, 3/))) then
      call getLatticePoints(cellVec, latVec, recVec2p, cutoff, posExtension=1, negExtension=1)
      allocate(rCellVec(3, size(cellVec, dim=2)))
      do ii = 1, size(rCellVec, dim=2)
        rCellVec(:,ii) = matmul(latVec, cellVec(:,ii))
      end do
    else if (all(shape(latVec) == (/3, 1/))) then
      ! Helical
      call getHelicalPoints(cellVec, rCellVec, latVec, cutoff)
    else
      call error("Miss-shaped cell vectors in getCellTranslations.")
    end if

  end subroutine getCellTranslations


  !> Returns a set which definitely contains all the points of a 3D grid which are nearer to the
  !> origin than a given distance.
  !> Without the onlyInside parameter, the returned set of lattice points shape a
  !> parallelepipedon. With the onlyInside parameter its not necessarily the case.
  !> Refine the algorithm with the help of a new routine which can calculate the minimal distance
  !> between two arbitary cells.
  subroutine getLatticePoints(latPoint, latVec, recVec2p, dist, posExtension, negExtension,&
      & onlyInside, reduceByInversion, withoutOrigin)

    !> Returns grid points in relative coords.
    real(dp), allocatable, intent(out) :: latPoint(:,:)

    !> Lattice vectors.
    real(dp), intent(in) :: latVec(:,:)

    !> Reciprocal lattice vectors in 2*pi units.
    real(dp), intent(in) :: recVec2p(:,:)

    !> Global cutoff for the diatomic interactions.
    real(dp), intent(in) :: dist

    !> Extend the set along the positive lattice vectors with that many additional lattice vectors.
    integer, intent(in), optional :: posExtension

    !> Same as posExtension for negative lattice vectors
    integer, intent(in), optional :: negExtension

    !> Return only those lattice points which are really not outside the given distance.
    logical,  intent(in), optional :: onlyInside

    !> whether to include time reversal symmetry when generating k-points
    logical,  intent(in), optional :: reduceByInversion

    !> whether to exclude the (0,0,0) point
    logical,  intent(in), optional :: withoutOrigin

    type(TLatPointIter) :: latPointGen
    real(dp), allocatable :: tmpLatPoint(:,:)
    logical :: tWithOrigin

    @:ASSERT(dist >= 0.0_dp)

    ! In order to make sure that (0,0,0) is the first point, the points that are generated by the
    ! iterator will not contain the origin, this is instead inserted as the first element later
    call TLatPointIter_init(latPointGen, latVec, recVec2p, dist, negExtension, posExtension,&
        & onlyInside, reduceByInversion, excludeOrigin=.true.)
    call latPointGen%getAllPoints(tmpLatPoint)
    if (present(withoutOrigin)) then
      tWithOrigin = .not. withoutOrigin
    else
      tWithOrigin = .true.
    end if
    if (tWithOrigin) then
      ! insert origin first
      allocate(latPoint(3, size(tmpLatPoint, dim=2) + 1))
      latPoint(:, 1) = 0.0_dp
      latPoint(:, 2:) = tmpLatPoint
    else
      call move_alloc(tmpLatPoint, latPoint)
    end if

  end subroutine getLatticePoints


  !> Cells in a helical arrangement
  subroutine getHelicalPoints(cellVec, rCellVec, latVec, cutoff)

    !> Returns cell translation vectors in relative coordinates.
    real(dp), allocatable, intent(out) :: cellVec(:, :)

    !> Returns cell translation vectors in absolute units.
    real(dp), allocatable, intent(out) :: rCellVec(:,:)

    !> Lattice vectors.
    real(dp), intent(in) :: latVec(:,:)

    !> Global cutoff for the diatomic interactions
    real(dp), intent(in) :: cutoff

    integer :: maxCells, ii

    ! cell extension along helix in +ve sense
    maxCells = ceiling(cutoff / latVec(1,1)) +1
    ! Total helix cells
    maxCells = 2 * maxCells
    allocate(cellVec(2,nint(latVec(3,1)) * (maxCells+1)))
    allocate(rCellVec(2,nint(latVec(3,1)) * (maxCells+1)))
    cellVec(:,:) = 0.0_dp
    ! Helix operation
    do ii = 2, maxCells, 2
      cellVec(1,ii) = real(ii / 2,dp)
      cellVec(1,ii+1) = real(-ii / 2,dp)
    end do
    ! c_n operation, duplicating helical points
    do ii = 2, nint(latVec(3,1))
      cellVec(1,(ii-1)*(maxCells+1)+1:(ii)*(maxCells+1)) = cellVec(1,:maxCells+1)
      cellVec(2,(ii-1)*(maxCells+1)+1:(ii)*(maxCells+1)) = ii - 1
    end do
    rCellVec(1,:) = latVec(1,1) * cellVec(1,:)
    rCellVec(2,:) = cellVec(2,:)

  end subroutine getHelicalPoints

  !> Fold coordinates back in the central cell.
  !>
  !> Throw away the integer part of the relative coordinates of every atom. If the resulting
  !> coordinate is very near to 1.0 (closer than 1e-12 in absolute length), fold it to 0.0 to make
  !> the algorithm more predictable and independent of numerical noise.
  subroutine foldCoordToUnitCell(coord, latVec, recVec2p)

    !> Contains the original coordinates on call and the folded ones on return.
    real(dp), intent(inout) :: coord(:,:)

    !> Lattice vectors (column format).
    real(dp), intent(in) :: latVec(:,:)

    !> Reciprocal vectors in units of 2pi (column format).
    real(dp), intent(in) :: recVec2p(:,:)


    !> Nr. of atoms in the cell.
    integer :: nAtom

    integer :: ii, jj
    real(dp) :: frac(3), frac2(3), tmp3(3), vecLen(3), thetaNew, thetaOld

    nAtom = size(coord, dim=2)

    @:ASSERT(size(coord, dim=1) == 3)
    @:ASSERT(all(shape(latVec) == (/3, 3/)) .or. all(shape(latVec) == (/3, 1/)))

    if (all(shape(latVec) == (/3, 3/))) then

      vecLen(:) = sqrt(sum(latVec(:,:)**2, dim=1))
      do ii = 1, nAtom
        do jj = 1, 3
          frac(jj) = dot_product(recVec2p(:,jj), coord(:,ii))
        end do
        tmp3(:) = coord(:,ii)
        frac2(:) = frac(:) - real(floor(frac(:)), dp)
        where (abs(vecLen*(1.0_dp - frac2)) < 1e-12_dp) frac2 = 0.0_dp
        coord(:, ii) = matmul(latVec, frac2)
      end do

    else if (all(shape(latVec) == (/3, 1/))) then

      do ii = 1, nAtom

        jj = floor(coord(3,ii)/latVec(1,1))
        ! want coordinate in eventual range 0..latVec(1,1) hence floor

        tmp3(:) = coord(:,ii)
        coord(3,ii) = coord(3,ii) - jj * latVec(1,1)
        call rotate3(coord(:,ii),-jj*latVec(2,1),zAxis)
        thetaOld = atan2(coord(2,ii),coord(1,ii))
        thetaNew = mod(thetaOld+2.0_dp*pi,2.0_dp*pi/latvec(3,1))
        call rotate3(coord(:,ii),-thetaOld+thetaNew,zAxis)

      end do

    else

      call error("Miss-shaped cell vectors in foldCoordToUnitCell.")

    end if

  end subroutine foldCoordToUnitCell


  !> Updates the neighbour list and the species arrays.
  subroutine updateNeighbourListAndSpecies(coord, species, img2CentCell, iCellVec, neigh, nAllAtom,&
      & coord0, species0, cutoff, rCellVec, symmetric, helicalBoundConds)

    !> Coordinates of all interacting atoms on exit
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Species of all interacting atoms on exit.
    integer, allocatable, intent(inout) :: species(:)

    !> Mapping on atoms in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Shift vector index for every interacting atom
    integer, allocatable, intent(inout) :: iCellVec(:)

    !> Updated neighbour list.
    type(TNeighbourList), intent(inout) :: neigh

    !> Number of all interacting atoms
    integer, intent(out) :: nAllAtom

    !> Coordinates of the atoms in the central cell
    real(dp), intent(in) :: coord0(:,:)

    !> Species of the atoms in the central cell
    integer, intent(in) :: species0(:)

    !> Cutoff until neighbourlist should be created
    real(dp), intent(in) :: cutoff

    !> Cell vector for the translated cells to consider.
    real(dp), intent(in) :: rCellVec(:,:)

    !> Whether the neighbour list should be symmetric or not (default)
    logical, intent(in), optional :: symmetric

    !> Helical translation and angle, if neccessary, along z axis
    real(dp), intent(in), optional :: helicalBoundConds(:,:)

    call updateNeighbourList(coord, img2CentCell, iCellVec, neigh, nAllAtom, coord0, cutoff,&
        & rCellVec, symmetric, helicalBoundConds)

    if (size(species) /= nAllAtom) then
      deallocate(species)
      allocate(species(nAllAtom))
    end if
    species(1:nAllAtom) = species0(img2CentCell(1:nAllAtom))

  end subroutine updateNeighbourListAndSpecies


  !> Updates the neighbour list according a given geometry.
  !> The neighbourlist for the given cutoff is calculated. Arrays are resized if necessary. The
  !> neighbour list determination is a simple N^2 algorithm, calculating the distance between the
  !> possible atom pairs.
  subroutine updateNeighbourList(coord, img2CentCell, iCellVec, neigh, nAllAtom, coord0, cutoff,&
      & rCellVec, symmetric, helicalBoundConds)

    !> Coordinates of the objects interacting with the objects in the central cell (on exit).
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Returns for all objects the index of an object in the central cell which the object is mapped
    !> on to.
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Returns the index of the translating superlattice vector for each object.
    integer, allocatable, intent(inout) :: iCellVec(:)

    !> Neighbourlist.
    type(TNeighbourList), intent(inout) :: neigh

    !> Returns the nr. of all objects (including those in the translated cells.)
    integer, intent(out) :: nAllAtom

    !> Coordinates of the objects in the central cell.
    real(dp), intent(in) :: coord0(:,:)

    !> Cutoff radius for the interactions.
    real(dp), intent(in) :: cutoff

    !> Absolute coordinates of the shifted supercells which could have interacting
    !> atoms with the central cell.
    real(dp), intent(in) :: rCellVec(:,:)

    !> Optional, whether the map should be symmetric (dftb default = .false.)
    logical, intent(in), optional :: symmetric

    !> Helical translation and angle, if neccessary, along z axis
    real(dp), intent(in), optional :: helicalBoundConds(:,:)

    ! Nr. of atoms in the system
    integer :: nAtom

    ! Max. nr. of atom without reallocation
    integer :: mAtom

    ! Max. nr. of neighbours without reallocation
    integer :: maxNeighbour

    ! Nr. of cell translation vectors
    integer :: nCellVec

    ! Square of the diatomic interaction cutoffs
    real(dp) :: cutoff2

    real(dp) :: dist2
    real(dp) :: rCell(3), rr(3)
    integer :: ii, iAtom1, oldIAtom1, iAtom2
    integer :: nn1, iAtom2End
    logical :: symm
    integer, allocatable :: indx(:)
    character(len=100) :: strError

    nAtom = size(neigh%nNeighbour, dim=1)
    mAtom = size(coord, dim=2)
    maxNeighbour = ubound(neigh%iNeighbour, dim=1)
    nCellVec = size(rCellVec, dim=2)

    @:ASSERT(nAtom <= mAtom)
    @:ASSERT(allocated(coord))
    @:ASSERT(size(coord, dim=1) == 3)
    @:ASSERT(allocated(img2CentCell))
    @:ASSERT(size(img2CentCell) == mAtom)
    @:ASSERT(allocated(iCellVec))
    @:ASSERT(size(iCellVec) == mAtom)
    @:ASSERT(size(neigh%iNeighbour, dim=2) == nAtom)
    @:ASSERT((size(coord0, dim=1) == 3) .and. size(coord0, dim=2) >= nAtom)
    @:ASSERT(cutoff >= 0.0_dp)

    symm = .false.
    if (present(symmetric)) then
      symm = symmetric
    end if
    neigh%cutoff = cutoff
    cutoff2 = cutoff**2
    nAllAtom = 0

    ! Clean arrays.
    !  (Every atom is the 0th neighbour of itself with zero distance square.)
    neigh%nNeighbour(:) = 0
    neigh%iNeighbour(:,:) = 0
    do ii = 1, nAtom
      neigh%iNeighbour(0, ii) = ii
    end do
    neigh%neighDist2(:,:) = 0.0_dp

    rCell(:) = 0.0_dp

    ! Loop over all possible neighbours for all atoms in the central cell.
    ! Only those neighbours are considered which map on atom with a higher
    ! or equal index in the central cell.
    ! Outer two loops: all atoms in all cells.
    ! Inner loop: all atoms in the central cell.
    lpCellVec: do ii = 1, nCellVec
      if (present(helicalBoundConds)) then
        ! helical structure
        rCell(:) = 0.0_dp
        rCell(3) = rCellVec(1, ii)
      else
        rCell(:) = rCellVec(:, ii)
      end if
      oldIAtom1 = 0
      lpIAtom1: do iAtom1 = 1, nAtom
        rr(:) = coord0(:, iAtom1) + rCell(:)
        if (symm) then
          iAtom2End = nAtom
        else
          iAtom2End = iAtom1
        end if
        if (present(helicalBoundConds)) then
          ! helical geometry
          if (size(helicalBoundConds,dim=1)==3) then
            ! an additional C rotation operation
            call rotate3(rr,2.0_dp*pi*rCellVec(2, ii)/helicalBoundConds(3,1), zAxis)
          end if
          ! helical operation, note nint() not floor() as roundoff can cause problems for floor
          ! here.
          call rotate3(rr,helicalBoundConds(2,1)*nint(rCellVec(1, ii)/helicalBoundConds(1,1)),&
              & zAxis)
        end if
        lpIAtom2: do iAtom2 = 1, iAtom2End
          !  If distance greater than cutoff -> skip
          dist2 = sum((coord0(:, iAtom2) - rr(:))**2)
          if (dist2 > cutoff2) then
            cycle
          end if
          ! New interacting atom -> append
          ! We need that before checking for interaction with dummy atom or
          ! with itself to make sure that atoms in the central cell are
          ! appended  exactly in the same order as found in the coord0 array.
          if (iAtom1 /= oldIAtom1) then
            nAllAtom = nAllAtom + 1
            if (nAllAtom > mAtom) then
              mAtom = incrmntOfArray(mAtom)
              call reallocateArrays1(img2CentCell, iCellVec, coord, mAtom)
            end if
            coord(:, nAllAtom) = rr(:)
            img2CentCell(nAllAtom) = iAtom1
            iCellVec(nAllAtom) = ii
            oldIAtom1 = iAtom1
          end if

          ! Check if atoms are not too close to each other
          if (dist2 < minNeighDist2) then
            if (ii == 1 .and. iAtom1 == iAtom2) then
              ! We calculated the distance between the same atom in the unit cell
              cycle
            else
99000         format ('Atoms ',I5,' and ',I5,' too close to each other!', ' (dist=',E13.6,')')
              write (strError, 99000) iAtom2, nAllAtom, sqrt(dist2)
              call warning(strError)
            end if
          end if

          neigh%nNeighbour(iAtom2) = neigh%nNeighbour(iAtom2) + 1
          if (neigh%nNeighbour(iAtom2) > maxNeighbour) then
            maxNeighbour = incrmntOfArray(maxNeighbour)
            call reallocateArrays3(neigh%iNeighbour, neigh%neighDist2, maxNeighbour)
          end if
          neigh%iNeighbour(neigh%nNeighbour(iAtom2), iAtom2) = nAllAtom
          neigh%neighDist2(neigh%nNeighbour(iAtom2), iAtom2) = dist2

        end do lpIAtom2
      end do lpIAtom1
    end do lpCellVec

    ! Sort neighbours for all atom by distance
    allocate(indx(maxNeighbour))
    do iAtom1 = 1, nAtom
      nn1 = neigh%nNeighbour(iAtom1)
      call index_heap_sort(indx(1:nn1), neigh%neighDist2(1:nn1, iAtom1), tolSameDist2)
      neigh%iNeighbour(1:nn1, iAtom1) = neigh%iNeighbour(indx(:nn1), iAtom1)
      neigh%neighDist2(1:nn1, iAtom1) = neigh%neighDist2(indx(:nn1), iAtom1)
    end do

    call reallocateArrays1(img2CentCell, iCellVec, coord, nAllAtom)

    ! check for atoms on top of each other
    do iAtom1 = 1, nAtom
      do nn1 = 1, neigh%nNeighbour(iAtom1)
        if (neigh%neighDist2(nn1, iAtom1) < minNeighDist) then
          iAtom2 = img2CentCell(neigh%iNeighbour(nn1, iAtom1))
          write (strError, "(A,I0,A,I0,A)") "Atoms ",iAtom1, " and ", iAtom2, " too close together"
          call error(strError)
        end if
      end do
    end do

  end subroutine updateNeighbourList


  !> Returns the nr. of neighbours for a given cutoff for all atoms.
  subroutine getNrOfNeighboursForAll(nNeighbourSK, neigh, cutoff)

    !> Contains the nr. of neighbours for each atom on exit.
    integer, intent(out) :: nNeighbourSK(:)

    !> Initialized neighbourlist
    type(TNeighbourList), intent(in) :: neigh

    !> Maximal neighbour distance to consider.
    real(dp),            intent(in) :: cutoff

    integer :: nAtom, iAtom

    nAtom = size(nNeighbourSK)

    @:ASSERT(size(neigh%iNeighbour, dim=2) == nAtom)
    @:ASSERT(size(neigh%nNeighbour) == nAtom)
    @:ASSERT(maxval(neigh%nNeighbour) <= size(neigh%iNeighbour, dim=1))
    @:ASSERT(all(shape(neigh%neighDist2) == shape(neigh%iNeighbour)))
    @:ASSERT(cutoff >= 0.0_dp)

    ! Get last interacting neighbour for given cutoff
    do iAtom = 1, nAtom
      nNeighbourSK(iAtom) = getNrOfNeighbours(neigh, cutoff, iAtom)
    end do

  end subroutine getNrOfNeighboursForAll


  !> Returns the nr. of neighbours for a given atom.
  function getNrOfNeighbours(neigh, cutoff, iAtom) result(nNeighbour)

    !> Intialised neihgborlist.
    type(TNeighbourList), intent(in) :: neigh

    !> Maximal neighbour distance to consider.
    real(dp),            intent(in) :: cutoff

    !> Index of the atom to get the nr. of neighbours for.
    integer, intent(in) :: iAtom

    !> Nr. of neighbours for the specified atom.
    integer :: nNeighbour

    character(len=100) :: strError

    @:ASSERT(cutoff >= 0.0_dp)
    @:ASSERT(iAtom <= size(neigh%nNeighbour))

    ! Issue warning, if cutoff is bigger as used for the neighbourlist.
    if (cutoff > neigh%cutoff) then
99010 format ('Cutoff (', E16.6, ') greater then last cutoff ', '(', E13.6,&
          & ') passed to updateNeighbourList!')
      write (strError, 99010) cutoff, neigh%cutoff
      call warning(strError)
    end if

    ! Get last interacting neighbour for given cutoff
    call bisection(nNeighbour, neigh%neighDist2(1:neigh%nNeighbour(iAtom), iAtom), cutoff**2,&
        & tolSameDist2)

  end function getNrOfNeighbours


  !> Reallocate arrays which depends on the maximal nr. of all atoms.
  subroutine reallocateArrays1(img2CentCell, iCellVec, coord, mNewAtom)

    !> array mapping images of atoms to originals in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index of unit cell containing atom
    integer, allocatable, intent(inout) :: iCellVec(:)

    !> coordinates of all atoms (actual and image)
    real(dp), allocatable, intent(inout) :: coord(:, :)

    !> maximum number of new atoms
    integer, intent(in) :: mNewAtom

    integer :: mAtom
    integer, allocatable :: tmpIntR1(:)
    real(dp), allocatable :: tmpRealR2(:, :)

    mAtom = size(img2CentCell)

    @:ASSERT(size(iCellVec) == mAtom)
    @:ASSERT(all(shape(coord) == (/ 3, mAtom /)))
    !@:ASSERT((mNewAtom > 0) .and. (mNewAtom > mAtom))
    @:ASSERT((mNewAtom > 0))
    mAtom = min(mAtom,mNewAtom)

    call move_alloc(img2CentCell, tmpIntR1)
    allocate(img2CentCell(mNewAtom))
    img2CentCell(:) = 0
    img2CentCell(:mAtom) = tmpIntR1(:mAtom)

    tmpIntR1(:) = iCellVec(:)
    deallocate(iCellVec)
    allocate(iCellVec(mNewAtom))
    iCellVec(:mAtom) = tmpIntR1(:mAtom)

    call move_alloc(coord, tmpRealR2)
    allocate(coord(3, mNewAtom))
    coord(:, :mAtom) = tmpRealR2(:, :mAtom)

  end subroutine reallocateArrays1


  !> Reallocate array which depends on the maximal nr. of neighbours.
  subroutine reallocateArrays3(iNeighbour, neighDist2, mNewNeighbour)

    !> list of neighbours
    integer, allocatable, intent(inout) :: iNeighbour(:, :)

    !> square of distances between atoms
    real(dp), allocatable, intent(inout) :: neighDist2(:,:)

    !> maximum number of new atoms
    integer, intent(in) :: mNewNeighbour

    integer :: mNeighbour, mAtom
    integer, allocatable :: tmpIntR2(:,:)
    real(dp), allocatable :: tmpRealR2(:,:)

    mNeighbour = ubound(iNeighbour, dim=1)
    mAtom = size(iNeighbour, dim=2)

    @:ASSERT(mNewNeighbour > 0 .and. mNewNeighbour > mNeighbour)
    @:ASSERT(all(shape(neighDist2) == shape(iNeighbour)))

    call move_alloc(iNeighbour, tmpIntR2)
    allocate(iNeighbour(0:mNewNeighbour, mAtom))
    iNeighbour(:,:) = 0
    iNeighbour(:mNeighbour, :mAtom) = tmpIntR2

    call move_alloc(neighDist2, tmpRealR2)
    allocate(neighDist2(0:mNewNeighbour, mAtom))
    neighDist2(:,:) = 0.0_dp
    neighDist2(:mNeighbour, :mAtom) = tmpRealR2

  end subroutine reallocateArrays3

  !> Calculate indexing array and number of elements in sparse arrays like the real space overlap
  subroutine getSparseDescriptor(iNeighbour, nNeighbourSK, img2CentCell, orb, iPair, sparseSize)

    !> Neighbours of each atom
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours of each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> Indexing for mapping image atoms to central cell
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Sparse array indexing for the start of atomic blocks in data structures
    integer, allocatable, intent(inout) :: iPair(:,:)

    !> Total number of elements in a sparse structure (ignoring extra indices like spin)
    integer, intent(out) :: sparseSize

    integer :: nAtom, mNeighbour
    integer :: ind, iAt1, nOrb1, iNeigh1, nOrb2

    nAtom = size(iNeighbour, dim=2)
    mNeighbour = size(iNeighbour, dim=1)

    @:ASSERT(allocated(iPair))
    @:ASSERT(size(iPair, dim=2) == nAtom)

    if (mNeighbour > size(iPair, dim=1)) then
      deallocate(iPair)
      allocate(iPair(0 : mNeighbour - 1, nAtom))
      iPair(:,:) = 0
    end if
    ind = 0
    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh1 = 0, nNeighbourSK(iAt1)
        iPair(iNeigh1, iAt1) = ind
        nOrb2 = orb%nOrbAtom(img2CentCell(iNeighbour(iNeigh1, iAt1)))
        ind = ind + nOrb1 * nOrb2
      end do
    end do
    sparseSize = ind

  end subroutine getSparseDescriptor


  !> Allocate (reallocate) space for the sparse hamiltonian and overlap matrix.
  subroutine reallocateHS_1(ham, over, iPair, iNeighbour, nNeighbourSK, orb, img2Centcell)

    !> Hamiltonian
    real(dp), allocatable, intent(inout):: ham(:)

    !> Overlap matrix
    real(dp), allocatable, intent(inout) :: over(:)

    !> Pair indexing array (specifying the offset for the interaction between atoms in the central
    !> cell and their neighbours)
    integer, allocatable, intent(inout) :: iPair(:,:)

    !> List of neighbours for each atom in the central cell. (Note: first index runs from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom in the central cell.
    integer, intent(in) :: nNeighbourSK(:)

    !> Orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> array mapping images of atoms to originals in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> nr. atoms in the central cell
    integer :: nAtom

    !> nr. of elements in the sparse H/S before and after resizing
    integer :: nOldElem, nElem

    !> nr. of max. possible neighbours (incl. itself)
    integer :: mNeighbour

    integer :: ind
    integer :: iAt1, iNeigh1, nOrb1

    nAtom = size(iNeighbour, dim=2)
    mNeighbour = size(iNeighbour, dim=1)
    nOldElem = size(ham, dim=1)

    @:ASSERT(allocated(ham))
    @:ASSERT(allocated(over))
    @:ASSERT(size(over) == nOldElem)
    @:ASSERT(allocated(iPair))
    @:ASSERT(size(iPair, dim=2) == nAtom)

    if (mNeighbour > size(iPair, dim=1)) then
      deallocate(iPair)
      allocate(iPair(0:mNeighbour-1, nAtom))
      iPair(:,:) = 0
    end if
    nElem = 0
    ind = 0
    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh1 = 0, nNeighbourSK(iAt1)
        iPair(iNeigh1, iAt1) = ind
        ind = ind + nOrb1 * orb%nOrbAtom(img2CentCell(iNeighbour(iNeigh1, iAt1)))
      end do
    end do
    nElem = ind
    if (nElem > nOldElem) then
      deallocate(ham)
      deallocate(over)
      allocate(ham(nElem))
      allocate(over(nElem))
      ham(:) = 0.0_dp
      over(:) = 0.0_dp
    end if

  end subroutine reallocateHS_1


  !> Allocate (reallocate) space for the sparse hamiltonian and overlap matrix.
  subroutine reallocateHS_2(ham, over, iPair, iNeighbour, nNeighbourSK, orb, img2CentCell)

    !> Hamiltonian.
    real(dp), allocatable, intent(inout) :: ham(:,:)

    !> Overlap matrix.
    real(dp), allocatable, intent(inout) :: over(:)

    !> Pair indexing array (specifying the offset for the interaction between atoms in the central
    !> cell and their neighbours).
    integer, allocatable, intent(inout) :: iPair(:,:)

    !> List of neighbours for each atom in the central cell. (Note: first index runs from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom in the central cell.
    integer, intent(in) :: nNeighbourSK(:)

    !> Orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Mapping on atoms in the central cell
    integer, intent(in) :: img2CentCell(:)


    !> nr. of spin blocks in the Hamiltonian
    integer :: nSpin

    !> nr. atoms in the central cell
    integer :: nAtom

    !> nr. of elements in the spare H/S
    integer :: nElem, nOldElem

    !> nr. of max. possible neighbours (incl. itself)
    integer :: mNeighbour

    integer :: ind
    integer :: iAt1, iNeigh1, nOrb1

    nAtom = size(iNeighbour, dim=2)
    mNeighbour = size(iNeighbour, dim=1)
    nSpin = size(ham, dim=2)
    nOldElem = size(ham, dim=1)

    @:ASSERT(allocated(ham))
    @:ASSERT(allocated(over))
    @:ASSERT(size(over) == nOldElem)
    @:ASSERT(allocated(iPair))
    @:ASSERT(size(iPair, dim=2) == nAtom)

    if (mNeighbour > size(iPair, dim=1)) then
      deallocate(iPair)
      allocate(iPair(0:mNeighbour-1, nAtom))
      iPair(:,:) = 0
    end if
    nElem = 0
    ind = 0
    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh1 = 0, nNeighbourSK(iAt1)
        iPair(iNeigh1, iAt1) = ind
        ind = ind +  nOrb1 * orb%nOrbAtom(img2CentCell(iNeighbour(iNeigh1,iAt1)))
      end do
    end do
    nElem = ind
    if (nElem > nOldElem) then
      deallocate(ham)
      deallocate(over)
      allocate(ham(nElem, nSpin))
      allocate(over(nElem))
      ham(:,:) = 0.0_dp
      over(:) = 0.0_dp
    end if

  end subroutine reallocateHS_2


  !> Allocate (reallocate) space for the sparse hamiltonian and overlap matrix.
  subroutine reallocateHS_Single(ham, iPair, iNeighbour, nNeighbourSK, orb, img2CentCell)

    !> Hamiltonian.
    real(dp), allocatable, intent(inout) :: ham(:)

    !> Pair indexing array (specifying the offset for the interaction between atoms in the central
    !> cell and their neigbhors).
    integer, allocatable, intent(inout) :: iPair(:,:)

    !> List of neighbours for each atom in the central cell. (Note: first index runs from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom in the central cell.
    integer, intent(in) :: nNeighbourSK(:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Mapping on atoms in the central cell.
    integer, intent(in) :: img2CentCell(:)


    !> nr. atoms in the central cell
    integer :: nAtom

    !> nr. of elements in the spare H/S before and after resizing
    integer :: nOldElem, nElem

    !> nr. of max. possible neighbours (incl. itself)
    integer :: mNeighbour

    integer :: ind
    integer :: iAt1, iNeigh1, nOrb1

    nAtom = size(iNeighbour, dim=2)
    mNeighbour = size(iNeighbour, dim=1)
    nOldElem = size(ham, dim=1)

    @:ASSERT(allocated(ham))
    @:ASSERT(allocated(iPair))
    @:ASSERT(size(iPair, dim=2) == nAtom)

    if (mNeighbour > size(iPair, dim=1)) then
      deallocate(iPair)
      allocate(iPair(0:mNeighbour-1, nAtom))
      iPair(:,:) = 0
    end if
    nElem = 0
    ind = 0
    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh1 = 0, nNeighbourSK(iAt1)
        iPair(iNeigh1, iAt1) = ind
        ind = ind +  nOrb1 * orb%nOrbAtom(img2CentCell(iNeighbour(iNeigh1,iAt1)))
      end do
    end do
    nElem = ind
    if (nElem > nOldElem) then
      deallocate(ham)
      allocate(ham(nElem))
      ham(:) = 0.0_dp
    end if

  end subroutine reallocateHS_Single


  !> Builds an atom offset array for the squared hamiltonain/overlap.
  subroutine buildSquaredAtomIndex(iAtomStart, orb)

    !> Returns the offset array for each atom.
    integer, intent(out) :: iAtomStart(:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    integer :: ind, iAt1
    integer :: nAtom

    nAtom = size(orb%nOrbAtom)

    @:ASSERT(all(shape(iAtomStart) == (/ nAtom + 1 /)))

    ind = 1
    do iAt1 = 1, nAtom
      iAtomStart(iAt1) = ind
      ind = ind + orb%nOrbAtom(iAt1)
    end do
    iAtomStart(nAtom+1) = ind

  end subroutine buildSquaredAtomIndex


  !> Creates a K-points sampling, equivalent to folding of a reciprocal point of a super lattice.
  !> The routine calculates those reciprocal lattice points of the super lattice, which are inside
  !> the Brillouin zone of the original lattice. The resulting points are then all shifted by
  !> sum(shift(i)*B(i)) where B(i) are the reciprocal lattice vectors of the super lattice.
  !> Finally, points equivalent by inversion are reduced, unless specified otherwise.
  subroutine getSuperSampling(coeffs, shifts, kPoints, kWeights, reduceByInversion)

    !> Coefficients of the lattice vectors in the linear combination for the super lattice vectors
    !> (should be integer values)
    real(dp), intent(in) :: coeffs(:,:)

    !> Shift of the grid along the three small reciprocal lattice vectors (between 0.0 and 1.0)
    real(dp), intent(in) :: shifts(:)

    !> Contains the kPoints on exit.
    real(dp), allocatable, intent(out) :: kPoints(:,:)

    !> Contains the weights of the kPoints on exit.
    real(dp), allocatable, intent(out) :: kWeights(:)

    !> If points equivalent by inversion should be reduced.
    logical, intent(in), optional :: reduceByInversion

    real(dp), allocatable :: allKPoints(:,:), allKWeights(:)
    logical, allocatable :: irreducible(:)
    logical :: tReduce
    real(dp) :: invCoeffs(3,3), rr(3)
    integer :: imgRange(2,3), itmp3(3)
    integer :: nAllKPoint, nKPoint
    integer :: i1, i2, i3
    type(TListRealR1) :: lr1

    real(dp), parameter :: tol = 1e-4_dp
    real(dp), parameter :: minLim = -tol, maxLim = 1.0_dp - tol

    @:ASSERT(all(shape(coeffs) == (/ 3, 3 /)))
    ! check they are integers
    @:ASSERT(all(coeffs - nint(coeffs) < epsilon(1.0_dp)))
    @:ASSERT(size(shifts) == 3)

    if (present(reduceByInversion)) then
      tReduce = reduceByInversion
    else
      tReduce = .true.
    end if

    ! Get the eight corners of the original (big) reciprocal unit cell as linear
    ! combination of the reciprocal lattice vectors (B) of the superlattice
    ! Note: b = B * N^T (b/B: rec.lat.vec. of lattice/superlattice)
    imgRange(:,:) = 0
    do i1 = 0, 1
      do i2 = 0, 1
        do i3 = 0, 1
          itmp3 = i1*nint(coeffs(1,:)) + i2*nint(coeffs(2,:)) + i3*nint(coeffs(3,:))
          imgRange(1,:) = min(itmp3, imgRange(1,:))
          imgRange(2,:) = max(itmp3, imgRange(2,:))
        end do
      end do
    end do
    ! Decrease by one to have the range [min, max)
    imgRange(2,:) = imgRange(2,:) - 1

    ! invCoeffs = (N^-1)^T
    call invert33(invCoeffs, coeffs)
    invCoeffs = transpose(invCoeffs)
    call init(lr1)

    do i1 = imgRange(1, 1), imgRange(2, 1)
      do i2 = imgRange(1, 2), imgRange(2, 2)
        do i3 = imgRange(1, 3), imgRange(2, 3)
          ! relative coordinate with respect to the original reciprocal lattice
          rr(:) = matmul(invCoeffs, real((/ i1, i2, i3 /), dp))
          if (all(rr >= minLim) .and. all(rr < maxLim)) then
            ! Add point + shift vector
            call append(lr1, rr + matmul(invCoeffs, shifts))
          end if
        end do
      end do
    end do

    nAllKPoint = len(lr1)
    if (abs(real(nAllKPoint,dp) - abs(determinant33(coeffs))) > tol) then
      call error("Monkhorst-Pack routine failed to find all K-points.")
    end if

    allocate(allKPoints(3, nAllKPoint))
    allocate(allKWeights(nAllKPoint))
    call asArray(lr1, allKPoints)
    call destruct(lr1)
    allKPoints = modulo(allKPoints, 1.0_dp)
    allKWeights = 1.0_dp / real(nAllKPoint, dp)

    ! Reduce by inversion if needed
    if (tReduce) then
      allocate(irreducible(nAllKPoint))
      irreducible(:) = .true.
      do i1 = 1, nAllKPoint
        if (.not. irreducible(i1)) then
          cycle
        end if
        rr(:) = modulo(-1.0_dp * allKPoints(:,i1), 1.0_dp)
        do i2 = i1 + 1, nAllKPoint
          if (.not. irreducible(i2)) then
            cycle
          end if
          if (all(abs(allKPoints(:,i2) - rr(:)) < tol)) then
            irreducible(i2) = .false.
            allKWeights(i1) = allKWeights(i1) + allKWeights(i2)
          end if
        end do
      end do
      nKPoint = count(irreducible)
      allocate(kPoints(3, nKPoint))
      allocate(kWeights(nKPoint))
      i1 = 1
      i2 = 1
      do while (i2 <= nKpoint)
        if (irreducible(i1)) then
          kPoints(:,i2) = allKPoints(:,i1)
          kWeights(i2) = allKWeights(i1)
          i2 = i2 + 1
        end if
        i1 = i1 + 1
      end do
    else
      allocate(kPoints(3, nAllKPoint))
      allocate(kWeights(nAllKPoint))
      kPoints(:,:) = allKPoints
      kWeights(:) = allKWeights
    end if

  end subroutine getSuperSampling


  !> convert fractional coordinates to cartesian
  subroutine fractionalCartesian(cartCoords,latvecs)

    !> fractional coordinates in unit cell on entry, cartesian on exit
    real(dp), intent(inout) :: cartCoords(:,:)

    !> periodic lattice vectors
    real(dp), intent(in) :: latvecs(3,3)

    @:ASSERT(size(cartCoords,dim=1) == 3)

    cartCoords = matmul(latvecs,cartCoords)

  end subroutine fractionalCartesian


  !> Cartesian to fractional coordinates in periodic geometry
  subroutine cartesianFractional(cartCoords,latvecs)

    !> cartesian coordinates on entry, fractional on exit
    real(dp), intent(inout) :: cartCoords(:,:)

    !> periodic lattice vectors
    real(dp), intent(in) :: latvecs(3,3)

    real(dp) :: invLatVecs(3,3)

    @:ASSERT(size(cartCoords,dim=1) == 3)

    call invert33(invLatVecs, latvecs)
    cartCoords = matmul(invLatvecs, cartCoords)

  end subroutine cartesianFractional


  !> Convert from cylindrical to Cartesian coordinate systems
  subroutine cyl2cart_vec(x,r)

    !> Cartesian coordinates
    real(dp), intent(out) :: x(3)

    !> Cylindrical coordinates stored as (radius, height, angle)
    real(dp), intent(in) :: r(3)

    x = 0.0_dp
    x(1) = r(1)*cos(r(3))
    x(2) = r(1)*sin(r(3))
    x(3) = r(2)

  end subroutine cyl2cart_vec


  !> Convert from Cartesian to cylindrical coordinate systems
  subroutine cart2cyl_vec(r,x)

    !> Cartesian coordinates
    real(dp), intent(out) :: r(3)

    !> Cylindrical coordinates stored as (radius, height, angle)
    real(dp), intent(in) :: x(3)

    r = 0.0_dp
    r(1) = sqrt(sum(x(:2)**2))
    r(2) = x(3)
    r(3) = atan2(x(2),x(1))

  end subroutine cart2cyl_vec


  !> Convert from cylindrical to Cartesian coordinate systems
  subroutine cyl2cart_array(x,r)

    !> Cartesian coordinates
    real(dp), intent(out) :: x(:,:)

    !> Cylindrical coordinates stored as (radius, height, angle)
    real(dp), intent(in) :: r(:,:)

  @:ASSERT(all(shape(x)==shape(r)))
  @:ASSERT(size(x,dim=1)==3)
    x = 0.0_dp
    x(1,:) = r(1,:)*cos(r(3,:))
    x(2,:) = r(1,:)*sin(r(3,:))
    x(3,:) = r(2,:)

  end subroutine cyl2cart_array

end module dftbp_periodic
