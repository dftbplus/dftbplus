!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains subroutines for the periodic boundary conditions and neighbour
!!* data
module periodic
  use assert
  use accuracy
  use constants, only : pi
  use message
  use sorting
  use bisect
  use linkedlist
  use simplealgebra, only : determinant33, invert33
  use commontypes
  use memman
  implicit none

  private

  public :: getCellTranslations, getLatticePoints, foldCoordToUnitCell
  public :: reallocateHS, buildSquaredAtomIndex
  public :: TNeighborList, init
  public :: updateNeighborList, updateNeighborListAndSpecies
  public :: getNrOfNeighbors, getNrOfNeighborsForAll
  public :: getImgRange, getSuperSampling
  public :: frac2cart, cart2frac

  interface reallocateHS
    module procedure reallocateHS_1
    module procedure reallocateHS_2
    module procedure reallocateHS_Single
  end interface

  !!* Initializes ADTs defined in this module
  interface init
    module procedure init_TNeighborList
  end interface

  interface frac2cart
    module procedure fractionalCartesian
  end interface

  interface cart2frac
    module procedure cartesianFractional
  end interface

  !!* Contains essential data for the neighborlist
  type TNeighborList
    integer, allocatable :: iNeighbor(:,:)  !* index of neighbor atoms
    integer, allocatable :: nNeighbor(:)    !* nr. of neighbors
    real(dp), allocatable :: neighDist2(:,:) !* temporary array for
                                                   !* neighbor distances
    real(dp)            :: cutoff          !* cutoff it was generated for
    logical             :: initialized = .false.
  end type TNeighborList

contains

  !!* Initializes a neighborlist instance.
  !!* @param neibhorlist   Neighborlist data.
  !!* @param nAtom         Nr. of atoms in the system.
  !!* @param nInitNeighbor Expected nr. of neighbors per atom.
  subroutine init_TNeighborList(neighborList, nAtom, nInitNeighbor)
    type(TNeighborList), intent(out) :: neighborList
    integer,             intent(in)  :: nAtom
    integer,             intent(in)  :: nInitNeighbor

    @:ASSERT(.not. neighborList%initialized)
    @:ASSERT(nAtom > 0)
    @:ASSERT(nInitNeighbor > 0)

    allocate(neighborList%nNeighbor(nAtom))
    allocate(neighborList%iNeighbor(0:nInitNeighbor, nAtom))
    allocate(neighborList%neighDist2(0:nInitNeighbor, nAtom))

    neighborList%cutoff = -1.0_dp
    neighborList%initialized = .true.

  end subroutine init_TNeighborList


  !!* Calculates the translation vectors for cells, which could contain atoms
  !!* interacting with any of the atoms in the central cell.
  !!* @param cellVec  Returns cell translation vectors in relative coordinates.
  !!* @param rCellVec Returns cell translation vectors in absolute units.
  !!* @param latVec   Lattice vectors
  !!* @param recVec2p Reciprocal lattice vectors in 2*pi units.
  !!* @param cutoff   Global cutoff for the diatomic interactions
  !!* @note This subroutine uses a simple guess to get the necessary translation
  !!*   vectors. This results in a set of vectors wich could for very asymmetric
  !!*   cells a large amount bigger than the real necessary one.
  subroutine getCellTranslations(cellVec, rCellVec, latVec, recVec2p, cutoff)
    real(dp), allocatable, intent(out) :: cellVec(:, :)
    real(dp), allocatable, intent(out) :: rCellVec(:,:)
    real(dp), intent(in) :: latVec(:,:)
    real(dp), intent(in) :: recVec2p(:,:)
    real(dp), intent(in) :: cutoff

    integer  :: ii

    @:ASSERT(all(shape(latVec) == [3, 3]))
    @:ASSERT(all(shape(recVec2p) == [3, 3]))
    @:ASSERT(cutoff >= 0.0_dp)

    call getLatticePoints(cellVec, latVec, recVec2p, cutoff, posExtension=1, &
        &negExtension=1)
    allocate(rCellVec(3, size(cellVec, dim=2)))
    do ii = 1, size(rCellVec, dim=2)
      rCellVec(:,ii) = matmul(latVec, cellVec(:,ii))
    end do

  end subroutine getCellTranslations



  subroutine getImgRange(imgRange, dist, recVec2p, posExt, negExt)
    integer, intent(out) :: imgRange(:,:)
    real(dp), intent(in) :: dist
    real(dp), intent(in) :: recVec2p(:,:)
    integer, intent(in) :: posExt
    integer, intent(in) :: negExt

    integer :: ii, iTmp

    @:ASSERT(all(shape(imgRange) == (/ 2, 3 /)))
    @:ASSERT(dist >= 0.0_dp)
    @:ASSERT(all(shape(recVec2p) == (/ 3, 3 /)))

    do ii = 1, 3
      iTmp = floor(dist * sqrt(sum(recVec2p(:, ii)**2)))
      imgRange(1, ii) = -(iTmp + negExt)
      imgRange(2, ii) = iTmp + posExt
    end do

  end subroutine getImgRange



  !!* Returns a set which definitely contains all the points of a 3D grid
  !!* which are nearer to the origin as a given distance.
  !!* @param latPoint     Returns grid points in relative coords.
  !!* @param latVec       Lattice vectors.
  !!* @param recVec2p     Reciprocal lattice vectors in 2*pi units.
  !!* @param dist         Global cutoff for the diatomic interactions.
  !!* @param posExtension Extend the set along the positive lattice vectors with
  !!*    that many additional lattice vectors.
  !!* @param negExtension Same as posExtension for negative lattice vectors
  !!* @param onlyInside   Return only those lattice points which are really
  !!*                     not outside the given distance.
  !!* @note Without the onlyInside parameter, the returned set of lattice points
  !!*   shape a parallelepipedon. With the onlyInside parameter its not
  !!*   necessarily the case.
  !!* @todo Refine the algorithm with the help of a new routine which can
  !!*   calculate the minimal distance between two arbitary cells.
  subroutine getLatticePoints(latPoint, latVec, recVec2p, dist, posExtension, &
      &negExtension, onlyInside, reduceByInversion, withoutOrigin)
    real(dp), allocatable, intent(out) :: latPoint(:,:)
    real(dp), intent(in) :: latVec(:,:)
    real(dp), intent(in) :: recVec2p(:,:)
    real(dp), intent(in) :: dist
    integer,  intent(in), optional :: posExtension
    integer,  intent(in), optional :: negExtension
    logical,  intent(in), optional :: onlyInside
    logical,  intent(in), optional :: reduceByInversion
    logical,  intent(in), optional :: withoutOrigin

    integer  :: imgRange(2, 3)
    integer  :: negExt, posExt
    integer  :: ii, jj, kk, ind, ii0, jj0, kk0
    real(dp) :: ri, rj, rk, rTmp, dist2
    logical  :: tAll, tOrig, tNoInv
    real(dp), allocatable :: tmpLatPoint(:,:)

    @:ASSERT(all(shape(latVec) == (/3, 3/)))
    @:ASSERT(all(shape(recVec2p) == (/3, 3/)))
    @:ASSERT(dist >= 0.0_dp)

    if (present(negExtension)) then
      negExt = negExtension
    else
      negExt = 0
    end if
    if (present(posExtension)) then
      posExt = posExtension
    else
      posExt = 0
    end if
    if (present(onlyInside)) then
      tAll = .not. onlyInside
    else
      tAll = .true.
    end if
    if (present(withoutOrigin)) then
      tOrig = .not. withoutOrigin
    else
      tOrig = .true.
    end if
    if (present(reduceByInversion)) then
      tNoInv = reduceByInversion
    else
      tNoInv = .false.
    end if

    call getImgRange(imgRange, dist, recVec2p, posExt, negExt)

    allocate(tmpLatPoint(3, product(sum(abs(imgRange), dim=1) + 1)))

    dist2 = dist**2
    if (tOrig) then
      tmpLatPoint(:,1) = [0.0_dp, 0.0_dp, 0.0_dp]
      ind = 1
    else
      ind = 0
    end if

    if (tNoInv) then
      ii0 = 0
    else
      ii0 = imgRange(1,1)
    end if
    iiLoop: do ii = ii0, imgRange(2, 1)
      ri = real(ii)

      if (tNoInv .and. ii == 0) then
        jj0 = 0
      else
        jj0 = imgRange(1, 2)
      end if
      jjLoop: do jj = jj0, imgRange(2, 2)
        rj = real(jj)

        if (tNoInv .and. ii == 0 .and. jj == 0) then
          kk0 = 0
        else
          kk0 = imgRange(1, 3)
        end if
        kkLoop: do kk = kk0, imgRange(2, 3)
          if (kk == 0 .and. jj == 0 .and. ii == 0) then
            cycle
          end if
          rk = real(kk)
          if (tAll) then
            ind = ind + 1
            tmpLatPoint(:,ind) = [ri, rj, rk]
          else
            rTmp = sum(matmul(latVec, reshape([ri, rj, rk], [3, 1]))**2)
            if (rTmp <= dist2) then
              ind = ind + 1
              tmpLatPoint(:,ind) = [ri, rj, rk]
            end if
          end if
        end do kkLoop
      end do jjLoop
    end do iiLoop

    if (ind == size(tmpLatPoint, dim=2)) then
      call move_alloc(tmpLatPoint, latPoint)
    else
      allocate(latPoint(3, ind))
      latPoint(:,:) = tmpLatPoint(:,:ind)
    end if

  end subroutine getLatticePoints



  !!* Fold coordinates back in the central cell
  !!* @param coord       Contains the original coordinates on call and
  !!*   the folded ones on return.
  !!* @param latVec      Lattice vectors (column format).
  !!* @param recVec2p    Reciprocal vectors in units of 2pi (column format).
  !!* @param invShift    Contains difference vectors old_coords - new_coords.
  !!* @desc
  !!*   Throw away the integer part of the relative coordinates of every
  !!*   atom. If the resulting coordinate is very near to 1.0 (closer
  !!*   than 1e-12 in absolute length), fold it to 0.0 to make the
  !!*   algorithm more predictable and independent of numeric noises.
  subroutine foldCoordToUnitCell(coord, latVec, recVec2p, invShift)
    real(dp), intent(inout) :: coord(:,:)
    real(dp), intent(in)    :: latVec(:,:)
    real(dp), intent(in)    :: recVec2p(:,:)
    real(dp), intent(out), optional :: invShift(:,:)

    integer  :: nAtom             ! Nr. of atoms in the cell.
    integer  :: ii, jj
    real(dp) :: frac(3), frac2(3), tmp3(3), vecLen(3)

    nAtom = size(coord, dim=2)

    @:ASSERT(size(coord, dim=1) == 3)
    @:ASSERT(all(shape(latVec) == (/3, 3/)))
    @:ASSERT(all(shape(recVec2p) == (/3, 3/)))
  #:call ASSERT_CODE
    if (present(invShift)) then
      ASSERT(all(shape(invShift) == shape(coord)))
    end if
  #:endcall ASSERT_CODE

    vecLen(:) = sqrt(sum(latVec(:,:)**2, dim=1))
    do ii = 1, nAtom
      do jj = 1, 3
        frac(jj) = dot_product(recVec2p(:,jj), coord(:,ii))
      end do
      tmp3(:) = coord(:,ii)
      frac2(:) = frac(:) - real(floor(frac(:)), dp)
      where (abs(vecLen*(1.0_dp - frac2)) < 1e-12_dp) frac2 = 0.0_dp
      coord(:, ii) = matmul(latVec, frac2)
      if (present(invShift)) then
        invShift(:,ii) = tmp3(:) - coord(:,ii)
      end if
    end do

  end subroutine foldCoordToUnitCell



  !!* Updates the neighbor list and the species arrays.
  !!* @param coord Coordinates of all interacting atoms on exit
  !!* @param species Species of all interacting atoms on exit.
  !!* @param img2CentCell Mapping on atoms in the central cell
  !!* @param iCellVec Shift vector index for every interacting atom
  !!* @param neig Updated neighbor list.
  !!* @param nAllAtom Number of all interacting atoms
  !!* @param coord0 Coordinates of the atoms in the central cell
  !!* @param species0 Species of the atoms in the central cell
  !!* @param cutoff Cutoff until neighborlist should be created
  !!* @param rCellVec Cell vector for the translated cells to consider.
  subroutine updateNeighborListAndSpecies(coord, species, img2CentCell, &
      &iCellVec, neigh, nAllAtom, coord0, species0, cutoff, rCellVec)
    real(dp), allocatable, intent(inout) :: coord(:,:)
    integer,  allocatable, intent(inout) :: species(:)
    integer,  allocatable, intent(inout) :: img2CentCell(:)
    integer,  allocatable, intent(inout) :: iCellVec(:)
    type(TNeighborList), intent(inout) :: neigh
    integer,  intent(out)              :: nAllAtom
    real(dp), intent(in)               :: coord0(:,:)
    integer,  intent(in)               :: species0(:)
    real(dp), intent(in)               :: cutoff
    real(dp), intent(in)               :: rCellVec(:,:)


    call updateNeighborList(coord, img2CentCell, iCellVec, neigh, nAllAtom, &
        &coord0, cutoff, rCellVec)
    if (size(species) < nAllAtom) then
      deallocate(species)
      allocate(species(nAllAtom))
    end if
    species(1:nAllAtom) = species0(img2CentCell(1:nAllAtom))

  end subroutine updateNeighborListAndSpecies



  !!* Updates the neighbor list according a given geometry.
  !!* @param coord1       Coordinates of the objects interacting with the
  !!*   objects in the central cell (on exit).
  !!* @param img2CentCell Returns for all objects the index of an object in
  !!*   the central cell which the object is mapped on.
  !!* @param iCellVec     Returns the index of the translating
  !!*   superlattice vector for each object.
  !!* @param neigh        Neighborlist.
  !!* @param nAllAtom     Returns the nr. of all objects (including those
  !!*   in the translated cells.)
  !!* @param coord0In     Coordinates of the objects in the central cell.
  !!* @param coord1In     Coorindates of the periodical images interacting
  !!*   with the objects in the centrall cell.
  !!* @param cutoff       Cutoff radius for the interactions.
  !!* @param rCellVec     Absolute coordinates of the shifted supercells which
  !!*   could have interacting atoms with the central cell.
  !!* @desc The neighborlist for the given cutoff is calculated. Arrays are
  !!*   resized if necessary. The neighbor list determination is a simple
  !!*   N^2 algorithm, calculating the distance between the possible atom pairs.
  subroutine updateNeighborList(coord, img2CentCell, iCellVec, neigh, &
      &nAllAtom, coord0, cutoff, rCellVec)
    real(dp), allocatable, intent(inout) :: coord(:,:)
    integer,  allocatable, intent(inout) :: img2CentCell(:)
    integer,  allocatable, intent(inout) :: iCellVec(:)
    type(TNeighborList), intent(inout) :: neigh
    integer,  intent(out)              :: nAllAtom
    real(dp), intent(in)               :: coord0(:,:)
    real(dp), intent(in)               :: cutoff
    real(dp), intent(in)               :: rCellVec(:,:)

    integer  :: nAtom             ! Nr. of atoms in the system
    integer  :: mAtom             ! Max. nr. of atom without reallocation
    integer  :: maxNeighbor       ! Max. nr. of neighbors without reallocation
    integer  :: nCellVec          ! Nr. of cell translation vectors

    !! Square of the diatomic interaction cutoffs
    real(dp) :: cutoff2

    real(dp) :: dist2
    real(dp) :: rCell(3), rr(3)
    integer  :: ii, iAtom1, oldIAtom1, iAtom2
    integer  :: nn1

    integer,  allocatable :: indx(:)
    character(len=100)    :: strError

    nAtom = size(neigh%nNeighbor, dim=1)
    mAtom = size(coord, dim=2)
    maxNeighbor = ubound(neigh%iNeighbor, dim=1)
    nCellVec = size(rCellVec, dim=2)

    @:ASSERT(nAtom <= mAtom)
    @:ASSERT(allocated(coord))
    @:ASSERT(size(coord, dim=1) == 3)
    @:ASSERT(allocated(img2CentCell))
    @:ASSERT(size(img2CentCell) == mAtom)
    @:ASSERT(allocated(iCellVec))
    @:ASSERT(size(iCellVec) == mAtom)
    @:ASSERT(size(neigh%iNeighbor, dim=2) == nAtom)
    @:ASSERT((size(coord0, dim=1) == 3) .and. size(coord0, dim=2) >= nAtom)
    @:ASSERT((size(rCellVec, dim=1) == 3))
    @:ASSERT(cutoff >= 0.0_dp)

    neigh%cutoff = cutoff
    cutoff2 = cutoff**2
    nAllAtom = 0

    !! Clean arrays.
    !! (Every atom is the 0th neighbor of itself with zero distance square.)
    neigh%nNeighbor(:) = 0
    neigh%iNeighbor(:,:) = 0
    do ii = 1, nAtom
      neigh%iNeighbor(0, ii) = ii
    end do
    neigh%neighDist2(:,:) = 0.0_dp


    !! Loop over all possible neighbors for all atoms in the central cell.
    !! Only those neighbors are considered which map on atom with a higher
    !! or equal index in the central cell.
    !! Outer two loops: all atoms in all cells.
    !! Inner loop: all atoms in the central cell.
    lpCellVec: do ii = 1, nCellVec
      rCell(:) = rCellVec(:, ii)
      oldIAtom1 = 0
      lpIAtom1: do iAtom1 = 1, nAtom
        rr(:) = coord0(:, iAtom1) + rCell(:)
        lpIAtom2: do iAtom2 = 1, iAtom1

          !! If distance greater than cutoff -> skip
          dist2 = sum((coord0(:, iAtom2) - rr(:))**2)
          if (dist2 > cutoff2) then
            cycle
          end if

          !! New interacting atom -> append
          !! We need that before checking for interaction with dummy atom or
          !! with itself to make sure that atoms in the central cell are
          !! appended  exactly in the same order as found in the coord0 array.
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

          !! Check if atoms are not too close to each other
          if (dist2 < minNeighDist2) then
            if (ii == 1 .and. iAtom1 == iAtom2) then
              !! We calculated the distance between the same atom in the
              !! unit cell
              cycle
            else
99000         format ('Atoms ',I5,' and ',I5,' too close to each other!', &
                  & ' (dist=',E13.6,')')
              write (strError, 99000) iAtom2, nAllAtom, sqrt(dist2)
              call warning(strError)
            end if
          end if

          neigh%nNeighbor(iAtom2) = neigh%nNeighbor(iAtom2) + 1
          if (neigh%nNeighbor(iAtom2) > maxNeighbor) then
            maxNeighbor = incrmntOfArray(maxNeighbor)
            call reallocateArrays3(neigh%iNeighbor, neigh%neighDist2, &
                & maxNeighbor)

          end if
          neigh%iNeighbor(neigh%nNeighbor(iAtom2), iAtom2) = nAllAtom
          neigh%neighDist2(neigh%nNeighbor(iAtom2), iAtom2) = dist2

        end do lpIAtom2
      end do lpIAtom1
    end do lpCellVec

    !! Sort neighbors for all atom by distance
    allocate(indx(maxNeighbor))
    do iAtom1 = 1, nAtom
      nn1 = neigh%nNeighbor(iAtom1)
      call index_heap_sort(indx(1:nn1), neigh%neighDist2(1:nn1, iAtom1),&
          &tolSameDist2)
      neigh%iNeighbor(1:nn1, iAtom1) = neigh%iNeighbor(indx(:nn1), iAtom1)
      neigh%neighDist2(1:nn1, iAtom1) = neigh%neighDist2(indx(:nn1), iAtom1)
    end do
    coord(:,nAllAtom+1:size(coord, dim=2)) = 0.0_dp

  end subroutine updateNeighborList



  !!* Returns the nr. of neighbors for a given cutoff for all atoms.
  !!* @param nNeighbor Contains the nr. of neighbors for each atom on exit.
  !!* @param neigh     Initialized neighborlist
  !!* @param cutoff    Maximal neighbor distance to consider.
  subroutine getNrOfNeighborsForAll(nNeighbor, neigh, cutoff)
    integer,             intent(out) :: nNeighbor(:)
    type(TNeighborList), intent(in)  :: neigh
    real(dp),            intent(in)  :: cutoff

    integer  :: nAtom, iAtom

    nAtom = size(nNeighbor)

    @:ASSERT(size(neigh%iNeighbor, dim=2) == nAtom)
    @:ASSERT(size(neigh%nNeighbor) == nAtom)
    @:ASSERT(maxval(neigh%nNeighbor) <= size(neigh%iNeighbor, dim=1))
    @:ASSERT(all(shape(neigh%neighDist2) == shape(neigh%iNeighbor)))
    @:ASSERT(cutoff >= 0.0_dp)

    !! Get last interacting neighbor for given cutoff
    do iAtom = 1, nAtom
      nNeighbor(iAtom) = getNrOfNeighbors(neigh, cutoff, iAtom)
    end do

  end subroutine getNrOfNeighborsForAll



  !!* Returns the nr. of neighbors for a given atom.
  !!* @param neigh   Intialized neihgborlist.
  !!* @param cutoff  Maximal neighbor distance to consider.
  !!* @param iAtom   Index of the atom to get the nr. of neighbors for.
  !!* @return        Nr. of neighbors for the specified atom.
  function getNrOfNeighbors(neigh, cutoff, iAtom) result(nNeighbor)
    type(TNeighborList), intent(in)  :: neigh
    real(dp),            intent(in)  :: cutoff
    integer,             intent(in)  :: iAtom
    integer :: nNeighbor

    character(len=100) :: strError

    @:ASSERT(cutoff >= 0.0_dp)
    @:ASSERT(iAtom <= size(neigh%nNeighbor))

    !! Issue warning, if cutoff is bigger as used for the neighborlist.
    if (cutoff > neigh%cutoff) then
99010 format ('Cutoff (', E16.6, ') greater then last cutoff ', &
          & '(', E13.6, ') passed to updateNeighborList!')
      write (strError, 99010) cutoff, neigh%cutoff
      call warning(strError)
    end if

    !! Get last interacting neighbor for given cutoff
    call bisection(nNeighbor, &
        &neigh%neighDist2(1:neigh%nNeighbor(iAtom), iAtom), cutoff**2, &
        &tolSameDist2)

  end function getNrOfNeighbors





  !!* Reallocate arrays which depends on the maximal nr. of all atoms.
  subroutine reallocateArrays1(img2CentCell, iCellVec, coord, mNewAtom)
    integer,  allocatable, intent(inout) :: img2CentCell(:)
    integer,  allocatable, intent(inout) :: iCellVec(:)
    real(dp), allocatable, intent(inout) :: coord(:, :)
    integer,  intent(in) :: mNewAtom

    integer               :: mAtom
    integer, allocatable :: tmpIntR1(:)
    real(dp), allocatable :: tmpRealR2(:, :)

    mAtom = size(img2CentCell)

    @:ASSERT(size(iCellVec) == mAtom)
    @:ASSERT(all(shape(coord) == (/ 3, mAtom /)))
    @:ASSERT((mNewAtom > 0) .and. (mNewAtom > mAtom))

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
    coord(:, :mAtom) = tmpRealR2

  end subroutine reallocateArrays1



  !!* Reallocate array which depends on the maximal nr. of neighbors.
  subroutine reallocateArrays3(iNeighbor, neighDist2, mNewNeighbor)
    integer, allocatable, intent(inout) :: iNeighbor(:, :)
    real(dp), allocatable, intent(inout) :: neighDist2(:,:)
    integer,  intent(in) :: mNewNeighbor

    integer               :: mNeighbor, mAtom
    integer,  allocatable :: tmpIntR2(:,:)
    real(dp), allocatable :: tmpRealR2(:,:)

    mNeighbor = ubound(iNeighbor, dim=1)
    mAtom = size(iNeighbor, dim=2)

    @:ASSERT(mNewNeighbor > 0 .and. mNewNeighbor > mNeighbor)
    @:ASSERT(all(shape(neighDist2) == shape(iNeighbor)))

    call move_alloc(iNeighbor, tmpIntR2)
    allocate(iNeighbor(0:mNewNeighbor, mAtom))
    iNeighbor(:,:) = 0
    iNeighbor(:mNeighbor, :mAtom) = tmpIntR2

    call move_alloc(neighDist2, tmpRealR2)
    allocate(neighDist2(0:mNewNeighbor, mAtom))
    neighDist2(:,:) = 0.0_dp
    neighDist2(:mNeighbor, :mAtom) = tmpRealR2

  end subroutine reallocateArrays3



  !!* Allocate (reallocate) space for the sparse hamiltonian and overlap
  !!* matrix.
  !!* @param ham       Hamiltonian.
  !!* @param over      Overlap matrix.
  !!* @param iPair     Pair indexing array (specifying the offset for the
  !!*   interaction between atoms in the central cell and their neigbhors).
  !!* @param iNeighbor List of neighbors for each atom in the central cell.
  !!*   (Note: first index runs from 0!)
  !!* @param nNeighbor Nr. of neighbors for each atom in the central cell.
  !!* @param orb Orbitals in the system.
  subroutine reallocateHS_1(ham, over, iPair, iNeighbor, nNeighbor, orb, &
      &img2Centcell)
    real(dp), allocatable, intent(inout):: ham(:)
    real(dp), allocatable, intent(inout) :: over(:)
    integer, allocatable, intent(inout) :: iPair(:,:)
    integer,  intent(in) :: iNeighbor(0:,:)
    integer,  intent(in) :: nNeighbor(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: img2CentCell(:)

    integer :: nAtom           ! nr. atoms in the central cell
    integer :: nElem, nOldElem ! nr. of elements in the spare H/S
    integer :: mNeighbor       ! nr. of max. possible neighbors (incl. itself)
    integer :: ind
    integer :: iAt1, iNeigh1, nOrb1

    nAtom = size(iNeighbor, dim=2)
    mNeighbor = size(iNeighbor, dim=1)
    nOldElem = size(ham, dim=1)

    @:ASSERT(allocated(ham))
    @:ASSERT(allocated(over))
    @:ASSERT(size(over) == nOldElem)
    @:ASSERT(allocated(iPair))
    @:ASSERT(size(iPair, dim=2) == nAtom)

    if (mNeighbor > size(iPair, dim=1)) then
      deallocate(iPair)
      allocate(iPair(0:mNeighbor-1, nAtom))
      iPair(:,:) = 0
    end if
    nElem = 0
    ind = 0
    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh1 = 0, nNeighbor(iAt1)
        iPair(iNeigh1, iAt1) = ind
        ind = ind + nOrb1 * orb%nOrbAtom(img2CentCell(iNeighbor(iNeigh1, iAt1)))
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



  !!* Allocate (reallocate) space for the sparse hamiltonian and overlap
  !!* matrix.
  !!* @param ham       Hamiltonian.
  !!* @param over      Overlap matrix.
  !!* @param iPair     Pair indexing array (specifying the offset for the
  !!*   interaction between atoms in the central cell and their neigbhors).
  !!* @param iNeighbor List of neighbors for each atom in the central cell.
  !!*   (Note: first index runs from 0!)
  !!* @param nNeighbor Nr. of neighbors for each atom in the central cell.
  !!* @param orb Orbitals in the system.
  !!* @param img2CentCell Mapping on atoms in the central cell
  subroutine reallocateHS_2(ham, over, iPair, iNeighbor, nNeighbor, orb, &
      &img2CentCell)
    real(dp), allocatable, intent(inout) :: ham(:,:)
    real(dp), allocatable, intent(inout) :: over(:)
    integer, allocatable, intent(inout) :: iPair(:,:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer,  intent(in) :: nNeighbor(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: img2CentCell(:)

    integer :: nSpin           ! nr. of spin blocks in the Hamiltonian
    integer :: nAtom           ! nr. atoms in the central cell
    integer :: nElem, nOldElem ! nr. of elements in the spare H/S
    integer :: mNeighbor       ! nr. of max. possible neighbors (incl. itself)
    integer :: ind
    integer :: iAt1, iNeigh1, nOrb1

    nAtom = size(iNeighbor, dim=2)
    mNeighbor = size(iNeighbor, dim=1)
    nSpin = size(ham, dim=2)
    nOldElem = size(ham, dim=1)

    @:ASSERT(allocated(ham))
    @:ASSERT(allocated(over))
    @:ASSERT(size(over) == nOldElem)
    @:ASSERT(allocated(iPair))
    @:ASSERT(size(iPair, dim=2) == nAtom)

    if (mNeighbor > size(iPair, dim=1)) then
      deallocate(iPair)
      allocate(iPair(0:mNeighbor-1, nAtom))
      iPair(:,:) = 0
    end if
    nElem = 0
    ind = 0
    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh1 = 0, nNeighbor(iAt1)
        iPair(iNeigh1, iAt1) = ind
        ind = ind +  nOrb1 * orb%nOrbAtom(img2CentCell(iNeighbor(iNeigh1,iAt1)))
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

  !!* Allocate (reallocate) space for the sparse hamiltonian and overlap
  !!* matrix.
  !!* @param ham       Hamiltonian.
  !!* @param iPair     Pair indexing array (specifying the offset for the
  !!*   interaction between atoms in the central cell and their neigbhors).
  !!* @param iNeighbor List of neighbors for each atom in the central cell.
  !!*   (Note: first index runs from 0!)
  !!* @param nNeighbor Nr. of neighbors for each atom in the central cell.
  !!* @param orb Information about the orbitals in the system.
  !!* @param img2CentCell Mapping on atoms in the central cell.
  subroutine reallocateHS_Single(ham, iPair, iNeighbor, nNeighbor, orb, &
      &img2CentCell)
    real(dp), allocatable, intent(inout) :: ham(:)
    integer, allocatable, intent(inout) :: iPair(:,:)
    integer,  intent(in) :: iNeighbor(0:,:)
    integer,  intent(in) :: nNeighbor(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: img2CentCell(:)

    integer :: nAtom           ! nr. atoms in the central cell
    integer :: nElem, nOldElem ! nr. of elements in the spare H/S
    integer :: mNeighbor       ! nr. of max. possible neighbors (incl. itself)
    integer :: ind
    integer :: iAt1, iNeigh1, nOrb1

    nAtom = size(iNeighbor, dim=2)
    mNeighbor = size(iNeighbor, dim=1)
    nOldElem = size(ham, dim=1)

    @:ASSERT(allocated(ham))
    @:ASSERT(allocated(iPair))
    @:ASSERT(size(iPair, dim=2) == nAtom)

    if (mNeighbor > size(iPair, dim=1)) then
      deallocate(iPair)
      allocate(iPair(0:mNeighbor-1, nAtom))
      iPair(:,:) = 0
    end if
    nElem = 0
    ind = 0
    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh1 = 0, nNeighbor(iAt1)
        iPair(iNeigh1, iAt1) = ind
        ind = ind +  nOrb1 * orb%nOrbAtom(img2CentCell(iNeighbor(iNeigh1,iAt1)))
      end do
    end do
    nElem = ind
    if (nElem > nOldElem) then
      deallocate(ham)
      allocate(ham(nElem))
      ham(:) = 0.0_dp
    end if

  end subroutine reallocateHS_Single



  !!* Builds an atom offset array for the squared hamiltonain/overlap.
  !!* @param iAtomStart  Returns the offset array for each atom.
  !!* @param orb Information about the orbitals in the system.
  subroutine buildSquaredAtomIndex(iAtomStart, orb)
    integer, intent(out) :: iAtomStart(:)
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

  !!* Creates a K-points sampling, equivalent to folding of a reciprocal point
  !!*   of a super lattice.
  !!* @param coeffs Coefficients of the lattice vectors in the linear
  !!*   combination for the super lattice vectors (should be integer values)
  !!* @param shifts Shift of the grid along the three small reciprocal lattice
  !!*   vectors (between 0.0 and 1.0)
  !!* @param latVecs Lattice vector of the original grid
  !!* @param recVecs2p Reciprocal lattice vectors in 2p units (inverse(latVecs))
  !!* @param kPoints Contains the kPoints on exit.
  !!* @param kWeights Contains the weights of the kPoints on exit.
  !!* @param reduceByInversion If points equivalent by inversion should be
  !!*   reduced.
  !!* @desc The routine calculates those reciprocal lattice points of the
  !!*   super lattice, which are inside the Brillouin zone of the original
  !!*   lattice. The resulting points are then all shifted by sum(shift(i)*B(i))
  !!*   where B(i) are the reciprocal lattice vectors of the super lattice.
  !!*   Finally, points equivalent by inversion are reduced, unless specified
  !!*   otherwise.
  subroutine getSuperSampling(coeffs, shifts, kPoints, kWeights, reduceByInversion)
    real(dp), intent(in) :: coeffs(:,:)
    real(dp), intent(in) :: shifts(:)
    real(dp), allocatable, intent(out) :: kPoints(:,:)
    real(dp), allocatable, intent(out) :: kWeights(:)
    logical, intent(in), optional :: reduceByInversion

    real(dp), allocatable :: allKPoints(:,:), allKWeights(:)
    logical, allocatable :: irreducible(:)
    logical :: tReduce
    real(dp) :: invCoeffs(3,3), rr(3)
    integer :: imgRange(2,3), itmp3(3)
    integer :: nAllKPoint, nKPoint
    integer :: i1, i2, i3
    type(listRealR1) :: lr1

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
          itmp3 = i1*nint(coeffs(1,:)) + i2*nint(coeffs(2,:)) &
              & + i3*nint(coeffs(3,:))
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


  subroutine fractionalCartesian(cartCoords,latvecs)
    real(dp), intent(inout) :: cartCoords(:,:)
    real(dp), intent(in)  :: latvecs(3,3)

    @:ASSERT(size(cartCoords,dim=1) == 3)

    cartCoords = matmul(latvecs,cartCoords)

  end subroutine fractionalCartesian

  subroutine cartesianFractional(cartCoords,latvecs)
    real(dp), intent(inout) :: cartCoords(:,:)
    real(dp), intent(in)  :: latvecs(3,3)

    real(dp) :: invLatVecs(3,3)


    @:ASSERT(size(cartCoords,dim=1) == 3)

    call invert33(invLatVecs, latvecs)
    cartCoords = matmul(invLatvecs, cartCoords)

  end subroutine cartesianFractional

end module periodic
