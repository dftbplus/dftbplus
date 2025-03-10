!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Contains subroutines for the periodic boundary conditions and neighbour data
module dftbp_dftb_periodic
  use dftbp_common_accuracy, only : dp, tolSameDist, tolSameDist2, minNeighDist, minNeighDist2
  use dftbp_common_constants, only : pi
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_memman, only : incrmntOfArray, TAlignedArray
  use dftbp_common_schedule, only : assembleChunks, distributeRangeInChunks, getChunkRanges
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_boundarycond, only : zAxis
#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : mpifx_win, mpifx_allreduceip, mpifx_allgather, MPI_MAX, MPI_LOR,&
      & MPIFX_SIZE_T
#:endif
  use dftbp_io_message, only : error, warning
  use dftbp_math_bisect, only : bisection
  use dftbp_math_quaternions, only : rotate3
  use dftbp_math_simplealgebra, only : determinant33, invert33
  use dftbp_math_sorting, only : index_heap_sort
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_latpointiter, only : TLatPointIter, TLatPointIter_init
  use dftbp_type_linkedlist, only : TListRealR1, len, init, append, asArray, destruct

  implicit none
  private

  public :: getCellTranslations, getLatticePoints
  public :: getSuperSampling
  public :: frac2cart, cart2frac
  public :: TNeighbourList, TNeighbourList_init, TSymNeighbourList
  public :: updateNeighbourList, updateNeighbourListAndSpecies, setNeighbourList
  public :: getNrOfNeighbours, getNrOfNeighboursForAll

#:if WITH_UNIT_TESTS
  ! NOTE: these entries are public only temporarily for unit testing purposes. Do not call them
  ! from the outside.
  public :: distributeAtoms, reallocateArrays2, allocateNeighbourArrays, fillNeighbourArrays
#:endif

  !> Contains essential data for the neighbourlist
  type TNeighbourList

    !> Number of neighbours
    integer, allocatable :: nNeighbour(:)

    !> Index of neighbour atoms
    integer, pointer :: iNeighbour(:,:) => null()

    !> Neighbour distances
    real(dp), pointer :: neighDist2(:,:) => null()

    !> Cutoff it was generated for
    real(dp) :: cutoff

    !> Initialised data
    logical :: initialized = .false.

    !> Whether the neighbour list has been set by an API call
    logical :: setExternally = .false.

    !> Whether memory should be allocated via MPI-windows (or directly via allocate() otherwise)
    !!
    !! Note: this variable cannot be inside the MPI block below with current code, as the parser
    !! uses this data structure before the MPI environment is set up (in addition to the main code
    !! use). This leads to a fall back to the usual distributed case if MPI is enabled, instead of
    !! the shared window in that case.
    logical, private :: useMpiWindows_ = .false.

    !> Memory allocated for the iNeighbour array
    integer, pointer, private :: iNeighbourMem_(:) => null()

    !> Memory allocated for the neighDist2 array
    real(dp), pointer, private :: neighDist2Mem_(:) => null()

  #:if WITH_MPI

    !> MPI shared memory window handler for iNeighbour
    type(mpifx_win), private :: iNeighbourWin_

    !> MPI shared memory window handler for neightDist2
    type(mpifx_win), private :: neighDist2Win_

  #:endif

  contains

    final :: TNeighbourList_final

  end type TNeighbourList


  !> Contains neighbour list instance and symmetry specific entries
  type TSymNeighbourList

    !> Neighbour list instance
    type(TNeighbourList), allocatable :: neighbourList

    !> Number of all interacting atoms, including periodic images
    integer :: nAllAtom

    !> Coordinates of all interacting atoms, including periodic images
    real(dp), allocatable :: coord(:,:)

    !> Species of all interacting atoms, including periodic images
    integer, allocatable :: species(:)

    !> Mapping of all atoms onto atoms in the central cell
    integer, allocatable :: img2CentCell(:)

    !> Shift vector index for every interacting atom, including periodic images
    integer, allocatable :: iCellVec(:)

    !> Sparse array indexing for the start of atomic blocks in data structures
    integer, allocatable :: iPair(:,:)

    !> Total number of elements in a sparse structure (ignoring extra indices like spin)
    integer :: sparseSize

  contains

    final :: TSymNeighbourList_final

  end type TSymNeighbourList

contains


  !> Initializes a neighbourlist instance.
  subroutine TNeighbourList_init(this, nAtom, nInitNeighbour)

    !> Neighbourlist data.
    type(TNeighbourList), intent(out) :: this

    !> Nr. of atoms in the system.
    integer, intent(in) :: nAtom

    !> Expected nr. of neighbours per atom.
    integer, intent(in) :: nInitNeighbour

    @:ASSERT(.not. this%initialized)
    @:ASSERT(nAtom > 0)
    @:ASSERT(nInitNeighbour > 0)

    allocate(this%nNeighbour(nAtom))

    this%cutoff = -1.0_dp
    this%initialized = .true.
    this%setExternally = .false.

  end subroutine TNeighbourList_init


  !> Deallocates MPI shared memory if required
  subroutine TNeighbourList_final(this)

    !> TNeighbourList instance.
    type(TNeighbourList), intent(inout) :: this

    if (this%useMpiWindows_) then
  #:if WITH_MPI
      if (associated(this%iNeighbourMem_)) call this%iNeighbourWin_%free()
      if (associated(this%neighDist2Mem_)) call this%neighDist2Win_%free()
  #:endif
    else
      if (associated(this%iNeighbourMem_)) deallocate(this%iNeighbourMem_)
      if (associated(this%neighDist2Mem_)) deallocate(this%neighDist2Mem_)
    end if

  end subroutine TNeighbourList_final


  !> Finalizes the symmetric neighbour-list instance.
  !!
  !! Workaround: Intel oneAPI 2021/22
  !! Without explicit deallocation, the oneAPI versions listed above do not correctly finalize the
  !! MPI windows.
  subroutine TSymNeighbourList_final(this)

    !> TSymNeighbourList instance
    type(TSymNeighbourList), intent(inout) :: this

    if (allocated(this%neighbourList)) deallocate(this%neighbourList)

  end subroutine TSymNeighbourList_final


  !> Calculates the translation vectors for cells, which could contain atoms interacting with any of
  !> the atoms in the central cell.
  !> This subroutine uses a simple guess to get the necessary translation vectors. This results in a
  !> set of vectors which could for very asymmetric cells a large amount bigger than the real
  !> necessary one.
  subroutine getCellTranslations(cellVec, rCellVec, latVec, recVec2p, cutoff)

    !> Returns cell translation vectors in relative coordinates.
    real(dp), allocatable, intent(out) :: cellVec(:,:)

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
  !> between two arbitrary cells.
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

    !> Whether to include time reversal symmetry when generating k-points
    logical,  intent(in), optional :: reduceByInversion

    !> Whether to exclude the (0,0,0) point
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


  !> Updates the neighbour list and the species arrays.
  subroutine updateNeighbourListAndSpecies(env, coord, species, img2CentCell, iCellVec, neigh,&
      & nAllAtom, coord0, species0, cutoff, rCellVec, errStatus, symmetric, helicalBoundConds)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

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

    !> Operation status, if an error needs to be returned
    type(TStatus), intent(out) :: errStatus

    !> Whether the neighbour list should be symmetric or not (default)
    logical, intent(in), optional :: symmetric

    !> Helical translation and angle, if necessary, along z axis
    real(dp), intent(in), optional :: helicalBoundConds(:,:)

    if (neigh%setExternally) then
      return
    end if

    call updateNeighbourList(coord, img2CentCell, iCellVec, neigh, nAllAtom, coord0, cutoff,&
        & rCellVec, errStatus, env, symmetric, helicalBoundConds)
    @:PROPAGATE_ERROR(errStatus)

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
      & rCellVec, errStatus, env, symmetric, helicalBoundConds)

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

    !> Operation status, if an error needs to be returned
    type(TStatus), intent(out) :: errStatus

    !> Environment settings
    type(TEnvironment), intent(in), optional :: env

    !> Optional, whether the map should be symmetric (dftb default = .false.)
    logical, intent(in), optional :: symmetric

    !> Helical translation and angle, if necessary, along z axis
    real(dp), intent(in), optional :: helicalBoundConds(:,:)

    ! Nr. of atoms in the system
    integer :: nAtom

    ! Max. nr. of atom without reallocation
    integer :: mAtom

    ! Max. nr. of neighbours without reallocation
    integer :: maxNeighbour, maxNeighbourLocal

    ! Nr. of cell translation vectors
    integer :: nCellVec

    ! Square of the diatomic interaction cutoffs
    real(dp) :: cutoff2

    real(dp) :: dist2
    real(dp) :: rCell(3), rr(3)
    integer :: ii, iAtom1, oldIAtom1, iAtom2, startAtom, endAtom
    integer :: nn1, iAtom2End, pairError(2)
    logical :: symm, isParallel, isSetupError
    real(dp), allocatable :: neighDist2(:,:)
    integer, allocatable :: indx(:), iNeighbour(:,:)
    character(len=100) :: strError

  #:if WITH_MPI
    logical :: isParallelSetupError
    integer, allocatable :: buffer(:,:)
  #:endif

    nAtom = size(neigh%nNeighbour, dim=1)
    mAtom = size(coord, dim=2)
    nCellVec = size(rCellVec, dim=2)

    @:ASSERT(nAtom <= mAtom)
    @:ASSERT(allocated(coord))
    @:ASSERT(size(coord, dim=1) == 3)
    @:ASSERT(allocated(img2CentCell))
    @:ASSERT(size(img2CentCell) == mAtom)
    @:ASSERT(allocated(iCellVec))
    @:ASSERT(size(iCellVec) == mAtom)
    @:ASSERT((size(coord0, dim=1) == 3) .and. size(coord0, dim=2) >= nAtom)
    @:ASSERT(cutoff >= 0.0_dp)

    symm = .false.
    if (present(symmetric)) then
      symm = symmetric
    end if
    neigh%cutoff = cutoff
    cutoff2 = cutoff**2
    nAllAtom = 0

    isParallel = .false.
  #:if WITH_MPI
    if (present(env)) then
      call distributeAtoms(env%mpi%nodeComm%rank, env%mpi%nodeComm%size, nAtom, startAtom, endAtom,&
          & isParallelSetupError)
      isParallel = .not. isParallelSetupError
    end if

    neigh%useMpiWindows_ = isParallel
  #:endif

    if (.not. isParallel) then
      startAtom = 1
      endAtom = nAtom
    end if

    maxNeighbour = max(nCellVec * nAtom / 4, 1)
    allocate(iNeighbour(1:maxNeighbour, startAtom:endAtom))
    allocate(neighDist2(1:maxNeighbour, startAtom:endAtom))

    ! Clean arrays.
    ! (Every atom is the 0th neighbour of itself with zero distance square.)
    neigh%nNeighbour(:) = 0
    iNeighbour(:,:) = 0
    neighDist2(:,:) = 0.0_dp

    ! No erroneous atom pair(s) at moment
    pairError(:) = 0

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
            call rotate3(rr, 2.0_dp * pi * rCellVec(2, ii) / helicalBoundConds(3, 1), zAxis)
          end if
          ! helical operation, note nint() not floor() as roundoff can cause problems for floor
          ! here.
          call rotate3(rr, helicalBoundConds(2,1) * nint(rCellVec(1, ii) / helicalBoundConds(1,1)),&
              & zAxis)
        end if
        lpIAtom2: do iAtom2 = 1, iAtom2End
          !  If distance greater than cutoff -> skip
          dist2 = sum((coord0(:, iAtom2) - rr)**2)
          if (dist2 > cutoff2) then
            cycle lpIAtom2
          end if
          ! New interacting atom -> append
          ! We need that before checking for interaction with dummy atom or
          ! with itself to make sure that atoms in the central cell are
          ! appended exactly in the same order as found in the coord0 array.
          if (iAtom1 /= oldIAtom1) then
            nAllAtom = nAllAtom + 1
            if (nAllAtom > mAtom) then
              mAtom = incrmntOfArray(mAtom)
              call reallocateArrays1(img2CentCell, iCellVec, coord, mAtom)
            end if
            coord(:, nAllAtom) = rr
            img2CentCell(nAllAtom) = iAtom1
            iCellVec(nAllAtom) = ii
            oldIAtom1 = iAtom1
          end if

          if (startAtom <= iAtom2 .and. iAtom2 <= endAtom) then
            ! Check if atoms are not too close to each other
            if (dist2 < minNeighDist2) then
              if (ii == 1 .and. iAtom1 == iAtom2) then
                ! We calculated the distance between the same atom in the unit cell
                cycle lpIAtom2
              else
                ! proximity problem
                pairError(:) = [iAtom1, iAtom2]
                exit lpCellVec
              end if
            end if

            neigh%nNeighbour(iAtom2) = neigh%nNeighbour(iAtom2) + 1
            if (neigh%nNeighbour(iAtom2) > maxNeighbour) then
              maxNeighbour = incrmntOfArray(maxNeighbour)
              call reallocateArrays2(iNeighbour, neighDist2, maxNeighbour)
            end if
            iNeighbour(neigh%nNeighbour(iAtom2), iAtom2) = nAllAtom
            neighDist2(neigh%nNeighbour(iAtom2), iAtom2) = dist2
          end if
        end do lpIAtom2
      end do lpIAtom1
    end do lpCellVec

    isSetupError = .false.
    if (pairError(1) /= 0) then
      isSetupError = .true.
    end if
  #:if WITH_MPI
    if (neigh%useMpiWindows_) then
      ! find if any of the processes in the node comm are in error state
      call mpifx_allreduceip(env%mpi%nodeComm, isSetupError, MPI_LOR)
      if (isSetupError) then
        ! synchronize over node and choose highest atom in erroneous pair
        allocate(buffer(2,env%mpi%nodeComm%size))
        call mpifx_allgather(env%mpi%nodeComm, pairError, buffer)
        pairError(:) = reshape(buffer(:,maxloc(buffer(1,:))),[2])
      end if
    end if
  #:endif
    if (isSetupError) then
      write(strError, "(A,I0,A,I0,A)") "Atoms ", pairError(1), " and ", pairError(2),&
          & " too close together"
      @:RAISE_ERROR(errStatus, -1, strError)
    end if

    call reallocateArrays1(img2CentCell, iCellVec, coord, nAllAtom)

    #:if WITH_MPI
    if (neigh%useMpiWindows_) then
      call mpifx_allreduceip(env%mpi%nodeComm, neigh%nNeighbour, MPI_MAX)
    end if
    #:endif

    maxNeighbour = maxval(neigh%nNeighbour(1:nAtom))
    maxNeighbourLocal = min(ubound(iNeighbour, dim=1), maxNeighbour)

    ! Sort neighbours for all atom by distance
    allocate(indx(maxNeighbourLocal))
    lpStoreAtoms: do iAtom1 = startAtom, endAtom
      nn1 = neigh%nNeighbour(iAtom1)
      call index_heap_sort(indx(1:nn1), neighDist2(1:nn1, iAtom1), tolSameDist2)
      iNeighbour(1:nn1, iAtom1) = iNeighbour(indx(1:nn1), iAtom1)
      neighDist2(1:nn1, iAtom1) = neighDist2(indx(1:nn1), iAtom1)

      if (nn1 < maxNeighbourLocal) then
        iNeighbour(nn1+1:maxNeighbourLocal, iAtom1) = 0
        neighDist2(nn1+1:maxNeighbourLocal, iAtom1) = 0.0_dp
      end if

    end do lpStoreAtoms

    call allocateNeighbourArrays(neigh, maxNeighbour, nAtom, env)
    call fillNeighbourArrays(neigh, iNeighbour, neighDist2, startAtom, endAtom, maxNeighbour,&
        & nAtom)

  end subroutine updateNeighbourList


  !> Explicitly set the neighbour list
  subroutine setNeighbourList(neigh, env, nNeighbour, iNeighbour, neighDist, cutOff, coord0,&
      & species0, coordNeighs, neigh2CentCell, rCellVec, nAllAtom, img2CentCell, iCellVec, coord,&
      & species)

    !> Neighbourlist instance
    type(TNeighbourList), intent(inout) :: neigh

    !> Environment
    type(TEnvironment), intent(inout) :: env

    !> Number of neighbours of an atom in the central cell
    integer, intent(in) :: nNeighbour(:)

    !> References to the neighbour atoms for an atom in the central cell
    integer, intent(in) :: iNeighbour(:,:)

    !> Distances to the neighbour atoms for an atom in the central cell
    real(dp), intent(in) :: neighDist(:,:)

    !> Cutoff distance used for this neighbour list
    real(dp), intent(in) :: cutOff

    !> Coordinates of the central cell
    real(dp), intent(in) :: coord0(:,:)

    !> Species of the atoms in the central cell
    integer, intent(in) :: species0(:)

    !> Coordinates of all neighbours
    real(dp), intent(in) :: coordNeighs(:,:)

    !> Mapping between neighbour reference and atom index in the central cell
    integer, intent(in) :: neigh2CentCell(:)

    real(dp), contiguous, target, intent(in) :: rCellVec(:,:)

    !> Nr. of all atoms in the system
    integer, intent(out) :: nAllAtom

    !> Mapping of the atoms into the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index of unit cell containing atom
    integer, allocatable, intent(inout) :: iCellVec(:)

    !> Coordinates of all atoms (may contain duplicates, access it only via neigh%iNeighbour)
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Species of all atoms (may contain duplicates, access it only via neigh%iNeighbour)
    integer, allocatable, intent(inout) :: species(:)

    type(TAlignedArray) :: dist2Mem
    real(dp), pointer :: rCellVecFlat(:), dist2(:)
    real(dp) :: diff(3)
    integer :: nMaxNeighbours, nCellVec, nAtom
    integer :: iAtom, iCell, iNeigh, iImage, iCellFlat, iAtFirst, iAtLast
    integer, allocatable :: indx(:)
    logical :: copyData

    nMaxNeighbours = maxval(nNeighbour)
    nAtom = size(nNeighbour)

    neigh%setExternally = .true.
  #:if WITH_MPI
    neigh%useMpiWindows_ = .true.
  #:endif

    call allocateNeighbourArrays(neigh, nMaxNeighbours, nAtom, env)

    neigh%nNeighbour(:) = nNeighbour(:)
    neigh%cutoff = cutOff

    nAllAtom = size(neigh2CentCell) + nAtom
    @:ASSERT(size(neigh2CentCell) == size(coordNeighs, dim=2))

    if (size(img2CentCell) /= nAllAtom) then
      call reallocateArrays1(img2CentCell, iCellVec, coord, nAllAtom)
    end if

    !> Prepend data for the atoms in the central cell
    coord(1:3, 1:nAtom) = coord0(1:3, 1:nAtom)
    img2CentCell(1:nAtom) = [(iAtom, iAtom = 1, nAtom)]
    iCellVec(1:nAtom) = 1

    if (nAtom < nAllAtom) then
      coord(1:3,  nAtom + 1:) = coordNeighs(1:3,:)
      img2CentCell(nAtom + 1:) = neigh2CentCell
      iCellVec(nAtom + 1:) = 0

      !> Now set iCellVec: Iterate over all cells, calculate the coordinates the atom would
      !> have there, and determine the cell in which this is very close to the actual coordinates
      nCellVec = size(rCellVec, dim=2)
      call dist2Mem%allocate(nCellVec)
      call dist2Mem%getArray(dist2)

      rCellVecFlat(1 : 3 * nCellVec) => rCellVec

      call distributeRangeInChunks(env, nAtom + 1, nAllAtom, iAtFirst, iAtLast)

      do iAtom = iAtFirst, iAtLast
        iImage = img2CentCell(iAtom)
        diff(:) = coord0(:,iImage) - coord(:,iAtom)
        do iCell = 1, nCellVec
          iCellFlat = 3 * (iCell - 1)
          ! Note that we need an explicit loop unrolling here to assist vectorization.
          ! Using the pointer rCellVecFlat helps the compiler to assume continuous memory access.
          dist2(iCell) = (diff(1) + rCellVecFlat(iCellFlat + 1))**2&
              & + (diff(2) + rCellVecFlat(iCellFlat + 2))**2&
              & + (diff(3) + rCellVecFlat(iCellFlat + 3))**2
        end do
        iCellVec(iAtom) = minloc(dist2, dim=1)
      end do

      call assembleChunks(env, iCellVec(nAtom + 1:))
    end if

    if (size(species) /= nAllAtom) then
      deallocate(species)
      allocate(species(nAllAtom))
    end if
    species(1:nAllAtom) = species0(img2CentCell(1:nAllAtom))

    copyData = .true.
  #:if WITH_MPI
    if (neigh%useMpiWindows_) then
      call neigh%iNeighbourWin_%lock_all()
      call neigh%neighDist2Win_%lock_all()
      copyData = env%mpi%nodeComm%lead
    end if
  #:endif

    !> This is done only for task 0 on the node due to MPI shared memory: Copy to the actual
    !> neighbour arrays.
    if (copyData) then
      @:ASSERT(nMaxNeighbours <= size(iNeighbour, dim=1))
      @:ASSERT(nMaxNeighbours <= size(neighDist, dim=1))
      @:ASSERT(size(iNeighbour, dim=2) == nAtom)
      @:ASSERT(size(neighDist, dim=2) == nAtom)

      allocate(indx(maxval(nNeighbour)))
      do iAtom = 1, nAtom
        call index_heap_sort(indx(1:nNeighbour(iAtom)), neighDist(1:nNeighbour(iAtom), iAtom),&
            & tolSameDist)
        !$OMP SIMD
        do iNeigh = 1, nNeighbour(iAtom)
          neigh%iNeighbour(iNeigh, iAtom) = iNeighbour(indx(iNeigh), iAtom) + nAtom
        end do
        !$OMP SIMD
        do iNeigh = 1, nNeighbour(iAtom)
          neigh%neighDist2(iNeigh, iAtom) = neighDist(indx(iNeigh), iAtom)**2
        end do

        if (nNeighbour(iAtom) < nMaxNeighbours) then
          !$OMP SIMD
          do iNeigh = nNeighbour(iAtom) + 1, nMaxNeighbours
            neigh%iNeighbour(iNeigh, iAtom) = 0
          end do
          !$OMP SIMD
          do iNeigh = nNeighbour(iAtom) + 1, nMaxNeighbours
            neigh%neighDist2(iNeigh, iAtom) = 0.0_dp
          end do
        end if

        neigh%iNeighbour(0, iAtom) = iAtom
        neigh%neighDist2(0, iAtom) = 0.0_dp
      end do
    end if

  #:if WITH_MPI
    if (neigh%useMpiWindows_) then
      call neigh%iNeighbourWin_%sync()
      call neigh%neighDist2Win_%sync()

      call neigh%iNeighbourWin_%unlock_all()
      call neigh%neighDist2Win_%unlock_all()
    end if
  #:endif

  end subroutine setNeighbourList


  !> Allocate arrays for type 'neigh'
  subroutine allocateNeighbourArrays(neigh, maxNeighbour, nAtom, env)

    !> Contains all neighbour information
    type(TNeighbourList), intent(inout) :: neigh

    !> Maximum number of neighbours an atom can have
    integer, intent(in) :: maxNeighbour

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Environment settings
    type(TEnvironment), intent(in), optional :: env

    integer :: dataLength
    #:if WITH_MPI
      integer(MPIFX_SIZE_T) :: longDataLength
    #:endif

    dataLength = (maxNeighbour + 1) * nAtom

    if (neigh%useMpiWindows_) then
    #:if WITH_MPI
      longDataLength = int(dataLength, kind=MPIFX_SIZE_T)
      if (associated(neigh%iNeighbourMem_)) call neigh%iNeighbourWin_%free()
      call neigh%iNeighbourWin_%allocate_shared(env%mpi%nodeComm, longDataLength,&
          & neigh%iNeighbourMem_)
      if (associated(neigh%neighDist2Mem_)) call neigh%neighDist2Win_%free()
      call neigh%neighDist2Win_%allocate_shared(env%mpi%nodeComm, longDataLength,&
          & neigh%neighDist2Mem_)
    #:endif
    else
      if (associated(neigh%iNeighbourMem_)) deallocate(neigh%iNeighbourMem_)
      allocate(neigh%iNeighbourMem_(dataLength))
      if (associated(neigh%neighDist2Mem_)) deallocate(neigh%neighDist2Mem_)
      allocate(neigh%neighDist2Mem_(dataLength))
    end if

    neigh%iNeighbour(0:maxNeighbour, 1:nAtom) => neigh%iNeighbourMem_(1:dataLength)
    neigh%neighDist2(0:maxNeighbour, 1:nAtom) => neigh%neighDist2Mem_(1:dataLength)

  end subroutine allocateNeighbourArrays


  !> Collect all neighbour data and copy to neighbour arrays
  subroutine fillNeighbourArrays(neigh, iNeighbour, neighDist2, startAtom, endAtom, maxNeighbour,&
      & nAtom)

    !> Contains all neighbour information
    type(TNeighbourList), intent(inout) :: neigh

    !> Array of neighbour indices
    integer, allocatable, intent(in) :: iNeighbour(:,:)

    !> Array if neighbour distances
    real(dp), allocatable, intent(in) :: neighDist2(:,:)

    !> First and last atomic indices for the current MPI process
    integer, intent(in) :: startAtom, endAtom

    !> Maximum number of neighbours an atom can have
    integer, intent(in) :: maxNeighbour

    !> Number of atoms
    integer, intent(in) :: nAtom

    integer :: ii

  #:if WITH_MPI
    integer :: maxNeighbourLocal
  #:endif

    if (neigh%useMpiWindows_) then
    #:if WITH_MPI
      maxNeighbourLocal = min(ubound(iNeighbour, dim=1), maxNeighbour)

      call neigh%iNeighbourWin_%lock_all()
      call neigh%neighDist2Win_%lock_all()

      neigh%iNeighbour(1:maxNeighbourLocal, startAtom:endAtom) =&
          & iNeighbour(1:maxNeighbourLocal, startAtom:endAtom)
      neigh%neighDist2(1:maxNeighbourLocal, startAtom:endAtom) =&
          & neighDist2(1:maxNeighbourLocal, startAtom:endAtom)

      if (maxNeighbourLocal < maxNeighbour) then
        neigh%iNeighbour(maxNeighbourLocal + 1 : maxNeighbour, startAtom:endAtom) = 0
        neigh%neighDist2(maxNeighbourLocal + 1 : maxNeighbour, startAtom:endAtom) = 0.0_dp
      end if

      do ii = 1, nAtom
        neigh%iNeighbour(0, ii) = ii
        neigh%neighDist2(0, ii) = 0.0_dp
      end do

      call neigh%iNeighbourWin_%sync()
      call neigh%neighDist2Win_%sync()

      call neigh%iNeighbourWin_%unlock_all()
      call neigh%neighDist2Win_%unlock_all()
    #:endif
    else
      neigh%iNeighbour(1:,:) = iNeighbour(1:maxNeighbour,:)
      neigh%neighDist2(1:,:) = neighDist2(1:maxNeighbour,:)

      do ii = 1, nAtom
        neigh%iNeighbour(0, ii) = ii
        neigh%neighDist2(0, ii) = 0.0_dp
      end do

    end if

  end subroutine fillNeighbourArrays


  !> Returns the nr. of neighbours for a given cutoff for all atoms.
  subroutine getNrOfNeighboursForAll(nNeighbourSK, neigh, cutoff)

    !> Contains the nr. of neighbours for each atom on exit.
    integer, intent(out) :: nNeighbourSK(:)

    !> Initialized neighbourlist
    type(TNeighbourList), intent(in) :: neigh

    !> Maximal neighbour distance to consider.
    real(dp), intent(in) :: cutoff

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

    ! Issue warning, if cutoff is bigger than used for the neighbourlist.
    if (cutoff > neigh%cutoff) then
99010 format ('Cutoff (', E16.6, ') greater than last cutoff ', '(', E13.6,&
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

    !> Array mapping images of atoms to originals in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index of unit cell containing atom
    integer, allocatable, intent(inout) :: iCellVec(:)

    !> Coordinates of all atoms (actual and image)
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Maximum number of new atoms
    integer, intent(in) :: mNewAtom

    integer :: mAtom
    integer, allocatable :: tmpIntR1(:)
    real(dp), allocatable :: tmpRealR2(:,:)

    mAtom = size(img2CentCell)

    @:ASSERT(size(iCellVec) == mAtom)
    @:ASSERT(all(shape(coord) == (/ 3, mAtom /)))
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
  subroutine reallocateArrays2(iNeighbour, neighDist2, mNewNeighbour)

    !> List of neighbours
    integer, allocatable, intent(inout) :: iNeighbour(:, :)

    !> Square of distances between atoms
    real(dp), allocatable, intent(inout) :: neighDist2(:,:)

    !> Maximum number of new atoms
    integer, intent(in) :: mNewNeighbour

    integer :: mNeighbour, startAtom, endAtom
    integer, allocatable :: tmpIntR2(:,:)
    real(dp), allocatable :: tmpRealR2(:,:)

    mNeighbour = size(iNeighbour, dim=1)
    startAtom = lbound(iNeighbour, dim=2)
    endAtom = ubound(iNeighbour, dim=2)

    @:ASSERT(mNewNeighbour > 0)
    @:ASSERT(all(shape(neighDist2) == shape(iNeighbour)))

    call move_alloc(iNeighbour, tmpIntR2)
    allocate(iNeighbour(mNewNeighbour, startAtom:endAtom))
    if (mNewNeighbour > mNeighbour) then
      iNeighbour(mNeighbour+1:,:) = 0
      iNeighbour(:mNeighbour,:) = tmpIntR2(:,:)
    else
      iNeighbour(:,:) = tmpIntR2(:mNewNeighbour,:)
    end if

    call move_alloc(neighDist2, tmpRealR2)
    allocate(neighDist2(mNewNeighbour, startAtom:endAtom))
    if (mNewNeighbour > mNeighbour) then
      neighDist2(mNeighbour+1:,:) = 0.0_dp
      neighDist2(:mNeighbour,:) = tmpRealR2(:,:)
    else
      neighDist2(:,:) = tmpRealR2(:mNewNeighbour,:)
    end if

  end subroutine reallocateArrays2


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


  !> Convert fractional coordinates to cartesian
  subroutine frac2cart(cartCoords,latvecs)

    !> Fractional coordinates in unit cell on entry, cartesian on exit
    real(dp), intent(inout) :: cartCoords(:,:)

    !> Periodic lattice vectors
    real(dp), intent(in) :: latvecs(3,3)

    @:ASSERT(size(cartCoords,dim=1) == 3)

    cartCoords = matmul(latvecs,cartCoords)

  end subroutine frac2cart


  !> Cartesian to fractional coordinates in periodic geometry
  subroutine cart2frac(cartCoords,latvecs)

    !> Cartesian coordinates on entry, fractional on exit
    real(dp), intent(inout) :: cartCoords(:,:)

    !> Periodic lattice vectors
    real(dp), intent(in) :: latvecs(3,3)

    real(dp) :: invLatVecs(3,3)

    @:ASSERT(size(cartCoords,dim=1) == 3)

    call invert33(invLatVecs, latvecs)
    cartCoords = matmul(invLatvecs, cartCoords)

  end subroutine cart2frac


  !> Computes a domain decomposition of n atoms
  subroutine distributeAtoms(mpiRank, mpiSize, nAtom, startAtom, endAtom, error)

    !> Current MPI rank
    integer, intent(in) :: mpiRank

    !> Number of MPI processes
    integer, intent(in) :: mpiSize

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Start and end atomic indices for the current MPI process
    integer, intent(out) :: startAtom, endAtom

    !> Whether an error occurred during domain decomposition
    logical, intent(out) :: error

    if (nAtom < mpiSize) then
      call warning("Cannot parallelize atomic distance computation: &
          &The number of MPI ranks must not be greater than the number of atoms.")
      error = .true.
      return
    end if

    call getChunkRanges(mpiSize, mpiRank, 1, nAtom, startAtom, endAtom)
    error = .false.

  end subroutine distributeAtoms

end module dftbp_dftb_periodic
