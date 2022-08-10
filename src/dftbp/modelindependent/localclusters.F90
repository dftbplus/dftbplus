!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains routines for producing local environments around atoms and convert between neighbour map
!> indexed information and local geometries
module dftbp_modelindependent_localclusters
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : Bohr__AA
  use dftbp_io_formatout, only : writeXYZFormat
  use dftbp_common_memman, only : incrmntOfArray
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_math_simplegeometry, only : dist2Segment, radialSort, uniqueCoords
  use iso_c_binding, only : c_int, c_double
  implicit none

#:include 'common.fypp'

  private
  public :: TLocalClusters, TLocalClusters_init

  type TLocalClusters

    private

    !> Atoms in central cell of the system
    integer :: nAtomCentCell

    !> Cutoff for interacting atoms around either sites or bonds in the neighbour map
    real(dp) :: envCutoff

    !> Cutoff for longest bond with interactions
    real(dp) :: bondCutoff

    !> Square of atomic surroundings to give rounded cylinder of radius of envCutoff around bonds
    real(dp) :: maxCutoff2

    !> Number of atoms surrounding each central atom choice
    integer, allocatable :: nAtomicHalo(:)

    !> Coordinates for neighbouring atomic sites around each central atom (which is placed at the
    !> origin)
    real(dp), allocatable :: haloCoords(:,:,:)

    !> Central cell number of atoms in the halo
    integer, allocatable :: haloAtomNos(:,:)

    ! Exposed cluster variables

    !> Start of atom clusters in geometry atomClusters
    integer(c_int), allocatable, public :: iStartAtCluster(:)

    !> Original atom numbers of atoms in clusters acording to central cell of structure
    integer(c_int), allocatable, public :: iAtCluster(:)

    !> Clusters around atoms, compressed over last index (3, sum(nAtsAroundAt(:)))
    real(c_double), allocatable, public :: atomClusters(:,:)

    !> Number of bonds that have surrounding clusters
    integer(c_int), public :: nBondClusters

    !> Start of specific bond clusters in the bondClusters array of coordinates
    integer(c_int), allocatable, public :: iStartBondCluster(:)

    !> Original atom numbers in central cell of structure
    integer(c_int), allocatable, public :: iBondCluster(:)

    !> Bonds around atoms, compressed over last index (3, sum(nAtsAroundBnd(:)))
    real(c_double), allocatable, public :: bondClusters(:,:)

  contains

    procedure :: updateGeometry
    procedure :: generateClusters

  end type TLocalClusters

contains

  !> Initialise local clusters of atoms around atomic sites
  subroutine TLocalClusters_init(this, nAtCentCell, envCutoff, bondCutoff)

    !> Instance
    type(TLocalClusters), intent(out) :: this

    !> Number of atoms in the central region cell
    integer, intent(in) :: nAtCentCell

    !> Cutoff distance for atoms surrounding bonds
    real(dp), intent(in) :: envCutoff

    !> Longest bond distance present
    real(dp), intent(in) :: bondCutoff

    this%envCutoff = envCutoff
    this%bondCutoff = bondCutoff
    ! Distance for a rounded cylinder centered on a bond of length bondCutOff which is contained
    ! within envCutoff of points along the bond
    this%maxCutoff2 = envCutoff**2 + (bondCutoff**2)/4.0_dp

    this%natomCentCell = nAtCentCell

  end subroutine TLocalClusters_init


  !> Updates local cluster geometries around atoms
  subroutine updateGeometry(this, coords, neighbourList, img2CentCell, species)

    !> Instance
    class(TLocalClusters), intent(inout) :: this

    !> Coordinates for the system
    real(dp), intent(in) :: coords(:,:)

    !> ADT for neighbour parameters
    type(TNeighbourList), intent(in) :: neighbourList

    !> Folds images into the central atom cell
    integer, intent(in) :: img2CentCell(:)

    !> Chemical species of atoms
    integer, intent(in) :: species(:)

    integer :: iAt1, iNeigh, iAt2, iAt2f, iSp1, iSp2
    integer :: initialAts
    integer, allocatable :: indx(:)

    initialAts = 5
    if (allocated(this%nAtomicHalo)) then
      initialAts = size(this%haloCoords, dim=2)
      deallocate(this%nAtomicHalo)
      deallocate(this%haloCoords)
      deallocate(this%haloAtomNos)
    end if
    allocate(this%nAtomicHalo(this%nAtomCentCell), source=0)
    allocate(this%haloCoords(3, initialAts, this%nAtomCentCell), source=0.0_dp)
    allocate(this%haloAtomNos(initialAts, this%nAtomCentCell), source=0)

    @:ASSERT(this%maxCutoff2 <= neighbourList%cutoff**2)

    do iAt1 = 1, this%nAtomCentCell
      iSp1 = species(iAt1)
      lpNeigh: do iNeigh = 0, neighbourList%nNeighbour(iAt1)
        iAt2 = neighbourList%iNeighbour(iNeigh, iAt1)
        if (sum((coords(:,iAt2)-coords(:,iAt1))**2) > this%maxCutoff2) then
          exit lpNeigh
        end if
        iAt2f = img2centcell(iAt2)
        iSp2 = species(iAt2f)
        this%nAtomicHalo(iAt1) = this%nAtomicHalo(iAt1) + 1
        if (any([this%nAtomicHalo(iAt1), this%nAtomicHalo(iAt2f)] >=&
            & size(this%haloAtomNos, dim=1))) then
          call resizeAtHalos(this%haloCoords, this%haloAtomNos)
        end if
        this%haloAtomNos(this%nAtomicHalo(iAt1), iAt1) = iAt2f
        this%halocoords(:, this%nAtomicHalo(iAt1), iAt1) = coords(:,iAt2) - coords(:,iAt1)
        if (iNeigh > 0) then
          this%nAtomicHalo(iAt2f) = this%nAtomicHalo(iAt2f) + 1
          this%haloAtomNos(this%nAtomicHalo(iAt2f), iAt2f) = iAt1
          this%halocoords(:, this%nAtomicHalo(iAt2f), iAt2f) = coords(:,iAt1) - coords(:,iAt2)
        end if
      end do lpNeigh
    end do

    ! sort by distance from central atom in each cluster
    do iAt1 = 1, this%nAtomCentCell
      allocate(indx(this%nAtomicHalo(iAt1)))
      call radialSort(this%halocoords(:, :this%nAtomicHalo(iAt1), iAt1), indx,&
          & this%nAtomicHalo(iAt1))
      this%haloCoords(:, :this%nAtomicHalo(iAt1), iAt1) =&
          & this%haloCoords(:, indx(:this%nAtomicHalo(iAt1)), iAt1)
      this%haloAtomNos(:this%nAtomicHalo(iAt1), iAt1) =&
          & this%haloAtomNos(indx(:this%nAtomicHalo(iAt1)), iAt1)
      deallocate(indx)
    end do

  end subroutine updateGeometry


  !> Generate clusters around atoms
  subroutine generateClusters(this, coords, nNeighbourSK, neighbourList, img2centcell, species,&
      & speciesName)

    !> Instance
    class(TLocalClusters), intent(inout) :: this

    !> Coordinates for the system
    real(dp), intent(in) :: coords(:,:)

    !> Number of neighbours of each real atom
    integer, intent(out) :: nNeighbourSK(:)

    !> ADT for neighbour parameters
    type(TNeighbourList), intent(in) :: neighbourList

    !> Folds images into the central atom cell
    integer, intent(in) :: img2CentCell(:)

    !> Chemical species
    integer, intent(in) :: species(:)

    !> label for each atomic chemical species
    character(*), intent(in), optional :: speciesName(:)

    call generateAtomClusters(this, species, speciesName)

    call generateBondClusters(this, coords, nNeighbourSK, neighbourList, img2centcell, species,&
        & speciesName)

  end subroutine generateClusters


  !> Generate packed storage of atomic clusters
  subroutine generateAtomClusters(this, species, speciesName)

    !> Instance
    type(TLocalClusters), intent(inout) :: this

    !> Chemical species
    integer, intent(in) :: species(:)

    !> label for each atomic chemical species
    character(*), intent(in), optional :: speciesName(:)

    integer :: iAt1, iAt2, nAt
    logical :: areClustersWritten

    areClustersWritten = present(speciesName)

    if (allocated(this%iStartAtCluster)) then
      if (size(this%iStartAtCluster) <= this%nAtomCentCell+1) then
        deallocate(this%iStartAtCluster)
      else
        this%iStartAtCluster(:) = 0
      end if
    end if
    if (.not.allocated(this%iStartAtCluster)) then
      allocate(this%iStartAtCluster(this%nAtomCentCell+1), source=0)
    end if
    this%iStartAtCluster(1) = 1
    do iAt1 = 1, this%nAtomCentCell
      ! count atoms in this cluster
      nAt = 0
      lpHalo1: do iAt2 = 1, this%nAtomicHalo(iAt1)
        ! Should use bisection if a large number of atoms is present in the halo, as already sorted
        if (sum(this%haloCoords(:, iAt2, iAt1)**2) <= this%envCutoff**2) then
          nAt = nAt + 1
        else
          exit lpHalo1
        end if
      end do lpHalo1
      this%iStartAtCluster(iAt1+1) = this%iStartAtCluster(iAt1) + nAt
    end do
    if (allocated(this%atomClusters)) then
      if (size(this%atomClusters, dim=2) < this%iStartAtCluster(this%nAtomCentCell+1)-1) then
        deallocate(this%atomClusters)
        deallocate(this%iAtCluster)
      end if
    end if
    if (.not.allocated(this%atomClusters)) then
      allocate(this%atomClusters(3,this%iStartAtCluster(this%nAtomCentCell+1)-1), source=0.0_dp)
      allocate(this%iAtCluster(this%iStartAtCluster(this%nAtomCentCell+1)-1), source=0)
    end if

    do iAt1 = 1, this%nAtomCentCell
      iAt2 = this%iStartAtCluster(iAt1+1) - this%iStartAtCluster(iAt1)
      this%atomClusters(:, this%iStartAtCluster(iAt1):this%iStartAtCluster(iAt1+1)-1) =&
          & this%haloCoords(:, :iAt2, iAt1)
      this%iAtCluster(this%iStartAtCluster(iAt1):this%iStartAtCluster(iAt1+1)-1)&
          & = this%haloAtomNos(:iAt2, iAt1)
    end do

    if (areClustersWritten) then
      do iAt1 = 1, this%nAtomCentCell
        iAt2 = this%iStartAtCluster(iAt1+1)-this%iStartAtCluster(iAt1)
        call writeXYZFormat("atomcluster.xyz", this%haloCoords(:, :iAt2, iAt1),&
            & species(this%haloAtomNos(:iAt2, iAt1)), speciesName, append = iAt1 > 1)
      end do
    end if

  end subroutine generateAtomClusters


  !> Generate clusters around bonds
  subroutine generateBondClusters(this, coords, nNeighbourSK, neighbourList, img2centcell, species,&
      & speciesName)

    !> Instance
    type(TLocalClusters), intent(inout) :: this

    !> Coordinates for the system
    real(dp), intent(in) :: coords(:,:)

    !> Number of neighbours of each real atom
    integer, intent(out) :: nNeighbourSK(:)

    !> ADT for neighbour parameters
    type(TNeighbourList), intent(in) :: neighbourList

    !> Folds images into the central atom cell
    integer, intent(in) :: img2CentCell(:)

    !> Chemical species
    integer, intent(in) :: species(:)

    !> label for each atomic chemical species
    character(*), intent(in), optional :: speciesName(:)

    integer :: iAt1, iNeigh, iAt2, iAt2f, ii, nAtClust, kk, nAt, iClust
    integer :: iBond, iStart, iEnd
    real(dp) :: p1(3), p2(3), shift(3)
    real(dp), allocatable :: xUnique(:,:)
    integer, allocatable :: iUnique(:), iAtInHaloGlobNumber(:)
    logical :: areClustersWritten

    areClustersWritten = present(speciesName)

    if (allocated(this%iStartBondCluster)) then
      deallocate(this%iStartBondCluster)
    end if
    if (.not.allocated(this%iStartBondCluster)) then
      ! count number of bonds present, hence number of bond clusters
      allocate(this%iStartBondCluster(sum(nNeighbourSK)+1), source=0)
    end if

    if (.not. allocated(this%bondClusters)) then
      ! smallest possible size, if environment radius is 0
      allocate(this%bondClusters(3, 2*sum(nNeighbourSK)))
      allocate(this%iBondCluster(2*sum(nNeighbourSK)))
    end if
    this%bondClusters(:,:) = 0.0_dp
    this%iBondCluster(:) = 0

    this%iStartBondCluster(1) = 1
    iClust = 1
    ! Generate clusters
    do iAt1 = 1, this%nAtomCentCell
      p1(:) = coords(:, iAt1)
      lpNeigh: do iNeigh = 1, neighbourList%nNeighbour(iAt1)
        if (neighbourList%neighDist2(iNeigh, iAt1) > this%bondCutoff**2) then
          exit lpNeigh
        end if
        iAt2 = neighbourList%iNeighbour(iNeigh, iAt1)
        iAt2f = img2centcell(iAt2)
        p2(:) = coords(:, iAt2)
        shift(:) = (p1 + p2)/2.0_dp

        ! store the lower index atom first
        p1 = coords(:,iAt1)-shift
        p2 = coords(:,iAt2)-shift

        ! count atoms in the cluster around this bond, before de-duplicating
        nAtClust = 2
        do ii = 2, this%nAtomicHalo(iAt1)
          if (dist2Segment(p1, p2, this%halocoords(:, ii, iAt1)+p1) <= this%envCutoff**2) then
            nAtClust = nAtClust + 1
          end if
        end do
        do ii = 2, this%nAtomicHalo(iAt2f)
          if (dist2Segment(p1, p2, this%halocoords(:, ii, iAt2f)+p2) <= this%envCutoff**2) then
            nAtClust = nAtClust + 1
          end if
        end do

        ! set mid-bond at the origin
        shift = (p1 + p2)/2.0_dp

        allocate(xUnique(3,nAtClust))
        allocate(iAtInHaloGlobNumber(nAtClust))

        kk = 1
        xUnique(:,kk) = (p1 - shift)
        iAtInHaloGlobNumber(kk) = this%haloAtomNos(1, iAt1)
        kk = kk + 1
        xUnique(:,kk) = (p2 - shift)
        iAtInHaloGlobNumber(kk) = this%haloAtomNos(1, iAt2f)
        kk = kk + 1
        do ii = 2, this%nAtomicHalo(iAt1)
          if (dist2Segment(p1, p2, this%halocoords(:, ii, iAt1)+p1) <= this%envCutoff**2) then
            xUnique(:,kk) = (this%halocoords(:, ii, iAt1)+p1-shift)
            iAtInHaloGlobNumber(kk) = this%haloAtomNos(ii, iAt1)
            kk = kk + 1
          end if
        end do
        do ii = 2, this%nAtomicHalo(iAt2f)
          if (dist2Segment(p1, p2, this%halocoords(:, ii, iAt2f)+p2) <= this%envCutoff**2) then
            xUnique(:,kk) = (this%halocoords(:, ii, iAt2f)+p2-shift)
            iAtInHaloGlobNumber(kk) = this%haloAtomNos(ii, iAt2f)
            kk = kk + 1
          end if
        end do

        ! de-duplicate coordinates
        call uniqueCoords(xUnique, iUnique)
        ! size of iUnique is now the number of atoms in that bond cluster
        nAt = size(iUnique)
        this%iStartBondCluster(iClust+1) = this%iStartBondCluster(iClust) + nAt

        call resizeBondHalos(this%bondClusters, this%iBondCluster, this%iStartBondCluster(iClust+1))

        this%bondClusters(:, this%iStartBondCluster(iClust):this%iStartBondCluster(iClust+1)-1) =&
            & xUnique(:,iUnique(:nAt))
        this%iBondCluster(this%iStartBondCluster(iClust):this%iStartBondCluster(iClust+1)-1) =&
            & iAtInHaloGlobNumber(iUnique(:nAt))

        deallocate(iUnique)
        deallocate(xUnique)
        deallocate(iAtInHaloGlobNumber)

        iClust = iClust + 1
      end do lpNeigh

    end do
    this%nBondClusters  = iClust -1

    if (areClustersWritten) then
      do iBond = 1, this%nBondClusters
        iStart =  this%iStartBondCluster(iBond)
        iEnd =  this%iStartBondCluster(iBond+1) -1
        call writeXYZFormat("bondcluster.xyz", this%bondClusters(:,iStart:iEnd),&
            & species(this%iBondCluster(iStart:iEnd)), speciesName, append = iBond > 1)
      end do
    end if

  end subroutine generateBondClusters


  !> Resize arrays holding atomic-centred cluster information, if needed
  subroutine resizeAtHalos(haloCoords, haloAtomNos)

    real(dp), allocatable, intent(inout) :: haloCoords(:,:,:)
    integer, allocatable, intent(inout) :: haloAtomNos(:,:)

    integer :: nNew, nOld, nAt

    real(dp), allocatable :: rTmp(:,:,:)
    integer, allocatable :: iTmp(:,:)

    nOld = size(haloAtomNos, dim=1)
    nNew = incrmntOfArray(nOld)
    nAt = size(haloAtomNos, dim=2)

    call move_alloc(haloCoords, rTmp)
    call move_alloc(haloAtomNos, iTmp)
    allocate(haloCoords(3, nNew, nAt), source=0.0_dp)
    allocate(haloAtomNos(nNew, nAt), source=0)
    haloCoords(:, :nOld, :) = rTmp
    haloAtomNos(:nOld, :) = iTmp

  end subroutine resizeAtHalos


  subroutine resizeBondHalos(clusterPos, clusterIndex, neededSize)

    real(dp), allocatable, intent(inout) :: clusterPos(:,:)
    integer, allocatable, intent(inout) :: clusterIndex(:)
    integer, intent(in) :: neededSize

    real(dp), allocatable :: rTmp(:,:)
    integer, allocatable :: iTmp(:)
    integer :: newSize

    integer :: oldSize
    oldSize = size(clusterIndex)
    if (oldSize < neededSize) then
      call move_alloc(clusterPos, rTmp)
      call move_alloc(clusterIndex, iTmp)
      newSize = incrmntOfArray(neededSize)
      allocate(clusterPos(3, newSize), source=0.0_dp)
      allocate(clusterIndex(newSize), source=0)
      clusterPos(:, :oldSize) = rTmp
      clusterIndex(:oldSize) = iTmp
    end if

  end subroutine resizeBondHalos

end module dftbp_modelindependent_localclusters
