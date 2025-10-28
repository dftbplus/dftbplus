!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> A cache for calculating molecule orbitals on a grid.
!! This object is responsible for reading in the eigenvectors from a specified file and passing the
!! appropriate eigenvectors to the molecule orbital calculator.
module waveplot_gridcache
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_file, only : closeFile, openFile, TFileDescr
  use dftbp_common_globalenv, only : stdOut
  use dftbp_io_message, only : error
  use dftbp_wavegrid, only : getValue, TMolecularOrbital
#:if WITH_MPI
  use dftbp_common_schedule, only : getStartAndEndIndex
#:endif
  implicit none

  private
  public :: TGridCache, TGridCache_init

  !> Contains the data for a grid cache.
  type TGridCache

    !> Molecular orbital calculator
    type(TMolecularOrbital), pointer :: molorb

    !> Grid vectors
    real(dp) :: gridVec(3,3)

    !> Origin of the grid
    real(dp) :: origin(3)

    !> Real eigenvectors
    real(dp), allocatable :: eigenvecReal(:,:)

    !> Complex eigenvectors
    complex(dp), allocatable :: eigenvecCmpl(:,:)

    !> Parameters of the levels
    integer, allocatable :: levelIndex(:,:)

    !> Eigenvectors read so far
    integer :: nReadEigVec

    !> Grids processed so far
    integer :: iGrid

    !> Nr. of grids to process
    integer :: nGrid

    !> Position in the cache
    integer :: cachePos

    !> Are we ready?
    logical :: isFinished

    !> File descriptor for eigenvec
    type(TFileDescr) :: fdEigVec

    !> Size of the eigenvectors
    integer :: nOrb

    !> Shape of the grid
    integer :: nPoints(3)

    !> Levels in the eigenvec file
    integer :: nAllLevel

    !> The k-points in the eigenv. file
    integer :: nAllKPoint

    !> Spins in the eigenvec. file
    integer :: nAllSpin

    !> Verbose?
    logical :: beVerbose

    !> Nr. of cached grids
    integer :: nCached

    !> Whether to enable GPU offloading.
    logical :: useGpu

    !> Cache for real grids
    real(dp), allocatable :: gridCacheReal(:,:,:,:)

    !> Cache for complex grids
    complex(dp), allocatable :: gridCacheCmpl(:,:,:,:)

    !> KPoints
    real(dp), allocatable :: kPoints(:,:)

    !> Are eigenvectors real
    logical :: isReal

    !> Initialised?
    logical :: isInitialised = .false.
  contains
    procedure :: TGridCache_next_real
    procedure :: TGridCache_next_cmpl
    generic :: next => TGridCache_next_real, TGridCache_next_cmpl
    procedure :: loadEigenvecs => TGridCache_loadEigenvecs

  end type TGridCache


contains

  !> Initialises a GridCache instance.
  !! Caveat: Level index is not allowed to contain duplicate entries!
  subroutine TGridCache_init(this, env, levelIndexAll, nOrb, nAllLevel, nAllKPoint, nAllSpin, nCached,&
      & nPoints, beVerbose, eigvecBin, gridVec, origin, kPointCoords, isReal, molorb, useGpu)

    !> Structure to initialise
    class(TGridCache), intent(out) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Contains indexes (spin, kpoint, state) of the levels which must be calculated
    integer, intent(in) :: levelIndexAll(:,:)

    !> Nr. of orbitals (elements) in an eigenvectors
    integer, intent(in) :: nOrb

    !> Nr. of levels, for which eigenvector is provided in the eigenvector file
    integer, intent(in) :: nAllLevel

    !> Nr. of k-points for which eigenvectors are provided.
    integer, intent(in) :: nAllKPoint

    !> Nr. of all spins for which eigenvectors are provided.
    integer, intent(in) :: nAllSpin

    !> Nr. of cached grids
    integer, intent(in) :: nCached

    !> Nr. of grid points along the grid vectors
    integer, intent(in) :: nPoints(:)

    !> If verbosity should turned on
    logical, intent(in) :: beVerbose

    !> Name of the binary eigenvector file
    character(len=*), intent(in) :: eigvecBin

    !> Grid vectors
    real(dp), intent(in) :: gridVec(:,:)

    !> Origin of the grid
    real(dp), intent(in) :: origin(:)

    !> Coordinates of the k-points
    real(dp), intent(in) :: kPointCoords(:,:)

    !> If grids and eigenvectors are real
    logical, intent(in) :: isReal

    !> Molecular orbital calculator
    type(TMolecularOrbital), pointer, intent(in) :: molorb

    !> Whether to enable GPU offloading
    logical, intent(in) :: useGpu

    !! Contains indexes (spin, kpoint, state) to be calculated by the current MPI process
    integer, allocatable :: levelIndex(:,:)

    !! Start and end level index for MPI parallelization, if applicable
    integer :: iLvlStart, iLvlEnd

    !! Auxiliary variables
    integer ::nAll
    integer :: iSpin, iKPoint, iLevel, ind, ii, iostat
    integer :: curVec(3)
    logical :: wasFound
  
  #:if WITH_MPI
    call getStartAndEndIndex(env, size(levelIndexAll, dim=2), iLvlStart, iLvlEnd)
  #:else
    iLvlStart = 1
    iLvlEnd = size(levelIndexAll, dim=2)
  #:endif

    levelIndex = levelIndexAll(:, iLvlStart:iLvlEnd)

    @:ASSERT(.not. this%isInitialised)
    @:ASSERT(size(levelIndex, dim=1) == 3)
    @:ASSERT(size(levelIndex, dim=2) > 0)
    @:ASSERT(minval(levelIndex) > 0)
    @:ASSERT(maxval(levelIndex(1, :)) <= nAllLevel)
    @:ASSERT(maxval(levelIndex(2, :)) <= nAllKPoint)
    @:ASSERT(maxval(levelIndex(3, :)) <= nAllSpin)
    @:ASSERT(associated(molorb))
    @:ASSERT(all(shape(gridVec) == [3, 3]))
    @:ASSERT(size(origin) == 3)
    @:ASSERT(size(nPoints) == 3)
    @:ASSERT(all(shape(kPointCoords) == [3, nAllKPoint]))


    this%molorb => molorb
    this%gridVec(:,:) = gridVec
    this%origin = origin
    this%nOrb = nOrb
    this%nAllLevel = nAllLevel
    this%nAllKPoint = nAllKPoint
    this%nAllSpin = nAllSpin
    this%beVerbose = beVerbose
    allocate(this%kPoints(3, nAllKPoint))
    this%kPoints(:,:) = 2.0_dp * pi * kPointCoords
    this%nCached = nCached
    this%isReal = isReal
    this%useGpu = useGpu
    this%nPoints = nPoints

    if (this%isReal) then
      allocate(this%eigenvecReal(this%nOrb, this%nCached))
    else
      allocate(this%eigenvecCmpl(this%nOrb, this%nCached))
    end if

    nAll = size(levelIndex, dim=2)
    allocate(this%levelIndex(3, nAll))
    ! Make sure, entries are correctly sorted in the list
    ind = 1
    do iSpin = 1, nAllSpin
      do iKPoint = 1, nAllKPoint
        do iLevel = 1, nAllLevel
          curVec = [iLevel, iKPoint, iSpin]
          wasFound = .false.
          lpLevelIndex: do ii = 1, size(levelIndex, dim=2)
            wasFound = all(levelIndex(:, ii) == curVec)
            if (wasFound) exit lpLevelIndex
          end do lpLevelIndex
          if (wasFound) then
            this%levelIndex(:, ind) = curVec
            ind = ind + 1
          end if
        end do
      end do
    end do
    @:ASSERT(ind == nAll + 1)
    this%nGrid = nAll
    this%iGrid = 1
    this%cachePos = 1
    this%nReadEigVec = 0
    this%isFinished = .false.
    call openFile(this%fdEigVec, eigvecBin, mode="rb", iostat=iostat)
    if (iostat /= 0) then
      call error("Can't open file '" // trim(eigvecBin) // "'.")
    end if
    read(this%fdEigVec%unit) ii

    this%isInitialised = .true.

  end subroutine TGridCache_init


  !> Returns the next entry from the cache (real version).
  subroutine TGridCache_next_real(this, gridValReal, levelIndex, isFinished)

    !> Gridcache instance
    class(TGridCache), intent(inout) :: this

    !> Contains the molecular orbital on the grid on exit
    real(dp), pointer :: gridValReal(:,:,:)

    !> Indices of the moleular orbital (spin, kpoint, level)
    integer, intent(out) :: levelIndex(:)

    !> If all orbitals had been processed.
    logical, intent(out) :: isFinished

    complex(dp), pointer, save :: gridValCmpl(:,:,:) => null()

    call localNext(this, gridValReal, gridValCmpl, levelIndex, isFinished)

  end subroutine TGridCache_next_real


  !> Returns the next entry from the cache (complex version).
  subroutine TGridCache_next_cmpl(this, gridValCmpl, levelIndex, isFinished)

    !> Gridcache instance
    class(TGridCache), intent(inout) :: this

    !> Contains the molecular orbital on the grid on exit
    complex(dp), pointer :: gridValCmpl(:,:,:)

    !> Indices of the molecular orbital (spin, kpoint, level)
    integer, intent(out) :: levelIndex(:)

    !> If all orbitals had been processed.
    logical, intent(out) :: isFinished

    real(dp), pointer, save :: gridValReal(:,:,:) => null()

    call localNext(this, gridValReal, gridValCmpl, levelIndex, isFinished)

  end subroutine TGridCache_next_cmpl
  
  !> Loads the Eigenvectors from disk
  subroutine TGridCache_loadEigenvecs(this, iEnd)
    !> Gridcache instance
    class(TGridCache), intent(inout) :: this
    !> Nr. of eigenvectors to read
    integer, intent(in) :: iEnd
    integer :: iStartAbs
    integer :: iSpin, iKPoint, iLevel
    integer :: ind, tmp
    iStartAbs = this%iGrid ! (1)
    ind = 1
    !print *, "Loading EV for index range", iStartAbs, "to", iStartAbs + iEnd - 1

    do while (ind <= iEnd)
      if (this%isReal) then
        read(this%fdEigVec%unit) this%eigenvecReal(:,ind)
      else
        read(this%fdEigVec%unit) this%eigenvecCmpl(:,ind)
      end if
      this%nReadEigVec = this%nReadEigVec + 1

      ! If eigenvec belongs to a level which must be plotted, keep it
      iSpin = (this%nReadEigVec - 1) / (this%nAllLevel * this%nAllKPoint) + 1
      tmp = mod(this%nReadEigVec - 1, this%nAllLevel * this%nAllKPoint)
      iKPoint = tmp / this%nAllLevel + 1
      iLevel = mod(tmp, this%nAllLevel) + 1
      if (all([iLevel, iKPoint, iSpin] == this%levelIndex(:,iStartAbs+ind-1))) then
        ind = ind + 1
        if (this%beVerbose) then
          write(stdout, "(I5,I7,I7,A8)") iSpin, iKPoint, iLevel, "read"
        end if
      end if
    end do
  end subroutine TGridCache_loadEigenvecs



  !> Working subroutine for the TGridCache_next_* subroutines.
  subroutine localNext(this, gridValReal, gridValCmpl, levelIndex, isFinished)

    !> Gridcache instance
    type(TGridCache), intent(inout), target :: this

    !> Contains the real grid onexit
    real(dp), pointer :: gridValReal(:,:,:)

    !> Contains the complex grid on exit
    complex(dp), pointer :: gridValCmpl(:,:,:)

    !> Level indexes of the processed orbital on exit
    integer, intent(out) :: levelIndex(:)

    !> If all orbitals had been processed.
    logical, intent(out) :: isFinished

    integer :: iEnd, iStartAbs, iEndAbs
    real(dp), pointer :: eigReal(:,:)
    complex(dp), pointer :: eigCmpl(:,:)


    @:ASSERT(this%isInitialised)
    @:ASSERT(.not. this%isFinished)

    !! Allocate the grid cache if not done yet
    if (this%isReal .and. .not. allocated(this%gridCacheReal)) then
      allocate(this%gridCacheReal(this%nPoints(1), this%nPoints(2), this%nPoints(3), this%nCached))
    else if (.not. this%isReal .and. .not. allocated(this%gridCacheCmpl)) then
      allocate(this%gridCacheCmpl(this%nPoints(1), this%nPoints(2), this%nPoints(3), this%nCached))
    end if


    ! We passed back everything from the cache, fill it with new grids
    if (mod(this%cachePos - 1, this%nCached) == 0) then
      this%cachePos = 1
      iStartAbs = this%iGrid
      iEnd = min(this%nGrid, this%iGrid + this%nCached - 1) - this%iGrid + 1
      iEndAbs = this%iGrid + iEnd - 1

      if (this%nReadEigVec < iEndAbs) then
        call this%loadEigenvecs(iEnd)
      end if

      ! Get molecular orbital for that eigenvector
      if (this%beVerbose) then
        write(stdout, "(/,A,/)") "Calculating grid"
      end if
      if (this%isReal) then
        eigReal => this%eigenvecReal(:, :iEnd)
        call getValue(this%molorb, eigReal, this%gridCacheReal(:,:,:,:iEnd), useGpu=this%useGpu)
      else
        eigCmpl => this%eigenvecCmpl(:, :iEnd)
        call getValue(this%molorb, eigCmpl, this%kPoints,&
            & this%levelIndex(2, iStartAbs:iEndAbs), this%gridCacheCmpl(:,:,:,:iEnd), useGpu=this%useGpu)
      end if
    end if

    ! Return the appropriate grid
    if (this%isReal) then
      gridValReal => this%gridCacheReal(:,:,:,this%cachePos)
    else
      gridValCmpl => this%gridCacheCmpl(:,:,:,this%cachePos)
    end if
    levelIndex(:) = this%levelIndex(:,this%iGrid)

    ! Increase cache and grid counters
    this%iGrid = this%iGrid + 1
    this%cachePos = this%cachePos + 1
    if (this%iGrid > this%nGrid) then
      this%isFinished = .true.
      call closeFile(this%fdEigVec)
    end if
    isFinished = this%isFinished

  end subroutine localNext

end module waveplot_gridcache
