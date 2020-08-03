!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module dftbp_degeneracyfind
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_message, only : error
  implicit none

  private
  public :: TDegeneracyFind

  type :: TDegeneracyFind

    private

    !> Ranges of degenerate groups of levels
    integer, allocatable :: levelRange(:,:)
    
    !> To which group of states any particular level belongs
    integer, allocatable :: degenerateGroup(:)

    !> Numerical tolerance for deciding degeneracy
    real(dp) :: tolerance

    ! redundant variables (could test from size), but makes code simpler to read in a few places:

    !> Number of groups of degenerate orbitals
    integer :: nGrp
    
  contains

    !> Initialises instance and set some optional parameters
    procedure :: init

    !> Set up degeneracy test on levels
    procedure :: degeneracyTest
    
    !> Are a pair of states in the same degenerate group
    procedure :: areDegenerate

    !> Return number of degenerate levels
    procedure :: degenerateGroups

    !> Returns range of levels in each degenerate group
    procedure :: degenerateRanges
    
  end type TDegeneracyFind

  !> a few times eps, just in case of minor symmetry breaking
  real(dp), parameter :: toleranceDefault = 128.0_dp*epsilon(0.0_dp)

contains


  !> Initialise the structure
  subroutine init(this, tolerance)

    !> Instance
    class(TDegeneracyFind), intent(out) :: this

    !> Tolerance for degeneracy testing
    real(dp), intent(in), optional :: tolerance

    if (present(tolerance)) then
      this%tolerance = tolerance
    else
      this%tolerance = toleranceDefault
    end if

  end subroutine init


  !> Set up bookeeping
  subroutine degeneracyTest(this, levels, tDegenerate)

    !> Instance
    class(TDegeneracyFind), intent(inout) :: this

    !> Eigenvalues
    real(dp), intent(in) :: levels(:)

    !> Are degenerate pairs present requiring transformation
    logical, intent(out), optional :: tDegenerate

    integer :: maxRange, nOrb
    
    nOrb = size(levels)

    if (allocated(this%levelRange)) then
      if (any(shape(this%levelRange) /= [2, nOrb])) then
        deallocate(this%levelRange)
      end if
    end if
    if (.not.allocated(this%levelRange)) then
      allocate(this%levelRange(2, nOrb))
    end if
    this%levelRange(:,:) = 0
    if (allocated(this%degenerateGroup)) then
      if (size(this%degenerateGroup) /= nOrb) then
        deallocate(this%degenerateGroup)
      end if
    end if
    if (.not.allocated(this%degenerateGroup)) then
      allocate(this%degenerateGroup(nOrb))
    end if

    call degeneracyRanges(this%levelRange, this%nGrp, Levels, this%tolerance,&
        & grpMembership=this%degenerateGroup)

    maxRange = maxval(this%levelRange(2,:this%nGrp) - this%levelRange(1,:this%nGrp)) + 1
    
    if (present(tDegenerate)) then
      tDegenerate = .false.
    end if
    if (maxRange == 1) then
      ! no transformations required
      ! also nGrp == nOrb
      return
    end if
    if (present(tDegenerate)) then
      tDegenerate = .true.
    end if

  end subroutine degeneracyTest


  !> Returns whether states are in the same degenerate group
  pure function areDegenerate(this, ii, jj)

    !> Instance
    class(TDegeneracyFind), intent(in) :: this

    !> First state
    integer, intent(in) :: ii

    !> second state
    integer, intent(in) :: jj

    !> Resulting test
    logical :: areDegenerate

    areDegenerate = (this%degenerateGroup(ii) == this%degenerateGroup(jj))

  end function areDegenerate


  !> Count of degenerate groups of levels
  pure function degenerateGroups(this)

    !> Instance
    class(TDegeneracyFind), intent(in) :: this

    !> Number of degenerate groups of levels
    integer :: degenerateGroups

    degenerateGroups = this%nGrp
    
  end function degenerateGroups


  !> Returns ranges of levels in degenerate groups (lower:upper, groups)
  pure function degenerateRanges(this)

    !> Instance
    class(TDegeneracyFind), intent(in) :: this

    !> Number of degenerate groups of levels
    integer :: degenerateRanges(2,this%nGrp)

    degenerateRanges(:,:) = this%levelRange(:2,:this%nGrp)
    
  end function degenerateRanges
  
  
  ! internal routines

  !> Find which groups of eigenvales are degenerate to within a tolerance
  !> Note, similar process is used in Casida excited state calculations, so should spin off as its
  !> own module at some point
  subroutine degeneracyRanges(levelRange, nGrp, Levels, tol, grpMembership, levelSubRange)

    !> Index array for lower and upper states in degenerate group
    integer, intent(out) :: levelRange(:,:)

    !> Number of degenerate groups
    integer, intent(out) :: nGrp

    !> Eigenvalues for degeneracy testing
    real(dp), intent(in) :: Levels(:)

    !> Tolerance for degeneracy testing
    real(dp), intent(in), optional :: tol

    !> Which group each level belongs to
    integer, intent(out), optional :: grpMembership(:)

    !> sub range of eigenvalues to process
    integer, intent(in), optional :: levelSubRange(2)
    
    integer :: ii, jj, nOrb, iStart, iEnd
    real(dp) :: localTol

    if (present(tol)) then
      localTol = tol
    else
      localTol = toleranceDefault
    end if
    nOrb = size(levels)
    levelRange(:,:) = 0
    
    if (present(levelSubRange)) then
      ! set states before group as not of interest
      iStart = levelSubRange(1)
      iEnd = levelSubRange(2)
    else
      ! full range
      iStart = 1
      iEnd = nOrb
    end if
  @:ASSERT(iStart <= iEnd)
    do ii = 1, iStart - 1
      levelRange(:, ii) = ii
    end do
    
    nGrp = iStart - 1

    do ii = 1, nGrp
      grpMembership(ii) = ii
    end do

    ii = iStart
    do while (ii <= iEnd)
      nGrp = nGrp + 1
      levelRange(1, nGrp) = ii
      grpMembership(ii) = nGrp
      do jj = ii + 1, iEnd
        ! assume sorted:
        if ( abs(levels(jj) - levels(jj-1)) > localTol) then
          exit
        end if
        grpMembership(jj) = nGrp
      end do
      ii = jj
      levelRange(2, nGrp) = jj - 1
    end do

    ! set states after grouped cases as not degenerate
    do ii = iEnd + 1, nOrb
      nGrp = nGrp + 1
      levelRange(:, nGrp) = ii
      grpMembership(ii) = nGrp
    end do

  end subroutine degeneracyRanges
  
end module dftbp_degeneracyfind
