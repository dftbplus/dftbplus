!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Container module for the Slater-Koster data
!>
!> This module contains the Slater-Koster tables. It decides, which one to call for which
!> species. It can be easily extended to contain different Slater-Koster schemes for different
!> species. At the moment, it handles only Slater-Koster data tabulated on an equidistant grid.
module dftbp_dftb_slakocont
  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_slakoeqgrid, only : TSlakoEqGrid, getSKIntegrals, getNIntegrals, getCutoff
#:if WITH_PLUGINS
  use dftbp_plugins_plugin, only: TPlugin
#:endif
  implicit none

  private
  public :: TSlakoCont, init
  public :: addTable, getMIntegrals, getCutoff, getSKIntegrals

  !> A specific Slater-Koster table implementation.
  type TSlaKo_
    integer :: iType = 0
    type(TSlakoEqGrid), allocatable :: pSlakoEqGrid
  end type TSlaKo_


  !> Container for Slater-Koster integrals for all pair-interactions
  type TSlakoCont
    private
    type(TSlaKo_), allocatable :: slakos(:,:)
    integer :: nAtom
    integer, allocatable :: species(:)
    integer :: mInt
    real(dp) :: cutoff
    logical :: tDataOK
    logical :: tInit = .false.
    logical :: isH = .false.
  #:if WITH_PLUGINS
    type(TPlugin), pointer, public :: plugin => null()
  #:endif
  end type TSlakoCont


  !> Initialises SlakoCont
  interface init
    module procedure SlakoCont_init
  end interface init


  !> Adds a Slater-Koster table for a given diatomic pair to the container.
  interface addTable
    module procedure SlakoCont_addTableEqGrid
  end interface addTable


  !> Returns the maximal number of integrals needed for the interactions.
  interface getMIntegrals
    module procedure SlakoCont_getMIntegrals
  end interface getMIntegrals


  !> Returns the cutoff for all interactions
  interface getCutoff
    module procedure SlakoCont_getCutoff
  end interface getCutoff


  !> Returns the Slater-Koster integrals for a given distance for a given species pair.
  interface getSKIntegrals
    module procedure SlakoCont_getSKIntegrals
  end interface getSKIntegrals

contains


  !> Initialises SlakoCont
  subroutine SlakoCont_init(this, nAtom, species, isH)

    !> SlakoCont instance
    type(TSlakoCont), intent(out) :: this

    !> Nr. of atoms in the system.
    integer, intent(in) :: nAtom

    !> Species index for each atom. Shape: [nAtom].
    integer, intent(in) :: species(:)

    !> Specified the container for Hamitonian (==.true.) or Overlap (==.false.).
    logical, intent(in) :: isH

    @:ASSERT(.not. this%tInit)
    @:ASSERT(size(species) == nAtom)

    this%nAtom = nAtom
    allocate(this%species(nAtom))
    this%species(:) = species(:)
    allocate(this%slakos(nAtom, nAtom))
    this%mInt = 0
    this%cutoff = 0.0_dp
    this%tDataOK = .false.
    this%tInit = .true.
    this%isH = isH

  end subroutine SlakoCont_init


  !> Adds a Slater-Koster table for a given diatomic species pair to the container.
  subroutine SlakoCont_addTableEqGrid(this, pTable, iSp1, iSp2)

    !> SlakoCont instance
    type(TSlakoCont), intent(inout) :: this

    !> Slater-Koster table to be added
    type(TSlakoEqGrid), allocatable, intent(inout) :: pTable

    !> Index of the first interacting species
    integer, intent(in) :: iSp1

    !> Index of the second interacting species
    integer, intent(in) :: iSp2

    integer :: iAt1, iAt2

    @:ASSERT(this%tInit)
    this%mInt = max(this%mInt, getNIntegrals(pTable))
    this%cutoff = max(this%cutoff, getCutoff(pTable))
    do iAt1 = 1, this%nAtom
      if (this%species(iAt1) /= iSp1) cycle
      do iAt2 = 1, this%nAtom
        if (this%species(iAt2) /= iSp2) cycle
        this%slakos(iAt2, iAt1)%iType = 1
        allocate(this%slakos(iAt2, iAt1)%pSlakoEqGrid, source=pTable)
      end do
    end do
    deallocate(pTable)
    this%tDataOK = all(this%slakos(:,:)%iType /= 0)

  end subroutine SlakoCont_addTableEqGrid


  !> Returns the maximal number of integrals needed for describing any of the interactions in the
  !> container
  !>
  !> This subroutine is "pure", so that it can be used to determine the size of static arrays.
  pure function SlakoCont_getMIntegrals(this) result(mInt)

    !> SlakoCont instance
    type(TSlakoCont), intent(in) :: this

    !> Max. number of integrals.
    integer :: mInt

    !! Pure procedures can not contain any I/O, therefore the following assertion is commented out
    !@:ASSERT(this%tInit)
    mInt = this%mInt

  end function SlakoCont_getMIntegrals


  !> Returns the cutoff for all interactions
  function SlakoCont_getCutoff(this) result(cutoff)

    !> SlakoCont instance
    type(TSlakoCont), intent(in) :: this

    !> Cutoff of interaction
    real(dp) :: cutoff

    @:ASSERT(this%tInit)
    cutoff = this%cutoff

  end function SlakoCont_getCutoff


  !> Returns the Slater-Koster integrals for a given distance for a given atom pair.
  subroutine SlakoCont_getSKIntegrals(this, sk, dist, atom1, atom2, sp1, sp2)

    !> SlakoCont instance
    type(TSlakoCont), intent(inout) :: this

    !> Contains the integrals on exit
    real(dp), intent(out) :: sk(:)

    !> Distance of the two atoms
    real(dp), intent(in) :: dist

    !> Index of the first atom
    integer, intent(in) :: atom1

    !> Index of the second atom
    integer, intent(in) :: atom2

    !> Index of the first interacting species.
    integer, intent(in) :: sp1

    !> Index of the second interacting species.
    integer, intent(in) :: sp2

    logical :: updated

  #:if WITH_PLUGINS
    if (associated(this%plugin)) then
      if (this%plugin%capabilities%provides_updateSKIntegrals) then
        updated = this%plugin%updateSKIntegrals(this%slakos(atom2, atom1)%pSlakoEqGrid%skTab, dist,&
            & atom1, atom2, sp1, sp2, this%isH, this%slakos(atom2, atom1)%pSlakoEqGrid%dist)
      end if
    end if
  #:endif

    @:ASSERT(this%tInit .and. this%tDataOK)
    call getSKIntegrals(this%slakos(atom2, atom1)%pSlakoEqGrid, sk, dist)

  end subroutine SlakoCont_getSKIntegrals

end module dftbp_dftb_slakocont
