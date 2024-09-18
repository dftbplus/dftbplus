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
  use dftbp_io_message, only : error
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
    integer :: nSpecies
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
  subroutine SlakoCont_init(this, nSpecies, isH)

    !> SlakoCont instance
    type(TSlakoCont), intent(out) :: this

    !> Nr. of species in the system.
    integer, intent(in) :: nSpecies

    !> Specified the container for Hamitonian (==.true.) or Overlap (==.false.).
    logical, intent(in) :: isH

    @:ASSERT(.not. this%tInit)

    this%nSpecies = nSpecies
    allocate(this%slakos(nSpecies, nSpecies))
    this%mInt = 0
    this%cutoff = 0.0_dp
    this%tDataOK = .false.
    this%tInit = .true.
    this%isH = isH

  end subroutine SlakoCont_init


  !> Adds a Slater-Koster table for a given diatomic pair to the container.
  subroutine SlakoCont_addTableEqGrid(this, pTable, iSp1, iSp2)

    !> SlakoCont instance
    type(TSlakoCont), intent(inout) :: this

    !> Slater-Koster table to be added
    type(TSlakoEqGrid), allocatable, intent(inout) :: pTable

    !> Index of the first interacting species
    integer, intent(in) :: iSp1

    !> Index of the second interacting species
    integer, intent(in) :: iSp2

    @:ASSERT(this%tInit)
    this%mInt = max(this%mInt, getNIntegrals(pTable))
    this%cutoff = max(this%cutoff, getCutoff(pTable))
    this%slakos(iSp2, iSp1)%iType = 1
    call move_alloc(pTable, this%slakos(iSp2, iSp1)%pSlakoEqGrid)
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


  !> Returns the Slater-Koster integrals for a given distance for a given species pair.
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

  #:if WITH_PLUGINS
    if (associated(this%plugin)) then
      if (this%plugin%capabilities%provides_getSKIntegrals) then
        if (.not. this%plugin%getSKIntegrals(this%slakos(sp2, sp1)%pSlakoEqGrid%skTab, dist, atom1,&
            & atom2, sp1, sp2, this%isH, this%slakos(sp2, sp1)%pSlakoEqGrid%dist)) then
          call error("Cannot fetch SK integrals from plugin")
        end if
      end if
    end if
  #:endif

    @:ASSERT(this%tInit .and. this%tDataOK)
    call getSKIntegrals(this%slakos(sp2, sp1)%pSlakoEqGrid, sk, dist)

  end subroutine SlakoCont_getSKIntegrals

end module dftbp_dftb_slakocont
