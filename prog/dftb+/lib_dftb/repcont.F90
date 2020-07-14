!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Container module for the repulsive data
!>
!> This module contains the repulsive functions. It decides, which one to call for which type
!> pairs. It can be easily extended to contain different repulsive schemes for different pairs. At
!> the moment, it handles only repulsive with spline interpolation.
module dftbp_repcont
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_repspline
  use dftbp_reppoly
  implicit none
  private

  public :: TRepCont, init
  public :: addRepulsive, getCutoff, getEnergy, getEnergyDeriv


  !> Types of repulsive currently supported
  integer, parameter :: typeRepInvalid = 0
  integer, parameter :: typeRepSpline = 1
  integer, parameter :: typeRepPoly = 2


  !> Contains repulsive types.
  type PRep_
    integer :: iType = typeRepInvalid
    type(TRepSpline), allocatable :: pRepSpline
    type(TRepPoly), allocatable :: pRepPoly
  end type PRep_


  !> Contains the repulsive interactions for the species pairs.
  type TRepCont
    private

    !> repulsive functions
    type(PRep_), allocatable :: repulsives(:,:)

    !> number of chemical species
    integer :: nSpecies

    !> Max. cutoff among all repulsives
    real(dp) :: cutoff

    !> All repulsives added.
    logical :: tDataOK

    !> Is structure initialised?
    logical :: tInit = .false.
  end type TRepCont


  !> Initialises the repulsive container.
  interface init
    module procedure RepCont_init
  end interface init


  !> Adds a new repulsive function for a given pair.
  interface addRepulsive
    module procedure RepCont_addRepSpline
    module procedure RepCont_addRepPoly
  end interface addRepulsive


  !> Returns global repulsive cutoff.
  interface getCutoff
    module procedure RepCont_getCutoff
  end interface getCutoff


  !> Returns the repulsive energy for a given distance and species pair.
  interface getEnergy
    module procedure RepCont_getEnergy
  end interface getEnergy


  !> Returns the repulsive gradient for a given distance and species pair.
  interface getEnergyDeriv
    module procedure RepCont_getEnergyDeriv
  end interface getEnergyDeriv

contains


  !> Initialises the repulsive container.
  subroutine RepCont_init(self, nSpecies)

    !> Repulsive container.
    type(TRepCont), intent(out) :: self

    !> Nr. of species.
    integer, intent(in) :: nSpecies

    @:ASSERT(.not. self%tInit)

    self%nSpecies = nSpecies
    allocate(self%repulsives(nSpecies, nSpecies))
    self%cutoff = 0.0_dp
    self%tDataOK = .false.
    self%tInit = .true.

  end subroutine RepCont_init


  !> Adds a spline repulsive function to the container for a given species pair.
  subroutine RepCont_addRepSpline(self, pRep, iSp1, iSp2)

    !> Repulsive container.
    type(TRepCont), intent(inout) :: self

    !> Repulsive function to add.
    type(TRepSpline), intent(in) :: pRep

    !> Nr. of the first interacting species.
    integer, intent(in) :: iSp1

    !> Nr. of the second interacting species.
    integer, intent(in) :: iSp2

    @:ASSERT(self%tInit)
    self%repulsives(iSp2, iSp1)%iType = typeRepSpline
    self%repulsives(iSp2, iSp1)%pRepSpline = pRep
    self%tDataOK = all(self%repulsives(:,:)%iType /= typeRepInvalid)
    self%cutoff = max(self%cutoff, getCutoff(pRep))

  end subroutine RepCont_addRepSpline


  !> Adds a polynomial repulsive function to the container for a given species pair.
  subroutine RepCont_addRepPoly(self, pRep, iSp1, iSp2)

    !> Repulsive container.
    type(TRepCont), intent(inout) :: self

    !> Repulsive function to add.
    type(TRepPoly), intent(in) :: pRep

    !> Nr. of the first interacting species.
    integer, intent(in) :: iSp1

    !> Nr. of the second interacting species.
    integer, intent(in) :: iSp2

    @:ASSERT(self%tInit)
    self%repulsives(iSp2, iSp1)%iType = typeRepPoly
    self%repulsives(iSp2, iSp1)%pRepPoly = pRep
    self%tDataOK = all(self%repulsives(:,:)%iType /= typeRepInvalid)
    self%cutoff = max(self%cutoff, getCutoff(pRep))

  end subroutine RepCont_addRepPoly


  !> Returns a global cutoff for all repulive functions.
  function RepCont_getCutoff(self) result(cutoff)

    !> Repulsive container.
    type(TRepCont), intent(in) :: self

    !> Global cutoff.
    real(dp) :: cutoff

    @:ASSERT(self%tInit .and. self%tDataOK)
    cutoff = self%cutoff

  end function RepCont_getCutoff


  !> Returns the repulsive energy for a given distance and species pair.
  subroutine RepCont_getEnergy(self, res, rr, sp1, sp2)

    !> Repulsive container.
    type(TRepCont), intent(in) :: self

    !> Energy contribution.
    real(dp), intent(out) :: res

    !> Distance between the atoms
    real(dp), intent(in) :: rr

    !> Type of the first interacting atom
    integer, intent(in) :: sp1

    !> Type of the second interacting atom
    integer, intent(in) :: sp2

    @:ASSERT(self%tInit .and. self%tDataOK)

    select case (self%repulsives(sp2, sp1)%iType)
    case(typeRepSpline)
      call getEnergy(self%repulsives(sp2, sp1)%pRepSpline, res, rr)
    case(typeRepPoly)
      call getEnergy(self%repulsives(sp2, sp1)%pRepPoly, res, rr)
    end select

  end subroutine RepCont_getEnergy


  !> Returns the repulsive gradient for a given distance and species pair.
  subroutine RepCont_getEnergyDeriv(self, res, xx, sp1, sp2)

    !> Repulsive container.
    type(TRepCont), intent(in) :: self

    !> Gradient on exit.
    real(dp), intent(out) :: res(:)

    !> Difference vector between the interacting atoms
    real(dp), intent(in) :: xx(:)

    !> Type of the first interacting atom
    integer, intent(in) :: sp1

    !> Type of the second interacting atom
    integer, intent(in) :: sp2

    @:ASSERT(self%tInit .and. self%tDataOK)
    @:ASSERT(size(res) == 3)
    @:ASSERT(size(xx) == 3)

    select case (self%repulsives(sp2, sp1)%iType)
    case(typeRepSpline)
      call getEnergyDeriv(self%repulsives(sp2, sp1)%pRepSpline, res, xx)
    case(typeRepPoly)
      call getEnergyDeriv(self%repulsives(sp2, sp1)%pRepPoly, res, xx)
    end select

  end subroutine RepCont_getEnergyDeriv

end module dftbp_repcont
