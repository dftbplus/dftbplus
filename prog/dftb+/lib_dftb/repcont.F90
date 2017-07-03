!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Container module for the repulsive data
!!* @desc This module contains the repulsive functions. It decides, which
!!* one to call for which type pairs. It can be easily extended to contain
!!* different repulsive schemes for different pairs. At the moment,
!!* it handles only repulsive with spline interpolation.
module repcont
  use assert
  use accuracy
  use repspline
  use reppoly
  implicit none
  private

  public :: ORepCont, init
  public :: addRepulsive, getCutoff, getEnergy, getEnergyDeriv

  integer, parameter :: typeRepInvalid = 0
  integer, parameter :: typeRepSpline = 1
  integer, parameter :: typeRepPoly = 2

  !!* Contains repulsive types.
  type PRep_
    integer :: iType = typeRepInvalid
    type(ORepSpline), allocatable :: pRepSpline
    type(ORepPoly), allocatable :: pRepPoly
  end type PRep_


  !!* Contains the repulsive interactions for the species pairs.
  type ORepCont
    private
    type(PRep_), allocatable :: repulsives(:,:)   ! repulsive functions
    integer :: nSpecies
    real(dp) :: cutoff                        ! Max. cutoff among all repulsives
    logical :: tDataOK                        ! All repulsives added.
    logical :: tInit = .false.
  end type ORepCont


  !!* Initialises the repulsive container.
  interface init
    module procedure RepCont_init
  end interface

  !!* Adds a new repulsive function for a given pair.
  interface addRepulsive
    module procedure RepCont_addRepSpline
    module procedure RepCont_addRepPoly
  end interface

  !!* Returns global repulsive cutoff.
  interface getCutoff
    module procedure RepCont_getCutoff
  end interface

  !!* Returns the repulsive energy for a given distance and species pair.
  interface getEnergy
    module procedure RepCont_getEnergy
  end interface

  !!* Returns the repulsive gradient for a given distance and species pair.
  interface getEnergyDeriv
    module procedure RepCont_getEnergyDeriv
  end interface


contains

  !!* Initialises the repulsive container.
  !!* @param self Repulsive container.
  !!* @param nSpecies Nr. of species.
  subroutine RepCont_init(self, nSpecies)
    type(ORepCont), intent(out) :: self
    integer, intent(in) :: nSpecies

    @:ASSERT(.not. self%tInit)

    self%nSpecies = nSpecies
    allocate(self%repulsives(nSpecies, nSpecies))
    self%cutoff = 0.0_dp
    self%tDataOK = .false.
    self%tInit = .true.

  end subroutine RepCont_init


  !!* Adds a spline repulsive function to the container for a given species
  !!* pair.
  !!* @param self Repulsive container.
  !!* @param pRep Repulsive function to add.
  !!* @param iSp1 Nr. of the first interacting species.
  !!* @param iSp2 Nr. of the second interacting species.
  subroutine RepCont_addRepSpline(self, pRep, iSp1, iSp2)
    type(ORepCont), intent(inout) :: self
    type(ORepSpline), intent(in) :: pRep
    integer, intent(in) :: iSp1, iSp2

    @:ASSERT(self%tInit)
    self%repulsives(iSp2, iSp1)%iType = typeRepSpline
    self%repulsives(iSp2, iSp1)%pRepSpline = pRep
    self%tDataOK = all(self%repulsives(:,:)%iType /= typeRepInvalid)
    self%cutoff = max(self%cutoff, getCutoff(pRep))

  end subroutine RepCont_addRepSpline



  !!* Adds a polynomial repulsive function to the container for a given species
  !!* pair.
  !!* @param self Repulsive container.
  !!* @param pRep Repulsive function to add.
  !!* @param iSp1 Nr. of the first interacting species.
  !!* @param iSp2 Nr. of the second interacting species.
  subroutine RepCont_addRepPoly(self, pRep, iSp1, iSp2)
    type(ORepCont), intent(inout) :: self
    type(ORepPoly), intent(in) :: pRep
    integer, intent(in) :: iSp1, iSp2

    @:ASSERT(self%tInit)
    self%repulsives(iSp2, iSp1)%iType = typeRepPoly
    self%repulsives(iSp2, iSp1)%pRepPoly = pRep
    self%tDataOK = all(self%repulsives(:,:)%iType /= typeRepInvalid)
    self%cutoff = max(self%cutoff, getCutoff(pRep))

  end subroutine RepCont_addRepPoly



  !!* Returns a global cutoff for all repulive functions.
  !!* @param self Repulsive container.
  !!* @return Global cutoff.
  function RepCont_getCutoff(self) result(cutoff)
    type(ORepCont), intent(in) :: self
    real(dp) :: cutoff

    @:ASSERT(self%tInit .and. self%tDataOK)
    cutoff = self%cutoff

  end function RepCont_getCutoff



  !!* Returns the repulsive energy for a given distance and species pair.
  !!* @param self Repulsive container.
  !!* @param res Energy contribution.
  !!* @param rr Distance between the atoms
  !!* @param sp1 Type of the first interacting atom
  !!* @param sp2 Type of the second interacting atom
  subroutine RepCont_getEnergy(self, res, rr, sp1, sp2)
    type(ORepCont), intent(in) :: self
    real(dp), intent(out) :: res
    real(dp), intent(in) :: rr
    integer, intent(in) :: sp1, sp2

    @:ASSERT(self%tInit .and. self%tDataOK)

    select case (self%repulsives(sp2, sp1)%iType)
    case(typeRepSpline)
      call getEnergy(self%repulsives(sp2, sp1)%pRepSpline, res, rr)
    case(typeRepPoly)
      call getEnergy(self%repulsives(sp2, sp1)%pRepPoly, res, rr)
    end select

  end subroutine RepCont_getEnergy



  !!* Returns the repulsive gradient for a given distance and species pair.
  !!* @param self Repulsive container.
  !!* @param res Gradient on exit.
  !!* @param xx Difference vector between the interacting atoms
  !!* @param sp1 Type of the first interacting atom
  !!* @param sp2 Type of the second interacting atom
  subroutine RepCont_getEnergyDeriv(self, res, xx, sp1, sp2)
    type(ORepCont), intent(in) :: self
    real(dp), intent(out) :: res(:)
    real(dp), intent(in) :: xx(:)
    integer, intent(in) :: sp1, sp2

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


end module repcont
