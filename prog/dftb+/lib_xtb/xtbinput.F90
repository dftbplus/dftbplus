!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Wrapper for xTB parametrisations
module dftbp_xtbinput
  use dftbp_accuracy, only : dp
  use dftbp_xtbparam, only : xtbCalculator, xtbParam, xtbGlobalParameter
  use dftbp_gfn1param, only : getGFN1Param, gfn1Globals
  implicit none
  private

  public :: setupGFN1Parameters
  public :: xtbInput


  type :: xtbInput
    type(xtbGlobalParameter) :: gPar
    type(xtbParam), allocatable :: param(:)
  end type xtbInput


contains


  !> Initialize parametrisation data with GFN1-xTB parameters
  subroutine setupGFN1Parameters(self, speciesNames)

    !> Instance of the xTB parametrisation data
    type(xtbInput), intent(inout) :: self

    !> Element symbols for all species
    character(len=*), intent(in) :: speciesNames(:)

    integer :: iSp, nSpecies

    nSpecies = size(speciesNames, dim=1)

    ! global parameters first
    self%gPar = gfn1Globals

    ! create a dummy array of parameters
    allocate(self%param(nSpecies))
    do iSp = 1, nSpecies
      self%param(iSp) = getGFN1Param(speciesNames(iSp))
    end do

  end subroutine setupGFN1Parameters


  !> Initialize parametrisation data from input
  subroutine setupXTBCalculator(self, input)

    !> Instance of the xTB parametrisation data
    type(xtbCalculator), intent(out) :: self

    !> Instance of the xTB parametrisation data
    type(xtbInput), intent(in) :: input

    integer :: iSp, nSpecies, mShell

    nSpecies = size(input%param, dim=1)
    mShell = maxval(input%param%nSh)

    ! global parameters first
    self%gPar = input%gPar

    ! allocate space
    call self%allocateSpecies(nSpecies)
    call self%allocateShells(mShell)

    self%nShell = input%param%nSh
    self%electronegativity = input%param%en
    self%atomicRad = input%param%radius
    self%eta = input%param%eta
    self%gam = input%param%gam
    self%alpha = input%param%alpha
    self%zeff = input%param%zeff
    self%xbond = input%param%xbond
    self%anisoRad = input%param%anisoRad
    self%dipoleKernel = input%param%dpolc
    self%quadrupoleKernel = input%param%qpolc
    self%valenceCN = input%param%valenceCN

    do iSp = 1, nSpecies
      self%lShell(:self%nShell(iSp), iSp) = input%param(iSp)%basis(:self%nShell(iSp))%l
      self%shellPoly(:self%nShell(iSp), iSp) = input%param(iSp)%basis(:self%nShell(iSp))%poly
      self%shellEta(:self%nShell(iSp), iSp) = input%param(iSp)%basis(:self%nShell(iSp))%lEta
      self%level(:self%nShell(iSp), iSp) = input%param(iSp)%basis(:self%nShell(iSp))%h
      self%kcn(:self%nShell(iSp), iSp) = input%param(iSp)%basis(:self%nShell(iSp))%kcn
      self%slaterExp(:self%nShell(iSp), iSp) = input%param(iSp)%basis(:self%nShell(iSp))%zeta
      self%valenceShell(:self%nShell(iSp), iSp) = input%param(iSp)%basis(:self%nShell(iSp))%valence == 1
    end do

  end subroutine setupXTBCalculator


end module dftbp_xtbinput
