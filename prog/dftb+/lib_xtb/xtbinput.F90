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
  use dftbp_constants, only : maxL
  use dftbp_coordinationnumber, only : cnInput
  use dftbp_gfn1param, only : getGFN1Param, gfn1Globals
  use dftbp_gtocont, only : TGaussCont
  use dftbp_gtoints, only : TGaussFunc
  use dftbp_slater, only : gaussFromSlater
  use dftbp_xtbparam, only : xtbCalculator, xtbParam, xtbBasis, xtbGlobalParameter
  implicit none
  private

  public :: setupGFN1Parameters, setupAtomEigVal, setupGaussCont
  public :: xtbInput


  !> Bundle of xTB related input
  type :: xtbInput

    !> Global parameters for xTB
    type(xtbGlobalParameter) :: gPar

    !> Element specific parameters for xTB
    type(xtbParam), allocatable :: param(:)

    !> Input for the coordination number container
    type(cnInput) :: cnInput

    !> Use weighting by exponent in the construction of H0
    logical :: tExponentWeighting

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


  subroutine setupAtomEigVal(param, atomEigVal)
    type(xtbParam), intent(in) :: param(:)
    real(dp), intent(out) :: atomEigVal(:, :)
    integer :: iSp1, iSh1

    do iSp1 = 1, size(atomEigVal, dim=1)
      do iSh1 = 1, param(iSp1)%nSh
        atomEigVal(iSh1, iSp1) = param(iSh1)%basis(iSh1)%h
      end do
    end do

  end subroutine setupAtomEigVal

  subroutine setupGaussCont(input, gaussCont)
    type(xtbInput), intent(in) :: input
    type(TGaussCont), intent(out) :: gaussCont

    call setupGaussBasis(input%param, gaussCont)

  end subroutine setupGaussCont

  subroutine setupGaussBasis(param, gaussCont)
    type(xtbParam), intent(in) :: param(:)
    type(TGaussCont), intent(inout) :: gaussCont
    type(TGaussFunc) :: cgto
    integer :: ortho(maxval(param%nSh))
    integer :: valSh(0:maxL)
    integer :: iSp1, iSh1


    do iSp1 = 1, size(param)
      ortho(:) = 0
      valSh(:) = 0
      lpExpand: do iSh1 = 1, param(iSp1)%nSh
        call cgtoFromBasis(cgto, param(iSp1)%basis(iSh1))
        if (valSh(cgto%l) == 0 .and. param(iSp1)%basis(iSh1)%valence == 1) then
          valSh(cgto%l) = iSh1
        else
          ortho(iSh1) = valSh(cgto%l)
        end if
        gaussCont%cgto(iSh1, iSp1) = cgto
      end do lpExpand

      lpOrtho: do iSh1 = 1, param(iSp1)%nSh
        if (ortho(iSh1) > 0) then
          call gsOrtho(gaussCont%cgto(ortho(iSh1), iSp1), gaussCont%cgto(iSh1, iSp1))
        end if
      end do lpOrtho

    end do

  end subroutine setupGaussBasis

  subroutine gsOrtho(val, pol)
    type(TGaussFunc), intent(inout) :: val
    type(TGaussFunc), intent(inout) :: pol
  end subroutine gsOrtho

  pure subroutine cgtoFromBasis(cgto, sto)
    type(xtbBasis), intent(in) :: sto
    type(TGaussFunc), intent(out) :: cgto
    integer :: info
    call gaussFromSlater(sto%nGauss, sto%n, sto%l, sto%zeta, cgto%alpha, &
        & cgto%coeff, .true., info)
    cgto%nPrim = sto%nGauss
    cgto%l = sto%l
  end subroutine cgtoFromBasis


end module dftbp_xtbinput
