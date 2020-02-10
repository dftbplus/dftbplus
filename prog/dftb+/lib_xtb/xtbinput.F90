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
  use dftbp_constants, only : maxL, pi
  use dftbp_coordinationnumber, only : cnInput
  use dftbp_gfn1param, only : getGFN1Param, gfn1Globals
  use dftbp_gtocont, only : TGTOCont, gaussInput
  use dftbp_gtoints, only : TGaussFunc
  use dftbp_repcont, only : ORepCont, addRepulsive, init
  use dftbp_repsimple, only : ORepSimple, TRepSimpleIn, init
  use dftbp_slater, only : gaussFromSlater
  use dftbp_xtbparam, only : xtbParam, xtbBasis, xtbGlobalParameter
  use dftbp_xtbcont, only : xtbCalculator
  implicit none
  private

  public :: setupGFN1Parameters, setupAtomEigVal, setupGaussCont, setupRepCont
  public :: setupXTBCalculator
  public :: xtbInput


  !> Bundle of xTB related input
  type :: xtbInput

    !> Global parameters for xTB
    type(xtbGlobalParameter) :: gPar

    !> Element specific parameters for xTB
    type(xtbParam), allocatable :: param(:)

    !> Input for the coordination number container
    type(cnInput) :: cnInput

    !> Input for the gaussian type orbital basisset
    type(gaussInput) :: gaussInput

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
  subroutine setupXTBCalculator(self, nAtom, input, latVecs)

    !> Instance of the xTB parametrisation data
    type(xtbCalculator), intent(out) :: self

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Instance of the xTB parametrisation data
    type(xtbInput), intent(in) :: input

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    integer :: iSp, nSpecies, mShell
    logical :: tPeriodic

    tPeriodic = present(latVecs)
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

    if (tPeriodic) then
      call self%cnCont%initialize(nAtom, input%cnInput, latVecs)
    else
      call self%cnCont%initialize(nAtom, input%cnInput)
    end if

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


  subroutine setupRepCont(input, repCont)
    type(xtbInput), intent(in) :: input
    type(ORepCont), intent(out) :: repCont
    integer :: iSp1, iSp2
    type(ORepSimple) :: rep
    real(dp) :: zeff, alpha
    integer :: nSpecies

    nSpecies = size(input%param)
    call init(repCont, nSpecies)
    do iSp1 = 1, nSpecies
      do iSp2 = 1, nSpecies
        zeff = input%param(iSp1)%zeff * input%param(iSp2)%zeff
        alpha = sqrt(input%param(iSp1)%alpha * input%param(iSp2)%alpha)
        call init(rep, TRepSimpleIn(zeff, alpha, 1.5_dp, 1.0_dp, 40.0_dp))
        call addRepulsive(repCont, rep, iSp1, iSp2)
      end do
    end do

  end subroutine setupRepCont


  subroutine setupGaussCont(input, gtoCont)
    type(xtbInput), intent(in) :: input
    type(TGTOCont), intent(inout) :: gtoCont

    call setupGaussBasis(input%param, gtoCont)

  end subroutine setupGaussCont


  !> Expand all STOs to Gaussian functions
  subroutine setupGaussBasis(param, gtoCont)

    !> xTB parametrisation data
    type(xtbParam), intent(in) :: param(:)

    !> Container for Gaussian type integrals
    type(TGTOCont), intent(inout) :: gtoCont

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
        gtoCont%cgto(iSh1, iSp1) = cgto
      end do lpExpand

      lpOrtho: do iSh1 = 1, param(iSp1)%nSh
        if (ortho(iSh1) > 0) then
          call gsOrtho(gtoCont%cgto(ortho(iSh1), iSp1), gtoCont%cgto(iSh1, iSp1))
        end if
      end do lpOrtho

    end do

  end subroutine setupGaussBasis


  !> Orthogonalize a polarisation shell against a valence shell
  subroutine gsOrtho(val, pol)

    !> Valence shell to orthogonalize against
    type(TGaussFunc), intent(in) :: val

    !> Polarisation shell to orthogonalize
    type(TGaussFunc), intent(inout) :: pol

    integer :: ii, jj
    real(dp) :: ss, ab

    @:ASSERT(val%l == pol%l)
    @:ASSERT(val%l == 0)

    ss = 0.0_dp
    do ii = 1, val%nPrim
      do jj = 1, pol%nPrim
        ab = 1.0_dp / (val%alpha(ii) + pol%alpha(jj))
        ss = ss + val%coeff(ii) * pol%coeff(jj) * sqrt(pi * ab)**3
      end do
    end do

    pol%alpha(pol%nPrim+1:pol%nPrim+val%nPrim) = val%alpha(:val%nPrim)
    pol%coeff(pol%nPrim+1:pol%nPrim+val%nPrim) = -ss*val%coeff(:val%nPrim)
    pol%nPrim = pol%nPrim + val%nPrim

    ss = 0.0_dp
    do ii = 1, pol%nPrim
      do jj = 1, pol%nPrim
        ab = 1.0_dp / (pol%alpha(ii) + pol%alpha(jj))
        ss = ss + pol%coeff(ii) * pol%coeff(jj) * sqrt(pi * ab)**3
      end do
    end do

    pol%coeff(:pol%nPrim) = pol%coeff(:pol%nPrim) / sqrt(ss)

  end subroutine gsOrtho


  !> Expand Slater type orbital to contracted Gaussian type orbital
  subroutine cgtoFromBasis(cgto, sto)

    !> xTB basisset information
    type(xtbBasis), intent(in) :: sto

    !> Contracted Gaussian type orbital
    type(TGaussFunc), intent(out) :: cgto

    integer :: info

    call gaussFromSlater(sto%nGauss, sto%n, sto%l, sto%zeta, cgto%alpha, &
        & cgto%coeff, .true., info)
    cgto%nPrim = sto%nGauss
    cgto%l = sto%l

  end subroutine cgtoFromBasis


end module dftbp_xtbinput
