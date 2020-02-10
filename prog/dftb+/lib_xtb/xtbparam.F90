!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> xTB parameter types
module dftbp_xtbparam
  use dftbp_accuracy, only : dp
  use dftbp_coordinationnumber, only : TCNCont
  use dftbp_periodic, only : TNeighbourList
  implicit none
  private

  public :: xtbParam, xtbBasis, xtbCalculator, xtbGlobalParameter

  !> xTB-basisset definition
  type xtbBasis

    !> principal quantum number
    integer :: n

    !> angular momentum of shell
    integer :: l

    !> shell polynomial
    real(dp) :: poly

    !> chemical hardness shell parameter
    real(dp) :: lEta

    !> atomic level
    real(dp) :: h

    !> coordination number dependence of the level
    real(dp) :: kcn

    !> Slater function exponent
    real(dp) :: zeta

    !> Number of Gaussian functions for expansion
    integer :: nGauss

    !> Is valence or polarization shell.
    integer :: valence

  end type

  !> xTB-parameters
  type :: xtbParam

    !> Number of shells
    integer :: nSh

    !> Metal
    integer :: metal

    !> Electronegativity
    real(dp) :: en

    !> Atomic radius
    real(dp) :: radius

    !> Chemical hardness
    real(dp) :: eta

    !> Third order Hubbard derivative
    real(dp) :: gam

    !> Repulsion parameter
    real(dp) :: alpha

    !> Repulsion parameter
    real(dp) :: zeff

    !> Halogen bond strength
    real(dp) :: xbond = 0.0_dp

    !> Anisotropy radius
    real(dp) :: anisoRad = 0.0_dp

    !> Dipole kernel
    real(dp) :: dpolc = 0.0_dp

    !> Quadrupole kernel
    real(dp) :: qpolc = 0.0_dp

    !> Valence coordination number
    real(dp) :: valenceCN = 0.0_dp

    !> Basis functions
    type(xtbBasis), allocatable :: basis(:)

  end type xtbParam


  !> Global parameters used in xTB methods
  type :: xtbGlobalParameter

    !> Scaling parameter for EN differences
    real(dp) :: kEnScale

    !> Repulsive polynomial parameter
    real(dp) :: krep

    !> Halogen bonding parameter
    real(dp) :: xbondRad

    !> Halogen bonding parameter
    real(dp) :: xbondExp

  end type xtbGlobalParameter


  !> xTB parametrisation data
  type :: xtbCalculator

    !> Number of species
    integer :: nSpecies

    !> Maximum number of shells per atom
    integer :: mShell

    !> Global parameters
    type(xtbGlobalParameter) :: gPar

    !> Coordination number container
    type(TCNCont) :: cnCont

    !> Number of shells per species
    integer, allocatable :: nShell(:)

    !> Angular momentum of each shell and species
    integer, allocatable :: lShell(:, :)

    !> Electronegativity for each species
    real(dp), allocatable :: electronegativity(:)

    !> Atomic radius for each species
    real(dp), allocatable :: atomicRad(:)

    !> Chemical hardness for each species
    real(dp), allocatable :: eta(:)

    !> Hubbard derivative for each species
    real(dp), allocatable :: gam(:)

    !> Repulsion parameter for each species
    real(dp), allocatable :: alpha(:)

    !> Effective nuclear charge for each species
    real(dp), allocatable :: zeff(:)

    !> Halogen bond strength for each species
    real(dp), allocatable :: xbond(:)

    !> Anisotropy radius for each species
    real(dp), allocatable :: anisoRad(:)

    !> Dipole kernel for each species
    real(dp), allocatable :: dipoleKernel(:)

    !> Quadrupole kernel for each species
    real(dp), allocatable :: quadrupoleKernel(:)

    !> Valence coordination number for each species
    real(dp), allocatable :: valenceCN(:)

    !> Shell polynomials for each shell and species
    real(dp), allocatable :: shellPoly(:, :)

    !> Chemical hardness parameter for each shell and species
    real(dp), allocatable :: shellEta(:, :)

    !> Atomic level for each shell and species
    real(dp), allocatable :: level(:, :)

    !> coordination number dependence of the level
    real(dp), allocatable :: kcn(:, :)

    !> Slater exponent of each shell and species
    real(dp), allocatable :: slaterExp(:, :)

    !> Gaussian exponents of each shell and species
    real(dp), allocatable :: gaussExp(:, :, :)

    !> Gaussian contraction coefficents of each shell and species
    real(dp), allocatable :: gaussCoeff(:, :, :)

    !> Status of valence character for each shell and species
    logical, allocatable :: valenceShell(:, :)

    !> Pair parameters for each species
    real(dp), allocatable :: pairParam(:, :)

  contains

    !> Allocate space for species dependent parameters
    procedure :: allocateSpecies

    !> Allocate space for basisset dependent parameters
    procedure :: allocateShells

    !> update internal copy of coordinates
    procedure :: updateCoords

    !> update internal copy of lattice vectors
    procedure :: updateLatVecs

    !> get real space cutoff
    procedure :: getRCutoff

  end type xtbCalculator


contains

  !> Allocate space for species dependent parameters
  subroutine allocateSpecies(self, nSpecies)

    !> Instance of the xTB parametrisation data
    class(xtbCalculator), intent(inout) :: self

    !> Number of species
    integer, intent(in) :: nSpecies

    @:ASSERT(nSpecies > 0)

    self%nSpecies = nSpecies

    allocate(self%nShell(nSpecies))
    allocate(self%electronegativity(nSpecies))
    allocate(self%atomicRad(nSpecies))
    allocate(self%eta(nSpecies))
    allocate(self%gam(nSpecies))
    allocate(self%alpha(nSpecies))
    allocate(self%zeff(nSpecies))
    allocate(self%xbond(nSpecies))
    allocate(self%anisoRad(nSpecies))
    allocate(self%dipoleKernel(nSpecies))
    allocate(self%quadrupoleKernel(nSpecies))
    allocate(self%valenceCN(nSpecies))
    allocate(self%pairParam(nSpecies, nSpecies))

  end subroutine allocateSpecies


  !> Allocate space for basisset dependent parameters
  subroutine allocateShells(self, mShell)

    !> Instance of the xTB parametrisation data
    class(xtbCalculator), intent(inout) :: self

    !> Maximum number of shells
    integer, intent(in) :: mShell

    @:ASSERT(mShell > 0)
    @:ASSERT(self%nSpecies > 0)

    self%mShell = mShell

    allocate(self%lShell(mShell, self%nSpecies))
    allocate(self%lShell(mShell, self%nSpecies))
    allocate(self%shellPoly(mShell, self%nSpecies))
    allocate(self%shellEta(mShell, self%nSpecies))
    allocate(self%level(mShell, self%nSpecies))
    allocate(self%kcn(mShell, self%nSpecies))
    allocate(self%slaterExp(mShell, self%nSpecies))
    allocate(self%valenceShell(mShell, self%nSpecies))

  end subroutine allocateShells


  !> Update internal stored coordinates
  subroutine updateCoords(self, neighList, img2CentCell, coords, species0)

    !> data structure
    class(xtbCalculator), intent(inout) :: self

    !> list of neighbours to atoms
    type(TNeighbourList), intent(in) :: neighList

    !> image to central cell atom index
    integer, intent(in) :: img2CentCell(:)

    !> atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    call self%cnCont%updateCoords(neighList, img2CentCell, coords, species0)

  end subroutine updateCoords


  !> update internal copy of lattice vectors
  subroutine updateLatVecs(self, latVecs)

    !> data structure
    class(xtbCalculator), intent(inout) :: self

    !> lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    call self%cnCont%updateLatVecs(latVecs)

  end subroutine updateLatVecs


  !> Distance cut off for calculator
  function getRCutoff(self) result(cutoff)

    !> data structure
    class(xtbCalculator), intent(inout) :: self

    !> resulting cutoff
    real(dp) :: cutoff

    cutoff = self%cnCont%getRCutoff()

  end function getRCutoff


  !> Serialize xTB-parameters to TOML
  subroutine xtbParamTOML(unit, param, gfn, root)

    !> Unit for IO
    integer, intent(in) :: unit

    !> Instance
    type(xtbParam), intent(in) :: param

    !> GFN level
    integer, intent(in), optional :: gfn

    !> TOML table name
    character(len=*), intent(in), optional :: root

    !> TOML key-value format
    character(len=*), parameter :: tomlf = '(a,1x,"=",1x,g0)'

    integer :: iSh

    if (present(gfn)) then
      write(unit, tomlf) "GFN", gfn
    end if
  #:for key in ("nSh", "metal", "en", "radius", "eta", "gam", "alpha", "zeff", "xbond", "anisoRad", "dpolc", "qpolc", "valenceCN")
    write(unit, tomlf) "${key}$", param%${key}$
  #:endfor

    do iSh = 1, param%nSh
      if (present(root)) then
        write(unit, '("[[",a,".",a,"]]")') root, "basis"
      else
        write(unit, '("[[",a,"]]")') root, "basis"
      end if
    #:for key in ("n", "l", "poly", "lEta", "h", "kcn", "zeta", "nGauss", "valence")
      write(unit, tomlf) "${key}$", param%basis(iSh)%${key}$
    #:endfor
    end do

  end subroutine xtbParamTOML

end module dftbp_xtbparam
