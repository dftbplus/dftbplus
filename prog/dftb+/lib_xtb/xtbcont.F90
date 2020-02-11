!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides a high-level xTB container, which carries all associated containers
module dftbp_xtbcont
  use dftbp_accuracy, only : dp
  use dftbp_coordinationnumber, only : TCNCont
  use dftbp_gtocont, only : TGTOCont
  use dftbp_periodic, only : TNeighbourList
  use dftbp_xtbparam, only : xtbGlobalParameter
  implicit none
  private

  public :: xtbCalculator

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

    !> Gaussian type orbital basis set container
    type(TGTOCont) :: gtoCont

    !> Number of shells per species
    integer, allocatable :: nShell(:)

    !> Angular momentum of each shell and species
    integer, allocatable :: lShell(:, :)

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

    !> Status of valence character for each shell and species
    logical, allocatable :: valenceShell(:, :)

    !> Reference occupation numbers
    real(dp), allocatable :: referenceN0(:, :)

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
    allocate(self%shellPoly(mShell, self%nSpecies))
    allocate(self%shellEta(mShell, self%nSpecies))
    allocate(self%level(mShell, self%nSpecies))
    allocate(self%kcn(mShell, self%nSpecies))
    allocate(self%slaterExp(mShell, self%nSpecies))
    allocate(self%valenceShell(mShell, self%nSpecies))
    allocate(self%referenceN0(mShell, self%nSpecies))

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

    cutoff = max(self%cnCont%getRCutoff(), self%gtoCont%getRCutoff())

  end function getRCutoff


end module dftbp_xtbcont
