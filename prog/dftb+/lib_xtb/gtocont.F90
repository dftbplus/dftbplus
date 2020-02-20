!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module dftbp_gtocont
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_constants, only: pi
  use dftbp_gtoints
  use dftbp_slater
  implicit none

  public :: TGTOCont, gaussInput


  !> Input to generate the integral container
  type :: gaussInput

    !> Electronegativity scaling factor
    real(dp) :: kENScale

    !> Default for weighting the core Hamiltonian elements with the STO exponents
    logical :: tExponentWeighting

    !> Electronegativity for each species
    real(dp), allocatable :: electronegativity(:)

    !> Atomic radii for each species
    real(dp), allocatable :: atomicRad(:)

    !> Atomic energy for each shell and *species*
    real(dp), allocatable :: atomEnergy(:, :)

    !> Shell polynomials for each shell and species
    real(dp), allocatable :: shellPoly(:, :)

    !> Slater exponent of each shell and species
    real(dp), allocatable :: slaterExp(:, :)

    !> Status of valence character for each shell and species
    logical, allocatable :: valenceShell(:, :)
  end type gaussInput


  !> Gaussian Integral Container
  type, extends(gaussInput) :: TGTOCont
    !> Contracted Gaussian basis functions
    type(TGaussFunc), allocatable :: cgto(:, :)

    !> Nr. of atoms in the system
    integer :: nAtom

    !> Nr. of species in the system
    integer :: nSpecies

    !> Maximum number of shells per species in the system
    integer :: mShell

    integer :: maxInt = 0
    real(dp) :: cutoff = 0.0_dp

    !> Atomic energy for each shell and *atom*
    real(dp), allocatable :: selfEnergy(:, :)

    !> Shell pair parameters
    real(dp), allocatable :: h0Scale(:, :, :, :)

  contains

    !> Construct the integral container from input data
    procedure :: initialize

    procedure :: getOverlapIntegrals

    procedure :: getMaxIntegrals

    procedure :: getRCutoff
  end type TGTOCont


contains


  !> Construct the integral container from input data
  subroutine initialize(self, nAtom, nSpecies, mShell, input)

    !> Gaussian integral container
    class(TGTOCont), intent(out) :: self

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Nr. of species in the system
    integer, intent(in) :: nSpecies

    !> Maximum number of shells per species in the system
    integer, intent(in) :: mShell

    !> Input to initialize integral container
    type(gaussInput), intent(in) :: input

    self%nAtom = nAtom
    self%mShell = mShell
    self%nSpecies = nSpecies
    allocate(self%cgto(mShell, nSpecies))
    allocate(self%selfEnergy(mShell, nAtom))
    allocate(self%h0Scale(mShell, mShell, nSpecies, nSpecies))

    self%gaussInput = input

  end subroutine initialize


  !> Returns the cutoff for all interactions
  pure function getRCutoff(self) result(cutoff)

    !> SlakoCont instance
    class(TGTOCont), intent(in) :: self

    !> Cutoff of interaction
    real(dp) :: cutoff

    cutoff = self%cutoff

  end function getRCutoff


  !> Returns the maximal number of integrals needed for describing any of the
  !  interactions in the container
  pure function getMaxIntegrals(self) result(maxInt)

    !> Gaussian Integral Container
    class(TGTOCont), intent(in) :: self

    !> Max. number of integrals.
    integer :: maxInt

    maxInt = self%maxInt

  end function getMaxIntegrals


  subroutine getOverlapIntegrals(self, ints, vec, dist, iSp1, iSh1, iSp2, iSh2)

    !> Gaussian Integral Container
    class(TGTOCont), intent(in) :: self

    !> Contains the integrals on exit
    real(dp), intent(out) :: ints(:, :)

    !> Distance of the two atoms
    real(dp), intent(in) :: vec(3)

    !> Distance of the two atoms
    real(dp), intent(in) :: dist

    !> Index of the first interacting species
    integer, intent(in) :: iSp1

    !> Index of the first interacting shell
    integer, intent(in) :: iSh1

    !> Index of the second interacting species
    integer, intent(in) :: iSp2

    !> Index of the second interacting shell
    integer, intent(in) :: iSh2

    @:ASSERT(allocated(self%cgto))

    call shellPairOverlapIntegral(self%cgto(iSh1, iSp1), self%cgto(iSh2, iSp2), &
        & vec, dist, ints)

  end subroutine getOverlapIntegrals


end module dftbp_gtocont
