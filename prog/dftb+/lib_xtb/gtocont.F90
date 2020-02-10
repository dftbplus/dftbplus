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
  end type gaussInput


  !> Gaussian Integral Container
  type :: TGTOCont
    type(TGaussFunc), allocatable :: cgto(:, :)
    integer :: moment = 0
    integer :: mShell = 0
    integer :: nSpecies = 0
    integer :: maxInt = 0
    real(dp) :: cutoff = 0.0_dp
    real(dp), allocatable :: slaterExp(:, :)
    real(dp), allocatable :: poly(:, :)
  contains
    !> Construct the integral container from input data
    procedure :: initialize
    procedure :: getGaussIntegrals
    procedure :: getMaxIntegrals
    procedure :: getCutoff
  end type TGTOCont


contains


  !> Construct the integral container from input data
  subroutine initialize(self, mShell, nSpecies, input)

    !> Gaussian integral container
    class(TGTOCont), intent(out) :: self

    !> Nr. of species in the system.
    integer, intent(in) :: nSpecies

    !> Maximum number of shells per species in the system.
    integer, intent(in) :: mShell

    !> Input to initialize integral container
    type(gaussInput), intent(in) :: input

    self%mShell = mShell
    self%nSpecies = nSpecies
    allocate(self%cgto(mShell, nSpecies), source=TGaussFunc())

  end subroutine initialize


  subroutine setGaussFunction(self, iSp1, iSh1, ang1)
    !> Gaussian Integral Container
    class(TGTOCont), intent(inout) :: self
    !> Species
    integer, intent(in) :: iSp1
    integer, intent(in) :: iSh1
    integer, intent(in) :: ang1
  end subroutine setGaussFunction

  !> Returns the cutoff for all interactions
  pure function getCutoff(self) result(cutoff)
    !> SlakoCont instance
    class(TGTOCont), intent(in) :: self
    !> Cutoff of interaction
    real(dp) :: cutoff

    cutoff = self%cutoff

  end function getCutoff

  !> Returns the maximal number of integrals needed for describing any of the
  !  interactions in the container
  pure function getMaxIntegrals(self) result(maxInt)
    !> SlakoCont instance
    class(TGTOCont), intent(in) :: self
    !> Max. number of integrals.
    integer :: maxInt

    maxInt = self%maxInt

  end function getMaxIntegrals

  subroutine getGaussIntegrals(self, ints, vec, dist, iSp1, iSh1, iSp2, iSh2)
    !> Gaussian Integral Container
    class(TGTOCont), intent(in) :: self
    !> Contains the integrals on exit
    real(dp), intent(out) :: ints(:)
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

  end subroutine getGaussIntegrals

end module dftbp_gtocont
