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
  implicit none

  public :: TGaussCont, init


  !> Gaussian Integral Container
  type :: TGaussCont
    private
    type(TGaussFunc), allocatable :: cgto(:)
    integer :: moment = 0
    integer :: nSpecies = 0
    integer :: maxInt = 0
    real(dp) :: cutoff = 0.0_dp
  contains
    procedure :: getGaussIntegrals
    procedure :: getMaxIntegrals
    procedure :: getCutoff
  end type TGaussCont


  interface init
    module procedure :: GaussContInitialize
  end interface init


contains


  subroutine GaussContInitialize(self, nSpecies)
    !> Gaussian Integral Container
    type(TGaussCont), intent(out) :: self
    !> Nr. of species in the system.
    integer, intent(in) :: nSpecies

    self%nSpecies = nSpecies
    allocate(self%cgto(nSpecies), source=TGaussFunc())

  end subroutine GaussContInitialize

  !> Returns the cutoff for all interactions
  pure function getCutoff(self) result(cutoff)
    !> SlakoCont instance
    class(TGaussCont), intent(in) :: self
    !> Cutoff of interaction
    real(dp) :: cutoff

    cutoff = self%cutoff

  end function getCutoff

  !> Returns the maximal number of integrals needed for describing any of the
  !  interactions in the container
  pure function getMaxIntegrals(self) result(maxInt)
    !> SlakoCont instance
    class(TGaussCont), intent(in) :: self
    !> Max. number of integrals.
    integer :: maxInt

    maxInt = self%maxInt

  end function getMaxIntegrals

  subroutine getGaussIntegrals(self, ints, vec, dist, iSp1, iSp2)
    !> Gaussian Integral Container
    class(TGaussCont), intent(in) :: self
    !> Contains the integrals on exit
    real(dp), intent(out) :: ints(:)
    !> Distance of the two atoms
    real(dp), intent(in) :: vec(3)
    !> Distance of the two atoms
    real(dp), intent(in) :: dist
    !> Index of the first interacting species.
    integer, intent(in) :: iSp1
    !> Index of the second interacting species.
    integer, intent(in) :: iSp2

    @:ASSERT(allocated(self%cgto))

    call shellPairOverlapIntegral(self%cgto(iSp1), self%cgto(iSp2), vec, dist, &
        & ints)

  end subroutine getGaussIntegrals

end module dftbp_gtocont
