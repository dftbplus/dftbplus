!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines to deal with HSD species lists
module dftbp_specieslist
  use dftbp_accuracy, only : dp
  use dftbp_constants, only : elementSymbol
  use dftbp_hsdutils, only : getChildValue, getChild
  use dftbp_xmlf90, only : fnode
  implicit none
  private

  public :: readSpeciesList


  !> Generic wrapper for species specific data
  interface readSpeciesList
    module procedure :: readSpeciesListReal
  end interface readSpeciesList


contains


  !> Read a list of real valued species data
  subroutine readSpeciesListReal(node, speciesNames, array, default, conv)

    !> Node to process
    type(fnode), pointer :: node

    !> Data array to read
    real(dp), intent(out) :: array(:)

    !> Names of all species
    character(len=*), intent(in) :: speciesNames(:)

    !> Data array to read
    real(dp), intent(in), optional :: default(:)

    !> Conversion factor
    real(dp), intent(in), optional :: conv

    type(fnode), pointer :: child
    real(dp) :: fact, dummy
    integer :: iSp

    if (present(conv)) then
      fact = 1.0_dp / conv
    else
      fact = 1.0_dp
    end if

    if (present(default)) then
      do iSp = 1, size(speciesNames)
        call getChildValue(node, speciesNames(iSp), array(iSp), default=default(iSp)*fact)
      end do
    else
      do iSp = 1, size(speciesNames)
        call getChildValue(node, speciesNames(iSp), array(iSp))
      end do
    end if

    do iSp = 1, size(elementSymbol)
      call getChild(node, elementSymbol(iSp), child, requested=.false.)
      if (associated(child)) then
        call getChildValue(child, "", dummy)
      end if
    end do

    if (present(conv)) then
        array(:) = array * conv
    end if

  end subroutine readSpeciesListReal


end module dftbp_specieslist
