!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Routines to deal with HSD species lists
module dftbp_dftbplus_specieslist
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : elementSymbol
  use dftbp_common_status, only : TStatus
  use dftbp_extlibs_xmlf90, only : fnode
  use dftbp_io_hsdutils, only : getChildValue, getChild
  implicit none

  private
  public :: readSpeciesList


  !> Generic wrapper for species specific data
  interface readSpeciesList
    module procedure :: readSpeciesListReal
    module procedure :: readSpeciesListInt
  end interface readSpeciesList


contains


  !> Read a list of real valued species data
  subroutine readSpeciesListReal(node, speciesNames, array, errStatus, default, conv)

    !> Node to process
    type(fnode), pointer :: node

    !> Names of all species
    character(len=*), intent(in) :: speciesNames(:)

    !> Data array to read
    real(dp), intent(out) :: array(:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

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
        call getChildValue(node, speciesNames(iSp), array(iSp), errStatus,&
            & default=default(iSp)*fact)
        @:PROPAGATE_ERROR(errStatus)
      end do
    else
      do iSp = 1, size(speciesNames)
        call getChildValue(node, speciesNames(iSp), array(iSp), errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end do
    end if

    do iSp = 1, size(elementSymbol)
      call getChild(node, elementSymbol(iSp), child, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(child)) then
        call getChildValue(child, "", dummy, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    end do

    if (present(conv)) then
        array(:) = array * conv
    end if

  end subroutine readSpeciesListReal


  !> Read a list of integer valued species data
  subroutine readSpeciesListInt(node, speciesNames, array, errStatus, default)

    !> Node to process
    type(fnode), pointer :: node

    !> Names of all species
    character(len=*), intent(in) :: speciesNames(:)

    !> Data array to read
    integer, intent(out) :: array(:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !> Data array to read
    integer, intent(in), optional :: default(:)

    type(fnode), pointer :: child
    integer :: dummy
    integer :: iSp

    if (present(default)) then
      do iSp = 1, size(speciesNames)
        call getChildValue(node, speciesNames(iSp), array(iSp), errStatus, default=default(iSp))
        @:PROPAGATE_ERROR(errStatus)
      end do
    else
      do iSp = 1, size(speciesNames)
        call getChildValue(node, speciesNames(iSp), array(iSp), errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end do
    end if

    do iSp = 1, size(elementSymbol)
      call getChild(node, elementSymbol(iSp), child, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(child)) then
        call getChildValue(child, "", dummy, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    end do

  end subroutine readSpeciesListInt


end module dftbp_dftbplus_specieslist
