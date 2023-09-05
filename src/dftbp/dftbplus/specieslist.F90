!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines to deal with HSD species lists
module dftbp_dftbplus_specieslist
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : elementSymbol
  use dftbp_common_unitconversion, only : TUnit
  use dftbp_extlibs_xmlf90, only : char, fnode, string
  use dftbp_io_hsdutils, only : getChild, getChildValue
  use dftbp_io_hsdutils2, only : convertUnitHsd
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
  subroutine readSpeciesListReal(node, speciesNames, array, default, units)

    !> Node to process
    type(fnode), pointer :: node

    !> Data array to read
    real(dp), intent(inout) :: array(:)

    !> Names of all species
    character(len=*), intent(in) :: speciesNames(:)

    !> Optional default values of data array to be read
    real(dp), intent(in), optional :: default(:)

    !> Conversion factor
    type(TUnit), intent(in), optional :: units(:)

    type(fnode), pointer :: child
    type(string) :: modifier
    integer :: iSp

    if (present(default)) then
      if (present(units)) then
        do iSp = 1, size(speciesNames)
          call getChildValue(node, speciesNames(iSp), array(iSp), default=default(iSp),&
              & modifier=modifier, child=child, isDefaultExported=.false.)
          call convertUnitHsd(char(modifier), units, child, array(iSp))
        end do
      else
        do iSp = 1, size(speciesNames)
          call getChildValue(node, speciesNames(iSp), array(iSp), default=default(iSp),&
              & isDefaultExported=.false.)
        end do
      end if
    else
      if (present(units)) then
        do iSp = 1, size(speciesNames)
          call getChildValue(node, speciesNames(iSp), array(iSp), modifier=modifier, child=child)
          call convertUnitHsd(char(modifier), units, child, array(iSp))
        end do
      else
        do iSp = 1, size(speciesNames)
          call getChildValue(node, speciesNames(iSp), array(iSp))
        end do
      end if
    end if

  end subroutine readSpeciesListReal


  !> Read a list of integer valued species data
  subroutine readSpeciesListInt(node, speciesNames, array, default)

    !> Node to process
    type(fnode), pointer :: node

    !> Data array to read
    integer, intent(out) :: array(:)

    !> Names of all species
    character(len=*), intent(in) :: speciesNames(:)

    !> Data array to read
    integer, intent(in), optional :: default(:)

    integer :: iSp

    if (present(default)) then
      do iSp = 1, size(speciesNames)
        call getChildValue(node, speciesNames(iSp), array(iSp), default=default(iSp),&
            & isDefaultExported = .false.)
      end do
    else
      do iSp = 1, size(speciesNames)
        call getChildValue(node, speciesNames(iSp), array(iSp))
      end do
    end if

  end subroutine readSpeciesListInt


end module dftbp_dftbplus_specieslist
