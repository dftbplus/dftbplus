!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module dftbp_dftbplus_input_fileaccess
  use dftbp_common_file, only : fileAccessValues
  use dftbp_extlibs_xmlf90, only : fnode
  use dftbp_io_charmanip, only : tolower, unquote
  use dftbp_io_hsdutils, only : detailedError, getChild, getChildValue, setChildValue
  use dftbp_type_linkedlist, only : asArray, destruct, init, len, TListString
  implicit none

  private
  public :: readBinaryAccessTypes

contains


  !> Reads in the file acess types
  subroutine readBinaryAccessTypes(node, accessTypes)

    !> Parent note which should contain the "BinaryAccessTypes" subnode
    type(fnode), pointer, intent(in) :: node

    !> Read and write access types on exit (defaulting to ["stream", "stream"])
    character(*), intent(out) :: accessTypes(:)

    type(fnode), pointer :: child
    type(TListString) :: stringList
    integer :: ii

    @:ASSERT(size(accessTypes) == 2)

    call getChild(node, "BinaryAccessTypes", child, requested=.false.)
    if (.not. associated(child)) then
      call setChildValue(node, "BinaryAccessTypes", ["stream"], child=child)
    end if
    call init(stringList)
    call getChildValue(child, "", stringList)
    if (len(stringList) < 1 .or. len(stringList) > 2) then
      call detailedError(child, "BinaryAccessTypes needs one or two arguments")
    end if
    call asArray(stringList, accessTypes(1 : len(stringList)))
    if (len(stringList) == 1) then
      accessTypes(2) = accessTypes(1)
    end if
    call destruct(stringList)
    accessTypes(:) = tolower(unquote(accessTypes))
    do ii = 1, size(accessTypes)
      if (.not. any(accessTypes(ii) == fileAccessValues)) then
        call detailedError(child, "Invalid file access type '" // trim(accessTypes(ii)) // "'")
      end if
    end do

  end subroutine readBinaryAccessTypes

end module dftbp_dftbplus_input_fileaccess
