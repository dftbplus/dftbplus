!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#! (TYPE, RANK, NAME) tuple for all data types
#:set DATA_TYPES = [('integer', 1, 'I1'), ('integer', 2, 'I2'), ('integer', 3, 'I3'),&
  & ('real(dp)', 1, 'R1'), ('real(dp)', 2, 'R2'), ('real(dp)', 3, 'R3'),&
  & ('complex(dp)', 1, 'C1'), ('complex(dp)', 2, 'C2'), ('complex(dp)', 3, 'C3')]

#! Convert an integer into a fortran dimension specification
#:def ranksuffix(RANK)
$:'' if RANK == 0 else ':' + ',:' * (RANK - 1)
#:enddef ranksuffix

#! Fortran array slices, using START and END for the ranges
#:def slice(RANK, START, END)
$:''.join(["%s(%i):%s(%i), " % (START, (s+1), END, (s+1)) for s in range(RANK-1)])&
    & + "%s(%i):%s(%i)" % (START, RANK, END, RANK)
#:enddef slice

#! Fortran sizing statement from an array dimension
#:def sizes(RANK, ARRAY)
$:''.join(["size(%s, dim=%i), " % (ARRAY, (s+1)) for s in range(RANK-1)])&
    & + "size(%s, dim=%i)" % (ARRAY, RANK)
#:enddef sizes

#! Fortran array slices, using START and END for the ranges
#:def sizerange(RANK, START, END)
$:''.join(["%s(%i)-%s(%i)+1, " % (END, (s+1), START, (s+1)) for s in range(RANK-1)])&
    & + "%s(%i)-%s(%i)+1" % (END, RANK, START, RANK)
#:enddef sizerange

!> Interface module to allow read only arrays by wrapping pointers, the target can then be kept
!> privately by creating module/type by read only access given elsewhere, tracking whether it has
!> been updated by the whatever can access the private target
module dftbp_wrappointer
  use dftbp_accuracy, only : dp
  implicit none

#:for TYPE, RANK, NAME in DATA_TYPES

  !> Target array for ${TYPE}$ (${RANK}$ dim)
  type :: TWrappedTarget${NAME}$

    private

    !> Count of whether the target has been updated by whatever has access to it
    integer :: updates = 0

    !> data payload
    ${TYPE}$, allocatable :: content(${ranksuffix(RANK)}$)

  contains

    !> Mechanism to update contents of target, re-sizing if needed
    procedure :: update => update${NAME}$

  end type TWrappedTarget${NAME}$


  !> Wrapped pointer to target
  type TWrappedPointer${NAME}$

    private

    !> The real pointer
    type(TWrappedTarget${NAME}$), pointer :: payload => null()

    !> Expected state of the target, updating to match the actual state once its been seen
    integer :: expectedUpdate = 0

  contains

    !> Has the contents of the target been seen since it was last updated
    procedure :: tainted => tainted${NAME}$

    !> Mark the target as seen
    procedure :: untaint => untaint${NAME}$

    !> Access contents of the target
    generic :: data => data${NAME}$, dataSlice${NAME}$
    procedure :: data${NAME}$, dataSlice${NAME}$

    !> Shape of the target
    procedure :: shape => shape${NAME}$

  end type TWrappedPointer${NAME}$

#:endfor

  !> Intialise the target and pointer pair
  interface TWrapped_associate
  #:for _, _, NAME in DATA_TYPES
    module procedure TWrapped${NAME}$_associate
  #:endfor
  end interface TWrapped_associate


  private
#:for _, _, NAME in DATA_TYPES
  public :: TWrappedTarget${NAME}$, TWrappedPointer${NAME}$
#:endfor
  public :: TWrapped_associate

contains

  ! Wrapped target datatype functions

#:for TYPE, RANK, NAME in DATA_TYPES

  !> Update the wrapped array
  subroutine update${NAME}$(this, array)

    !> Instance of target structure
    class(TWrappedTarget${NAME}$), intent(inout) :: this

    !> Data to update target
    ${TYPE}$, intent(in) :: array(${ranksuffix(RANK)}$)

    if (allocated(this%content)) then
      if (any(shape(this%content) /= shape(array))) then
        deallocate(this%content)
      end if
    end if
    this%content = array

    this%updates = this%updates + 1

  end subroutine update${NAME}$


  ! Wrapped pointer datatype functions

  !> Check whether the internal state of the data matches the expected state
  pure function tainted${NAME}$(this) result(tainted)

    !> Instance
    class(TWrappedPointer${NAME}$), intent(in) :: this

    !> True if internal state does not match expected state
    logical :: tainted

    tainted = .true.
    if (this%expectedUpdate == this%payload%updates) then
      tainted = .false.
    end if

  end function tainted${NAME}$


  !> Clear taint marker for data - this has been seen by the client routine with access to the
  !> wrapped pointer
  subroutine untaint${NAME}$(this)

    !> Instance
    class(TWrappedPointer${NAME}$), intent(inout) :: this

    this%expectedUpdate = this%payload%updates

  end subroutine untaint${NAME}$


  !> Return a copy of the whole target array
  pure function data${NAME}$(this) result(data)

    !> Instance
    class(TWrappedPointer${NAME}$), intent(in) :: this

    ${TYPE}$ :: data(${sizes(RANK, "this%payload%content")}$)

    data(${ranksuffix(RANK)}$) = this%payload%content

  end function data${NAME}$


  !> Return slice of target array
  pure function dataSlice${NAME}$(this, start, end) result(data)

    !> Instance
    class(TWrappedPointer${NAME}$), intent(in) :: this

  #:if RANK == 1

    !> Slice index start
    integer, intent(in) :: start

    !> Slice index end
    integer, intent(in) :: end

    ${TYPE}$ :: data(end-start+1)

    data(:) = this%payload%content(start:end)

  #:else

    !> Slice index start
    integer, intent(in) :: start(:)

    !> Slice index end
    integer, intent(in) :: end(:)

    ${TYPE}$ :: data(${sizerange(RANK, "start", "end")}$)

    data(${ranksuffix(RANK)}$) = this%payload%content(${slice(RANK, "start", "end")}$)

  #:endif

  end function dataSlice${NAME}$


  !> Shape of the target array
  pure function shape${NAME}$(this) result(dataShape)

    !> Instance
    class(TWrappedPointer${NAME}$), intent(in) :: this

    integer :: dataShape(${RANK}$)

    dataShape(:) = shape(this%payload%content)

  end function shape${NAME}$

#:endfor


  ! combined mutual initialisation routines over both target and pointer for specific types

#:for _, _, NAME in DATA_TYPES

  !> Load reference to internal data into externally readable structure
  subroutine TWrapped${NAME}$_associate(internalFacing, externalFacing)

    !> Internally facing target
    type(TWrappedTarget${NAME}$), target, intent(in) :: internalFacing

    !> Externally facing pointer
    type(TWrappedPointer${NAME}$), intent(out) :: externalFacing

    externalFacing%payload => internalFacing

    ! mark internal data as unseen so far by whatever is accessing the pointer
    externalFacing%expectedUpdate = internalFacing%updates - 1

  end subroutine TWrapped${NAME}$_associate

#:endfor

end module dftbp_wrappointer
