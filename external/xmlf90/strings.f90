!*******************************************************************************
! module STRINGS
! Mart Rentmeester, Mart.Rentmeester@sci.kun.nl
! http://nn-online.sci.kun.nl/fortran
! Version 1.0
!*******************************************************************************
! This is a stripped down version of the original m_strings module containing
! only the functions necessary for xml-parsing

module xmlf90_strings
  implicit none
  private

  type string
    private
    integer                 :: len = 0
    integer                 :: size = 0
    character, allocatable  :: chars(:)
  end type string

  character, parameter :: blank = ' '


  !     GENERIC PROCEDURE INTERFACE DEFINITIONS

  !---- LEN interface
  interface len
    module procedure len_s
  end interface

  !---- Conversion (to CHAR) procedure interfaces
  interface char
    module procedure s_to_c    ! string to character
    module procedure s_to_slc  ! string to specified length character
  end interface

  !---- ASSIGNMENT interfaces
  interface assignment(=)
    module procedure assign_s_to_s  ! string = string
    module procedure assign_s_to_c  ! character = string
    module procedure assign_c_to_s  ! string = character
  end interface

  !---- PREPEND_TO_STRING interface
  interface prepend_to_string
    module procedure prepend_to_string_c
    module procedure prepend_to_string_s
  end interface

  !---- APPEND_TO_STRING interface
  interface append_to_string
    module procedure append_to_string_c
    module procedure append_to_string_s
  end interface

  !---- ==  .eq. comparison operator interfaces
  interface operator(==)
    module procedure s_eq_s  ! string == string
    module procedure s_eq_c  ! string == character
    module procedure c_eq_s  ! character == string
  end interface

  !---- /=  .ne. comparison operator interfaces
  interface operator(/=)
    module procedure s_ne_s  ! string /= string
    module procedure s_ne_c  ! string /= character
    module procedure c_ne_s  ! character /= string
  end interface

  !---- TRIM interface
  interface len_trim
    module procedure len_trim_s
  end interface

  !---- LEN_TRIM interface
  interface trim
    module procedure trim_s
  end interface


  !---- Publically accessible entities
  public :: string
  public :: assignment(=)
  public :: char, len, len_trim, trim
  public :: prepend_to_string,append_to_string
  public :: operator(==),operator(/=)
  public :: resize_string

contains

  !*****************************************************************************
  !     LEN
  !*****************************************************************************

  elemental function len_s(s)
    type(string), intent(in)  :: s
    integer                   :: len_s

    len_s = s%len

  end function len_s



  !*****************************************************************************
  !     STRING_SIZE
  !*****************************************************************************

  elemental function string_size(s)
    type(string), intent(in)  :: s
    integer                   :: string_size

    string_size = s%size

  end function string_size



  !*****************************************************************************
  !     CHAR
  !*****************************************************************************
  !     Returns the characters of string as an automatically sized character

  pure function s_to_c(s)
    type(string),intent(in)   :: s
    character(len(s))         :: s_to_c

    s_to_c = transfer(s%chars(1:len(s)),s_to_c)

  end function s_to_c



  !*****************************************************************************
  !     Returns the character of fixed length, length, containing the characters
  !     of string either padded with blanks or truncated on the right to fit
  pure function s_to_slc(s,length)
    type(string),intent(in)  :: s
    integer, intent(in)      :: length
    character(length)        :: s_to_slc

    integer                  :: lc

    lc = min(len(s),length)
    s_to_slc(1:lc) = transfer(s%chars(1:lc),s_to_slc)
    !     Result longer than string: padding needed
    if (lc < length) s_to_slc(lc+1:length) = blank
    
  end function s_to_slc



  !*****************************************************************************
  ! Assign a string value to a string variable overriding default assignement.
  ! Reallocates string variable to size of string value and copies characters.

   elemental subroutine assign_s_to_s(var,expr)
    type(string), intent(inout)  :: var
    type(string), intent(in)   :: expr

    if (allocated(var%chars)) then
      if (var%size < expr%len) then
        deallocate(var%chars)
      end if
    end if
    if (.not. allocated(var%chars)) then
      allocate(var%chars(1:expr%len))
      var%size = expr%len
    end if
    var%len = expr%len
    var%chars(1:var%len) = expr%chars(1:var%len)

  end subroutine assign_s_to_s



  !*****************************************************************************
  ! Assign a string value to a character variable.
  ! If the string is longer than the character truncate the string on the right.
  ! If the string is shorter the character is blank padded on the right.

  elemental subroutine assign_s_to_c(var,expr)
    character(*), intent(out)  :: var
    type(string), intent(in)   :: expr

    integer                    :: i,lc,ls

    lc = len(var)
    ls = min(len(expr),lc)
    var(1:ls) = transfer(expr%chars(1:ls),var(1:ls))
    do i=ls+1,lc
      var(i:i) = blank
    enddo

  end subroutine assign_s_to_c



  !*****************************************************************************
  !     Assign a character value to a string variable.
  !     Disassociates the string variable from its current value, allocates new
  !     space to hold the characters and copies them from the character value
  !     into this space.

   elemental subroutine assign_c_to_s(var,expr)
    type(string), intent(inout) :: var
    character(*), intent(in)   :: expr

    integer                    :: i,lc

    lc = len(expr)
    if (allocated(var%chars)) then
      if (var%size < lc) then
        deallocate(var%chars)
      end if
    end if
    if (.not. allocated(var%chars)) then
      allocate(var%chars(1:lc))
      var%size = lc
    end if
    var%len = lc
    !var%chars(1:lc) = (/ (expr(i:i), i=1,lc) /)
    ! Workaround for gfortran 4.6
    do i = 1, lc
      var%chars(i:i) = expr(i:i)
    end do

  end subroutine assign_c_to_s



  !*****************************************************************************
  !     PREPEND_TO_STRING
  !*****************************************************************************

  pure subroutine prepend_to_string_s(s1,s2)
    type(string), intent(inout)  :: s1
    type(string), intent(in)     :: s2

    integer                      :: i,ls1,ls2
    character, allocatable :: ss(:)

    ls1 = len(s1)
    ls2 = len(s2)
    if (ls2 == 0) return
    if (ls1+ls2 > string_size(s1)) then
      allocate(ss(ls1+ls2))
      do i=1,ls2
        ss(i) = s2%chars(i)
      enddo
      do i=1,ls1
        ss(ls2+i) = s1%chars(i)
      enddo
      call move_alloc(ss, s1%chars)

      s1%len = ls1 + ls2
      s1%size = s1%len
    else
      do i=ls1,1,-1
        s1%chars(ls2+i) = s1%chars(i)
      enddo
      do i=1,ls2
        s1%chars(i) = s2%chars(i)
      enddo
      s1%len = ls1 + ls2
    endif

  end subroutine prepend_to_string_s


  
  !*****************************************************************************

  pure subroutine prepend_to_string_c(s,c)
    type(string), intent(inout)  :: s
    character(*), intent(in)     :: c

    integer                      :: i,ls,lc
    character, allocatable       :: ss(:)

    ls = len(s)
    lc = len(c)
    if (lc == 0) return
    if (ls+lc > string_size(s)) then
      allocate(ss(ls+lc))
              
      do i=1,lc
        ss(i) = c(i:i)
      enddo
      do i=1,ls
        ss(lc+i) = s%chars(i)
      enddo

      call move_alloc(ss, s%chars)

      s%len = ls + lc
      s%size = s%len
    else
      do i=ls,1,-1
        s%chars(lc+i) = s%chars(i)
      enddo
      do i=1,lc
        s%chars(i) = c(i:i)
      enddo
      s%len = ls + lc
    endif

  end subroutine prepend_to_string_c


  
  !*****************************************************************************
  !     APPEND_TO_STRING
  !*****************************************************************************

  pure subroutine append_to_string_s(s1,s2)
    type(string), intent(inout)  :: s1
    type(string), intent(in)     :: s2
    integer                      :: i,ls1,ls2

    character, allocatable       :: ss(:)

    ls1 = len(s1)
    ls2 = len(s2)
    if (ls2 == 0) return
    if (ls1+ls2 > string_size(s1)) then
      allocate(ss(ls1+ls2))

      do i=1,ls1
        ss(i) = s1%chars(i)
      enddo
      do i=ls1+1,ls1+ls2
        ss(i) = s2%chars(i-ls1)
      enddo

      call move_alloc(ss, s1%chars)

      s1%len = ls1 + ls2
      s1%size = s1%len
    else
      do i=ls1+1,ls1+ls2
        s1%chars(i) = s2%chars(i-ls1)
      enddo
      s1%len = ls1 + ls2
    endif

  end subroutine append_to_string_s



  !*****************************************************************************

  pure subroutine append_to_string_c(s,c)
    type(string), intent(inout)  :: s
    character(*), intent(in)     :: c

    integer                      :: i,ls,lc
    character, allocatable       :: ss(:)

    ls = len(s)
    lc = len(c)
    if (lc == 0) return
    if (ls+lc > string_size(s)) then
      allocate(ss(ls+lc))
              
      do i=1,ls
        ss(i) = s%chars(i)
      enddo
      do i=ls+1,ls+lc
        ss(i) = c(i-ls:i-ls)
      enddo

      call move_alloc(ss, s%chars)

      s%len = ls + lc
      s%size = s%len
    else
      do i=ls+1,ls+lc
        s%chars(i) = c(i-ls:i-ls)
      enddo
      s%len = ls + lc
    endif

  end subroutine append_to_string_c



  !*****************************************************************************
  !     LEN_TRIM
  !*****************************************************************************

  elemental function len_trim_s(s)
    type(string), intent(in)  :: s
    integer                   :: len_trim_s

    if (len(s) == 0) then
      len_trim_s = 0
      return
    endif
    do len_trim_s = len(s),1,-1
      if (s%chars(len_trim_s) /= blank) return
    end do

  end function len_trim_s



  !*****************************************************************************
  !     TRIM
  !*****************************************************************************

  pure function trim_s(s)
    type(string), intent(in)  :: s
    character(len_trim(s))    :: trim_s

    integer                   :: i

    do i=1,len(trim_s)
      trim_s(i:i) = s%chars(i)
    enddo

  end function trim_s

  !*****************************************************************************
  !     ==
  !*****************************************************************************
  ! string == string

  elemental function s_eq_s(s1,s2)
    type(string), intent(in)  :: s1,s2
    logical                   :: s_eq_s

    integer                   :: l1,l2

    l1 = len(s1)
    l2 = len(s2)
    if (l1 > l2) then
      s_eq_s = all(s1%chars(1:l2) == s2%chars) .and.  &
          all(s1%chars(l2+1:l1) == blank)
    elseif (l1 < l2) then
      s_eq_s = all(s1%chars == s2%chars(1:l1)) .and.  &
          all(blank == s2%chars(l1+1:l2))
    else
      s_eq_s = all(s1%chars == s2%chars)
    endif

  end function s_eq_s


  
  !*****************************************************************************
  ! string == character

  elemental function s_eq_c(s,c)
    type(string), intent(in)  :: s
    character(*), intent(in)  :: c

    logical                   :: s_eq_c
    integer                   :: i,ls,lc

    ls = len(s)
    lc = len(c)
    do i=1,min(ls,lc)
      if (s%chars(i) /= c(i:i)) then
        s_eq_c = .false.
        return
      endif
    enddo
    if ((ls > lc) .and. any(s%chars(lc+1:ls) /= blank)) then
      s_eq_c = .false.
    elseif ((ls < lc) .and. (blank /= c(ls+1:lc))) then
      s_eq_c = .false.
    else
      s_eq_c = .true.
    endif

  end function s_eq_c



  !*****************************************************************************
  ! character == string

  elemental function c_eq_s(c,s)
    character(*), intent(in)  :: c
    type(string), intent(in)  :: s
    
    logical                   :: c_eq_s
    integer                   :: i,lc,ls

    lc = len(c)
    ls = len(s)
    do i=1,min(lc,ls)
      if (c(i:i) /= s%chars(i)) then
        c_eq_s = .false.
        return
      endif
    enddo
    if ((lc > ls) .and. (c(ls+1:lc) /= blank)) then
      c_eq_s = .false.
    elseif ((lc < ls) .and. any(blank /= s%chars(lc+1:ls) ) )then
      c_eq_s = .false.
    else
      c_eq_s = .true.
    endif

  end function c_eq_s



  !*****************************************************************************
  !     /=
  !*****************************************************************************
  ! string /= string

  elemental function s_ne_s(s1,s2)
    type(string), intent(in)  :: s1,s2
    logical                   :: s_ne_s

    integer                   :: l1,l2

    l1 = len(s1)
    l2 = len(s2)
    if (l1 > l2) then
      s_ne_s = any(s1%chars(1:l2) /= s2%chars) .or.  &
          any(s1%chars(l2+1:l1) /= blank)
    elseif (l1 < l2) then
      s_ne_s = any(s1%chars /= s2%chars(1:l1)) .or. &
          any(blank /= s2%chars(l1+1:l2))
    else
      s_ne_s = any(s1%chars /= s2%chars)
    endif

  end function s_ne_s



  !*****************************************************************************
  ! string /= character

  elemental function s_ne_c(s,c)
    type(string), intent(in)  :: s
    character(*), intent(in)  :: c
    logical                   :: s_ne_c
    
    integer                   :: i,ls,lc

    ls = len(s)
    lc = len(c)
    do i=1,min(ls,lc)
      if (s%chars(i) /= c(i:i) )then
        s_ne_c = .true.
        return
      endif
    enddo
    if ((ls > lc) .and. any(s%chars(ls+1:lc) /= blank)) then
      s_ne_c = .true.
    elseif ((ls < lc) .and. blank /= c(ls+1:lc)) then
      s_ne_c = .true.
    else
      s_ne_c = .false.
    endif

  end function s_ne_c



  !*****************************************************************************
  ! character /= string

  elemental function c_ne_s(c,s)
    character(*), intent(in)  :: c
    type(string), intent(in)  :: s
    logical                   :: c_ne_s

    integer                   :: i,lc,ls

    lc = len(c)
    ls = len(s)
    do i=1,min(lc,ls)
      if (c(i:i) /= s%chars(i)) then
        c_ne_s = .true.
        return
      endif
    enddo
    if ((lc > ls) .and. c(ls+1:lc) /= blank) then
      c_ne_s = .true.
    elseif ((lc < ls) .and. any(blank /= s%chars(lc+1:ls))) then
      c_ne_s = .true.
    else
      c_ne_s = .false.
    endif

  end function c_ne_s



  !*****************************************************************************
  !     RESIZE_STRING procedure
  !*****************************************************************************

  !*** return code
  !*** n < 0  --> deallocate?

  pure subroutine resize_string(s,newsize)
    type(string), intent(inout)     :: s
    integer, intent(in)             :: newsize

    character, allocatable          :: c(:)
    integer                         :: i

    if (newsize <= 0) return

    if (allocated(s%chars)) then
      allocate(c(newsize))
      i = min(newsize,s%len)
      c(1:i) = s%chars(1:i)
      call move_alloc(c, s%chars)
      s%len = i
      s%size = newsize
    else
      s%size = newsize
      s%len = 0
      allocate(s%chars(s%size))
    endif

  end subroutine resize_string



end module xmlf90_strings
