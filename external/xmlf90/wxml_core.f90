module xmlf90_wxml_core

use xmlf90_wxml_buffer
use xmlf90_wxml_elstack
use xmlf90_wxml_dictionary

implicit none

logical, private, save  :: pcdata_advance_line_default = .false.
logical, private, save  :: pcdata_advance_space_default = .false.

integer, private, parameter ::  sp = selected_real_kind(6,30)
integer, private, parameter ::  dp = selected_real_kind(14,100)

private

type :: xmlf_t
   integer            :: lun
   type(buffer_t)     :: buffer
   type(elstack_t)    :: stack
   type(wxml_dictionary_t) :: dict
   logical            :: start_tag_closed
   logical            :: root_element_output
   logical            :: indenting_requested
end type xmlf_t

public :: xmlf_t
public :: xml_OpenFile, xml_NewElement, xml_EndElement, xml_Close
public :: xml_AddPcdata, xml_AddAttribute, xml_AddXMLDeclaration
public :: xml_AddComment, xml_AddCdataSection

public :: xml_AddArray
interface xml_AddArray
   module procedure  xml_AddArray_integer,  &
                xml_AddArray_real_dp, xml_AddArray_real_sp
end interface
private :: xml_AddArray_integer,  xml_AddArray_real_dp, xml_AddArray_real_sp

private :: get_unit
private :: add_eol
private :: write_attributes


integer, private, parameter  :: COLUMNS = 80

CONTAINS

!-------------------------------------------------------------------
subroutine xml_OpenFile(filename, xf, indent)
character(len=*), intent(in)  :: filename
type(xmlf_t), intent(inout)   :: xf
logical, intent(in), optional :: indent

integer :: iostat

call get_unit(xf%lun,iostat)
if (iostat /= 0) stop "cannot open file"
open(unit=xf%lun, file=filename, form="formatted", status="replace", &
     action="write", position="rewind") ! , recl=65536)

call reset_elstack(xf%stack)
call reset_dict(xf%dict)
call reset_buffer(xf%buffer)

xf%start_tag_closed = .true.
xf%root_element_output = .false.

xf%indenting_requested = .false.
if (present(indent)) then
   xf%indenting_requested = indent
endif
end subroutine xml_OpenFile

!-------------------------------------------------------------------
subroutine xml_AddXMLDeclaration(xf,encoding)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in), optional :: encoding

if (present(encoding)) then
   call add_to_buffer("<?xml version=""1.0"" encoding=""" &
                     // trim(encoding) // """ ?>", xf%buffer)
else
   call add_to_buffer("<?xml version=""1.0"" ?>", xf%buffer)
endif
end subroutine xml_AddXMLDeclaration

!-------------------------------------------------------------------
subroutine xml_AddComment(xf,comment)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in)  :: comment

call close_start_tag(xf,">")
call add_eol(xf)
call add_to_buffer("<!--", xf%buffer)
call add_to_buffer(comment, xf%buffer)
call add_to_buffer("-->", xf%buffer)
end subroutine xml_AddComment

!-------------------------------------------------------------------
subroutine xml_AddCdataSection(xf,cdata)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in)  :: cdata

call close_start_tag(xf,">")
call add_to_buffer("<![CDATA[", xf%buffer)
call add_to_buffer(cdata, xf%buffer)
call add_to_buffer("]]>", xf%buffer)
end subroutine xml_AddCdataSection

!-------------------------------------------------------------------
subroutine xml_NewElement(xf,name)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in)  :: name

if (is_empty(xf%stack)) then
   if (xf%root_element_output) stop "two root elements"
   xf%root_element_output = .true.
endif

call close_start_tag(xf,">")
call push_elstack(name,xf%stack)
call add_eol(xf)
call add_to_buffer("<" // trim(name),xf%buffer)
xf%start_tag_closed = .false.
call reset_dict(xf%dict)

end subroutine xml_NewElement
!-------------------------------------------------------------------
subroutine xml_AddPcdata(xf,pcdata,space,line_feed)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in)  :: pcdata
logical, intent(in), optional  :: space
logical, intent(in), optional  :: line_feed

logical :: advance_line , advance_space
integer :: n, i, jmax
integer, parameter   :: chunk_size = 128

advance_line = pcdata_advance_line_default 
if (present(line_feed)) then
   advance_line = line_feed
endif

advance_space = pcdata_advance_space_default 
if (present(space)) then
   advance_space = space
endif

if (is_empty(xf%stack)) then
   stop "pcdata outside element content"
endif

call close_start_tag(xf,">")

if (advance_line) then
   call add_eol(xf)
   advance_space = .false.
else
   if (xf%indenting_requested) then
      if ((len(xf%buffer) + len_trim(pcdata) + 1) > COLUMNS ) then
         call add_eol(xf)
         advance_space = .false.
      endif
   endif
endif
if (advance_space) call add_to_buffer(" ",xf%buffer)
if (len(xf%buffer) > 0) call dump_buffer(xf,lf=.false.)
!
! We bypass the buffer for the bulk of the dump
!
n = len(pcdata)
!print *, "writing pcdata of length: ", n
i = 1
do
   jmax = min(i+chunk_size-1,n)
!   print *, "writing chunk: ", i, jmax
   write(unit=xf%lun,fmt="(a)",advance="no") pcdata(i:jmax)
   if (jmax == n) exit
   i = jmax + 1
enddo
end subroutine xml_AddPcdata

!-------------------------------------------------------------------
subroutine xml_AddAttribute(xf,name,value)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in)  :: name
character(len=*), intent(in)  :: value

if (is_empty(xf%stack)) then
   stop "attributes outside element content"
endif

if (xf%start_tag_closed)  then
   stop "attributes outside start tag"
endif
if (has_key(xf%dict,name)) then
   stop "duplicate att name"
endif

call add_key_to_dict(trim(name),xf%dict)
call add_value_to_dict(trim(value),xf%dict)

end subroutine xml_AddAttribute

!-----------------------------------------------------------
subroutine xml_EndElement(xf,name)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in)  :: name

character(len=100)  :: current

if (is_empty(xf%stack)) then
   stop "Out of elements to close"
endif

call get_top_elstack(xf%stack,current)
if (current /= name) then
   print *, "current, name: ", trim(current), " ", trim(name)
   stop
endif
if (.not. xf%start_tag_closed)  then                ! Empty element
   if (len(xf%dict) > 0) call write_attributes(xf)
   call add_to_buffer(" />",xf%buffer)
   xf%start_tag_closed = .true.
else
   call add_eol(xf)
   call add_to_buffer("</" // trim(name) // ">", xf%buffer)
endif
call pop_elstack(xf%stack,current)

end subroutine xml_EndElement

!----------------------------------------------------------------

subroutine xml_Close(xf)
type(xmlf_t), intent(in)   :: xf

write(unit=xf%lun,fmt="(a)") char(xf%buffer)
close(unit=xf%lun)

end subroutine xml_Close

!==================================================================
!-------------------------------------------------------------------
subroutine get_unit(lun,iostat)

! Get an available Fortran unit number

integer, intent(out)  :: lun
integer, intent(out)  :: iostat

integer :: i
logical :: unit_used

do i = 10, 99
   lun = i
   inquire(unit=lun,opened=unit_used)
   if (.not. unit_used) then
      iostat = 0
      return
   endif
enddo
iostat = -1
lun = -1
end subroutine get_unit

!----------------------------------------------------------
subroutine add_eol(xf)
type(xmlf_t), intent(inout)   :: xf

integer :: indent_level
character(len=100), parameter  ::  blanks =  ""

indent_level = len(xf%stack) - 1
write(unit=xf%lun,fmt="(a)") char(xf%buffer)
call reset_buffer(xf%buffer)

if (xf%indenting_requested) &
   call add_to_buffer(blanks(1:indent_level),xf%buffer)

end subroutine add_eol
!------------------------------------------------------------
subroutine dump_buffer(xf,lf)
type(xmlf_t), intent(inout)   :: xf
logical, intent(in), optional :: lf

if (present(lf)) then
   if (lf) then
      write(unit=xf%lun,fmt="(a)",advance="yes") char(xf%buffer)
   else
      write(unit=xf%lun,fmt="(a)",advance="no") char(xf%buffer)
   endif
else
   write(unit=xf%lun,fmt="(a)",advance="no") char(xf%buffer)
endif
call reset_buffer(xf%buffer)

end subroutine dump_buffer

!------------------------------------------------------------
subroutine close_start_tag(xf,s)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in)  :: s

if (.not. xf%start_tag_closed)  then
   if (len(xf%dict) > 0)  call write_attributes(xf)
   call add_to_buffer(s, xf%buffer)
   xf%start_tag_closed = .true.
endif

end subroutine close_start_tag

!-------------------------------------------------------------
subroutine write_attributes(xf)
type(xmlf_t), intent(inout)   :: xf

integer  :: i, status, size
character(len=100)  :: key, value

do i = 1, len(xf%dict)
   call get_key(xf%dict,i,key,status)
   call get_value(xf%dict,key,value,status)
   key = adjustl(key)
   value = adjustl(value)
   size = len_trim(key) + len_trim(value) + 4
   if ((len(xf%buffer) + size) > COLUMNS) call add_eol(xf)
   call add_to_buffer(" ", xf%buffer)
   call add_to_buffer(trim(key), xf%buffer)
   call add_to_buffer("=", xf%buffer)
   call add_to_buffer("""",xf%buffer)
   call add_to_buffer(trim(value), xf%buffer)
   call add_to_buffer("""", xf%buffer)
enddo

end subroutine write_attributes

!---------------------------------------------------------------
    subroutine xml_AddArray_integer(xf,a,format)
      type(xmlf_t), intent(inout)         :: xf
      integer, intent(in), dimension(:)   :: a
      character(len=*), intent(in), optional  :: format

      call close_start_tag(xf,">")
      if (len(xf%buffer) > 0) call dump_buffer(xf,lf=.true.)
      if (present(format)) then
         write(xf%lun,format) a
      else
         write(xf%lun,"(6(i12))") a
      endif
    end subroutine xml_AddArray_integer

!-------------------------------------------------------------------
    subroutine xml_AddArray_real_dp(xf,a,format)
      type(xmlf_t), intent(inout)         :: xf
      real(kind=dp), intent(in), dimension(:)   :: a
      character(len=*), intent(in), optional  :: format

      call close_start_tag(xf,">")
      if (len(xf%buffer) > 0) call dump_buffer(xf,lf=.true.)
      if (present(format)) then
         write(xf%lun,format) a
      else
         write(xf%lun,"(4(es20.12))") a
      endif
    end subroutine xml_AddArray_real_dp

!------------------------------------------------------------------
    subroutine xml_AddArray_real_sp(xf,a,format)
      type(xmlf_t), intent(inout)         :: xf
      real(kind=sp), intent(in), dimension(:)   :: a
      character(len=*), intent(in), optional  :: format

      call close_start_tag(xf,">")
      if (len(xf%buffer) > 0) call dump_buffer(xf,lf=.true.)
      if (present(format)) then
         write(xf%lun,format) a
      else
         write(xf%lun,"(4(es20.12))") a
      endif
    end subroutine xml_AddArray_real_sp

end module xmlf90_wxml_core

