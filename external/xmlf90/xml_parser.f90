module xmlf90_xml_parser

  
!
! Basic module to parse XML in the SAX spirit.
!

use xmlf90_buffer
use xmlf90_reader
use xmlf90_fsm
use xmlf90_dictionary
use xmlf90_debug
use xmlf90_xml_error
use xmlf90_elstack          ! For element nesting checks
use xmlf90_entities
!

implicit none

private

!
!  XML file handle
!
type, public :: xml_t
private
      type(file_buffer_t)  :: fb
      type(fsm_t)          :: fx
      character(len=200)   :: path_mark
end type xml_t

!
public :: xml_parse
public :: open_xmlfile, close_xmlfile
public :: endfile_xmlfile, rewind_xmlfile
public :: eof_xmlfile, sync_xmlfile
public :: xml_char_count
public :: xml_path, xml_mark_path, xml_get_path_mark
public :: xml_name, xml_attributes

CONTAINS  !=============================================================

subroutine open_xmlfile(fname,fxml,iostat,record_size)
character(len=*), intent(in)      :: fname
integer, intent(out)              :: iostat
type(xml_t), intent(out)          :: fxml
integer, intent(in), optional     :: record_size

call open_file(fname,fxml%fb,iostat,record_size)
call init_fsm(fxml%fx)
fxml%path_mark = ""

end subroutine open_xmlfile
!-------------------------------------------------------------------------

subroutine rewind_xmlfile(fxml)
type(xml_t), intent(inout) :: fxml

call rewind_file(fxml%fb)
call reset_fsm(fxml%fx)
fxml%path_mark = ""

end subroutine rewind_xmlfile

!-----------------------------------------
subroutine endfile_xmlfile(fxml)
type(xml_t), intent(inout) :: fxml

call mark_eof_file(fxml%fb) 

end subroutine endfile_xmlfile

!-----------------------------------------
subroutine close_xmlfile(fxml)
type(xml_t), intent(inout) :: fxml

call close_file_buffer(fxml%fb)
call reset_fsm(fxml%fx)         ! just in case
fxml%path_mark = ""             ! ""

end subroutine close_xmlfile

!-----------------------------------------
subroutine sync_xmlfile(fxml,iostat)
type(xml_t), intent(inout) :: fxml
integer, intent(out)       :: iostat

call sync_file(fxml%fb,iostat)
! Do not reset fx: that's the whole point of synching.

end subroutine sync_xmlfile

!----------------------------------------------------
function eof_xmlfile(fxml) result (res)
type(xml_t), intent(in)          :: fxml
logical                          :: res

res = eof_file(fxml%fb)

end function eof_xmlfile
!
!----------------------------------------------------
!----------------------------------------------------
function xml_char_count(fxml) result (nc)
type(xml_t), intent(in)          :: fxml
integer                          :: nc

nc = nchars_processed(fxml%fb)

end function xml_char_count
!
!----------------------------------------------------
!

subroutine xml_parse(fxml, begin_element_handler,    &
                           end_element_handler,      &
                           pcdata_chunk_handler,     &
                           comment_handler,          &
                           xml_declaration_handler,  &
                           cdata_section_handler,     &
                           sgml_declaration_handler, &
                           error_handler,            &
                           signal_handler,            &
                           verbose,                  &
                           empty_element_handler)

type(xml_t), intent(inout), target  :: fxml

optional                            :: begin_element_handler
optional                            :: end_element_handler
optional                            :: pcdata_chunk_handler
optional                            :: comment_handler
optional                            :: xml_declaration_handler
optional                            :: sgml_declaration_handler
optional                            :: cdata_section_handler
optional                            :: error_handler
optional                            :: signal_handler
logical, intent(in), optional       :: verbose
optional                            :: empty_element_handler

interface
   subroutine begin_element_handler(name,attributes)
   use xmlf90_dictionary
   character(len=*), intent(in)     :: name
   type(dictionary_t), intent(in)   :: attributes
   end subroutine begin_element_handler

   subroutine end_element_handler(name)
   character(len=*), intent(in)     :: name
   end subroutine end_element_handler

   subroutine pcdata_chunk_handler(chunk)
   character(len=*), intent(in) :: chunk
   end subroutine pcdata_chunk_handler

   subroutine comment_handler(comment)
   character(len=*), intent(in) :: comment
   end subroutine comment_handler

   subroutine xml_declaration_handler(name,attributes)
   use xmlf90_dictionary
   character(len=*), intent(in)     :: name
   type(dictionary_t), intent(in)   :: attributes
   end subroutine xml_declaration_handler

   subroutine sgml_declaration_handler(sgml_declaration)
   character(len=*), intent(in) :: sgml_declaration
   end subroutine sgml_declaration_handler

   subroutine cdata_section_handler(cdata)
   character(len=*), intent(in) :: cdata
   end subroutine cdata_section_handler

   subroutine error_handler(error_info)
   use xmlf90_xml_error
   type(xml_error_t), intent(in)            :: error_info
   end subroutine error_handler

   subroutine signal_handler(code)
   logical, intent(out) :: code
   end subroutine signal_handler

   subroutine empty_element_handler(name,attributes)
   use xmlf90_dictionary
   character(len=*), intent(in)     :: name
   type(dictionary_t), intent(in)   :: attributes
   end subroutine empty_element_handler

end interface

character(len=1)     :: c
integer              :: iostat, status

character(len=150)   :: message
integer              :: signal

type(buffer_t)       :: translated_pcdata
type(buffer_t)       :: name, oldname, dummy

logical              :: have_begin_handler, have_end_handler, &
                        have_pcdata_handler, have_comment_handler, &
                        have_xml_declaration_handler, &
                        have_sgml_declaration_handler, &
                        have_cdata_section_handler, have_empty_handler, &
                        have_error_handler, have_signal_handler

logical              :: pause_signal

type(xml_error_t)            :: error_info
type(file_buffer_t), pointer :: fb
type(fsm_t), pointer         :: fx

have_begin_handler = present(begin_element_handler)
have_end_handler = present(end_element_handler)
have_pcdata_handler = present(pcdata_chunk_handler)
have_comment_handler = present(comment_handler)
have_xml_declaration_handler = present(xml_declaration_handler)
have_sgml_declaration_handler = present(sgml_declaration_handler)
have_cdata_section_handler = present(cdata_section_handler)
have_error_handler = present(error_handler)
have_signal_handler = present(signal_handler)
have_empty_handler = present(empty_element_handler)

fb => fxml%fb
fx => fxml%fx
if (present(verbose)) then
   debug = verbose                 ! For m_converters
   fx%debug = verbose              ! job-specific flag
endif

if (fx%debug) print *, " Entering xml_parse..."

!---------------------------------------------------------------------
do
      call get_character(fb,c,iostat)

      if (iostat /= 0) then          ! End of file...
         if (.not. is_empty(fx%element_stack)) then
            call build_error_info(error_info, &
                 "Early end of file.", &
                 line(fb),column(fb),fx%element_stack,SEVERE_ERROR_CODE)
            if (have_error_handler) then
               call error_handler(error_info)
            else
               call default_error_handler(error_info)
            endif
         endif
         call endfile_xmlfile(fxml)  ! Mark it as eof
         EXIT
      endif

      call evolve_fsm(fx,c,signal)

      if (fx%debug) print *, c, " ::: ", trim(fx%action)

      if (signal == END_OF_TAG) then
         !
         ! We decide whether we have ended an opening tag or a closing tag
         !
         if (fx%context == OPENING_TAG) then
            name = fx%element_name

            if (fx%debug) print *, "We have found an opening tag"
            if (fx%root_element_seen) then
               if (name .equal. fx%root_element_name) then
                  call build_error_info(error_info, &
                  "Duplicate root element: " // str(name), &
                  line(fb),column(fb),fx%element_stack,SEVERE_ERROR_CODE)
                  if (have_error_handler) then
                     call error_handler(error_info)
                  else
                     call default_error_handler(error_info)
                  endif
               endif
               if (is_empty(fx%element_stack)) then
                  call build_error_info(error_info, &
                  "Opening tag beyond root context: " // str(name), &
                  line(fb),column(fb),fx%element_stack,SEVERE_ERROR_CODE)
                  if (have_error_handler) then
                     call error_handler(error_info)
                  else
                     call default_error_handler(error_info)
                  endif
               endif
            else
               fx%root_element_name = name
               fx%root_element_seen = .true.
            endif
            call push_elstack(name,fx%element_stack)
            if (have_begin_handler) &
                call begin_element_handler(str(name),fx%attributes)

         else if (fx%context == CLOSING_TAG) then
            name = fx%element_name
         
            if (fx%debug) print *, "We have found a closing tag"
            if (is_empty(fx%element_stack)) then
               call build_error_info(error_info, &
                  "Nesting error: End tag: " // str(name) //  &
                  " does not match -- too many end tags", &
                  line(fb),column(fb),fx%element_stack,SEVERE_ERROR_CODE)
               if (have_error_handler) then
                  call error_handler(error_info)
               else
                  call default_error_handler(error_info)
               endif
            else
               call get_top_elstack(fx%element_stack,oldname)
               if (oldname .equal. name) then
                  call pop_elstack(fx%element_stack,oldname)
                  if (have_end_handler) call end_element_handler(str(name))
!!                  call pop_elstack(fx%element_stack,oldname)
               else
                  call build_error_info(error_info, &
                       "Nesting error: End tag: " // str(name) //  &
                       ". Expecting end of : " // str(oldname), &
                       line(fb),column(fb),fx%element_stack,SEVERE_ERROR_CODE)
                  if (have_error_handler) then
                     call error_handler(error_info)
                  else
                     call default_error_handler(error_info)
                  endif
               endif
            endif
         else if (fx%context == SINGLE_TAG) then
            name = fx%element_name

            if (fx%debug) print *, "We have found a single (empty) tag: ", &
                                char(name)
            !
            ! Push name on to stack to reveal true xpath
            !
            call push_elstack(name,fx%element_stack)
            if (have_empty_handler) then
               if (fx%debug) print *, "--> calling empty_element_handler."
               call empty_element_handler(str(name),fx%attributes)
               call pop_elstack(fx%element_stack,dummy)
            else
               if (have_begin_handler) then
                  if (fx%debug) print *, "--> calling begin_element_handler..."
                  call begin_element_handler(str(name),fx%attributes)
               endif
               call pop_elstack(fx%element_stack,dummy)
               if (have_end_handler) then
                  if (fx%debug) print *, "--> ... and end_element_handler."
                  call end_element_handler(str(name))
               endif
            endif
!!            call pop_elstack(fx%element_stack,dummy)

         else if (fx%context == CDATA_SECTION_TAG) then

            if (fx%debug) print *, "We found a CDATA section"
            if (is_empty(fx%element_stack)) then
               if (fx%debug) print *, &
                   "... Warning: CDATA section outside element context"
            else
               if (have_cdata_section_handler) then
                  call cdata_section_handler(str(fx%pcdata))
               else
                  if (have_pcdata_handler) &
                   call pcdata_chunk_handler(str(fx%pcdata))
               endif
            endif

         else if (fx%context == COMMENT_TAG) then

            if (fx%debug) print *, "We found a comment tag"
            if (have_comment_handler)  &
                 call comment_handler(str(fx%pcdata))

         else if (fx%context == SGML_DECLARATION_TAG) then

            if (fx%debug) print *, "We found an sgml declaration"
            if (have_sgml_declaration_handler)  &
                      call sgml_declaration_handler(str(fx%pcdata))

         else if (fx%context == XML_DECLARATION_TAG) then

            if (fx%debug) print *, "We found an XML declaration"
            name = fx%element_name
            if (have_xml_declaration_handler)  &
                      call xml_declaration_handler(str(name),fx%attributes)

         else

            ! do nothing

         endif

      else if (signal == CHUNK_OF_PCDATA) then

         if (fx%debug) print *, "We found a chunk of PCDATA"
         if (is_empty(fx%element_stack)) then
            if (fx%debug) print *, "... Warning: PCDATA outside element context"
               ! Just a warning
               call build_error_info(error_info, &
                  "PCDATA outside of element context", &
                  line(fb),column(fb),fx%element_stack,WARNING_CODE)
               if (have_error_handler) then
                  call error_handler(error_info)
               else
                  call default_error_handler(error_info)
               endif
         else
            !
            ! Replace entities by their value
            !
            call entity_filter(fx%pcdata,translated_pcdata,status,message)
            if (status < 0) then
               call build_error_info(error_info, message, &
                  line(fb),-status,fx%element_stack,SEVERE_ERROR_CODE)
               if (have_error_handler) then
                  call error_handler(error_info)
               else
                  call default_error_handler(error_info)
               endif
            else if (status > 0) then
               ! Just a warning
               call build_error_info(error_info, message, &
                  line(fb),status,fx%element_stack,WARNING_CODE)
               if (have_error_handler) then
                  call error_handler(error_info)
               else
                  call default_error_handler(error_info)
               endif
            else
               if (have_pcdata_handler) &
                   call pcdata_chunk_handler(str(translated_pcdata))
            endif
         endif

      else if (signal == EXCEPTION) then
         call build_error_info(error_info, fx%action, &
                  line(fb),column(fb),fx%element_stack,SEVERE_ERROR_CODE)
         if (have_error_handler) then
            call error_handler(error_info)
         else
            call default_error_handler(error_info)
         endif
      else
         ! QUIET, do nothing
      endif
      if (signal /= QUIET) then
         if (have_signal_handler) then
            call signal_handler(pause_signal)
            if (pause_signal) exit
         endif
      endif

enddo

end subroutine xml_parse

!-----------------------------------------
subroutine xml_path(fxml,path)
type(xml_t), intent(in) :: fxml
character(len=*), intent(out)  :: path

call get_elstack_signature(fxml%fx%element_stack,path)

end subroutine xml_path

!-----------------------------------------
subroutine xml_mark_path(fxml,path)
!
! Marks the current path
!
type(xml_t), intent(inout) :: fxml
character(len=*), intent(out)  :: path

call get_elstack_signature(fxml%fx%element_stack,fxml%path_mark)
path = fxml%path_mark

end subroutine xml_mark_path

!-----------------------------------------
subroutine xml_get_path_mark(fxml,path)
!
! Returns the currently markd path (or an empty string if there are no marks)
!
type(xml_t), intent(in)        :: fxml
character(len=*), intent(out)  :: path

path = fxml%path_mark

end subroutine xml_get_path_mark

!-----------------------------------------
subroutine xml_name(fxml,name)
type(xml_t), intent(in) :: fxml
character(len=*), intent(out)  :: name

name = char(fxml%fx%element_name)

end subroutine xml_name
!-----------------------------------------
subroutine xml_attributes(fxml,attributes)
type(xml_t), intent(in) :: fxml
type(dictionary_t), intent(out)  :: attributes

attributes = fxml%fx%attributes

end subroutine xml_attributes

end module xmlf90_xml_parser




