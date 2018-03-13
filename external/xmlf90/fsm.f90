module xmlf90_fsm
!
use xmlf90_buffer
use xmlf90_dictionary
use xmlf90_charset
use xmlf90_entities
use xmlf90_elstack

implicit none

private

type, public :: fsm_t
      !
      ! Contains information about the "finite state machine"
      ! Some of the components (marked *) could at this point be made into
      ! saved module variables.
      ! 
      !
      integer              :: state
      integer              :: context
      integer              :: nbrackets             !*
      integer              :: nlts                  !*
      character(len=1)     :: quote_char            !*
      type(buffer_t)       :: buffer                !*
      type(buffer_t)       :: element_name
      type(dictionary_t)   :: attributes
      type(buffer_t)       :: pcdata
      type(elstack_t)      :: element_stack
      logical              :: root_element_seen
      type(buffer_t)       :: root_element_name
      character(len=150)   :: action
      logical              :: debug
end type fsm_t

public :: init_fsm, reset_fsm, evolve_fsm

!
! State parameters
!
integer, parameter, public   ::  ERROR = -1
integer, parameter, public   ::  INIT = 1         
integer, parameter, private  ::  START_TAG_MARKER = 2
integer, parameter, private  ::  END_TAG_MARKER = 3
integer, parameter, private  ::  IN_NAME = 4
integer, parameter, private  ::  WHITESPACE_IN_TAG = 5
integer, parameter, private  ::  IN_PCDATA = 6
integer, parameter, private  ::  SINGLETAG_MARKER = 7
integer, parameter, private  ::  CLOSINGTAG_MARKER = 8
integer, parameter, private  ::  IN_COMMENT = 9
integer, parameter, private  ::  IN_ATT_NAME = 10
integer, parameter, private  ::  IN_ATT_VALUE = 11
integer, parameter, private  ::  EQUAL = 12
integer, parameter, private  ::  SPACE_AFTER_EQUAL = 13
integer, parameter, private  ::  SPACE_BEFORE_EQUAL = 14
integer, parameter, private  ::  START_QUOTE = 15
integer, parameter, private  ::  END_QUOTE = 16
integer, parameter, private  ::  BANG = 17
integer, parameter, private  ::  BANG_HYPHEN = 18
integer, parameter, private  ::  ONE_HYPHEN = 19
integer, parameter, private  ::  TWO_HYPHEN = 20
integer, parameter, private  ::  QUESTION_MARK = 21
integer, parameter, private  ::  START_XML_DECLARATION = 22
integer, parameter, private  ::  IN_SGML_DECLARATION = 23
integer, parameter, private  ::  IN_CDATA_SECTION = 24
integer, parameter, private  ::  ONE_BRACKET = 25
integer, parameter, private  ::  TWO_BRACKET = 26
integer, parameter, private  ::  CDATA_PREAMBLE = 27
integer, parameter, private  ::  IN_PCDATA_AT_EOL = 30
!
! Context parameters
!
integer, parameter, public   ::  OPENING_TAG  = 100
integer, parameter, public   ::  CLOSING_TAG  = 110
integer, parameter, public   ::  SINGLE_TAG   = 120
integer, parameter, public   ::  COMMENT_TAG  = 130
integer, parameter, public   ::  XML_DECLARATION_TAG  = 140
integer, parameter, public   ::  SGML_DECLARATION_TAG  = 150
integer, parameter, public   ::  CDATA_SECTION_TAG  = 160
integer, parameter, public   ::  NULL_CONTEXT          = 200
!
! Signal parameters
!
integer, parameter, public   ::  QUIET             = 1000
integer, parameter, public   ::  END_OF_TAG        = 1100
integer, parameter, public   ::  CHUNK_OF_PCDATA   = 1200
integer, parameter, public   ::  EXCEPTION         = 1500

CONTAINS

!------------------------------------------------------------
! Initialize once and for all the derived types (Fortran90 restriction)
!
subroutine init_fsm(fx) 
type(fsm_t), intent(inout)   :: fx

 fx%state = INIT
 call setup_xml_charsets()
 fx%context = NULL_CONTEXT
 call init_elstack(fx%element_stack)
 fx%root_element_seen = .false.
 fx%debug = .false.
 fx%action = ""
 call init_buffer(fx%buffer)
 call init_buffer(fx%element_name)
 call init_buffer(fx%pcdata)
 call init_buffer(fx%root_element_name)
 call init_dict(fx%attributes)
end subroutine init_fsm
!------------------------------------------------------------
subroutine reset_fsm(fx) 
type(fsm_t), intent(inout)   :: fx

 fx%state = INIT
 call setup_xml_charsets()
 fx%context = NULL_CONTEXT
 call reset_elstack(fx%element_stack)
 fx%action = ""
 fx%root_element_seen = .false.
 call reset_buffer(fx%buffer)
 call reset_buffer(fx%element_name)
 call reset_buffer(fx%pcdata)
 call reset_buffer(fx%root_element_name)
 call reset_dict(fx%attributes)
end subroutine reset_fsm

!------------------------------------------------------------
subroutine evolve_fsm(fx,c,signal)
!
! Finite-state machine evolution rules for XML parsing.
!
type(fsm_t), intent(inout)      :: fx    ! Internal state
character(len=1), intent(in)    :: c
integer, intent(out)            :: signal

!
! Reset signal
!
signal = QUIET
!

if (.not. (c .in. valid_chars)) then
!
!      Let it pass (in case the underlying encoding is UTF-8)
!      But this chars in a name will cause havoc
!
!      signal = EXCEPTION
!      fx%state = ERROR
!      fx%action = trim("Not a valid character in simple encoding: "//c)
!      RETURN
endif

select case(fx%state)

 case (INIT)
      if (c == "<") then
         fx%state = START_TAG_MARKER
         if (fx%debug) fx%action = ("Starting tag")
      else if (c == ">") then
         fx%state = ERROR
         fx%action = ("Ending tag without being in one!")
      else
         if (fx%debug) fx%action = ("Reading garbage chars")
      endif

 case (START_TAG_MARKER)
      if (c == ">") then
         fx%state = ERROR
         fx%action = ("Tag empty!")
      else if (c == "<") then
         fx%state = ERROR
         fx%action = ("Double opening of tag!!")
      else if (c == "/") then
         fx%state = CLOSINGTAG_MARKER
         if (fx%debug) fx%action = ("Starting endtag: ")
         fx%context = CLOSING_TAG
      else if (c == "?") then
         fx%state = START_XML_DECLARATION
         if (fx%debug) fx%action = ("Starting XML declaration ")
         fx%context = XML_DECLARATION_TAG
      else if (c == "!") then
         fx%state = BANG
         if (fx%debug) fx%action = ("Saw ! -- comment or SGML declaration expected...")
      else if (c .in. whitespace) then
         fx%state = ERROR
         fx%action = ("Cannot have whitespace after <")
      else if (c .in. initial_name_chars) then
         fx%context = OPENING_TAG
         fx%state = IN_NAME
         call add_to_buffer(c,fx%buffer)
         if (fx%debug) fx%action = ("Starting to read name in tag")
      else 
         fx%state = ERROR
         fx%action = ("Illegal initial character for name")
      endif


 case (BANG)
      if (c == "-") then
         fx%state = BANG_HYPHEN
         if (fx%debug) fx%action = ("Almost ready to start comment ")
      else if (c .in. uppercase_chars) then
         fx%state = IN_SGML_DECLARATION
         fx%nlts = 0
         fx%nbrackets = 0
         if (fx%debug) fx%action = ("SGML declaration ")
         fx%context = SGML_DECLARATION_TAG
         call add_to_buffer(c,fx%buffer)
      else if (c == "[") then
         fx%state = CDATA_PREAMBLE
         if (fx%debug) fx%action = ("Declaration with [ ")
         fx%context = CDATA_SECTION_TAG
      else
         fx%state = ERROR
         fx%action = ("Wrong character after ! ")
      endif

 case (CDATA_PREAMBLE)
      ! We assume a CDATA[ is forthcoming, we do not check
      if (c == "[") then
         fx%state = IN_CDATA_SECTION
         if (fx%debug) fx%action = ("About to start reading CDATA contents")
      else if (c == "]") then
         fx%state = ERROR
         fx%action = ("Unexpected ] in CDATA preamble")
      else
         if (fx%debug) fx%action = ("Reading CDATA preamble")
      endif

 case (IN_CDATA_SECTION)
      if (c == "]") then
         fx%state = ONE_BRACKET
         if (fx%debug) fx%action = ("Saw a ] in CDATA section")
      else
         call add_to_buffer(c,fx%buffer)
         if (fx%debug) fx%action = ("Reading contents of CDATA section")
      endif

 case (ONE_BRACKET)
      if (c == "]") then
         fx%state = TWO_BRACKET
         if (fx%debug) fx%action = ("Maybe finish a CDATA section")
      else
         fx%state = IN_CDATA_SECTION
         call add_to_buffer("]",fx%buffer)
         if (fx%debug) fx%action = ("Continue reading contents of CDATA section")
      endif

 case (TWO_BRACKET)
      if (c == ">") then
         fx%state = END_TAG_MARKER
         signal = END_OF_TAG
         if (fx%debug) fx%action = ("End of CDATA section")
         fx%pcdata = fx%buffer    ! Not quite the same behavior
                                  ! as pcdata... (not filtered)
         call reset_buffer(fx%buffer)
      else
         fx%state = IN_CDATA_SECTION
         call add_to_buffer("]",fx%buffer)
         if (fx%debug) fx%action = ("Continue reading contents of CDATA section")
      endif

 case (IN_SGML_DECLARATION)
      if (c == "<") then
         fx%nlts = fx%nlts + 1
         call add_to_buffer("<",fx%buffer)
         fx%action = "Read an intermediate < in SGML declaration"
      else if (c == "[") then
         fx%nbrackets = fx%nbrackets + 1
         call add_to_buffer("[",fx%buffer)
         fx%action = "Read a [ in SGML declaration"
      else if (c == "]") then
         fx%nbrackets = fx%nbrackets - 1
         call add_to_buffer("]",fx%buffer)
         fx%action = "Read a ] in SGML declaration"
      else if (c == ">") then
         if (fx%nlts == 0) then
            if (fx%nbrackets == 0) then
               fx%state = END_TAG_MARKER
               signal  = END_OF_TAG
               if (fx%debug) fx%action = ("Ending SGML declaration tag")
               fx%pcdata = fx%buffer       ! Same behavior as pcdata
               call reset_buffer(fx%buffer)
            else
               fx%state = ERROR
               fx%action = ("Unmatched ] in SGML declaration")
            endif
         else
            fx%nlts = fx%nlts -1
            call add_to_buffer(">",fx%buffer)
            fx%action = "Read an intermediate > in SGML declaration"
         endif
      else
         if (fx%debug) fx%action = ("Keep reading SGML declaration")
         call add_to_buffer(c,fx%buffer)
      endif

 case (BANG_HYPHEN)
      if (c == "-") then
         fx%state = IN_COMMENT
         fx%context = COMMENT_TAG
         if (fx%debug) fx%action = ("In comment ")
      else
         fx%state = ERROR
         fx%action = ("Wrong character after <!- ")
      endif

 case (START_XML_DECLARATION)
      if (c .in. initial_name_chars) then
         fx%state = IN_NAME
         call add_to_buffer(c,fx%buffer)
         if (fx%debug) fx%action = ("Starting to read name in XML declaration")
      else
         fx%state = ERROR
         fx%action = "Wrong character after ? in start of XML declaration"
      endif

 case (CLOSINGTAG_MARKER)
      if (c == ">") then
         fx%state = ERROR
         fx%action = ("Closing tag empty!")
      else if (c == "<") then
         fx%state = ERROR
         fx%action = ("Double opening of closing tag!!")
      else if (c == "/") then
         fx%state = ERROR
         fx%action = ("Syntax error (<//)")
      else if (c .in. whitespace) then
         fx%state = ERROR
         fx%action = ("Cannot have whitespace after </")
      else if (c .in. initial_name_chars) then
         fx%state = IN_NAME
         if (fx%debug) fx%action = ("Starting to read name inside endtag")
         call add_to_buffer(c,fx%buffer)
      else 
         fx%state = ERROR
         fx%action = ("Illegal initial character for name")
      endif

 case (IN_NAME)
      if (c == "<") then
         fx%state = ERROR
         fx%action = ("Starting tag within tag")
      else if (c == ">") then
         fx%state = END_TAG_MARKER
         signal  = END_OF_TAG
         if (fx%debug) fx%action = ("Ending tag")
!         call set_element_name(fx%buffer,fx%element)
         fx%element_name = fx%buffer
         call reset_buffer(fx%buffer)
         call reset_dict(fx%attributes)
      else if (c == "/") then
         if (fx%context /= OPENING_TAG) then
            fx%state = ERROR
            fx%action = ("Single tag did not open as start tag")
         else 
            fx%state = SINGLETAG_MARKER
            fx%context = SINGLE_TAG
            if (fx%debug) fx%action = ("Almost ending single tag")
!            call set_element_name(fx%buffer,fx%element)
            fx%element_name = fx%buffer
            call reset_buffer(fx%buffer)
            call reset_dict(fx%attributes)
         endif
      else if (c .in. whitespace) then
         fx%state = WHITESPACE_IN_TAG
         if (fx%debug) fx%action = ("Ending name chars")
!            call set_element_name(fx%buffer,fx%element)
         fx%element_name = fx%buffer
         call reset_buffer(fx%buffer)
         call reset_dict(fx%attributes)
      else if (c .in. name_chars) then
         if (fx%debug) fx%action = ("Reading name chars in tag")
         call add_to_buffer(c,fx%buffer)
      else
         fx%state = ERROR
         fx%action = ("Illegal character for name")
      endif

 case (IN_ATT_NAME)
      if (c == "<") then
         fx%state = ERROR
         fx%action = ("Starting tag within tag")
      else if (c == ">") then
         fx%state = ERROR
         fx%action = ("Ending tag in the middle of an attribute")
      else if (c == "/") then
         fx%state = ERROR
         fx%action = ("Ending tag in the middle of an attribute")
      else if (c .in. whitespace) then
         fx%state = SPACE_BEFORE_EQUAL  
         if (fx%debug) fx%action = ("Whitespace after attr. name (specs?)")
         call add_key_to_dict(fx%buffer,fx%attributes)
         call reset_buffer(fx%buffer)
      else if ( c == "=" ) then
         fx%state = EQUAL
         if (fx%debug) fx%action = ("End of attr. name")
         call add_key_to_dict(fx%buffer,fx%attributes)
         call reset_buffer(fx%buffer)
      else if (c .in. name_chars) then
         if (fx%debug) fx%action = ("Reading attribute name chars")
         call add_to_buffer(c,fx%buffer)
      else
         fx%state = ERROR
         fx%action = ("Illegal character for attribute name")
      endif

 case (EQUAL)
      if ( (c == """") .or. (c == "'") ) then
         fx%state = START_QUOTE
         if (fx%debug) fx%action = ("Found beginning quote")
         fx%quote_char = c
      else if (c .in. whitespace) then
         fx%state = SPACE_AFTER_EQUAL
         if (fx%debug) fx%action = ("Whitespace after equal sign...")
      else
         fx%state = ERROR
         fx%action = ("Must use quotes for attribute values")
      endif

 case (SPACE_BEFORE_EQUAL)
      if ( c == "=" ) then
         fx%state = EQUAL
         if (fx%debug) fx%action = ("Equal sign")
      else if (c .in. whitespace) then
         if (fx%debug) fx%action = ("More whitespace before equal sign...")
      else
         fx%state = ERROR
         fx%action = ("Must use equal sign for attribute values")
      endif

 case (SPACE_AFTER_EQUAL)
      if ( c == "=" ) then
         fx%state = ERROR
         fx%action = ("Duplicate Equal sign")
      else if (c .in. whitespace) then
         if (fx%debug) fx%action = ("More whitespace after equal sign...")
      else  if ( (c == """") .or. (c == "'") ) then
         fx%state = START_QUOTE
         fx%quote_char = c
         if (fx%debug) fx%action = ("Found beginning quote")
      else
         fx%state = ERROR
         fx%action = ("Must use quotes for attribute values")
      endif

 case (START_QUOTE)
      if (c == fx%quote_char) then
         fx%state = END_QUOTE
         if (fx%debug) fx%action = ("Emtpy attribute value...")
         call add_value_to_dict(fx%buffer,fx%attributes)
         call reset_buffer(fx%buffer)
      else if (c == "<") then
         fx%state = ERROR
         fx%action = ("Attribute value cannot contain <")
      else   ! actually allowed chars in att values... Specs: No "<"        
         fx%state = IN_ATT_VALUE
         if (fx%debug) fx%action = ("Starting to read attribute value")
         call add_to_buffer(c,fx%buffer)
      endif

 case (IN_ATT_VALUE)
      if (c == fx%quote_char) then
         fx%state = END_QUOTE
         if (fx%debug) fx%action = ("End of attribute value")
         call add_value_to_dict(fx%buffer,fx%attributes)
         call reset_buffer(fx%buffer)
      else if (c == "<") then
         fx%state = ERROR
         fx%action = ("Attribute value cannot contain <")
      else if ( (c == char(10)) ) then
         fx%state = ERROR
!
!        Aparently other whitespace is allowed...
!
         fx%action = ("No newline allowed in attr. value (specs?)")
      else        ! all other chars allowed in attr value
         if (fx%debug) fx%action = ("Reading attribute value chars")
         call add_to_buffer(c,fx%buffer)
      endif

 case (END_QUOTE)
      if ((c == """") .or. (c == "'")) then
         fx%state = ERROR
         fx%action = ("Duplicate end quote")
      else if (c .in. whitespace) then
         fx%state = WHITESPACE_IN_TAG
         if (fx%debug) fx%action = ("Space in between attributes or to end of tag")
      else if (c == "<") then
         fx%state = ERROR
         fx%action = ("Starting tag within tag")
      else if (c == ">") then
         if (fx%context == XML_DECLARATION_TAG) then
            fx%state = ERROR
            fx%action = "End of XML declaration without ?"
         else
            fx%state = END_TAG_MARKER
            signal  = END_OF_TAG
            if (fx%debug) fx%action = ("Ending tag after some attributes")
         endif
      else if (c == "/") then
         if (fx%context /= OPENING_TAG) then
            fx%state = ERROR
            fx%action = ("Single tag did not open as start tag")
         else 
            fx%state = SINGLETAG_MARKER
            fx%context = SINGLE_TAG
            if (fx%debug) fx%action = ("Almost ending single tag after some attributes")
         endif
      else if (c == "?") then
         if (fx%context /= XML_DECLARATION_TAG) then
            fx%state = ERROR
            fx%action = "Wrong lone ? in tag"
         else
            fx%state = QUESTION_MARK
            if (fx%debug) fx%action = ("About to end XML declaration")
         endif
      else   
         fx%state = ERROR
         fx%action = ("Must have some whitespace after att. value")
      endif


 case (WHITESPACE_IN_TAG)
      if ( c .in. whitespace) then
         if (fx%debug) fx%action = ("Reading whitespace in tag")
      else if (c == "<") then
         fx%state = ERROR
         fx%action = ("Starting tag within tag")
      else if (c == ">") then
         if (fx%context == XML_DECLARATION_TAG) then
            fx%state = ERROR
            fx%action = "End of XML declaration without ?"
         else
            fx%state = END_TAG_MARKER
            signal  = END_OF_TAG
            if (fx%debug) fx%action = ("End whitespace in tag")
         endif
      else if (c == "/") then
         if (fx%context /= OPENING_TAG) then
            fx%state = ERROR
            fx%action = ("Single tag did not open as start tag")
         else 
            fx%state = SINGLETAG_MARKER
            fx%context = SINGLE_TAG
            if (fx%debug) fx%action = ("End whitespace in single tag")
         endif
      else if (c .in. initial_name_chars) then
         fx%state = IN_ATT_NAME
         if (fx%debug) fx%action = ("Starting Attribute name in tag")
         call add_to_buffer(c,fx%buffer)
      else if (c == "?") then
         if (fx%context /= XML_DECLARATION_TAG) then
            fx%state = ERROR
            fx%action = "Wrong lone ? in tag"
         else
            fx%state = QUESTION_MARK
            if (fx%debug) fx%action = ("About to end XML declaration after whitespace")
         endif
      else
         fx%state = ERROR
         fx%action = ("Illegal initial character for attribute")
      endif

 case (QUESTION_MARK)
      if (c == ">") then
         fx%state = END_TAG_MARKER
         signal  = END_OF_TAG
         if (fx%debug) fx%action = ("End of XML declaration tag")
      else
         fx%state = ERROR
         fx%action = "No > after ? in XML declaration tag"
      endif

 case (IN_COMMENT)
      !
      ! End of comment is  "-->", and  ">" can appear inside comments
      !
      if (c == "-") then
         fx%state = ONE_HYPHEN
         if (fx%debug) fx%action = ("Saw - in Comment")
      else
         if (fx%debug) fx%action = ("Reading comment")
         call add_to_buffer(c,fx%buffer)
      endif

 case (ONE_HYPHEN)
      if (c == "-") then
         fx%state = TWO_HYPHEN
         if (fx%debug) fx%action = ("About to end comment")
      else
         fx%state = IN_COMMENT
         if (fx%debug) fx%action = ("Keep reading comment after -: ")
         call add_to_buffer("-",fx%buffer)
         call add_to_buffer(c,fx%buffer)
      endif

 case (TWO_HYPHEN)
      if (c == ">") then
         fx%state = END_TAG_MARKER
         signal  = END_OF_TAG
         if (fx%debug) fx%action = ("End of Comment")
         fx%pcdata = fx%buffer                  ! Same behavior as pcdata
         call reset_buffer(fx%buffer)
      else
         fx%state = ERROR
         fx%action = ("Cannot have -- in comment")
      endif

 case (SINGLETAG_MARKER)

      if (c == ">") then
         fx%state = END_TAG_MARKER
         signal  = END_OF_TAG
         if (fx%debug) fx%action = ("Ending tag")
         ! We have to call begin_element AND end_element
      else 
         fx%state = ERROR
         fx%action = ("Wrong ending of single tag")
      endif

 case (IN_PCDATA)
      if (c == "<") then
         fx%state = START_TAG_MARKER
         signal = CHUNK_OF_PCDATA
         if (fx%debug) fx%action = ("End of pcdata -- Starting tag")
         fx%pcdata = fx%buffer
         call reset_buffer(fx%buffer)
      else if (c == ">") then
         fx%state = ERROR
         fx%action = ("Ending tag without starting it!")
      else if  (c == char(10)) then
         fx%state = IN_PCDATA_AT_EOL
         signal = CHUNK_OF_PCDATA
         if (fx%debug) fx%action = ("Resetting PCDATA buffer at newline")
         call add_to_buffer(c,fx%buffer)
         fx%pcdata = fx%buffer
         call reset_buffer(fx%buffer)
      else
         call add_to_buffer(c,fx%buffer)
         if (fx%debug) fx%action = ("Reading chars outside tags")
         !
         ! Check whether we are close to the end of the buffer. 
         ! If so, make a chunk and reset the buffer
         if (c .in. whitespace) then
            if (buffer_nearly_full(fx%buffer)) then
               signal = CHUNK_OF_PCDATA
               if (fx%debug) fx%action = ("Resetting almost full PCDATA buffer")
               fx%pcdata = fx%buffer
               call reset_buffer(fx%buffer)
            endif
         endif
      endif

 case (IN_PCDATA_AT_EOL)
      !
      ! Avoid triggering an extra pcdata event
      !
      if (c == "<") then
         fx%state = START_TAG_MARKER
         if (fx%debug) fx%action = ("No more pcdata after eol-- Starting tag")
      else if (c == ">") then
         fx%state = ERROR
         fx%action = ("Ending tag without starting it!")
      else if  (c == char(10)) then
         fx%state = IN_PCDATA_AT_EOL
         signal = CHUNK_OF_PCDATA
         if (fx%debug) fx%action = ("Resetting PCDATA buffer at repeated newline")
         call add_to_buffer(c,fx%buffer)
         fx%pcdata = fx%buffer
         call reset_buffer(fx%buffer)
      else
         fx%state = IN_PCDATA
         call add_to_buffer(c,fx%buffer)
         if (fx%debug) fx%action = ("Reading chars outside tags")
         !
         ! Check whether we are close to the end of the buffer. 
         ! If so, make a chunk and reset the buffer
         if (c .in. whitespace) then
            if (buffer_nearly_full(fx%buffer)) then
               signal = CHUNK_OF_PCDATA
               if (fx%debug) fx%action = ("Resetting almost full PCDATA buffer")
               fx%pcdata = fx%buffer
               call reset_buffer(fx%buffer)
            endif
         endif
      endif



 case (END_TAG_MARKER)
!
      if (c == "<") then
         fx%state = START_TAG_MARKER
         if (fx%debug) fx%action = ("Starting tag")
      else if (c == ">") then
         fx%state = ERROR
         fx%action = ("Double ending of tag!")
!
!     We should make this whitespace in general (maybe not?
!     how about indentation in text chunks?)
!     See specs.
!
      else if (c == char(10)) then
        ! Ignoring LF after end of tag is probably non standard...

         if (fx%debug) &
            fx%action = ("---------Discarding newline after end of tag")

        !!!  New code for full compliance
        ! fx%state = IN_PCDATA_AT_EOL
        ! call add_to_buffer(c,fx%buffer)
        ! if (fx%debug) &
        !    fx%action = ("Found LF after end of tag. Emitting PCDATA event")
        ! signal = CHUNK_OF_PCDATA
        ! fx%pcdata = fx%buffer
        ! call reset_buffer(fx%buffer)
      else
         fx%state = IN_PCDATA
         call add_to_buffer(c,fx%buffer)
         if (fx%debug) fx%action = ("End of Tag. Starting to read PCDATA")
      endif

 case (ERROR)

      stop "Cannot continue after parsing errors!"

 end select

if (fx%state == ERROR) signal  = EXCEPTION

end subroutine evolve_fsm

end module xmlf90_fsm













