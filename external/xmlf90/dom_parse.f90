module xmlf90_dom_parse

  use xmlf90_dom_types
  use xmlf90_dom_element
  use xmlf90_dom_document
  use xmlf90_dom_node
!  use m_dom_namednodemap
  use xmlf90_dom_debug

  use xmlf90_flib_sax

  implicit none
  
  private

  public :: parsefile

  type(fnode), pointer, private, save  :: main
  type(fnode), pointer, private, save  :: current


CONTAINS

  subroutine begin_element_handler(name,attrs)

    character(len=*),   intent(in) :: name
    type(dictionary_t), intent(in) :: attrs
   
    type(fnode), pointer :: temp
    character(len=400)   :: attr_name, attr_value
    integer              :: status
    integer              :: i

    if (dom_debug) print *, "Adding node for element: ", name

    temp => createElement(name)
    current => appendChild(current,temp)
!
!   Add attributes
!
    ! WORKAROUND due to bug in Intel 11.0
    !do i = 1, len(attrs)
    do i = 1, number_of_entries(attrs)
       call get_name(attrs, i, attr_name, status)
       call get_value(attrs, attr_name, attr_value, status)
       if (dom_debug) print *, "Adding attribute: ", &
                       trim(attr_name), ":",trim(attr_value)
       call setAttribute(current,attr_name,attr_value)
    enddo

  end subroutine begin_element_handler

!---------------------------------------------------------

  subroutine end_element_handler(name)
    character(len=*), intent(in)     :: name

!!AG for IBM    type(fnode), pointer :: np

    if (dom_debug) print *, "End of element: ", name
!!AG for IBM    np => getParentNode(current)
!!AG for IBM    current => np
    current => getParentNode(current)
  end subroutine end_element_handler

!---------------------------------------------------------

  subroutine pcdata_chunk_handler(chunk)
    character(len=*), intent(in) :: chunk

    type(fnode), pointer :: temp, dummy
    
    if (dom_debug) print *, "Got PCDATA: |", chunk, "|"

    temp => createTextNode(chunk)
    dummy => appendChild(current,temp)

  end subroutine pcdata_chunk_handler

!---------------------------------------------------------

  subroutine comment_handler(comment)
    character(len=*), intent(in) :: comment

    type(fnode), pointer :: temp, dummy

    if (dom_debug) print *, "Got COMMENT: |", comment, "|"

    temp => createComment(comment)
    dummy => appendChild(current,temp)

  end subroutine comment_handler
!---------------------------------------------------------
  subroutine cdata_section_handler(chunk)
    character(len=*), intent(in) :: chunk

    type(fnode), pointer :: temp, dummy
    
    if (dom_debug) print *, "Got CDATA_SECTION: |", chunk, "|"

    temp => createCdataSection(chunk)
    dummy => appendChild(current,temp)

  end subroutine cdata_section_handler

!***************************************************
!   PUBLIC PROCEDURES
!***************************************************


  function parsefile(filename, verbose, sax_verbose)

    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: verbose
    logical, intent(in), optional :: sax_verbose

    type(fnode), pointer :: parsefile

    logical :: sax_debug = .false.

    type(xml_t) :: fxml
    integer :: iostat

    if (present(verbose)) then
       dom_debug = verbose
    endif

    if (present(sax_verbose)) then
       sax_debug = sax_verbose
    endif
    
    call open_xmlfile(filename, fxml, iostat)

    if (iostat /= 0) then
       write (*,*) "Cannot open file '" // trim(filename) // "'."
       stop
    endif

    main => createDocumentNode()
    current => main

    call xml_parse(fxml,  &
           begin_element_handler, end_element_handler, pcdata_chunk_handler, &
           comment_handler, cdata_section_handler=cdata_section_handler, &
           verbose = sax_debug)    
    call close_xmlfile(fxml)

    parsefile => main
    if (dom_debug) print *, "Number of allocated nodes: ", getNumberofAllocatedNodes()

  end function parsefile


END MODULE xmlf90_dom_parse
