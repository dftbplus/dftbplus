!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the HSD (Human readable Structured Data) parser.
!>
!> The HSD format is a more or less user friendly input format, which can be easily converted to a
!> simplified XML format. The parser returns a DOM-tree, which can be further processed. The
!> returned tree contains also information about the original name and position of the keywords in
!> the original HSD format, in order to enable user friendly error messages if inconsistent data are
!> detected during the processing of the DOM-tree.
!>
!> For the specification of the HSD format see the sample input
module dftbp_hsdparser
  use dftbp_assert
  use dftbp_message
  use dftbp_charmanip
  use dftbp_xmlutils
  use dftbp_xmlf90
  implicit none
  private


  !> Wrapper around the parsing function
  interface parseHSD
    module procedure parseHSD_file
    module procedure parseHSD_opened
  end interface parseHSD


  !> Wrapper around the HSD dumping
  interface dumpHSD
    module procedure dumpHSD_file
    module procedure dumpHSD_opened
  end interface dumpHSD

  ! Main token separator characters

  !> number of separator character strings
  integer, parameter :: nSeparator = 7

  !> XML includer
  character(len=*), parameter :: sIncludeXML = "<<!"

  !> include parsed material
  character(len=*), parameter :: sIncludeParsed = "<<+"

  !> include unparsed material
  character(len=*), parameter :: sIncludeUnparsed = "<<<"

  !> open for a single thing
  character(len=*), parameter :: sSingleOpen = "=  "

  !> open region
  character(len=*), parameter :: sOpen = "{  "

  !> close region
  character(len=*), parameter :: sClose = "}  "

  !> close region
  character(len=*), parameter :: sSingleClose = ";  "

  !> Collect together as an array
  character(len=*), parameter :: separators(nSeparator) = &
      &(/ sIncludeXML, sIncludeParsed, sIncludeUnparsed, &
      &sSingleOpen, sOpen, sClose, sSingleClose /)

  ! Other parsed characters

  !> open modifier
  character(len=*), parameter :: sModifierOpen = "["

  !> close modifier
  character(len=*), parameter :: sModifierClose = "]"

  !> comment mark
  character(len=*), parameter :: sComment = "#"

  ! Extension related stuff

  !> number of parser tag extensions
  integer, parameter :: nExtension = 5

  !> Extend with things, or halt
  character(len=*), parameter :: sExtendIfPresentOrDie = "+"

  !> Optionally extend if present
  character(len=*), parameter :: sExtendIfPresent = "?"

  !> Extend with 0 or more instances
  character(len=*), parameter :: sExtendIfPresentOrCreate = "*"

  !> create if missing
  character(len=*), parameter :: sCreateIfNotPresent = "/"

  !> replace if present, or create
  character(len=*), parameter :: sReplaceIfPresentOrCreate = "!"

  !> Collect together as an array
  character(len=*), parameter :: extensions(nExtension) = &
      &(/ sExtendIfPresentOrDie, sExtendIfPresent, sExtendIfPresentOrCreate, &
      &sCreateIfNotPresent, sReplaceIfPresentOrCreate /)

  ! Name and file descriptors for standard input/output

  !> Forbidden (even if quoted) characters in the iput
  character(len=*), parameter :: forbiddenChars = "<>"

  ! Attribute names

  !> Start of attribute
  character(len=*), parameter :: attrStart = "start"

  !> end of attribute
  character(len=*), parameter :: attrEnd = "end"

  !> file name
  character(len=*), parameter :: attrFile = "file"

  !> attribte name
  character(len=*), parameter :: attrName = "name"

  !> modifier label
  character(len=*), parameter :: attrModifier = "m"

  !> list label
  character(len=*), parameter :: attrList = "l"


  !> Length of a parsed line
  integer, parameter :: lc = 1024


  !> Maximal record lenght on output in characters (bytes).
  !> If text nodes bigger than that occur runtime error can be expected.
  integer, parameter :: MAXRECL = 1024 * 1024


  !> Name of the root tag
  character(len=lc) :: rootName


  !> Pointer to the top of the currently processed document
  !> (modified only in parseHSD and replaceTreeFromFile)
  type(fnode), pointer :: myDoc


  !> Format of the input line
  character(len=lc) :: lineFormat = ""

  public :: parseHSD, dumpHSD, newline
  public :: getNodeHSDName, getHSDPath
  public :: attrStart, attrEnd, attrFile, attrName, attrModifier, attrList

contains

  !> Parser HSD format from a file
  subroutine parseHSD_file(initRootName, file, xmlDoc)

    !> Name of the root tag, which should contain the parsed tree
    character(len=*), intent(in) :: initRootName

    !> Name of the file (used in error messages)
    character(len=*), intent(in) :: file

    !> DOM-tree of the parsed input on exit
    type(fnode), pointer :: xmlDoc

    integer :: fd
    integer :: iostat

    open(newunit=fd, file=file, iostat=iostat, status='old', action='read', recl=lc)
    if (iostat /= 0) then
      call parsingError("Error in opening file '" // trim(file) //"'.", file, -1)
    end if
    call parseHSD_opened(initRootName, fd, file, xmlDoc)
    close(fd, iostat=iostat)

  end subroutine parseHSD_file


  !> Parses HSD format from an already opened file
  subroutine parseHSD_opened(initRootName, fd, file, xmlDoc)

    !> Name of the root tag, which should contain the parsed tree
    character(len=*), intent(in) :: initRootName

    !> File descriptor of the open file containing the input
    integer, intent(in) :: fd

    !> Name of the file (used in error messages)
    character(len=*), intent(in), optional :: file

    !> DOM-tree of the parsed input on exit
    type(fnode), pointer :: xmlDoc

    type(fnode), pointer :: rootNode, dummy
    logical :: tFinished
    integer :: curLine
    character(len=lc) :: residual, curFile

    if (present(file)) then
      curFile = file
    else
      curFile = "???"
    end if

    if (len_trim(lineFormat) == 0) then
      lineFormat = "(A" // i2c(lc) // ")"
    end if
    rootName = tolower(initRootName(:min(lc, len(initRootName))))
    myDoc => createDocumentNode()
    rootNode => createElement(trim(rootName))
    dummy => appendChild(myDoc, rootNode)
    curLine = 0
    residual = ""
    tFinished = parse_recursive(rootNode, 0, residual, .false., fd, curFile, &
        &0, curLine, &
        &(/ .true., .true., .true., .true., .true., .true., .true. /), .false.)
    xmlDoc => myDoc
    myDoc => null()

  end subroutine parseHSD_opened


  !> Recursive parsing function for the HSD parser making the actual work
  recursive function parse_recursive(curNode, depth, residual, tRightValue, fd, curFile, fileDepth,&
      & curLine, parsedTypes, tNew) result (tFinished)

    !> Node which should contain parsed input
    type(fnode), pointer :: curNode

    !> Number of open blocks/assignments.
    integer, intent(in) :: depth

    !> Unparsed text from the previous line
    character(len=lc), intent(inout) :: residual

    !> Is next parsed token a right value of an assignment?
    logical, intent(in) :: tRightValue

    !> File descriptor of the input
    integer, intent(in) :: fd

    !> Name of the current input file
    character(len=lc), intent(in) :: curFile

    !> Number of open files
    integer, intent(in) :: fileDepth

    !> Number of current line in the current file
    integer, intent(inout) :: curLine

    !> True for those separators, which should be parsed
    logical, intent(in) :: parsedTypes(nSeparator)

    !> True, if parsing is done
    logical, intent(in) :: tNew

    logical :: tFinished

    character(len=lc) :: strLine, word

    type(fnode), pointer :: childNode, dummy
    type(string) :: buffer
    integer :: newFile
    integer :: iostat
    integer :: iType, sepPos
    integer :: newCurLine
    logical :: tTagClosed, tNewNodeCreated
    logical :: newParsedTypes(nSeparator)
    integer :: nTextLine
    integer :: iTmp
    integer :: nodetype

    tTagClosed = .false.
    tFinished = .false.
    tNewNodeCreated = .false.
    nTextLine = 0
    nodetype = 0

    lpMain: do while ((.not. tTagClosed) .and. (.not. tFinished))

      !! Read in next line or process residual from last line.
      if (len_trim(residual) /= 0) then
        strLine = adjustl(residual)
        residual = ""
      else
        read (fd, trim(lineFormat), iostat=iostat) strLine
        curLine = curLine + 1
        call convertWhitespaces(strLine)
        strLine = adjustl(strLine)
        !! If reading error (e.g. EOF) -> close current scope
        if (iostat /= 0) then
          tTagClosed = .true.
          if (depth /= 0) then
            call getAttribute(curNode, attrStart, buffer)
            call parsingError("Unexpected end of input (probably open node at &
                &line " // char(buffer) // " or after).", curFile, curLine)
          end if
          !! If outermost file, we are ready
          if (fileDepth == 0) then
            tFinished = .true.
          end if
          exit
        end if
      end if

      !! Remove comments
      iTmp = unquotedIndex(strLine, sComment)
      if (iTmp /= 0) then
        strLine = strLine(:iTmp-1)
      end if

      !! Get first occurance of any separator
      call getFirstOccurance(strLine, separators, parsedTypes, iType, sepPos)

      !! Handle various closing operators.
      select case (iType)
      case (6)
        !! Block closing char on level zero is invalid
        if (depth == 0) then
          call parsingError("Invalid block closing sign.", curFile, curLine)
        end if
        !! If block closing char is not first char of the line, text before it
        !! will be appended as text, and residual line reparsed in next cycle
        if (sepPos == 1) then
          tTagClosed = .true.
          residual = strLine(sepPos+1:)
        else
          iType = 0
          residual = strLine(sepPos:)
          strLine = strLine(:sepPos-1)
        end if
      case (7)
        if (tRightValue) then
          !! Single assignment is terminated by a semi-colon. Text after
          !! semicolon must be reparsed in next cycle.
          iType = 0
          residual = strLine(sepPos+1:)
          strLine = strLine(:sepPos-1)
          tTagClosed = .true.
        else
          call parsingError("Invalid assignment separator", curFile, curLine)
        end if
      end select

      !! Ignore empty lines
      if (len_trim(strLine) == 0) then
        cycle lpMain
      end if

      !! Check for forbidden characters in the current line
      if (.not. (iType == 1 .or. iType == 2 .or. iType == 3)) then
        call checkForbiddenChars(strLine, curFile, curLine)
      end if

      !! Process non-closing separators
      select case (iType)
      case(0)
        if (nodetype > 0) then
          call parsingError("Node already contains subnodes, no text content&
              & allowed any more", curFile, curLine)
        end if
        !! No separator found -> Add entire line as text
        !! If current node already contains text, prepend newline before
        !! appending. (Otherwise newlines in the text would get lost.)
        if (associated(curNode)) then
          nTextLine = nTextLine + 1
          if (nTextLine > 1) then
            childNode => createTextNode(newline // trim(strLine))
          else
            childNode => createTextNode(trim(strLine))
          end if
          dummy => appendChild(curNode, childNode)
        end if
        nodetype = -1

      case(1)
        !! XML inclusion
        call error("Mixed XML input in HSD input no longer supported")

      case(2, 3)
        !! File inclusion operator -> append content of new file to current node
        if (associated(curNode)) then
          if (sepPos /= 1) then
            call parsingError("Invalid character before file inclusion &
                &operator", curFile, curLine)
          end if
          strLine = adjustl(unquote(strLine))
          if (iType == 2) then
            word = adjustl(strLine(len(sIncludeParsed)+1:len_trim(strLine)))
          else
            word = adjustl(strLine(len(sIncludeUnparsed)+1:len_trim(strLine)))
          end if

          if (len_trim(word) == 0) then
            call parsingError("No file name specified after the inclusion &
                &operator.", curFile, curLine)
          end if

          open(newunit=newFile, file=trim(word), status='old', action='read', &
              &iostat=iostat, recl=lc)
          if (iostat /= 0) then
            call parsingError("Error in opening file '" // trim(word) // &
                &"'.", curFile, curLine)
          end if
          strLine = ""
          newCurLine = 0
          if (iType == 2) then
            !! Everything is parsed
            newParsedTypes = (/ .true., .true., .true., .true., .true., &
                &.true., .true. /)
          else
            !! Nothing is parsed
            newParsedTypes = (/ .false., .false., .false., .false., .false., &
                &.false., .false. /)
          end if
          tFinished = parse_recursive(curNode, 0, strLine, .false., newFile, &
              &word, fileDepth + 1, newCurLine, newParsedTypes, .false.)
          close(newFile, iostat=iostat)
        end if

      case(4)
        !! Assignment
        if (nodetype < 0) then
          call parsingError("Node already contains free text, no child nodes&
              & allowed any more", curFile, curLine)
        end if
        word = adjustl(strLine(:sepPos-1))
        strLine = adjustl(strLine(sepPos+1:))
        if (len_trim(word) == 0) then
          call parsingError("Missing field name on the left of the &
              &assignment.", curFile, curLine)
        elseif (len_trim(strLine) == 0) then
          call parsingError("Missing value on the right of the assignment.", &
              &curFile, curLine)
        elseif (nTextLine > 0) then
          call parsingError("Unparsed text before current node", curFile, &
              &curLine)
        end if
        if (associated(curNode)) then
          childNode => createChildNode(curNode, word, curLine, curFile)
          tNewNodeCreated = .true.
        else
          childNode => null()
        end if
        !! Only block opening/closing sign and single child separator are parsed
        newParsedTypes = (/ .false., .false., .false., .false., .true., &
            &.true., .true. /)
        tFinished = parse_recursive(childNode, depth+1, strLine, .true., fd, &
            &curFile, fileDepth, curLine, newParsedTypes, tNewNodeCreated)
        residual = strLine
        nodetype = 1

      case(5)
        if (nodetype < 0) then
          call parsingError("Node already contains free text, no child nodes&
              & allowed any more", curFile, curLine)
        end if
        !! Block opening sign
        word = adjustl(strLine(:sepPos-1))
        strLine = adjustl(strLine(sepPos+1:))
        if (nTextLine > 0) then
          call parsingError("Unparsed text before current node", curFile, &
              &curLine)
        end if
        if (associated(curNode)) then
          ! Currently node without name is allowed to support "= {" construct
          ! Should be turned to parsing error to deprecate that construct.
          if (len_trim(word) == 0) then
            childNode => curNode
            call setAttribute(curNode, attrList, "")
            !call parsingError("Node without name not allowed.", curFile,&
            !    & curLine)
          else
            childNode => createChildNode(curNode, word, curLine, curFile)
            tNewNodeCreated = .true.
          end if
        else
          childNode => null()
        end if
        newParsedTypes = (/ .true., .true., .true., .true., .true., .true., &
            &.true. /)
        tFinished = parse_recursive(childNode, depth+1, strLine, .false., &
            &fd, curFile, fileDepth, curLine, newParsedTypes, tNewNodeCreated)
        residual = strLine
        nodetype = 1

      end select

      tTagClosed = tTagClosed .or. tRightValue

    end do lpMain

    !! Set end attribute on tag end and normalise text nodes
    if (tTagClosed .and. associated(curNode)) then
      if (tNew) then
        call setAttribute(curNode, attrEnd, i2c(curLine))
      end if
      if (nTextLine > 1) then
        call normalize(curNode)
      end if
    end if

  end function parse_recursive


  !> Creates a child node with attributes related to the HSD input.
  function createChildNode(parentNode, childName, curLine, file) result(newChild)

    !> Parent node containing of the child to be created
    type(fnode), pointer :: parentNode

    !> Name of the new child
    character(len=lc), intent(in) :: childName

    !> Number of the current line
    integer, intent(in) :: curLine

    !> Name of the current file
    character(len=lc), intent(in) :: file

    !> Pointer to the new (appended) child node
    type(fnode), pointer :: newChild

    type(fnode), pointer :: dummy, sameChild
    character(len=lc) :: lowerName, truncName, modifier
    logical :: tModifier, tCreate
    integer :: pos1, pos2, iType
    integer :: ii

    truncName = childName
    lowerName = tolower(childName)

    !! Look for any extension operator
    iType = 0
    pos1 = 0
    do ii = 1, nExtension
      pos1 = len_trim(extensions(ii))
      if (lowerName(:pos1) == trim(extensions(ii))) then
        iType = ii
        exit
      end if
    end do

    !! Cut extension operator from the field name
    if (iType /= 0) then
      lowerName = lowerName(pos1+1:)
      truncName = truncName(pos1+1:)
    end if

    !! Look for modifier after field name
    tModifier = .false.
    pos1 = index(lowerName, sModifierOpen)
    pos2 = index(lowerName, sModifierClose)
    if (pos1 == 0) then
      if (pos2 /= 0) then
        call parsingError("Unbalanced modifier opening sign.", file, curLine)
      end if
    else
      if (pos2 /= len_trim(lowerName)) then
        call parsingError("Invalid character(s) after modifier closing sign.",&
            &file, curLine)
      end if
      !! Remove modifier from field name
      modifier = adjustl(truncName(pos1+1:pos2-1))
      lowerName = adjustl(lowerName(:pos1-1))
      truncName = adjustl(truncName(:pos1-1))
      tModifier = .true.
    end if

    !! Check if field name is nonempty
    if (len_trim(lowerName) == 0) then
      call parsingError("Missing field name", file, curLine)
    end if

    !! Create child according extension operator
    tCreate = .false.
    if (iType == 0) then
      !! Create and append new node
      tCreate = .true.
    else
      !! Look for already present node with the same name
      sameChild => getFirstChildByName(parentNode, trim(lowerName))
      if (associated(sameChild)) then
        !! We have found a block with the same name
        if (iType == 4) then
          newChild => null()
          return
        elseif (iType == 5) then
          dummy => removeChild(parentNode, sameChild)
          call destroyNode(sameChild)
          tCreate = .true.
        else
          newChild => sameChild
        end if
      else
        !! We did not found a child with the same name
        select case (iType)
        case(1)
          call parsingError("Containing block does not contain a(n) '" &
              &// trim(truncName) // "' block yet.", file, curLine)
        case(2)
          newChild => null()
          return
        case(3,4,5)
          tCreate = .true.
        end select
      end if
    end if

    !! Create and append the node
    if (tCreate) then
      newChild => createElement(trim(lowerName))
      dummy => appendChild(parentNode, newChild)
    end if

    !! Set useful attributes
    call setAttribute(newChild, attrStart, i2c(curLine))
    call setAttribute(newChild, attrName, trim(truncName))
    call setAttribute(newChild, attrFile, trim(file))
    if (tModifier) then
      call setAttribute(newChild, attrModifier, trim(modifier))
    end if

  end function createChildNode


  !> Checks for forbidden characters and issue error message, if any found.
  subroutine checkForbiddenChars(str, curFile, curLine)

    !> String to investigate
    character(len=*), intent(in) :: str

    !> Name of the current file
    character(len=*), intent(in) :: curFile

    !> Number of the current line in the current file
    integer, intent(in) :: curLine

    if (scan(str, forbiddenChars) /= 0) then
      call parsingError("Invalid character(s).", curFile, curLine)
    end if

  end subroutine checkForbiddenChars


  !> Issues a parsing error message containing file name and line number.
  subroutine parsingError(message, file, line)

    !> Parsing error message
    character(len=*), intent(in) :: message

    !> Name of the current file
    character(*), intent(in) :: file

    !> Number of current line
    integer, intent(in) :: line

    character(len=lc) :: msgArray(2)

    !! Watch out to trunk away enough from the file name to prevent overflow
    if (len_trim(file) > lc - 40) then
      write (msgArray(1), 9991) trim(file(1:lc-40)), line
    else
      write (msgArray(1), 9991) trim(file), line
    end if
9991 format("HSD parser error: File '",A,"', Line",I5,".")
    write (msgArray(2), "(A)") trim(message(:min(lc, len(message))))
    call error(msgArray)

  end subroutine parsingError


  !> Replaces the tree
  !>
  !> The solution with access to a global module variable is not very
  !> elegant, but it saves the deep cloning of the parsed document.
  subroutine replaceTreeFromFile(curNode, file)

    !> Node, which should contained the parsed children from file
    type(fnode), pointer :: curNode

    !> File to parse
    character(len=*), intent(in) :: file

    type(fnode), pointer :: newDoc, rootNode

    newDoc => parsefile(file)
    call removeSpace(newDoc)
    call normalize(newDoc)
    rootNode => getLastChildByName(newDoc, trim(rootName))
    if (.not. associated(rootNode)) then
      call parsingError("File '" // file // "' does not contain '" // trim(rootName) // "' node.",&
          & file, -1)
    else
      call destroyNode(myDoc)
      myDoc => newDoc
      curNode => rootNode
    end if

  end subroutine replaceTreeFromFile


  !> Dumps a HSD tree in a file.
  subroutine dumpHSD_file(myDoc, file, subnode)

    !> The DOM tree
    type(fnode), pointer :: myDoc

    !> Name of the file
    character(len=*), intent(in) :: file

    !> Whether passed node is an arbitrary node within the tree (or the tree top node otherwise).
    !> Default: .false.
    logical, optional, intent(in) :: subnode

    integer :: fd
    integer :: iostat
    character(len=lc) :: fileName

    open(newunit=fd, file=file, iostat=iostat, status='replace', action='write', recl=MAXRECL)
    if (iostat /= 0) then
      fileName = file
      call parsingError("Error in opening file for the HSD output.", fileName, -1)
    end if
    call dumpHSD_opened(myDoc, fd, subnode)
    close(fd)

  end subroutine dumpHSD_file


  !> Dumps a DOM-tree representing a HSD input in HSD format to an opened file.
  subroutine dumpHSD_opened(myDoc, fd, subnode)

    !> The DOM tree
    type(fnode), pointer :: myDoc

    !> File descriptor for an open file where output should go.
    integer, intent(in) :: fd

    !> Whether passed node is an arbitrary node within the tree (or the tree top node otherwise).
    !> Default: .false.
    logical, optional, intent(in) :: subnode

    type(fnode), pointer :: rootNode
    type(fnode), pointer :: child
    type(string) :: buffer
    logical :: subnode_

    if (present(subnode)) then
      subnode_ = subnode
    else
      subnode_ = .false.
    end if
    if (subnode_) then
      rootNode => myDoc
    else
      rootNode => getFirstChild(myDoc)
    end if
    if (.not. associated(rootNode)) then
      return
    end if
    child => getFirstChild(rootNode)
    do while (associated(child))
      call dumpHSD_recursive(child, 0, fd, .false., buffer)
      child => getNextSibling(child)
    end do

  end subroutine dumpHSD_opened


  !> Recursive workhorse for the dumpHSD routine.
  recursive subroutine dumpHSD_recursive(node, indent, fd, tRightValue, buffer)

    !> Node to dump
    type(fnode),      pointer :: node

    !> Current indentation level
    integer, intent(in) :: indent

    !> File descriptor for an open file where output should go
    integer, intent(in) :: fd

    !> Is current node the right hand side of an assignment?
    logical,          intent(in) :: tRightValue

    !> Buffer for storing temporary strings
    type(string),     intent(inout) :: buffer

    type(fnode), pointer :: child, attr
    logical :: tOpenBlock

    !! Text nodes are printed including trailing newline. No further processing
    if (getNodeType(node) == TEXT_NODE) then
      call getNodeValue(node, buffer)
      write (fd, "(A)", advance="yes") trim2(char(buffer))
      return
    end if

    !! If left value -> indent
    if (.not. tRightValue) then
      write (fd, "(A)", advance="no") repeat(" ", indent)
    end if

    !! Get and print tag name
    attr => getAttributeNode(node, attrName)
    if (associated(attr)) then
      call getNodeValue(attr, buffer)
    else
      call getNodeName(node, buffer)
    end if
    write (fd, "(A)", advance="no") char(buffer)

    !! Get and print modifier
    attr => getAttributeNode(node, attrModifier)
    if (associated(attr)) then
      call getNodeValue(attr, buffer)
      if (len(buffer) > 0) then
        write (fd, "(' ',A,A,A)", advance="no") trim(sModifierOpen), &
            &char(buffer), trim(sModifierClose)
      end if
    end if

    child => getFirstChild(node)
    if (tRightValue) then
      if (associated(child)) then
        !! Children exist -> open block
        write (fd, "(' ',A)", advance="yes") trim(sOpen)
        tOpenBlock = .true.
      else
        !! No children -> Open and close block immediately and return
        write (fd, "(' ',A,A)", advance="yes") trim(sOpen), trim(sClose)
        return
      end if
    else
      !! We are the left hand side of an assignment -> write assignment sign
      write (fd, "(' ',A,' ')", advance="no") trim(sSingleOpen)
      attr => getAttributeNode(node, attrList)
      if (associated(child)) then
        tOpenBlock = .false.
        if (associated(attr) .or. associated(getNextSibling(child))) then
          !! RHS has many children or is signed as block -> open a block
          tOpenBlock = .true.
        elseif (getNodeType(child) == TEXT_NODE) then
          !! RHS's child is text -> open block if it contains several tokens
          call getNodeValue(child, buffer)
          tOpenBlock = (unquotedScan(trim2(char(buffer)), whiteSpaces) /= 0)
        end if
        if (tOpenBlock) then
          write (fd, "(A)", advance="yes") trim(sOpen)
        end if
      else
        !! No child, write an empty block and return
        write (fd, "(A,A)", advance="yes") trim(sOpen), trim(sClose)
        return
      end if
    end if

    !! Process children
    do while (associated(child))
      if (tOpenBlock) then
        call dumpHSD_recursive(child, indent+2, fd, .false., buffer)
      else
        call dumpHSD_recursive(child, indent, fd, .true., buffer)
      end if
      child => getNextSibling(child)
    end do

    !! Dump closing sign
    if (tOpenBlock) then
      write (fd, "(A)", advance="no") repeat(" ", indent)
      write (fd, "(A)", advance="yes") trim(sClose)
    end if

  end subroutine dumpHSD_recursive


  !> Returns the name of a node, if present in pretty printing format.
  subroutine getNodeHSDName(node, name)

    !> Node to investigate.
    type(fnode), pointer :: node

    !> name of the node on return
    type(string), intent(inout) :: name

    call getAttribute(node, attrName, name)
    if (len(name) == 0) then
      call getNodeName(node, name)
    end if

  end subroutine getNodeHSDName


  !> Returns the path of a node, if possible in pretty printing format.
  subroutine getHSDPath(node, path, excludeRoot)

    !> Node to investigate
    type(fnode), pointer :: node

    !> String containing the path on return
    type(string), intent(inout) :: path

    !> If root node should be excluded
    logical, intent(in), optional :: excludeRoot

    type(fnode), pointer :: parent
    type(string) :: buffer
    logical :: inclRoot

    if (present(excludeRoot)) then
      inclRoot = .not. excludeRoot
    else
      inclRoot = .true.
    end if

    call getNodeHSDName(node, path)
    parent => getParentNode(node)
    if (associated(parent)) then
      call getHSDPath_recursive(parent, path, inclRoot, buffer)
    elseif (.not. inclRoot) then
      path = ""
    end if

  end subroutine getHSDPath


  !> Workhorse for the getHSDPath routine
  recursive subroutine getHSDPath_recursive(node, path, inclRoot, buffer)

    !> Node to look for
    type(fnode), pointer :: node

    !> String containing the path until now
    type(string), intent(inout) :: path

    !> If document root should be included
    logical, intent(in) :: inclRoot

    !> Buffer for strings (to avoid destruction at every call)
    type(string), intent(inout) :: buffer

    type(fnode), pointer :: parent

    parent => getParentNode(node)
    if (associated(parent) .or. inclRoot) then
      call prepend_to_string(path, "/")
      call getNodeHSDName(node, buffer)
      call prepend_to_string(path, buffer)
      if (associated(parent)) then
        call getHSDPath_recursive(parent, path, inclRoot, buffer)
      end if
    end if

  end subroutine getHSDPath_recursive

end module dftbp_hsdparser
