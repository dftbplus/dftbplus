!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Fills the derived type with the input parameters from an HSD or an XML file.
module parser_setup
  use dftbp_globalenv
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_constants
  use dftbp_typegeometryhsd
  use dftbp_hsdparser, only : dumpHSD, dumpHSDAsXML, getNodeHSDName
  use dftbp_hsdutils
  use dftbp_hsdutils2
  use dftbp_charmanip
  use dftbp_message
  use dftbp_linkedlist
  use dftbp_wrappedintr
  use dftbp_unitconversion
  use dftbp_periodic
  use dftbp_simplealgebra, only: cross3, determinant33
  use dftbp_dispersions
  use dftbp_slakocont
  use dftbp_slakoeqgrid
  use dftbp_repcont
  use dftbp_repspline
  use dftbp_reppoly
  use dftbp_commontypes
  use dftbp_oldskdata
  use dftbp_xmlf90
  use dftbp_wrappedintr
#:if WITH_TRANSPORT
  use libnegf_vars
  use helpsetupgeom
#:endif
  use inputdata_setup
  use inputconversion
  use oldcompat
  implicit none

  private

  public :: parseHsdInput, parserVersion


  ! Default file names

  !> Main HSD input file
  character(len=*), parameter :: hsdInputName = "dftb_in.hsd"

  !> XML input file
  character(len=*), parameter :: xmlInputName = "dftb_in.xml"

  !> Processed HSD input
  character(len=*), parameter :: hsdProcInputName = "dftb_pin.hsd"

  !> Processed  XML input
  character(len=*), parameter :: xmlProcInputName = "dftb_pin.xml"

  !> Tag at the head of the input document tree
  character(len=*), parameter :: rootTag = "dftb_in"


  !> Version of the current parser
  integer, parameter :: parserVersion = 6


  !> Version of the oldest parser for which compatibility is still maintained
  integer, parameter :: minVersion = 1


  !> Container type for parser related flags.
  type TParserFlags

    !> stop after parsing?
    logical :: tStop

    !> Continue despite unprocessed nodes
    logical :: tIgnoreUnprocessed

    !> XML output?
    logical :: tWriteXML

    !> HSD output?
    logical :: tWriteHSD
  end type TParserFlags


contains


  !> Parse input from an HSD/XML file
  subroutine parseHsdInput(input)

    !> Returns initialised input variables on exit
    type(inputData), intent(out) :: input

    type(fnode), pointer :: hsdTree
    type(fnode), pointer :: root, tmp, hamNode, analysisNode, child, dummy
    type(TParserflags) :: parserFlags
    logical :: tHSD, missing

    write(stdOut, "(/, A, /)") "***  Parsing and initializing"

    ! Read in the input
    call readHSDOrXML(hsdInputName, xmlInputName, rootTag, hsdTree, tHSD, &
        &missing)

    !! If input is missing return
    if (missing) then
      call error("No input file found.")
    end if

    write(stdout, '(A,1X,I0,/)') 'Parser version:', parserVersion
    if (tHSD) then
      write(stdout, "(A)") "Interpreting input file '" // hsdInputName // "'"
    else
      write(stdout, "(A)") "Interpreting input file '" // xmlInputName //  "'"
    end if
    write(stdout, "(A)") repeat("-", 80)

    ! Get the root of all evil ;-)
    call getChild(hsdTree, rootTag, root)

    ! Handle parser options
    call getChildValue(root, "ParserOptions", dummy, "", child=child, &
        &list=.true., allowEmptyValue=.true., dummyValue=.true.)
    call readParserOptions(child, root, parserFlags)

    ! Read in the different blocks

    ! Atomic geometry and boundary conditions
    call getChild(root, "Geometry", tmp)
    call readGeometry(tmp, input)

    call getChild(root, "Transport", child, requested=.false.)

  #:if WITH_TRANSPORT

    ! Read in transport and modify geometry if it is only a contact calculation
    if (associated(child)) then
      call readTransportGeometry(child, input%geom, input%transpar)
    else
      input%transpar%ncont=0
      allocate(input%transpar%contacts(0))
      ! set range of atoms in the 'device', as there are no contacts
      input%transpar%idxdevice(1) = 1
      input%transpar%idxdevice(2) = input%geom%nAtom
    end if

  #:else

    if (associated(child)) then
      call detailedError(child, "Program had been compiled without transport enabled")
    end if

  #:endif

    ! input data strucutre has been initialised
    input%tInitialized = .true.

    ! Issue warning about unprocessed nodes
    call warnUnprocessedNodes(root, parserFlags%tIgnoreUnprocessed)

    ! Dump processed tree in HSD and XML format
    if (tIoProc .and. parserFlags%tWriteHSD) then
      call dumpHSD(hsdTree, hsdProcInputName)
      write(stdout, '(/,/,A)') "Processed input in HSD format written to '" &
          &// hsdProcInputName // "'"
    end if
    if (tIoProc .and. parserFlags%tWriteXML) then
      call dumpHSDAsXML(hsdTree, xmlProcInputName)
      write(stdout, '(A,/)') "Processed input in XML format written to '" &
          &// xmlProcInputName // "'"
    end if

    ! Stop, if only parsing is required
    if (parserFlags%tStop) then
      call error("Keyword 'StopAfterParsing' is set to Yes. Stopping.")
    end if

    call destroyNode(hsdTree)
      
    write(stdout,*) 'Geometry processed. Job finished'

  end subroutine parseHsdInput


  !> Read in parser options (options not passed to the main code)
  subroutine readParserOptions(node, root, flags)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Root of the entire tree (in case it needs to be converted, for example because dftbp_of compability
    !> options)
    type(fnode), pointer :: root

    !> Contains parser flags on exit.
    type(TParserFlags), intent(out) :: flags

    integer :: inputVersion
    type(fnode), pointer :: child

    ! Check if input needs compatibility conversion.
    call getChildValue(node, "ParserVersion", inputVersion, parserVersion, &
        &child=child)
    if (inputVersion < 1 .or. inputVersion > parserVersion) then
      call detailedError(child, "Invalid parser version (" // i2c(inputVersion)&
          &// ")")
    elseif (inputVersion < minVersion) then
      call detailedError(child, &
          &"Sorry, no compatibility mode for parser version " &
          &// i2c(inputVersion) // " (too old)")
    elseif (inputVersion /= parserVersion) then
      write(stdout, "(A,I2,A,I2,A)") "***  Converting input from version ", &
          &inputVersion, " to version ", parserVersion, " ..."
      call convertOldHSD(root, inputVersion, parserVersion)
      write(stdout, "(A,/)") "***  Done."
    end if

    call getChildValue(node, "WriteHSDInput", flags%tWriteHSD, .true.)
    call getChildValue(node, "WriteXMLInput", flags%tWriteXML, .false.)
    if (.not. (flags%tWriteHSD .or. flags%tWriteXML)) then
      call detailedWarning(node, &
          &"WriteHSDInput and WriteXMLInput both turned off. You are not&
          & guaranteed" &
          &// newline // &
          &" to able to obtain the same results with a later version of the&
          & code!")
    end if
    call getChildValue(node, "StopAfterParsing", flags%tStop, .false.)

    call getChildValue(node, "IgnoreUnprocessedNodes", &
        &flags%tIgnoreUnprocessed, .false.)

  end subroutine readParserOptions


  !> Read in Geometry
  subroutine readGeometry(node, input)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Input structure to be filled
    type(inputData), intent(inout) :: input

    type(fnode), pointer :: value1, child
    type(string) :: buffer

    call getChildValue(node, "", value1, child=child)
    call getNodeName(value1, buffer)
    select case (char(buffer))
    case ("genformat")
      call readTGeometryGen(value1, input%geom)
    case default
      call setUnprocessed(value1)
      call readTGeometryHSD(child, input%geom)
    end select

  end subroutine readGeometry


#:if WITH_TRANSPORT
  !> Read geometry information for transport calculation
  subroutine readTransportGeometry(root, geom, transpar)

    !> Root node containing the current block
    type(fnode), pointer :: root

    !> geometry of the system, which may be modified for some types of calculation
    type(TGeometry), intent(inout) :: geom

    !> Parameters of the transport calculation
    type(TTransPar), intent(inout) :: transpar

    type(fnode), pointer :: pGeom, pDevice, pNode, pTask, pTaskType
    type(string) :: buffer, modif
    type(fnode), pointer :: pTmp, field
    type(fnodeList), pointer :: pNodeList
    integer :: ii, contact
    real(dp) :: acc, contactRange(2), lateralContactSeparation, plCutoff
    type(listInt) :: li
    type(WrappedInt1), allocatable :: iAtInRegion(:)
    real(dp), allocatable :: contVec(:,:)
    integer, allocatable :: nPLs(:)

    transpar%defined = .true.
    transpar%tPeriodic1D = .not. geom%tPeriodic
    !call getChild(pDevice, "ContactPLs", pTmp, requested=.false.)
    !if (associated(pTmp)) then
    !  call init(li)
    !  call getChildValue(pTmp, "", li)
    !  allocate(transpar%cblk(len(li)))
    !  call asArray(li,transpar%cblk)
    !  call destruct(li)
    !end if

    !! Note: we parse first the task because dftbp_we need to know it to define the
    !! mandatory contact entries. On the other hand we need to wait that
    !! contacts are parsed to resolve the name of the contact for task =
    !! contacthamiltonian
    call getChildValue(root, "Task", pTaskType, child=pTask, default='uploadcontacts')
    call getNodeName(pTaskType, buffer)

    if (char(buffer).ne."setupgeometry") then
      call getChild(root, "Device", pDevice)
      call getChildValue(pDevice, "AtomRange", transpar%idxdevice)
      !call getChild(pDevice, "FirstLayerAtoms", pTmp, requested=.false.)
      !call readFirstLayerAtoms(pTmp, transpar%PL, transpar%nPLs, transpar%idxdevice)
      !if (.not.associated(pTmp)) then
      !  call setChildValue(pDevice, "FirstLayerAtoms", transpar%PL)
      !end if
    end if
    
    call getChildren(root, "Contact", pNodeList)
    transpar%ncont = getLength(pNodeList)
    if (transpar%ncont < 2) then
      call detailedError(root, "At least two contacts must be defined")
    end if
    allocate(transpar%contacts(transpar%ncont))

    select case (char(buffer))
    case ("setupgeometry")
      
      call readContacts(pNodeList, transpar%contacts, geom, char(buffer), iAtInRegion, &
            contVec, nPLs)
      call getChildValue(pTask, "PLCutoff", plCutoff, 10.0_dp, modifier=modif, child=field)
      call convertByMul(char(modif), lengthUnits, field, plCutoff)
      call setupGeometry(geom, iAtInRegion, contVec, plCutoff, nPLs)

    case default

      call getNodeHSDName(pTaskType, buffer)
      call detailedError(pTask, "Invalid task '" // char(buffer) // "'")

   end select

   call destroyNodeList(pNodeList)

  end subroutine readTransportGeometry


  !> Read bias information, used in Analysis and Green's function eigensolver
  subroutine readContacts(pNodeList, contacts, geom, task, iAtInRegion, contVec, nPLs)
    type(fnodeList), pointer :: pNodeList
    type(ContactInfo), allocatable, dimension(:), intent(inout) :: contacts
    type(TGeometry), intent(in) :: geom
    character(*), intent(in) :: task
    type(WrappedInt1), allocatable, intent(out), optional :: iAtInRegion(:)
    real(dp), intent(out), allocatable, optional :: contVec(:,:)
    integer, intent(out), allocatable, optional :: nPLs(:)

    real(dp) :: contactLayerTol, vec(3)
    integer :: ii, jj
    type(fnode), pointer :: field, pNode, pTmp, pWide
    type(string) :: buffer, modif
    type(listReal) :: fermiBuffer, vecBuffer
    integer, allocatable :: tmpI1(:)

    if (present(iAtInRegion)) then
      allocate(iAtInRegion(size(contacts)+1))
    end if
    if (present(contVec)) then
      allocate(contVec(4,size(contacts)))
    end if 
    if (present(nPLs)) then
      allocate(nPLs(size(contacts)))
    end if

    do ii = 1, size(contacts)

      contacts(ii)%wideBand = .false.
      contacts(ii)%wideBandDos = 0.0_dp

      call getItem1(pNodeList, ii, pNode)
      call getChildValue(pNode, "Id", buffer, child=pTmp)
      buffer = tolower(trim(unquote(char(buffer))))
      if (len(buffer) > mc) then
        call detailedError(pTmp, "Contact id may not be longer than " // i2c(mc) // " characters.")
      end if
      contacts(ii)%name = char(buffer)
      if (any(contacts(1:ii-1)%name == contacts(ii)%name)) then
        call detailedError(pTmp, "Contact id '" // trim(contacts(ii)%name) &
            &//  "' already in use")
      end if

      call getChildValue(pNode, "PLShiftTolerance", contactLayerTol, 1e-5_dp, modifier=modif,&
          & child=field)
      call convertByMul(char(modif), lengthUnits, field, contactLayerTol)

      if (task .eq. "setupgeometry") then
        call getChildValue(pNode, "NumPLsDefined", nPLs(ii), 2)    
        call getChildValue(pNode, "Atoms", buffer, child=pTmp, modifier=modif, multiple=.true.)
        call convAtomRangeToInt(char(buffer), geom%speciesNames, geom%species, pTmp, &
             iAtInRegion(ii)%data, ishift=char_to_int(char(modif)))
        call init(vecBuffer)
        call getChildValue(pNode, "ContactVector", vecBuffer, modifier=modif)
        if (len(vecBuffer).eq.3) then
           call asArray(vecBuffer, vec)
           call convertByMul(char(modif), lengthUnits, pNode, vec)
           contVec(1:3,ii) = vec
           contVec(4,ii) = contactLayerTol
           call destruct(vecBuffer)
        else
           call error("ContactVector must define three entries")
        end if   
      else
        call error("Invalid task for setpugeometry tool")
      end if

    end do

    contains

      function char_to_int(chr) result(ind)
        character(*), intent(in) :: chr
        integer :: ind
        if (trim(chr) .eq. "") then
          ind = 0
          return
        end if
        if (verify(chr,"+-0123456789") .ne. 0) then
          call error("Modifier in Atoms should be an integer number")   
        end if  
        read(chr,*) ind
      end function char_to_int

  end subroutine readContacts

#:endif

end module parser_setup
