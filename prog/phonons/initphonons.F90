!!* Contains the routines for initialising phonons.
module InitProgram
#include "allocate.h"
#include "assert.h"
  use mpiglobal
  use HSDParser, only : parseHSD, dumpHSD, dumpHSDAsXML
  use TokenReader
  use XMLUtils
  use HSDUtils
  use HSDUtils2
  use flib_dom
  use LinkedList
  use CharManip
  use Accuracy
  use periodic
  use Constants
  use TypeGeometryHSD
  use Message
  use FileId
  use UnitConversion
  use StringList
  use OldSKData
  use libnegf_vars
  use WrappedIntrinsics
  implicit none

  private
  save

  character(len=*), parameter :: version =  "0.01"
  character(len=*), parameter :: rootTag = "phonons"
  character(len=*), parameter :: hsdInput = "phonons_in.hsd"
  character(len=*), parameter :: hsdParsedInput = "phonons_pin.hsd"
  character(len=*), parameter :: xmlInput = "phonons_in.xml"
  character(len=*), parameter :: xmlParsedInput = "phonons_pin.xml" 

  public :: initProgramVariables, destructProgramVariables
  public :: TPdos

  type TPdos 
    type(WrappedInt1), allocatable :: iAtInRegion(:)
    character(lc), allocatable :: regionLabels(:)
  end type TPdos

  !! Variables from detailed.xml
  integer, public :: identity                  ! Identity of the run
  type(TGeometry), public :: geo               ! Geometry
  type(TPdos), public :: pdos                  ! Projected dos infos
  type(TTransPar), public :: transpar          ! Contains transport parameters
  type(TGDFTBtundos), public :: tundos         ! Contains transport parameters
  type(TGDFTBStructure), public :: str         ! Structure input

  !! Variables from the Option block
  logical, public :: tVerbose          ! If program should be verbose

  real(dp), allocatable, public :: atomicMasses(:)
  real(dp), allocatable, public :: dynMatrix(:,:)
  integer, allocatable, public :: iMovedAtoms(:)
  integer, public :: nMovedAtom
  real(dp), allocatable, public :: KPoint(:,:), KWeight(:)
  integer, public :: nKPoints,  nAtomUnitCell
  integer, pointer, public  :: Img2CentCell(:) !* nr. of original atom in centre
  type(TNeighborList), public :: neighborList
  integer, allocatable, target, public :: nNeighbor(:)    !* nr. of neighbors 
  real(dp), public :: cutoff
  real(dp), public :: TempMin, TempMax, TempStep
  integer, public :: selTypeModes
  integer, public :: order
  real(dp), public :: atTemperature
  logical, public :: tCompModes
  logical, public :: tPlotModes
  logical, public :: tAnimateModes
  logical, public :: tXmakeMol
  logical, public :: tTransport
  logical, public :: tPhonDispersion
  integer, allocatable, public :: modesToPlot(:)
  integer, public :: nModesToPlot
  integer, public :: nCycles
  integer, public, parameter :: nSteps = 10
  
  !! Locally created variables

  !! Version of the current parser
  integer, parameter :: parserVersion = 4

  !! Version of the oldest parser, for which compatibility is maintained
  integer, parameter :: minVersion = 4 

  !! Container type for parser related flags.
  type TParserFlags
    logical :: tStop                        ! stop after parsing?
    logical :: tIgnoreUnprocessed           ! Continue despite unprocessed nodes
    logical :: tWriteXML, tWriteHSD         ! XML or HSD output?
  end type TParserFlags

  integer, parameter :: ALLMODES = 1
  integer, parameter :: INPLANE = 2 
  integer, parameter :: OUTOFPLANE = 3 

contains

  !!* Initialise program variables
  subroutine initProgramVariables()

    type(fnode), pointer :: input, root, node, tmp
    type(fnode), pointer :: child, value
    type(string) :: buffer, buffer2, modif
    integer :: inputVersion
    integer :: ii, iSp1, iAt
    logical :: tHSD, reqMass
    real(dp), allocatable :: speciesMass(:)
    integer :: nDerivs
    type(TParserflags) :: parserFlags
  
    integer :: cubicType, quarticType
   
    
    !! Write header
    write (*, "(A)") repeat("=", 80)
    write (*, "(A)") "      PHONONS   " // version
    write (*, "(A,/)") repeat("=", 80)

    !! Read in input file as HSD or XML.
    call readHSDOrXML(hsdInput, xmlInput, rootTag, input, tHSD)
    if (tHSD) then
      write (*, "(A)") "Interpreting input file '" // hsdInput // "'"
    else
      write (*, "(A)") "Interpreting input file '" // xmlInput //  "'"
    end if
    write (*, "(A)") repeat("-", 80)
    call getChild(input, rootTag, root)

    !! Check if input version is the one, which we can handle
    !! Handle parser options
    call getChildValue(root, "ParserOptions", tmp, "", child=child, &
        &list=.true., allowEmptyValue=.true.)
    call readParserOptions(child, root, parserFlags)

   ! call getChildValue(root, "InputVersion", inputVersion, parserVersion)
   ! if (inputVersion /= parserVersion) then
   !   call error("Version of input (" // i2c(inputVersion) // ") and parser (" &
   !       &// i2c(parserVersion) // ") do not match")
   ! end if

    call getChild(root, "Geometry", tmp)
    call readGeometry(tmp, geo)

    ! Read Transport block
    ! This defines system partitioning 
    call getChild(root, "Transport", child, requested=.false.)
    if (associated(child)) then
      tTransport = .true.
      call readTransportGeometry(child, geo, transpar)
      call fill_TStructure(geo, str)
    else
      tTransport = .false.
    end if

    call getChildValue(root, "Atoms", buffer2, "1:-1", child=child)
    call convAtomRangeToInt(char(buffer2), geo%specieNames, geo%species, &
        &child, iMovedAtoms)
    nMovedAtom = size(iMovedAtoms)

    call getChild(root, "ComputeModes",child=node,requested=.false.)
    if (associated(node)) then
      tCompModes = .true.
      call getChild(root, "DisplayModes",child=node,requested=.false.)
      if (associated(node)) then
        tPlotModes = .true.
        call getChildValue(node, "PlotModes", buffer2, "1:-1", child=child, &
            &multiple=.true.)
        call convRangeToInt(char(buffer2), child, modesToPlot, 3 * nMovedAtom)
        nModesToPlot = size(modesToPlot)
        call getChildValue(node, "Animate", tAnimateModes, .true.)
        call getChildValue(node, "XMakeMol", tXmakeMol, .true.)
      else
        nModesToPlot = 0
        tPlotModes = .false.
        tAnimateModes = .false.
        tXmakeMol = .false.
      end if
    else
      tCompModes = .false.
    end if

    if (tAnimateModes.and.tXmakeMol) then
      nCycles = 1
    else
      nCycles = 3
    end if
!!!!!!!!!!!!
    ! Reading K-points for Phonon Dispersion calculation
    call getChild(root, "PhononDispersion", child=node, requested=.false.)
    if  (associated(node))  then
      tPhonDispersion = .true.
      call getChildValue(node, "nAtomUnitCell", nAtomUnitCell, 0)
      call getChild(node, "KPoints", child=value, requested=.true.)
      call readKPointsFile(value)
    else
      tPhonDispersion = .false.
    end if
!!!!!!!!!!
    
    ALLOCATE_(speciesMass,(geo%nSpecie))

    ! Read the atomic masses from SlaterKosterFiles or Masses
    call getChild(root, "SlaterKosterFiles", child=value,requested=.false.)
    if ( associated(value) ) then
      call readSKfiles(value, geo, speciesMass)
    else
      call getChild(root,"Masses",child=value)
      call readMasses(value, geo, speciesMass)
    endif
    
    ALLOCATE_(atomicMasses,(nMovedAtom))
    do iAt = 1, nMovedAtom
      atomicMasses(iAt) = speciesMass(geo%species(iMovedAtoms(iAt)))
    end do
    DEALLOCATE_(speciesMass)

    ! --------------------------------------------------------------------------------------
    ! Reading Hessian block parameters 
    ! --------------------------------------------------------------------------------------
    call getChild(root, "Hessian", child=child)
    ! cutoff used to cut out interactions
    call getChildValue(child, "Cutoff", cutoff, 9.45_dp, modifier=modif, child=value) 
    call convertByMul(char(modif), lengthUnits, value, cutoff)

    !selecting the type of modes you want to analysis
    call getchildValue(child, "ModeType", buffer, "all")
    select case(trim(char(buffer)))
    case("all")
      selTypeModes = ALLMODES 
    case("in-plane")
      selTypeModes = INPLANE
    case("out-of-plane")
      selTypeModes = OUTOFPLANE
    case default
      call detailedError(root,"You should specify type of modes: all, in-plane or out-of-plane")
    end select

    ! Reading the actual Hessian matrix
    call getChildValue(child, "Matrix", buffer, child=child)
    !call getNodeName(value, buffer)
    select case(trim(char(buffer)))
    case ("dftb")
      call readDftbHessian(value)
    case ("dynmatrix")
      call readDynMatrix(value)
    case ("cp2k")
      call readCp2kHessian(value)
    case default
      call detailedError(root,"Unkown Hessian type "//char(buffer))  
    end select

    ! --------------------------------------------------------------------------------------

    ! --------------------------------------------------------------------------------------
    ! Reading cubic forces
    ! --------------------------------------------------------------------------------------
    order = 2
    call getChild(root, "Cubic", child=child, requested=.false.)
    if (associated(child)) then
      order = 3
      call getChildValue(child, "Matrix", buffer, child=child)
      select case(trim(char(buffer)))
      case ("gaussian")
        cubicType = 1 
      case default
        call detailedError(root,"Unkown Cubic forces type "//char(buffer))  
      end select
    end if

    call buildNeighborList()
    !call cutDynMatrix()
    ! Hacking to remove in-plane or out-of-plane modes assuming 2D structure on the xy plane
    call selectModes()

    call getChildValue(root, "Analysis", tmp, "", child=child, list=.true., &
        &allowEmptyValue=.true., dummyValue=.true.)

    if (associated(tmp)) then
      call readAnalysis(child, geo, pdos, tundos, transpar, atTemperature)
    endif   


    !! Issue warning about unprocessed nodes
    call warnUnprocessedNodes(root, parserFlags%tIgnoreUnprocessed )

    !! Dump processed tree in HSD and XML format
    if (ioproc .and. parserFlags%tWriteHSD) then
      call dumpHSD(input, hsdParsedInput)
      write(*, '(/,/,A)') "Processed input in HSD format written to '" &
          &// hsdParsedInput // "'"
    end if
    if (ioproc .and. parserFlags%tWriteXML) then
      call dumpHSDAsXML(input, xmlParsedInput)
      write(*, '(A,/)') "Processed input in XML format written to '" &
          &// xmlParsedInput // "'"
    end if
    
    !call destroyNode(input)

    !! Stop, if only parsing is required
    if (parserFlags%tStop) then
      call error("Keyword 'StopAfterParsing' is set to Yes. Stopping.")
    end if



    
  end subroutine initProgramVariables

  !!* Destroy the program variables created in initProgramVariables
  subroutine destructProgramVariables()
    
    call destruct(geo)
    DEALLOCATE_(atomicMasses)
    DEALLOCATE_(dynMatrix)
    DEALLOCATE_(iMovedAtoms)
    DEALLOCATE_(modesToPlot)
    write (*, "(/,A)") repeat("=", 80)
    
  end subroutine destructProgramVariables


  !!* Read in parser options (options not passed to the main code)
  !!* @param node Node to get the information from
  !!* @param root Root of the entire tree (in the case it must be converted)
  !!* @param flags Contains parser flags on exit.
  subroutine readParserOptions(node, root, flags)
    type(fnode), pointer :: node
    type(fnode), pointer :: root
    type(TParserFlags), intent(out) :: flags

    integer :: inputVersion
    type(fnode), pointer :: child

    !! Check if input needs compatibility conversion.
    call getChildValue(node, "ParserVersion", inputVersion, parserVersion, &
        &child=child)
    if (inputVersion < 1 .or. inputVersion > parserVersion) then
      call detailedError(child, "Invalid parser version (" // i2c(inputVersion)&
          &// ")")
    elseif (inputVersion < minVersion) then
      call detailedError(child, &
          &"Sorry, no compatibility mode for parser version " &
          &// i2c(inputVersion) // " (too old)")
    !elseif (inputVersion /= parserVersion) then
    !  write(*, "(A,I2,A,I2,A)") "***  Converting input from version ", &
    !      &inputVersion, " to version ", parserVersion, " ..."
    !  call convertOldHSD(root, inputVersion, parserVersion)
    !  write(*, "(A,/)") "***  Done."
    end if

    call getChildValue(node, "WriteHSDInput", flags%tWriteHSD, .true.)
    call getChildValue(node, "WriteXMLInput", flags%tWriteXML, .false.)
    if (.not. (flags%tWriteHSD .or. flags%tWriteXML)) then
      call detailedWarning(node, &
          &"WriteHSDInput and WriteXMLInput both turned off. You won't &
          &eventually be " &
          &// newline // &
          &" able to obtain the same results with a later version of the code!")
    end if
    call getChildValue(node, "StopAfterParsing", flags%tStop, .false.)

    call getChildValue(node, "IgnoreUnprocessedNodes", &
        &flags%tIgnoreUnprocessed, .false.)

  end subroutine readParserOptions


  !!* Read in the geometry stored as xml in internal or gen format.
  !!* @param geonode Node containing the geometry
  !!* @param geo     Contains the geometry information on exit
  subroutine readGeometry(geonode, geo)
    type(fnode), pointer :: geonode
    type(TGeometry), intent(out) :: geo

    type(fnode), pointer :: child, value
    type(string) :: buffer

    call getChildValue(geonode, "", value, child=child)
    call getNodeName(value, buffer)
    select case (char(buffer))
    case ("genformat")
      call readTGeometryGen(value, geo)
      !call removeChildNodes(geonode)
      !call writeTGeometryHSD(geonode, geo)
    case default
      call setUnprocessed(value)
      call readTGeometryHSD(child, geo)
    end select
    call unstring(buffer)
    
  end subroutine readGeometry

  !!* Read geometry information for transport calculation
  subroutine readTransportGeometry(root, geom, tp)
    type(fnode), pointer :: root
    type(TGeometry), intent(inout) :: geom
    type(TTransPar), intent(inout) :: tp

    type(fnode), pointer :: pGeom, pDevice, pNode, pTask, pTaskType
    type(string) :: buffer, modif
    type(fnode), pointer :: pTmp, field
    type(fnodelist), pointer :: pNodeList
    !type(fnodeList), pointer :: pNodeList
    integer :: ii, contact
    real(dp) :: acc, contactRange(2), sep
    
    tp%defined = .true.
    tp%tPeriodic1D = .not. geom%tPeriodic
    call getChild(root, "Device", pDevice)
    call getChildValue(pDevice, "AtomRange", tp%idxdevice)
    call getChild(pDevice, "FirstLayerAtoms", pTmp, requested=.false.)
    call readFirstLayerAtoms(pTmp, tp%PL, tp%nPLs, tp%idxdevice)
    if (.not.associated(pTmp)) then
      call setChildValue(pDevice, "FirstLayerAtoms", tp%PL)
    end if
    call getChildren(root, "Contact", pNodeList)
    tp%ncont = getLength(pNodeList)
    if (tp%ncont < 2) then
      call detailedError(pGeom, "At least two contacts must be defined")
    end if
    ALLOCATE_(tp%contacts, (tp%ncont))
    !! Parse contact geometry
    call readContacts(pNodeList, tp%contacts, geom, (buffer .eq. "uploadcontacts"))

    !! Note: here "Task" is always like uplaodContacts
    !! => no need for all the mess

    !call destroyNodeList(pNodeList)

  end subroutine readTransportGeometry

  subroutine readFirstLayerAtoms(pnode, pls, npl, idxdevice, check)
    logical, optional :: check
    type(fnode), pointer, intent(in) :: pnode
    integer :: idxdevice(2)
    integer, allocatable :: pls(:)
    integer :: npl
  
    type(listInt) :: li
    logical :: checkidx

    checkidx = .true.
    if (present(check)) checkidx = check
    
    if (associated(pnode)) then
        call init(li)
        call getChildValue(pnode, "", li)
        npl = len(li)
        ALLOCATE_(pls, (npl))
        call asArray(li, pls)
        call destroy(li)
        if (checkidx) then
          if (any(pls < idxdevice(1) .or. &
                  pls > idxdevice(2))) then
             call detailedError(pnode, "First layer atoms must be between " &
               &// i2c(idxdevice(1)) // " &
               & and " // i2c(idxdevice(2)) // ".")
          end if
        end if
      else
         npl = 1
         ALLOCATE_(pls, (npl))
         pls = (/ 1 /)
      end if

  end subroutine readFirstLayerAtoms

   !!* Read bias information, used in Analysis and Green's function eigensolver 
  subroutine readContacts(pNodeList, contacts, geom, upload)
    type(ContactInfo), allocatable, dimension(:), intent(inout) :: contacts
    type(fnodeList), pointer :: pNodeList
    type(TGeometry), intent(in) :: geom
    logical, intent(in) :: upload

    real(dp) :: acc
    integer :: ncont, ii
    type(fnode), pointer :: field, pNode, pTmp, pWide
    type(string) :: buffer, modif

    ncont = size(contacts)
    do ii = 1,ncont
      contacts(ii)%wideBand = .false.
      contacts(ii)%wideBandDos = 0.0

      call getItem1(pNodeList, ii, pNode)
      call getChildValue(pNode, "Id", buffer, child=pTmp)
      buffer = tolower(trim(unquote(char(buffer))))
      if (len(buffer) > mc) then
        call detailedError(pTmp, "Contact id may not be longer than " &
            &// i2c(mc) // " characters.")
      end if
      contacts(ii)%name = char(buffer)
      if (any(contacts(1:ii-1)%name == contacts(ii)%name)) then
        call detailedError(pTmp, "Contact id '" // trim(contacts(ii)%name) &
            &//  "' already in use")
      end if
      call getChildValue(pNode,"ShiftAccuracy",acc, 1e-5_dp, modifier=modif&
          &, child=field)
      call convertByMul(char(modif), lengthUnits, field, acc)
      call getChildValue(pNode, "AtomRange", contacts(ii)%idxrange, child=pTmp)
      !if (acc > 1.0_dp) then
      !  contactVecs(:,ii) = (/ 0.d0, 0.d0, acc /)
      !  ginfo%transport%cdir(ii) = 3
      !else
      call getContactVectorII(contacts(ii)%idxrange, geom, ii, pTmp, acc, &
                              & contacts(ii)%lattice, contacts(ii)%dir)
      !endif
      contacts(ii)%length = sqrt(sum(contacts(ii)%lattice**2))

      ! Contact temperatures. A negative default is used so it is quite clear
      ! when the user sets a different value. In such a case 
      ! this overrides values defined in Filling block
      call getChildValue(pNode, "Temperature", contacts(ii)%kbT,&
                         &0.0_dp, modifier=modif, child=field)
      call convertByMul(char(modif), energyUnits, field, contacts(ii)%kbT)

      if (upload) then
        !call getChildValue(pNode, 'potential', contacts(ii)%potential,&
        !                    &0.0_dp, modifier=modif, child=field)
        !call convertByMul(char(modif), energyUnits, field, contacts(ii)%potential)
        contacts(ii)%potential = 0.d0

        call getChildValue(pNode, "wideBand", contacts(ii)%wideBand, .false.)
        if (contacts(ii)%wideBand) then
          call getChildValue(pNode, 'LevelSpacing', contacts(ii)%wideBandDos, &
                             &0.735_dp, modifier=modif, child=field) 
          call convertByMul(char(modif), energyUnits, field,&
                              &contacts(ii)%wideBandDos) 
          !WideBandApproximation is defined as energy spacing between levels
          !In the code the inverse value (Density of states) is used
          !Convert the negf input value. Default is 20.e eV
          contacts(ii)%wideBandDos = 1.d0 / contacts(ii)%wideBandDos
        end if
        !call getChildValue(pNode, "FermiLevel", contacts(ii)%eFermi,modifier=modif)
        !call convertByMul(char(modif), energyUnits, pNode, contacts(ii)%eFermi)
        contacts(ii)%eFermi=0.d0
      end if

    enddo

  end subroutine readContacts

      ! Sanity checking of atom ranges and returning contact vector and direction.
  subroutine getContactVectorII(atomrange, geom, id, pContact, acc, &
      &contactVec, contactDir)
    integer, intent(in) :: atomrange(2)
    type(TGeometry), intent(in) :: geom
    integer, intent(in) :: id
    type(fnode), pointer :: pContact
    real(dp), intent(in) :: acc
    real(dp), intent(out) :: contactVec(3)
    integer, intent(out) :: contactDir

    integer :: iStart, iStart2, iEnd
    logical :: mask(3)

    !! Sanity check for the atom ranges
    iStart = atomrange(1)
    iEnd = atomrange(2)
    if (iStart < 1 .or. iEnd < 1 .or. iStart > geom%nAtom &
        &.or. iEnd > geom%nAtom .or. iEnd < iStart) then
      call detailedError(pContact, "Invalid atom range '" // i2c(iStart) &
          &// " " // i2c(iEnd) // "', values should be between " // i2c(1) &
          &// " and " // i2c(geom%nAtom) // ".")
    end if
    if (mod(iEnd - iStart + 1, 2) /= 0) then
      call detailedError(pContact, "Nr. of atoms in the contact must be even")
    end if

    ! Determining contact vector
    iStart2 = iStart + (iEnd - iStart + 1) / 2
    contactVec = geom%coords(:,iStart) - geom%coords(:,iStart2)
    if (any(sqrt(sum(&
        &(geom%coords(:,iStart:iStart2-1) - geom%coords(:,iStart2:iEnd) &
        &- spread(contactVec, dim=2, ncopies=iStart2-iStart))**2, dim=1)) &
        &> acc)) then
      write(*,*) 'coords:', geom%coords(:,iStart)
      write(*,*) 'coords:', geom%coords(:,iStart2)
      write(*,*) 'Contact Vector:', contactVec(1:3)
      write(*,*) iStart,iStart2,iEnd
      write(*,*) 'X:'
      write(*,*) ((geom%coords(1,iStart:iStart2-1)&
          & - geom%coords(1,iStart2:iEnd)&
          & - spread(contactVec(1), dim=1, ncopies=iStart2-iStart)))
      write(*,*) 'Y:'
      write(*,*) ((geom%coords(2,iStart:iStart2-1)&
          & - geom%coords(2,iStart2:iEnd) &
          & - spread(contactVec(2), dim=1, ncopies=iStart2-iStart)))
      write(*,*) 'Z:'
      write(*,*) ((geom%coords(3,iStart:iStart2-1)&
          & - geom%coords(3,iStart2:iEnd) &
          &- spread(contactVec(3), dim=1, ncopies=iStart2-iStart)))
      call error("Contact " // i2c(id) &
          &// " does not consist of two rigidly shifted layers")
    end if

    ! Determine to which axes it is parallel.
!    mask = (abs(abs(contactVec)  - sqrt(sum(contactVec**2))) < acc)
!    if (count(mask) /= 1) then
!      call error("Contact vector " // i2c(id) // " not parallel to any&
!          & of the coordinate axis.")
!    end if
    ! Workaround for bug in Intel compiler (can not use index function)
    contactDir = 0
!    do while (.not. mask(contactDir))
!      contactDir = contactDir + 1
!    end do

  end subroutine getContactVectorII



  ! this is a utility subroutine to fill in negf structure
  subroutine fill_TStructure(geo, str)
    type(TGeometry), intent(inout), target :: geo
    type(TGDFTBStructure), intent(out) :: str 
   
    integer :: ii

    str%nAtom = geo%nAtom
    str%nSpecies = geo%nSpecie
    str%specie0 => geo%species
    str%x0 => geo%coords
    str%isperiodic = geo%tPeriodic
    if (str%isperiodic) then
      str%latVecs = geo%latVecs
    else
      str%latVecs = 0.0_dp 
    end if

    str%nel = 0.0
    str%tempElec = 0.0

    INITALLOCATE_PARR(str%iAtomStart, (str%nAtom+1))
    do ii = 1, str%nAtom+1
      str%iAtomStart(ii) = 1+3*(ii-1)
    end do

  end subroutine fill_TStructure


  subroutine readSKfiles(child, geo, speciesMass)
    type(fnode), pointer :: child
    type(TGeometry), intent(in) :: geo
    real(dp), dimension(:) :: speciesMass

    type(TOldSKData) :: skData
    type(listCharLc), allocatable :: skFiles(:)
    type(fnode), pointer :: value, child2
    type(string) :: buffer, buffer2
    character(lc) :: prefix, suffix, separator, elem1, elem2, strTmp, filename
    type(listString) :: lStr
    integer :: ii, iSp1
    logical :: tLower, tExist

    !! Slater-Koster files
    ALLOCATE_(skFiles, (geo%nSpecie))
    do iSp1 = 1, geo%nSpecie
        call init(skFiles(iSp1))
    end do

    call getChildValue(child, "", value)
    call getNodeName(value, buffer)

    select case(char(buffer))
    case ("type2filenames")
      call getChildValue(value, "Prefix", buffer2, "")
      prefix = unquote(char(buffer2))
      call getChildValue(value, "Suffix", buffer2, "")
      suffix = unquote(char(buffer2))
      call getChildValue(value, "Separator", buffer2, "")
      separator = unquote(char(buffer2))
      call getChildValue(value, "LowerCaseTypeName", tLower, .false.)
      do iSp1 = 1, geo%nSpecie
        if (tLower) then
          elem1 = tolower(geo%specieNames(iSp1))
        else
          elem1 = geo%specieNames(iSp1)
        end if
        strTmp = trim(prefix) // trim(elem1) // trim(separator) &
            &// trim(elem1) // trim(suffix)
        call append(skFiles(iSp1), strTmp)
        inquire(file=strTmp, exist=tExist)
        if (.not. tExist) then
          call detailedError(value, "SK file with generated name '" &
              &// trim(strTmp) // "' does not exist.")
        end if
      end do
    case default
      call setUnprocessed(value)
      do iSp1 = 1, geo%nSpecie
        strTmp = trim(geo%specieNames(iSp1)) // "-" &
            &// trim(geo%specieNames(iSp1))
        call init(lStr)
        call getChildValue(child, trim(strTmp), lStr, child=child2)
        ! We can't handle selected shells here (also not needed I guess)
        if (len(lStr) /= 1) then
          call detailedError(child2, "Incorrect number of Slater-Koster &
              &files")
        end if
        do ii = 1, len(lStr)
          call get(lStr, strTmp, ii)
          inquire(file=strTmp, exist=tExist)
          if (.not. tExist) then
            call detailedError(child2, "SK file '" // trim(strTmp) &
                &// "' does not exist'")
          end if
          call append(skFiles(iSp1), strTmp)
        end do
        call destroy(lStr)
      end do
    end select

    do iSp1 = 1, geo%nSpecie
      call get(skFiles(iSp1), fileName, 1)
      call readFromFile(skData, fileName, .true.)
      DEALLOCATE_PARR(skData%skHam)
      DEALLOCATE_PARR(skData%skOver)
      speciesMass(iSp1) = skData%mass      
    end do
    
    do iSp1 = 1, geo%nSpecie
      call destroy(skFiles(iSp1))
    end do
    DEALLOCATE_(skFiles)

  end subroutine readSKfiles

  subroutine readMasses(value, geo, speciesMass)
    type(fnode), pointer :: value
    type(TGeometry), intent(in) :: geo
    real(dp), dimension(:) :: speciesMass

    type(fnode), pointer :: child, child2
    type(string) :: modif 
    integer :: iSp
    character(lc) :: strTmp
    real(dp) :: mass, defmass

    do iSp = 1, geo%nSpecie
      defmass = getAtomicMass(trim(geo%specieNames(iSp)))
      call getChildValue(value, geo%specieNames(iSp), mass, defmass,& 
               &modifier=modif, child= child2)
      speciesMass(iSp) = mass * amu__au
    end do

  end subroutine readMasses
  
  subroutine  readKPointsFile(child)
    type(fnode),  pointer ::  child
    type(string) :: text

    call getFirstTextChild(child, text)
    call readKPointsFile_help(child, char(text))

  end subroutine  readKPointsFile

  subroutine readKPointsFile_help(child,text)
    type(fnode),  pointer ::  child
    character(len=*), intent(in) :: text
    integer :: iStart, iErr=0, ii, iOldStart
    real(dp), dimension(:), allocatable :: tmparray
    real(dp), dimension(:,:), allocatable :: KPoint2

    iStart = 1
    call getNextToken(text, nKPoints, iStart, iErr)

      allocate(tmparray(4))
      allocate(KPoint2(nKPoints,4))
      allocate(KPoint(nKPoints,3))
      allocate(KWeight(nKPoints))

    iErr = TOKEN_ERROR
    iOldStart = iStart
    iStart  = iOldStart

    do ii = 1, nKPoints
      call getNextToken(text, tmparray, iStart, iErr)
      KPoint2(ii,:) = tmparray(:)
    end do

    do ii = 1, nKPoints
        KPoint(ii,1:3)  = KPoint2(ii,1:3)
        KWeight(ii)  = KPoint2(ii,4)
    end do

    do ii = 1, nKPoints
      print*, KPoint(ii,1:3)
    enddo

  end subroutine readKPointsFile_help
      
  subroutine readDftbHessian(child)
    type(fnode), pointer :: child

    type(listRealR1) :: realBuffer
    integer :: iCount, jCount, ii, kk, jj, ll 
    integer :: nDerivs

    real, dimension(:,:), allocatable :: h 
    integer ::  n, j1, j2

    nDerivs = 3 * nMovedAtom
    ALLOCATE_(dynMatrix,(nDerivs,nDerivs))

    !The derivatives matrix must be stored as the following order: 
    
    ! For the x y z directions of atoms 1..n
    !   d^2 E        d^2 E       d^2 E       d^2 E        d^2 E         
    ! ---------- + --------- + --------- + ---------- + ---------- +...
    ! dx_1 dx_1    dy_1 dx_1   dz_1 dx_1   dx_2 dx_1    dy_2 dx_1   

!    call init(realBuffer)
!    call getChildValue(child, "", nDerivs, realBuffer)
!    if (len(realBuffer)/=nDerivs) then
!      call detailedError(child,"wrong number of derivatives supplied:" &
!          & // i2c(len(realBuffer)) // " supplied, " &
!          & // i2c(nDerivs) // " required.")
!    end if
!    call asArray(realBuffer, dynMatrix)
!    call destroy(realBuffer)

    open(unit=65, file='hessian.out', action='read')
    do ii = 1,  nDerivs
        read(65,'(4f16.10)') dynMatrix(ii,1:nDerivs)
    end do

    ! mass weight the Hessian matrix to get the dynamical matrix
    iCount = 0
    do ii = 1, nMovedAtom
      do kk = 1, 3
        iCount = iCount + 1
        jCount = 0
        do jj = 1, nMovedAtom
          do ll = 1, 3
            jCount = jCount + 1
            dynMatrix(jCount,iCount) = dynMatrix(jCount,iCount) &
                & / (sqrt(atomicMasses(ii)) * sqrt(atomicMasses(jj)))
          end do
        end do
      end do
    end do


  close(65)

  end subroutine readDftbHessian

  subroutine readDynMatrix(child)
    type(fnode), pointer :: child

    type(listRealR1) :: realBuffer
    integer :: iCount, jCount, ii, kk, jj, ll 
    integer :: nDerivs

    nDerivs = 3 * nMovedAtom
    ALLOCATE_(dynMatrix,(nDerivs,nDerivs))

    !The derivatives matrix must be stored as the following order: 
    
    ! For the x y z directions of atoms 1..n
    !   d^2 E        d^2 E       d^2 E       d^2 E        d^2 E         
    ! ---------- + --------- + --------- + ---------- + ---------- +...
    ! dx_1 dx_1    dy_1 dx_1   dz_1 dx_1   dx_2 dx_1    dy_2 dx_1   

    call init(realBuffer)
    call getChildValue(child, "", nDerivs, realBuffer)
    if (len(realBuffer)/=nDerivs) then
      call detailedError(child,"wrong number of derivatives supplied:" &
          & // i2c(len(realBuffer)) // " supplied, " &
          & // i2c(nDerivs) // " required.")
    end if
    call asArray(realBuffer, dynMatrix)
    call destroy(realBuffer)
  
  end subroutine readDynMatrix

  subroutine readCp2kHessian(child)
    type(fnode), pointer :: child

    type(listRealR1) :: realBuffer
    integer :: iCount, jCount, ii, kk, jj, ll 
    integer :: nDerivs, nBlocks

    real, dimension(:,:), allocatable :: HessCp2k
    integer ::  n, j1, j2,  p,  q

    nDerivs = 3 * nMovedAtom
    ALLOCATE_(dynMatrix,(nDerivs,nDerivs))

    !The derivatives matrix must be stored as the following order: 
    
    ! For the x y z directions of atoms 1..n
    !   d^2 E        d^2 E       d^2 E       d^2 E        d^2 E         
    ! ---------- + --------- + --------- + ---------- + ---------- +...
    ! dx_1 dx_1    dy_1 dx_1   dz_1 dx_1   dx_2 dx_1    dy_2 dx_1   

!    call init(realBuffer)
!    call getChildValue(child, "", nDerivs, realBuffer)
!    if (len(realBuffer)/=nDerivs) then
!      call detailedError(child,"wrong number of derivatives supplied:" &
!          & // i2c(len(realBuffer)) // " supplied, " &
!          & // i2c(nDerivs) // " required.")
!    end if
!    call asArray(realBuffer, dynMatrix)
!    call destroy(realBuffer)

!    n = nDerivs*nDerivs/4

    open(unit=65, file='hessian.cp2k', action='read')
    nBlocks = nDerivs/5.0
    
    allocate(HessCp2k(nDerivs*nBlocks,5))

    do  ii  = 1,  nDerivs*nBlocks
        read(65,*) HessCp2k(ii,1:5)
    end do

    do ii = 1,  nBlocks
        do  jj  = 1, nDerivs
            p = 1+5*(ii-1)
            q = 5*ii
          dynMatrix(jj,p:q) =  HessCp2k(jj + nDerivs*(ii-1),1:5)
        end do
    end do

    ! mass weight the Hessian matrix to get the dynamical matrix
    iCount = 0
    do ii = 1, nMovedAtom
      do kk = 1, 3
        iCount = iCount + 1
        jCount = 0
        do jj = 1, nMovedAtom
          do ll = 1, 3
            jCount = jCount + 1
            dynMatrix(jCount,iCount) = dynMatrix(jCount,iCount) &
                & / (sqrt(atomicMasses(ii)) * sqrt(atomicMasses(jj)))
          end do
        end do
      end do
    end do


  close(65)

  end subroutine readCp2kHessian

  subroutine selectModes()
   
    integer :: iCount, jCount, ii, jj, kk, ll
    
    select case ( selTypeModes )
    case(INPLANE)
      iCount = 0
      do ii = 1, nMovedAtom
        do kk = 1, 3
          iCount = iCount + 1
          jCount = 0
          do jj = 1, nMovedAtom
            do ll = 1, 3
              jCount = jCount + 1
              if (mod(iCount,3).eq.0 .or. mod(jCount,3).eq.0) then
                  dynMatrix(jCount,iCount) = 0.0
              end if
            end do
          end do
        end do
      end do
    case(OUTOFPLANE)
      iCount = 0
      do ii = 1, nMovedAtom
        do kk = 1, 3
          iCount = iCount + 1
          jCount = 0
          do jj = 1, nMovedAtom
            do ll = 1, 3
              jCount = jCount + 1
              if (mod(iCount,3).ne.0 .and. mod(jCount,3).ne.0) then
                  dynMatrix(jCount,iCount) = 0.0
              end if
            end do
          end do
        end do
      end do
    end select 

  end subroutine selectModes    

  !! Reads the Analysis block.
  subroutine readAnalysis(node, geo, pdos, tundos, transpar, atTemperature)
    type(fnode), pointer :: node, pnode
    type(TGeometry), intent(in) :: geo
    type(TPdos), intent(inout) :: pdos
    type(TGDFTBTunDos), intent(inout) :: tundos
    type(TTransPar), intent(inout) :: transpar
    real(dp) :: atTemperature, TempRange(2)

    type(fnode), pointer :: val, child, field
    type(string) :: modif 
    type(fnodeList), pointer :: children
   

    !call getChildValue(node, "ProjectStates", val, "", child=child, &
    !    &allowEmptyValue=.true., list=.true.)
    
    !if (associated(child)) then
      !call detailedError(node,"ProjectStates available only in TunnelingAndDOS")
      !call getChildren(child, "Region", children)
      !call readPDOSRegions(children, geo, pdos%iAtInRegion, pdos%regionLabels)
      !call destroyNodeList(children)
    !end if
   
    call getChild(node, "TunnelingAndDOS", child, requested=.false.)
    if (associated(child)) then
      if (.not.tTransport) then
        call detailedError(node, "Tunneling requires Transport block")
      end if
      call readTunAndDos(child, geo, tundos, transpar, maxval(transpar%contacts(:)%kbT) )
    endif

    call getChild(node, "Conductance", child, requested=.false.)
    if (associated(child)) then
      if (.not.tTransport) then
        call detailedError(node, "Conductance requires Transport block")
      end if
      call getChildValue(child, "TempRange", TempRange, modifier=modif,&
      & child=field)
      call convertByMul(char(modif), energyUnits, field, TempRange)
    
      call getChildValue(child, "TempStep", TempStep, modifier=modif,&
      & child=field)
      call convertByMul(char(modif), energyUnits, field, TempStep)
    
       TempMin = TempRange(1)
       TempMax = TempRange(2)
    endif
   
  end subroutine readAnalysis
  
  subroutine readPDOSRegions(children, geo, iAtInregion, regionLabels)
    type(fnodeList), pointer :: children
    type(TGeometry), intent(in) :: geo
    type(WrappedInt1), allocatable, intent(out) :: iAtInRegion(:)
    character(lc), allocatable, intent(out) :: regionLabels(:)

    integer :: nReg, iReg
    integer, allocatable :: tmpI1(:)
    type(fnode), pointer :: child, child2
    type(string) :: buffer
    character(lc) :: strTmp

    nReg = getLength(children)
    ALLOCATE_(regionLabels, (nReg))
    ALLOCATE_(iAtInRegion, (nReg))
    do iReg = 1, nReg
      call getItem1(children, iReg, child)
      call getChildValue(child, "Atoms", buffer, child=child2, &
          & multiple=.true.)
      call convAtomRangeToInt(char(buffer), geo%specieNames, &
          & geo%species, child2, tmpI1)
      iAtInRegion(iReg)%data = tmpI1      
      write(strTmp, "('region',I0)") iReg
      call getChildValue(child, "Label", buffer, trim(strTmp))
      regionLabels(iReg) = unquote(char(buffer))
    end do
    call unstring(buffer)
    
  end subroutine readPDOSRegions

  !!* Read Tunneling and Dos options from analysis block
  !!* tundos is the container to be filled
  !!* ncont is needed for contact option allocation
  subroutine readTunAndDos(root, geo, tundos, transpar, temperature)
    type(fnode), pointer :: root
    type(TGeometry), intent(in) :: geo
    type(TGDFTBTunDos), intent(inout) :: tundos
    type(TTransPar), intent(inout) :: transpar
    real(dp), intent(in) :: temperature

    type(fnode), pointer :: pTmp, field
    type(fnode), pointer :: pGeom, pDevice, pNode
    type(fnodeList), pointer :: pNodeList
    integer :: ii, jj, ind, ncont, nKT
    real(dp) :: eRange(2), eRangeDefault(2) 
    type(string) :: buffer, modif
    type(WrappedInt1), allocatable :: iAtInRegion(:)
    logical, allocatable :: tDirectionResInRegion(:)
    character(lc), allocatable :: regionLabelPrefixes(:)

    tundos%defined = .true.
    ncont = transpar%ncont
    call getChildValue(root, "Verbosity", tundos%verbose, 51)
    call getChildValue(root, "WriteLDOS", tundos%writeLDOS, .true.)
    call getChildValue(root, "WriteTunn", tundos%writeTunn, .true.)
    
    !call getChildValue(root, "Temperature", temperature, &
    !    modifier=modif, child=field,requested=.false. )
    !call convertByMul(char(modif), energyUnits, field, temperature)    
    ! Parsing of energy range
    ! If the calculation is in equilibrium (all potentials to 0.0)
    ! then an energy range and step must be specified (it is assumed
    ! that the user use this filed to calculate a DOS or T(E) )
    ! If the calculation is out of equilibrium, a default similar to
    ! GreensFunction RealAxisStep is set to ensure that the current 
    ! can be calculated without manually specify the energy parameters.

    ! Default meaningful: eRange= (0..10*kT]
    ! nKT is set to GreensFunction default, i.e. 10
    ! I avoid an explicit nKT option because I find it confusing here 
    ! (it makes sense only out of equilibrium)
    ! What matters is w*[nB(w;T1)-nB(w;T2)] that is finite lim w->0
    nKT = 10
    eRangeDefault(1) = 0.0001_dp
    eRangeDefault(2) = nKT * temperature 

    call getChildValue(root, "FreqRange", eRange, eRangeDefault, &
         & modifier=modif, child=field)
    call convertByMul(char(modif), energyUnits, field, eRange)

    if (eRange(1).le.0.d0) then
       call detailedError(root, "FreqRange must be > 0")
    end if 

    call getChildValue(root, "FreqStep", tundos%estep, 1.0e-5_dp,  &
       & modifier=modif, child=field)

    call convertByMul(char(modif), energyUnits, field, tundos%estep)

    ! Terminal currents
    call getChild(root, "TerminalCurrents", pTmp, requested=.false.)
      if (associated(pTmp)) then
        call getChildren(pTmp, "EmitterCollector", pNodeList)
        ALLOCATE_(tundos%ni, (getLength(pNodeList)))
        ALLOCATE_(tundos%nf, (getLength(pNodeList)))
        do ii = 1, getLength(pNodeList)
          call getItem1(pNodeList, ii, pNode)
          call getEmitterCollectorByName(pNode, tundos%ni(ii),&
              & tundos%nf(ii), transpar%contacts(:)%name)
        end do
        call destroyNodeList(pNodeList)
      else
        ALLOCATE_(tundos%ni, (ncont-1) )
        ALLOCATE_(tundos%nf, (ncont-1) )
        call setChild(root, "TerminalCurrents", pTmp)
        ind = 1
        do ii = 1, 1
          do jj = ii + 1, ncont
            call setChildValue(pTmp, "EmitterCollector", &
                &(/ transpar%contacts(ii)%name, transpar%contacts(jj)%name /))
            tundos%ni(ind) = ii
            tundos%nf(ind) = jj
            ind = ind + 1
          end do
        end do
      end if
      call getChildValue(root, "Delta", tundos%delta, &
          &1.0e-7_dp, modifier=modif, child=field)
      call convertByMul(char(modif), energyUnits, field, &
          &tundos%delta)
      call getChildValue(root, "BroadeningDelta", tundos%broadeningDelta, &
          &0.0_dp, modifier=modif, child=field)
      call convertByMul(char(modif), energyUnits, field, &
          &tundos%broadeningDelta)

      call getChildren(root, "Region", pNodeList)
      call readPDOSRegions(pNodeList, geo, iAtInRegion, regionLabelPrefixes)
      call destroyNodeList(pNodeList)

      call addAtomResolvedRegion(tundos%dosOrbitals, tundos%dosLabels)

      tundos%emin = eRange(1)
      tundos%emax = eRange(2)
    
      
    contains

      !! Adds one region with all the orbitals of the atoms in it.
      subroutine addAtomResolvedRegion(iOrbRegion, regionLabels)
        
        type(WrappedInt1), allocatable, intent(out) :: iOrbRegion(:)
        character(lc), allocatable, intent(out) :: regionLabels(:)
  
        integer :: nRegion, nAtomInRegion, iReg, ind, ii, jj, iAt
        integer :: nIndices
  
        nRegion = size(iAtInRegion)
        ALLOCATE_(iOrbRegion, (nRegion))
        ALLOCATE_(regionLabels, (nRegion))

        do iReg = 1, nRegion
          nAtomInRegion = size(iAtInRegion(iReg)%data) 
          nIndices = 3*nAtomInRegion
          ALLOCATE_(iOrbRegion(iReg)%data, (nIndices))
          ind = 1
          do ii = 1, nAtomInRegion
            iAt = iAtInRegion(iReg)%data(ii)
            do jj = 0, 2
              iOrbRegion(iReg)%data(ind) = 3*iAt - 2 + jj 
              ind = ind + 1
            end do
          end do
          regionLabels(iReg) = regionLabelPrefixes(iReg)
        end do
       
      end subroutine addAtomResolvedRegion

  end subroutine readTunAndDos
    
  ! Get contacts for terminal currents by name
  subroutine getEmitterCollectorByName(pNode, emitter, collector, contactNames)
    type(fnode), pointer :: pNode
    integer, intent(out) :: emitter, collector
    character(len=*), intent(in) :: contactNames(:)

    type(listString) :: lString
    character(len=mc) :: buffer
    integer :: ind
    logical :: tFound

    call init(lString)
    call getChildValue(pNode, "", lString)
    if (len(lString) /= 2) then
      call detailedError(pNode, "You must provide two contacts")
    end if
    call get(lString, buffer, 1)
    emitter = getContactByName(contactNames, buffer, pNode)
    call get(lString, buffer, 2)
    collector = getContactByName(contactNames, buffer, pNode)
    call destroy(lString)

  end subroutine getEmitterCollectorByName

  ! Getting the contact by name
  function getContactByName(contactNames, contName, pNode) result(contact)
    character(len=*), intent(in) :: contactNames(:)
    character(len=*), intent(in) :: contName
    type(fnode), pointer :: pNode
    integer :: contact

    logical :: tFound

    tFound = .false.
    do contact = 1, size(contactNames)
      tFound = (contactNames(contact) == contName)
      if (tFound) then
        exit
      end if
    end do
    if (.not. tFound) then
      call detailedError(pNode, "Invalid collector contact name '" &
          &// trim(contName) // "'")
    end if

  end function getContactByName

  ! Build a simple neighbor list. Currently does not work for periodic systems.
  ! Have to fix this important point
  subroutine buildNeighborList()
  
    integer ::  iAtom, jAtom, ii, jj, kk, PL1, PL2
    integer, parameter :: nInitNeighbor = 100  !* First guess for nr. of neighbors.
    real :: disAtom, dd(3) 
    integer :: nAllAtom
    real(dp) :: mCutoff
    real(dp), pointer :: coords(:,:), cellVec(:,:), rCellVec(:,:)
    integer, pointer :: iCellVec(:)
    
    call init(neighborList, geo%nAtom, nInitNeighbor)

    mCutoff = 1.0_dp * cutoff

    if (geo%tPeriodic) then
      !! Make some guess for the nr. of all interacting atoms
      nAllAtom = int((real(geo%nAtom, dp)**(1.0_dp/3.0_dp) + 3.0_dp)**3)
      INIT_PARR(cellVec)
      INIT_PARR(rCellVec)
      call getCellTranslations(cellVec, rCellVec, geo%latVecs, &
            &geo%recVecs2p, mCutoff)
      DEALLOCATE_PARR(cellVec)
    else
      nAllAtom = geo%nAtom
      INITALLOCATE_PARR(rCellVec, (3, 1))
      rCellVec(:, 1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
    end if

    INITALLOCATE_PARR(coords, (3, nAllAtom))
    INITALLOCATE_PARR(img2CentCell, (nAllAtom))
    INITALLOCATE_PARR(iCellVec, (nAllAtom))
      
    call updateNeighborList(coords, img2CentCell, iCellVec, neighborList, &
          &nAllAtom, geo%coords, mCutoff, rCellVec)
      
    DEALLOCATE_PARR(coords)
    DEALLOCATE_PARR(iCellVec)
    DEALLOCATE_PARR(rCellVec)
    
    ALLOCATE_(nNeighbor, (geo%nAtom))
    nNeighbor(:) = 0
  
    call getNrOfNeighborsForAll(nNeighbor, neighborList, mCutoff)

    !print*, 'cutoff=',neighborList%cutoff 
    !do iAtom = 1, nMovedAtom
    !  print*,'nNeig=',iAtom, nNeighbor(iAtom)
    !  print*,img2CentCell(neighborList%iNeighbor(1:nNeighbor(iAtom),iAtom))
    !  print*,neighborList%neighDist2(1:nNeighbor(iAtom),iAtom)
    !end do
    
    ! Check PLs
    do iAtom = 1, transpar%idxdevice(2)
      PL1 = getPL(iAtom)
      do jj = 1, nNeighbor(iAtom)
        jAtom = img2CentCell(neighborList%iNeighbor(jj,iAtom))
        if (jAtom > transpar%idxdevice(2)) cycle 
        PL2 = getPL(jAtom)
        if (.not.(PL1.eq.PL2 .or. PL1.eq.PL2+1 .or. PL1.eq.PL2-1)) then
          write(*,*) 'ERROR: PL size inconsistent with cutoff'
          stop
        end if
      end do 
    end do

  end subroutine buildNeighborList

  subroutine cutDynMatrix()  

    integer :: iAtom, jAtom, jj
    real(dp), allocatable :: dynMat2(:,:)

    ALLOCATE_(dynMat2, (3*nMovedAtom, 3*nMovedAtom))
    dynMat2 = 0.0_dp

    do iAtom = 1, geo%nAtom
       do jj = 1, nNeighbor(iAtom)
          jAtom = img2CentCell(neighborList%iNeighbor(jj, iAtom))
          if (neighborList%neighDist2(jj,iAtom) .le. cutoff**2) then
            dynMat2(3*(iAtom-1)+1:3*(iAtom-1)+3, 3*(jAtom-1)+1:3*(jAtom-1)+3) = &
                dynMatrix(3*(iAtom-1)+1:3*(iAtom-1)+3, 3*(jAtom-1)+1:3*(jAtom-1)+3) 
            dynMat2(3*(jAtom-1)+1:3*(jAtom-1)+3, 3*(iAtom-1)+1:3*(iAtom-1)+3) = &
                dynMatrix(3*(iAtom-1)+1:3*(iAtom-1)+3, 3*(jAtom-1)+1:3*(jAtom-1)+3) 
          end if
       end do
    end do

    dynMatrix = dynMat2 

    DEALLOCATE_(dynMat2)
  
  end subroutine cutDynMatrix

  function getPL(iAt) result(PL)
    integer, intent(in) :: iAt
    integer :: PL
 
    integer :: ii

    do ii = 1, transpar%nPLs-1
      if (iAt>=transpar%PL(ii) .and. iAt<transpar%PL(ii+1)) then
        PL = ii
      end if
    end do
    if ( iAt>=transpar%PL(transpar%nPLs) .and. iAt<=transpar%idxdevice(2)) then
      PL = transpar%nPLs
    endif
    if (iAt > transpar%idxdevice(2)) then
      PL = transpar%nPLs + 1
    endif 
  
  end function getPL


end module InitProgram
