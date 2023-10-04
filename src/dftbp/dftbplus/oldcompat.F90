!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Contains routines to convert HSD input for old parser to the current format.
!> Note: parserVersion is set in parser.F90
module dftbp_dftbplus_oldcompat
  use dftbp_common_accuracy, only : dp
  use dftbp_common_status, only : TStatus
  use dftbp_extlibs_xmlf90, only : fnodeList, fnode, removeChild, string, char, getLength,&
      & getNodeName, destroyNode, getItem1, destroyNodeList
  use dftbp_io_charmanip, only : i2c
  use dftbp_io_hsdutils, only : getChildValue, setChildValue, getChild, setChild, detailedWarning,&
      & detailedError, getChildren
  use dftbp_io_hsdutils2, only : getDescendant, setUnprocessed, setNodeName
  use dftbp_io_message, only : error
  use dftbp_io_xmlutils, only : removeChildNodes
  implicit none

  private
  public :: convertOldHSD

contains


  !> Converts an HSD input for an older parser to the current format
  subroutine convertOldHSD(root, oldVersion, curVersion, errStatus)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    !> Version number of the old parser
    integer, intent(in) :: oldVersion

    !> Version number of the current parser
    integer, intent(in) :: curVersion

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    integer :: version
    type(fnode), pointer :: ch1, par

    version = oldVersion
    do while (version < curVersion)
      select case(version)
      case(1)
        call convert_1_2(root, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        version = 2
      case(2)
        call convert_2_3(root, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        version = 3
      case (3)
        call convert_3_4(root, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        version = 4
      case (4)
        call convert_4_5(root, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        version = 5
      case (5)
        call convert_5_6(root, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        version = 6
      case (6)
        call convert_6_7(root, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        version = 7
      case (7)
        call convert_7_8(root, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        version = 8
      case (8)
        call convert_8_9(root, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        version = 9
      case (9)
        call convert_9_10(root, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        version = 10
      case (10)
        call convert_10_11(root, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        version = 11
      case (11)
        call convert_11_12(root, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        version = 12
      case (12)
        call convert_12_13(root, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        version = 13
      case (13)
        call convert_13_14(root, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        version = 14
      end select
    end do

    ! increase the parser version number in the tree - since the resulting dftb_pin would not work
    ! with the old parser as the options have changed to the new parser by now
    call getChildValue(root, "ParserOptions", ch1, errStatus, "", child=par, allowEmptyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    call setChildValue(par, "ParserVersion", version, errStatus, replace=.true.)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine convertOldHSD


  !> Converts input from version 1 to 2. (Version 2 introduced in August 2006)
  subroutine convert_1_2(root, errStatus)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child1, child2

    call getChild(root, "Geometry", child1, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child1)) then
      call setUnprocessed(child1)
      call getChild(child1, "SpeciesNames", child2, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(child2)) then
        call setUnprocessed(child2)
        call setNodeName(child2, "TypeNames")
      end if
    end if

  end subroutine convert_1_2


  !> Converts input from version 2 to 3. (Version 3 introduced in Nov. 2006)
  subroutine convert_2_3(root, errStatus)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: ch1, ch2, par
    logical :: tValue

    call getDescendant(root,&
        & "Driver/VelocityVerlet/Thermostat/Andersen/RescalingProbability",&
        & ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword renamed to 'ReselectProbability'.")
      call setNodeName(ch1, "ReselectProbability")
    end if

    call getDescendant(root,&
        & "Driver/VelocityVerlet/Thermostat/Andersen/RescaleIndividually",&
        & ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword renamed to 'ReselectIndividually'.")
      call setNodeName(ch1, "ReselectIndividually")
    end if

    call getDescendant(root, "Hamiltonian/DFTB/Variational", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(ch1, "", tValue, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch1)
      if (.not. tValue) then
        call detailedError(ch1, "Sorry, non-variational energy calculation is not supported any&
            & more!", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      else
        call detailedWarning(ch1, "Energy calculation is made only variational, option removed.")
        call destroyNode(ch1)
      end if
    end if

    call getDescendant(root, "Hamiltonian/DFTB/SCC", ch1, errStatus, parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(ch1, "", tValue, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch1)
      if (tValue) then
        call setChildValue(par, "OrbitalResolvedSCC", .true., errStatus, child=ch2)
        @:PROPAGATE_ERROR(errStatus)
        call setUnprocessed(ch2)
        call detailedWarning(ch2, "Calculations are not orbital resolved &
            &per default any more. Keyword 'OrbitalResolvedSCC' added.")
      end if
    end if

    call getDescendant(root, "Options/PrintEigenvectors", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'WriteEigenvectors'")
      call setNodeName(ch1, "WriteEigenvectors")
    end if

    call getDescendant(root, "Options/WriteTaggedOut", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'WriteAutotestTag'. &
          &Output file name changed to 'autotest.out'")
      call setNodeName(ch1, "WriteAutotestTag")
    end if

    call getDescendant(root, "Options/WriteBandDat", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'WriteBandOut'. &
          &Output file name changed to 'band.out'")
      call setNodeName(ch1, "WriteBandOut")
    end if

  end subroutine convert_2_3


  !> Converts input from version 3 to 4. (Version 4 introduced in Mar. 2010)
  subroutine convert_3_4(root, errStatus)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: node, node2, node3
    type(fnodeList), pointer :: children
    integer :: ii

    ! Replace range operator with short start:end syntax
    call getDescendant(root, "Driver/SteepestDescent/MovedAtoms", node, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call replaceRange(node, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getDescendant(root, "Driver/ConjugateGradient/MovedAtoms", node, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call replaceRange(node, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getDescendant(root, "Driver/SecondDerivatives/Atoms", node, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call replaceRange(node, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getDescendant(root, "Driver/VelocityVerlet/MovedAtoms", node, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call replaceRange(node, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getDescendant(root, "Hamiltonian/DFTB/SpinPolarisation/Colinear&
        &/InitialSpin", node, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(node)) then
      call getChildren(node, "AtomSpin", children)
      do ii = 1, getLength(children)
        call getItem1(children, ii, node2)
        call getChild(node2, "Atoms", node3, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call replaceRange(node3, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end do
      call destroyNodeList(children)
    end if

    call getDescendant(root, "Hamiltonian/DFTB/SpinPolarisation/Colinear&
        &/InitialSpin", node, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(node)) then
      call detailedWarning(node, "Keyword renamed to 'InitalSpins'.")
      call setNodeName(node, "InitialSpins")
    end if

  end subroutine convert_3_4

  !> Helper function for Range keyword in convert_3_4
  subroutine replaceRange(node, errStatus)

    !> node to process
    type(fnode), pointer :: node

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: node2
    integer :: bounds(2)

    if (associated(node)) then
      call getChild(node, "Range", node2, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(node2)) then
        call getChildValue(node2, "", bounds, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call removeChildNodes(node)
        call setChildValue(node, "", i2c(bounds(1)) // ":" // i2c(bounds(2)), errStatus,&
            & replace=.true.)
        @:PROPAGATE_ERROR(errStatus)
        call detailedWarning(node, "Specification 'Range { start end }' &
            &not supported any more, using 'start:end' instead")
      end if
    end if

  end subroutine replaceRange


  !> Converts input from version 4 to 5. (Version 5 introduced in Dec. 2014)
  subroutine convert_4_5(root, errStatus)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: ch1, ch2, ch3, par, dummy
    logical :: tVal

    call getDescendant(root, "Hamiltonian/DFTB/Eigensolver/Standard", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword renamed to 'QR'.")
      call setNodeName(ch1, "QR")
    end if

    call getDescendant(root, "Options/MullikenAnalysis", ch1, errStatus, parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(ch1, "", tVal, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch1, "Keyword moved to Analysis block.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getChildValue(root, "Analysis", dummy, errStatus, "", child=ch1, list=.true.,&
          & allowEmptyValue=.true., dummyValue=.true.)
      @:PROPAGATE_ERROR(errStatus)
      if (.not.associated(ch1)) then
        call setChild(root, "Analysis", ch1, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call setChildValue(ch1, "MullikenAnalysis", tVal, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch1)
    end if

    call getDescendant(root, "Options/AtomResolvedEnergies", ch1, errStatus, parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(par, "AtomResolvedEnergies", tVal, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch1, "Keyword moved to Analysis block.")
      dummy => removeChild(par,ch1)
      call destroyNode(ch1)
      call getChildValue(root, "Analysis", dummy, errStatus, "", child=ch1, list=.true.,&
          &allowEmptyValue=.true., dummyValue=.true.)
      @:PROPAGATE_ERROR(errStatus)
      if (.not.associated(ch1)) then
        call setChild(root, "Analysis", ch1, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call setChildValue(ch1, "AtomResolvedEnergies", tVal, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch1)
    end if

    call getDescendant(root, "Options/WriteEigenvectors", ch1, errStatus, parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(par, "WriteEigenvectors", tVal, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch1, "Keyword moved to Analysis block.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getChildValue(root, "Analysis", dummy, errStatus, "", child=ch1, list=.true.,&
          &allowEmptyValue=.true., dummyValue=.true.)
      @:PROPAGATE_ERROR(errStatus)
      if (.not.associated(ch1)) then
        call setChild(root, "Analysis", ch1, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call setChildValue(ch1, "WriteEigenvectors", tVal, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch1)
    end if

    call getDescendant(root, "Options/WriteBandOut", ch1, errStatus, parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(par, "WriteBandOut", tVal, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch1, "Keyword moved to Analysis block.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getChildValue(root, "Analysis", dummy, errStatus, "", child=ch1, list=.true.,&
          & allowEmptyValue=.true., dummyValue=.true.)
      @:PROPAGATE_ERROR(errStatus)
      if (.not.associated(ch1)) then
        call setChild(root, "Analysis", ch1, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call setChildValue(ch1, "WriteBandOut", tVal, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch1)
    end if

    call getDescendant(root, "Options/CalculateForces", ch1, errStatus, parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(par, "CalculateForces", tVal, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch1, "Keyword moved to Analysis block.")
      dummy => removeChild(par,ch1)
      call destroyNode(ch1)
      call getChildValue(root, "Analysis", dummy, errStatus, "", child=ch1, list=.true.,&
          &allowEmptyValue=.true., dummyValue=.true.)
      @:PROPAGATE_ERROR(errStatus)
      if (.not.associated(ch1)) then
        call setChild(root, "Analysis", ch1, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call setChildValue(ch1, "CalculateForces", tVal, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch1)
    end if

    call getDescendant(root, "Hamiltonian/DFTB", ch1, errStatus, parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call setChild(ch1, "Differentiation", ch2, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setChild(ch2, "FiniteDiff", ch3, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setChildValue(ch3, "Delta", 1.0e-2_dp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch2, "Adding legacy step size for finite difference&
          & differentiation")
    end if

    call getDescendant(root, "Hamiltonian/DFTB/SpinConstants", ch1, errStatus, parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call setChildValue(ch1, "ShellResolvedSpin", .true., errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

  end subroutine convert_4_5

  !> Converts input from version 5 to 6. (Version 6 introduced in May. 2018)
  subroutine convert_5_6(root, errStatus)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: ch1, ch2, ch3, ch4, par, dummy
    logical :: tVal
    real(dp) :: rTmp

    call getDescendant(root, "Analysis/Localise/PipekMezey/Tollerance", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'Tolerance'.")
      call setNodeName(ch1, "Tolerance")
    end if

    call getDescendant(root, "Analysis/Localise/PipekMezey/SparseTollerances", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'SparseTollerances'.")
      call setNodeName(ch1, "SparseTolerances")
    end if

    call getDescendant(root, "Hamiltonian/DFTB/DampXH", ch1, errStatus, parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(par, "DampXH", tVal, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getDescendant(root, "Hamiltonian/DFTB/DampXHExponent", ch2, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (tVal .neqv. associated(ch2)) then
        call error("Incompatible combinaton of DampXH and DampXHExponent")
      end if
      if (associated(ch2)) then
        call getChildValue(par, "DampXHExponent", rTmp, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call detailedWarning(ch1, "Keyword DampXH moved to HCorrection block")
      dummy => removeChild(par,ch1)
      call destroyNode(ch1)
      dummy => removeChild(par,ch2)
      call destroyNode(ch2)

      ! clean out any HCorrection entry
      call getDescendant(root, "Hamiltonian/DFTB/HCorrection", ch2, errStatus, parent=par)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch2)) then
        call detailedError(ch2, "HCorrection already present.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

      call getDescendant(root, "Hamiltonian/DFTB", ch2, errStatus, parent=par)
      @:PROPAGATE_ERROR(errStatus)
      call setChild(ch2, "HCorrection", ch3, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setChild(ch3, "Damping", ch4, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setChildValue(ch4, "Exponent", rTmp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch3, "Adding Damping to HCorrection")
    end if

  end subroutine convert_5_6


  !> Converts input from version 6 to 7. (Version 7 introduced in April 2019)
  subroutine convert_6_7(root, errStatus)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: ch1

    call getDescendant(root, "Hamiltonian/DFTB/OrbitalResolvedSCC", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'ShellResolvedSCC'.")
      call setNodeName(ch1, "ShellResolvedSCC")
    end if
    call handleD3Defaults(root, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call getDescendant(root, "Hamiltonian/DFTB/Eigensolver", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword renamed to 'Solver'.")
      call setNodeName(ch1, "Solver")
    end if

  end subroutine convert_6_7


  !> Converts input from version 7 to 8. (Version 8 introduced in October 2019)
  subroutine convert_7_8(root, errStatus)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: ch1, ch2, par
    logical :: tVal
    type(fnode), pointer :: pTaskType
    type(string) :: buffer

    call getDescendant(root, "Analysis/EigenvectorsAsTxt", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'EigenvectorsAsText'.")
      call setNodeName(ch1, "EigenvectorsAsText")
    end if

    call getDescendant(root, "Transport", ch1, errStatus, parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getDescendant(ch1, "Task", ch2, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (.not. associated(ch2)) then
        call setChildValue(ch1, "readBinaryContact", .false., errStatus, child=ch2, replace=.true.)
        @:PROPAGATE_ERROR(errStatus)
      else
        call getChildValue(ch1, "Task", pTaskType, errStatus, child=ch2)
        @:PROPAGATE_ERROR(errStatus)
        call getNodeName(pTaskType, buffer)
        select case (char(buffer))
        case ("contacthamiltonian")
          call setChildValue(ch1, "writeBinaryContact", .false., errStatus, child=ch2,&
              & replace=.true.)
          @:PROPAGATE_ERROR(errStatus)
        case ("uploadcontacts")
          call setChildValue(ch1, "readBinaryContact", .false., errStatus, child=ch2,&
              & replace=.true.)
          @:PROPAGATE_ERROR(errStatus)
        end select
      end if
    end if

    call getDescendant(root, "ParserOptions/WriteXMLInput", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(ch1, "", tVal, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch1)
      if (tVal) then
        call detailedWarning(ch1, "Sorry, XML export of the dftb_in.hsd is not supported any more&
            & so is removed")
      else
        call detailedWarning(ch1, "XML export option is removed.")
      end if
      call destroyNode(ch1)
    end if

  end subroutine convert_7_8


  !> Converts input from version 8 to 9. (Version 9 introduced in August 2020)
  subroutine convert_8_9(root, errStatus)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: ch1, ch2
    logical :: tVal1, tVal2

    ! If this is an electron dynamics restart, then remove keywords for the (un-needed) ground state
    ! calculation (unless the eigenvectors are required)
    call getDescendant(root, "ElectronDynamics/Restart", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(ch1, "", tVal1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch1)
      tVal2 = .false.
      ! Population projection requires eigenvectors, which are not currently stored in the restart
      ! file.
      call getDescendant(root, "ElectronDynamics/Populations", ch2, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch2)) then
        call getChildValue(ch2, "", tVal2, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call setUnprocessed(ch2)
      end if
      if (tVal1 .and. .not.tVal2) then
        call getDescendant(root, "Hamiltonian/DFTB/Filling", ch1, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        if (associated(ch1)) then
          call detailedWarning(ch1, "Restarted electronDynamics does not require Filling{}&
              & settings unless projected onto ground state")
          call destroyNode(ch1)
        end if
        call getDescendant(root, "Analysis", ch1, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        if (associated(ch1)) then
          call detailedWarning(ch1, "Restarted electronDynamics does not use the Analysis{} block")
          call destroyNode(ch1)
        end if
      end if
    end if

    call getDescendant(root, "Driver/lBFGS", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call setChildValue(ch1, "LineSearch", .true., errStatus, child=ch2)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch2)
      call detailedWarning(ch2, "Set 'LineSearch = Yes'")
      call setChildValue(ch1, "oldLineSearch", .true., errStatus, child=ch2)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch2)
      call detailedWarning(ch2, "Set 'oldLineSearch = Yes'")
    end if

  end subroutine convert_8_9


  !> Converts input from version 9 to 10. (Version 10 introduced in November 2021)
  subroutine convert_9_10(root, errStatus)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: ch1, ch2, ch3, ch4, par, dummy
    logical :: tVal1, tVal2

    call getDescendant(root, "ExcitedState/Casida", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(ch1, "WriteStatusArnoldi", tVal1, errStatus, default=.false., child=ch2)
      @:PROPAGATE_ERROR(errStatus)
      dummy => removeChild(ch1, ch2)
      call getChildValue(ch1, "TestArnoldi", tVal2, errStatus, default=.false., child=ch2)
      @:PROPAGATE_ERROR(errStatus)
      dummy => removeChild(ch1, ch2)
      call detailedWarning(ch1, "Keyword moved to Diagonaliser block.")
      call setUnprocessed(ch1)
      call setChild(ch1, "Diagonaliser", ch2, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch2)
      call setChild(ch2, "Arpack", ch3, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch3)
      call setChildValue(ch3, "WriteStatusArnoldi", tVal1, errStatus, child=ch4)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch4)
      call setChildValue(ch3, "TestArnoldi", tVal2, errStatus, child=ch4)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch4)
    end if

    ! move ConvergentSccOnly and ConvergentForces into a common keyword

    call getDescendant(root, "Hamiltonian/Dispersion/Ts/ConvergentSCCOnly", ch1, errStatus,&
        & parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentSCCOnly", tVal1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian", ch1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, errStatus, child=ch2)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Hamiltonian/Dispersion/Mbd/ConvergentSCCOnly", ch1, errStatus,&
        & parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentSCCOnly", tVal1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call getDescendant(root, "Hamiltonian", ch1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, errStatus, child=ch2)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/ConjugateGradient/ConvergentForcesOnly", ch1, errStatus,&
        & parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call getDescendant(root, "Hamiltonian", ch1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, errStatus, child=ch2)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/VelocityVerlet/ConvergentForcesOnly", ch1, errStatus,&
        & parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call getDescendant(root, "Hamiltonian", ch1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, errStatus, child=ch2)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/SteepestDescent/ConvergentForcesOnly", ch1, errStatus,&
        & parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call getDescendant(root, "Hamiltonian", ch1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, errStatus, child=ch2)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/gDiis/ConvergentForcesOnly", ch1, errStatus, parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call getDescendant(root, "Hamiltonian", ch1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, errStatus, child=ch2)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/LBfgs/ConvergentForcesOnly", ch1, errStatus, parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call getDescendant(root, "Hamiltonian", ch1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, errStatus, child=ch2)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/Fire/ConvergentForcesOnly", ch1, errStatus, parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call getDescendant(root, "Hamiltonian", ch1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, errStatus, child=ch2)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/SecondDerivatives/ConvergentForcesOnly", ch1, errStatus,&
        & parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call getDescendant(root, "Hamiltonian", ch1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, errStatus, child=ch2)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/Socket/ConvergentForcesOnly", ch1, errStatus, parent=par)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call getDescendant(root, "Hamiltonian", ch1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, errStatus, child=ch2)
      @:PROPAGATE_ERROR(errStatus)
      call setUnprocessed(ch2)
    end if

  end subroutine convert_9_10


  !> Converts input from version 10 to 11. (Version 11 introduced in April 2022)
  subroutine convert_10_11(root, errStatus)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: ch1, ch2

    call getDescendant(root, "Hamiltonian/DFTB/Solvation/GeneralizedBorn", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Set solvated field scaling (RescaleSolvatedFields) to No.")
      call setChildValue(ch1, "RescaleSolvatedFields", .false., errStatus, child=ch2,&
          & replace=.true.)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call getDescendant(root, "Hamiltonian/DFTB/Solvation/Cosmo", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Set solvated field scaling (RescaleSolvatedFields) to No.")
      call setChildValue(ch1, "RescaleSolvatedFields", .false., errStatus, child=ch2,&
          & replace=.true.)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call getDescendant(root, "Hamiltonian/DFTB/Solvation/Sasa", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Set solvated field scaling (RescaleSolvatedFields) to No.")
      call setChildValue(ch1, "RescaleSolvatedFields", .false., errStatus, child=ch2,&
          & replace=.true.)
      @:PROPAGATE_ERROR(errStatus)
    end if

  end subroutine convert_10_11


  !> Converts input from version 11 to 12. (Version 12 introduced in June 2022)
  subroutine convert_11_12(root, errStatus)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: ch1, ch2
    type(string) :: buffer

    call getDescendant(root, "Transport", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call getChildValue(root, "Transport/Task", ch1, errStatus, child=ch2,&
          & default='uploadcontacts')
      @:PROPAGATE_ERROR(errStatus)
      call getNodeName(ch1, buffer)
      if (char(buffer) /= "contacthamiltonian") then
      #:for LABEL in [("xTB"), ("DFTB")]
        call getDescendant(root, "Hamiltonian/${LABEL}$/Charge", ch1, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        if (associated(ch1)) then
          call setUnprocessed(ch1)
          call detailedWarning(ch1, "Device region charge cannot be set if contacts are present.")
          call destroyNode(ch1)
        end if
      #:endfor
      end if
    end if

  end subroutine convert_11_12


  !> Converts input from version 12 to 13. (Version 13 introduced in February 2023)
  subroutine convert_12_13(root, errStatus)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: ch1, ch2, par1
    integer :: maxIter
    logical :: isPerturb, isConvRequired
    real(dp) :: sccTol

    call getDescendant(root, "Analysis/Polarisability", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    isPerturb = associated(ch1)
    if (.not.isPerturb) then
      call getDescendant(root, "Analysis/ResponseKernel", ch1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      isPerturb = associated(ch1)
    end if

    if (isPerturb) then

      call getDescendant(root, "Analysis/Eta", ch1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch1)) then
        call detailedWarning(ch1, "Keyword renamed to 'PerturbEta'.")
        call setNodeName(ch1, "PerturbEta")
      end if

      call getDescendant(root, "Analysis/DegeneracyTolerance", ch1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch1)) then
        call detailedWarning(ch1, "Keyword renamed to 'PerturbDegenTol'.")
        call setNodeName(ch1, "PertubDegenTol")
      end if

      call getDescendant(root, "Hamiltonian/DFTB/MaxSCCIterations", ch1, errStatus, parent=par1)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch1)) then
        call getChildValue(par1, "MaxSCCIterations", maxIter, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call getDescendant(root, "Analysis", ch1, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call setChildValue(ch1, "MaxPerturbIter", maxIter, errStatus, child=ch2)
        @:PROPAGATE_ERROR(errStatus)
        call setUnprocessed(ch2)
      end if

      call getDescendant(root, "Hamiltonian/DFTB/ConvergentSCCOnly", ch1, errStatus, parent=par1)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch1)) then
        call getChildValue(par1, "ConvergentSCCOnly", isConvRequired, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call getDescendant(root, "Analysis", ch1, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call setChildValue(ch1, "ConvergedPerturb", isConvRequired, errStatus, child=ch2)
        @:PROPAGATE_ERROR(errStatus)
        call setUnprocessed(ch2)
      end if

      call getDescendant(root, "Hamiltonian/DFTB/SccTolerance", ch1, errStatus, parent=par1)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(ch1)) then
        call getChildValue(par1, "SccTolerance", sccTol, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call getDescendant(root, "Analysis", ch1, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call setChildValue(ch1, "PerturbSccTol", sccTol, errStatus, child=ch2)
        @:PROPAGATE_ERROR(errStatus)
        call setUnprocessed(ch2)
      end if

    end if

  end subroutine convert_12_13


  !> Converts input from version 13 to 14. (Version 14 introduced in August 2023)
  subroutine convert_13_14(root, errStatus)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: ch1

    call getDescendant(root, "Analysis/CalculateForces", ch1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword renamed to 'PrintForces'.")
      call setNodeName(ch1, "PrintForces")
    end if

  end subroutine convert_13_14

  !> Update values in the DftD3 block to match behaviour of v6 parser
  subroutine handleD3Defaults(root, errStatus)

    !> Root node of the HSD-tree
    type(fnode), pointer :: root

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: pD3, pDampMethod, pChild
    type(string) :: buffer

    call getDescendant(root, "Hamiltonian/DFTB/Dispersion/DftD3", pD3, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (.not. associated(pD3)) then
      return
    end if

    call useDftb3Default(pD3, "s6", 1.0_dp, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call useDftb3Default(pD3, "s8", 0.5883_dp, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(pD3, "Damping", pDampMethod, errStatus, default="BeckeJohnson", child=pChild)
    @:PROPAGATE_ERROR(errStatus)
    call setUnprocessed(pChild)
    call setUnprocessed(pDampMethod)
    call getNodeName(pDampMethod, buffer)

    select case (char(buffer))
    case ("beckejohnson")
      call useDftb3Default(pDampMethod, "a1", 0.5719_dp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call useDftb3Default(pDampMethod, "a2", 3.6017_dp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine handleD3Defaults


  !> Helper routine to update values in the DftD3 block to match behaviour of v6 parser
  subroutine useDftb3Default(root, option, default, errStatus)

    !> Root node of the HSD-tree
    type(fnode), pointer, intent(in) :: root

    !> Name of option inside the DftD3 block
    character(*), intent(in) :: option

    !> Default value to set
    real(dp), intent(in) :: default

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: pChild

    call getChild(root, option, pChild, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (.not. associated(pChild)) then
      call setChildValue(root, option, default, errStatus, child=pChild)
      @:PROPAGATE_ERROR(errStatus)
      call detailedWarning(pChild, "Using DFTB3 optimised default value for parameter " // option)
    end if
    call setUnprocessed(pChild)

  end subroutine useDftb3Default


end module dftbp_dftbplus_oldcompat
