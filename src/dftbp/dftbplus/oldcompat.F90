!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to convert HSD input for old parser to the current format.
!> Note: parserVersion is set in parser.F90
module dftbp_dftbplus_oldcompat
  use dftbp_common_accuracy, only : dp
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
  subroutine convertOldHSD(root, oldVersion, curVersion)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    !> Version number of the old parser
    integer, intent(in) :: oldVersion

    !> Version number of the current parser
    integer, intent(in) :: curVersion

    integer :: version
    type(fnode), pointer :: ch1, ch2, par

    version = oldVersion
    do while (version < curVersion)
      select case(version)
      case(1)
        call convert_1_2(root)
        version = 2
      case(2)
        call convert_2_3(root)
        version = 3
      case (3)
        call convert_3_4(root)
        version = 4
      case (4)
        call convert_4_5(root)
        version = 5
      case (5)
        call convert_5_6(root)
        version = 6
      case (6)
        call convert_6_7(root)
        version = 7
      case (7)
        call convert_7_8(root)
        version = 8
      case (8)
        call convert_8_9(root)
        version = 9
      case (9)
        call convert_9_10(root)
        version = 10
      case (10)
        call convert_10_11(root)
        version = 11
      case (11)
        call convert_11_12(root)
        version = 12
      end select
    end do

    ! increase the parser version number in the tree - since the resulting dftb_pin would not work
    ! with the old parser as the options have changed to the new parser by now
    call getChildValue(root, "ParserOptions", ch1, "", child=par, &
        &allowEmptyValue=.true.)
    call setChildValue(par, "ParserVersion", version, replace=.true.)

  end subroutine convertOldHSD


  !> Converts input from version 1 to 2. (Version 2 introduced in August 2006)
  subroutine convert_1_2(root)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    type(fnode), pointer :: child1, child2

    call getChild(root, "Geometry", child1, requested=.false.)
    if (associated(child1)) then
      call setUnprocessed(child1)
      call getChild(child1, "SpeciesNames", child2, requested=.false.)
      if (associated(child2)) then
        call setUnprocessed(child2)
        call setNodeName(child2, "TypeNames")
      end if
    end if

  end subroutine convert_1_2


  !> Converts input from version 2 to 3. (Version 3 introduced in Nov. 2006)
  subroutine convert_2_3(root)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    type(fnode), pointer :: ch1, ch2, par
    logical :: tValue

    call getDescendant(root, &
        &"Driver/VelocityVerlet/Thermostat/Andersen/RescalingProbability", &
        &ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword renamed to 'ReselectProbability'.")
      call setNodeName(ch1, "ReselectProbability")
    end if

    call getDescendant(root, &
        &"Driver/VelocityVerlet/Thermostat/Andersen/RescaleIndividually", &
        &ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword renamed to 'ReselectIndividually'.")
      call setNodeName(ch1, "ReselectIndividually")
    end if

    call getDescendant(root, "Hamiltonian/DFTB/Variational", ch1)
    if (associated(ch1)) then
      call getChildValue(ch1, "", tValue)
      call setUnprocessed(ch1)
      if (.not. tValue) then
        call detailedError(ch1, "Sorry, non-variational energy calculation &
            &is not supported any more!")
      else
        call detailedWarning(ch1, "Energy calculation is made only &
            &variational, option removed.")
        call destroyNode(ch1)
      end if
    end if

    call getDescendant(root, "Hamiltonian/DFTB/SCC", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(ch1, "", tValue)
      call setUnprocessed(ch1)
      if (tValue) then
        call setChildValue(par, "OrbitalResolvedSCC", .true., child=ch2)
        call setUnprocessed(ch2)
        call detailedWarning(ch2, "Calculations are not orbital resolved &
            &per default any more. Keyword 'OrbitalResolvedSCC' added.")
      end if
    end if

    call getDescendant(root, "Options/PrintEigenvectors", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'WriteEigenvectors'")
      call setNodeName(ch1, "WriteEigenvectors")
    end if

    call getDescendant(root, "Options/WriteTaggedOut", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'WriteAutotestTag'. &
          &Output file name changed to 'autotest.out'")
      call setNodeName(ch1, "WriteAutotestTag")
    end if

    call getDescendant(root, "Options/WriteBandDat", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'WriteBandOut'. &
          &Output file name changed to 'band.out'")
      call setNodeName(ch1, "WriteBandOut")
    end if

  end subroutine convert_2_3


  !> Converts input from version 3 to 4. (Version 4 introduced in Mar. 2010)
  subroutine convert_3_4(root)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    type(fnode),pointer :: node, node2, node3
    type(fnodeList), pointer :: children
    integer :: ii

    ! Replace range operator with short start:end syntax
    call getDescendant(root, "Driver/SteepestDescent/MovedAtoms", node)
    call replaceRange(node)
    call getDescendant(root, "Driver/ConjugateGradient/MovedAtoms", node)
    call replaceRange(node)
    call getDescendant(root, "Driver/SecondDerivatives/Atoms", node)
    call replaceRange(node)
    call getDescendant(root, "Driver/VelocityVerlet/MovedAtoms", node)
    call replaceRange(node)
    call getDescendant(root, "Hamiltonian/DFTB/SpinPolarisation/Colinear&
        &/InitialSpin", node)
    if (associated(node)) then
      call getChildren(node, "AtomSpin", children)
      do ii = 1, getLength(children)
        call getItem1(children, ii, node2)
        call getChild(node2, "Atoms", node3)
        call replaceRange(node3)
      end do
      call destroyNodeList(children)
    end if

    call getDescendant(root, "Hamiltonian/DFTB/SpinPolarisation/Colinear&
        &/InitialSpin", node)
    if (associated(node)) then
      call detailedWarning(node, "Keyword renamed to 'InitalSpins'.")
      call setNodeName(node, "InitialSpins")
    end if

  end subroutine convert_3_4

  !> Helper function for Range keyword in convert_3_4
  subroutine replaceRange(node)

    !> node to process
    type(fnode), pointer :: node

    type(fnode), pointer :: node2
    integer :: bounds(2)

    if (associated(node)) then
      call getChild(node, "Range", node2, requested=.false.)
      if (associated(node2)) then
        call getChildValue(node2, "", bounds)
        call removeChildNodes(node)
        call setChildValue(node, "", &
            &i2c(bounds(1)) // ":" // i2c(bounds(2)), replace=.true.)
        call detailedWarning(node, "Specification 'Range { start end }' &
            &not supported any more, using 'start:end' instead")
      end if
    end if

  end subroutine replaceRange


  !> Converts input from version 4 to 5. (Version 5 introduced in Dec. 2014)
  subroutine convert_4_5(root)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    type(fnode), pointer :: ch1, ch2, ch3, par, dummy
    logical :: tVal

    call getDescendant(root, "Hamiltonian/DFTB/Eigensolver/Standard", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword renamed to 'QR'.")
      call setNodeName(ch1, "QR")
    end if

    call getDescendant(root, "Options/MullikenAnalysis", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(ch1, "", tVal)
      call detailedWarning(ch1, "Keyword moved to Analysis block.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getChildValue(root, "Analysis", dummy, "", child=ch1, list=.true., &
          & allowEmptyValue=.true., dummyValue=.true.)
      if (.not.associated(ch1)) then
        call setChild(root, "Analysis", ch1)
      end if
      call setChildValue(ch1, "MullikenAnalysis", tVal)
      call setUnprocessed(ch1)
    end if

    call getDescendant(root, "Options/AtomResolvedEnergies", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(par, "AtomResolvedEnergies", tVal)
      call detailedWarning(ch1, "Keyword moved to Analysis block.")
      dummy => removeChild(par,ch1)
      call destroyNode(ch1)
      call getChildValue(root, "Analysis", dummy, "", child=ch1, list=.true., &
          &allowEmptyValue=.true., dummyValue=.true.)
      if (.not.associated(ch1)) then
        call setChild(root, "Analysis", ch1)
      end if
      call setChildValue(ch1,"AtomResolvedEnergies",tVal)
      call setUnprocessed(ch1)
    end if

    call getDescendant(root, "Options/WriteEigenvectors", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(par, "WriteEigenvectors", tVal)
      call detailedWarning(ch1, "Keyword moved to Analysis block.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getChildValue(root, "Analysis", dummy, "", child=ch1, list=.true., &
          &allowEmptyValue=.true., dummyValue=.true.)
      if (.not.associated(ch1)) then
        call setChild(root, "Analysis", ch1)
      end if
      call setChildValue(ch1, "WriteEigenvectors", tVal)
      call setUnprocessed(ch1)
    end if

    call getDescendant(root, "Options/WriteBandOut", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(par, "WriteBandOut", tVal)
      call detailedWarning(ch1, "Keyword moved to Analysis block.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getChildValue(root, "Analysis", dummy, "", child=ch1, list=.true., &
          & allowEmptyValue=.true., dummyValue=.true.)
      if (.not.associated(ch1)) then
        call setChild(root, "Analysis", ch1)
      end if
      call setChildValue(ch1, "WriteBandOut", tVal)
      call setUnprocessed(ch1)
    end if

    call getDescendant(root, "Options/CalculateForces", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(par, "CalculateForces", tVal)
      call detailedWarning(ch1, "Keyword moved to Analysis block.")
      dummy => removeChild(par,ch1)
      call destroyNode(ch1)
      call getChildValue(root, "Analysis", dummy, "", child=ch1, list=.true., &
          &allowEmptyValue=.true., dummyValue=.true.)
      if (.not.associated(ch1)) then
        call setChild(root, "Analysis", ch1)
      end if
      call setChildValue(ch1, "CalculateForces", tVal)
      call setUnprocessed(ch1)
    end if

    call getDescendant(root, "Hamiltonian/DFTB", ch1, parent=par)
    if (associated(ch1)) then
      call setChild(ch1, "Differentiation", ch2)
      call setChild(ch2, "FiniteDiff", ch3)
      call setChildValue(ch3, "Delta", 1.0e-2_dp)
      call detailedWarning(ch2, "Adding legacy step size for finite difference&
          & differentiation")
    end if

    call getDescendant(root, "Hamiltonian/DFTB/SpinConstants", ch1, parent=par)
    if (associated(ch1)) then
      call setChildValue(ch1, "ShellResolvedSpin", .true.)
    end if

  end subroutine convert_4_5

  !> Converts input from version 5 to 6. (Version 6 introduced in May. 2018)
  subroutine convert_5_6(root)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    type(fnode), pointer :: ch1, ch2, ch3, ch4, par, dummy
    logical :: tVal
    real(dp) :: rTmp

    call getDescendant(root, "Analysis/Localise/PipekMezey/Tollerance", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'Tolerance'.")
      call setNodeName(ch1, "Tolerance")
    end if

    call getDescendant(root, "Analysis/Localise/PipekMezey/SparseTollerances", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'SparseTollerances'.")
      call setNodeName(ch1, "SparseTolerances")
    end if

    call getDescendant(root, "Hamiltonian/DFTB/DampXH", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(par, "DampXH", tVal)
      call getDescendant(root, "Hamiltonian/DFTB/DampXHExponent", ch2)
      if (tVal .neqv. associated(ch2)) then
        call error("Incompatible combinaton of DampXH and DampXHExponent")
      end if
      if (associated(ch2)) then
        call getChildValue(par, "DampXHExponent", rTmp)
      end if
      call detailedWarning(ch1, "Keyword DampXH moved to HCorrection block")
      dummy => removeChild(par,ch1)
      call destroyNode(ch1)
      dummy => removeChild(par,ch2)
      call destroyNode(ch2)

      ! clean out any HCorrection entry
      call getDescendant(root, "Hamiltonian/DFTB/HCorrection", ch2, parent=par)
      if (associated(ch2)) then
        call detailedError(ch2, "HCorrection already present.")
      end if

      call getDescendant(root, "Hamiltonian/DFTB", ch2, parent=par)
      call setChild(ch2, "HCorrection", ch3)
      call setChild(ch3, "Damping", ch4)
      call setChildValue(ch4, "Exponent", rTmp)
      call detailedWarning(ch3, "Adding Damping to HCorrection")
    end if

  end subroutine convert_5_6


  !> Converts input from version 6 to 7. (Version 7 introduced in April 2019)
  subroutine convert_6_7(root)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    type(fnode), pointer :: ch1

    call getDescendant(root, "Hamiltonian/DFTB/OrbitalResolvedSCC", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'ShellResolvedSCC'.")
      call setNodeName(ch1, "ShellResolvedSCC")
    end if
    call handleD3Defaults(root)

    call getDescendant(root, "Hamiltonian/DFTB/Eigensolver", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword renamed to 'Solver'.")
      call setNodeName(ch1, "Solver")
    end if

  end subroutine convert_6_7


  !> Converts input from version 7 to 8. (Version 8 introduced in October 2019)
  subroutine convert_7_8(root)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    type(fnode), pointer :: ch1, ch2, par
    logical :: tVal
    type(fnode), pointer :: pTaskType
    type(string) :: buffer

    call getDescendant(root, "Analysis/EigenvectorsAsTxt", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword converted to 'EigenvectorsAsText'.")
      call setNodeName(ch1, "EigenvectorsAsText")
    end if

    call getDescendant(root, "Transport", ch1, parent=par)
    if (associated(ch1)) then
      call getDescendant(ch1, "Task", ch2)
      if (.not. associated(ch2)) then
        call setChildValue(ch1, "readBinaryContact", .false., child=ch2, replace=.true.)
      else
        call getChildValue(ch1, "Task", pTaskType, child=ch2)
        call getNodeName(pTaskType, buffer)
        select case (char(buffer))
        case ("contacthamiltonian")
          call setChildValue(ch1, "writeBinaryContact", .false., child=ch2, replace=.true.)
        case ("uploadcontacts")
          call setChildValue(ch1, "readBinaryContact", .false., child=ch2, replace=.true.)
        end select
      end if
    end if

    call getDescendant(root, "ParserOptions/WriteXMLInput", ch1)
    if (associated(ch1)) then
      call getChildValue(ch1, "", tVal)
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
  subroutine convert_8_9(root)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    type(fnode), pointer :: ch1, ch2, ch3, ch4, par, dummy
    logical :: tVal1, tVal2
    type(fnode), pointer :: pTaskType
    type(string) :: buffer

    ! If this is an electron dynamics restart, then remove keywords for the (un-needed) ground state
    ! calculation (unless the eigenvectors are required)
    call getDescendant(root, "ElectronDynamics/Restart", ch1)
    if (associated(ch1)) then
      call getChildValue(ch1, "", tVal1)
      call setUnprocessed(ch1)
      tVal2 = .false.
      ! Population projection requires eigenvectors, which are not currently stored in the restart
      ! file.
      call getDescendant(root, "ElectronDynamics/Populations", ch2)
      if (associated(ch2)) then
        call getChildValue(ch2, "", tVal2)
        call setUnprocessed(ch2)
      end if
      if (tVal1 .and. .not.tVal2) then
        call getDescendant(root, "Hamiltonian/DFTB/Filling", ch1)
        if (associated(ch1)) then
          call detailedWarning(ch1, "Restarted electronDynamics does not require Filling{}&
              & settings unless projected onto ground state")
          call destroyNode(ch1)
        end if
        call getDescendant(root, "Analysis", ch1)
        if (associated(ch1)) then
          call detailedWarning(ch1, "Restarted electronDynamics does not use the Analysis{} block")
          call destroyNode(ch1)
        end if
      end if
    end if

    call getDescendant(root, "Driver/lBFGS", ch1)
    if (associated(ch1)) then
      call setChildValue(ch1, "LineSearch", .true., child=ch2)
      call setUnprocessed(ch2)
      call detailedWarning(ch2, "Set 'LineSearch = Yes'")
      call setChildValue(ch1, "oldLineSearch", .true., child=ch2)
      call setUnprocessed(ch2)
      call detailedWarning(ch2, "Set 'oldLineSearch = Yes'")
    end if

  end subroutine convert_8_9


  !> Converts input from version 9 to 10. (Version 10 introduced in November 2021)
  subroutine convert_9_10(root)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    type(fnode), pointer :: ch1, ch2, ch3, ch4, par, dummy
    logical :: tVal1, tVal2
    type(fnode), pointer :: pTaskType
    type(string) :: buffer

    call getDescendant(root, "ExcitedState/Casida", ch1)
    if (associated(ch1)) then
      call getChildValue(ch1, "WriteStatusArnoldi", tVal1, default=.false., child=ch2)
      dummy => removeChild(ch1, ch2)
      call getChildValue(ch1, "TestArnoldi", tVal2, default=.false., child=ch2)
      dummy => removeChild(ch1, ch2)
      call detailedWarning(ch1, "Keyword moved to Diagonaliser block.")
      call setUnprocessed(ch1)
      call setChild(ch1, "Diagonaliser", ch2)
      call setUnprocessed(ch2)
      call setChild(ch2, "Arpack", ch3)
      call setUnprocessed(ch3)
      call setChildValue(ch3, "WriteStatusArnoldi", tVal1, child=ch4)
      call setUnprocessed(ch4)
      call setChildValue(ch3, "TestArnoldi", tVal2, child=ch4)
      call setUnprocessed(ch4)
    end if

    ! move ConvergentSccOnly and ConvergentForces into a common keyword

    call getDescendant(root, "Hamiltonian/Dispersion/Ts/ConvergentSCCOnly", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentSCCOnly", tVal1)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian", ch1)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, child=ch2)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Hamiltonian/Dispersion/Mbd/ConvergentSCCOnly", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentSCCOnly", tVal1)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.")
      end if
      call getDescendant(root, "Hamiltonian", ch1)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, child=ch2)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/ConjugateGradient/ConvergentForcesOnly", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.")
      end if
      call getDescendant(root, "Hamiltonian", ch1)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, child=ch2)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/VelocityVerlet/ConvergentForcesOnly", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.")
      end if
      call getDescendant(root, "Hamiltonian", ch1)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, child=ch2)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/SteepestDescent/ConvergentForcesOnly", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.")
      end if
      call getDescendant(root, "Hamiltonian", ch1)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, child=ch2)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/gDiis/ConvergentForcesOnly", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.")
      end if
      call getDescendant(root, "Hamiltonian", ch1)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, child=ch2)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/LBfgs/ConvergentForcesOnly", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.")
      end if
      call getDescendant(root, "Hamiltonian", ch1)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, child=ch2)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/Fire/ConvergentForcesOnly", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.")
      end if
      call getDescendant(root, "Hamiltonian", ch1)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, child=ch2)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/SecondDerivatives/ConvergentForcesOnly", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.")
      end if
      call getDescendant(root, "Hamiltonian", ch1)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, child=ch2)
      call setUnprocessed(ch2)
    end if

    call getDescendant(root, "Driver/Socket/ConvergentForcesOnly", ch1, parent=par)
    if (associated(ch1)) then
      call getChildValue(par, "ConvergentForcesOnly", tVal1)
      call detailedWarning(ch1, "Keyword Moved to Hamiltonian {}.")
      dummy => removeChild(par, ch1)
      call destroyNode(ch1)
      call getDescendant(root, "Hamiltonian/ConvergentSCCOnly", ch3)
      if (associated(ch3)) then
        call detailedError(ch3, "ConvergentSCCOnly already present.")
      end if
      call getDescendant(root, "Hamiltonian", ch1)
      call setChildValue(ch1, "ConvergentSCCOnly", tVal1, child=ch2)
      call setUnprocessed(ch2)
    end if

  end subroutine convert_9_10


  !> Converts input from version 10 to 11. (Version 11 introduced in April 2022)
  subroutine convert_10_11(root)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    type(fnode), pointer :: ch1, ch2, ch3, ch4, par, dummy

    call getDescendant(root, "Hamiltonian/DFTB/Solvation/GeneralizedBorn", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Set solvated field scaling (RescaleSolvatedFields) to No.")
      call setChildValue(ch1, "RescaleSolvatedFields", .false., child=ch2, replace=.true.)
    end if

    call getDescendant(root, "Hamiltonian/DFTB/Solvation/Cosmo", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Set solvated field scaling (RescaleSolvatedFields) to No.")
      call setChildValue(ch1, "RescaleSolvatedFields", .false., child=ch2, replace=.true.)
    end if

    call getDescendant(root, "Hamiltonian/DFTB/Solvation/Sasa", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Set solvated field scaling (RescaleSolvatedFields) to No.")
      call setChildValue(ch1, "RescaleSolvatedFields", .false., child=ch2, replace=.true.)
    end if

  end subroutine convert_10_11


  !> Converts input from version 11 to 12. (Version 12 introduced in June 2022)
  subroutine convert_11_12(root)

    !> Root tag of the HSD-tree
    type(fnode), pointer :: root

    type(fnode), pointer :: ch1, ch2

    call getDescendant(root, "Driver/GeometryOptimization", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword renamed to 'GeometryOptimisation'.")
      call setNodeName(ch1, "GeometryOptimisation")
      call getDescendant(root, "Driver/GeometryOptimisation/Optimizer", ch2)
      if (associated(ch2)) then
        call detailedWarning(ch2, "Keyword renamed to 'Optimiser'.")
        call setNodeName(ch2, "Optimiser")
      end if
    end if

  #:for LABEL in [("Kick"), ("Laser")]
    call getDescendant(root, "ElectronDynamics/Perturbation/${LABEL}$/PolarizationDirection", ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword renamed to 'PolarisationDirection'.")
      call setNodeName(ch1, "PolarisationDirection")
    end if
    call getDescendant(root, "ElectronDynamics/Perturbation/${LABEL}$/ImagPolarizationDirection",&
        & ch1)
    if (associated(ch1)) then
      call detailedWarning(ch1, "Keyword renamed to 'ImagPolarisationDirection'.")
      call setNodeName(ch1, "ImagPolarisationDirection")
    end if
  #:endfor

  end subroutine convert_11_12


  !> Update values in the DftD3 block to match behaviour of v6 parser
  subroutine handleD3Defaults(root)

    !> Root node of the HSD-tree
    type(fnode), pointer :: root

    type(fnode), pointer :: pD3, pDampMethod, pChild
    type(string) :: buffer

    call getDescendant(root, "Hamiltonian/DFTB/Dispersion/DftD3", pD3)
    if (.not. associated(pD3)) then
      return
    end if

    call useDftb3Default(pD3, "s6", 1.0_dp)
    call useDftb3Default(pD3, "s8", 0.5883_dp)

    call getChildValue(pD3, "Damping", pDampMethod, default="BeckeJohnson", child=pChild)
    call setUnprocessed(pChild)
    call setUnprocessed(pDampMethod)
    call getNodeName(pDampMethod, buffer)

    select case (char(buffer))
    case ("beckejohnson")
      call useDftb3Default(pDampMethod, "a1", 0.5719_dp)
      call useDftb3Default(pDampMethod, "a2", 3.6017_dp)
    end select

  end subroutine handleD3Defaults


  !> Helper routine to update values in the DftD3 block to match behaviour of v6 parser
  subroutine useDftb3Default(root, option, default)

    !> Root node of the HSD-tree
    type(fnode), pointer, intent(in) :: root

    !> Name of option inside the DftD3 block
    character(*), intent(in) :: option

    !> Default value to set
    real(dp), intent(in) :: default

    type(fnode), pointer :: pChild

    call getChild(root, option, pChild, requested=.false.)
    if (.not. associated(pChild)) then
      call setChildValue(root, option, default, child=pChild)
      call detailedWarning(pChild, "Using DFTB3 optimised default value for parameter " // option)
    end if
    call setUnprocessed(pChild)

  end subroutine useDftb3Default


end module dftbp_dftbplus_oldcompat
