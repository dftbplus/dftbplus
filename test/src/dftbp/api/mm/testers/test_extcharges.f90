!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> example code for adding external charges to a water molecule calculation
program test_extcharges
  use, intrinsic :: iso_fortran_env, only : output_unit
  use dftbplus, only : convertAtomTypesToSpecies, dumpHsd, fnode, getDftbPlusApi, getDftbPlusBuild,&
      & getMaxAngFromSlakoFile, setChild, setChildValue, TDftbPlus, TDftbPlus_init, TDftbPlusInput
  ! Only needed for the internal test system
  use testhelpers, only : writeAutotestTag
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer, parameter :: nAtom = 3

  !> Number of external charges used in regression case
  integer, parameter :: nExtChrg = 2

  ! H2O coordinates, atomic units
  real(dp), parameter :: initialCoords(3, nAtom) = reshape([&
      & 0.000000000000000E+00_dp, -0.188972598857892E+01_dp,  0.000000000000000E+00_dp,&
      & 0.000000000000000E+00_dp,  0.000000000000000E+00_dp,  0.147977639152057E+01_dp,&
      & 0.000000000000000E+00_dp,  0.000000000000000E+00_dp, -0.147977639152057E+01_dp], [3, nAtom])

  ! Atomic number of each atom
  integer, parameter :: atomTypes(nAtom) = [8, 1, 1]

  ! list of atoms by atomic number
  character(2), parameter :: atomTypeNames(10) = [character(2) ::&
      & "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"]

  ! External charges (positions and charges, again atomic units). Extra charge at end of list to
  ! test reallocation.
  real(dp), parameter :: extCharges(4, nExtChrg + 1) = reshape([&
      &-0.94486343888717805E+00_dp,-0.94486343888717794E+01_dp,0.17007541899969201E+01_dp, 2.5_dp,&
      & 0.43463718188810203E+01_dp,-0.58581533211004997E+01_dp,0.26456176288841000E+01_dp, -1.9_dp,&
      & 0.43463718188810203E+01_dp,-0.58581533211004997E+01_dp,0.26456176288841000E+01_dp, 1.9_dp&
      &], [4, nExtChrg + 1])

  !> (optional variable to API call) sets widths of Gaussians to convolve with external charges
  !> (again set in a.u.). For this particular test, the first two choices are essentially still
  !> point charges:
  real(dp), parameter :: extChargeBlur(nExtChrg + 1) = [0.1_dp, 0.2_dp, 2.0_dp]

  character(100), parameter :: slakoFiles(2, 2) = reshape([character(100) :: &
      & "./O-O.skf", "./H-O.skf", "./O-H.skf", "./H-H.skf"], [2, 2])

  character(1), parameter :: maxAngNames(4) = ["s", "p", "d", "f"]

  call main_()

contains


  !! Main test routine
  !!
  !! All non-constant variables must be defined here to ensure that they are all explicitely
  !! deallocated before the program finishes  (avoiding residual memory that tools like valgrind notice).
  !!
  subroutine main_()

    type(TDftbPlus) :: dftbp
    type(TDftbPlusInput) :: input

    integer, allocatable :: species(:)
    character(2), allocatable :: speciesNames(:)
    real(dp) :: merminEnergy
    real(dp) :: coords(3, nAtom), gradients(3, nAtom), extChargeGrads(3, nExtChrg)
    real(dp) :: atomCharges(nAtom), cm5Charges(nAtom), atomMasses(nAtom)
    type(fnode), pointer :: pRoot, pGeo, pHam, pDftb, pMaxAng, pSlakos, pAnalysis, pCm5
    type(fnode), pointer :: pParserOpts

    character(:), allocatable :: DftbVersion
    integer :: major, minor, patch

    !integer :: devNull

    call getDftbPlusBuild(DftbVersion)
    write(*,*)'DFTB+ build: ' // "'" // trim(DftbVersion) // "'"
    call getDftbPlusApi(major, minor, patch)
    write(*,"(1X,A,1X,I0,'.',I0,'.',I0)")'API version:', major, minor, patch

    ! Note: setting the global standard output to /dev/null will also suppress run-time error messages
    !open(newunit=devNull, file="/dev/null", action="write")
    !call TDftbPlus_init(dftbp, outputUnit=devNull)
    call TDftbPlus_init(dftbp)

    ! You should provide the skfiles as found in the external/slakos/origin/mio-1-1/ folder. These can
    ! be downloaded with the utils/get_opt_externals script
    call dftbp%getEmptyInput(input)
    call input%getRootNode(pRoot)
    call setChild(pRoot, "Geometry", pGeo)
    call setChildValue(pGeo, "Periodic", .false.)

    ! Demonstrates how to convert the atom types if they are not numbered from 1 but use atomic
    ! numbers instead. The atomTypeNames array is optional, if not present, the resulting type names
    ! (which will have to be used at other places) will be X1, X2, etc.
    print "(A)", "Converting atom types"
    call convertAtomTypesToSpecies(atomTypes, species, speciesNames, atomTypeNames)

    call setChildValue(pGeo, "TypeNames", speciesNames)
    coords(:,:) = 0.0_dp
    call setChildValue(pGeo, "TypesAndCoordinates", reshape(species, [1, size(species)]), coords)
    call setChild(pRoot, "Hamiltonian", pHam)
    call setChild(pHam, "Dftb", pDftb)
    call setChildValue(pDftb, "Scc", .true.)
    call setChildValue(pDftb, "SccTolerance", 1e-12_dp)
    call setChild(pDftb, "MaxAngularMomentum", pMaxAng)

    ! read angular momenta from SK data
    call setChildValue(pMaxAng, speciesNames(1),&
        & maxAngNames(getMaxAngFromSlakoFile(slakoFiles(1, 1)) + 1))
    call setChildValue(pMaxAng, speciesNames(2),&
        & maxAngNames(getMaxAngFromSlakoFile(slakoFiles(2, 2)) + 1))

    ! set up locations for SK file data
    call setChild(pDftb, "SlaterKosterFiles", pSlakos)
    call setChildValue(pSlakos, "O-O", trim(slakoFiles(1, 1)))
    call setChildValue(pSlakos, "H-O", trim(slakoFiles(2, 1)))
    call setChildValue(pSlakos, "O-H", trim(slakoFiles(1, 2)))
    call setChildValue(pSlakos, "H-H", trim(slakoFiles(2, 2)))
    call setChild(pRoot, "Analysis", pAnalysis)
    call setChildValue(pAnalysis, "CalculateForces", .true.)
    call setChild(pAnalysis, "CM5", pCm5)
    call setChild(pRoot, "ParserOptions", pParserOpts)
    call setChildValue(pParserOpts, "ParserVersion", 5)

    print "(A)", 'Input tree in HSD format:'
    call dumpHsd(input%hsdTree, output_unit)

    ! convert input into settings for the DFTB+ calculator
    call dftbp%setupCalculator(input)

    ! replace QM atom coordinates
    coords(:,:) = initialCoords
    call dftbp%setGeometry(coords)

    ! add a single external charge
    call dftbp%setExternalCharges(extCharges(:3,:1), extCharges(4,:1), extChargeBlur(:1))

    ! get energy, charges and forces for the single external charge
    call dftbp%getEnergy(merminEnergy)
    call dftbp%getExtChargeGradients(extChargeGrads(:,:1))

    ! replace with a different number of external charges and re-compute
    call dftbp%setExternalCharges(extCharges(:3,:3), extCharges(4,:3), extChargeBlur(:3))

    ! get energy
    call dftbp%getEnergy(merminEnergy)

    ! Finally, external charges corresponding to the regression data
    call dftbp%setExternalCharges(extCharges(:3,:2), extCharges(4,:2), extChargeBlur(:2))

    ! get energy, charges and forces
    call dftbp%getEnergy(merminEnergy)
    call dftbp%getGradients(gradients)
    call dftbp%getExtChargeGradients(extChargeGrads)
    call dftbp%getGrossCharges(atomCharges)
    call dftbp%getCM5Charges(cm5Charges)
    call dftbp%getAtomicMasses(atomMasses)

    print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
    print "(A,3F15.10)", 'Obtained gross charges:', atomCharges
    print "(A,3F15.10)", 'Obtained CM5 charges:', cm5Charges
    print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)
    print "(A,3F15.10)", 'Obtained gradient of atom 2:', gradients(:,2)
    print "(A,3F15.10)", 'Obtained gradient of atom 3:', gradients(:,3)
    print "(A,3F15.10)", 'Obtained gradient of charge 1:', extChargeGrads(:,1)
    print "(A,3F15.10)", 'Obtained gradient of charge 2:', extChargeGrads(:,2)

    ! Write file for internal test system
    call writeAutotestTag(merminEnergy=merminEnergy, gradients=gradients, grossCharges=atomCharges,&
        & extChargeGradients=extChargeGrads, atomMasses=atomMasses, cm5Charges=cm5Charges)

  end subroutine main_

end program test_extcharges
