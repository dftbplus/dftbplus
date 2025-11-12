!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

program test_extefield
  use, intrinsic :: iso_fortran_env, only : output_unit
  use dftbplus, only : dumpHsd, fnode, getDftbPlusApi, getDftbPlusBuild, setChild, setChildValue,&
      & TDftbPlus, TDftbPlus_init, TDftbPlusInput
  ! Only needed for the internal test system
  use testhelpers, only : writeAutotestTag
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nAtom = 3

  ! H2O coordinates (atomic units)
  real(dp), parameter :: coords(3, nAtom) = reshape([&
      & 0.000000000000000E+00_dp, -0.188972598857892E+01_dp,  0.000000000000000E+00_dp,&
      & 0.000000000000000E+00_dp,  0.000000000000000E+00_dp,  0.147977639152057E+01_dp,&
      & 0.000000000000000E+00_dp,  0.000000000000000E+00_dp, -0.147977639152057E+01_dp], [3, nAtom])

  ! H2O atom types
  integer, parameter :: species(nAtom) = [1, 2, 2]

  call main_()

contains

  !! Main test routine
  !!
  !! All non-constant variables must be defined here to ensure that they are all explicitely
  !! deallocated before the program finishes (avoiding residual memory that tools like valgrind
  !! notice).
  !!
  subroutine main_()

    character(:), allocatable :: DftbVersion
    integer :: major, minor, patch, iCart, ii
    real(dp), parameter :: delta = 1.0E-5_dp
    real(dp) :: eFieldVec(3), eFieldStr
    type(TDftbPlus) :: dftbp
    type(TDftbPlusInput) :: input

    real(dp) :: merminEnergy, refMerminEnergy, atomCharges(nAtom), refDipole(3), dipole(3)
    type(fnode), pointer :: pRoot, pGeo, pHam, pDftb, pMaxAng, pSlakos, pType2Files, pAnalysis
    type(fnode), pointer :: pEField, pExtField, pParserOpts

    call getDftbPlusBuild(DftbVersion)
    write(*,*)'DFTB+ build: ' // "'" // trim(DftbVersion) // "'"
    call getDftbPlusApi(major, minor, patch)
    write(*,"(1X,A,1X,I0,'.',I0,'.',I0)")'API version:', major, minor, patch

    ! Note: setting the global standard output to /dev/null will also suppress run-time error
    ! messages
    !open(newunit=devNull, file="/dev/null", action="write")
    !call TDftbPlus_init(dftbp, outputUnit=devNull)
    call TDftbPlus_init(dftbp)

    call dftbp%getEmptyInput(input)
    call input%getRootNode(pRoot)
    call setChild(pRoot, "Geometry", pGeo)
    call setChildValue(pGeo, "Periodic", .false.)
    call setChildValue(pGeo, "TypeNames", ["O", "H"])
    call setChildValue(pGeo, "TypesAndCoordinates", reshape(species, [1, size(species)]), coords)
    call setChild(pRoot, "Hamiltonian", pHam)
    call setChild(pHam, "Dftb", pDftb)
    call setChildValue(pDftb, "Scc", .true.)
    call setChildValue(pDftb, "SccTolerance", 1e-12_dp)

    ! sub-block inside hamiltonian for the maximum angular momenta
    call setChild(pDftb, "MaxAngularMomentum", pMaxAng)
    ! explicitly set the maximum angular momenta for the species
    call setChildValue(pMaxAng, "O", "p")
    call setChildValue(pMaxAng, "H", "s")

    ! get the SK data
    ! You should provide the skfiles as found in the external/slakos/origin/mio-1-1/ folder. These
    ! can be downloaded with the utils/get_opt_externals script
    call setChild(pDftb, "SlaterKosterFiles", pSlakos)
    call setChild(pSlakos, "Type2FileNames", pType2Files)
    call setChildValue(pType2Files, "Prefix", "./")
    call setChildValue(pType2Files, "Separator", "-")
    call setChildValue(pType2Files, "Suffix", ".skf")

    call setChild(pDftb, "ElectricField", pEField)
    call setChild(pEField, "External", pExtField)
    call setChildValue(pExtField, "Direction", [1.0_dp, 0.0_dp, 0.0_dp])
    call setChildValue(pExtField, "Strength", 0.0_dp)

    !  set up analysis options
    call setChild(pRoot, "Analysis", pAnalysis)
    call setChildValue(pAnalysis, "WriteMulliken", .true.)

    call setChild(pRoot, "ParserOptions", pParserOpts)
    call setChildValue(pParserOpts, "ParserVersion", 14)

    print "(A)", 'Input tree in HSD format:'
    call dumpHsd(input%hsdTree, output_unit)

    ! initialise the DFTB+ calculator
    call dftbp%setupCalculator(input)
    call dftbp%setGeometry(coords)

    ! get results
    call dftbp%getEnergy(refMerminEnergy)
    call dftbp%getGrossCharges(atomCharges)
    refDipole(:) = -matmul(coords, atomCharges)

    ! Finite difference dipole moment from external field
    eFieldStr = delta
    do iCart = 1, 3
      dipole(iCart) = 0.0_dp
      do ii = -1, +1, 2
        eFieldVec(:) = 0.0_dp
        eFieldVec(iCart) = real(ii, dp)
        call dftbp%setExternalEfield(eFieldStr, eFieldVec)
        call dftbp%getEnergy(merminEnergy)
        dipole(iCart) = dipole(iCart) + real(ii, dp) * merminEnergy
      end do
    end do
    dipole(:) = dipole / (2.0_dp * delta)

    print "(A,F15.10)", 'Obtained Mermin Energy:', refMerminEnergy
    print "(A,3F15.10)", 'Obtained dipole:', refDipole
    print "(A,3F15.10)", 'Change in dipole wrt E field finite difference:', dipole-refDipole

    ! Write to file for internal test system
    call writeAutotestTag(merminEnergy=refMerminEnergy, groundDipole=dipole-refDipole)

  end subroutine main_

end program test_extefield
