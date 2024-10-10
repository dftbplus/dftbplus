!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

program test_timeprop
  use, intrinsic :: iso_fortran_env, only : output_unit
  use dftbp_common_constants, only : AA__Bohr, V_m__au
  use dftbplus
  ! Only needed for the internal test system
  use testhelpers, only : writeAutotestTag
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer, parameter :: nAtom = 12

  ! H2O coordinates (atomic units)
  real(dp), parameter :: initialCoords(3, nAtom) = reshape([&
       &    0.00000000_dp,      0.41727209_dp,      1.34035331_dp,&
       &    0.00000000_dp,      1.36277581_dp,      0.31264346_dp,&
       &    0.00000000_dp,      0.94549248_dp,     -1.02003049_dp,&
       &    0.00000000_dp,     -0.41727209_dp,     -1.32501204_dp,&
       &    0.00000000_dp,     -1.36277581_dp,     -0.29730219_dp,&
       &    0.00000000_dp,     -0.94549248_dp,      1.03537176_dp,&
       &    0.00000000_dp,      0.74550740_dp,      2.38862741_dp,&
       &    0.00000000_dp,      2.43471932_dp,      0.55254238_dp,&
       &    0.00000000_dp,      1.68922144_dp,     -1.82842183_dp,&
       &    0.00000000_dp,     -0.74550739_dp,     -2.37328614_dp,&
       &    0.00000000_dp,     -2.43471932_dp,     -0.53720110_dp,&
       &    0.00000000_dp,     -1.68922144_dp,      1.84376309_dp], [3, nAtom])

  ! H2O atom types
  integer, parameter :: species(nAtom) = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]

  real(dp), parameter :: timestep = 0.2_dp
  real(dp), parameter :: fstrength = 0.001_dp
  integer, parameter :: nsteps = 1000

  call main_()

contains

  !! Main test routine
  !!
  !! All non-constant variables must be defined here to ensure that they are all explicitely
  !! deallocated before the program finishes (avoiding residual memory that tools like valgrind notice).
  !!
  subroutine main_()

    type(TDftbPlus) :: dftbp
    type(TDftbPlusInput) :: input

    real(dp) :: coords(3, nAtom), merminEnergy, dipole(3, 1), energy, atomNetCharges(nAtom, 1)
    type(fnode), pointer :: pRoot, pGeo, pHam, pDftb, pMaxAng, pSlakos, pType2Files, pElecDyn
    type(fnode), pointer :: pPerturb, pKick

    character(:), allocatable :: DftbVersion
    integer :: major, minor, patch, istep, ii


    call getDftbPlusBuild(DftbVersion)
    write(*,*)'DFTB+ build: ' // "'" // trim(DftbVersion) // "'"
    call getDftbPlusApi(major, minor, patch)
    write(*,"(1X,A,1X,I0,'.',I0,'.',I0)")'API version:', major, minor, patch

    call TDftbPlus_init(dftbp)

    call dftbp%getEmptyInput(input)
    call input%getRootNode(pRoot)
    call setChild(pRoot, "Geometry", pGeo)
    call setChildValue(pGeo, "Periodic", .false.)
    call setChildValue(pGeo, "TypeNames", ["C", "H"])
    coords(:,:) = 0.0_dp
    call setChildValue(pGeo, "TypesAndCoordinates", reshape(species, [1, size(species)]), coords)
    call setChild(pRoot, "Hamiltonian", pHam)
    call setChild(pHam, "Dftb", pDftb)
    call setChildValue(pDftb, "Scc", .true.)
    call setChildValue(pDftb, "SccTolerance", 1e-10_dp)

    call setChild(pDftb, "MaxAngularMomentum", pMaxAng)
    call setChildValue(pMaxAng, "C", "p")
    call setChildValue(pMaxAng, "H", "s")

    call setChild(pDftb, "SlaterKosterFiles", pSlakos)
    call setChild(pSlakos, "Type2FileNames", pType2Files)
    call setChildValue(pType2Files, "Prefix", "./")
    call setChildValue(pType2Files, "Separator", "-")
    call setChildValue(pType2Files, "Suffix", ".skf")

    !  set up electron dynamics options
    call setChild(pRoot, "ElectronDynamics", pElecDyn)
    call setChildValue(pElecDyn, "Steps", nsteps)
    call setChildValue(pElecDyn, "TimeStep", timestep)
    call setChildValue(pElecDyn, "FieldStrength", fstrength*1.0e10_dp*V_m__au)
    call setChild(pElecDyn, "Perturbation", pPerturb)
    call setChild(pPerturb, "Kick", pKick)
    call setChildValue(pKick, "PolarisationDirection", "z")
    ! call setChildValue(pElecDyn, "WriteBondPopulation", .true.)

    print "(A)", 'Input tree in HSD format:'
    call dumpHsd(input%hsdTree, output_unit)

    ! initialise the DFTB+ calculator
    call dftbp%setupCalculator(input)

    ! Replace coordinates
    coords(:,:) = initialCoords*AA__Bohr
    call dftbp%setGeometry(coords)

    ! get ground state
    call dftbp%getEnergy(merminEnergy)

    call dftbp%initializeTimeProp(timestep, .false., .false.)

    do istep = 1, nsteps
      call dftbp%doOneTdStep(istep, dipole=dipole, energy=energy, atomNetCharges=atomNetCharges)
    end do

    print "(A,F15.10)", 'Final SCC Energy:', energy
    print "(A,3F15.10)", 'Final dipole:', (dipole(ii,1), ii=1,3)
    print "(A,100F15.10)", 'Final net atomic charges:', (atomNetCharges(ii,1), ii=1,nAtom)

    ! Write file for internal test system
    call writeAutotestTag(tdEnergy=energy, tdDipole=dipole, tdCharges=atomNetCharges)

  end subroutine main_

end program test_timeprop
