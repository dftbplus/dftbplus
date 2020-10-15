!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

program test_extpot
  use, intrinsic :: iso_fortran_env, only : output_unit
  use dftbplus
  use dftbp_constants, only : AA__Bohr
  use extchargepot
  ! Only needed for the internal test system
  use testhelpers, only : writeAutotestTag
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer, parameter :: nAtom = 3

  integer, parameter :: nExtChrg = 2

  ! H2O coordinates (atomic units)
  real(dp), parameter :: initialCoords(3, nAtom) = reshape([&
      & 0.000000000000000E+00_dp, -0.188972598857892E+01_dp,  0.000000000000000E+00_dp,&
      & 0.000000000000000E+00_dp,  0.000000000000000E+00_dp,  0.147977639152057E+01_dp,&
      & 0.000000000000000E+00_dp,  0.000000000000000E+00_dp, -0.147977639152057E+01_dp], [3, nAtom])

  ! H2O atom types
  integer, parameter :: species(nAtom) = [1, 2, 2]

  ! External charges (positions and charges, again atomic units)
  real(dp), parameter :: extCharges(4, nExtChrg) = reshape([&
      &-0.94486343888717805E+00_dp,-0.94486343888717794E+01_dp, 0.17007541899969201E+01_dp, 2.5_dp,&
      & 0.43463718188810203E+01_dp,-0.58581533211004997E+01_dp, 0.26456176288841000E+01_dp, -1.9_dp&
      &], [4, nExtChrg])

  type(TDftbPlus) :: dftbp
  type(TDftbPlusInput) :: input

  real(dp) :: merminEnergy
  real(dp) :: coords(3, nAtom), gradients(3, nAtom), extPot(nAtom), extPotGrad(3, nAtom)
  real(dp) :: atomCharges(nAtom), extChargeGrads(3, nExtChrg)
  type(fnode), pointer :: pRoot, pGeo, pHam, pDftb, pMaxAng, pSlakos, pType2Files, pAnalysis
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

  call dftbp%getEmptyInput(input)
  call input%getRootNode(pRoot)
  call setChild(pRoot, "Geometry", pGeo)
  call setChildValue(pGeo, "Periodic", .false.)
  call setChildValue(pGeo, "TypeNames", ["O", "H"])
  coords(:,:) = 0.0_dp
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
  ! You should provide the skfiles as found in the external/slakos/origin/mio-1-1/ folder. These can
  ! be downloaded with the utils/get_opt_externals script
  call setChild(pDftb, "SlaterKosterFiles", pSlakos)
  call setChild(pSlakos, "Type2FileNames", pType2Files)
  call setChildValue(pType2Files, "Prefix", "./")
  call setChildValue(pType2Files, "Separator", "-")
  call setChildValue(pType2Files, "Suffix", ".skf")

  !  set up analysis options
  call setChild(pRoot, "Analysis", pAnalysis)
  call setChildValue(pAnalysis, "CalculateForces", .true.)

  call setChild(pRoot, "ParserOptions", pParserOpts)
  call setChildValue(pParserOpts, "ParserVersion", 5)
  
  print "(A)", 'Input tree in HSD format:'
  call dumpHsd(input%hsdTree, output_unit)

  ! initialise the DFTB+ calculator
  call dftbp%setupCalculator(input)

  ! Replace coordinates
  coords(:,:) = initialCoords
  call dftbp%setGeometry(coords)

  ! add external point charges
  call getPointChargePotential(extCharges(1:3,:), extCharges(4,:), coords, extPot, extPotGrad)
  call dftbp%setExternalPotential(atomPot=extPot, potGrad=extPotGrad)

  ! get results
  call dftbp%getEnergy(merminEnergy)
  call dftbp%getGradients(gradients)
  call dftbp%getGrossCharges(atomCharges)
  call getPointChargeGradients(coords, atomCharges, extCharges(1:3,:), extCharges(4,:),&
      & extChargeGrads)

  print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
  print "(A,3F15.10)", 'Obtained gross charges:', atomCharges
  print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)
  print "(A,3F15.10)", 'Obtained gradient of atom 2:', gradients(:,2)
  print "(A,3F15.10)", 'Obtained gradient of atom 3:', gradients(:,3)
  print "(A,3F15.10)", 'Obtained gradient of charge 1:', extChargeGrads(:,1)
  print "(A,3F15.10)", 'Obtained gradient of charge 2:', extChargeGrads(:,2)

  call TDftbPlus_destruct(dftbp)

  ! Write file for internal test system
  call writeAutotestTag(merminEnergy=merminEnergy, gradients=gradients, grossCharges=atomCharges,&
      & extChargeGradients=extChargeGrads)

end program test_extpot
