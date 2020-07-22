!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!


!> Demonstrates the API with population dependant external potentials.
!>
!> Use it with the input in the test/api/mm/qdepextpot/ folder.
!>
program test_qdepextpot
  use dftbplus
  use extchargepot, only : getPointChargePotential
  use extchargepotgen, only : TExtChargePotGen, TExtChargePotGen_init
  ! Only needed for the internal test system
  use testhelpers, only : writeAutotestTag
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer, parameter :: nQmAtom = 3
  integer, parameter :: nExtCharge = 2

  ! coordinates in atomic units
  real(dp), parameter :: coords(3, nQmAtom) = reshape([&
      & 0.0000000000000000_dp, -1.8897259885789233_dp,  0.0000000000000000_dp,&
      & 0.0000000000000000_dp,  0.0000000000000000_dp,  1.4797763915205659_dp,&
      & 0.0000000000000000_dp,  0.0000000000000000_dp, -1.4797763915205659_dp], [3, 3])

  ! charges in atomic units
  real(dp), parameter :: extChargeCoords(3, nExtCharge) = reshape([&
      & -0.944863438887178_dp, -9.44863438887178_dp, 1.70075418999692_dp,&
      &  4.34637181888102_dp,  -5.85815332110050_dp, 2.64561762888410_dp], [3, 2])

  real(dp), parameter :: extCharges(nExtCharge) = [2.5_dp, -1.9_dp]


  type(TDftbPlus) :: dftbp
  type(TDftbPlusInput) :: input
  type(TExtChargePotGen) :: potGen

  real(dp), allocatable :: extPot(:), extPotGrad(:,:)
  real(dp) :: merminEnergy
  real(dp) :: gradients(3, nQmAtom), grossCharges(nQmAtom)

  character(:), allocatable :: DftbVersion
  integer :: major, minor, patch

  !integer :: devNull

  call getDftbPlusBuild(DftbVersion)
  write(*,*)'DFTB+ build: ' // "'" // trim(DftbVersion) // "'"
  call getDftbPlusApi(major, minor, patch)
  write(*,"(1X,A,1X,I0,'.',I0,'.',I0)")'API version:', major, minor, patch

  ! Pass the 1st external charge to dynamic potential generator from the extchargepotgen module,
  ! while the 2nd charge will be set as constant electrostatic potential
  call TExtChargePotGen_init(potGen, coords, extChargeCoords(:,1:1), extCharges(1:1))

  ! Note: setting the global standard output to /dev/null will also suppress run-time error messages
  !open(newunit=devNull, file="/dev/null", action="write")
  !call TDftbPlus_init(dftbp, outputUnit=devNull)

  ! initialise a calculation then read input from file
  call TDftbPlus_init(dftbp)
  call dftbp%getInputFromFile("dftb_in.hsd", input)
  call dftbp%setupCalculator(input)

  ! set up the above external potential generator for this calculation
  call dftbp%setQDepExtPotGen(potGen)

  ! add an extra fixed external charge
  allocate(extPot(nQmAtom))
  allocate(extPotGrad(3, nQmAtom))
  call getPointChargePotential(extChargeCoords(:,2:2), extCharges(2:2), coords, extPot, extPotGrad)
  call dftbp%setExternalPotential(extPot, extPotGrad)

  ! set the geometry from this program, replacing the dftb_in.hsd values
  call dftbp%setGeometry(coords)

  ! obtain energy and forces
  call dftbp%getEnergy(merminEnergy)
  call dftbp%getGradients(gradients)
  call dftbp%getGrossCharges(grossCharges)
  print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
  print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)
  print "(A,3F15.10)", 'Obtained gross charges:', grossCharges

  ! clean up
  call TDftbPlus_destruct(dftbp)

  ! Write file for internal test system
  call writeAutotestTag(merminEnergy=merminEnergy, gradients=gradients, grossCharges=grossCharges)

end program test_qdepextpot
