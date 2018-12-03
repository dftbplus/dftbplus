!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!


!> Demonstrates the API with population dependant external potentials.
!>
!> Use it with the input in the qdepextpot/ folder.
!>
program test_qdepextpot
  use dftbplus
  use extchargepot, only : getPointChargePotential
  use extchargepotgen, only : TExtChargePotGen, TExtChargePotGen_init
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nQmAtom = 3
  integer, parameter :: nExtCharge = 2

  real(dp), parameter :: coords(3, nQmAtom) = reshape([&
      & 0.0000000000000000_dp, -1.8897259885789233_dp,  0.0000000000000000_dp,&
      & 0.0000000000000000_dp,  0.0000000000000000_dp,  1.4797763915205659_dp,&
      & 0.0000000000000000_dp,  0.0000000000000000_dp, -1.4797763915205659_dp], [3, 3])

  real(dp), parameter :: extChargeCoords(3, nExtCharge) = reshape([&
      & -0.944863438887178_dp, -9.44863438887178_dp, 1.70075418999692_dp,&
      &  4.34637181888102_dp,  -5.85815332110050_dp, 2.64561762888410_dp], [3, 2])

  real(dp), parameter :: extCharges(nExtCharge) = [2.5_dp, -1.9_dp]


  type(TDftbPlus) :: dftbp
  type(TDftbPlusInput) :: input
  type(TExtChargePotGen) :: extChargePotGen

  real(dp), allocatable :: extPot(:), extPotGrad(:,:)
  real(dp) :: merminEnergy
  real(dp) :: gradients(3, 3)
  integer :: devNull

  ! Pass 1st external charge to dynamic potential generator, while 2nd will be set as
  ! constant electrostatic potential
  call TExtChargePotGen_init(extChargePotGen, coords, extChargeCoords(:,1:1), extCharges(1:1))

  open(newunit=devNull, file="/dev/null", action="write")
  !call TDftbPlus_init(dftbp, outputUnit=devNull)
  call TDftbPlus_init(dftbp)
  call dftbp%getInputFromFile("dftb_in.hsd", input)
  call dftbp%setupCalculator(input)
  call dftbp%setQDepExtPotGen(extChargePotGen)

  allocate(extPot(nQmAtom))
  allocate(extPotGrad(3, nQmAtom))
  call getPointChargePotential(extChargeCoords(:,2:2), extCharges(2:2), coords, extPot, extPotGrad)
  call dftbp%setExternalPotential(extPot, extPotGrad)
  call dftbp%setGeometry(coords)
  call dftbp%getEnergy(merminEnergy)
  call dftbp%getGradients(gradients)
  print "(A,F15.10)", 'Expected Mermin Energy:', -3.9854803392_dp
  print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
  print "(A,3F15.10)", 'Expected gradient of atom 1:', 0.017651363773_dp, -0.183137129803_dp,&
      & 0.003198251627_dp
  print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)

  call TDftbPlus_destruct(dftbp)

end program test_qdepextpot
