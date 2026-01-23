!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_dftb_coulomb
  use fortuno_serial, only : suite => serial_suite_item, test_list
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_boundarycond, only : boundaryCondsEnum, TBoundaryConds
  use dftbp_dftb_coulomb, only : TCoulombInput, TCoulomb, TCoulomb_init
  use dftbp_dftb_periodic, only : getCellTranslations, TNeighbourList, TNeighbourlist_init,&
      & updateNeighbourListAndSpecies
  use dftbp_math_simplealgebra, only : invert33, determinant33
  use dftbp_type_commontypes, only : TOrbitals
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains


  $:TEST("coulombderivcluster")
    type(TCoulombInput) :: input
    type(TCoulomb) :: coulomb
    type(TEnvironment) :: env
    type(TNeighbourList) :: neighbourList
    type(TOrbitals) :: orb
    type(TStatus) :: errStatus
    integer, parameter :: nAtom = 2, nSpecies = 1, nShell = 1
    integer, parameter :: species0(nAtom) = 1
    ! x aligned pair
    real(dp), parameter :: coords0(3,nAtom) = reshape([0.0,0.0,0.0,1.0,0.0,0.0],[3,nAtom])
    real(dp), parameter :: qAtom(nAtom) = [+1.0_dp, -1.0_dp]
    real(dp), parameter :: rCellVec(3,1) = 0.0_dp
    real(dp), parameter :: cutOff = 2.0_dp
    real(dp), allocatable :: coords(:,:), invR(:, :, :)
    integer, allocatable :: img2centcell(:), iCellVec(:), species(:)
    integer :: nAllAtom, iAt, iCart

    input%ewaldAlpha = 0.0_dp
    input%tolEwald = 1E-9_dp
    input%boundaryCond = boundaryCondsEnum%cluster

    ! Set up indexing for s-orbitals on the atoms
    allocate(orb%nShell(nSpecies))
    orb%nShell(:) = 1
    allocate(orb%nOrbSpecies(nSpecies))
    orb%nOrbSpecies(:) = 1
    allocate(orb%nOrbAtom(nAtom))
    orb%nOrbAtom(:) = 1
    orb%mShell = maxval(orb%nShell)
    orb%mOrb = maxval(orb%nOrbSpecies)
    orb%nOrbAtom(:) = orb%nOrbSpecies(1)
    orb%nOrb = sum(orb%nOrbAtom)
    allocate(orb%angShell(orb%mShell, nSpecies))
    orb%angShell(:,:) = 0
    allocate(orb%iShellOrb(orb%mOrb, nSpecies))
    orb%iShellOrb(:,:) = 1
    allocate(orb%posShell(orb%mShell+1, nSpecies))
    orb%posShell(1,:) = 1
    orb%posShell(2,:) = 2

    call TNeighbourlist_init(neighbourList, nAtom, 1)

    allocate(coords(3,nAtom))
    allocate(img2CentCell(nAtom))
    allocate(iCellVec(nAtom))
    allocate(species(nAtom))

    call updateNeighbourListAndSpecies(env, coords, species, img2CentCell, iCellVec, neighbourList,&
        & nAllAtom, coords0, species0, cutOff, rCellVec, errStatus)
    @:ASSERT(.not.errStatus%hasError())

    call TCoulomb_init(coulomb, input, env, nAtom)
    call coulomb%updateCoords(env, neighbourList, coords, species)

    call coulomb%updateCharges(env, orb, species, qAtom)

    allocate(invR(nAtom, 3, nAtom), source=0.0_dp)
    do iAt = 1, nAtom
      do iCart = 1, 3
        call coulomb%invRPrime(iCart, iAt, invR(:, iCart, iAt))
      end do
    end do
    write(*,*)
    do iAt = 1, nAtom
      write(*,"(2F20.12)")invR
    end do
    @:ASSERT(all(abs(invR - reshape([-1,+1,0,0,0,0,+1,-1,0,0,0,0],[nAtom, 3,nAtom]))&
        & < epsilon(1.0_dp)))

  $:END_TEST()


  $:TEST("coulombderivperiodic")
    type(TBoundaryConds) :: boundaryConds
    type(TCoulombInput) :: input
    type(TCoulomb) :: coulomb
    type(TEnvironment) :: env
    type(TNeighbourList) :: neighbourList
    type(TOrbitals) :: orb
    type(TStatus) :: errStatus
    integer, parameter :: nAtom = 2, nSpecies = 1, nShell = 1
    integer, parameter :: species0(nAtom) = 1
    ! x aligned pair
    !real(dp), parameter :: coords0(3,nAtom) = reshape([0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp],&
    !    & [3,nAtom])
    real(dp), parameter :: coords0(3,nAtom) = reshape([0.0_dp,0.0_dp,0.0_dp,1.1_dp,0.0_dp,0.0_dp],&
        & [3,nAtom])
    real(dp), parameter :: qAtom(nAtom) = [+1, -1], cutOff = 50.0_dp
    real(dp), allocatable :: coords(:,:), invR(:, :, :), cellVecs(:,:), rCellVecs(:,:)
    real(dp), parameter :: latVecs(3,3) = reshape([2,0,0,0,1,0,0,0,1],[3,3])
    integer, allocatable :: img2centcell(:), iCellVec(:), species(:)
    integer :: nAllAtom, iAt, iCart
    real(dp) :: invLatVecs(3,3), recVecs(3,3), vol, V(nAtom)

    input%ewaldAlpha = 0.0_dp
    input%tolEwald = 1E-9_dp
    input%boundaryCond = boundaryCondsEnum%pbc3d

    ! Set up indexing for s-orbitals on the atoms
    allocate(orb%nShell(nSpecies))
    orb%nShell(:) = 1
    allocate(orb%nOrbSpecies(nSpecies))
    orb%nOrbSpecies(:) = 1
    allocate(orb%nOrbAtom(nAtom))
    orb%nOrbAtom(:) = 1
    orb%mShell = maxval(orb%nShell)
    orb%mOrb = maxval(orb%nOrbSpecies)
    orb%nOrbAtom(:) = orb%nOrbSpecies(1)
    orb%nOrb = sum(orb%nOrbAtom)
    allocate(orb%angShell(orb%mShell, nSpecies))
    orb%angShell(:,:) = 0
    allocate(orb%iShellOrb(orb%mOrb, nSpecies))
    orb%iShellOrb(:,:) = 1
    allocate(orb%posShell(orb%mShell+1, nSpecies))
    orb%posShell(1,:) = 1
    orb%posShell(2,:) = 2

    call TNeighbourlist_init(neighbourList, nAtom, 1)

    allocate(coords(3,nAtom))
    allocate(img2CentCell(nAtom))
    allocate(iCellVec(nAtom))
    allocate(species(nAtom))

    invLatVecs(:,:) = latVecs
    call invert33(invLatVecs)
    invLatVecs(:,:) = transpose(invLatVecs)
    recVecs(:,:) = 2.0_dp * pi * invLatVecs
    vol = abs(determinant33(latVecs))
    write(*,*)'Volume', vol

    call TCoulomb_init(coulomb, input, env, nAtom)
    call coulomb%updateLatVecs(latVecs, recVecs, vol)

    call getCellTranslations(cellVecs, rCellVecs, latVecs, invLatVecs, cutOff, boundaryConds)

    call updateNeighbourListAndSpecies(env, coords, species, img2CentCell, iCellVec, neighbourList,&
        & nAllAtom, coords0, species0, cutOff, rCellVecs, errStatus)
    @:ASSERT(.not.errStatus%hasError())

    call coulomb%updateCoords(env, neighbourList, coords0, species0)
    call coulomb%updateCharges(env, orb, species, qAtom)
    call coulomb%updateShifts(env, orb, species, neighbourList%iNeighbour, img2CentCell)
    V(:) = 0.0_dp
    call coulomb%addShiftPerAtom(V)
    write(*,*)'V',V

    allocate(invR(nAtom, 3, nAtom), source=0.0_dp)
    do iAt = 1, nAtom
      do iCart = 1, 3
        call coulomb%invRPrime(env, iCart, iAt, invR(:, iCart, iAt))
      end do
    end do
    write(*,*)
    do iAt = 1, nAtom
      write(*,"(2F20.12)")invR
    end do

    #! @:ASSERT(all(abs(invR - reshape([-1,+1,0,0,0,0,+1,-1,0,0,0,0],[nAtom, 3,nAtom]))&
    !    & < epsilon(1.0_dp)))

    #!@:ASSERT(.false.)

  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("coulomb", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests


end module test_dftb_coulomb
