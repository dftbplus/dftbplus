!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"
module test_wavegrid_initmolorb
  use dftbp_wavegrid, only : TMolecularOrbital, TSlaterOrbital, TSpeciesBasis
  use dftbp_common_accuracy, only : dp
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_boundarycond, only : boundaryCondsEnum, TBoundaryConds, TBoundaryConds_init
  use dftbp_type_typegeometry, only : TGeometry
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: initMolorbH2O, initMolorbHchain

contains

  subroutine initOrbitalHydrogenS(sto, useRadialLut)
    type(TSlaterOrbital), intent(out) :: sto
    logical, intent(in), optional :: useRadialLut

    integer, parameter :: angMom = 0
    real(dp), parameter :: lutResolution = 0.01_dp
    real(dp), parameter :: cutoff = 6.0_dp
    real(dp), parameter :: alpha(3) = [0.5_dp, 1.0_dp, 2.0_dp]
    real(dp), parameter :: aa(3,3) = reshape([ &
        -2.2765228685565400_dp, 0.26641083435260687_dp, -7.9427553748566294E-003_dp, &
         17.453716731738609_dp, -5.4229751699602433_dp, 0.96370929548055750_dp, &
        -12.701455953438341_dp, -6.5568796727250120_dp, -0.85307020704514269_dp], [3,3])
    call sto%init(aa, alpha, angMom, lutResolution, cutoff, useRadialLut=useRadialLut)
  end subroutine initOrbitalHydrogenS


  subroutine initOrbitalOxygenS(sto, useRadialLut)
    type(TSlaterOrbital), intent(out) :: sto
    logical, intent(in), optional :: useRadialLut

    integer, parameter :: angMom = 0
    real(dp), parameter :: lutResolution = 0.01_dp
    real(dp), parameter :: cutoff = 6.0_dp
    real(dp), parameter :: alpha(4) = [0.5_dp, 1.26_dp, 3.17_dp, 8.0_dp]
    real(dp), parameter :: aa(3,4) = reshape([ &
        0.21323488915449521_dp, -0.031152441012403959_dp, 0.0011303933349092960_dp, -9.0596686234211106_dp, &
        4.0675254100925002_dp, -0.60689938592758674_dp, 10.444780695232501_dp, -3.1373217750406419_dp, &
        4.4644627824691749_dp, 8.9215282105208260_dp, 7.2210396633361826_dp, 16.146571373535430_dp], [3,4])
    call sto%init(aa, alpha, angMom, lutResolution, cutoff, useRadialLut=useRadialLut)
  end subroutine initOrbitalOxygenS


  subroutine initOrbitalOxygenP(sto, useRadialLut)
    type(TSlaterOrbital), intent(out) :: sto
    logical, intent(in), optional :: useRadialLut

    integer, parameter :: angMom = 1
    real(dp), parameter :: lutResolution = 0.01_dp
    real(dp), parameter :: cutoff = 6.0_dp
    real(dp), parameter :: alpha(4) = [0.5_dp, 1.26_dp, 3.17_dp, 8.0_dp]
    real(dp), parameter :: aa(3,4) = reshape([ &
        -0.021351405651207991_dp, 0.0028544859270132768_dp, -9.4141846289124166E-005_dp, 1.8517392789336220_dp, &
        -0.79114942586812875_dp, 0.10094277989615121_dp, 16.210706533770320_dp, -10.077615056451849_dp, &
        7.7615980276616314_dp, -1.7017045797631820_dp, -10.773616241206961_dp, -35.439076485248712_dp], [3,4])
    call sto%init(aa, alpha, angMom, lutResolution, cutoff, useRadialLut=useRadialLut)
  end subroutine initOrbitalOxygenP

  subroutine initSpeciesBasisH(speciesBasis, useRadialLut)
    type(TSpeciesBasis), intent(out) :: speciesBasis(1)
    logical, intent(in), optional :: useRadialLut

    ! Hydrogen
    speciesBasis(1)%atomicNumber = 1
    speciesBasis(1)%nOrb = 1
    allocate(speciesBasis(1)%stos(1))
    call initOrbitalHydrogenS(speciesBasis(1)%stos(1), useRadialLut=useRadialLut)

  end subroutine initSpeciesBasisH

  subroutine initSpeciesBasisH2O(speciesBasis, useRadialLut)
    type(TSpeciesBasis), intent(out) :: speciesBasis(2)
    logical, intent(in), optional :: useRadialLut
  
    ! Oxygen
    speciesBasis(1)%atomicNumber = 8
    speciesBasis(1)%nOrb = 2
    allocate(speciesBasis(1)%stos(2))
    call initOrbitalOxygenS(speciesBasis(1)%stos(1), useRadialLut=useRadialLut)
    call initOrbitalOxygenP(speciesBasis(1)%stos(2), useRadialLut=useRadialLut)
    
    ! Hydrogen 
    speciesBasis(2)%atomicNumber = 1
    speciesBasis(2)%nOrb = 1
    allocate(speciesBasis(2)%stos(1))
    call initOrbitalHydrogenS(speciesBasis(2)%stos(1), useRadialLut=useRadialLut)

  end subroutine initSpeciesBasisH2O


  subroutine initGeometryH2O(geometry)
    type(TGeometry), intent(out) :: geometry
    geometry%nAtom = 3
    geometry%tPeriodic = .false.
    geometry%tFracCoord = .false.
    allocate(geometry%species(3))
    geometry%species = [1, 2, 2]
    allocate(geometry%coords(3,3))
    geometry%coords(:,1) = [0.0_dp, -1.8897259885789233_dp, 0.0_dp]
    geometry%coords(:,2) = [0.0_dp, 0.0_dp, 1.4797763915205659_dp]
    geometry%coords(:,3) = [0.0_dp, 0.0_dp, -1.4797763915205659_dp]

    geometry%nSpecies = 2
    allocate(geometry%speciesNames(2))
    geometry%speciesNames = ["O", "H"]
    geometry%tHelical = .false.
  end subroutine initGeometryH2O


  subroutine initGeometryHchain(geometry)
    type(TGeometry), intent(out) :: geometry
    geometry%nAtom = 1
    geometry%tPeriodic = .true.
    geometry%tFracCoord = .false.
    allocate(geometry%species(1))
    geometry%species = [1]
    allocate(geometry%coords(3,1))
    geometry%coords(:,1) = 0.0_dp
    geometry%nSpecies = 1
    allocate(geometry%origin(3))
    geometry%origin = 0.0_dp
    allocate(geometry%latVecs(3,3))
    geometry%latVecs = 0.0_dp
    geometry%latVecs(1,1) = 188.97259885789234_dp
    geometry%latVecs(2,2) = 188.97259885789234_dp
    geometry%latVecs(3,3) = 1.5117807908631387_dp
    allocate(geometry%recVecs2p(3,3))
    geometry%recVecs2p = 0.0_dp
    geometry%recVecs2p(1,1) = 5.2917724899999999E-003_dp
    geometry%recVecs2p(2,2) = 5.2917724899999999E-003_dp
    geometry%recVecs2p(3,3) = 0.66147156125000006_dp
    allocate(geometry%speciesNames(1))
    geometry%speciesNames = ["H"]
    geometry%tHelical = .false.
  end subroutine initGeometryHchain


  
  subroutine initMolorbHchain(molorb, useRadialLut)
    type(TMolecularOrbital), intent(out) :: molorb
    logical, intent(in), optional :: useRadialLut
    
    real(dp), parameter :: gridOrigin(3) = [-5.0_dp, -5.0_dp, -5.0_dp]
    real(dp), parameter :: gridVecs(3,3) = reshape([ &
        0.1_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.1_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.1_dp], [3,3])
    type(TSpeciesBasis) :: speciesBasis(1)
    type(TGeometry) :: geometry
    type(TBoundaryConds) :: bconds
    type(TStatus) :: status

    call initSpeciesBasisH(speciesBasis, useRadialLut=useRadialLut)
    call initGeometryHchain(geometry)
    call TBoundaryConds_init(bconds, boundaryCondsEnum%pbc3d, errStatus=status)
    @:ASSERT(status%isOk())

    call molorb%init(geometry, bconds, speciesBasis, gridOrigin, gridVecs)
  end subroutine initMolorbHchain



  subroutine initMolorbH2O(molorb, useRadialLut)
    type(TMolecularOrbital), intent(out) :: molorb
    logical, intent(in), optional :: useRadialLut
    
    real(dp), parameter :: gridOrigin(3) = [-5.0_dp, -5.0_dp, -5.0_dp]
    real(dp), parameter :: gridVecs(3,3) = reshape([ &
        0.1_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 0.1_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.1_dp], [3,3])
    type(TSpeciesBasis) :: speciesBasis(2)
    type(TGeometry) :: geometry
    type(TBoundaryConds) :: bconds
    type(TStatus) :: status

    call initSpeciesBasisH2O(speciesBasis, useRadialLut=useRadialLut)
    call initGeometryH2O(geometry)
    call TBoundaryConds_init(bconds, boundaryCondsEnum%cluster, errStatus=status)
    @:ASSERT(status%isOk())

    call molorb%init(geometry, bconds, speciesBasis, gridOrigin, gridVecs)
  end subroutine initMolorbH2O

end module test_wavegrid_initmolorb
