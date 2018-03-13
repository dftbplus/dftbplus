!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

program test_cluster
  use, intrinsic :: iso_fortran_env, only : output_unit
  use dftbplus
  use dftbp_constants, only : AA__Bohr
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  
  real(dp), parameter :: initialCoords(3, 2) = reshape([&
      & 0.0000000000000000_dp, 0.0000000000000000_dp, 0.0000000000000000_dp,&
      & 2.5639291987021915_dp, 2.5639291987021915_dp, 2.5639291987021915_dp], [3, 2])

  real(dp), parameter :: initialLatVecs(3, 3) = reshape([&
      & 5.1278583974043830_dp, 5.1278583974043830_dp, 0.0000000000000000_dp,&
      & 0.0000000000000000_dp, 5.1278583974043830_dp, 5.1278583974043830_dp,&
      & 5.1278583974043830_dp, 0.0000000000000000_dp, 5.1278583974043830_dp], [3, 3])

  type(TDftbPlus) :: dftbp
  type(TDftbPlusInput) :: input

  real(dp) :: merminEnergy
  real(dp) :: coords(3, 2), latVecs(3, 3), gradients(3, 2)
  integer :: devNull
  type(fnode), pointer :: pRoot, pGeo, pHam, pDftb, pMaxAng, pSlakos, pOptions, pParserOpts

  ! Note: setting the global standard output to /dev/null will also suppress run-time error messages
  !open(newunit=devNull, file="/dev/null", action="write")
  !call TDftbPlus_init(dftbp, outputUnit=devNull)
  call TDftbPlus_init(dftbp)

  ! You should provide the dftb_in.hsd and skfiles as found in the
  ! test/prog/dftb+/non-scc/Si_2/ folder
  call dftbp%getEmptyInput(input)
  call input%getRootNode(pRoot)
  call setChild(pRoot, "Geometry", pGeo)
  call setChildValue(pGeo, "Periodic", .true.)
  call setChildValue(pGeo, "LatticeVectors", initialLatVecs)
  call setChildValue(pGeo, "TypeNames", ["Si"])
  coords(:,:) = 0.0_dp
  call setChildValue(pGeo, "TypesAndCoordinates", reshape([1, 1], [1, 2]), coords)
  call setChild(pRoot, "Hamiltonian", pHam)
  call setChild(pHam, "Dftb", pDftb)
  call setChildValue(pDftb, "Scc", .false.)
  call setChild(pDftb, "MaxAngularMomentum", pMaxAng)
  call setChildValue(pMaxAng, "Si", "p")
  call setChild(pDftb, "SlaterKosterFiles", pSlakos)
  call setChildValue(pSlakos, "Si-Si", "./external/slakos/origin/pbc-0-3/Si-Si.skf")
  call setChildValue(pDftb, "KPointsAndWeights", reshape([&
      &  0.25_dp,  0.25_dp, 0.25_dp, 1.00_dp,&
      & -0.25_dp,  0.25_dp, 0.25_dp, 1.00_dp,&
      &  0.25_dp, -0.25_dp, 0.25_dp, 1.00_dp,&
      & -0.25_dp, -0.25_dp, 0.25_dp, 1.00_dp], [4, 4]))
  call setChild(pRoot, "Options", pOptions)
  call setChildValue(pOptions, "CalculateForces", .true.)
  call setChild(pRoot, "ParserOptions", pParserOpts)
  call setChildValue(pParserOpts, "ParserVersion", 3)
  
  print *, 'Input tree in HSD format:'
  call dumpHsd(input%hsdTree, output_unit)
  

  call dftbp%setupCalculator(input)

  latVecs(:,:) = initialLatVecs
  coords(:,:) = initialCoords
  call dftbp%setGeometry(coords, latVecs)
  call dftbp%getEnergy(merminEnergy)
  call dftbp%getGradients(gradients)
  print "(A,F15.10)", 'Expected Mermin Energy:', -2.5933460731_dp
  print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
  print "(A,3F15.10)", 'Expected gradient of atom 1:', -0.010321385989_dp, -0.010321385989_dp,&
      & -0.010321385989_dp
  print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)

  latVecs(1, 1) = latVecs(1, 1) + 0.1_dp * AA__Bohr
  coords(1, 1) = coords(1, 1) + 0.1_dp * AA__Bohr
  call dftbp%setGeometry(coords, latVecs)
  call dftbp%getEnergy(merminEnergy)
  call dftbp%getGradients(gradients)
  print "(A,F15.10)", 'Expected Mermin Energy:', -2.5916977557_dp
  print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
  print "(A,3F15.10)", 'Expected gradient of atom 1:', 0.021290612216_dp, -0.010269102833_dp,&
      & -0.017111497265_dp
  print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)
  
  call TDftbPlus_destruct(dftbp)
  
end program test_cluster
