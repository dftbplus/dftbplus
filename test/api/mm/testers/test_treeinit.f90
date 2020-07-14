!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Test code which builds the DFTB input tree and then evaluates energy and forces. Example is a
!> periodic geometry with k-points.
program test_treeinit
  use, intrinsic :: iso_fortran_env, only : output_unit
  use dftbplus
  use dftbp_constants, only : AA__Bohr
  ! Only needed for the internal test system
  use testhelpers, only : writeAutotestTag
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  ! reference coordinates, (xyz,:nAtom) in atomic units
  real(dp), parameter :: initialCoords(3, 2) = reshape([&
      & 0.0000000000000000_dp, 0.0000000000000000_dp, 0.0000000000000000_dp,&
      & 2.5639291987021915_dp, 2.5639291987021915_dp, 2.5639291987021915_dp], [3, 2])

  ! lattice vectors in atomic units
  real(dp), parameter :: initialLatVecs(3, 3) = reshape([&
      & 5.1278583974043830_dp, 5.1278583974043830_dp, 0.0000000000000000_dp,&
      & 0.0000000000000000_dp, 5.1278583974043830_dp, 5.1278583974043830_dp,&
      & 5.1278583974043830_dp, 0.0000000000000000_dp, 5.1278583974043830_dp], [3, 3])

  ! DFTB+ calculation itself
  type(TDftbPlus) :: dftbp
  ! input settings
  type(TDftbPlusInput) :: input

  real(dp) :: merminEnergy
  real(dp) :: coords(3, 2), latVecs(3, 3), gradients(3, 2), stressTensor(3,3)

  ! pointers to the parts of the input tree that will be set
  type(fnode), pointer :: pRoot, pGeo, pHam, pDftb, pMaxAng, pSlakos, pOptions, pParserOpts

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

  ! You should provide the skfiles found in the external/slakos/origin/pbc-0-3/ folder. These can be
  ! downloaded with the utils/get_opt_externals script.

  ! initialise a DFTB input tree and populate entries which do not have relevant or appropriate
  ! default values
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
  call setChildValue(pSlakos, "Si-Si", "./Si-Si.skf")
  call setChildValue(pDftb, "KPointsAndWeights", reshape([&
      &  0.25_dp,  0.25_dp, 0.25_dp, 1.00_dp,&
      & -0.25_dp,  0.25_dp, 0.25_dp, 1.00_dp,&
      &  0.25_dp, -0.25_dp, 0.25_dp, 1.00_dp,&
      & -0.25_dp, -0.25_dp, 0.25_dp, 1.00_dp], [4, 4]))
  call setChild(pRoot, "Options", pOptions)
  call setChildValue(pOptions, "CalculateForces", .true.)
  call setChild(pRoot, "ParserOptions", pParserOpts)
  call setChildValue(pParserOpts, "ParserVersion", 3)

  ! print resulting input file, including defaults
  print *, 'Input tree in HSD format:'
  call dumpHsd(input%hsdTree, output_unit)
  print *

  ! parse the input for the DFTB+ instance
  call dftbp%setupCalculator(input)

  ! set the lattice vectors and coordinates in the document tree
  latVecs(:,:) = initialLatVecs
  coords(:,:) = initialCoords
  call dftbp%setGeometry(coords, latVecs)
  call dftbp%getEnergy(merminEnergy)
  call dftbp%getGradients(gradients)
  print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
  print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)

  ! make a small displacement in the lattice vectors and coordinates
  latVecs(1, 1) = latVecs(1, 1) + 0.1_dp * AA__Bohr
  coords(1, 1) = coords(1, 1) + 0.1_dp * AA__Bohr
  call dftbp%setGeometry(coords, latVecs)

  ! re-calculate energy and forces
  call dftbp%getEnergy(merminEnergy)
  call dftbp%getGradients(gradients)
  call dftbp%getStressTensor(stressTensor)
  print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
  print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)

  ! clean up
  call TDftbPlus_destruct(dftbp)

  ! Write file for internal test system
  call writeAutotestTag(merminEnergy=merminEnergy, gradients=gradients, stressTensor=stressTensor)

end program test_treeinit
