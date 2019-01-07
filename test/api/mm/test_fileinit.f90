!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Reads in a dftb_in.hsd file then overrides some of its supplied setting, before calculating
!> properties
program test_fileinit
  use dftbplus
  use dftbp_constants, only : AA__Bohr
  implicit none

  ! double precision real
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

  type(TDftbPlus) :: dftbp
  type(TDftbPlusInput) :: input

  real(dp) :: merminEnergy
  real(dp) :: coords(3, 2), latVecs(3, 3), gradients(3, 2)

  integer :: devNull

  ! Note: setting the global standard output to /dev/null to suppress output and run-time error
  ! messages
  open(newunit=devNull, file="/dev/null", action="write")
  call TDftbPlus_init(dftbp, outputUnit=devNull)

  !call TDftbPlus_init(dftbp)

  ! You should provide the dftb_in.hsd and skfile found in the test/prog/dftb+/non-scc/Si_2/ folder
  call dftbp%getInputFromFile("dftb_in.hsd", input)
  call dftbp%setupCalculator(input)

  ! a simple check for data values in the file
  if (dftbp%nrOfAtoms() /= size(coords, dim=2)) then
    print *, "Mismatch between expected number of atoms and dftb_in.hsd file"
    call TDftbPlus_destruct(dftbp)
    stop
  end if

  ! replace the lattice vectors and coordinates in the document tree
  latVecs(:,:) = initialLatVecs
  coords(:,:) = initialCoords
  call dftbp%setGeometry(coords, latVecs)

  ! evaluate energy and forces
  call dftbp%getEnergy(merminEnergy)
  call dftbp%getGradients(gradients)
  print "(A,F15.10)", 'Expected Mermin Energy:', -2.5933460731_dp
  print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
  print "(A,3F15.10)", 'Expected gradient of atom 1:', -0.010321385989_dp, -0.010321385989_dp,&
      & -0.010321385989_dp
  print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)

  ! make a small displacement in the lattice vectors and coordinates
  latVecs(1, 1) = latVecs(1, 1) + 0.1_dp * AA__Bohr
  coords(1, 1) = coords(1, 1) + 0.1_dp * AA__Bohr
  call dftbp%setGeometry(coords, latVecs)

  ! re-calculate energy and forces
  call dftbp%getEnergy(merminEnergy)
  call dftbp%getGradients(gradients)
  print "(A,F15.10)", 'Expected Mermin Energy:', -2.5916977557_dp
  print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
  print "(A,3F15.10)", 'Expected gradient of atom 1:', 0.021290612216_dp, -0.010269102833_dp,&
      & -0.017111497265_dp
  print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)

  ! clean up
  call TDftbPlus_destruct(dftbp)

end program test_fileinit
