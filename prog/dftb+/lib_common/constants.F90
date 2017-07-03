!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Contains a list of physical constants for the code
module constants
  use accuracy

  !! Natural constants

  real(dp), parameter :: pi = 3.14159265358979323846_dp !* $\pi$
  real(dp), parameter :: Bohr__AA = 0.529177249_dp      !* Bohr->Angstrom
  real(dp), parameter :: AA__Bohr = 1.0_dp / Bohr__AA   !* Angstrom->Bohr
  real(dp), parameter :: Hartree__eV = 27.2113845_dp        !* Hartre -> eV
  real(dp), parameter :: eV__Hartree = 1.0_dp / Hartree__eV !* eV->Hartree
  real(dp), parameter :: Boltzmann = 0.00000316681534524639_dp ! Hartree/K
  real(dp), parameter :: e__amu = 0.00054857990945_dp ! atomic mass -> a.m.u.
  real(dp), parameter :: amu__au = 1.0_dp/ e__amu ! a.m.u. -> atomic mass a.u.
  real(dp), parameter :: au__fs = 0.02418884326505_dp ! a.u. -> femtoseconds
  real(dp), parameter :: fs__au = 1.0_dp/au__fs  ! femtoseconds -> a.u.
  real(dp), parameter :: alpha_fs = 0.007297352568_dp ! fine structure constant
  real(dp), parameter :: c = 1.0_dp/alpha_fs ! speed of light in a.u. $\sim$137
  real(dp), parameter :: au__Coulomb = 1.60217653e-19_dp
  real(dp), parameter :: Coulomb__au = 1.0_dp / au__Coulomb
  real(dp), parameter :: Hartree__J = 4.3597441775e-18_dp ! Hartree -> Joule
  real(dp), parameter :: J__Hartree = 1.0_dp / Hartree__J
  real(dp), parameter :: hbar = 1.054571726e10-34_dp ! hbar in SI units
  real(dp), parameter :: gfac = 2.00231930436153_dp ! electron g factor
  real(dp), parameter :: mu_B = alpha_fs/2.0_dp ! Bohr magneton atomic CGS units

  !!* kcal/mol -> eV
  real(dp), parameter :: kcal_mol__eV = 4.33931346011e-2_dp

  !!* kcal/mol -> Hartree
  real(dp), parameter :: kcal_mol__Hartree = 0.0015946683859874898_dp

  !!* Hartree -> kcal/mol
  real(dp), parameter :: Hartree__kcal_mol = 1.0_dp / kcal_mol__Hartree

  !!* Rydberg -> m-1 codata 2006 R$_\infty$
  real(dp), parameter :: Rydberg__m = 10973731.568527 !

  !!* Hartree -> cm-1
  real(dp), parameter :: Hartree__cm = 2.0_dp * Rydberg__m / 100.0_dp

  !!* Debye -> a.u.
  real(dp), parameter :: Debye__au = 1.0e-16_dp*alpha_fs*au__fs*Coulomb__au*AA__Bohr**2
  real(dp), parameter :: au__Debye = 1.0_dp / Debye__au

  ! Conversion factors for Pascal
  real(dp), parameter :: pascal__au = J__Hartree / (1e10_dp * AA__Bohr)**3
  real(dp), parameter :: au__pascal = 1.0_dp / pascal__au

  ! Conversion factors for V/m
  real(dp), parameter :: V_m__au = 1e-10_dp * J__Hartree* au__Coulomb * Bohr__AA
  real(dp), parameter :: au__V_m = 1.0_dp / V_m__au

  !!* Golden mean plus one
  real(dp), parameter :: goldenMeanP1 = 1.6180339887498949_dp

  !!* Highest allowed angular momentum
  integer, parameter :: maxL = 3

  !!* Number of orbital names
  integer, parameter :: nOrbitalName = maxL + 1

  !!* Name of the orbitals
  character(len=1), parameter :: orbitalNames(nOrbitalName) = &
      &(/ "s", "p", "d", "f" /)

  !!* Names of the spin directions
  character(len=4), parameter :: spinName(2) = &
      &(/ "up  ", "down" /)
  character(len=1), parameter :: quaternionName(4) = &
      &(/ "q", "x", "y", "z" /)

end module constants
