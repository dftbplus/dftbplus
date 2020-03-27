!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains a list of physical constants for the code
module dftbp_constants
  use dftbp_accuracy

  !! Natural constants

  !> pi
  real(dp), parameter :: pi = 3.14159265358979323846_dp

  !> Bohr->Angstrom
  real(dp), parameter :: Bohr__AA = 0.529177249_dp

  !> Angstrom->Bohr
  real(dp), parameter :: AA__Bohr = 1.0_dp / Bohr__AA

  !> Hartre -> eV
  real(dp), parameter :: Hartree__eV = 27.2113845_dp

  !> eV->Hartree
  real(dp), parameter :: eV__Hartree = 1.0_dp / Hartree__eV

  !> Hartree/K
  real(dp), parameter :: Boltzmann = 0.00000316681534524639_dp

  !> atomic mass -> a.m.u.
  real(dp), parameter :: e__amu = 0.00054857990945_dp

  !> a.m.u. -> atomic mass a.u.
  real(dp), parameter :: amu__au = 1.0_dp/ e__amu

  !> a.u. -> femtoseconds
  real(dp), parameter :: au__fs = 0.02418884326505_dp

  !> femtoseconds -> a.u.
  real(dp), parameter :: fs__au = 1.0_dp/au__fs

  !> fine structure constant
  real(dp), parameter :: alpha_fs = 0.007297352568_dp

  !> speed of light in a.u. ~137
  real(dp), parameter :: c = 1.0_dp/alpha_fs

  !> atomic units -> SI
  real(dp), parameter :: au__Coulomb = 1.60217653e-19_dp

  !> SI -> a.u.
  real(dp), parameter :: Coulomb__au = 1.0_dp / au__Coulomb

  !> Hartree -> Joule
  real(dp), parameter :: Hartree__J = 4.3597441775e-18_dp

  !> J -> H
  real(dp), parameter :: J__Hartree = 1.0_dp / Hartree__J

  !> hbar in SI units
  real(dp), parameter :: hbar = 1.054571726e-34_dp

  !> electron g factor
  real(dp), parameter :: gfac = 2.00231930436153_dp

  !> Bohr magneton atomic CGS units
  real(dp), parameter :: mu_B = alpha_fs/2.0_dp


  !> kcal/mol -> eV
  real(dp), parameter :: kcal_mol__eV = 4.33931346011e-2_dp


  !> kcal/mol -> Hartree
  real(dp), parameter :: kcal_mol__Hartree = 0.0015946683859874898_dp


  !> Hartree -> kcal/mol
  real(dp), parameter :: Hartree__kcal_mol = 1.0_dp / kcal_mol__Hartree


  !> Rydberg -> m-1 codata 2006 R_infty
  real(dp), parameter :: Rydberg__m = 10973731.568527


  !> Hartree -> cm-1
  real(dp), parameter :: Hartree__cm = 2.0_dp * Rydberg__m / 100.0_dp


  !> Debye -> a.u.
  real(dp), parameter :: Debye__au = 1.0e-16_dp*alpha_fs*au__fs*Coulomb__au*AA__Bohr**2

  !> a.u. -> Debye
  real(dp), parameter :: au__Debye = 1.0_dp / Debye__au


  !> Conversion factors for Pascals to a.u.
  real(dp), parameter :: pascal__au = J__Hartree / (1e10_dp * AA__Bohr)**3

  !> Conversion factors for a.u. to Pascals
  real(dp), parameter :: au__pascal = 1.0_dp / pascal__au


  !> Conversion factors for V/m to a.u.
  real(dp), parameter :: V_m__au = 1e-10_dp * J__Hartree* au__Coulomb * Bohr__AA

  !>  a.u. to V/m
  real(dp), parameter :: au__V_m = 1.0_dp / V_m__au

  !> Avogadros constant
  real(dp), parameter :: avogadConst = 6.022140857e23_dp

  !> hartree to kj/mol
  real(dp), parameter :: Hartree__kJ_mol = (Hartree__J * avogadConst) / 1000.0_dp

  !> bohr to nm
  real(dp), parameter :: Bohr__nm = Bohr__AA / 10.0_dp

  !> a.u. to ps
  real(dp), parameter :: au__ps = au__fs / 1000.0_dp 

  !> Golden mean plus one
  real(dp), parameter :: goldenMeanP1 = 1.6180339887498949_dp


  !> Highest allowed angular momentum
  integer, parameter :: maxL = 3


  !> Number of shell names
  integer, parameter :: nShellName = maxL + 1


  !> Names of the atomic shells
  character(len=1), parameter :: shellNames(nShellName) = ["s", "p", "d", "f"]


  !> Names of the spin directions
  character(len=4), parameter :: spinName(2) = ["up  ", "down"]

  !> Names of quaternion directions
  character(len=1), parameter :: quaternionName(4) = ["q", "x", "y", "z"]

  !> Imaginary unit
  complex(dp), parameter :: imag = (0.0_dp,1.0_dp)

  !> Symbols of the periodic system of elements, up to 118
  character(len=2), private, parameter :: pse(1:118) = [&
      & 'h ','he',&
      & 'li','be','b ','c ','n ','o ','f ','ne',&
      & 'na','mg','al','si','p ','s ','cl','ar',&
      & 'k ','ca',&
      & 'sc','ti','v ','cr','mn','fe','co','ni','cu','zn',&
      &           'ga','ge','as','se','br','kr',&
      & 'rb','sr',&
      & 'y ','zr','nb','mo','tc','ru','rh','pd','ag','cd',&
      &           'in','sn','sb','te','i ','xe',&
      & 'cs','ba','la',&
      & 'ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb',&
      & 'lu','hf','ta','w ','re','os','ir','pt','au','hg',&
      &           'tl','pb','bi','po','at','rn',&
      & 'fr','ra','ac',&
      & 'th','pa','u ','np','pu','am','cm','bk','cf','es','fm','md','no',&
      & 'lr','rf','db','sg','bh','hs','mt','ds','rg','cn',&
      &           'nh','fl','mc','lv','ts','og' ]

contains


  !> get atomic number from element symbol.
  elemental function symbolToNumber(symbol) result(number)

    !> Element symbol
    character(len=*), intent(in) :: symbol

    !> Atomic number
    integer :: number

    integer, parameter :: offset = iachar('a') - iachar('A')

    character(len=2) :: lcSymbol
    integer :: i, j, k, l

    number = 0
    lcSymbol = '  '

    k = 0
    do j = 1, len_trim(symbol)
      if (k > 2) exit
      l = iachar(symbol(j:j))
      if (k >= 1 .and. l == iachar(' ')) exit
      if (k >= 1 .and. l == 9) exit
      if (l >= iachar('A') .and. l <= iachar('Z')) l = l + offset
      if (l >= iachar('a') .and. l <= iachar('z')) then
        k = k+1
        lcSymbol(k:k) = achar(l)
      end if
    end do

    do i = 1, size(pse)
      if (lcSymbol == pse(i)) then
        number = i
        exit
      end if
    end do

  end function symbolToNumber


end module dftbp_constants
