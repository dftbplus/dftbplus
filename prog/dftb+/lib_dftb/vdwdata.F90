!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains the Van der Waal radii from
!>
!> A. Bondi (1964). "van der Waals Volumes and
!> Radii". J. Phys. Chem. 68:441. doi:10.1021/j100785a001
!>
!> M. Mantina; A.C. Chamberlin; R. Valero; C.J. Cramer; D.G. Truhlar (2009). "Consistent van der
!> Waals Radii for the Whole Main Group". J. Phys. Chem. A. 113:5806–12. doi:10.1021/jp8111556
!>
module vdwdata
  use accuracy, only : dp
  use message, only : error
  use charmanip, only : tolower
  use constants
  implicit none
  private

  public :: getVdwData

  !> Contains van de Waals data (atomic number, chemical symbol, vdW radius)
  type TVdw

    !> Atomic number
    integer :: number

    !> Chemical symbol
    character(2) :: symbol

    !> vdW radii in pm
    integer :: radius

  end type TVdw


  type(TVdw), parameter :: database(118) = [&
      TVdw(1, "H ", 120), TVdw(2, "He", 140), TVdw(3, "Li", 182), TVdw(4, "Be", 153),&
      & TVdw(5, "B ", 192), TVdw(6, "C ", 170), TVdw(7, "N ", 155), TVdw(8, "O ", 152),&
      & TVdw(9, "F ", 147), TVdw(10, "Ne", 154), TVdw(11, "Na", 227), TVdw(12, "Mg", 173),&
      & TVdw(13, "Al", 184), TVdw(14, "Si", 210), TVdw(15, "P ", 180), TVdw(16, "S ", 180),&
      & TVdw(17, "Cl", 175), TVdw(18, "Ar", 188), TVdw(19, "K ", 275), TVdw(20, "Ca", 231),&
      & TVdw(21, "Sc", 211), TVdw(22, "Ti", -1), TVdw(23, "V ", -1), TVdw(24, "Cr", -1),&
      & TVdw(25, "Mn", -1), TVdw(26, "Fe", -1), TVdw(27, "Co", -1), TVdw(28, "Ni", 163),&
      & TVdw(29, "Cu", 140), TVdw(30, "Zn", 139), TVdw(31, "Ga", 187), TVdw(32, "Ge", 211),&
      & TVdw(33, "As", 185), TVdw(34, "Se", 190), TVdw(35, "Br", 185), TVdw(36, "Kr", 202),&
      & TVdw(37, "Rb", 303), TVdw(38, "Sr", 249), TVdw(39, "Y ", -1), TVdw(40, "Zr", -1),&
      & TVdw(41, "Nb", -1), TVdw(42, "Mo", -1), TVdw(43, "Tc", -1), TVdw(44, "Ru", -1),&
      & TVdw(45, "Rh", -1), TVdw(46, "Pd", 163), TVdw(47, "Ag", 172), TVdw(48, "Cd", 158),&
      & TVdw(49, "In", 193), TVdw(50, "Sn", 217), TVdw(51, "Sb", 206), TVdw(52, "Te", 206),&
      & TVdw(53, "I ", 198), TVdw(54, "Xe", 216), TVdw(55, "Cs", 343), TVdw(56, "Ba", 268),&
      & TVdw(57, "La", -1), TVdw(58, "Ce", -1), TVdw(59, "Pr", -1), TVdw(60, "Nd", -1),&
      & TVdw(61, "Pm", -1), TVdw(62, "Sm", -1), TVdw(63, "Eu", -1), TVdw(64, "Gd", -1),&
      & TVdw(65, "Tb", -1), TVdw(66, "Dy", -1), TVdw(67, "Ho", -1), TVdw(68, "Er", -1),&
      & TVdw(69, "Tm", -1), TVdw(70, "Yb", -1), TVdw(71, "Lu", -1), TVdw(72, "Hf", -1),&
      & TVdw(73, "Ta", -1), TVdw(74, "W ", -1), TVdw(75, "Re", -1), TVdw(76, "Os", -1),&
      & TVdw(77, "Ir", -1), TVdw(78, "Pt", 175), TVdw(79, "Au", 166), TVdw(80, "Hg", 155),&
      & TVdw(81, "Tl", 196), TVdw(82, "Pb", 202), TVdw(83, "Bi", 207), TVdw(84, "Po", 197),&
      & TVdw(85, "At", 202), TVdw(86, "Rn", 220), TVdw(87, "Fr", 348), TVdw(88, "Ra", 283),&
      & TVdw(89, "Ac", -1), TVdw(90, "Th", -1), TVdw(91, "Pa", -1), TVdw(92, "U ", 186),&
      & TVdw(93, "Np", -1), TVdw(94, "Pu", -1), TVdw(95, "Am", -1), TVdw(96, "Cm", -1),&
      & TVdw(97, "Bk", -1), TVdw(98, "Cf", -1), TVdw(99, "Es", -1), TVdw(100, "Fm", -1),&
      & TVdw(101, "Md", -1), TVdw(102, "No", -1), TVdw(103, "Lr", -1), TVdw(104, "Rf", -1),&
      & TVdw(105, "Db", -1), TVdw(106, "Sg", -1), TVdw(107, "Bh", -1), TVdw(108, "Hs", -1),&
      & TVdw(109, "Mt", -1), TVdw(110, "Ds", -1), TVdw(111, "Rg", -1), TVdw(112, "Cn", -1),&
      & TVdw(113, "Nh", -1), TVdw(114, "Fl", -1), TVdw(115, "Mc", -1), TVdw(116, "Lv", -1),&
      & TVdw(117, "Ts", -1), TVdw(118, "Og", -1)&
      & ]
  
  contains


  !> Returns distance, atomic number and van de Waals radius.
  subroutine getVdwData(name, radius, number, found)

    !> Name of the element for look for.
    character(len=*), intent(in) :: name

    !> Radius (in Bohr) at return.
    real(dp), intent(out) :: radius

    !> atomic number, if requested
    integer, optional, intent(out) :: number

    !> Flag for whether the element has been found. If this parameter is not specified and the
    !> element is not found, the program stops.
    logical, intent(out), optional :: found

    character(2) :: symbol
    integer :: ii

    symbol = trim(tolower(name))
    do ii = 1, size(database)
      if (trim(tolower(database(ii)%symbol)) == symbol) then
        if (database(ii)%radius > -1) then
          radius = 0.01_dp * real(database(ii)%radius, dp) * AA__Bohr
          if (present(number)) then
            number = database(ii)%number
          end if
          if (present(found)) then
            found = .true.
          end if
          return
        end if
      end if
    end do
    if (present(found)) then
      found = .false.
    else
      call error("VDW database search for element '" // trim(name) // "' failed")
    end if

  end subroutine getVdwData

end module vdwdata
