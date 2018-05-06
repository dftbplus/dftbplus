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
module vdwdata
  use accuracy, only : dp
  use message, only : error
  use charmanip, only : tolower
  use constants
  implicit none
  private

  public :: getVDWData

  !> Contains van de Waals data (atomic number, chemical symbol, vdW radius)
  type TVDW

    !> Atomic number
    integer :: number

    !> Chemical symbol
    character(2) :: symbol

    !> vdW radii in pm
    integer :: radius

  end type TVDW


  type(TVDW) :: database(118) = (/&
      TVDW(1,"H ",120), TVDW(2,"He",140), TVDW(3,"Li",182), TVDW(4,"Be",153), TVDW(5,"B ",192),&
      & TVDW(6,"C ",170), TVDW(7,"N ",155), TVDW(8,"O ",152), TVDW(9,"F ",147), TVDW(10,"Ne",154),&
      & TVDW(11,"Na",227), TVDW(12,"Mg",173), TVDW(13,"Al",184), TVDW(14,"Si",210),&
      & TVDW(15,"P ",180), TVDW(16,"S ",180), TVDW(17,"Cl",175), TVDW(18,"Ar",188),&
      & TVDW(19,"K ",275), TVDW(20,"Ca",231), TVDW(21,"Sc",211), TVDW(22,"Ti",-1),&
      & TVDW(23,"V ",-1), TVDW(24,"Cr",-1), TVDW(25,"Mn",-1), TVDW(26,"Fe",-1), TVDW(27,"Co",-1),&
      & TVDW(28,"Ni",163), TVDW(29,"Cu",140), TVDW(30,"Zn",139), TVDW(31,"Ga",187),&
      & TVDW(32,"Ge",211), TVDW(33,"As",185), TVDW(34,"Se",190), TVDW(35,"Br",185),&
      & TVDW(36,"Kr",202), TVDW(37,"Rb",303), TVDW(38,"Sr",249), TVDW(39,"Y ",-1),&
      & TVDW(40,"Zr",-1), TVDW(41,"Nb",-1), TVDW(42,"Mo",-1), TVDW(43,"Tc",-1), TVDW(44,"Ru",-1),&
      & TVDW(45,"Rh",-1), TVDW(46,"Pd",163), TVDW(47,"Ag",172), TVDW(48,"Cd",158),&
      & TVDW(49,"In",193), TVDW(50,"Sn",217), TVDW(51,"Sb",206), TVDW(52,"Te",206),&
      & TVDW(53,"I ",198), TVDW(54,"Xe",216), TVDW(55,"Cs",343), TVDW(56,"Ba",268),&
      & TVDW(57,"La",-1), TVDW(58,"Ce",-1), TVDW(59,"Pr",-1), TVDW(60,"Nd",-1), TVDW(61,"Pm",-1),&
      & TVDW(62,"Sm",-1), TVDW(63,"Eu",-1), TVDW(64,"Gd",-1), TVDW(65,"Tb",-1), TVDW(66,"Dy",-1),&
      & TVDW(67,"Ho",-1), TVDW(68,"Er",-1), TVDW(69,"Tm",-1), TVDW(70,"Yb",-1),&
      & TVDW(71,"Lu",-1), TVDW(72,"Hf",-1), TVDW(73,"Ta",-1), TVDW(74,"W ",-1), TVDW(75,"Re",-1),&
      & TVDW(76,"Os",-1), TVDW(77,"Ir",-1), TVDW(78,"Pt",175), TVDW(79,"Au",166),&
      & TVDW(80,"Hg",155), TVDW(81,"Tl",196), TVDW(82,"Pb",202), TVDW(83,"Bi",207),&
      & TVDW(84,"Po",197), TVDW(85,"At",202), TVDW(86,"Rn",220), TVDW(87,"Fr",348),&
      & TVDW(88,"Ra",283), TVDW(89,"Ac",-1), TVDW(90,"Th",-1), TVDW(91,"Pa",-1),&
      & TVDW(92,"U ",186), TVDW(93,"Np",-1), TVDW(94,"Pu",-1), TVDW(95,"Am",-1),&
      & TVDW(96,"Cm",-1), TVDW(97,"Bk",-1), TVDW(98,"Cf",-1), TVDW(99,"Es",-1), TVDW(100,"Fm",-1),&
      & TVDW(101,"Md",-1), TVDW(102,"No",-1), TVDW(103,"Lr",-1), TVDW(104,"Rf",-1),&
      & TVDW(105,"Db",-1), TVDW(106,"Sg",-1), TVDW(107,"Bh",-1), TVDW(108,"Hs",-1),&
      & TVDW(109,"Mt",-1), TVDW(110,"Ds",-1), TVDW(111,"Rg",-1), TVDW(112,"Cn",-1),&
      & TVDW(113,"Nh",-1), TVDW(114,"Fl",-1), TVDW(115,"Mc",-1), TVDW(116,"Lv",-1),&
      & TVDW(117,"Ts",-1), TVDW(118,"Og",-1)&
      & /)

  contains


  !> Returns distance, atomic number and van de Waals radius.
  subroutine getVDWData(name, radius, number, found)

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

  end subroutine getVDWData

end module vdwdata
