!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

#:set UNIT_CONVERTER_RANKS = [0, 1, 2]

!> Contains names of various units and their conversion factors to the corresponding unit used
!> internal in the code (atomic units).
module dftbp_common_unitconversion
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : AA__Bohr, Bohr__AA, J__Hartree, fs__au, c, Coulomb__au,&
      & V_m__au, pascal__au, kg__au, pi, Hartree__cm, eV__Hartree, au__fs, Debye__au, hbar,&
      & amu__au, Boltzmann, kcal_mol__Hartree, nm__Bohr
  use dftbp_io_charmanip, only : tolower
  implicit none

  private
  public :: TUnit
  public :: lengthUnits, inverseLengthUnits, energyUnits, forceUnits, timeUnits, freqUnits
  public :: volumeUnits, chargeUnits, EFieldUnits, BFieldUnits, pressureUnits, velocityUnits
  public :: dipoleUnits, massUnits, angularUnits, massDensityUnits
  public :: statusCodes, convertUnit


  !> Internal parameter for specifying the maximum unit name length
  integer, parameter :: maxUnitNameLen = 20


  !> Possible status codes emitted by public routines of this module
  type :: TStatusCodesEnum

    !> Everything OK
    integer :: ok = 0

    !> Unit name could not found in the unit array
    integer :: unitNotFound = 1

  end type TStatusCodesEnum

  !> Status code enum instance
  type(TStatusCodesEnum), parameter :: statusCodes = TStatusCodesEnum()


  !> Contains name of a unit and its conversion factor
  type TUnit
    private

    !> Name of the conversion factor
    character(maxUnitNameLen) :: name

    !> Conversion factor
    real(dp) :: conversionFact = 1.0_dp

    !> Whether the value should be inverted before multiplied with conversion factor
    logical :: invertValue = .false.

  contains

    procedure :: writeFormatted => TUnit_writeFormatted
    generic :: write(formatted) => writeFormatted

  end type TUnit


  interface convertUnit
    #:for RANK in UNIT_CONVERTER_RANKS
      module procedure convertUnitR${RANK}$
    #:endfor
  end interface convertUnit


  !> Length units
  type(TUnit), parameter :: lengthUnits(*) = [&
      & TUnit("angstrom            ", AA__Bohr),&
      & TUnit("aa                  ", AA__Bohr),&
      & TUnit("a                   ", AA__Bohr),&
      & TUnit("meter               ", 1.0e10_dp * AA__Bohr),&
      & TUnit("m                   ", 1.0e10_dp * AA__Bohr),&
      & TUnit("bohr                ", 1.0_dp),&
      & TUnit("pm                  ", 1.0e-2_dp * AA__Bohr),&
      & TUnit("picometer           ", 1.0e-2_dp * AA__Bohr),&
      & TUnit("au                  ", 1.0_dp)&
      & ]


  !> Inverse length units
  type(TUnit), parameter :: inverseLengthUnits(*) = [&
      & TUnit("1/angstrom          ", Bohr__AA),&
      & TUnit("1/aa                ", Bohr__AA),&
      & TUnit("1/a                 ", Bohr__AA),&
      & TUnit("1/meter             ", 1.0e-10_dp * Bohr__AA),&
      & TUnit("1/m                 ", 1.0e-10_dp * Bohr__AA),&
      & TUnit("1/bohr              ", 1.0_dp),&
      & TUnit("1/pm                ", 1.0e+2_dp * Bohr__AA),&
      & TUnit("1/picometer         ", 1.0e+2_dp * Bohr__AA),&
      & TUnit("1/au                ", 1.0_dp)&
      & ]


  !> Energy units
  type(TUnit), parameter :: energyUnits(*) = [&
      & TUnit("rydberg             ", 0.5_dp),&
      & TUnit("ry                  ", 0.5_dp),&
      & TUnit("electronvolt        ", eV__Hartree),&
      & TUnit("ev                  ", eV__Hartree),&
      & TUnit("kcal/mol            ", kcal_mol__Hartree),&
      & TUnit("kelvin              ", Boltzmann),&
      & TUnit("k                   ", Boltzmann),&
      & TUnit("cm^-1               ", 1.0_dp/Hartree__cm),&
      & TUnit("hartree             ", 1.0_dp),&
      & TUnit("ha                  ", 1.0_dp),&
      & TUnit("joule               ", J__Hartree),&
      & TUnit("j                   ", J__Hartree),&
      & TUnit("au                  ", 1.0_dp),&
      & TUnit("nm                  ", 2.0_dp * pi * c / nm__Bohr, invertValue=.true.)&
      & ]


  !> Force units
  type(TUnit), parameter :: forceUnits(*) = [&
      & TUnit("ev/angstrom         ", eV__Hartree / AA__Bohr),&
      & TUnit("ev/aa               ", eV__Hartree / AA__Bohr),&
      & TUnit("ev/a                ", eV__Hartree / AA__Bohr),&
      & TUnit("ha/bohr             ", 1.0_dp),&
      & TUnit("hartree/bohr        ", 1.0_dp),&
      & TUnit("j/m                 ", J__Hartree / (1.0e10_dp * AA__Bohr)),&
      & TUnit("joule/m             ", J__Hartree / (1.0e10_dp * AA__Bohr)),&
      & TUnit("au                  ", 1.0_dp)&
      & ]


  !> Time units
  type(TUnit), parameter :: timeUnits(*) = [&
      & TUnit("femtosecond         ", fs__au),&
      & TUnit("fs                  ", fs__au),&
      & TUnit("picosecond          ", 1e3_dp * fs__au),&
      & TUnit("ps                  ", 1e3_dp * fs__au),&
      & TUnit("second              ", 1e15_dp * fs__au),&
      & TUnit("s                   ", 1e15_dp * fs__au),&
      & TUnit("au                  ", 1.0_dp)&
      & ]


  !> Frequency units
  type(TUnit), parameter :: freqUnits(*) = [&
      & TUnit("au                  ", 1.0_dp),&
      & TUnit("hz                  ", 1e-15_dp * au__fs),&
      & TUnit("thz                 ", 1e-3_dp * au__fs),&
      & TUnit("cm^-1               ", 1e-8_dp * Bohr__AA * c),&
      & TUnit("ev                  ", eV__Hartree) &
      & ]


  !> Volume units
  type(TUnit), parameter :: volumeUnits(*) = [&
      & TUnit("angstrom^3          ", AA__Bohr**3),&
      & TUnit("aa^3                ", AA__Bohr**3),&
      & TUnit("a^3                 ", AA__Bohr**3),&
      & TUnit("meter^3             ", (1.0e10_dp * AA__Bohr)**3),&
      & TUnit("m^3                 ", (1.0e10_dp * AA__Bohr)**3),&
      & TUnit("picometer^3         ", (1.0e-2_dp * AA__Bohr)**3),&
      & TUnit("pm^3                ", (1.0e-2_dp * AA__Bohr)**3),&
      & TUnit("bohr^3              ", 1.0_dp),&
      & TUnit("au                  ", 1.0_dp)&
      & ]


  !> Charge units
  type(TUnit), parameter :: chargeUnits(*) = [&
      & TUnit("coulomb             ", Coulomb__au),&
      & TUnit("c                   ", Coulomb__au),&
      & TUnit("e                   ", 1.0_dp),&
      & TUnit("au                  ", 1.0_dp)&
      & ]


  !> Electric field units
  type(TUnit), parameter :: EFieldUnits(*) = [&
       & TUnit("v/m                 ", V_m__au),&
       & TUnit("v/a                 ", 1e10_dp * V_m__au),&
       & TUnit("v/aa                ", 1e10_dp * V_m__au),&
       & TUnit("v/angstrom          ", 1e10_dp * V_m__au),&
       & TUnit("au                  ", 1.0_dp)&
       & ]


  !> Magnetic field units (Atomic "Gaussian" CGS unit system!)
  type(TUnit), parameter :: BFieldUnits(*) = [&
      & TUnit("t                 ", 1.0E+24_dp*c/(hbar*Coulomb__au*AA__Bohr**2)),&
      & TUnit("tesla             ", 1.0E+24_dp*c/(hbar*Coulomb__au*AA__Bohr**2)),&
      & TUnit("gauss             ", 1.0E+20_dp*c/(hbar*Coulomb__au*AA__Bohr**2)),&
      & TUnit("g                 ", 1.0E+20_dp*c/(hbar*Coulomb__au*AA__Bohr**2)),&
      & TUnit("au                  ", 1.0_dp)&
      & ]


  !> Pressure units
  type(TUnit), parameter :: pressureUnits(*) = [&
      & TUnit("pa                  ", pascal__au),&
      & TUnit("au                  ", 1.0_dp)&
      & ]


  !> Velocity units
  type(TUnit), parameter :: velocityUnits(*) = [&
      & TUnit("au                  ", 1.0_dp),&
      & TUnit("m/s                 ", 1e10_dp * AA__Bohr / (1e15_dp * fs__au)),&
      & TUnit("pm/fs               ", 1e-2_dp * AA__Bohr / fs__au),&
      & TUnit("a/ps                ", AA__Bohr / (1e3_dp * fs__au)),&
      & TUnit("aa/ps               ", AA__Bohr / (1e3_dp * fs__au)),&
      & TUnit("angstrom/ps         ", AA__Bohr / (1e3_dp * fs__au))&
      & ]


  !> Dipole units
  type(TUnit), parameter :: dipoleUnits(*) = [&
      & TUnit("au                  ", 1.0_dp),&
      & TUnit("debye               ", Debye__au),&
      & TUnit("cm                  ", Coulomb__au * 1.0e10_dp * AA__Bohr),&
      & TUnit("coulombmeter        ", Coulomb__au * 1.0e10_dp * AA__Bohr)&
      & ]


  !> Mass units
  type(TUnit), parameter :: massUnits(*) = [&
      & TUnit("au                  ", 1.0_dp),&
      & TUnit("amu                 ", amu__au),&
      & TUnit("da                  ", amu__au),&
      & TUnit("dalton              ", amu__au),&
      & TUnit("kg                  ", kg__au),&
      & TUnit("g                   ", 1.0e+3_dp * kg__au)&
      & ]


  !> Angular units
  type(TUnit), parameter :: angularUnits(*) = [&
      & TUnit("degrees             ", pi / 180.0_dp),&
      & TUnit("deg                 ", pi / 180.0_dp),&
      & TUnit("radian              ", 1.0_dp),&
      & TUnit("rad                 ", 1.0_dp),&
      & TUnit("turns               ", 2.0_dp * pi),&
      & TUnit("gradians            ", pi / 200.0_dp)&
      & ]


  type(TUnit), parameter :: massDensityUnits(*) = [&
      & TUnit("kg/l                ", 1.0e+3_dp * kg__au / (1.0e10_dp * AA__Bohr)**3),&
      & TUnit("g/m^3               ", 1.0e+3_dp * kg__au / (1.0e10_dp * AA__Bohr)**3),&
      & TUnit("g/meter^3           ", 1.0e+3_dp * kg__au / (1.0e10_dp * AA__Bohr)**3),&
      & TUnit("g/l                 ", kg__au / (1.0e10_dp * AA__Bohr)**3),&
      & TUnit("kg/m^3              ", kg__au / (1.0e10_dp * AA__Bohr)**3),&
      & TUnit("kg/meter^3          ", kg__au / (1.0e10_dp * AA__Bohr)**3),&
      & TUnit("amu/aa^3            ", amu__au / AA__Bohr**3),&
      & TUnit("amu/angstrom^3      ", amu__au / AA__Bohr**3),&
      & TUnit("amu/a^3             ", amu__au / AA__Bohr**3),&
      & TUnit("amu/pm^3            ", amu__au/(1.0e-2_dp * AA__Bohr)**3),&
      & TUnit("amu/picometer^3     ", amu__au/(1.0e-2_dp * AA__Bohr)**3),&
      & TUnit("da/aa^3             ", amu__au / AA__Bohr**3),&
      & TUnit("da/angstrom^3       ", amu__au / AA__Bohr**3),&
      & TUnit("da/a^3              ", amu__au / AA__Bohr**3),&
      & TUnit("da/pm^3             ", amu__au/(1.0e-2_dp * AA__Bohr)**3),&
      & TUnit("da/picometer^3      ", amu__au/(1.0e-2_dp * AA__Bohr)**3),&
      & TUnit("dalton/aa^3         ", amu__au / AA__Bohr**3),&
      & TUnit("dalton/angstrom^3   ", amu__au / AA__Bohr**3),&
      & TUnit("dalton/a^3          ", amu__au / AA__Bohr**3),&
      & TUnit("dalton/pm^3         ", amu__au/(1.0e-2_dp * AA__Bohr)**3),&
      & TUnit("dalton/picometer^3  ", amu__au/(1.0e-2_dp * AA__Bohr)**3),&
      & TUnit("au                  ", 1.0_dp)]


contains


  !> User defined formatted output routine for the unit conversions
  subroutine TUnit_writeFormatted(this, unit, iotype, vlist, iostat, iomsg)

    !> Instance
    class(TUnit), intent(in) :: this

    !> Unit to write to
    integer, intent(in) :: unit

    !> IO type
    character(*), intent(in) :: iotype

    !> Output parameters (will be ignored)
    integer, intent(in) :: vlist(:)

    !> I/O status
    integer, intent(out) :: iostat

    !> Eventual error message
    character(*), intent(inout) :: iomsg

    if (this%invertValue) then
      write(unit,"(a20, a, e24.14, a)") this%name, ':', this%conversionFact, " / x"
    else
      write(unit,"(a20, a, e24.14, a)") this%name, ':', this%conversionFact, " * x"
    end if
    iostat = 0

  end subroutine TUnit_writeFormatted


#:for RANK in UNIT_CONVERTER_RANKS

  !> Applies unit conversion to a given value
  subroutine convertUnitR${RANK}$(units, unitName, val, status)

    !> Array of possible units
    type(TUnit), intent(in) :: units(:)

    !> Name of the unit from which the value should be converted
    character(*), intent(in) :: unitName

    !> Original value on entry, converted value on return
    real(dp), intent(inout) :: val${FORTRAN_ARG_DIM_SUFFIX(RANK)}$

    !> Status (statusCodes%OK or statusCodes%unitNotFound)
    integer, intent(out) :: status

    integer :: ind

    status = statusCodes%OK
    ind = getUnitIndex_(units, unitName)
    if (ind == 0) then
      status = statusCodes%unitNotFound
      return
    end if

    if (units(ind)%invertValue) then
      val = units(ind)%conversionFact / val
    else
      val = units(ind)%conversionFact * val
    end if

  end subroutine convertUnitR${RANK}$

#:endfor


  !> Returns the index of the unit with the given name or zero if not found.
  function getUnitIndex_(units, unitName) result(ind)
    type(TUnit), intent(in) :: units(:)
    character(*), intent(in) :: unitName
    integer :: ind

    character(maxUnitNameLen) :: unitNameLower

    unitNameLower = tolower(unitName)
    do ind = 1, size(units)
      if (units(ind)%name == unitNameLower) return
    end do
    ind = 0

  end function getUnitIndex_


end module dftbp_common_unitconversion
