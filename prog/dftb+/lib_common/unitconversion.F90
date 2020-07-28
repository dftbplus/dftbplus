!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains names of various units and their conversion factors to the corresponding unit used
!> internal in the code (atomic units).
module dftbp_unitconversion
  use dftbp_constants
  implicit none

  public


  !> Contains name of a unit and its conversion factor
  type unit

    !> label for conversion factor
    character(20) :: name

    !> its value
    real(dp) :: convertValue
  end type unit


  !> Number of length units
  integer, parameter :: nLengthUnit = 9


  !> Length units
  type(unit), parameter :: lengthUnits(nLengthUnit) = [&
      & unit("angstrom            ", AA__Bohr),&
      & unit("aa                  ", AA__Bohr),&
      & unit("a                   ", AA__Bohr),&
      & unit("meter               ", 1.0e10_dp * AA__Bohr),&
      & unit("m                   ", 1.0e10_dp * AA__Bohr),&
      & unit("bohr                ", 1.0_dp),&
      & unit("pm                  ", 1.0e-2_dp * AA__Bohr),&
      & unit("picometer           ", 1.0e-2_dp * AA__Bohr),&
      & unit("au                  ", 1.0_dp)&
      & ]


  !> Number of length units
  integer, parameter :: nInverseLengthUnit = 9


  !> Length units
  type(unit), parameter :: inverseLengthUnits(nInverseLengthUnit) = [&
      & unit("1/angstrom          ", Bohr__AA),&
      & unit("1/aa                ", Bohr__AA),&
      & unit("1/a                 ", Bohr__AA),&
      & unit("1/meter             ", 1.0e-10_dp * Bohr__AA),&
      & unit("1/m                 ", 1.0e-10_dp * Bohr__AA),&
      & unit("1/bohr              ", 1.0_dp),&
      & unit("1/pm                ", 1.0e+2_dp * Bohr__AA),&
      & unit("1/picometer         ", 1.0e+2_dp * Bohr__AA),&
      & unit("1/au                ", 1.0_dp)&
      & ]


  !> Number of energy units
  integer, parameter :: nEnergyUnit = 13


  !> Energy units
  type(unit), parameter :: energyUnits(nEnergyUnit) = [&
      & unit("rydberg             ", 0.5_dp),&
      & unit("ry                  ", 0.5_dp),&
      & unit("electronvolt        ", eV__Hartree),&
      & unit("ev                  ", eV__Hartree),&
      & unit("kcal/mol            ", kcal_mol__Hartree),&
      & unit("kelvin              ", Boltzmann),&
      & unit("k                   ", Boltzmann),&
      & unit("cm^-1               ", 1.0_dp/Hartree__cm),&
      & unit("hartree             ", 1.0_dp),&
      & unit("ha                  ", 1.0_dp),&
      & unit("joule               ", J__Hartree),&
      & unit("j                   ", J__Hartree),&
      & unit("au                  ", 1.0_dp)&
      & ]


  !> Number of force units
  integer, parameter :: nForceUnit = 8


  !> Force units
  type(unit), parameter :: forceUnits(nForceUnit) = [&
      & unit("ev/angstrom         ", eV__Hartree / AA__Bohr),&
      & unit("ev/aa               ", eV__Hartree / AA__Bohr),&
      & unit("ev/a                ", eV__Hartree / AA__Bohr),&
      & unit("ha/bohr             ", 1.0_dp),&
      & unit("hartree/bohr        ", 1.0_dp),&
      & unit("j/m                 ", J__Hartree / (1.0e10_dp * AA__Bohr)),&
      & unit("joule/m             ", J__Hartree / (1.0e10_dp * AA__Bohr)),&
      & unit("au                  ", 1.0_dp)&
      & ]


  !> Number of time units
  integer, parameter :: nTimeUnit = 7


  !> Time units
  type(unit), parameter :: timeUnits(nTimeUnit) = [&
      & unit("femtosecond         ", fs__au),&
      & unit("fs                  ", fs__au),&
      & unit("picosecond          ", 1e3_dp * fs__au),&
      & unit("ps                  ", 1e3_dp * fs__au),&
      & unit("second              ", 1e15_dp * fs__au),&
      & unit("s                   ", 1e15_dp * fs__au),&
      & unit("au                  ", 1.0_dp)&
      & ]


  !> Number of frequency units
  integer, parameter :: nFreqUnit = 4


  !> Frequency units
  type(unit), parameter :: freqUnits(nFreqUnit) = [&
      & unit("au                  ", 1.0_dp),&
      & unit("hz                  ", 1e-15_dp * au__fs),&
      & unit("thz                 ", 1e-3_dp * au__fs),&
      & unit("cm^-1               ", 1e-8_dp * Bohr__AA * c)&
      & ]


  !> Number of volume units
  integer, parameter :: nVolumeUnit = 9


  !> Volume units
  type(unit), parameter :: volumeUnits(nVolumeUnit) = [&
      & unit("angstrom^3          ", AA__Bohr**3),&
      & unit("aa^3                ", AA__Bohr**3),&
      & unit("a^3                 ", AA__Bohr**3),&
      & unit("meter^3             ", (1.0e10_dp * AA__Bohr)**3),&
      & unit("m^3                 ", (1.0e10_dp * AA__Bohr)**3),&
      & unit("picometer^3         ", (1.0e-2_dp * AA__Bohr)**3),&
      & unit("pm^3                ", (1.0e-2_dp * AA__Bohr)**3),&
      & unit("bohr^3              ", 1.0_dp),&
      & unit("au                  ", 1.0_dp)&
      & ]


  !> Number of charge units
  integer, parameter :: nChargeUnit = 4


  !> Volume units
  type(unit), parameter :: chargeUnits(nChargeUnit) = [&
      & unit("coulomb             ", Coulomb__au),&
      & unit("c                   ", Coulomb__au),&
      & unit("e                   ", 1.0_dp),&
      & unit("au                  ", 1.0_dp)&
      & ]


  !> Number of electric field units
  integer, parameter :: nEFieldUnit = 5


  !> Electric field units
  type(unit), parameter :: EFieldUnits(nEFieldUnit) = [&
       & unit("v/m                 ", V_m__au),&
       & unit("v/a                 ", 1e10_dp * V_m__au),&
       & unit("v/aa                ", 1e10_dp * V_m__au),&
       & unit("v/angstrom          ", 1e10_dp * V_m__au),&
       & unit("au                  ", 1.0_dp)&
       & ]


  !> Number of magnetic field units
  integer, parameter :: nBFieldUnit = 5


  !> Magnetic field units (Atomic "Gaussian" CGS unit system!)
  type(unit), parameter :: BFieldUnits(nBFieldUnit) = [&
      & unit("t                 ", 1.0E+24_dp*c/(hbar*Coulomb__au*AA__Bohr**2)),&
      & unit("tesla             ", 1.0E+24_dp*c/(hbar*Coulomb__au*AA__Bohr**2)),&
      & unit("gauss             ", 1.0E+20_dp*c/(hbar*Coulomb__au*AA__Bohr**2)),&
      & unit("g                 ", 1.0E+20_dp*c/(hbar*Coulomb__au*AA__Bohr**2)),&
      & unit("au                  ", 1.0_dp)&
      & ]


  !> Number of pressure units
  integer, parameter :: nPressureUnit = 2


  !> Pressure units
  type(unit), parameter :: pressureUnits(nPressureUnit) = [&
      & unit("pa                  ", pascal__au),&
      & unit("au                  ", 1.0_dp)&
      & ]


  !> Number of velocity units
  integer, parameter :: nVelocityUnit = 6


  !> Velocity units
  type(unit), parameter :: velocityUnits(nVelocityUnit) = [&
      & unit("au                  ", 1.0_dp),&
      & unit("m/s                 ", 1e10_dp * AA__Bohr / (1e15_dp * fs__au) ),&
      & unit("pm/fs               ", 1e-2_dp * AA__Bohr / fs__au),&
      & unit("a/ps                ", AA__Bohr / (1e3_dp * fs__au) ),&
      & unit("aa/ps               ", AA__Bohr / (1e3_dp * fs__au) ),&
      & unit("angstrom/ps         ", AA__Bohr / (1e3_dp * fs__au) )&
      & ]


  !> Number of dipole units
  integer, parameter :: nDipoleUnit = 4


  !> Dipole units
  type(unit), parameter :: dipoleUnits(nDipoleUnit) = [&
      & unit("au                  ", 1.0_dp),&
      & unit("debye               ", Debye__au ),&
      & unit("cm                  ", Coulomb__au*1.0e10_dp*AA__Bohr ),&
      & unit("coulombmeter        ", Coulomb__au*1.0e10_dp*AA__Bohr )&
      & ]


  !> Number of mass units
  integer, parameter :: nMassUnit = 6


  !> Mass units
  type(unit), parameter :: MassUnits(nMassUnit) = [&
      & unit("au                  ", 1.0_dp),&
      & unit("amu                 ", amu__au ),&
      & unit("da                  ", amu__au ),&
      & unit("dalton              ", amu__au ),&
      & unit("kg                  ", kg__au ),&
      & unit("g                   ", 1.0e+3_dp*kg__au )&
      & ]


  !> Number of angular units
  integer, parameter :: nAngularUnit = 6


  !> angular units
  type(unit), parameter :: angularUnits(nAngularUnit) = [&
      & unit("degrees             ", pi / 180.0_dp ),&
      & unit("deg                 ", pi / 180.0_dp ),&
      & unit("radian              ", 1.0_dp ),&
      & unit("rad                 ", 1.0_dp ),&
      & unit("turns               ", 2.0_dp * pi ),&
      & unit("gradians            ", pi / 200.0_dp )&
      & ]


  !> Number of mass density units
  integer, parameter :: nMassDensityUnit = 22


  type(unit), parameter :: massDensityUnits(nMassDensityUnit) = [&
      & unit("kg/l                ", 1.0e+3_dp*kg__au/(1.0e10_dp*AA__Bohr)**3),&
      & unit("g/m^3               ", 1.0e+3_dp*kg__au/(1.0e10_dp*AA__Bohr)**3),&
      & unit("g/meter^3           ", 1.0e+3_dp*kg__au/(1.0e10_dp*AA__Bohr)**3),&
      & unit("g/l                 ", kg__au/(1.0e10_dp*AA__Bohr)**3),&
      & unit("kg/m^3              ", kg__au/(1.0e10_dp*AA__Bohr)**3),&
      & unit("kg/meter^3          ", kg__au/(1.0e10_dp*AA__Bohr)**3),&
      & unit("amu/aa^3            ", amu__au/AA__Bohr**3),&
      & unit("amu/angstrom^3      ", amu__au/AA__Bohr**3),&
      & unit("amu/a^3             ", amu__au/AA__Bohr**3),&
      & unit("amu/pm^3            ", amu__au/(1.0e-2_dp*AA__Bohr)**3),&
      & unit("amu/picometer^3     ", amu__au/(1.0e-2_dp*AA__Bohr)**3),&
      & unit("da/aa^3             ", amu__au/AA__Bohr**3),&
      & unit("da/angstrom^3       ", amu__au/AA__Bohr**3),&
      & unit("da/a^3              ", amu__au/AA__Bohr**3),&
      & unit("da/pm^3             ", amu__au/(1.0e-2_dp*AA__Bohr)**3),&
      & unit("da/picometer^3      ", amu__au/(1.0e-2_dp*AA__Bohr)**3),&
      & unit("dalton/aa^3         ", amu__au/AA__Bohr**3),&
      & unit("dalton/angstrom^3   ", amu__au/AA__Bohr**3),&
      & unit("dalton/a^3          ", amu__au/AA__Bohr**3),&
      & unit("dalton/pm^3         ", amu__au/(1.0e-2_dp*AA__Bohr)**3),&
      & unit("dalton/picometer^3  ", amu__au/(1.0e-2_dp*AA__Bohr)**3),&
      & unit("au                  ", 1.0_dp)]


end module dftbp_unitconversion
