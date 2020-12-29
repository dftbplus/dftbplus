!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains the Van der Waal parameters for the UFF force field.
module dftbp_dispuffdata
  use dftbp_accuracy
  use dftbp_message, only : error
  use dftbp_charmanip, only : tolower
  use dftbp_constants
  implicit none
  private

  public :: getUffValues


  !> Contains UFF data (chemical symbol, vdW distance, depth of the LJ potential)
  type TUFF

    !> Chemical symbol
    character(2) :: symbol

    !> vdW minimum distance
    real(dp) :: distance

    !> depth of minimum
    real(dp) :: energy
  end type TUFF


  !> Values are in AA and kcal/mol
  !>
  type(TUFF), parameter :: database(*) = [&
      & TUFF("h ", 2.886_dp, 0.044_dp), TUFF("he", 2.362_dp, 0.056_dp),&
      & TUFF("li", 2.451_dp, 0.025_dp), TUFF("be", 2.745_dp, 0.085_dp),&
      & TUFF("b ", 4.083_dp, 0.180_dp), TUFF("c ", 3.851_dp, 0.105_dp),&
      & TUFF("n ", 3.660_dp, 0.069_dp), TUFF("o ", 3.500_dp, 0.060_dp),&
      & TUFF("f ", 3.364_dp, 0.050_dp), TUFF("ne", 3.243_dp, 0.042_dp),&
      & TUFF("na", 2.983_dp, 0.030_dp), TUFF("mg", 3.021_dp, 0.111_dp),&
      & TUFF("al", 4.499_dp, 0.505_dp), TUFF("si", 4.295_dp, 0.402_dp),&
      & TUFF("p ", 4.147_dp, 0.305_dp), TUFF("s ", 4.035_dp, 0.274_dp),&
      & TUFF("cl", 3.947_dp, 0.227_dp), TUFF("ar", 3.868_dp, 0.185_dp),&
      & TUFF("k ", 3.812_dp, 0.035_dp), TUFF("ca", 3.999_dp, 0.238_dp),&
      & TUFF("sc", 3.295_dp, 0.019_dp), TUFF("ti", 3.175_dp, 0.017_dp),&
      & TUFF("v ", 3.144_dp, 0.016_dp), TUFF("cr", 3.023_dp, 0.015_dp),&
      & TUFF("mn", 2.961_dp, 0.013_dp), TUFF("fe", 2.912_dp, 0.013_dp),&
      & TUFF("co", 2.872_dp, 0.014_dp), TUFF("ni", 2.834_dp, 0.015_dp),&
      & TUFF("cu", 3.495_dp, 0.005_dp), TUFF("zn", 2.763_dp, 0.124_dp),&
      & TUFF("ga", 4.383_dp, 0.415_dp), TUFF("ge", 4.280_dp, 0.379_dp),&
      & TUFF("as", 4.230_dp, 0.309_dp), TUFF("se", 4.205_dp, 0.291_dp),&
      & TUFF("br", 4.189_dp, 0.251_dp), TUFF("kr", 4.141_dp, 0.220_dp),&
      & TUFF("rb", 4.114_dp, 0.040_dp), TUFF("sr", 3.641_dp, 0.235_dp),&
      & TUFF("y ", 3.345_dp, 0.072_dp), TUFF("zr", 3.124_dp, 0.069_dp),&
      & TUFF("nb", 3.165_dp, 0.059_dp), TUFF("mo", 3.052_dp, 0.056_dp),&
      & TUFF("tc", 2.998_dp, 0.048_dp), TUFF("ru", 2.963_dp, 0.056_dp),&
      & TUFF("rh", 2.929_dp, 0.053_dp), TUFF("pd", 2.899_dp, 0.048_dp),&
      & TUFF("ag", 3.148_dp, 0.036_dp), TUFF("cd", 2.848_dp, 0.228_dp),&
      & TUFF("in", 4.463_dp, 0.599_dp), TUFF("sn", 4.392_dp, 0.567_dp),&
      & TUFF("sb", 4.420_dp, 0.449_dp), TUFF("te", 4.470_dp, 0.398_dp),&
      & TUFF("i ", 4.500_dp, 0.339_dp), TUFF("xe", 4.404_dp, 0.332_dp),&
      & TUFF("cs", 4.517_dp, 0.045_dp), TUFF("ba", 3.703_dp, 0.364_dp),&
      & TUFF("la", 3.522_dp, 0.017_dp), TUFF("ce", 3.556_dp, 0.013_dp),&
      & TUFF("pr", 3.606_dp, 0.010_dp), TUFF("nd", 3.575_dp, 0.010_dp),&
      & TUFF("pm", 3.547_dp, 0.009_dp), TUFF("sm", 3.520_dp, 0.008_dp),&
      & TUFF("eu", 3.493_dp, 0.008_dp), TUFF("gd", 3.368_dp, 0.009_dp),&
      & TUFF("tb", 3.451_dp, 0.007_dp), TUFF("dy", 3.428_dp, 0.007_dp),&
      & TUFF("ho", 3.409_dp, 0.007_dp), TUFF("er", 3.391_dp, 0.007_dp),&
      & TUFF("tm", 3.374_dp, 0.006_dp), TUFF("yb", 3.355_dp, 0.228_dp),&
      & TUFF("lu", 3.640_dp, 0.041_dp), TUFF("hf", 3.141_dp, 0.072_dp),&
      & TUFF("ta", 3.170_dp, 0.081_dp), TUFF("w ", 3.069_dp, 0.067_dp),&
      & TUFF("re", 2.954_dp, 0.066_dp), TUFF("os", 3.120_dp, 0.037_dp),&
      & TUFF("ir", 2.840_dp, 0.073_dp), TUFF("pt", 2.754_dp, 0.080_dp),&
      & TUFF("au", 3.293_dp, 0.039_dp), TUFF("hg", 2.705_dp, 0.385_dp),&
      & TUFF("tl", 4.347_dp, 0.680_dp), TUFF("pb", 4.297_dp, 0.663_dp),&
      & TUFF("bi", 4.370_dp, 0.518_dp), TUFF("po", 4.709_dp, 0.325_dp),&
      & TUFF("at", 4.750_dp, 0.284_dp), TUFF("rn", 4.765_dp, 0.248_dp),&
      & TUFF("fr", 4.900_dp, 0.050_dp), TUFF("ra", 3.677_dp, 0.404_dp),&
      & TUFF("ac", 3.478_dp, 0.033_dp), TUFF("th", 3.396_dp, 0.026_dp),&
      & TUFF("pa", 3.424_dp, 0.022_dp), TUFF("u ", 3.395_dp, 0.022_dp),&
      & TUFF("np", 3.424_dp, 0.019_dp), TUFF("pu", 3.424_dp, 0.016_dp),&
      & TUFF("am", 3.381_dp, 0.014_dp), TUFF("cm", 3.326_dp, 0.013_dp),&
      & TUFF("bk", 3.339_dp, 0.013_dp), TUFF("cf", 3.313_dp, 0.013_dp),&
      & TUFF("es", 3.299_dp, 0.012_dp), TUFF("fm", 3.286_dp, 0.012_dp),&
      & TUFF("md", 3.274_dp, 0.011_dp), TUFF("no", 3.248_dp, 0.011_dp),&
      & TUFF("lw", 3.236_dp, 0.011_dp)]

contains


  !> Returns distance and energy parameters of the UFF field.
  subroutine getUffValues(name, distance, energy, found)

    !> Name of the element for look for.
    character(len=*), intent(in) :: name

    !> Distance (in Bohr) at return.
    real(dp), intent(out) :: distance

    !> Energy (in Hartree) at return.
    real(dp), intent(out) :: energy

    !> Flags if the element has been found. If this parameter is not specified and the element is
    !> not found, the program stops.
    logical, intent(out), optional :: found

    character(2) :: symbol
    integer :: ii

    symbol = trim(tolower(name))
    do ii = 1, size(database)
      if (database(ii)%symbol == symbol) then
        distance = database(ii)%distance * AA__Bohr
        energy = database(ii)%energy * kcal_mol__Hartree
        if (present(found)) then
          found = .true.
        end if
        return
      end if
    end do
    if (present(found)) then
      found = .false.
    else
      call error("UFF database search for element '" // trim(name) // "' failed")
    end if

  end subroutine getUffValues

end module dftbp_dispuffdata
