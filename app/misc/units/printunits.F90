!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#! suffix and kinds for real types
#:set CONVERSIONS = [('Length', 'lengthUnits'), ('Inverse length', 'inverseLengthUnits'),&
  & ("Energy", "energyUnits"), ("Force", "forceUnits"), ("Time", "timeUnits"),&
  & ("Frequency", "freqUnits"), ("Volume", "volumeUnits"), ("Charge", "chargeUnits"),&
  & ("Elec. Field", "EFieldUnits"), ("Mag. Field", "BFieldUnits"), ("Pressure", "pressureUnits"),&
  & ("Velocity", "velocityUnits"), ("Elec. dipole", "dipoleUnits"), ("Mass", "massUnits"),&
  & ("Angular", "angularUnits"), ("Mass Density", "massDensityUnits")]

!> Printing out the conversion factors for the different units
program printunits
  use dftbp_common_unitconversion
  use dftbp_common_constants
  implicit none

  integer :: ii

  write(*,*)"Convert from unit to a.u. by multiplying with"
  #:for names, units in CONVERSIONS
    write(*,*)
    write(*,"(A)")"${names}$:"
    do ii = 1, size(${units}$)
      write(*,"(1X,A20,E24.15)") ${units}$(ii)%name, ${units}$(ii)%convertValue
    end do
  #:endfor

end program printunits
