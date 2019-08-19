!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *  
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program * 
!  under the LGPL licence.                                                *
!**************************************************************************
module gconstants
  
  use gprecision
  
  real(kind=dp), parameter    :: eovh = (1.05420882d-3)   ! A/H
  real(kind=dp), parameter    :: pi =  3.14159265358979323844_dp ! Greek p real
  real(kind=dp), parameter    :: HAR = 27.2113845_dp         ! H/eV
  real(kind=dp), parameter    :: ATU = 0.529177249_dp        ! a.u./Ang
  real(kind=dp), PARAMETER    :: Kb = (3.166830814d-6)    ! H/K
  
  COMPLEX(kind=dp), PARAMETER ::    j = (0.d0,1.d0)  ! CMPX unity
  
  
end module gconstants

