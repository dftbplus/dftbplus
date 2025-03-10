!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program *
!  under the LGPL licence.                                                *
!**************************************************************************
module dftbp_poisson_gradient


real(8), allocatable, save :: gr(:,:)     !3,NNDIM
real(8), allocatable, save :: hgrad(:,:)  !3,3*NNDIM
integer, allocatable, save :: conat(:)    !NNDIM+1
real(8), allocatable, save :: convec(:,:) !3,NNDIM
logical, save :: constr


end module dftbp_poisson_gradient
