!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *  
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program * 
!  under the LGPL licence.                                                *
!**************************************************************************
module gradient


real(8), allocatable, save :: gr(:,:)     !3,NNDIM
real(8), allocatable, save :: hgrad(:,:)  !3,3*NNDIM
integer, allocatable, save :: conat(:)    !NNDIM+1
real(8), allocatable, save :: convec(:,:) !3,NNDIM
logical, save :: constr


end module gradient
