!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains types for geometry control
module dftbp_geometry_control
  implicit none

  private
  public :: TGeomChanges


  !> Modifications allowed for the geometry
  type :: TGeomChanges

    !> Geometry optimisation needed?
    logical :: isGeoOpt =.false.

    !> Optimise coordinates inside unit cell (periodic)?
    logical :: tCoordOpt =.false.

    !> Optimise lattice constants?
    logical :: tLatOpt =.false.

    !> Fix angles between lattice vectors when optimising?
    logical :: tLatOptFixAng =.false.

    !> Fix length of specified lattice vectors when optimising?
    logical :: tLatOptFixLen(3) =.false.

    !> Optimise lattice isotropically
    logical :: tLatOptIsotropic =.false.

    !> Is this a MD calculation?
    logical :: tMD =.false.

    !> Barostat used if MD and periodic
    logical :: tBarostat = .false.

    !> Is this a finite difference derivatives calc?
    logical :: tDerivs = .false.

  end type TGeomChanges


end module dftbp_geometry_control
