!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for square dense matrix specification
module dftbp_type_densedescr
#:if WITH_SCALAPACK
  use dftbp_extlibs_scalapackfx, only : DLEN_
#:endif
  implicit none

  private
  public :: TDenseDescr


  type :: TDenseDescr

  #:if WITH_SCALAPACK
    !> BLACS specifier for the orbital sized matrix
    integer :: blacsOrbSqr(DLEN_)

    !> BLACS specifier for the reordered matrix
    integer :: blacsColumnSqr(DLEN_)

  #:endif

    !> Dense matrix indexing by the start of orbitals for each atom.
    !>
    !> Note: for Pauli matrix it contains the indexing of the left upper block only
    !>
    integer, allocatable :: iAtomStart(:)

    !> Dimension of the matrix
    integer :: fullSize

    !> Nr. of atomic orbitals represented in the matrix.
    !>
    !> Equals to fullSize for normal matrices and fullSize / 2 for Pauli matrices
    !>
    integer :: nOrb

    !> Whether atomic matrix represents a two-component Pauli matrix.
    logical :: t2Component

  end type TDenseDescr

end module dftbp_type_densedescr
