!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for square dense matrix specification
module densedescr
  use scalapackfx
  implicit none
  private

  public :: TDenseDescr


  type :: TDenseDescr

  #:if WITH_SCALAPACK
    !> BLACS specifier for the matrix
    integer :: blacsOrbSqr(DLEN_)
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

end module densedescr
