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

    !> Dense matrix indexing by the start of orbitals for each atom
    integer, allocatable :: iDenseStart(:)

    !> Dimension of the matrix
    integer :: fullSize

  end type TDenseDescr

end module densedescr
