#:include 'common.fypp'

module densedescr
  use scalapackfx
  implicit none
  private

  public :: TDenseDescr

  type :: TDenseDescr
  #:if WITH_SCALAPACK
    integer :: blacsOrbSqr(DLEN_)
  #:endif
    integer, allocatable :: iDenseStart(:)
    integer :: fullSize
  end type TDenseDescr

end module densedescr
