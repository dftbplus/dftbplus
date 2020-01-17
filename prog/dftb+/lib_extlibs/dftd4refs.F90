module dftbp_dftd4refs
  use dftbp_accuracy, only : dp
  implicit none
  private

  public :: sscale, clsq, secaiw, alphaiw, clsh, refn, refsys, refcn, refcovcn
  public :: hcount, ascale

  include 'dftd4_references.fh'

end module dftbp_dftd4refs
