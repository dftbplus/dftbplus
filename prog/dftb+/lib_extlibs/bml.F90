!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides interface to the BML library.
module dftbp_bml
  use bml
  implicit none
  private

  public :: bml_matrix_t
  public :: bml_element_real, bml_element_complex
  public :: bml_zero_matrix, bml_clear, bml_deallocate
  public :: bml_get_N, bml_get_M, bml_get_type
  public :: bml_set_row, bml_get_row
  public :: bml_transpose_triangle, bml_adjungate_triangle
  public :: bml_export_to_dense, bml_print_matrix

end module dftbp_bml
