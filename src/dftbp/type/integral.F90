!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Data types to handle overlap related integrals
module dftbp_type_integral
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TIntegral

  !> Container to store overlap related integrals
  type :: TIntegral

    !> Overlap integrals in atomic block sparse form
    real(dp), allocatable :: overlap(:)

    !> Real Hamiltonian integrals in atomic block sparse form
    real(dp), allocatable :: hamiltonian(:, :)

    !> Imaginary Hamiltonian integrals in atomic block sparse form
    real(dp), allocatable :: iHamiltonian(:, :)

  end type TIntegral

contains
end module dftbp_type_integral
