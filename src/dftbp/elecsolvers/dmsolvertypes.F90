!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Types for electronic solution of hamiltonian
module dftbp_elecsolvers_dmsolvertypes
  implicit none

  private
  public :: densityMatrixTypes

  !> Namespace for possible density matrix generation methods
  type :: TDensityMatrixTypesEnum

    !> No density matrix generated
    integer :: none = 0

    !> Electronic solver directly returns density matrix
    integer :: elecSolverProvided = 1

    !> Dense density matrix constructed from eigen values and vectors
    integer :: fromEigenVecs = 2

    !> GPU accellerated construction from eigen values and vectors
    integer :: magma_fromEigenVecs = 3

  end type TDensityMatrixTypesEnum


  !> Actual values for densityMatrixTypes.
  type(TDensityMatrixTypesEnum), parameter :: densityMatrixTypes = TDensityMatrixTypesEnum()

end module dftbp_elecsolvers_dmsolvertypes
