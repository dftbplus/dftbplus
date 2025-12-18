!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains types with electronic/hamiltonian interaction cutoff information.
module dftbp_type_eleccutoffs
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TCutoffs


  !> Interaction cutoff distances.
  type :: TCutoffs

    !> Cutoff for overlap and Hamiltonian according to SK-files
    real(dp) :: skCutOff

    !> Cutoff for overlap and Hamiltonian according to SK-files minus possible cutoff reduction
    real(dp) :: camCutOff

    !> Cutoff for real-space g-summation in CAM Hartree-Fock contributions
    real(dp), allocatable :: gSummationCutoff

    !> Number of unit cells along each supercell folding direction to substract from MIC
    !! Wigner-Seitz cell construction
    integer, allocatable :: wignerSeitzReduction

    !> Cutoff for truncated long-range Gamma integral
    real(dp), allocatable :: gammaCutoff

    !> Max. of all cutoffs
    real(dp) :: mCutOff

  end type TCutoffs

end module dftbp_type_eleccutoffs
