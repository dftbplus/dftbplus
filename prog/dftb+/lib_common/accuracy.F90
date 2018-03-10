!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains a list of constants for the control of precision of the
!! calculation, both for the fortran numerical model and defaults for the various algorithms in the
!! code.
!! Not all routines use the string length specifications to set their character string lengths.
module accuracy
  implicit none

  !> precision of the real data type
  integer, parameter :: dp = kind(1.0d0)

  !> precision of the complex data type
  integer, parameter :: cp = dp

  !> length of a short string
  integer, parameter :: sc = 10

  !> length of a medium length string
  integer, parameter :: mc = 50

  !> length of a long string
  integer, parameter :: lc = 200


  !> Real double precision - do not edit
  integer, parameter :: rdp = kind(0.0d0)

  !> Real single precision - do not edit
  integer, parameter :: rsp = kind(0.0)

  ! Program technical constants


  !> Length of the tail after the SK table at which point elements = 0
  real(dp), parameter :: distFudge = 1.0_dp


  !> Length of the tail after the SK table for the old extrapolation alg.
  real(dp), parameter :: distFudgeOld = 0.3_dp


  !> Desired tolerance for number of total electrons when finding electron chemical potential The
  !> Fermi level is searched to give the number of electrons as accurate as elecTol. If bisection
  !> ends and difference between nr. of electrons calculated/theoretical is bigger than elecTolMax,
  !> the program stops.
  real(dp), parameter :: elecTol = 1.0e-15_dp


  !> Maximal allowed tolerance for number of total electrons when finding Ef
  !> or when reading in charges from external file
  real(dp), parameter :: elecTolMax = 1.0e-7_dp


  !> Minimal temperature, temperatures below that are replaced by this value
  real(dp), parameter :: minTemp = 1.0e-8_dp


  !> Tolerance for atomic distances. Atoms closer than that are regarded to sit on the same
  !> positions. (Dummy atoms)
  real(dp), parameter :: tolSameDist = 1.0e-5_dp


  !> Tolerance for atomic square distances
  real(dp), parameter :: tolSameDist2 = tolSameDist**2


  !> Minimal distance between neihbors. (Neighbor distances smaller than that
  !> are meaningless because the parametrisation usually do not cover this

  !> region.)
  real(dp), parameter :: minNeighDist = 1.0e-2_dp


  !> Minimal square distance between neighbors
  real(dp), parameter :: minNeighDist2 = minNeighDist**2


  !> Cut-off value to calculate the short-range part of gamma_ab
  real(dp), parameter :: minShortGamma = 1.0e-10_dp


  !> Cut-off value to calculate the short-range part of Ewald sum
  real(dp), parameter :: minShortEwald = 1.0e-10_dp


  !> Tolerance for error in cut-off of short-range part of gamma_ab
  real(dp), parameter :: tolShortGamma = 1.0e-10_dp


  !> Tolerance for error in cut-off of short-range part of gamma_ab
  real(dp), parameter :: tolShortEwald = 1.0e-10_dp


  !> Minimum value for alpha in Ewald sum
  real(dp), parameter :: tolMinAlpha = 1.0e-4_dp


  !> Tolerance for minimum possible value of an atomic Hubbard U
  real(dp), parameter :: MinHubTol = 1.0e-6_dp


  !> Tolerance for minimum possible difference in values of Hubbard U
  real(dp), parameter :: MinHubDiff = 0.3125_dp*1.0e-5_dp


  !> Nr. of max. bisection steps
  integer, parameter :: nSearchIter = 30


  !> Exponential function treated as infinity with arguments higher than this
  !> (=-int(log(epsilon(1.0_8)))).
  real(dp), parameter :: mExpArg = 36.043653389117154_dp


  !> Tolerance for the error in the dispersion
  real(dp), parameter :: tolDispersion = 1.0e-9_dp


  !> Tolerance for the dispersion damping function being considered 1
  real(dp), parameter :: tolDispDamp = 1.0e-10_dp

end module accuracy
