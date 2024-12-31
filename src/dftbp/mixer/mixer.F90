!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set FLAVOURS = [('cmplx', 'complex', 'Cmplx'), ('real', 'real', 'Real')]

!> Provides a general mixer which contains the desired actual mixer.
module dftbp_mixer_mixer
  use dftbp_common_accuracy, only : dp
  implicit none

  private
#:for NAME, TYPE, LABEL in FLAVOURS
  public :: TMixer${LABEL}$, TMixer${LABEL}$_reset, TMixer${LABEL}$_mix
#:endfor
  public :: mixerTypes, TMixerInput




  type :: TMixerTypesEnum
    integer :: simple = 1
    integer :: anderson = 2
    integer :: broyden = 3
    integer :: diis = 4
  end type TMixerTypesEnum

  !> Contains mixer types
  type(TMixerTypesEnum), parameter :: mixerTypes = TMixerTypesEnum()

  !> Mixer specific Input
  type :: TMixerInput
    integer :: iMixSwitch = -1
    integer :: iGenerations = 0
    real(dp) :: mixParam = 0.0_dp
    logical :: tFromStart = .true.

    real(dp) :: broydenOmega0 = 0.01_dp
    real(dp) :: broydenMinWeight = 1.0_dp
    real(dp) :: broydenMaxWeight = 1.0e5_dp
    real(dp) :: broydenWeightFac = 1.0e-2_dp

    real(dp) :: andersonInitMixing = 0.01_dp
    integer :: andersonNrDynMix = 0
    real(dp), allocatable :: andersonDynMixParams(:,:)
    real(dp) :: andersonOmega0 = 1.0e-2_dp

    ! Beside the above mixer specific input, the
    ! Broyden Mixer needs to know the maximum iteration count.
    integer :: maxSccIter = -1
  end type TMixerInput

#:for NAME, TYPE, LABEL in FLAVOURS
  type, abstract :: TMixer${LABEL}$
    contains
        procedure(IInitialise${LABEL}$), deferred :: init
        
        procedure(IReset${LABEL}$), deferred :: reset

        procedure(IMix1D${LABEL}$), deferred :: mix1D
  end type TMixer${LABEL}$

  abstract interface 

    subroutine IInitialise${LABEL}$(this, mixerInp)
      import :: TMixer${LABEL}$, TMixerInput
      class(TMixer${LABEL}$), intent(out) :: this
      type(TMixerInput), intent(in) :: mixerInp
    end subroutine IInitialise${LABEL}$

    subroutine IReset${LABEL}$(this, nElem)
      import :: TMixer${LABEL}$
      class(TMixer${LABEL}$), intent(inout) :: this
      integer, intent(in) :: nElem
    end subroutine IReset${LABEL}$


    subroutine IMix1D${LABEL}$(this, qInpResult, qDiff)
      import :: TMixer${LABEL}$, dp
      class(TMixer${LABEL}$), intent(inout) :: this
      ${TYPE}$(dp), intent(inout) :: qInpResult(:)
      ${TYPE}$(dp), intent(in) :: qDiff(:)
    end subroutine IMix1D${LABEL}$

  end interface

  !> Does the actual mixing
  interface TMixer${LABEL}$_mix
    module procedure TMixer${LABEL}$_mix1D
    module procedure TMixer${LABEL}$_mix3D
    module procedure TMixer${LABEL}$_mix6D
  end interface TMixer${LABEL}$_mix
#:endfor



contains

#:for NAME, TYPE, LABEL in FLAVOURS
  !> Resets the mixer.
  subroutine TMixer${LABEL}$_reset(this, nElem)

    !> Mixer instance
    class(TMixer${LABEL}$), intent(inout) :: this

    !> Size of the vectors to mix
    integer, intent(in) :: nElem

    call this%reset(nElem)

  end subroutine TMixer${LABEL}$_reset

  !> Mixes two vectors.
  subroutine TMixer${LABEL}$_mix1D(this, qInpRes, qDiff)

    !> Mixer instance
    class(TMixer${LABEL}$), intent(inout) :: this

    !> Input vector on entry, result vector on exit
    ${TYPE}$(dp), intent(inout) :: qInpRes(:)

    !> Difference between input and output vectors (measure of lack of convergence)
    ${TYPE}$(dp), intent(in) :: qDiff(:)

    call this%mix1D(qInpRes, qDiff)

  end subroutine TMixer${LABEL}$_mix1D


  !> Mixes two 3D matrices.
  subroutine TMixer${LABEL}$_mix3D(this, qInpResSqr, qDiffSqr)

    !> Mixer instance
    class(TMixer${LABEL}$), intent(inout) :: this

    !> Input vector on entry, result vector on exit
    ${TYPE}$(dp), intent(inout), contiguous, target :: qInpResSqr(:,:,:)

    !> Difference between input and output vectors (measure of lack of convergence)
    ${TYPE}$(dp), intent(in), contiguous, target :: qDiffSqr(:,:,:)

    !! Difference between input and output vectors (1D pointer)
    ${TYPE}$(dp), pointer :: qDiff(:)

    !! Input vector on entry, result vector on exit (1D pointer)
    ${TYPE}$(dp), pointer :: qInpRes(:)

    qDiff(1:size(qDiffSqr)) => qDiffSqr
    qInpRes(1:size(qInpResSqr)) => qInpResSqr

    call TMixer${LABEL}$_mix1D(this, qInpRes, qDiff)

  end subroutine TMixer${LABEL}$_mix3D


  !> Mixes two 6D matrices.
  subroutine TMixer${LABEL}$_mix6D(this, qInpResSqr, qDiffSqr)

    !> Mixer instance
    class(TMixer${LABEL}$), intent(inout) :: this

    !> Input vector on entry, result vector on exit
    ${TYPE}$(dp), intent(inout), contiguous, target :: qInpResSqr(:,:,:,:,:,:)

    !> Difference between input and output vectors (measure of lack of convergence)
    ${TYPE}$(dp), intent(in), contiguous, target :: qDiffSqr(:,:,:,:,:,:)

    !!Difference between input and output vectors (1D pointer)
    ${TYPE}$(dp), pointer :: qDiff(:)

    !! Input vector on entry, result vector on exit (1D pointer)
    ${TYPE}$(dp), pointer :: qInpRes(:)

    qDiff(1:size(qDiffSqr)) => qDiffSqr
    qInpRes(1:size(qInpResSqr)) => qInpResSqr

    call TMixer${LABEL}$_mix1D(this, qInpRes, qDiff)

  end subroutine TMixer${LABEL}$_mix6D


#:endfor

end module dftbp_mixer_mixer
