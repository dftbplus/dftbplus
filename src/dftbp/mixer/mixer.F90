!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set FLAVOURS = [('cmplx', 'complex', 'Cmplx'), ('real', 'real', 'Real')]

!> Defines an abstract mixer as a base class, specifying the common interface for the
!> desired actual mixer and providing shared methods to flatten 6D and 3D data to 1D.
module dftbp_mixer_mixer
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TMixerReal, TMixerCmplx

#:for NAME, TYPE, LABEL in FLAVOURS
  type, abstract :: TMixer${LABEL}$
    contains
      procedure(IReset${LABEL}$), deferred :: reset

      !> The mixing implementation
      procedure(IMix1D${LABEL}$), deferred :: mix1D

      !> Flatten 6d and 3d mixing to 1d
      procedure :: mix3D => TMixer${LABEL}$_mix3D
      procedure :: mix6D => TMixer${LABEL}$_mix6D

      generic :: mix => mix1D, mix3D, mix6D
  end type TMixer${LABEL}$

  abstract interface

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
#:endfor

contains

#:for NAME, TYPE, LABEL in FLAVOURS
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

    call this%mix1D(qInpRes, qDiff)

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

    call this%mix1D(qInpRes, qDiff)

  end subroutine TMixer${LABEL}$_mix6D
#:endfor


end module dftbp_mixer_mixer
