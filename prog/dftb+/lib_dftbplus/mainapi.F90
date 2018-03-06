!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module mainapi
  use environment, only : TEnvironment
  use accuracy, only : dp
  use main, only : processGeometry
  use initprogram, only : initProgramVariables, destructProgramVariables, coord0, latVec,&
      & tCoordsChanged, tLatticeChanged, energy, derivs
  implicit none
  private

  public :: initProgramVariables, destructProgramVariables
  public :: updateGeometry
  public :: getEnergy, getGradients

contains


  subroutine updateGeometry(coords, latVecs)
    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(in), optional :: latVecs(:,:)
    
    coord0(:,:) = coords
    tCoordsChanged = .true.
    if (present(latVecs)) then
      latVec(:,:) = latVecs
      tLatticeChanged = .true.
    else
      tLatticeChanged = .false.
    end if

  end subroutine updateGeometry


  subroutine getEnergy(env, merminEnergy)
    type(TEnvironment), intent(inout) :: env
    real(dp), intent(out) :: merminEnergy

    call recalcGeometry(env)
    merminEnergy = energy%EMermin
    
  end subroutine getEnergy


  subroutine getGradients(env, gradients)
    type(TEnvironment), intent(inout) :: env
    real(dp), intent(out) :: gradients(:,:)

    call recalcGeometry(env)
    gradients(:,:) = derivs
    
  end subroutine getGradients


  subroutine recalcGeometry(env)
    type(TEnvironment), intent(inout) :: env

    logical :: tStopDriver, tStopScc
    
    if (tLatticeChanged .or. tCoordsChanged) then
      call processGeometry(env, 1, 1, .false., tStopDriver, tStopScc)
      tLatticeChanged = .false.
      tCoordsChanged = .false.
    end if

  end subroutine recalcGeometry

  
end module mainapi
