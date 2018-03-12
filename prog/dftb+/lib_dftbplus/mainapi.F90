!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module mainapi
  use environment, only : TEnvironment
  use assert
  use accuracy, only : dp
  use main, only : processGeometry
  use initprogram, only : initProgramVariables, destructProgramVariables, coord0, latVec,&
      & tCoordsChanged, tLatticeChanged, energy, derivs, TRefExtPot, refExtPot, tExtField, orb,&
      & nAtom, nSpin, q0, qOutput
  implicit none
  private

  public :: initProgramVariables, destructProgramVariables
  public :: setGeometry, setExternalPotential
  public :: getEnergy, getGradients, getGrossCharges

contains


  subroutine setGeometry(coords, latVecs)
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

  end subroutine setGeometry


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


  subroutine getGrossCharges(atomCharges)
    real(dp), intent(out) :: atomCharges(:)

    atomCharges(:) = sum(q0(:, :, 1) - qOutput(:, :, 1), dim=1)
    
  end subroutine getGrossCharges


  !> Sets up an external electrostatic potential.
  !>
  !> Sign convention: charge of electron is considered to be negative.
  !>
  subroutine setExternalPotential(atomPot, shellPot, potGrad)

    !> Atomic external potential
    real(dp), intent(in), optional :: atomPot(:)

    !> Shell resolved electrostatic potential
    real(dp), intent(in), optional :: shellPot(:,:)

    !> Gradient of the electrostatic potential
    real(dp), intent(in), optional :: potGrad(:,:)

    ! Using explicit allocation instead of F2003 automatic ones in order to stop eventual
    ! shape mismatches already at this point rather than later deep in the main code
    if (present(atomPot)) then
      if (.not. allocated(refExtPot%atomPot)) then
        allocate(refExtPot%atomPot(nAtom, nSpin))
      end if
      @:ASSERT(all(shape(atomPot) == [nAtom]))
      refExtPot%atomPot(:,1) = -atomPot
    end if
    if (present(shellPot)) then
      if (.not. allocated(refExtPot%shellPot)) then
        allocate(refExtPot%shellPot(orb%mShell, nAtom, nSpin))
      end if
      @:ASSERT(all(shape(shellPot) == [orb%mShell, nAtom]))
      refExtPot%shellPot(:,:,1) = -shellPot
    end if
    if (present(potGrad)) then
      if (.not. allocated(refExtPot%potGrad)) then
        allocate(refExtPot%potGrad(3, nAtom))
      end if
      @:ASSERT(all(shape(potGrad) == [3, nAtom]))
      refExtPot%potGrad(:,:) = -potGrad
    end if
    tExtField = .true.

  end subroutine setExternalPotential



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
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
