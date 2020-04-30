!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> main module for the DFTB+ API
module dftbp_mainapi
  use dftbp_environment, only : TEnvironment
  use dftbp_accuracy, only : dp
  use dftbp_main, only : processGeometry
  use dftbp_initprogram, only : initProgramVariables, destructProgramVariables, coord0, latVec,&
      & tCoordsChanged, tLatticeChanged, energy, derivs, TRefExtPot, refExtPot, tExtField, orb,&
      & nAtom, nSpin, q0, qOutput, dQAtomEx, sccCalc, isLinResp, tExtChrg, tForces, chrgForces,&
      & qDepExtPot, tStress, totalStress
  use dftbp_assert
  use dftbp_qdepextpotproxy, only : TQDepExtPotProxy
  use dftbp_message, only : error
  implicit none
  private

  public :: initProgramVariables, destructProgramVariables
  public :: setGeometry, setQDepExtPotProxy, setExternalPotential, setExternalCharges
  public :: getEnergy, getGradients, getExtChargeGradients, getGrossCharges, getStressTensor
  public :: nrOfAtoms

contains

  !> Sets up the atomic geometry
  subroutine setGeometry(coords, latVecs)

    !> atom coordinates
    real(dp), intent(in) :: coords(:,:)

    !> lattice vectors, if periodic
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


  !> Returns the free energy of the system at finite temperature
  subroutine getEnergy(env, merminEnergy)

    !> Instance
    type(TEnvironment), intent(inout) :: env

    !> Resulting energy
    real(dp), intent(out) :: merminEnergy

    call recalcGeometry(env)
    merminEnergy = energy%EMermin

  end subroutine getEnergy


  !> get forces on atoms
  subroutine getGradients(env, gradients)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> resulting gradients wrt atom positions
    real(dp), intent(out) :: gradients(:,:)

    if (.not. tForces) then
      call error("Forces not available, you must initialise your calculator&
          & with forces enabled.")
    end if

    call recalcGeometry(env)
    gradients(:,:) = derivs

  end subroutine getGradients


  !> get stress tensor for unit cell
  subroutine getStressTensor(env, stress)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> resulting gradients wrt atom positions
    real(dp), intent(out) :: stress(:,:)

    if (.not. tStress) then
      call error("Stress tensor not available, you must initialise your calculator with&
          & this property enabled.")
    end if

    call recalcGeometry(env)
    stress(:,:) = totalStress

  end subroutine getStressTensor


  !> get the gross (Mulliken projected) charges for atoms wrt neutral atoms
  subroutine getGrossCharges(env, atomCharges)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> resulting charges
    real(dp), intent(out) :: atomCharges(:)

    call recalcGeometry(env)
    atomCharges(:) = sum(q0(:, :, 1) - qOutput(:, :, 1), dim=1)

    !> Pass to the charges of the excited state if relevant
    if (isLinResp) then
      atomCharges(:) = atomCharges(:) + dQAtomEx(:)
    end if

  end subroutine getGrossCharges


  !> Sets up an external population independent electrostatic potential.
  !>
  !> Sign convention: charge of electron is considered to be positive.
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
      refExtPot%atomPot(:,1) = atomPot
    end if
    if (present(shellPot)) then
      if (.not. allocated(refExtPot%shellPot)) then
        allocate(refExtPot%shellPot(orb%mShell, nAtom, nSpin))
      end if
      @:ASSERT(all(shape(shellPot) == [orb%mShell, nAtom]))
      refExtPot%shellPot(:,:,1) = shellPot
    end if
    if (present(potGrad)) then
      if (.not. allocated(refExtPot%potGrad)) then
        allocate(refExtPot%potGrad(3, nAtom))
      end if
      @:ASSERT(all(shape(potGrad) == [3, nAtom]))
      refExtPot%potGrad(:,:) = potGrad
    end if
    tExtField = .true.

  end subroutine setExternalPotential


  !> Sets up a generator for external population dependant potentials
  subroutine setQDepExtPotProxy(extPotProxy)

    !> Generator for the external population dependant potential
    type(TQDepExtPotProxy), intent(in) :: extPotProxy

    qDepExtPot = extPotProxy

  end subroutine setQDepExtPotProxy


  !> Sets up external point charges
  subroutine setExternalCharges(chargeCoords, chargeQs, blurWidths)

    !> Coordiante of the external charges
    real(dp), intent(in) :: chargeCoords(:,:)

    !> Charges of the external point charges (sign convention: electron is negative)
    real(dp), intent(in) :: chargeQs(:)

    !> Widths of the Gaussian for each charge used for blurring (0.0 = no blurring)
    real(dp), intent(in), optional :: blurWidths(:)

    tExtChrg = .true.
    if (tForces) then
      if (.not. allocated(chrgForces)) then
        allocate(chrgForces(3, size(chargeQs)))
      end if
    end if
    call sccCalc%setExternalCharges(chargeCoords, chargeQs)

  end subroutine setExternalCharges


  !> Returns the gradient acting on the external point charges
  subroutine getExtChargeGradients(chargeGradients)

    !> Gradients
    real(dp), intent(out) :: chargeGradients(:,:)

    @:ASSERT(tForces .and. allocated(chrgForces))

    chargeGradients(:,:) = chrgForces

  end subroutine getExtChargeGradients


  !> Obtains number of atoms in the system
  function nrOfAtoms()
    integer :: nrOfAtoms

    nrOfAtoms = nAtom

  end function nrOfAtoms


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> re-evaluate the energy/forces if the geometry changes
  subroutine recalcGeometry(env)

    !> instance
    type(TEnvironment), intent(inout) :: env

    logical :: tStopDriver, tStopScc, tExitGeoOpt

    if (tLatticeChanged .or. tCoordsChanged) then
      call processGeometry(env, 1, 1, .false., tStopDriver, tStopScc, tExitGeoOpt)
      tLatticeChanged = .false.
      tCoordsChanged = .false.
    end if

  end subroutine recalcGeometry


end module dftbp_mainapi
