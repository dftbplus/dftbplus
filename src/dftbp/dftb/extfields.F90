!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> External fields
module dftbp_dftb_extfields
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_potentials, only : TPotentials
  use dftbp_io_message, only : warning
  implicit none

  public :: TEField, TElecFieldInput, addUpExternalField

  !> External field data structure
  type :: TEField

    ! homogeneous electric field

    !> field strength
    real(dp), allocatable :: EFieldStrength

    !> field direction
    real(dp) :: EfieldVector(3)

    !> Field is time dependent
    logical :: isTDEfield

    !> angular frequency
    real(dp) :: EfieldOmega

    !> phase of field at step 0
    integer :: EfieldPhase

    ! Derived properties from E field

    !> Instantaneous electric field
    real(dp) :: EField(3)

    !> Instantaneous magnitude of electric field
    real(dp) :: absEField

    ! potentials at atomic sites

    !> Atoms at which there is a external potential
    integer, allocatable :: atomicSites(:)

    !> Potential at the atomic sites
    real(dp), allocatable :: atomicPotential(:)

    !> is this an onsite (net charge related, so diagonal elements of hamiltonian) or not
    !> (gross/Mulliken charge related, so incluing off diagonal coupling)
    logical, allocatable :: atomicOnSite(:)

  end type TEField


  !> External homogeneous electric field
  type TElecFieldInput

    !> time dependent field in MD
    logical :: isTDEfield = .false.

    !> strength
    real(dp) :: EFieldStrength = 0.0_dp

    !> direction
    real(dp) :: EfieldVector(3)

    !> frequency of time dependent field
    real(dp) :: EfieldOmega

    !> relative phase of field
    integer :: EfieldPhase = 0

  end type TElecFieldInput

contains

  !> Sets up external electric or other fields
  subroutine addUpExternalField(eField, tPeriodic, neighbourList, nNeighbourSK, iCellVec,&
      & cellVec, deltaT, iGeoStep, coord0Fold, coord, potential)

    !> Whether an external field is present
    type(TEfield), intent(inout), allocatable :: eField

    !> Is this a periodic geometry
    logical, intent(in) :: tPeriodic

    !> Atomic neighbours
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for unit cells
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Time step in MD
    real(dp), intent(in) :: deltaT

    !> Number of the geometry step
    integer, intent(in) :: iGeoStep

    !> Atomic coordinates in central cell
    real(dp), allocatable, intent(in) :: coord0Fold(:,:)

    !> all coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Potential, including external potentials
    type(TPotentials), intent(inout) :: potential

    integer :: nAtom, iAt1, iAt2, iNeigh, ii
    logical :: isBoundaryCrossed

    if (allocated(eField)) then
      nAtom = size(nNeighbourSK)

      if (allocated(eField%EFieldStrength)) then
        eField%Efield(:) = eField%EFieldStrength * eField%EfieldVector
        if (eField%isTDEField) then
          eField%Efield(:) = eField%Efield * sin(eField%EfieldOmega * deltaT *&
              & real(iGeoStep + eField%EfieldPhase, dp))
        end if
        eField%absEfield = sqrt(sum(eField%Efield**2))
        if (tPeriodic) then
          isBoundaryCrossed = .false.
          do iAt1 = 1, nAtom
            do iNeigh = 1, nNeighbourSK(iAt1)
              iAt2 = neighbourList%iNeighbour(iNeigh, iAt1)
              ! overlap between atom in central cell and non-central cell
              if (iCellVec(iAt2) /= 0) then
                ! component of electric field projects onto vector between cells
                if (abs(dot_product(cellVec(:, iCellVec(iAt2)), eField%EfieldVector))&
                    & > epsilon(1.0_dp)) then
                  isBoundaryCrossed = .true.
                end if
              end if
            end do
          end do
          if (isBoundaryCrossed) then
            call warning("Interactions between atoms cross the saw-tooth discontinuity in the&
                & electric field")
          end if
          do iAt1 = 1, nAtom
            potential%extAtom(iAt1,1) = potential%extAtom(iAt1,1)&
                & + dot_product(coord0Fold(:, iAt1), eField%Efield)
          end do
        else
          do iAt1 = 1, nAtom
            potential%extAtom(iAt1,1) = potential%extAtom(iAt1,1) + dot_product(coord(:, iAt1),&
                & eField%Efield)
          end do
        end if
        potential%extGrad(:,:) = potential%extGrad + spread(eField%EField, 2, nAtom)
      end if

      if (allocated(eField%atomicSites)) then

        do ii = 1, size(eField%atomicSites)
          iAt1 = eField%atomicSites(ii)
          if (eField%atomicOnSite(ii)) then
            potential%extOnSiteAtom(iAt1,1) = potential%extOnSiteAtom(iAt1,1)&
                & + eField%atomicPotential(ii)
          else
            potential%extAtom(iAt1,1) = potential%extAtom(iAt1,1) + eField%atomicPotential(ii)
          end if
        end do

      end if

    end if

  end subroutine addUpExternalField


end module dftbp_dftb_extfields
