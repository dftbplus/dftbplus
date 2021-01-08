!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Helper routines for block charges (used in on-site corrections and DFTB+U)
module dftbp_blockpothelper
  use dftbp_accuracy, only : dp
  use dftbp_commontypes, only : TOrbitals
  use dftbp_assert
  implicit none

  private
  public :: appendBlockReduced

contains

  !> Adds block (dual) charges onto end of a 1D vector of reduced charges
  subroutine appendBlockReduced(input, equiv, orb, output, isSkew)

    !> unpacked data
    real(dp), intent(in) :: input(:,:,:,:)

    !> equivalences
    integer, intent(in) :: equiv(:,:,:,:)

    !> Information about the orbitals and their angular momenta
    type(TOrbitals), intent(in) :: orb

    !> 1D array with appended data
    real(dp), intent(inout) :: output(:)

    !> is skew symmetry required
    logical, optional, intent(in) :: isSkew

    integer :: nAtom, nSpin
    integer :: iS, iOrb1, iOrb2, iAt
    logical :: tSkew

    nAtom = size(input, dim=3)
    nSpin = size(input, dim=4)
    @:ASSERT(size(input, dim=1) == orb%mOrb)
    @:ASSERT(size(input, dim=2) == orb%mOrb)
    @:ASSERT(all(shape(equiv) == (/ orb%mOrb, orb%mOrb, nAtom, nSpin /)))

    if (present(isSkew)) then
      tSkew = isSkew
    else
      tSkew = .false.
    end if

    do iS = 1, nSpin
      do iAt = 1, nAtom
        do iOrb1 = 1, orb%nOrbAtom(iAt)
          do iOrb2 = 1, orb%nOrbAtom(iAt)
            if (equiv(iOrb1, iOrb2, iAt, iS) > 0) then
              if (tSkew) then
                output(equiv(iOrb1, iOrb2, iAt, iS)) = &
                    & 0.5_dp*( input(iOrb1, iOrb2, iAt, iS) &
                    &  - input(iOrb2, iOrb1, iAt, iS) )
              else
                output(equiv(iOrb1, iOrb2, iAt, iS)) = &
                    & 0.5_dp*( input(iOrb1, iOrb2, iAt, iS) &
                    &  + input(iOrb2, iOrb1, iAt, iS) )
              end if
            end if
          end do
        end do
      end do
    end do

  end subroutine appendBlockReduced

end module dftbp_blockpothelper
