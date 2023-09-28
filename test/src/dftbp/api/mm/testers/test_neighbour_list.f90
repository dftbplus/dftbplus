!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Explicitly sets the neighbour list via the API and checks the result
program test_neighbour_list
  use dftbplus
  ! Only needed for the internal test system
  use testhelpers, only : writeAutotestTag
  implicit none

  ! double precision real
  integer, parameter :: dp = kind(1.0d0)

  call main_()

contains


  !! Main test routine
  !!
  !! All non-constant variables must be defined here to ensure that they are all explicitely
  !! deallocated before the program finishes.
  !!
  subroutine main_()

    type(TDftbPlus) :: dftbp
    type(TDftbPlusInput) :: input

    real(dp) :: merminEnergy, cutoff, x, dist
    real(dp) :: gradients(3, 2)
    real(dp) :: coord(3, 4), neighDist(4, 2)
    integer :: img2CentCell(4), nNeighbour(2), iNeighbour(4, 2)

    character(:), allocatable :: DftbVersion
    integer :: major, minor, patch

    call getDftbPlusBuild(DftbVersion)
    write(*,*)'DFTB+ build: ' // "'" // trim(DftbVersion) // "'"
    call getDftbPlusApi(major, minor, patch)
    write(*,"(1X,A,1X,I0,'.',I0,'.',I0)")'API version:', major, minor, patch

    call TDftbPlus_init(dftbp)

    call dftbp%getInputFromFile("dftb_in.hsd", input)
    call dftbp%setupCalculator(input)
    call TDftbPlusInput_destruct(input)

    ! setup all data for the neighbour list
    cutoff = 6.0_dp
    x = 2.5639291987021915_dp
    coord(:,:) = reshape([-x, -x, x, x, -x, -x, -x, x, -x, x, x, x], shape(coord))
    img2CentCell(:) = [2, 2, 2, 2]
    nNeighbour(:) = [4, 0]
    iNeighbour(:,:) = reshape([1, 2, 3, 4, 0, 0, 0, 0], shape(iNeighbour))
    dist = sqrt(19.721198807872987_dp)
    neighDist(:,:) = reshape([dist, dist, dist, dist, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp],&
        & shape(neighDist))

    ! set the neighbour list
    call dftbp%setNeighbourList(nNeighbour, iNeighbour, neighDist, cutoff, coord, img2CentCell)

    ! evaluate energy and forces
    call dftbp%getEnergy(merminEnergy)
    call dftbp%getGradients(gradients)

    ! clean up
    call TDftbPlus_destruct(dftbp)

    ! Write file for internal test system
    call writeAutotestTag(merminEnergy=merminEnergy, gradients=gradients)

  end subroutine main_

end program test_neighbour_list
