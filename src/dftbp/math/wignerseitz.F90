!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains routines to construct Wigner-Seitz cells.
module dftbp_math_wignerseitz

  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error
  use dftbp_math_sorting, only : index_heap_sort

  implicit none
  private


  public :: generateWignerSeitzGrid


contains


  !> Calculates a grid of lattice vectors r that fall inside (and eventually on the surface of) the
  !! Wigner-Seitz supercell centered on the origin of the Bravais superlattice with primitive
  !! translations nk(1) * a1, nk(2) * a2 and nk(3) * a3.
  !!
  !! Adapted from the Wannier90 code: w90_postw90_common, d141f9f84dcd3ac54729b9e5874dabd451684237
  !! https://github.com/wannier-developers/wannier90
  !! Giovanni Pizzi et al 2020 J. Phys.: Condens. Matter 32 165902
  !! DOI: 10.1088/1361-648X/ab51ff
  subroutine generateWignerSeitzGrid(nk, latVecs, wsVectors, wsDistances)

    !> Supercell folding coefficients (diagonal elements, i.e. MP grid)
    integer, intent(in) :: nk(:)

    !> Real-space lattice vectors of periodic geometry
    real(dp), intent(in) :: latVecs(:,:)

    !> Wigner-Seitz grid points in units of lattice vectors
    integer, intent(out), allocatable :: wsVectors(:,:)

    !> Real-space lengths of Wigner-Seitz vectors in absolute coordinates
    real(dp), intent(out), allocatable, optional :: wsDistances(:)

    !! Number of Wigner-Seitz vectors
    integer :: nWsVectors

    !! Work buffer of Wigner-Seitz grid points in units of lattice vectors
    integer :: wsVectors_(3, 20 * nk(1) * nk(2) * nk(3))

    !! Work buffer of real-space lengths in absolute coordinates
    real(dp), allocatable :: wsDistances_(:)

    !! Auxiliary variables
    integer :: n1, n2, n3, i1, i2, i3, ii, icnt, ipol, jpol, ndiff(3), nind
    real(dp) :: mindist, adot(3,3), dist(125)
    real(dp), parameter :: eps = 1e-08_dp

    !! Index array for sorting operations
    integer, allocatable :: ind(:)

    !! True, if current point belongs to Wignez-Seitz cell
    logical :: tFound

    nind = 20 * product(nk)
    allocate(wsDistances_(nind))
    if (nind < 125) then
      nind = 125
    end if
    allocate(ind(nind))

    do ipol = 1, 3
      do jpol = 1, 3
        adot(ipol, jpol) = dot_product(latVecs(:, ipol), latVecs(:, jpol))
      end do
    end do

    ! Loop over grid points r on a unit cell that is 8 times larger than a primitive supercell.
    ! On exit, nWsVectors is the total number of grids points found in the Wigner-Seitz cell.
    nWsVectors = 0
    do n1 = 0, 4 * nk(1)
      do n2 = 0, 4 * nk(2)
        do n3 = 0, 4 * nk(3)
          ! Loop over the 5^3=125 points R. R=0 corresponds to i1=i2=i3=2 or icnt=63
          icnt = 0
          do i1 = 0, 4
            do i2 = 0, 4
              do i3 = 0, 4
                icnt = icnt + 1
                ! Calculate squared distance |r-R|^2
                ndiff(1) = n1 - i1 * nk(1)
                ndiff(2) = n2 - i2 * nk(2)
                ndiff(3) = n3 - i3 * nk(3)
                dist(icnt) = 0.0_dp
                do ipol = 1, 3
                  do jpol = 1, 3
                    dist(icnt) = dist(icnt)&
                        & + real(ndiff(ipol), dp) * adot(ipol, jpol) * real(ndiff(jpol), dp)
                  end do
                end do
              end do
            end do
          end do

          ! Sort the 125 vectors R by increasing value of |r-R|^2
          call index_heap_sort(ind(1:125), dist)
          dist(:) = dist(ind(1:125))

          ! Find all the vectors R with the (same) smallest |r-R|^2.
          ! If R=0 is one of them, the current point r belongs to the Wignez-Seitz cell.
          tFound = .false.
          ii = 1
          mindist = dist(1)
          do while (abs(dist(ii) - mindist) <= eps .and. ii < 125)
            if (ind(ii) == 63) tFound = .true.
            ii = ii + 1
          end do

          if (ii == 126) ii = 125
          if (tFound) then
            nWsVectors = nWsVectors + 1
            wsVectors_(1, nWsVectors) = n1 - 2 * nk(1)
            wsVectors_(2, nWsVectors) = n2 - 2 * nk(2)
            wsVectors_(3, nWsVectors) = n3 - 2 * nk(3)
          end if
        end do
      end do
    end do

    if (nWsVectors > nind) then
      call error('Wigner-Seitz Module: Too many Wigner-Seitz points.')
    end if

    !
    do ii = 1, nWsVectors
      wsDistances_(ii) = 0.0_dp
      do ipol = 1, 3
        do jpol = 1, 3
          wsDistances_(ii) = wsDistances_(ii) + real(wsVectors_(ipol, ii), dp) * adot(ipol, jpol)&
              & * real(wsVectors_(jpol, ii), dp)
        end do
      end do
      wsDistances_(ii) = sqrt(wsDistances_(ii))
    end do

    ! Sort Wigner-Seitz vectors by increasing length
    call index_heap_sort(ind(1:nWsVectors), wsDistances_(1:nWsVectors))
    wsDistances_(1:nWsVectors) = wsDistances_(ind(1:nWsVectors))
    wsVectors_(:,:nWsVectors) = wsVectors_(:, ind(:nWsVectors))

    ! Hand over results
    wsVectors = wsVectors_(:, :nWsVectors)
    if (present(wsDistances)) wsDistances = wsDistances_(:nWsVectors)

  end subroutine generateWignerSeitzGrid

end module dftbp_math_wignerseitz
