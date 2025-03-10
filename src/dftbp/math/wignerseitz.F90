!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'


!> Contains routines related to Wigner-Seitz cells.
module dftbp_math_wignerseitz

  use dftbp_common_accuracy, only : dp
  use dftbp_math_sorting, only : index_heap_sort
  use dftbp_dftb_periodic, only : frac2cart
  use dftbp_common_matrixappend, only : appendToArray1d_real, appendToArray2d_int

  implicit none
  private


  public :: generateWignerSeitzGrid


contains


  !> Determines the lattice vectors pointing to all primitive unit cells that are either located
  !! inside or on the surface of the Wigner-Seitz cell of the Born-von Karman supercell, defined by
  !! an unshifted Monkhorst-Pack-type k-point grid, that is, the diagonal coefficients nKpt(i) of
  !! the supercell folding matrix.
  !!
  !! The development of this routine was inspired by the pyscf code (under the Apache-2.0 License):
  !! Sun, Q., Berkelbach, T.C., Blunt, N.S., Booth, G.H., Guo, S., Li, Z., Liu, J., McClain, J.D.,
  !! Sayfutyarova, E.R., Sharma, S., Wouters, S. and Chan, G.K.-L. (2018),
  !! PySCF: the Python-based simulations of chemistry framework.
  !! WIREs Comput Mol Sci, 8: e1340.
  !! DOI: 10.1002/wcms.1340
  subroutine generateWignerSeitzGrid(nKpt, latVecs, wsVectors, wsSearchSize, distTol)

    !> Supercell folding coefficients (diagonal elements, i.e. MP-type grid)
    integer, intent(in) :: nKpt(:)

    !> Real-space lattice vectors of periodic geometry
    real(dp), intent(in) :: latVecs(:,:)

    !> Lattice vectors inside/on the surface of the WS-supercell in fractional coordinates
    integer, intent(out), allocatable :: wsVectors(:,:)

    !> Optional index range to search for primitive unit cells inside the WS-supercell
    integer, intent(in), optional :: wsSearchSize(:)

    !> Optional tolerance for comparing real-space distances
    real(dp), intent(in), optional :: distTol

    !> Actual tolerance for comparing real-space distances
    real(dp) :: distTol_

    !! Actual index range to search for primitive unit cells inside the WS-supercell
    integer :: wsSearchSize_(3)

    !> Lattice vectors inside/on the surface of the WS-supercell in Cartesian coordinates
    real(dp), allocatable :: wsVectorsCart(:,:)

    !! Work buffer of real-space lengths in Cartesian coordinates
    real(dp), allocatable :: wsLengths(:)

    !! Index array for sorting operations
    integer, allocatable :: ind(:)

    !! Temporary distance index that indicates the central cell
    integer :: iCenter

    !! Auxiliary variables
    integer :: n1, n2, n3, i1, i2, i3, ii, jj, diff(3), tmp3Int(3)
    real(dp) :: dotProducts(3, 3), diffMat(3, 3)

    !! Temporary storage for distances
    real(dp), allocatable :: dist(:)

    @:ASSERT(all(nKpt > 0))
    @:ASSERT(size(latVecs, dim=2) == 3)

    if (present(wsSearchSize)) then
      @:ASSERT(size(wsSearchSize, dim=1) == 3)
      wsSearchSize_(:) = wsSearchSize
    else
      wsSearchSize_(:) = 2
    end if

    if (present(distTol)) then
      @:ASSERT(distTol > 0.0_dp)
      distTol_ = distTol
    else
      distTol_ = 1.0e-08_dp
    end if

    ! Cache repeatedly used dot-products
    do jj = 1, size(latVecs, dim=2)
      do ii = 1, size(latVecs, dim=2)
        dotProducts(ii, jj) = dot_product(latVecs(:, ii), latVecs(:, jj))
      end do
    end do

    ! We brute-force search over a "large enough" supercell of the primitive Wigner-Seitz supercell
    ! to find all primitive unit cells inside/on the surface of the primitive WS-supercell.
    do n1 = -wsSearchSize_(1) * nKpt(1), wsSearchSize_(1) * nKpt(1)
      do n2 = -wsSearchSize_(2) * nKpt(2), wsSearchSize_(2) * nKpt(2)
        do n3 = -wsSearchSize_(3) * nKpt(3), wsSearchSize_(3) * nKpt(3)
          if (allocated(dist)) deallocate(dist)
          do i1 = -wsSearchSize_(1), wsSearchSize_(1)
            diff(1) = n1 - i1 * nKpt(1)
            do i2 = -wsSearchSize_(2), wsSearchSize_(2)
              diff(2) = n2 - i2 * nKpt(2)
              do i3 = -wsSearchSize_(3), wsSearchSize_(3)
                diff(3) = n3 - i3 * nKpt(3)
                diffMat(:,:) = spread(diff, 1, 3) * spread(diff, 2, 3)
                call appendToArray1d_real(dist, sum(diffMat * dotProducts))
                ! Remember the index that corresponds to the center, i.e. i1=i2=i3=0
                if (i1 == 0 .and. i2 == 0 .and. i3 == 0) iCenter = size(dist, dim=1)
              end do
            end do
          end do

          if (abs(dist(iCenter) - minval(dist)) <= distTol_) then
            tmp3Int(1) = n1; tmp3Int(2) = n2; tmp3Int(3) = n3
            call appendToArray2d_int(wsVectors, tmp3Int)
          end if

        end do
      end do
    end do

    ! Copy over translations in fractional coordinates, convert to Cartesian units and calculate
    ! real-space lengths of translations found
    allocate(wsVectorsCart(3, size(wsVectors, dim=2)))
    wsVectorsCart(:,:) = real(wsVectors, dp)
    call frac2cart(wsVectorsCart, latVecs)
    wsLengths = norm2(wsVectorsCart, dim=1)

    ! Sort Wigner-Seitz vectors by increasing Euclidean norm
    allocate(ind(size(wsVectors, dim=2)))
    call index_heap_sort(ind, wsLengths)
    wsVectors(:,:) = wsVectors(:, ind)

  end subroutine generateWignerSeitzGrid

end module dftbp_math_wignerseitz
