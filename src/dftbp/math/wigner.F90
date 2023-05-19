!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains routines to construct Wigner-Seitz cells.
module dftbp_math_wigner

  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error
  use dftbp_math_sorting, only : index_heap_sort

  implicit none
  private


  public :: generateWignerSeitzGrid, isInWignerSeitzGrid


contains

  !> Calculates a grid of points that fall inside of (and eventually on the surface of) the
  !! Wigner-Seitz supercell centered on the origin of the Bravais lattice with primitive
  !! translations nk1*a_1+nk2*a_2+nk3*a_3
  !!
  !! BUG FIX: in the case of the tetragonal cell with c>a (LSCO) the WS points are correctly
  !! determined, but there are a few points with the wrong degeneracies. To avoid this I search for
  !! points in the -2:2 replicas (5^2 replicas). I had the same problem in createkmap for the g0vec
  !! shift, and also there I have fixed it by extending the replicas to -2:2 instead of -1:1.
  subroutine generateWignerSeitzGrid(nk, latVecs, wsVectors, tExcludeSurface, wsDegeneracy,&
      & wsDistances)

    !> Supercell folding coefficients (diagonal elements)
    integer, intent(in) :: nk(:)

    !> Lattice vectors of (periodic) geometry
    real(dp), intent(in) :: latVecs(:,:)

    !> Wigner-Seitz grid points in units of lattice vectors
    integer, intent(out), allocatable :: wsVectors(:,:)

    !> True, if surface should be excluded (default: true)
    logical, intent(in), optional :: tExcludeSurface

    !> Degeneracy of the i-th Wigner-Seitz vector (weight is 1 / wsDegeneracy(iVec))
    integer, intent(out), allocatable, optional :: wsDegeneracy(:)

    !> Real-space lengths in absolute units
    real(dp), intent(out), allocatable, optional :: wsDistances(:)

    !! Number of Wigner-Seitz vectors
    integer :: nrr

    !! Work buffer of Wigner-Seitz grid points in units of lattice vectors
    integer :: wsVectors_(3, 20 * nk(1) * nk(2) * nk(3))

    !! Work buffer of degeneracy of the i-th Wigner-Seitz vector
    integer, allocatable :: wsDegeneracy_(:)

    !! Work buffer of real-space lengths in absolute units
    real(dp), allocatable :: wsDistances_(:)

    !> True, if surface should be excluded
    logical :: tExcludeSurface_

    !! Auxiliary variables
    integer :: n1, n2, n3, i1, i2, i3, ii, ipol, jpol, ndiff(3), nind
    real(dp) :: tot, mindist, adot(3,3), dist(125)
    real(dp), parameter :: eps = 1e-08_dp

    !! Index array for sorting operations
    integer, allocatable :: ind(:)

    !! True, if current point belongs to Wignez-Seitz cell
    logical :: tFound

    !! Number of highest equal distances
    integer :: nEqualDists

    if (present(tExcludeSurface)) then
      tExcludeSurface_ = tExcludeSurface
    else
      tExcludeSurface_ = .true.
    end if

    nind = 20 * product(nk)
    allocate(wsDegeneracy_(nind))
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

    ! Loop over grid points r on a unit cell that is 8 times larger than a
    ! primitive supercell. In the end nrr contains the total number of grids
    ! points that have been tFound in the Wigner-Seitz cell
    nrr = 0
    do n1 = 0, 4 * nk(1)
      do n2 = 0, 4 * nk(2)
        do n3 = 0, 4 * nk(3)
          ! Loop over the 27 points R. R=0 corresponds to i1=i2=i3=1, or icnt=14
          ! Loop over the 5^3 = 125 points R. R=0 corresponds to i1=i2=i3=2, or icnt=63
          ii = 0
          do i1 = 0, 4
            do i2 = 0, 4
              do i3 = 0, 4
                ii = ii + 1
                ! Calculate distance squared |r-R|^2
                ndiff(1) = n1 - i1 * nk(1)
                ndiff(2) = n2 - i2 * nk(2)
                ndiff(3) = n3 - i3 * nk(3)
                dist(ii) = 0.0_dp
                do ipol = 1, 3
                  do jpol = 1, 3
                    dist(ii) = dist(ii) + real(ndiff(ipol), dp) * adot(ipol, jpol)&
                        & * real(ndiff(jpol), dp)
                  end do
                end do
              end do
            end do
          end do

          ! Sort the  27 vectors R by increasing value of |r-R|^2
          ! Sort the 125 vectors R by increasing value of |r-R|^2

          ! NOTA BENE: hpsort really sorts the dist vector
          ! while the original subroutine by MVS did not. Therefore,
          ! dist(ind(i)) of the original version is here replaced by
          ! dist(i), while ind(i) is kept.

          call index_heap_sort(ind(1:125), dist)
          dist(:) = dist(ind(1:125))

          ! Find all the vectors R with the (same) smallest |r-R|^2;
          ! if R=0 is one of them, then the current point r belongs to
          ! Wignez-Seitz cell => set tFound to true
          tFound = .false.
          ii = 1
          mindist = dist(1)
          do while (abs(dist(ii) - mindist) .le. eps .and. ii .lt. 125)
            if (ind(ii) == 63) tFound = .true.
            ii = ii + 1
          end do

          if (ii .eq. 126) ii = 125
          if (tFound) then
            nrr = nrr + 1
            wsDegeneracy_(nrr) = ii - 1
            wsVectors_(1, nrr) = n1 - 2 * nk(1)
            wsVectors_(2, nrr) = n2 - 2 * nk(2)
            wsVectors_(3, nrr) = n3 - 2 * nk(3)
          end if
        end do
      end do
    end do

    ! Check the "sum rule"
    tot = 0.0_dp
    do ii = 1, nrr
      tot = tot + 1.0_dp / real(wsDegeneracy_(ii), dp)
    end do

    if (abs(tot - real(nk(1) * nk(2) * nk(3), dp)) .gt. eps) then
      call error('Wigner-Seitz Module: Weights do not add up to nk(1) * nk(2) * nk(3)')
    end if

    !@ JN it happens in 2d and 1d systems with small course grids. I've changed to 20**3
    !@ could calculate the max number of elements at the beginning
    ! Hopefully this will never happen, i.e., I think 2 * nk(1) * nk(2) * nk(3) is
    ! an upper bound to the number of lattice points in (or on
    ! the surface of) the Wigner-Seitz supercell
    if (nrr > nind) then
      call error('Wigner-Seitz Module: Too many Wigner-Seitz points, try to increase the bound&
          & 20 * product(nk)')
    end if

    ! Now sort the wigner-seitz vectors by increasing magnitude
    do ii = 1, nrr
      wsDistances_(ii) = 0.0_dp
      do ipol = 1, 3
        do jpol = 1, 3
          wsDistances_(ii) = wsDistances_(ii) + real(wsVectors_(ipol, ii), dp) * adot(ipol, jpol)&
              & * real(wsVectors_(jpol, ii), dp)
        end do
      end do
      wsDistances_(ii) = sqrt(wsDistances_(ii))
    end do

    call index_heap_sort(ind(1:nrr), wsDistances_(1:nrr))
    wsDistances_(1:nrr) = wsDistances_(ind(1:nrr))

    ! Now wsDistances_ is already sorted, but we still have to sort wsVectors_ and wsDegeneracy_
    wsDegeneracy_(:nrr) = wsDegeneracy_(ind(:nrr))
    wsVectors_(:,:nrr) = wsVectors_(:, ind(:nrr))

    if (tExcludeSurface_) then
      ! Remove most outer shell, to be on the save side
      nEqualDists = 1
      do ii = 1, nrr - 1
        if (abs(wsDistances_(nrr) - wsDistances_(nrr - ii)) < eps) then
          nEqualDists = nEqualDists + 1
        else
          exit
        end if
      end do
      if (nrr - nEqualDists > 0) nrr = nrr - nEqualDists
    end if

    ! Hand over results
    wsVectors = wsVectors_(:, :nrr)
    if (present(wsDegeneracy)) wsDegeneracy = wsDegeneracy_(:nrr)
    if (present(wsDistances)) wsDistances = wsDistances_(:nrr)

  end subroutine generateWignerSeitzGrid


  !> Searches given integer tripel in array of Wigner-Seitz grid points.
  pure function isInWignerSeitzGrid(vec, wsArray) result(tFound)

    !> Vector in relative coordinates to search for
    integer, intent(in) :: vec(:)

    !> Array of Wigner-Seitz grid points/vectors
    integer, intent(in) :: wsArray(:,:)

    !> True, if given vector was found in array
    logical :: tFound

    !! Iterates over all Wigner-Seitz vectors
    integer :: ii

    tFound = .false.

    do ii = 1, size(wsArray, dim=2)
      if (all(wsArray(:, ii) == vec)) then
        tFound = .true.
        return
      end if
    end do

  end function isInWignerSeitzGrid

end module dftbp_math_wigner
