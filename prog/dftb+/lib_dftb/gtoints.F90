!--------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations              !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                             !
!                                                                                !
!  See the LICENSE file for terms of usage and distribution.                     !
!--------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Actual implementation of the horizontal Obara--Saika recursion to evaluate
!  overlap and moment integrals over Gaussian type orbitals.
!
!  Since we are only looking at overlap related integrals here (no 1/r) we will
!  play a lot of tricks to make this integral evaluation as simple as possible.
!  First, we will work in one-dimensional primitive functions since we are allowed
!  to separate the overlap distribution in its Cartesian components and integrate
!  them separately, which saves a lot of logical overhead.
module dftbp_gtoints
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_constants, only: pi
  implicit none


  public :: TGaussFunc
  public :: shellPairOverlapIntegral, shellPairOverlapDeriv
  !public :: shellPairMomentsIntegral, shellPairMomentsDeriv
  private


  !> Allow up to d functions (l=2), for gradients l+1 is required.
  integer, parameter :: maxL = 3

  !> Maximal contraction depth of contracted Gaussian type orbitals.
  integer, parameter :: maxPrim = 7

  !> Wrapper for contracted Gaussian type orbitals.
  type :: TGaussFunc
    !> Angular momentum quantum number.
    integer :: l = -1
    !> Actual contraction depth.
    integer :: nPrim = 0
    !> Primitive exponents.
    real(dp) :: alpha(maxPrim) = 0.0_dp
    !> Contraction coefficients (containing normalization).
    real(dp) :: coeff(maxPrim) = 0.0_dp
  end type TGaussFunc


  real(dp), parameter :: sqrtpi = sqrt(pi)
  real(dp), parameter :: sqrtpi3 = sqrtpi**3

  real(dp), parameter :: small = 1.0e-9_dp

  integer, parameter :: lx(0:83) = [ &
      & 0, &
      & 1,0,0, &
      & 2,0,0,1,1,0, &
      & 3,0,0,2,2,1,0,1,0,1, &
      & 4,0,0,3,3,1,0,1,0,2,2,0,2,1,1, &
      & 5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1, &
      & 6,0,0,3,3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2]
  integer, parameter :: ly(0:83) = [ &
      & 0, &
      & 0,1,0, &
      & 0,2,0,1,0,1, &
      & 0,3,0,1,0,2,2,0,1,1, &
      & 0,4,0,1,0,3,3,0,1,2,0,2,1,2,1, &
      & 0,5,0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2, &
      & 0,6,0,3,0,3,1,0,0,1,5,5,2,0,0,2,4,4,2,1,3,1,3,2,1,4,1,2]
  integer, parameter :: lz(0:83) = [ &
      & 0, &
      & 0,0,1, &
      & 0,0,2,0,1,1, &
      & 0,0,3,0,1,0,1,2,2,1, &
      & 0,0,4,0,1,0,1,3,3,0,2,2,1,1,2, &
      & 0,0,5,0,2,0,3,2,3,0,1,0,1,4,4,3,1,1,1,2,2, &
      & 0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,1,1,4,2]

  integer, parameter :: lm(2, 0:6) = reshape(&
      & [0,0, 1,3, 4,9, 10,19, 20,34, 35,55, 56,83], shape(lm))

  real(dp), parameter :: cnorm(0:19) = [ &
      & 1.0_dp, &
      & 1.0_dp, 1.0_dp, 1.0_dp, &
      & 1.0_dp, 1.0_dp, 1.0_dp, sqrt(3.0_dp), sqrt(3.0_dp), sqrt(3.0_dp), &
      & 1.0_dp, 1.0_dp, 1.0_dp, sqrt(5.0_dp), sqrt(5.0_dp), sqrt(5.0_dp), &
      & sqrt(5.0_dp), sqrt(5.0_dp), sqrt(5.0_dp), sqrt(15.0_dp)]

contains


  !> evaluate overlap integrals over whole shell pair of contracted Gaussian
  !  type functions.
  subroutine shellPairOverlapIntegral(cgto1, cgto2, vec, dist, overlap)
    type(TGaussFunc), intent(in) :: cgto1
    type(TGaussFunc), intent(in) :: cgto2
    !> distance vector between atom 1 and 2
    real(dp), intent(in) :: vec(3)
    !> squared distance between atom 1 and 2
    real(dp), intent(in) :: dist
    real(dp), intent(out) :: overlap(:)
    real(dp) :: tmp(size(overlap, 1))
    real(dp) :: rho, orho, est, k12
    integer :: iPrim1, iPrim2

    overlap = 0.0_dp

    do iPrim1 = 1, cgto1%nPrim
      do iPrim2 = 1, cgto2%nPrim
        rho = cgto1%alpha(iPrim1) + cgto2%alpha(iPrim2)
        orho = 1.0_dp / rho
        est = dist * cgto1%alpha(iPrim1) * cgto2%alpha(iPrim2)* orho
        k12 = exp(-est) * sqrtpi3 * sqrt(orho) * orho &
            & * cgto1%coeff(iPrim1) * cgto2%coeff(iPrim2)

        call shellPrimOverlapIntegral(&
            & cgto1%l, cgto1%alpha(iPrim1), &
            & cgto2%l, cgto2%alpha(iPrim2), &
            & vec, tmp)

        overlap = overlap + tmp * k12
      end do
    end do
  end subroutine shellPairOverlapIntegral

  !> evaluate derivative of overlap integrals over whole shell pair of
  !  contracted Gaussian type functions.
  subroutine shellPairOverlapDeriv(cgto1, cgto2, vec, dist, overlap, doverlap)
    type(TGaussFunc), intent(in) :: cgto1
    type(TGaussFunc), intent(in) :: cgto2
    !> distance vector between atom 1 and 2
    real(dp), intent(in) :: vec(3)
    !> squared distance between atom 1 and 2
    real(dp), intent(in) :: dist
    real(dp), intent(out) :: overlap(:)
    real(dp), intent(out) :: doverlap(:, :)
    real(dp) :: tmp(size(overlap, 1))
    real(dp) :: dtmp(3, size(doverlap, 2))
    real(dp) :: rho, orho, est, k12
    integer :: iPrim1, iPrim2

    overlap = 0.0_dp
    doverlap = 0.0_dp

    do iPrim1 = 1, cgto1%nPrim
      do iPrim2 = 1, cgto2%nPrim
        rho = cgto1%alpha(iPrim1) + cgto2%alpha(iPrim2)
        orho = 1.0_dp / rho
        est = dist * cgto1%alpha(iPrim1) * cgto2%alpha(iPrim2)* orho
        k12 = exp(-est) * sqrtpi3 * sqrt(orho) * orho &
            & * cgto1%coeff(iPrim1) * cgto2%coeff(iPrim2)

        call shellPrimOverlapDeriv(&
            & cgto1%l, cgto1%alpha(iPrim1), &
            & cgto2%l, cgto2%alpha(iPrim2), &
            & vec, tmp, dtmp)

        overlap = overlap + tmp * k12
        doverlap = doverlap + dtmp * k12
      end do
    end do
  end subroutine shellPairOverlapDeriv


  !> evaluate overlap integrals over whole shell pair of primitive Gaussian type
  !  functions.
  subroutine shellPrimOverlapIntegral(l1, alpha1, l2, alpha2, vec, overlap)
    !> Angular momentum of shell on atom 1
    integer, intent(in) :: l1
    !> Angular momentum of shell on atom 2
    integer, intent(in) :: l2
    !> Exponent of primitive on atom 1
    real(dp), intent(in) :: alpha1
    !> Exponent of primitive on atom 2
    real(dp), intent(in) :: alpha2
    !> distance vector between atom 1 and 2
    real(dp), intent(in) :: vec(3)
    real(dp), intent(out) :: overlap(:)

    integer :: i, j, ij
    real(dp) :: tmp

    overlap = 0.0_dp

    ij = 0
    do i = lm(1, l1), lm(2, l1)
      do j = lm(1, l2), lm(2, l2)
        ij = ij+1
        call primOverlapIntegral(&
            & [lx(i), ly(i), lz(i)], alpha1, &
            & [lx(j), ly(j), lz(j)], alpha2, &
            & vec, tmp)
        overlap(ij) = tmp * cnorm(i) * cnorm(j)
      end do
    end do
  end subroutine shellPrimOverlapIntegral

  !> evaluate derivative of overlap integrals over whole shell pair of primitive
  !  Gaussian type functions.
  subroutine shellPrimOverlapDeriv(l1, alpha1, l2, alpha2, vec, overlap, doverlap)
    !> Angular momentum of shell on atom 1
    integer, intent(in) :: l1
    !> Angular momentum of shell on atom 2
    integer, intent(in) :: l2
    !> Exponent of primitive on atom 1
    real(dp), intent(in) :: alpha1
    !> Exponent of primitive on atom 2
    real(dp), intent(in) :: alpha2
    !> distance vector between atom 1 and 2
    real(dp), intent(in) :: vec(3)
    real(dp), intent(out) :: overlap(:)
    real(dp), intent(out) :: doverlap(:, :)

    integer :: i, j, ij
    real(dp) :: tmp, dtmp(3)

    overlap = 0.0_dp

    ij = 0
    do i = lm(1, l1), lm(2, l1)
      do j = lm(1, l2), lm(2, l2)
        ij = ij+1
        call primOverlapDeriv(&
            & [lx(i), ly(i), lz(i)], alpha1, &
            & [lx(j), ly(j), lz(j)], alpha2, &
            & vec, tmp, dtmp)
        overlap(ij) = tmp * cnorm(i) * cnorm(j)
        doverlap(:, ij) = dtmp * cnorm(i) * cnorm(j)
      end do
    end do
  end subroutine shellPrimOverlapDeriv


  !> primitive overlap integral over two Gaussian type orbitals.
  pure subroutine primOverlapIntegral(l1, alpha1, l2, alpha2, vec, overlap)
    !> Cartesian components of the angular momentum of primitive on atom 1
    integer, intent(in) :: l1(3)
    !> Cartesian components of the angular momentum of primitive on atom 2
    integer, intent(in) :: l2(3)
    !> Exponent of primitive on atom 1
    real(dp), intent(in) :: alpha1
    !> Exponent of primitive on atom 2
    real(dp), intent(in) :: alpha2
    !> distance vector between atom 1 and 2
    real(dp), intent(in) :: vec(3)
    !> primitive overlap integral
    real(dp), intent(out) :: overlap

    real(dp) :: rho, orho, val(3)
    real(dp) :: cfs1(0:maxL), cfs2(0:maxL), cfs(0:2*maxL)

    integer :: l

    val = 0.0_dp

    rho = alpha1 + alpha2
    orho = 1.0_dp / rho

    ! explicit loop unrolling by fypp
    #:for ixyz in (1, 2, 3)
      cfs = 0.0_dp
      cfs1(l1(${ixyz}$)) = 1.0_dp
      cfs2(l2(${ixyz}$)) = 1.0_dp
      call horizontalShift(alpha2*orho*vec(${ixyz}$), l1(${ixyz}$), cfs1)
      call horizontalShift(alpha1*orho*vec(${ixyz}$), l2(${ixyz}$), cfs2)
      call multiplyShifts(l1(${ixyz}$), cfs1, l2(${ixyz}$), cfs2, cfs)
      do l = 0, l1(${ixyz}$)+l2(${ixyz}$), 2
        if (cfs(l) > small) then
          val(${ixyz}$) = val(${ixyz}$) + cfs(l) * primOverlap1D(l, rho)
        end if
      end do
    #:endfor

    overlap = val(1) * val(2) * val(3)

  end subroutine primOverlapIntegral

  !> Derivative of primitive overlap integral over two Gaussian type orbitals.
  pure subroutine primOverlapDeriv(l1, alpha1, l2, alpha2, vec, &
      & overlap, doverlap)
    !> Cartesian components of the angular momentum of primitive on atom 1
    integer, intent(in) :: l1(3)
    !> Cartesian components of the angular momentum of primitive on atom 2
    integer, intent(in) :: l2(3)
    !> Exponent of primitive on atom 1
    real(dp), intent(in) :: alpha1
    !> Exponent of primitive on atom 2
    real(dp), intent(in) :: alpha2
    !> distance vector between atom 1 and 2
    real(dp), intent(in) :: vec(3)
    !> primitive overlap integral
    real(dp), intent(out) :: overlap
    !> derivative of primitive overlap integral w.r.t. interatomic distance
    real(dp), intent(out) :: doverlap(3)

    real(dp) :: rho, orho, tmp, val(3), grd(3)
    real(dp) :: cfs1(0:maxL), cfs2(0:maxL), cfs(0:2*maxL)
    real(dp) :: cfs1p(0:maxL), cfs1m(0:maxL), dcfs(0:2*maxL)

    integer :: l

    val = 0.0_dp
    grd = 0.0_dp

    rho = alpha1 + alpha2
    orho = 1.0_dp / rho

    ! explicit loop unrolling by fypp
    #:for ixyz in (1, 2, 3)
      cfs = 0.0_dp
      cfs1(l1(${ixyz}$)) = 1.0_dp
      cfs2(l2(${ixyz}$)) = 1.0_dp
      tmp = alpha2*orho*vec(${ixyz}$)
      call horizontalShift(tmp, l1(${ixyz}$), cfs1)
      call horizontalShift(alpha1*orho*vec(${ixyz}$), l2(${ixyz}$), cfs2)
      call multiplyShifts(l1(${ixyz}$), cfs1, l2(${ixyz}$), cfs2, cfs)

      dcfs = 0.0_dp
      cfs1p(l1(${ixyz}$)+1) = 2*alpha1
      call horizontalShift(tmp, l1(${ixyz}$)+1, cfs1p)
      if (l1(${ixyz}$) > 0) then
        cfs1m(l1(${ixyz}$)-1) = -l1(${ixyz}$)
        call horizontalShift(tmp, l1(${ixyz}$)-1, cfs1m)
        cfs1p = cfs1p + cfs1m
      end if
      call multiplyShifts(l1(${ixyz}$)+1, cfs1p, l2(${ixyz}$), cfs2, dcfs)

      do l = 0, l1(${ixyz}$)+l2(${ixyz}$)+1, 2
        if (cfs(l) > small .or. dcfs(l) > small) then
          tmp = primOverlap1D(l, rho)
          val(${ixyz}$) = val(${ixyz}$) + cfs(l) * tmp
          grd(${ixyz}$) = grd(${ixyz}$) + grd(l) * tmp
        end if
      end do
    #:endfor

    overlap = val(1) * val(2) * val(3)

    doverlap(1) = grd(1) * val(2) * val(3)
    doverlap(2) = val(1) * grd(2) * val(3)
    doverlap(3) = val(1) * val(2) * grd(3)

  end subroutine primOverlapDeriv


  !> Shift moments from one center to another.
  !
  !  Most used cases are unrolled by hand to avoid overhead from short loops.
  pure subroutine horizontalShift(ae, l, cfs)
    integer,intent(in) :: l
    real(dp), intent(in) :: ae
    real(dp), intent(inout) :: cfs(0:)
    integer :: i
    real(dp) :: aei
    select case(l)
    case(:0) ! s
      continue
    case(1) ! p
      cfs(0) = ae*cfs(1)
    case(2) ! d
      cfs(0) = ae*ae*cfs(2)
      cfs(1) = 2*ae*cfs(2)
    case(3) ! f
      cfs(0) = ae*ae*ae*cfs(3)
      cfs(1) = 3*ae*ae*cfs(3)
      cfs(2) = 3*ae*cfs(3)
    case(4) ! g
      cfs(0) = ae*ae*ae*ae*cfs(4)
      cfs(1) = 4*ae*ae*ae*cfs(4)
      cfs(2) = 6*ae*ae*cfs(4)
      cfs(3) = 4*ae*cfs(4)
    case default ! general case
      aei = 1.0_dp
      do i = l-1, 0, -1
        aei = aei * ae
        cfs(i) = aei * binomialCoefficient(l, i) * cfs(l)
      end do
    end select
  end subroutine horizontalShift

  !> Product of two horizontal shift vectors, results in the final coefficients
  !  for the expansion in primitive one-dimensional one-center Gaussian functions.
  pure subroutine multiplyShifts(l1, f1, l2, f2, prd)
    integer, intent(in) :: l1
    integer, intent(in) :: l2
    real(dp), intent(in) :: f1(0:)
    real(dp), intent(in) :: f2(0:)
    real(dp), intent(inout) :: prd(0:)
    integer :: i, j
    do i = 0, l1
      do j = 0, l2
        prd(i+j) = prd(i+j) + f1(i)*f2(j)
      end do
    end do
  end subroutine multiplyShifts

  !> binomial coefficient (n choose k)
  !
  !  n choose l can be expressed as Π(i=1→l) (n-l+i)/i,
  !  this formula should garantee that there is no problem with integer overflows.
  pure elemental function binomialCoefficient(n, k) result(b)
    integer,intent(in) :: n,k
    integer :: b
    integer :: i, l

    ! n choose n-k is always equal to n choose k,
    ! so if k greater n/2 we calculate n choose k instead of n choose n-k
    if ((2*k) > n) then
      l = n-k
    else
      l = k
    endif

    b = product([((n-l+i)/i, i = 1, k)])

  end function binomialCoefficient

  !> Evaluate one-dimensional one-center overlap integral (a Gaussian function).
  pure elemental function primOverlap1D(l, rho) result(s)
    implicit none
    integer, intent(in) :: l
    real(dp), intent(in) :: rho
    real(dp) :: s
    integer :: lh
    real(dp),parameter :: doubleFactorial(0:7) = & ! see OEIS A001147
      & [1._dp,1._dp,3._dp,15._dp,105._dp,945._dp,10395._dp,135135._dp]
    if (mod(l, 2) /= 0) then
      s = 0._dp
    else
      lh = l/2
      s = (0.5_dp/rho)**lh * doubleFactorial(lh)
    endif
  end function primOverlap1D

end module dftbp_gtoints
