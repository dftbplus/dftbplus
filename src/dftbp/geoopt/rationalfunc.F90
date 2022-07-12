!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Implementation of a rational function optimization procedure.
!>
!> The optimization problem is solved by determining the lowest eigensolution of the
!> augmented Hessian matrix. The eigenvector provides the optimal displacement for the
!> optimization. This rational function implementation will update the approximate Hessian
!> matrix using a BFGS-like update procedure.
module dftbp_geoopt_rationalfunc
  use dftbp_common_accuracy, only : dp
  use dftbp_geoopt_optimizer, only : TOptimizer, TOptimizerInput
  use dftbp_math_blasroutines, only : spmv
  use dftbp_math_eigensolver, only : syev => heev
  implicit none

  private
  public :: TRationalFuncInput, TRationalFunc, TRationalFunc_init


  !> Input for the rational function optimizer
  type, extends(TOptimizerInput) :: TRationalFuncInput

    !> Lower limit of diagonal Hessian elements
    real(dp) :: diagLimit = 1.0e-2_dp

  end type TRationalFuncInput


  !> Rational function optimization driver
  type, extends(TOptimizer) :: TRationalFunc

    !> Number of variables to optimize
    integer :: nvar

    !> Lower limit of diagonal Hessian elements
    real(dp) :: diagLimit

    !> Last gradient
    real(dp), allocatable :: gLast(:)

    !> Approximate Hessian matrix, will be updated every displacement step
    real(dp), allocatable :: hess(:)

    !> Scratch space for augmented Hessian using to solve rational function problem
    real(dp), allocatable :: aaug(:)

    !> Lowest eigenvector of rational function problem
    real(dp), allocatable :: uaug(:)

  contains

    !> Calculate displacement from gradient
    procedure :: step

  end type TRationalFunc

contains


  !> Create new rational function optimization driver
  subroutine TRationalFunc_init(this, input, nVar)

    !> Instance of the optimizer
    type(TRationalFunc), intent(out) :: this

    !> Input for the rational function optimizer
    type(TRationalFuncInput), intent(in) :: input

    !> Number of variables to optimize
    integer, intent(in) :: nVar

    integer :: ii, nvar1, npvar, npvar1

    this%nvar = nVar
    this%diagLimit = input%diagLimit

    nvar1  = this%nvar+1
    npvar = this%nvar*nvar1/2
    npvar1 = nvar1*(1+nvar1)/2
    allocate(this%gLast(nVar), source=0.0_dp)
    allocate(this%hess(npvar), source=0.0_dp)
    allocate(this%aaug(npvar1), source=0.0_dp)
    allocate(this%uaug(nvar1), source=0.0_dp)

    this%hess(:) = 0.0_dp
    do ii = 1, this%nvar
      this%hess(ii*(1+ii)/2) = 1.0_dp
    end do
  end subroutine TRationalFunc_init


  !> Calculate displacement from gradient
  subroutine step(this, val, grad, displ)

    !> Instance of geometry optimization driver
    class(TRationalFunc), intent(inout) :: this

    !> Current function value
    real(dp), intent(in) :: val

    !> Current gradient
    real(dp), intent(in) :: grad(:)

    !> Next displacement step
    real(dp), intent(out) :: displ(:)

    real(dp) :: eaug
    integer :: ii, nvar1, npvar, npvar1
    logical :: fail

    nvar1  = this%nvar+1
    npvar = this%nvar*nvar1/2
    npvar1 = nvar1*(1+nvar1)/2

    displ(:) = this%uaug(:this%nvar) / this%uaug(nvar1)
    call bfgsUpdate(grad, this%gLast, displ, this%diagLimit, this%hess)

    this%aaug(:) = [this%hess, grad, 0.0_dp]
    this%uaug(:) = [-grad, 1.0_dp]
    this%uaug(:) = this%uaug(:) / norm2(this%uaug)
    call davidson(nvar1, sqrt(epsilon(1.0_dp)), this%aaug, this%uaug, eaug, fail)
    displ(:) = this%uaug(:this%nvar) / this%uaug(nvar1)
    this%gLast(:) = grad

  end subroutine step


  !> Davidson iterative eigenvalue solver, solves for the first (lowest) eigenvalue only
  subroutine davidson(n, crite, Hp, C, e, fail)

    !> Eigenvalue convergence threshold
    real(dp), intent(in) :: crite

    !> Matrix to be diagonalized
    real(dp), intent(in) :: Hp(:)

    !> Eigenvectors
    real(dp), intent(inout) :: C(:)

    !> Eigenvalues
    real(dp), intent(out) :: e

    !> Failed to solve eigenvalue problem
    logical, intent(out) :: fail

    integer, parameter :: maxiter = 100
    integer, parameter :: initial_dyn_array_size = 10
    integer :: n, m, jold, ij, i, j, iold
    logical :: converged
    real(dp), allocatable :: lun1(:, :), lun2(:, :)
    real(dp) valn, uim, s, denerg
    real(dp), allocatable :: adiag(:), vecf1(:), vecf2(:), w(:)
    real(dp), allocatable :: Uaug(:, :), d(:)
    real(dp), allocatable :: av(:)

    fail = .true.

    n = size(C)
    allocate(adiag(n), vecf1(n), vecf2(n), w(n), av(maxiter*(maxiter+1)/2))

    allocate(lun1(n, initial_dyn_array_size), source=0.0_dp)
    allocate(lun2(n, initial_dyn_array_size), source=0.0_dp)

    ! H * C for initialization
    call mwrite(lun1, C, 1)
    call spmv(HP, C, vecf2, uplo='u')
    call mwrite(lun2, vecf2, 1)

    e = 0
    valn = 0
    converged = .false.

    do i = 1, n
      adiag(i) = HP(i*(i+1)/2)
    end do

    av(1) = dot_product(C, vecf2)
    ! done

    do m = 1, maxiter-1
      allocate(Uaug(m, m), d(m))

      do i = 1, m
        do j = 1, i
          ij = i*(i-1)/2+j
          Uaug(j, i) = av(ij)
          Uaug(i, j) = av(ij)
        end do
      end do
      call syev(Uaug, d, uplo='u', jobz='v')
      valn = d(1)

      vecf1 = 0.0_dp
      do i = 1, m
        w(:) = lun1(:, i)
        uim = Uaug(i, 1)
        vecf1(:) = uim * w + vecf1
      end do

      vecf2 = -valn * vecf1
      do i = 1, m
        w(:) = lun2(:, i)
        uim = Uaug(i, 1)
        vecf2(:) = uim * w + vecf2
      end do
      deallocate(Uaug, d)

      C(:) = vecf1

      vecf1(:) = vecf2/(valn-adiag)

      denerg = abs(valn - e)
      converged = abs(valn - e) < crite

      if (converged) then
        fail = .false.
        exit
      end if

      iold = m
      do jold = 1, iold
        w(:) = lun1(:, jold)
        s = -dot_product(w, vecf1)
        vecf1(:) = s * w + vecf1
      end do
      s = dot_product(vecf1, vecf1)
      if (s > 0.00000001_dp) then
        s = 1.0_dp /sqrt(s)
        vecf1 = vecf1 * s
        iold = iold + 1
        call mwrite(lun1, vecf1, jold)
      else
        fail = .false.
        exit
      end if

      ! H * C
      call spmv(HP, vecf1, vecf2)

      call mwrite(lun2, vecf2, m+1)

      do jold = 1, m
        w(:) = lun1(:, jold)
        av(m*(m+1)/2 + jold) = dot_product(w, vecf2)
      end do
      av((m+1)*(m+2)/2) = dot_product(vecf2, vecf1)

      ! increase expansion space and iterate further
      e = valn
    end do

  contains

    !> Write vector to storage, optionally reallocate storage to place new vectr
    pure subroutine mwrite(iwo, v, irec)
      !> Storage matrix
      real(dp), intent(inout), allocatable :: iwo(:, :)
      !> Vector to be stored
      real(dp), intent(in)  :: v(:)
      !> Record to store vector in
      integer,  intent(in)  :: irec

      real(dp), allocatable :: tmp(:, :)
      integer :: d2, dn, n
      n = size(iwo, 1)
      d2 = size(iwo, 2)
      if (irec > d2) then
        dn = d2 + d2/2 + 1
        allocate(tmp(n, dn))
        tmp(:, :d2) = iwo
        deallocate(iwo)
        call move_alloc(tmp, iwo)
      endif
      iwo(:, irec) = v
    end subroutine mwrite

  end subroutine davidson


  !> Perform BFGS-like update of packed Hessian matrix
  subroutine bfgsUpdate(gcurr, glast, displ, diagLimit, hess)

    !> Current gradient
    real(dp), intent(in) :: gcurr(:)

    !> Gradient from last optimization step
    real(dp), intent(in) :: glast(:)

    !> Displacement from last to current
    real(dp), intent(in) :: displ(:)

    !> Lower limit of diagonal Hessian elements
    real(dp), intent(in) :: diagLimit

    !> Approximate Hessian matrix
    real(dp), intent(inout) :: hess(:)

    integer  :: i, j, ij, ii, nn
    real(dp), allocatable :: svec(:), tvec(:)
    real(dp) :: ddtd, dds, ooddtd, oodds, sdds, tddtd
    real(dp), parameter :: thrs = 100*epsilon(0.0_dp)

    nn = size(gcurr)
    allocate(svec(nn), tvec(nn), source=0.0_dp)

    svec(:) = gcurr - glast
    call spmv(hess, displ, tvec)

    ! calculate scalar dxdx and jtdx
    ddtd = dot_product(tvec, displ)
    dds = dot_product(svec, displ)
    ooddtd = 1.0_dp / ddtd
    oodds = 1.0_dp / dds

    if (dds > thrs .and. ddtd > thrs) then
      !$omp parallel do default(none) shared(nn, oodds, ooddtd, svec, tvec, hess) &
      !$omp private(i, j, ii, ij, sdds, tddtd)
      do i = 1, nn
        ii = i*(i-1)/2
        sdds  = svec(i)*oodds
        tddtd = tvec(i)*ooddtd
        do j = 1, i
          ij = ii + j
          hess(ij) = hess(ij) + svec(j)*sdds - tvec(j)*tddtd
        end do
      end do
    end if

    ! limit diagonal to (0.01 slightly better than 0.001)
    do i = 1, nn
      ij = i*(i+1)/2
      if (abs(hess(ij)) < diagLimit) hess(ij) = diagLimit
    end do

  end subroutine bfgsUpdate


end module dftbp_geoopt_rationalfunc
