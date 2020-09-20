!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#:set FLAVOURS = [('Real', 'real', 'real'), ('Cmplx', 'complex', 'cmplx')]

!> Contains routines to calculate matrix determinants
module dftbp_determinant
  use dftbp_accuracy
  use dftbp_assert
  use dftbp_lapackroutines
#:if WITH_SCALAPACK
  use dftbp_mpifx
  use dftbp_scalapackfx
#:endif
  implicit none
  private

  public :: det

  interface det
  #:for SUFFIX, _, _ in FLAVOURS
    module procedure det${SUFFIX}$
  #:if WITH_SCALAPACK
      module procedure detScaLAPACK${SUFFIX}$
  #:endif
  #:endfor
  end interface det

contains

#:for SUFFIX, TYPE, CONVERT in FLAVOURS

  !> Determinant of a real matrix, matrix destroyed in process
  function det${SUFFIX}$(A) result(det)

    !> The matrix
    ${TYPE}$(dp), intent(inout) :: A(:,:)

    !> resulting determinant
    ${TYPE}$(dp) :: det

    integer, allocatable  :: ipiv(:)
    integer :: ii, n, exponent

    n = minval(shape(A))
    allocate(ipiv(n))

    call getrf(A,ipiv)

    det = ${CONVERT}$(1, kind=dp)
    exponent = 0
    do ii = 1, n
      if (ipiv(ii) /= ii) then
        det = -det * A(ii,ii)
      else
        det = det * A(ii,ii)
      end if
      if (det == 0.0_dp) then
        return
      end if
      do while (abs(det) > 2.0_dp)
        det = det / 2.0_dp
        exponent = exponent + 1
      end do
      do while (abs(det) < 0.5_dp)
        det = det * 2.0_dp
        exponent = exponent - 1
      end do
    end do
    det = det * 2.0_dp ** exponent

  end function det${SUFFIX}$

#:if WITH_SCALAPACK

  !> Determinant of a real matrix, matrix destroyed in process
  function detScaLAPACK${SUFFIX}$(A, descA, grid, myComm) result(det)

    !> The matrix
    ${TYPE}$(dp), intent(inout) :: A(:,:)

    !> Dense descriptor
    integer, intent(in) :: descA(DLEN_)

    !> BLACS grid involved in calculation
    type(blacsgrid), intent(in) :: grid

    !> Communicator for the region involved in the BLACS grid
    type(mpifx_comm), intent(in) :: myComm

    !> resulting determinant
    ${TYPE}$(dp) :: det

    integer, allocatable  :: ipiv(:)
    integer :: ii, jj, iLoc, jLoc, mm, nn
    logical :: tDiagBlock, tAnyDiag
    ${TYPE}$(dp) :: detLocal
    ${TYPE}$(dp), allocatable :: detBuffer(:)
    integer :: expLocal
    integer, allocatable :: expBuffer(:)

    @:ASSERT(grid%nProc == myComm%size)

    if (grid%iproc /= -1) then
      mm = descA(M_)
      nn = descA(N_)

      allocate(detBuffer(grid%nProc))
      allocate(expBuffer(grid%nProc))
      detBuffer = 0.0_dp
      expbuffer = 0.0_dp

      allocate(ipiv(min(mm,nn)))
      ipiv = 0
      call scalafx_pgetrf(A,descA,ipiv)

      ! note, this includes under-/over-flow protection similar to LINPACK routine dgedi.f
      detLocal = 1.0_dp
      expLocal = 0
      tAnyDiag = .false.
      lpLocal: do ii = 1, size(A,dim=2)

        ! Look for diagonal blocks
        jj = scalafx_indxl2g(ii, descA(NB_), grid%mycol, descA(CSRC_), grid%ncol)
        call scalafx_islocal(grid, descA, jj, jj, tDiagBlock, iLoc, jLoc)

        tAnyDiag = tAnyDiag .or. tDiagBlock

        if (tDiagBlock) then
          if (jj /= ipiv(ii)) then
            detLocal = -detLocal * A(iLoc,jLoc)
          else
            detLocal = detLocal * A(iLoc,jLoc)
          end if

          if (detLocal == 0.0_dp) then
            exit lpLocal
          end if

          do while (abs(detLocal) > 2)
            detLocal = detLocal / 2.0_dp
            expLocal = expLocal + 1
          end do
          do while (abs(detLocal) < 0.5_dp)
            detLocal = detLocal * 2.0_dp
            expLocal = expLocal - 1
          end do
        end if

      end do lpLocal

      if (tAnyDiag) then
        detBuffer(grid%iProc+1) = detLocal
        expBuffer(grid%iProc+1) = expLocal
      else
        ! node did not have any diagonal elements, so does not contribute to det
        detBuffer(grid%iProc+1) = 1.0_dp
        expBuffer(grid%iProc+1) = 0.0_dp
      end if

      ! now product the full det from the sub-processes
      call mpifx_allreduceip(myComm, detBuffer, MPI_SUM)
      call mpifx_allreduceip(myComm, expBuffer, MPI_SUM)

      detLocal = ${CONVERT}$(1, kind=dp)
      expLocal = 0
      lpTotal: do ii = 1, grid%nProc
        detLocal = detLocal * detBuffer(ii)
        expLocal = expLocal + expBuffer(ii)
        if (detLocal == 0.0_dp) then
          exit lpTotal
        end if
        do while (abs(detLocal) > 2)
          detLocal = detLocal / 2.0_dp
          expLocal = expLocal + 1
        end do
        do while (abs(detLocal) < 0.5_dp)
          detLocal = detLocal * 2.0_dp
          expLocal = expLocal - 1
        end do
      end do lpTotal

      det = detLocal * 2.0_dp ** expLocal

    end if

  end function detScaLAPACK${SUFFIX}$

#:endif

#:endfor

end module dftbp_determinant
