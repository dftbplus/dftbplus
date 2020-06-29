!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include "error.fypp"

#! (TYPE, RANK, NAME, DIM) tuple types which need to be tested for exact coherence across
#! processors
#:set EXACT_TYPES = [('real(dp)', '', 'R', '0'), ('real(dp)', '(:)', 'R', '1'),&
  & ('real(dp)', '(:,:)', 'R', '2'), ('real(dp)', '(:,:,:)', 'R', '3'),&
  & ('complex(dp)', '(:)', 'C', '1'), ('integer', '', 'I', '0'), ('integer', '(:)', 'I', '1'),&
  & ('logical', '', 'L', '0'), ('logical', '(:)', 'L', '1'), ('character(*)', '', 'S', '0')]


#! (TYPE, RANK, NAME, DIM) tuple types which need to be tested for approximate coherence across
#! processors
#:set APPROX_TYPES = [('real(dp)', '', 'R', '0'), ('real(dp)', '(:)', 'R', '1'),&
  & ('real(dp)', '(:,:)', 'R', '2'), ('real(dp)', '(:,:,:)', 'R', '3'),&
  & ('complex(dp)', '', 'C', '0'), ('complex(dp)', '(:)', 'C', '1')]

!> Contains MPI coherence tests across a comm world
module dftbp_coherence
  use dftbp_accuracy, only : dp, lc
  use dftbp_environment
#:if WITH_MPI
  use dftbp_mpifx
#:endif
  implicit none

  private
  public :: exactCoherence, toleranceCoherence, checkExactCoherence, checkToleranceCoherence

  !> Check for coherence of data across processor(s)
  interface exactCoherence
#:for _, _, NAME, DIM in EXACT_TYPES
    module procedure coherence${NAME}$${DIM}$
#:endfor
  end interface exactCoherence

  !> Check for coherence of data to a tolerance across processor(s)
  interface toleranceCoherence
#:for _, _, NAME, DIM in APPROX_TYPES
    module procedure approxCoherence${NAME}$${DIM}$
#:endfor
  end interface toleranceCoherence

  !> Check exact coherence of data across processor(s) with error handling
  interface checkExactCoherence
#:for _, _, NAME, DIM in EXACT_TYPES
    module procedure coherenceWithError${NAME}$${DIM}$
#:endfor
  end interface checkExactCoherence

  !> Check coherence of data across processor(s) to a tolerance, with error handling
  interface checkToleranceCoherence
#:for _, _, NAME, DIM in APPROX_TYPES
    module procedure approxCoherenceWithError${NAME}$${DIM}$
#:endfor
  end interface checkToleranceCoherence


contains

#:for TYPE, SHAPE, NAME, DIM in EXACT_TYPES

#:if WITH_MPI
  !> Comparison of data in global comm world
  function coherence${NAME}$${DIM}$(env, data) result(res)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Data to check for coherence
  #:if TYPE == 'character(*)' and DIM == '0'
    character(len=*), intent(in) :: data
  #:else
    ${TYPE}$, intent(in) :: data${SHAPE}$
  #:endif

    #:if DIM == '0'
    #:if TYPE == 'character(*)'
      character(len=len(data)) :: dataLocal
    #:else
      ${TYPE}$ :: dataLocal
    #:endif
    #:else
      ${TYPE}$, allocatable :: dataLocal${SHAPE}$
    #:endif

      !> Is the local data the same as the lead version?
      logical :: res

      logical :: resLocal

      resLocal = .false.
      dataLocal = data
      call mpifx_bcast(env%mpi%globalComm, dataLocal)
    #:if TYPE == 'logical'
      #:if DIM == '0'
      if (dataLocal .eqv. data) then
      #:else
      if (all(dataLocal .eqv. data)) then
      #:endif
        resLocal = .true.
      end if
    #:else
      #:if DIM == '0'
      if (dataLocal == data) then
      #:else
      if (all(dataLocal == data)) then
      #:endif
        resLocal = .true.
      end if
    #:endif
      call mpifx_allreduce(env%mpi%globalComm, resLocal, res, MPI_LAND)

  end function coherence${NAME}$${DIM}$

#:else
  !> Dummy serial comparison of data
  pure function coherence${NAME}$${DIM}$(env, data) result(res)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Data to check for coherence
  #:if TYPE == 'character(*)' and DIM == '0'
    character(len=*), intent(in) :: data
  #:else
    ${TYPE}$, intent(in) :: data${SHAPE}$
  #:endif

    logical :: res

    res = .true.

  end function coherence${NAME}$${DIM}$

#:endif

  !> Wrapper for exact coherence with error handling
  subroutine coherenceWithError${NAME}$${DIM}$(env, data, message, err)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Data to check for coherence
  #:if TYPE == 'character(*)' and DIM == '0'
    character(len=*), intent(in) :: data
  #:else
    ${TYPE}$, intent(in) :: data${SHAPE}$
  #:endif

    !> string detailing data
    character(len=*), intent(in) :: message

    !> Error code return, 0 if no problems
    integer, intent(out), optional :: err

    if (env%tAPICalculation) then
       if (.not. coherence${NAME}$${DIM}$(env, data)) then
         @:ERROR_HANDLING(err, -1, "Coherence failure in " //trim(adjustl(message))//&
             & " across nodes")
       end if
    end if

  end subroutine coherenceWithError${NAME}$${DIM}$

#:endfor

#:for TYPE, SHAPE, NAME, DIM in APPROX_TYPES

#:if WITH_MPI
  !> Comparison of data in global comm world
  function approxCoherence${NAME}$${DIM}$(env, data, tol) result(res)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Data to check for coherence
    ${TYPE}$, intent(in) :: data${SHAPE}$

    !> Tolerance for comparision, if absent use eps
    real(dp), intent(in), optional :: tol

    !> Is the local data the same as the lead version?
    logical :: res

    logical :: resLocal

    #:if DIM == '0'
      ${TYPE}$ :: dataLocal
    #:else
      ${TYPE}$, allocatable :: dataLocal${SHAPE}$
    #:endif

      real(dp) :: tol_

      if (present(tol)) then
        tol_ = tol
      else
        tol_ = epsilon(0.0_dp)
      end if

      resLocal = .false.
      dataLocal = data
      call mpifx_bcast(env%mpi%globalComm, dataLocal)
    #:if DIM == '0'
      if (abs(dataLocal - data) <= tol_) then
    #:else
      if (maxval(abs(dataLocal - data)) <= tol_) then
    #:endif
        resLocal = .true.
      end if
      call mpifx_allreduce(env%mpi%globalComm, resLocal, res, MPI_LAND)

  end function approxCoherence${NAME}$${DIM}$

#:else
  !> Dummy serial comparison of data
  pure function approxCoherence${NAME}$${DIM}$(env, data, tol) result(res)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Data to check for coherence
    ${TYPE}$, intent(in) :: data${SHAPE}$

    !> Tolerance for comparision, if absent use eps
    real(dp), intent(in), optional :: tol

    !> Is the data the same?
    logical :: res

    res = .true.

  end function approxCoherence${NAME}$${DIM}$

#:endif

  !> Wrapper for coherence within a specified tolerance, with error handling
  subroutine approxCoherenceWithError${NAME}$${DIM}$(env, data, message, tol, err)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Data to check for coherence
    ${TYPE}$, intent(in) :: data${SHAPE}$

    !> string detailing data
    character(len=*), intent(in) :: message

    !> Tolerance for comparision, if absent use eps
    real(dp), intent(in), optional :: tol

    !> Error code return, 0 if no problems
    integer, intent(out), optional :: err

    real(dp) :: tol_
    character(len=15) :: tol_str

    if (present(tol)) then
       tol_ = tol
    else
       tol_ = epsilon(0._dp)
    endif

    if (env%tAPICalculation) then
      if (.not. approxCoherence${NAME}$${DIM}$(env, data, tol_)) then
        Write(tol_str, '(E12.5)') tol_
        @:ERROR_HANDLING(err, -1, "Coherence failure in "//trim(adjustl(message))//" across nodes&
            & for a tolerance of: "//trim(adjustl(tol_str)))
      end if
    end if

  end subroutine approxCoherenceWithError${NAME}$${DIM}$

#:endfor

end module dftbp_coherence
