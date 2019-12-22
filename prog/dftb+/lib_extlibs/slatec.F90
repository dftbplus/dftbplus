!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Interfaces for the SLATEC routines needed in DFTB+
module dftbp_slatec
  use dftbp_accuracy, only : rdp
  implicit none
  private

  public :: zbesi, dgamit, dgamic


  interface

    !***BEGIN PROLOGUE  ZBESI
    !***PURPOSE  Compute a sequence of the Bessel functions I(a,z) for
    !            complex argument z and real nonnegative orders a=b,b+1,
    !            b+2,... where b>0.  A scaling option is available to
    !            help avoid overflow.
    !***LIBRARY   SLATEC
    !***CATEGORY  C10B4
    !***TYPE      COMPLEX (CBESI-C, ZBESI-C)
    !***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT, I BESSEL FUNCTIONS,
    !             MODIFIED BESSEL FUNCTIONS
    !***AUTHOR  Amos, D. E., (SNL)
    !***DESCRIPTION
    !
    !                    ***A DOUBLE PRECISION ROUTINE***
    !         On KODE=1, ZBESI computes an N-member sequence of complex
    !         Bessel functions CY(L)=I(FNU+L-1,Z) for real nonnegative
    !         orders FNU+L-1, L=1,...,N and complex Z in the cut plane
    !         -pi<arg(Z)<=pi where Z=ZR+i*ZI.  On KODE=2, CBESI returns
    !         the scaled functions
    !
    !            CY(L) = exp(-abs(X))*I(FNU+L-1,Z), L=1,...,N and X=Re(Z)
    !
    !         which removes the exponential growth in both the left and
    !         right half-planes as Z goes to infinity.
    !
    !         Input
    !           ZR     - DOUBLE PRECISION real part of argument Z
    !           ZI     - DOUBLE PRECISION imag part of argument Z
    !           FNU    - DOUBLE PRECISION initial order, FNU>=0
    !           KODE   - A parameter to indicate the scaling option
    !                    KODE=1  returns
    !                            CY(L)=I(FNU+L-1,Z), L=1,...,N
    !                        =2  returns
    !                            CY(L)=exp(-abs(X))*I(FNU+L-1,Z), L=1,...,N
    !                            where X=Re(Z)
    !           N      - Number of terms in the sequence, N>=1
    !
    !         Output
    !           CYR    - DOUBLE PRECISION real part of result vector
    !           CYI    - DOUBLE PRECISION imag part of result vector
    !           NZ     - Number of underflows set to zero
    !                    NZ=0    Normal return
    !                    NZ>0    CY(L)=0, L=N-NZ+1,...,N
    !           IERR   - Error flag
    !                    IERR=0  Normal return     - COMPUTATION COMPLETED
    !                    IERR=1  Input error       - NO COMPUTATION
    !                    IERR=2  Overflow          - NO COMPUTATION
    !                            (Re(Z) too large on KODE=1)
    !                    IERR=3  Precision warning - COMPUTATION COMPLETED
    !                            (Result has half precision or less
    !                            because abs(Z) or FNU+N-1 is large)
    !                    IERR=4  Precision error   - NO COMPUTATION
    !                            (Result has no precision because
    !                            abs(Z) or FNU+N-1 is too large)
    !                    IERR=5  Algorithmic error - NO COMPUTATION
    !                            (Termination condition not met)
    subroutine zbesi(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
      import rdp
      real(rdp), intent(in) :: zr
      real(rdp), intent(in) :: zi
      real(rdp), intent(in) :: fnu
      integer, intent(in) :: kode
      integer, intent(in) :: n
      real(rdp), dimension (n), intent(out) :: cyr
      real(rdp), dimension (n), intent(out) :: cyi
      integer, intent(out) :: nz
      integer, intent(out) :: ierr
    end subroutine zbesi


    !***BEGIN PROLOGUE  DGAMIT
    !***PURPOSE  Calculate Tricomi's form of the incomplete Gamma function.
    !***LIBRARY   SLATEC (FNLIB)
    !***CATEGORY  C7E
    !***TYPE      REAL(DP) (GAMIT-S, DGAMIT-D)
    !***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
    !             SPECIAL FUNCTIONS, TRICOMI
    !***AUTHOR  Fullerton, W., (LANL)
    !***DESCRIPTION
    !
    !   Evaluate Tricomi's incomplete Gamma function defined by
    !
    !   DGAMIT = X**(-A)/GAMMA(A) * integral from 0 to X of EXP(-T) *
    !              T**(A-1.)
    !
    !   for A .GT. 0.0 and by analytic continuation for A .LE. 0.0.
    !   GAMMA(X) is the complete gamma function of X.
    !
    !   DGAMIT is evaluated for arbitrary real values of A and for non-
    !   negative values of X (even though DGAMIT is defined for X .LT.
    !   0.0), except that for X = 0 and A .LE. 0.0, DGAMIT is infinite,
    !   which is a fatal error.
    !
    !   The function and both arguments are REAL(DP).
    function dgamit (a, x)
      import rdp
      real(rdp) dgamit
      real(rdp), intent(in) :: x
      real(rdp), intent(in) :: a
    end function dgamit


    !***BEGIN PROLOGUE  DGAMIC
    !***PURPOSE  Calculate the complementary incomplete Gamma function.
    !***LIBRARY   SLATEC (FNLIB)
    !***CATEGORY  C7E
    !***TYPE      DOUBLE PRECISION (GAMIC-S, DGAMIC-D)
    !***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
    !             SPECIAL FUNCTIONS
    !***AUTHOR  Fullerton, W., (LANL)
    !***DESCRIPTION
    !
    !   Evaluate the complementary incomplete Gamma function
    !
    !   DGAMIC = integral from X to infinity of EXP(-T) * T**(A-1.)  .
    !
    !   DGAMIC is evaluated for arbitrary real values of A and for non-
    !   negative values of X (even though DGAMIC is defined for X .LT.
    !   0.0), except that for X = 0 and A .LE. 0.0, DGAMIC is undefined.
    !
    !   The function and both arguments are REAL(DP).
    !
    !   A slight deterioration of 2 or 3 digits accuracy will occur when
    !   DGAMIC is very large or very small in absolute value, because log-
    !   arithmic variables are used.  Also, if the parameter A is very close
    !   to a negative INTEGER (but not a negative integer), there is a loss
    !   of accuracy, which is reported if the result is less than half
    !   machine precision.
    !
    ! REFERENCES  W. Gautschi, A computational procedure for incomplete
    !                 gamma functions, ACM Transactions on Mathematical
    !                 Software 5, 4 (December 1979), pp. 466-481.
    !               W. Gautschi, Incomplete gamma functions, Algorithm 542,
    !                 ACM Transactions on Mathematical Software 5, 4
    !                 (December 1979), pp. 482-489.
    function dgamic (a, x)
      import rdp
      real(rdp) dgamic
      real(rdp), intent(in) :: x
      real(rdp), intent(in) :: a
    end function dgamic

  end interface

end module dftbp_slatec
