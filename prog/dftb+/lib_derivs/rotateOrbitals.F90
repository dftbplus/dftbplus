!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!#:include 'common.fypp'

!> Module containing routines to make linear combinations of orbitals for degenerate perturbation
!> from a hermitian/symmetric matrix
module dftbp_degeneratePerturb
  use dftbp_accuracy, only : dp
  implicit none
  private

  public :: degenerateUnitary

interface degenerateUnitary
  module procedure realTransform
  module procedure cmplxTransform
end interface

contains

subroutine realTransform(matrixToProcess, ei, tol)

!> Matrix elements between (potentially degenerate) orbitals
real(dp), intent(inout) :: matrixToProcess(:,:)

!> Eigenvalues of states
real(dp), intent(in) :: ei(:)

!> Optional tolerance for comparisions between states
real(dp), intent(in), optional :: tol

real(dp), allocatable :: U(:,:)
integer :: nOrb

nOrb = size(ei)

allocate(U(nOrb,nOrb))

call realUnitary(U, matrixToProcess, ei, tol)

matrixToProcess(:,:) = matmul(transpose(U), matmul(matrixToProcess, U))

end subroutine


!> Transform degenerate states into combination that is orthogonal under the action of the matrix
subroutine realUnitary(U, matrixToProcess, ei, tol)

!> Unitary matrix to transform orbitals
real(dp), intent(out) :: U(:,:)

!> Matrix elements between (potentially degenerate) orbitals
real(dp), intent(in) :: matrixToProcess(:,:)

!> Eigenvalues of states
real(dp), intent(in) :: ei(:)

!> Tolerance for comparisions between states
real(dp), intent(in), optional :: tol

integer :: ii, nOrb, nGrp, iGrp, maxRange, iRange
integer, allocatable :: degeneracies(:), degenerateRange(:,:)
real(dp), allocatable :: subBlock(:,:), eigenvals(:)

U(:,:) = 0.0_dp
nOrb = size(ei)

do ii = 1, nOrb
  U(ii,ii) = 1.0_dp
end do


allocate(degenerateRange(2,nOrb))
call degeneracyRanges(degenerateRange, nGrp, Ei, tol)

maxRange = maxval(degenerateRange(2,:) - degenerateRange(1,:)) + 1

if (maxRange == 1) then
return
end if

allocate(subBlock(maxRange,maxRange))
allocate(eigenvals(maxRange))

do iGrp = 1, nGrp
subBlock(:,:) = 0.0_dp
eigenvals(:) = 0.0_dp
iRange = degenerateRange(2, iGrp) - degenerateRange(1, iGrp) + 1
subBlock(:iRange, :iRange) = matrixToProcess(degenerateRange(1, iGrp):degenerateRange(2, iGrp),&
& degenerateRange(1, iGrp):degenerateRange(2, iGrp))

call heev(subBlock(:iRange, :iRange), eigenvals, 'L', jobs='V')

U(degenerateRange(1, iGrp):degenerateRange(2, iGrp),&
& degenerateRange(1, iGrp):degenerateRange(2, iGrp)) = subBlock(:iRange, :iRange)

end do

end subroutine


subroutine cmplxTransform(matrixToProcess, ei, tol)

!> Matrix elements between (potentially degenerate) orbitals
complex(dp), intent(inout) :: matrixToProcess(:,:)

!> Eigenvalues of states
real(dp), intent(in) :: ei(:)

!> Optional tolerance for comparisions between states
real(dp), intent(in), optional :: tol

complex(dp), allocatable :: U(:,:)
integer :: nOrb

nOrb = size(ei)

allocate(U(nOrb,nOrb))

call cmplxUnitary(U, matrixToProcess, ei, tol)

matrixToProcess(:,:) = matmul(transpose(conjg(U)), matmul(matrixToProcess, U))

end subroutine


!> Transform degenerate states into combination that is orthogonal under the action of the matrix
subroutine cmplxUnitary(U, matrixToProcess, ei, tol)

!> Unitary matrix to transform orbitals
complex(dp), intent(out) :: U(:,:)

!> Matrix elements between (potentially degenerate) orbitals
complex(dp), intent(in) :: matrixToProcess(:,:)

!> Eigenvalues of states
real(dp), intent(in) :: ei(:)

!> Tolerance for comparisions between states
real(dp), intent(in), optional :: tol

integer :: ii, nOrb, nGrp, iGrp, maxRange, iRange
integer, allocatable :: degeneracies(:), degenerateRange(:,:)
complex(dp), allocatable :: subBlock(:,:)
real(dp), allocatable :: eigenvals(:)

U(:,:) = (0.0_dp,0.0_dp)
nOrb = size(ei)

do ii = 1, nOrb
  U(ii,ii) = (1.0_dp,0.0_dp)
end do


allocate(degenerateRange(2,nOrb))
call degeneracyRanges(degenerateRange, nGrp, Ei, tol)

maxRange = maxval(degenerateRange(2,:) - degenerateRange(1,:)) + 1

if (maxRange == 1) then
return
end if

allocate(subBlock(maxRange,maxRange))
allocate(eigenvals(maxRange))

do iGrp = 1, nGrp
subBlock(:,:) = (0.0_dp,0.0_dp)
eigenvals(:) = 0.0_dp
iRange = degenerateRange(2, iGrp) - degenerateRange(1, iGrp) + 1
subBlock(:iRange, :iRange) = matrixToProcess(degenerateRange(1, iGrp):degenerateRange(2, iGrp),&
& degenerateRange(1, iGrp):degenerateRange(2, iGrp))

call heev(subBlock(:iRange, :iRange), eigenvals, 'L', jobs='V')

U(degenerateRange(1, iGrp):degenerateRange(2, iGrp),&
& degenerateRange(1, iGrp):degenerateRange(2, iGrp)) = subBlock(:iRange, :iRange)

end do

end subroutine


!> Find which eigenvales are degenerate within tol
pure subroutine degeneracyRanges(degenerateRange, iGrp, Ei, tol)

integer, intent(out) :: degenerateRange(:,:)

integer, intent(out) :: iGrp

real(dp), intent(in) :: Ei(:)

real(dp), intent(in), optional :: tol


integer :: ii, jj, nOrb
real(dp) :: localTol

if (present(tol)) then
localTol = tol
else
localTol = epsilon(0.0_dp)
end if

nOrb = size(ei)

degenerateRange(:,:) = 0
iGrp = 0
ii = 1
do while (ii <= nOrb)
iGrp = iGrp + 1
degenerateRange(1, iGrp) = ii
do jj = ii + 1, nOrb

if ( ei(jj) - ei(jj-1) > localTol) then
  exit
end if
end do
ii = jj
degenerateRange(2, iGrp) = jj -1

end do

end subroutine

end module dftbp_rotateorbs

