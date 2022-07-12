!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

#! suffix and kinds for sparse with dense types
#:set ROUTINE_KINDS = [('real', 'r'), ('complex', 'c')]

!> DFTB+ sparse structure operations for a few blas matrix operations
module dftbp_math_sparseblas
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : imag, pi
  use dftbp_dftb_boundarycond, only : TBoundaryConditions
  use dftbp_type_commontypes, only : TOrbitals
#:if WITH_SCALAPACK
  use dftbp_common_blacsenv, only : TBlacsEnv
  use dftbp_extlibs_scalapackfx, only : blacsfx_gemr2d
  use dftbp_type_densedescr, only : TDenseDescr
#:endif
  implicit none
  private

  public :: symv, symm
#:if WITH_SCALAPACK
  public :: redist_sqr2rows, redist_rows2sqr

  interface redist_sqr2rows
    module procedure sqr2rows_real
    module procedure sqr2rows_complex
  end interface redist_sqr2rows

  interface redist_rows2sqr
    module procedure rows2sqr_real
    module procedure rows2sqr_complex
  end interface redist_rows2sqr

#:endif

  ! Rank 2 routines

  !> Routine for symmetric sparse matrix with dense vector multiply
  interface symv
  #:for _, suffix in ROUTINE_KINDS
    module procedure symv_${suffix}$_gamma
    module procedure symv_bc_${suffix}$_gamma
  #:endfor
    module procedure symv_kpt
    module procedure symv_bc_kpt
  end interface symv

  ! Rank 3 routines

  !> Routine for multiplication between a symmetric sparse matrix and a general
  !> dense matrix
  interface symm
  #:for _, suffix in ROUTINE_KINDS
    module procedure symm_${suffix}$_gamma
    module procedure symm_bc_${suffix}$_gamma
  #:endfor
    module procedure symm_kpt
    module procedure symm_bc_kpt
  end interface symm

contains

#:for typename, suffix in ROUTINE_KINDS

  !> Sparse Gamma point matrix with a ${typename}$ vector as a symv operation,
  !> y = alpha A x + beta * y
  subroutine symv_${suffix}$_gamma(y, A, x, iNeighbour, nNeighbour, img2CentCell, iSparseStart,&
      & iAtomStart, orb, alpha, beta)

    !> Vector on return
    ${typename}$(dp), intent(inout):: y(:)

    !> Sparse matrix
    real(dp), intent(in) :: A(:)

    !> Vector on entry
    ${typename}$(dp), intent(in) :: x(:)

    !> Atom neighbour list
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for atoms
    integer, intent(in) :: nNeighbour(:)

    !> Image to central cell mapping
    integer, intent(in) :: img2CentCell(:)

    !> Sparse indexing
    integer, intent(in) :: iSparseStart(0:,:)

    !> Dense indexing
    integer, intent(in) :: iAtomStart(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Scaling factor for A * x
    ${typename}$(dp), optional, intent(in) :: alpha

    !> Scaling factor for incoming y
    ${typename}$(dp), optional, intent(in) :: beta

    integer :: iOrig, ix, iy, jx, jy, iNeigh, nAtom, iAtom1, iAtom2, iAtom2f, nOrb1, nOrb2, ii
    real(dp) :: sqrTmp(orb%mOrb, orb%mOrb)
    ${typename}$(dp) :: alphaTmp

    @:ASSERT(size(x) == size(y))

    nAtom = size(iAtomStart)-1

    if (present(alpha)) then
      alphaTmp = alpha
    else
      alphaTmp = 1.0_dp
    end if

    if (present(beta)) then
      y(:) = beta * y
    else
      y(:) = 0.0_dp
    end if

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP &PRIVATE(iAtom1,nOrb1,ix,jx,iNeigh,sqrTmp,iAtom2,iAtom2f,nOrb2,iOrig,ii,iy,jy) &
    !$OMP & REDUCTION(+:y) SCHEDULE(RUNTIME)
    do iAtom1 = 1, nAtom
      ix = iAtomStart(iAtom1)
      jx = iAtomStart(iAtom1 + 1) -1
      nOrb1 = jx - ix + 1
      do iNeigh = 0, nNeighbour(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        sqrTmp(:,:) = 0.0_dp
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        iy = iAtomStart(iAtom2f)
        jy = iAtomStart(iAtom2f + 1) - 1
        nOrb2 = jy - iy + 1
        sqrTmp(:nOrb2,:nOrb1) = reshape(A(iOrig:iOrig+nOrb2*nOrb1-1), [nOrb2,nOrb1])
        ! symmetrize on-diagonal blocks just in case
        if (iAtom1 == iAtom2f) then
          do ii = 1, nOrb1
            sqrTmp(ii,ii+1:nOrb2) = sqrTmp(ii+1:nOrb2,ii)
          end do
        end if
        y(iy:jy) = y(iy:jy) + alphaTmp * matmul(sqrTmp(:nOrb2,:nOrb1),x(ix:jx))
        ! other triangle due to symmetry of matrix
        if (iAtom1 /= iAtom2f) then
          y(ix:jx) = y(ix:jx) + alphaTmp * matmul(transpose(sqrTmp(:nOrb2,:nOrb1)),x(iy:jy))
        end if
      end do
    end do
    !$OMP  END PARALLEL DO

  end subroutine symv_${suffix}$_gamma


  !> Sparse Gamma point matrix with a ${typename}$ vector as a symv operation, using specified
  !> boundary conditions, y = alpha A x + beta * y
  subroutine symv_bc_${suffix}$_gamma(y, A, x, bcs, iNeighbour, nNeighbour, img2CentCell,&
      & iSparseStart, iAtomStart, coords, species, orb, alpha, beta)

    !> Vector on return
    ${typename}$(dp), intent(inout):: y(:)

    !> Sparse matrix
    real(dp), intent(in) :: A(:)

    !> Vector on entry
    ${typename}$(dp), intent(in) :: x(:)

    !> Boundary conditions on the system
    type(TBoundaryConditions), intent(in) :: bcs

    !> Atom neighbour list
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for atoms
    integer, intent(in) :: nNeighbour(:)

    !> Image to central cell mapping
    integer, intent(in) :: img2CentCell(:)

    !> Sparse indexing
    integer, intent(in) :: iSparseStart(0:,:)

    !> Dense indexing
    integer, intent(in) :: iAtomStart(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coords(:,:)

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Scaling factor for A * x
    ${typename}$(dp), optional, intent(in) :: alpha

    !> Scaling factor for incoming y
    ${typename}$(dp), optional, intent(in) :: beta

    integer :: iOrig, ix, iy, jx, jy, iNeigh, nAtom, iAtom1, iAtom2, iAtom2f, nOrb1, nOrb2, ii
    real(dp) :: sqrTmp(orb%mOrb, orb%mOrb)
    ${typename}$(dp) :: alphaTmp

    @:ASSERT(size(x) == size(y))

    nAtom = size(iAtomStart)-1

    if (present(alpha)) then
      alphaTmp = alpha
    else
      alphaTmp = 1.0_dp
    end if

    if (present(beta)) then
      y(:) = beta * y
    else
      y(:) = 0.0_dp
    end if

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP &PRIVATE(iAtom1,nOrb1,ix,jx,iNeigh,sqrTmp,iAtom2,iAtom2f,nOrb2,iOrig,ii,iy,jy) &
    !$OMP & REDUCTION(+:y) SCHEDULE(RUNTIME)
    do iAtom1 = 1, nAtom
      ix = iAtomStart(iAtom1)
      jx = iAtomStart(iAtom1 + 1) -1
      nOrb1 = jx - ix + 1
      do iNeigh = 0, nNeighbour(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        sqrTmp(:,:) = 0.0_dp
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        iy = iAtomStart(iAtom2f)
        jy = iAtomStart(iAtom2f + 1) - 1
        nOrb2 = jy - iy + 1
        sqrTmp(:nOrb2,:nOrb1) = reshape(A(iOrig:iOrig+nOrb2*nOrb1-1), [nOrb2,nOrb1])
        call bcs%foldInDiatomicBlock(sqrTmp, iAtom1, iAtom2, coords, species, img2centCell, orb)
        ! symmetrize on-diagonal blocks just in case
        if (iAtom1 == iAtom2f) then
          do ii = 1, nOrb1
            sqrTmp(ii,ii+1:nOrb2) = sqrTmp(ii+1:nOrb2,ii)
          end do
        end if
        y(iy:jy) = y(iy:jy) + alphaTmp * matmul(sqrTmp(:nOrb2,:nOrb1),x(ix:jx))
        ! other triangle due to symmetry of matrix
        if (iAtom1 /= iAtom2f) then
          y(ix:jx) = y(ix:jx) + alphaTmp * matmul(transpose(sqrTmp(:nOrb2,:nOrb1)),x(iy:jy))
        end if
      end do
    end do
    !$OMP  END PARALLEL DO

  end subroutine symv_bc_${suffix}$_gamma

#:endfor


  !> Sparse matrix at specified k-point with a vector as a symv operation, y = alpha A x + beta * y
  subroutine symv_kpt(y, A, x, iNeighbour, nNeighbour, kPoint, iCellVec, cellVec, img2CentCell,&
      & iSparseStart, iAtomStart, orb, alpha, beta)

    !> Vector on return
    complex(dp), intent(inout):: y(:)

    !> Sparse matrix
    real(dp), intent(in) :: A(:)

    !> Vector on entry
    complex(dp), intent(in) :: x(:)

    !> Atom neighbour list
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for atoms
    integer, intent(in) :: nNeighbour(:)

    !> Relative coordinates of the k-point where the sparse matrix should be unfolded.
    real(dp), intent(in) :: kPoint(:)

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> Image to central cell mapping
    integer, intent(in) :: img2CentCell(:)

    !> Sparse indexing
    integer, intent(in) :: iSparseStart(0:,:)

    !> Dense indexing
    integer, intent(in) :: iAtomStart(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Scaling factor for A * x
    complex(dp), optional, intent(in) :: alpha

    !> Scaling factor for incoming y
    complex(dp), optional, intent(in) :: beta

    integer :: iOrig, ix, iy, jx, jy, iNeigh, nAtom, iAtom1, iAtom2, iAtom2f, nOrb1, nOrb2, ii, iVec
    complex(dp) :: sqrTmp(orb%mOrb, orb%mOrb), alphaTmp, phase
    real(dp) :: kPoint2p(3)

    @:ASSERT(size(x) == size(y))

    nAtom = size(iAtomStart)-1

    kPoint2p(:) = 2.0_dp * pi * kPoint

    if (present(alpha)) then
      alphaTmp = alpha
    else
      alphaTmp = 1.0_dp
    end if

    if (present(beta)) then
      y(:) = beta * y
    else
      y(:) = 0.0_dp
    end if

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP &PRIVATE(iAtom1,nOrb1,ix,jx,iNeigh,sqrTmp,iAtom2,iAtom2f,nOrb2,iOrig,ii,iy,jy,iVec,phase)&
    !$OMP &REDUCTION(+:y) SCHEDULE(RUNTIME)
    do iAtom1 = 1, nAtom
      ix = iAtomStart(iAtom1)
      jx = iAtomStart(iAtom1 + 1) -1
      nOrb1 = jx - ix + 1
      do iNeigh = 0, nNeighbour(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        sqrTmp(:,:) = 0.0_dp
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        iy = iAtomStart(iAtom2f)
        jy = iAtomStart(iAtom2f + 1) - 1
        nOrb2 = jy - iy + 1
        iVec = iCellVec(iAtom2)
        phase = exp(imag * dot_product(kPoint2p, cellVec(:, iVec)))
        sqrTmp(:nOrb2,:nOrb1) = phase * reshape(A(iOrig:iOrig+nOrb2*nOrb1-1), [nOrb2,nOrb1])
        ! hermitian symmetry for on-diagonal blocks, just in case
        if (iAtom1 == iAtom2f) then
          do ii = 1, nOrb1
            sqrTmp(ii,ii+1:nOrb2) = conjg(sqrTmp(ii+1:nOrb2,ii))
          end do
        end if
        y(iy:jy) = y(iy:jy) + alphaTmp * matmul(sqrTmp(:nOrb2,:nOrb1),x(ix:jx))
        ! other triangle due to symmetry of matrix
        if (iAtom1 /= iAtom2f) then
          y(ix:jx) = y(ix:jx) + alphaTmp*matmul(transpose(conjg(sqrTmp(:nOrb2,:nOrb1))),x(iy:jy))
        end if
      end do
    end do
    !$OMP  END PARALLEL DO

  end subroutine symv_kpt


  !> Sparse matrix at specified k-point with a vector as a symv operation, using specified boundary
  !> conditions, y = alpha A x + beta * y
  subroutine symv_bc_kpt(y, A, x, bcs, iNeighbour, nNeighbour, kPoint, iCellVec, cellVec,&
      & img2CentCell, iSparseStart, iAtomStart, coords, species, orb, alpha, beta)

    !> Vector on return
    complex(dp), intent(inout):: y(:)

    !> Sparse matrix
    real(dp), intent(in) :: A(:)

    !> Vector on entry
    complex(dp), intent(in) :: x(:)

    !> Boundary conditions on the system
    type(TBoundaryConditions), intent(in) :: bcs

    !> Atom neighbour list
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for atoms
    integer, intent(in) :: nNeighbour(:)

    !> Relative coordinates of the k-point where the sparse matrix should be unfolded.
    real(dp), intent(in) :: kPoint(:)

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> Image to central cell mapping
    integer, intent(in) :: img2CentCell(:)

    !> Sparse indexing
    integer, intent(in) :: iSparseStart(0:,:)

    !> Dense indexing
    integer, intent(in) :: iAtomStart(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coords(:,:)

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Scaling factor for A * x
    complex(dp), optional, intent(in) :: alpha

    !> Scaling factor for incoming y
    complex(dp), optional, intent(in) :: beta

    integer :: iOrig, ix, iy, jx, jy, iNeigh, nAtom, iAtom1, iAtom2, iAtom2f, nOrb1, nOrb2, ii, iVec
    complex(dp) :: sqrTmp(orb%mOrb, orb%mOrb), alphaTmp, phase
    real(dp) :: kPoint2p(2)

    @:ASSERT(size(x) == size(y))

    nAtom = size(iAtomStart)-1

    kPoint2p(:) = 2.0_dp * pi * kPoint

    if (present(alpha)) then
      alphaTmp = alpha
    else
      alphaTmp = 1.0_dp
    end if

    if (present(beta)) then
      y(:) = beta * y
    else
      y(:) = 0.0_dp
    end if

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP &PRIVATE(iAtom1,nOrb1,ix,jx,iNeigh,sqrTmp,iAtom2,iAtom2f,nOrb2,iOrig,ii,iy,jy,iVec,phase)&
    !$OMP &REDUCTION(+:y) SCHEDULE(RUNTIME)
    do iAtom1 = 1, nAtom
      ix = iAtomStart(iAtom1)
      jx = iAtomStart(iAtom1 + 1) -1
      nOrb1 = jx - ix + 1
      do iNeigh = 0, nNeighbour(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        sqrTmp(:,:) = 0.0_dp
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        iy = iAtomStart(iAtom2f)
        jy = iAtomStart(iAtom2f + 1) - 1
        nOrb2 = jy - iy + 1
        iVec = iCellVec(iAtom2)
        phase = exp(imag * dot_product(kPoint2p, cellVec(:, iVec)))
        sqrTmp(:nOrb2,:nOrb1) = reshape(A(iOrig:iOrig+nOrb2*nOrb1-1), [nOrb2,nOrb1])
        call bcs%foldInDiatomicBlock(sqrTmp, iAtom1, iAtom2, coords, species, img2centCell, orb)
        sqrTmp(:nOrb2,:nOrb1) = phase * sqrTmp
        ! hermitian symmetry for on-diagonal blocks, just in case
        if (iAtom1 == iAtom2f) then
          do ii = 1, nOrb1
            sqrTmp(ii,ii+1:nOrb2) = conjg(sqrTmp(ii+1:nOrb2,ii))
          end do
        end if
        y(iy:jy) = y(iy:jy) + alphaTmp * matmul(sqrTmp(:nOrb2,:nOrb1),x(ix:jx))
        ! other triangle due to symmetry of matrix
        if (iAtom1 /= iAtom2f) then
          y(ix:jx) = y(ix:jx) + alphaTmp*matmul(transpose(conjg(sqrTmp(:nOrb2,:nOrb1))),x(iy:jy))
        end if
      end do
    end do
    !$OMP  END PARALLEL DO

  end subroutine symv_bc_kpt


#:for typename, suffix in ROUTINE_KINDS

  !> Sparse Gamma point matrix with dense ${typename}$ matrix as a symm operation,
  !> symmetric matrix on 'l'eft or 'r'ight , where A is symmetric and sparse and B is general.
  subroutine symm_${suffix}$_gamma(C, side, A, B, iNeighbour, nNeighbour, img2CentCell,&
      & iSparseStart, iAtomStart, orb, alpha, beta)

    !> Dense matrix on return
    ${typename}$(dp), intent(inout):: C(:,:)

    !> Sparse matrix
    real(dp), intent(in) :: A(:)

    !> Side which is the sparse matrix. SIDE = 'L' or 'l' C := alpha*A*B + beta*C; SIDE = 'R' or 'r'
    !> C := alpha*B*A + beta*C
    character, intent(in) :: side

    !> Dense matrix on entry
    ${typename}$(dp), intent(in) :: B(:,:)

    !> Atom neighbour list
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for atoms
    integer, intent(in) :: nNeighbour(:)

    !> Image to central cell mapping
    integer, intent(in) :: img2CentCell(:)

    !> Sparse indexing
    integer, intent(in) :: iSparseStart(0:,:)

    !> Dense indexing
    integer, intent(in) :: iAtomStart(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Scaling factor for A * B
    ${typename}$(dp), optional, intent(in) :: alpha

    !> Scaling factor for incoming C
    ${typename}$(dp), optional, intent(in) :: beta

    integer :: iOrig, iB, iC, jB, jC, iNeigh, nAtom, iAtom1, iAtom2, iAtom2f, nOrb1, nOrb2, ii
    real(dp) :: sqrTmp(orb%mOrb, orb%mOrb)
    ${typename}$(dp) :: alphaTmp

    @:ASSERT(all(shape(C) == shape(B)))
    @:ASSERT(side == 'r' .or. side == 'R' .or. side == 'l' .or. side == 'L')

    nAtom = size(iAtomStart)-1

    if (present(alpha)) then
      alphaTmp = alpha
    else
      alphaTmp = 1.0_dp
    end if

    if (present(beta)) then
      C(:,:) = beta * C
    else
      C(:,:) = 0.0_dp
    end if

    select case(side)

    case ("L", "l")

      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP &PRIVATE(iAtom1,nOrb1,iB,jB,iNeigh,sqrTmp,iAtom2,iAtom2f,nOrb2,iOrig,ii,iC,jC) &
      !$OMP & REDUCTION(+:C) SCHEDULE(RUNTIME)
      do iAtom1 = 1, nAtom
        iB = iAtomStart(iAtom1)
        jB = iAtomStart(iAtom1 + 1) - 1
        nOrb1 = jB - iB + 1
        do iNeigh = 0, nNeighbour(iAtom1)
          iOrig = iSparseStart(iNeigh,iAtom1) + 1
          sqrTmp(:,:) = 0.0_dp
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          iC = iAtomStart(iAtom2f)
          jC = iAtomStart(iAtom2f + 1) - 1
          nOrb2 = jC - iC + 1
          sqrTmp(:nOrb2,:nOrb1) = reshape(A(iOrig:iOrig+nOrb2*nOrb1-1), [nOrb2,nOrb1])
          ! symmetrize on-diagonal blocks just in case
          if (iAtom1 == iAtom2f) then
            do ii = 1, nOrb1
              sqrTmp(ii,ii+1:nOrb2) = sqrTmp(ii+1:nOrb2,ii)
            end do
          end if
          C(iC:jC,:) = C(iC:jC,:) + alphaTmp*matmul(sqrTmp(:nOrb2,:nOrb1),B(iB:jB,:))
          ! other triangle due to symmetry of matrix
          if (iAtom1 /= iAtom2f) then
            C(iB:jB, :) = C(iB:jB,:) + alphaTmp*matmul(transpose(sqrTmp(:nOrb2,:nOrb1)),B(iC:jC,:))
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    case ("R", "r")

      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP &PRIVATE(iAtom1,nOrb1,iB,jB,iNeigh,sqrTmp,iAtom2,iAtom2f,nOrb2,iOrig,ii,iC,jC) &
      !$OMP & REDUCTION(+:C) SCHEDULE(RUNTIME)
      do iAtom1 = 1, nAtom
        iB = iAtomStart(iAtom1)
        jB = iAtomStart(iAtom1 + 1) - 1
        nOrb1 = jB - iB + 1
        do iNeigh = 0, nNeighbour(iAtom1)
          iOrig = iSparseStart(iNeigh,iAtom1) + 1
          sqrTmp(:,:) = 0.0_dp
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          iC = iAtomStart(iAtom2f)
          jC = iAtomStart(iAtom2f + 1) - 1
          nOrb2 = jC - iC + 1
          sqrTmp(:nOrb2,:nOrb1) = reshape(A(iOrig:iOrig+nOrb2*nOrb1-1), [nOrb2,nOrb1])
          ! symmetrize on-diagonal blocks just in case
          if (iAtom1 == iAtom2f) then
            do ii = 1, nOrb1
              sqrTmp(ii,ii+1:nOrb2) = sqrTmp(ii+1:nOrb2,ii)
            end do
          end if
          C(:,iC:jC) = C(:,iC:jC) + alphaTmp*matmul(B(:,iB:jB),transpose(sqrTmp(:nOrb2,:nOrb1)))
          ! other triangle due to symmetry of matrix
          if (iAtom1 /= iAtom2f) then
            C(:,iB:jB) = C(:,iB:jB) + alphaTmp*matmul(B(:,iC:jC),sqrTmp(:nOrb2,:nOrb1))
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    end select

  end subroutine symm_${suffix}$_gamma


  !> Sparse Gamma point matrix with dense ${typename}$ matrix as a symm operation, using specified
  !> boundary conditions. Symmetric matrix on 'l'eft or 'r'ight , where A is symmetric and sparse
  !> and B is general.
  subroutine symm_bc_${suffix}$_gamma(C, side, A, B, bcs, iNeighbour, nNeighbour, img2CentCell,&
      & iSparseStart, iAtomStart, coords, species, orb, alpha, beta)

    !> Dense matrix on return
    ${typename}$(dp), intent(inout):: C(:,:)

    !> Sparse matrix
    real(dp), intent(in) :: A(:)

    !> Side which is the sparse matrix. SIDE = 'L' or 'l' C := alpha*A*B + beta*C; SIDE = 'R' or 'r'
    !> C := alpha*B*A + beta*C
    character, intent(in) :: side

    !> Dense matrix on entry
    ${typename}$(dp), intent(in) :: B(:,:)

    !> Boundary conditions on the system
    type(TBoundaryConditions), intent(in) :: bcs

    !> Atom neighbour list
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for atoms
    integer, intent(in) :: nNeighbour(:)

    !> Image to central cell mapping
    integer, intent(in) :: img2CentCell(:)

    !> Sparse indexing
    integer, intent(in) :: iSparseStart(0:,:)

    !> Dense indexing
    integer, intent(in) :: iAtomStart(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coords(:,:)

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Scaling factor for A * B
    ${typename}$(dp), optional, intent(in) :: alpha

    !> Scaling factor for incoming C
    ${typename}$(dp), optional, intent(in) :: beta

    integer :: iOrig, iB, iC, jB, jC, iNeigh, nAtom, iAtom1, iAtom2, iAtom2f, nOrb1, nOrb2, ii
    real(dp) :: sqrTmp(orb%mOrb, orb%mOrb)
    ${typename}$(dp) :: alphaTmp

    @:ASSERT(all(shape(C) == shape(B)))
    @:ASSERT(side == 'r' .or. side == 'R' .or. side == 'l' .or. side == 'L')

    nAtom = size(iAtomStart)-1

    if (present(alpha)) then
      alphaTmp = alpha
    else
      alphaTmp = 1.0_dp
    end if

    if (present(beta)) then
      C(:,:) = beta * C
    else
      C(:,:) = 0.0_dp
    end if

    select case(side)

    case ("L", "l")

      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP &PRIVATE(iAtom1,nOrb1,iB,jB,iNeigh,sqrTmp,iAtom2,iAtom2f,nOrb2,iOrig,ii,iC,jC) &
      !$OMP & REDUCTION(+:C) SCHEDULE(RUNTIME)
      do iAtom1 = 1, nAtom
        iB = iAtomStart(iAtom1)
        jB = iAtomStart(iAtom1 + 1) - 1
        nOrb1 = jB - iB + 1
        do iNeigh = 0, nNeighbour(iAtom1)
          iOrig = iSparseStart(iNeigh,iAtom1) + 1
          sqrTmp(:,:) = 0.0_dp
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          iC = iAtomStart(iAtom2f)
          jC = iAtomStart(iAtom2f + 1) - 1
          nOrb2 = jC - iC + 1
          sqrTmp(:nOrb2,:nOrb1) = reshape(A(iOrig:iOrig+nOrb2*nOrb1-1), [nOrb2,nOrb1])
          call bcs%foldInDiatomicBlock(sqrTmp, iAtom1, iAtom2, coords, species, img2centCell, orb)
          ! symmetrize on-diagonal blocks just in case
          if (iAtom1 == iAtom2f) then
            do ii = 1, nOrb1
              sqrTmp(ii,ii+1:nOrb2) = sqrTmp(ii+1:nOrb2,ii)
            end do
          end if
          C(iC:jC,:) = C(iC:jC,:) + alphaTmp*matmul(sqrTmp(:nOrb2,:nOrb1),B(iB:jB,:))
          ! other triangle due to symmetry of matrix
          if (iAtom1 /= iAtom2f) then
            C(iB:jB, :) = C(iB:jB,:) + alphaTmp*matmul(transpose(sqrTmp(:nOrb2,:nOrb1)),B(iC:jC,:))
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    case ("R", "r")

      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP &PRIVATE(iAtom1,nOrb1,iB,jB,iNeigh,sqrTmp,iAtom2,iAtom2f,nOrb2,iOrig,ii,iC,jC) &
      !$OMP & REDUCTION(+:C) SCHEDULE(RUNTIME)
      do iAtom1 = 1, nAtom
        iB = iAtomStart(iAtom1)
        jB = iAtomStart(iAtom1 + 1) - 1
        nOrb1 = jB - iB + 1
        do iNeigh = 0, nNeighbour(iAtom1)
          iOrig = iSparseStart(iNeigh,iAtom1) + 1
          sqrTmp(:,:) = 0.0_dp
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          iC = iAtomStart(iAtom2f)
          jC = iAtomStart(iAtom2f + 1) - 1
          nOrb2 = jC - iC + 1
          sqrTmp(:nOrb2,:nOrb1) = reshape(A(iOrig:iOrig+nOrb2*nOrb1-1), [nOrb2,nOrb1])
          call bcs%foldInDiatomicBlock(sqrTmp, iAtom1, iAtom2, coords, species, img2centCell, orb)
          ! symmetrize on-diagonal blocks just in case
          if (iAtom1 == iAtom2f) then
            do ii = 1, nOrb1
              sqrTmp(ii,ii+1:nOrb2) = sqrTmp(ii+1:nOrb2,ii)
            end do
          end if
          C(:,iC:jC) = C(:,iC:jC) + alphaTmp*matmul(B(:,iB:jB),transpose(sqrTmp(:nOrb2,:nOrb1)))
          ! other triangle due to symmetry of matrix
          if (iAtom1 /= iAtom2f) then
            C(:,iB:jB) = C(:,iB:jB) + alphaTmp*matmul(B(:,iC:jC),sqrTmp(:nOrb2,:nOrb1))
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    end select

  end subroutine symm_bc_${suffix}$_gamma

#:endfor


  !> Sparse matrix at specified k-point with dense ${typename}$ matrix as a symm operation,
  !> symmetric matrix on 'l'eft or 'r'ight , where A is symmetric and sparse and B is general.
  subroutine symm_kpt(C, side, A, B, iNeighbour, nNeighbour, kPoint, iCellVec, cellVec,&
      & img2CentCell, iSparseStart, iAtomStart, orb, alpha, beta)

    !> Dense matrix on return
    complex(dp), intent(inout):: C(:,:)

    !> Sparse matrix
    real(dp), intent(in) :: A(:)

    !> Side which is the sparse matrix. SIDE = 'L' or 'l' C := alpha*A*B + beta*C; SIDE = 'R' or 'r'
    !> C := alpha*B*A + beta*C
    character, intent(in) :: side

    !> Dense matrix on entry
    complex(dp), intent(in) :: B(:,:)

    !> Atom neighbour list
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for atoms
    integer, intent(in) :: nNeighbour(:)

    !> Relative coordinates of the k-point where the sparse matrix should be unfolded.
    real(dp), intent(in) :: kPoint(:)

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> Image to central cell mapping
    integer, intent(in) :: img2CentCell(:)

    !> Sparse indexing
    integer, intent(in) :: iSparseStart(0:,:)

    !> Dense indexing
    integer, intent(in) :: iAtomStart(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Scaling factor for A * B
    complex(dp), optional, intent(in) :: alpha

    !> Scaling factor for incoming C
    complex(dp), optional, intent(in) :: beta

    integer :: iOrig, iB, iC, jB, jC, iNeigh, nAtom, iAtom1, iAtom2, iAtom2f, nOrb1, nOrb2, ii, iVec
    complex(dp) :: sqrTmp(orb%mOrb, orb%mOrb), alphaTmp, phase
    real(dp) :: kPoint2p(3)

    @:ASSERT(all(shape(C) == shape(B)))
    @:ASSERT(side == 'r' .or. side == 'R' .or. side == 'l' .or. side == 'L')

    nAtom = size(iAtomStart)-1

    kPoint2p(:) = 2.0_dp * pi * kPoint

    if (present(alpha)) then
      alphaTmp = alpha
    else
      alphaTmp = 1.0_dp
    end if

    if (present(beta)) then
      C(:,:) = beta * C
    else
      C(:,:) = 0.0_dp
    end if

    select case(side)

    case ("L", "l")

      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP &PRIVATE(nOrb1,iB,jB,iNeigh,sqrTmp,iAtom2,iAtom2f,nOrb2,iOrig,ii,iC,jC,iVec,&
      !$OMP &phase) REDUCTION(+:C) SCHEDULE(RUNTIME)
      do iAtom1 = 1, nAtom
        iB = iAtomStart(iAtom1)
        jB = iAtomStart(iAtom1 + 1) - 1
        nOrb1 = jB - iB + 1
        do iNeigh = 0, nNeighbour(iAtom1)
          iOrig = iSparseStart(iNeigh,iAtom1) + 1
          sqrTmp(:,:) = 0.0_dp
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          iC = iAtomStart(iAtom2f)
          jC = iAtomStart(iAtom2f + 1) - 1
          nOrb2 = jC - iC + 1
          iVec = iCellVec(iAtom2)
          phase = exp(imag * dot_product(kPoint2p, cellVec(:, iVec)))
          sqrTmp(:nOrb2,:nOrb1) = phase * reshape(A(iOrig:iOrig+nOrb2*nOrb1-1), [nOrb2,nOrb1])
          ! Hermitian symmetry on-diagonal blocks just in case
          if (iAtom1 == iAtom2f) then
            do ii = 1, nOrb1
              sqrTmp(ii,ii+1:nOrb2) = conjg(sqrTmp(ii+1:nOrb2,ii))
            end do
          end if
          C(iC:jC,:) = C(iC:jC,:) + alphaTmp*matmul(sqrTmp(:nOrb2,:nOrb1),B(iB:jB,:))
          ! other triangle due to symmetry of matrix
          if (iAtom1 /= iAtom2f) then
            C(iB:jB, :) = C(iB:jB,:)&
                & + alphaTmp * matmul(transpose(conjg(sqrTmp(:nOrb2,:nOrb1))),B(iC:jC,:))
          end if

        end do
      end do
      !$OMP END PARALLEL DO

    case ("R", "r")

      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP &PRIVATE(nOrb1,iB,jB,iNeigh,sqrTmp,iAtom2,iAtom2f,nOrb2,iOrig,ii,iC,jC,iVec,&
      !$OMP &phase) REDUCTION(+:C) SCHEDULE(RUNTIME)
      do iAtom1 = 1, nAtom
        iB = iAtomStart(iAtom1)
        jB = iAtomStart(iAtom1 + 1) - 1
        nOrb1 = jB - iB + 1
        do iNeigh = 0, nNeighbour(iAtom1)
          iOrig = iSparseStart(iNeigh,iAtom1) + 1
          sqrTmp(:,:) = 0.0_dp
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          iC = iAtomStart(iAtom2f)
          jC = iAtomStart(iAtom2f + 1) - 1
          nOrb2 = jC - iC + 1
          iVec = iCellVec(iAtom2)
          phase = exp(imag * dot_product(kPoint2p, cellVec(:, iVec)))
          sqrTmp(:nOrb2,:nOrb1) = phase * reshape(A(iOrig:iOrig+nOrb2*nOrb1-1), [nOrb2,nOrb1])
          ! Hermitian symmetry on-diagonal blocks just in case
          if (iAtom1 == iAtom2f) then
            do ii = 1, nOrb1
              sqrTmp(ii,ii+1:nOrb2) = conjg(sqrTmp(ii+1:nOrb2,ii))
            end do
          end if
          C(:,iC:jC) = C(:,iC:jC)&
              & + alphaTmp * matmul(B(:,iB:jB),transpose(conjg(sqrTmp(:nOrb2,:nOrb1))))
          ! other triangle due to symmetry of matrix
          if (iAtom1 /= iAtom2f) then
            C(:,iB:jB) = C(:,iB:jB) + alphaTmp*matmul(B(:,iC:jC),sqrTmp(:nOrb2,:nOrb1))
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    end select

  end subroutine symm_kpt


  !> Sparse matrix at specified k-point with dense ${typename}$ matrix as a symm operation, using
  !> specified boundary conditions. symmetric matrix on 'l'eft or 'r'ight , where A is symmetric and
  !> sparse and B is general.
  subroutine symm_bc_kpt(C, side, A, B, bcs, iNeighbour, nNeighbour, kPoint, iCellVec, cellVec,&
      & img2CentCell, iSparseStart, iAtomStart, coords, species, orb, alpha, beta)

    !> Dense matrix on return
    complex(dp), intent(inout):: C(:,:)

    !> Sparse matrix
    real(dp), intent(in) :: A(:)

    !> Side which is the sparse matrix. SIDE = 'L' or 'l' C := alpha*A*B + beta*C; SIDE = 'R' or 'r'
    !> C := alpha*B*A + beta*C
    character, intent(in) :: side

    !> Dense matrix on entry
    complex(dp), intent(in) :: B(:,:)

    !> Boundary conditions on the system
    type(TBoundaryConditions), intent(in) :: bcs

    !> Atom neighbour list
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for atoms
    integer, intent(in) :: nNeighbour(:)

    !> Relative coordinates of the k-point where the sparse matrix should be unfolded.
    real(dp), intent(in) :: kPoint(:)

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> Image to central cell mapping
    integer, intent(in) :: img2CentCell(:)

    !> Sparse indexing
    integer, intent(in) :: iSparseStart(0:,:)

    !> Dense indexing
    integer, intent(in) :: iAtomStart(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coords(:,:)

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Scaling factor for A * B
    complex(dp), optional, intent(in) :: alpha

    !> Scaling factor for incoming C
    complex(dp), optional, intent(in) :: beta

    integer :: iOrig, iB, iC, jB, jC, iNeigh, nAtom, iAtom1, iAtom2, iAtom2f, nOrb1, nOrb2, ii, iVec
    complex(dp) :: sqrTmp(orb%mOrb, orb%mOrb), alphaTmp, phase
    real(dp) :: kPoint2p(2)

    @:ASSERT(all(shape(C) == shape(B)))
    @:ASSERT(side == 'r' .or. side == 'R' .or. side == 'l' .or. side == 'L')

    nAtom = size(iAtomStart)-1

    kPoint2p(:) = 2.0_dp * pi * kPoint

    if (present(alpha)) then
      alphaTmp = alpha
    else
      alphaTmp = 1.0_dp
    end if

    if (present(beta)) then
      C(:,:) = beta * C
    else
      C(:,:) = 0.0_dp
    end if

    select case(side)

    case ("L", "l")

      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP &PRIVATE(nOrb1,iB,jB,iNeigh,sqrTmp,iAtom2,iAtom2f,nOrb2,iOrig,ii,iC,jC,iVec,&
      !$OMP &phase) REDUCTION(+:C) SCHEDULE(RUNTIME)
      do iAtom1 = 1, nAtom
        iB = iAtomStart(iAtom1)
        jB = iAtomStart(iAtom1 + 1) - 1
        nOrb1 = jB - iB + 1
        do iNeigh = 0, nNeighbour(iAtom1)
          iOrig = iSparseStart(iNeigh,iAtom1) + 1
          sqrTmp(:,:) = 0.0_dp
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          iC = iAtomStart(iAtom2f)
          jC = iAtomStart(iAtom2f + 1) - 1
          nOrb2 = jC - iC + 1
          iVec = iCellVec(iAtom2)
          phase = exp(imag * dot_product(kPoint2p, cellVec(:, iVec)))
          sqrTmp(:nOrb2,:nOrb1) = reshape(A(iOrig:iOrig+nOrb2*nOrb1-1), [nOrb2,nOrb1])
          call bcs%foldInDiatomicBlock(sqrTmp, iAtom1, iAtom2, coords, species, img2centCell, orb)
          sqrTmp(:nOrb2,:nOrb1) = phase * sqrTmp
          ! Hermitian symmetry on-diagonal blocks just in case
          if (iAtom1 == iAtom2f) then
            do ii = 1, nOrb1
              sqrTmp(ii,ii+1:nOrb2) = conjg(sqrTmp(ii+1:nOrb2,ii))
            end do
          end if
          C(iC:jC,:) = C(iC:jC,:) + alphaTmp*matmul(sqrTmp(:nOrb2,:nOrb1),B(iB:jB,:))
          ! other triangle due to symmetry of matrix
          if (iAtom1 /= iAtom2f) then
            C(iB:jB, :) = C(iB:jB,:)&
                & + alphaTmp * matmul(transpose(conjg(sqrTmp(:nOrb2,:nOrb1))),B(iC:jC,:))
          end if

        end do
      end do
      !$OMP END PARALLEL DO

    case ("R", "r")

      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP &PRIVATE(nOrb1,iB,jB,iNeigh,sqrTmp,iAtom2,iAtom2f,nOrb2,iOrig,ii,iC,jC,iVec,&
      !$OMP &phase) REDUCTION(+:C) SCHEDULE(RUNTIME)
      do iAtom1 = 1, nAtom
        iB = iAtomStart(iAtom1)
        jB = iAtomStart(iAtom1 + 1) - 1
        nOrb1 = jB - iB + 1
        do iNeigh = 0, nNeighbour(iAtom1)
          iOrig = iSparseStart(iNeigh,iAtom1) + 1
          sqrTmp(:,:) = 0.0_dp
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          iC = iAtomStart(iAtom2f)
          jC = iAtomStart(iAtom2f + 1) - 1
          nOrb2 = jC - iC + 1
          iVec = iCellVec(iAtom2)
          phase = exp(imag * dot_product(kPoint2p, cellVec(:, iVec)))
          sqrTmp(:nOrb2,:nOrb1) = reshape(A(iOrig:iOrig+nOrb2*nOrb1-1), [nOrb2,nOrb1])
          sqrTmp(:nOrb2,:nOrb1) = phase * sqrTmp
          call bcs%foldInDiatomicBlock(sqrTmp, iAtom1, iAtom2, coords, species, img2centCell, orb)
          ! Hermitian symmetry on-diagonal blocks just in case
          if (iAtom1 == iAtom2f) then
            do ii = 1, nOrb1
              sqrTmp(ii,ii+1:nOrb2) = conjg(sqrTmp(ii+1:nOrb2,ii))
            end do
          end if
          C(:,iC:jC) = C(:,iC:jC)&
              & + alphaTmp * matmul(B(:,iB:jB),transpose(conjg(sqrTmp(:nOrb2,:nOrb1))))
          ! other triangle due to symmetry of matrix
          if (iAtom1 /= iAtom2f) then
            C(:,iB:jB) = C(:,iB:jB) + alphaTmp*matmul(B(:,iC:jC),sqrTmp(:nOrb2,:nOrb1))
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    end select

  end subroutine symm_bc_kpt

#:if WITH_SCALAPACK

#:for VAR in [('real'),('complex')]

  !> Re-distributes data for square matrices between BLACS block cyclic data and whole global rows
  !> on each processor
  subroutine sqr2rows_${VAR}$(square, row, denseDesc, blacsEnv)

    !> Real matrix in block cyclic, last index over spin/kpts
    ${VAR}$(dp), allocatable, intent(in) :: square(:,:,:)

    !> Real matrix with individual rows on each processor, last index over spin/kpts
    ${VAR}$(dp), allocatable, intent(inout) :: row(:,:,:)

    !> Descriptors for dense matrices
    type(TDenseDescr), intent(in) :: denseDesc

    !> BLACS environment and information on grids
    type(TBlacsEnv), intent(in) :: blacsEnv

    integer :: iKS

    @:ASSERT(allocated(square) .eqv. allocated(row))

    if (.not.allocated(square)) then
      return
    end if

    #! same group choices over spin/k for the two grids, so same size on last index
    @:ASSERT(size(square, dim=3) == size(row, dim=3))

    #! Ensure these are over the same contexts
    @:ASSERT(blacsEnv%orbitalGrid%ctxt == blacsEnv%rowOrbitalGrid%ctxt)

    do iKS = 1, size(square, dim=3)
      call blacsfx_gemr2d(blacsEnv%nn, blacsEnv%nn, square(:,:,iKS), 1, 1,&
          & denseDesc%blacsOrbSqr, row(:,:,iKS), 1, 1, denseDesc%blacsColumnSqr,&
          & blacsEnv%orbitalGrid%ctxt)
    end do

  end subroutine sqr2rows_${VAR}$


  !> Re-distributes data for square matrices between BLACS whole global rows on each processor and
  !> block cyclic data
  subroutine rows2sqr_${VAR}$(row, square, denseDesc, blacsEnv)

    !> Real matrix with individual rows on each processor, last index over spin/kpts
    ${VAR}$(dp), allocatable, intent(in) :: row(:,:,:)

    !> Real matrix in block cyclic, last index over spin/kpts
    ${VAR}$(dp), allocatable, intent(inout) :: square(:,:,:)

    !> Descriptors for dense matrices
    type(TDenseDescr), intent(in) :: denseDesc

    !> BLACS environment and information on grids
    type(TBlacsEnv), intent(in) :: blacsEnv

    integer :: iKS

    @:ASSERT(allocated(row) .eqv. allocated(square))

    if (.not.allocated(row)) then
      return
    end if

    #! same group choices over spin/k for the two grids, so same size on last index
    @:ASSERT(size(row, dim=3) == size(square, dim=3))

    #! Ensure these are over the same contexts
    @:ASSERT(blacsEnv%orbitalGrid%ctxt == blacsEnv%rowOrbitalGrid%ctxt)

    do iKS = 1, size(square, dim=3)
      call blacsfx_gemr2d(blacsEnv%nn, blacsEnv%nn, row(:,:,iKS), 1, 1,&
          & denseDesc%blacsColumnSqr, square(:,:,iKS), 1, 1, denseDesc%blacsOrbSqr,&
          & blacsEnv%orbitalGrid%ctxt)
    end do

  end subroutine rows2sqr_${VAR}$

#:endfor

#:endif

end module dftbp_math_sparseblas
