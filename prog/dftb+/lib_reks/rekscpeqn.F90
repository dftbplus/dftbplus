!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> REKS and SI-SA-REKS formulation in DFTB as developed by Lee et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet appart from Gamma point.
!> * Orbital potentials or spin-orbit or external E-field does not work yet.
!> * Only for closed shell system.
!> * Onsite corrections are not included in this version
module dftbp_rekscpeqn

  use dftbp_accuracy
  use dftbp_blasroutines, only : gemm, gemv
  use dftbp_densedescr
  use dftbp_environment
  use dftbp_globalenv
  use dftbp_message
  use dftbp_orbitals
  use dftbp_periodic
  use dftbp_rekscommon
  use dftbp_reksgrad, only : getRmat, getZmat, getQ2mat

  implicit none

  private

  public :: CGgrad

  contains

  !> solve CP-REKS equation by using conjugate-gradient method
  subroutine CGgrad(env, denseDesc, neighbourList, nNeighbourSK, iSparseStart, &
      & img2CentCell, orb, XT, A1e, A1ePre, HxcSqrS, HxcSqrD, HxcHalfS, &
      & HxcHalfD, HxcSpS, HxcSpD, Fc, Fa, omega, SAweight, FONs, G1, GammaAO, &
      & SpinAO, LrGammaAO, overSqr, over, eigenvecs, fillingL, weight, &
      & ConvergeLimit, orderRmatL, getDenseAO, Lpaired, Nc, Na, maxIter, Glevel, &
      & reksAlg, tSaveMem, isRangeSep, ZT, RmatL, ZmatL, Q2mat)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb


    !> SA-REKS state vector
    real(dp), intent(in) :: XT(:)

    !> super A hessian matrix with one-electron term in front of orbital derivatives
    real(dp), allocatable, intent(in) :: A1e(:,:)

    !> preconditioner of super A hessian matrix with one-electron term in front of orbital
    !> derivatives
    real(dp), allocatable, intent(in) :: A1ePre(:,:)

    !> Hartree-XC kernel with dense form with same spin part
    real(dp), allocatable, intent(in) :: HxcSqrS(:,:,:,:)

    !> Hartree-XC kernel with dense form with different spin part
    real(dp), allocatable, intent(in) :: HxcSqrD(:,:,:,:)

    !> Hartree-XC kernel with half dense form with same spin part
    real(dp), allocatable, intent(in) :: HxcHalfS(:,:)

    !> Hartree-XC kernel with half dense form with different spin part
    real(dp), allocatable, intent(in) :: HxcHalfD(:,:)

    !> Hartree-XC kernel with sparse form with same spin part
    real(dp), allocatable, intent(in) :: HxcSpS(:,:)

    !> Hartree-XC kernel with sparse form with different spin part
    real(dp), allocatable, intent(in) :: HxcSpD(:,:)


    !> dense fock matrix for core orbitals
    real(dp), intent(in) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: Fa(:,:,:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1


    !> scc gamma integrals in AO basis
    real(dp), intent(in) :: GammaAO(:,:)

    !> spin W in AO basis
    real(dp), intent(in) :: SpinAO(:,:)

    !> long-range gamma integrals in AO basis
    real(dp), allocatable, intent(in) :: LrGammaAO(:,:)


    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Tolerance used in calculation of gradient with PCG and CG
    real(dp), intent(in) :: ConvergeLimit


    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> get dense AO index from sparse AO array
    integer, intent(in) :: getDenseAO(:,:)


    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Maximum iteration used in calculation of gradient with PCG and CG
    integer, intent(in) :: maxIter

    !> Algorithms to calculate analytic gradients
    integer, intent(in) :: Glevel

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> Save 'A' and 'Hxc' to memory in gradient calculation
    logical, intent(in) :: tSaveMem

    !> Whether to run a range separated calculation
    logical, intent(in) :: isRangeSep


    !> solution of A * Z = X equation with X is XT
    real(dp), intent(out) :: ZT(:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(out) :: RmatL(:,:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(out) :: ZmatL(:,:,:)

    !> auxiliary matrix in MO basis related to SA-REKS term
    real(dp), intent(out) :: Q2mat(:,:)


    real(dp), allocatable :: shift1e(:)          ! A1e * ZT
    real(dp), allocatable :: shift2e(:)          ! A2e * ZT
    real(dp), allocatable :: r0(:), r1(:)        ! residual vector
    real(dp), allocatable :: z0(:), z1(:)        ! PCG residual vector
    real(dp), allocatable :: p0(:), p1(:)        ! direction vector

    real(dp) :: rNorm1, rNorm2, alpha, beta, eps
    integer :: iter, superN

    superN = size(XT,dim=1)

    allocate(shift1e(superN),shift2e(superN))
    allocate(r0(superN),r1(superN))
    allocate(z0(superN),z1(superN))
    allocate(p0(superN),p1(superN))

    write(stdOut,"(A)") repeat("-", 82)

    ! initial guess for Z vector
    ! initial Z_initial = A_pre^{-1} * X
    if (tSaveMem) then
      ZT(:) = 0.0_dp
      call gemv(ZT, A1ePre, XT)
    else
      call shiftAY1ePre_(XT, Fc, Fa, omega, SAweight, FONs, G1, &
          & Nc, Na, Glevel, reksAlg, ZT)
    end if

    ! initial setting for r0, z0, p0 vectors
    ! 1-electron part
    if (tSaveMem) then
      shift1e(:) = 0.0_dp
      call gemv(shift1e, A1e, ZT)
    else
      call shiftAY1e_(ZT, Fc, Fa, omega, SAweight, FONs, G1, &
          & Nc, Na, reksAlg, shift1e)
    end if
    ! 2-electron part
    call getRmat(eigenvecs, ZT, fillingL, Nc, Na, reksAlg, RmatL)
    call getZmat(env, denseDesc, neighbourList, nNeighbourSK, &
        & iSparseStart, img2CentCell, orb, RmatL, HxcSqrS, HxcSqrD, &
        & HxcHalfS, HxcHalfD, HxcSpS, HxcSpD, overSqr, over, &
        & GammaAO, SpinAO, LrGammaAO, orderRmatL, getDenseAO, &
        & Lpaired, Glevel, tSaveMem, isRangeSep, ZmatL)
    call shiftAY2e_(ZmatL, eigenvecs, fillingL, weight, &
        & Nc, Na, reksAlg, shift2e)

    ! calculate r0 vector
    r0(:) = XT - shift1e - shift2e
    ! calculate z0, p0 vector
    if (tSaveMem) then
      z0(:) = 0.0_dp
      call gemv(z0, A1ePre, r0)
    else
      call shiftAY1ePre_(r0, Fc, Fa, omega, SAweight, FONs, G1, &
          & Nc, Na, Glevel, reksAlg, z0)
    end if
    p0(:) = z0

    iter = 0; eps = 0.0_dp
    write(stdOut,'(2x,a)') 'CG solver: Constructing Y initial guess'

    CGsolver: do iter = 1, maxIter

      ! Construct (A1e + A2e) * P
      ! 1-electron part
      if (tSaveMem) then
        shift1e(:) = 0.0_dp
        call gemv(shift1e, A1e, p0)
      else
        call shiftAY1e_(p0, Fc, Fa, omega, SAweight, FONs, G1, &
            & Nc, Na, reksAlg, shift1e)
      end if
      ! 2-electron part
      call getRmat(eigenvecs, p0, fillingL, Nc, Na, reksAlg, RmatL)
      call getZmat(env, denseDesc, neighbourList, nNeighbourSK, &
          & iSparseStart, img2CentCell, orb, RmatL, HxcSqrS, HxcSqrD, &
          & HxcHalfS, HxcHalfD, HxcSpS, HxcSpD, overSqr, over, &
          & GammaAO, SpinAO, LrGammaAO, orderRmatL, getDenseAO, &
          & Lpaired, Glevel, tSaveMem, isRangeSep, ZmatL)
      call shiftAY2e_(ZmatL, eigenvecs, fillingL, weight, &
          & Nc, Na, reksAlg, shift2e)

      ! compute step length
      rNorm1 = sum( r0(:)*z0(:) )
      if (rNorm1 == 0.0_dp) then
        alpha = 0.0_dp
      else
        alpha = rNorm1 / sum( p0(:)*(shift1e(:)+shift2e(:)) )
      end if
      ! update the approximate solution
      ZT(:) = ZT + alpha * p0
      ! update the residual
      r1(:) = r0 - alpha * (shift1e(:)+shift2e(:))
      ! solve new PCG residual
      if (tSaveMem) then
        z1(:) = 0.0_dp
        call gemv(z1, A1ePre, r1)
      else
        call shiftAY1ePre_(r1, Fc, Fa, omega, SAweight, FONs, G1, &
            & Nc, Na, Glevel, reksAlg, z1)
      end if
      ! compute a gradient correction factor
      rNorm2 = sum( r1(:)*z1(:) )
      if (rNorm2 == 0.0_dp) then
        beta = 0.0_dp
      else
        beta = rNorm2 / rNorm1
      end if
      p1(:) = z1 + beta * p0

      ! calculate square residual for current iteration
      eps = sum( r1(:)*r1(:) )

      ! show current iteration
      write(stdOut,'(2x,a,1x,i4,4x,a,F18.12)') &
          & 'CG solver: Iteration', iter, 'eps =', eps

      ! convergence check
      if (eps > ConvergeLimit) then
        ! update new variables to old variables
        r0(:) = r1
        z0(:) = z1
        p0(:) = p1
        if (iter == maxIter) then
          write(stdOut,'(2x,a,i4,a)') &
              & 'Warning! Maximum number of iterations (', maxIter, &
              & ') is exceeded in CG solver'
          call error("Increase the maximum number of iterations")
        end if
      else
        write(stdOut,'(2x,a,1x,i4,1x,a)') &
            & 'Convergence reached in CG solver after', iter, 'iterations'
        exit CGsolver
      end if

    end do CGsolver

    ! converged R, Z, Q2 value
    call getRmat(eigenvecs, ZT, fillingL, Nc, Na, reksAlg, RmatL)
    call getZmat(env, denseDesc, neighbourList, nNeighbourSK, &
        & iSparseStart, img2CentCell, orb, RmatL, HxcSqrS, HxcSqrD, &
        & HxcHalfS, HxcHalfD, HxcSpS, HxcSpD, overSqr, over, &
        & GammaAO, SpinAO, LrGammaAO, orderRmatL, getDenseAO, &
        & Lpaired, Glevel, tSaveMem, isRangeSep, ZmatL)
    call getQ2mat(eigenvecs, fillingL, weight, ZmatL, Q2mat)
    write(stdOut,'(2x,a)') 'CG solver: Calculating converged R, Z, Q2 matrix'
    write(stdOut,"(A)") repeat("-", 82)

  end subroutine CGgrad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate A1e * Y shift vectors without saving super A1e matrix
  subroutine shiftAY1e_(Y, Fc, Fa, omega, SAweight, FONs, G1, &
      & Nc, Na, reksAlg, shift1e)

    !> trial vector for soulution
    real(dp), intent(in) :: Y(:)

    !> dense fock matrix for core orbitals
    real(dp), intent(in) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: Fa(:,:,:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> computed (A1e * Y) shift vector
    real(dp), intent(out) :: shift1e(:)

    real(dp), allocatable :: tmpA(:)
    real(dp) :: e1, e2
    integer :: nOrb, Nv, superN, ij, pq, i, j, p, q

    superN = size(Y,dim=1)
    nOrb = size(Fc,dim=1)
    Nv = nOrb - Nc - Na

    allocate(tmpA(superN))

    shift1e(:) = 0.0_dp
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,p,q,e1,e2,tmpA) SCHEDULE(RUNTIME)
    do ij = 1, superN

      ! assign index i and j from ij
      call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)

      tmpA(:) = 0.0_dp
      do pq = 1, superN

        ! assign index p and q from pq
        call assignIndex(Nc, Na, Nv, reksAlg, pq, p, q)

        if (ij <= pq) then

          ! get lagrange multipliers with delta function
          if (i == p) then
            e1 = 0.0_dp; e2 = 0.0_dp;
            call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, q, &
                & 1, reksAlg, e1, e2)
            tmpA(pq) = tmpA(pq) + 0.5_dp*(e1 - e2)
          end if

          if (j == q) then
            e1 = 0.0_dp; e2 = 0.0_dp;
            call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, p, &
                & 2, reksAlg, e1, e2)
            tmpA(pq) = tmpA(pq) - 0.5_dp*(e1 - e2)
          end if

          if (i == q) then
            e1 = 0.0_dp; e2 = 0.0_dp;
            call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, p, &
                & 1, reksAlg, e1, e2)
            tmpA(pq) = tmpA(pq) - 0.5_dp*(e1 - e2)
          end if

          if (j == p) then
            e1 = 0.0_dp; e2 = 0.0_dp;
            call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, q, &
                & 2, reksAlg, e1, e2)
            tmpA(pq) = tmpA(pq) + 0.5_dp*(e1 - e2)
          end if

          ! SAweight(1) is equal to Wgss
          tmpA(pq) = tmpA(pq) - SAweight(1) * G1 * omega(ij) * omega(pq)

        else

          ! get lagrange multipliers with delta function
          if (p == i) then
            e1 = 0.0_dp; e2 = 0.0_dp;
            call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, p, q, j, &
                & 1, reksAlg, e1, e2)
            tmpA(pq) = tmpA(pq) + 0.5_dp*(e1 - e2)
          end if

          if (q == j) then
            e1 = 0.0_dp; e2 = 0.0_dp;
            call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, p, q, i, &
                & 2, reksAlg, e1, e2)
            tmpA(pq) = tmpA(pq) - 0.5_dp*(e1 - e2)
          end if

          if (p == j) then
            e1 = 0.0_dp; e2 = 0.0_dp;
            call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, p, q, i, &
                & 1, reksAlg, e1, e2)
            tmpA(pq) = tmpA(pq) - 0.5_dp*(e1 - e2)
          end if

          if (q == i) then
            e1 = 0.0_dp; e2 = 0.0_dp;
            call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, p, q, j, &
                & 2, reksAlg, e1, e2)
            tmpA(pq) = tmpA(pq) + 0.5_dp*(e1 - e2)
          end if

          ! SAweight(1) is equal to Wgss
          tmpA(pq) = tmpA(pq) - SAweight(1) * G1 * omega(ij) * omega(pq)

        end if

      end do

      shift1e(ij) = sum(tmpA(:)*Y(:))

    end do
!$OMP END PARALLEL DO

  end subroutine shiftAY1e_


  !> Calculate A1ePre * Y shift vectors without saving super A1ePre matrix
  subroutine shiftAY1ePre_(Y, Fc, Fa, omega, SAweight, FONs, G1, &
      & Nc, Na, Glevel, reksAlg, shift1ePre)

    !> trial vector for soulution
    real(dp), intent(in) :: Y(:)

    !> dense fock matrix for core orbitals
    real(dp), intent(in) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: Fa(:,:,:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Algorithms to calculate analytic gradients
    integer, intent(in) :: Glevel

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> computed (A1ePre * Y) shift vector
    real(dp), intent(out) :: shift1ePre(:)

    real(dp), allocatable :: tmpApre(:)
    real(dp) :: e1, e2
    integer :: nOrb, superN, Nv, ij, i, j

    superN = size(Y,dim=1)
    nOrb = size(Fc,dim=1)
    Nv = nOrb - Nc - Na

    allocate(tmpApre(superN))

    tmpApre(:) = 0.0_dp
    shift1ePre(:) = 0.0_dp
    do ij = 1, superN

      ! assign index i and j from ij
      call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)

      if (Glevel == 1) then

        ! get lagrange multipliers with delta function
        e1 = 0.0_dp; e2 = 0.0_dp;
        call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, j, &
            & 1, reksAlg, e1, e2)
        tmpApre(ij) = tmpApre(ij) + 0.5_dp*(e1 - e2)

        e1 = 0.0_dp; e2 = 0.0_dp;
        call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, i, &
            & 2, reksAlg, e1, e2)
        tmpApre(ij) = tmpApre(ij) - 0.5_dp*(e1 - e2)

        ! SAweight(1) is equal to Wgss
        tmpApre(ij) = tmpApre(ij) - SAweight(1) * G1 * omega(ij) * omega(ij)

      else if (Glevel == 2) then

        tmpApre(ij) = 1.0_dp

      end if

      ! check singularity for preconditioner
      if (abs(tmpApre(ij)) <= epsilon(1.0_dp)) then
        write(stdOut,'(A,f15.8)') " Current preconditioner value = ", tmpApre(ij)
        call error("A singularity exists in preconditioner for PCG, set Preconditioner = No")
      end if

      ! preconditioner part for CG
      tmpApre(ij) = 1.0_dp / tmpApre(ij)

      shift1ePre(ij) = tmpApre(ij) * Y(ij)

    end do

  end subroutine shiftAY1ePre_


  !> Calculate A2e * Y shift vectors
  subroutine shiftAY2e_(ZmatL, eigenvecs, fillingL, weight, &
      & Nc, Na, reksAlg, shift2e)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(in) :: ZmatL(:,:,:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> computed (A2e * Y) shift vector
    real(dp), intent(out) :: shift2e(:)

    real(dp), allocatable :: tmpMat(:,:)
    real(dp), allocatable :: Zmo(:,:)

    real(dp) :: fac
    integer :: nOrb, superN, Nv, iL, p, q, pq, Lmax

    superN = size(shift2e,dim=1)
    nOrb = size(fillingL,dim=1)
    Lmax = size(fillingL,dim=3)
    Nv = nOrb - Nc - Na

    allocate(tmpMat(nOrb,nOrb))
    allocate(Zmo(nOrb,nOrb))

    shift2e(:) = 0.0_dp
    do iL = 1, Lmax
      tmpMat(:,:) = 0.0_dp
      Zmo(:,:) = 0.0_dp
      call gemm(tmpMat, ZmatL(:,:,iL), eigenvecs(:,:,1))
      call gemm(Zmo, eigenvecs(:,:,1), tmpMat, transA='T')
      do pq = 1, superN
        ! assign index p and q from pq
        call assignIndex(Nc, Na, Nv, reksAlg, pq, p, q)
        fac = 2.0_dp * weight(iL) * Zmo(q,p) &
           & * (fillingL(p,1,iL) - fillingL(q,1,iL))
        shift2e(pq) = shift2e(pq) + fac
      end do
    end do

  end subroutine shiftAY2e_


end module dftbp_rekscpeqn
