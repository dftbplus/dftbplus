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
module dftbp_reksproperty

  use dftbp_accuracy
  use dftbp_blasroutines, only : gemm
  use dftbp_densitymatrix
  use dftbp_globalenv
  use dftbp_message
  use dftbp_sparse2dense
  use dftbp_rekscommon
  use dftbp_reksio
  use dftbp_reksvar, only : reksTypes

  implicit none

  private

  public :: getUnrelaxedDensMatAndTdp, getRelaxedDensMat, getRelaxedDensMatL
  public :: getDipoleIntegral, getDipoleMomentMatrix, getReksOsc

  contains

  !> Calculate unrelaxed density and transition density for target
  !> SA-REKS or SSR state (or L-th state)
  subroutine getUnrelaxedDensMatAndTdp(eigenvecs, overSqr, rhoSqrL, FONs, &
      & eigvecsSSR, Lpaired, Nc, Na, rstate, Lstate, reksAlg, tSSR, tTDP, &
      & unrelRhoSqr, unrelTdm)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> Dense density matrix for each microstate
    real(dp), intent(in) :: rhoSqrL(:,:,:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> eigenvectors from SA-REKS state
    real(dp), intent(in) :: eigvecsSSR(:,:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Target SSR state
    integer, intent(in) :: rstate

    !> Target microstate
    integer, intent(in) :: Lstate

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    logical, intent(in) :: tSSR

    !> Calculate transition dipole moments
    logical, intent(in) :: tTDP

    !> unrelaxed density matrix for target SSR or SA-REKS state
    real(dp), intent(out) :: unrelRhoSqr(:,:)

    !> unrelaxed transition density matrix between SSR or SA-REKS states
    real(dp), allocatable, intent(inout) :: unrelTdm(:,:,:)

    real(dp), allocatable :: rhoX(:,:,:)
    real(dp), allocatable :: rhoXdel(:,:,:)
    real(dp), allocatable :: tmpRho(:,:)
    real(dp), allocatable :: tmpMat(:,:)
    real(dp), allocatable :: tmpFilling(:,:)

    integer :: nOrb, nstates, nstHalf
    integer :: ii, ist, jst, kst, lst, ia, ib

    nOrb = size(eigenvecs,dim=1)
    nstates = size(eigvecsSSR,dim=1)
    nstHalf = nstates * (nstates - 1) / 2

    if (tSSR) then
      allocate(rhoX(nOrb,nOrb,nstates))
    else
      allocate(rhoX(nOrb,nOrb,1))
    end if
    if (Lstate == 0) then
      allocate(rhoXdel(nOrb,nOrb,nstHalf))
    end if
    allocate(tmpRho(nOrb,nOrb))
    allocate(tmpMat(nOrb,nOrb))
    allocate(tmpFilling(nOrb,nstates))

    ! core orbitals fillings
    tmpFilling(:,:) = 0.0_dp
    do ii = 1, Nc
      tmpFilling(ii,:) = 2.0_dp
    end do
    ! active orbitals fillings
    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getActiveFilling22_(FONs, Nc, tmpFilling)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    ! unrelaxed density matrix for SA-REKS or L-th state
    rhoX(:,:,:) = 0.0_dp
    if (tSSR) then
      do ist = 1, nstates
        call makeDensityMatrix(rhoX(:,:,ist), eigenvecs, tmpFilling(:,ist))
        call symmetrizeHS(rhoX(:,:,ist))
      end do
    else
      if (Lstate == 0) then
        call makeDensityMatrix(rhoX(:,:,1), eigenvecs, tmpFilling(:,rstate))
        call symmetrizeHS(rhoX(:,:,1))
      else
        ! find proper index for up+down in rhoSqrL
        if (Lstate <= Lpaired) then
          ii = Lstate
        else
          if (mod(Lstate,2) == 1) then
            ii = Lstate
          else
            ii = Lstate - 1
          end if
        end if
        rhoX(:,:,1) = rhoSqrL(:,:,1,ii)
      end if
    end if

    ! unrelaxed transition density matrix between SA-REKS states
    if (Lstate == 0) then
      rhoXdel(:,:,:) = 0.0_dp
      select case (reksAlg)
      case (reksTypes%noReks)
      case (reksTypes%ssr22)
        call getUnrelaxedTDM22_(eigenvecs, FONs, Nc, nstates, rhoXdel)
      case (reksTypes%ssr44)
        call error("SSR(4,4) is not implemented yet")
      end select
    end if

    ! Final unrelaxed density matrix for target state
    if (tSSR) then
      ! unrelRhoSqr is unrelaxed density matrix for target SSR state
      kst = 0
      unrelRhoSqr(:,:) = 0.0_dp
      do ist = 1, nstates
        do jst = ist, nstates
          if (ist == jst) then
            unrelRhoSqr(:,:) = unrelRhoSqr + eigvecsSSR(ist,rstate)**2 * rhoX(:,:,ist)
          else
            kst = kst + 1
            unrelRhoSqr(:,:) = unrelRhoSqr + 2.0_dp * eigvecsSSR(ist,rstate) * &
                & eigvecsSSR(jst,rstate) * rhoXdel(:,:,kst)
          end if
        end do
      end do
    else
      ! unrelRhoSqr is unrelaxed density matrix for target SA-REKS or L-th state
      unrelRhoSqr(:,:) = rhoX(:,:,1)
    end if

    ! Final unrelaxed transition density matrix between states
    if (tTDP .and. Lstate == 0) then
      if (tSSR) then
        ! unrelTdm is unrelaxed transition density matrix between SSR states
        unrelTdm(:,:,:) = 0.0_dp
        do lst = 1, nstHalf

          call getTwoIndices(nstates, lst, ia, ib, 1)

          kst = 0
          do ist = 1, nstates
            do jst = ist, nstates
              if (ist == jst) then
                unrelTdm(:,:,lst) = unrelTdm(:,:,lst) + rhoX(:,:,ist) &
                    & * eigvecsSSR(ist,ia) * eigvecsSSR(ist,ib)
              else
                kst = kst + 1
                ! <PPS|OSS> = <OSS|PPS>, etc
                unrelTdm(:,:,lst) = unrelTdm(:,:,lst) + rhoXdel(:,:,kst) &
                    & * ( eigvecsSSR(ist,ia) * eigvecsSSR(jst,ib) &
                    & + eigvecsSSR(jst,ia) * eigvecsSSR(ist,ib) )
              end if
            end do
          end do

        end do
      else
        ! unrelTdm is unrelaxed transition density matrix between SA-REKS states
        unrelTdm(:,:,:) = rhoXdel
      end if
    end if

    ! just calculate C^T*S*P*S*C = N, this will be diagonal.
    ! because, P = C*N*C^T, I = C^T*S*C, where
    ! P: density matrix, C: eigenvector, N: occupation number,
    ! T: transpose(real), I: identity matrix.
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, unrelRhoSqr, overSqr)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, overSqr, tmpMat)
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, tmpRho, eigenvecs)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, eigenvecs, tmpMat, transA='T')

    call printUnrelaxedFONs(tmpRho, rstate, Lstate, Nc, Na, tSSR)

  end subroutine getUnrelaxedDensMatAndTdp


  !> Calculate relaxed density for target SA-REKS or SSR state
  subroutine getRelaxedDensMat(eigenvecs, overSqr, unrelRhoSqr, ZT, omega, &
      & FONs, eigvecsSSR, SAweight, Rab, G1, Nc, Na, rstate, reksAlg, &
      & tSSR, tNAC, relRhoSqr)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> unrelaxed density matrix for target SSR or SA-REKS state
    real(dp), intent(in) :: unrelRhoSqr(:,:)

    !> solution of A * Z = X equation with X is XT
    real(dp), intent(in) :: ZT(:,:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> eigenvectors from SA-REKS state
    real(dp), intent(in) :: eigvecsSSR(:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> state-interaction term used in SSR gradients
    real(dp), allocatable, intent(in) :: Rab(:,:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Target SSR state
    integer, intent(in) :: rstate

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    logical, intent(in) :: tSSR

    !> Calculate nonadiabatic coupling vectors
    logical, intent(in) :: tNAC

    !> relaxed density matrix for target SSR or SA-REKS state
    real(dp), intent(out) :: relRhoSqr(:,:)

    real(dp), allocatable :: resRho(:,:)
    real(dp), allocatable :: resTdm(:,:,:)
    real(dp), allocatable :: tmpRho(:,:)
    real(dp), allocatable :: tmpMat(:,:)

    real(dp) :: fp, fq
    integer :: Nv, superN, nOrb, nstates, nstHalf
    integer :: ii, ist, jst, kst, pq, p, q

    superN = size(ZT,dim=1)
    nOrb = size(eigenvecs,dim=1)
    Nv = nOrb - Nc - Na
    nstates = size(eigvecsSSR,dim=1)
    nstHalf = nstates * (nstates - 1) / 2

    allocate(resRho(nOrb,nOrb))
    if (tSSR) then
      allocate(resTdm(nOrb,nOrb,nstHalf))
    end if
    allocate(tmpRho(nOrb,nOrb))
    allocate(tmpMat(nOrb,nOrb))

    ! a part of transition density matrix originating from the
    ! response of the orbital occupation numbers
    if (tSSR) then
      resTdm(:,:,:) = 0.0_dp
      select case (reksAlg)
      case (reksTypes%noReks)
      case (reksTypes%ssr22)
        call getResponseTDM22_(eigenvecs, FONs, SAweight, Rab, &
            & G1, Nc, resTdm)
      case (reksTypes%ssr44)
        call error("SSR(4,4) is not implemented yet")
      end select
    end if

    ! a part of relaxed density matrix originating from the
    ! response of the orbital occupation numbers with XT
    resRho(:,:) = 0.0_dp

    if (tNAC) then
      ii = rstate
    else
      ii = 1
    end if

    tmpRho(:,:) = 0.0_dp
    do pq = 1, superN
      ! assign index p and q from pq
      call assignIndex(Nc, Na, Nv, reksAlg, pq, p, q)
      ! assign average filling for pth orbital
      call assignFilling(FONs, SAweight, Nc, p, reksAlg, fp)
      ! assign average filling for qth orbital
      call assignFilling(FONs, SAweight, Nc, q, reksAlg, fq)
      tmpRho(p,q) = (fp - fq) * ZT(pq,ii)
    end do
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, tmpRho, eigenvecs, transB='T')
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, eigenvecs, tmpMat)

    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getResponseDM22_(eigenvecs, ZT(:,ii), tmpRho, omega, &
          & SAweight, G1, Nc, resRho)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    ! Final relaxed density matrix for target state
    relRhoSqr(:,:) = 0.0_dp
    if (tSSR) then
      kst = 0
      ! relRhoSqr is relaxed density matrix for target SSR state
      relRhoSqr(:,:) = unrelRhoSqr - 2.0_dp * resRho
      do ist = 1, nstates
        do jst = ist + 1, nstates
          kst = kst + 1
          relRhoSqr(:,:) = relRhoSqr + 2.0_dp * eigvecsSSR(ist,rstate) * &
              & eigvecsSSR(jst,rstate) * resTdm(:,:,kst)
        end do
      end do
    else
      ! relRhoSqr is relaxed density matrix for target SA-REKS state
      relRhoSqr(:,:) = unrelRhoSqr - 2.0_dp * resRho
    end if

    ! just calculate C^T*S*P*S*C = N, this will be diagonal.
    ! because, P = C*N*C^T, I = C^T*S*C, where
    ! P: density matrix, C: eigenvector, N: occupation number,
    ! T: transpose(real), I: identity matrix.
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, relRhosqr, overSqr)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, overSqr, tmpMat)
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, tmpRho, eigenvecs)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, eigenvecs, tmpMat, transA='T')

    call printRelaxedFONs(tmpRho, rstate, Nc, Na, tSSR)

  end subroutine getRelaxedDensMat


  !> Calculate relaxed density for target L-th microstate
  subroutine getRelaxedDensMatL(eigenvecs, rhoSqrL, overSqr, weight, &
      & SAweight, unrelRhoSqr, RmatL, ZT, omega, weightIL, G1, orderRmatL, &
      & Lpaired, Nc, Na, Lstate, reksAlg, relRhoSqr)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Dense density matrix for each microstate
    real(dp), intent(inout) :: rhoSqrL(:,:,:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> unrelaxed density matrix for target L-th state
    real(dp), intent(in) :: unrelRhoSqr(:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(in) :: RmatL(:,:,:,:)

    !> solution of A * Z = X equation with X is XT
    real(dp), intent(in) :: ZT(:,:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> modified weight of each microstate
    real(dp), intent(in) :: weightIL(:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Target microstate
    integer, intent(in) :: Lstate

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> relaxed density matrix for target L-th state
    real(dp), intent(out) :: relRhoSqr(:,:)

    real(dp), allocatable :: resRhoL(:,:)
    real(dp), allocatable :: tmpRho(:,:)
    real(dp), allocatable :: tmpMat(:,:)

    integer :: nOrb, tmpL, iL, Lmax

    nOrb = size(rhoSqrL,dim=1)
    Lmax = size(rhoSqrL,dim=4)

    allocate(resRhoL(nOrb,nOrb))
    allocate(tmpRho(nOrb,nOrb))
    allocate(tmpMat(nOrb,nOrb))

    ! rhoSqrL has (my_ud) component
    call qm2udL(rhoSqrL, Lpaired)

    ! resRhoL is response part for L-th state
    resRhoL(:,:) = 0.0_dp

    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getResponseDML22_(rhoSqrL, SAweight, ZT, omega, &
          & weightIL, G1, Lpaired, resRhoL)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    do iL = 1, Lmax
      ! find proper index for RmatL
      tmpL = orderRmatL(iL)
      resRhoL(:,:) = resRhoL - 2.0_dp * RmatL(:,:,tmpL,1) * weight(iL)
    end do

    ! relRhoSqr is relaxed density matrix for L-th state
    relRhoSqr(:,:) = unrelRhoSqr + resRhoL

    ! just calculate C^T*S*P*S*C = N, this will be diagonal.
    ! because, P = C*N*C^T, I = C^T*S*C, where
    ! P: density matrix, C: eigenvector, N: occupation number,
    ! T: transpose(real), I: identity matrix.
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, relRhosqr, overSqr)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, overSqr, tmpMat)
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, tmpRho, eigenvecs)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, eigenvecs, tmpMat, transA='T')

    call printRelaxedFONsL(tmpRho, Lstate, Nc, Na)

  end subroutine getRelaxedDensMatL


  !> Calculate dipole integral in DFTB formalism
  subroutine getDipoleIntegral(coord0, over, getAtomIndex, dipoleInt)

    !> central cell coordinates of atoms
    real(dp), intent(in) :: coord0(:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: over(:,:)

    !> get atom index from AO index
    integer, intent(in) :: getAtomIndex(:)

    !> dipole integral in DFTB formalism
    real(dp), intent(out) :: dipoleInt(:,:,:)

    real(dp), allocatable :: R(:,:)

    integer :: ii, mu, nu, nOrb

    nOrb = size(over,dim=1)

    allocate(R(nOrb,3))

    R(:,:) = 0.0_dp
    do mu = 1, nOrb
      ii = getAtomIndex(mu)
      R(mu,:) = coord0(:,ii)
    end do

    dipoleInt(:,:,:) = 0.0_dp
    do ii = 1, 3
      do mu = 1, nOrb
        do nu = 1, nOrb
          dipoleInt(mu,nu,ii) = 0.5_dp * over(mu,nu) * (R(mu,ii) + R(nu,ii))
        end do
      end do
    end do

  end subroutine getDipoleIntegral


  !> Calculate dipole moment using dipole integral and density matrix
  subroutine getDipoleMomentMatrix(Pmat, dipoleInt, dipole)

    !> density matrix related to dipole
    real(dp), intent(in) :: Pmat(:,:)

    !> dipole integral in DFTB formalism
    real(dp), intent(in) :: dipoleInt(:,:,:)

    !> resulting dipole moment
    real(dp), intent(out) :: dipole(:)

    integer :: ii

    dipole(:) = 0.0_dp
    do ii = 1, 3
      dipole(ii) = -sum(dipoleInt(:,:,ii)*Pmat(:,:))
    end do

  end subroutine getDipoleMomentMatrix


  !> get the oscillator strength between the states
  subroutine getReksOsc(tdp, energy)

    !> transition dipole moment between states
    real(dp), intent(in) :: tdp(:,:)

    !> energy of states
    real(dp), intent(in) :: energy(:)

    real(dp) :: osc
    integer :: ia, ib, ist, nstates, nstHalf

    nstates = size(energy,dim=1)
    nstHalf = size(tdp,dim=2)

    write(stdOut,*)
    write(stdOut,'(A)') " Oscillator Strength (au)"
    do ist = 1, nstHalf

      call getTwoIndices(nstates, ist, ia, ib, 1)

      osc = 2.0_dp / 3.0_dp * (energy(ib) - energy(ia)) * sum(tdp(:,ist)**2)

      write(stdOut,'(A4,I1,A6,I1,A5)',advance="no") " ( S", ia - 1, " <-> S", ib - 1, " ) : "
      write(stdOut,'(1(f12.6))') osc

    end do

  end subroutine getReksOsc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate filling information for active space in (2,2) case
  subroutine getActiveFilling22_(FONs, Nc, tmpFilling)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> temporary average filling for SA-REKS state
    real(dp), intent(inout) :: tmpFilling(:,:)

    real(dp) :: n_a, n_b
    integer :: a, b, nstates

    nstates = size(tmpFilling,dim=2)

    n_a = FONs(1,1)
    n_b = FONs(2,1)
    a = Nc + 1
    b = Nc + 2

    ! PPS state fillings
    tmpFilling(a,1) = n_a
    tmpFilling(b,1) = n_b
    ! OSS state fillings
    tmpFilling(a,2) = 1.0_dp
    tmpFilling(b,2) = 1.0_dp
    if (nstates == 3) then
      ! DES state fillings
      tmpFilling(a,3) = n_b
      tmpFilling(b,3) = n_a
    end if

  end subroutine getActiveFilling22_


  !> Calculate unrelaxed transition density between SA-REKS states in (2,2) case
  subroutine getUnrelaxedTDM22_(eigenvecs, FONs, Nc, nstates, rhoXdel)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of states
    integer, intent(in) :: nstates

    !> unrelaxed transition density matrix between SA-REKS states
    real(dp), intent(inout) :: rhoXdel(:,:,:)

    real(dp) :: n_a, n_b
    integer :: mu, nu, nOrb, a, b

    nOrb = size(eigenvecs,dim=1)

    n_a = FONs(1,1)
    n_b = FONs(2,1)
    a = Nc + 1
    b = Nc + 2

    do mu = 1, nOrb
      do nu = 1, nOrb
        rhoXdel(nu,mu,1) = eigenvecs(mu,a)*eigenvecs(nu,b) * &
            & (sqrt(n_a) - sqrt(n_b))
      end do
    end do
    if (nstates == 3) then
      do mu = 1, nOrb
        do nu = 1, nOrb
          rhoXdel(nu,mu,3) = eigenvecs(mu,a)*eigenvecs(nu,b) * &
              & (sqrt(n_a) + sqrt(n_b))
        end do
      end do
    end if

  end subroutine getUnrelaxedTDM22_


  !> Calculate response part of relaxed density for transition density contribution
  subroutine getResponseTDM22_(eigenvecs, FONs, SAweight, Rab, &
      & G1, Nc, resTdm)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> state-interaction term used in SSR gradients
    real(dp), intent(in) :: Rab(:,:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> a part of transition density matrix originating from the
    !> response of the orbital occupation numbers
    real(dp), intent(inout) :: resTdm(:,:,:)

    real(dp) :: n_a, n_b
    integer :: mu, nu, nOrb, a, b, nstHalf

    nstHalf = size(resTdm,dim=3)
    nOrb = size(eigenvecs,dim=1)

    n_a = FONs(1,1); n_b = FONs(2,1)
    a = Nc + 1; b = Nc + 2

    do mu = 1, nOrb
      do nu = 1, nOrb
        resTdm(nu,mu,1) = resTdm(nu,mu,1) + SAweight(1) * &
            & ( (n_a-1.0_dp)*sqrt(n_a) - (n_b-1.0_dp)*sqrt(n_b) ) * &
            & eigenvecs(mu,a) * eigenvecs(nu,b)
        resTdm(nu,mu,1) = resTdm(nu,mu,1) - 2.0_dp * G1 * &
            & Rab(1,2) * ( eigenvecs(mu,a) * eigenvecs(nu,a) - &
            & eigenvecs(mu,b) * eigenvecs(nu,b) )
      end do
    end do
    if (nstHalf == 3) then
      do mu = 1, nOrb
        do nu = 1, nOrb
          resTdm(nu,mu,3) = resTdm(nu,mu,3) + SAweight(1) * &
              & ( (n_a-1.0_dp)*sqrt(n_a) + (n_b-1.0_dp)*sqrt(n_b) ) * &
              & eigenvecs(mu,a) * eigenvecs(nu,b)
          resTdm(nu,mu,3) = resTdm(nu,mu,3) - 2.0_dp * G1 * &
              & Rab(2,3) * ( eigenvecs(mu,a) * eigenvecs(nu,a) - &
              & eigenvecs(mu,b) * eigenvecs(nu,b) )
        end do
      end do
    end if

  end subroutine getResponseTDM22_


  !> Calculate response part of relaxed density for density contribution
  subroutine getResponseDM22_(eigenvecs, ZT, tmpZ, omega, &
      & SAweight, G1, Nc, resRho)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> solution of A * Z = X equation with X is XT
    real(dp), intent(in) :: ZT(:)

    !> temporary matrix including ZT in MO basis
    real(dp), intent(in) :: tmpZ(:,:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> a part of relaxed density matrix originating from the
    !> response of the orbital occupation numbers with XT
    real(dp), intent(out) :: resRho(:,:)

    real(dp) :: tmpValue
    integer :: a, b, mu, nu, nOrb

    nOrb = size(eigenvecs,dim=1)

    a = Nc + 1
    b = Nc + 2

    tmpValue = sum(ZT(:)*omega(:))
    do mu = 1, nOrb
      do nu = 1, nOrb
        resRho(nu,mu) = tmpZ(mu,nu) - SAweight(1) * G1 * tmpValue * &
            & (eigenvecs(mu,a)*eigenvecs(nu,a) - eigenvecs(mu,b)*eigenvecs(nu,b))
      end do
    end do

  end subroutine getResponseDM22_


  !> Calculate response part of L-th relaxed density for density contribution
  subroutine getResponseDML22_(rhoSqrL, SAweight, ZT, omega, &
      & weightIL, G1, Lpaired, resRhoL)

    !> Dense density matrix for each microstate
    real(dp), intent(in) :: rhoSqrL(:,:,:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> solution of A * Z = X equation with X is XT
    real(dp), intent(in) :: ZT(:,:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> modified weight of each microstate
    real(dp), intent(in) :: weightIL(:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> response part of relaxed density matrix for target L-th state
    real(dp), intent(inout) :: resRhoL(:,:)

    real(dp) :: tmpValue
    integer :: tmpL, iL, Lmax

    Lmax = size(rhoSqrL,dim=4)

    tmpValue = sum(ZT(:,1)*omega(:))
    do iL = 1, Lmax

      ! find proper index for down spin in rhoSqrL
      if (iL <= Lpaired) then
        tmpL = iL
      else
        if (mod(iL,2) == 1) then
          tmpL = iL + 1
        else
          tmpL = iL - 1
        end if
      end if
      resRhoL(:,:) = resRhoL + SAweight(1) * G1 * tmpValue * &
          & weightIL(iL) * (rhoSqrL(:,:,1,iL) + rhoSqrL(:,:,1,tmpL))

    end do

  end subroutine getResponseDML22_


end module dftbp_reksproperty
