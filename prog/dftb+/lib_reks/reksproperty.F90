!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
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
! TODO
!> * Dispersion would be combined with REKS
module dftbp_reksproperty

  use dftbp_accuracy
  use dftbp_blasroutines, only : gemm
  use dftbp_densitymatrix
  use dftbp_globalenv
  use dftbp_mainio
  use dftbp_message
  use dftbp_sparse2dense
  use dftbp_rekscommon
  use dftbp_reksvar

  implicit none

  private

  public :: getUnrelaxedDMandTDP, getRelaxedDM, getRelaxedDML
  public :: getDipoleIntegral, getDipoleMomentMatrix, getReksOsc

  contains

  !> Calculate unrelaxed density and transition density for target
  !> SA-REKS or SSR state (or L-th state)
  ! TODO : variable name should be changed!
  subroutine getUnrelaxedDMandTDP(eigenvecs, self)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    real(dp), allocatable :: P_X_o(:,:,:)
    real(dp), allocatable :: P_X_del_o(:,:,:)
    real(dp), allocatable :: tmpRho(:,:)
    real(dp), allocatable :: tmpMat(:,:)
    real(dp), allocatable :: tmpFilling(:,:)

    real(dp) :: tmp
    integer :: nOrb, nstHalf
    integer :: ii, ist, jst, kst, lst, ia, ib

    nOrb = size(eigenvecs,dim=1)
    nstHalf = self%nstates * (self%nstates - 1) / 2

    if (self%useSSR == 1) then
      allocate(P_X_o(nOrb,nOrb,self%nstates))
    else
      allocate(P_X_o(nOrb,nOrb,1))
    end if
    if (self%Lstate == 0) then
      allocate(P_X_del_o(nOrb,nOrb,nstHalf))
    end if
    allocate(tmpRho(nOrb,nOrb))
    allocate(tmpMat(nOrb,nOrb))
    allocate(tmpFilling(nOrb,self%nstates))

    ! core orbitals fillings
    tmpFilling(:,:) = 0.0_dp
    do ii = 1, self%Nc
      tmpFilling(ii,:) = 2.0_dp
    end do
    ! active orbitals fillings
    if (self%tSSR22) then
      call getActiveFilling22_(self%FONs, self%Nc, tmpFilling)
    else if (self%tSSR44) then
      call error("SSR(4,4) is not implemented yet")
    end if

    ! unrelaxed density matrix for SA-REKS or L-th state
    P_X_o(:,:,:) = 0.0_dp
    if (self%useSSR == 1) then
      do ist = 1, self%nstates
        call makeDensityMatrix(P_X_o(:,:,ist), eigenvecs, tmpFilling(:,ist))
        call symmetrizeHS(P_X_o(:,:,ist))
      end do
    else
      if (self%Lstate == 0) then
        call makeDensityMatrix(P_X_o(:,:,1), eigenvecs, tmpFilling(:,self%rstate))
        call symmetrizeHS(P_X_o(:,:,ist))
      else
        ! find proper index for up+down in self%dm_L
        if (self%Lstate <= self%Lpaired) then
          ii = self%Lstate
        else
          if (mod(self%Lstate,2) == 1) then
            ii = self%Lstate
          else
            ii = self%Lstate - 1
          end if
        end if
        P_X_o(:,:,1) = self%rhoSqrL(:,:,1,ii)
      end if
    end if

    ! unrelaxed transition density matrix between SA-REKS states
    if (self%Lstate == 0) then
      P_X_del_o(:,:,:) = 0.0_dp
      if (self%tSSR22) then
        call getUnrelaxedTDM22_(eigenvecs, self%FONs, self%Nc, self%nstates, P_X_del_o)
      else if (self%tSSR44) then
        call error("SSR(4,4) is not implemented yet")
      end if
    end if

    ! Final unrelaxed density matrix for target state
    if (self%useSSR == 1) then
      ! self%unrelRhoSqr is unrelaxed density matrix for target SSR state
      kst = 0
      self%unrelRhoSqr(:,:) = 0.0_dp
      do ist = 1, self%nstates
        do jst = ist, self%nstates
          if (ist == jst) then
            self%unrelRhoSqr(:,:) = self%unrelRhoSqr + self%eigvecsSSR(ist,self%rstate)**2 * P_X_o(:,:,ist)
          else
            kst = kst + 1
            self%unrelRhoSqr(:,:) = self%unrelRhoSqr + 2.0_dp * self%eigvecsSSR(ist,self%rstate) * &
                & self%eigvecsSSR(jst,self%rstate) * P_X_del_o(:,:,kst)
          end if
        end do
      end do
    else
      ! self%unrelRhoSqr is unrelaxed density matrix for target SA-REKS or L-th state
      self%unrelRhoSqr(:,:) = P_X_o(:,:,1)
    end if

    ! Final unrelaxed transition density matrix between states
    if (self%tTDP .and. self%Lstate == 0) then
      if (self%useSSR == 1) then
        ! self%unrelTdm is unrelaxed transition density matrix between SSR states
        self%unrelTdm(:,:,:) = 0.0_dp
        do lst = 1, nstHalf

          ! (ia,ib) = (1,2) (1,3) (2,3) ...
          ! TODO
          tmp = ( dble(2.0_dp*self%nstates+1.0_dp) - dsqrt( (2.0_dp*self%nstates+ &
              & 1.0_dp)**2.0_dp - 8.0_dp*(self%nstates+lst) ) )/2.0_dp
          ia = int( tmp )
          if( (tmp - dble(ia)) < 1.0E-8_dp ) ia = ia - 1
          ib = ia**2/2 + ia/2 - self%nstates*ia + self%nstates + lst
          if( mod(ia,2) == 1 ) ib = ib + 1

          kst = 0
          do ist = 1, self%nstates
            do jst = ist, self%nstates
              if (ist == jst) then
                self%unrelTdm(:,:,lst) = self%unrelTdm(:,:,lst) + P_X_o(:,:,ist) &
                    & * self%eigvecsSSR(ist,ia) * self%eigvecsSSR(ist,ib)
              else
                kst = kst + 1
                ! <PPS|OSS> = <OSS|PPS>, etc
                self%unrelTdm(:,:,lst) = self%unrelTdm(:,:,lst) + P_X_del_o(:,:,kst) &
                    & * ( self%eigvecsSSR(ist,ia) * self%eigvecsSSR(jst,ib) &
                    & + self%eigvecsSSR(jst,ia) * self%eigvecsSSR(ist,ib) )
              end if
            end do
          end do

        end do
      else
        ! self%unrelTdm is unrelaxed transition density matrix between SA-REKS states
        self%unrelTdm(:,:,:) = P_X_del_o(:,:,:)
      end if
    end if

    ! just calculate C^T*S*P*S*C = N, this will be diagonal.
    ! because, P = C*N*C^T, I = C^T*S*C, where
    ! P: density matrix, C: eigenvector, N: occupation number,
    ! T: transpose(real), I: identity matrix.
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, self%unrelRhoSqr, self%overSqr)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, self%overSqr, tmpMat)
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, tmpRho, eigenvecs)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, eigenvecs, tmpMat, transA='T')

    call printUnrelaxedFONs(tmpRho, self%useSSR, self%rstate,&
        & self%Lstate, self%Nc, self%Na)

  end subroutine getUnrelaxedDMandTDP


  !> Calculate relaxed density for target SA-REKS or SSR state
  subroutine getRelaxedDM(eigenvecs, overSqr, unrelRhoSqr, ZT, omega, &
      & FONs, eigvecsSSR, SAweight, Rab, G1, Nc, Na, rstate, useSSR, &
      & tNAC, tSSR22, tSSR44, relRhoSqr)

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
    real(dp), intent(in) :: Rab(:,:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Target SSR state
    integer, intent(in) :: rstate

    !> Calculate SSR state (SI term is included)
    integer, intent(in) :: useSSR

    !> Calculate nonadiabatic coupling vectors
    logical, intent(in) :: tNAC

    !> Calculate DFTB/SSR(2,2) formalism
    logical, intent(in) :: tSSR22

    !> Calculate DFTB/SSR(4,4) formalism
    logical, intent(in) :: tSSR44

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
    if (useSSR == 1) then
      allocate(resTdm(nOrb,nOrb,nstHalf))
    end if
    allocate(tmpRho(nOrb,nOrb))
    allocate(tmpMat(nOrb,nOrb))

    ! a part of transition density matrix originating from the
    ! response of the orbital occupation numbers
    if (useSSR == 1) then
      resTdm(:,:,:) = 0.0_dp
      if (tSSR22) then
        call getResponseTDM22_(eigenvecs, FONs, SAweight, Rab, &
            & G1, Nc, resTdm)
      else if (tSSR44) then
        call error("SSR(4,4) is not implemented yet")
      end if
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
      call assignIndex(Nc, Na, Nv, tSSR22, tSSR44, pq, p, q)
      ! assign average filling for pth orbital
      call assignFilling(FONs, SAweight, Nc, p, tSSR22, tSSR44, fp)
      ! assign average filling for qth orbital
      call assignFilling(FONs, SAweight, Nc, q, tSSR22, tSSR44, fq)
      tmpRho(p,q) = (fp - fq) * ZT(pq,ii)
    end do
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, tmpRho, eigenvecs, transB='T')
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, eigenvecs, tmpMat)

    if (tSSR22) then
      call getResponseDM22_(eigenvecs, ZT(:,ii), tmpRho, omega, &
          & SAweight, G1, Nc, resRho)
    else if (tSSR44) then
      call error("SSR(4,4) is not implemented yet")
    end if

    ! Final relaxed density matrix for target state
    relRhoSqr(:,:) = 0.0_dp
    if (useSSR == 1) then
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

    call printRelaxedFONs(tmpRho, useSSR, rstate, Nc, Na)

  end subroutine getRelaxedDM


  !> Calculate relaxed density for target L-th microstate
  subroutine getRelaxedDML(eigenvecs, rhoSqrL, overSqr, weight, &
      & SAweight, unrelRhoSqr, RmatL, ZT, omega, weightIL, G1, orderRmatL, &
      & Lpaired, Nc, Na, Lstate, tSSR22, tSSR44, relRhoSqr)

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

    !> Calculate DFTB/SSR(2,2) formalism
    logical, intent(in) :: tSSR22

    !> Calculate DFTB/SSR(4,4) formalism
    logical, intent(in) :: tSSR44

    !> relaxed density matrix for target L-th state
    real(dp), intent(out) :: relRhoSqr(:,:)

    real(dp), allocatable :: resRhoL(:,:)
    real(dp), allocatable :: tmpRho(:,:)
    real(dp), allocatable :: tmpMat(:,:)

    integer :: nOrb, tmpL, ii, iL, Lmax

    nOrb = size(rhoSqrL,dim=1)
    Lmax = size(rhoSqrL,dim=4)

    allocate(resRhoL(nOrb,nOrb))
    allocate(tmpRho(nOrb,nOrb))
    allocate(tmpMat(nOrb,nOrb))

    ! resRhoL is response part for L-th state
    resRhoL(:,:) = 0.0_dp

    if (tSSR22) then
      call getResponseDML22_(rhoSqrL, SAweight, ZT, omega, &
          & weightIL, G1, Lpaired, resRhoL)
    else if (tSSR44) then
      call error("SSR(4,4) is not implemented yet")
    end if

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

  end subroutine getRelaxedDML


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
      dipole(ii) = -sum(dipoleInt(:,:,ii)*(Pmat(:,:)))
    end do

  end subroutine getDipoleMomentMatrix


  !> get the oscillator strength between the states
  subroutine getReksOsc(tdp, energy)

    !> transition dipole moment between states
    real(dp), intent(in) :: tdp(:,:)

    !> energy of states
    real(dp), intent(in) :: energy(:)

    real(dp) :: tmp, osc
    integer :: ia, ib, ist, nstates, nstHalf

    nstates = size(energy,dim=1)
    nstHalf = size(tdp,dim=2)

    write(stdOut,*)
    write(stdOut,'(A)') " Oscillator Strength (au)"
    do ist = 1, nstHalf

      ! ... (ia,ib) = (1,2) (1,3) (2,3) ...
      ! TODO
      tmp = ( dble(2.0_dp*nstates+1.0_dp) - dsqrt( (2.0_dp*nstates+ &
          & 1.0_dp)**2.0_dp - 8.0_dp*(nstates+ist) ) )/2.0_dp
      ia = int( tmp )
      if( (tmp - dble(ia)) < 1.0E-8_dp ) ia = ia - 1
      ib = ia**2/2 + ia/2 - nstates*ia + nstates + ist
      if( mod(ia,2) == 1 ) ib = ib + 1

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
    real(dp), intent(out) :: tmpFilling(:,:)

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
  subroutine getUnrelaxedTDM22_(eigenvecs, FONs, Nc, nstates, P_X_del_o)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of states
    integer, intent(in) :: nstates

    !> unrelaxed transition density matrix between SA-REKS states
    real(dp), intent(out) :: P_X_del_o(:,:,:)

    real(dp) :: n_a, n_b
    integer :: mu, nu, nOrb, a, b

    nOrb = size(eigenvecs,dim=1)

    n_a = FONs(1,1)
    n_b = FONs(2,1)
    a = Nc + 1
    b = Nc + 2

    do mu = 1, nOrb
      do nu = 1, nOrb
        P_X_del_o(nu,mu,1) = eigenvecs(mu,a)*eigenvecs(nu,b) * &
            & (dsqrt(n_a) - dsqrt(n_b))
      end do
    end do
    if (nstates == 3) then
      do mu = 1, nOrb
        do nu = 1, nOrb
          P_X_del_o(nu,mu,3) = eigenvecs(mu,a)*eigenvecs(nu,b) * &
              & (dsqrt(n_a) + dsqrt(n_b))
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
    real(dp), intent(out) :: resTdm(:,:,:)

    real(dp) :: n_a, n_b
    integer :: mu, nu, nOrb, a, b, nstHalf

    nstHalf = size(resTdm,dim=3)
    nOrb = size(eigenvecs,dim=1)

    n_a = FONs(1,1); n_b = FONs(2,1)
    a = Nc + 1; b = Nc + 2

    do mu = 1, nOrb
      do nu = 1, nOrb
        resTdm(nu,mu,1) = resTdm(nu,mu,1) + SAweight(1) * &
            & ( (n_a-1.0_dp)*dsqrt(n_a) - (n_b-1.0_dp)*dsqrt(n_b) ) * &
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
              & ( (n_a-1.0_dp)*dsqrt(n_a) + (n_b-1.0_dp)*dsqrt(n_b) ) * &
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
    real(dp), intent(out) :: resRhoL(:,:)

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
