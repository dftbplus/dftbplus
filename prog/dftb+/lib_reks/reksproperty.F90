!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

! TODO
!!!!#:include 'common.fypp'

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
  use dftbp_mainio
  use dftbp_message
  use dftbp_reksvar

  implicit none

  private

  public :: getUnrelaxedDMandTDP
  public :: getDipoleIntegral, getDipoleMomentMatrix, getReksOsc

  contains

  !> Calculate unrelaxed density and transition density for target
  !> SA-REKS or SSR state (or L-th state)
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
    integer :: ii, ist, jst, kst, lst, mu, nu, ia, ib

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
        tmpRho(:,:) = P_X_o(:,:,ist) + transpose(P_X_o(:,:,ist))
        do ii = 1, nOrb
          tmpRho(ii,ii) = tmpRho(ii,ii) - P_X_o(ii,ii,ist)
        end do
        P_X_o(:,:,ist) = tmpRho(:,:)
      end do
    else
      if (self%Lstate == 0) then
        call makeDensityMatrix(P_X_o(:,:,1), eigenvecs, tmpFilling(:,self%rstate))
        tmpRho(:,:) = P_X_o(:,:,1) + transpose(P_X_o(:,:,1))
        do ii = 1, nOrb
          tmpRho(ii,ii) = tmpRho(ii,ii) - P_X_o(ii,ii,1)
        end do
        P_X_o(:,:,1) = tmpRho(:,:)
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
        P_X_o(:,:,1) = self%dm_L(:,:,1,ii)
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
      ! self%P_m_o is unrelaxed density matrix for target SSR state
      kst = 0
      self%P_m_o(:,:) = 0.0_dp
      do ist = 1, self%nstates
        do jst = ist, self%nstates
          if (ist == jst) then
            self%P_m_o(:,:) = self%P_m_o(:,:) + self%eigvecsSSR(ist,self%rstate)**2 * P_X_o(:,:,ist)
          else
            kst = kst + 1
            self%P_m_o(:,:) = self%P_m_o(:,:) + 2.0_dp * self%eigvecsSSR(ist,self%rstate) * &
                & self%eigvecsSSR(jst,self%rstate) * P_X_del_o(:,:,kst)
          end if
        end do
      end do
    else
      ! self%P_m_o is unrelaxed density matrix for target SA-REKS or L-th state
      self%P_m_o(:,:) = P_X_o(:,:,1)
    end if

    ! Final unrelaxed transition density matrix between states
    if (self%tTDP .and. self%Lstate == 0) then
      if (self%useSSR == 1) then
        ! self%P_m_del_o is unrelaxed transition density matrix between SSR states
        self%P_m_del_o(:,:,:) = 0.0_dp
        do lst = 1, nstHalf

          ! (ia,ib) = (1,2) (1,3) (2,3) ...
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
                self%P_m_del_o(:,:,lst) = self%P_m_del_o(:,:,lst) + P_X_o(:,:,ist) &
                    & * self%eigvecsSSR(ist,ia) * self%eigvecsSSR(ist,ib)
              else
                kst = kst + 1
                ! <PPS|OSS> = <OSS|PPS>, etc
                self%P_m_del_o(:,:,lst) = self%P_m_del_o(:,:,lst) + P_X_del_o(:,:,kst) &
                    & * ( self%eigvecsSSR(ist,ia) * self%eigvecsSSR(jst,ib) &
                    & + self%eigvecsSSR(jst,ia) * self%eigvecsSSR(ist,ib) )
              end if
            end do
          end do

        end do
      else
        ! self%P_m_del_o is unrelaxed transition density matrix between SA-REKS states
        self%P_m_del_o(:,:,:) = P_X_del_o(:,:,:)
      end if
    end if

    ! just calculate C^T*S*P*S*C = N, this will be diagonal.
    ! because, P = C*N*C^T, I = C^T*S*C, where
    ! P: density matrix, C: eigenvector, N: occupation number,
    ! T: transpose(real), I: identity matrix.
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, self%P_m_o, self%over)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, self%over, tmpMat)
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, tmpRho, eigenvecs)
    tmpRho(:,:) = 0.0_dp
    call gemm(tmpRho, eigenvecs, tmpMat, transA='T')

    call printUnrelaxedFONs(tmpRho, self%useSSR, self%rstate,&
        & self%Lstate, self%Nc, self%Na)

  end subroutine getUnrelaxedDMandTDP


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

    write(*,*)
    write(*,'(A)') " Oscillator Strength (au)"
    do ist = 1, nstHalf

      ! ... (ia,ib) = (1,2) (1,3) (2,3) ...
      tmp = ( dble(2.0_dp*nstates+1.0_dp) - dsqrt( (2.0_dp*nstates+ &
          & 1.0_dp)**2.0_dp - 8.0_dp*(nstates+ist) ) )/2.0_dp
      ia = int( tmp )
      if( (tmp - dble(ia)) < 1.0E-8_dp ) ia = ia - 1
      ib = ia**2/2 + ia/2 - nstates*ia + nstates + ist
      if( mod(ia,2) == 1 ) ib = ib + 1

      osc = 2.0_dp / 3.0_dp * (energy(ib) - energy(ia)) * sum(tdp(:,ist)**2)

      write(*,'(A4,I1,A6,I1,A5)',advance="no") " ( S", ia - 1, " <-> S", ib - 1, " ) : "
      write(*,'(1(f12.6))') osc

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


end module dftbp_reksproperty
