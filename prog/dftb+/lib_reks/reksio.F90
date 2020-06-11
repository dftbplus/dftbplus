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
module dftbp_reksio

  use dftbp_accuracy
  use dftbp_constants, only : au__Debye
  use dftbp_globalenv
  use dftbp_message
  use dftbp_rekscommon, only : getTwoIndices, getSpaceSym
  use dftbp_reksvar, only : TReksCalc, reksTypes

  implicit none

  private

  public :: printReksMicrostates, printSaReksEnergy
  public :: printReksSAInfo, printReksSSRInfo
  public :: printReksGradInfo
  public :: printUnrelaxedFONs, printRelaxedFONs, printRelaxedFONsL
  public :: writeReksTDP, writeReksRelaxedCharge

  contains

  !> Print energy contribution for each microstate in SCC iteration
  subroutine printReksMicrostates(self, Erep)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    !> repulsive energy
    real(dp), intent(in) :: Erep

    integer :: iL

    write(stdOut,'(1x,A,5x,A,9x,A,9x,A,9x,A,8x,A,9x,A,8x,A)') &
        & "iL", "nonSCC", "SCC", "spin", "3rd", "fock", "Rep", "Total"
    do iL = 1, self%Lmax
      write(stdOut,'(I3,7(f13.8))',advance="no") iL, self%enLnonSCC(iL), &
          & self%enLscc(iL), self%enLspin(iL)
      if (self%t3rd) then
        write(stdOut,'(1(f13.8))',advance="no") self%enL3rd(iL)
      else
        write(stdOut,'(1(f13.8))',advance="no") 0.0_dp
      end if
      if (self%isRangeSep) then
        write(stdOut,'(1(f13.8))',advance="no") self%enLfock(iL)
      else
        write(stdOut,'(1(f13.8))',advance="no") 0.0_dp
      end if
      write(stdOut,'(2(f13.8))') Erep, self%enLtot(iL)
    end do

  end subroutine printReksMicrostates


  !> Print SA-REKS energy in SCC iteration
  subroutine printSaReksEnergy(self)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    integer :: ist

    write(stdOut,'(1x,A)') "SA-REKS state energies"
    do ist = 1, self%nstates
      if (mod(ist,5) == 0 .or. ist == self%nstates) then
        write(stdOut,"(I3,':',1x,1(f13.8),1x,'H')") ist, self%energy(ist)
      else
        write(stdOut,"(I3,':',1x,1(f13.8),1x,'H')",advance="no") ist, self%energy(ist)
      end if
    end do

  end subroutine printSaReksEnergy


  !> print SA-REKS result in standard output
  subroutine printReksSAInfo(self, Etotal)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    !> state-averaged energy
    real(dp), intent(in) :: Etotal

    select case (self%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call printReksSAInfo22_(Etotal, self%enLtot, self%energy, self%FONs, self%Efunction, self%Plevel)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine printReksSAInfo


  !> print SI-SA-REKS result in standard output
  subroutine printReksSSRInfo(self, Wab, tmpEn, StateCoup)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    !> converged Lagrangian values within active space
    real(dp), intent(in) :: Wab(:,:)

    !> SA-REKS energies
    real(dp), intent(in) :: tmpEn(:)

    !> state-interaction term between SA-REKS states
    real(dp), intent(in) :: StateCoup(:,:)

    select case (self%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call printReksSSRInfo22_(Wab, tmpEn, StateCoup, self%energy, self%eigvecsSSR, &
          & self%Na, self%tAllStates, self%tSSR)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine printReksSSRInfo


  !> print gradient results for REKS calculation
  subroutine printReksGradInfo(self, derivs)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    !> derivatives of energy wrt to atomic positions
    real(dp), intent(in) :: derivs(:,:)

    integer :: ist, ia, ib, nstHalf

    nstHalf = self%nstates * (self%nstates - 1) / 2

    write(stdOut,*)
    if (self%Efunction == 1) then

      write(stdOut,"(A)") repeat("-", 50)
      write(stdOut,"(A)") " Gradient Information"
      write(stdOut,"(A)") repeat("-", 50)
      write(stdOut,*) self%rstate, "state (single-state)"
      write(stdOut,'(3(f15.8))') derivs(:,:)
      write(stdOut,"(A)") repeat("-", 50)

    else

      if (self%tNAC) then

        write(stdOut,"(A)") repeat("-", 50)
        write(stdOut,"(A)") " Gradient Information"
        write(stdOut,"(A)") repeat("-", 50)
        do ist = 1, self%nstates
          write(stdOut,*) ist, "st state (SSR)"
          write(stdOut,'(3(f15.8))') self%SSRgrad(:,:,ist)
          if (ist == self%nstates) then
            write(stdOut,"(A)") repeat("-", 50)
          else
            write(stdOut,'(3(f15.8))')
          end if
        end do

!        write(stdOut,*) "AVG state"
!        write(stdOut,'(3(f15.8))') self%avgGrad(:,:)
!        write(stdOut,'(3(f15.8))')
!        do ist = 1, self%nstates
!          write(stdOut,*) ist, "st state (SA-REKS)"
!          write(stdOut,'(3(f15.8))') self%SAgrad(:,:,ist)
!          if (ist == self%nstates) then
!            write(stdOut,"(A)") repeat("-", 50)
!          else
!            write(stdOut,'(3(f15.8))')
!          end if
!        end do

        write(stdOut,"(A)") " Coupling Information"
        do ist = 1, nstHalf

          call getTwoIndices(self%nstates, ist, ia, ib, 1)

          write(stdOut,"(A)") repeat("-", 50)
          write(stdOut,'(" between ",I2," and ",I2," states")') ia, ib
          write(stdOut,"(A)") repeat("-", 50)
          write(stdOut,*) "g vector - difference gradient"
          write(stdOut,'(3(f15.8))') (self%SAgrad(:,:,ia) - self%SAgrad(:,:,ib)) * 0.5_dp
          write(stdOut,'(3(f15.8))')
          write(stdOut,*) "h vector - derivative coupling"
          write(stdOut,'(3(f15.8))') self%SIgrad(:,:,ist)
          write(stdOut,'(3(f15.8))')
          write(stdOut,*) "G vector - GDV"
          write(stdOut,'(3(f15.8))') self%nacG(:,:,ist)
          write(stdOut,'(3(f15.8))')
          write(stdOut,*) "H vector - DCV - non-adiabatic coupling"
          write(stdOut,'(3(f15.8))') self%nacH(:,:,ist)

        end do
        write(stdOut,"(A)") repeat("-", 50)

      else

        write(stdOut,"(A)") repeat("-", 50)
        write(stdOut,"(A)") " Gradient Information"
        write(stdOut,"(A)") repeat("-", 50)
        if (self%Lstate == 0) then
          if (self%tSSR) then
            write(stdOut,*) self%rstate, "state (SSR)"
          else
            write(stdOut,*) self%rstate, "state (SA-REKS)"
          end if
        else
          write(stdOut,*) self%Lstate, "microstate"
        end if
        write(stdOut,'(3(f15.8))') derivs(:,:)
        write(stdOut,"(A)") repeat("-", 50)

      end if

    end if
    write(stdOut,*)

  end subroutine printReksGradInfo


  !> print unrelaxed FONs for target state
  subroutine printUnrelaxedFONs(tmpRho, rstate, Lstate, Nc, Na, tSSR)

    !> Occupation number matrix
    real(dp), intent(in) :: tmpRho(:,:)

    !> Target SSR state
    integer, intent(in) :: rstate

    !> Target microstate
    integer, intent(in) :: Lstate

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    logical, intent(in) :: tSSR

    integer :: ii

    write(stdOut,*)
    if (tSSR) then
      write(stdOut,'(A25,I1,A1)',advance="no") " unrelaxed SSR FONs for S", &
          & rstate - 1, ":"
    else
      if (Lstate == 0) then
        write(stdOut,'(A29,I1,A1)',advance="no") " unrelaxed SA-REKS FONs for S", &
            & rstate - 1, ":"
      else
        write(stdOut,'(A20,I1,A12)',advance="no") " unrelaxed FONs for ", &
            & Lstate, " microstate:"
      end if
    end if
    do ii = 1, Na
      if (ii == Na) then
        write(stdOut,'(1(f10.6))') tmpRho(Nc+ii,Nc+ii)
      else
        write(stdOut,'(1(f10.6))',advance="no") tmpRho(Nc+ii,Nc+ii)
      end if
    end do

  end subroutine printUnrelaxedFONs


  !> print Relaxed FONs for target state
  subroutine printRelaxedFONs(tmpRho, rstate, Nc, Na, tSSR)

    !> Occupation number matrix
    real(dp), intent(in) :: tmpRho(:,:)

    !> Target SSR state
    integer, intent(in) :: rstate

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    logical, intent(in) :: tSSR

    integer :: ii

    if (tSSR) then
      write(stdOut,'(A23,I1,A1)',advance="no") " relaxed SSR FONs for S", &
          & rstate - 1, ":"
    else
      write(stdOut,'(A27,I1,A1)',advance="no") " relaxed SA-REKS FONs for S", &
          & rstate - 1, ":"
    end if
    do ii = 1, Na
      if (ii == Na) then
        write(stdOut,'(1(f10.6))') tmpRho(Nc+ii,Nc+ii)
      else
        write(stdOut,'(1(f10.6))',advance="no") tmpRho(Nc+ii,Nc+ii)
      end if
    end do
    write(stdOut,*)

  end subroutine printRelaxedFONs


  !> print Relaxed FONs for target L-th microstate
  subroutine printRelaxedFONsL(tmpRho, Lstate, Nc, Na)

    !> Occupation number matrix
    real(dp), intent(in) :: tmpRho(:,:)

    !> Target microstate
    integer, intent(in) :: Lstate

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    integer :: ii

    write(stdOut,'(A18,I1,A12)',advance="no") " relaxed FONs for ", &
        & Lstate, " microstate:"
    do ii = 1, Na
      if (ii == Na) then
        write(stdOut,'(1(f10.6))') tmpRho(Nc+ii,Nc+ii)
      else
        write(stdOut,'(1(f10.6))',advance="no") tmpRho(Nc+ii,Nc+ii)
      end if
    end do
    write(stdOut,*)

  end subroutine printRelaxedFONsL


  !> Write tdp.dat file with transidion dipole moment
  subroutine writeReksTDP(tdp)

    real(dp), intent(in) :: tdp(:,:)

    character(len=16), parameter :: fname = "tdp.dat"
    integer :: funit

    real(dp) :: tmp
    integer :: ia, ib, ist, nstates, nstHalf

    nstHalf = size(tdp,dim=2)

    tmp = 0.5_dp * (1.0_dp + sqrt(1.0_dp + 8.0_dp*real(nstHalf,dp)))
    nstates = nint(tmp)

    open(newunit=funit,file=fname,position="rewind",status="replace")
    write(funit,*)
    do ist = 1, nstHalf

      call getTwoIndices(nstates, ist, ia, ib, 1)

      write(funit,'(A4,I1,A8,I1,A2)') " < S", ia - 1, " | r | S", ib - 1, " >"
      write(funit,'(A)',advance="no") "Transition Dipole moment (au)    : "
      write(funit,'(3(f12.6))') tdp(:,ist)
      write(funit,'(A)',advance="no") "Transition Dipole moment (Debye) : "
      write(funit,'(3(f12.6))') tdp(:,ist) * au__Debye
      write(funit,*)

    end do
    close(funit)

  end subroutine writeReksTDP


  !> Write relaxed_charge.dat file with relaxed charges for target state
  subroutine writeReksRelaxedCharge(qOutput, q0, rstate, Lstate)

    !> Output electrons
    real(dp), intent(in) :: qOutput(:,:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> Target SSR state
    integer, intent(in) :: rstate

    !> Target microstate
    integer, intent(in) :: Lstate

    character(len=20), parameter :: fname = "relaxed_charge.dat"
    integer :: iAt, nAtom
    integer :: funit

    nAtom = size(qOutput,dim=2)

    open(newunit=funit,file=fname,position="rewind",status="replace")
    write(funit,'(A13,1X,F15.8,A4)') "total charge:", &
        & -sum(qOutput(:,:,1) - q0(:,:,1)), " (e)"
    write(funit,'(1X)')
    if (Lstate == 0) then
      write(funit,'(A9,I1,A18)') "relaxed S", rstate - 1, " atomic charge (e)"
    else
      write(funit,'(I3,A11,A18)') Lstate, " microstate", " atomic charge (e)"
    end if
    write(funit,'(3X,A18)') "atom        charge"
    do iAt = 1, nAtom
      write(funit,'(2X,I5,2X,F15.8)') iAt, -sum(qOutput(:,iAt,1) - q0(:,iAt,1))
    end do
    close(funit)

  end subroutine writeReksRelaxedCharge


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> print SA-REKS(2,2) result in standard output
  subroutine printReksSAInfo22_(Etotal, enLtot, energy, FONs, Efunction, Plevel)

    !> state-averaged energy
    real(dp), intent(in) :: Etotal

    !> total energy for each microstate
    real(dp), intent(in) :: enLtot(:)

    !> energy of states
    real(dp), intent(in) :: energy(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Minimized energy functional
    integer, intent(in) :: Efunction

    !> Print level in standard output file
    integer, intent(in) :: Plevel

    real(dp) :: n_a, n_b
    integer :: iL, Lmax, ist, nstates
    character(len=8) :: strTmp

    nstates = size(energy,dim=1)
    Lmax = size(enLtot,dim=1)

    n_a = FONs(1,1)
    n_b = FONs(2,1)

    write(stdOut,*) " "
    write(stdOut, "(A)") repeat("-", 50)
    if (Efunction == 1) then
      write(stdOut,'(A25,2x,F15.8)') " Final REKS(2,2) energy:", Etotal
      write(stdOut,*) " "
      write(stdOut,'(A46)') " State     Energy      FON(1)    FON(2)   Spin"
      write(strTmp,'(A)') "PPS"
      write(stdOut,'(1x,a4,1x,f13.8,1x,2(f10.6),2x,f4.2)') &
          & trim(strTmp), energy(1), n_a, n_b, 0.0_dp
    else if (Efunction == 2) then
      write(stdOut,'(A27,2x,F15.8)') " Final SA-REKS(2,2) energy:", Etotal
      write(stdOut,*) " "
      write(stdOut,'(A46)') " State     Energy      FON(1)    FON(2)   Spin"
      do ist = 1, nstates
        if (ist == 1) then
          write(strTmp,'(A)') "PPS"
          write(stdOut,'(1x,a4,1x,f13.8,1x,2(f10.6),2x,f4.2)') &
              & trim(strTmp), energy(1), n_a, n_b, 0.0_dp
        else if (ist == 2) then
          write(strTmp,'(A)') "OSS"
          write(stdOut,'(1x,a4,1x,f13.8,1x,2(f10.6),2x,f4.2)') &
              & trim(strTmp), energy(2), 1.0_dp, 1.0_dp, 0.0_dp
        else if (ist == 3) then
          write(strTmp,'(A)') "DES"
          write(stdOut,'(1x,a4,1x,f13.8,1x,2(f10.6),2x,f4.2)') &
              & trim(strTmp), energy(3), n_b, n_a, 0.0_dp
        end if
      end do
      write(strTmp,'(A)') "Trip"
      write(stdOut,'(1x,a4,1x,f13.8,1x,2(f10.6),2x,f4.2)') &
          & trim(strTmp), enLtot(5), 1.0_dp, 1.0_dp, 1.0_dp
    end if
    write(stdOut, "(A)") repeat("-", 50)

    if (Plevel >= 2) then
      write(stdOut,*) " "
      write(stdOut, "(A)") repeat("-", 25)
      write(stdOut,'(1x,A20,2x,F15.8)') " Microstate Energies"
      do iL = 1, Lmax
        write(stdOut,"(1x,'L =',1x,I2,':',1x,F13.8)") iL, enLtot(iL)
      end do
      write(stdOut, "(A)") repeat("-", 25)
    end if

  end subroutine printReksSAInfo22_


  !> print SI-SA-REKS(2,2) result in standard output
  subroutine printReksSSRInfo22_(Wab, tmpEn, StateCoup, energy, eigvecsSSR, &
      & Na, tAllStates, tSSR)

    !> converged Lagrangian values within active space
    real(dp), intent(in) :: Wab(:,:)

    !> SA-REKS energies
    real(dp), intent(in) :: tmpEn(:)

    !> state-interaction term between SA-REKS states
    real(dp), intent(in) :: StateCoup(:,:)

    !> energy of states
    real(dp), intent(in) :: energy(:)

    !> eigenvectors from SA-REKS state
    real(dp), intent(in) :: eigvecsSSR(:,:)

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Decide the energy states in SA-REKS
    logical, intent(in) :: tAllStates

    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    logical, intent(in) :: tSSR

    integer :: ist, jst, nstates, ia, ib, nActPair
    character(len=8) :: strTmp
    character(len=1) :: stA, stB

    nActPair = size(Wab,dim=1)
    nstates = size(energy,dim=1)

    write(stdOut,*)
    do ist = 1, nActPair

      call getTwoIndices(Na, ist, ia, ib, 1)

      call getSpaceSym(ia, stA)
      call getSpaceSym(ib, stB)

      write(stdOut,"(1x,'Lagrangian W',A1,A1,': ',2(f12.8))") &
          & trim(stA), trim(stB), Wab(ist,1), Wab(ist,2)

    end do

    write(stdOut,*)
    write(stdOut, "(A)") repeat("-", 50)
    if (.not. tAllStates) then
      write(stdOut,'(A)') " SSR: 2SI-2SA-REKS(2,2) Hamiltonian matrix"
      write(stdOut,'(15x,A3,11x,A3)') "PPS", "OSS"
    else
      write(stdOut,'(A)') " SSR: 3SI-2SA-REKS(2,2) Hamiltonian matrix"
      write(stdOut,'(15x,A3,11x,A3,11x,A3)') "PPS", "OSS", "DES"
    end if

    do ist = 1, nstates
      if (ist == 1) then
        write(strTmp,'(A)') "PPS"
      else if (ist == 2) then
        write(strTmp,'(A)') "OSS"
      else if (ist == 3) then
        write(strTmp,'(A)') "DES"
      end if
      write(stdOut,'(1x,a5,1x)',advance="no") trim(strTmp)
      do jst = 1, nstates
        if (ist == jst) then
          if (jst == nstates) then
            write(stdOut,'(1x,f13.8)') tmpEn(ist)
          else
            write(stdOut,'(1x,f13.8)',advance="no") tmpEn(ist)
          end if
        else
          if (jst == nstates) then
            write(stdOut,'(1x,f13.8)') StateCoup(ist,jst)
          else
            write(stdOut,'(1x,f13.8)',advance="no") StateCoup(ist,jst)
          end if
        end if
      end do
    end do
    write(stdOut, "(A)") repeat("-", 50)

    if (tSSR) then
      write(stdOut,*)
      write(stdOut, "(A)") repeat("-", 64)
      if (.not. tAllStates) then
        write(stdOut,'(A)') " SSR: 2SI-2SA-REKS(2,2) states"
        write(stdOut,'(19x,A4,7x,A7,4x,A7)') "E_n", "C_{PPS}", "C_{OSS}"
      else
        write(stdOut,'(A)') " SSR: 3SI-2SA-REKS(2,2) states"
        write(stdOut,'(19x,A4,7x,A7,4x,A7,4x,A7)') "E_n", "C_{PPS}", "C_{OSS}", "C_{DES}"
      end if
      do ist = 1, nstates
        if (.not. tAllStates) then
          write(stdOut,'(1x,A,I2,1x,f13.8,1x,f10.6,1x,f10.6)') &
              & "SSR state ", ist, energy(ist), eigvecsSSR(:,ist)
        else
          write(stdOut,'(1x,A,I2,1x,f13.8,1x,f10.6,1x,f10.6,1x,f10.6)') &
              & "SSR state ", ist, energy(ist), eigvecsSSR(:,ist)
        end if
      end do
      write(stdOut, "(A)") repeat("-", 64)
    end if

  end subroutine printReksSSRInfo22_


!  !> Calculate the weight of each microstate for current cycle, C_L
!  subroutine calcWeights(self)
!
!    !> data type for REKS
!    type(TReksCalc), intent(inout) :: self
!
!    select case (self%reksAlg)
!    case (reksTypes%noReks)
!    case (reksTypes%ssr22)
!      call getWeightL22_(self%FONs, self%delta, self%SAweight, self%weightL, self%weight)
!    case (reksTypes%ssr44)
!      call error("SSR(4,4) is not implemented yet")
!    end select
!
!  end subroutine calcWeights
!
!
!  !> Swap the active orbitals for feasible occupation in REKS
!  subroutine activeOrbSwap(self, eigenvecs)
!
!    !> data type for REKS
!    type(TReksCalc), intent(inout) :: self
!
!    !> eigenvectors
!    real(dp), intent(inout) :: eigenvecs(:,:)
!
!    select case (self%reksAlg)
!    case (reksTypes%noReks)
!    case (reksTypes%ssr22)
!      call MOswap22_(eigenvecs, self%SAweight, self%FONs, self%Efunction, self%Nc)
!    case (reksTypes%ssr44)
!      call error("SSR(4,4) is not implemented yet")
!    end select
!
!  end subroutine activeOrbSwap
!
!
!  !> Calculate filling for minimzed state with optimized FONs
!  subroutine getFilling(self, filling)
!
!    !> data type for REKS
!    type(TReksCalc), intent(inout) :: self
!
!    !> occupations (level)
!    real(dp), intent(out) :: filling(:)
!
!    select case (self%reksAlg)
!    case (reksTypes%noReks)
!    case (reksTypes%ssr22)
!      call getFilling22_(filling, self%SAweight, self%FONs, self%Efunction, self%Nc)
!    case (reksTypes%ssr44)
!      call error("SSR(4,4) is not implemented yet")
!    end select
!
!  end subroutine getFilling
!
!
!  !> Calculate the energy of SA-REKS states and averaged state
!  subroutine calcSaReksEnergy(self, energy)
!
!    !> data type for REKS
!    type(TReksCalc), intent(inout) :: self
!
!    !> Energy terms in the system
!    type(TEnergies), intent(inout) :: energy
!
!    integer :: ist
!
!    ! Compute the energy contributions for target SA-REKS state
!    ! electronic energy = nonSCC + scc + spin + 3rd + fock
!    energy%EnonSCC = sum(self%weightL(self%rstate,:)*self%enLnonSCC(:))
!    energy%Escc = sum(self%weightL(self%rstate,:)*self%enLscc(:))
!    energy%Espin = sum(self%weightL(self%rstate,:)*self%enLspin(:))
!    if (self%t3rd) then
!      energy%e3rd = sum(self%weightL(self%rstate,:)*self%enL3rd(:))
!    end if
!    if (self%isRangeSep) then
!      energy%Efock = sum(self%weightL(self%rstate,:)*self%enLfock(:))
!    end if
!
!    energy%Eelec = energy%EnonSCC + energy%Escc + energy%Espin + &
!        & energy%e3rd + energy%Efock
!
!    ! Compute the total energy for SA-REKS states
!    do ist = 1, self%nstates
!      self%energy(ist) = sum(self%weightL(ist,:)*self%enLtot(:))
!    end do
!
!    ! In this step Etotal becomes the energy of averaged state, not individual states
!    ! From this energy we can check the variational principle
!    energy%Etotal = 0.0_dp
!    do ist = 1, self%SAstates
!      energy%Etotal = energy%Etotal + self%SAweight(ist) * self%energy(ist)
!    end do
!
!  end subroutine calcSaReksEnergy
!
!
!  !> Make pseudo-fock operator with Hamiltonian of each microstate
!  !> and diagonalize the fock matrix
!  subroutine getFockandDiag(env, denseDesc, neighbourList, &
!      & nNeighbourSK, iSparseStart, img2CentCell, eigenvecs, &
!      & electronicSolver, eigen, self)
!
!    !> Environment settings
!    type(TEnvironment), intent(inout) :: env
!
!    !> Dense matrix descriptor
!    type(TDenseDescr), intent(in) :: denseDesc
!
!    !> neighbours to atoms
!    type(TNeighbourList), intent(in) :: neighbourList
!
!    !> Number of atomic neighbours
!    integer, intent(in) :: nNeighbourSK(:)
!
!    !> Index for atomic blocks in sparse data
!    integer, intent(in) :: iSparseStart(:,:)
!
!    !> map from image atom to real atoms
!    integer, intent(in) :: img2CentCell(:)
!
!    !> Eigenvectors on eixt
!    real(dp), intent(in) :: eigenvecs(:,:,:)
!
!    !> Electronic solver information
!    type(TElectronicSolver), intent(inout) :: electronicSolver
!
!    !> eigenvalues
!    real(dp), intent(out) :: eigen(:,:,:)
!
!    !> data type for REKS
!    type(TReksCalc), intent(inout) :: self
!
!    real(dp), allocatable :: orbFON(:)
!    real(dp), allocatable :: tmpOver(:,:)
!    real(dp), allocatable :: tmpMat(:,:)
!
!    integer :: ii, nOrb
!
!    nOrb = size(self%fockFc,dim=1)
!
!    allocate(orbFON(nOrb))
!    allocate(tmpOver(nOrb,nOrb))
!    allocate(tmpMat(nOrb,nOrb))
!
!    call getFockFcFa_(env, denseDesc, neighbourList, nNeighbourSK, &
!        & iSparseStart, img2CentCell, self%hamSqrL, self%hamSpL, self%weight, &
!        & self%fillingL, self%Nc, self%Na, self%Lpaired, self%isRangeSep, &
!        & orbFON, self%fockFc, self%fockFa)
!
!    call matAO2MO(self%fockFc, eigenvecs(:,:,1))
!    do ii = 1, self%Na
!      call matAO2MO(self%fockFa(:,:,ii), eigenvecs(:,:,1))
!    end do
!
!    call getPseudoFock_(self%fockFc, self%fockFa, orbFON, self%Nc, self%Na, self%fock)
!
!    call levelShifting_(self%fock, self%shift, self%Nc, self%Na)
!
!    ! Diagonalize the pesudo-Fock matrix
!    tmpOver(:,:) = 0.0_dp
!    do ii = 1, nOrb
!      tmpOver(ii,ii) = 1.0_dp
!    end do
!    tmpMat(:,:) = self%fock
!
!    eigen(:,1,1) = 0.0_dp
!    call diagDenseMtx(electronicSolver, 'V', tmpMat, tmpOver, eigen(:,1,1))
!    self%eigvecsFock(:,:) = tmpMat
!
!  end subroutine getFockandDiag
!
!
!  !> guess new eigenvectors from Fock eigenvectors
!  subroutine guessNewEigvecs(eigenvecs, eigvecsFock)
!
!    !> Eigenvectors on eixt
!    real(dp), intent(inout) :: eigenvecs(:,:)
!
!    !> eigenvectors from pesudo-fock matrix
!    real(dp), intent(in) :: eigvecsFock(:,:)
!
!    real(dp), allocatable :: tmpVec(:,:)
!    integer :: nOrb
!
!    nOrb = size(eigvecsFock,dim=1)
!
!    allocate(tmpVec(nOrb,nOrb))
!
!    tmpVec(:,:) = 0.0_dp
!    call gemm(tmpVec, eigenvecs, eigvecsFock)
!    eigenvecs(:,:) = tmpVec
!
!  end subroutine guessNewEigvecs
!
!
!  !> adjust the eigenvalues (eliminate shift values)
!  subroutine adjustEigenval(self, eigen)
!
!    !> data type for REKS
!    type(TReksCalc), intent(inout) :: self
!
!    !> eigenvalues
!    real(dp), intent(inout) :: eigen(:,:,:)
!
!    integer :: nOrb, ind, ii
!
!    nOrb = size(eigen,dim=1)
!
!    do ii = self%Nc + 1, self%Nc + self%Na
!      ind = ii - self%Nc
!      eigen(ii,1,1) = eigen(ii,1,1) - real(ind, dp) * self%shift
!    end do
!
!    do ii = self%Nc + self%Na + 1, nOrb
!      ind = self%Na + 1
!      eigen(ii,1,1) = eigen(ii,1,1) - real(ind, dp) * self%shift
!    end do
!
!  end subroutine adjustEigenval
!
!
!  !> Solve secular equation with coupling element between SA-REKS states
!  subroutine solveSecularEqn(env, denseDesc, neighbourList, &
!      & nNeighbourSK, iSparseStart, img2CentCell, electronicSolver, &
!      & eigenvecs, self)
!
!    !> Environment settings
!    type(TEnvironment), intent(inout) :: env
!
!    !> Dense matrix descriptor
!    type(TDenseDescr), intent(in) :: denseDesc
!
!    !> neighbours to atoms
!    type(TNeighbourList), intent(in) :: neighbourList
!
!    !> Number of atomic neighbours
!    integer, intent(in) :: nNeighbourSK(:)
!
!    !> Index for atomic blocks in sparse data
!    integer, intent(in) :: iSparseStart(:,:)
!
!    !> map from image atom to real atoms
!    integer, intent(in) :: img2CentCell(:)
!
!    !> Electronic solver information
!    type(TElectronicSolver), intent(inout) :: electronicSolver
!
!    !> Eigenvectors on eixt
!    real(dp), intent(in) :: eigenvecs(:,:,:)
!
!    !> data type for REKS
!    type(TReksCalc), intent(inout) :: self
!
!    real(dp), allocatable :: Wab(:,:)
!    real(dp), allocatable :: StateCoup(:,:)
!    real(dp), allocatable :: tmpOver(:,:)
!    real(dp), allocatable :: tmpState(:,:)
!    real(dp), allocatable :: tmpEigen(:)
!    real(dp), allocatable :: tmpEn(:)
!
!    integer :: ist, jst, nActPair
!
!    nActPair = self%Na * (self%Na - 1) / 2
!
!    allocate(Wab(nActPair,2))
!    allocate(StateCoup(self%nstates,self%nstates))
!    allocate(tmpOver(self%nstates,self%nstates))
!    allocate(tmpState(self%nstates,self%nstates))
!    allocate(tmpEigen(self%nstates))
!    allocate(tmpEn(self%nstates))
!
!    call getLagrangians_(env, denseDesc, neighbourList, nNeighbourSK, &
!        & iSparseStart, img2CentCell, eigenvecs(:,:,1), self%hamSqrL, &
!        & self%hamSpL, self%weight, self%fillingL, self%Nc, self%Na, &
!        & self%Lpaired, self%isRangeSep, Wab)
!
!    select case (self%reksAlg)
!    case (reksTypes%noReks)
!    case (reksTypes%ssr22)
!      call getStateCoup22_(Wab, self%FONs, StateCoup)
!    case (reksTypes%ssr44)
!      call error("SSR(4,4) is not implemented yet")
!    end select
!
!    ! diagonalize the state energies
!    ! obtain SSR energies & state-interaction term
!    tmpOver(:,:) = 0.0_dp
!    do ist = 1, self%nstates
!      tmpOver(ist,ist) = 1.0_dp
!    end do
!    tmpEigen(:) = 0.0_dp
!
!    tmpState(:,:) = 0.0_dp
!    do ist = 1, self%nstates
!      do jst = 1, self%nstates
!        if (ist == jst) then
!          tmpState(ist,jst) = self%energy(ist)
!        else
!          tmpState(ist,jst) = StateCoup(ist,jst)
!        end if
!      end do
!    end do
!
!    ! save state energies to print information
!    tmpEn(:) = self%energy
!    if (self%tSSR) then
!      call diagDenseMtx(electronicSolver, 'V', tmpState, tmpOver, tmpEigen)
!      self%eigvecsSSR(:,:) = tmpState
!      self%energy(:) = tmpEigen
!    end if
!
!    ! print state energies and couplings
!    call printReksSSRInfo(self, Wab, tmpEn, StateCoup)
!
!  end subroutine solveSecularEqn
!  !> Calculate filling of each microstate in REKS(2,2)
!  subroutine getFillingL22_(Nc, fillingL)
!    
!    !> Number of core orbitals
!    integer, intent(in) :: Nc
!
!    !> Filling for each microstate
!    real(dp), intent(out) :: fillingL(:,:,:)
!
!    integer :: iL, iSpin, ii, nSpin, Lmax
!
!    nSpin = size(fillingL,dim=2)
!    Lmax = size(fillingL,dim=3)
!
!    fillingL(:,:,:) = 0.0_dp
!
!    ! Filling of core orbitals
!    do iL = 1, Lmax
!      do iSpin = 1, nSpin
!        do ii = 1, Nc
!          fillingL(ii,iSpin,iL) = 1.0_dp
!        end do
!      end do
!    end do
!
!    ! Filling of active orbitals for REKS(2,2) case
!    ! 1 = a: up + down
!    fillingL(Nc+1,1,1) = 1.0_dp
!    fillingL(Nc+1,2,1) = 1.0_dp
!    ! 2 = b: up + down
!    fillingL(Nc+2,1,2) = 1.0_dp
!    fillingL(Nc+2,2,2) = 1.0_dp
!    ! 3 = a: up, b: down
!    fillingL(Nc+1,1,3) = 1.0_dp
!    fillingL(Nc+2,2,3) = 1.0_dp
!    ! 4 = a: down, b: up
!    fillingL(Nc+2,1,4) = 1.0_dp
!    fillingL(Nc+1,2,4) = 1.0_dp
!    ! 5 = a: up, b: up
!    fillingL(Nc+1,1,5) = 1.0_dp
!    fillingL(Nc+2,1,5) = 1.0_dp
!    ! 6 = a: down, b: down
!    fillingL(Nc+1,2,6) = 1.0_dp
!    fillingL(Nc+2,2,6) = 1.0_dp
!
!  end subroutine getFillingL22_
!
!
!  !> Make (2e,2o) weights, C_L used in SA-REKS
!  subroutine getWeightL22_(FONs, delta, SAweight, weightL, weight)
!
!    !> Fractional occupation numbers of active orbitals
!    real(dp), intent(in) :: FONs(:,:)
!
!    !> Smoothing factor used in FON optimization
!    real(dp), intent(in) :: delta
!
!    !> Weights used in state-averaging
!    real(dp), intent(in) :: SAweight(:)
!
!    !> Weight for each microstate per state
!    real(dp), intent(out) :: weightL(:,:)
!
!    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
!    real(dp), intent(out) :: weight(:)
!
!    integer :: iL, Lmax, ist, SAstates, nstates
!    real(dp) :: n_a, n_b, fac
!
!    Lmax = size(weightL,dim=2)
!    SAstates = size(SAweight,dim=1)
!    nstates = size(weightL,dim=1)
!
!    n_a = FONs(1,1)
!    n_b = FONs(2,1)
!
!    fac = getFactor(n_a, n_b, delta)
!
!    weightL(1,1) = 0.5_dp*n_a
!    weightL(1,2) = 0.5_dp*n_b
!    weightL(1,3) = fac
!    weightL(1,4) = fac
!    weightL(1,5) = -fac
!    weightL(1,6) = -fac
!
!    if(nstates >= 2) then
!      weightL(2,1) = 0.0_dp
!      weightL(2,2) = 0.0_dp
!      weightL(2,3) = 1.0_dp
!      weightL(2,4) = 1.0_dp
!      weightL(2,5) = -0.5_dp
!      weightL(2,6) = -0.5_dp
!    end if
!
!    if(nstates >= 3) then
!      weightL(3,1) = 0.5_dp*n_b
!      weightL(3,2) = 0.5_dp*n_a
!      weightL(3,3) = -fac
!      weightL(3,4) = -fac
!      weightL(3,5) = fac
!      weightL(3,6) = fac
!    end if
!
!    ! Decide which state will be optimized
!    ! SAstates = 1 -> PPS state is optimized
!    ! SAstates = 2 -> (PPS+OSS)/2 state (averaged state) is optimized
!    weight(:) = 0.0_dp
!    do iL = 1, Lmax
!      do ist = 1, SAstates
!        weight(iL) = weight(iL) + SAweight(ist)*weightL(ist,iL)
!      end do
!    end do
!
!  end subroutine getWeightL22_
!
!
!  !> Swap active orbitals when fa < fb in REKS(2,2) case
!  subroutine MOswap22_(eigenvecs, SAweight, FONs, Efunction, Nc)
!
!    !> eigenvectors
!    real(dp), intent(inout) :: eigenvecs(:,:)
!
!    !> Weights used in state-averaging
!    real(dp), intent(in) :: SAweight(:)
!
!    !> Fractional occupation numbers of active orbitals
!    real(dp), intent(in) :: FONs(:,:)
!
!    !> Minimized energy functional
!    integer, intent(in) :: Efunction
!
!    !> Number of core orbitals
!    integer, intent(in) :: Nc
!
!    real(dp), allocatable :: tmpMO(:)
!
!    real(dp) :: n_a, n_b, fa, fb
!    integer :: nOrb
!
!    nOrb = size(eigenvecs,dim=1)
!
!    n_a = FONs(1,1)
!    n_b = FONs(2,1)
!
!    allocate(tmpMO(nOrb))
!
!    if (Efunction == 1) then
!      ! REKS charge
!      fa = n_a * 0.5_dp
!      fb = n_b * 0.5_dp
!    else if (Efunction == 2) then
!      ! SA-REKS charge
!      fa = (n_a*SAweight(1) + SAweight(2)) * 0.5_dp
!      fb = (n_b*SAweight(1) + SAweight(2)) * 0.5_dp
!    end if
!
!    if (fa < fb) then
!      write(stdOut,'(A6,F9.6,A20,I4,A8,I4,A8)') " fa = ", fa, &
!          & ", MO swap between a(", Nc+1, ") and b(", Nc+2, ") occurs"
!      tmpMO(:) = eigenvecs(:,Nc+1)
!      eigenvecs(:,Nc+1) = eigenvecs(:,Nc+2)
!      eigenvecs(:,Nc+2) = tmpMO
!    end if
!
!  end subroutine MOswap22_
!
!
!  !> Calculate filling for minimzed state with optimized FONs in REKS(2,2)
!  subroutine getFilling22_(filling, SAweight, FONs, Efunction, Nc)
!
!    !> occupations (level)
!    real(dp), intent(out) :: filling(:)
!
!    !> Weights used in state-averaging
!    real(dp), intent(in) :: SAweight(:)
!
!    !> Fractional occupation numbers of active orbitals
!    real(dp), intent(in) :: FONs(:,:)
!
!    !> Minimized energy functional
!    integer, intent(in) :: Efunction
!
!    !> Number of core orbitals
!    integer, intent(in) :: Nc
!
!    real(dp) :: n_a, n_b
!    integer :: ii
!
!    n_a = FONs(1,1)
!    n_b = FONs(2,1)
!
!    filling(:) = 0.0_dp
!    do ii = 1, Nc
!      filling(ii) = 2.0_dp
!    end do
!    if (Efunction == 1) then
!      ! REKS charge
!      filling(Nc+1) = n_a
!      filling(Nc+2) = n_b
!    else if (Efunction == 2) then
!      ! SA-REKS charge
!      filling(Nc+1) = n_a*SAweight(1) + SAweight(2)
!      filling(Nc+2) = n_b*SAweight(1) + SAweight(2)
!    end if
!
!  end subroutine getFilling22_
!
!
!  !> Calculate Fc and Fa from Hamiltonian of each microstate
!  subroutine getFockFcFa_(env, denseDesc, neighbourList, nNeighbourSK, &
!      & iSparseStart, img2CentCell, hamSqrL, hamSpL, weight, fillingL, &
!      & Nc, Na, Lpaired, isRangeSep, orbFON, Fc, Fa)
!
!    !> Environment settings
!    type(TEnvironment), intent(inout) :: env
!
!    !> Dense matrix descriptor
!    type(TDenseDescr), intent(in) :: denseDesc
!
!    !> neighbours to atoms
!    type(TNeighbourList), intent(in) :: neighbourList
!
!    !> Number of atomic neighbours
!    integer, intent(in) :: nNeighbourSK(:)
!
!    !> Index for atomic blocks in sparse data
!    integer, intent(in) :: iSparseStart(:,:)
!
!    !> map from image atom to real atoms
!    integer, intent(in) :: img2CentCell(:)
!
!    !> state-averaged occupation numbers
!    real(dp), intent(inout) :: orbFON(:)
!
!    !> dense fock matrix for core orbitals
!    real(dp), intent(out) :: Fc(:,:)
!
!    !> dense fock matrix for active orbitals
!    real(dp), intent(out) :: Fa(:,:,:)
!
!    !> Dense Hamiltonian matrix for each microstate
!    real(dp), allocatable, intent(inout) :: hamSqrL(:,:,:,:)
!
!    !> Sparse Hamiltonian matrix for each microstate
!    real(dp), allocatable, intent(in) :: hamSpL(:,:,:)
!
!    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
!    real(dp), intent(in) :: weight(:)
!
!    !> Filling for each microstate
!    real(dp), intent(in) :: fillingL(:,:,:)
!
!    !> Number of core orbitals
!    integer, intent(in) :: Nc
!
!    !> Number of active orbitals
!    integer, intent(in) :: Na
!
!    !> Number of spin-paired microstates
!    integer, intent(in) :: Lpaired
!
!    !> Whether to run a range separated calculation
!    logical, intent(in) :: isRangeSep
!
!    real(dp), allocatable :: tmpHam(:,:)
!
!    integer :: iL, Lmax, nOrb
!
!    nOrb = size(Fc,dim=1)
!    Lmax = size(weight,dim=1)
!
!    if (.not. isRangeSep) then
!      allocate(tmpHam(nOrb,nOrb))
!    end if
!
!    call fockFON_(fillingL, weight, orbFON)
!
!    Fc(:,:) = 0.0_dp
!    Fa(:,:,:) = 0.0_dp
!    do iL = 1, Lmax
!
!      if (.not. isRangeSep) then
!        tmpHam(:,:) = 0.0_dp
!        ! convert from sparse to dense for hamSpL in AO basis
!        ! hamSpL has (my_ud) component
!        call env%globalTimer%startTimer(globalTimers%sparseToDense)
!        call unpackHS(tmpHam, hamSpL(:,1,iL), neighbourList%iNeighbour, nNeighbourSK, &
!            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
!        call env%globalTimer%stopTimer(globalTimers%sparseToDense)
!        call blockSymmetrizeHS(tmpHam, denseDesc%iAtomStart)
!      end if
!
!      ! compute the Fock operator with core, a, b orbitals in AO basis
!      if (isRangeSep) then
!        call fockFcAO_(hamSqrL(:,:,1,iL), weight, Lpaired, iL, Fc)
!        call fockFaAO_(hamSqrL(:,:,1,iL), weight, fillingL, orbFON, &
!            & Nc, Na, Lpaired, iL, Fa)
!      else
!        call fockFcAO_(tmpHam, weight, Lpaired, iL, Fc)
!        call fockFaAO_(tmpHam, weight, fillingL, orbFON, &
!            & Nc, Na, Lpaired, iL, Fa)
!      end if
!
!    end do
!
!  end subroutine getFockFcFa_
!
!
!  !> Calculate pseudo-fock matrix from Fc and Fa
!  subroutine getPseudoFock_(Fc, Fa, orbFON, Nc, Na, fock)
!
!    !> dense pseudo-fock matrix
!    real(dp), intent(out) :: fock(:,:)
!
!    !> dense fock matrix for core orbitals
!    real(dp), intent(in) :: Fc(:,:)
!
!    !> dense fock matrix for active orbitals
!    real(dp), intent(in) :: Fa(:,:,:)
!
!    !> state-averaged occupation numbers
!    real(dp), intent(in) :: orbFON(:)
!
!    !> Number of core orbitals
!    integer, intent(in) :: Nc
!
!    !> Number of active orbitals
!    integer, intent(in) :: Na
!
!    real(dp) :: res
!    integer :: ii, jj, ind1, ind2, nOrb
!
!    nOrb = size(fock,dim=1)
!
!    fock(:,:) = 0.0_dp
!    do ii = 1, Nc
!      do jj = ii, Nc
!        fock(jj,ii) = Fc(ii,jj)
!      end do
!      do jj = Nc + 1, Nc + Na
!        ind1 = jj - Nc
!        call fockFijMO_(res, Fc(ii,jj), Fa(ii,jj,ind1), &
!            & orbFON(ii), orbFON(jj))
!        fock(jj,ii) = res
!      end do
!      do jj = Nc + Na + 1, nOrb
!        fock(jj,ii) = Fc(ii,jj)
!      end do
!    end do
!
!    do jj = Nc + Na + 1, nOrb
!      do ii = Nc + 1, Nc + Na
!        ind1 = ii - Nc
!        call fockFijMO_(res, Fc(jj,ii), Fa(jj,ii,ind1), &
!            & orbFON(jj), orbFON(ii))
!        fock(jj,ii) = res
!      end do
!      do ii = jj, nOrb
!        fock(ii,jj) = Fc(jj,ii)
!      end do
!    end do
!
!    do ii = Nc + 1, Nc + Na
!      ind1 = ii - Nc
!      do jj = Nc + 1, Nc + Na
!        ind2 = jj - Nc
!        if (ii == jj) then
!          fock(jj,ii) = Fa(ii,jj,ind1)
!        else
!          call fockFijMO_(res, Fa(ii,jj,ind1), Fa(ii,jj,ind2), &
!              & orbFON(ii), orbFON(jj))
!          fock(jj,ii) = res
!        end if
!      end do
!    end do
!
!    call symmetrizeHS(fock)
!
!  end subroutine getPseudoFock_
!
!
!  !> Avoid changing the order of MOs
!  !> Required number of cycles increases as the number of shift increases
!  subroutine levelShifting_(fock, shift, Nc, Na)
!
!    !> dense pseudo-fock matrix
!    real(dp), intent(inout) :: fock(:,:)
!
!    !> Shift value in SCC cycle
!    real(dp), intent(in) :: shift
!
!    !> Number of core orbitals
!    integer, intent(in) :: Nc
!
!    !> Number of active orbitals
!    integer, intent(in) :: Na
!
!    integer :: nOrb, ind, ii
!
!    nOrb = size(fock,dim=1)
!
!    do ii = Nc + 1, Nc + Na
!      ind = ii - Nc
!      fock(ii,ii) = fock(ii,ii) + real(ind, dp) * shift
!    end do
!
!    do ii = Nc + Na + 1, nOrb
!      ind = Na + 1
!      fock(ii,ii) = fock(ii,ii) + real(ind, dp) * shift
!    end do
!
!  end subroutine levelShifting_
!
!
!  !> Calculate state-averaged FONs
!  subroutine fockFON_(fillingL, weight, orbFON)
!
!    !> state-averaged occupation numbers
!    real(dp), intent(out) :: orbFON(:)
!
!    !> Filling for each microstate
!    real(dp), intent(in) :: fillingL(:,:,:)
!
!    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
!    real(dp), intent(in) :: weight(:)
!
!    integer :: Lmax, iL
!
!    Lmax = size(weight,dim=1)
!
!    orbFON(:) = 0.0_dp
!    do iL = 1, Lmax
!      orbFON(:) = orbFON(:) + 0.5_dp * weight(iL) * &
!          & ( fillingL(:,1,iL) + fillingL(:,2,iL) )
!    end do
!
!  end subroutine fockFON_
!
!
!  !> Calculate fock matrix for core orbitals in AO basis
!  subroutine fockFcAO_(hamSqr, weight, Lpaired, iL, Fc)
!
!    !> dense fock matrix for core orbitals
!    real(dp), intent(inout) :: Fc(:,:)
!
!    !> Dense Hamiltonian matrix for each microstate
!    real(dp), intent(in) :: hamSqr(:,:)
!
!    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
!    real(dp), intent(in) :: weight(:)
!
!    !> Number of spin-paired microstates
!    integer, intent(in) :: Lpaired
!
!    !> current index in loop L
!    integer, intent(in) :: iL
!
!    if (iL <= Lpaired) then
!      Fc(:,:) = Fc + 0.5_dp * hamSqr * &
!          & ( weight(iL) + weight(iL) )
!    else
!      if (mod(iL,2) == 1) then
!        Fc(:,:) = Fc + 0.5_dp * hamSqr * &
!            & ( weight(iL) + weight(iL+1) )
!      else
!        Fc(:,:) = Fc + 0.5_dp * hamSqr * &
!            & ( weight(iL) + weight(iL-1) )
!      end if
!    end if
!
!  end subroutine fockFcAO_
!
!
!  !> Calculate fock matrix for active orbitals in AO basis
!  subroutine fockFaAO_(hamSqr, weight, fillingL, orbFON, Nc, Na, &
!      & Lpaired, iL, Fa)
!
!    !> dense fock matrix for active orbitals
!    real(dp), intent(inout) :: Fa(:,:,:)
!
!    !> Dense Hamiltonian matrix for each microstate
!    real(dp), intent(in) :: hamSqr(:,:)
!
!    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
!    real(dp), intent(in) :: weight(:)
!
!    !> Filling for each microstate
!    real(dp), intent(in) :: fillingL(:,:,:)
!
!    !> state-averaged occupation numbers
!    real(dp), intent(in) :: orbFON(:)
!
!    !> Number of core orbitals
!    integer, intent(in) :: Nc
!
!    !> Number of active orbitals
!    integer, intent(in) :: Na
!
!    !> Number of spin-paired microstates
!    integer, intent(in) :: Lpaired
!
!    !> current index in loop L
!    integer, intent(in) :: iL
!
!    integer :: ind, ind_a
!
!    do ind = 1, Na
!      ind_a = Nc + ind
!      if (iL <= Lpaired) then
!        Fa(:,:,ind) = Fa(:,:,ind) + 0.5_dp * fillingL(ind_a,1,iL) * &
!            & ( weight(iL) + weight(iL) ) * hamSqr / orbFON(ind_a)
!      else
!        if (mod(iL,2) == 1) then
!          Fa(:,:,ind) = Fa(:,:,ind) + 0.5_dp * fillingL(ind_a,1,iL) * &
!              & ( weight(iL) + weight(iL+1) ) * hamSqr / orbFON(ind_a)
!        else
!          Fa(:,:,ind) = Fa(:,:,ind) + 0.5_dp * fillingL(ind_a,1,iL) * &
!              & ( weight(iL) + weight(iL-1) ) * hamSqr / orbFON(ind_a)
!        end if
!      end if
!    end do
!
!  end subroutine fockFaAO_
!
!
!  !> Calculate pseudo-fock off-diagonal element in MO basis
!  subroutine fockFijMO_(res, fock_i, fock_j, f_i, f_j)
!
!    !> temporary pseudo-fock value
!    real(dp), intent(out) :: res
!
!    !> temporary Fc or Fa values
!    real(dp), intent(in) :: fock_i, fock_j
!
!    !> temporary orbFON values
!    real(dp), intent(in) :: f_i, f_j
!
!    real(dp) :: eps = 1.0E-3_dp
!
!    res = 0.0_dp
!    if (abs(f_j-f_i) .LT. eps) then
!      if (f_j >= f_i) then
!        res = ( f_j*fock_j - f_i*fock_i )
!      else
!        res = -( f_j*fock_j - f_i*fock_i )
!      end if
!    else
!      res = ( f_j*fock_j - f_i*fock_i ) / (f_j - f_i)
!    end if
!
!  end subroutine fockFijMO_
!
!
!  !> Calculate converged Lagrangian values
!  subroutine getLagrangians_(env, denseDesc, neighbourList, nNeighbourSK, &
!      & iSparseStart, img2CentCell, eigenvecs, hamSqrL, hamSpL, weight, &
!      & fillingL, Nc, Na, Lpaired, isRangeSep, Wab)
!
!    !> Environment settings
!    type(TEnvironment), intent(inout) :: env
!
!    !> Dense matrix descriptor
!    type(TDenseDescr), intent(in) :: denseDesc
!
!    !> neighbours to atoms
!    type(TNeighbourList), intent(in) :: neighbourList
!
!    !> Number of atomic neighbours
!    integer, intent(in) :: nNeighbourSK(:)
!
!    !> Index for atomic blocks in sparse data
!    integer, intent(in) :: iSparseStart(:,:)
!
!    !> map from image atom to real atoms
!    integer, intent(in) :: img2CentCell(:)
!
!    !> Eigenvectors on eixt
!    real(dp), intent(in) :: eigenvecs(:,:)
!
!    !> Dense Hamiltonian matrix for each microstate
!    real(dp), allocatable, intent(inout) :: hamSqrL(:,:,:,:)
!
!    !> Sparse Hamiltonian matrix for each microstate
!    real(dp), allocatable, intent(in) :: hamSpL(:,:,:)
!
!    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
!    real(dp), intent(in) :: weight(:)
!
!    !> Filling for each microstate
!    real(dp), intent(in) :: fillingL(:,:,:)
!
!    !> Number of core orbitals
!    integer, intent(in) :: Nc
!
!    !> Number of active orbitals
!    integer, intent(in) :: Na
!
!    !> Number of spin-paired microstates
!    integer, intent(in) :: Lpaired
!
!    !> Whether to run a range separated calculation
!    logical, intent(in) :: isRangeSep
!
!    !> converged Lagrangian values within active space
!    real(dp), intent(out) :: Wab(:,:)
!
!    real(dp), allocatable :: tmpHam(:,:)
!    real(dp), allocatable :: tmpHamL(:,:,:)
!
!    integer :: nOrb, iL, Lmax
!    integer :: ia, ib, ist, nActPair
!
!    nOrb = size(eigenvecs,dim=1)
!    Lmax = size(fillingL,dim=3)
!    nActPair = Na * (Na - 1) / 2
!
!    if (.not. isRangeSep) then
!      allocate(tmpHam(nOrb,nOrb))
!    end if
!    allocate(tmpHamL(nActPair,1,Lmax))
!
!    tmpHamL(:,:,:) = 0.0_dp
!    do ist = 1, nActPair
!
!      call getTwoIndices(Na, ist, ia, ib, 1)
!
!      do iL = 1, Lmax
!
!        if (isRangeSep) then
!          ! convert hamSqrL from AO basis to MO basis
!          ! hamSqrL has (my_ud) component
!          if (ist == 1) then
!            call matAO2MO(hamSqrL(:,:,1,iL), eigenvecs)
!          end if
!          tmpHamL(ist,1,iL) = hamSqrL(Nc+ia,Nc+ib,1,iL)
!        else
!          tmpHam(:,:) = 0.0_dp
!          ! convert from sparse to dense for hamSpL in AO basis
!          ! hamSpL has (my_ud) component
!          call env%globalTimer%startTimer(globalTimers%sparseToDense)
!          call unpackHS(tmpHam, hamSpL(:,1,iL), &
!              & neighbourList%iNeighbour, nNeighbourSK, &
!              & denseDesc%iAtomStart, iSparseStart, img2CentCell)
!          call env%globalTimer%stopTimer(globalTimers%sparseToDense)
!          call blockSymmetrizeHS(tmpHam, denseDesc%iAtomStart)
!          ! convert tmpHam from AO basis to MO basis
!          call matAO2MO(tmpHam, eigenvecs)
!          ! save F_{L,ab}^{\sigma} in MO basis
!          tmpHamL(ist,1,iL) = tmpHam(Nc+ia,Nc+ib)
!        end if
!
!      end do
!
!      ! calculate the Lagrangian eps_{ab} and state-interaction term
!      Wab(ist,1) = 0.0_dp
!      Wab(ist,2) = 0.0_dp
!      do iL = 1, Lmax
!        if (iL <= Lpaired) then
!          Wab(ist,1) = Wab(ist,1) + fillingL(Nc+ia,1,iL) * &
!              & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL) )
!          Wab(ist,2) = Wab(ist,2) + fillingL(Nc+ib,1,iL) * &
!              & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL) )
!        else
!          if (mod(iL,2) == 1) then
!            Wab(ist,1) = Wab(ist,1) + fillingL(Nc+ia,1,iL) * &
!                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL+1) )
!            Wab(ist,2) = Wab(ist,2) + fillingL(Nc+ib,1,iL) * &
!                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL+1) )
!          else
!            Wab(ist,1) = Wab(ist,1) + fillingL(Nc+ia,1,iL) * &
!                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL-1) )
!            Wab(ist,2) = Wab(ist,2) + fillingL(Nc+ib,1,iL) * &
!                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL-1) )
!          end if
!        end if
!      end do
!
!    end do
!
!  end subroutine getLagrangians_
!
!
!  !> calculate state-interaction terms between SA-REKS states in (2,2) case
!  subroutine getStateCoup22_(Wab, FONs, StateCoup)
!
!    !> converged Lagrangian values within active space
!    real(dp), intent(in) :: Wab(:,:)
!
!    !> Fractional occupation numbers of active orbitals
!    real(dp), intent(in) :: FONs(:,:)
!
!    !> state-interaction term between SA-REKS states
!    real(dp), intent(out) :: StateCoup(:,:)
!
!    real(dp) :: n_a, n_b
!    integer :: nstates
!
!    n_a = FONs(1,1)
!    n_b = FONs(2,1)
!    nstates = size(StateCoup,dim=1)
!
!    StateCoup(:,:) = 0.0_dp
!    StateCoup(1,2) = sqrt(n_a) * Wab(1,1) - sqrt(n_b) * Wab(1,1)
!    StateCoup(2,1) = StateCoup(1,2)
!    if (nstates == 3) then
!      StateCoup(2,3) = sqrt(n_a) * Wab(1,1) + sqrt(n_b) * Wab(1,1)
!      StateCoup(3,2) = StateCoup(2,3)
!    end if
!
!  end subroutine getStateCoup22_
!
!
!  !> Calculate factor from n_a, n_b, and delta for certain active orbital set
!  function getFactor(n_a, n_b, delta) result(factor)
!
!    !> Fractional occupation numbers of active orbitals
!    real(dp), intent(in) :: n_a, n_b
!
!    !> Smoothing factor used in FON optimization
!    real(dp), intent(in) :: delta
!
!    !> factor of n_a and n_b
!    real(dp) :: factor
!
!    factor = -0.5_dp*(n_a*n_b)**&
!        & (1.0_dp-0.5_dp*(n_a*n_b+delta)/(1.0_dp+delta))
!
!  end function getFactor


end module dftbp_reksio
