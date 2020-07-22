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
  subroutine printReksMicrostates(this, Erep)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> repulsive energy
    real(dp), intent(in) :: Erep

    integer :: iL

    write(stdOut,'(1x,A,5x,A,9x,A,9x,A,9x,A,8x,A,9x,A,8x,A)') &
        & "iL", "nonSCC", "SCC", "spin", "3rd", "fock", "Rep", "Total"
    do iL = 1, this%Lmax
      write(stdOut,'(I3,7(f13.8))',advance="no") iL, this%enLnonSCC(iL), &
          & this%enLscc(iL), this%enLspin(iL)
      if (this%t3rd) then
        write(stdOut,'(1(f13.8))',advance="no") this%enL3rd(iL)
      else
        write(stdOut,'(1(f13.8))',advance="no") 0.0_dp
      end if
      if (this%isRangeSep) then
        write(stdOut,'(1(f13.8))',advance="no") this%enLfock(iL)
      else
        write(stdOut,'(1(f13.8))',advance="no") 0.0_dp
      end if
      write(stdOut,'(2(f13.8))') Erep, this%enLtot(iL)
    end do

  end subroutine printReksMicrostates


  !> Print SA-REKS energy in SCC iteration
  subroutine printSaReksEnergy(this)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    integer :: ist

    write(stdOut,'(1x,A)') "SA-REKS state energies"
    do ist = 1, this%nstates
      if (mod(ist,5) == 0 .or. ist == this%nstates) then
        write(stdOut,"(I3,':',1x,1(f13.8),1x,'H')") ist, this%energy(ist)
      else
        write(stdOut,"(I3,':',1x,1(f13.8),1x,'H')",advance="no") ist, this%energy(ist)
      end if
    end do

  end subroutine printSaReksEnergy


  !> print SA-REKS result in standard output
  subroutine printReksSAInfo(this, Etotal)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> state-averaged energy
    real(dp), intent(in) :: Etotal

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call printReksSAInfo22_(Etotal, this%enLtot, this%energy, this%FONs, this%Efunction,&
          & this%Plevel)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine printReksSAInfo


  !> print SI-SA-REKS result in standard output
  subroutine printReksSSRInfo(this, Wab, tmpEn, StateCoup)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> converged Lagrangian values within active space
    real(dp), intent(in) :: Wab(:,:)

    !> SA-REKS energies
    real(dp), intent(in) :: tmpEn(:)

    !> state-interaction term between SA-REKS states
    real(dp), intent(in) :: StateCoup(:,:)

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call printReksSSRInfo22_(Wab, tmpEn, StateCoup, this%energy, this%eigvecsSSR, &
          & this%Na, this%tAllStates, this%tSSR)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine printReksSSRInfo


  !> print gradient results for REKS calculation
  subroutine printReksGradInfo(this, derivs)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> derivatives of energy wrt to atomic positions
    real(dp), intent(in) :: derivs(:,:)

    integer :: ist, ia, ib, nstHalf

    nstHalf = this%nstates * (this%nstates - 1) / 2

    write(stdOut,*)
    if (this%Efunction == 1) then

      write(stdOut,"(A)") repeat("-", 50)
      write(stdOut,"(A)") " Gradient Information"
      write(stdOut,"(A)") repeat("-", 50)
      write(stdOut,*) this%rstate, "state (single-state)"
      write(stdOut,'(3(f15.8))') derivs(:,:)
      write(stdOut,"(A)") repeat("-", 50)

    else

      if (this%tNAC) then

        write(stdOut,"(A)") repeat("-", 50)
        write(stdOut,"(A)") " Gradient Information"
        write(stdOut,"(A)") repeat("-", 50)
        do ist = 1, this%nstates
          write(stdOut,*) ist, "st state (SSR)"
          write(stdOut,'(3(f15.8))') this%SSRgrad(:,:,ist)
          if (ist == this%nstates) then
            write(stdOut,"(A)") repeat("-", 50)
          else
            write(stdOut,'(3(f15.8))')
          end if
        end do

!        write(stdOut,*) "AVG state"
!        write(stdOut,'(3(f15.8))') this%avgGrad(:,:)
!        write(stdOut,'(3(f15.8))')
!        do ist = 1, this%nstates
!          write(stdOut,*) ist, "st state (SA-REKS)"
!          write(stdOut,'(3(f15.8))') this%SAgrad(:,:,ist)
!          if (ist == this%nstates) then
!            write(stdOut,"(A)") repeat("-", 50)
!          else
!            write(stdOut,'(3(f15.8))')
!          end if
!        end do

        write(stdOut,"(A)") " Coupling Information"
        do ist = 1, nstHalf

          call getTwoIndices(this%nstates, ist, ia, ib, 1)

          write(stdOut,"(A)") repeat("-", 50)
          write(stdOut,'(" between ",I2," and ",I2," states")') ia, ib
          write(stdOut,"(A)") repeat("-", 50)
          write(stdOut,*) "g vector - difference gradient"
          write(stdOut,'(3(f15.8))') (this%SAgrad(:,:,ia) - this%SAgrad(:,:,ib)) * 0.5_dp
          write(stdOut,'(3(f15.8))')
          write(stdOut,*) "h vector - derivative coupling"
          write(stdOut,'(3(f15.8))') this%SIgrad(:,:,ist)
          write(stdOut,'(3(f15.8))')
          write(stdOut,*) "G vector - GDV"
          write(stdOut,'(3(f15.8))') this%nacG(:,:,ist)
          write(stdOut,'(3(f15.8))')
          write(stdOut,*) "H vector - DCV - non-adiabatic coupling"
          write(stdOut,'(3(f15.8))') this%nacH(:,:,ist)

        end do
        write(stdOut,"(A)") repeat("-", 50)

      else

        write(stdOut,"(A)") repeat("-", 50)
        write(stdOut,"(A)") " Gradient Information"
        write(stdOut,"(A)") repeat("-", 50)
        if (this%Lstate == 0) then
          if (this%tSSR) then
            write(stdOut,*) this%rstate, "state (SSR)"
          else
            write(stdOut,*) this%rstate, "state (SA-REKS)"
          end if
        else
          write(stdOut,*) this%Lstate, "microstate"
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


end module dftbp_reksio
