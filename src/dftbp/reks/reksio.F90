!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> REKS and SI-SA-REKS formulation in DFTB as developed by Lee et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet apart from Gamma point.
!> * Orbital potentials or spin-orbit or external E-field does not work yet.
!> * Only for closed shell system.
!> * Onsite corrections are not included in this version
module dftbp_reks_reksio
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : au__Debye
  use dftbp_common_file, only : TFileDescr, openFile, closeFile
  use dftbp_io_message, only : error
  use dftbp_reks_rekscommon, only : getTwoIndices, getSpaceSym
  use dftbp_reks_reksvar, only : TReksCalc, reksTypes

  implicit none

  private
  public :: printReksMicrostates, printSaReksEnergy
  public :: printReksSAInfo, printReksSSRInfo
  public :: printReksGradInfo
  public :: printUnrelaxedFONs, printRelaxedFONs, printRelaxedFONsL
  public :: writeReksTDP, writeReksRelaxedCharge

  contains

  !> Print energy contribution for each microstate in SCC iteration
  subroutine printReksMicrostates(this, Erep, output)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> repulsive energy
    real(dp), intent(in) :: Erep

    !> output for write processes
    integer, intent(in) :: output

    integer :: iL

    write(output,'(1x,A,5x,A,9x,A,9x,A,9x,A,10x,A,9x,A,10x,A,8x,A)') &
        & "iL", "nonSCC", "SCC", "spin", "3rd", "fock", "Rep", "Disp", "Total"
    do iL = 1, this%Lmax
      write(output,'(I3,7(f13.8))',advance="no") iL, this%enLnonSCC(iL), &
          & this%enLscc(iL), this%enLspin(iL)
      if (this%t3rd) then
        write(output,'(1(f13.8))',advance="no") this%enL3rd(iL)
      else
        write(output,'(1(f13.8))',advance="no") 0.0_dp
      end if
      if (this%isHybridXc) then
        write(output,'(1(f13.8))',advance="no") this%enLfock(iL)
      else
        write(output,'(1(f13.8))',advance="no") 0.0_dp
      end if
      write(output,'(1(f13.8))',advance="no") Erep
      if (this%isDispersion) then
        write(output,'(1(f13.8))',advance="no") this%enLdisp(iL)
      else
        write(output,'(1(f13.8))',advance="no") 0.0_dp
      end if
      write(output,'(1(f13.8))') this%enLtot(iL)
    end do

  end subroutine printReksMicrostates


  !> Print SA-REKS energy in SCC iteration
  subroutine printSaReksEnergy(this, output)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> output for write processes
    integer, intent(in) :: output

    integer :: ist

    write(output,'(1x,A)') "SA-REKS state energies"
    do ist = 1, this%nstates
      if (mod(ist,5) == 0 .or. ist == this%nstates) then
        write(output,"(I3,':',1x,1(f13.8),1x,'H')") ist, this%energy(ist)
      else
        write(output,"(I3,':',1x,1(f13.8),1x,'H')",advance="no") ist, this%energy(ist)
      end if
    end do

  end subroutine printSaReksEnergy


  !> print SA-REKS result in standard output
  subroutine printReksSAInfo(this, Eavg, output)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> Total energy for averaged state in REKS
    real(dp), intent(in) :: Eavg

    !> output for write processes
    integer, intent(in) :: output

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call printReksSAInfo22_(Eavg, this%enLtot, this%energy, this%FONs, this%Efunction,&
          & this%Plevel, output)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine printReksSAInfo


  !> print SI-SA-REKS result in standard output
  subroutine printReksSSRInfo(this, Wab, tmpEn, StateCoup, output)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> converged Lagrangian values within active space
    real(dp), intent(in) :: Wab(:,:)

    !> SA-REKS energies
    real(dp), intent(in) :: tmpEn(:)

    !> state-interaction term between SA-REKS states
    real(dp), intent(in) :: StateCoup(:,:)

    !> output for write processes
    integer, intent(in) :: output

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call printReksSSRInfo22_(output, Wab, tmpEn, StateCoup, this%energy, this%eigvecsSSR, &
          & this%Na, this%tAllStates, this%tSSR)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine printReksSSRInfo


  !> print gradient results for REKS calculation
  subroutine printReksGradInfo(this, derivs, output)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> derivatives of energy wrt to atomic positions
    real(dp), intent(in) :: derivs(:,:)

    !> output for write processes
    integer, intent(in) :: output

    integer :: ist, ia, ib, nstHalf
    character(3), parameter :: ordinals(6) = ['1st', '2nd', '3rd', '4th', '5th', '6th']

    nstHalf = this%nstates * (this%nstates - 1) / 2

    write(output,*)
    if (this%Efunction == 1) then

      write(output,"(A)") repeat("-", 50)
      write(output,"(A)") " Gradient Information"
      write(output,"(A)") repeat("-", 50)
      write(output,*) this%rstate, "state (single-state)"
      write(output,'(3(f15.8))') derivs(:,:)
      write(output,"(A)") repeat("-", 50)

    else

      if (this%tNAC) then

        write(output,"(A)") repeat("-", 50)
        write(output,"(A)") " Gradient Information"
        write(output,"(A)") repeat("-", 50)
        do ist = 1, this%nstates
          write(output,'(12X,A)') ordinals(ist) // " state (SSR)"
          write(output,'(3(f15.8))') this%SSRgrad(:,:,ist)
          if (ist == this%nstates) then
            write(output,"(A)") repeat("-", 50)
          else
            write(output,'(3(f15.8))')
          end if
        end do

        if (this%Plevel >= 2) then
          write(output,"(12X,A)") "Averaged state"
          write(output,'(3(f15.8))') this%avgGrad(:,:)
          write(output,'(3(f15.8))')
          do ist = 1, this%nstates
            write(output,'(10X,A)') ordinals(ist) // " state (SA-REKS)"
            write(output,'(3(f15.8))') this%SAgrad(:,:,ist)
            if (ist == this%nstates) then
              write(output,"(A)") repeat("-", 50)
            else
              write(output,'(3(f15.8))')
            end if
          end do
        end if

        write(output,"(A)") " Coupling Information"
        do ist = 1, nstHalf

          call getTwoIndices(this%nstates, ist, ia, ib, 1)

          write(output,"(A)") repeat("-", 50)
          write(output,'(" between ",I2," and ",I2," states")') ia, ib
          write(output,"(A)") repeat("-", 50)
          write(output,*) "g vector - difference gradient"
          write(output,'(3(f15.8))') (this%SAgrad(:,:,ia) - this%SAgrad(:,:,ib)) * 0.5_dp
          write(output,'(3(f15.8))')
          write(output,*) "h vector - derivative coupling"
          write(output,'(3(f15.8))') this%SIgrad(:,:,ist)
          write(output,'(3(f15.8))')
          write(output,*) "G vector - GDV"
          write(output,'(3(f15.8))') this%nacG(:,:,ist)
          write(output,'(3(f15.8))')
          write(output,*) "H vector - DCV - non-adiabatic coupling"
          write(output,'(3(f15.8))') this%nacH(:,:,ist)

        end do
        write(output,"(A)") repeat("-", 50)

      else

        write(output,"(A)") repeat("-", 50)
        write(output,"(A)") " Gradient Information"
        write(output,"(A)") repeat("-", 50)
        if (this%Lstate == 0) then
          if (this%tSSR) then
            write(output,*) this%rstate, "state (SSR)"
          else
            write(output,*) this%rstate, "state (SA-REKS)"
          end if
        else
          write(output,*) this%Lstate, "microstate"
        end if
        write(output,'(3(f15.8))') derivs(:,:)
        write(output,"(A)") repeat("-", 50)

      end if

    end if
    write(output,*)

  end subroutine printReksGradInfo


  !> print unrelaxed FONs for target state
  subroutine printUnrelaxedFONs(tmpRho, rstate, Lstate, Nc, Na, tSSR, output)

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

    !> output for write processes
    integer, intent(in) :: output

    integer :: ii

    write(output,*)
    if (tSSR) then
      write(output,'(A25,I1,A1)',advance="no") " unrelaxed SSR FONs for S", &
          & rstate - 1, ":"
    else
      if (Lstate == 0) then
        write(output,'(A29,I1,A1)',advance="no") " unrelaxed SA-REKS FONs for S", &
            & rstate - 1, ":"
      else
        write(output,'(A20,I1,A12)',advance="no") " unrelaxed FONs for ", &
            & Lstate, " microstate:"
      end if
    end if
    do ii = 1, Na
      if (ii == Na) then
        write(output,'(1(f10.6))') tmpRho(Nc+ii,Nc+ii)
      else
        write(output,'(1(f10.6))',advance="no") tmpRho(Nc+ii,Nc+ii)
      end if
    end do

  end subroutine printUnrelaxedFONs


  !> print Relaxed FONs for target state
  subroutine printRelaxedFONs(tmpRho, rstate, Nc, Na, tSSR, output)

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

    !> output for write processes
    integer, intent(in) :: output

    integer :: ii

    if (tSSR) then
      write(output,'(A23,I1,A1)',advance="no") " relaxed SSR FONs for S", &
          & rstate - 1, ":"
    else
      write(output,'(A27,I1,A1)',advance="no") " relaxed SA-REKS FONs for S", &
          & rstate - 1, ":"
    end if
    do ii = 1, Na
      if (ii == Na) then
        write(output,'(1(f10.6))') tmpRho(Nc+ii,Nc+ii)
      else
        write(output,'(1(f10.6))',advance="no") tmpRho(Nc+ii,Nc+ii)
      end if
    end do
    write(output,*)

  end subroutine printRelaxedFONs


  !> print Relaxed FONs for target L-th microstate
  subroutine printRelaxedFONsL(tmpRho, Lstate, Nc, Na, output)

    !> Occupation number matrix
    real(dp), intent(in) :: tmpRho(:,:)

    !> Target microstate
    integer, intent(in) :: Lstate

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> output for write processes
    integer, intent(in) :: output

    integer :: ii

    write(output,'(A18,I1,A12)',advance="no") " relaxed FONs for ", &
        & Lstate, " microstate:"
    do ii = 1, Na
      if (ii == Na) then
        write(output,'(1(f10.6))') tmpRho(Nc+ii,Nc+ii)
      else
        write(output,'(1(f10.6))',advance="no") tmpRho(Nc+ii,Nc+ii)
      end if
    end do
    write(output,*)

  end subroutine printRelaxedFONsL


  !> Write tdp.dat file with transidion dipole moment
  subroutine writeReksTDP(tdp)

    real(dp), intent(in) :: tdp(:,:)

    character(len=16), parameter :: fname = "tdp.dat"
    type(TFileDescr) :: fd
    real(dp) :: tmp
    integer :: ia, ib, ist, nstates, nstHalf

    nstHalf = size(tdp,dim=2)

    tmp = 0.5_dp * (1.0_dp + sqrt(1.0_dp + 8.0_dp*real(nstHalf,dp)))
    nstates = nint(tmp)

    call openFile(fd, fname, mode="w")
    write(fd%unit,*)
    do ist = 1, nstHalf

      call getTwoIndices(nstates, ist, ia, ib, 1)
      write(fd%unit,'(A4,I1,A8,I1,A2)') " < S", ia - 1, " | r | S", ib - 1, " >"
      write(fd%unit,'(A)',advance="no") "Transition Dipole moment (au)    : "
      write(fd%unit,'(3(f12.6))') tdp(:,ist)
      write(fd%unit,'(A)',advance="no") "Transition Dipole moment (Debye) : "
      write(fd%unit,'(3(f12.6))') tdp(:,ist) * au__Debye
      write(fd%unit,*)

    end do
    call closeFile(fd)

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
    type(TFileDescr) :: fd

    nAtom = size(qOutput,dim=2)

    call openFile(fd, fname, mode="w")
    write(fd%unit,'(A13,1X,F15.8,A4)') "total charge:", &
        & -sum(qOutput(:,:,1) - q0(:,:,1)), " (e)"
    write(fd%unit,'(1X)')
    if (Lstate == 0) then
      write(fd%unit,'(A9,I1,A18)') "relaxed S", rstate - 1, " atomic charge (e)"
    else
      write(fd%unit,'(I3,A11,A18)') Lstate, " microstate", " atomic charge (e)"
    end if
    write(fd%unit,'(3X,A18)') "atom        charge"
    do iAt = 1, nAtom
      write(fd%unit,'(2X,I5,2X,F15.8)') iAt, -sum(qOutput(:,iAt,1) - q0(:,iAt,1))
    end do
    call closeFile(fd)

  end subroutine writeReksRelaxedCharge


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> print SA-REKS(2,2) result in standard output
  subroutine printReksSAInfo22_(Eavg, enLtot, energy, FONs, Efunction, Plevel, output)

    !> Total energy for averaged state in REKS
    real(dp), intent(in) :: Eavg

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

    !> output for write processes
    integer, intent(in) :: output

    real(dp) :: n_a, n_b
    integer :: iL, Lmax, ist, nstates
    character(len=8) :: strTmp

    nstates = size(energy,dim=1)
    Lmax = size(enLtot,dim=1)

    n_a = FONs(1,1)
    n_b = FONs(2,1)

    write(output,*) " "
    write(output, "(A)") repeat("-", 50)
    if (Efunction == 1) then
      write(output,'(A25,2x,F15.8)') " Final REKS(2,2) energy:", Eavg
      write(output,*) " "
      write(output,'(A46)') " State     Energy      FON(1)    FON(2)   Spin"
      write(strTmp,'(A)') "PPS"
      write(output,'(1x,a4,1x,f13.8,1x,2(f10.6),2x,f4.2)') &
          & trim(strTmp), energy(1), n_a, n_b, 0.0_dp
    else if (Efunction == 2) then
      write(output,'(A27,2x,F15.8)') " Final SA-REKS(2,2) energy:", Eavg
      write(output,*) " "
      write(output,'(A46)') " State     Energy      FON(1)    FON(2)   Spin"
      do ist = 1, nstates
        if (ist == 1) then
          write(strTmp,'(A)') "PPS"
          write(output,'(1x,a4,1x,f13.8,1x,2(f10.6),2x,f4.2)') &
              & trim(strTmp), energy(1), n_a, n_b, 0.0_dp
        else if (ist == 2) then
          write(strTmp,'(A)') "OSS"
          write(output,'(1x,a4,1x,f13.8,1x,2(f10.6),2x,f4.2)') &
              & trim(strTmp), energy(2), 1.0_dp, 1.0_dp, 0.0_dp
        else if (ist == 3) then
          write(strTmp,'(A)') "DES"
          write(output,'(1x,a4,1x,f13.8,1x,2(f10.6),2x,f4.2)') &
              & trim(strTmp), energy(3), n_b, n_a, 0.0_dp
        end if
      end do
      write(strTmp,'(A)') "Trip"
      write(output,'(1x,a4,1x,f13.8,1x,2(f10.6),2x,f4.2)') &
          & trim(strTmp), enLtot(5), 1.0_dp, 1.0_dp, 1.0_dp
    end if
    write(output, "(A)") repeat("-", 50)

    if (Plevel >= 2) then
      write(output,*) " "
      write(output, "(A)") repeat("-", 25)
      write(output,'(1x,A20,2x,F15.8)') " Microstate Energies"
      do iL = 1, Lmax
        write(output,"(1x,'L =',1x,I2,':',1x,F13.8)") iL, enLtot(iL)
      end do
      write(output, "(A)") repeat("-", 25)
    end if

  end subroutine printReksSAInfo22_


  !> print SI-SA-REKS(2,2) result in standard output
  subroutine printReksSSRInfo22_(output, Wab, tmpEn, StateCoup, energy, eigvecsSSR, &
      & Na, tAllStates, tSSR)

    !> output for write processes
    integer, intent(in) :: output

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

    write(output,*)
    do ist = 1, nActPair

      call getTwoIndices(Na, ist, ia, ib, 1)

      call getSpaceSym(ia, stA)
      call getSpaceSym(ib, stB)

      write(output,"(1x,'Lagrangian W',A1,A1,': ',2(f12.8))") &
          & trim(stA), trim(stB), Wab(ist,1), Wab(ist,2)

    end do

    write(output,*)
    write(output, "(A)") repeat("-", 50)
    if (.not. tAllStates) then
      write(output,'(A)') " SSR: 2SI-2SA-REKS(2,2) Hamiltonian matrix"
      write(output,'(15x,A3,11x,A3)') "PPS", "OSS"
    else
      write(output,'(A)') " SSR: 3SI-2SA-REKS(2,2) Hamiltonian matrix"
      write(output,'(15x,A3,11x,A3,11x,A3)') "PPS", "OSS", "DES"
    end if

    do ist = 1, nstates
      if (ist == 1) then
        write(strTmp,'(A)') "PPS"
      else if (ist == 2) then
        write(strTmp,'(A)') "OSS"
      else if (ist == 3) then
        write(strTmp,'(A)') "DES"
      end if
      write(output,'(1x,a5,1x)',advance="no") trim(strTmp)
      do jst = 1, nstates
        if (ist == jst) then
          if (jst == nstates) then
            write(output,'(1x,f13.8)') tmpEn(ist)
          else
            write(output,'(1x,f13.8)',advance="no") tmpEn(ist)
          end if
        else
          if (jst == nstates) then
            write(output,'(1x,f13.8)') StateCoup(ist,jst)
          else
            write(output,'(1x,f13.8)',advance="no") StateCoup(ist,jst)
          end if
        end if
      end do
    end do
    write(output, "(A)") repeat("-", 50)

    if (tSSR) then
      write(output,*)
      write(output, "(A)") repeat("-", 64)
      if (.not. tAllStates) then
        write(output,'(A)') " SSR: 2SI-2SA-REKS(2,2) states"
        write(output,'(19x,A4,7x,A7,4x,A7)') "E_n", "C_{PPS}", "C_{OSS}"
      else
        write(output,'(A)') " SSR: 3SI-2SA-REKS(2,2) states"
        write(output,'(19x,A4,7x,A7,4x,A7,4x,A7)') "E_n", "C_{PPS}", "C_{OSS}", "C_{DES}"
      end if
      do ist = 1, nstates
        if (.not. tAllStates) then
          write(output,'(1x,A,I2,1x,f13.8,1x,f10.6,1x,f10.6)') &
              & "SSR state ", ist, energy(ist), eigvecsSSR(:,ist)
        else
          write(output,'(1x,A,I2,1x,f13.8,1x,f10.6,1x,f10.6,1x,f10.6)') &
              & "SSR state ", ist, energy(ist), eigvecsSSR(:,ist)
        end if
      end do
      write(output, "(A)") repeat("-", 64)
    end if

  end subroutine printReksSSRInfo22_


end module dftbp_reks_reksio
