!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module dftbp_transportio
  use dftbp_accuracy
  use dftbp_constants
  use dftbp_globalenv
  use dftbp_message
  use dftbp_orbitals
#:if WITH_TRANSPORT
  use libnegf_vars
#:endif
  implicit none

  private
  public :: writeShifts, readShifts, writeContactShifts, readContactShifts

contains

  !> Write the Hamiltonian self consistent shifts to file
  subroutine writeShifts(fShifts, orb, shiftPerL)
    !> filename where shifts are stored
    character(*), intent(in) :: fShifts

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> shifts organized per (shell , atom,  spin)
    real(dp), intent(in) :: shiftPerL(:,:,:)

    ! locals
    integer :: fdHS, nSpin, nAtom, ii, jj

    nSpin = size(shiftPerL, dim=3)
    nAtom = size(shiftPerL, dim=2)

    if (size(shiftPerL, dim=1) /= orb%mShell ) then
      call error("Internal error in writeshift: size(shiftPerL,1)")
    endif

    if (size(shiftPerL, dim=2) /= size(orb%nOrbAtom) ) then
      call error("Internal error in writeshift size(shiftPerL,2)")
    endif

    open(newunit=fdHS, file=trim(fShifts), form="formatted")
    write(fdHS, *) nAtom, orb%mShell, orb%mOrb, nSpin
    do ii = 1, nAtom
      write(fdHS, *) orb%nOrbAtom(ii), (shiftPerL(:,ii,jj), jj = 1, nSpin)
    end do

    close(fdHS)

    write(stdOut,*) ">> Shifts saved for restart in shifts.dat"

  end subroutine writeShifts


  !> Read the Hamiltonian potential shifts from file
  subroutine readShifts(fShifts, orb, nAtom, nSpin, shiftPerL)
    !> filename where shifts are stored
    character(*), intent(in) :: fShifts

    !> orbital information
    type(TOrbitals), intent(in) :: orb

    !> number of atoms and spin blocks
    integer, intent(in) :: nAtom, nSpin

    !> potential shifts (shell,atom,spin) charge/mag is used
    real(dp), intent(inout) :: shiftPerL(:,:,:)

    !Locals
    integer :: fdH, nAtomSt, nSpinSt, mOrbSt, mShellSt, ii, jj
    integer, allocatable :: nOrbAtom(:)

    shiftPerL(:,:,:) = 0.0_dp

    open(newunit=fdH, file=fShifts, form="formatted")
    read(fdH, *) nAtomSt, mShellSt, mOrbSt, nSpinSt

    if (nAtomSt /= nAtom .or. mShellSt /= orb%mShell .or. mOrbSt /= orb%mOrb) then
      call error("Shift upload error: Mismatch in number of atoms or max shell per atom.")
    end if
    if (nSpin /= nSpinSt) then
      call error("Shift upload error: Mismatch in number of spin channels.")
    end if

    allocate(nOrbAtom(nAtomSt))
    do ii = 1, nAtom
      read(fdH, *) nOrbAtom(ii), (shiftPerL(:,ii,jj), jj = 1, nSpin)
    end do

    close(fdH)

    if (any(nOrbAtom(:) /= orb%nOrbAtom(:))) then
      call error("Incompatible orbitals in the upload file!")
    end if

    deallocate(nOrbAtom)

  end subroutine readShifts


  !> Writes the contact potential shifts per shell (for transport)
  subroutine writeContactShifts(filename, orb, shiftPerL, charges, Ef)

    !> filename where shifts are written
    character(*), intent(in) :: filename

    !> orbital structure
    type(TOrbitals), intent(in) :: orb

    !> array of shifts per shell and spin
    real(dp), intent(in) :: shiftPerL(:,:,:)

    !> array of charges per shell and spin
    real(dp), intent(in) :: charges(:,:,:)

    !> Fermi level
    real(dp), intent(in) :: Ef(:)

    integer :: fdHS, nAtom, nSpin

    nSpin = size(shiftPerL, dim=3)
    nAtom = size(shiftPerL, dim=2)

    if (size(shiftPerL, dim=1) /= orb%mShell) then
      call error("Internal error in writeContactShifts: size(shiftPerL,1)")
    endif

    if (size(orb%nOrbAtom) /= nAtom) then
      call error("Internal error in writeContactShifts: size(shiftPerL,2)")
    endif

    if (all(shape(charges) /= (/ orb%mOrb, nAtom, nSpin /))) then
      call error("Internal error in writeContShift: shape(charges)")
    endif

    open(newunit=fdHS, file=trim(filename), form="formatted")
    write(fdHS, *) nAtom, orb%mShell, orb%mOrb, nSpin
    write(fdHS, *) orb%nOrbAtom
    write(fdHS, *) shiftPerL
    write(fdHS, *) charges
    if (nSpin == 2) then
      write(fdHS, *) 'Fermi level (up):', Ef(1), "H", Hartree__eV * Ef(1), 'eV'
      write(fdHS, *) 'Fermi level (down):', Ef(2), "H", Hartree__eV * Ef(2), 'eV'
    else
      write(fdHS, *) 'Fermi level :', Ef(1), "H", Hartree__eV * Ef(1), 'eV'
    end if

    close(fdHS)

    write(stdOut,*) 'shiftcont_'//trim(filename)//".dat written to file"

  end subroutine writeContactShifts


#:if WITH_TRANSPORT

  !> Read contact potential shifts from file
  subroutine readContactShifts(shiftPerL, charges, tp, orb, species)

    !> shifts for atoms in contacts
    real(dp), intent(out) :: shiftPerL(:,:)

    !> charges for atoms in contacts
    real(dp), intent(out) :: charges(:,:,:)

    !> transport parameters
    type(TTransPar), intent(in) :: tp

    !> atomic orbital parameters
    type(TOrbitals), intent(in) :: orb

    !> species of atoms in the system
    integer, intent(in) :: species(:)

    real(dp), allocatable :: shiftPerLSt(:,:,:), chargesSt(:,:,:)
    integer, allocatable :: nOrbAtom(:)
    integer :: nAtomSt, mShellSt, nContAtom, mOrbSt, nSpinSt, nSpin
    integer :: iCont, iStart, iEnd, ii
    integer :: fdH
    character(lc) :: strTmp
    logical :: iexist

    nSpin = size(charges, dim=3)

    if (size(shiftPerL,dim=2) /= size(charges, dim=2)) then
      call error("Mismatch between array charges and shifts")
    endif

    shiftPerL(:,:) = 0.0_dp
    charges(:,:,:) = 0.0_dp

    do iCont = 1, tp%ncont
      inquire(file="shiftcont_"// trim(tp%contacts(iCont)%name) // ".dat", exist = iexist)
      if (.not.iexist) then
        call error("Contact shift file shiftcont_"// trim(tp%contacts(iCont)%name) &
            &  // ".dat is missing"//new_line('a')//"Run ContactHamiltonian calculations first.")
      end if

      open(newunit=fdH, file="shiftcont_" // trim(tp%contacts(iCont)%name) // ".dat",&
          & form="formatted", status="OLD", action="READ")
      read(fdH, *) nAtomSt, mShellSt, mOrbSt, nSpinSt
      iStart = tp%contacts(iCont)%idxrange(1)
      iEnd = tp%contacts(iCont)%idxrange(2)
      nContAtom = iEnd - iStart + 1

      if (nAtomSt /= nContAtom) then
        call error("Upload Contacts: Mismatch in number of atoms.")
      end if
      if (mShellSt /= orb%mShell) then
        call error("Upload Contacts: Mismatch in max shell per atom.")
      end if
      if (mOrbSt /= orb%mOrb) then
        call error("Upload Contacts: Mismatch in orbitals per atom.")
      end if
      if (nSpin /= nSpinSt) then
        write(strTmp,"(A,I0,A,I0)")'Contact spin ',nSpinSt,'. Spin channels ',nSpin
        call error(trim(strTmp))
      end if

      allocate(nOrbAtom(nAtomSt))
      read(fdH, *) nOrbAtom
      allocate(shiftPerLSt(orb%mShell, nAtomSt, nSpin))
      read(fdH, *) shiftPerLSt(:,:,:)
      allocate(chargesSt(orb%mOrb, nAtomSt, nSpin))
      read(fdH, *) chargesSt
      close(fdH)

      if (any(nOrbAtom /= orb%nOrbAtom(iStart:iEnd))) then
        call error("Incompatible orbitals in the upload file!")
      end if

      !if (nSpin == 1) then
      shiftPerL(:,iStart:iEnd) = ShiftPerLSt(:,:,1)
      !else
      !  shiftPerL(:,iStart:iEnd) = sum(ShiftPerLSt, dim=3)
      !endif

      charges(:,iStart:iEnd,:) = chargesSt(:,:,:)
      deallocate(nOrbAtom)
      deallocate(shiftPerLSt)
      deallocate(chargesSt)
    end do

  end subroutine readContactShifts

#:endif

end module dftbp_transportio
