!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
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
  public :: writeShifts, readShifts, writeContactShifts
#:if WITH_TRANSPORT
  public :: readContactShifts
#:endif

  integer, parameter :: contactFormatVersion = 1

  character(len=*), parameter :: formatFermiWrite = "(1X,A,T20,F22.18,1X,A,F22.18,1X,A)"
  character(len=*), parameter :: formatFermiRead = "(T20, F22.18)"

contains

  !> Write the Hamiltonian self consistent shifts to file
  subroutine writeShifts(fShifts, orb, shiftPerL)

    !> filename where shifts are stored
    character(*), intent(in) :: fShifts

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> shifts organized per (shell , atom,  spin)
    real(dp), intent(in) :: shiftPerL(:,:,:)

    integer :: fdHS, nSpin, nAtom, ii, jj

    nSpin = size(shiftPerL, dim=3)
    nAtom = size(shiftPerL, dim=2)

    if (size(shiftPerL, dim=1) /= orb%mShell ) then
      call error("Internal error in writeshift: size(shiftPerL,1)")
    endif

    if (size(shiftPerL, dim=2) /= size(orb%nOrbAtom) ) then
      call error("Internal error in writeshift size(shiftPerL,2)")
    endif

    open(newunit=fdHS, file=fShifts, form="formatted")
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

    if (any(nOrbAtom /= orb%nOrbAtom)) then
      call error("Incompatible orbitals in the upload file!")
    end if

  end subroutine readShifts


  !> Writes the contact potential shifts per shell (for transport)
  subroutine writeContactShifts(filename, orb, shiftPerL, charges, Ef, blockShift, blockCharges,&
      & tWriteAscii)

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

    !> block shifts (for DFTB+U etc.)
    real(dp), allocatable, intent(in) :: blockShift(:,:,:,:)

    !> block charge populations
    real(dp), allocatable, intent(in) :: blockCharges(:,:,:,:)

    !> Should a text or binary file be saved
    logical, intent(in), optional :: tWriteAscii

    integer :: fdHS, nAtom, nSpin, iAt, iSp
    logical :: tAsciiFile

    nSpin = size(shiftPerL, dim=3)
    nAtom = size(shiftPerL, dim=2)

    tAsciiFile = .true.
    if (present(tWriteAscii)) then
      tAsciiFile = tWriteAscii
    end if

    if (tAsciiFile) then

      open(newunit=fdHS, file="shiftcont_" // trim(filename) // ".dat", form="formatted")

      ! now with a version number on the top of the file:
      write(fdHS, *) contactFormatVersion

      write(fdHS, *) nAtom, orb%mShell, orb%mOrb, nSpin, allocated(blockCharges)
      write(fdHS, *) orb%nOrbAtom
      write(fdHS, *) shiftPerL
      write(fdHS, *) charges

      if (allocated(blockCharges)) then
        do iSp = 1, nSpin
          do iAt = 1, nAtom
            write(fdHS, *) blockShift(:orb%nOrbAtom(iAt), :orb%nOrbAtom(iAt), iAt, iSp)
          end do
        end do
        do iSp = 1, nSpin
          do iAt = 1, nAtom
            write(fdHS, *) blockCharges(:orb%nOrbAtom(iAt), :orb%nOrbAtom(iAt), iAt, iSp)
          end do
        end do
      end if

      if (nSpin == 2) then
        write(fdHS, formatFermiWrite) 'Fermi level (up):', Ef(1), "H", Hartree__eV * Ef(1), 'eV'
        write(fdHS, formatFermiWrite) 'Fermi level (down):', Ef(2), "H", Hartree__eV * Ef(2), 'eV'
      else
        write(fdHS, formatFermiWrite) 'Fermi level :', Ef(1), "H", Hartree__eV * Ef(1), 'eV'
      end if

      write(stdOut,*) 'shiftcont_' // trim(filename) // '.dat written to file'

    else

      open(newunit=fdHS, file="shiftcont_" // trim(filename) // ".bin", form="unformatted")

      ! now with a version number on the top of the file:
      write(fdHS) contactFormatVersion

      write(fdHS) nAtom, orb%mShell, orb%mOrb, nSpin, allocated(blockCharges)
      write(fdHS) orb%nOrbAtom
      write(fdHS) shiftPerL
      write(fdHS) charges

      if (allocated(blockCharges)) then
        do iSp = 1, nSpin
          do iAt = 1, nAtom
            write(fdHS) blockShift(:orb%nOrbAtom(iAt), :orb%nOrbAtom(iAt), iAt, iSp)
          end do
        end do
        do iSp = 1, nSpin
          do iAt = 1, nAtom
            write(fdHS) blockCharges(:orb%nOrbAtom(iAt), :orb%nOrbAtom(iAt), iAt, iSp)
          end do
        end do
      end if

      write(fdHS) Ef(:)

      write(stdOut,*) 'shiftcont_' // trim(filename) // '.bin written to file'

    end if

    close(fdHS)


  end subroutine writeContactShifts


#:if WITH_TRANSPORT

  !> Read contact potential shifts from file
  subroutine readContactShifts(shiftPerL, charges, tp, orb, shiftBlockUp, blockUp)

    !> shifts for atoms in contacts
    real(dp), intent(out) :: shiftPerL(:,:)

    !> charges for atoms in contacts
    real(dp), intent(out) :: charges(:,:,:)

    !> transport parameters
    type(TTransPar), intent(inout) :: tp

    !> atomic orbital parameters
    type(TOrbitals), intent(in) :: orb

    !> uploded block per atom
    real(dp), allocatable, intent(inout) :: shiftBlockUp(:,:,:,:)

    !> uploaded block charges for atoms
    real(dp), allocatable, intent(inout) :: blockUp(:,:,:,:)

    integer :: fdH
    integer :: iCont
    logical :: tExist, tBlock, tAsciiFile
    character(lc) :: buffer
    integer :: iBuffer(5), fileVersion
    integer :: ii, iErr

    tBlock = .false.
    tAsciiFile = .not.tp%tReadBinShift

    shiftPerL(:,:) = 0.0_dp
    charges(:,:,:) = 0.0_dp

  @:ASSERT(allocated(blockUp) .eqv. allocated(shiftBlockUp))
    if (allocated(shiftBlockUp)) then
      shiftblockUp(:,:,:,:) = 0.0_dp
      blockUp(:,:,:,:) = 0.0_dp
      tBlock = .true.
    end if

    do iCont = 1, tp%ncont

      if (tAsciiFile) then

        inquire(file="shiftcont_" // trim(tp%contacts(iCont)%name) // ".dat", exist=tExist)
        if (.not. tExist) then
          call error("Contact shift file shiftcont_" // trim(tp%contacts(iCont)%name) &
              &  // ".dat is missing"// new_line('a') // "Run ContactHamiltonian calculations&
              & first.")
        end if

        open(newunit=fdH, file="shiftcont_" // trim(tp%contacts(iCont)%name) // ".dat",&
            & form="formatted", status="old", action="read")

        read(fdH, '(A)', iostat=iErr) buffer
        if (iErr /= 0) then
          call error("Error reading file contact file shiftcont_" // trim(tp%contacts(iCont)%name)&
              & // ".dat")
        endif

        ! count the integers in the first line and store them
        ! looking for either a version number or the old format start line
        do ii = 1, 5
          read(buffer, *, iostat=iErr) iBuffer(:ii)
          if (iErr /= 0) then
            exit
          end if
        end do

        select case (ii - 1)

        case (1)
          ! New file format with version number
          fileVersion = iBuffer(1)
          select case (fileVersion)
          case (1)
            ! Same format as old version (apart of the version number)
            call readContactShiftData1(fdH, orb, tp, iCont, shiftPerL, charges, shiftBlockUp,&
                & blockUp)
          case default
            write(buffer, "(I0)") fileVersion
            call error("Unknown contact version number in file shiftcont_" //&
                & trim(tp%contacts(iCont)%name) // ".dat : " // trim(buffer))
          end select

        case (4)

          ! Old version without format version number
          ! Re-read first line, since it contains relevant data instead of version number

          if (tBlock) then
            call error("File format of the " // trim(tp%contacts(iCont)%name) //&
                & " contact is too early for orbital potential support")
          end if

          if (.not. tp%contacts(iCont)%tFermiSet) then
            call error("File format of the " // trim(tp%contacts(iCont)%name) //&
                & " contact is too early to read the Fermi energy")
          end if

          rewind(fdH)
          call readContactShiftData1(fdH, orb, tp, iCont, shiftPerL, charges, shiftBlockUp, blockUp)

        case default

          write(stdOut,*) "Error reading file contact file shiftcont_" //&
              & trim(tp%contacts(iCont)%name) // ".dat"
          call error(trim(buffer))

        end select

      else

        inquire(file="shiftcont_" // trim(tp%contacts(iCont)%name) // ".bin", exist=tExist)
        if (.not. tExist) then
          call error("Contact shift file shiftcont_" // trim(tp%contacts(iCont)%name) &
              &  // ".bin is missing"// new_line('a') // "Run ContactHamiltonian calculations&
              & first.")
        end if

        open(newunit=fdH, file="shiftcont_" // trim(tp%contacts(iCont)%name) // ".bin",&
            & form="unformatted", status="old", action="read")

        read(fdH) fileVersion

        select case (fileVersion)
        case (1)
          ! From version 1, binary format supported, no need to check for earlier versions
          call readContactShiftData1(fdH, orb, tp, iCont, shiftPerL, charges, shiftBlockUp,&
              & blockUp)
        case default
          write(buffer, "(I0)") fileVersion
          call error("Unknown contact version number in file shiftcont_" //&
              & trim(tp%contacts(iCont)%name) // ".bin : " // trim(buffer))
        end select

      end if

      close(fdH)

    end do

  end subroutine readContactShifts


  !> Reads in contact shift data (format 1, ascii)
  subroutine readContactShiftData1(fdH, orb, tp, iCont, shiftPerL, charges, shiftBlockUp, blockUp)

    !> File handler
    integer, intent(in) :: fdH

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Transport parameters
    type(TTransPar), intent(inout) :: tp

    !> Contact to read
    integer, intent(in) :: iCont

    !> Shifts for atoms in contacts
    real(dp), intent(inout) :: shiftPerL(:,:)

    !> Charges for atoms in contacts
    real(dp), intent(inout) :: charges(:,:,:)

    !> uploded block per atom
    real(dp), allocatable, intent(inout) :: shiftBlockUp(:,:,:,:)

    !> uploaded block charges for atoms
    real(dp), allocatable, intent(inout) :: blockUp(:,:,:,:)

    real(dp), allocatable :: shiftPerLSt(:,:,:), chargesSt(:,:,:)
    integer, allocatable :: nOrbAtom(:)
    integer :: nAtomSt, mShellSt, nContAtom, mOrbSt, nSpinSt
    integer :: iStart, iEnd, iSpin, nSpin, iAt, ii
    character(lc) :: strTmp
    logical :: tAsciiFile

    nSpin = size(charges, dim=3)

    tAsciiFile = .not.tp%tReadBinShift

    if (tAsciiFile) then
      read(fdH, *) nAtomSt, mShellSt, mOrbSt, nSpinSt
    else
      read(fdH) nAtomSt, mShellSt, mOrbSt, nSpinSt
    end if
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
    if (nSpinSt /= nSpin) then
      write(strTmp,"(A,I0,A,I0)")'Contact spin ',nSpinSt,'. Expected spin channels ',nSpin
      call error(trim(strTmp))
    end if

    allocate(nOrbAtom(nAtomSt))
    allocate(shiftPerLSt(orb%mShell, nAtomSt, nSpin))
    allocate(chargesSt(orb%mOrb, nAtomSt, nSpin))
    if (tAsciiFile) then
      read(fdH, *) nOrbAtom
      read(fdH, *) shiftPerLSt
      read(fdH, *) chargesSt
    else
      read(fdH) nOrbAtom
      read(fdH) shiftPerLSt
      read(fdH) chargesSt
    end if

    if (any(nOrbAtom /= orb%nOrbAtom(iStart:iEnd))) then
      call error("Incompatible orbitals in the upload file!")
    end if

    if (allocated(shiftBlockUp)) then
      if (tAsciiFile) then
        do iSpin = 1, nSpin
          do ii = 0, iEnd-iStart
            iAt = iStart + ii
            read(fdH, *) shiftBlockUp(:orb%nOrbAtom(iAt), :orb%nOrbAtom(iAt), iAt, iSpin)
          end do
        end do
        do iSpin = 1, nSpin
          do ii = 0, iEnd-iStart
            iAt = iStart + ii
            read(fdH, *) blockUp(:orb%nOrbAtom(iAt),:orb%nOrbAtom(iAt),iAt,iSpin)
          end do
        end do
      else
        do iSpin = 1, nSpin
          do ii = 0, iEnd-iStart
            iAt = iStart + ii
            read(fdH) shiftBlockUp(:orb%nOrbAtom(iAt), :orb%nOrbAtom(iAt), iAt, iSpin)
          end do
        end do
        do iSpin = 1, nSpin
          do ii = 0, iEnd-iStart
            iAt = iStart + ii
            read(fdH) blockUp(:orb%nOrbAtom(iAt),:orb%nOrbAtom(iAt),iAt,iSpin)
          end do
        end do
      end if
    end if

    if (.not. tp%contacts(iCont)%tFermiSet) then
      if (tAsciiFile) then
        do iSpin = 1, nSpin
          read(fdH, formatFermiRead) tp%contacts(iCont)%eFermi(iSpin)
        end do
      else
        read(fdH) tp%contacts(iCont)%eFermi(:nSpin)
      end if
      tp%contacts(iCont)%tFermiSet = .true.
    end if

    !if (nSpin == 1) then
    shiftPerL(:,iStart:iEnd) = ShiftPerLSt(:,:,1)
    !else
    !  shiftPerL(:,iStart:iEnd) = sum(ShiftPerLSt, dim=3)
    !endif

    charges(:,iStart:iEnd,:) = chargesSt(:,:,:)

  end subroutine readContactShiftData1

#:endif

end module dftbp_transportio
