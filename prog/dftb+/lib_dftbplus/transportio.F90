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

      open(newunit=fdHS, file="shiftcont_"//trim(filename)//".dat", position="rewind",&
          & status="replace", form="formatted")

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

      write(stdOut,*) 'shiftcont_'//trim(filename)//'.dat written to file'

    else

      open(newunit=fdHS, file="shiftcont_"//trim(filename)//".bin", position="rewind",&
          & status="replace", form="unformatted")

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

      write(stdOut,*) 'shiftcont_'//trim(filename)//'.bin written to file'

    end if

    close(fdHS)

  end subroutine writeContactShifts


#:if WITH_TRANSPORT

  !> Read contact potential shifts from file
  subroutine readContactShifts(shiftPerL, charges, tp, orb, shiftBlockUp, blockUp,&
      & tReadAscii)

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

    !> Should a text or binary file be read?
    logical, intent(in), optional :: tReadAscii

    real(dp), allocatable :: shiftPerLSt(:,:,:), chargesSt(:,:,:)
    integer, allocatable :: nOrbAtom(:)
    integer :: nAtomSt, mShellSt, nContAtom, mOrbSt, nSpinSt, nSpin
    integer :: iCont, iStart, iEnd, ii, iAt, iSp

    integer :: fdH, iErr, iBuffer(5), fileVersion, iSpin, iEntries
    character(lc) :: strTmp
    logical :: iexist, tBlock, tAsciiFile
    character(lc) :: buffer, fileName
    character(sc) :: shortBuffer

  @:ASSERT(allocated(shiftBlockUp) .eqv. allocated(blockUp))

    tAsciiFile = .true.
    if (present(tReadAscii)) then
      tAsciiFile = tReadAscii
    end if

    tBlock = .false.

    nSpin = size(charges, dim=3)

    shiftPerL(:,:) = 0.0_dp
    charges(:,:,:) = 0.0_dp

    if (allocated(shiftBlockUp)) then
      shiftblockUp(:,:,:,:) = 0.0_dp
      blockUp(:,:,:,:) = 0.0_dp
    end if

    do iCont = 1, tp%ncont

      iStart = tp%contacts(iCont)%idxrange(1)
      iEnd = tp%contacts(iCont)%idxrange(2)
      nContAtom = iEnd - iStart + 1

      if (tAsciiFile) then
        fileName = "shiftcont_"// trim(tp%contacts(iCont)%name) // ".dat"
      else
        fileName = "shiftcont_"// trim(tp%contacts(iCont)%name) // ".bin"
      end if

      inquire(file = trim(fileName), exist = iexist)
      if (.not.iexist) then
        call error("Contact shift file shiftcont_" // trim(tp%contacts(iCont)%name) &
            &  // ".dat is missing"//new_line('a') // "Run ContactHamiltonian calculations first.")
      end if

      if (tAsciiFile) then
        open(newunit=fdH, file=trim(fileName), form="formatted", status="OLD", action = "READ",&
            & iostat=iErr)
      else
        open(newunit=fdH, file=trim(fileName), form="unformatted", status="OLD", action = "READ",&
            & iostat=iErr)
      end if
      if (iErr /= 0) then
        write(strTmp, *) "Failure to open contact shift file " // trim(fileName)
        call error(strTmp)
      end if
      rewind(fdH)

      if (tAsciiFile) then

        ! need to check for original unversioned format if ascii, so dump first line into a buffer
        ! and then post process to work out what this file is
        read(fdH,'(A)',iostat=iErr) buffer
        if (iErr/=0) then
          call error("Error reading file contact file "//trim(fileName))
        endif

        ! count the integers in the first line and store them
        ! looking for either a version number or the old format start line
        do iEntries = 1, 6
          read(buffer, *, iostat=iErr) iBuffer(:iEntries)
          if (iErr==-1) then
            exit
          end if
        end do
        iEntries = iEntries - 1

      else

        ! binary format only avalable from format 1 of the contact file
        iBuffer(:) = 0
        read(fdH) iBuffer(1)
        iEntries = 1

      end if

      ! version number of the format
      select case (iEntries)
      case(1)

        ! single value present - use it as the format version
        fileVersion = iBuffer(1)

        select case(fileVersion)
        case(1)

          if (tAsciiFile) then
            read(fdH, *) nAtomSt, mShellSt, mOrbSt, nSpinSt, tBlock
          else
            read(fdH) nAtomSt, mShellSt, mOrbSt, nSpinSt, tBlock
          end if

        case default

          write(shortBuffer, "(I0)")fileVersion
          call error("Unknown contact version number in file shiftcont_" //&
              & trim(tp%contacts(iCont)%name) // ".dat : " // trim(shortBuffer))

        end select

      case(4)

        ! old format, directly start of data at top of file
        fileVersion = 0
        nAtomSt = iBuffer(1)
        mShellSt = iBuffer(2)
        mOrbSt = iBuffer(3)
        nSpinSt = iBuffer(4)

      case default
        write(stdOut,*) "Error reading file contact file " // trim(fileName) //&
            & ", first line is unexpected:"
        call error(trim(buffer))
      end select

      if (fileVersion == 0 .and. allocated(shiftBlockUp)) then
        call error("Orbital potentials not supported for this contact file")
      end if


      select case(fileVersion)

      case(0,1)

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
        allocate(shiftPerLSt(orb%mShell, nAtomSt, nSpin))
        allocate(chargesSt(orb%mOrb, nAtomSt, nSpin))

        if (tAsciiFile) then
          read(fdH, *) nOrbAtom
          read(fdH, *) shiftPerLSt(:,:,:)
          read(fdH, *) chargesSt
        else
          read(fdH) nOrbAtom
          read(fdH) shiftPerLSt(:,:,:)
          read(fdH) chargesSt
        end if

        if (tBlock .neqv. allocated(blockUp)) then
          call error("Shift and orbital potential mismatch for "//trim(tp%contacts(iCont)%name))
        end if

        if (tBlock) then
          if (tAsciiFile) then
            do iSp = 1, size(shiftBlockUp,dim=4)
              do ii = 0, iEnd-iStart
                iAt = iStart + ii
                read(fdH, *) shiftBlockUp(:orb%nOrbAtom(iAt), :orb%nOrbAtom(iAt), iAt, iSp)
              end do
            end do
            do iSp = 1, size(blockUp,dim=4)
              do ii = 0, iEnd-iStart
                iAt = iStart + ii
                read(fdH, *) blockUp(:orb%nOrbAtom(iAt),:orb%nOrbAtom(iAt),iAt,iSp)
              end do
            end do
          else
            do iSp = 1, size(shiftBlockUp,dim=4)
              do ii = 0, iEnd-iStart
                iAt = iStart + ii
                read(fdH) shiftBlockUp(:orb%nOrbAtom(iAt), :orb%nOrbAtom(iAt), iAt, iSp)
              end do
            end do
            do iSp = 1, size(blockUp,dim=4)
              do ii = 0, iEnd-iStart
                iAt = iStart + ii
                read(fdH) blockUp(:orb%nOrbAtom(iAt),:orb%nOrbAtom(iAt),iAt,iSp)
              end do
            end do
          end if
        end if

        select case(fileVersion)
        case (0)
          ! unformatted Fermi energy line(s)
          if (.not. tp%contacts(iCont)%tFermiSet) then
            call error("File format of " // trim(fileName) //&
                & " is too early to read the Fermi energy")
          end if
        case (1)
          if (.not. tp%contacts(iCont)%tFermiSet) then
            if (tAsciiFile) then
              do iSpin = 1, nSpin
                read(fdH, formatFermiRead) tp%contacts(iCont)%eFermi(iSpin)
              end do
            else
              read(fdH) tp%contacts(iCont)%eFermi
            end if
            tp%contacts(iCont)%tFermiSet = .true.
          end if
        case default
          call error("Unknown format")
        end select

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

      end select

      close(fdH)

    end do

  end subroutine readContactShifts

#:endif

end module dftbp_transportio
