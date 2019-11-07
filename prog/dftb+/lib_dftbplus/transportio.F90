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
  subroutine writeContactShifts(filename, orb, shiftPerL, charges, Ef, blockShift, blockCharges)

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

    integer :: fdHS, nAtom, nSpin, iAt, iSp

    nSpin = size(shiftPerL, dim=3)
    nAtom = size(shiftPerL, dim=2)

    open(newunit=fdHS, file=trim(filename), form="formatted")

    ! now with a version number on the top of the file:
    write(fdHS, *) contactFormatVersion

    write(fdHS, *) nAtom, orb%mShell, orb%mOrb, nSpin, allocated(blockCharges)
    write(fdHS, *) orb%nOrbAtom
    write(fdHS, *) shiftPerL
    write(fdHS, *) charges

    if (allocated(blockShift)) then
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
  subroutine readContactShifts(shiftPerL, charges, tp, orb, species, shiftBlockUp, blockUp)

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

    !> uploded block per atom
    real(dp), allocatable, intent(inout) :: shiftBlockUp(:,:,:,:)

    !> uploaded block charges for atoms
    real(dp), allocatable, intent(inout) :: blockUp(:,:,:,:)

    real(dp), allocatable :: shiftPerLSt(:,:,:), chargesSt(:,:,:)
    integer, allocatable :: nOrbAtom(:)
    integer :: nAtomSt, mShellSt, nContAtom, mOrbSt, nSpinSt, nSpin
    integer :: iCont, iStart, iEnd, ii, iAt, iSp

    integer :: fdH, iErr, iBuffer(5), fileVersion
    character(lc) :: strTmp
    logical :: iexist, tBlock
    character(lc) :: buffer
    character(sc) :: shortBuffer

  @:ASSERT(allocated(shiftBlockUp) .eqv. allocated(blockUp))

    tBlock = .false.

    nSpin = size(charges, dim=3)

    shiftPerL(:,:) = 0.0_dp
    charges(:,:,:) = 0.0_dp

    if (allocated(shiftBlockUp)) then
      shiftblockUp(:,:,:,:) = 0.0_dp
      blockUp(:,:,:,:) = 0.0_dp
    end if

    do iCont = 1, tp%ncont

      inquire(file="shiftcont_"// trim(tp%contacts(iCont)%name) // ".dat", exist = iexist)
      if (.not.iexist) then
        call error("Contact shift file shiftcont_"// trim(tp%contacts(iCont)%name) &
            &  // ".dat is missing"//new_line('a')//"Run ContactHamiltonian calculations first.")
      end if

      open(newunit=fdH, file="shiftcont_" // trim(tp%contacts(iCont)%name) // ".dat",&
          & form="formatted", status="OLD", action="READ")

      read(fdH,'(A)',iostat=iErr) buffer
      if (iErr/=0) then
        call error("Error reading file contact file shiftcont_"// trim(tp%contacts(iCont)%name) //&
            & ".dat")
      endif

      ! count the integers in the first line and store them
      ! looking for either a version number or the old format start line
      do ii = 1, 5
        read(buffer, *, iostat=iErr) iBuffer(:ii)
        if (iErr==-1) then
          exit
        end if
      end do

      ! version number of format
      select case (ii-1)
      case(1)
        ! single value present - use it as the format version
        fileVersion = iBuffer(1)

        if (fileVersion > contactFormatVersion) then
          write(shortBuffer, "(I0)")fileVersion
          call error("Unknown contact version number in file shiftcont_" //&
              & trim(tp%contacts(iCont)%name) // ".dat : " // trim(shortBuffer))
        end if

        if (fileVersion == 1) then
          read(fdH, *) nAtomSt, mShellSt, mOrbSt, nSpinSt, tBlock
        end if

      case(4)
        ! old format, directly start of data
        nAtomSt = iBuffer(1)
        mShellSt = iBuffer(2)
        mOrbSt = iBuffer(3)
        nSpinSt = iBuffer(4)
        fileVersion = 0
      case default
        write(stdOut,*) "Error reading file contact file shiftcont_" //&
            & trim(tp%contacts(iCont)%name) // ".dat, first line is unexpected:"
        call error(trim(buffer))
      end select

      if (fileVersion == 0 .and. allocated(shiftBlockUp)) then
        call error("Orbital potentials not supported for this contact")
      end if

      select case(fileVersion)
      case(0,1)
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

        if (tBlock .neqv. allocated(blockUp)) then
          call error("Shift and orbital potential mismatch for "//trim(tp%contacts(iCont)%name))
        end if
        if (tBlock) then
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
        end if

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
      end select

    end do

  end subroutine readContactShifts

#:endif

end module dftbp_transportio
