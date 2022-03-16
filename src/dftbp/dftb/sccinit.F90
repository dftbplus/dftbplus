!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for initializing SCC part of the calculation
module dftbp_dftb_sccinit
  use dftbp_common_accuracy, only : dp, elecTolMax
  use dftbp_common_globalenv, only : stdOut
  use dftbp_io_message, only : error
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_multipole, only : TMultipole
  implicit none

  private
  public :: initQFromAtomChrg, initQFromShellChrg, initQFromFile, writeQToFile
  public :: initQFromUsrChrg

  !> version number for restart format, please increment if you change the interface.
  integer, parameter :: restartFormat = 6

contains


  !> Initialise the charge vector from the reference atomic charges
  subroutine initQFromAtomChrg(fOrb, qAtom, fRefShell, species, speciesNames, orb)

    !> The number of electrons per lm,atom,spin
    real(dp), intent(out) :: fOrb(:,:,:)

    !> The charges per atom.
    real(dp), intent(in) :: qAtom(:)

    !> Occupation of each shell in the neutral atom.
    real(dp), intent(in) :: fRefShell(:,:)

    !> List of chemical species for each atom
    integer, intent(in) :: species(:)

    !> Names of the species (for error messages)
    character(len=*), intent(in) :: speciesNames(:)

    !> Information about the orbitals.
    type(TOrbitals), intent(in) :: orb

    integer :: iAt, iSp, iSh1, nAtom, iShL, iShR, nSh
    real(dp) :: fShell, fAtomRes

    nAtom = size(orb%nOrbAtom)

    @:ASSERT(size(fOrb, dim=1) == orb%mOrb)
    @:ASSERT(size(fOrb, dim=2) == nAtom)
    @:ASSERT(size(fOrb, dim=3) >= 1)
    @:ASSERT(size(fRefShell, dim=1) >= orb%mShell)
    @:ASSERT(size(fRefShell, dim=2) == size(orb%nShell))
    @:ASSERT(size(species) == nAtom)
    @:ASSERT(size(speciesNames) == size(orb%nShell))

    fOrb = 0.0_dp
    ! fill degenerately over m for each shell l
    do iAt = 1, nAtom
      iSp = species(iAt)
      ! nr. of electrons = number of electrons in all shells - gross charge
      fAtomRes = sum(fRefShell(1:orb%nShell(iSp),iSp)) - qAtom(iAt)
      lpShell: do iSh1 = 1, orb%nShell(iSp)
        fShell = min(fAtomRes, real(2 * (2 * orb%angShell(iSh1, iSp) + 1), dp))
        iShL = orb%posShell(iSh1, iSp)
        iShR = orb%posShell(iSh1+1, iSp) - 1
        nSh = iShR-iShL+1
        fOrb(iShL:iShR, iAt, 1) = fShell / real(nSh, dp)
        fAtomRes = fAtomRes - fShell
        if (fAtomRes <= 0.0_dp) then
          exit lpShell
        end if
      end do lpShell
      if (fAtomRes > 1e-4_dp) then
        call error("Not enough orbitals on species '" // trim(speciesNames(iSp))&
            &// "' to hold specified atomic charges")
      end if
    end do

  end subroutine initQFromAtomChrg


  !> Initialise the charge vector from the reference atomic charges results in a set of charges
  !> appropriate for the neutral spin unpolarised atom reference system that DFTB assumes for
  !> SCC/spin extensions
  subroutine initQFromShellChrg(qq, qShell, species, orb)

    !> The charges per lm,atom,spin
    real(dp), intent(out) :: qq(:,:,:)

    !> The reference charges per shell
    real(dp), intent(in) :: qShell(:,:)

    !> List of chemical species for each atom
    integer, intent(in) :: species(:)

    !> Information about the orbitals
    type(TOrbitals), intent(in) :: orb

    integer :: iAt1, iSp1, iSh1, nAtom, iSh1l, iSh1r, nSh1

    nAtom = size(orb%nOrbAtom)

    @:ASSERT(size(qq, dim=1) == orb%mOrb)
    @:ASSERT(size(qq, dim=2) == nAtom)
    @:ASSERT(size(qq, dim=3) >= 1)
    @:ASSERT(size(qShell, dim=1) == orb%mShell)
    @:ASSERT(size(qShell, dim=2) == size(orb%angShell, dim=2))
    @:ASSERT(size(species) == nAtom)

    qq(:,:,:) = 0.0_dp
    ! fill degenerately over m for each shell l
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      do iSh1 = 1, orb%nShell(iSp1)
        iSh1l = orb%posShell(iSh1, iSp1)
        iSh1r = orb%posShell(iSh1+1, iSp1) - 1
        nSh1 = iSh1r-iSh1l+1
        qq(iSh1l:iSh1r, iAt1, 1) = qShell(iSh1, iSp1) / real(nSh1, dp)
      end do
    end do

  end subroutine initQFromShellChrg

  !> Initialise charge vector from user-defined reference atomic-charges.
  subroutine initQFromUsrChrg(qq, qAtShell, species, orb)

    !> The charges per lm,atom,spin
    real(dp), intent(out) :: qq(:,:,:)

    !> The reference charges per shell per Atom
    real(dp), intent(in) :: qAtShell(:,:)

    !> List of chemical species for each atom
    integer, intent(in) :: species(:)

    !> Information about the orbitals
    type(TOrbitals), intent(in) :: orb

    integer :: iAt1, iSp1, iSh1, nAtom, iSh1l, iSh1r, nSh1

    nAtom = size(orb%nOrbAtom)

    @:ASSERT(size(qq, dim=1) == orb%mOrb)
    @:ASSERT(size(qq, dim=2) == nAtom)
    @:ASSERT(size(qq, dim=3) >= 1)
    @:ASSERT(size(qAtShell, dim=1) == orb%mShell)
    @:ASSERT(size(qAtShell, dim=2) == nAtom)
    @:ASSERT(size(species) == nAtom)

    qq(:,:,:) = 0.0_dp

    ! fill degenerately over m for each shell l
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      do iSh1 = 1, orb%nShell(iSp1)
        iSh1l = orb%posShell(iSh1, iSp1)
        iSh1r = orb%posShell(iSh1+1, iSp1) - 1
        nSh1 = iSh1r-iSh1l+1
        qq(iSh1l:iSh1r, iAt1, 1) = qAtShell(iSh1, iAt1) / real(nSh1, dp)
      end do
    end do

  end subroutine initQFromUsrChrg


  !> Initialise the charge vector from a named external file. Check the total
  !> charge matches that expected for the calculation.
  !> Should test of the input, if the number of orbital charges per atom match the number from the
  !> angular momentum.
  subroutine initQFromFile(qq, fileName, tReadAscii, orb, qBlock, qiBlock, deltaRho,&
      & nAtInCentralRegion, magnetisation, nEl, multipoles)

    !> The charges per lm,atom,spin
    real(dp), intent(out) :: qq(:,:,:)

    !> The external file of charges for the orbitals, currently stored with each line containing the
    !> per-orbital charges in order of increasing m and l. Alternating lines give the spin case (if
    !> present)
    character(*), intent(in) :: fileName

    !> Should charges be read as ascii (cross platform, but potentially lower reproducibility) or
    !> binary files
    logical, intent(in) :: tReadAscii

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> block Mulliken population for LDA+U etc
    real(dp), intent(inout), allocatable :: qBlock(:,:,:,:)

    !> block Mulliken imagninary population for LDA+U and L.S
    real(dp), intent(inout), allocatable :: qiBlock(:,:,:,:)

    !> Full density matrix with on-diagonal adjustment
    real(dp), intent(inout), allocatable :: deltaRho(:)

    !> Number of atoms in central region (atoms outside this will have charges suplied from
    !> elsewhere)
    integer, intent(in) :: nAtInCentralRegion

    !> magnetisation checksum for regular spin polarization total magnetic moment
    real(dp), intent(in), optional :: magnetisation

    !> Nr. of electrons for each spin channel
    real(dp), intent(in), optional :: nEl

    !> Atomic multipoles, if relevant
    type(TMultipole), intent(inout), optional :: multipoles

    !> nr. of orbitals / atoms / spin channels
    integer :: nOrb, nAtom, nSpin

    !> error returned by the io commands
    integer :: iErr

    !> file unit number
    integer :: file

    !> total charge is present at the top of the file
    real(dp) :: CheckSum(size(qq, dim=3))

    integer :: iOrb, iAtom, iSpin, ii, nAtomInFile, nDipole, nQuadrupole
    integer :: fileFormat
    real(dp) :: sumQ

    !> present in the file itself
    logical :: tBlockPresent, tiBlockPresent, tRhoPresent

    !> requested to be re-loaded
    logical :: tBlock, tiBlock, tRho, isMultipolar

    character(len=120) :: error_string

    nAtom = size(qq, dim=2)
    nSpin = size(qq, dim=3)

    tBlock = allocated(qBlock)
    tiBlock = allocated(qiBlock)
    tRho = allocated(deltaRho)

    @:ASSERT(size(qq, dim=1) == orb%mOrb)
    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)

  #:block DEBUG_CODE

    if (present(magnetisation)) then
      @:ASSERT(nSpin==2)
    end if

    if (tBlock) then
      @:ASSERT(all(shape(qBlock) == (/orb%mOrb,orb%mOrb,nAtom,nSpin/)))
    end if

    if (tiBlock) then
      @:ASSERT(tBlock)
      @:ASSERT(all(shape(qiBlock) == shape(qBlock)))
    end if

    if (tRho) then
      @:ASSERT(size(deltaRho) == orb%nOrb*orb%nOrb*nSpin)
    end if

  #:endblock DEBUG_CODE

    if (tReadAscii) then
      open(newunit=file, file=trim(fileName)//'.dat', status='old', action='READ', iostat=iErr)
    else
      open(newunit=file, file=trim(fileName)//'.bin', status='old', action='READ',&
          & form='unformatted',iostat=iErr)
    end if
    if (iErr /= 0) then
      write(error_string, *) "Failure to open external file of charge data"
      call error(error_string)
    end if
    rewind(file)

    if (tReadAscii) then
      read(file, *, iostat=iErr)fileFormat
    else
      read(file, iostat=iErr)fileFormat
    end if
    if (iErr /= 0) then
      call error("Error during reading external file of charge data")
    end if
    isMultipolar = .false.
    select case(fileFormat)
    case(4)
      if (tReadAscii) then
        read(file, *, iostat=iErr)tBlockPresent, tiBlockPresent, tRhoPresent, iSpin, CheckSum
      else
        read(file, iostat=iErr)tBlockPresent, tiBlockPresent, tRhoPresent, iSpin, CheckSum
      end if
      nAtomInFile = nAtom
    case(5)
      if (tReadAscii) then
        read(file, *, iostat=iErr)tBlockPresent, tiBlockPresent, tRhoPresent, nAtomInFile, iSpin,&
            & CheckSum
      else
        read(file, iostat=iErr)tBlockPresent, tiBlockPresent, tRhoPresent, nAtomInFile, iSpin,&
            & CheckSum
      end if
    case(6)
      if (tReadAscii) then
        read(file, *, iostat=iErr)tBlockPresent, tiBlockPresent, tRhoPresent, isMultipolar,&
            & nAtomInFile, iSpin, CheckSum
      else
        read(file, iostat=iErr)tBlockPresent, tiBlockPresent, tRhoPresent, isMultipolar,&
            & nAtomInFile, iSpin, CheckSum
      end if
    case default
      call error("Incompatible file type for external charge data")
    end select
    if (iErr /= 0) then
      call error("Error during reading external file of charge data")
    end if
    if (nAtomInFile > nAtom) then
      call error("External charge file has more atoms than are present in the system")
    end if

    if (iSpin /= nSpin) then
      write(stdout, *) iSpin
      call error("Incorrect number of spins in restart file")
    end if

    qq(:,:,:) = 0.0_dp

    do iSpin = 1, nSpin
      do iAtom = 1, nAtomInFile
        nOrb = orb%nOrbAtom(iAtom)
        if (tReadAscii) then
          read (file, *, iostat=iErr) (qq(iOrb, iAtom, iSpin), iOrb = 1,nOrb)
        else
          read (file, iostat=iErr) (qq(iOrb, iAtom, iSpin), iOrb = 1,nOrb)
        end if
        if (iErr /= 0) then
          write (error_string, *) "Failure to read file of external charges"
          call error(error_string)
        end if
      end do
    end do

    if (any(abs(CheckSum(:) - sum(sum(qq(:,:nAtomInFile,:),dim=1),dim=1))>elecTolMax))then
      call error("Error during reading external file of charge data - internal checksum failure,&
          & probably a damaged file")
    end if
    sumQ = sum(qq(:,:,1))
    if (present(nEl)) then
      if (abs(nEl - sumQ) >= 1e-3_dp) then
        write(error_string, 99000) sumQ, nEl
99000   format ('External file of charges has a total charge:', F18.6,', instead of ',F18.6)
        call error(error_string)
      end if
    end if
    if (present(magnetisation)) then
      sumQ = sum(qq(:,:,2))
      if (abs(sumQ - magnetisation) >= 1e-3_dp) then
        write(error_string, 99010) sumQ, magnetisation
99010   format ('External file of charges has a total magnetisation:', F18.6,', instead of ',F18.6)
        call error(error_string)
      end if
    end if

    if (isMultipolar) then

      if (tReadAscii) then
        read(file, "(2I2)", iostat=iErr) nDipole, nQuadrupole
      else
        read(file, iostat=iErr)  nDipole, nQuadrupole
      end if
      if (allocated(multipoles%dipoleAtom)) then
        if (nDipole /= size(multipoles%dipoleAtom, dim=1)) then
          call error("Mismatch between expected and actual dipoles in charge restart")
        end if
      else
        if (nDipole /= 0) then
          call error("Mismatch between expected and actual dipoles in charge restart")
        end if
      end if
      if (allocated(multipoles%dipoleAtom)) then
        multipoles%dipoleAtom(:,:,:) = 0.0_dp
      end if
      if (allocated(multipoles%quadrupoleAtom)) then
        if (nQuadrupole /= size(multipoles%quadrupoleAtom, dim=1)) then
          call error("Mismatch between expected and actual quadrupoles in charge restart")
        end if
      else
        if (nQuadrupole /= 0) then
          call error("Mismatch between expected and actual quadrupoles in charge restart")
        end if
      end if
      if (allocated(multipoles%quadrupoleAtom)) then
        multipoles%quadrupoleAtom(:,:,:) = 0.0_dp
      end if
      if (tReadAscii) then
        do iSpin = 1, nSpin
          do iAtom = 1, nAtom
            do ii = 1, nDipole
              read (file, *, iostat=iErr) multipoles%dipoleAtom(ii,iAtom,iSpin)
              if (iErr /= 0) then
                write (error_string, *) "Failure to read file for atomic dipoles"
                call error(error_string)
              end if
            end do
          end do
        end do
        do iSpin = 1, nSpin
          do iAtom = 1, nAtom
            do ii = 1, nQuadrupole
              read (file, *, iostat=iErr) multipoles%quadrupoleAtom(ii,iAtom,iSpin)
              if (iErr /= 0) then
                write (error_string, *) "Failure to read file for atomic quadrupoles"
                call error(error_string)
              end if
            end do
          end do
        end do
      else
        if (nDipole > 0) then
          read(file, iostat=iErr)multipoles%dipoleAtom(:nDipole,:nAtom,:nSpin)
          if (iErr /= 0) then
            write(error_string, *)"Failure to read external dipoles from file for atomic multipoles"
            call error(error_string)
          end if
        end if
        if (nQuadrupole > 0) then
          read(file, iostat=iErr)multipoles%quadrupoleAtom(:nQuadrupole,:nAtom,:nSpin)
          if (iErr /= 0) then
            write(error_string, *) "Failure to read external quadrupoles from file for atomic&
                & multipoles"
            call error(error_string)
          end if
        end if
      end if
    end if

    if (tBlock) then
      qBlock(:,:,:,:) = 0.0_dp
      if (tBlockPresent) then
        do iSpin = 1, nSpin
          do iAtom = 1, nAtomInFile
            nOrb = orb%nOrbAtom(iAtom)
            do ii = 1, nOrb
              if (tReadAscii) then
                read (file, *, iostat=iErr) qBlock(1:nOrb, ii ,iAtom, iSpin)
              else
                read (file, iostat=iErr) qBlock(1:nOrb, ii ,iAtom, iSpin)
              end if
              if (iErr /= 0) then
                write (error_string, *) "Failure to read file for external block charges"
                call error(error_string)
              end if
            end do
          end do
        end do
      end if
      do iSpin = 1, nSpin
        do iAtom = 1, nAtomInFile
          nOrb = orb%nOrbAtom(iAtom)
          do ii = 1, nOrb
            qBlock(ii, ii ,iAtom, iSpin) = qq(ii ,iAtom, iSpin)
          end do
        end do
      end do
    end if

    if (tiBlock) then
      qiBlock(:,:,:,:) = 0.0_dp
      if (tiBlockPresent) then
        do iSpin = 1, nSpin
          do iAtom = 1, nAtomInFile
            nOrb = orb%nOrbAtom(iAtom)
            do ii = 1, nOrb
              if (tReadAscii) then
                read (file, *, iostat=iErr) qiBlock(1:nOrb, ii ,iAtom, iSpin)
              else
                read (file, iostat=iErr) qiBlock(1:nOrb, ii ,iAtom, iSpin)
              end if
              if (iErr /= 0) then
                write (error_string, *) "Failure to read file for external imaginary block charges"
                call error(error_string)
              end if
            end do
          end do
        end do
      end if
    end if
    ! need a checksum here

    if (tRho) then
      deltaRho(:) = 0.0_dp
      if (tRhoPresent) then
        do ii = 1, size(deltaRho)
          if (tReadAscii) then
            read (file, *, iostat=iErr) deltaRho(ii)
          else
            read (file, iostat=iErr) deltaRho(ii)
          end if
          if (iErr /= 0) then
            write (error_string, *) "Failure to read file for external imaginary block charges"
            call error(error_string)
          end if
        end do
      end if
    end if

    close(file)

  end subroutine initQFromFile


  !> Write the current charges to an external file
  subroutine writeQToFile(qq, fileName, tWriteAscii, orb, qBlock, qiBlock, deltaRhoIn,&
      & nAtInCentralRegion, multipoles)

    !> Array containing the charges
    real(dp), intent(in) :: qq(:,:,:)

    !> Name of the file to write the charges
    character(*), intent(in) :: fileName

    !> Write in a ascii format (T) or binary (F)
    logical, intent(in) :: tWriteAscii

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> block Mulliken population for LDA+U etc.
    real(dp), intent(in), allocatable :: qBlock(:,:,:,:)

    !> block Mulliken imagninary population for LDA+U and L.S
    real(dp), intent(in), allocatable :: qiBlock(:,:,:,:)

    !> Full density matrix with on-diagonal adjustment
    real(dp), intent(in), allocatable :: deltaRhoIn(:)

    !> Number of atoms in central region (atoms outside this will have charges suplied from
    !> elsewhere)
    integer, intent(in) :: nAtInCentralRegion

    !> Atomic multipoles, if relevant
    type(TMultipole), intent(in), optional :: multipoles

    character(len=120) :: error_string

    integer :: nAtom, nOrb, nSpin, nDipole, nQuadrupole
    integer :: iAtom, iOrb, iSpin, ii
    integer :: iErr, fd
    logical :: tqBlock, tqiBlock, tRho

    nAtom = nAtInCentralRegion
    nSpin = size(qq, dim=3)

    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin ==4)
    @:ASSERT(size(qq, dim=1) >= orb%mOrb)
    @:ASSERT(size(qq, dim=2) >= nAtInCentralRegion)

    tqBlock = allocated(qBlock)
    tqiBlock = allocated(qiBlock)
    tRho = allocated(deltaRhoIn)

  #:block DEBUG_CODE

    if (tqBlock) then
      @:ASSERT(all(shape(qBlock) >= [orb%mOrb,orb%mOrb,nAtom,nSpin]))
    end if

    if (tqiBlock) then
      @:ASSERT(allocated(qBlock))
      @:ASSERT(all(shape(qiBlock) == shape(qBlock)))
    end if

    if (tRho) then
      @:ASSERT(size(deltaRhoIn) == orb%nOrb*orb%nOrb*nSpin)
    end if

  #:endblock DEBUG_CODE

    if (tWriteAscii) then
      open(newunit=fd, file=trim(fileName)//'.dat', position="rewind", status="replace")
      write(fd, *, iostat=iErr) restartFormat
    else
      open(newunit=fd, file=trim(fileName)//'.bin', position="rewind", status="replace",&
          & form="unformatted")
      write(fd, iostat=iErr) restartFormat
    end if

    if (iErr /= 0) then
      write(error_string, *) "Failure to write file for external charges"
      call error(error_string)
    end if

    if (tWriteAscii) then
      write(fd, *, iostat=iErr) tqBlock, tqiBlock, tRho, present(multipoles), nAtom, nSpin,&
          & sum(sum(qq(:,:nAtom,:), dim=1), dim=1)
    else
      write(fd, iostat=iErr) tqBlock, tqiBlock, tRho, present(multipoles), nAtom,&
          & nSpin, sum(sum(qq(:,:nAtom,:), dim=1), dim=1)
    end if

    if (iErr /= 0) then
      write(error_string, *) "Failure to write file for external charges"
      call error(error_string)
    end if

    do iSpin = 1, nSpin
      do iAtom = 1, nAtom
        nOrb = orb%nOrbAtom(iAtom)
        if (tWriteAscii) then
          write(fd, *, iostat=iErr) (qq(iOrb, iAtom, iSpin), iOrb = 1, nOrb)
        else
          write(fd, iostat=iErr) (qq(iOrb, iAtom, iSpin), iOrb = 1, nOrb)
        end if
        if (iErr /= 0) then
          write(error_string, *) "Failure to write file for external charges"
          call error(error_string)
        end if
      end do
    end do

    if (present(multipoles)) then
      if (allocated(multipoles%dipoleAtom)) then
        nDipole = size(multipoles%dipoleAtom, dim=1)
        @:ASSERT(size(multipoles%dipoleAtom, dim=2) >= nAtom)
        @:ASSERT(nDipole > 0 .and. nDipole <= 3)
      else
        nDipole = 0
      end if
      if (allocated(multipoles%quadrupoleAtom)) then
        nQuadrupole = size(multipoles%quadrupoleAtom, dim=1)
        @:ASSERT(size(multipoles%quadrupoleAtom, dim=2) >= nAtom)
        @:ASSERT(nQuadrupole > 0 .and. nQuadrupole <= 6)
      else
        nQuadrupole =  0
      end if
      if (tWriteAscii) then
        write(fd, "(2I2)", iostat=iErr) nDipole, nQuadrupole
      else
        write(fd, iostat=iErr) nDipole, nQuadrupole
      end if
      if (iErr /= 0) then
        write(error_string, *) "Failure to write file for number of atomic multipoles"
        call error(error_string)
      end if
      if (tWriteAscii) then
        do iSpin = 1, nSpin
          do iAtom = 1, nAtom
            do ii = 1, nDipole
              write(fd, *, iostat=iErr) multipoles%dipoleAtom(ii, iAtom, iSpin)
              if (iErr /= 0) then
                write(error_string, "(A,I0)") "Failure to write dipole for atom ",iAtom
                call error(error_string)
              end if
            end do
          end do
        end do
        do iSpin = 1, nSpin
          do iAtom = 1, nAtom
            do ii = 1, nQuadrupole
              write(fd, "(3E20.12)", iostat=iErr) multipoles%quadrupoleAtom(ii, iAtom, iSpin)
              if (iErr /= 0) then
                write(error_string, "(A,I0)") "Failure to write quadrupole for atom ",iAtom
                call error(error_string)
              end if
            end do
          end do
        end do
      else
        if (nDipole > 0) then
          write(fd, iostat=iErr)multipoles%dipoleAtom(:nDipole,:nAtom,:nSpin)
          if (iErr /= 0) then
            write(error_string, *) "Failure to write atom dipoles"
            call error(error_string)
          end if
        end if
        if (nQuadrupole > 0) then
          write(fd, iostat=iErr)multipoles%quadrupoleAtom(:nQuadrupole,:nAtom,:nSpin)
          if (iErr /= 0) then
            write(error_string, *) "Failure to write atom quadrupoles"
            call error(error_string)
          end if
        end if
      end if
    end if

    if (tqBlock) then
      do iSpin = 1, nSpin
        do iAtom = 1, nAtom
          nOrb = orb%nOrbAtom(iAtom)
          do ii = 1, nOrb
            if (tWriteAscii) then
              write(fd, *, iostat=iErr) qBlock(1:nOrb, ii ,iAtom, iSpin)
            else
              write(fd, iostat=iErr) qBlock(1:nOrb, ii ,iAtom, iSpin)
            end if
            if (iErr /= 0) then
              write(error_string, *) "Failure to write file for external block charges"
              call error(error_string)
            end if
          end do
        end do
      end do
    end if

    if (tqiBlock) then
      do iSpin = 1, nSpin
        do iAtom = 1, nAtom
          nOrb = orb%nOrbAtom(iAtom)
          do ii = 1, nOrb
            if (tWriteAscii) then
              write(fd, *, iostat=iErr) qiBlock(1:nOrb, ii ,iAtom, iSpin)
            else
              write(fd, iostat=iErr) qiBlock(1:nOrb, ii ,iAtom, iSpin)
            end if
            if (iErr /= 0) then
              write(error_string, *) "Failure to write file for external block imaginary charges"
              call error(error_string)
            end if
          end do
        end do
      end do
    end if

    if (tRho) then
      do ii = 1, size(deltaRhoIn)
        if (tWriteAscii) then
          write(fd, *, iostat=iErr) deltaRhoIn(ii)
        else
          write(fd, iostat=iErr) deltaRhoIn(ii)
        end if
        if (iErr /= 0) then
          write(error_string, *) "Failure to write file for external density matrix"
          call error(error_string)
        end if
      end do
    end if

    close(fd)

  end subroutine writeQToFile

end module dftbp_dftb_sccinit
