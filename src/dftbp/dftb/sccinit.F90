!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for initializing SCC part of the calculation.
module dftbp_dftb_sccinit
  use dftbp_common_accuracy, only : dp, elecTolMax
  use dftbp_common_file, only : TFileDescr, openFile, closeFile
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_status, only : TStatus
  use dftbp_io_message, only : error
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_multipole, only : TMultipole
  use dftbp_dftb_densitymatrix, onLy : TDensityMatrix
  use dftbp_dftb_hybridxc, only : checkSupercellFoldingMatrix, hybridXcAlgo
  use dftbp_dftb_periodic, only : getSuperSampling
  implicit none

  private
  public :: initQFromAtomChrg, initQFromShellChrg, initQFromFile, writeQToFile
  public :: initQFromUsrChrg

  !> Version number for restart format, please increment if you change the interface.
  integer, parameter :: restartFormat = 8

contains


  !> Initialise the charge vector from the reference atomic charges.
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
  !! appropriate for the neutral spin unpolarised atom reference system that DFTB assumes for
  !! SCC/spin extensions.
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


  !> Initialise the charge vector from a named external file.
  !! Checks that the total charge matches the expected value for the calculation.
  !! Should test the input, if the number of orbital charges per atom matches the number from the
  !! angular momentum.
  subroutine initQFromFile(qq, fileName, tReadAscii, orb, qBlock, qiBlock, densityMatrix, tRealHS,&
      & errStatus, magnetisation, nEl, hybridXcAlg, coeffsAndShifts, multipoles)

    !> The charges per lm,atom,spin
    real(dp), intent(out) :: qq(:,:,:)

    !> The external file of charges for the orbitals, currently stored with each line containing the
    !! per-orbital charges in order of increasing m and l. Alternating lines give the spin case (if
    !! present)
    character(*), intent(in) :: fileName

    !> Should charges be read as ascii (cross platform, but potentially lower reproducibility) or
    !! binary files
    logical, intent(in) :: tReadAscii

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Block Mulliken population for LDA+U etc
    real(dp), intent(inout), allocatable :: qBlock(:,:,:,:)

    !> Block Mulliken imagninary population for LDA+U and L.S
    real(dp), intent(inout), allocatable :: qiBlock(:,:,:,:)

    !> Holds real and complex delta density matrices
    type(TDensityMatrix), intent(inout) :: densityMatrix

    !> True, if overlap and Hamiltonian are real-valued
    logical, intent(in) :: tRealHS

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !> Magnetisation checksum for regular spin polarization total magnetic moment
    real(dp), intent(in), optional :: magnetisation

    !> Nr. of electrons for each spin channel
    real(dp), intent(in), optional :: nEl

    !> Hybrid Hamiltonian construction algorithm
    integer, intent(in), optional :: hybridXcAlg

    !> Coefficients of the lattice vectors in the linear combination for the super lattice vectors
    !! (should be integer values) and shift of the grid along the three small reciprocal lattice
    !! vectors (between 0.0 and 1.0)
    real(dp), intent(inout), optional :: coeffsAndShifts(:,:)

    !> Atomic multipoles, if relevant
    type(TMultipole), intent(inout), optional :: multipoles

    !! Diagonal elements of supercell folding matrix
    integer :: supercellFoldingDiag(3)

    !! Nr. of orbitals / atoms / spin channels
    integer :: nOrb, nAtom, nSpin

    !! Error returned by the io commands
    integer :: iErr

    !! Total charge is present at the top of the file
    real(dp) :: checkSum(size(qq, dim=3))

    integer :: iOrb, iAtom, iSpin, ii, jj, kk, nAtomInFile, nDipole, nQuadrupole

    type(TFileDescr) :: file

    integer :: fileFormat
    real(dp) :: sumQ

    !! Requested to be re-loaded
    logical :: tBlock, tiBlock, tRho, tKpointInfo, isMultipolar

    !! Present in the file itself
    logical :: tBlockPresent, tiBlockPresent, tRhoPresent, tKpointInfoPresent

    character(len=120) :: error_string

    nAtom = size(qq, dim=2)
    nSpin = size(qq, dim=3)

    tBlock = allocated(qBlock)
    tiBlock = allocated(qiBlock)
    tRho = allocated(densityMatrix%deltaRhoIn) .or. allocated(densityMatrix%deltaRhoInCplx)&
        & .or. allocated(densityMatrix%deltaRhoInCplxHS)
    tKpointInfo = present(coeffsAndShifts)

    if (tRho .and. (.not. present(hybridXcAlg))) then
      call error("Missing hybrid xc-functional algorithm information.")
    end if

    @:ASSERT(size(qq, dim=1) == orb%mOrb)
    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)

  #:block DEBUG_CODE

    if (present(magnetisation)) then
      @:ASSERT(nSpin==2)
    end if

    if (tBlock) then
      @:ASSERT(all(shape(qBlock) == [orb%mOrb, orb%mOrb, nAtom, nSpin]))
    end if

    if (tiBlock) then
      @:ASSERT(tBlock)
      @:ASSERT(all(shape(qiBlock) == shape(qBlock)))
    end if

    if (tKpointInfo) then
      @:ASSERT(all(shape(coeffsAndShifts) == [3, 4]))
    end if

    if (allocated(densityMatrix%deltaRhoIn)) then
      @:ASSERT(size(densityMatrix%deltaRhoIn) == orb%nOrb * orb%nOrb * nSpin)
    end if

  #:endblock DEBUG_CODE

    if (tReadAscii) then
      call openFile(file, trim(fileName) // '.dat', mode="r", iostat=iErr)
    else
      call openFile(file, file=trim(fileName) // '.bin', mode="rb", iostat=iErr)
    end if
    if (iErr /= 0) then
      write(error_string, *) "Failure to open external file of charge data"
      call error(error_string)
    end if
    rewind(file%unit)

    if (tReadAscii) then
      read(file%unit, *, iostat=iErr) fileFormat
    else
      read(file%unit, iostat=iErr) fileFormat
    end if
    if (fileFormat < 7 .and. (.not. tRealHS) .and. tRho) then
      call error("Restart file belongs to a version that does not support general hybrid&
          & xc-functionals.")
    end if
    if (iErr /= 0) then
      call error("Error during reading external file of charge data")
    end if
    isMultipolar = .false.
    tKpointInfoPresent = .false.
    select case(fileFormat)
    case(4)
      if (tReadAscii) then
        read(file%unit, *, iostat=iErr) tBlockPresent, tiBlockPresent, tRhoPresent, iSpin, checkSum
      else
        read(file%unit, iostat=iErr) tBlockPresent, tiBlockPresent, tRhoPresent, iSpin, checkSum
      end if
      nAtomInFile = nAtom
    case(5)
      if (tReadAscii) then
        read(file%unit, *, iostat=iErr)tBlockPresent, tiBlockPresent, tRhoPresent, nAtomInFile,&
            & iSpin, checkSum
      else
        read(file%unit, iostat=iErr)tBlockPresent, tiBlockPresent, tRhoPresent, nAtomInFile, iSpin,&
            & checkSum
      end if
    case(6)
      if (tReadAscii) then
        read(file%unit, *, iostat=iErr)tBlockPresent, tiBlockPresent, tRhoPresent, isMultipolar,&
            & nAtomInFile, iSpin, checkSum
      else
        read(file%unit, iostat=iErr)tBlockPresent, tiBlockPresent, tRhoPresent, isMultipolar,&
            & nAtomInFile, iSpin, checkSum
      end if
    case(7, 8)
      if (tReadAscii) then
        read(file%unit, *, iostat=iErr) tBlockPresent, tiBlockPresent, tRhoPresent,&
            & tKpointInfoPresent, isMultipolar, nAtomInFile, iSpin, checkSum
      else
        read(file%unit, iostat=iErr) tBlockPresent, tiBlockPresent, tRhoPresent,&
            & tKpointInfoPresent, isMultipolar, nAtomInFile, iSpin, checkSum
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
          read(file%unit, *, iostat=iErr) (qq(iOrb, iAtom, iSpin), iOrb = 1,nOrb)
        else
          read(file%unit, iostat=iErr) (qq(iOrb, iAtom, iSpin), iOrb = 1,nOrb)
        end if
        if (iErr /= 0) then
          write(error_string, *) "Failure to read file of external charges"
          call error(error_string)
        end if
      end do
    end do

    if (any(abs(checkSum(:) - sum(sum(qq(:,:nAtomInFile,:),dim=1),dim=1)) > elecTolMax))then
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
        read(file%unit, "(2I2)", iostat=iErr) nDipole, nQuadrupole
      else
        read(file%unit, iostat=iErr)  nDipole, nQuadrupole
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
              read(file%unit, *, iostat=iErr) multipoles%dipoleAtom(ii,iAtom,iSpin)
              if (iErr /= 0) then
                write(error_string, *) "Failure to read file for atomic dipoles"
                call error(error_string)
              end if
            end do
          end do
        end do
        do iSpin = 1, nSpin
          do iAtom = 1, nAtom
            do ii = 1, nQuadrupole
              read(file%unit, *, iostat=iErr) multipoles%quadrupoleAtom(ii,iAtom,iSpin)
              if (iErr /= 0) then
                write(error_string, *) "Failure to read file for atomic quadrupoles"
                call error(error_string)
              end if
            end do
          end do
        end do
      else
        if (nDipole > 0) then
          read(file%unit, iostat=iErr)multipoles%dipoleAtom(:nDipole,:nAtom,:nSpin)
          if (iErr /= 0) then
            write(error_string, *)&
                & "Failure to read external dipoles from file for atomic multipoles"
            call error(error_string)
          end if
        end if
        if (nQuadrupole > 0) then
          read(file%unit, iostat=iErr)multipoles%quadrupoleAtom(:nQuadrupole,:nAtom,:nSpin)
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
                read(file%unit, *, iostat=iErr) qBlock(1:nOrb, ii, iAtom, iSpin)
              else
                read(file%unit, iostat=iErr) qBlock(1:nOrb, ii, iAtom, iSpin)
              end if
              if (iErr /= 0) then
                write(error_string, *) "Failure to read file for external block charges"
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
                read(file%unit, *, iostat=iErr) qiBlock(1:nOrb, ii, iAtom, iSpin)
              else
                read(file%unit, iostat=iErr) qiBlock(1:nOrb, ii, iAtom, iSpin)
              end if
              if (iErr /= 0) then
                write(error_string, *) "Failure to read file for external imaginary block charges"
                call error(error_string)
              end if
            end do
          end do
        end do
      end if
    end if
    ! need a checksum here

    if (allocated(densityMatrix%deltaRhoIn)) then
      densityMatrix%deltaRhoIn(:,:,:) = 0.0_dp
    end if
    ! In this case the size of deltaRho couldn't be known in advance, therefore deallocating
    if (allocated(densityMatrix%deltaRhoInCplxHS) .and. tKpointInfo .and. tKpointInfoPresent) then
      deallocate(densityMatrix%deltaRhoInCplxHS)
    end if
    if (allocated(densityMatrix%deltaRhoInCplx) .and. tKpointInfo .and. tKpointInfoPresent) then
      deallocate(densityMatrix%deltaRhoInCplx)
    end if

    if (tRho) then
      if (tKpointInfo) then
        coeffsAndShifts(:,:) = 0.0_dp
        if (tKpointInfoPresent) then
          if (tReadAscii) then
            read(file%unit, *, iostat=iErr) coeffsAndShifts
          else
            read(file%unit, iostat=iErr) coeffsAndShifts
          end if
          call checkSupercellFoldingMatrix(coeffsAndShifts, errStatus,&
              & supercellFoldingDiagOut=supercellFoldingDiag)
          if (hybridXcAlg == hybridXcAlgo%matrixBased) then
            call getSuperSampling(coeffsAndShifts(:,1:3), modulo(coeffsAndShifts(:,4), 1.0_dp),&
                & densityMatrix%kPointPrime, densityMatrix%kWeightPrime, reduceByInversion=.true.)
            allocate(densityMatrix%deltaRhoInCplx(orb%nOrb, orb%nOrb,&
                & size(densityMatrix%kPointPrime, dim=2) * nSpin))
          else
            allocate(densityMatrix%deltaRhoInCplxHS(orb%nOrb, orb%nOrb, supercellFoldingDiag(1),&
                & supercellFoldingDiag(2), supercellFoldingDiag(3), nSpin))
          end if
        end if
      end if
      if (tRhoPresent) then
        if (.not. tRealHS) then
          ! general k-point case
          if (hybridXcAlg == hybridXcAlgo%matrixBased) then
            ! matrix-multiplication based algorithm
            if (tReadAscii) then
              read(file%unit, *, iostat=iErr) densityMatrix%deltaRhoInCplx
            else
              read(file%unit, iostat=iErr) densityMatrix%deltaRhoInCplx
            end if
          else
            ! neighbor-list based algorithm
            if (tReadAscii) then
              read(file%unit, *, iostat=iErr) densityMatrix%deltaRhoInCplxHS
            else
              read(file%unit, iostat=iErr) densityMatrix%deltaRhoInCplxHS
            end if
          end if
        else
          ! cluster/Gamma-only case
          do ii = 1, size(densityMatrix%deltaRhoIn, dim=3)
            do jj = 1, size(densityMatrix%deltaRhoIn, dim=2)
              do kk = 1, size(densityMatrix%deltaRhoIn, dim=1)
                if (tReadAscii) then
                  read(file%unit, *, iostat=iErr) densityMatrix%deltaRhoIn(kk, jj, ii)
                else
                  read(file%unit, iostat=iErr) densityMatrix%deltaRhoIn(kk, jj, ii)
                end if
              end do
            end do
          end do
        end if
        if (iErr /= 0) then
          write(error_string, *) 'Failure to read file for delta density matrix.'
          call error(error_string)
        end if
      end if
    end if

    call closeFile(file)

  end subroutine initQFromFile


  !> Write the current charges to an external file
  subroutine writeQToFile(qq, fileName, tWriteAscii, orb, qBlock, qiBlock, densityMatrix, tRealHS,&
      & nAtInCentralRegion, hybridXcAlg, coeffsAndShifts, multipoles)

    !> Array containing the charges
    real(dp), intent(in) :: qq(:,:,:)

    !> Name of the file to write the charges
    character(len=*), intent(in) :: fileName

    !> Write in a ascii format (T) or binary (F)
    logical, intent(in) :: tWriteAscii

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Block Mulliken population for LDA+U etc.
    real(dp), intent(in), allocatable :: qBlock(:,:,:,:)

    !> Block Mulliken imagninary population for LDA+U and L.S
    real(dp), intent(in), allocatable :: qiBlock(:,:,:,:)

    !> Holds real and complex delta density matrices
    type(TDensityMatrix), intent(in) :: densityMatrix

    !> True, if overlap and Hamiltonian are real-valued
    logical, intent(in) :: tRealHS

    !> Number of atoms in central region (atoms outside this will have charges suplied from
    !> elsewhere)
    integer, intent(in) :: nAtInCentralRegion

    !> Hybrid Hamiltonian construction algorithm
    integer, intent(in), optional :: hybridXcAlg

    !> Coefficients of the lattice vectors in the linear combination for the super lattice vectors
    !! (should be integer values) and shift of the grid along the three small reciprocal lattice
    !! vectors (between 0.0 and 1.0)
    real(dp), intent(in), optional :: coeffsAndShifts(:,:)

    !> Atomic multipoles, if relevant
    type(TMultipole), intent(in), optional :: multipoles

    character(len=120) :: error_string

    integer :: nAtom, nOrb, nSpin, nDipole, nQuadrupole
    integer :: iAtom, iOrb, iSpin, ii, jj, kk
    integer :: iErr
    logical :: tqBlock, tqiBlock, tRho
    type(TFileDescr) :: fd

    nAtom = nAtInCentralRegion
    nSpin = size(qq, dim=3)

    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin ==4)
    @:ASSERT(size(qq, dim=1) >= orb%mOrb)
    @:ASSERT(size(qq, dim=2) >= nAtInCentralRegion)

    if (present(coeffsAndShifts)) then
      ! corresponds to Hybrid-DFTB case
      @:ASSERT(all(shape(coeffsAndShifts) == [3, 4]))
    end if

    tqBlock = allocated(qBlock)
    tqiBlock = allocated(qiBlock)
    tRho = allocated(densityMatrix%deltaRhoIn) .or. allocated(densityMatrix%deltaRhoInCplx)&
        & .or. allocated(densityMatrix%deltaRhoInCplxHS)

    if (tRho .and. (.not. allocated(densityMatrix%deltaRhoIn))&
        & .and. (.not. present(coeffsAndShifts))) then
      call error("Failure while writing restart file: This appears to be a cluster/Gamma-only&
          & calculations, but the associated density matrix is missing.")
    end if

    if (present(coeffsAndShifts) .and. (.not. tRealHS)&
        & .and. (.not. (allocated(densityMatrix%deltaRhoInCplx)&
        & .or. allocated(densityMatrix%deltaRhoInCplxHS)))) then
      call error("Failure while writing restart file: Supercell folding coefficients and shifts&
          & present, but the associated density matrix is missing.")
    end if

  #:block DEBUG_CODE

    if (tqBlock) then
      @:ASSERT(all(shape(qBlock) >= [orb%mOrb,orb%mOrb,nAtom,nSpin]))
    end if

    if (tqiBlock) then
      @:ASSERT(allocated(qBlock))
      @:ASSERT(all(shape(qiBlock) == shape(qBlock)))
    end if

  #:endblock DEBUG_CODE

    if (tWriteAscii) then
      call openFile(fd, trim(fileName) // '.dat', mode="w")
      write(fd%unit, *, iostat=iErr) restartFormat
    else
      call openFile(fd, trim(fileName) // '.bin', mode="wb")
      write(fd%unit, iostat=iErr) restartFormat
    end if

    if (iErr /= 0) then
      write(error_string, *) "Failure to write file for external charges"
      call error(error_string)
    end if

    if (tWriteAscii) then
      write(fd%unit, *, iostat=iErr) tqBlock, tqiBlock, tRho, present(coeffsAndShifts),&
          & present(multipoles), nAtom, nSpin, sum(sum(qq(:,:nAtom,:), dim=1), dim=1)
    else
      write(fd%unit, iostat=iErr) tqBlock, tqiBlock, tRho, present(coeffsAndShifts),&
          & present(multipoles), nAtom, nSpin, sum(sum(qq(:,:nAtom,:), dim=1), dim=1)
    end if

    if (iErr /= 0) then
      write(error_string, *) "Failure to write file for external charges"
      call error(error_string)
    end if

    do iSpin = 1, nSpin
      do iAtom = 1, nAtom
        nOrb = orb%nOrbAtom(iAtom)
        if (tWriteAscii) then
          write(fd%unit, *, iostat=iErr) (qq(iOrb, iAtom, iSpin), iOrb = 1, nOrb)
        else
          write(fd%unit, iostat=iErr) (qq(iOrb, iAtom, iSpin), iOrb = 1, nOrb)
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
        write(fd%unit, "(2I2)", iostat=iErr) nDipole, nQuadrupole
      else
        write(fd%unit, iostat=iErr) nDipole, nQuadrupole
      end if
      if (iErr /= 0) then
        write(error_string, *) "Failure to write file for number of atomic multipoles"
        call error(error_string)
      end if
      if (tWriteAscii) then
        do iSpin = 1, nSpin
          do iAtom = 1, nAtom
            do ii = 1, nDipole
              write(fd%unit, *, iostat=iErr) multipoles%dipoleAtom(ii, iAtom, iSpin)
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
              write(fd%unit, "(3E20.12)", iostat=iErr) multipoles%quadrupoleAtom(ii, iAtom, iSpin)
              if (iErr /= 0) then
                write(error_string, "(A,I0)") "Failure to write quadrupole for atom ",iAtom
                call error(error_string)
              end if
            end do
          end do
        end do
      else
        if (nDipole > 0) then
          write(fd%unit, iostat=iErr)multipoles%dipoleAtom(:nDipole,:nAtom,:nSpin)
          if (iErr /= 0) then
            write(error_string, *) "Failure to write atom dipoles"
            call error(error_string)
          end if
        end if
        if (nQuadrupole > 0) then
          write(fd%unit, iostat=iErr)multipoles%quadrupoleAtom(:nQuadrupole,:nAtom,:nSpin)
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
              write(fd%unit, *, iostat=iErr) qBlock(1:nOrb, ii, iAtom, iSpin)
            else
              write(fd%unit, iostat=iErr) qBlock(1:nOrb, ii, iAtom, iSpin)
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
              write(fd%unit, *, iostat=iErr) qiBlock(1:nOrb, ii, iAtom, iSpin)
            else
              write(fd%unit, iostat=iErr) qiBlock(1:nOrb, ii, iAtom, iSpin)
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
      ! Write k-point set information to file, if CAM calculation with k-points is present
      if (present(coeffsAndShifts)) then
        if (tWriteAscii) then
          write(fd%unit, *, iostat=iErr) coeffsAndShifts
        else
          write(fd%unit, iostat=iErr) coeffsAndShifts
        end if
      end if

      if (tRealHS) then
        ! cluster/Gamma-only case
        do ii = 1, size(densityMatrix%deltaRhoIn, dim=3)
          do jj = 1, size(densityMatrix%deltaRhoIn, dim=2)
            do kk = 1, size(densityMatrix%deltaRhoIn, dim=1)
              if (tWriteAscii) then
                write(fd%unit, *, iostat=iErr) densityMatrix%deltaRhoIn(kk, jj, ii)
              else
                write(fd%unit, iostat=iErr) densityMatrix%deltaRhoIn(kk, jj, ii)
              end if
            end do
          end do
        end do
      else
        ! general k-point case
        if (hybridXcAlg == hybridXcAlgo%matrixBased) then
          ! matrix-multiplication based algorithm
          if (tWriteAscii) then
            write(fd%unit, *, iostat=iErr) densityMatrix%deltaRhoInCplx
          else
            write(fd%unit, iostat=iErr) densityMatrix%deltaRhoInCplx
          end if
        else
          ! neighbor-list based algorithm
          if (tWriteAscii) then
            write(fd%unit, *, iostat=iErr) densityMatrix%deltaRhoInCplxHS
          else
            write(fd%unit, iostat=iErr) densityMatrix%deltaRhoInCplxHS
          end if
        end if
      end if
      if (iErr /= 0) then
        write(error_string, *) "Failure to write file for external density matrix"
        call error(error_string)
      end if
    end if

    call closeFile(fd)

  end subroutine writeQToFile

end module dftbp_dftb_sccinit
