!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Various I/O routines for the main program.
module mainio
  use mpifx
  use globalenv
  use environment
  use assert
  use accuracy
  use constants
  use periodic
  use commontypes
  use fifo
  use sparse2dense
  use blasroutines
  use charmanip, only : i2c
  use linkedlist
  use taggedoutput
  use fileid
  use spin, only : qm2ud
  use energies
  use xmlf90
  use hsdutils, only : writeChildValue
  use mdintegrator, only : OMdIntegrator, state
  use formatout
  use sccinit, only : writeQToFile
  use message
  use ipisocket
  implicit none
  private

  public :: writeEigvecs, writeProjEigvecs, getH, SetEigVecsTxtOutput
  public :: initOutputFile, writeAutotestTag, writeResultsTag, writeDetailedXml, writeBandOut
  public :: writeHessianOut
  public :: writeDetailedOut1, writeDetailedOut2, writeDetailedOut3, writeDetailedOut4
  public :: writeDetailedOut5
  public :: writeMdOut1, writeMdOut2, writeMdOut3
  public :: writeCharges, writeEigenvectors, writeProjectedEigenvectors
  public :: writeCurrentGeometry, writeFinalDriverStatus
  public :: writeHSAndStop, writeHS
  public :: printGeoStepInfo, printSccHeader, printSccInfo, printEnergies, printVolume
  public :: printPressureAndFreeEnergy, printMaxForce, printMaxLatticeForce
  public :: printMdInfo
  public :: receiveGeometryFromSocket


  ! output file names

  !> Ground state eigenvectors in text format
  character(*), parameter :: eigvecOut = "eigenvec.out"

  !> Ground state eigenvectors in binary format
  character(*), parameter :: eigvecBin = "eigenvec.bin"

  !> Format string for energy second derivative matrix
  character(len=*), parameter :: formatHessian = '(4f16.10)'

  !> Atomic geometries format
  character(len=*), parameter :: formatGeoOut = "(I5, F16.8, F16.8, F16.8)"

  !> Format for a single value with units
  character(len=*), parameter :: format1U = "(A, ':', T32, F18.10, T51, A)"

  !> Format for two values with units
  character(len=*), parameter :: format2U = "(A, ':', T32, F18.10, T51, A, T54, F16.4, T71, A)"

  !> Format for a single value using exponential notation with units
  character(len=*), parameter :: format1Ue = "(A, ':', T37, E13.6, T51, A)"

  !> Format for two using exponential notation values with units
  character(len=*), parameter :: format2Ue = "(A, ':', T37, E13.6, T51, A, T57, E13.6, T71,A)"

  !> Format for mixed decimal and exponential values with units
  character(len=*), parameter :: format1U1e =&
      & "(' ', A, ':', T32, F18.10, T51, A, T57, E13.6, T71, A)"

  ! Private module variables (suffixed with "_" for clarity)

  !> Should eigenvectors be written as text format in addition to  binary data
  logical :: EigVecsAsTxt_ = .false.


  !> Routines to get eigenvectors out of storage/memory
  interface getH
    module procedure getHreal
    module procedure getHcmplx
  end interface getH


  !> write eigenvectors to disc
  interface writeEigvecs
    module procedure writeRealEigvecs
    module procedure writeCplxEigvecs
  end interface writeEigvecs


  !> write eigenvector projections onto defined regions
  interface writeProjEigvecs
    module procedure writeProjRealEigvecs
    module procedure writeProjCplxEigvecs
  end interface writeProjEigvecs

contains


  !> Sets internal logical flag which controls whether to write a txt file for eigenvectors.
  subroutine SetEigVecsTxtOutput(tTxtWrite)

    !> Is a txt file written out, as well as the binary data file?
    logical, intent(in) :: tTxtWrite

    EigVecsAsTxt_ = tTxtWrite
  end subroutine SetEigVecsTxtOutput


  !> Write the real eigenvectors into text and binary output files.
  subroutine writeRealEigvecs(fdEigvec, runId, nAtom, nSpin, neighlist, &
      &nNeighbor, iAtomStart, iPair, img2CentCell, orb, species, speciesName, &
      &over, HSqrReal, SSqrReal, storeEigvecs, fileName)

    !> Fileid (file not yet opened) to use.
    integer, intent(in) :: fdEigvec

    !> Id of the current program run.
    integer, intent(in) :: runId

    !> Nr. of atoms in the system.
    integer, intent(in) :: nAtom

    !> Nr. of spin channels.
    integer, intent(in) :: nSpin

    !> Neighbor list.
    type(TNeighborList), intent(in) :: neighlist

    !> Nr. of neighbors for SK-interaction.
    integer, intent(in) :: nNeighbor(:)

    !> Positions of atoms int the dense matrices.
    integer, intent(in) :: iAtomStart(:)

    !> Positions of interactions in the sparse matrices.
    integer, intent(in) :: iPair(:,:)

    !> Mapping of atoms into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Species.
    integer, intent(in) :: species(:)

    !> Name of the species.
    character(mc), intent(in) :: speciesName(:)

    !> Sparse overlap matrix.
    real(dp), intent(in) :: over(:)

    !> Square Hamiltonian (or work array)
    real(dp), intent(inout) :: HSqrReal(:,:,:)

    !> Work array.
    real(dp), intent(inout) :: SSqrReal(:,:)

    !> If present, Hamiltonian(s) are fetched from this storage into HSqrReal, instead of using
    !> whatever is already there.
    type(OFifoRealR2), intent(inout), optional :: storeEigvecs(:)

    !> optional alternative file pre-fix
    character(len=*), intent(in), optional :: fileName

    character(lc) :: tmpStr
    integer :: iSpin, iSpin2, iAtom, iSp1, iSh1, iOrb, ang
    integer :: ii, jj
    real(dp), allocatable :: rVecTemp(:)

    @:ASSERT(nSpin == 1 .or. nSpin == 2)

    close(fdEigvec) ! just to be on the safe side
    ! Write eigenvalues in binary form
    if (present(fileName)) then
      write (tmpStr, "(A,A)") trim(fileName), ".bin"
      open(fdEigvec, file=tmpStr, action="write", status="replace", &
          &position="rewind", form="unformatted")
    else
      open(fdEigvec, file=eigvecBin, action="write", status="replace", &
          &position="rewind", form="unformatted")
    end if
    write (fdEigVec) runId
    do iSpin = 1, nSpin
      call getH(iSpin, HSqrReal, iSpin2, storeEigvecs)
      do ii = 1, size(HSqrReal, dim=2)
        write (fdEigvec) HSqrReal(:,ii,iSpin2)
      end do
    end do
    close(fdEigvec)

    if (EigVecsAsTxt_) then
      ! Write eigenvalues (together with Mulliken populations) in text form
      if (present(fileName)) then
        write (tmpStr, "(A,A)") trim(fileName), ".out"
        open(fdEigvec, file=tmpStr, action="write", status="replace", &
            &position="rewind")
      else
        open(fdEigvec, file=eigvecOut, action="write", status="replace", &
            &position="rewind")
      end if
      allocate(rVecTemp(size(HSqrReal, dim=1)))
      call unpackHS(SSqrReal, over, neighlist%iNeighbor, nNeighbor, &
          &iAtomStart, iPair, img2CentCell)
      do iSpin = 1, nSpin
        call getH(iSpin, HSqrReal, iSpin2, storeEigvecs)
        do ii = 1, orb%nOrb
          call hemv(rVecTemp, SSqrReal, HSqrReal(:,ii,iSpin2))
          write(fdEigvec, "('Eigenvector:',I4,4X,'(',A,')'/)") ii, &
              &trim(spinName(iSpin2))
          jj = 0
          do iAtom = 1, nAtom
            iSp1 = species(iAtom)
            do iSh1 = 1, orb%nShell(iSp1)
              ang = orb%angShell(iSh1, iSp1)
              if (iSh1 == 1) then
                write(tmpStr, "(I5,1X,A2,2X,A1)") iAtom, speciesName(iSp1), &
                    &orbitalNames(ang+1)
              else
                write(tmpStr, "(10X,A1)") orbitalNames(ang+1)
              end if
              do iOrb = 1, 2 * ang + 1
                jj = jj + 1
                write(fdEigvec,"(A,I1,T15,F12.6,3X,F12.6)") trim(tmpStr),&
                    &iOrb, HSqrReal(jj, ii, iSpin2), &
                    & HSqrReal(jj, ii, iSpin2) * rVecTemp(jj)
              end do
            end do
            write (fdEigvec,*)
          end do
        end do
      end do

      close(fdEigvec)
    end if

  end subroutine writeRealEigvecs


  !> Write the complex eigenvectors into text and binary output files.
  subroutine writeCplxEigvecs(fdEigvec, runId, nAtom, nSpin, neighlist, &
      &nNeighbor, cellVec, iCellVec, iAtomStart, iPair, img2CentCell, orb, &
      &species, speciesName, over, kpoint, HSqrCplx, SSqrCplx, storeEigvecs, &
      & fileName)

    !> Fileid (file not yet opened) to use.
    integer, intent(in) :: fdEigvec

    !> Id of the current program run.
    integer, intent(in) :: runId

    !> Nr. of atoms in the system.
    integer, intent(in) :: nAtom

    !> Nr. of spin channels.
    integer, intent(in) :: nSpin

    !> Neighbor list.
    type(TNeighborList), intent(in) :: neighlist

    !> Nr. of neighbors for SK-interaction.
    integer, intent(in) :: nNeighbor(:)

    !> Cell vectors of shifted cells.
    real(dp), intent(in) :: cellVec(:,:)

    !> Cell vector index of every atom.
    integer, intent(in) :: iCellVec(:)

    !> Positions of atoms int the dense matrices.
    integer, intent(in) :: iAtomStart(:)

    !> Positions of interactions in the sparse matrices.
    integer, intent(in) :: iPair(:,:)

    !> Mapping of atoms into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Species.
    integer, intent(in) :: species(:)

    !> Name of the species.
    character(mc), intent(in) :: speciesName(:)

    !> Sparse overlap matrix.
    real(dp), intent(in) :: over(:)

    !> KPoints.
    real(dp), intent(in) :: kpoint(:,:)

    !> Square Hamiltonian (or work array)
    complex(dp), intent(inout) :: HSqrCplx(:,:,:,:)

    !> Work array.
    complex(dp), intent(inout) :: SSqrCplx(:,:)

    !> If present, Hamiltonian(s) are fetched from this storage into HSqrCplx, instead of using
    !> whatever is already there.
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecs(:)

    !> optional alternative file prefix, to appear as "fileName".bin
    character(len=*), intent(in), optional :: fileName

    character(lc) :: tmpStr
    integer :: iSpin, iSpin2, iAtom, iSp1, iSh1, iOrb, ang, iK, iK2, nK
    integer :: ii, jj, nOrb, nSpinChannel
    complex(dp), allocatable :: cVecTemp(:), work(:,:)

    nK = size(kPoint, dim=2)
    close(fdEigvec) ! just to be on the safe side
    if (present(fileName)) then
      write (tmpStr, "(A,A)") trim(fileName), ".bin"
      open(fdEigvec, file=tmpStr, action="write", status="replace", &
          &position="rewind", form="unformatted")
    else
      open(fdEigvec, file=eigvecBin, action="write", status="replace", &
          &position="rewind", form="unformatted")
    end if
    write (fdEigVec) runId
    if (nSpin == 4) then
      nSpinChannel = 1
    else
      nSpinChannel = nSpin
    end if
    do iSpin = 1, nSpinChannel
      do iK = 1, nK
        call getH(iSpin, iK, HSqrCplx, iSpin2, iK2, storeEigvecs)
        do ii = 1, size(HSqrCplx, dim=2)
          write (fdEigvec) HSqrCplx(:,ii,iK2, iSpin2)
        end do
      end do
    end do

    close(fdEigvec)

    if (EigVecsAsTxt_) then
      ! Write eigenvalues (together with Mulliken populations) in text form
      if (present(fileName)) then
        write (tmpStr, "(A,A)") trim(fileName), ".out"
        open(fdEigvec, file=tmpStr, action="write", status="replace", &
            &position="rewind")
      else
        open(fdEigvec, file=eigvecOut, action="write", status="replace", &
            &position="rewind")
      end if
      write (fdEigvec,"(A/)") "Coefficients and Mulliken populations of the &
          &atomic orbitals"

      if (nSpin == 4) then
        write(fdEigvec,"(A/)")"   Atom   Orb  up spin coefficients        &
            &down spin coefficients         charge      x           y           z"
        allocate(cVecTemp(size(HSqrCplx, dim=1)))
        nOrb = size(HSqrCplx, dim=1) / 2
        allocate(work(nOrb,nOrb))
        do iK = 1, nK
          call getH(1, iK, HSqrCplx, iSpin2, iK2, storeEigvecs)
          SSqrCplx = 0.0_dp
          work = 0.0_dp
          call unpackHS(work, over, kPoint(:,iK), &
              & neighlist%iNeighbor, nNeighbor, iCellVec, cellVec, &
              & iAtomStart, iPair, img2CentCell)
          SSqrCplx(:nOrb,:nOrb) = work
          SSqrCplx(nOrb+1:,nOrb+1:) = work
          do ii = 1, 2*orb%nOrb
            call hemv(cVecTemp, SSqrCplx, HSqrCplx(:,ii,iK2,iSpin2))
            write(fdEigvec, "(A,I4,4X,A,I4)") "K-point: ", ik, &
                &"Eigenvector: ", ii
            jj = 0
            do iAtom = 1, nAtom
              iSp1 = species(iAtom)
              do iSh1 = 1, orb%nShell(iSp1)
                ang = orb%angShell(iSh1,iSp1)
                if (iSh1 == 1) then
                  write(tmpStr, "(I5,1X,A2,2X,A1)") iAtom, speciesName(iSp1), &
                      &orbitalNames(ang+1)
                else
                  write(tmpStr, "(10X,A1)") orbitalNames(ang+1)
                end if
                do iOrb = 1, 2*ang+1
                  jj = jj + 1
                  write(fdEigvec,&
                      &"(A,I1,T15,'(',F12.6,',',F12.6,')', &
                      & '(',F12.6,',',F12.6,')',3X,4F12.6)") &
                      &trim(tmpStr), iOrb, &
                      &real(HSqrCplx(jj, ii, iK2, iSpin2)), &
                      &aimag(HSqrCplx(jj, ii, iK2, iSpin2)), &
                      &real(HSqrCplx(jj+nOrb, ii, iK2, iSpin2)), &
                      &aimag(HSqrCplx(jj+nOrb, ii, iK2, iSpin2)), &
                      &real( conjg(HSqrCplx(jj, ii, iK2, iSpin2)) &
                      & * cVecTemp(jj) + &
                      & conjg(HSqrCplx(jj+nOrb, ii, iK2, iSpin2)) &
                      & * cVecTemp(jj+nOrb) ), &
                      &real( conjg(HSqrCplx(jj+nOrb, ii, iK2, iSpin2)) &
                      & * cVecTemp(jj) + &
                      & conjg(HSqrCplx(jj, ii, iK2, iSpin2)) &
                      & * cVecTemp(jj+nOrb) ), &
                      &aimag( conjg(HSqrCplx(jj, ii, iK2, iSpin2)) &
                      & * cVecTemp(jj+nOrb) - &
                      & conjg(HSqrCplx(jj+nOrb, ii, iK2, iSpin2)) &
                      & * cVecTemp(jj) ), &
                      & real( conjg(HSqrCplx(jj, ii, iK2, iSpin2)) &
                      & * cVecTemp(jj) - &
                      & conjg(HSqrCplx(jj+nOrb, ii, iK2, iSpin2)) &
                      & * cVecTemp(jj+nOrb) )
                end do
              end do
              write (fdEigvec,*)
            end do
          end do
        end do

        ! normal spin block structure
      else

        allocate(cVecTemp(size(HSqrCplx, dim=1)))
        do iSpin = 1, nSpin
          do iK = 1, nK
            call getH(iSpin, iK, HSqrCplx, iSpin2, iK2, storeEigvecs)
            call unpackHS(SSqrCplx, over, kPoint(:,iK), neighlist%iNeighbor, &
                &nNeighbor, iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
            do ii = 1, orb%nOrb
              call hemv(cVecTemp, SSqrCplx, HSqrCplx(:,ii,iK2,iSpin2))
              write(fdEigvec, "(A,I4,4X,A,I4,4X,'(',A,')'/)") "K-point: ", ik, &
                  &"Eigenvector: ", ii, trim(spinName(iSpin))
              jj = 0
              do iAtom = 1, nAtom
                iSp1 = species(iAtom)
                do iSh1 = 1, orb%nShell(iSp1)
                  ang = orb%angShell(iSh1,iSp1)
                  if (iSh1 == 1) then
                    write(tmpStr, "(I5,1X,A2,2X,A1)") iAtom, speciesName(iSp1), &
                        &orbitalNames(ang+1)
                  else
                    write(tmpStr, "(10X,A1)") orbitalNames(ang+1)
                  end if
                  do iOrb = 1, 2*ang+1
                    jj = jj + 1
                    write(fdEigvec,&
                        &"(A,I1,T15,'(',F12.6,',',F12.6,')',3X,F12.6)") &
                        &trim(tmpStr), iOrb, &
                        &real(HSqrCplx(jj, ii, iK2, iSpin2)), &
                        &aimag(HSqrCplx(jj, ii, iK2, iSpin2)), &
                        & real( conjg(HSqrCplx(jj, ii, iK2, iSpin2)) &
                        & * cVecTemp(jj))
                  end do
                end do
                write (fdEigvec,*)
              end do
            end do
          end do
        end do
      end if

      close(fdEigvec)

    end if

  end subroutine writeCplxEigvecs


  !> Write the projected eigenstates into text files
  subroutine writeProjRealEigvecs(filenames, fdProjEig, ei, nSpin, neighlist, &
      &nNeighbor, iAtomStart, iPair, img2CentCell, orb, &
      &over, HSqrReal, SSqrReal, iOrbRegion, storeEigvecs)

    !> List with filenames for each region
    type(listCharLc), intent(inout) :: filenames

    !> File unit IDs for each of the regions
    integer, intent(in) :: fdProjEig(:)

    !> eigenvalues
    real(dp), intent(in) :: ei(:,:,:)

    !> Nr. of spin channels
    integer, intent(in) :: nSpin

    !> Neighbor list
    type(TNeighborList), intent(in) :: neighlist

    !> Nr. of neighbors for SK-interaction
    integer, intent(in) :: nNeighbor(:)

    !> Positions of atoms int the dense matrices
    integer, intent(in) :: iAtomStart(:)

    !> Positions of interactions in the sparse matrices
    integer, intent(in) :: iPair(:,:)

    !> Mapping of atoms into the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Square Hamiltonian (or work array)
    real(dp), intent(inout) :: HSqrReal(:,:,:)

    !> Work array
    real(dp), intent(inout) :: SSqrReal(:,:)

    !> orbital number in each region
    type(listIntR1), intent(inout) :: iOrbRegion

    !> If present, Hamiltonian(s) are fetched from this storage into HSqrReal, instead of using
    !> whatever is already there
    type(OFifoRealR2), intent(inout), optional :: storeEigvecs(:)

    integer, allocatable :: iOrbs(:)
    integer :: iSpin, iSpin2, iLev, ii, nReg, dummy
    integer :: valshape(1)
    real(dp) :: qState
    real(dp), allocatable :: rVecTemp(:)
    character(lc) :: tmpStr

    nReg = len(iOrbRegion)
    @:ASSERT(len(filenames) == nReg)

    @:ASSERT(size(fdProjEig) == nReg)
    @:ASSERT(all(fdProjEig > 0))

    do ii = 1, nReg
      call get(filenames, tmpStr, ii)
      open(fdProjEig(ii), file=tmpStr, action="write", status="replace")
    end do

    allocate(rVecTemp(size(HSqrReal, dim=1)))
    call unpackHS(SSqrReal, over, neighlist%iNeighbor, nNeighbor, &
        &iAtomStart, iPair, img2CentCell)
    do iSpin = 1, nSpin
      do ii = 1, nReg
        if (nSpin <= 2) then
          write(fdProjEig(ii),*)' KPT',1,' SPIN ', iSpin
        else
          write(fdProjEig(ii),*)' KPT',1
        endif
      end do
      call getH(iSpin, HSqrReal, iSpin2, storeEigvecs)
      do iLev = 1, orb%nOrb
        call hemv(rVecTemp, SSqrReal, HSqrReal(:,iLev,iSpin2))
        rVecTemp = rVecTemp * HSqrReal(:,iLev,iSpin2)
        do ii = 1, nReg
          call elemShape(iOrbRegion, valshape, ii)
          allocate(iOrbs(valshape(1)))
          call intoArray(iOrbRegion, iOrbs, dummy, ii)
          qState = sum(rVecTemp(iOrbs))
          deallocate(iOrbs)
          write(fdProjEig(ii), "(f13.6,f10.6)")Hartree__eV*ei(iLev,1,iSpin),&
              & qState
        end do
      end do
      if (iSpin < nSpin) then
        do ii = 1, nReg
          write(fdProjEig(ii),*)
        end do
      end if
    end do

    do ii = 1, nReg
      close(fdProjEig(ii))
    end do

  end subroutine writeProjRealEigvecs


  !> Write the projected complex eigenstates into text files.
  subroutine writeProjCplxEigvecs(filenames, fdProjEig, ei, nSpin, neighlist, &
      & nNeighbor, cellVec, iCellVec, iAtomStart, iPair, img2CentCell, orb, &
      & over, kpoint, kweight, HSqrCplx, SSqrCplx, iOrbRegion, storeEigvecs)

    !> list of region names
    type(ListCharLc), intent(inout) :: filenames

    !> Fileid (file not yet opened) to use.
    integer, intent(in) :: fdProjEig(:)

    !> eigenvalues
    real(dp), intent(in) :: ei(:,:,:)

    !> Nr. of spin channels.
    integer, intent(in) :: nSpin

    !> Neighbor list.
    type(TNeighborList), intent(in) :: neighlist

    !> Nr. of neighbors for SK-interaction.
    integer, intent(in) :: nNeighbor(:)

    !> Cell vectors of shifted cells.
    real(dp), intent(in) :: cellVec(:,:)

    !> Cell vector index of every atom.
    integer, intent(in) :: iCellVec(:)

    !> Positions of first basis function of each atom in the dense matrices.
    integer, intent(in) :: iAtomStart(:)

    !> Positions of interactions in the sparse matrices.
    integer, intent(in) :: iPair(:,:)

    !> Mapping of atoms into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Sparse overlap matrix.
    real(dp), intent(in) :: over(:)

    !> KPoints.
    real(dp), intent(in) :: kpoint(:,:)

    !> KPoints weights
    real(dp), intent(in) :: kweight(:)

    !> Square Hamiltonian (or work array)
    complex(dp), intent(inout) :: HSqrCplx(:,:,:,:)

    !> Work array.
    complex(dp), intent(inout) :: SSqrCplx(:,:)

    !> orbital number in each region
    type(listIntR1), intent(inout) :: iOrbRegion

    !> If present, Hamiltonian(s) are fetched from this storage into HSqrReal, instead of using
    !> whatever is already there.
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecs(:)

    integer, allocatable :: iOrbs(:)
    integer :: iSpin, iSpin2, iK, iK2, nK, iLev, ii, nReg, dummy, nOrb
    integer :: valshape(1)
    real(dp) :: qState
    complex(dp), allocatable :: cVecTemp(:), work(:,:)
    character(lc) :: tmpStr

    nK = size(kPoint, dim=2)
    nReg = len(iOrbRegion)
    @:ASSERT(len(filenames) == nReg)
    @:ASSERT(size(kweight) == nK)

    @:ASSERT(size(fdProjEig) == nReg)
    @:ASSERT(all(fdProjEig > 0))

    do ii = 1, nReg
      call get(filenames, tmpStr, ii)
      open(fdProjEig(ii), file=tmpStr, action="write", status="replace")
    end do

    allocate(cVecTemp(size(HSqrCplx, dim=1)))

    if (nSpin <= 2) then

      do iSpin = 1, nSpin
        do iK = 1, nK

          do ii = 1, nReg
            write(fdProjEig(ii),*)'KPT ',iK,' SPIN ', iSpin, &
                &' KWEIGHT ', kweight(iK)
          end do

          call getH(iSpin, iK, HSqrCplx, iSpin2, iK2, storeEigvecs)
          call unpackHS(SSqrCplx, over, kPoint(:,iK), neighlist%iNeighbor, &
              & nNeighbor, iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
          do iLev = 1, orb%nOrb
            call hemv(cVecTemp, SSqrCplx, HSqrCplx(:,iLev,iK2,iSpin2))
            cVecTemp = conjg(HSqrCplx(:,iLev,iK2,iSpin2)) * cVecTemp
            do ii = 1, nReg
              call elemShape(iOrbRegion, valshape, ii)
              allocate(iOrbs(valshape(1)))
              call intoArray(iOrbRegion, iOrbs, dummy, ii)
              qState = real(sum(cVecTemp(iOrbs)), dp)
              write(fdProjEig(ii), "(f13.6,f10.6)") &
                  & Hartree__eV * ei(iLev,iK,iSpin), qState
              deallocate(iOrbs)
            end do
          end do
          if (iK < nK .or. iSpin < nSpin) then
            do ii = 1, nReg
              write(fdProjEig(ii),*)
            end do
          end if
        end do
      end do

    else

      iSpin = 1
      nOrb = orb%nOrb
      allocate(work(nOrb,nOrb))

      do iK = 1, nK

        do ii = 1, nReg
          write(fdProjEig(ii),*)'KPT ',iK, ' KWEIGHT ', kweight(iK)
        end do
        call getH(iSpin, iK, HSqrCplx, iSpin2, iK2, storeEigvecs)
        SSqrCplx = 0.0_dp
        work = 0.0_dp
        call unpackHS(work, over, kPoint(:,iK), neighlist%iNeighbor, &
            & nNeighbor, iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
        SSqrCplx(:nOrb,:nOrb) = work
        SSqrCplx(nOrb+1:,nOrb+1:) = work
        do iLev = 1, 2*nOrb
          call hemv(cVecTemp, SSqrCplx, HSqrCplx(:,iLev,iK2,iSpin2))
          do ii = 1, nReg
            call elemShape(iOrbRegion, valshape, ii)
            allocate(iOrbs(valshape(1)))
            call intoArray(iOrbRegion, iOrbs, dummy, ii)
            write(fdProjEig(ii), "(f13.6,4f10.6)") &
                & Hartree__eV * ei(iLev,iK,iSpin), &
                & real(sum( conjg(HSqrCplx(iOrbs, iLev, iK2, iSpin2)) &
                & * cVecTemp(iOrbs) + &
                & conjg(HSqrCplx(iOrbs+nOrb, iLev, iK2, iSpin2)) &
                & * cVecTemp(iOrbs+nOrb) )), &
                & real(sum( conjg(HSqrCplx(iOrbs+nOrb, iLev, iK2, iSpin2)) &
                & * cVecTemp(iOrbs) + &
                & conjg(HSqrCplx(iOrbs, iLev, iK2, iSpin2)) &
                & * cVecTemp(iOrbs+nOrb) )), &
                & aimag(sum( conjg(HSqrCplx(iOrbs, iLev, iK2, iSpin2)) &
                & * cVecTemp(iOrbs+nOrb) - &
                & conjg(HSqrCplx(iOrbs+nOrb, iLev, iK2, iSpin2)) &
                & * cVecTemp(iOrbs) )), &
                & real(sum( conjg(HSqrCplx(iOrbs, iLev, iK2, iSpin2)) &
                & * cVecTemp(iOrbs) - &
                & conjg(HSqrCplx(iOrbs+nOrb, iLev, iK2, iSpin2)) &
                & * cVecTemp(iOrbs+nOrb) ))
            deallocate(iOrbs)
          end do
        end do
        if (iK < nK) then
          do ii = 1, nReg
            write(fdProjEig(ii),*)
          end do
        end if

      end do

    end if

    do ii = 1, nReg
      close(fdProjEig(ii))
    end do

  end subroutine writeProjCplxEigvecs


  !> Routines to get eigenvectors out of storage/memory
  subroutine getHreal(iSpin, HSqrReal, iSpin2, storeEigvecs)

    !> required spin index
    integer, intent(in) :: iSpin

    !> Square Hamiltonian
    real(dp), intent(inout) :: HSqrReal(:,:,:)

    !> spin index in returned HSqrReal
    integer, intent(out) :: iSpin2

    !> If present, Hamiltonian(s) are fetched from this storage into HSqrReal, instead of using
    !> whatever is already there.
    type(OFifoRealR2), intent(inout), optional :: storeEigvecs(:)

    if (present(storeEigvecs)) then
      iSpin2 = 1
      call get(storeEigvecs(iSpin), HSqrReal(:,:,iSpin2))
    else
      iSpin2 = iSpin
    end if

  end subroutine getHreal


  !> Routines to get eigenvectors out of storage/memory
  subroutine getHcmplx(iSpin, iK, HSqrCplx, iSpin2, iK2, storeEigvecs)

    !> required spin index
    integer, intent(in) :: iSpin

    !> required kpoint index
    integer, intent(in) :: iK

    !> Square Hamiltonian
    complex(dp), intent(inout) :: HSqrCplx(:,:,:,:)

    !> spin index in returned HSqrCplx
    integer, intent(out) :: iSpin2

    !> kpoint index in returned HSqrCplx
    integer, intent(out) :: iK2

    !> If present, Hamiltonian(s) are fetched from this storage into HSqrCplx, instead of using
    !> whatever is already there.
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecs(:)

    if (present(storeEigvecs)) then
      iSpin2 = 1
      iK2 = 1
      call get(storeEigvecs(iSpin), HSqrCplx(:,:,iK2,iSpin2))
    else
      iSpin2 = iSpin
      iK2 = iK
    end if

  end subroutine getHcmplx

  !> Open an output file and return its unit number
  subroutine initOutputFile(fileName, fd)

    !> File name
    character(*), intent(in) :: fileName

    !> Associated file ID
    integer, intent(out) :: fd

    fd = getFileId()
    open(fd, file=fileName, action="write", status="replace")
    close(fd)

  end subroutine initOutputFile


  !> Write tagged output of data from the code at the end of the DFTB+ run, data being then used for
  !> regression testing
  subroutine writeAutotestTag(fd, fileName, tPeriodic, cellVol, tMulliken, qOutput, derivs,&
      & chrgForces, excitedDerivs, tStress, totalStress, pDynMatrix, freeEnergy, pressure,&
      & gibbsFree, endCoords, tLocalise, localisation)

    !> File ID to write to
    integer, intent(in) :: fd

    !> Name of output file
    character(*), intent(in) :: fileName

    !> Is the geometry periodic
    logical, intent(in) :: tPeriodic

    !> Unit cell volume if periodic (unreferenced otherwise)
    real(dp), intent(in) :: cellVol

    !> Are Mulliken charges to be output
    logical, intent(in) :: tMulliken

    !> Output Mulliken charges
    real(dp), intent(in) :: qOutput(:,:,:)

    !> Atomic derivatives (allocation status used as a flag)
    real(dp), allocatable, intent(in) :: derivs(:,:)

    !> Forces on external charges (allocation status used as a flag)
    real(dp), allocatable, intent(in) :: chrgForces(:,:)

    !> Excited state forces on atoms (allocation status used as a flag)
    real(dp), allocatable, intent(in) :: excitedDerivs(:,:)

    !> Should stresses be printed (assumes periodic)
    logical, intent(in) :: tStress

    !> Stress tensor
    real(dp), intent(in) :: totalStress(:,:)

    !> Hessian (dynamical) matrix
    real(dp), pointer, intent(in) ::  pDynMatrix(:,:)

    !> Mermin electronic free energy
    real(dp), intent(in) :: freeEnergy

    !> External pressure
    real(dp), intent(in) :: pressure

    !> Gibbs free energy (includes pV term)
    real(dp), intent(in) :: gibbsFree

    !> Final atomic coordinates
    real(dp), intent(in) :: endCoords(:,:)

    !> Has localisation of single particle states been applied
    logical, intent(in) :: tLocalise

    !> Localisation measure, if relevant
    real(dp), intent(in) :: localisation

    real(dp), allocatable :: qOutputUpDown(:,:,:)

    @:ASSERT(tPeriodic .eqv. tSress)

    open(fd, file=fileName, action="write", status="old", position="append")
    if (tPeriodic) then
      call writeTagged(fd, tag_volume, cellVol)
    end if
    if (tMulliken) then
      qOutputUpDown = qOutput
      call qm2ud(qOutputUpDown)
      call writeTagged(fd, tag_qOutput, qOutputUpDown(:,:,1))
    end if
    if (allocated(derivs)) then
      call writeTagged(fd, tag_forceTot, -derivs)
    end if
    if (allocated(chrgForces)) then
      call writeTagged(fd, tag_chrgForces, -chrgForces)
    end if
    if (allocated(excitedDerivs)) then
      call writeTagged(fd, tag_excForce, -excitedDerivs)
    end if
    if (tStress) then
      call writeTagged(fd, tag_stressTot, totalStress)
    end if
    if (associated(pDynMatrix)) then
      call writeTagged(fd, tag_HessianNum, pDynMatrix)
    end if
    call writeTagged(fd, tag_freeEgy, freeEnergy)
    if (pressure /= 0.0_dp) then
      call writeTagged(fd, tag_Gibbsfree, gibbsFree)
    end if
    call writeTagged(fd, tag_endCoord, endCoords)
    if (tLocalise) then
      call writeTagged(fd, tag_pmlocalise, localisation)
    end if
    close(fd)

  end subroutine writeAutotestTag


  !> Writes out machine readable data
  subroutine writeResultsTag(fd, fileName, derivs, chrgForces, tStress, totalStress,&
      & pDynMatrix, tPeriodic, cellVol)

    !> File ID to write to
    integer, intent(in) :: fd

    !> Name of output file
    character(*), intent(in) :: fileName

    !> Atomic derivatives (allocation status used as a flag)
    real(dp), allocatable, intent(in) :: derivs(:,:)

    !> Forces on external charges
    real(dp), allocatable, intent(in) :: chrgForces(:,:)

    !> Should stresses be printed (assumes periodic)
    logical, intent(in) :: tStress

    !> Stress tensor
    real(dp), intent(in) :: totalStress(:,:)

    !> Hessian (dynamical) matrix
    real(dp), pointer, intent(in) :: pDynMatrix(:,:)

    !> Is the geometry periodic
    logical, intent(in) :: tPeriodic

    !> Unit cell volume if periodic (unreferenced otherwise)
    real(dp), intent(in) :: cellVol

    @:ASSERT(tPeriodic .eqv. tSress)

    open(fd, file=fileName, action="write", status="replace")
    if (allocated(derivs)) then
      call writeTagged(fd, tag_forceTot, -derivs)
    end if
    if (allocated(chrgForces)) then
      call writeTagged(fd, tag_chrgForces, -chrgForces)
    end if
    if (tStress) then
      call writeTagged(fd, tag_stressTot, totalStress)
    end if
    if (associated(pDynMatrix)) then
      call writeTagged(fd, tag_HessianNum, pDynMatrix)
    end if
    if (tPeriodic) then
      call writeTagged(fd, tag_volume, cellVol)
    end if
    close(fd)

  end subroutine writeResultsTag


  !> Write XML format of derived results
  subroutine writeDetailedXml(runId, speciesName, species0, coord0Out, tPeriodic, latVec, tRealHS,&
      & nKPoint, nSpin, nStates, nOrb, kPoint, kWeight, filling, occNatural)

    !> Identifier for the run
    integer, intent(in) :: runId

    !> Labels for the atomic species
    character(*), intent(in) :: speciesName(:)

    !> Species numbers for central cell atoms
    integer, intent(in) :: species0(:)

    !> coordinates of atoms
    real(dp), intent(in) :: coord0Out(:,:)

    !> Periodic boundary conditions
    logical, intent(in) :: tPeriodic

    !> Lattice vectors if periodic
    real(dp), intent(in) :: latVec(:,:)

    !> Real Hamiltonian
    logical, intent(in) :: tRealHS

    !> Number of k-points present
    integer, intent(in) :: nKPoint

    !> Number of spin channels present
    integer, intent(in) :: nSpin

    !> Number of eigen states in the system / dimension of the Hamiltonian
    integer, intent(in) :: nStates

    !> Number of atomic orbitals (may not match nStates if non-collinear)
    integer, intent(in) :: nOrb

    !> k-points in the system
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights of the k-points
    real(dp), intent(in) :: kWeight(:)

    !> Filling of the eigenstates
    real(dp), intent(in) :: filling(:,:,:)

    !> Occupation numbers for natural orbitals
    real(dp), allocatable, target, intent(in) :: occNatural(:)

    type(xmlf_t) :: xf
    real(dp), allocatable :: bufferRealR2(:,:)
    integer :: ii, jj
    real(dp), pointer :: pOccNatural(:,:)

    call xml_OpenFile("detailed.xml", xf, indent=.true.)
    call xml_ADDXMLDeclaration(xf)
    call xml_NewElement(xf, "detailedout")
    call writeChildValue(xf, "identity", runId)
    call xml_NewElement(xf, "geometry")
    call writeChildValue(xf, "typenames", speciesName)
    call writeChildValue(xf, "typesandcoordinates", reshape(species0, [ 1, size(species0) ]),&
        & coord0Out)
    call writeChildValue(xf, "periodic", tPeriodic)
    if (tPeriodic) then
      call writeChildValue(xf, "latticevectors", latVec)
    end if
    call xml_EndElement(xf, "geometry")
    call writeChildValue(xf, "real", tRealHS)
    call writeChildValue(xf, "nrofkpoints", nKPoint)
    call writeChildValue(xf, "nrofspins", nSpin)
    call writeChildValue(xf, "nrofstates", nStates)
    call writeChildValue(xf, "nroforbitals", nOrb)
    allocate(bufferRealR2(4, nKPoint))
    bufferRealR2(1:3, :) = kPoint
    bufferRealR2(4,:) = kWeight
    call writeChildValue(xf, "kpointsandweights", bufferRealR2)
    call xml_NewElement(xf, "occupations")
    do ii = 1, nSpin
      call xml_NewElement(xf, "spin" // i2c(ii))
      do jj = 1, nKpoint
        call writeChildValue(xf, "k" // i2c(jj), filling(:, jj, mod(ii,3)))
      end do
      call xml_EndElement(xf, "spin" // i2c(ii))
    end do
    call xml_EndElement(xf, "occupations")
    if (allocated(occNatural)) then
      call xml_NewElement(xf, "excitedoccupations")
      call xml_NewElement(xf, "spin" // i2c(1))
      pOccNatural(1:size(occNatural), 1:1) => occNatural
      call writeChildValue(xf, "k" // i2c(1), pOccNatural)
      call xml_EndElement(xf, "spin" // i2c(1))
      call xml_EndElement(xf, "excitedoccupations")
    end if

    call xml_EndElement(xf, "detailedout")
    call xml_Close(xf)

  end subroutine writeDetailedXml


  !> Write the band structure data out
  subroutine writeBandOut(fd, fileName, eigen, filling, kWeight)

    !> File  ID
    integer, intent(in) :: fd

    !> Name of file to write to
    character(*), intent(in) :: fileName

    !> Eigenvalues for states, k-points and spin indices
    real(dp), intent(in) :: eigen(:,:,:)

    !> Fillings of the states
    real(dp), intent(in) :: filling(:,:,:)

    !> Weights of the k-points
    real(dp), intent(in) :: kWeight(:)

    integer :: iSpin, iK, iEgy

    open(unit=fd, file=fileName, action="write", status="replace")
    do iSpin = 1, size(eigen, dim=3)
      do iK = 1, size(eigen, dim=2)
        write(fd, *) 'KPT ', iK, ' SPIN ', iSpin, ' KWEIGHT ', kWeight(iK)
        do iEgy = 1, size(eigen, dim=1)
          write(fd, "(2f12.5)") Hartree__eV * eigen(iEgy, iK, iSpin), filling(iEgy, iK, iSpin)
        end do
        write(fd,*)
      end do
    end do
    close(fd)

  end subroutine writeBandOut


  !> Write the second derivative matrix
  subroutine writeHessianOut(fd, fileName, pDynMatrix)

    !> File ID
    integer, intent(in) :: fd

    !> File name
    character(*), intent(in) :: fileName

    !> Dynamical (Hessian) matrix
    real(dp), intent(in) :: pDynMatrix(:,:)

    integer :: ii

    open(unit=fd, file=fileName, action="write", status="replace")
    do ii = 1, size(pDynMatrix, dim=2)
      write(fd, formatHessian) pDynMatrix(:, ii)
    end do
    close(fd)
    write(stdOut, "(2A)") 'Hessian matrix written to ', fileName

  end subroutine writeHessianOut

  !> First group of data to go to detailed.out
  subroutine writeDetailedOut1(fd, fileName, tAppendDetailedOut, iDistribFn, nGeoSteps, iGeoStep,&
      & tMD, tDerivs, tCoordOpt, tLatOpt, iLatGeoStep, iSccIter, energy, diffElec, sccErrorQ,&
      & indMovedAtom, coord0Out, q0, qInput, qOutput, eigen, filling, orb, species,&
      & tDFTBU, tImHam, tPrintMulliken, orbitalL, qBlockOut, Ef, Eband, TS, E0, pressure, cellVol,&
      & tAtomicEnergy, tDispersion, tEField, tPeriodic, nSpin, tSpinOrbit, tScc)

    !> File  ID
    integer, intent(in) :: fd

    !> Name of file to write to
    character(*), intent(in) :: fileName

    !> Append to the end of the file or overwrite
    logical, intent(in) :: tAppendDetailedOut

    !> Electron distribution choice
    integer, intent(in) :: iDistribFn

    !> Total number of geometry steps
    integer, intent(in) :: nGeoSteps

    !> Current geometry step
    integer, intent(in) :: iGeoStep

    !> Is this an molecular dynamics run
    logical, intent(in) :: tMD

    !> Is this a finite difference derivative calculation
    logical, intent(in) :: tDerivs

    !> Are atomic coordinates being optimised?
    logical, intent(in) :: tCoordOpt

    !> Is the lattice being optimised?
    logical, intent(in) :: tLatOpt

    !> Which step of lattice optimisation is occuring
    integer, intent(in) :: iLatGeoStep

    !> Whic scc step is occuring
    integer, intent(in) :: iSccIter

    !> Energy terms in the system
    type(TEnergies), intent(in) :: energy

    !> Change in energy from previous SCC iteration
    real(dp), intent(in) :: diffElec

    !> Input/output charge error for SCC
    real(dp), intent(in) :: sccErrorQ

    !> Moving atoms
    integer, intent(in) :: indMovedAtom(:)

    !> Output atomic coordinates
    real(dp), intent(in) :: coord0Out(:,:)

    !> Reference atomic charges
    real(dp), intent(in) :: q0(:,:,:)

    !> Input atomic charges (if SCC)
    real(dp), intent(in) :: qInput(:,:,:)

    !> Output atomic charges (if SCC)
    real(dp), intent(in) :: qOutput(:,:,:)

    !> Eigenvalues/single particle states
    real(dp), intent(in) :: eigen(:,:,:)

    !> Occupation numbers
    real(dp), intent(in) :: filling(:,:,:)

    !> Type containing atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Chemical species of atoms
    integer, intent(in) :: species(:)

    !> Are orbital potentials being used
    logical, intent(in) :: tDFTBU

    !> Does the Hamiltonian have an imaginary component (spin-orbit, magnetic field, ...)
    logical, intent(in) :: tImHam

    !> Should Mulliken populations be printed
    logical, intent(in) :: tPrintMulliken

    !> Orbital angular momentum (if available)
    real(dp), allocatable, intent(in) :: orbitalL(:,:,:)

    !> Output block (dual) Mulliken charges
    real(dp), allocatable, intent(in) :: qBlockOut(:,:,:,:)

    !> Fermi level
    real(dp), intent(in) :: Ef(:)

    !> Band energy
    real(dp), intent(in) :: EBand(:)

    !> Electron entropy times temperature
    real(dp), intent(in) :: TS(:)

    !> Zero temperature extrapolated electron energy
    real(dp), intent(in) :: E0(:)

    !> External pressure
    real(dp), intent(in) :: pressure

    !> Unit cell volume
    real(dp), intent(in) :: cellVol

    !> Are atom resolved energies required
    logical, intent(in) :: tAtomicEnergy

    !> Are dispersion interactions included
    logical, intent(in) :: tDispersion

    !> Is there an external electric field
    logical, intent(in) :: tEfield

    !> Is the system periodic
    logical, intent(in) :: tPeriodic

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> Are spin orbit interactions present
    logical, intent(in) :: tSpinOrbit

    !> Is this a self consistent charge calculation
    logical, intent(in) :: tScc

    real(dp), allocatable :: qInputUpDown(:,:,:), qOutputUpDown(:,:,:), qBlockOutUpDown(:,:,:,:)
    real(dp) :: angularMomentum(3)
    integer :: ang
    integer :: nAtom, nLevel, nKPoint, nSpinHams, nMovedAtom
    integer :: iAt, iSpin, iEgy, iK, iSp, iSh, iOrb, kk
    logical :: tSpin

    character(*), parameter :: formatEigen = "(F14.8)"
    character(*), parameter :: formatFilling = "(F12.5)"

    nAtom = size(q0, dim=2)
    nLevel = size(eigen, dim=1)
    nKPoint = size(eigen, dim=2)
    nSpinHams = size(eigen, dim=3)
    nMovedAtom = size(indMovedAtom)
    tSpin = (nSpin == 2 .or. nSpin == 4)

    qInputUpDown = qInput
    call qm2ud(qInputUpDown)
    qOutputUpDown = qOutput
    call qm2ud(qOutputUpDown)
    if (allocated(qBlockOut)) then
      qBlockOutUpDown = qBlockOut
      call qm2ud(qBlockOutUpDown)
    end if

    if (iGeoStep == 0 .and. iSccIter == 1) then
      open(fd, file=fileName, status="replace", action="write")
    elseif (.not. tAppendDetailedOut) then
      close(fd)
      open(fd, file=fileName, status="replace", action="write")
    end if

    select case(iDistribFn)
    case(0)
      write(fd,*) 'Fermi distribution function'
    case(1)
      write(fd,*) 'Gaussian distribution function'
    case default
      write(fd,*) 'Methfessel-Paxton distribution function order', iDistribFn
    end select
    write(fd,*)

    if (nGeoSteps > 0) then
      if (tMD) then
        write(fd, "(A, I0)") "MD step: ", iGeoStep
      elseif (tDerivs) then
        write(fd, "(A, I0)") 'Difference derivative step: ', iGeoStep
      else
        if (tCoordOpt .and. tLatOpt) then
          write (fd, "(A, I0, A, I0)") "Geometry optimization step: ", iGeoStep,&
              & ", Lattice step: ", iLatGeoStep
        else
          write(fd, "(A, I0)") "Geometry optimization step: ", iGeoStep
        end if
      end if
    elseif (tScc) then
      ! Only written if scc is on, to be compatible with old output
      write(fd, "(A)") "Calculation with static geometry"
    end if
    write(fd, *)

    if (tSCC) then
      write(fd, "(/, A)") repeat("*", 80)
      write(fd, "(A5, A18, A18, A18)") "iSCC", " Total electronic ", "  Diff electronic ",&
          & "     SCC error    "
      write(fd, "(I5, E18.8, E18.8, E18.8, E18.8)") iSCCIter, energy%Eelec, diffElec, sccErrorQ
      write(fd, "(A)") repeat("*", 80)
      write(fd, *)
    end if

    if (nMovedAtom > 0 .and. .not. tDerivs) then
      write(fd, "(A)") "Coordinates of moved atoms (au):"
      do iAt = 1, nMovedAtom
        write(fd, formatGeoOut) indMovedAtom(iAt), coord0Out(:, indMovedAtom(iAt))
      end do
      write(fd, *)
    end if

    ! Write out atomic charges
    if (tPrintMulliken) then
      write(fd, "(A, F14.8)") " Net charge: ", sum(q0(:, :, 1) - qOutput(:, :, 1))
      write(fd, "(/,A)") " Net atomic charges (e)"
      write(fd, "(A5, 1X, A16)")" Atom", " Net charge"
      do iAt = 1, nAtom
        write(fd, "(I5, 1X, F16.8)") iAt, sum(q0(:, iAt, 1) - qOutput(:, iAt, 1))
      end do
      write(fd, *)
    end if

    lpSpinPrint: do iSpin = 1, size(eigen, dim=3)
      if (nSpin == 2) then
        write(fd, "(2A)") 'COMPONENT = ', trim(spinName(iSpin))
      else
        write(fd, "(2A)") 'COMPONENT = ', trim(quaternionName(iSpin))
      end if
      write(fd, "(/, A)") 'Eigenvalues /H'
      do iEgy = 1, size(eigen, dim=1)
        write(fd, formatEigen) (eigen(iEgy, iK, iSpin), iK = 1, nKPoint)
      end do
      write(fd, "(/, A)") 'Eigenvalues /eV'
      do iEgy = 1, size(eigen, dim=1)
        write(fd, formatEigen) (Hartree__eV * eigen(iEgy, iK, iSpin), iK = 1, nKPoint)
      end do
      write(fd, "(/, A)") 'Fillings'
      do iEgy = 1, nLevel
        write(fd, formatFilling) (filling(iEgy, iK, iSpin), iK = 1, nKPoint)
      end do
      write(fd, *)
    end do lpSpinPrint

    if (nSpin == 4) then
      if (tPrintMulliken) then
        do iSpin = 1, 4
          write(fd,"(3A, F16.8)") 'Nr. of electrons (', quaternionName(iSpin), '):',&
              & sum(qOutput(:,:, iSpin))
          write(fd, *)
          write(fd, "(/, 3A)") 'Atom populations (', quaternionName(iSpin), ')'
          write(fd, "(A5, 1X, A16)") " Atom", " Population"
          do iAt = 1, nAtom
            write(fd, "(1X, I5, 1X, F16.8)") iAt, sum(qOutput(:, iAt, iSpin))
          end do
          write(fd, "(/, 3A)") 'l-shell populations (', quaternionName(iSpin), ')'
          write(fd, "(A5, 1X, A3, 1X, A3, 1X, A16)") " Atom", "Sh.", "  l", " Population"
          do iAt = 1, nAtom
            iSp = species(iAt)
            do iSh = 1, orb%nShell(iSp)
              write(fd, "(I5, 1X, I3, 1X, I3, 1X, F16.8)") iAt, iSh, orb%angShell(iSh, iSp), &
                  & sum(qOutput(orb%posShell(iSh,iSp):orb%posShell(iSh+1, iSp) - 1, iAt, iSpin))
            end do
          end do
          write(fd,*)
          write(fd, "(/, 3A)") 'Orbital populations (', quaternionName(iSpin) ,')'
          write(fd, "(A5, 1X, A3, 1X, A3, 1X, A3, 1X, A16)") " Atom", "Sh.","  l","  m",&
              & " Population"
          do iAt = 1, nAtom
            iSp = species(iAt)
            do iSh = 1, orb%nShell(iSp)
              ang = orb%angShell(iSh, iSp)
              do kk = 0, 2 * ang
                write(fd, "(I5, 1X, I3, 1X, I3, 1X, I3, 1X, F16.8)") &
                    & iAt, iSh, ang, kk - ang, qOutput(orb%posShell(iSh, iSp) + kk, iAt, iSpin)
              end do
            end do
          end do
          write(fd, *)
        end do
      end if

      if (tDFTBU) then
        do iSpin = 1, 4
          write(fd, "(3A)") 'Block populations (', quaternionName(iSpin), ')'
          do iAt = 1, nAtom
            iSp = species(iAt)
            write(fd, "(A, 1X, I0)") 'Atom', iAt
            do iOrb = 1, orb%nOrbSpecies(iSp)
              write(fd, "(16F8.4)") qBlockOut(1:orb%nOrbSpecies(iSp), iOrb, iAt, iSpin)
            end do
            write(fd, *)
          end do
        end do
      end if

      if (tImHam .and. tPrintMulliken) then
        write(fd, "(/, A)") 'Electron angular momentum (mu_B/hbar)'
        write(fd, "(2X, A5, T10, A3, T14, A1, T20, A1, T35, A9)")&
            & "Atom", "Sh.", "l", "S", "Momentum"
        do iAt = 1, nAtom
          iSp = species(iAt)
          do iSh = 1, orb%nShell(iSp)
            write(fd, "(I5, 1X, I3, 1X, I3, 1X, F14.8, ' :', 3F14.8)") &
                & iAt, iSh, orb%angShell(iSh, iSp),&
                & 0.5_dp * sqrt(sum(sum(qOutput(orb%posShell(iSh, iSp)&
                & :orb%posShell(iSh + 1, iSp) - 1, iAt, 2:4), dim=1)**2)), &
                & -gfac * 0.25_dp * sum(qOutput(orb%posShell(iSh, iSp)&
                & :orb%posShell(iSh + 1, iSp) - 1, iAt, 2:4), dim=1)
          end do
        end do
        write(fd, "(/, A)") 'Orbital angular momentum (mu_B/hbar)'
        write(fd, "(2X, A5, T10, A3, T14, A1, T20, A1, T35, A9)")&
            & "Atom", "Sh.", "l", "L", "Momentum"
        do iAt = 1, nAtom
          iSp = species(iAt)
          do iSh = 1, orb%nShell(iSp)
            write(fd, "(I5, 1X, I3, 1X, I3, 1X, F14.8, ' :', 3F14.8)") &
                & iAt, iSh, orb%angShell(iSh, iSp), &
                & sqrt(sum(orbitalL(1:3, iSh, iAt)**2)), -orbitalL(1:3, iSh, iAt)
          end do
        end do

        write(fd, *)
        write(fd, "(A)") 'Total angular momentum (mu_B/hbar)'
        write(fd, "(2X, A5, T10, A3, T14, A1, T20, A1, T35, A9)")&
            & "Atom", "Sh.", "l", "J", "Momentum"
        angularMomentum(:) = 0.0_dp
        do iAt = 1, nAtom
          iSp = species(iAt)
          do iSh = 1, orb%nShell(iSp)
            write(fd, "(I5, 1X, I3, 1X, I3, 1X, F14.8, ' :', 3F14.8)") &
                & iAt, iSh, orb%angShell(iSh, iSp),&
                & sqrt(sum((orbitalL(1:3, iSh, iAt)&
                & + sum(0.5_dp * qOutput(orb%posShell(iSh, iSp)&
                & :orb%posShell(iSh + 1, iSp) - 1, iAt, 2:4), dim=1))**2)), &
                & -orbitalL(1:3, iSh, iAt) &
                & -gfac * 0.25_dp * sum(qOutput(orb%posShell(iSh, iSp)&
                & :orb%posShell(iSh + 1, iSp) - 1, iAt, 2:4), dim=1)
            angularMomentum(1:3) = angularMomentum(1:3) &
                & -orbitalL(1:3, iSh, iAt) &
                & -gfac * 0.25_dp * sum(qOutput(orb%posShell(iSh, iSp)&
                & :orb%posShell(iSh + 1, iSp) - 1, iAt, 2:4), dim=1)
          end do
        end do
        write(fd, *)
      end if
    else
      lpSpinPrint2: do iSpin = 1, nSpin
        if (tPrintMulliken) then
          write(fd, "(3A, F16.8)") 'Nr. of electrons (', trim(spinName(iSpin)), '):',&
              & sum(qOutputUpDown(:, :, iSpin))
          write(fd, "(3A)") 'Atom populations (', trim(spinName(iSpin)), ')'
          write(fd, "(A5, 1X, A16)") " Atom", " Population"
          do iAt = 1, nAtom
            write(fd, "(I5, 1X, F16.8)") iAt, sum(qOutputUpDown(:, iAt, iSpin))
          end do
          write(fd, *)
          write(fd, "(3A)") 'l-shell populations (', trim(spinName(iSpin)), ')'
          write(fd, "(A5, 1X, A3, 1X, A3, 1X, A16)")&
              & " Atom", "Sh.", "  l", " Population"
          do iAt = 1, nAtom
            iSp = species(iAt)
            do iSh = 1, orb%nShell(iSp)
              write(fd, "(I5, 1X, I3, 1X, I3, 1X, F16.8)") iAt, iSh, orb%angShell(iSh, iSp),&
                  & sum(qOutputUpDown(orb%posShell(iSh, iSp):orb%posShell(iSh + 1, iSp)-1, iAt,&
                  & iSpin))
            end do
          end do
          write(fd, *)
          write(fd, "(3A)") 'Orbital populations (', trim(spinName(iSpin)), ')'
          write(fd, "(A5, 1X, A3, 1X, A3, 1X, A3, 1X, A16)")&
              & " Atom", "Sh.", "  l", "  m", " Population"
          do iAt = 1, nAtom
            iSp = species(iAt)
            do iSh = 1, orb%nShell(iSp)
              ang = orb%angShell(iSh, iSp)
              do kk = 0, 2 * ang
                write(fd, "(I5, 1X, I3, 1X, I3, 1X, I3, 1X, F16.8)") &
                    &iAt, iSh, ang, kk - ang, qOutputUpDown(orb%posShell(iSh, iSp) + kk, iAt, iSpin)
              end do
            end do
          end do
          write(fd, *)
        end if
        if (tDFTBU) then
          write(fd, "(3A)") 'Block populations (', trim(spinName(iSpin)), ')'
          do iAt = 1, nAtom
            iSp = species(iAt)
            write(fd, "(A, 1X, I0)") 'Atom', iAt
            do iOrb = 1, orb%nOrbSpecies(iSp)
              write(fd, "(16F8.4)") qBlockOutUpDown(1:orb%nOrbSpecies(iSp), iOrb, iAt, iSpin)
            end do
          end do
          write(fd, *)
        end if
      end do lpSpinPrint2
    end if

    lpSpinPrint3: do iSpin = 1, nSpinHams
      if (nSpin == 2) then
        write(fd, "(A, 1X, A)") 'Spin ', trim(spinName(iSpin))
      end if
      write(fd, format2U) 'Fermi level', Ef(iSpin), "H", Hartree__eV * Ef(iSpin), 'eV'
      write(fd, format2U) 'Band energy', Eband(iSpin), "H", Hartree__eV * Eband(iSpin), 'eV'
      write(fd, format2U)'TS', TS(iSpin), "H", Hartree__eV * TS(iSpin), 'eV'
      write(fd, format2U) 'Band free energy (E-TS)', Eband(iSpin) - TS(iSpin), "H",&
          & Hartree__eV * (Eband(iSpin) - TS(iSpin)), 'eV'
      write(fd, format2U) 'Extrapolated E(0K)', E0(iSpin), "H", Hartree__eV * (E0(iSpin)), 'eV'
      if (tPrintMulliken) then
        if (nSpin == 2) then
          write(fd, "(3A, 2F16.8)") 'Input / Output electrons (', trim(spinName(iSpin)), '):',&
              & sum(qInputUpDown(:, :, iSpin)), sum(qOutputUpDown(:, :, iSpin))
        else
          write(fd, "(3A, 2F16.8)") 'Input / Output electrons (', quaternionName(iSpin), '):',&
              & sum(qInputUpDown(:, :, iSpin)), sum(qOutputUpDown(:, :, iSpin))
        end if
      end if
      write(fd, *)
    end do lpSpinPrint3

    write(fd, format2U) 'Energy H0', energy%EnonSCC, 'H', energy%EnonSCC * Hartree__eV, 'eV'

    if (tSCC) then
      write(fd, format2U) 'Energy SCC', energy%ESCC, 'H', energy%ESCC * Hartree__eV, 'eV'
      if (tSpin) then
        write(fd, format2U) 'Energy SPIN', energy%Espin, 'H', energy%Espin * Hartree__eV, 'eV'
      end if
      if (tDFTBU) then
        write(fd, format2U) 'Energy DFTB+U', energy%Edftbu, 'H',&
            & energy%Edftbu * Hartree__eV, 'eV'
      end if
    end if

    if (tSpinOrbit) then
      write(fd, format2U) 'Energy L.S', energy%ELS, 'H', energy%ELS * Hartree__eV, 'eV'
    end if

    if (tEfield) then
      write(fd, format2U) 'Energy ext. field', energy%Eext, 'H', energy%Eext * Hartree__eV, 'eV'
    end if

    write(fd, format2U) 'Total Electronic energy', energy%Eelec, 'H', &
        & energy%Eelec * Hartree__eV, 'eV'
    write(fd, format2U) 'Repulsive energy', energy%Erep, 'H', energy%Erep * Hartree__eV, 'eV'

    if (tDispersion) then
      write(fd, format2U) 'Dispersion energy', energy%eDisp, 'H',&
          & energy%eDisp * Hartree__eV, 'eV'
    end if

    write(fd, format2U) 'Total energy', energy%Etotal, 'H', energy%Etotal * Hartree__eV, 'eV'
    write(fd, format2U) 'Total Mermin free energy', energy%Etotal - sum(TS), 'H',&
        & (energy%Etotal - sum(TS)) * Hartree__eV, 'eV'
    if (tPeriodic .and. pressure /= 0.0_dp) then
      write(fd, format2U) 'Gibbs free energy', energy%Etotal - sum(TS) + cellVol * pressure,&
          & 'H', Hartree__eV * (energy%Etotal - sum(TS) + cellVol * pressure), 'eV'
    end if
    write(fd, *)

    if (tAtomicEnergy) then
      write(fd, "(A)") 'Atom resolved electronic energies '
      do iAt = 1, nAtom
        write(fd, "(I5, F16.8, A, F16.6, A)") iAt, energy%atomElec(iAt), ' H',&
            & Hartree__eV * energy%atomElec(iAt), ' eV'
      end do
      write(fd, *)

      write(fd, "(A)") 'Atom resolved repulsive energies '
      do iAt = 1, nAtom
        write(fd, "(I5, F16.8, A, F16.6, A)") iAt, energy%atomRep(iAt), ' H',&
            & Hartree__eV * energy%atomRep(iAt), ' eV'
      end do
      write(fd, *)
      write(fd, "(A)") 'Atom resolved total energies '
      do iAt = 1, nAtom
        write(fd, "(I5, F16.8, A, F16.6, A)") iAt, energy%atomTotal(iAt), ' H',&
            & Hartree__eV * energy%atomTotal(iAt), ' eV'
      end do
      write(fd, *)
    end if

  end subroutine writeDetailedOut1


  !> Second group of data for detailed.out
  subroutine writeDetailedOut2(fd, tScc, tConverged, tXlbomd, tLinResp, tGeoOpt, tMd, tPrintForces,&
      & tStress, tPeriodic, energy, totalStress, totalLatDeriv, derivs, chrgForces,&
      & indMovedAtom, cellVol, cellPressure, geoOutFile)

    !> File  ID
    integer, intent(in) :: fd

    !> Charge self consistent?
    logical, intent(in) :: tScc

    !> Has the SCC cycle converged?
    logical, intent(in) :: tConverged

    !> Is the extended Lagrangian in use for MD
    logical, intent(in) :: tXlbomd

    !> Is the Casida excited state in use?
    logical, intent(in) :: tLinResp

    !> Is the geometry being optimised
    logical, intent(in) :: tGeoOpt

    !> Is this a molcular dynamics run
    logical, intent(in) :: tMd

    !> Should forces be printed out?
    logical, intent(in) :: tPrintForces

    !> Is the stress tensor to be printed?
    logical, intent(in) :: tStress

    !> Is the geometry periodic
    logical, intent(in) :: tPeriodic

    !> Structure containing energy contributions
    type(TEnergies), intent(in) :: energy

    !> Stress tensor
    real(dp), intent(in) :: totalStress(:,:)

    !> Derivative with respect to lattice vectors
    real(dp), intent(in) :: totalLatDeriv(:,:)

    !> Energy derivative with respect to atomic coordinates
    real(dp), intent(in), allocatable :: derivs(:,:)

    !> Forces on external charges
    real(dp), intent(in), allocatable :: chrgForces(:,:)

    !> Index of moving atoms
    integer, intent(in) :: indMovedAtom(:)

    !> Unit cell volume
    real(dp), intent(in) :: cellVol

    !> Internal pressure in the unit cell
    real(dp), intent(in) :: cellPressure

    !> File for geometry output
    character(*), intent(in) :: geoOutFile

    integer :: iAt, ii

    if (tScc) then
      if (tConverged) then
        write(fd, "(A)") "SCC converged"
        write(fd, *)
      else
        if (.not. tXlbomd) then
          write(fd, "(A)") "SCC is NOT converged, maximal SCC iterations exceeded"
          write(fd, *)
        end if
      end if
    else
      write(fd, "(A)") "Non-SCC calculation"
      write(fd, *)
    end if

    ! only print excitation energy if 1) its been calculated and 2) its avaialable for a single
    ! state
    if (tLinResp .and. energy%Eexcited /= 0.0_dp) then
      write(fd, format2U) "Excitation Energy", energy%Eexcited, "H", &
          & Hartree__eV * energy%Eexcited, "eV"
      write(fd, *)
    end if

    if (tGeoOpt .or. tMd) then
      write(fd, "(3A)") "Full geometry written in ", trim(geoOutFile), ".{xyz|gen}"
      write(fd, *)
    end if

    if (tPrintForces) then
      write(fd, "(A)") 'Total Forces'
      do iAt = 1, size(derivs, dim=2)
        write(fd, "(3F20.12)") -derivs(:, iAt)
      end do
      write(fd, *)
      if (tStress .and. .not. tMd) then
        write(fd, "(A)") 'Total stress tensor'
        do ii = 1, 3
          write(fd, "(3F20.12)") totalStress(:, ii)
        end do
        write(fd, *)
        write(fd, "(A)") 'Total lattice derivs'
        do ii = 1, 3
          write(fd, "(3F20.12)") totalLatDeriv(:, ii)
        end do
        write(fd, *)
      end if

      write(fd, format1Ue) "Maximal derivative component", maxval(abs(derivs)), 'au'
      if (size(indMovedAtom) > 0) then
        write(fd, format1Ue) "Max force for moved atoms:",&
            & maxval(abs(derivs(:, indMovedAtom))), 'au'
      end if
      write(fd, *)

      if (allocated(chrgForces)) then
        write(fd, "(A)") "Forces on external charges"
        do ii = 1, size(chrgForces, dim=2)
          write(fd, "(3F20.12)") -chrgForces(:, ii)
        end do
        write(fd, *)
      end if

      if (tPeriodic .and. .not. tMd) then
        write(fd, format1Ue) 'Volume', cellVol, 'au^3'
        if (tStress) then
          write(fd, format2Ue)'Pressure', cellPressure, 'au', cellPressure * au__pascal, 'Pa'
        end if
        write(fd, *)
      end if
    end if

  end subroutine writeDetailedOut2


  !> Third group of data for detailed.out
  subroutine writeDetailedOut3(fd, tPrintForces, tSetFillingTemp, tPeriodic, tStress, totalStress,&
      & totalLatDeriv, energy, tempElec, pressure, cellPressure, tempIon)

    !> File ID
    integer, intent(in) :: fd

    !> Print forces on atoms
    logical, intent(in) :: tPrintForces

    !> If the electronic temperature is being set during the run
    logical, intent(in) :: tSetFillingTemp

    !> Is this a periodic geometry
    logical, intent(in) :: tPeriodic

    !> Should the stress tensor/lattice derivatives be printed?
    logical, intent(in) :: tStress

    !> Stress tensor
    real(dp), intent(in) :: totalStress(:,:)

    !> Energy derivatives with respect to lattice vectors
    real(dp), intent(in) :: totalLatDeriv(:,:)

    !> Data structure for energy components
    type(TEnergies), intent(in) :: energy

    !> electron temperature
    real(dp), intent(in) :: tempElec

    !> External pressure
    real(dp), intent(in) :: pressure

    !> Internal pressure in the unit cell
    real(dp), intent(in) :: cellPressure

    !> Atomic kinetic temperature
    real(dp), intent(in) :: tempIon

    integer :: ii

    if (tStress .and. tPrintForces) then
      write(fd, "(A)") 'Total stress tensor'
      do ii = 1, 3
        write(fd, "(3F20.12)") totalStress(:, ii)
      end do
      write(fd, *)
      write(fd, "(A)") 'Total lattice derivs'
      do ii = 1, 3
        write(fd, "(3F20.12)") totalLatDeriv(:, ii)
      end do
      write(fd, *)
    end if

    if (tSetFillingTemp) then
      write (fd, format2U) "Electronic Temperature", tempElec, 'au', tempElec * Hartree__eV,&
          & 'eV'
    end if
    write(fd, format1U) "MD Kinetic Energy", energy%EKin, "H"
    write(fd, format1U) "Total MD Energy", energy%EKin + energy%EMermin, "H"
    if (tPeriodic) then
      write(fd, format2Ue) 'Pressure', cellPressure, 'au', cellPressure * au__pascal, 'Pa'
      if (pressure /= 0.0_dp) then
        write(fd, format2U) 'Gibbs free energy including KE', energy%EGibbsKin, 'H',&
            & Hartree__eV * energy%EGibbsKin, 'eV'
      end if
    end if
    write(fd, format2U) "MD Temperature", tempIon, "H", tempIon / Boltzmann, "K"

  end subroutine writeDetailedOut3

  !> Fourth group of data for detailed.out
  subroutine writeDetailedOut4(fd, tMd, energy, tempIon)

    !> File ID
    integer, intent(in) :: fd

    !> Is this an MD calculation
    logical, intent(in) :: tMd

    !> Energy contributions
    type(TEnergies), intent(in) :: energy

    !> Atomic kinetic energy
    real(dp), intent(in) :: tempIon

    if (tMd) then
      write(fd, format1U) "MD Kinetic Energy", energy%Ekin, "H"
      write(fd, format2U) "Total MD Energy", energy%EMerminKin, "H",&
          & Hartree__eV * energy%EMerminKin, "eV"
      write(fd, format2U) "MD Temperature", tempIon, "H", tempIon / Boltzmann, "K"
      write(fd, *)
    end if

  end subroutine writeDetailedOut4

  !> Fifth group of data for detailed.out
  subroutine writeDetailedOut5(fd, tGeoOpt, tGeomEnd, tMd, tDerivs, tEField, absEField,&
      & dipoleMoment)

    !> File ID
    integer, intent(in) :: fd

    !> Is the geometry changing during the run
    logical, intent(in) :: tGeoOpt

    !> Did the geometry changes sucessfully complete
    logical, intent(in) :: tGeomEnd

    !> Is this a molecular dynamics run
    logical, intent(in) :: tMd

    !> Are finite difference derivatives being computed
    logical, intent(in) :: tDerivs

    !> Is there an external electric field
    logical, intent(in) :: tEField

    !> What is the external E field magnitude
    real(dp), intent(in) :: absEField

    !> What is the dipole moment (if available)
    real(dp), intent(in), optional :: dipoleMoment(:)

    if (tEfield) then
      write(fd, format1U1e) 'External E field', absEField, 'au', absEField * au__V_m, 'V/m'
    end if

    if (present(dipoleMoment)) then
      write(fd, "(A, 3F14.8, A)") 'Dipole moment:', dipoleMoment, ' au'
      write(fd, "(A, 3F14.8, A)") 'Dipole moment:', dipoleMoment * au__Debye, ' Debye'
      write(fd, *)
    end if

    if (tGeoOpt) then
      if (tGeomEnd) then
        write(fd, "(A)") "Geometry converged"
      else
        write(fd, "(A)") "!!! Geometry did NOT converge!"
      end if
    elseif (tMD) then
      if (tGeomEnd) then
        write(fd, "(A)") "Molecular dynamics completed"
      else
        write(fd, "(A)") "!!! Molecular dynamics terminated abnormally!"
      end if
    elseif (tDerivs) then
      if (tGeomEnd) then
        write(fd, "(A)") "Second derivatives completed"
      else
        write(fd, "(A)") "!!! Second derivatives terminated abnormally!"
      end if
    end if
    write(fd,*)
    close(fd)

  end subroutine writeDetailedOut5

  !> First group of output data during molecular dynamics
  subroutine writeMdOut1(fd, fileName, iGeoStep, pMdIntegrator)

    !> File ID
    integer, intent(in) :: fd

    !> File name
    character(*), intent(in) :: fileName

    !> Number of the current geometry step
    integer, intent(in) :: iGeoStep

    !> Molecular dynamics integrator
    type(OMdIntegrator), intent(in) :: pMdIntegrator

    if (iGeoStep == 0) then
      open(fd, file=fileName, status="replace", action="write")
    end if
    write(fd, "(A, 1X, I0)") "MD step:", iGeoStep
    call state(pMdIntegrator, fd)

  end subroutine writeMdOut1

  !> Second group of output data during molecular dynamics
  subroutine writeMdOut2(fd, tStress, tBarostat, tLinResp, tEField, tFixEf, tPrintMulliken,&
      & energy, latVec, cellVol, cellPressure, pressure, tempIon, absEField, qOutput, q0,&
      & dipoleMoment)

    !> File ID
    integer, intent(in) :: fd

    !> Is the stress tensor to be printed?
    logical, intent(in) :: tStress

    !> Is a barostat in use
    logical, intent(in) :: tBarostat

    !> Is linear response excitation being used
    logical, intent(in) :: tLinResp

    !> External electric field
    logical, intent(in) :: tEField

    !> Is the  Fermi level fixed
    logical, intent(in) :: tFixEf

    !> Should Mulliken charges be printed, hence total charge here
    logical, intent(in) :: tPrintMulliken

    !> energy contributions
    type(TEnergies), intent(in) :: energy

    !> Lattice vectors if periodic
    real(dp), intent(in) :: latVec(:,:)

    !> Unit cell volume
    real(dp), intent(in) :: cellVol

    !> Internal cell pressure
    real(dp), intent(in) :: cellPressure

    !> External applied pressure
    real(dp), intent(in) :: pressure

    !> Atomic kinetic energy
    real(dp), intent(in) :: tempIon

    !> magnitude of any applied electric field
    real(dp), intent(in) :: absEField

    !> Output atomic charges (if SCC)
    real(dp), intent(in) :: qOutput(:,:,:)

    !> Reference atomic charges
    real(dp), intent(in) :: q0(:,:,:)

    !> dipole moment if available
    real(dp), intent(in), optional :: dipoleMoment(:)

    integer :: ii

    if (tStress) then
      if (tBarostat) then
        write(fd, "(A)") 'Lattice vectors (A)'
        do ii = 1, 3
          write(fd, "(3E24.8)") latVec(:,ii) * Bohr__AA
        end do
        write(fd, format2Ue) 'Volume', cellVol, 'au^3', (Bohr__AA**3) * cellVol, 'A^3'
      end if
      write(fd, format2Ue) 'Pressure', cellPressure, 'au', cellPressure * au__pascal, 'Pa'
      if (pressure /= 0.0_dp) then
        write(fd, format2U) 'Gibbs free energy', energy%EGibbs, 'H',&
            & Hartree__eV * energy%EGibbs,'eV'
        write(fd, format2U) 'Gibbs free energy including KE', energy%EGibbsKin, 'H',&
            & Hartree__eV * energy%EGibbsKin, 'eV'
      end if
    end if
    if (tLinResp .and. energy%Eexcited /= 0.0_dp) then
      write (fd, format2U) "Excitation Energy", energy%Eexcited, "H",&
          & Hartree__eV * energy%Eexcited, "eV"
    end if
    write(fd, format2U) 'Potential Energy', energy%EMermin,'H', energy%EMermin * Hartree__eV, 'eV'
    write(fd, format2U) 'MD Kinetic Energy', energy%Ekin, 'H', energy%Ekin * Hartree__eV, 'eV'
    write(fd, format2U) 'Total MD Energy', energy%EMerminKin, 'H',&
        & energy%EMerminKin * Hartree__eV, 'eV'
    write(fd, format2U) 'MD Temperature', tempIon, 'au', tempIon / Boltzmann, 'K'
    if (tEfield) then
      write(fd, format1U1e) 'External E field', absEField, 'au', absEField * au__V_m, 'V/m'
    end if
    if (tFixEf .and. tPrintMulliken) then
      write(fd, "(A, F14.8)") 'Net charge: ', sum(q0(:, :, 1) - qOutput(:, :, 1))
    end if
    if (present(dipoleMoment)) then
      write(fd, "(A, 3F14.8, A)") 'Dipole moment:', dipoleMoment,  'au'
      write(fd, "(A, 3F14.8, A)") 'Dipole moment:', dipoleMoment * au__Debye,  'Debye'
    end if

  end subroutine writeMdOut2

  !> Third and final group of output data during molecular dynamics
  subroutine writeMdOut3(fd, fileName)

    !> File ID
    integer, intent(in) :: fd

    !> Output file name
    character(*), intent(in) :: fileName

    close(fd)
    write(stdOut, "(2A)") 'MD information accumulated in ', fileName

  end subroutine writeMdOut3


  !> Write out charges.
  subroutine writeCharges(fCharges, fdCharges, orb, qInput, qBlockIn, qiBlockIn)

    !> File name for charges to be written to
    character(*), intent(in) :: fCharges

    !> File descriptor for charge output
    integer, intent(in) :: fdCharges

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> input charges
    real(dp), intent(in) :: qInput(:,:,:)

    !> Block populations if present
    real(dp), intent(in), optional :: qBlockIn(:,:,:,:)

    !> Imaginary part of block populations if present
    real(dp), intent(in), optional :: qiBlockIn(:,:,:,:)

    call writeQToFile(qInput, fCharges, fdCharges, orb, qBlockIn, qiBlockIn)
    write(stdOut, "(A,A)") '>> Charges saved for restart in ', trim(fCharges)

  end subroutine writeCharges


  !> Writes Hamiltonian and overlap matrices and stops program execution.
  subroutine writeHSAndStop(tWriteHS, tWriteRealHS, tRealHS, over, neighborList, nNeighbor,&
      & iAtomStart, iPair, img2CentCell, kPoint, iCellVec, cellVec, ham, iHam)

    !> Write dense hamiltonian and overlap matrices
    logical, intent(in) :: tWriteHS

    !> write sparse hamitonian and overlap matrices
    logical, intent(in) :: tWriteRealHS

    !> Is the hamitonian real?
    logical, intent(in) :: tRealHS

    !> overlap in sparse storage
    real(dp), intent(in) :: over(:)

    !> atomic neighbours
    type(TNeighborList), intent(in) :: neighborList

    !> number of neighbours for each central cell atom
    integer, intent(in) :: nNeighbor(:)

    !> dense matrix indexing for atomic blocks
    integer, intent(in) :: iAtomStart(:)

    !> sparse matrix indexing for atomic blocks
    integer, intent(in) :: iPair(:,:)

    !> Image atoms to central cell
    integer, intent(in) :: img2CentCell(:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> index  for which unit cell an atom is in
    integer, intent(in) :: iCellVec(:)

    !> vectors to unit cells, in lattice constant units
    real(dp), intent(in) :: cellVec(:,:)

    !> sparse hamitonian
    real(dp), intent(in) :: ham(:,:)

    !> imaginary part of hamitonian (used if allocated)
    real(dp), allocatable, intent(in) :: iHam(:,:)

    real(dp), allocatable :: hamUpDown(:,:)
    integer :: nSpin

    nSpin = size(ham, dim=2)

    ! Sanity check, although this should have been caught in initprogram already.
    if (nSpin == 4) then
      call error('Internal error: Hamiltonian writing for Pauli-Hamiltoninan not implemented')
    end if

    hamUpDown = ham
    call qm2ud(hamUpDown)

    ! Write out matrices if necessary and quit.
    if (allocated(iHam)) then
      call writeHS(tWriteHS, tWriteRealHS, tRealHS, hamUpDown, over, neighborList%iNeighbor,&
          & nNeighbor, iAtomStart, iPair, img2CentCell, kPoint, iCellVec, cellVec, iHam=iHam)
    else
      call writeHS(tWriteHS, tWriteRealHS, tRealHS, hamUpDown, over, neighborList%iNeighbor,&
          & nNeighbor, iAtomStart, iPair, img2CentCell, kPoint, iCellVec, cellVec)
    end if
    write(stdOut, "(A)") "Hamilton/Overlap written, exiting program."
    stop

  end subroutine writeHSAndStop


  !> Invokes the writing routines for the Hamiltonian and overlap matrices.
  subroutine writeHS(tWriteHS, tWriteRealHS, tRealHS, ham, over, iNeighbor, nNeighbor, iAtomStart,&
      & iPair, img2CentCell, kPoint, iCellVec, cellVec, iHam)

    !> Should the hamiltonian and overlap be written out as dense matrices
    logical, intent(in) :: tWriteHS

    !> Should the (sparse) real space storage hamiltonian and overlap
    logical, intent(in) :: tWriteRealHS

    !> Is the hamiltonian real?
    logical, intent(in) :: tRealHS

    !> sparse hamitonian matrix
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Atomic neighbour data
    integer, intent(in) :: iNeighbor(0:,:)

    !> number of atomic neighbours for each atom
    integer, intent(in) :: nNeighbor(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iAtomStart(:)

    !> Index array for start of atomic block in sparse matrices
    integer, intent(in) :: iPair(0:,:)

    !> Index array for images of atoms
    integer, intent(in) :: img2CentCell(:)

    !> The kpoints in the system
    real(dp), intent(in) :: kPoint(:,:)

    !> Index from atom to which unit cell it belongs
    integer, intent(in) :: iCellVec(:)

    !> Vectors to specific unit cells
    real(dp), intent(in) :: cellVec(:,:)

    !> Imaginary part of the hamiltonian if present
    real(dp), intent(in), optional :: iHam(:,:)

    integer :: iS, nSpin

    nSpin = size(ham, dim=2)

    if (tWriteRealHS) then
      do iS = 1, nSpin
        call writeSparse("hamreal" // i2c(iS) // ".dat", ham(:,iS), iNeighbor, &
            &nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, cellVec)
        if (present(iHam)) then
          call writeSparse("hamimag" // i2c(iS) // ".dat", iHam(:,iS),&
              & iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell,iCellVec,&
              & cellVec)
        end if
      end do
      call writeSparse("overreal.dat", over, iNeighbor, &
          &nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, cellVec)
    end if
    if (tWriteHS) then
      if (tRealHS) then
        do iS = 1, nSpin
          call writeSparseAsSquare("hamsqr" // i2c(iS) // ".dat", ham(:,iS), &
              &iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell)
        end do
        call writeSparseAsSquare("oversqr.dat", over, iNeighbor, nNeighbor, &
            &iAtomStart, iPair, img2CentCell)
      else
        do iS = 1, nSpin
          call writeSparseAsSquare("hamsqr" // i2c(iS) // ".dat", ham(:,iS), &
              &kPoint, iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, &
              &iCellVec, cellVec)
        end do
        call writeSparseAsSquare("oversqr.dat", over, kPoint, iNeighbor, &
            &nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, cellVec)
      end if
    end if

  end subroutine writeHS


  !> Writes the eigenvectors to disc.
  subroutine writeEigenvectors(nSpin, fdEigvec, runId, nAtom, neighborList, nNeighbor, cellVec,&
      & iCellVec, iAtomStart, iPair, img2CentCell, species, speciesName, orb, kPoint, over,&
      & HSqrReal, SSqrReal, HSqrCplx, SSqrCplx, storeEigvecsReal, storeEigvecsCplx)

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> File ID for ground state eigenvectors
    integer, intent(in) :: fdEigvec

    !> Job ID for future identification
    integer, intent(in) :: runId

    !> Number of real atoms in the system
    integer, intent(in) :: nAtom

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index of start of atom blocks in dense matrices
    integer, intent(in) :: iAtomStart(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iPair(:,:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> label for each atomic chemical species
    character(*), intent(in) :: speciesName(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Storage for dense hamiltonian matrix
    real(dp), intent(inout), optional :: HSqrReal(:,:,:)

    !> Storage for dense overlap matrix
    real(dp), intent(inout), optional :: SSqrReal(:,:)

    !> Storage for dense hamitonian matrix (complex case)
    complex(dp), intent(inout), optional :: HSqrCplx(:,:,:,:)

    !> Storage for dense overlap matrix (complex case)
    complex(dp), intent(inout), optional :: SSqrCplx(:,:)
    type(OFifoRealR2), intent(inout), optional :: storeEigvecsReal(:)
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecsCplx(:)

    @:ASSERT(present(HSqrReal) .neqv. present(HSqrReal))
    @:ASSERT(present(SSqrReal) .neqv. present(SSqrReal))
    @:ASSERT(.not. present(storeEigvecsReal) .or. present(HSqrReal))
    @:ASSERT(.not. present(storeEigvecsCplx) .or. present(HSqrCplx))

    if (present(HSqrCplx)) then
      !> Complex Pauli-Hamiltonian without k-points
      call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, nNeighbor, cellVec, iCellVec,&
          & iAtomStart, iPair, img2CentCell, orb, species, speciesName, over, kPoint, HSqrCplx,&
          & SSqrCplx, storeEigvecsCplx)
    else
      !> Real Hamiltonian
      call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, nNeighbor, iAtomStart, iPair,&
          & img2CentCell, orb, species, speciesName, over, HSqrReal, SSqrReal, storeEigvecsReal)
    end if

  end subroutine writeEigenvectors


  !> Write projected eigenvectors.
  subroutine writeProjectedEigenvectors(regionLabels, fdProjEig, eigen, nSpin, neighborList,&
      & nNeighbor, cellVec, iCellVec, iAtomStart, iPair, img2CentCell, orb, over, kPoint, kWeight,&
      & iOrbRegion, HSqrReal, SSqrReal, HSqrCplx, SSqrCplx, storeEigvecsReal, storeEigvecsCplx)
    type(ListCharLc), intent(inout) :: regionLabels
    integer, intent(in) :: fdProjEig(:)
    real(dp), intent(in) :: eigen(:,:,:)

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Index of start of atom blocks in dense matrices
    integer, intent(in) :: iAtomStart(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iPair(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)
    type(ListIntR1), intent(inout) :: iOrbRegion

    !> Storage for dense hamiltonian matrix
    real(dp), intent(inout), optional :: HSqrReal(:,:,:)

    !> Storage for dense overlap matrix
    real(dp), intent(inout), optional :: SSqrReal(:,:)

    !> Storage for dense hamitonian matrix (complex case)
    complex(dp), intent(inout), optional :: HSqrCplx(:,:,:,:)

    !> Storage for dense overlap matrix (complex case)
    complex(dp), intent(inout), optional :: SSqrCplx(:,:)
    type(OFifoRealR2), intent(inout), optional :: storeEigvecsReal(:)
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecsCplx(:)

    @:ASSERT(present(HSqrReal) .neqv. present(HSqrReal))
    @:ASSERT(present(SSqrReal) .neqv. present(SSqrReal))
    @:ASSERT(.not. present(storeEigvecsReal) .or. present(HSqrReal))
    @:ASSERT(.not. present(storeEigvecsCplx) .or. present(HSqrCplx))


    if (present(SSqrCplx)) then
      call writeProjEigvecs(regionLabels, fdProjEig, eigen, nSpin, neighborList, nNeighbor,&
          & cellVec, iCellVec, iAtomStart, iPair, img2CentCell, orb, over, kpoint, kWeight,&
          & HSqrCplx, SSqrCplx, iOrbRegion, storeEigvecsCplx)
    else
      call writeProjEigvecs(regionLabels, fdProjEig, eigen, nSpin, neighborList, nNeighbor,&
          & iAtomStart, iPair, img2CentCell, orb, over, HSqrReal, SSqrReal, iOrbRegion,&
          & storeEigvecsReal)
    end if

  end subroutine writeProjectedEigenvectors


  !> Write current geometry to disc
  subroutine writeCurrentGeometry(geoOutFile, pCoord0Out, tLatOpt, tMd, tAppendGeo, tFracCoord,&
      & tPeriodic, tPrintMulliken, species0, speciesName, latVec, iGeoStep, iLatGeoStep, nSpin,&
      & qOutput, velocities)

    !>  file for geometry output
    character(*), intent(in) :: geoOutFile

    !> How central cell atoms are represented
    real(dp), intent(in) :: pCoord0Out(:,:)

    !> is the lattice being optimised?
    logical, intent(in) :: tLatOpt

    !> Is this a molecular dynamics calculation?
    logical, intent(in) :: tMd

    !> should the geometry be added to the end, or the file cleared first
    logical, intent(in) :: tAppendGeo

    !> are fractional GEN files expected
    logical, intent(in) :: tFracCoord

    !> Is the geometry periodic?
    logical, intent(in) :: tPeriodic

    !> should Mulliken charges be printed
    logical, intent(in) :: tPrintMulliken

    !> species of atoms in the central cell
    integer, intent(in) :: species0(:)

    !> label for each atomic chemical species
    character(*), intent(in) :: speciesName(:)

    !> lattice vectors
    real(dp), intent(in) :: latVec(:,:)

    !> current geometry step
    integer, intent(in) :: iGeoStep

    !> current lattice step
    integer, intent(in) :: iLatGeoStep

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> charges
    real(dp), intent(in), optional :: qOutput(:,:,:)

    !> atomic velocities
    real(dp), intent(in), optional :: velocities(:,:)

    real(dp), allocatable :: tmpMatrix(:,:)
    integer :: nAtom
    integer :: ii, jj
    character(lc) :: comment, fname

    nAtom = size(pCoord0Out, dim=2)

    fname = trim(geoOutFile) // ".gen"
    if (tPeriodic) then
      call writeGenFormat(fname, pCoord0Out, species0, speciesName, latVec, tFracCoord)
    else
      call writeGenFormat(fname, pCoord0Out, species0, speciesName)
    end if

    fname = trim(geoOutFile) // ".xyz"
    if (tLatOpt) then
      write (comment, "(A, I0, A, I0)") '** Geometry step: ', iGeoStep, ', Lattice step: ',&
          & iLatGeoStep
    elseif (tMD) then
      write(comment, "(A, I0)") 'MD iter: ', iGeoStep
    else
      write(comment,"(A, I0)") 'Geometry Step: ', iGeoStep
    end if

    if (tPrintMulliken) then
      ! For non-colinear spin without velocities write magnetisation into the velocity field
      if (nSpin == 4 .and. .not. present(velocities)) then
        allocate(tmpMatrix(3, nAtom))
        do jj = 1, nAtom
          do ii = 1, 3
            tmpMatrix(ii,jj) = sum(qOutput(:, jj, ii + 1))
          end do
        end do
        ! convert by the inverse of the scaling used in writeXYZFormat
        tmpMatrix(:,:) = tmpMatrix * au__fs / (1000_dp * Bohr__AA)
        call writeXYZFormat(fname, pCoord0Out, species0, speciesName,&
            & charges=sum(qOutput(:,:,1), dim=1), velocities=tmpMatrix, comment=comment,&
            & append=tAppendGeo)
      else
        call writeXYZFormat(fname, pCoord0Out, species0, speciesName,&
            & charges=sum(qOutput(:,:,1),dim=1), velocities=velocities, comment=comment,&
            & append=tAppendGeo)
      end if
    else
      call writeXYZFormat(fname, pCoord0Out, species0, speciesName, velocities=velocities,&
          & comment=comment, append=tAppendGeo)
    end if

  end subroutine writeCurrentGeometry


  !> Write out final status of the geometry driver.
  subroutine writeFinalDriverStatus(tGeoOpt, tGeomEnd, tMd, tDerivs)

    !> Is the geometry being optimised?
    logical, intent(in) :: tGeoOpt

    !> Has the optimisation terminated?
    logical, intent(in) :: tGeomEnd

    !> Is this a molecular dynamics calculation?
    logical, intent(in) :: tMd

    !> Are finite difference derivatives being calculated?
    logical, intent(in) :: tDerivs

    if (tGeoOpt) then
      if (tGeomEnd) then
        write(stdOut, "(/, A)") "Geometry converged"
      else
        call warning("!!! Geometry did NOT converge!")
      end if
    elseif (tMD) then
      if (tGeomEnd) then
        write(stdOut, "(/, A)") "Molecular dynamics completed"
      else
        call warning("!!! Molecular dynamics terminated abnormally!")
      end if
    elseif (tDerivs) then
      if (tGeomEnd) then
        write(stdOut, "(/, A)") "Second derivatives completed"
      else
        call warning("!!! Second derivatives terminated abnormally!")
      end if
    end if

  end subroutine writeFinalDriverStatus


  !> Prints geometry step information to standard out
  subroutine printGeoStepInfo(tCoordOpt, tLatOpt, iLatGeoStep, iGeoStep)

    !> Are coordinates being optimised
    logical, intent(in) :: tCoordOpt

    !> Is the lattice being optimised
    logical, intent(in) :: tLatOpt

    !> Which geometry step is this
    integer, intent(in) :: iGeoStep

    !> How many lattice optimisation steps have occurred
    integer, intent(in) :: iLatGeoStep

    write(stdOut, '(/, A)') repeat('-', 80)
    if (tCoordOpt .and. tLatOpt) then
      write(stdOut, "(/, A, I0, A, I0,/)") '***  Geometry step: ', iGeoStep, ', Lattice step: ',&
          & iLatGeoStep
    else
      write(stdOut, "(/, A, I0, /)") '***  Geometry step: ', iGeoStep
    end if

  end subroutine printGeoStepInfo


  !> Prints the line above the start of the SCC cycle data
  subroutine printSccHeader()

    write(stdOut, "(A5, A18, A18, A18)") "iSCC", " Total electronic ", "  Diff electronic ",&
        & "     SCC error    "

  end subroutine printSccHeader


  !> Prints info about scc convergence.
  subroutine printSccInfo(tDftbU, iSccIter, Eelec, diffElec, sccErrorQ)

    !> Are orbital potentials being used
    logical, intent(in) :: tDftbU

    !> Iteration count
    integer, intent(in) :: iSccIter

    !> electronic energy
    real(dp), intent(in) :: Eelec

    !> Difference in electronic energy between this iteration and the last
    real(dp), intent(in) :: diffElec

    !> Maximum charge difference between input and output
    real(dp), intent(in) :: sccErrorQ

    if (tDFTBU) then
      write(stdOut, "(I5,E18.8,E18.8,E18.8)") iSCCIter, Eelec, diffElec, sccErrorQ
    else
      write(stdOut, "(I5,E18.8,E18.8,E18.8)") iSCCIter, Eelec, diffElec, sccErrorQ
    end if

  end subroutine printSccInfo


  !> Prints current total energies
  subroutine printEnergies(energy)

    !> energy components
    type(TEnergies), intent(in) :: energy

    write(stdOut, *)
    write(stdOut, format2U) "Total Energy", energy%Etotal,"H", Hartree__eV * energy%Etotal,"eV"
    write(stdOut, format2U) "Total Mermin free energy", energy%EMermin, "H",&
        & Hartree__eV * energy%EMermin," eV"

  end subroutine printEnergies


  !> Prints cell volume.
  subroutine printVolume(cellVol)

    !> unit cell volume
    real(dp), intent(in) :: cellVol

    write(stdOut, format2Ue) 'Volume', cellVol, 'au^3', (Bohr__AA**3) * cellVol, 'A^3'
  end subroutine printVolume


  !> Prints pressure and free energy.
  subroutine printPressureAndFreeEnergy(pressure, cellPressure, EGibbs)

    !> applied external pressure
    real(dp), intent(in) :: pressure

    !> internal cell pressure
    real(dp), intent(in) :: cellPressure

    !> Gibbs free energy (E -TS_elec +pV)
    real(dp), intent(in) :: EGibbs

    write(stdOut, format2Ue) 'Pressure', cellPressure, 'au', cellPressure * au__pascal, 'Pa'
    if (abs(pressure) > epsilon(1.0_dp)) then
      write(stdOut, format2U) "Gibbs free energy", EGibbs, 'H', Hartree__eV * EGibbs, 'eV'
    end if

  end subroutine printPressureAndFreeEnergy


  !> Writes maximal force component.
  subroutine printMaxForce(maxForce)

    !> maximum of the atomic forces
    real(dp), intent(in) :: maxForce

    write(stdOut, "(A, ':', T30, E20.6)") "Maximal force component", maxForce

  end subroutine printMaxForce


  !> Print maximal lattice force component
  subroutine printMaxLatticeForce(maxLattForce)

    !> Maximum energy derivative with respect to lattice vectors
    real(dp), intent(in) :: maxLattForce

    write(stdOut, format1Ue) "Maximal Lattice force component", maxLattForce, 'au'

  end subroutine printMaxLatticeForce


  !> Prints out info about current MD step.
  subroutine printMdInfo(tSetFillingTemp, tEField, tPeriodic, tempElec, absEField, tempIon,&
      & cellPressure, pressure, energy)

    !> Is the electronic temperature set by the thermostat method?
    logical, intent(in) :: tSetFillingTemp

    !> Is an electric field being applied?
    logical, intent(in) :: tEFIeld

    !> Is the geometry periodic?
    logical, intent(in) :: tPeriodic

    !> Electronic temperature
    real(dp), intent(in) :: tempElec

    !> magnitude of applied electric field
    real(dp), intent(in) :: absEField

    !> Atomic kinetic energy
    real(dp), intent(in) :: tempIon

    !> Internal pressure
    real(dp), intent(in) :: cellPressure

    !> External pressure (applied)
    real(dp), intent(in) :: pressure

    !> data type for energy components and total
    type(TEnergies), intent(in) :: energy

    if (tSetFillingTemp) then
      write(stdOut, format2U) 'Electronic Temperature:', tempElec, 'H', tempElec / Boltzmann, 'K'
    end if
    if (tEfield) then
      write(stdOut, format1U1e) 'External E field', absEField, 'au', absEField * au__V_m, 'V/m'
    end if
    write(stdOut, format2U) "MD Temperature:", tempIon, "H", tempIon / Boltzmann, "K"
    write(stdOut, format2U) "MD Kinetic Energy", energy%Ekin, "H", Hartree__eV * energy%Ekin, "eV"
    write(stdOut, format2U) "Total MD Energy", energy%EMerminKin, "H",&
        & Hartree__eV * energy%EMerminKin, "eV"
    if (tPeriodic) then
      write(stdOut, format2Ue) 'Pressure', cellPressure, 'au', cellPressure * au__pascal, 'Pa'
      if (abs(pressure) < epsilon(1.0_dp)) then
        write(stdOut, format2U) 'Gibbs free energy including KE', energy%EGibbsKin, 'H',&
            & Hartree__eV * energy%EGibbsKin, 'eV'
      end if
    end if

  end subroutine printMdInfo


  !> Receives the geometry from socket communication.
  subroutine receiveGeometryFromSocket(env, socket, tPeriodic, coord0, latVecs, tCoordsChanged,&
      & tLatticeChanged, tStopDriver)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Socket communication object
    type(IpiSocketComm), allocatable, intent(in) :: socket

    !> Is the system periodic
    logical, intent(in) :: tPeriodic

    !> Coordinates for atoms
    real(dp), intent(inout) :: coord0(:,:)

    !> Lattice vectors for the unit cell (not referenced if not periodic)
    real(dp), intent(inout) :: latVecs(:,:)

    !> Have the atomic coordinates changed
    logical, intent(out) :: tCoordsChanged

    !> Have the lattice vectors changed
    logical, intent(out) :: tLatticeChanged

    !> Stop the geometry driver if true
    logical, intent(out) :: tStopDriver

    real(dp) :: tmpLatVecs(3, 3)

    @:ASSERT(env%tIoProc .eqv. allocated(socket))

    if (env%tIoProc) then
      call socket%receive(coord0, tmpLatVecs, tStopDriver)
    end if
    tCoordsChanged = .true.
    if (tPeriodic .and. .not. tStopDriver) then
      latVecs(:,:) = tmpLatVecs
    end if
    tLatticeChanged = tPeriodic
  #:if WITH_MPI
    call mpifx_bcast(env%mpi%all, coord0)
    call mpifx_bcast(env%mpi%all, latVecs)
    call mpifx_bcast(env%mpi%all, tCoordsChanged)
    call mpifx_bcast(env%mpi%all, tLatticeChanged)
    call mpifx_bcast(env%mpi%all, tStopDriver)
  #:endif

  end subroutine receiveGeometryFromSocket




end module mainio
