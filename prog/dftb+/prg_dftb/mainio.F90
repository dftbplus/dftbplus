!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Various I/O routines for the main program.
module mainio
  use io
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
  use spin, only : qm2ud, ud2qm
  use energies
  use xmlf90
  use hsdutils, only : writeChildValue
  use mdintegrator, only : OMdIntegrator, state
  use formatout
  implicit none
  private

  public :: writeEigvecs, writeProjEigvecs, getH, SetEigVecsTxtOutput
  public :: initOutputFile, writeAutotestTag, writeResultsTag, writeDetailedXml, writeBandOut
  public :: writeHessianOut
  public :: writeDetailedOut1, writeDetailedOut2, writeDetailedOut3, writeDetailedOut4
  public :: writeDetailedOut5
  public :: writeMdOut1, writeMdOut2, writeMdOut3
  public :: writeHS, writeGenGeometry
  public :: printGeoStepInfo, printSccHeader
  public :: format1U, format2U, format1Ue, format2Ue, format1U1e


  !> output file names
  character(*), parameter :: eigvecOut = "eigenvec.out"
  character(*), parameter :: eigvecBin = "eigenvec.bin"

  character(len=*), parameter :: formatHessian = '(4f16.10)'
  character(len=*), parameter :: formatEnergy = '(8f12.5)'
  character(len=*), parameter :: formatEigen = '(8f14.8)'
  character(len=*), parameter :: formatGeoOut = "(I5, F16.8, F16.8, F16.8)"
  character(len=*), parameter :: format1U = "(A, ':', T32, F18.10, T51, A)"
  character(len=*), parameter :: format2U = "(A, ':', T32, F18.10, T51, A, T54, F16.4, T71, A)"
  character(len=*), parameter :: format1Ue = "(A, ':', T37, E13.6, T51, A)"
  character(len=*), parameter :: format2Ue = "(A, ':', T37, E13.6, T51, A, T57, E13.6, T71,A)"
  character(len=*), parameter :: format1U1e =&
      & "(' ', A, ':', T32, F18.10, T51, A, T57, E13.6, T71, A)"

  !> Private module variables (suffixed with "_" for clarity)
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

    !> Is a txt file written out, as well as the binary data?
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

    !> optional alternative file pre-fix
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
            !qState = real(sum(cVecTemp(iOrbs)), dp) &
            !    & + real(sum(cVecTemp(iOrbs+nOrb)), dp)
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


  subroutine initOutputFile(fileName, fd)
    character(*), intent(in) :: fileName
    integer, intent(out) :: fd

    fd = getFileId()
    open(fd, file=fileName, action="write", status="replace")
    close(fd)

  end subroutine initOutputFile


  subroutine writeAutotestTag(fd, fileName, tPeriodic, cellVol, tMulliken, qOutput, totalDeriv,&
      & chrgForces, excitedDerivs, tStress, totalStress, pDynMatrix, freeEnergy, pressure,&
      & gibbsFree, endCoords, tLocalise, localisation)
    integer, intent(in) :: fd
    character(*), intent(in) :: fileName
    logical, intent(in) :: tPeriodic
    real(dp), intent(in) :: cellVol
    logical, intent(in) :: tMulliken
    real(dp), intent(inout) :: qOutput(:,:,:)
    real(dp), allocatable, intent(in) :: totalDeriv(:,:)
    real(dp), allocatable, intent(in) :: chrgForces(:,:)
    real(dp), allocatable, intent(in) :: excitedDerivs(:,:)
    logical, intent(in) :: tStress
    real(dp), intent(in) :: totalStress(:,:)
    real(dp), pointer, intent(in) ::  pDynMatrix(:,:)
    real(dp), intent(in) :: freeEnergy
    real(dp), intent(in) :: pressure
    real(dp), intent(in) :: gibbsFree
    real(dp), intent(in) :: endCoords(:,:)
    logical, intent(in) :: tLocalise
    real(dp), intent(in) :: localisation

    open(fd, file=fileName, action="write", status="old", position="append")
    if (tPeriodic) then
      call writeTagged(fd, tag_volume, cellVol)
    end if
    if (tMulliken) then
      call qm2ud(qOutput)
      call writeTagged(fd, tag_qOutput, qOutput(:,:,1))
      call ud2qm(qOutput)
    end if
    if (allocated(totalDeriv)) then
      call writeTagged(fd, tag_forceTot, -totalDeriv)
      if (allocated(chrgForces)) then
        call writeTagged(fd, tag_chrgForces, -chrgForces)
      end if
      if (allocated(excitedDerivs)) then
        call writeTagged(fd, tag_excForce, -excitedDerivs)
      end if
      if (tStress) then
        call writeTagged(fd, tag_stressTot, totalStress)
      end if
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



  subroutine writeResultsTag(fd, fileName, totalDeriv, chrgForces, tStress, totalStress,&
      & pDynMatrix, tPeriodic, cellVol)
    integer, intent(in) :: fd
    character(*), intent(in) :: fileName
    real(dp), allocatable, intent(in) :: totalDeriv(:,:)
    real(dp), allocatable, intent(in) :: chrgForces(:,:)
    logical, intent(in) :: tStress
    real(dp), intent(in) :: totalStress(:,:)
    real(dp), pointer, intent(in) :: pDynMatrix(:,:)
    logical, intent(in) :: tPeriodic
    real(dp), intent(in) :: cellVol

    open(fd, file=fileName, action="write", status="replace")
    if (allocated(totalDeriv)) then
      call writeTagged(fd, tag_forceTot, -totalDeriv)
      if (allocated(chrgForces)) then
        call writeTagged(fd, tag_chrgForces, -chrgForces)
      end if
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


  subroutine writeDetailedXml(runId, speciesName, species0, coord0Out, tPeriodic, latVec, tRealHS,&
      & nKPoint, nSpin, nStates, nOrb, kPoint, kWeight, filling, occNatural)
    integer, intent(in) :: runId
    character(*), intent(in) :: speciesName(:)
    integer, intent(in) :: species0(:)
    real(dp), intent(in) :: coord0Out(:,:)
    logical, intent(in) :: tPeriodic
    real(dp), intent(in) :: latVec(:,:)
    logical, intent(in) :: tRealHS
    integer, intent(in) :: nKPoint
    integer, intent(in) :: nSpin
    integer, intent(in) :: nStates
    integer, intent(in) :: nOrb
    real(dp), intent(in) :: kPoint(:,:)
    real(dp), intent(in) :: kWeight(:)
    real(dp), intent(in) :: filling(:,:,:)
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


  subroutine writeBandOut(fd, fileName, eigen, filling, kWeight)
    integer, intent(in) :: fd
    character(*), intent(in) :: fileName
    real(dp), intent(in) :: eigen(:,:,:)
    real(dp), intent(in) :: filling(:,:,:)
    real(dp), intent(in) :: kWeight(:)

    integer :: iSpin, iK, iEgy

    open(unit=fd, file=fileName, action="write", status="replace")
    do iSpin = 1, size(eigen, dim=3)
      do iK = 1, size(eigen, dim=2)
        write(fd, *) 'KPT ', iK, ' SPIN ', iSpin, ' KWEIGHT ', kWeight(iK)
        do iEgy = 1, size(eigen, dim=1)
          write(fd, formatEnergy) Hartree__eV * eigen(iEgy, iK, iSpin), filling(iEgy, iK, iSpin)
        end do
        write(fd,*)
      end do
    end do
    close(fd)

  end subroutine writeBandOut


  subroutine writeHessianOut(fd, fileName, pDynMatrix)
    integer, intent(in) :: fd
    character(*), intent(in) :: fileName
    real(dp), intent(in) :: pDynMatrix(:,:)

    integer :: ii

    open(unit=fd, file=fileName, action="write", status="replace")
    do ii = 1, size(pDynMatrix, dim=2)
      write(fd, formatHessian) pDynMatrix(:, ii)
    end do
    close(fd)
    write(stdOut, "(2A)") 'Hessian matrix written to ', fileName

  end subroutine writeHessianOut


  subroutine writeDetailedOut1(fd, fileName, tAppendDetailedOut, iDistribFn, nGeoSteps, iGeoStep,&
      & tMD, tDerivs, tCoordOpt, tLatOpt, iLatGeoStep, iSccIter, energy, diffElec, sccErrorQ,&
      & indMovedAtom, coord0Out, q0, qInput, qOutput, eigen, filling, orb, species,&
      & tDFTBU, tImHam, tPrintMulliken, orbitalL, qBlockOut, Ef, Eband, TS, E0, pressure, cellVol,&
      & tAtomicEnergy, tDispersion, tEField, tPeriodic, nSpin, tSpinOrbit, tScc)
    integer, intent(in) :: fd
    character(*), intent(in) :: fileName
    logical, intent(in) :: tAppendDetailedOut
    integer, intent(in) :: iDistribFn
    integer, intent(in) :: nGeoSteps
    integer, intent(in) :: iGeoStep
    logical, intent(in) :: tMD, tDerivs, tCoordOpt, tLatOpt
    integer, intent(in) :: iLatGeoStep, iSccIter
    type(TEnergies), intent(in) :: energy
    real(dp), intent(in) :: diffElec, sccErrorQ
    integer, intent(in) :: indMovedAtom(:)
    real(dp), intent(in) :: coord0Out(:,:)
    real(dp), intent(in) :: q0(:,:,:)
    real(dp), intent(inout) :: qInput(:,:,:), qOutput(:,:,:)
    real(dp), intent(in) :: eigen(:,:,:), filling(:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)
    logical, intent(in) :: tDFTBU, tImHam, tPrintMulliken
    real(dp), allocatable, intent(in) :: orbitalL(:,:,:)
    real(dp), allocatable, intent(inout) :: qBlockOut(:,:,:,:)
    real(dp), intent(in) :: Ef(:), EBand(:), TS(:), E0(:)
    real(dp), intent(in) :: pressure, cellVol
    logical, intent(in) :: tAtomicEnergy, tDispersion, tEfield, tPeriodic
    integer, intent(in) :: nSpin
    logical, intent(in) :: tSpinOrbit, tScc

    real(dp) :: angularMomentum(3)
    integer :: ang
    integer :: nAtom, nLevel, nKPoint, nSpinHams, nMovedAtom
    integer :: iAt, iSpin, iEgy, iK, iSp, iSh, iOrb, kk
    logical :: tSpin

    nAtom = size(q0, dim=2)
    nLevel = size(eigen, dim=1)
    nKPoint = size(eigen, dim=2)
    nSpinHams = size(eigen, dim=3)
    nMovedAtom = size(indMovedAtom)
    tSpin = (nSpin == 2 .or. nSpin == 4)

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

    !! Write out atomic charges
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
        write(fd, formatEnergy) (filling(iEgy, iK, iSpin), iK = 1, nKPoint)
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
      if (nSpin == 2) then
        call qm2ud(qOutput)
        if (tDFTBU) then
          call qm2ud(qBlockOut)
        end if
      end if
      lpSpinPrint2: do iSpin = 1, nSpin
        if (tPrintMulliken) then
          write(fd, "(3A, F16.8)") 'Nr. of electrons (', trim(spinName(iSpin)), '):',&
              & sum(qOutput(:, :, iSpin))
          write(fd, "(3A)") 'Atom populations (', trim(spinName(iSpin)), ')'
          write(fd, "(A5, 1X, A16)") " Atom", " Population"
          do iAt = 1, nAtom
            write(fd, "(I5, 1X, F16.8)") iAt, sum(qOutput(:, iAt, iSpin))
          end do
          write(fd, *)
          write(fd, "(3A)") 'l-shell populations (', trim(spinName(iSpin)), ')'
          write(fd, "(A5, 1X, A3, 1X, A3, 1X, A16)")&
              & " Atom", "Sh.", "  l", " Population"
          do iAt = 1, nAtom
            iSp = species(iAt)
            do iSh = 1, orb%nShell(iSp)
              write(fd, "(I5, 1X, I3, 1X, I3, 1X, F16.8)") iAt, iSh, orb%angShell(iSh, iSp),&
                  & sum(qOutput(orb%posShell(iSh, iSp):orb%posShell(iSh + 1, iSp)-1, iAt,&
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
                    &iAt, iSh, ang, kk - ang, qOutput(orb%posShell(iSh, iSp) + kk, iAt, iSpin)
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
              write(fd, "(16F8.4)") qBlockOut(1:orb%nOrbSpecies(iSp), iOrb, iAt, iSpin)
            end do
          end do
          write(fd, *)
        end if
      end do lpSpinPrint2
      if (nSpin == 2) then
        call ud2qm(qOutput)
        if (tDFTBU) then
          call ud2qm(qBlockOut)
        end if
      end if
    end if

    if (nSpin == 2) then
      call qm2ud(qOutput)
      call qm2ud(qInput)
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
              & sum(qInput(:, :, iSpin)), sum(qOutput(:, :, iSpin))
        else
          write(fd, "(3A, 2F16.8)") 'Input / Output electrons (', quaternionName(iSpin), '):',&
              & sum(qInput(:, :, iSpin)), sum(qOutput(:, :, iSpin))
        end if
      end if
      write(fd, *)
    end do lpSpinPrint3

    if (nSpin == 2) then
      call ud2qm(qOutput)
      call ud2qm(qInput)
    end if

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



  subroutine writeDetailedOut2(fd, tScc, tConverged, tXlbomd, tLinResp, tGeoOpt, tMd, tPrintForces,&
      & tStress, tPeriodic, energy, totalStress, totalLatDeriv, totalDeriv, chrgForces,&
      & indMovedAtom, cellVol, cellPressure, geoOutFile)
    integer, intent(in) :: fd
    logical, intent(in) :: tScc, tConverged, tXlbomd, tLinResp, tGeoOpt, tMd
    logical, intent(in) :: tPrintForces, tStress, tPeriodic
    type(TEnergies), intent(in) :: energy
    real(dp), intent(in) :: totalStress(:,:), totalLatDeriv(:,:)
    real(dp), intent(in), allocatable :: totalDeriv(:,:), chrgForces(:,:)
    integer, intent(in) :: indMovedAtom(:)
    real(dp), intent(in) :: cellVol, cellPressure
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
      do iAt = 1, size(totalDeriv, dim=2)
        write(fd, "(3F20.12)") -totalDeriv(:, iAt)
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

      write(fd, format1Ue) "Maximal derivative component", maxval(abs(totalDeriv)), 'au'
      if (size(indMovedAtom) > 0) then
        write(fd, format1Ue) "Max force for moved atoms:",&
            & maxval(abs(totalDeriv(:, indMovedAtom))), 'au'
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


  subroutine writeDetailedOut3(fd, tPrintForces, tSetFillingTemp, tPeriodic, tStress, totalStress,&
      & totalLatDeriv, energy, tempElec, pressure, cellPressure, kT)
    integer, intent(in) :: fd
    logical, intent(in) :: tPrintForces, tSetFillingTemp, tPeriodic, tStress
    real(dp), intent(in) :: totalStress(:,:), totalLatDeriv(:,:)
    type(TEnergies), intent(in) :: energy
    real(dp), intent(in) :: tempElec, pressure, cellPressure, kT

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
    write(fd, format2U) "MD Temperature", kT, "H", kT / Boltzmann, "K"

  end subroutine writeDetailedOut3


  subroutine writeDetailedOut4(fd, tMd, energy, kT)
    integer, intent(in) :: fd
    logical, intent(in) :: tMd
    type(TEnergies), intent(in) :: energy
    real(dp), intent(in) :: kT

    if (tMd) then
      write(fd, format1U) "MD Kinetic Energy", energy%Ekin, "H"
      write(fd, format2U) "Total MD Energy", energy%EMerminKin, "H",&
          & Hartree__eV * energy%EMerminKin, "eV"
      write(fd, format2U) "MD Temperature", kT, "H", kT / Boltzmann, "K"
      write(fd, *)
    end if

  end subroutine writeDetailedOut4


  subroutine writeDetailedOut5(fd, tGeoOpt, tGeomEnd, tMd, tDerivs, tEField, absEField,&
      & dipoleMoment)
    integer, intent(in) :: fd
    logical, intent(in) :: tGeoOpt, tGeomEnd, tMd, tDerivs, tEField
    real(dp), intent(in) :: absEField
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


  subroutine writeMdOut1(fd, fileName, iGeoStep, pMdIntegrator)
    integer, intent(in) :: fd
    character(*), intent(in) :: fileName
    integer, intent(in) :: iGeoStep
    type(OMdIntegrator), intent(in) :: pMdIntegrator

    if (iGeoStep == 0) then
      open(fd, file=fileName, status="replace", action="write")
    end if
    write(fd, "(A, 1X, I0)") "MD step:", iGeoStep
    call state(pMdIntegrator, fd)

  end subroutine writeMdOut1


  subroutine writeMdOut2(fd, tStress, tBarostat, tLinResp, tEField, tFixEf, tPrintMulliken,&
      & energy, latVec, cellVol, cellPressure, pressure, kT, absEField, qOutput, q0, dipoleMoment)
    integer, intent(in) :: fd
    logical, intent(in) :: tStress, tBarostat, tLinResp, tEField, tFixEf, tPrintMulliken
    type(TEnergies), intent(in) :: energy
    real(dp), intent(in) :: latVec(:,:)
    real(dp), intent(in) :: cellVol, cellPressure, pressure, kT, absEField
    real(dp), intent(in) :: qOutput(:,:,:), q0(:,:,:)
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
    write(fd, format2U) 'MD Temperature', kT, 'au', kT / Boltzmann, 'K'
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


  subroutine writeMdOut3(fd, fileName)
    integer, intent(in) :: fd
    character(*), intent(in) :: fileName

    close(fd)
    write(stdOut, "(2A)") 'MD information accumulated in ', fileName

  end subroutine writeMdOut3


  !> Invokes the writing routines for the Hamiltonian and overlap matrices.
  subroutine writeHS(tWriteHS, tWriteRealHS, tRealHS, ham, over, iNeighbor, nNeighbor, iAtomStart,&
      & iPair, img2CentCell, kPoint, iCellVec, cellVec, iHam)
    logical, intent(in) :: tWriteHS, tWriteRealHS, tRealHS
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:), img2CentCell(:)
    real(dp), intent(in) :: kPoint(:,:)
    integer, intent(in) :: iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
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


  !> Write out geometry in gen format if needed
  subroutine writeGenGeometry(tGeoOpt, tMd, tWriteRestart, tFracCoord, tPeriodic, geoOutFile,&
      & coords0, species0, speciesNames, latVecs)
    logical, intent(in) :: tGeoOpt, tMd, tWriteRestart, tFracCoord, tPeriodic
    character(*), intent(in) :: geoOutFile
    real(dp), intent(in) :: coords0(:,:)
    integer, intent(in) :: species0(:)
    character(*), intent(in) :: speciesNames(:)
    real(dp), intent(in) :: latVecs(:,:)

    character(lc) :: lcTmpLocal

    if (tGeoOpt .or. tMD) then
      if (tWriteRestart) then
        write(lcTmpLocal, "(A, A)") trim(geoOutFile), ".gen"
        call clearFile(trim(lcTmpLocal))
        if (tPeriodic) then
          call writeGenFormat(trim(lcTmpLocal), coords0, species0, speciesNames, latVecs,&
              & tFracCoord)
        else
          call writeGenFormat(trim(lcTmpLocal), coords0, species0, speciesNames)
        end if
      end if
    end if

  end subroutine writeGenGeometry


  !> Prints geometry step information
  subroutine printGeoStepInfo(tCoordOpt, tLatOpt, iLatGeoStep, iGeoStep)
    logical, intent(in) :: tCoordOpt, tLatOpt
    integer, intent(in) :: iGeoStep, iLatGeoStep
    
    write(stdOut, '(/, A)') repeat('-', 80)
    if (tCoordOpt .and. tLatOpt) then
      write(stdOut, "(/, A, I0, A, I0,/)") '***  Geometry step: ', iGeoStep, ', Lattice step: ',&
          & iLatGeoStep
    else
      write(stdOut, "(/, A, I0, /)") '***  Geometry step: ', iGeoStep
    end if

  end subroutine printGeoStepInfo


  subroutine printSccHeader()

    write(stdOut, "(A5, A18, A18, A18)") "iSCC", " Total electronic ", "  Diff electronic ",&
        & "     SCC error    "
  end subroutine printSccHeader
  

end module mainio
