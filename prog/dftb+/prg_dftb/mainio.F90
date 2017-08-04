!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Various I/O routines for the main program.
module mainio
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
  implicit none
  private

  public :: writeEigvecs, writeProjEigvecs, getH, SetEigVecsTxtOutput

  !> output file names
  character(*), parameter :: eigvecOut = "eigenvec.out"
  character(*), parameter :: eigvecBin = "eigenvec.bin"
  character(*), parameter :: regionOut = "region_"

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

end module mainio
