!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#! Note: This module contains preprocessor variable substitutions in subroutine names (${NAME}$)
#! which may break the documentation system. Make sure you preprocess this file before passing it
#! to a source code documentation tool.

#:include 'common.fypp'

!> Various I/O routines for the main program.
module mainio
#:if WITH_MPI
  use mpifx
#:endif
#:if WITH_SCALAPACK
  use scalapackfx
#:endif
  use globalenv
  use environment
  use densedescr
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
#:if WITH_SOCKETS
  use ipisocket
#:endif
  implicit none
  private

  public :: writeEigenvectors, writeRealEigvecs, writeCplxEigvecs
#:if WITH_SCALAPACK
  public :: writeRealEigvecsBinBlacs, writeRealEigvecsTxtBlacs
  public :: writeCplxEigvecsBinBlacs, writeCplxEigvecsTxtBlacs
#:else
  public :: writeRealEigvecsBinSerial, writeRealEigvecsTxtSerial
  public :: writeCplxEigvecsBinSerial, writeCplxEigvecsTxtSerial
#:endif
  public :: writeProjectedEigenvectors
  public :: initOutputFile, writeAutotestTag, writeResultsTag, writeDetailedXml, writeBandOut
  public :: writeHessianOut
  public :: writeDetailedOut1, writeDetailedOut2, writeDetailedOut3, writeDetailedOut4
  public :: writeDetailedOut5
  public :: writeMdOut1, writeMdOut2, writeMdOut3
  public :: writeCharges
  public :: writeCurrentGeometry, writeFinalDriverStatus
  public :: writeHSAndStop, writeHS
  public :: printGeoStepInfo, printSccHeader, printSccInfo, printEnergies, printVolume
  public :: printPressureAndFreeEnergy, printMaxForce, printMaxLatticeForce
  public :: printMdInfo
#:if WITH_SOCKETS
  public :: receiveGeometryFromSocket
#:endif

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


contains

  !> Writes the eigenvectors to disc.
  subroutine writeEigenvectors(env, fd, runId, neighborList, nNeighbor, cellVec, iCellVec,&
      & denseDesc, iPair, img2CentCell, species, speciesName, orb, kPoint, over, parallelKS,&
      & tPrintEigvecsTxt, eigvecsReal, SSqrReal, eigvecsCplx, SSqrCplx)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> File ID for ground state eigenvectors
    integer, intent(in) :: fd

    !> Job ID for future identification
    integer, intent(in) :: runId

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index of start of atom blocks in dense matrices
    type(TDenseDescr), intent(in) :: denseDesc

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

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Whether eigenvectors should be also written in text form
    logical, intent(in) :: tPrintEigvecsTxt

    !> Real eigenvectors (will be overwritten)
    real(dp), intent(inout), allocatable :: eigvecsReal(:,:,:)

    !> Storage for dense real overlap matrix
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    !> Complex eigenvectors (will be overwritten)
    complex(dp), intent(inout), allocatable :: eigvecsCplx(:,:,:)

    !> Storage for dense complex overlap matrix
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    @:ASSERT(allocated(eigvecsReal) .neqv. allocated(eigvecsCplx))
    @:ASSERT(allocated(SSqrReal) .neqv. allocated(SSqrCplx))

    if (allocated(eigvecsCplx)) then
      call writeCplxEigvecs(env, fd, runId, neighborList, nNeighbor, cellVec, iCellVec, denseDesc,&
          & iPair, img2CentCell, species, speciesName, orb, kPoint, over, parallelKS,&
          & tPrintEigvecsTxt, eigvecsCplx, SSqrCplx)
    else
      call writeRealEigvecs(env, fd, runId, neighborList, nNeighbor, denseDesc, iPair,&
          & img2CentCell, species, speciesName, orb, over, parallelKS, tPrintEigvecsTxt,&
          & eigvecsReal, SSqrReal)
    end if

  end subroutine writeEigenvectors


  !> Writes real eigenvectors
  subroutine writeRealEigvecs(env, fd, runId, neighborList, nNeighbor, denseDesc, iPair,&
      & img2CentCell, species, speciesName, orb, over, parallelKS, tPrintEigvecsTxt, eigvecsReal,&
      & SSqrReal, fileName)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> File ID for ground state eigenvectors
    integer, intent(in) :: fd

    !> Job ID for future identification
    integer, intent(in) :: runId

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Index of start of atom blocks in dense matrices
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iPair(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> label for each atomic chemical species
    character(*), intent(in) :: speciesName(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Whether eigenvectors should be also written in text form
    logical, intent(in) :: tPrintEigvecsTxt

    !> Real eigenvectors (will be overwritten)
    real(dp), intent(inout) :: eigvecsReal(:,:,:)

    !> Storage for dense real overlap matrix
    real(dp), intent(out) :: SSqrReal(:,:)

    !> optional alternative file prefix, to appear as "fileName".bin or "fileName".out
    character(len=*), intent(in), optional :: fileName

  #:if WITH_SCALAPACK
    call writeRealEigvecsBinBlacs(env, denseDesc, eigvecsReal, fd, runId, parallelKS,&
        & fileName=fileName)
    if (tPrintEigvecsTxt) then
      call writeRealEigvecsTxtBlacs(env, denseDesc, eigvecsReal, fd, parallelKS, orb, over,&
          & neighborList%iNeighbor, nNeighbor, iPair, img2CentCell, species, speciesName,&
          & fileName=fileName)
    end if
  #:else
    call writeRealEigvecsBinSerial(eigvecsReal, fd, runId, parallelKS, fileName=fileName)
    if (tPrintEigvecsTxt) then
      call writeRealEigvecsTxtSerial(fd, neighborList, nNeighbor, denseDesc, iPair, img2CentCell,&
          & orb, species, speciesName, over, parallelKS, eigvecsReal, SSqrReal, fileName=fileName)
    end if
  #:endif

  end subroutine writeRealEigvecs


  !> Writes complex eigenvectors.
  subroutine writeCplxEigvecs(env, fd, runId, neighborList, nNeighbor, cellVec, iCellVec,&
      & denseDesc, iPair, img2CentCell, species, speciesName, orb, kPoint, over, parallelKS,&
      & tPrintEigvecsTxt, eigvecsCplx, SSqrCplx, fileName)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> File ID for ground state eigenvectors
    integer, intent(in) :: fd

    !> Job ID for future identification
    integer, intent(in) :: runId

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Index of start of atom blocks in dense matrices
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iPair(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

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

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Whether eigenvectors should be also written in text form
    logical, intent(in) :: tPrintEigvecsTxt

    !> Complex eigenvectors (will be overwritten)
    complex(dp), intent(inout) :: eigvecsCplx(:,:,:)

    !> Storage for dense complex overlap matrix
    complex(dp), intent(out) :: SSqrCplx(:,:)

    !> optional alternative file prefix, to appear as "fileName".bin or "fileName".out
    character(len=*), intent(in), optional :: fileName

  #:if WITH_SCALAPACK
    call writeCplxEigvecsBinBlacs(env, denseDesc, eigvecsCplx, fd, runId, parallelKS,&
        & fileName=fileName)
    if (tPrintEigvecsTxt) then
      if (denseDesc%t2Component) then
        call writePauliEigvecsTxtBlacs(env, denseDesc, eigvecsCplx, fd, parallelKS, orb, over,&
            & kPoint, neighborList%iNeighbor, nNeighbor, iCellVec, cellVec, iPair, img2CentCell,&
            & species, speciesName, fileName=fileName)
      else
        call writeCplxEigvecsTxtBlacs(env, denseDesc, eigvecsCplx, fd, parallelKS, orb, over,&
            & kPoint, neighborList%iNeighbor, nNeighbor, iCellVec, cellVec, iPair, img2CentCell,&
            & species, speciesName, fileName=fileName)
      end if
    end if
  #:else
    call writeCplxEigvecsBinSerial(eigvecsCplx, fd, runId, parallelKS, fileName=fileName)
    if (tPrintEigvecsTxt) then
      if (denseDesc%t2Component) then
        call writePauliEigvecsTxtSerial(fd, neighborList, nNeighbor, denseDesc, iPair,&
            & img2CentCell, iCellVec, cellVec, orb, species, speciesName, over, parallelKS, kPoint,&
            & eigvecsCplx, SSqrCplx, fileName=fileName)
      else
        call writeCplxEigvecsTxtSerial(fd, neighborList, nNeighbor, denseDesc, iPair, img2CentCell,&
            & iCellVec, cellVec, orb, species, speciesName, over, parallelKS, kPoint, eigvecsCplx,&
            & SSqrCplx, fileName=fileName)
      end if
    end if
  #:endif

  end subroutine writeCplxEigvecs


#:for DTYPE, NAME in [('complex', 'Cplx'), ('real', 'Real')]

#:if WITH_SCALAPACK

  !> Write the real eigvectors into binary output file (BLACS version).
  subroutine write${NAME}$EigvecsBinBlacs(env, denseDesc, eigvecs, fd, runId, parallelKS, fileName)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Square Hamiltonian (or work array)
    ${DTYPE}$(dp), intent(in) :: eigvecs(:,:,:)

    !> Fileid (file not yet opened) to use.
    integer, intent(in) :: fd

    !> Id of the current program run.
    integer, intent(in) :: runId

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> optional alternative file prefix, to appear as "fileName".bin
    character(len=*), intent(in), optional :: fileName

    type(linecomm) :: collector
    ${DTYPE}$(dp), allocatable :: localEigvec(:)
    integer :: nOrb
    integer :: iKS, iGroup, iEig

    nOrb = denseDesc%fullSize
    allocate(localEigvec(nOrb))

    if (env%mpi%tGlobalMaster) then
      call prepareEigvecFileBin(fd, runId, fileName)
    end if
    call collector%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")

    ! Master process collects in the first run (iGroup = 0) the columns of the matrix in its own
    ! process group (as process group master) via the collector. In the subsequent runs it just
    ! receives the columns collected by the respective group masters. The number of available
    ! matrices (possible k and s indices) may differ for various process groups. Also note, that
    ! the (k, s) pairs are round-robin distributed between the process groups.

    masterOrSlave: if (env%mpi%tGlobalMaster) then
      do iKS = 1, parallelKS%maxGroupKS
        group: do iGroup = 0, env%mpi%nGroup - 1
          if (iKS > parallelKS%nGroupKS(iGroup)) then
            cycle group
          end if
          do iEig = 1, nOrb
            if (iGroup == 0) then
              call collector%getline_master(env%blacs%orbitalGrid, iEig, eigvecs(:,:,iKS),&
                  & localEigvec)
            else
              call mpifx_recv(env%mpi%interGroupComm, localEigvec, iGroup)
            end if
            write(fd) localEigvec
          end do
        end do group
      end do
    else
      do iKS = 1, parallelKS%nLocalKS
        do iEig = 1, nOrb
          if (env%mpi%tGroupMaster) then
            call collector%getline_master(env%blacs%orbitalGrid, iEig, eigvecs(:,:,iKS),&
                & localEigvec)
            call mpifx_send(env%mpi%interGroupComm, localEigvec, env%mpi%interGroupComm%masterrank)
          else
            call collector%getline_slave(env%blacs%orbitalGrid, iEig, eigvecs(:,:,iKS))
          end if
        end do
      end do
    end if masterOrSlave

    if (env%mpi%tGlobalMaster) then
      close(fd)
    end if

  end subroutine write${NAME}$EigvecsBinBlacs

#:else

  !> Writes ${DTYPE}$ eigenvectors in binary format.
  subroutine write${NAME}$EigvecsBinSerial(eigvecs, fd, runId, parallelKS, fileName)

    !> Square Hamiltonian (or work array)
    ${DTYPE}$(dp), intent(in) :: eigvecs(:,:,:)

    !> Fileid (file not yet opened) to use.
    integer, intent(in) :: fd

    !> Id of the current program run.
    integer, intent(in) :: runId

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> optional alternative file prefix, to appear as "fileName".bin
    character(len=*), intent(in), optional :: fileName

    integer :: iKS, iSpin
    integer :: ii

    call prepareEigvecFileBin(fd, runId, fileName)
    do iKS = 1, parallelKS%nLocalKS
      iSpin = parallelKS%localKS(2, iKS)
      do ii = 1, size(eigvecs, dim=2)
        write(fd) eigvecs(:,ii,iSpin)
      end do
    end do
    close(fd)

  end subroutine write${NAME}$EigvecsBinSerial

#:endif

#:endfor


#:if WITH_SCALAPACK

  !> Write the real eigvectors into human readible output file (BLACS version).
  subroutine writeRealEigvecsTxtBlacs(env, denseDesc, eigvecs, fd, parallelKS, orb, over,&
      & iNeighbor, nNeighbor, iSparseStart, img2CentCell, species, speciesName, fileName)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Square Hamiltonian (or work array)
    real(dp), intent(in) :: eigvecs(:,:,:)

    !> Fileid (file not yet opened) to use.
    integer, intent(in) :: fd

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Sparse overlap
    real(dp), intent(in) :: over(:)

    !> Neighbors of each atom
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom
    integer, intent(in) :: nNeighbor(:)

    !> Index array for sparse matrices
    integer, intent(in) :: iSparseStart(:,:)

    !> Mapping of atoms into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> Name of each species
    character(*), intent(in) :: speciesName(:)

    !> optional alternative file prefix, to appear as "fileName".bin
    character(len=*), intent(in), optional :: fileName

    type(linecomm) :: collector
    real(dp), allocatable :: localEigvec(:), localFrac(:)
    real(dp), allocatable :: globalS(:,:), globalFrac(:,:)
    integer :: nOrb, nAtom
    integer :: iKS, iS, iGroup, iEig

    nOrb = denseDesc%fullSize
    nAtom = size(nNeighbor)
    allocate(globalS(size(eigvecs, dim=1), size(eigvecs, dim=2)))
    allocate(globalFrac(size(eigvecs, dim=1), size(eigvecs, dim=2)))
    if (env%mpi%tGroupMaster) then
      allocate(localEigvec(nOrb))
      allocate(localFrac(nOrb))
    end if

    call collector%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")
    if (env%mpi%tGlobalMaster) then
      call prepareEigvecFileTxt(fd, .false., fileName)
    end if

    ! See comment about algorithm in routine write${NAME}$EigvecsBinBlacs

    masterOrSlave: if (env%mpi%tGlobalMaster) then
      ! Global master process
      do iKS = 1, parallelKS%maxGroupKS
        group: do iGroup = 0, env%mpi%nGroup - 1
          if (iKS > parallelKS%nGroupKS(iGroup)) then
            cycle group
          end if
          iS = parallelKS%groupKS(2, iKS, iGroup)
          if (iGroup == 0) then
            call unpackHSRealBlacs(env%blacs, over, iNeighbor, nNeighbor, iSparseStart,&
                & img2CentCell, denseDesc, globalS)
            call pblasfx_psymm(globalS, denseDesc%blacsOrbSqr, eigvecs(:,:,iKS),&
                & denseDesc%blacsOrbSqr, globalFrac, denseDesc%blacsOrbSqr)
            globalFrac(:,:) = globalFrac * eigvecs(:,:,iKS)
          end if
          do iEig = 1, nOrb
            if (iGroup == 0) then
              call collector%getline_master(env%blacs%orbitalGrid, iEig, eigvecs(:,:,iKS),&
                  & localEigvec)
              call collector%getline_master(env%blacs%orbitalGrid, iEig, globalFrac, localFrac)
            else
              call mpifx_recv(env%mpi%interGroupComm, localEigvec, iGroup)
              call mpifx_recv(env%mpi%interGroupComm, localFrac, iGroup)
            end if
            call writeSingleRealEigvecTxt(fd, localEigvec, localFrac, iS, iEig, orb, species,&
                & speciesName, nAtom)
          end do
        end do group
      end do
    else
      ! All processes except the global master process
      do iKS = 1, parallelKS%nLocalKS
        call unpackHSRealBlacs(env%blacs, over, iNeighbor, nNeighbor, iSparseStart, img2CentCell,&
            & denseDesc, globalS)
        call pblasfx_psymm(globalS, denseDesc%blacsOrbSqr, eigvecs(:,:,iKS), denseDesc%blacsOrbSqr,&
            & globalFrac, denseDesc%blacsOrbSqr)
        globalFrac(:,:) = globalFrac * eigvecs(:,:,iKS)
        do iEig = 1, nOrb
          if (env%mpi%tGroupMaster) then
            call collector%getline_master(env%blacs%orbitalGrid, iEig, eigvecs(:,:,iKS),&
                & localEigvec)
            call collector%getline_master(env%blacs%orbitalGrid, iEig, globalFrac, localFrac)
            call mpifx_send(env%mpi%interGroupComm, localEigvec, env%mpi%interGroupComm%masterrank)
            call mpifx_send(env%mpi%interGroupComm, localFrac, env%mpi%interGroupComm%masterrank)
          else
            call collector%getline_slave(env%blacs%orbitalGrid, iEig, eigvecs(:,:,iKS))
            call collector%getline_slave(env%blacs%orbitalGrid, iEig, globalFrac)
          end if
        end do
      end do
    end if masterOrSlave

    if (env%mpi%tGlobalMaster) then
      close(fd)
    end if

  end subroutine writeRealEigvecsTxtBlacs

#:else

    !> Writes real eigenvectors in text form.
  subroutine writeRealEigvecsTxtSerial(fd, neighlist, nNeighbor, denseDesc, iPair, img2CentCell,&
      & orb, species, speciesName, over, parallelKS, eigvecs, SSqr, fileName)

    !> Fileid (file not yet opened) to use.
    integer, intent(in) :: fd

    !> Neighbor list.
    type(TNeighborList), intent(in) :: neighlist

    !> Nr. of neighbors for SK-interaction.
    integer, intent(in) :: nNeighbor(:)

    !> Dense descriptor for H and S
    type(TDenseDescr), intent(in) :: denseDesc

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

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Square Hamiltonian (or work array)
    real(dp), intent(inout) :: eigvecs(:,:,:)

    !> Work array.
    real(dp), intent(out) :: SSqr(:,:)

    !> optional alternative file pre-fix
    character(len=*), intent(in), optional :: fileName

    real(dp), allocatable :: rVecTemp(:)
    integer :: nAtom
    integer :: iKS, iS, iEig

    nAtom = size(nNeighbor)
    call prepareEigvecFileTxt(fd, .false., fileName)
    allocate(rVecTemp(size(eigvecs, dim=1)))
    call unpackHS(SSqr, over, neighlist%iNeighbor, nNeighbor, denseDesc%iAtomStart, iPair,&
        & img2CentCell)
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      do iEig = 1, denseDesc%nOrb
        call hemv(rVecTemp, SSqr, eigvecs(:,iEig,iS))
        rVecTemp = rVecTemp * eigvecs(:,iEig,iS)
        call writeSingleRealEigvecTxt(fd, eigvecs(:,iEig,iS), rVecTemp, iS, iEig, orb, species,&
            & speciesName, nAtom)
      end do
    end do
    close(fd)

  end subroutine writeRealEigvecsTxtSerial

#:endif


#:if WITH_SCALAPACK

  !> Write the complex eigvectors into human readible output file (BLACS version).
  subroutine writeCplxEigvecsTxtBlacs(env, denseDesc, eigvecs, fdEigvec, parallelKS, orb, over,&
      & kPoints, iNeighbor, nNeighbor, iCellVec, cellVec, iSparseStart, img2CentCell, species,&
      & speciesName, fileName)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Square Hamiltonian (or work array)
    complex(dp), intent(in) :: eigvecs(:,:,:)

    !> Fileid (file not yet opened) to use.
    integer, intent(in) :: fdEigvec

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Sparse overlap
    real(dp), intent(in) :: over(:)

    !> Kpoints
    real(dp), intent(in) :: kPoints(:,:)

    !> Neighbors of each atom
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom
    integer, intent(in) :: nNeighbor(:)

    !> Cell vector index for each atom
    integer, intent(in) :: iCellVec(:)

    !> Cell vectors
    real(dp), intent(in) :: cellVec(:,:)

    !> Index array for sparse matrices
    integer, intent(in) :: iSparseStart(:,:)

    !> Mapping of atoms into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> Name of each species
    character(*), intent(in) :: speciesName(:)

    !> optional alternative file prefix, to appear as "fileName".bin
    character(len=*), intent(in), optional :: fileName

    type(linecomm) :: collector
    complex(dp), allocatable :: localEigvec(:)
    complex(dp), allocatable :: globalS(:,:), globalSDotC(:,:)
    real(dp), allocatable :: localFrac(:), globalFrac(:,:)
    integer :: nEigvec, nAtom
    integer :: iKS, iK, iS, iGroup, iEig

    nEigvec = denseDesc%nOrb
    nAtom = size(nNeighbor)
    allocate(globalS(size(eigvecs, dim=1), size(eigvecs, dim=2)))
    allocate(globalSDotC(size(eigvecs, dim=1), size(eigvecs, dim=2)))
    allocate(globalFrac(size(eigvecs, dim=1), size(eigvecs, dim=2)))
    if (env%mpi%tGroupMaster) then
      allocate(localEigvec(denseDesc%nOrb))
      allocate(localFrac(denseDesc%nOrb))
    end if

    call collector%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")
    if (env%mpi%tGlobalMaster) then
      call prepareEigvecFileTxt(fdEigvec, .false., fileName)
    end if

    if (env%mpi%tGlobalMaster) then
      do iKS = 1, parallelKS%maxGroupKS
        group: do iGroup = 0, env%mpi%nGroup - 1
          if (iKS > parallelKS%nGroupKS(iGroup)) then
            cycle group
          end if
          iK = parallelKS%groupKS(1, iKS, iGroup)
          iS = parallelKS%groupKS(2, iKS, iGroup)
          if (iGroup == 0) then
            call unpackHSCplxBlacs(env%blacs, over, kPoints(:,iK), iNeighbor, nNeighbor, iCellVec,&
                & cellVec, iSparseStart, img2CentCell, denseDesc, globalS)
            call pblasfx_phemm(globalS, denseDesc%blacsOrbSqr, eigvecs(:,:,iKS),&
                & denseDesc%blacsOrbSqr, globalSDotC, denseDesc%blacsOrbSqr)
            globalFrac(:,:) = real(conjg(eigvecs(:,:,iKS)) * globalSDotC)
          end if
          do iEig = 1, nEigvec
            if (iGroup == 0) then
              call collector%getline_master(env%blacs%orbitalGrid, iEig, eigvecs(:,:,iKS),&
                  & localEigvec)
              call collector%getline_master(env%blacs%orbitalGrid, iEig, globalFrac, localFrac)
            else
              call mpifx_recv(env%mpi%interGroupComm, localEigvec, iGroup)
              call mpifx_recv(env%mpi%interGroupComm, localFrac, iGroup)
            end if
            call writeSingleCplxEigvecTxt(fdEigvec, localEigvec, localFrac, iS, iK, iEig, orb,&
                & species, speciesName, nAtom)
          end do
        end do group
      end do
    else
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        call unpackHSCplxBlacs(env%blacs, over, kPoints(:,iK), iNeighbor, nNeighbor, iCellVec,&
            & cellVec, iSparseStart, img2CentCell, denseDesc, globalS)
        call pblasfx_phemm(globalS, denseDesc%blacsOrbSqr, eigvecs(:,:,iKS),&
            & denseDesc%blacsOrbSqr, globalSDotC, denseDesc%blacsOrbSqr)
        globalFrac(:,:) = real(conjg(eigvecs(:,:,iKS)) * globalSDotC)
        do iEig = 1, nEigvec
          if (env%mpi%tGroupMaster) then
            call collector%getline_master(env%blacs%orbitalGrid, iEig, eigvecs(:,:,iKS),&
                & localEigvec)
            call collector%getline_master(env%blacs%orbitalGrid, iEig, globalFrac, localFrac)
            call mpifx_send(env%mpi%interGroupComm, localEigvec, env%mpi%interGroupComm%masterrank)
            call mpifx_send(env%mpi%interGroupComm, localFrac, env%mpi%interGroupComm%masterrank)
          else
            call collector%getline_slave(env%blacs%orbitalGrid, iEig, eigvecs(:,:,iKS))
            call collector%getline_slave(env%blacs%orbitalGrid, iEig, globalFrac)
          end if
        end do
      end do
    end if

    if (env%mpi%tGlobalMaster) then
      close(fdEigvec)
    end if

  end subroutine writeCplxEigvecsTxtBlacs

#:else

    !> Writes complex eigenvectors in text form.
  subroutine writeCplxEigvecsTxtSerial(fd, neighlist, nNeighbor, denseDesc, iPair, img2CentCell,&
      & iCellVec, cellVec, orb, species, speciesName, over, parallelKS, kPoints, eigvecs, SSqr,&
      & fileName)

    !> Fileid (file not yet opened) to use.
    integer, intent(in) :: fd

    !> Neighbor list.
    type(TNeighborList), intent(in) :: neighlist

    !> Nr. of neighbors for SK-interaction.
    integer, intent(in) :: nNeighbor(:)

    !> Dense matrix descriptor for H and S
    type(TDenseDescr), intent(in) :: denseDesc

    !> Positions of interactions in the sparse matrices.
    integer, intent(in) :: iPair(:,:)

    !> Mapping of atoms into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Index of cell vector mapping the atom to the central cell
    integer, intent(in) :: iCellVec(:)

    !> Cell vector coordinates
    real(dp), intent(in) :: cellVec(:,:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Species.
    integer, intent(in) :: species(:)

    !> Name of the species.
    character(mc), intent(in) :: speciesName(:)

    !> Sparse overlap matrix.
    real(dp), intent(in) :: over(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> K point coordinates
    real(dp), intent(in) :: kPoints(:,:)

    !> Square Hamiltonian (or work array)
    complex(dp), intent(inout) :: eigvecs(:,:,:)

    !> Work array.
    complex(dp), intent(inout) :: SSqr(:,:)

    !> optional alternative file pre-fix
    character(len=*), intent(in), optional :: fileName

    complex(dp), allocatable :: cVecTemp(:)
    real(dp), allocatable :: fracs(:)
    integer :: nEigvecs, nAtom
    integer :: iKS, iK, iS, iEig

    nAtom = size(nNeighbor)
    call prepareEigvecFileTxt(fd, denseDesc%t2Component, fileName)
    allocate(cVecTemp(size(eigvecs, dim=1)))
    nEigvecs = size(eigvecs, dim=2)
    allocate(fracs(nEigvecs))

    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iS = parallelKS%localKS(2, iKS)
      call unpackHS(SSqr, over, kPoints(:,iK), neighlist%iNeighbor, nNeighbor, iCellVec,&
          & cellVec, denseDesc%iAtomStart, iPair, img2CentCell)
      do iEig = 1, nEigvecs
        call hemv(cVecTemp, SSqr, eigvecs(:,iEig,iKS))
        fracs(:) = real(conjg(eigvecs(:,iEig,iKS)) * cVecTemp)
        call writeSingleCplxEigvecTxt(fd, eigvecs(:,iEig,iKS), fracs, iS, iK, iEig, orb, species,&
            & speciesName, nAtom)
      end do
    end do
    close(fd)

  end subroutine writeCplxEigvecsTxtSerial

#:endif


#:if WITH_SCALAPACK

  !> Write the complex eigvectors into human readible output file (BLACS version).
  subroutine writePauliEigvecsTxtBlacs(env, denseDesc, eigvecs, fd, parallelKS, orb, over, kPoints,&
      & iNeighbor, nNeighbor, iCellVec, cellVec, iSparseStart, img2CentCell, species, speciesName,&
      & fileName)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Square Hamiltonian (or work array)
    complex(dp), intent(in) :: eigvecs(:,:,:)

    !> Fileid (file not yet opened) to use.
    integer, intent(in) :: fd

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Sparse overlap
    real(dp), intent(in) :: over(:)

    !> Kpoints
    real(dp), intent(in) :: kPoints(:,:)

    !> Neighbors of each atom
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom
    integer, intent(in) :: nNeighbor(:)

    !> Cell vector index for each atom
    integer, intent(in) :: iCellVec(:)

    !> Cell vectors
    real(dp), intent(in) :: cellVec(:,:)

    !> Index array for sparse matrices
    integer, intent(in) :: iSparseStart(:,:)

    !> Mapping of atoms into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> Name of each species
    character(*), intent(in) :: speciesName(:)

    !> optional alternative file prefix, to appear as "fileName".bin
    character(len=*), intent(in), optional :: fileName

    type(linecomm) :: collector
    real(dp), allocatable :: fracs(:,:)
    complex(dp), allocatable :: localEigvec(:), localSDotC(:)
    complex(dp), allocatable :: globalS(:,:), globalSDotC(:,:)
    integer :: nAtom, nOrb
    integer :: iKS, iK, iGroup, iEig

    nOrb = denseDesc%fullSize
    nAtom = size(nNeighbor)
    allocate(globalS(size(eigvecs, dim=1), size(eigvecs, dim=2)))
    allocate(globalSDotC(size(eigvecs, dim=1), size(eigvecs, dim=2)))
    if (env%mpi%tGroupMaster) then
      allocate(localEigvec(denseDesc%fullSize))
      allocate(localSDotC(denseDesc%fullSize))
      if (env%mpi%tGlobalMaster) then
        allocate(fracs(4, denseDesc%nOrb))
      end if
    end if

    call collector%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")
    if (env%mpi%tGlobalMaster) then
      call prepareEigvecFileTxt(fd, .true., fileName)
    end if

    ! See comment about algorithm in routine write${NAME}$EigvecsBinBlacs

    masterOrSlave: if (env%mpi%tGlobalMaster) then
      ! Global master process
      do iKS = 1, parallelKS%maxGroupKS
        group: do iGroup = 0, env%mpi%nGroup - 1
          if (iKS > parallelKS%nGroupKS(iGroup)) then
            cycle group
          end if
          iK = parallelKS%groupKS(1, iKS, iGroup)
          if (iGroup == 0) then
            call unpackSPauliBlacs(env%blacs, over, kPoints(:,iK), iNeighbor, nNeighbor, iCellVec,&
                & cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc, globalS)
            call pblasfx_phemm(globalS, denseDesc%blacsOrbSqr, eigvecs(:,:,iKS),&
                & denseDesc%blacsOrbSqr, globalSDotC, denseDesc%blacsOrbSqr)
          end if
          do iEig = 1, nOrb
            if (iGroup == 0) then
              call collector%getline_master(env%blacs%orbitalGrid, iEig, eigvecs(:,:,iKS),&
                  & localEigvec)
              call collector%getline_master(env%blacs%orbitalGrid, iEig, globalSDotC, localSDotC)
            else
              call mpifx_recv(env%mpi%interGroupComm, localEigvec, iGroup)
              call mpifx_recv(env%mpi%interGroupComm, localSDotC, iGroup)
            end if
            call getPauliFractions(localEigvec, localSDotC, fracs)
            call writeSinglePauliEigvecTxt(fd, localEigvec, fracs, iK, iEig, orb, species,&
                & speciesName, nAtom, denseDesc%nOrb)
          end do
        end do group
      end do
    else
      ! All processes except the global master process
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        call unpackSPauliBlacs(env%blacs, over, kPoints(:,iK), iNeighbor, nNeighbor, iCellVec,&
            & cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc, globalS)
        call pblasfx_phemm(globalS, denseDesc%blacsOrbSqr, eigvecs(:,:,iKS),&
            & denseDesc%blacsOrbSqr, globalSDotC, denseDesc%blacsOrbSqr)
        do iEig = 1, nOrb
          if (env%mpi%tGroupMaster) then
            call collector%getline_master(env%blacs%orbitalGrid, iEig, eigvecs(:,:,iKS),&
                & localEigvec)
            call collector%getline_master(env%blacs%orbitalGrid, iEig, globalSDotC, localSDotC)
            call mpifx_send(env%mpi%interGroupComm, localEigvec, env%mpi%interGroupComm%masterrank)
            call mpifx_send(env%mpi%interGroupComm, localSDotC, env%mpi%interGroupComm%masterrank)
          else
            call collector%getline_slave(env%blacs%orbitalGrid, iEig, eigvecs(:,:,iKS))
            call collector%getline_slave(env%blacs%orbitalGrid, iEig, globalSDotC)
          end if
        end do
      end do
    end if masterOrSlave

    if (env%mpi%tGlobalMaster) then
      close(fd)
    end if

  end subroutine writePauliEigvecsTxtBlacs

#:else

    !> Writes complex eigenvectors in text form.
  subroutine writePauliEigvecsTxtSerial(fd, neighlist, nNeighbor, denseDesc, iPair, img2CentCell,&
      & iCellVec, cellVec, orb, species, speciesName, over, parallelKS, kPoints, eigvecs, SSqr,&
      & fileName)

    !> Fileid (file not yet opened) to use.
    integer, intent(in) :: fd

    !> Neighbor list.
    type(TNeighborList), intent(in) :: neighlist

    !> Nr. of neighbors for SK-interaction.
    integer, intent(in) :: nNeighbor(:)

    !> Dense matrix descriptor for H and S
    type(TDenseDescr), intent(in) :: denseDesc

    !> Positions of interactions in the sparse matrices.
    integer, intent(in) :: iPair(:,:)

    !> Mapping of atoms into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Index of cell vector mapping the atom to the central cell
    integer, intent(in) :: iCellVec(:)

    !> Cell vector coordinates
    real(dp), intent(in) :: cellVec(:,:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Species.
    integer, intent(in) :: species(:)

    !> Name of the species.
    character(mc), intent(in) :: speciesName(:)

    !> Sparse overlap matrix.
    real(dp), intent(in) :: over(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> K point coordinates
    real(dp), intent(in) :: kPoints(:,:)

    !> Square Hamiltonian (or work array)
    complex(dp), intent(inout) :: eigvecs(:,:,:)

    !> Work array.
    complex(dp), intent(inout) :: SSqr(:,:)

    !> optional alternative file pre-fix
    character(len=*), intent(in), optional :: fileName

    complex(dp), allocatable :: cVecTemp(:)
    real(dp), allocatable :: fracs(:,:)
    integer :: nEigvecs, nAtom
    integer :: iKS, iK, iEig

    nAtom = size(nNeighbor)
    call prepareEigvecFileTxt(fd, denseDesc%t2Component, fileName)
    allocate(cVecTemp(size(eigvecs, dim=1)))
    nEigvecs = size(eigvecs, dim=2)
    allocate(fracs(4, nEigvecs / 2))

    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      call unpackSPauli(over, kPoints(:,iK), neighlist%iNeighbor, nNeighbor,&
          & denseDesc%iAtomStart, iPair, img2CentCell, iCellVec, cellVec, SSqr)
      do iEig = 1, nEigvecs
        call hemv(cVecTemp, SSqr, eigvecs(:,iEig,iKS))
        call getPauliFractions(eigvecs(:,iEig,iKS), cVecTemp, fracs)
        call writeSinglePauliEigvecTxt(fd, eigvecs(:,iEig,iKS), fracs, iK, iEig, orb,&
            & species, speciesName, nAtom, denseDesc%nOrb)
      end do
    end do
    close(fd)

  end subroutine writePauliEigvecsTxtSerial

#:endif


  !> Write projected eigenvectors.
  subroutine writeProjectedEigenvectors(env, regionLabels, fd, eigen, neighborList, nNeighbor,&
      & cellVec, iCellVec, denseDesc, iPair, img2CentCell, orb, over, kPoint, kWeight, iOrbRegion,&
      & parallelKS, eigvecsReal, workReal, eigvecsCplx, workCplx)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> File name prefix for each region
    type(ListCharLc), intent(inout) :: regionLabels

    !> File descriptor for each region
    integer, intent(in) :: fd(:)

    !> Eigenvalues
    real(dp), intent(in) :: eigen(:,:,:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Dense matrix descriptor for H and S
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iPair(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)
    type(ListIntR1), intent(inout) :: iOrbRegion

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Storage for eigenvectors (real)
    real(dp), intent(inout), allocatable :: eigvecsReal(:,:,:)

    !> Work space (real)
    real(dp), intent(inout), allocatable :: workReal(:,:)

    !> Storage for eigenvectors (complex)
    complex(dp), intent(inout), allocatable :: eigvecsCplx(:,:,:)

    !> Work space (complex)
    complex(dp), intent(inout), allocatable :: workCplx(:,:)

    @:ASSERT(allocated(eigvecsReal) .neqv. allocated(eigvecsCplx))
    @:ASSERT(allocated(workReal) .neqv. allocated(workCplx))

  #:if WITH_SCALAPACK
    if (allocated(eigvecsCplx)) then
      if (denseDesc%t2Component) then
        call writeProjPauliEigvecsBlacs(env, denseDesc, regionLabels, fd, iOrbRegion, eigen,&
            & eigvecsCplx, orb, parallelKS, kPoint, kWeight, over, neighborList, nNeighbor, iPair,&
            & img2CentCell, iCellVec, cellVec)
      else
        call writeProjCplxEigvecsBlacs(env, denseDesc, regionLabels, fd, iOrbRegion, eigen&
            &,eigvecsCplx, parallelKS, kPoint, kWeight, over, neighborList, nNeighbor, iPair,&
            & img2CentCell, iCellVec, cellVec)
      end if
    else
      call writeProjRealEigvecsBlacs(env, denseDesc, regionLabels, fd, iOrbRegion, eigen,&
          & eigvecsReal, parallelKS, over, neighborList, nNeighbor, iPair, img2CentCell)
    end if
  #:else
    if (allocated(eigvecsCplx)) then
      if (denseDesc%t2Component) then
        call writeProjPauliEigvecsSerial(regionLabels, fd, eigen, neighborList, nNeighbor, cellVec,&
            & iCellVec, denseDesc, iPair, img2CentCell, over, kpoint, kWeight, parallelKS,&
            & eigvecsCplx, workCplx, iOrbRegion)
      else
        call writeProjCplxEigvecsSerial(regionLabels, fd, eigen, neighborList, nNeighbor, cellVec,&
            & iCellVec, denseDesc, iPair, img2CentCell, over, kpoint, kWeight, parallelKS,&
            & eigvecsCplx, workCplx, iOrbRegion)
      end if
    else
      call writeProjRealEigvecsSerial(regionLabels, fd, eigen, neighborList, nNeighbor, denseDesc,&
          & iPair, img2CentCell, over, parallelKS, eigvecsReal, workReal, iOrbRegion)
    end if
  #:endif

  end subroutine writeProjectedEigenvectors


#:if WITH_SCALAPACK

  !> Write the real eigvectors into human readible output file (BLACS version).
  subroutine writeProjRealEigvecsBlacs(env, denseDesc, fileNames, fd, iOrbRegion, eigvals,&
      & eigvecs, parallelKS, over, neighborList, nNeighbor, iSparseStart, img2CentCell)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> List of region file names
    type(ListCharLc), intent(inout) :: fileNames

    !> File descriptor for a not yet opened file for  each region
    integer, intent(in) :: fd(:)

    !> orbital number in each region
    type(listIntR1), intent(inout) :: iOrbRegion

    !> Eigenvalues
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Eigenvectors
    real(dp), intent(in) :: eigvecs(:,:,:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Sparse overlap
    real(dp), intent(in) :: over(:)

    !> Neighbors of each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Nr. of neighbors for each atom
    integer, intent(in) :: nNeighbor(:)

    !> Index array for sparse matrices
    integer, intent(in) :: iSparseStart(:,:)

    !> Mapping of atoms into the central cell.
    integer, intent(in) :: img2CentCell(:)

    type(linecomm) :: collector
    real(dp), allocatable :: globalS(:,:), globalFrac(:,:), localFrac(:)
    integer :: nOrb, nReg
    integer :: iKS, iS, iGroup, iEig

    nReg = size(fd)
    nOrb = denseDesc%fullSize
    allocate(globalS(size(eigvecs, dim=1), size(eigvecs, dim=2)))
    allocate(globalFrac(size(eigvecs, dim=1), size(eigvecs, dim=2)))
    if (env%mpi%tGroupMaster) then
      allocate(localFrac(nOrb))
    end if

    if (env%mpi%tGlobalMaster) then
      call prepareProjEigvecFiles(fd, fileNames)
    end if
    call collector%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")

    ! See comment about algorithm in routine write${NAME}$EigvecsBinBlacs

    masterOrSlave: if (env%mpi%tGlobalMaster) then
      ! Global master process
      do iKS = 1, parallelKS%maxGroupKS
        group: do iGroup = 0, env%mpi%nGroup - 1
          if (iKS > parallelKS%nGroupKS(iGroup)) then
            cycle group
          end if
          iS = parallelKS%groupKS(2, iKS, iGroup)
          if (iGroup == 0) then
            call unpackHSRealBlacs(env%blacs, over, neighborList%iNeighbor, nNeighbor,&
                & iSparseStart, img2CentCell, denseDesc, globalS)
            call pblasfx_psymm(globalS, denseDesc%blacsOrbSqr, eigvecs(:,:,iKS),&
                & denseDesc%blacsOrbSqr, globalFrac, denseDesc%blacsOrbSqr)
            globalFrac(:,:) = eigvecs(:,:,iKS) * globalFrac
          end if
          call writeProjEigvecHeader(fd, iS)
          do iEig = 1, nOrb
            if (iGroup == 0) then
              call collector%getline_master(env%blacs%orbitalGrid, iEig, globalFrac, localFrac)
            else
              call mpifx_recv(env%mpi%interGroupComm, localFrac, iGroup)
            end if
            call writeProjEigvecData(fd, iOrbRegion, eigvals(iEig, 1, iS), localFrac)
          end do
          call writeProjEigvecFooter(fd)
        end do group
      end do
    else
      ! All processes except the global master process
      do iKS = 1, parallelKS%nLocalKS
        call unpackHSRealBlacs(env%blacs, over, neighborList%iNeighbor, nNeighbor, iSparseStart,&
            & img2CentCell, denseDesc, globalS)
        call pblasfx_psymm(globalS, denseDesc%blacsOrbSqr, eigvecs(:,:,iKS), denseDesc%blacsOrbSqr,&
            & globalFrac, denseDesc%blacsOrbSqr)
        globalFrac(:,:) = eigvecs(:,:,iKS) * globalFrac
        do iEig = 1, nOrb
          if (env%mpi%tGroupMaster) then
            call collector%getline_master(env%blacs%orbitalGrid, iEig, globalFrac, localFrac)
            call mpifx_send(env%mpi%interGroupComm, localFrac, env%mpi%interGroupComm%masterrank)
          else
            call collector%getline_slave(env%blacs%orbitalGrid, iEig, globalFrac)
          end if
        end do
      end do
    end if masterOrSlave

    if (env%mpi%tGlobalMaster) then
      call finishProjEigvecFiles(fd)
    end if

  end subroutine writeProjRealEigvecsBlacs

#:else

    !> Write the projected eigenstates into text files
  subroutine writeProjRealEigvecsSerial(fileNames, fd, eigvals, neighlist, nNeighbor, denseDesc,&
      & iPair, img2CentCell, over, parallelKS, eigvecs, work, iOrbRegion)

    !> List with fileNames for each region
    type(listCharLc), intent(inout) :: fileNames

    !> File unit IDs for each of the regions
    integer, intent(in) :: fd(:)

    !> eigenvalues
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Neighbor list
    type(TNeighborList), intent(in) :: neighlist

    !> Nr. of neighbors for SK-interaction
    integer, intent(in) :: nNeighbor(:)

    !> Dense matrix descriptor for H and S
    type(TDenseDescr), intent(in) :: denseDesc

    !> Positions of interactions in the sparse matrices
    integer, intent(in) :: iPair(:,:)

    !> Mapping of atoms into the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Square Hamiltonian (or work array)
    real(dp), intent(inout) :: eigvecs(:,:,:)

    !> Work array
    real(dp), intent(out) :: work(:,:)

    !> orbital number in each region
    type(listIntR1), intent(inout) :: iOrbRegion

    integer :: iKS, iS, iEig
    real(dp), allocatable :: rVecTemp(:)

    call prepareProjEigvecFiles(fd, fileNames)

    allocate(rVecTemp(size(eigvecs, dim=1)))
    call unpackHS(work, over, neighlist%iNeighbor, nNeighbor, denseDesc%iAtomStart, iPair,&
        & img2CentCell)
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      call writeProjEigvecHeader(fd, iS)
      do iEig = 1, denseDesc%nOrb
        call hemv(rVecTemp, work, eigvecs(:,iEig,iKS))
        rVecTemp(:) = rVecTemp * eigvecs(:,iEig,iKS)
        call writeProjEigvecData(fd, iOrbRegion, eigvals(iEig, 1, iS), rVecTemp)
      end do
      call writeProjEigvecFooter(fd)
    end do

  end subroutine writeProjRealEigvecsSerial

#:endif


#:if WITH_SCALAPACK

  !> Write the complex eigvectors into human readible output file (BLACS version).
  subroutine writeProjCplxEigvecsBlacs(env, denseDesc, fileNames, fd, iOrbRegion, eigvals,&
      & eigvecs, parallelKS, kPoints, kWeights, over, neighborList, nNeighbor, iSparseStart,&
      & img2CentCell, iCellVec, cellVec)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> List of region file names
    type(ListCharLc), intent(inout) :: fileNames

    !> File descriptor for a not yet opened file for  each region
    integer, intent(in) :: fd(:)

    !> orbital number in each region
    type(listIntR1), intent(inout) :: iOrbRegion

    !> Eigenvalues
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Eigenvectors
    complex(dp), intent(in) :: eigvecs(:,:,:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> K-points
    real(dp), intent(in) :: kPoints(:,:)

    !> Weights of the k-points
    real(dp), intent(in) :: kWeights(:)

    !> Sparse overlap
    real(dp), intent(in) :: over(:)

    !> Neighbors of each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Nr. of neighbors for each atom
    integer, intent(in) :: nNeighbor(:)

    !> Index array for sparse matrices
    integer, intent(in) :: iSparseStart(:,:)

    !> Mapping of atoms into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Cell vector index for each atom
    integer, intent(in) :: iCellVec(:)

    !> Cell vectors
    real(dp), intent(in) :: cellVec(:,:)

    type(linecomm) :: collector
    real(dp), allocatable :: globalFrac(:,:), localFrac(:)
    complex(dp), allocatable :: globalS(:,:), globalSDotC(:,:)
    integer :: nOrb, nReg
    integer :: iKS, iK, iS, iGroup, iEig

    nReg = size(fd)
    nOrb = denseDesc%fullSize
    allocate(globalS(size(eigvecs, dim=1), size(eigvecs, dim=2)))
    allocate(globalSDotC(size(eigvecs, dim=1), size(eigvecs, dim=2)))
    allocate(globalFrac(size(eigvecs, dim=1), size(eigvecs, dim=2)))
    if (env%mpi%tGroupMaster) then
      allocate(localFrac(nOrb))
    end if

    if (env%mpi%tGlobalMaster) then
      call prepareProjEigvecFiles(fd, fileNames)
    end if
    call collector%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")

    ! See comment about algorithm in routine write${NAME}$EigvecsBinBlacs

    masterOrSlave: if (env%mpi%tGlobalMaster) then
      ! Global master process
      do iKS = 1, parallelKS%maxGroupKS
        group: do iGroup = 0, env%mpi%nGroup - 1
          if (iKS > parallelKS%nGroupKS(iGroup)) then
            cycle group
          end if
          iK = parallelKS%groupKS(1, iKS, iGroup)
          iS = parallelKS%groupKS(2, iKS, iGroup)
          if (iGroup == 0) then
            call unpackHSCplxBlacs(env%blacs, over, kPoints(:,iK), neighborList%iNeighbor,&
                & nNeighbor, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, globalS)
            call pblasfx_phemm(globalS, denseDesc%blacsOrbSqr, eigvecs(:,:,iKS),&
                & denseDesc%blacsOrbSqr, globalSDotC, denseDesc%blacsOrbSqr)
            globalFrac(:,:) = real(globalSDotC * conjg(eigvecs(:,:,iKS)))
          end if
          call writeProjEigvecHeader(fd, iS, iK, kWeights(iK))
          do iEig = 1, nOrb
            if (iGroup == 0) then
              call collector%getline_master(env%blacs%orbitalGrid, iEig, globalFrac, localFrac)
            else
              call mpifx_recv(env%mpi%interGroupComm, localFrac, iGroup)
            end if
            call writeProjEigvecData(fd, iOrbRegion, eigvals(iEig, iK, iS), localFrac)
          end do
        end do group
        call writeProjEigvecFooter(fd)
      end do
    else
      ! All processes except the global master process
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        call unpackHSCplxBlacs(env%blacs, over, kPoints(:,iK), neighborList%iNeighbor, nNeighbor,&
            & iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, globalS)
        call pblasfx_phemm(globalS, denseDesc%blacsOrbSqr, eigvecs(:,:,iKS),&
            & denseDesc%blacsOrbSqr, globalSDotC, denseDesc%blacsOrbSqr)
        globalFrac(:,:) = real(conjg(eigvecs(:,:,iKS)) * globalSDotC)
        do iEig = 1, nOrb
          if (env%mpi%tGroupMaster) then
            call collector%getline_master(env%blacs%orbitalGrid, iEig, globalFrac, localFrac)
            call mpifx_send(env%mpi%interGroupComm, localFrac, env%mpi%interGroupComm%masterrank)
          else
            call collector%getline_slave(env%blacs%orbitalGrid, iEig, globalFrac)
          end if
        end do
      end do
    end if masterOrSlave

    if (env%mpi%tGlobalMaster) then
      call finishProjEigvecFiles(fd)
    end if

  end subroutine writeProjCplxEigvecsBlacs

#:else

  !> Write the projected complex eigenstates into text files.
  subroutine writeProjCplxEigvecsSerial(fileNames, fd, eigvals, neighlist, nNeighbor, cellVec,&
      & iCellVec, denseDesc, iPair, img2CentCell, over, kPoints, kWeights, parallelKS, eigvecs,&
      & work, iOrbRegion)

    !> list of region names
    type(ListCharLc), intent(inout) :: fileNames

    !> Fileid (file not yet opened) to use.
    integer, intent(in) :: fd(:)

    !> eigenvalues
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Neighbor list.
    type(TNeighborList), intent(in) :: neighlist

    !> Nr. of neighbors for SK-interaction.
    integer, intent(in) :: nNeighbor(:)

    !> Cell vectors of shifted cells.
    real(dp), intent(in) :: cellVec(:,:)

    !> Cell vector index of every atom.
    integer, intent(in) :: iCellVec(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Positions of interactions in the sparse matrices.
    integer, intent(in) :: iPair(:,:)

    !> Mapping of atoms into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Sparse overlap matrix.
    real(dp), intent(in) :: over(:)

    !> KPoints.
    real(dp), intent(in) :: kPoints(:,:)

    !> KPoints weights
    real(dp), intent(in) :: kWeights(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Eigen vectors
    complex(dp), intent(inout) :: eigvecs(:,:,:)

    !> Work array to unpack the overlap matrix
    complex(dp), intent(out) :: work(:,:)

    !> orbital number in each region
    type(listIntR1), intent(inout) :: iOrbRegion

    integer :: iKS, iS, iK, iEig, nOrb
    complex(dp), allocatable :: cVecTemp(:)

    nOrb = denseDesc%fullSize

    call prepareProjEigvecFiles(fd, fileNames)

    allocate(cVecTemp(size(eigvecs, dim=1)))
    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iS = parallelKS%localKS(2, iKS)
      call writeProjEigvecHeader(fd, iS, iK, kWeights(iK))
      call unpackHS(work, over, kPoints(:,iK), neighlist%iNeighbor, nNeighbor, iCellVec,&
          & cellVec, denseDesc%iAtomStart, iPair, img2CentCell)
      do iEig = 1, nOrb
        call hemv(cVecTemp, work, eigvecs(:,iEig,iKS))
        cVecTemp(:) = cVecTemp * conjg(eigvecs(:,iEig,iKS))
        call writeProjEigvecData(fd, iOrbRegion, eigvals(iEig, iK, iS), real(cVecTemp))
      end do
      call writeProjEigvecFooter(fd)
    end do

    call finishProjEigvecFiles(fd)

  end subroutine writeProjCplxEigvecsSerial

#:endif


#:if WITH_SCALAPACK

  !> Write the complex eigvectors into human readible output file (BLACS version).
  subroutine writeProjPauliEigvecsBlacs(env, denseDesc, fileNames, fd, iOrbRegion, eigvals,&
      & eigvecs, orb, parallelKS, kPoints, kWeights, over, neighborList, nNeighbor, iSparseStart,&
      & img2CentCell, iCellVec, cellVec)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> List of region file names
    type(ListCharLc), intent(inout) :: fileNames

    !> File descriptor for a not yet opened file for  each region
    integer, intent(in) :: fd(:)

    !> orbital number in each region
    type(listIntR1), intent(inout) :: iOrbRegion

    !> Eigenvalues
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Eigenvectors
    complex(dp), intent(in) :: eigvecs(:,:,:)

    !> Basis orbital information
    type(TOrbitals), intent(in) :: orb

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> K-points
    real(dp), intent(in) :: kPoints(:,:)

    !> Weights of the k-points
    real(dp), intent(in) :: kWeights(:)

    !> Sparse overlap
    real(dp), intent(in) :: over(:)

    !> Neighbors of each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Nr. of neighbors for each atom
    integer, intent(in) :: nNeighbor(:)

    !> Index array for sparse matrices
    integer, intent(in) :: iSparseStart(:,:)

    !> Mapping of atoms into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Cell vector index for each atom
    integer, intent(in) :: iCellVec(:)

    !> Cell vectors
    real(dp), intent(in) :: cellVec(:,:)

    type(linecomm) :: collector
    complex(dp), allocatable :: localSDotC(:), localEigvec(:)
    complex(dp), allocatable :: globalS(:,:), globalSDotC(:,:)
    real(dp), allocatable :: fracs(:,:)
    integer :: nOrb
    integer :: iKS, iK, iGroup, iEig

    nOrb = denseDesc%fullSize
    allocate(globalS(size(eigvecs, dim=1), size(eigvecs, dim=2)))
    allocate(globalSDotC(size(eigvecs, dim=1), size(eigvecs, dim=2)))
    if (env%mpi%tGroupMaster) then
      allocate(localEigvec(nOrb))
      allocate(localSDotC(nOrb))
      if (env%mpi%tGlobalMaster) then
        allocate(fracs(4, nOrb / 2))
      end if
    end if

    if (env%mpi%tGlobalMaster) then
      call prepareProjEigvecFiles(fd, fileNames)
    end if
    call collector%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")

    ! See comment about algorithm in routine write${NAME}$EigvecsBinBlacs

    masterOrSlave: if (env%mpi%tGlobalMaster) then
      ! Global master process
      do iKS = 1, parallelKS%maxGroupKS
        group: do iGroup = 0, env%mpi%nGroup - 1
          if (iKS > parallelKS%nGroupKS(iGroup)) then
            cycle group
          end if
          iK = parallelKS%groupKS(1, iKS, iGroup)
          if (iGroup == 0) then
            call unpackSPauliBlacs(env%blacs, over, kPoints(:,iK), neighborList%iNeighbor,&
                & nNeighbor, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
                & globalS)
            call pblasfx_phemm(globalS, denseDesc%blacsOrbSqr, eigvecs(:,:,iKS),&
                & denseDesc%blacsOrbSqr, globalSDotC, denseDesc%blacsOrbSqr)
          end if
          call writeProjEigvecHeader(fd, 1, iK, kWeights(iK))
          do iEig = 1, nOrb
            if (iGroup == 0) then
              call collector%getline_master(env%blacs%orbitalGrid, iEig, eigvecs(:,:,iKS),&
                  & localEigvec)
              call collector%getline_master(env%blacs%orbitalGrid, iEig, globalSDotC, localSDotC)
            else
              call mpifx_recv(env%mpi%interGroupComm, localEigvec, iGroup)
              call mpifx_recv(env%mpi%interGroupComm, localSDotC, iGroup)
            end if
            call getPauliFractions(localEigvec, localSDotC, fracs)
            call writeProjPauliEigvecData(fd, iOrbRegion, eigvals(iEig, iK, 1), fracs)
          end do
          call writeProjEigvecFooter(fd)
        end do group
      end do
    else
      ! All processes except the global master process
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        call unpackSPauliBlacs(env%blacs, over, kPoints(:,iK), neighborList%iNeighbor, nNeighbor,&
            & iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc, globalS)
        call pblasfx_phemm(globalS, denseDesc%blacsOrbSqr, eigvecs(:,:,iKS),&
            & denseDesc%blacsOrbSqr, globalSDotC, denseDesc%blacsOrbSqr)
        do iEig = 1, nOrb
          if (env%blacs%orbitalGrid%master) then
            call collector%getline_master(env%blacs%orbitalGrid, iEig, eigvecs(:,:,iKS),&
                & localEigvec)
            call collector%getline_master(env%blacs%orbitalGrid, iEig, globalSDotC, localSDotC)
            call mpifx_send(env%mpi%interGroupComm, localEigvec, env%mpi%interGroupComm%masterrank)
            call mpifx_send(env%mpi%interGroupComm, localSDotC, env%mpi%interGroupComm%masterrank)
          else
            call collector%getline_slave(env%blacs%orbitalGrid, iEig, eigvecs(:,:,iKS))
            call collector%getline_slave(env%blacs%orbitalGrid, iEig, globalSDotC)
          end if
        end do
      end do
    end if masterOrSlave

    if (env%mpi%tGlobalMaster) then
      call finishProjEigvecFiles(fd)
    end if

  end subroutine writeProjPauliEigvecsBlacs

#:else

  !> Write the projected complex eigenstates into text files.
  subroutine writeProjPauliEigvecsSerial(fileNames, fd, eigvals, neighlist, nNeighbor, cellVec,&
      & iCellVec, denseDesc, iPair, img2CentCell, over, kPoints, kWeights, parallelKS, eigvecs,&
      & work, iOrbRegion)

    !> list of region names
    type(ListCharLc), intent(inout) :: fileNames

    !> Fileid (file not yet opened) to use.
    integer, intent(in) :: fd(:)

    !> eigenvalues
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Neighbor list.
    type(TNeighborList), intent(in) :: neighlist

    !> Nr. of neighbors for SK-interaction.
    integer, intent(in) :: nNeighbor(:)

    !> Cell vectors of shifted cells.
    real(dp), intent(in) :: cellVec(:,:)

    !> Cell vector index of every atom.
    integer, intent(in) :: iCellVec(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Positions of interactions in the sparse matrices.
    integer, intent(in) :: iPair(:,:)

    !> Mapping of atoms into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Sparse overlap matrix.
    real(dp), intent(in) :: over(:)

    !> KPoints.
    real(dp), intent(in) :: kPoints(:,:)

    !> KPoints weights
    real(dp), intent(in) :: kWeights(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Eigenvectors
    complex(dp), intent(inout) :: eigvecs(:,:,:)

    !> Work array to unpack S.
    complex(dp), intent(out) :: work(:,:)

    !> orbital number in each region
    type(listIntR1), intent(inout) :: iOrbRegion

    complex(dp), allocatable :: cVecTemp(:)
    real(dp), allocatable :: fracs(:,:)
    integer :: nOrb
    integer :: iKS, iK, iEig

    nOrb = denseDesc%fullSize

    call prepareProjEigvecFiles(fd, fileNames)

    allocate(cVecTemp(size(eigvecs, dim=1)))
    allocate(fracs(4, nOrb / 2))

    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      call writeProjEigvecHeader(fd, 1, iK, kWeights(iK))
      call unpackSPauli(over, kPoints(:,iK), neighlist%iNeighbor, nNeighbor,&
          & denseDesc%iAtomStart, iPair, img2CentCell, iCellVec, cellVec, work)
      do iEig = 1, nOrb
        call hemv(cVecTemp, work, eigvecs(:,iEig,iKS))
        call getPauliFractions(eigvecs(:,iEig,iKS), cVecTemp, fracs)
        call writeProjPauliEigvecData(fd, iOrbRegion, eigvals(iEig, iK, 1), fracs)
      end do
      call writeProjEigvecFooter(fd)
    end do

    call finishProjEigvecFiles(fd)

  end subroutine writeProjPauliEigvecsSerial

#:endif


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
      if (size(excitedDerivs) > 0) then
        call writeTagged(fd, tag_excForce, -excitedDerivs)
      end if
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

    @:ASSERT(tPeriodic .eqv. tStress)

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
      & tAtomicEnergy, tDispersion, tEField, tPeriodic, nSpin, tSpinOrbit, tScc,&
      & invLatVec, kPoints)

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

    !> Is this a molecular dynamics run
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

    !> Reciprocal lattice vectors if periodic
    real(dp), intent(in) :: invLatVec(:,:)

    !> K-points if periodic
    real(dp), intent(in) :: kPoints(:,:)

    real(dp), allocatable :: qInputUpDown(:,:,:), qOutputUpDown(:,:,:), qBlockOutUpDown(:,:,:,:)
    real(dp) :: angularMomentum(3)
    integer :: ang
    integer :: nAtom, nLevel, nKPoint, nSpinHams, nMovedAtom
    integer :: iAt, iSpin, iEgy, iK, iSp, iSh, iOrb, kk
    logical :: tSpin

    character(*), parameter :: formatEigen = "(F14.8)"
    character(*), parameter :: formatFilling = "(F12.5)"
    character(lc) :: strTmp

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
          write(fd, "(A, I0, A, I0)") "Geometry optimization step: ", iGeoStep,&
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

    if (tPeriodic .and. tLatOpt) then
      do iK = 1, nKPoint
        if (iK == 1) then
          write(strTmp, "(A,':')") "K-points in absolute space"
        else
          write(strTmp, "(A)") ""
        end if
        write(fd, "(A,T28,I6,':',3F10.6)") trim(strTmp), iK, matmul(invLatVec,kPoints(:,iK))
      end do
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
      write(fd, "(A, F14.8)") " Total charge: ", sum(q0(:, :, 1) - qOutput(:, :, 1))
      write(fd, "(/,A)") " Atomic gross charges (e)"
      write(fd, "(A5, 1X, A16)")" Atom", " Charge"
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
      write(fd, format2U) "Electronic Temperature", tempElec, 'au', tempElec * Hartree__eV,&
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
    real(dp), intent(in), allocatable :: dipoleMoment(:)

    if (tEfield) then
      write(fd, format1U1e) 'External E field', absEField, 'au', absEField * au__V_m, 'V/m'
    end if

    if (allocated(dipoleMoment)) then
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
    real(dp), intent(inout), allocatable :: dipoleMoment(:)

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
      write(fd, format2U) "Excitation Energy", energy%Eexcited, "H",&
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
    if (allocated(dipoleMoment)) then
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
  subroutine writeCharges(fCharges, fdCharges, tWriteBinary, orb, qInput, qBlockIn, qiBlockIn)

    !> File name for charges to be written to
    character(*), intent(in) :: fCharges

    !> File descriptor for charge output
    integer, intent(in) :: fdCharges

    !> Charges should be output in binary (T) or ascii (F)
    logical, intent(in) :: tWriteBinary

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> input charges
    real(dp), intent(in) :: qInput(:,:,:)

    !> Block populations if present
    real(dp), intent(in), allocatable :: qBlockIn(:,:,:,:)

    !> Imaginary part of block populations if present
    real(dp), intent(in), allocatable :: qiBlockIn(:,:,:,:)

    if (allocated(qBlockIn)) then
      if (allocated(qiBlockIn)) then
        call writeQToFile(qInput, fCharges, fdCharges, tWriteBinary, orb, qBlockIn, qiBlockIn)
      else
        call writeQToFile(qInput, fCharges, fdCharges, tWriteBinary, orb, qBlockIn)
      end if
    else
      call writeQToFile(qInput, fCharges, fdCharges, tWriteBinary, orb)
    end if
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

    !> Dense matrix indexing for atomic blocks
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
    call writeHS(tWriteHS, tWriteRealHS, tRealHS, hamUpDown, over, neighborList%iNeighbor,&
        & nNeighbor, iAtomStart, iPair, img2CentCell, kPoint, iCellVec, cellVec, iHam)
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
    real(dp), intent(in), allocatable :: iHam(:,:)

    integer :: iS, nSpin

    nSpin = size(ham, dim=2)

    if (tWriteRealHS) then
      do iS = 1, nSpin
        call writeSparse("hamreal" // i2c(iS) // ".dat", ham(:,iS), iNeighbor, &
            &nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, cellVec)
        if (allocated(iHam)) then
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
    real(dp), intent(in), allocatable :: qOutput(:,:,:)

    !> atomic velocities
    real(dp), intent(in), allocatable :: velocities(:,:)

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
      write(comment, "(A, I0, A, I0)") '** Geometry step: ', iGeoStep, ', Lattice step: ',&
          & iLatGeoStep
    elseif (tMD) then
      write(comment, "(A, I0)") 'MD iter: ', iGeoStep
    else
      write(comment,"(A, I0)") 'Geometry Step: ', iGeoStep
    end if

    if (tPrintMulliken) then
      ! For non-colinear spin without velocities write magnetisation into the velocity field
      if (nSpin == 4 .and. .not. allocated(velocities)) then
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
      else if (allocated(velocities)) then
        call writeXYZFormat(fname, pCoord0Out, species0, speciesName,&
            & charges=sum(qOutput(:,:,1),dim=1), velocities=velocities, comment=comment,&
            & append=tAppendGeo)
      else
        call writeXYZFormat(fname, pCoord0Out, species0, speciesName,&
            & charges=sum(qOutput(:,:,1),dim=1), comment=comment, append=tAppendGeo)
      end if
    else if (allocated(velocities)) then
      call writeXYZFormat(fname, pCoord0Out, species0, speciesName, velocities=velocities,&
          & comment=comment, append=tAppendGeo)
    else
      call writeXYZFormat(fname, pCoord0Out, species0, speciesName, comment=comment,&
          & append=tAppendGeo)
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


#:if WITH_SOCKETS

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

    @:ASSERT(env%tGlobalMaster .eqv. allocated(socket))

    if (env%tGlobalMaster) then
      call socket%receive(coord0, tmpLatVecs, tStopDriver)
    end if
    tCoordsChanged = .true.
    if (tPeriodic .and. .not. tStopDriver) then
      latVecs(:,:) = tmpLatVecs
    end if
    tLatticeChanged = tPeriodic
  #:if WITH_MPI
    ! update all nodes with the received information
    call mpifx_bcast(env%mpi%globalComm, coord0)
    call mpifx_bcast(env%mpi%globalComm, latVecs)
    call mpifx_bcast(env%mpi%globalComm, tCoordsChanged)
    call mpifx_bcast(env%mpi%globalComm, tLatticeChanged)
    call mpifx_bcast(env%mpi%globalComm, tStopDriver)
  #:endif

  end subroutine receiveGeometryFromSocket

#:endif


  !> Prepares binary eigenvector file for writing.
  subroutine prepareEigvecFileBin(fd, runId, fileName)

    !> File descriptor to use
    integer, intent(in) :: fd

    !> Run id to write into the file header
    integer, intent(in) :: runId

    !> Name of the file
    character(*), intent(in), optional :: fileName

    character(lc) :: tmpStr

    if (present(fileName)) then
      write(tmpStr, "(A,A)") trim(fileName), ".bin"
      open(fd, file=tmpStr, action="write", status="replace", form="unformatted")
    else
      open(fd, file=eigvecBin, action="write", status="replace", form="unformatted")
    end if
    write(fd) runId

  end subroutine prepareEigvecFileBin


  !> Prepares text eigenvector file for writing.
  subroutine prepareEigvecFileTxt(fd, t2Component, fileName)

    !> File descriptor to use
    integer, intent(in) :: fd

    !> Whether eigenvectors present 2-component Pauli vectors
    logical, intent(in) :: t2Component

    !> Name of the file
    character(*), intent(in), optional :: fileName

    character(lc) :: tmpStr

    if (present(fileName)) then
      write(tmpStr, "(A,A)") trim(fileName), ".out"
      open(fd, file=tmpStr, action="write", status="replace", position="rewind")
    else
      open(fd, file=eigvecOut, action="write", status="replace", position="rewind")
    end if
    write(fd, "(A/)") "Coefficients and Mulliken populations of the atomic orbitals"
    if (t2Component) then
      write(fd,"(A/)")"   Atom   Orb  up spin coefficients       &
          & down spin coefficients         charge      x           y           z"
    end if

  end subroutine prepareEigvecFileTxt


  !> Writes a single real eigenvector into a file
  subroutine writeSingleRealEigvecTxt(fd, eigvec, fracs, iS, iEigvec, orb, species, speciesName,&
      & nAtom)

    !> File descriptor of open file
    integer, intent(in) :: fd

    !> Eigenvector to write
    real(dp), intent(in) :: eigvec(:)

    !> Fraction of each component in the eigenvector (c.S.c)
    real(dp), intent(in) :: fracs(:)

    !> Spin index of the eigenvector
    integer, intent(in) :: iS

    !> Index of the eigenvector
    integer, intent(in) :: iEigvec

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species for each atom
    integer, intent(in) :: species(:)

    !> Name of each species
    character(*), intent(in) :: speciesName(:)

    !> Number of atoms
    integer, intent(in) :: nAtom

    character(lc) :: tmpStr
    integer :: ind, ang
    integer :: iAt, iSp, iSh, iOrb

    write(fd, "('Eigenvector:',I4,4X,'(',A,')'/)") iEigvec, trim(spinName(iS))
    ind = 0
    do iAt = 1, nAtom
      iSp = species(iAt)
      do iSh = 1, orb%nShell(iSp)
        ang = orb%angShell(iSh, iSp)
        if (iSh == 1) then
          write(tmpStr, "(I5,1X,A2,2X,A1)") iAt, speciesName(iSp), orbitalNames(ang + 1)
        else
          write(tmpStr, "(10X,A1)") orbitalNames(ang + 1)
        end if
        do iOrb = 1, 2 * ang + 1
          ind = ind + 1
          write(fd, "(A,I1,T15,F12.6,3X,F12.6)") trim(tmpStr), iOrb, eigvec(ind), fracs(ind)
        end do
      end do
      write(fd,*)
    end do

  end subroutine writeSingleRealEigvecTxt


  !> Writes a single complex eigenvector into a file
  subroutine writeSingleCplxEigvecTxt(fd, eigvec, fracs, iS, iK, iEigvec, orb, species,&
      & speciesName, nAtom)

    !> File descriptor of open file
    integer, intent(in) :: fd

    !> Eigenvector to write
    complex(dp), intent(in) :: eigvec(:)

    !> Fraction of each basis function in the eigenvector (c.S.c)
    real(dp), intent(in) :: fracs(:)

    !> Spin index of the eigenvector
    integer, intent(in) :: iS

    !> K-point index of the eigenvector
    integer, intent(in) :: iK

    !> Index of the eigenvector
    integer, intent(in) :: iEigvec

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species for each atom
    integer, intent(in) :: species(:)

    !> Name of each species
    character(*), intent(in) :: speciesName(:)

    !> Number of atoms
    integer, intent(in) :: nAtom

    character(lc) :: tmpStr
    integer :: ind, ang
    integer :: iAt, iSp, iSh, iOrb

    write(fd, "(A,I4,4X,A,I4,4X,'(',A,')'/)") "K-point: ", iK, "Eigenvector: ", iEigvec,&
        & trim(spinName(iS))
    ind = 0
    do iAt = 1, nAtom
      iSp = species(iAt)
      do iSh = 1, orb%nShell(iSp)
        ang = orb%angShell(iSh, iSp)
        if (iSh == 1) then
          write(tmpStr, "(I5,1X,A2,2X,A1)") iAt, speciesName(iSp), orbitalNames(ang + 1)
        else
          write(tmpStr, "(10X,A1)") orbitalNames(ang + 1)
        end if
        do iOrb = 1, 2 * ang + 1
          ind = ind + 1
          write(fd, "(A,I1,T15,'(',F12.6,',',F12.6,')',3X,F12.6)") trim(tmpStr), iOrb,&
              & real(eigvec(ind)), aimag(eigvec(ind)), fracs(ind)
        end do
      end do
      write(fd,*)
    end do

  end subroutine writeSingleCplxEigvecTxt


  !> Writes a single Pauli two-component eigenvector into a file
  subroutine writeSinglePauliEigvecTxt(fd, eigvec, fracs, iK, iEigvec, orb, species, speciesName,&
      & nAtom, nOrb)

    !> File descriptor of open file
    integer, intent(in) :: fd

    !> Eigenvector to write
    complex(dp), intent(in) :: eigvec(:)

    !> Fraction of each orbital in the eigenvector, decomposed into 4 components.
    real(dp), intent(in) :: fracs(:,:)

    !> K-point index of the eigenvector
    integer, intent(in) :: iK

    !> Index of the eigenvector
    integer, intent(in) :: iEigvec

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species for each atom
    integer, intent(in) :: species(:)

    !> Name of each species
    character(*), intent(in) :: speciesName(:)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Number of orbitals
    integer, intent(in) :: nOrb

    character(lc) :: tmpStr
    integer :: ind, ang
    integer :: iAt, iSp, iSh, iOrb

    write(fd, "(A,I4,4X,A,I4)") "K-point: ", ik, "Eigenvector: ", iEigvec
    ind = 0
    do iAt = 1, nAtom
      iSp = species(iAt)
      do iSh = 1, orb%nShell(iSp)
        ang = orb%angShell(iSh,iSp)
        if (iSh == 1) then
          write(tmpStr, "(I5,1X,A2,2X,A1)") iAt, speciesName(iSp), orbitalNames(ang + 1)
        else
          write(tmpStr, "(10X,A1)") orbitalNames(ang + 1)
        end if
        do iOrb = 1, 2 * ang + 1
          ind = ind + 1
          write(fd, "(A,I1,T15,'(',F12.6,',',F12.6,')','(',F12.6,',',F12.6,')',3X,4F12.6)")&
              & trim(tmpStr), iOrb, real(eigvec(ind)), aimag(eigvec(ind)),&
              & real(eigvec(ind + nOrb)), aimag(eigvec(ind + nOrb)),&
              & fracs(:, ind)
        end do
      end do
      write(fd,*)
    end do

  end subroutine writeSinglePauliEigvecTxt


  !> Write projected real eigenvector data to disc
  subroutine writeProjEigvecData(fd, iOrbRegion, eigval, fracs)

    !> File descriptor for each region
    integer, intent(in) :: fd(:)

    !> List of orbital for each region
    type(listIntR1), intent(inout) :: iOrbRegion

    !> Eigenvalue for current eigenvector
    real(dp), intent(in) :: eigval

    !> Fraction of each orbital in the current eigenvector (c.S.c)
    real(dp), intent(in) :: fracs(:)

    integer, allocatable :: iOrbs(:)
    integer :: valShape(1)
    integer :: iReg, dummy

    do iReg = 1, size(fd)
      call elemShape(iOrbRegion, valshape, iReg)
      allocate(iOrbs(valshape(1)))
      call intoArray(iOrbRegion, iOrbs, dummy, iReg)
      write(fd(iReg), "(f13.6,f10.6)") Hartree__eV * eigval, sum(fracs(iOrbs))
      deallocate(iOrbs)
    end do

  end subroutine writeProjEigvecData


  !> Write projected real eigenvector data to disc (complex)
  subroutine writeProjPauliEigvecData(fd, iOrbRegion, eigval, fracs)

    !> File descriptor for each region
    integer, intent(in) :: fd(:)

    !> List of orbital for each region
    type(listIntR1), intent(inout) :: iOrbRegion

    !> Eigenvalue for current eigenvector
    real(dp), intent(in) :: eigval

    !> Fraction of each orbital in the eigenvector for the four Pauli components
    real(dp), intent(in) :: fracs(:,:)

    integer, allocatable :: iOrbs(:)
    integer :: valShape(1)
    integer :: iReg, dummy

    do iReg = 1, size(fd)
      call elemShape(iOrbRegion, valshape, iReg)
      allocate(iOrbs(valshape(1)))
      call intoArray(iOrbRegion, iOrbs, dummy, iReg)
      write(fd(iReg), "(f13.6,4f10.6)") Hartree__eV * eigval, sum(fracs(:,iOrbs), dim=2)
      deallocate(iOrbs)
    end do

  end subroutine writeProjPauliEigvecData


  !> Writes header for projected eigenvectors
  subroutine writeProjEigvecHeader(fd, iS, iK, kWeight)

    !> File descriptor for each region
    integer, intent(in) :: fd(:)

    !> Index fo current spin
    integer, intent(in) :: iS

    !> Index of current k-point
    integer, intent(in), optional :: iK

    !> Weight of current k-point
    real(dp), intent(in), optional :: kWeight

    integer :: iReg

    @:ASSERT(present(iK) .eqv. present(kWeight))

    do iReg = 1, size(fd)
      if (present(iK)) then
        write(fd(iReg), "(2(A,1X,I0,1X),A,1X,F12.8)") 'KPT', iK, 'SPIN', iS, 'KWEIGHT', kWeight
      else
        write(fd(iReg), "(A,1X,I0)") 'SPIN', iS
      end if
    end do

  end subroutine writeProjEigvecHeader


  !> Writes footer for projected eigenvectors
  subroutine writeProjEigvecFooter(fd)

    !> File descriptor for each region
    integer, intent(in) :: fd(:)

    integer :: iReg

    do iReg = 1, size(fd)
      write(fd(iReg), "(A)") ""
    end do

  end subroutine writeProjEigvecFooter


  !> Returns the fraction of each orbital in a 2-component Pauli-vector.
  subroutine getPauliFractions(eigvec, overDotEigvec, fracs)

    !> Pauli eigenvector
    complex(dp), intent(in) :: eigvec(:)

    !> Overlap times eigenvector
    complex(dp), intent(in) :: overDotEigvec(:)

    !> Fractions along the 4 (total, x, y and z) component. Shape (4, size(eigvec) / 2)
    real(dp), intent(out) :: fracs(:,:)

    integer :: nOrb
    integer :: iOrb

    nOrb = size(eigvec) / 2
    do iOrb = 1, nOrb
      fracs(1, iOrb) =  real(conjg(eigvec(iOrb)) * overDotEigvec(iOrb)&
          & + conjg(eigvec(iOrb + nOrb)) * overDotEigvec(iOrb + nOrb))
      fracs(2, iOrb) = real(conjg(eigvec(iOrb + nOrb)) * overDotEigvec(iOrb)&
          & + conjg(eigvec(iOrb)) * overDotEigvec(iOrb + nOrb))
      fracs(3, iOrb) = aimag(conjg(eigvec(iOrb)) * overDotEigvec(iOrb + nOrb)&
          & - conjg(eigvec(iOrb + nOrb)) * overDotEigvec(iOrb))
      fracs(4, iOrb) = real(conjg(eigvec(iOrb)) * overDotEigvec(iOrb)&
          & - conjg(eigvec(iOrb + nOrb)) * overDotEigvec(iOrb + nOrb))
    end do

  end subroutine getPauliFractions


  !> Prepare projected eigenvector file for each region
  subroutine prepareProjEigvecFiles(fd, fileNames)

    !> File descriptor for a not yet opened file for each region
    integer, intent(in) :: fd(:)

    !> List of region file names
    type(ListCharLc), intent(inout) :: fileNames

    integer :: iReg
    character(lc) :: tmpStr

    do iReg = 1, size(fd)
      call get(fileNames, tmpStr, iReg)
      open(fd(iReg), file=tmpStr, action="write", status="replace", form="formatted")
    end do

  end subroutine prepareProjEigvecFiles


  !> Finish projected eigenvector file for each region
  subroutine finishProjEigvecFiles(fd)

    !> File descriptor for opened file for each region
    integer, intent(in) :: fd(:)

    integer :: iReg

    do iReg = 1, size(fd)
      close(fd(iReg))
    end do

  end subroutine finishProjEigvecFiles


end module mainio
