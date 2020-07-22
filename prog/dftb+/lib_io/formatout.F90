!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines for formatted output of data
module dftbp_formatout
  use dftbp_globalenv
  use dftbp_environment
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_message
  use dftbp_constants
  use dftbp_lapackroutines, only: matinv
  use dftbp_sparse2dense
  implicit none
  private

  public :: clearFile, writeGenFormat, writeXYZFormat
  public :: printDFTBHeader
  public :: writeSparseAsSquare, writeSparse


  !> Clears contents of a file
  interface clearFile
    module procedure clearFile_fname
  end interface clearFile


  !> Writes geometry information in gen format to a file
  interface writeGenFormat
    module procedure writeGenFormat_fname
    module procedure writeGenFormat_fid
  end interface writeGenFormat


  !> Writes geometry information in xyz format to a file
  interface writeXYZFormat
    module procedure writeXYZFormat_fname
    module procedure writeXYZFormat_fid
  end interface writeXYZFormat


  !> Writes DFTB+ type sparse matrix in square form to disc
  interface writeSparseAsSquare
    module procedure writeSparseAsSquare_real
    module procedure writeSparseAsSquare_cplx
  end interface writeSparseAsSquare

contains


  !> Clears contents of file
  subroutine clearFile_fname(fileName)

    !> name of the file which should be cleared
    character(len=*), intent(in) :: fileName

    integer :: fd

    open(newunit=fd, file=fileName, status="replace", position="rewind")
    close(fd)

  end subroutine clearFile_fname


  !> A wrapper around writeGenFormat_fid to open a file first.
  subroutine writeGenFormat_fname(fileName, coord, species, speciesName, latVec, origin,&
      & tFracCoord, append)

    !> File name of the file which should be created
    character(len=*), intent(in) :: fileName

    !> Coordinates in atomic units
    real(dp), intent(in) :: coord(:,:)

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Name of the different species
    character(mc), intent(in) :: speciesName(:)

    !> Lattice vectors
    real(dp), intent(in), optional :: latVec(:,:)

    !> origin of coordinate system, if not molecular
    real(dp), intent(in), optional :: origin(:)

    !> Print out fractional coordinates?
    logical, intent(in), optional :: tFracCoord

    !> Whether geometry should be appended (default: it is overwritten)
    logical, intent(in), optional :: append

    integer :: fd

    logical :: append0

    if (present(append)) then
      append0 = append
    else
      append0 = .false.
    end if

    @:ASSERT((.not.(present(tFracCoord).neqv.present(latVec))) .or.(present(latVec)))
    @:ASSERT(present(latVec) .eqv. present(origin))

    if (append0) then
      open(newunit=fd, file=fileName, form="formatted", action="write", status="old",&
          & position="append")
    else
      open(newunit=fd, file=fileName, form="formatted", action="write", status="replace")
    end if
    call writeGenFormat_fid(fd, coord, species, speciesName, latVec, origin, tFracCoord)
    close(fd)

  end subroutine writeGenFormat_fname


  !> Writes coordinates in the famous GEN format to a file
  subroutine writeGenFormat_fid(fd, coord, species, speciesName, latVec, origin, tFracCoord)

    !> File id of an open file where output should be written
    integer, intent(in) :: fd

    !> Coordinates in atomic units
    real(dp),          intent(in) :: coord(:,:)

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Name of the different species
    character(mc), intent(in) :: speciesName(:)

    !> Lattice vectors
    real(dp), intent(in), optional :: latVec(:,:)

    !> origin of coordinate system, if not molecular
    real(dp), intent(in), optional :: origin(:)

    !> Print out fractional coordinates?
    logical, intent(in), optional :: tFracCoord

    integer :: nAtom, nSpecies
    character(mc) :: formatCoordinates
    integer :: ii
    logical :: tFractional, tHelical, tPeriodic
    real(dp) :: invLatVec(3,3)

    ! Format first line of a gen file
    character(len=*), parameter :: formatHead1 = '(1X,I0,1X,A2)'

    ! Format a vector
    character(len=*), parameter :: formatVec = '(3E20.10)'

    nAtom = size(coord, dim=2)
    nSpecies = maxval(species)

    @:ASSERT(size(coord, dim=1) == 3)
    @:ASSERT(size(species) == nAtom)
    @:ASSERT(size(speciesName) == nSpecies)
  #:block DEBUG_CODE
    if (present(latVec)) then
      @:ASSERT(all(shape(latVec) == [3,3]) .or. all(shape(latVec) == [2,1])&
          & .or. all(shape(latVec) == [3,1]))
      @:ASSERT(.not.(present(tFracCoord) .and.&
          & (all(shape(latVec) == [2, 1]) .or. all(shape(latVec) == [2, 1]) ) ) )
    end if
  #:endblock DEBUG_CODE
    @:ASSERT((.not.(present(tFracCoord).neqv.present(latVec))) .or.(present(latVec)))
    @:ASSERT(present(latVec) .eqv. present(origin))

    tPeriodic = .false.
    tFractional = .false.
    tHelical = .false.

    if (present(latVec)) then
      if (present(tFracCoord) ) then
        tFractional = tFracCoord
      end if
      if (all(shape(latVec) == [2, 1]) .or. all(shape(latVec) == [3,1])) then
        tHelical = .true.
      else
        tPeriodic = .true.
      end if
      if (tFractional) then
        write(fd, formatHead1) nAtom, "F"
      else if (tHelical) then
        write(fd, formatHead1) nAtom, "H"
      else if (tPeriodic) then
        write(fd, formatHead1) nAtom, "S"
      else
        call error("Unknown boundary conditions")
      end if
    else
      write(fd, formatHead1) nAtom, "C"
    end if

    do ii = 1, nSpecies
      write(fd, '(1X,A)', advance='NO')trim(speciesName(ii))
    end do
    write(fd,'(A)')

    write(formatCoordinates, '("(I5,2X,I",I0,",3E20.10)")') floor(log10(real(nSpecies)))+1
    if (tFractional) then
      invLatVec(:,:) = latVec
      call matinv(invLatVec)
      do ii = 1, nAtom
        write(fd, formatCoordinates) ii, species(ii), matmul(invLatVec,coord(:, ii) + origin)
      end do
    else if (tPeriodic .or. tHelical) then
      do ii = 1, nAtom
        write(fd, formatCoordinates) ii, species(ii), (coord(:, ii) + origin) * Bohr__AA
      end do
    else
      do ii = 1, nAtom
        write(fd, formatCoordinates) ii, species(ii), coord(:, ii) * Bohr__AA
      end do
    end if
    if (present(latVec)) then
      if (tHelical) then
        write(fd, formatVec) origin * Bohr__AA
        if (size(latvec,dim=1)==2) then
          write(fd, "(E20.10,F14.8)") latVec(1, 1) * Bohr__AA, latVec(2, 1) * 180.0_dp/pi
        else
          write(fd, "(E20.10,F14.8,1X,I0)") latVec(1, 1) * Bohr__AA, latVec(2, 1) * 180.0_dp/pi, &
              & nint(latVec(3,1))
        end if
      else if (tPeriodic) then
        if (tFractional) then
          write(fd, formatVec) matmul(invLatVec, origin)
        else
          write(fd, formatVec) origin * Bohr__AA
        end if
        do ii = 1, 3
          write(fd, formatVec) latVec(:, ii) * Bohr__AA
        end do
      else
        call error("Unknown boundary conditions")
      end if
    end if

  end subroutine writeGenFormat_fid


  !> Writes coordinates in the XYZ format
  subroutine writeXYZFormat_fname(fileName, coord, species, speciesName, charges, velocities,&
      & comment, append)

    !> File name of a file to be created
    character(len=*), intent(in) :: fileName

    !> Coordinates in atomic units
    real(dp), intent(in) :: coord(:,:)

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Name of the different species
    character(mc), intent(in) :: speciesName(:)

    !> Optional vector with charges for each atom.
    real(dp), intent(in), optional :: charges(:)

    !> Optional array of velocity vectors for each atom.
    real(dp), intent(in), optional :: velocities(:,:)

    !> Optional comment for line 2 of the file
    character(len=*), intent(in), optional :: comment

    !> Whether geometry should be appended (default: it is overwritten)
    logical, intent(in), optional :: append

    integer :: fd
    logical :: append0

    if (present(append)) then
      append0 = append
    else
      append0 = .false.
    end if

    if (append0) then
      open(newunit=fd, file=fileName, action="write", form="formatted", status="old",&
          & position="append")
    else
      open(newunit=fd, file=fileName, action="write", form="formatted", status="replace")
    end if
    call writeXYZFormat(fd, coord, species, speciesName, charges, velocities, comment)
    close(fd)

  end subroutine writeXYZFormat_fname


  !> Writes coordinates in the XYZ format with additional charges and vectors
  subroutine writeXYZFormat_fid(fd, coords, species, speciesNames, charges, velocities, comment)

    !> File id of an open file where output should be written
    integer, intent(in) :: fd

    !> Coordinates in atomic units
    real(dp), intent(in) :: coords(:,:)

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Name of the different species
    character(mc), intent(in) :: speciesNames(:)

    !> Optional vector with charges for each atom.
    real(dp), intent(in), optional :: charges(:)

    !> Optional array of velocity vectors for each atom.
    real(dp), intent(in), optional :: velocities(:,:)

    !> Optional comment for line 2 of the file
    character(len=*), intent(in), optional :: comment

    integer :: nAtom, nSpecies, ii

200 format(I5)
201 format(A5,3F16.8)
202 format(A5,6F16.8)
203 format(A5,4F16.8)
204 format(A5,7F16.8)

    nAtom = size(coords, dim=2)
    nSpecies = maxval(species)
    @:ASSERT(size(coords, dim=1) == 3)
    @:ASSERT(size(species) == nAtom)
    @:ASSERT(size(speciesNames) == nSpecies)
  #:block DEBUG_CODE
    if (present(charges)) then
      @:ASSERT(size(charges) == nAtom)
    end if
    if (present(velocities)) then
      @:ASSERT(all(shape(velocities) == (/ 3, nAtom /)))
    end if
  #:endblock DEBUG_CODE

    write(fd, 200) nAtom
    if (present(comment)) then
      write(fd, "(A)") trim(comment)
    elseif (present(velocities)) then
      write(fd, *) "Velocity in AA/ps"
    else
      write(fd, *) ""
    end if

    if (present(charges) .and. present(velocities)) then
      write(fd, 204) (trim(speciesNames(species(ii))), coords(:, ii) * Bohr__AA,&
          & charges(ii), velocities(:,ii) * Bohr__AA / au__fs * 1000.0_dp, ii = 1, nAtom)
    elseif (present(charges) .and. .not. present(velocities)) then
      write(fd, 203) (trim(speciesNames(species(ii))), coords(:, ii) * Bohr__AA,&
          & charges(ii), ii = 1, nAtom)
    elseif (.not. present(charges) .and. present(velocities)) then
      write(fd, 202) (trim(speciesNames(species(ii))), coords(:, ii) * Bohr__AA,&
          & velocities(:,ii) * Bohr__AA / au__fs * 1000.0_dp, ii = 1, nAtom)
    else
      write(fd, 201) (trim(speciesNames(species(ii))),&
          & coords(:, ii) * Bohr__AA, ii = 1, nAtom)
    end if

  end subroutine writeXYZFormat_fid


  !> Writes the greeting message of dftb+ on stdout
  subroutine printDFTBHeader(release, year)

    !> release version of the code
    character(len=*), intent(in) :: release

    !> release year
    integer, intent(in) :: year

    character, parameter :: verticalBar = '|'
    character, parameter :: horizontalBar = '='
    integer, parameter :: headerWidth = 80

    write(stdOut, '(2A,/,A)') verticalBar, repeat(horizontalBar, headerWidth - 1), verticalBar
    write(stdOut, '(3A)') verticalBar, '  DFTB+ ', trim(release)
    write(stdOut, '(A)') verticalBar
    write(stdOut, '(2A,I0,A)') verticalBar, '  Copyright (C) 2006 - ', year,&
        & '  DFTB+ developers group'
    write(stdOut, '(A,/,2A,/,A)') verticalBar, verticalBar, repeat(horizontalBar, headerWidth - 1),&
        & verticalBar
    write(stdOut, '(2A)') verticalBar,&
        & '  When publishing results obtained with DFTB+, please cite the following', verticalBar,&
        & '  reference:'
    write(stdOut, '(A)') verticalBar
    write(stdOut, '(2A)') verticalBar,&
        & '  * DFTB+, a software package for efficient approximate density functional',&
        & verticalBar,&
        & '    theory based atomistic simulations, J. Chem. Phys. 152, 124101 (2020).',&
        & verticalBar, '    [doi: 10.1063/1.5143190]'
    write(stdOut, '(A)') verticalBar
    write(stdOut, '(2A,2(/,2A))') verticalBar,&
        & '  You should also cite additional publications crediting the parametrization',&
        & verticalBar,&
        & '  data you use. Please consult the documentation of the SK-files for the', verticalBar,&
        & '  references.'
    write(stdOut, '(A,/,2A,/)') verticalBar, verticalBar, repeat(horizontalBar, headerWidth - 1)

  end subroutine printDFTBHeader


  !> Converts a sparse matrix to its square form and write it to a file.
  subroutine writeSparseAsSquare_real(env, fname, sparse, iNeighbour, nNeighbourSK, iAtomStart,&
      & iPair, img2CentCell)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Name of the file to write the matrix to.
    character(len=*), intent(in) :: fname

    !> Sparse matrix.
    real(dp), intent(in) :: sparse(:)

    !> Neighbour list index.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours.
    integer, intent(in) :: nNeighbourSK(:)

    !> Offset array in the square matrix.
    integer, intent(in) :: iAtomStart(:)

    !> Offset array in the sparse matrix
    integer, intent(in) :: iPair(0:,:)

    !> Pair indexing array.
    integer, intent(in) :: img2CentCell(:)

    !> Mapping of the atoms to the central cell.
    real(dp), allocatable :: square(:,:)
    character(mc) :: strForm
    integer :: fd, nOrb

    if (withMpi) then
      call error("Writing of HS not working with MPI yet")
    end if

    nOrb = iAtomStart(size(nNeighbourSK) + 1) - 1

    allocate(square(nOrb, nOrb))
    open(newunit=fd, file=fname, form="formatted", status="replace")
    write(fd, "(A1,A10,A10,A10,A10)") "#", "REAL", "NALLORB", "NKPOINT"
    write(fd, "(1X,L10,I10,I10,I10)") .true., nOrb, 1

    write (strForm, "(A,I0,A)") "(", nOrb, "ES24.15)"
    call unpackHS(square, sparse, iNeighbour, nNeighbourSK, iAtomStart, iPair, img2CentCell)
    call blockSymmetrizeHS(square, iAtomStart)
    write(fd, "(A1,A10,A10)") "#", "IKPOINT"
    write(fd, "(1X,I10,I10)") 1
    write(fd, "(A1,A)") "#", " MATRIX"
    write(fd, strForm) square
    close(fd)

  end subroutine writeSparseAsSquare_real


  !> Converts a sparse matrix to its square form and write it to a file.
  subroutine writeSparseAsSquare_cplx(env, fname, sparse, kPoints, iNeighbour, nNeighbourSK,&
      & iAtomStart, iPair, img2CentCell, iCellVec, cellVec)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Name of the file to write the matrix into.
    character(len=*), intent(in) :: fname

    !> Sparse matrix.
    real(dp), intent(in) :: sparse(:)

    !> List of k-points.
    real(dp), intent(in) :: kPoints(:,:)

    !> Neighbour list index.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours.
    integer, intent(in) :: nNeighbourSK(:)

    !> Offset array in the square matrix.
    integer, intent(in) :: iAtomStart(:)

    !> Pair indexing array.
    integer, intent(in) :: iPair(0:,:)

    !> Mapping of the atoms to the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Index of the cell translation vectors for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Cell translation vectors.
    real(dp), intent(in) :: cellVec(:,:)

    complex(dp), allocatable :: square(:,:)
    character(mc) :: strForm
    integer :: fd, nOrb, nKPoint
    integer :: iK

    if (withMpi) then
      call error("Writing of HS not working with MPI yet")
    end if

    nOrb = iAtomStart(size(nNeighbourSK) + 1) - 1
    nKPoint = size(kPoints, dim =2)

    allocate(square(nOrb, nOrb))
    open(newunit=fd, file=fname, form="formatted", status="replace")
    write(fd, "(A1,A10,A10,A10,A10)") "#", "REAL", "NALLORB", "NKPOINT"
    write(fd, "(1X,L10,I10,I10)") .false., nOrb, nKPoint

    write (strForm, "(A,I0,A)") "(", 2 * nOrb, "ES24.15)"
    do iK = 1, nKPoint
      call unpackHS(square, sparse, kPoints(:,iK), iNeighbour, nNeighbourSK, iCellVec, cellVec,&
          & iAtomStart, iPair, img2CentCell)
      call blockHermitianHS(square, iAtomStart)
      write(fd, "(A1,A10,A10)") "#", "IKPOINT"
      write(fd, "(1X,I10,I10)") iK
      write(fd, "(A1,A)") "#", " MATRIX"
      write(fd, strForm) square
    end do
    close(fd)

  end subroutine writeSparseAsSquare_cplx


  !> Writes a sparse matrix to a file.
  subroutine writeSparse(fname, sparse, iNeighbour, nNeighbourSK, iAtomStart, iPair, img2CentCell,&
      & iCellVec, cellVec)

    !> Name of the file to write the matrix to.
    character(len=*), intent(in) :: fname

    !> Sparse matrix.
    real(dp), intent(in) :: sparse(:)

    !> Neighbour list index.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours.
    integer, intent(in) :: nNeighbourSK(:)

    !> Offset array in the square matrix.
    integer, intent(in) :: iAtomStart(:)

    !> Pair indexing array.
    integer, intent(in) :: iPair(0:,:)

    !> Mapping of the atoms to the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Index of the cell translation vectors for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Cell translation vectors.
    real(dp), intent(in) :: cellVec(:,:)

    integer :: fd, nAtom
    integer :: iAt1, iAt2, iAt2f, iNeigh, iOrig, nOrb1, nOrb2
    character(mc) :: strForm

    if (.not. tIoProc) then
      return
    end if

    nAtom = size(nNeighbourSK)

    open(newunit=fd, file=fname, form="formatted", status="replace")
    write(fd, "(A1,A10)") "#", "NATOM"
    write(fd, "(1X,I10)") nAtom
    write(fd, "(A1,A10,A10,A10)") "#", "IATOM", "NNEIGH", "NORB"
    do iAt1 = 1, nAtom
      write(fd, "(1X,I10,I10,I10)") iAt1, nNeighbourSK(iAt1) + 1, iAtomStart(iAt1+1)&
          & - iAtomStart(iAt1)
    end do

    do iAt1 = 1, nAtom
      nOrb1 = iAtomStart(iAt1+1) - iAtomStart(iAt1)
      do iNeigh = 0, nNeighbourSK(iAt1)
        iOrig = iPair(iNeigh,iAt1) + 1
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        nOrb2 = iAtomStart(iAt2f+1) - iAtomStart(iAt2f)
        write(strForm, "(A,I0,A)") "(", nOrb2, "ES24.15)"
        write(fd, "(A1,A10,A10,A10,3A10)") "#", "IATOM1", "INEIGH", "IATOM2F", "ICELL(1)",&
            & "ICELL(2)", "ICELL(3)"
        write(fd, "(1X,I10,I10,I10,3I10)") iAt1, iNeigh, iAt2f, int(cellVec(:,iCellVec(iAt2)))
        write(fd, "(A1,A)") "#", " MATRIX"
        write(fd, strForm) sparse(iOrig:iOrig+nOrb1*nOrb2-1)
      end do
    end do
    close(fd)

  end subroutine writeSparse

end module dftbp_formatout
