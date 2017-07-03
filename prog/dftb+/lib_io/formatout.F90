!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains subroutines for formatted output of data
module formatout
  use assert
  use accuracy
  use fileid
  use constants
  use lapackroutines, only: matinv
  use sparse2dense
  implicit none
  private

  public :: clearFile, writeGenFormat, writeXYZFormat
  public :: printDFTBHeader
  public :: writeSparseAsSquare, writeSparse

  !!* Clears contents of a file
  interface clearFile
    module procedure clearFile_fname
  end interface


  !!* Writes geometry information in gen format
  interface writeGenFormat
    module procedure writeGenFormat_fname
    module procedure writeGenFormat_fid
  end interface

  !!* Writes geometry information in xyz format
  interface writeXYZFormat
    module procedure writeXYZFormat_fname
    module procedure writeXYZFormat_fid
  end interface

  !!* Writes sparse matrix in square form
  interface writeSparseAsSquare
    module procedure writeSparseAsSquare_real
    module procedure writeSparseAsSquare_cplx
  end interface



contains

  !!* Clears contents of file
  !!* @param fileName name of the file which should be cleared
  subroutine clearFile_fname(fileName)
    character(len=*), intent(in)   :: fileName

    integer, save :: fd = -1

    if (fd == -1) then
      fd = getFileId()
    end if
    open(fd, file=fileName, status="replace", position="rewind")
    close(fd)

  end subroutine clearFile_fname



  !!* A wrapper around writeGenFormat_fid.
  !!* @param fileName    File name of the file which should be created
  !!* @param coord       Coordinates in atomic units
  !!* @param species      Species of the atoms
  !!* @param speciesName  Name of the different species
  !!* @param latVec      Lattice vectors
  !!* @param tFracCoord  Print out fractional coordinates?
  subroutine writeGenFormat_fname(fileName, coord, species, speciesName, &
      & latVec, tFracCoord)
    character(len=*), intent(in)   :: fileName
    real(dp),         intent(in)   :: coord(:,:)
    integer,          intent(in)   :: species(:)
    character(mc),    intent(in)   :: speciesName(:)
    real(dp), intent(in), optional :: latVec(3,3)
    logical, intent(in), optional  :: tFracCoord

    integer, save :: fd = -1

    @:ASSERT((.not.(present(tFracCoord).neqv.present(latVec))) .or.(present(latVec)))

    if (fd == -1) then
      fd = getFileId()
    end if
    open(fd, file=fileName, position="append")
    call writeGenFormat(fd, coord, species, speciesName, latVec, tFracCoord)
    close(fd)

  end subroutine writeGenFormat_fname



  !!* Writes coordinates in the famous GEN format to a file
  !!* @param fd          File id of an open file where output should be written
  !!* @param coord       Coordinates in atomic units
  !!* @param species      Species of the atoms
  !!* @param speciesName  Name of the different species
  !!* @param latVec      Lattice vectors
  !!* @param tFracCoord  Print out fractional coordinates?
  subroutine writeGenFormat_fid(fd, coord, species, speciesName, latVec, &
      & tFracCoord)
    integer,           intent(in)  :: fd
    real(dp),          intent(in)  :: coord(:,:)
    integer,           intent(in)  :: species(:)
    character(mc),     intent(in)  :: speciesName(:)
    real(dp), intent(in), optional :: latVec(:,:)
    logical, intent(in), optional  :: tFracCoord

    integer :: nAtom, nSpecies
    character(6) :: formatSpecies
    integer :: ii, jj
    logical :: tFractional
    real(dp) :: invLatVec(3,3)

100 format(I5," ",A2)
101 format("(",I2.2,"A3)")
102 format(I5,I2,3E20.10)
103 format(3E20.10)

    nAtom = size(coord, dim=2)
    nSpecies = maxval(species)

    @:ASSERT(size(coord, dim=1) == 3)
    @:ASSERT(size(species) == nAtom)
    @:ASSERT(size(speciesName) == nSpecies)
  #:call ASSERT_CODE
    if (present(latVec)) then
      @:ASSERT(all(shape(latVec) == (/3, 3 /)))
    end if
  #:endcall ASSERT_CODE
    @:ASSERT((.not.(present(tFracCoord).neqv.present(latVec))) .or.(present(latVec)))

    tFractional = .false.
    if (present(latVec)) then
      if (present(tFracCoord) ) then
        tFractional = tFracCoord
      end if
      if (tFractional) then
        write(fd, 100) nAtom, "F"
      else
        write(fd, 100) nAtom, "S"
      end if
    else
      write(fd, 100) nAtom, "C"
    end if
    write(formatSpecies, 101) nSpecies
    write(fd, formatSpecies) (trim(speciesName(ii)), ii = 1, nSpecies)

    if (tFractional) then
      invLatVec(:,:) = latVec(:,:)
      call matinv(invLatVec)
      do ii = 1, nAtom
        write(fd, 102) ii, species(ii), matmul(invLatVec,coord(:, ii))
      end do
    else
      do ii = 1, nAtom
        write(fd, 102) ii, species(ii), (coord(jj, ii) * Bohr__AA, jj = 1, 3)
      end do
    end if
    if (present(latVec)) then
      write(fd, 103) 0.0_dp, 0.0_dp, 0.0_dp
      do ii = 1, 3
        write(fd, 103) (latVec(jj, ii) * Bohr__AA, jj = 1, 3)
      end do
    end if
  end subroutine writeGenFormat_fid



  !!* Writes coordinates in the XYZ format
  !!* @param fileName    File name of a file to be created
  !!* @param coord       Coordinates in atomic units
  !!* @param species      Species of the atoms
  !!* @param speciesName  Name of the different species
  !!* @param charges     Optional vector with charges for each atom.
  !!* @param velocities  Optional array of velocity vectors for each atom.
  !!* @param comment     Optional comment for line 2 of the file
  subroutine writeXYZFormat_fname(fileName, coord, species, speciesName, &
      &charges, velocities, comment)
    character(len=*), intent(in) :: fileName
    real(dp), intent(in) :: coord(:,:)
    integer,  intent(in) :: species(:)
    character(mc), intent(in) :: speciesName(:)
    real(dp), intent(in), optional :: charges(:)
    real(dp), intent(in), optional :: velocities(:,:)
    character(len=*), intent(in), optional :: comment

    integer, save :: fd = -1

    if (fd == -1) then
       fd = getFileId()
    end if
    open(fd, file=fileName, position="append")
    call writeXYZFormat(fd, coord, species, speciesName, charges, velocities, &
        &comment)
    close(fd)

  end subroutine writeXYZFormat_fname



  !!* Writes coordinates in the XYZ format with additional charges and vectors
  !!* @param fd          File id of an open file where output should be written
  !!* @param coord       Coordinates in atomic units
  !!* @param species      Species of the atoms
  !!* @param speciesName  Name of the different species
  !!* @param charges     Optional vector with charges for each atom.
  !!* @param velocities  Optional array of velocity vectors for each atom.
  !!* @param comment     Optional comment for line 2 of the file
  subroutine writeXYZFormat_fid(fd, coords, species, speciesNames, charges, &
      &velocities, comment)
    integer,  intent(in) :: fd
    real(dp), intent(in) :: coords(:,:)
    integer,  intent(in) :: species(:)
    character(mc), intent(in) :: speciesNames(:)
    real(dp), intent(in), optional :: charges(:)
    real(dp), intent(in), optional :: velocities(:,:)
    character(len=*), intent(in), optional :: comment

    integer :: nAtom, nSpecies
    integer :: ii, jj

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
  #:call ASSERT_CODE
    if (present(charges)) then
      @:ASSERT(size(charges) == nAtom)
    end if
    if (present(velocities)) then
      @:ASSERT(all(shape(velocities) == (/ 3, nAtom /)))
    end if
  #:endcall ASSERT_CODE

    write(fd, 200) nAtom
    if (present(comment)) then
      write(fd, "(A)") trim(comment)
    elseif (present(velocities)) then
      write(fd, *) "Velocity in AA/ps"
    else
      write(fd, *) ""
    end if

    if (present(charges) .and. present(velocities)) then
      write(fd, 204) (trim(speciesNames(species(ii))), &
          &coords(:, ii) * Bohr__AA, charges(ii), &
          &velocities(:,ii) * Bohr__AA / au__fs * 1000.0_dp, ii = 1, nAtom)
    elseif (present(charges) .and. .not. present(velocities)) then
      write(fd, 203) (trim(speciesNames(species(ii))), &
          &coords(:, ii) * Bohr__AA, &
          &charges(ii), ii = 1, nAtom)
    elseif (.not. present(charges) .and. present(velocities)) then
      write(fd, 202) (trim(speciesNames(species(ii))), &
          &coords(:, ii) * Bohr__AA, &
          &velocities(:,ii) * Bohr__AA / au__fs * 1000.0_dp, ii = 1, nAtom)
    else
      write(fd, 201) (trim(speciesNames(species(ii))), &
          & (coords(jj, ii) * Bohr__AA, jj = 1, 3), ii = 1, nAtom)
    end if

  end subroutine writeXYZFormat_fid



  !!* Writes the greeting message of dftb+ on stdout
  !!* @param revision Revision string from svn
  !!* @param headURL URL of the head (from svn)
  !!* @param parserVersion Version of the current parser
  subroutine printDFTBHeader(release, year)
    character(len=*), intent(in) :: release
    integer, intent(in) :: year

    character, parameter :: vbar = '|'
    character, parameter :: hbar = '='
    integer, parameter :: headerWidth = 80

    write(*, '(2A,/,A)') vbar, repeat(hbar, headerWidth - 1), vbar
    write(*, '(4A)') vbar, '  DFTB+ (Release ', release, ')'
    write(*, '(A)') vbar
    write(*, '(2A,I0,A)') vbar, '  Copyright (C) ', year, '  DFTB+ developers group'
    write(*, '(A,/,2A,/,A)') vbar, vbar, repeat(hbar, headerWidth - 1), vbar
    write(*, '(2A)') vbar,&
        & '  When publishing results obtained with DFTB+, please cite the following',&
        & vbar, '  reference:'
    write(*, '(A)') vbar
    write(*, '(2A)') vbar,'  * B. Aradi, B. Hourahine and T. Frauenheim,',&
        & vbar, '    DFTB+, a Sparse Matrix-Based Implementation of the DFTB Method,',&
        & vbar, '    J. Phys. Chem. A, 111 5678 (2007).  [doi: 10.1021/jp070186p]'
    write(*, '(A)') vbar
    write(*, '(2A,2(/,2A))') vbar,&
        & '  You should also cite additional publications crediting the parametrization',&
        & vbar,&
        & '  data you use. Please consult the documentation of the SK-files for the',&
        & vbar,&
        & '  references.'
    write(*, '(A,/,2A,/)') vbar, vbar, repeat(hbar, headerWidth - 1)

  end subroutine printDFTBHeader


  !!* Converts a sparse matrix to its square form and writes to a file.
  !!* @param fname Name of the file to write the matrix to.
  !!* @param sparse Sparse matrix.
  !!* @param iNeighbor Neighbor list index.
  !!* @param nNeighbor Number of neighbors.
  !!* @param iAtomStart Offset array in the square matrix.
  !!* @param iPair Pair indexing array.
  !!* @param img2CentCell Mapping of the atoms to the central cell.
  subroutine writeSparseAsSquare_real(fname, sparse, iNeighbor, nNeighbor, &
      &iAtomStart, iPair, img2CentCell)
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: sparse(:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:)

    real(dp), allocatable :: square(:,:)
    character(mc) :: strForm
    integer :: fd, nOrb

    nOrb = iAtomStart(size(nNeighbor) + 1) - 1

    allocate(square(nOrb, nOrb))
    fd = getFileId()
    open(fd, file=fname, form="formatted", status="replace")
    write(fd, "(A1,A10,A10,A10,A10)") "#", "REAL", "NALLORB", "NKPOINT"
    write(fd, "(1X,L10,I10,I10,I10)") .true., nOrb, 1

    write (strForm, "(A,I0,A)") "(", nOrb, "ES24.15)"
    call unpackHS(square, sparse, iNeighbor, nNeighbor, iAtomStart, iPair, &
        &img2CentCell)
    call blockSymmetrizeHS(square, iAtomStart)
    write(fd, "(A1,A10,A10)") "#", "IKPOINT"
    write(fd, "(1X,I10,I10)") 1
    write(fd, "(A1,A)") "#", " MATRIX"
    write(fd, strForm) square
    close(fd)

  end subroutine writeSparseAsSquare_real



  !!* Converts a sparse matrix to its square form and writes to a file.
  !!* @param fname Name of the file to write the matrix to.
  !!* @param sparse Sparse matrix.
  !!* @param kPoints List of k-points.
  !!* @param iNeighbor Neighbor list index.
  !!* @param nNeighbor Number of neighbors.
  !!* @param iAtomStart Offset array in the square matrix.
  !!* @param iPair Pair indexing array.
  !!* @param img2CentCell Mapping of the atoms to the central cell.
  !!* @param iCellVec Index of the cell translation vectors for each atom.
  !!* @param cellVec Cell translation vectors.
  subroutine writeSparseAsSquare_cplx(fname, sparse, kPoints, iNeighbor, &
      &nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, cellVec)
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: sparse(:)
    real(dp), intent(in) :: kPoints(:,:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)

    complex(dp), allocatable :: square(:,:)
    character(mc) :: strForm
    integer :: fd, nOrb, nKPoint
    integer :: iK

    nOrb = iAtomStart(size(nNeighbor) + 1) - 1
    nKPoint = size(kPoints, dim =2)

    allocate(square(nOrb, nOrb))
    fd = getFileId()
    open(fd, file=fname, form="formatted", status="replace")
    write(fd, "(A1,A10,A10,A10,A10)") "#", "REAL", "NALLORB", "NKPOINT"
    write(fd, "(1X,L10,I10,I10)") .false., nOrb, nKPoint

    write (strForm, "(A,I0,A)") "(", 2 * nOrb, "ES24.15)"
    do iK = 1, nKPoint
      call unpackHS(square, sparse, kPoints(:,iK), iNeighbor, nNeighbor, &
          &iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
      call blockHermitianHS(square, iAtomStart)
      write(fd, "(A1,A10,A10)") "#", "IKPOINT"
      write(fd, "(1X,I10,I10)") iK
      write(fd, "(A1,A)") "#", " MATRIX"
      write(fd, strForm) square
    end do
    close(fd)

  end subroutine writeSparseAsSquare_cplx



  !!* Writes a sparse matrix to a file.
  !!* @param fname Name of the file to write the matrix to.
  !!* @param sparse Sparse matrix.
  !!* @param iNeighbor Neighbor list index.
  !!* @param nNeighbor Number of neighbors.
  !!* @param iAtomStart Offset array in the square matrix.
  !!* @param iPair Pair indexing array.
  !!* @param img2CentCell Mapping of the atoms to the central cell.
  !!* @param iCellVec Index of the cell translation vectors for each atom.
  !!* @param cellVec Cell translation vectors.
  subroutine writeSparse(fname, sparse, iNeighbor, nNeighbor, iAtomStart, &
      &iPair, img2CentCell, iCellVec, cellVec)
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: sparse(:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)

    integer :: fd, nAtom
    integer :: iAt1, iAt2, iAt2f, iNeigh, iOrig, nOrb1, nOrb2
    character(mc) :: strForm

    nAtom = size(nNeighbor)

    fd = getFileId()
    open(fd, file=fname, form="formatted", status="replace")
    write(fd, "(A1,A10)") "#", "NATOM"
    write(fd, "(1X,I10)") nAtom
    write(fd, "(A1,A10,A10,A10)") "#", "IATOM", "NNEIGH", "NORB"
    do iAt1 = 1, nAtom
      write(fd, "(1X,I10,I10,I10)") iAt1, nNeighbor(iAt1) + 1, &
          &iAtomStart(iAt1+1) - iAtomStart(iAt1)
    end do

    do iAt1 = 1, nAtom
      nOrb1 = iAtomStart(iAt1+1) - iAtomStart(iAt1)
      do iNeigh = 0, nNeighbor(iAt1)
        iOrig = iPair(iNeigh,iAt1) + 1
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        nOrb2 = iAtomStart(iAt2f+1) - iAtomStart(iAt2f)
        write(strForm, "(A,I0,A)") "(", nOrb2, "ES24.15)"
        write(fd, "(A1,A10,A10,A10,3A10)") "#", "IATOM1", "INEIGH", "IATOM2F", &
            &"ICELL(1)", "ICELL(2)", "ICELL(3)"
        write(fd, "(1X,I10,I10,I10,3I10)") iAt1, iNeigh, iAt2f, &
            &int(cellVec(:,iCellVec(iAt2)))
        write(fd, "(A1,A)") "#", " MATRIX"
        write(fd, strForm) sparse(iOrig:iOrig+nOrb1*nOrb2-1)
      end do
    end do
    close(fd)

  end subroutine writeSparse



end module formatout

