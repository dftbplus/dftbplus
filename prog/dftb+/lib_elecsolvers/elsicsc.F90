!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines for converting from and to ELSI CSC format.
module dftbp_elsicsc
#:if WITH_MPI
  use dftbp_mpifx
#:endif
  use dftbp_accuracy, only : dp
  use dftbp_assert
  use dftbp_environment, only : TEnvironment
  use dftbp_periodic, only : TNeighbourList
  use dftbp_constants, only : pi
  use dftbp_message, only : error
  use dftbp_angmomentum, only : rotateZ
  use dftbp_commontypes, only : TOrbitals
  implicit none
  private

  public :: TElsiCsc, TElsiCsc_init

  #:set CONDITIONS = ['', 'Helical']

  !> Data needed to convert between DFTB+ sparse format and ELSI CSC format.
  type :: TElsiCsc
    private

    !> number of non-zero matrix elements on this processor
    integer, public :: nnzLocal

    !>  number of non-zero matrix elements in the whole sparse matrix
    integer, public :: nnzGlobal

    !> Number of local columns on this processor
    integer, public :: numColLocal

    !> Local column pointer in the CSC format on this processor
    integer, public, allocatable :: colPtrLocal(:)

    !> Local row index for non-zero elements (note: initialized at first pack2elsi() call)
    integer, public, allocatable :: rowIndLocal(:)

    !> On which column does this processor start its matrix
    integer :: colStartLocal

    !> On which column does this processor end its matrix
    integer :: colEndLocal

    !> Index for starting row of blocks in nzValLocal
    integer, allocatable :: blockRow(:,:)

    !> List of atoms with elements in the columns held locally
    integer, allocatable :: atomsInColumns(:)

    !> Count of the atoms with elements in the columns held locally
    integer :: nAtomsInColumns

  contains

    procedure, private :: TElsiCsc_convertPackedToElsiReal
    procedure, private :: TElsiCsc_convertPackedHelicalToElsiReal
    procedure, private :: TElsiCsc_convertPackedToElsiCmplx
    procedure, private :: TElsiCsc_convertPackedHelicalToElsiCmplx
    procedure, private :: TElsiCsc_convertElsiToPackedReal
    procedure, private :: TElsiCsc_convertElsiToPackedHelicalReal
    procedure, private :: TElsiCsc_convertElsiToPackedCmplx
    procedure, private :: TElsiCsc_convertElsiToPackedHelicalCmplx
    generic :: convertPackedToElsiReal => TElsiCsc_convertPackedToElsiReal,&
        & TElsiCsc_convertPackedHelicalToElsiReal
    generic :: convertPackedToElsiCmplx => TElsiCsc_convertPackedToElsiCmplx,&
        & TElsiCsc_convertPackedHelicalToElsiCmplx
    generic :: convertElsiToPackedReal => TElsiCsc_convertElsiToPackedReal,&
        & TElsiCsc_convertElsiToPackedHelicalReal
    generic :: convertElsiToPackedCmplx => TElsiCsc_convertElsiToPackedCmplx,&
        & TElsiCsc_convertElsiToPackedHelicalCmplx

  end type TElsiCsc

#:if WITH_ELSI

  #:set TYPES = ['real', 'complex']

  ! Internal routines addding blocks to ELSI matrices
  interface addBlock2Elsi
    #:for TYPE in TYPES
      module procedure addBlock2Elsi${TYPE}$
    #:endfor
  end interface addBlock2Elsi

  ! Internal routines copying blocks from  ELSI matrices
  interface cpElsi2Block
    #:for TYPE in TYPES
      module procedure cpElsi2Block${TYPE}$
    #:endfor
  end interface cpElsi2Block

#:endif


contains

  !> Initializes the ELSI CSC converter.
  subroutine TElsiCsc_init(this, env, nBasis, csrBlockSize, neighList, nNeighbourSK, iAtomStart,&
      & iSparseStart, img2CentCell)

    !> Sparse conversion instance
    type(TElsiCsc), intent(out) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Number of basis functions
    integer, intent(in) :: nBasis

    !> Number of csr rows on a given process.
    integer, intent(in) :: csrBlockSize

    !> Neighbour list for the atoms (First index from 0!)
    type(TNeighbourList), intent(in) :: neighList

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

  #:if WITH_ELSI

    integer :: numCol, iAtom1, iAtom2, iAtom2f, nAtom, ii, jj, iNeigh, iNext, nOrb1, nOrb2
    integer :: iOrig
    integer, allocatable :: blockList(:,:)
    logical, allocatable :: tRowTrans(:)

    numCol = nBasis

    this%colStartLocal = env%mpi%groupComm%rank * csrBlockSize + 1
    if (env%mpi%groupComm%rank /=  env%mpi%groupComm%size - 1) then
      this%numColLocal = csrBlockSize
    else
      this%numColLocal = numCol - (env%mpi%groupComm%size - 1) * csrBlockSize
    end if
    this%colEndLocal = this%colStartLocal + this%numColLocal - 1

    allocate(this%colPtrLocal(this%numColLocal + 1))

    call getColPtr(neighList%iNeighbour, nNeighbourSK, iAtomStart, iSparseStart, img2CentCell,&
        & this%colStartLocal, this%colEndLocal, this%nnzLocal, this%colPtrLocal)

    this%nnzGlobal = 0
    call mpifx_allreduce(env%mpi%groupComm, this%nnzLocal, this%nnzGlobal, MPI_SUM)

    call getBlockRow(neighList%iNeighbour, nNeighbourSK, iAtomStart, iSparseStart, img2CentCell,&
      & this%colStartLocal, this%colEndLocal, this%atomsInColumns, this%nAtomsInColumns,&
      & this%blockRow)

  #:else

    call error("Internal error: TElsiCsc_init() called despite missing ELSI support")

  #:endif

  end subroutine TElsiCsc_init

#:for BCS in CONDITIONS

  !> Convert sparse DFTB+ matrix to distributed CSC matrix format for ELSI calculations for real
  !> matrices
  !>
  !> NOTE: ELSI needs the full matrix (both triangles)
  !>
  subroutine TElsiCsc_convertPacked${BCS}$ToElsiReal(this, orig, iNeighbour, nNeighbourSK,&
      & iAtomStart, iSparseStart, img2CentCell, nzValLocal, orb&
    #:if BCS == 'Helical'
      &, species, coord&
    #:endif
    &)

    !> Sparse conversion instance
    class(TElsiCsc), intent(inout) :: this

    !> Sparse Hamiltonian
    real(dp), intent(in) :: orig(:)

    !> Neighbour list for each atom (first index from 0!).
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itthis).
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    !> Local non-zero elements
    real(dp), intent(out) :: nzValLocal(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

  #:if BCS == 'Helical'

    !> Species of each atom
    integer :: species(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coord(:,:)

  #:endif

  #:if WITH_ELSI

    integer :: nAtom
    integer :: iOrig, ii, jj, nOrb1, nOrb2
    integer :: iNeigh
    integer :: iAt, iAtom1, iAtom2, iAtom2f

    integer :: iNext
    integer, allocatable :: blockList(:,:)
    logical, allocatable :: tRowTrans(:)

    real(dp) :: tmpSqr(orb%mOrb,orb%mOrb)
  #:if BCS == 'Helical'
    real(dp) :: rotZ(orb%mOrb,orb%mOrb), theta
    integer  :: lShellVals(orb%mShell), iSh, iSp
  #:endif

    nAtom = size(iNeighbour, dim=2)

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(nNeighbourSK) == nAtom)
    @:ASSERT(size(iAtomStart) == nAtom + 1)

    allocate(blockList(nAtom,2))
    allocate(tRowTrans(nAtom))
    if (.not. allocated(this%rowIndLocal)) then
      allocate(this%rowIndLocal(this%nnzLocal))
      this%rowIndLocal(:) = 0
    end if
    nzValLocal(:) = 0.0_dp

    ! Offset in column belonging to transposed triangle
    blockList(:,2) = 1

  #:if BCS == 'Helical'
    lShellVals(:) = 0
    rotZ(:,:) = 0.0_dp
  #:endif

    ! loop over atoms relevant to this processor
    do iAt = 1, this%nAtomsInColumns
      iAtom1 = this%atomsInColumns(iAt)
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      ! Offset in current column
      blockList(:,1) = 0
      tRowTrans = .false.
      ! Starting index for column in DFTB+ sparse structure, because column probaly already contains
      ! blocks coming from transposing previously processed elements.
      iNext = blockList(iAtom1, 2)
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)

        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f+1) - jj

        if (blockList(iAtom2f,1) == 0) then
          blockList(iAtom2f,1) = iNext
          iNext = iNext + nOrb2
        end if

        if ( .not. isBlockInLocal(jj, jj + nOrb2 - 1, ii, ii + nOrb1 - 1, this%colStartLocal,&
            & this%colEndLocal)) then
          cycle
        end if

        tmpSqr(:nOrb2,:nOrb1) = reshape(orig(iOrig : iOrig + nOrb1 * nOrb2 - 1), [nOrb2, nOrb1])

      #:if BCS == 'Helical'
        iSp = species(iAtom2f)
        iSh = orb%nShell(iSp)
        lShellVals(:iSh) = orb%angShell(:iSh,iSp)
        theta = atan2(coord(2,iAtom2f),coord(1,iAtom2f)) - atan2(coord(2,iAtom2),coord(1,iAtom2))
        theta = mod(theta,2.0_dp*pi)
        call rotateZ(rotZ, lShellVals(:iSh), theta)
        tmpSqr(:nOrb2,:nOrb1) = matmul(rotZ(:nOrb2,:nOrb2),tmpSqr(:nOrb2,:nOrb1))
      #:endif

        call addBlock2Elsi(tmpSqr(:nOrb2,:nOrb1), this%colStartLocal, this%colEndLocal, jj, ii,&
            & this%colPtrLocal, blockList(iAtom2f, 1) - 1, nzValLocal, this%rowIndLocal)

        if (ii == jj) then
          cycle
        end if

        ! Because of folding of periodic images, it can happen that the transposed block has already
        ! been registered.
        if (.not. tRowTrans(iAtom2f)) then
          blockList(iAtom2f,2) = blockList(iAtom2f,2) + nOrb1
          tRowTrans(iAtom2f) = .true.
        end if

        call addBlock2Elsi(transpose(tmpSqr(:nOrb2,:nOrb1)), this%colStartLocal, this%colEndLocal,&
            & ii, jj, this%colPtrLocal, blockList(iAtom2f,2)-nOrb1-1, nzValLocal, this%rowIndLocal)

      end do
    end do

  #:else

    call error("Internal error: TElsiCsc_convertPackedToElsiReal() called despite missing ELSI&
        & support")

  #:endif

  end subroutine TElsiCsc_convertPacked${BCS}$ToElsiReal


  !> Convert sparse DFTB+ matrix to distributed CSC matrix format for ELSI calculations for complex
  !> matrices
  !>
  !> NOTE: ELSI needs the full matrix (both triangles)
  !>
  subroutine TElsiCsc_convertPacked${BCS}$ToElsiCmplx(this, orig, iNeighbour, nNeighbourSK,&
      & iAtomStart, iSparseStart, img2CentCell, kPoint, iCellVec, cellVec, nzValLocal, orb&
    #:if BCS == 'Helical'
      &, species, coord&
    #:endif
      &)

    !> Sparse conversion instance
    class(TElsiCsc), intent(inout) :: this

    !> Sparse Hamiltonian
    real(dp), intent(in) :: orig(:)

    !> Neighbour list for each atom (first index from 0!).
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itthis).
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    !> Current k-point
    real(dp), intent(in) :: kPoint(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Local non-zero elements
    complex(dp), intent(out) :: nzValLocal(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

  #:if BCS == 'Helical'

    !> Species of each atom
    integer :: species(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coord(:,:)

  #:endif

  #:if WITH_ELSI

    integer :: nAtom
    integer :: iOrig, ii, jj, nOrb1, nOrb2
    integer :: iNeigh
    integer :: iAt, iAtom1, iAtom2, iAtom2f

    integer :: iNext, iVec
    integer, allocatable :: blockList(:,:)
    logical, allocatable :: tRowTrans(:)
    real(dp) :: kPoint2p(size(kpoint)), tmpSqr(orb%mOrb,orb%mOrb)
    complex(dp) :: phase
  #:if BCS == 'Helical'
    real(dp) :: rotZ(orb%mOrb,orb%mOrb), theta
    integer  :: lShellVals(orb%mShell), iSh, iSp
  #:endif

    kPoint2p(:) = 2.0_dp * pi * kPoint(:)
    nAtom = size(iNeighbour, dim=2)

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(nNeighbourSK) == nAtom)
    @:ASSERT(size(iAtomStart) == nAtom + 1)

    allocate(blockList(nAtom,2))
    allocate(tRowTrans(nAtom))
    if (.not. allocated(this%rowIndLocal)) then
      allocate(this%rowIndLocal(this%nnzLocal))
      this%rowIndLocal(:) = 0
    end if
    nzValLocal(:) = cmplx(0,0,dp)

    ! Offset in column belonging to transposed triangle
    blockList(:,2) = 1

    ! loop over atoms relevant to this processor
    do iAt = 1, this%nAtomsInColumns
      iAtom1 = this%atomsInColumns(iAt)
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      ! Offset in current column
      blockList(:,1) = 0
      tRowTrans = .false.
      ! Starting index for column in DFTB+ sparse structure, because column probably already
      ! contains blocks coming from transposing previously processed elements.
      iNext = blockList(iAtom1, 2)
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)

        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f+1) - jj

        if (blockList(iAtom2f,1) == 0) then
          blockList(iAtom2f,1) = iNext
          iNext = iNext + nOrb2
        end if

        if (.not. isBlockInLocal(jj, jj + nOrb2 - 1, ii, ii + nOrb1 - 1, this%colStartLocal,&
            & this%colEndLocal)) then
          cycle
        end if

        iVec = iCellVec(iAtom2)
        phase = exp(cmplx(0, 1, dp) * dot_product(kPoint2p, cellVec(:, iVec)))

        tmpSqr(:nOrb2, :nOrb1) = reshape(orig(iOrig : iOrig + nOrb1 * nOrb2 - 1), [nOrb2, nOrb1])

      #:if BCS == 'Helical'
        iSp = species(iAtom2f)
        iSh = orb%nShell(iSp)
        lShellVals(:iSh) = orb%angShell(:iSh,iSp)
        theta = atan2(coord(2,iAtom2f),coord(1,iAtom2f)) - atan2(coord(2,iAtom2),coord(1,iAtom2))
        theta = mod(theta,2.0_dp*pi)
        call rotateZ(rotZ, lShellVals(:iSh), theta)
        tmpSqr(:nOrb2,:nOrb1) = matmul(rotZ(:nOrb2,:nOrb2),tmpSqr(:nOrb2,:nOrb1))
      #:endif

        call addBlock2Elsi(phase * tmpSqr(:nOrb2, :nOrb1), this%colStartLocal, this%colEndLocal,&
            & jj, ii, this%colPtrLocal, blockList(iAtom2f, 1) - 1, nzValLocal, this%rowIndLocal)

        if (ii == jj) then
          cycle
        end if

        ! other triangle, so Hermitian symmetry
        phase = conjg(phase)

        ! Because of folding of periodic images, it can happen that the transposed block has already
        ! been registered.
        if (.not. tRowTrans(iAtom2f)) then
          blockList(iAtom2f,2) = blockList(iAtom2f,2) + nOrb1
          tRowTrans(iAtom2f) = .true.
        end if

        call addBlock2Elsi(phase * transpose(tmpSqr(:nOrb2, :nOrb1)), this%colStartLocal,&
            & this%colEndLocal, ii, jj, this%colPtrLocal, blockList(iAtom2f,2) - nOrb1 - 1,&
            & nzValLocal, this%rowIndLocal)

      end do
    end do

  #:else

    call error("Internal error: TElsiCsc_convertPackedToElsiCmplx() called despite missing ELSI&
        & support")

  #:endif

  end subroutine TElsiCsc_convertPacked${BCS}$ToElsiCmplx


  !> Convert CSC matrix format into DFTB+ sparse format
  !>
  !> Note: primitive will not be set to zero on startup, and values are added to enable addition of
  !> spin components. Make sure, you set it to zero before invoking this routine the first time.
  subroutine TElsiCsc_convertElsiToPacked${BCS}$Real(this, iNeighbour, nNeighbourSK, orb,&
    #:if BCS == 'Helical'
      & species, coord,&
    #:endif
      & iAtomStart, iSparseStart, img2CentCell, nzval, primitive)

    !> Sparse conversion instance
    class(TElsiCsc), intent(inout) :: this

    !> Neighbour list for each atom (first index from 0!).
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itthis).
    integer, intent(in) :: nNeighbourSK(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

  #:if BCS == 'Helical'

    !> Species of each atom
    integer :: species(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coord(:,:)

  #:endif

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    !> Local non-zero elements
    real(dp), intent(in) :: nzval(:)

    !> Sparse Hamiltonian
    real(dp), intent(inout) :: primitive(:)

  #:if WITH_ELSI

    integer :: nAtom
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh
    integer :: iAt, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: tmpSqr(orb%mOrb,orb%mOrb)
  #:if BCS == 'Helical'
    real(dp) :: rotZ(orb%mOrb,orb%mOrb), theta
    integer  :: lShellVals(orb%mShell), iSh, iSp
  #:endif

    nAtom = size(iNeighbour, dim=2)

    tmpSqr(:,:) = 0.0_dp
  #:if BCS == 'Helical'
    lShellVals(:) = 0
    rotZ(:,:) = 0.0_dp
  #:endif

    ! loop over relevant atoms to pack back
    do iAt = 1, this%nAtomsInColumns
      iAtom1 = this%atomsInColumns(iAt)
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        nOrb2 = iAtomStart(iAtom2f+1) - jj

        if (.not. isBlockInLocal(jj, jj + nOrb2 - 1, ii, ii + nOrb1 - 1, this%colStartLocal,&
            & this%colEndLocal)) then
          cycle
        end if

        call cpElsi2Block(this%colStartLocal, this%colEndLocal, this%colPtrLocal, nzval, ii,&
            & this%blockRow(iNeigh, iAtom1), tmpSqr(1:nOrb2,1:nOrb1))

        ! Symmetrize the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = tmpSqr(kk+1:nOrb1, kk)
          end do
        end if

      #:if BCS == 'Helical'
        iSp = species(iAtom2f)
        iSh = orb%nShell(iSp)
        lShellVals(:iSh) = orb%angShell(:iSh,iSp)
        theta = atan2(coord(2,iAtom2),coord(1,iAtom2)) - atan2(coord(2,iAtom2f),coord(1,iAtom2f))
        theta = mod(theta,2.0_dp*pi)
        call rotateZ(rotZ, lShellVals(:iSh), theta)
        tmpSqr(:nOrb2,:nOrb1) =  matmul(rotZ(:nOrb2,:nOrb2),tmpSqr(:nOrb2,:nOrb1))
      #:endif

        primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1) = primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1)&
            & + reshape(tmpSqr(1:nOrb2,1:nOrb1), [nOrb1*nOrb2])

      end do
    end do

  #:else

    call error("Internal error: TElsiCsc_convertElsiToPackedReal() called despite missing ELSI&
        & support")

  #:endif

  end subroutine TElsiCsc_convertElsiToPacked${BCS}$Real


  !> Convert CSC matrix format into DFTB+ sparse format
  !>
  !> Note: primitive will not be set to zero on startup, and values are added to enable addition of
  !> spin components. Make sure, you set it to zero before invoking this routine the first time.
  subroutine TElsiCsc_convertElsiToPacked${BCS}$Cmplx(this, iNeighbour, nNeighbourSK, orb,&
    #:if BCS == 'Helical'
      & species, coord,&
    #:endif
      & iAtomStart, iSparseStart, img2CentCell, kPoint, kWeight, iCellVec, cellVec, nzval,&
      & primitive)

    !> Sparse conversion instance
    class(TElsiCsc), intent(in) :: this

    !> Neighbour list for each atom (first index from 0!).
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itthis).
    integer, intent(in) :: nNeighbourSK(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

  #:if BCS == 'Helical'

    !> Species of each atom
    integer :: species(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coord(:,:)

  #:endif

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    !> Current k-point
    real(dp), intent(in) :: kPoint(:)

    !> Weight for current k-points
    real(dp), intent(in) :: kWeight

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Local non-zero elements
    complex(dp), intent(in) :: nzval(:)

    !> Sparse Hamiltonian
    real(dp), intent(inout) :: primitive(:)

  #:if WITH_ELSI

    integer :: nAtom, iVec
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh
    integer :: iAt, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    complex(dp) :: tmpSqr(orb%mOrb,orb%mOrb)
    real(dp) :: kPoint2p(size(kpoint))
    complex(dp) :: phase
  #:if BCS == 'Helical'
    real(dp) :: rotZ(orb%mOrb,orb%mOrb), theta
    integer  :: lShellVals(orb%mShell), iSh, iSp
  #:endif

    kPoint2p(:) = 2.0_dp * pi * kPoint(:)
    nAtom = size(iNeighbour, dim=2)

    tmpSqr(:,:) = cmplx(0,0,dp)
  #:if BCS == 'Helical'
    lShellVals(:) = 0
    rotZ(:,:) = 0.0_dp
  #:endif

    ! loop over relevant atoms to pack back
    do iAt = 1, this%nAtomsInColumns
      iAtom1 = this%atomsInColumns(iAt)
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        nOrb2 = iAtomStart(iAtom2f+1) - jj

        if (.not. isBlockInLocal(jj, jj + nOrb2 - 1, ii, ii + nOrb1 - 1, this%colStartLocal,&
            & this%colEndLocal)) then
          cycle
        end if

        call cpElsi2Block(this%colStartLocal, this%colEndLocal, this%colPtrLocal, nzval, ii,&
            & this%blockRow(iNeigh, iAtom1), tmpSqr(1:nOrb2,1:nOrb1))

        iVec = iCellVec(iAtom2)
        phase = exp(cmplx(0, -1, dp) * dot_product(kPoint2p, cellVec(:, iVec)))

        ! Hermitian conjugate the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
          end do
        end if

      #:if BCS == 'Helical'
        iSp = species(iAtom2f)
        iSh = orb%nShell(iSp)
        lShellVals(:iSh) = orb%angShell(:iSh,iSp)
        theta = atan2(coord(2,iAtom2),coord(1,iAtom2)) - atan2(coord(2,iAtom2f),coord(1,iAtom2f))
        theta = mod(theta,2.0_dp*pi)
        call rotateZ(rotZ, lShellVals(:iSh), theta)
        tmpSqr(:nOrb2,:nOrb1) =  matmul(rotZ(:nOrb2,:nOrb2),tmpSqr(:nOrb2,:nOrb1))
      #:endif

        primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1) = primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1)&
            & + kWeight * real(phase * reshape(tmpSqr(1:nOrb2,1:nOrb1), [nOrb1*nOrb2]))

      end do
    end do

  #:else

    call error("Internal error: TElsiCsc_convertElsiToPackedCmplx() called despite missing ELSI&
        & support")

  #:endif

  end subroutine TElsiCsc_convertElsiToPacked${BCS}$Cmplx

#:endfor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:if WITH_ELSI


  !> Checks if atom block is part of local matrix (Elsi)
  pure logical function isBlockInLocal(rowStartBlock, rowEndBlock, colStartBlock, colEndBlock,&
      & colStartLocal, colEndLocal)

    !> Row of global matrix where block starts
    integer, intent(in) :: rowStartBlock

    !> Row of global matrix where block ends
    integer, intent(in) :: rowEndBlock

    !> Column of global matrix where block starts
    integer, intent(in) :: colStartBlock

    !> Column of global matrix where block ends
    integer, intent(in) :: colEndBlock

    !> Column of global matrix where local part starts
    integer, intent(in) :: colStartLocal

    !> Column of global matrix where local part ends
    integer, intent(in) :: colEndLocal

    isBlockInLocal = ((colStartLocal <= colStartBlock) .and. (colStartBlock <= colEndLocal))&
        & .or. ((colStartLocal <= colEndBlock) .and. (colEndBlock <= colEndLocal))&
        & .or. ((colStartBlock <= colStartLocal) .and. (colStartLocal <= colEndBlock))&
        & .or. ((colStartBlock <= colEndLocal) .and. (colEndLocal <= colEndBlock))&
        & .or. ((colStartLocal <= rowStartBlock) .and. (rowStartBlock <= colEndLocal))&
        & .or. ((colStartLocal <= rowEndBlock) .and. (rowEndBlock <= colEndLocal))&
        & .or. ((rowStartBlock <= colStartLocal) .and. (colStartLocal <= rowEndBlock))&
        & .or. ((rowStartBlock <= colEndLocal) .and. (colEndLocal <= rowEndBlock))

  end function isBlockInLocal


  !> Checks whether global column is in local matrix
  pure logical function isColumnInLocal(col, colStartLocal, colEndLocal)

    !> global column
    integer, intent(in) :: col

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStartLocal

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEndLocal

    if ((colStartLocal <= col) .and. (col <= colEndLocal)) then
      isColumnInLocal = .true.
    else
      isColumnInLocal = .false.
    end if

  end function isColumnInLocal


#:for TYPE in TYPES

  !> Add the content of a local matrix block to ELSI CSC format
  subroutine addBlock2Elsi${TYPE}$(loc, colStart, colEnd, ii, jj, colptr, rowOffset, val,&
      & rowIndLocal)

    !> Local matrix.
    ${TYPE}$(dp), intent(in) :: loc(:,:)

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStart

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEnd

    !> Starting row in the global matrix.
    integer, intent(in) :: ii

    !> Starting column in the global matrix
    integer, intent(in) :: jj

    !> column pointer
    integer, intent(in) :: colptr(:)

    !> index of the next element per column
    integer, intent(in) :: rowOffset

    !> values in CSC format
    ${TYPE}$(dp), intent(inout) :: val(:)

    !> row index pointer
    integer, intent(inout), optional :: rowIndLocal(:)

    integer :: j2, iloc, jloc
    integer :: jStart, jEnd

    jStart = max(jj, colStart) - colStart + 1
    jEnd = min(jj + size(loc, dim=2) - 1, colEnd) - colStart + 1
    jloc = max(jj, colStart) - jj + 1

    if (present(rowIndLocal)) then
      do j2 = jStart, jEnd
        do iloc = 1, size(loc, dim=1)
          val(colptr(j2) + rowOffset + iloc - 1) = val(colptr(j2) + rowOffset + iloc - 1)&
              & + loc(iloc, jloc)
          rowIndLocal(colptr(j2) + rowOffset + iloc - 1) = iloc + ii - 1
        end do
        jloc = jloc + 1
      end do
    else
      do j2 = jStart, jEnd
        do iloc = 1, size(loc, dim=1)
          val(colptr(j2) + rowOffset + iloc - 1) = val(colptr(j2) + rowOffset+iloc - 1)&
              & + loc(iloc, jloc)
        end do
        jloc = jloc + 1
      end do
    end if

  end subroutine addBlock2Elsi${TYPE}$


  !> Copies the content from the ELSI structure to block
  subroutine cpElsi2Block${TYPE}$(colStart, colEnd, colptr, nzval, jj, blockRow, loc)

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStart

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEnd

    !> column pointer
    integer, intent(in) :: colptr(:)

    !> non-zero values
    ${TYPE}$(dp), intent(in) :: nzval(:)

    !> Starting column in the global matrix
    integer, intent(in) :: jj

    !> Saves starting row of blocks in nzValLocal
    integer, intent(in) :: blockRow

    !> Local block of matrix.
    ${TYPE}$(dp), intent(out) :: loc(:,:)

    integer :: j2, iloc

    loc(:,:) = 0.0_dp

    do j2 = 1, size(loc, dim=2)
      if ( isColumnInLocal(jj + j2 - 1, colStart, colEnd )) then
        iloc = blockRow + colptr(jj-colStart+j2) - 1
        loc(:, j2) = nzval(iloc:iloc+size(loc,dim=1)-1)
      end if
    end do

  end subroutine cpElsi2Block${TYPE}$

#:endfor


  !> Creating colptr and nnz for CSC matrix format from packed format
  subroutine getColPtr(iNeighbour, nNeighbourSK, iAtomStart, iSparseStart,&
      & img2CentCell, colStartLocal, colEndLocal, nnzLocal, colPtrLocal)

    !> Neighbour list for each atom (first index from 0!).
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Atomic mapping indexes
    integer, intent(in) :: img2CentCell(:)

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStartLocal

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEndLocal

    !> Local number of  non-zero elements
    integer, intent(out) :: nnzLocal

    !> Local column pointer
    integer, intent(out) :: colPtrLocal(:)

    integer :: nAtom
    integer :: iOrig, ii, jj, nOrb1, nOrb2
    integer :: iNeigh
    integer :: iAtom1, iAtom2, iAtom2f, kk

    logical, allocatable :: blockList(:)

    nAtom = size(iNeighbour, dim=2)

    colPtrLocal(:) = 0
    nnzLocal = 0

    allocate(blockList(nAtom))

    ! Loop over all atom blocks
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      blockList(:) = .false.
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        nOrb2 = iAtomStart(iAtom2f+1) - jj

        if (blockList(iAtom2f)) then
          cycle
        else
          blockList(iAtom2f) = .true.
        end if

        if (isBlockInLocal(jj, jj+nOrb2-1, ii, ii+nOrb1-1, colStartLocal, colEndLocal)) then
          do kk = ii, ii + nOrb1 - 1
            if (isColumnInLocal(kk, colStartLocal, colEndLocal)) then
              nnzLocal = nnzLocal + nOrb2
              colPtrLocal(kk - colStartLocal + 2) = colPtrLocal(kk - colStartLocal + 2) + nOrb2
            end if
          end do
        end if

        if (ii == jj) then
          ! on the diagonal, can not be in other triangle
          cycle
        end if

        do kk = jj, jj + nOrb2 - 1
          if (isColumnInLocal(kk, colStartLocal, colEndLocal)) then
            nnzLocal = nnzLocal + nOrb1
            colPtrLocal(kk -colStartLocal + 2) = colPtrLocal(kk - colStartLocal + 2) + nOrb1
          end if
        end do
      end do
    end do
    colPtrLocal(1) = 1
    do ii = 2, size(colPtrLocal)
      colPtrLocal(ii) = colPtrLocal(ii-1) + colPtrLocal(ii)
    end do
    colPtrLocal(size(colPtrLocal)) = nnzLocal + 1

  end subroutine getColPtr


  !> Generates indexing for local rows of the CSC matrix
  subroutine getBlockRow(iNeighbour, nNeighbourSK, iAtomStart, iSparseStart, img2CentCell,&
      & colStart, colEnd, atomsInColumns, nAtomsInColumns, blockRow)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> On which column does this processor start its matrix
    integer, intent(in) :: colStart

    !> On which column does this processor end its matrix
    integer, intent(in) :: colEnd

    !> List of atoms with elements in the columns held locally
    integer, allocatable, intent(out) :: atomsInColumns(:)

    !> Count of the atoms with elements in the columns held locally
    integer, intent(out) :: nAtomsInColumns

    !> Index for starting row of blocks in nzValLocal
    integer, allocatable, intent(out) :: blockRow(:,:)

    integer, allocatable :: blockList(:,:)
    logical, allocatable :: tRowTrans(:)
    integer :: nAtom, nOrb1, nOrb2
    integer :: iAtom1, iAtom2, iAtom2f, iNeigh, iNext, iOrig, ii, jj

    nAtom = size(nNeighbourSK)

    allocate(blockRow(0 : size(iNeighbour, dim=1) - 1, nAtom))
    blockRow(:,:) = 0

    allocate(atomsInColumns(nAtom))
    atomsInColumns(nAtom) = 0
    nAtomsInColumns = 0

    allocate(blockList(nAtom, 2))
    allocate(tRowTrans(nAtom))
    blockList(:,:) = 0

    ! Offset in column belonging to transposed triangle
    blockList(:,2) = 1

    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1 + 1) - ii
      ! Offset in current column
      blockList(:,1) = 0
      tRowTrans(:) = .false.
      ! Starting index for column in DFTB+ sparse structure, because column probably already
      ! contains blocks coming from transposing previously processed elements.
      iNext = blockList(iAtom1, 2)
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)

        jj = iAtomStart(iAtom2f)
        nOrb2 = iAtomStart(iAtom2f + 1) - jj

        if (blockList(iAtom2f, 1) == 0) then
          blockList(iAtom2f, 1) = iNext
          iNext = iNext + nOrb2
        end if

        blockRow(iNeigh, iAtom1) = blockList(iAtom2f, 1)

        if (.not. isBlockInLocal(jj, jj + nOrb2 - 1, ii, ii + nOrb1 - 1, colStart, colEnd)) then
          cycle
        end if

        if (iAtom1 /= atomsInColumns(max(nAtomsInColumns,1))) then
          ! this atom is required for the local columns
          nAtomsInColumns = nAtomsInColumns + 1
          atomsInColumns(nAtomsInColumns) = iAtom1
        end if

        if (ii == jj) then
          ! on the diagonal, can not be in other triangle
          cycle
        end if

        ! Because of folding of periodic images, it can happen that the transposed block has already
        ! been registered.
        if (.not. tRowTrans(iAtom2f)) then
          blockList(iAtom2f,2) = blockList(iAtom2f,2) + nOrb1
          tRowTrans(iAtom2f) = .true.
        end if
      end do
    end do

  end subroutine getBlockRow

#:endif

end module dftbp_elsicsc
