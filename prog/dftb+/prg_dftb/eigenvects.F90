!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Module to wrap around the process of converting from a Hamiltonian and
!!* overlap in sparse form into eigenvectors
module eigenvects
  use assert
  use accuracy
  use eigensolver
  use sparse2dense
  use message
  use commontypes, only : TOrbitals
  use angmomentum
  implicit none

  public :: diagonalize
  private

  interface diagonalize
    module procedure realH
    module procedure cmplxH
    module procedure cmplx2Cmpnt
    module procedure cmplx2CmpntKpts
  end interface

contains

  !!* Diagonalizes a sparse represented Hamiltonian and overlap to give
  !!* the eigenvectors and values, as well as often the Cholesky factorized
  !!* overlap matrix (due to a side effect of lapack)
  !!* @param HSqrReal Large square matrix for the resulting eigenvectors
  !!* @param SSqrReal Large square matrix for the overlap workspace, often
  !!* overwritten with the  Cholesky factorized form.
  !!* @param eigen The eigenvalues of the matrices
  !!* @param ham The sparse represented Hamiltonian in real space
  !!* @param over The sparse represented overlap matrix in real space
  !!* @param iNeighbor List of atomic neighbors for each central cell atom
  !!* @param nNeighbor Number of atomic neighbors for each central cell atom
  !!* @param iAtomStart Indexing array for the large square matrices to relate
  !!* atom number to position in the matrix
  !!* @param img2CentCell Array to relate image atoms outside the central cell
  !!* to their real counterpart inside the cell
  !!* @param iSolver Choice of eigensolver, 4 lapack dense solvers currently
  !!* supported
  !!* @param jobz type of eigen-problem, either 'V'/'v' vectors or 'N'/'n'
  !!* eigenvalues only
  subroutine realH(HSqrReal, SSqrReal, eigen, ham, over, iNeighbor, nNeighbor, &
      &iAtomStart, iPair, img2CentCell, iSolver, jobz)
    real(dp), intent(out) :: HSqrReal(:,:)
    real(dp), intent(out) :: SSqrReal(:,:)
    real(dp), intent(out) :: eigen(:)
    real(dp), intent(in)  :: ham(:)
    real(dp), intent(in)  :: over(:)
    integer,  intent(in)  :: iNeighbor(0:,:)
    integer,  intent(in)  :: nNeighbor(:)
    integer,  intent(in)  :: iAtomStart(:)
    integer,  intent(in)  :: iPair(0:,:)
    integer,  intent(in)  :: img2CentCell(:)
    integer,  intent(in)  :: iSolver
    character, intent(in) :: jobz

    integer :: nOrb

    @:ASSERT(size(HSqrReal, dim=1) == size(HSqrReal, dim=2))
    @:ASSERT(all(shape(HSqrReal) == shape(SSqrReal)))
    @:ASSERT(size(HSqrReal, dim=1) == size(eigen))
    @:ASSERT(size(over) == size(ham))
    @:ASSERT(size(img2CentCell) >= maxval(iNeighbor))
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

    call unpackHS(HSqrReal,ham,iNeighbor,nNeighbor,iAtomStart,iPair, &
        &img2CentCell)
    !print *, "HSQRREAL:"
    !print *, HSqrReal
    call unpackHS(SSqrReal,over,iNeighbor,nNeighbor,iAtomStart,iPair, &
        & img2CentCell)

    select case(iSolver)
    case(1)
      call hegv(HSqrReal,SSqrReal,eigen,'L',jobz)
    case(2)
      call hegvd(HSqrReal,SSqrReal,eigen,'L',jobz)
    case(3)
      ! use subspace for all eigenstates
      nOrb = size(eigen)
      call gvr(HSqrReal,SSqrReal,eigen,'L',jobz,ilIn=1,iuIn=nOrb)
    case default
      call error('Unknown eigensolver')
    end select

  end subroutine realH



  !!* Diagonalizes a sparse represented Hamiltonian and overlap to give
  !!* the eigenvectors and values, as well as often the Cholesky factorized
  !!* overlap matrix (due to a side effect of lapack)
  !!* @param HSqrCplx Large square matrix for the resulting eigenvectors
  !!* @param SSqrCplx Large square matrix for the overlap workspace, often
  !!* overwritten with the  Cholesky factorized form.
  !!* @param eigen The eigenvalues of the matrices
  !!* @param ham The sparse represented Hamiltonian in real space
  !!* @param over The sparse represented overlap matrix in real space
  !!* @param kpoint The k-point to evaluate the phase factors for
  !!* @param iNeighbor List of atomic neighbors for each central cell atom
  !!* @param nNeighbor Number of atomic neighbors for each central cell atom
  !!* @param iCellVec Index of the cell translation vector for each atom.
  !!& @param cellVec Relative coordinates of the cell translation vectors.
  !!* @param iAtomStart Indexing array for the large square matrices to relate
  !!* atom number to position in the matrix
  !!* @param img2CentCell Array to relate image atoms outside the central cell
  !!* to their real counterpart inside the cell
  !!* @param iSolver Choice of eigensolver, 4 lapack dense solvers currently
  !!* supported
  !!* @param jobz type of eigen-problem, either 'V'/'v' vectors or 'N'/'n'
  !!* eigenvalues only
  subroutine cmplxH(HSqrCplx, SSqrCplx, eigen, ham, over, kpoint, iNeighbor, &
      &nNeighbor, iCellVec, cellVec, iAtomStart, iPair, img2CentCell, &
      &iSolver, jobz)
    complex(dp), intent(out) :: HSqrCplx(:,:)
    complex(dp), intent(out) :: SSqrCplx(:,:)
    real(dp), intent(out)    :: eigen(:)
    real(dp), intent(in)     :: ham(:)
    real(dp), intent(in)     :: over(:)
    real(dp), intent(in)     :: kpoint(3)
    integer,  intent(in)     :: iNeighbor(0:,:)
    integer,  intent(in)     :: nNeighbor(:)
    integer,  intent(in)     :: iCellVec(:)
    real(dp), intent(in)     :: cellVec(:,:)
    integer,  intent(in)     :: iAtomStart(:)
    integer,  intent(in)     :: iPair(0:,:)
    integer,  intent(in)     :: img2CentCell(:)
    integer,  intent(in)     :: iSolver
    character, intent(in)    :: jobz

    integer :: nOrb
    @:ASSERT(size(HSqrCplx, dim=1) == size(HSqrCplx, dim=2))
    @:ASSERT(all(shape(HSqrCplx) == shape(SSqrCplx)))
    @:ASSERT(size(HSqrCplx, dim=1) == size(eigen))
    @:ASSERT(size(over) == size(ham))
    @:ASSERT(size(img2CentCell) >= maxval(iNeighbor))
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

    call unpackHS(HSqrCplx, ham, kPoint, iNeighbor, nNeighbor, iCellVec, &
        &cellVec, iAtomStart, iPair, img2CentCell)
    call unpackHS(SSqrCplx, over, kPoint, iNeighbor, nNeighbor, iCellVec, &
        &cellVec, iAtomStart, iPair, img2CentCell)

    select case(iSolver)
    case(1)
      call hegv(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(2)
      call hegvd(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(3)
      ! use subspace for all eigenstates
      nOrb = size(eigen)
      call gvr(HSqrCplx,SSqrCplx,eigen,'L',jobz,ilIn=1,iuIn=nOrb)
    case default
      call error('Unknown eigensolver')
    end select

  end subroutine cmplxH

  !!* Diagonalizes a sparse represented Hamiltonian and overlap to give
  !!* the eigenvectors and values, as well as often the Cholesky factorized
  !!* overlap matrix (due to a side effect of lapack)
  !!* @param HSqrCplx Large square matrix for the resulting eigenvectors
  !!* @param SSqrCplx Large square matrix for the overlap workspace, often
  !!* overwritten with the  Cholesky factorized form.
  !!* @param eigen The eigenvalues of the matrices
  !!* @param ham The sparse represented Hamiltonian in real space
  !!* @param over The sparse represented overlap matrix in real space
  !!* @param iNeighbor List of atomic neighbors for each central cell atom
  !!* @param nNeighbor Number of atomic neighbors for each central cell atom
  !!* @param iAtomStart Indexing array for the large square matrices to relate
  !!* atom number to position in the matrix
  !!* @param img2CentCell Array to relate image atoms outside the central cell
  !!* to their real counterpart inside the cell
  !!* @param iSolver Choice of eigensolver, 4 lapack dense solvers currently
  !!* supported
  !!* @param jobz type of eigen-problem, either 'V'/'v' vectors or 'N'/'n'
  !!* eigenvalues only
  !!* @param xi optional spin orbit constants for each shell of each species
  !!* @param orb Contains information about the atomic orbitals in the system
  !!* @param species optional atomic species
  !!* @param iHam optional coefficients for imaginary part of the Hamiltonian
  subroutine cmplx2Cmpnt(HSqrCplx, SSqrCplx, eigen, ham, over, iNeighbor, &
      & nNeighbor, iAtomStart, iPair, img2CentCell, iSolver, jobz,xi,orb, &
      & species, iHam)

    complex(dp), intent(out) :: HSqrCplx(:,:)
    complex(dp), intent(out) :: SSqrCplx(:,:)
    real(dp), intent(out)    :: eigen(:)
    real(dp), intent(in)     :: ham(:,:)
    real(dp), intent(in)     :: over(:)
    integer,  intent(in)     :: iNeighbor(0:,:)
    integer,  intent(in)     :: nNeighbor(:)
    integer,  intent(in)     :: iAtomStart(:)
    integer,  intent(in)     :: iPair(0:,:)
    integer,  intent(in)     :: img2CentCell(:)
    integer,  intent(in)     :: iSolver
    character, intent(in)    :: jobz
    real(dp), intent(in), optional        :: xi(:,:)
    type(TOrbitals), intent(in), optional :: orb
    integer, intent(in), optional         :: species(:)
    real(dp), intent(in), optional        :: iHam(:,:)

    integer :: nOrb, nSpin, ii, jj, kk
    real(dp), allocatable :: work(:,:)
    logical :: tSpinOrb
    integer :: nAtom, nSpecies
    complex(dp), allocatable :: AtomZ(:,:,:)
    complex(dp), allocatable :: AtomPlus(:,:,:)
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)

    @:ASSERT(size(HSqrCplx, dim=1) == size(HSqrCplx, dim=2))
    @:ASSERT(all(shape(HSqrCplx) == shape(SSqrCplx)))
    @:ASSERT(size(HSqrCplx, dim=1) == size(eigen))
    @:ASSERT(size(over,dim=1) == size(ham,dim=1))
    @:ASSERT(size(img2CentCell) >= maxval(iNeighbor))
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(present(xi) .eqv. present(orb))
    @:ASSERT(present(xi) .eqv. present(species))
    @:ASSERT(.not.(present(xi).and.present(iHam)))

    nAtom = size(nNeighbor)
    nSpin = size(ham,dim=2)
    nOrb = size(eigen)

    tSpinOrb = present(xi)

    @:ASSERT(nSpin == 4)
    @:ASSERT(mod(nOrb,2)==0)
    nOrb = nOrb / 2
    ! for the moment, but will use S as workspace in the future
    allocate(work(nOrb,nOrb))

    SSqrCplx(:,:) = 0.0_dp
    HSqrCplx(:,:) = 0.0_dp

    work(:,:) = 0.0_dp
    call unpackHS(work,over,iNeighbor,nNeighbor,iAtomStart,iPair, &
        &img2CentCell)
    SSqrCplx(1:nOrb,1:nOrb) = work(1:nOrb,1:nOrb)
    SSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) = work(1:nOrb,1:nOrb)

    ! 1 0 charge part
    ! 0 1
    work(:,:) = 0.0_dp
    call unpackHS(work,ham(:,1),iNeighbor,nNeighbor,iAtomStart,iPair, &
        &img2CentCell)
    HSqrCplx(1:nOrb,1:nOrb) = 0.5_dp*work(1:nOrb,1:nOrb)
    HSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) = 0.5_dp*work(1:nOrb,1:nOrb)
    if (present(iHam)) then
      work(:,:) = 0.0_dp
      call unpackHS(work,iHam(:,1),iNeighbor,nNeighbor,iAtomStart,iPair, &
          &img2CentCell)
      HSqrCplx(1:nOrb,1:nOrb) = HSqrCplx(1:nOrb,1:nOrb) &
          & + 0.5_dp*cmplx(0,1,dp)*work(1:nOrb,1:nOrb)
      HSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) = &
          & HSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) &
          & + 0.5_dp*cmplx(0,1,dp)*work(1:nOrb,1:nOrb)
    end if

    ! 0 1 x part
    ! 1 0
    work(:,:) = 0.0_dp
    call unpackHS(work,ham(:,2),iNeighbor,nNeighbor,iAtomStart,iPair, &
        &img2CentCell)
    call blockSymmetrizeHS(work,iAtomStart)
    HSqrCplx(nOrb+1:2*nOrb,1:nOrb) = HSqrCplx(nOrb+1:2*nOrb,1:nOrb) &
        & + 0.5_dp * work(1:nOrb,1:nOrb)
    if (present(iHam)) then
      work(:,:) = 0.0_dp
      call unpackHS(work,iHam(:,2),iNeighbor,nNeighbor,iAtomStart,iPair, &
          &img2CentCell)
      call blockAntiSymmetrizeHS(work,iAtomStart)
      HSqrCplx(nOrb+1:2*nOrb,1:nOrb) = HSqrCplx(nOrb+1:2*nOrb,1:nOrb) &
          & + 0.5_dp *cmplx(0,1,dp)* work(1:nOrb,1:nOrb)
    end if

    ! 0 -i y part
    ! i  0
    work(:,:) = 0.0_dp
    call unpackHS(work,ham(:,3),iNeighbor,nNeighbor,iAtomStart,iPair, &
        &img2CentCell)
    call blockSymmetrizeHS(work,iAtomStart)
    HSqrCplx(nOrb+1:2*nOrb,1:nOrb) = HSqrCplx(nOrb+1:2*nOrb,1:nOrb) &
        & + cmplx(0.0,0.5,dp) * work(1:nOrb,1:nOrb)
    if (present(iHam)) then
      work(:,:) = 0.0_dp
      call unpackHS(work,iHam(:,3),iNeighbor,nNeighbor,iAtomStart,iPair, &
          &img2CentCell)
      call blockAntiSymmetrizeHS(work,iAtomStart)
      HSqrCplx(nOrb+1:2*nOrb,1:nOrb) = HSqrCplx(nOrb+1:2*nOrb,1:nOrb) &
          & - 0.5_dp * work(1:nOrb,1:nOrb)
    end if

    ! 1  0 z part
    ! 0 -1
    work(:,:) = 0.0_dp
    call unpackHS(work,ham(:,4),iNeighbor,nNeighbor,iAtomStart,iPair, &
        &img2CentCell)
    HSqrCplx(1:nOrb,1:nOrb) = HSqrCplx(1:nOrb,1:nOrb) &
        & + 0.5_dp * work(1:nOrb,1:nOrb)
    HSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) = &
        & HSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) &
        & - 0.5_dp * work(1:nOrb,1:nOrb)
    if (present(iHam)) then
      work(:,:) = 0.0_dp
      call unpackHS(work,iHam(:,4),iNeighbor,nNeighbor,iAtomStart,iPair, &
          &img2CentCell)
      HSqrCplx(1:nOrb,1:nOrb) = HSqrCplx(1:nOrb,1:nOrb) &
          & + 0.5_dp * cmplx(0,1,dp) * work(1:nOrb,1:nOrb)
      HSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) = &
          & HSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) &
          & - 0.5_dp * cmplx(0,1,dp) * work(1:nOrb,1:nOrb)
    end if

    if (tSpinOrb) then
      nSpecies = maxval(species(1:nAtom))
      allocate(AtomZ(orb%mOrb,orb%mOrb,nSpecies))
      AtomZ = 0.0_dp
      allocate(AtomPlus(orb%mOrb,orb%mOrb,nSpecies))
      AtomPlus = 0.0_dp
      allocate(Lz(orb%mOrb,orb%mOrb))
      allocate(Lplus(orb%mOrb,orb%mOrb))
      do ii = 1, nSpecies
        do jj = 1, orb%nShell(ii)
          Lz = 0.0_dp
          Lplus = 0.0_dp
          kk = orb%angShell(jj,ii)
          call loperators(Lplus(1:2*kk+1,1:2*kk+1),Lz(1:2*kk+1,1:2*kk+1),kk)
          AtomZ(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1, &
              & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,ii) &
              & = 0.5_dp*xi(jj,ii)*Lz(1:2*kk+1,1:2*kk+1)
          AtomPlus(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1, &
              & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,ii) &
              & = 0.5_dp*xi(jj,ii)*Lplus(1:2*kk+1,1:2*kk+1)
        end do
      end do
      deallocate(Lplus)
      deallocate(Lz)
      do ii = 1, nAtom
        jj = species(ii)
        HSqrCplx(iAtomStart(ii):iAtomStart(ii+1)-1, &
            & iAtomStart(ii):iAtomStart(ii+1)-1) = &
            & HSqrCplx(iAtomStart(ii):iAtomStart(ii+1)-1, &
            & iAtomStart(ii):iAtomStart(ii+1)-1) &
            & + AtomZ(1:orb%nOrbSpecies(jj),1:orb%nOrbSpecies(jj),jj)
        HSqrCplx(nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1, &
            & nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1) = &
            & HSqrCplx(nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1, &
            & nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1) &
            & - AtomZ(1:orb%nOrbSpecies(jj),1:orb%nOrbSpecies(jj),jj)
        HSqrCplx(nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1, &
            & iAtomStart(ii):iAtomStart(ii+1)-1) = &
            & HSqrCplx(nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1, &
            & iAtomStart(ii):iAtomStart(ii+1)-1) &
            & + AtomPlus(1:orb%nOrbSpecies(jj),1:orb%nOrbSpecies(jj),jj)
      end do
      deallocate(AtomZ)
      deallocate(AtomPlus)
    end if

    select case(iSolver)
    case(1)
      call hegv(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(2)
      call hegvd(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(3)
      call gvr(HSqrCplx,SSqrCplx,eigen,'L',jobz,ilIn=1,iuIn=2*nOrb)
    case default
      call error('Unknown eigensolver')
    end select

    deallocate(work)

  end subroutine cmplx2Cmpnt

  !!* Diagonalizes a sparse represented Hamiltonian and overlap to give
  !!* the eigenvectors and values, as well as often the Cholesky factorized
  !!* overlap matrix (due to a side effect of lapack)
  !!* @param HSqrCplx Large square matrix for the resulting eigenvectors
  !!* @param SSqrCplx Large square matrix for the overlap workspace, often
  !!* overwritten with the  Cholesky factorized form.
  !!* @param eigen The eigenvalues of the matrices
  !!* @param ham The sparse represented Hamiltonian in real space
  !!* @param over The sparse represented overlap matrix in real space
  !!* @param iNeighbor List of atomic neighbors for each central cell atom
  !!* @param nNeighbor Number of atomic neighbors for each central cell atom
  !!* @param iAtomStart Indexing array for the large square matrices to relate
  !!* atom number to position in the matrix
  !!* @param img2CentCell Array to relate image atoms outside the central cell
  !!* to their real counterpart inside the cell
  !!* @param iSolver Choice of eigensolver, 4 lapack dense solvers currently
  !!* supported
  !!* @param jobz type of eigen-problem, either 'V'/'v' vectors or 'N'/'n'
  !!* eigenvalues only
  !!* @param xi optional spin orbit constants for each shell of each species
  !!* @param orb Contains information about the atomic orbitals in the system
  !!* @param species optional atomic species
  !!* @param iHam optional coefficients for imaginary part of the Hamiltonian
  subroutine cmplx2CmpntKpts(HSqrCplx, SSqrCplx, eigen, ham, over, kpoint, &
      & iNeighbor, nNeighbor, iCellVec, cellVec, iAtomStart, iPair, &
      & img2CentCell, iSolver, jobz,xi,orb, species, iHam)
    complex(dp), intent(out) :: HSqrCplx(:,:)
    complex(dp), intent(out) :: SSqrCplx(:,:)
    real(dp), intent(out)    :: eigen(:)
    real(dp), intent(in)     :: ham(:,:)
    real(dp), intent(in)     :: over(:)
    real(dp), intent(in)     :: kpoint(3)
    integer,  intent(in)     :: iNeighbor(0:,:)
    integer,  intent(in)     :: nNeighbor(:)
    integer,  intent(in)     :: iCellVec(:)
    real(dp), intent(in)     :: cellVec(:,:)
    integer,  intent(in)     :: iAtomStart(:)
    integer,  intent(in)     :: iPair(0:,:)
    integer,  intent(in)     :: img2CentCell(:)
    integer,  intent(in)     :: iSolver
    character, intent(in)    :: jobz
    real(dp), intent(in), optional        :: xi(:,:)
    type(TOrbitals), intent(in), optional :: orb
    integer, intent(in), optional         :: species(:)
    real(dp), intent(in), optional        :: iHam(:,:)

    integer :: nOrb, nSpin, ii, jj, kk
    complex(dp), allocatable :: work(:,:)
    logical :: tSpinOrb
    integer :: nAtom, nSpecies
    complex(dp), allocatable :: AtomZ(:,:,:)
    complex(dp), allocatable :: AtomPlus(:,:,:)
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)

    @:ASSERT(size(HSqrCplx, dim=1) == size(HSqrCplx, dim=2))
    @:ASSERT(all(shape(HSqrCplx) == shape(SSqrCplx)))
    @:ASSERT(size(HSqrCplx, dim=1) == size(eigen))
    @:ASSERT(size(over,dim=1) == size(ham,dim=1))
    @:ASSERT(size(img2CentCell) >= maxval(iNeighbor))
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(present(xi) .eqv. present(orb))
    @:ASSERT(present(xi) .eqv. present(species))
    @:ASSERT(.not.(present(xi).and.present(iHam)))

    nAtom = size(nNeighbor)
    nSpin = size(ham,dim=2)
    nOrb = size(eigen)

    tSpinOrb = present(xi)

    @:ASSERT(nSpin == 4)
    @:ASSERT(mod(nOrb,2)==0)
    nOrb = nOrb / 2
     ! for the moment, but will use S as workspace in the future
    allocate(work(nOrb,nOrb))
    SSqrCplx(:,:) = 0.0_dp
    HSqrCplx(:,:) = 0.0_dp

    work(:,:) = 0.0_dp
    call unpackHS(work,over,kPoint, iNeighbor,nNeighbor,iCellVec, &
        & cellVec, iAtomStart,iPair, img2CentCell)
    SSqrCplx(1:nOrb,1:nOrb) = work(1:nOrb,1:nOrb)
    SSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) = work(1:nOrb,1:nOrb)

    ! 1 0 charge part
    ! 0 1
    work(:,:) = 0.0_dp
    call unpackHS(work,ham(:,1),kPoint,iNeighbor,nNeighbor,iCellVec, &
        & cellVec,iAtomStart,iPair, img2CentCell)
    HSqrCplx(1:nOrb,1:nOrb) = 0.5_dp*work(1:nOrb,1:nOrb)
    HSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) = 0.5_dp*work(1:nOrb,1:nOrb)
    if (present(iHam)) then
      work(:,:) = 0.0_dp
      call unpackHS(work,iHam(:,1),kPoint,iNeighbor,nNeighbor,iCellVec, &
        & cellVec,iAtomStart,iPair, img2CentCell)
      HSqrCplx(1:nOrb,1:nOrb) = HSqrCplx(1:nOrb,1:nOrb) &
          & + 0.5_dp*cmplx(0,1,dp)*work(1:nOrb,1:nOrb)
      HSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) = &
          & HSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) &
          & + 0.5_dp*cmplx(0,1,dp)*work(1:nOrb,1:nOrb)
    end if

    ! 0 1 x part
    ! 1 0
    work(:,:) = 0.0_dp
    call unpackHS(work,ham(:,2),kPoint,iNeighbor,nNeighbor,iCellVec, &
        & cellVec,iAtomStart,iPair, img2CentCell)
    do ii = 1, nOrb
      work(ii,ii+1:) = conjg(work(ii+1:,ii))
    end do

    HSqrCplx(nOrb+1:2*nOrb,1:nOrb) = HSqrCplx(nOrb+1:2*nOrb,1:nOrb) &
        & + 0.5_dp * work(1:nOrb,1:nOrb)
    if (present(iHam)) then
      work(:,:) = 0.0_dp
      call unpackHS(work,iHam(:,2),kPoint,iNeighbor,nNeighbor,iCellVec, &
         & cellVec,iAtomStart,iPair, img2CentCell)
      do ii = 1, nOrb
        work(ii,ii+1:) = -conjg(work(ii+1:,ii))
      end do
      HSqrCplx(nOrb+1:2*nOrb,1:nOrb) = HSqrCplx(nOrb+1:2*nOrb,1:nOrb) &
          & + 0.5_dp *cmplx(0,1,dp)* work(1:nOrb,1:nOrb)
    end if


    ! 0 -i y part
    ! i  0
    work(:,:) = 0.0_dp
    call unpackHS(work,ham(:,3),kPoint,iNeighbor,nNeighbor,iCellVec, &
        & cellVec,iAtomStart,iPair, img2CentCell)
    do ii = 1, nOrb
      work(ii,ii+1:) = conjg(work(ii+1:,ii))
    end do

    HSqrCplx(nOrb+1:2*nOrb,1:nOrb) = HSqrCplx(nOrb+1:2*nOrb,1:nOrb) &
        & + cmplx(0.0,0.5,dp) * work(1:nOrb,1:nOrb)
    if (present(iHam)) then
      work(:,:) = 0.0_dp
      call unpackHS(work,iHam(:,3),kPoint,iNeighbor,nNeighbor,iCellVec, &
          & cellVec,iAtomStart,iPair, img2CentCell)
      do ii = 1, nOrb
        work(ii,ii+1:) = -conjg(work(ii+1:,ii))
      end do
      do ii = 1, nOrb
        work(ii+1:,ii) = -conjg(work(ii,ii+1:))
      end do

      HSqrCplx(nOrb+1:2*nOrb,1:nOrb) = HSqrCplx(nOrb+1:2*nOrb,1:nOrb) &
          & - 0.5_dp * work(1:nOrb,1:nOrb)
    end if


    ! 1  0 z part
    ! 0 -1
    work(:,:) = 0.0_dp
    call unpackHS(work,ham(:,4),kPoint,iNeighbor,nNeighbor,iCellVec, &
        & cellVec,iAtomStart,iPair, img2CentCell)
    HSqrCplx(1:nOrb,1:nOrb) = HSqrCplx(1:nOrb,1:nOrb) &
        & + 0.5_dp * work(1:nOrb,1:nOrb)
    HSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) = &
        & HSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) &
        & - 0.5_dp * work(1:nOrb,1:nOrb)
    if (present(iHam)) then
      work(:,:) = 0.0_dp
      call unpackHS(work,iHam(:,4),kPoint,iNeighbor,nNeighbor,iCellVec, &
        & cellVec,iAtomStart,iPair, img2CentCell)
      HSqrCplx(1:nOrb,1:nOrb) = HSqrCplx(1:nOrb,1:nOrb) &
          & + 0.5_dp * cmplx(0,1,dp) * work(1:nOrb,1:nOrb)
      HSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) = &
          & HSqrCplx(nOrb+1:2*nOrb,nOrb+1:2*nOrb) &
          & - 0.5_dp * cmplx(0,1,dp) * work(1:nOrb,1:nOrb)
    end if

    if (tSpinOrb) then
      nSpecies = maxval(species(1:nAtom))
      allocate(AtomZ(orb%mOrb,orb%mOrb,nSpecies))
      AtomZ = 0.0_dp
      allocate(AtomPlus(orb%mOrb,orb%mOrb,nSpecies))
      AtomPlus = 0.0_dp
      allocate(Lz(orb%mOrb,orb%mOrb))
      allocate(Lplus(orb%mOrb,orb%mOrb))
      do ii = 1, nSpecies
        do jj = 1, orb%nShell(ii)
          Lz = 0.0_dp
          Lplus = 0.0_dp
          kk = orb%angShell(jj,ii)
          call loperators(Lplus(1:2*kk+1,1:2*kk+1),Lz(1:2*kk+1,1:2*kk+1),kk)
          AtomZ(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1, &
              & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,ii) &
              & = 0.5_dp*xi(jj,ii)*Lz(1:2*kk+1,1:2*kk+1)
          AtomPlus(orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1, &
              & orb%posShell(jj,ii):orb%posShell(jj+1,ii)-1,ii) &
              & = 0.5_dp*xi(jj,ii)*Lplus(1:2*kk+1,1:2*kk+1)
        end do
      end do
      deallocate(Lplus)
      deallocate(Lz)
      do ii = 1, nAtom
        jj = species(ii)
        HSqrCplx(iAtomStart(ii):iAtomStart(ii+1)-1, &
            & iAtomStart(ii):iAtomStart(ii+1)-1) = &
            & HSqrCplx(iAtomStart(ii):iAtomStart(ii+1)-1, &
            & iAtomStart(ii):iAtomStart(ii+1)-1) &
            & + AtomZ(1:orb%nOrbSpecies(jj),1:orb%nOrbSpecies(jj),jj)
        HSqrCplx(nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1, &
            & nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1) = &
            & HSqrCplx(nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1, &
            & nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1) &
            & - AtomZ(1:orb%nOrbSpecies(jj),1:orb%nOrbSpecies(jj),jj)
        HSqrCplx(nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1, &
            & iAtomStart(ii):iAtomStart(ii+1)-1) = &
            & HSqrCplx(nOrb+iAtomStart(ii):nOrb+iAtomStart(ii+1)-1, &
            & iAtomStart(ii):iAtomStart(ii+1)-1) &
            & + AtomPlus(1:orb%nOrbSpecies(jj),1:orb%nOrbSpecies(jj),jj)
      end do
      deallocate(AtomZ)
      deallocate(AtomPlus)
    end if

    select case(iSolver)
    case(1)
      call hegv(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(2)
      call hegvd(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(3)
      call gvr(HSqrCplx,SSqrCplx,eigen,'L',jobz,ilIn=1,iuIn=2*nOrb)
    case default
      call error('Unknown eigensolver')
    end select

    deallocate(work)

  end subroutine cmplx2CmpntKpts

end module eigenvects
