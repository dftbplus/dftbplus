!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to apply the sp2 method using the progress library.
!>
!> For details about the library and the method, see https://github.com/lanl/qmd-progress and
!> http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00154 respectively.
!>
module sp2progress
  use assert
  use constants
  use sparse2bml
  use periodic
  use densedescr
  use message
  use orbitals
  use commontypes, only : TParallelKS
  use dftbp_bml
  use dftbp_progress
  implicit none
  private

  public :: TSp2Solver
  public :: TSp2Solver_init


  !> Contains the internal parameters for the SP2 solver
  type :: TSp2Solver
    private

    !> Nr. of Z matrices generated (geometry step)
    integer :: iGenZ = 0

    !> Nr. of orbitals
    integer :: matrixDim

    !> Previous Z-matrices
    type(bml_matrix_t) :: zk(6)

    !> Z-matrix (congruence matrix)
    type(bml_matrix_t) :: zMat

    !> Settings for Z-matrix creation
    type(genZSPinp) :: zsp

    !> Settings for the SP2-algorithm
    type(sp2data_type) :: sp2

    !> Whether the solver had been initialized already.
    logical :: tInit = .false.

    !> Whether the Z-matrix has been initialized already.
    logical :: tInitZ = .false.
  contains
    procedure :: getDensity => TSp2Solver_getDensity
    procedure :: buildZMatrix => TSp2Solver_buildZMatrix
    final :: TSp2Solver_destruct
  end type TSp2Solver


contains

  !> Initializes an SP2 solver instance.
  subroutine TSp2Solver_init(this, matrixDim)

    !> Instance
    type(TSp2Solver), intent(out) :: this

    !> Matrix dimension (typically nr. of orbitals in the system)
    integer, intent(in) :: matrixDim

    this%matrixDim = matrixDim

    call prg_parse_sp2(this%sp2, "progress.in")
    if (this%sp2%mdim < 0) then
      this%sp2%mdim = this%matrixDim
    end if

    call prg_parse_zsp(this%zsp, "progress.in")
    if (this%zsp%mdim < 0) then
      this%zsp%mdim = this%matrixDim
    end if

    if (this%sp2%bml_type /= this%zsp%bml_type) then
      call error("SP2 and ZSP bml types differ")
    end if

    call bml_zero_matrix(this%zsp%bml_type, bml_element_real, dp, this%matrixDim, this%zsp%mdim,&
        & this%zMat)

    this%tInit = .true.

  end subroutine TSp2Solver_init


  !> Destructs an SP2 solver instance.
  subroutine TSp2Solver_destruct(this)

    !> Instance
    type(TSp2Solver), intent(inout) :: this

    integer :: ii

    if (this%tInit) then
      call bml_deallocate(this%zMat)
      if (this%tInitZ) then
        do ii = 1, size(this%zk)
          call bml_deallocate(this%zk(ii))
        end do
      end if
    end if

  end subroutine TSp2Solver_destruct


  !> Returns the density by the SP2 technique.
  !>
  !> Note: The method buildZMatrix() must be called for the given overlap matrix before this
  !> method is called.
  !>
  subroutine TSp2Solver_getDensity(this, ham, neighborList, nNeighbor, iSparseStart,&
      & img2CentCell, denseDesc, orb, nEl, parallelKS, rhoPrim, rhoSqrReal)

    !> Instance
    class(TSp2Solver), intent(inout) :: this

    !> Sparse Hamiltonian
    real(dp), intent(in) ::  ham(:,:)

    !> Neighbor list
    type(TNeighborList), intent(in) :: neighborList

    !> Nr. of neighbors for each atom.
    integer, intent(in) :: nNeighbor(:)

    !> Sparse matrix descriptor
    integer, intent(in) :: iSparseStart(:,:)

    !> Images of the atoms in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Orbital (basis) information
    type(TOrbitals), intent(in) :: orb

    !> Nr. of electrons
    real(dp), intent(in) :: nEl(:)

    !> Information about local k-points and spins
    type(TParallelKS), intent(in) :: parallelKS

    !> Sparse density matrix on exit.
    real(dp), intent(out) :: rhoPrim(:,:)

    !> Storage for the dense density matrix, if needed.
    real(dp), allocatable, intent(inout) :: rhoSqrReal(:,:,:)

    type(bml_matrix_t) :: bmlHam, bmlRho
    real(dp), allocatable :: buffer(:,:)
    real(dp) :: bandFilling, normFac
    integer :: iKS, iSpin, nSpin

    @:ASSERT(this%tInit)
    @:ASSERT(this%tInitZ)

    nSpin = size(ham, dim=2)

    ! SP2-solver needs idempotent density matrix with fillings 0 <= f <= 1.
    if (nSpin == 1) then
      normFac = 0.5_dp
    else
      normFac = 1.0_dp
    end if

    call bml_zero_matrix(this%sp2%bml_type, bml_element_real, dp, this%matrixDim, this%sp2%mdim,&
        & bmlHam)
    call bml_zero_matrix(this%sp2%bml_type, bml_element_real, dp, this%matrixDim, this%sp2%mdim,&
        & bmlRho)

    do iKS = 1, parallelKS%nLocalKS
      iSpin = parallelKS%localKS(2, iKS)
      bandFilling = nEl(iSpin) * normFac / real(this%matrixDim, dp)
      call foldToRealBml(ham(:,iSpin), neighborList%iNeighbor, nNeighbor, orb, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell, bmlHam)
      call getDensityBml(bmlHam, this%zMat, bandFilling, this%sp2, this%matrixDim, bmlRho)
      rhoPrim(:,iSpin) = 0.0_dp
      call unfoldFromRealBml(bmlRho, neighborList%iNeighbor, nNeighbor, orb, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell, rhoPrim(:,iSpin))
      if (allocated(rhoSqrReal)) then
        call bml_export_to_dense(bmlRho, buffer)
        rhoSqrReal(:,:,iSpin) = buffer
      end if
    end do

    call bml_deallocate(bmlHam)
    call bml_deallocate(bmlRho)

  end subroutine TSp2Solver_getDensity


  !> Computes the inverse overlap congruence transform.
  !>
  !> Computes the inverse overlap needed to orthogonalize the Hamiltonia before applying the SP2
  !> algorithm. This method must be called for the given overlap matrix once before getDensity()
  !> is called.
  !>
  subroutine TSp2Solver_buildZMatrix(this, over, neighborList, nNeighbor, iSparseStart,&
      & img2CentCell, denseDesc, orb)

    !> Instance.
    class(TSp2Solver), intent(inout) :: this

    !> Overlap matrix in the dftb+ compressed format.
    real(dp), intent(in) ::  over(:)

    !> Neighbor list for each atom.
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbors for each atom (incl. itself).
    integer, intent(in) :: nNeighbor(:)

    !> Indexing array for the sparse Hamiltonian.
    integer, intent(in) :: iSparseStart(:,:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    !> Dense matrix descriptor.
    type(TDenseDescr), intent(in) :: denseDesc

    !> Orbital (basis) information.
    type(TOrbitals), intent(in) :: orb

    type(bml_matrix_t) :: bmlOver

    @:ASSERT(this%tInit)

    associate (zsp => this%zsp)

      this%iGenZ = this%iGenZ + 1

      call bml_zero_matrix(zsp%bml_type, bml_element_real, dp, this%matrixDim, zsp%mdim, bmlOver)

      ! From dftb+ to bml format
      call foldToRealBml(over, neighborList%iNeighbor, nNeighbor, orb, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell, bmlOver)

      ! Congruence transformation.
      if (zsp%zsp) then
        ! Sparse Z matrix building
        call prg_buildzsparse(bmlOver, this%zMat, this%iGenZ, zsp%mdim, zsp%bml_type,&
            & this%zk(1), this%zk(2), this%zk(3), this%zk(4), this%zk(5), this%zk(6), zsp%nfirst,&
            & zsp%nrefi, zsp%nreff, zsp%numthresi, zsp%numthresf,&
            & zsp%integration, zsp%verbose)
      else
        ! Build Z matrix using diagonalization (usual method).
        call prg_buildzdiag(bmlOver, this%zMat, zsp%numthresf, zsp%mdim, zsp%bml_type)
      end if

      call bml_deallocate(bmlOver)

    end associate

    this%tInitZ = .true.

  end subroutine TSp2Solver_buildZMatrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the density in BML-format for a Hamiltonian in BML-format
  subroutine getDensityBml(bmlHam, zMat, bandFilling, sp2, matrixDim, bmlRho)

    !> Hamiltonian
    type(bml_matrix_t), intent(inout) :: bmlHam

    !> Z-matrix for the congruence transformation of the overlap
    type(bml_matrix_t), intent(inout) :: zMat

    !> Band filling factor
    real(dp), intent(in) :: bandFilling

    !> Settings for the SP2-algorithm
    type(sp2data_type), intent(in) :: sp2

    !> Dimension of the matrices
    integer, intent(in) :: matrixDim

    !> Density matrix
    type(bml_matrix_t), intent(out) :: bmlRho

    type(bml_matrix_t) :: bmlOrthoH, bmlOrthoRho
    character(100) :: errorStr

    call bml_zero_matrix(sp2%bml_type, bml_element_real, dp, matrixDim, sp2%mdim, bmlOrthoH)
    call bml_zero_matrix(sp2%bml_type, bml_element_real, dp, matrixDim, sp2%mdim, bmlOrthoRho)

    call prg_orthogonalize(bmlHam, zMat, bmlOrthoH, sp2%threshold, sp2%bml_type, sp2%verbose)

    ! Perform SP2 from progress
    if (sp2%flavor == "Basic") then
      call prg_sp2_basic(bmlOrthoH, bmlOrthoRho, sp2%threshold, bandFilling, sp2%minsp2iter,&
          & sp2%maxsp2iter, sp2%sp2conv, sp2%sp2tol, sp2%verbose)
    else if (sp2%flavor == "Alg1") then
      call prg_sp2_alg1(bmlOrthoH, bmlOrthoRho, sp2%threshold, bandFilling, sp2%minsp2iter,&
          & sp2%maxsp2iter, sp2%sp2conv, sp2%sp2tol, sp2%verbose)
    else if (sp2%flavor == "Alg2") then
      call prg_sp2_alg2(bmlOrthoH, bmlOrthoRho, sp2%threshold, bandFilling, sp2%minsp2iter,&
          & sp2%maxsp2iter, sp2%sp2conv, sp2%sp2tol, sp2%verbose)
    else
      write(errorStr, "(3A)") "Invalid SP2 flavor '", sp2%flavor, "'"
      call error(errorStr)
    end if

    call prg_deorthogonalize(bmlOrthoRho, zMat, bmlRho, sp2%threshold, sp2%bml_type, sp2%verbose)

    call bml_deallocate(bmlOrthoH)
    call bml_deallocate(bmlOrthoRho)

  end subroutine getDensityBml


end module sp2progress
