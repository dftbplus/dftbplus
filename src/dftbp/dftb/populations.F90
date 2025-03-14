!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Calculates various types of charge populations
!> To do: extend to other populations than Mulliken
module dftbp_dftb_populations
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  use dftbp_common_environment, only : TEnvironment, globalTimers
  use dftbp_common_schedule, only : distributeRangeWithWorkload, assembleChunks
  use dftbp_dftb_hybridxc, only : THybridXcFunc
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_sparse2dense, only : unpackHS
  use dftbp_math_matrixops, only : adjointLowerTriangle
  use dftbp_type_commontypes, only : TOrbitals, TParallelKS
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_integral, only : TIntegral
#:if WITH_SCALAPACK
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip
  use dftbp_extlibs_scalapackfx, only : DLEN_, NB_, MB_, CSRC_, RSRC_, scalafx_indxl2g,&
      & scalafx_getdescriptor, scalafx_addl2g, scalafx_addg2l
  use dftbp_math_bisect, only : bisection
#:endif
  implicit none

  private
  public :: mulliken, skewMulliken, denseMullikenReal
  public :: getChargePerShell, denseBlockMulliken
  public :: getOnsitePopulation, getAtomicMultipolePopulation
  public :: denseSubtractDensityOfAtomsReal, denseSubtractDensityOfAtomsCmplxNonperiodic,&
      & denseSubtractDensityOfAtomsCmplxPeriodic, denseSubtractDensityOfAtomsCmplxPeriodicGlobal,&
      & denseSubtractDensityOfAtomsNospinRealNonperiodicReks,&
      & denseSubtractDensityOfAtomsSpinRealNonperiodicReks
#:if WITH_SCALAPACK
  public :: denseMullikenRealBlacs
  public :: denseSubtractDensityOfAtomsRealNonperiodicBlacs,&
      & denseSubtractDensityOfAtomsRealPeriodicBlacs,&
      & denseSubtractDensityOfAtomsCmplxNonperiodicBlacs
#:endif


  !> Provides an interface to calculate Mulliken populations, either dual basis atomic block,
  !> orbitally resolved or atom resolved
  interface mulliken
    module procedure mullikenPerBlock
    module procedure mullikenPerOrbital
    module procedure mullikenPerOrbital_bvKDensityMatrix
    module procedure mullikenPerAtom
  end interface mulliken


  !> Provides an interface to calculate Mulliken populations for anti-symmetric density matrices
  interface skewMulliken
    module procedure skewMullikenPerBlock
  end interface skewMulliken


contains

  !> Calculate the Mulliken population for each atom in the system, by summing the individual
  !! orbital contriutions on each atom.
  subroutine mullikenPerAtom(env, qq, over, rho, orb, iNeighbour, nNeighbourSK, img2CentCell, iPair)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> The charge per atom
    real(dp), intent(inout) :: qq(:)

    !> Overlap matrix in packed format
    real(dp), intent(in) :: over(:)

    !> Density matrix in Packed format
    real(dp), intent(in) :: rho(:)

    !> Information about the orbitals
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours of each real atom (central cell)
    integer, intent(in) :: iNeighbour(0:,:)

    !> List of neighbours for each atom, starting at 0 for itself
    integer, intent(in) :: nNeighbourSK(:)

    !> Indexing array to convert images of atoms back into their number in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    real(dp), allocatable :: qPerOrbital(:,:)
    integer :: nAtom

    nAtom = size(orb%nOrbAtom)
    @:ASSERT(size(qq) == nAtom)
    @:ASSERT(size(over) == size(rho))

    allocate(qPerOrbital(orb%mOrb, nAtom))
    qPerOrbital(:,:) = 0.0_dp

    call mullikenPerOrbital(env, qPerOrbital, over, rho, orb, iNeighbour, nNeighbourSK,&
        & img2CentCell, iPair)

    qq(:) = qq + sum(qPerOrbital, dim=1)

  end subroutine mullikenPerAtom


  !> Calculate the Mulliken population for each orbital in the system using purely real-space
  !! overlap and density matrix values. Currently, Mulliken is transformed into real space sums
  !! over one triangle of real space extended matrices.
  !!
  !! To do: Add description of algorithm to programer manual / documentation.
  subroutine mullikenPerOrbital(env, qq, over, rho, orb, iNeighbour, nNeighbourSK, img2CentCell,&
      & iPair)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> The charge per orbital
    real(dp), intent(inout) :: qq(:,:)

    !> Overlap matrix in packed format
    real(dp), intent(in) :: over(:)

    !> Density matrix in Packed format
    real(dp), intent(in) :: rho(:)

    !> Information about the orbitals
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours of each real atom (central cell)
    integer, intent(in) :: iNeighbour(0:,:)

    !> List of neighbours for each atom, starting at 0 for itself
    integer, intent(in) :: nNeighbourSK(:)

    !> Indexing array to convert images of atoms back into their number in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    integer :: iOrig
    integer :: iIter, iNeigh
    integer :: nAtom, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: sqrTmp(orb%mOrb,orb%mOrb)
    real(dp) :: mulTmp(orb%mOrb**2)
    integer, allocatable :: iterIndices(:)

    nAtom = size(orb%nOrbAtom)

    @:ASSERT(all(shape(qq) == (/orb%mOrb, nAtom/)))
    @:ASSERT(size(over) == size(rho))

    call distributeRangeWithWorkload(env, 1, nAtom, nNeighbourSK, iterIndices)

    do iIter = 1, size(iterIndices)
      iAtom1 = iterIndices(iIter)
      nOrb1 = orb%nOrbAtom(iAtom1)
      do iNeigh = 0, nNeighbourSK(iAtom1)
        sqrTmp(:,:) = 0.0_dp
        mulTmp(:) = 0.0_dp
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        iOrig = iPair(iNeigh, iAtom1) + 1
        mulTmp(1:nOrb1*nOrb2) = over(iOrig:iOrig+nOrb1*nOrb2-1) * rho(iOrig:iOrig+nOrb1*nOrb2-1)
        sqrTmp(1:nOrb2,1:nOrb1) = reshape(mulTmp(1:nOrb1*nOrb2), (/nOrb2,nOrb1/))
        qq(1:nOrb1,iAtom1) = qq(1:nOrb1,iAtom1) + sum(sqrTmp(1:nOrb2,1:nOrb1), dim=1)
        ! Add contribution to the other triangle sum, using the symmetry, but only when off diagonal
        if (iAtom1 /= iAtom2f) then
          qq(1:nOrb2,iAtom2f) = qq(1:nOrb2,iAtom2f) + sum(sqrTmp(1:nOrb2,1:nOrb1), dim=2)
        end if
      end do
    end do

    call assembleChunks(env, qq)

  end subroutine mullikenPerOrbital


  !> Calculate the Mulliken population for each orbital in the system using purely real-space
  !! overlap and BvK density matrix values. Currently Mulliken is transformed into real space sums
  !! over one triangle of real space extended matrices.
  !!
  !! To do: Add description of algorithm to programer manual / documentation.
  subroutine mullikenPerOrbital_bvKDensityMatrix(qq, over, rho, orb, iNeighbour, nNeighbourSK,&
      & img2CentCell, iPair, iSquare, iCellVec, cellVecs, hybridXc)

    !> Mulliken orbital charges on output (mOrb, nAtom, nSpin)
    real(dp), intent(out) :: qq(:,:,:)

    !> Overlap matrix in packed format
    real(dp), intent(in) :: over(:)

    !> Real-space BvK density matrix in square format for every shift
    real(dp), intent(in) :: rho(:,:,:,:,:,:)

    !> Information about the orbitals
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours of each real atom (central cell)
    integer, intent(in) :: iNeighbour(0:,:)

    !> List of neighbours for each atom, starting at 0 for itself
    integer, intent(in) :: nNeighbourSK(:)

    !> Indexing array to convert images of atoms back into their number in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Indexing array for the sparse Hamiltonian/overlap
    integer, intent(in) :: iPair(0:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Shift vector index for every interacting atom, including periodic images
    integer, intent(in) :: iCellVec(:)

    !> Vectors to neighbouring unit cells in relative coordinates
    real(dp), intent(in) :: cellVecs(:,:)

    !> Data for hybrid xc-functional calculation
    class(THybridXcFunc), intent(in) :: hybridXc

    !! Real-space overlap shift in relative coordinates
    real(dp) :: overShift(3)

    !> Real-space overlap shift folded to BvK cell and translated to density matrix indices
    integer :: bvKIndex(3)

    !! Starting point for atomic block in sparse (packed) format
    integer :: iOrig

    !! Neighbour atom index
    integer :: iNeigh

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Atom indices
    integer :: iAtom1, iAtom2, iAtom2f

    !! Number of orbitals of atoms
    integer :: nOrb1, nOrb2

    !! Start and end index of atomic block atom1-atom2 in square matrices
    integer :: iAt1SqrStart, iAt1SqrEnd, iAt2SqrStart, iAt2SqrEnd

    !! Number of spin channels and spin index
    integer :: nSpin, iSpin

    !! Square, temporary storage
    real(dp) :: sqrTmp(orb%mOrb, orb%mOrb)

    nAtom0 = size(qq, dim=2)
    nSpin = size(qq, dim=3)

    qq(:,:,:) = 0.0_dp

    do iSpin = 1, nSpin
      do iAtom1 = 1, nAtom0
        nOrb1 = orb%nOrbAtom(iAtom1)
        iAt1SqrStart = iSquare(iAtom1)
        iAt1SqrEnd = iSquare(iAtom1 + 1) - 1
        do iNeigh = 0, nNeighbourSK(iAtom1)
          sqrTmp(:,:) = 0.0_dp
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          nOrb2 = orb%nOrbAtom(iAtom2f)
          iOrig = iPair(iNeigh, iAtom1) + 1
          ! get real-space overlap shift (relative coordinates)
          overShift(:) = cellVecs(:, iCellVec(iAtom2))
          bvKIndex(:) = hybridXc%foldToBvKIndex(overShift)
          iAt2SqrStart = iSquare(iAtom2f)
          iAt2SqrEnd = iSquare(iAtom2f + 1) - 1
          sqrTmp(1:nOrb2, 1:nOrb1) = reshape(over(iOrig:iOrig+nOrb1*nOrb2-1), [nOrb2, nOrb1])&
              & * rho(iAt2SqrStart:iAt2SqrEnd, iAt1SqrStart:iAt1SqrEnd, bvKIndex(1), bvKIndex(2),&
              & bvKIndex(3), iSpin)
          qq(1:nOrb1, iAtom1, iSpin) = qq(1:nOrb1, iAtom1, iSpin)&
              & + sum(sqrTmp(1:nOrb2, 1:nOrb1), dim=1)
          ! Add contribution to the other triangle sum using symmetry, but only when off-diagonal
          if (iAtom1 /= iAtom2f) then
            qq(1:nOrb2, iAtom2f, iSpin) = qq(1:nOrb2, iAtom2f, iSpin)&
                & + sum(sqrTmp(1:nOrb2, 1:nOrb1), dim=2)
          end if
        end do
      end do
    end do

  end subroutine mullikenPerOrbital_bvKDensityMatrix


  !> Calculate the Mulliken population for each element of the dual atomic blocks in the system
  !! using purely real-space overlap and density matrix values.
  !!
  !! To do: Add description of algorithm to programer manual / documentation.
  subroutine mullikenPerBlock(env, qq, over, rho, orb, iNeighbour, nNeighbourSK, img2CentCell,&
      & iPair)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> The charge per atom block
    real(dp), intent(inout) :: qq(:,:,:)

    !> Overlap matrix in packed format
    real(dp), intent(in) :: over(:)

    !> Density matrix in Packed format
    real(dp), intent(in) :: rho(:)

    !> Information about the orbitals
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours of each real atom (central cell)
    integer, intent(in) :: iNeighbour(0:,:)

    !> List of neighbours for each atom, starting at 0 for itself
    integer, intent(in) :: nNeighbourSK(:)

    !> Indexing array to convert images of atoms back into their number in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    integer :: iOrig
    integer :: iIter, iNeigh
    integer :: nAtom, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: STmp(orb%mOrb,orb%mOrb)
    real(dp) :: rhoTmp(orb%mOrb,orb%mOrb)
    integer, allocatable :: iterIndices(:)

    nAtom = size(orb%nOrbAtom)

    @:ASSERT(all(shape(qq) == (/orb%mOrb,orb%mOrb,nAtom/)))
    @:ASSERT(size(over) == size(rho))

    call distributeRangeWithWorkload(env, 1, nAtom, nNeighbourSK, iterIndices)

    do iIter = 1, size(iterIndices)
      iAtom1 = iterIndices(iIter)
      nOrb1 = orb%nOrbAtom(iAtom1)
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        iOrig = iPair(iNeigh,iAtom1) + 1
        sTmp(1:nOrb2,1:nOrb1) = reshape(over(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
        rhoTmp(1:nOrb2,1:nOrb1) = reshape(rho(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
        qq(1:nOrb1,1:nOrb1,iAtom1) = qq(1:nOrb1,1:nOrb1,iAtom1)&
            & + 0.5_dp*( matmul(transpose(rhoTmp(1:nOrb2,1:nOrb1)), sTmp(1:nOrb2,1:nOrb1) )&
            & + matmul(transpose(sTmp(1:nOrb2,1:nOrb1)), rhoTmp(1:nOrb2,1:nOrb1) ) )
        ! Add contribution to the other triangle sum, using the symmetry but only when off diagonal
        if (iAtom1 /= iAtom2f) then
          qq(1:nOrb2,1:nOrb2,iAtom2f) = qq(1:nOrb2,1:nOrb2,iAtom2f)&
              & + 0.5_dp*( matmul( rhoTmp(1:nOrb2,1:nOrb1), transpose(sTmp(1:nOrb2,1:nOrb1)) )&
              & + matmul( sTmp(1:nOrb2,1:nOrb1), transpose(rhoTmp(1:nOrb2,1:nOrb1)) ) )
        end if
      end do
    end do

    call assembleChunks(env, qq)

  end subroutine mullikenPerBlock


  !> Calculate the Mulliken population for each element of the dual atomic block orbital in the
  !! system using purely real-space overlap and density matrix values.
  !!
  !! To do: Add description of algorithm to programer manual / documentation.
  subroutine skewMullikenPerBlock(env, qq, over, rho, orb, iNeighbour, nNeighbourSK, img2CentCell,&
      & iPair)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> The charge per atom block
    real(dp), intent(inout) :: qq(:,:,:)

    !> Overlap matrix in packed format
    real(dp), intent(in) :: over(:)

    !> Density matrix in packed format
    real(dp), intent(in) :: rho(:)

    !> Information about the orbitals
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours of each real atom (central cell)
    integer, intent(in) :: iNeighbour(0:,:)

    !> List of neighbours for each atom, starting at 0 for itself
    integer, intent(in) :: nNeighbourSK(:)

    !> Indexing array to convert images of atoms back into their number in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    integer :: iOrig
    integer :: iIter, iNeigh
    integer :: nAtom, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: STmp(orb%mOrb,orb%mOrb)
    real(dp) :: rhoTmp(orb%mOrb,orb%mOrb)
    integer, allocatable :: iterIndices(:)

    nAtom = size(orb%nOrbAtom)

    @:ASSERT(all(shape(qq) == (/orb%mOrb,orb%mOrb,nAtom/)))
    @:ASSERT(size(over) == size(rho))

    call distributeRangeWithWorkload(env, 1, nAtom, nNeighbourSK, iterIndices)

    do iIter = 1, size(iterIndices)
      iAtom1 = iterIndices(iIter)
      nOrb1 = orb%nOrbAtom(iAtom1)
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        iOrig = iPair(iNeigh,iAtom1) + 1
        sTmp(1:nOrb2,1:nOrb1) = reshape(over(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
        rhoTmp(1:nOrb2,1:nOrb1) = reshape(rho(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
        qq(1:nOrb1,1:nOrb1,iAtom1) = qq(1:nOrb1,1:nOrb1,iAtom1)&
            & - 0.5_dp*( matmul(transpose(rhoTmp(1:nOrb2,1:nOrb1)), sTmp(1:nOrb2,1:nOrb1) )&
            & - matmul(transpose(sTmp(1:nOrb2,1:nOrb1)), rhoTmp(1:nOrb2,1:nOrb1) ) )
        ! Add contribution to the other triangle sum, using the symmetry but only when off diagonal
        if (iAtom1 /= iAtom2f) then
          qq(1:nOrb2,1:nOrb2,iAtom2f) = qq(1:nOrb2,1:nOrb2,iAtom2f)&
              & + 0.5_dp*( matmul( rhoTmp(1:nOrb2,1:nOrb1), transpose(sTmp(1:nOrb2,1:nOrb1)) )&
              & - matmul( sTmp(1:nOrb2,1:nOrb1), transpose(rhoTmp(1:nOrb2,1:nOrb1)) ) )
        end if
      end do
    end do

    call assembleChunks(env, qq)

  end subroutine skewMullikenPerBlock


#:if WITH_SCALAPACK

  !> Mulliken analysis with distributed, dense, square, real-space matrices.
  subroutine denseMullikenRealBlacs(env, parallelKS, denseDesc, rhoSqr, overSqr, qq)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> The k-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Square (lower triangular) spin polarized density matrix
    real(dp), intent(in) :: rhoSqr(:,:,:)

    !> Square (lower triangular) overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> Mulliken charges on output (mOrb, nAtom, nSpin)
    real(dp), intent(out) :: qq(:,:,:)

    integer :: iKS, iS, iK, iAt, iGlob, iLocCol
    integer :: nLocalCols, nLocalKS

    ! need distributed matrix descriptors
    integer :: desc(DLEN_), nn

    nn = denseDesc%fullSize
    call scalafx_getdescriptor(env%blacs%orbitalGrid, nn, nn, env%blacs%rowBlockSize,&
        & env%blacs%columnBlockSize, desc)

    qq(:,:,:) = 0.0_dp

    nLocalCols = size(rhoSqr, dim=2)
    nLocalKS = size(rhoSqr, dim=3)

    do iKS = 1, nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iS = parallelKS%localKS(2, iKS)
      do iLocCol = 1, nLocalCols
        iGlob = scalafx_indxl2g(iLocCol, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
            & env%blacs%orbitalGrid%ncol)
        ! search atom index that corresponds to the global matrix index
        call bisection(iAt, denseDesc%iAtomStart, iGlob)
        qq(iGlob - denseDesc%iAtomStart(iAt) + 1, iAt, iS)&
            & = sum(overSqr(:, iLocCol) * rhoSqr(:, iLocCol, iKS))
      end do
    end do

    ! distribute all charges to all nodes via a global summation
    call mpifx_allreduceip(env%mpi%globalComm, qq, MPI_SUM)

  end subroutine denseMullikenRealBlacs

#:endif


  !> Mulliken analysis with dense, square, real-space matrices.
  subroutine denseMullikenReal(rhoSqr, overSqr, iSquare, qq)

    !> Square (lower triangular) spin polarized density matrix
    real(dp), intent(in) :: rhoSqr(:,:,:)

    !> Square (lower triangular) overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> Atom positions in the row/column of square matrices
    integer, intent(in) :: iSquare(:)

    !> Mulliken charges on output (mOrb, nAtom, nSpin)
    real(dp), intent(out) :: qq(:,:,:)

    real(dp) :: tmpSqr(size(qq, dim=1), size(qq, dim=1))
    integer :: iSpin, iAtom1, iAtom2, iStart1, iEnd1, iStart2, iEnd2
    integer :: nAtom, nSpin, nOrb1, nOrb2

    nAtom = size(qq, dim=2)
    nSpin = size(qq, dim=3)

    qq(:,:,:) = 0.0_dp

    do iSpin = 1, nSpin
      do iAtom1 = 1, nAtom
        iStart1 = iSquare(iAtom1)
        iEnd1 = iSquare(iAtom1 + 1) - 1
        nOrb1 = iEnd1 - iStart1 + 1
        do iAtom2 = iAtom1, nAtom
          iStart2 = iSquare(iAtom2)
          iEnd2 = iSquare(iAtom2 + 1) - 1
          nOrb2 = iEnd2 - iStart2 + 1
          tmpSqr(1:nOrb2, 1:nOrb1) = &
              & rhoSqr(iStart2:iEnd2, iStart1:iEnd1, iSpin)&
              & * overSqr(iStart2:iEnd2, iStart1:iEnd1)
          qq(1:nOrb1, iAtom1, iSpin) = qq(1:nOrb1, iAtom1, iSpin)&
              & + sum(tmpSqr(1:nOrb2, 1:nOrb1), dim=1)
          if (iAtom1 /= iAtom2) then
            qq(1:nOrb2, iAtom2, iSpin) = qq(1:nOrb2, iAtom2, iSpin)&
                & + sum(tmpSqr(1:nOrb2, 1:nOrb1), dim=2)
          end if
        end do
      end do
    end do

  end subroutine denseMullikenReal


  !> Returns pre-factor of Mulliken populations, depending on the number of spin-channels.
  pure function populationScalingFactor(nSpin) result(scale)

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> Pre-factor of Mulliken populations
    real(dp) :: scale

    ! distinguish cases: spin-unpolarized, colinear-spin
    select case(nSpin)
    case(1)
      scale = 1.0_dp
    case(2)
      scale = 0.5_dp
    end select

  end function populationScalingFactor


  !> Subtracts superposition of atomic densities from dense, square, real-space density matrix.
  !!
  !! For spin-polarized calculations, q0 is distributed equally to alpha and beta density matrices.
  !! Note: For periodic systems, q0 is normalized by the overlap that includes periodic images.
  subroutine denseSubtractDensityOfAtomsReal(q0, iSquare, overSqr, rho)

    !> Reference atom populations
    real(dp), intent(in) :: q0(:,:,:)

    !> Atom positions in the row/column of square matrices
    integer, intent(in) :: iSquare(:)

    !> Square (lower triangular) overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> Spin polarized (lower triangular) density matrix
    real(dp), intent(inout) :: rho(:,:,:)

    integer :: nAtom, iAtom, nSpin, iStart, iEnd, iOrb, iSpin
    real(dp) :: scale

    nAtom = size(iSquare) - 1
    nSpin = size(q0, dim=3)

    scale = populationScalingFactor(nSpin)

    do iSpin = 1, nSpin
      do iAtom = 1, nAtom
        iStart = iSquare(iAtom)
        iEnd = iSquare(iAtom + 1) - 1
        do iOrb = 1, iEnd - iStart + 1
          rho(iStart+iOrb-1, iStart+iOrb-1, iSpin) = rho(iStart+iOrb-1, iStart+iOrb-1, iSpin)&
              & - scale * q0(iOrb, iAtom, 1) / overSqr(iStart+iOrb-1, iStart+iOrb-1)
        end do
      end do
    end do

  end subroutine denseSubtractDensityOfAtomsReal


  !> Subtracts superposition of atomic densities from dense, square, complex-valued density matrix.
  !!
  !! Note: q0 is equally distributed to alpha and beta density matrices.
  subroutine denseSubtractDensityOfAtomsCmplxNonperiodic(q0, iSquare, rho)

    !> Reference atom populations
    real(dp), intent(in) :: q0(:,:,:)

    !> Atom positions in the row/column of square matrices
    integer, intent(in) :: iSquare(:)

    !> Spin polarized density matrix to substract q0 from
    complex(dp), intent(inout) :: rho(:,:,:)

    integer :: nAtom, iAtom, nSpin, iStart, iEnd, iOrb, iSpin

    nAtom = size(iSquare) - 1
    nSpin = size(q0, dim=3)

    do iSpin = 1, nSpin
      do iAtom = 1, nAtom
        iStart = iSquare(iAtom)
        iEnd = iSquare(iAtom + 1) - 1
        do iOrb = 1, iEnd - iStart + 1
          rho(iStart+iOrb-1, iStart+iOrb-1, iSpin) = rho(iStart+iOrb-1, iStart+iOrb-1, iSpin)&
              & - q0(iOrb, iAtom, 1)
        end do
      end do
    end do

  end subroutine denseSubtractDensityOfAtomsCmplxNonperiodic


  !> Subtracts superposition of atomic densities from dense, square, dual-space density matrix.
  !!
  !! Note: q0 is normalized by the overlap that includes periodic images and distributed equally to
  !! alpha and beta density matrices for spin-polarized calculations.
  subroutine denseSubtractDensityOfAtomsCmplxPeriodic(env, ints, denseDesc, parallelKS,&
      & neighbourList, kPoint, nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, q0, rho)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Integral container
    type(TIntegral), intent(in) :: ints

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> The k-points and spins to be handled
    type(TParallelKS), intent(in) :: parallelKS

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> The k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Reference atom populations
    real(dp), intent(in) :: q0(:,:,:)

    !> Spin polarized density matrix for all k-points/spins
    complex(dp), intent(inout) :: rho(:,:,:)

    !! Temporarily stores square overlap matrix of current k-point
    complex(dp), allocatable :: SSqrCplx(:,:)

    integer :: nAtom, iAtom, iStart, iEnd, iOrb, iSpin, nSpin, iK, iKS
    real(dp) :: scale

    nAtom = size(denseDesc%iAtomStart) - 1
    nSpin = size(q0, dim=3)

    allocate(SSqrCplx(denseDesc%fullSize, denseDesc%fullSize))

    scale = populationScalingFactor(nSpin)

    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iSpin = parallelKS%localKS(2, iKS)
      ! Get full complex, square, k-space overlap and store for later q0 substraction
      call env%globalTimer%startTimer(globalTimers%sparseToDense)
      call unpackHS(SSqrCplx, ints%overlap, kPoint(:, iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
      do iAtom = 1, nAtom
        iStart = denseDesc%iAtomStart(iAtom)
        iEnd = denseDesc%iAtomStart(iAtom + 1) - 1
        do iOrb = 1, iEnd - iStart + 1
          rho(iStart+iOrb-1, iStart+iOrb-1, iKS) = rho(iStart+iOrb-1, iStart+iOrb-1, iKS)&
              & - scale * q0(iOrb, iAtom, 1) / SSqrCplx(iStart+iOrb-1, iStart+iOrb-1)
        end do
      end do
    end do

  end subroutine denseSubtractDensityOfAtomsCmplxPeriodic


  !> Subtracts superposition of atomic densities from dense, square, dual-space density matrix.
  !! (assumes that all dP(k) are present on every MPI process)
  !!
  !! Note: q0 is normalized by the overlap that includes periodic images and distributed equally to
  !! alpha and beta density matrices for spin-polarized calculations.
  subroutine denseSubtractDensityOfAtomsCmplxPeriodicGlobal(env, ints, denseDesc, neighbourList,&
      & kPoint, iKiSToiGlobalKS, nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, q0,&
      & rho)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Integral container
    type(TIntegral), intent(in) :: ints

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> The k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Composite index for mapping iK/iS --> iGlobalKS for arrays present at every MPI rank
    integer, intent(in) :: iKiSToiGlobalKS(:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Reference atom populations
    real(dp), intent(in) :: q0(:,:,:)

    !> Spin polarized density matrix for all k-points/spins
    complex(dp), intent(inout) :: rho(:,:,:)

    !! Temporarily stores square overlap matrix of current k-point
    complex(dp), allocatable :: SSqrCplx(:,:)

    integer :: nAtom, iAtom, iStart, iEnd, iOrb, iSpin, nSpin, iK
    real(dp) :: scale

    nAtom = size(denseDesc%iAtomStart) - 1
    nSpin = size(q0, dim=3)

    allocate(SSqrCplx(denseDesc%fullSize, denseDesc%fullSize))

    scale = populationScalingFactor(nSpin)

    do iK = 1, size(kPoint, dim=2)
      ! Get full complex, square, k-space overlap and store for later q0 substraction
      call env%globalTimer%startTimer(globalTimers%sparseToDense)
      call unpackHS(SSqrCplx, ints%overlap, kPoint(:, iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
      do iSpin = 1, nSpin
        do iAtom = 1, nAtom
          iStart = denseDesc%iAtomStart(iAtom)
          iEnd = denseDesc%iAtomStart(iAtom + 1) - 1
          do iOrb = 1, iEnd - iStart + 1
            rho(iStart+iOrb-1, iStart+iOrb-1, iKiSToiGlobalKS(iK, iSpin))&
                & = rho(iStart+iOrb-1, iStart+iOrb-1, iKiSToiGlobalKS(iK, iSpin))&
              & - scale * q0(iOrb, iAtom, 1) / SSqrCplx(iStart+iOrb-1, iStart+iOrb-1)
          end do
        end do
      end do
    end do

  end subroutine denseSubtractDensityOfAtomsCmplxPeriodicGlobal


  !> Subtracts superposition of atomic densities from dense, square, real-space density matrix.
  !! (spin-restricted version)
  subroutine denseSubtractDensityOfAtomsNospinRealNonperiodicReks(q0, iSquare, rho)

    !> Reference atom populations
    real(dp), intent(in) :: q0(:,:,:)

    !> Atom positions in the row/column of square matrices
    integer, intent(in) :: iSquare(:)

    !> Spin polarized (lower triangular) density matrix
    real(dp), intent(inout) :: rho(:,:,:)

    integer :: nAtom, iAtom, nSpin, iStart, iEnd, iOrb, iSpin

    nAtom = size(iSquare) - 1
    nSpin = size(rho, dim=3)

    do iSpin = 1, nSpin
      do iAtom = 1, nAtom
        iStart = iSquare(iAtom)
        iEnd = iSquare(iAtom + 1) - 1
        do iOrb = 1, iEnd - iStart + 1
          rho(iStart+iOrb-1, iStart+iOrb-1, iSpin) = rho(iStart+iOrb-1, iStart+iOrb-1, iSpin)&
              & - q0(iOrb, iAtom, iSpin)
        end do
      end do
    end do

  end subroutine denseSubtractDensityOfAtomsNospinRealNonperiodicReks


  !> Subtracts superposition of atomic densities from dense, square, real-space density matrix.
  !! (spin-unrestricted version)
  subroutine denseSubtractDensityOfAtomsSpinRealNonperiodicReks(q0, iSquare, rho, iSpin)

    !> Reference atom populations
    real(dp), intent(in) :: q0(:,:,:)

    !> Atom positions in the row/colum of square matrix
    integer, intent(in) :: iSquare(:)

    !> Spin polarized (lower triangular) matrix
    real(dp), intent(inout) :: rho(:,:,:)

    !> Spin index
    integer, intent(in) :: iSpin

    integer :: nAtom, iAtom, iStart, iEnd, iOrb

    nAtom = size(iSquare) - 1

    do iAtom = 1, nAtom
      iStart = iSquare(iAtom)
      iEnd = iSquare(iAtom + 1) - 1
      do iOrb = 1, iEnd - iStart + 1
        rho(iStart+iOrb-1, iStart+iOrb-1, iSpin) = rho(iStart+iOrb-1, iStart+iOrb-1, iSpin)&
            & - 0.5_dp * q0(iOrb, iAtom, 1)
      end do
    end do

  end subroutine denseSubtractDensityOfAtomsSpinRealNonperiodicReks


#:if WITH_SCALAPACK

  !> Subtracts superposition of atomic densities from distributed, dense, square, real-valued
  !! density matrix.
  !!
  !! For spin-polarized calculations, q0 is distributed equally to alpha and beta density matrices.
  subroutine denseSubtractDensityOfAtomsRealNonperiodicBlacs(env, parallelKS, q0, denseDesc, rho)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> The k-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Reference atom populations
    real(dp), intent(in) :: q0(:,:,:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Spin polarized (lower triangular) matrix
    real(dp), intent(inout) :: rho(:,:,:)

    integer :: nAtom, iKS, iAt, iOrbStart, nOrb, iOrb
    real(dp) :: tmp(size(q0, dim=1), size(q0, dim=1))
    real(dp) :: scale

    nAtom = size(q0, dim=2)

    scale = populationScalingFactor(size(q0, dim=3))

    do iKS = 1, parallelKS%nLocalKS
      do iAt = 1, nAtom
        iOrbStart = denseDesc%iAtomStart(iAt)
        nOrb = denseDesc%iAtomStart(iAt + 1) - iOrbStart
        tmp(:,:) = 0.0_dp
        do iOrb = 1, nOrb
          tmp(iOrb, iOrb) = -scale * q0(iOrb, iAt, 1)
        end do
        call scalafx_addl2g(env%blacs%orbitalGrid, tmp(1:nOrb, 1:nOrb), denseDesc%blacsOrbSqr,&
            & iOrbStart, iOrbStart, rho(:,:,iKS))
      end do
    end do

  end subroutine denseSubtractDensityOfAtomsRealNonperiodicBlacs


  !> Subtracts superposition of atomic densities from distributed, dense, square, complex-valued
  !! density matrix.
  !!
  !! For spin-polarized calculations, q0 is distributed equally to alpha and beta density matrices.
  subroutine denseSubtractDensityOfAtomsCmplxNonperiodicBlacs(env, parallelKS, q0, denseDesc, rho)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> The k-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Reference atom populations
    real(dp), intent(in) :: q0(:,:,:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Spin polarized (lower triangular) matrix
    complex(dp), intent(inout) :: rho(:,:,:)

    integer :: nAtom, iKS, iAt, iOrbStart, nOrb, iOrb
    complex(dp) :: tmp(size(q0, dim=1), size(q0, dim=1))

    nAtom = size(q0, dim=2)

    do iKS = 1, parallelKS%nLocalKS
      do iAt = 1, nAtom
        iOrbStart = denseDesc%iAtomStart(iAt)
        nOrb = denseDesc%iAtomStart(iAt + 1) - iOrbStart
        tmp(:,:) = (0.0_dp, 0.0_dp)
        do iOrb = 1, nOrb
          tmp(iOrb, iOrb) = -q0(iOrb, iAt, 1)
        end do
        call scalafx_addl2g(env%blacs%orbitalGrid, tmp(1:nOrb, 1:nOrb), denseDesc%blacsOrbSqr,&
            & iOrbStart, iOrbStart, rho(:,:,iKS))
      end do
    end do

  end subroutine denseSubtractDensityOfAtomsCmplxNonperiodicBlacs


  !> Subtracts superposition of atomic densities from distributed, dense, square, real-valued
  !! density matrix.
  !!
  !! Note: q0 is normalized by the overlap that includes periodic images and distributed equally to
  !! alpha and beta density matrices for spin-polarized calculations.
  subroutine denseSubtractDensityOfAtomsRealPeriodicBlacs(env, parallelKS, q0, denseDesc, overSqr,&
      & rho)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> The k-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Reference atom populations
    real(dp), intent(in) :: q0(:,:,:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Square (lower triangular) overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> Spin polarized (lower triangular) matrix
    real(dp), intent(inout) :: rho(:,:,:)

    integer :: nAtom, iKS, iS, iAt, iOrbStart, nOrb, iOrb
    real(dp) :: tmp(size(q0, dim=1), size(q0, dim=1)), tmpSqr(size(q0, dim=1), size(q0, dim=1))
    real(dp) :: scale

    nAtom = size(q0, dim=2)

    scale = populationScalingFactor(size(q0, dim=3))

    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      do iAt = 1, nAtom
        iOrbStart = denseDesc%iAtomStart(iAt)
        nOrb = denseDesc%iAtomStart(iAt + 1) - iOrbStart
        tmpSqr(:,:) = 0.0_dp
        call scalafx_addg2l(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, iOrbStart, iOrbStart,&
            & overSqr, tmpSqr(1:nOrb, 1:nOrb))
        call mpifx_allreduceip(env%mpi%groupComm, tmpSqr, MPI_SUM)
        tmp(:,:) = 0.0_dp
        do iOrb = 1, nOrb
          tmp(iOrb, iOrb) = -scale * q0(iOrb, iAt, 1) / tmpSqr(iOrb, iOrb)
        end do
        call scalafx_addl2g(env%blacs%orbitalGrid, tmp(1:nOrb, 1:nOrb), denseDesc%blacsOrbSqr,&
            & iOrbStart, iOrbStart, rho(:,:,iKS))
      end do
    end do

  end subroutine denseSubtractDensityOfAtomsRealPeriodicBlacs

#:endif


  !> Calculate the number of charges per shell given the orbital charges.
  subroutine getChargePerShell(qOrb, orb, species, chargePerShell, qRef)

    !> Charges in each orbital, for each atom and spin channel
    real(dp), intent(in) :: qOrb(:,:,:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> Resulting charges in atomic shells
    real(dp), intent(out) :: chargePerShell(:,:,:)

    !> Optional reference values
    real(dp), intent(in), optional :: qRef(:,:,:)

    real(dp), allocatable :: qq(:,:,:)
    integer :: iAt, iSp, iSh
    integer :: nAtom, nSpin

    if (present(qRef)) then
      qq = qOrb - qRef
    else
      qq = qOrb
    end if

    nAtom = size(chargePerShell, dim=2)
    nSpin = size(chargePerShell, dim=3)
    chargePerShell(:,:,:) = 0.0_dp
    do iAt = 1, nAtom
      iSp = species(iAt)
      do iSh = 1, orb%nShell(iSp)
        chargePerShell(iSh, iAt, 1:nSpin) = chargePerShell(iSh, iAt, 1:nSpin)&
            & + sum(qq(orb%posShell(iSh, iSp) : orb%posShell(iSh + 1, iSp) - 1, iAt, 1:nSpin),&
            & dim=1)
      end do
    end do

  end subroutine getChargePerShell


  !> Calculates the on-site Mulliken population for each atom.
  !!
  !! It returns the Mulliken population stemming from the orbitals on a given atom without any
  !! contributions due to the overlap with other atoms (net atomic population).
  subroutine getOnsitePopulation(rho, orb, iPair, qq)

    !> Density matrix in Packed format
    real(dp), intent(in) :: rho(:)

    !> Atomic species information
    type(TOrbitals), intent(in) :: orb

    !> Indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Onsite population per atom
    real(dp), intent(out) :: qq(:)

    integer :: iOrig, iAt
    integer :: nAtom, nOrb

    nAtom = size(orb%nOrbAtom)
    @:ASSERT(size(qq) == nAtom)

    do iAt = 1, nAtom
      nOrb = orb%nOrbAtom(iAt)
      iOrig = iPair(0, iAt) + 1
      ! Sum up the diagonal elements of the density matrix.
      qq(iAt) = sum(rho(iOrig : iOrig + nOrb * nOrb - 1 : nOrb + 1))
    end do

  end subroutine getOnsitePopulation


  !> Block mulliken analysis with dense, square, real-valued lower triangle matrices.
  subroutine denseBlockMulliken(rhoSqr, overSqr, iSquare, qq)

    !> Square (lower triangular) spin polarized density matrix
    real(dp), intent(in) :: rhoSqr(:,:,:)

    !> Square (lower triangular) overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> Atom positions in the row/column of square matrices
    integer, intent(in) :: iSquare(:)

    !> Mulliken block charges on output (mOrb, mOrb, nAtom, nSpin)
    real(dp), intent(out) :: qq(:,:,:,:)

    real(dp), allocatable :: tmpS(:,:)
    real(dp), allocatable :: tmpD(:,:)

    integer :: nAOs, nAtom, nSpin, ii, jj, iAt, iSpin, nOrb

    nSpin = size(rhoSqr, dim=3)
    nAtom = size(iSquare, dim=1) - 1
    nAOs = size(rhoSqr, dim=1)

    allocate(tmpS(nAOs,nAOs))
    allocate(tmpD(nAOs,nAOs))

    ! Symmetrize overlap
    tmpS(:,:) = overSqr + transpose(overSqr)
    do ii = 1, nAOs
      tmpS(ii,ii) = overSqr(ii,ii)
    end do

    qq(:,:,:,:) = 0.0_dp
    do iSpin = 1, nSpin

      ! Symmetrize density matrix for spin channel
      tmpD(:,:) = rhoSqr(:,:,iSpin) + transpose(rhoSqr(:,:,iSpin))
      do ii = 1, nAOs
        tmpD(ii,ii) = rhoSqr(ii,ii,iSpin)
      end do

      do iAt = 1, nAtom
        ii = iSquare(iAt)
        jj = iSquare(iAt+1)
        nOrb = jj - ii
        qq(:nOrb,:nOrb,iAt,iSpin) = matmul(tmpS(ii:jj-1,:), tmpD(:,ii:jj-1))
        qq(:nOrb,:nOrb,iAt,iSpin) = 0.5_dp * (qq(:nOrb,:nOrb,iAt,iSpin)&
            & + transpose(qq(:nOrb,:nOrb,iAt,iSpin)))
      end do

    end do

  end subroutine denseBlockMulliken


  !> Evaluate cumulative atomic multipole moments from block sparse density matrix and multipole
  !! integrals.
  subroutine getAtomicMultipolePopulation(mpAtom, mpintBra, mpintKet, rho, orb, iNeighbour,&
      & nNeighbourSK, img2CentCell, iPair)

    !> Cumulative atomic multipole moments
    real(dp), intent(inout) :: mpAtom(:,:,:)

    !> Multipole moment integral in packed format
    real(dp), intent(in) :: mpintBra(:,:), mpintKet(:,:)

    !> Density matrix in Packed format
    real(dp), intent(in) :: rho(:,:)

    !> Information about the orbitals
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours of each real atom (central cell)
    integer, intent(in) :: iNeighbour(0:,:)

    !> List of neighbours for each atom, starting at 0 for itself
    integer, intent(in) :: nNeighbourSK(:)

    !> Indexing array to convert images of atoms back into their number in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    integer :: nAtom, iAt1, iAt2, img, ind, iNeigh, iOrb1, iOrb2, nBlk, iBlk, iSpin, nSpin
    real(dp) :: pop1(size(mpAtom, 1)), pop2(size(mpAtom, 1))

    nSpin = size(rho, dim=2)
    nAtom = size(orb%nOrbAtom)
    mpAtom(:,:,:) = 0.0_dp

    do iSpin = 1, nSpin
      do iAt1 = 1, nAtom
        do iNeigh = 0, nNeighbourSK(iAt1)
          img = iNeighbour(iNeigh, iAt1)
          iAt2 = img2CentCell(img)
          ind = iPair(iNeigh,iAt1)
          nBlk = orb%nOrbAtom(iAt2)
          pop1(:) = 0.0_dp
          pop2(:) = 0.0_dp
          do iOrb1 = 1, orb%nOrbAtom(iAt1)
            do iOrb2 = 1, nBlk
              iBlk = ind + iOrb2 + nBlk*(iOrb1-1)
              pop1(:) = pop1 + mpintKet(:, iBlk) * rho(iBlk, iSpin)
              pop2(:) = pop2 + mpintBra(:, iBlk) * rho(iBlk, iSpin)
            end do
          end do
          mpAtom(:, iAt1, iSpin) = mpAtom(:, iAt1, iSpin) + pop1
          if (iAt1 /= iAt2) then
            mpAtom(:, iAt2, iSpin) = mpAtom(:, iAt2, iSpin) + pop2
          end if
        end do
      end do
    end do

  end subroutine getAtomicMultipolePopulation


end module dftbp_dftb_populations
