!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Module for evaluating dielectric functions
module dftbp_derivs_dielectric
  use dftbp_common_accuracy, only : dp, rsp
  use dftbp_common_constants, only : Hartree__eV, eV__Hartree, imag, pi, alpha_fs
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_file, only : TFileDescr, openFile, closeFile
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_status, only : TStatus
  use dftbp_derivs_fillings, only : filledOrEmpty, nIndependentHam, maximumFillings
  use dftbp_derivs_momentumelements, only : momentumElements
  use dftbp_dftb_densitymatrix, only : TDensityMatrix
  use dftbp_dftb_hybridxc, only : THybridXcFunc
  use dftbp_dftb_periodic, only : TAuxNeighbourList, TNeighbourList
  use dftbp_dftbplus_outputfiles, only : autotestTag, resultsTag
  use dftbp_io_taggedoutput, only : TTaggedWriter, tagLabels
  use dftbp_math_kramerskronig, only : kki2r
  use dftbp_math_simplealgebra, only : determinant33
  use dftbp_math_sorting, only : merge_sort
  use dftbp_type_commontypes, only : TOrbitals, TParallelKS
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_integral, only : TIntegral
  use dftbp_type_wrappedintr, only : TWrappedInt2, TWrappedReal1

  implicit none

  private
  public :: approxAtomDipole, dielectric, TDielectricSettings

  !> Input for dielectric calculations
  type :: TDielectricSettings

    !> Imaginary broadening factor for Kramer-Kronig and low frequency transition limit
    real(dp) :: eta = 1.0E-6_dp

    !> Grid spacing for the dielectric function
    real(dp) :: gridSpacing = 0.01_dp * eV__Hartree

    !> (Gaussian) broadening factor for the transitions
    real(dp) :: sigma = 0.1_dp * eV__Hartree

    !> Should atomic dipole, either approximated or calculated, be included
    logical :: isAtomicDipoleIncluded = .false.

  end type TDielectricSettings


  !> Ordering for symmetric 3x3 matrix elements
  integer, parameter :: quadOrder(2,6) = reshape([1,1,2,1,2,2,1,3,2,3,3,3], [2,6])
  character(1) :: labels(3) = ["x","y","z"]

contains


  !> Calculate dielectric function
  subroutine dielectric(env, settings, parallelKS, eigvals, filling, eigVecsCplx, ints,&
      & neighbourList, nNeighbourSK, symNeighbourList, nNeighbourCamSym, orb, denseDesc,&
      & iSparseStart, img2CentCell, kPoint, kWeight, rCellVecs, cellVec, iCellVec, latVecs,&
      & densityMatrix, hybridXc, taggedWriter, isAutotestWritten, isResultsTagWritten, dab,&
      & errStatus)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Calculations settings
    type(TDielectricSettings), intent(in) :: settings

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Occupations of orbitals
    real(dp), intent(in) :: filling(:,:,:)

    !> ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Integral container
    type(TIntegral), intent(in) :: ints

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> List of neighbouring atoms (symmetric version)
    type(TAuxNeighbourList), intent(in), allocatable :: symNeighbourList

    !> Symmetric neighbour list version of nNeighbourCam
    integer, intent(in), allocatable :: nNeighbourCamSym(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors to units cells in absolute units
    real(dp), intent(in) :: rCellVecs(:,:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Unit cell vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> Holds real and complex delta density matrices and pointers
    type(TDensityMatrix), intent(in) :: densityMatrix

    !> Hybrid functional, if in use
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Should the autotest tag file be writted to disc?
    logical, intent(in) :: isAutotestWritten

    !> Should the results tag file be writted to disc?
    logical, intent(in) :: isResultsTagWritten

    !> On-site atomic dipole matrix elements to second order in overlap
    real(dp), allocatable, intent(in) :: dab(:,:,:,:)

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

    integer :: iK, iKS, iS, nOrbs, nSpin, nKpts, nIndepHam
    integer :: iTrans, nLocCol, nLocRow, iFil, iEmp, ii, jj, kk
    real(dp) :: maxFill, cellVol, chiTmp(6)
    integer, allocatable :: nFilled(:,:), nEmpty(:,:)
    integer, allocatable :: nTrans(:), win(:)
    type(TWrappedInt2), allocatable :: getIA(:)
    type(TWrappedReal1), allocatable :: wij(:), fij(:)
    type(TFileDescr) :: fd

    real(dp) :: norm, lambda, widthOfNonZero
    integer :: nEgyPts, nNonZero
    real(dp) :: lowerEgy, upperEgy

    !! Small complex broadening
    complex(dp) :: eta

    real(dp), allocatable :: imChi(:,:), reChi(:,:)
    complex(dp), allocatable :: pElement(:,:,:)

    if (settings%sigma < epsilon(0.0_dp)) then
      @:RAISE_ERROR(errStatus, -1, "dielectric module: Broadening sigma too small")
    end if
    if (settings%gridSpacing < epsilon(0.0_dp)) then
      @:RAISE_ERROR(errStatus, -1, "dielectric module: GridSpacing too small")
    end if
    if (settings%eta < epsilon(0.0_dp)) then
      @:RAISE_ERROR(errStatus, -1, "dielectric module: Eta value too small")
    end if

    ! Broadening parameters
    norm = 1.0_dp / (settings%sigma * sqrt(2.0_dp * pi))
    lambda = -0.5 / (settings%sigma * settings%sigma)
    widthOfNonZero = 3.5_dp * settings%sigma
    nNonZero = nint(widthOfNonZero / settings%gridSpacing)
    eta = imag * settings%eta

    cellVol = abs(determinant33(latVecs))

    nOrbs = size(filling,dim=1)
    nKpts = size(filling,dim=2)
    nSpin = size(filling,dim=3)
    nLocCol = size(eigVecsCplx, dim=1)
    nLocRow = size(eigVecsCplx, dim=2)

    nIndepHam = nIndependentHam(nSpin)
    maxFill = maximumFillings(nSpin)
    call filledOrEmpty(nFilled, nEmpty, nIndepHam, nKpts, filling, nOrbs, maxFill)

    ! Count total number of transitions for locally stored k-points and spins
    allocate(nTrans(parallelKS%nLocalKS), source=0)
    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iS = parallelKS%localKS(2, iKS)
      nTrans(iKS) = nFilled(iS, iK) * (nOrbs - nEmpty(iS, iK) + 1)
    end do

    allocate(getIA(parallelKS%nLocalKS))
    allocate(wij(parallelKS%nLocalKS))
    allocate(fij(parallelKS%nLocalKS))

    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iS = parallelKS%localKS(2, iKS)
      allocate(getIA(iKS)%data(nTrans(iKS), 2))
      allocate(wij(iKS)%data(nTrans(iKS)))
      allocate(fij(iKS)%data(nTrans(iKS)))
      iTrans = 0
      do iFil = 1, nFilled(iS, iK)
        do iEmp = nEmpty(iS, iK), nOrbs
          iTrans = iTrans + 1
          getIA(iKS)%data(iTrans, :) = [iFil, iEmp]
          wij(iKS)%data(iTrans) = eigvals(iEmp, iK, iS) - eigvals(iFil, iK, iS)
          fij(iKS)%data(iTrans) = filling(iFil, iK, iS) - filling(iEmp, iK, iS)
        end do
      end do
      allocate(win(nTrans(iKS)))
      call merge_sort(win, wij(iKS)%data, 1.0E-3_dp*epsilon(0.0_rsp))
      getIA(iKS)%data(:,:) = getIA(iKS)%data(win, :)
      wij(iKS)%data(:) = wij(iKS)%data(win)
      fij(iKS)%data(:) = fij(iKS)%data(win)
      deallocate(win)
    end do

    lowerEgy = huge(1.0_dp)
    upperEgy = -huge(1.0_dp)
    do iKS = 1, parallelKS%nLocalKS
      lowerEgy = min(lowerEgy, wij(iKS)%data(1) - widthOfNonZero)
      upperEgy = max(upperEgy, wij(iKS)%data(nTrans(iKS)) + widthOfNonZero)
    end do
    lowerEgy = min(0.0_dp, lowerEgy)
    nEgyPts = ceiling((upperEgy-lowerEgy) / settings%gridSpacing)
    nEgyPts = nEgyPts + 2 * nNonZero

  #:if WITH_SCALAPACK

    @:ASSERT(.false.)

  #:endif

    call momentumElements(env, pElement, parallelKS, eigvals, eigVecsCplx, ints, nTrans, getIA,&
        & wij, neighbourList, nNeighbourSK, symNeighbourList, nNeighbourCamSym, orb, denseDesc,&
        & iSparseStart, img2CentCell, kPoint, kWeight, rCellVecs, cellVec, iCellVec, densityMatrix,&
        & hybridXc, dab, errStatus)

    allocate(imChi(6, 0:nEgyPts), source = 0.0_dp)

    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iS = parallelKS%localKS(2, iKS)

      do iTrans = 1, nTrans(iKS)

        do kk = 1, 6
          ii = quadOrder(1, kk)
          jj = quadOrder(2, kk)
          chiTmp(kk) = real(pElement(ii, iTrans, iKS),dp) * real(pElement(jj, iTrans, iKS),dp)&
              & + aimag(pElement(ii, iTrans, iKS)) * aimag(pElement(jj, iTrans, iKS))
        end do

        chiTmp(:) = 2.0_dp * pi * fij(iKS)%data(iTrans) * kWeight(iK) * chiTmp
        chiTmp(:) = real(chiTmp / (cellVol * (wij(iKS)%data(iTrans)**2 + eta)), dp)
        jj = ceiling((wij(iKS)%data(iTrans) - lowerEgy) / settings%gridSpacing)
        do ii = -nNonZero, nNonZero
          imChi(:, jj+ii) = imChi(:, jj+ii) + chiTmp*norm*exp(lambda * (ii*settings%gridSpacing)**2)
        end do

      end do

    end do

    deallocate(pElement)
    deallocate(getIA)
    deallocate(wij)
    deallocate(fij)

    allocate(reChi(6, 0:nEgyPts), source = 0.0_dp)
    call kki2r(settings%gridSpacing, lowerEgy, imChi, reChi, eta)

    ! Convert chi to epsilon
    reChi(:, :) = 1.0_dp + 4.0_dp * pi * reChi
    imChi(:, :) = 4.0_dp * pi * imChi

    call openFile(fd, "epsilon.dat", mode="w")

    do iK = 1, 2
      do kk = 1, 6
        if (iK == 1 .and. kk == 1) then
          write(fd%unit, "(A)", advance="NO")"#    omega"
          write(fd%unit, "(A)", advance="NO")"  "
        else
          write(fd%unit, "(A)", advance="NO")"    "
        end if
        if (iK == 1) then
          write(fd%unit, "(A)", advance="NO")"Re("
        else
          write(fd%unit, "(A)", advance="NO")"Im("
        end if
        ii = quadOrder(1, kk)
        jj = quadOrder(2, kk)
        write(fd%unit, "(A)", advance="NO")"eps_"//labels(ii)//labels(jj)//")"
      end do
    end do
    write(fd%unit,*)

    do ii = 1, nEgyPts
      write(fd%unit, "(F10.3, 12E14.6)")Hartree__eV * (lowerEgy + real(ii-1,dp) *&
          & settings%gridSpacing), reChi(:, ii), imChi(:, ii)
    end do
    call closeFile(fd)

    if (isAutotestWritten) then
      call openFile(fd, autotestTag, mode="a")
      call taggedWriter%write(fd%unit, tagLabels%dielectric, reChi(:,1) + imag * imChi(:,1))
      call taggedWriter%write(fd%unit, tagLabels%dielectricGrid, [lowerEgy, upperEgy,&
          & settings%gridSpacing])
      call closeFile(fd)
    end if
    if (isResultsTagWritten) then
      call openFile(fd, resultsTag, mode="a")
      call taggedWriter%write(fd%unit, tagLabels%dielectric, reChi + imag * imChi)
      call taggedWriter%write(fd%unit, tagLabels%dielectricGrid, [lowerEgy, upperEgy,&
          & settings%gridSpacing])
      call closeFile(fd)
    end if

  end subroutine dielectric


  !> Evaluate onsite dipole matrix elements using the approximation in the PRB of Sandu, doi:
  !! 10.1103/PhysvB.72.125105
  subroutine approxAtomDipole(over, nNeighbour, iNeighbour, iSparseStart, img2CentCell, orb,&
      & species, coord, dab)

    !> Overlap matrix
    real(dp), intent(in) :: over(:)

    !> Number of neighbours for overlap for each of the central cell atoms
    integer, intent(in) :: nNeighbour(:)

    !> Neighbour index
    integer, intent(in) :: iNeighbour(0:,:)

    !> Index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(0:,:)

    !> Image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species-list of atoms
    integer, intent(in) :: species(:)

    !> List of all atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> On-site atomic dipole matrix elements to second order in overlap
    real(dp), intent(out) :: dab(:,:,:,:)

    real(dp) :: tmpS21(orb%mOrb,orb%mOrb), tmpSdotS(orb%mOrb,orb%mOrb)
    integer :: nAtom0, iAt1, iAt2, iAt2f, iNeigh
    integer :: iSp1, iSp2, nOrb1, nOrb2, iOrig, iCart

    nAtom0 = size(nNeighbour)
    @:ASSERT(all(shape(dab) == [orb%mOrb, orb%mOrb, nAtom0, 3]))
    dab(:,:,:,:) = 0.0_dp

    do iAt1 = 1, nAtom0
      iSp1 = species(iAt1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      do iNeigh = 0, nNeighbour(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        nOrb2 = orb%nOrbSpecies(iSp2)
        iOrig = iSparseStart(iNeigh, iAt1) + 1

        tmpS21(:nOrb2, :nOrb1) = reshape(over(iOrig:iOrig+nOrb2*nOrb1-1), [nOrb2, nOrb1])

        tmpSdotS(:nOrb1,:nOrb1) = matmul(transpose(tmpS21(:nOrb2, :nOrb1)), tmpS21(:nOrb2, :nOrb1))
        do iCart = 1, 3
          dab(:nOrb1, :nOrb1, iAt1, iCart) = dab(:nOrb1, :nOrb1, iAt1, iCart)&
              & +(coord(iCart, iAt2) - coord(iCart, iAt1)) * tmpSdotS(:nOrb1,:nOrb1)
        end do

        tmpSdotS(:nOrb2,:nOrb2) = matmul(tmpS21(:nOrb2, :nOrb1), transpose(tmpS21(:nOrb2, :nOrb1)))
        do iCart = 1, 3
          dab(:nOrb2, :nOrb2, iAt2f, iCart) = dab(:nOrb2, :nOrb2, iAt2f, iCart)&
              & +(coord(iCart, iAt1) - coord(iCart, iAt2)) * tmpSdotS(:nOrb2,:nOrb2)
        end do

      end do
    end do

    dab(:,:,:,:) = 0.25_dp * dab

  end subroutine approxAtomDipole


end module dftbp_derivs_dielectric
