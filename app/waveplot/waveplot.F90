!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Program for plotting molecular orbitals as cube files.
program waveplot

  use dftbp_common_accuracy, only : dp
  use dftbp_common_globalenv, only : stdOut
  use dftbp_dftb_periodic, only : getCellTranslations, foldCoordToUnitCell
  use dftbp_io_charmanip, only : i2c
  use dftbp_io_fileid, only : getFileId
  use dftbp_math_simplealgebra, only : invert33
  use dftbp_type_linkedlist, only : TListInt, TListRealR1, len, init, append, asArray
  use dftbp_type_typegeometry, only : TGeometry

  use waveplot_grids, only : TGrid, TGrid_init, TGridData, TGridData_init, TRealTessY, &
      &TRealTessY_init, subgridsToGlobalGrid
  use waveplot_initwaveplot, only : TProgramVariables, TProgramVariables_init

#:if WITH_MPI
  use mpi, only : MPI_THREAD_FUNNELED
  use dftbp_common_mpienv, only : TMpiEnv, TMpiEnv_init
  use dftbp_extlibs_mpifx, only : mpifx_init_thread, mpifx_finalize
#:endif

  implicit none

  !> Instance containing all nessesary data
  type(TProgramVariables) :: wp

  !> Grid instances
  type(TGrid) :: totGrid
  type(TGrid), allocatable :: speciesGrids(:)

  !> Grid data instances
  type(TGridData) :: totGridDat, atomicGridDat
  type(TGridData), allocatable :: speciesGridsDat(:)

  !> Representation of real tesseral spherical harmonics
  type(TRealTessY), allocatable :: speciesRty(:)

  character(len=80) :: comments(2), fileName

  !> Variables needed for atom folding
  real(dp) :: invBoxVecs(3,3), recVecs2p(3,3)
  real(dp) :: cellMiddle(3), boxMiddle(3), frac(3), cubeCorner(3), coord(3)
  real(dp) :: mDist, dist
  real(dp), allocatable :: cellVec(:,:), rCellVec(:,:)
  integer :: ii, i1, i2, i3, iCell, iAtom
  type(TListRealR1) :: coordList
  type(TListInt) :: speciesList

  !> Variables needed in the loop over all states
  integer :: iLevel, iKPoint, iSpin
  integer :: iSpecies, iAng, ind

  !> If current level should be plotted
  logical :: tPlotLevel

  !> Arrays holding the volumetric grid data
  real(dp), allocatable :: totGridsDat(:,:,:,:)
  complex(dp), allocatable :: totGridsDatCplx(:,:,:,:,:,:)
  real(dp), allocatable, target ::  totChrg(:,:,:), atomicChrg(:,:,:), speciesChrg(:,:,:,:)
  real(dp), allocatable :: atomDensity(:,:,:)
  real(dp), allocatable, target :: totData(:,:,:), buffer(:,:,:), spinUp(:,:,:), spinDown(:,:,:)
  real(dp), allocatable, target :: copyBuffers(:,:,:,:), copyBuffersCplx(:,:,:,:)
  real(dp) :: sumTotChrg, sumAtomicChrg, sumChrg

  !> Pointers to arrays holding the volumetric grid data
  real(dp), pointer ::  pAtomicChrg(:,:,:), pSpeciesChrg(:,:,:)
  real(dp), pointer :: pTotData(:,:,:), pBuffer(:,:,:), pCopyBuffers(:,:,:,:)

  !> Auxillary variables
  integer :: jj, iL, iM
  integer :: maxAng
  integer, allocatable :: iKPointPrime(:), iLPrime(:), requiredKPoints(:), requiredLevels(:)
  integer, allocatable :: iSpinPrime(:), requiredSpins(:)
  integer :: kIndex, lIndex, sIndex, kPointCounter, levelCounter, spinCounter, KNum, LNum, SNum
  integer :: tmparray(1)

#:if WITH_MPI
  !> MPI environment, if compiled with mpifort
  type(TMpiEnv) :: mpiEnv

  ! As this is serial code, trap for run time execution
  ! on more than 1 processor with an mpi enabled build
  call mpifx_init_thread(requiredThreading=MPI_THREAD_FUNNELED)
  call TMpiEnv_init(mpiEnv)
  call mpiEnv%mpiSerialEnv()
#:endif

  ! Allocate resources
  call TProgramVariables_init(wp)
  write(stdout, "(/,A,/)") "Starting main program"

  ! Allocating buffer for general grid, total charge and spin up
  allocate(buffer(wp%option%nTotPoints(1), wp%option%nTotPoints(2), wp%option%nTotPoints(3)))
  buffer(:,:,:) = 0.0_dp
  pBuffer => buffer
  allocate(totData(wp%option%nTotPoints(1), wp%option%nTotPoints(2), wp%option%nTotPoints(3)))
  totData(:,:,:) = 0.0_dp
  pTotData => totData
  allocate(copyBuffers(wp%option%nTotPoints(1), wp%option%nTotPoints(2), wp%option%nTotPoints(3),&
      & size(wp%option%levelIndex, dim=2)))
  copyBuffers(:,:,:,:) = 0.0_dp
  pCopyBuffers => copyBuffers
  allocate(copyBuffersCplx(wp%option%nTotPoints(1), wp%option%nTotPoints(2),&
      & wp%option%nTotPoints(3), size(wp%option%levelIndex, dim=2)))
  copyBuffersCplx(:,:,:,:) = 0.0_dp

  if (wp%option%tCalcTotChrg) then
    allocate(totChrg(wp%option%nTotPoints(1), wp%option%nTotPoints(2), wp%option%nTotPoints(3)))
    totChrg(:,:,:) = 0.0_dp
    if (wp%option%tPlotTotSpin) then
      allocate(spinUp(wp%option%nTotPoints(1), wp%option%nTotPoints(2), wp%option%nTotPoints(3)))
      spinUp(:,:,:) = 0.0_dp
      allocate(spinDown(wp%option%nTotPoints(1), wp%option%nTotPoints(2), wp%option%nTotPoints(3)))
      spinDown(:,:,:) = 0.0_dp
    end if
  end if

  write (comments(1), "(A)") "Cube file generated by WAVEPLOT from data created by DFTB+"

  write(stdout, "(A)") "Origin"
  write(stdout, "(2X,3(F0.5,1X))") wp%option%origin(:)
  write(stdout, "(A)") "Total-Grid Origin"
  write(stdout, "(2X,3(F0.5,1X))") wp%option%totGridOrig
  write(stdout, "(A)") "Box"

  do jj = 1, 3
    write(stdout, "(2X,3(F0.5,1X))") wp%option%boxVecs(:, jj)
  end do

  write(stdout, "(A)") "Spatial resolution [1/Bohr]:"
  write(stdout, "(2X,3(F0.5,1X))") 1.0_dp / sqrt(sum(wp%internal%totGridVec**2, dim=1))
  write(stdout, *)

  ! Initialise total grid and volumetric data
  call TGrid_init(totGrid, wp%option%totGridOrig, wp%internal%totGridVec, wp%option%nTotPoints)
  call TGridData_init(totGridDat, totGrid, pBuffer, rwTabulationType=wp%option%rwTabulationType)

  ! Initialise atomic volumetric data
  allocate(speciesGrids(wp%xml%geo%nSpecies))
  allocate(speciesGridsDat(maxval(wp%internal%orbitalToSpecies)))
  allocate(speciesRty(maxval(wp%internal%orbitalToSpecies)))
  allocate(speciesChrg(wp%option%nSpPoints(1), wp%option%nSpPoints(2), wp%option%nSpPoints(3),&
      & maxval(wp%internal%orbitalToSpecies)))
  speciesChrg(:,:,:,:) = 0.0_dp

  ind = 1

  ! Initialise atomic grids and pretabulate orbital-species grids
  do iSpecies = 1, wp%xml%geo%nSpecies
    maxAng = maxval(wp%basis%basis(iSpecies)%angMoms)
    call TGrid_init(speciesGrids(iSpecies), wp%internal%speciesGridsOrigs(:, iSpecies),&
        & wp%internal%speciesGridsVecs(:, :, iSpecies), wp%option%nSpPoints)
    call TRealTessY_init(speciesRty(iSpecies), speciesGrids(iSpecies), maxAng)
    do iAng = wp%internal%molorb%iStos(iSpecies), wp%internal%molorb%iStos(iSpecies + 1) - 1
      iL = wp%internal%molorb%angMoms(iAng)
      do iM = - iL, iL
        pSpeciesChrg => speciesChrg(:,:,:, ind)
        call TGridData_init(speciesGridsDat(ind), speciesGrids(iSpecies), pSpeciesChrg,&
            & rwTabulationType=wp%option%rwTabulationType)
        call speciesGridsDat(ind)%tabulateBasis(speciesRty(iSpecies), iL, iM, &
            & wp%internal%molorb%rwfs(iAng)%rwf, wp%internal%molorb%rwfs(iAng)%exps, &
            & wp%internal%molorb%rwfs(iAng)%aa)
        ind = ind + 1
      end do
    end do
  end do

  ! Create density superposition of the atomic orbitals. Occupation is
  ! distributed equally on orbitals with the same angular momentum.
  if (wp%option%tCalcAtomDens) then

    allocate(atomicChrg(wp%option%nTotPoints(1), wp%option%nTotPoints(2), wp%option%nTotPoints(3)))
    atomicChrg(:,:,:) = 0.0_dp
    pAtomicChrg => atomicChrg

    ! Initialise volumetric data
    call TGridData_init(atomicGridDat, totGrid, pAtomicChrg,&
        & rwTabulationType=wp%option%rwTabulationType)

    ! Calculate total charge of atomic densities.
    call subgridsToGlobalGrid(atomicGridDat, speciesGridsDat, atomDensity, &
        & wp%internal%molorb%coords, wp%internal%orbitalOcc(:, 1), wp%internal%orbitalToAtom, &
        & wp%internal%orbitalToSpecies, wp%option%parallelRegionNum, .true., wp%option%gridInterType)

    sumAtomicChrg = sum(atomDensity) * wp%internal%gridVol

    if (wp%option%tVerbose) then
      write(stdout, "('Total charge of atomic densities:',F12.6,/)") sumAtomicChrg
    end if

    ! Plot total atomic density
    if (wp%option%tPlotAtomDens) then
      write (comments(2), "('Calc-Id:', I11, ', atomdens')") wp%xml%identity
      fileName = "wp-atomdens.cube"
      call writeCubeFile(wp%xml%geo, wp%atomicNumber%atomicNumbers, wp%internal%totGridVec,&
          & wp%option%totGridOrig, atomDensity, fileName, comments, wp%option%repeatBox)
      write(stdout, "(A)") "File '" // trim(fileName) // "' written"
    end if

  end if

  if (wp%option%tVerbose) then
    write(stdout, "(/,A5,' ',A6,' ',A6,' ',A7,' ',A11,' ',A11)") "Spin", "KPoint", "State",&
        & "Action", "Norm", "W. Occup."
  end if

  ! Fold in coordinates
  if (wp%option%tFoldCoords) then
    call invert33(invBoxVecs, wp%option%boxVecs)
    call invert33(recVecs2p, wp%xml%geo%latVecs)
    recVecs2p = reshape(recVecs2p, [3, 3], order=[2, 1])
    call foldCoordToUnitCell(wp%xml%geo%coords(:,:), wp%xml%geo%latVecs, recVecs2p)
  end if

  ! Fill the box with atoms
  if (wp%option%tFillBox) then

    ! Shift plotted region by integer lattice vectors, to have its center as close to the center
    ! of the lattice unit cell as possible.
    cellMiddle(:) = 0.5_dp * sum(wp%xml%geo%latVecs, dim=2)
    boxMiddle(:) = wp%option%origin(:) + 0.5_dp * sum(wp%option%boxVecs, dim=2)
    frac(:) = matmul(boxMiddle - cellMiddle, recVecs2p)
    wp%option%origin(:) = wp%option%origin(:) - matmul(wp%xml%geo%latVecs, real(anint(frac), dp))
    wp%option%totGridOrig(:) = wp%option%totGridOrig(:) - matmul(wp%xml%geo%latVecs,&
        & real(anint(frac), dp))

    ! We need all cells around, which could contain atoms in the sphere, drawn from the center of
    ! the unit cell, containing the entire plotted region
    mDist = 0.0_dp
    do i1 = 0, 1
      do i2 = 0, 1
        do i3 = 0, 1
          cubeCorner(:) = wp%option%origin(:) + matmul(wp%option%boxVecs, real([i1, i2, i3], dp))
          dist = sqrt(sum((cubeCorner - cellMiddle)**2))
          mDist = max(dist, mDist)
        end do
      end do
    end do

    call getCellTranslations(cellVec, rCellVec, wp%xml%geo%latVecs, recVecs2p, mDist)

    ! Check all atoms in the shifted cells, if they fall in the plotted region
    call init(coordList)
    call init(speciesList)

    do iCell = 1, size(rCellVec, dim=2)
      do iAtom = 1, wp%xml%geo%nAtom
        coord(:) = wp%xml%geo%coords(:, iAtom) + rCellVec(:, iCell)
        frac(:) = matmul(invBoxVecs, coord - wp%option%origin(:))
        if (all(frac > - 1e-04_dp) .and. all(frac < 1.0_dp + 1e-04_dp)) then
          call append(coordList, coord)
          call append(speciesList, wp%xml%geo%species(iAtom))
        end if
      end do
    end do

    deallocate(wp%xml%geo%coords)
    deallocate(wp%xml%geo%species)

    wp%xml%geo%nAtom = len(coordList)

    allocate(wp%xml%geo%coords(3, wp%xml%geo%nAtom))
    allocate(wp%xml%geo%species(wp%xml%geo%nAtom))

    call asArray(coordList, wp%xml%geo%coords)
    call asArray(speciesList, wp%xml%geo%species)

    deallocate(cellVec)
    deallocate(rCellVec)

  end if

  ! Calculate mapping from the values which are in the level index to array indices. The 'Prime'
  ! variables contain these mappings.
  allocate(iLPrime(size(wp%option%levelIndex, dim=2)))
  allocate(iKPointPrime(size(wp%option%levelIndex, dim=2)))
  allocate(iSpinPrime(size(wp%option%levelIndex, dim=2)))

  iLPrime(:) = 0
  iKPointPrime(:) = 0
  iSpinPrime(:) = 0

  levelCounter = 1
  kpointCounter = 1
  spinCounter = 1

  LNum = 1
  KNum = 1
  SNum = 1

  do ii = 1, size(wp%option%levelIndex, dim=2)
    if (ii .eq. 1) then
      iLPrime(ii) = 1
    else
      if (any(wp%option%levelIndex(1, 1:ii - 1) .eq. wp%option%levelIndex(1, ii))) then
        tmparray = findloc(wp%option%levelIndex(1,:), wp%option%levelIndex(1, ii))
        iLPrime(ii) = iLPrime(tmparray(1))
      else
        iLPrime(ii) = iLPrime(ii - 1) + 1
        LNum = LNum + 1
      end if
      levelCounter = levelCounter + 1
    end if
  end do

  do ii = 1, size(wp%option%levelIndex, dim=2)
    if (ii .eq. 1) then
      iKPointPrime(ii) = 1
    else
      if (any(wp%option%levelIndex(2, 1:ii - 1) .eq. wp%option%levelIndex(2, ii))) then
        tmparray = findloc(wp%option%levelIndex(2, :), wp%option%levelIndex(2, ii))
        iKPointPrime(ii) = iKPointPrime(tmparray(1))
      else
        iKPointPrime(ii) = iKPointPrime(ii - 1) + 1
        KNum = KNum + 1
      end if
      kpointCounter = kpointCounter + 1
    end if
  end do

  do ii = 1, size(wp%option%levelIndex, dim=2)
    if (ii .eq. 1) then
      iSpinPrime(ii) = 1
    else
      if (any(wp%option%levelIndex(3, 1:ii - 1) .eq. wp%option%levelIndex(3, ii))) then
        tmparray = findloc(wp%option%levelIndex(3, :), wp%option%levelIndex(3, ii))
        iSpinPrime(ii) = iSpinPrime(tmparray(1))
      else
        iSpinPrime(ii) = iSpinPrime(ii - 1) + 1
        SNum = SNum + 1
      end if
      spinCounter = spinCounter + 1
    end if
  end do

  ! Extract all required k-points, levels and spins which are required for the calculation from the
  ! levelIndex
  allocate(requiredKPoints(KNum))
  allocate(requiredLevels(LNum))
  allocate(requiredSpins(SNum))

  requiredKPoints(:) = 0
  requiredLevels(:) = 0
  requiredSpins(:) = 0

  requiredLevels(1) = wp%option%levelIndex(1, 1)
  requiredKPoints(1) = wp%option%levelIndex(2, 1)
  requiredSpins(1) = wp%option%levelIndex(3, 1)

  i1 = 2
  i2 = 2
  i3 = 2

  do ii = 2, size(wp%option%levelIndex, dim=2)
    if (.not. (any(wp%option%levelIndex(2, 1:ii - 1) .eq. wp%option%levelIndex(2, ii)))) then
      requiredKPoints(i1) = wp%option%levelIndex(2, ii)
      i1 = i1 + 1
    end if
    if (.not. (any(wp%option%levelIndex(1, 1:ii - 1) .eq. wp%option%levelIndex(1, ii)))) then
      requiredLevels(i2) = wp%option%levelIndex(1, ii)
      i2 = i2 + 1
    end if
    if (.not. (any(wp%option%levelIndex(3, 1:ii - 1) .eq. wp%option%levelIndex(3, ii)))) then
      requiredSpins(i3) = wp%option%levelIndex(3, ii)
      i3 = i3 + 1
    end if
  end do

  ! Calculate the molecular orbitals
  if (wp%xml%tRealHam) then
    call subgridsToGlobalGrid(totGridDat, speciesGridsDat, wp%internal%molorb%coords,&
        & wp%eig%eigvecsReal, wp%option%levelIndex, wp%internal%orbitalToAtom,&
        & wp%internal%orbitalToSpecies, wp%option%parallelRegionNum, pCopyBuffers, totGridsDat,&
        & wp%option%gridInterType, addDensities=.false.)
  else
    pCopyBuffers => copyBuffersCplx
    call subgridsToGlobalGrid(totGridDat, speciesGridsDat, wp%internal%molorb%coords,&
          & wp%eig%eigvecsCplx, wp%option%levelIndex, wp%internal%orbitalToAtom,&
          & wp%internal%orbitalToSpecies, wp%option%parallelRegionNum, pCopyBuffers,&
          & totGridsDatCplx, requiredKPoints, requiredLevels, requiredSpins,&
          & wp%internal%molorb%CellVec, wp%option%gridInterType,  addDensities=.false., &
          & kPointsandWeights=wp%xml%kPointsandWeight)
  end if

  ! Process the molecular orbitals and write them to the disc
  lpProcessStates: do ii = 1, size(wp%option%levelIndex, dim=2)

    iLevel = wp%option%levelIndex(1, ii)
    iKPoint = wp%option%levelIndex(2, ii)
    iSpin = wp%option%levelIndex(3, ii)

    ! The 'Index' variables are the indices of the arrays which contain the data for the i'th state
    kIndex = iKPointPrime(ii)
    lIndex = iLPrime(ii)
    sIndex = iSpinPrime(ii)

    ! Build charge if needed for total charge or if it was explicitely required
    tPlotLevel = any(wp%option%plottedSpins == iSpin) .and. any(wp%option%plottedKPoints ==&
        & iKPoint) .and. any(wp%option%plottedLevels == iLevel)

    if (wp%option%tCalcTotChrg .or. (tPlotLevel .and. (wp%option%tPlotChrg .or.&
        & wp%option%tPlotChrgDiff))) then

      if (wp%xml%tRealHam) then
        buffer(:,:,:) = totGridsDat(:,:,:,ii)**2
      else
        buffer(:,:,:) = abs(totGridsDatCplx(:,:,:, lIndex, kIndex, sIndex))**2
      end if

      if (wp%option%tCalcTotChrg) then
        totChrg(:,:,:) = totChrg(:,:,:) + wp%xml%occupations(iLevel, iKPoint, iSpin) * buffer(:,:,:)
      end if

      if (wp%option%tPlotTotSpin) then
        if (iSpin .eq. 1) then
          spinUp(:,:,:) = spinUp(:,:,:) + wp%xml%occupations(iLevel, iKPoint, iSpin) * buffer(:,:,:)
        else
          spinDown(:,:,:) = spinDown(:,:,:) + wp%xml%occupations(iLevel, iKPoint, iSpin) *&
              & buffer(:,:,:)
        end if
      end if

      sumChrg = sum(buffer) * wp%internal%gridVol

      if (wp%option%tVerbose) then
        write(stdout, "(I5,I7,I7,A8,F12.6,F12.6)") iSpin, iKPoint, iLevel, "calc", sumChrg,&
            & wp%xml%occupations(iLevel, iKPoint, iSpin)
      end if

    end if

    ! Build and dump desired properties of the current level
    if (tPlotLevel) then

      if (wp%option%tPlotChrg) then
        write (comments(2), "('Calc-Id:',I11,', Spin:',I2,', K-Point:',I6,', State:',I6,&
            & ', abs2')") wp%xml%identity, iSpin, iKPoint, iLevel
        fileName = "wp-" // i2c(iSpin) // "-" // i2c(iKPoint) // "-" //i2c(iLevel) // "-abs2.cube"
        call writeCubeFile(wp%xml%geo, wp%atomicNumber%atomicNumbers, wp%internal%totGridVec,&
            & wp%option%totGridOrig, buffer, fileName, comments, wp%option%repeatBox)
        write(stdout, "(A)") "File '" // trim(fileName) // "' written"
      end if

      ! Plot charge difference
      if (wp%option%tPlotChrgDiff) then
        buffer = buffer - (sumChrg / sumAtomicChrg) * atomDensity(:,:,:)
        write (comments(2), "('Calc-Id:',I11,', Spin:',I2,', K-Point:',I6,', State:',I6,&
            & ', abs2diff')") wp%xml%identity, iSpin, iKPoint, iLevel
        fileName = "wp-" // i2c(iSpin) // "-" // i2c(iKPoint) // "-" //i2c(iLevel) //&
            & "-abs2diff.cube"
        call writeCubeFile(wp%xml%geo, wp%atomicNumber%atomicNumbers, wp%internal%totGridVec,&
            & wp%option%totGridOrig, buffer, fileName, comments, wp%option%repeatBox)
        write(stdout, "(A)") "File '" // trim(fileName) // "' written"
      end if

      ! Plot real part of WFs
      if (wp%option%tPlotReal) then

        if (wp%xml%tRealHam) then
          buffer(:,:,:) = totGridsDat(:,:,:,ii)
        else
          buffer(:,:,:) = real(totGridsDatCplx(:,:,:, lIndex, kIndex, sIndex), dp)
        end if

        write (comments(2), "('Calc-Id:',I11,', Spin:',I2,', K-Point:',I6,', State:',I6,&
            & ', real')") wp%xml%identity, iSpin, iKPoint, iLevel
        fileName = "wp-" // i2c(iSpin) // "-" // i2c(iKPoint) // "-" //i2c(iLevel) // "-real.cube"
        call writeCubeFile(wp%xml%geo, wp%atomicNumber%atomicNumbers, wp%internal%totGridVec,&
            & wp%option%totGridOrig, buffer, fileName, comments, wp%option%repeatBox)
        write(stdout, "(A)") "File '" // trim(fileName) // "' written"

      end if

      ! Plot imaginary part of WFs
      if (wp%option%tPlotImag) then
        buffer(:,:,:) = aimag(totGridsDatCplx(:,:,:, lIndex, kIndex, sIndex))
        write (comments(2), "('Calc-Id:',I11,', Spin:',I2,', K-Point:',I6,', State:',I6,&
            & ', imag')") wp%xml%identity, iSpin, iKPoint, iLevel
        fileName = "wp-" // i2c(iSpin) // "-" // i2c(iKPoint) // "-" //i2c(iLevel) // "-imag.cube"
        call writeCubeFile(wp%xml%geo, wp%atomicNumber%atomicNumbers, wp%internal%totGridVec,&
            & wp%option%totGridOrig, buffer, fileName, comments, wp%option%repeatBox)
        write(stdout, "(A)") "File '" // trim(fileName) // "' written"
      end if

    end if

  end do lpProcessStates

  ! Dump total charge, if required
  if (wp%option%tCalcTotChrg) then
    sumTotChrg = sum(totChrg) * wp%internal%gridVol
  end if

  ! Plot total charge
  if (wp%option%tPlotTotChrg) then
    write (comments(2), "('Calc-Id:',I11,', abs2')") wp%xml%identity
    fileName = "wp-abs2.cube"
    call writeCubeFile(wp%xml%geo, wp%atomicNumber%atomicNumbers, wp%internal%totGridVec,&
        & wp%option%totGridOrig, totChrg, fileName, comments, wp%option%repeatBox)
    write(stdout, "(A)") "File '" // trim(fileName) // "' written"

    if (wp%option%tVerbose) then
      write(stdout, "(/,'Total charge:',F12.6,/)") sumTotChrg
    end if
  end if

  ! Dump total charge difference
  if (wp%option%tPlotTotDiff) then
    buffer = totChrg - (sumTotChrg / sumAtomicChrg) * atomDensity(:,:,:)
    write (comments(2), "('Calc-Id:',I11,', abs2diff')") wp%xml%identity
    fileName = 'wp-abs2diff.cube'
    call writeCubeFile(wp%xml%geo, wp%atomicNumber%atomicNumbers, wp%internal%totGridVec, &
        & wp%option%totGridOrig,&
        & buffer, fileName, comments, wp%option%repeatBox)
    write(stdout, "(A)") "File '" // trim(fileName) // "' written"
  end if

  ! Dump spin polarisation
  if (wp%option%tPlotTotSpin) then
    buffer = spinUp - spinDown
    write (comments(2), "('Calc-Id:',I11,', spinpol')") wp%xml%identity
    fileName = 'wp-spinpol.cube'
    call writeCubeFile(wp%xml%geo, wp%atomicNumber%atomicNumbers, wp%internal%totGridVec, &
        & wp%option%totGridOrig,&
        & buffer, fileName, comments, wp%option%repeatBox)
    write(stdout, "(A)") "File '" // trim(fileName) // "' written"
  end if
  

#:if WITH_MPI
  call mpifx_finalize()
#:endif


contains

  !> Writes a 3D function as cube file.
  subroutine writeCubeFile(geo, atomicNumbers, gridVecs, origin, gridVal, fileName, comments,&
      & repeatBox)

    !> Geometry information about the structure
    type(TGeometry), intent(in) :: geo

    !> Atomic numbers corresponding to each species
    integer, intent(in) :: atomicNumbers(:)

    !> Grid vectors
    real(dp), intent(in) :: gridVecs(:,:)

    !> Origin of the grid
    real(dp), intent(in) :: origin(:)

    !> Value of the 3D function in the grid points
    real(dp), intent(in) :: gridVal(:,:,:)

    !> Name of the file to create
    character(len=*), intent(in) :: fileName

    !> First two comment lines of the file
    character(len=*), intent(in), optional :: comments(:)

    !> How often the grid should be repeated along the direction of the grid vectors
    integer, intent(in), optional :: repeatBox(:)

    integer, parameter :: bufferSize = 6
    real(dp) :: buffer(bufferSize)
    character(len=*), parameter :: formBuffer = "(6E16.8)"
    integer :: rep(3)
    integer, save :: fd = -1
    integer :: ii, i1, i2, i3, ir1, ir2, ir3

    @:ASSERT(size(atomicNumbers) == size(geo%speciesNames))
    @:ASSERT(all(shape(gridVecs) == [3, 3]))
    @:ASSERT(size(origin) == 3)
    @:ASSERT(all(shape(gridVal) >= [0, 0, 0]))

  #:block DEBUG_CODE
    if (present(comments)) then
      @:ASSERT(size(comments) == 2)
    end if
    if (present(repeatBox)) then
      @:ASSERT(size(repeatBox) == 3)
      @:ASSERT(all(repeatBox > 0))
    end if
  #:endblock DEBUG_CODE

    if (present(repeatBox)) then
      rep(:) = repeatBox(:)
    else
      rep(:) = [1, 1, 1]
    end if

    if (fd == -1) then
      fd = getFileId()
    end if

    open(fd, file=fileName, action="write", status="replace")
    if (present(comments)) then
      write (fd, "(A)") trim(comments(1))
      write (fd, "(A)") trim(comments(2))
    else
      write (fd, "(A)") "Made by waveplot"
      write (fd, *)
    end if
    write (fd,"(I5,3F12.6)") geo%nAtom, origin(:)
    write (fd,"(I5,3F12.6)") rep(1) * size(gridVal, dim=1), gridVecs(:,1)
    write (fd,"(I5,3F12.6)") rep(2) * size(gridVal, dim=2), gridVecs(:,2)
    write (fd,"(I5,3F12.6)") rep(3) * size(gridVal, dim=3), gridVecs(:,3)
    do ii = 1, geo%nAtom
      write (fd, "(I5,4F12.6)") atomicNumbers(geo%species(ii)), 0.0_dp, &
          &geo%coords(:, ii)
    end do

    do ir1 = 1, rep(1)
      do i1 = 1, size(gridVal, dim=1)
        do ir2 = 1, rep(2)
          do i2 = 1, size(gridVal, dim=2)
            do ir3 = 1, rep(3)
              do i3 = 1, size(gridVal, dim=3)
                ii = mod(i3 - 1, bufferSize) + 1
                buffer(ii) = gridVal(i1, i2, i3)
                if (ii == bufferSize) then
                  write (fd,formBuffer) real(buffer)
                end if
              end do
              if (ii /= bufferSize) then
                write (fd, "(" // i2c(ii) // "E16.8)") real(buffer(:ii))
              end if
            end do
          end do
        end do
      end do
    end do
    close(fd)

  end subroutine writeCubeFile

end program waveplot
