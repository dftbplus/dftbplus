!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Program for plotting molecular orbitals as cube files.
program waveplot

  use dftbp_common_accuracy, only : dp
  use dftbp_common_globalenv, only : stdOut
  use dftbp_dftb_periodic, only : getCellTranslations
  use dftbp_io_charmanip, only : i2c
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
  real(dp) :: cellMiddle(3), boxMiddle(3), frac(3), cubeCorner(3), coord(3), shift(3)
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
  allocate(buffer(wp%opt%nTotPoints(1), wp%opt%nTotPoints(2), wp%opt%nTotPoints(3)))
  buffer(:,:,:) = 0.0_dp
  pBuffer => buffer
  allocate(totData(wp%opt%nTotPoints(1), wp%opt%nTotPoints(2), wp%opt%nTotPoints(3)))
  totData(:,:,:) = 0.0_dp
  pTotData => totData
  allocate(copyBuffers(wp%opt%nTotPoints(1), wp%opt%nTotPoints(2), wp%opt%nTotPoints(3),&
      & size(wp%opt%levelIndex, dim=2)))
  copyBuffers(:,:,:,:) = 0.0_dp
  pCopyBuffers => copyBuffers
  allocate(copyBuffersCplx(wp%opt%nTotPoints(1), wp%opt%nTotPoints(2),&
      & wp%opt%nTotPoints(3), size(wp%opt%levelIndex, dim=2)))
  copyBuffersCplx(:,:,:,:) = 0.0_dp

  if (wp%opt%tCalcTotChrg) then
    allocate(totChrg(wp%opt%nTotPoints(1), wp%opt%nTotPoints(2), wp%opt%nTotPoints(3)))
    totChrg(:,:,:) = 0.0_dp
    if (wp%opt%tPlotTotSpin) then
      allocate(spinUp(wp%opt%nTotPoints(1), wp%opt%nTotPoints(2), wp%opt%nTotPoints(3)))
      spinUp(:,:,:) = 0.0_dp
      allocate(spinDown(wp%opt%nTotPoints(1), wp%opt%nTotPoints(2), wp%opt%nTotPoints(3)))
      spinDown(:,:,:) = 0.0_dp
    end if
  end if

  write (comments(1), "(A)") "Cube file generated by WAVEPLOT from data created by DFTB+"

  write(stdout, "(A)") "Origin"
  write(stdout, "(2X,3(F0.5,1X))") wp%opt%origin(:)
  write(stdout, "(A)") "Total-Grid Origin"
  write(stdout, "(2X,3(F0.5,1X))") wp%opt%totGridOrig
  write(stdout, "(A)") "Box"

  do jj = 1, 3
    write(stdout, "(2X,3(F0.5,1X))") wp%opt%boxVecs(:, jj)
  end do

  write(stdout, "(A)") "Spatial resolution [1/Bohr]:"
  write(stdout, "(2X,3(F0.5,1X))") 1.0_dp / sqrt(sum(wp%loc%totGridVec**2, dim=1))
  write(stdout, *)

  ! Initialise total grid and volumetric data
  call TGrid_init(totGrid, wp%opt%totGridOrig, wp%loc%totGridVec, wp%opt%nTotPoints)
  call TGridData_init(totGridDat, totGrid, pBuffer, rwTabulationType=wp%opt%rwTabulationType)

  ! Initialise atomic volumetric data
  allocate(speciesGrids(wp%input%geo%nSpecies))
  allocate(speciesGridsDat(maxval(wp%loc%orbitalToSpecies)))
  allocate(speciesRty(maxval(wp%loc%orbitalToSpecies)))
  allocate(speciesChrg(wp%opt%nSpPoints(1), wp%opt%nSpPoints(2), wp%opt%nSpPoints(3),&
      & maxval(wp%loc%orbitalToSpecies)))
  speciesChrg(:,:,:,:) = 0.0_dp

  ind = 1

  ! Initialise atomic grids and pretabulate orbital-species grids
  do iSpecies = 1, wp%input%geo%nSpecies
    maxAng = maxval(wp%basis%basis(iSpecies)%angMoms)
    call TGrid_init(speciesGrids(iSpecies), wp%loc%speciesGridsOrigs(:, iSpecies),&
        & wp%loc%speciesGridsVecs(:, :, iSpecies), wp%opt%nSpPoints)
    call TRealTessY_init(speciesRty(iSpecies), speciesGrids(iSpecies), maxAng)
    do iAng = wp%loc%molorb%iStos(iSpecies), wp%loc%molorb%iStos(iSpecies + 1) - 1
      iL = wp%loc%molorb%angMoms(iAng)
      do iM = - iL, iL
        pSpeciesChrg => speciesChrg(:,:,:, ind)
        call TGridData_init(speciesGridsDat(ind), speciesGrids(iSpecies), pSpeciesChrg,&
            & rwTabulationType=wp%opt%rwTabulationType)
        call speciesGridsDat(ind)%tabulateBasis(speciesRty(iSpecies), iL, iM, &
            & wp%loc%molorb%rwfs(iAng)%rwf, wp%loc%molorb%rwfs(iAng)%exps, &
            & wp%loc%molorb%rwfs(iAng)%aa)
        ind = ind + 1
      end do
    end do
  end do

  ! Create density superposition of the atomic orbitals. Occupation is
  ! distributed equally on orbitals with the same angular momentum.
  if (wp%opt%tCalcAtomDens) then

    allocate(atomicChrg(wp%opt%nTotPoints(1), wp%opt%nTotPoints(2), wp%opt%nTotPoints(3)))
    atomicChrg(:,:,:) = 0.0_dp
    pAtomicChrg => atomicChrg

    ! Initialise volumetric data
    call TGridData_init(atomicGridDat, totGrid, pAtomicChrg,&
        & rwTabulationType=wp%opt%rwTabulationType)

    ! Calculate total charge of atomic densities.
    call subgridsToGlobalGrid(atomicGridDat, speciesGridsDat, atomDensity, &
        & wp%loc%molorb%coords, wp%loc%orbitalOcc(:, 1), wp%loc%orbitalToAtom, &
        & wp%loc%orbitalToSpecies, wp%opt%parallelRegionNum, .true., wp%opt%gridInterType)

    sumAtomicChrg = sum(atomDensity) * wp%loc%gridVol

    if (wp%opt%tVerbose) then
      write(stdout, "('Total charge of atomic densities:',F12.6,/)") sumAtomicChrg
    end if

    ! Plot total atomic density
    if (wp%opt%tPlotAtomDens) then
      write (comments(2), "('Calc-Id:', I11, ', atomdens')") wp%input%identity
      fileName = "wp-atomdens.cube"
      call writeCubeFile(wp%input%geo, wp%aNr%atomicNumbers, wp%loc%totGridVec,&
          & wp%opt%totGridOrig, atomDensity, fileName, comments, wp%opt%repeatBox)
      write(stdout, "(A)") "File '" // trim(fileName) // "' written"
    end if

  end if

  if (wp%opt%tVerbose) then
    write(stdout, "(/,A5,' ',A6,' ',A6,' ',A7,' ',A11,' ',A11)") "Spin", "KPoint", "State",&
        & "Action", "Norm", "W. Occup."
  end if

  ! Fold in coordinates
  if (wp%opt%tFoldCoords) then
    call invert33(invBoxVecs, wp%opt%boxVecs)
    call invert33(recVecs2p, wp%input%geo%latVecs)
    recVecs2p = reshape(recVecs2p, [3, 3], order=[2, 1])
    call wp%boundaryCond%foldCoordsToCell(wp%input%geo%coords(:,:), wp%input%geo%latVecs)
  end if

  ! Fill the box with atoms
  if (wp%opt%tFillBox) then

    ! Shift plotted region by integer lattice vectors, to have its center as close to the center
    ! of the lattice unit cell as possible.
    cellMiddle(:) = 0.5_dp * sum(wp%input%geo%latVecs, dim=2)
    boxMiddle(:) = wp%opt%origin(:) + 0.5_dp * sum(wp%opt%boxVecs, dim=2)
    ! Workaround for intel 2021 ICE, replacing matmul(boxMiddle - cellMiddle, recVecs2p)
    shift(:) = boxMiddle - cellMiddle
    frac(:) = matmul(shift, recVecs2p)
    wp%opt%origin(:) = wp%opt%origin(:) - matmul(wp%input%geo%latVecs, real(anint(frac), dp))
    wp%opt%totGridOrig(:) = wp%opt%totGridOrig(:)&
        & - matmul(wp%input%geo%latVecs, real(anint(frac), dp))
    ! We need all cells around, which could contain atoms in the sphere, drawn from the center of
    ! the unit cell, containing the entire plotted region
    mDist = 0.0_dp
    do i1 = 0, 1
      do i2 = 0, 1
        do i3 = 0, 1
          cubeCorner(:) = wp%opt%origin(:) + matmul(wp%opt%boxVecs, real([i1, i2, i3], dp))
          dist = sqrt(sum((cubeCorner - cellMiddle)**2))
          mDist = max(dist, mDist)
        end do
      end do
    end do

    call getCellTranslations(cellVec, rCellVec, wp%input%geo%latVecs, recVecs2p, mDist)

    ! Check all atoms in the shifted cells, if they fall in the plotted region
    call init(coordList)
    call init(speciesList)

    do iCell = 1, size(rCellVec, dim=2)
      do iAtom = 1, wp%input%geo%nAtom
        coord(:) = wp%input%geo%coords(:, iAtom) + rCellVec(:, iCell)
        frac(:) = matmul(invBoxVecs, coord - wp%opt%origin(:))
        if (all(frac > - 1e-04_dp) .and. all(frac < 1.0_dp + 1e-04_dp)) then
          call append(coordList, coord)
          call append(speciesList, wp%input%geo%species(iAtom))
        end if
      end do
    end do

    deallocate(wp%input%geo%coords)
    deallocate(wp%input%geo%species)

    wp%input%geo%nAtom = len(coordList)

    allocate(wp%input%geo%coords(3, wp%input%geo%nAtom))
    allocate(wp%input%geo%species(wp%input%geo%nAtom))

    call asArray(coordList, wp%input%geo%coords)
    call asArray(speciesList, wp%input%geo%species)

    deallocate(cellVec)
    deallocate(rCellVec)

  end if

  ! Calculate mapping from the values which are in the level index to array indices. The 'Prime'
  ! variables contain these mappings.
  allocate(iLPrime(size(wp%opt%levelIndex, dim=2)))
  allocate(iKPointPrime(size(wp%opt%levelIndex, dim=2)))
  allocate(iSpinPrime(size(wp%opt%levelIndex, dim=2)))

  iLPrime(:) = 0
  iKPointPrime(:) = 0
  iSpinPrime(:) = 0

  levelCounter = 1
  kpointCounter = 1
  spinCounter = 1

  LNum = 1
  KNum = 1
  SNum = 1

  do ii = 1, size(wp%opt%levelIndex, dim=2)
    if (ii .eq. 1) then
      iLPrime(ii) = 1
    else
      if (any(wp%opt%levelIndex(1, 1:ii - 1) .eq. wp%opt%levelIndex(1, ii))) then
        tmparray = findloc(wp%opt%levelIndex(1,:), wp%opt%levelIndex(1, ii))
        iLPrime(ii) = iLPrime(tmparray(1))
      else
        iLPrime(ii) = iLPrime(ii - 1) + 1
        LNum = LNum + 1
      end if
      levelCounter = levelCounter + 1
    end if
  end do

  do ii = 1, size(wp%opt%levelIndex, dim=2)
    if (ii .eq. 1) then
      iKPointPrime(ii) = 1
    else
      if (any(wp%opt%levelIndex(2, 1:ii - 1) .eq. wp%opt%levelIndex(2, ii))) then
        tmparray = findloc(wp%opt%levelIndex(2, :), wp%opt%levelIndex(2, ii))
        iKPointPrime(ii) = iKPointPrime(tmparray(1))
      else
        iKPointPrime(ii) = iKPointPrime(ii - 1) + 1
        KNum = KNum + 1
      end if
      kpointCounter = kpointCounter + 1
    end if
  end do

  do ii = 1, size(wp%opt%levelIndex, dim=2)
    if (ii .eq. 1) then
      iSpinPrime(ii) = 1
    else
      if (any(wp%opt%levelIndex(3, 1:ii - 1) .eq. wp%opt%levelIndex(3, ii))) then
        tmparray = findloc(wp%opt%levelIndex(3, :), wp%opt%levelIndex(3, ii))
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

  requiredLevels(1) = wp%opt%levelIndex(1, 1)
  requiredKPoints(1) = wp%opt%levelIndex(2, 1)
  requiredSpins(1) = wp%opt%levelIndex(3, 1)

  i1 = 2
  i2 = 2
  i3 = 2

  do ii = 2, size(wp%opt%levelIndex, dim=2)
    if (.not. (any(wp%opt%levelIndex(2, 1:ii - 1) .eq. wp%opt%levelIndex(2, ii)))) then
      requiredKPoints(i1) = wp%opt%levelIndex(2, ii)
      i1 = i1 + 1
    end if
    if (.not. (any(wp%opt%levelIndex(1, 1:ii - 1) .eq. wp%opt%levelIndex(1, ii)))) then
      requiredLevels(i2) = wp%opt%levelIndex(1, ii)
      i2 = i2 + 1
    end if
    if (.not. (any(wp%opt%levelIndex(3, 1:ii - 1) .eq. wp%opt%levelIndex(3, ii)))) then
      requiredSpins(i3) = wp%opt%levelIndex(3, ii)
      i3 = i3 + 1
    end if
  end do

  ! Calculate the molecular orbitals
  if (wp%input%tRealHam) then
    call subgridsToGlobalGrid(totGridDat, speciesGridsDat, wp%loc%molorb%coords,&
        & wp%eig%eigvecsReal, wp%opt%levelIndex, wp%loc%orbitalToAtom,&
        & wp%loc%orbitalToSpecies, wp%opt%parallelRegionNum, pCopyBuffers, totGridsDat,&
        & wp%opt%gridInterType, addDensities=.false.)
  else
    pCopyBuffers => copyBuffersCplx
    call subgridsToGlobalGrid(totGridDat, speciesGridsDat, wp%loc%molorb%coords,&
          & wp%eig%eigvecsCplx, wp%opt%levelIndex, wp%loc%orbitalToAtom,&
          & wp%loc%orbitalToSpecies, wp%opt%parallelRegionNum, pCopyBuffers,&
          & totGridsDatCplx, requiredKPoints, requiredLevels, requiredSpins,&
          & wp%loc%molorb%CellVec, wp%opt%gridInterType,  addDensities=.false., &
          & kPointsandWeights=wp%input%kPointsandWeight)
  end if

  ! Process the molecular orbitals and write them to the disc
  lpProcessStates: do ii = 1, size(wp%opt%levelIndex, dim=2)

    iLevel = wp%opt%levelIndex(1, ii)
    iKPoint = wp%opt%levelIndex(2, ii)
    iSpin = wp%opt%levelIndex(3, ii)

    ! The 'Index' variables are the indices of the arrays which contain the data for the i'th state
    kIndex = iKPointPrime(ii)
    lIndex = iLPrime(ii)
    sIndex = iSpinPrime(ii)

    ! Build charge if needed for total charge or if it was explicitely required
    tPlotLevel = any(wp%opt%plottedSpins == iSpin) .and. any(wp%opt%plottedKPoints ==&
        & iKPoint) .and. any(wp%opt%plottedLevels == iLevel)

    if (wp%opt%tCalcTotChrg .or. (tPlotLevel .and. (wp%opt%tPlotChrg .or.&
        & wp%opt%tPlotChrgDiff))) then

      if (wp%input%tRealHam) then
        buffer(:,:,:) = totGridsDat(:,:,:,ii)**2
      else
        buffer(:,:,:) = abs(totGridsDatCplx(:,:,:, lIndex, kIndex, sIndex))**2
      end if

      if (wp%opt%tCalcTotChrg) then
        totChrg(:,:,:) = totChrg(:,:,:) + wp%input%occupations(iLevel, iKPoint, iSpin) * buffer(:,:,:)
      end if

      if (wp%opt%tPlotTotSpin) then
        if (iSpin .eq. 1) then
          spinUp(:,:,:) = spinUp(:,:,:) + wp%input%occupations(iLevel, iKPoint, iSpin) * buffer(:,:,:)
        else
          spinDown(:,:,:) = spinDown(:,:,:) + wp%input%occupations(iLevel, iKPoint, iSpin) *&
              & buffer(:,:,:)
        end if
      end if

      sumChrg = sum(buffer) * wp%loc%gridVol

      if (wp%opt%tVerbose) then
        write(stdout, "(I5,I7,I7,A8,F12.6,F12.6)") iSpin, iKPoint, iLevel, "calc", sumChrg,&
            & wp%input%occupations(iLevel, iKPoint, iSpin)
      end if

    end if

    ! Build and dump desired properties of the current level
    if (tPlotLevel) then

      if (wp%opt%tPlotChrg) then
        write (comments(2), "('Calc-Id:',I11,', Spin:',I2,', K-Point:',I6,', State:',I6,&
            & ', abs2')") wp%input%identity, iSpin, iKPoint, iLevel
        fileName = "wp-" // i2c(iSpin) // "-" // i2c(iKPoint) // "-" //i2c(iLevel) // "-abs2.cube"
        call writeCubeFile(wp%input%geo, wp%aNr%atomicNumbers, wp%loc%totGridVec,&
            & wp%opt%totGridOrig, buffer, fileName, comments, wp%opt%repeatBox)
        write(stdout, "(A)") "File '" // trim(fileName) // "' written"
      end if

      ! Plot charge difference
      if (wp%opt%tPlotChrgDiff) then
        buffer = buffer - (sumChrg / sumAtomicChrg) * atomDensity(:,:,:)
        write (comments(2), "('Calc-Id:',I11,', Spin:',I2,', K-Point:',I6,', State:',I6,&
            & ', abs2diff')") wp%input%identity, iSpin, iKPoint, iLevel
        fileName = "wp-" // i2c(iSpin) // "-" // i2c(iKPoint) // "-" //i2c(iLevel) //&
            & "-abs2diff.cube"
        call writeCubeFile(wp%input%geo, wp%aNr%atomicNumbers, wp%loc%totGridVec,&
            & wp%opt%totGridOrig, buffer, fileName, comments, wp%opt%repeatBox)
        write(stdout, "(A)") "File '" // trim(fileName) // "' written"
      end if

      ! Plot real part of WFs
      if (wp%opt%tPlotReal) then

        if (wp%input%tRealHam) then
          buffer(:,:,:) = totGridsDat(:,:,:,ii)
        else
          buffer(:,:,:) = real(totGridsDatCplx(:,:,:, lIndex, kIndex, sIndex), dp)
        end if

        write (comments(2), "('Calc-Id:',I11,', Spin:',I2,', K-Point:',I6,', State:',I6,&
            & ', real')") wp%input%identity, iSpin, iKPoint, iLevel
        fileName = "wp-" // i2c(iSpin) // "-" // i2c(iKPoint) // "-" //i2c(iLevel) // "-real.cube"
        call writeCubeFile(wp%input%geo, wp%aNr%atomicNumbers, wp%loc%totGridVec,&
            & wp%opt%totGridOrig, buffer, fileName, comments, wp%opt%repeatBox)
        write(stdout, "(A)") "File '" // trim(fileName) // "' written"

      end if

      ! Plot imaginary part of WFs
      if (wp%opt%tPlotImag) then
        buffer(:,:,:) = aimag(totGridsDatCplx(:,:,:, lIndex, kIndex, sIndex))
        write (comments(2), "('Calc-Id:',I11,', Spin:',I2,', K-Point:',I6,', State:',I6,&
            & ', imag')") wp%input%identity, iSpin, iKPoint, iLevel
        fileName = "wp-" // i2c(iSpin) // "-" // i2c(iKPoint) // "-" //i2c(iLevel) // "-imag.cube"
        call writeCubeFile(wp%input%geo, wp%aNr%atomicNumbers, wp%loc%totGridVec,&
            & wp%opt%totGridOrig, buffer, fileName, comments, wp%opt%repeatBox)
        write(stdout, "(A)") "File '" // trim(fileName) // "' written"
      end if

    end if

  end do lpProcessStates

  ! Dump total charge, if required
  if (wp%opt%tCalcTotChrg) then
    sumTotChrg = sum(totChrg) * wp%loc%gridVol
  end if

  ! Plot total charge
  if (wp%opt%tPlotTotChrg) then
    write (comments(2), "('Calc-Id:',I11,', abs2')") wp%input%identity
    fileName = "wp-abs2.cube"
    call writeCubeFile(wp%input%geo, wp%aNr%atomicNumbers, wp%loc%totGridVec,&
        & wp%opt%totGridOrig, totChrg, fileName, comments, wp%opt%repeatBox)
    write(stdout, "(A)") "File '" // trim(fileName) // "' written"

    if (wp%opt%tVerbose) then
      write(stdout, "(/,'Total charge:',F12.6,/)") sumTotChrg
    end if
  end if

  ! Dump total charge difference
  if (wp%opt%tPlotTotDiff) then
    buffer = totChrg - (sumTotChrg / sumAtomicChrg) * atomDensity(:,:,:)
    write (comments(2), "('Calc-Id:',I11,', abs2diff')") wp%input%identity
    fileName = 'wp-abs2diff.cube'
    call writeCubeFile(wp%input%geo, wp%aNr%atomicNumbers, wp%loc%totGridVec, &
        & wp%opt%totGridOrig,&
        & buffer, fileName, comments, wp%opt%repeatBox)
    write(stdout, "(A)") "File '" // trim(fileName) // "' written"
  end if

  ! Dump spin polarisation
  if (wp%opt%tPlotTotSpin) then
    buffer = spinUp - spinDown
    write (comments(2), "('Calc-Id:',I11,', spinpol')") wp%input%identity
    fileName = 'wp-spinpol.cube'
    call writeCubeFile(wp%input%geo, wp%aNr%atomicNumbers, wp%loc%totGridVec, &
        & wp%opt%totGridOrig,&
        & buffer, fileName, comments, wp%opt%repeatBox)
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

    open(newunit=fd, file=fileName, action="write", status="replace")
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
