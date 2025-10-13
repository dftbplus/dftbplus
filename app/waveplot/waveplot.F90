!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Program for plotting molecular orbitals as cube files.
program waveplot
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment, TEnvironment_init
  use dftbp_common_file, only : closeFile, openFile, TFileDescr
  use dftbp_common_globalenv, only : destructGlobalEnv, initGlobalEnv, stdOut
  use dftbp_dftb_periodic, only : getCellTranslations
  use dftbp_io_charmanip, only : i2c
  use dftbp_io_message, only : error, warning
  use dftbp_math_simplealgebra, only : invert33
  use dftbp_type_linkedlist, only : append, asArray, init, len, TListInt, TListRealR1
  use dftbp_type_typegeometry, only : TGeometry
  use dftbp_wavegrid, only : getAtomicDensities, getTotalChrg
#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : MPI_LOR, MPI_SUM, mpifx_allreduceip, mpifx_bcast
#:endif
  use waveplot_initwaveplot, only : TProgramVariables, TProgramVariables_init

  implicit none

  !> Container of program variables
  type(TProgramVariables), target :: wp

  !> Environment settings
  type(TEnvironment) :: env

  !> Pointer to real-valued grid
  real(dp), pointer :: gridValReal(:,:,:)

  !> Pointer to complex-valued grid
  complex(dp), pointer :: gridValCmpl(:,:,:)

  !> Arrays holding the volumetric data
  real(dp), allocatable :: buffer(:,:,:), totChrg(:,:,:), atomicChrg(:,:,:,:), spinUp(:,:,:)

  !> Summation of all grid points
  real(dp) :: sumTotChrg, sumChrg, sumAtomicChrg

  !> Indices of current level, K-point and spin
  integer :: levelIndex(3)

  !> Auxiliary variables
  integer :: i1, ioStat, nEig, iEig, iLevel, iKPoint, iSpin
  logical :: isFinished, doPlotLevel, hasIoError, doRequireIndividual, doRepeatBox, doNeedCharge

  call initGlobalEnv()
  call TEnvironment_init(env)

  ! Allocate resources
  call TProgramVariables_init(wp, env)
  write(stdOut, "(/,A,/)") "Starting main program"

  ! Allocating buffer for general grid, total charge and spin up
  allocate(buffer(wp%opt%nPoints(1), wp%opt%nPoints(2), wp%opt%nPoints(3)))
  if (wp%opt%doCalcTotChrg) then
    allocate(totChrg(wp%opt%nPoints(1), wp%opt%nPoints(2), wp%opt%nPoints(3)), source=0.0_dp)
    if (wp%opt%doPlotTotSpin) then
      allocate(spinUp(wp%opt%nPoints(1), wp%opt%nPoints(2), wp%opt%nPoints(3)), source=0.0_dp)
    end if
  end if
  hasIoError = .false.

  ! Repeat boxes if necessary
  doRepeatBox = product(wp%opt%repeatBox) > 1
  if (doRepeatBox) then
    call expandToSupercell(wp)
  end if

  write(stdOut, "(A)") "Origin"
  write(stdOut, "(2X,3(F0.5,1X))") wp%opt%origin
  write(stdOut, "(A)") "Box"
  do i1 = 1, 3
    write(stdOut, "(2X,3(F0.5,1X))") wp%opt%boxVecs(:, i1)
  end do
  write(stdOut, "(A)") "Spatial resolution [1/Bohr]:"
  write(stdOut, "(2X,3(F0.5,1X))") 1.0_dp / norm2(wp%loc%gridVec, dim=1)
  write(stdOut, *)

  ! Create density superposition of the atomic orbitals. Occupation is distributed equally on
  ! orbitals with the same angular momentum.
  if (wp%opt%doCalcAtomDens) then
    call calcAtomicDensities(wp, env, hasIoError, atomicChrg, sumAtomicChrg)
  end if

  if (wp%opt%beVerbose) then
    write(stdOut, "(/,A5,' ',A6,' ',A6,' ',A7,' ',A11,' ',A11)") "Spin", "KPoint", "State",&
        & "Action", "Norm", "W. Occup."
  end if

  if (wp%opt%doFoldCoords) then
    call wp%boundaryCond%foldCoordsToCell(wp%input%geo%coords, wp%input%geo%latVecs)
  end if

  ! If requested, repeat the unit cell and keep all atoms that fall in the plotted region
  if (wp%opt%doFillBox) then
    call fillPlottedRegion(wp)
  end if
 
  doRequireIndividual = wp%opt%doPlotChrgDiff &
                   .or. wp%opt%doPlotReal &
                   .or. wp%opt%doPlotImag &
                   .or. wp%opt%doPlotTotSpin

  ! Wavegrid supports fast inplace accumulation for total charge.
  ! This avoids having to store all states in memory and can offer a
  ! significant speedup for large systems.
  if (wp%opt%doCalcTotChrg .and. .not. doRequireIndividual) then
      call calcTotChrgInplace(wp, totChrg)
  end if

  if (doRequireIndividual) then
    ! Calculate the molecular orbitals and write them to the disk
    isFinished = .false.
    lpStates: do while (.not. isFinished)
      ! Get the next grid and its parameters
      if (wp%input%isRealHam) then
        call wp%loc%grid%next(gridValReal, levelIndex, isFinished)
      else
        call wp%loc%grid%next(gridValCmpl, levelIndex, isFinished)
      end if

      ! Was the current state requested separately?
      iLevel = levelIndex(1); iKPoint = levelIndex(2); iSpin = levelIndex(3)
      doPlotLevel = any(wp%opt%plottedSpins   == iSpin) &
            & .and. any(wp%opt%plottedKPoints == iKPoint) &
            & .and. any(wp%opt%plottedLevels  == iLevel)

      ! Square the wavefunction if neccessary
      doNeedCharge = wp%opt%doCalcTotChrg .or. (doPlotLevel .and. (wp%opt%doPlotChrg .or. wp%opt%doPlotChrgDiff))
      if (doNeedCharge) then
        if (wp%input%isRealHam) then
          buffer(:,:,:) = gridValReal**2
        else
          buffer(:,:,:) = abs(gridValCmpl)**2
        end if

        sumChrg = sum(buffer) * wp%loc%gridVol
        if (wp%opt%beVerbose) then
          write(stdOut, "(I5,I7,I7,A8,F12.6,F12.6)") iSpin, iKPoint, iLevel, "calc", sumChrg,&
              & wp%input%occupations(iLevel, iKPoint, iSpin)
        end if

        if (wp%opt%doCalcTotChrg) then
          totChrg(:,:,:) = totChrg + wp%input%occupations(iLevel, iKPoint, iSpin) * buffer
        end if
      end if

      ! Save spin up density before processing first level for spin down
      if (wp%opt%doPlotTotSpin .and. (iSpin == 1)) then
        spinUp(:,:,:) = spinUp + wp%input%occupations(iLevel, iKPoint, iSpin) * buffer
      end if

      ! Build and dump desired properties of the current level
      if (doPlotLevel) then
        if (wp%opt%doPlotChrg) then
          call writePropertyToCube("charge", buffer, wp, levelIndex, ioStat=ioStat)
          hasIoError = hasIoError .or. ioStat /= 0
        end if

        if (wp%opt%doPlotChrgDiff) then
          buffer(:,:,:) = buffer - (sumChrg / sumAtomicChrg) * atomicChrg(:,:,:,1)
          call writePropertyToCube("chargediff", buffer, wp, levelIndex, ioStat=ioStat)
          hasIoError = hasIoError .or. ioStat /= 0
        end if

        if (wp%opt%doPlotReal) then
          if (wp%input%isRealHam) then
            buffer(:,:,:) = gridValReal
          else
            buffer(:,:,:) = real(gridValCmpl, dp)
          end if
          call writePropertyToCube("real", buffer, wp, levelIndex, ioStat=ioStat)
          hasIoError = hasIoError .or. ioStat /= 0
        end if

        if (wp%opt%doPlotImag) then
          buffer(:,:,:) = aimag(gridValCmpl)
          call writePropertyToCube("imag", buffer, wp, levelIndex, ioStat=ioStat)
          hasIoError = hasIoError .or. ioStat /= 0
        end if
      end if
    end do lpStates
  end if

#:if WITH_MPI
  call mpifx_allreduceip(env%mpi%globalComm, hasIoError, MPI_LOR)
#:endif

  if (env%tGlobalLead .and. hasIoError) then
    call error("At least one of the processes encountered an I/O error.")
  end if

  ! Dump total charge, if required
  if (wp%opt%doCalcTotChrg) then
  #:if WITH_MPI
    call mpifx_allreduceip(env%mpi%globalComm, totChrg, MPI_SUM)
  #:endif
    sumTotChrg = sum(totChrg) * wp%loc%gridVol
  end if
  if (env%tGlobalLead .and. wp%opt%doPlotTotChrg) then
    call writePropertyToCube("total_charge", totChrg, wp)
    if (wp%opt%beVerbose) then
      write(stdOut, "(/,'Total charge:',F12.6,/)") sumTotChrg
    end if
  end if

  ! Dump total charge difference
  if (env%tGlobalLead .and. wp%opt%doPlotTotDiff) then
    buffer(:,:,:) = totChrg - (sumTotChrg / sumAtomicChrg) * atomicChrg(:,:,:,1)
    call writePropertyToCube("total_chargediff", buffer, wp)
  end if

#:if WITH_MPI
  ! Collect spin polarisation
  if (wp%opt%doPlotTotSpin) then
    call mpifx_allreduceip(env%mpi%globalComm, spinUp, MPI_SUM)
  end if
#:endif

  if (env%tGlobalLead .and. wp%opt%doPlotTotSpin) then
    buffer(:,:,:) = 2.0_dp * spinUp - totChrg
    call writePropertyToCube("spinpol", buffer, wp)
  end if

  call env%destruct()
  call destructGlobalEnv()


contains

  !> Calculates atomic densities and their sum and saves to disk if requested.
  subroutine calcAtomicDensities(wp, env, hasIoError, atomicChrg, sumAtomicChrg)
    type(TProgramVariables), intent(in) :: wp
    type(TEnvironment), intent(in) :: env
    logical, intent(inout) :: hasIoError
    real(dp), allocatable, intent(out) :: atomicChrg(:,:,:,:)
    real(dp), intent(out) :: sumAtomicChrg
    real(dp), allocatable :: orbitalOcc(:,:)
    integer :: iAtom, iSpecies, iOrb, mAng, ind, ioStat, iL

    allocate(atomicChrg(wp%opt%nPoints(1), wp%opt%nPoints(2), wp%opt%nPoints(3), 1))
    allocate(orbitalOcc(wp%input%nOrb, 1))
    if (env%tGlobalLead) then
      ind = 1
      do iAtom = 1, wp%input%geo%nAtom
        iSpecies = wp%input%geo%species(iAtom)
        do iOrb = 1, size(wp%basis%basis(iSpecies)%orbitals)
          iL = wp%basis%basis(iSpecies)%orbitals(iOrb)%angMom
          mAng = 2 * iL + 1
          orbitalOcc(ind:ind + mAng - 1,1) = wp%basis%referenceOccupations(iOrb, iSpecies) &
              & / real(mAng, dp)
          ind = ind + mAng
        end do
      end do
      call getAtomicDensities(wp%loc%molorb, orbitalOcc, atomicChrg)
      sumAtomicChrg = sum(atomicChrg) * wp%loc%gridVol

      if (wp%opt%beVerbose) then
        write(stdOut, "('Total charge of atomic densities:',F12.6,/)") sumAtomicChrg
      end if
      if (wp%opt%doPlotAtomDens) then
        call writePropertyToCube("atomdens", atomicChrg(:,:,:,1), wp, iostat=ioStat)
        hasIoError = hasIoError .or. ioStat /= 0
      end if
    end if
  #:if WITH_MPI
    call mpifx_bcast(env%mpi%globalComm, atomicChrg)
    call mpifx_bcast(env%mpi%globalComm, sumAtomicChrg)
  #:endif
    deallocate(orbitalOcc)
  end subroutine calcAtomicDensities
  
  !> Calculates the total charge in-place, without storing all states in memory.
  subroutine calcTotChrgInplace(wp, totChrg)
    type(TProgramVariables), intent(inout) :: wp
    real(dp), intent(out) :: totChrg(:,:,:)
    real(dp), allocatable :: totChrg4d(:,:,:,:), eigCoeffs(:)

    if (wp%opt%beVerbose) then
      write(stdOut, "(A)") "Calculating total charge in-place."
    end if

    allocate(totChrg4d(wp%opt%nPoints(1), wp%opt%nPoints(2), wp%opt%nPoints(3), 1), source=0.0_dp)
    ! Get occupation by state
    nEig = wp%loc%grid%nCached
    call wp%loc%grid%loadEigenvecs(nEig)
    allocate(eigCoeffs(nEig))

    do iEig = 1, nEig
        levelIndex = wp%loc%grid%levelIndex(:, iEig)
        iLevel = levelIndex(1); iKPoint = levelIndex(2); iSpin = levelIndex(3)
        eigCoeffs(iEig) = wp%input%occupations(iLevel, iKPoint, iSpin)
    end do
    if (wp%input%isRealHam) then
      call getTotalChrg(wp%loc%molorb, wp%loc%grid%eigenvecReal, &
        & totChrg4d, eigCoeffs, wp%opt%useGPU)
    else
      call getTotalChrg(wp%loc%molorb, wp%loc%grid%eigenvecCmpl, &
        & wp%loc%grid%kPoints, wp%loc%grid%levelIndex(2,:), totChrg4d, eigCoeffs, wp%opt%useGPU)
    end if
    totChrg(:,:,:) = totChrg4d(:,:,:,1)
    deallocate(eigCoeffs)
    deallocate(totChrg4d)
  end subroutine calcTotChrgInplace



  !> Repeats the unit cell to a supercell as described in wp%opt%repeatBox.
  subroutine expandToSupercell(wp)
    type(TProgramVariables), intent(inout) :: wp
    integer :: nBox, i1, i2, i3, iAtom, ind
    real(dp) :: shift(3)
    real(dp), allocatable :: coords(:,:)
    integer, allocatable :: species(:)
    
    nBox = product(wp%opt%repeatBox)
    @:ASSERT(nBox > 1)
    wp%input%nOrb = wp%input%nOrb * nBox

    if (wp%opt%beVerbose) then
      write(stdOut, "(A,I1,A)") "Expanding unit cell to supercell with ", nBox, " boxes."
    end if

    ! If doFillBox is off, coordinates must be repeated here.
    ! Otherwise the part for filling with atoms will do that.
    if (.not. wp%opt%doFillBox) then
      ! Store old coordinates and species
      allocate(coords(3, size(wp%input%geo%coords, dim=2)))
      allocate(species(size(wp%input%geo%species)))
      coords(:,:) = wp%input%geo%coords
      species(:) = wp%input%geo%species

      deallocate(wp%input%geo%coords)
      deallocate(wp%input%geo%species)
      allocate(wp%input%geo%coords(3, nBox * wp%input%geo%nAtom))
      allocate(wp%input%geo%species(nBox * wp%input%geo%nAtom))
      ind = 0
      do i1 = 0, wp%opt%repeatBox(1) - 1
        do i2 = 0, wp%opt%repeatBox(2) - 1
          do i3 = 0, wp%opt%repeatBox(3) - 1
            shift(:) = matmul(wp%opt%boxVecs, real([i1, i2, i3], dp))
            do iAtom = 1, wp%input%geo%nAtom
              wp%input%geo%coords(:,ind+iAtom) = coords(:,iAtom) + shift
            end do
            wp%input%geo%species(ind+1:ind+wp%input%geo%nAtom) = species
            ind = ind + wp%input%geo%nAtom
          end do
        end do
      end do
      wp%input%geo%nAtom = nBox * wp%input%geo%nAtom
    end if
    do i1 = 1, 3
      wp%opt%boxVecs(:,i1) = wp%opt%boxVecs(:,i1) * real(wp%opt%repeatBox(i1), dp)
    end do
  end subroutine expandToSupercell
  
  !> Fills the plotted region with atoms from periodic images of the unit cell.
  subroutine fillPlottedRegion(wp)
    type(TProgramVariables), intent(inout) :: wp
    real(dp) :: invBoxVecs(3,3), recVecs2pi(3,3)
    real(dp) :: cellMiddle(3), boxMiddle(3), frac(3), cubeCorner(3), coord(3), shift(3)
    real(dp) :: mDist, dist
    integer :: i1, i2, i3, iCell, iAtom
    real(dp), allocatable :: fCellVec(:,:), rCellVec(:,:)
    type(TListRealR1) :: coordList
    type(TListInt) :: speciesList

    if (wp%opt%beVerbose) then
      write(stdOut, "(A)") "Filling plotted region with atoms from periodic images."
    end if

    ! Inverse box vectors and reciprocal lattice vectors of the unit cell
    call invert33(invBoxVecs, wp%opt%boxVecs)
    call invert33(recVecs2pi, wp%input%geo%latVecs)
    recVecs2pi = reshape(recVecs2pi, [3, 3], order=[2, 1])

    ! Shifting plotted region by integer lattice vectors, to have its center as close to the center
    ! of the lattice unit cell as possible.
    cellMiddle(:) = 0.5_dp * sum(wp%input%geo%latVecs, dim=2)
    boxMiddle(:) = wp%opt%origin + 0.5_dp * sum(wp%opt%boxVecs, dim=2)
    ! Workaround for intel 2021 ICE, replacing matmul(boxMiddle - cellMiddle, recVecs2pi)
    shift(:) = boxMiddle - cellMiddle
    frac(:) = matmul(shift, recVecs2pi)
    wp%opt%origin(:) = wp%opt%origin - matmul(wp%input%geo%latVecs, real(anint(frac), dp))
    wp%opt%gridOrigin(:) = wp%opt%gridOrigin - matmul(wp%input%geo%latVecs, real(anint(frac), dp))

    ! Determine how many unit cells to include
    ! We need all cells around, which could contain atoms in the sphere, drawn from the center of
    ! the unit cell, containing the entire plotted region
    mDist = 0.0_dp
    do i1 = 0, 1
      do i2 = 0, 1
        do i3 = 0, 1
          cubeCorner(:) = wp%opt%origin + matmul(wp%opt%boxVecs, real([i1, i2, i3], dp))
          dist = norm2(cubeCorner - cellMiddle)
          mDist = max(dist, mDist)
        end do
      end do
    end do
    ! Get all translation vectors that fall within this distance
    call getCellTranslations(fCellVec, rCellVec, wp%input%geo%latVecs, recVecs2pi, mDist)
  
    ! Loop over all atoms in the shifted cells and include them, if they fall in the plotted region
    call init(coordList)
    call init(speciesList)
    do iCell = 1, size(rCellVec, dim=2)
      do iAtom = 1, wp%input%geo%nAtom
        coord(:) = wp%input%geo%coords(:,iAtom) + rCellVec(:,iCell)
        frac(:) = matmul(invBoxVecs, coord - wp%opt%origin)
        if (all(frac > -1e-04_dp) .and. all(frac < 1.0_dp + 1e-04_dp)) then
          call append(coordList, coord)
          call append(speciesList, wp%input%geo%species(iAtom))
        end if
      end do
    end do

    ! Replace coordinates and species in wp with the new ones
    deallocate(wp%input%geo%coords)
    deallocate(wp%input%geo%species)
    wp%input%geo%nAtom = len(coordList)
    allocate(wp%input%geo%coords(3, wp%input%geo%nAtom))
    allocate(wp%input%geo%species(wp%input%geo%nAtom))
    call asArray(coordList, wp%input%geo%coords)
    call asArray(speciesList, wp%input%geo%species)
    deallocate(fCellVec)
    deallocate(rCellVec)
  end subroutine fillPlottedRegion

  

  !> Writes a 3D function as cube file.
  subroutine writeCubeFile(geo, atomicNumbers, gridVecs, origin, gridVal, fileName, comments,&
      & repeatBox, ioStat)

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

    !> Error status
    integer, intent(out), optional :: ioStat

    integer, parameter :: bufferSize = 6
    real(dp) :: buffer(bufferSize)
    character(len=*), parameter :: formBuffer = "(6E16.8)"
    integer :: rep(3)
    integer :: ii, i1, i2, i3, ir1, ir2, ir3
    type(TFileDescr) :: fd
    integer :: ioStat_

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
      rep(:) = repeatBox
    else
      rep(:) = [1, 1, 1]
    end if

    call openFile(fd, fileName, mode="w", iostat=ioStat_, parallelWriting=.true.)

    ! hand over the error status, if asked for
    if (present(ioStat)) ioStat = ioStat_

    if (ioStat_ /= 0) then
      call warning("Error while opening file '" // trim(fileName) // "'.")
      return
    end if
    if (present(comments)) then
      write(fd%unit, "(A)") trim(comments(1))
      write(fd%unit, "(A)") trim(comments(2))
    else
      write(fd%unit, "(A)") "Made by waveplot"
      write(fd%unit, *)
    end if
    write(fd%unit,"(I5,3F12.6)") geo%nAtom, origin
    write(fd%unit,"(I5,3F12.6)") rep(1) * size(gridVal, dim=1), gridVecs(:,1)
    write(fd%unit,"(I5,3F12.6)") rep(2) * size(gridVal, dim=2), gridVecs(:,2)
    write(fd%unit,"(I5,3F12.6)") rep(3) * size(gridVal, dim=3), gridVecs(:,3)
    do ii = 1, geo%nAtom
      write(fd%unit, "(I5,4F12.6)") atomicNumbers(geo%species(ii)), 0.0_dp, geo%coords(:, ii)
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
                  write(fd%unit, formBuffer) buffer
                end if
              end do
              if (ii /= bufferSize) then
                write(fd%unit, "(" // i2c(ii) // "E16.8)") buffer(:ii)
              end if
            end do
          end do
        end do
      end do
    end do

    call closeFile(fd)

  end subroutine writeCubeFile



  !> Helper routine to write a data grid to a cube file with standardized name and comments.
  !! Error handling depends on the presence of ioStat.
  !! If not passed and the writing fails, an error will be raised.
  subroutine writePropertyToCube(plotType, data, wp, levelIndex, ioStat)
    !> Determines file suffix and comment
    character(len=*), intent(in) :: plotType
    real(dp), intent(in) :: data(:,:,:)
    type(TProgramVariables), intent(in) :: wp
    !> Optional tuple of (level, kpoint, spin)
    integer, intent(in), optional :: levelIndex(3)
    integer, intent(out), optional :: ioStat

    character(len=128) :: comments(2), fileName, fileSuffix, commentSuffix

    integer :: iSpin, iKPoint, iLevel, ioStatTmp

    comments(1) = "Cube file generated by WAVEPLOT from data created by DFTB+"

    select case (trim(plotType))
    case ("charge")          ; fileSuffix = "-abs2.cube"       ; commentSuffix = ", abs2"
    case ("chargediff")      ; fileSuffix = "-abs2diff.cube"   ; commentSuffix = ", abs2diff"
    case ("real")            ; fileSuffix = "-real.cube"       ; commentSuffix = ", real"
    case ("imag")            ; fileSuffix = "-imag.cube"       ; commentSuffix = ", imag"
    case ("atomdens")        ; fileSuffix = "wp-atomdens.cube" ; commentSuffix = ", atomdens"
    case ("total_charge")    ; fileSuffix = "wp-abs2.cube"     ; commentSuffix = ", abs2"
    case ("total_chargediff"); fileSuffix = "wp-abs2diff.cube" ; commentSuffix = ", abs2diff"
    case ("spinpol")         ; fileSuffix = "wp-spinpol.cube"  ; commentSuffix = ", spinpol"
    case default
      call error("Unknown plotType in write_property_to_cube: " // trim(plotType))
    end select

    if (present(levelIndex)) then
      iLevel = levelIndex(1); iKPoint = levelIndex(2); iSpin = levelIndex(3)
      fileName = "wp-" // i2c(iSpin) // "-" // i2c(iKPoint) // "-" // i2c(iLevel) // trim(fileSuffix)
      write(comments(2), "('Calc-Id:',I11,', Spin:',I2,', K-Point:',I6,', State:',I6, A)") &
        & wp%input%identity, iSpin, iKPoint, iLevel, trim(commentSuffix)
    else
      fileName = trim(fileSuffix)
      write(comments(2), "('Calc-Id:',I11, A)") wp%input%identity, trim(commentSuffix)
    end if

    call writeCubeFile(wp%input%geo, wp%aNr%atomicNumbers, wp%loc%gridVec, wp%opt%gridOrigin, &
        & data, fileName, comments=comments, repeatBox=wp%opt%repeatBox, ioStat=ioStatTmp)

    if(present(ioStat)) then
        ioStat = ioStatTmp
    end if

    if (ioStatTmp == 0) then
      if (env%tGlobalLead) write(stdOut, "(A)") "File '" // trim(fileName) // "' written"
    else ! either error or warn depending on ioStat presence
      if (present(ioStat)) then
        call warning("Error while writing file '" // trim(fileName) // "'.")
      else
        call error("Error while writing file '" // trim(fileName) // "'.")
      end if
    end if
    

  end subroutine writePropertyToCube

end program waveplot
