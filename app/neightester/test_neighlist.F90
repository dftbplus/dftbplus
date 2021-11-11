module test_neighlist
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment, TEnvironment_init
  use dftbp_common_globalenv, only : initGlobalEnv, destructGlobalEnv
  use dftbp_dftb_neighbouriter, only : TNeighbourIter, TNeighbourIterFact
  use dftbp_dftb_neighbourmap, only : TNeighbourMap
  use dftbp_dftb_repcont, only : TRepCont, getEnergy
  use dftbp_dftb_splinepolyrep, only : TSplinePolyRep
  use dftbp_dftb_staticneighiter, only : TStaticNeighIterFact, TStaticNeighIterFact_init
  use dftbp_dftb_staticneighmap, only : TStaticNeighMap, TStaticNeighMap_init
  use dftbp_dftbplus_hsdhelpers, only : parseHsdInput
  use dftbp_dftbplus_initprogram, only : TDftbPlusMain
  use dftbp_dftbplus_inputdata, only : TInputData
  use dftbp_dftbplus_main, only : runDftbPlus, handleCoordinateChange, handleLatticeChange
  use dftbp_dftb_periodic, only : TNeighbourList, TNeighbourList_init, updateNeighbourList,&
      & getCellTranslations
  use dftbp_math_lapackroutines, only : matinv
  implicit none

  integer, parameter :: nRepeat = 1000
  real(dp), parameter :: cutoff = 20.0_dp
  integer, parameter :: chunkSize = 1000

contains

  subroutine main()

    type(TEnvironment) :: env
    type(TDftbPlusMain), target :: dpmain
    type(TNeighbourList), target :: neighList
    real(dp), allocatable :: coords(:,:)
    real(dp) :: resArr, resArrSum, resIter, resIterSum, resMap, resMapSum

    type(TStaticNeighIterFact), allocatable :: staticNeighIterFact
    class(TNeighbourIterFact), allocatable :: neighIterFact
    type(TStaticNeighMap), allocatable :: staticNeighMap
    class(TNeighbourMap), allocatable :: neighMap
    type(TRepCont), pointer :: pRepCont

    real :: tArr1, tArr2, tArrSum, tIter1, tIter2, tIterSum, tMap1, tMap2, tMapSum
    integer :: nAtom
    integer :: iRep, errStatus

    call initialize(env, dpmain)

    nAtom = size(dpmain%coord0, dim=2)

    call handleLatticeChange(dpmain%latVec, dpmain%scc, dpmain%tblite, dpmain%tStress, dpmain%extPressure,&
        & dpmain%cutOff%mCutoff, dpmain%dispersion, dpmain%solvation, dpmain%cm5Cont, dpmain%recVec,&
        & dpmain%invLatVec, dpmain%cellVol, dpmain%recCellVol, dpmain%extLatDerivs, dpmain%cellVec,&
        & dpmain%rCellVec)

    call handleCoordinateChange(env, dpmain%coord0, dpmain%latVec, dpmain%invLatVec, dpmain%species0,&
        & dpmain%cutoff, dpmain%orb, dpmain%tPeriodic, dpmain%tHelical, dpmain%scc, dpmain%tblite,&
        & dpmain%repulsive, dpmain%dispersion,dpmain%solvation, dpmain%thirdOrd, dpmain%rangeSep,&
        & dpmain%reks, dpmain%img2CentCell, dpmain%iCellVec, dpmain%neighbourList, dpmain%nAllAtom,&
        & dpmain%coord0Fold, dpmain%coord,dpmain%species, dpmain%rCellVec, dpmain%nNeighbourSk,&
        & dpmain%nNeighbourLC, dpmain%ints, dpmain%H0, dpmain%rhoPrim, dpmain%iRhoPrim,&
        & dpmain%ERhoPrim, dpmain%iSparseStart, dpmain%cm5Cont, errStatus)


    allocate(staticNeighIterFact)
    call TStaticNeighIterFact_init(staticNeighIterFact, dpmain%neighbourList)
    call move_alloc(staticNeighIterFact, neighIterFact)

    allocate(staticNeighMap)
    call TStaticNeighMap_init(staticNeighMap, dpmain%neighbourList)
    call move_alloc(staticNeighMap, neighMap)

    select type (rep => dpmain%repulsive)
    type is (TSplinePolyRep)
      pRepCont => rep%twoBodyCont_
    class default
      error stop "Can not repulsive handle class"
    end select

    tArrSum = 0.0_dp
    tIterSum = 0.0_dp
    tMapSum = 0.0_dp
    resArrSum = 0.0_dp
    resIterSum = 0.0_dp
    resMapSum = 0.0_dp

    do iRep = 1, nRepeat
      call cpu_time(tArr1)
      call iterateOverNeighArray(dpmain%neighbourList, pRepCont, dpmain%img2CentCell, dpmain%coord,&
          & dpmain%species, resArr)
      call cpu_time(tArr2)
      resArrSum = resArrSum + resArr
      tArrSum = tArrSum + (tArr2 - tArr1)

      call cpu_time(tIter1)
      call iterateOverNeighIter(neighIterFact, pRepCont, dpmain%coord, dpmain%species,&
          & dpmain%img2CentCell, cutoff, nAtom, chunkSize, resIter)
      call cpu_time(tIter2)
      resIterSum = resIterSum + resIter
      tIterSum = tIterSum + (tIter2 - tIter1)

      call cpu_time(tMap1)
      call iterateOverNeighMap(neighMap, pRepCont, dpmain%coord, dpmain%species, cutoff, nAtom,&
          & resMap)
      call cpu_time(tMap2)
      resMapSum = resMapSum + resMap
      tMapSum = tMapSum + (tMap2 - tMap1)
    end do

    print *, "CUTOFF:", cutoff
    print *, "REPETITIONS:", nRepeat
    print *, "RES ARRAY:", resArrSum, tArrSum / real(nRepeat, dp)
    print *, "RES ITER :", resIterSum, tIterSum / real(nRepeat, dp)
    print *, "RES Map :", resMapSum, tMapSum / real(nRepeat, dp)
    print *, "TIMEDIFF/ITER:", (tIterSum - tArrSum) / real(nRepeat, dp),&
        & (tMapSum - tArrSum) / real(nRepeat, dp)

    call finalize(env, dpmain)

  end subroutine main


  subroutine initialize(env, main)
    type(TEnvironment), intent(out) :: env
    type(TDftbPlusMain), intent(out) :: main

    type(TInputData), allocatable :: input

    call initGlobalEnv()
    allocate(input)
    call parseHsdInput(input)
    call TEnvironment_init(env)
    call main%initProgramVariables(input, env)
    deallocate(input)

  end subroutine initialize


  subroutine finalize(env, main)
    type(TEnvironment), intent(inout) :: env
    type(TDftbPlusMain), intent(inout) :: main

    call main%destructProgramVariables()
    call env%destruct()
    call destructGlobalEnv()

  end subroutine finalize


  subroutine buildNeighList(coord0, latVecs, cutoff, neighList, coords)
    real(dp), intent(in) :: coord0(:,:), latVecs(:,:)
    real(dp), intent(in) :: cutoff
    type(TNeighbourList), intent(out) :: neighList
    real(dp), allocatable, intent(out) :: coords(:,:)

    real(dp), allocatable :: cellVecs(:,:), rCellVecs(:,:)
    integer, allocatable :: img2CentCell(:), iCellVec(:)
    real(dp) :: recVecs2p(3, 3)
    integer :: nAtom, nAllAtom

    recVecs2p(:,:) = latVecs
    call matinv(recVecs2p)
    recVecs2p = transpose(recVecs2p)
    call getCellTranslations(cellVecs, rCellVecs, latVecs, recVecs2p, cutOff)
    nAtom = size(coord0, dim=2)
    allocate(coords(3, nAtom), img2CentCell(nAtom), iCellVec(nAtom))
    call TNeighbourList_init(neighList, size(coord0, dim=2), 100)
    call updateNeighbourList(coords, img2CentCell, iCellVec, neighList, nAllAtom, coord0, cutoff,&
      & rCellVecs)

  end subroutine buildNeighList


  subroutine iterateOverNeighArray(neighList, repCont, img2CentCell, coords, species, res)
    type(TNeighbourList), intent(in) :: neighList
    type(TRepCont), intent(in) :: repCont
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    real(dp), intent(out) :: res

    integer :: nAtom, iAt1, iAt2, iAt2f, iNeigh
    real(dp) :: resAtom(size(neighlist%nNeighbour))
    real(dp) :: vect(3), dist, intermed

    nAtom = size(neighList%nNeighbour)
    resAtom(:) = 0.0_dp
    !$omp parallel do default(none) private(iAt1, iNeigh, iAt2, iAt2f, intermed, vect, dist)&
    !$omp& shared(neighList, repCont, img2CentCell, coords, species, resAtom, nAtom)
    do iAt1 = 1, nAtom
      do iNeigh = 1, min(4, neighList%nNeighbour(iAt1))
        iAt2 = neighList%iNeighbour(iNeigh, iAt1)
        vect(:) = coords(:,iAt1) - coords(:,iAt2)
        dist = sqrt(sum(vect**2))
        call getEnergy(repCont, intermed, dist, species(iAt1), species(iAt2))
        resAtom(iAt1) = resAtom(iAt1) + 0.5_dp * intermed
        !if (iAt2f /= iAt1) then
        !  resAtom(iAt2f) = resAtom(iAt2f) + 0.5_dp * intermed
        !end if
        !res = res + 1.0_dp / sqrt(neighList%neighDist2(iNeigh, iAt1))
      end do
    end do
    !$omp end parallel do
    res = sum(resAtom)
  end subroutine iterateOverNeighArray


  subroutine iterateOverNeighIter(neighIterFact, repCont, coords, species, img2CentCell, cutoff, nAtom, chunkSize, res)
    class(TNeighbourIterFact), intent(in) :: neighIterFact
    type(TRepCont), intent(in) :: repCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:), img2CentCell(:)
    real(dp), intent(in) :: cutoff
    integer, intent(in) :: nAtom
    integer, intent(in) :: chunkSize
    real(dp), intent(out) :: res

    class(TNeighbourIter), allocatable :: neighIter
    integer :: iNeighs(chunkSize)
    real(dp) :: distances2(chunkSize), vect(3), dist, intermed
    integer :: iAt1, iAt2, iAt2f, iNeigh, nNeigh
    real(dp) :: resAtom(nAtom)

    resAtom(:) = 0.0_dp
    !$omp parallel default(none)&
    !$omp& firstprivate(neighIter)&
    !$omp& private(iAt1, nNeigh, iNeigh, iAt2, iAt2f, dist, iNeighs, distances2, vect, intermed)&
    !$omp& shared(cutoff, nAtom, chunkSize, neighIterFact, repCont, resAtom, img2CentCell, coords, species)
    call neighIterFact%getIterator(cutoff, neighIter)

    !$omp do
    do iAt1 = 1, nAtom
      call neighIter%start(iAt1, cutoff, .false., chunkSize)
      do
        call neighIter%get(nNeigh, iNeighs=iNeighs, distances2=distances2)
        do iNeigh = 1, min(4, nNeigh)
          iAt2 = iNeighs(iNeigh)
          vect(:) = coords(:,iAt1) - coords(:,iAt2)
          dist = sqrt(sum(vect**2))
          call getEnergy(repCont, intermed, dist, species(iAt1), species(iAt2))
          resAtom(iAt1) = resAtom(iAt1) + 0.5_dp * intermed
          !if (iAt2f /= iAt1) then
          !  resAtom(iAt2f) = resAtom(iAt2f) + 0.5_dp * intermed
          !end if
          !res = res + 1.0_dp / sqrt(distances2(iNeigh))
        end do
        exit
        !if (nNeigh < chunkSize) exit
      end do
    end do
    !$omp end do
    deallocate(neighIter)
    !$omp end parallel
    res = sum(resAtom)

  end subroutine iterateOverNeighIter


  subroutine iterateOverNeighMap(neighMap, repCont, coords, species, cutoff, nAtom, res)
    class(TNeighbourMap), intent(in) :: neighMap
    type(TRepCont), intent(in) :: repCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    real(dp), intent(in) :: cutoff
    integer, intent(in) :: nAtom
    real(dp), intent(out) :: res

    integer, allocatable :: iNeighs(:)
    real(dp), allocatable :: distances2(:)
    integer :: iAt1, iAt2, iNeigh
    real(dp) :: vect(3), dist, intermed
    real(dp) :: resAtom(nAtom)

    resAtom(:) = 0.0_dp

    !$omp parallel do default(none)&
    !$omp& firstprivate(distances2, iNeighs)&
    !$omp& private(iAt1, iNeigh, iAt2, vect, dist, intermed)&
    !$omp& shared(nAtom, cutoff, neighMap, coords, species, repCont, resAtom)
    do iAt1 = 1, nAtom
      call neighMap%getNeighbours(iAt1, cutoff, includeSelf=.false., iNeighs=iNeighs,&
          & distances2=distances2)
      do iNeigh = 1, min(4, size(iNeighs))
          iAt2 = iNeighs(iNeigh)
          vect(:) = coords(:,iAt1) - coords(:,iAt2)
          dist = sqrt(sum(vect**2))
          call getEnergy(repCont, intermed, dist, species(iAt1), species(iAt2))
          resAtom(iAt1) = resAtom(iAt1) + 0.5_dp * intermed
          !res = res + 1.0_dp / sqrt(distances2(iNeigh))
      end do
    end do
    !$omp end parallel do
    res = sum(resAtom)

  end subroutine iterateOverNeighMap


end module test_neighlist


program test_neighlist_program
  use test_neighlist, only : main
  implicit none

  call main()

end program test_neighlist_program