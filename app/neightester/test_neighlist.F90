#:include "common.fypp"

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
#:if WITH_OMP
  use omp_lib, only : omp_get_max_threads
#:endif
  implicit none

  integer, parameter :: maxAtom = 512
  integer, parameter :: maxNeigh = 1000
  integer, parameter :: nRepeat = 10000
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
    integer :: cArr1, cArr2, cIter1, cIter2, cMap1, cMap2, countRate
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

    resArrSum = 0.0_dp
    call system_clock(count=cArr1, count_rate=countRate)
    do iRep = 1, nRepeat
      call iterateOverNeighArray(dpmain%neighbourList, pRepCont, dpmain%img2CentCell, dpmain%coord,&
          & dpmain%species, resArr)
      resArrSum = resArrSum + resArr
    end do
    call system_clock(count=cArr2)
    tArrSum = real(cArr2 - cArr1, dp) / real(countRate, dp)

    resIterSum = 0.0_dp
    call system_clock(count=cIter1)
    do iRep = 1, nRepeat
      call iterateOverNeighIter(neighIterFact, pRepCont, dpmain%coord, dpmain%species,&
          & dpmain%img2CentCell, cutoff, nAtom, chunkSize, resIter)
      resIterSum = resIterSum + resIter
    end do
    call system_clock(count=cIter2)
    tIterSum = real(cIter2 - cIter1, dp) / real(countRate, dp)

    resMapSum = 0.0_dp
    call system_clock(count=cMap1)
    do iRep = 1, nRepeat
      call iterateOverNeighMap(neighMap, pRepCont, dpmain%coord, dpmain%species, cutoff, nAtom,&
          & resMap)
      resMapSum = resMapSum + resMap
    end do
    call system_clock(count=cMap2)
    tMapSum = real(cMap2 - cMap1, dp) / real(countRate, dp)

    print "(a, 1x, f5.2)", "CUTOFF:", cutoff
    print "(a, 1x, i0)", "MAX ATOM:", maxAtom
    print "(a, 1x, i0)", "MAX NEIGHS:", maxNeigh
    print "(a, 1x, i0)", "REPETITIONS:", nRepeat
    #:if WITH_OMP
      print "(a, 1x, i0)", "OMP_THREADS:", omp_get_max_threads()
    #:endif
    print "(a, 1x, es16.8, 1x, es10.2, 1x, es10.2, 1x, es10.2, 1x, es10.2)", "ARRAY", resArrSum,&
        & tArrSum, tArrSum / real(nRepeat, dp), 0.0_dp, 0.0_dp
    print "(a, 1x, es16.8, 1x, es10.2, 1x, es10.2, 1x, es10.2, 1x, es10.2)", "ITER ", resIterSum,&
        & tIterSum, tIterSum / real(nRepeat, dp),  (tIterSum - tArrSum) / real(nRepeat, dp),&
        & (tIterSum - tArrSum) / (real(nRepeat, dp) * min(maxAtom, nAtom))
    print "(a, 1x, es16.8, 1x, es10.2, 1x, es10.2, 1x, es10.2, 1x, es10.2)", "MAP  ", resMapSum,&
        & tMapSum, tMapSum / real(nRepeat, dp), (tMapSum - tArrSum) / real(nRepeat, dp),&
        & (tMapSum - tArrSum) / (real(nRepeat, dp) * min(maxAtom, nAtom))

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
    do iAt1 = 1, min(maxAtom, nAtom)
      do iNeigh = 1, min(maxNeigh, neighList%nNeighbour(iAt1))
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
    do iAt1 = 1, min(maxAtom, nAtom)
      call neighIter%start(iAt1, cutoff, .false., chunkSize)
      do
        call neighIter%get(nNeigh, iNeighs=iNeighs, distances2=distances2)
        do iNeigh = 1, min(maxNeigh, nNeigh)
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
    do iAt1 = 1, min(maxAtom, nAtom)
      call neighMap%getNeighbours(iAt1, cutoff, includeSelf=.false., iNeighs=iNeighs,&
          & distances2=distances2)
      do iNeigh = 1, min(maxNeigh, size(iNeighs))
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