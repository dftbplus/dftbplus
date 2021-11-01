module test_neighlist
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment, TEnvironment_init
  use dftbp_common_globalenv, only : initGlobalEnv, destructGlobalEnv
  use dftbp_dftb_neighbouriter, only : TNeighbourIter, TNeighbourIterFact
  use dftbp_dftb_staticneighiter, only : TStaticNeighIterFact, TStaticNeighIterFact_init
  use dftbp_dftbplus_hsdhelpers, only : parseHsdInput
  use dftbp_dftbplus_initprogram, only : TDftbPlusMain
  use dftbp_dftbplus_inputdata, only : TInputData
  use dftbp_dftbplus_main, only : runDftbPlus
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
    type(TDftbPlusMain) :: dpmain
    type(TNeighbourList), target :: neighList
    real(dp), allocatable :: coords(:,:)
    real(dp) :: resArr, resArrSum, resIter, resIterSum

    type(TStaticNeighIterFact), allocatable :: staticNeighIterFact
    class(TNeighbourIterFact), allocatable :: neighIterFact
    real :: tArr1, tArr2, tArrSum, tIter1, tIter2, tIterSum
    integer :: nAtom
    integer :: iRep

    call initialize(env, dpmain)

    nAtom = size(dpmain%coord0, dim=2)
    call buildNeighList(dpmain%coord0, dpmain%latVec, cutoff, neighList, coords)

    allocate(staticNeighIterFact)
    call TStaticNeighIterFact_init(staticNeighIterFact, neighList)
    call move_alloc(staticNeighIterFact, neighIterFact)

    tArrSum = 0.0_dp
    tIterSum = 0.0_dp
    resArrSum = 0.0_dp
    resIterSum = 0.0_dp

    do iRep = 1, nRepeat
      call cpu_time(tArr1)
      call iterateOverNeighArray(neighList, coords, resArr)
      call cpu_time(tArr2)
      resArrSum = resArrSum + resArr
      tArrSum = tArrSum + (tArr2 - tArr1)

      call cpu_time(tIter1)
      call iterateOverNeighIter(neighIterFact, cutoff, coords, nAtom, chunkSize, resIter)
      call cpu_time(tIter2)
      resIterSum = resIterSum + resIter
      tIterSum = tIterSum + (tIter2 - tIter1)
    end do

    print *, "CUTOFF:", cutoff
    print *, "REPETITIONS:", nRepeat
    print *, "RES ARRAY:", resArrSum, tArrSum
    print *, "RES ITER :", resIterSum, tIterSum
    print *, "TIMEDIFF/ATOM:", (tIterSum - tArrSum) / real(nRepeat * nAtom, dp)

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


  subroutine iterateOverNeighArray(neighList, coords, res)
    type(TNeighbourList), intent(in) :: neighList
    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(out) :: res

    integer :: nAtom, iAt1, iAt2, iNeigh

    nAtom = size(neighList%nNeighbour)
    res = 0.0_dp
    do iAt1 = 1, nAtom
      do iNeigh = 1, neighList%nNeighbour(iAt1)
        iAt2 = neighList%iNeighbour(iNeigh, iAt1)
        res = res + 1.0_dp / sqrt(neighList%neighDist2(iNeigh, iAt1))
      end do
    end do

  end subroutine iterateOverNeighArray


  subroutine iterateOverNeighIter(neighIterFact, cutoff, coords, nAtom, chunkSize, res)
    class(TNeighbourIterFact), intent(in) :: neighIterFact
    real(dp), intent(in) :: cutoff
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nAtom
    integer, intent(in) :: chunkSize
    real(dp), intent(out) :: res

    class(TNeighbourIter), allocatable :: neighIter
    integer :: iNeighs(chunkSize)
    real(dp) :: distances2(chunkSize)
    integer :: iAt1, iAt2, iNeigh, nNeigh

    res = 0.0_dp
    neighIter = neighIterFact%getIterator(cutoff)
    do iAt1 = 1, nAtom
      call neighIter%start(iAt1, cutoff, .false., chunkSize)
      do
        call neighIter%get(nNeigh, iNeighs=iNeighs, distances2=distances2)
        do iNeigh = 1, nNeigh
          iAt2 = iNeighs(iNeigh)
          res = res + 1.0_dp / sqrt(distances2(iNeigh))
        end do
        if (nNeigh < chunkSize) exit
      end do
    end do

  end subroutine iterateOverNeighIter


end module test_neighlist


program test_neighlist_program
  use test_neighlist, only : main
  implicit none

  call main()

end program test_neighlist_program