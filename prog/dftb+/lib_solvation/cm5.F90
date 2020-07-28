!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implementation of the charge model 5 (CM5)
module dftbp_cm5
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_atomicrad, only : getAtomicRad
  use dftbp_blasroutines, only : gemv
  use dftbp_constants, only : AA__Bohr, symbolToNumber
  use dftbp_periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use dftbp_simplealgebra, only : determinant33
  implicit none
  private

  public :: TChargeModel5, TCM5Input, TChargeModel5_init


  !> Charge model 5 input data
  type :: TCM5Input

    !> Real space cutoff
    real(dp) :: rCutoff

    !> Global parameter
    real(dp) :: alpha

    !> Atomic radii
    real(dp), allocatable :: atomicRad(:)

  end type TCM5Input


  !> Charge model 5 class
  type :: TChargeModel5
    private

    !> number of atoms
    integer :: nAtom

    !> lattice vectors if periodic
    real(dp) :: latVecs(3, 3)

    !> is this periodic
    logical :: tPeriodic

    !> are the coordinates current?
    logical :: tCoordsUpdated

    !> Real space cutoff
    real(dp) :: rCutoff

    !> Global parameter
    real(dp) :: alpha

    !> Atomic radii
    real(dp), allocatable :: atomicRad(:)

    !> Pair parameters
    real(dp), allocatable :: pairParam(:, :)

    !> CM5 correction
    real(dp), allocatable, public :: cm5(:)

    !> Derivative of CM5 correction w.r.t. coordinates
    real(dp), allocatable, public :: dcm5dr(:, :, :)

    !> Derivative of CM5 correction w.r.t. strain deformations
    real(dp), allocatable, public :: dcm5dL(:, :, :)

  contains

    !> update internal copy of coordinates
    procedure :: updateCoords

    !> update internal copy of lattice vectors
    procedure :: updateLatVecs

    !> get real space cutoff
    procedure :: getRCutoff

    !> get charge contributions
    procedure :: addCharges

    !> get force contributions
    procedure :: addGradients

    !> get stress tensor contributions
    procedure :: addSigma

  end type TChargeModel5


  !> Get pair parameters for two species (Dzz')
  interface getPairParameter
    module procedure :: getPairParameterSymbol
    module procedure :: getPairParameterNumber
  end interface getPairParameter


  ! Charge model 5 atomwise parameters (Dz in paper)
  real(dp), parameter :: pairParameters(1:118) = [ &
      & 0.0056_dp,-0.1543_dp, 0.0000_dp, 0.0333_dp,-0.1030_dp,-0.0446_dp, &
      &-0.1072_dp,-0.0802_dp,-0.0629_dp,-0.1088_dp, 0.0184_dp, 0.0000_dp, &
      &-0.0726_dp,-0.0790_dp,-0.0756_dp,-0.0565_dp,-0.0444_dp,-0.0767_dp, &
      & 0.0130_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
      & 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
      &-0.0512_dp,-0.0557_dp,-0.0533_dp,-0.0399_dp,-0.0313_dp,-0.0541_dp, &
      & 0.0092_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
      & 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
      &-0.0361_dp,-0.0393_dp,-0.0376_dp,-0.0281_dp,-0.0220_dp,-0.0381_dp, &
      & 0.0065_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
      & 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
      & 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
      & 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
      & 0.0000_dp, 0.0000_dp,-0.0255_dp,-0.0277_dp,-0.0265_dp,-0.0198_dp, &
      &-0.0155_dp,-0.0269_dp, 0.0046_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
      & 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
      & 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
      & 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
      & 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp,-0.0179_dp,-0.0195_dp, &
      &-0.0187_dp,-0.0140_dp,-0.0110_dp,-0.0189_dp]


contains


  !> Initialize generalized charge model 5 from geometry
  subroutine TChargeModel5_init(this, input, nAtom, speciesNames, tDerivs, latVecs)

    !> Initialised instance at return
    type(TChargeModel5), intent(out) :: this

    !> Charge model 5 input data
    type(TCM5Input), intent(in) :: input

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    !> Setup container to evaluate derivatives
    logical, intent(in) :: tDerivs

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    integer :: nSpecies
    integer :: iSp1, iSp2

    @:ASSERT(allocated(input%atomicRad))

    nSpecies = size(speciesNames)
    this%tPeriodic = present(latVecs)
    if (this%tPeriodic) then
      call this%updateLatVecs(LatVecs)
    end if
    this%nAtom = nAtom

    allocate(this%cm5(nAtom))
    if (tDerivs) then
      allocate(this%dcm5dr(3, nAtom, nAtom))
      allocate(this%dcm5dL(3, 3, nAtom))
    end if

    this%rCutoff = input%rCutoff

    allocate(this%atomicRad(nSpecies))
    allocate(this%pairParam(nSpecies, nSpecies))
    this%atomicRad(:) = input%atomicRad

    this%alpha = input%alpha
    do iSp1 = 1, nSpecies
      do iSp2 = 1, nSpecies
        this%pairParam(iSp1, iSp2) = getPairParameter(speciesNames(iSp1), speciesnames(iSp2))
      end do
      this%pairParam(iSp1, iSp1) = 0.0_dp
    end do

    this%tCoordsUpdated = .false.

  end subroutine TChargeModel5_init


  !> Update internal stored coordinates
  subroutine updateCoords(this, neighList, img2CentCell, coords, species0)

    !> data structure
    class(TChargeModel5), intent(inout) :: this

    !> list of neighbours to atoms
    type(TNeighbourList), intent(in) :: neighList

    !> image to central cell atom index
    integer, intent(in) :: img2CentCell(:)

    !> atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    integer, allocatable :: nNeigh(:)

    allocate(nNeigh(this%nAtom))
    call getNrOfNeighboursForAll(nNeigh, neighList, this%rCutoff)
    if (allocated(this%dcm5dr) .and. allocated(this%dcm5dL)) then
      call getCorrectionDerivs(this, nNeigh, neighList%iNeighbour, img2CentCell, &
          & neighList%neighDist2, species0, coords)
    else
      call getCorrection(this, nNeigh, neighList%iNeighbour, img2CentCell, &
          & neighList%neighDist2, species0, coords)
    end if

    this%tCoordsUpdated = .true.

  end subroutine updateCoords


  !> update internal copy of lattice vectors
  subroutine updateLatVecs(this, latVecs)

    !> data structure
    class(TChargeModel5), intent(inout) :: this

    !> lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    @:ASSERT(this%tPeriodic)
    @:ASSERT(all(shape(latvecs) == shape(this%latvecs)))

    this%latVecs(:,:) = latVecs

    this%tCoordsUpdated = .false.

  end subroutine updateLatVecs


  !> get charge corrections
  subroutine addCharges(this, charges)

    !> data structure
    class(TChargeModel5), intent(inout) :: this

    !> energy contributions for each atom
    real(dp), intent(inout) :: charges(:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(allocated(this%cm5))
    @:ASSERT(size(charges) == this%nAtom)

    charges(:) = charges + this%cm5

  end subroutine addCharges


  !> get force contributions
  subroutine addGradients(this, dEdcm5, gradients)

    !> data structure
    class(TChargeModel5), intent(inout) :: this

    !> Derivative w.r.t. CM5 correction
    real(dp), intent(in) :: dEdcm5(:)

    !> gradient contributions for each atom
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(allocated(this%dcm5dr))
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    call gemv(gradients, this%dcm5dr, dEdcm5, beta=1.0_dp)

  end subroutine addGradients


  !> get stress tensor contributions
  subroutine addSigma(this, dEdcm5, stress)

    !> data structure
    class(TChargeModel5), intent(inout) :: this

    !> Derivative w.r.t. CM5 correction
    real(dp), intent(in) :: dEdcm5(:)

    !> Stress tensor contributions (not volume scaled)
    real(dp), intent(inout) :: stress(:,:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(allocated(this%dcm5dL))
    @:ASSERT(all(shape(stress) == [3, 3]))

    call gemv(stress, this%dcm5dL, dEdcm5, beta=1.0_dp)

  end subroutine addSigma


  !> Distance cut off for charge model
  function getRCutoff(this) result(cutoff)

    !> data structure
    class(TChargeModel5), intent(inout) :: this

    !> resulting cutoff
    real(dp) :: cutoff

    cutoff = this%rCutoff

  end function getRCutoff


  !> Calculate CM5 correction for this geometry
  subroutine getCorrection(this, nNeighbour, iNeighbour, img2CentCell, neighDist2, &
      & species, coords)

    !> data structure
    type(TChargeModel5), intent(inout) :: this

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist
    integer, intent(in) :: iNeighbour(0:, :)

    !> Square distances of the neighbours
    integer, intent(in) :: img2CentCell(:)

    !> Square distances of the neighbours
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> Current atomic positions
    real(dp), intent(in) :: coords(:, :)

    integer :: iAt1, iSp1, iNeigh, iAt2, iAt2f, iSp2
    real(dp) :: dist, dEr, p12, p21

    this%cm5(:) = 0.0_dp

    do iAt1 = 1, this%nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeighbour(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        if (iSp1 == iSp2) cycle  ! includes iAt1 == iAt2f case
        dist = sqrt(neighDist2(iNeigh, iAt1))
        p12 = this%pairParam(iSp1, iSp2)
        p21 = this%pairParam(iSp2, iSp1)

        dEr = exp(-this%alpha*(dist-this%atomicRad(iSp1)-this%atomicRad(iSp2)))

        this%cm5(iAt1) = this%cm5(iAt1) + dEr * p12
        this%cm5(iAt2f) = this%cm5(iAt2f) + dEr * p21

      end do
    end do

  end subroutine getCorrection


  !> Calculate CM5 correction for this geometry
  subroutine getCorrectionDerivs(this, nNeighbour, iNeighbour, img2CentCell, &
      & neighDist2, species, coords)

    !> data structure
    type(TChargeModel5), intent(inout) :: this

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist
    integer, intent(in) :: iNeighbour(0:, :)

    !> Square distances of the neighbours
    integer, intent(in) :: img2CentCell(:)

    !> Square distances of the neighbours
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> Current atomic positions
    real(dp), intent(in) :: coords(:, :)

    integer :: iAt1, iSp1, iNeigh, iAt2, iAt2f, iSp2
    real(dp) :: dist, vec(3), dEr, dGr(3), dSr(3, 3), p12, p21

    this%cm5(:) = 0.0_dp
    this%dcm5dr(:, :, :) = 0.0_dp
    this%dcm5dL(:, :, :) = 0.0_dp

    do iAt1 = 1, this%nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeighbour(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        if (iSp1 == iSp2) cycle  ! includes iAt1 == iAt2f case
        dist = sqrt(neighDist2(iNeigh, iAt1))
        vec(:) = coords(:, iAt1) - coords(:, iAt2)
        p12 = this%pairParam(iSp1, iSp2)
        p21 = this%pairParam(iSp2, iSp1)

        dEr = exp(-this%alpha*(dist-this%atomicRad(iSp1)-this%atomicRad(iSp2)))
        dGr = dEr * this%alpha * vec/dist
        dSr = spread(dGr, 1, 3) * spread(vec, 2, 3)

        this%cm5(iAt1) = this%cm5(iAt1) + dEr * p12
        this%cm5(iAt2f) = this%cm5(iAt2f) + dEr * p21

        this%dcm5dr(:, iAt1, iAt1) = this%dcm5dr(:, iAt1, iAt1) - dGr * p12
        this%dcm5dr(:, iAt2f, iAt2f) = this%dcm5dr(:, iAt2f, iAt2f) + dGr * p21
        this%dcm5dr(:, iAt1, iAt2f) = this%dcm5dr(:, iAt1, iAt2f) - dGr * p21
        this%dcm5dr(:, iAt2f, iAt1) = this%dcm5dr(:, iAt2f, iAt1) + dGr * p12

        this%dcm5dL(:, :, iAt1) = this%dcm5dL(:, :, iAt1) + dSr * p12
        this%dcm5dL(:, :, iAt2f) = this%dcm5dL(:, :, iAt2f) + dSr * p21

      end do
    end do

  end subroutine getCorrectionDerivs


  !> Get pair parameter (Dzz') for species with a given symbols
  elemental function getPairParameterSymbol(symbol1, symbol2) result(pairPar)

    !> Element symbol
    character(len=*), intent(in) :: symbol1

    !> Element symbol
    character(len=*), intent(in) :: symbol2

    !> Pair parameter
    real(dp) :: pairPar

    pairPar = getPairParameter(symbolToNumber(symbol1), symbolToNumber(symbol2))

  end function getPairParameterSymbol


  !> Get pair parameter (Dzz') for species with a given atomic numbers
  elemental function getPairParameterNumber(number1, number2) result(pairPar)

    !> Atomic number
    integer, intent(in) :: number1

    !> Atomic number
    integer, intent(in) :: number2

    !> Pair parameter
    real(dp) :: pairPar

    if (number1 > 0 .and. number1 <= size(pairParameters, dim=1) .and. &
        & number2 > 0 .and. number2 <= size(pairParameters, dim=1)) then
      if (number1 == 1 .and. number2 == 6) then
        pairPar = 0.0502_dp
      else if (number1 == 6 .and. number2 == 1) then
        pairPar = -0.0502_dp
      else if (number1 == 1 .and. number2 == 7) then
        pairPar = 0.1747_dp
      else if (number1 == 7 .and. number2 == 1) then
        pairPar = -0.1747_dp
      else if (number1 == 1 .and. number2 == 8) then
        pairPar = 0.1671_dp
      else if (number1 == 8 .and. number2 == 1) then
        pairPar = -0.1671_dp
      else if (number1 == 6 .and. number2 == 7) then
        pairPar = 0.0556_dp
      else if (number1 == 7 .and. number2 == 6) then
        pairPar = -0.0556_dp
      else if (number1 == 6 .and. number2 == 8) then
        pairPar = 0.0234_dp
      else if (number1 == 8 .and. number2 == 6) then
        pairPar = -0.0234_dp
      else if (number1 == 7 .and. number2 == 8) then
        pairPar = -0.0346_dp
      else if (number1 == 8 .and. number2 == 7) then
        pairPar = 0.0346_dp
      else
        pairPar = pairParameters(number1) - pairParameters(number2)
      end if
    else
      pairPar = 0.0_dp
    end if

  end function getPairParameterNumber


end module dftbp_cm5
