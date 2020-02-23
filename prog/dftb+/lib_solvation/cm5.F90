!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implementation of the charge model 5 (CM5)
module dftbp_cm5
  use dftbp_accuracy, only : dp
  use dftbp_blasroutines, only : gemv
  use dftbp_constants, only : AA__Bohr, symbolToNumber
  use dftbp_periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use dftbp_simplealgebra, only : determinant33
  implicit none
  private

  public :: TChargeModel5, init


  !> Charge model 5 class
  type :: TChargeModel5
    private

    !> number of atoms
    integer :: nAtom = 0

    !> lattice vectors if periodic
    real(dp) :: latVecs(3, 3) = 0.0_dp

    !> is this periodic
    logical :: tPeriodic

    !> are the coordinates current?
    logical :: tCoordsUpdated = .false.

    !> Real space cutoff
    real(dp) :: rCutoff = 0.0_dp

    !> Global parameter
    real(dp) :: alpha = 2.4740_dp / AA__Bohr

    !> Atomic radii
    real(dp), allocatable :: atomicRad(:)

    !> Pair parameters
    real(dp), allocatable :: pairParam(:, :)

    !> CM5 correction
    real(dp), allocatable :: cm5(:)

    !> Derivative of CM5 correction w.r.t. coordinates
    real(dp), allocatable :: dcm5dr(:, :, :)

    !> Derivative of CM5 correction w.r.t. strain deformations
    real(dp), allocatable :: dcm5dL(:, :, :)

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
    procedure :: addStress

  end type TChargeModel5


  !> Initialize charge model 5 from geometry
  interface init
    module procedure :: initialize
  end interface init


  !> Get atomic radius for a species
  interface getAtomicRadius
    module procedure :: getAtomicRadiusSymbol
    module procedure :: getAtomicRadiusNumber
  end interface getAtomicRadius


  !> Covalent radii
  !>
  !> based on "Atomic Radii of the Elements," M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
  !> in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
  !> edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
  !> corrected Nov. 17, 2010 for the 92nd edition.
  real(dp), parameter :: atomicRadii(1:118) = AA__Bohr * [ &
      & 0.32_dp, 0.37_dp, 1.30_dp, 0.99_dp, 0.84_dp, 0.75_dp, 0.71_dp, 0.64_dp, &
      & 0.60_dp, 0.62_dp, 1.60_dp, 1.40_dp, 1.24_dp, 1.14_dp, 1.09_dp, 1.04_dp, &
      & 1.00_dp, 1.01_dp, 2.00_dp, 1.74_dp, 1.59_dp, 1.48_dp, 1.44_dp, 1.30_dp, &
      & 1.29_dp, 1.24_dp, 1.18_dp, 1.17_dp, 1.22_dp, 1.20_dp, 1.23_dp, 1.20_dp, &
      & 1.20_dp, 1.18_dp, 1.17_dp, 1.16_dp, 2.15_dp, 1.90_dp, 1.76_dp, 1.64_dp, &
      & 1.56_dp, 1.46_dp, 1.38_dp, 1.36_dp, 1.34_dp, 1.30_dp, 1.36_dp, 1.40_dp, &
      & 1.42_dp, 1.40_dp, 1.40_dp, 1.37_dp, 1.36_dp, 1.36_dp, 2.38_dp, 2.06_dp, &
      & 1.94_dp, 1.84_dp, 1.90_dp, 1.88_dp, 1.86_dp, 1.85_dp, 1.83_dp, 1.82_dp, &
      & 1.81_dp, 1.80_dp, 1.79_dp, 1.77_dp, 1.77_dp, 1.78_dp, 1.74_dp, 1.64_dp, &
      & 1.58_dp, 1.50_dp, 1.41_dp, 1.36_dp, 1.32_dp, 1.30_dp, 1.30_dp, 1.32_dp, &
      & 1.44_dp, 1.45_dp, 1.50_dp, 1.42_dp, 1.48_dp, 1.46_dp, 2.42_dp, 2.11_dp, &
      & 2.01_dp, 1.90_dp, 1.84_dp, 1.83_dp, 1.80_dp, 1.80_dp, 1.73_dp, 1.68_dp, &
      & 1.68_dp, 1.68_dp, 1.65_dp, 1.67_dp, 1.73_dp, 1.76_dp, 1.61_dp, 1.57_dp, &
      & 1.49_dp, 1.43_dp, 1.41_dp, 1.34_dp, 1.29_dp, 1.28_dp, 1.21_dp, 1.22_dp, &
      & 1.36_dp, 1.43_dp, 1.62_dp, 1.75_dp, 1.65_dp, 1.57_dp]


  !> Get pair parameters for two species
  interface getPairParameter
    module procedure :: getPairParameterSymbol
    module procedure :: getPairParameterNumber
  end interface getPairParameter


  ! Charge model 5 atomwise parameters
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
  subroutine initialize(self, nAtom, species0, speciesNames, rCutoff, latVecs)

    !> Initialised instance at return
    type(TChargeModel5), intent(out) :: self

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    !> Real space cutoff
    real(dp), intent(in) :: rCutoff

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    integer :: iSp1, iSp2

    self%tPeriodic = present(latVecs)
    if (self%tPeriodic) then
      call self%updateLatVecs(LatVecs)
    end if
    self%nAtom = nAtom

    allocate(self%cm5(nAtom))
    allocate(self%dcm5dr(3, nAtom, nAtom))
    allocate(self%dcm5dL(3, 3, nAtom))

    self%rCutoff = rCutoff

    allocate(self%atomicRad(size(species0)))
    allocate(self%pairParam(size(species0), size(species0)))

    do iSp1 = 1, size(species0)
      self%atomicRad(iSp1) = getAtomicRadius(speciesNames(iSp1))
      do iSp2 = 1, size(species0)
        self%pairParam(iSp1, iSp2) = getPairParameter(speciesNames(iSp1), speciesnames(iSp2))
      end do
      self%pairParam(iSp1, iSp1) = 0.0_dp
    end do

    self%tCoordsUpdated = .false.

  end subroutine initialize


  !> Update internal stored coordinates
  subroutine updateCoords(self, neighList, img2CentCell, coords, species0)

    !> data structure
    class(TChargeModel5), intent(inout) :: self

    !> list of neighbours to atoms
    type(TNeighbourList), intent(in) :: neighList

    !> image to central cell atom index
    integer, intent(in) :: img2CentCell(:)

    !> atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    integer, allocatable :: nNeigh(:)

    allocate(nNeigh(self%nAtom))
    call getNrOfNeighboursForAll(nNeigh, neighList, self%rCutoff)
    call getCorrection(self, nNeigh, neighList%iNeighbour, img2CentCell, neighList%neighDist2, species0, coords)

    self%tCoordsUpdated = .true.

  end subroutine updateCoords


  !> update internal copy of lattice vectors
  subroutine updateLatVecs(self, latVecs)

    !> data structure
    class(TChargeModel5), intent(inout) :: self

    !> lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    @:ASSERT(self%tPeriodic)
    @:ASSERT(all(shape(latvecs) == shape(self%latvecs)))

    self%latVecs(:,:) = latVecs

    self%tCoordsUpdated = .false.

  end subroutine updateLatVecs


  !> get charge corrections
  subroutine addCharges(self, charges)

    !> data structure
    class(TChargeModel5), intent(inout) :: self

    !> energy contributions for each atom
    real(dp), intent(inout) :: charges(:)

    @:ASSERT(self%tCoordsUpdated)
    @:ASSERT(size(charges) == self%nAtom)

    charges(:) = charges + self%cm5

  end subroutine addCharges


  !> get force contributions
  subroutine addGradients(self, dEdcm5, gradients)

    !> data structure
    class(TChargeModel5), intent(inout) :: self

    !> Derivative w.r.t. CM5 correction
    real(dp), intent(in) :: dEdcm5(:)

    !> gradient contributions for each atom
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(self%tCoordsUpdated)
    @:ASSERT(all(shape(stress) == [3, self%nAtom]))

    call gemv(gradients, self%dcm5dr, dEdcm5, beta=1.0_dp)

  end subroutine addGradients


  !> get stress tensor contributions
  subroutine addStress(self, dEdcm5, stress)

    !> data structure
    class(TChargeModel5), intent(inout) :: self

    !> Derivative w.r.t. CM5 correction
    real(dp), intent(in) :: dEdcm5(:)

    !> Stress tensor contributions (not volume scaled)
    real(dp), intent(inout) :: stress(:,:)

    @:ASSERT(self%tCoordsUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))

    call gemv(stress, self%dcm5dL, dEdcm5, beta=1.0_dp)

  end subroutine addStress


  !> Distance cut off for charge model
  function getRCutoff(self) result(cutoff)

    !> data structure
    class(TChargeModel5), intent(inout) :: self

    !> resulting cutoff
    real(dp) :: cutoff

    cutoff = self%rCutoff

  end function getRCutoff


  !> Calculate CM5 correction for this geometry
  subroutine getCorrection(self, nNeighbour, iNeighbour, img2CentCell, neighDist2, species, coords)

    !> data structure
    class(TChargeModel5), intent(inout) :: self

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

    self%cm5(:) = 0.0_dp
    self%dcm5dr(:, :, :) = 0.0_dp
    self%dcm5dL(:, :, :) = 0.0_dp

    do iAt1 = 1, self%nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeighbour(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        if (iSp1 == iSp2) cycle  ! includes iAt1 == iAt2f case
        dist = sqrt(neighDist2(iNeigh, iAt1))
        vec(:) = coords(:, iAt1) - coords(:, iAt2)
        p12 = self%pairParam(iSp1, iSp2)
        p21 = self%pairParam(iSp2, iSp1)

        dEr = exp(-self%alpha*(dist-self%atomicRad(iSp1)-self%atomicRad(iSp2)))
        dGr = dEr * self%alpha * vec/dist
        dSr = spread(dGr, 1, 3) * spread(vec, 2, 3)

        self%cm5(iAt1) = self%cm5(iAt1) + dEr * p12
        self%cm5(iAt2f) = self%cm5(iAt2f) + dEr * p21

        self%dcm5dr(:, iAt1, iAt1) = self%dcm5dr(:, iAt1, iAt1) - dGr * p12
        self%dcm5dr(:, iAt2f, iAt2f) = self%dcm5dr(:, iAt2f, iAt2f) + dGr * p21
        self%dcm5dr(:, iAt1, iAt2f) = self%dcm5dr(:, iAt1, iAt2f) - dGr * p21
        self%dcm5dr(:, iAt2f, iAt1) = self%dcm5dr(:, iAt2f, iAt1) + dGr * p12

        self%dcm5dL(:, :, iAt1) = self%dcm5dL(:, :, iAt1) + dSr * p12
        self%dcm5dL(:, :, iAt2f) = self%dcm5dL(:, :, iAt2f) + dSr * p21

      end do
    end do

  end subroutine getCorrection


  !> Get atomic radius for species with a given symbol
  elemental function getAtomicRadiusSymbol(symbol) result(radius)

    !> Element symbol
    character(len=*), intent(in) :: symbol

    !> atomic radius
    real(dp) :: radius

    radius = getAtomicRadius(symbolToNumber(symbol))

  end function getAtomicRadiusSymbol


  !> Get atomic radius for species with a given atomic number
  elemental function getAtomicRadiusNumber(number) result(radius)

    !> Atomic number
    integer, intent(in) :: number

    !> atomic radius
    real(dp) :: radius

    if (number > 0 .and. number <= size(atomicRadii, dim=1)) then
      radius = AtomicRadii(number)
    else
      radius = -1.0_dp
    end if

  end function getAtomicRadiusNumber


  !> Get pair parameter for species with a given symbols
  elemental function getPairParameterSymbol(symbol1, symbol2) result(pairPar)

    !> Element symbol
    character(len=*), intent(in) :: symbol1

    !> Element symbol
    character(len=*), intent(in) :: symbol2

    !> Pair parameter
    real(dp) :: pairPar

    pairPar = getPairParameter(symbolToNumber(symbol1), symbolToNumber(symbol2))

  end function getPairParameterSymbol


  !> Get pair parameter for species with a given atomic numbers
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
