!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Coordination number implementation
module dftbp_coordnumber
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_blasroutines, only : gemv
  use dftbp_constants, only : pi, AA__Bohr, symbolToNumber
  use dftbp_message, only : error
  use dftbp_periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use dftbp_simplealgebra, only : determinant33
  implicit none
  private

  public :: TCNCont, TCNInput, cnType, init
  public :: getElectronegativity, getCovalentRadius


  !> Get electronegativity for a species
  interface getElectronegativity
    module procedure :: getElectronegativitySymbol
    module procedure :: getElectronegativityNumber
  end interface getElectronegativity


  !> Get atomic radius for a species
  interface getCovalentRadius
    module procedure :: getCovalentRadiusSymbol
    module procedure :: getCovalentRadiusNumber
  end interface getCovalentRadius


  !> Covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
  !> 188-197), values for metals decreased by 10%.
  real(dp), parameter :: CovalentRadii(1:118) = [ &
    & 0.32_dp,0.46_dp, & ! H,He
    & 1.20_dp,0.94_dp,0.77_dp,0.75_dp,0.71_dp,0.63_dp,0.64_dp,0.67_dp, & ! Li-Ne
    & 1.40_dp,1.25_dp,1.13_dp,1.04_dp,1.10_dp,1.02_dp,0.99_dp,0.96_dp, & ! Na-Ar
    & 1.76_dp,1.54_dp, & ! K,Ca
    &                 1.33_dp,1.22_dp,1.21_dp,1.10_dp,1.07_dp, & ! Sc-
    &                 1.04_dp,1.00_dp,0.99_dp,1.01_dp,1.09_dp, & ! -Zn
    &                 1.12_dp,1.09_dp,1.15_dp,1.10_dp,1.14_dp,1.17_dp, & ! Ga-Kr
    & 1.89_dp,1.67_dp, & ! Rb,Sr
    &                 1.47_dp,1.39_dp,1.32_dp,1.24_dp,1.15_dp, & ! Y-
    &                 1.13_dp,1.13_dp,1.08_dp,1.15_dp,1.23_dp, & ! -Cd
    &                 1.28_dp,1.26_dp,1.26_dp,1.23_dp,1.32_dp,1.31_dp, & ! In-Xe
    & 2.09_dp,1.76_dp, & ! Cs,Ba
    &         1.62_dp,1.47_dp,1.58_dp,1.57_dp,1.56_dp,1.55_dp,1.51_dp, & ! La-Eu
    &         1.52_dp,1.51_dp,1.50_dp,1.49_dp,1.49_dp,1.48_dp,1.53_dp, & ! Gd-Yb
    &                 1.46_dp,1.37_dp,1.31_dp,1.23_dp,1.18_dp, & ! Lu-
    &                 1.16_dp,1.11_dp,1.12_dp,1.13_dp,1.32_dp, & ! -Hg
    &                 1.30_dp,1.30_dp,1.36_dp,1.31_dp,1.38_dp,1.42_dp, & ! Tl-Rn
    & 2.01_dp,1.81_dp, & ! Fr,Ra
    &         1.67_dp,1.58_dp,1.52_dp,1.53_dp,1.54_dp,1.55_dp,1.49_dp, & ! Ac-Am
    &         1.49_dp,1.51_dp,1.51_dp,1.48_dp,1.50_dp,1.56_dp,1.58_dp, & ! Cm-No
    &                 1.45_dp,1.41_dp,1.34_dp,1.29_dp,1.27_dp, & ! Lr-
    &                 1.21_dp,1.16_dp,1.15_dp,1.09_dp,1.22_dp, & ! -Cn
    &                 1.36_dp,1.43_dp,1.46_dp,1.58_dp,1.48_dp,1.57_dp] & ! Nh-Og
    & * AA__Bohr * 4.0_dp / 3.0_dp


  !> Pauling electronegativities, used for the covalent coordination number.
  real(dp), parameter :: paulingEN(1:118) = [ &
    & 2.20_dp,3.00_dp, & ! H,He
    & 0.98_dp,1.57_dp,2.04_dp,2.55_dp,3.04_dp,3.44_dp,3.98_dp,4.50_dp, & ! Li-Ne
    & 0.93_dp,1.31_dp,1.61_dp,1.90_dp,2.19_dp,2.58_dp,3.16_dp,3.50_dp, & ! Na-Ar
    & 0.82_dp,1.00_dp, & ! K,Ca
    &                 1.36_dp,1.54_dp,1.63_dp,1.66_dp,1.55_dp, & ! Sc-
    &                 1.83_dp,1.88_dp,1.91_dp,1.90_dp,1.65_dp, & ! -Zn
    &                 1.81_dp,2.01_dp,2.18_dp,2.55_dp,2.96_dp,3.00_dp, & ! Ga-Kr
    & 0.82_dp,0.95_dp, & ! Rb,Sr
    &                 1.22_dp,1.33_dp,1.60_dp,2.16_dp,1.90_dp, & ! Y-
    &                 2.20_dp,2.28_dp,2.20_dp,1.93_dp,1.69_dp, & ! -Cd
    &                 1.78_dp,1.96_dp,2.05_dp,2.10_dp,2.66_dp,2.60_dp, & ! In-Xe
    & 0.79_dp,0.89_dp, & ! Cs,Ba
    &         1.10_dp,1.12_dp,1.13_dp,1.14_dp,1.15_dp,1.17_dp,1.18_dp, & ! La-Eu
    &         1.20_dp,1.21_dp,1.22_dp,1.23_dp,1.24_dp,1.25_dp,1.26_dp, & ! Gd-Yb
    &                 1.27_dp,1.30_dp,1.50_dp,2.36_dp,1.90_dp, & ! Lu-
    &                 2.20_dp,2.20_dp,2.28_dp,2.54_dp,2.00_dp, & ! -Hg
    &                 1.62_dp,2.33_dp,2.02_dp,2.00_dp,2.20_dp,2.20_dp, & ! Tl-Rn
    ! only dummies below
    & 1.50_dp,1.50_dp, & ! Fr,Ra
    &         1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp, & ! Ac-Am
    &         1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp, & ! Cm-No
    &                 1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp, & ! Rf-
    &                 1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp, & ! Rf-Cn
    &                 1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp ] ! Nh-Og


  !> Possible counting functions for calculating coordination numbers
  type :: TCNCountEnum

    !> Counting function not specified
    integer :: invalid = 0

    !> Original DFT-D3 coordination number
    integer :: exp = 1

    !> Faster decaying error function CN, better for dense systems
    integer :: erf = 2

    !> Error function CN with covalency correction
    integer :: cov = 3

    !> Particular long-ranged version of the DFT-D3 coordination number
    integer :: gfn = 4

  end type TCNCountEnum

  !> Enumerator for different coordination number types
  type(TCNCountEnum), parameter :: cnType = TCNCountEnum()


  !> Input to generate coordination number container
  type :: TCNInput

    !> Coordination number type
    integer :: cnType

    !> Real space cutoff
    real(dp) :: rCutoff

    !> Upper bound for coordination number
    real(dp) :: maxCN

    !> Covalent radii
    real(dp), allocatable :: covRad(:)

    !> Electronegativity
    real(dp), allocatable :: en(:)

  end type TCNInput


  !> Coordination number container
  type :: TCNCont
    private

    !> number of atoms
    integer :: nAtom = 0

    !> lattice vectors if periodic
    real(dp) :: latVecs(3, 3) = 0.0_dp

    !> is this periodic
    logical :: tPeriodic

    !> are the coordinates current?
    logical :: tCoordsUpdated = .false.

    !> Use covalency correction from EN
    logical :: tENScale

    !> Real space cutoff
    real(dp) :: rCutoff = 0.0_dp

    !> Steepness of the counting function
    real(dp) :: kcn

    !> Cut coordination number
    logical :: tCutCN

    !> Upper bound for coordination number
    real(dp) :: maxCN

    !> Covalent radii
    real(dp), allocatable :: covRad(:)

    !> Electronegativity
    real(dp), allocatable :: en(:)

    !> Coordination number
    real(dp), public, allocatable :: cn(:)

    !> Derivative of coordination numbers w.r.t. coordinates
    real(dp), public, allocatable :: dcndr(:, :, :)

    !> Derivative of coordination numbers w.r.t. strain deformations
    real(dp), public, allocatable :: dcndL(:, :, :)

    !> Counting function for CN
    procedure(countFunction), nopass, pointer :: countFunc => null()

    !> Derivative of counting function w.r.t. distance
    procedure(countFunction), nopass, pointer :: countDeriv => null()

  contains

    !> update internal copy of coordinates
    procedure :: updateCoords

    !> update internal copy of lattice vectors
    procedure :: updateLatVecs

    !> get real space cutoff
    procedure :: getRCutoff

    !> get force contributions
    procedure :: addGradients

    !> get stress tensor contributions
    procedure :: addStress

  end type TCNCont


  abstract interface
    !> Abstract interface for the counting function (and its derivative)
    pure function countFunction(k, r, r0)
      import :: dp

      !> Constant for counting function
      real(dp), intent(in) :: k

      !> Actual distance
      real(dp), intent(in) :: r

      !> Critical distance
      real(dp), intent(in) :: r0

      !> Value of the counting function in the range of [0,1]
      real(dp) :: countFunction
    end function countFunction
  end interface

  !> Initialize container from geometry
  interface init
    module procedure :: initialize
  end interface init


contains


  !> Initialize coordination number container
  subroutine initialize(this, input, nAtom, latVecs)

    !> Initialised instance at return
    type(TCNCont), intent(out) :: this

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Input for container
    type(TCNInput), intent(in) :: input

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    @:ASSERT(allocated(input%covRad))
    @:ASSERT(allocated(input%en))

    this%tPeriodic = present(latVecs)
    if (this%tPeriodic) then
      call this%updateLatVecs(LatVecs)
    end if
    this%nAtom = nAtom

    allocate(this%cn(nAtom))
    allocate(this%dcndr(3, nAtom, nAtom))
    allocate(this%dcndL(3, 3, nAtom))

    this%rCutoff = input%rCutoff

    select case(input%cnType)
    case default
      call error("Fatal programming error in dftbp_coordnumber!")
    case(cnType%erf)
      this%countFunc => erfCount
      this%countDeriv => derfCount
      this%kcn = 7.5_dp
      this%tENScale = .false.
    case(cnType%exp)
      this%countFunc => expCount
      this%countDeriv => dexpCount
      this%kcn = 16.0_dp
      this%tENScale = .false.
    case(cnType%cov)
      this%countFunc => erfCount
      this%countDeriv => derfCount
      this%kcn = 7.5_dp
      this%tENScale = .true.
    case(cnType%gfn)
      this%countFunc => gfnCount
      this%countDeriv => dgfnCount
      this%kcn = 10.0_dp
      this%tENScale = .false.
    end select

    this%tCutCN = input%maxCN > 0.0_dp
    this%maxCN = input%maxCN

    allocate(this%covRad(size(input%covRad)))
    allocate(this%en(size(input%en)))
    this%covRad(:) = input%covRad
    this%en(:) = input%en

    this%tCoordsUpdated = .false.

  end subroutine initialize


  !> Update internal stored coordinates
  subroutine updateCoords(this, neighList, img2CentCell, coords, species0)

    !> data structure
    class(TCNCont), intent(inout) :: this

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
    call getCoordinationNumber(this%nAtom, coords, species0, nNeigh, &
        & neighList%iNeighbour, neighList%neighDist2, img2CentCell, &
        & this%covRad, this%en, this%tENScale, this%kcn, this%countFunc, &
        & this%countDeriv, this%cn, this%dcndr, this%dcndL)

    if (this%tCutCN) then
      call cutCoordinationNumber(this%nAtom, this%cn, this%dcndr, this%dcndL, &
          & this%maxCN)
    end if

    this%tCoordsUpdated = .true.

  end subroutine updateCoords


  !> update internal copy of lattice vectors
  subroutine updateLatVecs(this, latVecs)

    !> data structure
    class(TCNCont), intent(inout) :: this

    !> lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    @:ASSERT(this%tPeriodic)
    @:ASSERT(all(shape(latvecs) == shape(this%latvecs)))

    this%latVecs(:,:) = latVecs

    this%tCoordsUpdated = .false.

  end subroutine updateLatVecs


  !> get force contributions
  subroutine addGradients(this, dEdcn, gradients)

    !> data structure
    class(TCNCont), intent(inout) :: this

    !> Derivative w.r.t. CM5 correction
    real(dp), intent(in) :: dEdcn(:)

    !> gradient contributions for each atom
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(this%tCoordsUpdated)

    call gemv(gradients, this%dcndr, dEdcn, beta=1.0_dp)

  end subroutine addGradients


  !> get stress tensor contributions
  subroutine addStress(this, dEdcn, stress)

    !> data structure
    class(TCNCont), intent(inout) :: this

    !> Derivative w.r.t. CM5 correction
    real(dp), intent(in) :: dEdcn(:)

    !> Stress tensor contributions (not volume scaled)
    real(dp), intent(inout) :: stress(:,:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))

    call gemv(stress, this%dcndL, dEdcn, beta=1.0_dp)

  end subroutine addStress


  !> Distance cut off for coordination number
  function getRCutoff(this) result(cutoff)

    !> data structure
    class(TCNCont), intent(inout) :: this

    !> resulting cutoff
    real(dp) :: cutoff

    cutoff = this%rCutoff

  end function getRCutoff


  !> Actual implementation of the coordination number calculation
  subroutine getCoordinationNumber(nAtom, coords, species, nNeighbour, iNeighbour,&
      & neighDist2, img2CentCell, covalentRadius, electronegativity, tENscale, &
      & kcn, countFunc, countDeriv, cn, dcndr, dcndL)

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Coordinates of the atoms (including images)
    real(dp), intent(in) :: coords(:, :)

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist.
    integer, intent(in) :: iNeighbour(0:, :)

    !> Square distances of the neighbours.
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Mapping into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Use covalent coordination number by weighting with EN difference.
    logical, intent(in) :: tENscale

    !> Covalent Radii for counting function
    real(dp), intent(in) :: covalentRadius(:)

    !> Electronegativities for covalancy scaling
    real(dp), intent(in) :: electronegativity(:)

    !> Steepness of the counting function for the fractional coordination number.
    real(dp), intent(in) :: kcn

    !> Counting function for CN
    procedure(countFunction) :: countFunc

    !> Derivative of counting function w.r.t. distance
    procedure(countFunction) :: countDeriv

    !> Error function coordination number.
    real(dp), intent(out) :: cn(:)

    !> Derivative of the CN with respect to the Cartesian coordinates.
    real(dp), intent(out) :: dcndr(:, :, :)

    !> Derivative of the CN with respect to strain deformations.
    real(dp), intent(out) :: dcndL(:, :, :)

    !> EN scaling parameter
    real(dp), parameter :: k4 = 4.10451_dp
    real(dp), parameter :: k5 = 19.08857_dp
    real(dp), parameter :: k6 = 2*11.28174_dp**2

    integer :: iAt1, iSp1, iAt2, iAt2f, iSp2, iNeigh
    real(dp) :: r2, r1, rc, vec(3), countf, countd(3), stress(3, 3), dEN

    cn(:) = 0.0_dp
    dcndr(:, :, :) = 0.0_dp
    dcndL(:, :, :) = 0.0_dp

    !$omp parallel do default(none) schedule(runtime) &
    !$omp reduction(+:cn, dcndr, dcndL) &
    !$omp shared(nAtom, species, nNeighbour, iNeighbour, coords, img2CentCell) &
    !$omp shared(neighDist2, covalentRadius, tENScale, electronegativity, kcn) &
    !$omp private(iAt1, iSp1, iAt2, vec, iAt2f, iSp2, r2, r1, rc, dEN, countf) &
    !$omp private(countd, stress)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeighbour(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        vec(:) = coords(:, iAt1) - coords(:, iAt2)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        r2 = neighDist2(iNeigh, iAt1)
        r1 = sqrt(r2)

        rc = covalentRadius(iSp1) + covalentRadius(iSp2)

        if (tENScale) then
          dEN = abs(electronegativity(iSp1) - electronegativity(iSp2))
          dEN = k4 * exp(-(dEN + k5)**2 / k6)
        else
          dEN = 1.0_dp
        end if

        countf = dEN * countFunc(kcn, r1, rc)
        countd = dEN * countDeriv(kcn, r1, rc) * vec / r1

        cn(iAt1) = cn(iAt1) + countf
        if (iAt1 /= iAt2f) then
          cn(iAt2f) = cn(iAt2f) + countf
        end if

        dcndr(:, iAt1, iAt1) = dcndr(:, iAt1, iAt1) + countd
        dcndr(:, iAt2f, iAt2f) = dcndr(:, iAt2f, iAt2f) - countd
        dcndr(:, iAt1, iAt2f) = dcndr(:, iAt1, iAt2f) + countd
        dcndr(:, iAt2f, iAt1) = dcndr(:, iAt2f, iAt1) - countd

        stress = spread(countd, 1, 3) * spread(vec, 2, 3)

        dcndL(:, :, iAt1) = dcndL(:, :, iAt1) - stress
        if (iAt1 /= iAt2f) then
          dcndL(:, :, iAt2f) = dcndL(:, :, iAt2f) - stress
        end if

      end do
    end do
    !$omp end parallel do
  end subroutine getCoordinationNumber


  !> Error function counting function for coordination number contributions.
  pure function erfCount(k, r, r0) result(count)

    !> Steepness of the counting function.
    real(dp), intent(in) :: k

    !> Current distance.
    real(dp), intent(in) :: r

    !> Cutoff radius.
    real(dp), intent(in) :: r0

    real(dp) :: count

    count = 0.5_dp * (1.0_dp + erf(-k*(r-r0)/r0))

  end function erfCount


  !> Derivative of the counting function w.r.t. the distance.
  pure function derfCount(k, r, r0) result(count)

    !> Steepness of the counting function.
    real(dp), intent(in) :: k

    !> Current distance.
    real(dp), intent(in) :: r

    !> Cutoff radius.
    real(dp), intent(in) :: r0

    real(dp), parameter :: sqrtpi = sqrt(pi)

    real(dp) :: count

    count = -k/sqrtpi/r0*exp(-k**2*(r-r0)**2/r0**2)

  end function derfCount


  !> Exponential counting function for coordination number contributions.
  pure function expCount(k, r, r0) result(count)

    !> Steepness of the counting function.
    real(dp), intent(in) :: k

    !> Current distance.
    real(dp), intent(in) :: r

    !> Cutoff radius.
    real(dp), intent(in) :: r0

    real(dp) :: count

    count =1.0_dp/(1.0_dp+exp(-k*(r0/r-1.0_dp)))

  end function expCount


  !> Derivative of the counting function w.r.t. the distance.
  pure function dexpCount(k, r, r0) result(count)

    !> Steepness of the counting function.
    real(dp), intent(in) :: k

    !> Current distance.
    real(dp), intent(in) :: r

    !> Cutoff radius.
    real(dp), intent(in) :: r0

    real(dp) :: count
    real(dp) :: expterm

    expterm = exp(-k*(r0/r-1._dp))

    count = (-k*r0*expterm)/(r**2*((expterm+1._dp)**2))

  end function dexpCount


  !> Exponential counting function for coordination number contributions.
  pure function gfnCount(k, r, r0) result(count)

    !> Steepness of the counting function.
    real(dp), intent(in) :: k

    !> Current distance.
    real(dp), intent(in) :: r

    !> Cutoff radius.
    real(dp), intent(in) :: r0

    real(dp) :: count

    count = expCount(k, r, r0) * expCount(2*k, r, r0+2)

  end function gfnCount


  !> Derivative of the counting function w.r.t. the distance.
  pure function dgfnCount(k, r, r0) result(count)

    !> Steepness of the counting function.
    real(dp), intent(in) :: k

    !> Current distance.
    real(dp), intent(in) :: r

    !> Cutoff radius.
    real(dp), intent(in) :: r0

    real(dp) :: count

    count = dexpCount(k, r, r0) * expCount(2*k, r, r0+2) &
        &  + expCount(k, r, r0) * dexpCount(2*k, r, r0+2)

  end function dgfnCount


  !> Cutoff function for large coordination numbers
  pure subroutine cutCoordinationNumber(nAtom, cn, dcndr, dcndL, maxCN)

   !> number of atoms
    integer, intent(in) :: nAtom

    !> on input coordination number, on output modified CN
    real(dp), intent(inout) :: cn(:)

    !> on input derivative of CN w.r.t. cartesian coordinates,
    !> on output derivative of modified CN
    real(dp), intent(inout), optional :: dcndr(:, :, :)

    !> on input derivative of CN w.r.t. strain deformation,
    !> on output derivative of modified CN
    real(dp), intent(inout), optional :: dcndL(:, :, :)

    !> maximum CN (not strictly obeyed)
    real(dp), intent(in), optional :: maxCN

    real(dp) :: cnmax
    integer :: iAt

    if (present(maxCN)) then
      cnmax = maxCN
    else
      cnmax = 4.5_dp
    end if

    if (cnmax <= 0.0_dp) return

    if (present(dcndL)) then
      do iAt = 1, nAtom
        dcndL(:, :, iAt) = dcndL(:, :, iAt) * dCutCN(cn(iAt), cnmax)
      end do
    end if

    if (present(dcndr)) then
      do iAt = 1, nAtom
        dcndr(:, :, iAt) = dcndr(:, :, iAt) * dCutCN(cn(iAt), cnmax)
      end do
    end if

    do iAt = 1, nAtom
      cn(iAt) = cutCN(cn(iAt), cnmax)
    end do

  end subroutine cutCoordinationNumber


  !> Cutting function for the coordination number.
  elemental function cutCN(cn, cut) result(cnp)

    !> Current coordination number.
    real(dp), intent(in) :: cn

    !> Cutoff for the CN, this is not the maximum value.
    real(dp), intent(in) :: cut

    !> Cuting function vlaue
    real(dp) :: cnp

    cnp = log(1.0_dp + exp(cut)) - log(1.0_dp + exp(cut - cn))

  end function cutCN


  !> Derivative of the cutting function w.r.t. coordination number
  elemental function dCutCN(cn, cut) result(dcnpdcn)

    !> Current coordination number.
    real(dp), intent(in) :: cn

    !> Cutoff for the CN, this is not the maximum value.
    real(dp), intent(in) :: cut

    !> Derivative of the cutting function
    real(dp) :: dcnpdcn

    dcnpdcn = exp(cut)/(exp(cut) + exp(cn))

  end function dCutCn


  !> Populate covalent radius field from speciesNames
  subroutine covRadFromSpecies(this, speciesNames)

    !> Instance of the CN input data
    class(TCNInput), intent(inout) :: this

    !> Element symbols for all species
    character(len=*), intent(in) :: speciesNames(:)

    integer :: iSp

    @:ASSERT(.not.allocated(this%covRad))

    allocate(this%covRad(size(speciesNames)))
    do iSp = 1, size(speciesNames)
      this%covRad(iSp) = getCovalentRadius(speciesNames(iSp))
    end do

  end subroutine covRadFromSpecies


  !> Get covalent radius for species with a given symbol
  elemental function getCovalentRadiusSymbol(symbol) result(radius)

    !> Element symbol
    character(len=*), intent(in) :: symbol

    !> atomic radius
    real(dp) :: radius

    radius = getCovalentRadius(symbolToNumber(symbol))

  end function getCovalentRadiusSymbol


  !> Get covalent radius for species with a given atomic number
  elemental function getCovalentRadiusNumber(number) result(radius)

    !> Atomic number
    integer, intent(in) :: number

    !> atomic radius
    real(dp) :: radius

    if (number > 0 .and. number <= size(CovalentRadii, dim=1)) then
      radius = CovalentRadii(number)
    else
      radius = -1.0_dp
    end if

  end function getCovalentRadiusNumber


  !> Populate electronegativity field from speciesNames
  subroutine enFromSpecies(this, speciesNames)

    !> Instance of the CN input data
    class(TCNInput), intent(inout) :: this

    !> Element symbols for all species
    character(len=*), intent(in) :: speciesNames(:)

    integer :: iSp

    @:ASSERT(.not.allocated(this%en))

    allocate(this%en(size(speciesNames)))
    do iSp = 1, size(speciesNames)
      this%en(iSp) = getElectronegativity(speciesNames(iSp))
    end do

  end subroutine enFromSpecies


  !> Get electronegativity for species with a given symbol
  elemental function getElectronegativitySymbol(symbol) result(en)

    !> Element symbol
    character(len=*), intent(in) :: symbol

    !> atomic EN
    real(dp) :: en

    en = getElectronegativity(symbolToNumber(symbol))

  end function getElectronegativitySymbol


  !> Get electronegativity for species with a given atomic number
  elemental function getElectronegativityNumber(number) result(en)

    !> Atomic number
    integer, intent(in) :: number

    !> atomic EN
    real(dp) :: en

    if (number > 0 .and. number <= size(paulingEN, dim=1)) then
      en = paulingEN(number)
    else
      en = -1.0_dp
    end if

  end function getElectronegativityNumber


end module dftbp_coordnumber
