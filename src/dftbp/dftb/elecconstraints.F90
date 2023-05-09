!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module to impose constraints on the electronic ground state.
module dftbp_dftb_elecconstraints
  use dftbp_common_accuracy, only : dp
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_typegeometry, only : TGeometry
  use dftbp_extlibs_xmlf90, only : fnode, string, char, getLength, getItem1, fnodeList
  use dftbp_type_wrappedintr, only : TWrappedInt1, TWrappedReal1, TWrappedReal2
  use dftbp_geoopt_package, only : TOptimizer, TOptimizerInput, createOptimizer
  use dftbp_io_hsdutils, only : getChildValue, getChildren, getSelectedAtomIndices
  use dftbp_dftbplus_input_geoopt, only : readOptimizerInput
  implicit none

  private
  public :: TElecConstraint, TElecConstraint_init, TElecConstraintInput
  public :: readElecConstraintInput


  type TElecConstraintInput

    !> Optimiser input choice
    class(TOptimizerInput), allocatable :: optimiser

    !> Group of atoms in a constraint
    type(TWrappedInt1), allocatable :: atomGrp(:)

    !> Constraint targets for atom groups
    real(dp), allocatable :: atomNc(:)

    !> Direction of constraint in (q,m) space
    type(TWrappedReal1), allocatable :: atomSpinDir(:)

    !> Derivative tolerance for constraint
    real(dp) :: constrTol

    !> Number of iterations for enforcing constraint
    integer :: nConstrIter

    !> True, if converged micro-iterations are required
    logical :: isConstrConvRequired

  end type TElecConstraintInput


  type TElecConstraint

    !> Value of the constraint
    real(dp), allocatable :: Nc(:)

    !> Potential
    real(dp), allocatable :: Vc(:)

    ! Weighting function for constrain

    !> Atom(s) involved in each constrain
    type(TWrappedInt1), allocatable :: wAt(:)

    !> Atomic orbital(s) involved in each constrain
    type(TWrappedInt1), allocatable :: wAtOrb(:)

    !> Atomic orbital charge/spin quaternion involved in each constrain
    type(TWrappedReal2), allocatable :: wAtSpin(:)

    !> General optimiser
    class(TOptimizer), allocatable :: potOpt

    !> Derivative tolerance for constraint
    real(dp) :: constrTol

    !> Number of iterations for enforcing constraint
    integer :: nConstrIter

    !> True, if converged micro-iterations are required
    logical :: isConstrConvRequired

  contains

    procedure getConstrainShift
    procedure propagateConstraints
    procedure getMaxIter

  end type TElecConstraint


contains


  !> General entry point to read constraint on the electronic ground state.
  subroutine readElecConstraintInput(node, geo, input, tSpinPol)

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Control structure to be filled
    type(TElecConstraintInput), intent(out) :: input

    !> True, if this is a spin polarized calculation
    logical, intent(in) :: tSpinPol

    type(fnode), pointer :: val, child1, child2, child3
    type(fnodeList), pointer :: children
    type(string) :: buffer
    integer :: iConstr, nConstr

    call getChildValue(node, "Optimiser", child1, "FIRE")
    call readOptimizerInput(child1, input%optimiser)

    call getChildValue(node, "ConstrTolerance", input%constrTol, 1.0e-08_dp)
    call getChildValue(node, "MaxConstrIterations", input%nConstrIter, 100)
    call getChildValue(node, "ConvergentConstrOnly", input%isConstrConvRequired, .true.)

    call getChildValue(node, "Regions", val, "", child=child1, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)

    ! Read specification for regions of atoms
    call getChildren(child1, "Atoms", children)
    nConstr = getLength(children)

    allocate(input%atomGrp(nConstr))
    allocate(input%atomNc(nConstr))
    allocate(input%atomSpinDir(nConstr))

    do iConstr = 1, nConstr
      call getItem1(children, iConstr, child2)
      call getChildValue(child2, "Domain", buffer, child=child3, multiple=.true.)
      call getSelectedAtomIndices(child3, char(buffer), geo%speciesNames, geo%species,&
          & input%atomGrp(iConstr)%data)
      call getChildValue(child2, "Population", input%atomNc(iConstr))
      ! Functionality currently restricted to charges
      if (tSpinPol) then
        allocate(input%atomSpinDir(iConstr)%data(2))
        input%atomSpinDir(iConstr)%data(1) = 1.0_dp
      else
        allocate(input%atomSpinDir(iConstr)%data(1))
        input%atomSpinDir(iConstr)%data(1) = 1.0_dp
      end if
    end do

  end subroutine readElecConstraintInput


  !> Initialises the constraints structure.
  subroutine TElecConstraint_init(this, input, orb)

    !> Constrain structure instance
    type(TElecConstraint), intent(out) :: this

    !> Input data structure
    type(TElecConstraintInput), intent(inout) :: input

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    integer :: iConstr, nConstr, ii, jj, iAt, nOrb, nSpin

    nConstr = size(input%atomGrp)

    allocate(this%Vc(nConstr))
    ! should enable optional initialization of Vc from input
    this%Vc(:) = 0.0_dp

    call createOptimizer(input%optimiser, nConstr, this%potOpt)

    this%nConstrIter = input%nConstrIter
    this%isConstrConvRequired = input%isConstrConvRequired
    this%constrTol = input%constrTol
    this%Nc = input%AtomNc

    allocate(this%wAt(nConstr))
    allocate(this%wAtOrb(nConstr))
    allocate(this%wAtSpin(nConstr))

    ! Allocate + initialize arrays and build index mappings
    do iConstr = 1, nConstr
      ! Count orbitals subject to constraints
      nOrb = 0
      do ii = 1, size(input%atomGrp(iConstr)%data)
        iAt = input%atomGrp(iConstr)%data(ii)
        nOrb = nOrb + orb%nOrbAtom(iAt)
      end do
      allocate(this%wAt(iConstr)%data(nOrb))
      allocate(this%wAtOrb(iConstr)%data(nOrb))
      nSpin = size(input%atomSpinDir(iConstr)%data)
      allocate(this%wAtSpin(iConstr)%data(nOrb,nSpin))
      this%wAt(iConstr)%data(:) = 0
      this%wAtOrb(iConstr)%data(:) = 0
      nOrb = 0
      do ii = 1, size(input%atomGrp(iConstr)%data)
        iAt = input%atomGrp(iConstr)%data(ii)
        do jj = 1, orb%nOrbAtom(iAt)
          this%wAt(iConstr)%data(nOrb+jj) = iAt
          this%wAtOrb(iConstr)%data(nOrb+jj) = jj
          this%wAtSpin(iConstr)%data(nOrb+jj,:) = input%atomSpinDir(iConstr)%data
        end do
        nOrb = nOrb + orb%nOrbAtom(iAt)
      end do
    end do

  end subroutine TElecConstraint_init


  !> Returns maximum number of iterations for constraint driver.
  pure function getMaxIter(this) result(maxIter)

    !> Class instance
    class(TElecConstraint), intent(in) :: this

    !> Obtained maximum number of iterations
    integer :: maxIter

    maxIter = this%nConstrIter

  end function getMaxIter


  !> Applies electronic constraints to system.
  subroutine propagateConstraints(this, qq, energy, deltaW, dWdVcMax, tConverged)

    !> Class instance
    class(TElecConstraint), intent(inout) :: this

    !> Mulliken populations
    real(dp), intent(in) :: qq(:,:,:)

    !> Energy
    real(dp), intent(in) :: energy

    !> Contribution to free energy functional from constraint(s)
    real(dp), intent(out) :: deltaW

    !> Maximum derivative of energy functional with respect to Vc
    real(dp), intent(out) :: dWdVcMax

    !> Gradient convergence achieved
    logical, intent(out) :: tConverged

    !! Derivative of energy functional with respect to Vc
    real(dp), allocatable :: dWdVc(:)

    !! Potential displacement proposed by optimizer
    real(dp) :: potDisplace(size(this%Vc))

    !! Iterates over constraints
    integer :: iConstr

    !! Number of constraints requested by the user
    integer :: nConstr

    nConstr = size(this%wAt)
    allocate(dWdVc(nConstr))
    dWdVc(:) = 0.0_dp
    deltaW = 0.0_dp

    do iConstr = 1, nConstr
      call getConstrainEnergyAndPotQ(deltaW, dWdVc(iConstr), this%Vc(iConstr), this%Nc(iConstr),&
          & this%wAt(iConstr)%data, this%wAtOrb(iConstr)%data, this%wAtSpin(iConstr)%data,&
          & qq)
    end do

    ! Optimizers set up to minimize, therefore sign change in total energy and gradients
    ! call next(this%potOpt, -(energy + deltaW), -dWdVc, this%Vc, tConverged)
    call this%potOpt%step(energy + deltaW, -dWdVc, potDisplace)
    this%Vc(:) = this%Vc + potDisplace

    dWdVcMax = maxval(abs(dWdVc))

    ! In this case dWdVc is equivalent to the condition itself,
    ! so we can use it to measure convergence.
    tConverged = dWdVcMax < this%constrTol

  end subroutine propagateConstraints


  !> Calculate artificial potential to realize constraint on atomic charge.
  subroutine getConstrainEnergyAndPotQ(deltaW, dWdV, Vc, Nc, wAt, wOrb, wSp, qq)

    !> Free energy contribution from current contraint
    real(dp), intent(inout) :: deltaW

    !> Derivative of free energy with respect to potential
    real(dp), intent(out) :: dWdV

    !> Potential / Lagrange multiplier
    real(dp), intent(in) :: Vc

    !> Target population
    real(dp), intent(in) :: Nc

    !> Atom(s) involved in current constrain
    integer, intent(in) :: wAt(:)

    !> Orbital(s) involved in current constrain
    integer, intent(in) :: wOrb(:)

    !> Spin(s) involved in current constrain
    real(dp), intent(in) :: wSp(:,:)

    !> Mulliken populations
    real(dp), intent(in) :: qq(:,:,:)

    !! Number of spin channels to be constrained
    integer :: nSpin

    !! Index of spin channel
    integer :: iSpin

    !! Index of atomic orbital
    integer :: iW

    !! Present population
    real(dp) :: wn

    nSpin = size(wSp, dim=2)
    wn = 0.0_dp
    do iSpin = 1, nSpin
      do iW = 1, size(wAt)
        wn = wn + wSp(iW, iSpin) * qq(wOrb(iW), wAt(iW), iSpin)
      end do
    end do

    dWdV = wn - Nc
    deltaW = deltaW + Vc * dWdV

  end subroutine getConstrainEnergyAndPotQ


  !> Get total shift of all constraints.
  subroutine getConstrainShift(this, shift)

    !> Class instance
    class(TElecConstraint), intent(inout) :: this

    !> Total shift of all constraints
    real(dp), intent(out) :: shift(:,:,:,:)

    !! Iterates over constraints
    integer :: iConstr

    !! Number of constraints requested by the user
    integer :: nConstr

    shift(:,:,:,:) = 0.0_dp
    nConstr = size(this%wAt)

    do iConstr = 1, nConstr
      call getConstrainShiftQ(shift, this%Vc(iConstr), this%wAt(iConstr)%data,&
          & this%wAtOrb(iConstr)%data, this%wAtSpin(iConstr)%data)
    end do

  end subroutine getConstrainShift


  !> Get shift for atomic charge constraint.
  subroutine getConstrainShiftQ(shift, Vc, wAt, wOrb, wSp)

    !> Shift to which contribution is appended
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> Potential / Lagrange multiplier
    real(dp), intent(in) :: Vc

    !> Atom(s) involved in current constrain
    integer, intent(in) :: wAt(:)

    !> Orbital(s) involved in current constrain
    integer, intent(in) :: wOrb(:)

    !> Spin(s) involved in current constrain
    real(dp), intent(in) :: wSp(:,:)

    !! Number of spin channels to be constrained
    integer :: nSpin

    !! Index of spin channel
    integer :: iSpin

    !! Index of atomic orbital
    integer :: iW

    nSpin = size(wSp, dim=2)

    do iSpin = 1, nSpin
      do iW = 1, size(wAt)
        shift(wOrb(iW), wOrb(iW), wAt(iW), iSpin) = shift(wOrb(iW), wOrb(iW), wAt(iW), iSpin)&
            & + Vc * wSp(iW, iSpin)
      end do
    end do

  end subroutine getConstrainShiftQ

end module dftbp_dftb_elecconstraints
