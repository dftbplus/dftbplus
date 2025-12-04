!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

!> Module to impose constraints on the electronic ground state.
module dftbp_dftb_elecconstraints
  use dftbp_common_accuracy, only : dp, hugeIterations
  use dftbp_dftbplus_input_geoopt, only : readOptimizerInput
  use dftbp_extlibs_xmlf90, only : char, destroyNodeList, fnode, fnodeList, getItem1, getLength,&
      & string
  use dftbp_geoopt_package, only : createOptimizer, TOptimizer, TOptimizerInput
  use dftbp_io_hsdutils, only : detailedError, getChild, getChildren, getChildValue,&
      & getSelectedAtomIndices
  use dftbp_io_hsdutils2, only : renameChildren
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_typegeometry, only : TGeometry
  use dftbp_type_wrappedintr, only : TWrappedInt1, TWrappedReal2
  implicit none (type, external)

  private
  public :: TElecConstraint, TElecConstraint_init, TElecConstraintInp
  public :: readElecConstraintInput


  !> Represents the input for a single Mulliken population constraint
  type TMullikenConstrInp

    !> Atoms the constraint applies to
    integer, allocatable :: atoms(:)

    !> Constraint values (either single value = total or as as many as atoms = individual)
    real(dp), allocatable :: constrValues(:)

    !> Spin channel factors (1.0 for channels included in population calculations, 0.0 for rest)
    real(dp), allocatable :: spinChannelFactors(:)

    !> Whether constraint values are charges (otherwise they are populations)
    logical :: constrValuesAreCharges = .false.

  end type TMullikenConstrInp


  !> Contains input data for electronic constraints.
  type TElecConstraintInp

    !> Input for Mulliken population constraints
    type(TMullikenConstrInp), allocatable :: mullikenConstrs(:)

    !> Optimiser input choice
    class(TOptimizerInput), allocatable :: optimiser

    !> Derivative (absolute) tolerance for constraint
    real(dp) :: constrTol

    !> Number of iterations for enforcing constraint
    integer :: nConstrIter

    !> True, if converged micro-iterations are required
    logical :: isConvRequired

  end type TElecConstraintInp


  !> Collects data to apply all Mulliken constraints.
  type TMullikenConstr
    private

    !> Number of constraints
    integer :: nConstr

    !> Target value of the constraints
    real(dp), allocatable :: Nc(:)

    !> Atom(s) involved in each constraint
    type(TWrappedInt1), allocatable :: wAt(:)

    !> Atomic orbital(s) involved in each constraint
    type(TWrappedInt1), allocatable :: wAtOrb(:)

    !> Atomic orbital qm-values involved in each constraint
    !> ([q] for spin unpolarized case, [q, m] for colinear spin)
    type(TWrappedReal2), allocatable :: wAtSpin(:)

  contains

    !> Calculates Mulliken population constraints
    procedure :: getEnergyAndPot => TMullikenConstr_getEnergyAndPot

    !> Adds shift to Mulliken population constraints
    procedure :: addShift => TMullikenConstr_addShift

    !> Get nr. of constraints
    procedure :: getNConstr => TMullikenConstr_getNConstr

  end type TMullikenConstr


  !> Represents electronic constraints.
  type TElecConstraint
    private

    !> General optimiser
    class(TOptimizer), allocatable :: optimizer

    !> True, if converged micro-iterations are required
    logical :: isConvRequired = .false.

    !> Mulliken population constraints
    type(TMullikenConstr) :: mullikenConstr

    !> Number of constraints
    integer :: nConstr = 0

    !> Potential created by the constraints
    real(dp), allocatable :: Vc(:)

    !> Contribution to free energy functional from constraints
    real(dp), allocatable :: deltaW(:)

    !> Derivative of energy functional with respect to Vc
    real(dp), allocatable :: dWdVc(:)

    !> Derivative tolerance for constraint
    real(dp) :: constrTol = 0.0_dp

    !> Number of iterations for enforcing constraint
    integer :: nConstrIter = 0

  contains

    !> Returns constraining potential (Vc)
    procedure :: getConstraintShift => TElecConstraint_getConstraintShift

    !> Update constraints for current state of the system
    procedure :: propagateConstraints => TElecConstraint_propagateConstraints

    !> Maximum number of iterations for driver
    procedure :: getMaxIter => TElecConstraint_getMaxIter

    !> Free energy (W) for system in constraining potential(s)
    procedure :: getFreeEnergy => TElecConstraint_getFreeEnergy

    !> Maximum component of free energy derivative wrt. constraint potential(s)
    procedure :: getMaxEnergyDerivWrtVc => TElecConstraint_getMaxEnergyDerivWrtVc

    !> Reset optimizer
    procedure :: resetOptimizer => TElecConstraint_resetOptimizer

    !> Whether constraints require converged micro-iterations
    procedure :: requiresConvergence => TElecConstraint_requiresConvergence

  end type TElecConstraint

contains


  !> Reads electronic constraint input from HSD.
  subroutine readElecConstraintInput(node, geo, isSpinPol, is2Component, input)

    !> Input structure to be filled
    type(TElecConstraintInp), intent(out) :: input

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> True, if this is a spin polarized calculation
    logical, intent(in) :: isSpinPol

    !> Is this a calculation with Pauli wavefunctions
    logical, intent(in) :: is2Component

    type(fnode), pointer :: constrContainer, dummyNode, child1

    call renameChildren(node, "Optimizer", "Optimiser")
    call getChildValue(node, "Optimiser", child1, "FIRE")
    call readOptimizerInput(child1, input%optimiser)

    call getChildValue(node, "ConstrTolerance", input%constrTol, 1.0e-08_dp)
    call getChildValue(node, "MaxConstrIterations", input%nConstrIter, 100)
    call getChildValue(node, "ConvergentConstrOnly", input%isConvRequired, .true.)

    call getChildValue(node, "Constraints", dummyNode, "", child=constrContainer,&
        & allowEmptyValue=.true., dummyValue=.true., list=.true.)
    call readMullikenConstraintInputs(constrContainer, geo, isSpinPol, is2Component,&
        & input%mullikenConstrs)

  end subroutine readElecConstraintInput


  !> Reads Mulliken constraint inputs from HSD.
  subroutine readMullikenConstraintInputs(constrContainer, geo, isSpinPol, is2Component, inputs)

    !> Node containing all constraints
    type(fnode), pointer, intent(in) :: constrContainer

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> True, if this is a spin polarized calculation
    logical, intent(in) :: isSpinPol

    !> Is this a calculation with Pauli wavefunctions
    logical, intent(in) :: is2Component

    !> Array of input structures (depending on number of Mulliken constraints defined)
    type(TMullikenConstrInp), allocatable, intent(out) :: inputs(:)

    type(fnodeList), pointer :: constrNodes
    type(fnode), pointer :: constrNode, child1
    type(fnode), pointer :: totalPopNode, populationsNode, totalChargeNode, chargesNode
    type(string) :: buffer
    real(dp) :: rTmp
    integer :: iConstrInp, nConstrInp, nAssociated

    call getChildren(constrContainer, "MullikenPopulation", constrNodes)
    if (.not. associated(constrNodes)) return

    nConstrInp = getLength(constrNodes)
    allocate(inputs(nConstrInp))

    do iConstrInp = 1, nConstrInp
      associate(input => inputs(iConstrInp))
        call getItem1(constrNodes, iConstrInp, constrNode)
        call getChildValue(constrNode, "Atoms", buffer, child=child1, multiple=.true.)
        call getSelectedAtomIndices(child1, char(buffer), geo%speciesNames, geo%species,&
            & input%atoms)

        call getChild(constrNode, "Populations", populationsNode, requested=.false.)
        call getChild(constrNode, "TotalPopulation", totalPopNode, requested=.false.)
        call getChild(constrNode, "Charges", chargesNode, requested=.false.)
        call getChild(constrNode, "TotalCharge", totalChargeNode, requested=.false.)
        nAssociated = count([associated(populationsNode), associated(totalPopNode), &
            & associated(chargesNode), associated(totalChargeNode)])
        if (nAssociated /= 1) then
          call detailedError(constrNode, "You must specify exactly one and only one of the options&
              & Populations, TotalPopulation, Charges or TotalCharge")
        end if
        input%constrValuesAreCharges = associated(chargesNode) .or. associated(totalChargeNode)
        if (associated(populationsNode) .or. associated(chargesNode)) then
          allocate(input%constrValues(size(input%atoms)), source=0.0_dp)
          if (associated(populationsNode)) then
            call getChildValue(populationsNode, "", input%constrValues)
          else
            call getChildValue(chargesNode, "", input%constrValues)
          end if
        else
          if (associated(totalPopNode)) then
            call getChildValue(totalPopNode, "", rTmp)
          else
            call getChildValue(totalChargeNode, "", rTmp)
          end if
          input%constrValues = [rTmp]
        end if

        ! Functionality currently restricted to charge channel only
        if (isSpinPol) then
          ! [q, m] representation
          if (is2Component) then
            input%spinChannelFactors = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
          else
            input%spinChannelFactors = [1.0_dp, 0.0_dp]
          end if
        else
          input%spinChannelFactors = [1.0_dp]
        end if

      end associate
    end do

    call destroyNodeList(constrNodes)

  end subroutine readMullikenConstraintInputs


  !> Initialises the constraints structure.
  subroutine TElecConstraint_init(this, input, orb, q0)

    !> Constraint structure instance
    type(TElecConstraint), intent(out) :: this

    !> Input data structure
    type(TElecConstraintInp), intent(inout) :: input

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Reference Mulliken populations
    real(dp), intent(in) :: q0(:,:,:)

    call TMullikenConstr_init(this%mullikenConstr, input%mullikenConstrs, orb, q0)
    ! Currently we only support Mulliken constraints
    this%nConstr = this%mullikenConstr%getNConstr()

    ! We should enable optional initialization of Vc from input at some point.
    allocate(this%Vc(this%nConstr), source=0.0_dp)
    allocate(this%dWdVc(this%nConstr), source=0.0_dp)
    allocate(this%deltaW(this%nConstr), source=0.0_dp)

    call createOptimizer(input%optimiser, this%nConstr, this%optimizer)

    if (input%nConstrIter == -1) then
      this%nConstrIter = hugeIterations
    else
      this%nConstrIter = input%nConstrIter
    end if
    this%isConvRequired = input%isConvRequired
    this%constrTol = input%constrTol

  end subroutine TElecConstraint_init


  !> Returns maximum number of iterations for constraint driver.
  pure function TElecConstraint_getMaxIter(this) result(maxIter)

    !> Class instance
    class(TElecConstraint), intent(in) :: this

    !> Obtained maximum number of iterations
    integer :: maxIter

    maxIter = this%nConstrIter

  end function TElecConstraint_getMaxIter


  !> Returns total contribution to free energy functional from constraints.
  pure function TElecConstraint_getFreeEnergy(this) result(deltaWTotal)

    !> Class instance
    class(TElecConstraint), intent(in) :: this

    !> Summed up contribution to free energy functional from constraints
    real(dp) :: deltaWTotal

    deltaWTotal = sum(this%deltaW)

  end function TElecConstraint_getFreeEnergy


  !> Returns the maximum absolute value of dW/dVc across all constraints
  pure function TElecConstraint_getMaxEnergyDerivWrtVc(this) result(dWdVcMax)

    !> Class instance
    class(TElecConstraint), intent(in) :: this

    !> Maximum derivative of energy functional with respect to Vc
    real(dp) :: dWdVcMax

    dWdVcMax = maxval(abs(this%dWdVc))

  end function TElecConstraint_getMaxEnergyDerivWrtVc


  !> Applies electronic constraints to system.
  subroutine TElecConstraint_propagateConstraints(this, qq, energy, tConverged)

    !> Class instance
    class(TElecConstraint), intent(inout) :: this

    !> Mulliken populations
    real(dp), intent(in) :: qq(:,:,:)

    !> Energy
    real(dp), intent(in) :: energy

    !> Gradient convergence achieved
    logical, intent(out) :: tConverged

    ! Summed up contribution to free energy functional from constraints
    real(dp) :: deltaWTotal

    ! Maximum derivative of energy functional with respect to Vc
    real(dp) :: dWdVcMax

    ! Potential displacement proposed by optimizer
    real(dp) :: potDisplace(size(this%Vc))

    call this%mullikenConstr%getEnergyAndPot(this%Vc, qq, this%deltaW, this%dWdVc)

    ! Sum up all free energy contributions
    deltaWTotal = this%getFreeEnergy()

    ! Get maximum derivative of energy functional with respect to Vc
    dWdVcMax = this%getMaxEnergyDerivWrtVc()

    call this%optimizer%step(energy + deltaWTotal, -this%dWdVc, potDisplace)
    this%Vc(:) = this%Vc + potDisplace

    ! In this case dWdVc is equivalent to the condition itself,
    ! so we can use it to measure convergence.
    tConverged = dWdVcMax < this%constrTol

  end subroutine TElecConstraint_propagateConstraints


  !> Get total shift of all constraints.
  subroutine TElecConstraint_getConstraintShift(this, shift)

    !> Class instance
    class(TElecConstraint), intent(inout) :: this

    !> Total shift of all constraints
    real(dp), intent(out) :: shift(:,:,:,:)

    shift(:,:,:,:) = 0.0_dp
    call this%mullikenConstr%addShift(this%Vc(1:this%nConstr), shift)

  end subroutine TElecConstraint_getConstraintShift


  !> Resets the optimizer.
  subroutine TElecConstraint_resetOptimizer(this)

    !> Class instance
    class(TElecConstraint), intent(inout) :: this

    call this%optimizer%reset()

  end subroutine TElecConstraint_resetOptimizer


  !> Whether constraints require converged micro-iterations.
  pure function TElecConstraint_requiresConvergence(this) result(reqConv)

    !> Class instance
    class(TElecConstraint), intent(in) :: this

    !> True, if converged micro-iterations are required
    logical :: reqConv

    reqConv = this%isConvRequired

  end function TElecConstraint_requiresConvergence


  !> Initializes constraint helper arrays from Mulliken constraints.
  subroutine TMullikenConstr_init(this, inputs, orb, q0)

    !> Class instance
    type(TMullikenConstr), intent(out) :: this

    !> Mulliken population constraint inputs
    type(TMullikenConstrInp), intent(in) :: inputs(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Reference population for each orbital in the system
    real(dp), intent(in) :: q0(:,:,:)

    integer, allocatable :: atoms(:)
    integer :: nConstrInputs, nAllConstr, nConstrVal, nOrb, nOrbAtom, nSpin
    integer :: iConstrInp, constrInd, iConstrVal, iAt, iOrb, ii

    ! Count number of constraints
    nAllConstr = 0
    nConstrInputs = size(inputs)
    do iConstrInp = 1, nConstrInputs
      nAllConstr = nAllConstr + size(inputs(iConstrInp)%constrValues)
    end do

    allocate(this%Nc(nAllConstr))
    allocate(this%wAt(nAllConstr))
    allocate(this%wAtOrb(nAllConstr))
    allocate(this%wAtSpin(nAllConstr))

    ! Allocate + initialize arrays and build index mappings
    constrInd = 0
    do iConstrInp = 1, nConstrInputs
      associate (input => inputs(iConstrInp))
        nConstrVal = size(input%constrValues)
        do iConstrVal = 1, nConstrVal
          constrInd = constrInd + 1
          ! If only one constraint value had been specified, it represents the sum over all atoms in
          ! the constraint. Otherwise they are values for individual constraints on single atoms.
          if (nConstrVal == 1) then
            atoms = input%atoms
          else
            atoms = [input%atoms(iConstrVal)]
          end if
          nOrb = sum(orb%nOrbAtom(atoms))
          allocate(this%wAt(constrInd)%data(nOrb), source=0)
          allocate(this%wAtOrb(constrInd)%data(nOrb), source=0)
          nSpin = size(input%spinChannelFactors)
          allocate(this%wAtSpin(constrInd)%data(nOrb, nSpin), source=0.0_dp)

          this%Nc(constrInd) = input%constrValues(iConstrVal)
          if (input%constrValuesAreCharges) then
            this%Nc(constrInd) = sum(q0(:, atoms, :)) - this%Nc(constrInd)
          end if
          nOrb = 0
          do ii = 1, size(atoms)
            iAt = atoms(ii)
            nOrbAtom = orb%nOrbAtom(iAt)
            do iOrb = 1, nOrbAtom
              this%wAt(constrInd)%data(nOrb + iOrb) = iAt
              this%wAtOrb(constrInd)%data(nOrb + iOrb) = iOrb
              this%wAtSpin(constrInd)%data(nOrb + iOrb, :) = input%spinChannelFactors
            end do
            nOrb = nOrb + nOrbAtom
          end do
        end do
      end associate
    end do

    this%nConstr = nAllConstr

  end subroutine TMullikenConstr_init


  !> Calculate Mulliken population constraints.
  subroutine TMullikenConstr_getEnergyAndPot(this, Vc, qq, deltaW, dWdV)

    !> Class instance
    class(TMullikenConstr), intent(in) :: this

    !> Constraint potential
    real(dp), intent(in) :: Vc(:)

    !> Mulliken populations
    real(dp), intent(in) :: qq(:,:,:)

    !> Contribution to free energy functional from constraints
    real(dp), intent(out) :: deltaW(:)

    !> Derivative of energy functional with respect to Vc
    real(dp), intent(out) :: dWdV(:)

    integer :: iConstr

    @:ASSERT(size(Vc) == this%nConstr)
    @:ASSERT(size(deltaW) == this%nConstr)
    @:ASSERT(size(dWdV) == this%nConstr)

    do iConstr = 1, this%nConstr
      call getConstraintEnergyAndPotQ_(Vc(iConstr), this%Nc(iConstr), this%wAt(iConstr)%data,&
          & this%wAtOrb(iConstr)%data, this%wAtSpin(iConstr)%data, qq, deltaW(iConstr),&
          & dWdV(iConstr))
    end do

  end subroutine TMullikenConstr_getEnergyAndPot


  !> Add shift to Mulliken population constraints.
  subroutine TMullikenConstr_addShift(this, Vc, shift)

    !> Class instance
    class(TMullikenConstr), intent(in) :: this

    !> Constraint potential
    real(dp), intent(in) :: Vc(:)

    !> Shift to be updated
    real(dp), intent(inout) :: shift(:,:,:,:)

    integer :: iConstr

    @:ASSERT(size(Vc) == this%nConstr)

    do iConstr = 1, this%nConstr
      call addConstraintsShiftQ_(shift, Vc(iConstr), this%wAt(iConstr)%data,&
          & this%wAtOrb(iConstr)%data, this%wAtSpin(iConstr)%data)
    end do

  end subroutine TMullikenConstr_addShift


  !> Returns number of constraints.
  function TMullikenConstr_getNConstr(this) result(nConstr)

    !> Class instance
    class(TMullikenConstr), intent(in) :: this

    !> Number of constraints
    integer :: nConstr

    nConstr = this%nConstr

  end function TMullikenConstr_getNConstr


  !> Calculate artificial potential to realize constraint(s) on atomic charge.
  subroutine getConstraintEnergyAndPotQ_(Vc, Nc, wAt, wOrb, wSp, qq, deltaW, dWdV)

    !> Constraint potential value
    real(dp), intent(in) :: Vc

    !> Constraint population
    real(dp), intent(in) :: Nc

    !> Atoms involved in the constraint
    integer, intent(in) :: wAt(:)

    !> Orbitals involved in the constraint
    integer, intent(in) :: wOrb(:)

    !> Spin channel factors
    real(dp), intent(in) :: wSp(:,:)

    !> Mulliken populations
    real(dp), intent(in) :: qq(:,:,:)

    !> Contribution to free energy functional from constraints
    real(dp), intent(out) :: deltaW

    !> Derivative of energy functional with respect to Vc
    real(dp), intent(out) :: dWdV

    integer :: nSpin, iSpin, iW
    real(dp) :: wn

    nSpin = size(wSp, dim=2)
    wn = 0.0_dp
    do iSpin = 1, nSpin
      do iW = 1, size(wAt)
        wn = wn + wSp(iW, iSpin) * qq(wOrb(iW), wAt(iW), iSpin)
      end do
    end do

    dWdV = wn - Nc
    deltaW = Vc * dWdV

  end subroutine getConstraintEnergyAndPotQ_


  !> Get shift for atomic charge constraint.
  subroutine addConstraintsShiftQ_(shift, Vc, wAt, wOrb, wSp)

    !> Shift vector to modify
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> Constraint potential
    real(dp), intent(in) :: Vc

    !> Atoms involved in the constraint
    integer, intent(in) :: wAt(:)

    !> Orbitals involved in the constraint
    integer, intent(in) :: wOrb(:)

    !> Spin channel factors
    real(dp), intent(in) :: wSp(:,:)

    integer :: nSpin, iSpin, iW

    nSpin = size(wSp, dim=2)

    do iSpin = 1, nSpin
      do iW = 1, size(wAt)
        shift(wOrb(iW), wOrb(iW), wAt(iW), iSpin) = shift(wOrb(iW), wOrb(iW), wAt(iW), iSpin)&
            & + Vc * wSp(iW, iSpin)
      end do
    end do

  end subroutine addConstraintsShiftQ_

end module dftbp_dftb_elecconstraints
