!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Proviedes DFTB+ API for MM-type high level access
module dftbp_mmapi
  use iso_fortran_env, only : output_unit
  use dftbp_accuracy
  use dftbp_globalenv
  use dftbp_environment
  use dftbp_mainapi
  use dftbp_parser
  use dftbp_hsdutils
  use dftbp_inputdata_module
  use dftbp_xmlf90
  implicit none
  private

  public :: TDftbPlus
  public :: TDftbPlus_init, TDftbPlus_destruct
  public :: TDftbPlusInput
  public :: getMaxAngFromSlakoFile
  public :: getPointChargePotential, getPointChargeGradients

  integer :: nDftbPlusCalc = 0


  type :: TDftbPlusInput
    type(fnode), pointer :: hsdTree => null()
  contains
    procedure :: getRootNode => TDftbPlusInput_getRootNode
  end type TDftbPlusInput

  
  type :: TDftbPlus
    private
    type(TEnvironment) :: env
  contains
    procedure, nopass :: getInputFromFile => TDftbPlus_getInputFromFile
    procedure, nopass :: getEmptyInput => TDftbPlus_getEmptyInput
    procedure :: setupCalculator => TDftbPlus_setupCalculator
    procedure :: setGeometry => TDftbPlus_setGeometry
    procedure :: setExternalPotential => TDftbPlus_setExternalPotential
    procedure :: getEnergy => TDftbPlus_getEnergy
    procedure :: getGradients => TDftbPlus_getGradients
    procedure :: getGrossCharges => TDftbPlus_getGrossCharges
  end type TDftbPlus



contains

  subroutine TDftbPlusInput_getRootNode(this, root)
    class(TDftbPlusInput), intent(in) :: this
    type(fnode), pointer, intent(out) :: root

    if (.not. associated(this%hsdTree)) then
      print *, 'ERROR: input has not been created yet!'
      stop
    end if
    call getChild(this%hsdTree, rootTag, root)

  end subroutine TDftbPlusInput_getRootNode
    

  subroutine TDftbPlus_init(this, outputUnit, mpiComm)
    type(TDftbPlus), intent(out) :: this
    integer, intent(in), optional :: outputUnit
    integer, intent(in), optional :: mpiComm

    integer :: stdOut

    if (present(outputUnit)) then
      stdOut = outputUnit
    else
      stdOut = output_unit
    end if

    if (nDftbPlusCalc /= 0) then
      write(stdOut, "(A)") "Error: only one DFTB+ instance is allowed"
      stop
    end if
    nDftbPlusCalc = 1

    call initGlobalEnv(outputUnit=outputUnit, mpiComm=mpiComm)
    call TEnvironment_init(this%env)

  end subroutine TDftbPlus_init


  subroutine TDftbPlus_destruct(this)
    type(TDftbPlus), intent(inout) :: this

    call this%env%destruct()
    call destructGlobalEnv()
    nDftbPlusCalc = 0

  end subroutine TDftbPlus_destruct


  subroutine TDftbPlus_getInputFromFile(fileName, input)
    character(*), intent(in) :: fileName
    type(TDftbPlusInput), intent(out) :: input

    call readHsdFile(fileName, input%hsdTree)

  end subroutine TDftbPlus_getInputFromFile


  subroutine TDftbPlus_getEmptyInput(input)
    type(TDftbPlusInput), intent(out) :: input

    type(fnode), pointer :: root, dummy

    input%hsdTree => createDocumentNode()
    root => createElement(rootTag)
    dummy => appendChild(input%hsdTree, root)
    
  end subroutine TDftbPlus_getEmptyInput


  subroutine TDftbPlus_setupCalculator(this, input)
    class(TDftbPlus), intent(inout) :: this
    type(TDftbPlusInput), intent(inout) :: input

    type(TParserFlags) :: parserFlags
    type(inputData) :: inpData

    call parseHsdTree(input%hsdTree, inpData, parserFlags)
    call initProgramVariables(inpData, this%env)

  end subroutine TDftbPlus_setupCalculator


  subroutine TDftbPlus_setGeometry(this, coords, latVecs)
    class(TDftbPlus), intent(inout) :: this
    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(in), optional :: latVecs(:,:)

    call setGeometry(coords, latVecs)

  end subroutine TDftbPlus_setGeometry


  subroutine TDftbPlus_setExternalPotential(this, atomPot, shellPot, potGrad)
    class(TDftbPlus), intent(inout) :: this
    real(dp), intent(in), optional :: atomPot(:)
    real(dp), intent(in), optional :: shellPot(:,:)
    real(dp), intent(in), optional :: potGrad(:,:)

    call setExternalPotential(atomPot, shellPot, potGrad)

  end subroutine TDftbPlus_setExternalPotential

  
  subroutine TDftbPlus_getEnergy(this, merminEnergy)
    class(TDftbPlus), intent(inout) :: this
    real(dp), intent(out) :: merminEnergy

    call getEnergy(this%env, merminEnergy)
    
  end subroutine TDftbPlus_getEnergy


  subroutine TDftbPlus_getGradients(this, gradients)
    class(TDftbPlus), intent(inout) :: this
    real(dp), intent(out) :: gradients(:,:)

    call getGradients(this%env, gradients)
    
  end subroutine TDftbPlus_getGradients


  subroutine TDftbPlus_getGrossCharges(this, atomCharges)
    class(TDftbPlus), intent(inout) :: this
    real(dp), intent(out) :: atomCharges(:)

    call getGrossCharges(atomCharges)

  end subroutine TDftbPlus_getGrossCharges


  function getMaxAngFromSlakoFile(slakoFile) result(maxAng)
    character(*), intent(in) :: slakoFile
    integer :: maxAng

    real(dp) :: dr
    integer :: nGridPoints, nShells
    integer :: fd

    open(newunit=fd, file=slakoFile, status="old", action="read")
    read(fd, *) dr, nGridPoints, nShells
    close(fd)
    maxAng = nShells - 1

  end function getMaxAngFromSlakoFile


  subroutine getPointChargePotential(coordsMm, chargesMm, coordsQm, extPot, extPotGrad)
    real(dp), intent(in) :: coordsMm(:,:)
    real(dp), intent(in) :: chargesMm(:)
    real(dp), intent(in) :: coordsQm(:,:)
    real(dp), intent(out) :: extPot(:)
    real(dp), intent(out) :: extPotGrad(:,:)

    real(dp) :: atomPosQm(3), atomPosMm(3)
    real(dp) :: chargeMm, dist
    integer :: nAtomQm, nAtomMm
    integer :: iAtQm, iAtMm

    nAtomQm = size(coordsQm, dim=2)
    nAtomMm = size(coordsMm, dim=2)
    extPot(:) = 0.0_dp
    extPotGrad(:,:) = 0.0_dp
    do iAtQm = 1, nAtomQm
      atomPosQm(:) = coordsQm(:, iAtQm)
      do iAtMm = 1, nAtomMm
        atomPosMm(:) = coordsMm(1:3, iAtMm)
        chargeMm = chargesMm(iAtMm)
        dist = sqrt(sum((atomPosQm - atomPosMm)**2))
        extPot(iAtQm) = extPot(iAtQm) + chargeMm / dist
        extPotGrad(:, iAtQm) = extPotGrad(:, iAtQm) - chargeMm * (atomPosQm - atomPosMm) / dist**3
      end do
    end do

  end subroutine getPointChargePotential


  subroutine getPointChargeGradients(coordsQm, chargesQm, coordsMm, chargesMm, gradients)
    real(dp), intent(in) :: coordsQm(:,:)
    real(dp), intent(in) :: chargesQm(:)
    real(dp), intent(in) :: coordsMm(:,:)
    real(dp), intent(in) :: chargesMm(:)
    real(dp), intent(out) :: gradients(:,:)

    real(dp) :: atomPosQm(3), atomPosMm(3)
    real(dp) :: chargeQm, chargeMm, dist
    integer :: nAtomQm, nAtomMm
    integer :: iAtQm, iAtMm

    gradients(:,:) = 0.0_dp
    nAtomQm = size(coordsQm, dim=2)
    nAtomMm = size(coordsMm, dim=2)
    do iAtMm = 1, nAtomMm
      atomPosMm(:) = coordsMm(:, iAtMm)
      chargeMm = chargesMm(iAtMm)
      do iAtQm = 1, nAtomQm
        atomPosQm(:) = coordsQm(:, iAtQm)
        chargeQm = chargesQm(iAtQm)
        dist = sqrt(sum((atomPosQm - atomPosMm)**2))
        gradients(:, iAtMm) = gradients(:, iAtMm) &
            & - chargeQm * chargeMm * (atomPosMm - atomPosQm) / dist**3
      end do
    end do
    
  end subroutine getPointChargeGradients


end module dftbp_mmapi
