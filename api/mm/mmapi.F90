!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Proviedes DFTB+ API for MM-type high level access
module mmapi
  use accuracy
  use iso_fortran_env, only : output_unit
  use globalenv
  use environment
  use mainapi
  use parser
  use hsdutils
  use inputData_module
  use xmlf90
  implicit none
  private

  public :: TDftbPlus
  public :: TDftbPlus_init, TDftbPlus_destruct
  public :: TDftbPlusInput

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


end module mmapi
