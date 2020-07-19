!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains the C-API of DFTB+.
module dftbp_capi
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env
  use dftbp_accuracy, only : dp
  use dftbp_mmapi, only : TDftbPlus, TDftbPlus_init, TDftbPlus_destruct, TDftbPlusInput
  use dftbp_qdepextpotgenc, only :&
      & getExtPotIfaceC, getExtPotGradIfaceC, TQDepExtPotGenC, TQDepExtPotGenC_init
  implicit none
  private

  !> DFTB+ input tree
  type, bind(C) :: c_DftbPlusInput
    type(c_ptr) :: pDftbPlusInput
  end type c_DftbPlusInput

  !> DFTB+ calculation
  type, bind(C) :: c_DftbPlus
    type(c_ptr) :: instance
  end type c_DftbPlus


  !> Simple extension around the TDftbPlus with some additional variables for the C-API.
  type, extends(TDftbPlus) :: TDftbPlusC
    private
    integer :: outputUnit
    logical :: tOutputOpened
  end type TDftbPlusC


contains


  !> Returns the current API version
  subroutine c_DftbPlus_api(major, minor, patch) bind(C, name='dftbp_api')

    !> makor.minor.patch
    integer(c_int), intent(out) :: major, minor, patch

    major = ${APIMAJOR}$
    minor = ${APIMINOR}$
    patch = ${APIPATCH}$

  end subroutine c_DftbPlus_api


  !> Initialises a DFTB+ calculation with output sent to to some location
  subroutine c_DftbPlus_init(handler, outputFileName) bind(C, name='dftbp_init')

    !> DFTB+ handler
    type(c_DftbPlus), intent(out) :: handler

    !> output location
    type(c_ptr), value, intent(in) :: outputFileName

    type(TDftbPlusC), pointer :: instance
    character(c_char), pointer :: pOutputFileName
    character(:), allocatable :: fortranFileName

    allocate(instance)
    if (c_associated(outputFileName)) then
      call c_f_pointer(outputFileName, pOutputFileName)
      fortranFileName = fortranChar(pOutputFileName)
      open(newunit=instance%outputUnit, file=fortranFileName, action="write")
      instance%tOutputOpened = .true.
    else
      instance%outputUnit = output_unit
      instance%tOutputOpened = .false.
    end if
    call TDftbPlus_init(instance%TDftbPlus, outputUnit=instance%outputUnit)
    handler%instance = c_loc(instance)

  end subroutine c_DftbPlus_init


  !> finalises a DFTB+ instance
  subroutine c_DftbPlus_final(handler) bind(C, name='dftbp_final')

    !> DFTB+ handler
    type(c_DftbPlus), intent(inout) :: handler

    !> the specific instance to be finalised
    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    call TDftbPlus_destruct(instance%TDftbPlus)
    if (instance%tOutputOpened) then
      close(instance%outputUnit)
    end if
    deallocate(instance)
    handler%instance = c_null_ptr

  end subroutine c_DftbPlus_final


  !> Read input for DFTB+ from a specified file
  subroutine c_DftbPlus_getInputFromFile(handler, fileName, inputHandler)&
      & bind(C, name='dftbp_get_input_from_file')

    !> handler for the input
    type(c_DftbPlus), intent(inout) :: handler

    !> file to read
    character(c_char), intent(in) :: fileName(*)

    !> handler for the resulting input
    type(c_DftbPlusInput), intent(out) :: inputHandler

    type(TDftbPlusC), pointer :: instance
    type(TDftbPlusInput), pointer :: pDftbPlusInput
    character(:), allocatable :: fortranFileName

    allocate(pDftbPlusInput)
    fortranFileName = fortranChar(fileName)
    call c_f_pointer(handler%instance, instance)
    call instance%getInputFromFile(fortranFileName, pDftbPlusInput)
    inputHandler%pDftbPlusInput = c_loc(pDftbPlusInput)

  end subroutine c_DftbPlus_getInputFromFile


  !> process a document tree to get settings for the calculation
  subroutine c_DftbPlus_processInput(handler, inputHandler)&
      & bind(C, name='dftbp_process_input')

    !> handler for the calculation instance
    type(c_DftbPlus), intent(inout) :: handler

    !> input tree handler
    type(c_DftbPlusInput), intent(inout) :: inputHandler

    type(TDftbPlusC), pointer :: instance
    type(TDftbPlusInput), pointer :: pDftbPlusInput

    call c_f_pointer(handler%instance, instance)
    call c_f_pointer(inputHandler%pDftbPlusInput, pDftbPlusInput)
    call instance%setupCalculator(pDftbPlusInput)

  end subroutine c_DftbPlus_processInput


  !> set an external potential on the DFTB+ calculation
  subroutine c_DftbPlus_setExternalPotential(handler, extPot, extPotGrad)&
      & bind(C, name='dftbp_set_external_potential')

    !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> externally set potential
    real(c_double), intent(in) :: extPot(*)

    !> gradient of the potential wrt to atom positions
    type(c_ptr), value, intent(in) :: extPotGrad

    type(TDftbPlusC), pointer :: instance
    real(c_double), pointer :: pExtPotGrad(:,:)
    integer :: nAtom

    call c_f_pointer(handler%instance, instance)
    nAtom = instance%nrOfAtoms()

    if (c_associated(extPotGrad)) then
      call c_f_pointer(extPotGrad, pExtPotGrad, [3, nAtom])
    else
      pExtPotGrad => null()
    end if
    call instance%setExternalPotential(extPot(1:nAtom), pExtPotGrad)

  end subroutine c_DftbPlus_setExternalPotential


  !> register a generator for an external potential
  subroutine c_DftbPlus_registerExtPotGenerator(handler, refPtr, extPotFunc, extPotGradFunc)&
      & bind(C, name='dftbp_register_ext_pot_generator')

    !> handler for the potential
    type(c_DftbPlus), intent(inout) :: handler

    !> pointer to the C routine for the external potential
    type(c_ptr), value, intent(in) :: refPtr

    !> function for the external potential
    type(c_funptr), value, intent(in) :: extPotFunc

    !> function for the gradient of the potential
    type(c_funptr), value, intent(in) :: extPotGradFunc

    type(TDftbPlusC), pointer :: instance
    type(TQDepExtPotGenC) :: extPotGenC
    procedure(getExtPotIfaceC), pointer :: pExtPotFunc
    procedure(getExtPotGradIfaceC), pointer :: pExtPotGradFunc

    call c_f_procpointer(extPotFunc, pExtPotFunc)
    call c_f_procpointer(extPotGradFunc, pExtPotGradFunc)
    call c_f_pointer(handler%instance, instance)
    call TQDepExtPotGenC_init(extPotGenC, refPtr, pExtPotFunc, pExtPotGradFunc)
    call instance%setQDepExtPotGen(extPotGenC)

  end subroutine c_DftbPlus_registerExtPotGenerator


  !> set/replace the coordinates in a DFTB+ calculation instance
  subroutine c_DftbPlus_setCoords(handler, coords) bind(C, name='dftbp_set_coords')

    !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> coordinates, (xyz, :nAtom)
    real(c_double), intent(in) :: coords(3,*)

    type(TDftbPlusC), pointer :: instance
    integer :: nAtom

    call c_f_pointer(handler%instance, instance)
    nAtom = instance%nrOfAtoms()
    call instance%setGeometry(coords(:, 1:nAtom))

  end subroutine c_DftbPlus_setCoords


  !> Set both the coordinates and lattice vectors
  subroutine c_DftbPlus_setCoordsAndLatticeVecs(handler, coords, latVecs)&
      & bind(C, name='dftbp_set_coords_and_lattice_vecs')

    !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> coordinates, row major format (xyz, :nAtom)
    real(c_double), intent(in) :: coords(3,*)

    !> lattice vectors, row major format
    real(c_double), intent(in) :: latvecs(3, *)

    type(TDftbPlusC), pointer :: instance
    integer :: nAtom

    call c_f_pointer(handler%instance, instance)
    nAtom = instance%nrOfAtoms()
    call instance%setGeometry(coords(:, 1:nAtom), latVecs(:, 1:3))

  end subroutine c_DftbPlus_setCoordsAndLatticeVecs


  !> Set the coordinates and lattice vectors with an origin
  subroutine c_DftbPlus_setCoordsLatticeVecsOrigin(handler, coords, latVecs, origin)&
      & bind(C, name='dftbp_set_coords_lattice_origin')

    !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> coordinates, row major format (xyz, :nAtom)
    real(c_double), intent(in) :: coords(3,*)

    !> lattice vectors, row major format
    real(c_double), intent(in) :: latvecs(3, *)

    !> coordinate origin
    real(c_double), intent(in) :: origin(3)

    type(TDftbPlusC), pointer :: instance
    integer :: nAtom

    call c_f_pointer(handler%instance, instance)
    nAtom = instance%nrOfAtoms()
    call instance%setGeometry(coords(:, 1:nAtom), latVecs(:, 1:3), origin(1:3))

  end subroutine c_DftbPlus_setCoordsLatticeVecsOrigin


  !> Obtain nr. of atoms.
  function c_DftbPlus_nrOfAtoms(handler) result(nAtom) bind(C, name='dftbp_get_nr_atoms')
    type(c_DftbPlus), intent(inout) :: handler
    integer(c_int) :: nAtom

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    nAtom = instance%nrOfAtoms()

  end function c_DftbPlus_nrOfAtoms


  !> Obtain the DFTB+ energy
  subroutine c_DftbPlus_getEnergy(handler, merminEnergy) bind(C, name='dftbp_get_energy')

    !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> resulting energy
    real(c_double), intent(out) :: merminEnergy

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    call instance%getEnergy(merminEnergy)

  end subroutine c_DftbPlus_getEnergy


  !> Obtain the gradients wrt DFTB atom positions
  subroutine c_DftbPlus_getGradients(handler, gradients) bind(C, name='dftbp_get_gradients')

    !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> gradients, row major format
    real(c_double), intent(out) :: gradients(3, *)

    type(TDftbPlusC), pointer :: instance
    integer :: nAtom

    call c_f_pointer(handler%instance, instance)
    nAtom = instance%nrOfAtoms()
    call instance%getGradients(gradients(:, 1:nAtom))

  end subroutine c_DftbPlus_getGradients


  !> Obtain the stress tensor of the periodic system
  subroutine c_DftbPlus_getStressTensor(handler, stresstensor)&
      & bind(C, name='dftbp_get_stress_tensor')

    !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> gradients, row major format
    real(c_double), intent(out) :: stresstensor(3, 3)

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)

    call instance%getStressTensor(stresstensor(:, 1:3))

  end subroutine c_DftbPlus_getStressTensor


  !> Obtain gross (Mulliken) charges for atoms wrt to neutral references
  subroutine c_DftbPlus_getGrossCharges(handler, atomCharges)&
      & bind(C, name='dftbp_get_gross_charges')

    !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> resulting atomic charges
    real(c_double), intent(out) :: atomCharges(*)

    type(TDftbPlusC), pointer :: instance
    integer :: nAtom

    call c_f_pointer(handler%instance, instance)
    nAtom = instance%nrOfAtoms()
    call instance%getGrossCharges(atomCharges(1:nAtom))

  end subroutine c_DftbPlus_getGrossCharges


  !> Converts a 0-char terminated C-type string into a Fortran string.
  function fortranChar(cstring, maxlen)

    !> C-type string as array
    character(kind=c_char), intent(in) :: cstring(*)

    !> Maximal string lenght. If C-string is longer, it will be chopped.
    integer, intent(in), optional  :: maxlen

    !> Resulting Fortran string
    character(:, kind=c_char), allocatable :: fortranChar

    integer :: ii, maxlen0

    if (present(maxlen)) then
      maxlen0 = maxlen
    else
      maxlen0 = huge(maxlen0) - 1
    end if

    do ii = 1, maxlen0
      if (cstring(ii) == c_null_char) then
        exit
      end if
    end do
    allocate(character(ii - 1) :: fortranChar)
    fortranChar = transfer(cstring(1 : ii - 1), fortranChar)

  end function fortranChar


end module dftbp_capi
