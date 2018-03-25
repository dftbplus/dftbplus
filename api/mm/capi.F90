!> Contains the C-API to DFTB+
module dftbp_capi
  use, intrinsic :: iso_c_binding
  use dftbp_accuracy, only : dp
  use dftbp_mmapi
  implicit none
  private

  type, bind(C) :: TDftbPlusInputC
    type(c_ptr) :: pDftbPlusInput
  end type TDftbPlusInputC

  type, bind(C) :: TDftbPlusC
    type(c_ptr) :: pDftbPlus
  end type TDftbPlusC
  

contains

  subroutine TDftbPlusC_init(wrapper) bind(C, name='dftbp_init')
    type(TDftbPlusC), intent(out) :: wrapper

    type(TDftbPlus), pointer :: pDftbPlus

    allocate(pDftbPlus)
    call TDftbPlus_init(pDftbPlus)
    wrapper%pDftbPlus = c_loc(pDftbPlus)

  end subroutine TDftbPlusC_init

  
  subroutine TDftbPlusC_destruct(wrapper) bind(C, name='dftbp_destruct')
    type(TDftbPlusC), intent(inout) :: wrapper

    type(TDftbPlus), pointer :: pDftbPlus

    call c_f_pointer(wrapper%pDftbPlus, pDftbPlus)
    call TDftbPlus_destruct(pDftbPlus)
    deallocate(pDftbPlus)
    wrapper%pDftbPlus = c_null_ptr
    
  end subroutine TDftbPlusC_destruct


  subroutine TDftbPlusC_getInputFromFile(wrapper, fileName, inputWrapper)&
      & bind(C, name='dftbp_get_input_from_file')
    type(TDftbPlusC), intent(inout) :: wrapper
    character(c_char), intent(in) :: fileName(*)
    type(TDftbPlusInputC), intent(out) :: inputWrapper

    type(TDftbPlus), pointer :: pDftbPlus
    type(TDftbPlusInput), pointer :: pDftbPlusInput
    character(:), allocatable :: fortranFileName

    allocate(pDftbPlusInput)
    fortranFileName = fortranChar(fileName)
    call c_f_pointer(wrapper%pDftbPlus, pDftbPlus)
    call pDftbPlus%getInputFromFile(fortranFileName, pDftbPlusInput)
    inputWrapper%pDftbPlusInput = c_loc(pDftbPlusInput)
    
  end subroutine TDftbPlusC_getInputFromFile


  subroutine TDftbPlusC_setupCalculator(wrapper, inputWrapper)&
      & bind(C, name='dftbp_setup_calculator')
    type(TDftbPlusC), intent(inout) :: wrapper
    type(TDftbPlusInputC), intent(inout) :: inputWrapper

    type(TDftbPlus), pointer :: pDftbPlus
    type(TDftbPlusInput), pointer :: pDftbPlusInput

    call c_f_pointer(wrapper%pDftbPlus, pDftbPlus)
    call c_f_pointer(inputWrapper%pDftbPlusInput, pDftbPlusInput)
    call pDftbPlus%setupCalculator(pDftbPlusInput)
    
  end subroutine TDftbPlusC_setupCalculator


  subroutine TDftbPlusC_setCoords(wrapper, coords) bind(C, name='dftbp_set_coords')
    type(TDftbPlusC), intent(inout) :: wrapper
    real(c_double), intent(in) :: coords(3,*)

    type(TDftbPlus), pointer :: pDftbPlus
    integer :: nAtom

    call c_f_pointer(wrapper%pDftbPlus, pDftbPlus)
    nAtom = pDftbPlus%nrOfAtoms()
    call pDftbPlus%setGeometry(coords(:, 1:nAtom))

  end subroutine TDftbPlusC_setCoords


  subroutine TDftbPlusC_setCoordsAndLatVecs(wrapper, coords, latVecs)&
      & bind(C, name='dftbp_set_coords_and_latvecs')
    type(TDftbPlusC), intent(inout) :: wrapper
    real(c_double), intent(in) :: coords(3,*)
    real(c_double), intent(in) :: latvecs(3, *)

    type(TDftbPlus), pointer :: pDftbPlus
    integer :: nAtom

    call c_f_pointer(wrapper%pDftbPlus, pDftbPlus)
    nAtom = pDftbPlus%nrOfAtoms()
    call pDftbPlus%setGeometry(coords(:, 1:nAtom), latVecs(:, 1:3))

  end subroutine TDftbPlusC_setCoordsAndLatVecs


  subroutine TDftbPlusC_getEnergy(wrapper, merminEnergy) bind(C, name='dftbp_get_energy')
    type(TDftbPlusC), intent(inout) :: wrapper
    real(c_double), intent(out) :: merminEnergy

    type(TDftbPlus), pointer :: pDftbPlus

    call c_f_pointer(wrapper%pDftbPlus, pDftbPlus)
    call pDftbPlus%getEnergy(merminEnergy)

  end subroutine TDftbPlusC_getEnergy


  subroutine TDftbPlusC_getGradients(wrapper, gradients) bind(C, name='dftbp_get_gradients')
    type(TDftbPlusC), intent(inout) :: wrapper
    real(c_double), intent(out) :: gradients(3, *)

    type(TDftbPlus), pointer :: pDftbPlus
    integer :: nAtom

    call c_f_pointer(wrapper%pDftbPlus, pDftbPlus)
    nAtom = pDftbPlus%nrOfAtoms()
    call pDftbPlus%getGradients(gradients(:, 1:nAtom))

  end subroutine TDftbPlusC_getGradients


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
