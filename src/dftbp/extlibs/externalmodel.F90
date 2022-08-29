!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Connects to an external energy/hamiltonian model or provides a dummy external model if no
!> external library is linked
module dftbp_externalmodel
#:if WITH_EXTERNALMODEL
  use iso_c_binding, only : c_int, c_char, c_bool, c_null_char, c_ptr, c_double, c_loc
#:endif
  use dftbp_io_clang, only : fortranChar
  use dftbp_io_message, only : error, warning
  use dftbp_common_accuracy, only : dp, mc
  use dftbp_common_modelTypes, only : modelTypes
  use dftbp_common_status, only : TStatus
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: TExtModelProvides, getExtModelCapabilities, externalModel_init, TExternalModel

  !> Return type for capabilities of external electronic structure/energy model
  type TExtModelProvides

    !> Name of external model
    character(:), allocatable :: modelName

    !> Is a hamiltonian provided
    logical :: gives_ham = .false.

    !> Is an overlap provdided, or should a diagonal matrix be substituted
    logical :: gives_overlap = .false.

    !> Are double-counting energies provided (i.e. non-bandstructure)
    logical :: gives_dc_energy = .false.

    !> Is the model self-consistently evaluated
    logical :: is_self_consistent = .false.

    !> Maximum number of spin channels supported by the model
    integer :: max_spin_channels = 0

    !> Are results for a subset of the atoms provided (allowing MPI parallism on the DFTB+ side)
    logical :: gives_atomic_subset = .false.

    !> Does the model require a communicator?
    logical :: is_mpi = .false.

    !> Internal model to cover parts not included in the external model (also dictates if total
    !> energies/forces/etc. can be calcualted)
    integer :: intModel = modelTypes%none

    !> Is an internal model used for repulsive/double counting?
    logical :: is_dc_internal = .false.

  end type TExtModelProvides


  !> External model settings and state
  type TExternalModel



  #:if WITH_EXTERNALMODEL

    !> Cut-off for hamiltonian interactions with atoms
    real(dp), allocatable :: cutoffs(:)

    !> C pointer to internal state of the model (assume that model manages this)
    type(c_ptr) :: modelState


  #:else

  #:endif

  contains

    !> Fills up the atomic orbital data structure
    procedure :: orbData

  #:if WITH_EXTERNALMODEL

    !> Clean up model, including invoking it's finalisation process
    final :: deleteCModel

  #:endif

  end type TExternalModel


  !> Current version for model input API and ABI requirements
  type :: TModelAPBI
    !> Semantic version of interface, both program and binary
    integer :: major, minor, patch
  end type TModelAPBI

  !> Current state of the interface standard, remember to update versionand document if changing
  !> this!
  type(TModelAPBI), parameter :: currentInterface = TModelAPBI(0,1,0)

#:if WITH_EXTERNALMODEL

  !> ctype for receiving external declared capabilities
  type, bind(c) :: extModelAbilities

    logical(c_bool) :: gives_ham
    logical(c_bool) :: gives_overlap
    logical(c_bool) :: gives_dc_energy
    logical(c_bool) :: requires_self_consistency
    integer(c_int) :: max_spin_channels
    logical(c_bool) :: gives_atom_subset
    logical(c_bool) :: requires_mpi

  end type extModelAbilities

  !> C code interface
  interface

    !> Check external model API version
    subroutine model_interface_apbi(major,minor,patch) bind(C, name='dftbp_model_apbi')
      import :: c_int
      implicit none
      integer(c_int), intent(out) :: major
      integer(c_int), intent(out) :: minor
      integer(c_int), intent(out) :: patch
    end subroutine model_interface_apbi

    !> External model declared capabilities
    subroutine model_provides(modelname, capabilities)&
        & bind(C, name='dftbp_provided_with')
      import :: extModelAbilities
      import :: c_char
      implicit none
      character(c_char), intent(out) :: modelname(*)
      type(extModelAbilities), intent(out) :: capabilities
    end subroutine model_provides

    !> Initialise external model for calculation
    function model_init(nspecies, speciesnames, cutoffs, modelstate, errString)&
        & bind(C, name='initialise_model_for_dftbp')
      import :: c_int, c_ptr, c_char, c_double
      implicit none
      integer(c_int), intent(in) :: nspecies
      type(c_ptr), target, intent(in) :: speciesnames(nspecies)
      real(c_double), intent(out) :: cutoffs(nspecies)
      type(c_ptr), target, intent(out) :: modelstate
      character(c_char), intent(out) :: errString(*)
      integer(c_int) :: model_init
    end function model_init

    !> Update external model for state of calculation
    function model_update(modelstate, errString) bind(C, name='update_model_for_dftbp')
      import :: c_ptr, c_int, c_char
      implicit none
      type(c_ptr), target, intent(in) :: modelstate
      character(c_char), intent(out) :: errString(*)
      integer(c_int) :: model_update
    end function model_update

    function model_destroy(modelstate, errString) bind(C, name='cleanup_model_for_dftbp')
      import :: c_ptr, c_int, c_char
      implicit none
      type(c_ptr), target, intent(inout) :: modelstate
      character(c_char), intent(out) :: errString(*)
      integer(c_int) :: model_destroy
    end function model_destroy

  end interface

#:endif

  !> Buffer size on the Fortran side
  integer, parameter :: nBufChar = 256

contains

  !> What are the capabilities of the attached external model
  subroutine getExtModelCapabilities(modelProvides, status)

    !> Status of operation
    type(TStatus), intent(out) :: status

    !> Capabilities of externally provided hamiltonian/energy model
    type(TExtModelProvides), intent(out) :: modelProvides

  #:if WITH_EXTERNALMODEL

    !> Structure for returned model capabilities
    type(extModelAbilities) :: capabilities

    !> Buffer on Fortran side, expecting a null termination somewhere inside of this, or throws an
    !> error
    character(nBufChar) :: modelname = " "

    integer :: major, minor, patch

    ! Test declared API/ABI of external model
    call model_interface_apbi(major, minor, patch)
    if (major /= currentInterface%major) then
      @:RAISE_FORMATTED_ERROR(status, -1, '(A,I0,A,I0,A)', 'External model has a breaking mismatch,&
          & version:',major,'. DFTB+ expects ',currentInterface%major,'.')
    end if
    if (minor < currentInterface%minor) then
      ! need to check for use of interface extensions after version x.0.z within the major release
      ! group
    end if
    if (patch /= currentInterface%patch) then
      ! Harmless
    end if

    call model_provides(modelname, capabilities)
    if (.not.isCStringOK(modelname, nBufChar)) then
      @:RAISE_ERROR(status, -1, "External model name string does not fit in buffer")
    end if
    modelProvides%modelName = trim(modelname)

    modelProvides%gives_ham = capabilities%gives_ham
    modelProvides%gives_overlap = capabilities%gives_overlap
    modelProvides%gives_dc_energy = capabilities%gives_dc_energy
    modelProvides%is_self_consistent = capabilities%requires_self_consistency
    modelProvides%max_spin_channels = capabilities%max_spin_channels
    modelProvides%gives_atomic_subset = capabilities%gives_atom_subset
    modelProvides%is_mpi = capabilities%requires_mpi

  #:else

    modelProvides%modelName = "Dummy Model"

    @:RAISE_ERROR(status, -1, "Dummy external model present, non-functioning calculation")

  #:endif

  end subroutine getExtModelCapabilities


  !> Initialise external model for current calculation
  subroutine externalModel_init(this, speciesNames, maxCutOff, orb, status)

    !> External model state and settings
    type(TExternalModel), intent(out) :: this

    !> labels of atomic species
    character(mc), intent(in) :: speciesNames(:)

    !> Cut-off for generating neighbour maps in the geometry
    real(dp), intent(out) :: maxCutOff

    !> Atomic orbital basis information
    type(TOrbitals), intent(out) :: orb

    !> Status of operation
    type(TStatus), intent(out) :: status

    integer :: iStatus, ii
    character(nBufChar) :: errString = " "
    integer :: nspecies
    character(mc+1), allocatable, target :: stringArray(:)

  #:if WITH_EXTERNALMODEL
    type(c_ptr), dimension(size(speciesnames)) :: speciesPtrs

    ! set up list of species labels to pass in
    nspecies = size(speciesNames)
    allocate(stringArray(nspecies))
    do ii = 1, nspecies
      stringArray(ii) = trim(speciesNames(ii)) // c_null_char
      speciesPtrs(ii) = c_loc(stringArray(ii))
    end do

    ! cutoff for hamiltonian interactions
    allocate(this%cutoffs(nspecies))

    iStatus = model_init(nspecies, speciesPtrs, this%cutoffs, this%modelState, errString)

    if (any(this%cutoffs < 0.0_dp)) then
      @:RAISE_ERROR(status, -2, "External model returned geometric cutoff less than zero")
    end if
    maxCutOff = maxval(this%cutoffs)

    if (iStatus /= 0) then
      if (.not.isCStringOK(errString, nBufChar)) then
        @:RAISE_ERROR(status, -1, "External model name string does not fit in buffer")
      end if
      @:RAISE_ERROR(status, iStatus, trim(errString))
    end if

    iStatus = model_update(this%modelState, errString)
    if (iStatus /= 0) then
      if (.not.isCStringOK(errString, nBufChar)) then
        @:RAISE_ERROR(status, -1, "External model name string does not fit in buffer")
      end if
      @:RAISE_ERROR(status, iStatus, trim(errString))
    end if

  #:else

    @:RAISE_ERROR(status, -1, "Dummy model, should have stopped before geting to here")

  #:endif

  end subroutine externalModel_init


#:if WITH_EXTERNALMODEL
  !> Checks if string has a null termination withing the expected length
  function isCStringOK(string, nBufferChar)

    !> String to check
    character(c_char), intent(in) :: string(*)

    !> Expected max length of string
    integer, intent(in) :: nBufferChar

    !> Is the string within the length and null terminated
    logical isCStringOK

    integer :: ii

    isCStringOK = .false.
    lpBufCheck: do ii = 1, nBufferChar
      if (string(ii) == c_null_char) then
        isCStringOK = .true.
        exit lpBufCheck
      end if
    end do lpBufCheck

  end function isCStringOK
#:endif

  subroutine orbData(this, orb, species, status)

    !> Instance
    class(TExternalModel), intent(in) :: this

    !> Atomic orbital information
    type(TOrbitals), intent(out) :: orb

    !> Species list for each atom
    integer, intent(in) :: species(:)

    !> Status of operation
    type(TStatus), intent(out) :: status

   #:if WITH_EXTERNALMODEL

   #:else
    @:RAISE_ERROR(status, -1, "Dummy model, should have stopped before geting to here")
   #:endif

  end subroutine orbData


#:if WITH_EXTERNALMODEL

  !> Ask C model to clean up its data structure
  subroutine deleteCModel(this)

    !> Instance
    type(TExternalModel), intent(inout) :: this

    integer :: iStatus, ii
    character(nBufChar) :: errString = " "

    call warning("Cleaning up external model")

    iStatus = model_destroy(this%modelState, errString)
    if (iStatus /= 0) then
      ! take things down completely, as there isn't an error return above this
      call error(errString)
    end if

  end subroutine deleteCModel

#:endif

end module dftbp_externalmodel
