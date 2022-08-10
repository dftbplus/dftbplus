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
module dftbp_extlibs_externalmodel
  use dftbp_common_clang, only : c_free, c_strlen
  use dftbp_common_accuracy, only : dp, mc
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_common_status, only : TStatus
  use dftbp_common_globalenv, only : stdOut
  use dftbp_dftb_nonscc, only : buildOrthogonalS
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_modelindependent_localclusters, only : TLocalClusters
  use iso_c_binding, only : c_int, c_char, c_bool, c_null_char, c_ptr, c_double, c_intptr_t, c_loc,&
      & c_f_pointer
  implicit none

  private
  public :: TExtModelProvides, TExtModel, getExtModelCapabilities, externalModel_init

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

    !> Order of derivatives returned by the model
    integer :: model_derivative_order = 0

    !> Is the model self-consistently evaluated
    logical :: is_self_consistent = .false.

    !> Number of spin channels (0 for spin free)
    integer :: nSpin = 0

    !> Internal model to cover parts not included in the external model (also dictates if total
    !> energies/forces/etc. can be calcualted)
    integer :: intModel = hamiltonianTypes%none

    !> Is an internal model used for repulsive/double counting?
    logical :: is_dc_internal = .false.

  end type TExtModelProvides


  !> Type for instance of model
  type TExtModel

    !> External model's capabilities
    type(TExtModelProvides) :: capabilities

    !> Distance cutoff for interactions and local environment
    real(dp) :: interactionCutoff

    !> Distance cutoff for surrounding environment
    real(dp) :: environmentCutoff

    !> Distance atom atom sites such that environmentCutoff can be checked around bonds
    real(dp) :: maxCutoff

    !> Number of shells for each atom
    integer, allocatable :: nShellsOnSpecies(:)

    !> Angular momentum of each atomic shell
    integer, allocatable :: shellLValues(:,:)

    !> Reference (neutral) occupations for atomic shells
    real(dp), allocatable :: shellOccs(:,:)

    !> C pointer to internal state of the model (assume that the model internally manages whatever
    !> is on the end of this pointer)
    integer(c_intptr_t) :: modelState

    !> Index in H/overlap structure for start of each atomic onsite block of integrals
    integer(c_int), allocatable :: atomClusterIndex(:)

    !> Index in H/overlap structure for start of each bond block of integrals
    integer(c_int), allocatable :: bondClusterIndex(:)

  contains

    !> Update external model for geometric changes
    procedure :: update

    !> Generate overlap/hamiltonian
    procedure :: buildSH0

    !> remove any external model
    final :: cleanup_model

  end type TExtModel

  !> ctype for receiving external capabilities
  type, bind(C) :: extModelAbilities

    logical(c_bool) :: gives_ham
    logical(c_bool) :: gives_overlap
    logical(c_bool) :: gives_dc_energy
    integer(c_int) :: model_derivative_order
    logical(c_bool) :: requires_self_consistency
    integer(c_int) :: spinchannels

  end type extModelAbilities

  !> C code interfaces
  interface

    !> External model declared API/ABI level (for any checks)
    subroutine APBI(major, minor, patch) bind(C, name='dftbp_model_apbi')
      import :: c_int
      implicit none
      !> Major version of release (potentially breaking changes)
      integer(c_int), intent(out) :: major
      !> Minor version (revised on non-breaking extensions)
      integer(c_int), intent(out) :: minor
      !> Patch level (transparent changes)
      integer(c_int), intent(out) :: patch
    end subroutine APBI

    !> External model declared capabilities
    subroutine model_provides(modelname, capabilities) bind(C, name='dftbp_provided_with')
      import :: extModelAbilities
      import :: c_char
      implicit none
      character(c_char), intent(out) :: modelname
      type(extModelAbilities), intent(out) :: capabilities
    end subroutine model_provides

    !> Initialise external model for calculation
    function model_init(nspecies, speciesnames, interactionCutoff, environmentCutoff,&
        & nShellsOnSpecies, shellLValues, shellOccs, modelstate, msgString)&
        & bind(C, name='initialise_model_for_dftbp')
      import :: c_int, c_ptr, c_char, c_double, c_intptr_t
      implicit none
      !> Number of distinct atomic species
      integer(c_int), intent(in) :: nspecies
      !> Labels for the atomic species
      type(c_ptr), target, intent(in) :: speciesnames(nspecies)
      !> Distance over which atoms interact
      real(c_double), intent(out) :: interactionCutoff
      !> Surrounding environment which affects an interaction
      real(c_double), intent(out) :: environmentCutoff
      !> Number of atomic shells on each species
      type(c_ptr), target, intent(out) :: nShellsOnSpecies
      !> L values for atomic shells
      type(c_ptr), target, intent(out) :: shellLValues
      !> Reference occupations for atomic shells in the neutral ground state
      type(c_ptr), target, intent(out) :: shellOccs
      !> Internal state of the model, including any stored data passed in from DFTB+
      integer(c_intptr_t) :: modelstate
      !> Error return string
      type(c_ptr), intent(out) :: msgString
      !> Is the model initialised?
      integer(c_int) :: model_init
    end function model_init

    !> Update external model for calculation state (geometric changes)
    function model_update(modelstate, speciesOfAtoms, nAtomClusters, indexAtomClusters,&
        & atomClusters, atClustersGlobalNos, nBondClusters, indexBondClusters, bondClusters,&
        & bndClustersGlobalNos, atomClusterIndex, bondClusterIndex, msgString)&
        & bind(C, name='update_model_for_dftbp')
      import :: c_int, c_ptr, c_char, c_intptr_t, c_double
      implicit none
      !> Internal state of model, opaque to us
      integer(c_intptr_t) :: modelstate
      !> Species of atoms in original structure
      integer(c_int) :: speciesOfAtoms
      !> Number of atomic clusters
      integer(c_int) :: nAtomClusters
      !> Index array for atomic clusters
      integer(c_int) :: indexAtomClusters
      !> Atomic clusters
      real(c_double) :: atomClusters
      !> Global structure numbers of atoms in atomic clusters
      integer(c_int) :: atClustersGlobalNos
      !> Number of bond clusters
      integer(c_int) :: nBondClusters
      !> Index array for bond clusters
      integer(c_int) :: indexBondClusters
      !> Indexing for sparse matrices holding overlap/hamiltonian onsite elements
      integer(c_int) :: atomClusterIndex
      !> Bond clusters
      real(c_double) :: bondClusters
      !> Global structure numbers of atoms in bond clusters
      integer(c_int) :: bndClustersGlobalNos
      !> Indexing for sparse matrices holding overlap/hamiltonian dimer elements
      integer(c_int) :: bondClusterIndex
      !> Any returned error string
      type(c_ptr), intent(out) :: msgString
      !> Model state after operation
      integer(c_int) :: model_update
    end function model_update

    !> Retrieve predictions from external model
    function model_predictions(modelstate, h0, over, iSparseStart, nMaxStart, msgString)&
        & bind(C, name='predict_model_for_dftbp')
      import :: c_int, c_double, c_ptr, c_intptr_t
      implicit none
      !> Internal state of model, opaque to us
      integer(c_intptr_t) :: modelstate
      !> Hamiltonian
      real(c_double), intent(inout) :: h0
      !> Overlap (if accessed)
      real(c_double), intent(inout) :: over
      !> Indexing for sparse structures
      integer(c_int), intent(in) :: iSparseStart
      !> Indexing for sparse structures
      integer(c_int), intent(in) :: nMaxStart
      !> Any returned error string
      type(c_ptr), intent(out) :: msgString
      !> Model state after operation
      integer(c_int) :: model_predictions
    end function model_predictions

    !> Clean up the C external model
    subroutine cleanup(modelstate) bind(C, name='cleanup_model_for_dftbp')
      import :: c_int, c_double, c_ptr, c_intptr_t
      implicit none
      !> Internal state of model, opaque to us
      integer(c_intptr_t) :: modelstate
    end subroutine cleanup

  end interface

  !> Buffer size on the Fortran side
  integer, parameter :: nBufChar = 256

contains

  !> What are the capabilities of the attached external model
  subroutine getExtModelCapabilities(modelProvides, status)

    !> Capabilities of externally provided hamiltonian/energy model
    type(TExtModelProvides), intent(out) :: modelProvides

    !> Status of operation
    type(TStatus), intent(out) :: status

    !> Structure for returned model capabilities
    type(extModelAbilities) :: capabilities

    !> Buffer on Fortran side, expecting a null termination somewhere inside of this, or throws an
    !> error
    character(nBufChar) :: modelname = " "

    integer :: major, minor, patch

    call apbi(major, minor, patch)

  #:block DEBUG_CODE
    write(stdOut,'(1X,A,I0,A,I0,A,I0)')'External model API/ABI : ',major,'.',minor,'.',patch
  #:endblock DEBUG_CODE
    if (.not.(major == 0 .and. minor == 1)) then
      @:RAISE_ERROR(status, -1, "External model API/ABI non compliant (expecting 0.1.X)")
    end if

    call model_provides(modelname, capabilities)
    if (.not.isCStringOK(modelname, nBufChar)) then
      @:RAISE_ERROR(status, -1, "External model name string does not fit in buffer")
    end if
    modelProvides%modelName = trim(modelname)
    write(stdOut,'(1X,A,A,A)')'External model used : "', modelProvides%modelName,'"'
    modelProvides%gives_ham = capabilities%gives_ham
    modelProvides%gives_overlap = capabilities%gives_overlap
    modelProvides%gives_dc_energy = capabilities%gives_dc_energy
    modelProvides%model_derivative_order = capabilities%model_derivative_order
    modelProvides%is_self_consistent = capabilities%requires_self_consistency
    modelProvides%nSpin = capabilities%spinchannels

  end subroutine getExtModelCapabilities


  !> Initialise external model for current calculation
  subroutine externalModel_init(this, capabilities, speciesNames, status)

    !> Instance of external model
    type(TExtModel), allocatable, intent(out) :: this

    !> Structure for returned model capabilities
    type(TExtModelProvides), intent(in) :: capabilities

    !> labels of atomic species
    character(mc), intent(in) :: speciesNames(:)

    !> Status of operation
    type(TStatus), intent(out) :: status

    integer :: iStatus, ii, iSp, iSh
    type(c_ptr) :: msgString
    integer(c_int) :: nspecies
    type(c_ptr), dimension(size(speciesnames)) :: speciesPtrs
    character(mc+1), allocatable, target :: stringArray(:)
    real(c_double) :: interactionCutoff
    real(c_double) :: environmentCutoff
    type(c_ptr) :: cptr_nshells, cptr_shells, cptr_shellOccs
    integer, pointer :: fptr_nShells(:), fptr_shells(:)
    real(dp), pointer :: fptr_shellOccs(:)

    allocate(this)

    this%capabilities = capabilities
    @:PROPAGATE_ERROR(status)

    nspecies = size(speciesNames)
    allocate(stringArray(nspecies))

    do ii = 1, nspecies
      stringArray(ii) = trim(speciesNames(ii)) // c_null_char
      speciesPtrs(ii) = c_loc(stringArray(ii))
    end do

    iStatus = model_init(nspecies, speciesPtrs, interactionCutoff, environmentCutoff, cptr_nShells,&
        & cptr_shells, cptr_shellOccs, this%modelState, msgString)

    if (iStatus /= 0) then
      @:RAISE_ERROR(status, iStatus, "To do : from initialise_model_for_dftbp")
    end if

    this%interactionCutoff = interactionCutoff
    ! scale to allow for a cylinder around the bond between interacting atoms
    this%environmentCutoff = environmentCutoff
    this%maxCutoff = sqrt(environmentCutoff**2 + (interactionCutoff**2)/4.0_dp)

    allocate(this%nShellsOnSpecies(nspecies))
    call c_f_pointer(cptr_nShells, fptr_nShells, [nSpecies])
    this%nShellsOnSpecies(:) = fptr_nShells

    call c_f_pointer(cptr_shells, fptr_shells, [sum(this%nShellsOnSpecies)])
    allocate(this%shellLValues(maxval(this%nShellsOnSpecies), nspecies))
    this%shellLValues(:,:) = 0
    ii = 1
    do iSp = 1, nSpecies
      do iSh = 1, this%nShellsOnSpecies(iSp)
        this%shellLValues(iSh,iSp) = fptr_shells(ii)
        ii = ii + 1
      end do
    end do

    call c_f_pointer(cptr_shellOccs, fptr_shellOccs, [sum(this%nShellsOnSpecies)])
    allocate(this%shellOccs(maxval(this%nShellsOnSpecies), nspecies))
    this%shellOccs(:,:) = 0.0_dp
    ii = 1
    do iSp = 1, nSpecies
      do iSh = 1, this%nShellsOnSpecies(iSp)
        this%shellOccs(iSh,iSp) = fptr_shellOccs(ii)
        ii = ii + 1
      end do
    end do

    ! clean up
    fptr_nShells => null()
    call c_free(cptr_nShells)
    fptr_shells => null()
    call c_free(cptr_shells)
    fptr_shellOccs => null()
    call c_free(cptr_shellOccs)

  end subroutine externalModel_init


  !> Updates instance of external model for change in the geometry, and hence local environments
  !> around atoms/bonds
  subroutine update(this, species, clusters, iSparseStart, nNeighbourSK, status)

    !> Instance of external model
    class(TExtModel), intent(inout) :: this

    !> Chemical species of atoms in global structure
    integer, intent(in) :: species(:)

    !> Local geometric environments around atoms and bonds
    type(TLocalClusters), allocatable :: clusters

    !> Index vector, where the interaction between two atoms starts in the sparse format.
    integer, intent(in) :: iSparseStart(0:,:)

    !> Number of neighbours of each real atom
    integer, intent(in) :: nNeighbourSK(:)

    !> Status of operation
    type(TStatus), intent(out) :: status

    integer :: iStatus
    type(c_ptr) :: msgString
    integer :: iClust, iAt, iNeigh, nAtoms

    if (allocated(clusters)) then

      @:ASSERT(allocated(clusters%atomClusters).eqv.allocated(clusters%iStartAtCluster))
      @:ASSERT(allocated(clusters%atomClusters))
      write(stdOut,*)'Clusters set up'

      nAtoms = size(clusters%iStartAtCluster)-1
      if (allocated(this%atomClusterIndex)) then
        if (size(this%atomClusterIndex) < nAtoms) then
          deallocate(this%atomClusterIndex)
        end if
      end if
      if (.not. allocated(this%atomClusterIndex)) then
        allocate(this%atomClusterIndex(nAtoms), source=0)
      end if
      do iAt = 1, nAtoms
        this%atomClusterIndex(iAt) = iSparseStart(0, iAt)
      end do

      if (allocated(this%bondClusterIndex)) then
        deallocate(this%bondClusterIndex)
      end if
      allocate(this%bondClusterIndex(clusters%nBondClusters), source=0)
      iClust = 1
      do iAt = 1, size(clusters%iStartAtCluster)-1
        lpNeigh: do iNeigh = 1, nNeighbourSK(iAt)
          this%bondClusterIndex(iClust) = iSparseStart(iNeigh, iAt)
          iClust = iClust + 1
        end do lpNeigh
      end do

      iStatus = model_update(this%modelState, species(1), size(clusters%iStartAtCluster)-1,&
          & clusters%iStartAtCluster(1), clusters%atomClusters(1,1), clusters%iAtCluster(1),&
          & clusters%nBondClusters, clusters%iStartBondCluster(1), clusters%bondClusters(1,1),&
          & clusters%iBondCluster(1), this%atomClusterIndex(1), this%bondClusterIndex(1), msgString)

      if (iStatus /= 0) then
        @:RAISE_ERROR(status, iStatus, "To do : from update_model_for_dftbp")
      end if

    else

      @:RAISE_ERROR(status, iStatus, "Internal error, local cluster geometries not allocated")

    end if

  end subroutine update


  !> Generates electronic structure part of total model
  subroutine buildSH0(this, env, H0, over, nAtom, species, iSparseStart, orb, status)

    !> Instance of external model
    class(TExtModel), intent(in) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Non-self-consistent hamiltonian matrix
    real(dp), target, intent(inout) :: H0(:)

    !> Overlap matrix (needed if referenced by the external model)
    real(dp), intent(out) :: over(:)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Index vector, where the interaction between two atoms starts in the sparse format.
    integer, intent(in) :: iSparseStart(:,:)

    !> Information about the orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Status of operation
    type(TStatus), intent(out) :: status

    integer :: iStatus
    type(c_ptr) :: msgString

    H0(:) = 0.0_dp
    over(:) = 0.0_dp

    iStatus = model_predictions(this%modelState, h0(1), over(1), iSparseStart(1,1),&
        & size(iSparseStart,dim=1), msgString)

    if (iStatus /= 0) then
      @:RAISE_ERROR(status, iStatus, "To do : from predict_model_for_dftbp")
    end if

    if (.not.this%capabilities%gives_overlap) then
      call buildOrthogonalS(env, over, nAtom, species, iSparseStart, orb)
    end if

  end subroutine buildSH0


  !> Checks if string has a null termination within the expected length
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


  !> Cleans up the external model when this type goes out of scope
  subroutine cleanup_model(this)

    !> Instance
    type(TExtModel), intent(inout) :: this

    call cleanup(this%modelState)

  end subroutine cleanup_model

end module dftbp_extlibs_externalmodel
