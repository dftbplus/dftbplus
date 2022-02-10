!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains the C-API of DFTB+.
module dftbp_capi
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env
  use dftbp_common_accuracy, only : dp
  use dftbp_common_globalenv, only : instanceSafeBuild
  use dftbp_dftbplus_qdepextpotgenc, only :&
      & getExtPotIfaceC, getExtPotGradIfaceC, TQDepExtPotGenC, TQDepExtPotGenC_init
  use dftbp_mmapi, only :&
      & TDftbPlus, TDftbPlus_init, TDftbPlus_destruct, TDftbPlusInput, TDftbPlusAtomList
  use dftbp_type_linkedlist, only : TListString, append, init, destruct
  implicit none
  private

  !> DFTB+ atom list structure
  type, bind(C) :: c_DftbPlusAtomList
    type(c_ptr) :: pDftbPlusAtomList
  end type c_DftbPlusAtomList

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


  !> Queries, whether library is instance safe (allowing for concurrent DFTB+ instances)
  function c_DftbPlus_is_instance_safe() result(instSafe) bind(C, name='dftbp_is_instance_safe')

    !> Whether library was built ensuring instance safety (multiple instances supported)
    logical(c_bool) :: instSafe

    instSafe = instanceSafeBuild

  end function c_DftbPlus_is_instance_safe


  !> Initialises a DFTB+ calculation with output sent to to some location
  subroutine c_DftbPlus_init(handler, outputFileName) bind(C, name='dftbp_init')

    !> DFTB+ handler
    type(c_DftbPlus), intent(out) :: handler

    !> output location
    type(c_ptr), value, intent(in) :: outputFileName

    type(TDftbPlusC), pointer :: instance

    allocate(instance)
    call handleOutputFileName(outputFileName, instance%outputUnit, instance%tOutputOpened)
    call TDftbPlus_init(instance%TDftbPlus, outputUnit=instance%outputUnit)
    handler%instance = c_loc(instance)

  end subroutine c_DftbPlus_init


  !> Initialises a DFTB+ calculator (MPI-version)
  subroutine c_DftbPlus_init_mpi(handler, outputFileName, mpiComm) bind(C, name='dftbp_init_mpi')

    !> DFTB+ handler
    type(c_DftbPlus), intent(out) :: handler

    !> output location
    type(c_ptr), value, intent(in) :: outputFileName

    !> MPI-communicator id
    integer(c_int), value, intent(in) :: mpiComm

    type(TDftbPlusC), pointer :: instance

    allocate(instance)
    call handleOutputFileName(outputFileName, instance%outputUnit, instance%tOutputOpened)
    call TDftbPlus_init(instance%TDftbPlus, outputUnit=instance%outputUnit, mpiComm=mpiComm)
    handler%instance = c_loc(instance)

  end subroutine c_DftbPlus_init_mpi


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
    type(c_DftbPlusInput), intent(inout) :: inputHandler

    type(TDftbPlusC), pointer :: instance
    type(TDftbPlusInput), pointer :: pDftbPlusInput
    character(:), allocatable :: fortranFileName

    allocate(pDftbPlusInput)
    fortranFileName = fortranChar(fileName)
    call c_f_pointer(handler%instance, instance)
    call instance%getInputFromFile(fortranFileName, pDftbPlusInput)
    inputHandler%pDftbPlusInput = c_loc(pDftbPlusInput)

  end subroutine c_DftbPlus_getInputFromFile


  !> Process a document tree to get settings for the calculation
  subroutine c_DftbPlus_processInput(handler, inputHandler) bind(C, name='dftbp_process_input')

    !> handler for the calculation instance
    type(c_DftbPlus), intent(inout) :: handler

    !> input tree handler
    type(c_DftbPlusInput), intent(inout) :: inputHandler

    type(TDftbPlusC), pointer :: instance
    type(TDftbPlusInput), pointer :: pDftbPlusInput
    type(TDftbPlusAtomList), pointer :: pDftbPlusAtomList

    call c_f_pointer(handler%instance, instance)
    call c_f_pointer(inputHandler%pDftbPlusInput, pDftbPlusInput)
    call instance%setupCalculator(pDftbPlusInput)

  end subroutine c_DftbPlus_processInput


  !> Get electrostatic potential from DFTB+ at specified points
  subroutine c_DftbPlus_get_elstat_potential(handler, nLocations, pot, locations)&
      & bind(C, name='dftbp_get_elstat_potential')

    !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Number of requested points
    integer, value, intent(in) :: nLocations

    !> Resulting potentials
    real(dp), intent(out) :: pot(*)

    !> Sites to calculate potential
    real(c_double), intent(in) :: locations(3,*)

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)

    call instance%getElStatPotential(pot(1:nLocations), locations(:, 1:nLocations))

  end subroutine c_DftbPlus_get_elstat_potential


  !> Set an external potential on the DFTB+ calculation
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


  !> Register a generator for an external potential
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


  !> register DM exporting callback
  subroutine c_DftbPlus_registerDMCallback(handler, callback, aux_ptr)&
      & bind(C, name='dftbp_register_dm_callback')

    !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> callback function for DM export
    type(c_funptr), value :: callback

    !> pointer to a context object for the callback
    type(c_ptr), value :: aux_ptr

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    call instance%registerDMCallback(callback, aux_ptr)

  end subroutine c_DftbPlus_registerDMCallback

  !> register S exporting callback
  subroutine c_DftbPlus_registerSCallback(handler, callback, aux_ptr)&
      & bind(C, name='dftbp_register_s_callback')

    !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> callback function for DM export
    type(c_funptr), value :: callback

    !> pointer to a context object for the callback
    type(c_ptr), value :: aux_ptr

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    call instance%registerSCallback(callback, aux_ptr)

  end subroutine c_DftbPlus_registerSCallback

  !> register H exporting callback
  subroutine c_DftbPlus_registerHCallback(handler, callback, aux_ptr)&
      & bind(C, name='dftbp_register_h_callback')

    !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> callback function for DM export
    type(c_funptr), value :: callback

    !> pointer to a context object for the callback
    type(c_ptr), value :: aux_ptr

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    call instance%registerHCallback(callback, aux_ptr)

  end subroutine c_DftbPlus_registerHCallback

  !> Set/replace the coordinates in a DFTB+ calculation instance
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

  !> Obtain nr. of spin channels.
  function c_DftbPlus_nrOfSpin(handler) result(nSpin) bind(C, name='dftbp_get_nr_spin')
    type(c_DftbPlus), intent(inout) :: handler
    integer(c_int) :: nSpin

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    nSpin = instance%nrOfSpin()
  end function c_DftbPlus_nrOfSpin

  !> Obtain nr. of k-points.
  function c_DftbPlus_nrOfKPoints(handler) result(nKpoints) bind(C, name='dftbp_get_nr_kpoints')
    type(c_DftbPlus), intent(inout) :: handler
    integer(c_int) :: nKpoints

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    nKpoints = instance%nrOfKPoints()
  end function c_DftbPlus_nrOfKPoints
  
  !> Obtain nr. of (k-point,spin chanel) pairs in current process group.
  function c_DftbPlus_nrOfLocalKS(handler) result(nLocalKS) bind(C, name='dftbp_get_nr_local_ks')
    type(c_DftbPlus), intent(inout) :: handler
    integer(c_int) :: nLocalKS

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    nLocalKS = instance%nrOfLocalKS()

  end function c_DftbPlus_nrOfLocalKS
  
  !> Get (k-point,spin chanel) pairs in current process group, returns number of pairs
  function c_DftbPlus_getLocalKS(handler, localKS) result(nLocalKS) bind(C, name='dftbp_get_local_ks')
  !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler
    integer(c_int), intent(out) :: localKS(2, *)

    type(TDftbPlusC), pointer :: instance
    integer :: nLocalKS

    call c_f_pointer(handler%instance, instance)
    nLocalKS = instance%nrOfLocalKS()

    call instance%getLocalKS(localKS(:,1:nLocalKS))

  end function c_DftbPlus_getLocalKS


  !> Obtain size of the basis set (NxN)
  function c_DftbPlus_getBasisSize(handler) result(BasisSize) bind(C, name='dftbp_get_basis_size')
  !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    type(TDftbPlusC), pointer :: instance
    integer :: BasisSize

    call c_f_pointer(handler%instance, instance)
    
    BasisSize = instance%getBasisSize()

  end function c_DftbPlus_getBasisSize
  

  !> TODO
  function c_DftbPlus_isHSReal(handler) result(HSReal) bind(C, name='dftbp_is_hs_real')
  !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    type(TDftbPlusC), pointer :: instance
    logical :: HSReal

    call c_f_pointer(handler%instance, instance)
    
    HSReal = instance%isHSReal()

  end function c_DftbPlus_isHSReal
  
  
  !> TODO
  function c_DftbPlus_getOverlap(handler, blacs_descr) result(data_ptr) bind(C, name='dftbp_get_overlap')
  !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler
    type(c_ptr), value :: blacs_descr

    type(c_ptr) :: data_ptr
    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    
    data_ptr = instance%getOverlap(blacs_descr)

  end function c_DftbPlus_getOverlap
  
  
  !> TODO
  function c_DftbPlus_getHamiltonian(handler, blacs_descr) result(data_ptr) bind(C, name='dftbp_get_hamiltonian')
  !> handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler
    type(c_ptr), value :: blacs_descr

    type(c_ptr) :: data_ptr
    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    
    data_ptr = instance%getHamiltonian(blacs_descr)

  end function c_DftbPlus_getHamiltonian


  !> Obtain nr. of k-points.
  function c_DftbPlus_nrOfKPoints(handler) result(nKPoints) bind(C, name='dftbp_nr_kpoints')
    type(c_DftbPlus), intent(inout) :: handler
    integer(c_int) :: nKPoints

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    nKPoints = instance%nrOfKPoints()

  end function c_DftbPlus_nrOfKPoints


  !> Obtain masses of atoms in the DFTB+ calculation.
  subroutine c_DftbPlus_get_atom_masses(handler, masses) bind(C, name='dftbp_get_masses')
    type(c_DftbPlus), intent(inout) :: handler
    real(c_double), intent(out) :: masses(*)

    type(TDftbPlusC), pointer :: instance
    integer :: nAtom

    call c_f_pointer(handler%instance, instance)
    nAtom = instance%nrOfAtoms()
    call instance%getAtomicMasses(masses(:nAtom))

  end subroutine c_DftbPlus_get_atom_masses


  !> Obtain nr. basis functions at each atoms in the DFTB+ calculation.
  subroutine c_DftbPlus_get_atom_nr_basis(handler, nOrbitals) bind(C, name='dftbp_get_nr_orbitals')
    type(c_DftbPlus), intent(inout) :: handler
    integer(c_int), intent(out) :: nOrbitals(*)

    type(TDftbPlusC), pointer :: instance
    integer :: nAtom

    call c_f_pointer(handler%instance, instance)
    nAtom = instance%nrOfAtoms()
    call instance%getNOrbitalsOnAtoms(nOrbitals(:nAtom))

  end subroutine c_DftbPlus_get_atom_nr_basis


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

    !> Maximal string length. If C-string is longer, it will be chopped.
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


  !> Handles the optional output file name (which should be a NULL-ptr if not present)
  subroutine handleOutputFileName(outputFileName, outputUnit, tOutputOpened)
    type(c_ptr), intent(in) :: outputFileName
    integer, intent(out) :: outputUnit
    logical, intent(out) :: tOutputOpened

    character(c_char), pointer :: pOutputFileName
    character(:), allocatable :: fortranFileName

    if (c_associated(outputFileName)) then
      call c_f_pointer(outputFileName, pOutputFileName)
      fortranFileName = fortranChar(pOutputFileName)
      open(newunit=outputUnit, file=fortranFileName, action="write")
      tOutputOpened = .true.
    else
      outputUnit = output_unit
      tOutputOpened = .false.
    end if

  end subroutine handleOutputFileName


end module dftbp_capi
