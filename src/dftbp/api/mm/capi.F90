!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains the C-API of DFTB+.
module dftbp_capi
  use, intrinsic :: iso_c_binding, only : c_associated, c_bool, c_char, c_double, c_f_pointer,&
      & c_f_procpointer, c_funptr, c_int, c_loc, c_null_char, c_null_ptr, c_ptr
  use, intrinsic :: iso_fortran_env, only : output_unit
  use dftbp_capicallback, only: dmhs_callback_c_wrapper, set_dmhs_callback_c_wrapper, TCAuxWrapper
  use dftbp_common_accuracy, only : dp
  use dftbp_common_file, only : openFile, TFileDescr
  use dftbp_common_globalenv, only : instanceSafeBuild
  use dftbp_dftbplus_qdepextpotgenc, only : getExtPotGradIfaceC, getExtPotIfaceC, TQDepExtPotGenC,&
      & TQDepExtPotGenC_init
  use dftbp_mmapi, only : TDftbPlus, TDftbPlus_destruct, TDftbPlus_init, TDftbPlusAtomList,&
      & TDftbPlusInput
  implicit none
  private


  !> DFTB+ input tree
  type, bind(C) :: c_DftbPlusInput
    type(c_ptr) :: pDftbPlusInput = c_null_ptr
  end type c_DftbPlusInput


  !> DFTB+ calculation
  type, bind(C) :: c_DftbPlus
    type(c_ptr) :: instance = c_null_ptr
  end type c_DftbPlus


  !> Simple extension around the TDftbPlus with some additional variables for the C-API.
  type, extends(TDftbPlus) :: TDftbPlusC
    private
    type(TFileDescr) :: outputFile
  contains
    final :: TDftbPlusC_final
  end type TDftbPlusC


contains

  !> Finalises a DFTB+ input instance
  subroutine c_DftbPlusInput_final(handler) bind(C, name='dftbp_input_final')

    !> DFTB+ handler
    type(c_DftbPlusInput), intent(inout) :: handler

    !> The specific instance to be finalised
    type(TDftbPlusInput), pointer :: instance

    if (.not. c_associated(handler%pDftbPlusInput)) return
    call c_f_pointer(handler%pDftbPlusInput, instance)
    deallocate(instance)
    handler%pDftbPlusInput = c_null_ptr

  end subroutine c_DftbPlusInput_final


  !> Returns the current API version
  subroutine c_DftbPlus_api(major, minor, patch) bind(C, name='dftbp_api')

    !> major.minor.patch
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

    !> Output location
    type(c_ptr), value, intent(in) :: outputFileName

    type(TDftbPlusC), pointer :: instance

    allocate(instance)
    call handleOutputFile_(outputFileName, instance%outputFile)
    call TDftbPlus_init(instance%TDftbPlus, outputUnit=instance%outputFile%unit)
    handler%instance = c_loc(instance)

  end subroutine c_DftbPlus_init


  !> Initialises a DFTB+ calculator (MPI-version)
  subroutine c_DftbPlus_init_mpi(handler, outputFileName, mpiComm) bind(C, name='dftbp_init_mpi')

    !> DFTB+ handler
    type(c_DftbPlus), intent(out) :: handler

    !> Output location
    type(c_ptr), value, intent(in) :: outputFileName

    !> MPI-communicator id
    integer(c_int), value, intent(in) :: mpiComm

    type(TDftbPlusC), pointer :: instance

    allocate(instance)
    call handleOutputFile_(outputFileName, instance%outputFile)
    call TDftbPlus_init(instance%TDftbPlus, outputUnit=instance%outputFile%unit, mpiComm=mpiComm)
    handler%instance = c_loc(instance)

  end subroutine c_DftbPlus_init_mpi


  !> Finalises a DFTB+ instance
  subroutine c_DftbPlus_final(handler) bind(C, name='dftbp_final')

    !> DFTB+ handler
    type(c_DftbPlus), intent(inout) :: handler

    !> The specific instance to be finalised
    type(TDftbPlusC), pointer :: instance

    if (.not. c_associated(handler%instance)) return
    call c_f_pointer(handler%instance, instance)
    deallocate(instance)
    handler%instance = c_null_ptr

  end subroutine c_DftbPlus_final


  !> Read input for DFTB+ from a specified file
  subroutine c_DftbPlus_getInputFromFile(handler, fileName, inputHandler)&
      & bind(C, name='dftbp_get_input_from_file')

    !> Handler for the input
    type(c_DftbPlus), intent(inout) :: handler

    !> File to read
    character(c_char), intent(in) :: fileName(*)

    !> Handler for the resulting input
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

    !> Handler for the calculation instance
    type(c_DftbPlus), intent(inout) :: handler

    !> Input tree handler
    type(c_DftbPlusInput), intent(inout) :: inputHandler

    type(TDftbPlusC), pointer :: instance
    type(TDftbPlusInput), pointer :: pDftbPlusInput

    call c_f_pointer(handler%instance, instance)
    call c_f_pointer(inputHandler%pDftbPlusInput, pDftbPlusInput)
    call instance%setupCalculator(pDftbPlusInput)

  end subroutine c_DftbPlus_processInput


  !> Get electrostatic potential from DFTB+ at specified points
  subroutine c_DftbPlus_get_elstat_potential(handler, nLocations, pot, locations)&
      & bind(C, name='dftbp_get_elstat_potential')

    !> Handler for the calculation
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

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Externally set potential
    real(c_double), intent(in) :: extPot(*)

    !> Gradient of the potential wrt to atom positions
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

    !> Handler for the potential
    type(c_DftbPlus), intent(inout) :: handler

    !> Pointer to the C routine for the external potential
    type(c_ptr), value, intent(in) :: refPtr

    !> Function for the external potential
    type(c_funptr), value, intent(in) :: extPotFunc

    !> Function for the gradient of the potential
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


  !> Register a density matrix exporting callback
  subroutine c_DftbPlus_registerDMCallback(handler, callback, aux_ptr)&
      & bind(C, name='dftbp_register_dm_callback')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Callback function for DM export
    type(c_funptr), value :: callback

    !> Pointer to a context object for the callback
    type(c_ptr), value :: aux_ptr

    type(TDftbPlusC), pointer :: instance
    class(*), pointer :: wrapper

    call c_f_pointer(handler%instance, instance)

    allocate(TCAuxWrapper :: wrapper)
    select type(wrapper)
    type is (TCAuxWrapper)
      wrapper%auxPtr = aux_ptr
      wrapper%callback = callback
    end select
    call instance%registerDMCallback(dmhs_callback_c_wrapper, wrapper)

  end subroutine c_DftbPlus_registerDMCallback


  !> Register overlap matrix exporting callback
  subroutine c_DftbPlus_registerSCallback(handler, callback, aux_ptr)&
      & bind(C, name='dftbp_register_s_callback')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Callback function for DM export
    type(c_funptr), value :: callback

    !> Pointer to a context object for the callback
    type(c_ptr), value :: aux_ptr

    type(TDftbPlusC), pointer :: instance
    class(*), pointer :: wrapper

    call c_f_pointer(handler%instance, instance)

    allocate(TCAuxWrapper :: wrapper)
    select type(wrapper)
    type is (TCAuxWrapper)
      wrapper%auxPtr = aux_ptr
      wrapper%callback = callback
    end select

    call instance%registerSCallback(dmhs_callback_c_wrapper, wrapper)

  end subroutine c_DftbPlus_registerSCallback


  !> Register overlap matrix importing callback
  subroutine c_DftbPlus_registerSetSCallback(handler, callback, aux_ptr)&
      & bind(C, name='dftbp_register_set_s_callback')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Callback function for DM export
    type(c_funptr), value :: callback

    !> Pointer to a context object for the callback
    type(c_ptr), value :: aux_ptr

    type(TDftbPlusC), pointer :: instance
    class(*), pointer :: wrapper

    call c_f_pointer(handler%instance, instance)

    allocate(TCAuxWrapper :: wrapper)
    select type(wrapper)
    type is (TCAuxWrapper)
      wrapper%auxPtr = aux_ptr
      wrapper%callback = callback
    end select

    call instance%registerSetSCallback(set_dmhs_callback_c_wrapper, wrapper)

  end subroutine c_DftbPlus_registerSetSCallback


  !> Register hamiltonian exporting callback
  subroutine c_DftbPlus_registerHCallback(handler, callback, aux_ptr)&
      & bind(C, name='dftbp_register_h_callback')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Callback function for DM export
    type(c_funptr), value :: callback

    !> Pointer to a context object for the callback
    type(c_ptr), value :: aux_ptr

    type(TDftbPlusC), pointer :: instance
    class(*), pointer :: wrapper

    call c_f_pointer(handler%instance, instance)

    allocate(TCAuxWrapper :: wrapper)
    select type(wrapper)
    type is (TCAuxWrapper)
      wrapper%auxPtr = aux_ptr
      wrapper%callback = callback
    end select
    call instance%registerHCallback(dmhs_callback_c_wrapper, wrapper)

  end subroutine c_DftbPlus_registerHCallback


  !> Register hamiltonian importing callback
  subroutine c_DftbPlus_registerSetHCallback(handler, callback, aux_ptr)&
      & bind(C, name='dftbp_register_set_h_callback')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Callback function for DM export
    type(c_funptr), value :: callback

    !> Pointer to a context object for the callback
    type(c_ptr), value :: aux_ptr

    type(TDftbPlusC), pointer :: instance
    class(*), pointer :: wrapper

    call c_f_pointer(handler%instance, instance)

    allocate(TCAuxWrapper :: wrapper)
    select type(wrapper)
    type is (TCAuxWrapper)
      wrapper%auxPtr = aux_ptr
      wrapper%callback = callback
    end select
    call instance%registerSetHCallback(set_dmhs_callback_c_wrapper, wrapper)

  end subroutine c_DftbPlus_registerSetHCallback


  !> Set/replace the coordinates in a DFTB+ calculation instance
  subroutine c_DftbPlus_setCoords(handler, coords) bind(C, name='dftbp_set_coords')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Coordinates, (xyz, :nAtom)
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

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Coordinates, row major format (xyz, :nAtom)
    real(c_double), intent(in) :: coords(3,*)

    !> Lattice vectors, row major format
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

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Coordinates, row major format (xyz, :nAtom)
    real(c_double), intent(in) :: coords(3,*)

    !> Lattice vectors, row major format
    real(c_double), intent(in) :: latvecs(3, *)

    !> Coordinate origin
    real(c_double), intent(in) :: origin(3)

    type(TDftbPlusC), pointer :: instance
    integer :: nAtom

    call c_f_pointer(handler%instance, instance)
    nAtom = instance%nrOfAtoms()
    call instance%setGeometry(coords(:, 1:nAtom), latVecs(:, 1:3), origin(1:3))

  end subroutine c_DftbPlus_setCoordsLatticeVecsOrigin


  !> Set the neighbour list instead of computing it in DFTB+
  subroutine c_DftbPlus_setNeighbourList(handler, nAllAtom, nMaxNeighbours, nNeighbour,&
      & iNeighbour, neighDist, cutOff, coordNeighbours, neighbour2CentCell)&
      & bind(C, name='dftbp_set_neighbour_list')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Total number of neighbour atoms
    integer(c_int), value, intent(in) :: nAllAtom

    !> Maximum number of neighbours an atom in the central cell can have
    integer(c_int), value, intent(in) :: nMaxNeighbours

    !> Number of neighbours for each atom in the central cell
    integer(c_int), intent(in) :: nNeighbour(*)

    !> References to the neighbour atoms for an atom in the central cell
    integer(c_int), intent(in) :: iNeighbour(nMaxNeighbours, *)

    !> Distances to the neighbour atoms for an atom in the central cell
    real(c_double), intent(in) :: neighDist(nMaxNeighbours, *)

    !> Cutoff distance used for this neighbour list
    real(c_double), value, intent(in) :: cutOff

    !> Coordinates of all neighbours
    real(c_double), intent(in) :: coordNeighbours(3, *)

    !> Mapping between neighbour reference and atom index in the central cell
    integer(c_int), intent(in) :: neighbour2CentCell(*)

    type(TDftbPlusC), pointer :: instance
    integer :: nAtom

    call c_f_pointer(handler%instance, instance)
    nAtom = instance%nrOfAtoms()
    call instance%setNeighbourList(nNeighbour(1:nAtom), iNeighbour(:,1:nAtom),&
        & neighDist(:,1:nAtom), cutOff, coordNeighbours(:,1:nAllAtom),&
        & neighbour2CentCell(1:nAllAtom))

  end subroutine c_DftbPlus_setNeighbourList


  !> Obtain nr. of atoms.
  function c_DftbPlus_nrOfAtoms(handler) result(nAtom) bind(C, name='dftbp_get_nr_atoms')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Number of atoms in the system
    integer(c_int) :: nAtom

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    nAtom = instance%nrOfAtoms()

  end function c_DftbPlus_nrOfAtoms


  !> Obtain nr. of spin channels.
  function c_DftbPlus_nrOfSpin(handler) result(nSpin) bind(C, name='dftbp_get_nr_spin')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Number of spin channels in the calculation
    integer(c_int) :: nSpin

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    nSpin = instance%nrOfSpin()

  end function c_DftbPlus_nrOfSpin


  !> Obtain nr. of (k-point,spin chanel) pairs in current process group.
  function c_DftbPlus_nrOfLocalKS(handler) result(nLocalKS) bind(C, name='dftbp_get_nr_local_ks')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Total number of k and spin local to this processor group
    integer(c_int) :: nLocalKS

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    nLocalKS = instance%nrOfLocalKS()

  end function c_DftbPlus_nrOfLocalKS


  !> Get (k-point,spin chanel) pairs in current process group, returns number of pairs
  function c_DftbPlus_getLocalKS(handler, localKS) result(nLocalKS)&
      & bind(C, name='dftbp_get_local_ks')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> K and spin local to this processor group
    integer(c_int), intent(out) :: localKS(2, *)

    !> Number of local K&S values on this COMM
    integer(c_int) :: nLocalKS

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    nLocalKS = instance%nrOfLocalKS()

    call instance%getLocalKS(localKS(:,1:nLocalKS))

  end function c_DftbPlus_getLocalKS


  !> Queries weights of k-points
  subroutine c_DftbPlus_getKWeights(handler, kweights)  bind(C, name='dftbp_get_kweights')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> K-point weights for all k-points
    real(c_double), intent(out) :: kweights(*)

    type(TDftbPlusC), pointer :: instance
    integer :: nkpts

    call c_f_pointer(handler%instance, instance)
    nkpts = instance%nrOfKPoints()

    call instance%getKWeights(kweights(1:nkpts))

  end subroutine c_DftbPlus_getKWeights


  !> Obtain total size of the basis set
  function c_DftbPlus_getBasisSize(handler) result(BasisSize) bind(C, name='dftbp_get_basis_size')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Total number of spatial basis functions
    integer(c_int) :: BasisSize

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)

    BasisSize = instance%getBasisSize()

  end function c_DftbPlus_getBasisSize


  !> Whether the system is described with real matrices
  function c_DftbPlus_isHSReal(handler) result(HSReal) bind(C, name='dftbp_is_hs_real')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Is this a real matrix calculation (typically a molecule without spin-orbit)
    logical(kind=C_BOOL) :: HSReal

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)

    HSReal = instance%isHSReal()

  end function c_DftbPlus_isHSReal


  !> Obtain total nr. of k-points.
  function c_DftbPlus_nrOfKPoints(handler) result(nKPoints) bind(C, name='dftbp_nr_kpoints')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Total number of k-points
    integer(c_int) :: nKPoints

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    nKPoints = instance%nrOfKPoints()

  end function c_DftbPlus_nrOfKPoints


  !> Obtain masses of atoms in the DFTB+ calculation.
  subroutine c_DftbPlus_get_atom_masses(handler, masses) bind(C, name='dftbp_get_masses')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Masses of atoms
    real(c_double), intent(out) :: masses(*)

    type(TDftbPlusC), pointer :: instance
    integer :: nAtom

    call c_f_pointer(handler%instance, instance)
    nAtom = instance%nrOfAtoms()
    call instance%getAtomicMasses(masses(:nAtom))

  end subroutine c_DftbPlus_get_atom_masses


  !> Obtain nr. basis functions at each atoms in the DFTB+ calculation.
  subroutine c_DftbPlus_get_atom_nr_basis(handler, nOrbitals) bind(C, name='dftbp_get_nr_orbitals')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Number of atomic orbitals (basis fns.) for each atom
    integer(c_int), intent(out) :: nOrbitals(*)

    type(TDftbPlusC), pointer :: instance
    integer :: nAtom

    call c_f_pointer(handler%instance, instance)
    nAtom = instance%nrOfAtoms()
    call instance%getNOrbitalsOnAtoms(nOrbitals(:nAtom))

  end subroutine c_DftbPlus_get_atom_nr_basis


  !> Retrieve the maximum cutoff distance that is being used for all interactions
  function c_DftbPlus_getCutOff(handler) result(cutOff) bind(C, name='dftbp_get_cutoff')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Cutoff distance
    real(dp) :: cutOff

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    cutOff = instance%getCutOff()

  end function c_DftbPlus_getCutOff


  !> Obtain the DFTB+ energy
  subroutine c_DftbPlus_getEnergy(handler, merminEnergy) bind(C, name='dftbp_get_energy')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Resulting energy
    real(c_double), intent(out) :: merminEnergy

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    call instance%getEnergy(merminEnergy)

  end subroutine c_DftbPlus_getEnergy


  !> Obtain the gradients wrt DFTB atom positions
  subroutine c_DftbPlus_getGradients(handler, gradients) bind(C, name='dftbp_get_gradients')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Gradients, row major format
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

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Gradients, row major format
    real(c_double), intent(out) :: stresstensor(3, 3)

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)

    call instance%getStressTensor(stresstensor(:, 1:3))

  end subroutine c_DftbPlus_getStressTensor


  !> Obtain gross (Mulliken) charges for atoms wrt to neutral references
  subroutine c_DftbPlus_getGrossCharges(handler, atomCharges)&
      & bind(C, name='dftbp_get_gross_charges')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Resulting atomic charges
    real(c_double), intent(out) :: atomCharges(*)

    type(TDftbPlusC), pointer :: instance
    integer :: nAtom

    call c_f_pointer(handler%instance, instance)
    nAtom = instance%nrOfAtoms()
    call instance%getGrossCharges(atomCharges(1:nAtom))

  end subroutine c_DftbPlus_getGrossCharges


  !> Obtain CM5 charges
  subroutine c_DftbPlus_getCM5Charges(handler, atomCharges)&
      & bind(C, name='dftbp_get_cm5_charges')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Resulting atomic charges
    real(c_double), intent(out) :: atomCharges(*)

    !> F pointer of input arguments
    type(TDftbPlusC), pointer :: instance

    integer :: nAtom

    !> Translate c to f objects
    call c_f_pointer(handler%instance, instance)

    nAtom = instance%nrOfAtoms()
    call instance%getCM5Charges(atomCharges(1:nAtom))

  end subroutine c_DftbPlus_getCM5Charges


  !> Obtain reference atomic charge, the effective Z for the valence orbitals
  subroutine c_DftbPlus_getRefCharges(handler, refCharges)&
      & bind(C, name='dftbp_get_ref_charges')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Resulting atomic reference charges
    real(c_double), intent(out) :: refCharges(*)

    !> F pointer of input arguments
    type(TDftbPlusC), pointer :: instance

    integer :: nAtom

    ! translate c to f objects
    call c_f_pointer(handler%instance, instance)

    nAtom = instance%nrOfAtoms()
    call instance%getRefCharges(refCharges(1:nAtom))

  end subroutine c_DftbPlus_getRefCharges


  !> Set reference atomic charge, the effective Z for the valence orbitals
  subroutine c_DftbPlus_setRefCharge(handler, refCharges)&
      & bind(C, name='dftbp_set_ref_charges')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Provided atomic reference charges
    real(c_double), intent(in) :: refCharges(*)

    !> F pointer of input arguments
    type(TDftbPlusC), pointer :: instance

    integer :: nAtom

    ! translate c to f objects
    call c_f_pointer(handler%instance, instance)

    nAtom = instance%nrOfAtoms()
    call instance%setRefCharges(refCharges(1:nAtom))

  end subroutine c_DftbPlus_setRefCharge


  !> Finalizer for TDftbPlusC
  subroutine TDftbPlusC_final(this)

    !> Instance
    type(TDftbPlusC), intent(inout) :: this

    ! Note: Fortran finalizes all components of a child class instance (TDftbPlusC) first and only
    ! then the components of its parent (TDftbPlus). TDftbPlusC contains a descriptor connected to
    ! an open file, whose unit had been passed to and stored by TDftbPlus. When TDftbPlusC is
    ! finalized, the file is closed, so TDftbPlus will try to write the timings to an invalid unit
    ! when finalized aftewards. Therefore, we call TDftbPlus_destruct explicitely before
    ! finalization of TDftbPlusC happens.
    !
    call TDftbPlus_destruct(this%TDftbPlus)

  end subroutine TDftbPlusC_final


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


  !> Returns a unit for an opened output file.
  !!
  !! If outputFileName is associated, a file with that name will be created (and returned),
  !! otherwise the output_unit is returned (and no file is created)
  !!
  subroutine handleOutputFile_(outputFileName, outputFile)

    !> File name to open
    type(c_ptr), intent(in) :: outputFileName

    !> Resulting descriptor for the file handle
    type(TFileDescr), intent(out) :: outputFile

    character(c_char), pointer :: pOutputFileName
    character(:), allocatable :: fortranFileName

    if (c_associated(outputFileName)) then
      call c_f_pointer(outputFileName, pOutputFileName)
      fortranFileName = fortranChar(pOutputFileName)
      call openFile(outputFile, fortranFileName, mode="w")
    else
      call outputFile%connectToUnit(output_unit)
    end if

  end subroutine handleOutputFile_


end module dftbp_capi
