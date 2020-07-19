!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

! Tests the subroutine 'setSpeciesAndDependents' which can be used to reset quantities that depend
! on the global ordering of atoms species, i.e.
!
! species(1:nAtoms) = (/1,2,2,1,...,1/) -> (/1,2,2,1,...,2/),
!
! following an MD step. A situation that arises in molecular dynamics packages that use domain
! decomposition as a parallelisation scheme.
!
! In this tester the lead process reads in geometry and broadcasts to all other processes. Geometry
! is not distributed upon entry to or inside DFTB+ (note that it will be checked though if compiled
! in Debug mode). nAtoms conserved for a given hsd_tree, to reset, call destruct and re-initialise.
!
! Forces are not used to update atomic positions, rather two sets of atomic positions are read from
! file and forces are computed for them
!
! NOTE: count_substrings function requires no leading white space for species names in structure.gen

#! Fypp preprocessing required as this can be compiled either as serial or MPI parallel
#:set WITH_MPI = defined('WITH_MPI')

program test_setSpeciesAndDependents
  use, intrinsic :: iso_fortran_env, only: output_unit, REAL64, IOSTAT_END
#:if WITH_MPI
  use mpi
#:endif
  use dftbp_mmapi, only: TDftbPlus_init, TDftbPlus_destruct, TDftbPlus, TDftbPlusInput
  use dftbp_hsdapi, only: fnode, getChild, getChildren, setChild, getChildValue, setChildValue
  use dftbp_hsdapi, only: dumpHsd
  use testhelpers, only: writeAutotestTag
  implicit none

  !> Precision and unit conversion
  integer, parameter :: dp = REAL64
  real(dp), parameter :: Bohr__AA = 0.529177249_dp
  real(dp), parameter :: AA__Bohr = 1.0_dp / Bohr__AA

  !> DFTB Objects
  type(TDftbPlus) :: dftb
  type(TDftbPlusInput) :: hsd_tree
  character(*), parameter :: dftb_fname = 'dftb_in.hsd'

  !> Type for containing geometrical information
  !> One can write your own type, rather than copy DTFB+
  type TGeometry
    integer :: nAtom
    logical :: tPeriodic
    logical :: tFracCoord
    integer, allocatable :: species(:)
    real(dp), allocatable :: coords(:,:)
    integer :: nSpecies
    real(dp), allocatable :: origin(:)
    real(dp), allocatable :: latVecs(:,:)
    real(dp), allocatable :: recVecs2p(:,:)
    character(50), allocatable :: speciesNames(:)
  end type TGeometry

  !> MD indexing type
  type MDstatus_type
    integer :: initial_step
    integer :: final_step
  end type MDstatus_type

#:if WITH_MPI
  !> MPI variables
  integer, parameter :: requiredThreading = MPI_THREAD_FUNNELED
  integer :: providedThreading, rank, np, ierr
  integer, parameter :: lead_id = 0
  logical :: IO
#:else
  !> Dummy variables for serial code
  integer, parameter :: MPI_COMM_WORLD = 0
  integer, parameter :: lead_id = 0
  logical :: IO = .true.
#:endif

  ! Local
  integer, parameter :: nSteps = 2
  integer :: nAtoms

  type(MDstatus_type):: MDstatus
  type(TGeometry) :: geo
  integer :: imd
  real(dp) :: merminEnergy
  character(100) :: imd_lab
  character(100) :: fname
  real(dp), allocatable :: gradients(:,:), stressTensor(:,:), grossCharges(:)

#:if WITH_MPI
  ! Initialise MPI environment
  call mpi_init_thread(requiredThreading, providedThreading, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, np, ierr)
  IO = (rank == lead_id)
#:endif

  nAtoms = get_number_of_atoms('structure_1.gen')
  MDstatus = MDstatus_type(1, nSteps) ! initialise type
  allocate(gradients(3,nAtoms))
  allocate(grossCharges(nAtoms))

  do imd = MDstatus%initial_step, MDstatus%final_step

    if(IO) then
      write(*,*) 'MD step ', imd, 'of ', MDstatus%final_step
    end if

    ! For every step, read in atomic positions from file
    ! Equivalent to obtaining data externally, e.g. from an MD package
    write(imd_lab, '(I0)') imd
    fname = 'structure_'//trim(imd_lab)//'.gen'
    call read_in_geo(trim(fname), geo)
  #:if WITH_MPI
    call broadcast_geometry(MPI_COMM_WORLD, geo)
  #:endif

    if (geo%nAtom /= nAtoms) then
      write(*,*) 'Error: Number of atoms not conserved between MD steps'
    #:if WITH_MPI
      call mpi_abort(MPI_COMM_WORLD, 0, ierr)
    #:else
      stop
    #:endif
    endif

    if(imd == MDstatus%initial_step) then
    #:if WITH_MPI
      call TDftbPlus_init(dftb, output_unit, MPI_COMM_WORLD)
    #:else
      call TDftbPlus_init(dftb, output_unit)
    #:endif

      call initialise_dftbplus_tree(geo, dftb, hsd_tree)

      ! Dump hsd tree to file "fort.1" if debugging
      ! call dumpHsd(hsd_tree%hsdTree, 001)
      call dftb%setupCalculator(hsd_tree)
    endif

    ! Update coordinates and lattice vectors
    call dftb%setGeometry(geo%coords, geo%latVecs)

    ! Update species order every step
    call dftb%setSpeciesAndDependents(geo%speciesNames, geo%species)

    ! Do a total energy calculation
    call dftb%getEnergy(merminEnergy)
    call dftb%getGradients(gradients)
    call dftb%getGrossCharges(grossCharges)

    if (geo%tPeriodic) then
      if (.not.allocated(stressTensor)) then
        allocate(stressTensor(3,3))
      end if
      call dftb%getStressTensor(stressTensor)
    else
      if (allocated(stressTensor)) then
        deallocate(stressTensor)
      end if
    end if

    ! call output_forces_per_process(gradients, imd_lab)

  enddo

  ! Clean up
  call TDftbPlus_destruct(dftb)

  ! Write file for internal test system, using the last structure that was run
  call writeAutotestTag(merminEnergy=merminEnergy, gradients=gradients, stressTensor=stressTensor,&
      & grossCharges=grossCharges)

#:if WITH_MPI
  call mpi_finalize(ierr)
#:endif

contains

  !--------------------------------------------
  ! Utility routines based on DL_POLY V4.10 API
  !--------------------------------------------

  !> Initialise DFTB tree
  subroutine initialise_dftbplus_tree(geo, dftb, hsd_tree, input_fname)
    implicit none

    !> Geometry for DFTB+
    type(tGeometry), intent( In ) :: geo

    !> DFTB+ calculator
    type(TDftbPlus), intent( InOut ) :: dftb

    !> DFTB+ Input tree
    type(TDftbPlusInput), intent( Out ) :: hsd_tree

    !> Optional DFTB+ input file name
    character(*), intent( In ), optional :: input_fname

    ! Local data

    ! DFTB+ input file name
    character(100) :: fname

    ! Pointers to the parts of the input tree that will be set
    type(fnode), pointer :: pRoot, pGeo, pAnalysis

    ! "Does geometry already exist in DTFB+ input?" (== "replace geometry in HSD tree?")
    logical :: replace_geometry

    if(present(input_fname)) Then
      fname = input_fname
    else
      fname = dftb_fname
    endif

    ! Read dftb input file into hsd_tree
    call dftb%getInputFromFile(fname, hsd_tree)

    ! Set structure data, retain rest
    call hsd_tree%getRootNode(pRoot)

    call setChild(pRoot, "Geometry", pGeo, replace=replace_geometry)

    call setChildValue(pGeo, "Periodic", geo%tPeriodic)

    if (geo%tPeriodic) then
      call setChildValue(pGeo, "LatticeVectors", geo%latVecs)
    end if

    call setChildValue(pGeo, "TypeNames", geo%speciesNames)

    ! See DFTB+ subroutine in lib_type/typegeometryhsd.F90
    call setChildValue(pGeo, "TypesAndCoordinates",&
        & reshape(geo%species, (/ 1, size(geo%species) /)), geo%coords)

    ! Always compute forces
    call setChild(pRoot, "Analysis", pAnalysis, replace=.True.)
    call setChildValue(pAnalysis, "CalculateForces", .True.)

  end subroutine initialise_dftbplus_tree

#:if WITH_MPI

  !> Output forces per process
  subroutine output_forces_per_process(gradients, imd_lab)
    implicit none

    !> atomic gradients = -forces
    real(dp), allocatable, intent(in) :: gradients(:,:)

    !> Label string for file
    character(*), intent(in) :: imd_lab

    integer :: io_unit, ia
    character(100) :: rank_lab

    write(rank_lab, '(I0)') rank

    open(newunit=io_unit, file= 'forces_MD' // trim(imd_lab) // '_rank' // trim(rank_lab) // '.dat')
    write(io_unit,*) '# One line of header'
    do ia = 1, size(gradients,2)
      write(io_unit, "(I5, 3F20.12)") ia, -1.0_dp * gradients(:,ia)
    enddo

    close(io_unit)

  end subroutine output_forces_per_process

#:endif


  !> Count number of substrings on a line that are separated by whitespace
  function count_substrings(input_line) result(n_substrings)
    implicit none

    !> string to search
    character(*), intent(in) :: input_line

    !> Resulting count
    integer :: n_substrings

    character(len(input_line)) :: line
    character(1), parameter :: whitespace = ' '
    integer :: ii, length

    ! Add 1 prefixed and 1 trailing white space so counter works correctly
    line = ' ' // input_line
    length = len(trim(adjustl(line))) + 1

    n_substrings = 0
    do ii = 2, length
      if( line(ii:ii) /= whitespace .and. line(ii-1:ii-1) == whitespace)then
        n_substrings = n_substrings + 1
      endif
    enddo

  end function count_substrings


  !> Get number of atoms from DFTB+ geometry file
  function get_number_of_atoms(fname, headerLines) result(nAtoms)
    implicit none

    !> File name to read
    character(*), intent(in) :: fname

    !> Number of header lines to discard on reading
    integer, optional, intent(in) :: headerLines

    !> Resulting number of atoms
    integer :: nAtoms

    ! Boundary condition label
    character(1) :: boundary

    integer :: io_unit, ii

    if (IO) then
      open(newunit=io_unit,file=fname)

      ! Skip (additional) header lines at top of file
      if(present(headerLines))then
        do ii = 1, headerLines
          read(io_unit,*)
        enddo
      endif

      read(io_unit,*) nAtoms, boundary

      close(io_unit)
    endif

#:if WITH_MPI
    call mpi_bcast(nAtoms, 1, MPI_INTEGER, lead_id, MPI_COMM_WORLD, ierr)
#:endif

  end function get_number_of_atoms


  !> Read DFTB+ geometry file (on lead_id processor and broadcast, if MPI parallel)
  subroutine read_in_geo(fname, geo, headerLines)
    implicit none

    !> Name of file to read
    character(*), intent(in) :: fname

    !> Resulting geometry object
    type(TGeometry), intent(out) :: geo

    !> Lines of header material to discard before reading
    integer, optional, intent(in) :: headerLines

    integer :: ia, ii, nSpecies
    character(1) :: boundary

    ! note, this assumes 2 character chemical symbols (DFTB+ itself does not require this limit)
    character(2), allocatable :: species_list(:)

    character(100) :: line2
    integer :: io_unit

    if (IO) then

      !Initialise geo components
      geo%nAtom = 0
      geo%nSpecies = 0

      ! Open file
      open(newunit=io_unit, file=fname)

      ! Skip (optional) header lines at top of file
      if(present(headerLines))then
        do ii = 1, headerLines
          read(io_unit,*)
        enddo
      endif

      ! Assign geometry data from file in gen format

      ! Line 1
      read(io_unit,*) geo%nAtom, boundary

      if(boundary=='S')then
        geo%tPeriodic = .true.
        geo%tFracCoord = .false.
      elseif(boundary=='F')then
        geo%tPeriodic = .true.
        geo%tFracCoord = .true.
      elseif(boundary=='C')then
        geo%tPeriodic = .false.
        geo%tFracCoord = .false.
      endif

      ! Line 2
      read(io_unit, '(A)') line2
      nSpecies = count_substrings(line2)
      backspace io_unit

      ! assumes max 2 character atomic species names
      allocate(species_list(nSpecies))
      species_list = ' '
      read(io_unit,*) species_list(:)

      geo%nSpecies = nSpecies
      if(.not. allocated(geo%speciesNames)) then
        allocate(geo%speciesNames(geo%nSpecies))
      end if
      geo%speciesNames = ' '
      do ia=1,geo%nSpecies
        geo%speciesNames(ia) = species_list(ia)
      enddo

      ! Lines 3 -> 3 + Natoms
      if(.not. allocated(geo%species)) then
        allocate(geo%species(geo%nAtom))
      end if
      if(.not. allocated(geo%coords)) then
        allocate(geo%coords(3, geo%nAtom))
      end if
      geo%species(:) = 0
      geo%coords(:,:) = 0
      do ia = 1, geo%nAtom
        read(io_unit,*) ii, geo%species(ia), geo%coords(1:3,ia)
      enddo
      geo%coords(:,:) = geo%coords * AA__Bohr

      if(boundary /= 'C')then
        ! Origin, which we discard, instead of passing onto DFTB+
        read(io_unit,*)

        ! Lattice vectors in rows in .gen BUT as columns internally
        if(.not. allocated(geo%latVecs)) allocate(geo%latVecs(3,3))
        do ia=1,3
          read(io_unit,*) geo%latVecs(:,ia)
        enddo
        geo%latVecs(:,:) = geo%latVecs * AA__Bohr
      endif

      close(io_unit)
    endif

  end subroutine read_in_geo


#:if WITH_MPI
  !> Broadcast geometry from lead_id to all processes
  subroutine broadcast_geometry(comm, geo)

    !> MPI communicator
    integer, intent(in) :: comm

    !> Geometry: Filled on lead_id, uninitialised on all other processes
    type(TGeometry), intent(inout) :: geo

    integer :: ierr, lenSN

    call mpi_bcast(geo%nAtom, 1, MPI_INTEGER, lead_id, comm, ierr)
    call mpi_bcast(geo%tPeriodic, 1, MPI_LOGICAL, lead_id, comm, ierr)
    call mpi_bcast(geo%tFracCoord, 1, MPI_LOGICAL, lead_id, comm, ierr)
    call mpi_bcast(geo%nSpecies, 1, MPI_INTEGER, lead_id, comm, ierr)

    if(.not. allocated(geo%speciesNames)) then
      allocate(geo%speciesNames(geo%nSpecies))
    endif

    lenSN = len(geo%speciesNames)
    call mpi_bcast(geo%speciesNames, lenSN * geo%nSpecies, MPI_CHARACTER, lead_id, comm, ierr)

    if(.not. allocated(geo%species)) then
      allocate(geo%species(geo%nAtom))
    endif

    if(.not. allocated(geo%coords)) then
      allocate(geo%coords(3, geo%nAtom))
    endif

    call mpi_bcast(geo%species, geo%nAtom, MPI_INTEGER, lead_id, comm, ierr)
    call mpi_bcast(geo%coords, 3 * geo%nAtom, MPI_DOUBLE_PRECISION, lead_id, comm, ierr)

    if (geo%tPeriodic) then
      if(.not. allocated(geo%latVecs)) then
        allocate(geo%latVecs(3,3))
      endif
      call mpi_bcast(geo%latVecs, 9, MPI_DOUBLE_PRECISION, lead_id, comm, ierr)
    endif

    call mpi_barrier(comm, ierr)

  end subroutine broadcast_geometry

#:endif

end program test_setSpeciesAndDependents
