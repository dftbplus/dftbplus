!---------------------------------------------------------------------
! Test 'setSpeciesAndDependents' when order of species changes
!--------------------------------------------------------------------
!
! Test subroutine 'setSpeciesAndDependents' which can
! be used to reset quantities that depend on the global ordering of
! atoms species. I.e.
!
!  species(1:nAtoms) = (/1,2,2,1,...,1/) -> (/1,2,2,1,...,2/),
!
! following an MD step. A situation that arises in molecular dynamics
! packages that use domain decomposition as a parallelisation scheme.
!
! Master process reads in geometry and broadcasts to all other processes.
! Geometry is not distributed upon entry to or in DFTB+.
! nAtoms conserved.
!
! Forces are not used to update atomic positions, rather two sets
! of atomic positions are read from file and forces are computed
! for them
!
! NOTE: count_substrings function requires no leading white space
! for species names in structure.gen
!
! Author: A. Buccheri. 2020. University of Bristol.
! alexanderbuccheri@googlemail.com
!---------------------------------------------------------------------

#! Fypp preprocessing from command line is required
#:set WITH_MPI = defined('WITH_MPI')

program test_setSpeciesAndDependents
  use, intrinsic :: iso_fortran_env, only: output_unit, REAL64, IOSTAT_END
#:if WITH_MPI
  use mpi
#:endif
  use dftbp_mmapi,  only: TDftbPlus_init, TDftbPlus_destruct, TDftbPlus, TDftbPlusInput
  use dftbp_hsdapi, only: fnode, getChild, getChildren, setChild, getChildValue, setChildValue
  use dftbp_hsdapi, only: dumpHsd
  use testhelpers,  only: writeAutotestTag
  implicit none

  !> Precision and unit conversion
  integer, parameter :: dp = REAL64
  real(dp), parameter :: Bohr__AA = 0.529177249_dp
  real(dp), parameter :: AA__Bohr = 1.0_dp / Bohr__AA

  !> DFTB Objects
  type(TDftbPlus)      :: dftb
  type(TDftbPlusInput) :: hsd_tree
  Character(len=11), parameter :: dftb_fname = 'dftb_in.hsd'

  !> Type for containing geometrical information
  !  One can write their own type, rather than copy DTFB+'s
  type TGeometry
    integer :: nAtom
    logical :: tPeriodic
    logical :: tFracCoord
    integer,  allocatable :: species(:)
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

  !> MPI
#:if WITH_MPI
  integer, parameter :: requiredThreading = MPI_THREAD_FUNNELED
  integer            :: providedThreading, rank, np, ierr
  integer, parameter :: master_id = 0
  logical            :: IO
#:else
  integer, parameter :: MPI_COMM_WORLD = 0
  integer, parameter :: master_id = 0
  logical            :: IO = .true.
#:endif
  
  !Local
  integer, parameter :: Nsteps = 2
  integer            :: nAtoms 
  
  type(MDstatus_type):: MDstatus
  type(TGeometry)    :: geo
  integer            :: imd
  real(dp)           :: merminEnergy
  character(len=2)   :: imd_lab
  character(len=100) :: fname
  real(dp), allocatable :: gradients(:,:), stressTensor(:,:), grossCharges(:)

  !Initialise MPI environment
#:if WITH_MPI
  call mpi_init_thread(requiredThreading, providedThreading, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, np, ierr)
  IO = (rank == master_id)
#:endif
  
  nAtoms = get_number_of_atoms('structure_1.gen')
  MDstatus = MDstatus_type(1, Nsteps)
  allocate(gradients(3,nAtoms))
  allocate(grossCharges(nAtoms))

  do imd = MDstatus%initial_step, MDstatus%final_step
     
    if(IO) then
      write(*,*) 'MD step ', imd, 'of ', MDstatus%final_step
    end if

    !For every step, read in atomic positions from file
    !Equivalent to obtaining data from MD package
    write(imd_lab, '(I2)') imd
    fname = 'structure_'//trim(adjustl(imd_lab))//'.gen'
    call read_in_geo(trim(adjustl(fname)), geo)
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
      Call TDftbPlus_init(dftb, output_unit, MPI_COMM_WORLD)
  #:else
      Call TDftbPlus_init(dftb, output_unit)
  #:endif
      Call initialise_dftbplus_tree(geo, dftb, hsd_tree)
      !Dump hsd tree to fort.1 if debugging 
      !Call dumpHsd(hsd_tree%hsdTree, 001)
      Call dftb%setupCalculator(hsd_tree)
    endif

    !Update coordinates and lattice vectors
    call dftb%setGeometry(geo%coords, geo%latVecs)

    !Update species order every step
    call dftb%setSpeciesAndDependents(geo%speciesNames, geo%species)

    !Do a total energy calculation
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

    !call output_forces_per_process(gradients, imd_lab)
    
  enddo

  !Clean up
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
  ! A. Buccheri 2020
  !--------------------------------------------

  !> @brief Initialise DFTB tree
  !!
  !! @param[in]   geo         Geometry for DFTB+
  !! @param[out]  hsd_tree    DFTB+ Input tree
  !! @param[in]   input_fname Optional DFTB+ input file name
  !
  Subroutine initialise_dftbplus_tree(geo, dftb, hsd_tree, input_fname)
    !Arguments
    Type(tGeometry), Intent( In    ) :: geo
    Type(TDftbPlus),          Intent( InOut ) :: dftb
    Type(TDftbPlusInput),     Intent(   Out ) :: hsd_tree
    Character(Len=*),         Intent( In    ),   Optional :: input_fname

    !Local data
    !> DFTB+ input file name
    Character(Len=100) :: fname
    !> Pointers to the parts of the input tree that will be set
    Type(fnode), Pointer :: pRoot, pGeo, pAnalysis
    !> "Does geometry already exist in DTFB+ input?" (== "replace geometry in HSD tree?")
    Logical  :: replace_geometry

    If(Present(input_fname)) Then
       fname = input_fname
    Else
       fname = dftb_fname
    EndIf

    !Read dftb input file into hsd_tree
    Call dftb%getInputFromFile(fname, hsd_tree)

    !Set structure data, retain rest
    Call hsd_tree%getRootNode(pRoot)
    Call setChild(pRoot, "Geometry", pGeo, replace=replace_geometry)
    Call setChildValue(pGeo, "Periodic", geo%tPeriodic)
    Call setChildValue(pGeo, "LatticeVectors", geo%latVecs)
    Call setChildValue(pGeo, "TypeNames", geo%speciesNames)

    !See DFTB+ subroutine in lib_type/typegeometryhsd.F90
    Call setChildValue(pGeo, "TypesAndCoordinates", &
         Reshape(geo%species, (/ 1, Size(geo%species) /)), geo%coords)

    !Always compute forces
    Call setChild(pRoot, "Analysis", pAnalysis, replace=.True.)
    Call setChildValue(pAnalysis, "CalculateForces", .True.)

  End Subroutine initialise_dftbplus_tree

#:if WITH_MPI
  !> Output forces per process
  subroutine output_forces_per_process(gradients, imd_lab)
    implicit none

    real(dp), allocatable, intent(in) :: gradients(:,:)
    character(len=2), intent(in) :: imd_lab
    integer :: io_unit, ia
    character(len=2) :: rank_lab

    write(rank_lab, '(I2)') rank

    open(newunit=io_unit, file='forces_MD'//trim(adjustl(imd_lab))&
        //'_rank'//trim(adjustl(rank_lab))//'.dat')
    write(io_unit, *) '# One line of header'
    do ia =1,size(gradients,2)
      write(io_unit, "(I5, 3F20.12)") ia, -1._dp * gradients(:,ia)
    enddo

    close(io_unit)

  end subroutine output_forces_per_process
#:endif 
  
  !> Count number of substrings on a line separated by whitespace
  function count_substrings(input_line) result(n_substrings)
    implicit none
    character(len=*), intent(in) :: input_line
    character(len=len(input_line))  :: line
    character(len=1), parameter  :: whitespace = ' '
    integer :: i, length, n_substrings

    ! Add 1 prefixed and 1 trailing white space so counter
    ! works correctly
    line = ' '//input_line
    length = len(trim(adjustl(line))) + 1

    n_substrings = 0
    do i = 2,length
       if( line(i:i) /= whitespace .and. line(i-1:i-1) == whitespace)then
          n_substrings = n_substrings + 1
       endif
    enddo

  end function count_substrings


  !> Get number of atoms from DFTB+ geometry file
  function get_number_of_atoms(fname, header) result(nAtoms)
    implicit none

    character(len=*), intent(in)  :: fname
    integer, optional, intent(in) :: header
    integer :: nAtoms

    character(len=1) :: boundary
    integer :: io_unit , i

    if(IO) then
       open(newunit=io_unit,file=fname)
       
       if(present(header))then
          do i=1,header
             read(io_unit,*)
          enddo
       endif
       
       read(io_unit,*) nAtoms, boundary
       
       close(io_unit)
    endif
#:if WITH_MPI
    call mpi_bcast(nAtoms, 1, MPI_INTEGER, master_id, MPI_COMM_WORLD, ierr)
#:endif
  end function get_number_of_atoms
  

  !> Read DFTB+ geometry file onto master_id
  subroutine read_in_geo(fname, geo, header)
    implicit none

    character(len=*), intent(in)  :: fname
    type(TGeometry),  intent(out) :: geo
    integer, optional, intent(in) :: header
        
    integer          :: ia,i,Nspecies
    character(len=1) :: boundary
    character(len=2),   allocatable :: species_list(:)
    character(len=100) :: line2
    integer :: io_unit 

    if(IO) then
       !Initialise geo components
       geo%nAtom = 0
       geo%nSpecies = 0
       
       ! Open file
       open(newunit=io_unit,file=fname)
       
       !Skip nlines of header
       if(present(header))then
          do i=1,header
             read(io_unit,*)
          enddo
       endif
       
       !Assign geo data from file
       !Line 1
       read(io_unit,*) geo%nAtom, boundary
       
       if(boundary=='S')then
          geo%tPeriodic = .true.
          geo%tFracCoord = .false.
       elseif(boundary=='F')then
          geo%tPeriodic  = .true.
          geo%tFracCoord = .true.
       elseif(boundary=='C')then
          geo%tPeriodic = .false.
          geo%tFracCoord = .false.
       endif
       
       !Line 2.
       read(io_unit, '(A)') line2
       Nspecies = count_substrings(line2)
       backspace io_unit
       
       allocate(species_list(Nspecies))
       species_list=' '
       read(io_unit,*) species_list(:)
       
       geo%nSpecies = Nspecies
       if(.not. allocated(geo%speciesNames)) allocate(geo%speciesNames(geo%nSpecies))
       geo%speciesNames = ' '
       do ia=1,geo%nSpecies
          geo%speciesNames(ia) = species_list(ia)
       enddo
       
       !Lines 3 -> 3+Natoms
       if(.not. allocated(geo%species))   allocate(geo%species(geo%nAtom))
       if(.not. allocated(geo%coords))    allocate(geo%coords(3, geo%nAtom))
       geo%species = 0
       geo%coords = 0
       do ia=1,geo%nAtom
          read(io_unit,*) i,geo%species(ia),geo%coords(1:3,ia)
          geo%coords(1:3,ia) = geo%coords(1:3,ia) * AA__Bohr
       enddo
       
       if(boundary /= 'C')then
          !Origin
          read(io_unit, *)
          
          !Lattice vectors in rows in .gen BUT as columns internally
          if(.not. allocated(geo%latVecs)) allocate(geo%latVecs(3,3))
          do ia=1,3
             read(io_unit,*) geo%latVecs(:,ia)
          enddo
          geo%latVecs(:,:) =  geo%latVecs(:,:) * AA__Bohr
       endif
       
       close(io_unit)
    endif
    
  end subroutine read_in_geo

#:if WITH_MPI
  !> Broadcast geometry from master_id to all processes 
  subroutine broadcast_geometry(comm, geo)

    !> MPI communicator
    integer, intent(in) :: comm
    !> Geometry: Filled on master_id, uninitialised on all other processes
    type(TGeometry),  intent(inout) :: geo

    integer :: ierr, lenSN, ia
    
    call mpi_bcast(geo%nAtom,      1, MPI_INTEGER, master_id, comm, ierr)
    call mpi_bcast(geo%tPeriodic,  1, MPI_LOGICAL, master_id, comm, ierr)
    call mpi_bcast(geo%tFracCoord, 1, MPI_LOGICAL, master_id, comm, ierr)
    call mpi_bcast(geo%nSpecies,   1, MPI_INTEGER, master_id, comm, ierr)

    if(.not. allocated(geo%speciesNames)) then
       allocate(geo%speciesNames(geo%nSpecies))
    endif

    lenSN = len(geo%speciesNames)
    call mpi_bcast(geo%speciesNames, lenSN * geo%nSpecies, MPI_CHARACTER, master_id, comm, ierr)

    if(.not. allocated(geo%species)) then
       allocate(geo%species(geo%nAtom))
    endif

    if(.not. allocated(geo%coords)) then
       allocate(geo%coords(3, geo%nAtom))
    endif

    call mpi_bcast(geo%species, geo%nAtom, MPI_INTEGER, master_id, comm, ierr)
    call mpi_bcast(geo%coords, 3 * geo%nAtom, MPI_DOUBLE_PRECISION, master_id, comm, ierr)

    if (geo%tPeriodic) then
       if(.not. allocated(geo%latVecs)) then
          allocate(geo%latVecs(3,3))
       endif
       call mpi_bcast(geo%latVecs, 9, MPI_DOUBLE_PRECISION, master_id, comm, ierr)
    endif

    call mpi_barrier(comm, ierr)
    
  end subroutine broadcast_geometry
#:endif
  
end program test_setSpeciesAndDependents
