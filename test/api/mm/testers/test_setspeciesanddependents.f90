!---------------------------------------------------------------------
! Test 'setSpeciesAndDependents' when order of species changes
!--------------------------------------------------------------------
!
! Test API subroutine 'setSpeciesAndDependents' which can
! be used to reset quantities that depend on the global ordering of 
! atomic species. 
!
! If species(1:nAtoms) = (/1,2,2,1,...,1/) -> (/1,2,2,1,...,2/),  
! following an MD step, reset equivalency relations and partial charges.
! This situation can arise in molecular dynamics packages that use domain
! decomposition as an MPI parallelisation scheme.
!
! All processes read in geometry and all processes contain all atoms.
! nAtoms is conserved.
!
! Forces are not used to update atomic positions, rather two sets
! of atomic positions are read from files and forces are computed
! for each of them, sequentially. Forces of the 2nd geometry are evaluated
! by the test framework. 
!
! If 'setSpeciesAndDependents' is not called, species data from structure1.gen
! will be used for the structure2.gen calculation, resulting in erroneous forces. 
!
! Author: A. Buccheri. 2020. University of Bristol.
! alexanderbuccheri@googlemail.com
!---------------------------------------------------------------------

program test_setSpeciesAndDependents
  use, intrinsic :: iso_fortran_env, only: output_unit, REAL64, IOSTAT_END
  use mpi 

  use dftbp_mmapi,  only: TDftbPlus_init, TDftbPlus_destruct, TDftbPlus, TDftbPlusInput
  use dftbp_hsdapi, only: fnode, getChild, getChildren, setChild, getChildValue, &
                          setChildValue, dumpHsd
  use testhelpers,  only: writeAutotestTag
  implicit none

  !> Precision and unit conversion
  integer, parameter :: dp = REAL64
  real(dp), parameter :: Bohr__AA = 0.529177249_dp
  real(dp), parameter :: AA__Bohr = 1.0_dp / Bohr__AA
  
  !> DFTB Objects
  type(TDftbPlus)      :: dftb
  type(TDftbPlusInput) :: hsd_tree, input
  Character(len=11), parameter :: dftb_fname = 'dftb_in.hsd'
  !> Print more data for debugging
  Logical              :: suppressIO = .true.

  !> Type containing geometrical information
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
  integer, parameter :: requiredThreading = MPI_THREAD_FUNNELED
  integer            :: providedThreading, rank, np, ierr
  integer, parameter :: master_id = 0
  logical            :: IO

  !Local
  integer, parameter :: nAtoms = 568
  integer, parameter :: Nsteps = 2  

  type(MDstatus_type):: MDstatus
  type(TGeometry)    :: geo
  integer            :: imd, ia, i
  real(dp)           :: merminEnergy
  character(len=2)   :: imd_lab
  character(len=100) :: fname
  real(dp), allocatable :: gradients(:,:)
  
  !Initialise MPI environment
  call mpi_init_thread(requiredThreading, providedThreading, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, np, ierr)
  IO = (rank == master_id)

  MDstatus = MDstatus_type(1, Nsteps)
  allocate(gradients(3,nAtoms))
  
  do imd = MDstatus%initial_step, MDstatus%final_step
     if(IO) write(*,*) 'MD step ', imd, 'of ', MDstatus%final_step

     !For every step, read in atomic positions from file
     !Equivalent to obtaining data from MD package 
     write(imd_lab, '(I2)') imd 
     fname = 'structure_'//trim(adjustl(imd_lab))//'.gen'
     call read_in_geo(trim(adjustl(fname)), geo, header=1)
     
     if(imd == MDstatus%initial_step) then
        Call TDftbPlus_init(dftb, output_unit, MPI_COMM_WORLD)
        Call initialise_dftbplus_tree(geo, dftb, hsd_tree)
        !Dump hsd tree to fort.1
        if(.not. suppressIO) Call dumpHsd(hsd_tree%hsdTree, 001)
        Call dftb%setupCalculator(hsd_tree)
     endif

     !Update coordinates and lattice vectors
     call dftb%setGeometry(geo%coords, geo%latVecs)
     
     !Update species order every step
     call dftb%setSpeciesAndDependents(&
          geo%speciesNames, geo%species) 
     
     !Do a total energy calculation
     call dftb%getEnergy(merminEnergy)     
     call dftb%getGradients(gradients)
     if(.not. suppressIO)then
        call output_forces_per_process(gradients, imd_lab)
     endif
  enddo

  !Clean up
  call TDftbPlus_destruct(dftb)

  ! Write file for internal test system                                                                    
  call writeAutotestTag(merminEnergy=merminEnergy, gradients=gradients) 

  call mpi_finalize(ierr)

  
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

    Type(tGeometry), Intent( In    ) :: geo
    Type(TDftbPlus),          Intent( InOut ) :: dftb
    Type(TDftbPlusInput),     Intent(   Out ) :: hsd_tree
    Character(Len=*),         Intent( In    ),   Optional :: input_fname

    !> DFTB+ input file name 
    Character(Len=100) :: fname
    !> Pointers to the parts of the input tree that will be set                                
    Type(fnode), Pointer :: pRoot, pGeo, pOptions, pParserOpts, pDftb, pKpoints, pAnalysis
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

    !More info for debugging 
    If(.not. suppressIO) Then
       !Call setChildValue(pAnalysis, "WriteBandOut", .False.)
       !Call setChildValue(pAnalysis, "MullikenAnalysis", .True.)
       Call setChild(pRoot, "ParserOptions", pParserOpts, replace=.True.)
       Call setChildValue(pParserOpts, "ParserVersion", 8)
       Call setChildValue(pParserOpts, 'WriteResultsTag', .True.)
       Call setChildValue(pParserOpts, "WriteDetailedOut", .True.)
       Call setChildValue(pParserOpts, "WriteHSDInput", .True.)
    Endif
      
  End Subroutine initialise_dftbplus_tree

  
  !> Output forces per process
  subroutine output_forces_per_process(gradients, imd_lab)
    implicit none
    
    real(dp), allocatable, intent(in) :: gradients(:,:)
    character(len=2), intent(in) :: imd_lab 
    integer :: io_unit, ia
    character(len=2) :: rank_lab

    !Ensure each MPI process uses a different unit when writing 
    io_unit = int(100 * rank+1)
    write(rank_lab, '(I2)') rank 

    open(unit=io_unit, file='forces_MD'//trim(adjustl(imd_lab))&
         //'_rank'//trim(adjustl(rank_lab))//'.dat')
    write(io_unit, *) '# One line of header'
    do ia =1,size(gradients,2)
       write(io_unit, "(I5, 3F20.12)") ia, -1._dp * gradients(:,ia)
    enddo

    close(io_unit)
    
  end subroutine output_forces_per_process

  
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

  
  !> Read DFTB+ geometry file 
  subroutine read_in_geo(fname, geo, header)
    implicit none
 
    type(TGeometry),  intent(out) :: geo
    character(len=*), intent(in)    :: fname
    integer          :: ia,i,Nspecies
    character(len=1) :: boundary
    character(len=2),   allocatable :: species_list(:)
    character(len=100) :: line2
    integer, optional, intent(in) :: header

    !Initialise geo components
    geo%nAtom = 0
    geo%nSpecies = 0
   
    !Open file
    open(unit=100,file=fname)
    
    !Skip nlines of header
    if(present(header))then
       do i=1,header
          read(100,*)
       enddo
    endif
    
    !Assign geo data from file
    !Line 1
    read(100,*) geo%nAtom, boundary
   
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
    read(100, '(A)') line2
    Nspecies = count_substrings(line2) 
    backspace 100
    
    allocate(species_list(Nspecies))
    species_list=' '
    read(100,*) species_list(:)
    
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
       read(100,*) i,geo%species(ia),geo%coords(1:3,ia)
       geo%coords(1:3,ia) = geo%coords(1:3,ia) * AA__Bohr
    enddo

    if(boundary /= 'C')then
       !Origin
       read(100, *) 
       
       !Lattice vectors in rows in .gen BUT as columns internally
       if(.not. allocated(geo%latVecs)) allocate(geo%latVecs(3,3))
       do ia=1,3
          read(100,*) geo%latVecs(:,ia)
       enddo
       geo%latVecs(:,:) =  geo%latVecs(:,:) * AA__Bohr
    endif
    
    close(100)
       
  end subroutine read_in_geo

end program test_setSpeciesAndDependents
