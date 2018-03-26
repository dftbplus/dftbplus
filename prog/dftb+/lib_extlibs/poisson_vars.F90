module poisson_vars
  use Accuracy, only : dp 
  use CommonTypes
  implicit none
  private

  public :: TPoissonInfo 
  public :: TPoissonStructure
  public :: TSKdata
  
  type TPoissonStructure
    integer               :: nAtom          ! number of atoms in central cell 
    integer               :: nSpecies       ! number of species
    integer, pointer      :: specie0(:)     ! type of the atoms (nAtom)
    integer, pointer      :: iatomstart(:)  ! atom START pos for squared H/S
    real(dp), pointer     :: x0(:,:)        ! coordinates in central cell
    real(dp)              :: nel            ! total number of electrons
    real(dp)              :: latVecs(3,3)   ! lattice vectors
    real(dp)              :: tempElec       ! electron temperature
    logical               :: isperiodic     ! tells whether the system is periodic    
  end type TPoissonStructure

  type TSKdata
    type(TOrbitals), pointer :: orb        !* Information about orbitals
    real(dp), pointer    :: hubbU(:,:)     !* Hubbard Us (orbital, atom)
    real(dp)             :: mCutoff        !* longest pair interaction
  end type TSKdata
    
  type TPoissonInfo
 
    ! Informations for the Poisson section
    integer            :: verbose              ! verbosity level of the library
    logical            :: defined = .false.  ! solve or not Poisson
    real(dp)           :: poissBox(3)          ! Poisson box
    real(dp)           :: poissGrid(3)         ! Minimal grid spacing
    logical            :: foundBox             ! (.false. for periodic systems!)
    real(dp)           :: maxRadAtomDens       ! Maximum radius of atom density
    real(dp)           :: poissAcc             ! solution accuracy
    logical            :: bulkBC               ! use bulk potential as BC  
    logical            :: readBulkPot          ! read bulk potential from file
    character(1)       :: localBCType          ! activates local BC mode (C|S)
    real(dp)           :: bufferLocBC          ! buffer spacing of local BC 
    integer            :: overrideBC(6)        ! forced BC in each direction
    integer            :: overrBulkBC(6)       ! forced BC on bulk potential
    logical            :: savePotential = .false. ! save the potential on a file
    logical            :: solveTwice  = .false.  ! Recompute poisson after D.M.    
    integer            :: maxPoissIter         ! maximum number of poisson iter
    character(1)       :: gateType             ! Planar| Cyilindrical Gate
    integer            :: gatedir              ! gate direction
    real(dp)           :: gatePot              ! Gate potential
    real(dp)           :: gateLength_l         ! Gate length along transport dir 
    real(dp)           :: gateLength_t         ! Gate length in transv. dir
    real(dp)           :: insLength            ! Insulator length
    real(dp)           :: insRad               ! Radius of insulator
    real(dp)           :: gateRad              ! Radius of gate
    real(dp)           :: eps_r                ! Insulator relative dielect.
    real(dp)           :: dr_eps               ! Buffer layer between dielect
                                               ! and vacuum
    real(dp)           :: bufferBox            ! Box buffer inside the contact region
    logical :: exactRenorm                     ! Use new numerical renormalization
                                               ! volume (preserves total charge)
    logical :: cutoffcheck 
    integer :: maxNumNodes                     ! Number of Nodes for parallel P. 

    character(:), allocatable :: scratch

  end type TPoissonInfo


end module poisson_vars
