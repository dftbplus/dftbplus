module libnegf_vars
  use Accuracy, only : dp, mc, lc
  use CommonTypes
  use WrappedIntrinsics
  use xmlf90                                                                !DAR                 
  implicit none
  private

  public :: TNEGFStructure
  public :: TNEGFTunDos
  public :: TNEGFGreenDensInfo
  public :: TTransPar
  public :: ContactInfo
  public :: TElph

  type TNEGFStructure
    integer               :: nAtom          ! number of atoms in central cell 
    integer, pointer      :: iatomstart(:)  ! atom START pos for squared H/S
  end type TNEGFStructure


  !> Options for electron-phonon model
  type TElPh
    logical :: defined = .false. !True if filled up with info from an input block
    integer :: model = 0  !Specify which model in input: 1=elastic
    real(dp), allocatable, dimension(:) :: coupling
    integer  :: scba_niter = 0
    integer, allocatable, dimension(:) :: orbsperatm !List of orbital per atom
                                            ! for model = 2 (atom block local)
  end type TElPh

  !Structure for contact information in transport calculation
  type ContactInfo
    integer :: idxrange(2) ! Beginning (1) and end (2) of contact
    character(mc) :: name
    ! Note: I explicitely define a contact id because with multiple definition
    ! of contacts in the input file relying on contact ordering to assign an
    ! integer can be inconsistent
    real(dp) :: shiftAccuracy   = 0.0! Accuracy of rigid layer shift
    integer :: dir              = 0
    real(dp) :: length          = 0.0
    real(dp) :: lattice(3)    ! Lattice vectors
    real(dp) :: potential       = 0.0
    ! for colinear spin we may need two fermi level (up and down)
    real(dp) :: eFermi(2)       = (/0.0_dp, 0.0_dp /)
    real(dp) :: kbT             = 0.0
    ! Is it a contact in wide band approximation?
    logical :: wideBand         = .false.
    real(dp) :: wideBandDos     = 0.0
    character(lc) :: output   ! Filename for contact infos (shiftcont_) TO BE
                              !MOVED?
    !DAR begin - Write/Read SE, GF
    logical :: tWriteSelfEnergy = .false.                                   
    logical :: tReadSelfEnergy = .false.                                    
    logical :: tWriteSurfaceGF = .false.                                   
    logical :: tReadSurfaceGF = .false.                                    
    !DAR end
  end type ContactInfo


  !Options for Landauer (Tunneling and Dos) calculation
  type TNEGFTunDos
    logical            :: defined = .false.     ! true only if filling block is
                                               ! defined
    integer            :: verbose              ! verbosity level of the library
    integer            :: gSpin                ! spin degeneracy (used in transmission and current integration)
    real(dp)           :: emin                 ! Min integration energy
                                               ! (possible to define them different for colinear spin calculation)
    real(dp)           :: emax                 ! Max integration energy
    real(dp)           :: estep                ! Energy step
    real(dp)           :: delta                ! Delta for Green's function
    real(dp)           :: broadeningDelta                ! An additional broadening delta for DOS and tunneling
    integer, allocatable, dimension(:)  :: ni  !emitter contact(s)
    integer, allocatable, dimension(:)  :: nf  !collector contact(s)
    type(WrappedInt1), allocatable :: dosOrbitals(:)
    character(lc), allocatable :: dosLabels(:)
    logical :: writeLDOS = .false.  ! write DOS on separate files
    logical :: writeTunn = .false.  ! write tunneling on separate files
    real(dp), dimension(:), allocatable :: kbT ! contact temperatures 
    type(Telph) :: elph
    type(Telph) :: bp
  end type TNEGFTunDos

  type TNEGFGreenDensInfo
    ! Information for Green's function charge density section
    logical            :: defined = .false.     ! true only if filling block is
    integer            :: verbose              ! verbosity level of the library
    ! Fermi level for closed system calculation. If a coliner spin calculation 
    ! is defined, two values are needed (up and down)
    real(dp)   :: oneFermi(2) = (/0.0_dp, 0.0_dp /) ! unique Fermi closed systems
    real(dp)           :: delta                ! delta function in G.F.
    integer            :: nP(3)                ! Number of points in contour
    real(dp)           :: enLow                ! Lowest energy for contour int
    integer            :: nkT                  ! N. of kT for Fermi dist
    integer            :: nPoles               ! N. of poles included in contour
    logical            :: doGreenDens = .false.! use or not Green solver
    logical            :: saveSGF              ! save SGF on files
    logical            :: readSGF              ! read SGF from files
    logical            :: doLocalCurr = .false.! Calculate or not local J
    !! There is an independent definition of pls as in closed system 
    !! green's calculation a separate definition may be used
    integer            :: nPLs = 0             ! N. of principal layers
    integer, allocatable   :: PL(:)            ! PL indeces (starting atom)
    integer            :: gSpin                ! spin degeneracy (used in charge integration)
    real(dp), dimension(:), allocatable :: kbT ! contact temperatures 
    type(Telph) :: elph
    type(Telph) :: bp
  end type TNEGFGreenDensInfo
  
  
  ! Options from Transport section (geometry and task)  
  type TTransPar
     
    logical :: defined = .false.   ! True if the corresponding input block exists
    type(ContactInfo), dimension(:), allocatable :: contacts
    
    !DAR begin - type TTransPar new items
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical :: tNoGeometry = .false.
    logical :: tWriteDFTB = .false.
    logical :: tReadDFTB = .false.
    logical :: tOrthonormal = .false.
    logical :: tOrthonormalDevice = .false.
    logical :: tModel = .false.
    logical :: tRead_negf_in = .false.
    integer :: NumStates = 0
    integer, allocatable :: cblk(:)
    character(lc) :: FileName
    logical :: tManyBody =.false.
    logical :: tElastic =.true.
    logical :: tDephasingVE = .false.
    logical :: tDephasingBP = .false.
    logical :: tZeroCurrent = .false.
    integer :: MaxIter = 1000
    logical :: tWriteDOS = .false.
    logical :: tWrite_ldos = .false.
    logical :: tWrite_negf_params = .false.
    logical :: tDOSwithS =.true.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !DAR end
    
    integer :: ncont = 0
    !! Start and end index of device region
    integer :: idxdevice(2)
    integer :: nPLs =1                ! N. of principal layers
    integer, allocatable, dimension(:) :: PL   ! PL indeces (starting atom)
    ! False: run the full OBC calculation / True: upload contact phase
    logical :: taskUpload = .false.
    ! Index of contact for contact hamiltonian task, if any
    integer :: taskContInd = 0
    !! Not from input file
    logical :: tPeriodic1D = .false.
    
  end type TTransPar
    


end module libnegf_vars
