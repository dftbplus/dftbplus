module libnegf_vars
  use Accuracy, only : dp, mc, lc
  use CommonTypes
  use WrappedIntrinsics
  use xmlf90                                                                !DAR                 
  implicit none
  private

  public :: TNEGFTunDos
  public :: TNEGFGreenDensInfo
  public :: TTransPar
  public :: ContactInfo
  public :: TElph



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
    integer, allocatable  :: PL(:)             ! PL indeces (starting atom)
    integer            :: gSpin                ! spin degeneracy (used in charge integration)
    real(dp), allocatable :: kbT(:)            ! contact temperatures 
    type(Telph) :: elph
    type(Telph) :: bp
  end type TNEGFGreenDensInfo
  
  
  ! Options from Transport section (geometry and task)  
  type TTransPar
     
    logical :: defined = .false.   ! True if the corresponding input block exists
    type(ContactInfo), dimension(:), allocatable :: contacts
    ! number of contacts
    integer :: ncont = 0
    !! Start and end index of device region
    integer :: idxdevice(2)
    ! N. of principal layers
    integer :: nPLs =1                
    ! PL indeces (starting atom)
    integer, allocatable, dimension(:) :: PL   
    ! False: run the full OBC calculation / True: upload contact phase
    logical :: taskUpload = .false.
    ! Index of contact for contact hamiltonian task, if any
    integer :: taskContInd = 0
    !! Not from input file
    logical :: tPeriodic1D = .false.
    
    
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
    
  end type TTransPar
  
  contains 

  subroutine copyGreenDens(gIN, gOUT)
    type(TNEGFGreenDensInfo), intent(in) :: gIN
    type(TNEGFGreenDensInfo), intent(out) :: gOUT

    integer :: isz

    !copy greendens
    gOUT%defined    = gIN%defined
    gOUT%verbose    = gIN%verbose
    gOUT%oneFermi   = gIN%oneFermi   
    gOUT%delta      = gIN%delta      
    gOUT%np         = gIN%np         
    gOUT%enLow      = gIN%enLow      
    gOUT%nkT        = gIN%nkT        
    gOUT%nPoles     = gIN%nPoles     
    gOUT%doGreenDens= gIN%doGreenDens
    gOUT%saveSGF    = gIN%saveSGF    
    gOUT%readSGF    = gIN%readSGF    
    gOUT%doLocalCurr= gIN%doLocalCurr
    gOUT%nPLs       = gIN%nPLs       
    gOUT%gspin      = gIN%gspin      

    if (allocated(gIN%kbT)) then
      isz = size(gIN%kbT)
      allocate(gOUT%kbT(isz))
      gOUT%kbT = gIN%kbT
    end if
    if (allocated(gIN%PL)) then
      isz = size(gIN%PL)
      allocate(gOUT%PL(isz))
      gOUT%PL = gIN%PL
    end if
    ! copy greendens%elph
    if (allocated(gIN%elph%coupling)) then
      isz = size(gIN%elph%coupling)
      allocate(gOUT%elph%coupling(isz))
      gOUT%elph%coupling = gIN%elph%coupling
    end if
    if (allocated(gIN%elph%orbsperatm)) then
      isz = size(gIN%elph%orbsperatm)
      allocate(gOUT%elph%orbsperatm(isz))
      gOUT%elph%orbsperatm = gIN%elph%orbsperatm
    end if
    gOUT%elph%defined   = gIN%elph%defined  
    gOUT%elph%model     = gIN%elph%model    
    gOUT%elph%scba_niter = gIN%elph%scba_niter

    ! copy greendens%bp
    if (allocated(gIN%bp%coupling)) then
      isz = size(gIN%bp%coupling)
      allocate(gOUT%bp%coupling(isz))
      gOUT%bp%coupling = gIN%bp%coupling
    end if
    if (allocated(gIN%bp%orbsperatm)) then
      isz = size(gIN%bp%orbsperatm)
      allocate(gOUT%bp%orbsperatm(isz))
      gOUT%bp%orbsperatm = gIN%bp%orbsperatm
    end if
    gOUT%bp%defined   = gIN%bp%defined  
    gOUT%bp%model     = gIN%bp%model    
    gOUT%bp%scba_niter = gIN%bp%scba_niter
 
  end subroutine copyGreenDens  

  subroutine copyTunDos(gIN, gOUT)
    type(TNEGFTunDos), intent(in) :: gIN
    type(TNEGFTunDos), intent(out) :: gOUT

    integer :: isz

    ! copy tundos
    if (allocated(gIN%ni)) then
      isz = size(gIN%ni)
      allocate(gOUT%ni(isz))
      gOUT%ni = gIN%ni
    end if
    if (allocated(gIN%nf)) then
      isz = size(gIN%nf)
      allocate(gOUT%nf(isz))
      gOUT%nf = gIN%nf
    end if
    if (allocated(gIN%kbT)) then
      isz = size(gIN%kbT)
      allocate(gOUT%kbT(isz))
      gOUT%kbT = gIN%kbT
    end if
    if (allocated(gIN%dosOrbitals)) then
      isz = size(gIN%dosOrbitals)
      allocate(gOUT%dosOrbitals(isz))
      gOUT%dosOrbitals = gIN%dosOrbitals
    end if
    if (allocated(gIN%dosLabels)) then
      isz = size(gIN%dosLabels)
      allocate(gOUT%dosLabels(isz))
      gOUT%dosLabels = gIN%dosLabels
    end if
     
    gOUT%defined   = gIN%defined
    gOUT%verbose   = gIN%verbose
    gOUT%gspin     = gIN%gspin      
    gOUT%emin      = gIN%emin       
    gOUT%emax      = gIN%emax       
    gOUT%estep     = gIN%estep      
    gOUT%delta     = gIN%delta      
    gOUT%writeLDOS = gIN%writeLDOS  
    gOUT%writeTUNN = gIN%writeTUNN  
    gOUT%broadeningDelta = gIN%broadeningDelta

    ! copy tundos%elph
    if (allocated(gIN%elph%coupling)) then
      isz = size(gIN%elph%coupling)
      allocate(gOUT%elph%coupling(isz))
      gOUT%elph%coupling = gIN%elph%coupling
    end if
    if (allocated(gIN%elph%orbsperatm)) then
      isz = size(gIN%elph%orbsperatm)
      allocate(gOUT%elph%orbsperatm(isz))
      gOUT%elph%orbsperatm = gIN%elph%orbsperatm
    end if
    gOUT%elph%defined   = gIN%elph%defined  
    gOUT%elph%model     = gIN%elph%model    
    gOUT%elph%scba_niter = gIN%elph%scba_niter

    ! copy tundos%bp
    if (allocated(gIN%bp%coupling)) then
      isz = size(gIN%bp%coupling)
      allocate(gOUT%bp%coupling(isz))
      gOUT%bp%coupling = gIN%bp%coupling
    end if
    if (allocated(gIN%bp%orbsperatm)) then
      isz = size(gIN%bp%orbsperatm)
      allocate(gOUT%bp%orbsperatm(isz))
      gOUT%bp%orbsperatm = gIN%bp%orbsperatm
    end if
    gOUT%bp%defined   = gIN%bp%defined  
    gOUT%bp%model     = gIN%bp%model    
    gOUT%bp%scba_niter = gIN%bp%scba_niter

  end subroutine copyTunDos

  subroutine copyTranspar(tIN, tOUT)
    type(TTransPar), intent(in) :: tIN    
    type(TTransPar), intent(out) :: tOUT  
    
    integer :: isz

    tOUT%defined     = tIN%defined    
    tOUT%ncont       = tIN%ncont      
    tOUT%idxdevice   = tIN%idxdevice  
    tOUT%nPls        = tIN%nPls       
    tOUT%taskUpload  = tIN%taskUpload 
    tOUT%taskContInd = tIN%taskContInd
    tOUT%tPeriodic1D = tIN%tPeriodic1D
    if (allocated(tIN%PL)) then
      isz = size(tIN%PL)
      allocate(tOUT%PL(isz))
      tOUT%PL = tIN%PL
    end if
    if (allocated(tIN%contacts)) then
      isz = size(tIN%contacts)
      allocate(tOUT%contacts(isz))
      tOUT%contacts = tIN%contacts
    end if
        
  end subroutine copyTranspar


end module libnegf_vars
