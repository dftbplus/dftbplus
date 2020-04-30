!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

module libnegf_vars
  use dftbp_accuracy, only : dp, mc, lc
  use dftbp_commontypes
  use dftbp_wrappedintr
  use dftbp_xmlf90
  implicit none
  private

  public :: TNEGFTunDos
  public :: TNEGFGreenDensInfo
  public :: TTransPar
  public :: ContactInfo
  public :: TElph



  !> Options for electron-phonon model
  type TElPh

    !> True if filled up with info from an input block
    logical :: defined = .false.

    !> Specify which model in input, 1-3 elastic models: diagonal, block diagonal and overlap masked
    integer :: model = 0

    !> Coupling strength (?)
    real(dp), allocatable :: coupling(:)

    !> Iterations for self-consistent Born approximation
    integer  :: scba_niter = 0

    !> List of orbital per atom for models = (2,3)
    integer, allocatable :: orbsperatm(:)

  end type TElPh


  !Structure for contact information in a transport calculation
  type ContactInfo

    ! Beginning (1) and end (2) of contact in atoms (?)
    integer :: idxrange(2)

    !> Note: a contact id is specifically defined because, with multiple definition of contacts in
    !> the input file, relying on contact ordering to assign an integer can be inconsistent
    character(mc) :: name

    !> Accuracy of rigid layer shift
    real(dp) :: shiftAccuracy = 0.0_dp

    integer :: dir = 0

    real(dp) :: length = 0.0_dp

    !> Lattice vectors
    real(dp) :: lattice(3)

    real(dp) :: potential = 0.0_dp

    !> for colinear spin we may need two Fermi levels (up and down)
    real(dp) :: eFermi(2) = [0.0_dp, 0.0_dp]

    !> contact temperature
    real(dp) :: kbT = 0.0_dp

    ! Is it a contact in the wide band approximation?
    logical :: wideBand = .false.

    real(dp) :: wideBandDos = 0.0_dp

    !> Filename for contact infos (shiftcont_) TO BE MOVED?
    character(lc) :: output

    logical :: tWriteSelfEnergy = .false.
    logical :: tReadSelfEnergy = .false.
    logical :: tWriteSurfaceGF = .false.
    logical :: tReadSurfaceGF = .false.

  end type ContactInfo


  !> Options for Landauer (Tunneling and DOS) calculation
  type TNEGFTunDos

    !> true only if filling block is defined
    logical            :: defined = .false.

    !> verbosity level of the library
    integer            :: verbose

    !> spin degeneracy (used in transmission and current integration)
    integer            :: gSpin

    !> Min integration energy (possible to define them different for colinear spin calculation)
    real(dp) :: emin

    !> Max integration energy
    real(dp)           :: emax

    !> Energy step
    real(dp)           :: estep

    !> Delta for Green's function
    real(dp)           :: delta

    !> An additional broadening delta for DOS and tunneling
    real(dp)           :: broadeningDelta

    !> emitter contact(s)
    integer, allocatable  :: ni(:)

    !> collector contact(s)
    integer, allocatable  :: nf(:)

    type(TWrappedInt1), allocatable :: dosOrbitals(:)

    character(lc), allocatable :: dosLabels(:)

    !> write DOS on separate files
    logical :: writeLDOS = .false.

    !> write tunneling on separate files
    logical :: writeTunn = .false.

    !> contact temperatures
    real(dp), allocatable :: kbT(:)

    !> Electron-phonon coupling
    type(Telph) :: elph

    !> Buttiker Probe for dephasing
    type(Telph) :: bp

  end type TNEGFTunDos


  !> Information for Green's function charge density calculation
  type TNEGFGreenDensInfo

    !> true only if filling block is defined
    logical            :: defined = .false.

    !> verbosity level of the library
    integer            :: verbose

    !> Fermi level for closed system calculation. If a coliner spin calculation is defined, two
    !> values are needed (up and down) unique Fermi closed systems
    real(dp)   :: oneFermi(2) = [0.0_dp, 0.0_dp]

    !> delta function in G.F.
    real(dp)           :: delta

    !> Number of points in contour
    integer            :: nP(3)

    !> Lowest energy for contour int
    real(dp)           :: enLow

    !> Number of kT for Fermi dist
    integer            :: nkT

    !> Number of poles included in contour
    integer            :: nPoles

    !> use or not Green solver
    logical            :: doGreenDens = .false.

    !> save SGF in files
    logical            :: saveSGF

    !> read SGF from files
    logical            :: readSGF

    !> Calculate or not the local J. There is an independent definition of principal layers (pls),
    !> since in a closed system Green's calculation a separate definition may be used
    logical            :: doLocalCurr = .false.

    !> Number of principal layers
    integer            :: nPLs = 0

    !> PL indices (starting atom)
    integer, allocatable  :: PL(:)

    !> spin degeneracy (used in charge integration)
    integer            :: gSpin

    !> contact temperatures
    real(dp), allocatable :: kbT(:)

    !> Electron-phonon coupling
    type(Telph) :: elph

    !> Buttiker Probe for dephasing
    type(Telph) :: bp

  end type TNEGFGreenDensInfo


  !> Options from Transport section (geometry and task)
  type TTransPar

    !> True if the corresponding input block exists
    logical :: defined = .false.

    !> Contacts in the system
    type(ContactInfo), allocatable :: contacts(:)

    !> Number of contacts
    integer :: ncont = 0

    !> Start and end index of device region
    integer :: idxdevice(2)

    !> Number of principal layers
    integer :: nPLs =1

    !> PL indices (starting atom)
    integer, allocatable, dimension(:) :: PL

    !> False: run the full OBC calculation / True: upload contact phase
    logical :: taskUpload = .false.

    !> Index of contact for contact hamiltonian task, if any
    integer :: taskContInd = 0

    !> Not from input file
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

  !> Copies contents of a Green function density calculation structure
  subroutine copyGreenDens(gIN, gOUT)

    !> Original structure
    type(TNEGFGreenDensInfo), intent(in) :: gIN

    !> Duplicate
    type(TNEGFGreenDensInfo), intent(out) :: gOUT

    integer :: isz

    ! copy greendens
    gOUT%defined = gIN%defined
    gOUT%verbose = gIN%verbose
    gOUT%oneFermi = gIN%oneFermi
    gOUT%delta = gIN%delta
    gOUT%np = gIN%np
    gOUT%enLow = gIN%enLow
    gOUT%nkT = gIN%nkT
    gOUT%nPoles = gIN%nPoles
    gOUT%doGreenDens= gIN%doGreenDens
    gOUT%saveSGF = gIN%saveSGF
    gOUT%readSGF = gIN%readSGF
    gOUT%doLocalCurr= gIN%doLocalCurr
    gOUT%nPLs = gIN%nPLs
    gOUT%gspin = gIN%gspin

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
    gOUT%elph%defined = gIN%elph%defined
    gOUT%elph%model = gIN%elph%model
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
    gOUT%bp%defined = gIN%bp%defined
    gOUT%bp%model = gIN%bp%model
    gOUT%bp%scba_niter = gIN%bp%scba_niter

  end subroutine copyGreenDens


  !> Copies contents of a Landauer calculation structure
  subroutine copyTunDos(gIN, gOUT)

    !> Original structure
    type(TNEGFTunDos), intent(in) :: gIN

    !> Duplicate
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

    gOUT%defined = gIN%defined
    gOUT%verbose = gIN%verbose
    gOUT%gspin = gIN%gspin
    gOUT%emin = gIN%emin
    gOUT%emax = gIN%emax
    gOUT%estep = gIN%estep
    gOUT%delta = gIN%delta
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
    gOUT%elph%defined = gIN%elph%defined
    gOUT%elph%model = gIN%elph%model
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
    gOUT%bp%defined = gIN%bp%defined
    gOUT%bp%model = gIN%bp%model
    gOUT%bp%scba_niter = gIN%bp%scba_niter

  end subroutine copyTunDos


  !> Copies contents of a transport calculation structure
  subroutine copyTranspar(tIN, tOUT)

    !> Original structure
    type(TTransPar), intent(in) :: tIN

    !> Duplicate
    type(TTransPar), intent(out) :: tOUT

    integer :: isz

    tOUT%defined = tIN%defined
    tOUT%ncont = tIN%ncont
    tOUT%idxdevice = tIN%idxdevice
    tOUT%nPls = tIN%nPls
    tOUT%taskUpload = tIN%taskUpload
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
