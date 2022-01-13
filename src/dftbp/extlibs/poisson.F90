!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

#:assert not (WITH_POISSON and INSTANCE_SAFE_BUILD)

!> Interface to libPoisson routines
!>
!> NOTE: THIS MODULE IS NOT MULTI-INSTANCE SAFE
!>
module dftbp_extlibs_poisson
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  use dftbp_common_environment, only : TEnvironment, globalTimers
  use dftbp_common_globalenv, only : stdOut
  use dftbp_io_message, only : error
  use dftbp_poisson_poisson, only : poiss_savepotential, poiss_updcoords, active_id, natoms,&
      & verbose, bufferBox, deltaR_max, DoCilGate, DoGate, dR_cont, dr_eps, eps_r, fixed_renorm,&
      & FoundBox, Gate, GateDir, GateLength_l, GateLength_t, id0, InitPot, localBC, MaxPoissIter,&
      & numprocs, overrBulkBC, overrideBC, OxLength, period, ReadBulk, Rmin_Gate, Rmin_Ins,&
      & SavePot, scratchfolder, uhubb, lmax, izp, dQmat, init_poissbox, mudpack_drv,&
      & poiss_freepoisson, set_scratch, init_structure, init_skdata, init_charges, init_defaults,&
      & set_temperature, set_ncont, set_cluster, set_mol_indeces, set_dopoisson, set_poissonbox,&
      & set_poissongrid, set_accuracy, set_verbose, check_biasdir, check_poisson_box,&
      & check_parameters, check_localbc, check_contacts, write_parameters, poiss_getlatvecs
  use dftbp_type_commontypes, only : TOrbitals
#:if WITH_MPI
  use dftbp_poisson_poisson, only : global_comm, poiss_mpi_init, poiss_mpi_split
  use libmpifx_module, only : mpifx_barrier, mpifx_bcast
#:endif
#:if WITH_TRANSPORT
  use dftbp_poisson_poisson, only : ncont, set_cont_indeces, set_contdir, set_fermi,&
      & set_potentials, set_builtin
  use dftbp_transport_negfvars, only : TTransPar
#:endif
  implicit none

  private
  public :: withPoisson
  public :: TPoissonInfo, TPoissonStructure
  public :: TPoissonInput, TPoisson, TPoisson_init

  logical, parameter :: withPoisson = ${FORTRAN_LOGICAL(WITH_POISSON)}$


#:if WITH_POISSON

  !> Contains some part of the Poisson solvers internal data
  !>
  !> Note: As much of the data of the solver is kept in global variables, the solver is not
  !> multi-instance safe!
  !>
  type :: TPoisson
    private

    ! Stores last calculated shell resolved potential
    real(dp), allocatable :: shellPot_(:,:)

    ! Stores the shift vector to use in order to upload the shell potential
    real(dp), allocatable :: shellPotUpload_(:,:)

  contains

    !> Updates the lattice vectors.
    procedure, nopass :: updateLatVecs

    !> Updates the coordinates
    procedure, nopass :: updateCoords

    !> Updates the charges.
    procedure :: updateCharges

    !> Adds the potential calculated with the last charges.
    procedure :: addPotentials

    !> Returns the gradients corresponding to the last charges.
    procedure :: getGradients

    !> Saves the internal potential of the Poisson solver.
    procedure, nopass :: savePotential

    ! Finalizer
    final :: finalize_

  end type TPoisson

#:else

  type :: TPoisson
  end type TPoisson

#:endif


#:if WITH_POISSON

  !> Geometry of the atoms for the Poisson solver
  type :: TPoissonStructure

    !> number of atoms in central cell
    integer :: nAtom

    !> number of species
    integer :: nSpecies

    !> type of the atoms (nAtom)
    integer, allocatable :: specie0(:)

    !> coordinates in central cell
    real(dp), allocatable :: x0(:,:)

    !> lattice vectors
    real(dp), allocatable :: latVecs(:,:)

    !> tells whether the system is periodic
    logical :: isperiodic

  end type TPoissonStructure

#:else

  type :: TPoissonStructure
  end type TPoissonStructure

#:endif


#:if WITH_POISSON

  !> Information for the Poisson solver
  type TPoissonInfo

    !> verbosity level of the library
    integer :: verbose

    !> solve or not Poisson
    logical :: defined = .false.

    !> Poisson box
    real(dp) :: poissBox(3)

    !> Minimal grid spacing
    real(dp) :: poissGrid(3)

    !> Has a containing box been identified for the structure? (.false. for periodic systems!)
    logical :: foundBox

    !> Maximum radius of atom charge density
    real(dp) :: maxRadAtomDens

    !> Required solution accuracy
    real(dp) :: poissAcc

    !> use bulk potential as boundary condition
    logical :: bulkBC

    !> read bulk potential from file
    logical :: readBulkPot

    !> activates local boundary conditions mode (C|S)
    character(1) :: localBCType

    !> buffer spacing of local boundary condition
    real(dp) :: bufferLocBC

    !> forced boundary conditions in each direction
    integer :: overrideBC(6)

    !> forced boundary conditions on bulk potential
    integer :: overrBulkBC(6)

    !> save the potential on a file
    logical :: savePotential = .false.

    !> maximum number of poisson iter
    integer :: maxPoissIter

    !> Planar or cylindrical gate
    character(1) :: gateType

    !> gate direction
    integer :: gatedir

    !> Gate potential
    real(dp) :: gatePot

    !> Gate length along transport direction
    real(dp) :: gateLength_l

    !> Gate length in transverse direction
    real(dp) :: gateLength_t

    !> Insulator length
    real(dp) :: insLength

    !> Radius of insulator
    real(dp) :: insRad

    !> Radius of gate
    real(dp) :: gateRad

    !> Insulator relative dielect.
    real(dp) :: eps_r

    !> Buffer layer between dielectric and vacuum
    real(dp) :: dr_eps

    !> Box buffer inside the contact region
    real(dp) :: bufferBox

    !> Use new numerical renormalization volume (preserves total charge)
    logical :: numericNorm

    !> Check whether density cut off fits into the PLs
    logical :: cutoffcheck

    !> Number of Nodes for parallel P.
    integer :: maxNumNodes

    !> scratch folder name
    character(:), allocatable :: scratch

  end type TPoissonInfo

#:else

  type :: TPoissonInfo
  end type TPoissonInfo

#:endif


#:if WITH_POISSON

  !> Contains the main input for the Poisson-solver
  type :: TPoissonInput

    !> Geometrical structure related input data
    type(TPoissonStructure) :: poissonStruct

    !> Poisson solver specific parameters
    type(TPoissonInfo) :: poissonInfo

  #:if WITH_TRANSPORT

    !> Optional transport parameters (only when used in transport calculations)
    type(TTransPar), allocatable :: transPar

  #:endif

    !> Hubbard U values (uncontracted). Shape: [mShell, nSpecies]
    real(dp), allocatable :: hubbU(:,:)

    !> Optional shift vectors to upload shifts for atoms which are not calculated explicitely.
    !> Typically used in transport calculations to upload the contact atoms. Shape: [mShell, nAtom].
    real(dp), allocatable :: shellPotUpload(:,:)

  end type TPoissonInput

#:else

  type :: TPoissonInput
  end type TPoissonInput

#:endif


#:if WITH_POISSON

  ! Nr. of active Poisson-solver instances
  integer :: nInstances_ = 0

#:endif


contains

#:if WITH_POISSON

  !> Initialized Poisson solver
  !>
  !> Note: Only single instance should exist, as part of the internal state is stored in
  !>     global variables.
  !>
  subroutine TPoisson_init(this, input, env, orb, success)

    !> Instance
    type(TPoisson), intent(out) :: this

    !> Input parameters for the solver
    type(TPoissonInput), intent(inout) :: input

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> data structure with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Success of initialisation
    logical, intent(out) :: success

    integer :: nAtom

    if (nInstances_ > 0) then
      call error("Internal error: There exists already an instance of PoissonSolver")
    end if
    nInstances_ = 1

    nAtom = size(orb%nOrbAtom)
    allocate(this%shellPot_(orb%mShell, nAtom))

    call move_alloc(input%shellPotUpload, this%shellPotUpload_)

  #:if WITH_TRANSPORT
    call poiss_init_(env, input%poissonStruct, orb, input%hubbU, input%poissonInfo, input%transpar,&
        & success)
  #:else
    call poiss_init_(env, input%poissonStruct, orb, input%hubbU, input%poissonInfo, success)
  #:endif

  end subroutine TPoisson_init

#:else

  subroutine TPoisson_init()
  end subroutine TPoisson_init

#:endif


#:if WITH_POISSON

  !> Invokes the finalization of the instance.
  subroutine finalize_(this)
    type(TPoisson), intent(inout) :: this

    call poiss_destroy_()
    nInstances_ = nInstances_ - 1

  end subroutine finalize_


  !> Updates the lattice vectors.
  !>
  !> Note: The Poisson solver can not handle lattice vectors updates. This routine just checks
  !> whether the new lattice vectors are the ones used at the initialization and triggers
  !> an error, if it is not the case.
  !>
  subroutine updateLatVecs(latVecs)

    !> New lattice vectors. Shape: [3, 3]
    real(dp), intent(in) :: latVecs(:,:)

    real(dp) :: origLatVecs(3, 3)
    logical :: latticeChanged

    ! Poisson solver stores init data only on the lead node
    if (active_id) then
      call poiss_getlatvecs(origLatVecs)
      latticeChanged = any(abs(latVecs - origLatVecs) > 1e-10_dp)
    end if

    ! To warranty a proper stop, distribute result and call error() on all processes if needed
  #:if WITH_MPI
    call mpifx_bcast(global_comm, latticeChanged)
  #:endif
    if (latticeChanged) then
      call error("Internal error: Poisson solver can not handle lattice vector changes")
    end if

  end subroutine updateLatVecs


  !> Updates the coordinates
  subroutine updateCoords(coords)

    !> New coordinates. Shape: [3, nAtom]
    real(dp), intent(in) :: coords(:,:)

    call poiss_updcoords(coords)

  end subroutine updateCoords


  !> Updates the charges.
  subroutine updateCharges(this, env, qOrb, q0)

    !> Instance.
    class(TPoisson), intent(inout) :: this

    !> Environment
    type(TEnvironment), intent(inout) :: env

    !> Orbital resolved charges. Shape: [mOrb, nAtom]
    real(dp), intent(in) :: qOrb(:,:)

    !> Orbital resolved reference charges. Shape: [mOrb, nAtom].
    real(dp), intent(in), optional :: q0(:,:,:)

    call poiss_updcharges_(env, qOrb, q0)
    if (allocated(this%shellPotUpload_)) then
      this%shellPot_(:,:) = this%shellPotUpload_
    else
      this%shellPot_(:,:) = 0.0_dp
    end if
    call poiss_getshift_(env, this%shellPot_)

  end subroutine updateCharges


  !> Adds the potential calculated with the last charges.
  subroutine addPotentials(this, shellPot)

    !> Instance.
    class(TPoisson), intent(in) :: this

    !> Potential to increase by the potential corresponding to the last trasmitted charges.
    real(dp), intent(inout) :: shellPot(:,:)

    shellPot(:,:) = shellPot + this%shellPot_

  end subroutine addPotentials


  !> Returns the gradients corresponding to the last charges.
  subroutine getGradients(this, env, gradients)

    !> Instance.
    class(TPoisson), intent(inout) :: this

    !> Environment
    type(TEnvironment), intent(inout) :: env

    !> Gradients
    real(dp), intent(out) :: gradients(:,:)

    real(dp), allocatable :: dummyArray(:,:)

    allocate(dummyArray, mold=this%shellPot_)
    call poiss_getshift_(env, dummyArray, gradients)

  end subroutine getGradients


  !> Saves the internal potential of the Poisson solver.
  subroutine savePotential(env)

    !> Environment
    type(TEnvironment), intent(inout) :: env

    call poiss_savepotential(env)

  end subroutine savePotential

#:endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:if WITH_POISSON

  !> Initialise gDFTB environment and variables
#:if WITH_TRANSPORT
  subroutine poiss_init_(env, structure, orb, hubbU, poissoninfo, transpar, initinfo)
#:else
  subroutine poiss_init_(env, structure, orb, hubbU, poissoninfo, initinfo)
#:endif

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> initialisation choices for poisson solver
    Type(TPoissonStructure), intent(in) :: structure

    !> data structure with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Hubbard Us stored as (orbital, atom)
    real(dp), intent(in) :: hubbU(:,:)

    !> Solver settings
    Type(TPoissonInfo), intent(in) :: poissoninfo

  #:if WITH_TRANSPORT
    !> Transport parameters
    Type(TTransPar), intent(in) :: transpar
  #:endif

    !> Success of initialisation
    logical, intent(out) :: initinfo


    ! local variables
    integer :: iErr

    iErr = 0
    initinfo = .true.

  #:if WITH_MPI
    call poiss_mpi_init(env%mpi%globalComm)
    call poiss_mpi_split(min(poissoninfo%maxNumNodes, env%mpi%globalComm%size))
    call mpifx_barrier(env%mpi%globalComm, iErr)
  #:endif

    write(stdOut,*)
    write(stdOut,*) 'Poisson Initialisation:'
    write(stdOut,'(a,i0,a)') ' Poisson parallelized on ', numprocs, ' node(s)'
    write(stdOut,*)

    ! Directory for temporary files
    call set_scratch(poissoninfo%scratch)

    if (id0) then
      ! only use a scratch folder on the lead node
      call create_directory_(trim(scratchfolder),iErr)
    end if

    if (active_id) then
      ! processors over which the right hand side of the Poisson equation is parallelised

      call init_structure(structure%nAtom, structure%nSpecies, structure%specie0, structure%x0,&
          & structure%latVecs, structure%isperiodic)

      call init_skdata(orb%nShell, orb%angShell, hubbU, iErr)

      if (iErr.ne.0) then
        call error("I am sorry... cannot proceed. orbital shells should be in the order s,p,d")
      endif

      call init_charges()

      ! Initialise renormalization factors for grid projection

      if (iErr.ne.0) then
        call poiss_destroy_()
        initinfo = .false.
        return
      endif

      !init default values for parameters
      call init_defaults()

      call set_verbose(poissonInfo%verbose)

      call set_temperature(0.0_dp)

    #:if WITH_TRANSPORT
      !-----------------------------------------------------------------------------
      ! GP: verify if the calculation is on an open system (cluster=false) or not
      !
      !     note: in previous versions only ncont was checked but this is not
      !     sufficient as it does not take into account a contact calculation
      !     with poisson solver, where the transport block is defined
      !-----------------------------------------------------------------------------
      if (transpar%defined .and. transpar%taskUpload) then
        !init the number of contacts as in dftb+
        call set_ncont(transpar%ncont)
      else
        call set_ncont(0)
      endif
      call set_cluster(.false.)
      if (.not.transpar%defined) then
        call set_cluster(.true.)
      endif
      if (transpar%defined .and. (.not. transpar%taskUpload)) then
        call set_cluster(.true.)
      endif

      !-----------------------------------------------------------------------------+
      ! TRANSPORT PARAMETER NEEDED FOR POISSON (contact partitioning)
      !-----------------------------------------------------------------------------+
      call set_mol_indeces(transpar%idxdevice(1:2), structure%natom)
      call set_cont_indeces(transpar%contacts(1:ncont)%idxrange(1), 1)
      call set_cont_indeces(transpar%contacts(1:ncont)%idxrange(2), 2)
      call set_contdir(transpar%contacts(1:ncont)%dir)
      call set_fermi(transpar%contacts(1:ncont)%eFermi(1))
      call set_potentials(transpar%contacts(1:ncont)%potential)
      call set_builtin()
    #:else
      call set_ncont(0)
      call set_cluster(.true.)
      call set_mol_indeces([1,structure%natom], structure%natom)
    #:endif

      call set_dopoisson(poissoninfo%defined)
      call set_poissonbox(poissoninfo%poissBox)
      call set_poissongrid(poissoninfo%poissGrid)
      call set_accuracy(poissoninfo%poissAcc)

      FoundBox = poissoninfo%foundBox
      InitPot = poissoninfo%bulkBC
      ReadBulk = poissoninfo%readBulkPot
      MaxPoissIter = poissoninfo%maxPoissIter
      select case (poissoninfo%localBCType)
      case('G')
        localBC = 0
      case('C')
        localBC = 1
      case('S')
        localBC = 2
      end select
      deltaR_max = poissoninfo%maxRadAtomDens

      ! if deltaR_max > 0 is a radius cutoff, if < 0 a tolerance
      if (deltaR_max < 0.0_dp) then
        write(stdOut,*) "Atomic density tolerance: ", -deltaR_max
        deltaR_max = getAtomDensityCutoff_(-deltaR_max, uhubb)
      end if

      write(stdOut,*) "Atomic density cutoff: ", deltaR_max, "a.u."

    #:if WITH_TRANSPORT
      if (ncont /= 0 .and. poissoninfo%cutoffcheck) then
        call checkDensityCutoff_(deltaR_max, transpar%contacts(:)%length)
      end if
    #:endif

      dR_cont = poissoninfo%bufferLocBC
      bufferBox = poissoninfo%bufferBox
      SavePot = poissoninfo%savePotential
      overrideBC = poissoninfo%overrideBC
      overrBulkBC = poissoninfo%overrBulkBC
      !-----------------------------------------------------------------------------+
      ! Gate settings
      DoGate=.false.
      DoCilGate=.false.
      if (poissoninfo%gateType.eq.'C') then
        DoCilGate = .true.
      end if
      if (poissoninfo%gateType.eq.'P') then
        DoGate = .true.
      end if

      Gate = poissoninfo%gatePot
      GateLength_l = poissoninfo%gateLength_l
      GateLength_t = poissoninfo%gateLength_t

      ! Planar gate must be along y
      GateDir = poissoninfo%gatedir

      OxLength = poissoninfo%insLength
      Rmin_Gate = poissoninfo%gateRad
      Rmin_Ins = poissoninfo%insRad
      eps_r = poissoninfo%eps_r
      dr_eps = poissoninfo%dr_eps

      ! Use fixed analytical renormalization (approximate) or numerical
      fixed_renorm = .not.(poissoninfo%numericNorm)

      ! Performs parameters checks
      call check_biasdir(iErr)
      if  (iErr /= 0) then
        call error("Unable to build box for Poisson solver")
      end if
      call check_poisson_box(iErr)
      if  (iErr /= 0) then
        call error("Unable to build box for Poisson solver")
      end if
      if (any(overrideBC.ne.0)) then
        period = .false.
      end if
      call check_parameters()
      call check_localbc()
      call write_parameters()
      call check_contacts(iErr)
      if  (iErr /= 0) then
        call error("Unable to build contact potentials for Poisson solver")
      end if

      !-----------------------------------------------------------------------------+
      ! Parameters from old PAR.in not yet parsed:
      !
      ! PoissPlane
      ! DoTip,tip_atom,base_atom1,base_atom2
      !-----------------------------------------------------------------------------+

      write(stdOut,'(79(">"))')

    endif

  end subroutine poiss_init_


  subroutine create_directory_(dirName, iErr)

    character(*), intent(in) :: dirName

    integer, intent(out) :: iErr

    integer :: cstat
    character(len=255) :: cmsg

    iErr = -999
    cstat = -999
    cmsg  = "notfilled"

    call execute_command_line("mkdir "//trim(dirName), exitstat=iErr, cmdstat=cstat,&
        & cmdmsg=cmsg)
    if (iErr /= 0) then
      write (stdOut,*) 'error status of mkdir: ', iErr
      write (stdOut,*) 'command status: ', cstat
      write (stdOut,*) "command msg:  ", trim(cmsg)
    end if

  end subroutine create_directory_


  !> Release gDFTB varibles in Poisson library
  subroutine poiss_destroy_()

    if (active_id) then
      write(stdOut,'(A)')
      write(stdOut,'(A)') 'Release Poisson Memory:'
      call poiss_freepoisson()
    endif

  end subroutine poiss_destroy_


  !> Interface subroutine to call Poisson
  subroutine poiss_getshift_(env, V_L_atm,grad_V)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> potential for each shell at atom sites
    real(dp), intent(inout) :: V_L_atm(:,:)

    !> Gradient of potential
    real(dp), intent(out), optional :: grad_V(:,:)

    real(dp) :: fakegrad(3,1)

    integer :: ierr
    integer :: PoissFlag

    ! The subroutine gives the atom shifts.

    ! Poiss Flag is a control flag:
    ! 0 - potential in SCC
    ! 1 - atomic shift component of gradient

    ! This information is needed by all nodes
    PoissFlag=0
    if(present(grad_V)) then
      PoissFlag=1
    end if

    call env%globalTimer%startTimer(globalTimers%poisson)
    if (active_id) then

      select case(PoissFlag)
      case(0)
        if (verbose.gt.30) then
          write(stdOut,*)
          write(stdOut,'(80("="))')
          write(stdOut,*) '                       SOLVING POISSON EQUATION         '
          write(stdOut,'(80("="))')
        end if
        call init_PoissBox(iErr)
        if (iErr /= 0) then
          call error("Failure during initialisation of the Poisson box")
        end if
        call mudpack_drv(env, PoissFlag,V_L_atm,fakegrad)
      case(1)
        call mudpack_drv(env, PoissFlag,V_L_atm,grad_V)
      end select

      if (verbose.gt.30) then
        write(stdOut,'(80("*"))')
      end if
    end if

  #:if WITH_MPI
    call mpifx_barrier(global_comm, ierr)
    select case(PoissFlag)
      ! Note: V_L_atm and grad_V are allocated for all processes in dftb+.F90
    case(0)
      call mpifx_bcast(global_comm, V_L_atm)
    case(1)
      call mpifx_bcast(global_comm, V_L_atm)
      call mpifx_bcast(global_comm, grad_V)
    end select
    call mpifx_barrier(global_comm)
  #:endif

    call env%globalTimer%stopTimer(globalTimers%poisson)

  end subroutine poiss_getshift_


  !> Interface subroutine to overload Mulliken charges stored in libPoisson
  subroutine poiss_updcharges_(env, q, q0)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> populations
    real(dp), intent(in) :: q(:,:)

    !> reference charges
    real(dp), intent(in), optional :: q0(:,:,:)

    integer :: nsh, l, i, o, orb
    real(dp) :: Qtmp

    call env%globalTimer%startTimer(globalTimers%poisson)
    if (active_id) then

      if (size(q, dim=2).ne.natoms) then
        call error('ERROR in udpcharges size of q')
      end if

      if (present(q0)) then
        do i = 1, natoms
          nsh = lmax(izp(i))+1
          orb=0
          do l = 1, nsh
            Qtmp = 0.0_dp
            do o= 1,2*l-1
              orb = orb + 1
              ! - is for negative electron charge density
              Qtmp = Qtmp + q0(orb,i,1)-q(orb,i)
            enddo
            dQmat(l,i) = Qtmp
          enddo
        enddo
      else
        do i = 1, natoms
          nsh = lmax(izp(i))+1
          orb=0
          do l = 1, nsh
            Qtmp = 0.0_dp
            do o= 1,2*l-1
              orb = orb + 1
              ! - is for negative electron charge density
              Qtmp = Qtmp -q(orb,i)
            enddo
            dQmat(l,i) = Qtmp
          enddo
        enddo
      end if

    endif
    call env%globalTimer%stopTimer(globalTimers%poisson)

  end subroutine poiss_updcharges_


  !> Calculates the atom density cutoff from the density tolerance.
  function getAtomDensityCutoff_(denstol, uhubb) result(res)

    !> Density tolerance.
    real(dp), intent(in) :: denstol

    !> List of atomic Hubbard U values.
    real(dp), intent(in) :: uhubb(:,:)

    !> Maximal atom density cutoff.
    real(dp) :: res

    integer :: typ, nsh
    real(dp) :: tau

    res = 0.0_dp
    do typ = 1, size(uhubb,2)
      nsh = lmax(typ)+1
      tau = 3.2_dp * minval(uhubb(1:nsh,typ))
      res = max(res, -log(8.0_dp * pi / tau**3 * denstol) / tau)
    end do

  end function getAtomDensityCutoff_


#:if WITH_TRANSPORT

  !> Checks whether density cutoff fits into the PLs and stop if not.
  subroutine checkDensityCutoff_(rr, pllens)

    !> Density cutoff.
    real(dp), intent(in) :: rr

    !> Lengths of the principal layers in the contacts.
    real(dp), intent(in) :: pllens(:)

    integer :: ii

    ! GP In Poisson both the contact layers are used, which is the reason why we
    ! have a factor 2 in front of pllens
    do ii = 1, size(pllens)
      if (rr > 2.0_dp * pllens(ii) + 1e-12_dp) then
        write(stdOut,"(A,I0,A)") "!!! ERROR: Atomic density cutoff incompatible with the principle&
            & layer width in contact ", ii, "."
        write(stdOut,"(A,G10.3,A,G10.3,A)") "  (", rr, ">", pllens(ii), ")"
        call error("Either enlarge PL width in the contact or increase AtomDensityCutoff or&
            & AtomDensityTolerance.")
      end if
    end do

  end subroutine checkDensityCutoff_

#:endif

#:endif

end module dftbp_extlibs_poisson
