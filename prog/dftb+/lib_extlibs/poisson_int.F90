!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

!> Interface to libPoisson routines
module poisson_init
  use dftbp_accuracy, only : dp
  use dftbp_constants, only : pi
  use dftbp_commontypes, only : TOrbitals
  use dftbp_globalenv, only : stdOut
  use dftbp_message
#:if WITH_MPI
  use libmpifx_module
#:endif
  use poisson
  use libnegf_vars, only : TTransPar
  use system_calls, only: create_directory
  implicit none
  private

  public :: poiss_init
  public :: poiss_updcharges, poiss_getshift
  public :: poiss_destroy
  public :: poiss_updcoords
  public :: poiss_savepotential
  public :: TPoissonInfo
  public :: TPoissonStructure


  !> Geometry of the atoms for the Poisson solver
  type TPoissonStructure

    !> number of atoms in central cell
    integer :: nAtom

    !> number of species
    integer :: nSpecies

    !> type of the atoms (nAtom)
    integer, pointer :: specie0(:)

    !> atom START pos for squared H/S
    integer, pointer :: iatomstart(:)

    !> coordinates in central cell
    real(dp), pointer :: x0(:,:)

    !> total number of electrons
    real(dp) :: nel

    !> lattice vectors
    real(dp) :: latVecs(3,3)

    !> electron temperature
    real(dp) :: tempElec

    !> tells whether the system is periodic
    logical :: isperiodic

  end type TPoissonStructure


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

    !> Recompute poisson after D.M.
    logical :: solveTwice  = .false.

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


contains

  !> Initialise gDFTB environment and variables
#:if WITH_MPI
  subroutine poiss_init(structure, orb, hubbU, poissoninfo, transpar, mpicomm, initinfo)
#:else
  subroutine poiss_init(structure, orb, hubbU, poissoninfo, transpar, initinfo)
#:endif
    !> initialisation choices for poisson solver
    Type(TPoissonStructure), intent(in) :: structure

    !> data structure with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Hubbard Us stored as (orbital, atom)
    real(dp), intent(in) :: hubbU(:,:)

    !> Solver settings
    Type(TPoissonInfo), intent(in) :: poissoninfo

    !> Transport parameters
    Type(TTransPar), intent(in) :: transpar
#:if WITH_MPI
    !> MPI details
    Type(mpifx_comm), intent(in) :: mpicomm
#:endif
    !> Success of initialisation
    logical, intent(out) :: initinfo

    ! local variables
    integer :: i, iErr

    iErr = 0
    initinfo = .true.

  #:if WITH_MPI
    call poiss_mpi_init(mpicomm)
    call poiss_mpi_split(min(poissoninfo%maxNumNodes, mpicomm%size))
    call mpifx_barrier(mpicomm, iErr)
  !#:else
    !call error("The Poisson solver currently requires MPI parallelism to be enabled")
  #:endif

    write(stdOut,*)
    write(stdOut,*) 'Poisson Initialisation:'
    write(stdOut,'(a,i0,a)') ' Poisson parallelized on ',numprocs,' node(s)'
    write(stdOut,*)

    ! notify solver of standard out unit
    call set_stdout(stdOut)

    ! Directory for temporary files
    call set_scratch(poissoninfo%scratch)

    if (id0) then
      ! only use a scratch folder on the master node
      call create_directory(trim(scratchfolder),iErr)
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
        call poiss_destroy()
        initinfo = .false.
        return
      endif

      !init default values for parameters
      call init_defaults()

      call set_verbose(poissonInfo%verbose)

      call set_temperature(0.0_dp)

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
        deltaR_max = getAtomDensityCutoff(-deltaR_max, uhubb)
      end if

      write(stdOut,*) "Atomic density cutoff: ", deltaR_max, "a.u."

      if (ncont /= 0 .and. poissoninfo%cutoffcheck) then
        call checkDensityCutoff(deltaR_max, transpar%contacts(:)%length)
      end if

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

  end subroutine poiss_init


  !> Release gDFTB varibles in poisson library
  subroutine poiss_destroy()

    if (active_id) then
      write(stdOut,'(A)')
      write(stdOut,'(A)') 'Release Poisson Memory:'
      call poiss_freepoisson()
    endif

  end subroutine poiss_destroy


  !> Interface subroutine to call Poisson
  subroutine poiss_getshift(V_L_atm,grad_V)

    !> potential for each shell at atom sites
    real(dp), intent(inout) :: V_L_atm(:,:)

    !> Gradient of potential
    real(dp), intent(out), optional :: grad_V(:,:)

    real(dp) :: fakegrad(3,1)

    integer :: ndim, ierr, array_size
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
        call mudpack_drv(PoissFlag,V_L_atm,fakegrad)
      case(1)
        call mudpack_drv(PoissFlag,V_L_atm,grad_V)
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

  end subroutine poiss_getshift


  !> Interface subroutine to overload Mulliken charges stored in libPoisson
  subroutine poiss_updcharges(q,q0)

    !> populations
    real(dp), intent(in) :: q(:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:)

    integer :: nsh, l, i, o, orb
    real(dp) :: Qtmp

    if (active_id) then

      if (size(q, dim=2).ne.natoms) then
        call error('ERROR in udpcharges size of q')
      end if

      do i = 1, natoms
        nsh = lmax(izp(i))+1
        orb=0
        do l = 1, nsh
          Qtmp = 0.0_dp
          do o= 1,2*l-1
            orb = orb + 1
            ! - is for negative electron charge density
            Qtmp = Qtmp + q0(orb,i)-q(orb,i)
          enddo
          dQmat(l,i) = Qtmp
        enddo
      enddo

    endif

  end subroutine poiss_updcharges


  !> Calculates the atom density cutoff from the density tolerance.
  function getAtomDensityCutoff(denstol, uhubb) result(res)

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

  end function getAtomDensityCutoff


  !> Checks whether density cutoff fits into the PLs and stop if not.
  subroutine checkDensityCutoff(rr, pllens)

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

  end subroutine checkDensityCutoff


end module poisson_init
