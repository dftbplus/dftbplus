! **************************************************************************
! *                INTERFACE of gDFTB for DFTB+
! *
! *  call poiss_init(structure,skdata,gdftbinfo)
! *  
! *  call poiss_supercell(skdata)
! *
! *  call poiss_getshift(V_atm,grad_V)
! *
! *  call poiss_updcharges(qmulli_orb)
! *
! *  call poiss_updcoords(x0)
! *
! *  call poiss_density(miter,HH,SS,DensMat,EnMat)
! *
! *  call poiss_current(icycle,HH,SS)
! *
! *  call poiss_destroy()
! *  ----------------------------------------------------------------------
! *  To do:
! *
! *  1. set structure variables from DFTB+:  DONE
! *
! *  2. set matrix conversion utilities in module mat_conv: 
! *
! *  3. change structure of poiss_density and gdftb_current
! *     to accept in input DFTB matrix fomats and convert to gDFTB format
! *
! * ---------------------------------------------------------------------- 
! * Known Bugs & Problems:
! * 
! *
! *
! * ---------------------------------------------------------------------- 
#:include "common.fypp"


module poisson_int

  use gprecision
  use gconstants, only : Kb, pi
  use gallocation
  use parameters
  use parcheck
  use structure, only : natoms, nbeweg, x, boxsiz, period, &
                        period_dir, izp, dQmat, ntypes, &
                        inversebox, gamma_summind, buildsupercell, &
                        lmax, uhubb
  use poisson, only : init_poissbox,mudpack_drv
  use mpi_poisson
  use poisson_vars, only : TPoissonInfo, TPoissonStructure, TSKdata
  use libnegf_vars, only : TTransPar
  use CommonTypes, only : TOrbitals
#:if WITH_MPI
  use libmpifx_module 
#:endif
  use gclock
  use system_calls, only: create_directory
  implicit none
  private

  public :: poiss_init,  poiss_getshift 
  public :: poiss_destroy, poiss_supercell, poiss_savepotential
  public :: poiss_freepoisson, poiss_updcoords, poiss_updcharges
  !public :: poiss_mulliken 
  ! ---------------------------------------------------------------------------
  contains
    
!------------------------------------------------------------------------------
! Init gDFTB environment and varibles
!------------------------------------------------------------------------------
 subroutine poiss_init(structure, skdata, poissoninfo, transpar,&
        & mpicomm, initinfo)
      
  implicit none

  Type(TPoissonStructure), intent(IN) :: structure
  Type(TSKdata), intent(IN) :: skdata
  Type(TPoissonInfo) :: poissoninfo
  Type(TTransPar) :: transpar
  Type(mpifx_comm), intent(in) :: mpicomm
  logical, intent(OUT) :: initinfo

  ! local variables
  integer :: i,error  

  error = 0 
  initinfo = .true.

#:if WITH_MPI
  if (poissoninfo%maxNumNodes>mpicomm%size) then
    poissoninfo%maxNumNodes = mpicomm%size
  end if   
  call poiss_mpi_init(mpicomm)
  call poiss_mpi_split(poissoninfo%maxNumNodes)
  call mpi_barrier(mpicomm, error)
  print*,'poisson_int:',global_comm%rank, id, id0, active_id
#:endif

  if (id0) then
    write(*,*)
    write(*,*) 'Poisson Initializations'
    write(*,'(a,i0,a)') 'Poisson parallelized on ',numprocs,' node(s)'
  end if 

  allocate(character(len=len(poissoninfo%scratch))::scratchfolder)
  scratchfolder = trim(poissoninfo%scratch)
  if (id0) then
    call create_directory(trim(scratchfolder),error)
  end if
  
  if (active_id) then

    call init_structure(structure)
    
    call init_skdata(skdata,error)
 
    call init_charges()
 
    !! Initialize renormalization factors for grid projection
 
    if(error.ne.0) then 
        call poiss_destroy()
        initinfo = .false.; return
    endif       
 
    call init_defaults()            !init default values for parameters
 
 
    verbose = poissonInfo%verbose
 
    Temp = 0.0_dp
 
    !-----------------------------------------------------------------------------
    ! GP: verify if the calculation is on an open system (cluster=false) or not
    !
    !     note: in previous versions only ncont was checked but this is not 
    !     sufficient as it does not take into account a contact calculation 
    !     with poisson solver, where the transport block is defined
    !-----------------------------------------------------------------------------  
    if (transpar%defined .and. transpar%taskUpload) then
      ncont = transpar%ncont         !init the number of contacts as in dftb+
    endif
    cluster = .false.
    if (.not.transpar%defined) then
      cluster = .true.
    endif
    if (transpar%defined .and. (.not. transpar%taskUpload)) then
      cluster = .true.
    endif
 
    !-----------------------------------------------------------------------------+
    ! TRANSPORT PARAMETER NEEDED FOR POISSON (contact partitioning)
    !-----------------------------------------------------------------------------+
    
    if (cluster) then
      ! TODO: GP Does it make sense to have a specific direction with no contacts? I commented it out
      ! because I don't know where this is initialized. Idem for EFermi, I don't think is needed
      !if( associated( transpar%cdir )) then
      !   contdir(1) = transportinfo%cdir(1) 
      !else
      contdir(1) = 0
      !end if    
      !if( associated(transportinfo%eFermi) ) then
      !   EFermi(1) = transportinfo%eFermi(1) 
      !else
      EFermi(1) = 0.d0
      !endif
      ! Note: moved here iatm initialization for cluster case, consistent with not 
      ! cluster case
      iatm(1)=1 
      iatm(2)=natoms
    else
      iatm(1:2) = transpar%idxdevice(1:2)
      iatc(1,1:ncont) = transpar%contacts(1:ncont)%idxrange(1)
      iatc(2,1:ncont) = transpar%contacts(1:ncont)%idxrange(2)
      iatc(3,1:ncont) = 0
      contdir(1:ncont) = transpar%contacts(1:ncont)%dir
      mu(1:ncont) = transpar%contacts(1:ncont)%potential
      ! Note: the poisson machinery is not compatible with built in 
      ! potential with colinear spin calculations
      EFermi(1:ncont) = transpar%contacts(1:ncont)%eFermi(1)
    endif
 
    if (.not.cluster) then
      do i = 1,ncont
        mu(i) = mu(i) + EFermi(i) -  minval(EFermi(1:ncont))
      enddo
    end if
 
    DoPoisson = poissoninfo%defined
    PoissBox(1,1) = poissoninfo%poissBox(1)
    PoissBox(2,2) = poissoninfo%poissBox(2)
    PoissBox(3,3) = poissoninfo%poissBox(3)
    dmin(1:3) = poissoninfo%poissGrid(1:3)
    FoundBox = poissoninfo%foundBox 
    PoissAcc = poissoninfo%poissAcc
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
      if (id0) write(*,*) "Atomic density tolerance: ", -deltaR_max
      deltaR_max = getAtomDensityCutoff(-deltaR_max, uhubb)
    end if
    if (id0)  write(*,*) "Atomic density cutoff: ", deltaR_max,"a.u."
 
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
    if (poissoninfo%gateType.eq.'C')   DoCilGate = .true.
    if (poissoninfo%gateType.eq.'P')   DoGate = .true.
    
    Gate = poissoninfo%gatePot
    GateLength_l = poissoninfo%gateLength_l
    GateLength_t = poissoninfo%gateLength_t
    GateDir = poissoninfo%gatedir   ! Planar gate must be along y
    OxLength = poissoninfo%insLength
    Rmin_Gate = poissoninfo%gateRad
    Rmin_Ins = poissoninfo%insRad  
    eps_r = poissoninfo%eps_r
    dr_eps = poissoninfo%dr_eps
 
    ! Use fixed analytical renormalization (approximate) or new numerical one??
    fixed_renorm = .not.(poissoninfo%exactRenorm)
    ! Performs parameters checks    
    !!if(id0.and.verbose.gt.40) call echo_init()
 
    call check_biasdir()
    call check_poisson_box()
    if (any(overrideBC.ne.0)) period = .false.
    call check_parameters() 
    call check_localbc()
    call write_parameters()
    call check_contacts()
 
    !-----------------------------------------------------------------------------+
    ! Special hacking to allow one-dimensional 
    ! periodic calculations using Poisson
    ! ('cluster' means 'without contacts')
    !if(period.and.cluster.and.DoPoisson) then
    !   write(*,'(A)') 'Supercell vectors for Hamiltonian:'
    !   write(*,*) boxsiz(:,:)
    !   if(contdir(1).eq.0) then
    !     period_dir(:) = .true.      
    !   else
    !     period_dir(:) = .false.
    !     period_dir(contdir(1)) = .true.
    !     write(*,'(A)') 'The system has 1D periodicity for Poisson'
    !   endif  
    !end if
 
 
    !-----------------------------------------------------------------------------+
    ! Parameters from old PAR.IN not yet parsed:
    !
    ! PoissPlane
    ! DoTip,tip_atom,base_atom1,base_atom2
    !-----------------------------------------------------------------------------+
 
    if (id0) write(*,'(79(">"))')
  
  endif 

end subroutine poiss_init

!--------------------------------------------------------------------------
!----------------------------------------------------------------------------
! Release gDFTB varibles
!----------------------------------------------------------------------------
subroutine poiss_destroy()

  if (active_id) then
    if (id0) write(*,'(A)')
    if (id0) write(*,'(A)') 'Release Poisson Memory:'
    call poiss_freepoisson()
    !if(allocated(ind)) call log_gdeallocate(ind) 
    if(allocated(x)) call log_gdeallocate(x)
    if(allocated(izp)) call log_gdeallocate(izp) 
    if(allocated(dQmat)) call log_gdeallocate(dQmat)
    if(allocated(uhubb)) call log_gdeallocate(uhubb)
    if(allocated(lmax)) call log_gdeallocate(lmax)
    if(allocated(renorm)) call log_gdeallocate(renorm)
    if (id0) then
      call writePeakInfo(6)
      call writeMemInfo(6)
    endif
  endif 
  
  deallocate(scratchfolder)

end subroutine poiss_destroy

!----------------------------------------------------------------------------
! INTERFACE subroutine to call Poisson 
!----------------------------------------------------------------------------
subroutine poiss_getshift(V_L_atm,grad_V)

  real(kind=dp), DIMENSION(:,:), intent(inout) :: V_L_atm
  real(kind=dp), DIMENSION(:,:), optional :: grad_V

  real(kind=dp), DIMENSION(3,1) :: fakegrad

  integer :: ndim, ierr, array_size
  integer :: PoissFlag

  ! The subroutine gives the atom shifts.  

  ! Poiss Flag is a control flag:
  ! 0 - potential in SCC
  ! 1 - atomic shift component of gradient

  !! This iformation is needed by all nodes
  PoissFlag=0  
  if(present(grad_V)) PoissFlag=1

  if (active_id) then

    ! Sets a local copy of V_orb for shift calculations. 
    ! This will change soon by passing directly V_orb
    !ndim = ind(natoms+1)
    !write(*,*) 'Calling Poisson Solver'
    !write(*,*) 'cluster box:',cluster
    !write(*,*) 'periodic box:',period
    !write(*,*) 'periodic dir:',period_dir
    !write(*,*) 'contact dir:',contdir(1:ncont)
    !write(*,*) 'bias dir:',BiasDir 
    !write(*,*) 'Do Gate:',DoGate, GateDir
    !write(*,*) 'Do Cil Gate:',DoCilGate
    !write(*,*) 'Do Tip:',DoTip 
    !write(*,*) 'Use bulk potential:',InitPot
    !write(*,*) 'Read old bulk pot:',ReadBulk
    !write(*,*) 'Local BC:',LocalBC

    ! .....................................................
    select case(PoissFlag)
    case(0)
      if (id0.and.verbose.gt.30) then
        write(*,*)
        write(*,'(80("="))')
        write(*,*) '                       SOLVING POISSON EQUATION         '
        write(*,'(80("="))') 
      endif
      call init_PoissBox
      call mudpack_drv(PoissFlag,V_L_atm,fakegrad)
    case(1)
      call mudpack_drv(PoissFlag,V_L_atm,grad_V)
    end select

    if (id0.and.verbose.gt.30) write(*,'(80("*"))')

  endif

#:if WITH_MPI
  call mpifx_barrier(global_comm, ierr)
  select case(PoissFlag)
    !! Note: V_L_atm and grad_V are allocated for all processes in dftb+.F90
  case(0)
    call mpifx_bcast(global_comm, V_L_atm)
  case(1) 
    call mpifx_bcast(global_comm, V_L_atm)
    call mpifx_bcast(global_comm, grad_V)
  end select
  call mpifx_barrier(global_comm)
#:endif

end subroutine poiss_getshift


! -----------------------------------------------------------------------------
subroutine poiss_savepotential()
  
  real(kind=dp), DIMENSION(3,1) :: fakegrad
  real(kind=dp), DIMENSION(1,1) :: fakeshift
  integer :: PoissFlag, ndim

  if (active_id) then

    if(.not.SavePot) return
 
    if (id0) write(*,"('>> Saving Poisson output in potential.dat')")
 
    PoissFlag=2
    !ndim = ind(natoms+1)-1
    
    call  mudpack_drv(PoissFlag,fakeshift,fakegrad)
  
  endif

end subroutine poiss_savepotential


! -----------------------------------------------------------------------------
subroutine poiss_freepoisson()

  real(kind=dp), DIMENSION(3,1) :: fakegrad
  real(kind=dp), DIMENSION(1,1) :: fakeshift
  integer :: PoissFlag, ndim 

  if (active_id) then
    PoissFlag=3
    call  mudpack_drv(PoissFlag,fakeshift,fakegrad)  
  endif

end subroutine poiss_freepoisson

! -----------------------------------------------------------------------------
!  FILL UP Structure Parameters
! -----------------------------------------------------------------------------
subroutine init_structure(struct)

  Type(TPoissonStructure), intent(IN) :: struct
  !  From Tstructure:
  !  integer               :: nAtom          ! number of Atoms in central cell 
  !  integer               :: nSpecies       ! number of Species
  !  integer, pointer      :: specie0(:)     ! type of each Atoms (nAtom)
  !  integer, pointer      :: iatomstart(:)  ! atom START pos for squared H/S
  !  integer               :: nImages(3)     ! number of image cells
  !  real(dp), pointer     :: x0(:,:)        ! coordinates in central cell
  !  real(dp)              :: latVecs(3,3)   ! lattice vectors
  !  real(dp)              :: tempElec       ! electron temperature 
  !  real(dp)              :: q0(:)          ! neutral atom charge (nSpecies)   
  integer :: i,j
  real(dp) :: side(3)
  integer :: dir(3)

  if (active_id) then

    natoms=struct%nAtom
    nbeweg=0
    ntypes=struct%nSpecies
     
    dir = (/ 1, 2, 3 /)
    boxsiz(1:3,1:3) = 0_dp
    ! the latVec structure is :
    !   
    !  [e_x1, e_x2, e_x3]
    !  [e_y1, e_y2, e_y3] 
    !  [e_z1, e_z2, e_z3]
 
    ! computes the three sides
    ! side(:) = sqrt(sum(struct%latVecs,dim=1)**2)
 
    boxsiz(:,:) = transpose(struct%latVecs(:,:))
 
    do i = 1, 3 
       if (abs(boxsiz(i,i))<1.0d-4) then
         do j = i+1, 3      
           if (abs(boxsiz(j,i))>1.d-4) then
             side(:)=boxsiz(i,:)
             boxsiz(i,:)=boxsiz(j,:)
             boxsiz(j,:)=side(:)
           endif  
         enddo
       endif
       if(boxsiz(i,i)<0.d0) boxsiz(i,:)=-boxsiz(i,:)
    enddo
    !do i = 1, 3
    !   mask(:) = (abs(abs(struct%latVecs(:,i)) - side(i)) < 1.0d-8)
    !   j = minval(dir,mask)
    !   boxsiz(j,j) = side(i) 
    !enddo
 
    call inversebox()
  
    period=struct%isperiodic
 
    !call log_gallocate(ind,natoms+1)
    !ind(1:natoms+1)=struct%iatomstart(1:natoms+1)-1
  
    call log_gallocate(x,3,natoms)
    x(1:3,1:natoms)=struct%x0(1:3,1:natoms)
 
    call log_gallocate(izp,natoms)   
    izp(1:natoms)=struct%specie0(1:natoms)
 
 
    telec=struct%tempElec
 
    !open(99, file = 'str.gen') 
    !write(99,*) natoms  
    !do i = 1, natoms
    !   write(99,'(i4,i3,f14.8,f14.8,f14.8)') i,izp(i),x(:,i)*0.529177d0
    !end do
    !close(99)

  endif 
  
end subroutine init_structure

!------------------------------------------------------------------------------
subroutine init_charges()
  integer :: nshells

  if (active_id) then
    nshells = maxval(lmax)+1 
    call log_gallocate(dQmat,nshells,natoms)
  endif

end subroutine init_charges


!------------------------------------------------------------------------------
subroutine init_skdata(skdata,err)

  Type(TSKdata), intent(IN) :: skdata
  integer :: i,j,err,nshells

  if (active_id) then

    ! set maximum angular momentum per specie
    call log_gallocate(lmax,ntypes)

    !checks that all atoms have shells in s,p,d sequence
    err=0
    do i=1,ntypes
       do j=1,skdata%orb%nShell(i)
          if (skdata%orb%angShell(j,i).ne.j-1) err=1  
       enddo 
       lmax(i)=skdata%orb%angShell(skdata%orb%nShell(i),i)
    enddo
    if (err.ne.0 .and. id0) then
       write(*,*) 'I am sorry... cannot preceed with gdftb'
       write(*,*) 'orbital shells should be in the order s,p,d'
       return
    endif
 
    ! set Hubbard parameters
    nshells = maxval(lmax)+1
 
    call log_gallocate(uhubb,nshells,ntypes)
 
    do i = 1,ntypes
      do j = 1,nshells
         uhubb(j,i) = skdata%hubbU(j,i)
      enddo   
    end do

  endif
  
end subroutine init_skdata

!------------------------------------------------------------------------------
subroutine echo_init()

  integer :: i

  if (active_id) then

    write(*,*) 'natoms=',natoms
    write(*,*) 'ntypes=',ntypes
    write(*,*) 'is Periodic=',period

    !write(*,*) 'coordinates='
    !do i=1,natoms
    !   write(*,'(i5,i3,3(f9.4))') i,izp(i),x(1,i),x(2,i),x(3,i)
    !enddo

    !write(*,*) 'ind='
    !write(*,*) ind(1:natoms+1)
    write(*,*) '  qzero,   uhubb,   lmax' 

    do i=1,ntypes  
       write(*,'(f9.4,f9.4,i6)') uhubb(1:2,i),lmax(i)
    enddo

  endif

end subroutine echo_init

!------------------------------------------------------------------------------
! Utility subroutine to build gDFTB supercell
!------------------------------------------------------------------------------
subroutine poiss_supercell(skdata)

  Type(TSKdata), intent(IN) :: skdata
  !  From TSKdata:  
  !  real(dp)              :: mCutoff        !* longest pair interaction 
  real(dp) :: slkcutoff
  
  if (active_id) then
    slkcutoff=skdata%mCutoff
    call gamma_summind(slkcutoff)
                         
    call buildsupercell() !(computes ss_natoms and allocates inside ss_x, ss_izp)
                          ! period_dir should have been computed (check_biasdir)
  endif

end subroutine poiss_supercell
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! INTERFACE subroutine to update the coordinates
!------------------------------------------------------------------------------
subroutine poiss_updcoords(x0)

  real(dp), dimension(:,:), intent(IN) :: x0

  if (active_id) then
    x(1:3,1:natoms)=x0(1:3,1:natoms)
    do_renorm = .true.
  endif
  ! Check about periodic images !!
  
end subroutine poiss_updcoords

!----------------------------------------------------------------
! INTERFACE subroutine to overlaod Mulliken charges 
!----------------------------------------------------------------
subroutine poiss_updcharges(q,q0)
  real(dp), dimension(:,:), intent(IN) :: q, q0

  integer :: nsh, l, i, o, orb 
  real(dp) :: Qtmp

  if (active_id) then

    if (size(q, dim=2).ne.natoms) stop 'ERROR in udpcharges size of q'
 
    do i = 1, natoms
      nsh = lmax(izp(i))+1
      orb=0
      do l = 1, nsh 
      Qtmp = 0.d0
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
  !! \param denstol  Density tolerance.
  !! \param uhubb  List of atomic Hubbard U values.
  !! \return  Maximal atom density cutoff.
  function getAtomDensityCutoff(denstol, uhubb) result(res)
    real(dp), intent(in) :: denstol, uhubb(:,:)
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
  !! \param rr  Density cutoff.
  !! \param pllens  Lengths of the principal layers in the contacts.
  subroutine checkDensityCutoff(rr, pllens)
    real(dp), intent(in) :: rr, pllens(:)

    integer :: ii

    ! GP In Poisson both the contact layers are used, which is the reason why we
    ! have a factor 2 in fron of pllens
    do ii = 1, size(pllens)
      if (rr > 2.0_dp * pllens(ii) + 1e-12_dp) then
        if (id0) then
          write(*,"(A,I0,A)") "!!! ERROR: Atomic density cutoff incompatible with&
              & PL width in conctact ", ii, "."
          write(*,"(A,G10.3,A,G10.3,A)") "  (", rr, ">", pllens(ii), ")"
          write(*,"(A)") "Either enlarge PL width in the contact or increase&
              & AtomDensityCutoff or AtomDensityTolerance."
        end if    
        stop
      end if
    end do

  end subroutine checkDensityCutoff
  
  

! ---------------------------------------------------------------
!subroutine poiss_mulliken(csrDens, csrOver, kweight, q_atoms)
!   type(z_CSR) :: csrDens, csrOver
!   real(dp) :: kweight
!   real(dp), dimension(:) :: q_atoms     
!
!   real(dp), dimension(:), allocatable :: qmulli
!
!   call log_gallocate(qmulli,csrOver%nrow)
!
!   qmulli=0.d0
!
!   call mulliken(csrDens,csrOver,qmulli)
!
!   q_atoms = q_atoms + qmat*kweight
!
!   call log_gdeallocate(qmulli)
!
!end subroutine poiss_mulliken     

end module poisson_int
