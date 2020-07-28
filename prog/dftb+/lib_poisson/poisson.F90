!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program * 
!  under the LGPL licence.                                                *
!**************************************************************************
#:include "common.fypp"
#:include "error.fypp"

module poisson

  use dftbp_constants, only : pi, hartree__eV, Bohr__AA
  use gallocation
  use parameters
  use structure
  use parcheck                
  use gewald
  use bulkpot
  use fancybc
  use mpi_poisson
  use dftbp_globalenv, only : stdOut
  use dftbp_accuracy, only : lc, dp
  use dftbp_environment, only : TEnvironment, globalTimers
  implicit none
  private

  !> Verbosity threashold
  integer, parameter :: VBT=30

  ! from parameters
  public :: MAXNCONT
  public :: verbose, biasdir, gatedir, contdir, ncont, ni, nf
  public :: iatc, iatm, ncdim, mbound_end, maxiter, localBC, poissBC
  public :: overrideBC, overrBulkBC, maxpoissiter
  public :: Temp, telec, deltaR_max, LmbMax, gate
  public :: GateLength_l, GateLength_t, OxLength, Efermi
  public :: bias_dEf, Rmin_Ins, Rmin_Gate, dr_eps
  public :: eps_r, cntr_gate, tip_atom, base_atom1, base_atom2, tipbias
  public :: DOS, delta, racc, PoissAcc, dmin, PoissBox, PoissBounds
  public :: PoissPlane,mu,cntr_cont,R_cont,dR_cont,x0,y0,z0,bufferBox
  public :: etb ,cluster ,SavePot ,SaveNNList,SaveHS,FictCont,FoundBox,DoPoisson
  public :: Readold,InitPot,DoGate,DoCilGate,DoTip,ReadBulk,mixed
  public :: do_renorm, fixed_renorm
  public :: scratchfolder
  public :: init_defaults
  public :: set_verbose, set_scratch, set_contdir, set_cluster
  public :: set_fermi, set_potentials, set_builtin
  public :: set_temperature, set_ncont, set_mol_indeces, set_cont_indeces
  public :: set_dopoisson, set_poissonbox, set_poissongrid, set_accuracy  
  ! from structure
  public :: init_structure, init_charges, init_skdata
  public :: natoms, dQmat, lmax, uhubb, izp, period
  ! from parcheck
  public :: check_poisson_box, check_contacts, check_localbc
  public :: check_parameters, write_parameters, check_biasdir
  ! from mpi_poisson
#:if WITH_MPI 
  public :: poiss_mpi_init, poiss_mpi_split, global_comm, poiss_comm
#:endif
  public :: id0, active_id, numprocs
  ! from io


  public :: init_poissbox, mudpack_drv, save_pot, rho, n_alpha
  public :: poiss_updcoords, poiss_savepotential, poiss_freepoisson

 contains

 !------------------------------------------------------------------------------
 subroutine poiss_freepoisson(env)

   !> Environment settings
   type(TEnvironment), intent(inout) :: env

   real(kind=dp), DIMENSION(3,1) :: fakegrad
   real(kind=dp), DIMENSION(1,1) :: fakeshift
   integer :: PoissFlag

   if (active_id) then
     PoissFlag=3
     call  mudpack_drv(env, PoissFlag,fakeshift,fakegrad)
   endif
 
   if(allocated(x)) call log_gdeallocate(x)
   if(allocated(izp)) call log_gdeallocate(izp) 
   if(allocated(dQmat)) call log_gdeallocate(dQmat)
   if(allocated(uhubb)) call log_gdeallocate(uhubb)
   if(allocated(lmax)) call log_gdeallocate(lmax)
   if(allocated(renorm)) call log_gdeallocate(renorm)
   if (id0) then
     call writePoissPeakInfo(stdOut)
     call writePoissMemInfo(stdOut)
   endif
   
 end subroutine poiss_freepoisson

 ! -----------------------------------------------------------------------------
 subroutine poiss_savepotential(env)

   !> Environment settings
   type(TEnvironment), intent(inout) :: env

   real(kind=dp), DIMENSION(3,1) :: fakegrad
   real(kind=dp), DIMENSION(1,1) :: fakeshift
   integer :: PoissFlag

   call env%globalTimer%startTimer(globalTimers%poisson)

   if (active_id) then
 
     if(.not.SavePot) return
  
     write(stdOut,"('>> Saving Poisson output in potential.dat')")
  
     PoissFlag=2
     
     call  mudpack_drv(env, PoissFlag,fakeshift,fakegrad)
   
   endif

   call env%globalTimer%stopTimer(globalTimers%poisson)

 end subroutine poiss_savepotential
  
 !------------------------------------------------------------------------------
 subroutine poiss_updcoords(x0)
   real(dp), dimension(:,:), intent(in) :: x0

   if (active_id) then
     x = x0
     do_renorm = .true.
   endif

 end subroutine poiss_updcoords

 ! -----------------------------------------------------------------------------
 subroutine init_poissbox(iErr)

  !> Error code, 0 if no problems
  integer, intent(out), optional :: iErr

  integer :: i,m,s,f
  real(kind=dp) :: bound(MAXNCONT) 
  real(kind=dp) :: tmp,Lx, xmax, xmin 
  integer :: tmpdir(3)
  character, parameter :: dir(3) = ['x', 'y', 'z']

  if (present(iErr)) then
    iErr = 0
  end if

  !*******************************************************************************
  ! 1. Set Poisson Box 
  !*******************************************************************************

  !  e.g.: x direction
  !                     +++++++++++++++++++++++
  !            ---------+---------------------+---------
  ! <--- contdir(1)=-1 |+C1|    MOLECULE   |C2+| contdir(2)=1 --->   
  !            ---------+---------------------+---------
  !                     +++++++++++++++++++++++
  !  ------------------------------------------------------------------>x    
  
  ! The first check is used to decide the position of the molecular atoms
  ! with respect to contact 1. 
  ! The idea is to define the Poisson box from the atoms closer to the 
  ! molecule, e.g. :
  !           +++++++++++++++++++
  !  ---------+-----------------+---------
  !     C1   |+|    MOLECULE   |+|   C2 
  !  ---------+-----------------+---------
  !           +++++++++++++++++++
  ! 
  ! The box is cut in the middle between the contact and molecule atoms 
  ! closer to the boundary. Note that in this way it's easier to  reload
  ! chunks of periodic potentials from bulk calculations (GP: patch from Stan)
  ! If a buffer length for the box is defined, then the box boundary are shifted
  ! inside the contact region, according to the buffer parameter 
  ! (can be used to reproduce old default settings)

  do m = 1, ncont
     f=abs(contdir(m))
     if (contdir(m).gt.0) then
       xmax = maxval(x(f,iatm(1):iatm(2)))
       xmin = minval(x(f,iatc(3,m):iatc(2,m)))
       if (xmin > xmax) then
         bound(m) = 0.5_dp * (xmax + xmin) + bufferBox
       else
         @:FORMATTED_ERROR_HANDLING(iErr, -1, "(A,I0)",&
             & "Device and contact atoms overlap at contact", m)
       end if
     else                          
       xmin = minval(x(f,iatm(1):iatm(2)))
       xmax = maxval(x(f,iatc(3,m):iatc(2,m)))
       if (xmin > xmax) then
         bound(m) = 0.5_dp * (xmax + xmin) - bufferBox
       else
         @:FORMATTED_ERROR_HANDLING(iErr, -2, "(A,I0)",&
             & "Device and contact atoms overlap at contact", m)
       end if
     end if
  end do

  tmpdir(:) = 0
  !check if there are two or more contacts in the same direction
  do m=1,ncont
     f=abs(contdir(m))
     do s=m+1,ncont
        if( f.eq.abs(contdir(s)) .and. contdir(m).eq.-contdir(s) ) then                 
           PoissBox(f,f)= abs(bound(s) -  bound(m)) !set PoissonBox to the direction of m and s
           
           PoissBounds(f,2) = max(bound(s),bound(m))
           PoissBounds(f,1) = min(bound(s),bound(m)) 
           
           tmpdir(f)=1
        endif
        if (contdir(m).eq.contdir(s).and.bound(s).ne.bound(m)) then
          @:ERROR_HANDLING(iErr, -3, 'Contacts in the same direction must be aligned')
        endif
     enddo
     ! Adjust PoissonBox if there are no facing contacts
     if(tmpdir(f).eq.0) then

        if (contdir(m).gt.0) then     
           !fparm((f*2)-1) = bound(m) - PoissBox(f,f)
           PoissBounds(f,1) = min( bound(m) - PoissBox(f,f), &
                              minval(x(f,1:iatm(2)))-2.d0*deltaR_max )
           PoissBounds(f,2) = bound(m)
        else                          
           !fparm((f*2))   = bound(m) + PoissBox(f,f)
           PoissBounds(f,2) = max( bound(m) + PoissBox(f,f), &
                              maxval(x(f,1:iatm(2)))+2.d0*deltaR_max )
           PoissBounds(f,1) = bound(m)
        end if
        PoissBox(f,f) = PoissBounds(f,2)-PoissBounds(f,1)
        tmpdir(f)=1
     endif
     
  enddo
  
  !------------------------------------------------------
  ! Range of the remaining spatial variables 
  !------------------------------------------------------
  do i = 1, 3
     if (period_dir(i)) then
        PoissBox(i,i) = boxsiz(i,i)
     endif
  end do

  ! set bounds (PoissBox)
  ! Set PoissonBox in the directions where there are no facing contacts
  do i = 1, 3
     
     f=iatm(2)
     if (tmpdir(i) .eq. 0) then
        
        tmp=(maxval(x(i,1:f))-minval(x(i,1:f))+2.d0*deltaR_max)

        if(any(localBC.gt.0)) then
           tmp = tmp - 2.d0*deltaR_max + 2.d0*maxval(dR_cont)
        endif
        
        if (.not.period_dir(i).and.PoissBox(i,i).le.tmp) PoissBox(i,i) = tmp
           
        Lx= maxval(x(i,1:f))+minval(x(i,1:f))
        
        PoissBounds(i,2) = ( Lx+PoissBox(i,i) )/2.d0 
        PoissBounds(i,1) = ( Lx-PoissBox(i,i) )/2.d0 

     end if
     
  end do

  !-------------------------------
  ! Checking Poisson Box
  !---- ---------------------------
  do i=1,3  
    if(PoissBox(i,i) .le. 0.0_dp) then
      @:FORMATTED_ERROR_HANDLING(iErr, -4, '(A,A)', 'PoissBox negative along ', dir(i))
    end if
  enddo


  !-------------------------------
  ! Checking for gate dimension
  !-------------------------------
  if (DoGate) then
     biasdir = abs(contdir(1))
     if (((PoissBox(gatedir,gatedir))/2.d0).le.Rmin_Gate) then
       @:ERROR_HANDLING(iErr, -5, 'Gate Distance too large')
     end if
  endif
  
  if (DoCilGate) then
    
     biasdir = abs(contdir(1))

     if (abs(bound(2)-bound(1)).le.(OxLength+dr_eps)) then
       @:ERROR_HANDLING(iErr, -6, 'Gate insulator is longer than Poisson box!')
     end if
     
     do i = 1,3
        if (i.eq.biasdir) then
          cycle
        end if
        if (((PoissBox(i,i))/2.d0).le.Rmin_Gate) then
          @:ERROR_HANDLING(iErr, -7, 'Gate transversal section is bigger than Poisson box!')
        end if
      end do
  end if

  !---------------------------------------
  ! Gate geometry mid point definition
  !---------------------------------------
  if (DoGate.or.DoCilGate) then
     do i = 1,3
        cntr_gate(i) = ( PoissBounds(i,2) + PoissBounds(i,1) )/2.d0
     end do
  end if

 end subroutine init_poissbox

!> This subroutine is a driver for the mudpack (c) solver (see mud3.f) *
subroutine mudpack_drv(env, SCC_in, V_L_atm, grad_V, iErr)

  !> Environment settings
  type(TEnvironment), intent(inout) :: env

  !> Control flag:
  integer, intent(in) :: SCC_in

  !> Outputs of subroutine "shift_Ham"
  real(kind=dp), intent(inout) :: V_L_atm(:,:)

  !> Output of subroutine "grad_V"
  real(kind=dp), intent(out) :: grad_V(:,:)

  !> Error return
  integer, intent(out), optional :: iErr

 integer, parameter :: GetPOT=0            !potential in SCC
 integer, parameter :: GetGRAD=1           !atomic shift component of gradient
 integer, parameter :: SavePOT=2           !SAVE CHR AND POTENTIAL 
 integer, parameter :: CLEAN=3             !DEALLOCATE PHI AND WORK

!Internal variables

 real(kind=dp) :: dlx,dly,dlz
 integer :: i,m,s
 integer :: worksize, ncycles, err
 integer :: iparm(23),mgopt(4)

 real(kind=dp) :: fparm(8)

 real(kind=dp), SAVE, ALLOCATABLE, DIMENSION (:,:,:) :: phi
 real(kind=dp), ALLOCATABLE, DIMENSION (:) :: work
 real(kind=dp), SAVE, ALLOCATABLE, DIMENSION (:,:,:) :: rhs 

 Type(super_array), ALLOCATABLE, DIMENSION(:) :: bulk

 integer :: isx, jsy, ksz
 integer, save :: niter = 1

 integer :: na,nb,nc, cont_mem
 character(50) :: BCinfo
 !e.g.: tmpdir (1) = 0 if there aren't two or more contacts in the same "x direction"

 if (present(iErr)) then
   iErr = 0
 end if

 iparm = 0

 if (SCC_in.ne.3) then

    do i = 1,3
       fparm(2*i) =   PoissBounds(i,2) 
       fparm(2*i-1) = PoissBounds(i,1)
    end do
            
    !**********************************************************************************
    ! 2. Setting number of equally spaced grid points in each interval defining the box 
    !**********************************************************************************
    ! number of point per axis: nx = ixp*(2**(iex-1))+1
    ! ixp+1 are the number of points in the coarsest grid 
    ! iparm(8) = ixp
    ! iparm(11)= iex
    ! iparm(14)= nx

    iparm(8) = 2   
    iparm(9) = 2
    iparm(10) = 2
        
    if (cluster.and.period) then
       if (period_dir(1)) iparm(8) = 3
       if (period_dir(2)) iparm(9) = 3
       if (period_dir(3)) iparm(10) = 3
    endif    

    do i = 1,50 
       iparm(11) = i
       iparm(14) = iparm(8)*(2**(iparm(11) - 1)) + 1 
       if (((fparm(2) - fparm(1))/(iparm(14) - 1)).le.dmin(1)) then
          dlx = (fparm(2) - fparm(1))/(iparm(14) - 1) 
          exit 
       end if
    end do
    do i = 1,50 
       iparm(12) = i
       iparm(15) = iparm(9)*(2**(iparm(12) - 1)) + 1 
       if (((fparm(4) - fparm(3))/(iparm(15) -1)).le.dmin(2)) then
          dly = (fparm(4) - fparm(3))/(iparm(15) -1) 
          exit 
       end if
    end do
    do i = 1,50 
       iparm(13) = i
       iparm(16) = iparm(10)*(2**(iparm(13) - 1)) + 1 
       if (((fparm(6) - fparm(5))/(iparm(16) - 1)).le.dmin(3)) then
          dlz = (fparm(6) - fparm(5))/(iparm(16) - 1)
          exit
       end if
    end do

    if (niter.eq.1.and.verbose.gt.30) then
       write(stdOut,'(73("-"))')
       write(stdOut,*) "Poisson Box internally adjusted:"
       write(stdOut,'(a,f12.5,f12.5,a11,l3)') ' x range=',PoissBounds(1,1)*Bohr__AA,&
            PoissBounds(1,2)*Bohr__AA,'; Periodic:',period_dir(1)
       write(stdOut,'(a,f12.5,f12.5,a11,l3)') ' y range=',PoissBounds(2,1)*Bohr__AA,&
            PoissBounds(2,2)*Bohr__AA,'; Periodic:',period_dir(2)
       write(stdOut,'(a,f12.5,f12.5,a11,l3)') ' z range=',PoissBounds(3,1)*Bohr__AA,&
            PoissBounds(3,2)*Bohr__AA,'; Periodic:',period_dir(3)

       write(stdOut,*) 'Mesh details:' 
       write(stdOut,'(a,f10.3,a,i4,a,f9.5)') ' Lx=',PoissBox(1,1)*Bohr__AA,&
            '  nx=',iparm(14),'   dlx=',dlx*Bohr__AA
       write(stdOut,'(a,f10.3,a,i4,a,f9.5)') ' Ly=',PoissBox(2,2)*Bohr__AA,&
            '  ny=',iparm(15),'   dly=',dly*Bohr__AA
       write(stdOut,'(a,f10.3,a,i4,a,f9.5)') ' Lz=',PoissBox(3,3)*Bohr__AA,&
            '  nz=',iparm(16),'   dlz=',dlz*Bohr__AA
       write(stdOut,'(73("-"))')
    endif

    !--------------------------
    ! Relaxation method used
    !--------------------------
    
    iparm(19) = 0           !0 => Gauss-Siedel method 
    isx = 0                 !        
    jsy = 0                 ! => relaxation method num. 0
    ksz = 0                 ! 
    iparm(20) = 0           !To set only for methods num. 8,9,10
        
    worksize = (iparm(14) + 2)*(iparm(15) + 2)*(iparm(16) +2)*(12+isx+jsy+ksz)
    iparm(21) = worksize 

    x0=(fparm(2)+fparm(1))/2.d0
    y0=(fparm(4)+fparm(3))/2.d0
    z0=(fparm(6)+fparm(5))/2.d0


    ! Special settings for mud3sp solver --------------------------------
    !if(cluster.and.period) then
    !   iparm(19) = 0 ! Gauss-Siedel 
    !   worksize = 7*(iparm(14) + 2)*(iparm(15) + 2)*(iparm(16) +2)/2
    !   iparm(20) = worksize 
    !endif
    ! -------------------------------------------------------------------
    
 end if !(SCC_in.ne.3)
 !**********************************************************************************
 ! 3. Allocate space for phi,rhs and work
 !**********************************************************************************
 if (id0 .and. .not.allocated(phi)) then
    call log_gallocate(phi,iparm(14),iparm(15),iparm(16))
    phi(:,:,:) = 0.d0
 end if
 if (id0 .and. .not.allocated(rhs)) then
    call log_gallocate(rhs,iparm(14),iparm(15),iparm(16))
    rhs(:,:,:) = 0.d0
 end if

 !**********************************************************************************
 ! 4. solve Poisson or compute atomic potentials
 !**********************************************************************************
 select case(SCC_in)

 !/////////////////////////////////////////////////////////////////
 case(GetPOT)     !Poisson called in order to calculate potential in SCC
 !/////////////////////////////////////////////////////////////////

   !**********************************************************************************
   ! 5.  Setting boundary conditions, iparm(2..7)
   !**********************************************************************************
   ! Default Setting for cluster or PBC in transverse directions
   ! flags BC on x=xa plane  !  0 = Periodic  
   ! flags BC on x=xb plane  !  1 = Dirichlet  
   ! flags BC on y=yc plane  !  2 = Mixed
   iparm(2) = 2 - 2*booltoint(period_dir(1))
   iparm(3) = 2 - 2*booltoint(period_dir(1))
   iparm(4) = 2 - 2*booltoint(period_dir(2))
   iparm(5) = 2 - 2*booltoint(period_dir(2))
   iparm(6) = 2 - 2*booltoint(period_dir(3))
   iparm(7) = 2 - 2*booltoint(period_dir(3))
    
   ! PoissBC = 0 for Dirichlet, 1 for Neumann   
   do m = 1, ncont
      select case(contdir(m))   
      case(-1)
         iparm(2) = 1+PoissBC(m)                         
      case(1)           
         iparm(3) = 1+PoissBC(m)                   
      case(-2)
         iparm(4) = 1+PoissBC(m)
      case(2)           
         iparm(5) = 1+PoissBC(m)
      case(-3)
         iparm(6) = 1+PoissBC(m)
      case(3)
         iparm(7) = 1+PoissBC(m)
      end select
    end do

    ! PLANAR GATE 
    if (DoGate) then
       select case(GateDir)
       case(1)     
         iparm(2) = 2                  
         iparm(3) = 2 
       case(2)      
         iparm(4) = 2                  
         iparm(5) = 2 
       case(3)
         iparm(6) = 2
         iparm(7) = 2
       end select         
    end if

    ! SETTING PBC or CLUSTER CALCULATIONS 
    ! 3D Periodic Boundary cond are handled fixing 2 faces with Ewalds sums
    if (cluster) then      
       iparm(2) = 1 - booltoint(period_dir(1))   ! if periodic => 0    
       iparm(3) = 1 - booltoint(period_dir(1))   ! else        => 1   
       iparm(4) = 1 - booltoint(period_dir(2))                 
       iparm(5) = 1 - booltoint(period_dir(2))                 
       iparm(6) = 1 - booltoint(period_dir(3))                 
       iparm(7) = 1 - booltoint(period_dir(3))
       
       if (period) then
          ! override Dirichlet on 2 opposite faces with smallest area
          na = iparm(14)
          nb = iparm(15)
          nc = iparm(16)
          overrideBC( 2*minloc( (/nb*nc,na*nc,na*nb/),1 )-1 ) = 1
          overrideBC( 2*minloc( (/nb*nc,na*nc,na*nb/),1 )   ) = 1
       end if

    end if        

    if(overrideBC(1).ne.0) iparm(2) = overrideBC(1)
    if(overrideBC(2).ne.0) iparm(3) = overrideBC(2)
    if(overrideBC(3).ne.0) iparm(4) = overrideBC(3)
    if(overrideBC(4).ne.0) iparm(5) = overrideBC(4)
    if(overrideBC(5).ne.0) iparm(6) = overrideBC(5)
    if(overrideBC(6).ne.0) iparm(7) = overrideBC(6)


    if(id0.and.niter.eq.1.and.verbose.gt.VBT) then
      write(stdOut,*) 'Boundary Conditions:' 
      BCinfo = 'x: '//boundary2string(iparm(2),mixed(1))
      if (overrideBC(1).ne.0) BCinfo = trim(BCinfo)//' (overridden)'
      BCinfo = trim(BCinfo)//'  '//boundary2string(iparm(3),mixed(2))
      if (overrideBC(2).ne.0) BCinfo = trim(BCinfo)//' (overridden)'
      write(stdOut,*) trim(BCinfo)

      BCinfo = 'y: '//boundary2string(iparm(4),mixed(3))
      if (overrideBC(3).ne.0) BCinfo = trim(BCinfo)//' (overridden)'
      BCinfo = trim(BCinfo)//'  '//boundary2string(iparm(5),mixed(4))
      if (overrideBC(4).ne.0) BCinfo = trim(BCinfo)//' (overridden)'
      write(stdOut,*) trim(BCinfo)

      BCinfo = 'z: '//boundary2string(iparm(6),mixed(5))
      if (overrideBC(5).ne.0) BCinfo = trim(BCinfo)//' (overridden)' 
      BCinfo = trim(BCinfo)//'  '//boundary2string(iparm(7),mixed(6))
      if (overrideBC(6).ne.0) BCinfo = trim(BCinfo)//' (overridden)'
      write(stdOut,*) trim(BCinfo)
    endif

    !------------------------------------------
    ! Setting initial guess for potential 
    !------------------------------------------
    
    if (niter.eq.1) then
       iparm(17) = 0           !no initial guess is provided
    else
       iparm(17) = 1           !previous phi provided as initial guess 
    end if
    
    !------------------------------------------------------------------------
    ! Max num of cycles executed between the finest and coarsest grid levels 
    !------------------------------------------------------------------------
    
    iparm(18) = MaxPoissIter
    
    !---------------------------------------------
    ! maximum relative error tolerance accepted 
    !---------------------------------------------
    
    fparm(7) = PoissAcc                         

   !-----------------------------------
   ! Selection for multigrid options  
   !-----------------------------------
   
    mgopt(1) = 0                            !Default multigrid with w-cycles               
                                    
    !--------------------------------------------------------------------------
    if (cluster.and.period .and. niter.eq.1) then
      call env%globalTimer%startTimer(globalTimers%poissonEwald)
      if (id0) then
        call set_phi_periodic(phi,iparm,fparm)
      end if
      call env%globalTimer%stopTimer(globalTimers%poissonEwald)
    end if
    !--------------------------------------------------------------------------
    if (ncont.gt.0) then

      allocate(bulk(ncont))
       call create_phi_bulk(bulk,iparm,dlx,dly,dlz,cont_mem)

       if(InitPot.and.id0.and.niter.eq.1.and.verbose.gt.VBT) then
         write(stdOut,*) 'Bulk Potential Info:'
         do m = 1,ncont 
            write(stdOut,*) 'contact',m
            if (bulk(m)%doEwald) then
               write(stdOut,*) 'BC = all periodic solved with Ewalds on two planes'  
            endif
            write(stdOut,*) 'x: '//boundary2string(bulk(m)%iparm(2))//'  '&
                            //boundary2string(bulk(m)%iparm(3))
            write(stdOut,*) 'y: '//boundary2string(bulk(m)%iparm(4))//'  '&
                            //boundary2string(bulk(m)%iparm(5))  
            write(stdOut,*) 'z: '//boundary2string(bulk(m)%iparm(6))//'  '&
                            //boundary2string(bulk(m)%iparm(7)) 
         enddo
       endif

       if(id0.and.niter.eq.1.and.verbose.gt.VBT) then
         write(stdOut,*) 'Memory required for Poisson:', & 
                          (2*size(phi)+iparm(21)+cont_mem)*8.d0/1d6,'Mb'
         
         write(stdOut,'(73("-"))')
       endif

       ! -----------------------------------------------------------------------
       ! Bulk boundary conditions setup 
       ! -----------------------------------------------------------------------
       if (InitPot) then
          
         if (ReadBulk) then
           !Read old bulk potential
           call env%globalTimer%startTimer(globalTimers%poissonBulkRead)
           call readbulk_pot(bulk)
           call env%globalTimer%stopTimer(globalTimers%poissonBulkRead)
         else
           !Compute bulk potential
           call env%globalTimer%startTimer(globalTimers%poissonBulkCalc)
           call compbulk_pot(bulk,iparm,fparm)
           ReadBulk=.true.
           call env%globalTimer%stopTimer(globalTimers%poissonBulkCalc)
         end if

       else  
         if(id0.and.verbose.gt.VBT) write(stdOut,*) 'No bulk potential'
       endif
       
    else
       ! allocate fake bulk to avoid problems    
       allocate(bulk(0))
    endif

    !write(stdOut,*) 'debug bulk potential'
    !do m=1,ncont
    !  call write_super_array(bulk(m)) 
    !  call save_bulkpot(bulk,m)
    !enddo
    
    ! -----------------------------------------------------------------------
    ! Simple boundary conditions setup
    ! -----------------------------------------------------------------------

    do m=1,ncont

      na = bulk(m)%iparm(14)
      nb = bulk(m)%iparm(15)
      nc = bulk(m)%iparm(16)
  
      bulk(m)%val(1:na,1:nb,1:nc) = bulk(m)%val(1:na,1:nb,1:nc) + mu(m)

      if (contdir(m).gt.0) then        
         s = iparm(abs(contdir(m))+13) 
      else
         s = 1
      end if
      
      ! Dirichlet BC
      if (id0) then
        select case(abs(contdir(m)))        
        case(1)
          do i = 1, nb 
            phi(s,1:na,i) = bulk(m)%val(1:na,i,1)
          enddo   
        case(2)     
          do i = 1, nb 
             phi(i,s,1:na) = bulk(m)%val(1:na,i,1)
          enddo   
        case(3)
          do i = 1, nb 
             phi(1:na,i,s) = bulk(m)%val(1:na,i,1)
          enddo   
        end select
      endif

    enddo
    
    !*********************************************************************************
    ! Charge density evaluation on the grid points 
    !*********************************************************************************

    call set_rhs(env, iparm,fparm,dlx,dly,dlz,rhs,bulk)


    !*********************************************************************************
    !---------------------------------------------------------------------------------- 

    !**********************************************************************************
    ! 6.  Solving Poisson equation                 
    !**********************************************************************************
    if (id0) then

      call log_gallocate(work,worksize)

      call env%globalTimer%startTimer(globalTimers%poissonSoln)
      do i = 0,1
        iparm(1) = i

        if (DoGate) then
           call mud3(iparm,fparm,work,coef_gate,bndyc,rhs,phi,mgopt,err)
        elseif (DoCilGate) then
           call mud3(iparm,fparm,work,coef_cilgate,bndyc,rhs,phi,mgopt,err)
        elseif (DoTip) then
           call mud3(iparm,fparm,work,coef_tip,bndyc,rhs,phi,mgopt,err)
        elseif (cluster.and.period) then
           call mud3(iparm,fparm,work,coef,bndyc,rhs,phi,mgopt,err)
        else
           call mud3(iparm,fparm,work,coef,bndyc,rhs,phi,mgopt,err)
        end if
        
        worksize = iparm(22)

        if (err.ne.0.and.err.ne.9) then
          if(err.gt.0) then
            @:FORMATTED_ERROR_HANDLING(iErr, err, "(A,I0)", 'Fatal Error in poisson solver:', err)
          end if
        end if
        if (err.eq.9) then
          call log_gdeallocate(work)
          call log_gallocate(work,worksize)
        end if
      end do

      call env%globalTimer%stopTimer(globalTimers%poissonSoln)

      if (err.lt.-1) then
        write(stdOut,*) 'Non-fatal Error in poisson solver:',err
      endif

      ncycles = iparm(23)
      call log_gdeallocate(work)
    end if

    if (id0) then
      if (verbose.gt.30) then 
        write(stdOut,'(1x,73("-"))') 
        write(stdOut,*) 'Relative Poisson Error =',fparm(8)
        write(stdOut,'(a,i3,a,i3)') ' Number of cycles executed =',ncycles,'/',iparm(18)
        write(stdOut,'(1x,73("-"))') 
        flush(6)
      end if   

      if (err.eq.-1 .or. ncycles.eq.iparm(18)) then
        @:ERROR_HANDLING(iErr, -1, 'Convergence in Poisson solver not obtained')
      end if

    end if

    !--------------------------------------------
    ! Shift of the Hamiltonian matrix elements 
    !--------------------------------------------
    if (id0) then
      call env%globalTimer%startTimer(globalTimers%poissonShifts)
      call shift_Ham(iparm,fparm,dlx,dly,dlz,phi,bulk,V_L_atm)
      call env%globalTimer%stopTimer(globalTimers%poissonShifts)
    end if
    
    call destroy_phi_bulk(bulk)

    deallocate(bulk,stat=err)

 !//////////////////////////////////////////////////////////////////////
 case(GetGRAD)    ! Poisson called in order to calculate atomic shift gradient 
 !//////////////////////////////////////////////////////////////////////
   
   if (id0) call gradient_V(phi,iparm,fparm,dlx,dly,dlz,grad_V)
   
   
 !////////////////////////////////////////////////////////////////////////
 case(SavePOT)    ! Poisson called in order to save potential and charge density
 !////////////////////////////////////////////////////////////////////////

    if (id0) call save_pot(iparm,fparm,dlx,dly,dlz,phi,rhs)
   
 !///////////////////////////////////////////
 case(CLEAN)       ! Deallocate Poisson variables
 !///////////////////////////////////////////      
       
   if (allocated(phi)) call log_gdeallocate(phi)
   if (allocated(rhs)) call log_gdeallocate(rhs)
   niter = 0

 end select

  niter=niter+1

end subroutine Mudpack_drv

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine set_rhs(env, iparm, fparm, dlx, dly, dlz, rhs, bulk)

  !> Environment settings
  type(TEnvironment), intent(inout) :: env

  integer :: iparm(23)
  real(kind=dp) :: fparm(8),dlx,dly,dlz
  real(kind=dp), dimension(:,:,:) :: rhs
  type(super_array) :: bulk(:)

 #:if WITH_MPI
  !---------------------------------------------------------------------
  ! MPI parallelization of the r.h.s. assembly (charge density)
  ! This is done slicing the grid along the z-direction, iparm(16)
  ! Parallelization is done along z since this corresponds to the slowest
  ! index in rhs(:,:,:), therefore gather is done on a contiguus vector
  !
  !---------------------------------------------------------------------
  integer :: i, npid, iparm_tmp(23)
  real(kind=dp) :: fparm_tmp(8)
  real(kind=dp), ALLOCATABLE, DIMENSION (:,:,:) :: rhs_par
  integer, ALLOCATABLE, DIMENSION (:) :: istart,iend,dim_rhs
 
  if (numprocs > 1) then
 
    call log_gallocate(dim_rhs,numprocs)
    call log_gallocate(istart,numprocs)
    call log_gallocate(iend,numprocs)
 
    ! z points per processor
    npid = int(iparm(16)/numprocs)
 
    ! set start/end and size handled by each processor
    do i = 1,numprocs
      istart(i) = (i-1)*npid+1
      if (i .ne. numprocs) then
         iend(i) = i*npid
      else
         iend(i) = iparm(16)
      endif
      dim_rhs(i) = iparm(14)*iparm(15)*(iend(i)-istart(i)+1)
    end do
    ! Define a subproblem with appropriate iparm_tmp
    iparm_tmp = iparm
    fparm_tmp = fparm
 
    iparm_tmp(16) = iend(id+1) - istart(id+1) + 1
    fparm_tmp(5) = fparm(5) + (istart(id+1)-1)*dlz
    fparm_tmp(6) = fparm(5) + (iend(id+1)-1)*dlz
 
    ! Allocate space for local rhs.
    call log_gallocate(rhs_par,iparm_tmp(14),iparm_tmp(15),iparm_tmp(16))
    !---------------------------------------------------------------------

    call env%globalTimer%startTimer(globalTimers%poissonDensity)

    !! Set a renormalization volume for grid projection
 
    if (do_renorm) then
      call renormalization_volume(iparm_tmp,fparm_tmp,dlx,dly,dlz,fixed_renorm)
      do_renorm = .false.
    endif
 
    call charge_density(iparm_tmp,fparm_tmp,dlx,dly,dlz,rhs_par)
 
    if (DoGate) then
       call gate_bound(iparm_tmp,fparm_tmp,dlx,dly,dlz,rhs_par)
    endif
 
    if (DoCilGate) then
       call cilgate_bound(iparm_tmp,fparm_tmp,dlx,dly,dlz,rhs_par)
    endif
 
    if (DoTip) then
       call tip_bound(iparm_tmp,fparm_tmp,dlx,dly,dlz,rhs_par)
    endif

    call env%globalTimer%stopTimer(globalTimers%poissonDensity)

    ! gather all partial results on lead node 0
    call mpifx_gatherv(poiss_comm, rhs_par, rhs, dim_rhs)
 
    call log_gdeallocate(rhs_par)
    call log_gdeallocate(dim_rhs)
    call log_gdeallocate(istart)
    call log_gdeallocate(iend)
 
  else
 
 #:endif

    if (do_renorm) then
      call renormalization_volume(iparm,fparm,dlx,dly,dlz,fixed_renorm)
      do_renorm = .false.
    endif
   
    call charge_density(iparm,fparm,dlx,dly,dlz,rhs)

    if (DoGate) then
      call gate_bound(iparm,fparm,dlx,dly,dlz,rhs)
    endif
 
    if (DoCilGate) then
      call cilgate_bound(iparm,fparm,dlx,dly,dlz,rhs)
    endif
 
    if (DoTip) then
      call tip_bound(iparm,fparm,dlx,dly,dlz,rhs)
    endif

 #:if WITH_MPI
  endif
 #:endif

  if (any(localBC.gt.0)) then
    call local_bound(iparm,fparm,x,rhs,bulk)
  endif 

end subroutine set_rhs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine renormalization_volume(iparm, fparm, dlx, dly, dlz, fixed)

  implicit none 

  integer :: iparm(23) 
  real(kind=dp) :: fparm(8)
  real(kind=dp) :: dlx, dly, dlz
  logical, intent(in) :: fixed

  !Internal variables

  real(kind=dp) :: xi(3), deltaR
  integer :: i, j, k, atom, rag(3)
  integer :: ragx, ragy, ragz, npx, npy, npz

  real(kind=dp) :: dl(3), xmin(3), xmax(3)
  integer :: imin(3),imax(3), ii, jj, kk, l, nsh
  real(kind=dp), allocatable, dimension(:,:) :: tau
  ! Perform almost the same operations of charge_density, estimates
  ! for each atom the renormalization factors (inverse of exponential
  ! spherical shell volume), to cancel out 
  ! discretization and truncation errors (this runs on contact atoms 
  ! as well)

  dl(1)=dlx; dl(2)=dly; dl(3)=dlz;

  do i = 1,3
    if (period_dir(i)) then
      rag(i) = 1
    else
      rag(i) = 0    
    end if
  end do

  ! define aliases
  ragx = rag(1); ragy = rag(2); ragz = rag(3)
  npx = iparm(14)-ragx; npy = iparm(15)-ragy; npz = iparm(16)-ragz

  ! Strategy: along periodic direction we need to avoid the double counting
  ! introduced by the last grid-point that should fold on the first.
  ! for this reason on periodic dir we compute the point mod(np-1).

  tau = 3.2d0 * uhubb

  if (.not.(allocated(renorm))) call log_gallocate(renorm,maxval(lmax)+1,natoms)
  renorm = 0.0

  if (fixed) then
    ! This implement the old formula 8pi/tau**3
    do atom = 1, natoms
       nsh = lmax(izp(atom))+1
      do l = 1, nsh
        renorm(l,atom) = (8.d0*Pi)/(tau(l,izp(atom))**3)
      enddo
    enddo
  else

    do atom = 1, natoms !istart(id+1), iend(id+1)

      nsh = lmax(izp(atom))+1

      ! Set boundaries of a box around the atom 
      xmin(:) = x(:,atom) - deltaR_max
      xmax(:) = x(:,atom) + deltaR_max

      ! Cut atomic box out of PoissonBox
      do i=1,3
        imin(i) = nint( (xmin(i) - fparm(2*i-1))/dl(i) ) !+ 1 
        imax(i) = nint( (xmax(i) - fparm(2*i-1))/dl(i) ) !+ 1
        ! GP In non periodic directions DO NOT cut at boundaries. In this way we
        ! ensure the right renormalization volume also for atoms lying close
        ! to boundaries. The lines commented below are the corresponding ones 
        ! from the projection subroutine
        ! ----------------------------------------------------------
        ! in non periodic directions cut at PoissonBox boundaries 
        !if (.not.period_dir(i)) then
        !  imin(i) = max( 0, imin(i) )
        !  imax(i) = min( iparm(13+i)-1, imax(i) )
        !endif
        ! -----------------------------------------------------------
      enddo

      ! compute the charge density on the grid points
      do i = imin(1),imax(1)      
        xi(1) = fparm(1) + i*dlx
        ii=modulo(i, npx) + 1        
        do j = imin(2),imax(2)         
          xi(2) = fparm(3) + j*dly 
          jj=modulo(j, npy) + 1            
          do k = imin(3),imax(3)                       
            xi(3) = fparm(5) + k*dlz
            kk=modulo(k, npz) + 1  

            !Compute distance from atom
            deltaR = sqrt(dot_product(xi-x(:,atom),xi-x(:,atom)))
            ! add charge density contrib
            do l = 1, nsh
              renorm(l,atom) = renorm(l,atom) + exp(-tau(l,izp(atom))*deltaR) *dlx*dly*dlz
            enddo

          end do
        end do
      end do  
    end do
  endif

  !Inverting total volume to calculate the renormalization factors
  renorm = 1.d0/renorm

end subroutine renormalization_volume
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine charge_density(iparm,fparm,dlx,dly,dlz,rhs)

 implicit none 

 integer :: iparm(23)
 real(kind=dp) :: fparm(8)
 real(kind=dp) :: dlx,dly,dlz
 real(kind=dp) :: rhs(:,:,:)

 !Internal variables

 real(kind=dp) :: xi(3),deltaR
 integer :: i,j,k,atom, rag(3)
 integer :: ragx, ragy, ragz, npx, npy, npz

 real(kind=dp) :: tmp,dl(3),xmin(3),xmax(3)
 integer :: imin(3),imax(3), ii, jj, kk, l, nsh
 real(kind=dp), allocatable, dimension(:,:) :: tau

 rhs(:,:,:)=0.d0

 dl(1)=dlx; dl(2)=dly; dl(3)=dlz;
 
 do i = 1,3
    if (period_dir(i)) then
       rag(i) = 1
    else
       rag(i) = 0    
    end if
 end do

 ! define aliases
 ragx = rag(1); ragy = rag(2); ragz = rag(3)
 npx = iparm(14)-ragx; npy = iparm(15)-ragy; npz = iparm(16)-ragz

 
 !call distribute_atoms(1, natoms, 1, istart, iend, dims, displ)

 ! Strategy: along periodic direction we need to avoid the double counting
 ! introduced by the last grid-point that should fold on the first.
 ! for this reason on periodic dir we compute the point mod(np-1).
 
 ! Set for each shell a renormalization term which correctly 
 ! takes into account cutoff of exponential as well
 tau = 3.2d0 * uhubb

 do atom = 1, natoms !istart(id+1), iend(id+1)
    tmp=0.d0
    ! project the atom on the primitive cell
    !patom=mod(atom-1,natoms)+1    
    nsh = lmax(izp(atom))+1
    
    ! Set boundaries of a box around the atom 
    xmin(:) = x(:,atom) - deltaR_max
    xmax(:) = x(:,atom) + deltaR_max
    
    ! Cut atomic box out of PoissonBox
    do i=1,3
       imin(i) = nint( (xmin(i) - fparm(2*i-1))/dl(i) ) !+ 1 
       imax(i) = nint( (xmax(i) - fparm(2*i-1))/dl(i) ) !+ 1
       ! in non periodic directions cut at PoissonBox boundaries 
       if (.not.period_dir(i)) then
          imin(i) = max( 0, imin(i) )
          imax(i) = min( iparm(13+i)-1, imax(i) )
       endif
    enddo
    
    ! compute the charge density on the grid points
    do i = imin(1),imax(1)
       
       xi(1) = fparm(1) + i*dlx
       ii=modulo(i, npx) + 1  
       
       do j = imin(2),imax(2)
          
          xi(2) = fparm(3) + j*dly 
          jj=modulo(j, npy) + 1  
          
          do k = imin(3),imax(3)
                       
             xi(3) = fparm(5) + k*dlz
             kk=modulo(k, npz) + 1  

             !Compute distance from atom
             deltaR = sqrt(dot_product(xi-x(:,atom),xi-x(:,atom)))
             ! add charge density contrib
             do l = 1, nsh
               !tmp=3.2d0*uhubb(l,izp(atom))
               rhs(ii,jj,kk) = rhs(ii,jj,kk) - 4.d0 * Pi * renorm(l,atom) *&
                   & dQmat(l,atom) * exp(-tau(l,izp(atom))*deltaR)
               !note: rhs contains -4*pi in normalization term: 
               !Volume renormalization is in renorm and includes discretization and
               !truncation error, or it is the usual 8pi/tau**3 depending on options
             enddo

          end do
       end do
    end do
 end do


 ! Copy the ragged periodic copy
 if (period_dir(1)) rhs(npx+ragx,:,:) = rhs(1,:,:)
 if (period_dir(2)) rhs(:,npy+ragy,:) = rhs(:,1,:)
 if (period_dir(3)) rhs(:,:,npz+ragz) = rhs(:,:,1)
      
end subroutine charge_density


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function boundary2string(typ,local) result(str)
   integer :: typ
   character(10) :: str
   logical, optional :: local
              
   select case(typ)
   case(0)
       str=' Periodic '
   case(1)
       str='Dirichlet'
   case(2)
      str='  Neumann '
      if (present(local)) then
        if (local) then
           str='    Mixed '
        endif
      endif
   end select

end function boundary2string


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine shift_Ham(iparm,fparm,dlx,dly,dlz,phi,phi_bulk,V_atm)
 
  integer, intent(in) :: iparm(23)
  real(kind=dp), intent(in) :: fparm(8)
  real(kind=dp), intent(in) :: dlx,dly,dlz
  real(kind=dp), intent(in) :: phi(:,:,:)
  type(super_array), intent(in) :: phi_bulk(:)
  real(kind=dp), dimension(:,:), intent(inout) :: V_atm
 
  !Internal variables
 
  integer :: i,j,k,atm,f,s,pl
  real(kind=dp) :: g,Norm,xi(3),deltaR, vol, V_tmp, expgr
  integer :: m,a,b,c
 
  real(dp) :: dl(3), xmin(3), xmax(3), xhlp(3), dla, dlb, dlc
  integer :: imin(3), imax(3), n_cell(3), ii, jj, kk, rag(3)
  integer :: ncx,ncy,ncz, npx, npy, npz, nsh,l
 
  dl(1)=dlx; dl(2)=dly; dl(3)=dlz;
 
  do i = 1,3
     if (period_dir(i)) then
        n_cell(i) = nint(PoissBox(i,i)/boxsiz(i,i)) + 2 
        rag(i) = 1
     else
        n_cell(i) = 0
        rag(i) = 0    
     end if
  end do
 
  ! NOTE: DO NOT INITIALIZE V_atm = 0 here
  
  ! define aliases
  ncx = n_cell(1); ncy = n_cell(2); ncz = n_cell(3)
  npx = iparm(14)-rag(1); npy = iparm(15)-rag(2); npz = iparm(16)-rag(3)
 
  atoms: do atm = 1, iatm(2)  
 
     xhlp(:)=x(:,atm)
     do while (xhlp(1).lt.fparm(1))
        xhlp(1)=xhlp(1)+PoissBox(1,1)
     enddo
     do while (xhlp(1).gt.fparm(2))
        xhlp(1)=xhlp(1)-PoissBox(1,1)
     enddo
     do while (xhlp(2).lt.fparm(3))
        xhlp(2)=xhlp(2)+PoissBox(2,2)
     enddo
     do while (xhlp(2).gt.fparm(4))
        xhlp(2)=xhlp(2)-PoissBox(2,2)
     enddo
     do while (xhlp(3).lt.fparm(5))
        xhlp(3)=xhlp(3)+PoissBox(3,3)
     enddo
     do while (xhlp(3).gt.fparm(6))
        xhlp(3)=xhlp(3)-PoissBox(3,3)
     enddo
 
     ! Set boundaries of a box around the atom 
     xmin(:) = xhlp(:) - deltaR_max
     xmax(:) = xhlp(:) + deltaR_max
     
     ! Cut out box out of PoissonBox imin imax start from 0
     do i=1,3
        imin(i) = nint( (xmin(i) - fparm(2*i-1))/dl(i) ) !+ 1 
        imax(i) = nint( (xmax(i) - fparm(2*i-1))/dl(i) ) !+ 1
        if (.not.period_dir(i)) then
           imin(i) = max( 0, imin(i) )
           imax(i) = min( iparm(13+i)-1, imax(i) )
        endif
     enddo
        
     nsh = lmax(izp(atm))+1     
     shells: do l = 1, nsh
        g = 3.2d0*uhubb(l,izp(atm))
        vol = dlx*dly*dlz 
        Norm=0.d0
        V_tmp = 0.d0
                 
        do i = imin(1),imax(1)
           
           xi(1) = fparm(1) + i*dlx
           ii=mod(i + ncx*npx, npx) + 1  
           
           do j = imin(2),imax(2)
              
              xi(2) = fparm(3) + j*dly 
              jj=mod(j + ncy*npy, npy) + 1       
              
              do k = imin(3),imax(3)
                 !Compute point coordinates
              
                 xi(3) = fparm(5) + k*dlz
                 kk=mod(k + ncz*npz, npz) + 1
 
                 ! Compute distance from atom 
                 deltaR = sqrt(dot_product(xi-xhlp, xi-xhlp))
                 if (deltaR.gt.deltaR_max) then
                    cycle
                 else
                   expgr = exp(-g*deltaR)
                   V_tmp = V_tmp - phi(ii,jj,kk)*expgr*vol
                   !(Normalization g**3*dlx*dly*dlz/(8.0*pi) included in Norm)
                   Norm = Norm + expgr*vol
                 end if
              end do
           end do
        end do        
        
        !-------------------------------------------------------------------------------
        ! Hamiltonian shifts correction with contacts bulk-potential or bias 
        !-------------------------------------------------------------------------------
        do m = 1, ncont
 
           s=sign(1,contdir(m))
           f=(s-1)/2
 
           a=phi_bulk(m)%a
           b=phi_bulk(m)%b
           c=phi_bulk(m)%c
 
           dla=phi_bulk(m)%dla
           dlb=phi_bulk(m)%dlb
           dlc=phi_bulk(m)%dlc
 
           vol = dla*dlb*dlc 
 
           if(abs( xhlp(c)-fparm( 2*c + f ) ).lt.deltaR_max) then
           !We add the shift due to 2nd PL, too. This is consistent
           !with the calculation of chrge density which uses both
           !contact PLs 
           do pl = 1,2
              do i = 1,phi_bulk(m)%iparm(14)
                 do j = 1,phi_bulk(m)%iparm(15)
                    do k = 2,phi_bulk(m)%iparm(16)
                       
                       xi(a) = fparm( 2*a-1 ) +   (i - 1)*dla 
                       xi(b) = fparm( 2*b-1 ) +   (j - 1)*dlb    
                       !This line is modified to include 2nd PL 
                       xi(c) = fparm( 2*c+f ) + s*(k - 1)*dlc +s*(pl-1)*phi_bulk(m)%L_PL
                       
 
                       ! This is relying on the fact that the grid of the contact bulk 
                       ! potential is ordered from the device outwards. 
 
                       deltaR = sqrt(dot_product(xi-xhlp, xi-xhlp))
                       if (deltaR.gt.deltaR_max) then
                          cycle
                       else    
                          expgr = exp(-g*deltaR)
                          V_tmp = V_tmp - phi_bulk(m)%val(i,j,k)*expgr*vol
                          Norm = Norm + expgr*vol
                       end if
                    end do
                 end do
              end do
           end do
 
           end if
 
        end do
        !------------------------------------------------------------------------------ 
        V_atm(l,atm)= V_tmp/Norm
 
     end do shells
  end do atoms
 
  ! biasshift here 
  do m = 1, ncont
     do atm=iatc(3,m),iatc(2,m)
        do l = 1, lmax(izp(atm))+1 
           V_atm(l,atm) = V_atm(l,atm) - mu(m) 
        end do
     end do
  end do
 
end subroutine shift_Ham

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine gradient_V(phi,iparm,fparm,dlx,dly,dlz,grad_V)

  real(kind=dp), intent(in) :: phi(:,:,:)
  integer, intent(in) :: iparm(23) 
  real(kind=dp), intent(in) :: fparm(8)
  real(kind=dp), intent(in) :: dlx,dly,dlz
  real(kind=dp), intent(inout), dimension(:,:) :: grad_V
  
 
  integer :: i,j,k,atm, l, nsh
  real(kind=dp) :: xi(3),g,deltaR,tmp_Gr(3)
 
  real(dp) :: dl(3), xmin(3), xmax(3), tmp, Norm
  integer :: imin(3), imax(3), n_cell(3), ii, jj, kk
  integer :: ncx,ncy,ncz, ragx, ragy, ragz, npx, npy, npz, rag(3)
 
  grad_V = 0.0_dp
 
  dl(1)=dlx; dl(2)=dly; dl(3)=dlz;
 
  do i = 1,3
     if (period_dir(i)) then
        n_cell(i) = nint(PoissBox(i,i)/boxsiz(i,i)) + 2
        rag(i) = 1
     else
        n_cell(i) = 0
        rag(i) = 0
     end if
  end do
 
  ! define aliases
  ragx = rag(1); ragy = rag(2); ragz = rag(3)
  ncx = n_cell(1); ncy = n_cell(2); ncz = n_cell(3)
  npx = iparm(14)-ragx; npy = iparm(15)-ragy; npz = iparm(16)-ragz
  
  do atm = iatm(1),iatm(2)
     
     ! Set boundaries of a box around the atom 
     xmin(:) = x(:,atm) - deltaR_max
     xmax(:) = x(:,atm) + deltaR_max
     
     ! Cut out box out of PoissonBox
     do i=1,3
        imin(i) = int( (xmin(i) - fparm(2*i-1))/dl(i) ) !+ 1 
        imax(i) = int( (xmax(i) - fparm(2*i-1))/dl(i) ) !+ 1
        if (.not.period_dir(i)) then
           imin(i) = max( 0, imin(i) )
           imax(i) = min( iparm(13+i)-1, imax(i) )
        endif
     enddo
                
     nsh = lmax(izp(atm))+1
     do l = 1, nsh                      
       g = 3.2d0*uhubb(l,izp(atm))
       tmp_Gr = 0.d0
       Norm = 0.d0
 
       do i = imin(1),imax(1)
          
          xi(1) = fparm(1) + i*dlx
          ii=mod(i + ncx*npx, npx) + 1  
          
          do j = imin(2),imax(2)
             
             xi(2) = fparm(3) + j*dly 
             jj=mod(j + ncy*npy, npy) + 1       
             
             do k = imin(3),imax(3)
                !Compute point coordinates
                
                xi(3) = fparm(5) + k*dlz
                kk=mod(k + ncz*npz, npz) + 1
                
                ! Compute distance from atom 
                
                deltaR = sqrt(dot_product(xi(:)-x(:,atm), xi(:)-x(:,atm)))
                if (deltaR.gt.deltaR_max .or. deltaR.eq.0.d0) then
                   cycle
                else
 
                   tmp= exp(-g*deltaR)/deltaR
                   tmp_Gr = tmp_Gr + tmp*(xi(:)-x(:,atm))*phi(ii,jj,kk)
                   
                   Norm = Norm + exp(-g*deltaR) 
                   !(Normalization g**3*dlx*dly*dlz/(8.0*pi) included in Norm)
                end if
             end do
          end do
       end do
 
       grad_V(:,atm)= grad_V(:,atm)+ tmp_gr(:)*g*dQmat(l,atm)/Norm
                  
     end do 
  end do

end subroutine gradient_V


subroutine save_pot(iparm,fparm,dlx,dly,dlz,phi,rhs)

  integer, intent(in) :: iparm(23) 
  real(kind=dp), intent(in) :: fparm(8)
  real(kind=dp), intent(in) :: dlx,dly,dlz
  real(kind=dp), intent(in) :: phi(:,:,:)
  real(kind=dp), intent(in) :: rhs(:,:,:)

  integer :: i,j,k,nx_fix,ny_fix,nz_fix,FixDir, fp
  real(kind=dp) :: xi,yj,zk
  real(kind=dp) :: z_min_gate, z_max_gate, z_min_ox, z_max_ox 

  FixDir=int(PoissPlane(1)) 

  if (verbose.gt.70) then 
    if(id0) write(stdOut,'(1x,a)') 'Saving charge density and potential ...'
  endif

  ! Saving 3D potential and density
  !--------------------------------------------  
  if (id0.and.(FixDir.eq.0)) then
     open(newunit=fp,file='box3d.dat')
     write(fp,*) iparm(14),iparm(15),iparm(16)
     close(fp)

     open(newunit=fp,file='Xvector.dat')
     do i = 1,iparm(14)  
        xi = fparm(1) + (i - 1)*dlx
        xi = xi*Bohr__AA
        write(fp,'(E17.8)',ADVANCE='NO') xi 
     end do
     close(fp)

     open(newunit=fp,file='Yvector.dat')
     do j = 1,iparm(15)
        yj = fparm(3) + (j - 1)*dly  
        yj = yj*Bohr__AA
        write(fp,'(E17.8)',ADVANCE='NO') yj
     end do
     close(fp)

     open(newunit=fp,file='Zvector.dat')
     do k = 1,iparm(16)  
        zk = fparm(5) + (k - 1)*dlz
        zk = zk*Bohr__AA
        write(fp,'(E17.8)',ADVANCE='NO') zk 
     end do
     close(fp)
     
     open(newunit=fp,file='potential.dat') 
     do i = 1,iparm(14)  
        do j = 1,iparm(15)
           do k = 1,iparm(16)
              write(fp,'(E17.8)') phi(i,j,k)*hartree__eV
           end do
        end do
     end do
     close(fp)
     
     open(newunit=fp,file='charge_density.dat') 
     do i = 1,iparm(14)  
        do j = 1,iparm(15)
           do k = 1,iparm(16)
              write(fp,'(E17.8)') rhs(i,j,k)/(-4.d0*Pi)
           end do
        end do
     end do
     close(fp)

  endif
  
  !--------------------------------------------
  ! Saving just on a plane of the Poisson box
  !--------------------------------------------
  
  if (id0.and.(FixDir.ne.0)) then   

    select case(FixDir)
    case(1)
      nx_fix = nint(((fparm(2)-fparm(1))*PoissPlane(2))/dlx) + 1            
      open(newunit=fp,file='pot2D.dat') 
      do j = 1,iparm(15)
        yj = fparm(3) + (j - 1)*dly 
        do k = 1,iparm(16)
           zk = fparm(5) + (k - 1)*dlz
           write(fp,'(E17.8,E17.8,E17.8)') yj*Bohr__AA, zk*Bohr__AA, phi(nx_fix&
               &,j,k)*hartree__eV
           
         end do
       end do
       close(fp)
       
       open(newunit=fp,file='chr2D.dat') 
       do j = 1,iparm(15)
         yj = fparm(3) + (j - 1)*dly 
         do k = 1,iparm(16)
           zk = fparm(5) + (k - 1)*dlz
           write(fp,'(E17.8,E17.8,E17.8)') yj*Bohr__AA, zk*Bohr__AA, rhs(nx_fix&
               &,j,k)/(-4.0*4.0*atan(1.d0))
           
         end do
       end do
       close(fp)
       
       open(newunit=fp,file='box2d.dat')
       write(fp,'(I6,I6)') iparm(15),iparm(16)
       close(fp)
       
     case(2)
       ny_fix = nint(((fparm(4)-fparm(3))*PoissPlane(2))/dly) + 1            
       open(newunit=fp,file='pot2D.dat') 
       do i = 1,iparm(14)
         xi = fparm(1) + (i - 1)*dlx  
         do k = 1,iparm(16)
           zk = fparm(5) + (k - 1)*dlz
           write(fp,'(E17.8,E17.8,E17.8)') xi*Bohr__AA, zk*Bohr__AA, phi(i,ny_fix,k)*hartree__eV
         end do
       end do
       close(fp)
       
       open(newunit=fp,file='chr2D.dat') 
       do i = 1,iparm(14)
         xi = fparm(1) + (i - 1)*dlx  
         do k = 1,iparm(16)
           zk = fparm(5) + (k - 1)*dlz
           write(fp,'(E17.8,E17.8,E17.8)')  xi*Bohr__AA, zk*Bohr__AA, rhs(i&
               &,ny_fix,k)/(-4.0*4.0*atan(1.d0))
           
         end do
       end do
       close(fp)
       
       open(fp,file='box2d.dat')
       write(fp,'(I6,I6)') iparm(14),iparm(16)
       close(fp)
                
     case(3)
       nz_fix = nint(((fparm(6)-fparm(5))*PoissPlane(2))/dlz) + 1            
       open(fp,file='pot2D.dat') 
       do i = 1,iparm(14)
         xi = fparm(1) + (i - 1)*dlx  
         do j = 1,iparm(15)
           yj = fparm(3) + (j - 1)*dly
           write(fp,'(E17.8,E17.8,E17.8)') xi*Bohr__AA, yj*Bohr__AA, phi(i,j&
               &,nz_fix)*hartree__eV
           
         end do
       end do
       close(fp)
       
       open(newunit=fp,file='chr2D.dat') 
       do i = 1,iparm(14)
         xi = fparm(1) + (i - 1)*dlx  
         do j = 1,iparm(15)
           yj = fparm(3) + (j - 1)*dly
           write(fp,'(E17.8,E17.8,E17.8)') xi*Bohr__AA, yj*Bohr__AA, rhs(i,j&
               &,nz_fix)/(-4.0*4.0*atan(1.d0))
           
         end do
       end do
       close(fp)        
       
       open(newunit=fp,file='box2d.dat')
       write(fp,'(I6,I6)') iparm(14),iparm(15)
       close(fp)
       
     end select
   end if
   
   !-----------------------------
   ! Saving gate geometry data
   !-----------------------------
   
   if (id0.and.DoCilGate) then             
     z_min_gate = cntr_gate(biasdir) - GateLength_l/2.d0  
     z_max_gate = cntr_gate(biasdir) + GateLength_l/2.d0
     z_min_ox = cntr_gate(biasdir) - OxLength/2.d0  
     z_max_ox = cntr_gate(biasdir) + OxLength/2.d0  
     open(newunit=fp,file='gate.dat')
     write(fp,'(i2)') biasdir
     write(fp,'(E17.8,E17.8)') z_min_gate*Bohr__AA,z_max_gate*Bohr__AA
     write(fp,'(E17.8,E17.8)') z_min_ox*Bohr__AA,z_max_ox*Bohr__AA
     write(fp,'(E17.8,E17.8)') Rmin_Gate*Bohr__AA,Rmin_Ins*Bohr__AA
     write(fp,'(E17.8,E17.8)') cntr_gate(1)*Bohr__AA,cntr_gate(2)*Bohr__AA,cntr_gate(3)*Bohr__AA
     close(fp)
   end if

   if (id0.and.DoGate) then             
     open(newunit=fp,file='gate.dat')
     write(fp,'(i2)') gatedir, biasdir

     z_min_gate = cntr_gate(biasdir) - GateLength_l/2.d0  
     z_max_gate = cntr_gate(biasdir) + GateLength_l/2.d0
     write(fp,'(E17.8,E17.8)') z_min_gate*Bohr__AA,z_max_gate*Bohr__AA
     
     do i=1,3
       if (i.ne.gatedir .and. i.ne.biasdir) exit
     enddo
     z_min_gate = cntr_gate(i) - GateLength_t/2.d0  
     z_max_gate = cntr_gate(i) + GateLength_t/2.d0
     write(fp,'(E17.8,E17.8)') z_min_gate*Bohr__AA,z_max_gate*Bohr__AA
     
     z_min_ox = cntr_gate(gatedir) - OxLength/2.d0  
     z_max_ox = cntr_gate(gatedir) + OxLength/2.d0  
     write(fp,'(E17.8,E17.8)') z_min_ox*Bohr__AA,z_max_ox*Bohr__AA
     write(fp,'(E17.8,E17.8)') Rmin_Gate*Bohr__AA,Rmin_Ins*Bohr__AA
     write(fp,'(E17.8,E17.8)') cntr_gate(1)*Bohr__AA,cntr_gate(2)*Bohr__AA,cntr_gate(3)*Bohr__AA
     close(fp)
   end if
     
 end subroutine save_pot
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ! Perform a standard gamma-functional calculation for periodic systems 
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine set_phi_periodic(phi, iparm, fparm)
   real(dp), dimension(:,:,:) :: phi
   integer, dimension(:) :: iparm
   real(dp), dimension(:) :: fparm

   real(dp) :: distR(3), deltaQ, uhatm, sh_pot, lng_pot
   real(dp) :: basis(3,3), recbasis(3,3), vol, alpha, tol
   integer :: i,j,k,l, nx,ny,nz, atom, nsh, ml, stepx, stepy, stepz
   real(dp) :: xi, yj, zk, dlx, dly, dlz

   dlx = (fparm(2) - fparm(1))/(iparm(14) - 1) 
   dly = (fparm(4) - fparm(3))/(iparm(15) - 1) 
   dlz = (fparm(6) - fparm(5))/(iparm(16) - 1)

   basis(1,1) = fparm(2)-fparm(1)
   basis(1,2) = 0.0_dp
   basis(1,3) = 0.0_dp
   basis(2,1) = 0.0_dp
   basis(2,2) = fparm(4)-fparm(3)
   basis(2,3) = 0.0_dp  
   basis(3,1) = 0.0_dp
   basis(3,2) = 0.0_dp
   basis(3,3) = fparm(6)-fparm(5)

   CALL REZVOL(basis,recbasis,vol)     
   ! choose good convergence parameter alpha
   alpha = getalpha(basis)

   nx=iparm(14)
   ny=iparm(15)
   nz=iparm(16)

   ml = minloc( (/ny*nz,nx*nz,nx*ny/), 1 )

   select case(ml)   
   case(1)
      stepx = nx - 1  
      stepy = 1 
      stepz = 1
   case(2)
      stepx = 1  
      stepy = ny - 1
      stepz = 1       
   case(3)
      stepx = 1 
      stepy = 1
      stepz = nz - 1      
   end select
   
   tol = 1.0d-5

   do k = 1, nz, stepz
     zk = fparm(1) + (k-1)*dlz
     do j = 1, ny, stepy
       yj = fparm(1) + (k-1)*dly
       do i = 1, nx, stepx
         xi = fparm(1) + (k-1)*dlx
         do atom = iatm(1), iatm(2)

            distR(1) = xi - x(1,atom)
            distR(2) = yj - x(2,atom)
            distR(3) = zk - x(3,atom)

            nsh = lmax(izp(atom))+1 
            
            ! Compute L-independent part:
            call long_pot(distR,basis,recbasis,alpha,vol,tol,lng_pot)
            ! total atomic charge
            deltaQ = sum(dQmat(1:nsh,atom) )
            phi(i,j,k) = phi(i,j,k) + deltaQ*lng_pot
            
            ! compute L-dependent part:
            do l = 1, nsh
               deltaQ = dQmat(l,atom) 
               uhatm = uhubb(l,izp(atom)) 
               
               call short_pot(distR,basis,uhatm,deltaQ,tol,sh_pot)
               
               !OMP CRITICAL
               phi(i,j,k) = phi(i,j,k) - sh_pot
               !OMP END CRITICAL
            end do
                             
         end do
       end do
     end do
   end do


 end subroutine set_phi_periodic



 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 ! CALL - BACK functions used by libmesh poisson
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 real(dp) function rho(xx1,yy1,zz1)
   
   real(dp) :: xx1, yy1, zz1
   real(dp) :: xx, yy, zz
   real(dp) :: deltaR, tau 
   integer :: atom, nsh, l
   
   
   xx=xx1+PoissBounds(1,1)
   yy=yy1+PoissBounds(2,1)
   zz=zz1+PoissBounds(3,1)
   
   rho = 0.d0
   if (.not.period) then
      
      do atom = 1,natoms   
         
         deltaR = sqrt((xx - x(1,atom))**2 + (yy - x(2,atom))**2 + (zz - x(3,atom))**2) 
         if (deltaR.gt.deltaR_max) then
            cycle
         else
         
           nsh = lmax(izp(atom))+1
           do l = 1, nsh                      
             tau = 3.2d0*uhubb(l,izp(atom))
            
             rho = rho - 0.5d0* dQmat(l,atom) * (tau)**3 * exp(-tau*deltaR)
             ! written this way rho contains -4*pi factor 
           end do
         end if
         
      end do
      
   else !Periodic structure 
      
      write(stdOut,*) 'periodic poisson not yet implemented'
      
   end if
   
   rho=-rho/4.d0/Pi  !now rho is really charge density
   
 end function rho
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 real(dp) function n_alpha(atom,xx1,yy1,zz1)
        
   integer :: atom
   real(dp) :: xx1,yy1,zz1
   
   real(dp) :: xx,yy,zz,g
   real(dp) :: deltaR
   integer :: l, nsh
   
   xx=xx1+PoissBounds(1,1)
   yy=yy1+PoissBounds(2,1)
   zz=zz1+PoissBounds(3,1)
   
   n_alpha=0.d0
   
   deltaR = sqrt((xx - x(1,atom))**2 + (yy - x(2,atom))**2 + (zz - x(3,atom))**2) 
   if (deltaR.lt.deltaR_max) then
      
      nsh = lmax(izp(atom))+1
      do l = 1, nsh
         g = 3.2d0*uhubb(l,izp(atom))             
         n_alpha  = (g**3)*exp(-g*deltaR)
         n_alpha  = n_alpha/8.d0/Pi     
      enddo
   end if
   
 end function n_alpha
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function booltoint(bool)
   logical bool
   integer booltoint

   booltoint = 0
   if (bool) booltoint = 1

 end function booltoint

end module poisson



!!$   subroutine libmesh_drv(ndim,V_orb,V_atm)
!!$
!!$      integer :: ndim
!!$      real(dp) :: V_orb(*),V_atm(*)
!!$
!!$      !local variables:
!!$      integer :: i,nn,orb,atm,m,periodic(3)
!!$      real(dp) :: Lx,Ly,Lz,bc_pot_value(2)
!!$      real(dp), ALLOCATABLE, DIMENSION(:) :: X_atm,Y_atm,Z_atm
!!$    
!!$      nn=natoms
!!$
!!$      V_orb(1:ndim)=0.d0
!!$
!!$      call log_gallocate(X_atm,nn)
!!$      call log_gallocate(Y_atm,nn)
!!$      call log_gallocate(Z_atm,nn)
!!$
!!$      do i=1,nn
!!$         X_atm(i)=x(1,i)-PoissBounds(1,1);
!!$         Y_atm(i)=x(2,i)-PoissBounds(2,1);
!!$         Z_atm(i)=x(3,i)-PoissBounds(3,1);
!!$      enddo
!!$
!!$      Lx=PoissBox(1,1)
!!$      Ly=PoissBox(2,2)
!!$      Lz=PoissBox(3,3)
!!$
!!$      bc_pot_value(1)=mu(1)
!!$      bc_pot_value(2)=mu(2)
!!$
!!$      do i=1,3
!!$         if(period_dir(i)) then 
!!$            periodic(i)=1; 
!!$         else 
!!$            periodic(i)=0; 
!!$         endif
!!$      enddo
!!$
!!$      if(id0.and.verbose.gt.VBT) call set_clock
!!$      
!!$      call libmesh_entry(nn, ndim, ind, Lx, Ly, Lz, X_atm, Y_atm, Z_atm, &
!!$                         bc_pot_value,periodic,V_orb)
!!$
!!$      if(id0.and.verbose.gt.VBT) call write_message('libmesh done')      
!!$      if(id0.and.verbose.gt.VBT) call write_clock(stdOut)  
!!$
!!$      call log_gdeallocate(X_atm)
!!$      call log_gdeallocate(Y_atm)
!!$      call log_gdeallocate(Z_atm)
!!$
!!$      ! Add up all orbital shifts to get atom shift
!!$      V_atm(1:nn) = 0.d0
!!$      do atm = iatm(1),iatm(2)
!!$         do orb = ind(atm)+1,ind(atm+1) 
!!$            V_atm(atm) = V_atm(atm) + V_orb(orb)
!!$         end do
!!$         V_atm(atm) = V_atm(atm)/(ind(atm+1)-ind(atm))
!!$      end do
!!$      
!!$      ! biasshift here so it is like ExtShift 
!!$      do m = 1, ncont
!!$         do atm=iatc(3,m),iatc(2,m)
!!$            do orb = ind(atm) + 1,ind(atm + 1) 
!!$               V_orb(orb) = -mu(m)
!!$               V_atm(atm) = V_atm(atm) + V_orb(orb)
!!$            end do
!!$            V_atm(atm) = V_atm(atm)/(ind(atm+1)-ind(atm))
!!$         end do
!!$      end do
!!$      
!!$      return
!!$
!!$   end subroutine libmesh_drv
