!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *  
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program * 
!  under the LGPL licence.                                                *
!**************************************************************************
module parcheck

  use gprecision
  use parameters
  use structure, only : natoms, x, boxsiz, period, period_dir
  use mpi_poisson, only : id0, numprocs
  use std_io

implicit none
private

 public :: check_poisson_box, check_contacts, check_localbc
 public :: check_parameters, write_parameters, check_biasdir

contains
 ! ---------------------------------------------------------------------------
 ! Perform parameters checks
 ! ---------------------------------------------------------------------------
 subroutine check_poisson_box()

   integer i
   

   if(period) then
      !cheks that the latt vect are directed along x,y,z
      !if (boxsiz(1,2).ne.0.d0 .or. boxsiz(1,3).ne. 0.d0 .or. boxsiz(2,3).ne. 0.d0) then
      !   stop 'ERROR: Supercell box is not compatible with Poisson solver'
      !end if   
      !
      if (.not.FoundBox) then
         PoissBox(:,:)=boxsiz(:,:)
         if (verbose > 30) then
            write(stdOut,*) 'Box for Poisson not Found: Set equal to supercell box' 
            do i=1,3
               write(stdOut,'(a,i1,a,f20.10)') " L(",i,")",boxsiz(i,i)
            end do
         end if
      else
         do i = 1, 3
            if (period_dir(i)) then
               PoissBox(i,i) = boxsiz(i,i)
            endif
         enddo
      end if
   else
      if (.not.FoundBox) then
         if (verbose > 30) then
            write(stdOut,*) 'Box for Poisson not Found'
         end if
         stop 'ERROR: No way to build box for Poisson'
      end if
   end if


   !if(period.and.DoGate) then
   !   if (verbose > 30) then
   !      if (id0) write(stdOut,*) 'Periodic system is not compatible with Gate'
   !      if (id0) write(stdOut,*) 'Periodicity set to F'
   !   end if
   !   period=.false.
   !endif

 end subroutine check_poisson_box

 
 subroutine check_biasdir()
  
   integer i,m

   if (.not.cluster) then
     !-OLD Bias,BiasDir compatibility -----------------
     if (ncont.eq.2) then
       !if(mu(ni(1)).eq.0.d0) mu(ni(1))=0.d0
       !if(mu(nf(1)).eq.0.d0) mu(nf(1))=bias
       if (contdir(1).ne.contdir(2)) then
         if (localBC(1).eq.0) then
           write(stdOut,*) 'ERROR: local BC should be used when &
               & contacts are in different directions'
           stop        
         endif
         if(DoCilGate) then
           write(stdOut,*) 'ERROR: contacts must be in the same &
               & direction when using cylindrical gate'
           stop
         endif

       end if
     endif

     ! ------------------------------------------------
     ! Find oriented direction of contact m 
     do m = 1, ncont
       i = contdir(m)
       if ( x(i,iatc(1,m)).lt.x(i,1) ) then
         contdir(m) = -contdir(m)
       end if
     end do
   endif
   !---------------------------------------------------
   ! Find periodic directions (where there are no contacts)
   if (period) then
      period_dir(:)=.true.
      do i=1,ncont
         period_dir(abs(contdir(i)))=.false.
      enddo
      if (DoGate) period_dir(gatedir)=.false.
      do i=1,3
         if (overrideBC(2*i).ne.0) then
            period_dir(i)=.false.
         endif
      enddo
   else
      period_dir(:)=.false.
   end if
   
   if (period_dir(3) .and. numprocs>1) then
     write(stdOut,*) 'ERROR: periodicity along z is incompatible with &
                & grid parallelization strategy'
     stop
   end if

 end subroutine check_biasdir

 subroutine check_parameters()

   use gconstants, only : pi, Kb


   if(OxLength.lt.GateLength_l) OxLength=GateLength_l
   if(Rmin_ins.lt.Rmin_gate) Rmin_ins=Rmin_gate
   if(dr_eps.lt.0.5d0) dr_eps = 0.5d0


   ! Check nPoles ----------------------------------------
   ! if (Temp.gt.0.and.nPoles.eq.0) nPoles=3
   ! LmbMax=2.d0*pi*kb*Temp*nPoles

 end subroutine check_parameters

   !-------WRITE DOWN FEEDBACK ABOUT INPUT FILE -------------------!
 subroutine write_parameters()

   integer i
   
   if (id0.and.verbose.gt.40) then
      
      !if (DoTransport.or.DoGreenDens) then
      !   write(stdOut,*) 'Conversion factor Ang/a.u.',a_u
      !   write(stdOut,*) 'Conversion factor eV/a.u.',hartree
      !   write(stdOut,*) "Delta for Green's function=",delta*hartree 
      !end if
 
      !if (DoGreenDens) then
      !   write(stdOut,*) 'Contour integration parameters:'
      !   write(stdOut,'(a4,4(i4))') 'Np=',Np(1),Np(2),Np(3)
      !   write(stdOut,*) 'nPoles=',nPoles,' LmbMax=',LmbMax*hartree
      !   write(stdOut,*) 'N_omega=',N_omega
      !   write(stdOut,*) "ReadOld Surface Green=",Readold
      !end if 

      !if (DoTransport) then
      !   write(stdOut,*) 'Temperature for electronic distribution=',Temp
      !endif

      if (DoPoisson) then
         write(stdOut,'(a,3(f10.4),a)') ' Input PoissonBox=',PoissBox(1,1)*a_u,PoissBox(2,2)*a_u,&
              PoissBox(3,3)*a_u, '  A'
         write(stdOut,*) 'PoissAcc=',PoissAcc
         if(initPot) then
            write(stdOut,*) 'Bulk Boundary Potential:    Yes'
         else
            write(stdOut,*) 'Bulk Boundary Potential:    No'
         endif
   
         write(stdOut,*) 'Atomic cutoff radius=', deltaR_max*a_u,'A'
         
         if (DoGate) then
            write(stdOut,*) 'Gate: Planar'
            write(stdOut,*) 'Gate bias=',gate*hartree,'V'
            write(stdOut,*) 'Gate length l=',GateLength_l*a_u,'A'
            write(stdOut,*) 'Gate length t=',GateLength_t*a_u,'A'
            write(stdOut,*) 'Gate distance=',Rmin_Gate*a_u,'A'           
         endif
         if (DoCilGate) then
            write(stdOut,*) 'Gate: Cylindrical'
            write(stdOut,*) 'Gate bias=',gate*hartree,'V'
            write(stdOut,*) 'Gate length=',GateLength_l*a_u,'A'
            write(stdOut,*) 'Oxide length=',OxLength*a_u,'A'
            write(stdOut,*) 'Inner gate radius=',Rmin_Gate*a_u,'A'
            write(stdOut,*) 'Inner oxide radius=',Rmin_Ins*a_u,'A'    
            write(stdOut,*) 'Dielectric constant of gate insulator=',eps_r
            write(stdOut,*) 'Smoothing of eps_r=',(eps_r-1.d0)/(dr_eps*a_u)
         end if
         if (any(localBC.gt.0)) then
            do i=1,ncont
               if (localBC(i).eq.1) write(stdOut,*) 'Local Boundary Conditions= Circular'
               if (localBC(i).eq.2) write(stdOut,*) 'Local Boundary Conditions= Squared'
               write(stdOut,'(a9,i2,a2,f8.3,a1)') ' dR_cont(',i,')=',dR_cont(i)*a_u,'A'
            enddo
         endif
         write(stdOut,*)
      endif

      
   end if

 end subroutine write_parameters

 subroutine check_localbc()
   
   integer :: m, err

   err = 0
   mixed = .false.

   if (any(localBC.gt.0)) then
      do m=1,ncont
         PoissBC(m)=1 ! Sets Mixed BC on this contact
         select case(contdir(m))   
         case(-1)
            if (overrideBC(1).eq.2) err = 1
            if (overrideBC(1).eq.1) then
               err = 2
               overrideBC(1) = 2
            endif
            mixed(1) = .true.
         case(1) 
            if (overrideBC(2).eq.2) err = 1
            if (overrideBC(2).eq.1) then
               err = 2
               overrideBC(2) = 2
            endif   
            mixed(2) = .true.
         case(-2)
            if (overrideBC(3).eq.2) err = 1
            if (overrideBC(3).eq.1) then
               err = 2
               overrideBC(3) = 2
            endif
            mixed(3) = .true.
         case(2)           
            if (overrideBC(4).eq.2) err = 1
            if (overrideBC(4).eq.1) then
               err = 2
               overrideBC(4) = 2
            endif
            mixed(4) = .true.
         case(-3)
            if (overrideBC(5).eq.2) err = 1
            if (overrideBC(5).eq.1) then
               err = 2
               overrideBC(5) = 2
            endif
            mixed(5) = .true.
         case(3)
            if (overrideBC(6).eq.2) err = 1
            if (overrideBC(6).eq.1) then
               err = 2
               overrideBC(6) = 2
            endif
            mixed(6) = .true.
         end select

      enddo
   endif
   if (err.eq.1 .and. id0) then
      write(stdOut,*)
      write(stdOut,*) 'ERROR: BoundaryRegion{} sets a Dirichlet BC'
      write(stdOut,*) '       not compatible with overrideDefaultBC=Neumann'
   endif
   if (err.eq.2 .and. id0) then
      write(stdOut,*)
      write(stdOut,*) 'WARNING: '
      write(stdOut,*) 'BoundaryRegion{} sets a mixed Neumann/Dirichlet BC'
      write(stdOut,*) 'User setting OverrideDefaultBC = Dirichlet'
      write(stdOut,*) 'has been disregarded !'
      write(stdOut,*)
   endif
   

 end subroutine check_localbc
   !--- WRITE INFOS ABOUT THE CONTACT STRUCTURES ---------------
 subroutine check_contacts()

   integer i,ncdim_max
    
   if(cluster) then
      if (id0) then
          write(stdOut,'(1x,a)',advance='NO') 'System Type: UNCONTACTED '
          if(period) then
             write(stdOut,'(a)',advance='NO') 'PERIODIC '
          else
             write(stdOut,'(a)',advance='NO') 'CLUSTER '
          endif
          write(stdOut,*) 'STRUCTURE'
          write(stdOut,'(1x,a,I3)') 'periodicity direction:',contdir(1)
      endif  
   endif


   if(id0.and.verbose.gt.50) then      
      write(stdOut,'(a)') 'CENTRAL REGION'
      write(stdOut,'(1x,a,2I6)') 'Atom start - end = ',iatm(1), iatm(2)
      write(stdOut,*)
   endif

   
   if(.not.cluster) then

      do i=1,ncont

         ! CHECK CONTACTS  
         if(iatc(3,i).lt.iatc(1,i)) then 
            iatc(3,i)=iatc(1,i)
         endif

         ! Initialize end of contact blocks ------------------
         !ncdim(i)=ind(iatc(2,i)+1)-ind(iatc(1,i))
         !mbound_end(i)=ind(iatc(3,i))-ind(iatc(1,i))

         !if (iatc(3,i).eq.iatc(1,i)) then
         !   mbound_end(i)=0
         !else
         !   mbound_end(i)=ind(iatc(3,i)+1)-ind(iatc(1,i))   
         !endif
         ! ---------------------------------------------------

         write(stdOut,'(a,I3)') 'CONTACT #',i
         write(stdOut,'(1x,a,2I6)') 'Atom start - end = ',iatc(3,i), iatc(2,i)
         write(stdOut,'(1x,a,I3)') 'direction:',contdir(i)
         write(stdOut,*) 'Fermi Level=',Efermi(i)*hartree,'eV'
         write(stdOut,*) 'mu=',mu(i)*hartree,'V'
         write(stdOut,*) 

      end do !ncont
   endif !cluster


   ncdim_max=maxval(ncdim)

   !-----CHECK IF SYSTEM STRUCTURE IS CONSISTENT-------!
   if (.not.cluster) then
      do i=1,ncont
         if (iatc(1,i).lt.iatm(2)) then 
            write(stdOut,*) 'ERROR: The contacts MUST be defined after the scattering region'
            stop
         endif
      enddo
   endif
   if ((iatm(2)-iatm(1)+1).gt.natoms) then
      write(stdOut,*) 'ERROR: The number of atoms in the scattering region is higer'
      write(stdOut,*) '       than the total number of atoms'
      stop
   endif
   if (DoGate) then
      if (gatedir.ne.2) then 
         write(stdOut,*) "ERROR: gate direction must be along y"
         stop
      endif
      if(any(abs(contdir(:)).eq.gatedir)) then
         write(stdOut,*) "ERROR: gate direction along contacts!?"
         stop
      endif
   endif

   !    Checks if the number of mouvable atoms is set correctly

 end subroutine check_contacts

 end module parcheck
