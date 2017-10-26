module parcheck

  use gprecision
  use parameters
  use structure, only : natoms, nbeweg, x, boxsiz, period, period_dir
  use mpi_poisson, only : id0, numprocs

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
      !   STOP 'ERROR: Supercell box is not compatible with Poisson solver'
      !end if   
      !
      if (.not.FoundBox) then
         PoissBox(:,:)=boxsiz(:,:)
         if (id0.and.verbose > 30) then
            if (id0) write(*,*) 'Box for Poisson not Found: Set equal to supercell box' 
            do i=1,3
               if (id0) write(*,'(a,i1,a,f20.10)') " L(",i,")",boxsiz(i,i)
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
         if (id0 .and. verbose > 30) then
            write(*,*) 'Box for Poisson not Found'
         end if
         STOP 'ERROR: No way to build box for Poisson'
      end if
   end if


   !if(period.and.DoGate) then
   !   if (verbose > 30) then
   !      if (id0) write(*,*) 'Periodic system is not compatible with Gate'
   !      if (id0) write(*,*) 'Periodicity set to F'
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
       if ( contdir(1).ne.contdir(2) ) then
         if(localBC(1).eq.0) then
           if (id0) write(*,*) 'ERROR: local BC should be used when &
               & contacts are in different directions'
           STOP        
         endif
         if(DoCilGate) then
           if(id0) write(*,*) 'ERROR: contacts must be in the same &
               & direction when using cylindrical gate'
           STOP
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
     if(id0) write(*,*) 'ERROR: periodicity along z is incompatible with &
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
      !   write(*,*) 'Conversion factor Ang/a.u.',a_u
      !   write(*,*) 'Conversion factor eV/a.u.',hartree
      !   write(*,*) "Delta for Green's function=",delta*hartree 
      !end if
 
      !if (DoGreenDens) then
      !   write(*,*) 'Contour integration parameters:'
      !   write(*,'(a4,4(i4))') 'Np=',Np(1),Np(2),Np(3)
      !   write(*,*) 'nPoles=',nPoles,' LmbMax=',LmbMax*hartree
      !   write(*,*) 'N_omega=',N_omega
      !   write(*,*) "ReadOld Surface Green=",Readold
      !end if 

      !if (DoTransport) then
      !   write(*,*) 'Temperature for electronic distribution=',Temp
      !endif

      if (DoPoisson) then
         write(*,'(a,3(f10.4),a)') ' Input PoissonBox=',PoissBox(1,1)*a_u,PoissBox(2,2)*a_u,&
              PoissBox(3,3)*a_u, '  A'
         write(*,*) 'PoissAcc=',PoissAcc
         if(initPot) then
            write(*,*) 'Bulk Boundary Potential:    Yes'
         else
            write(*,*) 'Bulk Boundary Potential:    No'
         endif
   
         write(*,*) 'Atomic cutoff radius=', deltaR_max*a_u,'A'
         
         if(DoGate) write(*,*) 'Gate: Planar'
         if(DoCilGate) write(*,*) 'Gate: Cylindrical'
         
         if (DoGate) then
            write(*,*) 'Gate bias=',gate*hartree,'V'
            write(*,*) 'Gate length l=',GateLength_l*a_u,'A'
            write(*,*) 'Gate length t=',GateLength_t*a_u,'A'
            write(*,*) 'Gate distance=',Rmin_Gate*a_u,'A'           
         endif
         if (DoCilGate) then
            write(*,*) 'Gate bias=',gate*hartree,'V'
            write(*,*) 'Gate length=',GateLength_l*a_u,'A'
            write(*,*) 'Oxide length=',OxLength*a_u,'A'
            write(*,*) 'Inner gate radius=',Rmin_Gate*a_u,'A'
            write(*,*) 'Inner oxide radius=',Rmin_Ins*a_u,'A'    
            write(*,*) 'Dielectric constant of gate insulator=',eps_r
            write(*,*) 'Smoothing of eps_r=',(eps_r-1.d0)/(dr_eps*a_u)
         end if
         if (any(localBC.gt.0)) then
            do i=1,ncont
               if (localBC(i).eq.1) write(*,*) 'Local Boundary Conditions= Circular'
               if (localBC(i).eq.2) write(*,*) 'Local Boundary Conditions= Squared'
               write(*,'(a9,i2,a2,f8.3,a1)') ' dR_cont(',i,')=',dR_cont(i)*a_u,'A'
            enddo
         endif
         write(*,*)
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
      write(*,*)
      write(*,*) 'ERROR: BoundaryRegion{} sets a Dirichlet BC'
      write(*,*) '       not compatible with overrideDefaultBC=Neumann'
   endif
   if (err.eq.2 .and. id0) then
      write(*,*)
      write(*,*) 'WARNING: '
      write(*,*) 'BoundaryRegion{} sets a mixed Neumann/Dirichlet BC'
      write(*,*) 'User setting OverrideDefaultBC = Dirichlet'
      write(*,*) 'has been disregarded !'
      write(*,*)
   endif
   

 end subroutine check_localbc
   !--- WRITE INFOS ABOUT THE CONTACT STRUCTURES ---------------
 subroutine check_contacts()

   integer i,ncdim_max
    
   if(cluster) then
      if (id0) then
          write(*,'(1x,a)',advance='NO') 'System Type: UNCONTACTED '
          if(period) then
             write(*,'(a)',advance='NO') 'PERIODIC '
          else
             write(*,'(a)',advance='NO') 'CLUSTER '
          endif
          write(*,*) 'STRUCTURE'
          write(*,'(1x,a,I3)') 'periodicity direction:',contdir(1)
          write(*,*) 'Fermi Level=',Efermi(1)*hartree,'eV'
      endif  
   endif


   if(id0.and.verbose.gt.50) then      
      write(*,'(a)') 'CENTRAL REGION'
      write(*,'(1x,a,2I6)') 'Atom start - end = ',iatm(1), iatm(2)
      write(*,*)
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

         write(*,*) '(poiss_init) CONTACT INFO #',i
         write(*,'(a,I3)') 'CONTACT #',i
         write(*,'(1x,a,2I6)') 'Atom start - end = ',iatc(3,i), iatc(2,i)
         write(*,'(1x,a,I3)') 'direction:',contdir(i)
         write(*,*) 'Fermi Level=',Efermi(i)*hartree,'eV'
         write(*,*) 'mu=',mu(i)*hartree,'V'
         write(*,*) 

      end do !ncont
   endif !cluster


   ncdim_max=maxval(ncdim)

   !-----CHECK IF SYSTEM STRUCTURE IS CONSISTENT-------!
   if (.not.cluster) then
      do i=1,ncont
         if (iatc(1,i).lt.iatm(2)) then 
            if(id0) write(*,*) 'ERROR: The contacts MUST be defined after the scattering region'
            STOP
         endif
      enddo
   endif
   if ((iatm(2)-iatm(1)+1).gt.natoms) then
      if(id0) write(*,*) 'ERROR: The number of atoms in the scattering region is higer'
      if(id0) write(*,*) '       than the total number of atoms'
      STOP
   endif
   if(DoGate) then
      if(gatedir.ne.2) then 
         if (id0) write(*,*) "ERROR: gate direction must be along y"
         STOP
      endif
      if(any(abs(contdir(:)).eq.gatedir)) then
         if (id0) write(*,*) "ERROR: gate direction along contacts!?"
         STOP
      endif
   endif

   !    Checks if the number of mouvable atoms is set correctly
   if(nbeweg > natoms) nbeweg=natoms 
   if(nbeweg > (iatm(2)-iatm(1)+1)) nbeweg=iatm(2)-iatm(1)+1

 end subroutine check_contacts

 end module parcheck
