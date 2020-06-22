!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program * 
!  under the LGPL licence.                                                *
!**************************************************************************

#:include "error.fypp"

module parcheck

  use dftbp_accuracy, only : lc, dp
  use dftbp_constants
  use dftbp_message, only : warning
  use parameters
  use structure, only : natoms, x, boxsiz, period, period_dir
  use mpi_poisson, only : id0, numprocs
  use dftbp_globalenv, only : stdOut

implicit none
private

 public :: check_poisson_box, check_contacts, check_localbc
 public :: check_parameters, write_parameters, check_biasdir

 !> Verbosity threashold
 integer, parameter :: VBT=30

contains
 ! ---------------------------------------------------------------------------
 ! Perform parameters checks
 ! ---------------------------------------------------------------------------
 subroutine check_poisson_box(iErr)

   integer, intent(out), optional :: iErr

   integer i

   if (present(iErr)) then
     iErr = 0
   end if

   if(period) then
      !checks that the latt vect are directed along x,y,z
      !if (boxsiz(1,2).ne.0.d0 .or. boxsiz(1,3).ne. 0.d0 .or. boxsiz(2,3).ne. 0.d0) then
      !   call error('ERROR: Supercell box is not compatible with Poisson solver')
      !end if   
      !
      if (.not.FoundBox) then
         PoissBox(:,:)=boxsiz(:,:)
         if (verbose > VBT) then
            call warning('Box for Poisson not Found: Setting equal to supercell box')
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
        if (verbose > VBT) then
          call warning('Box for Poisson not Found')
        end if
        @:ERROR_HANDLING(iErr, -1, 'No way to build box for Poisson')
      end if
   end if


   !if(period.and.DoGate) then
   !   if (verbose > VBT) then
   !      if (id0) write(stdOut,*) 'Periodic system is not compatible with Gate'
   !      if (id0) write(stdOut,*) 'Periodicity set to F'
   !   end if
   !   period=.false.
   !endif

 end subroutine check_poisson_box

 
 subroutine check_biasdir(iErr)

   integer, intent(out), optional :: iErr

   integer i,m

   if (present(iErr)) then
     iErr = 0
   end if

   if (.not.cluster) then
     !-OLD Bias,BiasDir compatibility -----------------
     if (ncont.eq.2) then
       !if(mu(ni(1)).eq.0.d0) mu(ni(1))=0.d0
       !if(mu(nf(1)).eq.0.d0) mu(nf(1))=bias
       if (contdir(1).ne.contdir(2)) then
         if (localBC(1).eq.0) then
           @:ERROR_HANDLING(iErr, -1,&
               & 'Local BC should be used when contacts are in different directions')
         endif
         if(DoCilGate) then
           @:ERROR_HANDLING(iErr, -2,&
               & 'Contacts must be in the same direction when using cylindrical gate')
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
     @:ERROR_HANDLING(iErr, -3,&
         & 'Periodicity along z is incompatible with grid parallelization strategy')
   end if

 end subroutine check_biasdir

 subroutine check_parameters()

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
      !   write(stdOut,*) 'Conversion factor Ang/a.u.',Bohr__AA
      !   write(stdOut,*) 'Conversion factor eV/a.u.',hartree__eV
      !   write(stdOut,*) "Delta for Green's function=",delta*hartree__eV
      !end if
 
      !if (DoGreenDens) then
      !   write(stdOut,*) 'Contour integration parameters:'
      !   write(stdOut,'(a4,4(i4))') 'Np=',Np(1),Np(2),Np(3)
      !   write(stdOut,*) 'nPoles=',nPoles,' LmbMax=',LmbMax*hartree__eV
      !   write(stdOut,*) 'N_omega=',N_omega
      !   write(stdOut,*) "ReadOld Surface Green=",Readold
      !end if 

      !if (DoTransport) then
      !   write(stdOut,*) 'Temperature for electronic distribution=',Temp
      !endif

      if (DoPoisson) then
        write(stdOut,'(a,3(f10.4),a)') ' Input PoissonBox=',PoissBox(1,1)*Bohr__AA,&
            & PoissBox(2,2)*Bohr__AA, PoissBox(3,3)*Bohr__AA, '  A'
         write(stdOut,*) 'PoissAcc=',PoissAcc
         if(initPot) then
            write(stdOut,*) 'Bulk Boundary Potential:    Yes'
         else
            write(stdOut,*) 'Bulk Boundary Potential:    No'
         endif
   
         write(stdOut,*) 'Atomic cutoff radius=', deltaR_max*Bohr__AA,'A'
         
         if (DoGate) then
            write(stdOut,*) 'Gate: Planar'
            write(stdOut,*) 'Gate bias=',gate*hartree__eV,'V'
            write(stdOut,*) 'Gate length l=',GateLength_l*Bohr__AA,'A'
            write(stdOut,*) 'Gate length t=',GateLength_t*Bohr__AA,'A'
            write(stdOut,*) 'Gate distance=',Rmin_Gate*Bohr__AA,'A'
         endif
         if (DoCilGate) then
            write(stdOut,*) 'Gate: Cylindrical'
            write(stdOut,*) 'Gate bias=',gate*hartree__eV,'V'
            write(stdOut,*) 'Gate length=',GateLength_l*Bohr__AA,'A'
            write(stdOut,*) 'Oxide length=',OxLength*Bohr__AA,'A'
            write(stdOut,*) 'Inner gate radius=',Rmin_Gate*Bohr__AA,'A'
            write(stdOut,*) 'Inner oxide radius=',Rmin_Ins*Bohr__AA,'A'
            write(stdOut,*) 'Dielectric constant of gate insulator=',eps_r
            write(stdOut,*) 'Smoothing of eps_r=',(eps_r-1.d0)/(dr_eps*Bohr__AA)
         end if
         if (any(localBC.gt.0)) then
            do i=1,ncont
               if (localBC(i).eq.1) write(stdOut,*) 'Local Boundary Conditions= Circular'
               if (localBC(i).eq.2) write(stdOut,*) 'Local Boundary Conditions= Squared'
               write(stdOut,'(a9,i2,a2,f8.3,a1)') ' dR_cont(',i,')=',dR_cont(i)*Bohr__AA,'A'
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
 subroutine check_contacts(iErr)

   integer, intent(out), optional :: iErr

   integer i,ncdim_max

   if (present(iErr)) then
     iErr = 0
   end if

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
         write(stdOut,*) 'Fermi Level=',Efermi(i)*hartree__eV,'eV'
         write(stdOut,*) 'mu=',mu(i)*hartree__eV,'V'
         write(stdOut,*) 

      end do !ncont
   endif !cluster


   ncdim_max=maxval(ncdim)

   !-----CHECK IF SYSTEM STRUCTURE IS CONSISTENT-------!
   if (.not.cluster) then
     do i=1,ncont
       if (iatc(1,i).lt.iatm(2)) then
         @:ERROR_HANDLING(iErr, -1,&
             & 'The contacts MUST be defined after the scattering region')
       end if
     enddo
   endif
   if ((iatm(2)-iatm(1)+1).gt.natoms) then
     @:ERROR_HANDLING(iErr, -2, 'The number of atoms in the scattering region is higer than&
         & the total number of atoms')
   endif
   if (DoGate) then
     if (gatedir.ne.2) then
       @:ERROR_HANDLING(iErr, -3, 'Gate direction must be along y')
     end if
     if(any(abs(contdir(:)).eq.gatedir)) then
       @:ERROR_HANDLING(iErr, -4, 'Gate direction along contacts!?')
     end if
   endif

 end subroutine check_contacts


end module parcheck
