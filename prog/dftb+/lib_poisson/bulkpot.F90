!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *  
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program * 
!  under the LGPL licence.                                                *
!**************************************************************************

#:include "error.fypp"

Module bulkpot
  
 use dftbp_accuracy, only : dp
 use dftbp_constants
 use gallocation
 use parameters
 use structure
 use mpi_poisson
 use gewald

 implicit none

 private

 public :: super_array, create_super_array,destroy_super_array
 public :: create_phi_bulk,destroy_phi_bulk,readbulk_pot,compbulk_pot
 public :: save_bulkpot, write_super_array

 type super_array  
   integer :: a,b,c
   integer :: iparm(23)
   real(kind=dp) :: fparm(8)
   real(kind=dp) :: dla,dlb,dlc
   integer :: ibsize
   integer :: natm_PL
   real(kind=dp) :: L_PL
   real(kind=dp), DIMENSION (:,:,:), ALLOCATABLE :: val 
   real(kind=dp), DIMENSION (:,:,:), ALLOCATABLE :: rhs
   logical :: doEwald
 end type  super_array



contains

 !%--------------------------------------------------------------------------
 subroutine create_super_array(SA,na,nb,nc)

   type(super_array) :: SA
   integer :: na,nb,nc

   call log_gallocate(SA%val,na,nb,nc)

   SA%ibsize=na*nb*nc

 end subroutine create_super_array
 
 !%--------------------------------------------------------------------------
 subroutine destroy_super_array(SA)

   type(super_array) :: SA

   call log_gdeallocate(SA%val)

 end subroutine destroy_super_array
 !%--------------------------------------------------------------------------
   

 subroutine write_super_array(SA)

   type(super_array) :: SA

   if (id0) then
     write(*,*) SA%a,SA%b,SA%c
     write(*,*) SA%dla,SA%dlb,SA%dlc
     write(*,*) 'size',SA%ibsize
     write(*,*) 'iparm',SA%iparm
     write(*,*) 'fparm',SA%fparm
     write(*,*) 'natm_PL',SA%natm_PL
     write(*,*) 'L_PL',SA%L_PL  
     write(*,*) 'rhs',size(SA%rhs)
     write(*,*) 'val',size(SA%val)
   endif

 end subroutine write_super_array
 
 !%--------------------------------------------------------------------------  
 subroutine create_phi_bulk(phi_bulk,iparm,dlx,dly,dlz,cont_mem)

 type(super_array) :: phi_bulk(:)
 integer :: iparm(23)
 integer :: cont_mem,na,nb,nc,m
 real(kind=dp) :: dlx,dly,dlz
 integer :: nstart, nlast, num_p, i

 cont_mem=0

 do m=1,ncont

   nstart = iatc(3,m)
   nlast  = iatc(2,m)
   phi_bulk(m)%natm_PL = (nlast-nstart+1)/2
   phi_bulk(m)%doEwald = .false.

   select case(abs(contdir(m)))    
     ! contacts are internally oriented along z.
   case(1)
     phi_bulk(m)%a = 2   !(y)
     phi_bulk(m)%b = 3   !(z)
     phi_bulk(m)%c = 1   !(x)
     phi_bulk(m)%dla = dly
     phi_bulk(m)%dlb = dlz
     phi_bulk(m)%dlc = dlx

     ! set the bulk contact periodicity in x
     phi_bulk(m)%iparm(2)= iparm(4)    !Copy from device
     phi_bulk(m)%iparm(3)= iparm(5)    !     "     "         
     phi_bulk(m)%iparm(4)= iparm(6)    !     "     "   
     phi_bulk(m)%iparm(5)= iparm(7)    !     "     "  
     phi_bulk(m)%iparm(6)= 0           ! Periodic in  c
     phi_bulk(m)%iparm(7)= 0           ! Periodic in  c

     if(overrBulkBC(3).gt.-1) phi_bulk(m)%iparm(2)=overrBulkBC(3)
     if(overrBulkBC(4).gt.-1) phi_bulk(m)%iparm(3)=overrBulkBC(4)
     if(overrBulkBC(5).gt.-1) phi_bulk(m)%iparm(4)=overrBulkBC(5)
     if(overrBulkBC(6).gt.-1) phi_bulk(m)%iparm(5)=overrBulkBC(6)
     !------------------------------------------------------
     phi_bulk(m)%iparm(8)  = iparm(9)  ! iyp 
     phi_bulk(m)%iparm(9)  = iparm(10) ! izp
     phi_bulk(m)%iparm(10) = iparm(8)  ! ixp
   
     phi_bulk(m)%iparm(11) = iparm(12) !highest grid in a  (y)
     phi_bulk(m)%iparm(12) = iparm(13) !highest grid in b  (z)

     ! define the PL periodicity     
     phi_bulk(m)%L_PL = abs(x(1,nstart+phi_bulk(m)%natm_PL)-x(1,nstart))
     ! find the right grid
     do i = 1,50 
        phi_bulk(m)%iparm(13) = i      !highest grid in c  (x)
        num_p = phi_bulk(m)%iparm(10) *(2**(i - 1)) + 1 
        if (phi_bulk(m)%L_PL/(num_p - 1).le.dmin(1)) then
           phi_bulk(m)%dlc = phi_bulk(m)%L_PL/(num_p - 1) 
           exit 
        end if
     end do
     phi_bulk(m)%iparm(14) = iparm(15) !# grid points in a (y)
     phi_bulk(m)%iparm(15) = iparm(16) !# grid points in b (z)
     phi_bulk(m)%iparm(16) = num_p     !# grid points in c (x)
     
   case(2)
     phi_bulk(m)%a = 3  !(z)
     phi_bulk(m)%b = 1  !(x)
     phi_bulk(m)%c = 2  !(y)
     phi_bulk(m)%dla = dlz
     phi_bulk(m)%dlb = dlx
     phi_bulk(m)%dlc = dly
     ! set the bulk contact periodicity in y
     phi_bulk(m)%iparm(2)= iparm(6)    !Copy from device
     phi_bulk(m)%iparm(3)= iparm(7)    !     "     "         
     phi_bulk(m)%iparm(4)= iparm(2)    !     "     "   
     phi_bulk(m)%iparm(5)= iparm(3)    !     "     "  
     phi_bulk(m)%iparm(6)= 0           ! Periodic in  c
     phi_bulk(m)%iparm(7)= 0           ! Periodic in  c

     if(overrBulkBC(5).gt.-1) phi_bulk(m)%iparm(2)=overrBulkBC(5)
     if(overrBulkBC(6).gt.-1) phi_bulk(m)%iparm(3)=overrBulkBC(6)
     if(overrBulkBC(1).gt.-1) phi_bulk(m)%iparm(4)=overrBulkBC(1)
     if(overrBulkBC(2).gt.-1) phi_bulk(m)%iparm(5)=overrBulkBC(2)
     !------------------------------------------------------
     phi_bulk(m)%iparm(8)  = iparm(10) ! izp 
     phi_bulk(m)%iparm(9)  = iparm(8)  ! ixp
     phi_bulk(m)%iparm(10) = iparm(9)  ! iyp
   
     phi_bulk(m)%iparm(11) = iparm(13) !highest grid in a  (z)
     phi_bulk(m)%iparm(12) = iparm(11) !highest grid in b  (x)

     ! define the PL periodicity    
     phi_bulk(m)%L_PL = abs(x(2,nstart+phi_bulk(m)%natm_PL)-x(2,nstart))
     ! find the right grid
     do i = 1,50 
        phi_bulk(m)%iparm(13) = i      !highest grid in c  (y)
        num_p = phi_bulk(m)%iparm(10) *(2**(i - 1)) + 1 
        if (phi_bulk(m)%L_PL/(num_p - 1).le.dmin(2)) then
           phi_bulk(m)%dlc = phi_bulk(m)%L_PL/(num_p - 1)
           exit 
        end if
     end do
     
     phi_bulk(m)%iparm(14) = iparm(16) !# grid points in a (z)
     phi_bulk(m)%iparm(15) = iparm(14) !# grid points in b (x)
     phi_bulk(m)%iparm(16) = num_p     !# grid points in c (y)
           
   case(3)
     phi_bulk(m)%a = 1   !(x)
     phi_bulk(m)%b = 2   !(y)
     phi_bulk(m)%c = 3   !(z)
     phi_bulk(m)%dla = dlx
     phi_bulk(m)%dlb = dly
     phi_bulk(m)%dlc = dlz
     ! set the bulk contact periodicity in z
     phi_bulk(m)%iparm(2)= iparm(2)    !Copy from device
     phi_bulk(m)%iparm(3)= iparm(3)    !     "     "         
     phi_bulk(m)%iparm(4)= iparm(4)    !     "     "   
     phi_bulk(m)%iparm(5)= iparm(5)    !     "     " 
     phi_bulk(m)%iparm(6)= 0           ! Periodic in  c
     phi_bulk(m)%iparm(7)= 0           ! Periodic in  c

     if(overrBulkBC(1).gt.-1) phi_bulk(m)%iparm(2)=overrBulkBC(1)
     if(overrBulkBC(2).gt.-1) phi_bulk(m)%iparm(3)=overrBulkBC(2)
     if(overrBulkBC(3).gt.-1) phi_bulk(m)%iparm(4)=overrBulkBC(3)
     if(overrBulkBC(4).gt.-1) phi_bulk(m)%iparm(5)=overrBulkBC(4) 
     !------------------------------------------------------
     phi_bulk(m)%iparm(8)  = iparm(8)  ! ixp 
     phi_bulk(m)%iparm(9)  = iparm(9)  ! iyp
     phi_bulk(m)%iparm(10) = iparm(10) ! izp

     phi_bulk(m)%iparm(11) = iparm(11) !highest grid in a  (x) 
     phi_bulk(m)%iparm(12) = iparm(12) !highest grid in b  (y)

     ! define the PL periodicity    
     phi_bulk(m)%L_PL = abs(x(3,nstart+phi_bulk(m)%natm_PL)-x(3,nstart))
     ! find the right grid
     do i = 1,50 
        phi_bulk(m)%iparm(13) = i      !highest grid in c  (z)
        num_p = phi_bulk(m)%iparm(10) *(2**(i - 1)) + 1 
        if (phi_bulk(m)%L_PL/(num_p - 1).le.dmin(3)) then
           phi_bulk(m)%dlc = phi_bulk(m)%L_PL/(num_p - 1) 
           exit 
        end if
     end do

     phi_bulk(m)%iparm(14) = iparm(14) !# grid points in a (x)
     phi_bulk(m)%iparm(15) = iparm(15) !# grid points in b (y)
     phi_bulk(m)%iparm(16) = num_p     !# grid points in c (z)

   end select

   na= phi_bulk(m)%iparm(14)               !# grid points in a
   nb= phi_bulk(m)%iparm(15)               !# grid points in b
   nc= phi_bulk(m)%iparm(16)               !# grid points in c

   ! CASE ALL PERIODIC:
   ! Choose the smallest area where compute Ewald sums 
   ! Override periodic BC with Dirichlet
   if ( all(phi_bulk(m)%iparm(2:7).eq.0) ) then
     phi_bulk(m)%doEwald=.true.
     if(nb*nc .lt. na*nc) then
        phi_bulk(m)%iparm(2)= 1           !Dirichlet in a 
        phi_bulk(m)%iparm(3)= 1           !     "     "  a
     else 
        phi_bulk(m)%iparm(4)= 1           !Dirichlet in b 
        phi_bulk(m)%iparm(5)= 1           !     "     "  b
     endif
   end if

   ! CASE PERIODIC & NEUMANN:
   ! Override with Dirichlet.
   if ( sum(phi_bulk(m)%iparm(2:7)).eq.8 ) then
        phi_bulk(m)%iparm(2)= 1           !Dirichlet in a 
        phi_bulk(m)%iparm(3)= 1           !     "     "  a
        phi_bulk(m)%iparm(4)= 1           !Dirichlet in b 
        phi_bulk(m)%iparm(5)= 1           !     "     "  b     
   endif


   !------------------------------------------------------
   phi_bulk(m)%iparm(17) = 0            ! no initial guess
   phi_bulk(m)%iparm(18) = 50 !iparm(18)    ! # of iterations
   phi_bulk(m)%iparm(19) = 0            ! Gauss-Siedel
   phi_bulk(m)%iparm(20) = 7*(na+2)*(nb+2)*(nc+2)/2
   
   call log_gallocate(phi_bulk(m)%rhs,na,nb,nc)
   
   !write(*,*) '%LOC rhs=',%LOC(phi_bulk(m)%rhs)
   cont_mem = na*nb*nc          
   
   cont_mem = cont_mem+na*nb*nc 
   
   phi_bulk(m)%ibsize = cont_mem 
   
   call log_gallocate(phi_bulk(m)%val,na,nb,nc)
   
   phi_bulk(m)%val(1:na,1:nb,1:nc)=0.d0
   
   !write(*,*) '%LOC val=',%LOC(phi_bulk(m)%val)
   
   !if(id0.and.verbose.gt.80) then
   !   write(*,*) 'Bulk Potential Contact #',m,nstart,nlast
   !   write(*,*) 'N(a)=',phi_bulk(m)%iparm(14),'dla=',phi_bulk(m)%dla*Bohr__AA
   !   write(*,*) 'N(b)=',phi_bulk(m)%iparm(15),'dlb=',phi_bulk(m)%dlb*Bohr__AA
   !   write(*,*) 'N(c)=',phi_bulk(m)%iparm(16),'dlc=',phi_bulk(m)%dlc*Bohr__AA
   !endif

  enddo

end subroutine create_phi_bulk

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine destroy_phi_bulk(phi_bulk)

  type(super_array) :: phi_bulk(:)
  integer :: m

  do m=1,ncont 
    call log_gdeallocate(phi_bulk(m)%val)
    call log_gdeallocate(phi_bulk(m)%rhs)
  enddo

end subroutine destroy_phi_bulk

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%--------------------------------------------------------------------------  
Subroutine readbulk_pot(phi_bulk, iErr)
  type(super_array) :: phi_bulk(:)

  integer, intent(out), optional :: iErr

  integer :: i,j,k,m
  character(2) :: m_id
  real(kind=dp) :: tmp_dbl

  integer :: a,b,c, fp
  logical :: lex

  if (present(iErr)) then
    iErr = 0
  end if

  do m = 1, ncont
   
    write(m_id,'(i2.2)') m
    inquire(file='contacts/BulkPot_'//m_id//'.dat',exist=lex)
    if (.not.lex) then
      @:ERROR_HANDLING(iErr, -1, 'File contacts/BulkPot_'//m_id//'.dat not found')
    else
      open(newunit=fp,file='contacts/BulkPot_'//m_id//'.dat',form='formatted')
    endif   
    
    read(fp,*) a,b,c
    
    if (a.ne.phi_bulk(m)%iparm(14) .or. &
      b.ne.phi_bulk(m)%iparm(15) .or. &
      c.ne.phi_bulk(m)%iparm(16)) then
      if(id0) write(*,*) 'Warning: incompatible BulkPot: will be recomputed'   
      ReadBulk = .false.
      close(fp)
      return
    endif
    
    do i = 1,phi_bulk(m)%iparm(14)  
       do j = 1,phi_bulk(m)%iparm(15)
          do k = 1,phi_bulk(m)%iparm(16)
    
             read(fp,*) tmp_dbl
             phi_bulk(m)%val(i,j,k) = tmp_dbl
    
          end do
       end do
    end do
    
    close(fp)

  enddo

end subroutine  readbulk_pot
!%--------------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine compbulk_pot(phi_bulk,iparm,fparm)
  type(super_array) :: phi_bulk(:)
  integer :: iparm(23)
  real(kind=dp) :: fparm(8)

  call compbulk_pot_mud(phi_bulk,iparm,fparm)

end subroutine compbulk_pot

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine compbulk_pot_ewald(phi_bulk,m)
  type(super_array) :: phi_bulk(:)
  integer :: m

  !local variables
  real(kind=dp), allocatable, dimension(:,:,:) :: phi_bulk_PAR 
  integer :: i,j,k,ibsize,atom,a,b,c,na,nb,nc,ff,stepa,stepb
  real(kind=dp) :: yj,zk,xi
  real(kind=dp) :: basis(3,3), recbasis(3,3)
  real(kind=dp) :: alpha, vol, tol
  character(2) :: m_id
  real(kind=dp) :: distR(3), deltaQ, uhatm, sh_pot, lng_pot

  integer :: npid,istart,iend, nsh,l

  ! set tolerance for convergence
  tol = 1.0d-5

  !write(*,*) " " 
  !write(*,*) 'BULK-POTENTIAL CPU =', id            

  ! Ewald sum initialization 
  ! get reciprocal lattice vectors and cell volume   
  basis(1,1) = phi_bulk(m)%fparm(2)-phi_bulk(m)%fparm(1)
  basis(1,2) = 0.0_dp
  basis(1,3) = 0.0_dp
  basis(2,1) = 0.0_dp
  basis(2,2) = phi_bulk(m)%fparm(4)-phi_bulk(m)%fparm(3)
  basis(2,3) = 0.0_dp  
  basis(3,1) = 0.0_dp
  basis(3,2) = 0.0_dp
  basis(3,3) = phi_bulk(m)%fparm(6)-phi_bulk(m)%fparm(5)
  
  CALL REZVOL(basis,recbasis,vol)     
  ! choose good convergence parameter alpha
  alpha = getalpha(basis)

  istart = iatc(3,m)
  iend = iatc(2,m)

  !npid = int( (iatc(2,m)-iatc(3,m)+1)/numprocs )
  !istart = iatc(3,m)+id*npid
  !if (id.ne.(numprocs-1)) then
  !   iend = iatc(3,m)+(id+1)*npid-1
  !else
  !   iend = iatc(2,m)
  !endif

  na=phi_bulk(m)%iparm(14)
  nb=phi_bulk(m)%iparm(15)
  nc=phi_bulk(m)%iparm(16)
  a=phi_bulk(m)%a
  b=phi_bulk(m)%b
  c=phi_bulk(m)%c
  ibsize=phi_bulk(m)%ibsize

  ! decide which area is smaller and compute ewalds only on that
  ! (Dirichlet should have been set before) 
  if(nb*nc .lt. na*nc) then
      stepa=na-1 ! Ewalds is computed on 1 and na
      stepb=1
  else
      stepa=1
      stepb=nb-1 ! Ewalds is computed on 1 and nb
  endif 

  call log_gallocate( phi_bulk_PAR,na,nb,nc)

  phi_bulk_PAR(:,:,:) = 0.d0

  !OMP PARALLEL DO  FIRSTPRIVATE(basis,tol,alpha,vol) &
  !OMP& PRIVATE(distR,deltaQ,uhatm,sh_pot,lng_pot,xi,yj,zk)
  do k = 1,nc

     zk = phi_bulk(m)%fparm(5) + (k-1)*phi_bulk(m)%dlc 

     do i = 1, na, stepa
        do j = 1, nb, stepb
           
           xi = phi_bulk(m)%fparm(1)  + (i-1)*phi_bulk(m)%dla
           yj = phi_bulk(m)%fparm(3)  + (j-1)*phi_bulk(m)%dlb
                         
           do atom = istart,iend            
              
              distR(1) = xi - x(a,atom)
              distR(2) = yj - x(b,atom)
              distR(3) = zk - x(c,atom)
              nsh = lmax(izp(atom))+1 

              ! Compute L-independent part:
              call long_pot(distR,basis,recbasis,alpha,vol,tol,lng_pot)
              ! total atomic charge
              deltaQ = sum(dQmat(1:nsh,atom)) 
              phi_bulk_PAR(i,j,k) = phi_bulk_PAR(i,j,k) + deltaQ*lng_pot
              
              ! compute L-dependent part:
              do l = 1, nsh
                 deltaQ = dQmat(l,atom) 
                 uhatm = uhubb(l,izp(atom)) 

                 call short_pot(distR,basis,uhatm,deltaQ,tol,sh_pot)
              
                 !OMP CRITICAL
                 phi_bulk_PAR(i,j,k) = phi_bulk_PAR(i,j,k) - sh_pot
                 !OMP END CRITICAL
              enddo 
           end do
        end do
     end do
  end do
  !OMP END PARALLEL DO
 
  phi_bulk(m)%val(:,:,:)=phi_bulk_PAR(:,:,:)
  
  call log_gdeallocate(phi_bulk_PAR)
  

end subroutine  compbulk_pot_ewald
!!$
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine save_bulkpot(phi_bulk,m)

  type(super_array) :: phi_bulk(:)
  integer :: m,i,j,k, fp, fp2
  character(2) :: m_id
  real(dp) :: xk

  write(m_id,'(i2.2)') m
  
  open(newunit=fp,file='contacts/BulkPot_'//m_id//'.dat',form='formatted')
 
  write(fp,'(3(i5))') phi_bulk(m)%iparm(14), &
                      phi_bulk(m)%iparm(15), &
                      phi_bulk(m)%iparm(16)
  
  do i = 1,phi_bulk(m)%iparm(14)  
     do j = 1,phi_bulk(m)%iparm(15)
        do k = 1,phi_bulk(m)%iparm(16)
           write(fp,*) phi_bulk(m)%val(i,j,k)
        end do
     end do
  end do
  
  close(fp)
 
  open(newunit=fp,file='contacts/Xvector_'//m_id//'.dat')
  do k = 1,phi_bulk(m)%iparm(14) 
     xk = phi_bulk(m)%fparm(1)+(k-1)*phi_bulk(m)%dla
     write(fp,'(E17.8)',ADVANCE='NO') xk * Bohr__AA
  enddo
  close(fp)
  open(newunit=fp,file='contacts/Yvector_'//m_id//'.dat')
  do k = 1,phi_bulk(m)%iparm(15)
     xk = phi_bulk(m)%fparm(3)+(k-1)*phi_bulk(m)%dlb
     write(fp,'(E17.8)',ADVANCE='NO') xk * Bohr__AA
  enddo
  close(fp)
  open(newunit=fp,file='contacts/Zvector_'//m_id//'.dat')
  do k = 1,phi_bulk(m)%iparm(16)    
     xk = phi_bulk(m)%fparm(5)+(k-1)*phi_bulk(m)%dlc
     write(fp,'(E17.8)',ADVANCE='NO') xk * Bohr__AA
  enddo
  close(fp)
  open(newunit=fp,file='contacts/box3d_'//m_id//'.dat') 
  write(fp,*) phi_bulk(m)%iparm(14),phi_bulk(m)%iparm(15),phi_bulk(m)%iparm(16) 
  close(fp)
 
end subroutine save_bulkpot


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine compbulk_pot_mud(phi_bulk,iparm,fparm, iErr)
  type(super_array) :: phi_bulk(:)
  integer :: iparm(23)
  real(kind=dp) :: fparm(8)
  integer, intent(out), optional :: iErr
 
  integer :: m,err,mgopt(4), a, b, c, i, cont
  real(dp), allocatable, dimension(:) :: work

  if (present(iErr)) then
    iErr = 0
  end if

  do m =1, ncont

    ! set parameters for mudpack

    a= phi_bulk(m)%a
    b= phi_bulk(m)%b
    c= phi_bulk(m)%c
    
    phi_bulk(m)%fparm(1) = fparm(2*a-1)
    phi_bulk(m)%fparm(2) = fparm(2*a)
    phi_bulk(m)%fparm(3) = fparm(2*b-1)
    phi_bulk(m)%fparm(4) = fparm(2*b)
    if (contdir(m).gt.0) then
       phi_bulk(m)%fparm(5) = fparm(2*c)
       phi_bulk(m)%fparm(6) = fparm(2*c)+phi_bulk(m)%L_PL
    else
       phi_bulk(m)%fparm(5) = fparm(2*c-1)-phi_bulk(m)%L_PL
       phi_bulk(m)%fparm(6) = fparm(2*c-1)
    endif
    
    phi_bulk(m)%fparm(7) = PoissAcc      ! Desired accuracy

    !if(id0.and.verbose.gt.80) then
    !   write(*,*)
    !   write(*,*) 'Bulk potential, contact',m
    !   if (phi_bulk(m)%doEwald) then
    !       write(*,*) 'BC = all periodic solved with Ewalds on two planes'  
    !   endif  
    !   write(*,*) 'X(a)=',phi_bulk(m)%fparm(1)*Bohr__AA,phi_bulk(m)%fparm(2)*Bohr__AA, &
    !   &   boundary2string(phi_bulk(m)%iparm(2)), boundary2string(phi_bulk(m)%iparm(3))
    !   write(*,*) 'X(b)=',phi_bulk(m)%fparm(3)*Bohr__AA,phi_bulk(m)%fparm(4)*Bohr__AA, &
    !   &   boundary2string(phi_bulk(m)%iparm(4)), boundary2string(phi_bulk(m)%iparm(5))
    !   write(*,*) 'X(c)=',phi_bulk(m)%fparm(5)*Bohr__AA,phi_bulk(m)%fparm(6)*Bohr__AA, &
    !   &   boundary2string(phi_bulk(m)%iparm(6)), boundary2string(phi_bulk(m)%iparm(7))
    !   write(*,*) 'L(c)=',phi_bulk(m)%L_PL*Bohr__AA
    !   write(*,*) 'na=',phi_bulk(m)%iparm(14)
    !   write(*,*) 'nb=',phi_bulk(m)%iparm(15)
    !   write(*,*) 'nc=',phi_bulk(m)%iparm(16)
    !endif

    ! call Ewald sums to set Dirichlet BC on two faces
    if (phi_bulk(m)%doEwald) then
      call compbulk_pot_ewald(phi_bulk,m)
    endif

    ! set charge density (rhs of poisson)
    cont= m
    call set_bulk_rhs(phi_bulk,cont)

    ! solve poisson for bulkpot
    
    call log_gallocate(work,phi_bulk(m)%iparm(20))
    
    mgopt(1) = 0

    do i = 0,1  
      phi_bulk(m)%iparm(1)=i

      call mud3sp( phi_bulk(m)%iparm, phi_bulk(m)%fparm, work, &
                &  bulk_cofx, bulk_cofy, bulk_cofz, & 
                &  bulk_bndyc,phi_bulk(m)%rhs,phi_bulk(m)%val,mgopt,err )
      
      if (err.ne.0 .and. err.ne.9) then
        @:FORMATTED_ERROR_HANDLING(iErr, -1, '(A,I0)', 'Poisson solver error n=', err)
      endif
      if(err.eq.9) then
         call log_gdeallocate(work)
         call log_gallocate(work,phi_bulk(m)%iparm(21))          
      endif
    enddo

    call log_gdeallocate(work)

    if (phi_bulk(m)%iparm(22).eq.phi_bulk(m)%iparm(18)) then
      @:FORMATTED_ERROR_HANDLING(iErr, -2, '(A,E12.4)', 'Bulk potential not converged! Error:',&
          & phi_bulk(m)%fparm(8))
    endif

    if (id0) call save_bulkpot(phi_bulk,cont)

  end do

end Subroutine compbulk_pot_mud
!%--------------------------------------------------------------------------
Subroutine bulk_bndyc(kbdy,xory,yorz,alfa,gbdy)

  integer :: kbdy
  real(kind=dp) :: xory,yorz,alfa,gbdy

  alfa = 0.d0
  gbdy = 0.d0

end subroutine bulk_bndyc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine bulk_cofx(x,cxx,cx,cex)
  real(kind=dp) :: x,cxx,cx,cex
  
  cxx = 1.d0
  cx = 0.d0          
  cex = 0.d0    

end subroutine

Subroutine bulk_cofy(y,cyy,cy,cey)
  real(kind=dp) :: y,cyy,cy,cey
  
  cyy = 1.d0
  cy = 0.d0          
  cey = 0.d0    

end subroutine

Subroutine bulk_cofz(z,czz,cz,cez)
  real(kind=dp) :: z,czz,cz,cez
  
  czz = 1.d0
  cz = 0.d0          
  cez = 0.d0    

end subroutine

Subroutine bulk_coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)

  real(kind=dp) :: x,y,z,cxx,cyy,czz,cx,cy,cz,ce
  
  cxx = 1.d0
  cyy = 1.d0
  czz = 1.d0
  cx = 0.d0          
  cy = 0.d0  
  cz = 0.d0
  ce = 0.d0    

end subroutine bulk_coef
!%-------------------------------------------------------------------
subroutine set_bulk_rhs(phi_bulk,cont)

  type(super_array) :: phi_bulk(:)  
  integer :: cont

  real(kind=dp) :: tmp,dl(3),xmin(3),xmax(3)
  real(kind=dp) :: xi(3), deltaR
  real(kind=dp) :: amin, bmin, cmin

  integer :: imin(3),imax(3), m, ii, jj, kk, na, nb, nc
  integer :: c(3), i,j,k, atom, dir,nsh,l
  integer :: rag(3)

  m=cont

  dir=sign(1,contdir(m))
  phi_bulk(m)%rhs = 0.d0  
  dl(1)=phi_bulk(m)%dla; dl(2)=phi_bulk(m)%dlb; dl(3)=phi_bulk(m)%dlc;

  c(1)=phi_bulk(m)%a
  c(2)=phi_bulk(m)%b
  c(3)=phi_bulk(m)%c
  
  if (period) then
     rag(1) = 1
     rag(2) = 1
  else
     rag(1) = 0
     rag(2) = 0
  endif
  rag(3) = 1
  
  na=phi_bulk(m)%iparm(14) - rag(1)
  nb=phi_bulk(m)%iparm(15) - rag(2)
  nc=phi_bulk(m)%iparm(16) - rag(3)

  amin = phi_bulk(m)%fparm(1)
  bmin = phi_bulk(m)%fparm(3)
  cmin = phi_bulk(m)%fparm(5)
 

  do atom = iatc(3,m),iatc(3,m)+phi_bulk(m)%natm_PL-1
    nsh = lmax(izp(atom))+1
    do l = 1, nsh
      tmp=3.2d0*uhubb(l,izp(atom))   
      ! Set boundaries of a box around the atom 
      ! Reference system is rotated to a,b,c
      do i = 1,3
         xmin(i) = x(c(i),atom) - deltaR_max
         xmax(i) = x(c(i),atom) + deltaR_max
      enddo
      
      ! Cut out a box around the atom
      do i=1,3
        imin(i) = nint( (xmin(i) - phi_bulk(m)%fparm(2*i-1))/dl(i) ) !+ 1 
        imax(i) = nint( (xmax(i) - phi_bulk(m)%fparm(2*i-1))/dl(i) ) !+ 1
        ! In non periodic directions cut at PoissonBox boundaries 
        if ( i.eq.3 .or. period_dir(c(i)) ) then
        else
          imin(i) = max( 0, imin(i) )
          imax(i) = min( phi_bulk(m)%iparm(13+i)-1, imax(i) )
        endif
      enddo
     
      ! compute the charge density on the grid points
      do i = imin(1),imax(1)
        
        xi(c(1)) = amin + i*dl(1)
        ! in periodic directions this folds to the Box
        ii=modulo(i, na) + 1  
        
        do j = imin(2),imax(2)
           
           xi(c(2)) = bmin + j*dl(2)
           ! in periodic directions this folds to the Box
           jj=modulo(j, nb) + 1       
           
           do k = imin(3),imax(3)
             
             xi(c(3)) = cmin + k*dl(3)
             
             ! in periodic directions this folds to the Box
             kk=modulo(k, nc) + 1
             
             ! if contdir<1 then reverse the point along c
             ! in order to have the interface contact/device in np=1
             kk=(1-dir)/2*(nc+rag(3)+1) + dir*kk 
             
             !Compute distance from atom 
             !Reference rotation is performed
             deltaR = sqrt(dot_product(xi-x(:,atom), xi-x(:,atom)))
             ! add charge density contrib
             
             phi_bulk(m)%rhs(ii,jj,kk) = phi_bulk(m)%rhs(ii,jj,kk) &
                  -0.5d0* dQmat(l,atom)* (tmp**3)*exp(-tmp*deltaR)
             !rhs contains -4*pi in normalization term: -4pi*tau^3/(8pi) = -0.5*tau^3
             
           end do

        end do
      end do
    
    end do
  end do
     
  if (period) then
    phi_bulk(m)%rhs(na+rag(1),:,:) = phi_bulk(m)%rhs(1,:,:) 
    phi_bulk(m)%rhs(:,nb+rag(2),:) = phi_bulk(m)%rhs(:,1,:)
  endif

  if (dir.gt.0) phi_bulk(m)%rhs(:,:,nc+rag(3)) = phi_bulk(m)%rhs(:,:,1)
  if (dir.lt.0) phi_bulk(m)%rhs(:,:,1) = phi_bulk(m)%rhs(:,:,nc+rag(3))

end subroutine set_bulk_rhs
!%--------------------------------------------------------------------------


subroutine rec_pot(r,uhatm,basis,recbasis,vol,tol,nit,potential, iErr)
  real(dp) ::  r(3)
  real(dp) ::  uhatm
  real(dp) ::  basis(3,3), recbasis(3,3),  vol, tol
  integer :: nit
  real(dp) ::  potential
  integer, intent(out), optional :: iErr

  real(dp) ::  G(3),help 
  real(dp) :: lastshell,butlast,err,uhatm2
  integer nrezi, nreal, nmax, nmin
  integer i,j,k

  if (present(iErr)) then
    iErr = 0
  end if

  nmax = 100
  nmin = 2

  !evaluate reciprocal space term ( sum over G <> 0) ...  
  !/* sum over G until tolerance is reached */
  nrezi = 1
  err=1.0d8
  lastshell = 0.d0
  butlast=0.d0  
  uhatm2=10.24*uhatm*uhatm
  potential = 0.d0

  do WHILE (   (nrezi .le. nmax) .and. &
              ((nrezi .le. nmin).or.(err.gt.tol))  )

    lastshell = 0.d0  
  
    do i=-nrezi,nrezi,2*nrezi
       do j=-nrezi,nrezi
          do k=-nrezi,nrezi
             
             G(1)=i*recbasis(1,1)+j*recbasis(2,1)+k*recbasis(3,1)
             G(2)=i*recbasis(1,2)+j*recbasis(2,2)+k*recbasis(3,2)
             G(3)=i*recbasis(1,3)+j*recbasis(2,3)+k*recbasis(3,3)

             help = G(1)*G(1)+G(2)*G(2)+G(3)*G(3)
             help = exp(-help/(4.d0*uhatm2))/help
             help = cos(G(1)*r(1)+G(2)*r(2)+G(3)*r(3)) * help               

             !help = G(1)*G(1)+G(2)*G(2)+G(3)*G(3)
             !help2= (uhatm2/help + 1)*(uhatm2/help + 1)*help
             !help = (uhatm2/help)*(uhatm2/help)/help2
             !help = cos(G(1)*r(1)+G(2)*r(2)+G(3)*r(3)) * help
             
             lastshell = lastshell + help       
             
          end do
       end do
    end do
    do j=-nrezi,nrezi,2*nrezi
       do i=-nrezi+1,nrezi-1
          do k=-nrezi,nrezi
             
             G(1)=i*recbasis(1,1)+j*recbasis(2,1)+k*recbasis(3,1)
             G(2)=i*recbasis(1,2)+j*recbasis(2,2)+k*recbasis(3,2)
             G(3)=i*recbasis(1,3)+j*recbasis(2,3)+k*recbasis(3,3)

             help = G(1)*G(1)+G(2)*G(2)+G(3)*G(3)
             help = exp(-help/(4.d0*uhatm2))/help
             help = cos(G(1)*r(1)+G(2)*r(2)+G(3)*r(3)) * help 
             
             !help = G(1)*G(1)+G(2)*G(2)+G(3)*G(3)
             !help2= (uhatm2/help+1)*(uhatm2/help+1)*help
             !help = (uhatm2/help)*(uhatm2/help)/help2
             !help = cos(G(1)*r(1)+G(2)*r(2)+G(3)*r(3)) * help
             
             lastshell = lastshell + help       
             
          end do
       end do
    end do
    do k=-nrezi,nrezi,2*nrezi
       do i=-nrezi+1,nrezi-1
          do j=-nrezi+1,nrezi-1
             
             G(1)=i*recbasis(1,1)+j*recbasis(2,1)+k*recbasis(3,1)
             G(2)=i*recbasis(1,2)+j*recbasis(2,2)+k*recbasis(3,2)
             G(3)=i*recbasis(1,3)+j*recbasis(2,3)+k*recbasis(3,3)

             help = G(1)*G(1)+G(2)*G(2)+G(3)*G(3)
             help = exp(-help/(4.d0*uhatm2))/help
             help = cos(G(1)*r(1)+G(2)*r(2)+G(3)*r(3)) * help 
             
             lastshell = lastshell + help       
             
          end do
       end do
    end do

    potential = potential + lastshell
    err=abs(lastshell/potential)
    !butlast=lastshell
    nrezi = nrezi + 1
  end do
  potential=(4.d0*Pi*potential)/vol
  nit=nrezi-1

  ! Halt if tolerance not reached
  if ( err .gt. tol ) then
    @:FORMATTED_ERROR_HANDLING(iErr, -1, '(A,E12.4,1X,A,1X,E12.4)',&
        & 'Tolerance in rec_pot not reached in reciprocal space:', err, 'vs', tol)
  end if

end subroutine

function boundary2string(typ) result(str)
   integer :: typ
   character(3) :: str
              
   select case(typ)
   case(0)
       str=' P '
   case(1)
       str=' D '
   case(2)
       str=' N '     
   end select

end function boundary2string


end Module bulkpot
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

