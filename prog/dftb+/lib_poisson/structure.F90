!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *  
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program * 
!  under the LGPL licence.                                                *
!**************************************************************************
module structure

  use dftbp_accuracy, only : dp
  use gallocation
  use mpi_poisson
   
  implicit none
  private
  
  integer,  public, save :: natoms
  integer, allocatable, public, save :: izp(:)       ! specie
  integer,  public, save :: nlat(3)
  integer,  public, save :: ntypes
  
  real(kind=dp), allocatable, public, save :: x(:,:)
 
  real(kind=dp), allocatable, public, save :: ss_x(:,:)  !supercell coords
  integer,  public, save :: ss_natoms                !supercell atoms
  integer, allocatable, public, save :: ss_izp(:)
  integer, public, save :: ss_f(3)
  
  integer, public, allocatable, save :: lmax(:)        !(MAXTYP)
  real(dp), public, allocatable, save :: uhubb(:,:)    !(NORB,MAXTYP)
  
  
  real(kind=dp),  public, save :: boxsiz(3,3), xnullvec(3)
  logical,  public, save :: period
  logical,  public, save :: period_dir(3)
  
  
  real(kind=dp), allocatable, public, save :: dQmat(:,:) ! nshells,natoms
  
  character(3), public, save :: atnames(92)
  
  !! Renormalization volumes: to ensure charge neutrality
  real(kind=dp), public, allocatable :: renorm(:,:)


  public :: find_ntypes, buildsupercell
  public :: shortvertice
  public :: gamma_summind
  public :: init_structure, init_charges, init_skdata

  contains

  ! -----------------------------------------------------------------------------
  !  FILL UP Structure 
  ! -----------------------------------------------------------------------------
  subroutine init_structure(st_nAtom, st_nSpecies, st_specie0, st_x0, &
              st_latVecs, st_isperiodic)

    integer, intent(in)   :: st_nAtom          ! number of Atoms in central cell 
    integer, intent(in)   :: st_nSpecies       ! number of Species
    integer, intent(in)   :: st_specie0(:)     ! type of each Atoms (nAtom)
    real(dp), intent(in)  :: st_x0(:,:)        ! coordinates in central cell
    real(dp), intent(in)  :: st_latVecs(3,3)   ! lattice vectors
    logical, intent(in)   :: st_isperiodic     ! periodic structure

    integer :: i,j
    real(dp) :: side(3)
    integer :: dir(3)
 
    if (active_id) then
 
      natoms=st_nAtom
      ntypes=st_nSpecies
       
      dir = (/ 1, 2, 3 /)
      boxsiz(1:3,1:3) = 0_dp
      ! the latVec structure is :
      !   
      !  [e_x1, e_x2, e_x3]
      !  [e_y1, e_y2, e_y3] 
      !  [e_z1, e_z2, e_z3]
  
      ! computes the three sides
      ! side(:) = sqrt(sum(st_latVecs,dim=1)**2)
  
      boxsiz(:,:) = transpose(st_latVecs(:,:))
  
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

      period=st_isperiodic

      call log_gallocate(x,3,natoms)
      x(1:3,1:natoms)=st_x0(1:3,1:natoms)
  
      call log_gallocate(izp,natoms)   
      izp(1:natoms)=st_specie0(1:natoms)
  
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
  subroutine init_skdata(nShell, angShell, hubbU, err)
    integer, intent(in) :: nShell(:)
    integer, intent(in) :: angShell(:,:)
    real(dp), intent(in) :: hubbU(:,:)
    integer, intent(out) :: err

    integer :: i,j,nshells

    err=0

    if (active_id) then

      ! set maximum angular momentum per specie
      call log_gallocate(lmax,ntypes)
 
      !checks that all atoms have shells in s,p,d sequence
      do i = 1, ntypes
         do j= 1, nShell(i)
            if (angShell(j,i).ne.j-1) then
               err=1
            end if   
         enddo 
         lmax(i) = angShell(nShell(i),i)
      enddo
      if (err.ne.0) then
         return
      endif
  
      ! set Hubbard parameters
      nshells = maxval(lmax)+1
  
      call log_gallocate(uhubb,nshells,ntypes)
  
      do i = 1,ntypes
        do j = 1,nshells
           uhubb(j,i) = hubbU(j,i)
        enddo   
      end do
 
    endif
  
  end subroutine init_skdata

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
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

  ! ------------------------------------------------------------------
  subroutine find_ss(M,ii,jj,kk)

    implicit none
    
    integer :: M,ii,jj,kk
    integer :: maxf

    maxf=maxval(ss_f)

    ii=mod(M,maxf)+1
    jj=mod(M/maxf,maxf)+1
    kk=mod(M/maxf/maxf,maxf)+1
    
  end subroutine find_ss
  ! ------------------------------------------------------------------
  subroutine find_ntypes(err)

    integer :: n,err

    err=0
    ntypes = 1
    do n = 1,natoms 
       if (izp(n) > ntypes) ntypes = izp(n)
    end do
    if (ntypes > 20) then
       err=1
    end if
    
  end subroutine find_ntypes


  ! ----------------------------------------------------------------------
  ! find the lenght of the shortest vertex in your supercell 
  ! ---------------------------------------------------------------------

  real(dp) function shortvertice()

    real(dp) yhlp, testbox(6)
  
    testbox(1)=boxsiz(1,1)**2+boxsiz(1,2)**2+boxsiz(1,3)**2
    testbox(2)=boxsiz(2,1)**2+boxsiz(2,2)**2+boxsiz(2,3)**2
    testbox(3)=boxsiz(3,1)**2+boxsiz(3,2)**2+boxsiz(3,3)**2
    
    testbox(4)=(boxsiz(1,1)-boxsiz(2,1))**2+(boxsiz(1,2)-boxsiz(2,2))**2+ & 
         (boxsiz(1,3)-boxsiz(2,3))**2
    testbox(5)=(boxsiz(1,1)-boxsiz(3,1))**2+(boxsiz(1,2)-boxsiz(3,2))**2+ &
         (boxsiz(1,3)-boxsiz(3,3))**2
    testbox(6)=(boxsiz(3,1)-boxsiz(2,1))**2+(boxsiz(3,2)-boxsiz(2,2))**2+ & 
         (boxsiz(3,3)-boxsiz(2,3))**2
    
    yhlp=MIN(testbox(1),testbox(2),testbox(3),testbox(4),testbox(5),testbox(6))
    shortvertice=sqrt(yhlp)

  end function shortvertice 


   !
   !get the three summation limits for the matrix construction in
   !the gamma point approximation
   !output: nlat(3)
   !
   subroutine gamma_summind(slkcutoff)
     
     implicit none
     real(dp) slkcutoff
     real(dp) u(3),v(3),w(3),helpv(3),l,lu,lv,lw 
     
     !get vectors and lengths
     u(1) = boxsiz(1,1); u(2) = boxsiz(1,2); u(3) = boxsiz(1,3)
     v(1) = boxsiz(2,1); v(2) = boxsiz(2,2); v(3) = boxsiz(2,3)
     w(1) = boxsiz(3,1); w(2) = boxsiz(3,2); w(3) = boxsiz(3,3)
     lu = sqrt(u(1)**2 + u(2)**2 + u(3)**2)
     lv = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
     lw = sqrt(w(1)**2 + w(2)**2 + w(3)**2)
   
     !see, whether length of u is shorter than projected side of v
     CALL CROSS(u,v,helpv)
     l = min(lu,sqrt(helpv(1)**2 + helpv(2)**2 + helpv(3)**2)/lv)
     
     !see, whether length of u is shorter than projected side of w
     CALL CROSS(u,w,helpv)
     l = min(lu,sqrt(helpv(1)**2 + helpv(2)**2 + helpv(3)**2)/lw)
     
     !set nlat(1) according to slkcutoff. At least 1.
     nlat(1) = max(int(2.d0*slkcutoff/l),1)
     
     !for v,w do the same as for u...
     CALL CROSS(v,u,helpv)
     l = min(lv,sqrt(helpv(1)**2 + helpv(2)**2 + helpv(3)**2)/lu)
     
     CALL CROSS(v,w,helpv)
     l = min(lv,sqrt(helpv(1)**2 + helpv(2)**2 + helpv(3)**2)/lw)
     
     nlat(2) = max(int(2.d0*slkcutoff/l),1)
     
     CALL CROSS(w,u,helpv)
     l = min(lw,sqrt(helpv(1)**2 + helpv(2)**2 + helpv(3)**2)/lu)
     
     CALL CROSS(w,v,helpv)
     l = min(lw,sqrt(helpv(1)**2 + helpv(2)**2 + helpv(3)**2)/lv)
     
     nlat(3) = max(int(2.d0*slkcutoff/l),1)
     
   end subroutine gamma_summind

   subroutine CROSS( A, B, C) 
     IMPLICIT NONE
     REAL(kind=dp) ::  A(3), B(3), C(3)
     
     C(1)=A(2)*B(3)-A(3)*B(2)
     C(2)=A(3)*B(1)-A(1)*B(3)
     C(3)=A(1)*B(2)-A(2)*B(1)
   END subroutine CROSS

   !-------------------------------------------------------------------      
   ! This section builds the Super Structure for periodic systems
   ! 
   subroutine buildsupercell() 
 
     implicit none
     
     integer :: ijk(9),algn,nu,nv,nw,i,j,k,n

     
     if( .not.any(period_dir) ) period=.false.
     
     !
     ! Check boxsize according to the cutoff in SLK-data:
     !    the shortest distance of two vertices has to be twice the shortest
     !    cutoff
     !
     ! --------------------------------------------------------------
     ss_f(:)=1
     
     if (period) then
        
        do i=1,3
           if (period_dir(i)) ss_f(i)=(2*nlat(i)+1)  
        enddo
        ss_natoms=ss_f(1)*ss_f(2)*ss_f(3)*natoms 
        
     else
        
        ss_natoms=natoms
        ss_f(:)=1
    
     endif

     if (.not.allocated(ss_x)) call log_gallocate(ss_x,3,ss_natoms)
     if (.not.allocated(ss_izp)) call log_gallocate(ss_izp,ss_natoms)
     
     ijk(1)=0;  ijk(2)=-1; ijk(3)=1;  ijk(4)=-2;  ijk(5)=2;  
     ijk(6)=-3; ijk(7)=3;  ijk(8)=-4; ijk(9)=4; 
     
     algn=1
     do i=1,ss_f(1)
        do j=1,ss_f(2)
           do k=1,ss_f(3)
              do n=1,natoms
                 
                 nu=ijk(i); nv=ijk(j); nw=ijk(k); 
                 ss_x(:,algn)=x(:,n)+nu*boxsiz(1,:)+nv*boxsiz(2,:)+nw*boxsiz(3,:)
                 ss_izp(algn)=izp(n)
                 algn=algn+1 
                 
              enddo
           enddo
        enddo
     enddo
     
   end subroutine buildsupercell
   
   
 end module structure
 
