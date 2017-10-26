module structure

  use gprecision
  use gallocation

   
  implicit none
  private
  
  integer,  public, save :: natoms, nbeweg
  integer, ALLOCATABLE, public, save :: izp(:)       ! specie
  !integer, ALLOCATABLE, public, save :: ind(:)       ! iAtomStart(N) starting from 1
  integer,  public, save :: nlat(3)
  integer,  public, save :: ntypes
  
  integer,  public, save :: ss_natoms                !supercell atoms
  integer, ALLOCATABLE, public, save :: ss_izp(:)
  integer, public, save :: ss_f(3)
  
  integer, public, allocatable, save :: lmax(:)        !(MAXTYP)
  real(dp), public, allocatable, save :: uhubb(:,:)    !(NORB,MAXTYP)
  
 
  real(kind=dp), ALLOCATABLE, public, save :: x(:,:)
  real(kind=dp), ALLOCATABLE, public, save :: ss_x(:,:)  !supercell coords
  
  real(kind=dp),  public, save :: boxsiz(3,3), xinvbox(3,3), xnullvec(3)
  real(kind=dp),  public, save :: nel
  real(kind=dp), ALLOCATABLE, public, save :: dQmat(:,:) ! nshells,natoms
  real(kind=dp), ALLOCATABLE,  public, save :: xm(:)     ! mass of atomic specie (ntypes)
  real(kind=dp), ALLOCATABLE,  public, save :: xmrc(:)   ! reciprocal mass       (ntypes)
  
  logical,  public, save :: period
  logical,  public, save :: period_dir(3)
  
  character(3), public, save :: atnames(92)

  public :: find_ntypes, buildsupercell, inversebox
  public :: shortvertice, difback
  public :: coordback, gamma_summind


  contains
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

  ! ------------------------------------------------------------
  ! Calculate the inverse matrix (xinvbox) of the basis (boxsiz)
  ! ------------------------------------------------------------
  subroutine inversebox

    real(dp) :: xhlp

    xhlp=-boxsiz(1,3)*boxsiz(2,2)*boxsiz(3,1)+boxsiz(1,2)*boxsiz(2,3)*boxsiz(3,1)
    xhlp=xhlp+boxsiz(1,3)*boxsiz(2,1)*boxsiz(3,2)
    xhlp=xhlp-boxsiz(1,1)*boxsiz(2,3)*boxsiz(3,2)
    xhlp=xhlp-boxsiz(1,2)*boxsiz(2,1)*boxsiz(3,3)
    xhlp=xhlp+boxsiz(1,1)*boxsiz(2,2)*boxsiz(3,3)
    
    xinvbox(1,1)=(-boxsiz(2,3)*boxsiz(3,2)+boxsiz(2,2)*boxsiz(3,3))/xhlp
    xinvbox(2,1)=(boxsiz(2,3)*boxsiz(3,1)-boxsiz(2,1)*boxsiz(3,3))/xhlp
    xinvbox(3,1)=(-boxsiz(2,2)*boxsiz(3,1)+boxsiz(2,1)*boxsiz(3,2))/xhlp
    
    xinvbox(1,2)=(boxsiz(1,3)*boxsiz(3,2)-boxsiz(1,2)*boxsiz(3,3))/xhlp
    xinvbox(2,2)=(-boxsiz(1,3)*boxsiz(3,1)+boxsiz(1,1)*boxsiz(3,3))/xhlp
    xinvbox(3,2)=(boxsiz(1,2)*boxsiz(3,1)-boxsiz(1,1)*boxsiz(3,2))/xhlp
    
    xinvbox(1,3)=(-boxsiz(1,3)*boxsiz(2,2)+boxsiz(1,2)*boxsiz(2,3))/xhlp
    xinvbox(2,3)=(boxsiz(1,3)*boxsiz(2,1)-boxsiz(1,1)*boxsiz(2,3))/xhlp
    xinvbox(3,3)=(-boxsiz(1,2)*boxsiz(2,1)+boxsiz(1,1)*boxsiz(2,2))/xhlp

  end subroutine inversebox
  
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

  ! ------------------------------------------------------
  ! Put a difference vector back into supercell
  ! ------------------------------------------------------

  subroutine difback(dif)

     implicit REAL(dp) (A-H,O-Z)
     real(dp) dif(3)
     real(dp) xx1,xy1,xz1
     
     xx1=dif(1)*xinvbox(1,1)+dif(2)*xinvbox(2,1)+dif(3)*xinvbox(3,1)
     xy1=dif(1)*xinvbox(1,2)+dif(2)*xinvbox(2,2)+dif(3)*xinvbox(3,2)
     xz1=dif(1)*xinvbox(1,3)+dif(2)*xinvbox(2,3)+dif(3)*xinvbox(3,3)
     
     if(xx1>0.5) xx1=xx1-1.d0
     if(xx1<-0.5) xx1=xx1+1.d0
     if(xy1>0.5) xy1=xy1-1.d0
     if(xy1<-0.5) xy1=xy1+1.d0
     if(xz1>0.5) xz1=xz1-1.d0
     if(xz1<-0.5) xz1=xz1+1.d0    
     
     dif(1)=xx1*boxsiz(1,1)+xy1*boxsiz(2,1)+xz1*boxsiz(3,1)
     dif(2)=xx1*boxsiz(1,2)+xy1*boxsiz(2,2)+xz1*boxsiz(3,2)
     dif(3)=xx1*boxsiz(1,3)+xy1*boxsiz(2,3)+xz1*boxsiz(3,3)

   end subroutine difback

   !
   ! Put Atom at (x,y,z) back into supercell
   !
   subroutine coordback(xx,yy,zz)
     
     implicit REAL(dp) (A-H,O-Z)
     real(dp) xx,yy,zz
     real(dp) xx1,xy1,xz1

     xx1=(xx-xnullvec(1))*xinvbox(1,1)+(yy-xnullvec(2))*xinvbox(2,1)+(zz-xnullvec(3))*xinvbox(3,1)
     xy1=(xx-xnullvec(1))*xinvbox(1,2)+(yy-xnullvec(2))*xinvbox(2,2)+(zz-xnullvec(3))*xinvbox(3,2)
     xz1=(xx-xnullvec(1))*xinvbox(1,3)+(yy-xnullvec(2))*xinvbox(2,3)+(zz-xnullvec(3))*xinvbox(3,3)
     
     if(xx1>0.5) xx1=xx1-1.d0
     if(xx1<-0.5) xx1=xx1+1.d0
     if(xy1>0.5) xy1=xy1-1.d0
     if(xy1<-0.5) xy1=xy1+1.d0
     if(xz1>0.5) xz1=xz1-1.d0
     if(xz1<-0.5) xz1=xz1+1.d0
     
     xx=xx1*boxsiz(1,1)+xy1*boxsiz(2,1)+xz1*boxsiz(3,1)+xnullvec(1)
     yy=xx1*boxsiz(1,2)+xy1*boxsiz(2,2)+xz1*boxsiz(3,2)+xnullvec(2)
     zz=xx1*boxsiz(1,3)+xy1*boxsiz(2,3)+xz1*boxsiz(3,3)+xnullvec(3)
 
   end subroutine coordback

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
     !call log_gallocate(ss_fat,ss_natoms)
     !call log_gallocate(ss_ind,ss_natoms+1)
     
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
                 !ss_atf(algn)=i+j*maxf+k*maxf*maxf !repres. in maxf base
                 algn=algn+1 
                 
              enddo
           enddo
        enddo
     enddo
     
   end subroutine buildsupercell
   
   
 end module structure
 
