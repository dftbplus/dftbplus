c
c    file mudcom.f
c
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c .                                                             .
c .                  copyright (c) 1999 by UCAR                 .
c .                                                             .
c .       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH       .
c .                                                             .
c .                      all rights reserved                    .
c .                                                             .
c .                                                             .
c .                      MUDPACK version 5.0                    .
c .                                                             .
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c... author and specialist
c
c         John C. Adams (National Center for Atmospheric Research)
c         email: johnad@ucar.edu, phone: 303-497-1213

c... For MUDPACK 5.0 information, visit the website:
c    (http://www.scd.ucar.edu/css/software/mudpack)
c
c... purpose
c
c    mudcom.f is a common subroutines file containing subroutines
c    called by some or all of the real two- and three-dimensional
c    mudpack solvers.  mudcom.f must be loaded with any real mudpack
c    solver.
c

c----+|----------------------------------------------------------------|
c    set phif,rhsf input in arrays which include
c    virtual boundaries for phi (for all 2-d real codes)
c----+|----------------------------------------------------------------|
      subroutine swk2(nfx,nfy,phif,rhsf,phi,rhs)
      implicit none
      integer nfx,nfy,i,j
      double precision phif(nfx,nfy),rhsf(nfx,nfy)
      double precision phi(0:nfx+1,0:nfy+1),rhs(nfx,nfy)
      do j=1,nfy
        do i=1,nfx
          phi(i,j) = phif(i,j)
          rhs(i,j) = rhsf(i,j)
        end do
      end do
c----+|----------------------------------------------------------------|
c    set virtual boundaries in phi to zero
c----+|----------------------------------------------------------------|
      do j=0,nfy+1
        phi(0,j) = 0.d0
        phi(nfx+1,j) = 0.d0
      end do
      do i=0,nfx+1
        phi(i,0) = 0.d0
        phi(i,nfy+1) = 0.d0
      end do
      return
      end subroutine swk2

c----+|----------------------------------------------------------------|
c    transfer fine grid to coarse grid
c----+|----------------------------------------------------------------|
      subroutine trsfc2(nx,ny,phi,rhs,ncx,ncy,phic,rhsc)
      implicit none
      integer nx,ny,ncx,ncy,i,j,ic,jc
      double precision phi(0:nx+1,0:ny+1),rhs(nx,ny)
      double precision phic(0:ncx+1,0:ncy+1),rhsc(ncx,ncy)
c
c----+|----------------------------------------------------------------|
c    set virtual boundaries in phic to zero
c----+|----------------------------------------------------------------|
      do jc=0,ncy+1
        phic(0,jc) = 0.d0
        phic(ncx+1,jc) = 0.d0
      end do
      do ic=0,ncx+1
        phic(ic,0) = 0.d0
        phic(ic,ncy+1) = 0.d0
      end do

      if (ncx.lt.nx .and. ncy.lt.ny) then
c----+|----------------------------------------------------------------|
c    coarsening in both x and y
c----+|----------------------------------------------------------------|
        do jc=1,ncy
          j = jc+jc-1
          do ic=1,ncx
            i = ic+ic-1
            phic(ic,jc) = phi(i,j)
            rhsc(ic,jc) = rhs(i,j)
          end do
        end do
      else if (ncx.lt.nx .and. ncy.eq.ny) then
c----+|----------------------------------------------------------------|
c    coarsening in x only
c----+|----------------------------------------------------------------|
        do jc=1,ncy
          j = jc
          do ic=1,ncx
            i = ic+ic-1
            phic(ic,jc) = phi(i,j)
            rhsc(ic,jc) = rhs(i,j)
          end do
        end do
      else
c----+|----------------------------------------------------------------|
c    coarsening in y only
c----+|----------------------------------------------------------------|
        do jc=1,ncy
          j = jc+jc-1
          do ic=1,ncx
            i = ic
            phic(ic,jc) = phi(i,j)
            rhsc(ic,jc) = rhs(i,j)
          end do
        end do
      end if
      return
      end subroutine trsfc2

c----+|----------------------------------------------------------------|
c    restrict fine grid residual in resf to coarse grid in rhsc
c    using full weighting for all real 2d codes
c----+|----------------------------------------------------------------|
      subroutine res2(nx,ny,resf,ncx,ncy,rhsc,nxa,nxb,nyc,nyd)
      implicit none
      integer nx,ny,ncx,ncy,nxa,nxb,nyc,nyd
      integer i,j,ic,jc,im1,ip1,jm1,jp1,ix,jy
      double precision resf(nx,ny),rhsc(ncx,ncy)

c----+|----------------------------------------------------------------|
c    set x,y coarsening integer subscript scales
c----+|----------------------------------------------------------------|
      ix = 1
      if (ncx.eq.nx) ix = 0
      jy = 1
      if (ncy.eq.ny) jy = 0
c----+|----------------------------------------------------------------|
c    restrict on interior
c----+|----------------------------------------------------------------|
      if (ncy.lt.ny .and. ncx.lt.nx) then
c----+|----------------------------------------------------------------|
c    coarsening in both directions
c----+|----------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(i,j,ic,jc), SHARED(resf,rhsc,ncx,ncy)
        do jc=2,ncy-1
          j = jc+jc-1
          do ic=2,ncx-1
            i = ic+ic-1
            rhsc(ic,jc) = (resf(i-1,j-1)+resf(i+1,j-1)+resf(i-1,j+1)+
     +                resf(i+1,j+1)+2.d0*(resf(i-1,j)+resf(i+1,j)+
     +                resf(i,j-1)+resf(i,j+1))+4.d0*resf(i,j))*.0625d0
          end do
        end do
!$OMP END PARALLEL DO
      else if (ncy.eq.ny) then
c----+|----------------------------------------------------------------|
c    no coarsening in y but coarsening in x
c----+|----------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(i,j,ic,jc), SHARED(resf,rhsc,ncx,ncy)
        do jc=2,ncy-1
          j = jc
          do ic=2,ncx-1
            i = ic+ic-1
            rhsc(ic,jc) = (resf(i-1,j-1)+resf(i+1,j-1)+resf(i-1,j+1)+
     +                 resf(i+1,j+1)+2.d0*(resf(i-1,j)+resf(i+1,j)+
     +                 resf(i,j-1)+resf(i,j+1))+4.d0*resf(i,j))*.0625d0
          end do
        end do
!$OMP END PARALLEL DO
      else
c----+|----------------------------------------------------------------|
c    no coarsening in x but coarsening in y
c----+|----------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(i,j,ic,jc), SHARED(resf,rhsc,ncx,ncy)
        do jc=2,ncy-1
          j = jc+jc-1
          do ic=2,ncx-1
            i = ic
            rhsc(ic,jc) = (resf(i-1,j-1)+resf(i+1,j-1)+resf(i-1,j+1)+
     +                 resf(i+1,j+1)+2.d0*(resf(i-1,j)+resf(i+1,j)+
     +                 resf(i,j-1)+resf(i,j+1))+4.d0*resf(i,j))*.0625d0
          end do
        end do
!$OMP END PARALLEL DO
      end if
c----+|----------------------------------------------------------------|
c    set residual on boundaries
c----+|----------------------------------------------------------------|
      do jc=1,ncy,ncy-1
c----+|----------------------------------------------------------------|
c       y=yc,yd boundaries
c----+|----------------------------------------------------------------|
        j = jc+jy*(jc-1)
        jm1 = max0(j-1,2)
        jp1 = min0(j+1,ny-1)
        if (j.eq.1 .and. nyc.eq.0) jm1 = ny-1
        if (j.eq.ny .and. nyc.eq.0) jp1 = 2
c----+|----------------------------------------------------------------|
c       y=yc,yd and x=xa,xb cornors
c----+|----------------------------------------------------------------|
        do ic=1,ncx,ncx-1
          i = ic+ix*(ic-1)
          im1 = max0(i-1,2)
          ip1 = min0(i+1,nx-1)
          if (i.eq.1 .and. nxa.eq.0) im1 = nx-1
          if (i.eq.nx .and. nxa.eq.0) ip1 = 2
          rhsc(ic,jc) = (resf(im1,jm1)+resf(ip1,jm1)+resf(im1,jp1)+
     +                   resf(ip1,jp1)+2.d0*(resf(im1,j)+resf(ip1,j)+
     +                resf(i,jm1)+resf(i,jp1))+4.d0*resf(i,j))*.0625d0
        end do
c----+|----------------------------------------------------------------|
c       set y=yc,yd interior edges
c----+|----------------------------------------------------------------|
        do ic=2,ncx-1
          i = ic+ix*(ic-1)
          rhsc(ic,jc) = (resf(i-1,jm1)+resf(i+1,jm1)+resf(i-1,jp1)+
     +                 resf(i+1,jp1)+2.d0*(resf(i-1,j)+resf(i+1,j)+
     +                 resf(i,jm1)+resf(i,jp1))+4.d0*resf(i,j))*.0625d0
        end do
      end do
c----+|----------------------------------------------------------------|
c     set x=xa,xb interior edges
c----+|----------------------------------------------------------------|
      do ic=1,ncx,ncx-1
        i = ic+ix*(ic-1)
        im1 = max0(i-1,2)
        ip1 = min0(i+1,nx-1)
        if (i.eq.1 .and. nxa.eq.0) im1 = nx-1
        if (i.eq.nx .and. nxa.eq.0) ip1 = 2
 
        do jc=2,ncy-1
          j = jc+jy*(jc-1)
          rhsc(ic,jc) = (resf(im1,j-1)+resf(ip1,j-1)+resf(im1,j+1)+
     +                resf(ip1,j+1)+2.d0*(resf(im1,j)+resf(ip1,j)+
     +                resf(i,j-1)+resf(i,j+1))+4.d0*resf(i,j))*.0625d0
        end do
      end do
c----+|----------------------------------------------------------------|
c    set coarse grid residual zero on specified boundaries
c----+|----------------------------------------------------------------|
      if (nxa.eq.1) then
        do jc=1,ncy
          rhsc(1,jc) = 0.d0
        end do
        end if
        if (nxb.eq.1) then
        do jc=1,ncy
          rhsc(ncx,jc) = 0.d0
        end do
      end if
      if (nyc.eq.1) then
        do ic=1,ncx
          rhsc(ic,1) = 0.d0
        end do
      end if
      if (nyd.eq.1) then
        do ic=1,ncx
          rhsc(ic,ncy) = 0.d0
        end do
      end if
      return
      end subroutine res2 

c----+|----------------------------------------------------------------|
c     prolon2 modified from rgrd2u 11/20/97
c----+|----------------------------------------------------------------|
      subroutine prolon2(ncx,ncy,p,nx,ny,q,nxa,nxb,nyc,nyd,intpol)
      implicit none
      integer ncx,ncy,nx,ny,intpol,nxa,nxb,nyc,nyd
      double precision p(0:ncx+1,0:ncy+1),q(0:nx+1,0:ny+1)
      integer i,j,jc,ist,ifn,jst,jfn,joddst,joddfn
      ist = 1
      ifn = nx
      jst = 1
      jfn = ny
      joddst = 1
      joddfn = ny
      if (nxa.eq.1) then
        ist = 2
      end if
      if (nxb.eq.1) then
        ifn = nx-1
      end if
      if (nyc.eq.1) then
        jst = 2
        joddst = 3
      end if
      if (nyd.eq.1) then
        jfn = ny-1
        joddfn = ny-2
      end if

      if (intpol.eq.1 .or. ncy.lt.4) then
c----+|----------------------------------------------------------------|
c       linearly interpolate in y
c----+|----------------------------------------------------------------|
        if (ncy .lt. ny) then
c----+|----------------------------------------------------------------|
c         ncy grid is an every other point subset of ny grid
c         set odd j lines interpolating in x and then set even
c         j lines by averaging odd j lines
c----+|----------------------------------------------------------------|
          do j=joddst,joddfn,2
            jc = j/2+1
            call prolon1(ncx,p(0,jc),nx,q(0,j),nxa,nxb,intpol)
          end do
          do j=2,jfn,2
            do i=ist,ifn
              q(i,j) = 0.5d0*(q(i,j-1)+q(i,j+1))
            end do
          end do
c----+|----------------------------------------------------------------|
c         set periodic virtual boundaries if necessary
c----+|----------------------------------------------------------------|
          if (nyc.eq.0) then
            do i=ist,ifn
              q(i,0) = q(i,ny-1)
              q(i,ny+1) = q(i,2)
            end do
          end if
          return
        else
c----+|----------------------------------------------------------------|
c         ncy grid equals ny grid so interpolate in x only
c----+|----------------------------------------------------------------|
          do j=jst,jfn
            jc = j
            call prolon1(ncx,p(0,jc),nx,q(0,j),nxa,nxb,intpol)
          end do
c----+|----------------------------------------------------------------|
c         set periodic virtual boundaries if necessary
c----+|----------------------------------------------------------------|
          if (nyc.eq.0) then
            do i=ist,ifn
              q(i,0) = q(i,ny-1)
              q(i,ny+1) = q(i,2)
            end do
          end if
          return
        end if
      else
c----+|----------------------------------------------------------------|
c       cubically interpolate in y
c----+|----------------------------------------------------------------|
        if (ncy .lt. ny) then
c----+|----------------------------------------------------------------|
c         set every other point of ny grid by interpolating in x
c----+|----------------------------------------------------------------|
          do j=joddst,joddfn,2
            jc = j/2+1
            call prolon1(ncx,p(0,jc),nx,q(0,j),nxa,nxb,intpol)
          end do
c----+|----------------------------------------------------------------|
c         set deep interior of ny grid using values just
c         generated and symmetric cubic interpolation in y
c----+|----------------------------------------------------------------|
          do j=4,ny-3,2
            do i=ist,ifn
           q(i,j)=(-q(i,j-3)+9.d0*(q(i,j-1)+q(i,j+1))-q(i,j+3))*.0625d0
            end do
          end do
c----+|----------------------------------------------------------------|
c         interpolate from q at j=2 and j=ny-1
c----+|----------------------------------------------------------------|
          if (nyc.ne.0) then
c----+|----------------------------------------------------------------|
c           asymmetric formula near nonperiodic y boundaries
c----+|----------------------------------------------------------------|
            do i=ist,ifn
              q(i,2)=(5.d0*q(i,1)+15.d0*q(i,3)-5.d0*q(i,5)+
     +               q(i,7))*.0625d0
              q(i,ny-1)=(5.d0*q(i,ny)+15.d0*q(i,ny-2)-5.d0*q(i,ny-4)+
     +                  q(i,ny-6))*.0625d0
            end do
          else
c----+|----------------------------------------------------------------|
c           periodicity in y alows symmetric formula near bndys
c----+|----------------------------------------------------------------|
            do i=ist,ifn
              q(i,2) = (-q(i,ny-2)+9.d0*(q(i,1)+q(i,3))-q(i,5))*.0625d0
              q(i,ny-1)=(-q(i,ny-4)+9.d0*(q(i,ny-2)+q(i,ny))
     +                     -q(i,3))*.0625d0
              q(i,ny+1) = q(i,2)
              q(i,0) = q(i,ny-1)
            end do
          end if
          return
        else
c----+|----------------------------------------------------------------|
c         ncy grid is equals ny grid so interpolate in x only
c----+|----------------------------------------------------------------|
          do j=jst,jfn
            jc = j
            call prolon1(ncx,p(0,jc),nx,q(0,j),nxa,nxb,intpol)
          end do
c----+|----------------------------------------------------------------|
c         set periodic virtual boundaries if necessary
c----+|----------------------------------------------------------------|
          if (nyc.eq.0) then
            do i=ist,ifn
              q(i,0) = q(i,ny-1)
              q(i,ny+1) = q(i,2)
            end do
          end if
          return
        end if
      end if
      end subroutine prolon2

c----+|----------------------------------------------------------------|
c    11/20/97  modification of rgrd1u.f for mudpack
c----+|----------------------------------------------------------------|
      subroutine prolon1(ncx,p,nx,q,nxa,nxb,intpol)
      implicit none
      integer intpol,nxa,nxb,ncx,nx,i,ic,ist,ifn,ioddst,ioddfn
      double precision p(0:ncx+1),q(0:nx+1)
      ist = 1
      ioddst = 1
      ifn = nx
      ioddfn = nx
      if (nxa.eq.1) then
        ist = 2
        ioddst = 3
      end if
      if (nxb.eq.1) then
        ifn = nx-1
        ioddfn = nx-2
      end if

      if (intpol.eq.1 .or. ncx.lt.4) then
c----+|----------------------------------------------------------------|
c       linear interpolation in x
c----+|----------------------------------------------------------------|
        if (ncx .lt. nx) then
c----+|----------------------------------------------------------------|
c         every other point of nx grid is ncx grid
c----+|----------------------------------------------------------------|
          do i=ioddst,ioddfn,2
            ic = (i+1)/2
            q(i) = p(ic)
          end do
          do i=2,ifn,2
            q(i) = 0.5d0*(q(i-1)+q(i+1))
          end do
        else
c----+|----------------------------------------------------------------|
c         nx grid equals ncx grid
c----+|----------------------------------------------------------------|
          do i=ist,ifn
            q(i) = p(i)
          end do
        end if
c----+|----------------------------------------------------------------|
c       set virtual end points if periodic
c----+|----------------------------------------------------------------|
        if (nxa.eq.0) then
          q(0) = q(nx-1)
          q(nx+1) = q(2)
        end if
        return
      else
c----+|----------------------------------------------------------------|
c       cubic interpolation in x
c----+|----------------------------------------------------------------|
        if (ncx .lt. nx) then
          do i=ioddst,ioddfn,2
            ic = (i+1)/2
            q(i) = p(ic)
          end do
c----+|----------------------------------------------------------------|
c         set deep interior with symmetric formula
c----+|----------------------------------------------------------------|
          do i=4,nx-3,2
            q(i)=(-q(i-3)+9.d0*(q(i-1)+q(i+1))-q(i+3))*.0625d0
          end do
c----+|----------------------------------------------------------------|
c         interpolate from q at i=2 and i=nx-1
c----+|----------------------------------------------------------------|
          if (nxa.ne.0) then
c----+|----------------------------------------------------------------|
c           asymmetric formula near nonperiodic bndys
c----+|----------------------------------------------------------------|
            q(2)=(5.d0*q(1)+15.d0*q(3)-5.d0*q(5)+q(7))*.0625d0
            q(nx-1)=(5.d0*q(nx)+15.d0*q(nx-2)-5.d0*q(nx-4)+
     +              q(nx-6))*.0625d0
          else
c----+|----------------------------------------------------------------|
c           periodicity in x alows symmetric formula near bndys
c----+|----------------------------------------------------------------|
            q(2) = (-q(nx-2)+9.d0*(q(1)+q(3))-q(5))*.0625d0
            q(nx-1) = (-q(nx-4)+9.d0*(q(nx-2)+q(nx))-q(3))*.0625d0
            q(nx+1) = q(2)
            q(0) = q(nx-1)
          end if
          return
        else
c----+|----------------------------------------------------------------|
c         ncx grid equals nx grid
c----+|----------------------------------------------------------------|
          do i=ist,ifn
            q(i) = p(i)
          end do
          if (nxa.eq.0) then
            q(0) = q(nx-1)
            q(nx+1) = q(2)
          end if
          return
        end if
      end if
      end subroutine prolon1


c----+|----------------------------------------------------------------|
c     add coarse grid correction in phic to fine grid approximation
c     in phif using linear or cubic interpolation
c----+|----------------------------------------------------------------|
      subroutine cor2(nx,ny,phif,ncx,ncy,phic,nxa,nxb,nyc,nyd,intpol,
     +                phcor)
      implicit none
      integer i,j,nx,ny,ncx,ncy,nxa,nxb,nyc,nyd,intpol,ist,ifn,jst,jfn
      double precision phif(0:nx+1,0:ny+1),phic(0:ncx+1,0:ncy+1)
      double precision phcor(0:nx+1,0:ny+1)
      do j=0,ny+1
      do i=0,nx+1
        phcor(i,j) = 0.d0
      end do
      end do
c----+|----------------------------------------------------------------|
c     lift correction in phic to fine grid in phcor
c----+|----------------------------------------------------------------|
      call prolon2(ncx,ncy,phic,nx,ny,phcor,nxa,nxb,nyc,nyd,intpol)
c----+|----------------------------------------------------------------|
c     add correction in phcor to phif on nonspecified boundaries
c----+|----------------------------------------------------------------|
      ist = 1
      ifn = nx
      jst = 1
      jfn = ny
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.1) jst = 2
      if (nyd.eq.1) jfn = ny-1
      do j=jst,jfn
        do i=ist,ifn
          phif(i,j) = phif(i,j) + phcor(i,j)
        end do
      end do
c----+|----------------------------------------------------------------|
c    add periodic points if necessary
c----+|----------------------------------------------------------------|
      if (nyc.eq.0) then
        do i=ist,ifn
          phif(i,0) = phif(i,ny-1)
          phif(i,ny+1) = phif(i,2)
        end do
      end if
      if (nxa.eq.0) then
        do j=jst,jfn
          phif(0,j) = phif(nx-1,j)
          phif(nx+1,j) = phif(2,j)
        end do
      end if
      end subroutine cor2

c----+|----------------------------------------------------------------|
c     use second order approximation in u to estimate (second order)
c     third and fourth partial derivatives in the x and y direction
c     non-symmetric difference formula (derived from the  routine
c     finpdf,findif) are used at and one point in from mixed boundaries.
c----+|----------------------------------------------------------------|
      subroutine pde2(nx,ny,u,i,j,ux3,ux4,uy3,uy4,nxa,nyc)
      implicit none
      integer nx,ny,i,j,nxa,nyc
      double precision u(nx,ny),dlx,dly,dlxx,dlyy,tdlx3,tdly3,dlx4,dly4
      common/pde2com/dlx,dly,dlxx,dlyy,tdlx3,tdly3,dlx4,dly4
      double precision ux3,ux4,uy3,uy4

c----+|----------------------------------------------------------------|
c     x partial derivatives
c----+|----------------------------------------------------------------|
      if (nxa.ne.0) then
c----+|----------------------------------------------------------------|
c       non periodic in x
c----+|----------------------------------------------------------------|
        if(i.gt.2 .and. i.lt.nx-1) then
          ux3 = (-u(i-2,j)+2.d0*u(i-1,j)-2.d0*u(i+1,j)+u(i+2,j))/tdlx3
          ux4 = (u(i-2,j)-4.d0*u(i-1,j)+6.d0*u(i,j)-4.d0*u(i+1,j)+
     +          u(i+2,j))/dlx4
        else if (i.eq.1) then
          ux3 = (-5.d0*u(1,j)+18.d0*u(2,j)-24.d0*u(3,j)+14.d0*u(4,j)-
     +           3.d0*u(5,j))/tdlx3
          ux4 = (3.d0*u(1,j)-14.d0*u(2,j)+26.d0*u(3,j)-24.d0*u(4,j)+
     +           11.d0*u(5,j)-2.d0*u(6,j))/dlx4
        else if (i.eq.2) then
          ux3 = (-3.d0*u(1,j)+10.d0*u(2,j)-12.d0*u(3,j)+6.d0*u(4,j)-
     +          u(5,j))/tdlx3
          ux4 = (2.d0*u(1,j)-9.d0*u(2,j)+16.d0*u(3,j)-14.d0*u(4,j)+
     +           6.d0*u(5,j)-u(6,j))/dlx4
        else if (i.eq.nx-1) then
          ux3 = (u(nx-4,j)-6.d0*u(nx-3,j)+12.d0*u(nx-2,j)-
     +          10.d0*u(nx-1,j)+3.d0*u(nx,j))/tdlx3
          ux4 = (-u(nx-5,j)+6.d0*u(nx-4,j)-14.d0*u(nx-3,j)+
     +         16.d0*u(nx-2,j)-9.d0*u(nx-1,j)+2.d0*u(nx,j))/dlx4
        else if (i.eq.nx) then
          ux3 = (3.d0*u(nx-4,j)-14.d0*u(nx-3,j)+24.d0*u(nx-2,j)-
     +           18.d0*u(nx-1,j)+5.d0*u(nx,j))/tdlx3
          ux4 = (-2.d0*u(nx-5,j)+11.d0*u(nx-4,j)-24.d0*u(nx-3,j)+
     +           26.d0*u(nx-2,j)-14.d0*u(nx-1,j)+3.d0*u(nx,j))/dlx4
        end if

      else
c----+|----------------------------------------------------------------|
c       periodic in x
c----+|----------------------------------------------------------------|
        if(i.gt.2 .and. i.lt.nx-1) then
          ux3 = (-u(i-2,j)+2.d0*u(i-1,j)-2.d0*u(i+1,j)+u(i+2,j))/tdlx3
          ux4 = (u(i-2,j)-4.d0*u(i-1,j)+6.d0*u(i,j)-4.d0*u(i+1,j)+
     +           u(i+2,j))/dlx4
        else if (i.eq.1) then
          ux3 = (-u(nx-2,j)+2.d0*u(nx-1,j)-2.d0*u(2,j)+u(3,j))/tdlx3
          ux4 = (u(nx-2,j)-4.d0*u(nx-1,j)+6.d0*u(1,j)-4.d0*u(2,j)+
     +          u(3,j))/dlx4
        else if (i.eq.2) then
          ux3 = (-u(nx-1,j)+2.d0*u(1,j)-2.d0*u(3,j)+u(4,j))/(tdlx3)
          ux4 = (u(nx-1,j)-4.d0*u(1,j)+6.d0*u(2,j)-4.d0*u(3,j)+
     +          u(4,j))/dlx4
        else if (i.eq.nx-1) then
          ux3 = (-u(nx-3,j)+2.d0*u(nx-2,j)-2.d0*u(1,j)+u(2,j))/tdlx3
          ux4 = (u(nx-3,j)-4.d0*u(nx-2,j)+6.d0*u(nx-1,j)-4.d0*u(1,j)+
     +           u(2,j))/dlx4
        else if (i.eq.nx) then
          ux3 = (-u(nx-2,j)+2.d0*u(nx-1,j)-2.d0*u(2,j)+u(3,j))/tdlx3
          ux4 = (u(nx-2,j)-4.d0*u(nx-1,j)+6.d0*u(nx,j)-4.d0*u(2,j)+
     +          u(3,j))/dlx4
        end if
      end if

c----+|----------------------------------------------------------------|
c     y partial derivatives
c----+|----------------------------------------------------------------|
      if (nyc.ne.0) then
c----+|----------------------------------------------------------------|
c       non periodic in y
c----+|----------------------------------------------------------------|
        if (j.gt.2 .and. j.lt.ny-1) then
          uy3 = (-u(i,j-2)+2.d0*u(i,j-1)-2.d0*u(i,j+1)+u(i,j+2))/tdly3
          uy4 = (u(i,j-2)-4.d0*u(i,j-1)+6.d0*u(i,j)-4.d0*u(i,j+1)+
     +            u(i,j+2))/dly4
        else if (j.eq.1) then
          uy3 = (-5.d0*u(i,1)+18.d0*u(i,2)-24.d0*u(i,3)+14.d0*u(i,4)-
     +              3.d0*u(i,5))/tdly3
          uy4 = (3.d0*u(i,1)-14.d0*u(i,2)+26.d0*u(i,3)-24.d0*u(i,4)+
     +             11.d0*u(i,5)-2.d0*u(i,6))/dly4
        else if (j.eq.2) then
          uy3 = (-3.d0*u(i,1)+10.d0*u(i,2)-12.d0*u(i,3)+6.d0*u(i,4)-
     +            u(i,5))/tdly3
          uy4 = (2.d0*u(i,1)-9.d0*u(i,2)+16.d0*u(i,3)-14.d0*u(i,4)+
     +             6.d0*u(i,5)-u(i,6))/dly4
        else if (j.eq.ny-1) then
          uy3 = (u(i,ny-4)-6.d0*u(i,ny-3)+12.d0*u(i,ny-2)-
     +            10.d0*u(i,ny-1)+3.d0*u(i,ny))/tdly3
          uy4 = (-u(i,ny-5)+6.d0*u(i,ny-4)-14.d0*u(i,ny-3)+
     +            16.d0*u(i,ny-2)-9.d0*u(i,ny-1)+2.d0*u(i,ny))/dly4
        else if (j.eq.ny) then
          uy3 = (3.d0*u(i,ny-4)-14.d0*u(i,ny-3)+24.d0*u(i,ny-2)-
     +             18.d0*u(i,ny-1)+5.d0*u(i,ny))/tdly3
          uy4 = (-2.d0*u(i,ny-5)+11.d0*u(i,ny-4)-24.d0*u(i,ny-3)+
     +             26.d0*u(i,ny-2)-14.d0*u(i,ny-1)+3.d0*u(i,ny))/dly4
        end if

      else
c----+|----------------------------------------------------------------|
c       periodic in y
c----+|----------------------------------------------------------------|
        if (j.gt.2 .and. j.lt.ny-1) then
          uy3 = (-u(i,j-2)+2.d0*u(i,j-1)-2.d0*u(i,j+1)+u(i,j+2))/tdly3
          uy4 = (u(i,j-2)-4.d0*u(i,j-1)+6.d0*u(i,j)-4.d0*u(i,j+1)+
     +            u(i,j+2))/dly4
        else if (j.eq.1) then
          uy3 = (-u(i,ny-2)+2.d0*u(i,ny-1)-2.d0*u(i,2)+u(i,3))/tdly3
          uy4 = (u(i,ny-2)-4.d0*u(i,ny-1)+6.d0*u(i,1)-4.d0*u(i,2)+
     +            u(i,3))/dly4
        else if (j.eq.2) then
          uy3 = (-u(i,ny-1)+2.d0*u(i,1)-2.d0*u(i,3)+u(i,4))/(tdly3)
          uy4 = (u(i,ny-1)-4.d0*u(i,1)+6.d0*u(i,2)-4.d0*u(i,3)+
     +            u(i,4))/dly4
        else if (j.eq.ny-1) then
          uy3 = (-u(i,ny-3)+2.d0*u(i,ny-2)-2.d0*u(i,1)+u(i,2))/tdly3
          uy4 = (u(i,ny-3)-4.d0*u(i,ny-2)+6.d0*u(i,ny-1)-4.d0*u(i,1)+
     +             u(i,2))/dly4
        else if (j.eq.ny) then
          uy3 = (-u(i,ny-2)+2.d0*u(i,ny-1)-2.d0*u(i,2)+u(i,3))/tdly3
          uy4 = (u(i,ny-2)-4.d0*u(i,ny-1)+6.d0*u(i,ny)-4.d0*u(i,2)+
     +            u(i,3))/dly4
        end if
      end if
      return
      end subroutine pde2

c----+|----------------------------------------------------------------|
c     set phif,rhsf input in arrays which include
c     virtual boundaries for phi (for all 2-d real codes)
c----+|----------------------------------------------------------------|
      subroutine swk3(nfx,nfy,nfz,phif,rhsf,phi,rhs)
      implicit none
      integer nfx,nfy,nfz,i,j,k
      double precision phif(nfx,nfy,nfz),rhsf(nfx,nfy,nfz)
      double precision phi(0:nfx+1,0:nfy+1,0:nfz+1),rhs(nfx,nfy,nfz)

      do k=1,nfz
        do j=1,nfy
          do i=1,nfx
            phi(i,j,k) = phif(i,j,k)
            rhs(i,j,k) = rhsf(i,j,k)
          end do
        end do
      end do

c----+|----------------------------------------------------------------|
c    set virtual boundaries in phi to zero
c----+|----------------------------------------------------------------|
      do k=0,nfz+1
        do j=0,nfy+1
          phi(0,j,k) = 0.d0
          phi(nfx+1,j,k) = 0.d0
        end do
      end do
      do k=0,nfz+1
        do i=0,nfx+1
          phi(i,0,k) = 0.d0
          phi(i,nfy+1,k) = 0.d0
        end do
      end do
      do j=0,nfy+1
        do i=0,nfx+1
          phi(i,j,0) = 0.d0
          phi(i,j,nfz+1) = 0.d0
        end do
      end do
      return
      end subroutine swk3

c----+|----------------------------------------------------------------|
c     transfers fine grid to coarse grid
c----+|----------------------------------------------------------------|
      subroutine trsfc3(nx,ny,nz,phi,rhs,ncx,ncy,ncz,phic,rhsc)
      implicit none
      integer nx,ny,nz,ncx,ncy,ncz,i,j,k,ic,jc,kc,ix,jy,kz
      double precision phi(0:nx+1,0:ny+1,0:nz+1),rhs(nx,ny,nz)
      double precision phic(0:ncx+1,0:ncy+1,0:ncz+1),rhsc(ncx,ncy,ncz)

c----+|----------------------------------------------------------------|
c    set virtual boundaries in phic to zero
c----+|----------------------------------------------------------------|
      do kc=0,ncz+1
        do jc=0,ncy+1
          phic(0,jc,kc) = 0.d0
          phic(ncx+1,jc,kc) = 0.d0
        end do
      end do
      do kc=0,ncz+1
        do ic=0,ncx+1
          phic(ic,0,kc) = 0.0
          phic(ic,ncy+1,kc) = 0.0
        end do
      end do
      do jc=0,ncy+1
        do ic=0,ncx+1
          phic(ic,jc,0) = 0.0
          phic(ic,jc,ncz+1) = 0.0
        end do
      end do

c----+|----------------------------------------------------------------|
c     coarsening in x,y,z 
c     (no coarsening for equal grids)
c----+|----------------------------------------------------------------|
      ix = 1
      if (ncx.eq.nx) ix = 0
      jy = 1
      if (ncy.eq.ny) jy = 0
      kz = 1
      if (ncz.eq.nz) kz = 0

!$OMP PARALLEL DO PRIVATE(k,kc,j,jc,i,ic) 
      do kc=1,ncz
        k = kc+kz*(kc-1)
        do jc=1,ncy
          j = jc+jy*(jc-1)
          do ic=1,ncx
            i = ic+ix*(ic-1)
            phic(ic,jc,kc) = phi(i,j,k)
            rhsc(ic,jc,kc) = rhs(i,j,k)
          end do
        end do
      end do
!$OMP END PARALLEL DO
      
      return
      end subroutine trsfc3

c----+|----------------------------------------------------------------|
c     restrict fine grid residual in resf to coarse grid in rhsc
c     using full weighting
c----+|----------------------------------------------------------------|
      subroutine res3(nx,ny,nz,resf,ncx,ncy,ncz,rhsc,
     +                nxa,nxb,nyc,nyd,nze,nzf)
      implicit none
      integer nx,ny,nz,ncx,ncy,ncz,nxa,nxb,nyc,nyd,nze,nzf
      integer ix,jy,kz,i,j,k,ic,jc,kc,im1,ip1,jm1,jp1,km1,kp1
      double precision rm,rk,rp
      double precision resf(nx,ny,nz),rhsc(ncx,ncy,ncz)

c----+|----------------------------------------------------------------|
c    set x,y,z coarsening integer subscript scales
c----+|----------------------------------------------------------------|
      ix = 1
      if (ncx.eq.nx) ix = 0
      jy = 1
      if (ncy.eq.ny) jy = 0
      kz = 1
      if (ncz.eq.nz) kz = 0
c----+|----------------------------------------------------------------|
c     restrict on interior
c     coarsening in x,y,z
c----+|----------------------------------------------------------------|
      if (ncz.lt.nz .and. ncy.lt.ny .and. ncx.lt.nx) then
!$OMP PARALLEL DO PRIVATE(i,j,k,ic,jc,kc,rm,rk,rp)
!$OMP& SHARED(resf,rhsc,ncx,ncy,ncz)
        do kc=2,ncz-1
          k = kc+kc-1
          do jc=2,ncy-1
            j = jc+jc-1
            do ic=2,ncx-1
              i = ic+ic-1
c----+|----------------------------------------------------------------|
c             weight on k-1,k,k+1 z planes in rm,rk,rp
c----+|----------------------------------------------------------------|
              rm=(resf(i-1,j-1,k-1)+resf(i+1,j-1,k-1)+
     +           resf(i-1,j+1,k-1)+resf(i+1,j+1,k-1)+
     +           2.d0*(resf(i-1,j,k-1)+resf(i+1,j,k-1)+
     +           resf(i,j-1,k-1)+resf(i,j+1,k-1))+
     +           4.d0*resf(i,j,k-1))*.0625d0
              rk=(resf(i-1,j-1,k)+resf(i+1,j-1,k)+resf(i-1,j+1,k)+
     +           resf(i+1,j+1,k)+2.d0*(resf(i-1,j,k)+resf(i+1,j,k)+
     +           resf(i,j-1,k)+resf(i,j+1,k))+4.d0*resf(i,j,k))*.0625d0
              rp=(resf(i-1,j-1,k+1)+resf(i+1,j-1,k+1)+
     +           resf(i-1,j+1,k+1)+resf(i+1,j+1,k+1)+
     +           2.d0*(resf(i-1,j,k+1)+resf(i+1,j,k+1)+resf(i,j-1,k+1)+
     +           resf(i,j+1,k+1))+4.d0*resf(i,j,k+1))*.0625d0
c----+|----------------------------------------------------------------|
c             weight in z direction for final result
c----+|----------------------------------------------------------------|
              rhsc(ic,jc,kc) = 0.25d0*(rm+2.d0*rk+rp)
            end do
          end do
        end do
!$OMP END PARALLEL DO
      else
c----+|----------------------------------------------------------------|
c       allow for noncoarsening in any of x,y,z
c----+|----------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(i,j,k,ic,jc,kc,rm,rk,rp)
!$OMP& SHARED(ix,jy,kz,resf,rhsc,ncx,ncy,ncz)
        do kc=2,ncz-1
          k = kc+kz*(kc-1)
          do jc=2,ncy-1
            j = jc+jy*(jc-1)
            do ic=2,ncx-1
              i = ic+ix*(ic-1)
c----+|----------------------------------------------------------------|
c             weight on k-1,k,k+1 z planes in rm,rk,rp
c----+|----------------------------------------------------------------|
              rm=(resf(i-1,j-1,k-1)+resf(i+1,j-1,k-1)+
     +           resf(i-1,j+1,k-1)+resf(i+1,j+1,k-1)+
     +           2.d0*(resf(i-1,j,k-1)+resf(i+1,j,k-1)+resf(i,j-1,k-1)+
     +           resf(i,j+1,k-1))+4.d0*resf(i,j,k-1))*.0625d0
              rk=(resf(i-1,j-1,k)+resf(i+1,j-1,k)+resf(i-1,j+1,k)+
     +           resf(i+1,j+1,k)+2.d0*(resf(i-1,j,k)+resf(i+1,j,k)+
     +           resf(i,j-1,k)+resf(i,j+1,k))+4.d0*resf(i,j,k))*.0625d0
              rp=(resf(i-1,j-1,k+1)+resf(i+1,j-1,k+1)+resf(i-1,j+1,k+1)+
     +           resf(i+1,j+1,k+1)+2.d0*(resf(i-1,j,k+1)+
     +           resf(i+1,j,k+1)+resf(i,j-1,k+1)+resf(i,j+1,k+1))+
     +           4.d0*resf(i,j,k+1))*.0625d0
c----+|----------------------------------------------------------------|
c             weight in z direction for final result
c----+|----------------------------------------------------------------|
              rhsc(ic,jc,kc) = 0.25*(rm+2.d0*rk+rp)
            end do
          end do
        end do
!$OMP END PARALLEL DO
      end if

c----+|----------------------------------------------------------------|
c     set residual on boundaries
c----+|----------------------------------------------------------------|
      do ic=1,ncx,ncx-1
c----+|----------------------------------------------------------------|
c       x=xa and x=xb
c----+|----------------------------------------------------------------|
        i = ic+ix*(ic-1)
        im1 = max0(i-1,2)
        ip1 = min0(i+1,nx-1)
        if (i.eq.1 .and. nxa.eq.0) im1 = nx-1
        if (i.eq.nx .and. nxb.eq.0) ip1 = 2
c----+|----------------------------------------------------------------|
c       (y,z) interior
c----+|----------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(j,k,jc,kc,rm,rk,rp)
!$OMP& SHARED(kz,jy,ic,im1,i,ip1,resf,rhsc,ncy,ncz)
        do kc=2,ncz-1
          k = kc+kz*(kc-1)
          do jc=2,ncy-1
            j = jc+jy*(jc-1)
            rm=(resf(im1,j-1,k-1)+resf(ip1,j-1,k-1)+resf(im1,j+1,k-1)+
     +         resf(ip1,j+1,k-1)+2.d0*(resf(im1,j,k-1)+resf(ip1,j,k-1)+
     +         resf(i,j-1,k-1)+resf(i,j+1,k-1))+
     +         4.d0*resf(i,j,k-1))*.0625d0
            rk=(resf(im1,j-1,k)+resf(ip1,j-1,k)+resf(im1,j+1,k)+
     +         resf(ip1,j+1,k)+2.d0*(resf(im1,j,k)+resf(ip1,j,k)+
     +         resf(i,j-1,k)+resf(i,j+1,k))+4.d0*resf(i,j,k))*.0625d0
            rp=(resf(im1,j-1,k+1)+resf(ip1,j-1,k+1)+resf(im1,j+1,k+1)+
     +         resf(ip1,j+1,k+1)+2.d0*(resf(im1,j,k+1)+resf(ip1,j,k+1)+
     +         resf(i,j-1,k+1)+resf(i,j+1,k+1))+
     +         4.d0*resf(i,j,k+1))*.0625d0
            rhsc(ic,jc,kc) = 0.25d0*(rm+2.d0*rk+rp)
          end do
        end do
!$OMP END PARALLEL DO
c----+|----------------------------------------------------------------|
c       x=xa,xb and y=yc,yd interior edges
c----+|----------------------------------------------------------------|
        do jc=1,ncy,ncy-1
          j = jc+jy*(jc-1)
          jm1 = max0(j-1,2)
          jp1 = min0(j+1,ny-1)
          if (j.eq.1 .and. nyc.eq.0) jm1 = ny-1
          if (j.eq.ny .and. nyc.eq.0) jp1 = 2
          do kc=2,ncz-1
            k = kc+kz*(kc-1)
            rm=(resf(im1,jm1,k-1)+resf(ip1,jm1,k-1)+resf(im1,jp1,k-1)+
     +        resf(ip1,jp1,k-1)+2.d0*(resf(im1,j,k-1)+resf(ip1,j,k-1)+
     +        resf(i,jm1,k-1)+resf(i,jp1,k-1))+4.d0*resf(i,j,k-1))
     +        *.0625d0
            rk=(resf(im1,jm1,k)+resf(ip1,jm1,k)+resf(im1,jp1,k)+
     +        resf(ip1,jp1,k)+2.d0*(resf(im1,j,k)+resf(ip1,j,k)+
     +        resf(i,jm1,k)+resf(i,jp1,k))+4.d0*resf(i,j,k))
     +        *.0625d0
            rp=(resf(im1,jm1,k+1)+resf(ip1,jm1,k+1)+resf(im1,jp1,k+1)+
     +        resf(ip1,jp1,k+1)+2.d0*(resf(im1,j,k+1)+resf(ip1,j,k+1)+
     +        resf(i,jm1,k+1)+resf(i,jp1,k+1))+4.d0*resf(i,j,k+1))
     +        *.0625d0
            rhsc(ic,jc,kc) = 0.25d0*(rm+2.d0*rk+rp)
          end do
c----+|----------------------------------------------------------------|
c         x=xa,xb; y=yc,yd; z=ze,zf cornors
c----+|----------------------------------------------------------------|
          do kc=1,ncz,ncz-1
            k = kc+kz*(kc-1)
            km1 = max0(k-1,2)
            kp1 = min0(k+1,nz-1)
            if (k.eq.1 .and. nze.eq.0) km1 = nz-1
            if (k.eq.nz .and. nzf.eq.0) kp1 = 2
            rm=(resf(im1,jm1,km1)+resf(ip1,jm1,km1)+resf(im1,jp1,km1)+
     +         resf(ip1,jp1,km1)+2.d0*(resf(im1,j,km1)+resf(ip1,j,km1)+
     +         resf(i,jm1,km1)+resf(i,jp1,km1))+
     +         4.d0*resf(i,j,km1))*.0625d0
            rk=(resf(im1,jm1,k)+resf(ip1,jm1,k)+resf(im1,jp1,k)+
     +         resf(ip1,jp1,k)+2.d0*(resf(im1,j,k)+resf(ip1,j,k)+
     +         resf(i,jm1,k)+resf(i,jp1,k))+4.d0*resf(i,j,k))*.0625d0
            rp=(resf(im1,jm1,kp1)+resf(ip1,jm1,kp1)+resf(im1,jp1,kp1)+
     +         resf(ip1,jp1,kp1)+2.d0*(resf(im1,j,kp1)+resf(ip1,j,kp1)+
     +         resf(i,jm1,kp1)+resf(i,jp1,kp1))+
     +         4.d0*resf(i,j,kp1))*.0625d0
            rhsc(ic,jc,kc) = 0.25d0*(rm+2.d0*rk+rp)
          end do
        end do
c----+|----------------------------------------------------------------|
c       x=xa,xb and z=ze,zf edges
c----+|----------------------------------------------------------------|
        do kc=1,ncz,ncz-1
          k = kc+kz*(kc-1)
          km1 = max0(k-1,2)
          kp1 = min0(k+1,nz-1)
          if (k.eq.1 .and. nze.eq.0) km1 = nz-1
          if (k.eq.nz .and. nzf.eq.0) kp1 = 2
          do jc=2,ncy-1
            j = jc+jy*(jc-1)
            rm=(resf(im1,j-1,km1)+resf(ip1,j-1,km1)+resf(im1,j+1,km1)+
     +         resf(ip1,j+1,km1)+2.d0*(resf(im1,j,km1)+resf(ip1,j,km1)+
     +         resf(i,j-1,km1)+resf(i,j+1,km1))+
     +         4.d0*resf(i,j,km1))*.0625d0
            rk=(resf(im1,j-1,k)+resf(ip1,j-1,k)+resf(im1,j+1,k)+
     +         resf(ip1,j+1,k)+2.d0*(resf(im1,j,k)+resf(ip1,j,k)+
     +         resf(i,j-1,k)+resf(i,j+1,k))+4.d0*resf(i,j,k))*.0625d0
            rp=(resf(im1,j-1,kp1)+resf(ip1,j-1,kp1)+resf(im1,j+1,kp1)+
     +         resf(ip1,j+1,kp1)+2.d0*(resf(im1,j,kp1)+resf(ip1,j,kp1)+
     +         resf(i,j-1,kp1)+resf(i,j+1,kp1))+
     +         4.d0*resf(i,j,kp1))*.0625d0
            rhsc(ic,jc,kc) = 0.25d0*(rm+2.d0*rk+rp)
          end do
        end do
      end do
c----+|----------------------------------------------------------------|
c     y boundaries y=yc and y=yd
c----+|----------------------------------------------------------------|
      do jc=1,ncy,ncy-1
        j = jc+jy*(jc-1)
        jm1 = max0(j-1,2)
        jp1 = min0(j+1,ny-1)
        if (j.eq.1 .and. nyc.eq.0) jm1 = ny-1
        if (j.eq.ny .and. nyd.eq.0) jp1 = 2
c----+|----------------------------------------------------------------|
c       (x,z) interior
c----+|----------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(i,k,ic,kc,rm,rk,rp)
!$OMP& SHARED(ix,kz,jc,jm1,j,jp1,resf,rhsc,ncx,ncz)
        do kc=2,ncz-1
          k = kc+kz*(kc-1)
          do ic=2,ncx-1
            i = ic+ix*(ic-1)
            rm=(resf(i-1,jm1,k-1)+resf(i+1,jm1,k-1)+resf(i-1,jp1,k-1)+
     +        resf(i+1,jp1,k-1)+2.d0*(resf(i-1,j,k-1)+resf(i+1,j,k-1)+
     +        resf(i,jm1,k-1)+resf(i,jp1,k-1))+
     +        4.d0*resf(i,j,k-1))*.0625d0
            rk=(resf(i-1,jm1,k)+resf(i+1,jm1,k)+resf(i-1,jp1,k)+
     +        resf(i+1,jp1,k)+2.d0*(resf(i-1,j,k)+resf(i+1,j,k)+
     +        resf(i,jm1,k)+resf(i,jp1,k))+4.d0*resf(i,j,k))*.0625d0
            rp=(resf(i-1,jm1,k+1)+resf(i+1,jm1,k+1)+resf(i-1,jp1,k+1)+
     +        resf(i+1,jp1,k+1)+2.d0*(resf(i-1,j,k+1)+resf(i+1,j,k+1)+
     +        resf(i,jm1,k+1)+resf(i,jp1,k+1))+
     +        4.d0*resf(i,j,k+1))*.0625d0
            rhsc(ic,jc,kc) = 0.25d0*(rm+2.d0*rk+rp)
          end do
        end do
!$OMP END PARALLEL DO
c----+|----------------------------------------------------------------|
c       y=yc,yd and z=ze,zf edges
c----+|----------------------------------------------------------------|
        do kc=1,ncz,ncz-1
          k = kc+kz*(kc-1)
          km1 = max0(k-1,2)
          kp1 = min0(k+1,nz-1)
          if (k.eq.1 .and. nze.eq.0) km1 = nz-1
          if (k.eq.nz .and. nzf.eq.0) kp1 = 2
c----+|----------------------------------------------------------------|
c         interior in x
c----+|----------------------------------------------------------------|
          do ic=2,ncx-1
            i = ic+ix*(ic-1)
            rm=(resf(i-1,jm1,km1)+resf(i+1,jm1,km1)+resf(i-1,jp1,km1)+
     +        resf(i+1,jp1,km1)+2.d0*(resf(i-1,j,km1)+resf(i+1,j,km1)+
     +        resf(i,jm1,km1)+resf(i,jp1,km1))+
     +        4.d0*resf(i,j,km1))*.0625d0
            rk=(resf(i-1,jm1,k)+resf(i+1,jm1,k)+resf(i-1,jp1,k)+
     +        resf(i+1,jp1,k)+2.d0*(resf(i-1,j,k)+resf(i+1,j,k)+
     +        resf(i,jm1,k)+resf(i,jp1,k))+4.d0*resf(i,j,k))*.0625d0
            rp=(resf(i-1,jm1,kp1)+resf(i+1,jm1,kp1)+resf(i-1,jp1,kp1)+
     +        resf(i+1,jp1,kp1)+2.d0*(resf(i-1,j,kp1)+resf(i+1,j,kp1)+
     +        resf(i,jm1,kp1)+resf(i,jp1,kp1))+
     +        4.d0*resf(i,j,kp1))*.0625d0
            rhsc(ic,jc,kc) = 0.25d0*(rm+2.d0*rk+rp)
          end do
        end do
      end do
c----+|----------------------------------------------------------------|
c       z=ze,zf boundaries
c----+|----------------------------------------------------------------|
      do kc=1,ncz,ncz-1
        k = kc+kz*(kc-1)
        km1 = max0(k-1,2)
        kp1 = min0(k+1,nz-1)
        if (k.eq.1 .and. nze.eq.0) km1 = nz-1
        if (k.eq.nz .and. nzf.eq.0) kp1 = 2
c----+|----------------------------------------------------------------|
c         (x,y) interior
c----+|----------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(i,j,ic,jc,rm,rk,rp)
!$OMP& SHARED(ix,jy,kc,km1,k,kp1,resf,rhsc,ncx,ncz)
        do jc=2,ncy-1
          j = jc+jy*(jc-1)
          do ic=2,ncx-1
            i = ic+ix*(ic-1)
            rm=(resf(i-1,j-1,km1)+resf(i+1,j-1,km1)+resf(i-1,j+1,km1)+
     +        resf(i+1,j+1,km1)+2.d0*(resf(i-1,j,km1)+resf(i+1,j,km1)+
     +        resf(i,j-1,km1)+resf(i,j+1,km1))+
     +        4.d0*resf(i,j,km1))*.0625d0
            rk=(resf(i-1,j-1,k)+resf(i+1,j-1,k)+resf(i-1,j+1,k)+
     +        resf(i+1,j+1,k)+2.d0*(resf(i-1,j,k)+resf(i+1,j,k)+
     +        resf(i,j-1,k)+resf(i,j+1,k))+
     +        4.d0*resf(i,j,k))*.0625d0
            rp=(resf(i-1,j-1,kp1)+resf(i+1,j-1,kp1)+resf(i-1,j+1,kp1)+
     +        resf(i+1,j+1,kp1)+2.d0*(resf(i-1,j,kp1)+resf(i+1,j,kp1)+
     +        resf(i,j-1,kp1)+resf(i,j+1,kp1))+
     +        4.d0*resf(i,j,kp1))*.0625d0
            rhsc(ic,jc,kc) = 0.25d0*(rm+2.d0*rk+rp)
          end do
        end do
!$OMP END PARALLEL DO
      end do
c----+|----------------------------------------------------------------|
c     set coarse grid residual to zero at specified boundaries
c----+|----------------------------------------------------------------|
      if (nxa.eq.1) then
        ic = 1
        do kc=1,ncz
          do jc=1,ncy
            rhsc(ic,jc,kc) = 0.0
          end do
        end do
      end if
      if (nxb.eq.1) then
        ic = ncx
        do kc=1,ncz
          do jc=1,ncy
            rhsc(ic,jc,kc) = 0.0
          end do
        end do
      end if
      if (nyc.eq.1) then
        jc = 1
        do kc=1,ncz
          do ic=1,ncx
            rhsc(ic,jc,kc) = 0.0
          end do
        end do
      end if
      if (nyd.eq.1) then
        jc = ncy
        do kc=1,ncz
          do ic=1,ncx
            rhsc(ic,jc,kc) = 0.d0
          end do
        end do
      end if
      if (nze.eq.1) then
        kc = 1
        do jc=1,ncy
          do ic=1,ncx
            rhsc(ic,jc,kc) = 0.d0
          end do
        end do
      end if
      if (nzf.eq.1) then
        kc = ncz
        do jc=1,ncy
          do ic=1,ncx
            rhsc(ic,jc,kc) = 0.d0
          end do
        end do
      end if
      return
      end subroutine res3

c----+|----------------------------------------------------------------|
c     prolon3 modified from prolon2 11/25/97
c----+|----------------------------------------------------------------|
      subroutine prolon3(ncx,ncy,ncz,p,nx,ny,nz,q,nxa,nxb,nyc,nyd,
     +                   nze,nzf,intpol)
      implicit none
      integer ncx,ncy,ncz,nx,ny,nz,intpol,nxa,nxb,nyc,nyd,nze,nzf
      double precision p(0:ncx+1,0:ncy+1,0:ncz+1),
     +q(0:nx+1,0:ny+1,0:nz+1)
      integer i,j,k,kc,ist,ifn,jst,jfn,kst,kfn,koddst,koddfn
      ist = 1
      ifn = nx
      jst = 1
      jfn = ny
      kst = 1
      kfn = nz
      koddst = 1
      koddfn = nz
      if (nxa.eq.1) then
        ist = 2
      end if
      if (nxb.eq.1) then
        ifn = nx-1
      end if
      if (nyc.eq.1) then
        jst = 2
      end if
      if (nyd.eq.1) then
        jfn = ny-1
      end if
      if (nze.eq.1) then
        kst = 2
        koddst = 3
      end if
      if (nzf.eq.1) then
        kfn = nz-1
        koddfn = nz-2
      end if

c----+|----------------------------------------------------------------|
c     linearly interpolate in z
c----+|----------------------------------------------------------------|
      if (intpol.eq.1 .or. ncz.lt.4) then
        if (ncz .lt. nz) then
c----+|----------------------------------------------------------------|
c         ncz grid is an every other point subset of nz grid
c         set odd k planes interpolating in x&y and then set even
c         k planes by averaging odd k planes
c----+|----------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(k,kc)
          do k=koddst,koddfn,2
            kc = k/2+1
            call prolon2(ncx,ncy,p(0,0,kc),nx,ny,q(0,0,k),nxa,nxb,nyc,
     +                     nyd,intpol)
          end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(i,j,k)
          do k=2,kfn,2
            do j=jst,jfn
              do i=ist,ifn
                q(i,j,k) = 0.5d0*(q(i,j,k-1)+q(i,j,k+1))
              end do
            end do
          end do
!$OMP END PARALLEL DO
c----+|----------------------------------------------------------------|
c         set periodic virtual boundaries if necessary
c----+|----------------------------------------------------------------|
          if (nze.eq.0) then
            do j=jst,jfn
              do i=ist,ifn
                q(i,j,0) = q(i,j,nz-1)
                q(i,j,nz+1) = q(i,j,2)
              end do
            end do
          end if
          return
        else
c----+|----------------------------------------------------------------|
c       ncz grid is equals nz grid so interpolate in x&y only
c----+|----------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(k,kc)
          do k=kst,kfn
            kc = k
            call prolon2(ncx,ncy,p(0,0,kc),nx,ny,q(0,0,k),nxa,nxb,nyc,
     +                     nyd,intpol)
          end do
!$OMP END PARALLEL DO
c----+|----------------------------------------------------------------|
c       set periodic virtual boundaries if necessary
c----+|----------------------------------------------------------------|
          if (nze.eq.0) then
!$OMP PARALLEL DO PRIVATE(i,j)
            do j=jst,jfn
              do i=ist,ifn
                q(i,j,0) = q(i,j,nz-1)
                q(i,j,nz+1) = q(i,j,2)
              end do
            end do
!$OMP END PARALLEL DO
          end if
          return
        end if
      else
c----+|----------------------------------------------------------------|
c       cubically interpolate in z
c----+|----------------------------------------------------------------|
        if (ncz .lt. nz) then
c----+|----------------------------------------------------------------|
c       set every other point of nz grid by interpolating in x&y
c----+|----------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(k,kc)
          do k=koddst,koddfn,2
            kc = k/2+1
            call prolon2(ncx,ncy,p(0,0,kc),nx,ny,q(0,0,k),nxa,nxb,nyc,
     +                     nyd,intpol)
          end do
!$OMP END PARALLEL DO
c----+|----------------------------------------------------------------|
c       set deep interior of nz grid using values just
c       generated and symmetric cubic interpolation in z
c----+|----------------------------------------------------------------|
          do k=4,nz-3,2
            do j=jst,jfn
              do i=ist,ifn
                q(i,j,k)=(-q(i,j,k-3)+9.d0*(q(i,j,k-1)+q(i,j,k+1))-
     +                   q(i,j,k+3))*.0625d0
              end do
            end do
          end do
c----+|----------------------------------------------------------------|
c         interpolate from q at k=2 and k=nz-1
c----+|----------------------------------------------------------------|
          if (nze.ne.0) then
c----+|----------------------------------------------------------------|
c           asymmetric formula near nonperiodic z boundaries
c----+|----------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(i,j)
            do j=jst,jfn
              do i=ist,ifn
                q(i,j,2)=(5.d0*q(i,j,1)+15.d0*q(i,j,3)-5.d0*q(i,j,5)+
     +                     q(i,j,7))*.0625d0
                q(i,j,nz-1)=(5.d0*q(i,j,nz)+15.d0*q(i,j,nz-2)-
     +                       5.d0*q(i,j,nz-4)+q(i,j,nz-6))*.0625d0
              end do
            end do
!$OMP END PARALLEL DO
          else
c----+|----------------------------------------------------------------|
c         periodicity in y alows symmetric formula near bndys
c----+|----------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(i,j)
            do j=jst,jfn
              do i=ist,ifn
                q(i,j,2) = (-q(i,j,nz-2)+9.d0*(q(i,j,1)+q(i,j,3))-
     +                       q(i,j,5))*.0625d0
                q(i,j,nz-1)=(-q(i,j,nz-4)+9.d0*(q(i,j,nz-2)+q(i,j,nz))-
     +                          q(i,j,3))*.0625d0
                q(i,j,nz+1) = q(i,j,2)
                q(i,j,0) = q(i,j,nz-1)
              end do
            end do
!$OMP END PARALLEL DO
          end if
          return
        else
c----+|----------------------------------------------------------------|
c         ncz grid is equals nx grid so interpolate in x&y only
c----+|----------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(k,kc)
          do k=kst,kfn
            kc = k
            call prolon2(ncx,ncy,p(0,0,kc),nx,ny,q(0,0,k),nxa,nxb,nyc,
     +                   nyd,intpol)
          end do
!$OMP END PARALLEL DO
c----+|----------------------------------------------------------------|
c    set periodic virtual boundaries if necessary
c----+|----------------------------------------------------------------|
          if (nze.eq.0) then
!$OMP PARALLEL DO PRIVATE(i,j)
            do j=jst,jfn
              do i=ist,ifn
                q(i,j,0) = q(i,j,nz-1)
                q(i,j,nz+1) = q(i,j,2)
              end do
            end do
!$OMP END PARALLEL DO
          end if
          return
        end if
      end if
      end subroutine prolon3

c----+|----------------------------------------------------------------|
c     add coarse grid correction in phic to fine grid approximation
c     in phif using linear or cubic interpolation
c----+|----------------------------------------------------------------|
      subroutine cor3(nx,ny,nz,phif,ncx,ncy,ncz,phic,nxa,nxb,nyc,nyd,
     +                nze,nzf,intpol,phcor)
      implicit none
      integer nx,ny,nz,ncx,ncy,ncz,nxa,nxb,nyc,nyd,nze,nzf,intpol
      integer i,j,k,ist,ifn,jst,jfn,kst,kfn
      double precision phif(0:nx+1,0:ny+1,0:nz+1),
     +phic(0:ncx+1,0:ncy+1,0:ncz+1)
      double precision phcor(0:nx+1,0:ny+1,0:nz+1)

      phcor = 0.d0

c----+|----------------------------------------------------------------|
c    lift correction in phic to fine grid in phcor
c----+|----------------------------------------------------------------|
      call prolon3(ncx,ncy,ncz,phic,nx,ny,nz,phcor,nxa,nxb,nyc,nyd,
     +             nze,nzf,intpol)

c----+|----------------------------------------------------------------|
c    add correction in phcor to phif on nonspecified boundaries
c----+|----------------------------------------------------------------|
      ist = 1
      ifn = nx
      jst = 1
      jfn = ny
      kst = 1
      kfn = nz
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.1) jst = 2
      if (nyd.eq.1) jfn = ny-1
      if (nze.eq.1) kst = 2
      if (nzf.eq.1) kfn = nz-1
!$OMP PARALLEL DO PRIVATE(i,j,k)
      do k=kst,kfn
        do j=jst,jfn
          do i=ist,ifn
            phif(i,j,k) = phif(i,j,k) + phcor(i,j,k)
          end do
        end do
      end do
!$OMP END PARALLEL DO

c----+|----------------------------------------------------------------|
c    add periodic points if necessary
c----+|----------------------------------------------------------------|
      if (nze.eq.0) then
        do j=jst,jfn
          do i=ist,ifn
            phif(i,j,0) = phif(i,j,nz-1)
            phif(i,j,nz+1) = phif(i,j,2)
          end do
        end do
      end if
      if (nyc.eq.0) then
        do k=kst,kfn
          do i=ist,ifn
            phif(i,0,k) = phif(i,ny-1,k)
            phif(i,ny+1,k) = phif(i,2,k)
          end do
        end do
      end if
      if (nxa.eq.0) then
        do k=kst,kfn
          do j=jst,jfn
            phif(0,j,k) = phif(nx-1,j,k)
            phif(nx+1,j,k) = phif(2,j,k)
          end do
        end do
      end if
      end subroutine cor3

c----+|----------------------------------------------------------------|
c    set virtual periodic boundaries from interior values
c    in three dimensions (for all 3-d solvers)
c----+|----------------------------------------------------------------|
      subroutine per3vb(nx,ny,nz,phi,nxa,nyc,nze)
      implicit none
      integer nx,ny,nz,nxa,nyc,nze,j,k,i
      double precision phi(0:nx+1,0:ny+1,0:nz+1)
      if (nxa.eq.0) then
        do k=1,nz
          do j=1,ny
            phi(0,j,k) = phi(nx-1,j,k)
            phi(nx,j,k) = phi(1,j,k)
            phi(nx+1,j,k) = phi(2,j,k)
          end do
        end do
      end if
      if (nyc.eq.0) then
        do k=1,nz
          do i=1,nx
            phi(i,0,k) = phi(i,ny-1,k)
            phi(i,ny,k) = phi(i,1,k)
            phi(i,ny+1,k) = phi(i,2,k)
          end do
        end do
      end if
      if (nze.eq.0) then
        do j=1,ny
          do i=1,nx
            phi(i,j,0) = phi(i,j,nz-1)
            phi(i,j,nz) = phi(i,j,1)
            phi(i,j,nz+1) = phi(i,j,2)
          end do
        end do
      end if
      return
      end subroutine per3vb

c----+|----------------------------------------------------------------|
c     compute mixed partial derivative approximations
c----+|----------------------------------------------------------------|
      subroutine pde2cr(nx,ny,u,i,j,ux3y,uxy3,ux2y2)
      implicit none
      integer nx,ny,i,j,n1,n2,n3,n4,m1,m2,m3,m4
      double precision u(nx,ny),ux3y,uxy3,ux2y2
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      double precision xa,xb,yc,yd,tolmax,relmax
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
      double precision dlx,dly,dyox,dxoy,dlx2,dly2,dlxx,dlxy,dlyy,dlxy2,
     +             dlxy4,dxxxy4,dxyyy4,dxxyy,tdlx3,tdly3,dlx4,dly4,
     +             dlxxx,dlyyy
      common/com2dcr/dyox,dxoy,dlx2,dly2,dlxy,dlxy2,dlxy4,
     +               dxxxy4,dxyyy4,dxxyy,dlxxx,dlyyy
      common/pde2com/dlx,dly,dlxx,dlyy,tdlx3,tdly3,dlx4,dly4
      n1=ny-1
      n2=ny-2
      n3=ny-3
      n4=ny-4
      m1=nx-1
      m2=nx-2
      m3=nx-3
      m4=nx-4

      if (i.eq.1) then

      if ((j.gt.2.and.j.lt.ny-1)) then
c----+|----------------------------------------------------------------|
c      x=xa, yinterior
c----+|----------------------------------------------------------------|
        ux3y=(5.d0*u(1,j-1)-18.d0*u(2,j-1)+24.d0*u(3,j-1)-14.d0*u(4,j-1)
     +       +3.d0*u(5,j-1)-5.d0*u(1,j+1)+18.d0*u(2,j+1)-24.d0*u(3,j+1)
     +       +14.d0*u(4,j+1)-3.d0*u(5,j+1))/dxxxy4
        uxy3=(3.d0*u(1,j-2)-4.d0*u(2,j-2)+u(3,j-2)
     +       -6.d0*u(1,j-1)+8.d0*u(2,j-1)-2.d0*u(3,j-1)
     +       +6.d0*u(1,j+1)-8.d0*u(2,j+1)+2.d0*u(3,j+1)
     +       -3.d0*u(1,j+2)+4.d0*u(2,j+2)-u(3,j+2))/dxyyy4
      else if (j.eq.1) then
c----+|----------------------------------------------------------------|
c       (xa,yc)
c----+|----------------------------------------------------------------|
        ux3y=(15.d0*u(1,1)-54.d0*u(2,1)+72.d0*u(3,1)-42.d0*u(4,1)
     +       +9.d0*u(5,1)-20.d0*u(1,2)+72.d0*u(2,2)-96.d0*u(3,2)
     +       +56.d0*u(4,2)-12.d0*u(5,2)+5.d0*u(1,3)-18.d0*u(2,3)
     +       +24.d0*u(3,3)-14.d0*u(4,3)+3.d0*u(5,3))/dxxxy4
        uxy3=(15.d0*u(1,1)-20.d0*u(2,1)+5.d0*u(3,1)
     +       -54.d0*u(1,2)+72.d0*u(2,2)-18.d0*u(3,2)
     +       +72.d0*u(1,3)-96.d0*u(2,3)+24.d0*u(3,3)
     +       -42.d0*u(1,4)+56.d0*u(2,4)-14.d0*u(3,4)
     +       +9.d0*u(1,5)-12.d0*u(2,5)+3.d0*u(3,5))/dxyyy4
        ux2y2=(4.d0*u(1,1)-10.d0*u(2,1)+8.d0*u(3,1)-2.d0*u(4,1)
     +       -10.d0*u(1,2)+25.d0*u(2,2)-20.d0*u(3,2)+5.d0*u(4,2)
     +       +8.d0*u(1,3)-20.d0*u(2,3)+16.d0*u(3,3)-4.d0*u(4,3)
     +       -2.d0*u(1,4)+5.d0*u(2,4)-4.d0*u(3,4)+u(4,4))/dxxyy
      else if (j.eq.2) then
c----+|----------------------------------------------------------------|
c       (xa,yc+dly)
c----+|----------------------------------------------------------------|
        ux3y=(5.d0*u(1,1)-18.d0*u(2,1)+24.d0*u(3,1)-14.d0*u(4,1)
     +       +3.d0*u(5,1)-5.d0*u(1,3)+18.d0*u(2,3)-24.d0*u(3,3)
     +       +14.d0*u(4,3)-3.d0*u(5,3))/dxxxy4
        uxy3=(9.d0*u(1,1)-12.d0*u(2,1)+3.d0*u(3,1)
     +       -30.d0*u(1,2)+40.d0*u(2,2)-10.d0*u(3,2)
     +       +36.d0*u(1,3)-48.d0*u(2,3)+12.d0*u(3,3)
     +       -18.d0*u(1,4)+24.d0*u(2,4)-6.d0*u(3,4)
     +       +3.d0*u(1,5)-4.d0*u(2,5)+u(3,5))/dxyyy4
      else if (j.eq.ny-1) then
c----+|----------------------------------------------------------------|
c       x=xa,y=yd-dly
c----+|----------------------------------------------------------------|
        ux3y=(5.d0*u(1,j-1)-18.d0*u(2,j-1)+24.d0*u(3,j-1)
     +       -14.d0*u(4,j-1)+3.d0*u(5,j-1)-5.d0*u(1,j+1)+18.d0*u(2,j+1)
     +       -24.d0*u(3,j+1)+14.d0*u(4,j+1)-3.d0*u(5,j+1))
        uxy3=(5.d0*u(1,n2)-18.d0*u(2,n2)+24.d0*u(3,n2)
     +       -14.d0*u(4,n2)+3.d0*u(5,n2)-5.d0*u(1,ny)+18.d0*u(2,ny)
     +       -24.d0*u(3,ny)+14.d0*u(4,ny)-3.d0*u(5,ny))/dxyyy4
      else if (j.eq.ny) then
c----+|----------------------------------------------------------------|
c       x=xa, y=yd
c----+|----------------------------------------------------------------|
        ux3y=(-5.d0*u(1,n2)+18.d0*u(2,n2)-24.d0*u(3,n2)+14.d0*u(4,n2)
     +       -3.d0*u(5,n2)+20.d0*u(1,n1)-72.d0*u(2,n1)+96.d0*u(3,n1)
     +       -56.d0*u(4,n1)+12.d0*u(5,n1)-15.d0*u(1,ny)+54.d0*u(2,ny)
     +       -72.d0*u(3,ny)+42.d0*u(4,ny)-9.d0*u(5,ny))/dxxxy4
        uxy3=(-9.d0*u(1,n4)+12.d0*u(2,n4)-3.d0*u(3,n4)
     +       +42.d0*u(1,n3)-56.d0*u(2,n3)+14.d0*u(3,n3)
     +       -72.d0*u(1,n2)+96.d0*u(2,n2)-24.d0*u(3,n2)
     +       +54.d0*u(1,n1)-72.d0*u(2,n1)+18.d0*u(3,n1)
     +       -15.d0*u(1,ny)+20.d0*u(2,ny)-5.d0*u(3,ny))/dxyyy4
        ux2y2=(-2.d0*u(1,n3)+5.d0*u(2,n3)-4.d0*u(3,n3)+u(4,n3)
     +       +8.d0*u(1,n2)-20.d0*u(2,n2)+16.d0*u(3,n2)-4.d0*u(4,n2)
     +       -10.d0*u(1,n1)+25.d0*u(2,n1)-20.d0*u(3,n1)+5.d0*u(4,n1)
     +       +4.d0*u(1,ny)-10.d0*u(2,ny)+8.d0*u(3,ny)-2.d0*u(4,ny))
     +       /dxxyy
      end if

      else if (i.eq.2) then

        if ((j.gt.2.and.j.lt.ny-1)) then
c----+|----------------------------------------------------------------|
c         x=xa+dlx, y interior
c----+|----------------------------------------------------------------|
          ux3y=(3.d0*u(1,j-1)-10.d0*u(2,j-1)+12.d0*u(3,j-1)
     +        -6.d0*u(4,j-1)+u(5,j-1)-3.d0*u(1,j+1)+10.d0*u(2,j+1)
     +        -12.d0*u(3,j+1)+6.d0*u(4,j+1)-u(5,j+1))/dxxxy4
          uxy3=(u(1,j-2)-u(3,j-2)-2*u(1,j-1)+2.d0*u(3,j-1)
     +        +2.d0*u(1,j+1)-2.d0*u(3,j+1)-u(1,j+2)+u(3,j+2))/dxyyy4
        else if (j.eq.1) then
c----+|----------------------------------------------------------------|
c         x=xa+dlx, y=yc
c----+|----------------------------------------------------------------|
          ux3y=(9.d0*u(1,1)-30.d0*u(2,1)+36.d0*u(3,1)-18.d0*u(4,1)
     +         +3.d0*u(5,1)-12.d0*u(1,2)+40.d0*u(2,2)-48.d0*u(3,2)
     +         +24.d0*u(4,2)-4.d0*u(5,2)+3.d0*u(1,3)-10.d0*u(2,3)
     +         +12.d0*u(3,3)-6.d0*u(4,3)+u(5,3))/dxxxy4
          uxy3=(5.d0*u(1,1)-5.d0*u(3,1)-18.d0*u(1,2)+18.d0*u(3,2)
     +         +24.d0*u(1,3)-24.d0*u(3,3)-14.d0*u(1,4)
     +         +14.d0*u(3,4)+3.d0*u(1,5)-3.d0*u(3,5))/dxyyy4

        else if (j.eq.2) then
c----+|----------------------------------------------------------------|
c         x=xa+dlx,y=yc+dly
c----+|----------------------------------------------------------------|
          ux3y=(3.d0*u(1,1)-10.d0*u(2,1)+12.d0*u(3,1)-6.d0*u(4,1)+u(5,1)
     +         -3.d0*u(1,3)+10.d0*u(2,3)-12.d0*u(3,3)+6.d0*u(4,3)
     +         -u(5,3))/dxxxy4
          uxy3=(3.d0*u(1,1)-3.d0*u(3,1)-10.d0*u(1,2)+10.d0*u(3,2)
     +         +12.d0*u(1,3)-12.d0*u(3,3)-6.d0*u(1,4)+6.d0*u(3,4)
     +         +u(1,5)-u(3,5))/dxyyy4

        else if (j.eq.ny-1) then
c----+|----------------------------------------------------------------|
c       x=xa+dlx,y=yd-dly
c----+|----------------------------------------------------------------|
          ux3y=(3.d0*u(1,n2)-10.d0*u(2,n2)+12.d0*u(3,n2)-6.d0*u(4,n2)
     +         +u(5,n2)-3.d0*u(1,ny)+10.d0*u(2,ny)-12.d0*u(3,ny)
     +         +6.d0*u(4,ny)-u(5,ny))/dxxxy4
          uxy3=(-u(1,n4)+u(3,n4)+6.d0*u(1,n3)-6.d0*u(3,n3)
     +        -12.d0*u(1,n2)+12.d0*u(3,n2)+10.d0*u(1,n1)-10.d0*u(3,n1)
     +        -3.d0*u(1,ny)+3.d0*u(3,ny))/dxyyy4

        else if (j.eq.ny) then
c----+|----------------------------------------------------------------|
c         x=xa+dlx,y=yd
c----+|----------------------------------------------------------------|
          ux3y=(-3.d0*u(1,n2)+10.d0*u(2,n2)-12.d0*u(3,n2)+6.d0*u(4,n2)
     +         -u(5,n2)+12.d0*u(1,n1)-40.d0*u(2,n1)+48.d0*u(3,n1)
     +         -24.d0*u(4,n1)+4*u(5,n1)-9.d0*u(1,ny)+30.d0*u(2,ny)
     +         -36.d0*u(3,ny)+18.d0*u(4,ny)-3.d0*u(5,ny))/dxxxy4
          uxy3=(-3.d0*u(1,n4)+3.d0*u(3,n4)+14.d0*u(1,n3)-14.d0*u(3,n3)
     +         -24.d0*u(1,n2)+24.d0*u(3,n2)+18.d0*u(1,n1)-18.d0*u(3,n1)
     +         -5.d0*u(1,ny)+5.d0*u(3,ny))/dxyyy4
        end if

      else if (i.gt.2 .and. i.lt.nx-1) then

        if (j.eq.1) then
c----+|----------------------------------------------------------------|
c         y=yc,x interior
c----+|----------------------------------------------------------------|
          ux3y=(3.d0*u(i-2,1)-6.d0*u(i-1,1)+6.d0*u(i+1,1)-3.d0*u(i+2,1)
     +         -4.d0*u(i-2,2)+8.d0*u(i-1,2)-8.d0*u(i+1,2)+4.d0*u(i+2,2)
     +         +u(i-2,3)-2.d0*u(i-1,3)+2.d0*u(i+1,3)-u(i+2,3))/dxxxy4
          uxy3=(5.d0*u(i-1,1)-5.d0*u(i+1,1)-18.d0*u(i-1,2)
     +         +18.d0*u(i+1,2)+24.d0*u(i-1,3)-24.d0*u(i+1,3)
     +         -14.d0*u(i-1,4)+14.d0*u(i+1,4)+3.d0*u(i-1,5)
     +         -3.d0*u(i+1,5))/dxyyy4
        else if (j.eq.2) then
c----+|----------------------------------------------------------------|
c         y=yc+dly,x interior
c----+|----------------------------------------------------------------|
          ux3y=(u(i-2,1)-2.d0*u(i-1,1)+2.d0*u(i+1,1)-u(i+2,1)
     +         -u(i-2,3)+2.d0*u(i-1,3)-2.d0*u(i+1,3)+u(i+2,3))/dxxxy4
          uxy3=(u(i-1,1)-u(i+1,1)-2.d0*u(i-1,2)+2.d0*u(i+1,2)
     +         +2.d0*u(i-1,4)-2.d0*u(i+1,4)-u(i-1,5)+u(i+1,5))/dxyyy4
        else if (j.eq.ny-1) then
c----+|----------------------------------------------------------------|
c         y=yd-dly, x interior
c----+|----------------------------------------------------------------|
          ux3y=(u(i-2,n2)-2.d0*u(i-1,n2)+2.d0*u(i+1,n2)-u(i+2,n2)
     +        -u(i-2,ny)+2.d0*u(i-1,ny)-2.d0*u(i+1,ny)+u(i+2,ny))/dxxxy4
          uxy3=(-u(i-1,n4)+u(i+1,n4)+6.d0*u(i-1,n3)-6.d0*u(i+1,n3)
     +        -12.d0*u(i-1,n2)+12.d0*u(i+1,n2)+10.d0*u(i-1,n1)
     +        -10.d0*u(i+1,n1)-3.d0*u(i-1,ny)+3.d0*u(i+1,ny))/dxyyy4
        else if (j.eq.ny) then
c----+|----------------------------------------------------------------|
c        y=yd, x interior
c----+|----------------------------------------------------------------|
         ux3y=(-u(i-2,n2)+2.d0*u(i-1,n2)-2.d0*u(i+1,n2)+u(i+2,n2)
     +       +4.d0*u(i-2,n1)-8.d0*u(i-1,n1)+8.d0*u(i+1,n1)
     +       -4.d0*u(i+2,n1)-3.d0*u(i-2,ny)+6.d0*u(i-1,ny)
     +       -6.d0*u(i+1,ny)+3.d0*u(i+2,ny))/dxxxy4
         uxy3=(-3.d0*u(i-1,n4)+3.d0*u(i+1,n4)+14.d0*u(i-1,n3)
     +       -14.d0*u(i+1,n3)-24.d0*u(i-1,n2) +24.d0*u(i+1,n2)
     +       +18.d0*u(i-1,n1)-18.d0*u(i+1,n1)-5.d0*u(i-1,ny)
     +       +5.d0*u(i+1,ny))/dxyyy4
        end if

      else if (i.eq.nx-1) then

        if ((j.gt.2.and.j.lt.ny-1)) then
c----+|----------------------------------------------------------------|
c         x=xb-dlx,y interior
c----+|----------------------------------------------------------------|
          ux3y=(-u(m4,j-1)+6.d0*u(m3,j-1)-12.d0*u(m2,j-1)
     +         +10.d0*u(m1,j-1)-3.d0*u(nx,j-1)+u(m4,j+1)
     +         -6.d0*u(m3,j+1)+12.d0*u(m2,j+1)-10.d0*u(m1,j+1)
     +         +3.d0*u(nx,j+1)) /dxxxy4
          uxy3=(u(m2,j-2)-u(nx,j-2)-2.d0*u(m2,j-1)+2.d0*u(nx,j-1)
     +        +2.d0*u(m2,j+1)-2.d0*u(nx,j+1)-u(m2,j+2)+u(nx,j+2))/dxyyy4
        else if (j.eq.1) then
c----+|----------------------------------------------------------------|
c       x=xb-dlx, y=yc
c----+|----------------------------------------------------------------|
          ux3y=(-3.d0*u(m4,1)+18.d0*u(m3,1)-36.d0*u(m2,1)+30.d0*u(m1,1)
     +         -9.d0*u(nx,1)+4.d0*u(m4,2)-24.d0*u(m3,2)+48.d0*u(m2,2)
     +         -40.d0*u(m1,2)+12.d0*u(nx,2)-u(m4,3)+6.d0*u(m3,3)
     +         -12.d0*u(m2,3)+10.d0*u(m1,3)-3.d0*u(nx,3))/dxxxy4
          uxy3=(5.d0*u(m2,1)-5.d0*u(nx,1)-18.d0*u(m2,2)+18.d0*u(nx,2)
     +         +24.d0*u(m2,3)-24.d0*u(nx,3)-14.d0*u(m2,4)+14.d0*u(nx,4)
     +         +3.d0*u(m2,5)-3.d0*u(nx,5))/dxyyy4
        else if (j.eq.2) then
c----+|----------------------------------------------------------------|
c         x=xb-dlx,y=yc+dly
c----+|----------------------------------------------------------------|
          ux3y=(-u(m4,1)+6.d0*u(m3,1)-12.d0*u(m2,1)+10.d0*u(m1,1)
     +        -3.d0*u(nx,1)+u(m4,3)-6.d0*u(m3,3)+12.d0*u(m2,3)
     +        -10.d0*u(m1,3)+3.d0*u(nx,3))/dxxxy4
          uxy3=(3.d0*u(m2,1)-3.d0*u(nx,1)-10.d0*u(m2,2)+10.d0*u(nx,2)
     +        +12.d0*u(m2,3)-12.d0*u(nx,3)-6.d0*u(m2,4)+6.d0*u(nx,4)
     +        +u(m2,5)-u(nx,5)) / dxyyy4
        else if (j.eq.ny-1) then
c----+|----------------------------------------------------------------|
c         x=xb-dlx,y=yd-dly
c----+|----------------------------------------------------------------|
          ux3y=(-u(m4,n2)+6.d0*u(m3,n2)-12.d0*u(m2,n2)+10.d0*u(m1,n2)
     +        -3.d0*u(nx,n2)+u(m4,ny)-6.d0*u(m3,ny)+12.d0*u(m2,ny)
     +        -10.d0*u(m1,ny)+3.d0*u(nx,ny))/dxxxy4
          uxy3=(-u(m2,n4)+u(nx,n4)+6*u(m2,n3)-6.d0*u(nx,n3)
     +        -12.d0*u(m2,n2)+12.d0*u(nx,n2)+10.d0*u(m2,n1)
     +        -10.d0*u(nx,n1)-3.d0*u(m2,ny)+3.d0*u(nx,ny)) / dxyyy4
        else if (j.eq.ny) then
c----+|----------------------------------------------------------------|
c         x=xb.dlx,y=yd
c----+|----------------------------------------------------------------|
          ux3y=(u(m4,n2)-6.d0*u(m3,n2)+12.d0*u(m2,n2)-10.d0*u(m1,n2)
     +        +3.d0*u(nx,n2)-4.d0*u(m4,n1)+24.d0*u(m3,n1)-48.d0*u(m2,n1)
     +        +40.d0*u(m1,n1)-12.d0*u(nx,n1)+3.d0*u(m4,ny)
     +        -18.d0*u(m3,ny)+36.d0*u(m2,ny)-30.d0*u(m1,ny)
     +        +9.d0*u(nx,ny))/ dxxxy4
          uxy3=(-3.d0*u(m2,n4)+3.d0*u(nx,n4)+14.d0*u(m2,n3)
     +         -14.d0*u(nx,n3)-24.d0*u(m2,n2)+24.d0*u(nx,n2)
     +         +18.d0*u(m2,n1)-18.d0*u(nx,n1)
     +         -5.d0*u(m2,ny)+5.d0*u(nx,ny)) / dxyyy4
        end if

      else if (i.eq.nx) then

        if ((j.gt.2.and.j.lt.ny-1)) then
c----+|----------------------------------------------------------------|
c         x=xb,y interior
c----+|----------------------------------------------------------------|
          ux3y=(-3.d0*u(m4,j-1)+14.d0*u(m3,j-1)-24.d0*u(m2,j-1)
     +         +18.d0*u(m1,j-1)-5.d0*u(nx,j-1)+3.d0*u(m4,j+1)
     +         -14.d0*u(m3,j+1)+24.d0*u(m2,j+1)-18.d0*u(m1,j+1)
     +         +5.d0*u(nx,j+1)) / dxxxy4
          uxy3=(-u(m2,j-2)+4.d0*u(m1,j-2)-3.d0*u(nx,j-2)
     +         +2.d0*u(m2,j-1)-8.d0*u(m1,j-1)+6.d0*u(nx,j-1)
     +         -2.d0*u(m2,j+1)+8.d0*u(m1,j+1)-6.d0*u(nx,j+1)
     +         +u(m2,j+2)-4.d0*u(m1,j+2)+3.d0*u(nx,j+2)) / dxyyy4
        else if (j.eq.1) then
c----+|----------------------------------------------------------------|
c         x=xb,y=yc
c----+|----------------------------------------------------------------|
          ux3y=(-9.d0*u(m4,1)+42.d0*u(m3,1)-72.d0*u(m2,1)+54.d0*u(m1,1)
     +         -15.d0*u(nx,1)+12.d0*u(m4,2)-56.d0*u(m3,2)+96.d0*u(m2,2)
     +         -72.d0*u(m1,2)+20.d0*u(nx,2)-3.d0*u(m4,3)+14.d0*u(m3,3)
     +         -24.d0*u(m2,3)+18.d0*u(m1,3)-5.d0*u(nx,3))/dxxxy4
          uxy3=(-5.d0*u(m2,1)+20.d0*u(m1,1)-15.d0*u(nx,1)
     +         +18.d0*u(m2,2)-72.d0*u(m1,2)+54.d0*u(nx,2)
     +         -24.d0*u(m2,3)+96.d0*u(m1,3)-72.d0*u(nx,3)
     +         +14.d0*u(m2,4)-56.d0*u(m1,4)+42.d0*u(nx,4)
     +         -3.d0*u(m2,5)+12.d0*u(m1,5)-9.d0*u(nx,5)) / dxyyy4
          ux2y2=(-2.d0*u(m3,1)+8.d0*u(m2,1)-10.d0*u(m1,1)+4.d0*u(nx,1)
     +         +5.d0*u(m3,2)-20.d0*u(m2,2)+25.d0*u(m1,2)-10.d0*u(nx,2)
     +         -4.d0*u(m3,3)+16.d0*u(m2,3)-20.d0*u(m1,3)+8.d0*u(nx,3)
     +         +u(m3,4)-4.d0*u(m2,4)+5.d0*u(m1,4)-2.d0*u(nx,4)) / dxxyy
        else if (j.eq.2) then
c----+|----------------------------------------------------------------|
c         x=xb,y=yc+dly
c----+|----------------------------------------------------------------|
          ux3y=(-3.d0*u(m4,1)+14.d0*u(m3,1)-24.d0*u(m2,1)+18.d0*u(m1,1)
     +        -5.d0*u(nx,1)+3.d0*u(m4,3)-14.d0*u(m3,3)+24.d0*u(m2,3)
     +        -18.d0*u(m1,3)+5.d0*u(nx,3))/ dxxxy4
          uxy3=(-3.d0*u(m2,1)+12.d0*u(m1,1)-9.d0*u(nx,1)
     +        +10.d0*u(m2,2)-40.d0*u(m1,2)+30.d0*u(nx,2)
     +        -12.d0*u(m2,3)+48.d0*u(m1,3)-36.d0*u(nx,3)
     +        +6.d0*u(m2,4)-24.d0*u(m1,4)+18.d0*u(nx,4)
     +        -u(m2,5)+4.d0*u(m1,5)-3.d0*u(nx,5)) / dxyyy4
        else if (j.eq.ny-1) then
c----+|----------------------------------------------------------------|
c         x=xb,y=yd-dly
c----+|----------------------------------------------------------------|
          ux3y=(-3.d0*u(m4,n2)+14.d0*u(m3,n2)-24.d0*u(m2,n2)
     +        +18.d0*u(m1,n2)-5.d0*u(nx,n2)+3.d0*u(m4,ny)-14.d0*u(m3,ny)
     +        +24.d0*u(m2,ny)-18.d0*u(m1,ny)+5.d0*u(nx,ny)) / dxxxy4
          uxy3=(u(m2,n4)-4.d0*u(m1,n4)+3.d0*u(nx,n4)
     +        -6.d0*u(m2,n3)+24.d0*u(m1,n3)-18.d0*u(nx,n3)
     +        +12.d0*u(m2,n2)-48.d0*u(m1,n2)+36.d0*u(nx,n2)
     +        -10.d0*u(m2,n1)+40.d0*u(m1,n1)-30.d0*u(nx,n1)
     +        +3.d0*u(m2,ny)-12.d0*u(m1,ny)+9.d0*u(nx,ny)) / dxyyy4
        else if (j.eq.ny) then
c----+|----------------------------------------------------------------|
c         x=xb,y=yd
c----+|----------------------------------------------------------------|
          ux3y=(3.d0*u(m4,n2)-14.d0*u(m3,n2)+24.d0*u(m2,n2)
     +         -18.d0*u(m1,n2)+5.d0*u(nx,n2)-12.d0*u(m4,n1)
     +         +56.d0*u(m3,n1)-96.d0*u(m2,n1)+72.d0*u(m1,n1)
     +         -20.d0*u(nx,n1)+9.d0*u(m4,ny)-42.d0*u(m3,ny)
     +         +72.d0*u(m2,ny)-54.d0*u(m1,ny)+15.d0*u(nx,ny)) / dxxxy4
          uxy3=(3.d0*u(m2,n4)-12.d0*u(m1,n4)+9.d0*u(nx,n4)
     +         -14.d0*u(m2,n3)+56.d0*u(m1,n3)-42.d0*u(nx,n3)
     +         +24.d0*u(m2,n2)-96.d0*u(m1,n2)+72.d0*u(nx,n2)
     +         -18.d0*u(m2,n1)+72.d0*u(m1,n1)-54.d0*u(nx,n1)
     +         +5.d0*u(m2,ny)-20.d0*u(m1,ny)+15.d0*u(nx,ny)) / dxyyy4
          ux2y2=(u(m3,n3)-4.d0*u(m2,n3)+5.d0*u(m1,n3)-2.d0*u(nx,n3)
     +         -4.d0*u(m3,n2)+16.d0*u(m2,n2)-20.d0*u(m1,n2)
     +         +8.d0*u(nx,n2)+5.d0*u(m3,n1)-20.d0*u(m2,n1)
     +         +25.d0*u(m1,n1)-10.d0*u(nx,n1)-2.d0*u(m3,ny)
     +         +8.d0*u(m2,ny)-10.d0*u(m1,ny)+4.d0*u(nx,ny)) / dxxyy
        end if
      end if

      return
      end subroutine pde2cr
c----+|----------------------------------------------------------------|

c----+|----------------------------------------------------------------|
c     estimate third and fourth partial derivatives in x,y,z
c----+|----------------------------------------------------------------|
      subroutine pde3(nx,ny,nz,u,i,j,k,ux3,ux4,uy3,uy4,uz3,uz4,
     +                nxa,nyc,nze)
      implicit none
      integer nx,ny,nz,i,j,k,nxa,nyc,nze
      double precision u(nx,ny,nz)
      double precision dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,
     +dlx4,dly4,dlz4
      common/pde3com/dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,
     +               dlx4,dly4,dlz4
      double precision ux3,ux4,uy3,uy4,uz3,uz4
c----+|----------------------------------------------------------------|
c    x,y partial derivatives
c----+|----------------------------------------------------------------|
      call p3de2(nx,ny,u(1,1,k),i,j,ux3,ux4,uy3,uy4,nxa,nyc)
c----+|----------------------------------------------------------------|
c    z partial derivatives
c----+|----------------------------------------------------------------|
      if (nze.ne.0) then
c----+|----------------------------------------------------------------|
c    nonperiodic in z
c----+|----------------------------------------------------------------|
      if(k.gt.2 .and. k.lt.nz-1) then
      uz3=(-u(i,j,k-2)+2.d0*u(i,j,k-1)-2.d0*u(i,j,k+1)+
     +     u(i,j,k+2))/tdlz3
      uz4=(u(i,j,k-2)-4.d0*u(i,j,k-1)+6.d0*u(i,j,k)-4.d0*u(i,j,k+1)+
     +         u(i,j,k+2))/dlz4
      else if (k.eq.1) then
      uz3=(-5.d0*u(i,j,1)+18.d0*u(i,j,2)-24.d0*u(i,j,3)+14.d0*u(i,j,4)-
     +          3.d0*u(i,j,5))/tdlz3
      uz4 = (3.d0*u(i,j,1)-14.d0*u(i,j,2)+26.d0*u(i,j,3)-
     +      24.d0*u(i,j,4)+11.d0*u(i,j,5)-2.d0*u(i,j,6))/dlz4
      else if (k.eq.2) then
      uz3 = (-3.d0*u(i,j,1)+10.d0*u(i,j,2)-12.d0*u(i,j,3)+
     +      6.d0*u(i,j,4)-u(i,j,5))/tdlz3
      uz4 = (2.d0*u(i,j,1)-9.d0*u(i,j,2)+16.d0*u(i,j,3)-
     +      14.d0*u(i,j,4)+6.d0*u(i,j,5)-u(i,j,6))/dlz4
      else if (k.eq.nz-1) then
      uz3 = (u(i,j,nz-4)-6.d0*u(i,j,nz-3)+12.d0*u(i,j,nz-2)-10.d0*
     +           u(i,j,nz-1)+3.d0*u(i,j,nz))/tdlz3
      uz4 = (-u(i,j,nz-5)+6.d0*u(i,j,nz-4)-14.d0*u(i,j,nz-3)+16.d0*
     +            u(i,j,nz-2)-9.d0*u(i,j,nz-1)+2.d0*u(i,j,nz))/dlz4
      else if (k.eq.nz) then
      uz3 = (3.d0*u(i,j,nz-4)-14.d0*u(i,j,nz-3)+24.d0*u(i,j,nz-2)
     +      -18.d0*u(i,j,nz-1)+5.d0*u(i,j,nz))/tdlz3
      uz4 = (-2.d0*u(i,j,nz-5)+11.d0*u(i,j,nz-4)-24.d0*u(i,j,nz-3)+
     +       26.d0*u(i,j,nz-2)-14.d0*u(i,j,nz-1)+3.d0*u(i,j,nz))/dlz4
      end if
      else
c----+|----------------------------------------------------------------|
c    periodic in z so use symmetric formula even "near" z boundaies
c----+|----------------------------------------------------------------|
      if(k.gt.2 .and. k.lt.nz-1) then
      uz3=(-u(i,j,k-2)+2.d0*u(i,j,k-1)-2.d0*u(i,j,k+1)+
     +     u(i,j,k+2))/tdlz3
      uz4=(u(i,j,k-2)-4.d0*u(i,j,k-1)+6.d0*u(i,j,k)-4.d0*u(i,j,k+1)+
     +     u(i,j,k+2))/dlz4
      else if (k.eq.1) then
      uz3 = (-u(i,j,nz-2)+2.d0*u(i,j,nz-1)-2.d0*u(i,j,2)+
     +       u(i,j,3))/tdlz3
      uz4 = (u(i,j,nz-2)-4.d0*u(i,j,nz-1)+6.d0*u(i,j,1)-4.d0*u(i,j,2)+
     +       u(i,j,3))/dlz4
      else if (k.eq.2) then
      uz3 = (-u(i,j,nz-1)+2.d0*u(i,j,1)-2.d0*u(i,j,3)+u(i,j,4))/(tdlz3)
      uz4 = (u(i,j,nz-1)-4.d0*u(i,j,1)+6.d0*u(i,j,2)-4.d0*u(i,j,3)+
     +       u(i,j,4))/dlz4
      else if (k.eq.nz-1) then
      uz3 = (-u(i,j,nz-3)+2.d0*u(i,j,nz-2)-2.d0*u(i,j,1)+u(i,j,2))/
     +      tdlz3
      uz4 = (u(i,j,nz-3)-4.d0*u(i,j,nz-2)+6.d0*u(i,j,nz-1)
     +      -4.d0*u(i,j,1)+u(i,j,2))/ dlz4
      else if (k.eq.nz) then
      uz3 = (-u(i,j,nz-2)+2.d0*u(i,j,nz-1)-2.d0*u(i,j,2)+
     +      u(i,j,3))/tdlz3
      uz4 = (u(i,j,nz-2)-4.d0*u(i,j,nz-1)+6.d0*u(i,j,nz)-4.d0*u(i,j,2)+
     +       u(i,j,3))/dlz4
      end if
      end if
      return
      end subroutine pde3

      subroutine p3de2(nx,ny,u,i,j,ux3,ux4,uy3,uy4,nxa,nyc)
c----+|----------------------------------------------------------------|
c    third and fourth partial derivatives in x and y
c----+|----------------------------------------------------------------|
      implicit none
      integer nx,ny,i,j,nxa,nyc,l
      double precision u(nx,ny)
      double precision dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,
     +dlx4,dly4,dlz4
      common/pde3com/dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,
     +               dlx4,dly4,dlz4
      double precision ux3,ux4,uy3,uy4
      l=ny
c----+|----------------------------------------------------------------|
c    x partial derivatives
c----+|----------------------------------------------------------------|
      call p3de1(nx,u(1,j),i,ux3,ux4,nxa)
c----+|----------------------------------------------------------------|
c    y partial derivatives
c----+|----------------------------------------------------------------|
      if (nyc.ne.0) then
c----+|----------------------------------------------------------------|
c    not periodic in y
c----+|----------------------------------------------------------------|
      if (j.gt.2 .and. j.lt.ny-1) then
      uy3 = (-u(i,j-2)+2.d0*u(i,j-1)-2.d0*u(i,j+1)+u(i,j+2))/tdly3
      uy4 = (u(i,j-2)-4.d0*u(i,j-1)+6.d0*u(i,j)-4.d0*u(i,j+1)+
     +      u(i,j+2))/ dly4
      else if (j.eq.1) then
      uy3 = (-5.d0*u(i,1)+18.d0*u(i,2)-24.d0*u(i,3)+14.d0*u(i,4)-
     +        3.d0*u(i,5))/tdly3
      uy4 = (3.d0*u(i,1)-14.d0*u(i,2)+26.d0*u(i,3)-24.d0*u(i,4)+
     +       11.d0*u(i,5)-2.d0*u(i,6))/dly4
      else if (j.eq.2) then
      uy3 = (-3.d0*u(i,1)+10.d0*u(i,2)-12.d0*u(i,3)+6.d0*u(i,4)-
     +      u(i,5))/ tdly3
      uy4 = (2.d0*u(i,1)-9.d0*u(i,2)+16.d0*u(i,3)-14.d0*u(i,4)+
     +       6.d0*u(i,5)-u(i,6))/dly4
      else if (j.eq.ny-1) then
      uy3 = (u(i,l-4)-6.d0*u(i,l-3)+12.d0*u(i,l-2)-10.d0*u(i,l-1)+
     +       3.d0*u(i,l))/tdly3
      uy4 = (-u(i,l-5)+6.d0*u(i,l-4)-14.d0*u(i,l-3)+16.d0*u(i,l-2)-
     +       9.d0*u(i,l-1)+2.d0*u(i,l))/dly4
      else if (j.eq.ny) then
      uy3 = (3.d0*u(i,l-4)-14.d0*u(i,l-3)+24.d0*u(i,l-2)-18.d0*u(i,l-1)+
     +      5.d0*u(i,l))/tdly3
      uy4 = (-2.d0*u(i,l-5)+11.d0*u(i,l-4)-24.d0*u(i,l-3)+
     +      26.d0*u(i,l-2)-14.d0*u(i,l-1)+3.d0*u(i,l))/dly4
      end if
      else
c----+|----------------------------------------------------------------|
c    periodic in y
c----+|----------------------------------------------------------------|
      if (j.gt.2 .and. j.lt.ny-1) then
      uy3 = (-u(i,j-2)+2.d0*u(i,j-1)-2.d0*u(i,j+1)+u(i,j+2))/tdly3
      uy4 = (u(i,j-2)-4.d0*u(i,j-1)+6.d0*u(i,j)-4.d0*u(i,j+1)+
     +      u(i,j+2))/dly4
      else if (j.eq.1) then
      uy3 = (-u(i,l-2)+2.d0*u(i,l-1)-2.d0*u(i,2)+u(i,3))/tdly3
      uy4 = (u(i,l-2)-4.d0*u(i,l-1)+6.d0*u(i,1)-4.d0*u(i,2)+
     +      u(i,3))/dly4 
      else if (j.eq.2) then
      uy3 = (-u(i,l-1)+2.d0*u(i,1)-2.d0*u(i,3)+u(i,4))/(tdly3)
      uy4 = (u(i,l-1)-4.d0*u(i,1)+6.d0*u(i,2)-4.d0*u(i,3)+u(i,4))/dly4
      else if (j.eq.ny-1) then
      uy3 = (-u(i,l-3)+2.d0*u(i,l-2)-2.d0*u(i,1)+u(i,2))/tdly3
      uy4 = (u(i,l-3)-4.d0*u(i,l-2)+6.d0*u(i,l-1)-4.d0*u(i,1)+u(i,2))/
     +        dly4
      else if (j.eq.ny) then
      uy3 = (-u(i,l-2)+2.d0*u(i,l-1)-2.d0*u(i,2)+u(i,3))/tdly3
      uy4 = (u(i,l-2)-4.d0*u(i,l-1)+6.d0*u(i,l)-4.d0*u(i,2)
     +      +u(i,3))/dly4
      end if
      end if
      return
      end subroutine p3de2

      subroutine p3de1(nx,u,i,ux3,ux4,nxa)
c----+|----------------------------------------------------------------|
c    third and fourth derivatives in x
c----+|----------------------------------------------------------------|
      implicit none
      integer nx,i,nxa,k
      double precision u(nx)
      double precision dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,
     +dlx4,dly4,dlz4
      common/pde3com/dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,
     +               dlx4,dly4,dlz4
      double precision ux3,ux4
      k = nx
      if (nxa.ne.0) then
c----+|----------------------------------------------------------------|
c    nonperiodic in x
c----+|----------------------------------------------------------------|
      if(i.gt.2 .and. i.lt.nx-1) then
      ux3 = (-u(i-2)+2.d0*u(i-1)-2.d0*u(i+1)+u(i+2))/tdlx3
      ux4 = (u(i-2)-4.d0*u(i-1)+6.d0*u(i)-4.d0*u(i+1)+u(i+2))/dlx4
      else if (i.eq.1) then
      ux3 = (-5.d0*u(1)+18.d0*u(2)-24.d0*u(3)+14.d0*u(4)
     +      -3.d0*u(5))/tdlx3
      ux4 = (3.d0*u(1)-14.d0*u(2)+26.d0*u(3)-24.d0*u(4)+11.d0*u(5)
     +      -2.d0*u(6))/dlx4
      else if (i.eq.2) then
      ux3 = (-3.d0*u(1)+10.d0*u(2)-12.d0*u(3)+6.d0*u(4)-u(5))/tdlx3
      ux4 = (2.d0*u(1)-9.d0*u(2)+16.d0*u(3)-14.d0*u(4)+6.d0*u(5)-
     +      u(6))/dlx4
      else if (i.eq.nx-1) then
      ux3 = (u(k-4)-6.d0*u(k-3)+12.d0*u(k-2)-10.d0*u(k-1)
     +      +3.d0*u(k))/tdlx3
      ux4 = (-u(k-5)+6.d0*u(k-4)-14.d0*u(k-3)+16.d0*u(k-2)-9.d0*u(k-1)+
     +      2.d0*u(k))/dlx4
      else if (i.eq.nx) then
      ux3 = (3.d0*u(k-4)-14.d0*u(k-3)+24.d0*u(k-2)-18.d0*u(k-1)+
     +       5.d0*u(k))/tdlx3
      ux4 = (-2.d0*u(k-5)+11.d0*u(k-4)-24.d0*u(k-3)+26.d0*u(k-2)-
     +       14.d0*u(k-1)+3.d0*u(k))/dlx4
      end if
      else
c----+|----------------------------------------------------------------|
c    periodic in x
c----+|----------------------------------------------------------------|
      if(i.gt.2 .and. i.lt.nx-1) then
      ux3 = (-u(i-2)+2.d0*u(i-1)-2.d0*u(i+1)+u(i+2))/tdlx3
      ux4 = (u(i-2)-4.d0*u(i-1)+6.d0*u(i)-4.d0*u(i+1)+u(i+2))/dlx4
      else if (i.eq.1) then
      ux3 = (-u(k-2)+2.d0*u(k-1)-2.d0*u(2)+u(3))/tdlx3
      ux4 = (u(k-2)-4.d0*u(k-1)+6.d0*u(1)-4.d0*u(2)+u(3))/dlx4
      else if (i.eq.2) then
      ux3 = (-u(k-1)+2.d0*u(1)-2.d0*u(3)+u(4))/(tdlx3)
      ux4 = (u(k-1)-4.d0*u(1)+6.d0*u(2)-4.d0*u(3)+u(4))/dlx4
      else if (i.eq.nx-1) then
      ux3 = (-u(k-3)+2.d0*u(k-2)-2.d0*u(1)+u(2))/tdlx3
      ux4 = (u(k-3)-4.d0*u(k-2)+6.d0*u(k-1)-4.d0*u(1)+u(2))/dlx4
      else if (i.eq.nx) then
      ux3 = (-u(k-2)+2.d0*u(k-1)-2.d0*u(2)+u(3))/tdlx3
      ux4 = (u(k-2)-4.d0*u(k-1)+6.d0*u(k)-4.d0*u(2)+u(3))/dlx4
      end if
      end if
      return
      end subroutine p3de1

c----+|----------------------------------------------------------------|
c----+|----------------------------------------------------------------|
c    factri and factrip are:
c    subroutines called by any real mudpack solver which uses line
c    relaxation(s) within multigrid iteration.  these subroutines do
c    a vectorized factorization of m simultaneous tridiagonal systems
c    of order n arising from nonperiodic or periodic discretizations
c----+|----------------------------------------------------------------|
c----+|----------------------------------------------------------------|
c    factor the m simultaneous tridiagonal systems of order n
c----+|----------------------------------------------------------------|
      subroutine factri(m,n,a,b,c)
      implicit none
      integer m,n,i,j
      double precision a(n,m),b(n,m),c(n,m)
      do i=2,n
      do j=1,m
        a(i-1,j) = a(i-1,j)/b(i-1,j)
        b(i,j) = b(i,j)-a(i-1,j)*c(i-1,j)
       end do
      end do
      return
      end subroutine factri

c----+|----------------------------------------------------------------|
c    factor the m simultaneous "tridiagonal" systems of order n
c    from discretized periodic system (leave out periodic n point)
c    (so sweeps below only go from i=1,2,...,n-1) n > 3 is necessary
c----+|----------------------------------------------------------------|
      subroutine factrp(m,n,a,b,c,d,e,ssm)
      implicit none
      integer m,n,i,j
      double precision a(n,m),b(n,m),c(n,m),d(n,m),e(n,m),ssm(m)
      do j=1,m
        d(1,j) = a(1,j)
      end do
      do i=2,n-2
        do j=1,m
          a(i,j) = a(i,j)/b(i-1,j)
          b(i,j) = b(i,j)-a(i,j)*c(i-1,j)
          d(i,j) = -a(i,j)*d(i-1,j)
        end do
      end do
c----+|----------------------------------------------------------------|
c    correct computation of last d element
c----+|----------------------------------------------------------------|
      do j=1,m
        d(n-2,j) = c(n-2,j)+d(n-2,j)
      end do
      do j=1,m
        e(1,j) = c(n-1,j)/b(1,j)
      end do
      do i=2,n-3
        do j=1,m
          e(i,j) = -e(i-1,j)*c(i-1,j)/b(i,j)
        end do
      end do
      do j=1,m
        e(n-2,j) = (a(n-1,j)-e(n-3,j)*c(n-3,j))/b(n-2,j)
      end do
c----+|----------------------------------------------------------------|
c    compute  inner product (e,d) for each j in ssm(j)
c----+|----------------------------------------------------------------|
      do j=1,m
        ssm(j) = 0.
      end do
      do i=1,n-2
        do j=1,m
          ssm(j) = ssm(j)+e(i,j)*d(i,j)
        end do
      end do
c----+|----------------------------------------------------------------|
c     set last diagonal element
c----+|----------------------------------------------------------------|
      do j=1,m
        b(n-1,j) = b(n-1,j)-ssm(j)
      end do
      return
      end subroutine factrp

c----+|----------------------------------------------------------------|
c    transpose n by n real matrix
c----+|----------------------------------------------------------------|
      subroutine transp(n,amat)
      implicit none
      integer n,i,j
      double precision amat(n,n),temp
      do i=1,n-1
        do j=i+1,n
          temp = amat(i,j)
          amat(i,j) = amat(j,i)
          amat(j,i) = temp
        end do
      end do
      return
      end subroutine transp

c----+|----------------------------------------------------------------|
c     ????????????
c----+|----------------------------------------------------------------|
      subroutine sgfa (a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info                                                
      double precision a(lda,1)                                                             
      double precision t                                                                    
      integer isfmax,j,k,kp1,l,nm1

      info = 0                                                                  
      nm1 = n - 1                                                               
      if (nm1 .ge. 1) then  
        do k = 1, nm1                                                          
          kp1 = k + 1                                                            
          l = isfmax(n-k+1,a(k,k),1) + k - 1
          ipvt(k) = l                                                            
          if (a(l,k) /= 0.d0) then                                        
            if (l /= k) then
              t = a(l,k)                                                       
              a(l,k) = a(k,k)                                                  
              a(k,k) = t                                                       
            end if                                                     
            t = -1.d0/a(k,k)                                                   
            call sscl(n-k,t,a(k+1,k),1)
            do j = kp1, n                                                    
              t = a(l,j)                                                       
              if (l /= k) then                                           
                a(l,j) = a(k,j)                                               
                a(k,j) = t                                                    
              end if                                                         
              call sxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
            end do                                                            
          else                                                               
            info = k                                                            
          end if                                                               
        end do                                                                  
      end if                                                                  
      ipvt(n) = n                                                               
      if (a(n,n) .eq. 0.d0) info = n                                           
      return                                                                    
      end subroutine sgfa                                                                       
                                                                                
c----+|----------------------------------------------------------------|
c----+|----------------------------------------------------------------|
      subroutine sgsl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job                                                 
      double precision a(lda,1),b(1)                                                        
      double precision sdt,t
      integer k,kb,l,nm1                                                        

      nm1 = n - 1

      if (job .eq. 0) then                                                  
        if (nm1 .ge. 1) then                                               
          do k = 1, nm1                                                       
            l = ipvt(k)                                                         
            t = b(l)                                                            
            if (l .ne. k) then                                              
              b(l) = b(k)                                                      
              b(k) = t                                                         
            end if                                                            
            call sxpy(n-k,t,a(k+1,k),1,b(k+1),1)
          end do                                                               
        end if                                                              
        do kb = 1, n                                                        
          k = n + 1 - kb                                                      
          b(k) = b(k)/a(k,k)                                                  
          t = -b(k)                                                           
          call sxpy(k-1,t,a(1,k),1,b(1),1)
        end do
        return
      end if                                                                 
      do k = 1, n                                                         
        t = sdt(k-1,a(1,k),1,b(1),1)
        b(k) = (b(k) - t)/a(k,k)                                            
      end do
      if (nm1 .ge. 1) then
        do kb = 1, nm1                                                      
          k = n - kb                                                          
          b(k) = b(k) + sdt(n-k,a(k+1,k),1,b(k+1),1)
          l = ipvt(k)                                                         
          if (l .ne. k) then                                              
            t = b(l)                                                         
            b(l) = b(k)                                                      
            b(k) = t                                                         
          end if                                                          
        end do                                                               
      end if
      return                                                                    
      end subroutine sgsl                                                                     
                                                                                
c----+|----------------------------------------------------------------|
      double precision function sdt(n,sx,incx,sy,incy)
      double precision sx(1),sy(1),stemp                                                    
      integer i,incx,incy,ix,iy,m,mp1,n                                         
      
      stemp = 0.d0                                                             
      sdt = 0.d0
      
      if(n.le.0)return                                                          
      
      if (incx.ne.1 .or. incy.ne.1) then                                       
        ix = 1                                                                    
        iy = 1                                                                    
        if(incx.lt.0)ix = (-n+1)*incx + 1                                         
        if(incy.lt.0)iy = (-n+1)*incy + 1                                         
        do i = 1,n                                                             
          stemp = stemp + sx(ix)*sy(iy)                                           
          ix = ix + incx                                                          
          iy = iy + incy                                                          
        end do
        sdt = stemp
        return
      end if

      m = mod(n,5) 

      if( m .ne. 0 ) then                                                   
        do i = 1,m                                                             
          stemp = stemp + sx(i)*sy(i)                                             
        end do
        if (n .lt. 5) then 
          sdt = stemp
          return
        end if    
      end if
      
      mp1 = m + 1                                                               
      do i = mp1,n,5                                                         
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +                     
     +  sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)        
      end do
      sdt = stemp
      return                                                                    
      end function sdt
c----+|----------------------------------------------------------------|
                                                                                
c----+|----------------------------------------------------------------|
      integer function isfmax(n,sx,incx)
      double precision sx(1),smax                                                           
      integer i,incx,ix,n 

      isfmax = 0
      if (n.lt.1) return                                                     
      isfmax = 1
      if (n.eq.1) return 

      if (incx.ne.1) then                                                     
        ix = 1                                                                    
        smax = abs(sx(1))                                                         
        ix = ix + incx                                                            
        do i = 2,n                                                             
          if (abs(sx(ix)).gt.smax) then                                        
            isfmax = i
            smax = abs(sx(ix))
          end if            
          ix = ix + incx                                                         
        end do
        return
      end if

      smax = abs(sx(1))                                                         
      do i = 2,n                                                             
         if(abs(sx(i)).gt.smax) then                                        
           isfmax = i
           smax = abs(sx(i))                                                      
         end if
      end do         
      return                                                                    
      end function isfmax                                                                       

c----+|----------------------------------------------------------------|
      subroutine sxpy(n,sa,sx,incx,sy,incy)
      double precision sx(1),sy(1),sa                                                       
      integer i,incx,incy,ix,iy,m,mp1,n                                         
      
      if(n.le.0)return                                                          
      if (sa .eq. 0.d0) return                                                   
      
      if(incx.ne.1 .or. incy.ne.1) then                                      
        ix = 1                                                                    
        iy = 1                                                                    
        if(incx.lt.0)ix = (-n+1)*incx + 1                                         
        if(incy.lt.0)iy = (-n+1)*incy + 1                                         
        do i = 1,n                                                             
          sy(iy) = sy(iy) + sa*sx(ix)                                             
          ix = ix + incx                                                          
          iy = iy + incy                                                          
        end do                                                                  
        return                                                                    
      end if

      m = mod(n,4)                                                              
      if( m .ne. 0 ) then                                                   
        do i = 1,m                                                             
          sy(i) = sy(i) + sa*sx(i)                                                
        end do
        if( n .lt. 4 ) return 
      end if

      mp1 = m + 1

      do i = mp1,n,4                                                         
        sy(i) = sy(i) + sa*sx(i)                                                
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)                                    
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)                                    
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)                                    
      end do                                                                  
      return                                                                    
      end subroutine sxpy                                                                      

c----+|----------------------------------------------------------------|
      subroutine sscl(n,sa,sx,incx)
      double precision sa,sx(1)                                                             
      integer i,incx,m,mp1,n,nincx                                              

      if(n.le.0)return

      if(incx.ne.1) then                                                     
        nincx = n*incx                                                            
        do i = 1,nincx,incx                                                    
          sx(i) = sa*sx(i)                                                        
        end do                                                                 
        return 
      end if

      m = mod(n,5)                                                              
      if( m .ne. 0 ) then                                                   
        do i = 1,m                                                             
          sx(i) = sa*sx(i)                                                        
        end do
        if( n .lt. 5 ) return                                                     
      end if
      mp1 = m + 1                                                               
      do i = mp1,n,5                                                         
        sx(i) = sa*sx(i)                                                        
        sx(i + 1) = sa*sx(i + 1)                                                
        sx(i + 2) = sa*sx(i + 2)                                                
        sx(i + 3) = sa*sx(i + 3)                                                
        sx(i + 4) = sa*sx(i + 4)                                                
      end do                                                                  
      return                                                                    
      end subroutine sscl                                                                      

