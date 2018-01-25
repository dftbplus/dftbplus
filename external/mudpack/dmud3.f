c
c     file mud3.f
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1999 by UCAR                 .
c  .                                                             .
c  .       University Corporation for Atmospheric Research       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                      MUDPACK  version 5.0                   .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c ... author and specialist
c
c          John C. Adams (National Center for Atmospheric Research)
c          email: johnad@ucar.edu, phone: 303-497-1213

c ... For MUDPACK 5.0 information, visit the website:
c     (http://www.scd.ucar.edu/css/software/mudpack)
c
c ... purpose (see mud3.d for details)
c
c     mud3 attempts to produce a second order finite difference
c     approximation to the three dimensional nonseparable elliptic
c     partial differential equation of the form:
c
c       cxx(x,y,z)*pxx + cyy(x,y,z)*pyy + czz(x,y,z)*pzz +
c
c       cx(x,y,z)*px + cy(x,y,z)*py + cz(x,y,z)*pz +
c
c       ce(x,y,z)*p(x,y,z) = r(x,y,z)
c
c ... documentation and test files
c
c     see the documentation file "mud3.d" for a complete discussion
c     of how to use subroutine mud3.  file "tmud3.f" is a test/driver
c     sample program illustrating use of mud3
c
c ... required MUDPACK files
c
c     mudcom.f, mud3ln.f, mud3pn.f
c
c----+|-----------------------------------------------------------------------|
      subroutine mud3(iparm,fparm,work,coef,bndyc,rhs,phi,mgopt,ierror)
      implicit none
      integer iparm(23),mgopt(4),ierror
      double precision work(*),phi(*),rhs(*),fparm(8)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      double precision xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/fmud3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,ktzbgn,nxk,nyk,nzk,ngrid,
     +        klevel,kcur,kps
      common/mud3c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +ktzbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
      integer sxy,sxz,syz,kpspace,kcspace,int
      integer m,lxy,lxz,lyz,isx,jsy,ksz,iisx,jjsy,kksz
      integer nx,ny,nz,itx,ity,itz,k,kb,kkb,kk,ip,ic
      external coef,bndyc
      data int / 0 /
      save int
      ierror = 1
      intl = iparm(1)    ! set and check intl on ALL calls
      if (intl*(intl-1).ne.0) return
      if (int.eq.0) then
      int = 1
      if (intl.ne.0) return  ! very first call is not intl=0
      end if
      ierror = 0
c
c----+|-----------------------------------------------------------------------|
c     set input parameters from iparm,fparm internally
c----+|-----------------------------------------------------------------------|
c
      intl = iparm(1)
      nxa = iparm(2)
      nxb = iparm(3)
      nyc = iparm(4)
      nyd = iparm(5)
      nze = iparm(6)
      nzf = iparm(7)
c
c----+|-----------------------------------------------------------------------|
c     set grid size params
c----+|-----------------------------------------------------------------------|
c
      ixp = iparm(8)
      jyq = iparm(9)
      kzr = iparm(10)
      iex = iparm(11)
      jey = iparm(12)
      kez = iparm(13)
c
c----+|-----------------------------------------------------------------------|
c     set number of subgrids for mg cycling
c----+|-----------------------------------------------------------------------|
c
      ngrid = max0(iex,jey,kez)
      nfx = iparm(14)
      nfy = iparm(15)
      nfz = iparm(16)

      iguess = iparm(17)
      maxcy = iparm(18)
      method = iparm(19)
      meth2 = iparm(20)
      nwork = iparm(21)
c
c----+|-----------------------------------------------------------------------|
c     set floating point params
c----+|-----------------------------------------------------------------------|
c
      xa = fparm(1)
      xb = fparm(2)
      yc = fparm(3)
      yd = fparm(4)
      ze = fparm(5)
      zf = fparm(6)
      tolmax = fparm(7)
c
c----+|-----------------------------------------------------------------------|
c     set multigrid option parameters or defaults 
c----+|-----------------------------------------------------------------------|
      if (mgopt(1) .eq. 0) then
        kcycle = 2
        iprer = 2
        ipost = 1
        intpol = 3
      else
        kcycle = mgopt(1)
        iprer = mgopt(2)
        ipost = mgopt(3)
        intpol = mgopt(4)
      end if

      if (intl .eq. 0) then  ! intialization call
c----+|-----------------------------------------------------------------------|
c       check input arguments
c----+|-----------------------------------------------------------------------|
        ierror = 2   ! check boundary condition flags
        if (max0(nxa,nxb,nyc,nyd,nze,nzf).gt.2) return
        if (min0(nxa,nxb,nyc,nyd,nze,nzf).lt.0) return
        if (nxa.eq.0.and.nxb.ne.0) return
        if (nxa.ne.0.and.nxb.eq.0) return
        if (nyc.eq.0.and.nyd.ne.0) return
        if (nyc.ne.0.and.nyd.eq.0) return
        if (nze.eq.0.and.nzf.ne.0) return
        if (nze.ne.0.and.nzf.eq.0) return
        ierror = 3   ! check grid sizes
        if (ixp.lt.2) return
        if (jyq.lt.2) return
        if (kzr.lt.2) return
        ierror = 4
        ngrid = max0(iex,jey,kez)
        if (iex.lt.1) return
        if (jey.lt.1) return
        if (kez.lt.1) return
        if (ngrid.gt.50) return
        ierror = 5
        if (nfx.ne.ixp*2**(iex-1)+1) return
        if (nfy.ne.jyq*2**(jey-1)+1) return
        if (nfz.ne.kzr*2**(kez-1)+1) return
        ierror = 6
        if (iguess*(iguess-1).ne.0) return
        ierror = 7
        if (maxcy.lt.1) return
        ierror = 8
        if (method.lt.0 .or. method.gt.10) return
        if (meth2.lt.0 .or. meth2.gt.3) return
        ierror = 9
c
c----+|-----------------------------------------------------------------------|
c       compute required work space length
c----+|-----------------------------------------------------------------------|
        m = method
        isx = 0
        if ((m-1)*(m-4)*(m-5)*(m-7).eq.0) then
          isx = 3
          if (nxa.eq.0) then
            isx = 5
          end if
        end if
        jsy = 0
        if ((m-2)*(m-4)*(m-6)*(m-7).eq.0) then
          jsy = 3
          if (nyc.eq.0) then
            jsy = 5
          end if
        end if
        ksz = 0
        if ((m-3)*(m-5)*(m-6)*(m-7).eq.0) then
          ksz = 3
          if (nze.eq.0) then
            ksz = 5
          end if
        end if
c
c----+|-----------------------------------------------------------------------|
c     set scales for planar relaxation
c----+|-----------------------------------------------------------------------|
c
        iisx = 0
        jjsy = 0
        kksz = 0
        lxy = 0
        lxz = 0
        lyz = 0
        nx = nfx
        ny = nfy
        nz = nfz
        if (method.eq.8) then
          lxy = 1
          if (meth2.eq.1.or.meth2.eq.3) then
            iisx = 3
            if (nxa.eq.0) iisx = 5
          end if
          if (meth2.eq.2.or.meth2.eq.3) then
            jjsy = 3
            if (nyc.eq.0) jjsy = 5
          end if
        end if
        if (method.eq.9) then
          lxz = 1
          if (meth2.eq.1.or.meth2.eq.3) then
            iisx = 3
            if (nxa.eq.0) iisx = 5
          end if
          if (meth2.eq.2.or.meth2.eq.3) then
            kksz = 3
            if (nze.eq.0) kksz = 5
            end if
        end if
        if (method.eq.10) then
          lyz = 1
          if (meth2.eq.1.or.meth2.eq.3) then
            jjsy = 3
            if (nyc.eq.0) jjsy = 5
          end if
          if (meth2.eq.2.or.meth2.eq.3) then
            kksz = 3
            if (nze.eq.0) kksz = 5
          end if
        end if
        sxy = 0
        sxz = 0
        syz = 0
c
c----+|-----------------------------------------------------------------------|
c       set subgrid sizes
c----+|-----------------------------------------------------------------------|
        do k=1,ngrid
          nxk(k) = ixp*2**(max0(k+iex-ngrid,1)-1)+1
          nyk(k) = jyq*2**(max0(k+jey-ngrid,1)-1)+1
          nzk(k) = kzr*2**(max0(k+kez-ngrid,1)-1)+1
        end do
        kps = 1
        do kb=1,ngrid
          k = ngrid-kb+1
          kpspace = 0
          kcspace = 0
c----+|-----------------------------------------------------------------------|
c         set spacers for planar relaxation
c----+|-----------------------------------------------------------------------|
          if (method.gt.7) then
            do kkb=2,k
              kk = k-kkb+1
              nx = nxk(kk)
              ny = nyk(kk)
              nz = nzk(kk)
              if (method.eq.8) nz = nzk(k)
              if (method.eq.9) ny = nyk(k)
              if (method.eq.10) nx = nxk(k)
              kpspace = kpspace + (nx+2)*(ny+2)*(nz+2)
              kcspace = kcspace + 8*nx*ny*nz
            end do
          end if
          nx = nxk(k)
          ny = nyk(k)
          nz = nzk(k)
c----+|-----------------------------------------------------------------------|
c         set pointers
c----+|-----------------------------------------------------------------------|
          kpbgn(k) = kps
          kcbgn(k) = kpbgn(k)+(nx+2)*(ny+2)*(nz+2)+kpspace
          ktxbgn(k) = kcbgn(k) + 8*nx*ny*nz + kcspace
          ktybgn(k) = ktxbgn(k) + isx*nx*ny*nz
          ktzbgn(k) = ktybgn(k) + jsy*nx*ny*nz
          kps = ktzbgn(k) + ksz*nx*ny*nz
c----+|-----------------------------------------------------------------------|
c         sum space in case planar relaxation in xy or xz or yz
c----+|-----------------------------------------------------------------------|
          sxy = sxy + (6+iisx+jjsy)*nx*ny + (nx+2)*(ny+2)
          sxz = sxz + (6+iisx+kksz)*nx*nz + (nx+2)*(nz+2)
          syz = syz + (6+jjsy+kksz)*ny*nz + (ny+2)*(nz+2)
        end do
c----+|-----------------------------------------------------------------------|
c       set and check minimal work space
c----+|-----------------------------------------------------------------------|
        nx = nxk(ngrid)
        ny = nyk(ngrid)
        nz = nzk(ngrid)
        iparm(22)=kps+max0((nx+2)*(ny+2)*(nz+2),lxy*sxy,lxz*sxz,lyz*syz)
        lwork = iparm(22)
        if (lwork .gt. nwork) return
        ierror = 10   ! check solution region
        if (xb.le.xa .or. yd.le.yc .or. zf.le.ze) return
        ierror = 11
        if (tolmax .lt. 0.d0) return
        ierror = 12   ! multigrid parameters
        if (kcycle.lt.0) return
        if (min0(iprer,ipost).lt.1) return
        if ((intpol-1)*(intpol-3).ne.0) return
        if (max0(kcycle,iprer,ipost).gt.2) then
          ierror = -5   ! inefficient multigrid cycling
        end if
        if (ierror .gt. 0) ierror = 0   ! no fatal errors
c
c----+|-----------------------------------------------------------------------|
c       discretize pde at each grid level
c----+|-----------------------------------------------------------------------|
c
        do kb=1,ngrid
          k = ngrid-kb+1
          nx = nxk(k)
          ny = nyk(k)
          nz = nzk(k)
          ip = kpbgn(k)
          ic = kcbgn(k)
          itx = ktxbgn(k)
          ity = ktybgn(k)
          itz = ktzbgn(k)
          call dismd3(nx,ny,nz,work(ic),work(itx),
     +        work(ity),work(itz),bndyc,coef,work,ierror)

c----+|-----------------------------------------------------------------------|
c         discretize for planar coarsening
c----+|-----------------------------------------------------------------------|
          if (method.gt.7) then
            do kkb=2,k
              kk = k-kkb+1
              ip = ip+(nx+2)*(ny+2)*(nz+2)
              ic = ic + 8*nx*ny*nz
              nx = nxk(kk)
              ny = nyk(kk)
              nz = nzk(kk)
              if (method.eq.8) nz = nzk(k)
              if (method.eq.9) ny = nyk(k)
              if (method.eq.10) nx = nxk(k)
              call dismd3(nx,ny,nz,work(ic),work(itx),work(ity),
     +             work(itz),bndyc,coef,work,ierror)
            end do
          end if
        end do
        return
      end if   ! end of intl=0 initialization call block
      nx = nfx
      ny = nfy
      nz = nfz

      call mud31(nx,ny,nz,rhs,phi,coef,bndyc,work)

      iparm(23) = itero
c
c----+|-----------------------------------------------------------------------|
c     set final computed maximum relative difference
c     or flag convergence failure
c----+|-----------------------------------------------------------------------|
      if (ierror .le. 0) then
        if (tolmax.gt.0.d0) then
          fparm(8) = relmax
          if (relmax.gt.tolmax .and. ierror.eq.0) then
            ierror = -1
          end if
        end if
      end if
      return
      end subroutine mud3


c----+|-----------------------------------------------------------------------|
c     Solver
c----+|-----------------------------------------------------------------------|
      subroutine mud31(nx,ny,nz,rhsf,phif,coef,bndyc,wk)
      implicit none
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      double precision xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,ktzbgn,nxk,nyk,nzk,ngrid,
     +        klevel,kcur,kps
      integer nx,ny,nz,ip,ic,ir,icc,irc,ncx,ncy,ncz,ipc
      integer i,j,k,kb,jk,kk,kkb,ijk,iter
      double precision phif(nx,ny,nz),rhsf(nx,ny,nz),wk(*),phmax
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/fmud3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/mud3c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +ktzbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
      external coef,bndyc
      nx = nxk(ngrid)
      ny = nyk(ngrid)
      nz = nzk(ngrid)
      ip = kpbgn(ngrid)
      ic = kcbgn(ngrid)
      ir = ic+7*nx*ny*nz
c
c----+|-----------------------------------------------------------------------|
c     set phif,rhsf in wk
c----+|-----------------------------------------------------------------------|
      call swk3(nx,ny,nz,phif,rhsf,wk(ip),wk(ir))
c
c----+|-----------------------------------------------------------------------|
c     no initial guess at finest grid level
c----+|-----------------------------------------------------------------------|
      if (iguess.eq.0) then
        do kb=2,ngrid
          k = ngrid-kb+1
          nx = nxk(k+1)
          ny = nyk(k+1)
          nz = nzk(k+1)
          ip = kpbgn(k+1)
          ic = kcbgn(k+1)
          ir = kcbgn(k+1)+7*nx*ny*nz
          ncx = nxk(k)
          ncy = nyk(k)
          ncz = nzk(k)
          ipc = kpbgn(k)
          icc = kcbgn(k)
          irc = icc+7*ncx*ncy*ncz
c
c----+|-----------------------------------------------------------------------|
c         transfers down to all grid levels
c----+|-----------------------------------------------------------------------|
          call trsfc3(nx,ny,nz,wk(ip),wk(ir),ncx,ncy,ncz,
     +                wk(ipc),wk(irc))
c
c----+|-----------------------------------------------------------------------|
c         transfers down for planar coarsening
c----+|-----------------------------------------------------------------------|
          if (method.gt.7) then
            do kkb=1,k-1
              kk = k-kkb+1
              ipc = ip + (nx+2)*(ny+2)*(nz+2)
              icc = ic + 8*nx*ny*nz
              ncx = nxk(kk)
              ncy = nyk(kk)
              ncz = nzk(kk)
              if (method.eq.8) then
                ncz = nzk(k+1)
              else if (method.eq.9) then
                ncy = nyk(k+1)
              else
                ncx = nxk(k+1)
              end if
              irc = icc+7*ncx*ncy*ncz
              call trsfc3(nx,ny,nz,wk(ip),wk(ir),ncx,ncy,ncz,
     +                    wk(ipc),wk(irc))
              nx = ncx
              ny = ncy
              nz = ncz
              ip = ipc
              ir = irc
              ic = icc
            end do
          end if
        end do
c
c----+|-----------------------------------------------------------------------|
c       Adjust right hand side at all grid levels in case
c       rhs or specified b.c. in phi or gbdy changed
c----+|-----------------------------------------------------------------------|
        do kb=1,ngrid
          k = ngrid-kb+1
          nx = nxk(k)
          ny = nyk(k)
          nz = nzk(k)
          ip = kpbgn(k)
          ic = kcbgn(k)
          call adjmd3(nx,ny,nz,wk(ip),wk(ic),bndyc,coef)
c
c----+|-----------------------------------------------------------------------|
c         Adjust for planar grid levels if necessary
c----+|-----------------------------------------------------------------------|
          if (method.gt.7) then
            do kkb=2,k
              kk = k-kkb+1
              ip = ip+(nx+2)*(ny+2)*(nz+2)
              ic = ic + 8*nx*ny*nz
              nx = nxk(kk)
              ny = nyk(kk)
              nz = nzk(kk)
              if (method.eq.8) then
              nz = nzk(k)
              else if (method.eq.9) then
              ny = nyk(k)
              else
              nx = nxk(k)
              end if
              call adjmd3(nx,ny,nz,wk(ip),wk(ic),bndyc,coef)
            end do
          end if
        end do
c
c----+|-----------------------------------------------------------------------|
c       Execute one full multigrid cycle
c       Lift or prolong approximation from k to k+1
c----+|-----------------------------------------------------------------------|
        do k=1,ngrid-1
          kcur = k
          call kcymd3(wk)
          nx = nxk(k+1)
          ny = nyk(k+1)
          nz = nzk(k+1)
          ip = kpbgn(k+1)
          ipc = kpbgn(k)
          ncx = nxk(k)
          ncy = nyk(k)
          ncz = nzk(k)
          call prolon3(ncx,ncy,ncz,wk(ipc),nx,ny,nz,wk(ip),
     +                 nxa,nxb,nyc,nyd,nze,nzf,intpol)
        end do

      else
c----+|-----------------------------------------------------------------------|
c       We have an initial guess at finest grid level
c       Adjust rhs at finest grid level
c----+|-----------------------------------------------------------------------|
        nx = nxk(ngrid)
        ny = nyk(ngrid)
        nz = nzk(ngrid)
        ip = kpbgn(ngrid)
        ic = kcbgn(ngrid)
        call adjmd3(nx,ny,nz,wk(ip),wk(ic),bndyc,coef)
      end if
c
c----+|-----------------------------------------------------------------------|
c     Execute maxcy more multigrid k cycles from finest level
c     check for convergence on the maximum relative difference
c----+|-----------------------------------------------------------------------|
      kcur = ngrid
      do iter=1,maxcy
        itero = iter
        call kcymd3(wk)
        if (tolmax.gt.0.d0) then
          relmax = 0.d0
          phmax = 0.d0
          do k=1,nfz
            kk = k*(nfx+2)*(nfy+2)
            do j=1,nfy
              jk = kk+j*(nfx+2)
              do i=1,nfx
                ijk = jk+i+1
                phmax = dmax1(phmax,abs(wk(ijk)))
                relmax = dmax1(relmax,abs(wk(ijk)-phif(i,j,k)))
                phif(i,j,k) = wk(ijk)
              end do
            end do
          end do
          if (phmax.gt.0.d0) relmax = relmax/phmax
          if (relmax.le.tolmax) return
        end if
      end do
c
c----+|-----------------------------------------------------------------------|
c     Set final iterate in phif
c----+|-----------------------------------------------------------------------|
      do k=1,nfz
        kk = k*(nfx+2)*(nfy+2)
        do j=1,nfy
          jk = kk+j*(nfx+2)
          do i=1,nfx
            ijk = jk+i+1
            phif(i,j,k) = wk(ijk)
          end do
        end do
      end do
      return
      end subroutine mud31

c
c----+|-----------------------------------------------------------------------|
c     perform multigrid k-cycle at kcur level
c     kcycle = 1 corresponds to v cycles
c     kcycle = 2 corresponds to w cycles
c----+|-----------------------------------------------------------------------|
      subroutine kcymd3(wk)
      implicit none
      double precision wk(*)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      integer kpbgn,kcbgn,ktxbgn,ktybgn,ktzbgn,nxk,nyk,nzk,ngrid,
     +        klevel,kcur,kps
      common/mud3c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +ktzbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
      integer nx,ny,nz,ncx,ncy,ncz
      integer kount(50),ip,ic,ipc,irc,nrel,l
      klevel = kcur
c
c----+|-----------------------------------------------------------------------|
c     pre-relax at current finest grid level
c----+|-----------------------------------------------------------------------|
      do l = 1,iprer
        call relmd3(wk)
      end do
c
c----+|-----------------------------------------------------------------------|
c      Iff at coarsest level post-relax
c----+|-----------------------------------------------------------------------|
c     restrict residual to kcur-1
c----+|-----------------------------------------------------------------------|
      if (kcur .eq. 1) go to 2

      nx = nxk(klevel)
      ny = nyk(klevel)
      nz = nzk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      ipc = kpbgn(klevel-1)
      ncx = nxk(klevel-1)
      ncy = nyk(klevel-1)
      ncz = nzk(klevel-1)
      irc = kcbgn(klevel-1)+7*ncx*ncy*ncz
c----+|-----------------------------------------------------------------------|
c     Full weighting with residual restriction
c----+|-----------------------------------------------------------------------|
      call resmd3(nx,ny,nz,wk(ip),wk(ic),ncx,ncy,ncz,wk(ipc),wk(irc),
     +          wk(kps))
c
c----+|-----------------------------------------------------------------------|
c     Set counter for grid levels to zero
c----+|-----------------------------------------------------------------------|
      do l = 1,kcur
        kount(l) = 0
      end do
c
c----+|-----------------------------------------------------------------------|
c     Set new level and continue k-cycling
c----+|-----------------------------------------------------------------------|
      klevel = kcur-1
      nrel = iprer
c
c----+|-----------------------------------------------------------------------|
c     Kcycle control point
c----+|-----------------------------------------------------------------------|
    1 continue
c
c----+|-----------------------------------------------------------------------|
c      post-relax when kcur revisited
c----+|-----------------------------------------------------------------------|
      if (klevel .eq. kcur) go to 2
c
c----+|-----------------------------------------------------------------------|
c     Counts "hit" at current level
c----+|-----------------------------------------------------------------------|
      kount(klevel) = kount(klevel)+1
c
c----+|-----------------------------------------------------------------------|
c     Relax at current level
c----+|-----------------------------------------------------------------------|
      do l = 1,nrel
        call relmd3(wk)
      end do
      
c----+|-----------------------------------------------------------------------|
c     Kcycle(iprer,ipost) complete at klevel
c     Inject correction to finer grid
c----+|-----------------------------------------------------------------------|
      if (kount(klevel) .eq. kcycle+1) then
        nx = nxk(klevel+1)
        ny = nyk(klevel+1)
        nz = nzk(klevel+1)
        ip = kpbgn(klevel+1)
        ncx = nxk(klevel)
        ncy = nyk(klevel)
        ncz = nzk(klevel)
        ipc = kpbgn(klevel)
        call cor3(nx,ny,nz,wk(ip),ncx,ncy,ncz,wk(ipc),
     +            nxa,nxb,nyc,nyd,nze,nzf,intpol,wk(kps))
c
c----+|-----------------------------------------------------------------------|
c       Reset counter to zero at klevel
c----+|-----------------------------------------------------------------------|
        kount(klevel) = 0
c
c----+|-----------------------------------------------------------------------|
c       Ascend to next higher level and set to post-relax there
c----+|-----------------------------------------------------------------------|
        klevel = klevel+1
        nrel = ipost
        go to 1
      else
c----+|-----------------------------------------------------------------------|
c       Kcycle not complete so descend unless at coarsest
c----+|-----------------------------------------------------------------------|
        if (klevel .gt. 1) then
          nx = nxk(klevel)
          ny = nyk(klevel)
          nz = nzk(klevel)
          ip = kpbgn(klevel)
          ic = kcbgn(klevel)
          ncx = nxk(klevel-1)
          ncy = nyk(klevel-1)
          ncz = nzk(klevel-1)
          irc = kcbgn(klevel-1)+7*ncx*ncy*ncz
          ipc = kpbgn(klevel-1)
          call resmd3(nx,ny,nz,wk(ip),wk(ic),ncx,ncy,ncz,wk(ipc),
     +                 wk(irc),wk(kps))
c
c----+|-----------------------------------------------------------------------|
c         Pre-relax at next coarser level
c----+|-----------------------------------------------------------------------|
          klevel = klevel-1
          nrel = iprer
          go to 1
        else
c----+|-----------------------------------------------------------------------|
c         Post-relax at coarsest level (klevel=1)
c----+|-----------------------------------------------------------------------|
          do l = 1,ipost
            call relmd3(wk)
          end do
c----+|-----------------------------------------------------------------------|
c         Inject correction to grid level 2
c----+|-----------------------------------------------------------------------|
          ipc = kpbgn(1)
          ncx = nxk(1)
          ncy = nyk(1)
          ncz = nzk(1)
          ip = kpbgn(2)
          nx = nxk(2)
          ny = nyk(2)
          nz = nzk(2)
          call cor3(nx,ny,nz,wk(ip),ncx,ncy,ncz,wk(ipc),
     +                nxa,nxb,nyc,nyd,nze,nzf,intpol,wk(kps))
c
c----+|-----------------------------------------------------------------------|
c         Set to post-relax at level 2
c----+|-----------------------------------------------------------------------|
c
          nrel = ipost
          klevel = 2
          go to 1
        end if
      end if
c
    2 continue
c
c----+|-----------------------------------------------------------------------|
c     post-relax at kcur level
c----+|-----------------------------------------------------------------------|
      do l = 1,ipost
        call relmd3(wk)
      end do
      return
      end subroutine kcymd3

c----+|-----------------------------------------------------------------------|
c     compute fully weighted residual restriction in rhsc
c----+|-----------------------------------------------------------------------|
      subroutine resmd3(nx,ny,nz,phi,cof,ncx,ncy,ncz,phic,rhsc,resf)
      implicit none
      integer nx,ny,nz,ncx,ncy,ncz
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      double precision phi(0:nx+1,0:ny+1,0:nz+1)
      double precision phic(0:ncx+1,0:ncy+1,0:ncz+1)
      double precision rhsc(ncx,ncy,ncz),resf(nx,ny,nz),cof(nx,ny,nz,8)
      integer ic,jc,kc,i,j,k
c
c     initialize phic to zero
c
      do kc=0,ncz+1
        do jc=0,ncy+1
          do ic=0,ncx+1
            phic(ic,jc,kc) = 0.d0
          end do
        end do
      end do
c
c     compute fine grid residual
c
!$OMP PARALLEL DO PRIVATE(i,j,k), SHARED(resf,cof,phi,nx,ny,nz)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            resf(i,j,k) =  cof(i,j,k,8)-(
     +                     cof(i,j,k,1)*phi(i-1,j,k)+
     +                     cof(i,j,k,2)*phi(i+1,j,k)+
     +                     cof(i,j,k,3)*phi(i,j-1,k)+
     +                     cof(i,j,k,4)*phi(i,j+1,k)+
     +                     cof(i,j,k,5)*phi(i,j,k-1)+
     +                     cof(i,j,k,6)*phi(i,j,k+1)+
     +                     cof(i,j,k,7)*phi(i,j,k))
          end do
        end do
      end do
!$OMP END PARALLEL DO
c
c     restrict resf to interior coarse mesh in rhsc
c     using fully weighted residual restriction in 3-d
c
      call res3(nx,ny,nz,resf,ncx,ncy,ncz,rhsc,nxa,nxb,nyc,nyd,nze,nzf)
      return
      end subroutine resmd3

c
c----+|-----------------------------------------------------------------------|
c     discretize the 3-d elliptic pde
c----+|-----------------------------------------------------------------------|
      subroutine dismd3(nx,ny,nz,cof,tx,ty,tz,bndyc,coef,wk,ier)
      implicit none
      integer nx,ny,nz,ier
      double precision cof(nx,ny,nz,8)
      double precision tx(nx,ny,nz,*),ty(ny,nx,nz,*),tz(nz,nx,ny,*)
      double precision wk(*)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      double precision xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,ktzbgn,nxk,nyk,nzk,ngrid,
     +        klevel,kcur,kps
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/fmud3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/mud3c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +ktzbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
      double precision dlx,dly,dlz,dlx2,dly2,dlz2,dlxx,dlyy,dlzz,cmin,
     +              cemax,alfmax
      double precision cxx,cyy,czz,cx,cy,cz,ce,alfa,gbdy,x,y,z,c1,c2,c3,
     +              c4,c5,c6
      integer i,j,k,l,ist,ifn,jst,jfn,kst,kfn,kbdy
      integer nxny,nxnz,nynz,im1,jm1,km1
      external bndyc,coef
c
c----+|-----------------------------------------------------------------------|
c     set current grid increments
c----+|-----------------------------------------------------------------------|
      dlx = (xb-xa)/(nx-1)
      dlx2 = dlx+dlx
      dlxx = dlx*dlx
      dly = (yd-yc)/(ny-1)
      dly2 = dly+dly
      dlyy = dly*dly
      dlz = (zf-ze)/(nz-1)
      dlz2 = dlz+dlz
      dlzz = dlz*dlz
      cmin = 1.d0
      cemax = 0.d0
c
c----+|-----------------------------------------------------------------------|
c     set x,y,z subscript limits to bypass specified boundaries
c     when calling coef or bndyc
c----+|-----------------------------------------------------------------------|
      jst = 1
      jfn = ny
      ist = 1
      ifn = nx
      kst = 1
      kfn = nz
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.1) jst = 2
      if (nyd.eq.1) jfn = ny-1
      if (nze.eq.1) kst = 2
      if (nzf.eq.1) kfn = nz-1

      do k=kst,kfn
        z = ze+(k-1)*dlz
        do j=jst,jfn
          y = yc+(j-1)*dly
          do i=ist,ifn
            x = xa+(i-1)*dlx
            call coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
            cmin = dmin1(cmin,cxx,cyy,czz)
            cemax = dmax1(abs(ce),cemax)
c
c----+|-----------------------------------------------------------------------|
c     check if pde is "hyperbolic" at finest grid level
c----+|-----------------------------------------------------------------------|
            if (klevel.eq.ngrid) then
              if ((abs(cx)*dlx .gt. abs(cxx+cxx))) ier = -4
              if ((abs(cy)*dly .gt. abs(cyy+cyy))) ier = -4
              if ((abs(cz)*dlz .gt. abs(czz+czz))) ier = -4
            end if
c
c----+|-----------------------------------------------------------------------|
c     adjust second order coefficients so that pde is not "hyperbolic"
c     this is especially possible on coarser grids if there are non-zero
c     first order terms
c----+|-----------------------------------------------------------------------|
            cxx = dmax1(cxx,abs(cx)*dlx*0.5)
            cyy = dmax1(cyy,abs(cy)*dly*0.5)
            czz = dmax1(czz,abs(cz)*dlz*0.5)
            c1 = cxx/dlxx-cx/dlx2
            c2 = cxx/dlxx+cx/dlx2
            c3 = cyy/dlyy-cy/dly2
            c4 = cyy/dlyy+cy/dly2
            c5 = czz/dlzz-cz/dlz2
            c6 = czz/dlzz+cz/dlz2
            cof(i,j,k,1) = c1
            cof(i,j,k,2) = c2
            cof(i,j,k,3) = c3
            cof(i,j,k,4) = c4
            cof(i,j,k,5) = c5
            cof(i,j,k,6) = c6
            cof(i,j,k,7) = ce-(c1+c2+c3+c4+c5+c6)
          end do
        end do
      end do
c
c----+|-----------------------------------------------------------------------|
c     set nonfatal error flag if ellipticity test fails
c----+|-----------------------------------------------------------------------|
      if (ier .ne. -4) then
        if (cmin.le.0.d0) ier = -2
      end if
      alfmax = 0.d0
c
c----+|-----------------------------------------------------------------------|
c     adjust equation at mixed b.c.
c----+|-----------------------------------------------------------------------|
      if (nxa.eq.2) then
        kbdy = 1
        x = xa
        i = 1
        do k=kst,kfn
          z = ze+(k-1)*dlz
          do j=jst,jfn
            y = yc+(j-1)*dly
            call bndyc(kbdy,y,z,alfa,gbdy)
            alfmax = dmax1(abs(alfa),alfmax)
            c1 = cof(i,j,k,1)
            cof(i,j,k,1) = 0.d0
            cof(i,j,k,2) = cof(i,j,k,2)+c1
            cof(i,j,k,7) = cof(i,j,k,7)+dlx2*alfa*c1
          end do
        end do
      end if

      if (nxb.eq.2) then
        kbdy = 2
        x = xb
        i = nx
        do k=kst,kfn
          z = ze+(k-1)*dlz
          do j=jst,jfn
            y = yc+(j-1)*dly
            call bndyc(kbdy,y,z,alfa,gbdy)
            alfmax = dmax1(abs(alfa),alfmax)
            c2 = cof(i,j,k,2)
            cof(i,j,k,1) = cof(i,j,k,1)+c2
            cof(i,j,k,2) = 0.d0
            cof(i,j,k,7) = cof(i,j,k,7)-dlx2*alfa*c2
          end do
        end do
      end if

      if (nyc.eq.2) then
        kbdy = 3
        y = yc
        j = 1
        do k=kst,kfn
          z = ze+(k-1)*dlz
          do i=ist,ifn
            x = xa+(i-1)*dlx
            call bndyc(kbdy,x,z,alfa,gbdy)
            alfmax = dmax1(abs(alfa),alfmax)
            c3 = cof(i,j,k,3)
            cof(i,j,k,3) = 0.d0
            cof(i,j,k,4) = cof(i,j,k,4)+c3
            cof(i,j,k,7) = cof(i,j,k,7)+dly2*alfa*c3
          end do
        end do
      end if

      if (nyd.eq.2) then
        kbdy = 4
        y = yd
        j = ny
        do k=kst,kfn
        z = ze+(k-1)*dlz
        do i=ist,ifn
            x = xa+(i-1)*dlx
            call bndyc(kbdy,x,z,alfa,gbdy)
            alfmax = dmax1(abs(alfa),alfmax)
            c4 = cof(i,j,k,4)
            cof(i,j,k,3) = cof(i,j,k,3)+c4
            cof(i,j,k,4) = 0.d0
            cof(i,j,k,7) = cof(i,j,k,7)-dly2*c4*alfa
          end do
        end do
      end if

      if (nze.eq.2) then
        kbdy = 5
        z = ze
        k = 1
        do j=jst,jfn
          y = yc+(j-1)*dly
          do i=ist,ifn
            x = xa+(i-1)*dlx
            call bndyc(kbdy,x,y,alfa,gbdy)
            alfmax = dmax1(abs(alfa),alfmax)
            c5 = cof(i,j,k,5)
            cof(i,j,k,5) = 0.d0
            cof(i,j,k,6) = cof(i,j,k,6)+c5
            cof(i,j,k,7) = cof(i,j,k,7)+dlz2*c5*alfa
          end do
        end do
      end if

      if (nzf.eq.2) then
        kbdy = 6
        z = zf
        k = nz
        do j=jst,jfn
          y = yc+(j-1)*dly
          do i=ist,ifn
            x = xa+(i-1)*dlx
            call bndyc(kbdy,x,y,alfa,gbdy)
            alfmax = dmax1(abs(alfa),alfmax)
            c6 = cof(i,j,k,6)
            cof(i,j,k,5) = cof(i,j,k,5)+c6
            cof(i,j,k,6) = 0.d0
            cof(i,j,k,7) = cof(i,j,k,7)-dlz2*c6*alfa
          end do
        end do
      end if
c
c----+|-----------------------------------------------------------------------|
c     flag continuous singular elliptic pde if detected
c----+|-----------------------------------------------------------------------|
      if (ier .ne. -4) then
        if (cemax.eq.0.d0.and.alfmax.eq.0.d0) then
          if (nxa.eq.0.or.(nxa.eq.2.and.nxb.eq.2)) then
            if (nyc.eq.0.or.(nyc.eq.2.and.nyd.eq.2)) then
              if (nze.eq.0.or.(nze.eq.2.and.nzf.eq.2)) then
                ier = -3
              end if
            end if
          end if
        end if
      end if
c
c----+|-----------------------------------------------------------------------|
c     reset cof for specified b.c.
c----+|-----------------------------------------------------------------------|
      if (nxa.eq.1) then
        i = 1
        do j=1,ny
          do k=1,nz
            do l=1,7
              cof(i,j,k,l) = 0.d0
            end do
            cof(i,j,k,7) = 1.d0
          end do
        end do
      end if
      if (nxb.eq.1) then
        i = nx
        do k=1,nz
          do j=1,ny
            do l=1,7
              cof(i,j,k,l) = 0.d0
            end do
            cof(i,j,k,7) = 1.d0
          end do
        end do
      end if
      if (nyc.eq.1) then
        j = 1
        do k=1,nz
          do i=1,nx
            do l=1,7
              cof(i,j,k,l) = 0.d0
            end do
            cof(i,j,k,7) = 1.d0
          end do
        end do
      end if
      if (nyd.eq.1) then
        j = ny
        do i=1,nx
          do k=1,nz
            do l=1,7
              cof(i,j,k,l) = 0.d0
            end do
            cof(i,j,k,7) = 1.d0
          end do
        end do
      end if
      if (nze.eq.1) then
        k = 1
        do j=1,ny
          do i=1,nx
            do l=1,7
              cof(i,j,k,l) = 0.d0
            end do
            cof(i,j,k,7) = 1.d0
          end do
        end do
      end if
      if (nzf.eq.1) then
        k = nz
        do j=1,ny
          do i=1,nx
            do l=1,7
              cof(i,j,k,l) = 0.d0
            end do
            cof(i,j,k,7) = 1.d0
          end do
        end do
      end if
      if (method*(method-8)*(method-9)*(method-10) .eq. 0) return
c
c----+|-----------------------------------------------------------------------|
c     set,factor tridiagonal matrices for line relaxation
c----+|-----------------------------------------------------------------------|
      if ((method-1)*(method-4)*(method-5)*(method-7).eq.0) then
c
c----+|-----------------------------------------------------------------------|
c     line relaxation in x used
c----+|-----------------------------------------------------------------------|
c     set "non-periodic" tridiagonal matrices in tx and factor
c----+|-----------------------------------------------------------------------|
        if (nxa.ne.0) then
          do i=1,nx
            im1 = max0(i-1,1)
            do k=1,nz
              do j=1,ny
              tx(im1,j,k,1) = cof(i,j,k,1)
              tx(i,j,k,2) = cof(i,j,k,7)
              tx(i,j,k,3) = cof(i,j,k,2)
              end do
            end do
          end do
          nynz = ny*nz
          call factri(nynz,nx,tx(1,1,1,1),tx(1,1,1,2),tx(1,1,1,3))
        else
c----+|-----------------------------------------------------------------------|
c     set "periodic" tridiagonal matrices in tx and factor when nx > 3
c----+|-----------------------------------------------------------------------|
          if (nx .gt. 3) then
            do k=1,nz
              do j=1,ny
                do i=1,nx-1
                  tx(i,j,k,1) = cof(i,j,k,1)
                  tx(i,j,k,2) = cof(i,j,k,7)
                  tx(i,j,k,3) = cof(i,j,k,2)
                end do
              end do
            end do
            nynz = ny*nz
            call factrp(nynz,nx,tx,tx(1,1,1,2),tx(1,1,1,3),tx(1,1,1,4),
     +                  tx(1,1,1,5),wk(kps))
          end if
        end if
      end if

      if ((method-2)*(method-4)*(method-6)*(method-7).eq.0) then
c
c----+|-----------------------------------------------------------------------|
c     line relaxation in y used
c----+|-----------------------------------------------------------------------|
c     set "non-periodic" tridiagonal matrices and factor
c----+|-----------------------------------------------------------------------|
        if (nyc.ne.0) then
          do j=1,ny
            jm1 = max0(j-1,1)
            do k=1,nz
              do i=1,nx
                ty(jm1,i,k,1) = cof(i,j,k,3)
                ty(j,i,k,2) = cof(i,j,k,7)
                ty(j,i,k,3) = cof(i,j,k,4)
              end do
            end do
          end do
          nxnz = nx*nz
          call factri(nxnz,ny,ty(1,1,1,1),ty(1,1,1,2),ty(1,1,1,3))
        else
c----+|-----------------------------------------------------------------------|
c     set and factor periodic "tridiagonal" matrices when ny > 3
c----+|-----------------------------------------------------------------------|
          if (ny .gt. 3) then
            do k=1,nz
              do i=1,nx
                do j=1,ny-1
                  ty(j,i,k,1) = cof(i,j,k,3)
                  ty(j,i,k,2) = cof(i,j,k,7)
                  ty(j,i,k,3) = cof(i,j,k,4)
                end do
              end do
            end do
            nxnz = nx*nz
            call factrp(nxnz,ny,ty,ty(1,1,1,2),ty(1,1,1,3),ty(1,1,1,4),
     +                  ty(1,1,1,5),wk(kps))
          end if
        end if
      end if

      if ((method-3)*(method-5)*(method-6)*(method-7).eq.0) then
c
c----+|-----------------------------------------------------------------------|
c     line relaxation in z used
c----+|-----------------------------------------------------------------------|
c     set and factor non-periodic tridiagonal matrices
c----+|-----------------------------------------------------------------------|
        if (nze.ne.0) then
          do k=1,nz
            km1 = max0(k-1,1)
            do j=1,ny
              do i=1,nx
                tz(km1,i,j,1) = cof(i,j,k,5)
                tz(k,i,j,2) = cof(i,j,k,7)
                tz(k,i,j,3) = cof(i,j,k,6)
              end do
            end do
          end do
          nxny = nx*ny
          call factri(nxny,nz,tz(1,1,1,1),tz(1,1,1,2),tz(1,1,1,3))
        else
c----+|-----------------------------------------------------------------------|
c     set and factor periodic tridiagonal matrices when nz > 3
c----+|-----------------------------------------------------------------------|
          if (nz .gt. 3) then
            do j=1,ny
              do i=1,nx
                do k=1,nz-1
                  tz(k,i,j,1) = cof(i,j,k,5)
                  tz(k,i,j,2) = cof(i,j,k,7)
                  tz(k,i,j,3) = cof(i,j,k,6)
                end do
              end do
            end do
            nxny = nx*ny
            call factrp(nxny,nz,tz(1,1,1,1),tz(1,1,1,2),tz(1,1,1,3),
     +                  tz(1,1,1,4),tz(1,1,1,5),wk(kps))
          end if
        end if
      end if
      return
      end subroutine dismd3

c
c----+|-----------------------------------------------------------------------|
c     adjust rhs for solution in cof(i,j,k,8) on non-initial calls
c     (i.e., values in cof have not changed)
c----+|-----------------------------------------------------------------------|
      subroutine adjmd3(nx,ny,nz,phi,cof,bndyc,coef)
      implicit none
      integer nx,ny,nz
      double precision phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,8)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      double precision xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,ktzbgn,nxk,nyk,nzk,ngrid,
     +        klevel,kcur,kps
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/fmud3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/mud3c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +ktzbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
      double precision dlx,dly,dlz,dlx2,dly2,dlz2,dlxx,dlyy,dlzz
      double precision cxx,cyy,czz,cx,cy,cz,ce,alfa,gbdy,x,y,z,c1,c2,
     +c3,c4,c5,c6
      integer i,j,k,ist,ifn,jst,jfn,kst,kfn,kbdy
      external bndyc,coef
c
c----+|-----------------------------------------------------------------------|
c     set current grid increments
c----+|-----------------------------------------------------------------------|
      dlx = (xb-xa)/(nx-1)
      dlx2 = dlx+dlx
      dlxx = dlx*dlx
      dly = (yd-yc)/(ny-1)
      dly2 = dly+dly
      dlyy = dly*dly
      dlz = (zf-ze)/(nz-1)
      dlz2 = dlz+dlz
      dlzz = dlz*dlz
c
c----+|-----------------------------------------------------------------------|
c     set x,y,z subscript limits for calls to coef,bndyc
c----+|-----------------------------------------------------------------------|
      jst = 1
      jfn = ny
      ist = 1
      ifn = nx
      kst = 1
      kfn = nz
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.1) jst = 2
      if (nyd.eq.1) jfn = ny-1
      if (nze.eq.1) kst = 2
      if (nzf.eq.1) kfn = nz-1
c
c----+|-----------------------------------------------------------------------|
c     adjust for derivative b.c.
c----+|-----------------------------------------------------------------------|
      if (nxa.eq.2) then
        kbdy=1
        x=xa
        i = 1
        do k=kst,kfn
          z=ze+(k-1)*dlz
          do j=jst,jfn
            y=yc+(j-1)*dly
            call bndyc(kbdy,y,z,alfa,gbdy)
            call coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
            cxx = dmax1(cxx,abs(cx)*dlx*0.5)
            c1 = cxx/dlxx-cx/dlx2
            cof(i,j,k,8) = cof(i,j,k,8)+dlx2*c1*gbdy
          end do
        end do
      end if
      if (nxb.eq.2) then
        kbdy=2
        x=xb
        i = nx
        do k=kst,kfn
          z=ze+(k-1)*dlz
          do j=jst,jfn
            y=yc+(j-1)*dly
            call bndyc(kbdy,y,z,alfa,gbdy)
            call coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
            cxx = dmax1(cxx,abs(cx)*dlx*0.5)
            c2 = cxx/dlxx+cx/dlx2
            cof(i,j,k,8) = cof(i,j,k,8)-dlx2*c2*gbdy
          end do
        end do
      end if
      if (nyc.eq.2) then
        kbdy = 3
        y=yc
        j = 1
        do k=kst,kfn
          z=ze+(k-1)*dlz
          do i=ist,ifn
            x=xa+(i-1)*dlx
            call bndyc(kbdy,x,z,alfa,gbdy)
            call coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
            cyy = dmax1(cyy,abs(cy)*dly*0.5)
            c3 = cyy/dlyy-cy/dly2
            cof(i,j,k,8) = cof(i,j,k,8)+dly2*c3*gbdy
          end do
        end do
      end if
      if (nyd.eq.2) then
        kbdy = 4
        y=yd
        j = ny
        do k=kst,kfn
          z=ze+(k-1)*dlz
          do i=ist,ifn
            x=xa+(i-1)*dlx
            call bndyc(kbdy,x,z,alfa,gbdy)
            call coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
            cyy = dmax1(cyy,abs(cy)*dly*0.5)
            c4 = cyy/dlyy+cy/dly2
            cof(i,j,k,8) = cof(i,j,k,8)-dly2*c4*gbdy
          end do
        end do
      end if
      if (nze.eq.2) then
        kbdy = 5
        k = 1
        z=ze
        do j=jst,jfn
          y=yc+(j-1)*dly
          do i=ist,ifn
            x=xa+(i-1)*dlx
            call bndyc(kbdy,x,y,alfa,gbdy)
            call coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
            czz = dmax1(czz,abs(cz)*dlz*0.5)
            c5 = czz/dlzz-cz/dlz2
            cof(i,j,k,8) = cof(i,j,k,8)+dlz2*c5*gbdy
          end do
        end do
      end if
      if (nzf.eq.2) then
        kbdy = 6
        z=zf
        k = nz
        do j=jst,jfn
          y=yc+(j-1)*dly
          do i=ist,ifn
            x=xa+(i-1)*dlx
            call bndyc(kbdy,x,y,alfa,gbdy)
            call coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
            czz = dmax1(czz,abs(cz)*dlz*0.5)
            c6 = czz/dlzz+cz/dlz2
            cof(i,j,k,8) = cof(i,j,k,8)-dlz2*c6*gbdy
          end do
        end do
      end if
c
c----+|-----------------------------------------------------------------------|
c     set specified b.c.
c----+|-----------------------------------------------------------------------|
      if (nxa.eq.1) then
        i = 1
        do j=1,ny
          do k=1,nz
            cof(i,j,k,8) = phi(i,j,k)
          end do
        end do
      end if
      if (nxb.eq.1) then
        i = nx
        do j=1,ny
          do k=1,nz
            cof(i,j,k,8) = phi(i,j,k)
          end do
        end do
      end if
      if (nyc.eq.1) then
        j = 1
        do k=1,nz
          do i=1,nx
            cof(i,j,k,8) = phi(i,j,k)
          end do
        end do
      end if
      if (nyd.eq.1) then
        j = ny
        do k=1,nz
          do i=1,nx
            cof(i,j,k,8) = phi(i,j,k)
          end do
        end do
      end if
      if (nze.eq.1) then
        k = 1
        do j=1,ny
          do i=1,nx
            cof(i,j,k,8) = phi(i,j,k)
          end do
        end do
      end if
      if (nzf.eq.1) then
        k = nz
        do j=1,ny
          do i=1,nx
            cof(i,j,k,8) = phi(i,j,k)
          end do
        end do
      end if
      return
      end subroutine adjmd3

c----+|-----------------------------------------------------------------------|
c     use point or line relaxation in the x and/or y and/or z
c     or planar relaxation in the x,y or x,z or y,z planes
c----+|-----------------------------------------------------------------------|
      subroutine relmd3(wk)
      implicit none
      double precision wk(*)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      double precision xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,ktzbgn,nxk,nyk,nzk,ngrid,
     +        klevel,kcur,kps
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/fmud3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/mud3c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +ktzbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
      integer nx,ny,nz,ip,ic,m,itx,ity,itz
      nx = nxk(klevel)
      ny = nyk(klevel)
      nz = nzk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)

c----+|-----------------------------------------------------------------------|
c     gauss-seidel pointwise red/black relaxation
c----+|-----------------------------------------------------------------------|
      if (method.eq.0) then
        call relmd3p(nx,ny,nz,wk(ip),wk(ic))
        return
      end if

      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
      itz = ktzbgn(klevel)
      m = method
c----+|-----------------------------------------------------------------------|
c     check for line relaxation(s) (in combinations)
c----+|-----------------------------------------------------------------------|
c----+|-----------------------------------------------------------------------|
c     line - x relaxation
c     replace by point if x-periodic and nx=3
c----+|----------------------------------------------------------------|------|
      if ((m-1)*(m-4)*(m-5)*(m-7) .eq. 0 ) then
        if (nxa .ne. 0 .or. nx .gt. 3) then
          itx = ktxbgn(klevel)
          call slxmd3(nx,ny,nz,wk(ip),wk(ic),wk(itx),wk(kps),nxa,nyc,
     +                nze)
        else
          call relmd3p(nx,ny,nz,wk(ip),wk(ic))
        end if
        if (method .eq. 1) return
      end if
c
c----+|-----------------------------------------------------------------------|
c     line - y relaxation
c     replace by point if y-periodic and ny=3
c----+|-----------------------------------------------------------------------|
      if ((m-2)*(m-4)*(m-6)*(m-7) .eq. 0 ) then
        if (nyc .ne. 0 .or. ny .gt. 3) then
          ity = ktybgn(klevel)
          call slymd3(nx,ny,nz,wk(ip),wk(ic),wk(ity),wk(kps),nxa,nyc,
     +                nze)
        else
          call relmd3p(nx,ny,nz,wk(ip),wk(ic))
        end if
        if ((m-2)*(m-4) .eq. 0) return
      end if
c
c----+|-----------------------------------------------------------------------|
c     line - z relaxation
c     replace by point if z-periodic and nz=3
c----+|-----------------------------------------------------------------------|
      if ((m-3)*(m-5)*(m-6)*(m-7) .eq. 0 ) then
        if (nze .ne. 0 .or. nz .gt. 3) then
          itz = ktzbgn(klevel)
          call slzmd3(nx,ny,nz,wk(ip),wk(ic),wk(itz),wk(kps),nxa,nyc,
     +                nze)
        else
          call relmd3p(nx,ny,nz,wk(ip),wk(ic))
        end if
        return
      end if
c
c----+|-----------------------------------------------------------------------|
c     planar relaxation
c----+|-----------------------------------------------------------------------|
      if (method.eq.8) then
        call planxy(wk)
      else if (method.eq.9) then
        call planxz(wk)
      else if (method.eq.10) then
        call planyz(wk)
      end if
      return
      end subroutine relmd3

c
c----+|-----------------------------------------------------------------------|
c     gauss-seidel point relaxation with red/black ordering
c     in three dimensions for nonseparable pde
c----+|-----------------------------------------------------------------------|
      subroutine relmd3p(nx,ny,nz,phi,cof)
      implicit none
      integer nx,ny,nz
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      double precision phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,8)
      integer i,j,k,nper

      integer, external  :: omp_get_num_threads, omp_get_thread_num
c
c----+|-----------------------------------------------------------------------|
c     set periodic b.c. indicator
c----+|-----------------------------------------------------------------------|
      nper = nxa*nyc*nze
      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c----+|-----------------------------------------------------------------------|
c     relax in order:
c     (1) red (x,y) on odd z planes
c     (2) black (x,y) on even z planes
c     (3) black (x,y) on odd z planes
c     (4) red (x,y) on even z planes
c
c----+|-----------------------------------------------------------------------|
c     red (x,y) points on odd z planes
c----+|-----------------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(i,j,k), SHARED(phi,cof,nx,ny,nz)
      do k=1,nz,2
        do j=1,ny,2
          do i=1,nx,2
            phi(i,j,k) = (cof(i,j,k,8) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
          end do
        end do
        do j=2,ny,2
          do i=2,nx,2
            phi(i,j,k) = (cof(i,j,k,8) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
          end do
        end do
      end do
!$OMP END PARALLEL DO

      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c----+|-----------------------------------------------------------------------|
c     black (x,y) points on even z planes
c----+|-----------------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(i,j,k), SHARED(phi,cof,nx,ny,nz)
      do k=2,nz,2
        do j=2,ny,2
          do i=1,nx,2
            phi(i,j,k) = (cof(i,j,k,8) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
          end do
        end do
        do j=1,ny,2
          do i=2,nx,2
            phi(i,j,k) = (cof(i,j,k,8) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
          end do
        end do
      end do
!$OMP END PARALLEL DO

      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c----+|-----------------------------------------------------------------------|
c     black (x,y) points on odd z planes
c----+|-----------------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(i,j,k), SHARED(phi,cof,nx,ny,nz)
      do k=1,nz,2
        do j=2,ny,2
          do i=1,nx,2
            phi(i,j,k) = (cof(i,j,k,8) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
          end do
        end do
        do j=1,ny,2
          do i=2,nx,2
            phi(i,j,k) = (cof(i,j,k,8) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
          end do
        end do
      end do
!$OMP END PARALLEL DO

      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c----+|-----------------------------------------------------------------------|
c     red (x,y) points on even z planes
c----+|-----------------------------------------------------------------------|
!$OMP PARALLEL DO PRIVATE(i,j,k), SHARED(phi,cof,nx,ny,nz)
      do k=2,nz,2
        do j=1,ny,2
          do i=1,nx,2
            phi(i,j,k) = (cof(i,j,k,8) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
          end do
        end do
        do j=2,ny,2
          do i=2,nx,2
            phi(i,j,k) = (cof(i,j,k,8) - (
     +                    cof(i,j,k,1)*phi(i-1,j,k)+
     +                    cof(i,j,k,2)*phi(i+1,j,k)+
     +                    cof(i,j,k,3)*phi(i,j-1,k)+
     +                    cof(i,j,k,4)*phi(i,j+1,k)+
     +                    cof(i,j,k,5)*phi(i,j,k-1)+
     +                    cof(i,j,k,6)*phi(i,j,k+1)))
     +                   /cof(i,j,k,7)
          end do
        end do
      end do
!$OMP END PARALLEL DO
c
c----+|-----------------------------------------------------------------------|
c     final set of periodic virtual boundaries if necessary
c----+|-----------------------------------------------------------------------|
      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
      return
      end subroutine relmd3p

