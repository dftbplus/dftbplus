c
c     file muh3.f
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
c
c ... author and specialist
c
c          John C. Adams (National Center for Atmospheric Research)
c          email: johnad@ucar.edu, phone: 303-497-1213

c ... For MUDPACK 5.0 information, visit the website:
c     (http://www.scd.ucar.edu/css/software/mudpack)
c
c ... purpose (see muh3.d for details)
c
c     (muh3 is a hybrid multigrid/direct method solver)
c     muh3 attempts to produce a second order finite difference
c     approximation to the three dimensional nonseparable elliptic
c     partial differential equation of the form:
c
c       cxx(x,y,z)*pxx + cyy(x,y,z)*pyy + czz(x,y,z)*pzz +
c
c       cx(x,y,z)*px + cy(x,y,z)*py + cz(x,y,z)*pz +
c
c       ce(x,y,z)*p(x,y,z) = r(x,y,z)
c
c
c ... documentation and test files
c
c     see the documentation file "muh3.d" for a complete discussion
c     of how to use subroutine muh3.  file "tmuh3.f" is a test/driver
c     sample program illustrating use of muh3
c
c ... required MUDPACK files
c
c     mudcom.f, mud3ln.f, mud3pn.f
c
c
      subroutine muh3(iparm,fparm,wk,iwk,coef,bndyc,rhs,phi,mopt,ierror)
      implicit none
      integer iwk(*),iparm(23),mopt(4),ierror
      double precision wk(*),phi(*),rhs(*),fparm(8)
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
      integer ibeta,ialfa,izmat,idmat
      common/muh3c/ibeta,ialfa,izmat,idmat
      integer sxy,sxz,syz,kpspace,kcspace,int
      integer m,lxy,lxz,lyz,isx,jsy,ksz,iisx,jjsy,kksz
      integer nx,ny,nz,nxny,itx,ity,itz,k,kb,kkb,kk,ip,ic
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
c     set input parameters from iparm,fparm internally
c
      intl = iparm(1)
      nxa = iparm(2)
      nxb = iparm(3)
      nyc = iparm(4)
      nyd = iparm(5)
      nze = iparm(6)
      nzf = iparm(7)
c
c     set grid size params
c
      ixp = iparm(8)
      jyq = iparm(9)
      kzr = iparm(10)
      iex = iparm(11)
      jey = iparm(12)
      kez = iparm(13)
c
c     set number of subgrids for mg cycling
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
c     set floating point params
c
      xa = fparm(1)
      xb = fparm(2)
      yc = fparm(3)
      yd = fparm(4)
      ze = fparm(5)
      zf = fparm(6)
      tolmax = fparm(7)
c
c     set multigrid option parameters
c
      kcycle = mopt(1)
      if (kcycle .eq. 0) then
c
c     use default settings
c
      kcycle = 2
      iprer = 2
      ipost = 1
      intpol = 3
      else
      iprer = mopt(2)
      ipost = mopt(3)
      intpol = mopt(4)
      end if
      if (intl .eq. 0) then  ! intialization call
c
c     check input arguments
c
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
c       periodic size check
      if (nxa.eq.0 .and. ixp.lt.3) return
      if (nyc.eq.0 .and. jyq.lt.3) return
      if (nze.eq.0 .and. kzr.lt.3) return
      ierror = 4
      ngrid = max0(iex,jey,kez)
      if (iex.lt.1) return
      if (jey.lt.1) return
      if (kez.lt.1) return
      ierror = 13
      if ((iex-1)*(jey-1)*(kez-1) .eq. 0) then
        if (maxcy.gt.1) return
        if (iguess.eq.1) return
      end if
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
c     compute required work space length
c
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
c     set scales for planar relaxation
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
c     set subgrid sizes
c
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
        if (method.gt.7) then
c
c     set spacers for planar relaxation
c
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
c
c     set pointers
c
        kpbgn(k) = kps
        kcbgn(k) = kpbgn(k)+(nx+2)*(ny+2)*(nz+2)+kpspace
        ktxbgn(k) = kcbgn(k) + 8*nx*ny*nz + kcspace
        ktybgn(k) = ktxbgn(k) + isx*nx*ny*nz
        ktzbgn(k) = ktybgn(k) + jsy*nx*ny*nz
        kps = ktzbgn(k) + ksz*nx*ny*nz
c
c     sum space in case planar relaxation in xy or xz or yz
c
        sxy = sxy + (6+iisx+jjsy)*nx*ny + (nx+2)*(ny+2)
        sxz = sxz + (6+iisx+kksz)*nx*nz + (nx+2)*(nz+2)
        syz = syz + (6+jjsy+kksz)*ny*nz + (ny+2)*(nz+2)
      end do
c
c     set pointers for direct at coarse grid
c
      nx = ixp+1
      ny = jyq+1
      nz = kzr+1
      ibeta = kps+1
      nxny = nx*ny
      if (nze .ne. 0) then
c     nonperiodic z b.c.
        ialfa = ibeta + nxny*nxny*nz
        izmat = ialfa
        idmat = ialfa
        kps = ialfa+nxny*nxny*nz
      else
c     periodic z b.c.
        ialfa = ibeta + nxny*nxny*nz
        izmat = ialfa+nxny*nxny*nz
        idmat = izmat+nxny*nxny*(nz-2)
        kps = idmat+nxny*nxny*(nz-2)
      end if
c
c     set and check minimal work space
c
      nx = nxk(ngrid)
      ny = nyk(ngrid)
      nz = nzk(ngrid)
      iparm(22)=kps+max0((nx+2)*(ny+2)*(nz+2),lxy*sxy,lxz*sxz,lyz*syz)
      if (iparm(21) .lt. iparm(22)) return
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
c     discretize pde at each grid level
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
        klevel = k
        call dismh3(nx,ny,nz,wk(ic),wk(itx),
     +    wk(ity),wk(itz),bndyc,coef,wk,iwk,ierror)
        if (method.gt.7) then
c
c      discretize for planar coarsening
c
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
            call dismh3(nx,ny,nz,wk(ic),wk(itx),
     +        wk(ity),wk(itz),bndyc,coef,wk,iwk,ierror)
          end do
        end if
      end do
      return
      end if   ! end of intl=0 initialization call block
      nx = nfx
      ny = nfy
      nz = nfz
      call muh31(nx,ny,nz,rhs,phi,coef,bndyc,wk,iwk)
c
c     return if singular pde detected
c
      if (ierror .eq. 14) return
      iparm(23) = itero
      if (ierror .le. 0) then
      if (tolmax.gt.0.d0) then
c
c     set final computed maximum relative difference
c
        fparm(8) = relmax
        if (relmax.gt.tolmax .and. ierror.eq.0) then
c
c     flag convergence failure
c
          ierror = -1
        end if
      end if
      end if
      return
      end

      subroutine muh31(nx,ny,nz,rhsf,phif,coef,bndyc,wk,iwk)
      implicit none
      integer iwk(*)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      double precision xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,ktzbgn,nxk,nyk,nzk,ngrid,
     +        klevel,kcur,kps
      integer nx,ny,nz,ip,ic,ir,irc,ncx,ncy,ncz,ipc,icc
      integer i,j,k,kb,jk,kk,kkb,ijk,iter
      double precision phif(nx,ny,nz),rhsf(nx,ny,nz),wk(*),phmax
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/fmud3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/mud3c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +ktzbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
      integer ibeta,ialfa,izmat,idmat
      common/muh3c/ibeta,ialfa,izmat,idmat
      external coef,bndyc
      nx = nxk(ngrid)
      ny = nyk(ngrid)
      nz = nzk(ngrid)
      ip = kpbgn(ngrid)
      ic = kcbgn(ngrid)
      ir = ic+7*nx*ny*nz
c
c     set phif,rhsf in wk
c
      call swk3(nx,ny,nz,phif,rhsf,wk(ip),wk(ir))
      if (iguess.eq.0) then
c
c     no initial guess at finest grid level!
c
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
c     transfer down to all grid levels
c
        call trsfc3(nx,ny,nz,wk(ip),wk(ir),ncx,ncy,ncz,
     +                wk(ipc),wk(irc))
        if (method.gt.7) then
c
c     transfer down for planar coarsening
c
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
c     adjust right hand side at all grid levels in case
c     rhs or specified b.c. in phi or gbdy changed
c
      do kb=1,ngrid
        k = ngrid-kb+1
        nx = nxk(k)
        ny = nyk(k)
        nz = nzk(k)
        ip = kpbgn(k)
        ic = kcbgn(k)
        call adjmh3(nx,ny,nz,wk(ip),wk(ic),bndyc,coef)
c
c     adjust for planar grid levels if necessary
c
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
            call adjmh3(nx,ny,nz,wk(ip),wk(ic),bndyc,coef)
          end do
        end if
      end do
c
c     execute one full multigrid cycle
c
      do k=1,ngrid-1
        kcur = k
        call kcymh3(wk,iwk)
        nx = nxk(k+1)
        ny = nyk(k+1)
        nz = nzk(k+1)
        ip = kpbgn(k+1)
        ipc = kpbgn(k)
        ncx = nxk(k)
        ncy = nyk(k)
        ncz = nzk(k)
c
c     lift or prolong approximation from k to k+1
c
        call prolon3(ncx,ncy,ncz,wk(ipc),nx,ny,nz,wk(ip),
     +                 nxa,nxb,nyc,nyd,nze,nzf,intpol)
      end do

      else
c
c     adjust rhs at finest grid level only
c
      nx = nxk(ngrid)
      ny = nyk(ngrid)
      nz = nzk(ngrid)
      ip = kpbgn(ngrid)
      ic = kcbgn(ngrid)
      call adjmh3(nx,ny,nz,wk(ip),wk(ic),bndyc,coef)
      end if
c
c     execute maxcy more multigrid k cycles from finest level
c
      kcur = ngrid
      do iter=1,maxcy
      itero = iter
      call kcymh3(wk,iwk)
      if (tolmax.gt.0.d0) then
c
c      error control
c
        relmax = 0.d0
        phmax = 0.d0
        do k=1,nfz
          kk = k*(nfx+2)*(nfy+2)
          do j=1,nfy
            jk = kk+j*(nfx+2)
            do i=1,nfx
            ijk = jk+i+1
            phmax = max(phmax,abs(wk(ijk)))
            relmax = max(relmax,abs(wk(ijk)-phif(i,j,k)))
            phif(i,j,k) = wk(ijk)
            end do
          end do
        end do
c
c     set maximum relative difference and check for convergence
c
        if (phmax.gt.0.d0) relmax = relmax/phmax
        if (relmax.le.tolmax) return
      end if
      end do
c
c     set final iterate in phif
c
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
      end

      subroutine kcymh3(wk,iwk)
c
c     perform multigrid k-cycle at kcur level
c     kcycle = 1 corresponds to v cycles
c     kcycle = 2 corresponds to w cycles
c
      implicit none
      double precision wk(*)
      integer iwk(*)
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
      integer ibeta,ialfa,izmat,idmat
      common/muh3c/ibeta,ialfa,izmat,idmat
      integer nx,ny,nz,ncx,ncy,ncz,nxny
      integer kount(50),ip,ic,ipc,irc,nrel,l
      klevel = kcur
      if (kcur .eq. 1) then
      nx = nxk(1)
      ny = nyk(1)
      nz = nzk(1)
      nxny = nx*ny
      ip = kpbgn(1)
      ic = kcbgn(1)
c
c     solve at coarsest grid level with direct method and return
c
      if (nze .ne. 0) then
        call dir3(nxny,nx,ny,nz,wk(ip),wk(ic),wk(ibeta),wk(ialfa),
     +              wk(kps),iwk)
        return
      else
        call dir3p(nxny,nx,ny,nz,wk(ip),wk(ic),wk(ibeta),wk(ialfa),
     +               wk(izmat),wk(idmat),wk(kps),iwk)
        return
      end if
      end if
c
c     pre-relax at current finest grid level > 1
c
      do l = 1,iprer
      call relmh3(wk)
      end do
c
c     restrict residual to kcur-1
c
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
c
c     use full weighting with residual restriction
c
      call resmh3(nx,ny,nz,wk(ip),wk(ic),ncx,ncy,ncz,wk(ipc),wk(irc),
     +            wk(kps))
c
c     set counter for grid levels to zero
c
      do l = 1,kcur
      kount(l) = 0
      end do
c
c    set new level and continue k-cycling
c
      klevel = kcur-1
      nrel = iprer
c
c     kcycle control point
c
    1 continue
c
c     post-relax when kcur revisited
c
      if (klevel .eq. kcur) go to 2
c
c     count "hit" at current level
c
      kount(klevel) = kount(klevel)+1
c
c     relax at current level or use direct method if klevel=1
c
      if (klevel .gt. 1) then
      do l = 1,nrel
        call relmh3(wk)
      end do
      else
      nx = nxk(1)
      ny = nyk(1)
      nz = nzk(1)
      nxny = nx*ny
      ip = kpbgn(1)
      ic = kcbgn(1)
      if (nze .ne. 0) then
        call dir3(nxny,nx,ny,nz,wk(ip),wk(ic),wk(ibeta),wk(ialfa),
     +              wk(kps),iwk)
      else
        call dir3p(nxny,nx,ny,nz,wk(ip),wk(ic),wk(ibeta),wk(ialfa),
     +               wk(izmat),wk(idmat),wk(kps),iwk)
      end if
c
c     insure direct method is not called again at coarse level
c
      kount(1) = kcycle+1
      end if
      if (kount(klevel) .eq. kcycle+1) then
c
c     kcycle(iprer,ipost) complete at klevel
c     inject correction to finer grid
c
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
c     reset counter to zero at klevel
c
      kount(klevel) = 0
c
c     ascend to next higher level and set to post-relax there
c
      klevel = klevel+1
      nrel = ipost
      go to 1
      else
c
c     kcycle not complete so descend unless at coarsest
c
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
        call resmh3(nx,ny,nz,wk(ip),wk(ic),ncx,ncy,ncz,wk(ipc),
     +              wk(irc),wk(kps))
c
c     pre-relax at next coarser level
c
        klevel = klevel-1
        nrel = iprer
        go to 1
      else
c
c     direct method at coarsest level
c
        nx = nxk(1)
        ny = nyk(1)
        nz = nzk(1)
        nxny = nx*ny
        ip = kpbgn(1)
        ic = kcbgn(1)
        if (nze .ne. 0) then
          call dir3(nxny,nx,ny,nz,wk(ip),wk(ic),wk(ibeta),wk(ialfa),
     +                wk(kps),iwk)
        else
          call dir3p(nxny,nx,ny,nz,wk(ip),wk(ic),wk(ibeta),wk(ialfa),
     +                 wk(izmat),wk(idmat),wk(kps),iwk)
        end if
c
c     inject correction to grid level 2
c
        ipc = kpbgn(1)
        ncx = nxk(1)
        ncy = nyk(1)
        ncz = nzk(1)
        ip = kpbgn(2)
        nx = nxk(2)
        ny = nyk(2)
        nz = nzk(2)
        call cor3(nx,ny,nz,wk(ip),ncx,ncy,ncz,wk(ipc),
     +              nxa,nxb,nyc,nyd,nze,nzf,intpol,wk(kps))
c
c     set to post-relax at level 2
c
        nrel = ipost
        klevel = 2
        go to 1
      end if
      end if
    2 continue
c
c     post-relax at kcur level
c
      do l = 1,ipost
      call relmh3(wk)
      end do
      return
      end

      subroutine resmh3(nx,ny,nz,phi,cof,ncx,ncy,ncz,phic,rhsc,resf)
c
c     compute fully weighted residual restriction in rhsc
c
      implicit none
      integer nx,ny,nz,ncx,ncy,ncz,ic,jc,kc,i,j,k
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      double precision phi(0:nx+1,0:ny+1,0:nz+1)
      double precision phic(0:ncx+1,0:ncy+1,0:ncz+1)
      double precision rhsc(ncx,ncy,ncz),resf(nx,ny,nz),cof(nx,ny,nz,8)
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
c
c     restrict resf to interior coarse mesh in rhsc
c     using fully weighted residual restriction in 3-d
c
      call res3(nx,ny,nz,resf,ncx,ncy,ncz,rhsc,nxa,nxb,nyc,nyd,nze,nzf)
      return
      end

      subroutine dismh3(nx,ny,nz,cof,tx,ty,tz,bndyc,coef,wk,iwk,ier)
c
c     discretize the 3-d elliptic pde
c
      implicit none
      integer nx,ny,nz,ier
      double precision cof(nx,ny,nz,8)
      double precision tx(nx,ny,nz,*),ty(ny,nx,nz,*),tz(nz,nx,ny,*),
     +wk(*)
      integer iwk(*)
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
      integer ibeta,ialfa,izmat,idmat
      common/muh3c/ibeta,ialfa,izmat,idmat
      double precision dlx,dly,dlz,dlx2,dly2,dlz2,dlxx,dlyy,dlzz,cmin,
     +cemax,alfmax
      double precision cxx,cyy,czz,cx,cy,cz,ce,alfa,gbdy,x,y,z,c1,c2,c3,
     +c4,c5,c6
      integer i,j,k,l,ist,ifn,jst,jfn,kst,kfn,kbdy
      integer nxny,nxnz,nynz,im1,jm1,km1
      external bndyc,coef
c
c     set current grid increments
c
      dlx = (xb-xa)/(nx-1)
      dlx2 = dlx+dlx
      dlxx = dlx*dlx
      dly = (yd-yc)/(ny-1)
      dly2 = dly+dly
      dlyy = dly*dly
      dlz = (zf-ze)/(nz-1)
      dlz2 = dlz+dlz
      dlzz = dlz*dlz
      cmin = 1.0
      cemax = 0.d0
c
c     set x,y,z subscript limits to bypass specified boundaries
c     when calling coef or bndyc
c
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
          cmin = min(cmin,cxx,cyy,czz)
          cemax = max(abs(ce),cemax)
c
c     check if pde is "hyperbolic" at finest grid level
c
          if (klevel.eq.ngrid) then
            if ((abs(cx)*dlx .gt. abs(cxx+cxx))) ier = -4
            if ((abs(cy)*dlx .gt. abs(cyy+cyy))) ier = -4
            if ((abs(cz)*dlx .gt. abs(czz+czz))) ier = -4
          end if
c
c     adjust second order coefficients so that pde is not "hyperbolic"
c     this is especially possible on coarser grids if there are non-zero
c     first order terms
c
          cxx = max(cxx,abs(cx)*dlx*0.5d0)
          cyy = max(cyy,abs(cy)*dly*0.5d0)
          czz = max(czz,abs(cz)*dlz*0.5d0)
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
      if (ier .ne. -4) then
c
c     set nonfatal error flag if ellipticity test fails
c
      if (cmin.le.0.d0) ier = -2
      end if
      alfmax = 0.d0
c
c     adjust equation at mixed b.c.
c
      if (nxa.eq.2) then
      kbdy = 1
      x = xa
      i = 1
      do k=kst,kfn
        z = ze+(k-1)*dlz
        do j=jst,jfn
          y = yc+(j-1)*dly
          call bndyc(kbdy,y,z,alfa,gbdy)
          alfmax = max(abs(alfa),alfmax)
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
          alfmax = max(abs(alfa),alfmax)
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
          alfmax = max(abs(alfa),alfmax)
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
          alfmax = max(abs(alfa),alfmax)
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
          alfmax = max(abs(alfa),alfmax)
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
          alfmax = max(abs(alfa),alfmax)
          c6 = cof(i,j,k,6)
          cof(i,j,k,5) = cof(i,j,k,5)+c6
          cof(i,j,k,6) = 0.d0
          cof(i,j,k,7) = cof(i,j,k,7)-dlz2*c6*alfa
        end do
      end do
      end if
c
c     flag continuous singular elliptic pde if detected
c     a fatal error for muh3
c
      if (cemax.eq.0.d0.and.alfmax.eq.0.0) then
      if (nxa.eq.0.or.(nxa.eq.2.and.nxb.eq.2)) then
        if (nyc.eq.0.or.(nyc.eq.2.and.nyd.eq.2)) then
          if (nze.eq.0.or.(nze.eq.2.and.nzf.eq.2)) then
            ier =14
          end if
        end if
      end if
      end if
c
c     reset cof for specified b.c.
c
      if (nxa.eq.1) then
      i = 1
      do j=1,ny
        do k=1,nz
          do l=1,7
            cof(i,j,k,l) = 0.d0
          end do
          cof(i,j,k,7) = 1.0
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
          cof(i,j,k,7) = 1.0
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
          cof(i,j,k,7) = 1.0
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
          cof(i,j,k,7) = 1.0
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
          cof(i,j,k,7) = 1.0
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
          cof(i,j,k,7) = 1.0
        end do
      end do
      end if
c     if (klevel.eq.1 .and. method.lt.8) then
      if (klevel.eq.1) then

c
c     set coarse grid resolution
c
      nx = ixp+1
      ny = jyq+1
      nz = kzr+1
      nxny = nx*ny
      if (nze .ne. 0) then
c
c     factor non-periodic block matrix
c
        call lud3(nxny,nx,ny,nz,cof,wk(ibeta),wk(ialfa),iwk,nxa,nyc)
        return
      else
c
c     factor periodic block matrix
c
        do k=1,nz-1
          call setbet3(nxny,nx,ny,nz,cof,wk(ibeta),k,nxa,nyc)
          call setalf3(nxny,nx,ny,nz,cof,wk(ialfa),k)
        end do
        call lud3p(nxny,nx,ny,nz,cof,wk(ibeta),wk(ialfa),wk(izmat),
     +               wk(idmat),iwk,nxa,nyc)
        return
      end if
      end if
      if (method*(method-8)*(method-9)*(method-10) .eq. 0) return
c
c     set,factor tridiagonal matrices for line relaxation
c
      if ((method-1)*(method-4)*(method-5)*(method-7).eq.0) then
c
c     line relaxation in x used
c
      if (nxa.ne.0) then
c
c     set non-periodic tridiagonal matrices in tx and factor
c
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
        if (nx .gt. 3) then
c
c     set "periodic" tridiagonal matrices in tx and factor when nx > 3
c
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
c     line relaxation in y used
c
      if (nyc.ne.0) then
c
c     set non-periodic tridiagonal matrices and factor
c
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
        if (ny .gt. 3) then
c
c     set and factor periodic "tridiagonal" matrices when ny > 3
c
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
c     line relaxation in z used
c
      if (nze.ne.0) then
c
c     set and factor non-periodic tridiagonal matrices
c
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
        if (nz .gt. 3) then
c
c     set and factor periodic "tridiagonal matrices when nz > 3"
c
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
      end

      subroutine adjmh3(nx,ny,nz,phi,cof,bndyc,coef)
c
c     adjust rhs for solution in cof(i,j,k,8) on non-initial calls
c     (i.e., values in cof have not changed)
c
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
      integer ibeta,ialfa,izmat,idmat
      common/muh3c/ibeta,ialfa,izmat,idmat
      double precision dlx,dly,dlz,dlx2,dly2,dlz2,dlxx,dlyy,dlzz
      double precision cxx,cyy,czz,cx,cy,cz,ce,alfa,gbdy,x,y,z,c1,c2,c3,
     +c4,c5,c6
      integer i,j,k,ist,ifn,jst,jfn,kst,kfn,kbdy
      external bndyc,coef
c
c     set current grid increments
c
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
c     set x,y,z subscript limits for calls to coef,bndyc
c
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
c     adjust for derivative b.c.
c
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
          cxx = max(cxx,abs(cx)*dlx*0.5d0)
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
          cxx = max(cxx,abs(cx)*dlx*0.5d0)
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
          cyy = max(cyy,abs(cy)*dly*0.5d0)
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
          cyy = max(cyy,abs(cy)*dly*0.5d0)
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
          czz = max(czz,abs(cz)*dlz*0.5d0)
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
          czz = max(czz,abs(cz)*dlz*0.5d0)
          c6 = czz/dlzz+cz/dlz2
          cof(i,j,k,8) = cof(i,j,k,8)-dlz2*c6*gbdy
        end do
      end do
      end if
c
c     set specified b.c.
c
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
      end

      subroutine relmh3(wk)
c
c     use point or line relaxation in the x and/or y and/or z
c     or planar relaxation in the x,y or x,z or y,z planes
c
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
      integer ibeta,ialfa,izmat,idmat
      common/muh3c/ibeta,ialfa,izmat,idmat
      integer nx,ny,nz,ip,ic,m,itx,ity,itz
      nx = nxk(klevel)
      ny = nyk(klevel)
      nz = nzk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      if (method.eq.0) then
c
c     gauss-seidel pointwise red/black relaxation
c
      call relmh3p(nx,ny,nz,wk(ip),wk(ic))
      return
      end if
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
      itz = ktzbgn(klevel)
      m = method
c
c     check for line relaxation(s) (in combinations)
c
      if ((m-1)*(m-4)*(m-5)*(m-7) .eq. 0 ) then
c
c     line - x relaxation
c
      if (nxa .ne. 0 .or. nx .gt. 3) then
       itx = ktxbgn(klevel)
       call slxmd3(nx,ny,nz,wk(ip),wk(ic),wk(itx),wk(kps),nxa,nyc,nze)
      else
c
c     replace by point if x-periodic and nx=3
c
        call relmh3p(nx,ny,nz,wk(ip),wk(ic))
      end if
      if (method .eq. 1) return
      end if

      if ((m-2)*(m-4)*(m-6)*(m-7) .eq. 0 ) then
c
c     line - y relaxation
c
      if (nyc .ne. 0 .or. ny .gt. 3) then
       ity = ktybgn(klevel)
       call slymd3(nx,ny,nz,wk(ip),wk(ic),wk(ity),wk(kps),nxa,nyc,nze)
      else
c
c     replace by point if y-periodic and ny=3
c
        call relmh3p(nx,ny,nz,wk(ip),wk(ic))
      end if
      if ((m-2)*(m-4) .eq. 0) return
      end if

      if ((m-3)*(m-5)*(m-6)*(m-7) .eq. 0 ) then
c
c     line - z relaxation
c
      if (nze .ne. 0 .or. nz .gt. 3) then
       itz = ktzbgn(klevel)
       call slzmd3(nx,ny,nz,wk(ip),wk(ic),wk(itz),wk(kps),nxa,nyc,nze)
      else
c
c     replace by point if z-periodic and nz=3
c
      call relmh3p(nx,ny,nz,wk(ip),wk(ic))
      end if
      return
      end if
c
c     planar relaxation
c
      if (method.eq.8) then
c
c     planar relaxation in the x,y plane
c
      call planxy(wk)
      else if (method.eq.9) then
c
c     planar relaxation in the x,z plane
c
      call planxz(wk)
      else if (method.eq.10) then
c
c     planar relaxation in the y,z plane
c
      call planyz(wk)
      end if
      return
      end

      subroutine relmh3p(nx,ny,nz,phi,cof)
c
c     gauss-seidel point relaxation with red/black ordering
c     in three dimensions for nonseparable pde
c
      implicit none
      integer nx,ny,nz,i,j,k,nper
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      double precision phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,8)
c
c     set periodic b.c. indicator
c
      nper = nxa*nyc*nze
c
c     set periodic boundaries as necessary
c
      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     relax in order:
c     (1) red (x,y) on odd z planes
c     (2) black (x,y) on even z planes
c     (3) black (x,y) on odd z planes
c     (4) red (x,y) on even z planes
c
!$OMP PARALLEL DO PRIVATE(i,j,k), SHARED(phi,cof,nx,ny,nz)
      do k=1,nz,2
c
c     red (x,y) points on odd z planes
c
      do i=1,nx,2
        do j=1,ny,2
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
      do i=2,nx,2
        do j=2,ny,2
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
c    black (x,y) points on even z planes
c
!$OMP PARALLEL DO PRIVATE(i,j,k), SHARED(phi,cof,nx,ny,nz)
      do k=2,nz,2
      do i=1,nx,2
        do j=2,ny,2
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
      do i=2,nx,2
        do j=1,ny,2
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
c     black (x,y) points on odd z planes
c
!$OMP PARALLEL DO PRIVATE(i,j,k), SHARED(phi,cof,nx,ny,nz)
      do k=1,nz,2
      do i=1,nx,2
        do j=2,ny,2
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
      do i=2,nx,2
        do j=1,ny,2
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
c     red (x,y) points on even z planes
c
!$OMP PARALLEL DO PRIVATE(i,j,k), SHARED(phi,cof,nx,ny,nz)
      do k=2,nz,2
      do i=1,nx,2
        do j=1,ny,2
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
      do i=2,nx,2
        do j=2,ny,2
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
c     final set of periodic virtual boundaries if necessary
c
      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
      return
      end

      subroutine lud3(nxny,nx,ny,nz,cof,beta,alfa,index,nxa,nyc)
c
c     decompose nonperiodic block coefficient matrix
c     coming from 3-d  PDE discretization
c
      implicit none
      integer nxny,nx,ny,nz,nxa,nyc
      double precision cof(nx,ny,nz,8),beta(nxny,nxny,*)
      double precision alfa(nxny,nxny,*)
      integer index(nxny,nz)
      integer iz,i1,kcur,ii,ll,il,jl,km1
      iz = 0
      i1 = 1
c
c     set and factor umat(1) in beta(1)
c
      kcur = 1
      call setbet3(nxny,nx,ny,nz,cof,beta,kcur,nxa,nyc)
      call sgfa(beta,nxny,nxny,index,iz)
      do kcur=2,nz
      call setalf3(nxny,nx,ny,nz,cof,alfa,kcur)
      call transp(nxny,alfa(1,1,kcur))
      km1 = kcur-1
      do ll=1,nxny
        call sgsl(beta(1,1,km1),nxny,nxny,index(1,km1),
     +              alfa(1,ll,kcur),i1)
      end do
      call transp(nxny,alfa(1,1,kcur))
c
C     set beta(kcur)=beta(kcur)-alfa(kcur)*gama(kcur-1)
c
      call setbet3(nxny,nx,ny,nz,cof,beta,kcur,nxa,nyc)
      do ii=1,nxny
        do ll=1,nxny
          jl = ((ll-1)/nx) + 1
          il = ll - (jl-1)*nx
          beta(ii,ll,kcur)=beta(ii,ll,kcur)-alfa(ii,ll,kcur)*
     +                       cof(il,jl,kcur-1,6)
        end do
      end do
c
c     factor current beta for next pass
c
      iz = 0
      call sgfa(beta(1,1,kcur),nxny,nxny,index(1,kcur),iz)
      end do
c
C     matrix factorization now complete
c
      return
      end

      subroutine dir3(nxny,nx,ny,nz,phi,cof,beta,alfa,phxy,index)
c
c     direct solve at coarsest grid
c
      implicit none
      integer nxny,nx,ny,nz,index(nxny,ny)
      double precision phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,8)
      double precision phxy(nxny)
      double precision beta(nxny,nxny,*),alfa(nxny,nxny,*)
c     forward sweep
      call for3(nxny,nx,ny,nz,phi,cof(1,1,1,8),alfa)
c     backward sweep
      call bkw3(nxny,nx,ny,nz,phi,cof,beta,phxy,index)
      return
      end

      subroutine for3(nxny,nx,ny,nz,phi,frhs,alfa)
c
c     forward sweep
c
      implicit none
      integer nxny,nx,ny,nz,i,j,k,ii,ll,il,jl
      double precision phi(0:nx+1,0:ny+1,0:nz+1)
      double precision frhs(nx,ny,nz),alfa(nxny,nxny,*)
      double precision sum
      do k=1,nz
      do j=1,ny
        do i=1,nx
          phi(i,j,k) = frhs(i,j,k)
        end do
      end do
      end do
      do k=2,nz
c
C     SET PHI(k)=PHI(k)-alfa(k)*PHI(k-1)
c
      do ii=1,nxny
        j = ((ii-1)/nx)+1
        i = ii - (j-1)*nx
        sum=0.d0
        do ll=1,nxny
          jl = ((ll-1)/nx)+1
          il = ll - (jl-1)*nx
          sum = sum+alfa(ii,ll,k)*phi(il,jl,k-1)
        end do
        phi(i,j,k) = phi(i,j,k)-sum
      end do
      end do
      return
      end

      subroutine bkw3(nxny,nx,ny,nz,phi,cof,beta,phxy,index)
C                                                                               
C     backward sweep
c
      implicit none
      integer nxny,nx,ny,nz,index(nxny,nz)
      double precision phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,8)
      double precision beta(nxny,nxny,*),phxy(nxny)
      integer iz,i,j,k,kb,ij
C
C     solve beta*phi(nz)=phi(nz)
C
      iz = 0
      do j=1,ny
      do i=1,nx
        ij = (j-1)*nx+i
        phxy(ij) = phi(i,j,nz)
      end do
      end do
      call sgsl(beta(1,1,nz),nxny,nxny,index(1,nz),phxy,iz)
      do j=1,ny
      do i=1,nx
        ij = (j-1)*nx+i
        phi(i,j,nz) = phxy(ij)
      end do
      end do
      do kb=2,nz
      k = nz-kb+1
C
C     set current rhs
C
      do j=1,ny
        do i=1,nx
          phi(i,j,k) = phi(i,j,k) - cof(i,j,k,6)*phi(i,j,k+1)
        end do
      end do
C
C     solve beta(k)*phi(k)=phi(k)
C
      do j=1,ny
        do i=1,nx
          ij = (j-1)*nx+i
          phxy(ij) = phi(i,j,k)
        end do
      end do
      call sgsl(beta(1,1,k),nxny,nxny,index(1,k),phxy,iz)
      do j=1,ny
        do i=1,nx
          ij = (j-1)*nx+i
          phi(i,j,k) = phxy(ij)
        end do
      end do
      end do
      return
      end

      subroutine lud3p(nxny,nx,ny,nz,cof,beta,alfa,zmat,dmat,index,
     +                 nxa,nyc)
C                                                                               
C     DO LUD DECOMPOSITION OF BLOCK periodic TRIDIAGONAL MATRIX EQUATION
c
      implicit none
      integer nxny,nx,ny,nz,index(nxny,nz),nxa,nyc
      double precision cof(nx,ny,nz,8),alfa(nxny,nxny,*),
     +beta(nxny,nxny,*)
      double precision dmat(nxny,nxny,*),zmat(nxny,nxny,*),sum
      integer iz,kcur,ii,jj,ll,i1,km1,il,jl,i,j,ij
      kcur = 1
c
c     set dmat(1)=alfa(1)
c
      call setalf3(nxny,nx,ny,nz,cof,alfa,kcur)
      do ii=1,nxny
      do ll=1,nxny
        dmat(ii,ll,1) = alfa(ii,ll,kcur)
      end do
      end do
c
c     factor umat(1) in beta(1)
c
      call setbet3(nxny,nx,ny,nz,cof,beta,kcur,nxa,nyc)
      iz = 0
      call sgfa(beta(1,1,1),nxny,nxny,index(1,1),iz)
      do kcur=2,nz-2
c
c     solve transpose of lmat(kcur)umat(kcur-1)=alfa(kcur) in alfa(kcur)
c
      call setalf3(nxny,nx,ny,nz,cof,alfa,kcur)
c
c     transpose alfa
c
      call transp(nxny,alfa(1,1,kcur))
c
c     solve transpose equation
c
      km1 = kcur-1
      i1 = 1
      do ll=1,nxny
        call sgsl(beta(1,1,km1),nxny,nxny,index(1,km1),alfa(1,ll,kcur)
     +              ,i1)
      end do
c
c     re-transpose solution in alfa
c
      call transp(nxny,alfa(1,1,kcur))
c
C     set umat(kcur) = beta(kcur) - alfa(kcur)*gama(kcur-1) in beta(kcur)
c
      call setbet3(nxny,nx,ny,nz,cof,beta,kcur,nxa,nyc)
      do ii=1,nxny
        do ll=1,nxny
          jl = ((ll-1)/nx)+1
          il = ll-(jl-1)*nx
          beta(ii,ll,kcur)=beta(ii,ll,kcur)-alfa(ii,ll,kcur)*
     +                       cof(il,jl,kcur-1,6)
        end do
      end do
c
c     factor current beta(1,1,kcur) for next pass
c
      call sgfa(beta(1,1,kcur),nxny,nxny,index(1,kcur),iz)
c
c     set dmat(kcur) = -alfa(kcur)*dmat(kcur-1)
c
      do ii=1,nxny
        do jj=1,nxny
          dmat(ii,jj,kcur) = 0.d0
          do ll=1,nxny
            dmat(ii,jj,kcur)=dmat(ii,jj,kcur)-alfa(ii,ll,kcur)*
     +                         dmat(ll,jj,kcur-1)
          end do
        end do
      end do
      if (kcur .eq. nz-2) then
c
c     adjust dmat(nz-2) = gama(nz-2)-alfa(nz-2)*dmat(nz-2)
c
        do j=1,ny
          do i=1,nx
            ij = (j-1)*nx+i
            dmat(ij,ij,kcur) = cof(i,j,kcur,6)+dmat(ij,ij,kcur)
          end do
        end do
      end if
      end do
c
c     final phase with periodic factorization
c
c     solve transpose of zmat(1) beta(1) = gama(nz-1)
c
      do ii=1,nxny
      j = ((ii-1)/nx)+1
      i = ii-(j-1)*nx
      do jj=1,nxny
        zmat(jj,ii,1) = 0.d0
      end do
      zmat(ii,ii,1) = cof(i,j,nz-1,6)
      end do
      do ll=1,nxny
      call sgsl(beta(1,1,1),nxny,nxny,index(1,1),zmat(1,ll,1),i1)
      end do
      call transp(nxny,zmat(1,1,1))
c
c     solve transpose of zmat(k) umat(k) = -zmat(k-1) gama(k-1)
c
      do kcur = 2,nz-3
      do ii=1,nxny
        do jj=1,nxny
          j = ((jj-1)/nx)+1
          i = jj-(j-1)*nx
          zmat(jj,ii,kcur) = -zmat(ii,jj,kcur-1)*cof(i,j,kcur-1,6)
        end do
      end do
      do ll=1,nxny
        call sgsl(beta(1,1,kcur),nxny,nxny,index(1,kcur),
     +               zmat(1,ll,kcur),i1)
      end do
      call transp(nxny,zmat(1,1,kcur))
      end do
c
c     solve transpose of zmat(nz-2)umat(nz-2)=alfa(nz-1)-zmat(nz-3)gama(nz-3)
c
      call setalf3(nxny,nx,ny,nz,cof,alfa,nz-1)
      do ii=1,nxny
      do jj=1,nxny
        j = ((jj-1)/nx)+1
        i = jj-(j-1)*nx
        zmat(jj,ii,nz-2) = alfa(ii,jj,nz-1)-zmat(ii,jj,nz-3)*
     +                   cof(i,j,nz-3,6)
      end do
      end do
      do ll=1,nxny
      call sgsl(beta(1,1,nz-2),nxny,nxny,index(1,nz-2),
     +             zmat(1,ll,nz-2),i1)
      end do
      call transp(nxny,zmat(1,1,nz-2))
c
c     set umat(nz-1) = beta(nz-1)-(zmat(1)*dmat(1)+...+zmat(nz-2)*dmat(nz-2))
c     in beta(nz-1)
c
      call setbet3(nxny,nx,ny,nz,cof,beta,nz-1,nxa,nyc)
      do ii=1,nxny
      do jj=1,nxny
        sum = 0.d0
        do kcur=1,nz-2
          do ll=1,nxny
            sum = sum + zmat(ii,ll,kcur)*dmat(ll,jj,kcur)
          end do
        end do
        beta(ii,jj,nz-1) = beta(ii,jj,nz-1) - sum
      end do
      end do
c
c     factor bmat(nz-1) for backward sweep
c
      call sgfa(beta(1,1,nz-1),nxny,nxny,index(1,nz-1),iz)
c
c     lud is now complete
c
      return
      end

      subroutine dir3p(nxny,nx,ny,nz,phi,cof,beta,alfa,zmat,dmat,
     +                 phxy,index)
c
C     direct solution of periodic block matrix at coarsest level
c
      implicit none
      integer nxny,nx,ny,nz,index(nxny,nz)
      double precision phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,8)
      double precision beta(nxny,nxny,*),alfa(nxny,nxny,*)
      double precision zmat(nxny,nxny,*), dmat(nxny,nxny,*)
      double precision phxy(nxny)
c     forward sweep for periodic
      call for3p(nxny,nx,ny,nz,phi,cof(1,1,1,8),alfa,zmat)
c     backward sweep for periodic
      call bkw3p(nxny,nx,ny,nz,phi,cof,beta,dmat,phxy,index)
      return
      end

      subroutine for3p(nxny,nx,ny,nz,phi,frhs,alfa,zmat)
C                                                                               
C     DO FORWARD SWEEP PHASE OF DIRECT METHOD IN SOLUTION OF
C     periodic FACTORED BLOCK TRIDIAGONAL SYSTEM
C                                                                               
      implicit none
      integer nxny,nx,ny,nz
      double precision phi(0:nx+1,0:ny+1,0:nz+1)
      double precision frhs(nx,ny,nz)
      double precision alfa(nxny,nxny,*),zmat(nxny,nxny,*),sum
      integer i,j,k,kcur,ii,ll,il,jl
C     SET RHS IN PHI ADJUSTED AT DIRICHLET BOUNDARIES                           
      do k=1,nz-1
      do j=1,ny
        do i=1,nx
          phi(i,j,k) = frhs(i,j,k)
        end do
      end do
      end do
      do kcur=2,nz-2
      do ii=1,nxny
        j = ((ii-1)/nx)+1
        i = ii-(j-1)*nx
        sum = 0.d0
        do ll=1,nxny
          jl = ((ll-1)/nx)+1
          il = ll - (jl-1)*nx
          sum=sum+alfa(ii,ll,kcur)*phi(il,jl,kcur-1)
        end do
        phi(i,j,kcur) = phi(i,j,kcur)-sum
      end do
      end do
c
c     solve:
c     zmat(1)*phi(1)+...+zmat(nz-2)*phi(nz-2) + phi(nz-1) = f(nz-1)
c
      do ii=1,nxny
      j = ((ii-1)/nx)+1
      i = ii-(j-1)*nx
      sum = 0.d0
      do k=1,nz-2
        do ll=1,nxny
          jl = ((ll-1)/nx)+1
          il = ll - (jl-1)*nx
          sum = sum + zmat(ii,ll,k)*phi(il,jl,k)
        end do
      end do
      phi(i,j,nz-1) = phi(i,j,nz-1) - sum
      end do
      return
      end

      subroutine bkw3p(nxny,nx,ny,nz,phi,cof,beta,dmat,phxy,index)
C                                                                               
C     DO BACKWARD SWEEP
c
      implicit none
      integer nxny,nx,ny,nz,index(nxny,nz)
      double precision phi(0:nx+1,0:ny+1,0:nz+1),cof(nx,ny,nz,8)
      double precision beta(nxny,nxny,nz),dmat(nxny,nxny,nz)
      double precision phxy(nxny),sum
      integer iz,ii,ll,il,jl,k,kb,i,j,ij
C
C     solve beta(nz-1)*phi(nz-1) = phi(nz-1)
C
      iz = 0
      do j=1,ny
      do i=1,nx
        ij = (j-1)*nx+i
        phxy(ij) = phi(i,j,nz-1)
      end do
      end do

      call sgsl(beta(1,1,nz-1),nxny,nxny,index(1,nz-1),phxy,iz)
      do j=1,ny
      do i=1,nx
        ij = (j-1)*nx+i
        phi(i,j,nz-1) = phxy(ij)
      end do
      end do
c
c     solve beta(nz-2)*phi(nz-2) = phi(nz-2)-dmat(nz-2)*phi(nz-1)
c
      do ii=1,nxny
      j = ((ii-1)/nx)+1
      i = ii-(j-1)*nx
      sum = 0.d0
      do ll=1,nxny
        jl = ((ll-1)/nx)+1
        il = ll - (jl-1)*nx
        sum = sum + dmat(ii,ll,nz-2)*phi(il,jl,nz-1)
      end do
      phi(i,j,nz-2) = phi(i,j,nz-2) - sum
      end do
      do j=1,ny
      do i=1,nx
        ij = (j-1)*nx+i
        phxy(ij) = phi(i,j,nz-2)
      end do
      end do
      call sgsl(beta(1,1,nz-2),nxny,nxny,index(1,nz-2),phxy,iz)
      do j=1,ny
      do i=1,nx
        ij = (j-1)*nx+i
        phi(i,j,nz-2) = phxy(ij)
      end do
      end do
c
c     solve beta(k)*phi(k) = phi(k) - gama(k)*phi(k+1)-dmat(k)*phi(nz-1)
c     k=nz-3,...,1
c
      do kb=4,nz
      k = nz-kb+1
      do ii=1,nxny
        j = ((ii-1)/nx)+1
        i = ii-(j-1)*nx
        sum = 0.d0
        do ll=1,nxny
          jl = ((ll-1)/nx)+1
          il = ll - (jl-1)*nx
          sum = sum+dmat(ii,ll,k)*phi(il,jl,nz-1)
        end do
        phi(i,j,k) = phi(i,j,k) - sum - cof(i,j,k,6)*phi(i,j,k+1)
      end do
      do j=1,ny
        do i=1,nx
          ij = (j-1)*nx+i
          phxy(ij) = phi(i,j,k)
        end do
      end do
      call sgsl(beta(1,1,k),nxny,nxny,index(1,k),phxy,iz)
      do j=1,ny
        do i=1,nx
          ij = (j-1)*nx+i
          phi(i,j,k) = phxy(ij)
        end do
      end do
      end do
c
c     set k=nz by periodicity
c
      do j=1,ny
      do i=1,nx
        phi(i,j,nz) = phi(i,j,1)
      end do
      end do
      return
      end

      subroutine setbet3(nxny,nx,ny,nz,cof,beta,kcur,nxa,nyc)
c
c     set diagonal matrix on block for k = kcur
c
      implicit none
      integer nxny,nx,ny,nz,nxa,nyc
      double precision cof(nx,ny,nz,8),beta(nxny,nxny,*)
      integer k,kcur,ii,ll,ij,jj,i,j
      k = kcur
      do ii = 1,nxny
      do ll=1,nxny
        beta(ii,ll,k)=0.d0
      end do
      end do
      do j=1,ny
      do i=1,nx
        ij = (j-1)*nx+i
        beta(ij,ij,k) = cof(i,j,k,7)
      end do
      end do
      do i=2,nx
      do j=1,ny
        ij = (j-1)*nx+i
        beta(ij,ij-1,k) = cof(i,j,k,1)
      end do
      end do
      do i=1,nx-1
      do j=1,ny
        ij = (j-1)*nx+i
        beta(ij,ij+1,k) = cof(i,j,k,2)
      end do
      end do
      do j=2,ny
      do i=1,nx
        ij = (j-1)*nx+i
        beta(ij,ij-nx,k) = cof(i,j,k,3)
      end do
      end do
      do j=1,ny-1
      do i=1,nx
        ij = (j-1)*nx+i
        beta(ij,ij+nx,k) = cof(i,j,k,4)
      end do
      end do
c
c     adjust for periodic b.c in x and/or y as necessary
c
      if (nxa .eq. 0) then
      do j=1,ny
        ii = (j-1)*nx+1
        jj = (j-1)*nx+nx-1
        beta(ii,jj,k) = cof(1,j,k,1)
        ii = (j-1)*nx+nx
        jj = (j-1)*nx+2
        beta(ii,jj,k) = cof(nx,j,k,2)
      end do
      end if

      if (nyc.eq.0) then
      do i=1,nx
        ii = i
        jj = (ny-2)*nx+i
        beta(ii,jj,k) = cof(i,1,k,3)
        ii = (ny-1)*nx+i
        jj = nx + i
        beta(ii,jj,k) = cof(i,ny,k,4)
      end do
      end if
      return
      end

      subroutine setalf3(nxny,nx,ny,nz,cof,alfa,kcur)
C                                                                               
C     set block subdiagonal matrix for k = kcur
C                                                                               
      implicit none
      integer nxny,nx,ny,nz,kcur
      double precision cof(nx,ny,nz,8),alfa(nxny,nxny,*)
      integer k,j,i,ij,ll
      k = kcur
      do j=1,ny
      do i=1,nx
        ij = (j-1)*nx+i
        do ll=1,nxny
          alfa(ij,ll,k)=0.d0
        end do
        alfa(ij,ij,k)=cof(i,j,k,5)
      end do
      end do
      return
      end
