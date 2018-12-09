c
c     file mud3pn.f
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1999 by UCAR                 .
c  .                                                             .
c  .       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                      MUDPACK version 5.0                    .
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
c ... purpose
c
c     mud3pn.f contains subroutines for planar relaxation in the x or
c     y or z direction.  This file must be loaded with any of the real
c     3-d mudpack solvers except mud3sp.
c
      subroutine planxy(wk)
c
c     planar relaxation on (x,y) planes in z direction
c
      implicit none
      real wk(*)
      integer iparm(16),mgo(4)
      real fparm(6)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,
     +ktzbgn,nxk,nyk,nzk,ngrid,klevel,kcur,kps
      integer nx,ny,nz,k,kb,kset,kl,ixy,icxy,ip,ic,isx,jsy,kst,kfn
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/fmud3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/mud3c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +ktzbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
c
c     set integer parameters for xy plane
c
      iparm(1) = 1
      iparm(2) = nxa
      iparm(3) = nxb
      iparm(4) = nyc
      iparm(5) = nyd
      iparm(6) = ixp
      iparm(7) = jyq
      iparm(8) = max0(klevel+iex-ngrid,1)
      iparm(9) = max0(klevel+jey-ngrid,1)
      iparm(10) = nxk(klevel)
      iparm(11) = nyk(klevel)
      iparm(12) = 1
      iparm(13) = 1
c
c     set relaxation method for 2-d
c
      iparm(14) = meth2

      fparm(1) = xa
      fparm(2) = xb
      fparm(3) = yc
      fparm(4) = yd
      fparm(5) = 0.0
c
c     set 2-d multigrid options
c
      mgo(1) = kcycle
      mgo(2) = iprer
      mgo(3) = ipost
      mgo(4) = intpol
      isx = 0
      if (meth2.eq.1.or.meth2.eq.3) then
	isx = 3
	if (nxa.eq.0) isx = 5
      end if
      jsy = 0
      if (meth2.eq.2.or.meth2.eq.3) then
	jsy = 3
	if (nyc.eq.0) jsy = 5
      end if
      nz = nzk(klevel)
      kst = 2
      if (nze.ne.1) kst = 1
      kfn = nz-1
      if (nzf.ne.1) kfn = nz
      do k=kst,kfn
	kset = k
	ixy = kps
	ip = kpbgn(klevel)
	ic = kcbgn(klevel)
	do kb=1,klevel
	  kl = klevel-kb+1
	  nx = nxk(kl)
	  ny = nyk(kl)
c
c     set 2-d coefficient pointer
c
	  icxy = ixy + (nx+2)*(ny+2)
c
c     transfer coefficients from 3-d to 2-d for current kset
c
	  call trscxy(kset,nx,ny,nz,wk(ic),wk(icxy),wk(ip))
	  ixy = ixy+(6+isx+jsy)*nx*ny + (nx+2)*(ny+2)
	  ip = ip + (nx+2)*(ny+2)*(nz+2)
	  ic = ic + 8*nx*ny*nz
	end do
c
c     set 2-d solution for current kset for passage to mup2
c
	nx = nxk(klevel)
	ny = nyk(klevel)
	ip = kpbgn(klevel)
	ixy = 0
	call setpxy(kset,nx,ny,nz,wk(ip),wk(kps),ixy)
c
c     solve on current z plane with full 2-d multigrid cycling
c
	call mup2(iparm,fparm,wk(kps),mgo)
c
c     reset approx from mup2 in wk(ip)
c
	ixy = 1
	call setpxy(kset,nx,ny,nz,wk(ip),wk(kps),ixy)
      end do
c
c     set periodic virtual boundaries if necessary
c
      if (nxa*nyc*nze.eq.0) then
	nx = nxk(klevel)
	ny = nyk(klevel)
	ip = kpbgn(klevel)
	call per3vb(nx,ny,nz,wk(ip),nxa,nyc,nze)
      end if
      return
      end

      subroutine trscxy(kset,nx,ny,nz,cof,cofxy,phi)
c
c     transfer coefficients from 3-d to 2-d
c     to allow 2-d relaxation
c
      implicit none
      integer kset,nx,ny,nz,i,j,k
      real cof(nx,ny,nz,8),cofxy(nx,ny,6)
      real phi(0:nx+1,0:ny+1,0:nz+1)
      k = kset
      do j=1,ny
	do i=1,nx
	  cofxy(i,j,1) = cof(i,j,k,1)
	  cofxy(i,j,2) = cof(i,j,k,2)
	  cofxy(i,j,3) = cof(i,j,k,3)
	  cofxy(i,j,4) = cof(i,j,k,4)
	  cofxy(i,j,5) = cof(i,j,k,7)
	  cofxy(i,j,6) = cof(i,j,k,8)-(cof(i,j,k,5)*phi(i,j,k-1)+
     +                                 cof(i,j,k,6)*phi(i,j,k+1))
	end do
      end do
      return
      end

      subroutine setpxy(kset,nx,ny,nz,phi,phxy,ixy)
      implicit none
      integer kset,nx,ny,nz,i,j,ixy
      real phi(0:nx+1,0:ny+1,0:nz+1),phxy(0:nx+1,0:ny+1)
      if (ixy.eq.0) then
	do j=0,ny+1
	  do i=0,nx+1
	    phxy(i,j) = phi(i,j,kset)
	  end do
	end do
      else
	do j=0,ny+1
	  do i=0,nx+1
	    phi(i,j,kset) = phxy(i,j)
	  end do
	end do
      end if
      return
      end

      subroutine planxz(wk)
c
c     planar relaxation on (x,z) planes in y direction
c
      implicit none
      real wk(*)
      integer iparm(16),mgo(4)
      real fparm(6)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,
     +ktzbgn,nxk,nyk,nzk,ngrid,klevel,kcur,kps
      integer nx,ny,nz,j,kb,jset,kl,ixz,icxz,ip,ic,isx,ksz,jst,jfn
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/fmud3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/mud3c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +ktzbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
c
c     set integer parameters for xz plane
c
      iparm(1) = 1
      iparm(2) = nxa
      iparm(3) = nxb
      iparm(4) = nze
      iparm(5) = nzf
      iparm(6) = ixp
      iparm(7) = kzr
      iparm(8) = max0(klevel+iex-ngrid,1)
      iparm(9) = max0(klevel+kez-ngrid,1)
      iparm(10) = nxk(klevel)
      iparm(11) = nzk(klevel)
      iparm(12) = 1
      iparm(13) = 1
c
c     set relaxation method for 2-d
c
      iparm(14) = meth2

      fparm(1) = xa
      fparm(2) = xb
      fparm(3) = ze
      fparm(4) = zf
      fparm(5) = 0.0
c
c     set 2-d multigrid options
c
      mgo(1) = kcycle
      mgo(2) = iprer
      mgo(3) = ipost
      mgo(4) = intpol
      isx = 0
      if (meth2.eq.1.or.meth2.eq.3) then
	isx = 3
	if (nxa.eq.0) isx = 5
      end if
      ksz = 0
      if (meth2.eq.2.or.meth2.eq.3) then
	ksz = 3
	if (nze.eq.0) ksz = 5
      end if
      ny = nyk(klevel)
      jst = 2
      if (nyc.ne.1) jst = 1
      jfn = ny-1
      if (nyd.ne.1) jfn = ny
      do j=jst,jfn
	jset = j
	ixz = kps
	ip = kpbgn(klevel)
	ic = kcbgn(klevel)
	do kb=1,klevel
	  kl = klevel-kb+1
	  nx = nxk(kl)
	  nz = nzk(kl)
c
c     set 2-d coefficient pointer
c
	  icxz = ixz + (nx+2)*(nz+2)
c
c     transfer coefficients from 3-d to 2-d for current jset
c
	  call trscxz(jset,nx,ny,nz,wk(ic),wk(icxz),wk(ip))
	  ixz = ixz+(6+isx+ksz)*nx*nz + (nx+2)*(nz+2)
	  ip = ip + (nx+2)*(ny+2)*(nz+2)
	  ic = ic + 8*nx*ny*nz
	end do
c
c     set 2-d solution for current jset for passage to mup2
c
	nx = nxk(klevel)
	nz = nzk(klevel)
	ip = kpbgn(klevel)
	ixz = 0
	call setpxz(jset,nx,ny,nz,wk(ip),wk(kps),ixz)
c
c     solve on current y plane with full 2-d multigrid cycling
c
	call mup2(iparm,fparm,wk(kps),mgo)
	ixz = 1
	call setpxz(jset,nx,ny,nz,wk(ip),wk(kps),ixz)
      end do
c
c     set periodic virtual boundaries if necessary
c
      if (nxa*nyc*nze.eq.0) then
	nx = nxk(klevel)
	nz = nzk(klevel)
	ip = kpbgn(klevel)
	call per3vb(nx,ny,nz,wk(ip),nxa,nyc,nze)
      end if
      return
      end

      subroutine trscxz(jset,nx,ny,nz,cof,cofxz,phi)
c
c     transfer coefficients from 3-d to 2-d
c     to allow 2-d relaxation
c
      implicit none
      integer jset,nx,ny,nz,i,j,k
      real cof(nx,ny,nz,8),cofxz(nx,nz,6)
      real phi(0:nx+1,0:ny+1,0:nz+1)
      j = jset
      do k=1,nz
	do i=1,nx
	  cofxz(i,k,1) = cof(i,j,k,1)
	  cofxz(i,k,2) = cof(i,j,k,2)
	  cofxz(i,k,3) = cof(i,j,k,5)
	  cofxz(i,k,4) = cof(i,j,k,6)
	  cofxz(i,k,5) = cof(i,j,k,7)
	  cofxz(i,k,6) = cof(i,j,k,8)-(cof(i,j,k,3)*phi(i,j-1,k)+
     +                                 cof(i,j,k,4)*phi(i,j+1,k))
	end do
      end do
      return
      end

      subroutine setpxz(jset,nx,ny,nz,phi,phxz,ixz)
      implicit none
      integer jset,nx,ny,nz,i,k,ixz
      real phi(0:nx+1,0:ny+1,0:nz+1),phxz(0:nx+1,0:nz+1)
      if (ixz.eq.0) then
	do k=0,nz+1
	  do i=0,nx+1
	    phxz(i,k) = phi(i,jset,k)
	  end do
	end do
      else
	do k=0,nz+1
	  do i=0,nx+1
	    phi(i,jset,k) = phxz(i,k)
	  end do
	end do
      end if
      return
      end

      subroutine planyz(wk)
c
c     planar relaxation on (y,z) planes in x direction
c
      implicit none
      real wk(*)
      integer iparm(16),mgo(4)
      real fparm(6)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,
     +ktzbgn,nxk,nyk,nzk,ngrid,klevel,kcur,kps
      integer nx,ny,nz,i,kb,iset,kl,iyz,icyz,ip,ic,jsy,ksz,ist,ifn
      common/imud3/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,meth2,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/fmud3/xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/mud3c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +ktzbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
c
c     set integer parameters for yz plane
c
      iparm(1) = 1
      iparm(2) = nyc
      iparm(3) = nyd
      iparm(4) = nze
      iparm(5) = nzf
      iparm(6) = jyq
      iparm(7) = kzr
      iparm(8) = max0(klevel+jey-ngrid,1)
      iparm(9) = max0(klevel+kez-ngrid,1)
      iparm(10) = nyk(klevel)
      iparm(11) = nzk(klevel)
      iparm(12) = 1
      iparm(13) = 1
c
c     set relaxation method for 2-d
c
      iparm(14) = meth2

      fparm(1) = yc
      fparm(2) = yd
      fparm(3) = ze
      fparm(4) = zf
      fparm(5) = 0.0
c
c     set 2-d multigrid options
c
      mgo(1) = kcycle
      mgo(2) = iprer
      mgo(3) = ipost
      mgo(4) = intpol
      jsy = 0
      if (meth2.eq.1.or.meth2.eq.3) then
	jsy = 3
	if (nyc.eq.0) jsy = 5
      end if
      ksz = 0
      if (meth2.eq.2.or.meth2.eq.3) then
	ksz = 3
	if (nze.eq.0) ksz = 5
      end if
      nx = nxk(klevel)
      ist = 2
      if (nxa.ne.1) ist = 1
      ifn = nx-1
      if (nxb.ne.1) ifn = nx
      do i=ist,ifn
	iset = i
	iyz = kps
	ip = kpbgn(klevel)
	ic = kcbgn(klevel)
	do kb=1,klevel
	  kl = klevel-kb+1
	  ny = nyk(kl)
	  nz = nzk(kl)
c
c     set 2-d coefficient pointer
c
	  icyz = iyz + (ny+2)*(nz+2)
c
c     transfer coefficients from 3-d to 2-d for current iset
c
	  call trscyz(iset,nx,ny,nz,wk(ic),wk(icyz),wk(ip))
	  iyz = iyz+(6+jsy+ksz)*ny*nz + (ny+2)*(nz+2)
	  ip = ip + (nx+2)*(ny+2)*(nz+2)
	  ic = ic + 8*nx*ny*nz
	end do
c
c     set 2-d solution for current iset for passage to mup2
c
	ny = nyk(klevel)
	nz = nzk(klevel)
	ip = kpbgn(klevel)
	iyz = 0
	call setpyz(iset,nx,ny,nz,wk(ip),wk(kps),iyz)
c
c     solve on current y plane with full 2-d multigrid cycling
c
	call mup2(iparm,fparm,wk(kps),mgo)
	iyz = 1
	call setpyz(iset,nx,ny,nz,wk(ip),wk(kps),iyz)
      end do
c
c     set periodic virtual boundaries if necessary
c
      if (nxa*nyc*nze.eq.0) then
	ny = nyk(klevel)
	nz = nzk(klevel)
	ip = kpbgn(klevel)
	call per3vb(nx,ny,nz,wk(ip),nxa,nyc,nze)
      end if
      return
      end

      subroutine trscyz(iset,nx,ny,nz,cof,cofyz,phi)
c
c     transfer coefficients from 3-d to 2-d
c     to allow 2-d relaxation
c
      implicit none
      integer iset,nx,ny,nz,i,j,k
      real cof(nx,ny,nz,8),cofyz(ny,nz,6)
      real phi(0:nx+1,0:ny+1,0:nz+1)
      i = iset
      do k=1,nz
	do j=1,ny
	  cofyz(j,k,1) = cof(i,j,k,3)
	  cofyz(j,k,2) = cof(i,j,k,4)
	  cofyz(j,k,3) = cof(i,j,k,5)
	  cofyz(j,k,4) = cof(i,j,k,6)
	  cofyz(j,k,5) = cof(i,j,k,7)
	  cofyz(j,k,6) = cof(i,j,k,8)-(cof(i,j,k,1)*phi(i-1,j,k)+
     +                                 cof(i,j,k,2)*phi(i+1,j,k))
	end do
      end do
      return
      end

      subroutine setpyz(iset,nx,ny,nz,phi,phyz,iyz)
      implicit none
      integer iset,nx,ny,nz,iyz,j,k
      real phi(0:nx+1,0:ny+1,0:nz+1),phyz(0:ny+1,0:nz+1)
      if (iyz.eq.0) then
	do k=0,nz+1
	  do j=0,ny+1
	    phyz(j,k) = phi(iset,j,k)
	  end do
	end do
      else
	do k=0,nz+1
	  do j=0,ny+1
	    phi(iset,j,k) = phyz(j,k)
	  end do
	end do
      end if
      return
      end

      subroutine mup2(iparm,fparm,work,mgopt)
c
c     modification of mud2 for planar with mud3
c     whenever mup2 is called coefficients from discretization
c     have already been set but matrices for line (if flagged
c     have not been set)
c
      implicit none
      integer iparm,mgopt
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real fparm,xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer iw,k,kb,nx,ny,ic,itx,ity
      dimension iparm(16),fparm(6),mgopt(4)
      real work(*)
      common/imup2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fmup2/xa,xb,yc,yd,tolmax,relmax
      common/mup2c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
      intl = 1
      nxa = iparm(2)
      nxb = iparm(3)
      nyc = iparm(4)
      nyd = iparm(5)
      ixp = iparm(6)
      jyq = iparm(7)
      iex = iparm(8)
      jey = iparm(9)
      ngrid = max0(iex,jey)
      nfx = iparm(10)
      nfy = iparm(11)
      iguess = iparm(12)
      maxcy = iparm(13)
      method = iparm(14)
      nwork = iparm(15)
      kcycle = mgopt(1)
      iprer = mgopt(2)
      ipost = mgopt(3)
      intpol = mgopt(4)
      xa = fparm(1)
      xb = fparm(2)
      yc = fparm(3)
      yd = fparm(4)
      tolmax = fparm(5)
      isx = 0
      jsy = 0
      if ((method-1)*(method-3).eq.0) then
	isx = 3
	if (nxa.eq.0) isx = 5
      end if
      if ((method-2)*(method-3).eq.0) then
	jsy = 3
	if (nyc.eq.0) jsy = 5
      end if
      kps = 1
      do k=1,ngrid
c       set subgrid sizes
	  nxk(k) = ixp*2**(max0(k+iex-ngrid,1)-1)+1
	  nyk(k) = jyq*2**(max0(k+jey-ngrid,1)-1)+1
	  nx = nxk(k)
	  ny = nyk(k)
	  kps = kps+(nx+2)*(ny+2)+nx*ny*(6+isx+jsy)
	end do
c
c     set work space pointers and discretize pde at each grid level
c
	iw = 1
	do kb=1,ngrid
	  k = ngrid-kb+1
	  nx = nxk(k)
	  ny = nyk(k)
	  kpbgn(k) = iw
	  kcbgn(k) = kpbgn(k)+(nx+2)*(ny+2)
	  ktxbgn(k) = kcbgn(k)+6*nx*ny
	  ktybgn(k) = ktxbgn(k)+isx*nx*ny
	  iw = ktybgn(k)+jsy*nx*ny
	  ic = kcbgn(k)
	  itx = ktxbgn(k)
	  ity = ktybgn(k)
	  klevel = k
	  call dismp2(nx,ny,work(ic),work(itx),work(ity),work(kps))
	end do
      call mup21(work)
      return
      end

      subroutine mup21(wk)
      implicit none
      real wk(*)
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer iter
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imup2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/mup2c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
c
c     cycle from ngrid level
c
      kcur = ngrid
      do iter=1,maxcy
	itero = iter
	call kcymp2(wk)
      end do
      return
      end

      subroutine kcymp2(wk)
c
c     execute multigrid k cycle from kcur grid level
c     kcycle=1 for v cycles, kcycle=2 for w cycles
c
      implicit none
      real wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ip,ic,ipc,irc,itx,ity,ncx,ncy,l,nrel
      real xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      common/imup2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fmup2/xa,xb,yc,yd,tolmax,relmax
      common/mup2c/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
      integer kount(50)
      klevel = kcur
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
c
c     prerelax at current finest grid level
c
      do l=1,iprer
	call relmp2(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      if (kcur .eq. 1) go to 5
c
c     restrict residual to kcur-1 level
c
      ipc = kpbgn(klevel-1)
      ncx = nxk(klevel-1)
      ncy = nyk(klevel-1)
      irc = kcbgn(klevel-1)+5*ncx*ncy
      call resmp2(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),wk(kps))
c
c    set counter for grid levels to zero
c
      do l = 1,kcur
	kount(l) = 0
      end do
c
c    set new grid level and continue k-cycling
c
      klevel = kcur-1
      nrel = iprer
c
c   kcycle control point
c
   10 continue
c
c      post relax when kcur revisited
c
      if (klevel .eq. kcur) go to 5
c
c   count hit at current level
c
      kount(klevel) = kount(klevel)+1
c
c   relax at current level
c
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
      do l=1,nrel
	call relmp2(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      if (kount(klevel) .eq. kcycle+1) then
c
c     kcycle complete at klevel
c
	ipc = ip
	ip = kpbgn(klevel+1)
	ncx = nxk(klevel)
	ncy = nyk(klevel)
	nx = nxk(klevel+1)
	ny = nyk(klevel+1)
c
c    inject correction to finer grid
c
	call cor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd,
     +            intpol,wk(kps))
c
c    reset counter to zero
c
	kount(klevel) = 0
c
c     ascend to next higher level and set to postrelax there
c
	klevel = klevel+1
	nrel = ipost
	go to 10
      else
	if (klevel .gt. 1) then
c
c    kcycle not complete so descend unless at coarsest grid
c
	  ipc = kpbgn(klevel-1)
	  ncx = nxk(klevel-1)
	  ncy = nyk(klevel-1)
	  irc = kcbgn(klevel-1)+5*ncx*ncy
	  call resmp2(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),
     +                wk(kps))
c
c     prerelax at next coarser level
c
	  klevel = klevel-1
	  nrel = iprer
	  go to 10
	else
c
c    postrelax at coarsest level
c
	  do l=1,ipost
	    call relmp2(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
	  end do
	  ipc = ip
	  ip = kpbgn(2)
	  ncx = nxk(1)
	  ncy = nyk(1)
	  nx = nxk(2)
	  ny = nyk(2)
c
c    inject correction to level 2
c
	call cor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd,
     +            intpol,wk(kps))
c
c     set to postrelax at level 2
c
	  nrel = ipost
	  klevel = 2
	  go to 10
	end if
      end if
    5 continue
c
c     post relax at current finest grid level
c
      nx = nxk(kcur)
      ny = nyk(kcur)
      ip = kpbgn(kcur)
      ic = kcbgn(kcur)
      itx = ktxbgn(kcur)
      ity = ktybgn(kcur)
      do l=1,ipost
	call relmp2(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      return
      end

      subroutine dismp2(nx,ny,cof,tx,ty,sum)
c
c     set tridiagonal matrices for line relaxation if necessary
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      integer nx,ny,i,j,im1,jm1
      real tx(nx,ny,*),ty(ny,nx,*),cof(nx,ny,6),sum(*)
      common/imup2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      if (method.eq.1.or.method.eq.3) then
	if (nxa.ne.0) then
c
c    nonperiodic x line relaxation
c
	  do i=1,nx
	    im1 = max0(i-1,1)
	    do j=1,ny
	      tx(im1,j,1) = cof(i,j,1)
	      tx(i,j,2) = cof(i,j,5)
	      tx(i,j,3) = cof(i,j,2)
	    end do
	  end do
	  call factri(ny,nx,tx(1,1,1),tx(1,1,2),tx(1,1,3))
	else
c
c     periodic x line relaxation
c
	  if (nx .gt. 3) then
c
c     set and factor iff nx > 3
c
	    do i=1,nx-1
	      do j=1,ny
		tx(i,j,1) = cof(i,j,1)
		tx(i,j,2) = cof(i,j,5)
		tx(i,j,3) = cof(i,j,2)
	      end do
	    end do
	    call factrp(ny,nx,tx,tx(1,1,2),tx(1,1,3),tx(1,1,4),
     +                  tx(1,1,5),sum)
	  end if
	end if
      end if

      if (method.eq.2.or.method.eq.3) then
	if (nyc.ne.0) then
c
c     nonperiodic y line relaxation
c
	  do j=1,ny
	    jm1 = max0(j-1,1)
	    do i=1,nx
	      ty(jm1,i,1) = cof(i,j,3)
	      ty(j,i,2) = cof(i,j,5)
	      ty(j,i,3) = cof(i,j,4)
	    end do
	  end do
	  call factri(nx,ny,ty(1,1,1),ty(1,1,2),ty(1,1,3))
	else
c
c      periodic y line relaxation
c
	  if (ny .gt. 3) then
c
c     set and factor iff ny > 3
c
	    do j=1,ny-1
	      do i=1,nx
		ty(j,i,1) = cof(i,j,3)
		ty(j,i,2) = cof(i,j,5)
		ty(j,i,3) = cof(i,j,4)
	      end do
	    end do
	    call factrp(nx,ny,ty,ty(1,1,2),ty(1,1,3),ty(1,1,4),
     +                  ty(1,1,5),sum)
	  end if
	end if
      end if
      return
      end

      subroutine resmp2(nx,ny,phi,ncx,ncy,phic,rhsc,cof,resf)
c
c     restrict residual from fine to coarse mesh using fully weighted
c     residual restriction
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ncx,ncy,i,j,ic,jc
      common/imup2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real rhsc(ncx,ncy),resf(nx,ny)
      real phi(0:nx+1,0:ny+1),phic(0:ncx+1,0:ncy+1)
      real cof(nx,ny,6)
c
c     set phic zero
c
      do jc=0,ncy+1
	do ic=0,ncx+1
	  phic(ic,jc) = 0.0
	end do
      end do
c
c     compute residual on fine mesh in resf
c
!$OMP PARALLEL DO SHARED(resf,cof,phi,nx,ny) PRIVATE(i,j)
      do j=1,ny
	do i=1,nx
	  resf(i,j) =  cof(i,j,6)-(
     +             cof(i,j,1)*phi(i-1,j)+
     +             cof(i,j,2)*phi(i+1,j)+
     +             cof(i,j,3)*phi(i,j-1)+
     +             cof(i,j,4)*phi(i,j+1)+
     +             cof(i,j,5)*phi(i,j))
	end do
      end do
!$OMP END PARALLEL DO
c
c     restrict resf to coarse mesh in rhsc
c
      call res2(nx,ny,resf,ncx,ncy,rhsc,nxa,nxb,nyc,nyd)
      return
      end

      subroutine relmp2(nx,ny,phi,cof,tx,ty,sum)
c
c     relaxation for mud2
c
      implicit none
      integer nx,ny
      real phi(*),cof(*),tx(*),ty(*),sum(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imup2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      if (method.eq.0) then                ! point relaxation
	call relmp2p(nx,ny,phi,cof)
      else if (method.eq.1) then           ! line x relaxation
	call slxmp2(nx,ny,phi,cof,tx,sum)
      else if (method.eq.2) then           ! line y relaxation
	call slymp2(nx,ny,phi,cof,ty,sum)
      else if (method.eq.3) then           ! line x&y relaxation
	call slxmp2(nx,ny,phi,cof,tx,sum)
	call slymp2(nx,ny,phi,cof,ty,sum)
      end if
      return
      end

      subroutine relmp2p(nx,ny,phi,cof)
c
c     Gauss-Seidel red/black point relaxation
c
      implicit none
      integer nx,ny,i,j
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imup2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real phi(0:nx+1,0:ny+1),cof(nx,ny,6)
c
c    periodic adjustment bypass block
c
      if (nxa*nyc.ne.0) then
c
c     relax on red grid points
c
!$OMP PARALLEL DO SHARED(cof,phi,nx,ny) PRIVATE(i,j)
	do i=1,nx,2
	  do j=1,ny,2
	    phi(i,j) = (cof(i,j,6) -
     +                 (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                  cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                  cof(i,j,5)
	  end do
	end do
!$OMP END PARALLEL DO
c
!$OMP PARALLEL DO SHARED(cof,phi,nx,ny) PRIVATE(i,j)
	do i=2,nx,2
	  do j=2,ny,2
	    phi(i,j) = (cof(i,j,6) -
     +                 (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                  cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                  cof(i,j,5)
	  end do
	end do
!$OMP END PARALLEL DO
c
c     relax on black grid points
c
c
!$OMP PARALLEL DO SHARED(cof,phi,nx,ny) PRIVATE(i,j)
	do i=1,nx,2
	  do j=2,ny,2
	    phi(i,j) = (cof(i,j,6) -
     +                 (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                  cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                  cof(i,j,5)
	  end do
	end do
!$OMP END PARALLEL DO
c
!$OMP PARALLEL DO SHARED(cof,phi,nx,ny) PRIVATE(i,j)
	do i=2,nx,2
	  do j=1,ny,2
	    phi(i,j) = (cof(i,j,6) -
     +                 (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                  cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                  cof(i,j,5)
	  end do
	end do
!$OMP END PARALLEL DO
	return
      end if
c
c    set periodic virtual boundaries
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
c
c     relax on red grid points
c
!$OMP PARALLEL DO SHARED(cof,phi,nx,ny) PRIVATE(i,j)
      do i=1,nx,2
	do j=1,ny,2
	  phi(i,j) = (cof(i,j,6) -
     +               (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                cof(i,j,5)
	end do
      end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO SHARED(cof,phi,nx,ny) PRIVATE(i,j)
      do i=2,nx,2
	do j=2,ny,2
	  phi(i,j) = (cof(i,j,6) -
     +               (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                cof(i,j,5)
	end do
      end do
!$OMP END PARALLEL DO
c
c    ensure periodic virtual boundary red points are set
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
c
c     relax on black grid points
c
!$OMP PARALLEL DO SHARED(cof,phi,nx,ny) PRIVATE(i,j)
      do i=1,nx,2
	do j=2,ny,2
	  phi(i,j) = (cof(i,j,6) -
     +               (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                cof(i,j,5)
	end do
      end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO SHARED(cof,phi,nx,ny) PRIVATE(i,j)
      do i=2,nx,2
	do j=1,ny,2
	  phi(i,j) = (cof(i,j,6) -
     +               (cof(i,j,1)*phi(i-1,j)+cof(i,j,2)*phi(i+1,j) +
     +                cof(i,j,3)*phi(i,j-1)+cof(i,j,4)*phi(i,j+1)))/
     +                cof(i,j,5)
	end do
      end do
!$OMP END PARALLEL DO
c
c     final set of periodic virtual boundaries
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
      return
      end

      subroutine slxmp2(nx,ny,phi,cof,tx,sum)
c
c     line relaxation in the x direction (periodic or nonperiodic)
c
      implicit none
      integer nx,ny,i,ib,j
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imup2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real phi(0:nx+1,0:ny+1),cof(nx,ny,6),tx(nx,ny,*),sum(ny)
c
c     replace line x with point gauss-seidel if
c     x direction is periodic and nx = 3 (coarsest)
c
      if (nxa .eq. 0 .and. nx .eq. 3) then
	call relmp2p(nx,ny,phi,cof)
	return
      end if
c
c     set periodic y virtual boundary if necessary
c
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if

      if (nxa.ne.0) then
c
c     x direction not periodic
c
c     sweep odd j lines
c
!$OMP PARALLEL DO SHARED(cof,phi,tx,nx,ny) PRIVATE(i,ib,j)
	do j=1,ny,2
	  do i=1,nx
	    phi(i,j)=cof(i,j,6)-cof(i,j,3)*phi(i,j-1)-cof(i,j,4)*
     +               phi(i,j+1)
	  end do
c
c     forward sweep
c
	  do i=2,nx
	    phi(i,j) = phi(i,j)-tx(i-1,j,1)*phi(i-1,j)
	  end do
c
c     backward sweep
c
	  phi(nx,j) = phi(nx,j)/tx(nx,j,2)
	  do ib=2,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j))/tx(i,j,2)
	  end do
	end do
!$OMP END PARALLEL DO
c
c     sweep even j lines forward and back
c
!$OMP PARALLEL DO SHARED(cof,phi,tx,nx,ny) PRIVATE(i,j,ib)
	do j=2,ny,2
	  do i=1,nx
	    phi(i,j)=cof(i,j,6)-cof(i,j,3)*phi(i,j-1)-cof(i,j,4)*
     +               phi(i,j+1)
	  end do
	  do i=2,nx
	    phi(i,j) = phi(i,j)-tx(i-1,j,1)*phi(i-1,j)
	  end do
	  phi(nx,j) = phi(nx,j)/tx(nx,j,2)
	  do ib=2,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j))/tx(i,j,2)
	  end do
	end do
!$OMP END PARALLEL DO
      else
c
c     x direction periodic
c
	do j=1,ny
	  sum(j) = 0.0
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
c
c      sweep odd lines forward and back
c
!$OMP PARALLEL DO SHARED(sum,cof,phi,tx,nx,ny) PRIVATE(i,j,ib)
	do j=1,ny,2
	  do i=1,nx-1
	    phi(i,j)=cof(i,j,6)-cof(i,j,3)*phi(i,j-1)-cof(i,j,4)*
     +               phi(i,j+1)
	  end do
c
c     forward sweep
	  do i=2,nx-2
	    phi(i,j) = phi(i,j)-tx(i,j,1)*phi(i-1,j)
	  end do
	  do i=1,nx-2
	    sum(j) = sum(j)+tx(i,j,5)*phi(i,j)
	  end do
	  phi(nx-1,j) = phi(nx-1,j)-sum(j)
c
c     backward sweep
c
	  phi(nx-1,j) = phi(nx-1,j)/tx(nx-1,j,2)
	  phi(nx-2,j) = (phi(nx-2,j)-tx(nx-2,j,4)*phi(nx-1,j))/
     +                   tx(nx-2,j,2)
	  do ib=4,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j)-tx(i,j,4)*
     +                 phi(nx-1,j))/tx(i,j,2)
	  end do
	end do
!$OMP END PARALLEL DO
c
c     set periodic and virtual points for j odd
c
	do j=1,ny,2
	  phi(nx,j) = phi(1,j)
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
c
c     sweep even j lines
c
!$OMP PARALLEL DO SHARED(sum,cof,phi,tx,nx,ny) PRIVATE(i,j,ib)
	do j=2,ny,2
	  do i=1,nx-1
	    phi(i,j)=cof(i,j,6)-cof(i,j,3)*phi(i,j-1)-cof(i,j,4)*
     +               phi(i,j+1)
	  end do
c
c     forward sweep
c
	  do i=2,nx-2
	    phi(i,j) = phi(i,j)-tx(i,j,1)*phi(i-1,j)
	  end do
	  do i=1,nx-2
	    sum(j) = sum(j)+tx(i,j,5)*phi(i,j)
	  end do
	  phi(nx-1,j) = phi(nx-1,j)-sum(j)
c
c     backward sweep
c
	  phi(nx-1,j) = phi(nx-1,j)/tx(nx-1,j,2)
	  phi(nx-2,j) = (phi(nx-2,j)-tx(nx-2,j,4)*phi(nx-1,j))/
     +                   tx(nx-2,j,2)
	  do ib=4,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j)-tx(i,j,4)*
     +                 phi(nx-1,j))/tx(i,j,2)
	  end do
	end do
!$OMP END PARALLEL DO
c
c     set periodic and virtual points for j even
c
	do j=2,ny,2
	  phi(nx,j) = phi(1,j)
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
c
c     set periodic y virtual boundaries if necessary
c
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
      return
      end

      subroutine slymp2(nx,ny,phi,cof,ty,sum)
      implicit none
      integer nx,ny,i,j,jb
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imup2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real phi(0:nx+1,0:ny+1),cof(nx,ny,6),ty(ny,nx,*),sum(nx)
c
c     replace line y with point gauss-seidel if
c     y direction is periodic and ny = 3
c
      if (nyc .eq. 0 .and. ny .eq. 3) then
	call relmp2p(nx,ny,phi,cof)
	return
      end if
c
c      set periodic and virtual x boundaries if necessary
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx,j) = phi(1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if

      if (nyc.ne.0) then
c
c     y direction not periodic
c
c
c     sweep odd x lines
c
!$OMP PARALLEL DO SHARED(cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
	do i=1,nx,2
	  do j=1,ny
	    phi(i,j)=cof(i,j,6)-cof(i,j,1)*phi(i-1,j)-cof(i,j,2)*
     +               phi(i+1,j)
	  end do
c
c     forward sweep thru odd x lines
c
	  do j=2,ny
	    phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
	  end do
c
c      backward sweep
c
	  phi(i,ny) = phi(i,ny)/ty(ny,i,2)
	  do jb=2,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1))/ty(j,i,2)
	  end do
	end do
!$OMP END PARALLEL DO
c
c     forward sweep even x lines
c
!$OMP PARALLEL DO SHARED(cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
	do i=2,nx,2
	  do j=1,ny
	    phi(i,j)=cof(i,j,6)-cof(i,j,1)*phi(i-1,j)-cof(i,j,2)*
     +               phi(i+1,j)
	  end do
	  do j=2,ny
	    phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
	  end do
c
c      backward sweep
c
	  phi(i,ny) = phi(i,ny)/ty(ny,i,2)
	  do jb=2,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1))/ty(j,i,2)
	  end do
	end do
!$OMP END PARALLEL DO
      else
c
c     y direction periodic
c
	do i=1,nx
	  sum(i) = 0.0
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny) = phi(i,1)
	  phi(i,ny+1) = phi(i,2)
	end do
c
!$OMP PARALLEL DO SHARED(cof,phi,ty,sum,nx,ny) PRIVATE(i,j,jb)
	do i=1,nx,2
	  do j=1,ny-1
	    phi(i,j)=cof(i,j,6)-cof(i,j,1)*phi(i-1,j)-cof(i,j,2)*
     +               phi(i+1,j)
	  end do
	  do j=2,ny-2
	    phi(i,j) = phi(i,j)-ty(j,i,1)*phi(i,j-1)
	  end do
	  do j=1,ny-2
	    sum(i) = sum(i)+ty(j,i,5)*phi(i,j)
	  end do
	  phi(i,ny-1) = phi(i,ny-1)-sum(i)
c
c     backward sweep
c
	  phi(i,ny-1) = phi(i,ny-1)/ty(ny-1,i,2)
	  phi(i,ny-2) = (phi(i,ny-2)-ty(ny-2,i,4)*phi(i,ny-1))/
     +                   ty(ny-2,i,2)
	  do jb=4,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1)-ty(j,i,4)*
     +                  phi(i,ny-1))/ty(j,i,2)
	  end do
	end do
!$OMP END PARALLEL DO
c
c       set odd periodic and virtual y boundaries
c
	do i=1,nx,2
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny) = phi(i,1)
	  phi(i,ny+1) = phi(i,2)
	end do
c
c     forward sweep even x lines
c
c
!$OMP PARALLEL DO SHARED(sum,cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
	do i=2,nx,2
	  do j=1,ny-1
	    phi(i,j)=cof(i,j,6)-cof(i,j,1)*phi(i-1,j)-cof(i,j,2)*
     +               phi(i+1,j)

	  end do
	  do j=2,ny-2
	    phi(i,j) = phi(i,j)-ty(j,i,1)*phi(i,j-1)
	  end do
	  do j=1,ny-2
	    sum(i) = sum(i)+ty(j,i,5)*phi(i,j)
	  end do
	  phi(i,ny-1) = phi(i,ny-1)-sum(i)
c
c     backward sweep
c
	  phi(i,ny-1) = phi(i,ny-1)/ty(ny-1,i,2)
	  phi(i,ny-2) = (phi(i,ny-2)-ty(ny-2,i,4)*phi(i,ny-1))/
     +                   ty(ny-2,i,2)
	  do jb=4,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1)-ty(j,i,4)*
     +                  phi(i,ny-1))/ty(j,i,2)
	  end do
	end do
!$OMP END PARALLEL DO
c
c       set even periodic and virtual y boundaries
c
	do i=2,nx,2
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny) = phi(i,1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
c
c      set periodic and virtual x boundaries if necessary
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      return
      end

