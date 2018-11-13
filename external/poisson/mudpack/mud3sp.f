c
c     file mud3sp.f
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
c ... purpose (see mud3sp.d for details)
c
c     mud3sp attempts to produce a second order finite difference
c     approximation to the three dimensional SEPARABLE elliptic
c     partial differential equation of the form:
c
c       cxx(x)*pxx + cx(x)*px + cex(x)*p(x,y,z) +
c
c       cyy(y)*pyy + cy(y)*py + cey(y)*p(x,y,z) +
c
c       czz(z)*pzz + cz(z)*pz + cez(z)*p(x,y,z) = r(x,y,z)
c
c     SEPARABILITY means:
c
c       cxx,cx,cex depend only on x
c       cyy,cy,cey depend only on y
c       czz,cz,cez depend only on z
c
c     For example, LaPlace's equation in Cartesian coordinates is separable.
c     Nonseparable elliptic PDEs can be approximated with muh3 or mud3.
c
c
c ... documentation and test files
c
c     see the documentation file "mud3sp.d" for a complete discussion
c     of how to use subroutine mud3sp.  file "tmud3sp.f" is a test/driver
c     sample program illustrating use of mud3sp
c
c ... required MUDPACK files
c
c     mudcom.f
c
      subroutine mud3sp(iparm,fparm,work,cfx,cfy,cfz,bndyc,rhs,phi,
     +                  mgopt,ierror)
      implicit none
      integer iparm(22),mgopt(4),ierror
      real work(*),phi(*),rhs(*),fparm(8)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3sp/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/fmud3sp/xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer kpbgn,krbgn,kcxbgn,kcybgn,kczbgn,nxk,nyk,nzk,ngrid,
     +        klevel,kcur,kps
      common/mud3spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +kczbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
      integer int,nx,ny,nz,k,kb,icx,icy,icz
      external cfx,cfy,cfz,bndyc
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
      nwork = iparm(20)
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
      kcycle = mgopt(1)
      if (kcycle .eq. 0) then
c
c     use default settings
c
	kcycle = 2
	iprer = 2
	ipost = 1
	intpol = 3
      else
	iprer = mgopt(2)
	ipost = mgopt(3)
	intpol = mgopt(4)
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
	if (method.ne.0) return
	ierror = 9
c
c     set subgrid sizes and pointers
c
	kps = 1
	do kb=1,ngrid
	  k = ngrid-kb+1
	  nxk(k) = ixp*2**(max0(k+iex-ngrid,1)-1)+1
	  nyk(k) = jyq*2**(max0(k+jey-ngrid,1)-1)+1
	  nzk(k) = kzr*2**(max0(k+kez-ngrid,1)-1)+1
	  nx = nxk(k)
	  ny = nyk(k)
	  nz = nzk(k)
	  kpbgn(k) = kps
	  krbgn(k) = kpbgn(k)+(nx+2)*(ny+2)*(nz+2)
	  kcxbgn(k) = krbgn(k)+nx*ny*nz
	  kcybgn(k) = kcxbgn(k) + 3*nx
	  kczbgn(k) = kcybgn(k) + 3*ny
	  kps = kczbgn(k) + 3*nz
	end do
c
c     set and check minimal work space
c
	nx = nxk(ngrid)
	ny = nyk(ngrid)
	nz = nzk(ngrid)
	iparm(21)=kps+(nx+2)*(ny+2)*(nz+2)
	lwork = iparm(21)
	if (lwork .gt. nwork) return
	ierror = 10   ! check solution region
	if (xb.le.xa .or. yd.le.yc .or. zf.le.ze) return
	ierror = 11
	if (tolmax .lt. 0.0) return
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
	  icx = kcxbgn(k)
	  icy = kcybgn(k)
	  icz = kczbgn(k)
	  call dismd3sp(nx,ny,nz,work(icx),work(icy),work(icz),
     +    bndyc,cfx,cfy,cfz,ierror)
	end do
	return
      end if   ! end of intl=0 initialization call block
      nx = nfx
      ny = nfy
      nz = nfz
      call mud3sp1(nx,ny,nz,rhs,phi,cfx,cfy,cfz,bndyc,work)
      iparm(22) = itero
      if (ierror .le. 0) then
	if (tolmax.gt.0.0) then
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

      subroutine mud3sp1(nx,ny,nz,rhsf,phif,cfx,cfy,cfz,bndyc,wk)
      implicit none
      integer nx,ny,nz
      real rhsf(nx,ny,nz),phif(nx,ny,nz),wk(*)
      common/imud3sp/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer kpbgn,krbgn,kcxbgn,kcybgn,kczbgn,nxk,nyk,nzk,ngrid,
     +        klevel,kcur,kps
      common/mud3spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +kczbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
      integer ip,ir,irc,ncx,ncy,ncz,ipc
      integer i,j,k,kb,jk,kk,ijk,iter
      real phmax
      common/fmud3sp/xa,xb,yc,yd,ze,zf,tolmax,relmax
      external cfx,cfy,cfz,bndyc
      nx = nxk(ngrid)
      ny = nyk(ngrid)
      nz = nzk(ngrid)
      ip = kpbgn(ngrid)
      ir = krbgn(ngrid)
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
	  ir = krbgn(k+1)
	  ncx = nxk(k)
	  ncy = nyk(k)
	  ncz = nzk(k)
	  ipc = kpbgn(k)
	  irc = krbgn(k)
c
c     transfer down to all grid levels
c
	  call trsfc3(nx,ny,nz,wk(ip),wk(ir),ncx,ncy,ncz,
     +                wk(ipc),wk(irc))
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
	  ir = krbgn(k)
	  call adjmd3sp(nx,ny,nz,wk(ip),wk(ir),bndyc,cfx,cfy,cfz)
	end do
c
c     execute one full multigrid cycle
c
	do k=1,ngrid-1
	  kcur = k
	  call kcymd3sp(wk)
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
	ir = krbgn(ngrid)
	call adjmd3sp(nx,ny,nz,wk(ip),wk(ir),bndyc,cfx,cfy,cfz)
      end if
c
c     execute maxcy more multigrid k cycles from finest level
c
      kcur = ngrid
      do iter=1,maxcy
	itero = iter
	call kcymd3sp(wk)
	if (tolmax.gt.0.0) then
c
c      error control
c
	  relmax = 0.0
	  phmax = 0.0
	  do k=1,nfz
	    kk = k*(nfx+2)*(nfy+2)
	    do j=1,nfy
	      jk = kk+j*(nfx+2)
	      do i=1,nfx
		ijk = jk+i+1
		phmax = amax1(phmax,abs(wk(ijk)))
		relmax = amax1(relmax,abs(wk(ijk)-phif(i,j,k)))
		phif(i,j,k) = wk(ijk)
	      end do
	    end do
	  end do
c
c     set maximum relative difference and check for convergence
c
	  if (phmax.gt.0.0) relmax = relmax/phmax
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

      subroutine kcymd3sp(wk)
c
c     perform multigrid k-cycle at kcur level
c     kcycle = 1 corresponds to v cycles
c     kcycle = 2 corresponds to w cycles
c
      implicit none
      real wk(*)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3sp/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      integer kpbgn,krbgn,kcxbgn,kcybgn,kczbgn,nxk,nyk,nzk,ngrid,
     +        klevel,kcur,kps
      common/mud3spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +kczbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
      integer nx,ny,nz,ncx,ncy,ncz
      integer kount(50),ip,ir,icx,icy,icz,ipc,irc,nrel,l
      klevel = kcur
c
c     pre-relax at current finest grid level
c
      do l = 1,iprer
	call relmd3sp(wk)
      end do
c
c     if at coarsest level post-relax
c
      if (kcur .eq. 1) go to 2
c
c     restrict residual to kcur-1
c
      nx = nxk(klevel)
      ny = nyk(klevel)
      nz = nzk(klevel)
      ip = kpbgn(klevel)
      ir = krbgn(klevel)
      icx = kcxbgn(klevel)
      icy = kcybgn(klevel)
      icz = kczbgn(klevel)
      ipc = kpbgn(klevel-1)
      ncx = nxk(klevel-1)
      ncy = nyk(klevel-1)
      ncz = nzk(klevel-1)
      irc = krbgn(klevel-1)
c
c     use full weighting with residual restriction
c
      call resmd3sp(nx,ny,nz,wk(ip),wk(ir),wk(icx),wk(icy),wk(icz),
     +              ncx,ncy,ncz,wk(ipc),wk(irc),wk(kps))
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
c     relax at current level
c
      do l = 1,nrel
	call relmd3sp(wk)
      end do
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
	  ir = krbgn(klevel)
	  icx = kcxbgn(klevel)
	  icy = kcybgn(klevel)
	  icz = kczbgn(klevel)
	  ncx = nxk(klevel-1)
	  ncy = nyk(klevel-1)
	  ncz = nzk(klevel-1)
	  irc = krbgn(klevel-1)
	  ipc = kpbgn(klevel-1)
	  call resmd3sp(nx,ny,nz,wk(ip),wk(ir),wk(icx),wk(icy),wk(icz),
     +                  ncx,ncy,ncz,wk(ipc),wk(irc),wk(kps))
c
c     pre-relax at next coarser level
c
	  klevel = klevel-1
	  nrel = iprer
	  go to 1
	else
c
c     post-relax at coarsest level (klevel=1)
c
	  do l = 1,ipost
	    call relmd3sp(wk)
	  end do
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
	call relmd3sp(wk)
      end do
      return
      end

      subroutine resmd3sp(nx,ny,nz,phi,rhs,cofx,cofy,cofz,ncx,ncy,ncz,
     +                    phic,rhsc,resf)
c
c     compute fully weighted residual restriction in rhsc
c
      implicit none
      integer nx,ny,nz,ncx,ncy,ncz
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3sp/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      real phi(0:nx+1,0:ny+1,0:nz+1),phic(0:ncx+1,0:ncy+1,0:ncz+1)
      real rhs(nx,ny,nz),rhsc(ncx,ncy,ncz),resf(nx,ny,nz)
      real cofx(nx,3),cofy(ny,3),cofz(nz,3)
      integer ic,jc,kc,i,j,k,ist,ifn,jst,jfn,kst,kfn
c
c     initialize phic to zero
c
      do kc=0,ncz+1
	do jc=0,ncy+1
	  do ic=0,ncx+1
	    phic(ic,jc,kc) = 0.0
	  end do
	end do
      end do
c
c     intialize residual to zero and set limits
c
      do k=1,nz
	do j=1,ny
	  do i=1,nx
	    resf(i,j,k) = 0.0
	  end do
	end do
      end do
c
c     set loop limits
c
      ist = 1
      if (nxa.eq.1) ist = 2
      ifn = nx
      if (nxb.eq.1) ifn = nx-1
      jst = 1
      if (nyc.eq.1) jst = 2
      jfn = ny
      if (nyd.eq.1) jfn = ny-1
      kst = 1
      if (nze.eq.1) kst = 2
      kfn = nz
      if (nzf.eq.1) kfn = nz-1
c
c     compute fine grid residual
c
!$OMP PARALLEL DO PRIVATE(i,j,k), SHARED(phi,nx,ny,nz)
!$OMP+SHARED(ist,ifn,jst,jfn,kst,kfn,cofx,cofy,cofz,rhs)
      do k=kst,kfn
	do j=jst,jfn
	  do i=ist,ifn
	    resf(i,j,k) =  rhs(i,j,k)-(
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1) +
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))*phi(i,j,k)  )
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
      end

      subroutine dismd3sp(nx,ny,nz,cofx,cofy,cofz,bndyc,cfx,cfy,cfz,ier)
c
c     discretize the 3-d elliptic pde
c
      implicit none
      integer nx,ny,nz,ier
      real cofx(nx,3),cofy(ny,3),cofz(nz,3)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer kpbgn,krbgn,kcxbgn,kcybgn,kczbgn,nxk,nyk,nzk,ngrid,
     +        klevel,kcur,kps
      common/mud3spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +kczbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
      common/imud3sp/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/fmud3sp/xa,xb,yc,yd,ze,zf,tolmax,relmax
      real dlx,dly,dlz,dlx2,dly2,dlz2,dlxx,dlyy,dlzz,cmin,cemax,alfmax
      real cxx,cyy,czz,cx,cy,cz,cex,cey,cez,alfa,gbdy,x,y,z,c1,c2
      integer i,j,k,ist,ifn,jst,jfn,kst,kfn,kbdy
      external bndyc,cfx,cfy,cfz
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
      cemax = 0.0
c
c     set x,y,z subscript limits to bypass specified boundaries
c     when calling cfx,cfy,cfz or bndyc
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
      do i=ist,ifn
	x = xa+(i-1)*dlx
	call cfx(x,cxx,cx,cex)
	cmin = amin1(cmin,cxx)
	cemax = amax1(abs(cex),cemax)
c
c     check if pde is "hyperbolic" at finest grid level
c
	if (klevel.eq.ngrid) then
	  if ((abs(cx)*dlx .gt. abs(cxx+cxx))) ier = -4
	end if
c
c     adjust second order coefficients so that pde is not "hyperbolic"
c     this is especially possible on coarser grids if there are non-zero
c     first order terms
c
	cxx = amax1(cxx,abs(cx)*dlx*0.5)
	cofx(i,1) = cxx/dlxx-cx/dlx2
	cofx(i,2) = cxx/dlxx+cx/dlx2
	cofx(i,3) = cex-(cofx(i,1)+cofx(i,2))
      end do
      do j=jst,jfn
	y = yc+(j-1)*dly
	call cfy(y,cyy,cy,cey)
	cmin = amin1(cmin,cyy)
	cemax = amax1(abs(cey),cemax)
c
c     check if pde is "hyperbolic" at finest grid level
c
	if (klevel.eq.ngrid) then
	  if ((abs(cy)*dly .gt. abs(cyy+cyy))) ier = -4
	end if
c
c     adjust second order coefficients so that pde is not "hyperbolic"
c     this is especially possible on coarser grids if there are non-zero
c     first order terms
c
	cyy = amax1(cyy,abs(cy)*dly*0.5)
	cofy(j,1) = cyy/dlyy-cy/dly2
	cofy(j,2) = cyy/dlyy+cy/dly2
	cofy(j,3) = cey-(cofy(j,1)+cofy(j,2))
      end do
      do k=kst,kfn
	z = ze+(k-1)*dlz
	call cfz(z,czz,cz,cez)
	cmin = amin1(cmin,czz)
	cemax = amax1(abs(cez),cemax)
c
c     check if pde is "hyperbolic" at finest grid level
c
	if (klevel.eq.ngrid) then
	  if ((abs(cz)*dlz .gt. abs(czz+czz))) ier = -4
	end if
c
c     adjust second order coefficients so that pde is not "hyperbolic"
c     this is especially possible on coarser grids if there are non-zero
c     first order terms
c
	czz = amax1(czz,abs(cz)*dlz*0.5)
	cofz(k,1) = czz/dlzz-cz/dlz2
	cofz(k,2) = czz/dlzz+cz/dlz2
	cofz(k,3) = cez-(cofz(k,1)+cofz(k,2))
      end do
c
c     set nonfatal error flag if ellipticity test fails
c
      if (cmin.le.0.0) ier = -2
      alfmax = 0.0
c
c     adjust equation at mixed b.c.
c
      if (nxa.eq.2) then
	kbdy = 1
	i = 1
	c1 = cofx(i,1)
	cofx(i,1) = 0.0
	cofx(i,2) = cofx(i,2)+c1
	y = yc+dly
	z = ze+dlz
c
c     compute constant coefficient alfa
c
	call bndyc(kbdy,y,z,alfa,gbdy)
	alfmax = amax1(alfmax,abs(alfa))
	cofx(i,3) = cofx(i,3)+dlx2*alfa*c1
      end if
      if (nxb.eq.2) then
	kbdy = 2
	i = nx
	y = yc+dly
	z = ze+dlz
c
c     compute constant coefficient alfa
c
	call bndyc(kbdy,y,z,alfa,gbdy)
	c2 = cofx(i,2)
	cofx(i,1) = cofx(i,1)+c2
	cofx(i,2) = 0.0
	cofx(i,3) = cofx(i,3)-dlx2*alfa*c2
	alfmax = amax1(abs(alfa),alfmax)
      end if
      if (nyc.eq.2) then
	kbdy = 3
	j = 1
	x = xa+dlx
	z = ze+dlz
c
c     compute constant coefficient alfa
c
	call bndyc(kbdy,x,z,alfa,gbdy)
	c1 = cofy(j,1)
	cofy(j,1) = 0.0
	cofy(j,2) = cofy(j,2) + c1
	cofy(j,3) = cofy(j,3) + dly2*alfa*c1
	alfmax = amax1(abs(alfa),alfmax)
      end if
      if (nyd.eq.2) then
	kbdy = 4
	j = ny
	x = xa+dlx
	z = ze+dlz
c
c     compute constant coefficient alfa
c
	call bndyc(kbdy,x,z,alfa,gbdy)
	c2 = cofy(j,2)
	cofy(j,2) = 0.0
	cofy(j,1) = cofy(j,1) + c2
	cofy(j,3) = cofy(j,3) - dly2*c2*alfa
	alfmax = amax1(abs(alfa),alfmax)
      end if
      if (nze.eq.2) then
	kbdy = 5
	k = 1
	x = xa+dlx
	y = yc+dly
c
c     compute constant coefficient alfa
c
	call bndyc(kbdy,x,y,alfa,gbdy)
	c1 = cofz(k,1)
	cofz(k,1) = 0.0
	cofz(k,2) = cofz(k,2) + c1
	cofz(k,3) = cofz(k,3) + dlz2*alfa*c1
	alfmax = amax1(abs(alfa),alfmax)
      end if
      if (nzf.eq.2) then
	kbdy = 6
	k = nz
	x = xa+dlx
	y = yc+dly
c
c     compute constant coefficient alfa
c
	call bndyc(kbdy,x,y,alfa,gbdy)
	c2 = cofz(k,2)
	cofz(k,2) = 0.0
	cofz(k,1) = cofz(k,1) + c2
	cofz(k,3) = cofz(k,3) - dlz2*alfa*c2
	alfmax = amax1(abs(alfa),alfmax)
      end if
c
c     flag continuous singular elliptic pde if detected
c
      if (cemax.eq.0.0.and.alfmax.eq.0.0) then
	if (nxa.eq.0.or.(nxa.eq.2.and.nxb.eq.2)) then
	  if (nyc.eq.0.or.(nyc.eq.2.and.nyd.eq.2)) then
	    if (nze.eq.0.or.(nze.eq.2.and.nzf.eq.2)) then
	      ier = -3
	    end if
	  end if
	end if
      end if
      return
      end

      subroutine adjmd3sp(nx,ny,nz,phi,rhs,bndyc,cfx,cfy,cfz)
c
c     adjust rhs for solution in cof(i,j,k,8) on non-initial calls
c     (i.e., values in cof have not changed)
c
      implicit none
      integer nx,ny,nz
      real phi(0:nx+1,0:ny+1,0:nz+1),rhs(nx,ny,nz)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer kpbgn,krbgn,kcxbgn,kcybgn,kczbgn,nxk,nyk,nzk,ngrid,
     +        klevel,kcur,kps
      common/mud3spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +kczbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
      common/imud3sp/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/fmud3sp/xa,xb,yc,yd,ze,zf,tolmax,relmax
      real dlx,dly,dlz,dlx2,dly2,dlz2,dlxx,dlyy,dlzz
      real cxx,cyy,czz,cx,cy,cz,cex,cey,cez,alfa,gbdy,x,y,z,c1,c2
      integer i,j,k,ist,ifn,jst,jfn,kst,kfn,kbdy
      external bndyc,cfx,cfy,cfz
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
c     set x,y,z subscript limits for calls to cfx,cfy,cfz,bndyc
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
	call cfx(x,cxx,cx,cex)
	cxx = amax1(cxx,abs(cx)*dlx*0.5)
	c1 = cxx/dlxx-cx/dlx2
	do k=kst,kfn
	  z=ze+(k-1)*dlz
	  do j=jst,jfn
	    y=yc+(j-1)*dly
	    call bndyc(kbdy,y,z,alfa,gbdy)
	    rhs(i,j,k) = rhs(i,j,k)+dlx2*c1*gbdy
	  end do
	end do
      end if
      if (nxb.eq.2) then
	kbdy=2
	x=xb
	i = nx
	call cfx(x,cxx,cx,cex)
	cxx = amax1(cxx,abs(cx)*dlx*0.5)
	c2 = cxx/dlxx+cx/dlx2
	do k=kst,kfn
	  z=ze+(k-1)*dlz
	  do j=jst,jfn
	    y=yc+(j-1)*dly
	    call bndyc(kbdy,y,z,alfa,gbdy)
	    rhs(i,j,k) = rhs(i,j,k)-dlx2*c2*gbdy
	  end do
	end do
      end if
      if (nyc.eq.2) then
	kbdy = 3
	y=yc
	j = 1
	call cfy(y,cyy,cy,cey)
	cyy = amax1(cyy,abs(cy)*dly*0.5)
	c1 = cyy/dlyy-cy/dly2
	do k=kst,kfn
	  z=ze+(k-1)*dlz
	  do i=ist,ifn
	    x=xa+(i-1)*dlx
	    call bndyc(kbdy,x,z,alfa,gbdy)
	    rhs(i,j,k) = rhs(i,j,k)+dly2*c1*gbdy
	  end do
	end do
      end if
      if (nyd.eq.2) then
	kbdy = 4
	y=yd
	j = ny
	call cfy(y,cyy,cy,cey)
	cyy = amax1(cyy,abs(cy)*dly*0.5)
	c2 = cyy/dlyy+cy/dly2
	do k=kst,kfn
	  z=ze+(k-1)*dlz
	  do i=ist,ifn
	    x=xa+(i-1)*dlx
	    call bndyc(kbdy,x,z,alfa,gbdy)
	    rhs(i,j,k) = rhs(i,j,k)-dly2*c2*gbdy
	  end do
	end do
      end if
      if (nze.eq.2) then
	kbdy = 5
	k = 1
	z=ze
	call cfz(z,czz,cz,cez)
	czz = amax1(czz,abs(cz)*dlz*0.5)
	c1 = czz/dlzz-cz/dlz2
	do j=jst,jfn
	  y=yc+(j-1)*dly
	  do i=ist,ifn
	    x=xa+(i-1)*dlx
	    call bndyc(kbdy,x,y,alfa,gbdy)
	    rhs(i,j,k) = rhs(i,j,k)+dlz2*c1*gbdy
	  end do
	end do
      end if
      if (nzf.eq.2) then
	kbdy = 6
	z=zf
	k = nz
	call cfz(z,czz,cz,cez)
	czz = amax1(czz,abs(cz)*dlz*0.5)
	c2 = czz/dlzz+cz/dlz2
	do j=jst,jfn
	  y=yc+(j-1)*dly
	  do i=ist,ifn
	    x=xa+(i-1)*dlx
	    call bndyc(kbdy,x,y,alfa,gbdy)
	    rhs(i,j,k) = rhs(i,j,k)-dlz2*c2*gbdy
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
	    rhs(i,j,k) = phi(i,j,k)
	  end do
	end do
      end if
      if (nxb.eq.1) then
	i = nx
	do j=1,ny
	  do k=1,nz
	    rhs(i,j,k) = phi(i,j,k)
	  end do
	end do
      end if
      if (nyc.eq.1) then
	j = 1
	do k=1,nz
	  do i=1,nx
	    rhs(i,j,k) = phi(i,j,k)
	  end do
	end do
      end if
      if (nyd.eq.1) then
	j = ny
	do k=1,nz
	  do i=1,nx
	    rhs(i,j,k) = phi(i,j,k)
	  end do
	end do
      end if
      if (nze.eq.1) then
	k = 1
	do j=1,ny
	  do i=1,nx
	    rhs(i,j,k) = phi(i,j,k)
	  end do
	end do
      end if
      if (nzf.eq.1) then
	k = nz
	do j=1,ny
	  do i=1,nx
	    rhs(i,j,k) = phi(i,j,k)
	  end do
	end do
      end if
      return
      end

      subroutine relmd3sp(wk)
c
c     use point or line relaxation in the x and/or y and/or z
c     or planar relaxation in the x,y or x,z or y,z planes
c
      implicit none
      real wk(*)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer kpbgn,krbgn,kcxbgn,kcybgn,kczbgn,nxk,nyk,nzk,ngrid,
     +        klevel,kcur,kps
      common/mud3spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +kczbgn(50),nxk(50),nyk(50),nzk(50),ngrid,klevel,kcur,kps
      common/imud3sp/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/fmud3sp/xa,xb,yc,yd,ze,zf,tolmax,relmax
      integer nx,ny,nz,ip,ir,icx,icy,icz
      nx = nxk(klevel)
      ny = nyk(klevel)
      nz = nzk(klevel)
      ip = kpbgn(klevel)
      ir = krbgn(klevel)
      icx = kcxbgn(klevel)
      icy = kcybgn(klevel)
      icz = kczbgn(klevel)
c
c     gauss-seidel pointwise red/black relaxation
c
      call relmd3spp(nx,ny,nz,wk(ip),wk(ir),wk(icx),wk(icy),wk(icz))
      return
      end

      subroutine relmd3spp(nx,ny,nz,phi,rhs,cofx,cofy,cofz)
c
c     gauss-seidel point relaxation with red/black ordering
c     in three dimensions for nonseparable pde
c     relax in order:
c     (1) red (x,y) on odd z planes
c     (2) black (x,y) on even z planes
c     (3) black (x,y) on odd z planes
c     (4) red (x,y) on even z planes
c
      implicit none
      integer nx,ny,nz
      real phi(0:nx+1,0:ny+1,0:nz+1),rhs(nx,ny,nz)
      real cofx(nx,3),cofy(ny,3),cofz(nz,3)
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      common/imud3sp/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nfx,nfy,nfz,iguess,maxcy,method,nwork,lwork,itero,
     +kcycle,iprer,ipost,intpol
      integer i,j,k,nper,ist,ifn,jst,jfn,kst,kfn
c
c     set periodic b.c. indicator
c
      nper = nxa*nyc*nze
c
c     set loop limits to avoid specified boundaries
c     in red/black sweeps
c
      ist = 1
      if (nxa.eq.1) ist = 3
      ifn = nx
      if (nxb.eq.1) ifn = nx-1
      jst = 1
      if (nyc.eq.1) jst = 3
      jfn = ny
      if (nyd.eq.1) jfn = ny-1
      kst = 1
      if (nze.eq.1) kst = 3
      kfn = nz
      if (nzf.eq.1) kfn = nz-1
c
c     set periodic boundaries if necessary
c
      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c   red (x,y) on odd z planes
c
!$OMP PARALLEL DO SHARED(ist,ifn,jst,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz)
!$OMP+PRIVATE(i,j,k)
      do k=kst,kfn,2
	do i=ist,ifn,2
	  do j=jst,jfn,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
	  end do
	end do
	do i=2,ifn,2
	  do j=2,jfn,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
	  end do
	end do
      end do
!$OMP END PARALLEL DO
      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c   black (x,y) or even z planes
c
!$OMP PARALLEL DO SHARED(ist,ifn,jst,jfn,kfn,phi,rhs,cofx,cofy,cofz)
!$OMP+PRIVATE(i,j,k)
      do k=2,kfn,2
	do i=ist,ifn,2
	  do j=2,jfn,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
	  end do
	end do
	do i=2,ifn,2
	  do j=jst,jfn,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
	  end do
	end do
      end do
!$OMP END PARALLEL DO
      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c   black (x,y) on odd z planes
c
!$OMP PARALLEL DO SHARED(ist,ifn,jfn,kst,kfn,phi,rhs,cofx,cofy,cofz)
!$OMP+PRIVATE(i,j,k)
      do k=kst,kfn,2
	do i=ist,ifn,2
	  do j=2,jfn,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
	  end do
	end do
	do i=2,ifn,2
	  do j=jst,jfn,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
	  end do
	end do
      end do
!$OMP END PARALLEL DO
      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c   red(x,y) on even z planes
c
!$OMP PARALLEL DO SHARED(ist,ifn,jst,jfn,kfn,phi,rhs,cofx,cofy,cofz)
!$OMP+PRIVATE(i,j,k)
      do k=2,kfn,2
	do i=ist,ifn,2
	  do j=jst,jfn,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
	  end do
	end do
	do i=2,ifn,2
	  do j=2,jfn,2
	    phi(i,j,k) = (rhs(i,j,k) - (
     +      cofx(i,1)*phi(i-1,j,k)+cofx(i,2)*phi(i+1,j,k) +
     +      cofy(j,1)*phi(i,j-1,k)+cofy(j,2)*phi(i,j+1,k) +
     +      cofz(k,1)*phi(i,j,k-1)+cofz(k,2)*phi(i,j,k+1)))/
     +      (cofx(i,3)+cofy(j,3)+cofz(k,3))
	  end do
	end do
      end do
!$OMP END PARALLEL DO
      if (nper.eq.0) call per3vb(nx,ny,nz,phi,nxa,nyc,nze)
      return
      end

