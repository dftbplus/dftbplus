c
c     file muh3.d
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
c
c     For MUDPACK 5.0 information, visit the website:
c     (http://www.scd.ucar.edu/css/software/mudpack)
c
c ... file muh3.d
c
c     contains documentation for:
c     subroutine muh3(iparm,fparm,wk,iw,coef,bndyc,rhs,phi,mgopt,ierror)
c     a sample fortran driver is file "tmuh3.f".
c
c ... author and specialist
c
c          john c. adams (national center for atmospheric research)
c          email: johnad@ucar.edu, phone: 303-497-1213
c
c ... required mudpack files
c
c     mudcom.f, mud3ln.f, mud3pn.f
c
c ... purpose
c
c     the "hybrid" multigrid/direct method code muh3 approximates the
c     same 3-d nonseparable elliptic pde as the mudpack solver mud3.
c     the basic algorithm is modified by using block banded gaussian
c     elimination in place of relaxation whenever the coarsest subgrid is
c     encountered within multigrid cycling.  use of the direct method at
c     the coarsest grid level gives muh3 at least two advantages over mud3:
c
c     (1) improved convergence rates
c
c         the use of a direct method at the coarsest grid level can
c         improve convergence rates at a small additional computational
c         cost if the coarse grid parameters ixp=iparm(8),jyq=iparm(9),
c         kzr=iparm(10) are small relative to the fine grid parameters
c         nx=iparm(14),ny=iparm(15),nz=iparm(16).  this is especially true
c         in the presence of certain boundary conditions.  for example,
c         if all boundary conditions are neuman (pure derivative) and/or
c         periodic then mud3 may fail to converge.  muh3 should handle
c         these boundary conditions with the expected multigrid efficiency
c         (see tmuh3.f).  in all cases, muh3 should give convergence
c         rates which equal or exceed mud3.
c
c
c     (2) more resolution choices
c
c         muh3 allows more grid size flexibility by "relaxing" the
c         constraint on the coarse grid parameters that ixp,jyq,kzr
c         be "very" small (2 or 3 for mud3) for effective error
c         error reduction within multigrid cycling.  convergence
c         rates will not deteoriate with larger values for ixp,
c         jyq,kzr.
c
c *** caution
c
c     because of the very large computational and storage
c     requirements, the three-dimensional dimensional direct
c     method costs can overwhelm the multigrid cycling costs
c     if the coarsest grid is not small relative to the finest
c     solution grid.  this is a user decision set by the choice
c     of coarse and fine grid parameters (see iparm(8) through
c     iparm(16) and iparm(21) descriptions)
c
c     subroutine muh3 automatically discretizes and attempts to compute
c     the second order finite difference approximation to a three-
c     dimensional linear nonseparable elliptic partial differential
c     equation on a box.  the approximation is generated on a uniform
c     grid covering the box (see mesh description below).  boundary
c     conditions may be any combination of mixed, specified (Dirchlet)
c     or periodic.  the form of the pde solved is . . .
c
c          cxx(x,y,z)*pxx + cyy(x,y,z)*pyy + czz(z,y,z)*pzz +
c
c          cx(x,y,z)*px + cy(x,y,z)*py + cz(x,y,z)*pz +
c
c          ce(x,y,z)*p(x,y,z) = r(x,y,z)
c
c     here cxx,cyy,czz,cx,cy,cz,ce are the known real coefficients
c     of the pde; pxx,pyy,pzz,px,py,pz are the second and first
c     partial derivatives of the unknown solution function p(x,y,z)
c     with respect to the independent variables x,y,z; r(x,y,z) is
c     is the known real right hand side of the elliptic pde.  cxx,cyy
c     and czz should be positive for all (x,y,z) in the solution region.
c                                                                               
c
c ... argument differences with mud3.f
c
c     the input and output arguments of muh3 are almost identical to the
c     arguments of mud3 (see mud3.d) with the following exceptions:
c
c     (1) let mx=ixp+1, my=jyq+1, mz=kzr+1 (the coarsest grid
c         resolutions, ixp=iparm(8), jyq=iparm(9), kzr=iparm(10))
c         then the work space vector "wk" requires
c
c               mx*my*mz*(2*mx*my+1))           (nze.ne.0)
c
c         additional words of storage if periodic boundary conditions
c         are not flagged in the z direction or
c
c               mx*my*(mz*(4*mx*my+1))          (nze=0)
c
c         additional words of storage if periodic boundary conditions are
c         flagged in the z direction (nze = 0).  the extra work space is
c         used for a direct solution with gaussian elimination whenever the
c         coarsest grid is encountered within multigrid cycling.
c
c     (2) an integer work space iwk of length at least mx*my*mz words
c         must be provided.  the length of iwk is not checked!
c
c     (3) kzr > 2 if nze=0, jyq > 2 if nyc=0, ixp > 2 if nxe = 0 are
c         required (i.e., the coarsest grid must contain at least four
c         points in any direction with periodic boundary conditions,
c         see the expanded meaning of ierror=3).
c
c *** (4) it is no longer necessary that ixp,jyq,kzr be 2 or 3 for
c         effective error reduction with multigrid iteration.  there
c         is no reduction in convergence rates when larger values for
c         ixp,jyq,kzr are used .  this provides additional flexibility
c         in choosing grid size.  in many cases muh3 provides more
c         robust convergence than mud3.  it can be used in place of
c         mud3 for all nonsingular problems (see (5) below).
c
c     (5) iguess = iparm(17) = 1 (flagging an initial guess) or
c         maxcy = iparm(18) > 1 (setting more than one multigrid
c         cycle) are not allowed if muh3 becomes a full direct method
c         by choosing iex=jey=kez=1.  this conflicting combination
c         of input arguments for multigrid iteration and a full
c         direct method set the fatal error flag
c
c              ierror = 13
c
c         iguess = 0 and maxcy = 1 are required when muh3 becomes a
c         full direct method.  ordinarily (see *** caution above) this
c         should not happen except when testing with very coarse resolution.
c
c
c     (6) if a "singular" pde is detected (see ierror=-3 description in
c         mud3.d, ce(x,y) = 0.0 for all x,y and the boundary conditions
c         are a combination of periodic and/or pure derivatives) then muh3
c         sets the fatal error flag
c
c              ierror = 14
c
c         the direct method utilized by muh3 would likely cause a near
c         division by zero in the singular case.  mud3 can be tried for
c         singular problems.
c
c ... grid size considerations
c
c     (1) flexibility
c
c         muh3 should be used in place of mud3 whenever grid size
c         requirements do not allow choosing ixp,jyq,kzr to be 2 or 3.
c
c     (2) additional work space (see (1) under "arguments differences") is
c         required by muh3 to implement gaussian elimination at the coarsest
c         grid level.  this may limit the size of ixp,jyq,kzr.
c
c     (3) operation counts
c                                                              k
c         for simplicity, assume p=ixp=jyq=kzr and n=nx=ny=nz=2 *p.
c         gaussian elimination requires o(p**7) operations for solution
c         on the coarsest subgrid while multigrid iteration is a o(n**3)
c         algorithm.  consequently the storage and computational
c         requirements for the 3-d direct method will dominate the
c         calculation if p is "large."  note that  o(p**7)=:o(n**3)
c         whenever k =: (4/3)*log2(p) grid levels are used in cycling.
c         larger values mean the direct method will dominate the
c         calculation.  smaller values for k mean the direct method
c         will only marginally add to the cost of multigrid iteration
c         alone.
c
c *** the remaining documentation is almost identical to mud3.d
c     except for the modifications already indicated.
c
c
c ... mesh description . . .
c
c     the approximation is generated on a uniform nx by ny by nz grid.
c     the grid is superimposed on the rectangular solution region
c
c          [xa,xb] x [yc,yd] x [ze,zf].
c
c     let
c
c          dlx = (xb-xa)/(nx-1), dly = (yd-yc)/(ny-1), dlz = (zf-ze)/(nz-1)
c
c     be the uniform grid increments in the x,y,z  directions. then
c
c          xi=xa+(i-1)*dlx,  yj=yc+(j-1)*dly, zk = ze+(k-1)*dlz
c
c     for i=1,...,nx; j=1,...,ny; k=1,...,nz  denote the x,y,z uniform
c     mesh points.
c
c
c ... language                                                                  
c
c     fortran90/fortran77
c
c
c ... portability                                                               
c
c     mudpack5.0 software has been compiled and tested with fortran77
c     and fortran90 on a variety of platforms.
c                                                                               
c ... methods
c
c     details of the methods employeed by the solvers in mudpack are given
c     in [1,9].  [1,2,9] contain performance measurements on a variety of
c     elliptic pdes (see "references" in the file "readme").  in summary:
c
c *** discretization and solution (second-order solvers) (see [1])
c
c     the pde and boundary conditions are automatically discretized at all
c     grid levels using second-order finite difference formula.  diagonal
c     dominance at coarser grid levels is maintained in the presence of
c     nonzero first-order terms by adjusting the second-order coefficient
c     when necessary.  the resulting block tri-diagonal linear system is
c     approximated using multigrid iteration [10,11,13,15,16,18].  version
c     5.0 of mudpack uses only fully weighted residual restriction.  defaults
c     include cubic prolongation and w(2,1) cycles.  these can be overridden
c     with  selected multigrid options (see "mgopt").  error control based on
c     maximum relative differences is available. full multigrid cycling (fmg)
c     or cycling beginning or restarting at the finest grid level can be
c     selected. a menu of relaxation methods including gauss-seidel point,
c     line relaxation(s) (in any combination of directions) and planar
c     relaxation (for three-dimensional anisotropic problems) are provided.
c     all methods use ordering based on alternating points (red/black),
c     lines, or planes for cray vectorization and improved convergence
c     rates [14].
c
c *** higher order solution (fourth-order solvers) (see [9,19,21])
c
c     if the multigrid cycling results in a second-order estimate (i.e.,
c     discretization level error is reached) then this can be improved to a
c     fourth-order estimate using the technique of "deferred corrections."
c     the values in the solution array are used to generate a fourth-order
c     approximation to the truncation error.  second-order finite difference
c     formula are used to approximate third and fourth partial derivatives
c     of the solution function [3].  the truncation error estimate is
c     transferred down to all grid levels using weighted averaging where
c     it serves as a new right hand side.  the default multigrid options
c     are used to compute the fourth-order correction term which is added
c     to the original solution array.
c
c
c ... references (partial)
c
c     for a complete list see "references" in the mudpack information and
c     directory file "readme" or visit the mudpack web site
c     (http://www.scd.ucar.edu/css/software/mudpack)
c
c     [1] J. Adams, "MUDPACK: Multigrid Fortran Software for the Efficient
c     Solution of Linear Elliptic Partial Differential Equations,"
c     Applied Math. and Comput. vol.34, Nov 1989, pp.113-146.
c
c     [2] J. Adams,"FMG Results with the Multigrid Software Package MUDPACK,"
c     proceedings of the fourth Copper Mountain Conference on Multigrid, SIAM,
c     1989, pp.1-12.
c     .
c     .
c     .
c     [7] J. Adams, R. Garcia, B. Gross, J. Hack, D. Haidvogel, and V. Pizzo,
c     "Applications of Multigrid Software in the Atmospheric Sciences,"
c      Mon. Wea. Rev.,vol. 120 # 7, July 1992, pp. 1447-1458.
c     .
c     .
c     .
c     [9] J. Adams, "Recent Enhancements in MUDPACK, a Multigrid Software
c     package for Elliptic Partial Differential Equations," Applied Math.
c     and Comp., 1991, vol. 43, May 1991, pp. 79-94.
c
c     [10]J. Adams, "MUDPACK-2: Multigrid Software for Approximating
c     Elliptic Partial Differential Equations on Uniform Grids with
c     any Resolution," Applied Math. and Comp., 1993, vol. 53, February
c     1993, pp. 235-249
c     .
c     .
c     .
c
c ... argument description
c                                                                               
c
c **********************************************************************
c *** input arguments *************************************************
c **********************************************************************
c
c
c ... iparm                                                                     
c
c          an integer vector of length 23 used to efficiently pass
c          integer arguments.  iparm is set internally in muh3
c          and defined as follows . . .                                         
c
c
c ... intl=iparm(1)
c
c          an initialization argument.  intl=0  must be input
c          on an initial call. in this case input arguments will
c          be checked for errors and the elliptic partial differential
c          equation and boundary conditions will be discretized using
c          second order finite difference formula.
c
c ***      an approximation is not generated after an intl=0 call!
c          muh3 should be called with intl=1 to approximate the elliptic
c          pde discretized by the intl=0 call.  intl=1 should also
c          be input if muh3 has been called earlier and only the
c          values in in rhs (see below) or gbdy (see bndyc below)
c          or phi (see below) have changed.  this will bypass
c          redundant pde discretization and argument checking
c          and save computational time.  some examples of when
c          intl=1 calls should be used are:
c
c          (0) after a intl=0 argument checking and discretization call
c
c          (1) muh3 is being recalled for additional accuracy.  in
c              this case iguess=iparm(12)=1 should also be used.
c
c          (2) muh3 is being called every time step in a time dependent
c              problem (see discussion below) where the elliptic operator
c              does not depend on time.
c
c          (3) muh3 is being used to solve the same elliptic equation
c              for several different right hand sides (iguess=0 should
c              probably be used for each new righthand side).
c
c          intl = 0 must be input before calling with intl = 1 when any
c          of the following conditions hold:
c
c          (a) the initial call to muh3

c          (b) any of the integer arguments other than iguess=iparm(12)
c              or maxcy=iparm(13) or mgopt have changed since the previous
c              call.
c
c          (c) any of the floating point arguments other than tolmax=
c              fparm(5) have changed since the previous call
c
c          (d) any of the coefficients input by coef (see below) have
c              changed since the previous call
c
c          (e) any of the "alfa" coefficients input by bndyc (see below)
c              have changed since the previous call.
c
c          if any of (a) through (e) are true then the elliptic pde
c          must be discretized or rediscretized.  if none of (a)
c          through (e) holds, calls can be made with intl=1.
c          incorrect calls with intl=1 will produce erroneous results.
c  ***     the values set in the saved work space "wk" (see below) with
c          an intl=0 call must be preserved with subsequent intl=1 calls.
c
c          MUDPACK software performance should be monitored for intl=1
c          calls.  The intl=0 discretization call performance depends
c          primarily on the efficiency or lack of efficiency of the
c          user provided subroutines for pde coefficients and
c          boundary conditions.
c
c
c ... nxa=iparm(2)
c
c          flags boundary conditions on the (y,z) plane x=xa                    
c
c        = 0 if p(x,y,z) is periodic in x on [xa,xb]
c            (i.e., p(x+xb-xa,y,z) = p(x,y,z) for all x,y,z)
c
c        = 1 if p(xa,y,z) is specified  (this must be input thru phi(1,j,k))
c
c        = 2 if there are mixed derivative boundary conditions at x=xa
c            (see "bndyc" description below where kbdy = 1)
c
c
c ... nxb=iparm(3)
c
c          flags boundary conditions on the (y,z) plane x=xb                    
c
c        = 0 if p(x,y,z) is periodic in x on [xa,xb]
c            (i.e., p(x+xb-xa,y,z) = p(x,y,z) for all x,y,z)
c
c        = 1 if p(xb,y,z) is specified (this must be input thru phi(nx,j,k))
c
c        = 2 if there are mixed derivative boundary conditions at x=xb
c            (see "bndyc" description below where kbdy = 2)
c
c
c ... nyc=iparm(4)
c
c          flags boundary conditions on the (x,z) plane y=yc                    
c
c        = 0 if p(x,y,z) is periodic in y on [yc,yd]
c            (i.e., p(x,y+yd-yc,z) = p(x,y,z) for all x,y,z)

c        = 1 if p(x,yc,z) is specified (this must be input thru phi(i,1,k))

c        = 2 if there are mixed derivative boundary conditions at y=yc
c            (see "bndyc" description below where kbdy = 3)
c
c
c ... nyd=iparm(5)
c
c          flags boundary conditions on the (x,z) plane y=yd                    
c
c        = 0 if p(x,y,z) is periodic in y on [yc,yd]
c            (i.e., p(x,y+yd-yc,z) = p(x,y,z) for all x,y,z)

c        = 1 if p(x,yd,z) is specified (this must be input thru phi(i,ny,k))

c        = 2 if there are mixed derivative boundary conditions at y=yd
c            (see "bndyc" description below where kbdy = 4)
c
c
c ... nze=iparm(6)
c
c          flags boundary conditions on the (x,y) plane z=ze                    
c
c        = 0 if p(x,y,z) is periodic in z on [ze,zf]
c            (i.e., p(x,y,z+zf-ze) = p(x,y,z) for all x,y,z

c        = 1 if p(x,y,ze) is specified (this must be input thru phi(i,j,1))

c        = 2 if there are mixed derivative boundary conditions at z=ze
c            (see "bndyc" description below where kbdy = 5)
c
c
c ... nzf=iparm(7)
c
c          flags boundary conditions on the (x,y) plane z=zf                    
c
c        = 0 if p(x,y,z) is periodic in z on [ze,zf]
c            (i.e., p(x,y,z+zf-ze) = p(x,y,z) for all x,y,z

c        = 1 if p(x,y,zf) is specified (this must be input thru phi(i,j,nz))

c        = 2 if there are mixed derivative boundary conditions at z=zf
c            (see "bndyc" description below where kbdy = 6)
c
c
c *** grid size arguments
c
c
c ... ixp = iparm(8)
c
c          an integer greater than one which is used in defining the number
c          of grid points in the x direction (see nx = iparm(14)).  "ixp+1"
c          is the number of points on the coarsest x grid visited during
c          multigrid cycling.  ixp should be chosen as small as possible
c          within grid size requirements.
c
c
c ... jyq = iparm(9)
c
c          an integer greater than one which is used in defining the number
c          of grid points in the y direction (see ny = iparm(15)).  "jyq+1"
c          is the number of points on the coarsest y grid visited during
c          multigrid cycling.  jyq should be chosen as small as possible
c          within grid size requirements.
c
c
c ... kzr = iparm(10)
c
c          an integer greater than one which is used in defining the number
c          of grid points in the z direction (see nz = iparm(16)).  "kzr+1"
c          is the number of points on the coarsest z grid visited during
c          multigrid cycling.  kzr should be chosen as small as possible
c          within grid size requirements.
c
c
c ... iex = iparm(11)
c
c          a positive integer exponent of 2 used in defining the number
c          of grid points in the x direction (see nx = iparm(14)).
c          iex .le. 50 is required.  for efficient multigrid cycling,
c          iex should be chosen as large as possible and ixp=iparm(8)
c          as small as possible within grid size constraints when
c          defining nx = iparm(14).
c
c
c ... jey = iparm(12)
c
c          a positive integer exponent of 2 used in defining the number
c          of grid points in the y direction (see ny = iparm(15)).
c          jey .le. 50 is required.  for efficient multigrid cycling,
c          jey should be chosen as large as possible and jyq=iparm(9)
c          as small as possible within grid size constraints when
c          defining ny = iparm(15).
c
c
c ... kez = iparm(13)
c
c          a positive integer exponent of 2 used in defining the number
c          of grid points in the z direction (see nz = iparm(16)).
c          kez .le. 50 is required.  for efficient multigrid cycling,
c          kez should be chosen as large as possible and kzr=iparm(10)
c          as small as possible within grid size constraints when
c          defining nz = iparm(16).
c
c
c ... nx = iparm(14)
c
c          the number of equally spaced grid points in the interval [xa,xb]
c          (including the boundaries).  nx must have the form
c
c               nx = ixp*(2**(iex-1)) + 1
c
c          where ixp = iparm(8), iex = iparm(11).
c
c
c ... ny = iparm(15)
c
c          the number of equally spaced grid points in the interval [yc,yd]
c          (including the boundaries).  ny must have the form:
c
c               ny = jyq*(2**(jey-1)) + 1
c
c          where jyq = iparm(9), jey = iparm(12).
c
c
c ... nz = iparm(16)
c
c          the number of equally spaced grid points in the interval [ze,zf]
c          (including the boundaries).  nz must have the form
c
c               nz = kzr*(2**(kez-1)) + 1
c
c          where kzr = iparm(10), kez = iparm(13)
c
c
c *** note
c
c     let g be the nx by ny by nz fine grid on which the approximation is
c     generated and let n = max0(iex,jey,kez).  in mudpack, multigrid
c     cycling is implemented on the ascending chain of grids
c
c         g(1) < ... < g(k) < ... < g(n) = g.
c
c     each g(k) (k=1,...,n) has mx(k) by my(k) by mz(k) grid points
c     given by:
c
c         mx(k) = ixp*[2**(max0(iex+k-n,1)-1)] + 1
c
c         my(k) = jyq*[2**(max0(jey+k-n,1)-1)] + 1
c
c         mz(k) = kzr*[2**(max0(kez+k-n,1)-1)] + 1
c
c     additionally muh3 implements a direct method whenever g(1) is
c     encountered.
c
c ... iguess=iparm(17)
c
c          = 0 if no initial guess to the pde is provided
c              and/or full multigrid cycling beginning at the
c              coarsest grid level is desired.
c
c          = 1 if an initial guess to the pde at the finest grid
c              level is provided in phi (see below).  in this case
c              cycling beginning or restarting at the finest grid
c              is initiated.
c
c *** comments on iguess = 0 or 1 . . .
c
c
c     setting iguess=0 forces full multigrid or "fmg" cycling.  phi
c     must be initialized at all grid points.  it can be set to zero at
c     non-Dirchlet grid points if nothing better is available.
c
c     if iguess = 1 then the values input in phi are an initial guess to the
c     pde at the finest grid level where cycling begins.  this option should
c     be used only if a "very good" initial guess is available (as, for
c     example, when restarting from a previous iguess=0 call).
c
c *** time dependent problems . . .
c
c     assume we are solving an elliptic pde every time step in a
c     marching problem of the form:
c
c          l(p(t)) = r(t)
c
c     where the differential operator "l" has no time dependence,
c     "p(t)" is the solution and "r(t)" is the right hand side at
c     current time "t". let "dt" be the increment between time steps.
c     then p(t) can be used as an initial guess to p(t+dt) with
c     intl = 1 when solving
c
c          l(p(t+dt)) = r(t+dt).
c
c     after the first two time steps, rather than continue, it would
c     be better to define the "correction" term:
c
c          e(t,dt) = p(t+dt) - p(t)
c
c     this clearly satisfies the equation
c
c          l(e(t,dt)) = r(t+dt) - r(t).
c
c     this should be solved with iguess = 0 and intl = 1. boundary
c     conditions for e(t,dt) are obtained from the boundary conditions
c     for p(t) by subtracting given values at t from given values at
c     t+dt. for example if
c
c          d(p(t))/dx = f(t), d(p(t+dt))/dx = f(t+dt)
c
c     at some x boundary then e(t,dt) satisfies the derivative
c     boundary condition
c
c          d(e(t,dt))/dx = f(t+dt) - f(t).
c
c     e(t,dt) can be preset to 0.0 (at non-Dirchlet points) or (if p(t-dt)
c     is saved) to p(t)-p(t-dt).  with iguess = 0, these values will serve
c     as an initial guess to e(t,dt) at the coarsest grid level.  this
c     approach has the advantage that a full sequence of multigrid cycles,
c     beginning at the coarsest grid level, is invoked every time step in
c     solving for e(t,dt).  a few digits of accuracy in e(t,dt), which is
c     ordinarily much smaller than p(t), will yield several more digits of
c     accuracy in the final approximation:
c
c          p(t+dt) = p(t) + e(t,dt).
c
c     using this approach to integrate in time will give more accuracy
c     then using p(t) as an initial guess to p(t+dt) for all time steps.
c     it does require additional storage.
c
c     if the differential operator "l" has time dependence (either thru
c     the coefficients in the pde or the coefficients in the derivative
c     boundary conditions) then use p(t) as an initial guess to p(t+dt)
c     when solving
c
c          l(t+dt)(p(t+dt)) = r(t+dt)
c
c     with intl = 0 for all time steps (the discretization must be repeated
c     for each new "t"). either iguess = 0 (p(t) will then be an initial
c     guess at the coarsest grid level where cycles will commence) or
c     iguess = 1 (p(t) will then be an initial guess at the finest grid
c     level where cycles will remain fixed) can be tried.
c
c
c ... maxcy = iparm(18)
c
c          the exact number of cycles executed between the finest
c          (nx by ny by nz) and the coarsest ((ixp+1) by (jyq+1) by
c          (kzr+1)) grid levels when tolmax=fparm(7)=0.0 (no error
c          control). when tolmax=fparm(7).gt.0.0 is input (error control)
c          then maxcy is a limit on the number of cycles between the
c          finest and coarsest grid levels.  in any case, at most
c          maxcy*(iprer+ipost) relaxation sweeps are performed at the
c          finest grid level (see iprer=mgopt(2),ipost=mgopt(3) below)
c          when multigrid iteration is working "correctly" only a few
c          cycles are required for convergence.  large values for maxcy
c          should not be required.
c
c
c ... method = iparm(19)
c
c          this sets the method of relaxation (all relaxation
c          schemes in mudpack use red/black type ordering)
c
c          = 0 for gauss-seidel pointwise relaxation
c
c          = 1 for line relaxation in the x direction
c
c          = 2 for line relaxation in the y direction
c
c          = 3 for line relaxation in the z direction
c
c          = 4 for line relaxation in the x and y direction
c
c          = 5 for line relaxation in the x and z direction
c
c          = 6 for line relaxation in the y and z direction
c
c          = 7 for line relaxation in the x,y and z direction
c
c          = 8 for x,y planar relaxation
c
c          = 9 for x,z planar relaxation
c
c          =10 for y,z planar relaxation
c
c ***  choice of method
c
c      this is very important for efficient convergence.  in some cases
c      experimentation may be required.
c
c      let fx represent the quantity cxx(x,y,z)/dlx**2 over the solution box
c
c      let fy represent the quantity cyy(x,y,z)/dly**2 over the solution box
c
c      let fz represent the quantity czz(x,y,z)/dlz**2 over the solution box
c
c      (0) if fx,fy,fz are roughly the same size and do not vary too
c          much choose method = 0.  if this fails try method = 7.
c
c      (1) if fx is much greater then fy,fz and fy,fz are roughly the same
c          size choose method = 1
c
c      (2) if fy is much greater then fx,fz and fx,fz are roughly the same
c          size choose method = 2
c
c      (3) if fz is much greater then fx,fy and fx,fy are roughly the same
c          size choose method = 3
c
c      (4) if fx,fy are roughly the same and both are much greater then fz
c          try method = 4.  if this fails try method = 8
c
c      (5) if fx,fz are roughly the same and both are much greater then fy
c          try method = 5.  if this fails try method = 9
c
c      (6) if fy,fz are roughly the same and both are much greater then fx
c          try method = 6.  if this fails try method = 10
c
c      (7) if fx,fy,fz vary considerably with none dominating try method = 7
c
c      (8) if fx and fy are considerably greater then fz but not necessarily
c          the same size and method=4 fails try method = 8
c
c      (9) if fx and fz are considerably greater then fy but not necessarily
c          the same size and method=5 fails try method = 9
c
c      (10)if fy and fz are considerably greater then fx but not necessarily
c          the same size and method=6 fails try method = 10
c
c
c ... meth2 = iparm(20) determines the method of relaxation used in the planes
c             when method = 8 or 9 or 10.
c
c
c          as above, let fx,fy,fz represent the quantities cxx/dlx**2,
c          cyy/dly**2,czz/dlz**2 over the box.
c
c          (if method = 8)
c
c          = 0 for gauss-seidel pointwise relaxation
c              in the x,y plane for each fixed z
c          = 1 for line relaxation in the x direction
c              in the x,y plane for each fixed z
c          = 2 for line relaxation in the y direction
c              in the x,y plane for each fixed z
c          = 3 for line relaxation in the x and y direction
c              in the x,y plane for each fixed z
c
c          (1) if fx,fy are roughly the same and vary little choose meth2 = 0
c          (2) if fx is much greater then fy choose meth2 = 1
c          (3) if fy is much greater then fx choose meth2 = 2
c          (4) if none of the above or meth2 = 0 fails choose meth2 = 3
c
c          (if method = 9)
c
c          = 0 for gauss-seidel pointwise relaxation with red/black ordering
c              in the x,z plane for each fixed y
c          = 1 for simultaneous line relaxation in the x direction
c              of the x,z plane for each fixed y
c          = 2 for simultaneous line relaxation in the z direction
c              of the x,z plane for each fixed y
c          = 3 for simultaneous line relaxation in the x and z direction
c              of the x,z plane for each fixed y
c
c          (1) if fx,fz are roughly the same and vary little choose meth2 = 0
c          (2) if fx is much greater then fz choose meth2 = 1
c          (3) if fz is much greater then fx choose meth2 = 2
c          (4) if none of the above or meth2 = 0 fails choose meth2 = 3
c
c          (if method = 10)
c
c          = 0 for gauss-seidel pointwise relaxation with red/black ordering
c              in the y,z plane for each fixed x
c          = 1 for simultaneous line relaxation in the y direction
c              of the y,z plane for each fixed x
c          = 2 for simultaneous line relaxation in the z direction
c              of the y,z plane for each fixed x
c          = 3 for simultaneous line relaxation in the y and z direction
c              of the y,z plane for each fixed x
c
c          (1) if fy,fz are roughly the same and vary little choose meth2 = 0
c          (2) if fy is much greater then fz choose meth2 = 1
c          (3) if fz is much greater then fy choose meth2 = 2
c          (4) if none of the above or meth2 = 0 fails choose meth2 = 3
c
c
c ***  note that planar relaxation implements full two-dimensional multigrid
c      cycling for each plane visited during three dimensional multigrid
c      cycling.  Consequently it can be computationally expensive.
c
c ... length = iparm(21)
c
c          the length of the work space provided in vector wk.
c
c          let isx = 3 if method = 1,4,5 or 7 and nxa.ne.0
c          let isx = 5 if method = 1,4,5 or 7 and nxa.eq.0
c          let isx = 0 if method has any other value
c
c          let jsy = 3 if method = 2,4,6 or 7 and nyc.ne.0
c          let jsy = 5 if method = 2,4,6 or 7 and nyc.eq.0
c          let jsy = 0 if method has any other value
c
c          let ksz = 3 if method = 3,5,6 or 7 and nze.ne.0
c          let ksz = 5 if method = 3,5,6 or 7 and nze.eq.0
c          let ksz = 0 if method has any other value
c
c
c
c          let ls = (nx+2)*(ny+2)*(nz+2)*(10+isx+jsy+ksz)
c
c          let mx = ixp+1; my = jyq+1; mz = kzr+1.  the block gaussian
c          elimination at the coarsest mx by my by mz grid level requires
c
c             ld = mx*my*mz*(2*mx*my+1))           (nze.ne.0)
c
c          words of storage if z boundary conditions are not periodic or
c
c             ld = mx*my*(mz*(4*mx*my+1))          (nze=0)
c
c          words of storage if z boundary conditions are periodic.
c          if ixp,jyq,kzr are not the same, this quantity is
c          minimized if they are chosen so that kzr > max0(ixp,jyq).
c
c          finally
c
c             length = ls + ld

c          will usually but not always suffice.  the exact minimal length depends,
c          in a complex way, on the grid size arguments and method chosen.
c  ***     it can be predetermined for the current input arguments by calling
c          muh3 with length=iparm(21)=0 and printing iparm(22) or (in f90)
c          dynamically allocating the work space using the value in iparm(22)
c          in a subsequent muh3 call.
c
c ... fparm                                                                     
c
c          a floating point vector of length 8 used to efficiently
c          pass floating point arguments.  fparm is set internally
c          in muh3 and defined as follows . . .
c
c
c ... xa=fparm(1), xb=fparm(2)
c
c          the range of the x independent variable. xa must                     
c          be less than xb                                                      
c
c
c ... yc=fparm(3), yd=fparm(4)
c
c          the range of the y independent variable.  yc must                    
c          be less than yd.                                                     
c
c
c ... ze=fparm(5), zf=fparm(6)
c
c          the range of the z independent variable. ze must                     
c          be less than zf.                                                     
c
c
c ... tolmax = fparm(7)
c
c          when input positive, tolmax is a maximum relative error tolerance
c          used to terminate the relaxation iterations. assume phi1(i,j,k)
c          and phi2(i,j,k) are the last two computed approximations at all
c          grid points of the finest grid level. if we define
c
c              phdif = max(abs(phi2(i,j,k)-phi1(i,j,k))) for all i,j,k
c
c          and
c
c              phmax = max(abs(phi2(i,j,k))) for all i,j,k
c
c          then "convergence" is considered to have occurred if and only if
c
c              phdif/phmax < tolmax.
c
c
c          if tolmax=fparm(7)=0.0 is input then there is no error control
c          and maxcy cycles from the finest grid level are executed. maxcy
c          is a limit which cannot be exceeded even with error control.
c     ***  calls with tolmax=0.0, when appropriate because of known
c          convergence behavior, are more efficient than calls with tolmax
c          positive (i.e., if possible DO NOT use error control!).
c
c
c ... wk
c
c          a one dimensional array that must be provided for work space.        
c          see length = iparm(21). the values in wk must be preserved
c          if muh3 is called again with intl=iparm(1).ne.0 or if muh34
c          is called to improve accuracy.
c
c
c ... iwk
c
c          an integer vector dimensioned of length at least
c
c             (ixp+1)*(jyq+1)*(kzr+1)
c
c          in the routine calling muh3. the length of iwk is not
c          checked!  if iwk has length too small then undetectable
c          undetectable errors will result.
c
c
c ... bndyc                                                                     
c
c          a subroutine with arguments (kbdy,xory,yorz,alfa,gbdy).
c          which are used to input mixed boundary conditions to muh3.
c          the boundaries are numbered one thru six and the form of
c          conditions are described below.
c                                                                               
c
c     (1) the kbdy=1 boundary
c
c     this is the (y,z) plane x=xa where nxa=iparm(2) = 2 flags
c     a mixed boundary condition of the form
c
c        dp/dx + alfxa(y,z)*p(xa,y,z) = gbdxa(y,z)
c
c     in this case kbdy=1,xory=y,yorz=z will be input to bndyc and
c     alfa,gbdy corresponding to alfxa(y,z),gbdxa(y,z) must be returned.
c
c
c     (2) the kbdy=2 boundary
c
c     this is the (y,z) plane x=xb where nxb=iparm(3) = 2 flags
c     a mixed boundary condition of the form
c
c        dp/dx + alfxb(y,z)*p(xb,y,z) = gbdxb(y,z)
c
c     in this case kbdy=2,xory=y,yorz=z will be input to bndyc and
c     alfa,gbdy corresponding to alfxb(y,z),gbdxb(y,z) must be returned.
c
c
c     (3) the kbdy=3 boundary
c
c     this is the (x,z) plane y=yc where nyc=iparm(4) = 2 flags
c     a mixed boundary condition of the form
c
c        dp/dy + alfyc(x,z)*p(x,yc,z) = gbdyc(x,z)
c
c     in this case kbdy=3,xory=x,yorz=z will be input to bndyc and
c     alfa,gbdy corresponding to alfyc(x,z),gbdyc(x,z) must be returned.
c
c
c     (4) the kbdy=4 boundary
c
c     this is the (x,z) plane y=yd where nyd=iparm(5) = 2 flags
c     a mixed boundary condition of the form
c
c        dp/dy + alfyd(x,z)*p(x,yd,z) = gbdyd(x,z)
c
c     in this case kbdy=4,xory=x,yorz=z will be input to bndyc and
c     alfa,gbdy corresponding to alfyd(x,z),gbdyd(x,z) must be returned.
c
c                                                                               
c     (5) the kbdy=5 boundary
c
c     this is the (x,y) plane z=ze where nze=iparm(6) = 2 flags
c     a mixed boundary condition of the form
c
c        dp/dz + alfze(x,y)*p(x,y,ze) = gbdze(x,y)
c
c     in this case kbdy=5,xory=x,yorz=y will be input to bndyc and
c     alfa,gbdy corresponding to alfze(x,y),gbdze(x,y) must be returned.
c
c
c     (6) the kbdy=6 boundary
c
c     this is the (x,y) plane z=zf where nzf=iparm(7) = 2 flags
c     a mixed boundary condition of the form
c
c        dp/dz + alfzf(x,y)*p(x,y,zf) = gbdzf(x,y)
c
c     in this case kbdy=6,xory=x,yorz=y will be input to bndyc and
c     alfa,gbdy corresponding to alfzf(x,y),gbdzf(x,y) must be returned.
c
c                                                                               
c *** alfxa,alfyc,alfze nonpositive and alfxb,alfyd,alfze nonnegative
c     will help maintain matrix diagonal dominance during discretization
c     aiding convergence.
c
c *** bndyc must provide the mixed boundary condition
c     values in correspondence with those flagged in iparm(2)
c     thru iparm(7).  if all boundaries are specified then
c     muh3 will never call bndyc.  even then it must be entered
c     as a dummy subroutine. bndyc must be declared
c     external in the routine calling muh3.  the actual
c     name chosen may be different.
c
c
c ... coef
c
c         a subroutine with arguments (x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
c         which provides the known real coefficients for the elliptic pde
c         at any grid point (x,y,z).  the name chosen in the calling routine
c         may be different where the coefficient routine must be declared
c         external.
c
c ... rhs                                                                       
c
c          an array dimensioned nx by ny by nz which contains
c          the given right hand side values on the uniform 3-d mesh.            
c          rhs(i,j,k) = r(xi,yj,zk) for i=1,...,nx and j=1,...,ny               
c          and k=1,...,nz.
c
c ... phi                                                                       
c
c          an array dimensioned nx by ny by nz .  on input phi must
c          contain specified boundary values and an initial guess
c          to the solution if flagged (see iguess=iparm(17)=1).  for
c          example, if nyd=iparm(5)=1 then phi(i,ny,k) must be set
c          equal to p(xi,yd,zk) for i=1,...,nx and k=1,...,nz prior to
c          calling muh3.  the specified values are preserved by muh3.
c
c ***      if no initial guess is given (iguess=0) then phi must still
c          be initialized at non-Dirchlet grid points (this is not
c          checked). these values are projected down and serve as an initial
c          guess to the pde at the coarsest grid level.  set phi to 0.0 at
c          nonDirchlet grid points if nothing better is available.
c
c
c ... mgopt
c
c           an integer vector of length 4 which allows the user to select
c           among various multigrid options.  if mgopt(1)=0 is input then
c           a default set of multigrid arguments (chosen for robustness)
c           will be internally selected and the remaining values in mgopt
c           will be ignored.  if mgopt(1) is nonzero then the arguments
c           in mgopt are set internally and defined as follows:  (see the
c           basic coarse grid correction algorithm below)
c
c
c     kcycle = mgopt(1)
c
c            = 0 if default multigrid options are to be used
c
c            = 1 if v cycling is to be used (the least expensive per cycle)
c
c            = 2 if w cycling is to be used (the default)
c
c            > 2 if more general k cycling is to be used
c             *** warning--values larger than 2 increase
c                 the execution time per cycle considerably and
c                 result in the nonfatal error ierror = -5
c                 which indicates inefficient multigrid cycling.
c
c    iprer = mgopt(2)
c
c           the number of "pre-relaxation" sweeps executed before the
c           residual is restricted and cycling is invoked at the next
c           coarser grid level (default value is 2 whenever mgopt(1)=0)
c
c    ipost = mgopt(3)
c
c           the number of "post relaxation" sweeps executed after cycling
c           has been invoked at the next coarser grid level and the residual
c           correction has been transferred back (default value is 1
c           whenever mgopt(1)=0).
c
c *** if iprer, ipost, or (especially) kcycle is greater than 2
c     than inefficient multigrid cycling has probably been chosen and
c     the nonfatal error (see below) ierror = -5 will be set.  note
c     this warning may be overridden by any other nonzero value
c     for ierror.
c
c   intpol = mgopt(4)
c
c          = 1 if multilinear prolongation (interpolation) is used to
c              transfer residual corrections and the pde approximation
c              from coarse to fine grids within full multigrid cycling.
c
c          = 3 if multicubic prolongation (interpolation) is used to
c              transfer residual corrections and the pde approximation
c              from coarse to fine grids within full multigrid cycling.
c              (this is the default value whenever mgopt(1)=0).
c
c *** the default values (2,2,1,3) in the vector mgopt were chosen for
c     robustness.  in some cases v(2,1) cycles with linear prolongation will
c     give good results with less computation (especially in two-dimensions).
c     this  was the default and only choice in an earlier version of mudpack
c     (see [1]) and can be set with the integer vector (1,2,1,1) in mgopt.
c
c *** the schedules for one full multigrid cycle (iguess=0) using v(2,1)
c     cycles and w(2,1) cycles are depicted for a four level grid below.
c     the number of relaxation sweeps when each grid is visited are indicated.
c     the "*" stands for prolongation of the full approximation and the "."
c     stands for transfer of residuals and residual corrections within the
c     coarse grid correction algorithm (see below).  the "d" at level 1
c     indicates a direct method is used
c
c     one fmg with v(2,1) cycles:
c
c
c     ------------------------------2-----------------1------     level 4
c                                  * .               .
c                                 *   .             .
c     ---------------2-----------1-----2-----------1---------     level 3
c                   * .         .       .         .
c                  *   .       .         .       .
c     ------2-----1-----2-----1-----------2-----1------------     level 2
c          * .   .       .   .             .   .
c         *   . .         . .               . .
c     ---d-----d-----------d-----------------d---------------     level 1
c
c
c     one fmg with w(2,1) cycles:
c
c     ------------------------2---------------------------1--     level 4
c                            * .                         .
c     ----------2-----------1---2-----------3-----------1----     level 3
c              * .         .     .         . .         .
c     ----2---1---2---3---1-------2---3---1---2---3---1------     level 2
c        * . .     . . . .         . . . .     . . . .
c     --d---d-------d---d-----------d---d-------d---d--------     level 1
c
c
c     the form of the "recursive" coarse grid correction cycling used
c     when kcycle.ge.0 is input is described below in pseudo-algorithmic
c     language.  it is implemented non-recursively in fortran in mudpack.
c *** this algorithim is modified with the hybrid solvers which use
c     a direct method whenever grid level 1 is encountered.
c
c     algorithm cgc(k,l(k),u(k),r(k),kcycle,iprer,ipost,iresw,intpol)
c
c *** approximately solve l(k)*u(k) = r(k) using multigrid iteration
c *** k is the current grid level
c *** l(k) is the discretized pde operator at level k
c *** u(k) is the initial guess at level k
c *** r(k) is the right hand side at level k
c *** i(k,k-1) is the restriction operator from level k to level k-1
c *** (the form of i(k,k-1) depends on iresw)
c *** i(k-1,k) is the prolongation operator from level k-1 to level k
c *** (the form of i(k-1,k) depends on intpol)
c
c     begin algorithm cgc
c
c ***   pre-relax at level k
c
c     . do (i=1,iprer)
c
c     . .  relax(l(k),u(k),r(k))
c
c     . end do
c
c     . if (k > 1) then
c
c ***     restrict the residual from level k to level k-1
c
c     . . r(k-1) = i(k,k-1)(r(k)-l(k)*u(k))
c
c     . . kount = 0
c
c     . . repeat
c
c ***     solve for the residual correction at level k-1 in u(k-1)
c ***     using algorithm cgc "kcycle" times (this is the recursion)
c
c     . . . kount = kount+1
c
c     . . . invoke cgc(k-1,l(k-1),u(k-1),r(k-1),kcycle,iprer,ipost,iresw)
c
c
c     . . until (kount.eq.kcycle)
c
c ***     transfer residual correction in u(k-1) to level k
c ***     with the prolongation operator and add to u(k)
c
c     . . u(k) = u(k) + i(k-1,k)(u(k-1))
c
c     . end if
c
c ***   post relax at level k
c
c     . do (i=1,ipost)
c
c     . . relax(l(k),u(k),r(k))
c
c     . end do
c
c     . return
c
c     end algorithm cgc
c
c ***********************************************************************
c ****output arguments**************************************************
c ***********************************************************************
c
c
c ... iparm(22)
c
c          on output iparm(22) contains the actual work space length
c          required for the current grid sizes and method.  this value
c          will be computed and returned even if iparm(21) is less then
c          iparm(22) (see ierror=9).
c
c
c ... iparm(23)
c
c          if error control is selected (tolmax = fparm(7) .gt. 0.0) then
c          on output iparm(23) contains the actual number of cycles executed
c          between the coarsest and finest grid levels in obtaining the
c          approximation in phi.  the quantity (iprer+ipost)*iparm(23) is
c          the number of relaxation sweeps performed at the finest grid level.
c
c
c ... fparm(8)                                                                  
c
c          on output fparm(8) contains the final computed maximum relative
c          difference between the last two iterates at the finest grid level.
c          fparm(8) is computed only if there is error control (tolmax.gt.0.)
c          assume phi1(i,j,k) and phi2(i,j,k) are the last two computed
c          values for phi(i,j,k) at all points of the finest grid level.
c          if we define
c
c             phdif = max(abs(phi2(i,j,k)-phi1(i,j,k))) for all i,j,k
c
c          and
c
c             phmax = max(abs(phi2(i,j,k))) for all i,j,k
c
c          then
c
c             fparm(8) = phdif/phmax
c
c          is returned whenever phmax.gt.0.0.  in the degenerate case
c          phmax = 0.0, fparm(8) = phdif is returned.
c
c
c
c ... wk
c
c          on output wk contains intermediate values that must not be
c          destroyed if muh3 is to be called again with iparm(1)=1 or
c          if muh34 is to be called to improve the estimate to fourth
c          order.
c
c ... phi
c
c          on output phi(i,j,k) contains the approximation to
c          p(xi,yj,zk) for all mesh points i=1,...,nx; j=1,...,ny;
c          k=1,...,nz.  the last computed iterate in phi is returned
c          even if convergence is not obtained (ierror=-1)
c
c ... ierror
c
c          for intl=iparm(1)=0 initialization calls, ierror is an
c          error flag that indicates invalid input arguments when
c          returned positive and nonfatal warnings when returned
c          negative.  argument checking and discretization
c          is bypassed for intl=1 calls which can only return
c          ierror = -1 or 0 or 1.
c
c
c     non-fatal warnings * * *
c
c
c     =-5 if kcycle=mgopt(1) is greater than 2. values larger than 2 results
c         in an algorithm which probably does far more computation than
c         necessary.  kcycle = 1 (v cycles)  or kcycle=2 (w cycles) should
c         suffice for most problems.  ierror = -5 is also set if either
c         iprer = mgopt(2) or ipost=mgopt(3) is greater than 2.  the
c         ierror=-5 flag is overridden by any other fatal or non-fatal
c         error.
c
c     =-4 if there are dominant nonzero first order terms in the pde which
c         make it "hyperbolic" at the finest grid level. numerically, this
c         happens if:
c
c              abs(cx)*dlx > 2.*abs(cxx)   (dlx = (xb-xa)/(nx-1))
c
c                         (or)
c
c              abs(cy)*dly > 2.*abs(cyy)   (dly = (yd-yc)/(ny-1))
c
c
c         at some fine grid point (xi,yj).  if an adjustment is not made the
c         condition can lead to a matrix coming from the discretization
c         which is not diagonally dominant and divergence is possible. since
c         the condition is "likely" at coarser grid levels for pde's with
c         nonzero first order terms, the adjustments (actually first order
c         approximations to the pde)
c
c
c             cxx = amax1(cxx,0.5*abs(cx)*dx)
c
c                          (and)
c
c             cyy = amax1(cyy,0.5*abs(cy)*dy)
c
c
c         (here dx,dy are the x,y mesh sizes of the subgrid)
c
c         are made to preserve convergence of multigrid iteration. if made
c         at the finest grid level, it can lead to convergence to an
c         erroneous solution (flagged by ierror = -4).  a possible remedy
c         is to increase resolution. the ierror = -4 flag overrides the
c         nonfatal ierror = -5 flag.
c
c
c     =-3  if the continuous elliptic pde is singular.  this means the
c          boundary conditions are periodic or pure derivative at all
c          boundaries and ce(x,y) = 0.0 for all x,y.  a solution is still
c          attempted but convergence may not occur due to ill-conditioning
c          of the linear system coming from the discretization.  the
c          ierror = -3 flag overrides the ierror=-4 and ierror=-5 nonfatal
c          flags.
c
c
c     =-2  if the pde is not elliptic (i.e., cxx*cyy.le.0.0 for some (xi,yj))
c          in this case a solution is still attempted although convergence
c          may not occur due to ill-conditioning of the linear system.
c          the ierror = -2 flag overrides the ierror=-5,-4,-3 nonfatal
c          flags.
c
c
c     =-1  if convergence to the tolerance specified in tolmax=fparm(5)>0.
c          is not obtained in maxcy=iparm(13) multigrid cycles between the
c          coarsest (ixp+1,jyq+1) and finest (nx,ny) grid levels.
c          in this case the last computed iterate is still returned.
c          the ierror = -1 flag overrides all other nonfatal flags
c
c
c     no errors * * *
c
c     = 0
c
c     fatal argument errors * * *
c
c     = 1 if intl=iparm(1) is not 0 on initial call or not 0 or 1
c         on subsequent calls
c
c     = 2 if any of the boundary condition flags nxa,nxb,nyc,nyd,nze,nzf
c         in iparm(2) through iparm(7)is not 0,1 or 2 or if
c         (nxa,nxb) or (nyc,nyd) or (nze,nzf) are not pairwise zero.
c
c     = 3 if mino(ixp,jyq,kzr) < 2 (ixp=iparm(8),jyq=iparm(9),kzr=iparm(10))
c         or if ixp<3 when nxa=0 or jyq<3 when nyc=0 or kzr<3 when nze=0.
c
c     = 4 if min0(iex,jey,kez) < 1 (iex=iparm(11),jey=iparm(12),kez=iparm(13))
c         or if max0(iex,jey,kez) > 50
c
c     = 5 if nx.ne.ixp*2**(iex-1)+1 or if ny.ne.jyq*2**(jey-1)+1 or
c         if nz.ne.kzr*2**(kez-1)+1 (nx=iparm(14),ny=iparm(15),nz=iparm(16))
c
c     = 6 if iguess = iparm(17) is not equal to 0 or 1
c
c     = 7 if maxcy = iparm(18) < 1 (large values for maxcy should not be used)
c
c     = 8 if method = iparm(19) is less than 0 or greater than 10
c
c     = 9 if length = iparm(20) is too small (see iparm(21) on output
c         for minimum required work space length)
c
c     =10 if xa > xb or yc > yd or ze > zf
c         (xa=fparm(1),xb=fparm(2),yc=fparm(3),yd=fparm(4),ze=fparm(5),zf=fparm(6))
c
c     =11 if tolmax = fparm(7) < 0.0
c
c     errors in setting multigrid options * * * (see also ierror=-5)
c
c     =12 if kcycle = mgopt(1) < 0 or
c         if iprer = mgopt(2) < 1 or
c         if ipost = mgopt(3) < 1 or
c         if intpol = mgopt(4) is not 1 or 3
c
c     =13 if iex=jey=kez=1 (full direct method) and iguess=1 or maxcy > 1
c
c     =14 if the elliptic pde is singular (see ierror=-3 in mud3.d)
c
c
c *********************************************************
c *********************************************************
c
c     end of muh3 documentation
c
c **********************************************************
c **********************************************************
c
c
