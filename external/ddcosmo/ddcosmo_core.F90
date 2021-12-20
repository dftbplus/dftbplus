!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  COPYRIGHT (C) 2015 by Filippo Lipparini, Benjamin Stamm, Paolo Gatto        !
!  Eric Cancès, Yvon Maday, Jean-Philip Piquemal, Louis Lagardère and          !
!  Benedetta Mennucci.                                                         !
!                             ALL RIGHT RESERVED.                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! A modular implementation of COSMO using a domain decomposition linear scaling
! strategy.
!
! This code is governed by the LGPL license and abiding by the rules of
! distribution of free software.
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Lesser General Public License for more details.
module ddcosmo_core
  implicit none
  private

  public :: TDomainDecomposition, TDomainDecompositionInput, TDomainDecomposition_init
  public :: hsnorm, calcv, intrhs, prtsph, adjrhs, wghpot, ddupdate, fdoka, fdokb, fdoga


  integer, parameter :: dp = selected_real_kind(15)
  integer, parameter :: ndiis=25, iout=6, nngmax=100
  real(dp), parameter :: zero=0._dp, pt5=0.5_dp, one=1._dp, two=2._dp, four=4._dp
  real(dp), parameter :: se = -1.0_dp
  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
  real(dp), parameter :: sq2 = sqrt(2.0_dp)


  !> Quantities contained in control file
  type :: TDomainDecompositionInput

    !> Max angular momentum of spherical harmonics basis
    integer :: lmax

    !> Threshold for iterative solver
    real(dp) :: conv

    !> Regularization parameters
    real(dp) :: eta

  end type TDomainDecompositionInput


  !> Cosmo calculator
  type :: TDomainDecomposition

    !> Printing flag
    integer :: iprint

    !> Max angular momentum of spherical harmonics basis
    integer :: lmax

    !> Desired number of Lebedev integration points
    integer :: ngrid

    !> Threshold for iterative solver
    real(dp) :: conv

    !> 1) compute forces ; 0) do not compute forces
    integer :: igrad

    !> Regularization parameters
    real(dp) :: eta

    !> Number of spheres/atoms
    integer :: nat

    !> Number of integration points on cavity's boundary
    integer :: ncav

    !> Number of basis functions, i.e., spherical harmonics
    integer :: nylm

    !> Workspaces
    logical :: grad
    integer, allocatable :: inl(:), nl(:)
    real(dp), allocatable :: rvdw(:), xyz(:, :), ccav(:, :)
    real(dp), allocatable :: w(:), grid(:, :), basis(:, :)
    real(dp), allocatable :: fact(:), facl(:), facs(:)
    real(dp), allocatable :: fi(:, :), ui(:, :), zi(:, :, :)

  end type TDomainDecomposition


contains


  !> Initialize data structure
  !
  !  allocate the various arrays needed for ddcosmo,
  !  assemble the cavity and the various associated geometrical quantities.
  subroutine TDomainDecomposition_init(self, input, rvdw, wang, grid)

    character(len=*), parameter :: source = 'ddcosmo_core::TDomainDecomposition_init'
    type(TDomainDecomposition), intent(out) :: self
    type(TDomainDecompositionInput), intent(in) :: input
    real(dp), intent(in) :: rvdw(:)
    real(dp), intent(in) :: wang(:), grid(:, :)

    integer :: iat, jat, i, ii, lnl, l, ind, m, igrid, inear, jnear
    integer :: istatus
    real(dp) :: fac, fl, ffl, fnorm, d2, r2, v(3), vv, t, xt, swthr

    real(dp), allocatable :: vcos(:), vsin(:), vplm(:)

    self%nat = size(rvdw)
    self%ngrid = size(wang)
    self%iprint = 0
    self%lmax = input%lmax
    self%conv = input%conv
    self%eta = input%eta
    self%igrad = 1

    ! allocate:
    self%grad   = self%igrad /= 0
    self%nylm   = (self%lmax+1)*(self%lmax+1)
    allocate(self%rvdw(self%nat), self%xyz(3, self%nat), self%w(self%ngrid), &
      & self%grid(3, self%ngrid), self%basis(self%nylm, self%ngrid), &
      & self%inl(self%nat+1), self%nl(self%nat*nngmax), &
      & self%fi(self%ngrid, self%nat), self%ui(self%ngrid, self%nat), &
      & self%fact(max(2*self%lmax+1, 2)), self%facl(self%nylm), &
      & self%facs(self%nylm))
    if (self%grad) allocate(self%zi(3, self%ngrid, self%nat))

    ! precompute various quantities (diagonal blocks of ddcosmo matrix,
    ! normalization factors for spherical harmonics, factorials)
    self%fact(1) = 1.0_dp
    self%fact(2) = 1.0_dp
    do i = 3, 2*self%lmax + 1
      self%fact(i) = real(i-1, dp)*self%fact(i-1)
    end do

    do l = 0, self%lmax
      ind = l*l + l + 1
      fl  = (2.0_dp*real(l, dp) + 1.0_dp)/(4.0_dp*pi)
      ffl = sqrt(fl)
      self%facl(ind-l:ind+l) = fl
      self%facs(ind) = ffl
      do m = 1, l
        fnorm = sq2*ffl*sqrt(self%fact(l-m+1)/self%fact(l+m+1))
        if (mod(m, 2) == 1) fnorm = -fnorm
        self%facs(ind+m) = fnorm
        self%facs(ind-m) = fnorm
      end do
    end do

    ! copy the angular grid:
    self%grid = grid
    ! scaling because the weights are normalised
    self%w = wang * four*pi

    ! set the centers and radii of the spheres:
    self%rvdw(:) = rvdw

    ! build a basis of spherical harmonics at the gridpoints:
    allocate(vplm(self%nylm), vcos(self%lmax+1), vsin(self%lmax+1))
    !$omp parallel do default(none) schedule(runtime) &
    !$omp shared(self) private(i, vplm, vcos, vsin)
    do i = 1, self%ngrid
      call ylmbas(self, self%grid(:, i), self%basis(:, i), vplm, vcos, vsin)
    end do
    deallocate(vplm, vcos, vsin)

    if (self%iprint >= 4) then
      call prtsph(self, 'facs', 1, 0, self%facs)
      call prtsph(self, 'facl', 1, 0, self%facl)
      call prtsph(self, 'basis', self%ngrid, 0, self%basis)
      call ptcart(self, 'grid', 3, 0, self%grid)
      call ptcart(self, 'weights', 1, 0, self%w)
    end if

  end subroutine TDomainDecomposition_init


  subroutine ddupdate(self, xyz)

    type(TDomainDecomposition), intent(inout) :: self

    real(dp), intent(in) :: xyz(:, :)

    self%xyz(:, :) = xyz

    call mknnl(self)
    call mkfiui(self)
    call mkccav(self)

  end subroutine ddupdate


  ! build neighbors list (CSR format)
  !
  !  \\  jat  |
  ! iat   \\  |  1   2   3   4   5   6
  ! -----------------------------------
  !         1 |      x       x   x
  !         2 |  x       x       x   x
  !         3 |      x       x       x
  !         4 |  x       x       x   x
  !         5 |  x   x       x
  !         6 |      x   x   x
  !
  !
  !  inl =  1,       4,          8,      11,         15,      18, 21        pointer to 1st neighbor
  !
  !         |        |           |        |           |        |
  !         v        V           V        V           V        V
  !
  !         1| 2| 3| 4| 5| 6| 7| 8| 9|10|11|12|13|14|15|16|17|18|19|20
  !
  !  nl  =  2, 4, 5, 1, 3, 5, 6, 2, 4, 6, 1, 3, 5, 6, 1, 2, 4, 2, 3, 4     neighbors list
  subroutine mknnl(self)

    type(TDomainDecomposition), intent(inout) :: self

    integer :: ii, lnl, iat, jat
    real(dp) :: v(3), d2, r2

    character(len=*), parameter :: f1000 = "(t3, 'neighbours of sphere ', i6)"
    character(len=*), parameter :: f1010 = "(t5, 12i6)"

    ! index of nl
    ii  = 1
    lnl = 0
    do iat = 1, self%nat
      self%inl(iat) = lnl + 1
      do jat = 1, self%nat
        if (iat /= jat) then
          v(:) = self%xyz(:, iat) - self%xyz(:, jat)
          d2 = v(1)**2 + v(2)**2 + v(3)**2
          r2 = (self%rvdw(iat) + self%rvdw(jat))**2
          if (d2 <= r2) then
            self%nl(ii) = jat
            ii  = ii + 1
            lnl = lnl + 1
          end if
        end if
      end do
    end do
    self%inl(self%nat+1) = lnl+1

    if (self%iprint >= 4) then
      write(iout, *) '   inl:'
      write(iout, '(10i6)') self%inl(1:self%nat+1)
      write(iout, *)
      do iat = 1, self%nat
        write(iout, f1000) iat
        write(iout, f1010) self%nl(self%inl(iat):self%inl(iat+1)-1)
      end do
      write(iout, *)
    end if

  end subroutine mknnl


  ! Define :
  !
  !   N_i = list of neighbors of i-sphere [ excluding i-sphere ]
  !
  !            | r_i + rho_i s_n - r_j |
  !   t_n^ij = -------------------------
  !                      rho_j
  !
  !   fi(n, i) =    sum    \chi(t_n^ij)
  !             j \in N_i
  !
  ! Notice that the derivative of fi(n, i) wrt to r_k is (possibly) nonzero
  ! when either k = i, or k \in N_j .
  subroutine mkfiui(self)

    type(TDomainDecomposition), intent(inout) :: self

    integer :: iat, i, ii, jat
    real(dp) :: v(3), vv, t, xt, swthr, fac


    ! build arrays fi, ui, zi
    self%fi = 0.0_dp
    self%ui = 0.0_dp
    if (self%grad) self%zi = 0.0_dp

    !$omp parallel do default(none) schedule(runtime) collapse(2) &
    !$omp shared(self) private(iat, i, ii, jat, v, vv, t, xt, swthr, fac)
    ! loop over spheres
    do iat = 1, self%nat

      ! loop over integration points
      do i = 1, self%ngrid

        ! loop over neighbors of i-sphere
        do ii = self%inl(iat), self%inl(iat+1)-1

          ! neighbor's number
          jat = self%nl(ii)

          ! compute t_n^ij
          v(:) = self%xyz(:, iat) + self%rvdw(iat)*self%grid(:, i) - self%xyz(:, jat)
          vv = v(1)**2 + v(2)**2 + v(3)**2
          vv = sqrt(vv)
          t = vv/self%rvdw(jat)

          ! compute \chi(t_n^ij)
          xt = fsw(t, se, self%eta)

          ! upper bound of switch region
          swthr = 1.0_dp + (se + 1._dp)*self%eta / 2._dp

          ! t_n^ij belongs to switch region
          if (self%grad .and. (t < swthr .and. t > swthr-self%eta)) then
            fac = dfsw(t, se, self%eta) / self%rvdw(jat)

            ! accumulate for zi
            self%zi(:, i, iat) = self%zi(:, i, iat) + fac*v(:)/vv
          end if

          ! accumulate for fi
          self%fi(i, iat) = self%fi(i, iat) + xt
        end do

        ! compute ui
        if (self%fi(i, iat) <= 1.0_dp)  self%ui(i, iat) = 1.0_dp - self%fi(i, iat)
      end do
    end do

    if (self%iprint >= 4) then
      call ptcart(self, 'fi', self%nat, 0, self%fi)
      call ptcart(self, 'ui', self%nat, 0, self%ui)
    end if

  end subroutine mkfiui


  ! build cavity array
  subroutine mkccav(self)

    type(TDomainDecomposition), intent(inout) :: self

    integer :: iat, i, ii
    character(len=*), parameter :: f1100 = "(t3, i8, 3f14.6)"

    if (allocated(self%ccav)) deallocate(self%ccav)

    ! initialize number of cavity points
    self%ncav=0

    ! loop over spheres
    do iat = 1, self%nat

      ! loop over integration points
      do i = 1, self%ngrid

        ! positive contribution from integration point
        if (self%ui(i, iat) > 0.0_dp) then

          ! accumulate
          self%ncav = self%ncav + 1

        end if
      end do
    end do

    ! allocate cavity array
    allocate(self%ccav(3, self%ncav))

    ! initialize cavity array index
    ii = 0

    ! loop over spheres
    do iat = 1, self%nat

      ! loop over integration points
      do i = 1, self%ngrid

        ! positive contribution from integration point
        if (self%ui(i, iat) > 0.0_dp) then

          ! advance cavity array index
          ii = ii + 1

          ! store point
          self%ccav(:, ii) = self%xyz(:, iat) + self%rvdw(iat)*self%grid(:, i)

        end if
      end do
    end do

    if (self%iprint >= 4) then
      write(iout, *) '   external cavity points:'
      do ii = 1, self%ncav
        write(iout, f1100) ii, self%ccav(:, ii)
      end do
      write(iout, *)
    end if

  end subroutine mkccav


  !> Switch function (5th degree polynomial)
  elemental function fsw(t, s, eta)

    real(dp) :: fsw
    real(dp), intent(in) :: t
    real(dp), intent(in) :: s
    real(dp), intent(in) :: eta

    real(dp) :: a, b, flow, x
    real(dp), parameter :: f6=6.0_dp, f10=10._dp, f12=12._dp, f15=15._dp

    ! shift :
    ! s =  0   =>   t - eta/2  [ CENTERED ]
    ! s =  1   =>   t - eta    [ EXTERIOR ]
    ! s = -1   =>   t          [ INTERIOR ]

    ! apply shift
    x = t - (s + 1._dp)*eta / 2._dp

    ! lower bound of switch region
    flow = 1.0_dp - eta

    ! define switch function \chi
    if (x >= 1.0_dp) then
      fsw = 0.0_dp
    elseif (x <= flow) then
      fsw = 1.0_dp
    else
      a = f15*eta - f12
      b = f10*eta*eta - f15*eta + f6
      fsw = ((x-1.0_dp)*(x-1.0_dp)*(1.0_dp-x)*(f6*x*x + a*x + b)) / (eta**5)
    end if

  end function fsw


  !> First derivative of switch function
  elemental function dfsw(t, s, eta)

    real(dp), intent(in) :: t
    real(dp), intent(in) :: s
    real(dp), intent(in) :: eta
    real(dp) :: dfsw

    real(dp) :: flow, x
    real(dp), parameter :: f30=30.0_dp

    ! shift :
    ! s =  0   =>   t - eta/2  [ CENTERED ]
    ! s =  1   =>   t - eta    [ EXTERIOR ]
    ! s = -1   =>   t          [ INTERIOR ]
    !
    ! apply shift
    x = t - (s + 1._dp)*eta / 2._dp

    ! lower bound of switch region
    flow = 1.0_dp - eta

    ! define derivative of switch function \chi
    if (x >= 1.0_dp) then
      dfsw = 0.0_dp
    elseif (x <= flow) then
      dfsw = 0.0_dp
    else
      dfsw = f30*(1.0_dp-x)*(x-1.0_dp)*(x-1.0_dp+eta)*(x-1.0_dp+eta)/(eta**5)
    end if

  end function dfsw


  !> Dump an array (ngrid, ncol) or just a column.
  subroutine ptcart(self, label, ncol, icol, x)

    type(TDomainDecomposition), intent(in) :: self
    character(len=*), intent(in) :: label
    integer, intent(in) :: ncol, icol
    real(dp), intent(in) :: x(self%ngrid, ncol)

    integer :: ig, noff, nprt, ic, j
    character(len=*), parameter :: f1000 = "(1x, i5, f14.8)"
    character(len=*), parameter :: f1010 = "(6x, 5i14)"
    character(len=*), parameter :: f1020 = "(1x, i5, 5f14.8)"

    ! print header :
    if (ncol == 1) then
      write (iout, '(3x, a, 1x, "(column ", i4, ")")') label, icol
    else
      write (iout, '(3x, a)') label
    end if

    ! print entries :
    if (ncol == 1) then
      do ig = 1, self%ngrid
        write(iout, f1000) ig, x(ig, 1)
      end do

    else
      noff = mod (ncol, 5)
      nprt = max(ncol - noff, 0)
      do ic = 1, nprt, 5
        write(iout, f1010) (j, j = ic, ic+4)
        do ig = 1, self%ngrid
          write(iout, f1020) ig, x(ig, ic:ic+4)
        end do
      end do
      write (iout, f1010) (j, j = nprt+1, nprt+noff)
      do ig = 1, self%ngrid
        write(iout, f1020) ig, x(ig, nprt+1:nprt+noff)
      end do
    end if

  end subroutine ptcart


  !> Dump an array (nylm, ncol) or just a column.
  subroutine prtsph(self, label, ncol, icol, x)

    type(TDomainDecomposition), intent(in) :: self
    character (len=*), intent(in) :: label
    integer, intent(in) :: ncol, icol
    real(dp), intent(in) :: x(self%nylm, ncol)

    integer :: l, m, ind, noff, nprt, ic, j
    character(len=*), parameter :: f1000 = "(1x, i3, i4, f14.8)"
    character(len=*), parameter :: f1010 = "(8x, 5i14)"
    character(len=*), parameter :: f1020 = "(1x, i3, i4, 5f14.8)"

    ! print header :
    if (ncol == 1) then
      write (iout, '(3x, a, 1x, "(column ", i4, ")")') label, icol
    else
      write (iout, '(3x, a)') label
    end if

    ! print entries :
    if (ncol == 1) then
      do l = 0, self%lmax
        ind = l*l + l + 1
        do m = -l, l
          write(iout, f1000) l, m, x(ind+m, 1)
        end do
      end do

    else
      noff = mod(ncol, 5)
      nprt = max(ncol - noff, 0)
      do ic = 1, nprt, 5
        write(iout, f1010) (j, j = ic, ic+4)
        do l = 0, self%lmax
          ind = l*l + l + 1
          do m = -l, l
            write(iout, f1020) l, m, x(ind+m, ic:ic+4)
          end do
        end do
      end do
      write (iout, f1010) (j, j = nprt+1, nprt+noff)
      do l = 0, self%lmax
        ind = l*l + l + 1
        do m = -l, l
          write(iout, f1020) l, m, x(ind+m, nprt+1:nprt+noff)
        end do
      end do
    end if

  end subroutine prtsph


  !> Integrate against spherical harmonics
  subroutine intrhs(self, iat, x, xlm)

    type(TDomainDecomposition), intent(in) :: self
    integer, intent(in) :: iat
    real(dp), intent(in) :: x(:) ! [self%ngrid]
    real(dp), intent(inout) :: xlm(:) ! [self%nylm]

    integer :: ig

    ! initialize
    xlm = 0.0_dp

    ! accumulate over integration points
    do ig = 1, self%ngrid
      xlm = xlm + self%basis(:, ig)*self%w(ig)*x(ig)
    end do

    ! printing
    if (self%iprint >= 5) then
      call ptcart(self, 'pot', 1, iat, x)
      call prtsph(self, 'vlm', 1, iat, xlm)
    end if

  end subroutine intrhs


  !> Compute spherical harmonics
  pure subroutine ylmbas(self, x, basloc, vplm, vcos, vsin)

    type(TDomainDecomposition), intent(in) :: self
    real(dp), intent(in) :: x(:) ! [3]
    real(dp), intent(out) :: basloc(:), vplm(:) ! [self%nylm]
    real(dp), intent(out) :: vcos(:), vsin(:) ! [self%lmax+1]

    integer :: l, m, ind
    real(dp) :: cthe, sthe, cphi, sphi, plm

    ! get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
    ! coordinates of x.

    ! evaluate cos(theta) ; sin(theta)
    cthe = x(3)
    sthe = sqrt(1.0_dp - cthe*cthe)

    ! evalutate cos(phi) ; sin(phi)
    if (sthe /= 0.0_dp) then
      cphi = x(1)/sthe
      sphi = x(2)/sthe
    else
      cphi = 0.0_dp
      sphi = 0.0_dp
    end if

    ! evaluate cos(m*phi) and sin(m*phi) arrays. notice that this is
    ! pointless if z = 1, as the only non vanishing terms will be the
    ! ones with m=0.
    if(sthe /= 0.0_dp) then
      call trgev(self, cphi, sphi, vcos, vsin)
    else
      vcos = 1.0_dp
      vsin = 0.0_dp
    end if

    ! evaluate the generalized legendre polynomials
    call polleg(self, cthe, sthe, vplm)

    ! now build the spherical harmonics. we will distinguish m=0,
    ! m>0 and m<0:
    do l = 0, self%lmax
      ind = l**2 + l + 1

      ! m = 0
      basloc(ind) = self%facs(ind)*vplm(ind)

      do m = 1, l

        plm = vplm(ind+m)

        ! m > 0
        basloc(ind+m) = self%facs(ind+m)*plm*vcos(m+1)

        ! m < 0
        basloc(ind-m) = self%facs(ind-m)*plm*vsin(m+1)

      end do
    end do

  end subroutine ylmbas


  !> Compute first derivatives of spherical harmonics
  pure subroutine dbasis(self, x, basloc, dbsloc, vplm, vcos, vsin)

    type(TDomainDecomposition), intent(in) :: self
    real(dp), intent(in) :: x(:) ! [3]
    real(dp), intent(inout) :: basloc(:), vplm(:) ! [self%nylm]
    real(dp), intent(inout) :: dbsloc(:, :) ! [3, self%nylm]
    real(dp), intent(inout) :: vcos(:), vsin(:) ! [self%lmax+1]

    integer :: l, m, ind, VC, VS
    real(dp) :: cthe, sthe, cphi, sphi, plm, fln, pp1, pm1, pp
    real(dp) :: et(3), ep(3)

    ! get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
    ! coordinates of x.
    cthe = x(3)
    sthe = sqrt(1.0_dp - cthe*cthe)

    if (sthe /= 0.0_dp) then
      ! not (NORTH or SOUTH pole)
      cphi = x(1)/sthe
      sphi = x(2)/sthe
    else
      ! NORTH or SOUTH pole
      cphi = 1.0_dp
      sphi = 0.0_dp
    end if

    ! evaluate the derivatives of theta and phi:
    et(1) = cthe*cphi
    et(2) = cthe*sphi
    et(3) = -sthe

    if(sthe /= 0.0_dp) then
      ! not (NORTH or SOUTH pole)
      ep(1) = -sphi/sthe
      ep(2) = cphi/sthe
      ep(3) = 0.0_dp
    else
      ! NORTH or SOUTH pole
      ep(1) = 0.0_dp
      ep(2) = 1.0_dp
      ep(3) = 0.0_dp
    end if

    ! evaluate cos(m*phi) and sin(m*phi) arrays. notice that this is
    ! pointless if z = 1, as the only non vanishing terms will be the
    ! ones with m=0.

    if (sthe /= 0.0_dp) then
      ! not (NORTH or SOUTH pole)
      call trgev(self, cphi, sphi, vcos, vsin)
    else
      ! NORTH or SOUTH pole
      vcos(:) = 1.0_dp
      vsin(:) = 0.0_dp
    end if
    VC=0.0_dp
    VS=cthe

    ! evaluate the generalized legendre polynomials.
    call polleg(self, cthe, sthe, vplm)

    ! now build the spherical harmonics. we will distinguish m=0,
    ! m>0 and m<0:

    basloc = 0.0_dp
    dbsloc = 0.0_dp
    do l = 0, self%lmax
      ind = l*l + l + 1
      ! m = 0
      fln = self%facs(ind)
      basloc(ind) = fln*vplm(ind)
      if (l > 0) then
        dbsloc(:, ind) = fln*vplm(ind+1)*et(:)
      else
        dbsloc(:, ind) = 0.0_dp
      end if
      !$omp simd
      do m = 1, l
        fln = self%facs(ind+m)
        plm = fln*vplm(ind+m)
        pp1 = 0.0_dp
        if (m < l) pp1 = -0.5_dp*vplm(ind+m+1)
        pm1 = 0.5_dp*(real(l+m, dp)*real(l-m+1, dp)*vplm(ind+m-1))
        pp  = pp1 + pm1

        ! m > 0
        basloc(ind+m) = plm*vcos(m+1)

        if (sthe /= 0.0_dp) then
          ! not (NORTH or SOUTH pole)
          dbsloc(:, ind+m) = -fln*pp*vcos(m+1)*et(:) - real(m, dp)*plm*vsin(m+1)*ep(:)
        else
          ! NORTH or SOUTH pole
          dbsloc(:, ind+m) = -fln*pp*vcos(m+1)*et(:) - fln*pp*ep(:)*VC
        end if

        ! m < 0
        basloc(ind-m) = plm*vsin(m+1)

        if (sthe /= 0.0_dp) then
          ! not (NORTH or SOUTH pole)
          dbsloc(:, ind-m) = -fln*pp*vsin(m+1)*et(:) + real(m, dp)*plm*vcos(m+1)*ep(:)
        else
          ! NORTH or SOUTH pole
          dbsloc(:, ind-m) = -fln*pp*vsin(m+1)*et(:) - fln*pp*ep(:)*VS
        end if

      end do
    end do

  end subroutine dbasis


  !> Compute the l, m associated legendre polynomial for -1 <= x <= 1 using the
  !  recurrence formula
  !
  !  (l-m)p(l, m) = x(2l-1)p(l-1, m) - (l+m-1)p(l-2, m)
  pure subroutine polleg(self, x, y, plm)

    type(TDomainDecomposition), intent(in) :: self
    real(dp), intent(in) :: x, y
    real(dp), intent(inout) :: plm(:) ! [self%nylm]

    integer :: m, ind, l, ind2
    real(dp) :: fact, pmm, somx2, pmm1, pmmo, pll, fm, fl

    fact  = 1.0_dp
    pmm   = 1.0_dp
    somx2 = y

    do m = 0, self%lmax
      ind      = (m + 1)*(m + 1)
      plm(ind) = pmm
      if(m == self%lmax) return
      fm = real(m, dp)
      pmm1 = x*(2.0_dp*fm + 1.0_dp)*pmm
      ind2 = ind + 2*m + 2
      plm(ind2) = pmm1
      pmmo = pmm
      do l = m+2, self%lmax
        fl = real(l, dp)
        pll   = (x*(2.0_dp*fl - 1.0_dp)*pmm1 - (fl + fm - 1.0_dp)*pmm)/(fl - fm)
        ind = l*l + l + 1
        plm(ind+m) = pll
        pmm  = pmm1
        pmm1 = pll
      end do
      pmm  = -pmmo*fact*somx2
      fact = fact + 2.0_dp
    end do

  end subroutine polleg


  !> Service routine for computation of spherical harmonics
  pure subroutine trgev(self, x, y, cx, sx)

    type(TDomainDecomposition), intent(in) :: self
    real(dp), intent(in) :: x, y
    real(dp), intent(inout) :: cx(:), sx(:) ! [max((self%lmax+1), 2)]

    integer :: m

    cx(1) = 1.0_dp
    sx(1) = 0.0_dp
    cx(2) = x
    sx(2) = y
    do m = 3, self%lmax+1
      cx(m) = 2.0_dp*x*cx(m-1) - cx(m-2)
      sx(m) = 2.0_dp*x*sx(m-1) - sx(m-2)
    end do

  end subroutine trgev


  !> Compute
  !
  !  sum   4pi/(2l+1) t^l * Y_l^m(s) * sigma_l^m
  !  l, m
  !
  !  which is need to compute action of COSMO matrix L.
  pure function intmlp(self, t, sigma, basloc)

    type(TDomainDecomposition), intent(in) :: self
    real(dp), intent(in) :: t
    real(dp), intent(in) :: sigma(:), basloc(:) ! [self%nylm]
    real(dp) :: intmlp

    integer :: l, ind
    real(dp) :: tt, ss, fac

    ! initialize t^l
    tt = 1.0_dp

    ! initialize
    ss = 0.0_dp

    ! loop over l
    do l = 0, self%lmax

      ind = l*l + l + 1

      ! update factor 4pi / (2l+1) * t^l
      fac = tt / self%facl(ind)

      ! contract over l, m and accumulate
      ss = ss + fac * dot_product(basloc(ind-l:ind+l), sigma(ind-l:ind+l))

      ! update t^l
      tt = tt*t
    end do

    ! redirect
    intmlp = ss

  end function intmlp


  !> Weigh potential at cavity points by characteristic function "ui"
  subroutine wghpot(self, phi, g)

    type(TDomainDecomposition), intent(in) :: self
    real(dp), intent(in)  :: phi(:) ! [self%ncav]
    real(dp), intent(out) :: g(:, :) ! [self%ngrid, self%nat]

    integer :: iat, ig, ic

    !> Initialize
    ic = 0
    g(:, :) = 0._dp

    ! loop over spheres
    do iat = 1, self%nat

      ! loop over points
      do ig = 1, self%ngrid

        ! nonzero contribution from point
        if (self%ui(ig, iat) /= 0.0_dp) then

          ! advance cavity point counter
          ic = ic + 1

          ! weight by (negative) characteristic function
          g(ig, iat) = -self%ui(ig, iat) * phi(ic)

        end if

      end do
    end do

  end subroutine wghpot


  !> Compute H-norm
  pure subroutine hsnorm(self, u, unorm)

    type(TDomainDecomposition), intent(in) :: self
    real(dp), intent(in) :: u(:) !< [nylm]
    real(dp), intent(inout) :: unorm

    integer :: l, m, ind
    real(dp) :: fac

    ! initialize
    unorm = 0.0_dp

    ! loop over l
    do l = 0, self%lmax

      ! first index associated to l
      ind = l*l + l + 1

      ! scaling factor
      fac = 1.0_dp/(1.0_dp + real(l, dp))

      ! loop over m
      do m = -l, l

        ! accumulate
        unorm = unorm + fac*u(ind+m)*u(ind+m)

      end do
    end do

    ! the much neglected square root
    unorm = sqrt(unorm)

  end subroutine hsnorm


  !> Compute
  !
  !   v_l^m = v_l^m +
  !
  !               4 pi           l
  !     sum  sum  ---- (t_n^ji)  Y_l^m(s_n^ji) W_n^ji [ \xi_j ]_n
  !      j    n   2l+1
  !
  ! which is related to the action of the adjont COSMO matrix L^* in the following
  ! way. Set
  !
  !   [ \xi_j ]_n = sum  Y_l^m(s_n) [ s_j ]_l^m
  !                 l, m
  !
  ! then
  !
  !   v_l^m = -   sum    (L^*)_ij s_j
  !             j \ne i
  !
  ! The auxiliary quantity [ \xi_j ]_l^m needs to be computed explicitly.
  pure subroutine adjrhs(self, iat, xi, vlm, basloc, vplm, vcos, vsin)

    type(TDomainDecomposition), intent(in) :: self
    integer, intent(in) :: iat
    real(dp), intent(in) :: xi(:, :) ! [self%ngrid, self%nat]
    real(dp), intent(inout) :: vlm(:) ! [self%nylm]
    real(dp), intent(inout) :: basloc(:), vplm(:) ! [self%nylm]
    real(dp), intent(inout) :: vcos(:), vsin(:) ! [self%lmax+1]

    integer :: ij, jat, ig, l, ind, m
    real(dp) :: vji(3), vvji, tji, sji(3), xji, oji, fac, ffac, t

    ! loop over neighbors of i-sphere
    do ij = self%inl(iat), self%inl(iat+1)-1

      ! j-sphere is neighbor
      jat = self%nl(ij)

      ! loop over integration points
      do ig = 1, self%ngrid

        ! compute t_n^ji = | r_j + \rho_j s_n - r_i | / \rho_i
        vji  = self%xyz(:, jat) + self%rvdw(jat)*self%grid(:, ig) - self%xyz(:, iat)
        vvji = vji(1)**2 + vji(2)**2 + vji(3)**2
        vvji = sqrt(vvji)
        tji  = vvji/self%rvdw(iat)

        ! point is INSIDE i-sphere (+ transition layer)
        if (tji < (1.0_dp + (se+1.0_dp)/2.0_dp*self%eta)) then

          ! compute s_n^ji
          sji = vji/vvji

          ! compute \chi(t_n^ji)
          xji = fsw(tji, se, self%eta)

          ! compute W_n^ji
          if (self%fi(ig, jat) > 1.0_dp) then
            oji = xji/self%fi(ig, jat)
          else
            oji = xji
          end if

          ! compute Y_l^m(s_n^ji)
          call ylmbas(self, sji, basloc, vplm, vcos, vsin)

          ! initialize (t_n^ji)^l
          t = 1.0_dp

          ! compute w_n * xi(n, j) * W_n^ji
          fac = self%w(ig) * xi(ig, jat) * oji

          ! loop over l
          do l = 0, self%lmax
            ind  = l*l + l + 1

            ! compute 4pi / (2l+1) * (t_n^ji)^l * w_n * xi(n, j) * W_n^ji
            ffac = fac*t/self%facl(ind)

            ! loop over m
            do m = -l, l
              vlm(ind+m) = vlm(ind+m) + ffac*basloc(ind+m)
            end do

            ! update (t_n^ji)^l
            t = t*tji
          end do

        end if
      end do
    end do

  end subroutine adjrhs


  subroutine header(self)

    type(TDomainDecomposition), intent(in) :: self

    1000 format(/, &
      ' An implementation of COSMO using a domain decomposition linear scaling strategy.', /)
    1010 format(' Parameters:', /, &
      '   number of grid points:                  ', 8x, i8, /,   &
      '   number of spheres:                      ', 8x, i8, /,   &
      '   lmax for the spherical harmonics basis: ', 8x, i8, /,   &
      '   convergence threshold:                  ', 8x, d8.1, /, &
      '   regularization parameter (eta):         ', 8x, f8.3/)

    if (self%iprint > 0) then

      write(iout, 1000)
      write(iout, 1010) self%ngrid, self%nat, self%lmax, self%conv, self%eta

      if (self%igrad == 1)  write(iout, 1013)
      1013   format(' Compute forces.'//)

    end if

  end subroutine header


  !> Compute the first part of <S, L^(x)X>
  pure subroutine fdoka(self, iat, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, fx)

    type(TDomainDecomposition), intent(in) :: self
    integer, intent(in) :: iat
    real(dp), intent(in) :: sigma(:, :) ! [self%nylm, self%nat]
    real(dp), intent(in) :: xi(:) ! [self%ngrid]
    real(dp), intent(inout) :: basloc(:), vplm(:) ! [self%nylm]
    real(dp), intent(inout) :: dbsloc(:, :) ! [3, self%nylm]
    real(dp), intent(inout) :: vcos(:), vsin(:) ! [self%lmax+1]
    real(dp), intent(inout) :: fx(:) ! [3]

    integer :: ig, ij, jat, l, ind, m
    real(dp) :: vvij, tij, xij, oij, t, fac, fl, f1, f2, f3, beta, tlow, thigh
    real(dp) :: vij(3), sij(3), alp(3), va(3)

    tlow  = 1.0_dp - 0.5_dp*(1.0_dp - se)*self%eta
    thigh = 1.0_dp + 0.5_dp*(1.0_dp + se)*self%eta

    do ig = 1, self%ngrid
      va = 0.0_dp
      do ij = self%inl(iat), self%inl(iat+1) - 1
        jat = self%nl(ij)
        vij = self%xyz(:, iat) + self%rvdw(iat)*self%grid(:, ig) - self%xyz(:, jat)
        vvij = vij(1)**2 + vij(2)**2 + vij(3)**2
        vvij = sqrt(vvij)
        tij = vvij/self%rvdw(jat)

        if (tij >= thigh) cycle

        sij = vij/vvij
        call dbasis(self, sij, basloc, dbsloc, vplm, vcos, vsin)
        alp = 0.0_dp
        t = 1.0_dp
        do l = 1, self%lmax
          ind = l*l + l + 1
          fl = real(l, dp)
          fac = t/self%facl(ind)
          do m = -l, l
            f2 = fac*sigma(ind+m, jat)
            f1 = f2*fl*basloc(ind+m)
            alp(:) = alp(:) + f1*sij(:) + f2*dbsloc(:, ind+m)
          end do
          t = t*tij
        end do
        beta = intmlp(self, tij, sigma(:, jat), basloc)
        xij = fsw(tij, se, self%eta)
        if (self%fi(ig, iat) > 1.0_dp) then
          oij = xij/self%fi(ig, iat)
          f2  = -oij/self%fi(ig, iat)
        else
          oij = xij
          f2  = 0.0_dp
        end if
        f1 = oij/self%rvdw(jat)
        va(:) = va(:) + f1*alp(:) + beta*f2*self%zi(:, ig, iat)
        if (tij > tlow) then
          f3 = beta*dfsw(tij, se, self%eta)/self%rvdw(jat)
          if (self%fi(ig, iat) > 1.0_dp) f3 = f3/self%fi(ig, iat)
          va(:) = va(:) + f3*sij(:)
        end if
      end do
      fx = fx - self%w(ig)*xi(ig)*va(:)
    end do

  end subroutine fdoka


  !> Compute the the second part of <S, L^(x)X>
  pure subroutine fdokb(self, iat, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, fx)

    type(TDomainDecomposition), intent(in) :: self
    integer, intent(in) :: iat
    real(dp), intent(in) :: sigma(:, :) ! [self%nylm, self%nat]
    real(dp), intent(in) :: xi(:, :) ! [self%ngrid, self%nat]
    real(dp), intent(inout) :: basloc(:), vplm(:) ! [self%nylm]
    real(dp), intent(inout) :: dbsloc(:, :) ! [3, self%nylm]
    real(dp), intent(inout) :: vcos(:), vsin(:) ! [self%lmax+1]
    real(dp), intent(inout) :: fx(:) ! [3]

    integer :: ig, ji, jat, l, ind, m, jk, kat
    logical :: proc
    real(dp) :: vvji, tji, xji, oji, t, fac, fl, f1, f2, beta, di, tlow, thigh
    real(dp) :: b, g1, g2, vvjk, tjk, f, xjk
    real(dp) :: vji(3), sji(3), alp(3), vb(3), vjk(3), sjk(3), vc(3)



    tlow  = 1.0_dp - 0.5_dp*(1.0_dp - se)*self%eta
    thigh = 1.0_dp + 0.5_dp*(1.0_dp + se)*self%eta

    do ig = 1, self%ngrid
      vb = 0.0_dp
      vc = 0.0_dp
      do ji = self%inl(iat), self%inl(iat+1) - 1
        jat = self%nl(ji)
        vji = self%xyz(:, jat) + self%rvdw(jat)*self%grid(:, ig) - self%xyz(:, iat)
        vvji = vji(1)**2 + vji(2)**2 + vji(3)**2
        vvji = sqrt(vvji)
        tji = vvji/self%rvdw(iat)

        if (tji > thigh) cycle

        sji = vji/vvji
        call dbasis(self, sji, basloc, dbsloc, vplm, vcos, vsin)

        alp = 0.0_dp
        t   = 1.0_dp
        do l = 1, self%lmax
          ind = l*l + l + 1
          fl = real(l, dp)
          fac = t/self%facl(ind)
          do m = -l, l
            f2 = fac*sigma(ind+m, iat)
            f1 = f2*fl*basloc(ind+m)
            alp = alp + f1*sji + f2*dbsloc(:, ind+m)
          end do
          t = t*tji
        end do
        xji = fsw(tji, se, self%eta)
        if (self%fi(ig, jat) > 1.0_dp) then
          oji = xji/self%fi(ig, jat)
        else
          oji = xji
        end if
        f1 = oji/self%rvdw(iat)
        vb = vb + f1*alp*xi(ig, jat)
        if (tji > tlow) then
          beta = intmlp(self, tji, sigma(:, iat), basloc)
          if (self%fi(ig, jat) > 1.0_dp) then
            di  = 1.0_dp/self%fi(ig, jat)
            fac = di*xji
            proc = .false.
            b    = 0.0_dp
            do jk = self%inl(jat), self%inl(jat+1) - 1
              kat = self%nl(jk)
              vjk = self%xyz(:, jat) + self%rvdw(jat)*self%grid(:, ig) - self%xyz(:, kat)
              vvjk = vjk(1)**2 + vjk(2)**2 + vjk(3)**2
              vvjk = sqrt(vvjk)
              tjk = vvjk/self%rvdw(kat)
              if (kat /= iat) then
                if (tjk <= thigh) then
                  proc = .true.
                  sjk = vjk/vvjk
                  call ylmbas(self, sjk, basloc, vplm, vcos, vsin)
                  g1 = intmlp(self, tjk, sigma(:, kat), basloc)
                  xjk = fsw(tjk, se, self%eta)
                  b = b + g1*xjk
                end if
              end if
            end do
            if (proc) then
              g1 = di*di*dfsw(tji, se, self%eta)/self%rvdw(iat)
              g2 = g1*xi(ig, jat)*b
              vc = vc + g2*sji
            end if
          else
            di = 1.0_dp
            fac = 0.0_dp
          end if
          f2 = (1.0_dp-fac)*di*dfsw(tji, se, self%eta)/self%rvdw(iat)
          vb = vb + f2*xi(ig, jat)*beta*sji
        end if
      end do
      fx = fx + self%w(ig)*(vb - vc)
    end do

  end subroutine fdokb


  !> Compute the U^(x)\Phi contribution to <S, g^(x)>
  pure subroutine fdoga(self, iat, xi, phi, fx)

    type(TDomainDecomposition), intent(in) :: self
    integer, intent(in) :: iat
    real(dp), intent(in) :: xi(:, :), phi(:, :) ! [self%ngrid, self%nat]
    real(dp), intent(inout) :: fx(:) ! [3]

    integer :: ig, ji, jat
    real(dp) :: vvji, tji, fac, swthr
    real(dp) :: alp(3), vji(3), sji(3)

    do ig = 1, self%ngrid
      alp = 0.0_dp
      if (self%ui(ig, iat) > 0.0_dp .and. self%ui(ig, iat) <  1.0_dp) then
        alp = alp + phi(ig, iat)*xi(ig, iat)*self%zi(:, ig, iat)
      end if
      do ji = self%inl(iat), self%inl(iat+1) - 1
        jat = self%nl(ji)
        vji = self%xyz(:, jat) + self%rvdw(jat)*self%grid(:, ig) - self%xyz(:, iat)
        vvji = vji(1)**2 + vji(2)**2 + vji(3)**2
        vvji = sqrt(vvji)
        tji = vvji/self%rvdw(iat)
        swthr = 1.0_dp + (se + 1._dp)*self%eta / 2._dp
        if (tji < swthr .and. tji > swthr-self%eta .and. self%ui(ig, jat) > 0.0_dp) then
          sji = vji/vvji
          fac = - dfsw(tji, se, self%eta)/self%rvdw(iat)
          alp = alp + fac*phi(ig, jat)*xi(ig, jat)*sji
        end if
      end do
      fx = fx - self%w(ig)*alp
    end do

  end subroutine fdoga


  !> Auxiliary routine for COSMO action
  !  compute
  !
  !   \Phi(n) =
  !
  !                       4 pi           l
  !     sum  W_n^ij  sum  ---- (t_n^ij)  Y_l^m(s_n^ij) [ \sigma_j ]_l^m
  !      j           l, m  2l+1
  !
  ! which is related to the action of the COSMO matrix L in the following
  ! way :
  !
  !   -   sum    L_ij \sigma_j = sum  w_n Y_l^m(s_n) \Phi(n)
  !     j \ne i                   n
  !
  ! This second step is performed by routine "intrhs".
  pure subroutine calcv(self, first, iat, pot, sigma, basloc, vplm, vcos, vsin)

    type(TDomainDecomposition), intent(in) :: self
    logical, intent(in) :: first
    integer, intent(in) :: iat
    real(dp), intent(in) :: sigma(:, :) ! [self%nylm, self%nat]
    real(dp), intent(inout) :: pot(:) ! [self%ngrid]
    real(dp), intent(inout) :: basloc(:) ! [self%nylm]
    real(dp), intent(inout) :: vplm(:) ! [self%nylm]
    real(dp), intent(inout) :: vcos(:) ! [self%lmax+1]
    real(dp), intent(inout) :: vsin(:) ! [self%lmax+1]

    integer :: its, ij, jat
    real(dp) :: vij(3), sij(3)
    real(dp) :: vvij2, vvij, tij, xij, oij, stslm, stslm2, stslm3

    ! initialize
    pot(:) = 0.0_dp

    ! if 1st iteration of Jacobi method, then done!
    if (first)  return

    ! loop over grid points
    do its = 1, self%ngrid

      ! contribution from integration point present
      if (self%ui(its, iat) < 1.0_dp) then

        ! loop over neighbors of i-sphere
        do ij = self%inl(iat), self%inl(iat+1)-1

          ! neighbor is j-sphere
          jat = self%nl(ij)

          ! compute t_n^ij = | r_i + \rho_i s_n - r_j | / \rho_j
          vij  = self%xyz(:, iat) + self%rvdw(iat)*self%grid(:, its) - self%xyz(:, jat)
          vvij2 = vij(1)**2 + vij(2)**2 + vij(3)**2
          vvij = sqrt(vvij2)
          tij  = vvij / self%rvdw(jat)

          ! point is INSIDE j-sphere
          if (tij < 1.0_dp) then

            ! compute s_n^ij = (r_i + \rho_i s_n - r_j) / | ... |
            sij = vij / vvij

            ! compute \chi(t_n^ij)
            xij = fsw(tij, se, self%eta)

            ! compute W_n^ij
            if (self%fi(its, iat) > 1.0_dp) then

              oij = xij / self%fi(its, iat)

            else

              oij = xij

            end if

            ! compute Y_l^m(s_n^ij)
            call ylmbas(self, sij, basloc, vplm, vcos, vsin)

            ! accumulate over j, l, m
            pot(its) = pot(its) + oij * intmlp(self, tij, sigma(:, jat), basloc)

          end if
        end do
      end if
    end do

  end subroutine calcv


end module ddcosmo_core
