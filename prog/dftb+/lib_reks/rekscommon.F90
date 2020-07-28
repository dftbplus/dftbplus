!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> REKS and SI-SA-REKS formulation in DFTB as developed by Lee et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet appart from Gamma point.
!> * Orbital potentials or spin-orbit or external E-field does not work yet.
!> * Only for closed shell system.
!> * Onsite corrections are not included in this version
module dftbp_rekscommon

  use dftbp_accuracy
  use dftbp_blasroutines, only : gemm
  use dftbp_densedescr
  use dftbp_message
  use dftbp_reksvar, only : TReksCalc, reksTypes

  implicit none

  private

  public :: checkGammaPoint
  public :: getTwoIndices
  public :: qm2udL, ud2qmL
  public :: qmExpandL, udExpandL
  public :: matAO2MO, matMO2AO
  public :: getSpaceSym
  public :: findShellOfAO
  public :: assignIndex, assignEpsilon, assignFilling

  !> Swap from charge/magnetisation to up/down in REKS
  !> It converts my_qm to my_ud in the representation
  interface qm2udL
    module procedure qm2ud2
    module procedure qm2ud3
  end interface qm2udL

  !> Swap from up/down to charge/magnetisation in REKS
  !> It converts my_ud to my_qm in the representation
  interface ud2qmL
    module procedure ud2qm2
    module procedure ud2qm3
  end interface ud2qmL

  !> Set correct charge/magnetization in REKS
  !> It converts my_qm to qm in the representation
  interface qmExpandL
    module procedure qmExpand3
    module procedure qmExpand4
  end interface qmExpandL

  !> Set correct up/down in REKS
  !> It converts my_ud to ud in the representation
  interface udExpandL
    module procedure udExpand3
    module procedure udExpand4
  end interface udExpandL

  contains

  !> Check whether the cell size is proper to the Gamma point
  !> calculation or not, and set several convenient variables
  subroutine checkGammaPoint(denseDesc, iNeighbour, nNeighbourSK,&
      & iPair, img2CentCell, over, this)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Map from images of atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    integer :: mu, nu, nAtom, nOrb, nAtomSparse
    integer :: iAtom1, iAtom2, iAtom2f, iNeigh1, iOrig1
    integer :: nOrb1, nOrb2, ii, jj, kk, ll

    nAtom = size(denseDesc%iAtomStart,dim=1) - 1
    nOrb = size(this%overSqr,dim=1)

    nAtomSparse = 0
    do iAtom1 = 1, nAtom ! mu
      nAtomSparse = nAtomSparse + nNeighbourSK(iAtom1) + 1
    end do

    deallocate(this%getDenseAtom)
    allocate(this%getDenseAtom(nAtomSparse,2))

    ll = 1
    this%getDenseAO(:,:) = 0
    this%getDenseAtom(:,:) = 0
    do iAtom1 = 1, nAtom ! mu in A atom
      ii = denseDesc%iAtomStart(iAtom1)
      nOrb1 = denseDesc%iAtomStart(iAtom1 + 1) - ii
      do iNeigh1 = 0, nNeighbourSK(iAtom1) ! nu in B atom
        iOrig1 = iPair(iNeigh1,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh1, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = denseDesc%iAtomStart(iAtom2f)
        nOrb2 = denseDesc%iAtomStart(iAtom2f + 1) - jj

        do kk = 1, nOrb1*nOrb2 ! mu and nu
          if (nOrb1 == 1 .or. nOrb2 == 1) then
            ! Read vertical
            mu = denseDesc%iAtomStart(iAtom1) - 1 + mod(kk,nOrb1)
            if (mod(kk,nOrb1) == 0) mu = mu + nOrb1
            nu = denseDesc%iAtomStart(iAtom2f) + kk/nOrb1
            if (mod(kk,nOrb1) == 0) nu = nu - 1
          else
            ! Read horizontal
            mu = denseDesc%iAtomStart(iAtom1) + kk/nOrb2
            if (mod(kk,nOrb2) == 0) mu = mu - 1
            nu = denseDesc%iAtomStart(iAtom2f) - 1 + mod(kk,nOrb2)
            if (mod(kk,nOrb2) == 0) nu = nu + nOrb2
          end if
          ! Find inconsistent index between dense and sparse
          ! It means that current lattice is not proper to Gamma point calculation
          ! TODO : add the condition of Gamma point using nKpoint and Kpoints?
          if (abs(this%overSqr(mu,nu)-over(iOrig1+kk-1)) >= epsilon(1.0_dp)) then
            call error("Inconsistent maching exists between sparse and dense")
          end if
          this%getDenseAO(iOrig1+kk-1,1) = mu
          this%getDenseAO(iOrig1+kk-1,2) = nu
        end do

        this%getDenseAtom(ll,1) = iAtom1  ! A atom
        this%getDenseAtom(ll,2) = iAtom2f ! B atom
        ll = ll + 1
      end do
    end do

    do mu = 1, nOrb
      do iAtom1 = 1, nAtom
        if (mu > denseDesc%iAtomStart(iAtom1)-1 .and.&
            & mu <= denseDesc%iAtomStart(iAtom1+1)-1) then
          this%getAtomIndex(mu) = iAtom1
        end if
      end do
    end do

  end subroutine checkGammaPoint


  !> Calculate two indices from single index with size of index
  subroutine getTwoIndices(sizeInd, curInd, ind1, ind2, option)

    !> size of index
    integer, intent(in) :: sizeInd

    !> current index (single index)
    integer, intent(in) :: curInd

    !> generated two indices
    integer, intent(inout) :: ind1, ind2

    !> option for generating two indices
    !> 1: exclude diagonal elements, 2: include diagonal elements
    integer, intent(in) :: option

    real(dp) :: tmp

    if (option == 1) then

      ! when sizeInd = 3, (ind1,ind2) = (1,2) (1,3) (2,3)
      tmp = ( real(2.0_dp*sizeInd+1.0_dp, dp) - sqrt( (2.0_dp*sizeInd+ &
          & 1.0_dp)**2.0_dp - 8.0_dp*(sizeInd+curInd) ) )/2.0_dp
      ind1 = int( real(tmp, dp) )
      if ( real(tmp, dp)-real(ind1, dp) <= epsilon(1.0_dp) ) then
        ind1 = ind1 - 1
      end if
      ind2 = ind1**2/2 + ind1/2 - sizeInd*ind1 + sizeInd + curInd
      if (mod(ind1,2) == 1) then
        ind2 = ind2 + 1
      end if

    else if (option == 2) then

      ! when sizeInd = 3, (ind1,ind2) = (1,1) (1,2) (1,3) (2,2) (2,3) (3,3)
      tmp = ( real(2.0_dp*sizeInd+3.0_dp, dp) - sqrt( (2.0_dp*sizeInd+ &
          & 3.0_dp)**2.0_dp - 8.0_dp*(sizeInd+curInd) ) )/2.0_dp
      ind1 = int( real(tmp, dp) )
      ind2 = ind1**2/2 - ind1/2 - sizeInd*ind1 + sizeInd + curInd

    end if

  end subroutine getTwoIndices


  !> Converts charge/magnetization set into up/down in REKS
  subroutine qm2ud2(x, Lpaired)

    !> Array of data, second index spin, third index Lmax
    real(dp), intent(inout) :: x(:,:,:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    integer :: iL, Lmax

    Lmax = size(x,dim=3)

    do iL = 1, Lmax
      if (iL <= Lpaired) then
        x(:,1,iL) = x(:,1,iL) * 0.5_dp
      else
        if (mod(iL,2) == 1) then
          x(:,1,iL) = (x(:,1,iL) + x(:,1,iL+1)) * 0.5_dp
        else
          x(:,1,iL) = x(:,1,iL-1) - x(:,1,iL)
        end if
      end if
    end do

  end subroutine qm2ud2


  !> Converts charge/magnetization set into up/down in REKS
  subroutine qm2ud3(x, Lpaired)

    !> Array of data, third index spin, fourth index Lmax
    real(dp), intent(inout) :: x(:,:,:,:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    integer :: iL, Lmax

    Lmax = size(x,dim=4)

    do iL = 1, Lmax
      if (iL <= Lpaired) then
        x(:,:,1,iL) = x(:,:,1,iL) * 0.5_dp
      else
        if (mod(iL,2) == 1) then
          x(:,:,1,iL) = (x(:,:,1,iL) + x(:,:,1,iL+1)) * 0.5_dp
        else
          x(:,:,1,iL) = x(:,:,1,iL-1) - x(:,:,1,iL)
        end if
      end if
    end do

  end subroutine qm2ud3


  !> Converts up/down set into charge/magnetization in REKS
  subroutine ud2qm2(x, Lpaired)

    !> Array of data, second index spin, third index Lmax
    real(dp), intent(inout) :: x(:,:,:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    integer :: iL, Lmax

    Lmax = size(x,dim=3)

    do iL = 1, Lmax
      if (iL <= Lpaired) then
        x(:,1,iL) = x(:,1,iL) * 2.0_dp
      else
        if (mod(iL,2) == 1) then
          x(:,1,iL) = x(:,1,iL) + x(:,1,iL+1)
        else
          x(:,1,iL) = x(:,1,iL-1) - 2.0_dp * x(:,1,iL)
        end if
      end if
    end do

  end subroutine ud2qm2


  !> Converts up/down set into charge/magnetization in REKS
  subroutine ud2qm3(x, Lpaired)

    !> Array of data, third index spin, fourth index Lmax
    real(dp), intent(inout) :: x(:,:,:,:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    integer :: iL, Lmax

    Lmax = size(x,dim=4)

    do iL = 1, Lmax
      if (iL <= Lpaired) then
        x(:,:,1,iL) = x(:,:,1,iL) * 2.0_dp
      else
        if (mod(iL,2) == 1) then
          x(:,:,1,iL) = x(:,:,1,iL) + x(:,:,1,iL+1)
        else
          x(:,:,1,iL) = x(:,:,1,iL-1) - 2.0_dp * x(:,:,1,iL)
        end if
      end if
    end do

  end subroutine ud2qm3


  !> Decide a correct charge/magnetization in REKS
  subroutine qmExpand3(x, Lpaired)

    !> Array of data, third index spin, fourth index Lmax
    real(dp), intent(inout) :: x(:,:,:,:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    integer :: iL, Lmax

    Lmax = size(x,dim=4)

    do iL = 1, Lmax
      if (iL <= Lpaired) then
        x(:,:,2,iL) = 0.0_dp
      else
        if (mod(iL,2) == 1) then
          x(:,:,2,iL) = x(:,:,1,iL+1)
        else
          x(:,:,1,iL) = x(:,:,1,iL-1)
          x(:,:,2,iL) = -x(:,:,2,iL-1)
        end if
      end if
    end do

  end subroutine qmExpand3


  !> Decide a correct charge/magnetization in REKS
  subroutine qmExpand4(x, Lpaired)

    !> Array of data, fourth index spin, fifth index Lmax
    real(dp), intent(inout) :: x(:,:,:,:,:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    integer :: iL, Lmax

    Lmax = size(x,dim=5)

    do iL = 1, Lmax
      if (iL <= Lpaired) then
        x(:,:,:,2,iL) = 0.0_dp
      else
        if (mod(iL,2) == 1) then
          x(:,:,:,2,iL) = x(:,:,:,1,iL+1)
        else
          x(:,:,:,1,iL) = x(:,:,:,1,iL-1)
          x(:,:,:,2,iL) = -x(:,:,:,2,iL-1)
        end if
      end if
    end do

  end subroutine qmExpand4


  !> Decide a correct up/down in REKS
  subroutine udExpand3(x, Lpaired)

    !> array of data, third index spin, fourth index Lmax
    real(dp), intent(inout) :: x(:,:,:,:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    integer :: iL, Lmax

    Lmax = size(x,dim=4)

    do iL = 1, Lmax
      if (iL <= Lpaired) then
        x(:,:,2,iL) = x(:,:,1,iL)
      else
        if (mod(iL,2) == 1) then
          x(:,:,2,iL) = x(:,:,1,iL+1)
        else
          x(:,:,2,iL) = x(:,:,1,iL-1)
        end if
      end if
    end do

  end subroutine udExpand3


  !> Decide a correct up/down in REKS
  subroutine udExpand4(x, Lpaired)

    !> array of data, fourth index spin, fifth index Lmax
    real(dp), intent(inout) :: x(:,:,:,:,:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    integer :: iL, Lmax

    Lmax = size(x,dim=5)

    do iL = 1, Lmax
      if (iL <= Lpaired) then
        x(:,:,:,2,iL) = x(:,:,:,1,iL)
      else
        if (mod(iL,2) == 1) then
          x(:,:,:,2,iL) = x(:,:,:,1,iL+1)
        else
          x(:,:,:,2,iL) = x(:,:,:,1,iL-1)
        end if
      end if
    end do

  end subroutine udExpand4


  !> Convert the matrix from AO basis to MO basis
  subroutine matAO2MO(mat, eigenvecs)

    !> matrix converted from AO basis to MO basis
    real(dp), intent(inout) :: mat(:,:)

    !> eigenvectors
    real(dp), intent(in) :: eigenvecs(:,:)

    real(dp), allocatable :: tmpMat(:,:)
    integer :: nOrb

    nOrb = size(eigenvecs,dim=1)

    allocate(tmpMat(nOrb,nOrb))

    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, mat, eigenvecs)
    call gemm(mat, eigenvecs, tmpMat, transA='T')

  end subroutine matAO2MO


  !> Convert the matrix from MO basis to AO basis
  subroutine matMO2AO(mat, eigenvecs)

    !> matrix converted from MO basis to AO basis
    real(dp), intent(inout) :: mat(:,:)

    !> eigenvectors
    real(dp), intent(in) :: eigenvecs(:,:)

    real(dp), allocatable :: tmpMat(:,:)
    integer :: nOrb

    nOrb = size(eigenvecs,dim=1)

    allocate(tmpMat(nOrb,nOrb))

    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, mat, eigenvecs, transB='T')
    call gemm(mat, eigenvecs, tmpMat)

  end subroutine matMO2AO


  !> find proper string for active orbital
  subroutine getSpaceSym(ii, st)

    !> index for active space
    integer, intent(in) :: ii

    !> string for active space
    character(len=1), intent(inout) :: st

    if (ii == 1) then
      write(st,'(A1)') "a"
    else if (ii == 2) then
      write(st,'(A1)') "b"
    else if (ii == 3) then
      write(st,'(A1)') "c"
    else if (ii == 4) then
      write(st,'(A1)') "d"
    end if

  end subroutine getSpaceSym


  !> Find shell of index alpha with respect to mu (reference)
  subroutine findShellOfAO(al, mu, getAtomIndex, iSquare, iSpA, facP, facD)

    !> input AO index
    integer, intent(in) :: al

    !> reference AO index (standard of atom)
    integer, intent(in) :: mu

    !> get atom index from AO index
    integer, intent(in) :: getAtomIndex(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> output shell (s,p,d) for input AO index
    integer, intent(out) :: iSpA

    !> check whether al is included in p or d orbitals
    real(dp), intent(out) :: facP, facD

    integer :: iAtM, iAtA, iShA

    ! set the orbital shell for al index w.r.t mu
    iAtM = getAtomIndex(mu)
    iAtA = getAtomIndex(al)
    if (iAtA == iAtM) then
      iShA = al - iSquare(iAtA) + 1
      if (iShA == 1) then
        iSpA = 1
      else if (iShA > 1 .and. iShA <= 4) then
        iSpA = 2
        facP = 1.0_dp
      else if (iShA > 4 .and. iShA <= 9) then
        iSpA = 3
        facD = 1.0_dp
      end if
    else
      iSpA = 0
    end if

  end subroutine findShellOfAO


  !> Assign index in terms of dense form from super matrix form
  subroutine assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Number of vacant orbitals
    integer, intent(in) :: Nv

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> index for super matrix form
    integer, intent(in) :: ij

    !> index for dense form
    integer, intent(out) :: i, j

    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call assignIndex22_(Nc, Na, Nv, ij, i, j)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine assignIndex


  !> Assign converged epsilon value from fock matrix
  subroutine assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, t, &
      & chk, reksAlg, e1, e2)

    !> dense fock matrix for core orbitals
    real(dp), intent(in) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: Fa(:,:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> MO index for converged fock matrix
    integer, intent(in) :: i, j, t

    !> choice of calculations for converged fock matrix
    integer, intent(in) :: chk

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> output multiplier from converged fock matrix
    real(dp), intent(out) :: e1, e2

    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call assignEpsilon22_(Fc, Fa, SAweight, FONs, Nc, i, j, &
          & t, chk, e1, e2)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine assignEpsilon


  !> Assign average filling for i-th orbital
  subroutine assignFilling(FONs, SAweight, Nc, i, reksAlg, fi)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> orbital index
    integer, intent(in) :: i

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> output filling from fractional occupation numbers
    real(dp), intent(out) :: fi

    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call assignFilling22_(FONs, SAweight, Nc, i, fi)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine assignFilling


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Assign index in terms of dense form from super matrix form in REKS(2,2)
  subroutine assignIndex22_(Nc, Na, Nv, ij, i, j)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Number of vacant orbitals
    integer, intent(in) :: Nv

    !> index for super matrix form
    integer, intent(in) :: ij

    !> index for dense form
    integer, intent(out) :: i, j

    ! (i,j) = (core,active)
    if (ij <= Nc*Na) then
      i = mod(ij,Nc)
      j = ij / Nc + 1 + Nc
      if (i == 0) then
        i = i + Nc
        j = j - 1
      end if
    ! (i,j) = (core,vacant)
    else if (ij > Nc*Na .and. ij <= Nc*Na+Nc*Nv) then
      i = (ij-Nc*Na) / Nv + 1
      j = mod(ij-Nc*Na,Nv) + Nc + Na
      if (j == Nc+Na) then
        i = i - 1
        j = j + Nv
      end if
    ! (i,j) = (active,active)
    else if (ij > Nc*Na+Nc*Nv .and. ij <= Nc*Na+Nc*Nv+Na*(Na-1)*0.5_dp) then
      i = Nc + 1
      j = Nc + Na
    ! (i,j) = (active,vacant)
    else
      i = ( ij-Nc*Na-Nc*Nv-Na*(Na-1)*0.5_dp ) / Nv + 1 + Nc
      j = mod( int(ij-Nc*Na-Nc*Nv-Na*(Na-1)*0.5_dp) ,Nv) + Nc + Na
      if (j == Nc+Na) then
        i = i - 1
        j = j + Nv
      end if
    end if

  end subroutine assignIndex22_


!  !> Assign index in terms of super matrix form from dense form in REKS(2,2)
!  !> for only i < j case
!  subroutine assignIndexInverse22_(ij,i,j,Nc,Na,Nv)
!
!    integer, intent(out) :: ij
!    integer, intent(in) :: i, j
!    integer, intent(in) :: Nc, Na, Nv
!
!    if (i <= Nc) then
!      ! (i,j) = (core,active)
!      if (j > Nc .and. j <= Nc+Na) then
!        ij = i + (j-Nc-1) * Nc
!      end if
!      ! (i,j) = (core,vacant)
!      if (j > Nc+Na) then
!        ij = (i-1) * Nv + Nc*Na + j - Nc - Na
!      end if
!    end if
!    if (i > Nc .and. i <= Nc+Na) then
!      ! (i,j) = (active,active)
!      if (j > Nc .and. j <= Nc+Na) then
!        ij = (Na+Nv)*Nc + 1
!      end if
!      ! (i,j) = (active,vacant)
!      if (j > Nc+Na) then
!        ij = (i-1-Nc) * Nv + Nc*Na + Nc*Nv + Na*(Na-1)/DBLE(2) + j - Nc - Na
!      end if
!    end if
!
!  end subroutine assignIndexInverse22_


  !> Assign converged epsilon value from fock matrix in REKS(2,2)
  subroutine assignEpsilon22_(Fc, Fa, SAweight, FONs, Nc, i, j, &
      & t, chk, e1, e2)

    !> dense fock matrix for core orbitals
    real(dp), intent(in) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: Fa(:,:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> MO index for converged fock matrix
    integer, intent(in) :: i, j, t

    !> choice of calculations for converged fock matrix
    integer, intent(in) :: chk

    !> output multiplier from converged fock matrix
    real(dp), intent(out) :: e1, e2

    real(dp) :: n_a, n_b

    n_a = FONs(1,1)
    n_b = FONs(2,1)

    if (chk == 1) then

      if (i <= Nc) then
        e1 = Fc(j,t) * 2.0_dp
      else if (i == Nc + 1) then
        e1 = Fa(j,t,1) * (SAweight(1)*n_a+SAweight(2))
      else if (i == Nc + 2) then
        e1 = Fa(j,t,2) * (SAweight(1)*n_b+SAweight(2))
      else
        e1 = Fc(j,t) * 0.0_dp
      end if

      if (j <= Nc) then
        e2 = Fc(j,t) * 2.0_dp
      else if (j == Nc + 1) then
        e2 = Fa(j,t,1) * (SAweight(1)*n_a+SAweight(2))
      else if (j == Nc + 2) then
        e2 = Fa(j,t,2) * (SAweight(1)*n_b+SAweight(2))
      else
        e2 = Fc(j,t) * 0.0_dp
      end if

    else if (chk == 2) then

      if (i <= Nc) then
        e1 = Fc(i,t) * 2.0_dp
      else if (i == Nc + 1) then
        e1 = Fa(i,t,1) * (SAweight(1)*n_a+SAweight(2))
      else if (i == Nc + 2) then
        e1 = Fa(i,t,2) * (SAweight(1)*n_b+SAweight(2))
      else
        e1 = Fc(i,t) * 0.0_dp
      end if

      if (j <= Nc) then
        e2 = Fc(i,t) * 2.0_dp
      else if (j == Nc + 1) then
        e2 = Fa(i,t,1) * (SAweight(1)*n_a+SAweight(2))
      else if (j == Nc + 2) then
        e2 = Fa(i,t,2) * (SAweight(1)*n_b+SAweight(2))
      else
        e2 = Fc(i,t) * 0.0_dp
      end if

    else if (chk == 3) then

      if (i <= Nc) then
        e1 = Fc(t,i) * 2.0_dp
      else if (i == Nc + 1) then
        e1 = Fa(t,i,1) * (SAweight(1)*n_a+SAweight(2))
      else if (i == Nc + 2) then
        e1 = Fa(t,i,2) * (SAweight(1)*n_b+SAweight(2))
      else
        e1 = Fc(t,i) * 0.0_dp
      end if

      if (j <= Nc) then
        e2 = Fc(t,i) * 2.0_dp
      else if (j == Nc + 1) then
        e2 = Fa(t,i,1) * (SAweight(1)*n_a+SAweight(2))
      else if (j == Nc + 2) then
        e2 = Fa(t,i,2) * (SAweight(1)*n_b+SAweight(2))
      else
        e2 = Fc(t,i) * 0.0_dp
      end if

    end if

  end subroutine assignEpsilon22_


  !> Assign average filling for i-th orbital in REKS(2,2)
  subroutine assignFilling22_(FONs, SAweight, Nc, i, fi)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> orbital index
    integer, intent(in) :: i

    !> output filling from fractional occupation numbers
    real(dp), intent(out) :: fi

    real(dp) :: n_a, n_b

    n_a = FONs(1,1)
    n_b = FONs(2,1)

    if (i <= Nc) then
      fi = 1.0_dp
    else if (i == Nc + 1) then
      fi = (SAweight(1)*n_a + SAweight(2)) * 0.5_dp
    else if (i == Nc + 2) then
      fi = (SAweight(1)*n_b + SAweight(2)) * 0.5_dp
    else
      fi = 0.0_dp
    end if

  end subroutine assignFilling22_


end module dftbp_rekscommon
