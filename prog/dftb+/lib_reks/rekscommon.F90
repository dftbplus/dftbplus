!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

! TODO
!!!!#:include 'common.fypp'

!> REKS and SI-SA-REKS formulation in DFTB as developed by Lee et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet appart from Gamma point.
!> * Orbital potentials or spin-orbit or external E-field does not work yet.
!> * Only for closed shell system.
!> * Onsite corrections are not included in this version
! TODO
!> * Dispersion would be combined with REKS
module dftbp_rekscommon

  use dftbp_accuracy
!  use dftbp_assert
  use dftbp_blasroutines, only : gemm
  use dftbp_densedescr
  use dftbp_message
  use dftbp_reksvar

  implicit none

  private

  public :: checkGammaPoint
  public :: qm2udL, ud2qmL
  public :: qmExpandL
  public :: matAO2MO, matMO2AO
  public :: getSpaceSym

  !> Swap from charge/magnetisation to up/down in REKS
  interface qm2udL
    module procedure qm2ud2
    module procedure qm2ud3
  end interface qm2udL

  !> Swap from up/down to charge/magnetisation in REKS
  interface ud2qmL
    module procedure ud2qm2
    module procedure ud2qm3
  end interface ud2qmL

  !> Set correct charge/magnetization in REKS
  interface qmExpandL
    module procedure qmExpand3
    module procedure qmExpand4
  end interface qmExpandL

  contains

  !> Check whether the cell size is proper to the Gamma point
  !> calculation or not, and set several convenient variables
  subroutine checkGammaPoint(denseDesc, iNeighbour, nNeighbourSK,&
      & iPair, img2CentCell, over, reks)

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
    type(TReksCalc), intent(inout) :: reks

    integer :: mu, nu, nAtom, nOrb, nAtomSparse
    integer :: iAtom1, iAtom2, iAtom2f, iNeigh1, iOrig1
    integer :: nOrb1, nOrb2, ii, jj, kk, ll

    nAtom = size(denseDesc%iAtomStart,dim=1) - 1
    nOrb = size(reks%over,dim=1)

    nAtomSparse = 0
    do iAtom1 = 1, nAtom ! mu
      nAtomSparse = nAtomSparse + nNeighbourSK(iAtom1) + 1
    end do

    deallocate(reks%getDenseAtom)
    allocate(reks%getDenseAtom(nAtomSparse,2))

    ll = 1
    reks%getDenseAO(:,:) = 0
    reks%getDenseAtom(:,:) = 0
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
          if (reks%over(mu,nu) /= over(iOrig1+kk-1)) then
            call error("Inconsistent maching exists between sparse and dense.")
          end if
          reks%getDenseAO(iOrig1+kk-1,1) = mu
          reks%getDenseAO(iOrig1+kk-1,2) = nu
        end do

        reks%getDenseAtom(ll,1) = iAtom1  ! A atom
        reks%getDenseAtom(ll,2) = iAtom2f ! B atom
        ll = ll + 1
      end do
    end do

    do mu = 1, nOrb
      do iAtom1 = 1, nAtom
        if (mu > denseDesc%iAtomStart(iAtom1)-1 .and.&
            & mu <= denseDesc%iAtomStart(iAtom1+1)-1) then
          reks%getAtomIndex(mu) = iAtom1
        end if
      end do
    end do

  end subroutine checkGammaPoint


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

    ! ... use external blas library for matrix multiplication
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat,mat,eigenvecs,1.0_dp,0.0_dp,'N','N')
    call gemm(mat,eigenvecs,tmpMat,1.0_dp,0.0_dp,'T','N')

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
    call gemm(tmpMat,mat,eigenvecs,1.0_dp,0.0_dp,'N','T')
    call gemm(mat,eigenvecs,tmpMat,1.0_dp,0.0_dp,'N','N')

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


end module dftbp_rekscommon
