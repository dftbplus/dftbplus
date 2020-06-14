!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Construct Pipek-Mezey localised orbitals, either for molecules/gamma point periodic, or for each
!> k-point separately. Note that for the k-point case these are NOT localised Wannier functions as
!> each k-point is localised independently.
module dftbp_pmlocalisation
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_globalenv, only : stdOut
  use dftbp_blasroutines
  use dftbp_sparse2dense, only :unpackHS
  use dftbp_sorting
  use dftbp_message
  use dftbp_periodic, only : TNeighbourList
  implicit none
  private

  public :: TPipekMezeyInp, TPipekMezey, initialise


  !> Input data for Pipek-Mezey calculator
  type :: TPipekMezeyInp

    !> halting tolerance for localisation
    real(dp) :: tolerance

    !> number of localisation iterations
    integer :: maxIter

    !> Optional tolerances for element neglect if instead using a sparse version of Pipek-Mezey
    !> localisation
    real(dp), allocatable :: sparseTols(:)

  end type TPipekMezeyInp


  !> Pipek-Mezey localisation calculator
  type :: TPipekMezey
    private

    !> tolerances for element neglect if instead using a sparse version of Pipek-Mezey localisation
    real(dp), allocatable :: sparseTols(:)

    !> halting tolerance for localisation
    real(dp) :: tolerance

    !> number of localisation iterations
    integer :: maxIter

  contains
    private

    !> Performs Pipek-Mezey localisation for real case
    procedure :: calcCoeffsReal

    !> Performs Pipek-Mezey localisation for complex
    procedure :: calcCoeffsKPoint

    !> evaluate the localisation measure for real cases
    procedure, nopass :: getLocalisationReal

    !> evaluate the localisation measure for complex cases
    procedure, nopass :: getLocalisationKPoint

    !> Performs localisation on orbitals
    generic, public :: calcCoeffs => calcCoeffsReal, calcCoeffsKPoint

    !> Value of the localisation measure for orbitals
    generic, public :: getLocalisation => getLocalisationReal, getLocalisationKPoint

  end type TPipekMezey


  !> Initializes calculator instance.
  interface initialise
    module procedure TPipekMezey_initialise
  end interface initialise


contains

  !> Initialises calculator instance.
  subroutine TPipekMezey_initialise(this, input)

    !> Instance.
    type(TPipekMezey), intent(out) :: this

    !> Input data.
    type(TPipekMezeyInp), intent(inout) :: input

    this%tolerance = input%tolerance
    this%maxIter = input%maxIter
    call move_alloc(input%sparseTols, this%sparseTols)
    if (allocated(this%sparseTols)) then
      if (any(this%sparseTols < epsilon(0.0_dp))) then
        call error('Tolerances for sparse Pipek-Mezey localisation too small.')
      end if
    end if

  end subroutine TPipekMezey_initialise


  !> Performs Pipek-Mezey localisation for a molecule.
  subroutine calcCoeffsReal(this, ci, SSqrReal, iAtomStart)

    !> Instance
    class(TPipekMezey), intent(in) :: this

    !> wavefunction coefficients
    real(dp), intent(inout) :: ci(:,:)

    !> overlap matrix in square form
    real(dp), intent(in) :: SSqrReal(:,:)

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    integer :: ii

    if (allocated(this%sparseTols)) then
      do ii = 1, size(this%sparseTols)
        call PipekMezeySuprtRegion_real(ci, SSqrReal, iAtomStart, this%tolerance, this%maxIter,&
            & this%sparseTols(ii))
      end do
    else
      call pipekMezeyOld_real(ci, SSqrReal, iAtomStart, this%tolerance, this%maxIter)
    end if

  end subroutine calcCoeffsReal


  !> Performs Pipek-Mezey localisation for a periodic system at a specified k-point.
  subroutine calcCoeffsKPoint(this, ci, SSqrCplx, over, kPoint, neighbourList, nNeighbourSK,&
      & iCellVec, cellVec, iAtomStart, iPair, img2CentCell)

    !> Instance.
    class(TPipekMezey), intent(in) :: this

    !> wavefunction coefficients
    complex(dp), intent(inout) :: ci(:,:)

    !> overlap matrix, used as workspace
    complex(dp), intent(inout) :: SSqrCplx(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> current k-point
    real(dp), intent(in) :: kPoint(:)

    !> neighbour list
    type(TNeighbourList), intent(in) :: neighbourList

    !> number of neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> list of which image cells atoms outside the central cell fall into
    integer, intent(in) :: iCellVec(:)

    !> vectors to the image cells
    real(dp), intent(in) :: cellVec(:,:)

    !> index for the square matrices
    integer, intent(in) :: iAtomStart(:)

    !> index for the sparse matrices
    integer, intent(in) :: iPair(0:,:)

    !> index array back to central cell
    integer, intent(in) :: img2CentCell(:)

    call PipekMezeyOld_kpoint(ci, SSqrCplx, over, kPoint, neighbourList%iNeighbour, nNeighbourSK,&
        & iCellVec, cellVec, iAtomStart, iPair, img2CentCell, this%tolerance, this%maxIter)

  end subroutine calcCoeffsKPoint


  !> Localisation value of square of Mulliken charges summed over all levels
  function getLocalisationReal(ci, SSqrReal, iAtomStart) result(locality)

    !> wavefunction coefficients
    real(dp), intent(in) :: ci(:,:)

    !> overlap matrix
    real(dp), intent(in) :: SSqrReal(:,:)

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Calculated locality
    real(dp) :: locality

    locality = PipekMezyLocality_real(ci, SSqrReal, iAtomStart)

  end function getLocalisationReal



  !> Localisation value of square of Mulliken charges summed over all levels for each k-point.
  function getLocalisationKPoint(ci, SSqrCplx, over, kpoint, neighbourList, nNeighbourSK,&
      & iCellVec, cellVec, iAtomStart, iPair, img2CentCell)  result (locality)

    !> wavefunction coefficients
    complex(dp), intent(in) :: ci(:,:)

    !> overlap matrix, used as workspace
    complex(dp), intent(inout) :: SSqrCplx(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> current k-point
    real(dp), intent(in) :: kpoint(:)

    !> neighbour list
    type(TNeighbourList), intent(in) :: neighbourList

    !> number of neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> list of which image cells atoms outside the central cell fall into
    integer, intent(in) :: iCellVec(:)

    !> vectors to the image cells
    real(dp), intent(in) :: cellVec(:,:)

    !> index for the square matrices
    integer, intent(in) :: iAtomStart(:)

    !> index for the sparse matrices
    integer, intent(in) :: iPair(0:,:)

    !> index array back to central cell
    integer, intent(in) :: img2CentCell(:)

    !> Locality for current k-point
    real(dp) :: locality

    locality = PipekMezyLocality_kpoint(ci, SSqrCplx, over, kpoint, neighbourList%iNeighbour,&
        & nNeighbourSK, iCellVec, cellVec, iAtomStart, iPair, img2CentCell)

  end function getLocalisationKPoint




! Private functions below



  !> Performs conventional Pipek-Mezey localisation for a molecule given the square overlap matrix
  !> using iterative sweeps over each pair of orbitals
  subroutine PipekMezeyOld_real(ci, S, iAtomStart, pipekTol, mIter)

    !> wavefunction coefficients
    real(dp), intent(inout) :: ci(:,:)

    !> overlap matrix in square form
    real(dp), intent(in) :: S(:,:)

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> tolerance for halting localisation iterations
    real(dp), intent(in) :: pipekTol

    !> maximum number of iterations to use
    integer, intent(in), optional :: mIter

    integer :: iLev1, iLev2, nLev
    integer :: iAtom1, nAtom, nIter
    integer :: iOrb1, iOrb2, nOrb
    integer :: iIter
    real(dp) :: Ast, Bst, C4A, AB
    real(dp) :: sina, cosa, Pst, Pss, Ptt
    real(dp), allocatable :: Sci1(:), Sci2(:,:), ciTmp1(:), ciTmp2(:)

    real(dp) :: alpha, alphaMax, conv
    real(dp) :: alphalast = 1.0_dp
    logical :: tConverged = .false.

    @:ASSERT(size(ci,dim=1)>=size(ci,dim=2))
    @:ASSERT(size(ci,dim=1)==size(S,dim=1))
    @:ASSERT(size(S,dim=1)==size(S,dim=2))

    nOrb = size(ci,dim=1)
    nLev = size(ci,dim=2)
    nAtom = size(iAtomStart) -1

    @:ASSERT(iAtomStart(nAtom+1)-1 == nOrb)

    if (present(mIter)) then
      nIter = mIter
    else
      nIter = 20
    end if

    allocate(Sci1(nOrb))
    allocate(Sci2(nOrb,2:nLev))

    allocate(ciTmp1(nOrb))
    allocate(ciTmp2(nOrb))

    lpLocalise: do iIter = 1, nIter
      alphamax = 0.0_dp
      ! Sweep over all pairs of levels
      write(stdout, *)'Iter', iIter
      do iLev1 = 1, nLev

        if (iLev1 < nLev) then
          call symm(Sci2(1:nOrb,iLev1+1:nLev),'l',S,ci(:,iLev1+1:nLev),'L')
        else
          call hemv(Sci2(1:nOrb,nLev),S,ci(:,nLev))
        end if

        do iLev2 = iLev1 +1, nLev

          call hemv(Sci1,S,ci(:,iLev1))

          Ast = 0.0_dp
          Bst = 0.0_dp
          do iAtom1 = 1, nAtom
            iOrb1 = iAtomStart(iAtom1)
            iOrb2 = iAtomStart(iAtom1+1) - 1
            Pst = 0.5_dp * (sum(ci(iOrb1:iOrb2,iLev1)*Sci2(iOrb1:iOrb2,iLev2))&
                & + sum(ci(iOrb1:iOrb2,iLev2)*Sci1(iOrb1:iOrb2)))
            Pss = sum ( ci(iOrb1:iOrb2,iLev1)*Sci1(iOrb1:iOrb2) )
            Ptt = sum ( ci(iOrb1:iOrb2,iLev2)*Sci2(iOrb1:iOrb2,iLev2) )
            Ast = Ast + Pst * Pst - 0.25_dp * (Pss - Ptt) * (Pss - Ptt)
            Bst = Bst + Pst * (Pss - Ptt)
          end do

          AB = Ast * Ast + Bst * Bst
          if (abs(AB)>0.0_dp) then
            C4A = -Ast / sqrt(AB)
            alpha = 0.25_dp * acos(C4A)
            if (Bst <= 0.0) then
              alpha = -alpha
            end if
          else
            alpha = 0.0_dp
          end if

          if (alphamax < abs(alpha)) then
            alphamax = alpha
          end if

          ! now we have to mix the two orbitals
          SINA=SIN(alpha)
          COSA=COS(alpha)
          ciTmp1 = ci(:,iLev1)
          ciTmp2 = ci(:,iLev2)
          ci(:,iLev1) = COSA*ciTmp1 + SINA*ciTmp2
          ci(:,iLev2) = COSA*ciTmp2 - SINA*ciTmp1

        end do
      end do

      conv = abs(alphamax) - abs(alphalast)
      if (iIter > 2 .and. ((abs(conv)<pipekTol) .or. alphamax == 0.0)) then
        tConverged = .true.
        exit
      end if
      alphalast = alphamax

    end do lpLocalise

    if (.not.tConverged) then
      call warning("Exceeded iterations in Pipek-Mezey localisation!")
    end if

  end subroutine PipekMezeyOld_real


  !> performs Pipek-Mezey localisation for a molecule given the square overlap matrix, using a
  subroutine PipekMezeySuprtRegion_real(ci, S, iAtomStart, convergence, mIter, RegionTol)

    !> support region for each molecular orbital
    real(dp), intent(inout) :: ci(:,:)

    !> wavefunction coefficients
    real(dp), intent(in) :: S(:,:)

    !> overlap matrix in square form
    integer, intent(in) :: iAtomStart(:)

    !> Atom offset for the squared Hamiltonian
    real(dp), intent(in) :: convergence

    !> tolerance for halting localisation iterations
    integer, intent(in), optional :: mIter

    !> maximum number of iterations to use
    real(dp), intent(in) :: RegionTol

    !> Number of electrons per site to consider it within a domain.
    integer :: iLev1, iLev2, nLev
    integer :: iAtom1, nAtom, nIter
    integer :: iOrb1, iOrb2, nOrb
    integer :: iIter
    real(dp) :: Ast, Bst, C4A, AB
    real(dp) :: sina, cosa, Pst, Pss, Ptt
    real(dp), allocatable :: Sci1(:,:), Sci2(:,:), ciTmp1(:), ciTmp2(:)
    integer, allocatable :: LevAtAtom(:,:), nLevAtAtom(:)
    integer, allocatable :: SitesLev(:,:), nSitesLev(:)
    integer, allocatable :: LevPairs(:)
    integer :: ii, jj, kk, ll, ij, iLev, nLevPairs
    logical :: tPair, tPresent

    integer, allocatable :: oldSites(:,:)
    integer :: nOldSites(2)

    real(dp) :: alpha, alphaMax, conv
    real(dp) :: alphalast = 1.0_dp
    real(dp) :: rCount
    logical :: tConverged = .false.

    real(dp) :: Localisation, oldLocalisation
    integer, allocatable :: union(:)

    write(stdout, *)'Pipek Mezey localisation'

    Localisation = PipekMezyLocality_real(ci,S,iAtomStart)
    write(stdout, *)'Initial', Localisation

    @:ASSERT(size(ci,dim=1)>=size(ci,dim=2))
    @:ASSERT(size(ci,dim=1)==size(S,dim=1))
    @:ASSERT(size(S,dim=1)==size(S,dim=2))

    nOrb = size(ci,dim=1)
    nLev = size(ci,dim=2)
    nAtom = size(iAtomStart) -1

    @:ASSERT(iAtomStart(nAtom+1)-1 == nOrb)

    if (present(mIter)) then
      nIter = mIter
    else
      nIter = 20
    end if

    allocate(Sci1(nOrb,2))
    allocate(Sci2(nOrb,nLev))

    allocate(ciTmp1(nOrb))
    allocate(ciTmp2(nOrb))

    allocate(oldSites(nAtom,2))

    allocate(LevAtAtom(nLev,nAtom))
    allocate(nLevAtAtom(nAtom))
    LevAtAtom = 0
    nLevAtAtom = 0
    allocate(SitesLev(nAtom,nLev))
    allocate(nSitesLev(nLev))
    SitesLev = 0
    nSitesLev = 0
    allocate(LevPairs(nLev))
    LevPairs = 0

    allocate(union(2*nAtom))

    ! make Mulliken charges for each level
    call symm(Sci2,'L',S,ci(:,1:nLev),'L')
    Sci2 = Sci2 * ci(:,1:nLev)

    nSitesLev(:) = 0
    SitesLev(:,:) = 0
    do iLev1 = 1, nLev
      do iAtom1 = 1, nAtom
        iOrb1 = iAtomStart(iAtom1)
        iOrb2 = iAtomStart(iAtom1+1) - 1
        if (sum(abs(Sci2(iOrb1:iOrb2,iLev1)))>=RegionTol) then
          nSitesLev(iLev1) = nSitesLev(iLev1) + 1
          SitesLev(nSitesLev(iLev1),iLev1) = iAtom1
        end if
      end do
    end do

    nLevAtAtom(:) = 0
    LevAtAtom(:,:) = 0
    do iLev1 = 1, nLev
      do jj = 1, nSitesLev(iLev1)
        nLevAtAtom(SitesLev(jj,iLev1)) = nLevAtAtom(SitesLev(jj,iLev1)) + 1
        LevAtAtom(nLevAtAtom(SitesLev(jj,iLev1)),SitesLev(jj,iLev1)) = iLev1
      end do
    end do

    lpLocalise: do iIter = 1, nIter
      alphamax = 0.0_dp

      write(stdout, "(' Iter:',I0,', tol:',E10.2)")iIter,RegionTol
      rCount = 0.0

      do iLev1 = 1, nLev

        if (real(iLev1)/real(nLev) > rCount) then
          write(stdout, "(1X,I0,'%')")int(100*real(iLev1)/real(nLev))
          rCount = rCount + 0.1 ! every 10%
        end if

        nLevPairs = 0
        LevPairs(:) = 0
        do ii = 1, nSitesLev(iLev1)
          iAtom1 = SitesLev(ii,iLev1)
          do jj = 1, nLevAtAtom(iAtom1)
            if (LevAtAtom(jj,iAtom1) > iLev1) then
              tPair = .true.
              do kk = 1, nLevPairs
                if (LevPairs(kk) == LevAtAtom(jj,iAtom1)) then  !already have it
                  tPair = .false.
                  exit
                end if
              end do
              if (tPair) then
                nLevPairs = nLevPairs + 1
                LevPairs(nLevPairs) = LevAtAtom(jj,iAtom1)
              end if
            end if
          end do
        end do

        ! Sweep over all pairs of levels with shared regions of charge
        do ii = 1, nLevPairs
          iLev2 = LevPairs(ii)

          call hemv(Sci1(:,1),S,ci(:,iLev1))
          call hemv(Sci1(:,2),S,ci(:,iLev2))

          Ast = 0.0_dp
          Bst = 0.0_dp

          ! Find atomic sites that appear in both localised orbitals
          union = 0
          union(:nSitesLev(iLev1)) = SitesLev(:nSitesLev(iLev1),iLev1)
          union(nSitesLev(iLev1)+1:nSitesLev(iLev1)+nSitesLev(iLev2)) =&
              & SitesLev(:nSitesLev(iLev2),iLev2)
          call heap_sort(union(:nSitesLev(iLev1)+nSitesLev(iLev2)))
          kk = unique(union,nSitesLev(iLev1)+nSitesLev(iLev2))

          do ll = 1, kk
            iAtom1 = union(ll)
            !  regions
            iOrb1 = iAtomStart(iAtom1)
            iOrb2 = iAtomStart(iAtom1+1) - 1

            Pst = 0.5_dp * (sum(ci(iOrb1:iOrb2,iLev1)*Sci1(iOrb1:iOrb2,2))&
                & + sum(ci(iOrb1:iOrb2,iLev2)*Sci1(iOrb1:iOrb2,1)))
            Pss = sum ( ci(iOrb1:iOrb2,iLev1)*Sci1(iOrb1:iOrb2,1) )
            Ptt = sum ( ci(iOrb1:iOrb2,iLev2)*Sci1(iOrb1:iOrb2,2) )
            Ast = Ast + Pst * Pst - 0.25_dp * (Pss - Ptt) * (Pss - Ptt)
            Bst = Bst + Pst * (Pss - Ptt)
          end do

          AB = Ast * Ast + Bst * Bst
          if (abs(AB)>0.0_dp) then
            C4A = -Ast / sqrt(AB)
            alpha = 0.25_dp * acos(C4A)
            if (Bst <= 0.0) then
              alpha = -alpha
            end if
          else
            alpha = 0.0_dp
          end if

          if (alphamax < abs(alpha)) then
            alphamax = alpha
          end if

          ! now we have to mix the two orbitals
          SINA=SIN(alpha)
          COSA=COS(alpha)
          ciTmp1 = ci(:,iLev1)
          ciTmp2 = ci(:,iLev2)
          ci(:,iLev1) = COSA*ciTmp1 + SINA*ciTmp2
          ci(:,iLev2) = COSA*ciTmp2 - SINA*ciTmp1

          nOldSites(1) = nSitesLev(iLev1)
          oldSites(1:nOldSites(1),1) = SitesLev(1:nOldSites(1),iLev1)
          nOldSites(2) = nSitesLev(iLev2)
          oldSites(1:nOldSites(2),2) = SitesLev(1:nOldSites(2),iLev2)

          ! need to update index information now
          jj = iLev1
          call hemv(Sci2(:,jj),S,ci(:,jj))
          Sci2(:,jj) = Sci2(:,jj) * ci(:,jj)
          nSitesLev(jj) = 0
          SitesLev(:,jj) = 0
          do kk = 1, nOldSites(1)
            iAtom1 = oldSites(kk,1)
            iOrb1 = iAtomStart(iAtom1)
            iOrb2 = iAtomStart(iAtom1+1) - 1
            if (sum(abs(Sci2(iOrb1:iOrb2,jj)))>=RegionTol) then
              nSitesLev(jj) = nSitesLev(jj) + 1
              SitesLev(nSitesLev(jj),jj) = iAtom1
            end if
          end do
          jj = iLev2
          call hemv(Sci2(:,jj),S,ci(:,jj))
          Sci2(:,jj) = Sci2(:,jj) * ci(:,jj)
          nSitesLev(jj) = 0
          SitesLev(:,jj) = 0
          do kk = 1, nOldSites(2)
            iAtom1 = oldSites(kk,2)
            iOrb1 = iAtomStart(iAtom1)
            iOrb2 = iAtomStart(iAtom1+1) - 1
            if (sum(abs(Sci2(iOrb1:iOrb2,jj)))>=RegionTol) then
              nSitesLev(jj) = nSitesLev(jj) + 1
              SitesLev(nSitesLev(jj),jj) = iAtom1
            end if
          end do

          do jj = 1, nOldSites(1)
            iAtom1 = oldSites(jj,1) ! was a site of level1
            do kk = 1, nLevAtAtom(iAtom1)
              iLev = LevAtAtom(kk,iAtom1)
              if (iLev == iLev1) then ! this was a level in the atom
                tPresent = .false. ! is it still present at that site?
                do ij = 1, nSitesLev(jj)
                  if (iAtom1 == SitesLev(ij,iLev1)) then
                    tPresent = .true.
                    exit
                  end if
                end do
                if (.not.tPresent) then
                  LevAtAtom(kk:nLevAtAtom(iAtom1)-1,iAtom1) =&
                      & LevAtAtom(kk+1:nLevAtAtom(iAtom1),iAtom1)
                  nLevAtAtom(iAtom1) = nLevAtAtom(iAtom1) -1
                end if
              end if
            end do
          end do

          do jj = 1, nOldSites(2)
            iAtom1 = oldSites(jj,2) ! was a site of level2
            do kk = 1, nLevAtAtom(iAtom1)
              iLev = LevAtAtom(kk,iAtom1)
              if (iLev == iLev2) then ! this was a level in the atom
                tPresent = .false. ! is it still present at that site?
                do ij = 1, nSitesLev(jj)
                  if (iAtom1 == SitesLev(ij,iLev2)) then
                    tPresent = .true.
                    exit
                  end if
                end do
                if (.not.tPresent) then
                  LevAtAtom(kk:nLevAtAtom(iAtom1)-1,iAtom1) =&
                      & LevAtAtom(kk+1:nLevAtAtom(iAtom1),iAtom1)
                  nLevAtAtom(iAtom1) = nLevAtAtom(iAtom1) -1
                end if
              end if
            end do
          end do

        end do
      end do

      oldLocalisation = Localisation
      Localisation = PipekMezyLocality_real(ci,S,iAtomStart)
      write(stdout, "(A,F12.6,1X,A,E20.12)")'Current localisation ',Localisation,&
          & 'change ',Localisation-oldLocalisation

      conv = abs(alphamax) - abs(alphalast)
      if (iIter > 2 .and. ((abs(conv)<convergence) .or. alphamax == 0.0)) then
        write(stdout, *)'Converged on rotation angle'
        tConverged = .true.
        exit
      end if

      conv = abs(Localisation-oldLocalisation)
      if (abs(conv)<convergence) then
        write(stdout, *)'Converged on localization value.'
        tConverged = .true.
        exit
      end if

      alphalast = alphamax
      write(stdout, "(' max(alpha)',E10.2)")alphamax

    end do lpLocalise

    Localisation = PipekMezyLocality_real(ci,S,iAtomStart)
    write(stdout, *)'Final',Localisation

    if (.not.tConverged) then
      write(stdout, *)alphamax
      call warning("Exceeded iterations in Pipek-Mezey localisation!")
    end if

  end subroutine PipekMezeySuprtRegion_real


  !> Localisation value of square of Mulliken charges summed over all levels
  function PipekMezyLocality_real(ci,S,iAtomStart) result(PipekMezyLocality)

    !> Localisation
    real(dp) :: PipekMezyLocality

    !> wavefunction coefficients
    real(dp), intent(in) :: ci(:,:)

    !> overlap matrix
    real(dp), intent(in) :: S(:,:)

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    real(dp), allocatable :: Sci(:,:)
    integer :: nAtom, iAtom, iOrbStart, iOrbEnd, nOrb, nLev

    nAtom = size(iAtomStart) -1
    nOrb = size(ci,dim=1)
    nLev = size(ci,dim=2)

    allocate(Sci(nOrb,nLev))

    PipekMezyLocality = 0.0_dp

    call symm(Sci,'L',S,ci,'L')

    Sci = ci * Sci
    do iAtom = 1, nAtom
      iOrbStart = iAtomStart(iAtom)
      iOrbEnd = iAtomStart(iAtom+1) - 1
      PipekMezyLocality = PipekMezyLocality + sum(sum(Sci(iOrbStart:iOrbEnd,1:nLev),dim=1)**2)
    end do

  end function PipekMezyLocality_real


  !> Localisation value of square of Mulliken charges at a k-point
  function PipekMezyLocality_kpoint(ci, S, over, kpoint, iNeighbour, nNeighbourSK, iCellVec,&
      & cellVec, iAtomStart, iPair, img2CentCell)  result (PipekMezyLocality)

    !> wavefunction coefficients
    complex(dp), intent(in) :: ci(:,:)

    !> overlap matrix, used as workspace
    complex(dp), intent(inout) :: S(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> current k-point
    real(dp), intent(in) :: kpoint(:)

    !> neighbour list
    integer, intent(in) :: iNeighbour(0:,:)

    !> number of neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> list of which image cells atoms outside the central cell fall into
    integer, intent(in) :: iCellVec(:)

    !> vectors to the image cells
    real(dp), intent(in) :: cellVec(:,:)

    !> index for the square matrices
    integer, intent(in) :: iAtomStart(:)

    !> index for the sparse matrices
    integer, intent(in) :: iPair(0:,:)

    !> index array back to central cell
    integer, intent(in) :: img2CentCell(:)

    !> Locality for each k-point
    real(dp) :: PipekMezyLocality

    complex(dp), allocatable :: Sci(:,:)
    real(dp), allocatable :: tmp(:,:)
    integer :: nAtom, iAtom, iOrbStart, iOrbEnd, nOrb, iLev, nLev

    @:ASSERT(size(ci,dim=1)>=size(ci,dim=2))

    nAtom = size(iAtomStart) -1
    nOrb = size(ci,dim=1)
    nLev = size(ci,dim=2)

    @:ASSERT(all(shape(kpoint) == [3]))
    @:ASSERT(all(shape(S) == [nOrb,nOrb]))

    allocate(Sci(nOrb,nLev))
    allocate(tmp(nAtom,nLev))

    PipekMezyLocality = 0.0_dp

    tmp = 0.0_dp

    call unpackHS(S, over, kPoint, iNeighbour, nNeighbourSK, iCellVec, cellVec, iAtomStart,&
        & iPair, img2CentCell)

    call hemm(Sci,'L',S,ci,'L')
    Sci = conjg(ci) * Sci

    do iLev = 1, nLev
      do iAtom = 1, nAtom
        iOrbStart = iAtomStart(iAtom)
        iOrbEnd = iAtomStart(iAtom+1) - 1
        tmp(iAtom, iLev) = tmp(iAtom, iLev) + sum(real(Sci(iOrbStart:iOrbEnd,iLev)))
      end do
    end do
    PipekMezyLocality = sum(tmp**2)

  end function PipekMezyLocality_kpoint


  !> Performs conventional Pipek-Mezey localisation for a supercell using iterative sweeps over each
  !> pair of orbitals for a particular k and spin sub-matrix
  subroutine PipekMezeyOld_kpoint(ci, S, over, kpoint, iNeighbour, nNeighbourSK, iCellVec, cellVec,&
      & iAtomStart, iPair, img2CentCell, convergence, mIter)

    !> wavefunction coefficients
    complex(dp), intent(inout) :: ci(:,:)

    !> overlap matrix, used as workspace
    complex(dp), intent(inout) :: S(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> current k-point
    real(dp), intent(in) :: kpoint(:)

    !> neighbour list
    integer, intent(in) :: iNeighbour(0:,:)

    !> number of neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> list of which image cells atoms outside the central cell fall into
    integer, intent(in) :: iCellVec(:)

    !> vectors to the image cells
    real(dp), intent(in) :: cellVec(:,:)

    !> index for the square matrices
    integer, intent(in) :: iAtomStart(:)

    !> index for the sparse matrices
    integer, intent(in) :: iPair(0:,:)

    !> index array back to central cell
    integer, intent(in) :: img2CentCell(:)

    !> tolerance for halting localisation iterations
    real(dp), intent(in) :: convergence

    !> maximum number of iterations to use
    integer, intent(in), optional :: mIter

    integer :: iLev1, iLev2, nLev
    integer :: iAtom1, nAtom, nIter
    integer :: iOrb1, iOrb2, nOrb
    integer :: iIter, iLoc(1), ii
    real(dp) :: Ast, Bst, C4A, AB
    real(dp) :: sina, cosa
    complex(dp) :: Pst, Pss, Ptt
    complex(dp), allocatable :: Sci1(:), Sci2(:,:)
    complex(dp), allocatable :: ciTmp1(:), ciTmp2(:)
    complex(dp) :: phase
    complex(dp), parameter :: im = (0.0_dp,1.0_dp)

    real(dp) :: alpha, alphaMax, conv
    real(dp) :: alphalast = 1.0_dp
    logical :: tConverged

    @:ASSERT(size(ci,dim=1)>=size(ci,dim=2))

    tConverged = .false.

    nOrb = size(ci,dim=1)
    nLev = size(ci,dim=2)
    nAtom = size(iAtomStart) -1

    @:ASSERT(iAtomStart(nAtom+1)-1 == nOrb)
    @:ASSERT(all(shape(kpoint) == [3]))
    @:ASSERT(all(shape(S) == [nOrb,nOrb]))

    if (present(mIter)) then
      nIter = mIter
    else
      nIter = 20
    end if

    allocate(Sci1(nOrb))
    allocate(Sci2(nOrb,nLev))

    allocate(ciTmp1(nOrb))
    allocate(ciTmp2(nOrb))

    S = cmplx(0,0,dp)
    call unpackHS(S, over, kPoint, iNeighbour, nNeighbourSK, iCellVec, cellVec, iAtomStart, iPair,&
        & img2CentCell)

    lpLocalise: do iIter = 1, nIter

      write(stdout, *)'Iter', iIter

      alphamax = 0.0_dp

      ! sweep over all pairs of levels at that k-point
      do iLev1 = 1, nLev

        if (iLev1 < nLev) then
          call hemm(Sci2(1:nOrb,iLev1+1:nLev),'l',S,ci(:,iLev1+1:nLev), 'L')
        else
          call hemv(Sci2(1:nOrb,nLev),S,ci(:,nLev))
        end if

        do iLev2 = iLev1 +1, nLev

          call hemv(Sci1,S,ci(:,iLev1))

          Ast = 0.0_dp
          Bst = 0.0_dp
          do iAtom1 = 1, nAtom
            iOrb1 = iAtomStart(iAtom1)
            iOrb2 = iAtomStart(iAtom1+1) - 1
            Pst = 0.5_dp * (sum(ci(iOrb1:iOrb2,iLev1)*Sci2(iOrb1:iOrb2,iLev2))&
                & +sum(ci(iOrb1:iOrb2,iLev2)*Sci1(iOrb1:iOrb2)) )
            Pss = sum(ci(iOrb1:iOrb2,iLev1)*Sci1(iOrb1:iOrb2))
            Ptt = sum(ci(iOrb1:iOrb2,iLev2)*Sci2(iOrb1:iOrb2,iLev2))
            Ast = Ast + abs(Pst)**2 - 0.25_dp * abs(Pss - Ptt)**2
            Bst = Bst + real(Pst * (Pss - Ptt),dp)
          end do

          AB = Ast * Ast + Bst * Bst
          if (abs(AB)>0.0_dp) then
            C4A = -Ast / sqrt(AB)
            alpha = 0.25_dp * acos(C4A)
            if (Bst <= 0.0) then
              alpha = -alpha
            end if
          else
            alpha = 0.0_dp
          end if

          if (alphamax < abs(alpha)) then
            alphamax = alpha
          end if

          ! now we have to mix the two orbitals
          SINA=SIN(alpha)
          COSA=COS(alpha)
          ciTmp1 = ci(:,iLev1)
          ciTmp2 = ci(:,iLev2)
          ci(:,iLev1) = COSA*ciTmp1 + SINA*ciTmp2
          ci(:,iLev2) = COSA*ciTmp2 - SINA*ciTmp1

        end do
      end do

      conv = abs(alphamax) - abs(alphalast)
      if (iIter > 2 .and. ((abs(conv)<convergence) .or. alphamax == 0.0)) then
        tConverged = .true.
        exit
      end if
      alphalast = alphamax

    end do lpLocalise

    !write(stdout, *)'Localisations at each k-point'
    !write(stdout, "(6E12.4)") &
    !    & PipekMezyLocality_kpoint(ci, S, over, kpoint, kweights, &
    !    & iNeighbour, nNeighbourSK, iCellVec, cellVec, iAtomStart, iPair, &
    !    & img2CentCell)
    !write(stdout, "(1X,A,E12.4)")'Total', &
    !    & sum(PipekMezyLocality_kpoint(ci, S, over, kpoint, kweights, &
    !    & iNeighbour, nNeighbourSK, iCellVec, cellVec, iAtomStart, iPair, &
    !    & img2CentCell))

    if (.not.tConverged) then
      call warning("Exceeded iterations in Pipek-Mezey localisation!")
    end if

    ! Choose phase to make largest Bloch state element real for this k-point - dumb approch, as
    ! phases should be set globally for levels and k-points
    do iLev1 = 1, nLev
      iLoc = maxloc(abs(ci(:,iLev1)))
      ii = iLoc(1)
      phase = exp(-im * atan2(aimag(ci(ii,iLev1)), real(ci(ii,iLev1))))
      ci(:,iLev1) = phase * ci(:,iLev1)
    end do

  end subroutine PipekMezeyOld_kpoint

end module dftbp_pmlocalisation
