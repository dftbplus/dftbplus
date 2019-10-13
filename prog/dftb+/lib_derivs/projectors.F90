!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Dense projectors of states onto the empty manifold for BLACS matrices
module dftbp_projectBLACS
  use dftbp_accuracy, only : dp
  use dftbp_message
  use dftbp_commontypes
  use dftbp_blasroutines
  use dftbp_periodic, only : TNeighbourList
  use dftbp_densedescr
#:if WITH_SCALAPACK
  use dftbp_scalapackfx
  use dftbp_blacsenv
#:endif

  implicit none

  public

contains

#:if WITH_SCALAPACK

  !> Real projector onto conduction band
  subroutine projEmptyReal(myBlacs, desc, parallelKS, SSqrReal, nFilled, eigvecs, over,&
      & orb, species, iNeighbor, nNeighbor, iAtomStart, iPair,  img2CentCell, Pc)

    !> BLACS grid involved in calculation
    type(blacsgrid), intent(in) :: myBlacs

    !> Descriptor for dense matrix
    integer, intent(in) :: desc(DLEN_)

    !> K-points and spins to be handled
    type(TParallelKS), intent(in) :: parallelKS

    real(dp), intent(inout) :: SSqrReal(:,:)
    integer, intent(in) :: nFilled(:)
    real(dp), intent(in) :: eigvecs(:,:,:)
    real(dp), intent(in) :: over(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(0:,:)
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(out) :: Pc(:,:,:)

    integer :: iS, iKS
    integer :: nOrbs, nSparse
    integer :: ii, jj, iGlob, jGlob
    real(dp) :: SPc(size(Pc,dim=1),size(Pc,dim=2),size(Pc,dim=3))
    
    Pc = 0.0_dp
    SPc = 0.0_dp

    ! Pc = c . cT
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      call pblasfx_psyrk(eigvecs(:,:,iS), desc, Pc(:,:,iS), desc, kk=nFilled(iS))
    end do

    ! fill in other triangle

    ! transpose and add
    SPc = Pc
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      call pblasfx_ptran(Pc(:,:,iS),desc,SPc(:,:,iS),desc,beta=1.0_dp)
    end do
    ! correct diagonal, as is double included in transpose above
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      do jj = 1, size(SPc,dim=2)
        jGlob = scalafx_indxl2g(jj,desc(NB_), myBlacs%mycol,desc(CSRC_), myBlacs%ncol)
        do ii = 1, size(SPc,dim=1)
          iGlob = scalafx_indxl2g(ii,desc(MB_), myBlacs%myrow,desc(RSRC_), myBlacs%nrow)
          if (iGlob == jGlob) then
            SPc(ii,jj,iS) = SPc(ii,jj,iS) * 0.5_dp
          end if
        end do
      end do
    end do

    ! Pc = matmul(SSqr,Pc)
    call unpackhs_parallel_real(myBlacs, over, iNeighbor, nNeighbor, iAtomStart, iPair,&
        & img2CentCell, desc, SSqrReal)
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      call pblasfx_psymm(SSqrReal, desc, SPc(:,:,iS), desc, Pc(:,:,iS), desc)
    end do

    ! diag(Pc) = diag(Pc) - 1
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      do jj = 1, size(Pc,dim=2)
        jGlob = scalafx_indxl2g(jj,desc(NB_), myBlacs%mycol,desc(CSRC_), myBlacs%ncol)
        do ii = 1, size(Pc,dim=1)
          iGlob = scalafx_indxl2g(ii,desc(MB_), myBlacs%myrow,desc(RSRC_), myBlacs%nrow)
          if (iGlob == jGlob) then
            Pc(ii,jj,iS) = Pc(ii,jj,iS) - 1.0_dp
          end if
        end do
      end do
    end do
    
    Pc = -Pc ! now have (1 - S.Pv)
    
  end subroutine projEmptyReal
  
  !> Complex projector onto conduction band
  subroutine projEmptyCplx(myBlacs, desc, SSqr, nFilled, eigvecs, over, kpoint, iCellVec,&
      & cellVec, orb, species, iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, Pc)
    type(blacsgrid), intent(in) :: myBlacs
    integer, intent(in) :: desc(DLEN_)
    complex(dp), intent(inout) :: SSqr(:,:)
    integer, intent(in) :: nFilled(:)
    complex(dp), intent(in) :: eigvecs(:,:,:)
    real(dp), intent(in) :: over(:)
    real(dp), intent(in) :: kpoint(3)
    integer, intent(in) :: iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(0:,:)
    integer, intent(in) :: img2CentCell(:)
    complex(dp), intent(out) :: Pc(:,:,:)

    integer :: iS, iKS
    integer :: nOrbs, nSparse
    integer :: ii, jj, iGlob, jGlob
    complex(dp) :: SPc(size(Pc,dim=1),size(Pc,dim=2),size(Pc,dim=3))
    
    Pc = 0.0_dp
    SPc = 0.0_dp
    
    ! Pc = c . cT
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      call pblasfx_pherk(eigvecs(:,:,iS), desc, Pc(:,:,iS), desc, kk=nFilled(iS))
    end do
    
    ! fill in other triangle
    
    ! Hermitian transpose and add
    SPc = Pc
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      call pblasfx_ptranc(Pc(:,:,iS),desc,SPc(:,:,iS),desc,beta=(1.0_dp,0.0_dp))
    end do
    ! correct diagonal, as is double included in transpose above
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      do jj = 1, size(SPc,dim=2)
        jGlob = scalafx_indxl2g(jj,desc(NB_), myBlacs%mycol,desc(CSRC_), myBlacs%ncol)
        do ii = 1, size(SPc,dim=1)
          iGlob = scalafx_indxl2g(ii,desc(MB_), myBlacs%myrow,desc(RSRC_), myBlacs%nrow)
          if (iGlob == jGlob) then
            SPc(ii,jj,iS) = SPc(ii,jj,iS) * 0.5_dp
          end if
        end do
      end do
    end do
    
    ! Pc = matmul(SSqr,Pc)
    call unpackhs_parallel_cplx(myBlacs, over, kPoint, iNeighbor, nNeighbor, iCellVec, cellVec,&
        & iAtomStart, iPair, img2CentCell, desc, SSqr)
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      call pblasfx_phemm(SSqr, desc, SPc(:,:,iS), desc, Pc(:,:,iS), desc)
    end do
    
    ! diag(Pc) = diag(Pc) - 1
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      do jj = 1, size(Pc,dim=2)
        jGlob = scalafx_indxl2g(jj,desc(NB_), myBlacs%mycol,desc(CSRC_), myBlacs%ncol)
        do ii = 1, size(Pc,dim=1)
          iGlob = scalafx_indxl2g(ii,desc(MB_), myBlacs%myrow,desc(RSRC_), myBlacs%nrow)
          if (iGlob == jGlob) then
            Pc(ii,jj,iS) = Pc(ii,jj,iS) - 1.0_dp
          end if
        end do
      end do
    end do
    
    Pc = -Pc ! now have (1 - S.Pv)
    
  end subroutine projEmptyCplx
  
  !> Pauli matrix projection
  subroutine projEmptyPauli(myBlacs, desc, SSqr, nFilled, eigvecs, over, kpoint, iCellVec,&
      & cellVec, orb, mOrb, species, iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, Pc)
    type(blacsgrid), intent(in) :: myBlacs
    integer, intent(in) :: desc(DLEN_)
    complex(dp), intent(inout) :: SSqr(:,:)
    integer, intent(in) :: nFilled(:)
    complex(dp), intent(in) :: eigvecs(:,:,:)
    real(dp), intent(in) :: over(:)
    real(dp), intent(in) :: kpoint(3)
    integer, intent(in) :: iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: mOrb
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(0:,:)
    integer, intent(in) :: img2CentCell(:)
    complex(dp), intent(out) :: Pc(:,:,:)
    
    integer :: iS, iKS
    integer :: nOrbs, nSparse
    integer :: ii, jj, iGlob, jGlob
    complex(dp) :: SPc(size(Pc,dim=1),size(Pc,dim=2),size(Pc,dim=3))
    
    Pc = 0.0_dp
    SPc = 0.0_dp
    
    ! Pc = c . cT
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      call pblasfx_pherk(eigvecs(:,:,iS), desc, Pc(:,:,iS), desc, kk=nFilled(iS))
    end do
    
    ! fill in other triangle
    
    ! Hermitian transpose and add
    SPc = Pc
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      call pblasfx_ptranc(Pc(:,:,iS),desc,SPc(:,:,iS),desc,beta=(1.0_dp,0.0_dp))
    end do
    ! correct diagonal, as is double included in transpose above
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      do jj = 1, size(SPc,dim=2)
        jGlob = scalafx_indxl2g(jj,desc(NB_), myBlacs%mycol,desc(CSRC_), myBlacs%ncol)
        do ii = 1, size(SPc,dim=1)
          iGlob = scalafx_indxl2g(ii,desc(MB_), myBlacs%myrow,desc(RSRC_), myBlacs%nrow)
          if (iGlob == jGlob) then
            SPc(ii,jj,iS) = SPc(ii,jj,iS) * 0.5_dp
          end if
        end do
      end do
    end do
    
    ! Pc = matmul(SSqr,Pc)
    call unpacks_parallel_pauli(myBlacs, over, kPoint, iNeighbor, nNeighbor, iCellVec, cellVec,&
        & iAtomStart, iPair, img2CentCell, mOrb, desc, SSqr)
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      call pblasfx_phemm(SSqr, desc, SPc(:,:,iS), desc, Pc(:,:,iS), desc)
    end do
    
    ! diag(Pc) = diag(Pc) - 1
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)
      do jj = 1, size(Pc,dim=2)
        jGlob = scalafx_indxl2g(jj,desc(NB_), myBlacs%mycol,desc(CSRC_), myBlacs%ncol)
        do ii = 1, size(Pc,dim=1)
          iGlob = scalafx_indxl2g(ii,desc(MB_), myBlacs%myrow,desc(RSRC_), myBlacs%nrow)
          if (iGlob == jGlob) then
            Pc(ii,jj,iS) = Pc(ii,jj,iS) - 1.0_dp
          end if
        end do
      end do
    end do
    
    Pc = -Pc ! now have (1 - S.Pv)
    
  end subroutine projEmptyPauli

#:else

  !> Real projector onto conduction band
  subroutine projEmptyReal( Pc, ci, over, iNeighbor, nNeighbor, iPair, img2CentCell, denseDesc)

    !> Projector onto conduction band states
    real(dp), intent(out) :: Pc(:,:,:)

    !> Ground state wavefunctions
    real(dp), intent(in) :: ci(:,:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iPair(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    real(dp), allocatable :: SSqrReal(:,:)


    ! Pc = matmul(ci(:,:nFilled(1),1),transpose(ci(:,:nFilled(1),1)))
    Pc = 0.0_dp
    do iSpin = 1, nSpin
      call herk(Pc(:,:,iSpin),ci(:,:nFilled(1),iSpin))
      do ii = 1, size(ci,dim=2)
        Pc(ii,ii+1:,iSpin) = Pc(ii+1:,ii,iSpin)
      end do
    end do

    allocate(SSqrReal(size(Pc,dim=1), size(Pc,dim=2))

    SSqrReal = 0.0_dp
    call unpackHS(SSqrReal,over,neighborList%iNeighbor, nNeighbor, iAtomStart,iPair,img2CentCell)
    do ii = 1, size(ci,dim=2)
      SSqrReal(ii,ii+1:) = SSqrReal(ii+1:,ii)
    end do

    write(*,*)'Forming Pc and PcT'
    do iSpin = 1, nSpin
      call gemm(Pct(:,:,iSpin),Pc(:,:,iSpin),SSqrReal,transA='T')
    end do
    !PcT = matmul(Pc,SSqrReal)
    allocate(arrayTmp(size(ci,dim=1),size(ci,dim=2)))
    !Pc = matmul(SSqrReal,Pc)
    do iSpin = 1, nSpin
      call gemm(arrayTmp,SSqrReal,Pc(:,:,iSpin),transA='T')
      Pc(:,:,iSpin) = arrayTmp
    end do
    deallocate(arrayTmp)

    do iSpin = 1, nSpin
      do ii = 1, size(ci,dim=2)
        Pc(ii,ii,iSpin) = Pc(ii,ii,iSpin) - 1.0_dp
        PcT(ii,ii,iSpin) = PcT(ii,ii,iSpin) - 1.0_dp
      end do
    end do

  end subroutine projEmptyReal



#:endif

end module projectors
