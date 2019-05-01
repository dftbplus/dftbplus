!> Module for linear respose derivative calculations
!!
!! \details only singlet perturbation sums over states expression
!! implemented in parallel for the moment
module LRderivs
# include "assert.h"
# include "allocate.h"
  use accuracy
  use constants
  use message
  use CommonTypes
  use Potentials
  use scc
  use OrbitalEquiv
  use populations
  use spin
  use mainio
  use shift
  use mixer
  use projectors
  use parallel
  use finiteT
  use scalapackfx_extlib
  use mpiglobal
  use mpifx_extlib
  use timer_module
  implicit none

  private
  public :: perturbation
  
  interface perturbation
    module procedure perturbStat
    module procedure perturbDyn
  end interface perturbation
  
  !interface sternheimer
  !  module procedure sternheimerSingletDense
  !end interface sternheimer
  
  complex(dp), parameter :: eta = (0.0_dp,1.0E-8_dp) ! small complex
  ! value for frequency dependent cases
  
  character(len=1), parameter :: direction(3) = ['x','y','z']
  
contains
  
  subroutine perturbStat(mympi, allproc, grpproc, gridAtom, &
      & groupKS, nGroups, groupsize, desc, filling, nEl, SSqrReal, eigvals, &
      & eigvecs, ham, over, orb, nAtom, species, speciesnames, iNeighbor, &
      & nNeighbor, iAtomStart, iPair, img2CentCell, coord, nSCCIter, sccTol, &
      & nMixElements, nIneqOrb, iEqOrbitals, kT, Ef, W, pChrgMixer)
    type(mpifx_comm), intent(in)   :: mympi
    type(blacsgrid), intent(in)    :: allproc
    type(blacsgrid), intent(in)    :: grpproc
    type(blacsgrid), intent(in)    :: gridAtom
    integer, intent(in)            :: groupKS(:,:)
    integer, intent(in)            :: nGroups, groupsize
    integer, intent(in)            :: desc(DLEN_)
    real(dp), intent(in)           :: filling(:,:,:)
    real(dp), intent(in)           :: nEl(:)
    real(dp), intent(inout)        :: SSqrReal(:,:)
    real(dp), intent(in)           :: eigvals(:,:,:)
    real(dp), intent(in)           :: eigvecs(:,:,:)
    real(dp), intent(in)           :: ham(:,:), over(:)
    type(TOrbitals), pointer       :: orb
    integer, intent(in)            :: nAtom, species(:)
    character(mc), intent(in)      :: speciesnames(:)
    integer, intent(in)            :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in)            :: iAtomStart(:), iPair(0:,:)
    integer, intent(in)            :: img2CentCell(:)
    real(dp), intent(in)           :: coord(:,:)
    integer, intent(in)            :: nSCCIter
    real(dp), intent(in)           :: sccTol
    integer, intent(in)            :: nMixElements, nIneqOrb
    integer, intent(in)            :: iEqOrbitals(:,:,:)
    real(dp), intent(in)           :: kT, Ef(:)
    real(dp), intent(in)           :: W(:,:,:)
    type(omixer), pointer, intent(inout) :: pChrgMixer
    
    type(timer) :: myclock_all
    integer     :: iS, iK, iKS, iAt, iCart, iSCC, iLev, iSh, iSp
    integer     :: nSpin, nKS, nOrbs, nSparse
    integer, allocatable  :: nFilled(:), nEmpty(:)
    
    integer :: ii, jj, iGlob, jGlob
    integer :: iSCCIter
    logical :: tStopSCC
    
    ! derivatives of hamiltonian, eigenvalues and eigenvectors
    real(dp), allocatable :: dham(:,:)
    real(dp) :: drho(size(over),size(ham, dim=2))
    real(dp) :: drhoExtra(size(over),size(ham, dim=2))
    real(dp) :: dqIn(orb%mOrb,nAtom,size(ham, dim=2))
    real(dp) :: dqOut(orb%mOrb,nAtom,size(ham, dim=2))
    real(dp) :: dqInpRed(nMixElements), dqOutRed(nMixElements)
    real(dp) :: dqDiffRed(nMixElements), sccErrorQ
    real(dp) :: dqPerShell(orb%mOrb,nAtom,size(ham, dim=2))
    real(dp) :: tmpFill(size(filling,dim=1))
    type(TPotentials) :: dpotential
    
    real(dp) :: dci(size(eigvecs,dim=1),size(eigvecs,dim=2), &
        & size(eigvecs,dim=3))
    real(dp) :: deigvals(size(eigvals,dim=1),size(eigvals,dim=2), &
        & size(eigvals,dim=3))
    real(dp) :: work(size(eigvecs,dim=1),size(eigvecs,dim=2), &
        & size(eigvecs,dim=3))
    real(dp) :: work2(size(eigvecs,dim=1),size(eigvecs,dim=2), &
        & size(eigvecs,dim=3))
    
    real(dp) :: nF(size(ham, dim=2)), dEf(size(ham, dim=2))

    logical :: tSCC, tSpin, tMetallic, tConverged
    
    nSpin = size(ham, dim=2)
    nKS = size(groupKS, dim=2)
    nSparse = size(ham,dim=1)
    nOrbs = size(filling,dim=1)
        
    tSCC = (nSCCIter > 1)
    tSpin = (nSpin > 1)
    
    ASSERT(size(eigvecs,dim=3) == nSpin)
    ASSERT(all(shape(coord) == [3,nAtom]))
    
    call create(dpotential,orb,nAtom,nSpin)
    
    call myclock_all%start(1)
    
    ALLOCATE_(nFilled,(nSpin))
    ALLOCATE_(nEmpty,(nSpin))
    
    nFilled = -1
    do iS = 1, nSpin
      do iLev = 1, nOrbs
        if ( filling(iLev,1,iS) <= epsilon(1.0_dp) ) then
          nFilled(iS) = iLev - 1
          exit
        end if
      end do
      if (nFilled(iS) < 0) then
        nFilled(iS) = nOrbs
      end if
    end do
    nEmpty = -1
    do iS = 1, nSpin
      do iLev = 1, nOrbs
        if ( abs( filling(iLev,1,iS) - real(3-nSpin,dp) ) &
            & >= epsilon(1.0_dp)) then
          nEmpty(iS) = iLev
          exit
        end if
      end do
      if (nEmpty(iS) < 0) then
        nEmpty(iS) = 1
      end if
    end do
    
    tMetallic = (.not.all(nFilled == nEmpty -1))
    
    if (ioproc) then
      write(*,*)'Fully or partly filled states',nFilled
      write(*,*)'Fully or partly empty states',nEmpty
      if (tMetallic) then
        write(*,*)'Metallic system'
      else
        write(*,*)'Non-metallic system'
      end if      
    end if
    
    if (tMetallic) then
      
      ! Number of electrons at the Fermi energy
      do iS = 1, nSpin
        nf(iS) = 0.0_dp
        do ii = nEmpty(iS), nFilled(iS)
          nf(iS) = nf(iS) &
              & + deltamn(Ef(iS),eigvals(ii,1,iS),kT)
        end do
        nf(iS) = real(3-nSpin,dp)*nf(iS)
      end do
      if (ioproc) then
        write(*,*)'Electrons at the Fermi energy Nf:',nF
      end if
    end if
    
    ALLOCATE_(dham,(nSparse,nSpin))
    
    do iCart = 1, 3 ! polarization direction
      
      dqIn = 0.0_dp
      dqOut = 0.0_dp
      ! derivative of E.x
      dpotential%extAtom = 0.0_dp
      do iAt = 1, nAtom
        dpotential%extAtom(iAt,1) = coord(iCart,iAt)
      end do
      
      dpotential%extBlock = 0.0_dp
      call total_shift(dpotential%extBlock, dpotential%extAtom, orb, species)
      
      if (tSCC .and. ioproc) then
        call reset(pChrgMixer, nMixElements)
        dqInpRed = 0.0_dp
        dqPerShell = 0.0_dp
      end if
      
      iSCCIter = 1
      tStopSCC = .false.
      lpSCC: do while (iSCCiter <= nSCCIter)
        
        dpotential%intAtom = 0.0_dp
        dpotential%intBlock = 0.0_dp
        dpotential%intShell = 0.0_dp
        
        if (tSCC .and. iSCCiter>1) then
          call updateCharges_SCC(allproc, gridAtom, dqIn, orb, species, &
              & iNeighbor, img2CentCell)
          call getShiftPerAtom(dpotential%intAtom)
          call getShiftPerL(dpotential%intShell)
          
          if (tSpin) then
            call addSpinShift(dpotential%intShell,dqPerShell,species,orb,W)
          end if
          
          !if (tDFTBU) then
          !  call shift_DFTBU(dpotential%orbitalBlock,qBlockIn,specie,orb, &
          !      & 2, UJ, nUJ, niUJ, iUJ)
          !end if
          
        end if
        
        call total_shift(dpotential%intShell,dpotential%intAtom, orb,species)
        call total_shift(dpotential%intBlock,dpotential%intShell, orb,species)
        
        dpotential%intBlock = dpotential%intBlock + dpotential%extBlock
        
        dham = 0.0_dp
        call add_shift(dham, over, nNeighbor, iNeighbor, species, orb, iPair, &
            & nAtom, img2CentCell, dpotential%intBlock)
        
        if (nSpin == 2) then
          dham(:,:) = 2.0_dp *  dham(:,:)
        end if
        call qm2ud(dham)
        
        drho = 0.0_dp
        
        if (grpproc%iproc /= -1) then
          
          do iKS = 1, nKS
            iK = groupKS(1, iKS)
            iS = groupKS(2, iKS)
            
            call unpackhs_parallel_real(grpproc, dham(:,iS), iNeighbor, &
                & nNeighbor, iAtomStart, iPair, img2CentCell, desc, &
                & work(:,:,iKS))
            
            ! dH times c_i
            call pblasfx_psymm(work(:,:,iKS), desc, eigvecs(:,:,iKS), desc, &
                & work2(:,:,iKS), desc) !, mm=nFilled(iS))
            
            ! c_i times dH times c_i
            call pblasfx_pgemm(eigvecs(:,:,iKS), desc, work2(:,:,iKS), desc, &
                & work(:,:,iKS), desc, transa="T")
            
            ! derivative of eigenvalues stored diagonal of matrix work,
            ! from <c|h'|c>            
            do jj = 1, size(work,dim=2)
              jGlob = scalafx_indxl2g(jj,desc(NB_), grpproc%mycol,desc(CSRC_),&
                  & grpproc%ncol)
              do ii = 1, size(work,dim=1)
                iGlob = scalafx_indxl2g(ii,desc(MB_),grpproc%myrow,desc(RSRC_),&
                    & grpproc%nrow)
                
                !if (iGlob == jGlob) then ! derivative of eigenvalues
                !  write(*,*)iGlob,work(ii,jj,iKS)*Hartree__eV
                !end if
                
                ! weight with inverse of energy differences
                work(ii,jj,iKS) = work(ii,jj,iKS) * &
                    &invDiff(eigvals(jGlob,1,iS),eigvals(iGlob,1,iS),Ef(iS),kT)&
                    & * theta(eigvals(jGlob,1,iS),eigvals(iGlob,1,iS),kT)
              end do
            end do
            
            ! Derivatives of states
            call pblasfx_pgemm(eigvecs(:,:,iKS), desc, work(:,:,iKS), desc, &
                & work2(:,:,iKS), desc)
            
            ! Form derivative of density matrix and symmetrize
            call pblasfx_pgemm(work2(:,:,iKS), desc,eigvecs(:,:,iKS), desc, &
                & work(:,:,iKS), desc, transb="T", kk=nFilled(iS))
            work2(:,:,iKS) = work(:,:,iKS)
            
            call pblasfx_ptran(work(:,:,iKS),desc,work2(:,:,iKS),desc, &
                & beta=1.0_dp)
            
            call packrho_parallel_real(grpproc, desc, work2(:,:,iKS), &
                & iNeighbor, nNeighbor, orb%mOrb, iAtomStart, iPair, &
                & img2CentCell, drho(:,iS) )
            
          end do
          
          
        end if
        
        call blacsfx_gsum(allproc, drho, rdest=-1, cdest=-1)
        
        if (tMetallic.and.grpproc%iproc /= -1) then
          ! correct for Fermi level shift for q=0 fields
          drhoextra = 0.0_dp
          do iKS = 1, nKS
            iK = groupKS(1, iKS)
            iS = groupKS(2, iKS)
            
            dqOut = 0.0_dp
            call mulliken(dqOut(:,:,iS), over, drho(:,iS), orb, &
                & iNeighbor, nNeighbor, img2CentCell, iPair)
            
            dEf(iS) = -sum(dqOut(:,:,iS))/nF(iS)
            
            if (abs(dEf(iS)) > epsilon(1.0_dp)) then ! Fermi level changes,
              ! so need to correct for change in the number of charges
              work(:,:,iKS) = 0.0_dp
              do jj = 1, size(work,dim=2)
                jGlob = scalafx_indxl2g(jj,desc(NB_), &
                    & grpproc%mycol,desc(CSRC_), grpproc%ncol)
                if (jGlob >= nEmpty(iS) .and. jGlob <= nFilled(iS)) then
                  work(:,jj,iKS) = eigvecs(:,jj,iKS) * &
                      & deltamn(eigVals(jGlob,1,iKS),Ef(iS),kT)*dEf(iS)
                end if
              end do
              call pblasfx_pgemm(work(:,:,iKS), desc,eigvecs(:,:,iKS), desc, &
                  & work2(:,:,iKS), desc, transb="T",alpha=2.0_dp)
              call packrho_parallel_real(grpproc, desc, work2(:,:,iKS), &
                  & iNeighbor, nNeighbor, orb%mOrb, iAtomStart, iPair, &
                  & img2CentCell, drhoextra(:,iS) )
            end if
            
          end do
          
          call blacsfx_gsum(allproc, drhoextra, rdest=-1, cdest=-1)
          drho = drho + drhoextra
        end if
        
        if (nSpin == 1) then
          drho = 2.0_dp * drho
        end if
        
        call ud2qm(drho)
        
        dqOut = 0.0_dp
        do iS = 1, nSpin
          call mulliken(dqOut(:,:,iS), over, drho(:,iS), orb, &
              & iNeighbor, nNeighbor, img2CentCell, iPair)
          !if (tDFTBU) then
          !  qBlockOut(:,:,:,iS) = 0.0_dp
          !  call mulliken(qBlockOut(:,:,:,iS), over, drho(:,iS), &
          !      &orb, iNeighbor, nNeighbor, img2CentCell, iPair)
          !end if
        end do
        
        if (tSCC) then
          dqOutRed = 0.0_dp
          call OrbitalEquiv_reduce(dqOut, iEqOrbitals, orb, &
              & dqOutRed(:nIneqOrb))
          !if (tDFTBU) then
          !  call AppendBlock_reduce( qBlockOut, iEqBlockDFTBU, orb, &
          !      & qOutRed )
          !end if
          
          dqDiffRed(:) = dqOutRed(:) - dqInpRed(:)
          sccErrorQ = maxval(abs(dqDiffRed))
          
          if (ioproc) then
            write(*,*)'Iter',iSCCIter,'Error',sccErrorQ
          end if
          tConverged = (sccErrorQ < sccTol)
          
          if ((.not. tConverged) .and. iSCCiter /= nSCCiter) then
            if (iSCCIter == 1) then
              dqIn(:,:,:) = dqOut(:,:,:)
              dqInpRed(:) = dqOutRed(:)
              !if (tDFTBU) then
              !  qBlockIn(:,:,:,:) = qBlockOut(:,:,:,:)
              !end if
            else
              ! Mixing only done on master node as mixer may require IO.
              if (ioproc) then
                call mix(pChrgMixer, dqInpRed, dqDiffRed)
              end if
              call mpifx_bcast(mympi, dqInpRed, root=ioprocid)
              
              call OrbitalEquiv_expand(dqInpRed(:nIneqOrb), iEqOrbitals, &
                  & orb, dqIn)
              !if (tDFTBU) then
              !  qBlockIn = 0.0_dp
              !  call Block_expand( qInpRed ,iEqBlockDFTBU, orb, &
              !      & qBlockIn, specie0, nUJ, niUJ, iUJ, orbEquiv=iEqOrbitals )
              !end if
            end if
          end if
        else
          tConverged = .true.
        end if
        
        if (tConverged) then
          exit lpSCC
        end if
        
        if (tSpin) then
          dqPerShell = 0.0_dp
          do iAt = 1, nAtom
            iSp = species(iAt)
            do iSh = 1, orb%nShell(iSp)
              dqPerShell(iSh,iAt,:nSpin) = dqPerShell(iSh,iAt,:nSpin) +&
                  & sum(dqIn(orb%posShell(iSh,iSp): &
                  & orb%posShell(iSh+1,iSp)-1,iAt,:nSpin),dim=1)
            end do
          end do
          
        end if
        
        iSCCIter = iSCCIter +1
        
      end do lpSCC
      
      if (ioproc) then
        if (iCart==1) then
          write(*,*)'Polarisability'
        end if
        do ii = 1, 3
          write(*,"(E20.12)",advance ='no') &
              & -sum(sum(dqOut(:,:nAtom,1),dim=1)*coord(ii,:nAtom))
        end do
        write(*,*)
      end if
      
    end do
    
    call destroy(dpotential)
    
    DEALLOCATE_(dham)
    DEALLOCATE_(nFilled)
    DEALLOCATE_(nEmpty)
    
    call myclock_all%stop()
    call myclock_all%print_times(msg="linear response total", fp=stdout)

  end subroutine perturbStat
  
  subroutine perturbDyn(mympi, allproc, grpproc, gridAtom, &
      & groupKS, nGroups, groupsize, desc, filling, nEl, SSqrReal, eigvals, &
      & eigvecs, ham, over, orb, nAtom, species, speciesnames, iNeighbor, &
      & nNeighbor, iAtomStart, iPair, img2CentCell, coord, nSCCIter, sccTol, &
      & nMixElements, nIneqOrb, iEqOrbitals, kT, Ef, W, pChrgMixer, omega)
    type(mpifx_comm), intent(in)   :: mympi
    type(blacsgrid), intent(in)    :: allproc
    type(blacsgrid), intent(in)    :: grpproc
    type(blacsgrid), intent(in)    :: gridAtom
    integer, intent(in)            :: groupKS(:,:)
    integer, intent(in)            :: nGroups, groupsize
    integer, intent(in)            :: desc(DLEN_)
    real(dp), intent(in)           :: filling(:,:,:)
    real(dp), intent(in)           :: nEl(:)
    real(dp), intent(inout)        :: SSqrReal(:,:)
    real(dp), intent(in)           :: eigvals(:,:,:)
    real(dp), intent(in)           :: eigvecs(:,:,:)
    real(dp), intent(in)           :: ham(:,:), over(:)
    type(TOrbitals), pointer       :: orb
    integer, intent(in)            :: nAtom, species(:)
    character(mc), intent(in)      :: speciesnames(:)
    integer, intent(in)            :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in)            :: iAtomStart(:), iPair(0:,:)
    integer, intent(in)            :: img2CentCell(:)
    real(dp), intent(in)           :: coord(:,:)
    integer, intent(in)            :: nSCCIter
    real(dp), intent(in)           :: sccTol
    integer, intent(in)            :: nMixElements, nIneqOrb
    integer, intent(in)            :: iEqOrbitals(:,:,:)
    real(dp), intent(in)           :: kT, Ef(:)
    real(dp), intent(in)           :: W(:,:,:)
    type(omixer), pointer, intent(inout) :: pChrgMixer
    real(dp), intent(in)           :: omega(:)

    type(timer) :: myclock_all
    integer     :: iS, iK, iKS, iAt, iCart, iSCC, iLev, iSh, iSp
    integer     :: nSpin, nKS, nOrbs, nSparse
    integer, allocatable  :: nFilled(:), nEmpty(:)
    
    integer :: ii, jj, iGlob, jGlob, iSignOmega, iOmega, nOmega
    integer :: iSCCIter
    logical :: tStopSCC
    
    ! derivatives of hamiltonian, eigenvalues and eigenvectors
    real(dp), allocatable :: dham(:,:)
    real(dp) :: drho(size(over),size(ham, dim=2))
    real(dp) :: drhoExtra(size(over),size(ham, dim=2))
    real(dp) :: dqIn(orb%mOrb,nAtom,size(ham, dim=2))
    real(dp) :: dqOut(orb%mOrb,nAtom,size(ham, dim=2))
    real(dp) :: dqInpRed(nMixElements), dqOutRed(nMixElements)
    real(dp) :: dqDiffRed(nMixElements), sccErrorQ, polarisability(3,3)
    real(dp) :: dqPerShell(orb%mOrb,nAtom,size(ham, dim=2))
    real(dp) :: tmpFill(size(filling,dim=1))
    type(TPotentials) :: dpotential
    
    real(dp) :: dci(size(eigvecs,dim=1),size(eigvecs,dim=2), &
        & size(eigvecs,dim=3))
    real(dp) :: deigvals(size(eigvals,dim=1),size(eigvals,dim=2), &
        & size(eigvals,dim=3))
    real(dp) :: work(size(eigvecs,dim=1),size(eigvecs,dim=2), &
        & size(eigvecs,dim=3))
    real(dp) :: work2(size(eigvecs,dim=1),size(eigvecs,dim=2), &
        & size(eigvecs,dim=3))
    complex(dp) :: cwork1(size(eigvecs,dim=1),size(eigvecs,dim=2))
    complex(dp) :: tmpCeigvecs(size(eigvecs,dim=1),size(eigvecs,dim=2))    
    complex(dp) :: cwork2(size(eigvecs,dim=1),size(eigvecs,dim=2))    

    real(dp) :: nF(size(ham, dim=2)), dEf(size(ham, dim=2))

    logical :: tSCC, tSpin, tMetallic, tConverged
    
    nSpin = size(ham, dim=2)
    nKS = size(groupKS, dim=2)
    nSparse = size(ham,dim=1)
    nOrbs = size(filling,dim=1)
    nOmega = size(omega)
        
    tSCC = (nSCCIter > 1)
    tSpin = (nSpin > 1)
    
    ASSERT(size(eigvecs,dim=3) == nSpin)
    ASSERT(all(shape(coord) == [3,nAtom]))

    if (ioproc) then
      do iCart = 1, 3
        open(1001+iCart,file="pol_"//direction(iCart),position="rewind", &
            & status="replace")
      end do
    end if
    
    call create(dpotential,orb,nAtom,nSpin)
    
    call myclock_all%start(1)
    
    ALLOCATE_(nFilled,(nSpin))
    ALLOCATE_(nEmpty,(nSpin))
    
    nFilled = -1
    do iS = 1, nSpin
      do iLev = 1, nOrbs
        if ( filling(iLev,1,iS) <= epsilon(1.0_dp) ) then
          nFilled(iS) = iLev - 1
          exit
        end if
      end do
      if (nFilled(iS) < 0) then
        nFilled(iS) = nOrbs
      end if
    end do
    nEmpty = -1
    do iS = 1, nSpin
      do iLev = 1, nOrbs
        if ( abs( filling(iLev,1,iS) - real(3-nSpin,dp) ) &
            & >= epsilon(1.0_dp)) then
          nEmpty(iS) = iLev
          exit
        end if
      end do
      if (nEmpty(iS) < 0) then
        nEmpty(iS) = 1
      end if
    end do
    
    tMetallic = (.not.all(nFilled == nEmpty -1))
    
    if (ioproc) then
      write(*,*)'Fully or partly filled states',nFilled
      write(*,*)'Fully or partly empty states',nEmpty
      if (tMetallic) then
        write(*,*)'Metallic system'
      else
        write(*,*)'Non-metallic system'
      end if      
    end if
    
    if (tMetallic) then
      
      ! Number of electrons at the Fermi energy
      do iS = 1, nSpin
        nf(iS) = 0.0_dp
        do ii = nEmpty(iS), nFilled(iS)
          nf(iS) = nf(iS) &
              & + deltamn(Ef(iS),eigvals(ii,1,iS),kT)
        end do
        nf(iS) = real(3-nSpin,dp)*nf(iS)
      end do
      if (ioproc) then
        write(*,*)'Electrons at the Fermi energy Nf:',nF
      end if
    end if
    
    ALLOCATE_(dham,(nSparse,nSpin))
    
    do iCart = 1, 3 ! polarization direction
      
      dqIn = 0.0_dp
      dqOut = 0.0_dp
      ! derivative of E.x
      dpotential%extAtom = 0.0_dp
      do iAt = 1, nAtom
        dpotential%extAtom(iAt,1) = coord(iCart,iAt)
      end do
      
      dpotential%extBlock = 0.0_dp
      call total_shift(dpotential%extBlock, dpotential%extAtom, orb, species)

      do iOmega = 1, nOmega
        
        if (tSCC .and. ioproc) then
          call reset(pChrgMixer, nMixElements)          
        end if
        if (iOmega == 1) then
          dqInpRed = 0.0_dp
          dqPerShell = 0.0_dp
        end if
        
        if (ioproc) then
          write(*,*)'At',omega(iOmega),'au',omega(iOmega)*Hartree__eV,'eV'
        end if
        
        iSCCIter = 1
        tStopSCC = .false.
        lpSCC: do while (iSCCiter <= nSCCIter)
          
          dpotential%intAtom = 0.0_dp
          dpotential%intBlock = 0.0_dp
          dpotential%intShell = 0.0_dp
          
          if (tSCC .and. ( iSCCiter>1 .or. iOmega > 1)) then
            call updateCharges_SCC(allproc, gridAtom, dqIn, orb, species, &
                & iNeighbor, img2CentCell)
            call getShiftPerAtom(dpotential%intAtom)
            call getShiftPerL(dpotential%intShell)
            
            if (tSpin) then
              call addSpinShift(dpotential%intShell,dqPerShell,species,orb,W)
            end if
            
            !if (tDFTBU) then
            !  call shift_DFTBU(dpotential%orbitalBlock,qBlockIn,specie,orb, &
            !      & 2, UJ, nUJ, niUJ, iUJ)
            !end if
            
          end if
          
          call total_shift(dpotential%intShell,dpotential%intAtom, orb,species)
          call total_shift(dpotential%intBlock,dpotential%intShell, orb,species)
          
          dpotential%intBlock = dpotential%intBlock + dpotential%extBlock
          
          dham = 0.0_dp
          call add_shift(dham, over, nNeighbor, iNeighbor, species, orb, iPair, &
              & nAtom, img2CentCell, dpotential%intBlock)
          
          if (nSpin == 2) then
            dham(:,:) = 2.0_dp *  dham(:,:)
          end if
          call qm2ud(dham)
          
          drho = 0.0_dp
          
          if (grpproc%iproc /= -1) then
            
            do iKS = 1, nKS
              iK = groupKS(1, iKS)
              iS = groupKS(2, iKS)
              
              call unpackhs_parallel_real(grpproc, dham(:,iS), iNeighbor, &
                  & nNeighbor, iAtomStart, iPair, img2CentCell, desc, &
                  & work(:,:,iKS))
              
              ! dH times c_i
              call pblasfx_psymm(work(:,:,iKS), desc, eigvecs(:,:,iKS), desc, &
                  & work2(:,:,iKS), desc) !, mm=nFilled(iS))
              
              ! c_i times dH times c_i
              call pblasfx_pgemm(eigvecs(:,:,iKS), desc, work2(:,:,iKS), desc, &
                  & work(:,:,iKS), desc, transa="T")
              ! derivative of eigenvalues stored diagonal of matrix work,
              ! from <c|h'|c>            
              
              tmpCeigvecs = eigvecs(:,:,iKS) ! complex copy of eigenvectors
              
              do iSignOmega = -1, 1, 2 ! loop over positive and negative frequencies
                
                do jj = 1, size(work,dim=2)
                  jGlob = scalafx_indxl2g(jj,desc(NB_),grpproc%mycol,desc(CSRC_),&
                      & grpproc%ncol)
                  do ii = 1, size(work,dim=1)
                    iGlob = scalafx_indxl2g(ii,desc(MB_),grpproc%myrow,desc(RSRC_),&
                        & grpproc%nrow)
                    
                    ! weight with inverse of energy differences
                    cwork1(ii,jj) = work(ii,jj,iKS) * &
                        & ( theta(eigvals(jGlob,1,iS),Ef(iS),kT) &
                        & - theta(eigvals(iGlob,1,iS),Ef(iS),kT) ) &
                        & * theta(eigvals(jGlob,1,iS),eigvals(iGlob,1,iS),kT) &
                        & / (eigvals(jGlob,1,iS)-eigvals(iGlob,1,iS) &
                        &    + iSignOmega * omega(iOmega) + eta)
                  end do
                end do
                
                ! Derivatives of eigen-state vectors
                call pblasfx_pgemm(tmpCeigvecs, desc, cwork1, desc, &
                    & cwork2, desc)
                
                ! Form derivative of density matrix and then hermitian symmetrise
                call pblasfx_pgemm(cwork2, desc,tmpCeigvecs, desc, &
                    & cwork1, desc, transb="T", kk=nFilled(iS))
                
                cwork2 = cwork1
                
                call pblasfx_ptranc(cwork2,desc,cwork1,desc, &
                    & alpha=(0.5_dp,0.0_dp), beta=(0.5_dp,0.0_dp))
                
                call packrho_parallel_real(grpproc, desc, real(cwork1), &
                    & iNeighbor, nNeighbor, orb%mOrb, iAtomStart, iPair, &
                    & img2CentCell, drho(:,iS) )
                
              end do
              
            end do
            
          end if
          
          call blacsfx_gsum(allproc, drho, rdest=-1, cdest=-1)
          
          if (tMetallic.and.grpproc%iproc /= -1) then
            ! correct for Fermi level shift for q=0 fields
            drhoextra = 0.0_dp
            do iKS = 1, nKS
              iK = groupKS(1, iKS)
              iS = groupKS(2, iKS)
              
              dqOut = 0.0_dp
              call mulliken(dqOut(:,:,iS), over, drho(:,iS), orb, &
                  & iNeighbor, nNeighbor, img2CentCell, iPair)
              
              dEf(iS) = -sum(dqOut(:,:,iS))/nF(iS)
              
              if (abs(dEf(iS)) > epsilon(1.0_dp)) then ! Fermi level changes,
                ! so need to correct for change in the number of charges
                work(:,:,iKS) = 0.0_dp
                do jj = 1, size(work,dim=2)
                  jGlob = scalafx_indxl2g(jj,desc(NB_), &
                      & grpproc%mycol,desc(CSRC_), grpproc%ncol)
                  if (jGlob >= nEmpty(iS) .and. jGlob <= nFilled(iS)) then
                    work(:,jj,iKS) = eigvecs(:,jj,iKS) * &
                        & deltamn(eigVals(jGlob,1,iKS),Ef(iS),kT)*dEf(iS)
                  end if
                end do
                call pblasfx_pgemm(work(:,:,iKS), desc,eigvecs(:,:,iKS), desc, &
                    & work2(:,:,iKS), desc, transb="T",alpha=2.0_dp)
                call packrho_parallel_real(grpproc, desc, work2(:,:,iKS), &
                    & iNeighbor, nNeighbor, orb%mOrb, iAtomStart, iPair, &
                    & img2CentCell, drhoextra(:,iS) )
              end if
              
            end do
            
            call blacsfx_gsum(allproc, drhoextra, rdest=-1, cdest=-1)
            drho = drho + drhoextra
          end if
          
          if (nSpin == 1) then
            drho = 2.0_dp * drho
          end if
          
          call ud2qm(drho)
          
          dqOut = 0.0_dp
          do iS = 1, nSpin
            call mulliken(dqOut(:,:,iS), over, drho(:,iS), orb, &
                & iNeighbor, nNeighbor, img2CentCell, iPair)
            !if (tDFTBU) then
            !  qBlockOut(:,:,:,iS) = 0.0_dp
            !  call mulliken(qBlockOut(:,:,:,iS), over, drho(:,iS), &
            !      &orb, iNeighbor, nNeighbor, img2CentCell, iPair)
            !end if
          end do
          
          if (tSCC) then
            dqOutRed = 0.0_dp
            call OrbitalEquiv_reduce(dqOut, iEqOrbitals, orb, &
                & dqOutRed(:nIneqOrb))
            !if (tDFTBU) then
            !  call AppendBlock_reduce( qBlockOut, iEqBlockDFTBU, orb, &
            !      & qOutRed )
            !end if
            
            dqDiffRed(:) = dqOutRed(:) - dqInpRed(:)
            sccErrorQ = maxval(abs(dqDiffRed))
            
            if (ioproc) then
              write(*,*)'Iter',iSCCIter,'Error',sccErrorQ
            end if
            tConverged = (sccErrorQ < sccTol)
            
            if ((.not. tConverged) .and. iSCCiter /= nSCCiter) then
              if (iSCCIter == 1 .and. iOmega == 1) then
                dqIn(:,:,:) = dqOut(:,:,:)
                dqInpRed(:) = dqOutRed(:)
                !if (tDFTBU) then
                !  qBlockIn(:,:,:,:) = qBlockOut(:,:,:,:)
                !end if
              else
                ! Mixing only done on master node as mixer may require IO.
                if (ioproc) then
                  call mix(pChrgMixer, dqInpRed, dqDiffRed)
                end if
                call mpifx_bcast(mympi, dqInpRed, root=ioprocid)
                
                call OrbitalEquiv_expand(dqInpRed(:nIneqOrb), iEqOrbitals, &
                    & orb, dqIn)
                !if (tDFTBU) then
                !  qBlockIn = 0.0_dp
                !  call Block_expand( qInpRed ,iEqBlockDFTBU, orb, &
                !      & qBlockIn, specie0, nUJ, niUJ, iUJ, orbEquiv=iEqOrbitals )
                !end if
              end if
            end if
          else
            tConverged = .true.
          end if
          
          if (tConverged) then
            exit lpSCC
          end if
          
          if (tSpin) then
            dqPerShell = 0.0_dp
            do iAt = 1, nAtom
              iSp = species(iAt)
              do iSh = 1, orb%nShell(iSp)
                dqPerShell(iSh,iAt,:nSpin) = dqPerShell(iSh,iAt,:nSpin) +&
                    & sum(dqIn(orb%posShell(iSh,iSp): &
                    & orb%posShell(iSh+1,iSp)-1,iAt,:nSpin),dim=1)
              end do
            end do
            
          end if
          
          iSCCIter = iSCCIter +1
          
        end do lpSCC
        
        if (ioproc) then
          if (iCart==1) then
            write(*,*)'Polarisability'
          end if
          do ii = 1, 3
            write(*,"(E20.12)",advance ='no') &
                & -sum(sum(dqOut(:,:nAtom,1),dim=1)*coord(ii,:nAtom))
            polarisability(ii,iCart) = &
                & -sum(sum(dqOut(:,:nAtom,1),dim=1)*coord(ii,:nAtom))
          end do
          write(*,*)
          write(1001+icart,"(5E20.12)")omega(iOmega)*Hartree__eV, &
              & polarisability(:,iCart), sqrt(sum((dqOut(:,:,1)*dqOut(:,:,1))))
          call flush(1001+icart)
        end if
        
      end do
      
    end do
    
    if (ioproc) then
      do ii = 1, 3
        close(1001+icart)
      end do
    end if
    call destroy(dpotential)
    DEALLOCATE_(dham)
    DEALLOCATE_(nFilled)
    DEALLOCATE_(nEmpty)
    
    call myclock_all%stop()
    call myclock_all%print_times(msg="linear response total", fp=stdout)

  end subroutine perturbDyn
  
!  subroutine sternheimerSingletDense(mympi, allproc, grpproc, gridAtom, &
!      & groupKS, nGroups, groupsize, desc, filling, nEl, SSqrReal, eigvals, &
!      & eigvecs, ham, over, orb, nAtom, species, speciesnames, tol, iNeighbor, &
!      & nNeighbor, iAtomStart, iPair, img2CentCell, coord, nSCCIter, &
!      & nMixElements, omega)
!    type(mpifx_comm), intent(in)   :: mympi
!    type(blacsgrid), intent(in)    :: allproc
!    type(blacsgrid), intent(in)    :: grpproc
!    type(blacsgrid), intent(in)    :: gridAtom
!    integer, intent(in)            :: groupKS(:,:)
!    integer, intent(in)            :: nGroups, groupsize
!    integer, intent(in)            :: desc(DLEN_)
!    real(dp), intent(in)           :: filling(:,:,:)
!    real(dp), intent(in)           :: nEl(:)
!    real(dp), intent(inout)        :: SSqrReal(:,:)
!    real(dp), intent(in)           :: eigvals(:,:,:)
!    real(dp), intent(in)           :: eigvecs(:,:,:)
!    real(dp), intent(in)           :: ham(:,:), over(:)
!    type(TOrbitals), pointer       :: orb
!    integer, intent(in)            :: nAtom, species(:)
!    character(mc), intent(in)      :: speciesnames(:)
!    real(dp), intent(in)           :: tol
!    integer, intent(in)            :: iNeighbor(0:,:), nNeighbor(:)
!    integer, intent(in)            :: iAtomStart(:), iPair(0:,:)
!    integer, intent(in)            :: img2CentCell(:)
!    real(dp), intent(in)           :: coord(:,:)
!    integer, intent(in)            :: nSCCIter, nMixElements
!    real(dp), intent(in), optional :: omega(:)
!    
!    type(timer) :: myclock_all
!    integer     :: iS, iKS, iAt, iCart, iSCC
!    integer     :: nSpin, nKS, nOrbs, nSparse
!    integer, allocatable  :: nFilled(:)
!    
!    integer :: ii, jj, iGlob, jGlob
!    integer :: iSCCIter
!    logical :: tStopSCC
!    
!    ! derivatives of hamiltonian, eigenvalues and eigenvectors
!    real(dp), allocatable :: dham(:,:), deigvals(:,:,:), dci(:,:,:)
!    real(dp) :: drho(size(over),size(ham, dim=2))
!    real(dp) :: dqIn(orb%mOrb,nAtom,size(ham, dim=2))
!    real(dp) :: dqOut(orb%mOrb,nAtom,size(ham, dim=2))
!    real(dp) :: dqInpRed(nMixElements), dqOutRed(nMixElements)
!    type(TPotentials)     :: dpotential
!    
!    real(dp), allocatable :: ciTmp(:,:,:)
!    
!    real(dp), allocatable :: Pc(:,:,:) ! projector onto empty states
!    
!    type(OMixer), pointer :: pChrgMixer    !* Charge mixer
!    logical :: tSCC
!    
!    nSpin = size(ham, dim=2)
!    nKS = size(groupKS, dim=2)
!    nSparse = size(ham,dim=1)
!    
!    tSCC = (nSCCIter > 1)
!    
!    ASSERT(size(eigvecs,dim=3) == nSpin)
!    ASSERT(all(shape(coord) == [3,nAtom]))
!    
!    call create(dpotential,orb,nAtom,nSpin)
!    
!    call myclock_all%start(1)
!    
!    ALLOCATE_(nFilled,(nSpin))
!    do iS = 1, nSpin
!      nFilled(iS) = floor(nEl(iS)/real(3-nSpin,dp))
!    end do
!    write(*,*)'Filled',nFilled
!    
!    ALLOCATE_(dham,(nSparse,nSpin))
!    ALLOCATE_(Pc,(size(eigvecs,dim=1),size(eigvecs,dim=2),nSpin))
!    ALLOCATE_(ciTmp,(size(eigvecs,dim=1),size(eigvecs,dim=2),nSpin))
!    ALLOCATE_(dci,(size(eigvecs,dim=1),size(eigvecs,dim=2),nSpin))
!    
!    Pc = 0.0_dp
!    
!    call projEmptyReal(grpproc, groupKS, desc, SSqrReal, nFilled, eigvecs, &
!        & over, orb, species, iNeighbor, nNeighbor, iAtomStart, iPair, &
!        & img2CentCell, Pc)
!    
!    !! test that this is a projection onto the empty states
!    !ciTmp = 0.0_dp
!    !do iKS = 1, nKS
!    !  iS = groupKS(2, iKS)
!    !  call pblasfx_pgemm(Pc(:,:,iS), desc, eigvecs(:,:,iS), desc, &
!    !      & ciTmp(:,:,iS), desc, transa="T")
!    !
!    !  call writeEigvecs(200, 101, mympi,allProc,grpproc, groupKS, &
!    !      & nGroups,groupsize, over, desc, iNeighbor, &
!    !      & nNeighbor, iAtomStart, iPair, img2CentCell, dci(:,:,iS), &
!    !      & SSqrReal, ciTmp, nAtom, orb, species, speciesNames )
!    !end do
!    
!    do iCart = 1, 3 ! polarization direction
!      
!      dci = 1.0_dp
!      dqIn = 0.0_dp
!      
!      ! derivative of E.x
!      dpotential%extAtom = 0.0_dp
!      do iAt = 1, nAtom
!        dpotential%extAtom(iAt,1) = coord(iCart,iAt)
!      end do
!      dpotential%extBlock = 0.0_dp
!      call total_shift(dpotential%extBlock, dpotential%extAtom, orb, species)
!      
!      if (tSCC) then
!        call reset(pChrgMixer, nMixElements)
!        dqInpRed = 0.0_dp
!      end if
!      
!      iSCCIter = 1
!      tStopSCC = .false.
!      lpSCC: do while (iSCCiter <= nSCCIter)
!        
!        dpotential%intAtom = 0.0_dp
!        dpotential%intBlock = 0.0_dp
!        dpotential%intShell = 0.0_dp
!        
!        call updateCharges_SCC(allproc, gridAtom, dqIn, orb, species, &
!            & iNeighbor, img2CentCell)
!        call getShiftPerAtom(dpotential%intAtom)
!        call getShiftPerL(dpotential%intShell)
!        
!        call total_shift(dpotential%intShell,dpotential%intAtom, orb,species)
!        call total_shift(dpotential%intBlock,dpotential%intShell, orb,species)
!        
!        dpotential%intBlock = dpotential%intBlock + dpotential%extBlock
!        
!        dham = 0.0_dp
!        call add_shift(dham, over, nNeighbor, iNeighbor, species, orb, iPair, &
!            & nAtom, img2CentCell, dpotential%intBlock)
!        
!        if (nSpin == 2) then
!          dham(:,:) = 2.0_dp *  dham(:,:)
!        end if
!        call qm2ud(dham)
!        deigvals = 0.0_dp
!        drho = 0.0_dp
!        
!        
!      end do lpSCC
!      
!    end do
!    
!    DEALLOCATE_(Pc)
!    
!    call destroy(dpotential)
!    
!    DEALLOCATE_(dci)
!    DEALLOCATE_(ciTmp)
!    DEALLOCATE_(dham)
!
!    call myclock_all%stop()
!    call myclock_all%print_times(msg="linear response total", fp=stdout)
!
!    DEALLOCATE_(nFilled)
!
!  end subroutine sternheimerSingletDense
  
end module LRderivs
