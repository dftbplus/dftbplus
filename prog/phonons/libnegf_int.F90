! **************************************************************************
! *                INTERFACE of LIBNEGF for DFTB+
! *
! *  call negf_init(structure,transpar,greendens,tundos,mpicomm,initinfo)
! *  
! *  call negf_density(miter,spin,nkpoint,HH,SS,DensMat,EnMat)
! *
! *  call negf_current(HH,SS,spin,kpoint,wght,tunn,ledos,currents)
! *
! *  call negf_destroy()
! *  ----------------------------------------------------------------------
! *
! *  ---------------------------------------------------------------------- 
! * Known Bugs & Problems:
! * 
! *
! *
! * ---------------------------------------------------------------------- 
module libnegf_int
  
  use ln_precision
  use ln_constants
  use ln_allocation
  use libnegf_vars, only : TGDFTBstructure, TGDFTBTunDos, TTransPar
  use ln_structure
  use mpi_globals
  use mat_def 
  use lib_param, only : Tnegf, set_defaults, set_phph
  use libnegf
  use ln_extract
  use mat_conv 
  use commonTypes, only : TOrbitals
  use libmpifx_module 
  use sparsekit_drv, only : dns2csr, csr2dns, nzdrop, zprint_csrdns
  use initprogram, only : TempMin, TempMax, TempStep
  use phph

  implicit none
  private

  Type(Tnegf), target, save, public :: negf

  public :: negf_init, negf_destroy
 
  ! This initializes the partitioned structure 
  public :: negf_init_str

  ! initialize phonon-phonon interactions (call after negf_init_str)
  public :: negf_init_phph

  ! INTERFACE SUBROUTINE TO DRIVE PHONON CALCULATIONS:  
  public :: calc_phonon_current


  ! direct calls to compute phonon current
  public :: negf_phonon_current 

  ! direct calls to compute thermal conductance
  public :: thermal_conductance

  ! wrapped functions passing dftb matrices. Needed for parallel
  !public :: calc_current
  !public :: calcdensity_green
  !public :: calcEdensity_green
  !public :: calcPDOS_green
  !public :: local_currents
  
  ! interface csr matrices. The pattern structure of csrHam
  ! is defined by negf_init_csr 
  public :: negf_init_csr
  type(z_CSR) :: csrHam

  type(z_CSR) :: csrOver ! need to be removed
 

  contains

!------------------------------------------------------------------------------
! Init gDFTB environment and variables
!------------------------------------------------------------------------------
  subroutine negf_init(structure, transpar, tundos, mpicomm, initinfo)
    Type(TGDFTBstructure), intent(IN) :: structure
    Type(TTranspar), intent(IN) :: transpar
    Type(TGDFTBTunDos), intent(IN) :: tundos
    Type(mpifx_comm), intent(in) :: mpicomm
    logical, intent(OUT) :: initinfo
    

    ! local variables
    integer :: i,error, l, ncont, nc_vec(1), j, nldos
    integer, dimension(:), allocatable :: sizes  
    ! string needed to hold processor name
    character(:), allocatable :: hostname
    type(lnParams) :: parms
   
    error = 0 
    initinfo = .true.       

    call negf_mpi_init(mpicomm)
     
    if (transpar%defined) then
      ncont = transpar%ncont 
    else
      ncont = 0
    endif
      
    ! ------------------------------------------------------------------------------
    ! check that GS and contacts are created  (one day this will be automatic)
    ! ------------------------------------------------------------------------------
    if(transpar%defined .and. ncont.gt.0) then
                                
        open(199,file='GS/test',IOSTAT=error)
        if (error.ne.0) then 
          call mpifx_get_processor_name(hostname)
          write(*,*) "ERROR: please create a directory called GS on "//trim(hostname)
          initinfo = .false.; return  
        end if   
        close(199)
        open(199,file='contacts/test',IOSTAT=error)
        if (error.ne.0) then 
          call mpifx_get_processor_name(hostname)
          write(*,*) "ERROR: please create a directory called contacts "//trim(hostname)
          initinfo = .false.; return  
        end if         
        close(199)
    
    end if
    ! ------------------------------------------------------------------------------
    !! Set defaults and fill up the parameter structure with them
    call init_negf(negf)
    call get_params(negf, parms)

    ! ------------------------------------------------------------------------------
    ! This parameter is used to set the averall drop threshold in libnegf
    ! It affects especially transmission that is not accurate more than 
    ! this value.
    call set_drop(1.d-20)

 

    ! ------------------------------------------------------------------------------
    !                        SETTING CONTACT TEMPERATURES 
    ! ------------------------------------------------------------------------------
    ! If no contacts => no transport,
    if (transpar%defined) then
      
      do i = 1, ncont
        if (transpar%contacts(i)%kbT .ge. 0.0_dp) then
          parms%kbT(i) = transpar%contacts(i)%kbT
        else
          parms%kbT(i) = structure%tempElec
        end if
      enddo

      ! set parameters for wide band approximations
      do i=1, ncont
         parms%FictCont(i) = transpar%contacts(i)%wideBand
         parms%contact_DOS(i) = transpar%contacts(i)%wideBandDOS
      enddo

    end if

    ! ------------------------------------------------------------------------------
    !                    SETTING TRANSMISSION PARAMETERS
    ! ------------------------------------------------------------------------------

    if (tundos%defined) then

      parms%verbose = tundos%verbose
      parms%delta = tundos%delta      ! delta for G.F.
      parms%dos_delta = tundos%broadeningDelta

      l = size(tundos%ni)    
      parms%ni(1:l) = tundos%ni(1:l)
      parms%nf(1:l) = tundos%nf(1:l)

      parms%Emin =  tundos%Emin
      parms%Emax =  tundos%Emax
      parms%Estep = tundos%Estep
      !checks if the energy interval is appropriate
      !if (id0.and.mumin.lt.mumax) then
             
      !if (negf%Emin.gt.mumin-10*negf%kbT .or. negf%Emax.lt.mumin-10*negf%kbT) then
      !  write(*,*) 'WARNING: the interval Emin..Emax is smaller than the bias window'
      !  write(*,*) 'Emin=',negf%emin*negf%eneconv, 'Emax=',negf%emax*negf%eneconv
      !  write(*,*) 'kT=',negf%kbT*negf%eneconv    
      !  write(*,*) 'Suggested interval:', &
      !        (mumin-10*negf%kbT)*negf%eneconv,(mumax+10*negf%kbT)*negf%eneconv
      ! endif
      !endif

    endif
    
    ! Energy conversion only affects output units. 
    ! The library writes energies as (E * negf%eneconv) 
    parms%eneconv = 1.d0

    parms%isSid = .true.

    if (allocated(sizes)) call log_deallocate(sizes)

    call set_params(negf,parms)

    if (tundos%defined) then
      ! setting of intervals and indeces for projected DOS
      nldos = size(tundos%dosOrbitals)
      call init_ldos(negf, nldos)
      do i = 1, nldos
         call set_ldos_indexes(negf, i, tundos%dosOrbitals(i)%data)   
      end do 
    end if

  end subroutine negf_init
  
  !------------------------------------------------------------------------------
  subroutine negf_destroy()

    write(*,'(A)') 'Release Negf Memory:'
    call destruct(csrHam)
    call destruct(csrOver)
    call destroy_negf(negf)
    call writePeakInfo(6)    
    call writeMemInfo(6)

  end subroutine negf_destroy

  !------------------------------------------------------------------------------
  subroutine negf_init_str(structure, transpar, iNeigh, nNeigh, img2CentCell)
    Type(TTranspar), intent(IN) :: transpar
    Type(TGDFTBStructure), intent(IN) :: structure
    Integer, intent(in) :: nNeigh(:)
    Integer, intent(in) :: img2CentCell(:)
    Integer, intent(in) :: iNeigh(0:,:)
    
    Integer, allocatable :: PL_end(:), cont_end(:), surf_end(:), cblk(:), ind(:)
    Integer, allocatable :: atomst(:), plcont(:)   
    integer, allocatable :: minv(:,:)
    Integer :: natoms, ncont, nbl, iatm1, iatm2, iatc1, iatc2
    integer :: i, m, i1, j1

    iatm1 = transpar%idxdevice(1)
    iatm2 = transpar%idxdevice(2)

    ncont = transpar%ncont
    nbl = transpar%nPLs
    if (nbl.eq.0) STOP 'Internal ERROR: nbl = 0 ?!'
    natoms = structure%nAtom

    call log_allocate(PL_end,nbl)
    call log_allocate(atomst,nbl+1)
    call log_allocate(plcont,nbl)
    call log_allocate(cblk,ncont)
    call log_allocate(cont_end,ncont)
    call log_allocate(surf_end,ncont)
    call log_allocate(ind,natoms+1)
    call log_allocate(minv,nbl,ncont)

    ind(1:natoms+1)=structure%iatomstart(1:natoms+1) - 1

    do i = 1, ncont
       cont_end(i) = ind(transpar%contacts(i)%idxrange(2)+1)
       surf_end(i) = ind(transpar%contacts(i)%idxrange(1))
    enddo
     
    if (transpar%defined) then
      do i = 1, nbl-1
        PL_end(i) = ind(transpar%PL(i+1))
      enddo
      atomst(1:nbl) = transpar%PL(1:nbl)
      PL_end(nbl) = ind(transpar%idxdevice(2)+1)
      atomst(nbl+1) = iatm2 + 1
    endif

    ! For every contact finds the min-max atom indeces among
    ! the atoms in the central region interacting with contact 
    if (transpar%defined .and. ncont.gt.0) then

       minv = 0

       do m = 1, transpar%nPLs
          ! Loop over all PL atoms      
          do i = atomst(m), atomst(m+1)-1 

             ! Loop over all contacts
             do j1 = 1, ncont

                iatc1 = transpar%contacts(j1)%idxrange(1)
                iatc2 = transpar%contacts(j1)%idxrange(2)

                i1 = minval(img2CentCell(iNeigh(1:nNeigh(i),i)), &
                     mask = (img2CentCell(iNeigh(1:nNeigh(i),i)).ge.iatc1 .and. & 
                     img2CentCell(iNeigh(1:nNeigh(i),i)).le.iatc2) )
                if (i1.ge.iatc1 .and. i1.le.iatc2) then
                   minv(m,j1) = j1
                endif

             end do
          end do  
       end do


       do j1 = 1, ncont      
          if (count(minv(:,j1).eq.j1).gt.1) then
             write(*,*) 'ERROR: contact',j1,'interacts with more than one PL'
             write(*,*) minv(:,j1)
             stop
          end if
          do m = 1, transpar%nPLs
             if (minv(m,j1).eq.j1) cblk(j1) = m
          end do
       end do

       if (negf%verbose.gt.50) then
          write(*,*) ' Structure info:'
          write(*,*) ' Nconts:', ncont
          write(*,*) ' Number of PLs:',nbl 
          write(*,*) ' Interacting PLs:',cblk(1:ncont)
          write(*,*)
       endif

    end if     

    call kill_Tstruct(negf%str)
 
    call create_Tstruct(ncont, nbl, PL_end, cont_end, surf_end, cblk, negf%str)

    call log_deallocate(PL_end)
    call log_deallocate(plcont)
    call log_deallocate(atomst)
    call log_deallocate(cblk)
    call log_deallocate(cont_end)
    call log_deallocate(surf_end)
    call log_deallocate(ind)
    call log_deallocate(minv)
  
  end subroutine negf_init_str

   
  subroutine negf_init_phph(negf, order)
     type(Tnegf) :: negf
     integer, intent(in) :: order

     select case (order)
     case(3, 34)
       print*,'Init cubic phonon-phonon interactions'
       call set_phph(negf, 3, 'cubic.dat')
     case(4)
       print*,'Init quartic phonon-phonon interactions'
       call set_phph(negf, 4, 'quartic.dat')
     end select 

  end subroutine negf_init_phph

  !------------------------------------------------------------------------------
  ! INTERFACE subroutine to call phonon current computation
  !------------------------------------------------------------------------------    
  subroutine calc_phonon_current(mpicomm, DynMat, tunnTot, ldosTot, currTot, &
                        & writeTunn, writeLDOS)  
    type(mpifx_comm), intent(in) :: mpicomm
    real(dp), intent(in) :: DynMat(:,:)
    real(dp), allocatable, intent(inout) :: tunnTot(:,:), ldosTot(:,:)
    real(dp), allocatable, intent(inout) :: currTot(:)
    logical, intent(in) :: writeLDOS
    logical, intent(in) :: writeTunn
 
    ! locals 
    real(dp), pointer    :: tunnMat(:,:)=>null()
    real(dp), pointer    :: ldosMat(:,:)=>null()
    real(dp), pointer    :: currVec(:)=>null()    
    integer :: ii, jj, iK, nK, err, nnz, ntemp
    real(dp), allocatable :: kPoints(:,:), kWeights(:)
    type(z_DNS) :: zDynMat
    real(dp), allocatable :: tunnSKRes(:,:,:), ldosSKRes(:,:,:)
    real(dp) :: cutoff,TT1,emin,emax,estep, kappa
    type(unit) :: HessianUnits, HeatCurrUnits, HeatCondUnits

    nK = 1
    call log_allocate(kPoints, 3, nK)
    call log_allocate(kWeights, nK)
    kPoints = 0.0_dp
    kWeights = 1.0_dp
    cutoff = 1d-30
    HessianUnits%name = "H"
    HeatCurrUnits%name = "W"
    HeatCondUnits%name = "W/K"

    do iK = 1, nK

      call create(zDynMat, size(DynMat,1), size(DynMat,2) )
      zDynMat%val = DynMat

      nnz = nzdrop(zDynMat,cutoff)

      call create(csrHam, zDynMat%nrow, zDynMat%ncol, nnz)

      call dns2csr(zDynMat, csrHam)
   
      call destroy(zDynMat)


      call negf_phonon_current(csrHam, iK, kWeights(iK), &
            tunnMat, ldosMat, currVec)

      if(.not.allocated(currTot)) then
         allocate(currTot(size(currVec)), stat=err)
         if (err/=0) STOP 'Allocation error (currTot)'
         currTot = 0.0_dp
       endif
       currTot = currTot + currVec

       call add_partial_results(tunnMat, tunnTot, tunnSKRes, iK, nK) 
       
       call add_partial_results(ldosMat, ldosTot, ldosSKRes, iK, nK) 

       ! MPI Reduce K dependent stuff
       if (associated(ldosMat).and.(nK.gt.1)) then
         call mpifx_allreduceip(mpicomm, ldosSKRes(:,:,iK), MPI_SUM)
       end if
       if (associated(tunnMat).and.(nK.gt.1)) then
         call mpifx_allreduceip(mpicomm, tunnSKRes(:,:,iK), MPI_SUM)
       end if
      
    end do

    if (associated(negf%tunn_mat)) call log_deallocatep(negf%tunn_mat)
    if (associated(negf%ldos_mat)) call log_deallocatep(negf%ldos_mat) 
    if (associated(negf%currents)) call log_deallocatep(negf%currents)

    call mpifx_allreduceip(mpicomm, currTot, MPI_SUM)

    ! converts from internal atomic units into W
    currTot = currTot * convertHeatCurrent(HessianUnits, HeatCurrUnits)
   
    if (id0) then 
      do ii= 1, size(currTot) 
        write(*,'(1x,a,i3,i3,a,ES14.5,a,a)') &
             & ' contacts: ',negf%ni(ii),negf%nf(ii), &
             & ' current: ', currTot(ii),' ',HeatCurrUnits%name
      enddo
    endif

    if (allocated(tunnTot)) then
      call mpifx_allreduceip(mpicomm, tunnTot, MPI_SUM)
      ! Write Total tunneling on a separate file (optional)
      if (id0 .and. writeTunn) then
        call write_file(negf, tunnTot, tunnSKRes, 'tunneling', kpoints, kWeights)
        if (allocated(tunnSKRes)) deallocate(tunnSKRes)
      endif 
    else
      allocate(tunnTot(0,0))
    endif

    if (allocated(ldosTot)) then
      call mpifx_allreduceip(mpicomm, ldosTot, MPI_SUM)
      ! Multiply density by 2w 
      do ii = 1, size(ldosTot,1)
        ldosTot(ii,:) = ldosTot(ii,:) * 2.d0*(negf%emin + negf%estep*(ii-1))
      end do
      ! Write Total localDOS on a separate file (optional)
      if (id0 .and. writeLDOS) then
        call write_file(negf, ldosTot, ldosSKRes, 'localdos', kpoints, kWeights)
        if (allocated(ldosSKRes)) deallocate(ldosSKRes)
      end if
    else
      allocate(ldosTot(0,0))
    endif
    
    call log_deallocate(kPoints)
    call log_deallocate(kWeights)

    if (allocated(tunnTot)) then

      open(unit=50,file='conductance.dat',action='write')
      emin = negf%Emin*negf%eneconv
      emax = negf%Emax*negf%eneconv
      estep = negf%Estep*negf%eneconv
 
      ntemp=nint((TempMax-TempMin)/TempStep)
 
      do ii = 1, size(TunnTot,2)
        write(50,*) '# T [K]', 'Thermal Conductance [W/K]'  
        do jj = 1, ntemp
          TT1 = TempMin + TempStep*(jj-1)
          kappa = thermal_conductance(TunnTot(:,ii),TT1,emin,emax,estep) 
          kappa = kappa * convertHeatConductance(HessianUnits,HeatCondUnits)
          write(50,*)  TT1/kb, kappa  
        end do
 
      end do
 
    else
      allocate(tunnTot(0,0))
    end if

  end subroutine calc_phonon_current

  
  
  subroutine negf_phonon_current(HH,qpoint,wght,tunn,ledos,currents)

    type(z_CSR), intent(in) :: HH
    integer, intent(in) :: qpoint        ! kp index 
    real(dp), intent(in) :: wght      ! kp weight 
    real(dp), dimension(:,:), pointer :: tunn
    real(dp), dimension(:,:), pointer :: ledos
    real(dp), dimension(:), pointer :: currents 

    integer :: i
    
    negf%kpoint = qpoint
    negf%wght = wght
    
    if (associated(negf%currents)) then
      call log_deallocatep(negf%currents)
      negf%currents=> null()
    endif
    currents=>null() 
    if (associated(negf%tunn_mat)) then
      call log_deallocatep(negf%tunn_mat)
      negf%tunn_mat=> null()
    endif
    tunn=>null() 
    if (associated(negf%ldos_mat)) then
      call log_deallocatep(negf%ldos_mat)
      negf%ldos_mat=> null()
    endif
    ledos=>null()
    
    if (id0.and.negf%verbose.gt.30) then
      write(*,*)
      write(*,'(73("="))')
      write(*,*) '                   COMPUTING TUNNELING AND DOS          '
      write(*,'(73("="))') 
      write(*,*)
    endif

    call pass_HS(negf,HH)

    !call printH(negf%H)

    call compute_phonon_current(negf)

    if (associated(negf%currents)) then
      currents => negf%currents
    else
      stop 'INTERNAL ERROR: negf%currents not associated'      
    endif

    if (associated(negf%tunn_mat)) then
      tunn => negf%tunn_mat
    endif

    if (associated(negf%ldos_mat)) then
       ledos => negf%ldos_mat
    endif

  end subroutine negf_phonon_current

  subroutine printH(H)
    type(z_CSR) :: H

    type(z_DNS) :: tmp
    integer :: ii, jj
    real(dp) :: maxv

    call create(tmp, H%nrow, H%ncol)

    call csr2dns(H,tmp)

    maxv = maxval(abs(tmp%val)) 
    print*,'maxval= ',maxv
    print*,'Normalized Dynamical Matrix:'

    do ii = 1, tmp%nrow, 3
       do jj = 1, tmp%ncol, 3
          write(*,'(F8.4)',advance='no') real(tmp%val(ii,jj))/maxv
       end do 
       write(*,*)
    end do

    call destroy(tmp)
 
  end subroutine printH

  !----------------------------------------------------------------------------
  subroutine add_partial_results(pMat, pTot, pSKRes, iK, nK)
    real(dp), pointer :: pMat(:,:)
    real(dp), allocatable :: pTot(:,:)
    real(dp), allocatable :: pSKRes(:,:,:)
    integer, intent(in) :: iK, nK
      
    integer :: err 

    if (associated(pMat)) then
      if(.not.allocated(pTot)) then 
        allocate(pTot(size(pMat,1), size(pMat,2)), stat=err)
        if (err/=0) STOP 'Allocation error (tunnTot)'
        pTot = 0.0_dp
      endif
      pTot = pTot + pMat
   
      if (nK.gt.1) then
        if(.not.allocated(pSKRes)) then 
          allocate(pSKRes(size(pMat,1), size(pMat,2), nK), stat=err)
          if (err/=0) STOP 'Allocation error (tunnSKRes)'
          pSKRes = 0.0_dp
        endif
        pSKRes(:,:,iK) = pMat(:,:)
      end if
    endif

  end subroutine add_partial_results 
  !----------------------------------------------------------------------------
  
  subroutine negf_init_csr(iAtomStart, iNeighbor, nNeighbor, img2CentCell, orb)
    integer, intent(in) :: iAtomStart(:)    
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    type(TOrbitals), intent(in) :: orb

    if (allocated(csrHam%nzval)) call destroy(csrHam)
    call init(csrHam, iAtomStart, iNeighbor, nNeighbor, &
        &img2CentCell, orb) 
    if (allocated(csrOver%nzval)) call destroy(csrOver)    
    call init(csrOver, csrHam)
     
  end subroutine negf_init_csr


  subroutine write_file(negf, pTot, pSKRes, filename, kpoints, kWeights)
    type(TNegf) :: negf
    real(dp), intent(in) :: pTot(:,:)
    real(dp), intent(in) :: pSKRes(:,:,:)
    character(*), intent(in) :: filename
    real(dp), intent(in) :: kPoints(:,:)
    real(dp), intent(in) :: kWeights(:)

    integer :: ii, jj, iK, nK

    nK = size(kPoints,2)

    open(65000,file=trim(filename)//'.dat')
    if (trim(filename).eq.'tunneling') then
      write(65000,*)  '# Energy [H]', '  Tunneling' 
    else  
      write(65000,*)  '# Energy [H]', '  LDOS'
    endif 
    do ii=1,size(pTot,1)
      write(65000,'(es20.8)',ADVANCE='NO') (negf%Emin+(ii-1)*negf%Estep)*negf%eneconv
      do jj=1,size(pTot,2)
        write(65000,'(es20.8)',ADVANCE='NO') pTot(ii,jj)
      enddo
      write(65000,*)
    enddo
    close(65000)

    if (nK.gt.1) then
      open(65000,file=trim(filename)//'_kpoints.dat')
      write(65000,*)  '# NKpoints = ', nK
      write(65000,*)  '# Energy [eV], <k1 k2 k3 weight> '
      write(65000,'(A1)', ADVANCE='NO') '# '
      do iK = 1,nK
        write(65000,'(es15.5, es15.5, es15.5, es15.5)', ADVANCE='NO') kpoints(:,iK),&
                                                                      & kWeights(iK) 
      end do
      write(65000,*)
      do ii=1,size(pSKRes(:,:,1),1)
        write(65000,'(f20.8)',ADVANCE='NO') (negf%Emin+(ii-1)*negf%Estep)*negf%eneconv
        do jj=1,size(pSKRes(:,:,1),2)
          do iK = 1,nK
            write(65000,'(es20.8)',ADVANCE='NO') pSKRes(ii,jj, iK)
          enddo
          write(65000,*)
        enddo
      enddo
      close(65000)
    end if

  end subroutine write_file
  !----------------------------------------------------------------------------
  ! DEBUG routine dumping H and S on file in Matlab format
  !----------------------------------------------------------------------------
  subroutine negf_dumpHS(HH,SS)
    type(z_CSR), intent(in) :: HH, SS

    write(*,*) 'Dumping H and S on files...'    
    open(1121, file='HH.dat')
    write(1121,*) '% Size =',HH%nrow, HH%ncol
    write(1121,*) '% Nonzeros =',HH%nnz
    write(1121,*) '% '
    write(1121,*) 'zzz = ['
    call printcsr(1121,HH)
    write(1121,*) ']'
    close(1121) 
    open(1121, file='SS.dat')
    write(1121,*) '% Size =',SS%nrow, SS%ncol
    write(1121,*) '% Nonzeros =',SS%nnz
    write(1121,*) '% '
    write(1121,*) 'zzz = ['
    call printcsr(1121,SS)
    write(1121,*) ']'
    close(1121) 
  end subroutine negf_dumpHS
  
    
end module libnegf_int
       
