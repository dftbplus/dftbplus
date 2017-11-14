  ! **************************************************************************
! *                INTERFACE of LIBNEGF for DFTB+
! *
! *  call negf_init(transpar, greendens, tundos, mpicomm, itempElec, initinfo)
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
module negf_int
  
  use accuracy
  use libnegf_vars  !, only : TNEGFstructure, TGDFTBTunDos, &
                    !       & TNEGFGreenDensInfo, TTransPar, Telph
  use libnegf
  use mat_conv 
  use commonTypes, only : TOrbitals
  use libmpifx_module

  use FormatOut                                                          
  
  implicit none
  private

  Type(Tnegf), target, public :: negf

  public :: negf_init      ! general library initializations
  public :: negf_init_str  ! passing structure parameters
  public :: negf_init_dephasing
  public :: negf_init_elph                                               
  public :: negf_destroy
  
  ! wrapped functions passing dftb matrices. Needed for parallel
  public :: calcdensity_green
  public :: calcEdensity_green
  public :: calcPDOS_green
  public :: calc_current
  public :: local_currents
  
  ! interface csr matrices. The pattering must be predefined 
  ! using negf_init_csr 
  public :: negf_init_csr
  type(z_CSR) :: csrHam, csrOver

  ! non wrapped direct calls  
  private :: negf_density, negf_current, negf_ldos

  contains

!------------------------------------------------------------------------------
! Init gDFTB environment and variables
!------------------------------------------------------------------------------
  subroutine negf_init(transpar, greendens, tundos, mpicomm,&
        & tempElec, initinfo)
    
    Type(TTranspar), intent(IN) :: transpar
    Type(TNEGFGreenDensInfo), intent(IN) :: greendens
    Type(TNEGFTunDos), intent(IN) :: tundos
    Type(mpifx_comm), intent(in) :: mpicomm
    real(dp), intent(in) :: tempElec 
    logical, intent(OUT) :: initinfo
    
    ! local variables
    real(dp), dimension(:), allocatable :: pot, eFermi
    integer :: i,error, l, ncont, nc_vec(1), j, nldos
    integer, dimension(:), allocatable :: sizes  
    ! string needed to hold processor name
    character(:), allocatable :: hostname
    type(lnParams) :: params

    
    
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
    call get_params(negf, params)
    
    ! ------------------------------------------------------------------------------
    !This must be different for different initialisations, to be separated
    !Higher between transport and greendens is taken, temporary
    if (tundos%defined .and. greendens%defined) then                                          
       if (tundos%verbose.gt.greendens%verbose) then 
          params%verbose = tundos%verbose
       else
          params%verbose = greendens%verbose
       endif
    else
       if (tundos%defined) params%verbose = tundos%verbose
       if (greendens%defined) params%verbose = greendens%verbose
    end if
    ! ------------------------------------------------------------------------------
    ! This parameter is used to set the averall drop threshold in libnegf
    ! It affects especially transmission that is not accurate more than 
    ! this value.
    call set_drop(1.d-20)

 
    ! ------------------------------------------------------------------------------
    ! Assign spin degenracy and check consistency between different input blocks
    if (tundos%defined .and. greendens%defined) then
      if (tundos%gSpin .ne. greendens%gSpin) then
        write(*,*) "ERROR: spin degeneracy is not consistent between different input blocks"
        initinfo = .false.; return 
      else
        params%g_spin = real(tundos%gSpin) ! Spin degeneracy
      end if
    else if (tundos%defined) then
        params%g_spin = real(tundos%gSpin) ! Spin degeneracy
    else if (greendens%defined) then
        params%g_spin = real(greendens%gSpin) ! Spin degeneracy
    end if
  
    ! Setting contact temperatures
    do i = 1, ncont
      if (greendens%defined) then
        params%kbT_dm(i) = greendens%kbT(i)
      else        
        params%kbT_dm(i) = tempElec
      end if  
      ! Make sure low temperatures (< 10K) converted to 0.
      ! This avoid numerical problems with contour integration
      if (params%kbT_dm(i) < 3.0e-5_dp) then
        params%kbT_dm(i) = 0.0_dp
      end if
      if (tundos%defined) then
        params%kbT_t(i) = tundos%kbT(i)
      else
        params%kbT_t(i) = tempElec
      end if    
    end do

    ! ------------------------------------------------------------------------------
    !            SETTING ELECTROCHEMICAL POTENTIALS INCLUDING BUILT-IN 
    ! ------------------------------------------------------------------------------
    ! Fermi level is given by the contacts. If no contacts => no transport,
    ! Then Fermi is defined by the Green solver
    if (transpar%defined) then
      pot = transpar%contacts(1:ncont)%potential
      eFermi = transpar%contacts(1:ncont)%eFermi(1) 
      do i = 1,ncont
        ! Built-in potential to equilibrate Fermi levels
        pot(i) = pot(i) + eFermi(i) - minval(eFermi(1:ncont))

        ! set parameters for wide band approximations
        params%FictCont(i) = transpar%contacts(i)%wideBand
        params%contact_DOS(i) = transpar%contacts(i)%wideBandDOS
      
        if (id0) then
          write(*,*) '(negf_init) CONTACT INFO #',i
            if (params%FictCont(i)) then
              write(*,*) 'FICTICIOUS CONTACT '
              write(*,*) 'DOS: ', params%contact_DOS(i)
            end if
          write(*,*) 'Temperature (DM): ', params%kbT_dm(i)
          write(*,*) 'Temperature (Current): ', params%kbT_t(i)
          write(*,*) 'Potential (with built-in): ', pot(i)
          write(*,*) 'eFermi: ', eFermi(i) 
          write(*,*) 
        endif 

      enddo

      ! Define electrochemical potentials
      params%mu(1:ncont) = eFermi(1:ncont) - pot(1:ncont)
      if (id0) write(*,*) 'Electro-chemical potentials: ', params%mu(1:ncont)
      deallocate(pot)

    else !transpar not defined 
      params%mu(1) = greendens%oneFermi(1) 
    end if


    ! ------------------------------------------------------------------------------
    !                  SETTING COUNTOUR INTEGRATION PARAMETERS
    ! ------------------------------------------------------------------------------
    if (greendens%defined) then
      params%Ec = greendens%enLow           ! lowest energy
      params%Np_n(1:2) = greendens%nP(1:2)  ! contour npoints
      params%Np_real = greendens%nP(3)      ! real axis points
      params%n_kt = greendens%nkt           ! n*kT for Fermi
 
      !Read G.F. from very first iter
      if (greendens%readSGF .and. .not.greendens%saveSGF) params%readOldSGF=0  
      !compute G.F. at every iteration 
      if (.not.greendens%readSGF .and. .not.greendens%saveSGF) params%readOldSGF=1  
      !Default Write on first iter
      if (.not.greendens%readSGF .and. greendens%saveSGF) params%readOldSGF=2  

      if(any(params%kbT_dm.gt.0) .and. greendens%nPoles.eq.0) then
         STOP 'ERROR: Number of Poles = 0 but T > 0' 
      else
         params%n_poles = greendens%nPoles
      end if
      if(all(params%kbT_dm.eq.0)) then 
        params%n_poles = 0
      end if

    end if

    ! ------------------------------------------------------------------------------
    !! Setting the delta: priority on Green Solver, if present
    !! dos_delta is used by libnegf to smoothen T(E) and DOS(E)
    !! and is currently set in tunneling
    params%dos_delta = tundos%broadeningDelta
    if (tundos%defined) then
      params%delta = tundos%delta      ! delta for G.F.
    end if
    if (greendens%defined) then
      params%delta = greendens%delta   ! delta for G.F.
    end if
    
    ! ------------------------------------------------------------------------------
    !                    SETTING TRANSMISSION PARAMETERS
    ! ------------------------------------------------------------------------------

    if (tundos%defined) then

      l = size(tundos%ni)    
      params%ni(1:l) = tundos%ni(1:l)
      params%nf(1:l) = tundos%nf(1:l)

      ! setting of intervals and indeces for projected DOS
      nldos = size(tundos%dosOrbitals)
      call init_ldos(negf, nldos)
      do i = 1, nldos
         call set_ldos_indexes(negf, i, tundos%dosOrbitals(i)%data)   
      end do 
      
      params%Emin =  tundos%Emin
      params%Emax =  tundos%Emax
      params%Estep = tundos%Estep

    endif
    
    ! Energy conversion only affects output units. 
    ! The library writes energies as (E * negf%eneconv) 
    params%eneconv = HAR 

    if (allocated(sizes)) deallocate(sizes)

    call set_params(negf,params)

    !--------------------------------------------------------------------------
    !DAR begin - negf_init - TransPar to negf
    !--------------------------------------------------------------------------
    if (transpar%defined) then
      negf%tNoGeometry = transpar%tNoGeometry
      negf%tOrthonormal = transpar%tOrthonormal
      negf%tOrthonormalDevice = transpar%tOrthonormalDevice
      negf%NumStates = transpar%NumStates
      negf%tManyBody = transpar%tManyBody
      negf%tElastic = transpar%tElastic
      negf%tZeroCurrent = transpar%tZeroCurrent
      negf%MaxIter = transpar%MaxIter
      negf%tranas%out%tWriteDOS = transpar%tWriteDOS
      negf%tWrite_ldos = transpar%tWrite_ldos
      negf%tWrite_negf_params = transpar%tWrite_negf_params
      negf%tranas%out%tDOSwithS = transpar%tDOSwithS
      allocate(negf%tranas%cont(ncont))
      negf%tranas%cont(:)%name = transpar%contacts(:)%name
      negf%tranas%cont(:)%tWriteSelfEnergy = transpar%contacts(:)%tWriteSelfEnergy
      negf%tranas%cont(:)%tReadSelfEnergy = transpar%contacts(:)%tReadSelfEnergy
      negf%tranas%cont(:)%tWriteSurfaceGF = transpar%contacts(:)%tWriteSurfaceGF
      negf%tranas%cont(:)%tReadSurfaceGF = transpar%contacts(:)%tReadSurfaceGF
    end if

    ! Defined outside transpar%defined ... HAS TO BE FIXED
    negf%tDephasingVE = transpar%tDephasingVE
    negf%tDephasingBP = transpar%tDephasingBP

    if(id0.and.negf%tWrite_negf_params) call check_negf_params
    
    if((.not.negf%tElastic).and.(.not.negf%tManyBody)) then
         write(*,*)'Current is not calculated!'
         write(*,*)'Choose "Elastic = Yes" or "ManyBody = Yes"!'
         write(*,*)'Program is terminated!'
         stop
    end if 


  end subroutine negf_init
  !-----------------------------------------------------------------------------

  subroutine negf_init_dephasing(tundos)
    Type(TNEGFTunDos), intent(in) :: tundos

    if(negf%tDephasingVE) call negf_init_elph(tundos%elph)     
    if(negf%tDephasingBP) call negf_init_bp(tundos%bp)        
    
  end subroutine negf_init_dephasing
  !-----------------------------------------------------------------------------


  subroutine negf_init_elph(elph)
  
    type(TElPh), intent(in) :: elph
    
    if (id0) write(*,*)

    if (elph%model .eq. 1) then
       if (id0) write(*,*) 'Setting local fully diagonal (FD) elastic dephasing model'
       call set_elph_dephasing(negf, elph%coupling, elph%scba_niter)
    else if (elph%model .eq. 2) then
       if (id0) write(*,*) 'Setting local block diagonal (BD) elastic dephasing model'
       call set_elph_block_dephasing(negf, elph%coupling, elph%orbsperatm, &
            elph%scba_niter)
    else if (elph%model .eq. 3) then
       if (id0) write(*,*) 'Setting overlap mask (OM) block diagonal elastic dephasing model'
       call set_elph_s_dephasing(negf, elph%coupling, elph%orbsperatm, &
            elph%scba_niter)
    else
       write(*,*) "ERROR: el-ph model is not supported"
    endif   
    
  end subroutine negf_init_elph

!!$ !----------------------------------------------------------------------------
!!$ subroutine negf_destroy_elph()
!!$
!!$    call destroy_elph_model(negf)
!!$
!!$ end subroutine negf_destroy_elph
!!$ !----------------------------------------------------------------------------
!!$

  subroutine negf_init_bp(elph)
  
    type(TElPh), intent(in) :: elph
    
    if (id0) write(*,*)
       
    if (elph%model .eq. 1) then
       if (id0) then
         write(*,*) 'Setting local fully diagonal (FD) BP dephasing model'
         write(*,*) 'coupling=',elph%coupling
       end if
       call set_bp_dephasing(negf, elph%coupling)
    else if (elph%model .eq. 2) then
       if (id0) then 
         write(*,*) 'Setting local block diagonal (BD) BP dephasing model'
         write(*,*) 'NOT IMPLEMENTED! INTERRUPTED!'
       end if
       stop
    else if (elph%model .eq. 3) then
       if (id0) then
         write(*,*) 'Setting overlap mask (OM) block diagonal BP dephasing model'
         write(*,*) 'NOT IMPLEMENTED! INTERRUPTED!'
       end if
       stop
    else
       write(*,*) "ERROR: BP model is not supported"
    endif

  end subroutine negf_init_bp
 

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
  
  !------------------------------------------------------------------------------
  subroutine negf_destroy()

    call destruct(csrHam)
    call destruct(csrOver)
    call destroy_negf(negf)
    !call writePeakInfo(6)                                                 !DAR
    !call writeMemInfo(6)                                                  !DAR

  end subroutine negf_destroy

  !------------------------------------------------------------------------------
  subroutine negf_init_str(structure, transpar, greendens, iNeigh, nNeigh, img2CentCell)
    Type(TNEGFStructure), intent(in) :: structure
    Type(TTranspar), intent(in) :: transpar
    Type(TNEGFGreenDensInfo) :: greendens
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

    if (transpar%defined) then
       nbl = transpar%nPLs
    else if (greendens%defined) then
       nbl = greendens%nPLs
    endif
    if (nbl.eq.0) STOP 'Internal ERROR: nbl = 0 ?!'
    natoms = structure%nAtom

    allocate(PL_end(nbl))
    allocate(atomst(nbl+1))
    allocate(plcont(nbl))
    allocate(cblk(ncont))
    allocate(cont_end(ncont))
    allocate(surf_end(ncont))
    allocate(ind(natoms+1))
    allocate(minv(nbl,ncont))

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
    else if (greendens%defined) then
       do i = 1, nbl-1
          PL_end(i) = ind(greendens%PL(i+1))
       enddo
       atomst(1:nbl) = greendens%PL(1:nbl)
      PL_end(nbl) = ind(natoms+1)
      atomst(nbl+1) = natoms + 1
    endif

    ! For every contact finds the min-max atom indeces among
    ! the atoms in the central region interacting with contact 
    if (transpar%defined .and. ncont.gt.0) then

      if(.not.transpar%tNoGeometry) then         !DAR

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
             write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             write(*,*) 'WARNING: contact',j1,'interacts with more than one PL'
             write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             write(*,*) minv(:,j1)
          end if
          do m = 1, transpar%nPLs
             if (minv(m,j1).eq.j1) cblk(j1) = m
          end do
       end do

      else

         cblk=transpar%cblk              !DAR!!!

      end if   
  
       if (id0) then
          write(*,*) ' Structure info:'
          write(*,*) ' Number of PLs:',nbl 
          write(*,*) ' PLs coupled to contacts:',cblk(1:ncont)       !DAR
          write(*,*)
       endif

    end if

    call init_structure(negf, ncont, cont_end, surf_end, nbl, PL_end, cblk)

    deallocate(PL_end)
    deallocate(plcont)
    deallocate(atomst)
    deallocate(cblk)
    deallocate(cont_end)
    deallocate(surf_end)
    deallocate(ind)
    deallocate(minv)
  
  end subroutine negf_init_str

  !------------------------------------------------------------------------------
  ! INTERFACE subroutine to call Density matrix computation
  !------------------------------------------------------------------------------
  subroutine negf_density(miter,spin,nkpoint,HH,SS,mu,DensMat,EnMat)
    
    integer, intent (in) :: miter          ! SCC step (used in SGF)
    integer, intent (in) :: spin           ! spin component (SGF)
    integer, intent (in) :: nkpoint        ! nk point (used in SGF)
    type(z_CSR), intent(in) :: HH          ! Hamiltonian
    type(z_CSR), intent(in) :: SS          ! Overlap
    real(dp), intent(in) :: mu(:)       
    type(z_CSR), optional :: DensMat   ! Density matrix (See NOTE)
    type(z_CSR), optional :: EnMat     ! Energy weighted DM (See NOTE)

    type(lnParams) :: params
    integer :: ncont

    call get_params(negf, params)
      
    params%iteration = miter
    params%kpoint = nkpoint
    params%spin = spin
    params%DorE='N'
    ncont=negf%str%num_conts
    params%mu(1:ncont) = mu(1:ncont)
    
    if(present(DensMat)) then
       params%DorE = 'D'
       call set_params(negf,params)
       call pass_DM(negf,rho=DensMat)
       if (id0.and.params%verbose.gt.30) then
         write(*,'(73("="))')
         write(*,*) '                    COMPUTING DENSITY MATRIX      '
         write(*,'(73("="))') 
       endif
    endif
    if(present(EnMat)) then
       params%DorE = 'E'
       call set_params(negf,params)
       call pass_DM(negf,rhoE=EnMat)
       if (id0.and.params%verbose.gt.30) then
         write(*,'(73("="))')
         write(*,*) '                 COMPUTING E-WEIGHTED DENSITY MATRIX '
         write(*,'(73("="))') 
       endif
    endif
    if (present(DensMat).and.present(EnMat)) then
       params%DorE  = 'B'
       call set_params(negf,params)
       stop 'UNSUPPORTED CASE in negf_density'
    endif

    if (params%DorE.eq.'N') return 
    ! ---------------------------------------------

    call pass_HS(negf,HH,SS)
    
    call compute_density_dft(negf)
   
    call destroy_matrices(negf)
    
    if (id0.and.params%verbose.gt.30) write(*,'(73("*"))')

  end subroutine negf_density

  !------------------------------------------------------------------------------
  ! INTERFACE subroutine to call ldos computation
  !------------------------------------------------------------------------------
  subroutine negf_ldos(HH,SS,spin,kpoint,wght,ledos)
    type(z_CSR), intent(in) :: HH, SS
    integer, intent(in) :: spin      ! spin index 
    integer, intent(in) :: kpoint        ! kp index 
    real(dp), intent(in) :: wght      ! kp weight 
    real(dp), dimension(:,:), pointer :: ledos
    type(lnParams) :: params

    call get_params(negf, params)

    if (id0) then
      write(*,*)
      write(*,'(73("="))')
      write(*,*) '                   COMPUTING  LOCAL  DOS          '
      write(*,'(73("="))') 
      write(*,*)
    end if

    params%spin = spin
    params%kpoint = kpoint
    params%wght = wght

    if (associated(negf%ldos_mat)) then
      call log_deallocatep(negf%ldos_mat)
      negf%ldos_mat=> null()
    endif
    ledos=>null()
    
    call pass_HS(negf,HH,SS)

    call compute_ldos(negf)
    
    call destroy_matrices(negf)

    call associate_ldos(negf, ledos)

  end subroutine negf_ldos
  

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

  !------------------------------------------------------------------------------
  ! ALEX: DAR routines to set H and S have been moved here 
  !------------------------------------------------------------------------------
  subroutine prepare_HS(H_dev,S_dev,HH,SS)
    real(dp), dimension(:,:) :: H_dev, S_dev    
    type(z_CSR), intent(inout) :: HH, SS   

    if (negf%tOrthonormal) call Orthogonalization(H_dev, S_dev)

    if (negf%tOrthonormalDevice) call Orthogonalization_dev(H_dev, S_dev)

    if (negf%tOrthonormal.or.negf%tOrthonormalDevice) then
      call MakeHHSS(H_dev,S_dev,HH,SS)     
    end if

    !!if(negf%tManyBody) call MakeHS_dev

  end subroutine
  !------------------------------------------------------------------------------
         
  !------------------------------------------------------------------------------
  ! INTERFACE subroutine to call current computation
  !------------------------------------------------------------------------------    
  subroutine negf_current(HH,SS,spin,kpoint,wght,tunn,ledos,currents)
    
    type(z_CSR), intent(in) :: HH, SS   
    integer, intent(in) :: spin      ! spin index 
    integer, intent(in) :: kpoint        ! kp index 
    real(dp), intent(in) :: wght      ! kp weight 
    real(dp), dimension(:,:), pointer :: tunn
    real(dp), dimension(:,:), pointer :: ledos
    real(dp), dimension(:), pointer :: currents 

    type(lnParams) :: params   

    call get_params(negf, params)

    params%spin = spin
    params%kpoint = kpoint
    params%wght = wght
    
    call set_params(negf, params)

    call pass_HS(negf,HH,SS)

    !DAR begin
    !if(negf%tWrite_ldos) call WriteLDOS
    
    if (id0.and.negf%verbose.gt.30) then
       write(*,*)
       write(*,'(80("="))')
       write(*,*) '                          LibNEGF: Current calculation' 
       write(*,'(80("="))') 
       !write(*,*)
    endif

    print*,'spin:',spin,'kpoint:',kpoint
    call compute_current(negf)  
   
    !call write_tunneling_and_dos(negf)

    ! Associate internal negf pointers to local pointers 
    call associate_current(negf, currents)
    if (.not.associated(currents)) STOP 'Internal error: currVec not associated'     
    print*,'size currVec',size(currents)

    call associate_ldos(negf, ledos)  
    call associate_transmission(negf, tunn)

    !DAR begin
    if (id0.and.negf%verbose.gt.30) then
       write(*,*)
       write(*,'(80("="))')
       write(*,*) '                           LibNEGF: Current finished'  
       write(*,'(80("="))') 
       !write(*,*)
    endif
    !DAR end
    
  end subroutine negf_current



  !> Calculates density matrix with Green's functions
  !!
  subroutine calcdensity_green(iSCCIter, mpicomm, groupKS, ham, over, &
      & iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, &
      & cellVec, orb, nEl, tempElec, kPoints, kWeights, &
      & rho, Eband, Ef, E0, TS, mu)
    
    integer, intent(in) :: iSCCIter
    type(mpifx_comm), intent(in) :: mpicomm
    integer, intent(in) :: groupKS(:,:)
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb 
    real(dp), intent(in) :: nEl(:), tempElec, kPoints(:,:), kWeights(:)
    real(dp), intent(out) :: rho(:,:)
    real(dp), intent(out) :: Eband(:), Ef(:), E0(:), TS(:)

    integer :: nSpin, nKPoint, nKS, iK, iS, iKS, ncont
    real(dp) :: prefac, EBandTmp, EfTmp, TSTmp, E0Tmp
    type(z_CSR) :: csrDens
    !! We need this now for different fermi levels in colinear spin
    !! Note: the spin polirized does not work with
    !! built-int potentials (the unpolarized does) in the poisson
    real(dp), intent(in) :: mu(:,:)
    type(lnParams) :: params

    call get_params(negf, params)

    nKS = size(groupKS, dim=2)
    nKPoint = size(kPoints, dim=2)
    nSpin = size(ham, dim=2)
    ncont = size(mu,1)
    rho = 0.0_dp
    ncont = size(mu,1)

    do iKS = 1, nKS
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS)
      
      call foldToCSR(csrHam, ham(:,iS), kPoints(:,iK), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)
      call foldToCSR(csrOver, over, kPoints(:,ik), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)
      
      call negf_density(iSCCIter, iS, iKS, csrHam, csrOver, mu(:,iS), DensMat=csrDens)
      
      ! NOTE:
      ! unfold adds up to rho the csrDens(k) contribution
      !
      call unfoldFromCSR(rho(:,iS), csrDens, kPoints(:,iK), &
          kWeights(iK), iAtomStart, iPair, iNeighbor, nNeighbor, &
          img2CentCell, iCellVec, cellVec, orb)

      call destruct(csrDens)

      !! Set some fake energies:
      Eband(iS) = 0.0_dp
      Ef(iS) = 0.0_dp
      TS(iS) = 0.0_dp
      E0(iS) = 0.0_dp

    end do

    do iS = 1, nSpin
      ! In place all-reduce of the density matrix
      call mpifx_allreduceip(mpicomm, rho(:,iS), MPI_SUM)
    end do
  
  end subroutine calcdensity_green

  !> Calculates E-density matrix with Green's functions
  !!
  subroutine calcEdensity_green(iSCCIter, mpicomm, groupKS, ham, over, &
      & iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, &
      & cellVec, orb, nEl, tempElec, kPoints, kWeights, rhoE, mu)
    
    integer, intent(in) :: iSCCIter
    type(mpifx_comm), intent(in) :: mpicomm
    integer, intent(in) :: groupKS(:,:)
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb   !Needs only orb%nOrbAtom, orb%mOrb    
    real(dp), intent(in) :: nEl(:), tempElec, kPoints(:,:), kWeights(:)
    real(dp), intent(out) :: rhoE(:)

    integer :: nSpin, nKPoint, nKS, iK, iS, iKS, ncont
    real(dp) :: prefac, EBandTmp, EfTmp, TSTmp, E0Tmp
    type(z_CSR) :: csrEDens
    !! We need this now for different fermi levels in colinear spin
    !! Note: the spin polirized does not work with
    !! built-int potentials (the unpolarized does) in the poisson
    !! I do not set the fermi because it seems that in libnegf it is 
    !! not really needed
    real(dp), intent(in) :: mu(:,:)
    type(lnParams) :: params

    call get_params(negf, params)

    nKS = size(groupKS, dim=2)
    nKPoint = size(kPoints, dim=2)
    nSpin = size(ham, dim=2)
    ncont = size(mu,1)
    rhoE = 0.0_dp

    
    do iKS = 1, nKS
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS)
      
      call foldToCSR(csrHam, ham(:,iS), kPoints(:,iK), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)
      call foldToCSR(csrOver, over, kPoints(:,ik), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)
      
      call negf_density(iSCCIter, iS, iKS, csrHam, csrOver, mu(:,iS), EnMat=csrEDens)

      ! NOTE:
      ! unfold adds up to rhoEPrim the csrEDens(k) contribution
      !
      call unfoldFromCSR(rhoE, csrEDens, kPoints(:,iK), &
          kWeights(iK), iAtomStart, iPair, iNeighbor, nNeighbor, &
          img2CentCell, iCellVec, cellVec, orb)

      call destruct(csrEDens)
      
    end do

    ! In place all-reduce of the density matrix
    call mpifx_allreduceip(mpicomm, rhoE, MPI_SUM)
    
  end subroutine calcEdensity_green

  !------------------------------------------------------------------------------
  ! INTERFACE subroutine to call current computation
  !------------------------------------------------------------------------------
  subroutine calcPDOS_green(mpicomm, groupKS, ham, over, &
      & iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, &
      & cellVec, orb, nEl, tempElec, kPoints, kWeights, ldosTot, writeLDOS)
    integer, intent(in) :: groupKS(:,:)
    type(mpifx_comm), intent(in) :: mpicomm
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: nEl(:), tempElec, kPoints(:,:), kWeights(:)
    real(dp), allocatable, intent(inout) :: ldosTot(:,:)
    logical, intent(in) :: writeLDOS

    real(dp), pointer    :: ldosMat(:,:)=>null()
    real(dp), allocatable :: ldosSKRes(:,:,:)
    integer :: iKS, iK, iS, nKS, nKPoint, nSpin, ii, jj, err
    type(lnParams) :: params

    call get_params(negf, params)
    
    nKS = size(groupKS, dim=2)
    nKPoint = size(kPoints, dim=2)
    nSpin = size(ham, dim=2)

    do iKS = 1, nKS
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS)
      
      call foldToCSR(csrHam, ham(:,iS), kPoints(:,iK), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)
     
      call foldToCSR(csrOver, over, kPoints(:,ik), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)

      call negf_ldos(csrHam, csrOver, iS, iK, kWeights(iK), ldosMat)
     
      call add_partial_results(mpicomm, ldosMat, ldosTot, ldosSKRes, iKS, nKS)  
      
    end do

    if (allocated(ldosTot)) then
      
      ! Write Total localDOS on a separate file (optional)
      if (id0 .and. writeLDOS) then
        open(65000,file='localDOS.dat')
        do ii=1,size(ldosTot,1)
          write(65000,*) (params%Emin+(ii-1)*params%Estep) * HAR, &
              (ldosTot(ii,jj), jj=1,size(ldosTot,2))
        enddo
        close(65000)
      endif
    else
      allocate(ldosTot(0,0))
    endif
    
    
  end subroutine calcPDOS_green

      
  
  !------------------------------------------------------------------------------
  ! INTERFACE subroutine to call current computation
  !------------------------------------------------------------------------------
  subroutine calc_current(mpicomm, groupKS, ham, over, &
      & iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, &
      & cellVec, orb, nEl, tempElec, kPoints, kWeights, tunnTot,&
      & ldosTot, currTot, writeTunn, writeLDOS, mu)
    
    integer, intent(in) :: groupKS(:,:)
    type(mpifx_comm), intent(in) :: mpicomm
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: nEl(:), tempElec, kPoints(:,:), kWeights(:)
    real(dp), allocatable, intent(inout) :: tunnTot(:,:), ldosTot(:,:)
    real(dp), allocatable :: tunnSKRes(:,:,:), ldosSKRes(:,:,:)
    real(dp), allocatable, intent(inout) :: currTot(:)
    logical, intent(in) :: writeLDOS
    logical, intent(in) :: writeTunn
    !! We need this now for different fermi levels in colinear spin
    !! Note: the spin polirized does not work with
    !! built-int potentials (the unpolarized does) in the poisson
    !! I do not set the fermi because it seems that in libnegf it is 
    !! not really needed
    real(dp), intent(in) :: mu(:,:)
    
    real(dp), pointer    :: tunnMat(:,:)=>null()
    real(dp), pointer    :: ldosMat(:,:)=>null()
    real(dp), pointer    :: currVec(:)=>null()    
    integer :: iKS, iK, iS, nKS, nKPoint, nSpin, ii, jj, err, ncont
    type(unit) :: unitOfEnergy        ! Set the units of H
    type(unit) :: unitOfCurrent       ! Set desired units for Jel
    type(lnParams) :: params

    integer :: i,j,k,NumStates,icont                                       !DAR
    real(dp), dimension(:,:), allocatable :: H_all, S_all
    call get_params(negf, params)
 
    unitOfEnergy%name = "H"
    unitOfCurrent%name = "A"

    nKS = size(groupKS, dim=2)
    nKPoint = size(kPoints, dim=2)
    nSpin = size(ham, dim=2)
    ncont = size(mu,1)

    if (id0.and.params%verbose.gt.30) then
      write(*,*)
      write(*,'(80("="))')
      write(*,*) '                            COMPUTATION OF TRANSPORT         '
      write(*,'(80("="))') 
      write(*,*)
    end if

    do iKS = 1, nKS 
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS)

      params%mu(:ncont) = mu(:ncont,iS) 
     
      call set_params(negf, params)
      
      if (negf%NumStates.eq.0) negf%NumStates=csrHam%ncol

      !*** ORTOGONALIZATIONS ***
      ! THIS MAKES SENSE ONLY FOR A REAL MATRIX: k=0 && collinear spin
      ! TEMPORARY HACKING: we save H, S on files as square-matrices
      !                    then reload them and perfom Lowdin tranf.  
      if (all(kPoints(:,iK) .eq. 0.0_dp) .and. &
        &(negf%tOrthonormal .or. negf%tOrthonormalDevice)) then
            
          call writeSparseAsSquare('H_dftb.mtr', ham(:,iS)*HAR, iNeighbor, nNeighbor, &
                                       &iAtomStart, iPair, img2CentCell)
          write(*,"(' Hamiltonian is written to the file ',A)")trim('H_dftb.mtr')
          
          call writeSparseAsSquare('S_dftb.mtr', over, iNeighbor, nNeighbor, &
                                        &iAtomStart, iPair, img2CentCell)
          write(*,"(' Overlap is written to the file ',A)")trim('S_dftb.mtr')
       
          NumStates = negf%NumStates
        
          write(*,"(' Hamiltonian is red from the file ',A)")trim('H_dftb.mtr')
          open(11,file='H_dftb.mtr',action="read")
          read(11,*);read(11,*);read(11,*);read(11,*);read(11,*)
        
          if(.not.allocated(H_all)) allocate(H_all(NumStates,NumStates))
          H_all=0.0_dp
          do i = 1, NumStates
             read(11,*) H_all(i,1:NumStates)
          end do
          close(11)
          H_all = H_all/HAR
          write(*,"(' Overlap is red from the file ',A)")trim('S_dftb.mtr')
          open(11,file='S_dftb.mtr',action="read")
          read(11,*);read(11,*);read(11,*);read(11,*);read(11,*)
          
          if (.not.allocated(S_all)) allocate(S_all(NumStates,NumStates))
          S_all=0.0_dp
          do i = 1, NumStates
             read(11,*) S_all(i,1:NumStates)
          end do
          close(11)
           
          call prepare_HS(H_all,S_all,csrHam,csrOver)         

      else
  
        call foldToCSR(csrHam, ham(:,iS), kPoints(:,iK), iAtomStart, &
           &iPair, iNeighbor, nNeighbor, img2CentCell, &
           &iCellVec, cellVec, orb)
     
        call foldToCSR(csrOver, over, kPoints(:,ik), iAtomStart, &
           &iPair, iNeighbor, nNeighbor, img2CentCell, &
           &iCellVec, cellVec, orb)

      end if     

      call negf_current(csrHam, csrOver, iS, iK, kWeights(iK), &
           tunnMat, ldosMat, currVec)

      !-------------------------------------------------------------------------
      if(.not.allocated(currTot)) then
        allocate(currTot(size(currVec)), stat=err)
        if (err/=0) STOP 'Allocation error (currTot)'
        currTot = 0.0_dp
      endif
      currTot = currTot + currVec

      !GUIDE: tunnMat libNEGF output stores contact Tunneling T(iE, i->j) 
      !       tunnMat is MPI distributed on energy points (0.0 on other nodes)
      !       tunnTot MPI gather partial results and accumulate k-summation  
      !       tunnSKRes stores tunneling for all k-points and spin: T(iE, i->j, iSK)  
      call add_partial_results(mpicomm, tunnMat, tunnTot, tunnSKRes, iKS, nKS) 
      
      call add_partial_results(mpicomm, ldosMat, ldosTot, ldosSKRes, iKS, nKS) 

    end do
   
    call mpifx_allreduceip(mpicomm, currTot, MPI_SUM)

    call mpifx_barrier(mpicomm) 

    ! converts from internal atomic units into A
    currTot = currTot * convertCurrent(unitOfEnergy, unitOfCurrent) 
    
    if (id0) then 
      do ii=1, size(currTot)
        write(*,'(1x,a,i3,i3,a,ES14.5,a,a)') &
             & ' contacts: ',params%ni(ii),params%nf(ii), &
             & ' current: ', currTot(ii),' ',unitOfCurrent%name
      enddo
    endif
    
    if (allocated(tunnTot)) then
      ! Write Total tunneling on a separate file (optional)
      if (id0 .and. writeTunn) then
        call write_file(negf, tunnTot, tunnSKRes, 'tunneling', &
                         & groupKS, kpoints, kWeights)
      endif 
      if (allocated(tunnSKRes)) deallocate(tunnSKRes)
    else
      ! needed to avoid some segfault     
      allocate(tunnTot(0,0))
    endif

    !!DAR begin - print tunn_mat_bp
    !#if (negf%tZeroCurrent) then                                                
    !#  call mpifx_allreduceip(mpicomm, negf%tunn_mat, MPI_SUM)        
    !#  if (id0) call write_file(negf, negf%tunn_mat, tunnSKRes, 'tunneling_bp', &
    !#        & groupKS, kpoints, kWeights)
    !#  if (allocated(tunnSKRes)) deallocate(tunnSKRes)
    !#end if
    !DAR end

    if (allocated(ldosTot)) then
      ! Write Total localDOS on a separate file (optional)
      if (id0 .and. writeLDOS) then
        call write_file(negf, ldosTot, ldosSKRes, 'localDOS', &
                         & groupKS, kpoints, kWeights)
      end if
      if (allocated(ldosSKRes)) deallocate(ldosSKRes)
    else
      ! needed to avoid some segfault     
      allocate(ldosTot(0,0))
    end if
   
  end subroutine calc_current

  !----------------------------------------------------------------------------
  !   utility to allocate and sum partial results
  !----------------------------------------------------------------------------
  subroutine add_partial_results(mpicomm, pMat, pTot, pSKRes, iK, nK)
    type(mpifx_comm), intent(in) :: mpicomm
    real(dp), pointer :: pMat(:,:)
    real(dp), allocatable :: pTot(:,:)
    real(dp), allocatable :: pSKRes(:,:,:)
    integer, intent(in) :: iK, nK
    real(dp), allocatable :: tmpMat(:,:)
      
    integer :: err 

    if (associated(pMat)) then
      allocate(tmpMat(size(pMat,1), size(pMat,2)), stat=err)
      if (err/=0) STOP 'Allocation error (tmpMat)'
      tmpMat = pMat
      call mpifx_allreduceip(mpicomm, tmpMat, MPI_SUM)
      if(.not.allocated(pTot)) then 
        allocate(pTot(size(pMat,1), size(pMat,2)), stat=err)
        if (err/=0) STOP 'Allocation error (tunnTot)'
        pTot = 0.0_dp
      end if
      pTot = pTot + tmpMat
   
      if (nK > 1) then
        if(.not.allocated(pSKRes)) then 
          allocate(pSKRes(size(pMat,1), size(pMat,2), nK), stat=err)
          if (err/=0) STOP 'Allocation error (tunnSKRes)'
          pSKRes = 0.0_dp
        endif
        pSKRes(:,:,iK) = tmpMat(:,:)
      end if  
      deallocate(tmpMat)
    end if
  end subroutine add_partial_results
  ! ----------------------------------------------------------------------------
  !    utility to write tunneling or ldos on files 
  ! ----------------------------------------------------------------------------
  subroutine write_file(negf, pTot, pSKRes, filename, groupKS, kpoints, kWeights)
    type(TNegf) :: negf
    real(dp), intent(in) :: pTot(:,:)
    real(dp), intent(in) :: pSKRes(:,:,:)
    character(*), intent(in) :: filename
    integer, intent(in) :: groupKS(:,:)
    real(dp), intent(in) :: kPoints(:,:)
    real(dp), intent(in) :: kWeights(:)

    integer :: ii, jj, nKS, iKS, nK, nS, iK, iS
    type(lnParams) :: params

    call get_params(negf, params)

    nKS = size(groupKS, dim=2)
    nK = size(kpoints, dim=2)
    nS = nKS/nK 

    open(65000,file=trim(filename)//'.dat')
    do ii=1,size(pTot,1)
      write(65000,'(f20.8)',ADVANCE='NO') (params%Emin+ii*params%Estep) * HAR
      do jj=1,size(pTot,2)
        !write(65000,'(es20.8)',ADVANCE='NO') pTot(ii,jj)
        write(65000,'(f20.8)',ADVANCE='NO') pTot(ii,jj)                     !DAR
      enddo
      write(65000,*)
    enddo
    close(65000)

    if (nKS.gt.1) then
      open(65000,file=trim(filename)//'_kpoints.dat')
      write(65000,*)  '# NKpoints = ', nK
      write(65000,*)  '# NSpin = ', nS
      write(65000,*)  '# Energy [eV], <spin k1 k2 k3 weight> '
      write(65000,'(A1)', ADVANCE='NO') '# '
      do iKS = 1,nKS
        iK = groupKS(1,iKS)
        iS = groupKS(2,iKS)    
        write(65000,'(i5.2)', ADVANCE='NO') iS
        write(65000,'(es15.5, es15.5, es15.5, es15.5)', ADVANCE='NO') kpoints(:,iK),&
                                                                      & kWeights(iK) 
      end do
      write(65000,*)
      do ii=1,size(pSKRes(:,:,1),1)
        write(65000,'(f20.8)',ADVANCE='NO') (params%Emin+ii*params%Estep) * HAR
        do jj=1,size(pSKRes(:,:,1),2)
          do iKS = 1,nKS
            write(65000,'(es20.8)',ADVANCE='NO') pSKRes(ii,jj, iKS)
          enddo
          write(65000,*)
        enddo
      enddo
      close(65000)
    end if

  end subroutine write_file

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! THIS is a first version of local current computation. 
  ! It has been placed here since it depends on internal representations of DFTB
  !       
  ! NOTE: Limited to non-periodic systems             s !!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine local_currents(mpicomm, groupKS, ham, over, &
      & iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, &
      & cellVec, orb, kPoints, kWeights, coord, dumpDens, chempot)
    
    type(mpifx_comm), intent(in) :: mpicomm
    integer, intent(in) :: groupKS(:,:)
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: kPoints(:,:), kWeights(:)
    real(dp), intent(in) :: coord(:,:)
    logical, intent(in) :: dumpDens
    !! We need this now for different fermi levels in colinear spin
    !! Note: the spin polirized does not work with
    !! built-int potentials (the unpolarized does) in the poisson
    !! I do not set the fermi because it seems that in libnegf it is 
    !! not really needed
    real(dp), intent(in) :: chempot(:,:)
    
    integer :: n,m, mu, nu, nAtom, irow, nrow, ncont
    integer :: nKS, nKPoint, nSpin, iK, iKS, iS, inn, startn, endn, morb
    real(dp), dimension(:), allocatable :: Im 
    integer, dimension(:,:), allocatable :: nneig
    integer, dimension(:), allocatable :: nn
    integer, parameter :: NMAX=40
    complex(dp) :: c1,c2
    character(1) :: sp
    integer :: iSCCiter=2
    type(z_CSR) :: csrDens, csrEDens
    type(lnParams) :: params

    call get_params(negf, params)

    if (id0.and.params%verbose.gt.30) then
      write(*,*)
      write(*,'(73("="))')
      write(*,*) '                   COMPUTING LOCAL CURRENTS          '
      write(*,'(73("="))') 
      write(*,*)
    endif

    nKS = size(groupKS, dim=2)
    nKPoint = size(kPoints, dim=2)
    nSpin = size(ham, dim=2)
    nAtom = size(orb%nOrbAtom)
    
    do iKS = 1, nKS
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS)
      
      ! We need to recompute Rho and RhoE .....
      call foldToCSR(csrHam, ham(:,iS), kPoints(:,1), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)
      call foldToCSR(csrOver, over, kPoints(:,1), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)

      call negf_density(iSCCIter, iS, iKS, csrHam, csrOver, chempot(:,iS), DensMat=csrDens)

      call negf_density(iSCCIter, iS, iKS, csrHam, csrOver, chempot(:,iS), EnMat=csrEDens)

      call mpifx_allreduceip(mpicomm, csrDens%nzval, MPI_SUM)
      
      call mpifx_allreduceip(mpicomm, csrEDens%nzval, MPI_SUM)
      ! Gather is done here 

      !if (dumpDens) then
      !  open(65000, file='dens.csr',form='formatted')
      !  call print_mat(65000, csrDens, .true.)
      !  close(65000)
      !end if
      
      allocate(nneig(nAtom,NMAX))
      allocate(nn(nAtom))
      call symmetrize_neiglist(nAtom,img2CentCell,iNeighbor,nNeighbor,coord&
          &,nneig,nn)
      
      
      if (iS .eq. 1) sp = 'u'
      if (iS .eq. 2) sp = 'd'
      open(207,file='lcurrents_'//sp//'.dat')

      do m = 1, nAtom

        allocate(Im(nn(m)))
        Im(:) = 0.0_dp

        mOrb = orb%nOrbAtom(m)
        iRow = iAtomStart(m)

        write(207,'(I5,3(F12.6),I3)',advance='NO') m, coord(:,m), nn(m) 

        do inn = 1, nn(m)
          n = nneig(m,inn)
          startn = iAtomStart(n)
          endn = startn + orb%nOrbAtom(n) - 1

          ! tracing orbitals of atoms  n  m
          ! More efficient without getel ?    
          do mu = iRow, iRow+mOrb-1
            do nu = startn, endn
              c1=conjg(getel(csrDens,mu,nu))
              c2=conjg(getel(csrEDens,mu,nu))
              Im(inn)=Im(inn)+aimag(getel(csrHam,mu,nu)*c1 - &
                  &getel(csrOver,mu,nu)*c2)
            enddo
          enddo
          ! pi-factor  comes from  Gn = rho * pi 
          write(207,'(I5,ES17.8)',advance='NO') n, 2.d0*params%g_spin*pi*eovh*Im(inn)

        enddo

        write(207,*)

        deallocate(Im) 

      enddo
    enddo

    close(207)
    call destruct(csrDens)
    call destruct(csrEDens)
    deallocate(nneig)
    deallocate(nn)

  end subroutine local_currents

  !--------------------------------------------------------------
  ! Neighbor is non-symmetric: i->j  j>i
  ! The following routine symmetrizes the neighbor list 
  !
  subroutine symmetrize_neiglist(nAtom,img2CentCell,iNeighbor,nNeighbor,coord&
      &,nneig,nn)
    integer, intent(in) :: nAtom
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    real(dp), intent(in) :: coord(:,:)
    integer, dimension(:,:) :: nneig
    integer, dimension(:) :: nn    
    
    real(dp) :: dist, dr(3)
    integer :: m, n, inn, ii, jj, morb
    integer, parameter :: NMAX=40
    
    nn=0
    nneig=0

    do m = 1, nAtom
      do inn = 1, nNeighbor(m)
        if (inn.gt.NMAX) exit
        n = img2CentCell(iNeighbor(inn,m))
        if(nn(m).le.NMAX-1) then
          nn(m)=nn(m)+1 
          nneig(m,nn(m))=n
          ! sort by distances
          jj = nn(m)
          dr(:) = coord(:,m)-coord(:,nneig(m,jj))
          dist= dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
          do ii = nn(m)-1,1,-1
            dr(:) = coord(:,m)-coord(:,nneig(m,ii))
            if (dist.lt.dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)) then
              morb=nneig(m,ii)     
              nneig(m,ii)=nneig(m,jj)
              nneig(m,jj)=morb
              jj = jj - 1 
            endif
          enddo
          ! -----------------
        endif
        if(nn(n).le.NMAX-1) then
          nn(n)=nn(n)+1 
          nneig(n,nn(n))=m
          ! sort by distances
          jj = nn(n)
          dr(:) = coord(:,n)-coord(:,nneig(n,jj))
          dist= dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3) 
          do ii = nn(n)-1,1,-1
            dr(:) = coord(:,n)-coord(:,nneig(n,ii))
            if (dist.lt.dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)) then
              morb=nneig(n,ii)     
              nneig(n,ii)=nneig(n,jj)
              nneig(n,jj)=morb
              jj=jj-1
            endif
          enddo
          ! -----------------
        endif
      enddo
    enddo
    !--------------------------------------------------------------
    !nn = 36 
    !print*,'nneig:'
    !do m =1,nAtom
    !   write(*,'(8(i6))',advance='NO') m
    !   do jj=1,nn(m)
    !     dr(:) = coord(:,m)-coord(:,nneig(m,jj))
    !     dist= sqrt(dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)) 
    !     write(*,'(i5,ES15.6)',advance='NO') nneig(m,nn(jj)),dist
    !   enddo 
    !   write(*,*) 
    !enddo
    !print*,'done'
    !--------------------------------------------------------------

  end subroutine symmetrize_neiglist

  !-----------------------------------------------------------------------------
  !DAR begin - ReadDFTB, ReadModel, test_negf, orthogonalization, etc.
  !-----------------------------------------------------------------------------

  subroutine ReadDFTB(H,S)
    real(dp), dimension(:,:), allocatable :: H, S 
    integer :: i,j,k,l,m,n,NumStates

    NumStates = negf%NumStates

    if (id0) write(*,"(' Hamiltonian is red from the file ',A)")trim('H_dftb.mtr')
    open(11,file='H_dftb.mtr',action="read")
    read(11,*);read(11,*);read(11,*);read(11,*);read(11,*)
    if(.not.allocated(H)) allocate(H(NumStates,NumStates))
    H=0.0_dp
    do i=1,NumStates
       read(11,*) H(i,1:NumStates)
    end do
    close(11)
    H=H/HAR
    if (id0) write(*,"(' Overlap is red from the file ',A)")trim('S_dftb.mtr')
    open(11,file='S_dftb.mtr',action="read")
    read(11,*);read(11,*);read(11,*);read(11,*);read(11,*)
    if(.not.allocated(S)) allocate(S(NumStates,NumStates))
    S=0.0_dp
    do i=1,NumStates
       read(11,*) S(i,1:NumStates)
    end do
    close(11)

  end subroutine ReadDFTB

  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------

  subroutine ReadModel(H,S)
    real(dp), dimension(:,:), allocatable :: H, S 
    integer :: i,j,k,l,m,n,NumStates

    NumStates = negf%NumStates

    if (id0) write(*,"(' Hamiltonian is red from the file ',A)")trim('H.mtr')
    open(11,file='H.mtr',action="read")
    if(.not.allocated(H)) allocate(H(NumStates,NumStates))
    H=0.0_dp
    do i=1,NumStates
      read(11,*) H(i,1:NumStates)
    end do
    close(11)
    H=H/HAR
    if(.not.allocated(S)) allocate(S(NumStates,NumStates))
    S=0.0_dp
    do i=1,NumStates
      S(i,i)=1._dp
    end do

  end subroutine ReadModel

  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------

  subroutine ReadLibNEGF
   
    write(*,*) 'Import from the negf.in file'
    call read_negf_in(negf)
    call negf_partition_info(negf)

  end subroutine ReadLibNEGF

  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------

  subroutine MakeHHSS(H_all, S_all, HH, SS)
    real(dp), dimension(:,:), intent(in) :: H_all, S_all
    type(z_CSR), intent(inout) :: HH, SS
    integer :: i,j,k,l,m,n,NumStates, nnz

    NumStates = negf%NumStates
    
    nnz=0
    do i=1,NumStates
      do j=1,NumStates   
        if ((i.eq.j).or.(abs(H_all(i,j)).gt.0.00001)) then
          nnz = nnz+1
        end if
      end do
    end do

    call destroy(HH)
    call create(HH, NumStates, NumStates, nnz) 
    call destroy(SS)
    call create(SS, NumStates, NumStates, nnz) 

    HH%rowpnt(1)=1
    nnz=0
    do i=1,NumStates
       k=0
       do j=1,NumStates   
          if((i.eq.j).or.(abs(H_all(i,j)).gt.0.00001)) then
             k=k+1
             nnz=nnz+1            
             if(i.eq.j) then
                HH%nzval(nnz)= H_all(i,j)
                SS%nzval(nnz)= S_all(i,j)
             else
                HH%nzval(nnz)= H_all(i,j)
                SS%nzval(nnz)= S_all(i,j)
             end if
             HH%colind(nnz)=j
             !print *, HH%nnz, negf%H_all(i,j),  HH%nzval(HH%nnz), SS%nzval(HH%nnz)   !debug
          end if
       end do
       HH%rowpnt(i+1)=HH%rowpnt(i)+k
    end do

    SS%colind=HH%colind
    SS%rowpnt=HH%rowpnt

  end subroutine MakeHHSS
  
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  

  subroutine WriteLDOS
    
    if (id0.and.negf%verbose.gt.30) then
    write(*,*)
    write(*,'(80("="))')
    write(*,*) '                   LibNEGF: equilibrium (!) LDOS calculation'  
    write(*,'(80("="))') 
    write(*,*)
    endif
    call compute_ldos(negf)
    if (id0.and.negf%verbose.gt.30) then
    write(*,*)
    write(*,'(80("="))')
    write(*,*) '                     LibNEGF: equilibrium (!) LDOS finished'   
    write(*,'(80("="))') 
    !write(*,*)
    endif

  end subroutine WriteLDOS

  !-----------------------------------------------------------------------------
  
  subroutine check_negf_params

    use globals, only : LST
    
    Integer :: ncont, nbl, ii, jj, ist, iend
    Integer, dimension(:), allocatable :: PL_end, cont_end, surf_end
    character(32) :: tmp
    character(LST) :: file_re_H, file_im_H, file_re_S, file_im_S
    integer, dimension(:), allocatable :: cblk

    open(10,file='log',action="write")
    
    write(10,"('negf%H%nnz               = ',I6)")negf%H%nnz
    write(10,"('negf%H%nrow              = ',I6)")negf%H%nrow
    write(10,"('negf%H%ncol              = ',I6)")negf%H%ncol
    write(10,"('negf%H%nzval             = ',10000000F8.4)")negf%H%nzval
    write(10,"('negf%H%colind            = ',10000000I6)")negf%H%colind
    write(10,"('negf%H%rowpnt            = ',10000000I6)")negf%H%rowpnt
    write(10,"('negf%isSid               = ',L6)")negf%isSid
    write(10,"('negf%S%nnz               = ',I6)")negf%S%nnz
    write(10,"('negf%S%nrow              = ',I6)")negf%S%nrow
    write(10,"('negf%S%ncol              = ',I6)")negf%S%ncol
    write(10,"('negf%S%nzval             = ',10000000F8.4)")negf%S%nzval
    write(10,"('negf%S%colind            = ',10000000I6)")negf%S%colind
    write(10,"('negf%S%rowpnt            = ',10000000I6)")negf%S%rowpnt
    write(10,"('negf%str%num_conts       = ',I6)")negf%str%num_conts
    write(10,"('negf%str%num_PLs         = ',I6)")negf%str%num_PLs
    write(10,"('negf%str%active_cont     = ',I6)")negf%str%active_cont
    write(10,"('negf%str%mat_B_start     = ',I6,I6)")negf%str%mat_B_start
    write(10,"('negf%str%mat_C_start     = ',I6,I6)")negf%str%mat_C_start
    write(10,"('negf%str%mat_C_end       = ',I6,I6)")negf%str%mat_C_end
    write(10,"('negf%str%cblk            = ',I6,I6)")negf%str%cblk
    write(10,"('negf%str%cont_dim        = ',I6,I6)")negf%str%cont_dim
    write(10,"('negf%str%mat_PL_start    = ',10000000I6)")negf%str%mat_PL_start
    write(10,"('negf%str%mat_PL_end      = ',10000000I6)")negf%str%mat_PL_end
    write(10,"('negf%str%central_dim     = ',I6)")negf%str%central_dim
    write(10,"('negf%str%total_dim       = ',I6)")negf%str%total_dim
    write(10,"('negf%Ec, negf%Ev         = ',3F8.4)")negf%Ec, negf%Ev
    write(10,"('negf%DeltaEc, %DeltaEv   = ',3F8.4)")negf%DeltaEc, negf%DeltaEv
    write(10,"('negf%Emin, %Emax, %Estep = ',3F8.4)")negf%Emin, negf%Emax, negf%Estep
    write(10,"('negf%kbT_dm              = ',10000000E12.4)")negf%kbT_dm
    write(10,"('negf%kbT_t               = ',10000000E12.4)")negf%kbT_t
    write(10,"('negf%mu_n                = ',10000000F8.4)")negf%mu_n
    write(10,"('negf%mu_p                = ',10000000F8.4)")negf%mu_p
    write(10,"('negf%mu                  = ',10000000F8.4)")negf%mu
    write(10,"('negf%delta               = ',10000000E12.4)")negf%delta

    write(10,*)negf%wght
    write(10,*)negf%Np_n(1:2)
    write(10,*)negf%Np_p(1:2)
    write(10,*)negf%Np_real(1)
    write(10,*)negf%n_kt
    write(10,*)negf%n_poles
    write(10,*)negf%g_spin
    write(10,*)negf%nLDOS
    write(10,*)negf%LDOS(1)%indexes

    close (10)

    write(*,*)
    write(*,"(' negf parameters are written to the log file')")
        
  end subroutine check_negf_params
  !-----------------------------------------------------------------------------
 
  !-----------------------------------------------------------------------------
  subroutine orthogonalization(H,S)

    real(dp), dimension(:,:) :: H, S 
    integer :: i,m,n1_first,n1_last,n2_first,n2_last
    integer :: INFO, N
    real(dp), allocatable :: A(:,:), WORK(:), W(:)
    real(dp), allocatable :: B(:,:),C(:,:)

    N=size(H,1)
    if (N .ne. negf%NumStates) then
      if (id0) write(*,*) 'orthogonalization: negf init NumStates error'
      stop
    end if

    allocate(A(N,N),WORK(3*N),W(N))
    allocate(B(N,N),C(N,N))
    W=0.0_dp
    WORK=0.0_dp
    A=0.0_dp
    B=0.0_dp
    C=0.0_dp
    
    A=S

    call DSYEV('V','U',N,A,N,W,WORK,3*N,INFO )

    !print  *,'U matrix, Eigenvectors for S diagonalization'
    !do i=1,N
    !   write(*,*)A(i,1:N)
    !end do

    !print *,'U matrix unitarity check'
    !B=matmul(transpose(A),A)
    !do i=1,N
    !   write(*,*)B(i,1:N)
    !end do

    B=matmul(transpose(A),matmul(S,A))
    
    do i=1,N
      B(i,i)=1.0_dp/sqrt(B(i,i))
    end do

    !C=sqrt(S)
    C=matmul(A,matmul(B,transpose(A)))

    !print *,'sqrt(S) inverted'
    !do i=1,N
    !  write(*,*) C(i,1:N)
    !end do

    !print *,'S unity check'
    !B=matmul(transpose(C),matmul(S,C))
    !do i=1,N
    !   write(*,*) B(i,1:N)
    !end do

    !print *,'H_dftb before orthogonalization'
    !do i=1,N
    !   write(*,*) H(i,1:N)
    !end do

    H = matmul(transpose(C),matmul(H,C))

    !print *,'H_dftb_orth before replacement'
    !do i=1,N
    !   write(*,*) H(i,1:N)
    !end do

    ! COPY THE FIRST CONTACT PL ONTO THE SECOND 
    do m = 1, negf%str%num_conts
       n1_first = negf%str%mat_B_start(m)
       n1_last = (negf%str%mat_C_end(m)+negf%str%mat_B_start(m))/2
       n2_first = n1_last + 1
       n2_last = negf%str%mat_C_end(m)
       !print *,n1_first,n1_last,n2_first,n2_last
       H(n2_first:n2_last,n2_first:n2_last) = H(n1_first:n1_last,n1_first:n1_last)
    end do

    !print *,'H_dftb_orth after replacement'
    !do i=1,N
    !   write(*,*) H(i,1:N)
    !end do

    S = 0.0_dp
    do i=1,N
      S(i,i) = 1.0_dp
    end do

    !Save H_dftb_orth.mtr to file
    !open(12,file='H_dftb_orth.mtr',action="write")
    !do i = 1,N
    !  write(12,*) H(i,1:N)*HAR
    !end do
    !close(12)

    write(*,"(' Löwdin orthogonalization is done! ')")
    
  end subroutine orthogonalization

  !-----------------------------------------------------------------------------

  subroutine orthogonalization_dev(H, S)

    real(dp), dimension(:,:) :: H, S 
    integer :: i,m,n1_first,n1_last,n2_first,n2_last
    integer :: INFO, N, N2
    real(dp), allocatable :: A(:,:), WORK(:), W(:)
    real(dp), allocatable :: B(:,:), U(:,:), C(:,:)


    N=size(H,1)
    if (N .ne. negf%NumStates) then
      if (id0) write(*,*) 'orthogonalization: negf init NumStates error'
      stop
    end if
    
    N2=negf%str%central_dim

    allocate(A(N2,N2),WORK(3*N2),W(N2))
    allocate(B(N2,N2), U(N2,N2))
    W=0.0_dp
    WORK=0.0_dp
    
    A = S(1:N2,1:N2)
    U = A

    call DSYEV('V','U',N2,U,N2,W,WORK,3*N2,INFO )
 
    !print  *,'U matrix, Eigenvectors for S diagonalization'
    !do i=1,N2
    !   write(*,*) U(i,1:N2)
    !end do

    !U matrix unitarity check
    !B = matmul(transpose(U),U)
    !do i=1,N2
    !  if (abs(B(i,i)-1.0_dp) > 1.d-9 .or. any(abs(B(i,1:i-1)) > 1.d-9 ) then
    !    print*, 'ERROR: U is not unitary', B(i,:)    
    !    stop
    !  end if
    !end do

    B = matmul(transpose(U),matmul(A,U))

    do i=1,N2
      B(i,i)=1.0_dp/sqrt(B(i,i))
    end do


    ! Now A = S^-1/2 
    A = matmul(U,matmul(B,transpose(U)))
    !print *,'sqrt(S) inverted'
    !do i=1,N2
    !  write(*,*) A(i,1:N2)
    !end do

    deallocate(U, B)
    
    allocate(C(N,N))
    C=0.0_dp
    do i = N2+1,N
      C(i,i)=1.0_dp
    end do
    C(1:N2,1:N2) = A

    deallocate(A)

    !C=sqrt(S) big matrix
    !print *,'C=sqrt(S) big matrix'
    !do i=1,N
    !   write(*,*) C(i,1:N)
    !end do
        
    !print *,'H_dftb before orthogonalization'
    !do i=1,N
    !   write(*,*) H(i,1:N)
    !end do

    H = matmul(transpose(C),matmul(H,C))
    S = matmul(transpose(C),matmul(S,C))

    !print *,'H_dftb_orth before replacement'
    !do i=1,N
    !   write(*,*) H(i,1:N)
    !end do

    ! COPY THE FIRST CONTACT PL ONTO THE SECOND 
    do m=1,negf%str%num_conts
       n1_first=negf%str%mat_B_start(m)
       n1_last=(negf%str%mat_C_end(m)+negf%str%mat_B_start(m))/2
       n2_first=n1_last+1
       n2_last=negf%str%mat_C_end(m)
       !print *,n1_first,n1_last,n2_first,n2_last
       H(n2_first:n2_last,n2_first:n2_last) = H(n1_first:n1_last,n1_first:n1_last)
       S(n2_first:n2_last,n2_first:n2_last) = S(n1_first:n1_last,n1_first:n1_last)
    end do

    !print *,'H_dftb_orth after replacement'
    !do i=1,N
    !   write(*,*) H(i,1:N)
    !end do

    !print *,'S_dftb_orth after replacement'
    !do i=1,N
    !   write(*,*) S(i,1:N)
    !end do

    !Save H_dftb_orth.mtr to file
    !open(12,file='H_dftb_orth.mtr',action="write")
    !do i=1,N
    !   write(12,*) H(i,1:N)*HAR
    !end do
    !close(12)

    !Save S_dftb_orth.mtr to file
    !open(12,file='S_dftb_orth.mtr',action="write")
    !do i=1,N
    !   write(12,*) S(i,1:N)
    !end do
    !close(12)

    write(*,"(' Löwdin orthogonalization for device only is done! ')")
    
  end subroutine orthogonalization_dev

  !-----------------------------------------------------------------------------
  ! ALEX: this routine has been commented out for the moment. 
  !       Reading of DFTB Hamiltonian should be moved inside DFTB
  !       Reading of model Hamiltonian should be a separate libNEGF driver
  !
  !DAR begin - negf_current_nogeom
  !-----------------------------------------------------------------------------
  !subroutine negf_current_nogeom(mpicomm, tundos)   
  !  
  !  type(z_CSR) :: HH, SS
  !  type(mpifx_comm), intent(in) :: mpicomm
  !  real(dp), dimension(:,:), pointer :: tunn
  !  real(dp), dimension(:,:), pointer :: ledos
  !  real(dp), dimension(:), pointer :: currents
  !  type(unit) :: unitOfEnergy        ! Set the units of H
  !  type(unit) :: unitOfCurrent       ! Set desired units for Jel
  !
  !  integer :: i,j,k,l,m,n,NumStates,icont
  !  real(8), dimension(:), allocatable :: coupling   
  !  real(dp), allocatable :: H(:,:),S(:,:)
  !
  !  unitOfEnergy%name = "H"
  !  unitOfCurrent%name = "A"
  !  
  !  if (id0) then
  !  write(*,*)
  !  write(*,'(80("="))')
  !  write(*,*) '                            COMPUTATION OF TRANSPORT         '
  !  write(*,'(80("="))') 
  !  write(*,*)
  !  write(*,*)'Transport is started'
  !  write(*,"(' Number of States = ',I0)")negf%NumStates
  !  endif
  !  
  !  negf%kbT=0.00001_dp
  !  negf%g_spin=2
  !
  !  if(negf%tReadDFTB) call ReadDFTB
  !    
  !  if(negf%tModel) call ReadModel
  !
  !  if(negf%tOrthonormal) call Orthogonalization
  !
  !  if(negf%tOrthonormalDevice) call Orthogonalization_dev
  !
  !  if(negf%tElastic.and.(negf%tReadDFTB.or.negf%tModel.or.negf%tOrthonormal.or.negf%tOrthonormalDevice)) &
  !       call MakeHHSS(HH,SS)
  !
  !  if(negf%tManyBody) call MakeHS_dev
  !
  !  if(negf%tRead_negf_in) call ReadLibNEGF
  !
  !  call pass_HS(negf,HH,SS)
  ! 
  !
  !  if(negf%tWrite_ldos) call WriteLDOS
  !  
  !  if (id0.and.negf%verbose.gt.30) then
  !    write(*,*)
  !    write(*,'(80("="))')
  !    write(*,*) '                          LibNEGF: Current calculation'  
  !    write(*,'(80("="))') 
  !    write(*,*)
  !  endif
  !
  !  call compute_current(negf)    
  !
  !  call associate_current(negf, currents)
  !  call associate_ldos(negf, ledos)
  !  call associate_transmission(negf, tunn)
  !
  !  if (id0.and.negf%verbose.gt.30) then
  !    write(*,*)
  !    write(*,'(80("="))')
  !    write(*,*) '                           LibNEGF: Current finished'   
  !    write(*,'(80("="))') 
  !    write(*,*)
  !  endif      
  !
  !  !-------------------------------------------------------------------------   
  !  !WriteSelfEnergy / WriteSurfaceGF
  !  !-------------------------------------------------------------------------
  !  
  !  do icont=1,negf%str%num_conts
  !    if(negf%tranas%cont(icont)%tWriteSelfEnergy) &
  !       call mpifx_allreduceip(mpicomm, negf%tranas%cont(icont)%SelfEnergy, MPI_SUM)
  !  end do
  !     
  !  if (id0) then
  !    do icont=1,negf%str%num_conts
  !      if(negf%tranas%cont(icont)%tWriteSelfEnergy) then
  !        open(14,form="unformatted",file=trim(negf%tranas%cont(icont)%name)//'-SelfEnergy.mgf' &
  !             ,action="write")
  !        do i = 1, size(negf%en_grid)
  !           write(14)real(negf%en_grid(i)%Ec),negf%tranas%cont(icont)%SelfEnergy(:,:,i)
  !        end do
  !        close(14)
  !        write(*,"('    The retarded contact self-energy is written into the file ',A)") &
  !              trim(negf%tranas%cont(icont)%name)//'-SelfEnergy.mgf'
  !      end if
  !    end do
  !  end if
  !
  !  do icont=1,negf%str%num_conts
  !     if(negf%tranas%cont(icont)%tWriteSurfaceGF) &
  !          call mpifx_allreduceip(mpicomm, negf%tranas%cont(icont)%SurfaceGF, MPI_SUM)
  !  end do
  !
  !  if (id0) then
  !    do icont=1,negf%str%num_conts
  !      if (negf%tranas%cont(icont)%tWriteSurfaceGF) then
  !        open(14,form="unformatted",file=trim(negf%tranas%cont(icont)%name)//'-SurfaceGF.mgf' &
  !               ,action="write")
  !        do i = 1, size(negf%en_grid)
  !          write(14)real(negf%en_grid(i)%Ec),negf%tranas%cont(icont)%SurfaceGF(:,:,i)
  !        end do
  !        close(14)
  !        write(*,"('    The retarded contact self-energy is written into the file ',A)") &
  !                trim(negf%tranas%cont(icont)%name)//'-SurfaceGF.mgf'
  !      end if
  !    end do
  !  end if
  !
  !  !-------------------------------------------------------------------------
  !
  !  call mpifx_allreduceip(mpicomm, currents, MPI_SUM)
  !
  !  currents = currents * convertCurrent(unitOfEnergy, unitOfCurrent) 
  !  
  !  if (id0) then 
  !    do i=1, size(currents)
  !      write(*,'(1x,a,i3,i3,a,ES14.5,a,a)') &
  !           & ' contacts: ',negf%ni(i),negf%nf(i), &
  !           & ' current: ', currents(i),' ',unitOfCurrent%name
  !    enddo
  !  endif
  !
  !  call mpifx_allreduceip(mpicomm, tunn, MPI_SUM)
  !  if (id0 .and. tundos%writeTunn) then  
  !     open(65000,file='tunneling.dat')
  !     do i=1,size(tunn,1)
  !        write(65000,'(f20.8)',ADVANCE='NO') (negf%Emin+i*negf%Estep)*HAR
  !        do j=1,size(tunn,2)
  !           write(65000,'(f20.8)',ADVANCE='NO') tunn(i,j)
  !        enddo
  !        write(65000,*)
  !     enddo
  !     close(65000)
  !  endif
  !
  !  if(negf%tZeroCurrent) then                                                
  !     call mpifx_allreduceip(mpicomm, negf%tunn_mat_bp, MPI_SUM)        
  !  if (id0 .and. tundos%writeTunn) then  
  !     open(65000,file='tunneling_bp.dat')
  !     do i=1,size(negf%tunn_mat_bp,1)
  !        write(65000,'(f20.8)',ADVANCE='NO') (negf%Emin+i*negf%Estep)*HAR
  !        do j=1,size(negf%tunn_mat_bp,2)
  !           write(65000,'(f20.8)',ADVANCE='NO') negf%tunn_mat_bp(i,j)
  !        enddo
  !        write(65000,*)
  !     enddo
  !     close(65000)
  !  endif 
  !  endif
  !
  !  if ((.not.allocated(negf%inter)).and.(.not.negf%tDephasingBP)) then
  !  call mpifx_allreduceip(mpicomm, ledos, MPI_SUM)
  !  if (id0 .and. tundos%writeLDOS) then
  !     open(65000,file='localDOS.dat')
  !     do i=1,size(ledos,1)
  !        write(65000,'(f20.8)',ADVANCE='NO') (negf%Emin+i*negf%Estep)*HAR
  !        do j=1,size(ledos,2)
  !           write(65000,'(f20.8)',ADVANCE='NO') ledos(i,j)
  !        enddo
  !        write(65000,*)
  !     enddo
  !     close(65000)
  !  end if
  !  end if
  !
  !  if (negf%tWrite_ldos) call mpifx_allreduceip(mpicomm, ledos, MPI_SUM)
  !  if (id0 .and. negf%tWrite_ldos) then
  !     open(65000,file='localDOS.dat')
  !     do i=1,size(ledos,1)
  !        write(65000,'(f20.8)',ADVANCE='NO') (negf%Emin+i*negf%Estep)*HAR
  !        do j=1,size(ledos,2)
  !           write(65000,'(f20.8)',ADVANCE='NO') ledos(i,j)
  !        enddo
  !        write(65000,*)
  !     enddo
  !     close(65000)
  !  end if    
  ! 
  !  if (id0) print*,'calculation of current done'                                     
  !      
  !end subroutine negf_current_nogeom
  !-----------------------------------------------------------------------------
  !DAR end
  !-----------------------------------------------------------------------------
  !DAR end
  !-----------------------------------------------------------------------------
    
end module negf_int
       
