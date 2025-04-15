!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

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
#:include "common.fypp"

module phonons_libnegfint
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_file, only : closeFile, openFile, TFileDescr
  use dftbp_common_globalenv, only : tIoProc
  use dftbp_extlibs_negf, only : associate_ldos, associate_lead_currents, associate_transmission,&
      & COMP_SGF, COMPSAVE_SGF, compute_phonon_current, convertHeatConductance, convertHeatCurrent,&
      & create, create_scratch, csr2dns, DELTA_MINGO, DELTA_SQ, DELTA_W, destroy, destroy_matrices,&
      & destroy_negf, dns2csr, get_params, getel, init_contacts, init_ldos, init_negf,&
      & init_structure, kb, lnParams, nzdrop, pass_DM, pass_hs, printcsr, READ_SGF,&
      & set_bp_dephasing, set_drop, set_elph_block_dephasing, set_elph_dephasing,&
      & set_elph_s_dephasing, set_ldos_indexes, set_params, set_readoldDMsgf, set_scratch,&
      & set_tun_indexes, thermal_conductance, Tnegf, units, writememinfo, writepeakinfo, z_CSR,&
      & z_DNS
  use dftbp_io_message, only : error
  use dftbp_transport_matconv, only : destruct
  use dftbp_transport_negfvars, only : TNEGFTunDos, TTranspar
  use dftbp_type_commontypes, only : TOrbitals
#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_comm, mpifx_reduceip
  use dftbp_extlibs_negf, only : negf_mpi_init
#:endif
  use phonons_initphonons, only : modeEnum, TempMax, TempMin, TempStep
  implicit none
  private

  Type(Tnegf), target, public :: negf
  ! Workaround: ifort 17, ifort 16
  ! Passing negf for pointer dummy arguments fails despite target attribute, so pointer is needed
  type(TNegf), pointer :: pNegf

  public :: negf_init, negf_destroy

  !> Initialize tunneling projection on specific modes
  public :: init_tun_proj

  ! This initializes the partitioned structure
  public :: negf_init_str

  ! initialize phonon-phonon interactions (call after negf_init_str)
  !public :: negf_init_phph

  ! INTERFACE SUBROUTINE TO DRIVE PHONON CALCULATIONS:
  public :: calc_phonon_current

  ! direct calls to compute phonon current
  public :: negf_phonon_current

  ! interface csr matrices. The pattern structure of csrHam is defined by negf_init_csr.

  type(z_CSR), target :: csrHam
  type(Z_CSR), pointer :: pCsrHam => null()


contains


  !> Init gDFTB environment and variables
  subroutine negf_init(env, transpar, tundos, initinfo)

    !> Environment settings, suplying the global comm world
    type(TEnvironment), intent(in) :: env

    !> Parameters for the transport calculation
    Type(TTranspar), intent(in) :: transpar

    !> Parameters for tuneling and density of states evaluation
    Type(TNEGFTunDos), intent(in) :: tundos

    !> Initialization flag
    logical, intent(out) :: initinfo


    ! local variables
    integer :: i, l, ncont, nc_vec(1), nldos
    integer, dimension(:), allocatable :: sizes
    ! string needed to hold processor name
    character(:), allocatable :: hostname
    type(lnParams) :: parms

    initinfo = .true.

    pNegf=>negf

  #:if WITH_MPI
    call negf_mpi_init(env%mpi%globalComm, tIOproc)
  #:endif

    if (transpar%defined) then
      ncont = transpar%ncont
    else
      ncont = 0
    endif

    ! Set defaults and fill up the parameter structure with them
    call init_negf(negf)
    call init_contacts(negf, ncont)
    call get_params(negf, parms)
    call set_scratch(negf, ".")

    !                        SETTING CONTACT TEMPERATURES
    ! If no contacts => no transport,
    if (transpar%defined) then

      do i = 1, ncont
        if (transpar%contacts(i)%kbT .ge. 0.0_dp) then
          parms%kbT_t(i) = transpar%contacts(i)%kbT
          parms%kbT_dm(i) = transpar%contacts(i)%kbT
        end if
      enddo

      ! set parameters for wide band approximations
      do i=1, ncont
         parms%FictCont(i) = transpar%contacts(i)%wideBand
         parms%contact_DOS(i) = transpar%contacts(i)%wideBandDOS
       enddo

    end if

    ! This parameter is used to set the averall drop threshold in libnegf
    ! It affects especially transmission that is not accurate more than
    ! this value.
    call set_drop(1.d-20)

    !                    SETTING TRANSMISSION PARAMETERS
    if (tundos%defined) then
      parms%verbose = tundos%verbose
      select case (tundos%deltaModel)
      case(DELTA_SQ)
        parms%deltaModel = DELTA_SQ
      case(DELTA_W)
        parms%deltaModel = DELTA_W
      case(DELTA_MINGO)
        parms%deltaModel = DELTA_MINGO
      case default
        call error('Internal error deltaModel not properly set')
      end select
      parms%Wmax = tundos%Wmax
      parms%delta = tundos%delta      ! delta for G.F.
      parms%dos_delta = tundos%broadeningDelta

      l = size(tundos%ni)
      parms%ni(1:l) = tundos%ni(1:l)
      parms%nf(1:l) = tundos%nf(1:l)

      parms%Emin =  tundos%Emin
      parms%Emax =  tundos%Emax
      parms%Estep = tundos%Estep
      ! set SGF reload to compute
      parms%readOldDM_SGFs = COMP_SGF
      parms%readOldT_SGFs = COMP_SGF
    endif

    ! Energy conversion only affects output units.
    ! The library writes energies as (E * negf%eneconv)
    parms%eneconv = 1.d0

    parms%isSid = .true.

    if (allocated(sizes)) then
      deallocate(sizes)
    end if

    call set_params(negf,parms)

    ! set indeces for projected DOS
    if (tundos%defined) then
      nldos = size(tundos%dosOrbitals)
      call init_ldos(negf, nldos)
      do i = 1, nldos
         call set_ldos_indexes(negf, i, tundos%dosOrbitals(i)%data)
      end do
    end if

  end subroutine negf_init


  subroutine init_tun_proj(selTypeModes, nAtoms)
    integer, intent(in) :: selTypeModes
    integer, intent(in) :: nAtoms

    integer, dimension(:), allocatable :: indx
    integer :: i, l
    ! set indeces for projecter transmission
    ! This assumes derivatives are ordered x,y,z
    ! That transport is along z
    ! That 2d structures are on x-z
    select case(selTypeModes)
    case(modeEnum%ALLMODES)
      allocate(indx(3*nAtoms))
      do i = 1, 3*nAtoms
        indx(i) = i
      end do
    case(modeEnum%XX)
      allocate(indx(nAtoms))
      l = 1
      do i = 1, 3*nAtoms, 3
        indx(l) = i
        l = l + 1
      end do
    case(modeEnum%YY,modeEnum%OUTOFPLANE)
      allocate(indx(nAtoms))
      l = 1
      do i = 2, 3*nAtoms, 3
        indx(l) = i
        l = l + 1
      end do
    case(modeEnum%ZZ,modeEnum%LONGITUDINAL)
      allocate(indx(nAtoms))
      l = 1
      do i = 3, 3*nAtoms, 3
        indx(l) = i
        l = l + 1
      end do
    case(modeEnum%TRANSVERSE)
      allocate(indx(2*nAtoms))
      l = 1
      do i = 1, 3*nAtoms, 3
        indx(l) = i     ! X
        l = l + 1
        indx(l) = i+1   ! Y
        l = l + 1
      end do
    case(modeEnum%INPLANE)
      allocate(indx(2*nAtoms))
      l = 1
      do i = 1, 3*nAtoms, 3
        indx(l) = i     ! X
        l = l + 1
        indx(l) = i+2   ! Z
        l = l + 1
      end do
    end select

    call set_tun_indexes(negf, indx)

  end subroutine init_tun_proj


  subroutine negf_destroy(output)

    !> Output for write processes
    integer, intent(in) :: output

    write(output, *)
    write(output, *) 'Release NEGF memory:'
    call destruct(csrHam)
    call destroy_negf(negf)
    call writePeakInfo(output)
    call writeMemInfo(output)

  end subroutine negf_destroy


  subroutine negf_init_str(output, nAtoms, transpar, iNeigh, nNeigh, img2CentCell)

    !> Output for write processes
    integer, intent(in) :: output

    !> Number of atoms
    integer, intent(in) :: nAtoms

    !> Transport calculation parameters
    Type(TTranspar), intent(in) :: transpar

    !> Neighbours of each atom
    Integer, intent(in) :: iNeigh(0:,:)

    !> Number of neighbours for each atom
    Integer, intent(in) :: nNeigh(:)

    !> Mapping from image atoms to central cell
    Integer, intent(in) :: img2CentCell(:)

    Integer, allocatable :: PL_end(:), cont_end(:), surf_start(:), surf_end(:), cblk(:), ind(:)
    Integer, allocatable :: atomst(:), plcont(:)
    integer, allocatable :: minv(:,:)
    Integer :: ncont, nbl, iatm1, iatm2, iatc1, iatc2
    integer :: i, m, i1, j1

    iatm1 = transpar%idxdevice(1)
    iatm2 = transpar%idxdevice(2)

    ncont = transpar%ncont
    nbl = transpar%nPLs
    if (nbl.eq.0) then
      call error('Internal ERROR: nbl = 0 ?!')
    end if

    allocate(PL_end(nbl))
    allocate(atomst(nbl+1))
    allocate(plcont(nbl))
    allocate(cblk(ncont))
    allocate(cont_end(ncont))
    allocate(surf_start(ncont))
    allocate(surf_end(ncont))
    allocate(ind(natoms+1))
    allocate(minv(nbl,ncont))

    do i = 1, nAtoms + 1
      ind(i) = 3*(i-1)
    end do

    do i = 1, ncont
       cont_end(i) = ind(transpar%contacts(i)%idxrange(2)+1)
       surf_start(i) = ind(transpar%contacts(i)%idxrange(1)) + 1
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
             write(output,*) 'Contact',j1,'interacts with more than one PL:'
             write(output,*) 'PLs:',minv(:,j1)
             call error('check cutoff value or PL size')
          end if
          do m = 1, transpar%nPLs
             if (minv(m,j1).eq.j1) cblk(j1) = m
          end do
       end do


       write(output,*)
       write(output,*) ' Structure info:'
       write(output,*) ' Number of PLs:',nbl
       write(output,*) ' PLs coupled to contacts:',cblk(1:ncont)
       write(output,*)

    end if

    call init_structure(negf, ncont, surf_start, surf_end, cont_end, nbl, PL_end, cblk)

    deallocate(PL_end)
    deallocate(plcont)
    deallocate(atomst)
    deallocate(cblk)
    deallocate(cont_end)
    deallocate(surf_end)
    deallocate(ind)
    deallocate(minv)

  end subroutine negf_init_str


  !subroutine negf_init_phph(negf, order)
  !   type(Tnegf) :: negf
  !   integer, intent(in) :: order
  !
  !   select case (order)
  !   case(3, 34)
  !     print*,'Init cubic phonon-phonon interactions'
  !     call set_phph(negf, 3, 'cubic.dat')
  !   case(4)
  !     print*,'Init quartic phonon-phonon interactions'
  !     call set_phph(negf, 4, 'quartic.dat')
  !   end select
  !
  !end subroutine negf_init_phph


  !> INTERFACE subroutine to call phonon current computation
  subroutine calc_phonon_current(env, DynMat, tunnMat, ldosMat, &
                        & currLead, conductance, twriteTunn, twriteLDOS)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> The dynamical matrix of the system
    real(dp), intent(in) :: DynMat(:,:)

    !> Matrix of tunnelling amplitudes at each energy from contacts
    real(dp), allocatable, intent(inout) :: tunnMat(:,:)

    !> Density of states for each energy and region of projection
    real(dp), allocatable, intent(inout) :: ldosMat(:,:)

    !> Current into/out of contacts
    real(dp), allocatable, intent(inout) :: currLead(:)

    !> Thermal conductance
    real(dp), allocatable, intent(inout) :: conductance(:, :)

    !> Should tunneling data be written
    logical, intent(in) :: tWriteTunn

    !> Should DOS data be written
    logical, intent(in) :: tWriteLDOS

    ! locals
    real(dp), allocatable :: tunnSKRes(:,:,:), ldosSKRes(:,:,:)
    real(dp), pointer    :: tunnPMat(:,:)=>null()
    real(dp), pointer    :: ldosPMat(:,:)=>null()
    real(dp), pointer    :: currPVec(:)=>null()
    type(TFileDescr) :: fd
    integer :: ii, jj, iK, nK, err, nnz, ntemp
    real(dp), allocatable :: kPoints(:,:), kWeights(:)
    type(z_DNS) :: zDynMat
    real(dp) :: cutoff,TT1,emin,emax,estep, kappa
    type(units) :: HessianUnits, HeatCurrUnits, HeatCondUnits
    type(lnParams) :: params
    character(15) :: filename

    ! Sets the Gamma point calculation
    nK = 1
    allocate(kPoints(3, nK))
    allocate(kWeights(nK))
    kPoints = 0.0_dp
    kWeights = 1.0_dp
    cutoff = 1d-30
    HessianUnits%name = "H"
    HeatCurrUnits%name = "W"
    HeatCondUnits%name = "W/K"
    pCsrHam => csrHam

    call get_params(negf, params)

    do iK = 1, nK

      call create(zDynMat, size(DynMat,1), size(DynMat,2) )
      zDynMat%val = DynMat

      nnz = nzdrop(zDynMat,cutoff)
      call create(csrHam, zDynMat%nrow, zDynMat%ncol, nnz)

      call dns2csr(zDynMat, csrHam)

      call destroy(zDynMat)

      call negf_phonon_current(pCsrHam, iK, kWeights(iK), &
            tunnPMat, ldosPMat, currPVec)

      if(.not.allocated(currLead)) then
         allocate(currLead(size(currPVec)), stat=err)
         if (err/=0) then
            call error('Allocation error (currTot)')
         end if
         currLead = 0.0_dp
       endif
       currLead(:) = currLead + currPVec

     #:if WITH_MPI
       call add_partial_results(env%mpi%groupComm, tunnPMat, tunnMat, tunnSKRes, iK, nK)
       call add_partial_results(env%mpi%groupComm, ldosPMat, ldosMat, ldosSKRes, iK, nK)
     #:else
       call add_partial_results(tunnPMat, tunnMat, tunnSKRes, iK, nK)
       call add_partial_results(ldosPMat, ldosMat, ldosSKRes, iK, nK)
     #:endif

    end do

    ! MPI Reduce k-dependent stuff
  #:if WITH_MPI
    call mpifx_reduceip(env%mpi%groupComm, currLead, MPI_SUM)
    call mpifx_reduceip(env%mpi%interGroupComm, currLead, MPI_SUM)
    call add_k_results(env%mpi%interGroupComm, tunnMat, tunnSKRes )
    call add_k_results(env%mpi%interGroupComm, ldosMat, ldosSKRes )
  #:endif

    ! converts from internal atomic units into W
    currLead = currLead * convertHeatCurrent(HessianUnits, HeatCurrUnits)


    if (tIOProc) then
      do ii= 1, size(currLead)
        write(*,'(1x,a,i3,i3,a,ES14.5,a,a)') &
             & ' contacts: ', params%ni(ii), params%nf(ii), &
             & ' current: ', currLead(ii),' ',HeatCurrUnits%name
      enddo
    endif

    if (allocated(tunnMat)) then
      ntemp=nint((TempMax-TempMin)/TempStep)
      allocate(conductance(ntemp,size(tunnMat,2)+1))
      emin = negf%Emin*negf%eneconv
      emax = negf%Emax*negf%eneconv
      estep = negf%Estep*negf%eneconv
      do ii = 1, size(tunnMat,2)
        do jj = 1, ntemp
          TT1 = TempMin + TempStep*(jj-1)
          kappa = thermal_conductance(tunnMat(:,ii),TT1,emin,emax,estep)
          kappa = kappa * convertHeatConductance(HessianUnits,HeatCondUnits)
          conductance(jj,1) = TT1/kb
          conductance(jj,ii+1) = kappa
        end do
      end do
      ! Write Total tunneling on a separate file (optional)
      if (tIOProc .and. twriteTunn) then
        filename = 'transmission'
        call write_file(negf, tunnMat, tunnSKRes, filename, kpoints, kWeights)

        call openFile(fd, "conductance.dat", mode="w")
        do ii = 1, size(tunnMat,2)
          write(fd%unit, *) '# T [K]', 'Thermal Conductance [W/K]'
          do jj = 1, ntemp
            write(fd%unit, *) conductance(jj,1), conductance(jj,ii+1)
          end do
        end do
        call closeFile(fd)
      end if

    else
      allocate(tunnMat(0,0))
    endif

    if (allocated(tunnSKRes)) then
      deallocate(tunnSKRes)
    end if

    if (allocated(ldosMat)) then
      ! Multiply density by 2w
      do ii = 1, size(ldosMat,1)
        ldosMat(ii,:) = ldosMat(ii,:) * 2.d0*(negf%emin + negf%estep*(ii-1))
      end do
      ! Write Total localDOS on a separate file (optional)
      if (tIOProc .and. twriteLDOS) then
        filename = 'localdos'
        call write_file(negf, ldosMat, ldosSKRes, filename, kpoints, kWeights)
      end if
    else
      allocate(ldosMat(0,0))
    endif
    if (allocated(ldosSKRes)) then
      deallocate(ldosSKRes)
    end if

    deallocate(kPoints)
    deallocate(kWeights)

  end subroutine calc_phonon_current



  subroutine negf_phonon_current(HH, qpoint, wght, tunn, ledos, currents)

    !> Hessian
    type(z_CSR), pointer, intent(in) :: HH

    !> Phonon kpoint (kp index)
    integer, intent(in) :: qpoint

    !> Phonon k-weight (kp weight)
    real(dp), intent(in) :: wght

    !> Tunneling
    real(dp), dimension(:,:), pointer :: tunn

    !> Local or projected dos
    real(dp), dimension(:,:), pointer :: ledos

    !> Heat currents
    real(dp), dimension(:), pointer :: currents

    type(lnParams) :: params

    call get_params(negf, params)

    params%ikpoint = qpoint
    params%kwght = wght

    call set_params(negf, params)

    call pass_HS(negf,HH)

    call compute_phonon_current(negf)

    call associate_ldos(pNEgf, ledos)
    call associate_transmission(pNegf, tunn)

    call associate_lead_currents(pNegf, currents)
    if (.not.associated(currents)) then
      call error('Internal error: currVec not associated')
    end if

  end subroutine negf_phonon_current

  subroutine printH(fu, H)
    integer, intent(in) :: fu
    type(z_CSR), intent(in) :: H

    type(z_DNS) :: tmp
    integer :: ii, jj
    real(dp) :: maxv

    call create(tmp, H%nrow, H%ncol)

    call csr2dns(H,tmp)

    maxv = maxval(abs(tmp%val))
    write(fu, *) 'Normalized Dynamical Matrix:'

    do ii = 1, tmp%nrow, 3
       do jj = 1, tmp%ncol, 3
          write(fu,'(F8.4)',advance='no') real(tmp%val(ii,jj))/maxv
       end do
       write(fu,*)
    end do

    call destroy(tmp)

  end subroutine printH


  !> Utility to allocate and sum partial results from different channels
#:if WITH_MPI
  subroutine add_partial_results(mpicomm, pMat, matTot, matSKRes, iK, nK)

    !> MPI communicator
    type(mpifx_comm), intent(in) :: mpicomm
#:else
  subroutine add_partial_results(pMat, matTot, matSKRes, iK, nK)
#:endif

    !> Pointer to matrix of data
    real(dp), intent(in), pointer :: pMat(:,:)

    !> Sum total
    real(dp), allocatable, intent(inout) :: matTot(:,:)

    !> Sum resolved by k
    real(dp), allocatable, intent(inout)  :: matSKRes(:,:,:)

    !> Particular k-point
    integer, intent(in) :: iK

    !> Number of k-points
    integer, intent(in) :: nK

    #:if WITH_MPI
    real(dp), allocatable :: tmpMat(:,:)
    #:endif

    integer :: err

    if (associated(pMat)) then
    #:if WITH_MPI
      allocate(tmpMat(size(pMat,dim=1), size(pMat,dim=2)), stat=err)

      if (err /= 0) then
        call error('Allocation error (tmpMat)')
      end if

      tmpMat(:,:) = pMat
      call mpifx_reduceip(mpicomm, tmpMat, MPI_SUM)
    #:endif
      if(.not.allocated(matTot)) then
        allocate(matTot(size(pMat,dim=1), size(pMat,dim=2)), stat=err)

        if (err /= 0) then
          call error('Allocation error (tunnTot)')
        end if

        matTot(:,:) = 0.0_dp
      end if
    #:if WITH_MPI
      matTot(:,:) = matTot + tmpMat
    #:else
      matTot(:,:) = matTot + pMat
    #:endif

      if (nK > 1) then
        if (.not.allocated(matSKRes)) then
          allocate(matSKRes(size(pMat,dim=1), size(pMat,dim=2), nK), stat=err)

          if (err/=0) then
            call error('Allocation error (tunnSKRes)')
          end if

          matSKRes(:,:,:) = 0.0_dp
        endif
      #:if WITH_MPI
        matSKRes(:,:,iK) = tmpMat
      #:else
        matSKRes(:,:,iK) = pMat
      #:endif
      end if

    #:if WITH_MPI
      deallocate(tmpMat)
    #:endif

    end if

  end subroutine add_partial_results



#:if WITH_MPI

  !> Utility to sum up partial results over K communicator
  subroutine add_k_results(kcomm, mat, matSKRes)

    !> MPI communicator
    type(mpifx_comm), intent(in) :: kcomm

    !> Sum total
    real(dp), allocatable, intent(inout) :: mat(:,:)

    !> Sum resolved by k
    real(dp), allocatable, intent(inout)  :: matSKRes(:,:,:)

    if (allocated(mat)) then
      call mpifx_reduceip(kcomm, mat, MPI_SUM)
    endif

    if (allocated(matSKRes)) then
      call mpifx_reduceip(kcomm, matSKRes, MPI_SUM)
    endif

  end subroutine add_k_results

#:endif


  !----------------------------------------------------------------------------
  ! init_csr: is needed only if H and S are stored in dftb+ format
  !----------------------------------------------------------------------------
  !subroutine negf_init_csr(iAtomStart, iNeighbor, nNeighbor, img2CentCell, orb)

  !  !> Start of orbitals for each atom
  !  integer, intent(in) :: iAtomStart(:)

  !  !> neighbours of each atom
  !  integer, intent(in) :: iNeighbor(0:,:)

  !  !> number of neighbours for each atom
  !  integer, intent(in) :: nNeighbor(:)

  !  !> mapping from image atoms to central cell
  !  integer, intent(in) :: img2CentCell(:)

  !  !> atomic orbital information
  !  type(TOrbitals), intent(in) :: orb

  !  pCsrHam => csrHam
  !  pCsrOver => csrOver
  !  if (allocated(csrHam%nzval)) then
  !    call destroy(csrHam)
  !  end if
  !  call init(csrHam, iAtomStart, iNeighbor, nNeighbor, img2CentCell, orb)
  !  if (allocated(csrOver%nzval)) then
  !    call destroy(csrOver)
  !  end if
  !  call init(csrOver, csrHam)

  !end subroutine negf_init_csr


  !> Write the transmission or local dos to a file
  subroutine write_file(negf, pTot, pKRes, filename, kpoints, kWeights)

    !> Contains input data, runtime quantities and output data
    type(TNegf) :: negf

    !> Total data to be written
    real(dp), intent(in) :: pTot(:,:)

    !> The k-point resolved data, if allocated
    real(dp), allocatable, intent(in) :: pKRes(:,:,:)

    !> File to print out to
    character(*), intent(in) :: filename

    !> The k-points for the system
    real(dp), intent(in) :: kPoints(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeights(:)

    type(TFileDescr) :: fd
    integer :: ii, jj, iK, nK

    nK = size(kPoints,2)
    call openFile(fd, trim(filename) // '.dat', mode="w")
    if (trim(filename).eq.'transmission') then
      write(fd%unit, *)  '# Energy [H]', '  Transmission'
    else
      write(fd%unit, *)  '# Energy [H]', '  LDOS'
    endif
    do ii=1,size(pTot,1)
      write(fd%unit, '(es20.8)', advance='no') (negf%Emin+(ii-1)*negf%Estep)*negf%eneconv
      do jj=1,size(pTot,2)
        write(fd%unit, '(es20.8)', advance='no') pTot(ii,jj)
      enddo
      write(fd%unit, *)
    enddo
    call closeFile(fd)

    if (nK.gt.1) then
      call openFile(fd, trim(filename) // '_kpoints.dat', mode="w")
      write(fd%unit, *)  '# NKpoints = ', nK
      write(fd%unit, *)  '# Energy [eV], <k1 k2 k3 weight> '
      write(fd%unit, '(A1)', advance="no") '# '
      do iK = 1,nK
        write(fd%unit, '(es15.5, es15.5, es15.5, es15.5)', advance="no") &
            & kpoints(:,iK), kWeights(iK)
      end do
      write(fd%unit, *)
      do ii = 1, size(pKRes(:,:,1),1)
        write(fd%unit, '(f20.8)',advance="no") (negf%Emin+(ii-1)*negf%Estep)*negf%eneconv
        do jj=1,size(pKRes(:,:,1),2)
          do iK = 1,nK
            write(fd%unit, '(es20.8)',advance="no") pKRes(ii,jj, iK)
          enddo
          write(fd%unit, *)
        enddo
      enddo
      call closeFile(fd)
    end if

  end subroutine write_file


  !> DEBUG routine dumping H and S on file in Matlab format
  subroutine negf_dumpHS(output, HH,SS)

    !> Output for write processes
    integer, intent(in) :: output

    type(z_CSR), intent(in) :: HH, SS

    type(TFileDescr) :: fd

    write(output,*) 'Dumping H and S on files...'
    call openFile(fd, 'HH.dat', mode="w")
    write(fd%unit, *) '% Size =',HH%nrow, HH%ncol
    write(fd%unit, *) '% Nonzeros =',HH%nnz
    write(fd%unit, *) '% '
    write(fd%unit, *) 'zzz = ['
    call printcsr(fd%unit, HH)
    write(fd%unit, *) ']'
    call closeFile(fd)

    call openFile(fd, 'SS.dat', mode="w")
    write(fd%unit, *) '% Size =',SS%nrow, SS%ncol
    write(fd%unit, *) '% Nonzeros =',SS%nnz
    write(fd%unit, *) '% '
    write(fd%unit, *) 'zzz = ['
    call printcsr(fd%unit, SS)
    write(fd%unit, *) ']'
    call closeFile(fd)
  end subroutine negf_dumpHS


end module phonons_libnegfint
