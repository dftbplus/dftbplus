!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Interface to LIBNEGF for DFTB+
module negf_int

  use accuracy
  use constants
  use libnegf_vars
  use libnegf, only : convertcurrent, eovh, getel, lnParams, negf_mpi_init, pass_DM, Tnegf, unit
  use libnegf, only : z_CSR
  use libnegf, only : associate_lead_currents, associate_ldos, associate_transmission
  use libnegf, only : associate_current, compute_current, compute_density_dft, compute_ldos
  use libnegf, only : create, create_scratch, destroy
  use libnegf, only : destroy_matrices, destroy_negf, get_params, init_contacts, init_ldos
  use libnegf, only : init_negf, init_structure, log_deallocatep, pass_hs, set_bp_dephasing
  use libnegf, only : set_drop, set_elph_block_dephasing, set_elph_dephasing, set_elph_s_dephasing
  use libnegf, only : set_ldos_indexes, set_params, set_scratch, writememinfo, writepeakinfo
  use libnegf, only : read_negf_in, negf_partition_info, printcsr
  use mat_conv
  use sparse2dense
  use densedescr
  use commonTypes, only : TOrbitals
  use libmpifx_module
  use FormatOut
  use globalenv
  use message
  use solvertypes

  implicit none
  private

  type(TNegf), target, public :: negf
  ! Workaround: ifort 17, ifort 16
  ! Passing negf for pointer dummy arguments fails despite target attribute, so pointer is needed
  type(TNegf), pointer :: pNegf

  !> general library initializations
  public :: negf_init

  !> passing structure parameters
  public :: negf_init_str

  public :: negf_init_dephasing
  public :: negf_init_elph
  public :: negf_destroy

  !> wrapped functions passing dftb matrices. Needed for parallel
  public :: calcdensity_green

  !> Calculate the energy weighted density from the GF
  public :: calcEdensity_green

  !> Calculate the partial density of states from the GF
  public :: calcPDOS_green

  !> Calculate the total current
  public :: calc_current

  !> calculate local currents
  public :: local_currents

  !> interface csr matrices. The pattering must be predefined using negf_init_csr
  public :: negf_init_csr

  !> compressed sparse row hamiltonian
  type(z_CSR), target :: csrHam
  type(Z_CSR), pointer :: pCsrHam => null()

  !> compressed sparse row overlap
  type(z_CSR), target :: csrOver
  type(Z_CSR), pointer :: pCsrOver => null()

  !> non wrapped direct calls
  private :: negf_density, negf_current, negf_ldos

  contains

  !> Init gDFTB environment and variables
  !>
  !> Note: mpicomm should be the global commworld here
  subroutine negf_init(transpar, greendens, tundos, mpicomm, tempElec, solver)

    Type(TTranspar), intent(in) :: transpar
    Type(TNEGFGreenDensInfo), intent(in) :: greendens
    Type(TNEGFTunDos), intent(in) :: tundos
    Type(mpifx_comm), intent(in) :: mpicomm
    real(dp), intent(in) :: tempElec

    !> Which solver call is used in the main code
    integer, intent(in) :: solver

    ! local variables
    real(dp), allocatable :: pot(:), eFermi(:)
    integer :: i, l, ncont, nc_vec(1), j, nldos
    integer, allocatable :: sizes(:)
    type(lnParams) :: params

    ! Workaround: ifort 16
    ! Pointer must be set within a subroutine. Initialization at declaration fails.
    pNegf => negf
    call negf_mpi_init(mpicomm)

    if (transpar%defined) then
      ncont = transpar%ncont
    else
      ncont = 0
    endif

    ! ------------------------------------------------------------------------------
    ! Set defaults and fill up the parameter structure with them
    call init_negf(negf)
    call init_contacts(negf, ncont)
    call set_scratch(negf, ".")

    if (tIoProc .and. transpar%defined .and. solver == solverGF) then
      call create_scratch(negf)
    end if

    call get_params(negf, params)

    ! ------------------------------------------------------------------------------
    ! This must be different for different initialisations, to be separated
    ! Higher between transport and greendens is taken, temporary
    if (tundos%defined .and. greendens%defined) then
       if (tundos%verbose.gt.greendens%verbose) then
          params%verbose = tundos%verbose
       else
          params%verbose = greendens%verbose
       endif
    else
      if (tundos%defined) then
        params%verbose = tundos%verbose
      end if
      if (greendens%defined) then
        params%verbose = greendens%verbose
      end if
    end if
    ! ------------------------------------------------------------------------------
    ! This parameter is used to set the averall drop threshold in libnegf
    ! It affects especially transmission that is not accurate more than
    ! this value.
    call set_drop(1.0e-20_dp)


    ! ------------------------------------------------------------------------------
    ! Assign spin degenracy and check consistency between different input blocks
    if (tundos%defined .and. greendens%defined) then
      if (tundos%gSpin /= greendens%gSpin) then
        call error("spin degeneracy is not consistent between different input blocks")
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
      if (tundos%defined) then
        params%kbT_t(i) = tundos%kbT(i)
      else
        params%kbT_t(i) = tempElec
      end if
    end do

    if (ncont == 0) then
      if (greendens%defined) then
        params%kbT_dm(1) = greendens%kbT(1)
      else
        params%kbT_dm(1) = tempElec
      end if
    end if

    ! Make sure low temperatures (< 10K) converted to 0.
    ! This avoid numerical problems with contour integration
    do i = 1, size(params%kbT_dm)
      if (params%kbT_dm(i) < 3.0e-5_dp) then
        params%kbT_dm(i) = 0.0_dp
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

        write(stdOut,*) '(negf_init) CONTACT INFO #',i

        if (params%FictCont(i)) then
          write(stdOut,*) 'FICTITIOUS CONTACT '
          write(stdOut,*) 'DOS: ', params%contact_DOS(i)
        end if
        write(stdOut,*) 'Temperature (DM): ', params%kbT_dm(i)
        write(stdOut,*) 'Temperature (Current): ', params%kbT_t(i)
        write(stdOut,*) 'Potential (with built-in): ', pot(i)
        write(stdOut,*) 'eFermi: ', eFermi(i)
        write(stdOut,*)

      enddo

      ! Define electrochemical potentials
      params%mu(1:ncont) = eFermi(1:ncont) - pot(1:ncont)
      write(stdOut,*) 'Electro-chemical potentials: ', params%mu(1:ncont)
      write(stdOut,*)
      deallocate(pot)

    else
      params%mu(1) = greendens%oneFermi(1)
    end if

    ! ------------------------------------------------------------------------------
    !                  SETTING COUNTOUR INTEGRATION PARAMETERS
    ! ------------------------------------------------------------------------------
    if (greendens%defined) then
      params%Ec = greendens%enLow           ! lowest energy
      params%Np_n(1:2) = greendens%nP(1:2)  ! contour npoints
      params%n_kt = greendens%nkt           ! n*kT for Fermi

      ! Real-axis points.
      ! Override to 0 if bias is 0.0
      params%Np_real = 0
      if (ncont > 0) then
        if (any(abs(params%mu(2:ncont)-params%mu(1)) > 1.0e-10_dp)) then
           params%Np_real = greendens%nP(3)  ! real axis points
        end if
      end if

      !Read G.F. from very first iter
      if (greendens%readSGF .and. .not.greendens%saveSGF) then
        params%readOldSGF=0
      end if
      !compute G.F. at every iteration
      if (.not.greendens%readSGF .and. .not.greendens%saveSGF) then
        params%readOldSGF=1
      end if
      !Default Write on first iter
      if (.not.greendens%readSGF .and. greendens%saveSGF) then
        params%readOldSGF=2
      end if

      if(any(params%kbT_dm > 0) .and. greendens%nPoles == 0) then
        call error("Number of Poles = 0 but T > 0")
      else
         params%n_poles = greendens%nPoles
      end if
      if(all(params%kbT_dm.eq.0)) then
        params%n_poles = 0
      end if

      write(stdOut,*) 'Density Matrix Parameters'
      if (.not.transpar%defined) then
        write(stdOut,*) 'Temperature (DM): ', params%kbT_dm(1)
        write(stdOut,*) 'eFermi: ', params%mu(1)
      end if
      write(stdOut,*) 'Contour Points: ', params%Np_n(1:2)
      write(stdOut,*) 'Number of poles: ', params%N_poles
      write(stdOut,*) 'Real-axis points: ', params%Np_real(1)
      write(stdOut,*)

    end if

    ! ------------------------------------------------------------------------------
    ! Setting the delta: priority on Green Solver, if present
    ! dos_delta is used by libnegf to smoothen T(E) and DOS(E)
    ! and is currently set in tunneling
    if (tundos%defined) then
      params%dos_delta = tundos%broadeningDelta
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

      ! setting of intervals and indices for projected DOS
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
    params%eneconv = Hartree__eV

    if (allocated(sizes)) then
      deallocate(sizes)
    end if

    call set_params(negf,params)

    !--------------------------------------------------------------------------
    ! DAR begin - negf_init - TransPar to negf
    !--------------------------------------------------------------------------
    if (transpar%defined) then
      !negf%tNoGeometry = transpar%tNoGeometry
      negf%tOrthonormal = transpar%tOrthonormal
      negf%tOrthonormalDevice = transpar%tOrthonormalDevice
      negf%NumStates = transpar%NumStates
      negf%tManyBody = transpar%tManyBody
      negf%tElastic = transpar%tElastic
      negf%tZeroCurrent = transpar%tZeroCurrent
      negf%MaxIter = transpar%MaxIter
      negf%trans%out%tWriteDOS = transpar%tWriteDOS
      negf%tWrite_ldos = transpar%tWrite_ldos
      negf%tWrite_negf_params = transpar%tWrite_negf_params
      negf%trans%out%tDOSwithS = transpar%tDOSwithS
      negf%cont(:)%name = transpar%contacts(:)%name
      negf%cont(:)%tWriteSelfEnergy = transpar%contacts(:)%tWriteSelfEnergy
      negf%cont(:)%tReadSelfEnergy = transpar%contacts(:)%tReadSelfEnergy
      negf%cont(:)%tWriteSurfaceGF = transpar%contacts(:)%tWriteSurfaceGF
      negf%cont(:)%tReadSurfaceGF = transpar%contacts(:)%tReadSurfaceGF
    end if

    ! Defined outside transpar%defined ... HAS TO BE FIXED
    negf%tDephasingVE = transpar%tDephasingVE
    negf%tDephasingBP = transpar%tDephasingBP


    if((.not.negf%tElastic).and.(.not.negf%tManyBody)) then
         write(stdOut, *)'Current is not calculated!'
         call error('Choose "Elastic = Yes" or "ManyBody = Yes"!')
    end if


  end subroutine negf_init


  !> Initialise dephasing effects
  subroutine negf_init_dephasing(tundos)

    !> density of states in tunnel region
    Type(TNEGFTunDos), intent(in) :: tundos

    if(negf%tDephasingVE) then
      call negf_init_elph(tundos%elph)
    end if

    if(negf%tDephasingBP) then
      call negf_init_bp(tundos%bp)
    end if

  end subroutine negf_init_dephasing


  !> Initialise electron-phonon coupling model
  subroutine negf_init_elph(elph)

    !> el-ph coupling structure
    type(TElPh), intent(in) :: elph

    write(stdOut,*)
    select case(elph%model)
    case(1)
      write(stdOut,*) 'Setting local fully diagonal (FD) elastic dephasing model'
      call set_elph_dephasing(negf, elph%coupling, elph%scba_niter)
    case(2)
      write(stdOut,*) 'Setting local block diagonal (BD) elastic dephasing model'
      call set_elph_block_dephasing(negf, elph%coupling, elph%orbsperatm, elph%scba_niter)
    case(3)
      write(stdOut,*) 'Setting overlap mask (OM) block diagonal elastic dephasing model'
      call set_elph_s_dephasing(negf, elph%coupling, elph%orbsperatm, elph%scba_niter)
    case default
      call error("This electron-phonon model is not supported")
    end select

  end subroutine negf_init_elph


  !> Initialise Buttiker Probe dephasing
  subroutine negf_init_bp(elph)

    !> el-ph coupling structure
    type(TElPh), intent(in) :: elph

    write(stdOut,*)
    select case(elph%model)
    case(1)
      write(stdOut,*) 'Setting local fully diagonal (FD) BP dephasing model'
      !write(stdOut,*) 'coupling=',elph%coupling
      call set_bp_dephasing(negf, elph%coupling)
    case(2)
      write(stdOut,*) 'Setting local block diagonal (BD) BP dephasing model'
      call error('NOT IMPLEMENTED! INTERRUPTED!')
    case(3)
      write(stdOut,*) 'Setting overlap mask (OM) block diagonal BP dephasing model'
      call error('NOT IMPLEMENTED! INTERRUPTED!')
    case default
      call error("BP model is not supported")
    end select

  end subroutine negf_init_bp

  !> Initialise compressed sparse row matrices
  subroutine negf_init_csr(iAtomStart, iNeighbor, nNeighbor, img2CentCell, orb)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    type(TOrbitals), intent(in) :: orb

    pCsrHam => csrHam
    pCsrOver => csrOver
    if (allocated(csrHam%nzval)) then
      call destroy(csrHam)
    end if
    call init(csrHam, iAtomStart, iNeighbor, nNeighbor, img2CentCell, orb)
    if (allocated(csrOver%nzval)) then
      call destroy(csrOver)
    end if
    call init(csrOver, csrHam)

  end subroutine negf_init_csr

  !> Destroy (module stored!) CSR matrices
  subroutine negf_destroy()

    call destruct(csrHam)
    call destruct(csrOver)
    call destroy_negf(negf)

    write(stdOut, *)
    write(stdOut, *) 'Release NEGF memory:'
    !if (tIoProc) then
    !  call writePeakInfo(6)
    !  call writeMemInfo(6)
    !end if
    call writePeakInfo(stdOut)
    call writeMemInfo(stdOut)

  end subroutine negf_destroy

  !------------------------------------------------------------------------------
  subroutine negf_init_str(denseDescr, transpar, greendens, iNeigh, nNeigh, img2CentCell)
    Type(TDenseDescr), intent(in) :: denseDescr
    Type(TTranspar), intent(in) :: transpar
    Type(TNEGFGreenDensInfo) :: greendens
    Integer, intent(in) :: nNeigh(:)
    Integer, intent(in) :: img2CentCell(:)
    Integer, intent(in) :: iNeigh(0:,:)

    Integer, allocatable :: PL_end(:), cont_end(:), surf_end(:), cblk(:), ind(:)
    Integer, allocatable :: atomst(:), plcont(:)
    integer, allocatable :: minv(:,:)
    Integer :: natoms, ncont, nbl, iatc1, iatc2, iatm2
    integer :: i, m, i1, j1, info
    integer, allocatable :: inRegion(:)

    iatm2 = transpar%idxdevice(2)
    ncont = transpar%ncont
    nbl = 0

    if (transpar%defined) then
       nbl = transpar%nPLs
    else if (greendens%defined) then
       nbl = greendens%nPLs
    endif

    if (nbl.eq.0) then
      call error('Internal ERROR: nbl = 0 ?!')
    end if

    natoms = size(denseDescr%iatomstart) - 1

    call check_pls(transpar, greendens, natoms, iNeigh, nNeigh, img2CentCell, info)

    allocate(PL_end(nbl))
    allocate(atomst(nbl+1))
    allocate(plcont(nbl))
    allocate(cblk(ncont))
    allocate(cont_end(ncont))
    allocate(surf_end(ncont))
    allocate(ind(natoms+1))
    allocate(minv(nbl,ncont))

    ind(:) = DenseDescr%iatomstart(:) - 1
    minv = 0
    cblk = 0

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

    if (transpar%defined .and. ncont.gt.0) then

      if(.not.transpar%tNoGeometry) then

       ! For each PL finds the min atom index among the atoms in each contact
       ! At the end the array minv(iPL,iCont) can have only one value != 0
       ! for each contact and this is the interacting PL
       ! NOTE: the algorithm works with the asymmetric neighbor-map of dftb+
       !       because atoms in contacts have larger indices than in the device
       do m = 1, transpar%nPLs
          ! Loop over all PL atoms
          do i = atomst(m), atomst(m+1)-1

             ! Loop over all contacts
             do j1 = 1, ncont

                iatc1 = transpar%contacts(j1)%idxrange(1)
                iatc2 = transpar%contacts(j1)%idxrange(2)

                i1 = minval(img2CentCell(iNeigh(1:nNeigh(i),i)), &
                    & mask = (img2CentCell(iNeigh(1:nNeigh(i),i)).ge.iatc1 .and. &
                    & img2CentCell(iNeigh(1:nNeigh(i),i)).le.iatc2) )

                if (i1.ge.iatc1 .and. i1.le.iatc2) then
                    minv(m,j1) = j1
                endif

             end do
          end do
       end do


       do j1 = 1, ncont

         if (all(minv(:,j1) == 0)) then
           write(stdOut,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(stdOut,*) 'WARNING: contact',j1,' does no interact with any PL '
           write(stdOut,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           minv(1,j1) = j1
         end if

         if (count(minv(:,j1).eq.j1).gt.1) then
           write(stdOut,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(stdOut,*) 'ERROR: contact',j1,' interacts with more than one PL'
           write(stdOut,*) '       check structure and increase PL size         '
           write(stdOut,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           call error("")
         end if

         do m = 1, transpar%nPLs
           if (minv(m,j1).eq.j1) then
             cblk(j1) = m
           end if
         end do

       end do

      else

         cblk=transpar%cblk

      end if

      write(stdOut,*) ' Structure info:'
      write(stdOut,*) ' Number of PLs:',nbl
      write(stdOut,*) ' PLs coupled to contacts:',cblk(1:ncont)
      write(stdOut,*)

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
  ! Subroutine to check the PL definitions
  !------------------------------------------------------------------------------
  subroutine check_pls(transpar, greendens, natoms, iNeigh, nNeigh, img2CentCell, info)
    Type(TTranspar), intent(in) :: transpar
    Type(TNEGFGreenDensInfo) :: greendens
    integer, intent(in) :: natoms
    Integer, intent(in) :: iNeigh(0:,:)
    Integer, intent(in) :: nNeigh(:)
    Integer, intent(in) :: img2CentCell(:)
    !> error message
    integer, intent(out) :: info

    integer :: nbl, iatm1, iatm2, iats, iate
    integer :: mm, nn, ii, kk
    integer, allocatable :: atomst(:)

    ! The contacts have been already checked
    ! Here checks the PL definition and contact/device interactions

    iatm1 = transpar%idxdevice(1)
    iatm2 = transpar%idxdevice(2)

    nbl = 0

    if (transpar%defined) then
       nbl = transpar%nPLs
    else if (greendens%defined) then
       nbl = greendens%nPLs
    endif

    if (nbl.eq.0) then
      call error('Internal ERROR: nbl = 0 ?!')
    end if

    allocate(atomst(nbl+1))

    if (transpar%defined) then
      atomst(1:nbl) = transpar%PL(1:nbl)
      atomst(nbl+1) = iatm2 + 1
    else if (greendens%defined) then
      atomst(1:nbl) = greendens%PL(1:nbl)
      atomst(nbl+1) = natoms + 1
    endif

    info = 0
    do mm = 1, nbl-1
       do nn = mm+1, nbl
         iats = atomst(nn)
         iate = atomst(nn+1)-1
         do ii = atomst(mm), atomst(mm+1)-1
            kk = maxval( img2CentCell(iNeigh(1:nNeigh(ii),ii)), &
               mask = (img2CentCell(iNeigh(1:nNeigh(ii),ii)).ge.iats .and. &
               img2CentCell(iNeigh(1:nNeigh(ii),ii)).le.iate) )
         end do
         if (nn .gt. mm+1 .and. kk .ge. iats .and. kk .le. iate) then
           write(stdOut,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(stdOut,*) 'WARNING: PL ',mm,' interacts with PL',nn
           write(stdOut,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           info = mm
         end if
       end do
    end do

    deallocate(atomst)

  end subroutine check_pls


  !------------------------------------------------------------------------------
  ! INTERFACE subroutine to call Density matrix computation
  !------------------------------------------------------------------------------
  subroutine negf_density(miter,spin,nkpoint,HH,SS,mu,DensMat,EnMat)

    integer, intent (in) :: miter          ! SCC step (used in SGF)
    integer, intent (in) :: spin           ! spin component (SGF)
    integer, intent (in) :: nkpoint        ! nk point (used in SGF)
    type(z_CSR), pointer, intent(in) :: HH          ! Hamiltonian
    type(z_CSR), pointer, intent(in) :: SS          ! Overlap
    real(dp), intent(in) :: mu(:)
    type(z_CSR), pointer, intent(in), optional :: DensMat   ! Density matrix (See NOTE)
    type(z_CSR), pointer, intent(in), optional :: EnMat     ! Energy weighted DM (See NOTE)

    type(lnParams) :: params
    integer :: nn

    call get_params(negf, params)

    params%iteration = miter
    params%kpoint = nkpoint
    params%spin = spin
    params%DorE='N'
    nn=size(mu,1)
    params%mu(1:nn) = mu(1:nn)

    if(present(DensMat)) then
       params%DorE = 'D'
       call set_params(negf,params)
       call pass_DM(negf,rho=DensMat)
    endif
    if(present(EnMat)) then
       params%DorE = 'E'
       call set_params(negf,params)
       call pass_DM(negf,rhoE=EnMat)
    endif
    if (present(DensMat).and.present(EnMat)) then
       params%DorE  = 'B'
       call set_params(negf,params)
       stop 'UNSUPPORTED CASE in negf_density'
    endif

    if (params%DorE.eq.'N') then
      return
    end if

    call pass_HS(negf,HH,SS)

    call compute_density_dft(negf)

    call destroy_matrices(negf)

  end subroutine negf_density

  !------------------------------------------------------------------------------
  ! INTERFACE subroutine to call ldos computation
  !------------------------------------------------------------------------------
  subroutine negf_ldos(HH,SS,spin,kpoint,wght,ledos)
    type(z_CSR), pointer, intent(in) :: HH, SS
    integer, intent(in) :: spin      ! spin index
    integer, intent(in) :: kpoint        ! kp index
    real(dp), intent(in) :: wght      ! kp weight
    real(dp), dimension(:,:), pointer :: ledos
    type(lnParams) :: params

    call get_params(negf, params)

    params%spin = spin
    params%kpoint = kpoint
    params%wght = wght

    call pass_HS(negf,HH,SS)

    call compute_ldos(negf)

    call destroy_matrices(negf)

    call associate_ldos(pNegf, ledos)

  end subroutine negf_ldos


  !> DEBUG routine to dump H and S as a file in Matlab format
  subroutine negf_dumpHS(HH,SS)

    !> hamiltonian in CSR format
    type(z_CSR), intent(in) :: HH

    !> Overlap in CSR format
    type(z_CSR), intent(in) :: SS

    integer :: fdUnit

    write(stdOut, *) 'Dumping H and S in files...'

    open(newunit=fdUnit, file='HH.dat')
    write(fdUnit, *) '% Size =',HH%nrow, HH%ncol
    write(fdUnit, *) '% Nonzeros =',HH%nnz
    write(fdUnit, *) '% '
    write(fdUnit, *) 'zzz = ['
    call printcsr(fdUnit, HH)
    write(fdUnit, *) ']'
    close(fdUnit)

    open(newunit=fdUnit, file='SS.dat')
    write(fdUnit, *) '% Size =',SS%nrow, SS%ncol
    write(fdUnit, *) '% Nonzeros =',SS%nnz
    write(fdUnit, *) '% '
    write(fdUnit, *) 'zzz = ['
    call printcsr(fdUnit, SS)
    write(fdUnit, *) ']'
    close(fdUnit)

  end subroutine negf_dumpHS

  !------------------------------------------------------------------------------
  ! ALEX: DAR routines to set H and S have been moved here
  !------------------------------------------------------------------------------
  subroutine prepare_HS(H_dev,S_dev,HH,SS)
    real(dp), dimension(:,:) :: H_dev, S_dev
    type(z_CSR), intent(inout) :: HH, SS

    if (negf%tOrthonormal) then
      write(stdOut, "(' Lowdin orthogonalization for the whole system ')")
      call Orthogonalization(H_dev, S_dev)
    end if

    if (negf%tOrthonormalDevice) then
      write(stdOut, "(' Lowdin orthogonalization for device-only')")
      call Orthogonalization_dev(H_dev, S_dev)
    end if

    if (negf%tOrthonormal.or.negf%tOrthonormalDevice) then
      call MakeHHSS(H_dev,S_dev,HH,SS)
    end if

    !if(negf%tManyBody) call MakeHS_dev

  end subroutine
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  ! INTERFACE subroutine to call current computation
  !------------------------------------------------------------------------------
  subroutine negf_current(HH,SS,spin,kpoint,wght,tunn,curr,ledos,currents)

    type(z_CSR), pointer, intent(in) :: HH, SS
    integer, intent(in) :: spin      ! spin index
    integer, intent(in) :: kpoint        ! kp index
    real(dp), intent(in) :: wght      ! kp weight
    real(dp), dimension(:,:), pointer :: tunn
    real(dp), dimension(:,:), pointer :: curr
    real(dp), dimension(:,:), pointer :: ledos
    real(dp), dimension(:), pointer :: currents

    type(lnParams) :: params

    call get_params(negf, params)

    params%spin = spin
    params%kpoint = kpoint
    params%wght = wght

    call set_params(negf, params)

    call pass_HS(negf,HH,SS)

    call compute_current(negf)

    ! Associate internal negf arrays to local pointers
    call associate_ldos(pNegf, ledos)
    call associate_transmission(pNegf, tunn)
    call associate_current(pNegf, curr)

    call associate_lead_currents(pNegf, currents)
    if (.not.associated(currents)) then
      call error('Internal error: currVec not associated')
    end if

 end subroutine negf_current



  !> Calculates density matrix with Green's functions
 subroutine calcdensity_green(iSCCIter, mpicomm, groupKS, ham, over, iNeighbor, nNeighbor,&
     & iAtomStart, iPair, img2CentCell, iCellVec, cellVec, orb, kPoints, kWeights, mu, rho, Eband,&
     & Ef, E0, TS)

    integer, intent(in) :: iSCCIter
    type(mpifx_comm), intent(in) :: mpicomm
    integer, intent(in) :: groupKS(:,:)
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: kPoints(:,:), kWeights(:)
    real(dp), intent(in) :: mu(:,:)
    real(dp), intent(out) :: rho(:,:)
    real(dp), intent(out) :: Eband(:), Ef(:), E0(:), TS(:)

    integer :: nSpin, nKS, iK, iS, iKS
    type(z_CSR), target :: csrDens
    type(z_CSR), pointer :: pCsrDens

    pCsrDens => csrDens

    call negf_mpi_init(mpicomm)

    ! We need this now for different fermi levels in colinear spin
    ! Note: the spin polirized does not work with
    ! built-int potentials (the unpolarized does) in the poisson
    nKS = size(groupKS, dim=2)
    nSpin = size(ham, dim=2)
    rho = 0.0_dp

    write(stdOut, *)
    write(stdOut, '(80("="))')
    write(stdOut, *) '                         COMPUTING DENSITY MATRIX      '
    write(stdOut, '(80("="))')


    do iKS = 1, nKS
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS)

      write(stdOut,*) 'k-point',iK,'Spin',iS

      call foldToCSR(csrHam, ham(:,iS), kPoints(:,iK), iAtomStart, iPair, iNeighbor, nNeighbor,&
          & img2CentCell, iCellVec, cellVec, orb)
      call foldToCSR(csrOver, over, kPoints(:,ik), iAtomStart, iPair, iNeighbor, nNeighbor,&
          & img2CentCell, iCellVec, cellVec, orb)

      call negf_density(iSCCIter, iS, iKS, pCsrHam, pCsrOver, mu(:,iS), DensMat=pCsrDens)

      ! NOTE:
      ! unfold adds up to rho the csrDens(k) contribution
      !
      call unfoldFromCSR(rho(:,iS), csrDens, kPoints(:,iK), kWeights(iK), iAtomStart, iPair,&
          & iNeighbor, nNeighbor, img2CentCell, iCellVec, cellVec, orb)

      call destruct(csrDens)

      ! Set some fake energies:
      Eband(iS) = 0.0_dp
      Ef(iS) = mu(1,iS)
      TS(iS) = 0.0_dp
      E0(iS) = 0.0_dp

    end do

    do iS = 1, nSpin
      ! In place all-reduce of the density matrix
      call mpifx_allreduceip(mpicomm, rho(:,iS), MPI_SUM)
    end do

    write(stdOut,'(80("="))')
    write(stdOut,*)

  end subroutine calcdensity_green

  !> Calculates E-density matrix with Green's functions
  subroutine calcEdensity_green(iSCCIter, mpicomm, groupKS, ham, over, iNeighbor, nNeighbor,&
      & iAtomStart, iPair, img2CentCell, iCellVec, cellVec, orb, kPoints, kWeights, mu, rhoE)

    integer, intent(in) :: iSCCIter
    type(mpifx_comm), intent(in) :: mpicomm
    integer, intent(in) :: groupKS(:,:)
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb   !Needs only orb%nOrbAtom, orb%mOrb
    real(dp), intent(in) :: kPoints(:,:), kWeights(:)
    real(dp), intent(in) :: mu(:,:)
    real(dp), intent(out) :: rhoE(:)

    integer :: nSpin, nKS, iK, iS, iKS
    type(z_CSR), target :: csrEDens
    type(z_CSR), pointer :: pCsrEDens

    pCsrEDens => csrEDens

    call negf_mpi_init(mpicomm)
    ! We need this now for different fermi levels in colinear spin
    ! Note: the spin polirized does not work with
    ! built-int potentials (the unpolarized does) in the poisson
    ! I do not set the fermi because it seems that in libnegf it is
    ! not really needed

    nKS = size(groupKS, dim=2)
    nSpin = size(ham, dim=2)
    rhoE = 0.0_dp

    write(stdOut, *)
    write(stdOut, '(80("="))')
    write(stdOut, *) '                     COMPUTING E-WEIGHTED DENSITY MATRIX '
    write(stdOut, '(80("="))')

    do iKS = 1, nKS
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS)

      write(stdOut,*) 'k-point',iK,'Spin',iS

      call foldToCSR(csrHam, ham(:,iS), kPoints(:,iK), iAtomStart, iPair, iNeighbor, nNeighbor,&
          & img2CentCell, iCellVec, cellVec, orb)
      call foldToCSR(csrOver, over, kPoints(:,ik), iAtomStart, iPair, iNeighbor, nNeighbor,&
          & img2CentCell, iCellVec, cellVec, orb)

      call negf_density(iSCCIter, iS, iKS, pCsrHam, pCsrOver, mu(:,iS), EnMat=pCsrEDens)

      ! NOTE:
      ! unfold adds up to rhoEPrim the csrEDens(k) contribution
      !
      call unfoldFromCSR(rhoE, csrEDens, kPoints(:,iK), kWeights(iK), iAtomStart, iPair, iNeighbor,&
          & nNeighbor, img2CentCell, iCellVec, cellVec, orb)

      call destruct(csrEDens)

    end do

    ! In place all-reduce of the density matrix
    call mpifx_allreduceip(mpicomm, rhoE, MPI_SUM)

    write(stdOut,'(80("="))')
    write(stdOut,*)

  end subroutine calcEdensity_green

  !> Calculate the partial density of states
  subroutine calcPDOS_green(mpicomm, groupKS, ham, over, iNeighbor, nNeighbor, iAtomStart, iPair,&
      & img2CentCell, iCellVec, cellVec, orb, kPoints, kWeights, ldosTot, writeLDOS)
    integer, intent(in) :: groupKS(:,:)
    type(mpifx_comm), intent(in) :: mpicomm
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: kPoints(:,:), kWeights(:)
    real(dp), allocatable, intent(inout) :: ldosTot(:,:)
    logical, intent(in) :: writeLDOS

    real(dp), pointer    :: ldosMat(:,:)=>null()
    real(dp), allocatable :: ldosSKRes(:,:,:)
    integer :: iKS, iK, iS, nKS, nKPoint, nSpin, ii, jj, err
    type(lnParams) :: params
    integer :: fdUnit

    call negf_mpi_init(mpicomm)

    call get_params(negf, params)

    nKS = size(groupKS, dim=2)
    nKPoint = size(kPoints, dim=2)
    nSpin = size(ham, dim=2)

    write(stdOut, *)
    write(stdOut, '(80("="))')
    write(stdOut, *) '                        COMPUTING  LOCAL  DOS          '
    write(stdOut, '(80("="))')
    write(stdOut, *)

    do iKS = 1, nKS
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS)

      write(stdOut,*) 'k-point',iK,'Spin',iS

      call foldToCSR(csrHam, ham(:,iS), kPoints(:,iK), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)

      call foldToCSR(csrOver, over, kPoints(:,ik), iAtomStart, &
          &iPair, iNeighbor, nNeighbor, img2CentCell, &
          &iCellVec, cellVec, orb)

      call negf_ldos(pCsrHam, PCsrOver, iS, iK, kWeights(iK), ldosMat)

      call add_partial_results(mpicomm, ldosMat, ldosTot, ldosSKRes, iKS, nKS)

    end do

    if (allocated(ldosTot)) then

      ! Write Total localDOS on a separate file (optional)
      if (tIoProc .and. writeLDOS) then
        open(newunit=fdUnit, file='localDOS.dat')
        do ii=1,size(ldosTot,1)
          write(fdUnit,*) (params%Emin+(ii-1)*params%Estep) * Hartree__eV, &
              (ldosTot(ii,jj), jj=1,size(ldosTot,2))
        enddo
        close(fdUnit)
      endif

    else

      allocate(ldosTot(0,0))

    endif


  end subroutine calcPDOS_green



  !> Calculate the current
  subroutine calc_current(mpicomm, groupKS, ham, over, iNeighbor, nNeighbor, iAtomStart, iPair,&
      & img2CentCell, iCellVec, cellVec, orb, kPoints, kWeights, tunnMat, currMat, ldosMat,&
      & currLead, writeTunn, writeLDOS, mu)

    integer, intent(in) :: groupKS(:,:)
    type(mpifx_comm), intent(in) :: mpicomm
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)
    integer, intent(in) :: img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: kPoints(:,:), kWeights(:)
    real(dp), allocatable, intent(inout) :: tunnMat(:,:)
    real(dp), allocatable, intent(inout) :: currMat(:,:)
    real(dp), allocatable, intent(inout) :: ldosMat(:,:)
    real(dp), allocatable, intent(inout) :: currLead(:)
    logical, intent(in) :: writeLDOS
    logical, intent(in) :: writeTunn
    ! We need this now for different fermi levels in colinear spin
    ! Note: the spin polirized does not work with
    ! built-int potentials (the unpolarized does) in the poisson
    ! I do not set the fermi because it seems that in libnegf it is
    ! not really needed
    real(dp), intent(in) :: mu(:,:)

    real(dp), allocatable :: tunnSKRes(:,:,:), currSKRes(:,:,:), ldosSKRes(:,:,:)
    real(dp), pointer    :: tunnPMat(:,:)=>null()
    real(dp), pointer    :: currPMat(:,:)=>null()
    real(dp), pointer    :: ldosPMat(:,:)=>null()
    real(dp), pointer    :: currPVec(:)=>null()
    integer :: iKS, iK, iS, nKS, ii, err, ncont
    type(unit) :: unitOfEnergy        ! Set the units of H
    type(unit) :: unitOfCurrent       ! Set desired units for Jel
    type(lnParams) :: params

    integer :: i, j, k, NumStates, icont
    real(dp), dimension(:,:), allocatable :: H_all, S_all
    character(:), allocatable :: filename

    call negf_mpi_init(mpicomm)

    call get_params(negf, params)

    unitOfEnergy%name = "H"
    unitOfCurrent%name = "A"

    nKS = size(groupKS, dim=2)
    ncont = size(mu,1)

    if (params%verbose.gt.30) then
      write(stdOut, *)
      write(stdOut, '(80("="))')
      write(stdOut, *) '                            COMPUTATION OF CURRENT         '
      write(stdOut, '(80("="))')
      write(stdOut, *)
    end if

    do iKS = 1, nKS
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS)

      write(stdOut,*) 'Spin',iS,'k-point',iK,'k-weight',kWeights(iK)

      params%mu(1:ncont) = mu(1:ncont,iS)

      call set_params(negf, params)

      if (negf%NumStates.eq.0) then
        negf%NumStates=csrHam%ncol
      end if

      !*** ORTHOGONALIZATIONS ***
      ! THIS MAKES SENSE ONLY FOR A REAL MATRICES, i.e. k==0 && collinear spin
      if (all(kPoints(:,iK) .eq. 0.0_dp) .and. &
         (negf%tOrthonormal .or. negf%tOrthonormalDevice)) then

        NumStates = negf%NumStates

        if (.not.allocated(H_all)) then
          allocate(H_all(NumStates,NumStates))
        end if
        if (.not.allocated(S_all)) then
          allocate(S_all(NumStates,NumStates))
        end if

        call unpackHS(H_all, ham(:,iS), iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell)
        call blockSymmetrizeHS(H_all, iAtomStart)

        call unpackHS(S_all, over, iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell)
        call blockSymmetrizeHS(S_all, iAtomStart)

        call prepare_HS(H_all,S_all,csrHam,csrOver)

      else

        call foldToCSR(csrHam, ham(:,iS), kPoints(:,iK), iAtomStart, iPair, iNeighbor, nNeighbor,&
            & img2CentCell, iCellVec, cellVec, orb)

        call foldToCSR(csrOver, over, kPoints(:,ik), iAtomStart, iPair, iNeighbor, nNeighbor,&
            & img2CentCell, iCellVec, cellVec, orb)

      end if

      call negf_current(pCsrHam, pCsrOver, iS, iK, kWeights(iK), &
                       & tunnPMat, currPMat, ldosPMat, currPVec)

      if(.not.allocated(currLead)) then
        allocate(currLead(size(currPVec)), stat=err)
        if (err /= 0) then
          call error('Allocation error (currTot)')
        end if
        currLead = 0.0_dp
      endif
      currLead = currLead + currPVec

      !GUIDE: tunnPMat libNEGF output stores Transmission, T(iE, i->j)
      !       tunnPMat is MPI distributed on energy points (0.0 on other nodes)
      !       tunnMat MPI gather partial results and accumulate k-summation
      !       currPMat stores contact current I_i(iE)
      !       tunnSKRes stores tunneling for all k-points and spin: T(iE, i->j, iSK)
      call add_partial_results(mpicomm, tunnPMat, tunnMat, tunnSKRes, iKS, nKS)

      call add_partial_results(mpicomm, currPMat, currMat, currSKRes, iKS, nKS)

      call add_partial_results(mpicomm, ldosPMat, ldosMat, ldosSKRes, iKS, nKS)

    end do

    call mpifx_allreduceip(mpicomm, currLead, MPI_SUM)

    call mpifx_barrier(mpicomm)

    ! converts from internal atomic units into A
    currLead = currLead * convertCurrent(unitOfEnergy, unitOfCurrent)

    do ii = 1, size(currLead)
      write(stdOut, *)
      write(stdOut, '(1x,a,i3,i3,a,ES14.5,a,a)') ' contacts: ',params%ni(ii),params%nf(ii),&
          & ' current: ', currLead(ii),' ',unitOfCurrent%name
    enddo

    ! Write Total transmission, T(E), on a separate file (optional)
    if (allocated(tunnMat)) then
      filename = 'transmission'
      if (tIOProc .and. writeTunn) then
        call write_file(negf, tunnMat, tunnSKRes, filename, groupKS, kpoints, kWeights)
      end if
      if (allocated(tunnSKRes)) then
        deallocate(tunnSKRes)
      end if
    else
      ! needed to avoid some segfault
      allocate(tunnMat(0,0))
    end if

    ! Write Total lead current, I_i(E), on a separate file (optional)
    if (allocated(currMat)) then
      filename = 'current'
      if (tIOProc .and. writeTunn) then
        call write_file(negf, currMat, currSKRes, filename, groupKS, kpoints, kWeights)
      end if
      if (allocated(currSKRes)) then
        deallocate(currSKRes)
      end if
    else
      ! needed to avoid some segfault
      allocate(currMat(0,0))
    endif

    if (allocated(ldosMat)) then
      ! Write Total localDOS on a separate file (optional)
      if (tIoProc .and. writeLDOS) then
        call write_file(negf, ldosMat, ldosSKRes, 'localDOS', groupKS, kpoints, kWeights)
      end if
      if (allocated(ldosSKRes)) then
        deallocate(ldosSKRes)
      end if
    else
      ! needed to avoid some segfault
      allocate(ldosMat(0,0))
    end if

  end subroutine calc_current

  !----------------------------------------------------------------------------
  !   utility to allocate and sum partial results
  !----------------------------------------------------------------------------
  subroutine add_partial_results(mpicomm, pMat, matTot, matSKRes, iK, nK)
    type(mpifx_comm), intent(in) :: mpicomm
    real(dp), intent(in), pointer :: pMat(:,:)
    real(dp), allocatable :: matTot(:,:)
    real(dp), allocatable :: matSKRes(:,:,:)
    integer, intent(in) :: iK, nK

    real(dp), allocatable :: tmpMat(:,:)
    integer :: err

    if (associated(pMat)) then
      allocate(tmpMat(size(pMat,1), size(pMat,2)), stat=err)

      if (err /= 0) then
        call error('Allocation error (tmpMat)')
      end if

      tmpMat = pMat
      call mpifx_allreduceip(mpicomm, tmpMat, MPI_SUM)

      if(.not.allocated(matTot)) then
        allocate(matTot(size(pMat,1), size(pMat,2)), stat=err)

        if (err /= 0) then
          call error('Allocation error (tunnTot)')
        end if

        matTot = 0.0_dp
      end if
      matTot = matTot + tmpMat

      if (nK > 1) then
        if (.not.allocated(matSKRes)) then
          allocate(matSKRes(size(pMat,1), size(pMat,2), nK), stat=err)

          if (err/=0) then
            call error('Allocation error (tunnSKRes)')
          end if

          matSKRes = 0.0_dp
        endif
        matSKRes(:,:,iK) = tmpMat(:,:)
      end if

      deallocate(tmpMat)

    end if

  end subroutine add_partial_results


  !> utility to write tunneling or ldos on files
  subroutine write_file(negf, matTot, matSKRes, filename, groupKS, kpoints, kWeights)

    !> Contains input data, runtime quantities and output data
    type(TNegf) :: negf

    !> results to print if allocated
    real(dp), intent(in), allocatable :: matTot(:,:)

    !> k- and spin-resolved quantities, if allocated
    real(dp), intent(in), allocatable :: matSKRes(:,:,:)

    !> file to print out to
    character(*), intent(in) :: filename

    !> local k-points and spins on this processor
    integer, intent(in) :: groupKS(:,:)

    !> k-points
    real(dp), intent(in) :: kPoints(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeights(:)

    integer :: ii, jj, nKS, iKS, nK, nS, iK, iS, fdUnit
    type(lnParams) :: params

    call get_params(negf, params)

    nKS = size(groupKS, dim=2)
    nK = size(kpoints, dim=2)
    nS = nKS/nK

    open(newunit=fdUnit, file=trim(filename)//'.dat')
    do ii=1,size(matTot, dim=1)
      write(fdUnit,'(F20.6)',ADVANCE='NO') (params%Emin+(ii-1)*params%Estep) * Hartree__eV
      do jj=1,size(matTot, dim=2)
        write(fdUnit,'(ES20.8)',ADVANCE='NO') matTot(ii,jj)
      enddo
      write(fdUnit,*)
    enddo
    close(fdUnit)

    if (nKS > 1) then

      open(newunit=fdUnit, file=trim(filename)//'_kpoints.dat')
      write(fdUnit,*)'# NKpoints = ', nK
      write(fdUnit,*)'# NSpin = ', nS
      write(fdUnit,*)'# Energy [eV], <spin k1 k2 k3 weight> '
      write(fdUnit,'(A1)', ADVANCE='NO') '# '

      do iKS = 1, nKS
        iK = groupKS(1,iKS)
        iS = groupKS(2,iKS)
        write(fdUnit,'(i5.2)', ADVANCE='NO') iS
        write(fdUnit,'(es15.5, es15.5, es15.5, es15.5)', ADVANCE='NO') kpoints(:,iK), kWeights(iK)
      end do
      write(fdUnit,*)

      if (allocated(matSKRes)) then
        do ii = 1, size(matSKRes(:,:,:), dim=1)
          write(fdUnit,'(f20.6)',ADVANCE='NO') (params%Emin+(ii-1)*params%Estep) * Hartree__eV
          do jj = 1, size(matSKRes(:,:,:), dim=2)
            do iKS = 1,nKS
              write(fdUnit,'(es20.8)',ADVANCE='NO') matSKRes(ii,jj, iKS)
            enddo
            write(fdUnit,*)
          enddo
        enddo
      end if
      close(fdUnit)

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
    ! We need this now for different fermi levels in colinear spin
    ! Note: spin polarized does not work with
    ! built-int potentials (the unpolarized does) in the poisson
    ! I do not set the fermi because it seems that in libnegf it is
    ! not really needed
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
    type(z_CSR), target :: csrDens, csrEDens
    type(z_CSR), pointer :: pCsrDens, pCsrEDens
    type(lnParams) :: params
    integer :: fdUnit

    pCsrDens => csrDens
    pCsrEDens => csrEDens

    call negf_mpi_init(mpicomm)

    call get_params(negf, params)

    write(stdOut, *)
    write(stdOut, '(80("="))')
    write(stdOut, *) '                        COMPUTING LOCAL CURRENTS          '
    write(stdOut, '(80("="))')
    write(stdOut, *)

    nKS = size(groupKS, dim=2)
    nKPoint = size(kPoints, dim=2)
    nSpin = size(ham, dim=2)
    nAtom = size(orb%nOrbAtom)

    do iKS = 1, nKS
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS)

      write(stdOut,*) 'k-point',iK,'Spin',iS

      ! We need to recompute Rho and RhoE .....
      call foldToCSR(csrHam, ham(:,iS), kPoints(:,1), iAtomStart, iPair, iNeighbor, nNeighbor,&
          & img2CentCell, iCellVec, cellVec, orb)
      call foldToCSR(csrOver, over, kPoints(:,1), iAtomStart, iPair, iNeighbor, nNeighbor,&
          & img2CentCell, iCellVec, cellVec, orb)

      call negf_density(iSCCIter, iS, iKS, pCsrHam, pCsrOver, chempot(:,iS), DensMat=pCsrDens)

      call negf_density(iSCCIter, iS, iKS, pCsrHam, pCsrOver, chempot(:,iS), EnMat=pCsrEDens)

      call mpifx_allreduceip(mpicomm, csrDens%nzval, MPI_SUM)

      call mpifx_allreduceip(mpicomm, csrEDens%nzval, MPI_SUM)

      allocate(nneig(size(iNeighbor,2), NMAX))
      allocate(nn(size(iNeighbor,2)))
      call symmetrize_neiglist(nAtom, iNeighbor, nNeighbor, img2CentCell, coord, nneig, nn)


      if (iS .eq. 1) then
        sp = 'u'
      end if
      if (iS .eq. 2) then
        sp = 'd'
      end if

      open(newUnit = fdUnit, file = 'lcurrents_'//sp//'.dat')

      do m = 1, nAtom
        allocate(Im(nn(m)))
        Im(:) = 0.0_dp

        mOrb = orb%nOrbAtom(m)
        iRow = iAtomStart(m)

        write(fdUnit,'(I5,3(F12.6),I3)',advance='NO') m, coord(:,m), nn(m)

        do inn = 1, nn(m)
          n = nneig(m,inn)
          startn = iAtomStart(img2CentCell(n))
          endn = startn + orb%nOrbAtom(n) - 1
          ! tracing orbitals of atoms  n  m
          ! More efficient without getel ?
          do mu = iRow, iRow+mOrb-1
            do nu = startn, endn
              c1 = conjg(getel(csrDens,mu,nu))
              c2 = conjg(getel(csrEDens,mu,nu))
              Im(inn) = Im(inn) + aimag(getel(csrHam,mu,nu)*c1 - getel(csrOver,mu,nu)*c2)
            enddo
          enddo
          ! pi-factor  comes from  Gn = rho * pi
          write(fdUnit,'(I5,ES17.8)',advance='NO') n, 2.0_dp*params%g_spin*pi*eovh*Im(inn)

        enddo

        write(fdUnit,*)

        deallocate(Im)

      enddo
    enddo

    close(fdUnit)
    call destruct(csrDens)
    call destruct(csrEDens)
    deallocate(nneig)
    deallocate(nn)

  end subroutine local_currents

  !--------------------------------------------------------------
  ! Neighbor is non-symmetric: i->j  j>i
  ! The following routine symmetrizes the neighbor list
  !
  subroutine symmetrize_neiglist(nAtom,iNeighbor,nNeighbor,img2CentCell,coord,nneig,nn)
    integer, intent(in) :: nAtom
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: coord(:,:)
    integer, dimension(:,:), intent(out) :: nneig
    integer, dimension(:), intent(out) :: nn

    real(dp) :: dist, dr(3)
    integer :: m, n, inn, ii, jj, morb
    integer, parameter :: NMAX=40

    nn=0
    nneig=0

    do m = 1, nAtom
      do inn = 1, nNeighbor(m)
        if (inn.gt.NMAX) then
          exit
        end if
        n = img2CentCell(iNeighbor(inn,m))
        !n = iNeighbor(inn,m)
        if(nn(m).lt.NMAX) then
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
        if(nn(n).lt.NMAX) then
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
    integer :: i,j,k,l,m,n,NumStates, fdUnit

    NumStates = negf%NumStates

    write(stdOut,"(' Hamiltonian is read from the file ',A)")trim('H_dftb.mtr')


    open(newunit=fdUnit, file='H_dftb.mtr', action="read")
    read(fdUnit,*);read(fdUnit,*);read(fdUnit,*);read(fdUnit,*);read(fdUnit,*)
    if(.not.allocated(H)) allocate(H(NumStates,NumStates))
    H = 0.0_dp
    do i=1,NumStates
       read(fdUnit,*) H(i,1:NumStates)
    end do
    close(fdUnit)

    H = H * eV__Hartree
    write(stdOut,"(' Overlap is read from the file ',A)")trim('S_dftb.mtr')
    open(newunit=fdUnit, file='S_dftb.mtr', action="read")
    read(fdUnit,*);read(fdUnit,*);read(fdUnit,*);read(fdUnit,*);read(fdUnit,*)
    if(.not.allocated(S)) allocate(S(NumStates,NumStates))
    S=0.0_dp
    do i=1,NumStates
       read(fdUnit,*) S(i,1:NumStates)
    end do
    close(fdUnit)

  end subroutine ReadDFTB

  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------

  subroutine ReadModel(H,S)
    real(dp), dimension(:,:), allocatable :: H, S
    integer :: i,j,k,l,m,n,NumStates, fdUnit

    NumStates = negf%NumStates

    write(stdOut,"(' Hamiltonian is read from the file ',A)")trim('H.mtr')
    open(newunit = fdUnit,file='H.mtr',action="read")
    if(.not.allocated(H)) allocate(H(NumStates,NumStates))
    H=0.0_dp
    do i=1,NumStates
      read(fdUnit,*) H(i,1:NumStates)
    end do
    close(fdUnit)
    H=H * eV__Hartree
    if(.not.allocated(S)) allocate(S(NumStates,NumStates))
    S=0.0_dp
    do i=1,NumStates
      S(i,i)=1._dp
    end do

  end subroutine ReadModel

  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------

  subroutine ReadLibNEGF

    write(stdOut, *) 'Import from the negf.in file'
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
        if ((i.eq.j).or.(abs(H_all(i,j)).gt.0.00001_dp)) then
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
          if((i.eq.j).or.(abs(H_all(i,j)).gt.0.00001_dp)) then
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

    if (negf%verbose.gt.30) then
      write(stdOut, *)
      write(stdOut, '(80("="))')
      write(stdOut, *) '                   LibNEGF: equilibrium (!) LDOS calculation'
      write(stdOut, '(80("="))')
      write(stdOut, *)
    endif

    call compute_ldos(negf)

    if (negf%verbose.gt.30) then
      write(stdOut, *)
      write(stdOut, '(80("="))')
      write(stdOut, *) '                     LibNEGF: equilibrium (!) LDOS finished'
      write(stdOut, '(80("="))')
    endif

  end subroutine WriteLDOS

  !-----------------------------------------------------------------------------

  subroutine check_negf_params()

    use globals, only : LST

    Integer :: ncont, nbl, ii, jj, ist, iend, fdUnit
    Integer, dimension(:), allocatable :: PL_end, cont_end, surf_end
    character(32) :: tmp
    character(LST) :: file_re_H, file_im_H, file_re_S, file_im_S
    integer, dimension(:), allocatable :: cblk

    open(newunit=fdUnit, file='log', action="write")

    write(fdUnit,"('negf%H%nnz               = ',I6)") negf%H%nnz
    write(fdUnit,"('negf%H%nrow              = ',I6)") negf%H%nrow
    write(fdUnit,"('negf%H%ncol              = ',I6)") negf%H%ncol
    write(fdUnit,"('negf%H%nzval             = ',10000000F8.4)") negf%H%nzval
    write(fdUnit,"('negf%H%colind            = ',10000000I6)") negf%H%colind
    write(fdUnit,"('negf%H%rowpnt            = ',10000000I6)") negf%H%rowpnt
    write(fdUnit,"('negf%isSid               = ',L6)") negf%isSid
    write(fdUnit,"('negf%S%nnz               = ',I6)") negf%S%nnz
    write(fdUnit,"('negf%S%nrow              = ',I6)") negf%S%nrow
    write(fdUnit,"('negf%S%ncol              = ',I6)") negf%S%ncol
    write(fdUnit,"('negf%S%nzval             = ',10000000F8.4)") negf%S%nzval
    write(fdUnit,"('negf%S%colind            = ',10000000I6)") negf%S%colind
    write(fdUnit,"('negf%S%rowpnt            = ',10000000I6)") negf%S%rowpnt
    write(fdUnit,"('negf%str%num_conts       = ',I6)") negf%str%num_conts
    write(fdUnit,"('negf%str%num_PLs         = ',I6)") negf%str%num_PLs
    write(fdUnit,"('negf%str%active_cont     = ',I6)") negf%str%active_cont
    write(fdUnit,"('negf%str%mat_B_start     = ',I6,I6)") negf%str%mat_B_start
    write(fdUnit,"('negf%str%mat_C_start     = ',I6,I6)") negf%str%mat_C_start
    write(fdUnit,"('negf%str%mat_C_end       = ',I6,I6)") negf%str%mat_C_end
    write(fdUnit,"('negf%str%cblk            = ',I6,I6)") negf%str%cblk
    write(fdUnit,"('negf%str%cont_dim        = ',I6,I6)") negf%str%cont_dim
    write(fdUnit,"('negf%str%mat_PL_start    = ',10000000I6)") negf%str%mat_PL_start
    write(fdUnit,"('negf%str%mat_PL_end      = ',10000000I6)") negf%str%mat_PL_end
    write(fdUnit,"('negf%str%central_dim     = ',I6)") negf%str%central_dim
    write(fdUnit,"('negf%str%total_dim       = ',I6)") negf%str%total_dim
    write(fdUnit,"('negf%Ec, negf%Ev         = ',3F8.4)") negf%Ec, negf%Ev
    write(fdUnit,"('negf%DeltaEc, %DeltaEv   = ',3F8.4)") negf%DeltaEc, negf%DeltaEv
    write(fdUnit,"('negf%Emin, %Emax, %Estep = ',3F8.4)") negf%Emin, negf%Emax, negf%Estep
    write(fdUnit,"('negf%kbT_dm              = ',10000000E12.4)") negf%cont(:)%kbT_dm
    write(fdUnit,"('negf%kbT_t               = ',10000000E12.4)") negf%cont(:)%kbT_t
    write(fdUnit,"('negf%mu_n                = ',10000000F8.4)") negf%cont(:)%mu_n
    write(fdUnit,"('negf%mu_p                = ',10000000F8.4)") negf%cont(:)%mu_p
    write(fdUnit,"('negf%mu                  = ',10000000F8.4)") negf%cont(:)%mu
    write(fdUnit,"('negf%delta               = ',10000000E12.4)") negf%delta
    write(fdUnit,"('negf%Np_real             = ',10000000I6)") negf%Np_real
    write(fdUnit,"('negf%n_kt                = ',I6)") negf%n_kt
    write(fdUnit,"('negf%n_poles             = ',I6)") negf%n_poles
    write(fdUnit,"('negf%wght                = ',E12.4)") negf%wght
    write(fdUnit,"('negf%g_spin              = ',E12.4)") negf%g_spin
    write(fdUnit,"('negf%nLDOS               = ',I6)") negf%nLDOS
    do ii = 1, negf%nLDOS
      write(fdUnit,"('negf%LDOS(',I4,')%indexes   = ',10000000I6)") ii, negf%LDOS(ii)%indexes
    end do
    write(fdUnit,"('negf%Np_n                = ',10000000I6)") negf%Np_n
    write(fdUnit,"('negf%Np_p                = ',10000000I6)") negf%Np_p

    close (fdUnit)

    write(stdOut, *)
    write(stdOut, "(' negf parameters are written to the log file')")

  end subroutine check_negf_params


  !-----------------------------------------------------------------------------
  subroutine orthogonalization(H,S)

    real(dp), dimension(:,:) :: H, S
    integer :: i,m,n1_first,n1_last,n2_first,n2_last
    integer :: INFO, N
    real(dp), allocatable :: A(:,:), WORK(:), W(:)
    real(dp), allocatable :: B(:,:),C(:,:)

    N=size(H,1)

    if (N /= negf%NumStates) then
      call error('orthogonalization: negf init NumStates error')
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

    B(:,:) = matmul(transpose(A), matmul(S, A))

    do i=1,N
      B(i,i) = 1.0_dp / sqrt(B(i,i))
    end do

    !C=sqrt(S)
    C(:,:) = matmul(A, matmul(B, transpose(A)))

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

    H(:,:) = matmul(transpose(C), matmul(H, C))

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
    !  write(12,*) H(i,1:N)* Hartree__eV
    !end do
    !close(12)
  end subroutine orthogonalization

  !-----------------------------------------------------------------------------

  subroutine orthogonalization_dev(H, S)

    real(dp), dimension(:,:) :: H, S
    integer :: i,m,n1_first,n1_last,n2_first,n2_last
    integer :: INFO, N, N2
    real(dp), allocatable :: A(:,:), WORK(:), W(:)
    real(dp), allocatable :: B(:,:), U(:,:), C(:,:)


    N=size(H,1)
    if (N /= negf%NumStates) then
      call error('orthogonalization: negf init NumStates error')
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
    !   write(12,*) H(i,1:N)*Hartree__eV
    !end do
    !close(12)

    !Save S_dftb_orth.mtr to file
    !open(12,file='S_dftb_orth.mtr',action="write")
    !do i=1,N
    !   write(12,*) S(i,1:N)
    !end do
    !close(12)

  end subroutine orthogonalization_dev

end module negf_int
