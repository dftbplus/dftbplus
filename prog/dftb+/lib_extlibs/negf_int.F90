!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
#:include "common.fypp"

!> Interface to LIBNEGF for DFTB+
module negf_int
  use libnegf_vars
  use libnegf, only : convertcurrent, eovh, getel, lnParams, pass_DM, Tnegf, units
#:if WITH_MPI
  use libnegf, only : negf_mpi_init
#:endif
  use libnegf, only : z_CSR, READ_SGF, COMP_SGF, COMPSAVE_SGF
  use libnegf, only : associate_lead_currents, associate_ldos, associate_transmission
  use libnegf, only : associate_current, compute_current, compute_density_dft, compute_ldos
  use libnegf, only : create, create_scratch, destroy, set_readoldDMsgf
  use libnegf, only : destroy_matrices, destroy_negf, get_params, init_contacts, init_ldos
  use libnegf, only : init_negf, init_structure, pass_hs, set_bp_dephasing
  use libnegf, only : set_drop, set_elph_block_dephasing, set_elph_dephasing, set_elph_s_dephasing
  use libnegf, only : set_ldos_indexes, set_params, set_scratch, writememinfo, writepeakinfo
  use libnegf, only : printcsr
  use dftbp_accuracy
  use dftbp_environment
  use dftbp_constants
  use dftbp_matconv
  use dftbp_sparse2dense
  use dftbp_densedescr
  use dftbp_commontypes, only : TOrbitals
  use dftbp_formatout
  use dftbp_globalenv, only : stdOut, tIOproc
  use dftbp_message
  use dftbp_elecsolvertypes, only : electronicSolverTypes
  use dftbp_linkedlist
  use dftbp_periodic
  use dftbp_assert
  use dftbp_eigensolver
#:if WITH_MPI
  use dftbp_mpifx
#:endif
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

  !> initialisation of dephasing effects
  public :: negf_init_dephasing

  !> electron-phonon coupling initialisation
  public :: negf_init_elph

  !> clean up CSR matrices and the NEGF structure
  public :: negf_destroy

  !> wrapped functions passing dftb matrices. Needed for parallel
  public :: calcdensity_green

  !> Calculate the energy weighted density from the GF
  public :: calcEdensity_green

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

  !> Format for two values with units
  character(len=*), parameter :: format2U = "(1X,A, ':', T32, F18.10, T51, A, T54, F16.4, T71, A)"

  contains

  !> Init gDFTB environment and variables
  subroutine negf_init(transpar, env, greendens, tundos, tempElec, solver)

    !> Parameters for the transport calculation
    Type(TTranspar), intent(in) :: transpar

    !> Environment settings, suplying the global comm world
    type(TEnvironment), intent(in) :: env

    !> Parameters for the Green's function calculation
    Type(TNEGFGreenDensInfo), intent(in) :: greendens

    !> parameters for tuneling and density of states evaluation
    Type(TNEGFTunDos), intent(in) :: tundos

    !> Electronic temperature
    real(dp), intent(in) :: tempElec

    !> Which solver call is used in the main code
    integer, intent(in) :: solver

    ! local variables
    real(dp), allocatable :: pot(:), eFermi(:)
    integer :: i, l, ncont, nldos
    integer, allocatable :: sizes(:)
    type(lnParams) :: params

    ! Workaround: ifort 16
    ! Pointer must be set within a subroutine. Initialization at declaration fails.
    pNegf => negf
#:if WITH_MPI
    call negf_mpi_init(env%mpi%globalComm)
#:endif

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

    if (tIoProc .and. greendens%saveSGF ) then
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

        write(stdOut,"(1X,A,I0,A)") '(negf_init) CONTACT INFO #', i,&
            & ' "'//trim(transpar%contacts(i)%name)//'"'

        if (params%FictCont(i)) then
          write(stdOut,*) 'FICTITIOUS CONTACT '
          write(stdOut,*) 'DOS: ', params%contact_DOS(i)
        end if
        write(stdOut,*) 'Temperature (DM): ', params%kbT_dm(i)
        write(stdOut,*) 'Temperature (Current): ', params%kbT_t(i)
        if (transpar%contacts(i)%tFermiSet) then
          write(stdOut,format2U)'Potential (with built-in)', pot(i), 'H', Hartree__eV*pot(i), 'eV'
          write(stdOut,format2U)'eFermi', eFermi(i), 'H', Hartree__eV*eFermi(i), 'eV'
        end if
        write(stdOut,*)

        ! Define electrochemical potentials
        params%mu(i) = eFermi(i) - pot(i)

        if (transpar%contacts(i)%tFermiSet) then
          write(stdOut,format2U)'Electro-chemical potentials', params%mu(i), 'H',&
              & Hartree__eV*params%mu(i), 'eV'
          write(stdOut,*)
        end if

      enddo

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
      ! Override to 0 real points if bias is 0.0
      params%Np_real = 0
      if (ncont > 0) then
        if (any(abs(params%mu(2:ncont)-params%mu(1)) > 1.0e-10_dp)) then
           params%Np_real = greendens%nP(3)  ! real axis points
        end if
      end if
     
      ! Setting for Read/Write Surface GFs. 
      ! NOTE: for the moment in tunneling and dos SGF are always   
      ! recomputed because bias may change points and errors are easy

      ! Read G.F. from very first iteration
      if (greendens%readSGF) then
        params%readOldDM_SGFs = READ_SGF
      end if
      ! Compute G.F. at every iteration
      if (.not.greendens%readSGF .and. .not.greendens%saveSGF) then
        params%readOldDM_SGFs = COMP_SGF
      end if
      ! Default Write on first iter
      if (.not.greendens%readSGF .and. greendens%saveSGF) then
        params%readOldDM_SGFs = COMPSAVE_SGF
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
      if (params%readOldDM_SGFs==0) then
        write(stdOut,*) 'Read Existing SGFs: Yes '
      else 
        write(stdOut,*) 'Read Existing SGFs: No, option ', params%readOldDM_SGFs
      end if
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
      
      ! For the moment tunneling and ldos SGFs are always recomputed 
      params%readOldT_SGFs = COMP_SGF

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

    !> Start of orbitals for each atom
    integer, intent(in) :: iAtomStart(:)

    !> neighbours of each atom
    integer, intent(in) :: iNeighbor(0:,:)

    !> number of neighbours for each atom
    integer, intent(in) :: nNeighbor(:)

    !> mapping from image atoms to central cell
    integer, intent(in) :: img2CentCell(:)

    !> atomic orbital information
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


  !> Initialise the structures for the libNEGF library
  subroutine negf_init_str(denseDescr, transpar, greendens, iNeigh, nNeigh, img2CentCell)

    !> Dense matrix information
    Type(TDenseDescr), intent(in) :: denseDescr

    !> transport calculation parameters
    Type(TTranspar), intent(in) :: transpar

    !> Green's function calculational parameters
    Type(TNEGFGreenDensInfo) :: greendens

    !> number of neighbours for each atom
    Integer, intent(in) :: nNeigh(:)

    !> mapping from image atoms to central cell
    Integer, intent(in) :: img2CentCell(:)

    !> neighbours of each atom
    Integer, intent(in) :: iNeigh(0:,:)

    integer, allocatable :: PL_end(:), cont_end(:), surf_start(:), surf_end(:), cblk(:)
    integer, allocatable :: ind(:), atomst(:), plcont(:)
    integer, allocatable :: minv(:,:)
    Integer :: natoms, ncont, nbl, iatc1, iatc2, iatm2
    integer :: i, m, i1, j1, info

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
    allocate(surf_start(ncont))
    allocate(ind(natoms+1))
    allocate(minv(nbl,ncont))

    ind(:) = DenseDescr%iatomstart(:) - 1
    minv = 0
    cblk = 0

    do i = 1, ncont
       cont_end(i) = ind(transpar%contacts(i)%idxrange(2)+1)
       surf_start(i) = ind(transpar%contacts(i)%idxrange(1))+1
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
           write(stdOut,*) 'WARNING: contact',j1,' does not interact with any PL'
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


  !> Subroutine to check the principal layer (PL) definitions
  subroutine check_pls(transPar, greenDens, nAtoms, iNeigh, nNeigh, img2CentCell, info)

    !> transport calculation parameters
    type(TTranspar), intent(in) :: transPar

    !> Green's function calculational parameters
    type(TNEGFGreenDensInfo) :: greenDens

    !> Number of atoms in the system
    integer, intent(in) :: nAtoms

    !> neighbours of each atom
    integer, intent(in) :: iNeigh(0:,:)

    !> number of neighbours for each atom
    integer, intent(in) :: nNeigh(:)

    !> mapping from image atoms to central cell
    integer, intent(in) :: img2CentCell(:)

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


  !> Interface subroutine to call density matrix computation
  subroutine negf_density(miter,spin,nkpoint,HH,SS,mu,DensMat,EnMat)

    !> SCC step (used in SGF)
    integer, intent (in) :: miter

    !> spin component (SGF)
    integer, intent (in) :: spin

    !> nk point (used in SGF)
    integer, intent (in) :: nkpoint

    !> Hamiltonian
    type(z_CSR), pointer, intent(in) :: HH

    !> Overlap
    type(z_CSR), pointer, intent(in) :: SS

    !> chemical potential
    real(dp), intent(in) :: mu(:)

    !> Density matrix
    type(z_CSR), pointer, intent(in), optional :: DensMat

    !> Energy weighted DM
    type(z_CSR), pointer, intent(in), optional :: EnMat

    type(lnParams) :: params
    integer :: nn

    call get_params(negf, params)

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
       call error('UNSUPPORTED CASE in negf_density')
    endif

    if (params%DorE.eq.'N') then
      return
    end if

    call pass_HS(negf,HH,SS)

    call compute_density_dft(negf)

    call destroy_matrices(negf)

  end subroutine negf_density


  !> Interface subroutine to call ldos computation
  subroutine negf_ldos(HH, SS, spin, kpoint, wght, ledos)

    !> hamiltonian in CSR format
    type(z_CSR), pointer, intent(in) :: HH

    !> overlap  in CSR format
    type(z_CSR), pointer, intent(in) :: SS

    !> spin index
    integer, intent(in) :: spin

    !> kpoint index
    integer, intent(in) :: kpoint

    !> kpoint weight                
    real(dp), intent(in) :: wght

    !> local DOS
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


  !> Debug routine to dump H and S as a file in Matlab format
  !>
  !> NOTE: This routine is not MPI-aware, call it only on MPI-lead!
  !>
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


  !> Routines to setup orthogonalised H and S have been moved here
  subroutine prepare_HS(H_dev,S_dev,HH,SS)

    !> hamiltonian in dense format
    real(dp), dimension(:,:) :: H_dev

    !> overlap in dense format
    real(dp), dimension(:,:) :: S_dev

    !> hamiltonian in CSR format
    type(z_CSR), intent(inout) :: HH

    !> overlap in CSR format
    type(z_CSR), intent(inout) :: SS

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


  !> Interface subroutine to call calculation of currents
  subroutine negf_current(HH, SS, spin, kpoint, wght, tunn, curr, ledos, currents)

    !> hamiltonian
    type(z_CSR), pointer, intent(in) :: HH

    !> overlap
    type(z_CSR), pointer, intent(in) :: SS

    !> spin index
    integer, intent(in) :: spin

    !> kpoint index
    integer, intent(in) :: kpoint

    !> kpoint weight
    real(dp), intent(in) :: wght

    !> Tunneling amplitudes
    real(dp), dimension(:,:), pointer :: tunn

    !> current magnitudes
    real(dp), dimension(:,:), pointer :: curr

    !> local density of states
    real(dp), dimension(:,:), pointer :: ledos

    !> current directions
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
  subroutine calcdensity_green(iSCCIter, env, groupKS, ham, over, iNeighbor, nNeighbor, iAtomStart,&
      & iPair, img2CentCell, iCellVec, cellVec, orb, kPoints, kWeights, mu, rho, Eband, Ef, E0, TS)

    !> SCC iteration
    integer, intent(in) :: iSCCIter

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> kpoint and spin descriptor
    integer, intent(in) :: groupKS(:,:)

    !> hamiltonian matrix
    real(dp), intent(in) :: ham(:,:)

    !> overlap matrix
    real(dp), intent(in) :: over(:)

    !> neighbours of atoms
    integer, intent(in) :: iNeighbor(0:,:)

    !> number of neighbours
    integer, intent(in) :: nNeighbor(:)

    !> dense indexing for orbitals
    integer, intent(in) :: iAtomStart(:)

    !> sparse indexing for orbitals
    integer, intent(in) :: iPair(0:,:)

    !> map from image to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> index for unit cell an atom is associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors to unit cells
    real(dp), intent(in) :: cellVec(:,:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> k-points
    real(dp), intent(in) :: kPoints(:,:)

    !> k-point weights
    real(dp), intent(in) :: kWeights(:)

    !> chemical potentials of reservoirs
    real(dp), intent(in) :: mu(:,:)

    !> density matrix
    real(dp), intent(out) :: rho(:,:)

    !> band energy
    real(dp), intent(out) :: Eband(:)

    !> Fermi energy
    real(dp), intent(out) :: Ef(:)

    !> zero temperature (extrapolated) electronic energy
    real(dp), intent(out) :: E0(:)

    !> Electron entropy
    real(dp), intent(out) :: TS(:)

    integer :: nSpin, nKS, iK, iS, iKS
    type(z_CSR), target :: csrDens
    type(z_CSR), pointer :: pCsrDens

    pCsrDens => csrDens

#:if WITH_MPI
    call negf_mpi_init(env%mpi%groupComm)
#:endif
    !Decide what to do with surface GFs.
    !sets readOldSGF: if it is 0 or 1 it is left so
    if (negf%readOldDM_SGFs.eq.COMPSAVE_SGF) then
      if(iSCCIter.eq.1) then
        call set_readOldDMsgf(negf, COMPSAVE_SGF)  ! compute and save SGF on files
      else
        call set_readOldDMsgf(negf, READ_SGF)  ! read from files
      endif
    endif
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

    #:if WITH_MPI
      if (env%mpi%nGroup == 1) then
        write(stdOut,*) 'k-point',iK,'Spin',iS
      end if
    #:else
      write(stdOut,*) 'k-point',iK,'Spin',iS
    #:endif

      call foldToCSR(csrHam, ham(:,iS), kPoints(:,iK), iAtomStart, iPair, iNeighbor, nNeighbor,&
          & img2CentCell, iCellVec, cellVec, orb)
      call foldToCSR(csrOver, over, kPoints(:,ik), iAtomStart, iPair, iNeighbor, nNeighbor,&
          & img2CentCell, iCellVec, cellVec, orb)

      call negf_density(iSCCIter, iS, iK, pCsrHam, pCsrOver, mu(:,iS), DensMat=pCsrDens)

      ! NOTE:
      ! unfold adds up to rho the csrDens(k) contribution
      call unfoldFromCSR(rho(:,iS), csrDens, kPoints(:,iK), kWeights(iK), iAtomStart, iPair,&
          & iNeighbor, nNeighbor, img2CentCell, iCellVec, cellVec, orb)

      call destruct(csrDens)

      ! Set some fake energies:
      Eband(iS) = 0.0_dp
      Ef(iS) = mu(1,iS)
      TS(iS) = 0.0_dp
      E0(iS) = 0.0_dp

    end do

#:if WITH_MPI
    do iS = 1, nSpin
      ! In place all-reduce of the density matrix
      call mpifx_allreduceip(env%mpi%groupComm, rho(:,iS), MPI_SUM)
    end do
    call mpifx_allreduceip(env%mpi%interGroupComm, rho, MPI_SUM)
#:endif

    ! Now SGFs can be read unless not stored 
    if (negf%readOldDM_SGFs.ne.COMP_SGF) then
      call set_readOldDMsgf(negf, READ_SGF)  ! read from files
    end if

    write(stdOut,'(80("="))')
    write(stdOut,*)

  end subroutine calcdensity_green


  !> Calculates energy-weighted density matrix with Green's functions
  subroutine calcEdensity_green(iSCCIter, env, groupKS, ham, over, iNeighbor, nNeighbor,&
      & iAtomStart, iPair, img2CentCell, iCellVec, cellVec, orb, kPoints, kWeights, mu, rhoE)

    !> SCC iteration
    integer, intent(in) :: iSCCIter

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> kpoint and spin descriptor
    integer, intent(in) :: groupKS(:,:)

    !> hamiltonian matrix
    real(dp), intent(in) :: ham(:,:)

    !> overlap matrix
    real(dp), intent(in) :: over(:)

    !> neighbours of atoms
    integer, intent(in) :: iNeighbor(0:,:)

    !> number of neighbours
    integer, intent(in) :: nNeighbor(:)

    !> dense indexing for orbitals
    integer, intent(in) :: iAtomStart(:)

    !> sparse indexing for orbitals
    integer, intent(in) :: iPair(0:,:)

    !> map from image to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> index for unit cell an atom is associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors to unit cells
    real(dp), intent(in) :: cellVec(:,:)

    !> atomic orbital information Needs only orb%nOrbAtom, orb%mOrb
    type(TOrbitals), intent(in) :: orb

    !> k-points
    real(dp), intent(in) :: kPoints(:,:)

    !> k-point weights
    real(dp), intent(in) :: kWeights(:)

    !> chemical potentials
    real(dp), intent(in) :: mu(:,:)

    !> Energy weighted density matrix
    real(dp), intent(out) :: rhoE(:)

    integer :: nSpin, nKS, iK, iS, iKS
    type(z_CSR), target :: csrEDens
    type(z_CSR), pointer :: pCsrEDens

    pCsrEDens => csrEDens

#:if WITH_MPI
    call negf_mpi_init(env%mpi%groupComm)
#:endif
    !Decide what to do with surface GFs.
    !sets readOldSGF: if it is 0 or 1 it is left so
    if (negf%readOldDM_SGFs.eq.COMPSAVE_SGF) then
      if(iSCCIter.eq.1) then
        call set_readOldDMsgf(negf, COMPSAVE_SGF)  ! compute and save SGF on files
      else
        call set_readOldDMsgf(negf, READ_SGF)  ! read from files
      endif
    endif

    ! We need this now for different Fermi energies in colinear spin
    ! Note: the spin polarised does not work with
    ! built-int potentials (the unpolarized does) in the Poisson
    ! I do not set the Fermi because it seems that in libnegf it is
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

      call negf_density(iSCCIter, iS, iK, pCsrHam, pCsrOver, mu(:,iS), EnMat=pCsrEDens)

      ! NOTE:
      ! unfold adds up to rhoEPrim the csrEDens(k) contribution
      !
      call unfoldFromCSR(rhoE, csrEDens, kPoints(:,iK), kWeights(iK), iAtomStart, iPair, iNeighbor,&
          & nNeighbor, img2CentCell, iCellVec, cellVec, orb)

      call destruct(csrEDens)

    end do

    ! In place all-reduce of the energy-weighted density matrix
#:if WITH_MPI
    call mpifx_allreduceip(env%mpi%groupComm, rhoE, MPI_SUM)
    call mpifx_allreduceip(env%mpi%interGroupComm, rhoE, MPI_SUM)
#:endif

    ! Now SGFs can be read unless not stored 
    if (negf%readOldDM_SGFs.ne.COMP_SGF) then
      call set_readOldDMsgf(negf, READ_SGF)  ! read from files
    end if

    write(stdOut,'(80("="))')
    write(stdOut,*)

  end subroutine calcEdensity_green


  !> Calculate the current and optionally density of states
  subroutine calc_current(env, groupKS, ham, over, iNeighbor, nNeighbor, iAtomStart, iPair,&
      & img2CentCell, iCellVec, cellVec, orb, kPoints, kWeights, tunnMat, currMat, ldosMat,&
      & currLead, writeTunn, tWriteLDOS, regionLabelLDOS, mu)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> kpoint and spin descriptor
    integer, intent(in) :: groupKS(:,:)

    !> hamiltonian matrix
    real(dp), intent(in) :: ham(:,:)

    !> overlap matrix
    real(dp), intent(in) :: over(:)

    !> neighbours of atoms
    integer, intent(in) :: iNeighbor(0:,:)

    !> number of neighbours
    integer, intent(in) :: nNeighbor(:)

    !> dense indexing for orbitals
    integer, intent(in) :: iAtomStart(:)

    !> sparse indexing for orbitals
    integer, intent(in) :: iPair(0:,:)

    !> map from image to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> index for unit cell an atom is associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors to unit cells
    real(dp), intent(in) :: cellVec(:,:)

    !> atomic orbital information Needs only orb%nOrbAtom, orb%mOrb
    type(TOrbitals), intent(in) :: orb

    !> k-points
    real(dp), intent(in) :: kPoints(:,:)

    !> k-point weights
    real(dp), intent(in) :: kWeights(:)

    !> matrix of tunnelling amplitudes at each energy from contacts
    real(dp), allocatable, intent(inout) :: tunnMat(:,:)

    !> matrix of current at each energy from currents
    real(dp), allocatable, intent(inout) :: currMat(:,:)

    !> density of states for each energy and region of projection
    real(dp), allocatable, intent(inout) :: ldosMat(:,:)

    !> current into/out of contacts
    real(dp), allocatable, intent(inout) :: currLead(:)

    !> should tunneling data be written
    logical, intent(in) :: writeTunn

    !> should DOS data be written
    logical, intent(in) :: tWriteLDOS

    !> labels for DOS projected regions
    character(lc), allocatable, intent(in) :: regionLabelLDOS(:)

    !> We need this now for different fermi levels in colinear spin
    !> Note: the spin polarised does not work with
    !> built-int potentials (the unpolarized does) in the poisson
    !> I do not set the fermi because it seems that in libnegf it is
    !> not really needed
    real(dp), intent(in) :: mu(:,:)

    real(dp), allocatable :: tunnSKRes(:,:,:), currSKRes(:,:,:), ldosSKRes(:,:,:)
    real(dp), pointer    :: tunnPMat(:,:)=>null()
    real(dp), pointer    :: currPMat(:,:)=>null()
    real(dp), pointer    :: ldosPMat(:,:)=>null()
    real(dp), pointer    :: currPVec(:)=>null()
    integer :: iKS, iK, iS, nKS, nS,  nTotKS, ii, err, ncont
    type(units) :: unitOfEnergy        ! Set the units of H
    type(units) :: unitOfCurrent       ! Set desired units for Jel
    type(lnParams) :: params

    integer :: NumStates
    real(dp), dimension(:,:), allocatable :: H_all, S_all
    character(:), allocatable :: filename

#:if WITH_MPI
    call negf_mpi_init(env%mpi%groupComm)
#:endif
    call get_params(negf, params)

    unitOfEnergy%name = "H"
    unitOfCurrent%name = "A"

    ! groupKS is local, hence nKS il local
    nKS = size(groupKS, dim=2)
    nS=size(ham,2)
    if (nS>2) then
      nS=1
    end if
    nTotKS = nS * size(kpoints, dim=2)
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
      ! THIS MAKES SENSE ONLY FOR REAL MATRICES, i.e. k==0 && collinear spin
      if (all(kPoints(:,iK) == 0.0_dp) .and. (negf%tOrthonormal .or. negf%tOrthonormalDevice)) then

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

      call negf_current(pCsrHam, pCsrOver, iS, iK, kWeights(iK), tunnPMat, currPMat, ldosPMat,&
          & currPVec)

      if(.not.allocated(currLead)) then
        allocate(currLead(size(currPVec)), stat=err)
        if (err /= 0) then
          call error('Allocation error (currTot)')
        end if
        currLead(:) = 0.0_dp
      endif
      currLead(:) = currLead + currPVec

      !GUIDE: tunnPMat libNEGF output stores Transmission, T(iE, i->j)
      !       tunnPMat is MPI distributed on energy points (0.0 on other nodes)
      !       tunnMat MPI gather partial results and accumulate k-summation
      !       currPMat stores contact current I_i(iE)
      !       tunnSKRes stores tunneling for all k-points and spin: T(iE, i->j, iSK)

#:if WITH_MPI
      ii = (iS-1)*size(kpoints,2) + iK
      call add_partial_results(env%mpi%groupComm, tunnPMat, tunnMat, tunnSKRes, ii, nTotKS)
      call add_partial_results(env%mpi%groupComm, currPMat, currMat, currSKRes, ii, nTotKS)
      call add_partial_results(env%mpi%groupComm, ldosPMat, ldosMat, ldosSKRes, ii, nTotKS)
#:else
      call add_partial_results(tunnPMat, tunnMat, tunnSKRes, iKS, nKS)
      call add_partial_results(currPMat, currMat, currSKRes, iKS, nKS)
      call add_partial_results(ldosPMat, ldosMat, ldosSKRes, iKS, nKS)
#:endif

    end do

#:if WITH_MPI
    call mpifx_reduceip(env%mpi%groupComm, currLead, MPI_SUM)
    call mpifx_reduceip(env%mpi%interGroupComm, currLead, MPI_SUM)
    call add_ks_results(env%mpi%interGroupComm, tunnMat, tunnSKRes)
    call add_ks_results(env%mpi%interGroupComm, currMat, currSKRes)
    call add_ks_results(env%mpi%interGroupComm, ldosMat, ldosSKRes)
#:endif

    ! converts from internal atomic units into amperes
    currLead(:) = currLead * convertCurrent(unitOfEnergy, unitOfCurrent)

    do ii = 1, size(currLead)
      write(stdOut, *)
      write(stdOut, '(1x,a,i3,i3,a,ES14.5,a,a)') ' contacts: ',params%ni(ii),params%nf(ii),&
          & ' current: ', currLead(ii),' ',unitOfCurrent%name
    enddo

    ! Write Total transmission, T(E), on a separate file (optional)
    if (allocated(tunnMat)) then
      filename = 'transmission'
      if (tIOProc .and. writeTunn) then
        call write_file(negf, tunnMat, tunnSKRes, filename, nS, kpoints, kWeights)
      end if
    else
      ! needed to avoid some segfault
      allocate(tunnMat(0,0))
    end if
    if (allocated(tunnSKRes)) then
      deallocate(tunnSKRes)
    end if

    ! Write Total lead current, I_i(E), on a separate file (optional)
    if (allocated(currMat)) then
      filename = 'current'
      if (tIOProc .and. writeTunn) then
        call write_file(negf, currMat, currSKRes, filename, nS, kpoints, kWeights)
      end if
    else
      ! needed to avoid some segfault
      allocate(currMat(0,0))
    endif
    if (allocated(currSKRes)) then
      deallocate(currSKRes)
    end if

    if (allocated(ldosMat)) then
      ! Write Total localDOS on a separate file (optional)
      if (tIoProc .and. tWriteLDOS) then
    @:ASSERT(allocated(regionLabelLDOS))
        call write_files(negf, ldosMat, ldosSKRes, nS, kpoints, kWeights, regionLabelLDOS)
      end if
    else
        ! needed to avoid some segfault
        allocate(ldosMat(0,0))
    end if
    if (allocated(ldosSKRes)) then
      deallocate(ldosSKRes)
    end if
    
  end subroutine calc_current


  !> utility to allocate and sum partial results from different channels
#:if WITH_MPI
  subroutine add_partial_results(mpicomm, pMat, matTot, matSKRes, iK, nK)

    !> MPI communicator
    type(mpifx_comm), intent(in) :: mpicomm
#:else
  subroutine add_partial_results(pMat, matTot, matSKRes, iK, nK)
#:endif

    !> pointer to matrix of data
    real(dp), intent(in), pointer :: pMat(:,:)

    !> sum total
    real(dp), allocatable, intent(inout) :: matTot(:,:)

    !> k-resolved sum
    real(dp), allocatable, intent(inout)  :: matSKRes(:,:,:)

    !> particular k-point
    integer, intent(in) :: iK

    !> number of k-points
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
  !> utility to sum up partial results over SK communicator
  subroutine add_ks_results(kscomm, mat, matSKRes)

    !> MPI communicator
    type(mpifx_comm), intent(in) :: kscomm

    !> sum total
    real(dp), allocatable, intent(inout) :: mat(:,:)

    !> k-resolved sum
    real(dp), allocatable, intent(inout)  :: matSKRes(:,:,:)

    if (allocated(mat)) then
      call mpifx_reduceip(kscomm, mat, MPI_SUM)
    endif

    if (allocated(matSKRes)) then
      call mpifx_reduceip(kscomm, matSKRes, MPI_SUM)
    endif

  end subroutine add_ks_results
#:endif


  !> utility to write tunneling files
  subroutine write_file(negf, matTot, matSKRes, filename, nS, kpoints, kWeights)

    !> Contains input data, runtime quantities and output data
    type(TNegf) :: negf

    !> results to print if allocated
    real(dp), intent(in), allocatable :: matTot(:,:)

    !> k- and spin-resolved quantities, if allocated
    real(dp), intent(in), allocatable :: matSKRes(:,:,:)

    !> file to print out to
    character(*), intent(in) :: filename

    !> number of spins 
    integer, intent(in) :: nS

    !> k-points
    real(dp), intent(in) :: kPoints(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeights(:)

    integer :: ii, jj, nK, iK, iS, iKS, fdUnit
    type(lnParams) :: params

    call get_params(negf, params)

    nK = size(kpoints, dim=2)

    open(newunit=fdUnit, file=trim(filename)//'.dat')
    do ii = 1, size(matTot, dim=1)
      write(fdUnit,'(F20.6)',ADVANCE='NO') (params%Emin+(ii-1)*params%Estep) * Hartree__eV
      do jj = 1, size(matTot, dim=2)
        write(fdUnit,'(ES20.8)',ADVANCE='NO') matTot(ii,jj)
      enddo
      write(fdUnit,*)
    enddo
    close(fdUnit)

    if (nK*nS > 1) then

      open(newunit=fdUnit, file=trim(filename)//'_kpoints.dat')
      write(fdUnit,*)'# NKpoints = ', nK
      write(fdUnit,*)'# NSpin = ', nS
      write(fdUnit,*)'# Energy [eV], <spin k1 k2 k3 weight> '
      write(fdUnit,'(A1)', ADVANCE='NO') '# '

      ! iKS = 1 2 3 4 5 6 7 8 9 10
      ! iK=groupKS(1,iKS), iS=groupKS(2,iKS)
      ! iK  = 1 2 3 4 5 1 2 3 4 5  
      ! iS  = 1 1 1 1 1 2 2 2 2 2 
      do iS = 1, nS
        do iK = 1, nK
          write(fdUnit,'(i5.2)', ADVANCE='NO') iS
          write(fdUnit,'(es15.5, es15.5, es15.5, es15.5)', ADVANCE='NO') kpoints(:,iK), kWeights(iK)
        end do
      end do  
      write(fdUnit,*)

      if (allocated(matSKRes)) then
        do ii = 1, size(matSKRes(:,:,:), dim=1)
          write(fdUnit,'(f20.6)',ADVANCE='NO') (params%Emin+(ii-1)*params%Estep) * Hartree__eV
          do jj = 1, size(matSKRes(:,:,:), dim=2)
            do iKS = 1, nK*nS
              write(fdUnit,'(es20.8)',ADVANCE='NO') matSKRes(ii,jj, iKS)
            enddo
            write(fdUnit,*)
          enddo
        enddo
      end if
      close(fdUnit)

    end if

  end subroutine write_file


  !> utility to write tunneling/ldos files with names from labels
  subroutine write_files(negf, matTot, matSKRes, nS, kpoints, kWeights, regionLabels)

    !> Contains input data, runtime quantities and output data
    type(TNegf) :: negf

    !> results to print if allocated
    real(dp), intent(in) :: matTot(:,:)

    !> k- and spin-resolved quantities, if allocated
    real(dp), intent(in), allocatable :: matSKRes(:,:,:)

    !> number of spins
    integer, intent(in) :: nS

    !> k-points
    real(dp), intent(in) :: kPoints(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeights(:)

    !> Labels for the separate regions
    character(lc), intent(in) :: regionLabels(:)

    integer :: ii, jj, nKS, iKS, nK, iK, iS, fdUnit
    type(lnParams) :: params

    call get_params(negf, params)

    nK = size(kpoints, dim=2)
    nKS = nK*nS

    do jj=1,size(matTot, dim=2)
      open(newunit=fdUnit, file=trim(regionLabels(jj))//'.dat')
      write(fdUnit,"(A)")'# Energy / eV     States / e'
      do ii=1,size(matTot, dim=1)
        write(fdUnit,'(F12.6)',ADVANCE='NO') (params%Emin+(ii-1)*params%Estep) * Hartree__eV
        write(fdUnit,'(ES20.8)') matTot(ii,jj)
      enddo
      close(fdUnit)
    enddo

    if (allocated(matSKRes)) then
      if (nKS > 1) then
        do jj = 1, size(matSKRes(:,:,:), dim=2)
          open(newunit=fdUnit, file=trim(regionLabels(jj))//'_kpoints.dat')
          write(fdUnit,"(A,I0)")'# NKpoints = ', nK
          write(fdUnit,"(A,I1)")'# NSpin = ', nS
          write(fdUnit,"(A)")'# <spin k1 k2 k3 weight> '
          write(fdUnit,'(A1)', ADVANCE='NO') '# '
          do iS = 1, nS
            do iK = 1, nK
              write(fdUnit,'(i5.1)', ADVANCE='NO') iS
              write(fdUnit,'(es15.5, es15.5, es15.5, es15.5)', ADVANCE='NO') kpoints(:,iK),&
                  & kWeights(iK)
            end do
          end do
          write(fdUnit,*)

          do ii = 1, size(matSKRes(:,:,:), dim=1)
            write(fdUnit,'(f20.6)',ADVANCE='NO') (params%Emin+(ii-1)*params%Estep) * Hartree__eV
            do iKS = 1,nKS
              write(fdUnit,'(es20.8)',ADVANCE='NO') matSKRes(ii,jj, iKS)
            enddo
            write(fdUnit,*)
          enddo
          close(fdUnit)
        enddo
      end if
    end if

  end subroutine write_files


  !> THIS is a first version of local current computation.
  ! It has been placed here since it depends on internal representations of DFTB
  !
  ! NOTE: Limited to non-periodic systems
  subroutine local_currents(env, groupKS, ham, over, neighbourList, nNeighbour, skCutoff,&
      & iAtomStart, iPair, img2CentCell, iCellVec, cellVec, rCellVec, orb, kPoints, kWeights,&
      & coord0, species0, speciesName, chempot, testArray)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> kpoint and spin descriptor
    integer, intent(in) :: groupKS(:,:)

    !> Hamiltonian and Overlap matrices 
    real(dp), intent(in) :: ham(:,:), over(:)

    !> Neighbor list container
    type(TNeighbourList), intent(in) :: neighbourList

    !> nNeighbor list for SK interactions
    integer, intent(in) :: nNeighbour(:)

    !> SK interaction cutoff
    real(dp), intent(in) :: skCutoff

    !> Indices of staring atom block and Pairs 
    integer, intent(in) :: iAtomStart(:), iPair(0:,:)

    !> map of atoms to central cell
    integer, intent(in) :: img2CentCell(:), iCellVec(:)

    !> translation vectors to lattice cells in units of lattice constants
    real(dp), allocatable, intent(in) :: cellVec(:,:)

    !> Vectors to unit cells in absolute units
    real(dp), allocatable, intent(in) :: rCellVec(:,:)

    !> Orbital descriptor 
    type(TOrbitals), intent(in) :: orb

    !> k-points and weights 
    real(dp), intent(in) :: kPoints(:,:), kWeights(:)

    !> central cell coordinates (folded to central cell) 
    real(dp), intent(in) :: coord0(:,:)

    !> Species indices for atoms in central cell
    integer, intent(in) :: species0(:)

    !> Species Names (as in gen file)
    character(*), intent(in) :: speciesName(:)

    ! We need this now for different fermi levels in colinear spin
    ! Note: spin polarized does not work with
    ! built-in potential (the unpolarized does) in the poisson
    ! I do not set the fermi because it seems that in libnegf it is
    ! not really needed
    real(dp), intent(in) :: chempot(:,:)

    !> Array passed back to main for autotests (will become the output)
    real(dp), allocatable, intent(out) :: testArray(:,:)


    ! Local stuff ---------------------------------------------------------
    integer :: n0, nn, mm,  mu, nu, nAtom, irow
    integer :: nKS, nK, nSpin, iKS, iK, iS, iKgl, inn, startn, endn, morb
    real(dp), dimension(:,:,:), allocatable :: lcurr 
    real(dp) :: Im
    type(TNeighbourList) :: lc_neigh
    integer, dimension(:), allocatable :: lc_img2CentCell, lc_iCellVec, lc_species
    real(dp), dimension(:,:), allocatable :: lc_coord 
    integer :: lc_nAllAtom
    integer, parameter :: nInitNeigh=40
    complex(dp) :: c1,c2
    character(:), allocatable :: skp
    character(6) :: fmtstring
    integer :: iSCCiter
    type(z_CSR), target :: csrDens, csrEDens
    type(z_CSR), pointer :: pCsrDens, pCsrEDens
    type(lnParams) :: params
    integer :: fdUnit
    logical :: tPrint

    pCsrDens => csrDens
    pCsrEDens => csrEDens

#:if WITH_MPI
    call negf_mpi_init(env%mpi%groupComm)
#:endif
    call get_params(negf, params)

    !Decide what to do with surface GFs.
    !sets readOldSGF: if it is 0 or 1 it is left so
    if (negf%readOldDM_SGFs.eq.COMPSAVE_SGF) then
      if(iSCCIter.eq.1) then
        call set_readOldDMsgf(negf, COMPSAVE_SGF)  ! compute and save SGF on files
      else
        call set_readOldDMsgf(negf, READ_SGF)  ! read from files
      endif
    endif

    write(stdOut, *)
    write(stdOut, '(80("="))')
    write(stdOut, *) '                        COMPUTING LOCAL CURRENTS          '
    write(stdOut, '(80("="))')
    write(stdOut, *)

    nKS = size(groupKS, dim=2)
    nK = size(kPoints, dim=2)
    nSpin = size(ham, dim=2)
    nAtom = size(orb%nOrbAtom)
    call get_fmtstring(nK, skp, fmtstring)

    ! Create a symmetrized neighbour list extended to periodic cell in lc_coord 
    if (any(iCellVec.ne.1)) then
      lc_nAllAtom = int((real(nAtom, dp)**(1.0_dp/3.0_dp) + 3.0_dp)**3)
    else
      lc_nAllAtom = nAtom
    end if
    allocate(lc_coord(3, lc_nAllAtom))
    allocate(lc_species(lc_nAllAtom))
    allocate(lc_img2CentCell(lc_nAllAtom))
    allocate(lc_iCellVec(lc_nAllAtom))
    call init(lc_neigh, nAtom, nInitNeigh) 

    call updateNeighbourListAndSpecies(lc_coord, lc_species, lc_img2CentCell, lc_iCellVec, &
        & lc_neigh, lc_nAllAtom, coord0, species0, skCutoff, rCellVec, symmetric=.true.)

    allocate(lcurr(maxval(lc_neigh%nNeighbour), nAtom, nSpin))
    lcurr(:,:,:) = 0.0_dp

    looppKS: do iKS = 1, nKS
      iK = groupKS(1, iKS)
      iS = groupKS(2, iKS)

      write(stdOut,*) 'k-point',iK,'Spin',iS

      ! We need to recompute Rho and RhoE .....
      call foldToCSR(csrHam, ham(:,iS), kPoints(:,iK), iAtomStart, iPair, neighbourList%iNeighbour,&
          & nNeighbour, img2CentCell, iCellVec, CellVec, orb)
      call foldToCSR(csrOver, over, kPoints(:,iK), iAtomStart, iPair, neighbourList%iNeighbour, &
          & nNeighbour, img2CentCell, iCellVec, CellVec, orb)

      call negf_density(iSCCIter, iS, iK, pCsrHam, pCsrOver, chempot(:,iS), DensMat=pCsrDens)

      ! Unless SGFs are not stored, read them from file
      if (negf%readOldDM_SGFs.ne.COMP_SGF) then
         call set_readOldDMsgf(negf, READ_SGF) 
      end if   

      call negf_density(iSCCIter, iS, iK, pCsrHam, pCsrOver, chempot(:,iS), EnMat=pCsrEDens)

    #:if WITH_MPI
      ! Reduce on node 0 as group lead node
      call mpifx_reduceip(env%mpi%groupComm, csrDens%nzval, MPI_SUM)
      call mpifx_reduceip(env%mpi%groupComm, csrEDens%nzval, MPI_SUM)

      ! Each group lead node prints the local currents
      tPrint = env%mpi%groupComm%lead

    #:else

      tPrint = .true.

    #:endif

      if (tPrint) then
        ! print local currents
        iKgl = (iS-1) * nK + iK    
        write(skp, fmtstring) iKgl
        open(newUnit = fdUnit, file = 'lcurrents_'//skp//"_"//spin2ch(iS)//'.dat')

        ! loop on central cell atoms and write local currents to all other 
        ! interacting atoms within the cell and neighbour cells
        do mm = 1, nAtom

          mOrb = orb%nOrbAtom(mm)
          iRow = iAtomStart(mm)

          write(fdUnit,'(I5,3(F12.6),I4)',advance='NO') mm, lc_coord(:,mm), lc_neigh%nNeighbour(mm)

          do inn = 1, lc_neigh%nNeighbour(mm)
            nn = lc_neigh%iNeighbour(inn, mm) 
            n0 = lc_img2CentCell(nn)
            startn = iAtomStart(n0)
            endn = startn + orb%nOrbAtom(n0) - 1
            Im = 0.0_dp
            ! tracing orbitals of atoms  n  m
            ! More efficient without getel ?
            do mu = iRow, iRow+mOrb-1
              do nu = startn, endn
                c1 = conjg(getel(csrDens,mu,nu))
                c2 = conjg(getel(csrEDens,mu,nu))
                Im = Im + aimag(getel(csrHam,mu,nu)*c1 - getel(csrOver,mu,nu)*c2)
              end do
            end do
            ! pi-factor  comes from  Gn = rho * pi
            Im = Im * 2.0_dp*params%g_spin*pi*eovh*kWeights(iK)
            write(fdUnit,'(I5,ES17.8)',advance='NO') nn, Im 
            lcurr(inn, mm, iS) = lcurr(inn, mm, iS) + Im
          end do

          write(fdUnit,*)
        end do

        close(fdUnit)
      end if

      call destruct(csrDens)
      call destruct(csrEDens)

    end do looppKS

#:if WITH_MPI
    call mpifx_reduceip(env%mpi%interGroupComm, lcurr, MPI_SUM)
#:endif

    if (tIoProc) then
      allocate(testArray(maxval(lc_neigh%nNeighbour),nAtom*nSpin))
      testArray(:,:) = 0.0_dp
      ! Write the total current per spin channel  
      do iS = 1, nSpin
        open(newUnit = fdUnit, file = 'lcurrents_'//spin2ch(iS)//'.dat')
        do mm = 1, nAtom
          write(fdUnit,'(I5,3(F12.6),I4)',advance='NO') mm, lc_coord(:,mm), lc_neigh%nNeighbour(mm)
          do inn = 1, lc_neigh%nNeighbour(mm)
            write(fdUnit,'(I5,ES17.8)',advance='NO') lc_neigh%iNeighbour(inn, mm), lcurr(inn,mm,iS)
            testArray(inn,(iS-1)*nAtom+mm) = lcurr(inn,mm,iS) 
          end do
          write(fdUnit,*)
        end do
      end do
      close(fdUnit)
    end if
    deallocate(lcurr)

    if (tIoProc) then
      write(stdOut,*) 
      call writeXYZFormat("supercell.xyz", lc_coord, lc_species, speciesName)
      write(stdOut,*) " <<< supercell.xyz written on file"
    end if

  contains

    !> labels from spin channel number
    pure function spin2ch(iS) result(ch)

      !> channel number
      integer, intent(in) :: iS

      !> label
      character(1) :: ch

      character(1), parameter :: labels(2) = ['u', 'd']

      ch = labels(iS)

    end function spin2ch
    
    subroutine get_fmtstring(nK, skp, fmtstring)
      integer, intent(in) :: nK
      character(:), allocatable :: skp
      character(6) :: fmtstring
      integer :: nchars
          
      nchars = 3
      do while (nK/(10**nchars) > 1 )
        nchars = nchars + 1
      end do 
      allocate(character(len=nchars)::skp)    
      ! create fmtstring = '(In.n)'
      write(fmtstring, '( "(I",I1,".",I1,")" )') nchars, nchars

    end subroutine get_fmtstring

  end subroutine local_currents


  !> pack dense matrices into CSR format
  subroutine MakeHHSS(H_all, S_all, HH, SS)

    !> hamiltonian matrix
    real(dp), intent(in) :: H_all(:,:)

    !> overlap matrix
    real(dp), intent(in) :: S_all(:,:)

    !> hamiltonian in CSR
    type(z_CSR), intent(inout) :: HH

    !> overlap in CSR
    type(z_CSR), intent(inout) :: SS

    integer :: i, j, k, NumStates, nnz

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


  !> form orthogonal matrices via Lowdin transform for whole system
  subroutine orthogonalization(H,S)

    !> hamiltonian matrix
    real(dp), intent(inout) :: H(:,:)

    !> overlap matrix
    real(dp), intent(inout) :: S(:,:)

    integer :: i, m, n1_first, n1_last, n2_first, n2_last
    integer :: N
    real(dp), allocatable :: A(:,:), W(:)
    real(dp), allocatable :: B(:,:),C(:,:)

    N = size(H, dim=1)

    if (N /= negf%NumStates) then
      call error('orthogonalization: negf init NumStates error')
    end if

    allocate(A(N,N),W(N))
    allocate(B(N,N),C(N,N))
    W(:) = 0.0_dp
    A(:,:) = 0.0_dp
    B(:,:) = 0.0_dp
    C(:,:) = 0.0_dp

    A(:,:) = S

    call heev(A, W, 'L', 'V')
    !call DSYEV('V','U',N,A,N,W,WORK,3*N,INFO )

    !print  *,'U matrix, Eigenvectors for S diagonalization'
    !do i=1,N
    !   write(stdOut,*)A(i,1:N)
    !end do

    !print *,'U matrix unitarity check'
    !B=matmul(transpose(A),A)
    !do i=1,N
    !   write(stdOut,*)B(i,1:N)
    !end do

    B(:,:) = matmul(transpose(A), matmul(S, A))

    do i=1,N
      B(i,i) = 1.0_dp / sqrt(B(i,i))
    end do

    !C=sqrt(S)
    C(:,:) = matmul(A, matmul(B, transpose(A)))

    !print *,'sqrt(S) inverted'
    !do i=1,N
    !  write(stdOut,*) C(i,1:N)
    !end do

    !print *,'S unity check'
    !B=matmul(transpose(C),matmul(S,C))
    !do i=1,N
    !   write(stdOut,*) B(i,1:N)
    !end do

    !print *,'H_dftb before orthogonalization'
    !do i=1,N
    !   write(stdOut,*) H(i,1:N)
    !end do

    H(:,:) = matmul(transpose(C), matmul(H, C))

    !print *,'H_dftb_orth before replacement'
    !do i=1,N
    !   write(stdOut,*) H(i,1:N)
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
    !   write(stdOut,*) H(i,1:N)
    !end do

    S(:,:) = 0.0_dp
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


  !> Orthogonalise basis in device region
  subroutine orthogonalization_dev(H, S)

    !> hamilonian matrix
    real(dp), intent(inout) :: H(:,:)

    !> overlap matrix
    real(dp), intent(inout) :: S(:,:)


    integer :: i, m, n1_first, n1_last, n2_first, n2_last
    integer :: N, N2
    real(dp), allocatable :: A(:,:), W(:)
    real(dp), allocatable :: B(:,:), U(:,:), C(:,:)

    N = size(H, dim=1)
    if (N /= negf%NumStates) then
      call error('orthogonalization: negf init NumStates error')
    end if

    N2=negf%str%central_dim

    allocate(A(N2,N2),W(N2))
    allocate(B(N2,N2), U(N2,N2))
    W(:) = 0.0_dp

    A(:,:) = S(1:N2,1:N2)
    U(:,:) = A

    call heev(U, W, 'L', 'V')
    !call DSYEV('V','U',N2,U,N2,W,WORK,3*N2,INFO )

    !print  *,'U matrix, Eigenvectors for S diagonalization'
    !do i=1,N2
    !   write(stdOut,*) U(i,1:N2)
    !end do

    !U matrix unitarity check
    !B = matmul(transpose(U),U)
    !do i=1,N2
    !  if (abs(B(i,i)-1.0_dp) > 1.d-9 .or. any(abs(B(i,1:i-1)) > 1.d-9 ) then
    !    print*, 'ERROR: U is not unitary', B(i,:)
    !    stop
    !  end if
    !end do

    B(:,:) = matmul(transpose(U),matmul(A,U))

    do i = 1, N2
      B(i,i) = 1.0_dp / sqrt(B(i,i))
    end do


    ! Now A = S^-1/2
    A(:,:) = matmul(U,matmul(B,transpose(U)))
    !print *,'sqrt(S) inverted'
    !do i=1,N2
    !  write(stdOut,*) A(i,1:N2)
    !end do

    deallocate(U, B)

    allocate(C(N,N))
    C(:,:) = 0.0_dp
    do i = N2+1, N
      C(i,i) = 1.0_dp
    end do
    C(1:N2,1:N2) = A

    deallocate(A)

    !C=sqrt(S) big matrix
    !print *,'C=sqrt(S) big matrix'
    !do i=1,N
    !   write(stdOut,*) C(i,1:N)
    !end do

    !print *,'H_dftb before orthogonalization'
    !do i=1,N
    !   write(stdOut,*) H(i,1:N)
    !end do

    H(:,:) = matmul(transpose(C),matmul(H,C))
    S(:,:) = matmul(transpose(C),matmul(S,C))

    !print *,'H_dftb_orth before replacement'
    !do i=1,N
    !   write(stdOut,*) H(i,1:N)
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
    !   write(stdOut,*) H(i,1:N)
    !end do

    !print *,'S_dftb_orth after replacement'
    !do i=1,N
    !   write(stdOut,*) S(i,1:N)
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
