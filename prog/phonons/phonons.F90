!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program phonons 
  use dftbp_assert
  use dftbp_globalenv
  use dftbp_environment
  use dftbp_initphonons
  use dftbp_accuracy, only : dp, lc
  use dftbp_constants, only : Hartree__cm, Bohr__AA, Hartree__J, Hartree__eV, hbar, pi
  use dftbp_simplealgebra, only : invert33
  use dftbp_typegeometry
  use dftbp_eigensolver, only : heev
  use dftbp_message
  use dftbp_taggedoutput
  use libnegf_int
  use ln_structure
  implicit none

  type(TEnvironment) :: env
  logical :: tInitialized
  real(dp), allocatable :: tunnTot(:,:), ldosTot(:,:), conductance(:,:)
  real(dp), allocatable :: currLead(:)
  logical :: twriteLDOS
  logical :: twriteTunn
  type (TTaggedWriter) :: taggedWriter
  integer :: err, nProcs

  call initGlobalEnv()
  call printHeader()
  call TEnvironment_init(env)
  call initProgramVariables(env)
  call TTaggedWriter_init(taggedWriter)

  nProcs = 1
  write(stdOut,*)'Computing environment'
#:if WITH_MPI
  print*,env%mpi%globalComm%rank, env%mpi%globalComm%size, tIOProc, env%mpi%globalComm%lead
  nProcs = env%mpi%globalComm%size
#:else
  write(stdOut,*)'Not compiled with MPI enabled'
#:endif

  if (tCompModes) then
    if (nProcs > 1) then
      call error("Mode calculation is not parallel yet. Run just on 1 node")
      call destructProgramVariables()
    end if
    call ComputeModes()
  end if

  if (tPhonDispersion) then
    if (nProcs > 1) then
      call error("Phonon dispersion is not parallel yet. Run just on 1 node")
      call destructProgramVariables()
    end if
    call PhononDispersion(taggedWriter)
  end if

  if (tTransport) then
    twriteTunn = .true.
    twriteLDOS = .true.
    call negf_init(env, transpar, tundos, tInitialized)
    call init_tun_proj(selTypeModes, geo%nAtom)

    if (.not. tInitialized) then
      call error("libnegf not initialized")
    end if

    call negf_init_str(geo%nAtom, transpar, neighbourList%iNeighbour, nNeighbour, img2CentCell)
   
    call calc_phonon_current(env, dynMatrix, tunnTot, ldosTot, currLead, conductance, &
                        & twriteTunn, twriteLDOS)  

    if (tWriteTagged) then              
      call writeTaggedOut(taggedWriter, tunnTot, ldosTot, conductance)      
    end if

  end if

  call destructProgramVariables()
  call destructGlobalEnv()

  write(stdOut, "(A)") "Program completed successfully"

contains
    
  subroutine printHeader()
    write (stdOut, "(A)") repeat("=", 80)
    write (stdOut, "(30X,A)") "PHONONS" 
    write (stdOut, "(A,/)") repeat("=", 80)
    write (stdOut, "(A)") "" 
    write (stdOut, "(A)") "Version 0.1" 
    write (stdOut, "(A)") "A tool to compute phonon transmission in nanostructures based on Hessians" 
    write (stdOut, "(A)") "Authors: Alessandro Pecchia, Leonardo Medrano Sandonas" 
    write (stdOut, "(A)") "When using this code, please cite this work:" 
    write (stdOut, "(A)") "Leonardo Medrano Sandonas, Rafaei Gutierrez, Alessandro Pecchia,"
    write (stdOut, "(A)") "Alexander Croy, Gianaurelio Cuniberti, Quantum phonon transport in"
    write (stdOut, "(A)") "nanomaterials: combining atomistic with non-equilibrium Green's functions"
    write (stdOut, "(A)") "techniques, Entropy 21, 735 (2019)" 
    write (stdOut, "(A)") "" 
  end subroutine printHeader

  subroutine ComputeModes()

    integer  :: ii, jj, kk, ll, nAtom, iMode, jCount, iAt, iAtMoved, fu
    real(dp), allocatable :: eigenValues(:)
    real(dp), allocatable :: displ(:,:,:)
    character(lc) :: lcTmp, lcTmp2
    
    nAtom = geo%nAtom
 
    allocate(eigenValues(3*nMovedAtom))
 
    ! solve the eigenproblem
    if (tPlotModes) then
      write(stdOut,*) 'Computing vibrational frequencies and modes'
      call heev(dynMatrix,eigenValues,'U','V')
    else
      write(stdOut,*) 'Computing vibrational frequencies'
      call heev(dynMatrix,eigenValues,'U','N')
    end if
 
    ! take square root of modes (allowing for imaginary modes) and print
    eigenValues =  sign(sqrt(abs(eigenValues)),eigenValues)
    write(stdOut,*)'Vibrational modes (cm-1):'
    do ii = 1, 3 * nMovedAtom
      write(stdOut,'(f8.2)')eigenValues(ii)*Hartree__cm
    end do
    write(stdOut,*)
 
    if (tPlotModes) then
      write(stdOut,*)'Plotting eigenmodes:'
      write(stdOut,*)ModesToPlot(:)
      ! scale mode compoents on each atom by mass and then normalise total mode
      do ii = 1, nModesToPlot
        iMode = ModesToPlot(ii)
        jCount = 0
        do jj = 1, nMovedAtom
          do ll = 1, 3
            jCount = jCount + 1
            dynMatrix(jCount,iMode) = dynMatrix(jCount,iMode) &
                & /sqrt(atomicMasses(jj))
          end do
        end do
        dynMatrix(:,iMode) = dynMatrix(:,iMode) &
            & / sqrt(sum(dynMatrix(:,iMode)**2))
      end do
 
     ! Create displacment vectors for every atom in every mode.
      allocate(displ(3, nAtom, nModesToPlot))
      displ(:,:,:) = 0.0_dp
      do iAt = 1, nAtom
        if (any(iMovedAtoms == iAt)) then
          ! Index of atom in the list of moved atoms
          iAtMoved = minloc(abs(iMovedAtoms - iAt), 1)
          do ii = 1, nModesToPlot
            iMode = ModesToPlot(ii)
            displ(:,iAt, ii) =  dynMatrix(3*iAtMoved-2:3*iAtMoved, iMode)
          end do
        end if
      end do
 
 
      if (tAnimateModes) then
        do ii = 1, nModesToPlot
          iMode = ModesToPlot(ii)
          write(lcTmp,"('mode_',I0)")iMode
          write(lcTmp2, "(A,A)") trim(lcTmp), ".xyz"
          open(newunit=fu, file=trim(lcTmp2), position="rewind", status="replace")
          do kk = 1, nCycles
            do ll = 1, nSteps
              write(fu,*)nAtom
              write(fu,*)'Eigenmode',iMode,eigenValues(iMode)*Hartree__cm,'cm-1'
              do iAt = 1, nAtom
                write(fu,'(A3,T4,3F10.6)') &
                    & geo%speciesNames(geo%species(iAt)), &
                    & (geo%coords(:,iAt)&
                    & + cos(2.0_dp * pi * real(ll) / real(nSteps))&
                    & * displ(:,iAt,ii)) * Bohr__AA
              end do
            end do
          end do
          close(fu)
        end do
      else
        open(fu, file="modes.xyz", position="rewind", status="replace")
        do ii = 1, nModesToPlot
          iMode = ModesToPlot(ii)
          write(fu,*)nAtom
          write(fu,*)'Eigenmode',iMode,eigenValues(iMode)*Hartree__cm,'cm-1'
          if (tXmakeMol) then ! need to account for its non-standard xyz vector
            ! format:
            do iAt = 1, nAtom
              write(fu,'(A3,T4,3F10.6,A,3F10.6)') &
                  & geo%speciesNames(geo%species(iAt)), &
                  & geo%coords(:,iAt)* Bohr__AA, ' atom_vector ',&
                  & displ(:,iAt,ii)
            end do
          else ! genuine xyz format
            do iAt = 1, nAtom
              write(fu,'(A3,T4,6F10.6)') &
                  & geo%speciesNames(geo%species(iAt)), &
                  & geo%coords(:,iAt)* Bohr__AA, &
                  & displ(:,iAt,ii)
            end do
          end if
        end do
        close(fu)
      end if
      
    end if

    deallocate(eigenValues)

  end subroutine ComputeModes

  subroutine PhononDispersion(tWriter)
    type(TTaggedWriter) :: tWriter

    integer  :: ii, jj, kk, ll, nAtom,  iAtom,  iK, jAtom,  kAtom
    integer  :: i2, j2, k2, fu, ftag, nrep
    real(dp), allocatable :: eigenValues(:)
    real(dp)::  ModKPoint,  ModDeltaR
    character(lc) :: lcTmp, lcTmp2
    complex(dp), dimension(:,:), allocatable :: KdynMatrix
    real(dp) :: latVecs(3,3), invLatt(3,3) 
    real(dp) :: DeltaR(3), q(3), qold(3)
    complex(dp), parameter ::    j = (0.d0,1.d0)
    real(dp) :: unitsConv
      
    call setConversionUnits(unitsConv)

    write(stdOut,*) 'Supercell repetitions:' 
    write(stdOut,*) nCells(1),'x',nCells(2),'x',nCells(3) 

    latVecs(1,:) = geo%latVecs(1,:)/real(nCells(1),dp) 
    latVecs(2,:) = geo%latVecs(2,:)/real(nCells(2),dp) 
    latVecs(3,:) = geo%latVecs(3,:)/real(nCells(3),dp) 

    call invert33(invLatt, latVecs) 
    write(stdOut,*) 'reciprocal lattice vectors:'
    write(stdOut,*) 'b1:',invLatt(:,1)
    write(stdOut,*) 'b2:',invLatt(:,2)
    write(stdOut,*) 'b3:',invLatt(:,3)
    invLatt = transpose(invLatt) * 2.0_dp * pi

    allocate(KdynMatrix(3*nAtomUnitCell,3*nAtomUnitCell))
    allocate(eigenValues(3*nAtomUnitCell))

    write(stdOut,*) 'Computing Phonon Dispersion (units '//trim(outputUnits)//')'
    if (tIOProc) then
      open(newunit=fu, file='phononDispersion.dat', action='write')
      if (tWriteTagged) then
        open(newunit=ftag, file=autotestTag) 
      end if
    end if   

    qold(1) = dot_product(invLatt(:,1), kPoint(:,1))
    qold(2) = dot_product(invLatt(:,2), kPoint(:,1))
    qold(3) = dot_product(invLatt(:,3), kPoint(:,1))
    ModKpoint=0.0_dp

    do iK  = 1, nKPoints
      KdynMatrix(:,:) = 0.d0
      q(1) = dot_product(invLatt(:,1), kPoint(:,iK))
      q(2) = dot_product(invLatt(:,2), kPoint(:,iK))
      q(3) = dot_product(invLatt(:,3), kPoint(:,iK))

      write(stdOut,*) ' q:',q(:)
      do  iAtom = 1,  nAtomUnitCell
        do  jAtom = 1, nAtomUnitCell
          ! This loops over all periodic copies
          do  kAtom  = jAtom,  geo%nAtom,  nAtomUnitCell
            DeltaR(:) = geo%Coords(:,kAtom)-geo%Coords(:,jAtom)
            i2 = 3*(iAtom-1)
            j2 = 3*(jAtom-1)
            k2 = 3*(kAtom-1)
            KdynMatrix(i2+1:i2+3,j2+1:j2+3) = KdynMatrix(i2+1:i2+3,j2+1:j2+3) &
               & + dynMatrix(i2+1:i2+3,k2+1:k2+3)*exp(j*dot_product(q,DeltaR))
          end do
        end do
      end do
      
      ! solve the eigenproblem
      call heev(KdynMatrix,eigenValues,'U','N')
 
      ! take square root of modes (allowing for imaginary modes) and print
      eigenValues =  sign(sqrt(abs(eigenValues)),eigenValues)
 
      if (tIOProc) then
        ModKPoint = ModKPoint + sqrt(dot_product(q-qold,q-qold)) 
        qold = q 
        do ii = 1, 3*nAtomUnitCell
          write(fu,*) ModKPoint,  eigenValues(ii)*unitsConv
        end do

        if (tWriteTagged) then
          call tWriter%write(ftag, "kpoint", kPoint(:,iK))
          call tWriter%write(ftag, "bands", eigenValues)
        end if
      end if

    end do


  end subroutine PhononDispersion

  subroutine setConversionUnits(unitsConv)
    real(dp), intent(out) :: unitsConv

    select case(trim(outputUnits))  
    case("H")
      unitsConv = 1.0_dp    
    case("eV")
      unitsConv = Hartree__eV
    case("meV")    
      unitsConv = Hartree__eV*100.0_dp
    case("cm")
      unitsConv = Hartree__cm
    case("THz")
      unitsConv = Hartree__J/(hbar*2.0_dp*pi)*1e-12_dp
    case default
      unitsConv = 1.0_dp    
    end select     
  end subroutine setConversionUnits

  subroutine writeTaggedOut(tWriter, tunnTot, ldosTot, conductance)      
    type(TTaggedWriter) :: tWriter
    real(dp), dimension(:,:) :: tunnTot
    real(dp), dimension(:,:) :: ldosTot
    real(dp), dimension(:,:) :: conductance

    integer :: fu

    open(newunit=fu,file=autotestTag,form="formatted", status="replace")
   
    if (size(tunnTot) > 0) then 
      call tWriter%write(fu,"transmission",tunnTot)
    end if
    if (size(ldosTot) > 0) then 
      call tWriter%write(fu,"PDOS",ldosTot)
    end if
    if (size(conductance) > 0) then 
      call tWriter%write(fu,"conductance",conductance)
    end if

    close(fu)

  end subroutine writeTaggedOut

end program phonons 
