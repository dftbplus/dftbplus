!!* Program for calculating phonon transport from an Hessian

#:include 'common.fypp'

program phonons 
  use dftbp_assert
  use dftbp_globalenv
  use dftbp_environment
  use dftbp_initphonons
  use dftbp_accuracy, only : dp, lc
  use dftbp_constants, only : Hartree__cm, Bohr__AA, pi
  use dftbp_typegeometry
  use dftbp_eigensolver, only : heev
  use dftbp_message
  use libnegf_int
  use ln_structure
  implicit none

  type(TEnvironment) :: env
  logical :: tInitialized
  real(dp), allocatable :: tunnTot(:,:), ldosTot(:,:), currLead(:)
  logical :: twriteLDOS
  logical :: twriteTunn

  call initGlobalEnv()
  call printHeader()
  call TEnvironment_init(env)
  call initProgramVariables(env)

  if (tCompModes) then
    call ComputeModes()
  end if

  if (tPhonDispersion) then
    call PhononDispersion()
  end if

  if (tTransport) then
    twriteTunn = .true.
    twriteLDOS = .true.
    call negf_init(env, transpar, tundos, tInitialized)

    if (.not. tInitialized) then
      call error("libnegf not initialized")
    end if

    call negf_init_str(geo%nAtom, transpar, neighbourList%iNeighbour, nNeighbour, img2CentCell)
   
    call calc_phonon_current(env, dynMatrix, tunnTot, ldosTot, currLead, &
                        & twriteTunn, twriteLDOS)  
  end if

  call destructProgramVariables()

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
    write (stdOut, "(A)") "L. Medrano Sandonas et al. Quantum Phonon transport in Nanomaterials: ..., Entropy 735 (2019)" 
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

  subroutine PhononDispersion()

    integer  :: ii, jj, kk, ll, nAtom,  iAtom,  iK, jAtom,  kAtom
    integer  :: i2, j2, k2, fu
    real(dp), allocatable :: eigenValues(:)
    real(dp), allocatable :: DeltaR(:), q(:)
    real(dp)::  ModKPoint,  ModDeltaR
    character(lc) :: lcTmp, lcTmp2
    complex(dp), dimension(:,:), allocatable :: KdynMatrix
    complex(dp), parameter ::    j = (0.d0,1.d0)

    nAtom = geo%nAtom

    allocate(KdynMatrix(3*nAtomUnitCell,3*nAtomUnitCell))
    allocate(DeltaR(3))
    allocate(q(3))
    allocate(eigenValues(3*nAtomUnitCell))

    write(stdOut,*) 'Computing Phonon Dispersion'
    open(newunit=fu, file='PhononDispersion.dat', action='write')

    do iK  = 1, nKPoints
      KdynMatrix(:,:) = 0.d0
      q(:)  = KPoint(iK,:)
      do  iAtom = 1,  nAtomUnitCell
        do  jAtom = 1, nAtomUnitCell
          do  kAtom  = jAtom,  nAtom,  nAtomUnitCell
            DeltaR(:) = geo%Coords(:,kAtom)-geo%Coords(:,jAtom)
            i2 = 3*(iAtom-1)
            j2 = 3*(jAtom-1)
            k2 = 3*(kAtom-1)
            KdynMatrix(i2+1:i2+3,j2+1:j2+3) = KdynMatrix(i2+1:i2+3,j2+1:j2+3) &
               & + dynMatrix(i2+1:i2+3,k2+1:k2+3)*exp(j*dot_product(q,DeltaR))
          end do
        end do
      end do
      
      ModKPoint = dsqrt(q(1)**2.0+q(2)**2.0+q(3)**2.0)
 
      ! solve the eigenproblem
      call heev(KdynMatrix,eigenValues,'U','N')
 
      ! take square root of modes (allowing for imaginary modes) and print
      eigenValues =  sign(sqrt(abs(eigenValues)),eigenValues)
 
      do ii = 1, 3*nAtomUnitCell
        write(fu,*) iK*1.0,  eigenValues(ii)*Hartree__cm
      end do
    end do

  end subroutine PhononDispersion

end program phonons 
