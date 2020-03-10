!!* Program for calculating system normal modes from a Hessian
program phonons 
# include "allocate.h"
# include "assert.h"
  use mpiglobal
  use InitProgram
  use accuracy, only : dp, lc
  use constants, only : Hartree__cm, Bohr__AA, pi
  use TypeGeometry
  use eigensolver, only : heev
  use message
  use libnegf_int
  use ln_structure
  implicit none

  logical :: tInitialized
  real(dp), allocatable :: tunnTot(:,:), ldosTot(:,:), currTot(:)
  logical :: writeLDOS
  logical :: writeTunn

  call mpiglobal_init()

  !! Allocate resources
  call initProgramVariables()
  write (*, "(/,A,/)") "Starting main program"

  if (tCompModes) then
    call ComputeModes()
  end if

  if (tPhonDispersion) then
    call PhononDispersion()
  end if

  if (tTransport) then
    writeTunn = .true.
    writeLDOS = .true.
    call negf_init(str, transpar, tundos, mympi, tInitialized)
    if (.not. tInitialized) call error("libnegf not initialized")

    call negf_init_str(str, transpar, neighborList%iNeighbor, nNeighbor, img2CentCell)
   
    call negf_init_phph(negf, order)

    call calc_phonon_current(mympi, dynMatrix, tunnTot, ldosTot, currTot, &
                        & writeTunn, writeLDOS)  

  end if

  call destructProgramVariables()
 

contains

  subroutine ComputeModes()

    integer  :: ii, jj, kk, ll, nAtom, iMode, jCount, iAt, iAtMoved
    real(dp), allocatable :: eigenValues(:)
    real(dp), allocatable :: displ(:,:,:)
    character(lc) :: lcTmp, lcTmp2
  
    nAtom = geo%nAtom
 
    ALLOCATE_(eigenValues,(3*nMovedAtom))
 
    ! solve the eigenproblem
    if (tPlotModes) then
      write(*,*) 'Computing vibrational frequencies and modes'
      call heev(dynMatrix,eigenValues,'U','V')
    else
      write(*,*) 'Computing vibrational frequencies'
      call heev(dynMatrix,eigenValues,'U','N')
    end if
 
    ! take square root of modes (allowing for imaginary modes) and print
    eigenValues =  sign(sqrt(abs(eigenValues)),eigenValues)
    write(*,*)'Vibrational modes (cm-1):'
    do ii = 1, 3 * nMovedAtom
      write(*,'(f8.2)')eigenValues(ii)*Hartree__cm
    end do
    write(*,*)
 
    if (tPlotModes) then
      write(*,*)'Plotting eigenmodes:'
      write(*,*)ModesToPlot(:)
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
          open(123, file=trim(lcTmp2), position="rewind", status="replace")
          do kk = 1, nCycles
            do ll = 1, nSteps
              write(123,*)nAtom
              write(123,*)'Eigenmode',iMode,eigenValues(iMode)*Hartree__cm,'cm-1'
              do iAt = 1, nAtom
                write(123,'(A3,T4,3F10.6)') &
                    & geo%specieNames(geo%species(iAt)), &
                    & (geo%coords(:,iAt)&
                    & + cos(2.0_dp * pi * real(ll) / real(nSteps))&
                    & * displ(:,iAt,ii)) * Bohr__AA
              end do
            end do
          end do
          close(123)
        end do
      else
        open(123, file="modes.xyz", position="rewind", status="replace")
        do ii = 1, nModesToPlot
          iMode = ModesToPlot(ii)
          write(123,*)nAtom
          write(123,*)'Eigenmode',iMode,eigenValues(iMode)*Hartree__cm,'cm-1'
          if (tXmakeMol) then ! need to account for its non-standard xyz vector
            ! format:
            do iAt = 1, nAtom
              write(123,'(A3,T4,3F10.6,A,3F10.6)') &
                  & geo%specieNames(geo%species(iAt)), &
                  & geo%coords(:,iAt)* Bohr__AA, ' atom_vector ',&
                  & displ(:,iAt,ii)
            end do
          else ! genuine xyz format
            do iAt = 1, nAtom
              write(123,'(A3,T4,6F10.6)') &
                  & geo%specieNames(geo%species(iAt)), &
                  & geo%coords(:,iAt)* Bohr__AA, &
                  & displ(:,iAt,ii)
            end do
          end if
        end do
        close(123)
      end if
      
    end if

    DEALLOCATE_(eigenValues)

  end subroutine ComputeModes

  subroutine PhononDispersion()

    integer  :: ii, jj, kk, ll, nAtom,  iAtom,  iK, jAtom,  kAtom
    integer  :: i2, j2, k2
    real(dp), allocatable :: eigenValues(:)
    real(dp), allocatable :: DeltaR(:), q(:)
    real(dp)::  ModKPoint,  ModDeltaR
    character(lc) :: lcTmp, lcTmp2
    complex(dp), dimension(:,:), allocatable :: KdynMatrix
    COMPLEX(dp), PARAMETER ::    j = (0.d0,1.d0)

    nAtom = geo%nAtom

    Allocate(KdynMatrix(3*nAtomUnitCell,3*nAtomUnitCell))
    Allocate(DeltaR(3))
    Allocate(q(3))

    ALLOCATE_(eigenValues,(3*nAtomUnitCell))

    write(*,*) 'Computing Phonon Dispersion'
    open(unit=70, file='PhononDispersion.dat', action='write')

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
                                 + dynMatrix(i2+1:i2+3,k2+1:k2+3)*exp(j*dot_product(q,DeltaR))
            end do
        end do
      end do
      

    ModKPoint = dsqrt(q(1)**2.0+q(2)**2.0+q(3)**2.0)

    ! solve the eigenproblem
    call heev(KdynMatrix,eigenValues,'U','N')

    ! take square root of modes (allowing for imaginary modes) and print
    eigenValues =  sign(sqrt(abs(eigenValues)),eigenValues)
    !    write(70,*)'     K-Point:  ', q(1:3)
    !    write(70,*)'Vibrational modes (cm-1):'

    do ii = 1, 3*nAtomUnitCell
      write(70,*) iK*1.0,  eigenValues(ii)*Hartree__cm
    end do
  end do

  end subroutine PhononDispersion

end program phonons 
