!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Program for calculating system normal modes from a Hessian
program modes
  use dftbp_assert
  use dftbp_globalenv, only : stdOut
  use dftbp_initmodes
  use dftbp_accuracy, only : dp, lc
  use dftbp_constants, only : Hartree__cm, Bohr__AA, pi
  use dftbp_typegeometry
  use dftbp_eigensolver, only : heev
  use dftbp_taggedoutput
  use dftbp_message
  use dftbp_modeprojection
#:if WITH_MPI
  use dftbp_mpienv, only : TMpiEnv, TMpiEnv_init
  use dftbp_mpifx, only : mpifx_init_thread, mpifx_finalize
  use mpi, only : MPI_THREAD_FUNNELED
#:endif
  implicit none

  type(TTaggedWriter) :: taggedWriter

  integer :: ii, jj, kk, ll, iMode, iAt, iAtMoved, nAtom
  integer :: iCount, jCount
  real(dp), allocatable :: eigenValues(:)
  real(dp), allocatable :: displ(:,:,:)

  character(lc) :: lcTmp, lcTmp2
  integer :: fdUnit

#:if WITH_MPI
  !> MPI environment, if compiled with mpifort
  type(TMpiEnv) :: mpiEnv

  ! As this is serial code, trap for run time execution on more than 1 processor with an mpi enabled
  ! build
  call mpifx_init_thread(requiredThreading=MPI_THREAD_FUNNELED)
  call TMpiEnv_init(mpiEnv)
  call mpiEnv%mpiSerialEnv()
#:endif

  ! Allocate resources
  call initProgramVariables()
  write(stdout, "(/,A,/)") "Starting main program"

  allocate(eigenValues(3 * nMovedAtom))

  ! mass weight the Hessian matrix to get the dynamical matrix
  iCount = 0
  do ii = 1, nMovedAtom
    do kk = 1, 3
      iCount = iCount + 1
      jCount = 0
      do jj = 1, nMovedAtom
        do ll = 1, 3
          jCount = jCount + 1
          dynMatrix(jCount,iCount) = dynMatrix(jCount,iCount) &
              & / (sqrt(atomicMasses(ii)) * sqrt(atomicMasses(jj)))
        end do
      end do
    end do
  end do

  ! remove translations or rotations if neccessary
  call project(dynMatrix, tRemoveTranslate, tRemoveRotate, nDerivs, nMovedAtom, geo, atomicMasses)

  ! solve the eigenproblem
  if (tPlotModes) then
    call heev(dynMatrix,eigenValues,'U','V')
  else
    call heev(dynMatrix,eigenValues,'U','N')
  end if

  ! take square root of modes (allowing for imaginary modes) and print
  eigenValues =  sign(sqrt(abs(eigenValues)),eigenValues)
  write(stdout, *)'Vibrational modes (cm-1):'
  do ii = 1, 3 * nMovedAtom
    write(stdout, '(i5,f8.2)')ii,eigenValues(ii)*Hartree__cm
  end do
  write(stdout, *)

  call TTaggedWriter_init(taggedWriter)
  open(newunit=fdUnit, file="vibrations.tag", form="formatted", status="replace")
  call taggedWriter%write(fdUnit, "frequencies", eigenValues)

  if (tPlotModes) then
    call taggedWriter%write(fdUnit, "saved_modes", modesToPlot)
    write(stdout, *) "Writing eigenmodes to vibrations.tag"
    call taggedWriter%write(fdUnit, "eigenmodes", dynMatrix(:,ModesToPlot))

    write(stdout, *)'Plotting eigenmodes:'
    write(stdout, *)ModesToPlot(:)
    ! scale mode components on each atom by mass and then normalise total mode
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
    call taggedWriter%write(fdUnit, "eigenmodes_scaled", dynMatrix(:,ModesToPlot))
    close(fdUnit)

    ! Create displacment vectors for every atom in every mode.
    nAtom = geo%nAtom
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
        open(newunit=fdUnit, file=trim(lcTmp2), position="rewind", status="replace")
        do kk = 1, nCycles
          do ll = 1, nSteps
            write(fdUnit,*)nAtom
            write(fdUnit,*)'Eigenmode',iMode,eigenValues(iMode)*Hartree__cm,'cm-1'
            do iAt = 1, nAtom
              write(fdUnit,'(A3,T4,3F10.6)') &
                  & geo%speciesNames(geo%species(iAt)), &
                  & (geo%coords(:,iAt)&
                  & + cos(2.0_dp * pi * real(ll) / real(nSteps))&
                  & * displ(:,iAt,ii)) * Bohr__AA
            end do
          end do
        end do
        close(fdUnit)
      end do
    else
      open(newunit=fdUnit, file="modes.xyz", position="rewind", status="replace")
      do ii = 1, nModesToPlot
        iMode = ModesToPlot(ii)
        write(fdUnit,*)nAtom
        write(fdUnit,*)'Eigenmode',iMode,eigenValues(iMode)*Hartree__cm,'cm-1'
        if (tXmakeMol) then
          ! need to account for its non-standard xyz vector format:
          do iAt = 1, nAtom
            write(fdUnit,'(A3,T4,3F10.6,A,3F10.6)') &
                & geo%speciesNames(geo%species(iAt)), &
                & geo%coords(:,iAt)* Bohr__AA, ' atom_vector ',&
                & displ(:,iAt,ii)
          end do
        else
          ! genuine xyz format
          do iAt = 1, nAtom
            write(fdUnit,'(A3,T4,6F10.6)') &
                & geo%speciesNames(geo%species(iAt)), &
                & geo%coords(:,iAt)* Bohr__AA, &
                & displ(:,iAt,ii)
          end do
        end if
      end do
      close(fdUnit)
    end if

  end if

#:if WITH_MPI
  call mpifx_finalize()
#:endif

end program modes
