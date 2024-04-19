!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Program for calculating system normal modes from a Hessian.
program modes
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_constants, only : Hartree__cm, pi
  use dftbp_common_file, only : TFileDescr, closeFile, openFile
  use dftbp_common_globalenv, only : stdOut
  use dftbp_io_formatout, only : writeXYZFormat
  use dftbp_io_taggedoutput, only : TTaggedWriter, TTaggedWriter_init
  use dftbp_math_eigensolver, only : heev
  use modes_initmodes, only : dynMatrix, bornMatrix, bornDerivsMatrix, modesToPlot, geo,&
      & iMovedAtoms, nCycles, nDerivs, nModesToPlot, nMovedAtom, nSteps, tAnimateModes, tPlotModes,&
      & tEigenVectors, tRemoveRotate, tRemoveTranslate, atomicMasses, initProgramVariables
  use modes_modeprojection, only : project
#:if WITH_MPI
  use mpi, only : MPI_THREAD_FUNNELED
  use dftbp_common_mpienv, only : TMpiEnv, TMpiEnv_init
  use dftbp_extlibs_mpifx, only : mpifx_init_thread, mpifx_finalize
#:endif
  implicit none

  type(TTaggedWriter) :: taggedWriter

  integer :: ii, jj, kk, ll, iMode, iAt, iAtMoved, nAtom, nTrans
  integer :: iCount, jCount
  real(dp), allocatable :: eigenValues(:), eigenModesScaled(:,:), displ(:,:,:)
  real(dp), allocatable :: transDip(:), degenTransDip(:), transPol(:), degenTransPol(:)
  real(dp) :: zStar(3,3), dMu(3), zStarDeriv(3,3,3), dQ(3,3)

  character(lc) :: lcTmp, lcTmp2
  type(TFileDescr) :: fd
  logical :: isAppend

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
  if (tPlotModes) allocate(eigenModesScaled(3 * nMovedAtom, 3 * nMovedAtom))

  ! mass weight the Hessian matrix to get the dynamical matrix
  ! H_{ij} = \frac{\partial^2 \Phi}{\partial u_i \partial u_j}
  ! D_{ij} = \frac{H_{ij}}{\sqrt{m_i m_j}}
  !        = \frac{\partial^2 \Phi}{\partial w_i \partial w_j}
  ! where w_i = \sqrt{m_i} u_i
  iCount = 0
  do ii = 1, nMovedAtom
    do kk = 1, 3
      iCount = iCount + 1
      jCount = 0
      do jj = 1, nMovedAtom
        do ll = 1, 3
          jCount = jCount + 1
          dynMatrix(jCount, iCount) = dynMatrix(jCount, iCount)&
              & / (sqrt(atomicMasses(ii)) * sqrt(atomicMasses(jj)))
        end do
      end do
    end do
  end do

  ! remove translations or rotations if necessary
  call project(dynMatrix, tRemoveTranslate, tRemoveRotate, nDerivs, nMovedAtom, geo, atomicMasses)

  ! solve the eigenproblem
  if (tEigenVectors) then
    call heev(dynMatrix, eigenValues, "U", "V")
  else
    call heev(dynMatrix, eigenValues, "U", "N")
  end if

  ! save original eigenvectors
  if (allocated(eigenModesScaled)) eigenModesScaled(:,:) = dynMatrix

  ! take square root of eigenvalues of modes (allowing for imaginary modes)
  eigenValues(:) = sign(sqrt(abs(eigenValues)), eigenValues)

  call TTaggedWriter_init(taggedWriter)
  call openFile(fd, "vibrations.tag", mode="w")

  ! scale mode components on each atom by mass and then normalise total mode
  do ii = 1, nDerivs
    jCount = 0
    do jj = 1, nMovedAtom
      do ll = 1, 3
        jCount = jCount + 1
        dynMatrix(jCount, ii) = dynMatrix(jCount, ii) / sqrt(atomicMasses(jj))
      end do
    end do
    dynMatrix(:, ii) = dynMatrix(:, ii) / sqrt(sum(dynMatrix(:, ii)**2))
  end do

  nAtom = geo%nAtom
  allocate(displ(3, nAtom, nDerivs))
  displ(:,:,:) = 0.0_dp
  ! Create displacement vectors for every atom in every mode.
  do iAt = 1, nAtom
    if (any(iMovedAtoms == iAt)) then
      ! Index of atom in the list of moved atoms
      iAtMoved = minloc(abs(iMovedAtoms - iAt), 1)
      do ii = 1, nDerivs
        displ(:, iAt, ii) = dynMatrix(3 * iAtMoved - 2:3 * iAtMoved, ii)
      end do
    end if
  end do

  if (allocated(bornMatrix)) then
    allocate(transDip(nDerivs), source=0.0_dp)
    do jj = 1, nDerivs
      dMu(:) = 0.0_dp
      do ii = 1, nMovedAtom
        iAt = iMovedAtoms(ii)
        zStar(:,:) = reshape(bornMatrix(9 * (ii - 1) + 1:9 * ii), [3, 3])
        dMu(:) = dMu + matmul(zStar, displ(:, iAt, jj))
      end do
      if (eigenValues(jj) > epsilon(0.0_dp)) then
        transDip(jj) = transDip(jj) + sum(dMu**2)
      end if
    end do
    allocate(degenTransDip(nDerivs), source=0.0_dp)
    degenTransDip(1) = transDip(1)
    nTrans = 1
    do jj = 2, nDerivs
      ! test for energy degeneracy greater than printing cutoff:
      if (abs(eigenValues(jj) - eigenValues(jj - 1)) * Hartree__cm >= 1.0E-2_dp) then
        nTrans = nTrans + 1
      end if
      degenTransDip(nTrans) = degenTransDip(nTrans) + transDip(jj)
    end do
  end if

  if (allocated(bornDerivsMatrix)) then
    allocate(transPol(nDerivs), source=0.0_dp)
    do jj = 1, nDerivs
      dQ(:,:) = 0.0_dp
      do ii = 1, nMovedAtom
        iAt = iMovedAtoms(ii)
        zStarDeriv(:,:,:) = reshape(bornDerivsMatrix(27 * (ii - 1) + 1:27 * ii), [3, 3, 3])
        dQ(:,:) = dQ + reshape(matmul(reshape(zStarDeriv, [9, 3]),  displ(:, iAt, jj)), [3, 3])
      end do
      if (eigenValues(jj) > epsilon(0.0_dp)) then
        transPol(jj) = transPol(jj) + sum(dQ**2)
      end if
    end do
    allocate(degenTransPol(nDerivs), source=0.0_dp)
    degenTransPol(1) = transPol(1)
    nTrans = 1
    do jj = 2, nDerivs
      ! test for energy degeneracy greater than printing cutoff:
      if (abs(eigenValues(jj) - eigenValues(jj - 1)) * Hartree__cm >= 1.0E-2_dp) then
        nTrans = nTrans + 1
      end if
      degenTransPol(nTrans) = degenTransPol(nTrans) + transPol(jj)
    end do
  end if

  if (tPlotModes) then
    call taggedWriter%write(fd%unit, "saved_modes", modesToPlot)
    write(stdout, *) "Writing eigenmodes to vibrations.tag"
    call taggedWriter%write(fd%unit, "eigenmodes", dynMatrix(:, modesToPlot))
    write(stdout, *) "Plotting eigenmodes:"
    write(stdout, "(16I5)") modesToPlot(:)
    call taggedWriter%write(fd%unit, "eigenmodes_scaled", eigenModesScaled(:, modesToPlot))
    if (tAnimateModes) then
      do ii = 1, nModesToPlot
        iMode = modesToPlot(ii)
        write(lcTmp,"('mode_',I0,'.xyz')") iMode
        do kk = 1, nCycles
          do ll = 1, nSteps
            isAppend = (kk > 1 .or. ll > 1)
            write(lcTmp2, *) "Eigenmode", iMode, eigenValues(iMode) * Hartree__cm, "cm-1"
            call writeXYZFormat(lcTmp,&
                & geo%coords + cos(2.0_dp * pi * real(ll) / real(nSteps)) * displ(:,:, iMode),&
                & geo%species, geo%speciesNames, comment=trim(lcTmp2), append=isAppend)
          end do
        end do
      end do
    else
      lcTmp = "modes.xyz"
      do ii = 1, nModesToPlot
        isAppend = (ii > 1)
        iMode = modesToPlot(ii)
        write(lcTmp2, *) "Eigenmode", iMode, eigenValues(iMode) * Hartree__cm, "cm-1"
        call writeXYZFormat(lcTmp, geo%coords, geo%species, geo%speciesNames,&
            & vectors=displ(:,:, iMode), comment=trim(lcTmp2), append=isAppend)
      end do
    end if
  end if

  write(stdout, *) "Vibrational modes"
  if (allocated(bornMatrix) .and. allocated(bornDerivsMatrix)) then
    write(stdout, "(T7,A,T16,A,T28,A)") "freq.", "IR", "Polarisability"
    write(stdout, "(A,T7,A,T16,A,T28,A)") "Mode", "/ cm-1", "/ a.u.", "change / a.u."
  else if (allocated(bornMatrix)) then
    write(stdout, "(T7,A,T16,A)") "freq.", "IR"
    write(stdout, "(A,T7,A,T16,A)") "Mode", "/ cm-1", "/ a.u."
  else if (allocated(bornDerivsMatrix)) then
    write(stdout, "(T7,A,T16,A)") "freq.", "Polarisability"
    write(stdout, "(A,T7,A,T16,A)") "Mode", "/ cm-1", "change / a.u."
  else
    write(stdout, "(T7,A)") "freq."
    write(stdout, "(A,T7,A)") "Mode", "cm-1"
  end if
  if (allocated(bornMatrix) .and. allocated(bornDerivsMatrix)) then
    do ii = 1, 3 * nMovedAtom
      write(stdout, "(i5,f8.2,2E12.4)") ii, eigenValues(ii) * Hartree__cm, transDip(ii),&
          & transPol(ii)
    end do
  else if (allocated(bornMatrix)) then
    do ii = 1, 3 * nMovedAtom
      write(stdout, "(i5,f8.2,E12.4)") ii, eigenValues(ii) * Hartree__cm, transDip(ii)
    end do
  else if (allocated(bornDerivsMatrix)) then
    do ii = 1, 3 * nMovedAtom
      write(stdout, "(i5,f8.2,E12.4)") ii, eigenValues(ii) * Hartree__cm, transPol(ii)
    end do
  else
    do ii = 1, 3 * nMovedAtom
      write(stdout, "(i5,f8.2)") ii,eigenValues(ii) * Hartree__cm
    end do
  end if
  write(stdout, *)

  call taggedWriter%write(fd%unit, "frequencies", eigenValues)

  if (allocated(bornMatrix)) then
    call taggedWriter%write(fd%unit, "intensities", degenTransDip(:nTrans))
  end if

  if (allocated(bornDerivsMatrix)) then
    if (tRemoveTranslate .or. tRemoveRotate) then
      call taggedWriter%write(fd%unit, "scattering", degenTransPol(2:nTrans))
    else
      call taggedWriter%write(fd%unit, "scattering", degenTransPol(:nTrans))
    end if
  end if

  call closeFile(fd)

#:if WITH_MPI
  call mpifx_finalize()
#:endif

end program modes
