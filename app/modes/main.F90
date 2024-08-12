!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> The main routines for MODES.
module modes_main
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_constants, only : Hartree__cm, pi
  use dftbp_common_file, only : TFileDescr, closeFile, openFile
  use dftbp_common_globalenv, only : stdOut
  use dftbp_io_formatout, only : writeXYZFormat
  use dftbp_io_message, only : error
  use dftbp_io_taggedoutput, only : TTaggedWriter, TTaggedWriter_init
  use dftbp_math_eigensolver, only : heev, heevd, heevr
#:if WITH_MAGMA
  use dftbp_math_eigensolver, only : gpu_evd
#:endif
  use modes_initmodes, only : TModesMain
  use modes_inputdata, only : solverTypes
#:if WITH_MPI
  use modes_initmodes, only : setEigvecGaugeBlacs
  use modes_modeprojection, only : projectBlacs
#:else
  use modes_initmodes, only : setEigvecGauge
  use modes_modeprojection, only : project
#:endif
  use dftbp_common_environment, only : TEnvironment
  use dftbp_math_matrixops, only : adjointLowerTriangle
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_math_bisect, only : bisection
#:if WITH_SCALAPACK
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip, mpifx_bcast, mpifx_send
  use dftbp_extlibs_scalapackfx, only : scalafx_psyev, scalafx_psyevd, scalafx_psyevr,&
      & scalafx_indxl2g, RSRC_, MB_, CSRC_, NB_, linecomm
#:endif
  implicit none

  private
  public :: runModes


contains

  !> The main MODES program itself.
  subroutine runModes(this, env)

    !> Global variables
    type(TModesMain), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    type(TTaggedWriter) :: taggedWriter

    !! Eigenvectors of the dynamical matrix (scaled by mass and normalized)
    real(dp), allocatable :: eigenModes(:,:)

    integer :: ii, jj, kk, ll, iMode, iAt, nTrans
    real(dp), allocatable :: transDip(:), degenTransDip(:), transPol(:), degenTransPol(:)
    real(dp), allocatable :: eigenModesFull(:,:), eigenModesScaledFull(:,:)
    real(dp) :: zStar(3,3), dMu(3), zStarDeriv(3,3,3), dQ(3,3)

    character(lc) :: lcTmp, lcTmp2
    character(1) :: eigenSolverMode
    type(TFileDescr) :: fd
    logical :: isAppend

    if (this%tEigenVectors) then
      eigenSolverMode = "V"
    else
      eigenSolverMode = "N"
    end if

  #:if WITH_SCALAPACK
    ! remove translations or rotations if necessary
    call projectBlacs(env, this%denseDesc, this%dynMatrix, this%tRemoveTranslate,&
        & this%tRemoveRotate, this%nDerivs, this%nMovedAtom, this%geo, this%atomicMasses)

    select case(this%iSolver)
    case(solverTypes%qr)
      call scalafx_psyev(this%dynMatrix, this%denseDesc%blacsOrbSqr, this%eigen,&
          & this%eigenModesScaled, this%denseDesc%blacsOrbSqr, uplo="U", jobz=eigenSolverMode)
    case(solverTypes%divideAndConquer)
      call scalafx_psyevd(this%dynMatrix, this%denseDesc%blacsOrbSqr, this%eigen,&
          & this%eigenModesScaled, this%denseDesc%blacsOrbSqr, uplo="U", jobz=eigenSolverMode)
    case(solverTypes%relativelyRobust)
      call scalafx_psyevr(this%dynMatrix, this%denseDesc%blacsOrbSqr, this%eigen,&
          & this%eigenModesScaled, this%denseDesc%blacsOrbSqr, uplo="U", jobz=eigenSolverMode)
    case(solverTypes%magmaEvd)
      call error("Magma-solver is not supported together with MPI parallelization.")
    case default
      call error("Unknown eigensolver choice")
    end select

    deallocate(this%dynMatrix)
    call setEigvecGaugeBlacs(env, this%denseDesc, this%eigenModesScaled)
    eigenModes = this%eigenModesScaled
  #:else
    ! remove translations or rotations if necessary
    call project(this%dynMatrix, this%tRemoveTranslate, this%tRemoveRotate, this%nDerivs,&
        & this%nMovedAtom, this%geo, this%atomicMasses)

    select case(this%iSolver)
    case(solverTypes%qr)
      call heev(this%dynMatrix, this%eigen, "U", eigenSolverMode)
    case(solverTypes%divideAndConquer)
      call heevd(this%dynMatrix, this%eigen, "U", eigenSolverMode)
    case(solverTypes%relativelyRobust)
      call heevr(this%dynMatrix, this%eigen, "U", eigenSolverMode)
    case(solverTypes%magmaEvd)
    #:if WITH_MAGMA
      call gpu_evd(env%gpu%nGpu, this%dynMatrix, this%eigen, "U", eigenSolverMode)
    #:else
      call error("Magma-solver selected, but program was compiled without MAGMA")
    #:endif
    case default
      call error("Unknown eigensolver choice")
    end select

    if (this%tEigenVectors) then
      call setEigvecGauge(this%dynMatrix)
    end if
    call move_alloc(this%dynMatrix, eigenModes)

    ! save original eigenvectors
    if (allocated(this%eigenModesScaled)) this%eigenModesScaled(:,:) = eigenModes
  #:endif

    ! take square root of eigenvalues of modes (allowing for imaginary modes)
    this%eigen(:) = sign(sqrt(abs(this%eigen)), this%eigen)

    call TTaggedWriter_init(taggedWriter)
    call openFile(fd, "vibrations.tag", mode="w")

  #:if WITH_SCALAPACK
    call scaleNormalizeEigenmodesBlacs(env, this%denseDesc, this%atomicMasses, eigenModes)
    call displFromEigenmodesBlacs(env, this%denseDesc, this%iMovedAtoms, eigenModes, this%displ)
  #:else
    call scaleNormalizeEigenmodes(this%atomicMasses, eigenModes)
    call displFromEigenmodes(this%iMovedAtoms, eigenModes, this%displ)
  #:endif

    if (allocated(this%bornMatrix)) then
      allocate(transDip(this%nDerivs), source=0.0_dp)
      do jj = 1, this%nDerivs
        dMu(:) = 0.0_dp
        do ii = 1, this%nMovedAtom
          iAt = this%iMovedAtoms(ii)
          zStar(:,:) = reshape(this%bornMatrix(9 * (ii - 1) + 1:9 * ii), [3, 3])
          dMu(:) = dMu + matmul(zStar, this%displ(:, iAt, jj))
        end do
        if (this%eigen(jj) > epsilon(0.0_dp)) then
          transDip(jj) = transDip(jj) + sum(dMu**2)
        end if
      end do
      allocate(degenTransDip(this%nDerivs), source=0.0_dp)
      degenTransDip(1) = transDip(1)
      nTrans = 1
      do jj = 2, this%nDerivs
        ! test for energy degeneracy greater than printing cutoff:
        if (abs(this%eigen(jj) - this%eigen(jj - 1)) * Hartree__cm >= 1.0E-2_dp) then
          nTrans = nTrans + 1
        end if
        degenTransDip(nTrans) = degenTransDip(nTrans) + transDip(jj)
      end do
    end if

    if (allocated(this%bornDerivsMatrix)) then
      allocate(transPol(this%nDerivs), source=0.0_dp)
      do jj = 1, this%nDerivs
        dQ(:,:) = 0.0_dp
        do ii = 1, this%nMovedAtom
          iAt = this%iMovedAtoms(ii)
          zStarDeriv(:,:,:) = reshape(this%bornDerivsMatrix(27 * (ii - 1) + 1:27 * ii), [3, 3, 3])
          dQ(:,:) = dQ + reshape(matmul(reshape(zStarDeriv, [9, 3]),&
              & this%displ(:, iAt, jj)), [3, 3])
        end do
        if (this%eigen(jj) > epsilon(0.0_dp)) then
          transPol(jj) = transPol(jj) + sum(dQ**2)
        end if
      end do
      allocate(degenTransPol(this%nDerivs), source=0.0_dp)
      degenTransPol(1) = transPol(1)
      nTrans = 1
      do jj = 2, this%nDerivs
        ! test for energy degeneracy greater than printing cutoff:
        if (abs(this%eigen(jj) - this%eigen(jj - 1)) * Hartree__cm >= 1.0E-2_dp) then
          nTrans = nTrans + 1
        end if
        degenTransPol(nTrans) = degenTransPol(nTrans) + transPol(jj)
      end do
    end if

    if (this%tPlotModes) then
      call taggedWriter%write(fd%unit, "saved_modes", this%modesToPlot)
      write(stdout, *) "Writing eigenmodes to vibrations.tag"
    #:if WITH_SCALAPACK
      call getFullFromDistributed(env, this%denseDesc, eigenModes, eigenModesFull)
      call getFullFromDistributed(env, this%denseDesc, this%eigenModesScaled, eigenModesScaledFull)
    #:else
      eigenModesFull = eigenModes
      eigenModesScaledFull = this%eigenModesScaled
    #:endif
      if (env%tGlobalLead) then
        call taggedWriter%write(fd%unit, "eigenmodes", eigenModesFull(:, this%modesToPlot))
      end if
      write(stdout, *) "Plotting eigenmodes:"
      write(stdout, "(16I5)") this%modesToPlot
      if (env%tGlobalLead) then
        call taggedWriter%write(fd%unit, "eigenmodes_scaled",&
            & eigenModesScaledFull(:, this%modesToPlot))
      end if
      if (this%tAnimateModes) then
        do ii = 1, this%nModesToPlot
          iMode = this%modesToPlot(ii)
          write(lcTmp,"('mode_',I0,'.xyz')") iMode
          do kk = 1, this%nCycles
            do ll = 1, this%nSteps
              isAppend = (kk > 1 .or. ll > 1)
              write(lcTmp2, *) "Eigenmode", iMode, this%eigen(iMode) * Hartree__cm, "cm-1"
              call writeXYZFormat(lcTmp, this%geo%coords + cos(2.0_dp * pi * real(ll)&
                  & / real(this%nSteps)) * this%displ(:,:, iMode), this%geo%species,&
                  & this%geo%speciesNames, comment=trim(lcTmp2), append=isAppend)
            end do
          end do
        end do
      else
        lcTmp = "modes.xyz"
        do ii = 1, this%nModesToPlot
          isAppend = (ii > 1)
          iMode = this%modesToPlot(ii)
          write(lcTmp2, *) "Eigenmode", iMode, this%eigen(iMode) * Hartree__cm, "cm-1"
          call writeXYZFormat(lcTmp, this%geo%coords, this%geo%species, this%geo%speciesNames,&
              & vectors=this%displ(:,:, iMode), comment=trim(lcTmp2), append=isAppend)
        end do
      end if
    end if

    write(stdout, *) "Vibrational modes"
    if (allocated(this%bornMatrix) .and. allocated(this%bornDerivsMatrix)) then
      write(stdout, "(T7,A,T16,A,T28,A)") "freq.", "IR", "Polarisability"
      write(stdout, "(A,T7,A,T16,A,T28,A)") "Mode", "/ cm-1", "/ a.u.", "change / a.u."
    else if (allocated(this%bornMatrix)) then
      write(stdout, "(T7,A,T16,A)") "freq.", "IR"
      write(stdout, "(A,T7,A,T16,A)") "Mode", "/ cm-1", "/ a.u."
    else if (allocated(this%bornDerivsMatrix)) then
      write(stdout, "(T7,A,T16,A)") "freq.", "Polarisability"
      write(stdout, "(A,T7,A,T16,A)") "Mode", "/ cm-1", "change / a.u."
    else
      write(stdout, "(T7,A)") "freq."
      write(stdout, "(A,T7,A)") "Mode", "cm-1"
    end if
    if (allocated(this%bornMatrix) .and. allocated(this%bornDerivsMatrix)) then
      do ii = 1, this%nDerivs
        write(stdout, "(i5,f8.2,2E12.4)") ii, this%eigen(ii) * Hartree__cm, transDip(ii),&
            & transPol(ii)
      end do
    else if (allocated(this%bornMatrix)) then
      do ii = 1, this%nDerivs
        write(stdout, "(i5,f8.2,E12.4)") ii, this%eigen(ii) * Hartree__cm, transDip(ii)
      end do
    else if (allocated(this%bornDerivsMatrix)) then
      do ii = 1, this%nDerivs
        write(stdout, "(i5,f8.2,E12.4)") ii, this%eigen(ii) * Hartree__cm, transPol(ii)
      end do
    else
      do ii = 1, this%nDerivs
        write(stdout, "(i5,f8.2)") ii, this%eigen(ii) * Hartree__cm
      end do
    end if
    write(stdout, *)

    if (env%tGlobalLead) then
      call taggedWriter%write(fd%unit, "frequencies", this%eigen)
    end if

    if (allocated(this%bornMatrix)) then
      if (env%tGlobalLead) then
        call taggedWriter%write(fd%unit, "intensities", degenTransDip(:nTrans))
      end if
    end if

    if (allocated(this%bornDerivsMatrix)) then
      if (env%tGlobalLead) then
        if (this%tRemoveTranslate .or. this%tRemoveRotate) then
          call taggedWriter%write(fd%unit, "scattering", degenTransPol(2:nTrans))
        else
          call taggedWriter%write(fd%unit, "scattering", degenTransPol(:nTrans))
        end if
      end if
    end if

    call closeFile(fd)

  end subroutine runModes


#:if WITH_SCALAPACK
  !> Collects full, collected square matrix from distributed contributions.
  subroutine getFullFromDistributed(env, denseDesc, distrib, collected)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Distributed matrix (source)
    real(dp), intent(in) :: distrib(:,:)

    !> Collected, full square matrix (target)
    real(dp), intent(out), allocatable :: collected(:,:)

    !! Type for communicating a row or a column of a distributed matrix
    type(linecomm) :: collector

    !! Temporary storage for a single line of the collected, dense, square matrix
    real(dp), allocatable :: localLine(:)

    !! Auxiliary variables
    integer :: iDerivs, nDerivs

    nDerivs = denseDesc%fullSize
    allocate(localLine(nDerivs))
    if (env%tGlobalLead) allocate(collected(nDerivs, nDerivs), source=0.0_dp)

    call collector%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")

    if (env%tGlobalLead) then
      do iDerivs = 1, nDerivs
        call collector%getline_lead(env%blacs%orbitalGrid, iDerivs, distrib, localLine)
        collected(:, iDerivs) = localLine
      end do
    else
      do iDerivs = 1, nDerivs
        if (env%mpi%tGroupLead) then
          call collector%getline_lead(env%blacs%orbitalGrid, iDerivs, distrib, localLine)
          call mpifx_send(env%mpi%interGroupComm, localLine, env%mpi%interGroupComm%leadrank)
        else
          call collector%getline_follow(env%blacs%orbitalGrid, iDerivs, distrib)
        end if
      end do
    end if

  end subroutine getFullFromDistributed


  !> Scale mode components on each atom by mass and then normalizes the total mode.
  subroutine scaleNormalizeEigenmodesBlacs(env, denseDesc, atomicMasses, eigenModes)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Atomic masses to build dynamical matrix
    real(dp), intent(in) :: atomicMasses(:)

    !> Eigenvectors of the dynamical matrix (scaled by mass and normalized)
    real(dp), intent(inout) :: eigenModes(:,:)

    !! Number of derivatives
    integer :: nDerivs

    !! Auxiliary variables
    integer :: iCount, jCount, iDeriv, iAt
    real(dp), allocatable :: sqrtModes(:)

    ! scale mode components by mass of atoms
    do iCount = 1, size(eigenModes, dim=2)
      do jCount = 1, size(eigenModes, dim=1)
        iDeriv = scalafx_indxl2g(jCount, denseDesc%blacsOrbSqr(MB_),&
            & env%blacs%orbitalGrid%myrow, denseDesc%blacsOrbSqr(RSRC_),&
            & env%blacs%orbitalGrid%nrow)
        call bisection(iAt, denseDesc%iAtomStart, iDeriv)
        eigenModes(jCount, iCount) = eigenModes(jCount, iCount) / sqrt(atomicMasses(iAt))
      end do
    end do

    nDerivs = denseDesc%fullSize
    allocate(sqrtModes(nDerivs), source=0.0_dp)

    ! collect square roots of distributed eigenvectors
    do iCount = 1, size(eigenModes, dim=2)
      iDeriv = scalafx_indxl2g(iCount, denseDesc%blacsOrbSqr(NB_), env%blacs%orbitalGrid%mycol,&
          & denseDesc%blacsOrbSqr(CSRC_), env%blacs%orbitalGrid%ncol)
      sqrtModes(iDeriv) = sum(eigenModes(:, iCount)**2)
    end do

    call mpifx_allreduceip(env%mpi%globalComm, sqrtModes, MPI_SUM)
    sqrtModes(:) = sqrt(sqrtModes)

    ! normalize total modes
    do iCount = 1, size(eigenModes, dim=2)
      iDeriv = scalafx_indxl2g(iCount, denseDesc%blacsOrbSqr(NB_), env%blacs%orbitalGrid%mycol,&
          & denseDesc%blacsOrbSqr(CSRC_), env%blacs%orbitalGrid%ncol)
      eigenModes(:, iCount) = eigenModes(:, iCount) / sqrtModes(iDeriv)
    end do

  end subroutine scaleNormalizeEigenmodesBlacs


  !> Creates displacement vectors for every atom in every mode.
  subroutine displFromEigenmodesBlacs(env, denseDesc, iMovedAtoms, eigenModes, displ)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> List of atoms in dynamical matrix
    integer, intent(in) :: iMovedAtoms(:)

    !> Eigenvectors of the dynamical matrix (scaled by mass and normalized)
    real(dp), intent(in) :: eigenModes(:,:)

    !> Displacement vectors for every atom in every mode. Shape: [3, nAtom, nDerivs]
    real(dp), intent(out) :: displ(:,:,:)

    !! Type for communicating a row or a column of a distributed matrix
    type(linecomm) :: collector

    !! Temporary storage for a single line of the collected, dense, square matrix
    real(dp), allocatable :: localLine(:)

    !! Number of derivatives
    integer :: nDerivs

    !! Total number of atoms
    integer :: nAtom

    !! Auxiliary variables
    integer :: iAt, iAtMoved, iDerivs

    nAtom = size(displ, dim=2)
    nDerivs = size(displ, dim=3)

    displ(:,:,:) = 0.0_dp

    allocate(localLine(nDerivs))

    call collector%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")

    if (env%tGlobalLead) then
      ! runs over eigenvectors
      do iDerivs = 1, nDerivs
        call collector%getline_lead(env%blacs%orbitalGrid, iDerivs, eigenModes, localLine)
        do iAt = 1, nAtom
          if (any(iMovedAtoms == iAt)) then
            ! index of atom in the list of moved atoms
            iAtMoved = minloc(abs(iMovedAtoms - iAt), dim=1)
            displ(:, iAt, iDerivs) = localLine(3 * iAtMoved - 2:3 * iAtMoved)
          end if
        end do
      end do
    else
      do iDerivs = 1, nDerivs
        if (env%mpi%tGroupLead) then
          call collector%getline_lead(env%blacs%orbitalGrid, iDerivs, eigenModes, localLine)
          call mpifx_send(env%mpi%interGroupComm, localLine, env%mpi%interGroupComm%leadrank)
        else
          call collector%getline_follow(env%blacs%orbitalGrid, iDerivs, eigenModes)
        end if
      end do
    end if

    call mpifx_bcast(env%mpi%globalComm, displ)

  end subroutine displFromEigenmodesBlacs

#:else

  !> Scale mode components on each atom by mass and then normalizes the total mode.
  subroutine scaleNormalizeEigenmodes(atomicMasses, eigenModes)

    !> Atomic masses to build dynamical matrix
    real(dp), intent(in) :: atomicMasses(:)

    !> Eigenvectors of the dynamical matrix (scaled by mass and normalized)
    real(dp), intent(inout) :: eigenModes(:,:)

    !! Number of derivatives
    integer :: nDerivs

    !! Number of atoms which should be moved.
    integer :: nMovedAtom

    !! Auxiliary variables
    integer :: ii, jj, ll, jCount

    nDerivs = size(eigenModes, dim=2)
    nMovedAtom = nDerivs / 3

    do ii = 1, nDerivs
      jCount = 0
      do jj = 1, nMovedAtom
        do ll = 1, 3
          jCount = jCount + 1
          eigenModes(jCount, ii) = eigenModes(jCount, ii) / sqrt(atomicMasses(jj))
        end do
      end do
      eigenModes(:, ii) = eigenModes(:, ii) / sqrt(sum(eigenModes(:, ii)**2))
    end do

  end subroutine scaleNormalizeEigenmodes


  !> Creates displacement vectors for every atom in every mode.
  subroutine displFromEigenmodes(iMovedAtoms, eigenModes, displ)

    !> List of atoms in dynamical matrix
    integer, intent(in) :: iMovedAtoms(:)

    !> Eigenvectors of the dynamical matrix (scaled by mass and normalized)
    real(dp), intent(in) :: eigenModes(:,:)

    !> Displacement vectors for every atom in every mode. Shape: [3, nAtom, nDerivs]
    real(dp), intent(out) :: displ(:,:,:)

    !! Number of derivatives
    integer :: nDerivs

    !! Total number of atoms
    integer :: nAtom

    !! Auxiliary variables
    integer :: iAt, iAtMoved, ii

    nAtom = size(displ, dim=2)
    nDerivs = size(displ, dim=3)

    displ(:,:,:) = 0.0_dp

    do iAt = 1, nAtom
      if (any(iMovedAtoms == iAt)) then
        ! index of atom in the list of moved atoms
        iAtMoved = minloc(abs(iMovedAtoms - iAt), 1)
        do ii = 1, nDerivs
          displ(:, iAt, ii) = eigenModes(3 * iAtMoved - 2:3 * iAtMoved, ii)
        end do
      end if
    end do

  end subroutine displFromEigenmodes
#:endif

end module modes_main
