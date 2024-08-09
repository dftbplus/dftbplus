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
  use dftbp_io_taggedoutput, only : TTaggedWriter, TTaggedWriter_init
  use dftbp_math_eigensolver, only : heev
  use modes_initmodes, only : TModesMain
  use dftbp_common_environment, only : TEnvironment
  use modes_modeprojection, only : project
#:if WITH_SCALAPACK
  use dftbp_extlibs_scalapackfx, only : scalafx_psyev
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

    integer :: ii, jj, kk, ll, iMode, iAt, iAtMoved, nTrans
    integer :: jCount
    real(dp), allocatable :: transDip(:), degenTransDip(:), transPol(:), degenTransPol(:)
    real(dp) :: zStar(3,3), dMu(3), zStarDeriv(3,3,3), dQ(3,3)

    character(lc) :: lcTmp, lcTmp2
    type(TFileDescr) :: fd
    logical :: isAppend

    #:if WITH_SCALAPACK
    ! remove translations or rotations if necessary
    ! call projectBlacs()
    if (this%tEigenVectors) then
      call scalafx_psyev(this%dynMatrix, this%denseDesc%blacsOrbSqr, this%eigen,&
          & this%eigenModesScaled, this%denseDesc%blacsOrbSqr, uplo="L", jobz="V")
    else
      call scalafx_psyev(this%dynMatrix, this%denseDesc%blacsOrbSqr, this%eigen,&
          & this%eigenModesScaled, this%denseDesc%blacsOrbSqr, uplo="L", jobz="N")
    end if
  #:else
    ! remove translations or rotations if necessary
    call project(this%dynMatrix, this%tRemoveTranslate, this%tRemoveRotate, this%nDerivs,&
        & this%nMovedAtom, this%geo, this%atomicMasses)
    if (this%tEigenVectors) then
      call heev(this%dynMatrix, this%eigen, "U", "V")
    else
      call heev(this%dynMatrix, this%eigen, "U", "N")
    end if

    ! save original eigenvectors
    if (allocated(this%eigenModesScaled)) this%eigenModesScaled(:,:) = this%dynMatrix
  #:endif

    ! take square root of eigenvalues of modes (allowing for imaginary modes)
    this%eigen(:) = sign(sqrt(abs(this%eigen)), this%eigen)

    call TTaggedWriter_init(taggedWriter)
    call openFile(fd, "vibrations.tag", mode="w")

    ! ! scale mode components on each atom by mass and then normalise total mode
    ! do ii = 1, this%nDerivs
    !   jCount = 0
    !   do jj = 1, this%nMovedAtom
    !     do ll = 1, 3
    !       jCount = jCount + 1
    !       this%dynMatrix(jCount, ii) = this%dynMatrix(jCount, ii) / sqrt(this%atomicMasses(jj))
    !     end do
    !   end do
    !   this%dynMatrix(:, ii) = this%dynMatrix(:, ii) / sqrt(sum(this%dynMatrix(:, ii)**2))
    ! end do

    ! ! create displacement vectors for every atom in every mode.
    ! do iAt = 1, this%geo%nAtom
    !   if (any(this%iMovedAtoms == iAt)) then
    !     ! Index of atom in the list of moved atoms
    !     iAtMoved = minloc(abs(this%iMovedAtoms - iAt), 1)
    !     do ii = 1, this%nDerivs
    !       this%displ(:, iAt, ii) = this%dynMatrix(3 * iAtMoved - 2:3 * iAtMoved, ii)
    !     end do
    !   end if
    ! end do

    ! if (allocated(this%bornMatrix)) then
    !   allocate(transDip(this%nDerivs), source=0.0_dp)
    !   do jj = 1, this%nDerivs
    !     dMu(:) = 0.0_dp
    !     do ii = 1, this%nMovedAtom
    !       iAt = this%iMovedAtoms(ii)
    !       zStar(:,:) = reshape(this%bornMatrix(9 * (ii - 1) + 1:9 * ii), [3, 3])
    !       dMu(:) = dMu + matmul(zStar, this%displ(:, iAt, jj))
    !     end do
    !     if (this%eigen(jj) > epsilon(0.0_dp)) then
    !       transDip(jj) = transDip(jj) + sum(dMu**2)
    !     end if
    !   end do
    !   allocate(degenTransDip(this%nDerivs), source=0.0_dp)
    !   degenTransDip(1) = transDip(1)
    !   nTrans = 1
    !   do jj = 2, this%nDerivs
    !     ! test for energy degeneracy greater than printing cutoff:
    !     if (abs(this%eigen(jj) - this%eigen(jj - 1)) * Hartree__cm >= 1.0E-2_dp) then
    !       nTrans = nTrans + 1
    !     end if
    !     degenTransDip(nTrans) = degenTransDip(nTrans) + transDip(jj)
    !   end do
    ! end if

    ! if (allocated(this%bornDerivsMatrix)) then
    !   allocate(transPol(this%nDerivs), source=0.0_dp)
    !   do jj = 1, this%nDerivs
    !     dQ(:,:) = 0.0_dp
    !     do ii = 1, this%nMovedAtom
    !       iAt = this%iMovedAtoms(ii)
    !       zStarDeriv(:,:,:) = reshape(this%bornDerivsMatrix(27 * (ii - 1) + 1:27 * ii), [3, 3, 3])
    !       dQ(:,:) = dQ + reshape(matmul(reshape(zStarDeriv, [9, 3]), this%displ(:, iAt, jj)), [3, 3])
    !     end do
    !     if (this%eigen(jj) > epsilon(0.0_dp)) then
    !       transPol(jj) = transPol(jj) + sum(dQ**2)
    !     end if
    !   end do
    !   allocate(degenTransPol(this%nDerivs), source=0.0_dp)
    !   degenTransPol(1) = transPol(1)
    !   nTrans = 1
    !   do jj = 2, this%nDerivs
    !     ! test for energy degeneracy greater than printing cutoff:
    !     if (abs(this%eigen(jj) - this%eigen(jj - 1)) * Hartree__cm >= 1.0E-2_dp) then
    !       nTrans = nTrans + 1
    !     end if
    !     degenTransPol(nTrans) = degenTransPol(nTrans) + transPol(jj)
    !   end do
    ! end if

    ! if (this%tPlotModes) then
    !   call taggedWriter%write(fd%unit, "saved_modes", this%modesToPlot)
    !   write(stdout, *) "Writing eigenmodes to vibrations.tag"
    !   call taggedWriter%write(fd%unit, "eigenmodes", this%dynMatrix(:, this%modesToPlot))
    !   write(stdout, *) "Plotting eigenmodes:"
    !   write(stdout, "(16I5)") this%modesToPlot
    !   call taggedWriter%write(fd%unit, "eigenmodes_scaled",&
    !       & this%eigenModesScaled(:, this%modesToPlot))
      ! if (this%tAnimateModes) then
      !   do ii = 1, this%nModesToPlot
      !     iMode = this%modesToPlot(ii)
      !     write(lcTmp,"('mode_',I0,'.xyz')") iMode
      !     do kk = 1, this%nCycles
      !       do ll = 1, this%nSteps
      !         isAppend = (kk > 1 .or. ll > 1)
      !         write(lcTmp2, *) "Eigenmode", iMode, this%eigen(iMode) * Hartree__cm, "cm-1"
      !         call writeXYZFormat(lcTmp, this%geo%coords + cos(2.0_dp * pi * real(ll)&
      !             & / real(this%nSteps)) * this%displ(:,:, iMode), this%geo%species,&
      !             & this%geo%speciesNames, comment=trim(lcTmp2), append=isAppend)
      !       end do
      !     end do
      !   end do
      ! else
      !   lcTmp = "modes.xyz"
      !   do ii = 1, this%nModesToPlot
      !     isAppend = (ii > 1)
      !     iMode = this%modesToPlot(ii)
      !     write(lcTmp2, *) "Eigenmode", iMode, this%eigen(iMode) * Hartree__cm, "cm-1"
      !     call writeXYZFormat(lcTmp, this%geo%coords, this%geo%species, this%geo%speciesNames,&
      !         & vectors=this%displ(:,:, iMode), comment=trim(lcTmp2), append=isAppend)
      !   end do
      ! end if
    ! end if

    ! write(stdout, *) "Vibrational modes"
    ! if (allocated(this%bornMatrix) .and. allocated(this%bornDerivsMatrix)) then
    !   write(stdout, "(T7,A,T16,A,T28,A)") "freq.", "IR", "Polarisability"
    !   write(stdout, "(A,T7,A,T16,A,T28,A)") "Mode", "/ cm-1", "/ a.u.", "change / a.u."
    ! else if (allocated(this%bornMatrix)) then
    !   write(stdout, "(T7,A,T16,A)") "freq.", "IR"
    !   write(stdout, "(A,T7,A,T16,A)") "Mode", "/ cm-1", "/ a.u."
    ! else if (allocated(this%bornDerivsMatrix)) then
    !   write(stdout, "(T7,A,T16,A)") "freq.", "Polarisability"
    !   write(stdout, "(A,T7,A,T16,A)") "Mode", "/ cm-1", "change / a.u."
    ! else
    !   write(stdout, "(T7,A)") "freq."
    !   write(stdout, "(A,T7,A)") "Mode", "cm-1"
    ! end if
    ! if (allocated(this%bornMatrix) .and. allocated(this%bornDerivsMatrix)) then
    !   do ii = 1, this%nDerivs
    !     write(stdout, "(i5,f8.2,2E12.4)") ii, this%eigen(ii) * Hartree__cm, transDip(ii),&
    !         & transPol(ii)
    !   end do
    ! else if (allocated(this%bornMatrix)) then
    !   do ii = 1, this%nDerivs
    !     write(stdout, "(i5,f8.2,E12.4)") ii, this%eigen(ii) * Hartree__cm, transDip(ii)
    !   end do
    ! else if (allocated(this%bornDerivsMatrix)) then
    !   do ii = 1, this%nDerivs
    !     write(stdout, "(i5,f8.2,E12.4)") ii, this%eigen(ii) * Hartree__cm, transPol(ii)
    !   end do
    ! else
    !   do ii = 1, this%nDerivs
    !     write(stdout, "(i5,f8.2)") ii, this%eigen(ii) * Hartree__cm
    !   end do
    ! end if
    ! write(stdout, *)

    call taggedWriter%write(fd%unit, "frequencies", this%eigen)

    ! if (allocated(this%bornMatrix)) then
    !   call taggedWriter%write(fd%unit, "intensities", degenTransDip(:nTrans))
    ! end if

    ! if (allocated(this%bornDerivsMatrix)) then
    !   if (this%tRemoveTranslate .or. this%tRemoveRotate) then
    !     call taggedWriter%write(fd%unit, "scattering", degenTransPol(2:nTrans))
    !   else
    !     call taggedWriter%write(fd%unit, "scattering", degenTransPol(:nTrans))
    !   end if
    ! end if

    call closeFile(fd)

  end subroutine runModes

end module modes_main
