!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> REKS and SI-SA-REKS formulation in DFTB as developed by Lee et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet appart from Gamma point.
!> * Orbital potentials or spin-orbit or external E-field does not work yet.
!> * Only for closed shell system.
!> * Onsite corrections are not included in this version
module dftbp_reksinterface

  use dftbp_accuracy
  use dftbp_densedescr
  use dftbp_dispiface
  use dftbp_environment
  use dftbp_globalenv
  use dftbp_mainio
  use dftbp_nonscc
  use dftbp_orbitals
  use dftbp_periodic
  use dftbp_populations
  use dftbp_rangeseparated
  use dftbp_repcont
  use dftbp_repulsive
  use dftbp_scc
  use dftbp_slakocont
  use dftbp_sparse2dense
  use dftbp_stress
  use dftbp_taggedoutput, only : TTaggedWriter, tagLabels
  use dftbp_rekscommon
  use dftbp_rekscpeqn
  use dftbp_reksgrad
  use dftbp_reksproperty
  use dftbp_reksvar

  implicit none

  private

  public :: getReksGradients, getReksGradProperties
  public :: getReksStress

  contains

  !> Calculate SI-SA-REKS state gradient by solving CP-REKS equations
  subroutine getReksGradients(env, denseDesc, sccCalc, rangeSep, dispersion, &
      & neighbourList, nNeighbourSK, nNeighbourRep, iSparseStart, img2CentCell, &
      & orb, nonSccDeriv, skHamCont, skOverCont, pRepCont, coord, coord0, &
      & species, q0, eigenvecs, chrgForces, over, spinW, derivs, tWriteTagged, &
      & autotestTag, taggedWriter, self)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> Range separation contributions
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> dispersion interactions
    class(TDispersionIface), allocatable, intent(inout) :: dispersion

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Number of neighbours for each of the atoms closer than the repulsive cut-off
    integer, intent(in) :: nNeighbourRep(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> method for calculating derivatives of S and H0
    type(TNonSccDiff), intent(in) :: nonSccDeriv

    !> non-SCC hamiltonian information
    type(TSlakoCont), intent(in) :: skHamCont

    !> overlap information
    type(TSlakoCont), intent(in) :: skOverCont

    !> repulsive information
    type(TRepCont), intent(in) :: pRepCont

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> central cell coordinates of atoms
    real(dp), intent(inout) :: coord0(:,:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> forces on external charges
    real(dp), allocatable, intent(inout) :: chrgForces(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> spin constants
    real(dp), intent(in) :: spinW(:,:,:)

    !> derivatives of energy wrt to atomic positions
    real(dp), intent(out) :: derivs(:,:)

    !> print tag information
    logical, intent(in) :: tWriteTagged

    !> File name for regression data
    character(*), intent(in) :: autotestTag

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    real(dp), allocatable :: Qmat(:,:)
    integer :: ist, ia, ib, nstHalf, fac

    nstHalf = self%nstates * (self%nstates-1) / 2

    allocate(Qmat(orb%nOrb,orb%nOrb))

    ! get the periodic information
    if (self%tPeriodic) then
      call sccCalc%coulombCont%getPeriodicInfo(self%rVec, self%gVec, self%alpha, self%volume)
    end if

    call getHellmannFeynmanGradientL_(env, denseDesc, sccCalc, neighbourList, &
        & nNeighbourSK, nNeighbourRep, iSparseStart, img2CentCell, orb, &
        & nonSccDeriv, skHamCont, skOverCont, pRepCont, coord, species, q0, &
        & dispersion, rangeSep, chrgForces, eigenvecs, derivs, self)

    if (self%Efunction == 1) then
      call weightGradient(self%gradL, self%weight, derivs)
    else

      ! get REKS parameters used in CP-REKS and gradient equations
      call getReksParameters_(env, denseDesc, sccCalc, rangeSep, &
          & neighbourList, nNeighbourSK, iSparseStart, img2CentCell, &
          & eigenvecs, coord, species, over, spinW, self)

      call buildStateVectors_(env, denseDesc, neighbourList, nNeighbourSK, &
          & iSparseStart, img2CentCell, eigenvecs, orb, over, self)

      if (self%tNAC) then

        ! SA-REKS calculations
        do ist = 1, self%nstates
          if (ist /= self%SAstates) then

            call printBlankLine()
            write(stdOut,"(A)") repeat("-", 82)
            write(stdOut,'(1x,a,1x,I2,1x,a)') &
                & 'Solving CP-REKS equation for', ist, 'state vector...'

            ! solve CP-REKS equation for SA-REKS state
            ! save information about ZT, RmatL
            call solveCpReks_(env, denseDesc, neighbourList, nNeighbourSK, &
                & iSparseStart, img2CentCell, eigenvecs, over, orb, self, &
                & self%XT(:,ist), self%ZT(:,ist), self%RmatL(:,:,:,ist), &
                & self%ZmatL, self%Q1mat, self%Q2mat, .false.)
            Qmat(:,:) = self%Q1mat + self%Q2mat
            ! compute SA-REKS shift
            call SSRshift(eigenvecs, self%gradL, Qmat, self%Sderiv, &
                & self%ZT(:,ist), self%SAweight, self%weightL(ist,:), &
                & self%omega, self%weightIL, self%G1, &
                & denseDesc%iAtomStart, orb%mOrb, self%SAgrad(:,:,ist), 1)

          end if
        end do

        ! state-interaction calculations
        do ist = 1, nstHalf

          call getTwoIndices(self%nstates, ist, ia, ib, 1)
          call printBlankLine()
          write(stdOut,"(A)") repeat("-", 82)
          write(stdOut,'(1x,a,1x,I2,1x,a,1x,I2,1x,a)') &
              & 'Solving CP-REKS equation for SI between', ia, 'and', ib, 'state vectors...'

          ! solve CP-REKS equation for state-interaction term
          ! save information about ZTdel, tmpRL
          call solveCpReks_(env, denseDesc, neighbourList, nNeighbourSK, &
              & iSparseStart, img2CentCell, eigenvecs, over, orb, self, &
              & self%XTdel(:,ist), self%ZTdel(:,ist), self%tmpRL(:,:,:,ist), &
              & self%ZmatL, self%Q1mat, self%Q2mat, .false.)
          call SIshift(eigenvecs, self%gradL, self%Q1del(:,:,ist), &
              & self%Q2del(:,:,ist), self%Q1mat, self%Q2mat, Qmat, &
              & self%Sderiv, self%ZTdel(:,ist), self%SAweight, self%omega, &
              & self%weightIL, self%Rab, self%G1, denseDesc%iAtomStart, &
              & orb%mOrb, ist, self%nstates, self%SIgrad(:,:,ist))

        end do

      else

        if (self%useSSR == 1) then

          ! Convert necessary variables for SSR state from SA and SI terms
          call SaToSsrXT(self%XTdel, self%eigvecsSSR, self%rstate, self%XT)
          call SaToSsrWeight(self%Rab, self%weightIL, self%G1, &
              & self%eigvecsSSR, self%rstate, self%weightL)

          call printBlankLine()
          write(stdOut,"(A)") repeat("-", 82)
          write(stdOut,'(1x,a,1x,I2,1x,a)') &
              & 'Solving CP-REKS equation for', self%rstate, 'state vector...'

          ! solve CP-REKS equation for SSR state
          ! save information about ZT, RmatL, ZmatL, Q1mat, Q2mat
          call solveCpReks_(env, denseDesc, neighbourList, nNeighbourSK, &
              & iSparseStart, img2CentCell, eigenvecs, over, orb, self, &
              & self%XT(:,self%rstate), self%ZT(:,1), self%RmatL(:,:,:,1), &
              & self%ZmatL, self%Q1mat, self%Q2mat, .false.)

          ! add remaining SI component to SSR state
          call addSItoRQ(eigenvecs, self%RdelL, self%Q1del, self%Q2del, &
              & self%eigvecsSSR, self%rstate, self%RmatL(:,:,:,1), &
              & self%Q1mat, self%Q2mat, Qmat)
          fac = 2

        else

          call printBlankLine()
          write(stdOut,"(A)") repeat("-", 82)
          if (self%Lstate == 0) then
            write(stdOut,'(1x,a,1x,I2,1x,a)') &
                & 'Solving CP-REKS equation for', self%rstate, 'state vector...'
          else
            write(stdOut,'(1x,a,1x,I2,1x,a)') &
                & 'Solving CP-REKS equation for', self%Lstate, 'microstate vector...'
          end if

          ! solve CP-REKS equation for SA-REKS or L state
          ! save information about ZT, RmatL, ZmatL, Q1mat, Q2mat
          call solveCpReks_(env, denseDesc, neighbourList, nNeighbourSK, &
              & iSparseStart, img2CentCell, eigenvecs, over, orb, self, &
              & self%XT(:,1), self%ZT(:,1), self%RmatL(:,:,:,1), &
              & self%ZmatL, self%Q1mat, self%Q2mat, .false.)
          Qmat(:,:) = self%Q1mat + self%Q2mat
          fac = 1

        end if

        if (self%Lstate == 0) then
          ! compute SSR or SA-REKS shift
          call SSRshift(eigenvecs, self%gradL, Qmat, self%Sderiv, &
              & self%ZT(:,1), self%SAweight, self%weightL(self%rstate,:), &
              & self%omega, self%weightIL, self%G1, &
              & denseDesc%iAtomStart, orb%mOrb, self%SSRgrad(:,:,1), fac)
        else
          ! compute L shift
          call Lshift(eigenvecs, self%gradL, Qmat, self%Sderiv, &
              & self%ZT(:,1), self%SAweight, self%omega, self%weightIL, self%G1, &
              & denseDesc%iAtomStart, orb%mOrb, self%Lstate, self%SSRgrad(:,:,1))
        end if

      end if

      ! compute R*T shift & final SSR gradient (grad, nac)
      call getRTGradient_(env, sccCalc, denseDesc, neighbourList, &
          & nNeighbourSK, iSparseStart, img2CentCell, orb, q0, coord0, self)
      if (self%tNAC) then
        call getReksNACinfo_(tWriteTagged, autotestTag, taggedWriter, self)
      end if

      ! set final gradient to derivs in SA-REKS or SSR case
      call setReksGradients_(derivs, self)

    end if

    if (self%Plevel >= 1) then
      call printReksGradInfo(self, derivs)
    end if

  end subroutine getReksGradients


  !> Get the gradient-related properties such as relaxed density, charge, etc
  subroutine getReksGradProperties(env, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, eigenvecs, orb, &
      & iAtInCentralRegion, coord, coord0, over, rhoPrim, qOutput, &
      & q0, tDipole, dipoleMoment, chrgForces, self)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:,:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Atoms over which to sum the total energies
    integer, intent(in) :: iAtInCentralRegion(:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> central cell coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> sparse density matrix
    real(dp), intent(inout) :: rhoPrim(:,:)

    !> Output electrons
    real(dp), intent(inout) :: qOutput(:,:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> Print the dipole moment
    logical, intent(in) :: tDipole

    !> resulting dipole moment
    real(dp), allocatable, intent(inout) :: dipoleMoment(:)

    !> forces on external charges
    real(dp), allocatable, intent(inout) :: chrgForces(:,:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    ! calculate the relaxed density matrix, dipole moment,
    !           the forces of external charges.
    if (self%tRD) then

      if (self%Efunction > 1) then

        ! currently, relaxed density matrix is possible
        !            relaxed charge is possible
        !            relaxed dipole moment is possible
        !            force of external point charges for target state is possible
        if (self%tNAC) then
          ! solve CP-REKS equation for remaining SA-REKS state in nac case
          ! save information about Z_T
          call solveCpReks_(env, denseDesc, neighbourList, nNeighbourSK, &
              & iSparseStart, img2CentCell, eigenvecs, over, orb, self, &
              & self%XT(:,self%SAstates), self%ZT(:,self%SAstates), &
              & self%RmatL(:,:,:,self%SAstates), self%ZmatL, &
              & self%Q1mat, self%Q2mat, .true.)
          ! now, ZT has information about target SSR state
          call SaToSsrXT(self%ZTdel, self%eigvecsSSR, self%rstate, self%ZT)
        end if

        if (self%Lstate == 0) then
          ! get the relaxed density matrix for target SSR or SA-REKS state
          call getRelaxedDensMat(eigenvecs(:,:,1), self%overSqr, self%unrelRhoSqr, &
              & self%ZT, self%omega, self%FONs, self%eigvecsSSR, self%SAweight, &
              & self%Rab, self%G1, self%Nc, self%Na, self%rstate, self%useSSR, &
              & self%tNAC, self%tSSR22, self%tSSR44, self%relRhoSqr)
        else
          ! get the relaxed density matrix for L-th microstate
          call getRelaxedDensMatL(eigenvecs(:,:,1), self%rhoSqrL, self%overSqr, &
              & self%weight, self%SAweight, self%unrelRhoSqr, self%RmatL, &
              & self%ZT, self%omega, self%weightIL, self%G1, self%orderRmatL, &
              & self%Lpaired, self%Nc, self%Na, self%Lstate, self%tSSR22, &
              & self%tSSR44, self%relRhoSqr)
        end if

        call getMullikenPopFromRelaxedDensity_(env, denseDesc, neighbourList, &
            & nNeighbourSK, iSparseStart, img2CentCell, orb, over, rhoPrim, qOutput, self)

      end if

      call writeReksRelaxedCharge(qOutput, q0, self%rstate, self%Lstate)

      if (self%Efunction > 1) then
        ! Get the dipole moment for target SSR state
        ! Relaxed dipole moment requires gradient infomation
        ! tDipole = (total charge = 0.0) * (non-periodic system) * (mulliken)
        if (tDipole) then
          call getDipoleMoment_(qOutput, q0, coord, dipoleMoment, iAtInCentralRegion)
        end if
      end if

      if (self%tExtChrg) then

        call getExtChrgGradients(env, coord0, self%extCharges(1:3,:), &
            & qOutput, q0, self%extCharges(4,:), self%blurWidths, self%rVec, self%gVec, &
            & self%alpha, self%volume, self%tPeriodic, self%tBlur, chrgForces)

      end if

    end if

  end subroutine getReksGradProperties


  !> Calculates stress tensor and lattice derivatives.
  subroutine getReksStress(env, denseDesc, sccCalc, nonSccDeriv, &
      & skHamCont, skOverCont, pRepCont, neighbourList, nNeighbourSk, &
      & nNeighbourRep, species, img2CentCell, iSparseStart, orb, &
      & dispersion, coord, q0, invLatVec, cellVol, totalStress, &
      & totalLatDeriv, intPressure, self)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

!    !> Is 3rd order SCC being used
!    type(ThirdOrder), intent(inout), allocatable :: thirdOrd

    !> method for calculating derivatives of S and H0
    type(TNonSccDiff), intent(in) :: nonSccDeriv

    !> non-SCC hamitonian information
    type(TSlakoCont), intent(in) :: skHamCont

    !> overlap information
    type(TSlakoCont), intent(in) :: skOverCont

    !> repulsive information
    type(TRepCont), intent(in) :: pRepCont

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Number of neighbours for each of the atoms closer than the repulsive cut-off
    integer, intent(in) :: nNeighbourRep(:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> dispersion interactions
    class(TDispersionIface), allocatable, intent(inout) :: dispersion

    !> coordinates of all atoms
    real(dp), intent(in) :: coord(:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> inverse of the lattice vectors
    real(dp), intent(in) :: invLatVec(:,:)

    !> unit cell volume
    real(dp), intent(in) :: cellVol

    !> stress tensor
    real(dp), intent(out) :: totalStress(:,:)

    !> energy derivatives with respect to lattice vectors
    real(dp), intent(out) :: totalLatDeriv(:,:)

    !> internal pressure in cell
    real(dp), intent(out) :: intPressure

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    real(dp), allocatable :: tmpRhoSp(:,:)
    real(dp), allocatable :: tmpStress(:,:)

    real(dp) :: fac
    integer :: sparseSize, nSpin, tmpLu, tmpLd, iL

    sparseSize = size(self%getDenseAO,dim=1)
    nSpin = size(self%fillingL,dim=2)

    if (self%Efunction == 1) then

      allocate(tmpRhoSp(sparseSize,nSpin))
      allocate(tmpStress(3,3))

      do iL = 1, self%Lmax

        if (iL <= self%Lpaired) then
          tmpLu = iL
        else
          if (mod(iL,2) == 1) then
            tmpLu = iL
            tmpLd = iL + 1
            fac = 1.0
          else
            tmpLu = iL - 1
            tmpLd = iL
            fac = -1.0
          end if
        end if

        ! self%qOutputL has (qm) component
        call sccCalc%updateCharges(env, self%qOutputL(:,:,:,iL), &
            & q0, orb, species)
        call sccCalc%updateShifts(env, orb, species, &
            & neighbourList%iNeighbour, img2CentCell)

        call env%globalTimer%startTimer(globalTimers%denseToSparse)
        ! self%rhoSqrL has (my_qm) component
        tmpRhoSp(:,:) = 0.0_dp
        call packHS(tmpRhoSp(:,1), self%rhoSqrL(:,:,1,tmpLu), &
            & neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb, &
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        if (iL > self%Lpaired) then
          call packHS(tmpRhoSp(:,2), fac*self%rhoSqrL(:,:,1,tmpLd), &
              & neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb, &
              & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
        call env%globalTimer%stopTimer(globalTimers%denseToSparse)

        tmpStress(:,:) = 0.0_dp
        call getBlockStress(env, tmpStress, nonSccDeriv, tmpRhoSp, &
            & self%edmSpL(:,iL), skHamCont, skOverCont, coord, &
            & species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell, &
            & iSparseStart, orb, self%intBlockL(:,:,:,:,iL), cellVol)
        call sccCalc%addStressDc(tmpStress, env, species, &
            & neighbourList%iNeighbour, img2CentCell)

!        if (allocated(thirdOrd)) then
!          call thirdOrd%updateCharges(pSpecies0, neighbourList, &
!              & reks%qOutput_L(:,:,:,iL), reks%q0, img2CentCell, orb)
!          call thirdOrd%addStressDc(neighbourList, species, coord, &
!              & img2CentCell, cellVol, tmpStress)
!        end if

        self%elecStressL(:,:,iL) = tmpStress

      end do

      call weightGradient(self%elecStressL, self%weight, totalStress)

      if (allocated(dispersion)) then
        tmpStress(:,:) = 0.0_dp
        call dispersion%getStress(tmpStress)
        totalStress(:,:) = totalStress + tmpStress
      end if

      tmpStress(:,:) = 0.0_dp
      call getRepulsiveStress(tmpStress, coord, nNeighbourRep, &
          & neighbourList%iNeighbour, species, img2CentCell, pRepCont, cellVol)
      totalStress(:,:) = totalStress + tmpStress

      intPressure = (totalStress(1,1) + totalStress(2,2) + totalStress(3,3)) / 3.0_dp
      totalLatDeriv(:,:) = -cellVol * matmul(totalStress, invLatVec)

    end if

  end subroutine getReksStress


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate Hellmann-Feynman gradient term of each microstate in REKS
  subroutine getHellmannFeynmanGradientL_(env, denseDesc, sccCalc, neighbourList, &
      & nNeighbourSK, nNeighbourRep, iSparseStart, img2CentCell, orb, &
      & nonSccDeriv, skHamCont, skOverCont, pRepCont, coord, species, q0, &
      & dispersion, rangeSep, chrgForces, eigenvecs, derivs, self)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Number of neighbours for each of the atoms closer than the repulsive cut-off
    integer, intent(in) :: nNeighbourRep(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> method for calculating derivatives of S and H0
    type(TNonSccDiff), intent(in) :: nonSccDeriv

    !> non-SCC hamiltonian information
    type(TSlakoCont), intent(in) :: skHamCont

    !> overlap information
    type(TSlakoCont), intent(in) :: skOverCont

    !> repulsive information
    type(TRepCont), intent(in) :: pRepCont

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> dispersion interactions
    class(TDispersionIface), allocatable, intent(inout) :: dispersion

    !> Range separation contributions
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> forces on external charges
    real(dp), allocatable, intent(inout) :: chrgForces(:,:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> derivatives of energy wrt to atomic positions
    real(dp), intent(out) :: derivs(:,:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    real(dp), allocatable :: repDerivs(:,:)
    real(dp), allocatable :: dispDerivs(:,:)
    real(dp), allocatable :: lcDerivs(:,:,:)
!    integer, pointer :: pSpecies0(:)

    integer :: nAtom, iL

    nAtom = size(self%gradL,dim=2)
!    pSpecies0 => species(1:nat)

    allocate(repDerivs(3,nAtom))
    allocate(dispDerivs(3,nAtom))
    if (self%tRangeSep) then
      allocate(lcDerivs(3,nAtom,self%Lmax))
    end if

    ! hamSpL has (my_ud) component
    ! hamSqrL has (my_ud) component
    call env%globalTimer%startTimer(globalTimers%energyDensityMatrix)
    call getEnergyWeightedDensityL(env, denseDesc, neighbourList, &
        & nNeighbourSK, iSparseStart, img2CentCell, orb, self%hamSqrL, &
        & self%hamSpL, self%fillingL, eigenvecs(:,:,1), self%Lpaired, &
        & self%Efunction, self%tRangeSep, self%edmSpL)
    call env%globalTimer%stopTimer(globalTimers%energyDensityMatrix)

    ! rhoSpL has (my_qm) component
    ! intBlockL has (qm) component
    call derivative_blockL(env, derivs, nonSccDeriv, self%rhoSqrL, &
        & self%edmSpL, skHamCont, skOverCont, coord, species, &
        & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, &
        & iSparseStart, denseDesc%iAtomStart, orb, self%intBlockL, &
        & self%Lpaired, self%gradL, self%Hderiv, self%Sderiv)

    do iL = 1, self%Lmax

      derivs(:,:) = 0.0_dp

      ! qOutput_L has (qm) component
      call sccCalc%updateCharges(env, self%qOutputL(:,:,:,iL), &
          & q0, orb, species)
      call sccCalc%updateShifts(env, orb, species, &
          & neighbourList%iNeighbour, img2CentCell)
      if (self%tExtChrg) then
        chrgForces(:,:) = 0.0_dp
        call sccCalc%addForceDc(env, derivs, species, &
            & neighbourList%iNeighbour, img2CentCell, chrgForces)
      else
        call sccCalc%addForceDc(env, derivs, species, &
            & neighbourList%iNeighbour, img2CentCell)
      endif

!      if (allocated(thirdOrd)) then
!        call thirdOrd%updateCharges(pSpecies0, neighbourList, &
!            & self%qOutputL(:,:,:,iL), q0, img2CentCell, orb)
!        call thirdOrd%addGradientDc(neighbourList, species, &
!            & coord, img2CentCell, derivs)
!      end if

      if (self%tRangeSep) then
        ! deltaRhoSqrL has (my_ud) component
        lcDerivs(:,:,iL) = 0.0_dp
        call rangeSep%addLRGradients(lcDerivs(:,:,iL), nonSccDeriv, &
            & self%deltaRhoSqrL(:,:,:,iL), skHamCont, skOverCont, coord, &
            & species, orb, denseDesc%iAtomStart, self%overSqr, &
            & neighbourList%iNeighbour, nNeighbourSK)
      end if

      if (iL == 1) then

        dispDerivs(:,:) = 0.0_dp
        if (allocated(dispersion)) then
          call dispersion%addGradients(dispDerivs)
        end if

        repDerivs(:,:) = 0.0_dp
        call getERepDeriv(repDerivs, coord, nNeighbourRep, &
            & neighbourList%iNeighbour, species, pRepCont, img2CentCell)

      end if
      derivs(:,:) = derivs + repDerivs + dispDerivs

      self%gradL(:,:,iL) = self%gradL(:,:,iL) + derivs

    end do

    if(self%tRangeSep) then
      do iL = 1, self%Lmax
        if (iL <= self%Lpaired) then
          derivs(:,:) = lcDerivs(:,:,iL) + lcDerivs(:,:,iL)
        else
          if (mod(iL,2) == 1) then
            derivs(:,:) = lcDerivs(:,:,iL) + lcDerivs(:,:,iL+1)
          else
            derivs(:,:) = lcDerivs(:,:,iL) + lcDerivs(:,:,iL-1)
          end if
        end if
        self%gradL(:,:,iL) = self%gradL(:,:,iL) + derivs
      end do
    end if

  end subroutine getHellmannFeynmanGradientL_


  !> Set several REKS variables used in CP-REKS equations
  subroutine getReksParameters_(env, denseDesc, sccCalc, rangeSep, &
      & neighbourList, nNeighbourSK, iSparseStart, img2CentCell, &
      & eigenvecs, coord, species, over, spinW, self)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> Range separation contributions
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> spin constants
    real(dp), intent(in) :: spinW(:,:,:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    ! get gamma, spinW, gamma deriv, LR-gamma deriv, on-site constants
    ! LRgamma is already defined in 1st scc loop
    call getSccSpinLrPars(env, sccCalc, rangeSep, coord, species, &
        & neighbourList%iNeighbour, img2CentCell, denseDesc%iAtomStart, &
        & spinW, self%getAtomIndex, self%tRangeSep, self%GammaAO, &
        & self%GammaDeriv, self%SpinAO, self%LrGammaAO, self%LrGammaDeriv)

    ! get Hxc kernel -> (\mu,\nu|f_{Hxc}|\tau,\gam)
    call getHxcKernel(denseDesc%iAtomStart, self%getAtomIndex, self%getDenseAO, &
        & over, self%overSqr, self%GammaAO, self%SpinAO, self%LrGammaAO, &
        & self%tRangeSep, self%Glevel, self%Mlevel, self%HxcSpS, &
        & self%HxcSpD, self%HxcHalfS, self%HxcHalfD, self%HxcSqrS, self%HxcSqrD)

    ! get G1, weightIL, Omega, Rab values
    call getG1ILOmegaRab(env, denseDesc, neighbourList, nNeighbourSK, &
        & iSparseStart, img2CentCell, eigenvecs, self%hamSqrL, self%hamSpL, &
        & self%fockFa, self%fillingL, self%FONs, self%SAweight, self%enLtot, &
        & self%hess, self%Nc, self%Na, self%useSSR, self%tSSR22, self%tSSR44, &
        & self%tRangeSep, self%G1, self%weightIL, self%omega, self%Rab)

    ! get A1e or Aall values based on GradOpt
    call getSuperAMatrix(eigenvecs, self%HxcSqrS, self%HxcSqrD, self%fockFc, &
        & self%fockFa, self%omega, self%fillingL, self%weight, self%SAweight, &
        & self%FONs, self%G1, self%Lpaired, self%Nc, self%Na, self%Glevel, &
        & self%Mlevel, self%tSSR22, self%tSSR44, self%A1e, self%A1ePre, self%Aall)

  end subroutine getReksParameters_


  !> Calculate SA-REKS state vectors and state-interaction vectors
  subroutine buildStateVectors_(env, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, eigenvecs, orb, over, self)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    integer :: ia, ib, ist, nstHalf

    nstHalf = self%nstates * (self%nstates - 1) / 2

    if (self%Lstate == 0) then

      ! build X^T value for state X = PPS, OSS, etc
      call buildSaReksVectors(env, denseDesc, neighbourList, nNeighbourSK, &
          & iSparseStart, img2CentCell, eigenvecs, self%hamSqrL, self%hamSpL, &
          & self%fillingL, self%weightL, self%Nc, self%Na, self%rstate, &
          & self%useSSR, self%tSSR22, self%tSSR44, self%tRangeSep, self%XT)

      if (self%useSSR == 1) then

        ! get R^delta values between states X's
        call getRdel(eigenvecs, self%fillingL, self%FONs, self%Nc, &
            & self%nstates, self%tSSR22, self%tSSR44, self%RdelL)

        do ist = 1, nstHalf

          call getTwoIndices(self%nstates, ist, ia, ib, 1)

          ! get Z^delta values from R^delta values
          call getZmat(env, denseDesc, neighbourList, nNeighbourSK, &
              & iSparseStart, img2CentCell, orb, self%RdelL(:,:,:,ist), &
              & self%HxcSqrS, self%HxcSqrD, self%HxcHalfS, self%HxcHalfD, &
              & self%HxcSpS, self%HxcSpD, self%overSqr, over, self%GammaAO, &
              & self%SpinAO, self%LrGammaAO, self%orderRmatL, self%getDenseAO, &
              & self%Lpaired, self%Glevel, self%Mlevel, self%tRangeSep, self%ZdelL)

          ! build XTdel with Z^delta values
          call buildInteractionVectors(eigenvecs, self%ZdelL, self%fockFc, &
              & self%fockFa, self%FONs, self%fillingL, self%weight, &
              & self%SAweight, self%omega, self%Rab, self%G1, self%Nc, self%Na, &
              & ia, ib, self%tSSR22, self%tSSR44, self%XTdel(:,ist))

          ! get Q2^delta values from Z^delta values : MO index
          call getQ2mat(eigenvecs, self%fillingL, self%weight, &
              & self%ZdelL, self%Q2del(:,:,ist))

        end do

        ! get Q1^delta values : AO index
        call getQ1del(eigenvecs, self%fockFc, self%fockFa, self%FONs, self%SAweight, &
            & self%Nc, self%nstates, self%tSSR22, self%tSSR44, self%Q1del)

      end if

    else

      ! build X^T_L value
      call buildLstateVector(env, denseDesc, neighbourList, nNeighbourSK, &
          & iSparseStart, img2CentCell, eigenvecs, self%hamSqrL, self%hamSpL, &
          & self%fillingL, self%Nc, self%Na, self%Lstate, self%Lpaired, &
          & self%tSSR22, self%tSSR44, self%tRangeSep, self%XT(:,1))

    end if

  end subroutine buildStateVectors_


  !> Solve CP-REKS equations by using CG based algorithms or direct matrix multiplication
  subroutine solveCpReks_(env, denseDesc, neighbourList, nNeighbourSK, &
      & iSparseStart, img2CentCell, eigenvecs, over, orb, self, XT, ZT, &
      & RmatL, ZmatL, Q1mat, Q2mat, optionQMMM)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    !> SA-REKS state vector
    real(dp), intent(in) :: XT(:)

    !> solution of A * Z = X equation with X is XT
    real(dp), intent(inout) :: ZT(:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(inout) :: RmatL(:,:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(inout) :: ZmatL(:,:,:)

    !> auxiliary matrix in MO basis related to SA-REKS term
    real(dp), intent(inout) :: Q1mat(:,:)

    !> auxiliary matrix in MO basis related to SA-REKS term
    real(dp), intent(inout) :: Q2mat(:,:)

    !> option for relaxed properties (QM/MM) calculations
    logical, intent(in) :: optionQMMM

    if (self%Glevel == 1 .or. self%Glevel == 2) then

      ! solve linear equation A * Z = X using CG algorithm
      ! get converged R, Z, Q2 matrices & ZT vector
      ! RmatL : AO index, ZmatL : AO index, Q2mat : MO index
      if (optionQMMM) then
        write(stdOut,'(a)') &
            & ' Warning! For calculating relaxed density for SSR state,'
        write(stdOut,'(a)') &
            & '          run CG again to obtain ZT solution in (nac) case.'
      end if
      call CGgrad(env, denseDesc, neighbourList, nNeighbourSK, iSparseStart, &
          & img2CentCell, orb, XT, self%A1e, self%A1ePre, self%HxcSqrS, &
          & self%HxcSqrD, self%HxcHalfS, self%HxcHalfD, self%HxcSpS, self%HxcSpD, &
          & self%fockFc, self%fockFa, self%omega, self%SAweight, self%FONs, &
          & self%G1, self%GammaAO, self%SpinAO, self%LrGammaAO, self%overSqr, &
          & over, eigenvecs, self%fillingL, self%weight, self%Glimit, self%orderRmatL, &
          & self%getDenseAO, self%Lpaired, self%Nc, self%Na, self%CGmaxIter, self%Glevel, &
          & self%Mlevel, self%tRangeSep, self%tSSR22, self%tSSR44, ZT, RmatL, ZmatL, Q2mat)

    else if (self%Glevel == 3) then

      ! solve A * Z = X using direct matrix inversion
      call solveZT(self%Aall, XT, ZT)

      ! get direct R, Z, Q2 matrices
      if (.not. optionQMMM) then
        call getRmat(eigenvecs, ZT, self%fillingL, self%Nc, self%Na, &
            & self%tSSR22, self%tSSR44, RmatL)
        call getZmat(env, denseDesc, neighbourList, nNeighbourSK, &
            & iSparseStart, img2CentCell, orb, RmatL, &
            & self%HxcSqrS, self%HxcSqrD, self%HxcHalfS, self%HxcHalfD, &
            & self%HxcSpS, self%HxcSpD, self%overSqr, over, self%GammaAO, &
            & self%SpinAO, self%LrGammaAO, self%orderRmatL, self%getDenseAO, &
            & self%Lpaired, self%Glevel, self%Mlevel, self%tRangeSep, ZmatL)
        call getQ2mat(eigenvecs, self%fillingL, self%weight, ZmatL, Q2mat)
        write(stdOut,"(A)") repeat("-", 82)
      end if

    end if

    ! get Q1 matrix from converged ZT vector : MO index
    if (.not. optionQMMM) then
      call getQ1mat(ZT, self%fockFc, self%fockFa, self%SAweight, &
          & self%FONs, self%Nc, self%Na, self%tSSR22, self%tSSR44, Q1mat)
    else
      call printBlankLine()
    end if

  end subroutine solveCpReks_


  !> Calculate R*T contribution of gradient (2nd term)
  subroutine getRTGradient_(env, sccCalc, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, orb, q0, coord0, self)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> central cell coordinates of atoms
    real(dp), intent(inout) :: coord0(:,:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    call RTshift(env, sccCalc, denseDesc, neighbourList, nNeighbourSK, &
        & iSparseStart, img2CentCell, orb, coord0, self%Hderiv, self%Sderiv, &
        & self%rhoSqrL, self%overSqr, self%deltaRhoSqrL, self%qOutputL, &
        & q0, self%GammaAO, self%GammaDeriv, self%SpinAO, self%LrGammaAO, &
        & self%LrGammaDeriv, self%RmatL, self%RdelL, self%tmpRL, self%weight, &
        & self%extCharges, self%blurWidths, self%rVec, self%gVec, self%alpha, self%volume, &
        & self%getDenseAO, self%getDenseAtom, self%getAtomIndex, self%orderRmatL, &
        & self%Lpaired, self%SAstates, self%tNAC, self%tRangeSep, self%tExtChrg, &
        & self%tPeriodic, self%tBlur, self%SAgrad, self%SIgrad, self%SSRgrad)

  end subroutine getRTGradient_


  !> Get the information related to non-adiabatic coupling
  subroutine getReksNACinfo_(tWriteTagged, autotestTag, taggedWriter, self)

    !> print tag information
    logical, intent(in) :: tWriteTagged

    !> File name for regression data
    character(*), intent(in) :: autotestTag

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    integer :: fdTagged

    call weightGradient(self%gradL, self%weight, self%avgGrad)
    call getOtherSAgrad(self%avgGrad, self%tSSR22, self%tSSR44, self%SAgrad)
    call SaToSsrGradient(self%SAgrad, self%SIgrad, self%eigvecsSSR, self%SSRgrad)

    call getReksNAC(self%SAgrad, self%SIgrad, self%SSRgrad, self%eigvecsSSR, &
        & self%energy, self%nacG, self%nacH)

    if (tWriteTagged) then
      open(newUnit=fdTagged, file=autotestTag, position="append")
      ! nonadiabatic coupling vector has a phase, just check the value not sign
      call taggedWriter%write(fdTagged, tagLabels%nacH, abs(self%nacH))
      close(fdTagged)
    end if

  end subroutine getReksNACinfo_


  !> Set the final gradient for REKS
  subroutine setReksGradients_(derivs, self)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    !> derivatives of energy wrt to atomic positions
    real(dp), intent(out) :: derivs(:,:)

    ! when Efunc = 1, derivs is already defined in previous routine
    if (self%tNAC) then
      derivs(:,:) = self%SSRgrad(:,:,self%rstate)
    else
      derivs(:,:) = self%SSRgrad(:,:,1)
    end if

  end subroutine setReksGradients_


  !> Creates (relaxed) density and mulliken population for target state
  subroutine getMullikenPopFromRelaxedDensity_(env, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, orb, over, rhoPrim, qOutput, self)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> sparse density matrix
    real(dp), intent(inout) :: rhoPrim(:,:)

    !> Output electrons
    real(dp), intent(inout) :: qOutput(:,:,:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    rhoPrim(:,:) = 0.0_dp
    call env%globalTimer%startTimer(globalTimers%denseToSparse)
    call packHS(rhoPrim(:,1), 0.5_dp*(self%relRhoSqr + transpose(self%relRhoSqr)), &
        & neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb, &
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)
    call env%globalTimer%stopTimer(globalTimers%denseToSparse)

    qOutput(:,:,:) = 0.0_dp
    call mulliken(qOutput(:,:,1), over, rhoPrim(:,1), orb, &
        & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)

  end subroutine getMullikenPopFromRelaxedDensity_


  !> Calculates dipole moment.
  !> This is same routine with getDipoleMoment in main.F90
  subroutine getDipoleMoment_(qOutput, q0, coord, dipoleMoment, iAtInCentralRegion)

    !> electrons in orbitals
    real(dp), intent(in) :: qOutput(:,:,:)

    !> reference atomic charges
    real(dp), intent(in) :: q0(:,:,:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> resulting dipole moment
    real(dp), intent(inout) :: dipoleMoment(:)

    !> atoms in the central cell (or device region if transport)
    integer, intent(in) :: iAtInCentralRegion(:)

    integer :: nAtom, ii, iAtom

    nAtom = size(qOutput, dim=2)
    dipoleMoment(:) = 0.0_dp
    do ii = 1, size(iAtInCentralRegion)
      iAtom = iAtInCentralRegion(ii)
      dipoleMoment(:) = dipoleMoment(:)&
          & + sum(q0(:, iAtom, 1) - qOutput(:, iAtom, 1)) * coord(:,iAtom)
    end do

  end subroutine getDipoleMoment_


end module dftbp_reksinterface
