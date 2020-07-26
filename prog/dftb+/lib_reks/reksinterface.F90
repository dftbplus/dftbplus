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
  use dftbp_elecsolvers
  use dftbp_environment
  use dftbp_globalenv
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
  use dftbp_reksen
  use dftbp_reksgrad
  use dftbp_reksio
  use dftbp_reksproperty
  use dftbp_reksvar

  implicit none

  private

  public :: getStateInteraction, getReksEnProperties
  public :: getReksGradients, getReksGradProperties
  public :: getReksStress

  contains

  !> Calculate SSR state from SA-REKS states and state-interaction terms
  subroutine getStateInteraction(env, denseDesc, neighbourList, nNeighbourSK,&
      & iSparseStart, img2CentCell, coord, iAtInCentralRegion, eigenvecs,&
      & electronicSolver, eigen, qOutput, q0, tDipole, dipoleMoment, this)

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

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Atoms over which to sum the total energies
    integer, intent(in) :: iAtInCentralRegion(:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> eigenvalues
    real(dp), intent(inout) :: eigen(:,:,:)

    !> Output electrons
    real(dp), intent(in) :: qOutput(:,:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> calculate an electric dipole?
    logical, intent(in) :: tDipole

    !> resulting dipole moment
    real(dp), allocatable, intent(inout) :: dipoleMoment(:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    call adjustEigenval(this, eigen)

    if (this%Efunction > 1) then
      call solveSecularEqn(env, denseDesc, neighbourList, nNeighbourSK, &
          & iSparseStart, img2CentCell, electronicSolver, eigenvecs, this)
    else
      ! Get the dipole moment for single-state REKS case
      ! In this case dipole moment can be calculated w/o gradient result
      ! tDipole = (total charge = 0.0) * (non-periodic system) * (mulliken)
      if (tDipole) then
        call getDipoleMoment_(qOutput, q0, coord, dipoleMoment, iAtInCentralRegion)
      end if
    end if


  end subroutine getStateInteraction


  !> get the energy-related properties; unrelaxed density matrix,
  !> dipole integral, transition dipole, oscillator strength
  subroutine getReksEnProperties(eigenvecs, coord0, this)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:,:)

    !> central cell coordinates of atoms
    real(dp), intent(in) :: coord0(:,:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    real(dp), allocatable :: dipoleInt(:,:,:)

    integer :: ist, nstHalf, nOrb

    nOrb = size(eigenvecs,dim=1)
    nstHalf = this%nstates * (this%nstates - 1) / 2

    allocate(dipoleInt(nOrb,nOrb,3))

    ! Get the unrelaxed density matrix for SA-REKS or SSR state
    ! The matrix that used in this calculation is not relaxed density
    ! matrix, so this unrelaxed FONs are not equal to relaxed FONS,
    ! but this value is easy to calculate without the information of
    ! gradient. Therefore, we can easily guess the behavior of the states.
    if (this%nstates > 1) then

      call getUnrelaxedDensMatAndTdp(eigenvecs(:,:,1), this%overSqr, this%rhoSqrL, &
          & this%FONs, this%eigvecsSSR, this%Lpaired, this%Nc, this%Na, &
          & this%rstate, this%Lstate, this%reksAlg, this%tSSR, this%tTDP, &
          & this%unrelRhoSqr, this%unrelTdm)

      if (this%tTDP) then
        call getDipoleIntegral(coord0, this%overSqr, this%getAtomIndex, dipoleInt)
        ! Get the transition dipole moment between states
        ! For (SI-)SA-REKS dipole moment requires gradient info.
        ! But TDP use only zero-th part without gradient info.
        do ist = 1, nstHalf
          call getDipoleMomentMatrix(this%unrelTdm(:,:,ist), dipoleInt, this%tdp(:,ist))
        end do
        call writeReksTDP(this%tdp)
        call getReksOsc(this%tdp, this%energy)
      end if

    end if

  end subroutine getReksEnProperties


  !> Calculate SI-SA-REKS state gradient by solving CP-REKS equations
  subroutine getReksGradients(env, denseDesc, sccCalc, rangeSep, dispersion, &
      & neighbourList, nNeighbourSK, nNeighbourRep, iSparseStart, img2CentCell, &
      & orb, nonSccDeriv, skHamCont, skOverCont, pRepCont, coord, coord0, &
      & species, q0, eigenvecs, chrgForces, over, spinW, derivs, tWriteTagged, &
      & autotestTag, taggedWriter, this)

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
    type(TReksCalc), intent(inout) :: this

    real(dp), allocatable :: Qmat(:,:)
    integer :: ist, ia, ib, nstHalf, fac

    nstHalf = this%nstates * (this%nstates-1) / 2

    allocate(Qmat(orb%nOrb,orb%nOrb))

    ! get the periodic information
    if (this%tPeriodic) then
      call sccCalc%coulombCont%getPeriodicInfo(this%rVec, this%gVec, this%alpha, this%volume)
    end if

    call getHellmannFeynmanGradientL_(env, denseDesc, sccCalc, neighbourList, &
        & nNeighbourSK, nNeighbourRep, iSparseStart, img2CentCell, orb, &
        & nonSccDeriv, skHamCont, skOverCont, pRepCont, coord, species, q0, &
        & dispersion, rangeSep, chrgForces, eigenvecs, derivs, this)

    if (this%Efunction == 1) then
      call weightGradient(this%gradL, this%weight, derivs)
    else

      ! get REKS parameters used in CP-REKS and gradient equations
      call getReksParameters_(env, denseDesc, sccCalc, rangeSep, &
          & neighbourList, nNeighbourSK, iSparseStart, img2CentCell, &
          & eigenvecs, coord, species, over, spinW, this)

      call buildStateVectors_(env, denseDesc, neighbourList, nNeighbourSK, &
          & iSparseStart, img2CentCell, eigenvecs, orb, over, this)

      if (this%tNAC) then

        ! SA-REKS calculations
        do ist = 1, this%nstates
          if (ist /= this%SAstates) then

            write(stdOut,"(A)")
            write(stdOut,"(A)") repeat("-", 82)
            write(stdOut,'(1x,a,1x,I2,1x,a)') &
                & 'Solving CP-REKS equation for', ist, 'state vector...'

            ! solve CP-REKS equation for SA-REKS state
            ! save information about ZT, RmatL
            call solveCpReks_(env, denseDesc, neighbourList, nNeighbourSK, &
                & iSparseStart, img2CentCell, eigenvecs, over, orb, this, &
                & this%XT(:,ist), this%ZT(:,ist), this%RmatL(:,:,:,ist), &
                & this%ZmatL, this%Q1mat, this%Q2mat, optionQMMM=.false.)
            Qmat(:,:) = this%Q1mat + this%Q2mat
            ! compute SA-REKS shift
            call SSRshift(eigenvecs, this%gradL, Qmat, this%Sderiv, &
                & this%ZT(:,ist), this%SAweight, this%weightL(ist,:), &
                & this%omega, this%weightIL, this%G1, &
                & denseDesc%iAtomStart, orb%mOrb, this%SAgrad(:,:,ist), 1)

          end if
        end do

        ! state-interaction calculations
        do ist = 1, nstHalf

          call getTwoIndices(this%nstates, ist, ia, ib, 1)
          write(stdOut,"(A)")
          write(stdOut,"(A)") repeat("-", 82)
          write(stdOut,'(1x,a,1x,I2,1x,a,1x,I2,1x,a)') &
              & 'Solving CP-REKS equation for SI between', ia, 'and', ib, 'state vectors...'

          ! solve CP-REKS equation for state-interaction term
          ! save information about ZTdel, tmpRL
          call solveCpReks_(env, denseDesc, neighbourList, nNeighbourSK, &
              & iSparseStart, img2CentCell, eigenvecs, over, orb, this, &
              & this%XTdel(:,ist), this%ZTdel(:,ist), this%tmpRL(:,:,:,ist), &
              & this%ZmatL, this%Q1mat, this%Q2mat, optionQMMM=.false.)
          call SIshift(eigenvecs, this%gradL, this%Q1del(:,:,ist), &
              & this%Q2del(:,:,ist), this%Q1mat, this%Q2mat, Qmat, &
              & this%Sderiv, this%ZTdel(:,ist), this%SAweight, this%omega, &
              & this%weightIL, this%Rab, this%G1, denseDesc%iAtomStart, &
              & orb%mOrb, ist, this%nstates, this%SIgrad(:,:,ist))

        end do

      else

        if (this%tSSR) then

          ! Convert necessary variables for SSR state from SA and SI terms
          call SaToSsrXT(this%XTdel, this%eigvecsSSR, this%rstate, this%XT)
          call SaToSsrWeight(this%Rab, this%weightIL, this%G1, &
              & this%eigvecsSSR, this%rstate, this%weightL)

          write(stdOut,"(A)")
          write(stdOut,"(A)") repeat("-", 82)
          write(stdOut,'(1x,a,1x,I2,1x,a)') &
              & 'Solving CP-REKS equation for', this%rstate, 'state vector...'

          ! solve CP-REKS equation for SSR state
          ! save information about ZT, RmatL, ZmatL, Q1mat, Q2mat
          call solveCpReks_(env, denseDesc, neighbourList, nNeighbourSK, &
              & iSparseStart, img2CentCell, eigenvecs, over, orb, this, &
              & this%XT(:,this%rstate), this%ZT(:,1), this%RmatL(:,:,:,1), &
              & this%ZmatL, this%Q1mat, this%Q2mat, optionQMMM=.false.)

          ! add remaining SI component to SSR state
          call addSItoRQ(eigenvecs, this%RdelL, this%Q1del, this%Q2del, &
              & this%eigvecsSSR, this%rstate, this%RmatL(:,:,:,1), &
              & this%Q1mat, this%Q2mat, Qmat)
          fac = 2

        else

          write(stdOut,"(A)")
          write(stdOut,"(A)") repeat("-", 82)
          if (this%Lstate == 0) then
            write(stdOut,'(1x,a,1x,I2,1x,a)') &
                & 'Solving CP-REKS equation for', this%rstate, 'state vector...'
          else
            write(stdOut,'(1x,a,1x,I2,1x,a)') &
                & 'Solving CP-REKS equation for', this%Lstate, 'microstate vector...'
          end if

          ! solve CP-REKS equation for SA-REKS or L state
          ! save information about ZT, RmatL, ZmatL, Q1mat, Q2mat
          call solveCpReks_(env, denseDesc, neighbourList, nNeighbourSK, &
              & iSparseStart, img2CentCell, eigenvecs, over, orb, this, &
              & this%XT(:,1), this%ZT(:,1), this%RmatL(:,:,:,1), &
              & this%ZmatL, this%Q1mat, this%Q2mat, optionQMMM=.false.)
          Qmat(:,:) = this%Q1mat + this%Q2mat
          fac = 1

        end if

        if (this%Lstate == 0) then
          ! compute SSR or SA-REKS shift
          call SSRshift(eigenvecs, this%gradL, Qmat, this%Sderiv, &
              & this%ZT(:,1), this%SAweight, this%weightL(this%rstate,:), &
              & this%omega, this%weightIL, this%G1, &
              & denseDesc%iAtomStart, orb%mOrb, this%SSRgrad(:,:,1), fac)
        else
          ! compute L shift
          call Lshift(eigenvecs, this%gradL, Qmat, this%Sderiv, &
              & this%ZT(:,1), this%SAweight, this%omega, this%weightIL, this%G1, &
              & denseDesc%iAtomStart, orb%mOrb, this%Lstate, this%SSRgrad(:,:,1))
        end if

      end if

      ! compute R*T shift & final SSR gradient (grad, nac)
      call getRTGradient_(env, sccCalc, denseDesc, neighbourList, &
          & nNeighbourSK, iSparseStart, img2CentCell, orb, q0, coord0, this)
      if (this%tNAC) then
        call getReksNACinfo_(tWriteTagged, autotestTag, taggedWriter, this)
      end if

      ! set final gradient to derivs in SA-REKS or SSR case
      call setReksGradients_(derivs, this)

    end if

    if (this%Plevel >= 1) then
      call printReksGradInfo(this, derivs)
    end if

  end subroutine getReksGradients


  !> Get the gradient-related properties such as relaxed density, charge, etc
  subroutine getReksGradProperties(env, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, eigenvecs, orb, &
      & iAtInCentralRegion, coord, coord0, over, rhoPrim, qOutput, &
      & q0, tDipole, dipoleMoment, chrgForces, this)

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
    type(TReksCalc), intent(inout) :: this

    ! calculate the relaxed density matrix, dipole moment,
    !           the forces of external charges.
    if (this%tRD) then

      if (this%Efunction > 1) then

        ! currently, relaxed density matrix is possible
        !            relaxed charge is possible
        !            relaxed dipole moment is possible
        !            force of external point charges for target state is possible
        if (this%tNAC) then
          ! solve CP-REKS equation for remaining SA-REKS state in nac case
          ! save information about Z_T
          call solveCpReks_(env, denseDesc, neighbourList, nNeighbourSK, &
              & iSparseStart, img2CentCell, eigenvecs, over, orb, this, &
              & this%XT(:,this%SAstates), this%ZT(:,this%SAstates), &
              & this%RmatL(:,:,:,this%SAstates), this%ZmatL, &
              & this%Q1mat, this%Q2mat, optionQMMM=.true.)
          ! now, ZT has information about target SSR state
          call SaToSsrXT(this%ZTdel, this%eigvecsSSR, this%rstate, this%ZT)
        end if

        if (this%Lstate == 0) then
          ! get the relaxed density matrix for target SSR or SA-REKS state
          call getRelaxedDensMat(eigenvecs(:,:,1), this%overSqr, this%unrelRhoSqr, &
              & this%ZT, this%omega, this%FONs, this%eigvecsSSR, this%SAweight, &
              & this%Rab, this%G1, this%Nc, this%Na, this%rstate, this%reksAlg, &
              & this%tSSR, this%tNAC, this%relRhoSqr)
        else
          ! get the relaxed density matrix for L-th microstate
          call getRelaxedDensMatL(eigenvecs(:,:,1), this%rhoSqrL, this%overSqr, &
              & this%weight, this%SAweight, this%unrelRhoSqr, this%RmatL, &
              & this%ZT, this%omega, this%weightIL, this%G1, this%orderRmatL, &
              & this%Lpaired, this%Nc, this%Na, this%Lstate, this%reksAlg, &
              & this%relRhoSqr)
        end if

        call getMullikenPopFromRelaxedDensity_(env, denseDesc, neighbourList, &
            & nNeighbourSK, iSparseStart, img2CentCell, orb, over, rhoPrim, qOutput, this)

      end if

      call writeReksRelaxedCharge(qOutput, q0, this%rstate, this%Lstate)

      if (this%Efunction > 1) then
        ! Get the dipole moment for target SSR state
        ! Relaxed dipole moment requires gradient infomation
        ! tDipole = (total charge = 0.0) * (non-periodic system) * (mulliken)
        if (tDipole) then
          call getDipoleMoment_(qOutput, q0, coord, dipoleMoment, iAtInCentralRegion)
        end if
      end if

      if (this%tExtChrg) then

        call getExtChrgGradients(env, coord0, this%extCharges(1:3,:), &
            & qOutput, q0, this%extCharges(4,:), this%blurWidths, this%rVec, this%gVec, &
            & this%alpha, this%volume, this%tPeriodic, this%tBlur, chrgForces)

      end if

    end if

  end subroutine getReksGradProperties


  !> Calculates stress tensor and lattice derivatives.
  subroutine getReksStress(env, denseDesc, sccCalc, nonSccDeriv, &
      & skHamCont, skOverCont, pRepCont, neighbourList, nNeighbourSk, &
      & nNeighbourRep, species, img2CentCell, iSparseStart, orb, &
      & dispersion, coord, q0, invLatVec, cellVol, totalStress, &
      & totalLatDeriv, intPressure, this)

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

    !> non-SCC hamiltonian information
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
    type(TReksCalc), intent(inout) :: this

    real(dp), allocatable :: tmpRhoSp(:,:)
    real(dp), allocatable :: tmpStress(:,:)

    real(dp) :: fac
    integer :: sparseSize, nSpin, tmpLu, tmpLd, iL

    sparseSize = size(this%getDenseAO,dim=1)
    nSpin = size(this%fillingL,dim=2)

    if (this%Efunction == 1) then

      allocate(tmpRhoSp(sparseSize,nSpin))
      allocate(tmpStress(3,3))

      do iL = 1, this%Lmax

        if (iL <= this%Lpaired) then
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

        ! this%qOutputL has (qm) component
        call sccCalc%updateCharges(env, this%qOutputL(:,:,:,iL), q0, orb, species)
        call sccCalc%updateShifts(env, orb, species, &
            & neighbourList%iNeighbour, img2CentCell)

        call env%globalTimer%startTimer(globalTimers%denseToSparse)
        ! this%rhoSqrL has (my_qm) component
        tmpRhoSp(:,:) = 0.0_dp
        call packHS(tmpRhoSp(:,1), this%rhoSqrL(:,:,1,tmpLu), &
            & neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb, &
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        if (iL > this%Lpaired) then
          call packHS(tmpRhoSp(:,2), fac*this%rhoSqrL(:,:,1,tmpLd), &
              & neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb, &
              & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
        call env%globalTimer%stopTimer(globalTimers%denseToSparse)

        tmpStress(:,:) = 0.0_dp
        call getBlockStress(env, tmpStress, nonSccDeriv, tmpRhoSp, &
            & this%edmSpL(:,iL), skHamCont, skOverCont, coord, &
            & species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell, &
            & iSparseStart, orb, this%intBlockL(:,:,:,:,iL), cellVol)
        call sccCalc%addStressDc(tmpStress, env, species, &
            & neighbourList%iNeighbour, img2CentCell)

!        if (allocated(thirdOrd)) then
!          call thirdOrd%updateCharges(pSpecies0, neighbourList, &
!              & this%qOutput_L(:,:,:,iL), q0, img2CentCell, orb)
!          call thirdOrd%addStressDc(neighbourList, species, coord, &
!              & img2CentCell, cellVol, tmpStress)
!        end if

        this%elecStressL(:,:,iL) = tmpStress

      end do

      call weightGradient(this%elecStressL, this%weight, totalStress)

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
      & dispersion, rangeSep, chrgForces, eigenvecs, derivs, this)

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
    type(TReksCalc), intent(inout) :: this

    real(dp), allocatable :: repDerivs(:,:)
    real(dp), allocatable :: dispDerivs(:,:)
    real(dp), allocatable :: lcDerivs(:,:,:)
!    integer, pointer :: pSpecies0(:)

    integer :: nAtom, iL

    nAtom = size(this%gradL,dim=2)
!    pSpecies0 => species(1:nat)

    allocate(repDerivs(3,nAtom))
    allocate(dispDerivs(3,nAtom))
    if (this%isRangeSep) then
      allocate(lcDerivs(3,nAtom,this%Lmax))
    end if

    ! hamSpL has (my_ud) component
    ! hamSqrL has (my_ud) component
    call env%globalTimer%startTimer(globalTimers%energyDensityMatrix)
    call getEnergyWeightedDensityL(env, denseDesc, neighbourList, &
        & nNeighbourSK, iSparseStart, img2CentCell, orb, this%hamSqrL, &
        & this%hamSpL, this%fillingL, eigenvecs(:,:,1), this%Lpaired, &
        & this%Efunction, this%isRangeSep, this%edmSpL)
    call env%globalTimer%stopTimer(globalTimers%energyDensityMatrix)

    ! rhoSpL has (my_qm) component
    ! intBlockL has (qm) component
    call derivative_blockL(env, derivs, nonSccDeriv, this%rhoSqrL, &
        & this%edmSpL, skHamCont, skOverCont, coord, species, &
        & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, &
        & iSparseStart, denseDesc%iAtomStart, orb, this%intBlockL, &
        & this%Lpaired, this%gradL, this%Hderiv, this%Sderiv)

    do iL = 1, this%Lmax

      derivs(:,:) = 0.0_dp

      ! qOutput_L has (qm) component
      call sccCalc%updateCharges(env, this%qOutputL(:,:,:,iL), q0, orb, species)
      call sccCalc%updateShifts(env, orb, species, &
          & neighbourList%iNeighbour, img2CentCell)
      if (this%tExtChrg) then
        chrgForces(:,:) = 0.0_dp
        call sccCalc%addForceDc(env, derivs, species, &
            & neighbourList%iNeighbour, img2CentCell, chrgForces)
      else
        call sccCalc%addForceDc(env, derivs, species, &
            & neighbourList%iNeighbour, img2CentCell)
      endif

!      if (allocated(thirdOrd)) then
!        call thirdOrd%updateCharges(pSpecies0, neighbourList, &
!            & this%qOutputL(:,:,:,iL), q0, img2CentCell, orb)
!        call thirdOrd%addGradientDc(neighbourList, species, &
!            & coord, img2CentCell, derivs)
!      end if

      if (this%isRangeSep) then
        ! deltaRhoSqrL has (my_ud) component
        lcDerivs(:,:,iL) = 0.0_dp
        call rangeSep%addLRGradients(lcDerivs(:,:,iL), nonSccDeriv, this%deltaRhoSqrL(:,:,:,iL),&
            & skOverCont, coord, species, orb, denseDesc%iAtomStart, this%overSqr,&
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

      this%gradL(:,:,iL) = this%gradL(:,:,iL) + derivs

    end do

    if(this%isRangeSep) then
      do iL = 1, this%Lmax
        if (iL <= this%Lpaired) then
          derivs(:,:) = lcDerivs(:,:,iL) + lcDerivs(:,:,iL)
        else
          if (mod(iL,2) == 1) then
            derivs(:,:) = lcDerivs(:,:,iL) + lcDerivs(:,:,iL+1)
          else
            derivs(:,:) = lcDerivs(:,:,iL) + lcDerivs(:,:,iL-1)
          end if
        end if
        this%gradL(:,:,iL) = this%gradL(:,:,iL) + derivs
      end do
    end if

  end subroutine getHellmannFeynmanGradientL_


  !> Set several REKS variables used in CP-REKS equations
  subroutine getReksParameters_(env, denseDesc, sccCalc, rangeSep, &
      & neighbourList, nNeighbourSK, iSparseStart, img2CentCell, &
      & eigenvecs, coord, species, over, spinW, this)

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
    type(TReksCalc), intent(inout) :: this

    ! get gamma, spinW, gamma deriv, LR-gamma deriv, on-site constants
    ! LRgamma is already defined in 1st scc loop
    call getSccSpinLrPars(env, sccCalc, rangeSep, coord, species, &
        & neighbourList%iNeighbour, img2CentCell, denseDesc%iAtomStart, &
        & spinW, this%getAtomIndex, this%isRangeSep, this%GammaAO, &
        & this%GammaDeriv, this%SpinAO, this%LrGammaAO, this%LrGammaDeriv)

    ! get Hxc kernel -> (\mu,\nu|f_{Hxc}|\tau,\gam)
    call getHxcKernel(this%getDenseAO, over, this%overSqr, this%GammaAO, this%SpinAO,&
        & this%LrGammaAO, this%Glevel, this%tSaveMem, this%isRangeSep, this%HxcSpS, &
        & this%HxcSpD, this%HxcHalfS, this%HxcHalfD, this%HxcSqrS, this%HxcSqrD)

    ! get G1, weightIL, Omega, Rab values
    call getG1ILOmegaRab(env, denseDesc, neighbourList, nNeighbourSK, &
        & iSparseStart, img2CentCell, eigenvecs, this%hamSqrL, this%hamSpL, &
        & this%fockFa, this%fillingL, this%FONs, this%SAweight, this%enLtot, &
        & this%hess, this%Nc, this%Na, this%reksAlg, this%tSSR, &
        & this%isRangeSep, this%G1, this%weightIL, this%omega, this%Rab)

    ! get A1e or Aall values based on GradOpt
    call getSuperAMatrix(eigenvecs, this%HxcSqrS, this%HxcSqrD, this%fockFc, &
        & this%fockFa, this%omega, this%fillingL, this%weight, this%SAweight, &
        & this%FONs, this%G1, this%Lpaired, this%Nc, this%Na, this%Glevel, &
        & this%reksAlg, this%tSaveMem, this%A1e, this%A1ePre, this%Aall)

  end subroutine getReksParameters_


  !> Calculate SA-REKS state vectors and state-interaction vectors
  subroutine buildStateVectors_(env, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, eigenvecs, orb, over, this)

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
    type(TReksCalc), intent(inout) :: this

    integer :: ia, ib, ist, nstHalf

    nstHalf = this%nstates * (this%nstates - 1) / 2

    if (this%Lstate == 0) then

      ! build X^T value for state X = PPS, OSS, etc
      call buildSaReksVectors(env, denseDesc, neighbourList, nNeighbourSK, &
          & iSparseStart, img2CentCell, eigenvecs, this%hamSqrL, this%hamSpL, &
          & this%fillingL, this%weightL, this%Nc, this%Na, this%rstate, &
          & this%reksAlg, this%tSSR, this%isRangeSep, this%XT)

      if (this%tSSR) then

        ! get R^delta values between states X's
        call getRdel(eigenvecs, this%fillingL, this%FONs, this%Nc, &
            & this%nstates, this%reksAlg, this%RdelL)

        do ist = 1, nstHalf

          call getTwoIndices(this%nstates, ist, ia, ib, 1)

          ! get Z^delta values from R^delta values
          call getZmat(env, denseDesc, neighbourList, nNeighbourSK, &
              & iSparseStart, img2CentCell, orb, this%RdelL(:,:,:,ist), &
              & this%HxcSqrS, this%HxcSqrD, this%HxcHalfS, this%HxcHalfD, &
              & this%HxcSpS, this%HxcSpD, this%overSqr, over, this%GammaAO, &
              & this%SpinAO, this%LrGammaAO, this%orderRmatL, this%getDenseAO, &
              & this%Lpaired, this%Glevel, this%tSaveMem, this%isRangeSep, this%ZdelL)

          ! build XTdel with Z^delta values
          call buildInteractionVectors(eigenvecs, this%ZdelL, this%fockFc, &
              & this%fockFa, this%FONs, this%fillingL, this%weight, &
              & this%SAweight, this%omega, this%Rab, this%G1, this%Nc, this%Na, &
              & ia, ib, this%reksAlg, this%XTdel(:,ist))

          ! get Q2^delta values from Z^delta values : MO index
          call getQ2mat(eigenvecs, this%fillingL, this%weight, &
              & this%ZdelL, this%Q2del(:,:,ist))

        end do

        ! get Q1^delta values : AO index
        call getQ1del(eigenvecs, this%fockFc, this%fockFa, this%FONs, this%SAweight, &
            & this%Nc, this%nstates, this%reksAlg, this%Q1del)

      end if

    else

      ! build X^T_L value
      call buildLstateVector(env, denseDesc, neighbourList, nNeighbourSK, &
          & iSparseStart, img2CentCell, eigenvecs, this%hamSqrL, this%hamSpL, &
          & this%fillingL, this%Nc, this%Na, this%Lstate, this%Lpaired, &
          & this%reksAlg, this%isRangeSep, this%XT(:,1))

    end if

  end subroutine buildStateVectors_


  !> Solve CP-REKS equations by using CG based algorithms or direct matrix multiplication
  subroutine solveCpReks_(env, denseDesc, neighbourList, nNeighbourSK, &
      & iSparseStart, img2CentCell, eigenvecs, over, orb, this, XT, ZT, &
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
    type(TReksCalc), intent(inout) :: this

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

    if (this%Glevel == 1 .or. this%Glevel == 2) then

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
          & img2CentCell, orb, XT, this%A1e, this%A1ePre, this%HxcSqrS, &
          & this%HxcSqrD, this%HxcHalfS, this%HxcHalfD, this%HxcSpS, this%HxcSpD, &
          & this%fockFc, this%fockFa, this%omega, this%SAweight, this%FONs, &
          & this%G1, this%GammaAO, this%SpinAO, this%LrGammaAO, this%overSqr, &
          & over, eigenvecs, this%fillingL, this%weight, this%Glimit, this%orderRmatL, &
          & this%getDenseAO, this%Lpaired, this%Nc, this%Na, this%CGmaxIter, this%Glevel, &
          & this%reksAlg, this%tSaveMem, this%isRangeSep, ZT, RmatL, ZmatL, Q2mat)

    else if (this%Glevel == 3) then

      ! solve A * Z = X using direct matrix inversion
      call solveZT(this%Aall, XT, ZT)

      ! get direct R, Z, Q2 matrices
      if (.not. optionQMMM) then
        call getRmat(eigenvecs, ZT, this%fillingL, this%Nc, this%Na, &
            & this%reksAlg, RmatL)
        call getZmat(env, denseDesc, neighbourList, nNeighbourSK, &
            & iSparseStart, img2CentCell, orb, RmatL, &
            & this%HxcSqrS, this%HxcSqrD, this%HxcHalfS, this%HxcHalfD, &
            & this%HxcSpS, this%HxcSpD, this%overSqr, over, this%GammaAO, &
            & this%SpinAO, this%LrGammaAO, this%orderRmatL, this%getDenseAO, &
            & this%Lpaired, this%Glevel, this%tSaveMem, this%isRangeSep, ZmatL)
        call getQ2mat(eigenvecs, this%fillingL, this%weight, ZmatL, Q2mat)
        write(stdOut,"(A)") repeat("-", 82)
      end if

    end if

    ! get Q1 matrix from converged ZT vector : MO index
    if (.not. optionQMMM) then
      call getQ1mat(ZT, this%fockFc, this%fockFa, this%SAweight, &
          & this%FONs, this%Nc, this%Na, this%reksAlg, Q1mat)
    else
      write(stdOut,"(A)")
    end if

  end subroutine solveCpReks_


  !> Calculate R*T contribution of gradient (2nd term)
  subroutine getRTGradient_(env, sccCalc, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, orb, q0, coord0, this)

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
    type(TReksCalc), intent(inout) :: this

    call RTshift(env, sccCalc, denseDesc, neighbourList, nNeighbourSK, &
        & iSparseStart, img2CentCell, orb, coord0, this%Hderiv, this%Sderiv, &
        & this%rhoSqrL, this%overSqr, this%deltaRhoSqrL, this%qOutputL, &
        & q0, this%GammaAO, this%GammaDeriv, this%SpinAO, this%LrGammaAO, &
        & this%LrGammaDeriv, this%RmatL, this%RdelL, this%tmpRL, this%weight, &
        & this%extCharges, this%blurWidths, this%rVec, this%gVec, this%alpha, this%volume, &
        & this%getDenseAO, this%getDenseAtom, this%getAtomIndex, this%orderRmatL, &
        & this%Lpaired, this%SAstates, this%tNAC, this%isRangeSep, this%tExtChrg, &
        & this%tPeriodic, this%tBlur, this%SAgrad, this%SIgrad, this%SSRgrad)

  end subroutine getRTGradient_


  !> Get the information related to non-adiabatic coupling
  subroutine getReksNACinfo_(tWriteTagged, autotestTag, taggedWriter, this)

    !> print tag information
    logical, intent(in) :: tWriteTagged

    !> File name for regression data
    character(*), intent(in) :: autotestTag

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    integer :: fdTagged

    call weightGradient(this%gradL, this%weight, this%avgGrad)
    call getOtherSAgrad(this%avgGrad, this%reksAlg, this%SAgrad)
    call SaToSsrGradient(this%SAgrad, this%SIgrad, this%eigvecsSSR, this%SSRgrad)

    call getReksNAC(this%SAgrad, this%SIgrad, this%SSRgrad, this%eigvecsSSR, &
        & this%energy, this%nacG, this%nacH)

    if (tWriteTagged) then
      open(newUnit=fdTagged, file=autotestTag, position="append")
      ! nonadiabatic coupling vector has a phase, just check the value not sign
      call taggedWriter%write(fdTagged, tagLabels%nacH, abs(this%nacH))
      close(fdTagged)
    end if

  end subroutine getReksNACinfo_


  !> Set the final gradient for REKS
  subroutine setReksGradients_(derivs, this)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> derivatives of energy wrt to atomic positions
    real(dp), intent(out) :: derivs(:,:)

    ! when Efunc = 1, derivs is already defined in previous routine
    if (this%tNAC) then
      derivs(:,:) = this%SSRgrad(:,:,this%rstate)
    else
      derivs(:,:) = this%SSRgrad(:,:,1)
    end if

  end subroutine setReksGradients_


  !> Creates (relaxed) density and mulliken population for target state
  subroutine getMullikenPopFromRelaxedDensity_(env, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, orb, over, rhoPrim, qOutput, this)

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
    type(TReksCalc), intent(inout) :: this

    rhoPrim(:,:) = 0.0_dp
    call env%globalTimer%startTimer(globalTimers%denseToSparse)
    call packHS(rhoPrim(:,1), 0.5_dp*(this%relRhoSqr + transpose(this%relRhoSqr)), &
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
