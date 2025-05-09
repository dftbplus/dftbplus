!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

!> Contains the calculator for the short-range part of the Gamma-electrostatics.
module dftbp_dftb_shortgamma
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_schedule, only : getIndicesWithWorkload
  use dftbp_dftb_h5correction, only : TH5Correction, TH5Correction_init, TH5CorrectionInput
  use dftbp_dftb_periodic, only : getNrOfNeighbours, TNeighbourList
  use dftbp_dftb_shortgammafuncs, only : expGamma, expGammaCutOff, expGammaDamped,&
      & expGammaDampedPrime, expGammaPrime
  use dftbp_dftb_uniquehubbard, only : TUniqueHubbard
  use dftbp_type_commontypes, only : TOrbitals
#:if WITH_SCALAPACK
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip
#:endif
  implicit none

  private
  public :: TShortGammaInput, TShortGamma, TShortGamma_init, TShortGammaDamp


  !> Contains damping details.
  type :: TShortGammaDamp

    !> Flag for each species, whether damping should be applied to it. Shape: [nSpecies]
    logical, allocatable :: isDamped(:)

    !> Damping exponent
    real(dp) :: exponent

  end type TShortGammaDamp


  !> Input parameters for short gamma calculation.
  type :: TShortGammaInput

    !> Unique Hubbard U parameters
    type(TUniqueHubbard), allocatable :: hubbU

    !> Optional H5 correction
    type(TH5CorrectionInput), allocatable :: h5CorrectionInp

    !> Optional damping
    type(TShortGammaDamp), allocatable :: damping

  end type TShortGammaInput


  !> Contains the internal state of the short-gamma interaction calculator.
  type :: TShortGamma
    private

    !> Nr. of species
    integer :: nSpecies_

    !> Maximal shells per species
    integer :: mShell_

    !> Nr. of atoms
    integer :: nAtom_

    !> Contracted Hubbard U values
    type(TUniqueHubbard), allocatable :: hubbU_

    !> Cutoff for short range interaction. Shape [mHubbU, mHubbU, nSpecies, nSpecies]
    real(dp), allocatable :: shortCutoffs_(:,:,:,:)

    !> Number of neighbours. Shape: [mHubbU, mHubbU, nSpecies, nAtom]
    integer, allocatable :: nNeigh_(:,:,:,:)

    !> Net charges per shell
    real(dp), allocatable :: deltaQShell_(:,:)

    !> Net charges summed up over equivalent Hubbard U values
    real(dp), allocatable :: deltaQUniqU_(:,:)

    !> H5 correction calculator
    type(TH5Correction), allocatable :: h5Correction_

    !> Gamma damping calculator
    type(TShortGammaDamp), allocatable :: damping_

    !> Shell resolved shift
    real(dp), allocatable :: shiftShell_(:,:)

    !> Cache storing the short gamma interaction
    real(dp), allocatable :: shortGamma_(:,:,:,:)

  contains

    procedure :: getCutoff
    procedure :: updateCoords
    procedure :: updateCharges
    procedure :: updateShifts
    procedure :: getShiftPerShell
    procedure :: addAtomicMatrix
    procedure :: addAtomicMatrixCustomU
    procedure :: addGradientsDc
    procedure :: addStressDc
    procedure :: addGradientsDcXlbomd
    procedure :: addDerivativeMatrix

  end type TShortGamma


contains

  !> Initializes a TShortGamma instance.
  subroutine TShortGamma_init(this, input, orb)

    !> Initialized instance
    type(TShortGamma), intent(out) :: this

    !> Input parameters
    type(TShortGammaInput), intent(inout) :: input

    !> Basis orbital information
    type(TOrbitals), intent(in) :: orb

    this%nSpecies_ = size(orb%nOrbSpecies)
    this%nAtom_ = size(orb%nOrbAtom)
    this%mShell_ = orb%mShell

    call move_alloc(input%hubbU, this%hubbU_)
    allocate(this%shortCutOffs_(this%hubbU_%mHubbU, this%hubbU_%mHubbU, this%nSpecies_,&
        & this%nSpecies_))
    call setupShortGammaCutoffs_(this%hubbU_, this%shortCutoffs_)

    allocate(this%nNeigh_(this%hubbU_%mHubbU, this%hubbU_%mHubbU, this%nSpecies_, this%nAtom_))
    allocate(this%deltaQShell_(this%mShell_, this%nAtom_))
    allocate(this%deltaQUniqU_(this%hubbU_%mHubbU, this%nAtom_))

    if (allocated(input%h5CorrectionInp)) then
      allocate(this%h5Correction_)
      call TH5Correction_init(this%h5Correction_, input%h5CorrectionInp)
    end if

    if (allocated(input%damping)) then
      call move_alloc(input%damping, this%damping_)
    else
      allocate(this%damping_)
      allocate(this%damping_%isDamped(this%nSpecies_))
      this%damping_%isDamped(:) = .false.
      this%damping_%exponent = 0.0_dp
    end if

    allocate(this%shiftShell_(this%mShell_, this%nAtom_))
    allocate(this%shortGamma_(0, 0, 0, 0))

  end subroutine TShortGamma_init


  !> Returns the real space cutoff needed for the neighbour list.
  function getCutoff(this) result(cutoff)

    !> Instance
    class(TShortGamma), intent(in) :: this

    !> Cutoff until which neighbours are required by the instance
    real(dp) :: cutoff

    cutoff = maxval(this%shortCutoffs_)

  end function getCutoff


  !> Updates the coordinates.
  subroutine updateCoords(this, coords, species, neighList)

    !> Instance
    class(TShortGamma), intent(inout) :: this

    !> Coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Species of each atom. Shape: [nAtom]
    integer, intent(in) :: species(:)

    !> Neighbour list
    type(TNeighbourList), intent(in) :: neighList

    call updateNrOfNeighbours_(this%shortCutoffs_, species, this%hubbU_, neighList, this%nNeigh_)
    call updateShortGammaValues_(coords, species, neighList, this%nNeigh_, this%hubbU_,&
        & this%damping_, this%h5Correction_, this%shortGamma_)

  end subroutine updateCoords


  !> Updates the charges.
  subroutine updateCharges(this, orb, species, deltaQShell)

    !> Instance
    class(TShortGamma), intent(inout) :: this

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> Shell resolved net charges. Shape [mShell, nAtom]
    real(dp), intent(in) :: deltaQShell(:,:)

    this%deltaQShell_(:,:) = deltaQShell
    call this%hubbU_%sumOverUniqueU(deltaQShell, species, orb, this%deltaQUniqU_)

  end subroutine updateCharges


  !> Instructs the calculator to recalculate its shift vectors.
  subroutine updateShifts(this, env, orb, species, iNeighbour, img2CentCell)

    !> Instance
    class(TShortGamma), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> Index of the neighbours
    integer, intent(in) :: iNeighbour(0:,:)

    !> Mapping of the atoms into the central cell
    integer, intent(in) :: img2CentCell(:)

    call buildShifts_(env, orb, species, iNeighbour, img2CentCell, this%hubbU_, this%nNeigh_,&
        & this%shortGamma_, this%deltaQShell_, this%shiftShell_)

  end subroutine updateShifts


  !> Returns the shifts calculated.
  subroutine getShiftPerShell(this, shiftShell)

    !> Instance
    class(TShortGamma), intent(inout) :: this

    !> Shell resolved shift vectors. Shape: [mShell, nAtom]
    real(dp), intent(out) :: shiftShell(:,:)

    shiftShell(:,:) = this%shiftShell_

  end subroutine getShiftPerShell


  !> Adds the atom resolved interaction matrix.
  subroutine addAtomicMatrix(this, gammamat, iNeighbour, img2CentCell, iAt, jAt)

    !> Instance
    class(TShortGamma), intent(in) :: this

    !> Atomic interaction matrix
    real(dp), intent(inout) :: gammamat(:,:)

    !> Index of the neighbours
    integer, intent(in) :: iNeighbour(0:,:)

    !> Mapping of the atoms into the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Optional pair of atoms for which to evaluate and add interaction
    integer, intent(in), optional :: iAt, jAt

    integer :: iAt1, iAt2f, iNeigh

    if (present(jAt)) then
      @:ASSERT(present(iAt))
      @:ASSERT(iAt >= jAt)
      do iNeigh = 0, maxval(this%nNeigh_(:,:,:, jAt))
        iAt2f = img2CentCell(iNeighbour(iNeigh, jAt))
        if (iAt2f == iAt) then
          gammamat(iAt, jAt) = gammamat(iAt, jAt) - this%shortGamma_(1, 1, iNeigh, jAt)
        end if
      end do
    else
      do iAt1 = 1, this%nAtom_
        do iNeigh = 0, maxval(this%nNeigh_(:,:,:, iAt1))
          iAt2f = img2CentCell(iNeighbour(iNeigh, iAt1))
          gammamat(iAt2f, iAt1) = gammamat(iAt2f, iAt1) - this%shortGamma_(1, 1, iNeigh, iAt1)
        end do
      end do
    end if

  end subroutine addAtomicMatrix


  !> Adds the atom resolved interaction matrix using customized Us.
  !!
  !! Note: Although the customized U-values are customised passed as arguments, the neighbours
  !! considered during the calculation are determined by the U values which had been used to
  !! initialize this instance. If the two sets of U-values differ significantly, the results
  !! may be very inaccurate.
  subroutine addAtomicMatrixCustomU(this, gammamat, hubbU, species, coords, iNeighbour,&
      & img2CentCell)

    !> Instance
    class(TShortGamma), intent(in) :: this

    !> Atomic interaction matrix. Shape [nAtom, nAtom]
    real(dp), intent(inout) :: gammamat(:,:)

    !> Customized hubbard U values to use
    real(dp), intent(in) :: hubbU(:)

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> Coordinates of the atoms. Shape [3, nAtom]
    real(dp), intent(in) :: coords(:,:)

    !> Neighbour indices
    integer, intent(in) :: iNeighbour(0:,:)

    !> Mapping of the atoms into the central cell
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iAt2, iAt2f, iSp1, iSp2, iNeigh
    real(dp) :: dist

    do iAt1 = 1, this%nAtom_
      iSp1 = species(iAt1)
      do iNeigh = 0, maxval(this%nNeigh_(:,:,:, iAt1))
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        dist = sqrt(sum((coords(:,iAt1) - coords(:,iAt2))**2))
        gammamat(iAt2f, iAt1) = gammamat(iAt2f, iAt1) - expGamma(dist, hubbU(iSp2), hubbU(iSp1))
      end do
    end do

  end subroutine addAtomicMatrixCustomU


  !> Adds the "double counting" contribution of the gradient (the gradient of the field).
  subroutine addGradientsDc(this, force, species, coords, iNeighbour, img2CentCell)

    !> Resulting module variables
    class(TShortGamma), intent(in) :: this

    !> Force vector to add the short-range part of gamma contribution
    real(dp), intent(inout) :: force(:,:)

    !> List of the species for each atom
    integer, intent(in) :: species(:)

    !> Atom coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Index of neighbouring atoms for each atom
    integer, intent(in) :: iNeighbour(0:,:)

    !> Image of each atom in the central cell
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iAt2, iAt2f, iU1, iU2, iNeigh, ii, iSp1, iSp2
    real(dp) :: rab, tmpGammaPrime, u1, u2
    real(dp) :: tmpGamma

    ! some additional symmetry not used
    do iAt1 = 1, this%nAtom_
      iSp1 = species(iAt1)
      do iNeigh = 1, maxval(this%nNeigh_(:,:,:, iAt1))
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        rab = sqrt(sum((coords(:,iAt1) - coords(:,iAt2))**2))
        do iU1 = 1, this%hubbU_%nHubbU(species(iAt1))
          u1 = this%hubbU_%uniqHubbU(iU1, iSp1)
          do iU2 = 1, this%hubbU_%nHubbU(species(iAt2f))
            u2 = this%hubbU_%uniqHubbU(iU2, species(iAt2f))
            if (iNeigh <= this%nNeigh_(iU2, iU1, species(iAt2f), iAt1)) then
              if (this%damping_%isDamped(iSp1) .or. this%damping_%isDamped(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, this%damping_%exponent)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
                if (allocated(this%h5Correction_)) then
                  tmpGamma = expGamma(rab, u2, u1)
                  call this%h5Correction_%scaleShortGammaDeriv(tmpGamma, tmpGammaPrime, iSp1, iSp2,&
                      & rab)
                end if
              end if
              do ii = 1,3
                force(ii, iAt1) = force(ii, iAt1) - this%deltaQUniqU_(iU1, iAt1)&
                    & * this%deltaQUniqU_(iU2, iAt2f) * tmpGammaPrime&
                    & * (coords(ii, iAt1) - coords(ii, iAt2)) / rab
                force(ii, iAt2f) = force(ii, iAt2f) + this%deltaQUniqU_(iU1, iAt1)&
                    & * this%deltaQUniqU_(iU2, iAt2f) * tmpGammaPrime&
                    & * (coords(ii, iAt1) - coords(ii, iAt2)) / rab
              end do
            end if
          end do
        end do
      end do
    end do

  end subroutine addGradientsDc


  !> Adds the "double counting" contribution of the stress (due to the gradient of the field).
  subroutine addStressDc(this, stress, coords, species, volume, iNeighbour, img2CentCell)

    !> Instance
    class(TShortGamma), intent(in) :: this

    !> Stress tensor which value should be added to
    real(dp), intent(inout) :: stress(:,:)

    !> Coordinates of the atoms
    real(dp), intent(in) :: coords(:,:)

    !> Species of the atoms
    integer, intent(in) :: species(:)

    !> Volume of the system
    real(dp), intent(in) :: volume

    !> Index of neighbouring atoms for each atom
    integer, intent(in) :: iNeighbour(0:,:)

    !> Image of each atom in the central cell
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iAt2, iAt2f, iU1, iU2, iNeigh, ii, jj, iSp1, iSp2
    real(dp) :: rab, tmpGammaPrime, u1, u2
    real(dp) :: intermed(3), vect(3)
    real(dp) :: tmpGamma

    real(dp) :: st(3, 3)

    st(:,:) = 0.0_dp
    ! some additional symmetry not used
    do iAt1 = 1, this%nAtom_
      iSp1 = species(iAt1)
      do iNeigh = 1, maxval(this%nNeigh_(:,:,:, iAt1))
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        vect(:) = coords(:,iAt1) - coords(:,iAt2)
        rab = sqrt(sum((vect)**2))
        intermed(:) = 0.0_dp
        do iU1 = 1, this%hubbU_%nHubbU(species(iAt1))
          u1 = this%hubbU_%uniqHubbU(iU1, iSp1)
          do iU2 = 1, this%hubbU_%nHubbU(species(iAt2f))
            u2 = this%hubbU_%uniqHubbU(iU2, species(iAt2f))
            if (iNeigh <= this%nNeigh_(iU2, iU1, species(iAt2f), iAt1)) then
              if (this%damping_%isDamped(iSp1) .or. this%damping_%isDamped(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, this%damping_%exponent)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
                if (allocated(this%h5Correction_)) then
                  tmpGamma = expGamma(rab, u2, u1)
                  call this%h5Correction_%scaleShortGammaDeriv(tmpGamma, tmpGammaPrime, iSp1, iSp2,&
                      & rab)
                end if
              end if
              do ii = 1,3
                intermed(ii) = intermed(ii)&
                    & - this%deltaQUniqU_(iU1, iAt1) * this%deltaQUniqU_(iU2, iAt2f)&
                    & * tmpGammaPrime * vect(ii) / rab
              end do
            end if
          end do
        end do
        if (iAt2f /= iAt1) then
          do ii = 1, 3
            do jj = 1, 3
              st(jj,ii) = st(jj,ii) + (vect(jj) * intermed(ii) + intermed(jj) * vect(ii))
            end do
          end do
        else
          do ii = 1, 3
            do jj = 1, 3
              st(jj,ii) = st(jj,ii) + 0.5_dp * (vect(jj) * intermed(ii) + intermed(jj) * vect(ii))
            end do
          end do
        end if
      end do
    end do

    st(:,:) = -0.5_dp * st / volume
    stress(:,:) = stress + st

  end subroutine addStressDc


  !> Calculate the derivative of the short range contributions using the linearized XLBOMD
  !! formulation with auxiliary charges.
  subroutine addGradientsDcXlbomd(this, dQInShell, dQOutShell, coords, species, orb, iNeighbour,&
      & img2CentCell, force)

    !> Instance
    class(TShortGamma), intent(in) :: this

    !> Input charges
    real(dp), intent(in) :: dQInShell(:,:)

    !> Output charges
    real(dp), intent(in) :: dQOutShell(:,:)

    !> Chemical species
    integer, intent(in) :: species(:)

    !> Atom coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Neighbours around atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Image to real atom indexing
    integer, intent(in) :: img2CentCell(:)

    !> Term to add force contributions to
    real(dp), intent(inout) :: force(:,:)

    integer :: iAt1, iAt2, iAt2f, iU1, iU2, iNeigh, iSp1, iSp2
    real(dp) :: rab, tmpGammaPrime, u1, u2, prefac
    real(dp) :: contrib(3)
    real(dp), allocatable :: dQInUniqU(:,:), dQOutUniqU(:,:)

    allocate(dQInUniqU(this%hubbU_%mHubbU, this%nAtom_))
    call this%hubbU_%sumOverUniqueU(dQInShell, species, orb, dQInUniqU)
    allocate(dQOutUniqU(this%hubbU_%mHubbU, this%nAtom_))
    call this%hubbU_%sumOverUniqueU(dQOutShell, species, orb, dQOutUniqU)

    do iAt1 = 1, this%nAtom_
      iSp1 = species(iAt1)
      do iNeigh = 1, maxval(this%nNeigh_(:,:,:, iAt1))
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        rab = sqrt(sum((coords(:,iAt1) - coords(:,iAt2))**2))
        do iU1 = 1, this%hubbU_%nHubbU(species(iAt1))
          u1 = this%hubbU_%uniqHubbU(iU1, iSp1)
          do iU2 = 1, this%hubbU_%nHubbU(species(iAt2f))
            u2 = this%hubbU_%uniqHubbU(iU2, species(iAt2f))
            if (iNeigh <= this%nNeigh_(iU2, iU1, species(iAt2f), iAt1)) then
              if (this%damping_%isDamped(iSp1) .or. this%damping_%isDamped(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, this%damping_%exponent)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
              end if
              prefac = dQOutUniqU(iU1, iAt1) * dQInUniqU(iU2, iAt2f)&
                  & + dQInUniqU(iU1, iAt1) * dQOutUniqU(iU2, iAt2f)&
                  & - dQInUniqU(iU1, iAt1) * dQInUniqU(iU2, iAt2f)
              contrib(:) = prefac * tmpGammaPrime / rab  * (coords(:,iAt2) - coords(:,iAt1))
              force(:,iAt1) = force(:,iAt1) + contrib
              force(:,iAt2f) = force(:,iAt2f) - contrib
            end if
          end do
        end do
      end do
    end do

  end subroutine addGradientsDcXlbomd


  !> Adds the atom pair resolved gradients.
  subroutine addDerivativeMatrix(this, coords, species, iNeighbour, img2CentCell, derivMat)

    !> Instance
    class(TShortGamma), intent(in) :: this

    !> Coordinates of the atoms
    real(dp), intent(in) :: coords(:,:)

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> Index of neighbour atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Mapping of atoms into central cell
    integer, intent(in) :: img2CentCell(:)

    !> Atom-pair resolved derivative matrix. Shape: [nAtom, nAtom, 3]
    real(dp), intent(inout) :: derivMat(:,:,:)

    integer :: iAt1, iAt2, iAt2f, iSp1, iSp2, iNeigh, iU1, iU2, ii
    real(dp) :: rab, tmpGamma, tmpGammaPrime, u1, u2

    do iAt1 = 1, this%nAtom_
      iSp1 = species(iAt1)
      do iNeigh = 1, maxval(this%nNeigh_(:,:,:, iAt1))
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        rab = sqrt(sum((coords(:,iAt1) - coords(:,iAt2))**2))
        do iU1 = 1, this%hubbU_%nHubbU(species(iAt1))
          u1 = this%hubbU_%uniqHubbU(iU1, iSp1)
          do iU2 = 1, this%hubbU_%nHubbU(species(iAt2f))
            u2 = this%hubbU_%uniqHubbU(iU2, species(iAt2f))
            if (iNeigh <= this%nNeigh_(iU2, iU1, species(iAt2f), iAt1)) then
              if (this%damping_%isDamped(iSp1) .or. this%damping_%isDamped(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, this%damping_%exponent)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
                if (allocated(this%h5Correction_)) then
                  tmpGamma = expGamma(rab, u2, u1)
                  call this%h5Correction_%scaleShortGammaDeriv(tmpGamma, tmpGammaPrime, iSp1, iSp2,&
                      & rab)
                end if
              end if
              do ii = 1,3
                derivMat(iAt2f, iAt1, ii) = -tmpGammaPrime * (coords(ii, iAt1) - coords(ii, iAt2))&
                    & / rab
              end do
            end if
          end do
        end do
      end do
    end do

  end subroutine addDerivativeMatrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !> Sets up cutoffs for the short gamma interaction.
  subroutine setupShortGammaCutoffs_(hubb, shortCutoffs)
    type(TUniqueHubbard), intent(in) :: hubb
    real(dp), intent(out) :: shortCutoffs(:,:,:,:)

    integer :: nSpecies
    integer :: iSp1, iSp2, iU1, iU2

    nSpecies = size(hubb%nHubbU)

    shortCutoffs(:,:,:,:) = 0.0_dp
    do iSp1 = 1, nSpecies
      do iSp2 = iSp1, nSpecies
        do iU1 = 1, hubb%nHubbU(iSp1)
          do iU2 = 1, hubb%nHubbU(iSp2)
            shortCutoffs(iU2, iU1, iSp2, iSp1) =&
                & expGammaCutoff(hubb%uniqHubbU(iU2, iSp2), hubb%uniqHubbU(iU1, iSp1))
            shortCutoffs(iU1, iU2, iSp1, iSp2) = shortCutoffs(iU2, iU1, iSp2, iSp1)
          end do
        end do
      end do
    end do

  end subroutine setupShortGammaCutoffs_


  !> Updates the number of neighbours.
  subroutine updateNrOfNeighbours_(shortCutoffs, species, hubb, neighList, nNeigh)
    real(dp), intent(in) :: shortCutoffs(:,:,:,:)
    integer, intent(in) :: species(:)
    type(TUniqueHubbard), intent(in) :: hubb
    type(TNeighbourList), intent(in) :: neighList
    integer, intent(out) :: nNeigh(:,:,:,:)

    integer :: iAt, iSp, iU1, iU2

    nNeigh(:,:,:,:) = 0
    do iAt = 1, size(nNeigh, dim=4)
      do iSp = 1, size(nNeigh, dim=3)
        do iU1 = 1, hubb%nHubbU(species(iAt))
          do iU2 = 1, hubb%nHubbU(iSp)
            nNeigh(iU2, iU1, iSp, iAt) =&
                & getNrOfNeighbours(neighList, shortCutoffs(iU2, iU1, iSp, species(iAt)), iAt)
          end do
        end do
      end do
    end do

  end subroutine updateNrOfNeighbours_


  !> Set up the storage and internal values for the short range part of Gamma.
  subroutine updateShortGammaValues_(coords, species, neighList, nNeigh, hubb, damping,&
      & h5Correction, shortGamma)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    type(TNeighbourList), intent(in) :: neighList
    integer, intent(in) :: nNeigh(:,:,:,:)
    type(TUniqueHubbard), intent(in) :: hubb
    type(TShortGammaDamp), intent(in) :: damping
    type(TH5Correction), intent(in), optional :: h5Correction
    real(dp), allocatable, intent(inout) :: shortGamma(:,:,:,:)

    integer :: nAtom
    integer :: iAt1, iAt2, iU1, iU2, iNeigh, iSp1, iSp2
    real(dp) :: rab, u1, u2

    nAtom = size(nNeigh, dim=4)

    ! Reallocate shortGamma, if it does not contain enough neighbours
    if (size(shortGamma, dim=3) < maxval(nNeigh) + 1) then
      deallocate(shortGamma)
      allocate(shortGamma(hubb%mHubbU, hubb%mHubbU, 0 : maxval(nNeigh), nAtom))
    end if
    shortGamma(:,:,:,:) = 0.0_dp

    ! some additional symmetry not used, as the value of gamma for atoms
    ! interacting with themselves is the same for all atoms of the same species
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      do iNeigh = 0, maxval(nNeigh(:,:,:, iAt1))
        iAt2 = neighList%iNeighbour(iNeigh, iAt1)
        iSp2 = species(iAt2)
        rab = sqrt(sum((coords(:, iAt1) - coords(:, iAt2))**2))
        do iU1 = 1, hubb%nHubbU(species(iAt1))
          u1 = hubb%uniqHubbU(iU1, iSp1)
          do iU2 = 1, hubb%nHubbU(species(iAt2))
            u2 = hubb%uniqHubbU(iU2, iSp2)
            if (iNeigh <= nNeigh(iU2, iU1, iSp2, iAt1)) then
              if (damping%isDamped(iSp1) .or. damping%isDamped(iSp2)) then
                shortGamma(iU2, iU1, iNeigh, iAt1) = expGammaDamped(rab, u2, u1, damping%exponent)
              else
                shortGamma(iU2, iU1, iNeigh, iAt1) = expGamma(rab, u2, u1)
                if (present(h5Correction)) then
                  call h5Correction%scaleShortGamma(shortGamma(iU2, iU1, iNeigh, iAt1), iSp1, iSp2,&
                      & rab)
                end if
              end if
            end if
          end do
        end do
      end do
    end do

  end subroutine updateShortGammaValues_


  !> Builds the short range shell resolved part of the shift vector.
  subroutine buildShifts_(env, orb, species, iNeighbours, img2CentCell, hubb, nNeigh, shortGamma,&
      & deltaQShell, shiftShell)

    type(TEnvironment), intent(in) :: env
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbours(0:,:)
    integer, intent(in) :: img2CentCell(:)
    type(TUniqueHubbard), intent(in) :: hubb
    integer, intent(in) :: nNeigh(:,:,:,:)
    real(dp), intent(in) :: shortGamma(:,:,0:,:), deltaQShell(:,:)
    real(dp), intent(out) :: shiftShell(:,:)

    integer :: iAt1, iAt2f, iSp1, iSp2, iSh1, iSh2, iU1, iU2, iIter, iNeigh
    integer :: nAtom, maxShell1, maxShell2
    integer, allocatable :: maxNeigh(:)
    real(dp) :: shortGammaValue, deltaQShellValue
    integer, allocatable :: iterIndices(:)
    logical :: skipLoop

    nAtom = size(nNeigh, dim=4)

    allocate(maxNeigh(nAtom))
    do iAt1 = 1, nAtom
      maxNeigh(iAt1) = maxval(nNeigh(:,:,:,iAt1))
    end do

    skipLoop = .false.
    #:if WITH_SCALAPACK
      if (env%blacs%atomGrid%iProc /= -1) then
        call getIndicesWithWorkload(env%blacs%atomGrid%nproc, env%blacs%atomGrid%iproc,&
            & 1, nAtom, maxNeigh, iterIndices)
      else
        ! Do not calculate anything if process is not part of the atomic grid
        skipLoop = .true.
      end if
    #:else
      allocate(iterIndices(nAtom))
      iterIndices(:) = [(iIter, iIter = 1, nAtom)]
    #:endif

    shiftShell(:,:) = 0.0_dp

    if (.not. skipLoop) then
      do iIter = 1, size(iterIndices)
        iAt1 = iterIndices(iIter)
        iSp1 = species(iAt1)
        maxShell1 = orb%nShell(iSp1)
        do iSh1 = 1, maxShell1
          iU1 = hubb%iHubbU(iSh1, iSp1)
          deltaQShellValue = deltaQShell(iSh1, iAt1)
          do iNeigh = 0, maxNeigh(iAt1)
            iAt2f = img2CentCell(iNeighbours(iNeigh, iAt1))
            iSp2 = species(iAt2f)
            maxShell2 = orb%nShell(iSp2)
            do iSh2 = 1, maxShell2
              iU2 = hubb%iHubbU(iSh2, iSp2)
              shortGammaValue = shortGamma(iU2, iU1, iNeigh, iAt1)
              shiftShell(iSh1, iAt1) = shiftShell(iSh1, iAt1)&
                  & - shortGammaValue * deltaQShell(iSh2, iAt2f)
              if (iAt2f /= iAt1) then
                shiftShell(iSh2, iAt2f) = shiftShell(iSh2, iAt2f)&
                    & - shortGammaValue * deltaQShellValue
              end if
            end do
          end do
        end do
      end do
    end if

    #:if WITH_SCALAPACK
      call mpifx_allreduceip(env%mpi%groupComm, shiftShell, MPI_SUM)
    #:endif

    end subroutine buildShifts_

end module dftbp_dftb_shortgamma
