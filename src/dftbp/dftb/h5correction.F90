!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> H5 H-bond correction. Scales the gamma function at short-range for H-bond acceptor element pairs.
!> See http://dx.doi.org/10.1021/acs.jctc.7b00629 for details.
module dftbp_dftb_h5correction
  use dftbp_common_accuracy, only : dp, mc
  use dftbp_dftb_vdwdata, only : getVdwData
  use dftbp_io_message, only : warning
  implicit none

  private
  public :: TH5CorrectionInput, TH5Correction, TH5Correction_init


  !> Contains the input parameters for the H5 correction
  type :: TH5CorrectionInput

    !> Distance (r) scaling factor, if set to -1 is replaced with value from the paper.
    real(dp) :: rScale

    !> Widths (w) scaling factor, if set to -1 is replaced with value from the paper.
    real(dp) :: wScale

    !> Elementwise scaling factors, if set to -1 (or any negative value) species is not corrected.
    !> Shape: [nSpecies].
    real(dp), allocatable  :: elementParams(:)

    !> Name of the species (needed for vdW parameters based on species names). Shape: [nSpecies].
    character(mc), allocatable :: speciesNames(:)

  end type TH5CorrectionInput


  !> Internal data of the H5 correction
  type :: TH5Correction
    private

    ! distance scale of correction
    real(dp) :: rScale_

    ! width scale of correction
    real(dp) :: wScale_

    ! Pairwise H5 scaling factors (negative value indicates no scaling)
    real(dp), allocatable :: h5Scaling_(:,:)

    ! Pairwise sum of the Van der Waals radii (negative value indicates no scaling)
    real(dp), allocatable :: sumVdw_(:,:)

  contains

    procedure :: scaleShortGamma
    procedure :: scaleShortGammaDeriv

  end type TH5Correction


  ! Conversion from full-width-at-half-maximum to c 2.35482 == 2*sqrt(2*ln(2))
  real(dp), parameter :: gaussianWidthFactor_ = 0.5_dp / sqrt(2.0_dp * log(2.0_dp))

contains

  !> Initialization of a H5Corr instance.
  subroutine TH5Correction_init(this, input)

    !> Initialised instance at return.
    type(TH5Correction), intent(out) :: this

    !> Input parameters
    type(TH5CorrectionInput), intent(inout) :: input

    integer :: nSpecies

    this%rScale_ = input%rScale
    this%wScale_ = input%wScale

    nSpecies = size(input%speciesNames)
    allocate(this%sumVdw_(nSpecies, nSpecies))
    allocate(this%h5Scaling_(nSpecies, nSpecies))
    call getParams_(input%speciesNames, input%elementParams, this%h5Scaling_, this%sumVdw_)

  end subroutine TH5Correction_init


  !> Apply the correction to the short-range part of gamma function, returning modified shortGamma
  subroutine scaleShortGamma(this, shortGamma, iSp1, iSp2, rab)

    !> Instance of the correction
    class(TH5Correction), intent(in) :: this

    !> Short range gamma value
    real(dp), intent(inout) :: shortGamma

    !> Species of first atom in the pair
    integer, intent(in) :: iSp1

    !> Species of the second atom in the pair
    integer, intent(in) :: iSp2

    !> Separation between atoms
    real(dp), intent(in) :: rab

    real(dp) :: h5Scaling_, gauss, sumVdw_, fwhm, r0, cc

    sumVdw_ = this%sumVdw_(iSp2, iSp1)
    h5Scaling_ = this%h5Scaling_(iSp2, iSp1)
    if (h5Scaling_ > 0.0_dp) then
      fwhm = this%wScale_ * sumVdw_
      r0 = this%rScale_ * sumVdw_
      cc = fwhm * gaussianWidthFactor_
      gauss = exp(-0.5_dp * (rab - r0)**2 / cc**2) * h5Scaling_
      shortGamma = shortGamma * (1.0_dp + gauss) - gauss / rab
    end if

  end subroutine scaleShortGamma


  !> Apply the correction to the derivative of the short-range part of gamma function
  subroutine scaleShortGammaDeriv(this, shortGamma, shortGammaDeriv, iSp1, iSp2, rab)

    !> Instance of the correction
    class(TH5Correction), intent(in) :: this

    !> Short range gamma
    real(dp), intent(in) :: shortGamma

    !> Derivative of short range gamma to scale
    real(dp), intent(inout) :: shortGammaDeriv

    !> Species of first atom in the pair
    integer, intent(in) :: iSp1

    !> Species of second atom in the pair
    integer, intent(in) :: iSp2

    !> Separation of pair
    real(dp), intent(in) :: rab

    ! Local variables
    real(dp) :: h5Scaling_, gauss, sumVdw_, fwhm, r0, cc, dgauss, deriv1, deriv2

    sumVdw_ = this%sumVdw_(iSp2, iSp1)
    h5Scaling_ = this%h5Scaling_(iSp2, iSp1)

    if (h5Scaling_ > 0.0_dp) then
      fwhm = this%wScale_ * sumVdw_
      r0 = this%rScale_ * sumVdw_
      cc = fwhm * gaussianWidthFactor_
      gauss = exp(-0.5_dp * (rab - r0)**2 / cc**2) * h5Scaling_

      dgauss = -gauss * (rab - r0) / cc**2
      deriv1 = shortGamma * dgauss + shortGammaDeriv * (1.0_dp + gauss)
      deriv2 = dgauss / rab - (h5Scaling_ * exp(-0.5_dp * ((rab - r0)**2) / cc**2)) / rab**2
      shortGammaDeriv = deriv1 - deriv2
    end if

  end subroutine scaleShortGammaDeriv


  ! Get H5 parameters for all species pairs.
  subroutine getParams_(speciesNames, elementParams, h5Scaling_, sumVdw_)
    character(*), intent(in) :: speciesNames(:)
    real(dp), intent(in) :: elementParams(:)
    real(dp), intent(out) :: h5Scaling_(:,:)
    real(dp), intent(out) :: sumVdw_(:,:)

    character(mc) :: spName1, spName2
    integer :: iSpHeavy
    real(dp) :: vdwH, vdwHeavy
    logical :: tFoundRadius
    integer :: iSp1, iSp2

    call getVdwData("H", vdwH)

    h5Scaling_(:,:) = -1.0_dp
    sumVdw_(:,:) = -1.0_dp
    do iSp1 = 1, size(speciesNames)
      do iSp2 = 1, size(speciesNames)
        spName1 = speciesNames(iSp1)
        spName2 = speciesNames(iSp2)

        ! Scaling only needed if exactly one of the pair is a H-atom
        if (count([spName1, spName2] == "H") /= 1) then
          cycle
        end if

        if (spName1 /= "H") then
          iSpHeavy = iSp1
        else
          iSpHeavy = iSp2
        end if

        h5Scaling_(iSp2, iSp1) = elementParams(iSpHeavy)

        call getVdwData(speciesNames(iSpHeavy), vdwHeavy, found=tFoundRadius)
        if (.not. tFoundRadius .and. h5Scaling_(iSp2, iSp1) > 0.0_dp) then
          call warning("The van de  Waals radius for " // trim(speciesNames(iSpHeavy)) //&
              & " is required for the H5 correction but is not available. H-" //&
              & trim(speciesNames(iSpHeavy)) // " contributions therefore neglected.")
          h5Scaling_(iSp2, iSp1) = -1.0_dp
        end if

        if (h5Scaling_(iSp2, iSp1) > 0.0_dp) then
          sumVdw_(iSp2, iSp1) = vdwH + vdwHeavy
        end if
      end do
    end do

  end subroutine getParams_


end module dftbp_dftb_h5correction
