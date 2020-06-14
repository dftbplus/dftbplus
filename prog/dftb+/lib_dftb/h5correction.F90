!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> H5 H-bond correction. Scales the gamma function at short-range for H-bond acceptor element pairs.
!> See http://dx.doi.org/10.1021/acs.jctc.7b00629 for details.
module dftbp_h5correction
  use dftbp_accuracy
  use dftbp_vdwdata
  use dftbp_message, only : warning
  implicit none
  private

  public :: TH5Corr, H5Corr_init

  !> Internal data of the H5 correction
  type :: TH5Corr
    private

    !> distance scale of correction
    real(dp) :: rScale

    !> width scale of correction
    real(dp) :: wScale

    !> Pairwise H5 scaling factors (negative value indicates no scaling)
    real(dp), allocatable :: h5Scaling(:,:)

    !> Pairwise sum of the Van der Waals radii (negative value indicates no scaling)
    real(dp), allocatable :: sumVdw(:,:)

  contains

    procedure :: scaleShortGamma
    procedure :: scaleShortGammaDeriv

  end type TH5Corr


  ! Conversion from full-width-at-half-maximum to c 2.35482 == 2*sqrt(2*ln(2))
  real(dp), parameter :: gaussianWidthFactor = 0.5_dp / sqrt(2.0_dp * log(2.0_dp))

contains

  !> Initialization of a H5Corr instance.
  subroutine H5Corr_init(this, speciesNames, rr, ww, elementParams)

    !> Initialised instance at return.
    type(TH5Corr), intent(out) :: this

    !> Names of the species
    character(mc), allocatable, intent(in) :: speciesNames(:)

    !> r scaling factor, if set to -1 is replaced with value from the paper
    real(dp), intent(in) :: rr

    !> w scaling factor, if set to -1 is replaced with value from the paper
    real(dp), intent(in) :: ww

    !> Elementwise scaling factors, if set to -1 (or any negative value) species is not corrected
    real(dp), allocatable, intent(in) :: elementParams(:)

    integer :: nSpecies, iSp

    this%rScale = rr
    this%wScale = ww

    nSpecies = size(speciesNames)
    allocate(this%sumVdw(nSpecies, nSpecies))
    allocate(this%h5scaling(nSpecies, nSpecies))
    call getParams(speciesNames, elementParams, this%h5Scaling, this%sumVdw)

  end subroutine H5Corr_init


  !> Apply the correction to the short-range part of gamma function, returning modified shortGamma
  subroutine scaleShortGamma(this, shortGamma, iSp1, iSp2, rab)

    !> instance of the correction
    class(TH5Corr), intent(in) :: this

    !> short range gamma value
    real(dp), intent(inout) :: shortGamma

    !> species of first atom in the pair
    integer, intent(in) :: iSp1

    !> species of the second atom in the pair
    integer, intent(in) :: iSp2

    !> separation between atoms
    real(dp), intent(in) :: rab

    real(dp) :: h5Scaling, gauss, sumVdw, fwhm, r0, cc

    sumVdw = this%sumVdw(iSp2, iSp1)
    h5Scaling = this%h5Scaling(iSp2, iSp1)
    if (h5Scaling > 0.0_dp) then
      fwhm = this%wScale * sumVdw
      r0 = this%rScale * sumVdw
      cc = fwhm * gaussianWidthFactor
      gauss = exp(-0.5_dp * (rab - r0)**2 / cc**2) * h5Scaling
      shortGamma = shortGamma * (1.0_dp + gauss) - gauss / rab
    end if

  end subroutine scaleShortGamma


  !> Apply the correction to the derivative of the short-range part of gamma function
  subroutine scaleShortGammaDeriv(this, shortGamma, shortGammaDeriv, iSp1, iSp2, rab)

    !> instance of the correction
    class(TH5Corr), intent(in) :: this

    !> short range gamma
    real(dp), intent(in) :: shortGamma

    !> derivative of short range gamma to scale
    real(dp), intent(inout) :: shortGammaDeriv

    !> species of first atom in the pair
    integer, intent(in) :: iSp1

    !> species of second atom in the pair
    integer, intent(in) :: iSp2

    !> separation of pair
    real(dp), intent(in) :: rab

    ! Local variables
    real(dp) :: h5Scaling, gauss, sumVdw, fwhm, r0, cc, dgauss, deriv1, deriv2

    sumVdw = this%sumVdw(iSp2, iSp1)
    h5Scaling = this%h5Scaling(iSp2, iSp1)

    if (h5Scaling > 0.0_dp) then
      fwhm = this%wScale * sumVdw
      r0 = this%rScale * sumVdw
      cc = fwhm * gaussianWidthFactor
      gauss = exp(-0.5_dp * (rab - r0)**2 / cc**2) * h5Scaling

      dgauss = -gauss * (rab - r0) / cc**2
      deriv1 = shortGamma * dgauss + shortGammaDeriv * (1.0_dp + gauss)
      deriv2 = dgauss / rab - (h5Scaling * exp(-0.5_dp * ((rab - r0)**2) / cc**2)) / rab**2
      shortGammaDeriv = deriv1 - deriv2
    end if

  end subroutine scaleShortGammaDeriv


  !> Get H5 parameters for all species pairs.
  subroutine getParams(speciesNames, elementParams, h5Scaling, sumVdw)

    !> Names of the atom types
    character(*), intent(in) :: speciesNames(:)

    !> H5-parameters for the elements
    real(dp), intent(in) :: elementParams(:)

    !> Pair-specific scaling factor (negative if scaling should be ommitted)
    real(dp), intent(out) :: h5Scaling(:,:)

    !> Sum of vdW radii of the pair (negative if scaling should be ommitted)
    real(dp), intent(out) :: sumVdw(:,:)

    character(mc) :: spName1, spName2
    integer :: iSpHeavy
    real(dp) :: vdwH, vdwHeavy
    logical :: tFoundRadius
    integer :: iSp1, iSp2

    call getVdwData("H", vdwH)

    h5scaling(:,:) = -1.0_dp
    sumVdw(:,:) = -1.0_dp
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

        h5Scaling(iSp2, iSp1) = elementParams(iSpHeavy)

        call getVdwData(speciesNames(iSpHeavy), vdwHeavy, found=tFoundRadius)
        if (.not. tFoundRadius .and. h5Scaling(iSp2, iSp1) > 0.0_dp) then
          call warning("The van de  Waals radius for " // trim(speciesNames(iSpHeavy)) //&
              & " is required for the H5 correction but is not available. H-" //&
              & trim(speciesNames(iSpHeavy)) // " contributions therefore neglected.")
          h5Scaling(iSp2, iSp1) = -1.0_dp
        end if

        if (h5Scaling(iSp2, iSp1) > 0.0_dp) then
          sumVdw(iSp2, iSp1) = vdwH + vdwHeavy
        end if
      end do
    end do

  end subroutine getParams


end module dftbp_h5correction
