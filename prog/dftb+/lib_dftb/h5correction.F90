!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> H5 H-bond correction. Scales the gamma function at short-range for H-bond acceptor element pairs.
!> See http://dx.doi.org/10.1021/acs.jctc.7b00629 for details.
module h5correction
  use accuracy
  use vdwdata
  use message, only : warning
  implicit none
  private

  public :: H5Corr, H5Corr_init

  !> Internal data of the H5 correction
  type :: H5Corr

    !> distance scale of correction
    real(dp) :: rScale

    !> width scale of correction
    real(dp) :: wScale

    !> Nr. of species
    integer :: nSpecies

    !> Species names for H5 correction
    character(mc), allocatable :: speciesName(:)

    !> Elementwise parameters for species
    real(dp), allocatable :: elementPara(:)

  contains

    procedure :: writeH5Setup
    procedure :: getParams
    procedure :: scaleShortGamma
    procedure :: scaleShortGammaDeriv

  end type H5Corr


  ! Conversion from full-width-at-half-maximum to c 2.35482 == 2*sqrt(2*ln(2))
  real(dp), parameter :: gaussianWidthFactor = 0.5_dp / sqrt(2.0_dp * log(2.0_dp))

contains

  !> Initialization of a H5Corr instance.
  subroutine H5Corr_init(this, nType, typeNames, rr, ww, elePara)

    !> Initialised instance at return.
    type(H5Corr), intent(out) :: this

    !> Number of atom types/species
    integer, intent(in) :: nType

    !> Names of the species
    character(mc), allocatable, intent(in) :: typeNames(:)

    !> r scaling factor, if set to -1 is replaced with value from the paper
    real(dp), intent(in) :: rr

    !> w scaling factor, if set to -1 is replaced with value from the paper
    real(dp), intent(in) :: ww

    !> elementwise scaling factors, if set to -1 are replaced with values from the paper in case of
    !> (O,N,S), but zeroed otherwise
    real(dp), allocatable, intent(in) :: elePara(:)

    integer :: iSp

    this%nSpecies = nType
    this%speciesName = typeNames
    this%rScale = rr
    this%wScale = ww
    this%elementPara = elePara

    ! If value of the parameter is -1.0, it was not read from the input and a default value will be
    ! used
    if (this%rScale == -1.0_dp) then
      this%rScale = 0.714_dp
    end if

    ! If value of the parameter is -1.0, it was not read from the input and a default value will be
    ! used
    if (this%wScale == -1.0_dp) then
      this%wScale = 0.25_dp
    end if

    do iSp = 1, this%nSpecies
      if (this%elementPara(iSp) == -1.0_dp) then
        ! Default values taken from the paper for O,N and S, zero for other elements
        select case (this%speciesName(iSp))
        case ("O")
          this%elementPara(iSp) =  0.06_dp
        case ("N")
          this%elementPara(iSp) = 0.18_dp
        case ("S")
          this%elementPara(iSp) =  0.21_dp
        case default
          this%elementPara(iSp) = 0.0_dp
        end select
      end if
    end do

  end subroutine H5Corr_init


  !> Print all the parameters for debugging
  subroutine writeH5Setup(this)

    !> Instance of the correction
    class(H5Corr), intent(in) :: this

    integer :: iSp1
    integer :: h5unit

    open(newunit=h5unit, file='h5_debugging.dat')
    write(h5unit,*) "H5 setup:"
    write(h5unit,*) "   rScale = ", this%rScale
    write(h5unit,*) "   wScale = ", this%wScale
    write(h5unit,*) "H5 species parameters:"
    do iSp1 = 1, this%nSpecies
      write(h5unit,*) "   ", this%speciesName(iSp1), " = ", this%elementPara(iSp1)
    end do
    close(h5unit)

  end subroutine writeH5Setup

  
  !> Get H5 parameters for a pair of species, also returning a flag as to whether the correction is
  !> applied to the pair, and the species-specific parameters
  subroutine getParams(this, iSp1, iSp2, applyCorrection, h5Scaling, sumVdw)

    ! Arguments
    class(H5Corr), intent(in) :: this

    !> Input: species names
    integer, intent(in) :: iSp1

    !> Input: species names
    integer, intent(in) :: iSp2

    !> Output: flag indicating if the correction is applied to that pair
    logical, intent(out) :: applyCorrection

    !> Output: pair-specific scaling factor
    real(dp), intent(out) :: h5Scaling

    !> Output: Sum of vdW radii of the pair
    real(dp), intent(out) :: sumVdw

    character(mc) :: spName1, spName2
    integer :: iSpHeavy
    real(dp) :: vdwH, vdwHeavy
    logical :: tFoundRadii

    spName1 = this%speciesName(iSp1)
    spName2 = this%speciesName(iSp2)

    applyCorrection = .false.

    ! If the pair cannot make an H-bond, return immediately with applyCorrection = .false.
    if (all([spName1, spName2] /= "H")) then
      return
    end if

    ! No heavy atom present in the pair
    if (all([spName1, spName2] == "H")) then
      return
    end if

    ! If there is one hydrogen and one other atom, save the species of the heavy atom
    if (spName1 /= "H") then
      iSpHeavy = iSp1
    else
      iSpHeavy = iSp2
    end if

    call getVdwData("H", vdwH)
    call getVdwData(this%speciesName(iSpHeavy), vdwHeavy, found=tFoundRadii)

    h5Scaling = this%elementPara(iSpHeavy)

    if (.not.tFoundRadii .and. h5Scaling /= 0.0_dp) then
      call warning("The van de  Waals radius for " // trim(this%speciesName(iSpHeavy)) //&
          & " is not available, but is needed for the H5 correction, so is neglected")
    end if

    if (tFoundRadii .and. h5Scaling /= 0.0_dp) then
      applyCorrection = .true.
      sumVdw = vdwH + vdwHeavy
    else
      ! If there are no parameters for this H-bond, return
    end if

  end subroutine getParams


  !> Apply the correction to the short-range part of gamma function, returning modified shortGamma
  subroutine scaleShortGamma(this, shortGamma, iSp1, iSp2, rab)

    !> instance of the correction
    class(H5Corr), intent(in) :: this

    !> short range gamma value
    real(dp), intent(inout) :: shortGamma

    !> species of first atom in the pair
    integer, intent(in) :: iSp1

    !> species of the second atom in the pair
    integer, intent(in) :: iSp2

    !> separation between atoms
    real(dp), intent(in) :: rab

    real(dp) :: h5Scaling, gauss, sumVdw, fwhm, r0, cc
    logical :: applyCorrection

    call this%getParams(iSp1, iSp2, applyCorrection, h5Scaling, sumVdw)

    if (applyCorrection) then
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
    class(H5Corr), intent(in) :: this

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
    logical :: applyCorrection

    ! Get parameters for current pair of species
    call this%getParams(iSp1, iSp2, applyCorrection, h5Scaling, sumVdw)

    ! If applicable to the current pair, modify the gamma
    if (applyCorrection) then
      fwhm = this%wScale * sumVdw
      r0 = this%rScale * sumVdw
      cc = fwhm * gaussianWidthFactor
      gauss = exp(-0.5_dp * (rab - r0)**2 / cc**2) * h5Scaling

      dgauss = -gauss * (rab - r0) / cc**2
      deriv1 = shortGamma * dgauss + shortGammaDeriv * (1.0_dp + gauss)
      deriv2 = dgauss / rab - (h5Scaling * exp(-0.5_dp * ((rab - r0)**2 ) / cc**2 )) / rab**2
      shortGammaDeriv = deriv1 - deriv2
    end if

  end subroutine scaleShortGammaDeriv

end module h5correction
