!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Necessary parameters to perform DFT-D4 calculations.
module dftbp_dftd4param
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_constants, only: pi, AA__Bohr
  use dftbp_dftd4refs
  implicit none

  public :: DftD4Calculator, DispDftD4Inp, initializeCalculator
  private

  !> Maximum atomic number allowed in EEQ calculations
  integer, parameter :: max_element_eeq = 86

  !> Maximum atomic number allowed in D4 calculations
  integer, parameter :: max_element_d4 = 118

  !> Maximum allowed number of reference systems, arbitrary choice
  integer, parameter :: max_references = 7

  !> Number of frequencies used in Casimir-Polder integration
  integer, parameter :: imag_frequencies = 23

  !> Damping parameters for DFT-D4 calculation.
  type :: DispDftD4Inp

    !> scaling parameter for dipole-dipole coefficients.
    real(dp) :: s6 = 1.0_dp

    !> scaling parameter for dipole-quadrupole coefficients.
    real(dp) :: s8

    !> scaling parameter for quadrupole-quadrupole coefficients.
    real(dp) :: s10 = 0.0_dp

    !> scaling parameter for non-additive triple dipole coefficients.
    real(dp) :: s9

    !> scaling parameter for <r4>/<r2> expectation value based critical radii.
    real(dp) :: a1

    !> constant offset of critical radii.
    real(dp) :: a2

    !> exponent of for the zero-damping function used for non-addititive triple dipole
    !>  contributions.
    real(dp) :: alpha = 16.0_dp

    !> cutoff radius for dispersion interactions.
    real(dp) :: cutoffInter = 64.0_dp

    !> cutoff radius for CN counting function.
    real(dp) :: cutoffCount = 40.0_dp

    !> cutoff radius for real space Ewald summation.
    real(dp) :: cutoffEwald = 40.0_dp

    !> cutoff radius for three-body interactions.
    real(dp) :: cutoffThree = 40.0_dp

    !> Gaussian weighting factor for interpolation of dispersion coefficients.
    real(dp) :: weightingFactor = 6.0_dp

    !> maximum charge scaling height for partial charge extrapolation.
    real(dp) :: chargeScale = 3.0_dp

    !> charge scaling steepness for partial charge extrapolation.
    real(dp) :: chargeSteepness = 2.0_dp

    !> if > 0 -> manual setting for alpha
    real(dp) :: parEwald = 0.0_dp

    !> Ewald tolerance
    real(dp) :: tolEwald = 0.0_dp

    !> Net charge
    real(dp) :: nrChrg = 0.0_dp

  end type DispDftD4Inp


  type :: DftD4Calculator

    !> scaling parameter for dipole-dipole coefficients.
    real(dp) :: s6 = 1.0_dp

    !> scaling parameter for dipole-quadrupole coefficients.
    real(dp) :: s8

    !> scaling parameter for quadrupole-quadrupole coefficients.
    real(dp) :: s10 = 0.0_dp

    !> scaling parameter for non-additive triple dipole coefficients.
    real(dp) :: s9

    !> scaling parameter for <r4>/<r2> expectation value based critical radii.
    real(dp) :: a1

    !> constant offset of critical radii.
    real(dp) :: a2

    !> exponent of for the zero-damping function used for non-addititive
    !>  triple dipole contributions.
    real(dp) :: alpha = 16.0_dp

    !> Gaussian weighting factor for interpolation of dispersion coefficients.
    real(dp) :: wf

    !> maximum charge scaling height for partial charge extrapolation.
    real(dp) :: ga

    !> charge scaling steepness for partial charge extrapolation.
    real(dp) :: gc

    !> cutoff radius for dispersion interactions.
    real(dp) :: cutoffInter

    !> cutoff radius for CN counting function.
    real(dp) :: cutoffCount

    !> cutoff radius for real space Ewald summation.
    real(dp) :: cutoffEwald

    !> cutoff radius for three-body interactions.
    real(dp) :: cutoffThree

    !> net charge
    real(dp) :: nrChrg = 0.0_dp

    real(dp), allocatable :: chi(:)
    real(dp), allocatable :: gam(:)
    real(dp), allocatable :: rad(:)
    real(dp), allocatable :: kcn(:)

    real(dp), allocatable :: sqrtZr4r2(:)
    real(dp), allocatable :: electronegativity(:)
    real(dp), allocatable :: covalentRadius(:)
    real(dp), allocatable :: chemicalHardness(:)
    real(dp), allocatable :: effectiveNuclearCharge(:)

    integer, allocatable :: numberOfReferences(:)
    integer, allocatable :: atoms(:)
    integer, allocatable :: countNumber(:, :)
    real(dp), allocatable :: referenceCN(:, :)
    real(dp), allocatable :: referenceCharge(:, :)
    real(dp), allocatable :: referenceAlpha(:, :, :)
    real(dp), allocatable :: referenceC6(:, :, :, :)

  end type DftD4Calculator


  !> Element-specific electronegativity for the electronegativity equilibration charges used in
  !>  DFT-D4.
  real(dp), parameter :: chi(max_element_eeq) = [&
    & 1.23695041_dp, 1.26590957_dp, 0.54341808_dp, 0.99666991_dp, 1.26691604_dp, &
    & 1.40028282_dp, 1.55819364_dp, 1.56866440_dp, 1.57540015_dp, 1.15056627_dp, &
    & 0.55936220_dp, 0.72373742_dp, 1.12910844_dp, 1.12306840_dp, 1.52672442_dp, &
    & 1.40768172_dp, 1.48154584_dp, 1.31062963_dp, 0.40374140_dp, 0.75442607_dp, &
    & 0.76482096_dp, 0.98457281_dp, 0.96702598_dp, 1.05266584_dp, 0.93274875_dp, &
    & 1.04025281_dp, 0.92738624_dp, 1.07419210_dp, 1.07900668_dp, 1.04712861_dp, &
    & 1.15018618_dp, 1.15388455_dp, 1.36313743_dp, 1.36485106_dp, 1.39801837_dp, &
    & 1.18695346_dp, 0.36273870_dp, 0.58797255_dp, 0.71961946_dp, 0.96158233_dp, &
    & 0.89585296_dp, 0.81360499_dp, 1.00794665_dp, 0.92613682_dp, 1.09152285_dp, &
    & 1.14907070_dp, 1.13508911_dp, 1.08853785_dp, 1.11005982_dp, 1.12452195_dp, &
    & 1.21642129_dp, 1.36507125_dp, 1.40340000_dp, 1.16653482_dp, 0.34125098_dp, &
    & 0.58884173_dp, 0.68441115_dp, 0.56999999_dp, 0.56999999_dp, 0.56999999_dp, &
    & 0.56999999_dp, 0.56999999_dp, 0.56999999_dp, 0.56999999_dp, 0.56999999_dp, &
    & 0.56999999_dp, 0.56999999_dp, 0.56999999_dp, 0.56999999_dp, 0.56999999_dp, &
    & 0.56999999_dp, 0.87936784_dp, 1.02761808_dp, 0.93297476_dp, 1.10172128_dp, &
    & 0.97350071_dp, 1.16695666_dp, 1.23997927_dp, 1.18464453_dp, 1.14191734_dp, &
    & 1.12334192_dp, 1.01485321_dp, 1.12950808_dp, 1.30804834_dp, 1.33689961_dp, &
    & 1.27465977_dp]

  !> Element-specific chemical hardnesses for the electronegativity equilibration charges used in
  !>  DFT-D4.
  real(dp), parameter :: gam(max_element_eeq) = [&
    &-0.35015861_dp, 1.04121227_dp, 0.09281243_dp, 0.09412380_dp, 0.26629137_dp, &
    & 0.19408787_dp, 0.05317918_dp, 0.03151644_dp, 0.32275132_dp, 1.30996037_dp, &
    & 0.24206510_dp, 0.04147733_dp, 0.11634126_dp, 0.13155266_dp, 0.15350650_dp, &
    & 0.15250997_dp, 0.17523529_dp, 0.28774450_dp, 0.42937314_dp, 0.01896455_dp, &
    & 0.07179178_dp,-0.01121381_dp,-0.03093370_dp, 0.02716319_dp,-0.01843812_dp, &
    &-0.15270393_dp,-0.09192645_dp,-0.13418723_dp,-0.09861139_dp, 0.18338109_dp, &
    & 0.08299615_dp, 0.11370033_dp, 0.19005278_dp, 0.10980677_dp, 0.12327841_dp, &
    & 0.25345554_dp, 0.58615231_dp, 0.16093861_dp, 0.04548530_dp,-0.02478645_dp, &
    & 0.01909943_dp, 0.01402541_dp,-0.03595279_dp, 0.01137752_dp,-0.03697213_dp, &
    & 0.08009416_dp, 0.02274892_dp, 0.12801822_dp,-0.02078702_dp, 0.05284319_dp, &
    & 0.07581190_dp, 0.09663758_dp, 0.09547417_dp, 0.07803344_dp, 0.64913257_dp, &
    & 0.15348654_dp, 0.05054344_dp, 0.11000000_dp, 0.11000000_dp, 0.11000000_dp, &
    & 0.11000000_dp, 0.11000000_dp, 0.11000000_dp, 0.11000000_dp, 0.11000000_dp, &
    & 0.11000000_dp, 0.11000000_dp, 0.11000000_dp, 0.11000000_dp, 0.11000000_dp, &
    & 0.11000000_dp,-0.02786741_dp, 0.01057858_dp,-0.03892226_dp,-0.04574364_dp, &
    &-0.03874080_dp,-0.03782372_dp,-0.07046855_dp, 0.09546597_dp, 0.21953269_dp, &
    & 0.02522348_dp, 0.15263050_dp, 0.08042611_dp, 0.01878626_dp, 0.08715453_dp, &
    & 0.10500484_dp]

  !> Element-specific CN scaling constant for the electronegativity equilibration charges used in
  !>  DFT-D4.
  real(dp), parameter :: kcn(max_element_eeq) = [&
    & 0.04916110_dp, 0.10937243_dp,-0.12349591_dp,-0.02665108_dp,-0.02631658_dp, &
    & 0.06005196_dp, 0.09279548_dp, 0.11689703_dp, 0.15704746_dp, 0.07987901_dp, &
    &-0.10002962_dp,-0.07712863_dp,-0.02170561_dp,-0.04964052_dp, 0.14250599_dp, &
    & 0.07126660_dp, 0.13682750_dp, 0.14877121_dp,-0.10219289_dp,-0.08979338_dp, &
    &-0.08273597_dp,-0.01754829_dp,-0.02765460_dp,-0.02558926_dp,-0.08010286_dp, &
    &-0.04163215_dp,-0.09369631_dp,-0.03774117_dp,-0.05759708_dp, 0.02431998_dp, &
    &-0.01056270_dp,-0.02692862_dp, 0.07657769_dp, 0.06561608_dp, 0.08006749_dp, &
    & 0.14139200_dp,-0.05351029_dp,-0.06701705_dp,-0.07377246_dp,-0.02927768_dp, &
    &-0.03867291_dp,-0.06929825_dp,-0.04485293_dp,-0.04800824_dp,-0.01484022_dp, &
    & 0.07917502_dp, 0.06619243_dp, 0.02434095_dp,-0.01505548_dp,-0.03030768_dp, &
    & 0.01418235_dp, 0.08953411_dp, 0.08967527_dp, 0.07277771_dp,-0.02129476_dp, &
    &-0.06188828_dp,-0.06568203_dp,-0.11000000_dp,-0.11000000_dp,-0.11000000_dp, &
    &-0.11000000_dp,-0.11000000_dp,-0.11000000_dp,-0.11000000_dp,-0.11000000_dp, &
    &-0.11000000_dp,-0.11000000_dp,-0.11000000_dp,-0.11000000_dp,-0.11000000_dp, &
    &-0.11000000_dp,-0.03585873_dp,-0.03132400_dp,-0.05902379_dp,-0.02827592_dp, &
    &-0.07606260_dp,-0.02123839_dp, 0.03814822_dp, 0.02146834_dp, 0.01580538_dp, &
    &-0.00894298_dp,-0.05864876_dp,-0.01817842_dp, 0.07721851_dp, 0.07936083_dp, &
    & 0.05849285_dp]

  !> Element-specific charge widths for the electronegativity equilibration charges used in DFT-D4.
  real(dp), parameter :: rad(max_element_eeq) = [&
    & 0.55159092_dp, 0.66205886_dp, 0.90529132_dp, 1.51710827_dp, 2.86070364_dp, &
    & 1.88862966_dp, 1.32250290_dp, 1.23166285_dp, 1.77503721_dp, 1.11955204_dp, &
    & 1.28263182_dp, 1.22344336_dp, 1.70936266_dp, 1.54075036_dp, 1.38200579_dp, &
    & 2.18849322_dp, 1.36779065_dp, 1.27039703_dp, 1.64466502_dp, 1.58859404_dp, &
    & 1.65357953_dp, 1.50021521_dp, 1.30104175_dp, 1.46301827_dp, 1.32928147_dp, &
    & 1.02766713_dp, 1.02291377_dp, 0.94343886_dp, 1.14881311_dp, 1.47080755_dp, &
    & 1.76901636_dp, 1.98724061_dp, 2.41244711_dp, 2.26739524_dp, 2.95378999_dp, &
    & 1.20807752_dp, 1.65941046_dp, 1.62733880_dp, 1.61344972_dp, 1.63220728_dp, &
    & 1.60899928_dp, 1.43501286_dp, 1.54559205_dp, 1.32663678_dp, 1.37644152_dp, &
    & 1.36051851_dp, 1.23395526_dp, 1.65734544_dp, 1.53895240_dp, 1.97542736_dp, &
    & 1.97636542_dp, 2.05432381_dp, 3.80138135_dp, 1.43893803_dp, 1.75505957_dp, &
    & 1.59815118_dp, 1.76401732_dp, 1.63999999_dp, 1.63999999_dp, 1.63999999_dp, &
    & 1.63999999_dp, 1.63999999_dp, 1.63999999_dp, 1.63999999_dp, 1.63999999_dp, &
    & 1.63999999_dp, 1.63999999_dp, 1.63999999_dp, 1.63999999_dp, 1.63999999_dp, &
    & 1.63999999_dp, 1.47055223_dp, 1.81127084_dp, 1.40189963_dp, 1.54015481_dp, &
    & 1.33721475_dp, 1.57165422_dp, 1.04815857_dp, 1.78342098_dp, 2.79106396_dp, &
    & 1.78160840_dp, 2.47588882_dp, 2.37670734_dp, 1.76613217_dp, 2.66172302_dp, &
    & 2.82773085_dp]

  !> Element-specific chemical hardnesses for the charge scaling function used to extrapolate the C6
  !>  coefficients in DFT-D4.
  real(dp), parameter :: chemicalHardness(max_element_d4) = [ &
    & 0.47259288_dp, 0.92203391_dp, 0.17452888_dp, 0.25700733_dp, 0.33949086_dp, &
    & 0.42195412_dp, 0.50438193_dp, 0.58691863_dp, 0.66931351_dp, 0.75191607_dp, &
    & 0.17964105_dp, 0.22157276_dp, 0.26348578_dp, 0.30539645_dp, 0.34734014_dp, &
    & 0.38924725_dp, 0.43115670_dp, 0.47308269_dp, 0.17105469_dp, 0.20276244_dp, &
    & 0.21007322_dp, 0.21739647_dp, 0.22471039_dp, 0.23201501_dp, 0.23933969_dp, &
    & 0.24665638_dp, 0.25398255_dp, 0.26128863_dp, 0.26859476_dp, 0.27592565_dp, &
    & 0.30762999_dp, 0.33931580_dp, 0.37235985_dp, 0.40273549_dp, 0.43445776_dp, &
    & 0.46611708_dp, 0.15585079_dp, 0.18649324_dp, 0.19356210_dp, 0.20063311_dp, &
    & 0.20770522_dp, 0.21477254_dp, 0.22184614_dp, 0.22891872_dp, 0.23598621_dp, &
    & 0.24305612_dp, 0.25013018_dp, 0.25719937_dp, 0.28784780_dp, 0.31848673_dp, &
    & 0.34912431_dp, 0.37976593_dp, 0.41040808_dp, 0.44105777_dp, 0.05019332_dp, &
    & 0.06762570_dp, 0.08504445_dp, 0.10247736_dp, 0.11991105_dp, 0.13732772_dp, &
    & 0.15476297_dp, 0.17218265_dp, 0.18961288_dp, 0.20704760_dp, 0.22446752_dp, &
    & 0.24189645_dp, 0.25932503_dp, 0.27676094_dp, 0.29418231_dp, 0.31159587_dp, &
    & 0.32902274_dp, 0.34592298_dp, 0.36388048_dp, 0.38130586_dp, 0.39877476_dp, &
    & 0.41614298_dp, 0.43364510_dp, 0.45104014_dp, 0.46848986_dp, 0.48584550_dp, &
    & 0.12526730_dp, 0.14268677_dp, 0.16011615_dp, 0.17755889_dp, 0.19497557_dp, &
    & 0.21240778_dp, 0.07263525_dp, 0.09422158_dp, 0.09920295_dp, 0.10418621_dp, &
    & 0.14235633_dp, 0.16394294_dp, 0.18551941_dp, 0.22370139_dp, 0.00000000_dp, &
    & 0.00000000_dp, 0.00000000_dp, 0.00000000_dp, 0.00000000_dp, 0.00000000_dp, &
    & 0.00000000_dp, 0.00000000_dp, 0.00000000_dp, 0.00000000_dp, 0.00000000_dp, &
    & 0.00000000_dp, 0.00000000_dp, 0.00000000_dp, 0.00000000_dp, 0.00000000_dp, &
    & 0.00000000_dp, 0.00000000_dp, 0.00000000_dp, 0.00000000_dp, 0.00000000_dp, &
    & 0.00000000_dp, 0.00000000_dp, 0.00000000_dp]

  !> Effective nuclear charges from the def2-ECPs used for calculating the reference
  !>  polarizibilities for DFT-D4.
  real(dp), parameter :: effectiveNuclearCharge(max_element_d4) = [ &
    &   1,                                                 2,  & ! H-He
    &   3, 4,                               5, 6, 7, 8, 9,10,  & ! Li-Ne
    &  11,12,                              13,14,15,16,17,18,  & ! Na-Ar
    &  19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,  & ! K-Kr
    &   9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,  & ! Rb-Xe
    &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Cs-Lu
    &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26, & ! Hf-Rn
    !  just copy & paste from above
    &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Fr-Lr
    &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26] ! Rf-Og

  !> Covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197), values for
  !>  metals decreased by 10%.
  real(dp), parameter :: CovalentRadiusD3(max_element_d4) = [ &
    & 0.32_dp,0.46_dp, & ! H,He
    & 1.20_dp,0.94_dp,0.77_dp,0.75_dp,0.71_dp,0.63_dp,0.64_dp,0.67_dp, & ! Li-Ne
    & 1.40_dp,1.25_dp,1.13_dp,1.04_dp,1.10_dp,1.02_dp,0.99_dp,0.96_dp, & ! Na-Ar
    & 1.76_dp,1.54_dp, & ! K,Ca
    &                 1.33_dp,1.22_dp,1.21_dp,1.10_dp,1.07_dp, & ! Sc-
    &                 1.04_dp,1.00_dp,0.99_dp,1.01_dp,1.09_dp, & ! -Zn
    &                 1.12_dp,1.09_dp,1.15_dp,1.10_dp,1.14_dp,1.17_dp, & ! Ga-Kr
    & 1.89_dp,1.67_dp, & ! Rb,Sr
    &                 1.47_dp,1.39_dp,1.32_dp,1.24_dp,1.15_dp, & ! Y-
    &                 1.13_dp,1.13_dp,1.08_dp,1.15_dp,1.23_dp, & ! -Cd
    &                 1.28_dp,1.26_dp,1.26_dp,1.23_dp,1.32_dp,1.31_dp, & ! In-Xe
    & 2.09_dp,1.76_dp, & ! Cs,Ba
    &         1.62_dp,1.47_dp,1.58_dp,1.57_dp,1.56_dp,1.55_dp,1.51_dp, & ! La-Eu
    &         1.52_dp,1.51_dp,1.50_dp,1.49_dp,1.49_dp,1.48_dp,1.53_dp, & ! Gd-Yb
    &                 1.46_dp,1.37_dp,1.31_dp,1.23_dp,1.18_dp, & ! Lu-
    &                 1.16_dp,1.11_dp,1.12_dp,1.13_dp,1.32_dp, & ! -Hg
    &                 1.30_dp,1.30_dp,1.36_dp,1.31_dp,1.38_dp,1.42_dp, & ! Tl-Rn
    & 2.01_dp,1.81_dp, & ! Fr,Ra
    &         1.67_dp,1.58_dp,1.52_dp,1.53_dp,1.54_dp,1.55_dp,1.49_dp, & ! Ac-Am
    &         1.49_dp,1.51_dp,1.51_dp,1.48_dp,1.50_dp,1.56_dp,1.58_dp, & ! Cm-No
    &                 1.45_dp,1.41_dp,1.34_dp,1.29_dp,1.27_dp, & ! Lr-
    &                 1.21_dp,1.16_dp,1.15_dp,1.09_dp,1.22_dp, & ! -Cn
    &                 1.36_dp,1.43_dp,1.46_dp,1.58_dp,1.48_dp,1.57_dp] & ! Nh-Og
    & * AA__Bohr * 4.0_dp / 3.0_dp

  !> Pauling electronegativities, used for the covalent coordination number.
  real(dp), parameter :: paulingEN(max_element_d4) = [ &
    & 2.20_dp,3.00_dp, & ! H,He
    & 0.98_dp,1.57_dp,2.04_dp,2.55_dp,3.04_dp,3.44_dp,3.98_dp,4.50_dp, & ! Li-Ne
    & 0.93_dp,1.31_dp,1.61_dp,1.90_dp,2.19_dp,2.58_dp,3.16_dp,3.50_dp, & ! Na-Ar
    & 0.82_dp,1.00_dp, & ! K,Ca
    &                 1.36_dp,1.54_dp,1.63_dp,1.66_dp,1.55_dp, & ! Sc-
    &                 1.83_dp,1.88_dp,1.91_dp,1.90_dp,1.65_dp, & ! -Zn
    &                 1.81_dp,2.01_dp,2.18_dp,2.55_dp,2.96_dp,3.00_dp, & ! Ga-Kr
    & 0.82_dp,0.95_dp, & ! Rb,Sr
    &                 1.22_dp,1.33_dp,1.60_dp,2.16_dp,1.90_dp, & ! Y-
    &                 2.20_dp,2.28_dp,2.20_dp,1.93_dp,1.69_dp, & ! -Cd
    &                 1.78_dp,1.96_dp,2.05_dp,2.10_dp,2.66_dp,2.60_dp, & ! In-Xe
    & 0.79_dp,0.89_dp, & ! Cs,Ba
    &         1.10_dp,1.12_dp,1.13_dp,1.14_dp,1.15_dp,1.17_dp,1.18_dp, & ! La-Eu
    &         1.20_dp,1.21_dp,1.22_dp,1.23_dp,1.24_dp,1.25_dp,1.26_dp, & ! Gd-Yb
    &                 1.27_dp,1.30_dp,1.50_dp,2.36_dp,1.90_dp, & ! Lu-
    &                 2.20_dp,2.20_dp,2.28_dp,2.54_dp,2.00_dp, & ! -Hg
    &                 1.62_dp,2.33_dp,2.02_dp,2.00_dp,2.20_dp,2.20_dp, & ! Tl-Rn
    ! only dummies below
    & 1.50_dp,1.50_dp, & ! Fr,Ra
    &         1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp, & ! Ac-Am
    &         1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp, & ! Cm-No
    &                 1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp, & ! Rf-
    &                 1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp, & ! Rf-Cn
    &                 1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp,1.50_dp ] ! Nh-Og


  !> PBE0/def2-QZVP atomic <r⁴>/<r²> expectation values.
  real(dp), parameter :: r4r2(max_element_d4) = [ &
    &  8.0589_dp, 3.4698_dp, & ! H,He
    & 29.0974_dp,14.8517_dp,11.8799_dp, 7.8715_dp, &
    &  5.5588_dp, 4.7566_dp, 3.8025_dp, 3.1036_dp, & ! Li-Ne
    & 26.1552_dp,17.2304_dp,17.7210_dp,12.7442_dp, &
    &  9.5361_dp, 8.1652_dp, 6.7463_dp, 5.6004_dp, & ! Na-Ar
    & 29.2012_dp,22.3934_dp, & ! K,Ca
    &            19.0598_dp,16.8590_dp,15.4023_dp,12.5589_dp,13.4788_dp, & ! Sc-
    &            12.2309_dp,11.2809_dp,10.5569_dp,10.1428_dp, 9.4907_dp, & ! -Zn
    & 13.4606_dp,10.8544_dp, 8.9386_dp, 8.1350_dp, 7.1251_dp, 6.1971_dp, & ! Ga-Kr
    & 30.0162_dp,24.4103_dp, & ! Rb,Sr
    &            20.3537_dp,17.4780_dp,13.5528_dp,11.8451_dp,11.0355_dp, & ! Y-
    &            10.1997_dp, 9.5414_dp, 9.0061_dp, 8.6417_dp, 8.9975_dp, & ! -Cd
    & 14.0834_dp,11.8333_dp,10.0179_dp, 9.3844_dp, 8.4110_dp, 7.5152_dp, & ! In-Xe
    & 32.7622_dp,27.5708_dp, & ! Cs,Ba
    & 23.1671_dp,21.6003_dp,20.9615_dp,20.4562_dp,20.1010_dp,19.7475_dp,19.4828_dp, & ! La-Eu
    & 15.6013_dp,19.2362_dp,17.4717_dp,17.8321_dp,17.4237_dp,17.1954_dp,17.1631_dp, & ! Gd-Yb
    &            14.5716_dp,15.8758_dp,13.8989_dp,12.4834_dp,11.4421_dp, & ! Lu-
    &            10.2671_dp, 8.3549_dp, 7.8496_dp, 7.3278_dp, 7.4820_dp, & ! -Hg
    & 13.5124_dp,11.6554_dp,10.0959_dp, 9.7340_dp, 8.8584_dp, 8.0125_dp, & ! Tl-Rn
    & 29.8135_dp,26.3157_dp, & ! Fr,Ra
    & 19.1885_dp,15.8542_dp,16.1305_dp,15.6161_dp,15.1226_dp,16.1576_dp, 0.0000_dp, & ! Ac-Am
    &  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, & ! Cm-No
    &             0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, & ! Lr-
    &             0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 5.4929_dp, & ! -Cn
    &  6.7286_dp, 6.5144_dp,10.9169_dp,10.3600_dp, 9.4723_dp, 8.6641_dp]   ! Nh-Og

  integer :: iDummy
  real(dp), parameter :: sqrtZr4r2(max_element_d4) = &
    &  sqrt(0.5_dp*(r4r2*[(sqrt(real(iDummy, dp)), iDummy=1, max_element_d4)]))

contains

  !> charge scaling function
  pure elemental function zetaScale(a, c, qref, qmod) result(zeta)
    real(dp),intent(in) :: qmod, qref
    real(dp),intent(in) :: a, c
    real(dp) :: zeta

    if (qmod < 0.0_dp) then
      zeta = exp(a)
    else
      zeta = exp(a * (1.0_dp - exp(c * (1.0_dp - qref/qmod))))
    end if

  end function zetaScale


  !> numerical Casimir--Polder integration using trapezoidal rule
  pure function numIntegration(pol) result(trapzd)

    !> polarizabilities at imaginary frequencies
    real(dp), intent(in) :: pol(imag_frequencies)

    !> resulting integral
    real(dp) :: trapzd

    real(dp),parameter  :: freq(imag_frequencies) = [ &
      & 0.000001_dp, 0.050000_dp, 0.100000_dp, 0.200000_dp, 0.300000_dp, &
      & 0.400000_dp, 0.500000_dp, 0.600000_dp, 0.700000_dp, 0.800000_dp, &
      & 0.900000_dp, 1.000000_dp, 1.200000_dp, 1.400000_dp, 1.600000_dp, &
      & 1.800000_dp, 2.000000_dp, 2.500000_dp, 3.000000_dp, 4.000000_dp, &
      & 5.000000_dp, 7.500000_dp, 10.00000_dp ]

    !  just precalculate all weights and get the job done
    real(dp),parameter :: weights(imag_frequencies) = 0.5_dp * [ &
      & (freq (2) - freq (1)),  &
      & (freq (2) - freq (1)) + (freq (3) - freq (2)), &
      & (freq (3) - freq (2)) + (freq (4) - freq (3)), &
      & (freq (4) - freq (3)) + (freq (5) - freq (4)), &
      & (freq (5) - freq (4)) + (freq (6) - freq (5)), &
      & (freq (6) - freq (5)) + (freq (7) - freq (6)), &
      & (freq (7) - freq (6)) + (freq (8) - freq (7)), &
      & (freq (8) - freq (7)) + (freq (9) - freq (8)), &
      & (freq (9) - freq (8)) + (freq(10) - freq (9)), &
      & (freq(10) - freq (9)) + (freq(11) - freq(10)), &
      & (freq(11) - freq(10)) + (freq(12) - freq(11)), &
      & (freq(12) - freq(11)) + (freq(13) - freq(12)), &
      & (freq(13) - freq(12)) + (freq(14) - freq(13)), &
      & (freq(14) - freq(13)) + (freq(15) - freq(14)), &
      & (freq(15) - freq(14)) + (freq(16) - freq(15)), &
      & (freq(16) - freq(15)) + (freq(17) - freq(16)), &
      & (freq(17) - freq(16)) + (freq(18) - freq(17)), &
      & (freq(18) - freq(17)) + (freq(19) - freq(18)), &
      & (freq(19) - freq(18)) + (freq(20) - freq(19)), &
      & (freq(20) - freq(19)) + (freq(21) - freq(20)), &
      & (freq(21) - freq(20)) + (freq(22) - freq(21)), &
      & (freq(22) - freq(21)) + (freq(23) - freq(22)), &
      & (freq(23) - freq(22)) ]

    trapzd = sum(pol*weights)

  end function numIntegration


  subroutine initializeCalculator(calculator, input, nAtom, izp)

    !> Calculator
    type(DftD4Calculator), intent(inout) :: calculator

    !> Input
    type(DispDftD4Inp), intent(in) :: input

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Atomic number of every atom.
    integer, intent(in) :: izp(:)

    integer :: mAt
    integer :: iAt1, iZp1, iSec, iCN, iZp2, iRef1, iRef2
    integer :: cncount(0:18)
    real(dp) :: alpha(imag_frequencies), zEff1, c6, eta1
    real(dp) :: tmp_hq(max_references, max_element_d4)

    real(dp), parameter :: thopi = 3.0_dp/pi

    ! initialize charge model
    calculator%chi = chi
    calculator%gam = gam
    calculator%rad = rad
    calculator%kcn = kcn

    calculator%sqrtZr4r2 = sqrtZr4r2
    calculator%covalentRadius = CovalentRadiusD3
    calculator%electronegativity = paulingEN
    calculator%ChemicalHardness = ChemicalHardness
    calculator%EffectiveNuclearCharge = EffectiveNuclearCharge

    calculator%s6 = input%s6
    calculator%s8 = input%s8
    calculator%s10 = input%s10
    calculator%s9 = input%s9
    calculator%a1 = input%a1
    calculator%a2 = input%a2
    calculator%alpha = input%alpha

    calculator%nrChrg = input%nrChrg

    calculator%cutoffCount = input%cutoffCount
    calculator%cutoffInter = input%cutoffInter
    calculator%cutoffThree = input%cutoffThree
    calculator%cutoffEwald = input%cutoffEwald

    calculator%wf = input%weightingFactor
    calculator%ga = input%chargeScale
    calculator%gc = input%chargeSteepness

    mAt = maxval(izp)
    allocate(calculator%numberOfReferences(mAt), calculator%atoms(mAt), &
        & calculator%countNumber(max_references, mAt))
    calculator%numberOfReferences(:) = 0
    calculator%atoms(:) = 0
    calculator%countNumber(:, :) = 0
    allocate(calculator%referenceCN(max_references, mAt), &
        & calculator%referenceCharge(max_references, mAt), &
        & calculator%referenceAlpha(imag_frequencies, max_references, mAt), &
        & calculator%referenceC6(max_references, max_references, mAt, mAt))
    calculator%referenceCN(:, :) = 0.0_dp
    calculator%referenceCharge(:, :) = 0.0_dp
    calculator%referenceAlpha(:, :, :) = 0.0_dp
    calculator%referenceC6(:, :, :, :) = 0.0_dp

    calculator%referenceCharge(:,:) = clsq(:, :mAt)
    tmp_hq(:,:) = clsh

    ! evaluate α(Z) = 1/k·(α(ZkBn) - ζ(n, m) · n/m · α(Bm))
    ! α(Z) is referenceAlpha, α(ZkBn) is alphaiw, α(Bm) is secaiw
    ! 1/m is sscale, 1/k is ascale, n is hcount, ζ(n, m) is zetaScale
    do iAt1 = 1, nAtom
      cncount(:) = 0
      cncount(0) = 1
      iZp1 = izp(iAt1)
      if (calculator%atoms(iZp1) == 0) then
        calculator%numberOfReferences(iZp1) = refn(iZp1)
        do iRef1 = 1, calculator%numberOfReferences(iZp1)
          iSec = refsys(iRef1,iZp1)
          eta1 = calculator%gc * calculator%chemicalHardness(iSec)
          zEff1 = calculator%effectiveNuclearCharge(iSec)
          alpha = sscale(iSec) * secaiw(:,iSec) &
              & * zetaScale(calculator%ga, eta1, zEff1, tmp_hq(iRef1,iZp1) + zEff1)
          iCN = nint(refcn(iRef1,iZp1))
          calculator%referenceCN(iRef1,iZp1) = refcovcn(iRef1,iZp1)
          cncount(iCN) = cncount(iCN) + 1
          calculator%referenceAlpha(:,iRef1,iZp1) = &
              & max(0.0_dp, ascale(iRef1,iZp1) * (alphaiw(:,iRef1,iZp1) - hcount(iRef1,iZp1)*alpha))
        end do
        ! setup the number of Gaussian functions for the weighting in countNumber
        do iRef1 = 1, calculator%numberOfReferences(iZp1)
          iCN = cncount(nint(refcn(iRef1,iZp1)))
          calculator%countNumber(iRef1,iZp1) = iCN * (iCN + 1) / 2
        end do
      end if
      calculator%atoms(iZp1) = calculator%atoms(iZp1) + 1
    end do

    ! integrate C6 coefficients
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(iZp1, iZp2, iRef1, iRef2, alpha, c6) SHARED(calculator, mAt)
    !$OMP DO SCHEDULE(RUNTIME)
    do iZp1 = 1, mAt
      do iZp2 = 1, iZp1
        if (calculator%atoms(iZp1) > 0 .and. calculator%atoms(iZp2) > 0) then
          do iRef1 = 1, calculator%numberOfReferences(iZp1)
            do iRef2 = 1, calculator%numberOfReferences(iZp2)
              alpha = calculator%referenceAlpha(:,iRef1,iZp1)&
                  & * calculator%referenceAlpha(:,iRef2,iZp2)
              c6 = thopi * numIntegration(alpha)
              calculator%referenceC6(iRef2,iRef1,iZp2,iZp1) = c6
              calculator%referenceC6(iRef1,iRef2,iZp1,iZp2) = c6
            end do
          end do
        end if
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine initializeCalculator


end module dftbp_dftd4param
