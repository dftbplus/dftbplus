!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module dftbp_gfn1param
  use dftbp_accuracy, only : dp
  use dftbp_constants, only : eV => eV__Hartree, AA => AA__Bohr, symbolToNumber
  use dftbp_xtbparam, only : xtbParam, xtbBasis, xtbCalculator, xtbGlobalParameter
  implicit none
  private

  public :: getGFN1Param, gfn1Globals


  interface getGFN1Param
    module procedure :: getGFN1ParamNumber
    module procedure :: getGFN1ParamSymbol
  end interface getGFN1Param


  !> Maximum number of elements supported by GFN1-xTB
  integer, parameter :: gfn1Elem = 86


  type(xtbGlobalParameter), parameter :: gfn1Globals = xtbGlobalParameter(&
      & kEnScale = -7e-3_dp, &
      & krep = 1.5_dp, &
      & xbondRad = 1.3_dp, &
      & xbondExp = 0.44_dp)


  !> Atomic parameters for GFN1-xTB
  type(xtbParam), parameter :: gfn1Param(1:gfn1Elem) = [ &
    !                -EN----  -radius---  -eta-------  -gam-------  -alpha-----   -zeff------  -xbond-
    & xtbParam(2, 0, 2.20_dp, 0.32_dp*AA, 0.470099_dp, 0.000000_dp, 2.209700_dp,  1.116244_dp,-1.00_dp), & ! H
    & xtbParam(1, 0, 3.00_dp, 0.37_dp*AA, 1.441379_dp, 1.500000_dp, 1.382907_dp,  0.440231_dp, 0.00_dp), & ! He
    & xtbParam(2, 1, 0.98_dp, 1.30_dp*AA, 0.205342_dp, 1.027370_dp, 0.671797_dp,  2.747587_dp, 0.00_dp), & ! Li
    & xtbParam(2, 1, 1.57_dp, 0.99_dp*AA, 0.274022_dp, 0.900554_dp, 0.865377_dp,  4.076830_dp, 0.00_dp), & ! Be
    & xtbParam(2, 1, 2.04_dp, 0.84_dp*AA, 0.340530_dp, 1.300000_dp, 1.093544_dp,  4.458376_dp, 0.00_dp), & ! B
    & xtbParam(2, 0, 2.55_dp, 0.75_dp*AA, 0.479988_dp, 1.053856_dp, 1.281954_dp,  4.428763_dp, 0.00_dp), & ! C
    & xtbParam(2, 0, 3.04_dp, 0.71_dp*AA, 0.476106_dp, 0.042507_dp, 1.727773_dp,  5.498808_dp, 0.00_dp), & ! N
    & xtbParam(2, 0, 3.44_dp, 0.64_dp*AA, 0.583349_dp,-0.005102_dp, 2.004253_dp,  5.171786_dp, 0.00_dp), & ! O
    & xtbParam(2, 0, 3.98_dp, 0.60_dp*AA, 0.788194_dp, 1.615037_dp, 2.507078_dp,  6.931741_dp, 0.00_dp), & ! F
    & xtbParam(3, 0, 4.50_dp, 0.62_dp*AA, 0.612878_dp, 1.600000_dp, 3.038727_dp,  9.102523_dp, 0.00_dp), & ! Ne
    & xtbParam(2, 1, 0.93_dp, 1.60_dp*AA, 0.165908_dp, 1.200000_dp, 0.704472_dp, 10.591259_dp, 0.00_dp), & ! Na
    & xtbParam(2, 1, 1.31_dp, 1.40_dp*AA, 0.354151_dp, 1.100000_dp, 0.862629_dp, 15.238107_dp, 0.00_dp), & ! Mg
    & xtbParam(3, 1, 1.61_dp, 1.24_dp*AA, 0.221658_dp, 1.200000_dp, 0.929219_dp, 16.283595_dp, 0.00_dp), & ! Al
    & xtbParam(3, 0, 1.90_dp, 1.14_dp*AA, 0.438331_dp, 1.500000_dp, 0.948165_dp, 16.898359_dp, 0.00_dp), & ! Si
    & xtbParam(3, 0, 2.19_dp, 1.09_dp*AA, 0.798319_dp, 1.500000_dp, 1.067197_dp, 15.249559_dp, 0.00_dp), & ! P
    & xtbParam(3, 0, 2.58_dp, 1.04_dp*AA, 0.643959_dp, 1.500000_dp, 1.200803_dp, 15.100323_dp, 0.00_dp), & ! S
    & xtbParam(3, 0, 3.16_dp, 1.00_dp*AA, 0.519712_dp, 1.000000_dp, 1.404155_dp, 17.000000_dp, 0.00_dp), & ! Cl
    & xtbParam(3, 0, 3.50_dp, 1.01_dp*AA, 0.529906_dp, 0.829312_dp, 1.323756_dp, 17.153132_dp, 0.00_dp), & ! Ar
    & xtbParam(2, 1, 0.82_dp, 2.00_dp*AA, 0.114358_dp, 0.732923_dp, 0.581529_dp, 20.831436_dp, 0.00_dp), & ! K
    & xtbParam(3, 1, 1.00_dp, 1.74_dp*AA, 0.134187_dp, 1.116963_dp, 0.665588_dp, 19.840212_dp, 0.00_dp), & ! Ca
    & xtbParam(3,-1, 1.36_dp, 1.59_dp*AA, 0.778545_dp, 1.000000_dp, 0.841357_dp, 18.676202_dp, 0.00_dp), & ! Sc
    & xtbParam(3,-1, 1.54_dp, 1.48_dp*AA, 1.044998_dp, 0.739203_dp, 0.828638_dp, 17.084130_dp, 0.00_dp), & ! Ti
    & xtbParam(3,-1, 1.63_dp, 1.44_dp*AA, 0.985157_dp, 0.800000_dp, 1.061627_dp, 22.352532_dp, 0.00_dp), & ! V
    & xtbParam(3,-1, 1.66_dp, 1.30_dp*AA, 0.468100_dp, 0.800000_dp, 0.997051_dp, 22.873486_dp, 0.00_dp), & ! Cr
    & xtbParam(3, 1, 1.55_dp, 1.29_dp*AA, 0.609868_dp, 0.300000_dp, 1.019783_dp, 24.160655_dp, 0.00_dp), & ! Mn
    & xtbParam(3, 1, 1.83_dp, 1.24_dp*AA, 0.900000_dp, 0.500000_dp, 1.137174_dp, 25.983149_dp, 0.00_dp), & ! Fe
    & xtbParam(3, 1, 1.88_dp, 1.18_dp*AA, 0.426680_dp, 0.300000_dp, 1.188538_dp, 27.169215_dp, 0.00_dp), & ! Co
    & xtbParam(3, 1, 1.91_dp, 1.17_dp*AA, 0.367019_dp, 1.000000_dp, 1.399197_dp, 23.396999_dp, 0.00_dp), & ! Ni
    & xtbParam(3, 1, 1.90_dp, 1.22_dp*AA, 0.260192_dp, 0.237602_dp, 1.199230_dp, 29.000000_dp, 0.00_dp), & ! Cu
    & xtbParam(2, 1, 1.65_dp, 1.20_dp*AA, 0.209459_dp, 1.400000_dp, 1.145056_dp, 31.185765_dp, 0.00_dp), & ! Zn
    & xtbParam(3, 1, 1.81_dp, 1.23_dp*AA, 0.193302_dp, 1.400000_dp, 1.047536_dp, 33.128619_dp, 0.00_dp), & ! Ga
    & xtbParam(3, 0, 2.01_dp, 1.20_dp*AA, 0.800000_dp, 1.400000_dp, 1.129480_dp, 35.493164_dp, 0.00_dp), & ! Ge
    & xtbParam(3, 0, 2.18_dp, 1.20_dp*AA, 0.732367_dp, 1.300000_dp, 1.233641_dp, 36.125762_dp, 0.00_dp), & ! As
    & xtbParam(3, 0, 2.55_dp, 1.18_dp*AA, 0.714534_dp, 1.300000_dp, 1.270088_dp, 32.148852_dp, 0.00_dp), & ! Se
    & xtbParam(3, 0, 2.96_dp, 1.17_dp*AA, 0.732530_dp,-0.500000_dp, 1.153580_dp, 35.000000_dp, 0.38_dp), & ! Br
    & xtbParam(3, 0, 3.00_dp, 1.16_dp*AA, 0.820312_dp, 1.000000_dp, 1.335287_dp, 36.000000_dp, 0.00_dp), & ! Kr
    & xtbParam(2, 1, 0.82_dp, 2.15_dp*AA, 0.075735_dp, 1.500000_dp, 0.554032_dp, 39.653032_dp, 0.00_dp), & ! Rb
    & xtbParam(3, 1, 0.95_dp, 1.90_dp*AA, 0.122861_dp, 1.300000_dp, 0.657904_dp, 38.924904_dp, 0.00_dp), & ! Sr
    & xtbParam(3, 1, 1.22_dp, 1.76_dp*AA, 0.351290_dp, 1.400000_dp, 0.760144_dp, 39.000000_dp, 0.00_dp), & ! Y
    & xtbParam(3, 1, 1.33_dp, 1.64_dp*AA, 0.168219_dp, 0.581478_dp, 0.739520_dp, 36.521516_dp, 0.00_dp), & ! Zr
    & xtbParam(3, 1, 1.60_dp, 1.56_dp*AA, 0.175875_dp, 0.280147_dp, 0.895357_dp, 40.803132_dp, 0.00_dp), & ! Nb
    & xtbParam(3, 1, 2.16_dp, 1.46_dp*AA, 0.384677_dp, 0.041052_dp, 0.944064_dp, 41.939347_dp, 0.00_dp), & ! Mo
    & xtbParam(3, 1, 1.90_dp, 1.38_dp*AA, 0.405474_dp, 0.500000_dp, 1.028240_dp, 43.000000_dp, 0.00_dp), & ! Tc
    & xtbParam(3, 1, 2.20_dp, 1.36_dp*AA, 0.305394_dp, 0.001205_dp, 1.066144_dp, 44.492732_dp, 0.00_dp), & ! Ru
    & xtbParam(3, 1, 2.28_dp, 1.34_dp*AA, 0.293973_dp, 0.622690_dp, 1.131380_dp, 45.241537_dp, 0.00_dp), & ! Rh
    & xtbParam(3, 1, 2.20_dp, 1.30_dp*AA, 0.280766_dp, 0.500000_dp, 1.206869_dp, 42.105527_dp, 0.00_dp), & ! Pd
    & xtbParam(3, 1, 1.93_dp, 1.36_dp*AA, 0.472978_dp,-0.445675_dp, 1.058886_dp, 43.201446_dp, 0.00_dp), & ! Ag
    & xtbParam(2, 1, 1.69_dp, 1.40_dp*AA, 0.130828_dp, 1.362587_dp, 1.026434_dp, 49.016827_dp, 0.00_dp), & ! Cd
    & xtbParam(3, 1, 1.78_dp, 1.42_dp*AA, 0.132120_dp, 1.063557_dp, 0.898148_dp, 51.718417_dp, 0.00_dp), & ! In
    & xtbParam(3, 0, 1.96_dp, 1.40_dp*AA, 0.480655_dp,-0.321283_dp, 1.008192_dp, 54.503455_dp, 0.00_dp), & ! Sn
    & xtbParam(3, 0, 2.05_dp, 1.40_dp*AA, 0.564406_dp,-0.341503_dp, 0.982673_dp, 50.757213_dp, 0.00_dp), & ! Sb
    & xtbParam(3, 0, 2.10_dp, 1.37_dp*AA, 0.400301_dp, 0.894388_dp, 0.973410_dp, 49.215262_dp, 0.00_dp), & ! Te
    & xtbParam(3, 0, 2.66_dp, 1.36_dp*AA, 0.520472_dp,-0.500000_dp, 0.949181_dp, 53.000000_dp, 0.32_dp), & ! I
    & xtbParam(3, 0, 2.60_dp, 1.36_dp*AA, 0.935394_dp,-0.800000_dp, 1.074785_dp, 52.500985_dp, 0.00_dp), & ! Xe
    & xtbParam(2, 1, 0.79_dp, 2.38_dp*AA, 0.085110_dp, 1.500000_dp, 0.579919_dp, 65.029838_dp, 0.00_dp), & ! Cs
    & xtbParam(3, 1, 0.89_dp, 2.06_dp*AA, 0.137819_dp, 1.500000_dp, 0.606485_dp, 46.532974_dp, 0.00_dp), & ! Ba
    & xtbParam(3, 1, 1.10_dp, 1.94_dp*AA, 0.495969_dp, 1.500000_dp, 1.311200_dp, 48.337542_dp, 0.00_dp), & ! La
    & xtbParam(3, 1, 1.12_dp, 1.84_dp*AA, 0.350000_dp, 1.200000_dp, 0.839861_dp, 30.638143_dp, 0.00_dp), & ! Ce
    & xtbParam(3, 1, 1.13_dp, 1.90_dp*AA, 0.342306_dp, 1.200000_dp, 0.847281_dp, 34.130718_dp, 0.00_dp), & ! Pr
    & xtbParam(3, 1, 1.14_dp, 1.88_dp*AA, 0.334612_dp, 1.200000_dp, 0.854701_dp, 37.623294_dp, 0.00_dp), & ! Nd
    & xtbParam(3, 1, 1.13_dp, 1.86_dp*AA, 0.326917_dp, 1.200000_dp, 0.862121_dp, 41.115870_dp, 0.00_dp), & ! Pm
    & xtbParam(3, 1, 1.17_dp, 1.85_dp*AA, 0.319223_dp, 1.200000_dp, 0.869541_dp, 44.608445_dp, 0.00_dp), & ! Sm
    & xtbParam(3, 1, 1.20_dp, 1.83_dp*AA, 0.311529_dp, 1.200000_dp, 0.876961_dp, 48.101021_dp, 0.00_dp), & ! Eu
    & xtbParam(3, 1, 1.20_dp, 1.82_dp*AA, 0.303835_dp, 1.200000_dp, 0.884381_dp, 51.593596_dp, 0.00_dp), & ! Gd
    & xtbParam(3, 1, 1.22_dp, 1.81_dp*AA, 0.296140_dp, 1.200000_dp, 0.891801_dp, 55.086172_dp, 0.00_dp), & ! Tb
    & xtbParam(3, 1, 1.23_dp, 1.80_dp*AA, 0.288446_dp, 1.200000_dp, 0.899221_dp, 58.578748_dp, 0.00_dp), & ! Dy
    & xtbParam(3, 1, 1.24_dp, 1.79_dp*AA, 0.280752_dp, 1.200000_dp, 0.906641_dp, 62.071323_dp, 0.00_dp), & ! Ho
    & xtbParam(3, 1, 1.24_dp, 1.77_dp*AA, 0.273058_dp, 1.200000_dp, 0.914061_dp, 65.563899_dp, 0.00_dp), & ! Er
    & xtbParam(3, 1, 1.25_dp, 1.77_dp*AA, 0.265364_dp, 1.200000_dp, 0.921481_dp, 69.056474_dp, 0.00_dp), & ! Tm
    & xtbParam(3, 1, 1.10_dp, 1.78_dp*AA, 0.257669_dp, 1.200000_dp, 0.928901_dp, 72.549050_dp, 0.00_dp), & ! Yb
    & xtbParam(3, 1, 1.27_dp, 1.74_dp*AA, 0.249975_dp, 1.200000_dp, 0.936321_dp, 76.041625_dp, 0.00_dp), & ! Lu
    & xtbParam(3, 1, 1.30_dp, 1.64_dp*AA, 0.269977_dp, 0.847011_dp, 0.853744_dp, 55.222897_dp, 0.00_dp), & ! Hf
    & xtbParam(3, 1, 1.50_dp, 1.58_dp*AA, 0.239696_dp, 0.064592_dp, 0.971873_dp, 63.743065_dp, 0.00_dp), & ! Ta
    & xtbParam(3, 1, 2.36_dp, 1.50_dp*AA, 0.243663_dp,-0.014599_dp, 0.992643_dp, 74.000000_dp, 0.00_dp), & ! W
    & xtbParam(3, 1, 1.90_dp, 1.41_dp*AA, 0.362512_dp, 0.300000_dp, 1.132106_dp, 75.000000_dp, 0.00_dp), & ! Re
    & xtbParam(3, 1, 2.20_dp, 1.36_dp*AA, 0.354318_dp,-0.170295_dp, 1.118216_dp, 76.000000_dp, 0.00_dp), & ! Os
    & xtbParam(3, 1, 2.20_dp, 1.32_dp*AA, 0.290898_dp, 0.965726_dp, 1.245003_dp, 77.000000_dp, 0.00_dp), & ! Ir
    & xtbParam(3, 1, 2.28_dp, 1.30_dp*AA, 0.370447_dp, 1.092759_dp, 1.304590_dp, 78.000000_dp, 0.00_dp), & ! Pt
    & xtbParam(3, 1, 2.54_dp, 1.30_dp*AA, 0.496380_dp, 0.123512_dp, 1.293034_dp, 79.000000_dp, 0.00_dp), & ! Au
    & xtbParam(2, 1, 2.00_dp, 1.32_dp*AA, 0.334997_dp,-0.267745_dp, 1.181865_dp, 80.000000_dp, 0.00_dp), & ! Hg
    & xtbParam(2, 1, 1.62_dp, 1.44_dp*AA, 0.671316_dp, 0.936157_dp, 0.976397_dp, 81.000000_dp, 0.00_dp), & ! Tl
    & xtbParam(2, 0, 2.33_dp, 1.45_dp*AA, 1.000000_dp, 1.500000_dp, 0.988859_dp, 79.578302_dp, 0.00_dp), & ! Pb
    & xtbParam(2, 0, 2.02_dp, 1.50_dp*AA, 0.944879_dp, 0.877488_dp, 1.047194_dp, 83.000000_dp, 0.00_dp), & ! Bi
    & xtbParam(3, 0, 2.00_dp, 1.42_dp*AA, 1.091248_dp,-0.035874_dp, 1.013118_dp, 84.000000_dp, 0.00_dp), & ! Po
    & xtbParam(3, 0, 2.20_dp, 1.48_dp*AA, 1.264162_dp,-0.860502_dp, 0.964652_dp, 85.000000_dp, 0.22_dp), & ! At
    & xtbParam(3, 0, 2.20_dp, 1.46_dp*AA, 0.798170_dp,-0.838429_dp, 0.998641_dp, 86.000000_dp, 0.00_dp)]   ! Rn


  !> Slater basisset and shell related parameters for GFN1-xTB
  type(xtbBasis), parameter :: gfn1Basis(3, 1:gfn1Elem) = reshape([ &
    !                  -poly------  -leta------  -h-------------  -kcn---  -zeta------
    & xtbBasis( 1, 0,  0.000000_dp, 0.000000_dp,-10.923452_dp*eV, 6e-3_dp, 1.207940_dp, 4, 1), & ! H  1s
    & xtbBasis( 2, 0,  0.000000_dp, 0.000000_dp, -2.171902_dp*eV, 6e-3_dp, 1.993207_dp, 3, 0), & ! H  2s
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 1, 0,  8.084149_dp, 0.000000_dp,-22.121015_dp*eV, 6e-3_dp, 2.133698_dp, 4, 1), & ! He 1s
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 2, 0, -4.102845_dp, 0.000000_dp, -7.270105_dp*eV, 0e+0_dp, 0.743881_dp, 6, 1), & ! Li 2s
    & xtbBasis( 2, 1,  9.259276_dp,-0.772012_dp, -4.609277_dp*eV, 0e+0_dp, 0.541917_dp, 6, 1), & ! Li 2p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 2, 0,-12.991482_dp, 0.000000_dp, -9.541494_dp*eV, 0e+0_dp, 0.876888_dp, 6, 1), & ! Be 2s
    & xtbBasis( 2, 1, -1.308797_dp, 1.113005_dp, -5.812621_dp*eV, 0e+0_dp, 1.104598_dp, 6, 1), & ! Be 2p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 2, 0, -7.088823_dp, 0.000000_dp,-12.497913_dp*eV, 0e+0_dp, 1.667617_dp, 6, 1), & ! B  2s
    & xtbBasis( 2, 1,  0.655877_dp, 0.165643_dp, -7.604923_dp*eV, 0e+0_dp, 1.495078_dp, 6, 1), & ! B  2p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 2, 0, -7.082170_dp, 0.000000_dp,-13.587210_dp*eV, 6e-3_dp, 1.960324_dp, 6, 1), & ! C  2s
    & xtbBasis( 2, 1,  0.812216_dp,-0.471181_dp,-10.052785_dp*eV,-3e-3_dp, 1.832096_dp, 6, 1), & ! C  2p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 2, 0,-12.745585_dp, 0.000000_dp,-20.058000_dp*eV, 6e-3_dp, 2.050067_dp, 6, 1), & ! N  2s
    & xtbBasis( 2, 1, -1.428367_dp, 0.315090_dp,-12.889326_dp*eV,-3e-3_dp, 2.113682_dp, 6, 1), & ! N  2p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 2, 0,-13.729047_dp, 0.000000_dp,-23.398376_dp*eV, 6e-3_dp, 2.345365_dp, 6, 1), & ! O  2s
    & xtbBasis( 2, 1, -4.453341_dp, 0.374608_dp,-17.886554_dp*eV,-3e-3_dp, 2.153060_dp, 6, 1), & ! O  2p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 2, 0, -3.921613_dp, 0.000000_dp,-24.776175_dp*eV, 6e-3_dp, 2.968015_dp, 6, 1), & ! F  2s
    & xtbBasis( 2, 1,-11.422491_dp,-0.827352_dp,-17.274415_dp*eV,-3e-3_dp, 2.256959_dp, 6, 1), & ! F  2p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 2, 0, -2.115896_dp, 0.000000_dp,-31.167487_dp*eV, 6e-3_dp, 3.200000_dp, 6, 1), & ! Ne 2s
    & xtbBasis( 2, 1,-15.124326_dp,-3.892542_dp,-18.268975_dp*eV,-3e-3_dp, 2.294365_dp, 6, 1), & ! Ne 2p
    & xtbBasis( 3, 2,  0.000000_dp, 0.000000_dp,  1.487984_dp*eV,-5e-3_dp, 2.684436_dp, 4, 1), & ! Ne 3d
    & xtbBasis( 3, 0, 13.188489_dp, 0.000000_dp, -4.717569_dp*eV, 0e+0_dp, 0.819143_dp, 6, 1), & ! Na 3s
    & xtbBasis( 3, 1, 10.969376_dp,-3.004391_dp, -2.797054_dp*eV, 0e+0_dp, 0.628961_dp, 6, 1), & ! Na 3p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 3, 0,-19.219408_dp, 0.000000_dp, -9.970921_dp*eV, 0e+0_dp, 1.271287_dp, 6, 1), & ! Mg 3s
    & xtbBasis( 3, 1, 18.272922_dp, 0.674819_dp, -2.901013_dp*eV, 0e+0_dp, 0.797143_dp, 6, 1), & ! Mg 3p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 3, 0,-21.085827_dp, 0.000000_dp,-12.916245_dp*eV, 0e+0_dp, 1.497753_dp, 6, 1), & ! Al 3s
    & xtbBasis( 3, 1, 24.805127_dp, 0.503564_dp, -3.441043_dp*eV, 0e+0_dp, 1.232966_dp, 6, 1), & ! Al 3p
    & xtbBasis( 3, 2, 26.405814_dp, 0.000000_dp, -1.751415_dp*eV, 0e+0_dp, 0.606937_dp, 4, 1), & ! Al 3d
    & xtbBasis( 3, 0,-14.201582_dp, 0.000000_dp,-14.506128_dp*eV, 6e-3_dp, 1.521960_dp, 6, 1), & ! Si 3s
    & xtbBasis( 3, 1, -3.893343_dp,-5.925834_dp, -7.557337_dp*eV,-3e-3_dp, 1.609138_dp, 6, 1), & ! Si 3p
    & xtbBasis( 3, 2, 25.499221_dp, 0.000000_dp, -2.508113_dp*eV,-5e-3_dp, 1.168971_dp, 4, 1), & ! Si 3d
    & xtbBasis( 3, 0,-16.118985_dp, 0.000000_dp,-18.865587_dp*eV, 6e-3_dp, 1.993165_dp, 6, 1), & ! P  3s
    & xtbBasis( 3, 1, -2.241189_dp,-2.530875_dp, -9.386464_dp*eV,-3e-3_dp, 1.826973_dp, 6, 1), & ! P  3p
    & xtbBasis( 3, 2, 30.984577_dp, 0.000000_dp, -0.673989_dp*eV,-5e-3_dp, 1.293345_dp, 4, 1), & ! P  3d
    & xtbBasis( 3, 0,-16.989922_dp, 0.000000_dp,-23.819013_dp*eV, 6e-3_dp, 2.506934_dp, 6, 1), & ! S  3s
    & xtbBasis( 3, 1, -6.067779_dp,-1.678147_dp,-12.120136_dp*eV,-3e-3_dp, 1.992775_dp, 6, 1), & ! S  3p
    & xtbBasis( 3, 2, 16.248395_dp, 0.000000_dp, -1.711261_dp*eV,-5e-3_dp, 1.964867_dp, 4, 1), & ! S  3d
    & xtbBasis( 3, 0, -9.341919_dp, 0.000000_dp,-24.452163_dp*eV, 6e-3_dp, 2.847946_dp, 6, 1), & ! Cl 3s
    & xtbBasis( 3, 1, -8.499805_dp,-4.481841_dp,-12.883714_dp*eV,-3e-3_dp, 2.077534_dp, 6, 1), & ! Cl 3p
    & xtbBasis( 3, 2, 13.088867_dp, 0.000000_dp, -1.190095_dp*eV,-5e-3_dp, 1.932463_dp, 4, 1), & ! Cl 3d
    & xtbBasis( 3, 0, -0.082808_dp, 0.000000_dp,-31.395427_dp*eV, 6e-3_dp, 3.502323_dp, 6, 1), & ! Ar 3s
    & xtbBasis( 3, 1, -9.217948_dp,-1.450000_dp,-17.412901_dp*eV,-3e-3_dp, 2.287983_dp, 6, 1), & ! Ar 3p
    & xtbBasis( 3, 2, 12.204172_dp, 0.000000_dp, -1.119399_dp*eV,-5e-3_dp, 1.761181_dp, 4, 1), & ! Ar 3d
    & xtbBasis( 4, 0, 12.482844_dp, 0.000000_dp, -5.815562_dp*eV, 0e+0_dp, 0.841791_dp, 6, 1), & ! K  4s
    & xtbBasis( 4, 1, 22.323655_dp,-5.332978_dp, -3.747255_dp*eV, 0e+0_dp, 0.771618_dp, 6, 1), & ! K  4p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 4, 0,-11.421376_dp, 0.000000_dp, -7.979180_dp*eV, 0e+0_dp, 1.321845_dp, 6, 1), & ! Ca 4s
    & xtbBasis( 4, 1, 14.628284_dp,11.522018_dp, -2.517008_dp*eV, 0e+0_dp, 0.734954_dp, 6, 1), & ! Ca 4p
    & xtbBasis( 3, 2, 10.129602_dp, 0.000000_dp, -2.752355_dp*eV, 0e+0_dp, 0.947032_dp, 4, 1), & ! Ca 3d
    & xtbBasis( 3, 2,-36.027863_dp,-5.934820_dp, -7.172021_dp*eV, 5e-3_dp, 2.200000_dp, 4, 1), & ! Sc 3d
    & xtbBasis( 4, 0,  9.522966_dp, 0.000000_dp, -9.632943_dp*eV, 6e-3_dp, 1.532191_dp, 6, 1), & ! Sc 4s
    & xtbBasis( 4, 1, 44.183320_dp,-2.000000_dp, -0.696628_dp*eV,-3e-3_dp, 1.017366_dp, 6, 1), & ! Sc 4p
    & xtbBasis( 3, 2,-24.908650_dp,-7.388986_dp, -7.617343_dp*eV, 5e-3_dp, 1.941479_dp, 4, 1), & ! Ti 3d
    & xtbBasis( 4, 0, 24.879987_dp, 0.000000_dp, -7.948161_dp*eV, 6e-3_dp, 1.477526_dp, 6, 1), & ! Ti 4s
    & xtbBasis( 4, 1, 18.910954_dp,-1.500000_dp, -0.902143_dp*eV,-3e-3_dp, 1.063921_dp, 6, 1), & ! Ti 4p
    & xtbBasis( 3, 2,-29.197847_dp,-5.229338_dp, -6.677563_dp*eV, 5e-3_dp, 1.812440_dp, 4, 1), & ! V  3d
    & xtbBasis( 4, 0, -5.301066_dp, 0.000000_dp, -9.000000_dp*eV, 6e-3_dp, 1.345487_dp, 6, 1), & ! V  4s
    & xtbBasis( 4, 1, 22.945047_dp,-2.000000_dp, -0.108008_dp*eV,-3e-3_dp, 1.100000_dp, 6, 1), & ! V  4p
    & xtbBasis( 3, 2,-22.608167_dp, 0.786859_dp, -7.357172_dp*eV, 5e-3_dp, 1.915482_dp, 4, 1), & ! Cr 3d
    & xtbBasis( 4, 0, -2.432193_dp, 0.000000_dp, -7.024438_dp*eV, 6e-3_dp, 1.241910_dp, 6, 1), & ! Cr 4s
    & xtbBasis( 4, 1, 11.274054_dp,-2.500000_dp, -3.933133_dp*eV,-3e-3_dp, 1.130000_dp, 6, 1), & ! Cr 4p
    & xtbBasis( 3, 2,-25.016650_dp,10.544199_dp, -8.558648_dp*eV, 0e+0_dp, 2.016302_dp, 4, 1), & ! Mn 3d
    & xtbBasis( 4, 0,  1.025345_dp, 0.000000_dp, -6.149482_dp*eV, 0e+0_dp, 1.882798_dp, 6, 1), & ! Mn 4s
    & xtbBasis( 4, 1,  1.834626_dp,-2.500000_dp, -4.360801_dp*eV, 0e+0_dp, 1.270000_dp, 6, 1), & ! Mn 4p
    & xtbBasis( 3, 2,-22.920815_dp, 1.018896_dp, -9.705009_dp*eV, 0e+0_dp, 2.264485_dp, 4, 1), & ! Fe 3d
    & xtbBasis( 4, 0, -2.182723_dp, 0.000000_dp, -6.617863_dp*eV, 0e+0_dp, 1.382959_dp, 6, 1), & ! Fe 4s
    & xtbBasis( 4, 1, 11.769535_dp,-2.000000_dp, -4.595985_dp*eV, 0e+0_dp, 1.300000_dp, 6, 1), & ! Fe 4p
    & xtbBasis( 3, 2,-21.678930_dp, 0.222849_dp,-10.285239_dp*eV, 0e+0_dp, 2.279966_dp, 4, 1), & ! Co 3d
    & xtbBasis( 4, 0,  0.815250_dp, 0.000000_dp, -4.593686_dp*eV, 0e+0_dp, 1.925082_dp, 6, 1), & ! Co 4s
    & xtbBasis( 4, 1, 15.765732_dp,-2.000000_dp, -3.855768_dp*eV, 0e+0_dp, 1.350000_dp, 6, 1), & ! Co 4p
    & xtbBasis( 3, 2,-26.348820_dp, 1.282426_dp,-10.841022_dp*eV, 0e+0_dp, 2.356745_dp, 4, 1), & ! Ni 3d
    & xtbBasis( 4, 0, 15.160508_dp, 0.000000_dp, -8.687611_dp*eV, 0e+0_dp, 1.532263_dp, 6, 1), & ! Ni 4s
    & xtbBasis( 4, 1, 15.782685_dp,-2.000000_dp, -3.332933_dp*eV, 0e+0_dp, 1.350000_dp, 6, 1), & ! Ni 4p
    & xtbBasis( 3, 2,-21.142399_dp,-1.290373_dp,-11.114050_dp*eV, 0e+0_dp, 2.598287_dp, 4, 1), & ! Cu 3d
    & xtbBasis( 4, 0, -3.590501_dp, 0.000000_dp, -8.373193_dp*eV, 0e+0_dp, 1.583677_dp, 6, 1), & ! Cu 4s
    & xtbBasis( 4, 1,  7.413473_dp, 0.000000_dp, -4.419045_dp*eV, 0e+0_dp, 1.350000_dp, 6, 1), & ! Cu 4p
    & xtbBasis( 4, 0,-15.535695_dp, 0.000000_dp,-11.263459_dp*eV, 0e+0_dp, 1.722526_dp, 6, 1), & ! Zn 4s
    & xtbBasis( 4, 1,  4.061664_dp, 0.200991_dp, -4.666731_dp*eV, 0e+0_dp, 1.061945_dp, 6, 1), & ! Zn 4p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 4, 0,-14.584657_dp, 0.000000_dp,-13.273222_dp*eV, 0e+0_dp, 1.992354_dp, 6, 1), & ! Ga 4s
    & xtbBasis( 4, 2, 19.671655_dp, 1.000000_dp, -2.245112_dp*eV, 0e+0_dp, 0.712761_dp, 4, 1), & ! Ga 4d
    & xtbBasis( 4, 1,  9.375082_dp,-2.021175_dp, -4.859478_dp*eV, 0e+0_dp, 1.482052_dp, 6, 1), & ! Ga 4p
    & xtbBasis( 4, 0,-12.195371_dp, 0.000000_dp,-12.558286_dp*eV, 6e-3_dp, 2.172951_dp, 6, 1), & ! Ge 4s
    & xtbBasis( 4, 1,-11.374296_dp,-7.631942_dp, -8.035796_dp*eV,-3e-3_dp, 1.794495_dp, 6, 1), & ! Ge 4p
    & xtbBasis( 4, 2,  9.364108_dp,-1.300000_dp, -2.752271_dp*eV,-5e-3_dp, 0.769997_dp, 4, 1), & ! Ge 4d
    & xtbBasis( 4, 0,-17.489686_dp, 0.000000_dp,-17.515251_dp*eV, 6e-3_dp, 2.265106_dp, 6, 1), & ! As 4s
    & xtbBasis( 4, 1, -6.747956_dp,-0.335509_dp, -8.272706_dp*eV,-3e-3_dp, 1.986411_dp, 6, 1), & ! As 4p
    & xtbBasis( 4, 2, 17.858510_dp,-1.000000_dp, -1.245776_dp*eV,-5e-3_dp, 1.113511_dp, 4, 1), & ! As 4d
    & xtbBasis( 4, 0,-14.852299_dp, 0.000000_dp,-23.000000_dp*eV, 6e-3_dp, 3.044672_dp, 6, 1), & ! Se 4s
    & xtbBasis( 4, 1, -9.863477_dp,-3.213580_dp,-10.398968_dp*eV,-3e-3_dp, 2.098532_dp, 6, 1), & ! Se 4p
    & xtbBasis( 4, 2,  9.556181_dp,-2.500000_dp, -0.821804_dp*eV,-5e-3_dp, 1.863317_dp, 4, 1), & ! Se 4d
    & xtbBasis( 4, 0,-17.815502_dp, 0.000000_dp,-19.875752_dp*eV, 6e-3_dp, 2.886237_dp, 6, 1), & ! Br 4s
    & xtbBasis( 4, 1,-14.058044_dp,-1.440020_dp,-12.818655_dp*eV,-3e-3_dp, 2.190987_dp, 6, 1), & ! Br 4p
    & xtbBasis( 4, 2,  5.468245_dp,-1.000000_dp, -3.348113_dp*eV,-5e-3_dp, 1.789395_dp, 4, 1), & ! Br 4d
    & xtbBasis( 4, 0,-25.437273_dp, 0.000000_dp,-20.280017_dp*eV, 6e-3_dp, 2.828105_dp, 6, 1), & ! Kr 4s
    & xtbBasis( 4, 1,-12.813227_dp,-3.743296_dp,-15.200155_dp*eV,-3e-3_dp, 1.965472_dp, 6, 1), & ! Kr 4p
    & xtbBasis( 4, 2, 10.440712_dp, 0.000000_dp, -4.253986_dp*eV,-5e-3_dp, 1.512609_dp, 4, 1), & ! Kr 4d
    & xtbBasis( 5, 0, -7.450752_dp, 0.000000_dp, -7.616948_dp*eV, 0e+0_dp, 0.809529_dp, 6, 1), & ! Rb 5s
    & xtbBasis( 5, 1, 16.670533_dp,-5.181667_dp, -4.369842_dp*eV, 0e+0_dp, 0.950253_dp, 6, 1), & ! Rb 5p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 5, 0, -6.087125_dp, 0.000000_dp, -6.840171_dp*eV, 0e+0_dp, 1.458742_dp, 6, 1), & ! Sr 5s
    & xtbBasis( 5, 1,  2.115262_dp,-8.003590_dp, -3.338573_dp*eV, 0e+0_dp, 0.730658_dp, 6, 1), & ! Sr 5p
    & xtbBasis( 4, 2, 17.076466_dp, 0.000000_dp, -1.715680_dp*eV, 0e+0_dp, 1.028147_dp, 4, 1), & ! Sr 4d
    & xtbBasis( 4, 2,-28.061976_dp,-4.159186_dp, -5.731066_dp*eV, 0e+0_dp, 2.300000_dp, 4, 1), & ! Y  4d
    & xtbBasis( 5, 0, 10.950764_dp, 0.000000_dp, -8.748292_dp*eV, 0e+0_dp, 1.593058_dp, 6, 1), & ! Y  5s
    & xtbBasis( 5, 1, 45.679760_dp,-0.800000_dp, -0.838555_dp*eV, 0e+0_dp, 1.170000_dp, 6, 1), & ! Y  5p
    & xtbBasis( 4, 2,-22.240873_dp, 0.337914_dp, -6.771010_dp*eV, 0e+0_dp, 2.175661_dp, 4, 1), & ! Zr 4d
    & xtbBasis( 5, 0, 44.110231_dp, 0.000000_dp, -3.979156_dp*eV, 0e+0_dp, 1.665905_dp, 6, 1), & ! Zr 5s
    & xtbBasis( 5, 1, 25.863572_dp,-2.500000_dp, -3.954049_dp*eV, 0e+0_dp, 1.230000_dp, 6, 1), & ! Zr 5p
    & xtbBasis( 4, 2,-25.998052_dp, 0.638436_dp, -9.245726_dp*eV, 0e+0_dp, 2.092288_dp, 4, 1), & ! Nb 4d
    & xtbBasis( 5, 0, 15.379439_dp, 0.000000_dp, -9.268975_dp*eV, 0e+0_dp, 1.459971_dp, 6, 1), & ! Nb 5s
    & xtbBasis( 5, 1, 30.159730_dp,-2.000000_dp, -1.348707_dp*eV, 0e+0_dp, 1.200000_dp, 6, 1), & ! Nb 5p
    & xtbBasis( 4, 2,-22.556077_dp,-3.426221_dp, -8.176239_dp*eV, 0e+0_dp, 1.891236_dp, 4, 1), & ! Mo 4d
    & xtbBasis( 5, 0,  5.815301_dp, 0.000000_dp, -7.645737_dp*eV, 0e+0_dp, 1.827996_dp, 6, 1), & ! Mo 5s
    & xtbBasis( 5, 1, 14.527159_dp,-2.500000_dp, -3.802884_dp*eV, 0e+0_dp, 1.220000_dp, 6, 1), & ! Mo 5p
    & xtbBasis( 4, 2,-23.231470_dp, 2.642680_dp, -8.690050_dp*eV, 0e+0_dp, 2.120497_dp, 4, 1), & ! Tc 4d
    & xtbBasis( 5, 0, 24.977603_dp, 0.000000_dp, -5.089073_dp*eV, 0e+0_dp, 1.789115_dp, 6, 1), & ! Tc 5s
    & xtbBasis( 5, 1,  1.953838_dp,-2.000000_dp, -4.878724_dp*eV, 0e+0_dp, 1.250000_dp, 6, 1), & ! Tc 5p
    & xtbBasis( 4, 2,-23.099524_dp, 1.772831_dp,-10.960165_dp*eV, 0e+0_dp, 2.352683_dp, 4, 1), & ! Ru 4d
    & xtbBasis( 5, 0, 15.281981_dp, 0.000000_dp, -6.304229_dp*eV, 0e+0_dp, 1.883645_dp, 6, 1), & ! Ru 5s
    & xtbBasis( 5, 1,  1.340798_dp,-1.500000_dp, -5.569969_dp*eV, 0e+0_dp, 1.370000_dp, 6, 1), & ! Ru 5p
    & xtbBasis( 4, 2,-23.540560_dp, 3.782936_dp,-11.935915_dp*eV, 0e+0_dp, 2.436353_dp, 4, 1), & ! Rh 4d
    & xtbBasis( 5, 0, 10.450086_dp, 0.000000_dp, -4.883179_dp*eV, 0e+0_dp, 2.000000_dp, 6, 1), & ! Rh 5s
    & xtbBasis( 5, 1, 15.559547_dp,-2.500000_dp, -4.427854_dp*eV, 0e+0_dp, 1.470000_dp, 6, 1), & ! Rh 5p
    & xtbBasis( 4, 2,-23.290322_dp, 3.210802_dp,-12.059626_dp*eV, 0e+0_dp, 2.528954_dp, 4, 1), & ! Pd 4d
    & xtbBasis( 5, 0, 17.475085_dp, 0.000000_dp, -5.724219_dp*eV, 0e+0_dp, 2.073217_dp, 6, 1), & ! Pd 5s
    & xtbBasis( 5, 1, 21.621321_dp,-2.500000_dp, -2.575000_dp*eV, 0e+0_dp, 1.550000_dp, 6, 1), & ! Pd 5p
    & xtbBasis( 4, 2, -6.963262_dp,-1.477715_dp, -9.675945_dp*eV, 0e+0_dp, 2.720329_dp, 4, 1), & ! Ag 4d
    & xtbBasis( 5, 0,-12.856324_dp, 0.000000_dp, -5.723081_dp*eV, 0e+0_dp, 1.994885_dp, 6, 1), & ! Ag 5s
    & xtbBasis( 5, 1, 0.1871550_dp,-1.500000_dp, -3.273430_dp*eV, 0e+0_dp, 1.620000_dp, 6, 1), & ! Ag 5p
    & xtbBasis( 5, 0,-10.281188_dp, 0.000000_dp,-12.099216_dp*eV, 0e+0_dp, 1.980518_dp, 6, 1), & ! Cd 5s
    & xtbBasis( 5, 1, 6.2471240_dp,-0.775216_dp, -3.859493_dp*eV, 0e+0_dp, 1.191810_dp, 6, 1), & ! Cd 5p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 5, 0,-10.488459_dp, 0.000000_dp,-16.894094_dp*eV, 0e+0_dp, 2.226101_dp, 6, 1), & ! In 5s
    & xtbBasis( 5, 1, 19.136222_dp, 0.762515_dp, -3.502771_dp*eV, 0e+0_dp, 1.625926_dp, 6, 1), & ! In 5p
    & xtbBasis( 5, 2, 5.5843660_dp, 0.000000_dp, -3.650350_dp*eV, 0e+0_dp, 0.663076_dp, 4, 1), & ! In 5d
    & xtbBasis( 5, 0,-19.310676_dp, 0.000000_dp,-24.164818_dp*eV, 6e-3_dp, 2.474055_dp, 6, 1), & ! Sn 5s
    & xtbBasis( 5, 1, -5.460959_dp,-3.444851_dp, -7.640096_dp*eV,-3e-3_dp, 1.893755_dp, 6, 1), & ! Sn 5p
    & xtbBasis( 5, 2, 10.683419_dp,-1.500000_dp, -1.908531_dp*eV,-5e-3_dp, 1.547485_dp, 4, 1), & ! Sn 5d
    & xtbBasis( 5, 0,-17.310388_dp, 0.000000_dp,-20.650528_dp*eV, 6e-3_dp, 2.761687_dp, 6, 1), & ! Sb 5s
    & xtbBasis( 5, 1, -7.203718_dp,-1.459812_dp, -7.536020_dp*eV,-3e-3_dp, 2.076379_dp, 6, 1), & ! Sb 5p
    & xtbBasis( 5, 2, 10.096015_dp,-2.000000_dp, -2.185884_dp*eV,-5e-3_dp, 1.071094_dp, 4, 1), & ! Sb 5d
    & xtbBasis( 5, 0,-17.836704_dp, 0.000000_dp,-29.899753_dp*eV, 6e-3_dp, 2.880945_dp, 6, 1), & ! Te 5s
    & xtbBasis( 5, 1, -9.887978_dp, 0.137154_dp,-10.026096_dp*eV,-3e-3_dp, 2.254863_dp, 6, 1), & ! Te 5p
    & xtbBasis( 5, 2, 20.942979_dp,-2.000000_dp, -0.372055_dp*eV,-5e-3_dp, 1.724516_dp, 4, 1), & ! Te 5d
    & xtbBasis( 5, 0,-21.954071_dp, 0.000000_dp,-23.832631_dp*eV, 6e-3_dp, 3.117622_dp, 6, 1), & ! I  5s
    & xtbBasis( 5, 1,-10.823970_dp,-0.387987_dp,-11.604442_dp*eV,-3e-3_dp, 2.248195_dp, 6, 1), & ! I  5p
    & xtbBasis( 5, 2, 12.522287_dp,-1.500000_dp, -2.025327_dp*eV,-5e-3_dp, 1.831809_dp, 4, 1), & ! I  5d
    & xtbBasis( 5, 0,-22.530281_dp, 0.000000_dp,-21.969064_dp*eV, 6e-3_dp, 3.128524_dp, 6, 1), & ! Xe 5s
    & xtbBasis( 5, 1,-16.667114_dp,-3.435282_dp,-11.870978_dp*eV,-3e-3_dp, 2.316580_dp, 6, 1), & ! Xe 5p
    & xtbBasis( 5, 2,  8.021956_dp,-1.500000_dp, -2.697796_dp*eV,-5e-3_dp, 1.888452_dp, 4, 1), & ! Xe 5d
    & xtbBasis( 5, 0, -1.460631_dp, 0.000000_dp, -6.341379_dp*eV, 0e+0_dp, 0.779877_dp, 6, 1), & ! Cs 5s
    & xtbBasis( 5, 1, 15.879494_dp,-7.035550_dp, -3.944275_dp*eV, 0e+0_dp, 0.810404_dp, 6, 1), & ! Cs 5p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 5, 0, -5.468018_dp, 0.000000_dp, -6.452630_dp*eV, 0e+0_dp, 1.387083_dp, 6, 1), & ! Ba 5s
    & xtbBasis( 5, 1,  4.368854_dp,-8.801363_dp, -3.975353_dp*eV, 0e+0_dp, 0.532658_dp, 6, 1), & ! Ba 5p
    & xtbBasis( 4, 2, 14.328052_dp, 0.000000_dp, -2.305768_dp*eV, 0e+0_dp, 0.853415_dp, 4, 1), & ! Ba 4d
    & xtbBasis( 5, 2,-44.208463_dp,-6.396752_dp, -5.872226_dp*eV, 0e+0_dp, 3.000000_dp, 4, 1), & ! La 5d
    & xtbBasis( 6, 0, -3.988102_dp, 0.000000_dp, -6.500000_dp*eV, 0e+0_dp, 1.492677_dp, 6, 1), & ! La 6s
    & xtbBasis( 6, 1, 40.847293_dp,-1.500000_dp, -0.727921_dp*eV, 0e+0_dp, 1.350000_dp, 6, 1), & ! La 6p
    & xtbBasis( 5, 2,-36.440945_dp,-5.245538_dp, -5.032003_dp*eV, 0e+0_dp, 3.000000_dp, 4, 1), & ! Ce 5d
    & xtbBasis( 6, 0,  6.148475_dp, 0.000000_dp, -6.275363_dp*eV, 0e+0_dp, 1.553483_dp, 6, 1), & ! Ce 6s
    & xtbBasis( 6, 1, 42.873822_dp,-1.500000_dp,  0.291196_dp*eV, 0e+0_dp, 1.380859_dp, 6, 1), & ! Ce 6p
    & xtbBasis( 5, 2,-36.021673_dp,-5.064761_dp, -4.944984_dp*eV, 0e+0_dp, 2.992307_dp, 4, 1), & ! Pr 5d
    & xtbBasis( 6, 0,  7.806576_dp, 0.000000_dp, -6.271128_dp*eV, 0e+0_dp, 1.578839_dp, 6, 1), & ! Pr 6s
    & xtbBasis( 6, 1, 42.846148_dp,-1.500000_dp,  0.241817_dp*eV, 0e+0_dp, 1.385620_dp, 6, 1), & ! Pr 6p
    & xtbBasis( 5, 2,-35.602402_dp,-4.883984_dp, -4.857964_dp*eV, 0e+0_dp, 2.984614_dp, 4, 1), & ! Nd 5d
    & xtbBasis( 6, 0,  9.464678_dp, 0.000000_dp, -6.266893_dp*eV, 0e+0_dp, 1.604196_dp, 6, 1), & ! Nd 6s
    & xtbBasis( 6, 1, 42.818474_dp,-1.500000_dp,  0.192438_dp*eV, 0e+0_dp, 1.390381_dp, 6, 1), & ! Nd 6p
    & xtbBasis( 5, 2,-35.183130_dp,-4.703207_dp, -4.770945_dp*eV, 0e+0_dp, 2.976922_dp, 4, 1), & ! Pm 5d
    & xtbBasis( 6, 0, 11.122779_dp, 0.000000_dp, -6.262657_dp*eV, 0e+0_dp, 1.629552_dp, 6, 1), & ! Pm 6s
    & xtbBasis( 6, 1, 42.790801_dp,-1.500000_dp,  0.143059_dp*eV, 0e+0_dp, 1.395142_dp, 6, 1), & ! Pm 6p
    & xtbBasis( 5, 2,-34.763859_dp,-4.522429_dp, -4.683925_dp*eV, 0e+0_dp, 2.969229_dp, 4, 1), & ! Sm 5d
    & xtbBasis( 6, 0, 12.780881_dp, 0.000000_dp, -6.258422_dp*eV, 0e+0_dp, 1.654909_dp, 6, 1), & ! Sm 6s
    & xtbBasis( 6, 1, 42.763127_dp,-1.500000_dp,  0.093680_dp*eV, 0e+0_dp, 1.399903_dp, 6, 1), & ! Sm 6p
    & xtbBasis( 5, 2,-34.344587_dp,-4.341652_dp, -4.596906_dp*eV, 0e+0_dp, 2.961536_dp, 4, 1), & ! Eu 5d
    & xtbBasis( 6, 0, 14.438982_dp, 0.000000_dp, -6.254187_dp*eV, 0e+0_dp, 1.680265_dp, 6, 1), & ! Eu 6s
    & xtbBasis( 6, 1, 42.735454_dp,-1.500000_dp,  0.044301_dp*eV, 0e+0_dp, 1.404664_dp, 6, 1), & ! Eu 6p
    & xtbBasis( 5, 2,-33.925315_dp,-4.160875_dp, -4.509886_dp*eV, 0e+0_dp, 2.953843_dp, 4, 1), & ! Gd 5d
    & xtbBasis( 6, 0, 16.097083_dp, 0.000000_dp, -6.249952_dp*eV, 0e+0_dp, 1.705622_dp, 6, 1), & ! Gd 6s
    & xtbBasis( 6, 1, 42.707780_dp,-1.500000_dp, -0.005078_dp*eV, 0e+0_dp, 1.409425_dp, 6, 1), & ! Gd 6p
    & xtbBasis( 5, 2,-33.506044_dp,-3.980098_dp, -4.422867_dp*eV, 0e+0_dp, 2.946150_dp, 4, 1), & ! Tb 5d
    & xtbBasis( 6, 0, 17.755185_dp, 0.000000_dp, -6.245716_dp*eV, 0e+0_dp, 1.730979_dp, 6, 1), & ! Tb 6s
    & xtbBasis( 6, 1, 42.680106_dp,-1.500000_dp, -0.054457_dp*eV, 0e+0_dp, 1.414186_dp, 6, 1), & ! Tb 6p
    & xtbBasis( 5, 2,-33.086772_dp,-3.799321_dp, -4.335848_dp*eV, 0e+0_dp, 2.938457_dp, 4, 1), & ! Dy 5d
    & xtbBasis( 6, 0, 19.413286_dp, 0.000000_dp, -6.241481_dp*eV, 0e+0_dp, 1.756335_dp, 6, 1), & ! Dy 6s
    & xtbBasis( 6, 1, 42.652433_dp,-1.500000_dp, -0.103836_dp*eV, 0e+0_dp, 1.418947_dp, 6, 1), & ! Dy 6p
    & xtbBasis( 5, 2,-32.667501_dp,-3.618544_dp, -4.248828_dp*eV, 0e+0_dp, 2.930765_dp, 4, 1), & ! Ho 5d
    & xtbBasis( 6, 0, 21.071387_dp, 0.000000_dp, -6.237246_dp*eV, 0e+0_dp, 1.781692_dp, 6, 1), & ! Ho 6s
    & xtbBasis( 6, 1, 42.624759_dp,-1.500000_dp, -0.153215_dp*eV, 0e+0_dp, 1.423708_dp, 6, 1), & ! Ho 6p
    & xtbBasis( 5, 2,-32.248229_dp,-3.437767_dp, -4.161809_dp*eV, 0e+0_dp, 2.923072_dp, 4, 1), & ! Er 5d
    & xtbBasis( 6, 0, 22.729489_dp, 0.000000_dp, -6.233011_dp*eV, 0e+0_dp, 1.807048_dp, 6, 1), & ! Er 6s
    & xtbBasis( 6, 1, 42.597085_dp,-1.500000_dp, -0.202593_dp*eV, 0e+0_dp, 1.428469_dp, 6, 1), & ! Er 6p
    & xtbBasis( 5, 2,-31.828957_dp,-3.256989_dp, -4.074789_dp*eV, 0e+0_dp, 2.915379_dp, 4, 1), & ! Tm 5d
    & xtbBasis( 6, 0, 24.387590_dp, 0.000000_dp, -6.228775_dp*eV, 0e+0_dp, 1.832405_dp, 6, 1), & ! Tm 6s
    & xtbBasis( 6, 1, 42.569412_dp,-1.500000_dp, -0.251972_dp*eV, 0e+0_dp, 1.433230_dp, 6, 1), & ! Tm 6p
    & xtbBasis( 5, 2,-31.409686_dp,-3.076212_dp, -3.987770_dp*eV, 0e+0_dp, 2.907686_dp, 4, 1), & ! Yb 5d
    & xtbBasis( 6, 0, 26.045692_dp, 0.000000_dp, -6.224540_dp*eV, 0e+0_dp, 1.857761_dp, 6, 1), & ! Yb 6s
    & xtbBasis( 6, 1, 42.541738_dp,-1.500000_dp, -0.301351_dp*eV, 0e+0_dp, 1.437991_dp, 6, 1), & ! Yb 6p
    & xtbBasis( 5, 2,-30.990414_dp,-2.895435_dp, -3.900750_dp*eV, 0e+0_dp, 2.899993_dp, 4, 1), & ! Lu 5d
    & xtbBasis( 6, 0, 27.703793_dp, 0.000000_dp, -6.220305_dp*eV, 0e+0_dp, 1.883118_dp, 6, 1), & ! Lu 6s
    & xtbBasis( 6, 1, 42.514065_dp,-1.500000_dp, -0.350730_dp*eV, 0e+0_dp, 1.442752_dp, 6, 1), & ! Lu 6p
    & xtbBasis( 5, 2,-21.116286_dp,-1.485678_dp, -4.360558_dp*eV, 0e+0_dp, 2.466693_dp, 4, 1), & ! Hf 5d
    & xtbBasis( 6, 0, 15.014122_dp, 0.000000_dp, -5.910623_dp*eV, 0e+0_dp, 2.039390_dp, 6, 1), & ! Hf 6s
    & xtbBasis( 6, 1, 22.898249_dp,-1.500000_dp, -2.814338_dp*eV, 0e+0_dp, 1.450000_dp, 6, 1), & ! Hf 6p
    & xtbBasis( 5, 2,-23.077812_dp,-1.870583_dp, -9.232014_dp*eV, 0e+0_dp, 2.177327_dp, 4, 1), & ! Ta 5d
    & xtbBasis( 6, 0, 29.782424_dp, 0.000000_dp, -8.600553_dp*eV, 0e+0_dp, 1.692963_dp, 6, 1), & ! Ta 6s
    & xtbBasis( 6, 1, 36.420564_dp,-1.500000_dp, -0.252865_dp*eV, 0e+0_dp, 1.400000_dp, 6, 1), & ! Ta 6p
    & xtbBasis( 5, 2,-17.030630_dp, 0.130920_dp, -8.997799_dp*eV, 0e+0_dp, 2.300752_dp, 4, 1), & ! W  5d
    & xtbBasis( 6, 0, 35.195571_dp, 0.000000_dp, -2.878936_dp*eV, 0e+0_dp, 2.096013_dp, 6, 1), & ! W  6s
    & xtbBasis( 6, 1, 18.760746_dp,-1.500000_dp, -3.369287_dp*eV, 0e+0_dp, 1.400000_dp, 6, 1), & ! W  6p
    & xtbBasis( 5, 2,-23.115824_dp, 2.507095_dp, -7.858164_dp*eV, 0e+0_dp, 2.470782_dp, 4, 1), & ! Re 5d
    & xtbBasis( 6, 0, 23.560994_dp, 0.000000_dp, -6.430285_dp*eV, 0e+0_dp, 2.220548_dp, 6, 1), & ! Re 6s
    & xtbBasis( 6, 1, -0.067497_dp,-2.000000_dp, -5.165147_dp*eV, 0e+0_dp, 1.450000_dp, 6, 1), & ! Re 6p
    & xtbBasis( 5, 2,-19.564083_dp,-0.262294_dp,-10.716969_dp*eV, 0e+0_dp, 2.734340_dp, 4, 1), & ! Os 5d
    & xtbBasis( 6, 0, 24.928002_dp, 0.000000_dp, -3.655133_dp*eV, 0e+0_dp, 2.365840_dp, 6, 1), & ! Os 6s
    & xtbBasis( 6, 1, -4.330556_dp,-1.000000_dp, -7.060522_dp*eV, 0e+0_dp, 1.650000_dp, 6, 1), & ! Os 6p
    & xtbBasis( 5, 2,-21.172493_dp, 3.805255_dp,-12.054598_dp*eV, 0e+0_dp, 2.797508_dp, 4, 1), & ! Ir 5d
    & xtbBasis( 6, 0, 25.774929_dp, 0.000000_dp, -5.686006_dp*eV, 0e+0_dp, 2.274300_dp, 6, 1), & ! Ir 6s
    & xtbBasis( 6, 1, -0.704597_dp,-1.500000_dp, -6.208990_dp*eV, 0e+0_dp, 1.650000_dp, 6, 1), & ! Ir 6p
    & xtbBasis( 5, 2,-22.169385_dp, 0.996400_dp,-11.571582_dp*eV, 0e+0_dp, 2.807068_dp, 4, 1), & ! Pt 5d
    & xtbBasis( 6, 0, 38.415536_dp, 0.000000_dp, -7.184794_dp*eV, 0e+0_dp, 2.341428_dp, 6, 1), & ! Pt 6s
    & xtbBasis( 6, 1, -0.665483_dp,-2.000000_dp, -5.080419_dp*eV, 0e+0_dp, 1.650000_dp, 6, 1), & ! Pt 6p
    & xtbBasis( 5, 2,-11.067532_dp,-4.380921_dp,-10.047575_dp*eV, 0e+0_dp, 3.117733_dp, 4, 1), & ! Au 5d
    & xtbBasis( 6, 0,-11.443658_dp, 0.000000_dp, -6.530840_dp*eV, 0e+0_dp, 2.325119_dp, 6, 1), & ! Au 6s
    & xtbBasis( 6, 1, -5.119735_dp,-1.500000_dp, -3.296026_dp*eV, 0e+0_dp, 1.750000_dp, 6, 1), & ! Au 6p
    & xtbBasis( 6, 0, -6.581368_dp, 0.000000_dp,-12.452637_dp*eV, 0e+0_dp, 2.062597_dp, 6, 1), & ! Hg 6s
    & xtbBasis( 6, 1,  3.995243_dp,-4.204099_dp, -4.169731_dp*eV, 0e+0_dp, 1.721925_dp, 6, 1), & ! Hg 6p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 6, 0, -2.193199_dp, 0.000000_dp,-12.563376_dp*eV, 0e+0_dp, 2.647541_dp, 6, 1), & ! Tl 6s
    & xtbBasis( 6, 1,  0.060451_dp,-8.101017_dp, -5.131043_dp*eV, 0e+0_dp, 1.717991_dp, 6, 1), & ! Tl 6p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 6, 0,-10.874138_dp, 0.000000_dp,-14.496335_dp*eV, 6e-3_dp, 2.847707_dp, 6, 1), & ! Pb 6s
    & xtbBasis( 6, 1, -6.034796_dp,-7.925216_dp, -5.848584_dp*eV,-3e-3_dp, 2.068091_dp, 6, 1), & ! Pb 6p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 6, 0,-20.410234_dp, 0.000000_dp,-18.199529_dp*eV, 6e-3_dp, 2.895660_dp, 6, 1), & ! Bi 6s
    & xtbBasis( 6, 1, -9.424568_dp,-7.150589_dp, -6.735929_dp*eV,-3e-3_dp, 2.256279_dp, 6, 1), & ! Bi 6p
    & xtbBasis(-1,-1,  0.000000_dp, 0.000000_dp,  0.000000_dp*eV, 0e+0_dp, 0.000000_dp, 0, 0), & ! -- --
    & xtbBasis( 6, 0,-18.477865_dp, 0.000000_dp,-23.908422_dp*eV, 6e-3_dp, 3.150662_dp, 6, 1), & ! Po 6s
    & xtbBasis( 6, 1,-14.037423_dp,-3.955914_dp, -8.889548_dp*eV,-3e-3_dp, 2.382063_dp, 6, 1), & ! Po 6p
    & xtbBasis( 5, 2, 13.809093_dp, 0.000000_dp, -0.921251_dp*eV,-5e-3_dp, 1.241625_dp, 4, 1), & ! Po 5d
    & xtbBasis( 6, 0,-21.965390_dp, 0.000000_dp,-21.752193_dp*eV, 6e-3_dp, 3.516922_dp, 6, 1), & ! At 6s
    & xtbBasis( 6, 1,-12.804436_dp,-3.402676_dp,-10.031093_dp*eV,-3e-3_dp, 2.392024_dp, 6, 1), & ! At 6p
    & xtbBasis( 5, 2, 16.836546_dp, 0.000000_dp, -0.852571_dp*eV,-5e-3_dp, 1.380239_dp, 4, 1), & ! At 5d
    & xtbBasis( 6, 0,-22.139701_dp, 0.000000_dp,-18.381647_dp*eV, 6e-3_dp, 3.520683_dp, 6, 1), & ! Rn 6s
    & xtbBasis( 6, 1,-20.539955_dp,-2.380762_dp,-10.236606_dp*eV,-3e-3_dp, 2.535389_dp, 6, 1), & ! Rn 6p
    & xtbBasis( 5, 2, 17.249637_dp, 0.000000_dp, -0.973687_dp*eV,-5e-3_dp, 1.418875_dp, 4, 1)],& ! Rn 5d
    & shape(gfn1Basis))


contains


  !> Return GFN1-xTB parameters for species with given atomic number
  function getGFN1ParamSymbol(symbol) result(param)

    !> Atomic number of species
    character(len=*), intent(in) :: symbol

    !> GFN1-xTB parameters
    type(xtbParam) :: param

    integer :: iSh

    param = getGFN1Param(symbolToNumber(symbol))

  end function getGFN1ParamSymbol


  !> Return GFN1-xTB parameters for species with given atomic number
  function getGFN1ParamNumber(number) result(param)

    !> Atomic number of species
    integer, intent(in) :: number

    !> GFN1-xTB parameters
    type(xtbParam) :: param

    integer :: iSh

    if (number > 0 .and. number <= gfn1Elem) then
      param = gfn1Param(number)
      allocate(param%basis(param%nSh))
      do iSh = 1, param%nSh
        param%basis(iSh) = gfn1Basis(iSh, number)
        param%basis(iSh)%kcn = param%basis(iSh)%h * param%basis(iSh)%kcn
      end do
    else
      param = xtbParam(0, 0, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)
    end if

  end function getGFN1ParamNumber


end module dftbp_gfn1param
