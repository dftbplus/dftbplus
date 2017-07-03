!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains routines to write out various data structures in a comprehensive
!!* tagged format.
module taggedoutput
  use assert
  use accuracy, only : dp
  implicit none
  private

  public :: initTaggedWriter, writeTagged

  !!* Writes objects in standardized form to the output
  interface writeTagged
    module procedure writeTaggedRealR0
    module procedure writeTaggedRealR1
    module procedure writeTaggedRealR2
    module procedure writeTaggedRealR3
    module procedure writeTaggedRealR4
    module procedure writeTaggedComplexR0
    module procedure writeTaggedComplexR1
    module procedure writeTaggedComplexR2
    module procedure writeTaggedComplexR3
    module procedure writeTaggedComplexR4
    module procedure writeTaggedIntegerR0
    module procedure writeTaggedIntegerR1
    module procedure writeTaggedIntegerR2
    module procedure writeTaggedIntegerR3
    module procedure writeTaggedIntegerR4
    module procedure writeTaggedLogicalR0
    module procedure writeTaggedLogicalR1
    module procedure writeTaggedLogicalR2
    module procedure writeTaggedLogicalR3
    module procedure writeTaggedLogicalR4
  end interface

  integer, parameter :: lenLabel = 20


  !! Tag names (Should be shorter than lenLabel!)
  character(*), parameter, public :: tag_SCC        = 'scc'
  character(*), parameter, public :: tag_nSCC       = 'n_scc_iters'
  character(*), parameter, public :: tag_sccConv    = 'scc_convergence'
  character(*), parameter, public :: tag_spin       = 'n_spins'
  character(*), parameter, public :: tag_DFTBU      = 'dftb+u'
  character(*), parameter, public :: tag_nDFTBU     = 'dftb+u_functional'
  character(*), parameter, public :: tag_LS         = 'ls'
  character(*), parameter, public :: tag_LSdual     = 'ls_dual'
  character(*), parameter, public :: tag_species     = 'species'
  character(*), parameter, public :: tag_mAngSpecies = 'angular_momenta'
  character(*), parameter, public :: tag_initCoord  = 'init_coords'
  character(*), parameter, public :: tag_endCoord   = 'end_coords'
  character(*), parameter, public :: tag_kPoint     = 'k_points'
  character(*), parameter, public :: tag_kWeight    = 'k_weights'
  character(*), parameter, public :: tag_nSKGrid    = 'sk_n_gridpoints'
  character(*), parameter, public :: tag_skDist     = 'sk_grid_distances'
  character(*), parameter, public :: tag_skHam      = 'sk_hamiltonian'
  character(*), parameter, public :: tag_skOver     = 'sk_overlap'
  character(*), parameter, public :: tag_nRepGrid   = 'rep_n_gridpoints'
  character(*), parameter, public :: tag_repGrid    = 'rep_gripoints'
  character(*), parameter, public :: tag_repCoeff   = 'rep_coeffs'
  character(*), parameter, public :: tag_repCutoff  = 'rep_cutoffs'
  character(*), parameter, public :: tag_repExp     = 'rep_exponentials'
  character(*), parameter, public :: tag_atomEigVal = 'atomic_eigenvalues'
  character(*), parameter, public :: tag_hubbU      = 'hubbard_us'
  character(*), parameter, public :: tag_tempElec   = 'electronic_temp'
  character(*), parameter, public :: tag_distribFn  = 'electron_distrib_fn'
  character(*), parameter, public :: tag_nElUp      = 'n_up_electrons'
  character(*), parameter, public :: tag_nElDown    = 'n_down_electrons'
  character(*), parameter, public :: tag_eigenVal   = 'eigenvalues'
  character(*), parameter, public :: tag_egyBand    = 'band_energy'
  character(*), parameter, public :: tag_egyBandT0  = 'band_energy_t0'
  character(*), parameter, public :: tag_egyRep     = 'repulsive_energy'
  character(*), parameter, public :: tag_qOutput    = 'orbital_charges'
  character(*), parameter, public :: tag_forces     = 'forces_calculated'
  character(*), parameter, public :: tag_forceTot   = 'forces'
  character(*), parameter, public :: tag_forceBand  = 'electronic_forces'
  character(*), parameter, public :: tag_forceRep   = 'repulsive_forces'
  character(*), parameter, public :: tag_stressRep  = 'repulsive_stress'
  character(*), parameter, public :: tag_stressElec = 'electronic_stress'
  character(*), parameter, public :: tag_volume     = 'cell_volume'
  character(*), parameter, public :: tag_stressKE   = 'kinetic_stress'
  character(*), parameter, public :: tag_stressTot  = 'stress'
  character(*), parameter, public :: tag_pV         = 'pv'
  character(*), parameter, public :: tag_nNeighbor  = 'n_neighbors'
  character(*), parameter, public :: tag_iNeighbor  = 'i_neighbors'
  character(*), parameter, public :: tag_egySCC     = 'scc_energy'
  character(*), parameter, public :: tag_egySpin    = 'spin_energy'
  character(*), parameter, public :: tag_egyDFTBU   = 'dftb+u_energy'
  character(*), parameter, public :: tag_egyLS      = 'ls_energy'
  character(*), parameter, public :: tag_egyExt     = 'extfield_energy'
  character(*), parameter, public :: tag_egyTotal   = 'total_energy'
  character(*), parameter, public :: tag_entropy    = 'entropy'
  character(*), parameter, public :: tag_freeEgy    = 'mermin_energy'
  character(*), parameter, public :: tag_Gibbsfree  = 'gibbs_energy'
  character(*), parameter, public :: tag_filling    = 'fillings'
  character(*), parameter, public :: tag_egyTotElec = 'total_elec_energy'
  character(*), parameter, public :: tag_qOutputAt  = 'atomic_charges'
  character(*), parameter, public :: tag_qOutAtNet  = 'net_atomic_charges'
  character(*), parameter, public :: tag_chrgForces = 'forces_ext_charges'
  character(*), parameter, public :: tag_efermi     = 'fermi_level'
  character(*), parameter, public :: tag_dispersn   = 'dispersion'
  character(*), parameter, public :: tag_egyDispersn= 'dispersion_energy'
  character(*), parameter, public :: tag_egyDispAt  = 'atomic_dispn_energy'
  character(*), parameter, public :: tag_egyRepAt   = 'atomic_rep_energy'
  character(*), parameter, public :: tag_egySCCAt   = 'atomic_scc_energy'
  character(*), parameter, public :: tag_egySpinAt  = 'atomic_spin_energy'
  character(*), parameter, public :: tag_egyDFTBUAt = 'atomic_+u_energy'
  character(*), parameter, public :: tag_egyLSAt    = 'atomic_ls'
  character(*), parameter, public :: tag_egyExtAt   = 'atomic_extfield_energy'
  character(*), parameter, public :: tag_egyTotElAt = 'atomic_elec_energy'
  character(*), parameter, public :: tag_egyTotalAt = 'atomic_egyTotal'
  character(*), parameter, public :: tag_HessianNum = 'hessian_numerical'
  character(*), parameter, public :: tag_pmlocalise = 'pm_localisation'
  ! linear response Casida related tags:
  character(*), parameter, public :: tag_spExcEgy   = 'sp_exc_energies'
  character(*), parameter, public :: tag_spExcOsc   = 'sp_exc_oscillator'
  character(*), parameter, public :: tag_excEgy     = 'exc_energies_sqr'
  character(*), parameter, public :: tag_excOsc     = 'exc_oscillator'
  character(*), parameter, public :: tag_excCharges = 'exc_charges'
  character(*), parameter, public :: tag_excForce   = 'exc_forces'
  ! Ground state dipole
  character(*), parameter, public :: tag_dipole     = 'dipole'

  character(len=lenLabel) :: formReal
  character(len=lenLabel) :: formCmplx
  character(len=lenLabel) :: formInt
  character(len=lenLabel) :: formLogical

  logical :: initialized = .false.

contains


  subroutine initTaggedWriter()

    integer :: nDecDigit, nExpDigit, nChar, nField

    if (initialized) then
      return
    end if

    !! "-3.1234567E-123 ": nDec = 7, nExpDigit = 3, nChar = 16
    nExpDigit = ceiling(log(maxexponent(1.0_dp)/log(10.0))/log(10.0))
    nDecDigit = precision(1.0_dp)
    nChar = nDecDigit + nExpDigit + 6
    nField = 80 / nChar
    if (nField == 0) then
      nField = 1
    end if
99000 format ('(', I2.2, 'E', I2.2, '.', I2.2, 'E', I3.3, ')')
    write (formReal, 99000) &
        & nField, nChar, nDecDigit, nExpDigit
99010 format ('(', I2.2, '(2E', I2.2, '.', I2.2, 'E', I3.3, '))')
    write (formCmplx, 99010) &
        & nField/2, nChar, nDecDigit, nExpDigit

    !! "-12345 "
    nChar = digits(1) + 2
    nField = 80 / nChar
    if (nField == 0) then
      nField = 1
    end if
99020 format ('(', I2.2, 'I', I2.2, ')')
    write (formInt, 99020) nField, nChar

99030 format ('(40L2)')
    write (formLogical, 99030)

    initialized = .true.

  end subroutine initTaggedWriter



  subroutine writeTaggedRealR0(file, tag, value, optForm)
    integer,                    intent(in) :: file
    character(len=*),           intent(in) :: tag
    real(dp),                   intent(in) :: value
    character(len=*), optional, intent(in) :: optForm

    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formReal)
    end if

99040 format (A, ':real:0:')
    write (file, 99040) getLabel(tag)
    write (file, form) value

  end subroutine writeTaggedRealR0



  subroutine writeTaggedRealR1(file, tag, value, optForm)
    integer,                    intent(in) :: file
    character(len=*),           intent(in) :: tag
    real(dp),                   intent(in) :: value(:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formReal)
    end if


99050 format (A, ':real:1:', I0)
    write (file, 99050) getLabel(tag), size(value)
    write (file, form) (value(ii), ii = 1, size(value))
  end subroutine writeTaggedRealR1



  subroutine writeTaggedRealR2(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    real(dp),         intent(in) :: value(:,:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii, jj
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formReal)
    end if

99060 format (A, ':real:2:', I0, ',', I0)
    write (file, 99060) getLabel(tag), &
        & size(value, dim=1), size(value, dim=2)
    write (file, form) ((value(ii, jj), &
        & ii = 1, size(value, dim=1)), &
        & jj = 1, size(value, dim=2))
  end subroutine writeTaggedRealR2



  subroutine writeTaggedRealR3(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    real(dp),         intent(in) :: value(:,:,:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii, jj, kk
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formReal)
    end if

99070 format (A, ':real:3:', I0, ',', I0, ',', I0)
    write (file, 99070) getLabel(tag), &
        & size(value, dim=1), size(value, dim=2), size(value, dim=3)
    write (file, form) (((value(ii, jj, kk), &
        & ii = 1, size(value, dim=1)), &
        & jj = 1, size(value, dim=2)), &
        & kk = 1, size(value, dim=3))
  end subroutine writeTaggedRealR3



  subroutine writeTaggedRealR4(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    real(dp),         intent(in) :: value(:,:,:,:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii, jj, kk, ll
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formReal)
    end if


99080 format (A, ':real:4:', I0, ',', I0, ',', I0, ',', I0)
    write (file, 99080) getLabel(tag), &
        & size(value, dim=1), size(value, dim=2), size(value, dim=3), &
        & size(value, dim=4)
    write (file, form) ((((value(ii, jj, kk, ll), &
        & ii = 1, size(value, dim=1)), &
        & jj = 1, size(value, dim=2)), &
        & kk = 1, size(value, dim=3)), &
        & ll = 1, size(value, dim=4))
  end subroutine writeTaggedRealR4



  subroutine writeTaggedComplexR0(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    complex(dp),      intent(in) :: value
    character(len=*), optional, intent(in) :: optForm

    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formCmplx)
    end if

99090 format (A, ':complex:0:')
    write (file, 99090) getLabel(tag)
    write (file, form) value

  end subroutine writeTaggedComplexR0



  subroutine writeTaggedComplexR1(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    complex(dp),      intent(in) :: value(:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formCmplx)
    end if

99100 format (A, ':complex:1:', I0)
    write (file, 99100) getLabel(tag), size(value)
    write (file, form) (value(ii), ii = 1, size(value))
  end subroutine writeTaggedComplexR1



  subroutine writeTaggedComplexR2(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    complex(dp),      intent(in) :: value(:,:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii, jj
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formCmplx)
    end if

99110 format (A, ':complex:2:', I0, ',', I0)
    write (file, 99110) getLabel(tag), &
        & size(value, dim=1), size(value, dim=2)
    write (file, form) ((value(ii, jj), &
        & ii = 1, size(value, dim=1)), &
        & jj = 1, size(value, dim=2))
  end subroutine writeTaggedComplexR2



  subroutine writeTaggedComplexR3(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    complex(dp),      intent(in) :: value(:,:,:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii, jj, kk
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formCmplx)
    end if

99120 format (A, ':complex:3:', I0, ',', I0, ',', I0)
    write (file, 99120) getLabel(tag), &
        & size(value, dim=1), size(value, dim=2), size(value, dim=3)
    write (file, form) (((value(ii, jj, kk), &
        & ii = 1, size(value, dim=1)), &
        & jj = 1, size(value, dim=2)), &
        & kk = 1, size(value, dim=3))
  end subroutine writeTaggedComplexR3



  subroutine writeTaggedComplexR4(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    complex(dp),      intent(in) :: value(:,:,:,:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii, jj, kk, ll
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formCmplx)
    end if

99130 format (A, ':complex:4:', I0, ',', I0, ',', I0, ',', I0)
    write (file, 99130) getLabel(tag), &
        & size(value, dim=1), size(value, dim=2), size(value, dim=3), &
        & size(value, dim=4)
    write (file, form) ((((value(ii, jj, kk, ll), &
        & ii = 1, size(value, dim=1)), &
        & jj = 1, size(value, dim=2)), &
        & kk = 1, size(value, dim=3)), &
        & ll = 1, size(value, dim=4))
  end subroutine writeTaggedComplexR4



  subroutine writeTaggedIntegerR0(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    integer,          intent(in) :: value
    character(len=*), optional, intent(in) :: optForm

    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formInt)
    end if

99140 format (A, ':integer:0:')
    write (file, 99140) getLabel(tag)
    write (file, form) value

  end subroutine writeTaggedIntegerR0



  subroutine writeTaggedIntegerR1(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    integer,          intent(in) :: value(:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formInt)
    end if

99150 format (A, ':integer:1:', I0)
    write (file, 99150) getLabel(tag), size(value)
    write (file, form) (value(ii), ii = 1, size(value))
  end subroutine writeTaggedIntegerR1



  subroutine writeTaggedIntegerR2(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    integer,          intent(in) :: value(:,:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii, jj
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formInt)
    end if

99160 format (A, ':integer:2:', I0, ',', I0)
    write (file, 99160) getLabel(tag), &
        & size(value, dim=1), size(value, dim=2)
    write (file, form) ((value(ii, jj), &
        & ii = 1, size(value, dim=1)), &
        & jj = 1, size(value, dim=2))
  end subroutine writeTaggedIntegerR2



  subroutine writeTaggedIntegerR3(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    integer,          intent(in) :: value(:,:,:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii, jj, kk
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formInt)
    end if

99170 format (A, ':integer:3:', I0, ',', I0, ',', I0)
    write (file, 99170) getLabel(tag), &
        & size(value, dim=1), size(value, dim=2), size(value, dim=3)
    write (file, form) (((value(ii, jj, kk), &
        & ii = 1, size(value, dim=1)), &
        & jj = 1, size(value, dim=2)), &
        & kk = 1, size(value, dim=3))
  end subroutine writeTaggedIntegerR3



  subroutine writeTaggedIntegerR4(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    integer,          intent(in) :: value(:,:,:,:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii, jj, kk, ll
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formInt)
    end if

99180 format (A, ':integer:4:', I0, ',', I0, ',', I0, ',', I0)
    write (file, 99180) getLabel(tag), &
        & size(value, dim=1), size(value, dim=2), size(value, dim=3), &
        & size(value, dim=4)
    write (file, form) ((((value(ii, jj, kk, ll), &
        & ii = 1, size(value, dim=1)), &
        & jj = 1, size(value, dim=2)), &
        & kk = 1, size(value, dim=3)), &
        & ll = 1, size(value, dim=4))
  end subroutine writeTaggedIntegerR4



  subroutine writeTaggedLogicalR0(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    logical,          intent(in) :: value
    character(len=*), optional, intent(in) :: optForm

    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formLogical)
    end if

99190 format (A, ':logical:0:')
    write (file, 99190) getLabel(tag)
    write (file, form) value

  end subroutine writeTaggedLogicalR0



  subroutine writeTaggedLogicalR1(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    logical,          intent(in) :: value(:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formLogical)
    end if

99200 format (A, ':logical:1:', I0)
    write (file, 99200) getLabel(tag), size(value)
    write (file, form) (value(ii), ii = 1, size(value))
  end subroutine writeTaggedLogicalR1



  subroutine writeTaggedLogicalR2(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    logical,          intent(in) :: value(:,:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii, jj
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formLogical)
    end if

99210 format (A, ':logical:2:', I0, ',', I0)
    write (file, 99210) getLabel(tag), &
        & size(value, dim=1), size(value, dim=2)
    write (file, form) ((value(ii, jj), &
        & ii = 1, size(value, dim=1)), &
        & jj = 1, size(value, dim=2))
  end subroutine writeTaggedLogicalR2



  subroutine writeTaggedLogicalR3(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    logical,          intent(in) :: value(:,:,:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii, jj, kk
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formLogical)
    end if

99220 format (A, ':logical:3:', I0, ',', I0, ',', I0)
    write (file, 99220) getLabel(tag), &
        & size(value, dim=1), size(value, dim=2), size(value, dim=3)
    write (file, form) (((value(ii, jj, kk), &
        & ii = 1, size(value, dim=1)), &
        & jj = 1, size(value, dim=2)), &
        & kk = 1, size(value, dim=3))
  end subroutine writeTaggedLogicalR3



  subroutine writeTaggedLogicalR4(file, tag, value, optForm)
    integer,          intent(in) :: file
    character(len=*), intent(in) :: tag
    logical,          intent(in) :: value(:,:,:,:)
    character(len=*), optional, intent(in) :: optForm

    integer :: ii, jj, kk, ll
    character(len=20) :: form

    @:ASSERT(initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(formLogical)
    end if

99230 format (A, ':logical:4:', I0, ',', I0, ',', I0, ',', I0)
    write (file, 99230) getLabel(tag), &
        & size(value, dim=1), size(value, dim=2), size(value, dim=3), &
        & size(value, dim=4)
    write (file, form) ((((value(ii, jj, kk, ll), &
        & ii = 1, size(value, dim=1)), &
        & jj = 1, size(value, dim=2)), &
        & kk = 1, size(value, dim=3)), &
        & ll = 1, size(value, dim=4))
  end subroutine writeTaggedLogicalR4



  character(len=20) function getLabel(tag)
    character(len=*), intent(in) :: tag

    integer :: lentrim

    @:ASSERT(initialized)

    lentrim = len_trim(tag)
    if (lentrim >= lenLabel) then
      getLabel(:) = tag(1:lenLabel)
    else
      getLabel(1:lentrim) = tag(1:lentrim)
      getLabel(lentrim+1:lenLabel) = " "
    end if

  end function getLabel


end module taggedoutput
