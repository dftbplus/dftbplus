!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#! Template paramters: (data type, suffix name, data type name in the tagged output, format string)
#:set TEMPLATE_PARAMS = [('real(dp)', 'Real', 'real', 'formReal'),&
    & ('complex(dp)', 'Cplx', 'complex', 'formCmplx'),&
    & ('integer', 'Integer', 'integer', 'formInt'),&
    & ('logical', 'Logical', 'logical', 'formLogical')]

#! Maximal rank to include into the interface (string from 0 - scalar)
#:set MAX_RANK = 4


!> Contains routines to write out various data structures in a comprehensive tagged format.
module dftbp_taggedoutput
  use dftbp_assert
  use dftbp_accuracy, only : dp
  implicit none
  private

  public :: tagLabels
  public :: TTaggedWriter, TTaggedWriter_init


  !> Length of permissible tag labels. Tag names should be shorter than lenLabel!
  integer, parameter :: lenLabel = 20

  !> Max length of the format strings for individual items
  integer, parameter :: lenFormStr = 20


  !> Contains a writer to write data in tagged form.
  type :: TTaggedWriter
    private
    ! Format strings
    character(len=lenFormStr) :: formReal
    character(len=lenFormStr) :: formCmplx
    character(len=lenFormStr) :: formInt
    character(len=lenFormStr) :: formLogical
    logical :: initialized = .false.
  contains
    #:for _, SUFFIX, _, _ in TEMPLATE_PARAMS
      #:for RANK in range(MAX_RANK + 1)
        procedure, private :: write${SUFFIX}$${RANK}$ => TTaggedWriter_write${SUFFIX}$${RANK}$
        generic :: write => write${SUFFIX}$${RANK}$
      #:endfor
    #:endfor
  end type TTaggedWriter


  !> Enumeration of the possible tag labels
  type :: TTagLabelsEnum

    !> unit cell volume (periodic)
    character(lenLabel) :: volume = 'cell_volume'

    !> final geometry
    character(lenLabel) :: endCoord = 'end_coords'

    !> excitation energies in Casida formalism
    character(lenLabel) :: excEgy = 'exc_energies_sqr'

    !> excited state force contributions
    character(lenLabel) :: excForce = 'exc_forces'

    !> oscillator strength for excitations
    character(lenLabel) :: excOsc = 'exc_oscillator'

    !> Transition dipole moments for excitations
    character(lenLabel) :: excDipole = 'exc_transdip'

    !> nonadiabatic coupling vector, H
    character(lenLabel) :: nacH = 'coupling_vectors'

    !> ground state total forces
    character(lenLabel) :: forceTot = 'forces'

    !> forces on any external charges present
    character(lenLabel) :: chrgForces = 'forces_ext_charges'

    !> Fermi level(s)
    character(lenLabel) :: fermiLvl = 'fermi_level'

    !> number of electrons
    character(lenLabel) :: nElec = 'number_of_electrons'

    !> eigenvalues/single particle states
    character(lenLabel) :: eigvals = 'eigenvalues'

    !> filling of the eigenstates
    character(lenLabel) :: eigFill = 'filling'

    !> Gibbs free energy for finite pressure periodic systems
    character(lenLabel) :: gibbsFree = 'gibbs_energy'

    !> Gross atomic charges
    character(lenLabel) :: qOutAtGross  = 'gross_atomic_charges'

    !> Charge model 5 corrected atomic gross charges
    character(lenLabel) :: qOutAtCM5 = 'cm5_atomic_charges'

    !> numerically calculated second derivatives matrix
    character(lenLabel) :: hessianNum = 'hessian_numerical'

    !> final energy components after real-time propagation
    character(lenLabel) :: tdenergy = 'final_energy'

    !> final dipole moment vector after real-time propagation
    character(lenLabel) :: tddipole = 'final_dipole_moment'

    !> final negative gross atomic Mulliken charges after real-time propagation
    character(lenLabel) :: tdcharges = 'final_td_charges'

    !> final forces components after real-time (Ehrenfest) propagation
    character(lenLabel) :: ehrenforces = 'final_ehrenfest_forc'

    !> final geometry after real-time (Ehrenfest) propagation
    character(lenLabel) :: ehrencoords = 'final_ehrenfest_geom'

    !> final velocities after real-time (Ehrenfest) propagation
    character(lenLabel) :: ehrenvelos = 'final_ehrenfest_velo'

    !> final molecular orbitals occupations after real-time (Ehrenfest) propagation
    character(lenLabel) :: tdprojocc = 'final_td_proj_occ'

    !> Sum of bond populaion values (should be number of electrons)
    character(lenLabel) :: sumBondPopul = 'sum_bond_pops'

    !> total energy including electron TS contribution
    character(lenLabel) :: freeEgy = 'mermin_energy'

    !> Mulliken charges
    character(lenLabel) :: qOutput = 'orbital_charges'

    !> Pipek-Mezey localisation score of single particle levels
    character(lenLabel) :: pmlocalise = 'pm_localisation'

    !> total stress tensor for periodic geometries
    character(lenLabel) :: stressTot = 'stress'

    !> total tunneling vector
    character(lenLabel) :: tunn = 'total_tunneling'

    !> total projected DOS vector
    character(lenLabel) :: ldos = 'total_localdos'

    !> total bond currents
    character(lenLabel) :: localCurrents = 'local_currents'

    !> total internal energy
    character(lenLabel) :: egyTotal   = 'total_energy'

    !> total internal energy extrapolated to 0 K
    character(lenLabel) :: egy0Total   = 'extrapolated0_energy'

    !> Energy, which if differentiated gives - force
    character(lenLabel) :: egyForceRelated = 'forcerelated_energy'

    !> Internal electric field
    character(lenLabel) :: internField = 'internal_efield'

    !> External electric field
    character(lenLabel) :: externField = 'external_efield'

    !> two-electron addition/removal energies in ppRPA formalism
    character(lenLabel) :: egyppRPA = '2e_add-rem_energies'

  end type TTagLabelsEnum


  !> Enum containing the tag labels used.
  type(TTagLabelsEnum), parameter :: tagLabels = TTagLabelsEnum()


contains


  !> initialise writer
  subroutine TTaggedWriter_init(this)

    !> Instance
    type(TTaggedWriter), intent(out) :: this

    integer :: nDecDigit, nExpDigit, nChar, nField

    if (this%initialized) then
      return
    end if

    !! "-3.1234567E-123 ": nDec = 7, nExpDigit = 3, nChar = 16
    nExpDigit = ceiling(log(maxexponent(1.0_dp) / log(10.0)) / log(10.0))
    nDecDigit = precision(1.0_dp)
    nChar = nDecDigit + nExpDigit + 6
    nField = 80 / nChar
    if (nField == 0) then
      nField = 1
    end if

    write (this%formReal, "('(', I2.2, 'E', I2.2, '.', I2.2, 'E', I3.3, ')')") nField, nChar,&
        & nDecDigit, nExpDigit

    write (this%formCmplx, "('(', I2.2, '(2E', I2.2, '.', I2.2, 'E', I3.3, '))')") nField / 2,&
        & nChar, nDecDigit, nExpDigit

    nChar = digits(1) + 2
    nField = 80 / nChar
    if (nField == 0) then
      nField = 1
    end if
    write (this%formInt, "('(', I2.2, 'I', I2.2, ')')") nField, nChar
    write (this%formLogical, "('(40L2)')")

    this%initialized = .true.

  end subroutine TTaggedWriter_init


#:for DATA_TYPE, SUFFIX, DATA_TYPE_TAG_NAME, FORMAT_STRING in TEMPLATE_PARAMS
  #:for RANK in range(MAX_RANK + 1)

  !> Write tagged data (data type: ${DATA_TYPE}$)
  subroutine TTaggedWriter_write${SUFFIX}$${RANK}$(this, file, tag, data, optForm)

    !> Instance
    class(TTaggedWriter), intent(inout) :: this

    !> File ID
    integer, intent(in) :: file

    !> tag label
    character(len=*), intent(in) :: tag

    !> data to print
    ${DATA_TYPE}$, intent(in) :: data${FORTRAN_ARG_DIM_SUFFIX(RANK)}$

    !> optional formatting string
    character(len=*), optional, intent(in) :: optForm

    character(len=20) :: form

    @:ASSERT(this%initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(this%${FORMAT_STRING}$)
    end if
    #:if RANK
      call writeTaggedHeader(file, tag, '${DATA_TYPE_TAG_NAME}$', shape(data))
    #:else
      call writeTaggedHeader(file, tag, '${DATA_TYPE_TAG_NAME}$')
    #:endif
    write(file, form) data

  end subroutine TTaggedWriter_write${SUFFIX}$${RANK}$

  #:endfor
#:endfor


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !> Writes the tagged header.
  subroutine writeTaggedHeader(file, tag, dataType, dataShape)

    !> File id to write to
    integer, intent(in) :: file

    !> Tag name
    character(*), intent(in) :: tag

    !> Form string to use
    character(*), intent(in) :: dataType

    !> Original shape of the data
    integer, intent(in), optional :: dataShape(:)

    character(100) :: buffer

    if (present(dataShape)) then
      if (size(dataShape) == 1) then
        write(buffer, "(5A,I0,4A)") '("', getLabel(tag), ":", trim(dataType), ":", size(dataShape),&
            & ":", '",', 'I0', ')'
      else
        write(buffer, "(5A,I0,3A,I0,2A)") '("', getLabel(tag), ":", trim(dataType), ":",&
            & size(dataShape), ":", '",', 'I0,', size(dataShape) - 1, '(",",I0)', ')'
      end if
      write(file, buffer) dataShape
    else
      write(file, "(4A,I0,A)") getLabel(tag), ":", trim(dataType), ":", 0, ":"
    end if

  end subroutine writeTaggedHeader


  !> Extracts the label for a tag
  function getLabel(tag)

    !> relevant tag
    character(len=*), intent(in) :: tag

    !> Label
    character(len=20) :: getLabel

    integer :: lentrim

    lentrim = len_trim(tag)
    if (lentrim >= lenLabel) then
      getLabel(:) = tag(1:lenLabel)
    else
      getLabel(1:lentrim) = tag(1:lentrim)
      getLabel(lentrim+1:lenLabel) = " "
    end if

  end function getLabel


end module dftbp_taggedoutput
