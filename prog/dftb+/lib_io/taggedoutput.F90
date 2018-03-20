!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to write out various data structures in a comprehensive tagged format.
module taggedoutput
  use assert
  use accuracy, only : dp
  implicit none
  private

  public :: initTaggedWriter, writeTagged


  !> Writes objects in standardized form to the output
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
  end interface writeTagged

  !> Length of permissible tag labels. Tag names (Should be shorter than lenLabel!)
  integer, parameter :: lenLabel = 20

  !> unit cell volume (periodic)
  character(*), parameter, public :: tag_volume = 'cell_volume'

  !> final geometry
  character(*), parameter, public :: tag_endCoord = 'end_coords'

  !> excitation energies in Casida formalism
  character(*), parameter, public :: tag_excEgy = 'exc_energies_sqr'

  !> excited state force contributions
  character(*), parameter, public :: tag_excForce = 'exc_forces'

  !> oscillator strength for excitations
  character(*), parameter, public :: tag_excOsc = 'exc_oscillator'

  !> ground state total forces
  character(*), parameter, public :: tag_forceTot = 'forces'

  !> forces on any external charges present
  character(*), parameter, public :: tag_chrgForces = 'forces_ext_charges'

  !> Gibbs free energy for finite pressure periodic systems
  character(*), parameter, public :: tag_Gibbsfree = 'gibbs_energy'

  !> Gross atomic charges
  character(*), parameter, public :: tag_qOutAtGross  = 'gross_atomic_charges'

  !> numerically calculated second derivatives matrix
  character(*), parameter, public :: tag_HessianNum = 'hessian_numerical'

  !> total energy including electron TS contribution
  character(*), parameter, public :: tag_freeEgy = 'mermin_energy'

  !> Mulliken charges
  character(*), parameter, public :: tag_qOutput = 'orbital_charges'

  !> Pipek-Mezey localisation score of single particle levels
  character(*), parameter, public :: tag_pmlocalise = 'pm_localisation'

  !> total stress tensor for periodic geometries
  character(*), parameter, public :: tag_stressTot = 'stress'

  !> total internal energy
  character(*), parameter, public :: tag_egyTotal   = 'total_energy'

  !> Internal electric field
  character(*), parameter, public :: tag_internfield = 'internal_efield'

  !> External electric field
  character(*), parameter, public :: tag_externfield = 'external_efield'

  ! general format strings

  !> real string format
  character(len=lenLabel) :: formReal

  !> complex  string format
  character(len=lenLabel) :: formCmplx

  !> integer string format
  character(len=lenLabel) :: formInt

  !> logical string format
  character(len=lenLabel) :: formLogical


  !> is the write initialised? required to get relevant machine/compiler constants. Should all be
  !> stored as a derived type
  logical :: initialized = .false.

contains


  !> initialise writer
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


  !> write single real values
  subroutine writeTaggedRealR0(file, tag, value, optForm)

    !> File ID
    integer, intent(in) :: file

    !> tag label
    character(len=*), intent(in) :: tag

    !> value to print
    real(dp), intent(in) :: value

    !> optional formatting string
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


  !> write real vectors
  subroutine writeTaggedRealR1(file, tag, value, optForm)

    !> file ID to write to
    integer, intent(in) :: file

    !> tag label
    character(len=*), intent(in) :: tag

    !> data to write
    real(dp), intent(in) :: value(:)

    !> optional formatting string
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


  !> write real arrays
  subroutine writeTaggedRealR2(file, tag, value, optForm)
    integer, intent(in) :: file
    character(len=*), intent(in) :: tag
    real(dp), intent(in) :: value(:,:)

    !> optional formatting string
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


  !> write 3d real arrays
  subroutine writeTaggedRealR3(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    real(dp), intent(in) :: value(:,:,:)

    !> optional formatting string
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


  !> write 4d real arrays
  subroutine writeTaggedRealR4(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    real(dp), intent(in) :: value(:,:,:,:)

    !> optional formatting string
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


  !> single complex values
  subroutine writeTaggedComplexR0(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    complex(dp), intent(in) :: value

    !> optional formatting string
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


  !> complex vectors
  subroutine writeTaggedComplexR1(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    complex(dp), intent(in) :: value(:)

    !> optional formatting string
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


  !> complex arrays
  subroutine writeTaggedComplexR2(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    complex(dp), intent(in) :: value(:,:)

    !> optional formatting string
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


  !> complex 3d arrays
  subroutine writeTaggedComplexR3(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    complex(dp), intent(in) :: value(:,:,:)

    !> optional formatting string
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


  !> complex 4d arrays
  subroutine writeTaggedComplexR4(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    complex(dp), intent(in) :: value(:,:,:,:)

    !> optional formatting string
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


  !> write integer values
  subroutine writeTaggedIntegerR0(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    integer, intent(in) :: value

    !> optional formatting string
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


  !> write integer vectors
  subroutine writeTaggedIntegerR1(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    integer, intent(in) :: value(:)

    !> optional formatting string
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


  !> write integer arrays
  subroutine writeTaggedIntegerR2(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    integer, intent(in) :: value(:,:)

    !> optional formatting string
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


  !> write 3d integer arrays
  subroutine writeTaggedIntegerR3(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    integer, intent(in) :: value(:,:,:)

    !> optional formatting string
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


  !> write 4d integer arrays
  subroutine writeTaggedIntegerR4(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    integer, intent(in) :: value(:,:,:,:)

    !> optional formatting string
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


  !> write logical values
  subroutine writeTaggedLogicalR0(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    logical, intent(in) :: value

    !> optional formatting string
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


  !> write logical vectors
  subroutine writeTaggedLogicalR1(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    logical, intent(in) :: value(:)

    !> optional formatting string
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


  !> write logical arrays
  subroutine writeTaggedLogicalR2(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    logical, intent(in) :: value(:,:)

    !> optional formatting string
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


  !> write 3d logical arrays
  subroutine writeTaggedLogicalR3(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    logical, intent(in) :: value(:,:,:)

    !> optional formatting string
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


  !> write 4d logical arrays
  subroutine writeTaggedLogicalR4(file, tag, value, optForm)

    !> file id
    integer, intent(in) :: file

    !> tag name
    character(len=*), intent(in) :: tag

    !> data to write
    logical, intent(in) :: value(:,:,:,:)

    !> optional formatting string
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


  !> Extracts the label for a tag
  function getLabel(tag)

    !> relevant tag
    character(len=*), intent(in) :: tag

    !> Label
    character(len=20) :: getLabel

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
