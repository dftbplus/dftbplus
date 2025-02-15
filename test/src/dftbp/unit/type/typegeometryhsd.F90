!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_type_typegeometryhsd
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : AA__Bohr
  use dftbp_extlibs_xmlf90, only : fnode, createDocumentNode, createTextNode, destroyNode
  use dftbp_io_hsdutils, only : setChildValue
  use dftbp_type_typegeometryhsd, only : readTGeometryLammps, TGeometry
  use fortuno_serial, only : suite => serial_suite_item, test_list
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests


  ! Floating point parsing precision
  real(dp), parameter :: prec = 10.0_dp * epsilon(1.0_dp)

  ! New line character
  character(len=1), parameter :: nl = new_line('a')

  ! Node wrapper to trigger automatic destruction when going out of scope
  type :: TLammpsGeoNode
    type(fnode), pointer :: node => null()
  contains
    final :: TLammpsGeoNode_final
  end type TLammpsGeoNode

contains


  $:TEST("minimalCaseWithDefaultValues")
    type(TLammpsGeoNode) :: nw
    type(TGeometry) :: geo
    integer :: ii

    call TLammpsGeoNode_init(nw, "units real",&
        & "3 atoms" // nl //&
        & "1 atom types" // nl //&
        & "Masses" // nl //&
        & "1 1.0" // nl //&
        & "Atoms" // nl //&
        & "1 1" // nl //&
        & "2 1" // nl //&
        & "3 1")
    call readTGeometryLammps(nw%node, geo)
    @:ASSERT(geo%natom == 3)
    @:ASSERT(geo%nSpecies == 1)
    @:ASSERT(geo%speciesNames(1) == 'H')
    @:ASSERT(geo%tPeriodic)
    @:ASSERT(all(abs(geo%origin + 0.5_dp * AA__Bohr) < prec))
    @:ASSERT(all(geo%coords == 0.0_dp))
    do ii = 1, 3
      @:ASSERT(abs(geo%latVecs(ii,ii) - AA__Bohr) < prec)
    end do
  $:END_TEST()


  $:TEST("explicitValuesForAll")
    type(TLammpsGeoNode) :: nw
    type(TGeometry) :: geo
    integer :: ii

    call TLammpsGeoNode_init(nw, "boundary p p p" // nl //&
        & "atom_style atomic" // nl //&
        & "pair_style dftbplus" // nl //&
        & "units real",&
        & "3 atoms" // nl //&
        & "1 atom types" // nl //&
        & "0.0 1.0 xlo xhi" // nl //&
        & "0.0 2.0 ylo yhi" // nl //&
        & "0.0 3.0 zlo zhi" // nl //&
        & "Masses" // nl //&
        & "1 1.0" // nl //&
        & "Atoms" // nl //&
        & "1 1" // nl //&
        & "2 1" // nl //&
        & "3 1")
    call readTGeometryLammps(nw%node, geo)
    @:ASSERT(geo%natom == 3)
    @:ASSERT(geo%nSpecies == 1)
    @:ASSERT(geo%speciesNames(1) == 'H')
    @:ASSERT(geo%tPeriodic)
    @:ASSERT(all(geo%origin == 0.0_dp))
    do ii = 1, 3
      @:ASSERT(abs(geo%latVecs(ii,ii) - real(ii, dp) * AA__Bohr) < prec)
    end do
  $:END_TEST()


  $:TEST("differentOrder")
    type(TLammpsGeoNode) :: nw
    type(TGeometry) :: geo
    integer :: ii

    call TLammpsGeoNode_init(nw, "atom_style atomic" // nl //&
        & "units real" // nl //&
        & "boundary p p p" // nl //&
        & "pair_style dftbplus",&
        & "0.0 3.0 zlo zhi" // nl //&
        & "1 atom types" // nl //&
        & "0.0 2.0 ylo yhi" // nl //&
        & "3 atoms" // nl //&
        & "0.0 1.0 xlo xhi" // nl //&
        & "Atoms" // nl //&
        & "1 1" // nl //&
        & "2 1" // nl //&
        & "3 1" // nl //&
        & "Masses" // nl //&
        & "1 1.0")
    call readTGeometryLammps(nw%node, geo)
    @:ASSERT(geo%natom == 3)
    @:ASSERT(geo%nSpecies == 1)
    @:ASSERT(geo%speciesNames(1) == 'H')
    @:ASSERT(geo%tPeriodic)
    @:ASSERT(all(geo%origin == 0.0))
    @:ASSERT(all(geo%species == 1))
    do ii = 1, 3
      @:ASSERT(abs(geo%latVecs(ii,ii) - real(ii, dp) * AA__Bohr) < prec)
    end do
  $:END_TEST()


  $:TEST("atomsWithCoordinates")
    type(TLammpsGeoNode) :: nw
    type(TGeometry) :: geo
    integer :: ii, jj

    call TLammpsGeoNode_init(nw, "atom_style full" // nl //&
        & "units real",&
        & "3 atoms" // nl //&
        & "1 atom types" // nl //&
        & "Masses" // nl //&
        & "1 1.0" // nl //&
        & "Atoms" // nl //&
        & "1 0 1 99 11.0 12.0 13.0" // nl //&
        & "2 0 1 99 21.0 22.0 23.0" // nl //&
        & "3 0 1 99 31.0 32.0 33.0")
    call readTGeometryLammps(nw%node, geo)
    @:ASSERT(geo%natom == 3)
    @:ASSERT(geo%nSpecies == 1)
    @:ASSERT(geo%speciesNames(1) == 'H')
    @:ASSERT(geo%tPeriodic)
    @:ASSERT(all(geo%species == 1))
    do ii = 1, 3
      do jj = 1, 3
        @:ASSERT(abs(geo%coords(ii,jj) - real(ii + 10*jj, dp) * AA__Bohr) < prec)
      end do
    end do
  $:END_TEST()


  $:TEST("speciesNames")
    type(TLammpsGeoNode) :: nw
    type(TGeometry) :: geo
    integer :: ii

    call TLammpsGeoNode_init(nw, "units real",&
        & "6 atoms" // nl //&
        & "6 atom types" // nl //&
        & "Masses" // nl //&
        & "1 1.0" // nl //&
        & "2 16.0" // nl //&
        & "3 12.0" // nl //&
        & "4 4.0" // nl //&
        & "5 7.0" // nl //&
        & "6 130.0" // nl //&
        & "Atoms" // nl //&
        & "1 1" // nl //&
        & "2 2" // nl //&
        & "3 3" // nl //&
        & "4 4" // nl //&
        & "5 5" // nl //&
        & "6 6")
    call readTGeometryLammps(nw%node, geo)
    @:ASSERT(geo%nSpecies == 6)
    @:ASSERT(geo%speciesNames(1) == 'H')
    @:ASSERT(geo%speciesNames(2) == 'O')
    @:ASSERT(geo%speciesNames(3) == 'C')
    @:ASSERT(geo%speciesNames(4) == 'He')
    @:ASSERT(geo%speciesNames(5) == 'Li')
    @:ASSERT(geo%speciesNames(6) == 'Xe')
    do ii = 1, 6
      @:ASSERT(geo%species(ii) == ii)
    end do
  $:END_TEST()


  $:TEST("strangeFormatting")
    type(TLammpsGeoNode) :: nw
    type(TGeometry) :: geo
    integer :: ii

    call TLammpsGeoNode_init(nw, "# comment" // nl // nl //&
        & "  boundary     p p p #some comment" // nl // nl // nl //&
        & "atom_style atomic   # additional comment" // nl //&
        & "   units real # indeed!" // nl //&
        & " pair_style dftbplus",&
        & "   3  atoms  #another comment" // nl //&
        & "  1 atom    types  " // nl // nl //&
        & "  -1.0     1.0  xlo     xhi" // nl //&
        & "-2   2.0      ylo yhi    " // nl //&
        & " -3.  3.0 zlo zhi" // nl //&
        & "  Masses" // nl // nl //&
        & " 1    1.0" // nl //&
        & "Atoms   # atomic" // nl // nl//&
        & " 1 1 4 4 4 #ignore everything here" // nl //&
        & "  2 1 #ignore me" // nl // nl //&
        & "3 1")
    call readTGeometryLammps(nw%node, geo)
    @:ASSERT(geo%natom == 3)
    @:ASSERT(geo%nSpecies == 1)
    @:ASSERT(geo%tPeriodic)
    @:ASSERT(geo%speciesNames(1) == 'H')
    @:ASSERT(geo%species(1) == 1)
    @:ASSERT(geo%species(2) == 1)
    @:ASSERT(geo%species(3) == 1)
    @:ASSERT(all(geo%species == 1))
    @:ASSERT(all(abs(geo%coords(:,1) - 4.0_dp * AA__Bohr) < prec))
    @:ASSERT(all(geo%coords(:,2) == 0.0_dp))
    @:ASSERT(all(geo%coords(:,3) == 0.0_dp))
    do ii = 1, 3
      @:ASSERT(abs(geo%origin(ii) + real(ii, dp) * AA__Bohr) < prec)
      @:ASSERT(abs(geo%latVecs(ii,ii) - real(2 * ii, dp) * AA__Bohr) < prec)
    end do
  $:END_TEST()


  $:TEST("triclinicCell")
    type(TLammpsGeoNode) :: nw
    type(TGeometry) :: geo

    call TLammpsGeoNode_init(nw, "units real",&
        & "1 atoms" // nl //&
        & "1 atom types" // nl //&
        & "0.0 1.0 xlo xhi" // nl //&
        & "0.0 2.0 ylo yhi" // nl //&
        & "0.0 3.0 zlo zhi" // nl //&
        & "4.0 5.0 6.0 xy xz yz" // nl //&
        & "Masses" // nl //&
        & "1 1.0" // nl //&
        & "Atoms" // nl //&
        & "1 1")
    call readTGeometryLammps(nw%node, geo)
    geo%latVecs = geo%latVecs / AA__Bohr
    @:ASSERT(abs(geo%latVecs(1,1) - 1.0_dp) < prec)
    @:ASSERT(abs(geo%latVecs(2,1)) < prec)
    @:ASSERT(abs(geo%latVecs(3,1)) < prec)
    @:ASSERT(abs(geo%latVecs(1,2) - 4.0_dp) < prec)
    @:ASSERT(abs(geo%latVecs(2,2) - 2.0_dp) < prec)
    @:ASSERT(abs(geo%latVecs(3,2)) < prec)
    @:ASSERT(abs(geo%latVecs(1,3) - 5.0_dp) < prec)
    @:ASSERT(abs(geo%latVecs(2,3) - 6.0_dp) < prec)
    @:ASSERT(abs(geo%latVecs(3,3) - 3.0_dp) < prec)
  $:END_TEST()


  $:TEST("unitConversion")
    type(TLammpsGeoNode) :: nw
    type(TGeometry) :: geo
    integer :: ii

    call TLammpsGeoNode_init(nw, "units si",&
        & "1 atoms" // nl //&
        & "1 atom types" // nl //&
        & "0.0 1.0 xlo xhi" // nl //&
        & "-0.5 0.5 ylo yhi" // nl //&
        & "1.0 2.0 zlo zhi" // nl //&
        & "Masses" // nl //&
        & "1 1.99e-26" // nl //&
        & "Atoms" // nl //&
        & "1 1 1.0e-10 2.0e-10 3.0e-10")
    call readTGeometryLammps(nw%node, geo)
    @:ASSERT(geo%speciesNames(1) == 'C')
    do ii = 1, 3
      @:ASSERT(abs(geo%latVecs(ii,ii) - 1.0e10_dp * AA__Bohr) < prec)
      @:ASSERT(abs(geo%coords(ii,1) - real(ii, dp) * AA__Bohr) < prec)
    end do
  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("typegeometryhsd", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests


  !> Prepare the node for reading the geometry
  subroutine TLammpsGeoNode_init(this, text1, text2)

      !> The geometry node
      type(TLammpsGeoNode), intent(out) :: this

      !> The mock input texts
      character(len=*), intent(in) :: text1, text2

      type(fnode), pointer :: child1, child2

      this%node => createDocumentNode()
      child1 => createTextNode(text1)
      child2 => createTextNode(text2)
      call setChildValue(this%node, "CommandFile", child1)
      call setChildValue(this%node, "DataFile", child2)

  end subroutine TLammpsGeoNode_init


  !> Clean-up the node instance
  subroutine TLammpsGeoNode_final(this)

      !> The geometry node
      type(TLammpsGeoNode), intent(inout) :: this

      if (associated(this%node)) call destroyNode(this%node)

  end subroutine TLammpsGeoNode_final

end module test_type_typegeometryhsd
