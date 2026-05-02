!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!


!> Demonstrates API use with derivative responses
!!
!! Use it with the inputs in the test/src/dftbp/api/mm/testcases/deriv_responses/ folder.
!!
program test_responses
  use dftbplus, only : getDftbPlusApi, getDftbPlusBuild, TDftbPlus, TDftbPlus_init, TDftbPlusInput
  use dftbp_common_constants, only : AA__Bohr
  ! Only needed for the internal test system
  use testhelpers, only : writeAutotestTag
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  call main_()

contains


  !! Main test routine
  !!
  !! All non-constant variables must be defined here to ensure that they are all explicitely
  !! deallocated before the program finishes (avoiding residual memory that tools like valgrind
  !! notice).
  !!
  subroutine main_()

    type(TDftbPlus) :: dftbp
    type(TDftbPlusInput) :: input

    real(dp) :: merminEnergy, rNum
    real(dp), allocatable :: latVecs(:,:)
    real(dp), allocatable :: gradients(:, :), grossCharges(:), dqdx(:,:,:), dqdxExt(:,:,:)
    real(dp), allocatable :: dqdxTmp(:,:,:), dqdxExtTmp(:,:,:)
    integer, allocatable :: wrtAtoms(:), wrtCharges(:)

    character(:), allocatable :: DftbVersion
    integer :: major, minor, patch, iAt, jAt, nAtom, nExtCharges, iExt, nDeriv, iDeriv
    logical :: isPeriodic, internalTests(2)
    integer :: iTests

    !integer :: devNull

    call random_seed()

    call getDftbPlusBuild(DftbVersion)
    write(*,*)'DFTB+ build: ' // "'" // trim(DftbVersion) // "'"
    call getDftbPlusApi(major, minor, patch)
    write(*,"(1X,A,1X,I0,'.',I0,'.',I0)")'API version:', major, minor, patch

    ! initialise a calculation then read input from file
    call TDftbPlus_init(dftbp)

    ! Note: setting the global standard output to /dev/null will also suppress run-time error
    ! messages
    !open(newunit=devNull, file="/dev/null", action="write")
    !call TDftbPlus_init(dftbp, outputUnit=devNull)

    call dftbp%getInputFromFile("dftb_in.hsd", input)
    call dftbp%setupCalculator(input)

    ! obtain energy and forces
    call dftbp%getEnergy(merminEnergy)

    internalTests(:) = .true.
    iTests = 0

    nAtom = dftbp%nrOfAtoms()
    allocate(gradients(3,nAtom))
    allocate(grossCharges(nAtom))

    call dftbp%getGradients(gradients)
    call dftbp%getGrossCharges(grossCharges)
    print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
    print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)
    print "(A,3F15.10)", 'Obtained atomic gross charges:'
    do iAt = 1, nAtom
      print "(T22,F15.10)", grossCharges(iAt)
    end do

    call dftbp%getGeometry(latVecs = latVecs)
    isPeriodic = allocated(latVecs)

    if (.not. isPeriodic) then

      print *,'Derivatives testing for', nAtom, 'atoms'

      iTests = iTests + 1

      ! Calculate for all atoms
      allocate(dqdx(nAtom,3,nAtom))

      call dftbp%getChargeDerivatives(dqdx=dqdx)

      print *, "Coupled-perturbed derivatives of atom charges wrt all DFTB atoms (e/AA)"
      do iAt = 1, nAtom
        print "(A,I0)","dx_At", iAt
        do jAt = 1, nAtom
          print "(3F15.10)", dqdx(jAt, :, iAt) * AA__Bohr
        end do
      end do

      ! Now calculate for a random sub-set of atoms and check they match the previous calculation
      call random_number(rNum)
      nDeriv = ceiling(nAtom * rNum) ! 1 .. nAtom
      call randomSelection(wrtAtoms, nDeriv, nAtom)
      allocate(dqdxTmp(nAtom, 3, nDeriv))
      call dftbp%getChargeDerivatives(dqdx=dqdxTmp, wrtAtoms=wrtAtoms)

      do iDeriv = 1, nDeriv
        if (any(abs(dqdxTmp(:,:,iDeriv) - dqdx(:,:,wrtAtoms(iDeriv))) > epsilon(1.0_dp) )) then
          write(*,*)'Mismatch for atom', wrtAtoms(iDeriv)
          internalTests(iTests) = .false.
        else
          write(*,*)'Match for atom', wrtAtoms(iDeriv)
        end if
      end do

      if (.not.internalTests(iTests)) deallocate(dqdx)
      deallocate(wrtAtoms)
      deallocate(dqdxTmp)

    else

      print *, "Skipping atom derivative calculation as periodic"

    end if

    nExtCharges = dftbp%nrOfExternalCharges()

    if (nExtCharges > 0) then
      print *,'Derivatives testing for', nExtCharges, 'external charges'

      iTests = iTests + 1

      allocate(dqdxExt(nAtom,3,nExtCharges))

      call dftbp%getChargeDerivatives(dqdxExt=dqdxExt)

      print "(A)", "Coupled-perturbed derivatives of DFTB charges wrt to all external charge&
          & positions (e/AA)"
      do iExt = 1, nExtCharges
        print "(1X,A,I0)","/dx_Chrg", iExt
        do jAt = 1, nAtom
          print "(3F15.10)", dqdxExt(jAt, :, iExt) * AA__Bohr
        end do
      end do

      ! Now calculate for a random sub-set of atoms and check they match the previous calculation
      call random_number(rNum)
      nDeriv = ceiling(nExtCharges * rNum) ! 1 .. nExtCharges
      call randomSelection(wrtCharges, nDeriv, nExtCharges)
      allocate(dqdxExtTmp(nAtom, 3, nDeriv))
      call dftbp%getChargeDerivatives(dqdxExt=dqdxExtTmp, dxExtCharges=wrtCharges)

      do iDeriv = 1, nDeriv
        if (any(abs(dqdxExtTmp(:,:,iDeriv) - dqdxExt(:,:,wrtCharges(iDeriv))) > epsilon(1.0_dp) ))&
            & then
          write(*,*)'Mismatch for charge', wrtCharges(iDeriv)
          internalTests(iTests) = .false.
        else
          write(*,*)'Match for charge', wrtCharges(iDeriv)

        end if
      end do

      if (.not.internalTests(iTests)) deallocate(dqdxExt)
      deallocate(dqdxExtTmp)
      deallocate(wrtCharges)

    else

      print *, "Skipping external charge derivative calculation as no charges"

    end if

    ! Write file for internal test system
    call writeAutotestTag(merminEnergy=merminEnergy, gradients=gradients,&
        & grossCharges=grossCharges, dqdx=dqdx, dqdxExt=dqdxExt)

  end subroutine main_

  !> Generate a list of m random, non-repeating values in the range 1..n
  subroutine randomSelection(list, m, n)

    !> Resulting list of values
    integer, allocatable, intent(out) :: list(:)

    !> Number of values to select (m) from 1..m range
    integer, intent(in) :: m, n

    integer :: work(n), im, in, iRand
    real(dp) :: rNum

    allocate(list(m))
    do in = 1, n
      work(in) = in
    end do
    do im = 1, m
      call random_number(rNum)
      iRand = ceiling(rNum * (n - im + 1))
      list(im) = work(iRand)
      ! shift rest of numbers down
      if (im < m) work(iRand : n - im) = work(iRand+1 : n - im + 1)
    end do

  end subroutine randomSelection

end program test_responses
