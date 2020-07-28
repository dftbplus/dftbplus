!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Calculates the first and second derivatives of matrix elements
program skderivs
  use dftbp_assert
  use dftbp_globalenv, only : stdOut
  use dftbp_accuracy
  use dftbp_constants
  use dftbp_message
  use xmlf90_flib_dom
  use dftbp_hsdparser, only : parseHSD, dumpHSD, getNodeHSDName
  use dftbp_hsdutils
  use dftbp_hsdutils2
  use dftbp_charmanip
  use dftbp_linkedlist
  use dftbp_slakoeqgrid
  use dftbp_oldskdata
  use dftbp_fileid
#:if WITH_MPI
  use dftbp_mpienv
#:endif
  implicit none


  !> Contains the data necessary for the main program
  type TInputData
    type(TSlakoEqGrid), pointer :: skHam, skOver
    integer, allocatable :: iHam(:), iOver(:)
    real(dp) :: from, to, step
    logical :: value, first, second
    real(dp) :: displ
    character(lc) :: output
  end type TInputData


  !> input data for the calculation of the derivatives
  type(TInputData) :: inp

#:if WITH_MPI
  !> MPI environment, if compiled with mpifort
  type(TMpiEnv) :: mpi

  ! As this is serial code, trap for run time execution on more than 1 processor with an mpi enabled
  ! build
  call TMpiEnv_init(mpi)
  call mpi%mpiSerialEnv()
#:endif

  call parseHSDInput(inp, "skderivs_in.hsd", "skderivs_in")
  call main(inp)

contains


  !> Does the main job of calculating derivatives and writing them to disc
  subroutine main(inp)

    !> instance
    type(TInputData), intent(inout) :: inp

    real(dp), allocatable :: sk(:,:), ham(:,:), over(:,:)
    integer, allocatable :: fpHam(:), fpOver(:)
    character(lc) :: strTmp
    type(string) :: buffer
    integer :: ii, jj, nGrid
    real(dp) :: rr

    nGrid = floor((inp%to - inp%from) / inp%step) + 1
    allocate(sk(getNIntegrals(inp%skHam), -1:1))
    allocate(ham(getNIntegrals(inp%skHam), 0:2))
    allocate(over(getNIntegrals(inp%skOver), 0:2))
    allocate(fpHam(size(inp%iHam)))
    allocate(fpOver(size(inp%iOver)))

    write(stdout, "(A)") ""
    write(stdout, "(A)") "Following files will be created:"
    call resize_string(buffer, 1024)
    do ii = 1, size(inp%iHam)
      fpHam(ii) = getFileId()
      strTmp = trim(inp%output) // ".ham." // i2c(ii)
      open(fpHam(ii), file=strTmp, status="replace", position="rewind")
      write(stdout, "(2X,A)") trim(strTmp)
    end do
    do ii = 1, size(inp%iOver)
      fpOver(ii) = getFileId()
      strTmp = trim(inp%output) // ".ovr." // i2c(ii)
      open(fpOver(ii), file=strTmp, status="replace", position="rewind")
      write(stdout, "(2X,A)") trim(strTmp)
    end do

    do ii = 1, nGrid
      !! Calculate and write value, first and second derivatives
      rr = inp%from + real(ii-1, dp) * inp%step
      call getSKIntegrals(inp%skHam, sk(:,0), rr)
      call getSKIntegrals(inp%skHam, sk(:,-1), rr - inp%displ)
      call getSKIntegrals(inp%skHam, sk(:,1), rr + inp%displ)
      ham(:, 0) = sk(:,0)
      ham(:, 1) = (sk(:,1) - sk(:,-1)) / (2.0_dp * inp%displ)
      ham(:, 2) = (sk(:,1) + sk(:,-1) - 2.0_dp * sk(:,0)) / (inp%displ**2)
      call getSKIntegrals(inp%skOver, sk(:,0), rr)
      call getSKIntegrals(inp%skOver, sk(:,-1), rr - inp%displ)
      call getSKIntegrals(inp%skOver, sk(:,1), rr + inp%displ)
      over(:, 0) = sk(:,0)
      over(:, 1) = (sk(:,1) - sk(:,-1)) / (2.0_dp * inp%displ)
      over(:, 2) = (sk(:,1) + sk(:,-1) - 2.0_dp * sk(:,0)) / (inp%displ**2)

      !! Write out those integrals asked by the user
      do jj = 1, size(inp%iHam)
        write (strTmp, "(E23.15)") rr
        buffer = trim(strTmp)
        if (inp%value) then
          write (strTmp, "(E23.15)") ham(inp%iHam(jj), 0)
          call append_to_string(buffer, trim(strTmp))
        end if
        if (inp%first) then
          write (strTmp, "(E23.15)") ham(inp%iHam(jj), 1)
          call append_to_string(buffer, trim(strTmp))
        end if
        if (inp%second) then
          write (strTmp, "(E23.15)") ham(inp%iHam(jj), 2)
          call append_to_string(buffer, trim(strTmp))
        end if
        write (fpHam(jj), "(A)") char(buffer)
      end do

      do jj = 1, size(inp%iOver)
        write (strTmp, "(E23.15)") rr
        buffer = trim(strTmp)
        if (inp%value) then
          write (strTmp, "(E23.15)") over(inp%iOver(jj), 0)
          call append_to_string(buffer, trim(strTmp))
        end if
        if (inp%first) then
          write (strTmp, "(E23.15)") over(inp%iOver(jj), 1)
          call append_to_string(buffer, trim(strTmp))
        end if
        if (inp%second) then
          write (strTmp, "(E23.15)") over(inp%iOver(jj), 2)
          call append_to_string(buffer, trim(strTmp))
        end if
        write (fpOver(jj), "(A)") char(buffer)
      end do
    end do

    do jj = 1, size(inp%iHam)
      close(fpHam(jj))
    end do
    do jj = 1, size(inp%iOver)
      close(fpOver(jj))
    end do

  end subroutine main


  !> Parses the HSD input
  subroutine parseHSDInput(inp, hsdInputName, rootTag)

    !> parsed data
    type(TInputData), intent(out) :: inp

    !> file name for HSD input
    character(*), intent(in) :: hsdInputName

    !> name of the tag at the root of the tree
    character(*), intent(in) :: rootTag

    type(fNode), pointer :: hsdTree, root, child
    type(TOldSKData) :: skData12(1,1), skData21(1,1)
    character(lc) :: strTmp
    logical :: useOldInter
    type(string) :: buffer
    integer :: angShellOrdered(size(shellNames))
    type(TListIntR1) :: angShells(2)
    type(TListInt), allocatable :: lIntTmp
    real(dp), allocatable :: skHam(:,:), skOver(:,:)
    integer :: skInterMeth, nInt, nSpecies
    integer :: ii, jj

    do ii = 1, maxL+1
      angShellOrdered(ii) = ii - 1
    end do

    call parseHSD(rootTag, hsdInputName, root)

    write(stdout, "(A)") repeat("-", 80)
    write(stdout, "(A)") "Interpreting input file '" // hsdInputName // "'"

    call getChild(hsdTree, rootTag, root)

    !! First interaction
    call getChildValue(root, "File1", buffer)
    call readFromFile(skData12(1,1), unquote(char(buffer)), .false.)
    call getChildValue(root, "MaxAngularMomentum1", buffer, child=child)
    strTmp = unquote(char(buffer))
    call init(angShells(1))
    do jj = 1, size(shellNames)
      if (trim(strTmp) == trim(shellNames(jj))) then
        call append(angShells(1), angShellOrdered(:jj))
      end if
    end do
    if (len(angShells(1)) < 1) then
      call detailedError(child, "Invalid orbital name '" // &
          &trim(strTmp) // "'")
    end if

    call getChildValue(root, "File2", buffer, "")
    if (char(buffer) == "") then
      nSpecies = 1
    else
      nSpecies = 2
      call readFromFile(skData21(1,1), unquote(char(buffer)), .false.)
      call getChildValue(root, "MaxAngularMomentum2", buffer, child=child)
      strTmp = unquote(char(buffer))
      call init(angShells(2))
      do jj = 1, size(shellNames)
        if (trim(strTmp) == trim(shellNames(jj))) then
          call append(angShells(2), angShellOrdered(:jj))
        end if
      end do
      if (len(angShells(2)) < 1) then
        call detailedError(child, "Invalid orbital name '" // &
            &trim(strTmp) // "'")
      end if
    end if

    call getChildValue(root, "OldSKInterpolation", useOldInter)
    if (useOldInter) then
      skInterMeth = skEqGridOld
    else
      skInterMeth = skEqGridNew
    end if

    !! Create Slako tables
    nInt = getNSKIntegrals(angShells(1), angShells(nSpecies))
    allocate(skHam(size(skData12(1,1)%skHam, dim=1), nInt))
    allocate(skOver(size(skData12(1,1)%skOver, dim=1), nInt))
    if (nSpecies == 1) then
      call getFullTable(skHam, skOver, skData12, skData12, angShells(1), &
          &angShells(1))
    else
      call getFullTable(skHam, skOver, skData21, skData21, angShells(1), &
          &angShells(2))
    end if
    allocate(inp%skHam)
    allocate(inp%skOver)
    call init(inp%skHam, skData12(1,1)%dist, skHam, skInterMeth)
    call init(inp%skOver, skData12(1,1)%dist, skOver, skInterMeth)

    call getChildValue(root, "Step", inp%step)
    call getChildValue(root, "Value", inp%value)
    call getChildValue(root, "FirstDerivative", inp%first)
    call getChildValue(root, "SecondDerivative", inp%second)
    call getChildValue(root, "Displacement", inp%displ)
    call getChildValue(root, "Start", inp%from, child=child)
    if (inp%from - inp%displ < skData12(1, 1)%dist) then
      write(strTmp, "(A, F8.4)") "With given displacement, start point must be larger than ",&
          & skData12(1, 1)%dist + inp%displ
      call detailedError(child, trim(strTmp))
    end if
    call getChildValue(root, "End", inp%to, child=child)
    call getChildValue(root, "OutputPrefix", buffer)
    inp%output = unquote(char(buffer))

    allocate(lIntTmp)
    call init(lIntTmp)
    call getChildValue(root, "Hamiltonian", lIntTmp, child=child)
    allocate(inp%iHam(len(lIntTmp)))
    call asArray(lIntTmp, inp%iHam)
    if (any(inp%iHam < 1) .or. any(inp%iHam > nInt)) then
      call detailedError(child, "Integral index must be between 1 and " &
          &// i2c(nInt))
    end if
    deallocate(lIntTmp)
    allocate(lIntTmp)
    call init(lIntTmp)
    call getChildValue(root, "Overlap", lIntTmp)
    allocate(inp%iOver(len(lIntTmp)))
    call asArray(lIntTmp, inp%iOver)
    if (any(inp%iOver < 1) .or. any(inp%iover > nInt)) then
      call detailedError(child, "Integral index must be between 1 and " &
          &// i2c(nInt))
    end if

    !! Issue warning about unprocessed nodes
    call warnUnprocessedNodes(root)
    write(stdout, "(A)") "Done."
    write(stdout, "(A)") repeat("-", 80)

  end subroutine parseHSDInput


  !> Creates from the columns of the Slater-Koster files for A-B and B-A a full table for A-B,
  !> containing all integrals.
  subroutine getFullTable(skHam, skOver, skData12, skData21, angShells1, angShells2)

    !> Resulting table of H integrals
    real(dp), intent(out) :: skOver(:,:)

    !> Resulting table of S integrals
    real(dp), intent(out) :: skHam(:,:)

    !> Contains all SK files describing interactions for A-B
    type(TOldSKData), intent(in), target :: skData12(:,:)

    !> Contains all SK files describing interactions for B-A
    type(TOldSKData), intent(in), target :: skData21(:,:)

    !> Angular momenta to pick from the SK-files for species A
    type(TListIntR1), intent(inout) :: angShells1

    !> Angular momenta to pick from the SK-files for species B
    type(TListIntR1), intent(inout) :: angShells2

    integer :: ind, iSK1, iSK2, iSh1, iSh2, nSh1, nSh2, l1, l2, lMin, lMax, mm
    integer :: angShell1(maxL+1), angShell2(maxL+1)
    real(dp), pointer :: pHam(:,:), pOver(:,:)

    !! Maps (mm, l1, l2 ) onto an element in the SK table.
    !! l2 >= l1 (l1 = 0, 1, ...; l2 = 0, 1, ...), m <= l1.
    integer, parameter :: skMap(0:maxL, 0:maxL, 0:maxL) &
        &= reshape((/&
        &20, 0,  0,  0,  19,  0,  0,  0,  18,  0,  0,  0,  17,  0,  0,  0,&
        & 0, 0,  0,  0,  15, 16,  0,  0,  13, 14,  0,  0,  11, 12,  0,  0,&
        & 0, 0,  0,  0,   0,  0,  0,  0,   8,  9, 10,  0,   5,  6,  7,  0,&
        & 0, 0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   1,  2,  3,  4/),&
        &(/maxL + 1, maxL + 1, maxL + 1/))

    ind = 1
    do iSK1 = 1, len(angShells1)
      call intoArray(angShells1, angShell1, nSh1, iSK1)
      do iSh1 = 1, nSh1
        l1 = angShell1(iSh1)
        do iSK2 = 1, len(angShells2)
          call intoArray(angShells2, angShell2, nSh2, iSK2)
          do iSh2 = 1, nSh2
            l2 = angShell2(iSh2)
            if (l1 <= l2) then
              pHam => skData12(iSK2,iSK1)%skHam
              pOver => skData12(iSK2,iSK1)%skOver
              lMin = l1
              lMax = l2
            else
              pHam => skData21(iSK1,iSK2)%skHam
              pOver => skData21(iSK1,iSK2)%skOver
              lMin = l2
              lMax = l1
            end if
            do mm = 0, lMin
              !! Safety check, if array size are appropriate
              @:ASSERT(all(shape(skHam) >= (/ size(pHam, dim=1), ind /)))
              @:ASSERT(all(shape(skOver) >= (/ size(pOver, dim=1), ind /)))
              @:ASSERT(size(pHam, dim=1) == size(pOver, dim=1))
              skHam(:,ind) = pHam(:,skMap(mm,lMax,lMin))
              skOver(:,ind) = pOver(:,skMap(mm,lMax,lMin))
              ind = ind + 1
            end do
          end do
        end do
      end do
    end do

  end subroutine getFullTable


  !> Returns the nr. of Slater-Koster integrals necessary to describe the interactions between two
  !> species.
  function getNSKIntegrals(angShells1, angShells2) result(nInt)

    !> list of shells for species B
    type(TListIntR1), intent(inout) :: angShells2

    !> list of shells for species A
    type(TListIntR1), intent(inout) :: angShells1

    !> count of integrals
    integer :: nInt

    integer :: iSh1, iSh2, nSh1, nSh2, ang1, ang2
    integer :: angs1(maxL+1), angs2(maxL+1)

    call intoArray(angShells1, angs1, nSh1, 1)
    call intoArray(angShells2, angs2, nSh2, 1)
    nInt = 0
    do iSh1 = 1, nSh1
      ang1 = angs1(iSh1)
      do iSh2 = 1, nSh2
        ang2 = angs2(iSh2)
        nInt = nInt + min(ang1, ang2) + 1
      end do
    end do

  end function getNSKIntegrals

end program skderivs
