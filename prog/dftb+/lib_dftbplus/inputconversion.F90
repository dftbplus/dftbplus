!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to convert data (typically delivered by the parser) to the internal form.
module dftbp_inputconversion
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_commontypes
  use dftbp_wrappedintr
  use dftbp_periodic, only : buildSquaredAtomIndex
  implicit none
  private

  public :: transformPdosRegionInfo

contains

  !> Transforms the PDOS region as returned by the parser to the internal form.
  subroutine transformPdosRegionInfo(iAtInRegion, tShellResInRegion, regionLabelPrefixes, orb,&
      & specie0, iOrbRegion, regionLabels)

    !> Array of 1D arrays. Each 1D array represents one region. Elements in the 1D array gives the
    !> index of the atoms in the region
    type(TWrappedInt1), intent(in) :: iAtInregion(:)

    !> Whether a region must be treated shell resolved.
    logical, intent(in) :: tShellResInRegion(:)

    !> Label prefix for each region
    character(lc), intent(in) :: regionLabelPrefixes(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Atom type of each atom in the unit cell
    integer, intent(in) :: specie0(:)

    !> Array of 1D arrays, each containing the list of orbitals which must be summed up in the given
    !> region
    type(TWrappedInt1), allocatable, intent(out) :: iOrbRegion(:)

    !> File name for each region
    character(lc), allocatable, intent(out) :: regionLabels(:)

    integer, allocatable :: iAtomStart(:)
    integer :: nRegion, nRegionShell, nAtom
    integer :: iREg, iRegShell

    nRegion = size(iAtInRegion)
    nRegionShell = countRegionsWithShellResolution(iAtInRegion, tShellResInRegion, specie0, orb)
    nAtom = size(orb%nOrbAtom)
    allocate(iOrbRegion(nRegionShell))
    allocate(regionLabels(nRegionShell))
    allocate(iAtomStart(nAtom+1))
    call buildSquaredAtomIndex(iAtomStart, orb)

    iRegShell = 1
    do iReg = 1, nRegion
      if (tShellResInRegion(iReg)) then
        call addShellResolvedRegions(iAtInRegion(iReg)%data, regionLabelPrefixes(iReg), specie0,&
            & orb, iAtomStart, iOrbRegion, regionLabels, iRegShell)
      else
        call addAtomResolvedRegion(iAtInRegion(iReg)%data, regionLabelPrefixes(iReg), orb,&
            & iAtomStart, iOrbRegion, regionLabels, iRegShell)
      end if
    end do

  end subroutine transformPdosRegionInfo


  !> Determines nr. of regions when shell resolution is also considered.
  function countRegionsWithShellResolution(iAtInRegion, tShellResInRegion, specie0, orb)&
      & result(nRegShell)

    !> Array of 1D arrays. Each 1D array represents one region. Elements in the 1D array gives the
    !> index of the atoms in the region
    type(TWrappedInt1), intent(in) :: iAtInRegion(:)

    !> Whether a region must be treated shell resolved.
    logical, intent(in) :: tShellResInRegion(:)

    !> Atom type of each atom in the unit cell
    integer, intent(in) :: specie0(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Resulting count of shells in the region
    integer :: nRegShell

    integer :: iReg, firstAtom

    nRegShell = 0
    do iReg = 1, size(iAtInRegion)
      if (tShellResInRegion(iReg)) then
        firstAtom = iAtInRegion(iReg)%data(1)
        nRegShell = nRegShell + orb%nShell(specie0(firstAtom))
      else
        nRegShell = nRegShell + 1
      end if
    end do

  end function countRegionsWithShellResolution


  !> Adds regions according to the shells of the constituting atoms.
  subroutine addShellResolvedRegions(atomIndices, regionLabelPrefix, specie0, orb, iAtomStart,&
      & iOrbRegion, regionLabels, curReg)

    !> Elements in the 1D array gives the index of the atoms in the region
    integer, intent(in) :: atomIndices(:)

    !> Prefix for the label of the region
    character(*), intent(in) :: regionLabelPrefix

    !> Atom type of each atom in the unit cell
    integer, intent(in) :: specie0(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Start of the atom orbitals in the basis
    integer, intent(in) :: iAtomStart(:)

    !> Array of 1D arrays, each containing the list of orbitals which must be summed up in the given
    !> region
    type(TWrappedInt1), intent(inout) :: iOrbRegion(:)

    !> File name for each region
    character(*), intent(inout) :: regionLabels(:)

    !> Current region being processed
    integer, intent(inout) :: curReg

    integer :: nOrb, nAtom
    integer :: iAt, iSp, iSh
    integer :: ind, ii, jj

    nAtom = size(atomIndices)
    iSp = specie0(atomIndices(1))
    ! Make sure that all atom have the same type within the region.
    @:ASSERT(all(specie0(atomIndices) == iSp))
    do iSh = 1, orb%nShell(iSp)
      nOrb = nAtom * (orb%posShell(iSh + 1, iSp) - orb%posShell(iSh, iSp))
      ind = 1
      allocate(iOrbRegion(curReg)%data(nOrb))
      do ii = 1, nAtom
        iAt = atomIndices(ii)
        do jj = orb%posShell(iSh, iSp), orb%posShell(iSh + 1, iSp) - 1
          iOrbRegion(curReg)%data(ind) = iAtomStart(iAt) + jj - 1
          ind = ind + 1
        end do
      end do
      write(regionLabels(curReg), "(A,A,I0)") trim(regionLabelPrefix), ".", iSh
      curReg = curReg + 1
    end do

  end subroutine addShellResolvedRegions


  !> Adds one region with all the orbitals of the atoms in it.
  subroutine addAtomResolvedRegion(atomIndices, regionLabelPrefix, orb, iAtomStart, iOrbRegion,&
      & regionLabels, curReg)

    !> Elements in the 1D array gives the index of the atoms in the region
    integer, intent(in) :: atomIndices(:)

    !> Prefix for the label of the region
    character(*), intent(in) :: regionLabelPrefix

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Start of the atom orbitals in the basis
    integer, intent(in) :: iAtomStart(:)

    !> Array of 1D arrays, each containing the list of orbitals which must be summed up in the given
    !> region
    type(TWrappedInt1), intent(inout) :: iOrbRegion(:)

    !> File name for each region
    character(*), intent(inout) :: regionLabels(:)

    !> Current region being processed
    integer, intent(inout) :: curReg

    integer :: nOrb, ind, ii, jj, iAt

    nOrb = sum(orb%nOrbAtom(atomIndices))
    allocate(iOrbRegion(curReg)%data(nOrb))
    ind = 1
    do ii = 1, size(atomIndices)
      iAt = atomIndices(ii)
      do jj = 1, orb%nOrbAtom(iAt)
        iOrbRegion(curReg)%data(ind) = iAtomStart(iAt) + jj - 1
        ind = ind + 1
      end do
    end do
    regionLabels(curReg) = regionLabelPrefix
    curReg = curReg + 1

  end subroutine addAtomResolvedRegion


end module dftbp_inputconversion
