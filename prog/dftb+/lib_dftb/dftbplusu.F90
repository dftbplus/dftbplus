!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module containing various routines for DFTB+U calculations
!> Intended to be used with SCC switched on !
module dftbp_dftbplusu
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_message
  use dftbp_fileid
  use dftbp_commontypes
  use dftbp_spin
  implicit none
  private

  public :: getDftbUShift, AppendBlock_reduce, Block_expand
  public :: E_DFTBU, DFTBplsU_getOrbitalEquiv, DFTBU_blockIndx
  public :: plusUFunctionals

  
  !> Contains functional characteristics
  type :: TPlusUFuncHelper

    !> Fully localised limit
    integer :: fll = 1

    !> Pseudo self-interaction corrected functional
    integer :: pSic = 2

    !> Possible functional indices
    integer :: indices(2) = [1, 2]

    !> Functional names (to be used in output)
    character(sc) :: names(2) = [character(sc) :: 'FLL', 'pSIC']
    
  end type TPlusUFuncHelper


  !> Can be queried for functional indices and names
  type(TPlusUFuncHelper), parameter :: plusUFunctionals = TPlusUFuncHelper()


  !> Potential shift from LDA+U type potentials
  interface getDftbUShift
    module procedure shift_U
    module procedure shift_iU
  end interface getDftbUShift

contains


  !> Construct the Orbital contribution to the Hamiltonian
  !> Ref: Petukhov, Mazin, Chioncel, and Lichtenstein PHYSICAL REVIEW B 67 (15): 153106 APR 15 2003
  subroutine shift_U(shift, qBlock, species, orb, functional, UJ, nUJ, niUJ, iUJ)

    !> potential to augment
    real(dp), intent(inout) :: shift(:,:,:,:)

    !> block charges
    real(dp), intent(in) :: qBlock(:,:,:,:)

    !> list of the species for each atom
    integer, intent(in) :: species(:)

    !> Angular momentum information about the orbitals.
    type(TOrbitals), intent(in) :: orb

    !> choice of functional, so far FLL, pSIC
    integer, intent(in), optional :: functional

    !> list of U-J values for each species
    real(dp), intent(in) :: UJ(:,:)

    !> number of blocks in each case
    integer, intent(in) :: nUJ(:)

    !> number of shells in each block
    integer, intent(in) :: niUJ(:,:)

    !> shells in the block
    integer, intent(in) :: iUJ(:,:,:)

    integer :: nAtom, nSpin, iAt, iSp, iSpecies
    integer :: iFunctional
    integer :: iStart1, iEnd1, iStart2, iEnd2
    integer :: ii, jj, kk, ll, ik

    @:ASSERT(all(shape(shift)==shape(qBlock)))
    @:ASSERT(size(shift,dim=1)==orb%mOrb)
    @:ASSERT(size(shift,dim=2)==orb%mOrb)

    nAtom = size(shift,dim=3)
    nSpin = size(shift,dim=4)

    if (present(functional)) then
      iFunctional = functional
    else
      iFunctional = plusUFunctionals%fll
    end if

    @:ASSERT(any(plusUFunctionals%indices == iFunctional))

    if (iFunctional == plusUFunctionals%fll) then
      ! move empty states on affected orbitals upwards
      do iAt = 1, nAtom
        iSpecies = species(iAt)
        do ii = 1, nUJ(iSpecies)
          do jj = 1, niUJ(ii,iSpecies)
            iStart1 = orb%posShell(iUJ(jj,ii,iSpecies),iSpecies)
            iEnd1 = orb%posShell(iUJ(jj,ii,iSpecies)+1,iSpecies)-1
            do kk = iStart1, iEnd1
              shift(kk,kk,iAt,1) = shift(kk,kk,iAt,1) + 0.5_dp * UJ(ii,iSpecies)
            end do
          end do
        end do
      end do
    end if

    do iSp = 1, nSpin
      do iAt = 1, nAtom
        iSpecies = species(iAt)
        do ii = 1, nUJ(iSpecies)
          do jj = 1, niUJ(ii,iSpecies)
            iStart1 = orb%posShell(iUJ(jj,ii,iSpecies),iSpecies)
            iEnd1 = orb%posShell(iUJ(jj,ii,iSpecies)+1,iSpecies)-1
            do ik = 1, niUJ(ii,iSpecies)
              iStart2 = orb%posShell(iUJ(ik,ii,iSpecies),iSpecies)
              iEnd2 = orb%posShell(iUJ(ik,ii,iSpecies)+1,iSpecies)-1
              do kk = iStart1, iEnd1
                do ll = iStart2, iEnd2
                  ! factor of 1/2 as using qm not Pauli matrix coefficients
                  shift(ll,kk,iAt,iSp) = shift(ll,kk,iAt,iSp) &
                      & - UJ(ii,iSpecies) * 0.5_dp * qBlock(ll,kk,iAt,iSp)
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine Shift_U


  !> Construct the orbital contribution to the Hamiltonian
  !>
  !> Ref: Petukhov, Mazin, Chioncel, and Lichtenstein Physical Review B 67, 153106 (2003)
  subroutine shift_iU(shiftRe, shiftIm, qBlockR, qBlockI, species, orb, functional, UJ, nUJ, niUJ, &
      & iUJ)

    !> Real part of shift
    real(dp), intent(inout) :: shiftRe(:,:,:,:)

    !> imaginary part of shift
    real(dp), intent(inout) :: shiftIm(:,:,:,:)

    !> real part of block charges
    real(dp), intent(in) :: qBlockR(:,:,:,:)

    !> imaginary part of block charges
    real(dp), intent(in) :: qBlockI(:,:,:,:)

    !> list of the species for each atom
    integer, intent(in) :: species(:)

    !> Angular momentum information about the orbitals.
    type(TOrbitals), intent(in) :: orb

    !> choice of functional, so far FLL, pSIC (1,2)
    integer, intent(in), optional :: functional

    !> list of U-J values for each species
    real(dp), intent(in) :: UJ(:,:)

    !> number of +U blocks to calculate for each species
    integer, intent(in) :: nUJ(:)

    !> number of l values contained in each block for each species
    integer, intent(in) :: niUJ(:,:)

    !> list of l values in each block for each species
    integer, intent(in) :: iUJ(:,:,:)

    integer :: nAtom, nSpin, iAt, iSp, iSpecies
    integer :: iFunctional
    integer :: iStart1, iEnd1, iStart2, iEnd2
    integer :: ii, jj, kk, ll, ik

    @:ASSERT(all(shape(shiftRe)==shape(qBlockR)))
    @:ASSERT(all(shape(shiftIm)==shape(qBlockI)))
    @:ASSERT(all(shape(shiftRe)==shape(shiftIm)))
    @:ASSERT(size(shiftRe,dim=1)==orb%mOrb)
    @:ASSERT(size(shiftRe,dim=2)==orb%mOrb)

    nAtom = size(shiftRe,dim=3)
    nSpin = size(shiftRe,dim=4)

    ! should not get here without spin-orbit present (unless absorbing potentials get added in
    ! future)
    @:ASSERT(nSpin == 4)

    if (present(functional)) then
      iFunctional = functional
    else
      iFunctional = plusUFunctionals%fll
    end if

    @:ASSERT(any(plusUFunctionals%indices == iFunctional))

    if (iFunctional == plusUFunctionals%fll) then
      ! move empty states on affected orbitals upwards
      do iAt = 1, nAtom
        iSpecies = species(iAt)
        do ii = 1, nUJ(iSpecies)
          do jj = 1, niUJ(ii,iSpecies)
            iStart1 = orb%posShell(iUJ(jj,ii,iSpecies),iSpecies)
            iEnd1 = orb%posShell(iUJ(jj,ii,iSpecies)+1,iSpecies)-1
            do kk = iStart1, iEnd1
              shiftRe(kk,kk,iAt,1) = shiftRe(kk,kk,iAt,1) + 0.5_dp * UJ(ii,iSpecies)
            end do
          end do
        end do
      end do
    end if

    do iSp = 1, nSpin
      do iAt = 1, nAtom
        iSpecies = species(iAt)
        do ii = 1, nUJ(iSpecies)
          do jj = 1, niUJ(ii,iSpecies)
            iStart1 = orb%posShell(iUJ(jj,ii,iSpecies),iSpecies)
            iEnd1 = orb%posShell(iUJ(jj,ii,iSpecies)+1,iSpecies)-1
            do ik = 1, niUJ(ii,iSpecies)
              iStart2 = orb%posShell(iUJ(ik,ii,iSpecies),iSpecies)
              iEnd2 = orb%posShell(iUJ(ik,ii,iSpecies)+1,iSpecies)-1
              do kk = iStart1, iEnd1
                do ll = iStart2, iEnd2
                  ! factor of 1/2 as using qm not Pauli matrix coefficients
                  shiftRe(ll,kk,iAt,iSp) = shiftRe(ll,kk,iAt,iSp) &
                      & - UJ(ii,iSpecies) * 0.5_dp * qBlockR(ll,kk,iAt,iSp)
                  shiftIm(ll,kk,iAt,iSp) = shiftIm(ll,kk,iAt,iSp) &
                      & - UJ(ii,iSpecies) * 0.5_dp * qBlockI(ll,kk,iAt,iSp)
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine Shift_iU


  !> Calculates the energy contribution for the DFTB+U type functionals
  !>
  !> Note: factor of 0.5 in expressions as using double the Pauli spinors
  subroutine E_dftbU(egy, qBlock, species, orb, functional, UJ, nUJ, niUJ, iUJ, qiBlock)

    !> energy contribution
    real(dp), intent(out) :: egy(:)

    !> charge block populations
    real(dp), intent(in) :: qBlock(:,:,:,:)

    !> list of the species for each atom
    integer, intent(in) :: species(:)

    !> Angular momentum information about the orbitals.
    type(TOrbitals), intent(in) :: orb

    !> choice of functional, so far FLL, pSIC (1,2)
    integer, intent(in), optional :: functional

    !> list of U-J values for each species
    real(dp), intent(in) :: UJ(:,:)

    !> number of +U blocks to calculate for each species
    integer, intent(in) :: nUJ(:)

    !> number of l values contained in each block for each species
    integer, intent(in) :: niUJ(:,:)

    !> list of l values in each block for each species
    integer, intent(in) :: iUJ(:,:,:)

    !> optional skew population for L.S cases
    real(dp), intent(in), optional :: qiBlock(:,:,:,:)

    integer :: nAtom, nSpin, iAt, iSp, iSpecies
    integer :: iFunctional
    integer :: iStart1, iEnd1, iStart2, iEnd2
    integer :: ii, jj, kk, ll, ik
    real(dp) :: blockTmp(orb%mOrb,orb%mOrb)

    @:ASSERT(size(qBlock,dim=1)==orb%mOrb)
    @:ASSERT(size(qBlock,dim=2)==orb%mOrb)

    nAtom = size(qBlock,dim=3)
    nSpin = size(qBlock,dim=4)

  #:block DEBUG_CODE
    if (present(qiBlock)) then
      @:ASSERT(all(shape(qiBlock)==shape(qBlock)))
      @:ASSERT(nSpin == 4)
    end if
  #:endblock DEBUG_CODE

    @:ASSERT(size(egy)==nAtom)

    egy(:) = 0.0_dp

    if (present(functional)) then
      iFunctional = functional
    else
      iFunctional = plusUFunctionals%fll
    end if

    @:ASSERT(any(plusUFunctionals%indices == iFunctional))

    do iSp = 1, nSpin
      do iAt = 1, nAtom
        iSpecies = species(iAt)
        do ii = 1, nUJ(iSpecies)
          blockTmp(:,:) = 0.0_dp
          do jj = 1, niUJ(ii,iSpecies)
            iStart1 = orb%posShell(iUJ(jj,ii,iSpecies),iSpecies)
            iEnd1 = orb%posShell(iUJ(jj,ii,iSpecies)+1,iSpecies)-1
            do ik = 1, niUJ(ii,iSpecies)
              iStart2 = orb%posShell(iUJ(ik,ii,iSpecies),iSpecies)
              iEnd2 = orb%posShell(iUJ(ik,ii,iSpecies)+1,iSpecies)-1
              do kk = iStart1, iEnd1
                do ll = iStart2, iEnd2
                  blockTmp(ll,kk) = qBlock(ll,kk,iAt,iSp)
                end do
              end do
            end do
          end do
          ! factor of 1/2 as using qm not Pauli matrix coefficients, and another from the +U
          ! functional itself
          egy(iAt) = egy(iAt) - 0.25_dp * UJ(ii,iSpecies) * sum(blockTmp(:,:)**2)
        end do
      end do
    end do

    if (present(qiBlock)) then
      do iSp = 1, nSpin
        do iAt = 1, nAtom
          iSpecies = species(iAt)
          do ii = 1, nUJ(iSpecies)
            blockTmp(:,:) = 0.0_dp
            do jj = 1, niUJ(ii,iSpecies)
              iStart1 = orb%posShell(iUJ(jj,ii,iSpecies),iSpecies)
              iEnd1 = orb%posShell(iUJ(jj,ii,iSpecies)+1,iSpecies)-1
              do ik = 1, niUJ(ii,iSpecies)
                iStart2 = orb%posShell(iUJ(ik,ii,iSpecies),iSpecies)
                iEnd2 = orb%posShell(iUJ(ik,ii,iSpecies)+1,iSpecies)-1
                do kk = iStart1, iEnd1
                  do ll = iStart2, iEnd2
                    blockTmp(ll,kk) = qiBlock(ll,kk,iAt,iSp)
                  end do
                end do
              end do
            end do
            ! factor of 1/2 as using qm not Pauli matrix coefficients, and another from the +U
            ! functional itself
            egy(iAt) = egy(iAt) - 0.25_dp * UJ(ii,iSpecies) * sum(blockTmp(:,:)**2)
          end do
        end do
      end do
    end if

    ! only trace of the identity (charge) part of the density matrix appears in this term
    if (iFunctional == plusUFunctionals%fll) then
      do iAt = 1, nAtom
        iSpecies = species(iAt)
        do ii = 1, nUJ(iSpecies)
          blockTmp(:,:) = 0.0_dp
          do jj = 1, niUJ(ii,iSpecies)
            iStart1 = orb%posShell(iUJ(jj,ii,iSpecies),iSpecies)
            iEnd1 = orb%posShell(iUJ(jj,ii,iSpecies)+1,iSpecies)-1
            do ik = 1, niUJ(ii,iSpecies)
              iStart2 = orb%posShell(iUJ(ik,ii,iSpecies),iSpecies)
              iEnd2 = orb%posShell(iUJ(ik,ii,iSpecies)+1,iSpecies)-1
              do kk = iStart1, iEnd1
                do ll = iStart2, iEnd2
                  blockTmp(ll,kk) = qBlock(ll,kk,iAt,1)
                end do
              end do
            end do
          end do
          do jj = 1, orb%mOrb
            egy(iAt) = egy(iAt) + 0.5_dp * UJ(ii,iSpecies) * blockTmp(jj,jj)
          end do
        end do
      end do
    end if

  end subroutine E_dftbU


  !> Returns the equivalence between the orbitals in the DFTB+U interactions
  subroutine DFTBplsU_getOrbitalEquiv(equiv, orb, species, nUJ, niUJ, iUJ)

    !> The equivalence vector on return
    integer, intent(out) :: equiv(:,:,:)

    !> Information about the orbitals and their angular momenta
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> How many U-J for each species
    integer, intent(in) :: nUJ(:)

    !> number of l-values of U-J for each block
    integer, intent(in) :: niUJ(:,:)

    !> l-values of U-J for each block
    integer, intent(in) :: iUJ(:,:,:)

    integer :: nAtom, iCount, iSpin, nSpin
    integer :: iAt, iSp, ii, jj, kk, iStart, iEnd

    nAtom = size(equiv, dim=2)
    nSpin = size(equiv, dim=3)

    @:ASSERT(size(equiv, dim=1) == orb%mOrb)
    @:ASSERT(size(nUJ) == maxval(species))
    @:ASSERT(all(nUJ <= orb%mShell))
    @:ASSERT(size(niUJ,dim=2) == maxval(species))
    @:ASSERT(all(niUJ <= orb%mShell))
    @:ASSERT(size(iUJ,dim=3) == maxval(species))
    @:ASSERT(size(iUJ,dim=1) <= orb%mShell)
    @:ASSERT(all(iUJ <= orb%mShell))

    equiv(:,:,:) = 0

    ! set all atoms to be initially equivalent to themselves
    do iSpin = 1, nSpin
      do iAt = 1, nAtom
        iSp = species(iAt)
        equiv(1:orb%nOrbSpecies(iSp), iAt, iSpin) = iAt + (iSpin-1)*nAtom
      end do
    end do

    iCount = nSpin*nAtom
    ! set LDA+U blocks to be full of unique orbitals
    do iSpin = 1, nSpin
      do iAt = 1, nAtom
        iSp = species(iAt)
        do ii = 1, nUJ(iSp)
          do jj = 1, niUJ(ii,iSp)
            iStart = orb%posShell(iUJ(jj,ii,iSp),iSp)
            iEnd = orb%posShell(iUJ(jj,ii,iSp)+1,iSp)-1
            do kk = iStart, iEnd
              iCount = iCount + 1
              equiv(kk, iAt, iSpin) = iCount
            end do
          end do
        end do
      end do
    end do

  end subroutine DFTBplsU_getOrbitalEquiv


  !> Returns the index for packing the relevant parts of DFTB+U atomic blocks into a 1D array
  subroutine DFTBU_blockIndx(iEqBlockDFTBU, count, orb, species, nUJ, niUJ, iUJ)

    !> The mapping array on return
    integer, intent(out) :: iEqBlockDFTBU(:,:,:,:)

    !> Number of prior entries in 1D array holding regular charges
    integer, intent(in) :: count

    !> Information about the orbitals and their angular momenta
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> How many U-J for each species
    integer, intent(in) :: nUJ(:)

    !> number of l-values of U-J for each block
    integer, intent(in) :: niUJ(:,:)

    !> l-values of U-J for each block
    integer, intent(in) :: iUJ(:,:,:)

    integer :: nAtom, nSpin, iCount
    integer :: iAt, iSp, iSpecies
    integer :: iStart1, iEnd1, iStart2, iEnd2
    integer :: ii, jj, kk, ll, ik

    nAtom = size(iEqBlockDFTBU, dim=3)
    nSpin = size(iEqBlockDFTBU, dim=4)

    @:ASSERT(size(iEqBlockDFTBU, dim=1) == orb%mOrb)
    @:ASSERT(size(iEqBlockDFTBU, dim=2) == orb%mOrb)
    @:ASSERT(size(nUJ) == maxval(species))
    @:ASSERT(all(nUJ <= orb%mShell))
    @:ASSERT(size(niUJ,dim=2) == maxval(species))
    @:ASSERT(all(niUJ <= orb%mShell))
    @:ASSERT(size(iUJ,dim=3) == maxval(species))
    @:ASSERT(size(iUJ,dim=1) <= orb%mShell)
    @:ASSERT(all(iUJ <= orb%mShell))

    iEqBlockDFTBU(:,:,:,:) = 0

    iCount = count
    do iSp = 1, nSpin
      do iAt = 1, nAtom
        iSpecies = species(iAt)
        do ii = 1, nUJ(iSpecies)
          do jj = 1, niUJ(ii,iSpecies)
            iStart1 = orb%posShell(iUJ(jj,ii,iSpecies),iSpecies)
            iEnd1 = orb%posShell(iUJ(jj,ii,iSpecies)+1,iSpecies)-1
            do ik = 1, niUJ(ii,iSpecies)
              iStart2 = orb%posShell(iUJ(ik,ii,iSpecies),iSpecies)
              iEnd2 = orb%posShell(iUJ(ik,ii,iSpecies)+1,iSpecies)-1
              do kk = iStart1, iEnd1
                do ll = iStart2, iEnd2
                  if (ll > kk) then
                    iCount = iCount + 1
                    iEqBlockDFTBU(ll,kk,iAt,iSp) = iCount
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine DFTBU_blockIndx


  !> Adds DFTB+U blocks onto end of a 1D vector
  subroutine AppendBlock_reduce(input, equiv, orb, output, skew)

    !> unpacked data
    real(dp), intent(in) :: input(:,:,:,:)

    !> equivalences
    integer, intent(in) :: equiv(:,:,:,:)

    !> Information about the orbitals and their angular momenta
    type(TOrbitals), intent(in) :: orb

    !> 1D array with appended data
    real(dp), intent(inout) :: output(:)

    !> is skew symmetry required
    logical, optional, intent(in) :: skew

    integer :: nAtom, nSpin
    integer :: iS, iOrb1, iOrb2, iAt
    logical :: iSkew

    nAtom = size(input, dim=3)
    nSpin = size(input, dim=4)
    @:ASSERT(size(input, dim=1) == orb%mOrb)
    @:ASSERT(size(input, dim=2) == orb%mOrb)
    @:ASSERT(all(shape(equiv) == (/ orb%mOrb, orb%mOrb, nAtom, nSpin /)))

    if (present(skew)) then
      iSkew = skew
    else
      iSkew = .false.
    end if

    do iS = 1, nSpin
      do iAt = 1, nAtom
        do iOrb1 = 1, orb%nOrbAtom(iAt)
          do iOrb2 = 1, orb%nOrbAtom(iAt)
            if (equiv(iOrb1, iOrb2, iAt, iS) > 0) then
              if (iSkew) then
                output(equiv(iOrb1, iOrb2, iAt, iS)) = &
                    & 0.5_dp*( input(iOrb1, iOrb2, iAt, iS) &
                    &  - input(iOrb2, iOrb1, iAt, iS) )
              else
                output(equiv(iOrb1, iOrb2, iAt, iS)) = &
                    & 0.5_dp*( input(iOrb1, iOrb2, iAt, iS) &
                    &  + input(iOrb2, iOrb1, iAt, iS) )
              end if
            end if
          end do
        end do
      end do
    end do

  end subroutine AppendBlock_reduce


  !> Extract DFTB+U blocks from the end of a 1D vector
  subroutine Block_expand(input, blockEquiv, orb, output, species, nUJ, niUJ, iUJ, orbEquiv, skew)

    !> 1D array of packed data
    real(dp), intent(in) :: input(:)

    !> equivalences for blocks on atomic sites
    integer, intent(in) :: blockEquiv(:,:,:,:)

    !> Information about the orbitals and their angular momenta
    type(TOrbitals), intent(in) :: orb

    !> unpacked data
    real(dp), intent(out) :: output(:,:,:,:)

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> How many U-J for each species
    integer, intent(in) :: nUJ(:)

    !> number of l-values of U-J for each block
    integer, intent(in) :: niUJ(:,:)

    !> l-values of U-J for each block
    integer, intent(in) :: iUJ(:,:,:)

    !> equivalences for atoms
    integer, intent(in),optional :: orbEquiv(:,:,:)

    !> is skew symmetry required
    logical, optional, intent(in) :: skew

    integer :: nAtom, nSpin
    integer :: iAt, iSp, iSpecies
    integer :: iStart1, iEnd1, iStart2, iEnd2
    integer :: ii, jj, kk, ll, ik
    logical :: iSkew

    nAtom = size(output, dim=3)
    nSpin = size(output, dim=4)

    if (present(skew)) then
      iSkew = skew
    else
      iSkew = .false.
    end if

    @:ASSERT(size(output, dim=1) == orb%mOrb)
    @:ASSERT(size(output, dim=2) == orb%mOrb)
  #:block DEBUG_CODE
    if (present(orbEquiv)) then
      @:ASSERT(all(shape(orbEquiv) == (/ orb%mOrb, nAtom, nSpin /)))
    end if
  #:endblock DEBUG_CODE
    @:ASSERT(all(shape(blockEquiv) == shape(output)))

    output = 0.0_dp

    do iSp = 1, nSpin
      do iAt = 1, nAtom
        iSpecies = species(iAt)
        do ii = 1, nUJ(iSpecies)
          do jj = 1, niUJ(ii,iSpecies)
            iStart1 = orb%posShell(iUJ(jj,ii,iSpecies),iSpecies)
            iEnd1 = orb%posShell(iUJ(jj,ii,iSpecies)+1,iSpecies)-1
            do ik = 1, niUJ(ii,iSpecies)
              iStart2 = orb%posShell(iUJ(ik,ii,iSpecies),iSpecies)
              iEnd2 = orb%posShell(iUJ(ik,ii,iSpecies)+1,iSpecies)-1
              do kk = iStart1, iEnd1
                if (present(orbEquiv)) then
                  output(kk,kk,iAt,iSp) = input(orbEquiv(kk,iAt,iSp))
                end if
                do ll = iStart2, iEnd2
                  if (ll > kk) then
                    output(ll,kk,iAt,iSp) = input(blockEquiv(ll,kk,iAt,iSp))
                    if (iSkew) then
                      output(kk,ll,iAt,iSp) = -input(blockEquiv(ll,kk,iAt,iSp))
                    else
                      output(kk,ll,iAt,iSp) = input(blockEquiv(ll,kk,iAt,iSp))
                    end if
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine Block_expand

end module dftbp_dftbplusu
