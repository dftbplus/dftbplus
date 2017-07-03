!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Module containing various routines for DFTB+U calculations
!!* Intended to be used with SCC switched on !
module dftbplusu
  use assert
  use accuracy
  use message
  use fileid
  use commontypes
  use spin
  implicit none
  private

  public :: shift_DFTBU, AppendBlock_reduce, Block_expand
  public :: E_DFTBU, DFTBplsU_getOrbitalEquiv, DFTBU_blockIndx

  interface shift_DFTBU
    module procedure shift_U
    module procedure shift_iU
  end interface

contains

  !!* Construct the Orbital contribution to the Hamiltonian
  !!* @param H the sparse Hamiltonian to add the +U potential to
  !!* @param rho the density matrix
  !!* @param nNeigh number of surrounding neighbours for each atom
  !!* @param iNeigh list of surrounding neighbours for each atom
  !!* @param species list of the species for each atom
  !!* @param orb  Angular momentum information about the orbitals.
  !!* @param iPair indexing array for the Hamiltonian/density matrix
  !!* @param img2CentCell indexing array to fold from periodic image atom
  !!* numbers back to the central cell number
  !!* @param UJ list of U-J values for each species
  !!* @param nUBlocks number of +U blocks to calculate for each species
  !!* @param nLinBlock number of l values contained in each block for each
  !!* species
  !!* @param iLVals list of l values in each block for each species
  !!* @param functional choice of functional, so far FLL, pSIC (1,2)
  !!* @ref Petukhov, Mazin, Chioncel, and Lichtenstein PHYSICAL REVIEW B
  !!* 67 (15): 153106 APR 15 2003
  !!* @author B. Hourahine
  !!* @todo add other +U potentials
  subroutine shift_U(shift, qBlock, species, orb, functional, &
      & UJ, nUJ, niUJ, iUJ)
    real(dp), intent(inout)        :: shift(:,:,:,:)
    real(dp), intent(in)           :: qBlock(:,:,:,:)
    integer,  intent(in)           :: species(:)
    type(TOrbitals), intent(in)    :: orb
    integer,  intent(in), optional :: functional
    real(dp), intent(in)           :: UJ(:,:)
    integer, intent(in)            :: nUJ(:)
    integer, intent(in)            :: niUJ(:,:)
    integer, intent(in)            :: iUJ(:,:,:)

    integer     :: nAtom, nSpin, iAt, iSp, iSpecies
    integer     :: iFunctional
    integer     :: iStart1, iEnd1, iStart2, iEnd2
    integer     :: ii, jj, kk, ll, ik

    @:ASSERT(all(shape(shift)==shape(qBlock)))
    @:ASSERT(size(shift,dim=1)==orb%mOrb)
    @:ASSERT(size(shift,dim=2)==orb%mOrb)

    nAtom = size(shift,dim=3)
    nSpin = size(shift,dim=4)

    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)

    if (present(functional)) then
      iFunctional = functional
    else
      iFunctional = 1
    end if

    @:ASSERT(iFunctional==1 .or. iFunctional==2)

    if (iFunctional == 1) then
      do iAt = 1, nAtom
        iSpecies = species(iAt)
        do ii = 1, nUJ(iSpecies)
          do jj = 1, niUJ(ii,iSpecies)
            iStart1 = orb%posShell(iUJ(jj,ii,iSpecies),iSpecies)
            iEnd1 = orb%posShell(iUJ(jj,ii,iSpecies)+1,iSpecies)-1
            do kk = iStart1, iEnd1
              shift(kk,kk,iAt,1) = shift(kk,kk,iAt,1) &
                  & + 0.5_dp * UJ(ii,iSpecies)
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
                  shift(ll,kk,iAt,iSp) = shift(ll,kk,iAt,iSp) &
                      & - UJ(ii,iSpecies) * 0.5_dp*qBlock(ll,kk,iAt,iSp)
                  ! factor of 1/2 as using qm not Pauli matrix coefficients
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine Shift_U

  !!* Construct the Orbital contribution to the Hamiltonian
  !!* @param H the sparse Hamiltonian to add the +U potential to
  !!* @param rho the density matrix
  !!* @param nNeigh number of surrounding neighbours for each atom
  !!* @param iNeigh list of surrounding neighbours for each atom
  !!* @param species list of the species for each atom
  !!* @param orb  Angular momentum information about the orbitals.
  !!* @param iPair indexing array for the Hamiltonian/density matrix
  !!* @param img2CentCell indexing array to fold from periodic image atom
  !!* numbers back to the central cell number
  !!* @param UJ list of U-J values for each species
  !!* @param nUJ number of +U blocks to calculate for each species
  !!* @param ninUJ number of l values contained in each block for each
  !!* species
  !!* @param iUJ list of l values in each block for each species
  !!* @param functional choice of functional, so far FLL, pSIC (1,2)
  !!* @ref Petukhov, Mazin, Chioncel, and Lichtenstein PHYSICAL REVIEW B
  !!* 67 (15): 153106 APR 15 2003
  !!* @author B. Hourahine
  !!* @todo add other +U potentials
  subroutine shift_iU(shiftR, shiftI, qBlockR, qBlockI, species, orb, &
      & functional, UJ, nUJ, niUJ, iUJ)
    real(dp), intent(inout)        :: shiftR(:,:,:,:)
    real(dp), intent(inout)        :: shiftI(:,:,:,:)
    real(dp), intent(in)           :: qBlockR(:,:,:,:)
    real(dp), intent(in)           :: qBlockI(:,:,:,:)
    integer,  intent(in)           :: species(:)
    type(TOrbitals), intent(in)    :: orb
    integer,  intent(in), optional :: functional
    real(dp), intent(in)           :: UJ(:,:)
    integer, intent(in)            :: nUJ(:)
    integer, intent(in)            :: niUJ(:,:)
    integer, intent(in)            :: iUJ(:,:,:)

    integer     :: nAtom, nSpin, iAt, iSp, iSpecies
    integer     :: iFunctional
    integer     :: iStart1, iEnd1, iStart2, iEnd2
    integer     :: ii, jj, kk, ll, ik

    @:ASSERT(all(shape(shiftR)==shape(qBlockR)))
    @:ASSERT(all(shape(shiftI)==shape(qBlockI)))
    @:ASSERT(all(shape(shiftR)==shape(shiftI)))
    @:ASSERT(size(shiftR,dim=1)==orb%mOrb)
    @:ASSERT(size(shiftR,dim=2)==orb%mOrb)

    nAtom = size(shiftR,dim=3)
    nSpin = size(shiftR,dim=4)

    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)

    if (present(functional)) then
      iFunctional = functional
    else
      iFunctional = 1
    end if

    @:ASSERT(iFunctional==1 .or. iFunctional==2)

    if (iFunctional == 1) then
      do iAt = 1, nAtom
        iSpecies = species(iAt)
        do ii = 1, nUJ(iSpecies)
          do jj = 1, niUJ(ii,iSpecies)
            iStart1 = orb%posShell(iUJ(jj,ii,iSpecies),iSpecies)
            iEnd1 = orb%posShell(iUJ(jj,ii,iSpecies)+1,iSpecies)-1
            do kk = iStart1, iEnd1
              shiftR(kk,kk,iAt,1) = shiftR(kk,kk,iAt,1) &
                  & + 0.5_dp * UJ(ii,iSpecies)
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
                  shiftR(ll,kk,iAt,iSp) = shiftR(ll,kk,iAt,iSp) &
                      & - UJ(ii,iSpecies) * 0.5_dp*qBlockR(ll,kk,iAt,iSp)
                  shiftI(ll,kk,iAt,iSp) = shiftI(ll,kk,iAt,iSp) &
                      & - UJ(ii,iSpecies) * 0.5_dp*qBlockI(ll,kk,iAt,iSp)
                  ! factor of 1/2 as using qm not Pauli matrix coefficients
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine Shift_iU

  !!* Calculates the energy contribution for the DFTB+U type functionals
  !!* @param egy energy contribution
  !!* @param qBlock charge block populations
  !!* @param species list of the species for each atom
  !!* @param orb  Angular momentum information about the orbitals.
  !!* @param functional choice of functional, so far FLL, pSIC (1,2)
  !!* @param UJ list of U-J values for each species
  !!* @param nUJ number of +U blocks to calculate for each species
  !!* @param ninUJ number of l values contained in each block for each
  !!* species
  !!* @param iUJ list of l values in each block for each species
  !!* @param qiBlock optional skew population for L.S cases
  !!* @note factor of 0.5 in expressions as using double the Pauli spinors
  subroutine E_dftbU(egy, qBlock, species, orb, functional, &
      & UJ, nUJ, niUJ, iUJ, qiBlock)
    real(dp), intent(inout)        :: egy(:)
    real(dp), intent(in)           :: qBlock(:,:,:,:)
    integer,  intent(in)           :: species(:)
    type(TOrbitals), intent(in)    :: orb
    integer,  intent(in), optional :: functional
    real(dp), intent(in)           :: UJ(:,:)
    integer, intent(in)            :: nUJ(:)
    integer, intent(in)            :: niUJ(:,:)
    integer, intent(in)            :: iUJ(:,:,:)
    real(dp), intent(in), optional :: qiBlock(:,:,:,:)

    integer     :: nAtom, nSpin, iAt, iSp, iSpecies
    integer     :: iFunctional
    integer     :: iStart1, iEnd1, iStart2, iEnd2
    integer     :: ii, jj, kk, ll, ik
    real(dp)    :: blockTmp(orb%mOrb,orb%mOrb)

    @:ASSERT(size(qBlock,dim=1)==orb%mOrb)
    @:ASSERT(size(qBlock,dim=2)==orb%mOrb)

    nAtom = size(qBlock,dim=3)
    nSpin = size(qBlock,dim=4)

  #:call ASSERT_CODE
    if (present(qiBlock)) then
      @:ASSERT(all(shape(qiBlock)==shape(qBlock)))
      @:ASSERT(nSpin == 4)
    end if
  #:endcall ASSERT_CODE

    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)
    @:ASSERT(size(egy)==nAtom)

    if (present(functional)) then
      iFunctional = functional
    else
      iFunctional = 1
    end if

    @:ASSERT(iFunctional==1 .or. iFunctional==2)

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
          egy(iAt) = egy(iAt) - 0.25_dp * UJ(ii,iSpecies) *sum(blockTmp(:,:)**2)
          ! extra factor of 1/2 as using qm not Pauli matrix coefficients
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
            egy(iAt) = egy(iAt) - 0.25_dp*UJ(ii,iSpecies) * sum(blockTmp(:,:)**2)
          end do
        end do
      end do
    end if

    ! only trace of the identity (charge) part of the density matrix appears
    ! in this term
    if (iFunctional == 1) then
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


  !!* Returns the equivalence between the orbitals in the DFTB+U interactions
  !!* @param equiv   The equivalence vector on return
  !!* @param orb     Information about the orbitals and their angular momenta
  !!* @param species Species of each atom
  !!* @param  nUJ    How many U-J for each species
  !!* @param  niUJ   number of l-values of U-J for each block
  !!* @param  iUJ    l-values of U-J for each block
  subroutine DFTBplsU_getOrbitalEquiv(equiv, orb, species, nUJ, niUJ, iUJ)
    integer, intent(out) :: equiv(:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)
    integer, intent(in) :: nUJ(:)
    integer, intent(in) :: niUJ(:,:)
    integer, intent(in) :: iUJ(:,:,:)

    integer :: nAtom, iCount, iSpin, nSpin
    integer :: iAt, iSp, ii, jj, kk, iStart, iEnd

    nAtom = size(equiv, dim=2)
    nSpin = size(equiv, dim=3)

    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)
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

  !!* Returns the index for packing the relevant parts of DFTB+U atomic blocks
  !!* into a 1D array
  !!* @param iEqBlockDFTBU The mapping array on return
  !!* @param count Number of prior entries in 1D array holding regular
  !!* charges
  !!* @param orb     Information about the orbitals and their angular momenta
  !!* @param species Species of each atom
  !!* @param  nUJ    How many U-J for each species
  !!* @param  niUJ   number of l-values of U-J for each block
  !!* @param  iUJ    l-values of U-J for each block
  subroutine DFTBU_blockIndx(iEqBlockDFTBU, count, orb, species, nUJ, niUJ, iUJ)
    integer, intent(out)        :: iEqBlockDFTBU(:,:,:,:)
    integer, intent(in)         :: count
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: species(:)
    integer, intent(in)         :: nUJ(:)
    integer, intent(in)         :: niUJ(:,:)
    integer, intent(in)         :: iUJ(:,:,:)

    integer :: nAtom, nSpin, iCount
    integer :: iAt, iSp, iSpecies
    integer :: iStart1, iEnd1, iStart2, iEnd2
    integer :: ii, jj, kk, ll, ik

    nAtom = size(iEqBlockDFTBU, dim=3)
    nSpin = size(iEqBlockDFTBU, dim=4)
    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)

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

  !!* Adds DFTB+U blocks onto end of a 1D vector
  subroutine AppendBlock_reduce(input, equiv, orb, output, skew)
    real(dp), intent(in) :: input(:,:,:,:)
    integer, intent(in) :: equiv(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(inout) :: output(:)
    logical, optional, intent(in) :: skew

    integer :: nAtom, nSpin
    integer :: iS, iOrb1, iOrb2, iAt
    logical :: iSkew

    nAtom = size(input, dim=3)
    nSpin = size(input, dim=4)
    @:ASSERT(size(input, dim=1) == orb%mOrb)
    @:ASSERT(size(input, dim=2) == orb%mOrb)
    @:ASSERT(all(shape(equiv) == (/ orb%mOrb, orb%mOrb, nAtom, nSpin /)))
    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)

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

  subroutine Block_expand(input, blockEquiv, orb, output, &
      & species, nUJ, niUJ, iUJ, orbEquiv, skew)
    real(dp), intent(in)        :: input(:)
    integer, intent(in)         :: blockEquiv(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(out)       :: output(:,:,:,:)
    integer, intent(in)         :: species(:)
    integer, intent(in)         :: nUJ(:)
    integer, intent(in)         :: niUJ(:,:)
    integer, intent(in)         :: iUJ(:,:,:)
    integer, intent(in),optional :: orbEquiv(:,:,:)
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

    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)
    @:ASSERT(size(output, dim=1) == orb%mOrb)
    @:ASSERT(size(output, dim=2) == orb%mOrb)
  #:call ASSERT_CODE
    if (present(orbEquiv)) then
      @:ASSERT(all(shape(orbEquiv) == (/ orb%mOrb, nAtom, nSpin /)))
    end if
  #:endcall ASSERT_CODE
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

end module dftbplusu
