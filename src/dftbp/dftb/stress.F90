!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines to calculate contributions to the stress tensor
module dftbp_dftb_stress
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_schedule, only : distributeRangeWithWorkload, assembleChunks
  use dftbp_dftb_nonscc, only : TNonSccDiff
  use dftbp_dftb_slakocont, only : TSlakoCont
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: getKineticStress, getNonSCCStress, getBlockStress, getBlockiStress

contains

  !> The kinetic contribution to the stress tensor
  subroutine getKineticStress(st, mass, species, velo, cellVol)

    !> stress tensor
    real(dp), intent(out) :: st(:,:)

    !> particle masses
    real(dp), intent(in) :: mass(:)

    !> Species of atoms in the central cell.
    integer, intent(in) :: species(:)

    !> particle velocities
    real(dp), intent(in) :: velo(:,:)

    !> cell volume.
    real(dp), intent(in) :: cellVol

    integer :: ii, jj, iAtom, nAtom

    nAtom = size(species)
    @:ASSERT(all(shape(st) == [3, 3]))
    @:ASSERT(all(shape(velo) == (/3, nAtom/)))
    @:ASSERT(size(mass) == nAtom)

    st(:,:) = 0.0_dp
    do iAtom = 1, nAtom
      do ii = 1, 3
        do jj = 1, 3
          st(jj,ii) = st(jj,ii) + mass(iAtom) * velo(jj,iAtom) * velo(ii,iAtom)
        end do
      end do
    end do

    st(:,:) = st(:,:) / cellVol

  end subroutine getKineticStress


  !> The stress tensor contributions from the non-SCC energy
  subroutine getNonSCCStress(env, st, derivator, DM, EDM, skHamCont, skOverCont, coords, species,&
      & iNeighbour, nNeighbourSK, img2CentCell, iPair, orb, cellVol)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> stress tensor
    real(dp), intent(out) :: st(:,:)

    !> Derivative calculator for (H0,S)
    class(TNonSccDiff), intent(in) :: derivator

    !> density matrix in packed format
    real(dp), intent(in) :: DM(:)

    !> energy-weighted density matrix in packed format
    real(dp), intent(in) :: EDM(:)

    !> Container for SK Hamiltonian integrals
    type(TSlakoCont), intent(in) :: skOverCont

    !> Container for SK overlap integrals
    type(TSlakoCont), intent(in) :: skHamCont

    !> list of all atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> list of all atomic species
    integer, intent(in) :: species(:)

    !> neighbour list for atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> number of neighbours of each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> indexing array for periodic image atoms
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Information about the shells and orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> cell volume.
    real(dp), intent(in) :: cellVol

    integer :: iOrig, ii, jj
    integer :: nAtom, iIter, iNeigh, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: sqrDMTmp(orb%mOrb,orb%mOrb), sqrEDMTmp(orb%mOrb,orb%mOrb)
    real(dp) :: hPrimeTmp(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: vect(3), intermed(3)
    integer, allocatable :: iterIndices(:)

    @:ASSERT(all(shape(st) == [3, 3]))
    @:ASSERT(size(DM,dim=1)==size(EDM,dim=1))

    nAtom = size(orb%nOrbAtom)
    st(:,:) = 0.0_dp

    call distributeRangeWithWorkload(env, 1, nAtom, nNeighbourSK, iterIndices)

    do iIter = 1, size(iterIndices)
      iAtom1 = iterIndices(iIter)
      nOrb1 = orb%nOrbAtom(iAtom1)
      ! loop from 1 as no contribution from the atom itself
      do iNeigh = 1, nNeighbourSK(iAtom1)
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        iOrig = iPair(iNeigh,iAtom1)
        sqrDMTmp(:,:) = 0.0_dp
        sqrEDMTmp(:,:) = 0.0_dp
        hPrimeTmp(:,:,:) = 0.0_dp
        sPrimeTmp(:,:,:) = 0.0_dp
        sqrDMTmp(1:nOrb2,1:nOrb1) = reshape(DM(iOrig+1:iOrig+nOrb1*nOrb2), (/nOrb2,nOrb1/))
        sqrEDMTmp(1:nOrb2,1:nOrb1) = reshape(EDM(iOrig+1:iOrig+nOrb1*nOrb2), (/nOrb2,nOrb1/))
        call derivator%getFirstDeriv(hPrimeTmp, skHamCont, coords, species, iAtom1, iAtom2, orb)
        call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, species, iAtom1, iAtom2, orb)
        do ii = 1, 3
          ! note factor of 2 for implicit summation over lower triangle of density matrix:
          intermed(ii) =  2.0_dp * (sum(sqrDMTmp(1:nOrb2,1:nOrb1)* hPrimeTmp(1:nOrb2,1:nOrb1,ii))&
              & - sum(sqrEDMTmp(1:nOrb2,1:nOrb1) * sPrimeTmp(1:nOrb2,1:nOrb1,ii)) )
        end do
        vect(:) = coords(:,iAtom1) - coords(:,iAtom2)
        if (iAtom1/=iAtom2f) then
          do ii = 1, 3
            do jj = 1, 3
              st(jj,ii) = st(jj,ii) + 2.0_dp*intermed(jj) * vect(ii)
            end do
          end do
        else
          do ii = 1, 3
            do jj = 1, 3
              st(jj,ii) = st(jj,ii) + intermed(jj) * vect(ii)
            end do
          end do
        end if
      end do
    end do

    call assembleChunks(env, st)

    st(:,:) = -0.5_dp * st(:,:) / cellVol

  end subroutine getNonSCCStress


  !> The stress tensor contributions from a potential
  subroutine getBlockStress(env, st, derivator, DM, EDM, skHamCont, skOverCont, coords, species,&
      & iNeighbour, nNeighbourSK, img2CentCell, iPair, orb, shift, cellVol)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> stress tensor
    real(dp), intent(out) :: st(:,:)

    !> density matrix in packed format
    class(TNonSccDiff), intent(in) :: derivator

    !> energy-weighted density matrix in packed format
    real(dp), intent(in) :: DM(:,:)

    !> Container for SK Hamiltonian integrals
    real(dp), intent(in) :: EDM(:)

    !> Container for SK overlap integrals
    type(TSlakoCont), intent(in) :: skHamCont

    !> list of all atomic coordinates
    type(TSlakoCont), intent(in) :: skOverCont

    !> list of all atomic species
    real(dp), intent(in) :: coords(:,:)

    !> neighbour list for atoms
    integer, intent(in) :: species(:)

    !> number of neighbours of each atom
    integer, intent(in) :: iNeighbour(0:,:)

    !> number of real atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> indexing array for periodic image atoms
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Information about the orbitals
    type(TOrbitals), intent(in) :: orb

    !> block shift from the potential
    real(dp), intent(in) :: shift(:,:,:,:)

    !> cell volume.
    real(dp), intent(in) :: cellVol

    integer :: iOrig, iSpin, nSpin, ii, jj
    integer :: nAtom, iIter, iNeigh, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, iSp1, iSp2
    real(dp) :: sqrDMTmp(orb%mOrb,orb%mOrb), sqrEDMTmp(orb%mOrb,orb%mOrb)
    real(dp) :: hPrimeTmp(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: shiftSprime(orb%mOrb,orb%mOrb)
    real(dp) :: vect(3), intermed(3)
    integer, allocatable :: iterIndices(:)

    nAtom = size(orb%nOrbAtom)
    nSpin = size(shift,dim=4)

    @:ASSERT(all(shape(st) == [3, 3]))
    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)
    @:ASSERT(size(DM,dim=1)==size(EDM,dim=1))
    @:ASSERT(size(shift,dim=1)==orb%mOrb)
    @:ASSERT(size(shift,dim=2)==orb%mOrb)
    @:ASSERT(size(shift,dim=3)==nAtom)
    @:ASSERT(size(DM,dim=2)==nSpin)

    st(:,:) = 0.0_dp

    call distributeRangeWithWorkload(env, 1, nAtom, nNeighbourSK, iterIndices)

    do iIter = 1, size(iterIndices)
      iAtom1 = iterIndices(iIter)
      iSp1 = species(iAtom1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      do iNeigh = 1, nNeighbourSK(iAtom1)
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        iSp2 = species(iAtom2f)
        if (iAtom1 /= iAtom2) then
          nOrb2 = orb%nOrbSpecies(iSp2)
          iOrig = iPair(iNeigh,iAtom1) + 1
          sqrDMTmp(1:nOrb2,1:nOrb1) = reshape(DM(iOrig:iOrig+nOrb1*nOrb2-1,1),(/nOrb2,nOrb1/))
          sqrEDMTmp(1:nOrb2,1:nOrb1) = reshape(EDM(iOrig:iOrig+nOrb1*nOrb2-1),(/nOrb2,nOrb1/))
          call derivator%getFirstDeriv(hPrimeTmp, skHamCont, coords, species, iAtom1, iAtom2, orb)
          call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, species, iAtom1, iAtom2, orb)

          intermed(:) = 0.0_dp
          do ii = 1, 3
            ! note factor of 2 for implicit summation over lower triangle of density matrix:
            intermed(ii) = 2.0_dp * (sum(sqrDMTmp(1:nOrb2,1:nOrb1)*hPrimeTmp(1:nOrb2,1:nOrb1,ii))&
                & - sum(sqrEDMTmp(1:nOrb2,1:nOrb1)*sPrimeTmp(1:nOrb2,1:nOrb1,ii)))
          end do

          do iSpin = 1, nSpin
            do ii = 1, 3
              shiftSprime(1:nOrb2,1:nOrb1) = 0.5_dp * ( matmul(sPrimeTmp(1:nOrb2,1:nOrb1,ii),&
                  & shift(1:nOrb1,1:nOrb1,iAtom1,iSpin) )&
                  & + matmul(shift(1:nOrb2,1:nOrb2,iAtom2f,iSpin), sPrimeTmp(1:nOrb2,1:nOrb1,ii)) )
              ! again factor of 2 from lower triangle sum of DM
              intermed(ii) = intermed(ii) + 2.0_dp * (sum(shiftSprime(1:nOrb2,1:nOrb1) *&
                  & reshape(DM(iOrig:iOrig+nOrb1*nOrb2-1,iSpin),(/nOrb2,nOrb1/)) ) )
            end do
          end do

          vect(:) = coords(:,iAtom1) - coords(:,iAtom2)
          if (iAtom1/=iAtom2f) then
            do ii = 1, 3
              do jj = 1, 3
                st(jj,ii) = st(jj,ii) + 2.0_dp * intermed(jj) * vect(ii)
              end do
            end do
          else
            do ii = 1, 3
              do jj = 1, 3
                st(jj,ii) = st(jj,ii) + intermed(jj) * vect(ii)
              end do
            end do
          end if

        end if
      end do
    end do

    call assembleChunks(env, st)

    st(:,:) = -0.5_dp * st(:,:) / cellVol

  end subroutine getBlockStress


  !> The stress tensor contributions from a complex potential
  subroutine getBlockiStress(env, st, derivator, DM, iDM, EDM, skHamCont, skOverCont, coords,&
      & species, iNeighbour, nNeighbourSK, img2CentCell, iPair, orb, shift, iShift, cellVol)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> stress tensor
    real(dp), intent(out) :: st(:,:)

    !> density matrix in packed format
    class(TNonSccDiff), intent(in) :: derivator

    !> imaginary part of density matrix in packed format
    real(dp), intent(in) :: DM(:,:)

    !> energy-weighted density matrix in packed format
    real(dp), intent(in) :: iDM(:,:)

    !> Container for SK Hamiltonian integrals
    real(dp), intent(in) :: EDM(:)

    !> Container for SK overlap integrals
    type(TSlakoCont), intent(in) :: skHamCont

    !> list of all atomic coordinates
    type(TSlakoCont), intent(in) :: skOverCont

    !> list of all atomic species
    real(dp), intent(in) :: coords(:,:)

    !> neighbour list for atoms
    integer, intent(in) :: species(:)

    !> number of neighbours of each atom
    integer, intent(in) :: iNeighbour(0:,:)

    !> number of real atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> indexing array for periodic image atoms
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Information about the orbitals
    type(TOrbitals), intent(in) :: orb

    !> block shift from the potential
    real(dp), intent(in) :: shift(:,:,:,:)

    !> imaginary block shift from the potential
    real(dp), intent(in) :: iShift(:,:,:,:)

    !> cell volume.
    real(dp), intent(in) :: cellVol

    integer :: iOrig, iSpin, nSpin, ii, jj
    integer :: nAtom, iIter, iNeigh, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, iSp1, iSp2
    real(dp) :: sqrDMTmp(orb%mOrb,orb%mOrb), sqrEDMTmp(orb%mOrb,orb%mOrb)
    real(dp) :: hPrimeTmp(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: shiftSprime(orb%mOrb,orb%mOrb)
    real(dp) :: vect(3), intermed(3)
    integer, allocatable :: iterIndices(:)

    nAtom = size(orb%nOrbAtom)
    nSpin = size(shift,dim=4)

    @:ASSERT(all(shape(st) == [3, 3]))
    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)
    @:ASSERT(size(DM,dim=1)==size(EDM,dim=1))
    @:ASSERT(size(shift,dim=1)==orb%mOrb)
    @:ASSERT(size(shift,dim=2)==orb%mOrb)
    @:ASSERT(size(shift,dim=3)==nAtom)
    @:ASSERT(size(DM,dim=2)==nSpin)

    st = 0.0_dp

    call distributeRangeWithWorkload(env, 1, nAtom, nNeighbourSK, iterIndices)

    do iIter = 1, size(iterIndices)
      iAtom1 = iterIndices(iIter)
      iSp1 = species(iAtom1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      do iNeigh = 1, nNeighbourSK(iAtom1)
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        iSp2 = species(iAtom2f)
        if (iAtom1 /= iAtom2) then
          nOrb2 = orb%nOrbSpecies(iSp2)
          iOrig = iPair(iNeigh,iAtom1) + 1
          sqrDMTmp(1:nOrb2,1:nOrb1) = reshape(DM(iOrig:iOrig+nOrb1*nOrb2-1,1),(/nOrb2,nOrb1/))
          sqrEDMTmp(1:nOrb2,1:nOrb1) = reshape(EDM(iOrig:iOrig+nOrb1*nOrb2-1),(/nOrb2,nOrb1/))
          call derivator%getFirstDeriv(hPrimeTmp, skHamCont, coords, species, iAtom1, iAtom2, orb)
          call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, species, iAtom1, iAtom2, orb)

          intermed(:) = 0.0_dp
          do ii = 1, 3
            ! again factor of 2 from lower triangle sum of DM
            intermed(ii) = 2.0_dp * (sum(sqrDMTmp(1:nOrb2,1:nOrb1)*hPrimeTmp(1:nOrb2,1:nOrb1,ii))&
                & - sum(sqrEDMTmp(1:nOrb2,1:nOrb1)*sPrimeTmp(1:nOrb2,1:nOrb1,ii)))
          end do

          do iSpin = 1, nSpin
            do ii = 1, 3
              shiftSprime(1:nOrb2,1:nOrb1) = 0.5_dp *  (matmul(sPrimeTmp(1:nOrb2,1:nOrb1,ii),&
                  & shift(1:nOrb1,1:nOrb1,iAtom1,iSpin) )&
                  & + matmul(shift(1:nOrb2,1:nOrb2,iAtom2f,iSpin), sPrimeTmp(1:nOrb2,1:nOrb1,ii)) )
              ! again factor of 2 from lower triangle sum of DM
              intermed(ii) = intermed(ii) + 2.0_dp * (sum(shiftSprime(1:nOrb2,1:nOrb1) *&
                  & reshape(DM(iOrig:iOrig+nOrb1*nOrb2-1,iSpin),(/nOrb2,nOrb1/)) ) )
            end do
          end do

          do iSpin = 1, nSpin
            do ii = 1, 3
              shiftSprime(1:nOrb2,1:nOrb1) = 0.5_dp *  (matmul(sPrimeTmp(1:nOrb2,1:nOrb1,ii),&
                  & iShift(1:nOrb1,1:nOrb1,iAtom1,iSpin) )&
                  & + matmul(iShift(1:nOrb2,1:nOrb2,iAtom2f,iSpin), sPrimeTmp(1:nOrb2,1:nOrb1,ii)) )
              ! again factor of 2 from lower triangle sum of DM
              intermed(ii) = intermed(ii) + 2.0_dp * sum(shiftSprime(1:nOrb2,1:nOrb1) *&
                  & reshape(iDM(iOrig:iOrig+nOrb1*nOrb2-1,iSpin),(/nOrb2,nOrb1/)) )
            end do
          end do

          vect(:) = coords(:,iAtom1) - coords(:,iAtom2)
          if (iAtom1/=iAtom2f) then
            do ii = 1, 3
              do jj = 1, 3
                st(jj,ii) = st(jj,ii) + 2.0_dp * intermed(jj) * vect(ii)
              end do
            end do
          else
            do ii = 1, 3
              do jj = 1, 3
                st(jj,ii) = st(jj,ii) + intermed(jj) * vect(ii)
              end do
            end do
          end if
        end if
      end do
    end do

    call assembleChunks(env, st)

    st(:,:) = -0.5_dp * st(:,:) / cellVol

  end subroutine getBlockiStress


end module dftbp_dftb_stress
