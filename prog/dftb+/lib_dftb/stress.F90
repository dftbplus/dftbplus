!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Routines to calculate contributions to the stress tensor
module stress
  use assert
  use accuracy
  use nonscc, only : NonSccDiff
  use scc
  use commontypes
  use slakocont
  use repcont
  implicit none
  private

  public :: getRepulsiveStress, getKineticStress, getNonSCCStress,&
      & getBlockStress, getBlockiStress

contains

  !!* The stress tensor contribution from the repulsive energy term
  !!* @param st stress tensor
  !!* @param coords coordinates (x,y,z, all atoms including possible images)
  !!* @param nNeighbors Number of neighbors for atoms in the central cell
  !!* @param iNeighbors Index of neighbors for a given atom.
  !!* @param species Species of atoms in the central cell.
  !!* @param img2CentCell indexing array for periodic image atoms
  !!* @param repCont Container for repulsive potentials.
  !!* @param cellVol cell volume.
  subroutine getRepulsiveStress(st, coords, nNeighbors, iNeighbors, species, &
      & img2CentCell, repCont, cellVol)
    real(dp), intent(out) :: st(:,:)
    real(dp), intent(in)  :: coords(:,:)
    integer, intent(in)   :: nNeighbors(:)
    integer, intent(in)   :: iNeighbors(0:,:)
    integer, intent(in)   :: species(:)
    integer, intent(in)   :: img2CentCell(:)
    type(ORepCont), intent(in) :: repCont
    real(dp), intent(in)  :: cellVol

    integer :: iAt1, iNeigh, iAt2, iAt2f, ii
    real(dp) :: vect(3), intermed(3), prefac

    @:ASSERT(all(shape(st)==(/ 3, 3/)))

    st(:,:) = 0.0_dp

    do iAt1 = 1, size(nNeighbors)
      do iNeigh = 1, nNeighbors(iAt1)
        iAt2 = iNeighbors(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        vect(:) = coords(:,iAt1) - coords(:,iAt2)
        call getEnergyDeriv(repCont, intermed, vect, species(iAt1), &
            &species(iAt2))
        if (iAt1 == iAt2f) then
          prefac = 0.5_dp
        else
          prefac = 1.0_dp
        end if
        do ii = 1, 3
          st(:, ii) = st(:, ii) - prefac * intermed(:) * vect(ii)
        end do
      end do
    end do
    st = st / cellVol

  end subroutine getRepulsiveStress


  !!* The kinetic contribution to the stress tensor
  !!* @param st stress tensor
  !!* @param mass particle masses
  !!* @param species Species of atoms in the central cell.
  !!* @param velo particle velocities
  !!* @param cellVol cell volume.
  subroutine getKineticStress(st, mass, species, velo, cellVol)
    real(dp), intent(out) :: st(:,:)
    real(dp), intent(in)  :: mass(:)
    integer, intent(in)   :: species(:)
    real(dp), intent(in)  :: velo(:,:)
    real(dp), intent(in)  :: cellVol

    integer :: ii, jj, iAtom, nAtom

    nAtom = size(species)
    @:ASSERT(all(shape(st) == (/3, 3/)))
    @:ASSERT(all(shape(velo) == (/3, nAtom/)))
    @:ASSERT(size(mass) == nAtom)

    st(:,:) = 0.0_dp

    do iAtom = 1, nAtom
      do ii = 1, 3
        do jj = 1, 3
          st(jj,ii) = st(jj,ii) &
              & + mass(iAtom) * velo(jj,iAtom) * velo(ii,iAtom)
        end do
      end do
    end do

    st(:,:) = st(:,:) / cellVol

  end subroutine getKineticStress


  !!* The stress tensor contributions from the non-SCC energy
  !!* @param st stress tensor
  !!* @param DM density matrix in packed format
  !!* @param EDM energy-weighted density matrix in packed format
  !!* @param skHamCont Container for SK Hamiltonian integrals
  !!* @param skOverCont Container for SK overlap integrals
  !!* @param coords list of all atomic coordinates
  !!* @param species list of all atomic species
  !!* @param iNeighbor neighbor list for atoms
  !!* @param nNeighbor number of neighbors of each atom
  !!* @param img2CentCell indexing array for periodic image atoms
  !!* @param iPair indexing array for the Hamiltonian
  !!* @param orb Information about the shells and orbitals in the system.
  !!* @param cellVol cell volume.
  subroutine getNonSCCStress(st, derivator, DM, EDM, skHamCont, skOverCont,&
      & coords, species, iNeighbor, nNeighbor, img2CentCell, iPair, orb,&
      & cellVol)
    real(dp), intent(out) :: st(:,:)
    class(NonSccDiff), intent(in) :: derivator
    real(dp), intent(in) :: DM(:)
    real(dp), intent(in) :: EDM(:)
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iPair(0:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in)  :: cellVol

    integer   :: iOrig, ii, jj
    integer   :: nAtom, iNeigh, iAtom1, iAtom2, iAtom2f
    integer   :: nOrb1, nOrb2
    real(dp)  :: sqrDMTmp(orb%mOrb,orb%mOrb), sqrEDMTmp(orb%mOrb,orb%mOrb)
    real(dp)  :: hPrimeTmp(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp)  :: vect(3), intermed(3)

    @:ASSERT(all(shape(st)==(/3,3/)))
    @:ASSERT(size(DM,dim=1)==size(EDM,dim=1))

    nAtom = size(orb%nOrbAtom)
    st(:,:) = 0.0_dp

    do iAtom1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAtom1)
      !! loop from 1 as no contribution from the atom itself
      do iNeigh = 1, nNeighbor(iAtom1)
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        iOrig = iPair(iNeigh,iAtom1)
        sqrDMTmp(:,:) = 0.0_dp
        sqrEDMTmp(:,:) = 0.0_dp
        hPrimeTmp(:,:,:) = 0.0_dp
        sPrimeTmp(:,:,:) = 0.0_dp
        sqrDMTmp(1:nOrb2,1:nOrb1) = &
            & reshape(DM(iOrig+1:iOrig+nOrb1*nOrb2), (/nOrb2,nOrb1/))
        sqrEDMTmp(1:nOrb2,1:nOrb1) = &
            & reshape(EDM(iOrig+1:iOrig+nOrb1*nOrb2), (/nOrb2,nOrb1/))
        call derivator%getFirstDeriv(hPrimeTmp, skHamCont, coords, species,&
            & iAtom1, iAtom2, orb)
        call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, species,&
              & iAtom1, iAtom2, orb)
        do ii = 1, 3
          ! note factor of 2 for implicit summation over lower triangle of
          ! density matrix:
          intermed(ii) =  2.0_dp * (&
              & sum(sqrDMTmp(1:nOrb2,1:nOrb1)&
              &* hPrimeTmp(1:nOrb2,1:nOrb1,ii)) &
              &- sum(sqrEDMTmp(1:nOrb2,1:nOrb1)&
              &* sPrimeTmp(1:nOrb2,1:nOrb1,ii)) )
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

    st(:,:) = -0.5_dp * st(:,:) / cellVol

  end subroutine getNonSCCStress


  !!* The stress tensor contributions from a potential
  !!* @param st stress tensor
  !!* @param DM density matrix in packed format
  !!* @param EDM energy-weighted density matrix in packed format
  !!* @param skHamCont Container for SK Hamiltonian integrals
  !!* @param skOverCont Container for SK overlap integrals
  !!* @param coords list of all atomic coordinates
  !!* @param species list of all atomic species
  !!* @param iNeighbor neighbor list for atoms
  !!* @param nNeighbor number of neighbors of each atom
  !!* @param nAtom number of real atoms
  !!* @param img2CentCell indexing array for periodic image atoms
  !!* @param iPair indexing array for the Hamiltonian
  !!* @param orb Information about the orbitals
  !!* @param shift block shift from the potential
  !!* @param cellVol cell volume.
  subroutine getBlockStress(st, derivator, DM, EDM, skHamCont, skOverCont,&
      & coords, species, iNeighbor, nNeighbor, img2CentCell, iPair, orb, shift,&
      & cellVol)
    real(dp), intent(out) :: st(:,:)
    class(NonSccDiff), intent(in) :: derivator
    real(dp), intent(in) :: DM(:,:)
    real(dp), intent(in) :: EDM(:)
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iPair(0:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: shift(:,:,:,:)
    real(dp), intent(in)  :: cellVol

    integer   :: iOrig, iSpin, nSpin, ii, jj
    integer   :: nAtom, iNeigh, iAtom1, iAtom2, iAtom2f
    integer   :: nOrb1, nOrb2, iSp1, iSp2
    real(dp)  :: sqrDMTmp(orb%mOrb,orb%mOrb), sqrEDMTmp(orb%mOrb,orb%mOrb)
    real(dp)  :: hPrimeTmp(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp)  :: shiftSprime(orb%mOrb,orb%mOrb)
    real(dp)  :: vect(3), intermed(3)

    nAtom = size(orb%nOrbAtom)
    nSpin = size(shift,dim=4)

    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)
    @:ASSERT(all(shape(st)==(/3,3/)))
    @:ASSERT(size(DM,dim=1)==size(EDM,dim=1))
    @:ASSERT(size(shift,dim=1)==orb%mOrb)
    @:ASSERT(size(shift,dim=2)==orb%mOrb)
    @:ASSERT(size(shift,dim=3)==nAtom)
    @:ASSERT(size(DM,dim=2)==nSpin)

    st(:,:) = 0.0_dp

    do iAtom1 = 1, nAtom
      iSp1 = species(iAtom1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      do iNeigh = 1, nNeighbor(iAtom1)
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        iSp2 = species(iAtom2f)
        if (iAtom1 /= iAtom2) then
          nOrb2 = orb%nOrbSpecies(iSp2)
          iOrig = iPair(iNeigh,iAtom1) + 1
          sqrDMTmp(1:nOrb2,1:nOrb1) = &
              & reshape(DM(iOrig:iOrig+nOrb1*nOrb2-1,1),(/nOrb2,nOrb1/))
          sqrEDMTmp(1:nOrb2,1:nOrb1) = &
              & reshape(EDM(iOrig:iOrig+nOrb1*nOrb2-1),(/nOrb2,nOrb1/))
          call derivator%getFirstDeriv(hPrimeTmp, skHamCont, coords, species,&
              & iAtom1, iAtom2, orb)
          call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, species,&
              & iAtom1, iAtom2, orb)

          intermed(:) = 0.0_dp
          do ii = 1, 3
            ! note factor of 2 for implicit summation over lower triangle of
            ! density matrix:
            intermed(ii) = 2.0_dp * (&
                & sum(sqrDMTmp(1:nOrb2,1:nOrb1)*hPrimeTmp(1:nOrb2,1:nOrb1,ii))&
                &-sum(sqrEDMTmp(1:nOrb2,1:nOrb1)*sPrimeTmp(1:nOrb2,1:nOrb1,ii)))
          end do

          do iSpin = 1, nSpin
            do ii = 1, 3
              shiftSprime(1:nOrb2,1:nOrb1) = 0.5_dp *  ( &
                  & matmul(sPrimeTmp(1:nOrb2,1:nOrb1,ii), &
                  & shift(1:nOrb1,1:nOrb1,iAtom1,iSpin) ) &
                  & + matmul(shift(1:nOrb2,1:nOrb2,iAtom2f,iSpin), &
                  & sPrimeTmp(1:nOrb2,1:nOrb1,ii)) )
              ! again factor of 2 from lower triangle sum of DM
              intermed(ii) = intermed(ii) + 2.0_dp * (&
                  &sum(shiftSprime(1:nOrb2,1:nOrb1) * &
                  &reshape(DM(iOrig:iOrig+nOrb1*nOrb2-1,iSpin),(/nOrb2,nOrb1/))&
                  & ) )
            end do
          end do

          vect(:) = coords(:,iAtom1) - coords(:,iAtom2)
          if (iAtom1/=iAtom2f) then
            do ii = 1, 3
              do jj = 1, 3
                st(jj,ii) = st(jj,ii) &
                    & + 2.0_dp * intermed(jj) * vect(ii)
              end do
            end do
          else
            do ii = 1, 3
              do jj = 1, 3
                st(jj,ii) = st(jj,ii) &
                    & + intermed(jj) * vect(ii)
              end do
            end do
          end if

        end if
      enddo
    enddo

    st(:,:) = -0.5_dp * st(:,:) / cellVol

  end subroutine getBlockStress

  !!* The stress tensor contributions from a complex potential
  !!* @param st stress tensor
  !!* @param DM density matrix in packed format
  !!* @param iDM imaginary part of density matrix in packed format
  !!* @param EDM energy-weighted density matrix in packed format
  !!* @param skHamCont Container for SK Hamiltonian integrals
  !!* @param skOverCont Container for SK overlap integrals
  !!* @param coords list of all atomic coordinates
  !!* @param species list of all atomic species
  !!* @param iNeighbor neighbor list for atoms
  !!* @param nNeighbor number of neighbors of each atom
  !!* @param nAtom number of real atoms
  !!* @param img2CentCell indexing array for periodic image atoms
  !!* @param iPair indexing array for the Hamiltonian
  !!* @param orb Information about the orbitals
  !!* @param shift block shift from the potential
  !!* @param iShift imaginary block shift from the potential
  !!* @param cellVol cell volume.
  subroutine getBlockiStress(st, derivator, DM, iDM, EDM, skHamCont,&
      & skOverCont, coords, species, iNeighbor, nNeighbor, img2CentCell, iPair,&
      & orb, shift, iShift, cellVol)
    real(dp), intent(out) :: st(:,:)
    class(NonSccDiff), intent(in) :: derivator
    real(dp), intent(in) :: DM(:,:)
    real(dp), intent(in) :: iDM(:,:)
    real(dp), intent(in) :: EDM(:)
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iPair(0:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: shift(:,:,:,:)
    real(dp), intent(in) :: iShift(:,:,:,:)
    real(dp), intent(in)  :: cellVol

    integer   :: iOrig, iSpin, nSpin, ii, jj
    integer   :: nAtom, iNeigh, iAtom1, iAtom2, iAtom2f
    integer   :: nOrb1, nOrb2, iSp1, iSp2
    real(dp)  :: sqrDMTmp(orb%mOrb,orb%mOrb), sqrEDMTmp(orb%mOrb,orb%mOrb)
    real(dp)  :: hPrimeTmp(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp)  :: shiftSprime(orb%mOrb,orb%mOrb)
    real(dp)  :: vect(3), intermed(3)

    nAtom = size(orb%nOrbAtom)
    nSpin = size(shift,dim=4)

    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)
    @:ASSERT(all(shape(st)==(/3,3/)))
    @:ASSERT(size(DM,dim=1)==size(EDM,dim=1))
    @:ASSERT(size(shift,dim=1)==orb%mOrb)
    @:ASSERT(size(shift,dim=2)==orb%mOrb)
    @:ASSERT(size(shift,dim=3)==nAtom)
    @:ASSERT(size(DM,dim=2)==nSpin)

    st(:,:) = 0.0_dp

    do iAtom1 = 1, nAtom
      iSp1 = species(iAtom1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      do iNeigh = 1, nNeighbor(iAtom1)
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        iSp2 = species(iAtom2f)
        if (iAtom1 /= iAtom2) then
          nOrb2 = orb%nOrbSpecies(iSp2)
          iOrig = iPair(iNeigh,iAtom1) + 1
          sqrDMTmp(1:nOrb2,1:nOrb1) = &
              & reshape(DM(iOrig:iOrig+nOrb1*nOrb2-1,1),(/nOrb2,nOrb1/))
          sqrEDMTmp(1:nOrb2,1:nOrb1) = &
              & reshape(EDM(iOrig:iOrig+nOrb1*nOrb2-1),(/nOrb2,nOrb1/))
          call derivator%getFirstDeriv(hPrimeTmp, skHamCont, coords, species,&
              & iAtom1, iAtom2, orb)
          call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, species,&
              & iAtom1, iAtom2, orb)

          intermed(:) = 0.0_dp
          do ii = 1, 3
            ! again factor of 2 from lower triangle sum of DM
            intermed(ii) = 2.0_dp * ( &
                & sum(sqrDMTmp(1:nOrb2,1:nOrb1)*hPrimeTmp(1:nOrb2,1:nOrb1,ii))&
                &-sum(sqrEDMTmp(1:nOrb2,1:nOrb1)*sPrimeTmp(1:nOrb2,1:nOrb1,ii))&
                & )
          end do

          do iSpin = 1, nSpin
            do ii = 1, 3
              shiftSprime(1:nOrb2,1:nOrb1) = 0.5_dp *  ( &
                  & matmul(sPrimeTmp(1:nOrb2,1:nOrb1,ii), &
                  & shift(1:nOrb1,1:nOrb1,iAtom1,iSpin) ) &
                  & + matmul(shift(1:nOrb2,1:nOrb2,iAtom2f,iSpin), &
                  & sPrimeTmp(1:nOrb2,1:nOrb1,ii)) )
              ! again factor of 2 from lower triangle sum of DM
              intermed(ii) = intermed(ii) + 2.0_dp * ( &
                  &sum(shiftSprime(1:nOrb2,1:nOrb1) * &
                  &reshape(DM(iOrig:iOrig+nOrb1*nOrb2-1,iSpin),(/nOrb2,nOrb1/))&
                  & ) )
            end do
          end do

          do iSpin = 1, nSpin
            do ii = 1, 3
              shiftSprime(1:nOrb2,1:nOrb1) = 0.5_dp *  ( &
                  & matmul(sPrimeTmp(1:nOrb2,1:nOrb1,ii), &
                  & iShift(1:nOrb1,1:nOrb1,iAtom1,iSpin) ) &
                  & + matmul(iShift(1:nOrb2,1:nOrb2,iAtom2f,iSpin), &
                  & sPrimeTmp(1:nOrb2,1:nOrb1,ii)) )
              intermed(ii) = intermed(ii) + &
                  &sum(shiftSprime(1:nOrb2,1:nOrb1) * &
                  &reshape(iDM(iOrig:iOrig+nOrb1*nOrb2-1,iSpin),(/nOrb2,nOrb1/)) )
            end do
          end do

          vect(:) = coords(:,iAtom1) - coords(:,iAtom2)
          if (iAtom1/=iAtom2f) then
            do ii = 1, 3
              do jj = 1, 3
                st(jj,ii) = st(jj,ii) &
                    & + 2.0_dp * intermed(jj) * vect(ii)
              end do
            end do
          else
            do ii = 1, 3
              do jj = 1, 3
                st(jj,ii) = st(jj,ii) &
                    & + intermed(jj) * vect(ii)
              end do
            end do
          end if
        end if
      enddo
    enddo

    st(:,:) = -0.5_dp * st(:,:) / cellVol

  end subroutine getBlockiStress


end module stress
