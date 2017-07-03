!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Code to calculate forces for several different types of calculation
!!* (non-scc, scc, sDFTB etc)
module forces
  use assert
  use accuracy
  use nonscc, only : NonSccDiff
  use scc
  use commontypes
  use slakocont
  implicit none

  private

  public :: derivative_nonSCC, derivative_shift

  interface derivative_shift
    module procedure derivative_nonSCC !* derivatives without any shift
    module procedure derivative_block  !* derivatives with shift
    module procedure derivative_iBlock !* derivatives with complex shift

  end interface

contains

  !!* The non-SCC electronic force contribution for all atoms, calculated from
  !!* $F_\alpha = \sum_{\mu\nu} \rho_{\mu\nu}\frac{H^0_{\mu\nu}}{\partial
  !!* R_\alpha} - \rho^E_{\mu\xonu}\frac{S_{\mu\nu}}{\partial R_\alpha}$
  !!* with $\rho$ and $\rho^E$ being the density and energy-density
  !!* matrices
  !!*
  !!* @param deriv x,y,z derivatives for each real atom in the system
  !!* @param derivator Differentiatior for the non-scc components.
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
  !!*
  subroutine derivative_nonSCC(deriv, derivator, DM, EDM, skHamCont,&
      & skOverCont, coords, species, iNeighbor, nNeighbor, img2CentCell, iPair,&
      & orb)
    real(dp), intent(out) :: deriv(:,:)
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

    integer   :: iOrig, ii
    integer   :: nAtom, iNeigh, iAtom1, iAtom2, iAtom2f
    integer   :: nOrb1, nOrb2
    real(dp)  :: sqrDMTmp(orb%mOrb,orb%mOrb), sqrEDMTmp(orb%mOrb,orb%mOrb)
    real(dp)  :: hPrimeTmp(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3)

    @:ASSERT(size(deriv,dim=1) == 3)

    nAtom = size(orb%nOrbAtom)
    deriv(:,:) = 0.0_dp

    do iAtom1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAtom1)
      !! loop from 1 as no contribution from the atom itself
      do iNeigh = 1, nNeighbor(iAtom1)
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        if (iAtom1 /= iAtom2f) then
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
          ! note factor of 2 for implicit summation over lower triangle of
          ! density matrix:
          do ii = 1, 3
            deriv(ii,iAtom1) = deriv(ii,iAtom1) &
                &+ sum(sqrDMTmp(1:nOrb2,1:nOrb1)&
                &* 2.0_dp*hPrimeTmp(1:nOrb2,1:nOrb1,ii)) &
                &- sum(sqrEDMTmp(1:nOrb2,1:nOrb1)&
                &* 2.0_dp*sPrimeTmp(1:nOrb2,1:nOrb1,ii))
          end do
          !! Add contribution to the force from atom 1 onto atom 2f using the
          !! symmetry in the blocks, and also the skew symmetry in the
          !! derivatives
          do ii = 1, 3
            deriv(ii,iAtom2f) = deriv(ii,iAtom2f) &
                &- sum(sqrDMTmp(1:nOrb2,1:nOrb1) &
                &* 2.0_dp*hPrimeTmp(1:nOrb2,1:nOrb1,ii)) &
                &+ sum(sqrEDMTmp(1:nOrb2,1:nOrb1)&
                &* 2.0_dp*sPrimeTmp(1:nOrb2,1:nOrb1,ii))
          end do
        end if
      end do
    end do

  end subroutine derivative_nonSCC


  !!* The SCC and spin electronic force contribution for all atoms, calculated
  !!* from
  !!* $F_\alpha = \sum_{\mu\nu} \rho_{\mu\nu}\frac{H^0_{\mu\nu}}{\partial
  !!* R_\alpha} - (\rho^E_{\mu\nu}-\frac{H^1_{\mu\nu}}{S_{\mu\nu}}
  !!* \rho_{\mu\nu})\frac{S_{\mu\nu}}{\partial R_\alpha}$
  !!* with $\rho$ and $\rho^E$ being the density and energy-density
  !!* matrices
  !!*
  !!* @param deriv x,y,z derivatives for each real atom in the system
  !!* @param derivator  Differentiatior for the non-scc components.
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
  !!*
  subroutine derivative_block(deriv, derivator, DM, EDM, skHamCont, skOverCont,&
      & coords, species, iNeighbor, nNeighbor, img2CentCell, iPair, orb, shift)
    real(dp), intent(out) :: deriv(:,:)
    class(NonSccDiff), intent(in) :: derivator
    real(dp), intent(in) :: DM(:,:)
    real(dp), intent(in) :: EDM(:)
    type(OSlakoCont) :: skHamCont, skOverCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iPair(0:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: shift(:,:,:,:)

    integer  :: iOrig, iSpin, ii, nSpin, nAtom
    integer  :: iNeigh, iAtom1, iAtom2, iAtom2f, iSp1, iSp2
    integer  :: nOrb1, nOrb2

    real(dp) :: sqrDMTmp(orb%mOrb,orb%mOrb), sqrEDMTmp(orb%mOrb,orb%mOrb)
    real(dp) :: shiftSprime(orb%mOrb,orb%mOrb)
    real(dp) :: hPrimeTmp(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: derivTmp(3)

    nAtom = size(orb%nOrbAtom)
    nSpin = size(shift,dim=4)
    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin ==4)
    @:ASSERT(size(deriv,dim=1) == 3)
    @:ASSERT(size(deriv,dim=2)==nAtom)
    @:ASSERT(size(DM,dim=1)==size(EDM,dim=1))
    @:ASSERT(size(shift,dim=1)==orb%mOrb)
    @:ASSERT(size(shift,dim=2)==orb%mOrb)
    @:ASSERT(size(shift,dim=3)==nAtom)
    @:ASSERT(size(DM,dim=2)==nSpin)

    deriv(:,:) = 0.0_dp

    do iAtom1 = 1, nAtom
      iSp1 = species(iAtom1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      do iNeigh = 1, nNeighbor(iAtom1)
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        iSp2 = species(iAtom2f)
        if (iAtom1 /= iAtom2f) then
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

          derivTmp(:) = 0.0_dp
          ! note factor of 2 for implicit summation over lower triangle of
          ! density matrix:
          do ii = 1, 3
            derivTmp(ii) = 2.0_dp * (&
                & sum(sqrDMTmp(1:nOrb2,1:nOrb1)*hPrimeTmp(1:nOrb2,1:nOrb1,ii))&
                &-sum(sqrEDMTmp(1:nOrb2,1:nOrb1)*sPrimeTmp(1:nOrb2,1:nOrb1,ii)))
          end do

          do iSpin = 1, nSpin
            do ii = 1, 3
              shiftSprime(1:nOrb2,1:nOrb1) = 0.5_dp * ( &
                  & matmul(sPrimeTmp(1:nOrb2,1:nOrb1,ii), &
                  & shift(1:nOrb1,1:nOrb1,iAtom1,iSpin) ) &
                  & + matmul(shift(1:nOrb2,1:nOrb2,iAtom2f,iSpin), &
                  & sPrimeTmp(1:nOrb2,1:nOrb1,ii)) )
              ! again factor of 2 from lower triangle, cf published force
              ! expressions for SCC:
              derivTmp(ii) = derivTmp(ii) + 2.0_dp * ( &
                  &sum(shiftSprime(1:nOrb2,1:nOrb1) * &
                  &reshape(DM(iOrig:iOrig+nOrb1*nOrb2-1,iSpin),(/nOrb2,nOrb1/))&
                  & ) )
            end do
          end do

          ! forces from atom 1 on atom 2f and 2f onto 1
          deriv(:,iAtom1) = deriv(:,iAtom1) + derivTmp(:)
          deriv(:,iAtom2f) = deriv(:,iAtom2f) - derivTmp(:)

        end if
      enddo
    enddo

  end subroutine derivative_block


  !!* The SCC and spin electronic force contribution for all atoms, calculated
  !!* from
  !!* $F_\alpha = \sum_{\mu\nu} \rho_{\mu\nu}\frac{H^0_{\mu\nu}}{\partial
  !!* R_\alpha} - (\rho^E_{\mu\nu}-\frac{H^1_{\mu\nu}}{S_{\mu\nu}}
  !!* \rho_{\mu\nu})\frac{S_{\mu\nu}}{\partial R_\alpha}$
  !!* with $\rho$ and $\rho^E$ being the density and energy-density
  !!* matrices
  !!*
  !!* @param deriv x,y,z derivatives for each real atom in the system
  !!* @param derivator  Differentiatior for the non-scc components.
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
  !!*
  subroutine derivative_iBlock(deriv, derivator, DM, iDM, EDM, skHamCont,&
      & skOverCont,coords, species, iNeighbor, nNeighbor, img2CentCell, iPair,&
      & orb, shift, iShift)
    real(dp), intent(out) :: deriv(:,:)
    class(NonSccDiff), intent(in) :: derivator
    real(dp), intent(in) :: DM(:,:)
    real(dp), intent(in) :: iDM(:,:)
    real(dp), intent(in) :: EDM(:)
    type(OSlakoCont) :: skHamCont, skOverCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iPair(0:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: shift(:,:,:,:)
    real(dp), intent(in) :: iShift(:,:,:,:)

    integer  :: iOrig, iSpin, ii, nSpin, nAtom
    integer  :: iNeigh, iAtom1, iAtom2, iAtom2f, iSp1, iSp2
    integer  :: nOrb1, nOrb2

    real(dp) :: sqrDMTmp(orb%mOrb,orb%mOrb)
    real(dp)    :: sqrEDMTmp(orb%mOrb,orb%mOrb)
    complex(dp) :: shiftSprime(orb%mOrb,orb%mOrb)
    real(dp)    :: hPrimeTmp(orb%mOrb,orb%mOrb,3),sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp)    :: derivTmp(3)
    complex(dp), parameter :: i = (0.0_dp,1.0_dp)

    nAtom = size(orb%nOrbAtom)
    nSpin = size(shift,dim=4)
    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin ==4)
    @:ASSERT(size(deriv,dim=1) == 3)
    @:ASSERT(size(deriv,dim=2)==nAtom)
    @:ASSERT(size(DM,dim=1)==size(EDM,dim=1))
    @:ASSERT(size(DM,dim=2)==nSpin)
    @:ASSERT(all(shape(iDM)==shape(DM)))
    @:ASSERT(size(shift,dim=1)==orb%mOrb)
    @:ASSERT(size(shift,dim=2)==orb%mOrb)
    @:ASSERT(size(shift,dim=3)==nAtom)
    @:ASSERT(all(shape(iShift)==shape(shift)))

    deriv(:,:) = 0.0_dp

    do iAtom1 = 1, nAtom
      iSp1 = species(iAtom1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      do iNeigh = 1, nNeighbor(iAtom1)
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        iSp2 = species(iAtom2f)
        if (iAtom1 /= iAtom2f) then
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

          derivTmp(:) = 0.0_dp
          ! note factor of 2 for implicit summation over lower triangle of
          ! density matrix:
          do ii = 1, 3
            derivTmp(ii) = 2.0_dp * (&
                & sum(sqrDMTmp(1:nOrb2,1:nOrb1)*hPrimeTmp(1:nOrb2,1:nOrb1,ii))&
                &-sum(sqrEDMTmp(1:nOrb2,1:nOrb1)*sPrimeTmp(1:nOrb2,1:nOrb1,ii)))
          end do

          do iSpin = 1, nSpin
            do ii = 1, 3
              shiftSprime(1:nOrb2,1:nOrb1) = 0.5_dp * ( &
                  & matmul(sPrimeTmp(1:nOrb2,1:nOrb1,ii), &
                  & shift(1:nOrb1,1:nOrb1,iAtom1,iSpin) ) &
                  & + matmul(shift(1:nOrb2,1:nOrb2,iAtom2f,iSpin), &
                  & sPrimeTmp(1:nOrb2,1:nOrb1,ii)) )
              ! again factor of 2 from lower triangle sum of DM
              derivTmp(ii) = derivTmp(ii) &
                  & + 2.0_dp* ( real(sum(shiftSprime(1:nOrb2,1:nOrb1)&
                  & * reshape(DM(iOrig:iOrig+nOrb1*nOrb2-1,iSpin),&
                  & (/nOrb2,nOrb1/)))) )
            end do
          end do

          do iSpin = 1, nSpin
            do ii = 1, 3
              shiftSprime(1:nOrb2,1:nOrb1) = 0.5_dp *  ( &
                  & matmul(sPrimeTmp(1:nOrb2,1:nOrb1,ii), &
                  & ishift(1:nOrb1,1:nOrb1,iAtom1,iSpin) ) &
                  & + matmul(ishift(1:nOrb2,1:nOrb2,iAtom2f,iSpin), &
                  & sPrimeTmp(1:nOrb2,1:nOrb1,ii)) )
              derivTmp(ii) = derivTmp(ii)&
                  & + real(sum(shiftSprime(1:nOrb2,1:nOrb1) *&
                  & reshape(iDM(iOrig:iOrig+nOrb1*nOrb2-1,iSpin),&
                  & (/nOrb2,nOrb1/))))
            end do
          end do

          deriv(:,iAtom1) = deriv(:,iAtom1) + derivTmp(:)
          deriv(:,iAtom2f) = deriv(:,iAtom2f) - derivTmp(:)

        end if
      enddo
    enddo

  end subroutine derivative_iBlock

end module forces
