#:include "common.fypp"

!!* Dispersion with Many Body Dispersion method by Tkatchenko et al.
module manybodydisp
  use globalenv, only : stdOut, tIoProc
  use mpifx, only : mpifx_comm
  use accuracy
  use constants, only: pi, Bohr__AA, AA__Bohr, eV__Hartree, Hartree__eV
  use mbd
  use commontypes, only : TOrbitals
  use message, only: error
  implicit none
  private

  public:: TMbdInit, TMbd
  public:: MBDinit, MBDupdateCoords, MBDupdateLatVecs
  public:: MBDgetEnergy, MBDgetGradients, MBDcalculateCPA
  public:: MBDgetStress
  public :: onsiteMullikenPerAtom

contains

  !takes TMbd instance, TMbdInit instance, nAtom, species0,
  !referenceN0, nShell, species_name, latVecs, recVecs
  subroutine MBDinit(this, inp, mympi, nAtom, species0, referenceN0, nShell, species_name, latVecs)

    type(TMbd), intent(inout) :: this
    type(TMbdInit), intent(inout) :: inp
    type(mpifx_comm), intent(in) :: mympi

    integer, intent(in) :: nAtom
    integer, intent(in) :: species0(:)
    real(dp), intent(in) :: referenceN0(:,:)
    integer, intent(in) :: nShell(:)
    character(mc), intent(in) :: species_name(:)
    real(dp), intent(in), optional :: latVecs(:,:)

    character(len=2) :: spec_name
    integer :: i_flag, i_atom, i_spec, nprow, npcol, npmax, i


    !start filling up MBD object this
    if (present(latVecs)) then
      inp%latvecs = latVecs
    end if
    inp%mbd_stdout = stdOut
    inp%species = species0
    inp%species_names = species_name
    inp%mbd_intra_comm = mympi%id

    !set free atom charges
    allocate(inp%free_charges(nAtom))
    if (inp%mbd_debug .and. tIoProc) then
      write(stdOut,*) 'i_atom    free_charge(i_atom)'
    end if
    do i_atom = 1, nAtom
      i_spec = inp%species(i_atom)
      inp%free_charges(i_atom) = sum(referenceN0(1:nShell(i_spec), i_spec))
      if (inp%mbd_debug.and. tIoProc) then
        write(stdOut,*) i_atom, inp%free_charges(i_atom)
      end if
    end do

    if (inp%mbd_debug .and. tIoProc) then
      do i_atom = 1, nAtom
        write(stdOut,*) i_atom, species0(i_atom), species_name(species0(i_atom))
      end do
    end if

    ! call the actual MBDinit
    call TMbd_init(this, inp)

  end subroutine MBDinit


  !!* Notifies the objects about changed coordinates
  !!* @param this Object instance.
  !!* @param coords Current coordinates of the atoms
  subroutine MBDupdateCoords(this, coords0)
    type(TMbd), intent(inout) :: this
    real(dp), intent(in) :: coords0(:,:)

    !mbd_api coordinate convention is (3,n_atoms)
    call this%updateCoords(coords0)

  end subroutine MBDupdateCoords


  !!* Notifies the object about updated lattice vectors.
  !!* @param latVecs  New lattice vectors
  !!* @param recVecs  New reciprocal vectors
  subroutine MBDupdateLatVecs(this, latVecs, volume)
    type(TMbd), intent(inout) :: this
    real(dp), intent(in) :: latVecs(:,:)
    real(dp), intent(in) :: volume

    call this%updateLatVecs(latVecs, volume)

  end subroutine MBDupdateLatVecs


  !!* Returns the CPA ratios from a mulliken analysis.
  !!* @param this Object instance
  !!* @param
  subroutine MBDcalculateCPA(this, cpatmp)
    type(TMbd), intent(inout)  :: this
    real(dp), intent(inout)    :: cpatmp(:)

    integer :: i_atom

    call this%getScalingRatios(cpatmp)

  end subroutine MBDcalculateCPA


  !!* Returns the MBD energy due to the dispersion.
  !!* @param this Object instance
  !!* @param energy contains the total MBD energy on exit.
  subroutine MBDgetEnergy(this, energy)
    type(TMbd), intent(inout) :: this
    real(dp), intent(inout) :: energy

    call this%getEnergy(energy)

    if (tIoProc) then
      write(stdOut,*) 'MBD energy ', energy, ' ', energy * Hartree__eV
    end if

  end subroutine MBDgetEnergy


  subroutine MBDgetGradients(this, gradients)
    type(TMbd), intent(inout) :: this
    real(dp), intent(out) :: gradients(:,:)

    integer :: i_cart, i_atom

    call this%getGradients(gradients)

    if (tIoProc) then
      write(stdOut,*) '!!!!!!!!!!!!!!!CALCULATING MBD GRADIENTS!!!!!!!!!!!!!!!!'
      write(stdOut,*) gradients(:, :)
    end if

  end subroutine MBDgetGradients


  !!* Calculate finite-difference strain
  !!* @param this Object instance
  !!* @return Cutoff
  subroutine MBDgetStress(this, stress)
    type(TMbd), intent(inout) :: this
    real(dp), intent(inout) :: stress(:,:)

    call this%getStress(stress)

    if (tIoProc) then
      write(stdOut,*) '!!!!!!!!!!!!!!!CALCULATING MBD STRESS!!!!!!!!!!!!!!!!'
      write(stdOut,*) stress(:, :)
    endif

  end subroutine MBDgetStress


  !!* Calculate the ON-Site Mulliken population for each orbital in the system
  !!* using purely real-space overlap and density matrix values.
  !!* Currently Mulliken defined as
  !!* $q_a = \sum_k w_k\sum_{\mu on a} \rho_{\mu\mu}(k)$
  !!* but transformed into real space sums over one triangle of real space
  !!* extended matrices
  !!* @param qq The charge per orbital
  !!* @param rho Density matrix in Packed format
  !!* @param orb Information about the orbitals.
  !!* @param iNeighbor Number of neighbours of each real atom (central cell)
  !!* @param nNeighbor List of neighbours for each atom, starting at 0 for
  !!* itself
  !!* @param img2CentCell indexing array to convert images of atoms
  !!* back into their number in the central cell
  !!* @param iPair indexing array for the Hamiltonian
  !!* @todo add description of algorithm to programer manual / documentation.
  subroutine onsiteMullikenPerAtom(qq, rho, orb, iNeighbor, nNeighbor, img2CentCell, iPair)
    real(dp), intent(inout) :: qq(:)
    real(dp), intent(in) :: rho(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iPair(0:,:)

    integer   :: iOrig
    integer   :: iNeigh
    integer   :: nAtom, iAtom1, iAtom2, iAtom2f
    integer   :: nOrb1, nOrb2, iOrb
    real(dp)  :: sqrTmp(orb%mOrb,orb%mOrb)
    real(dp)  :: mulTmp(orb%mOrb**2)

    nAtom = size(orb%nOrbAtom)

    do iAtom1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAtom1)
      do iNeigh = 0, nNeighbor(iAtom1)
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        if (iAtom1==iAtom2) then
            sqrTmp(:,:) = 0.0_dp
            mulTmp(:) = 0.0_dp
            nOrb2 = orb%nOrbAtom(iAtom2f)
            iOrig = iPair(iNeigh,iAtom1) + 1
            mulTmp(1:nOrb1*nOrb2) = rho(iOrig:iOrig+nOrb1*nOrb2-1)
            sqrTmp(1:nOrb2,1:nOrb1) = &
                & reshape(mulTmp(1:nOrb1*nOrb2), (/nOrb2,nOrb1/))
            do iOrb=1, nOrb1
              qq(iAtom1) = qq(iAtom1) &
                &+ sqrTmp(iOrb,iOrb)
            enddo
        endif
      end do
    end do

  end subroutine onsiteMullikenPerAtom


end module manybodydisp
