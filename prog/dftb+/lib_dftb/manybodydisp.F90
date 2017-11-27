#:include "common.fypp"

!!* Dispersion with Many Body Dispersion method by Tkatchenko et al.
module manybodydisp
  use globalenv, only : stdOut, tIoProc
  use mpifx, only : mpifx_comm
  use accuracy
  use constants, only: pi, Bohr__AA, AA__Bohr, eV__Hartree, Hartree__eV
  use mbd, only: TMbdInit, TMbd
  use commontypes, only : TOrbitals

  implicit none
  private

  public:: TMbdInit, TMbd

#:if WITH_MBD
  public:: MBDinit, MBDupdateCoords, MBDupdateLatVecs
  public:: MBDgetEnergy, MBDaddGradients, MBDcalculateCPA
  ! Stress convention. Stress provided by MBD_API is in energy units per volume Positive stress
  ! expands the unit cell, negative stress contracts it DOUBLE CHECK WITH YOUR CODE IF IT FOLLOWS
  ! THE SAME CONVENTION (This is the convention used by DFTB+)
  public:: MBDgetStress
#:endif

  public :: onsiteMullikenPerAtom

contains

  !takes TMbd instance, TMbdInit instance, nAtom, species0,
  !referenceN0, nShell, species_name, latVecs, recVecs
  subroutine MBDinit(this, inp, mympi, nAtom, species0, referenceN0, &
          & nShell, species_name, latVecs, recVecs)
    use message, only: error
    use mbd_api, only: initialize => MBDInit

    type(TMbd), intent(inout) :: this
    type(TMbdInit), intent(inout) :: inp
    type(mpifx_comm), intent(in) :: mympi

    integer, intent(in) :: nAtom
    integer, intent(in) :: species0(:)
    real(dp), intent(in) :: referenceN0(:,:)
    integer, intent(in) :: nShell(:)
    character(mc), intent(in) :: species_name(:)
    real(dp), intent(in), optional :: latVecs(:,:), recVecs(:,:)

    character(len=2) :: spec_name
    integer :: i_flag, i_atom, i_spec, nprow, npcol, npmax, i

    !start filling up MBD object this
    this%n_atoms = nAtom
    this%n_species = len(species_name)
    this%is_periodic = present(latVecs)

    this%mbd_stdOut = stdOut

    allocate(this%species(this%n_atoms))
    allocate(this%species_name(this%n_species))
    allocate(this%free_charge(this%n_atoms))
    this%species = species0
    this%species_name = species_name

    !TODO this needs to go into API
    !satithisy dependencies of mbdvdw_interface
    this%mbd_myid = mympi%rank
    this%mbd_n_tasks = mympi%size
    this%mbd_intra_comm = mympi%id

    if (this%is_periodic) then
        this%latvecs = latVecs
        this%recvecs = recVecs
    endif

    !set free atom charges
    if (inp%mbd_debug.and. tIoProc)  write(stdOut,*) 'i_atom    free_charge(i_atom)'
    do i_atom=1, this%n_atoms
      i_spec = this%species(i_atom)
      this%free_charge(i_atom) = sum(referenceN0(1:nShell(i_spec), i_spec))
      if (inp%mbd_debug.and. tIoProc) write(stdOut,*) i_atom, this%free_charge(i_atom)
    enddo

    if (inp%mbd_debug.and. tIoProc) then
      do i_atom = 1, this%n_atoms
        write(stdOut,*) i_atom, species0(i_atom) ,species_name(species0(i_atom))
      end do
    endif

    !call the actual MBDinit
    call initialize(this, inp)

    this%tInit = .true.

  end subroutine MBDinit


  !!* Notifies the objects about changed coordinates
  !!* @param this Object instance.
  !!* @param coords Current coordinates of the atoms
  subroutine MBDupdateCoords(this,coords)
    use mbd_api, only: updatecoords => MBDupdateCoords
    type(TMbd), intent(inout) :: this
    real(dp), intent(inout) :: coords(:,:)

    !mbd_api coordinate convention is (3,n_atoms)
    call updatecoords(this, coords)

    this%coordsUpdated = .true.

  end subroutine MBDupdateCoords

  !!* Notifies the object about updated lattice vectors.
  !!* @param latVecs  New lattice vectors
  !!* @param recVecs  New reciprocal vectors
  subroutine MBDupdateLatVecs(this, latVecs, recVecs)
    use mbd_api, only: updatelatvecs => MBDupdateLatVecs
    type(TMbd), intent(inout) :: this
    real(dp), intent(in) :: latVecs(:,:), recVecs(:,:)

    call updatelatvecs(this, latVecs, recVecs)

  end subroutine MBDupdateLatVecs

  !!* Returns the CPA ratios from a mulliken analysis.
  !!* @param this Object instance
  !!* @param
  subroutine MBDcalculateCPA(this, cpatmp)
    use mbd_api, only: MBDcalculate_scaling_ratios
    type(TMbd), intent(inout)  :: this
    real(dp), intent(inout)    :: cpatmp(:)

    integer :: i_atom

    call MBDcalculate_scaling_ratios(this, cpatmp)

    if (this%mbd_debug .and. (tIoProc)) then
        write(stdOut,*) ' '
        write(stdOut,*) 'iAtom  CPA ratio'
        do i_atom=1, this%n_atoms
          write(stdOut,*) i_atom, cpatmp(i_atom)
        enddo
    endif
  end subroutine MBDcalculateCPA


  !!* Returns the MBD energy due to the dispersion.
  !!* @param this Object instance
  !!* @param energy contains the total MBD energy on exit.
  subroutine MBDgetEnergy(this, energy)
    use mbd_api, only: get_energy => MBDgetEnergy
    type(TMbd), intent(inout) :: this
    real(dp), intent(inout) :: energy

    if (this%mbd_debug .and. tIoProc) &
        write(stdOut,*) '!!!!!!!!!!!!!!!CALCULATING MBD ENERGY!!!!!!!!!!!!!!!!'

    call get_energy(this, energy)

    if (tIoProc) write(stdOut,*) 'MBD energy ', energy, ' ', energy*27.2114

  end subroutine MBDgetEnergy

  subroutine MBDaddGradients(this, gradients)
    use mbd_api, only: add_gradients => MBDaddGradients
    type(TMbd), intent(inout) :: this
    real(dp), intent(inout) :: gradients(:,:)

    integer :: i_cart, i_atom

    call add_gradients(this, gradients)

    if (this%mbd_debug.and. tIoProc) &
        write(stdOut,*) '!!!!!!!!!!!!!!!CALCULATING MBD GRADIENTS!!!!!!!!!!!!!!!!'
    if (this%mbd_debug.and. tIoProc) then
      write(stdOut,*) gradients(:, :)
    endif

  end subroutine MBDaddGradients

  !!* Calculate finite-difference strain
  !!* @param this Object instance
  !!* @return Cutoff
  subroutine MBDgetStress(this, stress, latVec, orig_vol)
    use mbd_api, only: get_stress => MBDgetStress
    type(TMbd), intent(inout) :: this
    real(dp), intent(inout) :: stress(:,:)
    real(dp), intent(in) :: latVec(:,:)
    real(dp), intent(in) :: orig_vol

    if (this%mbd_debug.and. tIoProc) &
        write(stdOut,*) '!!!!!!!!!!!!!!!CALCULATING MBD STRESS!!!!!!!!!!!!!!!!'

    call get_stress(this, stress, latVec, orig_vol)

    if (this%mbd_debug.and. tIoProc) then
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
