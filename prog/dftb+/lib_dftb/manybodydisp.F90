!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

!> Dispersion with Many Body Dispersion method by Tkatchenko et al.
module manybodydisp
  use assert
  use globalenv, only : stdOut, tIoProc
  use mpifx, only : mpifx_comm
  use accuracy
  use constants, only: pi, Bohr__AA, AA__Bohr, eV__Hartree, Hartree__eV
  use mbd_module, TMbdInit => mbd_input, TMbd => mbd_calculation
  use commontypes, only : TOrbitals
  use message, only: error
  implicit none
  private

  public:: TMbdInit, TMbd
  public:: MBDinit, MBDupdateCoords, MBDupdateLatVecs
  public:: MBDgetEnergy, MBDgetGradients, MBDcalculateCPA
  public:: MBDgetStress

contains

  !> Initialises the MBD use
  subroutine MBDinit(this, inp, mympi, nAtom, species0, referenceN0, nShell, species_name, latVecs)

    !> ADT holding the structure
    type(TMbd), intent(inout) :: this

    !> initial input to the library
    type(TMbdInit), intent(inout) :: inp

    !> MPI common details
    type(mpifx_comm), intent(in) :: mympi

    !> Number of (central cell) atoms
    integer, intent(in) :: nAtom

    !> central cell atomic species
    ! TODO why is this needed?
    integer, intent(in) :: species0(:)

    !> Reference free atom charges
    real(dp), intent(in) :: referenceN0(:,:)

    !> Number of atomic shells
    integer, intent(in) :: nShell(:)

    !> labels of the atoms
    character(mc), intent(in) :: species_name(:)

    !> Lattice vectors if periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    ! label for species names
    character(len=2) :: spec_name

    integer :: i_flag, i_atom, i_spec, nprow, npcol, npmax, i

    ! start filling up MBD object
    ! TODO why should the input hold the lattice vectors? These should come into
    ! TMbd via MBDupdateLatVecs
    ! if (present(latVecs)) then
    !   inp%lattice_vectors = latVecs
    ! end if
    inp%free_values = mbd_get_free_vdw_params(species_name, 'ts')

    ! set free atom charges
    ! TODO this should be stored within DFTB+, not the MBD library, which cares
    ! only about the ratios
    ! allocate(inp%free_charges(nAtom))
    ! if (inp%mbd_debug .and. tIoProc) then
    !   write(stdOut,*) 'i_atom    free_charge(i_atom)'
    ! end if
    ! do i_atom = 1, nAtom
    !   i_spec = inp%species(i_atom)
    !   inp%free_charges(i_atom) = sum(referenceN0(1:nShell(i_spec), i_spec))
    !   if (inp%mbd_debug.and. tIoProc) then
    !     write(stdOut,*) i_atom, inp%free_charges(i_atom)
    !   end if
    ! end do

    ! TODO
    ! I removed all output from the mbd library for the moment, so there is not
    ! place for mbd_debug. All messages are stored in mbd_calc%info, and I'll
    ! created some mechanism to retrieve them
    ! if (inp%mbd_debug .and. tIoProc) then
    !   do i_atom = 1, nAtom
    !     write(stdOut,*) i_atom, species0(i_atom), species_name(species0(i_atom))
    !   end do
    ! end if

    ! call the actual MBDinit
    call this%init(inp)

  end subroutine MBDinit


  !> Notifies the objects about changed coordinates
  subroutine MBDupdateCoords(this, coords0)

    !> Object instance.
    type(TMbd), intent(inout) :: this

    !> Current coordinates of the atoms
    real(dp), intent(in) :: coords0(:,:)

    ! mbd_api coordinate convention is (3,n_atoms)
    call this%update_coords(coords0)

  end subroutine MBDupdateCoords


  !> Notifies the object about updated lattice vectors.
  subroutine MBDupdateLatVecs(this, latVecs, volume)

    !> New lattice vectors
    type(TMbd), intent(inout) :: this

    !> New reciprocal vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> cell volume
    real(dp), intent(in) :: volume

    call this%update_lattice_vectors(latVecs)

  end subroutine MBDupdateLatVecs


  !> Returns the CPA ratios from a mulliken analysis.
  subroutine MBDcalculateCPA(this, cpatmp)

    !> Object instance
    type(TMbd), intent(inout)  :: this

    !> charge population ratio
    real(dp), intent(inout)    :: cpatmp(:)

    integer :: i_atom

    !TODO make this compatible with Hirshfeld and CPA-based methods
    !currently only CPA

    ! do i_atom=1, this%n_atoms
    !   this%scaling_ratio(i_atom) = 1.0_dp + (cpatmp(i_atom)-this%free_charge(i_atom))/ &
    !       real(name2nucleus(this%species_name(this%species(i_atom))))
    !   !passing 'volumes' into MBD module
    !   vefftsvdw(i_atom) = real(name2nucleus(this%species_name(this%species(i_atom)))) + &
    !       (cpatmp(i_atom) - this%free_charge(i_atom))
    !   !passing CPA ratios into sdc_recode module TS
    !   sedc_ts_veff_div_vfree = this%scaling_ratio
    ! enddo
    ! !! overwrite cpatmp for output
    ! cpatmp(:) = this%scaling_ratio

    ! TODO this needs to be called once cpatmp contains the Hirshfeld ratios
    call this%update_vdw_params_from_ratios(cpatmp)

    !RESCALING only if Jans code is included, which it currently is not
    !alpha_0_scaled = cpa_ratio*alpha_0_free
    !C6_scaled = (cpa_ratio**2)*C6_free
    !R_vdw_scaled = (cpa_ratio**(1.d0/3))*R_vdw_free

    !if (this%mbd_debug .and. (ioproc)) then
    !    write(stdout,*) ' '
    !    write(stdout,*) 'iAtom  scaling ratio'
    !    do i_atom=1, this%n_atoms
    !      write(stdout,*) i_atom, this%scaling_ratio(i_atom)
    !    enddo
    !write(stdout,*) ' '
    !write(stdout,*) 'alpha_0_free    C6_free     R_vdw_free'
    !do i_atom=1, n_atoms
    !  write(stdout,*) alpha_0_free(i_atom), C6_free(i_atom), R_vdw_free(i_atom)
    !enddo
    !write(stdout,*) ' '
    !write(stdout,*) 'alpha_0_scaled    C6_scaled     R_vdw_scaled'
    !do i_atom=1, n_atoms
    !  write(stdout,*) alpha_0_scaled(i_atom), C6_scaled(i_atom), R_vdw_scaled(i_atom)
    !enddo
    !endif

  end subroutine MBDcalculateCPA


  !> Returns the MBD energy due to the dispersion.
  subroutine MBDgetEnergy(this, energy)

    !> Object instance
    type(TMbd), intent(inout) :: this

    !> contains the total MBD energy on exit.
    real(dp), intent(inout) :: energy

    call this%get_energy(energy)

    if (tIoProc) then
      write(stdOut,*) 'MBD energy ', energy, ' ', energy * Hartree__eV
    end if

  end subroutine MBDgetEnergy

  !> Energy gradient
  subroutine MBDgetGradients(this, gradients)

    !> Object instance
    type(TMbd), intent(inout) :: this

    !> resulting gradient
    real(dp), intent(out) :: gradients(:,:)

    integer :: i_cart, i_atom

    call this%get_gradients(gradients)

    if (tIoProc) then
      write(stdOut,*) '!!!!!!!!!!!!!!!CALCULATING MBD GRADIENTS!!!!!!!!!!!!!!!!'
      write(stdOut,*) gradients(:, :)
    end if

  end subroutine MBDgetGradients


  !> Calculate finite-difference strain
  subroutine MBDgetStress(this, stress)

    !> Object instance
    type(TMbd), intent(inout) :: this

    !> Stress tensor
    real(dp), intent(inout) :: stress(:,:)

    call this%get_lattice_derivs(stress)

    if (tIoProc) then
      write(stdOut,*) '!!!!!!!!!!!!!!!CALCULATING MBD STRESS!!!!!!!!!!!!!!!!'
      write(stdOut,*) stress(:, :)
    endif

  end subroutine MBDgetStress

end module manybodydisp
