!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"
#:include "error.fypp"

!> excitations energies according to the particle-particle Random Phase Approximation
!> (doi:10.1063/1.4977928)
module dftbp_pprpa
  use dftbp_assert
  use dftbp_linrespcommon
  use dftbp_commontypes
  use dftbp_slakocont
  use dftbp_shortgamma
  use dftbp_accuracy
  use dftbp_constants, only : Hartree__eV
  use dftbp_scc, only : TScc
  use dftbp_blasroutines
  use dftbp_eigensolver
  use dftbp_message
  use dftbp_taggedoutput, only : TTaggedWriter, tagLabels
  use dftbp_sorting
  use dftbp_qm
  use dftbp_transcharges
  use dftbp_densedescr
  use dftbp_fileid
  implicit none
  private

  public :: ppRPAenergies, TppRPAcal

  !> Data type for pp-RPA calculations
  type :: TppRPAcal

    !> number of excitations to be found
    integer :: nExc

    !> symmetry of states being calculated
    character :: sym

    !> Coulomb part of atom resolved Hubbard U
    real(dp), allocatable :: hhubbard(:)

    !> Tamm Dancoff Approximation?
    logical :: tTDA

    !> virtual orbital constraint?
    logical :: tConstVir

    !> number of virtual orbitals
    integer :: nvirtual

  end type TppRPAcal

  !> Name of output file
  character(*), parameter :: excitationsOut = "ppRPA_ener.DAT"


contains


  !> This subroutine analytically calculates excitations energies
  !> based on Time Dependent DFRT
  subroutine ppRPAenergies(RPA, denseDesc, grndEigVecs, grndEigVal, sccCalc, SSqr, species0, rnel,&
      & iNeighbour, img2CentCell, orb, tWriteTagged, autotestTag, taggedWriter, err)

    !> Container for RPA calculation data
    type(TppRPAcal), allocatable, intent(in) :: RPA

    !> index vector for S and H matrices
    type(TDenseDescr), intent(in) :: denseDesc

    !> ground state MO-coefficients
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> ground state MO-energies
    real(dp), intent(in) :: grndEigVal(:,:)

    !> Self-consistent charge module settings
    type(TScc), intent(in) :: sccCalc

    !> square overlap matrix between basis functions, both triangles required
    real(dp), intent(in) :: SSqr(:,:)

    !> chemical species of each atom
    integer, intent(in) :: species0(:)

    !> real number of electrons in system
    real(dp), intent(in) :: rnel

    !> Atomic neighbour lists
    integer, intent(in) :: iNeighbour(0:,:)

    !> Mapping of atom number to central cell atom number
    integer, intent(in) :: img2CentCell(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> print tag information
    logical, intent(in) :: tWriteTagged

    !> File name for regression data
    character(*), intent(in) :: autotestTag

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Error code return, 0 if no problems
    integer, intent(out), optional :: err

    logical :: tSpin

    real(dp), allocatable :: stimc(:,:,:)
    real(dp), allocatable :: gamma_eri(:,:)

    character, allocatable :: symmetries(:)

    !> file handle for excitation energies
    integer:: fdExc

    integer, allocatable :: iAtomStart(:)
    integer :: nocc, nvir, nxoo, nxvv, dim_rpa
    integer :: norb, natom
    integer :: iSpin, isym
    integer :: nSpin
    integer :: fdTagged
    character :: sym

    integer :: homo = 0
    real(dp) :: eval_0 = 0.0_dp

    !> ppRPA eigenvalues (two-electron addition/removal energies)
    real(dp), allocatable :: pp_eval(:)

    !> ppRPA eigenvector
    real(dp), allocatable :: vr(:,:)

  @:ASSERT(allocated(RPA))

    if (present(err)) then
      err = 0
    end if

    nSpin = size(grndEigVal, dim=2)
    @:ASSERT(nSpin > 0 .and. nSpin <=2)

    norb = orb%nOrb
    natom = size(species0)
    tSpin = (nSpin == 2)

    allocate(iAtomStart(size(denseDesc%iAtomStart)))
    iAtomStart = denseDesc%iAtomStart


    ! Select symmetries to process
    if (.not. tSpin) then
      select case (RPA%sym)
      case ("B")
        allocate(symmetries(2))
        symmetries(:) = [ "S", "T" ]
      case ("S")
        allocate(symmetries(1))
        symmetries(:) = [ "S" ]
      case ("T")
        ! like B. Triplet calculation requires first singlet exc. energy
        allocate(symmetries(2))
        symmetries(:) = [ "S", "T" ]
      end select
    else
      @:ERROR_HANDLING(err, -1, "spin-unrestricted calculations currently not possible with pp-RPA")
    end if

    ! Allocation for general arrays
    allocate(gamma_eri(natom, natom))
    allocate(stimc(norb, norb, nSpin))

    if (tWriteTagged) then
      open(newUnit=fdTagged, file=autotestTag, position="append")
    end if

    ! Overlap times wave function coefficients - most routines in DFTB+ use lower triangle (would
    ! remove the need to symmetrize the overlap and ground state density matrix in the main code if
    ! this could be used everywhere in these routines)
    do iSpin = 1, nSpin
      call symm(stimc(:,:,iSpin), "L", SSqr, grndEigVecs(:,:,iSpin))
    end do

    nocc = nint(rnel)
    if (abs(rnel - real(nocc,dp)) > epsilon(0.0_dp)) then
      @:ERROR_HANDLING(err, -1, "Fractionally charged systems not possible with pp-RPA")
    end if
    if (mod(abs(nocc),2) == 1) then
      @:ERROR_HANDLING(err, -1, "Odd numbers of electrons not possible with pp-RPA")
    end if
    nocc = nocc / 2

    if ((.not. RPA%tConstVir) .or. ( (RPA%tConstVir) .and. (RPA%nVirtual > nOrb-nocc) )) then
      nvir = nOrb - nocc
    else
      nvir = RPA%nVirtual
    end if

    ! elements in a triangle plus the diagonal of the occ-occ and virt-virt blocks
    nxoo = (nocc * (nocc + 1)) / 2
    nxvv = (nvir * (nvir + 1)) / 2

    ! ground state coulombic interactions
    call sccCalc%getAtomicGammaMatU(gamma_eri, RPA%hHubbard, species0, iNeighbour, img2CentCell)

    ! excitation energies  output file
    open(newUnit=fdExc, file=excitationsOut, position="rewind", status="replace")
    write(fdExc,*)
    write(fdExc,'(5x,a,6x,a,14x,a,6x,a,12x,a)') 'w [eV]', 'Transitions', 'Weight', 'KS [eV]',&
        & 'Symm.'
    write(fdExc,*)
    write(fdExc,'(1x,80("="))')
    write(fdExc,*)

    do isym = 1, size(symmetries)

      sym = symmetries(isym)

      if (sym == "S") then
        if (.not. RPA%tTDa) then
          dim_rpa = nxvv + nxoo
        else
          dim_rpa = nxvv
        end if
      else
        if (.not. RPA%tTDa) then
          dim_rpa = nxvv + nxoo - nvir - nocc
        else
          dim_rpa = nxvv - nvir
        end if
      end if

      allocate(pp_eval(dim_rpa))
      allocate(vr(dim_rpa, dim_rpa))

      call buildAndDiagppRPAmatrix(RPA%tTDa, sym, grndEigVal(:,1), nocc, nvir, nxvv, nxoo,&
          & iAtomStart, gamma_eri, stimc, grndEigVecs, pp_eval, vr, err)

      call writeppRPAExcitations(RPA%tTDa, sym, grndEigVal(:,1), RPA%nExc, pp_eval, vr, nocc, nvir,&
          & nxvv, nxoo, fdExc, tWriteTagged, fdTagged, taggedWriter, eval_0, homo)

      deallocate(pp_eval)
      deallocate(vr)
    end do

    close(fdExc)

    if (tWriteTagged) then
      close(fdTagged)
    end if

  end subroutine ppRPAenergies


  !> Builds and diagonalizes the pp-RPA matrix
  subroutine buildAndDiagppRPAmatrix(tTDA, sym, eigVal, nocc, nvir, nxvv, nxoo, ind, gamma_eri,&
      & stimc, cc, pp_eval, vr, err)

    !> Tamm-Dancoff approximation?
    logical, intent(in) :: tTDA

    !> symmetry to calculate transitions
    character, intent(in) :: sym

    !> ground state MO-energies
    real(dp), intent(in) :: eigVal(:)

    !> number of occupied orbitals
    integer, intent(in)  :: nocc

    !> number of virtual orbitals
    integer, intent(in)  :: nvir

    !> number of virtual-virtual transitions
    integer, intent(in)  :: nxvv

    !> number of occupied-occupied transitions
    integer, intent(in)  :: nxoo

    !> indexing array for square matrices
    integer, intent(in)  :: ind(:)

    !> coulomb interaction
    real(dp), intent(in) :: gamma_eri(:,:)

    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: stimc(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: cc(:,:,:)

    !> pp-RPA eigenvalues
    real(dp), intent(out):: pp_eval(:)

    !> pp-RPA eigenvectors
    real(dp), intent(out) :: vr(:,:)

    !> Error code return, 0 if no problems
    integer, intent(out), optional :: err

    real(dp):: A_s(nxvv,nxvv), A_t(nxvv-nvir,nxvv-nvir)
    real(dp):: B_s(nxvv,nxoo), B_t(nxvv-nvir,nxoo-nocc)
    real(dp):: C_s(nxoo,nxoo), C_t(nxoo-nocc,nxoo-nocc)
    integer :: a, b, c, d, i, j, k, l, ab, cd, ij, kl
    integer :: ab_r, cd_r, ij_r, kl_r
    integer :: at1, at2, natom
    real(dp):: factor1, factor2
    logical :: updwn
    real(dp):: q_1(size(gamma_eri, dim=1))
    real(dp):: q_2(size(gamma_eri, dim=1))
    real(dp):: q_3(size(gamma_eri, dim=1))
    real(dp):: q_4(size(gamma_eri, dim=1))
    integer :: info, nRPA, nxvv_r, nxoo_r
    real(dp), allocatable :: work(:)
    real(dp), allocatable :: wi(:), vl(:,:)
    real(dp), allocatable :: PP(:,:)

    real(dp), parameter :: sqrtFact = 1.0_dp/sqrt(2.0_dp)

    if (present(err)) then
      err = 0
    end if

    natom  = size(gamma_eri, dim=1)
    nRPA   = size(pp_eval)
    nxoo_r = nxoo - nocc
    nxvv_r = nxvv - nvir

    allocate(work(5*nRPA))
    allocate(wi(nRPA))
    allocate(PP(nRPA, nRPA))

    ! spin-up channel only for closed-shell systems
    updwn = .true.

    if (sym == "S") then !------ singlets -------

      ! build Matrix A
      A_s(:,:) = 0.0_dp

      do ab = 1, nxvv

        call indxvv(nocc, ab, a, b)

        factor1 = 1.0_dp
        if (a == b) then
          factor1 = sqrtFact
        end if

        A_s(ab,ab) = A_s(ab,ab) + eigVal(a) + eigVal(b)

        do cd = ab, nxvv

          call indxvv(nocc, cd, c, d)

          factor2 = 1.0_dp
          if (c == d) then
            factor2 = sqrtFact
          end if

          q_1(:) = transq(a, c, ind, updwn, stimc, cc)
          q_2(:) = transq(b, d, ind, updwn, stimc, cc)
          q_3(:) = transq(a, d, ind, updwn, stimc, cc)
          q_4(:) = transq(b, c, ind, updwn, stimc, cc)


          !A_s(ab,cd) = A_s(ab,cd) + factor1 * factor2 * (&
          ! & dot_product(q_1, matmul(gamma_eri, q_2)) + dot_product(q_3, matmul(gamma_eri, q_4)))

          do at1 = 1, natom
            do at2 = 1, natom
              !optimize this: do at2 = at1, natom
              A_s(ab,cd) = A_s(ab,cd) + gamma_eri(at1,at2)*(q_1(at1)*q_2(at2) + q_3(at1)*q_4(at2))&
                  & * factor1 * factor2
            end do
          end do

          if (ab /= cd) then
            A_s(cd,ab) = A_s(cd,ab) + A_s(ab,cd)
          end if

        end do

      end do

      if (.not. tTDA) then
        ! build Matrix B
        B_s(:,:) = 0.0_dp

        do ab = 1, nxvv
          call indxvv(nocc, ab, a, b)
          factor1 = 1.0_dp
          if (a == b) then
            factor1 = sqrtFact
          end if
          do kl = 1, nxoo
            call indxoo(kl, k, l)
            factor2 = 1.0_dp
            if (k == l) then
              factor2 = sqrtFact
            end if

            q_1(:) = transq(k, a, ind, updwn, stimc, cc)
            q_2(:) = transq(l, b, ind, updwn, stimc, cc)
            q_3(:) = transq(l, a, ind, updwn, stimc, cc)
            q_4(:) = transq(k, b, ind, updwn, stimc, cc)


            do at1 = 1, natom
              do at2 = 1, natom
                B_s(ab,kl) = B_s(ab,kl) +gamma_eri(at1,at2)*(q_1(at1)*q_2(at2) + q_3(at1)*q_4(at2))&
                    & * factor1 * factor2
              end do
            end do

          end do
        end do

        ! build Matrix C
        C_s(:,:) = 0.0_dp

        do ij = 1, nxoo
          call indxoo(ij, i, j)
          factor1 = 1.0_dp
          if (i == j) then
            factor1 = sqrtFact
          end if

          C_s(ij,ij) = C_s(ij,ij) - eigVal(i) - eigVal(j)

          do kl = ij, nxoo
            call indxoo(kl, k, l)
            factor2 = 1.0_dp
            if (k == l) then
              factor2 = sqrtFact
            end if

            q_1(:) = transq(i, k, ind, updwn, stimc, cc)
            q_2(:) = transq(j, l, ind, updwn, stimc, cc)
            q_3(:) = transq(i, l, ind, updwn, stimc, cc)
            q_4(:) = transq(j, k, ind, updwn, stimc, cc)

            do at1 = 1, natom
              do at2 = 1, natom
                C_s(ij,kl) = C_s(ij,kl) + gamma_eri(at1,at2)&
                &                *(q_1(at1)*q_2(at2) + q_3(at1)*q_4(at2))&
                &                *factor1*factor2
              end do
            end do

            if (ij /= kl) then
              C_s(kl,ij) = C_s(kl,ij) + C_s(ij,kl)
            end if
          end do

        end do

        !build ppRPA matrix
        PP(1:nxvv,1:nxvv)   =  A_s(:,:)
        PP(1:nxvv,nxvv+1:)  =  B_s(:,:)
        PP(nxvv+1:,1:nxvv)  = -transpose(B_s)
        PP(nxvv+1:,nxvv+1:) = -C_s(:,:)
      else
        PP(:,:) = A_s(:,:)
      end if

    else !-------- triplets ----------

      ! build Matrix A
      A_t(:,:) = 0.0_dp
      ab_r = 0
      do ab = 1, nxvv
        call indxvv(nocc, ab, a, b)
        if (a == b) then
          cycle
        end if
        ab_r = ab_r + 1

        A_t(ab_r,ab_r) = A_t(ab_r,ab_r) + eigVal(a) + eigVal(b)

        cd_r = ab_r - 1
        do cd = ab, nxvv
          call indxvv(nocc, cd, c, d)
          if (c == d) then
            cycle
          end if
          cd_r = cd_r + 1

          q_1(:) = transq(a, c, ind, updwn, stimc, cc)
          q_2(:) = transq(b, d, ind, updwn, stimc, cc)
          q_3(:) = transq(a, d, ind, updwn, stimc, cc)
          q_4(:) = transq(b, c, ind, updwn, stimc, cc)

          do at1 = 1, natom
            do at2 = 1, natom
              !optimize this: do at2 = at1, natom
              A_t(ab_r,cd_r) = A_t(ab_r,cd_r) + gamma_eri(at1,at2)&
                 &              *(q_1(at1)*q_2(at2) - q_3(at1)*q_4(at2))
            end do
          end do

          if (ab_r /= cd_r) then
            A_t(cd_r,ab_r) = A_t(cd_r,ab_r) + A_t(ab_r,cd_r)
          end if
        end do

      end do

      if (.not. tTDA) then
        ! build Matrix B
        B_t(:,:) = 0.0_dp
        ab_r = 0
        do ab = 1, nxvv
          call indxvv(nocc, ab, a, b)
          if (a == b) then
            cycle
          end if
          ab_r = ab_r + 1

          kl_r = 0
          do kl = 1, nxoo
            call indxoo(kl, k, l)
            if (k == l) then
              cycle
            end if
            kl_r = kl_r + 1

            q_1(:) = transq(k, a, ind, updwn, stimc, cc)
            q_2(:) = transq(l, b, ind, updwn, stimc, cc)
            q_3(:) = transq(l, a, ind, updwn, stimc, cc)
            q_4(:) = transq(k, b, ind, updwn, stimc, cc)

            do at1 = 1, natom
              do at2 = 1, natom
                B_t(ab_r,kl_r) = B_t(ab_r,kl_r) + gamma_eri(at1,at2)&
                   &             *(q_1(at1)*q_2(at2) - q_3(at1)*q_4(at2))
              end do
            end do

          end do
        end do

        ! build Matrix C
        C_t(:,:) = 0.0_dp
        ij_r = 0
        do ij = 1, nxoo
          call indxoo(ij, i, j)
          if (i == j) then
            cycle
          end if
          ij_r = ij_r + 1

          C_t(ij_r,ij_r) = C_t(ij_r,ij_r) - eigVal(i) - eigVal(j)

          kl_r = ij_r - 1
          do kl = ij, nxoo

            call indxoo(kl, k, l)
            if (k == l) then
              cycle
            end if
            kl_r = kl_r + 1

            q_1(:) = transq(i, k, ind, updwn, stimc, cc)
            q_2(:) = transq(j, l, ind, updwn, stimc, cc)
            q_3(:) = transq(i, l, ind, updwn, stimc, cc)
            q_4(:) = transq(j, k, ind, updwn, stimc, cc)

            do at1 = 1, natom
              do at2 = 1, natom
                C_t(ij_r,kl_r) = C_t(ij_r,kl_r) + gamma_eri(at1,at2)&
                    & * (q_1(at1)*q_2(at2) - q_3(at1)*q_4(at2))
              end do
            end do

            if (ij_r /= kl_r) then
              C_t(kl_r,ij_r) = C_t(kl_r,ij_r) + C_t(ij_r,kl_r)
            end if

          end do

        end do

        !build ppRPA matrix
        PP(1:nxvv_r,1:nxvv_r)   =  A_t
        PP(1:nxvv_r,nxvv_r+1:)  =  B_t
        PP(nxvv_r+1:,1:nxvv_r)  = -transpose(B_t)
        PP(nxvv_r+1:,nxvv_r+1:) = -C_t
      else
        PP(:,:) = A_t
      end if
    end if

    !diagonalize ppRPA matrix
    call geev(PP, pp_eval, wi, vl, vr, info)

    if (info /= 0) then
      @:FORMATTED_ERROR_HANDLING(err, info, "(A,I0)", " Error with dgeev, info = ", info)
    end if

  end subroutine buildAndDiagppRPAmatrix


  !> write pp-RPA excitation energies in output file
  subroutine writeppRPAExcitations(tTDA, sym, eigVal, nexc, pp_eval, vr, nocc, nvir, nxvv, nxoo,&
      & fdExc, tWriteTagged, fdTagged, taggedWriter, eval_0, homo)

    !> Tamm-Dancoff approximation?
    logical, intent(in) :: tTDA

    !> symmetry to calculate transitions
    character, intent(in) :: sym

    !> ground state MO-energies
    real(dp), intent(in) :: eigVal(:)

    !> number of excitation energies to show
    integer, intent(in) :: nexc

    !> pp-RPA eigenvalues
    real(dp), intent(inout) :: pp_eval(:)

    !> pp-RPA eigenvectors
    real(dp), intent(in) :: vr(:,:)

    !> number of occupied orbitals
    integer, intent(in)  :: nocc

    !> number of virtual orbitals
    integer, intent(in)  :: nvir

    !> number of virtual-virtual transitions
    integer, intent(in)  :: nxvv

    !> number of occupied-occupied transitions
    integer, intent(in)  :: nxoo

    !> file unit for excitation energies
    integer, intent(in) :: fdExc

    !> print tag information
    logical, intent(in) :: tWriteTagged

    !> file unit for tagged output (> -1 for write out)
    integer, intent(in) :: fdTagged

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> first eigenvalue (considers singlet ground state of N-system)
    real(dp), intent(inout) :: eval_0

    !> HOMO of N-system
    integer, intent(inout) :: homo

    integer :: i, ii, aa, bb, diff_a, diff_b
    integer :: nRPA, nxoo_r, nxvv_r
    integer, allocatable :: e_ind(:)
    real(dp), allocatable :: norm(:)
    real(dp), allocatable :: wvr(:)
    integer, allocatable :: wvin(:)
    real(dp):: weight
    real(dp):: ks_ener_a, ks_ener_b

    nRPA = size(pp_eval)
    allocate(e_ind(nRPA))
    allocate(wvr(nRPA))
    allocate(wvin(nRPA))


    if (sym == "S") then
      nxoo_r = nxoo
      nxvv_r = nxvv
    else
      nxoo_r = nxoo - nocc
      nxvv_r = nxvv - nvir
    end if

    !neglect negative-norm eigenvectors
    if (.not. tTDA) then
      allocate(norm(nRPA))
      norm(:) = sum(vr(:nxvv_r,:)**2, dim=1)
      norm(:) = norm - sum(vr(nxvv_r+1:nxvv_r+nxoo_r,:)**2, dim=1)
    end if

    call index_heap_sort(e_ind, pp_eval)

    pp_eval = pp_eval(e_ind)

    ii = 0
    RPA: do i = 1, nRPA

      if (.not.tTDA) then
        if (norm(e_ind(i)) < 0.0_dp) then
          cycle RPA
        end if
      end if

      wvr(:) = vr(:,e_ind(i))**2
      wvr(:) = wvr / sqrt(sum(wvr**2))

      call index_heap_sort(wvin, wvr)
      wvin = wvin(size(wvin):1:-1)
      wvr = wvr(wvin)

      call indxvv(nocc, wvin(1), aa, bb)
      if (sym == "T") then
        aa = aa + 1
      end if

      weight = wvr(1)

      ii = ii + 1
      if ((ii == 1) .and. (sym == "S")) then
        eval_0 = pp_eval(i)
        homo = aa
        cycle
      endif

      diff_a = aa - homo - 1
      diff_b = bb - homo - 1
      ks_ener_a = Hartree__eV * (eigVal(aa) - eigVal(homo))
      ks_ener_b = Hartree__eV * (eigVal(bb) - eigVal(homo))

      !> double excitations
      if (diff_b /= -1) then
        write(fdExc,'(1x,f10.3,6x,a,i2,a,i2,5x,f6.3,6x,f7.3,a,f7.3,7x,a)')&
            & Hartree__eV * (pp_eval(i) - eval_0),&
            & 'HOMO -> LUMO +', diff_b,',', diff_a, weight, ks_ener_b, ',', ks_ener_a, sym
      else
        write(fdExc,'(1x,f10.3,6x,a,i2,8x,f6.3,6x,f7.3,15x,a)')&
            & Hartree__eV * (pp_eval(i) - eval_0),&
            & 'HOMO -> LUMO +', diff_a, weight, ks_ener_a, sym
      endif

      if ( ((ii > nExc) .and. (sym == "S")) .or. ((ii == nExc) .and. (sym == "T")) ) then
         exit
      end if

    end do RPA

    if (tWriteTagged) then
      call taggedWriter%write(fdTagged, tagLabels%egyppRPA, pp_eval(:nRPA))
    end if

  end subroutine writeppRPAExcitations


end module dftbp_pprpa
