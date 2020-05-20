!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Generalized Born solvation model.
module dftbp_born
  use dftbp_accuracy, only : dp
  use dftbp_blasroutines, only : hemv, gemv
  use dftbp_charges, only : getSummedCharges
  use dftbp_cm5, only : TChargeModel5, TCM5Input, TChargeModel5_init
  use dftbp_commontypes, only : TOrbitals
  use dftbp_constants, only : Hartree__eV
  use dftbp_environment, only : TEnvironment
  use dftbp_periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use dftbp_sasa, only : TSASACont, TSASAInput, TSASACont_init, writeSASAContInfo
  use dftbp_simplealgebra, only : determinant33
  use dftbp_solvation, only : TSolvation
  implicit none
  private

  public :: TGeneralizedBorn, TGBInput, TGeneralizedBorn_init
  public :: writeGeneralizedBornInfo


  !> Global parameters for the solvation
  type :: TGBParameters

    !> Energy shift to the reference system
    real(dp) :: freeEnergyShift = 0.0_dp

    !> Dielectric screening
    real(dp) :: keps = 1.0_dp

    !> Scaling factor for Born radii
    real(dp) :: bornScale = 1.0_dp

    !> Offset parameter for Born radii
    real(dp) :: bornOffset = 0.0_dp

    !> Onufriev--Bashford--Case correction to Born radii
    real(dp) :: obc(3) = [1.00_dp, 0.80_dp, 4.85_dp]

    !> Van-der-Waals radii
    real(dp), allocatable :: vdwRad(:)

  end type TGBParameters

  !> Input parameters to initialize generalized Born model
  type, extends(TGBParameters) :: TGBInput

    !> Dielectric descreening parameter
    real(dp), allocatable :: descreening(:)

    !> Real space cutoff
    real(dp) :: rCutoff = 0.0_dp

    !> Use charge model 5
    type(TCM5Input), allocatable :: cm5Input

    !> Input for solvent accessible surface area model
    type(TSASAInput), allocatable :: sasaInput

    !> Parameter for H-bond correction
    real(dp), allocatable :: hBondPar(:)

  end type TGBInput


  !> Data for the Generalized Born solvation model
  type, extends(TSolvation) :: TGeneralizedBorn
    private

    !> number of atoms
    integer :: nAtom = 0

    !> solvation free energy
    real(dp), allocatable :: energies(:)

    !> lattice vectors if periodic
    real(dp) :: latVecs(3, 3) = 0.0_dp

    !> Volume of the unit cell
    real(dp) :: volume = 0.0_dp

    !> Strain derivatives
    real(dp) :: sigma(3, 3) = 0.0_dp

    !> is this periodic
    logical :: tPeriodic

    !> are the coordinates current?
    logical :: tCoordsUpdated = .false.

    !> are the charges current?
    logical :: tChargesUpdated = .false.

    !> Real space cutoff
    real(dp) :: rCutoff = 0.0_dp

    !> Model parameters
    type(TGBParameters) :: param

    !> Correction to charges with CM5
    type(TChargeModel5), allocatable :: cm5

    !> Born shifts to the hamiltonian
    real(dp), allocatable :: shift(:)

    !> Charges
    real(dp), allocatable :: chargesPerAtom(:)

    !> Born radii
    real(dp), allocatable :: bornRad(:)

    !> Born matrix
    real(dp), allocatable :: bornMat(:, :)

    !> Pair descreening approximation radii
    real(dp), allocatable :: rho(:)

    !> Gradient of the Born radii
    real(dp), allocatable :: dbrdr(:, :, :)

    !> Strain derivative of the Born radii
    real(dp), allocatable :: dbrdL(:, :, :)

    !> Solvent accessible surface area model
    type(TSASACont), allocatable :: sasaCont

    !> Parameter for H-bond correction
    real(dp), allocatable :: hBondStrength(:)

  contains

    !> update internal copy of coordinates
    procedure :: updateCoords

    !> update internal copy of lattice vectors
    procedure :: updateLatVecs

    !> get real space cutoff
    procedure :: getRCutoff

    !> get energy contributions
    procedure :: getEnergies

    !> get force contributions
    procedure :: addGradients

    !> get stress tensor contributions
    procedure :: getStress

    !> Updates with changed charges for the instance
    procedure :: updateCharges

    !> Returns shifts per atom
    procedure :: getShifts
  end type TGeneralizedBorn


contains


  !> Initialize generalized Born model from input data
  subroutine TGeneralizedBorn_init(self, input, nAtom, species0, speciesNames, &
      & latVecs)

    !> Initialised instance at return
    type(TGeneralizedBorn), intent(out) :: self

    !> Specific input parameters for generalized Born
    type(TGBinput), intent(in) :: input

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    integer :: nSpecies
    integer :: iAt1, iSp1

    nSpecies = size(speciesNames)
    self%tPeriodic = present(latVecs)

    if (allocated(input%sasaInput)) then
       allocate(self%sasaCont)
       if (self%tPeriodic) then
         call TSASACont_init(self%sasaCont, input%sasaInput, nAtom, species0, &
             & speciesNames, latVecs)
       else
         call TSASACont_init(self%sasaCont, input%sasaInput, nAtom, species0, &
             & speciesNames)
       end if
    end if

    if (self%tPeriodic) then
      call self%updateLatVecs(LatVecs)
    end if
    self%nAtom = nAtom

    allocate(self%energies(nAtom))
    allocate(self%shift(nAtom))
    allocate(self%chargesPerAtom(nAtom))
    allocate(self%bornRad(nAtom))
    allocate(self%bornMat(nAtom, nAtom))
    allocate(self%rho(nSpecies))
    allocate(self%dbrdr(3, nAtom, nAtom))
    allocate(self%dbrdL(3, 3, nAtom))

    self%param = input%TGBParameters
    self%rho(:) = input%vdwRad(:) * input%descreening(:)

    if (allocated(self%sasaCont) .and. allocated(input%hBondPar)) then
      if (any(input%hBondPar /= 0.0_dp)) then
        allocate(self%hBondStrength(nAtom))
        do iAt1 = 1, nAtom
          iSp1 = species0(iAt1)
          self%hBondStrength(iAt1) = input%hBondPar(iSp1) / self%sasaCont%probeRad(iSp1)**2
        end do
      end if
    end if

    self%rCutoff = input%rCutoff

    if (allocated(input%cm5Input)) then
      allocate(self%cm5)
      if (self%tPeriodic) then
        call TChargeModel5_init(self%cm5, input%cm5Input, nAtom, speciesNames, &
           & .true., latVecs)
      else
        call TChargeModel5_init(self%cm5, input%cm5Input, nAtom, speciesNames, &
           & .true.)
      end if
    end if

    self%tCoordsUpdated = .false.
    self%tChargesUpdated = .false.

  end subroutine TGeneralizedBorn_init


  !> Print the solvation model used
  subroutine writeGeneralizedBornInfo(unit, solvation)

    !> Formatted unit for IO
    integer, intent(in) :: unit

    !> Solvation model
    type(TGeneralizedBorn), intent(in) :: solvation

    write(unit, '(a, ":", t30, es14.6)') "Dielectric constant", &
        & 1.0_dp/(solvation%param%keps + 1.0_dp)
    write(unit, '(a, ":", t30, es14.6, 1x, a, t50, es14.6, 1x, a)') &
        & "Free energy shift", solvation%param%freeEnergyShift, "H", &
        & Hartree__eV * solvation%param%freeEnergyShift, "eV"
    write(unit, '(a, ":", t30, a)') "Born radii integrator", "GBOBC"

    write(unit, '(a, ":", t30)', advance='no') "SASA model"
    if (allocated(solvation%sasaCont)) then
      write(unit, '(a)') "Yes"
      call writeSASAContInfo(unit, solvation%sasaCont)
    else
      write(unit, '(a)') "No"
    end if

    write(unit, '(a, ":", t30)', advance='no') "CM5 correction"
    if (allocated(solvation%cm5)) then
      write(unit, '(a)') "Yes"
    else
      write(unit, '(a)') "No"
    end if

    write(unit, '(a, ":", t30)', advance='no') "Hydrogen bond correction"
    if (allocated(solvation%hBondStrength)) then
      write(unit, '(a)') "Yes"
    else
      write(unit, '(a)') "No"
    end if
  end subroutine writeGeneralizedBornInfo


  !> Update internal stored coordinates
  subroutine updateCoords(self, env, neighList, img2CentCell, coords, species0)

    !> Data structure
    class(TGeneralizedBorn), intent(inout) :: self

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> List of neighbours to atoms
    type(TNeighbourList), intent(in) :: neighList

    !> Image to central cell atom index
    integer, intent(in) :: img2CentCell(:)

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    integer, allocatable :: nNeigh(:)

    if (allocated(self%sasaCont)) then
      call self%sasaCont%updateCoords(env, neighList, img2CentCell, coords, species0)
    end if

    allocate(nNeigh(self%nAtom))
    call getNrOfNeighboursForAll(nNeigh, neighList, self%rCutoff)
    call getBornRadii(self, nNeigh, neighList%iNeighbour, img2CentCell, &
        & neighList%neighDist2, species0, coords)
    call getBornMatrixCluster(self, coords)

    if (allocated(self%cm5)) then
      call self%cm5%updateCoords(neighList, img2CentCell, coords, species0)
    end if

    self%tCoordsUpdated = .true.
    self%tChargesUpdated = .false.

  end subroutine updateCoords


  !> Update internal copy of lattice vectors
  subroutine updateLatVecs(self, latVecs)

    !> Data structure
    class(TGeneralizedBorn), intent(inout) :: self

    !> Lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    @:ASSERT(self%tPeriodic)
    @:ASSERT(all(shape(latvecs) == shape(self%latvecs)))

    if (allocated(self%sasaCont)) then
      call self%sasaCont%updateLatVecs(latVecs)
    end if

    self%volume = abs(determinant33(latVecs))
    self%latVecs(:,:) = latVecs

    if (allocated(self%cm5)) then
      call self%cm5%updateLatVecs(LatVecs)
    end if

    self%tCoordsUpdated = .false.
    self%tChargesUpdated = .false.

  end subroutine updateLatVecs


  !> Get energy contributions
  subroutine getEnergies(self, energies)

    !> data structure
    class(TGeneralizedBorn), intent(inout) :: self

    !> energy contributions for each atom
    real(dp), intent(out) :: energies(:)

    @:ASSERT(self%tCoordsUpdated)
    @:ASSERT(self%tChargesUpdated)
    @:ASSERT(size(energies) == self%nAtom)

    if (allocated(self%sasaCont)) then
      call self%sasaCont%getEnergies(energies)
    else
      energies(:) = 0.0_dp
    end if

    energies(:) = energies + 0.5_dp * (self%shift * self%chargesPerAtom) &
       & + self%param%freeEnergyShift / real(self%nAtom, dp)

  end subroutine getEnergies


  !> Get force contributions
  subroutine addGradients(self, env, neighList, species, coords, img2CentCell, gradients)

    !> Data structure
    class(TGeneralizedBorn), intent(inout) :: self

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Neighbour list
    type(TNeighbourList), intent(in) :: neighList

    !> Specie for each atom
    integer, intent(in) :: species(:)

    !> Coordinate of each atom
    real(dp), intent(in) :: coords(:,:)

    !> Mapping of atoms to cetnral cell
    integer, intent(in) :: img2CentCell(:)

    !> Gradient contributions for each atom
    real(dp), intent(inout) :: gradients(:,:)

    real(dp) :: iAt1
    real(dp) :: sigma(3, 3)
    real(dp), allocatable :: dEdcm5(:)
    integer, allocatable :: nNeigh(:)
    real(dp), allocatable :: dhbds(:)

    @:ASSERT(self%tCoordsUpdated)
    @:ASSERT(self%tChargesUpdated)
    @:ASSERT(all(shape(gradients) == [3, self%nAtom]))

    if (allocated(self%sasaCont)) then
      call self%sasaCont%addGradients(env, neighList, species, coords, img2CentCell, gradients)
      if (allocated(self%hBondStrength)) then
        allocate(dhbds(self%nAtom))
        dhbds(:) = self%hBondStrength * self%chargesPerAtom**2
        call gemv(gradients, self%sasaCont%dsdr, dhbds, beta=1.0_dp)
        deallocate(dhbds)
      end if
    end if

    allocate(nNeigh(self%nAtom))
    sigma(:, :) = 0.0_dp
    self%energies(:) = 0.0_dp

    call getNrOfNeighboursForAll(nNeigh, neighList, self%rCutoff)
    call getBornEGCluster(self, coords, self%energies, gradients, sigma)

    if (allocated(self%cm5)) then
      allocate(dEdcm5(self%nAtom))
      dEdcm5(:) = 0.0_dp
      call hemv(dEdcm5, self%bornMat, self%chargesPerAtom)
      call self%cm5%addGradients(dEdcm5, gradients)
      call self%cm5%addSigma(dEdcm5, sigma)
    end if

    self%energies = self%energies + self%param%freeEnergyShift / real(self%nAtom, dp)

    if (self%tPeriodic) then
      self%sigma(:, :) = sigma
    end if

  end subroutine addGradients


  !> get stress tensor contributions
  subroutine getStress(self, stress)

    !> data structure
    class(TGeneralizedBorn), intent(inout) :: self

    !> Stress tensor contributions
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(self%tCoordsUpdated)
    @:ASSERT(self%tChargesUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))
    @:ASSERT(self%tPeriodic)
    @:ASSERT(self%volume > 0.0_dp)

    if (allocated(self%sasaCont)) then
      call self%sasaCont%getStress(stress)
    else
      stress(:, :) = 0.0_dp
    end if

    stress(:,:) = stress + self%sigma / self%volume

  end subroutine getStress


  !> Distance cut off for generalized Born calculations
  function getRCutoff(self) result(cutoff)

    !> data structure
    class(TGeneralizedBorn), intent(inout) :: self

    !> resulting cutoff
    real(dp) :: cutoff

    cutoff = self%rCutoff
    if (allocated(self%cm5)) then
      cutoff = max(cutoff, self%cm5%getRCutoff())
    end if

    if (allocated(self%sasaCont)) then
      cutoff = max(cutoff, self%sasaCont%getRCutoff())
    end if

  end function getRCutoff


  !> Updates with changed charges for the instance.
  subroutine updateCharges(self, env, species, neighList, qq, q0, img2CentCell, orb)

    !> Data structure
    class(TGeneralizedBorn), intent(inout) :: self

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Species, shape: [nAtom]
    integer, intent(in) :: species(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighList

    !> Orbital charges.
    real(dp), intent(in) :: qq(:,:,:)

    !> Reference orbital charges.
    real(dp), intent(in) :: q0(:,:,:)

    !> Mapping on atoms in central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    @:ASSERT(self%tCoordsUpdated)

    if (allocated(self%sasaCont)) then
      call self%sasaCont%updateCharges(env, species, neighList, qq, q0, img2CentCell, orb)
    end if

    call getSummedCharges(species, orb, qq, q0, dQAtom=self%chargesPerAtom)
    if (allocated(self%cm5)) then
      call self%cm5%addCharges(self%chargesPerAtom)
    end if

    if (allocated(self%sasaCont) .and. allocated(self%hBondStrength)) then
      self%shift(:) = 2.0_dp * self%sasaCont%sasa * self%hBondStrength * self%chargesPerAtom
    else
      self%shift(:) = 0.0_dp
    end if
    call hemv(self%shift, self%bornMat, self%chargesPerAtom, beta=1.0_dp)

    self%tChargesUpdated = .true.

  end subroutine updateCharges


  !> Returns shifts per atom
  subroutine getShifts(self, shiftPerAtom, shiftPerShell)

    !> Data structure
    class(TGeneralizedBorn), intent(inout) :: self

    !> Shift per atom
    real(dp), intent(out) :: shiftPerAtom(:)

    !> Shift per shell
    real(dp), intent(out) :: shiftPerShell(:,:)

    @:ASSERT(self%tCoordsUpdated)
    @:ASSERT(self%tChargesUpdated)
    @:ASSERT(size(shiftPerAtom) == self%nAtoms)
    @:ASSERT(size(shiftPerShell, dim=2) == self%nAtoms)

    if (allocated(self%sasaCont)) then
      call self%sasaCont%getShifts(shiftPerAtom, shiftPerShell)
    else
      shiftPerAtom(:) = 0.0_dp
      shiftPerShell(:,:) = 0.0_dp
    end if

    shiftPerAtom(:) = shiftPerAtom + self%shift

  end subroutine getShifts


  !> Calculate Born radii for a given geometry
  pure subroutine getBornRadii(self, nNeighbour, iNeighbour, img2CentCell, &
      & neighDist2, species, coords)

    !> data structure
    type(TGeneralizedBorn), intent(inout) :: self

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist
    integer, intent(in) :: iNeighbour(0:, :)

    !> Mapping into the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Square distances of the neighbours
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Central cell chemical species
    integer, intent(in) :: species(:)

    !> current atomic positions
    real(dp), intent(in) :: coords(:, :)

    integer :: iAt1, iSp1
    real(dp) :: br
    real(dp) :: dpsi
    real(dp) :: svdwi,vdwri
    real(dp) :: s1, v1, s2
    real(dp) :: arg, arg2, th, ch

    self%bornRad(:) = 0.0_dp

    call getPsi(self, nNeighbour, iNeighbour, img2CentCell, neighDist2, &
        & species, coords)

    do iAt1 = 1, self%nAtom
      iSp1 = species(iAt1)

      br = self%bornRad(iAt1)

      svdwi = self%param%vdwRad(iSp1) - self%param%bornOffset
      vdwri = self%param%vdwRad(iSp1)
      s1 = 1.0_dp/svdwi
      v1 = 1.0_dp/vdwri
      s2 = 0.5_dp*svdwi

      br = br*s2

      arg2 = br*(self%param%obc(3)*br-self%param%obc(2))
      arg = br*(self%param%obc(1)+arg2)
      arg2 = 2.0_dp*arg2+self%param%obc(1)+self%param%obc(3)*br*br

      th = tanh(arg)
      ch = cosh(arg)

      br = 1.0_dp/(s1-v1*th)
      ! Include GBMV2-like scaling
      br = self%param%bornScale*br

      dpsi = ch*(s1-v1*th)
      dpsi = s2*v1*arg2/(dpsi*dpsi)
      dpsi = self%param%bornScale*dpsi

      self%bornRad(iAt1) = br
      self%dbrdr(:, :, iAt1) = self%dbrdr(:, :, iAt1) * dpsi
      self%dbrdL(:, :, iAt1) = self%dbrdL(:, :, iAt1) * dpsi

    end do

  end subroutine getBornRadii


  !> Evaluate volume integrals, intermediate values are stored in Born radii fields
  pure subroutine getPsi(self, nNeighbour, iNeighbour, img2CentCell, &
      & neighDist2, species, coords)

    !> data structure
    class(TGeneralizedBorn), intent(inout) :: self

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist
    integer, intent(in) :: iNeighbour(0:, :)

    !> Mapping into the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Square distances of the neighbours
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Central cell chemical species
    integer, intent(in) :: species(:)

    !> current atomic positions
    real(dp), intent(in) :: coords(:, :)

    integer :: iAt1, iNeigh, iAt2, iAt2f, iSp1, iSp2
    logical :: tOvij,tOvji
    real(dp) :: vec(3),dist,rhoi,rhoj
    real(dp) :: gi,gj,ap,am,lnab,rhab,ab,dgi,dgj
    real(dp) :: dGr(3), dSr(3, 3)
    real(dp) :: rh1,rhr1,r24,rh2,r1,aprh1,r12
    real(dp) :: rvdwi,rvdwj
    real(dp), allocatable :: psi(:),dpsidr(:,:,:),dpsitr(:,:)

    allocate(psi(self%nAtom), dpsidr(3, self%nAtom, self%nAtom), dpsitr(3, self%nAtom))
    psi(:) = 0.0_dp
    dpsidr(:, :, :) = 0.0_dp
    dpsitr(:, :) = 0.0_dp

    do iAt1 = 1, self%nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeighbour(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        vec(:) = coords(:, iAt1) - coords(:, iAt2)
        dist = sqrt(neighDist2(iNeigh, iAt1))

        rhoi = self%rho(iSp1)
        rhoj = self%rho(iSp2)
        rvdwi = self%param%vdwRad(iSp1)
        rvdwj = self%param%vdwRad(iSp2)

        tOvij = dist < (rvdwi + rhoj)
        tOvji = dist < (rhoi + rvdwj)

        tOverlap: if (.not. tOvij .and. .not. tOvji) then ! ij do not overlap; ji do not overlap
          ! nonoverlaping spheres
          if(abs(rhoi-rhoj) < 1.e-8_dp) then
            ! equal reduced radiAt1
            r1 = 1.0_dp/dist
            ap = dist+rhoj
            am = dist-rhoj
            ab = ap*am
            rhab = rhoj/ab
            lnab = 0.5_dp*log(am/ap)*r1
            gi = rhab+lnab
            dgi = -2.0_dp*rhab/ab+(rhab-lnab)*r1*r1
            ! accumulate psi
            psi(iAt1) = psi(iAt1)+gi
            psi(iAt2f) = psi(iAt2f)+gi
            ! accumulate psi gradient
            dGr(:) = dgi*vec(:)
            dpsitr(:,iAt1) = dpsitr(:,iAt1)+dGr(:)
            dpsidr(:,iAt2f,iAt1) = dpsidr(:,iAt2f,iAt1)-dGr(:)
            dpsitr(:,iAt2f) = dpsitr(:,iAt2f)-dGr(:)
            dpsidr(:,iAt1,iAt2f) = dpsidr(:,iAt1,iAt2f)+dGr(:)
          else
            ! unequal reduced radiAt1
            ! ij contribution
            r1 = 1.0_dp/dist
            ap = dist+rhoj
            am = dist-rhoj
            ab = ap*am
            rhab = rhoj/ab
            lnab = 0.5_dp*log(am/ap)*r1
            gi = rhab+lnab
            dgi = -2.0_dp*rhab/ab+(rhab-lnab)*r1*r1
            ! ji contribution
            ap = dist+rhoi
            am = dist-rhoi
            ab = ap*am
            rhab = rhoi/ab
            lnab = 0.5_dp*log(am/ap)*r1
            gj = rhab+lnab
            dgj = -2.0_dp*rhab/ab+(rhab-lnab)*r1*r1
            ! accumulate psi
            psi(iAt1) = psi(iAt1)+gi
            psi(iAt2f) = psi(iAt2f)+gj
            ! accumulate psi gradient
            dGr(:) = dgi*vec(:)
            dpsitr(:,iAt1) = dpsitr(:,iAt1)+dGr(:)
            dpsidr(:,iAt2f,iAt1) = dpsidr(:,iAt2f,iAt1)-dGr(:)

            dGr(:) = dgj*vec(:)
            dpsitr(:,iAt2f) = dpsitr(:,iAt2f)-dGr(:)
            dpsidr(:,iAt1,iAt2f) = dpsidr(:,iAt1,iAt2f)+dGr(:)
          end if

        else if (.not. tOvij .and. tOvji) then tOverlap ! ij do not overlap; ji overlap

          ! ij contribution
          r1 = 1.0_dp/dist
          ap = dist+rhoj
          am = dist-rhoj
          ab = ap*am
          rhab = rhoj/ab
          lnab = 0.5_dp*log(am/ap)*r1
          gi = rhab+lnab
          dgi = -2.0_dp*rhab/ab+(rhab-lnab)*r1*r1
          ! accumulate psi
          psi(iAt1) = psi(iAt1)+gi
          ! accumulate psi gradient
          dGr(:) = dgi*vec(:)
          dpsitr(:,iAt1) = dpsitr(:,iAt1)+dGr(:)
          dpsidr(:,iAt2f,iAt1) = dpsidr(:,iAt2f,iAt1)-dGr(:)

          if((dist+rhoi) > rvdwj) then
            ! ji contribution
            r1 = 1.0_dp/dist
            r12 = 0.5_dp*r1
            r24 = r12*r12

            ap = dist+rhoi
            am = dist-rhoi
            rh1 = 1.0_dp/rvdwj
            rhr1 = 1.0_dp/ap
            aprh1 = ap*rh1
            lnab = log(aprh1)

            gj = rh1-rhr1+r12*(0.5_dp*am*(rhr1-rh1*aprh1)-lnab)

            dgj = rhr1*rhr1*(1.0_dp-0.25_dp*am*r1*(1.0_dp+aprh1*aprh1))+ &
              & rhoi*r24*(rhr1-rh1*aprh1)+r12*(r1*lnab-rhr1)
            dgj = dgj*r1
            ! accumulate psi
            psi(iAt2f) = psi(iAt2f)+gj
            ! accumulate psi gradient
            dGr(:) = dgj*vec(:)
            dpsitr(:,iAt2f) = dpsitr(:,iAt2f)-dGr(:)
            dpsidr(:,iAt1,iAt2f) = dpsidr(:,iAt1,iAt2f)+dGr(:)
          end if

        else if (tOvij .and. .not. tOvji) then ! ij overlap; ji do not overlap

          if((dist+rhoj) > rvdwi) then
            ! ij contribution
            r1 = 1.0_dp/dist
            r12 = 0.5_dp*r1
            r24 = r12*r12

            ap = dist+rhoj
            am = dist-rhoj
            rh1 = 1.0_dp/rvdwi
            rhr1 = 1.0_dp/ap
            aprh1 = ap*rh1
            lnab = log(aprh1)

            gi = rh1-rhr1+r12*(0.5_dp*am*(rhr1-rh1*aprh1)-lnab)

            dgi = rhr1*rhr1*(1.0_dp-0.25_dp*am*r1*(1.0_dp+aprh1*aprh1))+ &
              & rhoj*r24*(rhr1-rh1*aprh1)+r12*(r1*lnab-rhr1)
            dgi = dgi*r1
            ! accumulate psi
            psi(iAt1) = psi(iAt1)+gi
            ! accumulate psi gradient
            dGr(:) = dgi*vec(:)
            dpsitr(:,iAt1) = dpsitr(:,iAt1)+dGr(:)
            dpsidr(:,iAt2f,iAt1) = dpsidr(:,iAt2f,iAt1)-dGr(:)
          end if

          ! ji contribution
          ap = dist+rhoi
          am = dist-rhoi
          ab = ap*am
          rhab = rhoi/ab
          lnab = 0.5_dp*log(am/ap)*r1
          gj = rhab+lnab
          dgj = -2.0_dp*rhab/ab+(rhab-lnab)*r1*r1
          ! accumulate psi
          psi(iAt2f) = psi(iAt2f)+gj
          ! accumulate psi gradient
          dGr(:) = dgj*vec(:)
          dpsitr(:,iAt2f) = dpsitr(:,iAt2f)-dGr(:)
          dpsidr(:,iAt1,iAt2f) = dpsidr(:,iAt1,iAt2f)+dGr(:)

        else if (tOvij .and. tOvji) then tOverlap ! ij and ji overlap
          ! overlaping spheres
          if((dist+rhoj) > rvdwi) then
            ! ij contribution
            r1 = 1.0_dp/dist
            r12 = 0.5_dp*r1
            r24 = r12*r12

            ap = dist+rhoj
            am = dist-rhoj
            rh1 = 1.0_dp/rvdwi
            rhr1 = 1.0_dp/ap
            aprh1 = ap*rh1
            lnab = log(aprh1)

            gi = rh1-rhr1+r12*(0.5_dp*am*(rhr1-rh1*aprh1)-lnab)

            dgi = rhr1*rhr1*(1.0_dp-0.25_dp*am*r1*(1.0_dp+aprh1*aprh1))+ &
              & rhoj*r24*(rhr1-rh1*aprh1)+r12*(r1*lnab-rhr1)
            dgi = dgi*r1
            ! accumulate psi
            psi(iAt1) = psi(iAt1)+gi
            ! accumulate psi gradient
            dGr(:) = dgi*vec(:)
            dpsitr(:,iAt1) = dpsitr(:,iAt1)+dGr(:)
            dpsidr(:,iAt2f,iAt1) = dpsidr(:,iAt2f,iAt1)-dGr(:)
          end if

          if((dist+rhoi) > rvdwj) then
            ! ji contribution
            r1 = 1.0_dp/dist
            r12 = 0.5_dp*r1
            r24 = r12*r12

            ap = dist+rhoi
            am = dist-rhoi
            rh1 = 1.0_dp/rvdwj
            rhr1 = 1.0_dp/ap
            aprh1 = ap*rh1
            lnab = log(aprh1)

            gj = rh1-rhr1+r12*(0.5_dp*am*(rhr1-rh1*aprh1)-lnab)

            dgj = rhr1*rhr1*(1.0_dp-0.25_dp*am*r1*(1.0_dp+aprh1*aprh1))+ &
              & rhoi*r24*(rhr1-rh1*aprh1)+r12*(r1*lnab-rhr1)
            dgj = dgj*r1
            ! accumulate psi
            psi(iAt2f) = psi(iAt2f)+gj
            ! accumulate psi gradient
            dGr(:) = dgj*vec(:)
            dpsitr(:,iAt2f) = dpsitr(:,iAt2f)-dGr(:)
            dpsidr(:,iAt1,iAt2f) = dpsidr(:,iAt1,iAt2f)+dGr(:)
          end if

        end if tOverlap

      end do
    end do

    ! save Born radii
    self%bornRad(:) = psi
    ! save derivative of Born radii w.r.t. atomic positions
    self%dbrdr(:,:,:) = dpsidr
    ! save one-center terms
    do iAt1 = 1, self%nAtom
      self%dbrdr(:,iAt1,iAt1) = self%dbrdr(:,iAt1,iAt1) + dpsitr(:,iAt1)
    end do

  end subroutine getPsi


  !> compute Born matrix
  pure subroutine getBornMatrixCluster(self, coords0)

    !> data structure
    type(TGeneralizedBorn), intent(inout) :: self

    !> coordinates in the central cell
    real(dp), intent(in) :: coords0(:, :)

    integer :: iAt1, iAt2, iAt2f, iNeigh
    real(dp) :: aa, dist2, dd, expd, dfgb, fgb

    self%bornMat(:, :) = 0.0_dp

    do iAt1 = 1, self%nAtom
       do iAt2 = 1, iAt1-1
          dist2 = sum((coords0(:, iAt1) - coords0(:, iAt2))**2)

          aa = self%bornRad(iAt1)*self%bornRad(iAt2)
          dd = 0.25_dp*dist2/aa
          expd = exp(-dd)
          dfgb = 1.0_dp/(dist2+aa*expd)
          fgb = self%param%keps*sqrt(dfgb)
          self%bornMat(iAt1, iAt2) = self%bornMat(iAt1, iAt2) + fgb
          self%bornMat(iAt2, iAt1) = self%bornMat(iAt2, iAt1) + fgb
       end do
    end do

    !> self-energy part
    do iAt1 = 1, self%nAtom
       self%bornMat(iAt1, iAt1) = self%param%keps/self%bornRad(iAt1)
    end do

  end subroutine getBornMatrixCluster


  !> GB energy and gradient
  subroutine getBornEGCluster(self, coords, energies, gradients, sigma)

    !> data structure
    type(TGeneralizedBorn), intent(in) :: self

    !> Current atomic positions
    real(dp), intent(in) :: coords(:, :)

    !> Atom resolved energies
    real(dp), intent(out) :: energies(:)

    !> Molecular gradient
    real(dp), intent(inout) :: gradients(:, :)

    !> Strain derivative
    real(dp), intent(inout) :: sigma(:, :)

    integer :: iAt1, iAt2
    real(dp) :: aa, dist2, fgb, fgb2, qq, dd, expd, dfgb, dfgb2, dfgb3, ap, bp
    real(dp) :: grddbi,grddbj, vec(3), dGr(3), dSr(3, 3)
    real(dp), allocatable :: dEdbr(:)
    real(dp), allocatable :: derivs(:, :)

    allocate(dEdbr(self%nAtom), derivs(3, self%nAtom))

    derivs(:, :) = 0.0_dp
    energies(:) = 0.0_dp
    dEdbr(:) = 0.0_dp

    do iAt1 = 1, self%nAtom
       do iAt2 = 1, iAt1-1
          vec(:) = coords(:, iAt1) - coords(:, iAt2)
          dist2 = sum(vec**2)

          ! dielectric scaling of the charges
          qq = self%chargesPerAtom(iAt1)*self%chargesPerAtom(iAt2)
          aa = self%bornRad(iAt1)*self%bornRad(iAt2)
          dd = 0.25_dp*dist2/aa
          expd = exp(-dd)
          fgb2 = dist2+aa*expd
          dfgb2 = 1.0_dp/fgb2
          dfgb = sqrt(dfgb2)
          dfgb3 = dfgb2*dfgb*self%param%keps

          energies(iAt1) = energies(iAt1) + qq*self%param%keps*dfgb/2
          if (iAt1 /= iAt2) then
             energies(iAt2) = energies(iAt2) + qq*self%param%keps*dfgb/2
          end if

          ap = (1.0_dp-0.25_dp*expd)*dfgb3
          dGr = ap*vec
          derivs(:,iAt1) = derivs(:,iAt1) - dGr*qq
          derivs(:,iAt2) = derivs(:,iAt2) + dGr*qq

          dSr = spread(dGr, 1, 3) * spread(vec, 2, 3)
          if (iAt1 /= iAt2) then
             sigma = sigma + dSr
          else
             sigma = sigma + dSr/2
          end if

          bp = -0.5_dp*expd*(1.0_dp+dd)*dfgb3
          grddbi = self%bornRad(iAt2)*bp
          grddbj = self%bornRad(iAt1)*bp
          dEdbr(iAt1) = dEdbr(iAt1) + grddbi*qq
          if (iAt1 /= iAt2) then
             dEdbr(iAt2) = dEdbr(iAt2) + grddbj*qq
          end if

       end do
    end do

    gradients(:, :) = gradients + derivs

    !> self-energy part
    do iAt1 = 1, self%nAtom
       bp = 1.0_dp/self%bornRad(iAt1)
       qq = self%chargesPerAtom(iAt1)*bp
       energies(iAt1) = energies(iAt1) + 0.5_dp*self%chargesPerAtom(iAt1)*qq*self%param%keps
       grddbi = -0.5_dp*self%param%keps*qq*bp
       dEdbr(iAt1) = dEdbr(iAt1) + grddbi*self%chargesPerAtom(iAt1)
    end do

    !> contract with the Born radii derivatives
    call gemv(gradients, self%dbrdr, dEdbr, beta=1.0_dp)
    call gemv(sigma, self%dbrdL, dEdbr, beta=1.0_dp)

  end subroutine getBornEGCluster


end module dftbp_born
