!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module dftbp_solvation_cosmo
  use ddcosmo_core, only : ddupdate, TDomainDecomposition_init, wghpot, intrhs, fdoka, fdokb,&
      & fdoga
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi, Hartree__eV, Bohr__AA
  use dftbp_common_environment, only : TEnvironment
  use dftbp_dftb_charges, only : getSummedCharges
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_extlibs_ddcosmo, only : TDomainDecompositionInput, TDomainDecomposition, jacobi_diis,&
      & lx, lstarx, ldm1x, hnorm
  use dftbp_extlibs_lebedev, only : getAngGrid, gridSize
  use dftbp_io_charmanip, only : tolower
  use dftbp_io_message, only : error
  use dftbp_math_blasroutines, only : gemv
  use dftbp_solvation_sasa, only : TSASACont, TSASAInput, TSASACont_init, writeSASAContInfo
  use dftbp_solvation_solvation, only : TSolvation
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: TCosmo, TCosmoInput, TCosmo_init, writeCosmoInfo
  public :: TDomainDecompositionInput


  !> Input for the conductor like screening model
  type :: TCosmoInput

    !> Energy shift to the reference system
    real(dp) :: freeEnergyShift

    !> Dielectric constant
    real(dp) :: dielectricConst

    !> dielectric scaling factor on going between inside cavity (vacuum) and outside region
    real(dp) :: keps

    !> Grid for numerical integration of atomic surfaces
    integer :: gridSize

    !> Van-der-Waals radii
    real(dp), allocatable :: vdwRad(:)

    !> Input for the domain decomposition algorithm
    type(TDomainDecompositionInput) :: ddInput

    !> Input for solvent accessible surface area model
    type(TSASAInput), allocatable :: sasaInput

  end type TCosmoInput


  !> Conductor like screening model
  type, extends(TSolvation) :: TCosmo

    !> Number of atoms
    integer :: nAtom

    !> Domain decomposition COSMO solver
    type(TDomainDecomposition) :: ddCosmo

    !> Electrostatic potential phi(ncav)
    real(dp), allocatable :: phi(:)

    !> Psi vector psi(nylm, nAtom)
    real(dp), allocatable :: psi(:, :)

    !> ddcosmo solution sigma(nylm, nAtom)
    real(dp), allocatable :: sigma(:, :)

    !> ddcosmo adjoint solution s(nylm, nAtom)
    real(dp), allocatable :: s(:, :)

    !> Dielectric constant
    real(dp) :: dielectricConst

    !> scaling factor from dielectric inside cavity (vacuum) and outside region
    real(dp) :: keps

    !> Energy shift to the reference system
    real(dp) :: freeEnergyShift

    !> Van-der-Waal radii
    real(dp), allocatable :: vdwRad(:)

    !> Angular grid for surface integration
    real(dp), allocatable :: angGrid(:, :)

    !> Weights of grid points for surface integration
    real(dp), allocatable :: angWeight(:)

    !> Interaction matrix with surface charges jmat(ncav, nAtom)
    real(dp), allocatable :: jmat(:, :)

    !> Partial charges per atom
    real(dp), allocatable :: chargesPerAtom(:)

    !> Potential shift per atom
    real(dp), allocatable :: shiftsPerAtom(:)

    !> Charges are current
    logical :: tChargesUpdated

    !> Coordinates are current
    logical :: tCoordsUpdated

    !> Solvent accessible surface area model
    type(TSASACont), allocatable :: sasaCont

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

    !> Updates with changed charges for the instance.
    procedure :: updateCharges

    !> Returns shifts per atom
    procedure :: getShifts

    !> Is the electrostic field modified by this solvent model?
    procedure :: isEFieldModified

    !> Write cavity information
    procedure :: writeCosmoFile

  end type TCosmo


  real(dp), parameter :: fourpi = 4.0_dp * pi
  integer, parameter :: ndiis=25


contains


  subroutine TCosmo_init(this, input, nAtom, species0, speciesNames, latVecs)

    !> Instance of the solvation model
    type(TCosmo), intent(out) :: this

    !> Input to setup the solvation model
    type(TCosmoInput), intent(in) :: input

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    integer :: iat, stat

    this%tChargesUpdated = .false.
    this%tCoordsUpdated = .false.

    this%nAtom = nAtom
    this%dielectricConst = input%dielectricConst
    this%keps = input%keps
    this%freeEnergyShift = input%freeEnergyShift

    allocate(this%vdwRad(nAtom))
    do iat = 1, nAtom
      this%vdwRad(iat) = input%vdwRad(species0(iat))
    end do

    allocate(this%angGrid(3, gridSize(input%gridSize)))
    allocate(this%angWeight(gridSize(input%gridSize)))
    call getAngGrid(input%gridSize, this%angGrid, this%angWeight, stat)
    if (stat /= 0) then
      call error("Could not initialize angular grid for COSMO model")
    end if

    allocate(this%chargesPerAtom(nAtom))
    allocate(this%shiftsPerAtom(nAtom))

    call TDomainDecomposition_init(this%ddCosmo, input%ddInput, &
      this%vdwRad, this%angWeight, this%angGrid)

    if (allocated(input%sasaInput)) then
       allocate(this%sasaCont)
       call TSASACont_init(this%sasaCont, input%sasaInput, nAtom, species0, &
           & speciesNames, latVecs)
    end if

  end subroutine TCosmo_init


  subroutine writeCosmoInfo(unit, solvation)

    !> Formatted unit for IO
    integer, intent(in) :: unit

    !> Solvation model
    type(TCosmo), intent(in) :: solvation

    write(unit, '(a, ":", t30, es14.6)') "Dielectric constant", &
        & solvation%dielectricConst
    write(unit, '(a, ":", t30, es14.6, 1x, a, t50, es14.6, 1x, a)') &
        & "Free energy shift", solvation%freeEnergyShift, "H", &
        & Hartree__eV * solvation%freeEnergyShift, "eV"
    write(unit, '(a, ":", t30, es14.6, 1x, a, t50, i14, 1x, a)') "Grid points", &
        & solvation%nAtom*real(size(solvation%angWeight, dim=1), dp), "total", &
        & size(solvation%angWeight, dim=1), "per atom"
    write(unit, '(a, ":", t30, a)') "Solver", "domain decomposition"

    write(unit, '(a, ":", t30)', advance='no') "SASA model"
    if (allocated(solvation%sasaCont)) then
      write(unit, '(a)') "Yes"
      call writeSASAContInfo(unit, solvation%sasaCont)
    else
      write(unit, '(a)') "No"
    end if

  end subroutine writeCosmoInfo


  !> Update internal stored coordinates
  subroutine updateCoords(this, env, neighList, img2CentCell, coords, species0)

    !> Data structure
    class(TCosmo), intent(inout) :: this

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

    if (allocated(this%sasaCont)) then
      call this%sasaCont%updateCoords(env, neighList, img2CentCell, coords, species0)
    end if

    if (allocated(this%phi)) deallocate(this%phi)
    if (allocated(this%psi)) deallocate(this%psi)
    if (allocated(this%sigma)) deallocate(this%sigma)
    if (allocated(this%s)) deallocate(this%s)
    if (allocated(this%jmat)) deallocate(this%jmat)

    call ddupdate(this%ddcosmo, coords)

    allocate(this%phi(this%ddCosmo%ncav), this%psi(this%ddCosmo%nylm, this%nAtom))
    allocate(this%jmat(this%ddCosmo%ncav, this%nAtom))

    call getCoulombMatrix(coords(:, :this%nAtom), this%ddCosmo%ccav, this%jmat)

    this%tChargesUpdated = .false.
    this%tCoordsUpdated = .true.

  end subroutine updateCoords


  !> Update internal copy of lattice vectors
  subroutine updateLatVecs(this, latVecs)

    !> Data structure
    class(TCosmo), intent(inout) :: this

    !> Lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    if (allocated(this%sasaCont)) then
      call this%sasaCont%updateLatVecs(latVecs)
    end if

    this%tChargesUpdated = .false.
    this%tCoordsUpdated = .false.

  end subroutine updateLatVecs


  !> Get energy contributions
  subroutine getEnergies(this, energies)

    !> Data structure
    class(TCosmo), intent(inout) :: this

    !> Energy contributions for each atom
    real(dp), intent(out) :: energies(:)

    integer :: iat

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(size(energies) == this%nAtom)

    if (allocated(this%sasaCont)) then
      call this%sasaCont%getEnergies(energies)
    else
      energies(:) = 0.0_dp
    end if

    do iat = 1, size(energies)
      energies(iat) = this%keps * dot_product(this%sigma(:, iat), this%psi(:, iat)) &
         & + this%freeEnergyShift / real(this%nAtom, dp) + energies(iat)
    end do

  end subroutine getEnergies


  !> Get force contributions
  subroutine addGradients(this, env, neighList, species, coords, img2CentCell, gradients)

    !> Data structure
    class(TCosmo), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighList

    !> Specie for each atom.
    integer, intent(in) :: species(:)

    !> Coordinate of each atom.
    real(dp), intent(in) :: coords(:,:)

    !> Mapping of atoms to cetnral cell.
    integer, intent(in) :: img2CentCell(:)

    !> Gradient contributions for each atom
    real(dp), intent(inout) :: gradients(:,:)

    integer :: ii, iat, ig
    real(dp) :: xx(1), esolv
    real(dp), allocatable :: fx(:, :), zeta(:), ef1(:, :), ef2(:, :)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    if (allocated(this%sasaCont)) then
      call this%sasaCont%addGradients(env, neighList, species, coords, img2CentCell, gradients)
    end if

    allocate(fx(3, this%nAtom), zeta(this%ddCosmo%ncav), &
      & ef1(3, this%ddCosmo%ncav), ef2(3, this%nAtom))

    ! reset Psi
    call getPsi(this%chargesPerAtom, this%psi)

    call solveCosmoAdjoint(this%ddCosmo, this%psi, this%s, .true., &
      & accuracy=this%ddCosmo%conv*1e-3_dp)

    ! reset Phi
    call getPhi(this%chargesPerAtom, this%jmat, this%phi)

    ! now call the routine that computes the ddcosmo specific contributions
    ! to the forces.
    call forces(this%ddCosmo, this%keps, this%phi, this%sigma, this%s, fx)

    ! form the "zeta" intermediate
    call getZeta(this%ddCosmo, this%keps, this%s, zeta)

    ! 1. solute's electric field at the cav points times zeta:
    !    compute the electric field
    call efld(this%nAtom, this%chargesPerAtom, this%ddCosmo%xyz, this%ddCosmo%ncav, &
      & this%ddCosmo%ccav, ef1)

    ! contract it with the zeta intermediate
    ii = 0
    do iat = 1, this%nAtom
      do ig = 1, size(this%angWeight)
        if (this%ddCosmo%ui(ig, iat) > 0.0_dp) then
          ii = ii + 1
          fx(:, iat) = fx(:, iat) - zeta(ii)*ef1(:, ii)
        end if
      end do
    end do

    @:ASSERT(ii == this%ddCOSMO%ncav)

    ! 2. "zeta's" electric field at the nuclei times the charges.
    !    compute the "electric field"
    call efld(this%ddCosmo%ncav, zeta, this%ddCosmo%ccav, this%nAtom, &
      & this%ddCosmo%xyz, ef2)

    ! contract it with the solute's charges.
    do iat = 1, this%nAtom
      fx(:, iat) = fx(:, iat) - ef2(:, iat)*this%chargesPerAtom(iat)
    end do

    gradients(:, :) = gradients(:, :) - fx

  end subroutine addGradients


  !> Get stress tensor contributions
  subroutine getStress(this, stress)

    !> Data structure
    class(TCosmo), intent(inout) :: this

    !> Stress tensor contributions
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))

    if (allocated(this%sasaCont)) then
      call this%sasaCont%getStress(stress)
    else
      stress(:, :) = 0.0_dp
    end if

  end subroutine getStress


  !> Distance cut off for dispersion interactions
  function getRCutoff(this) result(cutoff)

    !> Data structure
    class(TCosmo), intent(inout) :: this

    !> Resulting cutoff
    real(dp) :: cutoff

    if (allocated(this%sasaCont)) then
      cutoff = this%sasaCont%getRCutoff()
    else
      cutoff = 0.0_dp
    end if

  end function getRCutoff


  !> Updates with changed charges for the instance.
  subroutine updateCharges(this, env, species, neighList, qq, q0, img2CentCell, orb)

    !> Data structure
    class(TCosmo), intent(inout) :: this

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

    real(dp) :: xx(1, 1)
    logical :: restart

    @:ASSERT(this%tCoordsUpdated)

    if (allocated(this%sasaCont)) then
      call this%sasaCont%updateCharges(env, species, neighList, qq, q0, img2CentCell, orb)
    end if

    restart = allocated(this%sigma)
    if (.not.allocated(this%sigma)) then
      allocate(this%sigma(this%ddCosmo%nylm, this%nAtom))
    end if

    call getSummedCharges(species, orb, qq, q0, dQAtom=this%chargesPerAtom)

    call getPhi(this%chargesPerAtom, this%jmat, this%phi)

    call solveCosmoDirect(this%ddCosmo, .true., this%phi, xx, this%sigma, restart)

    restart = allocated(this%s)
    if (.not.allocated(this%s)) then
      allocate(this%s(this%ddCosmo%nylm, this%nAtom))
    end if

    call getPsi(this%chargesPerAtom, this%psi)

    ! solve adjoint ddCOSMO equation to get full potential contributions
    call solveCosmoAdjoint(this%ddCosmo, this%psi, this%s, restart)

    this%tChargesUpdated = .true.

  end subroutine updateCharges


  !> Returns shifts per atom
  subroutine getShifts(this, shiftPerAtom, shiftPerShell)

    !> Data structure
    class(TCosmo), intent(inout) :: this

    !> Shift per atom
    real(dp), intent(out) :: shiftPerAtom(:)

    !> Shift per shell
    real(dp), intent(out) :: shiftPerShell(:,:)

    real(dp) :: xx(1)
    logical :: restart

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(size(shiftPerAtom) == this%nAtom)
    @:ASSERT(size(shiftPerShell, dim=2) == this%nAtom)

    if (allocated(this%sasaCont)) then
      call this%sasaCont%getShifts(shiftPerAtom, shiftPerShell)
    else
      shiftPerAtom(:) = 0.0_dp
      shiftPerShell(:,:) = 0.0_dp
    end if

    shiftPerAtom(:) = this%keps * this%sigma(1, :) * sqrt(fourpi) + shiftPerAtom

    ! we abuse Phi to store the unpacked and scaled value of s
    call getZeta(this%ddCosmo, this%keps, this%s, this%phi)
    ! and contract with the Coulomb matrix
    call gemv(shiftPerAtom, this%jmat, this%phi, alpha=-1.0_dp, &
      & beta=1.0_dp, trans='t')

  end subroutine getShifts


  !> Is the electrostic field modified by this solvent model?
  function isEFieldModified(this) result(isChanged)

    !> Data structure
    class(TCosmo), intent(in) :: this

    !> Has the solvent model changed the electrostatic environment
    logical :: isChanged

    isChanged = .false.

  end function isEFieldModified

  !> Evaluate the Coulomb interactions between the atomic sides (xyz) and the
  !> surface elements of the cavity (ccav).
  subroutine getCoulombMatrix(xyz, ccav, jmat)
    real(dp), intent(in) :: xyz(:, :)
    real(dp), intent(in) :: ccav(:, :)
    real(dp), intent(inout) :: jmat(:, :)

    integer :: ic, j
    real(dp) :: vec(3), d2, d

    jmat(:, :) = 0.0_dp
    !$omp parallel do default(none) schedule(runtime) collapse(2) &
    !$omp shared(ccav, xyz, jmat) private(ic, j, vec, d2, d)
    do ic = 1, size(ccav, 2)
      do j = 1, size(xyz, 2)
        vec(:) = ccav(:, ic) - xyz(:, j)
        d2 = vec(1)**2 + vec(2)**2 + vec(3)**2
        d = sqrt(d2)
        jmat(ic, j) = 1.0_dp / d
      end do
    end do

  end subroutine getCoulombMatrix


  !> Routine to compute the psi vector
  subroutine getPsi(charge, psi)
    real(dp), intent(in) :: charge(:)
    real(dp), intent(out) :: psi(:, :)

    integer :: iat
    real(dp) :: fac

    fac = sqrt(fourpi)
    psi(:,:) = 0.0_dp

    do iat = 1, size(charge)
      psi(1, iat) = fac*charge(iat)
    end do

  end subroutine getPsi


  !> Routine to compute the potential vector
  subroutine getPhi(charge, jmat, phi)
    real(dp), intent(in) :: charge(:)
    real(dp), intent(in) :: jmat(:, :)
    real(dp), intent(out) :: phi(:)

    @:ASSERT(size(jmat, 1) == size(phi))
    @:ASSERT(size(jmat, 2) == size(charge))

    phi(:) = 0.0_dp

    call gemv(phi, jmat, charge)

  end subroutine getPhi


  !> Wrapper for the linear solvers for COSMO equation
  !    L sigma = G
  !  This routine performs the following operations :
  !   - allocates memory for the linear solvers
  !   - if star is false and cart is true, assembles the right-hand side for the COSMO
  !     equations.
  !   - computes a guess for the solution (using the inverse diagonal);
  !   - calls the iterative solver;
  subroutine solveCosmoDirect(ddCosmo, cart, phi, glm, sigma, restart)

    !> Error source
    character(len=*), parameter :: source = 'cosmo::solveCosmoDirect'
    type(TDomainDecomposition), intent(in) :: ddCosmo

    !> true:  the right-hand side for the COSMO has to be assembled
    !         inside this routine and the unscaled potential at the
    !         external points of the cavity is provided in phi.
    !  false: the right-hand side for the COSMO equations is provided
    !         in glm.
    logical, intent(in) :: cart

    !> Contains the potential at the external cavity points if cart is true.
    !  phi is not referenced in any other case.
    real(dp), intent(in) :: phi(:)

    !> Contains the right-hand side for the COSMO equations if cart is false.
    !  glm is not referenced in any other case
    real(dp), intent(in) :: glm(:, :)

    !> The solution to the COSMO (adjoint) equations
    real(dp), intent(inout) :: sigma(:, :)

    !> Initial guess is provided on sigma
    logical, intent(in) :: restart

    integer :: iat, istatus, n_iter, info, c1, c2, cr
    real(dp) :: tol, r_norm
    logical :: ok

    real(dp), allocatable  :: g(:, :), rhs(:, :), work(:, :)

    ! parameters for the solver and matvec routine
    tol     = ddCosmo%conv
    n_iter  = 200

    ! DIRECT COSMO EQUATION L X = g

    ! allocate workspace for rhs
    allocate(rhs(ddCosmo%nylm, ddCosmo%nat), stat=istatus)
    if (istatus /= 0) then
      write(*, *) ' cosmo: [2] failed allocation'
    endif

    ! 1. RHS
    ! assemble rhs
    if (cart) then

      ! allocate workspace for weighted potential
      allocate(g(ddCosmo%ngrid, ddCosmo%nat) , stat=istatus)
      if (istatus /= 0) then
        write(*, *) ' cosmo: [3] failed allocation'
      endif

      ! weight the potential...
      call wghpot(ddCosmo, phi, g)

      ! ... and compute its multipolar expansion
      do iat = 1, ddCosmo%nat
        call intrhs(ddCosmo, iat, g(:, iat), rhs(:, iat))
      enddo

      ! deallocate workspace
      deallocate(g , stat=istatus)
      if (istatus /= 0) then
        write(*, *) 'cosmo: [1] failed deallocation'
      endif

    else
      ! no need to manipulate rhs
      rhs = glm
    end if

    ! 2. INITIAL GUESS
    if (.not.restart) then
      do iat = 1, ddCosmo%nat
        sigma(:, iat) = ddCosmo%facl(:)*rhs(:, iat)
      end do
    end if

    ! 3. SOLVER CALL
    ! Jacobi method :
    ! L X = (diag + offdiag) X = g   ==>    X = diag^-1 (g - offdiag X_guess)
    ! action of  diag^-1 :  ldm1x
    ! action of  offdiag :  lx
    call jacobi_diis(ddCosmo, ddCosmo%nat*ddCosmo%nylm, ddCosmo%iprint, &
      & ndiis, 4, tol, rhs, sigma, n_iter, ok, lx, ldm1x, hnorm)

    ! check solution
    if (.not.ok) then
      call error('direct ddCOSMO did not converge!')
      return
    endif

  end subroutine solveCosmoDirect


  !> Wrapper for the linear solvers for adjoint COSMO equation
  !>   L^* sigma = Psi
  !> This routine performs the following operations :
  !>  - allocates memory for the linear solvers
  !>  - computes a guess for the solution (using the inverse diagonal);
  !>  - calls the iterative solver;
  subroutine solveCosmoAdjoint(ddCosmo, psi, sigma, restart, accuracy)
    ! Error source
    character(len=*), parameter :: source = 'cosmo::solveCosmoAdjoint'
    type(TDomainDecomposition), intent(in) :: ddCosmo

    !> The psi vector. it is used as a right-hand side
    real(dp), intent(in) :: psi(:, :)

    !> The solution to the COSMO (adjoint) equations
    real(dp), intent(inout) :: sigma(:, :)

    !> Initial guess is provided on sigma
    logical, intent(in) :: restart

    !> Overwrite accuracy
    real(dp), intent(in), optional :: accuracy

    integer :: iat, istatus, n_iter, info, c1, c2, cr
    real(dp) :: tol, r_norm
    logical :: ok

    real(dp), allocatable  :: g(:, :), rhs(:, :), work(:, :)

    ! parameters for the solver and matvec routine
    if (present(accuracy)) then
      tol = accuracy
    else
      tol = ddCosmo%conv
    end if
    n_iter  = 200

    ! 1. INITIAL GUESS
    if (.not.restart) then
      do iat = 1, ddCosmo%nat
        sigma(:, iat) = ddCosmo%facl(:)*psi(:, iat)
      end do
    end if

    ! 2. SOLVER CALL
    ! Jacobi method : see above
    call jacobi_diis(ddCosmo, ddCosmo%nat*ddCosmo%nylm, ddCosmo%iprint, &
      & ndiis, 4, tol, psi, sigma, n_iter, ok, lstarx, ldm1x, hnorm)

    ! check solution
    if (.not.ok) then
      call error('adjoint ddCOSMO did not converge!')
      return
    endif

  end subroutine solveCosmoAdjoint


  !> Compute
  !
  ! \zeta(n, i) =
  !
  !  1/2 f(\eps) sum w_n U_n^i Y_l^m(s_n) [S_i]_l^m
  !              l, m
  !
  subroutine getZeta(ddCosmo, keps, s, zeta)
    type(TDomainDecomposition), intent(in) :: ddCosmo
    real(dp), intent(in) :: keps
    real(dp), intent(in) :: s(:, :) ! [ddCosmo%nylm, ddCosmo%nat]
    real(dp), intent(inout) :: zeta(:) ! [ddCosmo%ncav]

    integer :: its, iat, ii

    @:ASSERT(size(s, 1) == ddCosmo%nylm)
    @:ASSERT(size(s, 2) == ddCosmo%nat)
    @:ASSERT(size(zeta) == ddCosmo%ncav)

    ii = 0
    do iat = 1, ddCosmo%nat
      do its = 1, ddCosmo%ngrid
        if (ddCosmo%ui(its, iat) > 0.0_dp) then
          ii = ii + 1
          zeta(ii) = keps * ddCosmo%w(its) * ddCosmo%ui(its, iat) &
            & * dot_product(ddCosmo%basis(:, its), s(:, iat))
        end if
      end do
    end do

    @:ASSERT(ii == ddCOSMO%ncav)

  end subroutine getZeta


  !> Sample driver for the calculation of the ddCOSMO forces.
  subroutine forces(ddCosmo, keps, phi, sigma, s, fx)
    type(TDomainDecomposition), intent(in) :: ddCosmo
    real(dp), intent(in) :: keps
    real(dp), intent(in) :: phi(:)
    real(dp), intent(in) :: sigma(:, :)
    real(dp), intent(in) :: s(:, :)
    real(dp), intent(inout) :: fx(:, :)

    integer :: iat, ig, ii, c1, c2, cr
    real(dp) :: fep

    real(dp), allocatable :: xi(:, :), phiexp(:, :), zeta(:), ef(:, :)
    real(dp), allocatable :: basloc(:), dbsloc(:, :), vplm(:), vcos(:), vsin(:)

    @:ASSERT(size(phi) == ddCosmo%ncav)
    @:ASSERT(size(sigma, 1) == ddCosmo%nylm)
    @:ASSERT(size(sigma, 2) == ddCosmo%nat)
    @:ASSERT(size(s, 1) == ddCosmo%nylm)
    @:ASSERT(size(s, 2) == ddCosmo%nat)
    @:ASSERT(size(fx, 2) == ddCosmo%nat)

    allocate (xi(ddCosmo%ngrid, ddCosmo%nat), phiexp(ddCosmo%ngrid, ddCosmo%nat))
    allocate (basloc(ddCosmo%nylm), dbsloc(3, ddCosmo%nylm), vplm(ddCosmo%nylm), &
      & vcos(ddCosmo%lmax+1), vsin(ddCosmo%lmax+1))

    ! compute xi:
    !$omp parallel do default(none) collapse(2) schedule(runtime) &
    !$omp shared(ddCosmo, s, xi) private(iat, ig)
    do iat = 1, ddCosmo%nat
      do ig = 1, ddCosmo%ngrid
        xi(ig, iat) = dot_product(s(:, iat), ddCosmo%basis(:, ig))
      end do
    end do

    ! expand the potential on a sphere-by-sphere basis (needed for parallelism):
    ii = 0
    phiexp(:, :) = 0.0_dp
    do iat = 1, ddCosmo%nat
      do ig = 1, ddCosmo%ngrid
        if (ddCosmo%ui(ig, iat) > 0.0_dp) then
          ii = ii + 1
          phiexp(ig, iat) = phi(ii)
        end if
      end do
    end do

    fx(:, :) = 0.0_dp
    do iat = 1, ddCosmo%nat
      call fdoka(ddCosmo, iat, sigma, xi(:, iat), basloc, dbsloc, vplm, &
        & vcos, vsin, fx(:, iat))
      call fdokb(ddCosmo, iat, sigma, xi, basloc, dbsloc, vplm, vcos, vsin, &
        & fx(:, iat))
      call fdoga(ddCosmo, iat, xi, phiexp, fx(:, iat))
    end do

    deallocate (basloc, dbsloc, vplm, vcos, vsin)
    deallocate (xi, phiexp)

    ! scale the forces time the cosmo factor:
    fx  = keps*fx

  end subroutine forces


  !> Computes the electric field produced by the sources src (nsrc point charges
  !  with coordinates csrc) at the ntrg target points ctrg:
  subroutine efld(nsrc, src, csrc, ntrg, ctrg, ef)
    integer, intent(in) :: nsrc, ntrg
    real(dp), intent(in) :: src(:)
    real(dp), intent(in) :: csrc(:, :)
    real(dp), intent(in) :: ctrg(:, :)
    real(dp), intent(inout) :: ef(:, :)

    integer :: i, j
    real(dp) :: vec(3), r2, rr, r3, f
    real(dp), parameter :: zero=0.0_dp

    @:ASSERT(size(src) == size(csrc, 2))
    @:ASSERT(size(ef, 2) == size(ctrg, 2))

    ef(:, :) = 0.0_dp
    !$omp parallel do default(none) schedule(runtime) collapse(2) &
    !$omp reduction(+:ef) shared(ntrg, nsrc, ctrg, csrc, src) &
    !$omp private(j, i, f, vec, r2, rr, r3)
    do j = 1, ntrg
      do i = 1, nsrc
        vec(:) = ctrg(:, j) - csrc(:, i)
        r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
        rr = sqrt(r2)
        r3 = r2*rr
        f = src(i)/r3
        ef(:, j) = ef(:, j) + f*vec
      end do
    end do

  end subroutine efld


  !> Write a COSMO file output
  subroutine writeCosmoFile(this, unit, species0, speciesNames, coords0, energy)

    !> COSMO container
    class(TCosmo), intent(in) :: this

    !> Formatted unit for output
    integer, intent(in) :: unit

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Atomic coordinates
    real(dp), intent(in) :: coords0(:,:)

    !> Total energy
    real(dp), intent(in) :: energy

    integer :: ii, ig, iat
    real(dp) :: dielEnergy
    real(dp), allocatable :: phi(:), zeta(:), area(:)

    allocate(phi(this%ddCosmo%ncav), zeta(this%ddCosmo%ncav), area(this%ddCosmo%ncav))
    ! Reset potential on the cavity, note that the potential is expected in e/Å
    call getPhi(this%chargesPerAtom, this%jmat, phi)
    ii = 0
    do iat = 1, this%ddCosmo%nat
      do ig = 1, this%ddCosmo%ngrid
        if (this%ddCosmo%ui(ig, iat) > 0.0_dp) then
          ii = ii + 1
          ! Calculate surface charge per area
          zeta(ii) = this%ddCosmo%w(ig) * this%ddCosmo%ui(ig, iat) &
              & * dot_product(this%ddCosmo%basis(:, ig), this%s(:, iat))
          ! Save surface area in Ångström²
          area(ii) = this%ddCosmo%w(ig) * Bohr__AA**2 * this%vdwRad(iat)**2
        end if
      end do
    end do

    ! Dielectric energy is the energy on the dielectric continuum
    dielEnergy = this%keps * dot_product(zeta, phi)

    write(unit, '(a)') &
        & "$info", &
        & "prog.: dftb+"

    write(unit, '(a)') &
        & "$cosmo"
    write(unit, '(2x, a:, "=", g0)') &
        & "epsilon", this%dielectricConst

    write(unit, '(a)') &
        & "$cosmo_data"
    write(unit, '(2x, a:, "=", g0)') &
        & "fepsi", this%keps, &
        & "area", sum(area)

    write(unit, '(a)') &
        & "$coord_rad", &
        & "#atom   x                  y                  z             element  radius [A]"
    do iat = 1, size(coords0, 2)
      write(unit, '(i4, 3(1x, f18.14), 2x, a4, 1x, f9.5)') &
          & iat, coords0(:, iat), trim(tolower(speciesNames(species0(iat)))), &
          & this%vdwRad(iat)*Bohr__AA
    end do

    write(unit, '(a)') &
        & "$coord_car", &
        & "!BIOSYM archive 3", &
        & "PBC=OFF", &
        & "coordinates from COSMO calculation", &
        & "!DATE"
    do iat = 1, size(coords0, 2)
      write(unit, '(a, i0, t5, 3(1x, f14.9), 1x, "COSM 1", 2(6x, a2), 1x, f6.3)') &
          & trim(speciesNames(species0(iat))), iat, coords0(:, iat)*Bohr__AA, &
          & tolower(speciesNames(species0(iat))), speciesNames(species0(iat)), 0.0_dp
    end do
    write(unit, '(a)') &
        & "end", "end"

    write(unit, '(a)') &
        & "$screening_charge"
    write(unit, '(2x, a:, "=", g0)') &
        & "cosmo", sum(zeta), &
        & "correction", 0.0_dp, &
        & "total", sum(zeta)

    write(unit, '(a)') &
        & "$cosmo_energy"
    write(unit, '(2x, a:, "=", f21.10)') &
        & "Total energy [a.u.]            ", energy, &
        & "Total energy + OC corr. [a.u.] ", energy, &
        & "Total energy corrected [a.u.]  ", energy, &
        & "Dielectric energy [a.u.]       ", dielEnergy, &
        & "Diel. energy + OC corr. [a.u.] ", dielEnergy

    write(unit, '(a)') &
        & "$segment_information", &
        & "# n             - segment number", &
        & "# atom          - atom associated with segment n", &
        & "# position      - segment coordinates [a.u.]", &
        & "# charge        - segment charge (corrected)", &
        & "# area          - segment area [A**2]", &
        & "# potential     - solute potential on segment (A length scale)", &
        & "#", &
        & "#  n   atom              position (X, Y, Z)                   charge         area        charge/area     potential", &
        & "#", &
        & "#"

    ii = 0
    do iat = 1, this%ddCosmo%nat
      do ig = 1, this%ddCosmo%ngrid
        if (this%ddCosmo%ui(ig, iat) > 0.0_dp) then
          ii = ii + 1
          write(unit, '(2i5, 7(1x, f14.9))') &
              & ii, iat, this%ddCosmo%ccav(:, ii), &
              & zeta(ii), area(ii), zeta(ii)/area(ii), phi(ii)/Bohr__AA
        end if
      end do
    end do
  end subroutine writeCosmoFile


end module dftbp_solvation_cosmo
