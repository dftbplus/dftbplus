!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Data types for linear response routines
module dftbp_linresptypes
  use dftbp_accuracy
  implicit none

  public

  !> Data type for linear response internal settings
  type :: TLinResp

    !> Number of excitations to be calculated
    integer :: nExc

    !> state of interest (< 0 find brightest, 0 calculate all nexc states, > 0 that specific state)
    integer :: nStat

    !> is an energy window specified
    logical :: tEnergyWindow

    !> energy window for transitions above energy of nexc-th single particle transtion
    real(dp) :: energyWindow

    !> is an oscillator window specified
    logical :: tOscillatorWindow

    !> the window for transitions not included in nexc and energy window (if used)
    real(dp) :: oscillatorWindow

    !> onsite corrections (if in use)
    real(dp), allocatable :: onSiteMatrixElements(:,:,:,:)

    !> should transition charges be cached or evaluated on the fly?
    logical :: tCacheCharges

    !> Number of atoms
    integer :: nAtom

    !> number of electrons in system
    real(dp) :: nEl

    !> symmetry required singlet ('S'), triplet ("T") or both ("B")
    character :: symmetry

    !> ground state spin constants for each species
    real(dp), allocatable :: spinW(:)

    !> ground state Hubbard U values for each species
    real(dp), allocatable :: HubbardU(:)

    !> file for X+Y data
    integer :: fdXplusY = -1

    !> file unit if the coefficients for the excited states should be written to disc
    integer :: fdCoeffs = -1

    !> Add the ground state to the excited state transition density matrix when determining the
    !> natural orbitals
    logical :: tGrndState = .true.

    !> file unit for excited Mulliken populations?
    integer :: fdMulliken = -1

    !> File unit for single particle (KS) transitions if required
    integer :: fdTrans = -1

    !> File unit for single particle (KS) transition charges if required
    integer :: fdTransQ = -1

    !> File unit for single particle transition dipole strengths
    integer :: fdSPTrans = -1

    !> file handle for excitation energies
    integer :: fdExc = -1

    !> File unit for transition dipole data
    integer :: fdTradip = -1

    !> For calculations where the geometry changes, previous vectors for restarting the iterative
    !> eigensolver. Note: in the case of ARPACK this is the residual not the eigenvectors
    real(dp), allocatable :: oldEigenVectors(:,:)

    !> Should the density matrix be stored to disc?
    logical :: tWriteDensityMatrix

    ! ARPACK related

    !> write state of Arnoldi solver to disc
    logical :: tArnoldi

    !> file unit for Arnoldi solver file unit for tests on output of Arnoldi solver
    integer :: fdArnoldi = -1

    !> file unit for Arnoldi solver tests, if this is < 1 no tests are performed
    integer :: fdArnoldiDiagnosis = -1

    !> Is the data structure initialised?
    logical :: tInit = .false.

  end type TLinResp
  
end module dftbp_linresptypes
