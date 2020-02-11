!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> xTB parameter types
module dftbp_xtbparam
  use dftbp_accuracy, only : dp
  use dftbp_constants, only : maxL
  implicit none
  private

  public :: xtbParam, xtbBasis, xtbGlobalParameter

  !> xTB-basisset definition
  type xtbBasis

    !> principal quantum number
    integer :: n

    !> angular momentum of shell
    integer :: l

    !> shell polynomial
    real(dp) :: poly

    !> chemical hardness shell parameter
    real(dp) :: lEta

    !> atomic level
    real(dp) :: h

    !> coordination number dependence of the level
    real(dp) :: kcn

    !> Slater function exponent
    real(dp) :: zeta

    !> Number of Gaussian functions for expansion
    integer :: nGauss

    !> Is valence or polarization shell.
    integer :: valence

    !> Reference occupation
    real(dp) :: referenceN0

  end type


  !> xTB-parameters
  type :: xtbParam

    !> Number of shells
    integer :: nSh

    !> Metal
    integer :: metal

    !> Electronegativity
    real(dp) :: en

    !> Atomic radius
    real(dp) :: radius

    !> Chemical hardness
    real(dp) :: eta

    !> Third order Hubbard derivative
    real(dp) :: gam

    !> Repulsion parameter
    real(dp) :: alpha

    !> Repulsion parameter
    real(dp) :: zeff

    !> Halogen bond strength
    real(dp) :: xbond = 0.0_dp

    !> Anisotropy radius
    real(dp) :: anisoRad = 0.0_dp

    !> Dipole kernel
    real(dp) :: dpolc = 0.0_dp

    !> Quadrupole kernel
    real(dp) :: qpolc = 0.0_dp

    !> Valence coordination number
    real(dp) :: valenceCN = 0.0_dp

    !> Basis functions
    type(xtbBasis), allocatable :: basis(:)

  end type xtbParam


  !> Global parameters used in xTB methods
  type :: xtbGlobalParameter

    !> Scaling parameter for EN differences
    real(dp) :: kEnScale

    !> Repulsive polynomial parameter
    real(dp) :: krep

    !> Repulsive polynomial parameter
    real(dp) :: rrep

    !> Halogen bonding parameter
    real(dp) :: xbondRad

    !> Halogen bonding parameter
    real(dp) :: xbondExp

    !> Shell pair scaling parameters for the Hamiltonian
    real(dp) :: kH0Scale(0:maxL, 0:maxL)

    !> Shell pair scaling parameter for polarisation functions
    real(dp) :: kPol

  end type xtbGlobalParameter


contains


  !> Serialize xTB-parameters to TOML
  subroutine xtbParamTOML(unit, param, gfn, root)

    !> Unit for IO
    integer, intent(in) :: unit

    !> Instance
    type(xtbParam), intent(in) :: param

    !> GFN level
    integer, intent(in), optional :: gfn

    !> TOML table name
    character(len=*), intent(in), optional :: root

    !> TOML key-value format
    character(len=*), parameter :: tomlf = '(a,1x,"=",1x,g0)'

    integer :: iSh

    if (present(gfn)) then
      write(unit, tomlf) "GFN", gfn
    end if
  #:for key in ("nSh", "metal", "en", "radius", "eta", "gam", "alpha", "zeff", "xbond", "anisoRad", "dpolc", "qpolc", "valenceCN")
    write(unit, tomlf) "${key}$", param%${key}$
  #:endfor

    do iSh = 1, param%nSh
      if (present(root)) then
        write(unit, '("[[",a,".",a,"]]")') root, "basis"
      else
        write(unit, '("[[",a,"]]")') root, "basis"
      end if
    #:for key in ("n", "l", "poly", "lEta", "h", "kcn", "zeta", "nGauss", "valence")
      write(unit, tomlf) "${key}$", param%basis(iSh)%${key}$
    #:endfor
    end do

  end subroutine xtbParamTOML

end module dftbp_xtbparam
