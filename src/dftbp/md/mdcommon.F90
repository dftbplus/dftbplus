!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Common routines for MD calculations
module dftbp_md_mdcommon
  use dftbp_common_accuracy, only : dp, rsp
  use dftbp_common_constants, only : pi
  use dftbp_geometry_projectvectors, only : getCentreOfMass, getPrincipleAxes
  use dftbp_math_lapackroutines, only : gesv
  use dftbp_math_ranlux, only : getRandom, TRanlux
  use dftbp_math_simplealgebra, only : cross3
  implicit none

  private
  public :: TMDCommon, TMDCommon_init, restFrame, evalKT, rescaleToKT
  public :: evalKE, BoxMueller, MaxwellBoltzmann, TMDOutput


  !> Contains necessary data for the MD framework
  type TMDCommon

    !> Nr. of degrees of freedom
    real(dp) :: Nf

    !> Number of moving atoms
    integer :: nMovedAtom

    !> Total number of atoms
    integer :: nAllAtom

    !> Should atomic translation be transformed to rest frame?
    logical :: isTranslationRemoved

    !> Should atomic rotation be transformed to rest frame?
    logical :: isRotationRemoved

    !> Number of degrees of freedom for rotation
    integer :: nRotationDegrees = 0

  end type TMDCommon


  !> Output variables accumulated during MD
  type TMDOutput

    !> Are eigenvalues printed out at every write time?
    logical :: bandStructure = .false.

    !> Are 1st derivative data accumulated?
    logical :: printForces = .false.

    !> Are charge-related data accumulated?
    logical :: printCharges = .false.

    !> Are atom resolved energies printed?
    logical :: printAtomEnergies = .false.

  end type TMDOutput

  !> transform to co-moving frame if needed
  interface restFrame
    module procedure MDCommon_restFrame
  end interface


  !> evaluate thermal energy
  interface evalKT
    module procedure MDCommon_evalKT
  end interface


  !> rescale velocities to temperature
  interface rescaleToKT
    module procedure MDCommon_rescaleToKT
  end interface

contains


  !> Creates MD Framework.
  subroutine TMDCommon_init(sf, nMovedAtom, nAllAtom, isTranslationRemoved, mass, coords,&
      & isRotationRemoved)

    !> MD Framework instance.
    type(TMDCommon), intent(out) :: sf

    !> Number of moving atoms in the system
    integer, intent(in) :: nMovedAtom

    !> Total number of real atoms in the system
    integer, intent(in) :: nAllAtom

    !> If system should be transformed to translational rest frame
    logical, intent(in) :: isTranslationRemoved

    !> Particle masses
    real(dp), intent(in) :: mass(:)

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> If system should be transformed to rotational rest frame
    logical, intent(in) :: isRotationRemoved

    integer :: ii
    real(dp) :: cm(3), axes(3,3), moments(3)

    @:ASSERT(nMovedAtom <= nAllAtom)
    sf%nAllAtom = nAllAtom
    sf%nMovedAtom = nMovedAtom
    sf%isTranslationRemoved = isTranslationRemoved
    sf%isRotationRemoved = isRotationRemoved

    if (sf%isRotationRemoved) then
      cm(:) = getCentreOfMass(coords, mass)
      call getPrincipleAxes(axes, moments, coords, mass, cm)
      sf%nRotationDegrees = 0
      do ii = 1, 3
        if (moments(ii) < epsilon(0.0_rsp)) then
          ! zero moment of inertia - linear molecule, and this direction is along its axis
          cycle
        end if
        sf%nRotationDegrees = sf%nRotationDegrees + 1
      end do
    end if

    call updateNf(sf)

  end subroutine TMDCommon_init


  !> Shift velocities so the total velocity is 0
  subroutine  MDCommon_restFrame(sf, velocity, mass, coords)

    !> MD Framework instance.
    type(TMDCommon), intent(inout) :: sf

    !> Particle velocities
    real(dp), intent(inout) :: velocity(:,:)

    !> Particle masses
    real(dp), intent(in) :: mass(:)

    !> Particle coordinates
    real(dp), intent(in) :: coords(:,:)

    integer :: ii, iAt
    real(dp) :: mv(3), cm(3), axes(3,3), moments(3), iTmp(3), rTmp(3), vTmp(3), rMag
    real(dp) :: L(3), omega(3)

    if (sf%isTranslationRemoved) then
      ! calculate total momentum of the system
      mv(:) = sum(spread(mass(:), 1, 3) * velocity(:,:), dim=2)
      ! get the total velocity of the system
      mv(:) = mv(:) / sum(mass)
      ! and shift so that it is 0
      velocity(:,:) = velocity(:,:) - spread(mv(:), 2, size(mass))
    end if

    if (sf%isRotationRemoved) then

      cm(:) = getCentreOfMass(coords, mass)
      call getPrincipleAxes(axes, moments, coords, mass, cm)

      sf%nRotationDegrees = 0
      do ii = 1, 3
        if (moments(ii) < epsilon(0.0_rsp)) then
          ! zero moment of inertia - linear molecule, and this direction is along its axis
          cycle
        end if
        sf%nRotationDegrees = sf%nRotationDegrees + 1
      end do

      L(:) = 0.0_dp
      do iAt = 1, sf%nAllAtom
        rTmp(:) = coords(:, iAt) - cm
        L(:) = L + cross3(rTmp, velocity(:,iAt) * mass(iAt))
      end do

      omega(:) = 0.0_dp
      do ii = 1, 3
        if (moments(ii) < epsilon(0.0_rsp)) then
          ! zero moment of inertia - linear molecule, and this direction is along its axis
          cycle
        end if
        omega(:) = omega + axes(:,ii) * dot_product(L,axes(:,ii)) / moments(ii)
      end do

      do iAt = 1, sf%nAllAtom
        rTmp(:) = coords(:, iAt) - cm
        vTmp(:) = cross3(rTmp, omega)
        velocity(:,iAt) = velocity(:,iAt) + vTmp
      end do

    end if

    call updateNf(sf)

  end subroutine MDCommon_restFrame


  !> Calculate the kinetic temperature of an integrator.
  subroutine MDCommon_evalKT(sf, kT, velocity, mass)

    !> MD Framework instance.
    type(TMDCommon), intent(in) :: sf

    !> resulting thermal energy
    real(dp),intent(out) :: kT

    !> particle velocities
    real(dp), intent(in) :: velocity(:,:)

    !> particle masses
    real(dp), intent(in) :: mass(:)

    real(dp) :: kinE

    @:ASSERT(size(velocity,dim=2)==size(mass))
    @:ASSERT(size(velocity,dim=1)==3)

    call evalKE(kinE, velocity, mass)

    kT = 2.0_dp * kinE / sf%Nf

  end subroutine MDCommon_evalKT


  !> Rescales the velocities of a system to match the target thermal energy.
  subroutine MDCommon_rescaleTokT(sf, velocity, mass, kTtarget)

    !> MD Framework instance.
    type(TMDCommon), intent(in) :: sf

    !> particle velocities
    real(dp), intent(inout) :: velocity(:,:)

    !> particle masses
    real(dp), intent(in) :: mass(:)

    !> intended kinetic energy
    real(dp), intent(in) :: kTtarget

    real(dp) :: currentkT

    @:ASSERT(size(velocity,dim=2) == size(mass))
    @:ASSERT(size(velocity,dim=1) == 3)

    call evalkT(sf, currentkT, velocity, mass)
    velocity(:,:) = velocity(:,:) * sqrt(kTtarget/currentkT)

  end subroutine MDCommon_rescaleTokT


  !> Calculate the kinetic energy of the atoms
  subroutine evalKE(kinE, velocity, mass)

    !> resulting energy
    real(dp),intent(out) :: kinE

    !> particle velocities
    real(dp), intent(in) :: velocity(:,:)

    !> particle masses
    real(dp), intent(in) :: mass(:)

    @:ASSERT(size(velocity,dim=2) == size(mass))
    @:ASSERT(size(velocity,dim=1) == 3)

    kinE = 0.5_dp * sum(spread(mass(:),1,3) * velocity(:,:)**2)

  end subroutine evalKE


  !> Converts a uniform distribution into a Gaussian distribution.
  subroutine BoxMueller(eta1,eta2,u1,u2)

    !> number with Gaussian distribution
    real(dp), intent(out) :: eta1

    !> number with Gaussian distribution
    real(dp), intent(out) :: eta2

    !> number with uniform distribution
    real(dp), intent(in) :: u1

    !> number from uniform distribution
    real(dp), intent(in) :: u2

    real(dp) :: theta, a

    a = sqrt( -2.0_dp*log(u1) )
    theta = 2.0_dp*pi*u2

    eta1 = a * cos(theta)
    eta2 = a * sin(theta)

  end subroutine BoxMueller


  !> Draws an atom velocity from a Maxwell-Boltzmann distribution.
  subroutine MaxwellBoltzmann(velocity,mass,kT,pRanlux)

    !> resulting velocity
    real(dp), intent(out) :: velocity(3)

    !> atomic mass in a.u.
    real(dp), intent(in) :: mass

    !> system thermal energy
    real(dp), intent(in) :: kT

    !> Random number generator
    type(TRanlux), intent(inout) :: pRanlux

    real(dp) :: ranvals(7)
    real(dp) :: junk

    ! use the uniform distribution to get a normal (Gaussian) distribution
    ! and then scale to get a Maxwell-Boltzmann

    call getRandom(pRanlux, ranvals)

    call BoxMueller(velocity(1),velocity(2),ranvals(1),ranvals(2))
    call BoxMueller(velocity(3),junk,ranvals(3),ranvals(4))

    if (ranvals(5) < 0.5_dp) then
      velocity(1) = - velocity(1)
    end if
    if (ranvals(6) < 0.5_dp) then
      velocity(2) = - velocity(2)
    end if
    if (ranvals(7) < 0.5_dp) then
      velocity(3) = - velocity(3)
    end if

    velocity(:) = velocity(:) * sqrt(kT/mass)

  end subroutine MaxwellBoltzmann


  !> Updates the number of degrees of freedom
  subroutine updateNf(sf)

    !> MD Framework instance
    type(TMDCommon), intent(inout) :: sf

    integer :: nDegrees

    nDegrees = 3 * sf%nMovedAtom
    if (sf%isTranslationRemoved) then
      ! Translational motion is removed, total of 3n - 3 degrees of freedom
      nDegrees = nDegrees - 3
    end if
    if (sf%isRotationRemoved) then
      ! Rotational motion is removed, either 3 or 2 depending if system is linear
      nDegrees = nDegrees - sf%nRotationDegrees
    end if

    sf%Nf = real(nDegrees, dp)

  end subroutine updateNf

end module dftbp_md_mdcommon
