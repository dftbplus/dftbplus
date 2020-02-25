!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Common routines for MD calculations
module dftbp_mdcommon
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_constants
  use dftbp_ranlux
  implicit none
  private

  public :: TMDCommon, init, restFrame, evalKT, rescaleToKT
  public :: evalKE, BoxMueller, MaxwellBoltzmann


  !> Contains necessary data for the MD framework
  type TMDCommon

    !> Nr. of degrees of freedom
    real(dp) :: Nf

    !> Should transform to rest frame?
    logical :: tStationary
  end type TMDCommon


  !> initialise thermostat
  interface init
    module procedure MDCommon_init
  end interface


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
  subroutine MDCommon_init(sf, nMovedAtom, nAllAtom, tStationary)

    !> MD Framework instance.
    type(TMDCommon), intent(out) :: sf

    !> Number of moving atoms in the system
    integer, intent(in) :: nMovedAtom

    !> Total number of real atoms in the system
    integer, intent(in) :: nAllAtom

    !> If system should be transformed to rest frame.
    logical :: tStationary

    @:ASSERT(nMovedAtom <= nAllAtom)

    if (nMovedAtom /= nAllAtom .or. .not. tStationary) then
      ! there are fixed atoms present, all  moving atoms have 3 degrees of fd.
      sf%Nf = real(3 * nMovedAtom, dp)
    else
      ! translational motion is removed, total of 3n - 3 degrees of freedom
      sf%Nf = real(3 * (nMovedAtom - 1), dp)
    end if
    sf%tStationary = tStationary

  end subroutine MDCommon_init


  !> Shift velocities so the total velocity is 0
  subroutine  MDCommon_restFrame(sf, velocity, mass)

    !> MD Framework instance.
    type(TMDCommon), intent(in) :: sf

    !> Particle velocities
    real(dp), intent(inout) :: velocity(:,:)

    !> Particle masses
    real(dp), intent(in) :: mass(:)

    real(dp) :: mv(3)

    if (sf%tStationary) then
      ! calculate total momentum of the system
      mv(:) = sum(spread(mass(:), 1, 3) * velocity(:,:), dim=2)

      ! get the total velocity of the system
      mv(:) = mv(:) / sum(mass)

      ! and shift so that it is 0
      velocity(:,:) = velocity(:,:) - spread(mv(:), 2, size(mass))
    end if

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

end module dftbp_mdcommon
