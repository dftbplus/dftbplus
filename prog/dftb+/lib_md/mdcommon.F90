!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Common routines for MD calculations
module mdcommon
  use assert
  use accuracy
  use constants
  use ranlux
  implicit none
  private

  public :: OMDCommon, init, restFrame, evalKT, rescaleToKT
  public :: evalKE, BoxMueller, MaxwellBoltzmann

  !!* Contains necessary data for the MD framework
  type OMDCommon
!    private
    real(dp) :: Nf               ! Nr. of degrees of freedom
    logical :: tStationary       ! Always transform to rest frame
  end type OMDCommon

  interface init
    module procedure MDCommon_init
  end interface

  interface restFrame
    module procedure MDCommon_restFrame
  end interface

  interface evalKT
    module procedure MDCommon_evalKT
  end interface

  interface rescaleToKT
    module procedure MDCommon_rescaleToKT
  end interface


contains

  !!* Creates MD Framework.
  !!* @param sf MD Framework instance.
  !!* @param nMoverAtoms Number of moving atoms in the system
  !!* @param nAllAtom   Total number of real atoms in the system
  !!* @param tStationary  If system should be transformed to rest frame.
  subroutine MDCommon_init(sf, nMovedAtom, nAllAtom, tStationary)
    type(OMDCommon), intent(out) :: sf
    integer, intent(in) :: nMovedAtom
    integer, intent(in) :: nAllAtom
    logical             :: tStationary

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


  !!* Shift velocities so the total velocity is 0
  !!* @param sf MD Framework instance.
  !!* @param velocity particle velocities
  !!* @param mass particle masses
  subroutine  MDCommon_restFrame(sf, velocity, mass)
    type(OMDCommon), intent(in) :: sf
    real(dp), intent(inout) :: velocity(:,:)
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



  !!* Calculate the kinetic temperature of an integrator.
  !!* @param sf MD Framework instance.
  !!* @parameter kT resulting thermal energy
  !!* @param velocity particle velocities
  !!* @param mass particle masses
  subroutine MDCommon_evalKT(sf, kT, velocity, mass)
    type(OMDCommon), intent(in) :: sf
    real(dp),intent(out) :: kT
    real(dp), intent(in) :: velocity(:,:)
    real(dp), intent(in) :: mass(:)

    real(dp) :: kinE

    @:ASSERT(size(velocity,dim=2)==size(mass))
    @:ASSERT(size(velocity,dim=1)==3)

    call evalKE(kinE, velocity, mass)

    kT = 2.0_dp * kinE / sf%Nf

  end subroutine MDCommon_evalKT



  !!* Rescales the velocities of a system to match the target thermal energy.
  !!* @param sf MD Framework instance.
  !!* @param velocity particle velocities
  !!* @param mass particle masses
  !!* @param kTtarget intended kinetic energy
  subroutine MDCommon_rescaleTokT(sf, velocity, mass, kTtarget)
    type(OMDCommon), intent(in) :: sf
    real(dp), intent(inout) :: velocity(:,:)
    real(dp), intent(in)    :: mass(:)
    real(dp), intent(in)    :: kTtarget

    real(dp) :: currentkT

    @:ASSERT(size(velocity,dim=2) == size(mass))
    @:ASSERT(size(velocity,dim=1) == 3)

    call evalkT(sf, currentkT, velocity, mass)
    velocity(:,:) = velocity(:,:) * sqrt(kTtarget/currentkT)

  end subroutine MDCommon_rescaleTokT



  !!* Calculate the kinetic energy of the atoms
  !!* @parameter kinE resulting energy
  !!* @param velocity particle velocities
  !!* @param mass particle masses
  subroutine evalKE(kinE, velocity, mass)
    real(dp),intent(out) :: kinE
    real(dp), intent(in) :: velocity(:,:)
    real(dp), intent(in) :: mass(:)

    @:ASSERT(size(velocity,dim=2) == size(mass))
    @:ASSERT(size(velocity,dim=1) == 3)

    kinE = 0.5_dp * sum(spread(mass(:),1,3) * velocity(:,:)**2)

  end subroutine evalKE



  !!* Converts a uniform distribution into a Gaussian distribution.
  !!* @parameter eta1 number with Gaussian distribution
  !!* @parameter eta2 number with Gaussian distribution
  !!* @parameter u1 number with uniform distribution
  !!* @parameter u2 number from uniform distribution
  subroutine BoxMueller(eta1,eta2,u1,u2)
    real(dp), intent(out) :: eta1
    real(dp), intent(out) :: eta2
    real(dp), intent(in)  :: u1
    real(dp), intent(in)  :: u2

    real(dp) :: theta, a

    a = sqrt( -2.0_dp*log(u1) )
    theta = 2.0_dp*pi*u2

    eta1 = a * cos(theta)
    eta2 = a * sin(theta)

  end subroutine BoxMueller



  !!* Draws an atom velocity from a Maxwell-Boltzmann distribution.
  !!* @param velocity resulting velocity
  !!* @param mass atomic mass in a.u.
  !!* @param kT system thermal energy
  !!* @param pRanlux pointer to a random number generator
  subroutine MaxwellBoltzmann(velocity,mass,kT,pRanlux)
    real(dp), intent(out)  :: velocity(:)
    real(dp), intent(in)   :: mass
    real(dp), intent(in)   :: kT
    type(ORanlux), pointer :: pRanlux

    real(dp) :: ranvals(7)
    real(dp) :: junk

    @:ASSERT(size(velocity)==3)

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

end module mdcommon
