!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'


!> Contains workhorse routines for range-separated gamma functions.
module dftbp_dftb_rshgamma

  use dftbp_common_accuracy, only : dp, tolSameDist, MinHubDiff
  use dftbp_common_status, only : TStatus
  use dftbp_math_simplealgebra, only : cross3

  implicit none
  private

  public :: getCoulombTruncationCutoff
  public :: getCamAnalyticalGammaValue_workhorse
  public :: getHfAnalyticalGammaValue_workhorse, getdHfAnalyticalGammaValue_workhorse
  public :: getLrAnalyticalGammaValue_workhorse, getdLrAnalyticalGammaValue_workhorse
  public :: getddLrNumericalGammaValue_workhorse


contains

  !> Workhorse for calculating the general CAM gamma integral.
  subroutine getCoulombTruncationCutoff(latVecs, cutoff, errStatus, nK)

    !> Real-space lattice vectors of unit cell
    real(dp), intent(in) :: latVecs(:,:)

    !> Coulomb truncation cutoff
    real(dp), intent(out) :: cutoff

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !> Number of k-points along each direction
    integer, intent(in), optional :: nK(:)

    !! Actual number of k-points along each direction
    integer :: nK_(3)

    !! Real-space Born-von Karman supercell
    real(dp) :: bvkLatVecs(3, 3)

    !! Diagonals of triclinic cell
    real(dp) :: diag1(3), diag2(3), diag3(3)

    !! Normal vectors of three planes of the cell
    real(dp) :: n1(3), n2(3), n3(3)

    !! Point of intersection in absolute coordinates
    real(dp) :: intersect(3)

    !! Distances of sphere origin to cell planes
    real(dp) :: distances(3)

    !! Origin dummy
    real(dp) :: origin(3)

    !! Auxiliary variable
    integer :: ii

    @:ASSERT(size(latVecs, dim=1) == 3)
    @:ASSERT(size(latVecs, dim=2) == 3)

    origin(:) = 0.0_dp

    if (present(nK)) then
      @:ASSERT(size(nK) == 3)
      nK_(:) = nK
    else
      nK_(:) = 1
    end if

    do ii = 1, 3
      bvkLatVecs(:, ii) = nK_ * latVecs(:, ii)
    end do

    diag1(:) = bvkLatVecs(:, 1) + bvkLatVecs(:, 2) + bvkLatVecs(:, 3)
    diag2(:) = bvkLatVecs(:, 1) - bvkLatVecs(:, 2) + bvkLatVecs(:, 3)
    diag3(:) = bvkLatVecs(:, 2) - bvkLatVecs(:, 1) + bvkLatVecs(:, 3)

    call getDiagIntersectionPoint(bvkLatVecs(:, 2), bvkLatVecs(:, 2) + diag2, origin,&
        & origin + diag1, intersect, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! calculate normal vectors of three planes
    n1 = cross3(bvkLatVecs(:, 1), bvkLatVecs(:, 2))
    n1 = n1 / norm2(n1)

    n2 = cross3(bvkLatVecs(:, 3), bvkLatVecs(:, 1))
    n2 = n2 / norm2(n2)

    n3 = cross3(bvkLatVecs(:, 2), bvkLatVecs(:, 3))
    n3 = n3 / norm2(n3)

    distances(1) = abs(dot_product(intersect, n1))
    distances(2) = abs(dot_product(intersect, n2))
    distances(3) = abs(dot_product(intersect, n3))

    cutoff = minval(distances)

  contains

    !> Returns intersection point of two diagonals in a triclinic cell.
    subroutine getDiagIntersectionPoint(startDiag1, endDiag1, startDiag2, endDiag2, intersect1,&
        & errStatus)

      !> Start and end point of first diagonal
      real(dp), intent(in) :: startDiag1(:), endDiag1(:)

      !> Start and end point of second diagonal
      real(dp), intent(in) :: startDiag2(:), endDiag2(:)

      !> Points of intersection in absolute coordinates
      real(dp), intent(out) :: intersect1(3)

      !> Error status
      type(TStatus), intent(out) :: errStatus

      !> Points of intersection in absolute coordinates
      real(dp) :: intersect2(3)

      !! Deviation between to calculated intersection points
      real(dp) :: deviation

      !! Parameters of diagonal lines
      real(dp) :: tt, ss

      ! calculate parameters of diagonal lines
      tt = ((startDiag2(1) - startDiag1(1)) * (startDiag2(2) - endDiag2(2))&
          & - (startDiag2(2) - startDiag1(2)) * (startDiag2(1) - endDiag2(1)))&
          & / ((startDiag2(2) - endDiag2(2)) * (endDiag1(1) - startDiag1(1))&
          & - (startDiag2(1) - endDiag2(1)) * (endDiag1(2) - startDiag1(2)))
      ss = ((endDiag1(2) - startDiag1(2)) * (startDiag2(1) - startDiag1(1))&
          & - (endDiag1(1) - startDiag1(1)) * (startDiag2(2) - startDiag1(2)))&
          & / ((endDiag1(2) - startDiag1(2)) * (startDiag2(1) - endDiag2(1))&
          & - (endDiag1(1) - startDiag1(1)) * (startDiag2(2) - endDiag2(2)))

      ! check if third equation (z-component) is also satisfied
      ! (3 equations for x,y,z but only 2 variables ss and tt to determine)
      if ((tt * (endDiag1(3) - startDiag1(3)) + ss * (startDiag2(3) - endDiag2(3)))&
          & - (startDiag2(3) - startDiag1(3)) < 1.0e-10_dp) then
        if ((abs(tt) >= 0.0_dp .and. abs(tt) <= 1.0_dp)&
            & .and. (abs(ss) >= 0.0_dp .and. abs(ss) <= 1.0_dp)) then
          intersect1 = startDiag1 + tt * (endDiag1 - startDiag1)
          intersect2 = startDiag2 + ss * (endDiag2 - startDiag2)
          deviation = norm2(intersect2 - intersect1)
          if (deviation > 1.0e-10_dp) then
            @:RAISE_ERROR(errStatus, -1, "Error while calculating Coulomb truncation cutoff for&
                & given cell geometry.")
          end if
        else
          @:RAISE_ERROR(errStatus, -1, "Error while calculating Coulomb truncation cutoff for given&
              & cell geometry.")
        end if
      end if

    end subroutine getDiagIntersectionPoint

  end subroutine getCoulombTruncationCutoff


  !> Workhorse for calculating the general CAM gamma integral.
  function getCamAnalyticalGammaValue_workhorse(hubbu1, hubbu2, omega, camAlpha, camBeta, dist)&
      & result(gamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> CAM alpha and beta parameter
    real(dp), intent(in) :: camAlpha, camBeta

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting CAM gamma
    real(dp) :: gamma

    gamma =&
        & camAlpha * getHfAnalyticalGammaValue_workhorse(hubbu1, hubbu2, dist)&
        & + camBeta * getLrAnalyticalGammaValue_workhorse(hubbu1, hubbu2, omega, dist)

  end function getCamAnalyticalGammaValue_workhorse


  !> Workhorse for getHfAnalyticalGammaValue wrapper in hybrid module.
  function getHfAnalyticalGammaValue_workhorse(hubbu1, hubbu2, dist) result(gamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma
    real(dp) :: gamma

    real(dp) :: tmp, tau, tauA, tauB, distTauA, invDist

    tauA = 3.2_dp * hubbu1
    tauB = 3.2_dp * hubbu2

    if (dist < tolSameDist) then
      @:ASSERT(abs(tauA - tauB) < MinHubDiff)
    end if

    if (dist < tolSameDist) then
      ! on-site case
      tau = 0.5_dp * (tauA + tauB)
      gamma = tau * 0.3125_dp
    else
      invDist = 1.0_dp / dist
      ! off-site case, Ua == Ub
      if (abs(tauA - tauB) < MinHubDiff) then
        tauA = 0.5_dp * (tauA + tauB)
        distTauA = dist * tauA
        tmp = (distTauA**3 / 48.0_dp + 0.1875_dp * distTauA**2 + 0.6875_dp * distTauA + 1.0_dp)&
            & * exp(-tauA * dist) * invDist
        gamma = invDist - tmp
      ! off-site, Ua != Ub
      else
        gamma = invDist&
            & - getYGammaSubPart(tauA, tauB, dist, 0.0_dp)&
            & - getYGammaSubPart(tauB, tauA, dist, 0.0_dp)
      end if
    end if

  end function getHfAnalyticalGammaValue_workhorse


  !> Workhorse for getLrAnalyticalGammaValue wrapper in hybrid module.
  function getLrAnalyticalGammaValue_workhorse(hubbu1, hubbu2, omega, dist) result(gamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma
    real(dp) :: gamma

    real(dp) :: tauA, tauB, distTau, distTauA, invDist
    real(dp) :: prefac, tmp, tmp2, tau, omega2

    omega2 = omega**2

    tauA = 3.2_dp * hubbu1
    tauB = 3.2_dp * hubbu2

    if (dist < tolSameDist) then
      @:ASSERT(abs(tauA - tauB) < MinHubDiff)
    end if

    if (dist < tolSameDist) then
      ! on-site case
      tau = 0.5_dp * (tauA + tauB)
      tmp = 5.0_dp * tau**6 + 15.0_dp * tau**4 * omega2 - 5.0_dp * tau**2 * omega**4 + omega**6
      tmp = tmp * 0.0625_dp - omega * tau**5
      tmp = tmp * tau**3 / (tau**2 - omega2)**4
      gamma = tau * 0.3125_dp - tmp
    else
      invDist = 1.0_dp / dist
      ! off-site case, Ua == Ub
      if (abs(tauA - tauB) < MinHubDiff) then
        tauA = 0.5_dp * (tauA + tauB)
        distTauA = dist * tauA
        tmp2 = (distTauA**3 / 48.0_dp + 0.1875_dp * distTauA**2 +&
            & 0.6875_dp * distTauA + 1.0_dp) * exp(-tauA * dist) * invDist
        tmp = -tauA**8 / (tauA**2 - omega2)**4 * (tmp2 + exp(-tauA * dist) * &
            & (dist**2 * (3.0_dp * tauA**4 * omega**4 - 3.0_dp * tauA**6 * omega2 - &
            & tauA**2 * omega**6) + dist * (15.0_dp * tauA**3 * omega**4 - &
            & 21.0_dp * tauA**5 * omega2 - 3.0_dp * tauA * omega**6) + &
            & (15.0_dp * tauA**2 * omega**4 - 45.0_dp * tauA**4 * omega2 - &
            & 3.0_dp * omega**6)) / (48.0_dp * tauA**5))
        gamma = 1.0_dp/dist - tmp2 - (tauA**8 / (tauA**2 - omega2)**4 *&
            & exp(-omega * dist) * invDist + tmp)
      else
        ! off-site, Ua != Ub
        prefac = tauA**4 / (tauA**2 - omega2)**2
        prefac = prefac * tauB**4 / (tauB**2 - omega2)**2
        prefac = prefac * exp(-omega * dist) * invDist
        tmp = prefac&
            & - getYGammaSubPart(tauA, tauB, dist, omega)&
            & - getYGammaSubPart(tauB, tauA, dist, omega)
        tmp = invDist - tmp
        tmp = tmp&
            & - getYGammaSubPart(tauA, tauB, dist, 0.0_dp)&
            & - getYGammaSubPart(tauB, tauA, dist, 0.0_dp)
        gamma = tmp
      end if
    end if

  end function getLrAnalyticalGammaValue_workhorse


  !> Returns analytical derivative of full-range Hartree-Fock gamma.
  function getdHfAnalyticalGammaValue_workhorse(hubbu1, hubbu2, dist) result(dGamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting d gamma / d dist
    real(dp) :: dGamma

    real(dp) :: tauA, tauB
    real(dp) :: dTmp, invDist2

    tauA = 3.2_dp * hubbu1
    tauB = 3.2_dp * hubbu2

    if (dist < tolSameDist) then
      @:ASSERT(abs(tauA - tauB) < MinHubDiff)
    end if

    if (dist < tolSameDist) then
      ! on-site case
      dGamma = 0.0_dp
    else
      invDist2 = 1.0_dp / dist**2
      ! off-site case, Ua == Ub
      if (abs(tauA - tauB) < MinHubDiff) then
        tauA = 0.5_dp * (tauA + tauB)

        dTmp = &
            & (2.0_dp * dist * tauA**3 / 48.0_dp + 0.1875_dp * tauA**2 - invDist2)&
            & * exp(-tauA * dist) - (dist**2 * tauA**3 / 48.0_dp + 0.1875_dp * dist * tauA**2&
            & + 0.6875_dp * tauA + 1.0_dp / dist) * tauA * exp(-tauA * dist)

        dGamma = -invDist2 - dtmp

      ! off-site, Ua != Ub
      else
        dGamma = -invDist2&
            & - getdYGammaSubPart(tauA, tauB, dist, 0.0_dp)&
            & - getdYGammaSubPart(tauB, tauA, dist, 0.0_dp)
      end if
    end if

  end function getdHfAnalyticalGammaValue_workhorse


  !> Returns 1st analytical derivative of long-range gamma.
  function getdLrAnalyticalGammaValue_workhorse(hubbu1, hubbu2, omega, dist) result(dGamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting d gamma / d dist
    real(dp) :: dGamma

    real(dp) :: tauA, tauB, omega2, omegaDist, invDist
    real(dp) :: prefac, tmp, tmp2, dTmp, dTmp2

    tauA = 3.2_dp * hubbu1
    tauB = 3.2_dp * hubbu2

    if (dist < tolSameDist) then
      @:ASSERT(abs(tauA - tauB) < MinHubDiff)
    end if

    if (dist < tolSameDist) then
      ! on-site case
      dGamma = 0.0_dp
    else
      omega2 = omega**2
      omegaDist = omega * dist
      invDist = 1.0_dp / dist
      ! off-site case, Ua == Ub
      if (abs(tauA - tauB) < MinHubDiff) then
        tauA = 0.5_dp * (tauA + tauB)

        tmp = dist**2 * (3.0_dp * tauA**4 * omega**4 &
            &- 3.0_dp * tauA**6 * omega2 - tauA**2 * omega**6)&
            & + dist * (15.0_dp * tauA**3 * omega**4 &
            &- 21.0_dp * tauA**5 * omega2 - 3.0_dp * tauA * omega**6)&
            & + (15.0_dp * tauA**2 * omega**4 - 45.0_dp * tauA**4 * omega2 - 3.0_dp * omega**6)

        dTmp = 2.0_dp * dist * (3.0_dp * tauA**4 * omega**4&
            & - 3.0_dp * tauA**6 * omega2 - tauA**2 * omega**6)&
            & + (15.0_dp * tauA**3 * omega**4&
            & - 21.0_dp * tauA**5 * omega2 - 3.0_dp * tauA * omega**6)

        dtmp = (dtmp * exp(-tauA * dist) - tmp * tauA * exp(-tauA * dist)) / (48.0_dp * tauA**5)

        tmp2 = (dist**2 * tauA**3 / 48.0_dp + 0.1875_dp * dist * tauA**2 + 0.6875_dp * tauA&
            & + invDist) * exp(-tauA * dist)

        dTmp2 = &
            & (2.0_dp * dist * tauA**3 / 48.0_dp + 0.1875_dp * tauA**2 - invDist**2)&
            & * exp(-tauA * dist)&
            & -(dist**2 * tauA**3 / 48.0_dp + 0.1875_dp * dist * tauA**2&
            & + 0.6875_dp * tauA + invDist) * tauA * exp(-tauA * dist)

        dGamma = -invDist**2 - dtmp2&
            & + (tauA**8 / (tauA**2 - omega2)**4)&
            & * (dtmp + dtmp2 + omega * exp(-omegaDist) * invDist + exp(-omegaDist) * invDist**2)

      else
        ! off-site, Ua != Ub
        prefac = tauA**4 / (tauA**2 - omega2)**2
        prefac = prefac * tauB**4 / (tauB**2 - omega2)**2
        prefac = prefac * (-omega * exp(-omegaDist) * invDist - exp(-omegaDist) * invDist**2)
        dGamma = -invDist**2 - prefac&
            & + getdYGammaSubPart(tauA, tauB, dist, omega)&
            & + getdYGammaSubPart(tauB, tauA, dist, omega)&
            & - getdYGammaSubPart(tauA, tauB, dist, 0.0_dp)&
            & - getdYGammaSubPart(tauB, tauA, dist, 0.0_dp)
      end if
    end if

  end function getdLrAnalyticalGammaValue_workhorse


  !> Workhorse routine for getddLrNumericalGammaValue in hybrid module.
  function getddLrNumericalGammaValue_workhorse(hubbu1, hubbu2, omega, dist, delta) result(ddGamma)

    !> Hubbard U's
    real(dp), intent(in) :: hubbu1, hubbu2

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Delta for finite differences
    real(dp), intent(in) :: delta

    !> Numerical gamma derivative
    real(dp) :: ddGamma

    ddGamma = (getdLrAnalyticalGammaValue_workhorse(hubbu1, hubbu2, omega, dist + delta)&
        & - getdLrAnalyticalGammaValue_workhorse(hubbu1, hubbu2, omega, dist - delta))&
        & / (2.0_dp * delta)

  end function getddLrNumericalGammaValue_workhorse


  !> Returns the subexpression for the evaluation of the off-site Y-Gamma-integral.
  pure function getYGammaSubPart(tauA, tauB, dist, omega) result(yGamma)

    !> Decay constant site A
    real(dp), intent(in) :: tauA

    !> Decay constant site B
    real(dp), intent(in) :: tauB

    !> Separation of the sites A and B
    real(dp), intent(in) :: dist

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> Resulting off-site Y-Gamma-integral
    real(dp) :: yGamma

    real(dp) :: prefac, tmp

    tmp = (tauA - omega)
    tmp = tmp * (tauA + omega)
    prefac = tauA**2 / tmp
    tmp = (tauB**6 - 3.0_dp * tauA**2 * tauB**4 + 2.0_dp * omega**2 * tauB**4) / dist
    tmp = tmp * prefac * prefac / (tauA**2 - tauB**2)**3
    tmp = tauA * tauB**4 * 0.5_dp * prefac / (tauB**2 - tauA**2)**2 - tmp
    yGamma = tmp * exp(-tauA * dist)

  end function getYGammaSubPart


  !> Returns the derivative of the subexpression for the evaluation of the off-site
  !! Y-Gamma-integral. Note that tauA /= tauB.
  pure function getdYGammaSubPart(tauA, tauB, dist, omega) result(dYGammaSubPart)

    !> Decay constant site A
    real(dp), intent(in) :: tauA

    !> Decay constant site B
    real(dp), intent(in) :: tauB

    !> Separation of the sites A and B
    real(dp), intent(in) :: dist

    !> Range-separation parameter
    real(dp), intent(in) :: omega

    !> Resulting derivative of the subexpression
    real(dp) :: dYGammaSubPart

    !! Auxiliary variables
    real(dp) :: prefac, tmp, tmp2, dtmp, invDist

    invDist = 1.0_dp / dist
    tmp = tauA**2 - omega**2
    prefac = tauA**2 / tmp
    tmp = prefac * prefac / (tauA**2 - tauB**2)**3
    dtmp = tmp * (tauB**6 - 3.0_dp * tauA**2 * tauB**4 + 2.0_dp * omega**2 * tauB**4) * invDist**2
    tmp = tmp * (tauB**6 - 3.0_dp * tauA**2 * tauB**4 + 2.0_dp * omega**2 * tauB**4) * invDist
    tmp2 = tauA * tauB**4 * 0.5_dp * prefac / (tauB**2 - tauA**2)**2 - tmp

    dYGammaSubPart = (dtmp - tmp2 * tauA) * exp(-tauA * dist)

  end function getdYGammaSubPart

end module dftbp_dftb_rshgamma
