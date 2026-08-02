!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines related to finite electron temperature, including Fermi, Gaussian and
!! Methfessel-Paxton broadening functions.
!! To do: Add other methods, including possibly Pederson and Jackson method PRB 43, 7312
!! (1991). Also fix exact occupation for electron numers, using interpolation instead of bisection.
module dftbp_dftb_etemp
  use dftbp_common_accuracy, only : dp, elecTol, elecTolMax, mExpArg, rsp
  use dftbp_common_constants, only : pi
  use dftbp_io_message, only : error
  use dftbp_math_binarysearch, only : search_asc_real_gt, search_asc_real_geq
  use dftbp_math_errorfunction, only : erfcwrap
  use dftbp_math_factorial, only : fact
  use dftbp_math_hermite, only : hx
  use dftbp_math_sorting, only : index_heap_sort
  use dftbp_math_summation, only : kahan
  implicit none

  private
  public :: Efilling, electronFill, fillingTypes

  type :: TFillingTypesEnum

    !> Definition of a type of broadening function - Fermi-Dirac in this case
    integer :: Fermi = 0

    !> Definition of a type of broadening function - Gaussian in this case
    integer :: Gaussian = 1

    !> Definition of a type of broadening function - Methfessel-Paxton, for higher orders use
    !! Methfessel + n as a value
    integer :: Methfessel = 2

  end type TFillingTypesEnum

  !> Enumerated filling types.
  type(TFillingTypesEnum), parameter :: fillingTypes = TFillingTypesEnum()

  !> Twice the machine precision
  real(dp), parameter :: epsilon2 = 2.0_dp * epsilon(1.0_dp)

  !> Fermi-broadening function decays to effective zero (double precision)
  real(dp), parameter :: fermiWidth = 36.0_dp

  !> Gauss-broadening function decays to effective zero (double precision)
  real(dp), parameter :: gaussianWidth = 6.0_dp


contains


  !> Driver to calculate electron filling, the band-structure energy at T and extrapolated to T=0K,
  !! and the entropy of the electron energy for the Mermin free energy, returning band energy and
  !! entropy for each channel but common Fermi level
  !!
  !! Note: use slices of eigenvalues and other input arrays if you want a different Fermi energy and
  !! filling for different k-points and/or spins.
  !!
  !! Note: If no electrons are present, the Fermi energy is set to zero per default.
  subroutine Efilling(Ebs, Ef, TS, E0, filling, eigenvals, nElectrons, kT, kWeight, distrib)

    !> Band structure energy at T
    real(dp), intent(out) :: Ebs(:)

    !> Fermi energy for given distribution
    real(dp), intent(out) :: Ef

    !> Entropy
    real(dp), intent(out) :: TS(:)

    !> Band structure energy extrapolated to T=0K
    real(dp), intent(out) :: E0(:)

    !> Electron occupancies
    real(dp), intent(out) :: filling(:,:,:)

    !> The eigenvalues of the levels, 1st index is energy 2nd index is k-point and 3nd index is spin
    real(dp), intent(in) :: eigenvals(:,:,:)

    !> Number of electrons
    real(dp), intent(in) :: nElectrons

    !> Thermal energy in atomic units
    real(dp), intent(in) :: kT

    !> The k-point weightings
    real(dp), intent(in) :: kWeight(:)

    !> Choice of distribution functions, currently Fermi, Gaussian and Methfessle-Paxton supported.
    !! The flags is defined symbolically, so (Methfessel + 2) gives the 2nd order M-P scheme
    integer, intent(in) :: distrib

    real(dp) :: upperEf, lowerEf
    real(dp) :: nElec
    real(dp) :: nElecMax, nElecMin, maxEig, minEig

    @:ASSERT(all(shape(filling) == shape(eigenvals)))
    @:ASSERT(size(eigenvals, dim=3) == size(Ebs))
    @:ASSERT(size(eigenvals, dim=3) == size(TS))
    @:ASSERT(size(eigenvals, dim=3) == size(E0))
    @:ASSERT(nElectrons >= 0.0_dp)

    ! Not a tight enough bound ? :
    @:ASSERT(ceiling(nElectrons) <= 2 * size(eigenvals, dim=1) * size(eigenvals, dim=3))

    @:ASSERT(kT > 0.0_dp)
    @:ASSERT(size(kWeight) > 0)
    @:ASSERT(all(kWeight >= 0.0_dp))

    @:ASSERT(distrib >= fillingTypes%Fermi)

    ! If no electrons there, we are ready
    if (nElectrons < epsilon(1.0_dp)) then
      filling(:,:,:) = 0.0_dp
      Ebs(:) = 0.0_dp
      ! place the Fermi energy well below the lowest eigenvalue
      Ef = minval(eigenvals) - 1000.0_dp * (kT + epsilon(1.0_rsp))
      TS(:) = 0.0_dp
      E0(:) = 0.0_dp
      return
    end if

    if (size(filling,dim=1)*size(filling,dim=3) <= nElectrons) then
      ! place the Fermi energy well above the highest eigenvalue, as nOrbs * spin <= nElec
      Ef = maxval(eigenvals) + 1000.0_dp * (kT + epsilon(1.0_rsp))
      call electronFill(Ebs, filling, TS, E0, Ef, eigenvals, kT, distrib, kWeight)
      return
    end if

    ! For integer number of electrons, try middle gap for Ef
    if (abs(nElectrons - nint(nElectrons)) <= elecTol) then
      Ef = middleGap(eigenvals, kWeight, nElectrons)
      nElec = electronCount(Ef, eigenvals, kT, distrib, kWeight)
      if (abs(nElectrons - nElec) <= elecTolMax) then
        call electronFill(Ebs, filling, TS, E0, Ef, eigenvals, kT, distrib, kWeight)
        return
      end if
    end if

    ! find maximum and minimum possible value of Fermi Energy
    minEig = minval(eigenvals(1,:,:))
    maxEig = maxval(eigenvals(size(eigenvals, dim=1),:,:))
    ! Fermi level hopefully between highest and lowest eigenvalue
    upperEf = maxEig + 0.01_dp
    lowerEf = minEig - 0.01_dp

    ! but just to be on the safe side if the temperature is BIG compared to the bandwidth, or if the
    ! system has a fully filled band structure:
    nElecMax = electronCount(upperEf, eigenvals, kT, distrib, kWeight)
    nElecMin = electronCount(lowerEf, eigenvals, kT, distrib, kWeight)
    do while (nElecMin > nElectrons)
      lowerEf = 2.0_dp * (lowerEf - upperEf) + lowerEf
      nElecMin = electronCount(lowerEf, eigenvals, kT, distrib, kWeight)
    end do
    do while (nElecMax < nElectrons)
      upperEf = 2.0_dp * (upperEf - lowerEf) + lowerEf
      nElecMax = electronCount(upperEf, eigenvals, kT, distrib, kWeight)
    end do

    Ef = 0.5_dp * (upperEf + lowerEf)
    nElec = electronCount(Ef, eigenvals, kT, distrib, kWeight)

    ! Bisection as long as nr. electrons is not accurate enough or next change in the Fermi level
    ! would go below precision
    do while (abs(nElectrons - nElec) > elecTol &
        & .and. abs(upperEf - lowerEf) >= max(abs(Ef) * epsilon2, epsilon2))
      if ((nElecMax >= nElecMin) .eqv. (nElectrons >= nElec)) then
        lowerEf = Ef
        nElecMin = nElec
      else
        upperEf = Ef
        nElecMax = nElec
      end if
      Ef = 0.5_dp * (upperEf + lowerEf)
      nElec = electronCount(Ef, eigenvals, kT, distrib, kWeight)
    end do

    ! If number of electrons deviates from theoretical value too much: stop
    if (abs(nElectrons - nElec) > elecTolMax) then
      call error("Fermi level search did not converge.")
    end if

    nElec = electronCount(Ef, eigenvals, kT, distrib, kWeight)

    call electronFill(Ebs, filling, TS, E0, Ef, eigenvals, kT, distrib, kWeight)

    ! Re-scale to give exact number of electrons, this is a temporay hack
    if (nElectrons > epsilon2 .and. abs(nElectrons - nElec) > epsilon2) then
      filling(:,:,:) = filling * nElectrons / nElec
    end if

  end subroutine Efilling


  !> Calculates the number of electrons for a given Fermi energy and distribution function
  !! Totals are made using using Kahan–Babušk summation
  function electronCount(Ef,eigenvals,kT,distrib,kWeight)

    !> Electrons for this Fermi energy
    real(dp) :: electronCount

    !> Fermi energy for given distribution
    real(dp), intent(in) :: Ef

    !> The eigenvalues of the levels, 1st index is energy 2nd index is k-point and 3nd index is spin
    real(dp), intent(in) :: eigenvals(:,:,:)

    !> Thermal energy in atomic units
    real(dp), intent(in) :: kT

    !> Choice of distribution functions, currently Fermi, Gaussian and Methfessle-Paxton
    !! supported. The flags is defined symbolically, so (Methfessel + 2) gives the 2nd order M-P
    !! scheme
    integer, intent(in) :: distrib

    !> The k-point weightings
    real(dp), intent(in) :: kWeight(:)

    real(dp), allocatable :: A(:), hermites(:)
    integer ii, jj , kk, MPorder, iSpin
    integer :: nLastFully, nLastPartially
    real(dp) :: beta, occ, comp, x

    beta = 1.0_dp / kT
    electronCount = 0.0_dp
    comp = 0.0_dp

    if (distrib /= fillingTypes%Fermi) then

      MPorder = distrib - fillingTypes%Methfessel
      allocate(A(0:MPorder))
      allocate(hermites(0:2*MPorder))
      call Aweights(A,MPorder)
      do iSpin = 1, size(eigenvals,dim=3)
        do ii = 1, size(kWeight)

          call search_asc_real_geq(nLastFully, eigenvals(:, ii, iSpin), Ef - gaussianWidth * kT)
          nLastFully = max(nLastFully, 1)
          call search_asc_real_geq(nLastPartially, eigenvals(:, ii, iSpin), Ef + gaussianWidth * kT)
          nLastPartially = min(nLastPartially + 1, size(eigenvals, dim=1))

          ! Count fully filled states
          call kahan(electronCount, kWeight(ii) * real(nLastFully-1,dp), comp)

          do jj = nLastFully, nLastPartially
            x = (eigenvals(jj, ii, iSpin) - Ef) * beta
            call hX(hermites, MPorder*2, x)
            occ = 0.5_dp * erfcwrap(x)
            do kk = 1, MPorder
              occ = occ + A(kk) * hermites(2 * kk - 1) * exp(-x**2)
            end do
            call kahan(electronCount, occ * kWeight(ii), comp)
          end do

        end do
      end do

    else

      ! Fermi function
      do iSpin = 1, size(eigenvals, dim=3)
        do ii = 1, size(kWeight)

          call search_asc_real_geq(nLastFully, eigenvals(:, ii, iSpin), Ef - fermiWidth * kT)
          nLastFully = max(nLastFully, 1)
          call search_asc_real_geq(nLastPartially, eigenvals(:, ii, iSpin), Ef + fermiWidth * kT)
          nLastPartially = min(nLastPartially + 1, size(eigenvals, dim=1))

          ! Count fully filled states
          call kahan(electronCount, kWeight(ii) * real(nLastFully-1,dp), comp)

          do jj = nLastFully, nLastPartially
            x = ( eigenvals(jj, ii, iSpin) - Ef ) * beta
            ! Where the compiler does not handle inf gracefully, trap the exponential function for
            ! small input values
          #:if EXP_TRAP
            if (x <= mExpArg) then
              call kahan(electronCount, kWeight(ii) / (1.0_dp + exp(x)), comp)
            end if
          #:else
            call kahan(electronCount, kWeight(ii) / (1.0_dp + exp(x)), comp)
          #:endif
          end do

        end do
      end do

    end if

    electronCount = electronCount + comp

  end function electronCount


  !> Calculate filling and TS for the given eigenspectrum and distribution function and Fermi
  !! energy, for two spin channels
  !!
  !! Ref: G. Kresse and J. Furthm&uuml;ller, Phys. Rev. B vol 54, pp 11169 (1996).
  !! Ref: M. Methfessel and A. T. Paxton,, Phys. Rev. B vol 40, pp 3616 (1989).
  !! Ref: F. Wagner, Th. Laloyaux and M. Scheffler, Phys. Rev. B, vol 57 pp 2102 (1998).
  subroutine electronFill(Eband, filling, TS, E0, Ef, eigenvals, kT, distrib, kWeights)

    !> Band structure energy at T
    real(dp), intent(out) :: Eband(:)

    !> Electron occupancies
    real(dp), intent(out) :: filling(:,:,:)

    !> Entropy times temperature
    real(dp), intent(out) :: TS(:)

    !> Band structure energy extrapolated to T=0K
    real(dp), intent(out) :: E0(:)

    !> Fermi energy for given distribution
    real(dp), intent(in) :: Ef

    !> The eigenvalues of the levels, 1st index is energy 2nd index is k-point and 3nd index is spin
    real(dp), intent(in) :: eigenvals(:,:,:)

    !> Thermal energy in atomic units
    real(dp), intent(in) :: kT

    !> Choice of distribution functions, currently Fermi, Gaussian and Methfessle-Paxton
    !! supported. The MP order is defined symbolically, so (Methfessel + 2) gives the 2nd order M-P
    !! scheme
    integer, intent(in) :: distrib

    !> The k-point weightings
    real(dp), intent(in) :: kWeights(:)

    integer :: ii, iSpin, jj , kk, MPorder, nKpts, nLastFully, nLastPartially
    real(dp), allocatable :: A(:), hermites(:)
    real(dp) :: beta, compensateBND, compensateTS, occ, x

    @:ASSERT(size(filling, dim=3) == size(Eband))
    @:ASSERT(size(filling, dim=3) == size(TS))
    @:ASSERT(size(filling, dim=3) == size(E0))

    nKpts = size(kWeights)

    beta = 1.0_dp / kT

    Eband(:) = 0.0_dp
    TS(:) = 0.0_dp
    filling(:,:,:) = 0.0_dp
    E0(:) = 0.0_dp

    if (distrib /= fillingTypes%Fermi) then
      ! The Gaussian and Methfessel-Paxton broadening functions first

      MPorder = distrib - fillingTypes%Methfessel
      allocate(A(0:MPorder))
      allocate(hermites(0 : 2 * MPorder))
      call Aweights(A, MPorder)

      do iSpin = 1, size(eigenvals,dim=3)

        compensateBND = 0.0_dp
        compensateTS = 0.0_dp

        do ii = 1, nKpts

          call search_asc_real_geq(nLastFully, eigenvals(:, ii, iSpin), Ef - gaussianWidth * kT)
          nLastFully = max(nLastFully, 1)
          call search_asc_real_geq(nLastPartially, eigenvals(:, ii, iSpin), Ef + gaussianWidth * kT)
          nLastPartially = min(nLastPartially + 1, size(eigenvals, dim=1))

          filling(:nLastFully - 1, ii, iSpin) = 1.0_dp

          do jj = 1, nLastFully - 1
            call kahan(Eband(iSpin), kWeights(ii) * (filling(jj, ii, iSpin)&
                & * eigenvals(jj, ii, iSpin)), compensateBND)
          end do

          do jj = nLastFully, nLastPartially

            x = (eigenvals(jj, ii, iSpin) - Ef) * beta
            call hX(hermites, MPorder * 2, x)

            ! Gauusian broadened occupancy
            occ = 0.5_dp * erfcwrap(x)

            ! Gaussian broadening entropy
            call kahan(TS(iSpin), kWeights(ii) * 0.5_dp * exp(-x**2) / sqrt(pi), compensateTS)

            ! Methfessel-Paxton occupation sum
            do kk = 1, MPorder
              occ = occ + A(kk) * hermites(2 * kk - 1) * exp(-x**2)
            end do
            filling(jj, ii, iSpin) = occ

            ! Accumulate the band-structure energy
            call kahan(Eband(iSpin), kWeights(ii) * filling(jj, ii, iSpin)&
                & * eigenvals(jj, ii, iSpin), compensateBND)

            ! Methfessel-Paxton broadening entropy
            do kk = 1, MPorder
              call kahan(TS(iSpin), kWeights(ii) * 0.5_dp * A(kk) * hermites(2 * kk) * exp(-x**2),&
                  & compensateTS)
            end do

          end do

          filling(nLastPartially +1:, ii, iSpin) = 0.0_dp

        end do

        Eband(iSpin) = Eband(iSpin) + compensateBND
        TS(iSpin) = TS(iSpin) + compensateTS

      end do

      TS(:) = TS * kT
      E0(:) = (real(MPorder + 1, dp) * (Eband - TS) + Eband) / real(MPorder + 2, dp)

    else

      ! Fermi function
      do iSpin = 1, size(eigenvals, dim=3)

        compensateBND = 0.0_dp
        compensateTS = 0.0_dp

        do ii = 1, nKpts

          call search_asc_real_geq(nLastFully, eigenvals(:, ii, iSpin), Ef - fermiWidth * kT)
          nLastFully = max(nLastFully, 1)
          call search_asc_real_geq(nLastPartially, eigenvals(:, ii, iSpin), Ef + fermiWidth * kT)
          nLastPartially = min(nLastPartially + 1, size(eigenvals, dim=1))

          filling(:nLastFully - 1, ii, iSpin) = 1.0_dp

          do jj = 1, nLastFully - 1
            call kahan(Eband(iSpin), kWeights(ii) * (filling(jj, ii, iSpin)&
                & * eigenvals(jj, ii, iSpin)), compensateBND)
          end do

          do jj = nLastFully, nLastPartially
            x = (eigenvals(jj, ii, iSpin) - Ef) * beta
            ! Where the compiler does not handle inf gracefully, trap the exponential function for
            ! small values
          #:if EXP_TRAP
            if (x > mExpArg) then
              filling(jj, ii, iSpin) = 0.0_dp
            else
              filling(jj, ii, iSpin) = 1.0_dp / (1.0_dp + exp(x))
            endif
          #:else
            filling(jj, ii, iSpin) = 1.0_dp / (1.0_dp + exp(x))
          #:endif

            call kahan(Eband(iSpin), kWeights(ii) * (filling(jj, ii, iSpin)&
                & * eigenvals(jj, ii, iSpin)), compensateBND)

            if (filling(jj, ii, iSpin) > epsilon(0.0_dp) .and.&
                & filling(jj, ii, iSpin) < (1.0_dp - epsilon(1.0_dp))) then
              ! Fermi-Dirac entropy
              call kahan(TS(iSpin),&
                  & - kWeights(ii) * (filling(jj, ii, iSpin) * log(filling(jj, ii, iSpin))&
                  & + (1.0_dp - filling(jj, ii, iSpin)) * log(1.0_dp - filling(jj, ii, iSpin))),&
                  & compensateTS)
            end if

          end do

          filling(nLastPartially +1:, ii, iSpin) = 0.0_dp

        end do

        Eband(iSpin) = Eband(iSpin) + compensateBND
        TS(iSpin) = TS(iSpin) + compensateTS

      end do

      TS(:) = TS * kT
      E0(:) = Eband - 0.5_dp * TS

    end if

  end subroutine electronFill


  !> Calculate the weighting factors for the Methfessel-Paxton smearing scheme
  !!
  !! Ref: M. Methfessel and A. T. Paxton, Phys. Rev. B Vol 40, pp 3616 (1989)
  subroutine Aweights(A,n)

    !> Returned weighting values for the scheme
    real(dp), intent(out) :: A(0:)

    !> The required order to calculate A_n up to
    integer, intent(in) :: n

    real(dp) :: nbang(0:n)
    integer ii

    @:ASSERT(n >= 0)
    @:ASSERT(size(A) >= n)
    A(:) = 0.0_dp
    call fact(nbang, n)
    do ii = 0, n
      A(ii) = real((-1)**ii, dp) / (nbang(ii) * real(4**ii, dp) * sqrt(pi))
    end do

  end subroutine Aweights


  !> Middle gap position, assuming aufbau principle for the filling
  function middleGap(eigenvals, kWeight, nElectrons)

    !> Eigenvalues of states
    real(dp), intent(in) :: eigenvals(:,:,:)

    !> Weights of k-points
    real(dp), intent(in) :: kWeight(:)

    !> Number of electrons to fill in
    real(dp), intent(in) :: nElectrons

    !> Resulting mid gap position
    real(dp) :: middleGap

    integer, allocatable :: tmpIndx(:)
    integer :: size1, size2
    integer :: ind, iLev, iOrb, iKpt, iSpin, jOrb, jKpt, jSpin
    real(dp) :: nElec

    allocate(tmpIndx(size(eigenvals)))
    call index_heap_sort(tmpIndx, reshape(eigenvals, [size(eigenvals)]))
    size1 = size(eigenvals, dim=1)
    size2 = size(eigenvals, dim=2)
    ind = 1
    nElec = 0.0_dp
    do while (nElec < nElectrons)
      iLev = tmpIndx(ind)
      iOrb = mod(iLev - 1, size1) + 1
      iKpt = mod((iLev - 1) / size1, size2) + 1
      iSpin = (iLev - 1) / (size1 * size2) + 1
      nElec = nElec + kWeight(iKpt)
      ind = ind + 1
    end do

    ! just in case the system has all levels filled, but eventually this means Ef has to be above
    ! last eigenvalue:
    ind = min(size(eigenvals), ind)

    iLev = tmpIndx(ind)
    jOrb = mod(iLev - 1, size1) + 1
    jKpt = mod((iLev - 1) / size1, size2) + 1
    jSpin = (iLev - 1) / (size1 * size2) + 1
    middleGap = 0.5_dp * (eigenvals(jOrb, jKpt, jSpin) + eigenvals(iOrb, iKpt, iSpin))

  end function middleGap

end module dftbp_dftb_etemp
