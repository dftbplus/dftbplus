!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines related to finite electron temperature, including Fermi, Gaussian and
!> Methfessel-Paxton broadening functions.
!> To do: Add other methods, including possibly Pederson and Jackson method
!> PRB 43, 7312 (1991). Also fix exact occupation for electron numers, using
!> interpolation instead of bisection.
module dftbp_etemp
  use dftbp_assert
  use dftbp_accuracy, only : dp, elecTol, elecTolMax, mExpArg
  use dftbp_errorfunction
  use dftbp_message
  use dftbp_hermite
  use dftbp_sorting
  use dftbp_constants
  use dftbp_factorial, only : fact
  implicit none
  private

  public :: Efilling, electronFill, fillingTypes

  type :: TFillingTypesEnum

    !> Definition of a type of broadening function - Fermi-Dirac in this case
    integer :: Fermi = 0

    !> Definition of a type of broadening function - Gaussian in this case
    integer :: Gaussian = 1

    !> Definition of a type of broadening function - Methfessel-Paxton, for higher orders use
    !> Methfessel + n as a value
    integer :: Methfessel = 2

  end type TFillingTypesEnum

  !> Enumerated filling types.
  type(TFillingTypesEnum), parameter :: fillingTypes = TFillingTypesEnum()

  !> Twice the machine precision
  real(dp), parameter :: epsilon2 = 2.0_dp * epsilon(1.0_dp)

contains


  !> Driver to calculate electron filling, the band-structure energy at T and extrapolated to T=0K,
  !> and the entropy of the electron energy for the Mermin free energy, returning band energy and
  !> entropy for each channel but common Fermi level
  !>
  !> Note: use slices of eigenvalues and other input arrays if you want a different Fermi energy and
  !> filling for different k-points and/or spins.
  !>
  !> Note: If no electrons are present, the Fermi energy is set to zero per default.
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

    !> k-point weightings
    real(dp), intent(in) :: kWeight(:)

    !> Choice of distribution functions, currently Fermi, Gaussian and Methfessle-Paxton
    !> supported. The flags is defined symbolically, so (Methfessel + 2) gives the 2nd order M-P
    !> scheme
    integer, intent(in) :: distrib

    real(dp) :: upperEf, lowerEf
    real(dp) :: nElec
    real(dp) :: nElecMax, nElecMin, maxEig, minEig
    real(dp) :: EfOld

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
    ! Polish resulting root with Newton-Raphson type steps
    if (abs(nElectrons - nElec) > elecTol) then
      if (distrib == fillingTypes%Fermi) then ! only derivs for Fermi so far
        if (abs(derivElectronCount(Ef, eigenvals, kT, distrib, kWeight)) >= epsilon(1.0_dp)) then
          EfOld = Ef
          Ef = Ef - (electronCount(Ef, eigenvals, kT, distrib, kWeight) - nElectrons)&
              & / derivElectronCount(Ef, eigenvals, kT, distrib, kWeight)
          do while(abs(electronCount(EfOld, eigenvals, kT, distrib, kWeight) - nElectrons)&
              & > abs(electronCount(Ef, eigenvals, kT, distrib, kWeight) - nElectrons))
            if (abs(derivElectronCount(Ef, eigenvals, kT, distrib, kWeight))&
                & >= epsilon(1.0_dp)) then
              EfOld = Ef
              Ef = Ef - (electronCount(Ef, eigenvals, kT, distrib, kWeight) - nElectrons)&
                  & / derivElectronCount(Ef, eigenvals, kT, distrib, kWeight)
            else
              exit
            end if
          end do
          Ef = EfOld
        end if
      end if
    end if

    nElec = electronCount(Ef, eigenvals, kT, distrib, kWeight)
    call electronFill(Ebs,filling,TS,E0,Ef,eigenvals,kT,distrib,kWeight)

    ! re-scale to give exact number of electrons, this is a temporay hack
    if (nElec > epsilon(1.0_dp)) then
      filling(:,:,:) = filling * nElectrons / nElec
    end if

  end subroutine Efilling


  !> Calculates the number of electrons for a given Fermi energy and distribution function
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
    !> supported. The flags is defined symbolically, so (Methfessel + 2) gives the 2nd order M-P
    !> scheme
    integer, intent(in) :: distrib

    !> k-point weightings
    real(dp), intent(in) :: kWeight(:)

    integer :: MPorder
    real(dp) :: w
    real(dp), allocatable :: A(:)
    real(dp), allocatable :: hermites(:)
    integer i, j , k, l, ispin
    real(dp) :: occ, x

    w = 1.0_dp/kT
    electronCount=0.0_dp
    if (distrib /= fillingTypes%Fermi) then
      MPorder = distrib - fillingTypes%Methfessel - 1
      allocate(A(0:MPorder))
      allocate(hermites(0:2*MPorder))
      call Aweights(A,MPorder)
      do ispin = 1, size(eigenvals,dim=3)
        do i = 1, size(kWeight)
          do j = 1, size(eigenvals,dim=1)
            if (eigenvals(j,i,ispin)>(Ef-3.0_dp*w)) then
              exit
            else
              electronCount=electronCount+kWeight(i)
            end if
          end do
          do k = j, size(eigenvals,dim=1)
            if (eigenvals(k,i,ispin)>(Ef+3.0_dp*w)) then
              exit
            end if
            x = ( eigenvals(k,i,ispin) - Ef ) / kT
            call hX(hermites,MPorder*2,x)
            occ = 0.5_dp * erfcwrap(x)
            do l=1, MPorder
              occ = occ + A(l) * hermites(2*l-1) * exp(-x**2)
            end do
            electronCount = electronCount + occ * kWeight(i)
          end do
        end do
      end do
    else
      do ispin = 1, size(eigenvals,dim=3)
        do i = 1, size(kWeight)
          do j = 1, size(eigenvals,dim=1)
            x = ( eigenvals(j,i,ispin) - Ef ) / kT
            ! Where the compiler does not handle inf gracefully, trap the exponential function for
            ! small input values
          #:if EXP_TRAP
            if (x <= mExpArg) then
              electronCount = electronCount + kWeight(i)/(1.0_dp + exp(x))
            end if
          #:else
            electronCount = electronCount + kWeight(i)/(1.0_dp + exp(x))
          #:endif
          end do
        end do
      end do
    end if
  end function electronCount


  !> Calculates the derivative of the number of electrons for a given Fermi energy and distribution
  !> function
  !>
  !> To do: support Methfestle-Paxton
  function derivElectronCount(Ef,eigenvals,kT,distrib,kWeight)

    !> Derivative of electrons wrt to Ef
    real(dp) :: derivElectronCount

    !> Fermi energy for given distribution
    real(dp), intent(in) :: Ef

    !> The eigenvalues of the levels, 1st index is energy
    !> 2nd index is k-point and 3nd index is spin
    real(dp), intent(in) :: eigenvals(:,:,:)

    !> Thermal energy in atomic units
    real(dp), intent(in) :: kT

    !> Choice of distribution functions, currently
    integer, intent(in) :: distrib

    !> Fermi supported.
    real(dp), intent(in) :: kWeight(:)

    !> k-point weightings
    real(dp) :: w
    integer i, j, ispin
    real(dp) :: x

    w = 1.0_dp/kT
    derivElectronCount=0.0_dp
    if (distrib /= fillingTypes%Fermi) then
      call error("Fermi distribution only supported")
    else
      do ispin = 1, size(eigenvals,dim=3)
        do i = 1, size(kWeight)
          do j = 1, size(eigenvals,dim=1)
            x = ( eigenvals(j,i,ispin) - Ef ) * w
            if (x<10.0_dp) then
              ! Where the compiler does not handle inf gracefully, trap the exponential function for
              ! small input values
            #:if EXP_TRAP
              if (x <= mExpArg) then
                derivElectronCount = derivElectronCount + &
                    & (w*kWeight(i)) * (exp(x)/((1.0_dp + exp(x))**2))
              end if
            #:else
              derivElectronCount = derivElectronCount + &
                  & (w*kWeight(i)) * (exp(x)/((1.0_dp + exp(x))**2))
            #:endif
            end if
          end do
        end do
      end do
    end if
  end function derivElectronCount


  !> Calculate filling and TS for the given eigenspectrum and distribution function and Fermi
  !> energy, for two spin channels
  !>
  !> Ref: G. Kresse and J. Furthm&uuml;ller, Phys. Rev. B vol 54, pp 11169 (1996).
  !> Ref: M. Methfessel and A. T. Paxton,, Phys. Rev. B vol 40, pp 3616 (1989).
  !> Ref: F. Wagner, Th.\ Laloyaux and M. Scheffler, Phys. Rev. B, vol 57 pp 2102 (1998).
  subroutine electronFill(Eband, filling, TS, E0, Ef, eigenvals, kT, distrib, kWeights)

    !> Band structure energy at T
    real(dp), intent(out) :: Eband(:)

    !> Electron occupancies
    real(dp), intent(out) :: filling(:,:,:)

    !> Entropy * temperature
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
    !> supported. The flags is defined symbolically, so (Methfessel + 2) gives the 2nd order M-P
    !> scheme
    integer, intent(in) :: distrib

    !> k-point weightings
    real(dp), intent(in) :: kWeights(:)

    integer :: MPorder
    integer :: kpts
    real(dp) :: w
    real(dp), allocatable :: A(:)
    real(dp), allocatable :: hermites(:)
    integer :: i, j , k, l, iSpin
    real(dp) :: occ, x

    @:ASSERT(size(filling, dim=3) == size(Eband))
    @:ASSERT(size(filling, dim=3) == size(TS))
    @:ASSERT(size(filling, dim=3) == size(E0))

    kpts = size(kWeights)

    Eband(:) = 0.0_dp
    TS(:) = 0.0_dp
    filling(:,:,:) = 0.0_dp
    w = 1.0_dp / kT
    E0(:) = 0.0_dp

    ! The Gaussian and Methfessel-Paxton broadening functions first
    if (distrib /= fillingTypes%Fermi) then
      MPorder = distrib - fillingTypes%Methfessel -1
      allocate(A(0:MPorder))
      allocate(hermites(0 : 2 * MPorder))
      call Aweights(A, MPorder)
      do iSpin = 1, size(eigenvals,dim=3)
        do i = 1, kpts
          do j = 1, size(eigenvals,dim=1)
            if (eigenvals(j,i,iSpin)>(Ef-3.0_dp*w)) then
              exit
            else
              filling(j,i,iSpin)=1.0_dp
              Eband(iSpin) = Eband(iSpin) + eigenvals(j,i,iSpin)
            end if
          end do
          do k = j, size(eigenvals,dim=1)
            if (eigenvals(k, i, iSpin) > (Ef + 3.0_dp * w)) then
              exit
            end if
            x = (eigenvals(k,i,iSpin) - Ef) / kT
            call hX(hermites, MPorder * 2, x)
            ! Gauusian broadened occupancy
            occ = 0.5_dp * erfcwrap(x)
            ! Gaussian broadening entropy
            TS(iSpin) = TS(iSpin) + kWeights(i) * 0.5_dp * exp(-x**2) / sqrt(pi)
            ! Methfessel-Paxton occupation sum
            do l = 1, MPorder
              occ = occ + A(l) * hermites(2*l-1) * exp(-x**2)
            end do
            filling(k, i, iSpin) = occ
            ! Sum up the band-structure energy, including k-point weight where needed
            Eband(iSpin) = Eband(iSpin)&
                & + kWeights(i) * filling(k, i, iSpin) * eigenvals(k, i, iSpin)
            ! Methfessel-Paxton broadening entropy
            do l = 1, MPorder
              TS(iSpin) = TS(iSpin) &
                  & + kWeights(i) * 0.5_dp * A(l) * hermites(2 * l) * exp(-x**2)
            end do
          end do
        end do
      end do
      TS = TS * kT
      E0(:) = (real(MPorder + 1,dp) * (Eband - TS) + Eband) / real(MPorder + 2, dp)
    else
      do iSpin = 1, size(eigenvals, dim=3)
        do i = 1, kpts
          do j = 1, size(eigenvals, dim=1)
            x = (eigenvals(j, i, iSpin) - Ef) / kT
            ! Where the compiler does not handle inf gracefully, trap the exponential function for
            ! small values
          #:if EXP_TRAP
            if (x > mExpArg) then
              filling(j, i, iSpin) = 0.0_dp
            else
              filling(j, i, iSpin) = 1.0_dp / (1.0_dp + exp(x))
            endif
          #:else
            filling(j, i, iSpin) = 1.0_dp / (1.0_dp + exp(x))
          #:endif
            if (filling(j, i, iSpin) <= elecTol) then
              exit
            end if
            if (filling(j, i, iSpin) > epsilon(0.0_dp) .and.&
                & filling(j, i, iSpin) < (1.0_dp - epsilon(1.0_dp))) then
              ! Fermi-Dirac entropy :
              TS(iSpin) = TS(iSpin)&
                  & - kWeights(i) * (filling(j, i, iSpin)* log(filling(j, i, iSpin))&
                  & + (1.0_dp - filling(j, i, iSpin)) * log(1.0_dp - filling(j, i, iSpin)))
            end if
            Eband(iSpin) = Eband(iSpin) &
                & + kWeights(i) * (filling(j, i, iSpin) * eigenvals(j, i, iSpin))
          end do
        end do
      end do
      TS(:) = TS * kT
      E0(:) = Eband - 0.5_dp * TS
    end if

  end subroutine electronFill


  !> Calculate the weighting factors for the Methfessel-Paxton smearing scheme
  !>
  !> Ref: M. Methfessel and A. T. Paxton, Phys. Rev. B Vol 40, pp 3616 (1989)
  subroutine Aweights(A,n)

    !> returned weighting values for the scheme
    real(dp), intent(out) :: A(0:)

    !> the required order to calculate A_n up to
    integer, intent(in) :: n

    real(dp) :: nbang(0:n)
    integer i
    @:ASSERT(n>=0)
    @:ASSERT(size(A)>=n)
    A(:) = 0.0_dp
    call fact(nbang,n)
    do i = 0, n
      A(i) = real((-1)**i,dp)/(nbang(i)*real(4**i,dp)*sqrt(pi))
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

end module dftbp_etemp
