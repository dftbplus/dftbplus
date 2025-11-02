!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
!> Dispatches to either the GPU or CPU implementation, and contains the OMP parallel CPU implementation.
module dftbp_wavegrid_molorb_parallel
  use dftbp_wavegrid_molorb_types, only : TCalculationContext, TPeriodicParams, TSystemParams
  use dftbp_wavegrid_basis, only: TOrbital, realTessY
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error
#:if WITH_CUDA
  use dftbp_wavegrid_molorb_offloaded, only : evaluateCuda
#:endif
  implicit none
  private

  public :: evaluateParallel

contains


  !> Returns the values of several molecular orbitals on grids.
  !> This dispatches to either the CPU or GPU implementation, and handles total charge calculation using occupationVec if present.
  subroutine evaluateParallel(system, periodic, kIndexes, phases, orbitals, &
      & ctx, eigVecsReal, eigVecsCmpl,  valueReal, valueCmpl, occupationVec)

    !> System geometry and composition
    type(TSystemParams), intent(in) :: system

    !> Periodic boundary conditions data
    type(TPeriodicParams), intent(in) :: periodic

    !> Index of the k-points for each orbital in kPoints
    integer, intent(in) :: kIndexes(:)

    !> Phase factors for periodic images
    complex(dp), intent(in) :: phases(:,:)

    !> Basis set data in AoS format
    class(TOrbital), intent(in) :: orbitals(:)

    !> Calculation control flags
    type(TCalculationContext), intent(in) :: ctx

    !> Real eigenvectors, or null-array
    real(dp), intent(in) :: eigVecsReal(:,:)

    !> Complex eigenvectors, or null-array
    complex(dp), intent(in) :: eigVecsCmpl(:,:)

    !> Contains the real grid on exit
    real(dp), intent(out) :: valueReal(:,:,:,:)

    !> Contains the complex grid on exit
    complex(dp), intent(out) :: valueCmpl(:,:,:,:)

    !> Used for calculating total charge if ctx%calcTotalChrg is set.
    !> Done by summing the squared states weighted by occupationVec(nEig).
    !> Output valueReal and valueCmpl will collapse to one slice in the last dimension, (x,y,z,1).
    real(dp), intent(in), optional :: occupationVec(:)

    !! Variables for total charge calculation
    real(dp), allocatable :: coeffVecsReal(:,:)
    complex(dp), allocatable :: coeffVecsCmpl(:,:)

    ! To reduce repeated multiplications we bake the totalCharge occupations
    ! directly into the eigenvector coefficients.
    call prepareCoefficients(ctx, eigVecsReal, eigVecsCmpl, occupationVec, coeffVecsReal, coeffVecsCmpl)
    
    ! Dispatch to CPU / GPU implementation
    if (ctx%runOnGPU) then
      #:if WITH_CUDA
        print *, "Wavegrid: running on GPU using CUDA"
        call evaluateCuda(system, orbitals, periodic, kIndexes, phases, ctx, &
            & coeffVecsReal, coeffVecsCmpl, valueReal, valueCmpl)
      #:else
        call error("Wavegrid: GPU offloaded molorb requested, but compiled without CUDA support.")
      #:endif
    else ! CPU implementation
      #:if WITH_OMP
        print *, "Wavegrid: running OMP parallel on CPU"
      #:else
        print *, "Wavegrid: missing OMP, running serially on CPU"
      #:endif
      call evaluateOMP(system, orbitals, periodic, kIndexes, phases, ctx, &
            & coeffVecsReal, coeffVecsCmpl, valueReal, valueCmpl)
    end if

  end subroutine evaluateParallel


  subroutine evaluateOMP(system, orbitals, periodic, kIndexes, phases, ctx, &
      & eigVecsReal, eigVecsCmpl, valueReal, valueCmpl)

    !> System
    type(TSystemParams), intent(in) :: system

    !> Basis set
    class(TOrbital), intent(in) :: orbitals(:)

    !> Periodic boundary conditions
    type(TPeriodicParams), intent(in) :: periodic
    
    !> K-point indexes (for phase factor selection)
    integer, intent(in) :: kIndexes(:)

    !> Phase factors
    complex(dp), intent(in) :: phases(:, :)

    !> Calculation flags
    type(TCalculationContext), intent(in) :: ctx

    !> Real Eigenvectors
    real(dp), intent(in) :: eigVecsReal(:, :)

    !> Complex Eigenvectors (if not real)
    complex(dp), intent(in) :: eigVecsCmpl(:, :)

    !> Real output grid (if real input or total charge calculation)
    real(dp), intent(out) :: valueReal(:, :, :, :)

    !> Complex output grid (if complex input)
    complex(dp), intent(out) :: valueCmpl(:, :, :, :)

    !! Thread private variables
    integer ::  ind, iSpecies
    real(dp) :: xyz(3), diff(3), frac(3)
    real(dp) :: invR, r, rSq, val, radialVal
    !! Variables for inplace charge calculation
    real(dp), allocatable :: orbValsPerPointReal(:)
    complex(dp), allocatable :: orbValsPerPointCmpl(:)
    !! Loop Variables
    integer :: i1, i2, i3, iAtom, iOrb, iM, iL, iCell
    integer :: gridDimensions(4)

    if (ctx%isRealInput .or. ctx%calcTotalChrg) then
      valueReal = 0.0_dp
      gridDimensions = shape(valueReal)
    else
      valueCmpl = 0.0_dp
      gridDimensions = shape(valueCmpl)
    end if

    @:ASSERT(.not.(.not. ctx%isRealInput .and. ctx%calcAtomicDensity))
    @:ASSERT(.not.(ctx%calcTotalChrg .and. ctx%calcAtomicDensity))


    !$omp parallel private(i1, i2, i3, iCell, iAtom, iOrb, iL, iM, xyz, frac, diff, &
    !$omp&              r, invR, val, radialVal, ind, iSpecies, rSq, &
    !$omp&              orbValsPerPointReal, orbValsPerPointCmpl)
    if (ctx%calcTotalChrg) then
      if (ctx%isRealInput) then
        allocate(orbValsPerPointReal(size(eigVecsReal, dim=2)))
      else
        allocate(orbValsPerPointCmpl(size(eigVecsCmpl, dim=2)))
      end if
    end if

    !$omp do collapse(3)
    lpI3: do i3 = 1, gridDimensions(3)
      lpI2: do i2 = 1, gridDimensions(2)
        lpI1: do i1 = 1, gridDimensions(1)
            xyz(:) = system%origin(:) + real(i1 - 1, dp) * system%gridVecs(:, 1) &
                                    & + real(i2 - 1, dp) * system%gridVecs(:, 2) &
                                    & + real(i3 - 1, dp) * system%gridVecs(:, 3)

            ! Map grid coordinates into unit cell
            if (periodic%isPeriodic) then
              frac(:) = matmul(xyz, periodic%recVecs2pi)
              xyz(:) = matmul(periodic%latVecs, frac - real(floor(frac), dp))
            end if

            if (ctx%calcTotalChrg) then
              if (ctx%isRealInput) then
                orbValsPerPointReal(:) = 0.0_dp
              else
                orbValsPerPointCmpl(:) = 0.0_dp
              end if
            end if

            ! Get contribution from every atom in every cell for current point
            lpCell: do iCell = 1, size(system%coords, dim=3)
              ind = 0
              lpAtom: do iAtom = 1, size(system%coords, dim=2)
                iSpecies = system%species(iAtom)
                diff(:) = xyz - system%coords(:, iAtom, iCell)
                rSq = dot_product(diff, diff)

                lpOrb: do iOrb = system%iStos(iSpecies), system%iStos(iSpecies + 1) - 1
                  iL = orbitals(iOrb)%angMom
                  ! Calculate wave function only if atom is inside the cutoff
                  if (rSq > orbitals(iOrb)%cutoffSq) then
                    ind = ind + 2*iL + 1
                    cycle lpOrb
                  end if
                  r = sqrt(rSq)

                  radialVal = orbitals(iOrb)%getRadial(r)
                  
                  ! Only calculate inverse once
                  invR = 0.0_dp
                  if (r > epsilon(1.0_dp)) then
                    invR = 1.0_dp / r
                  end if 

                  lpM : do iM = -iL, iL
                    ind = ind + 1
                    val = radialVal * realTessY(iL, iM, diff, invR)
                    if (ctx%calcAtomicDensity) val = val * val

                    if (ctx%calcTotalChrg) then
                      if (ctx%isRealInput) then
                        orbValsPerPointReal(:) = orbValsPerPointReal(:) + val * eigVecsReal(ind, :)
                      else ! Complex
                        orbValsPerPointCmpl(:) = orbValsPerPointCmpl(:) + val &
                            & * phases(iCell, kIndexes(:)) * eigVecsCmpl(ind, :)
                      end if
                    else
                      if (ctx%isRealInput) then
                        valueReal(i1, i2, i3, :) = valueReal(i1, i2, i3, :) + val * eigVecsReal(ind, :)
                      else ! Complex
                        valueCmpl(i1, i2, i3, :) = valueCmpl(i1, i2, i3, :) + val &
                            & * phases(iCell, kIndexes(:)) *  eigVecsCmpl(ind, :)
                      end if
                    end if
                  end do lpM
                end do lpOrb
              end do lpAtom
            end do lpCell

            if (ctx%calcTotalChrg) then
              if (ctx%isRealInput) then
                valueReal(i1, i2, i3, 1) = sum(orbValsPerPointReal(:)**2)
              else ! Complex
                valueReal(i1, i2, i3, 1) = sum(abs(orbValsPerPointCmpl(:))**2)
              end if
            end if
        end do lpI1
      end do lpI2
    end do lpI3
    !$omp end do

    if (ctx%calcTotalChrg .and. ctx%isRealInput) then
        deallocate(orbValsPerPointReal)
    else if (ctx%calcTotalChrg) then
        deallocate(orbValsPerPointCmpl)
    end if
    !$omp end parallel
  end subroutine evaluateOMP



  !> Prepare coefficient vectors for calculation by, if required due to total charge calculation,
  !> scaling the eigenvectors by sqrt(occupationVec).
  subroutine prepareCoefficients(ctx, eigVecsReal, eigVecsCmpl, occupationVec, coeffVecReal, coeffVecCmpl)
    type(TCalculationContext), intent(in) :: ctx
    real(dp), intent(in) :: eigVecsReal(:,:)
    complex(dp), intent(in) :: eigVecsCmpl(:,:)
    real(dp), intent(in), optional :: occupationVec(:)
    real(dp), allocatable, intent(out) :: coeffVecReal(:,:)
    complex(dp), allocatable, intent(out) :: coeffVecCmpl(:,:)

    integer :: iEig

    allocate(coeffVecReal(size(eigVecsReal, dim=1), size(eigVecsReal, dim=2)))
    allocate(coeffVecCmpl(size(eigVecsCmpl, dim=1), size(eigVecsCmpl, dim=2)))

    if (ctx%calcTotalChrg) then
      if(ctx%isRealInput) then
        @:ASSERT(size(occupationVec) == size(eigVecsReal, dim=2))
      else
        @:ASSERT(size(occupationVec) == size(eigVecsCmpl, dim=2))
      end if

      do iEig = 1, size(eigVecsReal, dim=2)
        coeffVecReal(:, iEig) = eigVecsReal(:, iEig) * sqrt(occupationVec(iEig))
      end do
      do iEig = 1, size(eigVecsCmpl, dim=2)
        coeffVecCmpl(:, iEig) = eigVecsCmpl(:, iEig) * sqrt(occupationVec(iEig))
      end do
    else
      coeffVecReal = eigVecsReal
      coeffVecCmpl = eigVecsCmpl
    end if

  end subroutine prepareCoefficients

end module dftbp_wavegrid_molorb_parallel
