!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for linear response calculations using perturbation methods
module dftbp_perturbation
  use dftbp_accuracy, only : dp
  use dftbp_commontypes, only : TOrbitals
  use dftbp_rotateDegenerateOrbs, only : TDegeneracyTransform
  use dftbp_environment, only : TEnvironment
  use dftbp_sparse2dense
  use dftbp_densedescr
  use dftbp_finitethelper, only : dEfda, dEida
  use dftbp_scalapackfx
#:if WITH_MPI
  use dftbp_mpifx
#:endif
#:if WITH_SCALAPACK
  use dftbp_scalafxext
#:else
  use dftbp_blasroutines
#:endif
  implicit none

  private
  public :: perturbation

#:set ROUTINE_TYPES = [('real'), ('complex')]

  interface perturbation
    !module procedure sparse
  #:for NAME in ROUTINE_TYPES
    module procedure densePerturb${NAME}$
  #:endfor
  end interface perturbation

  !> small complex value for frequency dependent perturbation cases
  complex(dp), parameter :: eta = (0.0_dp,1.0E-8_dp)

contains

#:for NAME in ROUTINE_TYPES
  !> Evaluate perturbation for ${NAME}$ case of a dense matrix
  subroutine densePerturb${NAME}$(dEi, dFilling, dEf, env, dH, ei, psi, filling, orb,&
      & tMetallic, iFracFilled, tempElec, omega, dS, ePsi, dPsi, dRho)

    !> Perturbed eigenvalues
    real(dp), intent(out) :: dEi(:)

    !> derivative of occupations
    real(dp), intent(out) :: dFilling(:)

    !> derivative of Fermi energy
    real(dp), intent(out) :: dEf

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Perturbation matrix
    ${NAME}$(dp), intent(in) :: dH(:,:)

    !> unperturbed eigenvalues
    real(dp), intent(in) :: ei(:)

    !> Wavefunctions
    ${NAME}$(dp), intent(in) :: psi(:,:)

    !> unperturbed occupations
    real(dp), intent(in) :: filling(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Is this a metallic system, i.e. fractional occupations
    logical, intent(in) :: tMetallic

    !> range of fractionally filled states, ordered as first (partly) empty state and last (partly)
    !> filled state
    integer, intent(in):: iFracFilled(2)

    !> Electronic temperature
    real(dp), intent(in) :: tempElec

    !> Optional driving frequency
    ${NAME}$(dp), intent(in), optional :: omega

    !> Optional basis overlap perturbation
    ${NAME}$(dp), intent(in), optional :: dS(:,:)

    !> Optional energy weighted eigenvectors
    ${NAME}$(dp), intent(in), optional :: ePsi(:,:)

    !> Perturbed wavefunctions
    ${NAME}$(dp), intent(inout), optional :: dPsi(:,:)

    !> Perturbed density matrix
    ${NAME}$(dp), intent(inout), optional :: dRho(:,:)

    type(TDegeneracyTransform) :: transform
    logical :: tDegenerate
    integer :: iOrb, nOrb, iFilled, iEmpty
    ${NAME}$(dp), allocatable :: Umat(:,:), psiTilde(:,:), work(:,:)
  #:if NAME == 'complex'
    complex(dp), parameter :: zero = (0.0_dp, 0.0_dp)
  #:else
    real(dp), parameter :: zero = 0.0_dp
  #:endif

  @:ASSERT(present(dS) .eqv. present(ePsi))

    nOrb = size(psi,dim=2)

    call transform%init()

    allocate(Umat(nOrb,nOrb))

    ! form H' |c>
  #:if NAME == 'complex'
    call hemm(Umat, 'l', dH, psi)
  #:else
    call symm(Umat, 'l', dH, psi)
  #:endif

    if (present(dS)) then
      ! form extra (- e S') |c> term if overlap matrix is changing
  #:if NAME == 'complex'
      call hemm(Umat, 'l', dS, ePsi, alpha=(-1.0_dp,0.0_dp), beta=(1.0_dp,0.0_dp))
  #:else
      call symm(Umat, 'l', dS, ePsi, alpha=-1.0_dp, beta=1.0_dp)
  #:endif
    end if

    ! form <c| H' |c> or <c| (H' - e S') |c>
    allocate(work(nOrb,nOrb))
  #:if NAME == 'complex'
    work(:,:) = matmul(transpose(conjg(psi)), Umat)
  #:else
    work(:,:) = matmul(transpose(psi), Umat)
  #:endif
    Umat(:,:) = work

    call transform%generateUnitary(Umat, ei, tDegenerate)

    if (tDegenerate) then
      call transform%degenerateTransform(Umat)
      allocate(psiTilde(nOrb,nOrb))
      psiTilde = psi(:,:)
      call transform%applyUnitary(psiTilde)
    end if

    ! diagonal elements of Umat are now derivatives of eigenvalues
    do iOrb = 1, nOrb
      dEi(iOrb) = Umat(iOrb,iOrb)
    end do

    if (present(dS)) then
      ! needed for diagonal element
      if (tDegenerate) then
      #:if NAME == 'complex'
        call hemm(work, 'l', dS, psiTilde)
        work(:,:) = work * conjg(psiTilde)
      #:else
        call symm(work, 'l', dS, psiTilde)
        work(:,:) = work * psiTilde
      #:endif
      else
      #:if NAME == 'complex'
        call hemm(work, 'l', dS, psi)
        work(:,:) = work * conjg(psi)
      #:else
        call symm(work, 'l', dS, psi)
        work(:,:) = work * psi
      #:endif
      end if

      ! Form actual perturbation U matrix for eigenvectors by weighting the elements
      do iFilled = 1, nOrb
        do iEmpty = 1, nOrb
          if (iFilled == iEmpty) then
            Umat(iFilled, iFilled) = -0.5_dp * real(sum(work(:, iFilled)),dp)
          else
            if (.not.transform%degenerate(iFilled,iEmpty)) then
              Umat(iEmpty, iFilled) = Umat(iEmpty, iFilled) / (ei(iFilled) - ei(iEmpty))
            else
              Umat(iEmpty, iFilled) = zero
            end if
          end if
        end do
      end do

    else

      ! Form actual perturbation U matrix for eigenvectors by weighting the elements
      do iFilled = 1, nOrb
        do iEmpty = 1, nOrb
          if (iFilled == iEmpty) then
            ! element is already correct
          else
            if (.not.transform%degenerate(iFilled,iEmpty)) then
              Umat(iEmpty, iFilled) = Umat(iEmpty, iFilled) / (ei(iFilled) - ei(iEmpty))
            else
              Umat(iEmpty, iFilled) = zero
            end if
          end if
        end do
      end do

    end if

    ! calculate the derivatives of the eigenvectors
    if (tDegenerate) then
      work(:, :) = matmul(psiTilde, Umat)
    else
      work(:, :) = matmul(psi, Umat)
    end if

    if (present(dPsi)) then
      dPsi(:,:) = work
    end if

    dFilling(:) = 0.0_dp
    dEf = 0.0_dp
    if (tMetallic) then
      dEf = dEfda(filling, dEi)
      call dEida(dFilling, filling, dEi, tempElec)
    end if

    if (present(dRho)) then

      dRho(:,:) = zero

      do iOrb = 1, nOrb
        work(:,iOrb) = filling(iOrb) * work(:,iOrb)
      end do

      ! form the derivative of the density matrix
      if (tDegenerate) then
        dRho(:,:) = matmul(work, transpose(psiTilde)) + matmul(psiTilde, transpose(work))
      else
        dRho(:,:) = matmul(work, transpose(psi)) + matmul(psi, transpose(work))
      end if

      if (tMetallic) then

      @:ASSERT(iFracFilled(1) < iFracFilled(2))

        if (present(dPsi)) then
          work(:,iFracFilled(1):iFracFilled(2)) = dPsi(:,iFracFilled(1):iFracFilled(2))
        else
          ! re-calculate the derivatives of the eigenvectors, as they were overwritten above
          if (tDegenerate) then
            work(:,iFracFilled(1):iFracFilled(2)) = matmul(psiTilde,&
                & Umat(:,iFracFilled(1):iFracFilled(2)))
          else
            work(:,iFracFilled(1):iFracFilled(2)) = matmul(psi,&
                & Umat(:,iFracFilled(1):iFracFilled(2)))
          end if
        end if

        ! reuse Umat as a temporary workspace
        Umat(:,iFracFilled(1):iFracFilled(2)) = work(:,iFracFilled(1):iFracFilled(2))

        do iOrb = iFracFilled(1), iFracFilled(2)
          Umat(:,iOrb) = dFilling(iOrb) * Umat(:,iOrb)
        end do

      #:if NAME == 'complex'
        dRho(:,:) = dRho + 0.5_dp * (&
            & matmul(Umat(:, iFracFilled(1):iFracFilled(2)),&
            & transpose(conjg(work(:, iFracFilled(1):iFracFilled(2)))))&
            & + matmul(work(:, iFracFilled(1):iFracFilled(2)),&
            & transpose(conjg(Umat(:, iFracFilled(1):iFracFilled(2))))) )
      #:else
        dRho(:,:) = dRho + 0.5_dp * (&
            & matmul(Umat(:, iFracFilled(1):iFracFilled(2)),&
            & transpose(work(:, iFracFilled(1):iFracFilled(2))))&
            & + matmul(work(:, iFracFilled(1):iFracFilled(2)),&
            & transpose(Umat(:, iFracFilled(1):iFracFilled(2)))) )
      #:endif

      end if

    end if

  end subroutine densePerturb${NAME}$
#:endfor

end module dftbp_perturbation
