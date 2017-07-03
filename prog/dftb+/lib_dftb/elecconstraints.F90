!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module elecconstraints
  use assert
  use accuracy, only : dp
  use commontypes, only : TOrbitals
  use angmomentum, only : Loperators

  implicit none

  private
  public :: constrainQ, constrainS, constrainL, constrainJ, constrainMj

  interface constrainQ
    module procedure constrainQ_
  end interface

  interface constrainS
    module procedure constrainS_
  end interface

  interface constrainL
    module procedure constrainL_
  end interface

  interface constrainJ
    module procedure constrainJ_
  end interface

  interface constrainMj
    module procedure constrainMj_
  end interface

contains

  subroutine constrainQ_(shift, qIn, orb, species, conAt, conSh, Qtarget, &
      & V)
    real(dp), intent(inout)     :: shift(:,:,:,:)
    real(dp), intent(in)        :: qIn(:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: species(:)
    integer, intent(in)         :: conAt
    integer, intent(in)         :: conSh
    real(dp), intent(in)        :: Qtarget
    real(dp), intent(in)        :: V

    integer :: iOrb
    real(dp) :: Qshell

    Qshell = sum(qIn(orb%posShell(conSh,species(conAt)): &
        & orb%posShell(conSh+1,species(conAt))-1,conAt,1))

    ! Push q towards required value
    do iOrb = orb%posShell(conSh,species(conAt)), &
        & orb%posShell(conSh+1,species(conAt))-1
      shift(iOrb,iOrb,conAt,1) = shift(iOrb,iOrb,conAt,1) &
          & + V * 0.5_dp*(Qshell - Qtarget)
    end do

  end subroutine constrainQ_

  subroutine constrainS_(shift, qIn, orb, species, conAt, conSh, Starget, &
      & V, vec)
    real(dp), intent(inout)     :: shift(:,:,:,:)
    real(dp), intent(in)        :: qIn(:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: species(:)
    integer, intent(in)         :: conAt
    integer, intent(in)         :: conSh
    real(dp), intent(in)        :: Starget
    real(dp), intent(in)        :: V
    real(dp), intent(in)        :: vec(3)

    integer :: iOrb, nSpin, iSpin
    real(dp) :: Sshell(3), W, vecNorm(3)

    nSpin = size(shift,dim=4)

    vecNorm = vec / sqrt(sum(vec**2))

    Sshell = sum(qIn(orb%posShell(conSh,species(conAt)): &
        & orb%posShell(conSh+1,species(conAt))-1,conAt,2:4),dim=1)

    if (sqrt(sum(Sshell**2)) < 1.0E-8_dp) then
      Sshell = Sshell + 1.0E-8_dp*(/1,1,1/)
    end if

    vecNorm = Sshell  / sqrt(sum(Sshell**2))

    ! Push S towards required value

    w = V * 0.5_dp*(dot_product(Sshell,vecNorm) - Starget)

    do iSpin = 2, nSpin
      do iOrb = orb%posShell(conSh,species(conAt)), &
          & orb%posShell(conSh+1,species(conAt))-1
        shift(iOrb,iOrb,conAt,iSpin) = shift(iOrb,iOrb,conAt,iSpin) &
            & + w * vecNorm(iSpin-1)
      end do
    end do

  end subroutine constrainS_

  !!* @param shift block shift
  !!* @param qBlockSkew Antisymmetric Mulliken block populations for imaginary
  !!* coefficients of Pauli matrics
  !!* @param orb Information about the orbitals in the system.
  !!* @param species Species of the atoms
  !!* @param conAt Atom for constraint
  !!* @param conSh Shell for constraint
  !!* @param Ltarget value of L
  !!* @param V strength of constraint
  !!* @param vec direction of constrain
  subroutine constrainL_(iShift,qBlockSkew, orb, species, conAt, conSh, &
      & Ltarget, V, vec)
    real(dp), intent(inout)     :: iShift(:,:,:,:)
    real(dp), intent(in)        :: qBlockSkew(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: species(:)
    integer, intent(in)         :: conAt
    integer, intent(in)         :: conSh
    real(dp), intent(in)        :: Ltarget
    real(dp), intent(in)        :: V
    real(dp), intent(in)        :: vec(3)

    integer :: ii, iSp, iSh, iOrb, iStart, iEnd
    real(dp), allocatable :: SpeciesL(:,:,:)
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    real(dp), allocatable :: tmpBlock(:,:)
    real(dp) :: Lshell(3), W, vecNorm(3)

    complex(dp), parameter :: i = (0.0_dp,1.0_dp)

    vecNorm = vec / sqrt(sum(vec**2))

    allocate(SpeciesL(orb%mOrb,orb%mOrb,3))
    SpeciesL = 0.0_dp
    allocate(Lz(orb%mOrb,orb%mOrb))
    allocate(Lplus(orb%mOrb,orb%mOrb))

    iSp = species(conAt)
    Lz = 0.0_dp
    Lplus = 0.0_dp
    iSh = orb%angShell(conSh,iSp)
    call loperators(Lplus(1:2*iSh+1,1:2*iSh+1),Lz(1:2*iSh+1,1:2*iSh+1),iSh)
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1, &
        & orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,1) &
        & = aimag(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1, &
        & orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,2) &
        & = -real(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1, &
        & orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,3) &
        & = aimag(Lz(1:2*iSh+1,1:2*iSh+1))

    allocate(tmpBlock(orb%mOrb,orb%mOrb))

    Lshell = 0.0_dp

    iOrb = orb%nOrbSpecies(iSp)
    tmpBlock(:,:) = 0.0_dp
    tmpBlock(1:iOrb,1:iOrb) = qBlockSkew(1:iOrb,1:iOrb,conAt,1) ! identity part
    iStart = orb%posShell(conSh,iSp)
    iEnd = orb%posShell(conSh+1,iSp)-1
    do ii = 1, 3
      Lshell(ii) = &
          & - sum(SpeciesL(iStart:iEnd,iStart:iEnd,ii) &
          &  * transpose(tmpBlock(iStart:iEnd,iStart:iEnd)))
    end do

    if (sqrt(sum(Lshell**2)) < 1.0E-8_dp) then
      Lshell = Lshell + 1.0E-8_dp*(/1,1,1/)
    end if

    vecNorm = Lshell / sqrt(sum(Lshell**2))

    ! Push L towards required value

    w = V * 0.5_dp*(dot_product(lshell,vecNorm)-Ltarget)

    do ii = 1, 3
      iShift(iStart:iEnd,iStart:iEnd,conAt,1) = &
          & iShift(iStart:iEnd,iStart:iEnd,conAt,1) &
          & + w * vecNorm(ii) * SpeciesL(iStart:iEnd,iStart:iEnd,ii)
    end do

  end subroutine constrainL_


  subroutine constrainJ_(shift, qIn, iShift, qBlockSkew, orb, species, &
      & conAt, conSh, Jtarget, V, vec)
    real(dp), intent(inout)     :: shift(:,:,:,:)
    real(dp), intent(in)        :: qIn(:,:,:)
    real(dp), intent(inout)     :: iShift(:,:,:,:)
    real(dp), intent(in)        :: qBlockSkew(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: species(:)
    integer, intent(in)         :: conAt
    integer, intent(in)         :: conSh
    real(dp), intent(in)        :: Jtarget
    real(dp), intent(in)        :: V
    real(dp), intent(in)        :: vec(3)

    integer :: ii, iSp, iSh, iOrb, iStart, iEnd, nSpin, iSpin
    real(dp), allocatable :: SpeciesL(:,:,:)
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    real(dp), allocatable :: tmpBlock(:,:)
    real(dp) :: Lshell(3), Sshell(3), W, vecNorm(3)

    complex(dp), parameter :: i = (0.0_dp,1.0_dp)

    nSpin = size(shift,dim=4)

    vecNorm = vec / sqrt(sum(vec**2))

    allocate(SpeciesL(orb%mOrb,orb%mOrb,3))
    SpeciesL = 0.0_dp
    allocate(Lz(orb%mOrb,orb%mOrb))
    allocate(Lplus(orb%mOrb,orb%mOrb))

    iSp = species(conAt)
    Lz = 0.0_dp
    Lplus = 0.0_dp
    iSh = orb%angShell(conSh,iSp)
    call loperators(Lplus(1:2*iSh+1,1:2*iSh+1),Lz(1:2*iSh+1,1:2*iSh+1),iSh)
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1, &
        & orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,1) &
        & = aimag(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1, &
        & orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,2) &
        & = -real(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1, &
        & orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,3) &
        & = aimag(Lz(1:2*iSh+1,1:2*iSh+1))

    allocate(tmpBlock(orb%mOrb,orb%mOrb))

    Lshell = 0.0_dp

    iOrb = orb%nOrbSpecies(iSp)
    tmpBlock(:,:) = 0.0_dp
    tmpBlock(1:iOrb,1:iOrb) = qBlockSkew(1:iOrb,1:iOrb,conAt,1) ! identity part
    iStart = orb%posShell(conSh,iSp)
    iEnd = orb%posShell(conSh+1,iSp)-1
    do ii = 1, 3
      Lshell(ii) = &
          & - sum(SpeciesL(iStart:iEnd,iStart:iEnd,ii) &
          &  * transpose(tmpBlock(iStart:iEnd,iStart:iEnd)))
    end do

    Sshell = sum(qIn(orb%posShell(conSh,species(conAt)): &
        & orb%posShell(conSh+1,species(conAt))-1,conAt,2:4),dim=1)

    if ( sqrt(sum((lshell + 0.5_dp*Sshell)**2)) < 1.0E-8_dp) then
      Sshell = Sshell + 1.0E-8_dp*(/1,1,1/)
    end if

    vecNorm = (lshell + 0.5_dp*Sshell) / sqrt(sum((lshell + 0.5_dp*Sshell)**2))

    ! Push J towards required value

    w = V * 0.5_dp*(dot_product(lshell,vecNorm)+ &
        & 0.5_dp*dot_product(Sshell,vecNorm) -Jtarget)

    do ii = 1, 3
      iShift(iStart:iEnd,iStart:iEnd,conAt,1) = &
          & iShift(iStart:iEnd,iStart:iEnd,conAt,1) &
          & + w * vecNorm(ii) * SpeciesL(iStart:iEnd,iStart:iEnd,ii)
    end do


    do iSpin = 2, nSpin
      do iOrb = orb%posShell(conSh,species(conAt)), &
          & orb%posShell(conSh+1,species(conAt))-1
        shift(iOrb,iOrb,conAt,iSpin) = shift(iOrb,iOrb,conAt,iSpin) &
            & + w * vecNorm(iSpin-1)
      end do
    end do


  end subroutine constrainJ_

  subroutine constrainMj_(shift, qIn, iShift, qBlockSkew, orb, species, &
      & conAt, conSh, Jtarget, V, vec)
    real(dp), intent(inout)     :: shift(:,:,:,:)
    real(dp), intent(in)        :: qIn(:,:,:)
    real(dp), intent(inout)     :: iShift(:,:,:,:)
    real(dp), intent(in)        :: qBlockSkew(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in)         :: species(:)
    integer, intent(in)         :: conAt
    integer, intent(in)         :: conSh
    real(dp), intent(in)        :: Jtarget
    real(dp), intent(in)        :: V
    real(dp), intent(in)        :: vec(3)

    integer :: ii, iSp, iSh, iOrb, iStart, iEnd, nSpin, iSpin
    real(dp), allocatable :: SpeciesL(:,:,:)
    complex(dp), allocatable :: Lz(:,:)
    complex(dp), allocatable :: Lplus(:,:)
    real(dp), allocatable :: tmpBlock(:,:)
    real(dp) :: Lshell(3), Sshell(3), W, vecNorm(3)

    complex(dp), parameter :: i = (0.0_dp,1.0_dp)

    nSpin = size(shift,dim=4)

    vecNorm = vec / sqrt(sum(vec**2))

    allocate(SpeciesL(orb%mOrb,orb%mOrb,3))
    SpeciesL = 0.0_dp
    allocate(Lz(orb%mOrb,orb%mOrb))
    allocate(Lplus(orb%mOrb,orb%mOrb))

    iSp = species(conAt)
    Lz = 0.0_dp
    Lplus = 0.0_dp
    iSh = orb%angShell(conSh,iSp)
    call loperators(Lplus(1:2*iSh+1,1:2*iSh+1),Lz(1:2*iSh+1,1:2*iSh+1),iSh)
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1, &
        & orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,1) &
        & = aimag(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1, &
        & orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,2) &
        & = -real(Lplus(1:2*iSh+1,1:2*iSh+1))
    speciesL(orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1, &
        & orb%posShell(conSh,iSp):orb%posShell(conSh+1,iSp)-1,3) &
        & = aimag(Lz(1:2*iSh+1,1:2*iSh+1))

    allocate(tmpBlock(orb%mOrb,orb%mOrb))

    Lshell = 0.0_dp

    iOrb = orb%nOrbSpecies(iSp)
    tmpBlock(:,:) = 0.0_dp
    tmpBlock(1:iOrb,1:iOrb) = qBlockSkew(1:iOrb,1:iOrb,conAt,1) ! identity part
    iStart = orb%posShell(conSh,iSp)
    iEnd = orb%posShell(conSh+1,iSp)-1
    do ii = 1, 3
      Lshell(ii) = &
          & - sum(SpeciesL(iStart:iEnd,iStart:iEnd,ii) &
          &  * transpose(tmpBlock(iStart:iEnd,iStart:iEnd)))
    end do

    Sshell = sum(qIn(orb%posShell(conSh,species(conAt)): &
        & orb%posShell(conSh+1,species(conAt))-1,conAt,2:4),dim=1)

    ! Push J towards required value

    w = V * 0.5_dp*(dot_product(lshell,vecNorm)+ &
        & 0.5_dp*dot_product(Sshell,vecNorm) -Jtarget)

    do ii = 1, 3
      iShift(iStart:iEnd,iStart:iEnd,conAt,1) = &
          & iShift(iStart:iEnd,iStart:iEnd,conAt,1) &
          & + w * vecNorm(ii) * SpeciesL(iStart:iEnd,iStart:iEnd,ii)
    end do

    do iSpin = 2, nSpin
      do iOrb = orb%posShell(conSh,species(conAt)), &
          & orb%posShell(conSh+1,species(conAt))-1
        shift(iOrb,iOrb,conAt,iSpin) = shift(iOrb,iOrb,conAt,iSpin) &
            & + w * vecNorm(iSpin-1)
      end do
    end do

  end subroutine constrainMj_

end module elecconstraints
