!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Helper routines to handle input to solvation models
module dftbp_solvinput
  use dftbp_accuracy, only : dp
  use dftbp_born, only : TGeneralizedBorn, TGBInput, init
  use dftbp_solvation, only : TSolvation
  implicit none
  private

  public :: TSolvationInp, createSolvationModel, writeSolvationInfo

  !> Collection of solvation input data
  type :: TSolvationInp

    !> Generalized Born input data
    type(TGBInput), allocatable :: GBInp

  end type TSolvationInp

  !> Wrapper to create a solvation model from its respective input
  interface createSolvationModel
    module procedure :: createGeneralizedBornModel
  end interface

contains

  !> Print the solvation model used
  subroutine writeSolvationInfo(unit, solvation)

    !> Formatted unit for IO
    integer, intent(in) :: unit

    !> Solvation model
    class(TSolvation), intent(in) :: solvation

    select type(solvation)
    type is(TGeneralizedBorn)
      write(unit, '(a)') "Using generalized Born solvation model"
    end select

  end subroutine writeSolvationInfo

  !> Wrapper to create a generalized Born model
  subroutine createGeneralizedBornModel(solvation, input, nAtom, species0, speciesNames, latVecs)

    !> Generic solvation model
    class(TSolvation), allocatable, intent(out) :: solvation

    !> Input to setup the solvation model
    type(TGBInput), intent(in) :: input

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    type(TGeneralizedBorn), allocatable :: model

    allocate(model)

    call init(model, input, nAtom, species0, speciesNames, latVecs)

    call move_alloc(model, solvation)

  end subroutine createGeneralizedBornModel

end module dftbp_solvinput
