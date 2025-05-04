!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Helper routines to handle input to solvation models
module dftbp_solvation_solvinput
  use dftbp_common_accuracy, only : dp
  use dftbp_common_status, only : TStatus
  use dftbp_io_message, only : error
  use dftbp_solvation_born, only : TGeneralizedBorn, TGBInput, TGeneralizedBorn_init, writeGeneralizedBornInfo
  use dftbp_solvation_cosmo, only : TCosmo, TCosmoInput, TCosmo_init, writeCosmoInfo
  use dftbp_solvation_sasa, only : TSASACont, TSASAInput, TSASACont_init, writeSASAContInfo
  use dftbp_solvation_solvation, only : TSolvation
  implicit none

  private
  public :: TSolvationInp, createSolvationModel, writeSolvationInfo


  !> Collection of solvation input data
  type :: TSolvationInp

    !> Generalized Born input data
    type(TGBInput), allocatable :: GBInp

    !> Generalized Born input data
    type(TCosmoInput), allocatable :: cosmoInp

    !> Solvent accessible surface area input data
    type(TSASAInput), allocatable :: SASAInp

  end type TSolvationInp


  !> Wrapper to create a solvation model from its respective input
  interface createSolvationModel
    module procedure :: createGeneralizedBornModel
    module procedure :: createCosmoModel
    module procedure :: createSASAModel
  end interface


contains


  !> Print the solvation model used
  subroutine writeSolvationInfo(unit, solvation)

    !> Formatted unit for IO
    integer, intent(in) :: unit

    !> Solvation model
    class(TSolvation), intent(in) :: solvation

    write(unit, '(a, ":", t30)', advance='no') "Solvation"
    select type(solvation)
    type is(TGeneralizedBorn)
      if (solvation%isALPB()) then
        write(unit, '(a)') "analytical linearized Poisson-Boltzmann model"
      else
        write(unit, '(a)') "generalized Born model"
      end if
      call writeGeneralizedBornInfo(unit, solvation)

    type is(TCosmo)
      write(unit, '(a)') "conductor like screening model"
      call writeCosmoInfo(unit, solvation)

    type is(TSASACont)
      write(unit, '(a)') "surface area model"
      call writeSASAContInfo(unit, solvation)

    class default
      write(unit, '(a)') "internal error"
      call error("Unknown solvation model passed")

    end select

  end subroutine writeSolvationInfo


  !> Wrapper to create a generalized Born model
  subroutine createGeneralizedBornModel(solvation, input, nAtom, species0, speciesNames, errStatus,&
      & latVecs)

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

    !> Error status
    type(TStatus), intent(out) :: errStatus

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    type(TGeneralizedBorn), allocatable :: model

    allocate(model)

    call TGeneralizedBorn_init(model, input, nAtom, species0, speciesNames, errStatus, latVecs)
    @:PROPAGATE_ERROR(errStatus)

    call move_alloc(model, solvation)

  end subroutine createGeneralizedBornModel


  !> Wrapper to create a conductor like screening model
  subroutine createCosmoModel(solvation, input, nAtom, species0, speciesNames, errStatus, latVecs)

    !> Generic solvation model
    class(TSolvation), allocatable, intent(out) :: solvation

    !> Input to setup the solvation model
    type(TCosmoInput), intent(in) :: input

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    type(TCosmo), allocatable :: model

    allocate(model)

    call TCosmo_init(model, input, nAtom, species0, speciesNames, errStatus, latVecs)
    @:PROPAGATE_ERROR(errStatus)

    call move_alloc(model, solvation)

  end subroutine createCosmoModel


  !> Wrapper to create a generalized Born model
  subroutine createSASAModel(solvation, input, nAtom, species0, speciesNames, errStatus, latVecs)

    !> Generic solvation model
    class(TSolvation), allocatable, intent(out) :: solvation

    !> Input to setup the solvation model
    type(TSASAInput), intent(in) :: input

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    type(TSASACont), allocatable :: model

    allocate(model)

    call TSASACont_init(model, input, nAtom, species0, speciesNames, errStatus, latVecs)
    @:PROPAGATE_ERROR(errStatus)

    call move_alloc(model, solvation)

  end subroutine createSASAModel


end module dftbp_solvation_solvinput
