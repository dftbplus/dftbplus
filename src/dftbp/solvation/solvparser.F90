!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Fills the derived type with the input parameters from an HSD or an XML file.
module dftbp_solvation_solvparser
  use, intrinsic :: ieee_arithmetic, only : ieee_support_inf, ieee_value, ieee_positive_inf
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_atomicrad, only : getAtomicRad
  use dftbp_common_constants, only : Boltzmann, amu__au, kg__au, AA__Bohr
  use dftbp_common_filesystem, only : findFile, getParamSearchPath
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_unitconversion, only : lengthUnits, energyUnits, massUnits, &
      & massDensityUnits, inverseLengthUnits
  use dftbp_dftbplus_specieslist, only : readSpeciesList
  use dftbp_extlibs_lebedev, only : gridSize
  use dftbp_extlibs_xmlf90, only : fnode, string, char, getNodeName
  use dftbp_io_charmanip, only : tolower, unquote
  use dftbp_io_hsdutils, only : getChild, getChildValue, setChild, detailedError, &
      & detailedWarning
  use dftbp_io_hsdutils2, only : convertUnitHsd
  use dftbp_math_bisect, only : bisection
  use dftbp_solvation_born, only : TGBInput, fgbKernel
  use dftbp_solvation_cm5, only : TCM5Input
  use dftbp_solvation_cosmo, only : TCosmoInput, TDomainDecompositionInput
  use dftbp_solvation_gbsafile, only : readParamGBSA
  use dftbp_solvation_sasa, only : TSASAInput
  use dftbp_solvation_solvdata, only : getVanDerWaalsRadiusD3, getVanDerWaalsRadiusCosmo, &
      & getVanDerWaalsRadiusBondi
  use dftbp_solvation_solventdata, only : TSolventData, SolventFromName
  use dftbp_solvation_solvinput, only : TSolvationInp
  use dftbp_type_typegeometry, only : TGeometry
  implicit none

  private
  public :: readSolvation
  public :: readSolvGB, readSolvSASA, readCM5, readSolvCosmo


  real(dp), parameter :: ambientTemperature = 298.15_dp * Boltzmann


contains


  !> Reads in solvation related settings
  subroutine readSolvation(node, geo, input)

    !> Node to parse
    type(fnode), pointer :: node

    !> Geometry, including atomic information
    type(TGeometry), intent(in) :: geo

    !> Solvation data on exit
    type(TSolvationInp), intent(out) :: input

    type(fnode), pointer :: solvModel
    type(string) :: buffer

    call getChildValue(node, "", solvModel)
    call getNodeName(solvModel, buffer)

    select case (char(buffer))
    case default
      call detailedError(node, "Invalid solvation model name.")
    case ("generalizedborn")
      allocate(input%GBInp)
      call readSolvGB(solvModel, geo, input%GBInp)
    case ("cosmo")
      allocate(input%cosmoInp)
      call readSolvCosmo(solvModel, geo, input%cosmoInp)
    case ("sasa")
      allocate(input%SASAInp)
      call readSolvSASA(solvModel, geo, input%SASAInp)
    end select
  end subroutine readSolvation


  !> Reads in generalized Born related settings
  subroutine readSolvGB(node, geo, input)

    !> Node to process
    type(fnode), pointer :: node

    !> Geometry of the current system
    type(TGeometry), intent(in) :: geo

    !> Contains the input for the solvation module on exit
    type(TGBInput), intent(out) :: input

    type(TGBInput), allocatable :: defaults
    type(string) :: buffer, modifier
    type(fnode), pointer :: child, value1, field
    logical :: found, tHBondCorr, tALPB
    real(dp) :: temperature, shift, alphaALPB
    type(string), allocatable :: searchPath(:)
    type(TSolventData) :: solvent
    real(dp), parameter :: alphaDefault = 0.571412_dp
    character(len=:), allocatable :: paramFile, paramTmp

    if (geo%tPeriodic .or. geo%tHelical) then
      call detailedError(node, "Generalized Born model currently not available with the&
         & selected boundary conditions")
    end if

    call getChild(node, "ParamFile", value1, requested=.false.)
    if (associated(value1)) then
      allocate(defaults)
      call getChildValue(node, "ParamFile", buffer, "", child=child)
      paramFile = trim(unquote(char(buffer)))
      call getParamSearchPath(searchPath)
      call findFile(searchPath, paramFile, paramTmp)
      if (allocated(paramTmp)) call move_alloc(paramTmp, paramFile)
      write(stdOut, '(a)') "Reading GBSA parameter file '"//paramFile//"'"
      call readParamGBSA(paramFile, defaults, solvent, geo%speciesNames, node=child)
    else
      call readSolvent(node, solvent)
    end if

    call getChildValue(node, "ALPB", tALPB, .false.)
    if (tALPB) then
      call getChildValue(node, "Alpha", alphaALPB, alphaDefault)
      input%alpbet = alphaALPB / solvent%dielectricConstant
    else
      input%alpbet = 0.0_dp
    end if
    input%keps = (1.0_dp / solvent%dielectricConstant - 1.0_dp) / (1.0_dp + input%alpbet)
    input%dielectricConstant = solvent%dielectricConstant

    call getChildValue(node, "Kernel", buffer, "Still", child=child)
    select case(tolower(unquote(char(buffer))))
    case default
      call detailedError(child, "Unknown interaction kernel: "//char(buffer))
    case("still")
      input%kernel = fgbKernel%still
    case("p16")
      input%kernel = fgbKernel%p16
    end select

    ! shift value for the free energy (usually fitted)
    if (allocated(defaults)) then
      call getChildValue(node, "FreeEnergyShift", shift, defaults%freeEnergyShift, &
          & modifier=modifier, child=field)
    else
      call getChildValue(node, "FreeEnergyShift", shift, modifier=modifier, &
          & child=field)
    end if
    call convertUnitHsd(char(modifier), energyUnits, field, shift)

    ! temperature, influence depends on the reference state
    call getChildValue(node, "Temperature", temperature, ambientTemperature, &
        & modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), energyUnits, field, temperature)

    ! reference state for free energy calculation
    call readReferenceState(node, solvent, temperature, shift, input%freeEnergyShift)

    if (allocated(defaults)) then
      call getChildValue(node, "BornScale", input%bornScale, defaults%bornScale)
      call getChildValue(node, "BornOffset", input%bornOffset, defaults%bornOffset, &
          & modifier=modifier, child=field)
    else
      call getChildValue(node, "BornScale", input%bornScale)
      call getChildValue(node, "BornOffset", input%bornOffset, modifier=modifier, child=field)
    end if
    call convertUnitHsd(char(modifier), lengthUnits, field, input%bornOffset)
    call getChildValue(node, "OBCCorrection", input%obc, [1.00_dp, 0.80_dp, 4.85_dp])

    call getChild(node, "CM5", child, requested=.false.)
    if (associated(child)) then
      allocate(input%cm5Input)
      call readCM5(child, input%cm5Input, geo)
    end if

    call readVanDerWaalsRad(node, geo, input%vdwRad)

    allocate(input%descreening(geo%nSpecies))
    if (allocated(defaults)) then
      call getChildValue(node, "Descreening", value1, "Defaults", child=child)
    else
      call getChildValue(node, "Descreening", value1, child=child)
    end if
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(child, "Unknown method '"//char(buffer)//"' to generate descreening parameters")
    case("defaults")
      if (.not.allocated(defaults)) then
        call detailedError(child, "No defaults available for descreening parameters")
      end if
      call readSpeciesList(value1, geo%speciesNames, input%descreening, &
          & defaults%descreening)
    case("unity")
      input%descreening(:) = 1.0_dp
    case("values")
      call readSpeciesList(value1, geo%speciesNames, input%descreening)
    end select

    call getChildValue(node, "Cutoff", input%rCutoff, 35.0_dp * AA__Bohr, &
        & modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), lengthUnits, field, input%rCutoff)

    call getChild(node, "SASA", value1, requested=.false.)
    if (associated(value1) .or. allocated(defaults)) then
      allocate(input%sasaInput)
      if (.not.associated(value1)) then
        call setChild(node, "SASA", value1)
      end if
      if (allocated(defaults)) then
        call readSolvSASA(value1, geo, input%sasaInput, defaults%sasaInput%probeRad, &
            & defaults%sasaInput%surfaceTension)
      else
        call readSolvSASA(value1, geo, input%sasaInput)
      end if

      if (allocated(defaults)) then
        call getChildValue(node, "HBondCorr", tHBondCorr, &
            & allocated(defaults%hBondPar), child=child)
      else
        call getChildValue(node, "HBondCorr", tHBondCorr, child=child)
      end if

      if (tHBondCorr) then
        allocate(input%hBondPar(geo%nSpecies))
        if (allocated(defaults)) then
          call getChildValue(node, "HBondStrength", value1, "Defaults", child=child)
        else
          call getChildValue(node, "HBondStrength", value1, child=child)
        end if
        call getNodeName(value1, buffer)
        select case(char(buffer))
        case default
          call detailedError(child, "Unknown method '"//char(buffer)//"' to generate H-bond parameters")
        case("defaults")
          if (allocated(defaults)) then
            if (.not.allocated(defaults%hBondPar)) then
              call detailedError(child, "No defaults available for hydrogen bond strengths")
            end if
          else
            call detailedError(child, "No defaults available for hydrogen bond strengths")
          end if
          call readSpeciesList(value1, geo%speciesNames, input%hBondPar, &
              & defaults%hBondPar)
        case("values")
          call readSpeciesList(value1, geo%speciesNames, input%hBondPar)
        end select
      end if
    end if

  end subroutine readSolvGB


  !> Reads in conductor like screening model settings
  subroutine readSolvCosmo(node, geo, input)

    !> Node to process
    type(fnode), pointer :: node

    !> Geometry of the current system
    type(TGeometry), intent(in) :: geo

    !> Contains the input for the solvation module on exit
    type(TCosmoInput), intent(out) :: input

    type(string) :: buffer, modifier
    type(fnode), pointer :: child, value1, field
    real(dp) :: temperature, shift, radScale
    type(TSolventData) :: solvent

    if (geo%tPeriodic .or. geo%tHelical) then
      call detailedError(node, "COSMO solvation currently not available with the&
         & selected boundary conditions")
    end if

    call readSolvent(node, solvent)
    input%dielectricConst = solvent%dielectricConstant
    input%keps = 0.5_dp * (1.0_dp - 1.0_dp/solvent%dielectricConstant)

    ! shift value for the free energy (usually zero)
    call getChildValue(node, "FreeEnergyShift", shift, 0.0_dp, modifier=modifier, &
      & child=field)
    call convertUnitHsd(char(modifier), energyUnits, field, shift)

    ! temperature, influence depends on the reference state
    call getChildValue(node, "Temperature", temperature, ambientTemperature, &
        & modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), energyUnits, field, temperature)

    call readReferenceState(node, solvent, temperature, shift, input%freeEnergyShift)

    call readVanDerWaalsRad(node, geo, input%vdwRad)
    call getChildValue(node, "RadiiScaling", radScale)
    input%vdwRad(:) = input%vdwRad * radScale

    call readAngularGrid(node, input%gridSize)

    call getChildValue(node, "Solver", value1, "DomainDecomposition", child=child)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(child, "Unknown method '"//char(buffer)//"' to solve COSMO equation")
    case("domaindecomposition")
      call readDomainDecomposition(value1, input%ddInput)
    end select

    call getChild(node, "SASA", value1, requested=.false.)
    if (associated(value1)) then
      allocate(input%sasaInput)
      call readSolvSASA(value1, geo, input%sasaInput)
    end if

  end subroutine readSolvCosmo


  subroutine readDomainDecomposition(node, input)

    !> Node to process
    type(fnode), pointer :: node

    !> Input for the domain decomposition algorithm
    type(TDomainDecompositionInput), intent(out) :: input

    type(fnode), pointer :: child

    call getChildValue(node, "MaxMoment", input%lmax, child=child)
    call getChildValue(node, "Regularisation", input%eta, 0.2_dp, child=child)
    call getChildValue(node, "Accuracy", input%conv, child=child)

  end subroutine readDomainDecomposition


  !> Read input data for non-polar surface area solvation model.
  subroutine readSolvSASA(node, geo, input, probeRadDefault, surfaceTensionDefault)

    !> Node to process
    type(fnode), pointer :: node

    !> Geometry of the current system
    type(TGeometry), intent(in) :: geo

    !> Contains the input for the solvation module on exit
    type(TSASAInput), intent(out) :: input

    !> Default value for the probe radius
    real(dp), intent(in), optional :: probeRadDefault

    !> Default values for the surface tension
    real(dp), intent(in), optional :: surfaceTensionDefault(:)

    type(string) :: buffer, modifier
    type(fnode), pointer :: child, value1, field, dummy
    character(lc) :: errorStr
    real(dp), allocatable :: vdwRadDefault(:)

    if (geo%tPeriodic .or. geo%tHelical) then
      call detailedError(node, "SASA model currently not available with the&
         & selected boundary conditions")
    end if

    call getChildValue(node, "ProbeRadius", input%probeRad, probeRadDefault, &
        & modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), lengthUnits, field, input%probeRad)

    call getChildValue(node, "Smoothing", input%smoothingPar, 0.3_dp*AA__Bohr, &
        & modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), lengthUnits, field, input%smoothingPar)

    call getChildValue(node, "Tolerance", input%tolerance, 1.0e-6_dp, child=child)

    call readAngularGrid(node, input%gridSize, 230)

    call readVanDerWaalsRad(node, geo, input%vdwRad)

    allocate(input%surfaceTension(geo%nSpecies))
    if (present(surfaceTensionDefault)) then
      call getChildValue(node, "SurfaceTension", value1, "Defaults", child=child)
    else
      call getChildValue(node, "SurfaceTension", value1, child=child)
    end if
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(child, "Unknown method '"//char(buffer)//"' to generate surface tension")
    case("defaults")
      if (.not.present(surfaceTensionDefault)) then
        call detailedError(child, "No defaults available for surface tension values")
      end if
      call readSpeciesList(value1, geo%speciesNames, input%surfaceTension, &
          & surfaceTensionDefault)
    case("values")
      call readSpeciesList(value1, geo%speciesNames, input%surfaceTension)
    end select

    call getChildValue(node, "Offset", input%sOffset, 2.0_dp * AA__Bohr, &
        & modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), lengthUnits, field, input%sOffset)

  end subroutine readSolvSASA


  !> Read settings for charge model 5.
  subroutine readCM5(node, input, geo)

    !> Node to process
    type(fnode), pointer :: node

    !> Geometry of the current system
    type(TGeometry), intent(in) :: geo

    !> Contains the input for the CM5 module on exit
    type(TCM5Input), intent(out) :: input

    type(fnode), pointer :: value1, dummy, child, field
    type(string) :: buffer, modifier
    real(dp) :: conv
    real(dp), allocatable :: atomicRadDefault(:)

    call getChildValue(node, "Alpha", input%alpha, 2.474_dp/AA__Bohr, &
      & modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), inverseLengthUnits, field, input%alpha)

    conv = 1.0_dp
    allocate(input%atomicRad(geo%nSpecies))
    call getChildValue(node, "Radii", value1, "AtomicRadii", child=child)
    call getChild(value1, "", dummy, modifier=modifier)
    call convertUnitHsd(char(modifier), lengthUnits, child, conv)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(child, "Unknown method '"//char(buffer)//"' to generate radii")
    case("atomicradii")
      allocate(atomicRadDefault(geo%nSpecies))
      atomicRadDefault(:) = getAtomicRad(geo%speciesNames)
      call readSpeciesList(value1, geo%speciesNames, input%atomicRad, conv=conv, &
        & default=atomicRadDefault)
      deallocate(atomicRadDefault)
    case("values")
      call readSpeciesList(value1, geo%speciesNames, input%atomicRad, conv=conv)
    end select
    if (any(input%atomicRad <= 0.0_dp)) then
      call detailedError(value1, "Atomic radii must be positive for all species")
    end if
    input%atomicRad(:) = input%atomicRad * conv

    call getChildValue(node, "Cutoff", input%rCutoff, 30.0_dp, &
        & modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), lengthUnits, field, input%rCutoff)

  end subroutine readCM5


  subroutine readSolvent(node, solvent)

    !> Node to process
    type(fnode), pointer :: node

    !> Data associated with the solvent
    type(TSolventData), intent(out) :: solvent

    type(string) :: buffer, modifier
    type(fnode), pointer :: child, value1, field
    logical :: found

    call getChildValue(node, "Solvent", value1, child=child)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(child, "Invalid solvent method '" // char(buffer) // "'")
    case('fromname')
      call getChildValue(value1, "", buffer)
      call SolventFromName(solvent, unquote(char(buffer)), found)
      if (.not. found) then
        call detailedError(value1, "Invalid solvent " // char(buffer))
      end if
    case('fromconstants')
      call getChildValue(value1, "Epsilon", buffer)
      if (unquote(char(buffer)) == "Inf") then
         if (ieee_support_inf(solvent%dielectricConstant)) then
            solvent%dielectricConstant = ieee_value(solvent%dielectricConstant, ieee_positive_inf)
         else
            solvent%dielectricConstant = huge(solvent%dielectricConstant)
         end if
      else
         call getChildValue(value1, "Epsilon", solvent%dielectricConstant)
      end if
      call getChildValue(value1, "MolecularMass", solvent%molecularMass, &
        & modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), massUnits, field, solvent%molecularMass)
      call getChildValue(value1, "Density", solvent%density, modifier=modifier, &
        & child=field)
      call convertUnitHsd(char(modifier), massDensityUnits, field, solvent%density)
    end select

  end subroutine readSolvent


  !> Reference state for free energy calculation
  subroutine readReferenceState(node, solvent, temperature, shift, freeEnergyShift)

    !> Node to process
    type(fnode), pointer :: node

    !> Data associated with the solvent
    type(TSolventData), intent(in) :: solvent

    !> Temperature for calculation
    real(dp), intent(in) :: temperature

    !> Shift to free energy
    real(dp), intent(in) :: shift

    !> Free energy shift includings state specific terms
    real(dp), intent(out) :: freeEnergyShift

    type(string) :: state
    type(fnode), pointer :: child
    real(dp), parameter :: referenceDensity = kg__au/(1.0e10_dp*AA__Bohr)**3
    real(dp), parameter :: referenceMolecularMass = amu__au
    real(dp), parameter :: idealGasMolVolume = 24.79_dp

    call getChildValue(node, "State", state, "gsolv", child=child)
    select case(tolower(unquote(char(state))))
    case default
      call detailedError(child, "Unknown reference state: "//char(state))
    case("gsolv") ! just the bare shift
      freeEnergyShift = shift
    case("reference") ! gsolv=reference option in cosmotherm
      ! RT * ln(ideal gas mol volume) + ln(rho/M)
      freeEnergyShift = shift + temperature &
          & * (log(idealGasMolVolume * temperature / ambientTemperature) &
          & + log(solvent%density/referenceDensity * referenceMolecularMass/solvent%molecularMass))
    case("mol1bar")
      ! RT * ln(ideal gas mol volume)
      freeEnergyShift = shift + temperature &
          & * log(idealGasMolVolume * temperature / ambientTemperature)
    end select
  end subroutine readReferenceState


  subroutine readVanDerWaalsRad(node, geo, vdwRad)

    !> Node to process
    type(fnode), pointer :: node

    !> Geometry of the current system
    type(TGeometry), intent(in) :: geo

    !> Van-der-Waals Radii
    real(dp), allocatable, intent(out) :: vdwRad(:)

    type(string) :: buffer, modifier
    type(fnode), pointer :: child, value1, field, dummy
    real(dp) :: conv
    real(dp), allocatable :: vdwRadDefault(:)

    conv = 1.0_dp
    allocate(vdwRad(geo%nSpecies))
    call getChildValue(node, "Radii", value1, "vanDerWaalsRadiiD3", child=child)
    call getChild(value1, "", dummy, modifier=modifier)
    call convertUnitHsd(char(modifier), lengthUnits, child, conv)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(child, "Unknown method '"//char(buffer)//"' to generate radii")
    case("vanderwaalsradiid3")
      allocate(vdwRadDefault(geo%nSpecies))
      vdwRadDefault(:) = getVanDerWaalsRadiusD3(geo%speciesNames)
      call readSpeciesList(value1, geo%speciesNames, vdwRad, vdwRadDefault, &
        & conv=conv)
      deallocate(vdwRadDefault)
    case("vanderwaalsradiicosmo")
      allocate(vdwRadDefault(geo%nSpecies))
      vdwRadDefault(:) = getVanDerWaalsRadiusCosmo(geo%speciesNames)
      call readSpeciesList(value1, geo%speciesNames, vdwRad, vdwRadDefault, &
        & conv=conv)
      deallocate(vdwRadDefault)
    case("vanderwaalsradiibondi")
      allocate(vdwRadDefault(geo%nSpecies))
      vdwRadDefault(:) = getVanDerWaalsRadiusBondi(geo%speciesNames)
      call readSpeciesList(value1, geo%speciesNames, vdwRad, vdwRadDefault, &
        & conv=conv)
      deallocate(vdwRadDefault)
    case("values")
      call readSpeciesList(value1, geo%speciesNames, vdwRad, conv=conv)
    end select
    vdwRad(:) = vdwRad * conv

  end subroutine readVanDerWaalsRad


  subroutine readAngularGrid(node, angGrid, default)

    !> Node to process
    type(fnode), pointer :: node

    !> Grid identifier
    integer, intent(out) :: angGrid

    !> Default grid size
    integer, intent(in), optional :: default

    type(fnode), pointer :: child
    character(lc) :: errorStr
    integer :: gridPoints

    call getChildValue(node, "AngularGrid", gridPoints, default, child=child)
    angGrid = 0
    call bisection(angGrid, gridSize, gridPoints)
    if (angGrid == 0) then
      call detailedError(child, "Illegal number of grid points for numerical integration")
    end if
    if (gridSize(angGrid) /= gridPoints) then
      write(errorStr, '(a, *(1x, i0, 1x, a))') &
          & "No angular integration grid with", gridPoints, &
          & "points available, using",  gridSize(angGrid), "points instead"
      call detailedWarning(child, trim(errorStr))
    end if

  end subroutine readAngularGrid


end module dftbp_solvation_solvparser
