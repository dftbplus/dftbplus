!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> main module for the DFTB+ API
module dftbp_mainapi
  use dftbp_environment, only : TEnvironment
  use dftbp_accuracy, only : dp
  use dftbp_main, only : processGeometry
  use dftbp_initprogram, only : initProgramVariables, destructProgramVariables, coord0, latVec,&
      & tCoordsChanged, tLatticeChanged, energy, derivs, TRefExtPot, refExtPot, tExtField, orb,&
      & nAtom, nSpin, q0, qOutput, sccCalc, tExtChrg, tForces, chrgForces, qDepExtPot
  use dftbp_assert
  use dftbp_qdepextpotproxy, only : TQDepExtPotProxy
  implicit none
  private

  public :: initProgramVariables, destructProgramVariables
  public :: setGeometry, setQDepExtPotProxy, setExternalPotential, setExternalCharges, setShellResolvedCharges
  public :: getEnergy, getGradients, getExtChargeGradients, getGrossCharges
  public :: nrOfAtoms
  public :: updateDataDependentOnSpeciesOrdering, checkSpeciesNames

contains

  !> Sets up the atomic geometry
  subroutine setGeometry(coords, latVecs)

    !> atom coordinates
    real(dp), intent(in) :: coords(:,:)

    !> lattice vectors, if periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    @:ASSERT(size(coords,1) == 3)
    coord0(:,:) = coords
    tCoordsChanged = .true.
    if (present(latVecs)) then
      latVec(:,:) = latVecs
      tLatticeChanged = .true.
    else
      tLatticeChanged = .false.
    end if

  end subroutine setGeometry


  !> Returns the free energy of the system at finite temperature
  subroutine getEnergy(env, merminEnergy)

    !> Instance
    type(TEnvironment), intent(inout) :: env

    !> Resulting energy
    real(dp), intent(out) :: merminEnergy

    call recalcGeometry(env)
    merminEnergy = energy%EMermin

  end subroutine getEnergy


  !> get forces on atoms
  subroutine getGradients(env, gradients)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> resulting gradients wrt atom positions
    real(dp), intent(out) :: gradients(:,:)

    @:ASSERT(size(gradients,1) == 3)
    call recalcGeometry(env)
    gradients(:,:) = derivs

  end subroutine getGradients


  !> get the gross (Mulliken projected) charges for atoms wrt neutral atoms
  subroutine getGrossCharges(atomCharges)

    !> resulting charges
    real(dp), intent(out) :: atomCharges(:)

    atomCharges(:) = sum(q0(:, :, 1) - qOutput(:, :, 1), dim=1)

  end subroutine getGrossCharges


  !> Sets up an external population independent electrostatic potential.
  !>
  !> Sign convention: charge of electron is considered to be positive.
  !>
  subroutine setExternalPotential(atomPot, shellPot, potGrad)

    !> Atomic external potential
    real(dp), intent(in), optional :: atomPot(:)

    !> Shell resolved electrostatic potential
    real(dp), intent(in), optional :: shellPot(:,:)

    !> Gradient of the electrostatic potential
    real(dp), intent(in), optional :: potGrad(:,:)

    ! Using explicit allocation instead of F2003 automatic ones in order to stop eventual
    ! shape mismatches already at this point rather than later deep in the main code
    if (present(atomPot)) then
      if (.not. allocated(refExtPot%atomPot)) then
        allocate(refExtPot%atomPot(nAtom, nSpin))
      end if
      @:ASSERT(all(shape(atomPot) == [nAtom]))
      refExtPot%atomPot(:,1) = atomPot
    end if
    if (present(shellPot)) then
      if (.not. allocated(refExtPot%shellPot)) then
        allocate(refExtPot%shellPot(orb%mShell, nAtom, nSpin))
      end if
      @:ASSERT(all(shape(shellPot) == [orb%mShell, nAtom]))
      refExtPot%shellPot(:,:,1) = shellPot
    end if
    if (present(potGrad)) then
      if (.not. allocated(refExtPot%potGrad)) then
        allocate(refExtPot%potGrad(3, nAtom))
      end if
      @:ASSERT(all(shape(potGrad) == [3, nAtom]))
      refExtPot%potGrad(:,:) = potGrad
    end if
    tExtField = .true.

  end subroutine setExternalPotential


  !> Sets up a generator for external population dependant potentials
  subroutine setQDepExtPotProxy(extPotProxy)

    !> Generator for the external population dependant potential
    type(TQDepExtPotProxy), intent(in) :: extPotProxy

    qDepExtPot = extPotProxy

  end subroutine setQDepExtPotProxy


  !> Sets up external point charges
  subroutine setExternalCharges(chargeCoords, chargeQs, blurWidths)

    !> Coordiante of the external charges
    real(dp), intent(in) :: chargeCoords(:,:)

    !> Charges of the external point charges (sign convention: electron is negative)
    real(dp), intent(in) :: chargeQs(:)

    !> Widths of the Gaussian for each charge used for blurring (0.0 = no blurring)
    real(dp), intent(in), optional :: blurWidths(:)

    tExtChrg = .true.
    if (tForces) then
      if (.not. allocated(chrgForces)) then
        allocate(chrgForces(3, size(chargeQs)))
      end if
    end if
    call sccCalc%setExternalCharges(chargeCoords, chargeQs)

  end subroutine setExternalCharges


  !> Returns the gradient acting on the external point charges
  subroutine getExtChargeGradients(chargeGradients)

    !> Gradients
    real(dp), intent(out) :: chargeGradients(:,:)

    @:ASSERT(tForces .and. allocated(chrgForces))

    chargeGradients(:,:) = chrgForces

  end subroutine getExtChargeGradients


  !> Obtains number of atoms in the system
  function nrOfAtoms()
    integer :: nrOfAtoms

    nrOfAtoms = nAtom

  end function nrOfAtoms


  !> Check that the order of speciesName remains constant
  !> Keeping speciesNames constant avoids the need to reset:
  !> atomEigVal, referenceN0, speciesMass and SK parameters 
  function checkSpeciesNames(inputSpeciesName) result(tSpeciesNameChanged)
    use dftbp_initprogram, only: speciesName
    use dftbp_accuracy, only: mc

    !> Labels of atomic species from external program
    character(mc),  allocatable, intent(in) :: inputSpeciesName(:)
    !> Has speciesName changed?
    logical                                 :: tSpeciesNameChanged
    
    tSpeciesNameChanged = any(speciesName /= inputSpeciesName)
 
  end function checkSpeciesNames

  
  !> When order of atoms changes, update arrays containing atom type indices,
  !> and all subsequent dependencies. Updated data returned via use statements
  subroutine updateDataDependentOnSpeciesOrdering(env, inputSpecies)
    !TODO(Alex) These will go to the top of the module
    use dftbp_environment, only: TEnvironment
    use dftbp_accuracy,    only: mc
    use dftbp_initprogram, only: species0, mass, orb, iEqOrbitals, nIneqOrb, nMixElements, &
         denseDesc, HSqrCplx, SSqrCplx, eigVecsCplx, HSqrReal, SSqrReal, eigVecsReal, &
         q0, qShell0, qInput, qOutput, qInpRed, qOutRed, qDiffRed, nEl, speciesName

    !> dftb+ environment 
    type(TEnvironment),   intent(in) :: env
    !> types of the atoms (nAllAtom) 
    integer, allocatable, intent(in) :: inputSpecies(:)
    
    !Number of atoms must be conserved    
    @:ASSERT(size(inputSpecies) == nAtom) 
    species0 = inputSpecies
    mass =  updateAtomicMasses(species0)
    orb%nOrbAtom = updateAtomicOrbitals(species0)
    call updateEquivalencyRelations(species0, iEqOrbitals, nIneqOrb, nMixElements)
    call updateBLACSDecomposition(env, denseDesc)
    call reallocateHSArrays(env, denseDesc, HSqrCplx, SSqrCplx, eigVecsCplx, HSqrReal, SSqrReal, eigVecsReal)
    !If atomic order changes, partial charges need to be initialised, else wrong charge will be associated
    !with each atom
    call initialiseCharges(species0, speciesName, orb, nEl, iEqOrbitals, nIneqOrb, nMixElements,&
         q0, qShell0, qInput, qInpRed, qOutput, qOutRed, qDiffRed)
    
  end subroutine updateDataDependentOnSpeciesOrdering

  
  !> Set shell-resolved partial charges.
  !> If no argument is given, the final charges from the prior calculation are used 
  subroutine setShellResolvedCharges(qSeed)
    use dftbp_initprogram, only: qInput, qOutput, qInpRed, qOutRed, iEqOrbitals, nIneqOrb, orb 
    use dftbp_spin,        only: qm2ud,ud2qm
    use dftbp_orbitalequiv,only: OrbitalEquiv_reduce
    
    real(dp), allocatable, intent(in), optional :: qSeed(:, :, :)
    
    if(present(qSeed)) then
       @:ASSERT(size(qSeed,1) == size(qInput,1))
       @:ASSERT(size(qSeed,2) == size(qInput,2))
       @:ASSERT(size(qSeed,3) == size(qInput,3))
       qInput = qSeed
       !Fill qInpRed
       if (nSpin == 2) call qm2ud(qInput)
       call OrbitalEquiv_reduce(qInput, iEqOrbitals, orb, qInpRed(1:nIneqOrb))
       if (nSpin == 2) call ud2qm(qInput)
    else
       qInput  = qOutput
       qInpRed = qOutRed 
    endif
    
    qOutput = 0._dp
    qOutRed = 0._dp
 
  end subroutine setShellResolvedCharges

   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  !> Update order of nr. atomic orbitals for each atom, orb%nOrbAtom
  function updateAtomicOrbitals(species0) result(nOrbAtomReordered)
    use dftbp_initprogram, only: orb, nAtom
    !> Type of the atoms (nAtom)
    integer,         allocatable,  intent(in) :: species0(:)
    !> Nr. of orbitals for each atom (nAtom)  
    integer,         allocatable :: nOrbAtomReordered(:)

    @:ASSERT(nAtom == size(species0))
    allocate(nOrbAtomReordered(nAtom))
    nOrbAtomReordered(:) = orb%nOrbSpecies(species0(:)) 
  end function updateAtomicOrbitals

  
  !> Update atomic masses
  function updateAtomicMasses(species0) result(massReordered)
    use dftbp_initprogram, only: speciesMass, nAtom
    !> Type of the atoms (nAtom)
    integer,  allocatable,   intent(in) :: species0(:)
    !> List of atomic masses (nAtom)
    real(dp), allocatable :: massReordered(:)

    @:ASSERT(nAtom == size(speciesMass))
    allocate(massReordered(size(speciesMass)))
    massReordered = speciesMass(species0)
  end function updateAtomicMasses

  
  ! TODO(Alex) This contains copy and pasted code from initprogram.F90
  ! Could be adapted into a function that is also called in initprogram.F90
  !
  !> Update equivalency relations
  !> Required for partial charge initialisation
  subroutine updateEquivalencyRelations(species0, iEqOrbitals, nIneqOrb, nMixElements)
    use dftbp_initprogram,        only: tSccCalc, orb ,nAtom, nSpin, sccCalc, tDFTBU, onSiteElements
    use dftbp_spin,               only: Spin_getOrbitalEquiv
    use dftbp_orbitalequiv,       only: OrbitalEquiv_merge
    use dftbp_message,            only: error
    
    !> Type of the atoms (nAtom)
    integer,   allocatable, intent(in) :: species0(:)
    !> Orbital equivalence relations
    integer,   allocatable, intent(inout) :: iEqOrbitals(:,:,:)
    !> nr. of inequivalent orbitals
    integer,   intent(inout) :: nIneqOrb
    !> nr. of elements to go through the mixer - may contain reduced orbitals and also orbital blocks
    !> (if tDFTBU or onsite corrections)
    integer,   intent(inout) :: nMixElements
    ! Local data 
    integer,   allocatable  :: iEqOrbSCC(:,:,:), iEqOrbSpin(:,:,:)

    integer :: i,j
    
    if (tSccCalc) then
       if(.not. allocated(iEqOrbitals)) allocate(iEqOrbitals(orb%mOrb, nAtom, nSpin))
       if(.not. allocated(iEqOrbSCC))   allocate(iEqOrbSCC(orb%mOrb, nAtom, nSpin))
       call sccCalc%getOrbitalEquiv(orb, species0, iEqOrbSCC)
       
       if (nSpin == 1) then
          iEqOrbitals(:,:,:) = iEqOrbSCC(:,:,:)
       else
          allocate(iEqOrbSpin(orb%mOrb, nAtom, nSpin))
          call Spin_getOrbitalEquiv(orb, species0, iEqOrbSpin)
          call OrbitalEquiv_merge(iEqOrbSCC, iEqOrbSpin, orb, iEqOrbitals)
          deallocate(iEqOrbSpin)
       end if

       deallocate(iEqOrbSCC)
       nIneqOrb = maxval(iEqOrbitals)
       nMixElements = nIneqOrb

       if(tDFTBU) then
          call error("DFTB+U is not supported by DFTB+ API")
       elseif(allocated(onSiteElements)) then
          call error("User-specified onsite elements is not supported by DFTB+ API")
       endif
       
    elseif(.not. tSccCalc)then
       nIneqOrb = orb%nOrb
       nMixElements = 0
    end if

  end subroutine updateEquivalencyRelations


  !> Initialise partial charges
  subroutine initialiseCharges(species0, speciesName, orb, nElectrons, iEqOrbitals, nIneqOrb, nMixElements, &
       q0, qShell0, qInput, qInpRed, qOutput, qOutRed, qDiffRed, initialSpins, initialCharges)    

    use dftbp_initprogram, only: tMulliken, tSccCalc, tSpin, tReadChrg, tCustomOccAtoms, tMixBlockCharges, tRangeSep, &
                                 referenceN0, nSpin, nAtom  !, tImHam
    use dftbp_sccinit,     only: initQFromShellChrg, initQFromAtomChrg
    use dftbp_orbitalequiv,only: OrbitalEquiv_reduce
    use dftbp_spin,        only: qm2ud,ud2qm
    use dftbp_orbitals,    only: TOrbitals
    use dftbp_accuracy,    only: mc
    use dftbp_message,     only: error
    
    !> Type of the atoms (nAtom)
    integer,  allocatable, intent(in) :: species0(:)
    character(mc),  allocatable, intent(in) :: speciesName(:)
    !> Data type for atomic orbitals
    type(TOrbitals),       intent(in) :: orb
    !> Number of electrons
    real(dp),              intent(in) :: nElectrons(:)
    !> Orbital equivalence relations
    integer,  allocatable, intent(in) :: iEqOrbitals(:,:,:)
    !> nr. of inequivalent orbitals
    integer,               intent(in) :: nIneqOrb
    !> nr. of elements to go through the mixer - may contain reduced orbitals and also orbital blocks                
    !> (if tDFTBU or onsite corrections)                                           
    integer, intent(in) :: nMixElements

    !> Initial spins
    real(dp), allocatable, intent(in), optional :: initialSpins(:,:)
    !> Set of atom-resolved atomic charges 
    real(dp), allocatable, intent(in), optional :: initialCharges(:)
 
    !> reference neutral atomic occupations
    real(dp), allocatable, intent(inout) :: q0(:, :, :)
    !> shell resolved neutral reference
    real(dp), allocatable, intent(inout) :: qShell0(:,:)
    !> input charges (for potentials)
    real(dp), allocatable, intent(inout) :: qInput(:, :, :)
    !> output charges
    real(dp), allocatable, intent(inout) :: qOutput(:, :, :)
    !> input charges packed into unique equivalence elements
    real(dp), allocatable, intent(inout) :: qInpRed(:)
    !> output charges packed into unique equivalence elements
    real(dp), allocatable, intent(inout) :: qOutRed(:)
    !> charge differences packed into unique equivalence elements
    real(dp), allocatable, intent(inout) :: qDiffRed(:)

    !> TODO(Alex) Not global in initprogram. Currently hard-coded
    !> Is check-sum for charges read externally to be used?
    logical  :: tSkipChrgChecksum = .false.

    !Local variables 
    integer        :: iAt,iSp,iSh,ii,jj,   i,j

    if (.not. tMulliken .or. .not. tSccCalc) then
       return
    endif

    !Support to be added 
    if(tReadChrg) then
       call error("Reading charges from file is not supported by DFTB+ API")
    endif
    
    if(tMixBlockCharges) then
       call error("Mixing block charges is not supported by DFTB+ API")
    endif

    if(tRangeSep) then
       call error("Range-separated charges is not supported by DFTB+ API")
    endif

    if(tCustomOccAtoms) then
       call error("Custom occupation is not supported by DFTB+ API")
    endif

    @:ASSERT(size(species0) == nAtom) 
    !call error("Number of atoms in species array differ from current state of DFTB+")

    ! Initialise reference neutral atoms 
    if(.not. allocated(q0)) allocate(q0(orb%mOrb, nAtom, nSpin))
    q0(:,:,:) = 0.0_dp
    !referenceN0 consistent as long as speciesName remains constant
    call initQFromShellChrg(q0, referenceN0, species0, orb)

    !TODO(Alex) This does not appear to get used in any instance in DFTB+. Remove?
    ! Shell resolved neutral reference
    if(.not. allocated(qShell0)) allocate(qShell0(orb%mShell, nAtom))
    qShell0(:,:) = 0.0_dp

    do iAt = 1, nAtom
       iSp = species0(iAt)
       do iSh = 1, orb%nShell(iSp)
          qShell0(iSh,iAt) = sum(q0(orb%posShell(iSh,iSp):orb%posShell(iSh+1,iSp)-1,iAt,1))
       end do
    end do
    
    !Input charges for potentials 
    if(.not. allocated(qInput))  allocate(qInput(orb%mOrb, nAtom, nSpin))
    qInput(:,:,:) = 0.0_dp
    
    if (present(initialCharges) .and. allocated(initialCharges)) then
       call initQFromAtomChrg(qInput, initialCharges, referenceN0, species0, &
            & speciesName, orb)
    else
       qInput(:,:,:) = q0
    endif
    
    !Rescaling to ensure correct number of electrons in the system
    if (.not. tSkipChrgChecksum) then
       qInput(:,:,1) = qInput(:,:,1) *  sum(nElectrons) / sum(qInput(:,:,1))
    endif
    
    !Fill additional spin channels 
    select case (nSpin)
    case (1)
       write(*,*) 'case1'
       continue 
     
    case (2)
       if (present(initialSpins) .and. allocated(initialSpins)) then   
          do ii = 1, nAtom
             qInput(1:orb%nOrbAtom(ii),ii,2) = qInput(1:orb%nOrbAtom(ii),ii,1)&
                  * initialSpins(1,ii) / sum(qInput(1:orb%nOrbAtom(ii),ii,1))
          end do
       elseif(.not. tSkipChrgChecksum) then
          do ii = 1, nAtom
             qInput(1:orb%nOrbAtom(ii),ii,2) = qInput(1:orb%nOrbAtom(ii),ii,1)&
                  * (nElectrons(1)-nElectrons(2))/sum(qInput(:,:,1))
          end do
       end if
       
    case (4)
       if (tSpin) then
          if ((.not. present(initialSpins)) .or. (.not. allocated(initialSpins))) then
             call error("Missing initial spins!")
             
          elseif (any(shape(initialSpins)/=(/3,nAtom/))) then
             call error("Incorrect shape initialSpins array!")

          elseif (.not. tSkipChrgChecksum) then
             do ii = 1, nAtom
                do jj = 1, 3
                   qInput(1:orb%nOrbAtom(ii),ii,jj+1) = qInput(1:orb%nOrbAtom(ii),ii,1)&
                        & * initialSpins(jj,ii) / sum(qInput(1:orb%nOrbAtom(ii),ii,1))
                end do
             end do
          end if
       endif
    end select

    !Input charges packed into unique equivalence elements
    if(.not. allocated(qInpRed))   allocate(qInpRed(nMixElements))
    qInpRed = 0.0_dp
  
    !Swap from charge/magnetisation to up/down 
    if (nSpin == 2) call qm2ud(qInput)
    call OrbitalEquiv_reduce(qInput, iEqOrbitals, orb, qInpRed(1:nIneqOrb))
    !Converts up/down set back to charge/magnetization                                                              
    if (nSpin == 2) call ud2qm(qInput)
    
    !Only reinitialised, not assigned specific values
    if(.not. allocated(qOutput))   allocate(qOutput(orb%mOrb, nAtom, nSpin))
    if(.not. allocated(qOutRed))   allocate(qOutRed(nMixElements))
    if(.not. allocated(qDiffRed))  allocate(qDiffRed(nMixElements))
    qOutput  = 0._dp
    qOutRed  = 0.0_dp
    qDiffRed = 0.0_dp

  end subroutine initialiseCharges

  
  ! TODO(Alex) Not sure this actually requires calling as long as the total
  ! number of orbitals and processes remains constant - CHECK
  ! Will certainly change if nAtom changes 
  !
  !> Update dense matrix descriptor for H and S in BLACS decomposition
  subroutine updateBLACSDecomposition(env, denseDesc)
    use dftbp_initprogram, only: orb, nAtom, t2Component, getDenseDescCommon, getDenseDescBlacs
    use dftbp_environment, only: TEnvironment
    use dftbp_densedescr,  only: TDenseDescr

    !> Environment settings
    type(TEnvironment), intent(in)    :: env
    !> Dense matrix descriptor for H and S
    type(TDenseDescr),  intent(inout) :: denseDesc

    ! Specificaly, denseDesc uses orb%nOrbAtom
    call getDenseDescCommon(orb, nAtom, t2Component, denseDesc)
  #:if WITH_SCALAPACK
    call getDenseDescBlacs(env, env%blacs%rowBlockSize, env%blacs%columnBlockSize, denseDesc)
  #:endif
    
  end subroutine updateBLACSDecomposition


  ! TODO(Alex) Check if these quantities change per MD step 
  ! My assumption is that if the total number of orbitals and the number of MPI processes is
  ! fixed, then nLocalRows and nLocalCols will also stay constant
  ! Hence will change if nAtom changes 
  !
  !> Reassign Hamiltonian, overlap and eigenvector arrays
  subroutine reallocateHSArrays(env, denseDesc, HSqrCplx, SSqrCplx, eigVecsCplx, HSqrReal, SSqrReal ,eigVecsReal)
    use dftbp_initprogram, only: t2Component, tRealHS, parallelKS
    use dftbp_scalapackfx, only: scalafx_getlocalshape
    use dftbp_environment, only: TEnvironment
    use dftbp_densedescr,  only: TDenseDescr
    
    !> Environment instance
    type(TEnvironment), intent(in) :: env
    !> Dense matrix descriptor for H and S
    type(TDenseDescr),  intent(in) :: denseDesc 
    !> Square dense hamiltonian storage for cases with k-points
    complex(dp), allocatable, intent(inout) :: HSqrCplx(:,:)
    !> Square dense overlap storage for cases with k-points
    complex(dp), allocatable, intent(inout)  :: SSqrCplx(:,:)
    !> Complex eigenvectors
    complex(dp), allocatable :: eigvecsCplx(:,:,:)
    !> Square dense hamiltonian storage
    real(dp),    allocatable :: HSqrReal(:,:)
    !> Square dense overlap storage
    real(dp),    allocatable :: SSqrReal(:,:)
    !> Real eigenvectors
    real(dp),    allocatable :: eigvecsReal(:,:,:)
    ! Local variables 
    integer                  :: nLocalRows, nLocalCols, nLocalKS

    !Retrieved from index array for spin and k-point index
    nLocalKS = size(parallelKS%localKS, dim=2)
    
    !Get nLocalRows and nLocalCols 
    call scalafx_getlocalshape(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, nLocalRows, nLocalCols)
  
    if (t2Component .or. .not. tRealHS) then
       if( (size(HSqrCplx,1) /= nLocalRows) .or. (size(HSqrCplx,2) /= nLocalCols) )then
          write(*,*) 'Complex arrays reallocated - see mainapi.F90' 
          deallocate(HSqrCplx)
          deallocate(SSqrCplx)
          deallocate(eigVecsCplx)
          allocate(HSqrCplx(nLocalRows, nLocalCols))
          allocate(SSqrCplx(nLocalRows, nLocalCols))
          allocate(eigVecsCplx(nLocalRows, nLocalCols, nLocalKS))
       endif
       
    else
       if( (size(HSqrReal,1) /= nLocalRows) .or. (size(HSqrReal,2) /= nLocalCols) )then
          write(*,*) 'Real arrays reallocated - see mainapi.F90' 
          deallocate(HSqrReal)
          deallocate(SSqrReal)
          deallocate(eigVecsReal)
          allocate(HSqrReal(nLocalRows, nLocalCols))
          allocate(SSqrReal(nLocalRows, nLocalCols))
          allocate(eigVecsReal(nLocalRows, nLocalCols, nLocalKS))
       endif
       
    end if
    
  end subroutine reallocateHSArrays


  !> re-evaluate the energy/forces if the geometry changes
  subroutine recalcGeometry(env)

    !> instance
    type(TEnvironment), intent(inout) :: env

    logical :: tStopDriver, tStopScc, tExitGeoOpt

    if (tLatticeChanged .or. tCoordsChanged) then
      call processGeometry(env, 1, 1, .false., tStopDriver, tStopScc, tExitGeoOpt)
      tLatticeChanged = .false.
      tCoordsChanged = .false.
    end if

  end subroutine recalcGeometry


end module dftbp_mainapi
