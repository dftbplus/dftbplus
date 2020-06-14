!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *  
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program * 
!  under the LGPL licence.                                                *
!**************************************************************************
Module parameters
  
  use dftbp_accuracy, only : dp

  implicit none
  private

  public :: init_defaults
  public :: set_verbose, set_scratch, set_contdir, set_fermi, set_potentials
  public :: set_temperature, set_ncont, set_mol_indeces, set_cont_indeces
  public :: set_dopoisson, set_poissonbox, set_poissongrid, set_accuracy  
  public :: set_cluster, set_builtin

  interface
    real(dp) function DLAMCH(C)
      import :: dp
      character C
    end function DLAMCH
  end interface
  
  integer, public, parameter :: MAXNCONT=10

  integer,  public :: verbose
  integer,  public :: biasdir
  integer,  public :: gatedir
  integer,  public :: contdir(MAXNCONT)
  integer,  public :: ncont
  integer,  public :: ni(MAXNCONT)
  integer,  public :: nf(MAXNCONT)
  integer,  public :: iatc(3,MAXNCONT)
  integer,  public :: iatm(2)
  integer,  public :: ncdim(MAXNCONT)
  integer,  public :: mbound_end(MAXNCONT)
  integer,  public :: maxiter
  integer,  public :: localBC(MAXNCONT)
  integer,  public :: poissBC(MAXNCONT)
  integer,  public :: overrideBC(6)
  integer,  public :: overrBulkBC(6)
  integer,  public :: maxpoissiter
  
  real(kind=dp),  public :: Temp
  real(kind=dp),  public :: telec
  real(kind=dp),  public :: deltaR_max
  real(kind=dp),  public :: LmbMax
  real(kind=dp),  public :: gate
  real(kind=dp),  public :: GateLength_l, GateLength_t, OxLength
  real(kind=dp),  public :: Efermi(MAXNCONT)
  real(kind=dp),  public :: bias_dEf
  real(kind=dp),  public :: Rmin_Ins
  real(kind=dp),  public :: Rmin_Gate
  real(kind=dp),  public :: dr_eps
  real(kind=dp),  public :: eps_r
  real(kind=dp),  public :: cntr_gate(3)
  real(kind=dp),  public :: tip_atom(3)
  real(kind=dp),  public :: base_atom1(3),base_atom2(3)
  real(kind=dp),  public :: tipbias
  real(kind=dp),  public :: DOS(MAXNCONT)
  real(kind=dp),  public :: delta
  real(kind=dp),  public :: racc
  real(kind=dp),  public :: PoissAcc
  real(kind=dp),  public :: dmin(3)
  real(kind=dp),  public :: PoissBox(3,3)
  real(kind=dp),  public :: PoissBounds(3,2)
  real(kind=dp),  public :: PoissPlane(2)
  real(kind=dp),  public :: mu(MAXNCONT)
  real(kind=dp),  public :: cntr_cont(3,MAXNCONT)
  real(kind=dp),  public :: R_cont(MAXNCONT)
  real(kind=dp),  public :: dR_cont(MAXNCONT)
  real(kind=dp),  public :: x0
  real(kind=dp),  public :: y0
  real(kind=dp),  public :: z0
  real(kind=dp),  public :: bufferBox
  
  logical,  public :: etb 
  logical,  public :: cluster
  logical,  public :: SavePot
  logical,  public :: SaveNNList
  logical,  public :: SaveHS
  logical,  public :: FictCont(MAXNCONT)
  logical,  public :: FoundBox
  logical,  public :: DoPoisson
  logical,  public :: Readold
  logical,  public :: InitPot
  logical,  public :: DoGate
  logical,  public :: DoCilGate
  logical,  public :: DoTip
  logical,  public :: ReadBulk
  logical,  public :: mixed(6)
  !! Specify if the renormalization volume needs to be calculated
  !! Modified runtime to avoid multiple evaluation
  logical, public :: do_renorm
    !! Or you want just the "fixed" renormalization without 
    !! discretization error and cutoff compensation?
  logical, public :: fixed_renorm
  
  character(:), allocatable, public :: scratchfolder
  
  
  contains
  
  ! -----------------------------------------------------------------
  ! INIT DEFAULT VALUES 
  ! -----------------------------------------------------------------
  subroutine init_defaults()
    
    verbose=0
    biasdir=0
    gatedir=0
    contdir(:)=0
    ncont=0
    ni(:)=0; ni(1)=1;
    nf(:)=0; nf(1)=2;
    iatc(:,:)=0
    iatm(:)=0
    ncdim(:)=0
    mbound_end(:)=0
    maxiter=30    
    localBC=0
    poissBC=0
    overrideBC=0
    overrBulkBC=-1
    mixed = .false.
    maxpoissiter=60
    
    Temp=0.0_dp
    deltaR_max=9.0_dp
    LmbMax=0.50_dp
    gate=0.0_dp
    GateLength_l=0.0_dp
    GateLength_t=0.0_dp
    OxLength=0.0_dp
    Efermi(:)=0.0_dp
    bias_dEf=0.0_dp
    Rmin_Ins=0.0_dp
    Rmin_Gate=0.0_dp
    dr_eps=0.0_dp
    eps_r=0.0_dp
    cntr_gate(:)=0.0_dp
    tip_atom(:)=0.0_dp
    base_atom1(:)=0.0_dp
    base_atom2(:)=0.0_dp
    tipbias=0.0_dp
    DOS(:)=0.0_dp
    delta=1e-4_dp
    racc=DLAMCH('Precision')
    PoissAcc=1e-5_dp
    dmin(:)=0.50_dp
    PoissBox(:,:)=0.0_dp
    PoissPlane(:)=0.0_dp
    mu(:)=0.0_dp
    cntr_cont(:,:)=0.0_dp
    R_cont(:)=0.0_dp
    dR_cont(:)=1.0_dp
    bufferBox = 0.0_dp
    
    etb=.false.
    cluster=.false.
    SavePot=.false.
    SaveNNList=.false.
    SaveHS=.false.
    FictCont(:)=.false.
    DoPoisson=.true.
    Readold=.false.
    InitPot=.false.
    DoGate=.false.
    DoCilGate=.false.
    DoTip=.false.
    ReadBulk=.false.
    FoundBox=.false.
    do_renorm=.true.
    fixed_renorm=.true.

  end subroutine
 
  subroutine set_scratch(scratch)
    character(*) :: scratch      
    allocate(character(len=len(scratch))::scratchfolder)
    scratchfolder = trim(scratch)
  end subroutine set_scratch  

  subroutine set_verbose(verb)
    integer, intent(in) :: verb
    verbose = verb
  end subroutine set_verbose

  subroutine set_temperature(tt)
    real(dp), intent(in) :: tt
    temp = tt
  end subroutine set_temperature

  subroutine set_cluster(logic)
    logical, intent(in) :: logic
    cluster = logic
  end subroutine set_cluster

  subroutine set_mol_indeces(molinds,natoms)
    integer, intent(in) :: molinds(2)
    integer, intent(in) :: natoms
    if (cluster) then
      iatm(1) = 1
      iatm(2) = natoms
    else    
      iatm(1:2) = molinds(1:2)
    end if  
  end subroutine set_mol_indeces      

  subroutine set_cont_indeces(continds,ii)
    integer, intent(in) :: continds(:)
    integer, intent(in) :: ii
    integer :: nconts
    nconts = size(continds)
    iatc(ii,1:nconts) = continds(1:nconts)
  end subroutine set_cont_indeces      

  subroutine set_contdir(dir)
    integer, intent(in) :: dir(:)
    integer :: nconts
    if (cluster) then
      contdir(1) = 0
    else
      nconts = size(dir)
      contdir(1:nconts) = dir(1:nconts)
    end if  
  end subroutine set_contdir    

  subroutine set_potentials(pot)
    real(dp), intent(in) :: pot(:)
    integer :: nconts
    nconts = size(pot)
    mu(1:nconts) = pot(1:nconts)
  end subroutine set_potentials  

  subroutine set_fermi(fm)
    real(dp), intent(in) :: fm(:)
    integer :: nconts
    if (cluster) then
      EFermi(1) = 0.0_dp
    else    
      nconts = size(fm)
      EFermi(1:nconts) = fm(1:nconts)
    end if  
  end subroutine set_fermi  

  subroutine set_builtin()
    integer :: i    
    if (.not.cluster) then
      do i = 1,ncont
        mu(i) = mu(i) + EFermi(i) -  minval(EFermi(1:ncont))
      enddo
    end if
  end subroutine set_builtin      

  subroutine set_ncont(nc)
    integer, intent(in) :: nc
    ncont = nc
  end subroutine set_ncont
  
  subroutine set_dopoisson(logic)
    logical :: logic
    DoPoisson = logic
  end subroutine set_dopoisson

  subroutine set_poissonbox(Box)
    real(dp), intent(in) :: Box(:)
    integer i
    do i = 1, 3
      PoissBox(i,i)=box(i)
    end do
  end subroutine set_poissonbox

  subroutine set_poissongrid(Grid)
    real(dp), intent(in) :: Grid(:)
    dmin(1:3) = Grid(1:3)
  end subroutine set_poissongrid

  subroutine set_accuracy(Acc)
    real(dp), intent(in) :: Acc
    PoissAcc = Acc
  end subroutine set_accuracy

end module Parameters
