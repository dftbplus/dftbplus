Module parameters
  
  use gprecision
  use gconstants, only : HAR,ATU

  implicit none
  private

  public :: init_defaults
  
  interface
    real(dp) function DLAMCH(C)
      import :: dp
      character C
    end function DLAMCH
  end interface
  
integer, public, parameter :: MAXNCONT=10

  !! Renormalization volumes: to ensure charge neutrality
  real(kind=dp), public, save, allocatable, dimension(:,:) :: renorm

integer,  public, save :: verbose
integer,  public, save :: biasdir
integer,  public, save :: gatedir
integer,  public, save :: contdir(MAXNCONT)
integer,  public, save :: ncont
integer,  public, save :: ni(MAXNCONT)
integer,  public, save :: nf(MAXNCONT)
integer,  public, save :: iatc(3,MAXNCONT)
integer,  public, save :: iatm(2)
integer,  public, save :: ncdim(MAXNCONT)
integer,  public, save :: mbound_end(MAXNCONT)
integer,  public, save :: maxiter
integer,  public, save :: localBC(MAXNCONT)
integer,  public, save :: poissBC(MAXNCONT)
integer,  public, save :: overrideBC(6)
integer,  public, save :: overrBulkBC(6)
integer,  public, save :: maxpoissiter

real(kind=dp),  public, save :: hartree, a_u
real(kind=dp),  public, save :: Temp
real(kind=dp),  public, save :: telec
real(kind=dp),  public, save :: deltaR_max
real(kind=dp),  public, save :: LmbMax
real(kind=dp),  public, save :: gate
real(kind=dp),  public, save :: GateLength_l, GateLength_t, OxLength
real(kind=dp),  public, save :: Efermi(MAXNCONT)
real(kind=dp),  public, save :: bias_dEf
real(kind=dp),  public, save :: Rmin_Ins
real(kind=dp),  public, save :: Rmin_Gate
real(kind=dp),  public, save :: dr_eps
real(kind=dp),  public, save :: eps_r
real(kind=dp),  public, save :: cntr_gate(3)
real(kind=dp),  public, save :: tip_atom(3)
real(kind=dp),  public, save :: base_atom1(3),base_atom2(3)
real(kind=dp),  public, save :: tipbias
real(kind=dp),  public, save :: DOS(MAXNCONT)
real(kind=dp),  public, save :: delta
real(kind=dp),  public, save :: racc
real(kind=dp),  public, save :: PoissAcc
real(kind=dp),  public, save :: dmin(3)
real(kind=dp),  public, save :: PoissBox(3,3)
real(kind=dp),  public, save :: PoissBounds(3,2)
real(kind=dp),  public, save :: PoissPlane(2)
real(kind=dp),  public, save :: mu(MAXNCONT)
real(kind=dp),  public, save :: cntr_cont(3,MAXNCONT)
real(kind=dp),  public, save :: R_cont(MAXNCONT)
real(kind=dp),  public, save :: dR_cont(MAXNCONT)
real(kind=dp),  public, save :: x0
real(kind=dp),  public, save :: y0
real(kind=dp),  public, save :: z0
real(kind=dp),  public, save :: bufferBox

logical,  public, save :: etb 
logical,  public, save :: cluster
logical,  public, save :: SavePot
logical,  public, save :: SaveNNList
logical,  public, save :: SaveHS
logical,  public, save :: FictCont(MAXNCONT)
logical,  public, save :: FoundBox
logical,  public, save :: DoPoisson
logical,  public, save :: Readold
logical,  public, save :: InitPot
logical,  public, save :: DoGate
logical,  public, save :: DoCilGate
logical,  public, save :: DoTip
logical,  public, save :: ReadBulk
logical,  public, save :: mixed(6)
!! Specify if the renormalization volume needs to be calculated
!! Modified runtime to avoid multiple evaluation
logical, public, save :: do_renorm
!! Or you want just the "fixed" renormalization without 
!! discretization error and cutoff compensation?
logical, public, save :: fixed_renorm

interface set_global
     module procedure setglobal_d, setglobal_i, setglobal_l
end interface

contains


subroutine setglobal_i(param,value)
   integer :: param,value
   param=value
end subroutine

subroutine setglobal_d(param,value)
   real(kind=dp) :: param,value
   param=value
end subroutine

subroutine setglobal_l(param,value)
   logical :: param,value
   param=value
end subroutine

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

 hartree=HAR      ! this sets the conversion factors as set 
 a_u=ATU           ! in constants.F90
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

end module Parameters
