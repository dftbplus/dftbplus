module onscorrection
# include "assert.h"
  use Accuracy
  use CommonTypes
  use Message
  use nonscc, only : H0Sprime
  use SlakoCont  
  implicit none
  public
  
contains

  subroutine addOnsShift(potential,qBlock,orb,ons_en,species)
    real(dp), intent(inout) :: potential(:,:,:,:)
    real(dp), intent(in) :: qBlock(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: ons_en(:,:)
    integer, intent(in) :: species(:)
    
    integer :: nOrb, nAtom, nSpin
    integer :: iAt, iSp, iSpin, sigma, ud 
    integer :: mu, nu
    real(dp) :: ons_conv(orb%mOrb,orb%mOrb,2)
    real(dp) :: factor
!    real(dp) :: sfact
    
    nAtom = size(potential, dim=3)
    nSpin = size(potential, dim=4)
    
!    sfact = 0.4_dp
    
    do iSpin = 1,nSpin
      do iAt = 1, nAtom
        nOrb = orb%nOrbAtom(iAt)
        iSp  = species(iAt)
        call getOnsME(orb,iSp,ons_en,ons_conv)
        do mu = 1, nOrb-1
          do nu = mu +1, nOrb
            if (iSpin == 1) then
              factor = 1.0_dp
            else
              factor = -1.0_dp
            end if
            potential(mu,nu,iAt,iSpin) = potential(mu,nu,iAt,iSpin)&
                & + qBlock(mu,nu,iAt,iSpin)&
                &*( ons_conv(mu,nu,1) + factor*ons_conv(mu,nu,2) )

            potential(nu,mu,iAt,iSpin) = potential(mu,nu,iAt,iSpin)
          end do
        end do
        
      end do
    end do  
        
  end subroutine addOnsShift
 
  subroutine add2centShift(potential,qBlock,q,q0,orb,ons_en,species,&
      &nNeighbor,iNeighbor,over,iPair,hubbU,W)
    real(dp), intent(inout)        :: potential(:,:,:,:)
    real(dp), intent(in)           :: qBlock(:,:,:,:)
    real(dp), intent(in)           :: q(:,:,:), q0(:,:,:) 
    type(TOrbitals), intent(in)    :: orb
    real(dp), intent(in)           :: ons_en(:,:)
    integer, intent(in)            :: species(:)    
    integer, intent(in)            :: nNeighbor(:)
    integer, intent(in)            :: iNeighbor(0:,:)
    real(dp), intent(in)           :: over(:)
    integer, intent(in)            :: iPair(0:,:) 
    real(dp), intent(in)           :: hubbU(:,:)
    real(dp), intent(in), optional :: W(:,:,:)
    
    integer :: nOrb1, nOrb2, nAtom, nSpin
    integer :: iAt1, iAt2, iSpin, iNeigh, iOrig 
    integer :: ii, mu1, nu1, mu2, nu2
    real(dp) :: factor,factor1, factor2 
    real(dp) :: tmpS(orb%mOrb,orb%mOrb)
    real(dp) :: twoCInt(2)
    real(dp) :: deltaP( orb%mOrb,orb%mOrb,&
        &size(qBlock, dim=3),size(qBlock, dim=4) )
    logical :: dcount
    
    nAtom = size(qBlock, dim=3)
    nSpin = size(qBlock, dim=4)
    
    deltaP = qBlock    
    do ii =1, orb%mOrb
      deltaP(ii,ii,:,:) = q(ii,:,:) - q0(ii,:,:)
    end do
    
    do iSpin = 1, nSpin
      if (iSpin == 1) then
        factor = 1.0_dp
      else
        factor = -1.0_dp
      end if
      do iAt1 = 1, nAtom
        nOrb1 = orb%nOrbAtom(iAt1)
        do mu1 = 1, nOrb1
          do nu1 = mu1, nOrb1
            if (mu1 == nu1) then
              factor1 = 0.5_dp
            else
              factor1 = 1.0_dp
            end if
            
            do iNeigh = 1, nNeighbor(iAt1)
              iAt2 = iNeighbor(iNeigh, iAt1)
              nOrb2 = orb%nOrbAtom(iAt2)
              iOrig = iPair(iNeigh,iAt1)
              tmpS(1:nOrb2,1:nOrb1) = reshape( &
              & over(iOrig+1:iOrig+nOrb1*nOrb2),(/nOrb2,nOrb1/) )
              do mu2 = 1, nOrb2
                do nu2 = mu2, nOrb2
                  if (mu2 == nu2) then
                    factor2 = 0.5_dp
                  else
                    factor2 = 1.0_dp
                  end if
                  !avoiding double counting with SCC terms
                  dcount = ( (mu1 == nu1).and.(mu2 == nu2) )
                  if (.not. dcount) then 
                    if (present(W)) then
                      call get2CentInt(twoCInt,iAt1,mu1,nu1,iAt2,mu2,nu2,orb,&
                          &ons_en,species,tmpS,hubbU,W)
                    else
                      call get2CentInt(twoCInt,iAt1,mu1,nu1,iAt2,mu2,nu2,orb,&
                          &ons_en,species,tmpS,hubbU) 
                    end if      
                    
                    potential(mu1,nu1,iAt1,iSpin) = potential(mu1,nu1,iAt1,iSpin)&
                        & + factor2*deltaP(mu2,nu2,iAt2,iSpin)&
                        &*( twoCInt(1) + factor*twoCInt(2) )
                    potential(mu2,nu2,iAt2,iSpin) = potential(mu2,nu2,iAt2,iSpin)&
                        & + factor1*deltaP(mu1,nu1,iAt1,iSpin)&
                        &*( twoCInt(1) + factor*twoCInt(2) )  
                        
                  end if      
                end do
              end do
            end do
            potential(nu1,mu1,iAt1,iSpin) = potential(mu1,nu1,iAt1,iSpin)
          end do
        end do
      end do
    end do  
    
  end subroutine add2centShift 
  
  subroutine addRIshift(potential,q,q0,orb,ons_en,species)
    real(dp), intent(inout)        :: potential(:,:,:,:)
    real(dp), intent(in)           :: q(:,:,:), q0(:,:,:) 
    type(TOrbitals), intent(in)    :: orb
    real(dp), intent(in)           :: ons_en(:,:)
    integer, intent(in)            :: species(:)     
  
    integer :: nOrb, nAtom, nSpin
    integer :: iAt, iSpin, iSp, iOrb
    real(dp):: qDiff( orb%mOrb,&
        &size(potential, dim=3),size(potential, dim=4))
    real(dp):: factor, onethird
    real(dp) :: sfact
    
    nAtom = size(potential, dim=3)
    nSpin = size(potential, dim=4)    
    
    onethird = 1.0_dp/3.0_dp
    
    sfact = 0.4_dp
    
    do iOrb =1,orb%mOrb
      qDiff(iOrb,:,:) = q0(2,:,:) + q0(3,:,:) + q0(4,:,:) &
           &        - q(2,:,:) -  q(3,:,:) -  q(4,:,:)
    end do
    
    do iOrb =2,4
      qDiff(iOrb,:,:) = qDiff(iOrb,:,:) + &
           & 3.0_dp*( q(iOrb,:,:) - q0(iOrb,:,:) )
    end do
    
    do iSpin = 1,nSpin
      if (iSpin == 1) then
        factor = 1.0_dp
      else
        factor = -1.0_dp
      end if       
      do iAt = 1, nAtom
        nOrb = orb%nOrbAtom(iAt)
        iSp  = species(iAt)
        if (nOrb > 1) then
          do iOrb = 1,nOrb 
             potential(iOrb,iOrb,iAt,iSpin) = potential(iOrb,iOrb,iAt,iSpin)&
                 & + onethird*qDiff(iOrb,iAt,iSpin)&
                 & *( ons_en(iSp,3) + factor*ons_en(iSp,4) )  
          end do
        end if  
      end do
    end do     
  
  end subroutine addRIshift
  
  subroutine getEons(Eons,qBlock,q,q0,orb,ons_en,species)
    real(dp), intent(out)       :: Eons(:)
    real(dp), intent(in)        :: q(:,:,:), q0(:,:,:)
    real(dp), intent(in)        :: qBlock(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in)        :: ons_en(:,:)
    integer, intent(in)         :: species(:)
    
    integer :: nOrb, nAtom, nSpin
    integer :: iAt, iSp, iSpin
    integer :: mu, nu
    real(dp) :: ons_conv(orb%mOrb,orb%mOrb,2)
    real(dp) :: factor
    real(dp) :: Eri(size(Eons))
    
    nAtom = size(qBlock, dim=3)
    nSpin = size(qBlock, dim=4)
    
    Eons = 0.0_dp
    
    do iSpin = 1,nSpin
      do iAt = 1, nAtom
        nOrb = orb%nOrbAtom(iAt)
        if (nOrb == 1) cycle
        iSp  = species(iAt)
        call getOnsME(orb,iSp,ons_en,ons_conv)
        do mu = 1, nOrb-1
          do nu = mu +1, nOrb
            if (iSpin == 1) then
              factor = 1.0_dp
            else
              factor = -1.0_dp
            end if
            Eons(iAt) = Eons(iAt) + qBlock(mu,nu,iAt,iSpin)**2&
               &*( ons_conv(mu,nu,1) + factor*ons_conv(mu,nu,2) )

          end do
        end do
        
      end do
    end do
    
    call getEri(Eri,q,q0,orb,ons_en,species)
    Eons = Eons + Eri    
  
  end subroutine getEons
  
  subroutine getEri(Eri,q,q0,orb,ons_en,species)
    real(dp), intent(out)          :: Eri(:)
    real(dp), intent(in)           :: q(:,:,:), q0(:,:,:) 
    type(TOrbitals), intent(in)    :: orb
    real(dp), intent(in)           :: ons_en(:,:)
    integer, intent(in)            :: species(:)
    
    integer  :: nOrb, nAtom, nSpin
    integer  :: iAt, iSp, iSpin
    integer  :: mu, nu
    real(dp) :: factor, fact, onethird, onefifth  
    real(dp) :: qDiff(orb%mOrb,size(q, dim=2),size(q, dim=3))
        
    nAtom = size(q, dim=2)
    nSpin = size(q, dim=3)
    
    onethird = 1.0_dp/3.0_dp
    onefifth = 1.0_dp/5.0_dp
    
    Eri = 0.0_dp   
    qDiff = q - q0
    
    do iSpin = 1,nSpin
      do iAt = 1, nAtom
        nOrb = orb%nOrbAtom(iAt)
        iSp  = species(iAt)       
        if (nOrb > 1) then
          do mu = 2, 4
            do nu = mu, 4
              fact = -1.0_dp
              factor = -1.0_dp
              if (iSpin == 1) factor = 1.0_dp
              if (mu == nu) fact = 1.0_dp
              
              Eri(iAt) = Eri(iAt) + fact*onethird &
                 & *qDiff(mu,iAt,iSpin)*qDiff(nu,iAt,iSpin)&
                 & *( ons_en(iSp,3) + factor*ons_en(iSp,4) )

            end do
          end do
        end if
        !d-orbital (will be optimized)
        if (nOrb > 4) then
          do mu = 5, nOrb
            do nu = mu, nOrb
              fact = -1.0_dp
              factor = -1.0_dp
              if (iSpin == 1) factor = 1.0_dp
              if (mu == nu) fact = 2.0_dp
              
              Eri(iAt) = Eri(iAt) + fact*onefifth &
                 & *qDiff(mu,iAt,iSpin)*qDiff(nu,iAt,iSpin)&
                 & *( ons_en(iSp,9) + factor*ons_en(iSp,10) )

            end do
          end do
        end if
      end do
    end do
  
  end subroutine getEri
    
  subroutine getOnsME(orb,iSp,ons_en,ons_conv)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iSp
    real(dp), intent(in) :: ons_en(:,:)
    real(dp), intent(out) :: ons_conv(:,:,:)
    
    integer :: iOrb1, iOrb2, nOrb, ud
       
    ons_conv(:,:,:) = 0.0_dp
    
    nOrb = orb%nOrbSpecie(iSp)
    if (nOrb > 9) then
      call error("Linresp: onsite correction does not work&
          & for atoms containing 'f' orbitals yet.")
    end if
      
    do iOrb1 =1, nOrb-1
      do iOrb2 = iOrb1+1, nOrb
        do ud = 1,2
        
          if (iOrb1 == 1) then        
            ons_conv(iOrb1,iOrb2,ud) = ons_en(iSp,ud)
            if (iOrb2 > iOrb1 + 3) then
              ons_conv(iOrb1,iOrb2,ud)= ons_en(iSp,4+ud)
            endif
          !__________d orbital_______________
          elseif (iOrb1 > 4) then              
            ons_conv(iOrb1,iOrb2,ud) = ons_en(iSp,8+ud)
          !__________________________________
          else
            ons_conv(iOrb1,iOrb2,ud) = ons_en(iSp,2+ud)
            !for d orbitals:  
            if (iOrb2 > 4) then
              ons_conv(iOrb1,iOrb2,ud) = ons_en(iSp,6+ud)
            end if
          end if
          
          ons_conv(iOrb2,iOrb1,ud) = ons_conv(iOrb1,iOrb2,ud)
          
        end do  
      end do
    end do 
    
  end subroutine getOnsME
  
  subroutine Ons_blockIndx(iEqBlockOns,count,orb)
    integer, intent(out)        :: iEqBlockOns(:,:,:,:)
    integer, intent(in)         :: count
    type(TOrbitals), intent(in) :: orb

    integer :: nAtom, nSpin, nOrb, iCount
    integer :: iAt, iSp, iOrb1, iOrb2
    
    nAtom = size(iEqBlockOns, dim=3)
    nSpin = size(iEqBlockOns, dim=4)
    ASSERT(nSpin == 1 .or. nSpin == 2)

    ASSERT(size(iEqBlockOns, dim=1) == orb%mOrb)
    ASSERT(size(iEqBlockOns, dim=2) == orb%mOrb)
    
    iEqBlockOns(:,:,:,:) = 0

    iCount = count
    do iSp = 1, nSpin
      do iAt = 1, nAtom
        nOrb = orb%nOrbAtom(iAt)
        do iOrb1 = 1, nOrb - 1
	  do iOrb2 = iOrb1 + 1, nOrb
	    iCount = iCount + 1
	    iEqBlockOns(iOrb1, iOrb2, iAt, iSp) = iCount
	  end do
	end do
      end do
    end do
    
  end subroutine Ons_blockIndx
  
  subroutine Onsblock_reduce(input, equiv, orb, output)
    real(dp), intent(in) :: input(:,:,:,:)
    integer, intent(in) :: equiv(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(inout) :: output(:)

    integer :: nAtom, nSpin, nOrb
    integer :: iS, iOrb1, iOrb2, iAt
    integer :: pos
    
    nAtom = size(input, dim=3)
    nSpin = size(input, dim=4)
    ASSERT(size(input, dim=1) == orb%mOrb)
    ASSERT(size(input, dim=2) == orb%mOrb)
    ASSERT(all(shape(equiv) == (/ orb%mOrb, orb%mOrb, nAtom, nSpin /)))
    ASSERT(nSpin == 1 .or. nSpin == 2)
    
    do iS = 1, nSpin
      do iAt = 1, nAtom
        nOrb = orb%nOrbAtom(iAt)
        do iOrb1 = 1, nOrb
          do iOrb2 = 1, nOrb
            pos = equiv(iOrb1, iOrb2, iAt, iS)
            if (pos > 0) then
              output(pos) = 0.5_dp*( input(iOrb1, iOrb2, iAt, iS) &
                                &  + input(iOrb2, iOrb1, iAt, iS) )
            end if
          end do
        end do
      end do
    end do
    
  end subroutine OnsBlock_reduce
  
  subroutine Onsblock_expand(input,equiv,orb,output)
    real(dp), intent(in) :: input(:)
    integer, intent(in) :: equiv(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(inout) :: output(:,:,:,:)

    integer :: nAtom, nSpin, nOrb
    integer :: iS, iOrb1, iOrb2, iAt
    logical :: mask(size(input))
    integer :: pos
    
    nAtom = size(output, dim=3)
    nSpin = size(output, dim=4)
    ASSERT(size(output, dim=1) == orb%mOrb)
    ASSERT(size(output, dim=2) == orb%mOrb)
    ASSERT(all(shape(equiv) == (/ orb%mOrb, orb%mOrb, nAtom, nSpin /)))
    ASSERT(nSpin == 1 .or. nSpin == 2)
    
    mask(:) = .true.
    do iS = 1, nSpin
      do iAt = 1, nAtom
        nOrb = orb%nOrbAtom(iAt)
        do iOrb1 = 1, nOrb
          do iOrb2 = 1, nOrb
            pos = equiv(iOrb1, iOrb2, iAt, iS)
            if (pos > 0) then
              if (mask(pos)) then
                output(iOrb1, iOrb2, iAt, iS) = input(pos)
                output(iOrb2, iOrb1, iAt, iS) = input(pos)
                mask(pos) = .false.
              end if
            end if
          end do
        end do
      end do
    end do
    
  end subroutine Onsblock_expand
  
  subroutine get2CentInt(out,A,a1,a2,B,b1,b2,orb,ons,sp,S,U,W)
    real(dp), intent(out)          :: out(:)
    integer, intent(in)            :: A, a1, a2, B, b1, b2
    type(TOrbitals), intent(in)    :: orb    
    real(dp), intent(in)           :: ons(:,:)
    integer, intent(in)            :: sp(:)
    real(dp), intent(in)           :: S(:,:)
    real(dp), intent(in)           :: U(:,:)
    real(dp), intent(in), optional :: W(:,:,:)
    
    integer :: iSp_A, iSp_B, ishellA, ishellB
    real(dp) :: ons_A(orb%nOrbAtom(A),orb%nOrbAtom(A),2)
    real(dp) :: ons_B(orb%nOrbAtom(B),orb%nOrbAtom(B),2)
    real(dp) :: intA(2), intB(2)
    real(dp) :: Sfactor
       
    iSp_A = sp(A)
    iSp_B = sp(B)
    call getOnsME(orb,iSp_A,ons,ons_A)
    call getOnsME(orb,iSp_B,ons,ons_B)
        
    if (a1 == a2) then
      ishellA = orb%iShellOrb(a1,iSp_A)
      intA(:) = U(ishellA,iSp_A)
      if (present(W)) then
        intA(1) = intA(1) + W(ishellA,ishellA,iSp_A)
        intA(2) = intA(2) - W(ishellA,ishellA,iSp_A)
      end if
    else
      intA(:) = ons_A(a1,a2,:)
    end if
    
    if (b1 == b2) then
      ishellB = orb%iShellOrb(b1,iSp_B)
      intB(:) = U(ishellB,iSp_B)
      if (present(W)) then
        intB(1) = intB(1) + W(ishellB,ishellB,iSp_B)
        intB(2) = intB(2) - W(ishellB,ishellB,iSp_B)
      end if 
    else
      intB(:) = ons_B(b1,b2,:)
    end if
    
    Sfactor = S(b1,a1)*S(b2,a2) + S(b2,a1)*S(b1,a2)
    
    out = intA + intB
    out(:) = 0.25_dp*Sfactor*out(:)
  
  end subroutine get2CentInt

  subroutine addonsForce(deriv, DM, species, iNeighbor, nNeighbor, coords,&
      & iPair, img2CentCell, qblock, q, q0, orb, ons_en, skOverCont, skHamCont)
    real (dp), intent(inout)       :: deriv(:,:)
    real (dp), intent(in)          :: DM(:,:)
    integer, intent(in)            :: species(:)
    integer, intent(in)            :: iNeighbor(0:,:)
    integer, intent(in)            :: nNeighbor(:)
    real(dp), intent(in)           :: coords(:,:)
    integer, intent(in)            :: iPair(0:,:) 
    integer, intent(in)            :: img2CentCell(:)
    real(dp), intent(in)           :: qBlock(:,:,:,:)
    real(dp), intent(in)           :: q(:,:,:), q0(:,:,:)    
    type(TOrbitals), intent(in)    :: orb
    real(dp), intent(in)           :: ons_en(:,:)
    type(OSlakoCont), intent(in)   :: skOverCont, skHamCont

    integer  :: iOrig, iSpin, ii, nSpin, nAtom
    integer  :: iNeigh, iAtom1, iAtom2, iAtom2f, iSp1, iSp2
    integer  :: nOrb1, nOrb2, mu, nu
    real(dp) :: sqrDMTmp(orb%mOrb,orb%mOrb)
    real(dp) :: hPrimeTmp(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: ons_conv_1(orb%mOrb,orb%mOrb,2), ons_conv_2(orb%mOrb,orb%mOrb,2)
    real(dp) :: factor, fact, onethird, onefifth
    real(dp) :: Deltaq(orb%mOrb,size(q, dim=2),size(q, dim=3))
    real(dp) :: derivTmp(3)
    

    nAtom = size(orb%nOrbAtom)
    nSpin = size(qBlock, dim=4)
    Deltaq = q - q0
    
    onethird = 1.0_dp/3.0_dp
    onefifth = 1.0_dp/5.0_dp
    
    do iAtom1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAtom1)
      iSp1 = species(iAtom1)
      call getOnsME(orb,iSp1,ons_en,ons_conv_1)
      
!      print *, ">>------------"
!      print *, "atom 1", iAtom1
      
      do iNeigh = 1, nNeighbor(iAtom1)
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        if (iAtom1 == iAtom2f) then
          cycle
        end if
        
!        print *, "atom 2", iAtom2f
        
        iSp2  = species(iAtom2f)
        call getOnsME(orb,iSp2,ons_en,ons_conv_2)         
        iOrig = iPair(iNeigh,iAtom1) + 1         
        derivTmp(:) = 0.0_dp  
        call H0Sprime(hPrimeTmp, sPrimeTmp, skHamCont, skOverCont, coords, &
            &species, iAtom1, iAtom2, orb)    
           
        do iSpin = 1,nSpin
          sqrDMTmp(:,:) = 0.0_dp
          sqrDMTmp(1:nOrb2,1:nOrb1) = &
              & reshape(DM(iOrig:iOrig+nOrb1*nOrb2-1,iSpin), (/nOrb2,nOrb1/))
          if (iSpin == 1) then
            factor = 1.0_dp
          else
            factor = -1.0_dp
          end if 
          
          if (nOrb1 > 1) then
            do mu = 1, nOrb1-1
!              print *, "mu", mu
              do nu = mu+1, nOrb1
!                print *, "nu", nu
!                print *, "DeltaP(s_H,mu)", sqrDMTmp(1:nOrb2,mu)
!                print *, "DeltaP(s_H,nu)", sqrDMTmp(1:nOrb2,nu)
!                print *, "DeltaP(nu, S_H)", sqrDMTmp(nu,1:nOrb2)
!                print *, "S_prime(s_H,mu)", sPrimeTmp(1:nOrb2,mu,2)
!                print *, "S_prime(s_H,nu)", sPrimeTmp(1:nOrb2,nu,2)
                do ii = 1, 3
                  derivTmp(ii) = derivTmp(ii) + qBlock(mu,nu,iAtom1,iSpin)&
                     &*(sum(sqrDMTmp(1:nOrb2,mu)*sPrimeTmp(1:nOrb2,nu,ii))&
                     & +sum(sqrDMTmp(1:nOrb2,nu)*sPrimeTmp(1:nOrb2,mu,ii)))& 
                     &*( ons_conv_1(mu,nu,1) + factor*ons_conv_1(mu,nu,2) )
                end do                 
              end do  
            end do
          end if
          
          if (nOrb2 > 1) then
!            print *, "this should not be printed"
            do mu = 1, nOrb2-1
              do nu = mu+1, nOrb2            
                do ii = 1, 3
                  derivTmp(ii) = derivTmp(ii) + qBlock(mu,nu,iAtom2f,iSpin)&
                     &*(sum(sqrDMTmp(mu,1:nOrb1)*sPrimeTmp(nu,1:nOrb1,ii))&
                     & +sum(sqrDMTmp(nu,1:nOrb1)*sPrimeTmp(mu,1:nOrb1,ii)))& 
                     &*( ons_conv_2(mu,nu,1) + factor*ons_conv_2(mu,nu,2) )
                end do                
              end do  
            end do
          end if 
!          print *, "ons_part:", derivTmp(2)

          !RI correction (this part will be optimized)
          !p-orbitals
          if (nOrb1 > 1) then
            do mu = 2, 4
              do nu = mu,4
                fact = -1.0_dp
                if (mu == nu) fact = 1.0_dp
                do ii = 1, 3
                  derivTmp(ii) = derivTmp(ii) + fact*onethird&
                     & *(Deltaq(mu,iAtom1,iSpin) &
                     & *sum(sqrDMTmp(1:nOrb2,nu)*sPrimeTmp(1:nOrb2,nu,ii))&
                     & + Deltaq(nu,iAtom1,iSpin) &
                     & *sum(sqrDMTmp(1:nOrb2,mu)*sPrimeTmp(1:nOrb2,mu,ii)))&
                     & *( ons_en(iSp1,3) + factor*ons_en(iSp1,4))
                end do
              end do
            end do  
          end if
          if (nOrb2 > 1) then
!            print *, "this should not be printed (RI)"
            do mu = 2, 4
              do nu = mu,4
                fact = -1.0_dp
                if (mu == nu) fact = 1.0_dp
                do ii = 1, 3
                  derivTmp(ii) = derivTmp(ii) + fact*onethird&
                       & *(Deltaq(mu,iAtom2f,iSpin) &
                       & *sum(sqrDMTmp(nu,1:nOrb1)*sPrimeTmp(nu,1:nOrb1,ii))&
                       & + Deltaq(nu,iAtom2f,iSpin) &
                       & *sum(sqrDMTmp(mu,1:nOrb1)*sPrimeTmp(mu,1:nOrb1,ii)))&
                       & *( ons_en(iSp2,3) + factor*ons_en(iSp2,4))   
                end do
              end do
            end do
          end if  
!          print *, "ons+RI_part:", derivTmp(2)
          !d-orbitals
!           if (nOrb1 > 4) then
!             do mu = 5, nOrb1
!               do nu = mu, nOrb1
!                 fact = -1.0_dp
!                 if (mu == nu) fact = 2.0_dp
!                 do ii = 1, 3
!                   derivTmp(ii) = derivTmp(ii) + fact*onefifth&
!                      & *(Deltaq(mu,iAtom1,iSpin) &
!                      & *sum(sqrDMTmp(1:nOrb2,nu)*sPrimeTmp(1:nOrb2,nu,ii))&
!                      & + Deltaq(nu,iAtom1,iSpin) &
!                      & *sum(sqrDMTmp(1:nOrb2,mu)*sPrimeTmp(1:nOrb2,mu,ii)))&
!                      & *( ons_en(iSp1,9) + factor*ons_en(iSp1,10))
!                 end do
!               end do
!             end do
!           end if
!           if (nOrb2 > 4) then        
!             do mu = 5, nOrb2
!               do nu = mu, nOrb2
!                 fact = -1.0_dp
!                 if (mu == nu) fact = 2.0_dp
!                 do ii = 1, 3
!                   derivTmp(ii) = derivTmp(ii) + fact*onefifth&
!                      & *(Deltaq(mu,iAtom2f,iSpin) &
!                      & *sum(sqrDMTmp(nu,1:nOrb1)*sPrimeTmp(nu,1:nOrb1,ii))&
!                      & + Deltaq(nu,iAtom2f,iSpin) &
!                      & *sum(sqrDMTmp(mu,1:nOrb1)*sPrimeTmp(mu,1:nOrb1,ii)))&
!                      & *( ons_en(iSp2,9) + factor*ons_en(iSp2,10))
!                 end do
!               end do
!             end do
!           end if
                
        end do  

        deriv(:,iAtom1) = deriv(:,iAtom1) + derivTmp(:)
        deriv(:,iAtom2f) = deriv(:,iAtom2f) - derivTmp(:)
             
      end do     
    end do

  end subroutine addonsForce  
   
  
end module onscorrection