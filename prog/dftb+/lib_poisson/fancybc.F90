!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *  
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program * 
!  under the LGPL licence.                                                *
!**************************************************************************
module fancybc

  use dftbp_accuracy, only : dp
  use parameters
  use bulkpot , only : super_array
  use mpi_poisson
  
  implicit none
  private

  public :: bndyc, coef, coef_period, coef_gate, gate_bound
  public :: coef_tip ,tip_bound, coef_cilgate, cilgate_bound
  public :: coef_local, local_bound
  public :: cofx, cofy, cofz
  public :: mix_bndyc

contains
  
Subroutine bndyc(kbdy,xory,yorz,alfa,gbdy)

 integer :: kbdy
 real(kind=dp) :: xory,yorz,alfa,gbdy

 alfa = 0.d0
 gbdy = 0.d0

end subroutine bndyc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine cofx(x,cxx,cx,cex)
  real(dp) :: x, cxx, cx, cex

  cxx = 1.d0
  cx = 0.d0
  cex = 0.d0

end subroutine cofx
! ----------------------------------------------------------------------
Subroutine cofy(y,cyy,cy,cey)
  real(dp) :: y, cyy, cy, cey

  cyy = 1.d0
  cy = 0.d0
  cey = 0.d0

end subroutine cofy
! ----------------------------------------------------------------------
Subroutine cofz(z,czz,cz,cez)
  real(dp) :: z, czz, cz, cez

  czz = 1.d0
  cz = 0.d0
  cez = 0.d0

end subroutine cofz
! ----------------------------------------------------------------------

Subroutine coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)

  real(kind=dp) :: x,y,z,cxx,cyy,czz,cx,cy,cz,ce
  
   cxx = 1.d0
   cyy = 1.d0
   czz = 1.d0
   cx = 0.d0          
   cy = 0.d0  
   cz = 0.d0
   ce = 0.d0    

   if(any(localBC.gt.0)) call coef_local(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)

end subroutine coef
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine coef_period(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)

  real(kind=dp) :: x,y,z,cxx,cyy,czz,cx,cy,cz,ce
  
   cxx = 1.d0
   cyy = 1.d0
   czz = 1.d0
   cx = 0.d0          
   cy = 0.d0  
   cz = 0.d0
   ce = 0.d0    


end subroutine coef_period
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine coef_gate(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)

real(kind=dp) :: cxx,cyy,czz,cx,cy,cz,ce
real(kind=dp) :: x,y,z,x_x,y_y,z_z
real(kind=dp) :: x_min_gate,x_max_gate,y_min_gate,y_max_gate
real(kind=dp) :: z_min_gate,z_max_gate
integer :: i_x,i_y,i_z


   cxx = 1.d0
   cyy = 1.d0
   czz = 1.d0
   cx = 0.d0          
   cy = 0.d0  
   cz = 0.d0
   ce = 0.d0    

   if(any(localBC.gt.0)) call coef_local(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)

   select case (gatedir)
   case(1)
       x_x = y
       y_y = z
       z_z = x
       i_x = 2
       i_y = 3
       i_z = 1
   case(2)  
       x_x = x
       y_y = z
       z_z = y
       i_x = 1
       i_y = 3
       i_z = 2
   case(3)
       x_x = x
       y_y = y    
       z_z = z
       i_x = 1 
       i_y = 2 
       i_z = 3
   end select

     ! The gate extends on the centre of Poisson Box
     ! If Rmin_Gate is positive the gate is on +gatedir
     ! If Rmin_Gate is negative the gate is on -gatedir

        x_min_gate=cntr_gate(i_x)-GateLength_t/2.d0
        x_max_gate=cntr_gate(i_x)+GateLength_t/2.d0
        y_min_gate=cntr_gate(i_y)-GateLength_l/2.d0
        y_max_gate=cntr_gate(i_y)+GateLength_l/2.d0

        if (Rmin_Gate.gt.0.d0) then
           z_min_gate=cntr_gate(i_z)+Rmin_Gate
           z_max_gate=cntr_gate(i_z)+PoissBox(i_z,i_z)
        else
           z_min_gate=cntr_gate(i_z)-PoissBox(i_z,i_z)
           z_max_gate=cntr_gate(i_z)-Rmin_Gate           
        end if   

        if(x_x.ge.x_min_gate.and.x_x.le.x_max_gate) then
           if(y_y.ge.y_min_gate.and.y_y.le.y_max_gate) then
              if(z_z.ge.z_min_gate.and.z_z.le.z_max_gate) then

                   cxx = 0.d0
                   cyy = 0.d0
                   czz = 0.d0
                   cx = 0.d0          
                   cy = 0.d0  
                   cz = 0.d0
                   ce = 1.d0   

               end if
            end if
         end if

  return

end subroutine coef_gate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine gate_bound(iparm,fparm,dlx,dly,dlz,rhs)

 integer :: iparm(23)
 real(kind=dp) :: fparm(8)
 real(kind=dp) :: dlx,dly,dlz
 real(kind=dp) :: rhs(iparm(14),iparm(15),iparm(16))
 logical :: period 

 real(kind=dp) :: x_x,y_y,z_z,xi,yj,zk
 real(kind=dp) :: x_min_gate,x_max_gate,y_min_gate,y_max_gate
 real(kind=dp) :: z_min_gate,z_max_gate
 integer :: i,j,k,i_x,i_y,i_z
 
   select case (gatedir)
   case(1)
       i_x = 2
       i_y = 3
       i_z = 1
   case(2)  
       i_x = 1
       i_y = 3
       i_z = 2
   case(3)
       i_x = 1 
       i_y = 2 
       i_z = 3
   end select

   x_min_gate=cntr_gate(i_x)-GateLength_t/2.d0
   x_max_gate=cntr_gate(i_x)+GateLength_t/2.d0
   y_min_gate=cntr_gate(i_y)-GateLength_l/2.d0
   y_max_gate=cntr_gate(i_y)+GateLength_l/2.d0
          
   if (Rmin_Gate.gt.0.d0) then
      z_min_gate=cntr_gate(i_z)+Rmin_Gate
      z_max_gate=cntr_gate(i_z)+PoissBox(i_z,i_z)
   else
      z_min_gate=cntr_gate(i_z)-PoissBox(i_z,i_z)
      z_max_gate=cntr_gate(i_z)-Rmin_Gate           
   end if

   do i = 1,iparm(14)  
     do j = 1,iparm(15)
       do k = 1,iparm(16)

          xi = fparm(1) + (i - 1)*dlx
          yj = fparm(3) + (j - 1)*dly    
          zk = fparm(5) + (k - 1)*dlz

          select case (gatedir)
          case(1)
             x_x = yj
             y_y = zk
             z_z = xi
          case(2)  
             x_x = xi
             y_y = zk
             z_z = yj
          case(3)
             x_x = xi
             y_y = yj   
             z_z = zk
          end select

          ! The gate extends on the centre of Poisson Box
          ! If Rmin_Gate is positive the gate is on +gatedir
          ! If Rmin_Gate is negative the gate is on -gatedir
 
          if(z_z.ge.z_min_gate.and.z_z.le.z_max_gate) then
             if(x_x.ge.x_min_gate.and.x_x.le.x_max_gate) then
                if(y_y.ge.y_min_gate.and.y_y.le.y_max_gate) then

                   rhs(i,j,k)=gate
                   
                endif
             endif
          endif


       enddo
    enddo
 enddo

 end subroutine gate_bound
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine coef_cilgate(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
 
 real(kind=dp) :: cxx,cyy,czz,cx,cy,cz,ce
 real(kind=dp) :: x,y,z,a,b,a1,a2,b1,b2
 real(kind=dp) :: x_min_gate,x_max_gate,y_min_gate,y_max_gate
 real(kind=dp) :: z_min_gate,z_max_gate,d_cntr
 real(kind=dp) :: x_min_ox,x_max_ox,y_min_ox,y_max_ox
 real(kind=dp) :: z_min_ox,z_max_ox
 real(kind=dp) :: x_x,y_y,z_z,zero(3)
 integer :: i_x, i_y, i_z


  select case(biasdir)

     case(1)

        x_x = x
        i_x= 1
        y_y = y
        z_z = z
        i_y= 2
        i_z= 3
        zero=(/1.d0,0.d0,0.d0/)

     case(2)

        x_x = y
        i_x = 2
        y_y = x
        z_z = z
        i_y= 1
        i_z= 3
        zero=(/0.d0,1.d0,0.d0/)

     case(3)

        x_x = z
        i_x= 3
        y_y = x
        z_z = y
        i_y= 1
        i_z= 2
        zero=(/0.d0,0.d0,1.d0/)

     end select

     x_min_gate = cntr_gate(i_x) - GateLength_l/2.d0  
     x_max_gate = cntr_gate(i_x) + GateLength_l/2.d0
     x_min_ox = cntr_gate(i_x) - OxLength/2.d0  
     x_max_ox = cntr_gate(i_x) + OxLength/2.d0


     a = (eps_r - 1.d0)/dr_eps
     b = (eps_r + 1.d0)/2.d0 - a*Rmin_Ins
     
     a1 = -(eps_r - 1.d0)/dr_eps
     b1 = (eps_r + 1.d0)/2.d0 - a1*x_max_ox
     
     a2 = (eps_r - 1.d0)/dr_eps
     b2 = (eps_r + 1.d0)/2.d0 - a2*x_min_ox


     d_cntr = sqrt((y_y-cntr_gate(i_y))**2+(z_z-cntr_gate(i_z))**2)

     ! First checks whether (xyz) lies outside gate regions
     if (((x_x.lt.(x_min_ox-dr_eps/2.d0)).or.(x_x.gt.(x_max_ox+dr_eps/2.d0))).or. &
          (((x_x.ge.(x_min_ox-dr_eps/2.d0)).and.(x_x.le.(x_max_ox+dr_eps/2.d0))).and. &
          (d_cntr.lt.(Rmin_Ins-dr_eps/2.d0)))) then
     
        cxx = 1.d0
        cyy = 1.d0
        czz = 1.d0
        
        cx = 0.d0          
        cy = 0.d0
        cz = 0.d0
        
        ce = 0.d0

        if(any(localBC.gt.0)) call coef_local(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)

     else

        if ((d_cntr.lt.(Rmin_Ins+dr_eps/2.d0)).and. &
             ((x_x.le.(-d_cntr+x_max_ox+Rmin_Ins)).and.(x_x.ge.(d_cntr+x_min_ox-Rmin_Ins)))) then  

           cxx = a*d_cntr + b
           cyy = a*d_cntr + b
           czz = a*d_cntr + b
           
           cx = (1.d0-zero(1)) * a * ( x - cntr_gate(1) )/d_cntr           
           cy = (1.d0-zero(2)) * a * ( y - cntr_gate(2) )/d_cntr  
           cz = (1.d0-zero(3)) * a * ( z - cntr_gate(3) )/d_cntr  
           
           ce = 0.d0             
  
        else

           if (x_x.gt.(x_max_ox-dr_eps/2.d0)) then
              
              cxx = a1*x_x + b1 
              cyy = a1*x_x + b1
              czz = a1*x_x + b1
              
              cx = zero(1)*a1
              cy = zero(2)*a1
              cz = zero(3)*a1
              
              ce = 0.d0
              
           else
              
              if (x_x.lt.(x_min_ox+dr_eps/2.d0)) then
                 
                 cxx = a2*x_x + b2 
                 cyy = a2*x_x + b2
                 czz = a2*x_x + b2
                 
                 cx = zero(1)*a2
                 cy = zero(2)*a2
                 cz = zero(3)*a2
                 
                 ce = 0.d0 
                 
              else   
                 
                 if (((d_cntr.lt.Rmin_Gate)).or. &
                      ((d_cntr.ge.Rmin_Gate).and.((x_x.lt.x_min_gate).or.(x_x.gt.x_max_gate)))) then 
                    
                    cxx = eps_r           
                    cyy = eps_r
                    czz = eps_r
                    
                    cx = 0.d0
                    cy = 0.d0
                    cz = 0.d0
                    
                    ce = 0.d0
                    
                 else
                    
                    cxx = 0.d0
                    cyy = 0.d0
                    czz = 0.d0
                    
                    cx = 0.d0
                    cy = 0.d0
                    cz = 0.d0
                    
                    ce = 1.d0
                    
                 end if
                 
              end if
              
           end if
           
        end if
        
     end if


end subroutine coef_cilgate
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine cilgate_bound(iparm,fparm,dlx,dly,dlz,rhs)

 integer :: iparm(23)
 real(kind=dp) :: fparm(8)
 real(kind=dp) :: dlx,dly,dlz
 real(kind=dp) :: rhs(iparm(14),iparm(15),iparm(16))
 logical :: period 

 integer :: i,j,k
 real(kind=dp) :: xi,yj,zk
 real(kind=dp) :: x_min_gate,x_max_gate,y_min_gate,y_max_gate
 real(kind=dp) :: z_min_gate,z_max_gate,d_cntr
 real(kind=dp) :: x_min_ox,x_max_ox,y_min_ox,y_max_ox
 real(kind=dp) :: z_min_ox,z_max_ox


 do i = 1,iparm(14)  
    do j = 1,iparm(15)
       do k = 1,iparm(16)

         xi = fparm(1) + (i - 1)*dlx
         yj = fparm(3) + (j - 1)*dly    
         zk = fparm(5) + (k - 1)*dlz

         
         selectcase(biasdir)

           case(1)

             x_min_gate = cntr_gate(1) - GateLength_l/2.d0  
             x_max_gate = cntr_gate(1) + GateLength_l/2.d0
             x_min_ox = cntr_gate(1) - OxLength/2.d0  
             x_max_ox = cntr_gate(1) + OxLength/2.d0

             d_cntr = sqrt((yj-cntr_gate(2))**2+(zk-cntr_gate(3))**2)

             !if ( (d_cntr.gt.(Rmin_Ins-dr_eps/2.d0)).and.(xi.ge.x_min_ox).and.(xi.le.x_max_ox) ) then
             !   rhs(i,j,k) = 0.d0
             !endif

             if ( (d_cntr.ge.Rmin_Gate).and.(xi.ge.x_min_gate).and.(xi.le.x_max_gate) ) then                  
                rhs(i,j,k) = gate   
             endif

           case(2)

             y_min_gate = cntr_gate(2) - GateLength_l/2.d0  
             y_max_gate = cntr_gate(2) + GateLength_l/2.d0
             y_min_ox = cntr_gate(2) - OxLength/2.d0  
             y_max_ox = cntr_gate(2) + OxLength/2.d0

             d_cntr = sqrt((xi-cntr_gate(1))**2+(zk-cntr_gate(3))**2)

             !if ( (d_cntr.gt.(Rmin_Ins-dr_eps/2.d0)).and.(yj.ge.y_min_ox).and.(yj.le.y_max_ox) ) then
             !   rhs(i,j,k) = 0.d0
             !endif

             if ( (d_cntr.ge.Rmin_Gate).and.(yj.ge.y_min_gate).and.(yj.le.y_max_gate) ) then
                 rhs(i,j,k) = gate   
             endif

           case(3)

             z_min_gate = cntr_gate(3) - GateLength_l/2.d0  
             z_max_gate = cntr_gate(3) + GateLength_l/2.d0
             z_min_ox = cntr_gate(3) - OxLength/2.d0  
             z_max_ox = cntr_gate(3) + OxLength/2.d0

             d_cntr = sqrt((xi-cntr_gate(1))**2+(yj-cntr_gate(2))**2)

             !if ( (d_cntr.gt.(Rmin_Ins-dr_eps/2.d0)).and.(zk.ge.z_min_ox).and.(zk.le.z_max_ox) ) then
             !   rhs(i,j,k) = 0.d0
             !endif

             if ( (d_cntr.ge.Rmin_Gate).and.(zk.ge.z_min_gate).and.(zk.le.z_max_gate) ) then
                rhs(i,j,k) = gate   
             endif

          end select

       end do
    end do
 end do

end subroutine cilgate_bound
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine coef_tip(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)

 real(kind=dp) :: cxx,cyy,czz,cx,cy,cz,ce
 real(kind=dp) :: x_x,y_y,z_z
 real(kind=dp) :: x,y,z
 real(kind=dp) :: x_min_tip,x_max_tip,y_min_tip,y_max_tip
 real(kind=dp) :: z_min_tip,z_max_tip,tip_high
 real(kind=dp) :: y_min_plane,y_max_plane,z_min_plane,z_max_plane
 integer :: i_x,i_y,i_z

   cxx = 1.d0
   cyy = 1.d0
   czz = 1.d0
   cx = 0.d0          
   cy = 0.d0  
   cz = 0.d0
   ce = 0.d0    

  selectcase(gatedir)

     case(1)

        x_x = x
        y_y = y
        z_z = z
        i_x= 1
        i_y= 2
        i_z= 3

     case(2)

        x_x = y
        y_y = x
        z_z = z
        i_x= 2
        i_y= 1
        i_z= 3

     case(3)

        x_x = z
        y_y = y
        z_z = x
        i_x= 3
        i_y= 2
        i_z= 1

     
     end select

        
        x_min_tip=min( tip_atom(i_x),base_atom1(i_x) )
        x_max_tip=max( tip_atom(i_x),base_atom1(i_x) )
        tip_high=(x_max_tip-x_min_tip)      

        if (x_min_tip.eq.base_atom1(i_x)) then
          if (x_x.lt.x_min_tip) then
            cxx = 0.d0
            cyy = 0.d0
            czz = 0.d0
            cx = 0.d0
            cy = 0.d0
            cz = 0.d0
            ce = 1.d0
          endif
        end if
        if (x_min_tip.eq.tip_atom(i_x)) then
          if(x_x.gt.x_max_tip) then
            cxx = 0.d0
            cyy = 0.d0
            czz = 0.d0
            cx = 0.d0
            cy = 0.d0
            cz = 0.d0
            ce = 1.d0
          endif
        end if

        if(x_x.gt.x_min_tip.and.x_x.lt.x_max_tip) then

             y_min_tip=min(base_atom1(i_y),base_atom2(i_y))
             y_max_tip=max(base_atom1(i_y),base_atom2(i_y))
             z_min_tip=min(base_atom1(i_z),base_atom2(i_z))
             z_max_tip=max(base_atom1(i_z),base_atom2(i_z))

             y_max_plane=tip_atom(i_y)+(y_max_tip-tip_atom(i_y))*abs(x_x-tip_atom(i_x))/tip_high
             y_min_plane=tip_atom(i_y)+(y_min_tip-tip_atom(i_y))*abs(x_x-tip_atom(i_x))/tip_high
             z_max_plane=tip_atom(i_z)+(z_max_tip-tip_atom(i_z))*abs(x_x-tip_atom(i_x))/tip_high
             z_min_plane=tip_atom(i_z)+(z_min_tip-tip_atom(i_z))*abs(x_x-tip_atom(i_x))/tip_high

             if(y_y.gt.y_min_plane.and.y_y.lt.y_max_plane) then
                if(z_z.gt.z_min_plane.and.z_z.lt.z_max_plane) then

                   cxx = 0.d0
                   cyy = 0.d0
                   czz = 0.d0
                   cx = 0.d0          
                   cy = 0.d0  
                   cz = 0.d0
                   ce = 1.d0   

                endif
              endif
             
          endif

  return

end subroutine coef_tip
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine tip_bound(iparm,fparm,dlx,dly,dlz,rhs)

 integer :: iparm(23)
 real(kind=dp) :: fparm(8)
 real(kind=dp) :: dlx,dly,dlz
 real(kind=dp) :: rhs(iparm(14),iparm(15),iparm(16))

 real(kind=dp) :: x_x,y_y,z_z
 real(kind=dp) :: xi,yj,zk
 real(kind=dp) :: x_min_tip,x_max_tip,y_min_tip,y_max_tip
 real(kind=dp) :: z_min_tip,z_max_tip,tip_high
 real(kind=dp) :: y_min_plane,y_max_plane,z_min_plane,z_max_plane
 integer :: i,j,k,i_x,i_y,i_z

 do i = 1,iparm(14)  
    do j = 1,iparm(15)
       do k = 1,iparm(16)
          xi = fparm(1) + (i - 1)*dlx
          yj = fparm(3) + (j - 1)*dly    
          zk = fparm(5) + (k - 1)*dlz


          selectcase(gatedir)
             
          case(1)
             
             x_x = xi
             y_y = yj
             z_z = zk
             i_x= 1
             i_y= 2
             i_z= 3
             
          case(2)
             
             x_x = yj
             y_y = xi
             z_z = zk
             i_x= 2
             i_y= 1
             i_z= 3
             
          case(3)
             
             x_x = zk
             y_y = yj
             z_z = xi
             i_x= 3
             i_y= 2
             i_z= 1
             
             
          end select


          x_min_tip=min( tip_atom(i_x),base_atom1(i_x) )
          x_max_tip=max( tip_atom(i_x),base_atom1(i_x) )      
          tip_high=x_max_tip-x_min_tip
          

          if (x_min_tip.eq.base_atom1(i_x)) then
            if (x_x.lt.x_min_tip) then 
                rhs(i,j,k) = tipbias
            endif 
          endif
          if (x_min_tip.eq.tip_atom(i_x)) then
            if(x_x.gt.x_max_tip) then 
                 rhs(i,j,k) = tipbias
            endif
          endif

          if (x_x.gt.x_min_tip.and.x_x.lt.x_max_tip) then

             y_min_tip=min(base_atom1(i_y),base_atom2(i_y))
             y_max_tip=max(base_atom1(i_y),base_atom2(i_y))
             z_min_tip=min(base_atom1(i_z),base_atom2(i_z))
             z_max_tip=max(base_atom1(i_z),base_atom2(i_z))

             y_max_plane=tip_atom(i_y)+(y_max_tip-tip_atom(i_y))*abs(x_x-tip_atom(i_x))/tip_high
             y_min_plane=tip_atom(i_y)+(y_min_tip-tip_atom(i_y))*abs(x_x-tip_atom(i_x))/tip_high
             z_max_plane=tip_atom(i_z)+(z_max_tip-tip_atom(i_z))*abs(x_x-tip_atom(i_x))/tip_high
             z_min_plane=tip_atom(i_z)+(z_min_tip-tip_atom(i_z))*abs(x_x-tip_atom(i_x))/tip_high

             if(y_y.gt.y_min_plane.and.y_y.lt.y_max_plane) then
                if(z_z.gt.z_min_plane.and.z_z.lt.z_max_plane) then
                   
                  rhs(i,j,k) = tipbias

                endif
            endif
         endif
     
     enddo
   enddo
enddo

if(id0) close(90)
 
end subroutine tip_bound
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

Subroutine coef_local(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)

 real(kind=dp) :: x,y,z, cxx,cyy,czz,cx,cy,cz,ce
 
 
 integer :: m,i_x,i_y,i_z      ! m = number of contact
 real(kind=dp) :: xx,yy,zz,R
 

 cxx = 1.d0
 cyy = 1.d0
 czz = 1.d0
 cx = 0.d0          
 cy = 0.d0  
 cz = 0.d0
 ce = 0.d0    
 
 do m = 1, ncont
    
    selectcase( abs(contdir(m)) )
    case(1)           
       xx=x; yy=y; zz=z;
       i_x=1; i_y=2; i_z=3;
    case(2)
       xx=y; yy=z; zz=x;
       i_x=2; i_y=3; i_z=1;
    case(3)
       xx=z; yy=x; zz=y;
       i_x=3; i_y=1; i_z=2;
    end select
    
    ! The contacts extends on the face of Poisson Box
    if(abs(xx-cntr_cont(i_x,m)).le.1d-10) then  
       select case(localBC(m))
       case(1)
          R=sqrt((zz-cntr_cont(i_z,m))**2+(yy-cntr_cont(i_y,m))**2)
          if(R.le.R_cont(m)) then
             cxx = 0.d0
             cyy = 0.d0
             czz = 0.d0
             cx = 0.d0          
             cy = 0.d0  
             cz = 0.d0
             ce = 1.d0 
          endif
          
       case(2) 

          if(zz.ge.cntr_cont(i_z,m)-R_cont(m).and.zz.le.cntr_cont(i_z,m)+R_cont(m).and. &
             yy.ge.cntr_cont(i_y,m)-R_cont(m).and.yy.le.cntr_cont(i_y,m)+R_cont(m)) then
             
             cxx = 0.d0
             cyy = 0.d0
             czz = 0.d0
             cx = 0.d0          
             cy = 0.d0  
             cz = 0.d0
             ce = 1.d0   
          endif
          
       end select 

    endif

 enddo !m       

end subroutine coef_local
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine local_bound(iparm,fparm,x,rhs,phi_bulk)

 integer :: iparm(23)
 real(kind=dp) :: fparm(8)
 real(dp), dimension(:,:) :: x
 real(kind=dp) :: rhs(iparm(14),iparm(15),iparm(16))
 Type(super_array) :: phi_bulk(ncont)

 !Common variables         
 integer :: m,jmax,kmax,i_y,i_z,i,j,k                ! m = number of contact

 real(kind=dp) :: Rmin_x,Rmin_y,ymin,zmin,dy,dz,R,yj,zk
 real(kind=dp) :: x_min,y_min,x_max,y_max,p
 
    do m = 1, ncont
       
       x_min=minval(x(phi_bulk(m)%a,iatc(3,m):iatc(2,m))) !along a
       x_max=maxval(x(phi_bulk(m)%a,iatc(3,m):iatc(2,m)))
       y_min=minval(x(phi_bulk(m)%b,iatc(3,m):iatc(2,m))) !along b
       y_max=maxval(x(phi_bulk(m)%b,iatc(3,m):iatc(2,m)))

       if(contdir(m).gt.0) then 
          p=fparm(2*phi_bulk(m)%c)
          k=iparm(contdir(m)+13)
       else ! p is the position of the contact on the face of poisson box
          p=fparm(2*phi_bulk(m)%c - 1)
          k=1
       endif
       jmax=phi_bulk(m)%iparm(14)  ! along a
       kmax=phi_bulk(m)%iparm(15)  ! along b
       ymin=fparm(2*phi_bulk(m)%a - 1)
       zmin=fparm(2*phi_bulk(m)%b - 1)
       dy=phi_bulk(m)%dla
       dz=phi_bulk(m)%dlb
       i_y=phi_bulk(m)%a
       i_z=phi_bulk(m)%b

       !	R of contact  
       Rmin_x=(x_max-x_min)/2.d0
       Rmin_y=(y_max-y_min)/2.d0
       R_cont(m)=max(Rmin_x,Rmin_y)+dR_cont(m)

       select case ( abs(contdir(m)) )
       case(1)
          cntr_cont(:,m)=(/ p,(x_max-Rmin_x),(y_max-Rmin_y) /)  
       case(2)
          cntr_cont(:,m)=(/ (y_max-Rmin_y),p,(x_max-Rmin_x) /)
       case(3)
          cntr_cont(:,m)=(/ (x_max-Rmin_x),(y_max-Rmin_y),p /)
       end select

       do i = 1,jmax
          do j = 1,kmax
             
             yj = ymin + (i - 1)*dy
             zk = zmin + (j - 1)*dz
             
             select case(localBC(m))
             case(1)
                ! The contacts extends on the faces of Poisson Box
                R=sqrt((zk-cntr_cont(i_z,m))**2+(yj-cntr_cont(i_y,m))**2)
                if(R.le.R_cont(m)) then
                   select case ( abs(contdir(m)) )
                   case(1)
                      rhs(k,i,j) = phi_bulk(m)%val(i,j,1)
                   case(2)
                      rhs(j,k,i) = phi_bulk(m)%val(i,j,1)
                   case(3)
                      rhs(i,j,k) = phi_bulk(m)%val(i,j,1)
                   end select
                endif
             case(2)
                if(yj.ge.cntr_cont(i_y,m)-R_cont(m).and.yj.le.cntr_cont(i_y,m)+R_cont(m).and. &
                   zk.ge.cntr_cont(i_z,m)-R_cont(m).and.zk.le.cntr_cont(i_z,m)+R_cont(m)) then
                   select case ( abs(contdir(m)) )
                   case(1)
                      rhs(k,i,j) = phi_bulk(m)%val(i,j,1)
                   case(2)
                      rhs(j,k,i) = phi_bulk(m)%val(i,j,1)
                   case(3)
                      rhs(i,j,k) = phi_bulk(m)%val(i,j,1)
                   end select
                endif
             end select

          enddo
       enddo
 
    enddo !m


end subroutine local_bound


Subroutine mix_bndyc(kbdy,xory,yorz,alfa,gbdy)

 integer :: kbdy
 real(kind=dp) :: xory,yorz,alfa,gbdy

 integer :: m, ii, jj, i_z, i_y
 real(dp) :: R, yy, zz
 ! kbdy=1: xa; xory = y
 ! kbdy=2: xb; yorz = z
 ! kbdy=3: ya; xory = x
 ! kbdy=4: yb; yorz = z
 ! kbdy=5: za; xory = x
 ! kbdy=6: zb; yorz = y

 ! Standard Neumann BC 
 
 alfa = 0.d0
 gbdy = 0.d0
 
 if (mixed(kbdy)) then 
    
    yy=xory; zz=yorz;

    do m = 1, ncont
       !we should map yy, zz onto the grid points ii, jj 
       selectcase( abs(contdir(m)) )
       case(1) ! yy=y zz=z           
          i_y=2; i_z=3;
          !ii=nint((yy-PoissBounds(2,1))/phi_bulk(m)%dla)+1
          !jj=nint((zz-PoissBounds(3,1))/phi_bulk(m)%dlb)+1
       case(2) ! yy=x zz=z
          i_y=1; i_z=3;
          !ii=nint((yy-PoissBounds(1,1))/phi_bulk(m)%dla)+1
          !jj=nint((zz-PoissBounds(3,1))/phi_bulk(m)%dlb)+1
       case(3) ! yy=x zz=y
          i_y=1; i_z=2;
          !ii=nint((yy-PoissBounds(1,1))/phi_bulk(m)%dla)+1
          !jj=nint((zz-PoissBounds(2,1))/phi_bulk(m)%dlb)+1
       end select
       
       select case(localBC(m))
       case(1)
          
          R=sqrt((zz-cntr_cont(i_z,m))**2+(yy-cntr_cont(i_y,m))**2)
          if(R.le.R_cont(m)) then
             alfa = 1.d6
             gbdy = 1.d6 ! bulk_pot(m)%val(ii,jj,1)
          endif
          
       case(2) 
          
          if(zz.ge.cntr_cont(i_z,m)-R_cont(m).and.zz.le.cntr_cont(i_z,m)+R_cont(m).and. &
             yy.ge.cntr_cont(i_y,m)-R_cont(m).and.yy.le.cntr_cont(i_y,m)+R_cont(m)) then

             alfa = 1.d6
             gbdy = 1.d6 ! bulk_pot(m)%val(ii,jj,1)
             
          endif
          
       end select
       
    enddo !m
    
 end if 

end subroutine mix_bndyc

end module fancybc
