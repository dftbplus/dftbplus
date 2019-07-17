!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains functions to calculate the Mataga-Nishimoto-Ohno-Klopman damped interactions
module dftbp_mnok
  use dftbp_accuracy
  use dftbp_message
  implicit none

  private

contains

        !> Mataga-Nishimoto-Ohno-Klopman real space case
        pure function gamma_cluster_l(etaA, kappaAl, etaB, kappaBl, Rab) result(gamma)
                
                !> Fitted hardness for atom A
                real(dp), intent(in) :: etaA

                !> Scaling factor for shell l on atom A
                real(dp), intent(in) :: kappaAl

                !> Fitted hardness for atom B
                real(dp), intent(in) :: etaB
                
                !> Scaling factor for shell l' on atom B
                real(dp), intent(in) :: kappaBl

                !> Distance between atoms
                real(dp), intent(in) :: Rab

                !> Resulting interaction between atoms
                real(dp) :: gamma

                gamma = sqrt(1.0_dp / (Rab**2 + 1.0_dp / eta(etaA, kappaAl, etaB, kappaBl)**2 ))

        end function gamma_cluster_l

        
        !> effective hardness parameter for shells of atoms
        pure function eta(etaA, kappaAl, etaB, kappaBl)
                
                !> Fitted hardness for atom A
                real(dp), intent(in) :: etaA

                !> Scaling factor for shell l on atom A
                real(dp), intent(in) :: kappaAl

                !> Fitted hardness for atom B
                real(dp), intent(in) :: etaB
                
                !> Scaling factor for shell l' on atom B
                real(dp), intent(in) :: kappaBl

                eta = 0.5_dp * ( (1.0_dp + kappaAl) * etaB + (1.0_dp + kappaBl) )

        end function eta

end module

