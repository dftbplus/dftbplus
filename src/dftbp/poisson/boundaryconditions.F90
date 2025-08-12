!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

module dftbp_poisson_boundaryconditions

  !> Enumerator containing possible Poisson boundary conditions
  type, private :: TBCPoissonEnum_

    !> Unset
    integer :: unset = -1

    !> Periodic on face
    integer :: periodic = 0

    !> Potential specified at the cell edge
    integer :: dirichlet = 1

    !> Derivative specified at the cell edge
    integer :: neumann = 2

  end type TBCPoissonEnum_

  !> Actual instance of the boundary condition enumerator
  type(TBCPoissonEnum_), parameter, public :: poissonBCsEnum = TBCPoissonEnum_()

  ! Note, order corresponds to TBCPoissonEnum_
  character(10), parameter, public :: bcPoissonNames(-1:2) =&
      & [ character(10) :: "Unset", "Periodic", "Dirichlet", "Neumann" ]

end module dftbp_poisson_boundaryconditions
