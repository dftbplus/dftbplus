!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *  
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program * 
!  under the LGPL licence.                                                *
!**************************************************************************

#:include 'common.fypp'

module mpi_poisson

#:if WITH_MPI
  use libmpifx_module
  implicit none
  private  

  public :: mpifx_gatherv
  public :: poiss_mpi_init, poiss_mpi_split

  !! mpifx communicator as received by calling program
  type(mpifx_comm), public :: global_comm
  !! Local communicator (calculations of rhs and shifts)
  type(mpifx_comm), public :: poiss_comm
  integer, public ::  numprocs
  integer, public ::  id
  logical, public ::  id0
  logical, public :: active_id

  contains

  subroutine  poiss_mpi_init(comm)

      type(mpifx_comm), intent(in) :: comm

      global_comm = comm 

      !! If no split occurs, then the local and the global communicator 
      !! must be the same
      poiss_comm = comm
      id = comm%rank
      numprocs = comm%size
      
      ! i/o node that also actually solves the multigrid
      id0 = (id == 0)
      
      ! nodes parallelizing r.h.s. assembly
      active_id = (id == 0) 

  end subroutine poiss_mpi_init

  !-----------------------------------------------------------------------
  !   MPI PARALLELIZATION STRATEGY
  !----------------------------------------------------------------------- 
  !   case even number of processors:
  !   suppose maxnumproc = 3 => split into 3 groups
  !
  !   global_comm:         0  1  2  3  4  5  6  7 
  !   color      :         0  0  0  1  1  1  2  2 
  !   key        :         0  1  2  0  1  2  0  1
  !
  !   temp_comm  :         0  1  2  0  1  2  0  1 
  !   active:              Y  N  N  Y  N  N  Y  N
  !   id0:                 Y  N  N  N  N  N  N  N
  !
  !   poiss_comm :         0        1        2  
  !
  !   HENCE: it is possible to fit ngroups <= numprocs/2
  !       => numprocs/ngroups >= 2
  !
  !------------------------------------------------------------------------
  !   case odd number of processors:
  !   suppose numprocs = 5; maxnumprocs = 3
  !   global_comm:         0  1  2  3  4
  !   color      :         0  0  1  1  2
  !   temp_comm  :         0  1  0  1  0 
  !   
  !   HENCE: it is possible to fit ngroups <= (numprocs+1)/2
  !   HENCE: (numprocs+1)/ngroups >= 2
  !    
  !   EXCEPTION: if maxnumproc = numprocs => then poiss_comm = global_comm
  !  
  !   only active nodes are parallelizing the rhs assembly
  !   only id0 is actually solving multigrid  
  !------------------------------------------------------------------------
  subroutine poiss_mpi_split(maxnumprocs)
      integer, intent(in) :: maxnumprocs

      integer :: ngroups, npg, cf, color, key, ierr
      type(mpifx_comm) :: temp_comm
         
      cf=0
      if (mod(numprocs,2)>0) cf = 1
      ngroups = maxnumprocs
      ! number of processors per group
      npg = (numprocs+cf)/ngroups

      if (npg >= 2) then
        ! this correction is needed when division is not perfect
        ! but we know at least one more id=0 will fit 
        if (mod(numprocs+cf, ngroups) > 0) npg = npg + 1
      else !(npg<2)
        ! if npg < 2 then ngroups is limited to npg=2     
        ngroups = (numprocs+cf)/2 
        npg = 2 
      end if

      ! ovverride the exception case
      if (numprocs == maxnumprocs) npg = 1      
      
      color = global_comm%rank/npg
      key   = mod(global_comm%rank,npg)

      ! A temporary comm is created to get the active processes
      call global_comm%split(color, key, temp_comm, ierr)  
      active_id = (temp_comm%rank == 0) 
     
      ! now the poisson comm is created as the intercomm
      ! things works because only only active_id(s) will do   
      call global_comm%split(key, color, poiss_comm, ierr)  
      
      id = poiss_comm%rank
      numprocs = poiss_comm%size

  end subroutine poiss_mpi_split

#:else
  
  logical, parameter :: id0 = .true.
  integer, parameter :: id = 0
  integer, parameter :: numprocs = 1
  logical, parameter :: active_id = .true.
         
#:endif


end module mpi_poisson
