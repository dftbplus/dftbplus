
#:include 'common.fypp'

module mpi_poisson

#:if WITH_MPI
  use mpi
  use libmpifx_module
  implicit none
  private  

  !! mpifx communicator as received by calling program
  TYPE(mpifx_comm), PUBLIC :: global_comm
  !! Local communicator (obtained by splitting), if any
  TYPE(mpifx_comm), PUBLIC :: poiss_comm
  INTEGER, PUBLIC ::  numprocs
  INTEGER, PUBLIC ::  id
  LOGICAL, PUBLIC ::  id0
  logical, public :: active_id
  !INTEGER, PUBLIC :: poiss_comm 
  PUBLIC :: MPI_DOUBLE_PRECISION
  PUBLIC :: MPI_IN_PLACE
  PUBLIC :: MPI_DATATYPE_NULL
  PUBLIC :: MPI_SUM
  PUBLIC :: poiss_mpi_init
  PUBLIC :: poiss_mpi_split

  contains

  subroutine  poiss_mpi_init(comm)
      type(mpifx_comm), intent(in) :: comm

      global_comm = comm 
      !! If no split occurs, then the local and the global communicator 
      !! must be the same
      poiss_comm = comm
      id = comm%rank
      numprocs = comm%size
      ! It doesn't make difference (in theory) if we define
      ! id0 here or in the split as only one communicator
      ! ranks<maxnumprocs is active during poisson
      id0 = (id == 0)
      
  end subroutine poiss_mpi_init

  subroutine poiss_mpi_split(maxnumprocs)
      integer, intent(in) :: maxnumprocs

      integer :: irow, jcol, ierr
           
      irow = global_comm%rank/maxnumprocs
      jcol = mod(global_comm%rank,maxnumprocs)

      !! Only the rank in the first sub communicator are active
      !! This could maybe achieved in a more elegant way using 
      !! different groups and communicators for idle and active ids
      if (global_comm%rank < maxnumprocs) then
        active_id = .true.
      else 
        active_id = .false.
      endif
      
      call global_comm%split(irow,jcol,poiss_comm,ierr)  
      id = poiss_comm%rank
      numprocs = maxnumprocs
      !! id0 does not change, as we mean id=0 for first communicator only

  end subroutine poiss_mpi_split

#:else
  
  logical, parameter :: id0 = .true.
  integer, parameter :: id = 0, numprocs = 1
  logical, parameter :: active_id = .true.
         
#:endif


end module mpi_poisson
