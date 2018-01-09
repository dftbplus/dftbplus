module solvertypes
  implicit none     
  !> Solver types (used like enumerators)

  integer, parameter :: mixerSimple = 1 
  integer, parameter :: mixerAnderson = 2   
  integer, parameter :: mixerBroyden = 3  
  integer, parameter :: mixerDIIS = 4

  integer, parameter :: optSD = 1 
  integer, parameter :: optCG = 2   
  integer, parameter :: optDIIS = 3  
   
  integer, parameter :: forceOrig = 0 
  integer, parameter :: forceRediag = 1 
  integer, parameter :: forceDynT0 = 2   
  integer, parameter :: forceDynT = 3  

  integer, parameter :: solverQR = 1
  integer, parameter :: solverDAC = 2
  integer, parameter :: solverRR1 = 3
  integer, parameter :: solverRR2 = 4
  integer, parameter :: solverGF = 5
  integer, parameter :: onlyTransport = 6


  integer, parameter :: gammaf = 0 
  integer, parameter :: poisson = 1
 
end module solvertypes
