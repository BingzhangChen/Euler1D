!This module declares the state variables for phytoplankton in Eulerian fields
MODULE STATE_VARIABLES
use grid, only : nlev, Hz
IMPLICIT NONE

!State variables
integer, private   :: i
integer, parameter :: NZOO    = 1  !Number of zooplankton size classes
integer, parameter :: N_ESD   = 1  !Number of phytoplankton size classes
integer, parameter :: N_Topt  = 1   !Number of phytoplankton optimal temperature classes
integer, parameter :: N_alpha = 1   !Number of phytoplankton classes for alphaI (light)
integer, parameter :: NPHY = N_ESD * N_Topt * N_alpha  !Total Number of phytoplankton species

!Indexes for the state variables
integer, parameter :: iNO3 = 1 
integer, parameter :: iPHYC(NPHY) = [(iNO3 + i ,        i = 1, NPHY)]
integer, parameter :: iPHYN(NPHY) = [(iPHYC(NPHY) + i , i = 1, NPHY)]
integer, parameter :: iCHL(NPHY)  = [(iPHYN(NPHY) + i , i = 1, NPHY)]
integer, parameter :: iZOO(NZOO)  = [(iCHL(NPHY)  + i , i = 1, NZOO)]
integer, parameter :: iDET = iZOO(NZOO) + 1 
integer, parameter :: nvar = iDET !Total number of state variables
real               :: t(nvar, nlev) = 0.d0
real               :: Ntot = 0d0  !Total nitrogen in the domain
real, parameter :: MinSzoo = log(0.8d0)!Minimal zooplankton log ESD (micron)
real, parameter :: MaxSzoo = log(36d2) !Maximal zooplankton log ESD (micron)
real, parameter :: dZOOESD = (MaxSzoo - MinSzoo)/dble(NZOO)! ESD difference between adjacent zooplankton size class (log)
real, parameter :: MinSPHY = log(0.6d0)!Minimal phytoplankton log ESD (micron)
real, parameter :: MaxSPHY = log(2d2) !Maximal PHYplankton log ESD (micron)
real, parameter :: dPHYESD = (MaxSPHY - MinSPHY)/dble(NPHY) ! ESD difference between adjacent PHYplankton size class (log)
real, parameter :: dPHYLnV = 3.d0 * dPHYESD

!Log ESD of each PHY. size class
real            :: ESDPHY(NPHY) = 0d0

!Volume of each PHYplankton size class
real            :: VolPHY(NPHY) = 0d0

!Log ESD of each zoo. size class
real            :: ESDZOO(NZOO) = 0d0

!Volume of each zooplankton size class
real            :: VolZOO(NZOO) = 0d0

! Define the number of sinking tracers:
integer, parameter :: NVsinkterms =  1 + 3*NPHY ! DET plus all phytoplankton

! Define the index of sinking tracers in the Vars matrix: 
integer, parameter :: Windex(NVsinkterms) = [(iPHYC(i), i = 1, NPHY), (iPHYN(i), i = 1, NPHY), (iCHL(i), i = 1, NPHY), iDET]

!Sinking rate
real               :: ww(0:nlev, NVsinkterms) = 0.d0

!Labels for each state variable (for saving model output)
integer, parameter :: oNPP = nvar + 1   !Daily net primary production integrated over a whole day
integer, parameter :: oTEMP= oNPP + 1
integer, parameter :: oPAR = oTEMP+ 1
integer, parameter :: Nout = oPAR
real               :: Varout(Nout, nlev) = 0.d0
character(LEN=7)   :: Labelout(Nout)     = 'Unknown'

!Model choices
integer, parameter :: GMK98_simple         = 1 
integer, parameter :: GMK98_Topt           = 2 
integer, parameter :: GMK98_Size           = 3 
integer, parameter :: GMK98_Light          = 4 
integer, parameter :: GMK98_ToptLight      = 5 
integer, parameter :: GMK98_ToptSize       = 6 
integer, parameter :: GMK98_SizeLight      = 7 
integer, parameter :: GMK98_ToptSizeLight  = 8 

!Current model selection
integer :: Model_ID = 1
real, parameter :: NO3_min = 0.01 !Mimimal NO3 concentration
END MODULE STATE_VARIABLES
