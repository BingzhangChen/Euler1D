PROGRAM Euler1D
USE GRID
IMPLICIT NONE
real(4) :: start,finish

!Count time
call cpu_time(start) 

!Start the pseudo random number generator
call random_seed()

!Define grid
call setup_grid

!Initialization
call initialize

!Timestep
call Timestep

call cpu_time(finish)
Print '("Whole simulation time = ",f8.3," hours.")', (finish-start)/3600.0 

END program