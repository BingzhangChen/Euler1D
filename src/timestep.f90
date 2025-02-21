SUBROUTINE TIMESTEP
use forcing
use Time_setting
use state_variables
use grid
use NETCDF_IO
use Trait_functions,  only : PHY_C2Vol
implicit none

real,    parameter  :: cnpar      = 0.6d0
real,    parameter  :: Taur(nlev) = 1D12  !Relaxation time
real,    parameter  :: zero       = 0.d0  !Vectors of zero
real,    parameter  :: Vec0(nlev) = zero  !Vectors of zero
integer, parameter  :: mode0      = 0
integer, parameter  :: mode1      = 1
integer :: j = 0
real    :: par_save_freq = 0d0           !scratch variable for saving frequency of particles

!Benchmarking
real(4) :: dt1, dt2, dt3, dt4, t1, t2, t3, t4, t5

! 'START TIME STEPPING'
dt1 = 0.d0
dt2 = 0.d0
dt3 = 0.d0
dt4 = 0.d0

DO it = 1, Nstep+1

  call cpu_time(t1) 
  call update_time

  !For each time step, read in external environmental data%%%%%%%%%%%%%%%
  !Interpolate vertical profile of temperature at each timestep
  call time_interp(int(current_sec), N_time_temp, nlev, obs_time_temp, VTemp, Temp)

  !Calculate PAR
  call VERTICAL_LIGHT(current_DOY, sec_of_day, t(iChl,:))

  !Directly use the ROMS model output of Kv profiles
  call time_interp(int(current_sec), N_time_Kv, nlev+1, obs_time_Kv, VKv, Kv)

  !Start biology
  call cpu_time(t2) 
  dt1 = t2 - t1 + dt1 !The time for interpolating environmental data

  call BIOLOGY
  call cpu_time(t3) 
  dt2 = t3 - t2 + dt2 !The time for biology

  !Save the Eulerian output every day
  IF (mod(it, nsave) == 1) THEN

    !Update record 
    irec_Euler = irec_Euler + 1

    ! Add calculations of total nitrogen and save to Eulerian output files
    call Cal_total_N 

    write(6, 101) "Day", current_day, ": Total Nitrogen =", Ntot

    !Save data of Eulerian fields to the Euler.nc
    call write_Eulerfile(irec_Euler, current_day, current_hour)

  ENDIF

  call cpu_time(t4) 
  dt3 = t4 - t3 + dt3 !The time for saving data

  ! Diffusion
  Do j = 1,NVAR
     ! Zero flux at bottom
     call diff_center(nlev,dtsec,cnpar,1,Hz, Neumann, Neumann, &
                   zero, zero, Kv, Vec0,Vec0,Taur, t(j,:), t(j,:), t(j,:))
  Enddo

  ! Sinking:
  do j = 1,NVsinkterms
     SELECT CASE (bot_bound)
     case(Neumann)  ! closed at bottom (Conserve total N mass)
        call adv_center(nlev,dtsec,Hz,Hz,ww(:,j),1,1,zero,zero,    6,mode1,t(Windex(j),:))
     case(Dirichlet)
        ! Open bottom boundary
        call adv_center(nlev,dtsec,Hz,Hz,ww(:,j),1,2,zero,t(Windex(j),1),6,mode1,t(Windex(j),:))
     case default
        stop "The boundary conditions incorrect! STOP!"
     ENDSELECT
  enddo

  call cpu_time(t5) 
  dt4 = t5 - t4 + dt4  !The time for diffusion and sinking
ENDDO

print '("Environmental interpolation costs ",f8.3," hours.")', dt1/3600.0 
print '("Biology costs ",f8.3," hours.")', dt2/3600.0 
print '("Saving data costs ",f8.3," hours.")', dt3/3600.0 
print '("Diffusion and detritus sinking cost ",f8.3," hours.")', dt4/3600.0 

100 format(A4,I0,A3)
101 format(A3,1x, I0, 1x, A25, 1x, F10.4)
102 format(A5,I0,A3)
END SUBROUTINE TIMESTEP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Cal_total_N
use grid, only : nlev, Hz, Z_r
use state_variables, only : t, Ntot, iPHYC, iCHL, iPHYN, iZOO, iNO3, iDET, NZOO, NPHY
implicit none
integer :: k,i,j,m
real    :: Max_N = 0.d0

Ntot = 0d0
do k = 1, nlev
   do m = 1, NPHY
     if (t(iPHYN(m),k) .ne. t(iPHYN(m),k)) then
        write(6,*) "Phyto N", m, "is NaN at depth", Z_r(k)
        stop 
     endif

     if (t(iPHYN(m),k) < 0d0) then
        write(6,*) "Phyto N", m, "is negative at depth", Z_r(k)
        stop 
     endif
   enddo

   do m = 1, NZOO
     if (t(iZOO(m),k) .ne. t(iZOO(m),k)) then
        write(6,*) "ZOO", m, " is NaN at depth", Z_r(k)
        stop 
     endif

     if (t(iZOO(m),k) < 0d0) then
        write(6,*) "ZOO", m, " is negative at depth", Z_r(k)
        stop 
     endif
   enddo

   if (t(iDET,k) .ne. t(iDET,k)) then
      write(6,*) "DET is NaN at depth", Z_r(k)
      stop 
   endif

   if (t(iDET,k) < 0d0) then
      write(6,*) "DET is negative at depth", Z_r(k)
      stop 
   endif

   if (t(iNO3,k) .ne. t(iNO3,k)) then
      write(6,*) "NO3 is NaN at depth", Z_r(k)
      stop 
   endif

   if (t(iNO3,k) < 0d0) then
      write(6,*) "NO3 is negative at depth", Z_r(k)
      stop 
   endif

   !Update total N
   Ntot = Ntot + Hz(k)*(t(iNO3, k) + t(iDET, k))

   !Add total ZOOplankton N into total N
   do i = 1, NZOO
      Ntot = Ntot + t(iZOO(i), k) * Hz(k)
   enddo

   !Add total phytoplankton N into total N
   do i = 1, NPHY
      Ntot = Ntot + t(iPHYN(i),k) * Hz(k)
   enddo

enddo
END subroutine Cal_total_N