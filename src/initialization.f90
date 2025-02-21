SUBROUTINE INITIALIZE
USE params
USE state_variables
USE Time_setting
USE grid,             only: Z_w, hmax
USE Trait_functions,  only: PHY_ESD2C, PHY_C2Vol, phyto_sinking
USE NETCDF_IO
USE forcing
IMPLICIT NONE

integer    :: rc    = 0
integer    :: k     = 0
integer    :: i     = 0
integer    :: j     = 0
integer    :: ct    = 0
integer    :: j_    = 0
real       :: cff   = 0.d0
real       :: Z_avg = 0.d0
real       :: Vol   = 0.
real       :: QN    = 0.
real       :: NO3_Jan(nlev, 2) = 0.
integer, parameter :: namlst         = 20   !Unit time for namelist files
integer            :: AllocateStatus = 0
logical            :: exists         = .true.

character(LEN=20) :: NO3_file='BATS_NO3_Jan.dat'

!==========================================================

!Namelist definition of time settings
namelist /timelist/  NDay_Run, dtsec, nsave

!Namelist definition of model choice and parameters
namelist /paramlist/ Model_ID, mu0, aI0, KN, gmax, Kp, mz,GGE, unass, RDN, &
                     wDET, SDZoo, sigma

! Check whether the namelist file exists.
inquire (file='time.nml', iostat=rc)

if (rc /= 0) then
    write (6, '(a)') 'Error: namelist file time.nml does not exist.'
    stop
end if

!  open the namelist file and read initial paramter values
open(namlst,file='time.nml',status='old',action='read')
read(namlst,nml=timelist)
close(namlst)

!Total Number of time steps
Nstep = NDay_Run * INT(d_per_s)/INT(dtsec) 

!Calculate dtdays
dtdays = dtsec/d_per_s
write(6,'(A13,1x,1pe12.2,A12)') 'Timestepping: ', dtdays, 'of one day.'
write(6,'(A30,1x,I0)') 'Total number of simulation days: ', NDay_Run
!==========================================================
!Read parameter namelist
!Check whether the namelist file exists.
inquire (file='param.nml', exist = exists)

if (.not. exists) then
    write (6, '(a)') 'Error: namelist file Model.nml does not exist.'
    stop
end if

!  open the namelist file and read initial paramter values
open(namlst,file='param.nml',status='old',action='read')
read(namlst,nml=paramlist)
close(namlst)

WRITE(6,'(A15,1x,I1)') 'Select Model ID', Model_ID

!Calculate delta
delta = sigma**2 / 3.d0 / dPHYLnV**2

! Prepare forcing
! Prepare temperature for temporal interpolation
call extract_WOAtemp
  
! Prepare Kv
call extract_Kv

!Initialize ZOOplankton size
If (Model_ID .eq. GMK98_Size          .or. Model_ID .eq. GMK98_ToptSize .or.  &
    Model_ID .eq. GMK98_ToptSizeLight .or. Model_ID .eq. GMK98_SizeLight) then

  !Initialize phytoplankton size (logESD)
  ct = 0 !Count phytoplankton index
  DO k = 1, N_ESD
    do i = 1, N_Topt
      do j = 1, N_alpha

        ct = ct + 1
  
        ESDPHY(ct) = MinSPHY + dble(k-1)*dPHYESD
  
        !Write out the phytoplankton size
        write(6,1001) "PHY", ct, "ESD = ", exp(ESDPHY(ct)), "micron"
  
        !Compute volume of phytoplankton
        VolPHY(ct) = pi/6d0*exp(ESDPHY(ct))**3
      enddo
    enddo
  ENDDO

  !Initialize zooplankton size (logESD)
  do k = 1, NZOO
  
     ESDZOO(k) = MinSzoo + dble(k-1)*dZOOESD
  
     !Write out the zooplankton size
     write(6,1001) "ZOO", k, "ESD = ", exp(ESDZOO(k)), "micron"
  
     !Compute volume of zooplankton
     VolZOO(k) = pi/6d0*exp(ESDZOO(k))**3
  enddo
Else

  !Check if NZOO is consistent with Model_ID (i.e., without size classes, NZOO should be one)
  if(NZOO > 1) then
    stop "Number of zooplankton size classes should be ONE if size is not modelled!"
  endif

  ESDPHY(1) = 3.d0
  VolPHY(1) = pi/6d0*ESDPHY(1)**3  !Assume phytoplankton volume 3 micron
  ESDZOO(1) = 60.d0
  VolZOO(1) = pi/6d0*ESDZOO(1)**3  !Assume zooplankton volume 60 micron

Endif
 
! Initialize initial NO3 using WOA data:
! Read NO3 data that are already matched to the grid:
call Readcsv(NO3_file, nlev, 2, NO3_Jan)

do k = 1, nlev
  t(iNO3,k) = NO3_Jan(k,2)
enddo

do k = 1, NPHY
   t(iPHYN(k),:) = 0.1d0/dble(NPHY) !Assuming an initial condition of uniform biomass among different phyto. species
   t(iPHYC(k),:) = t(iPHYN(k),:)/16.d0*106.d0
   t(iCHL(k), :) = t(iPHYC(k),:)*12.d0/50.d0  !Unit: mgChl m-3
enddo
do k = 1, NZOO
   t(iZOO(k),:) = 0.1d0/dble(NZOO) !Assuming an initial condition of uniform biomass among different zoo. size classes
enddo

!Following Verity et al. AME (1996)
t(iDET,:) = .1d0

if (Model_ID .eq. GMK98_Topt          .or. Model_ID .eq. GMK98_ToptLight .or. &
    Model_ID .eq. GMK98_ToptSizeLight .or. Model_ID .eq. GMK98_ToptSize) then
   !Initialize phytoplankton optimal temperature (Topt) from a uniform distribution between 2 and 30 degree celcius
   !TODO
endif

if (Model_ID .eq. GMK98_Light         .or. Model_ID .eq. GMK98_ToptLight .or. &
    Model_ID .eq. GMK98_ToptSizeLight .or. Model_ID .eq. GMK98_SizeLight) then

    !Assign light traits
    !TODO
endif

!Initialize Varout
do k = 1, NVAR
   Varout(k,:) = t(k,:)
enddo

!Initialize time
it = 1
call update_time

!Save initial state to external file
call create_Eulerian_file
irec_Euler = 1
call write_Eulerfile(irec_Euler, current_day, current_hour)

!Sinking rate
!Initialize sinking rate (UNIT: m/s !):
ww(:,:) = 0d0
do k = 0,nlev-1
  !Phytoplankton sinking rate following Durante et al. JPR 2019 
  do i = 1, NPHY
    ww(k,i) = -phyto_sinking(VolPHY(i))/dble(d_per_s) 
  enddo

  !Detritus sinking rate (convert to UNIT: m/s)
  ww(k,NVsinkterms) = -wDET/dble(d_per_s) 
enddo
return

1001 format(A3, 1x, I0, 1x, A6, F8.2, 1x, A6)
END SUBROUTINE INITIALIZE