MODULE NETCDF_IO
USE NETCDF
USE Grid,            only : nlev,  Z_r, Z_w
USE State_variables, only : t, NPHY, NZOO, iNO3, iPHYC, iPHYN, iCHL, iZOO, iDET
USE State_variables, only : Nout, Varout, oNPP, oTEMP, oPAR
USE Time_setting,    only : it, Nstep
IMPLICIT NONE

private

public :: irec_Euler 
public :: create_Eulerian_file, write_Eulerfile

character (len=8),  public     ::  Euler_FNAME = 'Euler.nc'
character (len=3),  parameter  ::  Zr_NAME = 'Z_r'
character (len=3),  parameter  ::  Zw_NAME = 'Z_w'
character (len=4),  parameter  ::  hr_NAME = 'Hour'
character (len=3),  parameter  ::  DAY_NAME = 'Day'
character (len=3),  parameter  ::  DOY_NAME = 'DOY'
character (len=4),  parameter  ::  timestep_NAME = 'Step'   !For restart file
character (len=3),  parameter  ::  NO3_NAME   = 'NO3'
character (len=4),  parameter  ::  PC_NAME    = 'PHYC'
character (len=4),  parameter  ::  PN_NAME    = 'PHYN'
character (len=3),  parameter  ::  CHL_NAME   = 'CHL'
character (len=3),  parameter  ::  ZOO_NAME   = 'ZOO'
character (len=3),  parameter  ::  DET_NAME   = 'DET'
character (len=3),  parameter  ::  NPP_NAME   = 'NPP'
character (len=2),  parameter  ::  Kv_NAME    = 'Kv'
character (len=4),  parameter  ::  Temp_NAME  = 'Temp'
character (len=3),  parameter  ::  PAR_NAME   = 'PAR'
character (len=8),  parameter  ::  FZ_NAME    = 'FZ'

character (len=5),  parameter  ::  UNITS = 'units'

INTEGER :: Temp_varid, PAR_varid, Kv_varid, NPP_varid
INTEGER :: Zr_varid, Zw_varid, step_varid, DAY_varid, DOY_varid, hr_varid
INTEGER :: alpha_varid, CHL_varid, DET_varid, NO3_varid, PC_varid
INTEGER :: PN_varid, ZOO_varid

!Records of Euler
integer :: irec_Euler = 0

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Create Eulerian file for saving Eulerian outputs
SUBROUTINE create_Eulerian_file
IMPLICIT NONE
CHARACTER (LEN=10), PARAMETER :: REC_NAME = 'Time'
CHARACTER (LEN=10), PARAMETER :: UNIT_dist = 'm'

!Number of Dimensions for normal tracers
integer,  parameter :: NDIMS    = 2
integer,  parameter :: NDIM_PHY = 3 !Three traits, time, and depth
integer,  parameter :: NDIM_ZOO = 3
character (len=8)   :: date

INTEGER :: ncid, rec_dimid, NZr_dimid, NZw_dimid, NZOO_dimid, NPHY_dimid
INTEGER :: dimid_r(NDIMS), dimid_w(NDIMS), dimid_zoo(NDIM_ZOO), dimid_PHY(NDIM_PHY)

! Create Euler file (nc file)
CALL check(nf90_create(Euler_FNAME, nf90_clobber, ncid) )
  
! Define the dimensions. The record dimension is defined to have
! unlimited length - it can grow as needed.
CALL check(nf90_def_dim(ncid, Zr_NAME,  nlev,   NZr_dimid))
CALL check(nf90_def_dim(ncid, Zw_NAME,  nlev+1, NZw_dimid))
CALL check(nf90_def_dim(ncid, ZOO_NAME, NZOO,   NZOO_dimid))
CALL check(nf90_def_dim(ncid, PC_NAME,  NPHY,   NPHY_dimid))
CALL check(nf90_def_dim(ncid, PN_NAME,  NPHY,   NPHY_dimid))
CALL check(nf90_def_dim(ncid, CHL_NAME, NPHY,   NPHY_dimid))
CALL check(nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )

! The dimids array is used to pass the dimids of the dimensions of the netCDF variables. 
! In Fortran, the unlimited dimension must come last on the list of dimids.
dimid_r   = (/ NZr_dimid, rec_dimid /)
dimid_w   = (/ NZw_dimid, rec_dimid /)
dimid_zoo = (/NZOO_dimid, NZr_dimid, rec_dimid /)
dimid_PHY = (/NPHY_dimid, NZr_dimid, rec_dimid /)

! Define the variables of Eulerian fields
CALL check(nf90_def_var(ncid, Zr_NAME,  NF90_REAL, NZr_dimid, Zr_varid) )
CALL check(nf90_def_var(ncid, Zw_NAME,  NF90_REAL, NZw_dimid, Zw_varid) )
CALL check(nf90_def_var(ncid, DAY_NAME, NF90_INT,  rec_dimid, DAY_varid) )
CALL check(nf90_def_var(ncid, hr_NAME,  NF90_INT,  rec_dimid, hr_varid) )

! Define the netCDF variables for variables that need to be saved into the Euler.nc
CALL check( nf90_def_var(ncid, Temp_NAME,   NF90_REAL, dimid_r,   Temp_varid) )
CALL check( nf90_def_var(ncid, PAR_NAME,    NF90_REAL, dimid_r,   PAR_varid)  )
CALL check( nf90_def_var(ncid, Kv_NAME,     NF90_REAL, dimid_w,   Kv_varid)   )
CALL check( nf90_def_var(ncid, NPP_NAME,    NF90_REAL, dimid_r,   NPP_varid)  )
CALL check( nf90_def_var(ncid, NO3_NAME,    NF90_REAL, dimid_r,   NO3_varid)  )
CALL check( nf90_def_var(ncid, DET_NAME,    NF90_REAL, dimid_r,   DET_varid)  )
CALL check( nf90_def_var(ncid, PC_NAME,     NF90_REAL, dimid_PHY, PC_varid)   )
CALL check( nf90_def_var(ncid, PN_NAME,     NF90_REAL, dimid_PHY, PN_varid)   )
CALL check( nf90_def_var(ncid, CHL_NAME,    NF90_REAL, dimid_PHY, CHL_varid)  )
CALL check( nf90_def_var(ncid, ZOO_NAME,    NF90_REAL, dimid_zoo, ZOO_varid)  )

! Assign units attributes to the netCDF variables.
CALL date_and_time(DATE=date)
CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Date', date))
CALL check( nf90_put_att(ncid,   Zr_varid, UNITS, UNIT_dist))
CALL check( nf90_put_att(ncid,   Zw_varid, UNITS, UNIT_dist))
CALL check( nf90_put_att(ncid, DAY_varid, UNITS, 'days'))
CALL check( nf90_put_att(ncid, hr_varid, UNITS, 'hour'))
CALL check( nf90_put_att(ncid, NPP_varid, UNITS, 'mg C m-3 d-1'))
CALL check( nf90_put_att(ncid, Kv_varid, UNITS, 'm^2 s-1'))
CALL check( nf90_put_att(ncid, PAR_varid, UNITS, 'W m-2'))
CALL check( nf90_put_att(ncid, Temp_varid, UNITS, 'ºC'))
CALL check( nf90_put_att(ncid, NO3_varid, UNITS, 'mmol N m-3'))
CALL check( nf90_put_att(ncid, PN_varid, UNITS, 'mmol N m-3'))
CALL check( nf90_put_att(ncid, PC_varid, UNITS, 'mmol C m-3'))
CALL check( nf90_put_att(ncid, CHL_varid, UNITS, 'mg Chl m-3'))
CALL check( nf90_put_att(ncid, DET_varid, UNITS, 'mmol N m-3'))
CALL check( nf90_put_att(ncid, ZOO_varid, UNITS, 'mmol N m-3'))

! End define mode.
CALL check( nf90_enddef(ncid) )
  
! Write the Z coordinate data
CALL check( nf90_put_var(ncid, Zr_varid, Z_r) )
CALL check( nf90_put_var(ncid, Zw_varid, Z_w) )
CALL check( nf90_close(ncid) )
return
END SUBROUTINE create_Eulerian_file

SUBROUTINE write_Eulerfile(rec, day, hour)
use Forcing, only : Kv
IMPLICIT NONE
INTEGER, INTENT(IN)  :: rec  !The time index to be written
INTEGER, INTENT(IN)  :: day 
INTEGER, INTENT(IN)  :: hour 
INTEGER              :: i, j
INTEGER              :: ncid = 0
real ::  cff(NZOO, nlev) = 0.
real :: cff1(NPHY, nlev) = 0.

!Open the nc file for writing
CALL check(NF90_OPEN(Euler_FNAME, NF90_WRITE, ncid))

CALL check(NF90_INQ_VARID(ncid, Kv_NAME,  Kv_varid))     ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, Temp_NAME,Temp_varid))   ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, DAY_NAME, DAY_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, Hr_NAME,  Hr_varid))     ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, PAR_NAME, PAR_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, NO3_NAME, NO3_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, PC_NAME, PC_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, PN_NAME, PN_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, CHL_NAME, CHL_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, DET_NAME, DET_varid))    ! get variable IDs
CALL check(NF90_INQ_VARID(ncid, ZOO_NAME, ZOO_varid))    ! get variable IDs

!Add data into the Euler.nc
CALL check(NF90_PUT_VAR(ncid, DAY_varid, day, start = (/rec/)))
CALL check(NF90_PUT_VAR(ncid, Hr_varid, hour, start = (/rec/)))

CALL check(NF90_PUT_VAR(ncid, Kv_varid, Kv(:), start=[1,rec],           &
                                               count=[1+nlev,1]))

CALL check(NF90_PUT_VAR(ncid, PAR_varid, Varout(oPAR,:),start=[1,rec],  &
                                                        count=[nlev,1]))

CALL check(NF90_PUT_VAR(ncid, Temp_varid,Varout(oTemp,:),start=[1,rec],  &
                                                         count=[nlev,1]))
 
CALL check(NF90_PUT_VAR(ncid, DET_varid, t(iDET,:),start=[1,rec],  &
                                                   count=[nlev,1]))
CALL check(NF90_PUT_VAR(ncid, NO3_varid, t(iNO3,:),start=[1,rec],  &
                                                   count=[nlev,1]))

!save phytoplankton into a temporary matrix                                                   
do i = 1, NPHY                                                
  cff1(i,:) = t(iPHYC(i),:)
enddo

CALL check(NF90_PUT_VAR(ncid, PC_varid, cff1,  start=[1,   1,   rec],  &
                                               count=[NPHY,nlev,1  ]))

do i = 1, NPHY                                                
  cff1(i,:) = t(iPHYN(i),:)
enddo

CALL check(NF90_PUT_VAR(ncid, PN_varid, cff1,  start=[1,1,rec],  &
                                               count=[NPHY,nlev,1]))

do i = 1, NPHY                                                
  cff1(i,:) = t(iCHL(i),:)
enddo

CALL check(NF90_PUT_VAR(ncid, CHL_varid, cff1,start=[1,1,rec],  &
                                              count=[NPHY,nlev,1]))

!save zooplankton into a temporary matrix                                                   
do i = 1, NZOO                                                
  cff(i,:) = t(iZOO(i),:)
enddo

CALL check(NF90_PUT_VAR(ncid, ZOO_varid, cff, start=[1,1,rec],  &
                                              count=[NZOO,nlev,1]))

CALL check(NF90_PUT_VAR(ncid, NPP_varid, Varout(oNPP,:),start=[1,rec],  &
                                                        count=[nlev,1]))

! Close the file. This causes netCDF to flush all buffers and make
! sure your data are really written to disk.
CALL check(nf90_close(ncid))

RETURN
END SUBROUTINE write_Eulerfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check(status)
integer, intent (in) :: status

if(status /= nf90_noerr) then
  print *, trim(nf90_strerror(status))
  stop "Stopped"
end if
end subroutine check

END MODULE
