&notebook_obcforcing
!*******************************************************************************
!Large scale OGCM files in the Open Boundary Conditions:
!*******************************************************************************
!Set the following parameter to 1 to add the OGCM forcing (0 otherwise):
iobc_ogcm=0            !  1=on  0=off                                       IOBC_AF
obc_option=2           !  2="online interpolation"

! OGCM Name:
obc_ogcm_type='nemo_z' !  OGCM key

! Variables files list:
obcfile(1,1)='../../../WM/LIST_OGCM/list_var_T'
obcfile(2,1)='../../../WM/LIST_OGCM/list_var_S'
obcfile(3,1)='../../../WM/LIST_OGCM/list_var_U'
obcfile(4,1)='../../../WM/LIST_OGCM/list_var_V'
obcfile(5,1)='../../../WM/LIST_OGCM/list_var_SSH'

! Grid files list:
obcfile(6,1)='../../../WM/LIST_OGCM/list_grid_T'
obcfile(7,1)='../../../WM/LIST_OGCM/list_grid_U'
obcfile(8,1)='../../../WM/LIST_OGCM/list_grid_V'

! netcdf file with annual average of the very same variables (and same ogcm grid):
 obcfile(9,1)='none'
!obcfile(9,1)='quickinitial'
!obcfile(9,1)='../../../WM/OGCM_ANNUAL/saltem2012.nc'

!Add the Inverse Barometer if not contained in the ssh ogcm fied:
bi_onoff=1  ! 1=add the BI   0=otherwise   

! Files conventions:
ogcm_time_shift=0  ! netcdf time is the middle (0) or the end (1) of the sampling period

! Value of the offset to add to SSH in obc files
offset_sshobc=0.

/
