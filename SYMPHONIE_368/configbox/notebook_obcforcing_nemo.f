&notebook_obcforcing
!*******************************************************************************
!Large scale OGCM files in the Open Boundary Conditions:
!*******************************************************************************

! Nouveautes:
! v255: possibilite d'enchainer des grilles nemo differentes 
! Details dans:
! https://docs.google.com/document/d/1VMmzC7jXV7sycAxoFCUcOBxVbydbyXD3imgM0wAnAFQ/edit

!Set the following parameter to 1 to add the OGCM forcing (0 otherwise):
iobc_ogcm=0            !  1=on  0=off                                       IOBC_AF
obc_option=2           !  2="online interpolation"

! OGCM Name:
obc_ogcm_type='nemo_z' !  OGCM key

! Variables files list:
obcfile(1,1)='../../../REGION/LIST/list_var_T'
obcfile(2,1)='../../../REGION/LIST/list_var_S'
obcfile(3,1)='../../../REGION/LIST/list_var_U'
obcfile(4,1)='../../../REGION/LIST/list_var_V'
obcfile(5,1)='../../../REGION/LIST/list_var_SSH'

! Grid files list:
obcfile(6,1)='../../../REGION/LIST/list_grid_T'
obcfile(7,1)='../../../REGION/LIST/list_grid_U'
obcfile(8,1)='../../../REGION/LIST/list_grid_V'

! Add the expected effect of the tidal mixing on the temperature and salinity fields 
! if NEMO was not forced by the tide.
flag_ogcmtidemixing=1 !on=1, off=0

!Add the Inverse Barometer if not contained in the ssh ogcm fied:
bi_onoff=1  ! 1=add the BI   0=otherwise   

! Files conventions:
ogcm_time_shift=0  ! netcdf time is the middle (0) or the end (1) of the sampling period
ogcm_time_lag=0.   ! Time delay (seconds) to be added to the OGCM time
obctime_order=2    ! 2=linear time interpolation, 4=3rd degree polynomial time interpolation

! Value (in meters) of the offset to add to SSH in obc files 
offset_sshobc=0.

/
