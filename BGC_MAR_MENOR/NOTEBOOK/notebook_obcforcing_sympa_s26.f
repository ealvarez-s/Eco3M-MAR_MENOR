&notebook_obcforcing
!*******************************************************************************
!Large scale OGCM files in the Open Boundary Conditions:
!*******************************************************************************
!Set the following parameter to 1 to add the OGCM forcing (0 otherwise):
iobc_ogcm=0            !  1=on  0=off                                       IOBC_AF

! OGCM Name:
obc_ogcm_type='sympa' !  OGCM key

! Variables files list:
obcfile(1,1)='../../../BGC_MAR_MENOR/LIST/obc/list_ogcm_uvts'
obcfile(2,1)='../../../BGC_MAR_MENOR/LIST/obc/list_ogcm_grid'

! Add the expected effect of the tidal mixing on the temperature and salinity fields 
! if SYMPHONIE was not forced by the tide.
flag_ogcmtidemixing=0 !on=1, off=0

!Add the Inverse Barometer if not contained in the ssh ogcm fied:
bi_onoff=1  ! 1=add the BI   0=otherwise   

! Files conventions:
ogcm_time_shift=0  ! netcdf time is the middle (0) or the end (1) of the sampling period
ogcm_time_lag=0.   ! Time delay (seconds) to be added to the OGCM time
obctime_order=2    ! 2=linear time interpolation, 4=3rd degree polynomial time interpolation

! anciens parametres, ne pas modifier:
obc_option=2           !  2="online interpolation"
/
