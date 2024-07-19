&notebook_time
! https://docs.google.com/document/d/1og-v_nRLtCZ13-Qs0NrcNJZaCqwhHwwtDREZdP2dbgI/edit

! Enter the time for the start and the end of the simulation:
datesim(1:6,1)= 2009 , 03 , 01 , 00  , 00  , 00  ! Start time (yyyy mm dd hh mm ss)
datesim(1:6,2)= 2009 , 06 , 09 , 00  , 00  , 00  ! End   time (yyyy mm dd hh mm ss)

! Define datesim(1:6,2) from a maximum number of iterations of the internal mode:
iteration3d_max=-999 ! active if > 0 . Number of iterations from datesim(1:6,1) 

!-----------------------------------------------------------------------------------
initial=0                      ! 0> Standard procedure  >1 Start from the end the previous simulation

! Restart files (relevant if initial=1):
restartfileperiod=100.             ! Production periodicity (see restartfileunits) of the restart files (no file if <=0.)
restartfileunits='d'               ! units of restartfileperiod. 'd'=days, 'm'=minutes, 's'=seconds
flag_0status_option=0              ! 0=stop if field not found, 1=ignore and go on if field not found
flag_w_binary=1                    ! write  0= netcdf restart single file. 1= binary chanel9 restart files (one per mpi process)
flag_r_binary=1                    ! read   0= netcdf restart single file. 1= binary chanel9 restart files (one per mpi process)
restartdir_out1='restart_output/'  ! Directory name for restart output files
restartdir_out2='restart_outbis/'  ! Directory name for restart output files
restartdir_in='restart_input/'     ! Directory name for restart input files

! simple restart file (but initial=0)
simple_restart_file_txt='none'    ! simple restart file if not "none" (following link for details)
simple_restart_biofile_txt='none' ! same but for passive tracers
! https://docs.google.com/document/d/1og-v_nRLtCZ13-Qs0NrcNJZaCqwhHwwtDREZdP2dbgI/edit
!-----------------------------------------------------------------------------------
! Time step: 
! External mode: set the CFL input parameters of the barotropic momentum equation
   cfl_sshmax=3.           ! A priori maximum value for the sea surface height (m)  
   cfl_hsshmax=3.          ! A priori maximum value for H+SSH (m)
   cfl_umax=3.             ! A priori maximum value for the current (m/s)             
   cfl_reduce=1.8          ! A priori attenuation factor for the CFL first guess

! Internal mode:
! Option 1: give the ratio dt_internal/dt_external
   iteration2d_max_now=16  ! ratio dt_internal/dt_external   
! Option 2: give the internal step (iteration2d_max_now will be adjusted)
!  dti_fw=180.             ! internal time step in seconds

! dti_fw variable with time using the following file (if not 'none')
 variable_time_step_txt='none'
! details in:
! https://docs.google.com/document/d/1og-v_nRLtCZ13-Qs0NrcNJZaCqwhHwwtDREZdP2dbgI/edit

! ECO3M-S time step is dti_fw*modulo_biotimestep (ex: 1h if dti_fw=180s and modulo_biotimestep=20)
  modulo_biotimestep=20

!-----------------------------------------------------------------------------------
! SPINUP
  spinup_forcing=86400.    ! time in second of ramping up of certain forcings (tides, rivers, swell)

!-------------------------------------------------------------------------------------
! Simulation mode:
run_option=0    ! 0=Standard procedure  -1=Stop the simulation at the end of the initialization

/

Note about RESTARTFILEPERIOD:
If RESTARTFILEPERIOD=0 or RESTARTFILEPERIOD<0 no restart file
Otherwise RESTARTFILEPERIOD indicates the periodicity (in days) of the restart file production.
If RESTARTFILEPERIOD>=0 a restart file is created at the end of the simulation.
