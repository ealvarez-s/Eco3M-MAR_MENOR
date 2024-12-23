&notebook_time

! Enter the time for the start and the end of the simulation:
!datesim(1:6,1)= 2013 , 06 , 01 , 00  , 00  , 00 ! Start time (yyyy mm dd hh mm ss)
datesim(1:6,1)= 2011 , 08 , 01 , 00  , 00  , 00 ! Start time (yyyy mm dd hh mmss)

!datesim(1:6,2)= 2012 , 1 , 1 , 00  , 00  , 00  ! Start time (yyyy mm dd hh mm ss)
 datesim(1:6,2)= 2016 , 08 , 01 , 00  , 00  , 00  ! Start time (yyyy mm dd hh mm ss)
!datesim(1:6,2)= 2019 , 03 , 30 , 00  , 00  , 00  ! Start time (yyyy mm dd hh mm ss)

! Define datesim(1:6,2) from a maximum number of iterations of the internal mode:
!iteration3d_max=10 ! active if > 0 . Number of iterations from datesim(1:6,1)
!iteration3d_max=100 ! active if > 0 . Number of iterations from datesim(1:6,1)

!-----------------------------------------------------------------------------------
! Initial state:

! Attention, si tu changes des options qui entrainent l'allocation dynamique de
! nouveaux tableaux (par ex tu changes le schema de turbulence) utiliser les
! restart impose que tu crees d'abord un premier jeux de fichiers netcdf en mode "grid
! full". Vider les repertoire restart avant operation.
initial=2               ! 0> Standard procedure  >1 Start from the end the previous simulation
restartfileperiod=30.   ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXProduction periodicity (days) of the restart files (no file if <=0.)
restartfileunits='d'    ! units of restartfileperiod. 'd'=days, 'm'=minutes, 's'=seconds
flag_0status_option=0   ! 0=stop if field not found, 1=ignore and go on if field not found
flag_binary=1           ! 0= netcdf restart files. 1= old chanel9 restart files


!-----------------------------------------------------------------------------------
! Time step: 
! External mode: set the CFL input parameters of the barotropic momentum equation
   cfl_sshmax=3.           ! A priori maximum value for the sea surface height (m)  
   cfl_umax=3.             ! A priori maximum value for the current (m/s)             
   cfl_reduce=1.8  !1.5        ! A priori attenuation factor for the CFL first guess

! Internal mode:
! Option 1: give the ratio dt_internal/dt_external
   iteration2d_max_now=150 ! 26  ! ratio dt_internal/dt_external   
! Option 2: give the internal step (iteration2d_max_now will be adjusted)
!  dti_fw=180.             ! internal time step in seconds

! time stepping method for T,S, & baroclinic current:
timestep_type=1            ! 0=Leap-Frog 1=Forward-Backward

!-------------------------------------------------------------------------------------
! Simulation mode:
run_option=0  ! 0=Standard procedure  -1=Stop the simulation at the end of the initialization
! penser a remettre restartfileperiod>0

/

Note about RESTARTFILEPERIOD:
If RESTARTFILEPERIOD=0 or RESTARTFILEPERIOD<0 no restart file
Otherwise RESTARTFILEPERIOD indicates the periodicity (in days) of the restart file production.
If RESTARTFILEPERIOD>=0 a restart file is created at the end of the simulation.
