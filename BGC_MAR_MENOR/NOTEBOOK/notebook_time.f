&notebook_time

! Enter the time for the start and the end of the simulation:
datesim(1:6,1)= 2016 , 09 , 01 , 12  , 00  , 00  ! Start time (yyyy mm dd hh mm ss)
datesim(1:6,2)= 2016 , 09 , 15 , 12  , 00  , 00  ! End time (yyyy mm dd hh mm ss)

! Define datesim(1:6,2) from a maximum number of iterations of the internal mode:
 iteration3d_max=-999 ! active if > 0 . Number of iterations from datesim(1:6,1)
!iteration3d_max=100 ! active if > 0 . Number of iterations from datesim(1:6,1)

!-----------------------------------------------------------------------------------
! Initial state:

! Attention, si tu changes des options qui entrainent l'allocation dynamique de
! nouveaux tableaux (par ex tu changes le schema de turbulence) utiliser les
! restart impose que tu crees d'abord un premier jeux de fichiers netcdf en mode "grid
! full". Vider les repertoire restart avant operation.
initial=0               ! 0> Standard procedure  >1 Start from the end the previous simulation
restartfileperiod=30.    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXProduction periodicity (days) of the restart files (no file if <=0.)

!-----------------------------------------------------------------------------------
! Time step: 
! External mode: set the CFL input parameters of the barotropic momentum equation
   cfl_sshmax=3.           ! A priori maximum value for the sea surface height (m)  
   cfl_umax=3.             ! A priori maximum value for the current (m/s)             
   cfl_reduce=1.8  !1.5        ! A priori attenuation factor for the CFL first guess

! Internal mode:
! Option 1: give the ratio dt_internal/dt_external
   iteration2d_max_now=16 ! 30  ! ratio dt_internal/dt_external   
! Option 2: give the internal step (iteration2d_max_now will be adjusted)
!  dti_fw=180.             ! internal time step in seconds
   dti_fw=40.              ! internal time step in seconds

! time stepping method for T,S, & baroclinic current:
timestep_type=1            ! 0=Leap-Frog 1=Forward-Backward

!-------------------------------------------------------------------------------------
! Simulation mode:
run_option=0   ! 0=Standard procedure  -1=Stop the simulation at the end of the initialization
! penser a remettre restartfileperiod>0

/

Note about RESTARTFILEPERIOD:
If RESTARTFILEPERIOD=0 or RESTARTFILEPERIOD<0 no restart file
Otherwise RESTARTFILEPERIOD indicates the periodicity (in days) of the restart file production.
If RESTARTFILEPERIOD>=0 a restart file is created at the end of the simulation.
