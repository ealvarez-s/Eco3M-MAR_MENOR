&notebook_bathy

!...............
! Land-Sea mask and (TEMPORARY) H file:
! NOTE: the bathymetry contained in the following file will be smoothed by the model so 
! don't consider it as the true bathymetry
! texte250='../../../GLOBMED/BATHYMASK/bathycote_GLOBMED_curvigibraltar_1120_865.ijh'
! texte250='../../../GLOBMED/BATHYMASK/bathycote_GLOBMED_curvigibraltar_dardanelles_1120_865.ijh'
! texte250='../../../GLOBMED/BATHYMASK/bathy_GLOBMED_curvgib_dard_reparCors2fois_Sicil_1120_865.ijh'
!  texte250='../../../GLOBMED/BATHYMASK/bathy_GLOBMED_curvgib_dard_reparCors2fois_Sicil_Gib_dejaliss.ijh'
   texte250='../../../GLOBMED/BATHYMASK/bathy_GLOBMED_curvgib_dardnew_reparCors2fois_Sicil_Gib_dejaliss.ijh'
! texte250='none' ! if no file h=h1d
! h1d=100.        ! if no file h=h1d

! Remove secondary bassins
flag_remove_secondary_bassins=1 ! Note: High cost if large grid

!...............
! Bounds and Smoothing if option=1
! Run1
!ioption=1           ! 1=apply smoothing and threshold options, 0 otherwise
!h_inf=1.            ! h>h_inf(meters) (if ioption==1)
!h_inf_obc=10.       ! as h_inf but within a 10 points OBC layer
!rmax_=0.15          ! smoothing parameter "rmax". Default value 1/kmax if rmax<0
!hrmax=20.           ! (m) prevents division by 0 in r=abs(delta_H)/max(2H,hrmax) < rmax
!nsmooth=10          ! loop number in the smoothing "regular" routine
!flag_smooth_h_mask=0 ! masked h involved in smoothing if 1 

! Run2
!ioption=0            ! 1=apply smoothing and threshold options, 0 otherwise
 ioption=1            ! 1=apply smoothing and threshold options, 0 otherwise
 h_inf=1.             ! h>h_inf(meters) (if ioption==1)
 h_inf_obc=10.        ! as h_inf but within a 10 points OBC layer
!rmax_=0.1            ! smoothing parameter "rmax". Default value 1/kmax if rmax<0
 rmax_=2.            ! smoothing parameter "rmax". Default value 1/kmax if rmax<0
 hrmax=20.            ! (m) prevents division by 0 in r=abs(delta_H)/max(2H,hrmax) < rmax
 nsmooth=50 ! 5       ! loop number in the smoothing "regular" routine 
 flag_smooth_h_mask=1 ! masked h involved in smoothing if 1 (no masked h influence if 0)

!...............
! Wetdry parameters:
wetdry_cst1=1.              ! h for cancellation of baroclinic velocities
wetdry_cst2=0.1             ! h for cancellation of barotropic velocities
wetdry_cst3=0.1             ! h for cancellation of traceur tendencies 
mangrove_file_name='none'   ! magrove sea mask file
coef_diss_mangrove=0.1      ! Coef of the linear mangrove friction
linear_coef_mangrove=1      ! 0 for quadratic law or 1 for linear coef in mangrove

!...............
! Merge the bathymetry with the external ogcm bathy within a boundary sponge layer
  mergebathy_filename='nomerge' ! if no merge...
! mergebathy_filename='../../../GLOBMED/samplemercator/ext-PSY4V1R3_mesh_zgr.nc'
! Size of the boundary layer where the merge is applied:
  mergebathy_sponge=0  ! no merge
! mergebathy_sponge=30 ! merge over 30 grid nodes
!...............
/
&notebook_bathy2
! post smoothed bathymetry corrector file
! texte250='none'
  texte250='../../../GLOBMED/BATHYMASK/post_smoothed_bathy_corrector.ijDeltaH'
/
