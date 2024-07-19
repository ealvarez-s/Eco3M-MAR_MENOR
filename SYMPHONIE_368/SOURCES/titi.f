!...........................................................
! Declaration pour la subroutine my_outputs_zone1bioflux
      double precision ::        &
      ,zone1bioflux_w    =0.    &
      ,zone1waterflux_w   =0.

      integer ::              &
       zone1_nlayer   =0      &
      ,zone1_max      =0      &
      ,zone1_u_max    =0      &
      ,zone1_v_max    =0       

      real*4 ::                        &
             zone1_inv_dz      =0.01   &
            ,zone1_stretch_dz  =1.

      double precision , dimension(:,:) , allocatable ::    &
      ,zone1bioflux_glb      &
      ,zone1bioflux_u        &
      ,zone1bioflux_v        &
      ,zone1waterflux_glb     &
      ,zone1waterflux_u       &
      ,zone1waterflux_v        

      double precision , dimension(:,:) , allocatable ::    &
      ,zone1bioflux_glb_in      &
      ,zone1bioflux_u_in        &
      ,zone1bioflux_v_in        &
      ,zone1waterflux_glb_in     &
      ,zone1waterflux_u_in       &
      ,zone1waterflux_v_in        

      double precision , dimension(:,:) , allocatable ::    &
      ,zone1bioflux_glb_out      &
      ,zone1bioflux_u_out        &
      ,zone1bioflux_v_out        &
      ,zone1waterflux_glb_out     &
      ,zone1waterflux_u_out       &
      ,zone1waterflux_v_out   

      double precision , dimension(:) , allocatable ::      &
       zone1biocumul_glb      &
      ,zone1biocumul_loc      &
      ,zone1biomasst0         &
      ,zone1watercumul_glb     &
      ,zone1watercumul_loc     &
      ,zone1watermasst0         

      integer , dimension(:,:) , allocatable ::      &
       zone1_mask            &
      ,zone1_flux_u_node     &
      ,zone1_flux_v_node
!...........................................................
