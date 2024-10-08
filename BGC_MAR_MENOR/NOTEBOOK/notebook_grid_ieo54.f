&notebook_grid

kmax=20           ! dimension along the Ok vertical axis
   iglb=430       ! dimension along the Oi axis (global grid)
   jglb=530       ! dimension along the Oj axis (global grid)
grid_j0=530       ! j0=j(longi0,latit0)
dxb=17000.        ! Along Oi axis cellboxe size (m)
dyb=17000.        ! Along Oj axis cellboxe size (m)
  phi0=69.0       ! Reference latitude for the Mercator projection
latit0=69.0       ! latitude of grid point (I0,J0) relative to north pole
longi0=0.       
grid_i0=110        
rayonterre=6370949.    
northpole_lon=-0.836
northpole_lat=37.822
southpole_lon=-0.805 
southpole_lat=37.635

! fplan1_grid=.true.       ! constant coriolis (lat=phi0)
  fplan1_grid=.false.        

  fplan2_grid=0
! fplan2_grid=1       ! dx_t(:,:)=dxb dy_t(:,:)=dyb
! fplan2_grid=2       ! dx_t(:,:) & dy_t(:,:) are constant along Oi axis
! fplan2_grid=3       ! dx_t(:,:) & dy_t(:,:) are constant along Oj axis


!--- Longitudes & Latitudes provided by a netcdf file
  lonlatfile='nofile'     ! case "no file"
! OFFLINE MODE: use ../../../BGC_MAR_MENOR/OFFLINE/grid.nc file
! lonlatfile='../../../BGC_MAR_MENOR/OFFLINE/grid.nc' ! OFFLINE MODE

! lonlatfile='../../../BGC_MAR_MENOR/BATHYMASK/lonlat_4col.txt'   ! format ascii i,j,lon,lat
!lonlatfile='../../../MAR_MENOR/OFFLINE/grid.nc' ! OFFLINE MODE
 lonlatfile='../../../BGC_MAR_MENOR/OFFLINE/full_grid.nc' ! OFFLINE MODE

!*** MPI SECTION *******
nbdom_imax=6                       ! number of subdomains along the Oi axis
nbdom_jmax=9                       ! number of subdomains along the Oj axis
iperiodicboundary=.true.          ! periodic boundaries in the Oi direction (if .true.)
jperiodicboundary=.true.           ! periodic boundaries in the Oj direction (if .true.)
discard_lonlat_periodicity=0        ! if =1 periodicity (if any) do not apply on lon,lat

! Distribution of subdomains in the mpi space:
! .. The "default" option automatically gives a regular mpi map including 100% masked subdomains. Other
! distributions requires an input file whose name (other than 'default') is given in the next line:
 mpi_map_file_name='default'        
!  mpi_map_file_name='../../../BGC_MAR_MENOR/BATHYMASK/description_domaine36.next'

! .. Lost subdomains (if any) must be recovered when netcdf files will be created. 
! Lost areas coordinates and ranks in charge of plugging are listed in a file whose name
! (other than 'none') is given in the next line:
 mpi_hole_plugging='none'
!  mpi_hole_plugging='../../../BGC_MAR_MENOR/BATHYMASK/description_trous36.txt'
! Details: https://docs.google.com/document/d/1CeW2GhCUSjmx_f7oIeFICK6gX-7WrV3hTIMCVgQ4-24/edit
!**********

!...............
! Web of canals
! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit
  webcanals_list='none'
! webcanals_list='../../../BGC_MAR_MENOR/LIST/webcanals_list'
!...............

!...............
! Name of the temporary output files:
tmpdirname='tmp/'
!...............

!...............................................................................
! Offline procedure (ioffline=2 in notebook_offline):
! Case: grid provided by NEMO
!initialgridfile_txt='../../../BGC_MAR_MENOR/OFFLINE/mesh_mask.ps.nouv.nc'
!vert_axis_conv_direc='dw' ! vertical axis: 'up'=upward (roms,s) 'dw'=downward (nemo) 
!vert_axis_conv_start='w'  ! vertical axis: 'w'=start with a w  point 't'=start with a t point
!vert_axis_conv_end='t'    ! vertical axis: 'w'=end   with a w  point 't'=end   with a t point 
!hori_axis_conv_start='t'  ! horizontal axis: 'v'=start with a velocity point 't'=start with a t point 
!hori_axis_conv_end='v'    ! horizontal axis: 'v'=end   with a velocity point 't'=end   with a t point 

! Case: grid provided by S model
initialgridfile_txt='none' ! S case='none' because the grid file will be the one found in the OFFLINE files directory
vert_axis_conv_direc='up'  ! vertical axis: 'up'=upward (roms,s) 'dw'=downward (nemo) 
vert_axis_conv_start='w'   ! vertical axis: 'w'=start with a w  point 't'=start with a t point
vert_axis_conv_end='w'     ! vertical axis: 'w'=end   with a w  point 't'=end   with a t point 
hori_axis_conv_start='t'   ! horizontal axis: 'v'=start with a velocity point 't'=start with a t point 
hori_axis_conv_end='t'     ! horizontal axis: 'v'=end   with a velocity point 't'=end   with a t point 
!...............................................................................


/

! Ne pas calculer les domaines 100% masquEs: 
https://docs.google.com/document/d/1CeW2GhCUSjmx_f7oIeFICK6gX-7WrV3hTIMCVgQ4-24/edit#

! parametrer notebook_grid et notebook_ofline pour simulations bio offline:
https://docs.google.com/document/d/17xbkkh_KwMwHof8Q0iDcKnOcS2ZIt_ztgkW0EzswARo/edit?ts=565c83db

! Simulations 1DV:
https://docs.google.com/document/d/1K-pQHjy77t0viUND4B9HwXCRYufeM-QcStr3wqlZQI4/edit

! Simulation 2DV le long des axes Oj,Ok periodique dans la direcionOi:
https://docs.google.com/document/d/1RXyYF3TIOxWHwXQPr8KHbxBy7p_qW3q59jNUJjY3Nng/edit

!Note: a 360° grid is obtained with
!DXB=DYB=rayonterre*2.*pi/(iglb-2)*cos(latit0(radians))
!and iperiodicboundary=.true.

      program calculette_bon_rayon_cercle_complet
      implicit none
      double precision rayonterre,pi
      pi=acos(-1.)
      rayonterre=6370949.000000  
      iglb=848
      latit0=89.42*pi/180
      dxb=rayonterre*2.*pi/(iglb-2)*cos(latit0)
      write(6,*) dxb
      end


