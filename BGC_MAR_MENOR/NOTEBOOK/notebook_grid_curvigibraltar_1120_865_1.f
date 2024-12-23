&notebook_grid

! this namelist is read by the s26 routine main.F90 and by the preprocessing routine dom.F90

!--- Usual parameters of the "old" notebook_grid version:
iglb=1120 ! 2240          ! dimension along the Oi axis (global grid)
jglb=865  ! 1730          ! dimension along the Oj axis (global grid)
kmax=50   ! 50                    ! dimension along the Ok vertical axis
dxb=900.                  ! Along Oi axis cellboxe size (m)
dyb=900.                  ! Along Oj axis cellboxe size (m)
phi0=77.                  ! Reference latitude for the Mercator projection
longi0=20.8               ! longitude of grid point (I0,J0) relative to north pole
latit0=77.                ! latitude of grid point (I0,J0) relative to north pole
grid_i0=700 ! 160               ! i0=i(longi0,latit0)

grid_j0=1730 ! 706               ! j0=j(longi0,latit0)

rayonterre=6370949.       ! Earth radius (m)
northpole_lon=-6.5        ! longitude (° decimal) of the grid north pole
northpole_lat=54.7        ! latitude  (° decimal) of the grid north pole
southpole_lon=9999.       ! longitude (° decimal) of the grid south pole (antipode=9999.)
southpole_lat=9999.       ! latitude  (° decimal) of the grid south pole (antipode=9999.)

! fplan_grid=.true.       ! constant coriolis (lat=phi0) and constant resolution dx=dxb dy=dyb
! fplan_grid=.false.      ! 
  fplan1_grid=.false.
  fplan2_grid=0


!--- Longitudes & Latitudes provided by a netcdf file
! lonlatfile='nofile'     ! case "no file"
! lonlatfile='../../../GLOBMED/BATHYMASK/lonlat.nc' 
! lonlatfile='../../../GLOBMED/OFFLINE/grid.nc' ! format netcdf
  lonlatfile='../../../GLOBMED/BATHYMASK/lonlat_GLOBMED_curvigibraltar_1120_865_4col.txt'

!*** MPI SECTION *******
!nbdom_imax=24 ! 24                       ! number of subdomains along the Oi axis
!nbdom_jmax=22 ! 20                        ! number of subdomains along the Oj axis
nbdom_imax=34 ! 24                       ! number of subdomains along the Oi axis
nbdom_jmax=30 ! 20                        ! number of subdomains along the Oj axis
! iperiodicboundary=.true.          ! periodic boundaries in the Oi direction (if .true.)
iperiodicboundary=.false.           ! periodic boundaries in the Oi direction (if .true.)
jperiodicboundary=.false.           ! periodic boundaries in the Oj direction (if .true.)

! Distribution of subdomains in the mpi space:
! .. The "default" option automatically gives a regular mpi map including 100% masked subdomains. Other
! distributions requires an input file whose name (other than 'default') is given in the next line:
! mpi_map_file_name='default'        
! mpi_map_file_name='../../../GLOBMED/BATHYMASK/description_domaine_1120_865.next'
 mpi_map_file_name='../../../GLOBMED/BATHYMASK/description_domaine_34_30.next'

! .. Lost subdomains (if any) must be recovered when netcdf files will be created. 
! Lost areas coordinates and ranks in charge of plugging are listed in a file whose name
! (other than 'none') is given in the next line:
! mpi_hole_plugging='none'
! mpi_hole_plugging='../../../GLOBMED/BATHYMASK/description_trous_1120_865.txt'
 mpi_hole_plugging='../../../GLOBMED/BATHYMASK/description_trous_34_30.txt'
! Details: https://docs.google.com/document/d/1CeW2GhCUSjmx_f7oIeFICK6gX-7WrV3hTIMCVgQ4-24/edit
!**********

!...............................................................................
! Offline procedure (ioffline=2 in notebook_offline):
! Case: grid provided by NEMO
!initialgridfile_txt='../../../GLOBMED/OFFLINE/mesh_mask.ps.nouv.nc'
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


