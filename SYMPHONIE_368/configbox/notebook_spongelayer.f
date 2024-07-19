&notebook_spongelayer

! Define the parameters of the lateral nudging (or sponge) layer

!............................
! For tracers only:
! Sponge layer depending on distance ('dist') or on horizontal resolution ('dxdy') or full grid ('full')
 flag_sponge_txt='dist'
!flag_sponge_txt='full'
!flag_sponge_txt='dxdy'
! details of the 'dxdy' option in:
! https://docs.google.com/document/d/10A6XnRxEiuXqTIloBfFt_G-bFn-UHUw95E49TuttA0w/edit
!............................

! sponge_l is the size of the nudging/sponge layer for velocities (and tracers if flag_sponge_txt=='dist')
sponge_l=30     ! Width (in grid nodes) of the nudging/sponge layer              

! sponge_dx_critic is the resolution below which the restoring force is applied (relevant if flag_sponge_txt=='dxdy')
 sponge_dx_critic=2500.  !meters (relevant if flag_sponge_txt=='dxdy')

! sponge_dx_width is the decay scale of the hyperbolic tangent function (relevant if flag_sponge_txt=='dxdy')
 sponge_dx_width=1000. !meters (relevant if flag_sponge_txt=='dxdy')

! Barotropic velocity:
relax_ext=0.1   ! Nudging time scale (days) for barotropic velocity. 
                ! UNUSED IF relax_ext=0. 

! Baroclinic velocity:
relax_int=1.    ! Nudging time scale (days) for baroclinic velocity. 
                ! UNUSED IF relax_int=0.

! Temperature and Salinity:
! https://docs.google.com/document/d/1bt_mnhZCXMuNnw3JqymtFjInqG6nNHQT0vuXX0HtZps/edit
! Scheme:
relaxtype_ts=1  ! 2= low frequency filter
                ! 1= restoring force is current dependent
                ! 5= restoring force is based on relax_ts 

relax_ts=60.    ! Case relaxtype_ts=2: relax_ts is a nudging time scale (days)
                ! Case relaxtype_ts=1: the nudging time scale depends on velobc_u  velobc_v
                ! Case relaxtype_ts=6: the nudging time scale depends on vel_u     vel_v
                ! UNUSED IF relax_ts<0.
relax_lwf=10.   ! Case relaxtype_ts=2: relax_lwf is a time scale (days) separating long scales from short ones.

! Cancel sponge layer along specific boundaries if flagspoXX=0
flagspo_i1=1 ! i=1    boundary
flagspo_i2=1 ! i=iglb boundary 
flagspo_j1=1 ! j=1    boundary
flagspo_j2=1 ! j=jglb boundary


! T,S,u,v advection in spongelayer
! https://docs.google.com/document/d/1djA4xtqN5U-kf0c0ezLFCydew2Ln4VrcIBjUh8G-dAs/edit
flag_upwind_obc=0 ! upwind advection for T,S,u,v at the open borders if flag_upwind_obc=1

/
