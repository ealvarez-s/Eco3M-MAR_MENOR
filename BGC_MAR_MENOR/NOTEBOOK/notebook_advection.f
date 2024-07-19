&notebook_advection

! T & S horizontal advection scheme:
      iadvec_ts_hor=1      ! 1=Quickest   
                           ! 2=Quickest-2 
! T & S vertical advection scheme:
      iadvec_ts_ver=0      ! 0=C2
                           ! 1=Quickest   
                           ! 2=Quickest-2

! QUICKEST diffusivity:
! https://docs.google.com/document/d/1-4NdFRKCug7u4-reURBVsHOeHMtf5p_8RlCYDLwAZjA/edit
!quick_coef=0.125  ! 1/8 or 1/6, in factor of the CURV term of the QUICKEST scheme
 quick_coef=0.1667 ! 1/8 or 1/6, in factor of the CURV term of the QUICKEST scheme

! relevant if iadvec_ts_hor=1 or 2:
flag_ofactor=0         ! 1=apply orthogonality factor

! Negative diffusion
flag_negdif_hor=1            ! 1 or 0 = on or off
ratio_negdif_hor=1.          ! ratio antidiffusion/diffusion 
flag_negdif_ver=1            ! 1 or 0 = on or off
ratio_negdif_ver=1.          ! ratio antidiffusion/diffusion 
quick_filter_points=0        ! =0 =3 or =5 =number of points in the filter for vertical advection
                             ! Note: use quick_filter_points=0 with iadvec_ts_ver=0

! Define an "upwind" T,S advection" area according to bathymetry (upw_hrange1,upw_hrange2) in meters
      upw_hrange1=0.    ! 100%   upwind scheme si h<upw_hrange1
      upw_hrange2=1.    ! 100% standard scheme si h>upw_hrange2


! MOMENTUM
! Amplification factor of the bilaplacian viscosity
! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit
      cst_adv_vel=1.     ! Amplification factor of the bi-harmonic viscosity in INTERNAL mode
      biharm_2dfactor=1. ! Amplification factor of the bi-harmonic viscosity in EXTERNAL mode

! External mode: advection scheme of depth-averaged velocities in external mode
! https://docs.google.com/document/d/1U7DDLl9p1BZzG5NVgp4VWQ1Y6Y0if7Pt7bOjpr5GkfE/edit
      flag_adve2d=2 ! 0=c4+biharm (version1) , 1=upwind , 2=c4+biharm(version2)
! Internal mode: advection scheme for baroclinic velocities
      flag_adve3d=2 ! 0=c4+biharm (version1) , 1=upwind , 2=c4+biharm(version2)

! Advection scheme for biogeochemical and passive tracers
      iadvec_bio=3     ! 3 = QUICKEST

      flag_rmnegval=0     ! 0 = do not remove negative values 
                          ! 1 = use the negative values remover scheme 1 (vertical direction only)
                          ! 2 = use the negative values remover scheme 2 (tridimensional)

! Cancel advection of biogeochemical and passive tracers in lateral border
! layers. Define the size (grid point number) of the 4 lateral BL:
      ibl1_advbio=0 ! side "i=1"
      ibl2_advbio=0 ! side "i=iglb"
      jbl1_advbio=0 ! side "j=1"
      jbl2_advbio=0 ! side "j=jglb"

! Central points of upwind zones
! https://docs.google.com/document/d/1djA4xtqN5U-kf0c0ezLFCydew2Ln4VrcIBjUh8G-dAs/edit
 upwindzone_file='none'
!upwindzonefile='../../../MAR_MENOR/BATHYMASK/central_points_of_upwind_zones.txt'

! implicit horizontal advection settings:
! https://docs.google.com/document/d/1Z9G4Sl1fcbKqju7RdzyvMsD11_LcPNi3RdjQvF3t-kI/edit?usp=sharing
looplimit_hor=1

/
