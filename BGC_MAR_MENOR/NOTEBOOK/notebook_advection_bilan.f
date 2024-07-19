&notebook_advection

! T & S advection scheme:
      iadvec_ts=2      ! 1=Quickest (horizontal & vertical) 
                       ! 2=Quickest (horizontal) and Centred  (vertical)
                       ! 3=Quickest (horizontal) and Lax Wen. (vertical)

flag_ts_effectivedensity=0 ! If =1 QUICK diffusivity is applied on the T,S effective density

! Define an "upwind" T,S advection" area according to bathymetry (upw_hrange1,upw_hrange2) in meters
      upw_hrange1=2.    ! 100%   upwind scheme si h<upw_hrange1
      upw_hrange2=10.   ! 100% standard scheme si h>upw_hrange2

! Amplification factor of the bilaplacian viscosity
! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit
      cst_adv_vel=1.     ! Amplification factor of the bi-harmonic viscosity in INTERNAL mode
      biharm_2dfactor=1. ! Amplification factor of the bi-harmonic viscosity in EXTERNAL mode
!     biharm_c1=0.       ! Ratio of bi-harmonic viscosity scheme (scheme1=constant, scheme2=velocity gradient) : valeur ds OFFLINE_SOL_ET_TIDE3
      biharm_c1=0.01       ! Ratio of bi-harmonic viscosity scheme (scheme1=constant, scheme2=velocity gradient)
                         ! biharm_c1=1   : scheme 1: bi-harmonic viscosity is constant
                         ! biharm_c1=0   : scheme 2: bi-harmonic viscosity is "velocity-gradient" dependent
                         ! 0<biharm_c1<1 : linear mixing of scheme1 and scheme2

! External mode: advection scheme of depth-averaged velocities in external mode
! https://docs.google.com/document/d/1U7DDLl9p1BZzG5NVgp4VWQ1Y6Y0if7Pt7bOjpr5GkfE/edit
      flag_adve2d=0      ! 0=c4+biharm , 1=upwind

! Coupling external and internal modes with momentum advection
! https://docs.google.com/document/d/1U7DDLl9p1BZzG5NVgp4VWQ1Y6Y0if7Pt7bOjpr5GkfE/edit
      flag_timesplitting_adve_uv=0 ! 0=coupling with u'u' (only)
                                   ! 1=coupling with depth-averaged total advection (+ delta of barotropic advection)

! Advection scheme for biogeochemical and passive tracers
      iadvec_bio=3     ! 3 = QUICKEST

      substep_advbio=1 ! Maximum number of sub time steps of vertical explicit advection
                       ! for passive tracer 

      flag_rmnegval=1  ! O = do not remove negative values 
                       ! 1 = use the negative values remover scheme 1
                       ! 2 = use the negative values remover scheme 2

! Cancel advection of biogeochemical and passive tracers in lateral border
! layers. Define the size (grid point number) of the 4 lateral BL:
      ibl1_advbio=0 ! side "i=1"
      ibl2_advbio=0 ! side "i=iglb"
      jbl1_advbio=0 ! side "j=1"
      jbl2_advbio=0 ! side "j=jglb"

! Central points of upwind zones
! https://docs.google.com/document/d/1djA4xtqN5U-kf0c0ezLFCydew2Ln4VrcIBjUh8G-dAs/edit
 texte250='none'
!texte250='../../../REGION/BATHYMASK/central_points_of_upwind_zones.txt'

/
