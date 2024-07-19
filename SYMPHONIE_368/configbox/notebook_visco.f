&notebook_visco

!....................................
! Turbulence closure:
      iturbulence=1     ! Turbulence closure scheme: 
                        ! 0: Gaspar1990 
                        ! 1: K-Eps (Michaud et al, OS 2012)
                        ! 2: Kh=constant_Kh, Km=constant_Km
                        ! 3: Kh=Km=((karman*L)**2)*sqrt( (du/dz)**2+(dv/dz)**2 ) Prandlt
                        ! 4: "tke NEMO"
                        ! 5: "tke NEMO" neutral
      constant_km=1.e-6 ! Relevant if iturbulence=2
      constant_kh=1.e-7 ! Relevant if iturbulence=2

      emin=1.d-7        ! Minimum value of the Turbulent Kinetic Energy (TKE)     
      epsmin=1.d-12     ! Minimum value for EPSIPLON       
      tke_surf=1        ! TKE surf cond: 0=PROduction DIssipation Equilibrium, 1=Craig & Banner

! Stability functions (relevant if k-eps selected):
!     flag_tke_stab_func='kc94' ! kantha Clayson 1994
      flag_tke_stab_func='ca01' ! Canuto A 2001
!     flag_tke_stab_func='gp88' ! Galperin 1988

! relevant if Gaspar selected:
      ctke1=0.1 !
      ctke2=0.7 ! 
! Gaspar equivalent to Craig et Banner 1994 scheme (SM,B)=(0.39,16.6) (provided Lespilon=Lkm)
! Be careful, Gaspar omitted the karman constant in his formulation which partly explains the deviation of his coefficients 
! from those generally given in the literature. To find the equivalence, it is necessary to consider that the mixing length 
! of the other schemes (Lo) is the length of Gaspar (Lg) multiplied by Karman (km=0.41): Lo=km*Lg
! Examples:
! Gaspar equivalent to Craig et Banner 1994 scheme (SM,B)=(0.39,16.6):
! ctke1=0.225  !=0.55*karman  Km=Lo.q.SM=Lo.sqrt(2.tken).SM=Lg.km.sqrt(tken).sqrt(2).SM ! ==>ctke1=km*sqrt(2).SM=0.55*km=0.225
! ctke2=0.415  !=0.17/karman  epsilon=2.sqrt(2.tken).tken/(B.Lo)=(tken**3/2)*2.sqrt(2)/(B*km*Lg) ==> ctke2=2.sqrt(2)/(B*km)=0.415

!....................................
! At molecular scale: 
! Apel J. R., 1987, Principle of Ocean Physics, Table in Appendix Three:
      kmol_m=0.  !  1.049e-6   ! molecular viscosity
      kmol_h=0.  !  1.49e-7    ! molecular diffusivity for heat
      kmol_s=0.  !  1.5e-9     ! molecular diffusivity for salt


!....................................
! SURFACE LAYER:
      z0s=0.01                ! Default surface roughness (meters)
      momentum_input_depth=1. ! Default momentum input depth (meters)

!....................................
! BOTTOM LAYER:
!  Bottom roughness 
      z0b=0.001      ! Sea bottom roughness (meters)
! z0b from a netcdf file if available:
!texte250='../../../REGION/BATHYMASK/symphonie-symtools-comodo.nc' ! lon-lat-dependent from file
      z0b_land=0.01    ! Bottom roughness on land. Useful for flooding of land, ...
      zlevel_land=0.   ! (related to z0b_land) z levels (in meters) of sea land transition.
      z0b_rivers=0.01  ! Bed rivers bottom roughness 
! local values if available:
!https://docs.google.com/document/d/1cgBEAL0qvT7ernFu7FoxlNq8zwX7vhhwB5NNFeBgsb4/edit
!texte90='../../../REGION/BATHYMASK/z0b.ijvaldist'

! Macro roughness scale (effect of rocks on friction,....) in addition to small-scale friction
! https://docs.google.com/document/d/1QkYuJuw6InVor7g8NKgxRPjSJPrS6TLwC4kxOggG3c4/edit?usp=sharing
  flag_z0_macro=0      ! (=1 if used, =0 if unused))

! Avoid the singularity of the logarithmic layer by limiting the ratio (z(k=1)+h)/z0
  dz_over_z0_min=10.   ! dz_over_z0_min is the limit for (z(k=1)+h)/z0


!  Bottom Drag Coefficient:
      flag_linearfric=0      ! Linear bottom friction if flag_linearfric=1
      coef_linearfric=0.0001 ! (relevant if flag_linearfric=1) coef of the linear friction
      cdseuil=1.d-3          ! Minimum value of the bottom drag coefficient    
      cdb_2dh=2.d-3          ! Drag coefficient for 2DH (kmax=1) simulations

! coastal no slip condition: [coastal_viscosity*DX] is the horizontal viscosity in m**2/s (<=0 if free splip condition)
! Example: coastal_viscosity=0.01 and DX=1000 ==> horizontal viscosity = 10m**2/s
! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit?usp=sharing
coastal_viscosity=0.      ! du/dt=-(1/dy**2)*coastal_viscosity*u, dv/dt=-(1/DX)*coastal_viscosity*v (1<nu_no_slip<100)

!....................................
! Convection:
! https://docs.google.com/presentation/d/19f3SSwAQ4flJZG0lJAP-axdwOeSqtoRbucJojyTE_Ac/edit#slide=id.g57193dca64_0_11
      convect_yn=0    ! 1 or 2 or 3 = Quick convection methods, 0 otherwise                    
!....................................
! Time stepping:
      assel0=0.1 ! 0.05     ! Time filter coefficient for the Laplacian scheme (Marsaleix et al OM 2012)

/
