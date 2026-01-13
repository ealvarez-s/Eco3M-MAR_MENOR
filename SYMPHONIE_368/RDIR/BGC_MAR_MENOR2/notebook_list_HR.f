&notebook_list

!main points to check before starting the simulation:
!https://docs.google.com/document/d/1TRF8uYjsVen8EiY0rERHw3YxiofZmpZ6cH_sEhWiMwA/edit?usp=sharing

! DIRECTORY
directory='../../../BGC_MAR_MENOR/NOTEBOOK'             ! Directory of the notebooks

! TIME
nomfichier(1)='notebook_time.f'                   ! Departure/End time of the runs, time steps,...

! GRID
!nomfichier(2) ='notebook_grid_ieo23.f'            ! 24 total, 23 active (4x6)     
nomfichier(2) ='notebook_grid_ieo48.f'            ! 54 total, 48 active (6x9)
!nomfichier(2) ='notebook_grid_ieo62.f'            ! 72 total, 62 active (8x9)
!nomfichier(2) ='notebook_grid_ieo70.f'            ! 84 total, 70 active (7x12)
!nomfichier(2) ='notebook_grid_ieo72.f'            ! 84 total, 72 active (12x7)     
!nomfichier(2) ='notebook_grid_ieo92.f'            ! 108 total, 92 active (9x12)
!nomfichier(2) ='notebook_grid_ieo94.f'            ! 112 total, 94 active (8x14)
      
nomfichier(3) ='notebook_bathy.f'                 ! Land/Sea mask, bathymetry, wetdrying,...   
nomfichier(13)='notebook_vertcoord.f'             ! Vertical coordinate, sigma stretching,...

! FORCING
nomfichier(4) ='notebook_rivers'                  ! RIVERS
nomfichier(7) ='notebook_airseaflux_ecmwf_s26.f'  ! METEO
nomfichier(8) ='notebook_obcforcing_nemo.f'       ! OGCM
nomfichier(22)='notebook_wave.f'                  ! WAVES
!nomfichier(11)='notebook_tide'                   ! TIDES
!nomfichier(11)='notebook_tide_fes2014'     
!nomfichier(11)='notebook_tide_s26_fes2012        
 nomfichier(11)='notebook_tide_fes2012_v211'

! I/O
nomfichier(20)='notebook_offline.f'              ! Offline files
nomfichier(21)='notebook_graph'                  ! Outputs files for graph

! PHYSIC
nomfichier(5)='notebook_advection.f'              ! Advection schemes
nomfichier(9)='notebook_visco.f'                  ! Turbulence schemes
nomfichier(15)='notebook_optical.f'               ! Light attenuation
nomfichier(17)='notebook_eqstate.f'               ! Equations of state
nomfichier(14)='notebook_spongelayer.f'           ! OBC schemes, nudging layer
nomfichier(34)='notebook_nh.f'                    ! m0v0m

! TRACERS
nomfichier(10)='notebook_tracer.f'                ! Eulerian (passive)
nomfichier(12)='notebook_bio_lowSed_10highInput'  ! Eulerian (bio)
nomfichier(16)='notebook_drifter'                 ! Lagrangian

! BIO
nomfichier(23)='notebook_light'
nomfichier(24)='notebook_zooplankton_test14'
!nomfichier(24)='notebook_zooplankton_test17'
!nomfichier(24)='notebook_zooplankton_testNoGraz'
nomfichier(25)='notebook_phytoplankton_fay'
!nomfichier(25)='notebook_phytoplankton_fay3'
nomfichier(26)='notebook_bacteria_fay'
nomfichier(27)='notebook_remineralisation_test11'
nomfichier(28)='notebook_initpelagic'
nomfichier(29)='notebook_biobcforcing'
nomfichier(30)='notebook_benthic'
nomfichier(31)='notebook_oxygen2'

! OTHERS
nomfichier(18)='notebook_dateoutput'
nomfichier(19)='notebook_atlas_LEGOS'
nomfichier(33)='notebook_sedim.f'

! OASIS COUPLER
!nomfichier(32)='notebook_oasis'
nomfichier(32)='notebook_oasis_sym_sym'

/
