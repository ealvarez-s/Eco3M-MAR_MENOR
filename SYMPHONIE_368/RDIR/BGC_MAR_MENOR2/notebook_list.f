&notebook_list

!main points to check before starting the simulation:
!https://docs.google.com/document/d/1TRF8uYjsVen8EiY0rERHw3YxiofZmpZ6cH_sEhWiMwA/edit?usp=sharing

! DIRECTORY
directory='../../../BGC_MAR_MENOR/NOTEBOOK'      ! Directory of the notebooks

! TIME
nomfichier(1)='notebook_time_LR_LF.f'         ! Departure/End time of the runs, time steps,...

! GRID
nomfichier(2) ='notebook_grid24_LR.f'             ! 1 node
      
nomfichier(3) ='notebook_bathy_LR.f'              ! Land/Sea mask, bathymetry, wetdrying,...   
nomfichier(13)='notebook_vertcoord_LR.f'          ! Vertical coordinate, sigma stretching,...

! FORCING
!nomfichier(4) ='notebook_rivers_LR'              ! 10 RIVERS (legos)
nomfichier(4) ='notebook_rivers_SWAT'            ! 23 RIVERS (swat)      
!nomfichier(4) ='notebook_rivers_horario_LR'      ! 10 RIVERS in dana   
nomfichier(7) ='notebook_airseaflux_ecmwf_s26.f' ! METEO
nomfichier(8) ='notebook_obcforcing_nemo.f'       ! OGCM
!nomfichier(8) ='notebook_obcforcing_sympa_s26.f'  ! OGCM
nomfichier(22)='notebook_wave.f' ! WAVES
!nomfichier(11)='notebook_tide'                   ! TIDES
!nomfichier(11)='notebook_tide_fes2014'     
!nomfichier(11)='notebook_tide_s26_fes2012        
nomfichier(11)='notebook_tide_fes2012_v211'

!I/O
!nomfichier(20)='notebook_offline_LR_T35.f'              ! Offline files T35 (no rivers)
!nomfichier(20)='notebook_offline_LR.f'                  ! Offline files T37 (swat rivers)     
nomfichier(20)='notebook_offline_LR_ssh-mean.f'          ! Offline files T37 with ssh mean (swat rivers)
!nomfichier(20)='notebook_offline_LR_ssh-mean_runoffCNT.f'  ! Offline files T37 with runoff cnt and ssh mean (swat rivers)
!nomfichier(20)='notebook_offline_LR_danaOct25.f'          ! Offline files danaOct2025 + T35 (10 rivers biel)
nomfichier(21)='notebook_graph_LR'                       ! Outputs files for graph

! PHYSIC
nomfichier(5)='notebook_advection.f'              ! Advection schemes
nomfichier(9)='notebook_visco.f'                  ! Turbulence schemes
nomfichier(15)='notebook_optical.f'               ! Light attenuation
nomfichier(17)='notebook_eqstate.f'               ! Equations of state
nomfichier(14)='notebook_sponge_LR.f'        ! OBC schemes, nudging layer
nomfichier(34)='notebook_nh.f'                    ! m0v0m

! TRACERS
nomfichier(10)='notebook_tracer.f'                       ! Eulerian (passive)
!nomfichier(12)='notebook_bio_10highInputx10_CAMSco2'     ! co2 file (ppm) from CAMS monthly
!nomfichier(12)='notebook_bio_lowSed_01onlyInput_CAMSco2'  ! rios 1-3, [baja]
!nomfichier(12)='notebook_bio_lowSed_23realInput_chs'    ! real-date concentrations CHS
!nomfichier(12)='notebook_bio_lowSed_23realInput_swat'   ! real-date concentrations SWAT
!nomfichier(12)='notebook_bio_lowSed_23realInput_chla'   ! 
!nomfichier(12)='notebook_bio_lowSed_23realInput_plus'   ! 
!nomfichier(12)='notebook_bio_lowSed_23realInput_avSi'   ! 
!nomfichier(12)='notebook_bio_lowSed_23cntInput_250'   ! constant concentrations
!nomfichier(12)='notebook_bio_lowSed_23cntInput_500'   ! constant concentrations
!nomfichier(12)='notebook_bio_lowSed_23cntInput_5R'    ! constant concentrations
nomfichier(12)='notebook_bio_lowSed_23zeroInput'         ! all input to zero
!nomfichier(12)='notebook_bio_lowSed_23zeroInput_CAMSco2' ! co2 file (ppm) from CAMS monthly
!!nomfichier(12)='notebook_bio_lowSed_23highInput_noSi'  ! constant input but Si 
!!nomfichier(12)='notebook_bio_lowSed_23noRamblas'       ! zero input
nomfichier(16)='notebook_drifter' ! Lagrangian

! BIO
nomfichier(23)='notebook_light'
nomfichier(24)='notebook_zooplankton_test14'
!nomfichier(24)='notebook_zooplankton_test17'
!nomfichier(24)='notebook_zooplankton_testNoGraz'
!nomfichier(25)='notebook_phytoplankton_modified'
nomfichier(25)='notebook_phytoplankton_fay'
nomfichier(26)='notebook_bacteria_fay'
nomfichier(27)='notebook_remineralisation_test11'
nomfichier(28)='notebook_initpelagic'
nomfichier(29)='notebook_biobcforcing'
!nomfichier(30)='notebook_benthic'       ! Benthic=1
nomfichier(30)='notebook_benthic2'       ! Benthic=2
nomfichier(31)='notebook_oxygen2'

! OTHERS
nomfichier(18)='notebook_dateoutput'
nomfichier(19)='notebook_atlas_SWAT'
nomfichier(33)='notebook_sedim.f'

! OASIS COUPLER
!nomfichier(32)='notebook_oasis'
nomfichier(32)='notebook_oasis_sym_sym'

/
