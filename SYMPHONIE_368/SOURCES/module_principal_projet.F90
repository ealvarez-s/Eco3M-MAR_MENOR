! POUR Faire un restart independant du nombre de processeur
! l'ideal serait de faire une simulation (3D bien sûr) où
! j'ecrirais dans un fichier les tableaux qui ont ete rejetes par
! module_netcdfrestart et écrit par write(9)



















!______________________________________________________________________
! S model
! release S26.1 - last update: 01-07-15
!______________________________________________________________________
! Version Date      Description des modifications
!         26/12/02: VHZ_X VHZ_Y OMEGA_Z augmentés pour filiere forward
!         27/12/02: bienvenue à KLOW_X, KLOW_Y, KLOW_Z
!         27/01/03: KOUNTGRAPH devient real (nb: on a oublié de le faire
!c                                            dans symphonie2002.h, dommage...)
!c                   bienvenue à IOFFLINE, DFVOFL_Z, OFFLINEDT, OFFLINEFILE
!c                   bienvenue à KMAX_DOF, IDATE_OUTPUT, KDATE_OUTPUT
!c         05/02/03: binevenue à X14
!c         19/03/03: EVA_Z supprimé, PRECIPI_Z ajouté
!c         24/03/03: ajout de IHYBSIG
!c         09/07/03: ajout de TKE_SURF
!c         10/07/03: OBCFILE taille aggrandie de (8,3) à (9,3)
!c         25/07/03: bienvenue à Y1 & Y2
!c         28/07/03: bienvenue à GRH_OUT_MI, GRH_OUT_VAR, GRH_NB & GRH_TITRE
!c                   L1_SCA, L1_VEC, L2_SCA, L2_VEC, L3_SCA, L3_VEC
!c         29/07/03: bienvenue à VMI_X & VMI_Y tableaux pour faire des
!c                   moyennes temporelles des courants et SUM_MI
!c         02/08/03: ajout de NHYBSIG
!c         07/08/03: DIMENSIONs de ANYVAR3D aggrandies. Avant NR, maintenant 0:NR
!c         21/08/03: supression TIDESTEP1
!c         28/08/03: bienvenue à NC1 et à toute la famille "NEST"
!c         21/10/03: bienvenue à RELAX_R
!c         02/11/03: plus de memoire pour NRJ3D
!c         30/11/03: plus de memoire pour XY_X et XY_Y
!c         17/04/04: on ecrit les DIMENSIONs à partir de la colonne 23
!c         24/04/04: La zone de bouclage d'update des tableaux de forcage sur
!c                   zone eponge est mise en memoire dans des tableaux
!c         05/05/04: bienvenue forcage HF et temperature des rivieres evolutive
!c         06/05/04: bienvenue NCMIN_AIRSEA NCMIN_RIVER NCMIN_AF
!c         25/05/04: la nouvelle filiere d'imbrication permet de reduire
!c                   les DIMENSION de obcfile
!c         28/05/04: suppression relaxmin _int & _ext
!c         01/11/04: bienvenue à la pression atmopherique moyenne PSS_MEAN
!c         01/11/04: bienvenue à BI_ONOFF pour activer ou non le barometre
!c                   inversea aux obc
!c         13/11/04: bienvenue à EQS_STATE1 EQS_STATE2
!c         02/12/04: GAMA et BETA disparaissent. ANYV3D(i,j,k,2ou1) les
!c                   remplace.
!c         15/12/04: mask tem et sal aggrandis pour schema d'advection
!c         26/04/05: bienvenue à GIVE_CHANEL9 qui permet de demander un
!c                   fichier restart sans arreter le model_.
!c         27/04/05: NEST_ONOFF_DEMO
!c         03/06/05: bienvenue à T0SURF S0SURF RHO_TMP
!c         19/06/05: CST_ADV_HOR CST_ADV_VER
!c         19/07/05: augmentation taille memoire tableaux spo_i1_x et cie....
!c                   (on passe de 4 à 5)
!c         25/07/05: bienvenue à EQS_STATE3
!c         29/07/05: bienvenue à const8
!c         03/08/05: bienvenue à ZTAREFOBC_I ZTAREFOBC_J OBCFREEORFIX
!c         27/09/05: DIMENSIONs de ANYVARLR aggrandies suite à
!c                   debug dans nest_inout.F
!c         28/12/05: bienvenue à X1_R8 X2_R8 X3_R8 X4_R8
!c                   à IALBEDO
!c         17/01/06: bienvenue à IWAVE (active model_houle de Cléa Denamiel)
!c         20/01/06: bienvenue à ZTAOBC_CUM et à LAGRANGE_ZTA
!c         26/01/06: supression ANYVAR3D et aggrandissement ANYV3D(i,j,k,0:4)
!c         16/02/06: bienvenue a RELAX_TS,LIGHT_ATT1,LIGHT_ATT2,LIGHT_RAT1
!c                  ,LIGHT_RAT2
!c         14/03/06: reduction taille des tableaux OBC
!c         28/03/06: bienvenue à RAP_OBC
!c         07/04/06: bienvenue a SP_OR_DB (simple ou double precision)
!c         11/04/06: bienvenue a ANYV1D
!c         13/04/06: bienvenue a LNAME4 DIRGRAPH et TEXTE250
!c         18/04/06: anyvarlr et dfvofl_z force real4
!c                   REAL regroupe pour passer vite de simple a double
!c                   precision
!c         04/05/06: les tableaux de forcage HF barotropes sont forces
!c                   a la simple precision
!c         08/06/06: bienvenue a LV
!c         16/06/06: bienvenue CONVECT_YN
!c         02/07/06: schema tke forward entraine elimination TKENHZ_Z
!c         05/07/06: suppression ADHY et cie
!c         15/09/06: bienvenue à DIF3D2D_X et DIF3D2D_Y et CST_ADV_VEL
!c         04/12/06: bienvenue à VSUMSAVE
!c         19/01/07: CN: Convertion common en module et allocation dynamique
!          08/04/07: Ajout des tableaux dx_x et cie...
!          30/04/07: debug taille tableau grdobc3d_i
!          09/05/07: bienvenue à ADVMIX_X & ADVMIX_Y
!          31/05/07: bienvenue à ZTAREFSAVE_I VBRREFSAVE_I ZTAREFSAVE_J
!                    VBRREFSAVE_J
!          08/06/07: Bienvenue à HZDY_X et HZDX_Y
!          18/06/07: Bienvenue à DT_OBC
!          02/10/07: bienvenue à ANYVAR3D
!          28/11/07: ajout de nouveau tableaux
!          15/01/08: bienvenue à KSL_Z
!          17/01/08: IWAVE renomme IWVE
!          18/01/08: Bienvenue à ZTATOTIDE_Z, la somme de toutes les harmoniques
!                    et aussi Z1MEATIDE Z2MEATIDE ZTAMEATIDE
!          21/01/08: Bienvenue à OBCTYPE_TS et RELAX_BPC
!          24/02/08: Bienvenue à OBCTYPE_P et modif taille VMEAN_
!          31/03/08: Bienvenue à HSTEPMIN HSTEPMAX NBVSTEPMIN
!          07/04/08: ajout RELAXTYPE_TS pour distinguer la relaxation des
!                    traceur en mode classique du mode "densité inchangée"
!          18/04/08: bienvenue à KOUNT0_RST
!          13/05/08: extension memoire THZ_Z, SHZ_Z, TEM_Z, SAL_Z, ANYV3D
!          16/10/08: extension memoire   mask_x &   mask_y
!          01/12/08: introduction de tableau pour calculer l'interpolation
!                    "online" des champs meteo
!          06/12/08: Retour en arriere par rapport au point precedent
!          18/02/09: bienvenue à Z1 et Z2
!          23/02/09: bienvenue à DX_R et DY_R
!          24/03/09: parallelisation
!          27/04/09: parallelisation: RIVER_DOM attribue un numero de sous domaine
!                    à chaque fleuve
!          28/04/09: parallelisation: MSKI_Y MKSJ_X eviteront de sommer 2 fois
!                    les zones de chevauchement
!          19-05-09  bienvenue à ZTA_AVR_NEST_OUT, GRID_AREA
!          30-05-09  ajout SUM1GLB
!          31-05-09  MSKI_Y et MSKJ_X deviennent MASK_I_Z & MASK_J_Z
!          01-06-09  ajout SUM0GLB, SUM2GLB, SUM3GLB, SUM4GLB, SUM5GLB, SUM6GLB
!          04-06-09  ajout ZTATIDE_CUM, DLON, DLAT
!          05-06-09  modif dim ZTAMEATIDE
!          14-06-09  ajout RESTART_FILE_Y_OR_N
! 2009-2   07-09-09  ajout ogcm_rec_max
!          08-09-09  ajout ind
! 2009-3   03-10-09  ajout des increments pour ecriture generique
!          05-10-09  ajout de grid_areaglb
!          08-10-09  albedo devient une variable 2D
!          10-10-09  tem_c et sal_c dimensionnés de 0 à 2 pour le temps
!          08-11-09  ajout x_xhl_dim etc...
!          15-11-09  suppression difhor_x,difhor_y,difhor_z,difhor_r
!          06-12-09  ajout ipu,ipv etc..
! 2010.2   17-12-09  ajout deg2rad et rad2deg
!          19-12-09  dayindex_size est la taille effective du tableaux
!                    (c.a.d. nombre de valeurs non nulles)
!          21-12-09  modif dimensions parametres nodaux
!          22-12-09  ajout tide_year_min, tide_year_max ,kounttide_prev_rdv
!                   ,kounttide_next_rdv
!          28-12-09  suppression tidevel
! 2010.3   06-01-10  ajout wetmask_t, wetmask_x, wetmask_y
!          08-01-10  ajout pole_lon, pole_lat, gridrotcos_z, gridrotsin_z
!          12-01-10  ajout wetdry_cst2, wetdry_cst3
!          14-01-10  cfl_sshmax cfl_umax cfl_reduce
! 2010.5   29-01-10  drifter_t
! 2010.6   02-01-10  renomme lat_t lon_t
!          03-02-10  nc_or_bin_airsea, dimensio airseainfo
! 2010.7   15-02-10  modif dimension tableaux de vitesse
!          17-02-10  difvel renommé velbot3d2d
!          22-02-10  reparation d'un nom (h_fiver -> hriver) et ajout friver
!          24-02-10  ajout rap_wave t_wave_t hs_wave_t kx_wave_t ky_wave_t
!          25-02-10  dir_wave_t wavedt velstokes0m_u velstokes0m_v
!          02-02-10  dimensions xy_t
!          03-03-10  ajout ustokesvortex_t ustokesvortex_f vstokesvortex_t
!                    vstokesvortex_f stokesvortex3d2d_u stokesvortex3d2d_v
!          04-03-10  ajout sshstokes_t
! 2010.8   10-03-10  - ajout before2 before3
!                    - modif dimension vel_u vel_v tem_t sal_t
!          11-03-10  modif dimension velbar_u, velbar_v
!          12-03-10  tfc1 tfc2
!          18-03-10  - supprime obcdt
!                    - augmenter dimension wetmask
!          19-03-10  ajout heatrelax_w irelaxsst
!          22-03-10  incides time stepping desormais definis dans parameter
!                    ajout hssh_w
!          21-05-10  ajout stress des vagues
! 2010.9   31-05-10  Indexation generique maintenant declaree dans parameter
!          04-06-10  ajout dataperwavefile
!          13-06-10  ajout ale_selected
! 2010.10  18-06-10  modif sur declaration de pss_w, uwind_t, vwind_t
! 2010.11  27-06-10  ajout tideweight    
!          04-07-10  ajout de tableau specialement au post-traitement de la maree
!          15-07-10  des variables declarees dans module_forcage sont
!                    maintenant dans module_principal
!          16-07-10  temobc et salobc deviennent temobc et salobc
!          08-08-10  ajout alpha
! 2010.12  26-08-10  ajout meteo_lonstr meteo_lonend meteo_londlt meteo_latstr    
!                    meteo_latend meteo_latdlt
!          14-09-10  fluxbar_u et fluxbar_v passent de real*4 à double precision
!                    nouveaux tableaux pour parametrisation houle/courant
!          15-09-10  ajout wave_year_min
!          03-10-10  ajout istr,jstr,iend,jend 
!                    ajout velobl et cie...
! 2010.13  09-10-10  ajout advsmin, advsmax, advsstp
!          11-10-10  ajout cellboxfactor1 et cellboxfactor2
!          04-11-10  ajout year_now month_now day_now hour_now minute_now second_now   
!          09-11-10  ajout upwindriver_t tfalpha tfa0 tfa1 et cie... 
! 2010.14  24-11-10  ajout riverupwdist
! 2010.15  27-12-10  ajout tfd0 tfd1 tfd2 tfd3
!                    ajout kbbc_u kbbc_v kbbc_w
! 2010.17  15-01-11  extension dimension kbbc_w
! 2010.18  28-01-11  ajout flag_asselin (0 ou 1)
!          03-02-11  ajout ncid2
!          22-02-11  Retour en arrirere par rapport à kbbc
!          23-02-11  Extension dimension kmin_w
! 2010.22  26-04-11  Ajout grid_volumeglb
!          01-05-11  Ajout rap1 rap2
!          06-05-11  Ajout hsw_wave et foc_wave
! 2010.23  21-05-11  Ajout tab1_code tab2_code tab_code_glb
! 2010.25  02-02-12  nouvelles declarations pour subroutine drifter
!          03-04-12  ajout title_for_netcdf_files_txt
!          31-05-12  dmensions de grille lues dans name_list
! S26      26-09-12  ajout ofl_rotation
!          02-11-12  ajout initialgridfile_txt,vertical_axis_convention
!          27-01-13  etendre les dimensions des tableaux de vitesse suite
!                    a advection qdm o4
!          01-02-13  ajout noise_velbar
!          01-03-13  ajout ssh_avr_w (moyenne de la ssh sur cycle barotrope)
!          09-03-13  ajout subcycle_exchange et synchro
!          01-05-13  ajout flag_externalmode_lf_or_fb
!          04-05-13  ajout type_structured type_unstructured
!          12-05-13  augmneter les dimensions de xy_t pour module_wave
!          02-06-13  ajout grid_i0 grid_j0
!          23-06-13  ajout mangrove
!          17-07-13  ajout sqr_hoverg_u sqr_hoverg_v
!          19-07-13  ajout h_inf_obc
!          09-10-13  ajout de temlwf et sallwf
!          31-10-13  ajout ofl_readtime
!          24-11-13  ajout ogcm_readtime_next, ..._prev
!          14-02-14  mpi_neighbor_list passe dans module_principal pour faire_hotrestart
!                    Ajout tkeb & epsb
!          02-05-15  ajout ofl_sshtype
!          15-06-14  logical fplan_grid
!          26-06-14  modif allocation 4eme arg tridia_in(0:4) 
!          03-07-14  modif allocation 4eme arg tridia_in(-1:4) 
!          11-07-14  ajout ogcm_time_shift
!          14-07-14  ajout kp2 km2
!          16-07-14  ajout offline_init_status                   
!          23-07-14  nomfichier allocate dans main.F90
!          25-07-14  ajout ofl_surflux
!          30-07-14  ajout flag_ssr24avr
!          21-08-14  ajout veldtrodx_u veldtody_v
!          11-09-14  ioffline_prv en tete de declaration
!          18-09-14  ajout divflux0_u divflux0_v
!          01-10-14  ajout anyv3dint
!          27-10-14  ajout meteo_t0
!          11-11-14  redimensionnement hssh et fluxbar
!          24-11-14  ajout flag_sequoia
!          17-12-14  ajout tableaux globaux de latitude longitude
!          31-01-15  ajout h0_w
!          17-02-15  ajout departuredate
!          05-03-15  ajout ofl_tke tkeofl_w dpt_wave_t     
!          01-07-15  ajout flag_p0m_filter
!______________________________________________________________________

      module module_principal
      use module_parameter
      use module_biology

      character*1 ::                                         &
        txtslash='/'

      character*4 ::                                         &
        meteo_t0='file'                                       !norestart

      character*5 ::                                         &
        ofl_type='short'                                     & !norestart
       ,ofl_sshtype='real'                                     !norestart

      character                                              &
        texte90*200                                          &
       ,texte250*250                                         & !13/04/06
       ,texte60*60                                           &
       ,texte3*3                                             &
       ,texte4*4                                             &
       ,texte11*11                                           &
       ,directory*200                                        &
       ,directory_offline*200                                &
!      ,txtslash*1                                           &
       ,texte30*30                                           &
       ,dirgraph*250                                         & !11-09-09
       ,variablesdesc*250                                    & !07/02/07
       ,extmodtype*6                                         &
       ,obc_ogcm_type*30                                     & !05-11-09
       ,txt_wavemodeltype*60                                 & ! 24-02-10 
       ,verticalgridfile*200                                 & !31-05-10
       ,sigstepgridfile*200                                  & !31-05-10
       ,flag_humspec_dewpoint*8                              & !26-08-10
       ,flag_meteodata*20                                    & !26-08-10
       ,liste_ogcm_t*60                                      &
       ,liste_ogcm_s*60                                      &
       ,liste_ogcm_u*60                                      &
       ,liste_ogcm_v*60                                      &
       ,liste_ogcm_ssh*60                                    &
       ,liste_ogcm_uvts*60                                   &
       ,liste_ogcm_grid*60                                   &
       ,liste_ogcm_grid_t*60                                 &
       ,liste_ogcm_grid_u*60                                 &
       ,liste_ogcm_grid_v*60                                 &
       ,dynrestartfilename*100                               &
       ,model_name*60                                        &
       ,attribute_netcdf_files*90                            &
       ,biorestartfilename*100                               &   
       ,drifter_initial_mode*8                               &
       ,meteo_sealand_file*80                                &
       ,initialgridfile_txt*100                              &
       ,title_for_netcdf_files_txt*30                        &
       ,vert_axis_conv_direc*2                               & ! norestart
       ,vert_axis_conv_start*1                               & ! norestart
       ,vert_axis_conv_end*1                                 & ! norestart
       ,hori_axis_conv_start*1                               & ! norestart
       ,hori_axis_conv_end*1                                 & ! norestart
       ,flag_externalmode_lf_or_fb*2                         &
       ,lonlatfile*100                                       &
       ,mangrove_file_name*100                               & ! norestart
       ,mergebathy_filename*100                              &
       ,biobc_nextfilename*200                               & 
       ,biobc_prevfilename*200                               &
       ,mpi_map_file_name*100                                &
       ,mpi_hole_plugging*100                                &
       ,subroutinetitle*40                                   &
       ,subroutinedescription*200                            &
       ,filename_runoff*200                                  &
       ,pgfscheme*3                                          &
       ,tmpdirname*30                                        & !norestart
       ,restartdirname*30                                    & !norestart
       ,flag_tke_stab_func*4                                   !norestart

      character(len=6) :: cl_model_name = 'symph1'
      character(len=4) :: cl_grid_name = 'sym1'

      character,dimension(:,:),allocatable :: &
       nametide*100                           &
      ,obcfile*90    
      character,dimension(:),allocatable :: &
       tidenodalfile*100                    &
      ,texte80*200                          &     
      ,grh_titre*30                         & 
      ,rivername*40                         &
      ,airseafile*90                        &
      ,riverfile*90                         &
      ,month*9                              &
      ,offlinefile*90                       & 
      ,nest_path_in*90                      &
      ,nest_path_out*90                     &
      ,obc_hf_file*90                       &
      ,bioriv_file*90                       &
      ,airsea_standard_name*100             &
      ,airseabinreclist*30                  &
      ,nomfichier*200

      logical fplan_grid

! debut module ne pas toucher à cette ligne
      integer           & 
       ioffline_prv     ! en tete de declaration pour "faire_dyn_restart.F"     

      double precision,dimension(:,:),allocatable :: &
       analysis3dmatrix
      double precision,dimension(:,:,:),allocatable :: &
       c_wave_mode                         &
      ,pcoefmode_t                         &
      ,ucoefmode_t                         &
      ,vcoefmode_t                         &
      ,rhrefmode_t
      double precision,dimension(:,:,:,:),allocatable :: &
       uv_wmode_t          &
      ,analysis_p3d_t      &
      ,analysis_u3d_t      &
      ,analysis_v3d_t      &
      ,p3dmode_cos_t       &
      ,p3dmode_sin_t       &
      ,u3dmode_cos_t       &
      ,u3dmode_sin_t       &
      ,v3dmode_cos_t       &
      ,v3dmode_sin_t          

      real,dimension(:,:,:,:,:),allocatable :: &
       bio_relax_north                         &
      ,bio_relax_south                         &
      ,bio_relax_east                          &
      ,bio_relax_west

      double precision,dimension(:,:,:),allocatable ::  &
       analysetide_u                                    &
      ,analysetide_v                                    &
      ,analysetide_w                                    &
      ,cocosisi_f                                       & !norestart             
      ,gridrotcos_f                                     & 
      ,gridrotsin_f                                     &
      ,temref_z                                         & 
      ,salref_z             
      double precision,dimension(:,:),allocatable :: &
       tideanalysismatrix                            &
      ,iaveraged_in                                  & !norestart
      ,iaveraged_out                                 & !norestart
      ,iaverag1d_in                                  & !norestart
      ,iaverag1d_out                                 & !norestart
      ,light_kpar2_w                                   !norestart

      double precision,dimension(:,:,:),allocatable :: &
       tidepotential_w                       &
      ,sshtotide_w                           &
      ,sshtidecos_w                          &
      ,sshtidesin_w                          &
      ,veltidecos_u                          &
      ,veltidesin_u                          &
      ,veltidecos_v                          &
      ,veltidesin_v                          &
      ,potidecos_w                           &
      ,potidesin_w                           &
      ,veltidecosout_u                       &
      ,veltidesinout_u                       &
      ,veltidecosout_v                       &
      ,veltidesinout_v                       &
      ,sshtidecosout_w                       &   
      ,sshtidesinout_w       

      double precision,dimension(:,:),allocatable :: &
       passetide                           &
      ,z1meatide                           &
      ,z2meatide                           &
      ,ftide                               &
      ,utide                               &
      ,northflux_sumt_v                    & 
      ,southflux_sumt_v      

      real*4,dimension(:,:),allocatable :: &
       drifter_l                               ! mpirestart

      real*4,dimension(:),allocatable :: &
       drifter_send_nord                 & ! mpirestart
      ,drifter_recv_nord                 & ! mpirestart
      ,drifter_send_sud                  & ! mpirestart
      ,drifter_recv_sud                  & ! mpirestart
      ,drifter_send_est                  & ! mpirestart
      ,drifter_recv_est                  & ! mpirestart
      ,drifter_send_ouest                & ! mpirestart
      ,drifter_recv_ouest                & ! mpirestart
      ,drifter_send_nordest              & ! mpirestart
      ,drifter_send_nordouest            & ! mpirestart
      ,drifter_send_sudest               & ! mpirestart
      ,drifter_send_sudouest             & ! mpirestart
      ,drifter_recv_nordest              & ! mpirestart
      ,drifter_recv_nordouest            & ! mpirestart
      ,drifter_recv_sudest               & ! mpirestart
      ,drifter_recv_sudouest               ! mpirestart

      double precision,dimension(:,:,:),allocatable :: &
       tkeb_w                                          &
      ,tken_w                                          &
      ,tkea_w                                          &
      ,tkle_w                                          &
      ,tkll_w                                          &
      ,epsb_w                                          &
      ,epsn_w                                          &
      ,epsa_w                                          &
      ,gradssh_u                                       &
      ,gradssh_v                                       &
      ,pgf_u                                           &
      ,pgf_v

      double precision,dimension(:,:,:,:),allocatable ::               &
       dz_u                          &
      ,dz_v                          &
      ,dz_t                           
      double precision,dimension(:,:,:),allocatable ::               &
       stokesforces_u                   &
      ,stokesforces_v                   &
      ,depth_t                          &
      ,sigma_w                          &
      ,sigma_f                          &
      ,dsig_u                           &
      ,dsig_v                           &
      ,dsig_t                           &
      ,dsig_f                           &
      ,depth_u                          &
      ,depth_v                          &
      ,depth_w                          &
      ,depth_f                          &
      ,km_w                             &
      ,kh_w                             &
      ,rhpref_u                         &
      ,temref_u                         &
      ,salref_u                         &
      ,rhpref_v                         &
      ,temref_v                         &
      ,salref_v                         &
      ,rhp_t                            &                       
      ,rhpref_t                         &   
      ,rhcref_t                         &   
      ,temref_t                         &   
      ,salref_t                          
      double precision,dimension(:,:,:,:),allocatable ::               &
       vel_u                         &
      ,vel_v                         &
      ,presgrad_u                    &
      ,presgrad_v                    &
      ,omega_w                       &
      ,veldydz_u                     &
      ,veldxdz_v                     &
      ,veldtodx_u                    &
      ,veldtody_v                    &
      ,tem_t                         &
      ,sal_t                         &
      ,tridia_in                     & !norestart
      ,hr_z2lr_w                     &
      ,anyv3d                          !norestart

!     double precision,dimension(:,:),allocatable :: &
!      noise_velbar_u                                &
!     ,noise_velbar_v                                &
!     ,divflux0_u                                    &
!     ,divflux0_v                                    &
!      rhpzavr_w                   !norestart

      double precision,dimension(:,:,:),allocatable ::               &
       tridia_out                            &
      ,fluxbar_sumt_u                        &
      ,fluxbar_sumt_v                        & 
      ,velbar_u                              &
      ,uflux2d_f                             &
      ,vflux2d_f                             &
      ,vlxbar_f                              &
      ,velbar_v                              &
      ,vlybar_f                              &
      ,flux2d_u                              &
      ,flux2d_v                              &
      ,sponge_u                              &
      ,sponge_v                              &
      ,ssh_w                                 &
!     ,sshe_w                                &
      ,hssh_f                                &
      ,hz_u                                  &
      ,hz_v                                  &
      ,hz_w                                  &
      ,sponge_t                              &
      ,zeroleveldepth_u                      &
      ,zeroleveldepth_v                      &
      ,heatrelax_w                           &
      ,xy_u                                  & !norestart
      ,xy_v                                  & !norestart
      ,xy_t                                  & !norestart
      ,xy_f                                  & !norestart
      ,pss_w                                 &
      ,uwind_t                               &
      ,vwind_t                               &
      ,uwind100_t                            &
      ,vwind100_t                            &
      ,wstress_u                             &
      ,wstress_v                             &
      ,sshf_w                                &
      ,slhf_w                                &
!     ,ssh_int_u                             &
!     ,ssh_int_v                             &
      ,ssh_int_w                             &
!     ,ssh_avr_w                             &
!     ,streamf_f                             &
      ,ssr_w                                 &
      ,ssr24prv_w                            &
      ,snsf_w                                &
      ,precipi_w                             &
      ,ustokesvortex_t                       &
      ,ustokesvortex_f                       &
      ,vstokesvortex_t                       &
      ,vstokesvortex_f                       &
      ,velavr_u                              &
      ,velavr_v                              &
      ,timefilter_u                          &
      ,timefilter_v                          &
!     ,grdobc3d_j                            &
!     ,grdobc3d_i                            &
      ,teta0_t                               &
      ,teta2_t                               &
      ,q2_t                                  &
      ,teta2delta_t                          &
      ,q2delta_t                             &
      ,uwinddelta_t                          &
      ,vwinddelta_t                          &
      ,tetastar_t                            &
      ,ustar_t                               &
      ,qstar_t                               &
      ,vmi_u                                 &
      ,vmi_v                                 &
      ,frozenterm2d_u                        &
      ,frozenterm2d_v                        &
      ,frozenterm3d_u                        &
      ,frozenterm3d_v                        &
      ,fluxbar_u                             &
      ,fluxbar_v       

      double precision,dimension(:,:),allocatable ::               &
       lon_t                                  &
      ,lat_t                                  &
      ,lon_u                                  &
      ,lat_u                                  &
      ,lon_v                                  &
      ,lat_v                                  &
      ,lon_f                                  &
      ,lat_f                                  &
      ,globlon_t                              & !norestart
      ,globlat_t                              & !norestart
      ,globlon_u                              & !norestart
      ,globlat_u                              & !norestart
      ,globlon_v                              & !norestart
      ,globlat_v                              & !norestart
      ,fluxbarsum_ofl_u & ! restart if ioffline_prv>=1
      ,fluxbarsum_ofl_v & ! restart if ioffline_prv>=1
      ,dxdy_u                                 &
      ,dxdy_v                                 &
      ,dxdy_t                                 &
      ,dxdy_e_t                               &
      ,dx_e_f                                 &
      ,dy_e_f                                 &
      ,dx_u                                   &
      ,dx_v                                   &
      ,dx_t                                   &
      ,dx_f                                   &
      ,dy_u                                   &
      ,dy_v                                   &
      ,dy_t                                   &
      ,dy_f                                   &
      ,h_w                                    &
      ,h0_w                                   &
      ,h_u                                    &
      ,h_v                                    &
      ,h_f                                    &
      ,coriolis_t                             &
      ,coriolis_f                             &
      ,wetmask_u                              &
      ,wetmask_v                              &
      ,wetmask_t                              &
      ,rhpzavr_w                              &  !norestart
      ,q10_t                                  &
      ,teta10_t                               &
      ,fric_u                                 &
      ,fric_v                                 &
      ,fric_t                                 &
      ,velbot3d2d_u                           &
      ,velbot3d2d_v                           &
      ,cdb_t                                  &
      ,cdb_f                                  &
      ,xflux_t                                &
      ,xflux_f                                &
      ,yflux_t                                & 
      ,yflux_f                                &
      ,pres3d2d_u                             &
      ,pres3d2d_v                             &
      ,adve3d2d_u                             &
      ,adve3d2d_v                             &
      ,restoring3d2d_u                              &
      ,restoring3d2d_v                              &
      ,stokesforces3d2d_u                     &
      ,stokesforces3d2d_v                     &
      ,wstress_w                              &
      ,z0_w                                   &
      ,albedo_w                               &
      ,gridrotcos_t                           &
      ,gridrotsin_t                           &
      ,sshstokes_w                            &
      ,cwi_int_u                              &
      ,cwi_int_v                              &
      ,cwj_int_u                              &
      ,cwj_int_v                              &
      ,sqr_hoverg_u                           &
      ,sqr_hoverg_v                           &
      ,sshrefobc_i                            &
      ,vbrrefobc_i                            &
      ,sshrefobc_j                            &
      ,vbrrefobc_j                            &
      ,anyv1d                                 &
      ,gridcard                               &
      ,novector                               &
      ,proficard                              &
      ,airseainfo                             &
      ,airseadt                               &
      ,riverdt                                &
      ,river_t                                &
      ,riverflux                               

      integer,dimension(:,:,:),allocatable :: &
       mask_t       &
      ,mask_f       &
      ,mask_u       &
      ,mask_v       & 
      ,lonlat2ij_t  & ! norestart
      ,sodate         ! norestart

      integer,dimension(:,:),allocatable :: &
       kmin_u                               &
      ,kmin_v                               &
      ,kmin_w                               &
      ,ksl_t                                &
      ,glob_mask_mangrove                   &
      ,mask_mangrove_t                       
 

      real*4,dimension(:,:),allocatable :: &
       upwindriver_t          ! norestart 

      real*4,dimension(:,:,:,:),allocatable :: &
       temobc_t         &
      ,salobc_t         &
      ,velobc_u         &
      ,velobc_v         &
      ,temofl_t         & ! restart if ioffline_prv>=1
      ,salofl_t         & ! restart if ioffline_prv>=1
      ,dzofl_t          & ! restart if ioffline_prv>=1
      ,velofl_u         & ! restart if ioffline_prv>=1
      ,velofl_v         & ! restart if ioffline_prv>=1
      ,tkeofl_w         & ! restart if ioffline_prv>=1
      ,dfvofl_w         & ! restart if ioffline_prv==1
      ,uwindabl_t       &
      ,vwindabl_t

      real*4,dimension(:,:,:),allocatable :: &
       w0mofl_w           & ! restart if ioffline_prv==1
      ,ablheight_t        &
      ,wwindabl_w         &
      ,kz_abl_w

      real*4,dimension(:,:,:,:),allocatable ::                     &
       velstokes_u                                                 &
      ,velstokes_v

      real*4,dimension(:,:,:),allocatable ::  &
       velbarstokes_u                         &
      ,velbarstokes_v                         &
      ,nhp1_t                                 &
      ,nhp2_t                                 &
      ,temlwf_t                               &
      ,sallwf_t    

      real*4,dimension(:,:,:),allocatable ::                       &
       t_wave_t                                                    &
      ,hs_wave_t                                                   &
      ,hsw_wave_t                                                  &
      ,foc_wave_t                                                  &
      ,k_wave_t                                                    &
      ,kx_wave_t                                                   &
      ,ky_wave_t                                                   &
      ,twox_wave_t                                                 &
      ,twoy_wave_t                                                 &
      ,tawx_wave_t                                                 &
      ,tawy_wave_t                                                 &
      ,usf_wave_t                                                  &
      ,vsf_wave_t                                                  &  
      ,dir_wave_t                                                  &
      ,uss_wave_t                                                  &
      ,j_wave_t                                                    &
      ,vss_wave_t

      real*4,dimension(:,:),allocatable ::                         &
       ubw                                                         &
      ,fw                                                          &
      ,dpt_wave_t                                                  &
      ,wstresb_u                                                   &
      ,wstresb_v                                                   &
      ,ij2ww3_i                                                    &
      ,ij2ww3_j                                                    &
      ,ij2ww3_teta

      real*4,dimension(:,:),allocatable ::                         &
       slhf_aver_w                                                 &
      ,sshf_aver_w                                                 &
      ,snsf_aver_w                                                 &
      ,ssr_aver_w                                                  &
      ,precipi_aver_w                                              &
      ,wstress_aver_u                                              &
      ,wstress_aver_v


      real,dimension(:,:,:),allocatable ::   &
       anyvar3d           !norestart                         

      integer,dimension(:,:,:),allocatable ::   &
       anyv3dint           !norestart                         

      real*4,dimension(:,:,:),allocatable :: &
       sshobc_w                              &
      ,velbarobc_u                           &
      ,velbarobc_v                           &
      ,sshofl_w         & ! restart if ioffline_prv>=1
      ,velbarofl_u      & ! restart if ioffline_prv>=1
      ,velbarofl_v        ! restart if ioffline_prv>=1

      real,dimension(:,:),allocatable ::                     &
       anyvar2d           !norestart

      integer,dimension(:),allocatable ::                    &
       riverdir                                              & ! norestart
      ,l_river                                               & ! norestart
      ,river_no                                              & ! norestart
      ,rivertrc_inout                                        & ! norestart
      ,rivervel_inout                                        & ! norestart
      ,riverupwdist                                          & ! norestart
      ,nest_len_in                                           & ! norestart
      ,jpm                                                     ! norestart

      integer,dimension(:,:),allocatable ::                  &
       iriver                                                &
      ,jriver

      integer,dimension(:),allocatable ::                           &
       spo_i1_u                                                     & ! norestart
      ,spo_i2_u                                                     & ! norestart
      ,spo_j1_u                                                     & ! norestart
      ,spo_j2_u                                                     & ! norestart
      ,spo_i1_v                                                     & ! norestart
      ,spo_i2_v                                                     & ! norestart
      ,spo_j1_v                                                     & ! norestart
      ,spo_j2_v                                                     & ! norestart
      ,spo_i1_t                                                     & ! norestart
      ,spo_i2_t                                                     & ! norestart
      ,spo_j1_t                                                     & ! norestart
      ,spo_j2_t                                                     & ! norestart
      ,obcstreamf                                                   & ! norestart
      ,mask_i_w                                                     & ! norestart
      ,mask_i_v                                                     & ! norestart
      ,mask_i_u                                                     & ! norestart
      ,mask_j_w                                                     & ! norestart
      ,mask_j_u                                                     & ! norestart
      ,mask_j_v                                                       ! norestart

!     integer,dimension(:,:),allocatable :: &
!      l2ij_out_u                           &
!     ,l2ij_out_v                           &
!     ,l2ij_out_w                           &
!     ,l2ij_in_u                            &
!     ,l2ij_in_v                            &
!     ,l2ij_in_w                            &
!     ,fillmask_t

      integer,dimension(:),allocatable :: &
       dayindex                           &
      ,grh_out_var                        &  ! norestart
      ,varid                              &  ! norestart
      ,grh_nb                             &
      ,vardim                             &  ! norestart
      ,varstart                           &  ! norestart
      ,varcount                           &  ! norestart
      ,departuredate

      integer,dimension(:,:),allocatable :: &
       datesim                              & 
      ,dateobc                              &
      ,dateairsea                           &
      ,dateriver

      real*4,dimension(:),allocatable :: &
       cgridshift                             ! norestart

      real*4,dimension(:),allocatable ::   &
        runoff_dt
      real*4,dimension(:,:),allocatable :: &
        runoff_w                            

      double precision,dimension(:),allocatable ::   &
       frqtide                             &
      ,v0tide                              &
      ,ti0tide                             &
      ,equitide                            &
      ,tideanaweight                       &
      ,sshmeatide                          &
      ,airseafile_prvtime                  &
      ,airseafile_nextime                  & 
      ,ogcm_readtime_next                  &
      ,ogcm_readtime_prev                  &
      ,ogcm_period_prev                    &
      ,ogcm_period_next                    &
      ,timeweightobc                       &
      ,dtobc                               &
      ,dt_abl_max                          &
      ,zref_z                              &
      ,checkanyv3d                           ! norestart

      integer,dimension(:),allocatable ::   &
       tab1_code                      &   ! norestart
      ,tab2_code                      &   ! norestart
      ,tab_code_glb                   &   ! norestart
      ,ub2                            &   ! norestart
      ,lb2                            &   ! norestart
      ,ub3                            &   ! norestart
      ,lb3                            &   ! norestart
      ,ub4                            &   ! norestart
      ,lb4                            &   ! norestart
      ,nutide                         &
      ,nest_len_out                   &
      ,ind                            &
      ,airseafile_prvtrec             &
      ,airseafile_nextrec             &      
      ,mpi_neighbor_list              &
      ,ogcm_rec_next                  &
      ,ogcm_rec_prev                  &
      ,obcstatus                         ! norestart             

      real*4 ::                          &
       drifter_out_sampling              &  ! norestart
      ,obc2dtype                         &  ! norestart
      ,coef_diss_mangrove                &  ! norestart 
      ,expnum              =1.           &  ! norestart 
      ,cst_c0cub                         &  ! norestart 
      ,relativewind        =0.              ! norestart 
 
      integer,dimension(:),allocatable :: &
       drifter_send_order_nord            & ! mpirestart
      ,drifter_send_order_sud             & ! mpirestart
      ,drifter_send_order_est             & ! mpirestart
      ,drifter_send_order_ouest           & ! mpirestart
      ,drifter_send_order_out             & ! mpirestart
      ,drifter_send_order_nordest         & ! mpirestart
      ,drifter_send_order_nordouest       & ! mpirestart
      ,drifter_send_order_sudest          & ! mpirestart
      ,drifter_send_order_sudouest        & ! mpirestart
      ,zw_abl                               ! norestart

      integer                                                        &
       drifter_onoff       ! norestart

      double precision,dimension(:),allocatable ::   &
       river_s                                       &
      ,river_tmin                                    &
      ,river_tmax                                    &
      ,hriver                                        &
      ,friver                                        &
      ,realriver                                     &
      ,riverinfo                                     &
      ,river_timeref                                 &
      ,daysim                                        &
      ,tdate_output                                  &  !norestart
      ,obcinfo                                       &
      ,offlinedt                                     &
      ,pss_mean                                      &
      ,sshobc_cum                                    &
      ,nest_dt_in                                    &
      ,nest_dt_out                                   &
!     ,obcafdt                                       &
      ,wavedt                                        &
      ,obc_hf_dt                                     &
      ,glob_dte_lp                                   &
      ,glob_dte_lp_tmp

      integer ::                 &
       ww3_varmax                &
      ,ww3_type_grid             &
      ,type_unstructured         &
      ,type_structured           &
      ,vststep                =0 & ! norestart
      ,zprofile2d3d           =0 & ! norestart
      ,timestep_type          =1 & ! norestart
      ,timestep_leapfrog      =0 & ! norestart
      ,timestep_forwbckw      =1 & ! norestart
      ,bulk_core              =1 & ! norestart
      ,bulk_moon              =2 & ! norestart
      ,bulk_coare             =3 & ! norestart
      ,bulk_scheme            =1   ! norestart

      real*4,dimension(:),allocatable ::                           &
       freq_wave                   ! restartbin 

      double precision :: &
!      ti0tide            &
       tideana_spinup     & ! norestart
      ,tideana_delta      & ! norestart
      ,sshtide_cum        & 
      ,albedo_val         & ! norestart
      ,tfilterfb  =0.01   & ! norestart
      ,tfilterlf  =0.01     ! norestart

      integer ::                      &
       ktide                          &
      ,tideforces                     &
      ,tideana_yesno                  &   ! norestart
      ,ioffline                       &   ! norestart
      ,tideanalysis_count  =0         &   ! norestart
      ,ibl1_advbio         =0         &   ! norestart              
      ,ibl2_advbio         =0         &   ! norestart              
      ,jbl1_advbio         =0         &   ! norestart              
      ,jbl2_advbio         =0         &   ! norestart              
      ,ogcm_time_shift                &
      ,ieq1                           & 
      ,ieqimax                        & 
      ,jeq1                           & 
      ,jeqjmax                        &
      ,ieq1_jeq1                      &
      ,ieq1_jeqjmax                   & 
      ,ieqimax_jeq1                   &
      ,ieqimax_jeqjmax                &
      ,dim_varid                      &   ! norestart
      ,albedo_constant   =0           &   !notrestart
      ,albedo_apel1987   =1           &   !notrestart
      ,albedo_br1982     =2           &   !notrestart
      ,albedo_br1986     =3           &   !notrestart
      ,albedo_case       =-9999       &   !notrestart
      ,il_oasis_time                      !notrestart

      integer ::                   &
       loop1                       &
      ,loop2                       &
      ,loop3                       &
      ,nairsea                     &
      ,bi_onoff                    &
      ,nc_or_bin_airsea            &
      ,loop_netcdf                 &
      ,count_netcdfvar             &
      ,wavefile_prvtrec            &
      ,wavefile_nextrec            &
      ,ofl_rec_now                 &
      ,ofl_rec_max                 &
      ,flag_kz_enhanced    =0      & !norestart
      ,flag_dt_adjust      =0      & !norestart
      ,deltamax                    &
      ,tide_interpolation  =0      & !norestart
      ,tide_flagrotation   =0      & !norestart
      ,flag_ssr24avr       =0      & !norestart
      ,flag_abl            =0      & !norestart 
      ,flag_abl2           =0      & !norestart 
      ,flag_sequoia        =0      & !norestart 
      ,dimssr24prv                 &
      ,flag_meteo_average  =0      & !norestart
      ,flag_nemoffline     =0      & !norestart
      ,trc_id              =1      & !norestart
      ,vel_id              =2      & !norestart
      ,ssh_id              =3      & !norestart
      ,var_num             =3      & !norestart
      ,flag_p0m_filter     =0      & !norestart
      ,flag_refstate       =0      &
      ,flag_linearfric     =0      & !norestart
      ,flag_1dv            =0      & !norestart
      ,flag_smooth_h_mask  =0        !norestart

      double precision ::                                            &
       var_misval                                                    &
      ,un_r8                                                         &
      ,heure                                                         &
      ,suma                                                          &
      ,sumb                                                          &
      ,sum_mi                                                        &
      ,grid_area                                                     &
      ,grid_areaglb                                                  &
      ,grid_volumeglb                                                &
      ,sumarchive                                                    &
      ,lonmin                                                        &
      ,lonmax                                                        &
      ,latmin                                                        &
      ,latmax                                                        &
      ,wavefile_prvtime                                              &
      ,wavefile_nextime                                              &
      ,ofl_period_prev                                               &
      ,ofl_period_now                                                &
      ,ofl_period_next                                               &
      ,ofl_nextrec_time                                              &
      ,ofl_writime                                                   &
      ,ofl_readtime_next                                             &
      ,ofl_readtime_prev                                             &
      ,biobc_nextfiletime                                            &
      ,biobc_prevfiletime                                            & 
      ,stability_index                                               &
      ,iteration2d_max_r8                                            &
      ,iteration2d_upbound                                           & !norestart
      ,coef_linearfric          =0.                                  & !norestart
      ,gh                                                            & !norestart
      ,gm                                                            & !norestart
      ,nn                                                            & !norestart
      ,tke2overeps                                                   & !norestart
      ,tkeovereps                                                    & !norestart
      ,gravoverrho                                                   & !norestart
      ,sh                                                            & !norestart
      ,sm                                                            & !norestart
      ,heatfluxbias             =0.                                    !norestart

      double precision ::                                             &
       zero                                                           &
      ,un                                                             &
      ,cdseuil                                                        & !norestart
      ,cdb_2dh                                                        & !norestart
      ,small1                                                         &
      ,deci                                                           &
      ,decj                                                           &
      ,deck                                                           &
      ,rap1                                                           &
      ,rap2                                                           &
      ,rapi                                                           &
      ,rapj                                                           &
      ,rapk                                                           &
      ,rap                                                            &
      ,const0                                                         & !norestart
      ,const1                                                         & !norestart
      ,const2                                                         & !norestart
      ,const3                                                         & !norestart
      ,const4                                                         & !norestart
      ,const5                                                         & !norestart
      ,const6                                                         & !norestart
      ,const7                                                         & !norestart
      ,const8                                                         & !norestart
      ,const9                                                         & !norestart
      ,uv_10                                                          &
      ,karman                                        =0.4             & !norestart
      ,stefan                                        =5.67032e-8      & !norestart
      ,z2m                                           =2.              & !norestart
      ,z10m                                          =10.             & !norestart
      ,pss0                                          =1.e5            & !norestart
      ,boltz                                         =1.380658e-23    & !norestart
      ,planck                                        =6.6260755e-34   & !norestart
      ,avogadro                                      =6.0221367e+23   & !norestart
      ,ce                                                             &! CE is the Dalton number.
      ,cen                                                            &! cas de la neutralité
      ,ch                                                             &! CH is the Stanton number
      ,chn                                                            &! cas de la neutralité
      ,cd                                                             &! CD is the drag coefficient
      ,cdn                                                            &! cas de la neutralité
      ,z_2                                                            &! niveau 2  mètres
      ,z_10                                                           &! niveau 10 mètres
      ,q_0                                                            &! humidité spécifique à  0 mètres (surface mer)
      ,r_0                                                            &! rapport de mélange à 0 mètres (surface mer)
      ,pvs_0                                                          &! Pression Vapeur Saturante 0 metre
      ,psih_10                                                        &! fonctions universelles psi à 10m
      ,psih_2                                                         &! fonctions universelles psi à 2m
      ,phim                                                           &! fonctions universelles phi
      ,zl_10                                                          &! rapport z sur L à 10m
      ,zl_2                                                           &! rapport z sur L à  2m
      ,falpha                                                         &! coef fonction universelle instable (1-ALPHA Z/L)
      ,fbeta                                                          &! coef fonction universelle   stable (1+BETA  Z/L)
      ,ro                                                             &! densite de l'air
      ,cp_air                                                         &!C capacité calorifique de l'air
      ,lv                                                             &
      ,psim_10                                                        &
      ,psim_2                                                         &
      ,dte_lp                                                         &
      ,dte_fw                                                         &
      ,dti_lpbef                                                      &
      ,dti_lp                                                         &
      ,dti_lpmax                                                      &
      ,dti_lpsub                                                      &
      ,dti_fwsub                                                      &
      ,dti_fw                                                         &
      ,dti_bef                                                        &
      ,dti_now                                                        &
      ,dtiratio                                                       &
      ,fbtfiltercoef                                                  &
      ,assel0                                                         &!C coef d'asselin n0
      ,assel1                                                         &!C coef d'asselin n1
      ,assel2                                                         &!C coef d'asselin n2
      ,assel3                                                         &!C coef d'asselin n2
      ,wetdry_cst1                                                    &! norestart
      ,wetdry_cst2                                                    &! norestart
      ,wetdry_cst3                                                    &! norestart
      ,h_inf                                                          &
      ,h_inf_obc                                                      &
      ,h_sup                                                          &
      ,dist                                                           &
      ,dist0                                                          &
      ,dist1                                                          &
      ,dist2                                                          &
      ,dist3                                                          &
      ,sponge_l                                                       &
      ,t0surf                                                         &
      ,s0surf                                                         &
      ,deg2rad                                                        &
      ,rad2deg                                                        &
      ,zmin                                                           &
      ,zmax                                                           &
      ,coa                                                            &
      ,cob                                                            &
      ,coc                                                            &
      ,sol1                                                           &
      ,sol2                                                           &
      ,discri                                                         &
      ,profm1                                                         &
      ,profp1                                                         &
      ,hmax                                                           &
      ,hstepmin                                                   &
      ,hstepmax                                                   &
      ,grav                                                       &
      ,cfl_sshmax                                                 &    ! norestart
      ,cfl_umax                                                   &    ! norestart
      ,cfl_reduce                                                 &    ! norestart
      ,relax_es                                                   &
      ,relax_ext                                                  &
      ,relax_int                                                  &
      ,relax_ts                                                   &
      ,relax_bpc                                                  &    ! norestart
      ,z0s                                                            &
      ,z1                                                             &
      ,z2                                                             &
      ,z3                                                             &
      ,z4                                                             &
      ,z5                                                             &
      ,y0                                                             &
      ,y2                                                             &
      ,y3                                                             &
      ,y4                                                             &
      ,y5                                                             &
      ,x0                                                             &
      ,x1                                                             &
      ,x2                                                             &
      ,x3                                                             &
      ,x4                                                             &
      ,x5                                                             &
      ,x6                                                             &
      ,x7                                                             &
      ,x8                                                             &
      ,x9                                                             &
      ,x10                                                            &
      ,x11                                                            &
      ,x12                                                            &
      ,x13                                                            &
      ,x14                                                            &
      ,x20                                                            &
      ,x21                                                            &
      ,x22                                                            &
      ,x33                                                            &
      ,x44                                                       

! Constants involved in Moon Scheme
      double precision ::       &
       xmd                      & ! norestart
      ,xmv                      & ! norestart
      ,xrd                      & ! norestart
      ,xrv                      & ! norestart
      ,xcpd                     & ! norestart
      ,xcpv                     & ! norestart
      ,xcl                      & ! norestart
      ,celsius2kelvin   =273.15 & ! norestart
      ,xlvtt                    & ! norestart
      ,xestt                    & ! norestart
      ,xgamw                    & ! norestart
      ,xbetaw                   & ! norestart
      ,xalpw                    & ! norestart
      ,zrvsrdm1                 & ! norestart
      ,z0_u                     & ! norestart
      ,qsat_sea_z               & ! norestart
      ,airdensity               & ! norestart
      ,sst_kelvin               & ! norestart
      ,prs_atm_z                & ! norestart   
      ,tem_atm_z                & ! norestart   
      ,exner_atm_z              & ! norestart   
      ,delta_u                  & ! norestart    
      ,delta_t                  & ! norestart    
      ,delta_q                  & ! norestart    
      ,delta_u_n                & ! norestart    
      ,exner_sea_z              & ! norestart
      ,psifunctt                & ! norestart
      ,psifunctu                & ! norestart
      ,z0_q                     & ! norestart
      ,z0_t                     & ! norestart
      ,sst1000hpa_kelvin        & ! norestart
      ,visa                     & ! norestart
      ,charnock        =0.011   & ! norestart
      ,qsat_atm_z               & ! norestart   
      ,ustar_bef                &
      ,qstar_bef                &
      ,tetastar_bef

      double precision                                               &
       rayonterre                                                    &
      ,northpole_lon                                                 &
      ,northpole_lat                                                 &
      ,southpole_lon                                                 &
      ,southpole_lat                                                 &
      ,phi0                                                          &
      ,longi                                                         &
      ,latit                                                         &
      ,longi0                                                        &
      ,latit0                                                        &
      ,longi1                                                        &
      ,latit1                                                        &
      ,angle0                                                        &
      ,alp_t                                                         &
      ,alp_s                                                         &
      ,pi                                                            &
      ,rho                                                           &
      ,rhoair                                                        &
      ,valmax                                                        &
      ,vis                                                           &
      ,tkee1                                                         &
      ,tkee2                                                         &
      ,tkee3                                                         &
      ,tkeg2                                                         &
      ,tkeg3                                                         &
      ,tkeg4                                                         &
      ,tkeg5                                                         &
      ,tkeg6                                                         &
      ,tkeb1                                                         &
      ,tkeb2                                                         &
      ,ctke1                                                         &
      ,ctke2                                                         &
      ,t0                                                            &
      ,s0                                                            &
      ,cp                                                            &
      ,light_kpar1           & !norestart
      ,light_att1            & !norestart
      ,light_att2            & !norestart
      ,light_rat1            & !norestart
      ,light_rat2            & !norestart
      ,light_att2_val1       & !norestart
      ,light_att2_h1         & !norestart
      ,light_att2_val2       & !norestart
      ,light_att2_h2         & !norestart
      ,z0b                                                           &
      ,t0_base                                                       &
      ,s0_base                                                       &
      ,rho_base                                                      &
      ,alp_t_base                                                    &
      ,alp_s_base                                                    &
!     ,tfalpha                                                       &
!     ,tfa0                                                          &
!     ,tfa1                                                          &
!     ,tfa2                                                          &
!     ,tfa3                                                          &
      ,tfb0                                                          &
!     ,tfc0                                                          &
!     ,tfc1                                                          &
!     ,tfc2                                                          &
!     ,tfc3                                                          &
!     ,tfd0                                                          &
!     ,tfd1                                                          &
!     ,tfd2                                                          &
!     ,tfd3                                                          &
!     ,tfd4                                                          &
      ,meteo_lonmin                                                  &
      ,meteo_latmin                                                  &
      ,meteo_lonmax                                                  &
      ,meteo_latmax                                                  &
      ,meteo_resol                                                   &
      ,meteo_resol_u                                                 &
      ,meteo_resol_v                                                 &
      ,meteo_lonstr                                                  &
      ,meteo_lonend                                                  &
      ,meteo_londlt                                                  &
      ,meteo_latstr                                                  &
      ,meteo_latend                                                  &
      ,meteo_latdlt                                                  &
      ,var_lonmin                                                    &
      ,var_latmin                                                    &
      ,var_lonmax                                                    &
      ,var_latmax                                                    &
      ,ww3_lonmin                                                    &
      ,ww3_latmin                                                    &
      ,ww3_lonmax                                                    &
      ,ww3_latmax                                                    &
      ,ww3_dlon                                                      &
      ,ww3_dlat                                                      &
      ,tide_lonmin                                                   &
      ,tide_latmin                                                   &
      ,tide_lonmax                                                   &
      ,tide_latmax                                                   &
      ,tide_dlon                                                     &
      ,tide_dlat                                                     &
      ,dxb                                                           &
      ,dyb                                                           &
      ,dxa                                                           &
      ,dya                                                           &
!     ,dt_obc                                                        &
      ,hmin                                                          &
      ,h1d                                                           &
      ,dlon                                                          &
      ,dlat                                                          &
      ,epsi                                                          &
      ,lagrange_ssh                                                  &
      ,diffu                                                         &
      ,difnorm                                                       &
      ,lup                                                           &
      ,ldown                                                         &
      ,zup                                                           &
      ,zdown                                                         &
      ,rbase                                                         &
      ,small                                                         &
      ,small3                                                        &
      ,hgesig                                                        &
      ,pgesig                                                        &
      ,windfactor                                                    &
      ,xdtk_out                                                      &
      ,tfond                                                         &
      ,sfond                                                         &
      ,rfond                                                         &
      ,c1streamf                                                     &
      ,c2streamf                                                     &
      ,rampe                                                         &
      ,rampe_wind                                                    &
      ,y1                                                            &
      ,obc_hf_reset                                                  &
      ,rho_0d                                                        &
      ,tem_0d                                                        &
      ,sal_0d                                                        &
      ,rho_tmp                                                       &
      ,cst_adv_hor                                                   & !norestart
      ,cst_adv_ver                                                   & !norestart
      ,cst_adv_vel                                                   & !norestart
      ,ssh_avr_nest_out                                              &
      ,graphperiod                                                   & !norestart
      ,graph_nextime                                                 & 
      ,restartfileperiod                                             & !norestart
      ,tidenodal_prev_rdv                                            &
      ,tidenodal_next_rdv                                            &
      ,tideana_modulo                                                & ! norestart
      ,tideana_nextime                                               & 
      ,cellboxfactor1                                                &
      ,cellboxfactor2                                                &
      ,rap_wave

      integer ::                                                    &
       ihmax                                                        &
      ,jhmax                                                        &
      ,ihmin                                                        &
      ,jhmin                                                        &
      ,convect_yn                                                   &
      ,nbvstepmin                                                   &
      ,bio_relax_size

      real ::                  &
       unit_r4                 &
      ,x1_r4                   &
      ,x2_r4                   &
      ,x3_r4                   &
      ,discharge               &
      ,filval                  & 
      ,var_validmin            &                         
      ,var_validmax            &
      ,vdw_loc                 &
      ,vup_loc                 &
      ,zdw_loc                 &
      ,var_scalefactor         &
      ,inv_scalefactor         &
      ,var_addoffset           &
      ,zup_loc                 &
      ,hrmax                   &
      ,relax_bio               & !norestart
      ,kmol_m                  & !norestart
      ,kmol_h                  & !norestart
      ,kmol_s                  & !norestart
      ,inv_ekman_depth         & !norestart
      ,constant_kz               !norestart

      integer ::                                                    &
       i                                                            &
      ,j                                                            &
      ,k                                                            &
      ,flag3d                                                          &
      ,lrec                                                         &
      ,iteration3d                                                  &
      ,compt1                                                       &
      ,compt2                                                       &
      ,compt3                                                       &
      ,compt4                                                       &
      ,kount0                                                       &
      ,kount1                                                       &
      ,kount2                                                       &
      ,kount3                                                       &
      ,kount4                                                       &
      ,kount5                                                       &
      ,kount6                                                       &
      ,kount7                                                       &
      ,kount8                                                       &
      ,kount9                                                       &
      ,kountmod                                                     &
      ,kountrdv1                                                    &
      ,kountrdv2                                                    &
      ,kountrdv3                                                    &
      ,kountrdv4                                                    &
      ,substep_advbio                                               & ! norestart
      ,subcycle_exchange                                            &
      ,subcycle_onoff                                               &
      ,subcycle_synchro                                             &
      ,subcycle_modulo                                              &
      ,iteration3d_restart                                          &
      ,filvalshort
      integer                      &
           status                  & ! norestart
          ,forcedstatus            & ! norestart
          ,decision                & ! norestart
          ,ncid1                   & ! norestart
          ,ncid2                   & ! norestart
          ,dim_x_id                & ! norestart
          ,dim_y_id                & ! norestart
          ,dim_z_id                & ! norestart
          ,dim_t_id                & ! norestart
          ,dim_b_id                & ! norestart
          ,max_x                   & 
          ,max_y                   & 
          ,max_z                   & 
          ,max_time_counter        &
          ,max_meteo_time_counter  &
          ,var_id                  &
          ,var_nftype              &
          ,meteo_imax              &
          ,meteo_jmax              &
          ,meteo_kmax              &
          ,meteozoom_istr          &
          ,meteozoom_iend          &
          ,meteozoom_jstr          &
          ,meteozoom_jend          &
          ,meteofull_imax         &
          ,meteofull_jmax        &
          ,tide_imax            &
          ,tide_jmax           &
          ,tide_kmax           &
          ,tidezoom_istr       & !norestart
          ,tidezoom_iend       & !norestart
          ,tidezoom_jstr       & !norestart
          ,tidezoom_jend       & !norestart
          ,tidezoom_istr_t     & !norestart
          ,tidezoom_iend_t     & !norestart
          ,tidezoom_jstr_t     & !norestart
          ,tidezoom_jend_t     & !norestart
          ,tidezoom_istr_u     & !norestart
          ,tidezoom_iend_u     & !norestart
          ,tidezoom_jstr_u     & !norestart
          ,tidezoom_jend_u     & !norestart
          ,tidezoom_istr_v     & !norestart
          ,tidezoom_iend_v     & !norestart
          ,tidezoom_jstr_v     & !norestart
          ,tidezoom_jend_v     & !norestart
          ,tidefull_imax       &
          ,tidefull_jmax       &
          ,ww3_imax            &
          ,ww3_jmax            &
          ,ww3_kmax            &
          ,ww3_fmax            &
          ,ww3zoom_istr        &
          ,ww3zoom_iend        &
          ,ww3zoom_jstr        &
          ,ww3zoom_jend        &
          ,ww3full_imax        &
          ,ww3full_jmax
      integer ::                                                       &
           jour                &
          ,kstop               & !norestart
          ,istr                & !norestart
          ,jstr                & !norestart
          ,kstr                & !norestart
          ,tstr                & !norestart
          ,bstr                & !norestart
          ,iend                & !norestart
          ,jend                & !norestart
          ,kend                & !norestart
          ,tend                & !norestart
          ,bend                & !norestart
          ,dimend               &
          ,kpvwave               &
          ,give_chanel9           &
          ,i0                      &
          ,j0                       &
          ,iteration2d        =1     &
          ,iteration2d_max_now        &
          ,iteration2d_max_bef         &
          ,i2dh                         &
          ,l1                            &
          ,l1_sca                         &
          ,l1_vec                          &
          ,l2                               &
          ,l2_sca                            &
          ,l2_vec                             &
          ,l3                                  &
          ,l3_sca                              &
          ,l3_vec                              &
          ,len1
      integer ::                                                        &
           nc                                                           &
          ,nc1                                                          &
          ,initial                                                      &
          ,itimets                                                      &
          ,itimebio                                                     &
          ,iadvec_ts              &     ! norestart
          ,iturbulence            &     ! norestart
          ,istreamf                                                     &
          ,itime                                                        &
          ,kts                                                          &
          ,kuv                                                          &
          ,ko                                                           &
          ,jm1                                                          &
          ,jm2                                                          &
          ,jp1                                                          &
          ,jp2                                                          &
          ,im1                                                          &
          ,im2                                                          &
          ,ip1                                                          &
          ,ip2                                                          &
          ,iwind                                                        &
          ,ip                                                           &
          ,im                                                           &
          ,jp                                                           &
          ,jm                                                           &
          ,kp                                                           &
          ,km                                                           &
          ,ipnoc
      integer ::                  &
           i1                     & !norestart
          ,i2                     & !norestart
          ,i3                     & !norestart
          ,i4                     & !norestart
          ,i5                     & !norestart
          ,i6                     & !norestart
          ,i7                     & !norestart
          ,i8                     & !norestart
          ,i9                     & !norestart
          ,i10                    & !norestart
          ,j1                     & !norestart
          ,j2                     & !norestart
          ,j3                     & !norestart
          ,j4                     & !norestart
          ,j5                     & !norestart
          ,j6                     & !norestart
          ,j7                     & !norestart
          ,j8                     & !norestart
          ,j9                     & !norestart
          ,k0                     & !norestart
          ,k1                     & !norestart
          ,k2                     & !norestart
          ,k3                     & !norestart
          ,k4                     & !norestart
          ,k5                     & !norestart
          ,k6                     & !norestart
          ,k7                     & !norestart
          ,k8                     & !norestart
          ,k9                     & !norestart
          ,kr                     & !norestart
          ,kp1                    & !norestart           
          ,kp2                    & !norestart 
          ,km1                    & !norestart 
          ,km2                    & !norestart
          ,key                    & !norestart
          ,flag                     !norestart
      integer ::                                                     &
           itest                                                     &
          ,itest1                                                    &
          ,itest2                                                    &
          ,itest3                                                    &
          ,ioption                                                   &
          ,nsmooth                                                   &
          ,nriver                                                    &
          ,ian0                                                      &
          ,imois0                                                    &
          ,ijour0                                                    &
          ,iheure0                                                   &
          ,iminute0                                                  &
          ,iseconde0                                                 &
          ,iref                                                      &
          ,jref                                                      &
          ,lname1                                                    &
          ,lname2                                                    &
          ,lname3                                                    &
          ,lname4                                                    & !C13/04/06
          ,kmode                                                     &
          ,kmodemax                                                  &
          ,fgrid_or_wgrid                                            &
          ,fgrid_case                                                &
          ,wgrid_case                                                &
          ,typegrid                                                  &
          ,typegrid_monopole                                         &
          ,typegrid_file                                             &
          ,typegrid_bipole 
      integer ::                                                        &
           istar,                                                       &
           jstar,                                                       &
           istop,                                                       &
           jstop,                                                       &
           izoomin,                                                     &
           izoomax,                                                     &
           jzoomin,                                                     &
           jzoomax,                                                     &
           iairsea,                                                     &
           ialbedo,                                                     &
           airseaoption,                                                &
           iobc_f,                                                      &
           iobc_wv,                                                     &
           iobc_ogcm,                                                     &
           iobc_lr,                                                     &
           obc_option,                                                  &
           iobc_demo_wv,                                                &
           iarchive,                                                    &
           imodeltrc,                                                   &
           imodelbio,                                                   &
           multiple,                                                    &
           kvarmax,                                                     &
           igesig,                                                      &
           isigfile,                                                    &
           ksecu,                                                       &
           run_option                !norestart
      integer ::                                                       &
           kmaxtide           &
          ,kmaxtidep1         &
          ,kminserie          &
          ,nzctdmax           &
          ,nsctdmax           &
          ,ktctdmin           &
          ,ktctdmax           &
          ,kdtk_out           &
          ,k10                &
          ,k11                &
          ,nbinco             &
          ,nbequa             &
          ,nbsparse           &
          ,kland              &
          ,tidestep2          &
          ,tidestep3          &
          ,ncfilm_max         &
          ,in_out_tide        &
          ,kmax_dof           &
          ,k_in               &
          ,k_out              &
          ,used_unused_dom    &
          ,mergebathy_sponge  &
          ,id_tem             &
          ,id_tem2            &
          ,id_dtem            &
          ,id_sal             &
          ,id_sal2            &
          ,id_rhp             &
          ,id_rhop            &
          ,id_rhom            &
          ,id_rhf             &
          ,id_rhc             &
          ,id_rhpa            &
          ,id_rhpb            &
          ,id_z               &
          ,id_prs             &
          ,id_now             &
          ,id_aft             &
          ,id_eost            &
          ,id_eoss            &
          ,id_u_now           & ! noretstart
          ,id_v_now           & ! noretstart
          ,id_ffreq      =0   & ! noretstart
          ,id_coriolis1       & ! noretstart
          ,id_coriolis2       & ! noretstart
          ,id_tcn        =0   & ! noretstart
          ,id_scn        =7   & ! noretstart
          ,id_gradssh    =1   & ! noretstart
          ,iDtOvRhCp     =0

      integer ::                                                        &
           idate_output,                                                &
           ihybsig,                                                     &
           nhybsig,                                                     &
           tke_surf,                                                    &
           grh_out_mi,                                                  &
           nest_onoff_in,                                               &
           nest_onoff_demo,                                             &
           nest_full_in,                                                &
           nest_full_out,                                               &
           nest_onoff_out,                                              &
           ncmin_airsea,                                                &
           ncmin_river,                                                 &
           ogcm_rec_max,                                                  & ! norestart
           eos_author,                                                & ! norestart
           eos_comprs,                                                & ! norestart
           eos_linear,                                                & ! norestart
           obcfreeorfix,                                              &
           iwve,                                                      &
           wave_obc_type,                                             & !norestart
           dataperwavefile,                                           &
           sp_or_db ,                                                 &
           kbu,                                                         &
           kbumax,                                                      &
           n_element,                                                   &
           relaxtype_ts,                                                &
           obctype_ts,                                                  &
           obctype_p                                                    &
          ,restart_file_y_or_n                                          &
          ,ncid                                                         &
          ,x_xhl_dim                                                    &
          ,y_xhl_dim                                                    &
          ,z_xhl_dim                                                    &
          ,x_yhl_dim                                                    &
          ,y_yhl_dim                                                    &
          ,z_yhl_dim                                                    &
          ,x_zhl_dim                                                    &
          ,y_zhl_dim                                                    &
          ,z_zhl_dim                                                    &
          ,x_zl_dim                                                     &
          ,y_zl_dim                                                     &
          ,z_zl_dim                                                     &
          ,time_dim                                                     &
          ,i_t_dim                                                      &
          ,j_t_dim                                                      &
          ,k_t_dim                                                      &
          ,i_w_dim                                                      &
          ,j_w_dim                                                      &
          ,k_w_dim                                                      &
          ,i_u_dim                                                      &
          ,j_u_dim                                                      &
          ,k_u_dim                                                      &
          ,i_v_dim                                                      &
          ,j_v_dim                                                      &
          ,k_v_dim                                                      &
          ,i_f_dim                                                      &
          ,j_f_dim                                                      &
          ,k_f_dim                                                      &
          ,dayindex_size                                                &
          ,tide_year_min                                                &
          ,airsea_year_min                                              &
          ,tide_year_max                                                &
          ,airsea_year_max                                              &
          ,wave_year_min                                                &
          ,irelaxsst                                                    &
          ,removetide                             & ! norestart
          ,ofl_rotation=0                         & ! norestart
          ,ofl_rhp=0                              & ! norestart
          ,ofl_tke=0                              & ! norestart
          ,ofl_surflux                            & ! norestart
          ,ale_selected                                                 &
          ,flag_asselin                                                 &
!         ,flag_lptimefilter                                            &
          ,freq
      integer ::                                                     &
           year_now                                                  &
          ,month_now                                                 &
          ,day_now                                                   &
          ,hour_now                                                  &
          ,minute_now                                                &
          ,second_now                                                &
          ,nd_send_est                                               &
          ,nd_send_ouest                                             &
          ,nd_send_nord                                              &
          ,nd_send_sud                                               &
          ,nd_send_out                                               &
          ,nd_recv_est                                               &
          ,nd_recv_ouest                                             &
          ,nd_recv_nord                                              &
          ,nd_recv_sud                                               & 
          ,nd_send_sudouest                                          &
          ,nd_send_sudest                                            &
          ,nd_send_nordouest                                         &
          ,nd_send_nordest                                           &
          ,nd_recv_sudouest                                          &
          ,nd_recv_sudest                                            &
          ,nd_recv_nordouest                                         &
          ,nd_recv_nordest                                           &
          ,initial_main_status                                       &
          ,offline_init_status                                       &
          ,meteo_sealand_mask                                        &
          ,grid_i0                                                   &
          ,grid_j0                                                   &
          ,ifb                                                       &
          ,dim_airsea                                                &
          ,rhp_zavr_xy  =0     ! norestart

      integer ::          &
           ssr_id    =0   & ! norestart
          ,ir_id     =0   & ! norestart
          ,rain_id   =0   & ! norestart
          ,t2m_id    =0   & ! norestart
          ,t0m_id    =0   & ! norestart
          ,abl_id    =0   & ! norestart
          ,dp2m_id   =0   & ! norestart
          ,u10m_id   =0   & ! norestart
          ,v10m_id   =0   & ! norestart
          ,u100m_id  =0   & ! norestart
          ,v100m_id  =0   & ! norestart
          ,p0m_id    =0   & ! norestart
          ,ustrs_id  =0   & ! norestart
          ,vstrs_id  =0   & ! norestart
          ,slhf_id   =0   & ! norestart
          ,netir_id  =0   & ! norestart
          ,sshf_id   =0     ! norestart

      double precision ::                                            &
        som0                                                         &
       ,som2                                                         &
       ,b2delev                                                      &
       ,b2delev_glb                                                  &
       ,sum0                                                         &
       ,sum1                                                         &
       ,sum2                                                         &
       ,sum3                                                         &
       ,sum4                                                         &
       ,sum5                                                         &
       ,sum6                                                         &
       ,sum7                                                         &
       ,sum8                                                         &
       ,sum9                                                         &
       ,sum10                                                        &
       ,sum11                                                        &
       ,sum12                                                        &
       ,sum13                                                        &
       ,time0                                                        &
       ,time1                                                        &
       ,time2                                                        &
       ,small2                                                       &
       ,x1_r8                                                        &
       ,x2_r8                                                        &
       ,x3_r8                                                        &
       ,x4_r8                                                        &
       ,sum0glb                                                      &
       ,sum1glb                                                      &
       ,sum2glb                                                      &
       ,sum3glb                                                      &
       ,sum4glb                                                      &
       ,sum5glb                                                      &
       ,sum6glb                                                      &
       ,emin                                                         &
       ,epsmin                                                       &
       ,elapsedtime_out                                              &
       ,elapsedtime_rst                                              &
       ,elapsedtime_end                                        &  ! norestart
       ,elapsedtime_now                                        &
       ,elapsedtime_bef                                        &
       ,elapsedtime_aft                                        &
       ,cpu_seconds                                            &
       ,alpha                                                  &
       ,eos_pgfzref                                            & ! norestart
       ,eos_tkezref                                              ! norestart

! fin module ne pas toucher à cette ligne
!*************************************************************

contains

      subroutine principal_allocate
      implicit none
#ifdef parallele
      include 'mpif.h'
#endif
      integer :: rank, ierr
#ifdef synopsis
       call main_synopsis('principal_allocate' &
                         ,'allocate principal arrays')
#endif

      m_mi=  imax*mi_onoff
      n_mi=  jmax*mi_onoff
      nr_mi=(kmax-1)*mi_onoff+1

      nest_m=nest_dim0*(imax-1)+1
      nest_n=nest_dim0*(jmax-1)+1
      nest_r=nest_dim0*(kmax-1)+1

      imax_w=onoff_wave*imax+(1-onoff_wave)*1
      jmax_w=onoff_wave*jmax+(1-onoff_wave)*1

!! Extra area for multi-grid approaches
      deltamax=16
      if(deltamax<1)stop 'deltamax<1 in module_principal'

      ifb=0
      dim_airsea=8

!........................................................

      allocate(tab1_code                   (0:nbdom-1))  ; tab1_code=0
      allocate(tab2_code                   (0:nbdom-1))  ; tab2_code=0
      allocate(tab_code_glb                (0:nbdom-1))  ; tab_code_glb=0

      allocate(ub2                         (2)) ; ub2=0
      allocate(lb2                         (2)) ; lb2=0
      allocate(ub3                         (3)) ; ub3=0
      allocate(lb3                         (3)) ; lb3=0
      allocate(ub4                         (4)) ; ub4=0
      allocate(lb4                         (4)) ; lb4=0
      allocate(obcstatus                  (10)) ; obcstatus=0

      allocate(anyv3d          (-1:imax+2,-1:jmax+2,0:kmax+1,any1:any2)) ; anyv3d=0
      allocate(checkanyv3d     (                             any1:any2)) ; checkanyv3d=-9999.
! checkanyv3d sert A verifier que le tableau anyv3d est disponible

      allocate(dz_u                (1:imax+1,0:jmax+1,0:kmax+1,-1:2)   ) ; dz_u=0
      allocate(dz_v                (0:imax+1,1:jmax+1,0:kmax+1,-1:2)   ) ; dz_v=0
!     allocate(dz_u                (0:imax+1,0:jmax+1,0:kmax+1, 0:2)   ) ; dz_u=0
!     allocate(dz_v                (0:imax+1,0:jmax+1,0:kmax+1, 0:2)   ) ; dz_v=0
      allocate(dz_t                (0:imax+1,0:jmax+1,0:kmax+1,-1:2)   ) ; dz_t=0
      allocate(stokesforces_u      (2:imax,2:jmax-1,kmax)              ) ; stokesforces_u=0
      allocate(stokesforces_v      (2:imax-1,2:jmax,kmax)              ) ; stokesforces_v=0
      allocate(depth_t             (-1:imax+2,-1:jmax+2,0:kmax+1)      ) ; depth_t=0
!     allocate(sigma_w             (0:imax+1,0:jmax+1,0:kmax+1)        ) ; sigma_w=0
      allocate(dsig_t              (0:imax+1,0:jmax+1,0:kmax+1)        ) ; dsig_t=0
      allocate(dsig_u              (1:imax+1,0:jmax+1,1:kmax  )        ) ; dsig_u=0
      allocate(dsig_v              (0:imax+1,1:jmax+1,1:kmax  )        ) ; dsig_v=0
      allocate(depth_u             (1:imax+1,0:jmax+1,0:kmax+1)        ) ; depth_u=0
      allocate(depth_v             (0:imax+1,1:jmax+1,0:kmax+1)        ) ; depth_v=0
      allocate(depth_w             (0:imax+1,0:jmax+1,0:kmax+1)        ) ; depth_w=0
      allocate(km_w                (1:imax  ,1:jmax  ,1:kmax+1)        ) ; km_w=0
      allocate(kh_w                (1:imax  ,1:jmax  ,1:kmax+1)        ) ; kh_w=0
      allocate(presgrad_u          (2:imax  ,2:jmax-1,1:kmax  ,0:1)    ) ; presgrad_u=0
      allocate(rhpref_u            (1:imax+1,1:jmax  ,0:kmax+1)        ) ; rhpref_u=0
      allocate(temref_u            (1:imax+1,1:jmax  ,0:kmax+1)        ) ; temref_u=0
      allocate(salref_u            (1:imax+1,1:jmax  ,0:kmax+1)        ) ; salref_u=0
      allocate(rhpref_v            (1:imax  ,1:jmax+1,0:kmax+1)        ) ; rhpref_v=0
      allocate(temref_v            (1:imax  ,1:jmax+1,0:kmax+1)        ) ; temref_v=0
      allocate(salref_v            (1:imax  ,1:jmax+1,0:kmax+1)        ) ; salref_v=0
      allocate(presgrad_v          (2:imax-1,2:jmax  ,1:kmax  ,0:1)    ) ; presgrad_v=0
      allocate(rhp_t               (0:imax+1,0:jmax+1,0:kmax+1)        ) ; rhp_t=0
      allocate(rhpzavr_w           (0:imax+1,0:jmax+1)                 ) ; rhpzavr_w=0.
      allocate(vel_u               (0:imax+2,0:jmax+1,-1:kmax+2,-1:2)   ) ; vel_u=0
      allocate(vel_v               (0:imax+1,0:jmax+2,-1:kmax+2,-1:2)   ) ; vel_v=0
      allocate(omega_w             (0:imax+1,0:jmax+1,0:kmax+1,1)      ) ; omega_w=0
      allocate(veldydz_u           (0:imax+2,0:jmax+1,0:kmax+1,1)      ) ; veldydz_u=0
      allocate(veldxdz_v           (0:imax+1,0:jmax+2,0:kmax+1,1)      ) ; veldxdz_v=0
      allocate(veldtodx_u          (0:imax+2,0:jmax+1,0:kmax+1,1)      ) ; veldtodx_u=0
      allocate(veldtody_v          (0:imax+1,0:jmax+2,0:kmax+1,1)      ) ; veldtody_v=0
      allocate(tem_t               (-1:imax+2,-1:jmax+2,0:kmax+1,-1:2) ) ; tem_t=0
      allocate(rhpref_t            (-1:imax+2,-1:jmax+2,0:kmax+1) )      ; rhpref_t=0.
      allocate(rhcref_t            ( 1:imax  , 1:jmax  ,1:kmax  ) )      ; rhcref_t=0.
      allocate(temref_t            (-1:imax+2,-1:jmax+2,0:kmax+1) )      ; temref_t=0.
      allocate(salref_t            (-1:imax+2,-1:jmax+2,0:kmax+1) )      ; salref_t=0.
      allocate(temref_z            ( 0:imax+1, 0:jmax+1,1:kmax  ) )      ; temref_z=0.
      allocate(salref_z            ( 0:imax+1, 0:jmax+1,1:kmax  ) )      ; salref_z=0.
      allocate(zref_z              (                    1:kmax  ) )      ; zref_z=0.
      allocate(sal_t               (-1:imax+2,-1:jmax+2,0:kmax+1,-1:2) ) ; sal_t=0
      allocate(tridia_out          (0:imax+1,0:jmax+1,0:kmax+1)        ) ; tridia_out=0
      allocate(tridia_in           (0:imax+1,0:jmax+1,kmax+1,-1:4)     ) ; tridia_in=0
      allocate(fluxbar_sumt_u      (1:imax+1,1:jmax  ,0:1)             ) ; fluxbar_sumt_u=0
      allocate(fluxbar_sumt_v      (1:imax  ,1:jmax+1,0:1)             ) ; fluxbar_sumt_v=0
      allocate(lon_t               (-1:imax+2,-1:jmax+2)               ) ; lon_t=0
      allocate(lat_t               (-1:imax+2,-1:jmax+2)               ) ; lat_t=0
      allocate(lon_u               (0:imax+2,0:jmax+2)                 ) ; lon_u=0
      allocate(lat_u               (0:imax+2,0:jmax+2)                 ) ; lat_u=0
      allocate(lon_v               (0:imax+2,0:jmax+2)                 ) ; lon_v=0
      allocate(lat_v               (0:imax+2,0:jmax+2)                 ) ; lat_v=0
      allocate(lon_f               (0:imax+2,0:jmax+2)                 ) ; lon_f=0
      allocate(lat_f               (0:imax+2,0:jmax+2)                 ) ; lat_f=0
      allocate(velbar_u            (0:imax+2,0:jmax+1, -1:2)           ) ; velbar_u=0
      allocate(velbar_v            (0:imax+1,0:jmax+2, -1:2)           ) ; velbar_v=0
!     allocate(noise_velbar_u      (1:imax+1,1:jmax)                   ) ; noise_velbar_u=0
!     allocate(noise_velbar_v      (1:imax  ,1:jmax+1)                 ) ; noise_velbar_v=0
      allocate(gradssh_u           (2:imax,2:jmax-1,0:1)) ; gradssh_u=0.
      allocate(pgf_u               (2:imax,2:jmax-1,0:1)) ; pgf_u=0.
      allocate(gradssh_v           (2:imax-1,2:jmax,0:1)) ; gradssh_v=0.
      allocate(pgf_v               (2:imax-1,2:jmax,0:1)) ; pgf_v=0.
      allocate(sponge_u            (1:imax+1,1:jmax  ,2)               ) ; sponge_u=0
      allocate(sponge_v            (1:imax  ,1:jmax+1,2)               ) ; sponge_v=0
      allocate(fluxbar_u           (1:imax+1,1:jmax  ,0:1)             ) ; fluxbar_u=0
      allocate(fluxbar_v           (1:imax  ,1:jmax+1,0:1)             ) ; fluxbar_v=0
      allocate(ssh_w               (0:imax+1,0:jmax+1,-1:2)            ) ; ssh_w=0
      allocate(dxdy_u              (1:imax+1,0:jmax+1)                 ) ; dxdy_u=0
      allocate(dxdy_v              (0:imax+1,1:jmax+1)                 ) ; dxdy_v=0
      allocate(dxdy_t              (0:imax+1,0:jmax+1)                 ) ; dxdy_t=0
      allocate(dx_u                (1:imax+1,0:jmax+1)                 ) ; dx_u=0
      allocate(dx_v                (0:imax+1,1:jmax+1)                 ) ; dx_v=0
      allocate(dx_t                (0:imax+1,0:jmax+1)                 ) ; dx_t=0
      allocate(dx_f                (0:imax+1,0:jmax+1)                 ) ; dx_f=0
      allocate(dy_u                (1:imax+1,0:jmax+1)                 ) ; dy_u=0
      allocate(dy_v                (0:imax+1,1:jmax+1)                 ) ; dy_v=0
      allocate(dy_t                (0:imax+1,0:jmax+1)                 ) ; dy_t=0
      allocate(dy_f                (0:imax+1,0:jmax+1)                 ) ; dy_f=0
      allocate(hz_u                (1:imax+1,0:jmax+1,0:2)             ) ; hz_u=0
      allocate(hz_v                (0:imax+1,1:jmax+1,0:2)             ) ; hz_v=0
      allocate(hz_w                (0:imax+1,0:jmax+1,0:2)             ) ; hz_w=0
      allocate(h_w                 (0:imax+1,0:jmax+1)                 ) ; h_w=0
      allocate(h_u                 (1:imax+1,0:jmax+1)                 ) ; h_u=0
      allocate(h_v                 (0:imax+1,1:jmax+1)                 ) ; h_v=0
      allocate(h_f                 (0:imax+2,0:jmax+2)                 ) ; h_f=0
      allocate(sponge_t            (1:imax  ,1:jmax  ,1)               ) ; sponge_t=0
      allocate(zeroleveldepth_u    (1:imax+1,1:jmax  ,1)               ) ; zeroleveldepth_u=0
      allocate(zeroleveldepth_v    (1:imax  ,1:jmax+1,1)               ) ; zeroleveldepth_v=0
      allocate(heatrelax_w         (0:imax+1,0:jmax+1,1)               ) ; heatrelax_w=0
      allocate(coriolis_t          (0:imax+1,0:jmax+1)                 ) ; coriolis_t=0
      allocate(wetmask_u           (0:imax+1,0:jmax+1)                 ) ; wetmask_u=0
      allocate(wetmask_v           (0:imax+1,0:jmax+1)                 ) ; wetmask_v=0
      allocate(wetmask_t           (0:imax+1,0:jmax+1)                 ) ; wetmask_t=0
      allocate(xy_u                (0:imax+2,0:jmax+1,0:10)            ) ; xy_u=0
      allocate(xy_v                (0:imax+1,0:jmax+2,0:10)            ) ; xy_v=0
      allocate(xy_t                (-1:imax+2,-1:jmax+2,0:8)           ) ; xy_t=0
      allocate(xy_f                (0:imax+2,0:jmax+2,5)               ) ; xy_f=0
      allocate(ssr_w               (0:imax+1,0:jmax+1,0:2)             ) ; ssr_w=0
      allocate(snsf_w              (0:imax+1,0:jmax+1,0:2)             ) ; snsf_w=0
      allocate(precipi_w           (0:imax+1,0:jmax+1,0:2)             ) ; precipi_w=0
      allocate(pss_w               (0:imax+1,0:jmax+1,0:2)             ) ; pss_w=0
      allocate(uwind_t             (0:imax+1,0:jmax+1,0:2)             ) ; uwind_t=0
      allocate(vwind_t             (0:imax+1,0:jmax+1,0:2)             ) ; vwind_t=0
      allocate(fric_u              (0:imax+1,0:jmax+1)                 ) ; fric_u=0
      allocate(fric_v              (0:imax+1,0:jmax+1)                 ) ; fric_v=0
      allocate(fric_t              (0:imax+1,0:jmax+1)                 ) ; fric_t=0
      allocate(velbot3d2d_u        (0:imax+1,0:jmax+1)                 ) ; velbot3d2d_u=0
      allocate(velbot3d2d_v        (0:imax+1,0:jmax+1)                 ) ; velbot3d2d_v=0
      allocate(cdb_t               (0:imax+1,0:jmax+1)                 ) ; cdb_t=0
      allocate(xflux_t             (0:imax+1,0:jmax+1)                 ) ; xflux_t=0
      allocate(xflux_f             (0:imax+1,0:jmax+1)                 ) ; xflux_f=0
      allocate(yflux_t             (0:imax+1,0:jmax+1)                 ) ; yflux_t=0
      allocate(yflux_f             (0:imax+1,0:jmax+1)                 ) ; yflux_f=0
      allocate(frozenterm2d_u      (2:imax  ,2:jmax-1,-1:1)            ) ; frozenterm2d_u=0.
      allocate(frozenterm2d_v      (2:imax-1,2:jmax  ,-1:1)            ) ; frozenterm2d_v=0.
      allocate(frozenterm3d_u      (2:imax  ,2:jmax-1,-1:1)            ) ; frozenterm3d_u=0.
      allocate(frozenterm3d_v      (2:imax-1,2:jmax  ,-1:1)            ) ; frozenterm3d_v=0.
      allocate(pres3d2d_u          (0:imax+1,0:jmax+1)                 ) ; pres3d2d_u=0
      allocate(pres3d2d_v          (0:imax+1,0:jmax+1)                 ) ; pres3d2d_v=0
      allocate(adve3d2d_u          (0:imax+1,0:jmax+1)                 ) ; adve3d2d_u=0
      allocate(adve3d2d_v          (0:imax+1,0:jmax+1)                 ) ; adve3d2d_v=0
!     allocate(restoring3d2d_u     (2:imax  ,2:jmax-1)                 ) ; restoring3d2d_u=0
!     allocate(restoring3d2d_v     (2:imax-1,2:jmax  )                 ) ; restoring3d2d_v=0
      allocate(stokesforces3d2d_u  (0:imax+1,0:jmax+1)                 ) ; stokesforces3d2d_u=0
      allocate(stokesforces3d2d_v  (0:imax+1,0:jmax+1)                 ) ; stokesforces3d2d_v=0
      allocate(wstress_w           (0:imax+1,0:jmax+1)                 ) ; wstress_w=0
      allocate(z0_w                (0:imax+1,0:jmax+1)                 ) ; z0_w=0
      allocate(albedo_w            (0:imax+1,0:jmax+1)                 ) ; albedo_w=0
      allocate(gridrotcos_t        (0:imax+1,0:jmax+1)                 ) ; gridrotcos_t=0
      allocate(gridrotsin_t        (0:imax+1,0:jmax+1)                 ) ; gridrotsin_t=0
      allocate(light_kpar2_w       (imax,jmax)                         ) ; light_kpar2_w=1.
      allocate(sshstokes_w         (0:imax+1,0:jmax+1)                 ) ; sshstokes_w=0
!     allocate(divflux0_u          (2:imax,2:jmax-1)) ; divflux0_u=0.
!     allocate(divflux0_v          (2:imax-1,2:jmax)) ; divflux0_v=0.
!     allocate(ssh_int_u           (0:imax+1,0:jmax+1,2:2)             ) ; ssh_int_u=0
!     allocate(ssh_int_v           (0:imax+1,0:jmax+1,2:2)             ) ; ssh_int_v=0
      allocate(ssh_int_w           (0:imax+1,0:jmax+1,0:2)             ) ; ssh_int_w=0
!     allocate(ssh_avr_w           (0:imax+1,0:jmax+1,0:1)             ) ; ssh_avr_w=0
!     allocate(streamf_f           (0:imax+1,0:jmax+1,0:2)             ) ; streamf_f=0
      allocate(ustokesvortex_t     (0:imax+1,0:jmax+1,2)               ) ; ustokesvortex_t=0
      allocate(ustokesvortex_f     (0:imax+1,0:jmax+1,2)               ) ; ustokesvortex_f=0
      allocate(vstokesvortex_t     (0:imax+1,0:jmax+1,2)               ) ; vstokesvortex_t=0
      allocate(vstokesvortex_f     (0:imax+1,0:jmax+1,2)               ) ; vstokesvortex_f=0
      allocate(velavr_u            (0:imax+1,0:jmax+1,0:1)             ) ; velavr_u=0
      allocate(velavr_v            (0:imax+1,0:jmax+1,0:1)             ) ; velavr_v=0
      allocate(timefilter_u        (0:imax+1,0:jmax+1,0:1)             ) ; timefilter_u=0
      allocate(timefilter_v        (0:imax+1,0:jmax+1,0:1)             ) ; timefilter_v=0
!     allocate(grdobc3d_j          (jmax,kmax,2)                       ) ; grdobc3d_j=0
!     allocate(grdobc3d_i          (imax,kmax,2)                       ) ; grdobc3d_i=0
      allocate(cwi_int_u           (0:imax+1,2)                        ) ; cwi_int_u=0
      allocate(cwi_int_v           (0:imax+1,2)                        ) ; cwi_int_v=0
      allocate(cwj_int_u           (0:jmax+1,2)                        ) ; cwj_int_u=0
      allocate(cwj_int_v           (0:jmax+1,2)                        ) ; cwj_int_v=0
      allocate(sqr_hoverg_u             (2:jmax-1,2)                   ) ; sqr_hoverg_u=0.
      allocate(sqr_hoverg_v             (2:imax-1,2)                   ) ; sqr_hoverg_v=0.
      allocate(sshrefobc_i         (imax,2)                            ) ; sshrefobc_i=0
      allocate(vbrrefobc_i         (imax,2)                            ) ; vbrrefobc_i=0
      allocate(sshrefobc_j         (jmax,2)                            ) ; sshrefobc_j=0
      allocate(vbrrefobc_j         (jmax,2)                            ) ; vbrrefobc_j=0
      allocate(anyv1d              (kmax+1,2)                          ) ; anyv1d=0
      allocate(hr_z2lr_w           (0:nest_m+1,0:nest_n+1,nest_r,3)    ) ; hr_z2lr_w=0
      allocate(vmi_u               (0:m_mi+1,0:n_mi+1,nr_mi)           ) ; vmi_u=0
      allocate(vmi_v               (0:m_mi+1,0:n_mi+1,nr_mi)           ) ; vmi_v=0
      allocate(gridcard            (11,2)                              ) ; gridcard=0
      allocate(novector            (4,4)                               ) ; novector=0
      allocate(proficard           (0:nbincomax,3)                     ) ; proficard=0
      allocate(airseainfo          (dim_airsea,3)                      ) ; airseainfo=0
      allocate(airseadt            (dim_airsea,2)                      ) ; airseadt=0
      allocate(riverdt             (dim_river,2)                       ) ; riverdt=0
      allocate(river_s             (dim_river)                         ) ; river_s=0
      allocate(river_tmin          (dim_river)                         ) ; river_tmin=0
      allocate(river_tmax          (dim_river)                         ) ; river_tmax=0
      allocate(hriver              (dim_river)                         ) ; hriver=0
      allocate(friver              (dim_river)                         ) ; friver=0
      allocate(realriver           (dim_river)                         ) ; realriver=0
      allocate(riverinfo           (dim_river)                         ) ; riverinfo=0
      allocate(river_timeref       (dim_river)                         ) ; river_timeref=0
      allocate(river_t             (dim_river,0:2)                     ) ; river_t=0
      allocate(riverflux           (dim_river,0:2)                     ) ; riverflux=0
      allocate(daysim              (0:10+dim_airsea+3+dim_river)       ) ; daysim=0
      allocate(tdate_output        (dim_dof)                           ) ; tdate_output=0
      allocate(obcinfo             (3)                                 ) ; obcinfo=0
      allocate(offlinedt           (2)                                 ) ; offlinedt=0
      allocate(pss_mean            (0:2)                               ) ; pss_mean=0
      allocate(sshobc_cum          (0:2)                               ) ; sshobc_cum=0
      allocate(nest_dt_in          (2)                                 ) ; nest_dt_in=0
      allocate(nest_dt_out         (2)                                 ) ; nest_dt_out=0
!     allocate(obcafdt             (2)                                 ) ; obcafdt=0
      allocate(wavedt              (2)                                 ) ; wavedt=0
      allocate(obc_hf_dt           (2)                                 ) ; obc_hf_dt=0
      allocate(mpi_neighbor_list   (10)) ; mpi_neighbor_list=-999

!     allocate(mask_t              (-1:imax+2,-1:jmax+2,0:kmax+1)) ; mask_t=0
      i1=0-deltamax ; i2=imax+1+deltamax ; j1=0-deltamax ; j2=jmax+1+deltamax
      allocate(mask_t              (i1:i2,j1:j2,0:kmax+1)) ; mask_t=0

      allocate(mask_f              (0:imax+1,0:jmax+1,0:kmax+1))   ; mask_f=0
      allocate(mask_u              (0:imax+2,0:jmax+1,0:kmax+1))   ; mask_u=0
      allocate(mask_v              (0:imax+1,0:jmax+2,0:kmax+1))   ; mask_v=0

      allocate(ksl_t               (imax,jmax))         ; ksl_t=kmax
      allocate(upwindriver_t       (0:imax+1,0:jmax+1)) ; upwindriver_t=1.

!     allocate(anyvarlr            (0:mlr+1,0:nlr+1,nrlr-1,2)) ; anyvarlr=0

      allocate(temobc_t            (0:imax+1,0:jmax+1,kmax+1,0:2)) ; temobc_t=0
      allocate(salobc_t            (0:imax+1,0:jmax+1,kmax+1,0:2)) ; salobc_t=0
      allocate(velobc_u            (0:imax+1,0:jmax+1,kmax+1,0:2)) ; velobc_u=0
      allocate(velobc_v            (0:imax+1,0:jmax+1,kmax+1,0:2)) ; velobc_v=0

      allocate(velstokes_u         (0:imax+1,0:jmax+1,0:kmax+1,1:2)) ; velstokes_u=0
      allocate(velstokes_v         (0:imax+1,0:jmax+1,0:kmax+1,1:2)) ; velstokes_v=0
      allocate(velbarstokes_u      (0:imax+1,0:jmax+1,1:2))          ; velbarstokes_u=0
      allocate(velbarstokes_v      (0:imax+1,0:jmax+1,1:2))          ; velbarstokes_v=0

      allocate(anyv3dint           (-1:imax+2,-1:jmax+2,kmax+1)) ; anyv3dint=0
      allocate(anyvar3d            (-1:imax+2,-1:jmax+2,kmax+1)) ; anyvar3d=0
      allocate(anyvar2d            (-1:imax+2,-1:jmax+2))        ; anyvar2d=0

      allocate(sshobc_w            (0:imax+1,0:jmax+1,0:2)) ; sshobc_w=0
      allocate(velbarobc_u         (0:imax+1,0:jmax+1,0:2)) ; velbarobc_u=0
      allocate(velbarobc_v         (0:imax+1,0:jmax+1,0:2)) ; velbarobc_v=0

      allocate(riverdir            (dim_river))   ; riverdir=0
      allocate(l_river             (dim_river))   ; l_river=0
      allocate(river_no            (dim_river))   ; river_no=0
      allocate(rivertrc_inout      (dim_river))   ; rivertrc_inout=0
      allocate(rivervel_inout      (dim_river))   ; rivervel_inout=0
      allocate(riverupwdist        (dim_river))   ; riverupwdist=0
      allocate(iriver              (dim_river,3)) ; iriver=0
      allocate(jriver              (dim_river,3)) ; jriver=0

      allocate(jpm                 (12))  ; jpm=0
      allocate(nest_len_in         (0:1)) ; nest_len_in=0

      allocate(spo_i1_u            (5)) ; spo_i1_u=0
      allocate(spo_i2_u            (5)) ; spo_i2_u=0
      allocate(spo_j1_u            (5)) ; spo_j1_u=0
      allocate(spo_j2_u            (5)) ; spo_j2_u=0
      allocate(spo_i1_v            (5)) ; spo_i1_v=0
      allocate(spo_i2_v            (5)) ; spo_i2_v=0
      allocate(spo_j1_v            (5)) ; spo_j1_v=0
      allocate(spo_j2_v            (5)) ; spo_j2_v=0
      allocate(spo_i1_t            (5)) ; spo_i1_t=0
      allocate(spo_i2_t            (5)) ; spo_i2_t=0
      allocate(spo_j1_t            (5)) ; spo_j1_t=0
      allocate(spo_j2_t            (5)) ; spo_j2_t=0
      allocate(obcstreamf          (4)) ; obcstreamf=0

      allocate(mask_i_w            (-1:imax+2))   ; mask_i_w=0
      allocate(mask_i_v            (imax))   ; mask_i_v=0
      allocate(mask_i_u            (imax+1)) ; mask_i_u=0

      allocate(mask_j_w            (-1:jmax+2))   ; mask_j_w=0
      allocate(mask_j_u            (jmax))   ; mask_j_u=0
      allocate(mask_j_v            (jmax+1)) ; mask_j_v=0

      allocate(kmin_u              ( 0:imax+1, 0:jmax+1)) ; kmin_u=1
      allocate(kmin_v              ( 0:imax+1, 0:jmax+1)) ; kmin_v=1
      allocate(kmin_w              (-1:imax+2,-1:jmax+2)) ; kmin_w=1

      allocate(l2ij_out_u          (0:nest_dim1,2)) ; l2ij_out_u=0
      allocate(l2ij_out_v          (0:nest_dim2,2)) ; l2ij_out_v=0
      allocate(l2ij_out_w          (0:nest_dim3,2)) ; l2ij_out_w=0
      allocate(l2ij_in_u           (0:nest_dim4,2)) ; l2ij_in_u=0
      allocate(l2ij_in_v           (0:nest_dim5,2)) ; l2ij_in_v=0
      allocate(l2ij_in_w           (0:nest_dim6,2)) ; l2ij_in_w=0
      allocate(fillmask_t          (0:nest_dim7,4)) ; fillmask_t=0

      allocate(dayindex            (366*numbyears)) ; dayindex=0
      allocate(grh_out_var         (kmax_grh))      ; grh_out_var=0
!     allocate(varid               (kmax_grh))      ; varid=0
      allocate(grh_nb              (6))             ; grh_nb=0

      allocate(departuredate       (6))            ; departuredate=-999
      allocate(datesim             (6,0:10))       ; datesim=0
      allocate(dateobc             (6,3))          ; dateobc=0
      allocate(dateairsea          (6,dim_airsea)) ; dateairsea=0
      allocate(dateriver           (6,dim_river))  ; dateriver=0

      allocate(vardim              (4)) ; vardim=0
      allocate(varstart            (4)) ; varstart=0
      allocate(varcount            (4)) ; varcount=0

      allocate(cgridshift          (3)) ; cgridshift=0

      allocate(nest_len_out        (0:nest_max)) ; nest_len_out=0
      allocate(ind                 (4))          ; ind=0

      allocate(texte80             (20))           ; texte80='s' 
!     allocate(nomfichier          (40))           ; nomfichier='s'
      allocate(grh_titre           (dimgrh_titre)) ; grh_titre='s'
      allocate(rivername           (dim_river))    ; rivername='s'
      allocate(airseafile          (dim_airsea))   ; airseafile='s'
      allocate(riverfile           (dim_river))    ; riverfile='s'
      allocate(obcfile             (9,2))          ; obcfile='s'
      allocate(month               (12))           ; month='s'
      allocate(offlinefile         (7))            ; offlinefile='s'
      allocate(nest_path_in        (0:1))          ; nest_path_in='s'
      allocate(nest_path_out       (0:nest_max))   ; nest_path_out='s'
      allocate(obc_hf_file         (3))            ; obc_hf_file='s'
      allocate(bioriv_file         (dim_river))    ; bioriv_file='s'

      allocate(tidepotential_w     (0:imax+1,0:jmax+1,1)) ; tidepotential_w=0
      allocate( sshtotide_w        (0:imax+1,0:jmax+1,1)) ; sshtotide_w=0
      allocate(sshmeatide          (0:2))                 ; sshmeatide=0
      allocate(timeweightobc       (var_num))             ; timeweightobc=0
      allocate(dtobc               (var_num))             ; dtobc=1.
      allocate(runoff_dt           (2))                   ; runoff_dt=0.
      allocate(runoff_w            (dim_river,0:2))       ; runoff_w=0.


      end subroutine principal_allocate

!.....................................................................

      subroutine principal_allocate_tide
      implicit none
#ifdef synopsis
       call main_synopsis('principal_allocate_tide' &
                         ,'allocate tide arrays')
#endif

      allocate(analysetide_u     (0:imax+1,0:jmax+1,kmaxtide*2)) ; analysetide_u=0
      allocate(analysetide_v     (0:imax+1,0:jmax+1,kmaxtide*2)) ; analysetide_v=0
      allocate(analysetide_w     (0:imax+1,0:jmax+1,kmaxtide*2)) ; analysetide_w=0
      allocate(tideanalysismatrix(kmaxtide*2,kmaxtide*2)       ) ; tideanalysismatrix=0
      allocate(sshtidecos_w      (0:imax+1,0:jmax+1,kmaxtide)  ) ; sshtidecos_w=0
      allocate(sshtidesin_w      (0:imax+1,0:jmax+1,kmaxtide)  ) ; sshtidesin_w=0
      allocate(veltidecos_u      (0:imax+1,0:jmax+1,kmaxtide)  ) ; veltidecos_u=0
      allocate(veltidesin_u      (0:imax+1,0:jmax+1,kmaxtide)  ) ; veltidesin_u=0
      allocate(veltidecos_v      (0:imax+1,0:jmax+1,kmaxtide)  ) ; veltidecos_v=0
      allocate(veltidesin_v      (0:imax+1,0:jmax+1,kmaxtide)  ) ; veltidesin_v=0
      allocate(potidecos_w       (0:imax+1,0:jmax+1,kmaxtide)  ) ; potidecos_w=0
      allocate(potidesin_w       (0:imax+1,0:jmax+1,kmaxtide)  ) ; potidesin_w=0
      allocate(veltidecosout_u   (0:imax+1,0:jmax+1,kmaxtide)  ) ; veltidecosout_u=0
      allocate(veltidesinout_u   (0:imax+1,0:jmax+1,kmaxtide)  ) ; veltidesinout_u=0
      allocate(veltidecosout_v   (0:imax+1,0:jmax+1,kmaxtide)  ) ; veltidecosout_v=0
      allocate(veltidesinout_v   (0:imax+1,0:jmax+1,kmaxtide)  ) ; veltidesinout_v=0
      allocate(sshtidecosout_w   (0:imax+1,0:jmax+1,kmaxtide)  ) ; sshtidecosout_w=0
      allocate(sshtidesinout_w   (0:imax+1,0:jmax+1,kmaxtide)  ) ; sshtidesinout_w=0
      allocate(passetide         (kmaxtide,2)                  ) ; passetide=0
      allocate(z1meatide         (kmaxtide,2)                  ) ; z1meatide=0
      allocate(z2meatide         (kmaxtide,2)                  ) ; z2meatide=0
      allocate(frqtide           (kmaxtide+1)                  ) ; frqtide=0
      allocate(v0tide            (kmaxtide+1)                  ) ; v0tide=0
      allocate(ti0tide           (kmaxtide+1)                  ) ; ti0tide=0
      allocate(equitide          (kmaxtide)                    ) ; equitide=0
      allocate(tideanaweight     (kmaxtide)                    ) ; tideanaweight=0
      allocate(ftide             (kmaxtide+1,0:2)              ) ; ftide=1.
      allocate(utide             (kmaxtide+1,0:2)              ) ; utide=0
      allocate(nutide            (kmaxtide)                    ) ; nutide=0
      allocate(nametide          (kmaxtide,4)                  ) ; nametide='s' 
      allocate(tidenodalfile     (kmaxtide)                    ) ; tidenodalfile='s'

      end subroutine principal_allocate_tide

!...........................................................................

      end module module_principal
