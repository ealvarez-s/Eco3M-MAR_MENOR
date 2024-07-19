      subroutine set_parameters
!______________________________________________________________________
! SYMPHONIE ocean model
! release 367 last update: 07-03-23
!______________________________________________________________________

      use module_principal  ; use module_parallele ; use module_wave
      use module_airseaflux ; use module_optics    ; use module_s
      implicit none
      integer unit_

! notebook_spongelayer !20-05-18
      namelist/notebook_spongelayer/sponge_l,relax_ext,relax_int,relaxtype_ts     &
              ,relax_ts,relax_es,obcfreeorfix,relax_bpc,obctype_ts,obctype_p      &
              ,flag_upwind_obc,flag_sponge_txt,sponge_dx_critic,sponge_dx_width   &
              ,flagspo_i1,flagspo_i2,flagspo_j1,flagspo_j2 & !17-01-19
              ,relax_lwf !22-11-20

! notebook_eqstate
      namelist /notebook_eqstate/eos_author,eos_comprs,eos_pgfzref     &
               ,eos_linear,t0,s0,alp_t,alp_s,rho,eos_tkezref,pgfscheme &
               ,rhp_zavr_xy,flag_refstate,flag_steric_effect           &
               ,tem_validmin,tem_validmax,sal_validmin,sal_validmax !04-03-19

! notebook_offline:
      namelist /notebook_offline/ioffline,removetide,directory_offline &
      ,offlinefile,ofl_rotation,ofl_rhp,ofl_surflux,flag_kz_enhanced   &
      ,ofl_type,ofl_sshtype,ofl_tke,flag_maxbotstress,flag_offline_binary & !08-06-18 
      ,flag_ksloffline,ofl_reversedtime,ofl_bio !26-05-19!09-09-19

! notebook_vertcoord
      namelist/notebook_vertcoord/igesig,isigfile,hgesig,pgesig      &
      ,ihybsig,nhybsig,nbvstepmin,hstepmin,hstepmax,verticalgridfile &
      ,ale_selected,sigstepgridfile,fgrid_or_wgrid,vststep,flag_merged_levels & !03-05-17
      ,dzsurfmin,flag_z2dv_outputs,vqs_file  &     !28-10-18!17-04-20
      ,dz_vertical_incr_fact                 &     !07-05-20
      ,vqs_cst1,vqs_cst2,vqs_cst3                  !22-03-21

! notebook_wave (Part I)
      namelist/notebook_wave1/iwve

! notebook_obcforcing
      namelist/notebook_obcforcing/iobc_ogcm,obc_option,obc_ogcm_type & !14-10-15
              ,obcfile,bi_onoff,ogcm_time_shift                       &
              ,offset_sshobc,flag_ogcmtidemixing & !19-09-16 !11-12-19
              ,flag_ogcm_instab                  & !19-02-22
              ,ogcm_time_lag                     & !17-03-22
              ,obctime_order                       !29-10-22

#ifdef synopsis
       subroutinetitle='set_parameters'
       subroutinedescription= &
       'Notebooks reading and initial value of some constants & arrays'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!______________________________________________________________
! Version date      Description des modifications:
!         13/08/01: plus d'infos dans le fichier "messages"
!                   ajout de la lecture de la position des
!                   4 coins de la zone d'archivage pour la
!                   procedure d'imbrication des grille.
!         20/08/01: ajout de la lecture de la carte
!                   d'identité de la grille "basse résolution".
!         21/08/01: bouton d'option "basse résolution" IOBC_LR
!                   + deplacement de l'initialisation des constantes
!                   PI, G, etc, en début de programme.
!         24/08/01: l'initialisation de ICHANEL9 est maintenant faite
!                   dans time_step.F
!         13/09/01: bienvenue à SMALL1
!         25/10/01: si iobc_lr=0 on se fout que MLR et NLR soient mal definis
!         12/11/01: ajout de relaxmin_ext et de relaxmin_int
!         06/12/01: lecture des coef des fonctions de stabilité de Mellor
!         01/01/02: lecture du temps initial des parametres maree
!         15/01/02: bienvenue à IN_OUT_TIDE
!         18/01/02: bienvenue à NAMETIDE
!         24/01/02: bienvenue à EQUITIDE
!         26/02/02: retrait de relaxmin_ext et de relaxmin_int
!         28/02/02: lecture notebook_obcforcing securisée
!         05/03/02: bienvenue à small2
!         17/04/02: ajout securité sur taille KMAXTIDE et NTIDE
!         20/04/02: notebook_obcforcing compte 2 lignes de plus
!         01/07/02: bienvenue à C1STREAMF & C2STREAMF
!         23/08/02: mise à jour lecture de notebook_airseaflux
!                   pour nouvelle filiere "formule bulk"
!                   et lecture de notebook_streamf (bienvenue
!                   aux conditions radiatives pour streamf_matrix.F)
!         02/10/02: retour à T0=20.
!                   bienvenue à T0_BASE S0_BASE RHO_BASE ALP_T_BASE ALP_S_BASE
!         24/03/03: ajout lecture de notebook_advection et notebook_vertcoord
!                   C1STREAMF et C2STREAMF supprimés
!         25/03/03: Include de biology2004.h
!         05/06/03: correction bug. On lit NRLR et pas NR
!         09/07/03: bienvenue à TKE_SURF (choix obc PRODI ou CB)
!         10/07/03: dans notebook_obcforcing lire nom fichier record.dat,
!                   reset real de KOUNTGRAPH
!         24/07/03: lecture du nom de fichier des courant de marée dans
!                   notebook_tide
!         02/08/03: bienvenue NHYBSIG (iteration max pour calcul grille hybride)
!         21/08/03: ajout d'une possibilité de selectionner les forces de
!                   marée: ondes + tout potentiel
!                           ondes + potentiel astro
!                           ondes + potentiel charge
!                           ondes + aucun potentiel
!                   c.a.d. bienvenue à TIDEFORCES
!                   supression TIDESTEP1
!         26/08/03: les parametres de la couche eponge sont maintenant lus
!                   à part dans un notebook_spongelayer
!         09/01/04: cas IAIRSEA=3
!         25/05/04: la nouvelle filiere d'imbrication amene à reduire la
!                   longueur de notebook_obcforcing
!         28/07/04: Mise à jour lecture de notebook_tide
!         02/12/04: IWIND=0 si AIRSEAFLUX.NE.0
!         06/01/05: coed d'asselin réduit à 0.1
!         20/01/05: un "close(3)" oublié....
!         04/02/05: des messages à l'ecran qui indiquent que la lecture des
!                   notebook se passe bien.
!         12/04/05: valeurs à priori des paramètres de l'équation d'état plus
!                   proches des valeurs méditerranéennes
!         14/06/05: suppression des accents sur le message d'erreur
!         19/06/05: lecture des coef de moderation du limiteur de flux du
!                   schema d'advection CST_ADV_HOR & CST_ADV_VER
!         03/08/05: lecture de OBCFREEORFIX ajouté dans notebook_spongelayer
!         30/09/05: bienvenue à UNIT_R4 pour compiler avec -r8
!         28/12/05: bienvenue à IALBEDO
!         18/01/06: ABSORPTION plusieurs cas de Jerlov proposés (CAS I, II, III)
!         16/02/06: bienvenue à RELAX_TS
!                   et à notebook_optical dans lequel on lit les parametres
!                   d'attenuation de la lumiere dans l'eau.
!         07/04/06: SP_OR_DB permet de distinguer un fichier binaire simple
!                   precision (=4) d'un fichier double precision (=8)
!         04/05/06: bienvenu a ZERO pour pouvoir compiler en double precision
!                   avec f95 sur aeropc37
!         16/06/06: CONVECT_YN =1 ou 0 selon qu'on active ou pas la convection
!                   automaitique
!         19/07/06: de nouveau le coef d'asselin est est lu dans un notebook
!         15/09/06: Bienvenue à CST_ADV_VEL
!         06/03/07: commentaires dans messsages
!         26/03/07: initialisation BI_ONOFF (attention en plusieurs etapes)
!                   Reset de RAMPE
!         30/04/07: G est maintenant defini dans coriolis.F à partir de la latitude
!         09/05/07: small2 passe de 1D-10 à 1D-20
!         18/06/07: reset DT_OBC
!         27/12/07: le notebook_optical est desormais obligatoire
!         16/01/08: des noms de fichier pour l'analyse harmonique de la solution
!                   de symphonie
!         17/01/08: IWAVE renommé IWVE
!         21/01/08: ajout de lecture de RELAX_BPC et OBCTYPE_TS dans
!                   notebook_spongelayer
!         22/02/08: En attendant un posttraitement adapté, la version
!                   par defaut supprime le barometre inverse. Pour le
!                   reintroduire, agir directement dans le code (dans la
!                   presente routine). Agir egalement dans external_mode.F90
!         24/02/08: Lecture d'un nouveau parametre, OBCTYPE_P, dans
!                   notebook_spongelayer
!         31/03/08: lecture de HSTEPMIN, HSTEPMAX, NBVSTEPMIN dans
!                   notebook_vertcoord
!         07/04/08: ajout RELAXTYPE_TS pour distinguer la relaxation des
!                   traceur en mode classique du mode "densité inchangée"
!         10/04/08: Par default, la variable BIOBC_TYPE =1
!         11/04/08: Ajout test debug: on ne peut avoir RELAX_TS=0 et
!                   RELAX_BPC non nul etc etc...
!         08/12/08: lecture d'un nouveau notebook_obcforcing
!         24/12/08: IAIRSEADEMO est remplacé par AIRSEAOPTION
!         24/03/09: parallelisation
!         27/04/09: notebook_optical et notebook_spongelayer passe dans NOMFICHIER(15&14)
!         28/04/09: message d'avertissement et arret si NTIDE=0
!         29/04/09: MSKI_Y & MSKJ_X faciliteront les bilans en mode parallele
!         30/04/09: un stop si maree tugo et forcage par grille mere
!         31-05-09  MSKI_Y & MSKJ_X deviennent MASK_I_Z & MASK_J_Z
!         04-06-09  TYPEGRID remplace I2DH
! 2009.2  07-09-09  amenagement lecture notebook_obcforcing
! 2009.3  30-09-09  introduction de before, bef ,now ,aft ,after
!         02-10-09  introduction des increments d'indice génériques
!         15-10-09  echeance beforebefore
!         16-10-09  barometre inverse activé
!         05-11-09  modification notebook_obcforcing
!         15-11-09  tableaux z0 initialiser dans set_parameters.F
!         06-12-09  ajout ipu ipv etc...
! 2010.2  17-12-09  definition deg2rad rad2deg
!         21-12-09  modif dimension parametres nodaux
!         22-12-09  lecture nouveau notebook_tide
!         28-12-09  suppression tidevel
! 2010.3  15-01-10  la multiplication de cst_adv_vel par 0.25 ayant été supprimée
!                   on veille à ce que la valeur en entree de notebook_grid soit reduite
! 2010.4  21-01-10  modifs sur lecture notebook_obcforcing
! 2010.6  04-02-10  airseainfo(2) forcé à zéro. A la place lecture de airseainfo(3)
!         05-02-10  blindage division par sumarchive
! 2010.8  10-03-10  ajout before2 et before3
!         12-03-10  reset tfc1 et tfc2 au time filter 3L
!         17-03-10  pi=dacos(-1.d0)
!         19-03-10  ajout relax T de surface
!         22-03-10  indices de temps definis dans parameter
!         26-03-10  suite du point precedent: attention irelaxsst doit
!                   être traité avant nairsea
! 2010.8  09-05-10  seul le proc 0 ecrit messages
!         18-05-10  si iairsea<0 rappel vers Tsurf
! 2010.9  31-05-31  - les indexations generiques sont maintenant définis dans
!                   parameter
!                   - lecture de verticalgridfile dans notebok_vertcoord
!         13-06-10  Lecture ale_selected dans notebook_vertcoord
! 2010.10 24-06-10  debug s'adapte au fait que kmax remplace nr-1
! 2010.11 16-07-10  re-analyse maree: nouveau nom pour le fichier
! 2010.12 15-09-10  lecture du notebook_wave pour connaitre l'année minimum
!         18-09-10  Message d'alerte si iwve=1 sans etiquette stokes
!         19-09-10  Par defaut z0s=1cm
!         20-09-10  Possibilité de calcul en simple precision
!         27-09-10  securites sur le parametrage
! 2010.13 09-10-10  ajout advsmin advsmax advsstp pour detection du
!                   panache dans advection_scal
!         11-10-10  definition de cellboxfactor1 et cellboxfactor2
! 2010.14 04-11-10  Pour eviter division par zero grille sigma si kmax=1
!         05-11-10  Initialisation des FD time filter parameters
!         15-11-10  suppresion advsmin advsmax qui sont remplacés par upwindriver
! 2010.16 12-01-11  initialisation des coef du filtre temporel d'ordre élevé
! 2010.18 28-01-11  Possibilité de tester le vrai filtre d'Asselin
! 2010.19 13-04-11  Version adaptée à la lecture des champs MERCATOR PSY4V1R3
! 2010.20 15-04-11  kountgraph renommé graphperiod
! 2010.22 22-04-11  Differencier cas nemo psi3 et nemo psi4
!         28-04-11  Differencier cas nemo psi3 et nemo psi4 et sympa
! 2010.23 12-05-11  Initialisation dynrestartfilename
!         18-05-11  - Ajout d'une securité à la lecture de notebook_obcforcing
!                   - Une ligne de plus dans notebook_wave
!                   - bi_onoff lu dans notebook_obcforcing
!         25-05-11  debug filtre FD ALE
! 2010.24 02-09-11  Nom generique pour fichier de grille nemo mercator version p4
! 2010.25 03-02-12  Pas d'eponge si jperiodic
!         09-02-12  Pas d'eponge si iperiodic
!         20-02-12  Ne pas lire notebook_obcforcing is ogcm non utilise
!         25-02-12  allocation dynamique
!         27-03-12  ajout d'un fichier SSH au cas nemo_z
!         04-04-12  debug lecture notebook_airseaflux
!         21-06-12  8 fichiers dans le notebook_obcforcing cas nemo_z
!         29-06-12  notebook_airseaflux: le modele deduit lui-même le nbre d'echeances
!                   par fichier
! S26     20-09-12  - affichages ecran
!                   - rustine pour flux Samuel
!         05-02-13  eponges grille periodique
!         07-02-13  seul proc zero ecrit à l'écran
!         10-02-13  Champs glorys
!         04-04-13  lire le type de fichier ww3 dans le notebook_wave
!         08-04-13  Lire notebook_rivers avec appel set_rivers(0)
!         03-06-13  Plus besoin de lire le notebook_grid qui est maintenant lu dans
!                   module_principal, idem pour reset de pi, rad2deg, deg2rad
!         12-09-13  modifs notebook_tide
!         14-09-13  notebook_visco devient un namelist
!         11-11-13  stop si checkmpi pas compatible avec la config
!         29-11-13  notebook_obcforcing ne contient plus la periodicite
!                   et la date du premier champs.
!         01-12-13  lecture de notebook_eq_state passe dans set_parameters.F90
!         13-01-14  modif sur mask_i_w et cie....
!         24-03-14  ajout flag_nemoffline
!         23-04-14  debug cas notebook_obcforcing sympa
!         01-05-14  Savoir de suite si le status de la procedure offline
!                   pour aiguillage dans construction de grille
!         19-06-14  notebook_eqstate devient un namelist
!         22-06-14  affichage a l'ecran
!         09-07-14  notebook_offline: ne plus deduire la directory des fichiers du
!                   nom de la liste des fichiers offline
!         11-07-14  ajout ogcm_time_shift, dt_obc devient dtobc(:)
!         13-07-14  notebook_offline partitionne en 2 zones (namelist & fichier txt)
!         25-07-14  ajout ofl_surflux
!         12-10-14  ajout initialisation constantes cp_air lv pour pouvoir en disposer
!                   (cas academique) meme si on ne passe pas par formules bulk
!         01-11-14  ajout flag_kz_enhanced
!         03-11-14  ofl_type et ofl_sshtype defini dans notebook_list
!         04-12-14  notebook_vertcoord passe au format namelist
!         06-01-15  definition de mask_i_w, j sur imax+1 jmax+1
!         09-01-15  fichier offline format real
!         15-01-15  suppression d'un stop
!         22-01-15  suppression de lignes inutiles
!         10-02-15  mask_i_u et cie....
!         05-03-15  ajout ofl_tke pour archive de la tke dans fichiers "offline"
!         21-03-15  extension dimensions de mask_i_w et mask_j_w
!         25-03-15  ibl1_advbio,ibl2_advbio,jbl1_advbio,jbl2_advbio definisse
!                   une zone buffer sans advection pour la bio
!         19-04-15  Ajout de cas dans le detecteur d'incoherence des notebook
!         16-05-15  vite lire notebook_graph pour connaitre rapidement dim_varid
!         24-05-15  proprietes optiques initialisees dans initial_main apres
!                   initialisation de h_w
!         15-06-15  tester l'existence du repertoire de la reanalyse de la maree
!                   a la lecture de notebook_tide
!         23-06-15  ajout pgfscheme dans notebook_eqstate
!         10-07-15  expnum (notebook_advection) et rhp_zavr_xy (notebook_eq_state)
!         11-07-15  ajout turbulence echelle moleculaire kmol_m,kmol_h,kmol_s
!         15-07-15  ajout de inv_ekman_depth et z2d3dprofile dans noebook_visco
!         22-07-15  z0s (default surface roughness) defini dans notebook_visco
!         05-08-15  vst continue vs discontinue
!         09-08-15  flag_tke_stab_func permet de choisir les fonctions de stabilite
!                   de la turbulence (canuto, kantha clayson galperin)
!         14-10-15  notebook_obcforcing devient un namelist
!         11-11-15  ajout flag_refstate
!         16-11-15  flag_refstate=0 si kmax=1
!         28-01-16  nouveau filtre 2dt-Ndt
!         10-04-16  si kmax=1 pas d'Etat de reference pour rhp
!         15-04-16  kmaxtidep1=kmaxtide+1
!         24-04-16  - nouveau calcul pour cp
!                   - deplacer le reset de frqtide(kmaxtide+1)
!         26-04-16  Ajout constant_kz dans lecture de notebook_visco
!         13-09-16  pouvoir utiliser iwve=0 avec l'etiquette de compilation stokes
!         19-09-16  ajout offset_sshobc
!         20-11-16  simplification du notebook_advection
!         24-11-16  message A l'ecran
!         05-12-16  Ajout lecture flag_rmnegval dans notebook_advection.f
!         23-12-16  momentum_input_depth est la profondeur de penetration
!                   des inputs de QDM, defini dans notebook_visco
!         16-01-17  upw_hrange1,upw_hrange2 pour definir une zone d'advection upwind
!         25-02-17  if(iadvec_ts==0)upwindriver_t=0.
!         25-03-17  ajout flag_steric_effect dans notebook_eqstate
!         29-03-17  ajout d'un fichier different pour le LSA
!         10-04-17  debug reset de vststep
!         03-05-17  ajout flag_merged_levels dans notebook_vertcoord
!         14-05-17  - biharm_c1 biharm_c2 pour un coef viscositE basEe sur une
!                   proportion du gradient de vitesse et de la vitesse
!                   - biharm_2dfactor permet d'amplifier la viscositE du
!                   mode externe
!         10-09-17  z0b_land et zlevel_land dans notebook_visco
!         01-11-17  schema turbulence zeroequation
!         09-11-17  ajout flag_maxbotstress dans notebook_offline
!         25-01-18  ajout dzsurfmin dans notebook_vertcoord
!         30-01-18  reset dzsurfmin<0
!         20-04-18  call set_parameters_visco deplacE avant wave_reset pour que z0s
!                   soit initialisE avant reset de hsw
!         20-05-18  notebook_spongelayer format namelist
!         02-06-18  ajout zones upwind dans notebook_advection.f
!                   ajout subroutine set_parameters_postgrid
!         08-06-18  ajout flag_offline_binary
!         09-06-18  ajout d'une securite suie au changement de significatio de biharm_c1
!         21-06-18  ajout flag_sponge_txt,sponge_dx_critic
!         01-09-18  constant_kz remplacE par constant_kh et constant_km
!         06-09-18  ajout flag_adve2d,flag_timesplitting_adve_uv    
!         07-10-18  ajout flag_ts_effectivedensity
!         28-10-18  ajout flag_z2dv_outputs  
!         17-01-19  ajout flagspo_i1,flagspo_i2,flagspo_j1,flagspo_j2 
! v246    09-02-19  ajout nouveaux schemas de tke
!         12-02-19  ajout flag_ksloffline=0
! v247    23-02-19  STOP 'kmax=1 requires flag_merged_levels=0 in notebook_vertcoord' 
! v248    04-03-19  tem_validmin,tem_validmax,sal_validmin,sal_validmax !04-03-19
! v255    26-05-19  ajout ofl_reversedime 
! v258    17-07-19  - call set_parameters_visco commentE le 17-07-19 car repart dans initial_main 
!                   car z0_w peut avoir besoin de connaitre h_w
!                   - call wave_reset commentE car appele depuis initial_main apres que z0s ait
!                   ete initialise
! v259    09-09-19  ofl_bio
! v265    28-10-19  ajout d'une securite sur valeur de relax_ext et relax_int sortant de notebook_spongelayer
! v267    16-11-19  iwve/=0
! v269    11-12-19  flag_ogcmtidemixing
! v271    13-12-19  z0b_rivers           
! v273    29-01-19  debug ecriture ecran
!         04-02-20  if(typegrid==typegrid_file)call grid_lonlat2ij_initial
! v275    18-02-20  if(texte90/='none')call dragcoef_initial_z0_local !18-02-20
! v276    04-04-20  analyse harmonique 3D
! v278    17-04-20  ajout du fichier vqs_file dans notebook_vertcoord
! v279    25-04-20  hstepmax=-999
! v281    07-05-20  dz_vertical_incr_fact           
! v284    21-05-20  ajout flag_ofactor dans notebook_advection
! v285    23-05-20  utiliser kmerged plutot que kmin pour calcul upwzone0
! v287    18-07-20  utiliser restartdir_out1 et cie...
! v293    02-12-20  .or.iturbulence==6      
! v299    11-03-21  suppression cn_power
!         19-03-21  ajout flag_negdif
! v300    22-03-21  ajout timescale_negdif dans notebok_advection
!                   ajout vqs_cst1,vqs_cst2,vqs_cst3 dans notebook_vertcoord
! v301    29-04-21  ajout ratio_negdif_ver
! v303    21-07-21  ajouT ratio_bionegdif                                              !21-07-21
!         02-08-21  substep_advbio=1 !02-08-21
! v309    25-08-21  iadvec_ts_hor,iadvec_ts_ver    
! v310    20-10-21  flag_adve3d 
! v314    29-11-21  Ne pas prendre quickest2 si diffusion negative
! v316    20-12-21  valeurs par defaut de parametres tels que flag_ofactor=0 etc...
! v327    19-02-22  ajout flag_ogcm_instab   
! v337    17-03-22  +ogcm_time_lag, nul par defaut, est un decallage de temps ajoutE au temps de l'ogcm
! v339    24-03-22  rho_t est la densite utilisee par le gradient de pression, sur laquelle porte 
!                   l'analyse harmonique 
! v348    27-05-22  ajout flag_z0_macro et z0_0d_macro
! v349    28-06-22  reset dz_over_z0_min=10. 
!         29-06-22  ajout coastal_viscosity
! v350    18-07-22  lignes commentees 
!         20-07-22  ajout quick_filter_points et quick_coef dans notebook_advection
! v358    29-10-22  obctime_order    
! v362    03-01-23  looplimit_hor     
! v367    07-03-23  rho_t desormais allocatE dans module_principal
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! (°v°)
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! (°O°)
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      ! (°L°)
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     ! m°v°m
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

       ioffline=0 ; removetide=0 ; directory_offline='none'
       offlinefile='none' ; ofl_rotation=0 ; ofl_rhp=0 ; ofl_surflux=0
!      flag_kz_enhanced=0 ; ofl_type='short' ; ofl_sshtype='real'
       flag_kz_enhanced=0 ; ofl_type='real' ; ofl_sshtype='real' !09-01-15


      igesig=1 ; isigfile=0 ; hgesig=100. ; pgesig=100. ; ihybsig=0
      nhybsig=0 ; nbvstepmin=5 ; hstepmin=3. ; hstepmax=-999. !25-04-20
      verticalgridfile='none' ; ale_selected=0 ; sigstepgridfile='none'

      iwve=0


!$--------------------------------------------------------------------------------
!$
!$ Generic Modelling section:
!$
!$ time indexes:
!$ SYMPHONIE:
      if(before/=0)stop 'set_parameter erreur1'     !22-03-10
      if(bef/=0)stop 'set_parameter erreur2'
      if(now/=1)stop 'set_parameter erreur3'
      if(after/=2)stop 'set_parameter erreur4'
      if(aft/=2)stop 'set_parameter erreur5'
      if(before2ext/=2)stop 'set_parameter erreur6'
      if(before2/=-1)stop 'set_parameter erreur7'
      if(before3/=2)stop 'set_parameter erreur8'
!$ horizontal indexe increments relative to the centre of cell boxes:
!$ ROMS, POM & SYMPHONIE:
!     ipt=1     ! i forward  increment relative to "c" points
!     imt=0     ! i backward increment relative to "c" points
!     jpt=1     ! j forward  increment relative to "c" points
!     jmt=0     ! j backward increment relative to "c" points
!     ipu=0     ! i forward  increment relative to "x" points             !06-12-09
!     imu=1     ! i backward increment relative to "x" points
!     jpu=0     ! j forward  increment relative to "x" points
!     jmu=1     ! j backward increment relative to "x" points
!     ipv=0     ! i forward  increment relative to "y" points
!     imv=1     ! i backward increment relative to "y" points
!     jpv=0     ! j forward  increment relative to "y" points
!     jmv=1     ! j backward increment relative to "y" points
!$ vertical indexe increments relative to the centre of cell boxes:
!$ SYMPHONIE:
!     kpt=1     ! k upward   increment relative to "c" points
!     kmt=0     ! k downward increment relative to "c" points
!     kpw=0     ! k upward   increment relative to "z" points
!     kmw=1     ! k downward increment relative to "z" points
!$--------------------------------------------------------------------------------

#ifdef checkmpi
      if(par%subcycle/=1) then !------->  !11-11-13
       write(6,*)'checkmpi compilation key not possible' &
                ,' if par%subcycle/=1'
       stop 'STOP set_parameters.F90'
      endif                    !------->
#endif

      flag_nemoffline=0
      if(index(nomfichier(20),'another_model')/=0)flag_nemoffline=1

      write(dynrestartfilename,'(a,i0)') &
           trim(restartdirname)//'restart_output/chanel9_',par%rank

      cellboxfactor1=0.8 ; cellboxfactor2=1.2                         !11-10-10

!     if(ntide==0)stop 'mauvais parameter: choisir ntide=1 ou +'       !29/04/09
!     kountgraph=1.e10                                                !10/07/03
!     graphperiod=1.d10
      unit_r4=1.
      zero=0.
      un=1.
      un_r8=1.

!#MPI ! Pour ne pas sommer 2 fois la zone fantome !10-02-15
      mask_i_w=1
      mask_i_v=1
      mask_i_u=1
      mask_j_w=1
      mask_j_v=1
      mask_j_u=1
      if(obcstatus(ieq1)==0)then    !>>>>>>
        mask_i_w(-1:1)=0            !21-03-15
        mask_i_v(1)=0
        mask_i_u(1:2)=0
      endif                         !>>>>>>
      if(obcstatus(ieqimax)==0)then !>>>>>>
        mask_i_w(imax:imax+2)=0     !21-03-15
        mask_i_u(imax+1)=0
        mask_i_v(imax)=0
      endif                         !>>>>>>
      if(obcstatus(jeq1)==0)then    !>>>>>>
        mask_j_w(-1:1)=0            !21-03-15
        mask_j_v(1:2)=0
        mask_j_u(1)=0
      endif                         !>>>>>>
      if(obcstatus(jeqjmax)==0)then !>>>>>>
        mask_j_w(jmax:jmax+2)=0     !21-03-15
        mask_j_u(jmax)=0
        mask_j_v(jmax+1)=0
      endif                         !>>>>>>


!#MPI ! Pour ne pas sommer 2 fois la zone fantome ! Fin.


! Rest du pas de temps entre 2 états du champs de reference.
! Par defaut de dernier est pris tres grand de ce sorte que
! dv/dt borné. Prend sa vraie valeur dans routine update_obcforcingterms:
!     dt_obc=1.e10                                                     !18/06/07
      dtobc(:)=1.d10                                                   !11-07-14

      t0=   13.0037155                                                 !12/04/05
      s0=   38.4349442                                                 !12/04/05
      alp_t=2.0176679e-04                                              !12/04/05
      alp_s=7.5432844e-04                                              !12/04/05
      rho=  1.0290560e+03                                              !12/04/05
! Lecture de notebook_eqstate                                          !13/11/04
       rhp_zavr_xy=0
       flag_refstate=0                                                 !16-11-15
       flag_steric_effect=0                                            !25-03-17
       tem_validmin=-999. ; tem_validmax=999.                          !04-03-19
       sal_validmin=-999. ; sal_validmax=999.
       open(100,file=nomfichier(17)) ! 'notebook_eqstate
       read(100,nml=notebook_eqstate)
       close(100)
       if(kmax==1.and.flag_refstate/=0) then !00000>
!       flag_refstate=0 !10-04-16
!       if(par%rank==0)write(6,*)                            &
!        'kmax==1.and.flag_refstate/=0 inconsistent choice!' &
!       ,' Setting is changed for kmax=1 & flag_refstate=0 '  
       stop ' Err kmax==1.and.flag_refstate/=0 inconsistent choice'    !16-11-15
       endif                                 !00000>

       if(eos_author/=0.and.eos_author/=3) then !>>>>>
        write(6,*)'Error in notebook_eqstate!'
        write(6,*)'Possibles choices for eos_author are 0 or 3 only!'
        stop 'Stop in set_parameters'
       endif                                    !>>>>>

       t0_base=t0
       s0_base=s0
       rho_base=rho
       alp_t_base=alp_t
       alp_s_base=alp_s


      rhoair=     1.226     ! densité de l'air
! Choix de la valeur pour cp
!https://docs.google.com/document/d/1xJmiLSxhQIl4Z_1RH2h3jjf5AHUw88E00-Go4KOdNbA/edit
!     cp=         3950.     ! Capacite calorifique de l'air Estournel et al, 2009
      cp=         3992.     ! Capacite calorifique de l'air dEduit de figure 4 TEOS-10 !24-04-16
      cp_air=1000.5         ! capacité calorifique de l'air            !12-10-14
      lv=2.5e6              ! chaleur latente de vaporisation (J/kg)   !12-10-14

      small1=     1.d-10                                               !09-07-14
      small2=     1.d-20                                               !09/05/07

      bi_onoff=0                                                       !26/03/07
      rampe=0.                                                         !26/03/07

!____________________________________________________________

      vis=1.4e-6

!_________________________________________
      vststep=0                     !  05-08-15
      flag_merged_levels=0          !  03-05-17
      ihybsig=0 ; nbvstepmin=5 ; sigstepgridfile='none' !25-01-18
      dzsurfmin=-999.                                   !25-01-18
      vqs_file='none' 
      open(100,file=nomfichier(13)) ! 'notebook_eqstate
      read(100,nml=notebook_vertcoord)
      close(100)

      if(ale_selected==1.and.igesig<2) then !-debug->
      write(*,*)'Mauvaise configuration de:'
      write(*,'(a)')trim(nomfichier(13))
      write(*,*)'Si ALE selectionné alors selectionner partial_step'
      write(*,*)'c.a.d. igesig=2'
      stop ' STOP dans set_parameters.F'
      endif                                 !-debug->
      if(vststep==1.and.ihybsig==0) then    !-debug-> !10-04-17
       stop 'Err vststep==1.and.ihybsig==0 in set_parameters.F'
      endif                                 !-debug->

      if(flag_merged_levels==1.and.kmax==1) then !m°v°m>
      STOP 'kmax=1 requires flag_merged_levels=0 in notebook_vertcoord' !23-02-19
      endif


         i2dh=1                                                         !04-06-09
!_________________________________________
!      if(par%rank==0)write(6,'(a11,a60)')'lecture de ',nomfichier(2)
! Simulation 1D vertical I1D=0 sinon I1D=1
! Simulation 2D mode  externe seul      I2DH=0
! Simulation 3D modes externe & interne I2DH=1
!      open(unit=3,file=nomfichier(2))
!      do k=1,9
!        read(3,*)
!      enddo
!        i2dh=1                                                         !04-06-09
!        read(3,*)typegrid ! I2DH                                       !04-06-09
!        read(3,*)flag3d
!        read(3,*)i1,j1,k1
!!.....................................................!
! debugage du parameter:
!#MPI
!       if(i1.ne.iglb.or.j1.ne.jglb.or.k1.ne.kmax) then       !24-06-10
!        if(par%rank==0)write(6,*)
!        if(par%rank==0)write(6,*)'attention anomalie dans la definition'
!        if(par%rank==0)write(6,*)'des dimensions de la grille qui ne'
!        if(par%rank==0)write(6,*)'correspondent pas avec ceux de'
!        if(par%rank==0)write(6,'(a1,a60)')' ',nomfichier(2)
!        if(par%rank==0)write(6,*)'dans ce dernier:   ',i1,j1,k1
!        if(par%rank==0)write(6,*)'et dans namelist_parameter_grid: ',iglb,jglb,kmax
!        if(par%rank==0)write(6,*)'corriger puis relancer.'
!        stop ' dans set_parameters.f'
!        endif
!.....................................................!
!      close(3)
! debug:
!      if(flag3d.ne.1)stop 'flag3d=0 cas non prévu!'
!      if(i2dh.ne.0.and.i2dh.ne.1)stop 'i2dh mal initialisé'


!_________________________________________

! INITIALISATION DE LA MAILLE HORIZONTALE:
! lecture du fichier parametres_grille
!     open(unit=3,file=nomfichier(2))
!       read(3,*)dxb
!       read(3,*)dyb
!     close(3)
!     WRITE(6,*)'OK!'
      dxa=0.5*dxb
      dya=0.5*dyb

!debug
! commente le !15-01-15
!     if(dxb.gt.dyb+1.e-10.or.dxb.lt.dyb-1.e-10)then
!     if(par%rank==0)write(6,*)'dxb different de dyb pas encore validé'
!     stop 'set_parameter'
!     endif


!____________________________________________________________

! Lecture de notebook_airseaflux:                                      !23/08/02

      meteo_sealand_mask=0
      meteo_sealand_file='none'
      irelaxsst=0                                                      !26-03-10
      iairsea=0                                                        !23/08/02
      windfactor=1. ! modif 15/06/01                                   !23/08/02
      do 114 k=1,dim_airsea                                            !23/08/02
      airseainfo(k,1)=0.                                               !23/08/02
      airseainfo(k,2)=0.                                               !23/08/02
      airseainfo(k,3)=0.                                               !04-02-10
  114 continue                                                         !23/08/02

      airsea_year_min=9999              !22-12-09
      airsea_year_max=-9999

      if(par%rank==0)write(6,'(a11,a60)')'lecture de ',nomfichier(7)

      if(index(nomfichier(7),'_s26')/=0) then  !*****************>     !17-09-13
       call airseaflux_readnotebook                                    !17-09-13
      else                                     !*****************>
       stop 'dans set_parameters doute sur notebook_airseaflux'        !22-01-15
      endif                                    !*****************>

      if(iairsea>1)ifb=1
      call airseaflux_allocate
!____________________________________________________________
!      open(100,file=nomfichier(8) ! 'notebook_obcforcing
!      read(100,nml=notebook_obcforcing)
!      close(100)

! Lecture de notebook_obcforcing:
      iobc_f =0
      iobc_ogcm=0
      iobc_wv=0
      iarchive=0
      iobc_lr=0                                                        !25/05/04
      offset_sshobc=0.                                                 !19-09-16
      sumarchive=small1                                                !05-02-10
      ogcm_time_lag=0.                                                 !17-03-22
      obcfile(9,1)='none'                                              !11-12-19
      if(par%rank==0)write(6,'(a11,a60)')'lecture de ',nomfichier(8)
        open(100,file=nomfichier(8)) ! 'notebook_obcforcing            !14-10-15
        read(100,nml=notebook_obcforcing)
       close(100)
      if(obctime_order==4)call principal_allocate_obc !29-10-22
      if(iobc_ogcm==1) then

       if(obcfile(9,1)/='none'.and.flag_refstate==0) then
         write(6,*)'notebook_obcforcing & notebook_vertcoord ' &
       ,'inconsistent regarding the reference state method'    &
       ,'obcfile & flag_refstate= ',trim(obcfile(9,1)),flag_refstate !24-11-16
         stop 'Err 688 ref state inconsistent choice'
       endif

       if(obcfile(9,1)=='none'.and.flag_refstate==1) then
         write(6,*)'notebook_obcforcing & notebook_vertcoord ' &
       ,'inconsistent regarding the reference state method'    &
       ,'obcfile & flag_refstate= ',trim(obcfile(9,1)),flag_refstate
         stop 'Err 687 ref state inconsistent choice'
       endif

      endif
!____________________________________________________________          !01/01/02

      if(par%rank==0)write(6,'(a11,a60)')'lecture de ',nomfichier(14)
      sponge_l=30     
      relax_ext=0.1   
      relax_int=1.    
      relaxtype_ts=1  
      relax_ts=60.    
      relax_es=-999.  
      obcfreeorfix=0   
      relax_bpc=-9999. 
      obctype_ts=0   
      obctype_p=0  
      open(100,file=nomfichier(14)) ! 'notebook_spongelayer !20-05-18
      read(100,nml=notebook_spongelayer)
      close(100)
      relax_ext=max(relax_ext,0.) !28-10-19
      relax_int=max(relax_int,0.) !28-10-19
      if(jperiodicboundary)flagspo_j1=0 !14-01-19
      if(jperiodicboundary)flagspo_j2=0
      if(iperiodicboundary)flagspo_i1=0
      if(iperiodicboundary)flagspo_i2=0

!____________________________________________________________          !01/01/02

      kmaxtide=0
      if(par%rank==0)write(6,'(a11,a60)')'lecture de ',nomfichier(11)
      open(unit=3,file=nomfichier(11)) ! lecture preliminaire du notebook_tide
      read(3,*)
      read(3,*)k0
      if(k0==1) then !------------->
!      do k1=1,26 ; read(3,*) ; enddo

! Lire la valeur de tideana_yesno pour savoir si on alloue les tableaux de l'analyse harmonique 3d !04-04-20
       do k1=1,5 ; read(3,*) ; enddo
       read(3,*)tideana_yesno     
       do k1=1,20 ; read(3,*) ; enddo
! Si tideana_yesno=2 alors activer le flag de de l'analyse harmonique 3d
! puis remettre tideana_yesno=1
       if(tideana_yesno==2) then !m°v°m> !04-04-20
          flag_tide3d_analysis=1
          tideana_yesno=1
          if(par%rank==0)write(6,*)'flag_tide3d_analysis=',flag_tide3d_analysis
!         allocate(rho_t(0:imax+1,0:jmax+1,kmax)) !07-03-23
       endif                     !m°v°m> !04-04-20

  859  read(3,*,end=857)i
       if(i==1)kmaxtide=kmaxtide+1
       do k1=1,10 ; read(3,*) ; enddo ; goto 859
      endif          !------------->
  857 close(3)
      kmaxtidep1=kmaxtide+1 !15-04-16
      if(par%rank==0)write(6,*)'kmaxtide=',kmaxtide

      if(kmaxtide>0) then !tidetidetidetide>

      call principal_allocate_tide

      ktide=0
      open(unit=3,file=nomfichier(11)) ! relecture du notebook_tide
      do k1=1,4 ; read(3,*) ; enddo
      read(3,*)tideforces                                              !21/08/03
      read(3,*) ; read(3,*)
      read(3,*)tideana_yesno                                           !28/07/04
            if(tideana_yesno==2)tideana_yesno=1 !04-04-20
      read(3,*)tideana_spinup                                          !28/07/04
      read(3,*)tideana_delta                                           !28/07/04
      read(3,*) ; read(3,*)
      read(3,*)obc2dtype
      read(3,*) ; read(3,*)
      read(3,*)tide_interpolation
      read(3,*)texte30
                       if(texte30(1:9)=='alongaxis') then
                        tide_flagrotation=1
                       else
                        if(texte30(1:9)=='earthaxis') then
                         tide_flagrotation=0
                        else
                         stop 'tidal uv alongaxis or earthaxis?'
                        endif
                       endif
      read(3,*) ; read(3,*)
      read(3,*) ; read(3,'(a)')texte80(1) ; texte80(2)=texte80(1) ! forcing files directories for ssh & uv
      read(3,*) ; read(3,'(a)')texte80(3) ! directory detiding files
      read(3,*) ; read(3,'(a)')texte80(4) ! directory output (reanalysis files)
      read(3,*) ; read(3,'(a)')texte80(5) ! directory nodal factors files
      read(3,*)

! Tester l'existence du repertoire
       unit_=s_unit(7)
       open(unit=unit_,file=trim(texte80(4))//'/test',err=848) !15-06-15
       close(unit_)

  898 read(3,*,end=885)i ! 1=the following harmonic is selected, 0 otherwise

      if(i==1) then !oooooooooo>
       ktide=ktide+1

       do k1=1,6 !pmx>                                  !17-04-18 

! k1=1 ssh forcing file
! k1=2 uv  forcing file
! k1=3     detiting file
! k1=4     output reanalysis file
! k1=5 LSA forcing file
! k1=6 ssh biases correction file
         read(3,'(a)')texte90 ; k2=index(texte90,' ')-1

         if(k1==1)shortnametide(ktide)=texte90(1:3)

                    nametide(ktide,k1)=trim(texte80(k1))//texte90(1:k2)
         if (k1==5) nametide(ktide,k1)=trim(texte80(1))//texte90(1:k2)
         if (k1==6) nametide(ktide,k1)=trim(texte80(4))//texte90(1:k2)

         k2=len_trim(nametide(ktide,k1))
         if(nametide(ktide,k1)(k2-3:k2)=='none')nametide(ktide,k1)='none'
         if(par%rank==0)write(6,'(a)')nametide(ktide,k1)

       enddo     !pmx>
       if(shortnametide(ktide)(3:3)=='.')shortnametide(ktide)(3:3)=' '
       if(par%rank==0)write(6,'(a,a)')'Tide wave considered: ',trim(shortnametide(ktide))

       read(3,'(a)')texte90 ; k2=index(texte90,' ')-1
       tidenodalfile(ktide)=trim(texte80(5))//texte90(1:k2)
       if(par%rank==0)write(6,'(a)')tidenodalfile(ktide)
       read(3,*)nutide(ktide)
       read(3,*)equitide(ktide)
       read(3,*)
!      if(par%rank==0)write(6,*)'nutide equitide=',nutide(ktide),equitide(ktide)

      else          !oooooooooo>
       do k1=1,10 ; read(3,*) ; enddo ! LIRE VITE SI L'ONDE N'EST PAS SELECTIONNEE
      endif         !oooooooooo>
      goto 898

  885 close(3)

! Ajout d'une onde a tres basse frequence pour analyse 3D seulement: !24-04-16
!     frqtide(kmaxtidep1)=2.*pi/(365.25*10*86400.) ! 10 an
!     frqtide(kmaxtidep1)=2.*pi/(1000.*86400.) ! 10 an
!     frqtide(kmaxtidep1)=0.
      do ktide=kmaxtide+1,kmaxtidep1
       if(ktide==kmaxtide+1)frqtide(ktide)=0.
       if(ktide==kmaxtide+2)frqtide(ktide)=2.*pi/(17.*3600.)
      enddo

      endif               !tidetidetidetide>

      tide_year_min=9999              !22-12-09
      tide_year_max=-9999
      do ktide=1,kmaxtide
       open(unit=3,file=tidenodalfile(ktide))
        read(3,*)
        read(3,*)frqtide(ktide)
        frqtide(ktide)=frqtide(ktide)*deg2rad/3600.
        read(3,*)
        do loop1=1,9999
         read(3,*,end=800)k
         tide_year_min=min0(tide_year_min,k)
         tide_year_max=max0(tide_year_max,k)
        enddo
 800    close(3)
      enddo

!____________________________________________________________



!____________________________________________________________
! activation d'un sous programme simple pour le calcul de
! la tension du vent en fonction de la direction et de la
! vitesse du vent à 10 mètres IWIND=1 sinon IWIND=0
      iwind=1
      if(iairsea.ne.0)iwind=0                                          !02/12/04

      month(1) ='january'
      month(2) ='february'
      month(3) ='march'
      month(4) ='april'
      month(5) ='may'
      month(6) ='june'
      month(7) ='july'
      month(8) ='august'
      month(9) ='september'
      month(10)='october'
      month(11)='november'
      month(12)='december'

      sp_or_db=4
!     UN=1.
#ifdef doubleprecision
      sp_or_db=8
#endif

! Lecture du notebook_rivers
      call set_rivers(0) ! 08-04-13

! lecture de notebook_wave:
!     texte90=nomfichier(22)
      if(par%rank==0)write(6,'(a,a)')'set_parameters va lire: ',trim(nomfichier(22)) !29-01-20
!     open(unit=3,file=texte90)
!     read(3,*)
!     read(3,*)
!     read(3,*)iwve                  ! on / off
!     read(3,'(a)')texte60    !04-04-13
!     k=index(texte60,' ') ; txt_wavemodeltype=texte60(1:k-1)
!     read(3,*)wave_obc_type          !18-05-11
!     close(3)
      open(100,file=nomfichier(22)) ! notebook_wave
      read(100,nml=notebook_wave1)
      close(100)

! call set_parameters_visco commentE le 17-07-19 car repart dans initial_main car z0_w peut avoir besoin
! de connaitre h_w
!     call set_parameters_visco ! notebook_visco !20-04-18!17-07-19
! Note: set_parameters_visco appelE avant wave_reset pour reset hsw avec z0s

! Verification parametrage - debut :
       ksecu=0
#ifdef stokes
! call wave_reset commentE car appele depuis initial_main apres que z0s ait ete initialise
!      call wave_reset ! (including allocation of arrays) !13-09-16!17-07-19
#else
       if(iwve/=0) then !16-11-19
       write(6,*)'Stop: The wave-current effect computation' !22-06-14
       write(6,*)'requires the key compilation "stokes".'
       write(6,*)'Correct the make file and recompile.'
       ksecu=1
      endif
#endif
!      if(iwve==1)call wave_reset ! (including allocation of arrays)
       if(ksecu==1)stop ' STOP dans set_parameters.F90'
! Verification parametrage - Fin.


! Vite connaitre le status de la procedure offline: !01-05-14
       ofl_reversedtime=0 !26-05-19
       open(100,file=nomfichier(20)) ! 'notebook_offline
       read(100,nml=notebook_offline)
       close(100)
       call initial_main_inconsistencies3 !19-04-15

! Vite connaitre la valeur de dim_varid !16-05-15
       call initial_graph_notebook

#ifdef bidon
      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
      write(3,*)'-----------------------------------------------------'
      write(3,*)'subroutine set_parameters:'
      write(3,*)
      write(3,*)'horizontal resolution dxb=',dxb,' meters'
      write(3,*)'grid points in the first horizontal direction : ',imax!13/08/01
      write(3,*)'grid points in the second horizontal direction: ',jmax!13/08/01
      write(3,*)'grid points in the vertical direction         : ',kmax+1  !13/08/01
      if(flag3d.eq.0)write(3,*)'1d oz simulation'
      if(i2dh.eq.0.and.flag3d.eq.1)write(3,*)'2d oxy simulation'
      if(i2dh.eq.1.and.flag3d.eq.1)write(3,*)'3d oxyz simulation'
      if(iwind.eq.1)write(3,*)'subroutine windstress used'
      if(i2dh.eq.1) then
      write(3,*)'minimum bottom drag coef cdseuil=',cdseuil
      endif
      if(flag3d.eq.0)write(3,*)'1d simulation'
      if(istreamf.ge.1)write(3,*)'model_streamf activé'
      if(iairsea.eq.1)write(3,*)'forcage par un model_ meteo activé'

      write(3,*)
      write(3,*)'proprités optiques:'
      write(3,*)'light_kpar1=',light_kpar1,' m'
      write(3,*)'light_rat1=',light_rat1,'=ratio1'
      write(3,*)'light_kpar2_w=',light_kpar2_w,' m'
      write(3,*)'light_rat2=',light_rat2,'=ratio2'
      write(3,*)
      write(3,*)'proprietes schema advection t & s:'                   !06/03/07
      write(3,*)'iadvec_ts_hor=  ',iadvec_ts_hor !25-08-21
      write(3,*)'iadvec_ts_ver=  ',iadvec_ts_ver !25-08-21
      write(3,*)'cst_adv_hor=    ',cst_adv_hor
      write(3,*)'cst_adv_ver=    ',cst_adv_ver
      write(3,*)
      write(3,*)'proprietes schema advection u & v:'                   !06/03/07
      write(3,*)'cst_adv_vel=',cst_adv_vel

      close(3)
      endif                !#mpi-->>-->                       !09-05-10
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
#endif
#endif

      if(flag_asselin==1.and.ale_selected==1) then   !28-01-11
      write(*,*)'Asselin filter can not be selected if ALE selected'
      stop ' STOP in routine set_parameters'
      endif

      return
  848 write(6,'(a,a,a)')'Err 848: directory ',trim(texte80(4)) &
                       ,' not found'
      stop 'Stop 848 in set_parameters'
      end subroutine set_parameters

!.................................................................

      subroutine set_parameters_visco !10-09-17
      use module_principal  ; use module_parallele 
      implicit none

! notebook_visco
      namelist/notebook_visco/iturbulence,emin,epsmin,z0b,cdseuil & !14-09-13
                             ,tke_surf,convect_yn,assel0,cdb_2dh  &
                             ,texte250,kmol_m,kmol_h,kmol_s       & !24-05-15 !11-07-15
                             ,zprofile2d3d,inv_ekman_depth,z0s    & !15-07-15 !22-07-15
                             ,flag_tke_stab_func,ctke1,ctke2      & !09-08-15
                             ,flag_linearfric,coef_linearfric     & !19-01-16
                             ,constant_km                         & !01-09-18
                             ,constant_kh                         & !01-09-18
                             ,momentum_input_depth                & !23-12-26
                             ,z0b_land,zlevel_land                & !10-09-17
                             ,z0b_rivers                          & !13-12-19
                             ,texte90                             & !18-02-20
                             ,flag_z0_macro,dz_over_z0_min        & !27-05-22
                             ,coastal_viscosity                     !29-06-22

      if(par%rank==0)write(6,'(a11,a60)')'lecture de ',nomfichier(9)
      convect_yn=0 ! par defaut par de convection automatique          !16/06/06
      assel0=0.1   ! coef d'Asselin par defaut                         !19/07/06

! Schema de Gaspar avec les constantes recommandees par l'auteur:
      ctke1=0.1
      ctke2=0.7
! Schema de Gaspar avec les constantes qui permettent de retrouver l'equation
! de Craig et banner 1994 avec (SM,B)=(0.39,16.6)
!     ctke1=0.55  ! Km=l.q.SM=l.sqrt(2.tken).SM=l.sqrt(tken).sqrt(2).SM ==>ctke1=sqrt(2).SM=0.55
!     ctke2=0.17  ! epsilon=2.sqrt(2.tken).tken/(B.l)=(tken**3/2)*(1/l)*2.sqrt(2)/B ==> ctke2=2.sqrt(2)/B

      texte90='none'   !18-02-20
      texte250='none' ; z0b=-9999. ; kmol_m=0. ; kmol_h=0. ; kmol_s=0. !11-07-15
      zprofile2d3d=0  ; inv_ekman_depth=0.1    ; z0s=0.01  ; cst_c0cub=-999.
      flag_linearfric=0 ; coef_linearfric=0. ; constant_km=1.e-6 ; constant_kh=1.e-7
      flag_tke_stab_func='kc94' ! ou 'ca01' ou 'gp88'
      momentum_input_depth=1.   !23-12-16
      z0b_land=-999.            !10-09-17
      flag_z0_macro=0     !27-05-22
      dz_over_z0_min=10.  !28-06-22
      coastal_viscosity=0. !29-06-22
      open(3,file=nomfichier(9))    ! lecture du namelist "notebook_visco" !14-09-13
      read(3,nml=notebook_visco)
      close(3)
      if(z0b_land==-999.)z0b_land=z0b
       if(z0s==0.)stop 'Err 512 z0s=0 leads to a division by zero'
       if(flag_tke_stab_func=='kc94' &
      .or.flag_tke_stab_func=='gp88')cst_c0cub=0.5544**3 ! table 2 Warner et al OM 2005
       if(flag_tke_stab_func=='ca01')cst_c0cub=0.5270**3 ! table 2 Warner et al OM 2005
       if(cst_c0cub==-999.)stop 'Err 515 flag_tke_stab_func not defined'

! Bottom roughness:
! l'appel A dragcoef est desormais fait apres l'initialisation de la
! bathy !10-09-17
! Les lignes suivantes ont ete commentees le 10-09-17
      if(texte250/='none') then !>>>>
       call dragcoef_initial_z0_file
      else                      !>>>>
       if(z0b==-9999.)stop 'Err 595 set_parameters z0b==-9999'
       call dragcoef_initial_z0_h
      endif                     !>>>>

! Ajouter des valeurs ponctuelles (avec bulle d'influence)
      if(texte90/='none')call dragcoef_initial_z0_local !18-02-20

      if(    iturbulence==0         &
         .or.iturbulence==2         &
         .or.iturbulence==4         & !09-02-19
         .or.iturbulence==5         & !09-02-19
         .or.iturbulence==6         & !02-12-20
                       )call turbulence_gaspar_allocate_init !26-04-16
      if(iturbulence==1)call turbulence_k_eps_allocate_init
      if(iturbulence/=0.and.iturbulence/=1.and.iturbulence/=2  &
                                          .and.iturbulence/=3  & !01-11-17
                                          .and.iturbulence/=4  & !09-02-19
                                          .and.iturbulence/=5  & !09-02-19
                                          .and.iturbulence/=6) & !02-12-20
      stop 'turbulence sch. missing'

! coef d'Asselin: desormais dans time_step.F90
      end subroutine set_parameters_visco

!.................................................................

      subroutine set_parameters_postgrid
      use module_principal  ; use module_parallele ; use module_grid
      implicit none
! Lecture des notebook necessitant de connaitre les details de la grille 
! (par exemple si call latlontoij

! notebook_advection
      namelist/notebook_advection/iadvec_bio,cst_adv_hor           &
       ,cst_adv_ver,cst_adv_vel,substep_advbio                     &
       ,ibl1_advbio,ibl2_advbio,jbl1_advbio,jbl2_advbio            & !25-03-15
       ,expnum,flag_rmnegval                                       & !10-07-15!05-12-16
       ,upw_hrange1,upw_hrange2                                    & !16-01-17
!      ,biharm_c1,biharm_c2                                        & !14-05-17
       ,biharm_2dfactor                                            & !14-05-17
       ,texte250                                                   & !02-06-18
       ,flag_adve2d,flag_timesplitting_adve_uv                     & !06-09-18
       ,flag_ts_effectivedensity                                   & !07-10-18
       ,flag_ofactor                                               & !21-05-20
       ,upwindzone_file                                            & !22-05-20
       ,flag_ts_quicklim                                           & !22-05-20
       ,flag_negdif_hor                                            & !04-10-21
       ,flag_negdif_ver                                            & !04-10-21
       ,ratio_negdif_hor                                           & !04-10-21
       ,ratio_negdif_ver                                           & !04-10-21
       ,ratio_bionegdif                                            & !21-07-21
       ,iadvec_ts_hor,iadvec_ts_ver                                & !25-08-21
       ,flag_adve3d                                                & !20-10-21
       ,quick_coef                                                 & !20-07-22
       ,quick_filter_points                                        & !20-07-22
       ,looplimit_hor                                                !03-01-23

      itimets=0                                                      !24/03/03
      expnum=1.
      cst_adv_hor=1.
      cst_adv_ver=1.
      cst_adv_vel=1.
      biharm_2dfactor=1.  !20-12-21
      flag_negdif_ver=1   !20-12-21
      flag_negdif_hor=1   !20-12-21
      ratio_negdif_hor=1. !20-12-21
      ratio_negdif_ver=1. !29-04-21
      iadvec_ts_hor=1     !20-12-21
      iadvec_ts_ver=1     !20-12-21
      ratio_bionegdif=1.  !21-07-21
      texte250='none'     !02-06-18
      flag_rmnegval=0 ! Par default on ne s'oppose pas aux valeurs negatives !05-12-16
      flag_ofactor=0      !20-12-21
      substep_advbio=1    !02-08-21
      flag_adve2d=2       !20-12-21
      flag_adve3d=1       !20-12-21
      quick_coef=0.125 ! 1/8. !18-07-22
!     quick_coef=0.166 ! 1/6. !18-07-22
      quick_filter_points=3   !20-07-22
      if(par%rank==0)write(6,'(a11,a60)')'lecture de ',nomfichier(5)
      open(3,file=nomfichier(5))    ! lecture du namelist "notebook_advection" !29-04-14
      read(3,nml=notebook_advection)
      close(3)

!     if(flag_negdif_hor>0) then !m°v°m> !29-11-21
!       if(iadvec_ts_hor==iadvec_ts_hor_quickest2) then !>>>
!        stop 'STOP iadvec_ts_hor_quickest2 and flag_negdif_hor>0'
!       endif                                           !>>>
!     endif                      !m°v°m> !29-11-21

! commentees le !18-07-22
!     if(flag_negdif_ver>0) then !m°v°m> !29-11-21
!       if(iadvec_ts_ver==iadvec_ts_ver_c2) then !>>>
!        stop 'STOP iadvec_ts_ver_c2 and flag_negdif_ver>0'
!       endif                                           !>>>
!     endif                      !m°v°m> !29-11-21

      if(par%rank==0) then !0000000>
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
       write(3,*)'-----------------------------------------------------'
       write(3,*)'subroutine set_parameters:'
       write(3,*)
       write(3,*)'iadvec_ts_hor=   ',iadvec_ts_hor !25-08-21
       write(3,*)'iadvec_ts_ver=   ',iadvec_ts_ver !25-08-21
       write(3,*)'iadvec_bio=      ',iadvec_bio
       write(3,*)'substep_advbio=  ',substep_advbio
       write(3,*)'cst_adv_hor=     ',cst_adv_hor ! 0.=0%  1.=100% du limiteur
       write(3,*)'cst_adv_ver=     ',cst_adv_ver ! 0.=0%  1.=100% du limiteur
       write(3,*)'cst_adv_vel=     ',cst_adv_vel ! 0.=0%  1.=100% du limiteur
       write(3,*)'expnum=          ' ,expnum
       write(3,*)'biharm_2dfactor= ',biharm_2dfactor
       write(3,*)'flag_negdif_hor  ',flag_negdif_hor
       write(3,*)'flag_negdif_ver  ',flag_negdif_ver
       write(3,*)'ratio_negdif_hor ',ratio_negdif_hor
       write(3,*)'ratio_negdif_ver ',ratio_negdif_ver
       write(3,*)'ratio_bionegdif  ',ratio_bionegdif
       close(3)
      endif                !000000>

      if(flag_timesplitting_adve_uv==1) then !w°v°w>
       dbefore=-1
       dnow=0
       dafter=1
      endif                                  !w°v°w>


!...........
! upwindriver_t=0 <==> schema d'advection upwind !25-02-17
      if(iadvec_ts_hor==iadvec_ts_hor_upwind)upwindriver_t=0.
      

!     if(cst_adv_vel>2.)    &
!     stop 'Reduire cst_adv_vel dans notebook_advection' !14-01-10

! zones upwind: !02-06-18
      if(upwindzone_file/='none') then !pmx> !02-06-18

        open(unit=3,file=upwindzone_file)
   539  continue

        read(3,'(a)',end=540)texte30      
        backspace 3


        if(index(texte30,'3D')==0) then !--->
          read(3,*)x1,x2,i0      
        else                            !--->
          read(3,*)x1,x2,k0,i0   
        endif

        flag_stop=1
        if(index(texte30,'indice')/=0) then !>>>>
         i1=nint(x1)-par%timax(1) ; j1=nint(x2)-par%tjmax(1) ; flag_stop=0
        endif                               !>>>>
        if(index(texte30,'lonlat')/=0) then !>>>>
         x1=x1*deg2rad ; x2=x2*deg2rad
         call latlontoij(x1,x2,'loc') ; i1=nint(deci) ; j1=nint(decj) ; flag_stop=0
        endif                               !>>>>
        if(flag_stop==1) stop 'Err 8 unrecognized format'

        if(index(texte30,'3D')==0) then !-cas-upwindriver->

         do i=i1-i0,i1+i0 ; do j=j1-i0,j1+i0
          if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then !------->
                  upwindriver_t(i,j)=                              & !02-06-18
          max(min(upwindriver_t(i,j)                               &
                 ,sqrt(real(i-i1)**2+real(j-j1)**2)/(0.5*i0)-1.),0.)
          endif                                              !------->
         enddo            ; enddo

        endif                           !-cas-upwindriver->

        goto 539
   540  close(3)        
      endif                    !pmx>

      if(typegrid==typegrid_file)call grid_lonlat2ij_initial !04-02-20

      end subroutine set_parameters_postgrid

!.................................................................


      subroutine set_parameters_postkmin
      use module_principal  ; use module_parallele ; use module_grid
      implicit none

! Reset upwzone0_t
      upwzone0_t=1.

! zones upwind: !02-06-18
      if(upwindzone_file/='none') then !pmx> !02-06-18

        open(unit=3,file=upwindzone_file)
   539  continue

        read(3,'(a)',end=540)texte30      
        backspace 3

        if(index(texte30,'3D')==0) then !--->
          read(3,*)x1,x2,i0      
        else                            !--->
          read(3,*)x1,x2,k0,i0   
        endif

        flag_stop=1
        if(index(texte30,'indice')/=0) then !>>>>
         i1=nint(x1)-par%timax(1) ; j1=nint(x2)-par%tjmax(1) ; flag_stop=0
        endif                               !>>>>
        if(index(texte30,'lonlat')/=0) then !>>>>
         x1=x1*deg2rad ; x2=x2*deg2rad
         call latlontoij(x1,x2,'loc') ; i1=nint(deci) ; j1=nint(decj) ; flag_stop=0
        endif                               !>>>>
        if(flag_stop==1) stop 'Err 8 unrecognized format'

        if(index(texte30,'3D')/=0) then !-cas-upwzone0-> !21-05-20

        if(index(texte30,'surf')/=0) then !-surface-case->

         do i=i1-i0,i1+i0 ; do j=j1-i0,j1+i0
          if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then !------->
           do k=kmax,max(kmax-k0+1,1),-1
                    upwzone0_t(i,j,k)=                              & !21-05-20
            max(min(upwzone0_t(i,j,k)                               &
                   ,sqrt(real(i-i1)**2+real(j-j1)**2)/(0.5*i0)-1.),0.)
           enddo ! k loop
          endif                                              !------->
         enddo            ; enddo

        else                              !-bottom-case->

         do i=i1-i0,i1+i0 ; do j=j1-i0,j1+i0
          if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then !------->
!          do k=1,min(kmin_w(i,j)+k0-1,kmax)
           do k=1,min(kmerged_t(i,j)+k0-1,kmax)  !23-05-20
                     upwzone0_t(i,j,k)=                              & !21-05-20
             max(min(upwzone0_t(i,j,k)                               &
                    ,sqrt(real(i-i1)**2+real(j-j1)**2)/(0.5*i0)-1.),0.)
           enddo ! k loop
          endif                                              !------->
         enddo            ; enddo
        endif                             !-bottom-case->

        endif                           !-cas-upwzone0->



        goto 539
   540  close(3)        
      endif                    !pmx>

      end subroutine set_parameters_postkmin

!.................................................................
