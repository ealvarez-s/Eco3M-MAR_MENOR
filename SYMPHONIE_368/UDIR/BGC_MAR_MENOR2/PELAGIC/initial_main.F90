      subroutine initial_main
!______________________________________________________________________
! SYMPHONIE ocean model
! release S26 - last update: 29-10-18
!______________________________________________________________________

      use module_principal ; use module_wave ; use module_parallele
      use module_drifter ; use module_airseaflux ; use module_offline
      use module_modeanalysis ; use module_my_outputs ; use module_grid
      use module_curvgrdtoolbox ; use module_q
      implicit none
#ifdef synopsis
       subroutinetitle='initial_main'
       subroutinedescription='Driver of the initial state subroutines'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!______________________________________________________________________
! Version date      Descriptions des modifications:
!         09/08/01: amelioration de la procedure d'initialisation des
!                   fleuves.
!         10/08/01: dans le cas ISTREAMF=1 (model_streamf activé) on
!                   n'initialise pas le vent et les débits des fleuves
!         13/08/01: ajout de call model_out(0)
!         21/08/01: bouton IOBC_LR, appel à model_in
!         29/08/01: nouvelles conditions d'appel à MODELE_OUT(0)
!         09/10/01: creer grid3D.out: fichier comportant les caracteristiques
!                   3D de la grille: latitude longitude profondeur de chaque
!                   point de grille _Z
!         27/10/01: apparition d'un argument passé dans INIT_WITH_OBC
!         30/10/01: initialisation de T & S aux sources des fleuves
!         08/03/02: on n'oublie pas d'initialiser les points i=0 i=imax+1
!                   j=0 j=jmax+1 pour T et S
!         09/03/02: appel à tgt_const
!         14/03/02: applel à grid_out et tgt_const apres z_averaged
!         15/07/02: amenagements pour nouvelle filiere de calcul
!                   des flux atmospheriques (formules bulk)
!         23/08/02: mise à jour à partir de la version 2003. L'appel
!                   à initial_tide est réintroduit pour être promu à
!                   2003.
!         02/10/02: bienvenue à CALL INIT_STATE_EQ
!         07/10/02: appel à CALL BOUEES(0)
!         14/10/02: correction d'un bug introduit malencontreusement
!                   qq semaines avant la presente correction. L'argument
!                   ds CALL DTVAR_OBC avait mysterieusement disparu. IL
!                   est reintroduit
!         06/11/02: Position de sites géographiques sur la grille du
!                   model_: ajout de CALL ATLAS.
!         16/12/02: ajout de CALL CWAVE_INT pour initialiser les vitesses
!                   de propagation pour la condition aux limites du mode
!                   interne
!         27/12/02: vers une grille hybride: call mix_sig_step
!         27/01/03: appel filiere offline
!         28/01/03: appel à MIX_SIG_STEP
!         25/03/03: suppression CALL ANIMATION
!                   Choix d'appel sur mix_sig_step
!         05/06/03: appel à coriolis déplacé avant appel à bathycote.
!                   Motif: l'imbrication de la bathy "mère" dans bathycote
!                   suppose d'interpoller avec la filiere de model_in
!                   qui requiere les tableaux lat lon. CQFD.
!                   Passage d'un argument dans initial_sponge
!         28/07/03: appel à initial_graph
!         14/08/03: pas de KSERIE devant grid_out ou tgt_const
!         26/08/03: bienvenue à CALL NEST_INOUT
!         29/09/03: bienvenue à CALL ANALYSEHARMONIQUE(0)
!         01/10/03: ce dernier est en fait déplacé (avant CALL INIT_GRAPH)
!         09/01/04: cas IAIRSEA=3
!         21/04/04: si IOBC_AF=2 on utilise un fichier "particulier" à
!                   l'initialisation.
!         28/05/04: amenagement pour pouvoir faire quelques modifs dans
!                   les notebooks apres la lecture d'un fichier restart
!         29/06/04: des argument dans §ALL OBC_SCAL
!         08/07/04: introduction des routines "ondelettes" de Francis
!         27/07/04: CALL RIVER_UPD(1) doit être appelé plus tôt pour
!                   que les températures des fleuves variables avec le
!                   temps soient bien initialisées au moment de l'appel
!                   à set_rivers(3)
!         13/04/05: Relire notebook_nesting apres lecture fichier restart.
!                   Permet de redefinir les sorties pour de nouvelles
!                   filles qui n'avaient pas été prévue au départ du
!                   run de la maman.
!         21/04/05: appel à set_river(3) déplacé. Celui est appelé qu'à
!                   condition que river_upd soit egalement appelé
!         26/04/05: dans interrupteur on peut maintenant demander un
!                   fichier restart sans arreter le model_
!         22/07/05: initial_sponge appelé quelque soit sponge_l
!         29/07/05: Attention dans le cas d'une equation prognostique
!                   il ne faut pas appeler les obc sur T et S sinon
!                   T et S avancent artificiellement d'un pas de temps
!                   de trop. Pour assurer l'initialisation des variables
!                   hors domaine (advection upstream) initial_with_obc.F a
!                   été modifié en consequence.
!         14/02/06: un argument apparait dans l'appel à dragcoef
!         16/06/06: restart special biologie
!         25/07/06: lecture du notebook_offline avant initialisation du
!                   forcage atmospherique car forcage meteo reduit si
!                   filiere offline activee.
!         05/03/07: possibilite de sortir diverses grilles
!         26/03/07: ajout appel add_bi (initialisation du BI)
!         03/04/07: suppression definitive model_in et model_out
!         07/06/07: bornes min et max pour i & j passes en argument de
!                   z_averaged
!         27/08/07: Les appels  aux subroutines de "re-initialisation" apres
!                   hot_restart sont annules
!         02/10/07: Blocage en cas de demande simultanee de filiere
!                   offline et nesting
!                   Appel à offline(-1) pour connaitre rapidement la valeur
!                   de IOFFLINE
!         22/02/08: les lignes pour l'analyse de la sse ogcm (commentées par
!                   defaut) sont revues
!         18/04/08: KOUNT0_RST est le KOUNT correspondant à la date du dernier
!                   redemarrage.
!         04/12/08: Possibilite de stopper le run a la fin d'initier.
!                   Le parametre RUN_OPTION se regle dans notebook_time
! 2009.2  03-09-09: quick_initial passant par les tableaux obc on deplace
!                   l'appel à initial_with_obc
! 2009.3  30-09-09: Initialisation des nouveaux facteurs d'echelle verticale
!         05-10-09: deplacement de cellbox_thickness apres mix_sig_step
!         19-10-09: des arguments passés dans cellbox_thickness
!         12-11-09: creation d'une subroutine qui passe en revue le fichier
!                   variables.txt et qui stoppe le programme si anomalie
!         15-11-09: initier renommé initial_main
!         19-12-09: time_step renommé initial_time_step
! 2010.2  20-12-09: - subroutine bathycote renommée initial_mask_and_bathy
!                   - subroutine coriolis renommée initial_lonlat_dxdy_coriolis
! 2010.6  02-02-10: renomme lon_t lat_t
! 2010.7  04-03-10: initialiser le wet mask
! 2010.8  09-05-10  seul le proc 0 ecrit interrupteur
! 2010.9  29-05-10  suppression dsig_u dsig_v
!         06-06-10  model_wave renommé waveforcing
! 2010.10 23-06-10  hot_restart renommé dyn_restart et bio_restart
! 2010.14 18-11-10  archiver la grille à la fin de l'initialisation
! 2010.20 16-04-11  Calculs sur la base d'un temps en secondes
!         19-04-11  initial_time_step renommé time_step_initial
! 2010.23 24-05-11  initialiser dz_t(i,j,k,-1)
! 2010.25 01-02-12  routine drifter remplace routine bouee
!         21-02-12  modif argument dans dragcoef
!         23-02-12  appel a waveforcing remplace par wave_initial
!         09-03-12  afin de pouvoir regarder l'etat initial des vagues
!                   l'appel a graph_out quand iteration3d=0 est place
!                   dans model_3d
! 2010.25 08-06-12  use module_drifter
! S25.4   30-06-12  use module_airseaflux
!         24-07-12  iteration3d_restart remplace iteration3d_restart
! S26.1   24-11-12  module_offline
!         10-02-13  appel a airseaflux_driver
!         15-03-13  appel a my_outputs
!         01-06-13  initial_lonlat_dxdy_coriolis remplace par module_grid
!         13-09-13  detecteur de mauvais parametrages des notebook
!         13-10-13  offline_write_grid_2 devient offline_write_grid(1)
!         17-11-13 obc2dtype devient real
!         19-01-14 Sauver une copie de ioffline apres lecture du fichier restart
!         05-05-14 ajout d'une verification de coherence de notebook_grid
!         11-07-14 supprime nest_inout
!         23-07-14 suppression routine files_name, notebook_list lu dans routine main
!         22-08-14 Enlever le stop si advection upg et grille deformee
!         25-10-14 Ecrire le fichier grid.nc egalement dans OFFLINE si full proc
!         28-11-14 call z_levels deplace apres cell_box et argument 0 au lieu de 1
!         08-01-15 Pas besoin de passer par read_ogcm si ioffline=2
!         23-01-15 Des actions realisees dans main.F90 passe au debut de initial_main
!         17-02-15 Message d'erreur si start time > end time
!         02-03-15 iglb>=4 jglb>=4 sinon stop
!         19-04-15 ajout de cas dans le detecteur d'incoherence des notebook
!         16-05-15 ajout de cas dans le detecteur d'incoherence des notebook
!         24-05-15 appel a optics_notebook_optical apres connaissance de h_w
!         10-06-15 call initial_main_ww3gridfile: gives ww3 input grid file using S
!                  curvilineatr grid
!         08-07-15 suite point precedent, correction taille de boucle i
!         26-08-15 fichier interrupteur ecrit dans repertoire tmp
!         27-10-15 - initialisation des parametres de la maree avant l'interpolation 
!                  de l'ogcm pour mixing de l'interpolation
!                  - Calcul etat de reference pour PGF
!         09-11-15 reset des tableaux de reference a zero si methode etat de 
!                  reference non utilisee
!         11-11-15 ajout flag_refstate
!         18-11-15 suite du point precedent. Modif pour que rhcref soit calcule
!                  meme si flag_refstate==0                  
!         26-11-15 aide et avertissements au cas modele 1DV
!         27-11-15 Cas model 1DV : verifier la continuite
!         29-11-15 ajout use module_curvgrdtoolbox
!         23-08-16 alarmes debug cas 1DV
!         29-09-16 ajout clefs oasis
!         02-01-17 un commentaire
!         15-04-17 possibilite de definir z0b dans les reservoirs des fleuves
!         04-06-17 plus de cas dans initial_main_inconsistencies2
!         10-07-17 autorisation du k-epsilon en coordonnee s-z
!         18-08-17 if(flag_nh2d==1)call q_allocate
!         10-09-17 parametres de turbulence et frottement definis apres lon,lat,h,mask
!         03-10-17 ajout iteration2d_begin
!         07-11-17 ajout dtmultiple
!         17-11-17 amenagement cas ioffline=2
!         20-11-17 - ajout appel principal_allocate_last
!                  - ajout d'une verification sur valeur du parametre bulk_scheme
!         09-12-17 deplacement de cwave_init avant q_initial
!         20-04-18 call set_parameters_visco commentE car appelE dans set_parameters
!         14-05-18 call atlas est deplace dans initial_mask_and_bathy pour pouvoir
!                  utiliser glob_mask
!         02-06-18 ajout call set_parameters_postgrid
!         17-06-18 initial_main_bin2netcdf est un convertiseur de fichiers binaires en netcdf
!         05-07-18 nf_clobber+ NF_64BIT_DATA permet de creer des fichiers pour de grandes grilles
!         17-08-18 cas flag_nh3d_uv > 1 possible
!         29-08-18 if(flag_timesplitting==0.and.rhp_zavr_xy==1) STOP
!         03-10-18 ajout fplan1_grid fplan2_grid
!         16-10-18 modif sur test fplan2_grid
!         29-10-18 allocate module_nh avant call q_initial
!...............................................................................
!    _________                    .__                  .__             !    (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

      call initial_main_inconsistencies1 !02-03-15

      initial_main_status=0 !23-02-12

!     call initial_main_bin2netcdf !17-06-18

! Array allocation:
      call principal_allocate       !23-01-15

! Check that the RDIR/xxx/tmp directory is empty:
      call initial_main_checktmpdir !23-01-15

      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
      open(unit=3,file=trim(tmpdirname)//'messages')
      write(3,*)'               messages divers'
      write(3,*)'..............................................'
      close(3)
      endif                !#mpi-->>-->                       !09-05-10
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
#endif

! mise à zéro de certains tableaux (assure des départs propres dans le
! cas d'un ensemble de simulations à l'interieur d'un seul run)
      call reset

! initialisation de la maille horizontale, des constantes physiques...
      call set_parameters

! initialisation des latitudes et longitudes des points de grille       !25/06/03
! et le coef. de coriolis:
      call grid_driver ! lon, lat, dxdy, coriolis                       !01-06-13

! notebook  & initialisation necessitant que la grille soit connue (donc apres call grid_driver)
      call set_parameters_postgrid                                      !02-06-18

! initialisation du masque terre/mer + bathymetrie:
      call initial_mask_and_bathy                                       !20-12-09


! Attirail A grille curvi !02-01-17
!     call curvgrdtoolbox_driver(0)
!     call curvgrdtoolbox_driver(1)
!     call curvgrdtoolbox_driver(2)
!     stop 'cl'

! lecture notebook_visco (parametres de turbulence, frottement, etc...)
!     call set_parameters_visco  !10-09-17 ! commentE le 20-04-18 car appelE dans set_parameters

! Verification de la fonction de passage des lon,lat vers indices i,j
! La grille horizontale generalisee necessite mask_t donc appel place
! apres lecture de mask_t:
!     call grid_check_lonlat2ij !29-12-13

      call q_initial_beforetimestep !20171107

! configuration des schemas d'avancée dans le temps, pas de temps....
! définir les dates de début et de fin de simulation....
!     call initial_time_step                                            !19-12-09
      call time_step_initial                                            !19-04-11

! initialisation de la grille verticale (positionner les niveaux sigma)
      call sigma_levels ! attention cette routine peut modifier h_w et donc
                        ! l'appel à z_thickness doit être placé apres et non avant

! hybrid sigma/step grid generator: (compute kmin_w kmin_u kmin_v)
!     call sigstepgrid_driver(ihybsig)  passe dans sigma_levels

! initialisation de l'epaisseur de la colonne d'eau                     !27-09-10
      call z_thickness(0,2)
! computes the thickness of the cell boxes
      call cellbox_thickness(0,-1,2)                           !24-05-11
      call z_levels(0)                                         !28-11-14

! Initialisation des vitesses de propagation pour la condition aux limites ouvertes du mode interne
      call cwave_int !16/12/02 deplace avant q_initial le !09-12-17

      if(flag_nh3d>=1)call q_allocate       !18-08-17 !17-08-18
      call q_initial !29-10-18

! initialisation du coef de frottement sur le fond
      call dragcoef_initial !21-02-12

! Initial optical properties (ici car necessite connaissance de h_w)
      call optics_notebook_optical !24-05-15

! initialisations simples de tableaux:
! (placé vers la fin pour beneficier de l'initialisation
!  des parametres geographiques de la grille (z, sigma latitude
!  longitude etc...) pour initialiser la stratification)
      call quick_initial

! initialisation de la couche eponge:
      call initial_sponge(0)                                              !22/07/05

!     call z_averaged(1,imax+1,1,jmax+1,1)                             !07/06/07

! initialisation des drifter 3D                                        !01-02-12
      call drifter_initial

!c initialisation du model_ avec des champs de grande echelle:

! Initialisation des fleuves à partir de fichiers de débit:
      call river_upd(1)                                   ! déplacé le !27/07/04
      call set_rivers_reservoir_z0 !15-04-17

      call offline_inout(-1) ! Vite connaitre IOFFLINE                 !02/10/07
!     call nest_inout(0) ! lire notebook_nesting                       !26/08/03
!     call nest_inout(2) ! lire un fichier "maman" et l'interpoller    !30/08/03

! Initialisattion du modèle de houle
#ifdef key_oasis_symp_ww3
#else
        call wave_driver                                                 !04-03-13
#endif                                             

! Initialisation des parametres de la marée
      call initial_tide                                                !27-10-15

      if(iobc_ogcm==1)call read_ogcm_fields(1) !08-01-15

      call update_obcforcingterms(0)  ! interpole au temps t initial les !14/10/02
                                      ! echeances de l'analyse
      call initial_with_obc(1) ! initialise le model_ avec les tableaux obc !03-09-09

      call offline_inout(0) ! filiere offline                          !25/07/06

! Initialisation du forçage atmospherique à partir d'un fichier météo:
#ifdef  key_oasis_symp_surfex
#else
      call airseaflux_driver(1)        !10-02-13
#endif

! Ajout du BI à l'initialisation:
      call add_bi(0,0,0)                                               !26/03/07

      if(flag_nh3d==1)call quick_nh !17-08-18

! Initialisation de T et S aux sources des fleuves:
      call set_rivers(3)                                               !21/04/05

! Initialisation des parametres de l'equation d'etat:                  !02/10/02
      call initial_state_eq                                               !02/10/02

! Initialisation éventuelle d'un modèle biologique ou de transport sédimentaire
      call strada(0)

! Initialisation des parametres de la marée
!     call initial_tide                                                   !23/08/02

! Modes propres 3D, harmoniques ondes internes !26-02-13
      call modeanalysis_roadmap

! Initialisattion du modèle de houle
!     call wave_initial

!$ computes the wetdry mask and cancels airsea fluxes in dried aeras:   !04-03-10
      call wetdry_mask_airseafluxes

! Position de sites géographiques sur la grille du model_:             !06/11/02
!     call atlas  !commentE le 14-05-18 

      call initial_graph                                               !28/07/03

      call offline_write_grid(1)                                       !13-10-13
!     if(mpi_hole_plugging=='none')call offline_write_grid(2)          !25-10-14

      call my_outputs_driver ! allocater les tableaux avant lecture chanel9 !15-03-13

      call initial_main_ww3gridfile ! ww3 input grid file using S curvilineatr grid !10-06-15

! Reference state for PGF
      if(ioffline/=2)call pressure_gradient_initref  !27-10-15  !17-11-17
#ifdef bidonref
      if(flag_refstate==0) then     !>>>>        !18-11-15
       salref_t=0.  ; temref_t=0.  ; rhpref_t=0. !09-11-15
       salref_u=0.  ; temref_u=0.  ; rhpref_u=0. !09-11-15
       salref_v=0.  ; temref_v=0.  ; rhpref_v=0. !09-11-15
      endif                         !>>>>
#endif

! Ultimes allocations avant lecture restart !20-11-17
      call principal_allocate_last !20-11-17

! depart d'un etat issu d'une precedente simulation. Obligatoirement
! derniere action d'initial_main.F
      iteration3d_restart=kount0                                                !18/04/08
      elapsedtime_rst=elapsedtime_now                                  !16-04-11
      if(initial.eq.1)then  ! hothothot >                              !28/05/04
        call dyn_restart('r') ! départ d'un état précédent  !23-06-10
        iteration3d_restart=iteration3d
        elapsedtime_rst=elapsedtime_now                                !16-04-11
!       if(flag_nh2d==1)iteration2d_begin=iteration2d                  !03-10-17
      endif                 ! hothothot >

      if(initial.eq.2)then                                             !16/06/06
      if(imodelbio.eq.1.or.imodeltrc.eq.1)call bio_restart('r') !23-06-10
!        call InitNutDoxy
      endif

      initial_main_status=1 !23-02-12
      kstop=0
      kpvwave=0
      give_chanel9=0
      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
       open(unit=3,file=trim(tmpdirname)//'interrupteur')     !26-08-15
       write(3,'(i1,9x,i1,9x,i1)')kstop,kpvwave,give_chanel9            !26/04/05
       close(3)
      endif                !#mpi-->>-->                       !09-05-10

! Sauver une copie de ioffline apres lecture du fichier restart:
      ioffline_prv=ioffline !19-01-14

      call initial_main_inconsistencies2              !13-09-13

#ifdef checkmpi
      if(flag_1dv==1) then
         call check_mpi_1dv_obc !27-11-15
         call check_mpi_1dv_ts('I')
      endif
#endif

      end subroutine initial_main

!............................................................................

      subroutine initial_main_inconsistencies1 !02-03-15
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='initial_main_inconsistencies1'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


! Order of advection scheme combined to mpi deepness require iglb>=4 and jglb>=4
       if(iglb<4) &
       stop ' ERROR in specifying Oi dimension: iglb must be >=4'
       if(jglb<4) &
       stop ' ERROR in specifying Oj dimension: jglb must be >=4'

       if(iglb==4.and.jglb==4) then ! 1DV case !26-11-15

        flag_1dv=1 !27-11-15

        if(nbdom_imax/=1.or.nbdom_jmax/=1) then
         stop 'nbdom_imax/=1.or.nbdom_jmax/=1 & iglb==4.and.jglb==4' 
        endif

        if(.not.iperiodicboundary) & !23-08-16
        stop 'notebook_grid: IF 1DV CASE SET iperiodicboundary=.true.'
        if(.not.jperiodicboundary) & !23-08-16
        stop 'notebook_grid: IF 1DV CASE SET jperiodicboundary=.true.'

        if(.not.fplan1_grid) & !23-08-16!03-10-18
          stop 'notebook_grid: IF 1DV CASE SET fplan1_grid=.true.'
        if(fplan2_grid/=1) & !16-10-18
          stop 'notebook_grid: IF 1DV CASE SET fplan2_grid=1'

        write(6,'(a)')lonlatfile
        if(lonlatfile/='nofile') &
        stop 'notebook_grid: IF 1DV CASE SET lonlatfile=nofile'

       endif                        ! 1DV case


      end subroutine initial_main_inconsistencies1 !02-03-15

!............................................................................

      subroutine initial_main_inconsistencies2 !13-09-13
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='initial_main_inconsistencies2'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

        if(elapsedtime_now>elapsedtime_end) then !>>>>> !17-02-15
         write(10+par%rank,*)'Restart date error:' &
                             ,'elapsedtime_now>elapsedtime_end'
         write(10+par%rank,*)'elapsedtime_now=',elapsedtime_now
         write(10+par%rank,*)'elapsedtime_end=',elapsedtime_end
         call elapsedtime2date(elapsedtime_now,i1,i2,i3,i4,i5,i6)
         write(10+par%rank,*)'elapsedtime_now corresponding day:'
         write(10+par%rank,*)i1,i2,i3,i4,i5,i6
         call elapsedtime2date(elapsedtime_end,i1,i2,i3,i4,i5,i6)
         write(10+par%rank,*)'elapsedtime_end corresponding day:'
         write(10+par%rank,*)i1,i2,i3,i4,i5,i6
         call elapsedtime2date(0.,i1,i2,i3,i4,i5,i6)
         write(10+par%rank,*)'t0 time corresponding day:'
         write(10+par%rank,*)i1,i2,i3,i4,i5,i6

         stop'Start (or restart) time > elapsedtime_end'
        endif                                    !>>>>>

      do k=1,6
       if(departuredate(k)/=datesim(k,1)) then !------> !17-02-15
        write(6,*)'Do not change the departure date in case ' &
                 ,'of restart file!'
        stop 'initial_main_inconsistencies departuredate/=datesim(1)'
       endif                                   !------>
      enddo

      if(tideforces==4.and.removetide/=0) then !44444444444444>
        write(6,*)'In notebook_tide and notebook_obcforcing'
        write(6,*)'tideforces==4.and.removetide/=0', &
                  ' are inconsistent choices'
        stop ' STOP in initial_main_inconsistencies'
      endif                                    !44444444444444>


      if(obc2dtype<1.) then !ooooooooooooooooooooooo> !13-03-13
       if(kmaxtide>1) then   !11111111111>
        write(6,*)'In notebook_tide:'
        write(6,*)'kmaxtide>1 obc2dtype<1 are inconsistent choices'
!       stop ' STOP in initial_main_inconsistencies'
       endif                 !11111111111>
       if(iobc_ogcm/=0) then !22222222222>
        write(6,*)'In notebook_obcforcing:'
        write(6,*)'iobc_ogcm/=0 obc2dtype<1 are inconsistent choices'
        stop ' STOP in initial_main_inconsistencies'
       endif                 !22222222222>
       if(iairsea/=0)   then !33333333333>
        write(6,*)'In notebook_airseaflux:'
        write(6,*)'iairsea/=0 obc2dtype<1 are inconsistent choices'
        stop ' STOP in initial_main_inconsistencies'
       endif                 !33333333333>
      endif                 !ooooooooooooooooooooooo>

! detecteur de bugs debut:
!     if(kount.lt.0) then
      if(iteration3d<0) then
      write(6,*)
      write(6,*)'.....................................................'
      write(6,*)'la valeur initiale de kount, ',iteration3d  &
               ,', est négative.'
      write(6,*)'certains calculs risquent d''etre faux comme par'
      write(6,*)'exemple le calcul du modulo pour obtenir rap'
      write(6,*)'dans aladin.f. corriger kount initial puis recommencer'
      stop ' dans symphonie.f kount négatif'
      endif
! detecteur de bugs fin.

!___________________________________________________
! Blocage en cas de demande simultanee de filiere                      !02/10/07
! offline et nesting
      if(nest_onoff_out.ge.1.and.ioffline.eq.1) then
      write(6,*)
      write(6,*)'...................................'
      write(6,*)'demande simultanee filiere offline mode ecriture et'
      write(6,*)'filiere nesting mode ecriture pour fille impossible'
      write(6,*)'choisir l''un ou l''autre dans les fichiers'
      write(6,*)'notebooline et noteboo_nesting.'
      write(6,*)'en attendant...'
      stop ' dans initial_main.f'
      endif

       if(flag_merged_levels==1) then !pmxpmx>
        
!        if(iturbulence==1) then !ooo> ! commentE le 10-07-17
!         write(6,*)'POUR LE MOMENT flag_merged_levels=1' &
!                  ,' ne peut pas etre utilisE avec k-epsilon'
!         stop ' initial_main_inconsistencies2 ' ; call mpi_finalize !04-06-17
!        endif                   !ooo>

         if(kmol_h/=0.or.kmol_s/=0.) then !aaa>
          write(6,*)'POUR LE MOMENT flag_merged_levels=1' &
                   ,' ne peut pas etre utilisE avec'      &
                   , 'kmol_h/=0.or.kmol_s/=0.'
          stop ' initial_main_inconsistencies2 ' ; call mpi_finalize !04-06-17
         endif                            !aaa>

       endif                          !pmxpmx>

      if(bulk_scheme/=bulk_core.and.       & !20-11-17
         bulk_scheme/=bulk_coare.and.      &
         bulk_scheme/=bulk_moon)      stop &
         'Err560 initial_main_inconsistencies2 error on bulk_scheme'

      if(iteration2d_max_now/=0.and.flag_timesplitting==0) then !>>>
      write(6,*)'iteration2d_max_now,flag_timesplitting',iteration2d_max_now,flag_timesplitting
      stop 'Err 260 iteration2d_max_now/=0.and.flag_timesplitting==0'
      endif                                                     !>>>

      if(flag_timesplitting==0.and.rhp_zavr_xy==1) then !29-08-18
       write(6,*)' No z-averaged density if NO time-splitting '
       stop 'Err 583 flag_timesplitting==0.and.rhp_zavr_xy==1'
      endif

      end subroutine initial_main_inconsistencies2

!.......................................................................

      subroutine initial_main_checktmpdir
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='initial_main_checktmpdir'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0) then !--------------------->
      open(unit=3,file=trim(tmpdirname)//'check_tmp_is_empty')
      read(3,'(a)',end=46)texte90
   46 close(3)
      if(texte90=='tmp is not empty') then
      stop 'STOP tmp directory is not empty'
      endif
      open(unit=3,file=trim(tmpdirname)//'check_tmp_is_empty')
      write(3,'(a)')'tmp is not empty'
      close(3)
      endif                !--------------------->

      end subroutine initial_main_checktmpdir

!............................................................................

      subroutine initial_main_inconsistencies3 !19-04-15
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='initial_main_inconsistencies3'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(ioffline==2) then !*****************>
       if(initialgridfile_txt=='none') then !------> !05-05-14
        ksecu=0
        if(vert_axis_conv_direc/='up')ksecu=1
        if(vert_axis_conv_start/='w')ksecu=1
        if(vert_axis_conv_end/='w')ksecu=1
        if(hori_axis_conv_start/='t')ksecu=1
        if(hori_axis_conv_end/='t')ksecu=1
        if(ksecu==1) then !>>>>>o
         write(6,'(a,a,a)')           &
        ' Stop notebook_grid inconsistent with choice ' // &
        'initialgridfile_txt=none. If "read offline" files provided by ' // &
        'S26, use: '
         write(6,'(a)')'vert_axis_conv_direc=up'
         write(6,'(a)')'vert_axis_conv_start=w'
         write(6,'(a)')'vert_axis_conv_end=w'
         write(6,'(a)')'hori_axis_conv_start=t'
         write(6,'(a)')'hori_axis_conv_end=t'
         stop ' Stop in routine initial_main.F90'
        endif             !>>>>>o
       endif                                !------>
      endif                !*****************>

      if(initialgridfile_txt=='none'.and.flag_nemoffline==1) then !----> !05-05-14
       write(6,'(a,a)')' initialgridfile_txt=none in notebook_grid and' &
      ,' notebook_offline inconsistent choice. NEMO or S case?'
       stop ' Stop in routine initial_main.F90'
      endif                                                       !---->

      if(initialgridfile_txt/='none'.and.flag_nemoffline==0) then !----> !05-05-14
       write(6,'(a,a)')' initialgridfile_txt/=none in notebook_grid and' &
      ,' notebook_offline inconsistent choice. NEMO or S case?'
       stop ' Stop in routine initial_main.F90'
      endif                                                       !---->

      if(iobc_ogcm==1.and.ioffline==2) then !08-01-15
       write(6,*)'--------------------------------'
       write(6,*)'iobc_ogcm=1 & ioffline=2 is not logical.' &
        ,' Cancel the use of external fields in notebook_obcforcing.'
       stop 'Stop in subroutine initial_main_inconsistencies'
      endif

      if(lonlatfile/='nofile'.and.initialgridfile_txt/='none') then
       write(6,'(a,a)')'lonlatfile/=nofile & initialgridfile_txt/=none'&
      ,'is an inconsistent choicer: NEMO or S case?'
       stop 'Stop in subroutine initial_main_inconsistencies'
      endif

      end subroutine initial_main_inconsistencies3

!.......................................................................

      subroutine initial_main_inconsistencies4(rmax_)
      use module_principal ; use module_parallele
      implicit none
      real rmax_
#ifdef synopsis
       subroutinetitle='initial_main_inconsistencies3'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


      if(rmax_>0.3.and.ihybsig==1) then !>>>>>>>

       if(par%rank==0)                                          &
       write(6,*)' WARNING (and stop)!'                         &
       ,' ihybsig==1 builds the VST grid based on the respect'  &
       ,' of the rmax criteria which seems too big (rmax>0.3).' &
       ,' Fix notebook_vertcoord or notebook_bathy.'


#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
       stop ' STOP initial_main_inconsistencies4'

      endif                             !>>>>>>>


      end subroutine initial_main_inconsistencies4

!...............................................................................

      subroutine initial_main_ww3gridfile !10-06-15
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='initial_main_ww3gridfile'
       subroutinedescription=                                        &
        ' Produces the input grid files for the ww3 model using the' &
       ,' same curvilinear grid as S model'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif



! Pour le moment cette routine ne fonctionne pas en mode parallele
      if(nbdom==1) then !--------->

      open(unit=10,file=trim(tmpdirname)//'lon_input_ww3.txt')
      open(unit=11,file=trim(tmpdirname)//'lat_input_ww3.txt')
      open(unit=12,file=trim(tmpdirname)//'h_input_ww3.txt')
      do j=1,jmax
       write(10,'(9999(f10.5))')(lon_t(i,j)*rad2deg         ,i=1,imax) !08-07-15
       write(11,'(9999(f10.5))')(lat_t(i,j)*rad2deg         ,i=1,imax)
       write(12,'(9999(f10.5))')(  h_w(i,j)*mask_t(i,j,kmax),i=1,imax)
      enddo
      close(10)
      close(11)
      close(12)

      open(unit=3,file=trim(tmpdirname)//'obc_input_ww3.txt')
! OBC i=1
      i=1
      j0=0
      do j=1,jmax,5
       if(mask_t(i,j,kmax)==1)then
         j0=j0+1
         write(texte30,'(a3,i0)')'obc',j0
         write(3,'(f10.5,1x,f10.5,1x,3a)')lon_t(i,j)*rad2deg,lat_t(i,j)*rad2deg,'"',trim(texte30),'"'
       endif
      enddo

! OBC i=imax
      i=imax
      do j=1,jmax,5
       if(mask_t(i,j,kmax)==1)then
         j0=j0+1
         write(texte30,'(a3,i0)')'obc',j0
         write(3,'(f10.5,1x,f10.5,1x,3a)')lon_t(i,j)*rad2deg,lat_t(i,j)*rad2deg,'"',trim(texte30),'"'
       endif
      enddo

! OBC j=1
      j=1
      do i=1,imax,5
       if(mask_t(i,j,kmax)==1)then
         j0=j0+1
         write(texte30,'(a3,i0)')'obc',j0
         write(3,'(f10.5,1x,f10.5,1x,3a)')lon_t(i,j)*rad2deg,lat_t(i,j)*rad2deg,'"',trim(texte30),'"'
       endif
      enddo

! OBC j=jmax
      j=jmax
      do i=1,imax,5
       if(mask_t(i,j,kmax)==1)then
         j0=j0+1
         write(texte30,'(a3,i0)')'obc',j0
         write(3,'(f10.5,1x,f10.5,1x,3a)')lon_t(i,j)*rad2deg,lat_t(i,j)*rad2deg,'"',trim(texte30),'"'
       endif
      enddo

      close(3)

      endif !---------------------->

      end subroutine initial_main_ww3gridfile
!............................................................
      subroutine initial_main_bin2netcdf
      use module_principal
      use module_parallele
      use pnetcdf
      implicit none
      character*2 txt2
      integer :: compteur=0 , time_=0                      &
        ,imin_=9999 ,imax_=-9999 ,jmin_=9999 ,jmax_=-9999  &
        ,loop_=0,ishift_=0,jshift_=0,flag_netcdf=0         &
        ,looptime_,looptimeend_=1 ,imaxm2        &
        ,time_selection_=0,flag_grid_=0           &

! ICI ON CHOISIT SI ON SORT TOUTES LES ECHEANCES (flag_all_or_one_=0) 
! OU UNE ECHEANCE PARTICULIERE (flag_all_or_one_= numero de l'echeance)
        ,flag_all_or_one_=4

      real, dimension(:,:)   , allocatable :: anyvar2dbis
      real, dimension(:,:,:) , allocatable :: anyvar3dbis

      if(nbdom_jmax/=1) then !>>>
       stop 'nbdom_jmax/=1 seul nbdom_imax peut etre > 1'
      endif                  !>>>

! Pour pnetcdf il faut que toutes les coeurs soient au travail. On
! verifie donc que le nombre de fichier dans la liste est un multiple
! du nombre de coeurs
       flag_stop=0
       k=0
       open(unit=10,file='liste_in')
 801    read(10,'(a)',end=800)texte90
        k=k+1
        goto 801
 800   close(10)
       if(mod(k,nbdom)/=0)flag_stop=1
       if(flag_stop/=0) then !debug>
        write(6,*)' Le nombre de fichiers n''est ' &
                 ,'pas un multiple du nombre de coeurs'
       stop 'Err 809 initial_main_bin2netcdf'
       endif                 !debug>

      dim_varid=100
      filval=-32767.

      allocate(texte80(11))
      allocate(ub2                 (2)) ; ub2=0
      allocate(lb2                 (2)) ; lb2=0
      allocate(ub3                 (3)) ; ub3=0
      allocate(lb3                 (3)) ; lb3=0
      allocate(vardim              (4)) ; vardim=0
      allocate(varstart            (4)) ; varstart=0
      allocate(varcount            (4)) ; varcount=0
      allocate(cgridshift          (3)) ; cgridshift=0
      allocate(varid(dim_varid))        ; varid=0
      texte80(7)='real'
      texte80(11)='none'


      do loop_=1,4
      if(par%rank==0) then !>>>>
        if(loop_==1)write(6,*)'-------------------------loop_=',loop_ &
                      ,' Trouver min max en i et j & nombre d echeances'
        if(loop_==2)write(6,*)'-------------------------loop_=',loop_ &
                             ,' Construction du header'
        if(loop_==3)write(6,*)'-------------------------loop_=',loop_ &
                             ,' Ecrire filval'
        if(loop_==4)write(6,*)'-------------------------loop_=',loop_ &
                             ,' Ecrire la donnee'
      endif                !>>>>

! On ne peut pas faire tous les fichiers en meme temps (ca rame...)
! on les fait donc les uns apres les autres, soit la boucle looptime_
! A chaque fois on fait le fichier qui correspond A looptime_
      do looptime_=1,looptimeend_ ! looptimeend_=1 pour loop_=1 
                                  ! si loop_>1 et flag_all_or_one_==0 looptimeend_=max(time_)
                                  ! si loop_>1 et flag_all_or_one_/=0 looptimeend_=1

      if(loop_>1) then !pmx>
         if(flag_all_or_one_==0)time_selection_=looptime_
         if(flag_all_or_one_/=0)time_selection_=flag_all_or_one_
      endif            !pmx>

       if(loop_==1)loop_netcdf=-999
       if(loop_==2)loop_netcdf=0
       if(loop_==3)loop_netcdf=1
       if(loop_==4)loop_netcdf=1
       compteur=0


!     elapsedtime_bef=elapsedtime_now
!     elapsedtime_bef=-99999.

      open(unit=10,file='liste_in')

  825 continue
      elapsedtime_bef=-99999.

      do k=0,nbdom-1
      compteur=compteur+1
      if(loop_==2.and.compteur>nbdom)goto 904 ! la construction du header est terminee pour compteur>nbdom
      if(loop_==3.and.compteur>nbdom)goto 904 ! le reset par filval est termine pour compteur>nbdom
       read(10,'(a)',end=904)texte90
       if(index(texte90,'/grid_')/=0)flag_grid_=1
       if(k==par%rank)offline_binrec_name=trim(texte90)
      enddo
! A chaque ouverture de fichier le compteur de variable est remis A zero:
!     count_netcdfvar=0

      if(par%rank==0)write(6,'(a)')trim(offline_binrec_name)

      if(allocated(anyvar2d)) then
        deallocate(anyvar2d)
        deallocate(anyvar2dbis)
      endif
      if(allocated(anyvar3d)) then
        deallocate(anyvar3d)
        deallocate(anyvar3dbis)
        deallocate(anyv3d)
      endif

      time_=0
      open(unit=8,file=trim(offline_binrec_name),access='sequential',form='unformatted')
      read(8)iglb,jglb,imax,jmax,kmax

      allocate(anyvar3dbis(-1:imax+2,-1:jmax+2,0:kmax+1)) ; anyvar3dbis=0
      allocate(anyvar2dbis(-1:imax+2,-1:jmax+2))          ; anyvar2dbis=0

      if(loop_>1) then !resize>
        iglb=imax_-ishift_
        jglb=jmax_-jshift_
      endif            !resize

      if(loop_==3) then !reset file with filval>
! On ecrit nbdom "tranches verticales" de filval
        jmax=jglb
        imaxm2=(iglb-2)/nbdom
        if(par%rank+1==nbdom)imaxm2=imaxm2+iglb-(imaxm2*nbdom+2)
        imax=imaxm2+2
      endif             !reset file with filval> 

      allocate(anyv3d  (-1:imax+2,-1:jmax+2,0:kmax+1,any1:any2)) ; anyv3d=0
      allocate(anyvar3d(-1:imax+2,-1:jmax+2,0:kmax+1)) ; anyvar3d=0
      allocate(anyvar2d(-1:imax+2,-1:jmax+2))        ; anyvar2d=0

  826 continue

        read(8,end=903)texte80(1)      ! character*200 nom de la variable
        read(8)texte80(2)              ! character*200 unites
        read(8)texte80(3)              ! character*200 long name
        read(8)texte80(4)              ! character*200 standard name
        read(8)texte80(5)              ! character*200 TZYX ou TYX
        read(8)txt2                    ! character*2   _t _u _v _w
        read(8)i0    ! integer       i local --> i global
        read(8)j0    ! integer       j local --> j global
        read(8)elapsedtime_now ! real*8        temps ecoule en secondes
        read(8)year_now,month_now,day_now,hour_now,minute_now,second_now ! 6 integer pour la date
  

        if(loop_==1) then !1111111>
! A la fin i=imin_ sera i=1 et i=imax_ sera i=iglb
         imin_=min(imin_,1+i0) !i0 est par%timax(1)
         jmin_=min(jmin_,1+j0) !j0 est par%tjmax(1)
         imax_=max(imax_,imax+i0)
         jmax_=max(jmax_,jmax+j0)
        endif             !1111111>
        par%timax(1)=i0-ishift_
        par%tjmax(1)=j0-jshift_

      if(loop_==3) then !reset file with filval>
! On ecrit nbdom "tranches verticales" de filval
! Attention dans la ligne suivante les parentheses sont importantes car par%rank*(iglb-2)/nbdom ne donnerait pas
! le meme resultat que par%rank*((iglb-2)/nbdom). Il est important que le block (iglb-2)/nbdom soit calculE independement
        par%timax(1)=par%rank*((iglb-2)/nbdom)
        par%tjmax(1)=0
      endif             !reset file with filval> 

!--------!
        if(elapsedtime_now/=elapsedtime_bef) then !pmx>
        count_netcdfvar=0
! On ne passe ici seulement quand time change (quelque soit le nombre de champs par echeance)
! c'est donc l'endroit pour faire les operations d'ouverture, de definition du header
            time_=time_+1

!           if(loop_==1.and.time_==2)goto 825 ! Quand loop_=1 on veut seulement savoir l'emprise spatiale max, une seule echeance suffit donc.

! Si fichier A chaque echeance fermer le fichier precedent avant de commencer le suivant:
!           if(loop_netcdf==0.and.time_selection_==0.and.time_>1) then
!                                                                status=nfmpi_enddef(ncid)
!                                                                status=nfmpi_close(ncid)
!           endif
!           if(loop_netcdf==1.and.time_selection_==0.and.time_>1)status=nfmpi_close(ncid)

            elapsedtime_bef=elapsedtime_now
! decider quelle echeance on ecrit:
!           if(time_==time_selection_.or.time_selection_==0) then !--->
            if(time_==time_selection_) then !--->
! soit le temps correspond A l'echeance choisie soit on veut convertir toutes les echeances disponibles
              flag_netcdf=1
              write(texte90(1:8),'(i8)')day_now+100*month_now+10000*year_now
              write(texte90(9:15),'(i7)')1000000+second_now+100*minute_now+10000*hour_now
              write(texte90(9:9),'(a1)')'_'
              texte250=texte90(1:15)//'.nc'
            else                            !--->
              flag_netcdf=0
            endif                           !--->
            if(flag_grid_==1) then ! grid-case >
                 texte250='grid.nc' 
                 flag_netcdf=1
            endif                  ! grid-case >

            if(loop_netcdf==0.and.compteur<=nbdom.and.flag_netcdf==1) then !header>
                 status=nfmpi_create(par%comm2d,texte250,nf_clobber+ NF_64BIT_DATA, MPI_INFO_NULL,ncid) !05-07-18
                 if(status/=0)stop ' Err 983 nfmpi_create'
                 call netcdf_dim
            endif                                                          !header>
            if(loop_netcdf==1.and.flag_netcdf==1) then !data>
                 status=nfmpi_open(par%comm2d,texte250,nf_write+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid)
                 if(status/=0)stop ' Err 989 nfmpi_open'
            endif                                      !data>

            if(loop_==1.and.par%rank==0.and.compteur==nbdom) &
            write(6,'(10(i0,1x))')time_   &
                                  ,year_now,month_now,day_now  &
                                  ,hour_now,minute_now,second_now &
                                  ,looptime_,compteur

        endif                                     !pmx>

!--------!

        if(texte80(5)=='TYX'.or. &
           texte80(5)=='YX')  then  !2d2d>
             read(8)lb2 ! 2 integer
             read(8)ub2 ! 2 intger
             if(loop_==3) then  !oooo>
                read(8)anyvar2dbis  
                anyvar2d=filval    !reset file with filval
             else               !oooo>
                read(8)anyvar2d
             endif              !oooo>
        endif                       !2d2d
        if(texte80(5)=='TZYX'.or. &
           texte80(5)=='ZYX' ) then !3d3d>
             read(8)lb3 ! 3 integer
             read(8)ub3 ! 3 intger
             if(loop_==3) then  !oooo>
                read(8)anyvar3dbis  
                anyvar3d=filval    !reset file with filval
             else               !oooo>
                read(8)anyvar3d
             endif              !oooo>
        endif 


     if(flag_netcdf==1 ) then !>>>>> 

      if(loop_netcdf==0.and.compteur<=nbdom)call netcdf_main(txt2)  
      if(loop_netcdf==1)                    call netcdf_main(txt2)

      endif                                             !vvvv>

      goto 826
  903 close(8)

      goto 825


  904 continue
      close(10)

      if(loop_==1) then  !>>>

       if(flag_all_or_one_==0) then
          looptimeend_=time_ ! faire N passages si N echeances et si flag_all_or_one_=0 (faire toutes les echeances)
       else
          looptimeend_=1     ! faire 1 passage si flag_all_or_one_/=0 (Ecrire une seule echeance)
       endif

       call mpi_allreduce(imin_,i1,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       call mpi_allreduce(jmin_,j1,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       call mpi_allreduce(imax_,i2,1,mpi_integer,mpi_max,par%comm2d ,ierr)
       call mpi_allreduce(jmax_,j2,1,mpi_integer,mpi_max,par%comm2d ,ierr)
       imin_=i1
       jmin_=j1
       imax_=i2
       jmax_=j2

! A la fin i=imin_ sera i=1 et i=imax_ sera i=iglb
       ishift_=imin_-1
       jshift_=jmin_-1
       if(par%rank==0) then !oooo>
           write(6,*)'imin et max',i1,i2
           write(6,*)'jmin et max',j1,j2
           write(6,*)'imin et max shifte',imin_-ishift_,imax_-ishift_
           write(6,*)'jmin et max shifte',jmin_-jshift_,jmax_-jshift_
       endif                !oooo>


      endif              !>>>

       if(loop_netcdf==0) then
             status=nfmpi_enddef(ncid)
          if(status/=0)stop ' Err 1074 nfmpi_enddef'
             status=nfmpi_close(ncid)
          if(status/=0)stop ' Err 1076 nfmpi_close'
       endif
       if(loop_netcdf==1)status=nfmpi_close(ncid)
                      if(status/=0)stop ' Err 1079 nfmpi_close'

      enddo  ! boucle looptime_
      enddo  ! boucle loop_

      stop 'SUCCESSFULL CONVERSION'
      end subroutine initial_main_bin2netcdf
!............................................................
