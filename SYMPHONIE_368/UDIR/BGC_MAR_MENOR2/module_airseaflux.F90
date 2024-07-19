      module module_airseaflux
!______________________________________________________________________
! SYMPHONIE ocean model
! release S26.1 - last update: 28-10-18
!______________________________________________________________________
!  _________                    .__                  .__         !m°v°m
! /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____   !
! \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \  !
! /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/  !
!/_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  > !
!        \/\/          \/|__|        \/            \/        \/  !
!.................................................................
! Version date      Description des modifications
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         19/03/02: amenagement etat initial similaire e obc_af.F et model_in.F
!         11/07/02: amenagements pour filiere formules bulk. Le 3eme argument
!                   des tableaux (0 1 ou 2) n'est plus code en dur mais
!                   parametre avec I0 et I2 pour eviter les warning e la
!                   compilation (avec formules bulk le 3eme argument est
!                   reduit e 1)
!         29/01/03: securite pour eviter divisions par zero. PSS_Z est bornee.
!         19/03/03: ajout des precipitations
!         06/01/04: ajout option calcul du longwave flux a partir des parametres
!                   meteo et preciser option form='unformatted'
!         06/05/04: ajout d'une option permettant de commencer avant nc=1
!         01/11/04: Calcul de la pression atmospherique moyenne dans le but
!                   d'introduire une condition de barometre inverse sur obc de la
!                   sse.
!         07/03/05: Ecriture e l'ecran pour verifier que la lecture des fichiers
!                   se passe bien.
!         16/06/06: rationnalisation du calcul en mode offline
!         26/03/07: on prefere (pour le moment) que la pression atmospherique de
!                   reference (pour calcul du BI) soit une constante plutet qu'une
!                   moyenne sur le domaine numerique...
!         22/01/08: Le renouvellement des echeance via le modulo d'entier est
!                   supprime.
!         02/02/08: modif affichage ecran
!         03/12/08: introduction de l'interpolation online des fichiers meteo
! 2009.3  05-11-09  ecrire les fichiers temporaires dans repertoire tmp
! 2010.3  08-01-10  la rotation de grille est calcule avec les tableaux gridrotcos_z
!                   et gridrotsin_z qui remplacent desormais cos et sin de angle0
! 2010.6  02-02-10  renomme lon_t lat_t
!         04-02-10  Ajout d'une filiere de lecture de champs au format netcdf
!         10-02-10  Modif calcul bouchage des trous pour conservation de la
!                   parallelisation
! 2010.7  14-02-10  deplacement affectation valeur sur i0 i2
! 2010.8  11-03-10  une rampe sur le vent
!         19-03-10  rappel T de surface
! 2010.11 01-07-10  debug sur detecteur de depassement de grille
!         16-07-10  temobc & salobc renommes temobc & salobc
! 2010.12 26-08-10  adaptation aux fichiers ecmwf
! 2010.13 13-10-10  Ecriture reservee e proc 0
!         03-11-10  des arguments passes dans date_to_kount
! 2010.18 24-03-11  Version simulation "cascade" par Claude
!         26-03-11  lecture time_ compatible avec f95
! 2010.19 14-04-11  fermeture fichier netcdf
! 2010.20 19-04-11  Calculs sur la base d'un temps en secondes
! 2010.22 30-04-11  poursuite du point precedent
! 2010.24 14-11-11  suppression de la rampe sur le vent
! 2010.25 19-04-12  debug aiguillage dewpoint
! S25.4   02-07-12  mises a jour cas ecmwf et ajout cas arome
! S26     18-09-12  Possibilite de mask terre/mer "reel" pour le modele meteo
!         08-02-13  test etoffe sur Units pour decumul arome
!         07-03-13  ajout echange mpi
!         12-03-13  rustine pour rattraper bug de glorys
!         26-06-13  debug ecriture des listes binaires
!         02-01-14  affichages a l'ecran
!         03-01-14  debug cumul
!         07-02-14  debug cas glorys
!         13-03-14  debug cas glorys
!         08-04-14  acceleration liste binaire
!         09-04-14  ajout airseaflux_inquire_albedo
!                      et airseaflux_sphum_or_dewp
!         20-04-14  suite du point precedent
!         08-05-14  reduire le nombre d'actions dans le cas ioffline=2
!         15-05-14  compatibilite procedure "glorys" aux champs meteo du projet simed
!         11-07-14  rap_obc devient timeweightobc
!         30-07-14  Pour Caroline: prodedure automatisee depuis notebook_airseaflux
!                   de l'option moyenne 24h du flux solaire
!         31-07-14  Nouveaux echanges
!         10-08-14  Tous les procs font les listes binaires
!         04-10-14  - ssr24prv_w=0.
!                   - echange za pour le wind stress
!         27-10-14  meteo_t0 pour reset des variables cumulees
!         08-01-15  si ioffline=2 on considere maintenant pss pour
!                   convenance diagnostique de verification de conservation de ssh
!                   stop si la grille ocean n'est pas incluse a 100% dans la grille
!                   ecmwf
!         24-01-15  subroutine airseaflux_driver(case_) reorganisee
!         03-02-15  Bulk cas simu 2DH
!         27-02-15  Detection de listes desordonnees
!         10-03-15  ajout variable rsds
!         11-03-15  message ecran
!         19-03-15  detection de noms de listes invalides - correction d'un bug introduit
!                   lors d'une recente mise a jour
!         21-03-15  conservation mpi option ecmwf
!         23-03-15  - correction d'une boucle
!                   - nf_nowrite remplace par nf_share (Cyril)
!         24-03-15  ncid_ remplace ncid1
!         09-04-15  evaporation/precipition = terme puit/source ssh
!         23-04-15  echange za sur tableaux ij2meteo_i; ij2meteo_j
!         25-05-15  module_optics remplace module_albedo
!         06-06-15  relaxed sst deplace dans le driver
!         21-06-15  si le calcul de correspondance de grilles ne converge pas car
!                   la grille meteo est discontinue, un algo alternatif est possible
!         01-07-15  flag_p0m_filter permet d'activer le filtrage du bruit de la pression
!                   atmospherique d'ecmwf + correction d'un bug concernant ce meme point
!         12-10-15  ajout du cas 1DV
!         31-10-15  ajout wetmask_t sur omega_w(kmax+1)
!         24-11-15  ajout du vent dans la couche limite simplifiee
!         22-12-15  test mpi_proc_null
!         23-12-15  dimensions reduites pour champs de la CLA
!         14-02-16  - Possibilite temps au format real
!                   - La regle de decumul a EtE generalisee (un peu plus) pour Etre
!                   Egalement adaptEe aux fichier CLAS
!         10-03-16  cas 1DV glorys (Leo)
!         11-03-16  bug calcul moyenne ssr24h corrigE
!         31-03-16  Procedure de resynchronisation des fichiers pour s'adapter au cas test
!                   d'Alex sur le traitement de ssr_w (moyenne - non moyenne etc....)
!         12-04-16  - pas de post-traitement de la sst de ecmwf
!                   - Valeur de la pression atmospherique moyenne au niveau de la mer affinEe
!         18-04-16  Possibilite de repartir de listes binaires existantes dans tmp
!         19-04-16  Seul le rank0 teste l'existence des listes binaires dans tmp
!         19-04-16  modif mpi
!         24-04-16  - rho remplace par rhp+rho dans calculer omega_w(:,:,kmax+1)
!                   - choisir SSR plutot que SSRD si possible
!         25-04-16  verifier la coherence de la liste des fichiers meteo
!         28-04-16  - listes binaires: recl=540
!                   - ajustements pour lecture des fichiers MESO-NH
!         02-05-16  Choix des formules bulk
!         05-05-16  Mise A jour de la chaine d'interpolation si discontinuite dans la liste 
!                   des fichiers meteo (ECMWF, MESO-NH etc...)
!         07-05-16  lecture plus stricte de la liste des bouchetrous pour conservation mpi 
!                   en cas de successions ECMWF-MESO-NH
!         08-05-16  test relatif A formules bulk CORE
!         11-05-16  ajout lecture vent moyen de la ABL de meso-nh
!         15-05-16  verifier deallocation des tableaux avant allocation
!         19-05-16  Ajout d'un cas de lecture de l'humidite relative
!         21-06-16  Ajout d'une possibilite de lecture de la pression
!                   corrigee au niveau de la mer (champ nommE MSL)
!         08-11-16  ajout des champs arome
!         28-11-16  prevoir cas arome 4 dimensions
!         19-01-17  cas glorys IR et Net IR
!         20-01-17  Message Ecran
!         15-02-17  un stop commentE
!         12-03-17  Hors cas evident de la SABL, il vaut mieux ne pas utiliser le vent
!                   au dessus du continent
!         25-03-17  possibilite de parametrer l'effet sterique
!         27-03-17  ajout d'un lien http pour l'effet sterique
!         10-06-17  possibilite d'avoir une grille oceanique plus grande
!                   que la grille atmospherique
!         27-07-17  On remplace wetmask basE sur wetdry_cst3 par la formule
!                   basEe sur wetdry_cst1 (plus grand que wetdry_cst3)
!         20-11-17  Exiger MSL A la place de SP
!         25-01-18  wstress donnE par le modele meteo et formules bulk pour autres flux
!         11-04-18  - conversion humidite relative --> specificique par l'equation de Clausius-Clapeyron
!                   - reset varname         
!         05-05-18  - Procedure d'erret d'urgence pour occigen 
!                   - ajout meteo_enlarged_bounds
!                   - seul par%rank=0 lit liste binaire puis mpi_bcast
!         09-05-18  - Procedure d'erret d'urgence pour occigen 
!         14-06-18  fichier bouchetrou remplacE par tableau
!         10-07-18  verifier les depassements de memoire
!                   modifs pour passer le cas du modele meteo plus petit que le modele ocean
!         13-07-18  flag_meteo_land_plug permet de choisir si on prend les valeurs en terre ou 
!                   si on les remplace par la valeur en mer la plus proche
!         28-10-18  ajout meteo_cumul_modulo
!...............................................................................
      use module_principal ; use module_forcages ; use module_parallele
      use module_optics ; use module_atmboundlayer ; use module_s
      implicit none
      include 'netcdf.inc'

      real*4,dimension(:,:),allocatable :: meteo_lonlat_r4  &
                      ,ij2meteo_i,ij2meteo_j,ij2meteo_teta  &
                      ,meteo_var2

      double precision,dimension(:)  ,allocatable :: scalefct_all !26-01-22
      integer,dimension(:)  ,allocatable :: varcode &
                                           ,flag_cumul_all !26-01-22
      integer,dimension(:,:),allocatable :: meteoplugs !14-06-18

      integer :: flag_meteotime,flag_cumul,var_dims,var_type,lon_id,lat_id &
             ,tabdim(4),flag_lsm,meteo_imax_full,meteo_jmax_full &
             ,meteolandconvention,misspointmax=0,meteo_enlarged_bounds=20 &
             ,meteo_cumul_modulo=0     & !28-10-18
             ,autocumul_status=0       & !01-07-20
             ,nc_meteo_max=0             !26-01-22


      integer(kind=1) :: flag_gridoverflow       & !10-06-17
                        ,flag_refresh_interp=0     !26-01-22


!     integer :: ssr_id=0  ,ir_id=0,rain_id=0,t2m_id=0       &
!              ,dp2m_id=0,u10m_id=0,v10m_id=0,p0m_id=0       &
!              ,ustrs_id=0,vstrs_id=0,slhf_id=0,netir_id=0   &
!              ,sshf_id=0

      double precision lat_1_1_gridt,lat_1_1_gridu,lat_1_1_gridv &
                      ,lon_1_1_gridt,lon_1_1_gridu,lon_1_1_gridv &
                      ,lon_2_1_gridt,lat_1_2_gridt,scalefct

      real c_grid_shift

      character*10 :: meteovarname(7)
      character*100 :: meteo_grid_file='none'

contains

!..............................................................................

      subroutine airseaflux_driver(case_) !24-01-15
      integer case_
#ifdef synopsis
       subroutinetitle='airseaflux_driver'
       subroutinedescription= &
          'Principal driver of the subroutines computing the surface' &
       //' air/sea fluxes'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Note: if case_==1 initial state
!       if case_==2 iterative mode

!............
! Case of simple academic wind stress (no meteo file input):
        if(iwind==1.and.iairsea==0)call windstress

!............
! Case straightforward fluxes (no bulk):
        if(iairsea==1) then ! fluxes (no bulk)
         if(flag_meteodata/='glorys')call airseaflux_upd(case_)
         if(flag_meteodata=='glorys')call airseaflux_glorys(case_)
        endif               ! fluxes (no bulk)

!............
! Case where fluxes are computed using bulk formulae:
        if(iairsea==2.or.iairsea==3) then !-=-=-=-=-=->
         if(kmax==1)then !-2Dcase-2Dcase-> ! 2D run assumes SST=T0m since SST
          if(iteration3d/=0) then !>>>>             !03-02-15
           do j=1,jmax ; do i=1,imax                ! is required to compute windstress CD
            tem_t(i,j,kmax,1)=teta2_t(i,j,1)-273.15 !03-02-15
           enddo ; enddo
          endif                   !>>>>
         endif           !-2Dcase-2Dcase->
!          if(flag_meteodata/='ecmwf')call airseaflux_fbk(case_)
           if(    flag_meteodata=='ecmwf' &
              .or.flag_meteodata=='arome') then !pmxpmx> !08-11-16
                       call airseaflux_ecmwfbk(case_)
           else                                 !pmxpmx>
            stop 'huuum are you sure you want to use airseaflux_fbk?' !08-01-23
                       call airseaflux_fbk(case_)
           endif                                !pmxpmx>
        
        endif                             !-=-=-=-=-=->

!............
! Update albedo (if required):
        if(ialbedo==1)call optics_albedo_upd !24-05-15

! Surface vertical velocity deduced from evaporation and precipitation
        call airseaflux_omega_surf !09-04-15

! irelaxsst default value is 0. SST is relaxed toward temobc(kmax) if irelaxsst=1 !06-06-15
      if(irelaxsst==1)call airseaflux_relaxed_sst !06-06-15

      end subroutine airseaflux_driver

!..............................................................................

      subroutine airseaflux_omega_surf !09-04-15
      implicit none
#ifdef synopsis
       subroutinetitle='airseaflux_omega_surf'
       subroutinedescription=                                &
        'Surface vertical velocity deduced from evaporation' &
       ' and precipitation'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(flag_steric_effect==0) then !pjppjp>

! No steric effect:
       do j=1,jmax ; do i=1,imax
        omega_w(i,j,kmaxp1,1)=-mask_t(i,j,kmax)*( & !pmxpmx>
         slhf_w(i,j,1)/(lv*(rho+rhp_t(i,j,kmax)))*max(min( hz_w(i,j,1)/wetdry_cst1,1.),0.) & 
                           +precipi_w(i,j,1) & !18-02-22
                                                )   !pmxpmx>

       enddo ; enddo

      else                           !pjppjp> !25-03-17

! Parametrization of the steric effect through a vertical surface velocity. Details in:
! link 1:
! https://docs.google.com/document/d/1x9lWebKZzgpWGAlKypjulwtNvN-YiufIKb3EKN__yk8/edit
! link 2:
! https://docs.google.com/presentation/d/1zEqG6sk01MrGwlybE2PXkm9bzBZj4Z7Q9DUtBPutW1M/edit#slide=id.g201e91634d_0_100

      do j=1,jmax ; do i=1,imax

!       omega_evaprec_w(i,j,0)=omega_evaprec_w(i,j,1) ! commentee le 06-06-19

       omega_w(i,j,kmaxp1,1)=-mask_t(i,j,kmax)*( & !m°v°m>

! Precipitations
       precipi_w(i,j,1) &  !18-02-22 

      +( &                                           !-wetdry->

! Evaporation condensation
        slhf_w(i,j,1)/(lv*(rho+rhp_t(i,j,kmax)))                & !24-04-16

! Steric effect approximation: !25-03-17
       -( (1.-albedo_w(i,j))*ssr_w(i,j,1)   &
                           +snsf_w(i,j,1)   &
                           +slhf_w(i,j,1)   &
                           +sshf_w(i,j,1)   &
        )                                   &
       *(4.9762e-2+tem_t(i,j,kmax,1)*(-2.*7.5911e-3                   & ! Thermal coef expansion derived from equation 3 in :
                                      +3.*3.5187e-5*tem_t(i,j,kmax,1) & ! Friedrich and Levitus, 1972. JPO.
                                      +2.*3.7297e-5*sal_t(i,j,kmax,1))& ! An Approximation to the Equation of State for Sea Water
                  -sal_t(i,j,kmax,1)*3.0063e-3)                       & ! suitable for numerical ocean models. 
               /(cp*(rho+rhp_t(i,j,kmax))**2)  &


       )*max(min( hz_w(i,j,1)/wetdry_cst1,1.),0.)  & !-wetdry->

                                               )   !m°v°m>
      enddo ; enddo

      endif                          !pjppjp>

! NOTE: si pour Economiser une etape de calcul on calcule
! directement omega_evaprec_w(i,j,1) dans les lignes precedentes
! attention A:
!              1- remplacer omega_w(kmax+1) par omega_evaprec_w(1) dans
!                 module_external_mode et autres routines utilisant omega_w(kmax+1)
!              2- ...
! voir aussi https://docs.google.com/document/d/1fIzAG9mo_zvtTdVauu7yw-VkIqFRldqq4tF2b6ZOnzI/edit
      do j=1,jmax ; do i=1,imax !28-02-19
        omega_evaprec_w(i,j,0)=omega_evaprec_w(i,j,1)
        omega_evaprec_w(i,j,1)=omega_w(i,j,kmaxp1,1)
      enddo ; enddo

! Ce flag indique que omega_w(i,j,kmaxp1,1) est desormais non nul et que 
! les autres operations sur omega de surface (par exemple les rivieres ajoutees en surface)
! doivent incrementer omega_w sans perdre la valeur acquise dans la presente subroutine
      flag_omega_cumul=1 !18-02-21

      end subroutine airseaflux_omega_surf

!..............................................................................

      subroutine airseaflux_ecmwfbk(case_)
      implicit none
      integer case_,loop_,t_,tstr_,flag_rotate_,flag_phys1_    &
             ,flag_phys2_,flag_bulk_,flag_bulk_sstmeteo_
#ifdef synopsis
       subroutinetitle='airseaflux_ecmwfbk'
       subroutinedescription='ECMWF general driver'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


      flag_bulk_sstmeteo_=0
      flag_rotate_=0
      flag_phys1_=0
      flag_phys2_=0
      tstr_=2
      flag_bulk_=2



!...................................................!
! INITIAL PHASE:
      if(case_==1) then !- Initial phase -> !15-03-13
       tstr_=0
       flag_bulk_=1
       call airseaflux_initial
      endif             !> Etat initial >
!...................................................!




       do t_=tstr_,2,2 ! Si etat initial, tstr_=0 donc 2 passages, sinon tstr_=2 donc un seul passage

        loop_=ssr_id
        call airseaflux_decision(loop_,case_) ! donne feu vert pour lecture fichier !09-02-14
        if(decision==1) then !ooooo>
         call airseaflux_interp_driver(loop_,t_,loop_) ! > xy_t(:,:,1)
         do j=0,jmax+1 ; do i=0,imax+1
          ssr_w(i,j,0)=ssr_w(i,j,2)
          ssr_w(i,j,2)= xy_t(i,j,1)
         enddo; enddo
         if(flag_ssr24avr==1)call airseaflux_ssr24avr !30-07-14

!.....>
! Aide debug (pour par exemple aider a detecter les erreurs de decumul)
! Ici on detecte les valeurs negatives (le test ne porte pas stritement
! sur la valeur 0 pour ne pas confondre avec des erreurs de precision):
         flag_stop=0
         if(minval(ssr_w(:,:,2))<-1.) then !>>> !01-07-20
          flag_stop=1
          write(10+par%rank,*)'Erreur ssr negatif'
          write(10+par%rank,*)'minval(ssr_w(:,:,2)',minval(ssr_w(:,:,2))
          write(10+par%rank,*)'minloc(ssr_w(:,:,2)',minloc(ssr_w(:,:,2))
          write(10+par%rank,*)'indices globaux obtenus en ajoutant',par%timax(1),par%tjmax(1)
         endif                             !>>> !01-07-20
         call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
         if(k0/=0)stop 'negative SSR. See fort.xxx error files'
!.....>

        endif                !ooooo>
! commente le !04-10-14 car le vent intervient desormais dans le flux d'O2
!       if(ioffline==2)goto 185 !08-05-14

        loop_=ir_id
        call airseaflux_decision(loop_,case_) ! donne feu vert pour lecture fichier !09-02-14
        if(decision==1) then
         call airseaflux_interp_driver(loop_,t_,loop_) ! > xy_t(:,:,1)
         do j=0,jmax+1 ; do i=0,imax+1
          snsf_w(i,j,0)=snsf_w(i,j,2)
          snsf_w(i,j,2)=  xy_t(i,j,1)
         enddo; enddo
        endif

        loop_=rain_id
        call airseaflux_decision(loop_,case_) ! donne feu vert pour lecture fichier !09-02-14
        if(decision==1) then
         call airseaflux_interp_driver(loop_,t_,loop_) ! > xy_t(:,:,1)
         do j=0,jmax+1 ; do i=0,imax+1
          precipi_w(i,j,0)=precipi_w(i,j,2)
          precipi_w(i,j,2)=     xy_t(i,j,1)
         enddo; enddo
!.....> !commentE le 22-01-20 apres avoir constatE que les precip pouvaient etre negatives
! dans "nos" fichiers ecmwf "mensualisEs"
! Aide debug (pour par exemple aider a detecter les erreurs de decumul) !14-01-20
! Ici on detecte les valeurs negatives (le test ne porte pas stritement sur la valeur 0 pour
! ne pas confondre avec des erreurs de precision):
!        flag_stop=0
!        if(precipi_w(imax/2,jmax/2,2)*mask_t(imax/2,jmax/2,kmax)<-small1) then !>>> !14-01-20
!         flag_stop=1
!         write(10+par%rank,*)'Erreur precipi negatif'
!         write(10+par%rank,*)'i,j,loc',imax/2,jmax/2
!         write(10+par%rank,*)'i,j,glb',imax/2+par%timax(1),jmax/2+par%tjmax(1)
!         write(10+par%rank,*)'precipi_w',precipi_w(imax/2,jmax/2,2)
!        endif                                                             !>>> !14-01-20
!        call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
!        if(k0/=0)stop 'negative PRECIPI_W. See fort.xxx error files'
!.....>
        endif                !oooo>

        loop_=t2m_id
        call airseaflux_decision(loop_,case_) ! donne feu vert pour lecture fichier !09-02-14
        if(decision==1) then !ddddd>

          call airseaflux_interp_driver(loop_,t_,loop_) ! > xy_t(:,:,1)
          do j=0,jmax+1 ; do i=0,imax+1
           teta2_t(i,j,0)=teta2_t(i,j,2)
           teta2_t(i,j,2)=   xy_t(i,j,1)
          enddo; enddo
          flag_phys1_=1

         if(flag_abl==1) then !ablablabl>
! rembobiner la liste binaire de 1 (que le airseaflux_interp_driver precedent incrementa)
          airseafile_nextrec(loop_)=airseafile_nextrec(loop_)-1
          call airseaflux_interp_driver(t0m_id,t_,loop_) ! > xy_t(:,:,1)
          do j=0,jmax+1 ; do i=0,imax+1
           teta0_t(i,j,1)=teta0_t(i,j,2)
           teta0_t(i,j,2)=   xy_t(i,j,1)
          enddo; enddo

! rembobiner la liste binaire de 1 (que le airseaflux_interp_driver precedent incrementa)
          airseafile_nextrec(loop_)=airseafile_nextrec(loop_)-1
          call airseaflux_interp_driver(abl_id,t_,loop_) ! > xy_t(:,:,1)
          do j=0,jmax+1 ; do i=0,imax+1
           ablheight_t(i,j,1)=ablheight_t(i,j,2)
           ablheight_t(i,j,2)=       xy_t(i,j,1)
          enddo; enddo

          flag_bulk_sstmeteo_=1
         endif                !ablablabl>


        endif                 !ddddd>

        loop_=dp2m_id
        call airseaflux_decision(loop_,case_) ! donne feu vert pour lecture fichier !09-02-14
        if(decision==1) then
         call airseaflux_interp_driver(loop_,t_,loop_) ! > xy_t(:,:,1)
         do j=0,jmax+1 ; do i=0,imax+1
          q2_t(i,j,0)=q2_t(i,j,2)
          q2_t(i,j,2)=xy_t(i,j,1)
         enddo; enddo
         flag_phys2_=1
        endif
#ifdef checkmpi
      call airseaflux_checkmpi_q2
#endif

        loop_=u10m_id
        call airseaflux_decision(loop_,case_) ! donne feu vert pour lecture fichier !09-02-14
        if(decision==1) then !ddddd>

            call airseaflux_interp_driver(loop_,t_,loop_) ! > xy_t(:,:,1)
            do j=0,jmax+1 ; do i=0,imax+1
             uwind_t(i,j,0)=uwind_t(i,j,2)
             uwind_t(i,j,2)=   xy_t(i,j,1)
            enddo; enddo

         if(flag_abl==1) then !ablablabl>
! rembobiner la liste binaire de 1 (que le airseaflux_interp_driver precedent incrementa)
            airseafile_nextrec(loop_)=airseafile_nextrec(loop_)-1
            call airseaflux_interp_driver(u100m_id,t_,loop_) ! > xy_t(:,:,1)
            do j=0,jmax+1 ; do i=0,imax+1
             uwindabl_t(i,j,1,1)=uwindabl_t(i,j,1,2)
             uwindabl_t(i,j,1,2)=        xy_t(i,j,1)
            enddo; enddo
         endif                !ablablabl>

         flag_rotate_=1
        endif                !ddddd>

        loop_=v10m_id
        call airseaflux_decision(loop_,case_) ! donne feu vert pour lecture fichier !09-02-14
        if(decision==1) then !ddddd>

            call airseaflux_interp_driver(loop_,t_,loop_) ! > xy_t(:,:,1)
            do j=0,jmax+1 ; do i=0,imax+1
             vwind_t(i,j,0)=vwind_t(i,j,2)
             vwind_t(i,j,2)=   xy_t(i,j,1)
            enddo; enddo

         if(flag_abl==1) then !ablablabl>
! rembobiner la liste binaire de 1 (que le airseaflux_interp_driver precedent incrementa)
            airseafile_nextrec(loop_)=airseafile_nextrec(loop_)-1
            call airseaflux_interp_driver(v100m_id,t_,loop_) ! > xy_t(:,:,1)
            do j=0,jmax+1 ; do i=0,imax+1
             vwindabl_t(i,j,1,1)=vwindabl_t(i,j,1,2)
             vwindabl_t(i,j,1,2)=      xy_t(i,j,1)
            enddo; enddo
         endif                !ablablabl>

         flag_rotate_=1
        endif                !ddddd>

        if(flag_wstressbulk==0) then !oooooooooooooooooooo> ! wstress from meteo model !25-01-18

         loop_=ustrs_id
         call airseaflux_decision(loop_,case_) ! donne feu vert pour lecture fichier !09-02-14
         if(decision==1) then !ddddd>
            call airseaflux_interp_driver(loop_,t_,loop_) ! > xy_t(:,:,1)
            do j=0,jmax+1 ; do i=0,imax+1
             wstress_u(i,j,0)=wstress_u(i,j,2)
             wstress_u(i,j,2)=     xy_t(i,j,1) ! NOTE WSTRESS(:,:,2) SUR POINT T SANS ROTATION
            enddo; enddo
          flag_rotate_=1
         endif                !ddddd>

         loop_=vstrs_id
         call airseaflux_decision(loop_,case_) ! donne feu vert pour lecture fichier !09-02-14
         if(decision==1) then !ddddd>
            call airseaflux_interp_driver(loop_,t_,loop_) ! > xy_t(:,:,1)
            do j=0,jmax+1 ; do i=0,imax+1
             wstress_v(i,j,0)=wstress_v(i,j,2)
             wstress_v(i,j,2)=     xy_t(i,j,1) ! NOTE WSTRESS(:,:,2) SUR POINT T SANS ROTATION
            enddo; enddo
          flag_rotate_=1
         endif                !ddddd>
        endif                        !oooooooooooooooooooo>

        loop_=p0m_id
        call airseaflux_decision(loop_,case_) ! donne feu vert pour lecture fichier !09-02-14
        if(decision==1) then
         call airseaflux_interp_driver(loop_,t_,loop_) ! > xy_t(:,:,1)
         do j=0,jmax+1 ; do i=0,imax+1
          pss_w(i,j,0)=pss_w(i,j,2)
          pss_w(i,j,2)= xy_t(i,j,1)
         enddo; enddo
         pss_mean(0)=pss_mean(2) ; pss_mean(2)=101325. ! pression atmospherique moyenne au niveau de la mer !12-04-16
        endif

! Post-traitement (rotation, physique) apres que la boucle loop_ soit terminee...
! Les flags servent a ne faire le calcul qu'au moment des mises a jours des fichiers
! et non pas a chaque iteration du modele S
        if(flag_rotate_==1)call airseaflux_rotate     ! > uwind_t(:,:,2) vwind_t(:,:,2)
        if( flag_phys1_==1.and.      &
            bulk_scheme==bulk_core)  & !08-05-16
                           call airseaflux_physics(1) ! > teta2=temp pot
        if( flag_phys2_==1)call airseaflux_physics(2) ! > q2=hum. spec.

#ifdef checkmpi
      call airseaflux_checkmpi_teta2
#endif

! Si modele de couche limite atmospherique active alors calculer les flux de chaleur latente
! et sensible a partir de la SST du modele meteo. Note que ce calcul n'est fait que si
! flag_bulk_sstmeteo=1 et uniquement lorsqu'on lit de nouveaux champs (donc pas a chaque
! iteration du modele)
       if(flag_bulk_sstmeteo_==1) then
         call atmboundlayer_sstmeteo(flag_bulk_,t_)
         call atmboundlayer_windcfl
       endif

  185  continue ! goto test ioffline !08-05-14
       enddo ! t_

! Interpolation temporelle:

      loop_=ssr_id
      x2=(elapsedtime_now          -airseafile_prvtime(loop_))      &
        /(airseafile_nextime(loop_)-airseafile_prvtime(loop_))
! Bornage car l'attente de la synchro subcycling peut faire sortir (legerement) de l'intervale [O:1]
      x2=min(max(x2,0.d00),1.d00) !09-02-14
      x0=1.-x2
      do j=0,jmax+1 ; do i=0,imax+1
       ssr_w(i,j,1)=x0*ssr_w(i,j,0)+x2*ssr_w(i,j,2)
      enddo; enddo
! commente le !04-10-14 car le vent intervient desormais dans le flux d'O2
!     if(ioffline==2)return !08-05-14

      loop_=ir_id
      x2=(elapsedtime_now          -airseafile_prvtime(loop_))      &
        /(airseafile_nextime(loop_)-airseafile_prvtime(loop_))
! Bornage car l'attente de la synchro subcycling peut faire sortir (legerement) de l'intervale [O:1]
      x2=min(max(x2,0.d00),1.d00) !09-02-14
      x0=1.-x2
      do j=0,jmax+1 ; do i=0,imax+1
       snsf_w(i,j,1)=x0*snsf_w(i,j,0)+x2*snsf_w(i,j,2)
      enddo; enddo

      loop_=rain_id
      x2=(elapsedtime_now          -airseafile_prvtime(loop_))      &
        /(airseafile_nextime(loop_)-airseafile_prvtime(loop_))
! Bornage car l'attente de la synchro subcycling peut faire sortir (legerement) de l'intervale [O:1]
      x2=min(max(x2,0.d00),1.d00) !09-02-14
      x0=1.-x2
      do j=0,jmax+1 ; do i=0,imax+1
       precipi_w(i,j,1)=x0*precipi_w(i,j,0)+x2*precipi_w(i,j,2)
      enddo; enddo

      loop_=t2m_id
      x2=(elapsedtime_now          -airseafile_prvtime(loop_))      &
        /(airseafile_nextime(loop_)-airseafile_prvtime(loop_))
! Bornage car l'attente de la synchro subcycling peut faire sortir (legerement) de l'intervale [O:1]
      x2=min(max(x2,0.d00),1.d00) !09-02-14
      x0=1.-x2
      do j=0,jmax+1 ; do i=0,imax+1
       teta2_t(i,j,1)=x0*teta2_t(i,j,0)+x2*teta2_t(i,j,2)
      enddo; enddo

      loop_=dp2m_id
      x2=(elapsedtime_now          -airseafile_prvtime(loop_))      &
        /(airseafile_nextime(loop_)-airseafile_prvtime(loop_))
! Bornage car l'attente de la synchro subcycling peut faire sortir (legerement) de l'intervale [O:1]
      x2=min(max(x2,0.d00),1.d00) !09-02-14
      x0=1.-x2
      do j=0,jmax+1 ; do i=0,imax+1
       q2_t(i,j,1)=x0*q2_t(i,j,0)+x2*q2_t(i,j,2)
      enddo; enddo

      loop_=u10m_id
      x2=(elapsedtime_now          -airseafile_prvtime(loop_))      &
        /(airseafile_nextime(loop_)-airseafile_prvtime(loop_))
! Bornage car l'attente de la synchro subcycling peut faire sortir (legerement) de l'intervale [O:1]
      x2=min(max(x2,0.d00),1.d00) !09-02-14
      x0=1.-x2
      do j=0,jmax+1 ; do i=0,imax+1
       uwind_t(i,j,1)=x0*uwind_t(i,j,0)+x2*uwind_t(i,j,2)
      enddo; enddo

      loop_=v10m_id
      x2=(elapsedtime_now          -airseafile_prvtime(loop_))      &
        /(airseafile_nextime(loop_)-airseafile_prvtime(loop_))
! Bornage car l'attente de la synchro subcycling peut faire sortir (legerement) de l'intervale [O:1]
      x2=min(max(x2,0.d00),1.d00) !09-02-14
      x0=1.-x2
      do j=0,jmax+1 ; do i=0,imax+1
       vwind_t(i,j,1)=x0*vwind_t(i,j,0)+x2*vwind_t(i,j,2)
      enddo; enddo

! Si wstress n'est pas donnE par les formules bulk:
      if(flag_wstressbulk==0) then !pmx> !25-01-18

       loop_=ustrs_id 
       x2=(elapsedtime_now          -airseafile_prvtime(loop_))      &
         /(airseafile_nextime(loop_)-airseafile_prvtime(loop_))
! Bornage car l'attente de la synchro subcycling peut faire sortir (legerement) de l'intervale [O:1]
       x2=min(max(x2,0.d00),1.d00) !09-02-14
       x0=1.-x2
       do j=0,jmax+1 ; do i=0,imax+1
        xy_t(i,j,1)=x0*wstress_u(i,j,0)+x2*wstress_u(i,j,2) ! NOTE WSTRESS(:,:,0) & WSTRESS(:,:,2) sur points traceurs
       enddo; enddo
       loop_=vstrs_id 
       x2=(elapsedtime_now          -airseafile_prvtime(loop_))      &
         /(airseafile_nextime(loop_)-airseafile_prvtime(loop_))
! Bornage car l'attente de la synchro subcycling peut faire sortir (legerement) de l'intervale [O:1]
       x2=min(max(x2,0.d00),1.d00) !09-02-14
       x0=1.-x2
       do j=0,jmax+1 ; do i=0,imax+1
        xy_t(i,j,2)=x0*wstress_v(i,j,0)+x2*wstress_v(i,j,2) ! NOTE WSTRESS(:,:,0) & WSTRESS(:,:,2) sur points traceurs
       enddo; enddo
! Interpoler sur points U et V
       do j=0,jmax+1 ; do i=1,imax+1
        wstress_u(i,j,1)=0.5*(xy_t(i,j,1)+xy_t(i-1,j,1))
       enddo         ; enddo
       do j=1,jmax+1 ; do i=0,imax+1
        wstress_v(i,j,1)=0.5*(xy_t(i,j,2)+xy_t(i,j-1,2))
       enddo         ; enddo
        
      endif                        !pmx>

! Apres le vent (necessaire au calcul du flux de O2) il n'est pas utile d'aller
! plus loin dans le traitement des variables meteo dans le cas d'une simulation
! bio offline:
!     if(ioffline==2)return !04-10-14

      loop_=p0m_id
      x2=(elapsedtime_now          -airseafile_prvtime(loop_))      &
        /(airseafile_nextime(loop_)-airseafile_prvtime(loop_))
! Bornage car l'attente de la synchro subcycling peut faire sortir (legerement) de l'intervale [O:1]
      x2=min(max(x2,0.d00),1.d00) !09-02-14
      x0=1.-x2
      do j=0,jmax+1 ; do i=0,imax+1
       pss_w(i,j,1)=x0*pss_w(i,j,0)+x2*pss_w(i,j,2)
      enddo; enddo
      pss_mean(1)=x0*pss_mean(0)+x2*pss_mean(2)

! Pas necessaire d'aller au dela dans le cas ioffline=2
      if(ioffline==2)return !08-01-15

! Par defaut irelaxsst=0. Rappel de la SST si irelaxsst=1
!     if(irelaxsst==1)call airseaflux_relaxed_sst !06-06-15 deplace dans le driver

!     do loop_=1,1000
!     if(par%rank==0)write(6,*)'loop_=',loop_
!     teta2_t(:,:,1)=teta2_t(:,:,0) !BIDOUILLE
!        q2_t(:,:,1)=   q2_t(:,:,0) !BIDOUILLE

! Formules Bulk:
      if(flag_abl==1)call atmboundlayer_fieldcorrection ! Ajouter l'anomalie de Tair et Q produit par le modele de CLA


! CASE 1DV MODEL: les tableaux OBC constants sur (i,j)
      if(flag_1dv==1)call airseaflux_1dv 

      call bulk_formulae(flag_bulk_,0,1,relativewind)

! Anomalies de temperature et d'humidite de l'air resultant de l'incoherence
! des flux air/mer entre ceux calcules par le modele ocean et ceux du modele
! meteo
      
      if(flag_abl==1)call atmboundlayer_balance

!      write(6,*)'flag_abl=',flag_abl
!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      stop 'Chet'

      end subroutine airseaflux_ecmwfbk

!..............................................................................

      subroutine airseaflux_1dv !12-10-15
      implicit none

      do j=0,jmax+1 ; do i=0,imax+1
           ssr_w(i,j,1)=    ssr_w(2,2,1)
          snsf_w(i,j,1)=   snsf_w(2,2,1)
       precipi_w(i,j,1)=precipi_w(2,2,1)
         teta2_t(i,j,1)=  teta2_t(2,2,1)
            q2_t(i,j,1)=     q2_t(2,2,1)
         uwind_t(i,j,1)=  uwind_t(2,2,1)
         vwind_t(i,j,1)=  vwind_t(2,2,1)
           pss_w(i,j,1)=    pss_w(2,2,1)
      enddo; enddo

      end subroutine airseaflux_1dv

!..............................................................................

      subroutine airseaflux_1dv_glorys !10-03-16 Leo
      implicit none

             ssr_w(:,:,1)=    ssr_w(2,2,1)
            snsf_w(:,:,1)=   snsf_w(2,2,1)
         precipi_w(:,:,1)=precipi_w(2,2,1)
            slhf_w(:,:,1)=   slhf_w(2,2,1)
            sshf_w(:,:,1)=   sshf_w(2,2,1)
         wstress_u(:,:,1)=wstress_u(2,2,1)
         wstress_v(:,:,1)=wstress_v(2,2,1)

      end subroutine airseaflux_1dv_glorys
!..............................................................................

      subroutine airseaflux_physics(case_)
      implicit none
      integer case_
#ifdef synopsis
       subroutinetitle='airseaflux_physics'
       subroutinedescription= &
          'If necessary: deduces the specific humidity from the' &
       //' dewpoint, deduces the potential air temperature from' &
       //' the air temperaturer'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


      if(case_==2) then !222222222222222222>

      if(flag_humspec_dewpoint=='dewpoint')then  !dddddddddddddddddddd>
! on calcule l'humidite specifique !25-01-11 changement formules

       do j=0,jmax+1 ; do i=0,imax+1

        x1=6.1078*10.**((q2_t(i,j,2)-273.15)*7.5/(q2_t(i,j,2)-273.15+237.3))  ! pression vapeur saturante fct(point de rosee): unite mb
                                                                      !       ! verif tableau9 page 120 Queney
        x2=x1*100./pss_w(i,j,2)      ! e/p avec  passe de mb en Pascal
        q2_t(i,j,2)=0.622*x2/(1.-0.378*x2)                            ! expression 52 p115 Queney q=f(e/p)
!       x1 = 6.11 * exp(5417.7530 * ((1/273.16) - (1/q2_t(i,j,2))))   !  pression de vapeur saturante fction du pt de rosee
!       q2_t(i,j,2)=           &  ! humidite spec
!             0.622*x1          & ! 0.622 * Pression Vapeur Saturante
!             /pss_w(i,j,2)       !  Pression Surface
! Equivalent e:
!        q2_t(i,j,2)=const5*exp(-(5417.7530/q2_t(i,j,2)))/pss_w(i,j,2)

       enddo ; enddo
      endif                                      !dddddddddddddddddddd>

      if(flag_humspec_dewpoint=='humrelat')then  !huhuhuhuhuhuhuhuhuhu> !19-05-16

! humidite specifique = humidite relative fois humidite specifique saturante

#ifdef bidon
! Quelques constantes doivent etre initialisees:
       call initial_bulk_moon(0,0) ! Arguments sans importance

       do j=0,jmax+1 ; do i=0,imax+1
!       x1=0.98*exp(xalpw-xbetaw/teta2_t(i,j,2)  &
        x1=     exp(xalpw-xbetaw/teta2_t(i,j,2)  & ! sans 0.98 pour effet de sel
                      -xgamw*log(teta2_t(i,j,2)) )/pss_w(i,j,2)
! x2=Q saturante
        x2=(xrd/xrv)*x1/(1.+((xrd/xrv)-1.)*x1)
        q2_t(i,j,2)=x2*q2_t(i,j,2)

        if(i==imax/2.and.j==jmax/2.and.mask_t(i,j,kmax)==1)write(10+par%rank,*)q2_t(i,j,2),teta2_t(i,j,2),pss_w(i,j,2)

       enddo ; enddo
#endif
! Formule utilisant l'equation de Clausius-Clapeyron !11-04-18
       do j=0,jmax+1 ; do i=0,imax+1
          q2_t(i,j,2)=q2_t(i,j,2)*100.*exp(17.67*(teta2_t(i,j,2)-273.15)/(teta2_t(i,j,2)-29.65))/(0.263*pss_w(i,j,2))
       enddo         ; enddo

      endif                                      !huhuhuhuhuhuhuhuhuhu>

      return
      endif        !222222222222222222>

      if(case_==1) then !111111111111111111>

!       Transformer la temperature en temperature potentielle:
        const3=2.*grav*rhoair ! correction de pression hydrostatique e 2m
        do j=0,jmax+1 ; do i=0,imax+1
         teta2_t(i,j,2)=teta2_t(i,j,2)*(1.e5/(pss_w(i,j,2)-const3))**.286
        enddo ; enddo

! commente le !12-04-16
!      if(flag_abl==1) then !ablablabl>
!       do j=0,jmax+1 ; do i=0,imax+1
!        teta0_t(i,j,2)=teta0_t(i,j,2)*(1.e5/(pss_w(i,j,2)       ))**.286
!       enddo ; enddo
!      endif                !ablablabl>

      endif        !111111111111111111>

      end subroutine airseaflux_physics

!..............................................................................

      subroutine airseaflux_rotate
      implicit none
#ifdef synopsis
       subroutinetitle='airseaflux_rotate'
       subroutinedescription=                                       &
          'Computes the along-axis components of the wind from the' &
       //' the E-W S-N wind components'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      do j=0,jmax+1 ; do i=0,imax+1
       xy_t(i,j,1)=gridrotcos_t(i,j)*uwind_t(i,j,2)       &
                  -gridrotsin_t(i,j)*vwind_t(i,j,2)
       xy_t(i,j,2)=gridrotsin_t(i,j)*uwind_t(i,j,2)       &
                  +gridrotcos_t(i,j)*vwind_t(i,j,2)
      enddo ; enddo
      do j=0,jmax+1 ; do i=0,imax+1
       uwind_t(i,j,2)=xy_t(i,j,1)
       vwind_t(i,j,2)=xy_t(i,j,2)
      enddo ; enddo

      if(flag_abl==1) then !>>>>>>>

       do j=0,jmax+1 ; do i=0,imax+1
        xy_t(i,j,1)=gridrotcos_t(i,j)*uwindabl_t(i,j,1,2)       &
                   -gridrotsin_t(i,j)*vwindabl_t(i,j,1,2)
        xy_t(i,j,2)=gridrotsin_t(i,j)*uwindabl_t(i,j,1,2)       &
                   +gridrotcos_t(i,j)*vwindabl_t(i,j,1,2)
       enddo ; enddo
       do j=0,jmax+1 ; do i=0,imax+1
        uwindabl_t(i,j,1,2)=xy_t(i,j,1)
        vwindabl_t(i,j,1,2)=xy_t(i,j,2)
       enddo ; enddo

      endif                !>>>>>>>

      if(flag_wstressbulk==0) then !pmx> !25-01-18
! NOTE: WSTRESS(:,:,2) sur points traceurs 
       do j=0,jmax+1 ; do i=0,imax+1
        xy_t(i,j,1)=gridrotcos_t(i,j)*wstress_u(i,j,2)       &
                   -gridrotsin_t(i,j)*wstress_v(i,j,2)
        xy_t(i,j,2)=gridrotsin_t(i,j)*wstress_u(i,j,2)       &
                   +gridrotcos_t(i,j)*wstress_v(i,j,2)
       enddo ; enddo
       do j=0,jmax+1 ; do i=0,imax+1
        wstress_u(i,j,2)=xy_t(i,j,1)
        wstress_v(i,j,2)=xy_t(i,j,2)
       enddo ; enddo
      endif                        !pmx>

      end subroutine airseaflux_rotate
!..............................................................................

      subroutine airseaflux_varname(loop_)
      implicit none
      integer loop_
#ifdef synopsis
       subroutinetitle='airseaflux_varname'
       subroutinedescription= &
          'Defines the variable names used to inquire the variables' &
       //' in the netcdf meteo files.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(flag_meteodata=='glorys') then !ggggggggg>
          if(loop_==ustrs_id)then;meteovarname(1)='sozotaux';return;endif
          if(loop_==vstrs_id)then;meteovarname(1)='sometauy';return;endif
          if(loop_==ssr_id  )then;meteovarname(1)='soceshwf';meteovarname(2)='rsds';return;endif
          if(loop_==slhf_id )then;meteovarname(1)='socelatf';return;endif
          if(loop_==netir_id)then;meteovarname(1)='socelowf';return;endif
          if(loop_==sshf_id )then;meteovarname(1)='socesenf';return;endif
          if(loop_==rain_id )then;meteovarname(1)='sowaprec';return;endif
      stop ' STOP airseaflux_varname varname non trouve cas glorys'
      endif                             !ggggggggg>

      if(flag_meteodata/='glorys') then !ooooooooo>

      meteovarname='undefined' !11-04-18

        if(loop_==ssr_id) then
! SSR PRIORITAIRE:
         meteovarname(1)='SSR'   
         meteovarname(2)='ssr'   ! 24-04-16
         meteovarname(3)='SSRD'  ! Necessite calcul albedo
         meteovarname(4)='ssrd'  ! Necessite calcul albedo
! SSRD PRIORITAIRE:
!        meteovarname(1)='SSRD'  ! Necessite calcul albedo
!        meteovarname(2)='ssrd'  ! Necessite calcul albedo
!        meteovarname(3)='SSR'   
!        meteovarname(4)='ssr'   ! 24-04-16
         meteovarname(5)='rsds'
         meteovarname(6)='flsolaire' !08-11-16
         meteovarname(7)='SWDOWN'    !17-09-21
         return
        endif

        if(loop_==ustrs_id) then
          meteovarname(1)='EWSS'
         return
        endif
        if(loop_==vstrs_id) then
          meteovarname(1)='NSSS'
         return
        endif

        if(loop_==u10m_id) then
          meteovarname(1)='U10M' ; meteovarname(2)='u10m'
          meteovarname(3)='U10'  ; meteovarname(4)='u10'
          meteovarname(5)='uas' !14-09-22
         return
        endif
        if(loop_==v10m_id) then
          meteovarname(1)='V10M' ; meteovarname(2)='v10m'
          meteovarname(3)='V10'  ; meteovarname(4)='v10'
          meteovarname(5)='vas' !14-09-22
         return
        endif
        if(loop_==t2m_id) then
         meteovarname(1)='T2M' ; meteovarname(2)='t2m'
         meteovarname(3)='T2'  !17-09-21
         meteovarname(4)='tas' !14-09-22
         return
        endif
        if(loop_==dp2m_id) then
!        meteovarname(1)='hu2m'                        !19-05-16
         meteovarname(1)='hrel2m'                        !19-05-16
         meteovarname(2)='D2M' ; meteovarname(3)='d2m'
         meteovarname(4)='Q2M' ; meteovarname(5)='q2m' !28-04-16
         meteovarname(6)='Q2'   !17-09-21
         meteovarname(7)='huss' !13-09-22
         return
        endif
        if(loop_==p0m_id) then
         meteovarname(1)='MSL' ; meteovarname(2)='msl' !21-06-16
! CommentE le 20-11-17 pour inciter A demander MSL. Decommenter si
! necessaire
!        meteovarname(3)='SP'  ; meteovarname(4)='sp'
         meteovarname(3)='???'  ; meteovarname(4)='???'
         meteovarname(5)='pmer'
         meteovarname(6)='PSFC' !17-09-21
         meteovarname(7)='psl'  !14-09-22
         return
        endif
        if(loop_==ir_id) then
         meteovarname(1)='STRD' ; meteovarname(2)='strd'
         meteovarname(3)='fltherm' !08-11-16
         meteovarname(4)='GLW'     !17-09-21
         meteovarname(5)='rlds'    !12-09-22
         return
        endif
        if(loop_==rain_id) then
         meteovarname(1)='TP' ; meteovarname(2)='tp'
         meteovarname(3)='eau' !08-11-16
         meteovarname(4)='RAINNC' !17-09-21
         meteovarname(5)='pr'     !12-09-22
         return
        endif

        if(loop_==t0m_id) then
         meteovarname(1)='t0m' ; meteovarname(2)='T0M'
         meteovarname(3)='sst' ; meteovarname(4)='SST'
         return
        endif
        if(loop_==u100m_id) then
         meteovarname(1)='ublh' !11-05-16
         meteovarname(2)='u100'  ; meteovarname(3)='U100'
         meteovarname(4)='u100m' ; meteovarname(5)='U100M'
         return
        endif
        if(loop_==v100m_id) then
         meteovarname(1)='vblh' !11-05-16
         meteovarname(2)='v100'  ; meteovarname(3)='V100'
         meteovarname(4)='v100m' ; meteovarname(5)='V100M'
         return
        endif
        if(loop_==abl_id) then
         meteovarname(1)='blh' ; meteovarname(2)='BLH'
         return
        endif

      if(par%rank==0) then !>>>>>
      write(6,*)'ssr_id,u10m_id,v10m_id,t2m_id,dp2m_id,p0m_id' &
               ,'ir_id,rain_id,t0m_id,u100m_id,v100m_id,abl_id'&
               ,ssr_id,u10m_id,v10m_id,t2m_id,dp2m_id,p0m_id   &
               ,ir_id,rain_id,t0m_id,u100m_id,v100m_id,abl_id
      write(6,'(a,i0,3a,5a)')                           &
               'loop_=',loop_                           &
              ,' flag_meteodata=',trim(flag_meteodata)  &
              ,' meteovarname=',meteovarname(1:5)

      stop ' STOP airseaflux_varname liste binrec inexistante'
      endif                !>>>>>

      endif                             !ooooooooo>

      end subroutine airseaflux_varname

!..............................................................................

      subroutine airseaflux_inquire_var(loop_,ncid_)
      implicit none
      integer loop_,ncid_,varid_
#ifdef synopsis
       subroutinetitle='airseaflux_inquire_var'
       subroutinedescription= &
       'Inquires and reads the variable in the meteo netcdf files'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      flag_cumul=0 ! pas de cumul
      scalefct=1.  ! pas de chgt d'unites
      texte90=''   ! reset

      if( loop_==ir_id.or.   &
          loop_==ssr_id.or.  &
          loop_==netir_id) then                          !*****>

!      write(6,*)'loop_ ir_id ssr_id netir_id',loop_,ir_id,ssr_id,netir_id
!      write(6,'(a,a)')'meteovarname=',meteovarname(1)

                   status=nf_inq_varid(ncid_,meteovarname(1),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(2),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(3),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(4),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(5),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(6),varid_) !08-11-16
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(7),varid_) !17-09-21
      if(status/=0) then !----->
       write(6,*)'ir_id ssr_id netir_id loop_=' &
                 ,ir_id,ssr_id,netir_id,loop_
       write(6,'(8(a,1x))')trim(meteovarname(1)) &
                          ,trim(meteovarname(2)) &
                          ,trim(meteovarname(3)) &
                          ,trim(meteovarname(4)) &
                          ,trim(meteovarname(5)) &
                          ,trim(meteovarname(6)) &
                          ,trim(meteovarname(7)) &
                          ,'not found in the netcdf file' &
                          ,trim(texte80(1))
       stop 'Erreur nf_inq_varid airseaflux_inquire_var 1'
      endif              !----->

      status=nf_get_att_text(ncid_,varid_,'units',texte90);if(status/=0)stop 'echec attribut units pour ssr'
!      write(6,*)'Unites=',trim(texte90)

      flag_stop=1
      if(   texte90=='W m**-2 s'            &
        .or.texte90=='J m**-2') then !pmx> !14-09-22
                                flag_stop=0
                                flag_cumul=1 ! cumul
      endif                          !pmx>
      if(texte90=='W m-2') then !ooo> !14-09-22
                                flag_stop=0
                                flag_cumul=0 ! pas de cumul
      endif                     !ooo>

! Il semblerait (ca sera A reverifier) qu'Arome se trompe dans les unites et qu'il manque "s". 
! Du coup procedure d'exception afin de ne pas interferer avec le cas WRF:
       if(flag_meteodata=='arome'.and.texte90=='W m-2') then !>>>
        flag_stop=0
        flag_cumul=1 !17-09-21 
       endif                                                 !>>>

       call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
       if(k0/=0) then !pmx>
       write(10+par%rank,*)'Units not found in airseaflux_inquire_var'
       write(10+par%rank,*)'ir_id,ssr_id,netir_id,loop_',ir_id,ssr_id,netir_id,loop_
       write(10+par%rank,*)'texte90=',trim(texte90)
       stop 'Units not found airseaflux_inquire_var. See fort.xxx files'
       endif          !pmx>

      if(par%rank==0) then !RK0>
      if(loop_==ir_id)write(6,*)'ir_id flag_cumul=',flag_cumul
      if(loop_==ssr_id)write(6,*)'ssr_id flag_cumul=',flag_cumul
      if(loop_==netir_id)write(6,*)'netir_id flag_cumul=',flag_cumul
      endif                !RK0>

!      write(10+par%rank,*)loop_,flag_cumul
       
      endif                                                           !*****>

      if(loop_==rain_id) then                                         !*****>

!      write(6,*)'loop_ rain_id',loop_,rain_id
!      write(6,'(a,a)')'meteovarname=',meteovarname(1)

                   status=nf_inq_varid(ncid_,meteovarname(1),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(2),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(3),varid_) !08-11-16
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(4),varid_) !17-09-21
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(5),varid_) !12-09-22
      if(status/=0) then !>>>
        write(6,*)'rain_id name not found in name list ',meteovarname(:) &
         ,' in netcdf file ',trim(texte80(1))
        stop 'Erreur nf_inq_varid airseaflux_inquire_var rain_id'
      endif              !>>>
      status=nf_get_att_text(ncid_,varid_,'units',texte90);if(status/=0)stop 'echec attribut units pour ssr'

      flag_stop=1 !17-09-21
      if(texte90=='kg/m2') then       !------->
       flag_cumul=1 ;  scalefct=0.001 ! conversion en m
       flag_stop=0 !17-09-21
      endif                           !------->
      if(texte90=='kg m-2') then      !-------> !08-11-16
       flag_cumul=1 ;  scalefct=0.001 ! conversion en m
       flag_stop=0 !17-09-21
      endif                           !------->
      if(texte90=='m') then           !------->
       flag_cumul=1
       flag_stop=0 !17-09-21
      endif                           !------->
      if(texte90=='mm') then          !------->
       flag_cumul=1 ;  scalefct=0.001 ! conversion en m
       flag_stop=0 !17-09-21
      endif                           !------->
      if(texte90=='kg m-2 s-1') then  !------->
       flag_cumul=0 ;  scalefct=0.001 ! conversion en m/s
       flag_stop=0 !17-09-21
      endif                           !------->
      if(par%rank==0) then !pmx>
        write(6,*)'Unites precipitations: ',trim(texte90) ! 14-09-22
        write(6,*)'precipitations flag_cumul,scalefct=',flag_cumul,scalefct
      endif                !pmx>
      if(flag_stop==1)stop 'Err units inconnues pour precip' !17-09-21

!      write(6,*)'Unites=',trim(texte90)
!      write(6,*)'flag_cumul=',flag_cumul
!      write(6,*)'scalefct=',scalefct
!      if(loop_==netir_id)stop 'sof?'

!      write(10+par%rank,*)loop_,flag_cumul,scalefct

      endif                                                          !*****>

!     if(loop_==dp2m_id) then                                         !------->
!                  status=nf_inq_varid(ncid_,meteovarname(1),varid_)
!     if(status/=0)status=nf_inq_varid(ncid_,meteovarname(2),varid_)
!     if(status/=0)stop 'Erreur nf_inq_varid airseaflux_inquire_var 3'
!     status=nf_get_att_text(ncid_,varid_,'long_name',texte90) ; if(status/=0)stop 'echec attribut long_name dp2m'
!     flag_humspec_dewpoint='humspeci'
!     if(index(texte90,'dewpoint')/=0)flag_humspec_dewpoint='dewpoint'
!     endif                                                           !------->

       if(loop_==ustrs_id)flag_cumul=1 !25-01-18
       if(loop_==vstrs_id)flag_cumul=1

      end subroutine airseaflux_inquire_var

!..............................................................................

      subroutine airseaflux_interp_driver(var_id_,t_,var_id_file_)
      implicit none
      integer var_id_,t_,vstart1_,vstart2_,ncid_,var_id_file_
      double precision :: time1_=0.,time2_=0. !28-10-18
#ifdef synopsis
       subroutinetitle='airseaflux_interp_driver'
       subroutinedescription='Driver of interpolation subroutines'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Obtenir le nom du fichier netcdf, texte80, et les numeros d'echeances time_ et time2_
      call airseaflux_get_time_from_binrecfile(var_id_file_,t_,vstart1_,vstart2_,time1_,time2_)

      call airseaflux_varname(var_id_) ! donne meteovarname le nom de la variable du fichier netcdf

      call airseaflux_read_only(vstart1_,1)

! Cas de variables cumules:
      if(flag_cumul==1) then !1111111111111111>

      k0=0
      if(meteo_t0=='file'.and.vstart1_>vstart2_)k0=1
      if(meteo_t0=='time'.and.time1_>time2_)    k0=1
      if(meteo_t0=='user') then !pmx>
        if(modulo(nint(time2_/3600.),meteo_cumul_modulo)/=0)k0=1 ! note: oui c'est le temps precedent que l'on teste
      endif                     !pmx>
      if(meteo_t0=='auto') then !ooo>
         if(var_id_==ssr_id) then !m°v°m>
! note: si autocumul_status=1 le cumul se poursuit
!       si autocumul_status=0 le cumul repart de zero
            do j=1,meteo_jmax ; do i=1,meteo_imax
             meteo_cum(i,j)=meteo_var(i,j)
            enddo ; enddo
            call airseaflux_read_only(vstart2_,2)
            if(time1_-time2_<0)time2_=0. ; x1=scalefct/(time1_-time2_)
            autocumul_status=1 ; i=1 ; j=1
! pour eviter erreurs de precision on teste en suivant par rapport a -1
! plutot que par rapport a zero:
            if((meteo_cum(i,j)-meteo_var(i,j))*x1<-1.)autocumul_status=0
            call mpi_bcast(autocumul_status,1,mpi_integer,0,par%comm2d,ierr)
! Recharger meteo_var avec l'echeance 2 pour etre coherent avec l'etape A suivre:
            call airseaflux_read_only(vstart1_,1)
         endif                    !m°v°m>
         k0=autocumul_status
      endif                     !ooo>
      if(meteo_t0=='unli')k0=1  ! 'unlimited' le cumul n'est jamais remis A zero !17-09-21


! Debugging help (commente le !16-11-18 car defectueux)
!     if(meteo_t0/='user'.and.time1_>86400.) then !m°v°m> !28-10-18
!     write(10+par%rank,*)'Oups! time1_=',time1_,' >24h'
!     write(10+par%rank,*)'The ECMWF file contains more than 24h and' &
!     ,' strangely the decumulation mode is not "user". '             &
!     ,' Delete stop lines if I worry unnecessarily'
!     stop 'Warning notebook_airseaflux reported in fort.xx files'
!     endif                                   !m°v°m> !28-10-18


        if(k0==1) then !---------------->

            do j=1,meteo_jmax ; do i=1,meteo_imax
             meteo_cum(i,j)=meteo_var(i,j)
            enddo ; enddo
            call airseaflux_read_only(vstart2_,2)

            if(time1_-time2_==0.)stop ' Stop 1281 time1_-time2_==0.'

            x1=scalefct/(time1_-time2_)
            do j=1,meteo_jmax ; do i=1,meteo_imax
             meteo_var(i,j)=(meteo_cum(i,j)-meteo_var(i,j))*x1
            enddo ; enddo

        else           !---------------->

! Si le compteur temps a EtE remis A zero alors time1_-time2_<0 
! et on impose time2_=0. sinon cela signifie que le compteur temps
! continue d'etre incremente (cas des nouveaux fichiers CLAS) et donc
! il faut continuer de normaliser par time1- - time2_ !14-02-16
            if(time1_-time2_<0)time2_=0. !14-02-16
            if(time1_-time2_==0.) &
            stop ' Stop airseaflux_interp_driver time1_-time2_==0.' !14-02-16
            x1=scalefct/(time1_-time2_) !14-02-16

            do j=1,meteo_jmax ; do i=1,meteo_imax
             meteo_var(i,j)=meteo_var(i,j)*x1
            enddo ; enddo

        endif          !---------------->

      else                   !1111111111111111> !10-12-22
            do j=1,meteo_jmax ; do i=1,meteo_imax
             meteo_var(i,j)=meteo_var(i,j)*scalefct
            enddo ; enddo
      endif                  !1111111111111111>

      if(var_id_==p0m_id.and.flag_p0m_filter==1)call airseaflux_2dnoiseremover(var_id_) !01-07-15
      call airseaflux_interpolation(var_id_) ! calcule xy_t(i,j,1), le champs interpole

! La ligne suivante est deplacee e la fin de airseaflux_get_time_from_binrecfile....
!     airseafile_nextrec(var_id_)=airseafile_nextrec(var_id_)+1

      end subroutine airseaflux_interp_driver

!..............................................................................

      subroutine airseaflux_interpolation(var_id_)
      implicit none
      integer var_id_
#ifdef synopsis
       subroutinetitle='airseaflux_interpolation'
       subroutinedescription='Applies the bilinear interpolation scheme'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Bouchage des trous pour certaines variables (vent exclus):
      if(flag_meteo_land_plug==0) then !m°v°m> !13-07-18
       k0=0 ! Valeurs en terre prises en comptes
      else                             !m°v°m>
       k0=1 ! Valeurs en terre remplacees par valeurs en mer la plus proche
      endif                            !m°v°m>

! Bouchage des trous pour le vent seulement: !27-12-21
      if(flag_meteo_land_plug_wind==0) then !m°O°m> !27-12-21
       if(var_id_==u10m_id) k0=0 !Valeurs en terre prises en comptes
       if(var_id_==v10m_id) k0=0 !Valeurs en terre prises en comptes
      else                                  !m°O°m> !27-12-21
       if(var_id_==u10m_id) k0=1 !Valeurs en terre remplacees par valeurs en mer la plus proche
       if(var_id_==v10m_id) k0=1 !Valeurs en terre remplacees par valeurs en mer la plus proche
      endif                                 !m°O°m> !27-12-21

      if(flag_abl==1) then !pmxpmx> ! 12-03-17
! Hors cas evident de la SABL, il vaut mieux ne pas utiliser le vent
! au dessus du continent:
       if(var_id_==u10m_id) k0=0
       if(var_id_==v10m_id) k0=0
       if(var_id_==u100m_id)k0=0
       if(var_id_==v100m_id)k0=0
       if(var_id_==abl_id)  k0=0
      endif                !pmxpmx> ! 12-03-17

      if(k0==1)call appli_bouchetrou_meteo

! DEBUG: a la premiere iteration verifier les depassements de tableaux !10-07-18
      if(iteration3d==0) then !debug>
       flag_stop=0
       do j=0,jmax+1 ; do i=0,imax+1
        i1=int( ij2meteo_i(i,j) - (meteozoom_istr-1) ) ; j1=int( ij2meteo_j(i,j) - (meteozoom_jstr-1) )
        if(i1<1.or.i1>meteo_imax.or.j1<1.or.j1>meteo_jmax) then !ooo>
         flag_stop=1 
         write(10+par%rank,*)'par%rank ',par%rank,' i,j loc ',i,j,' i,j glob ',i+par%timax(1),j+par%tjmax(1) &
         ,' i1 j1 ',i1,j1,' meteozoom_istr,meteozoom_jstr ',meteozoom_istr,meteozoom_jstr
        endif                                                   !ooo>
       enddo         ; enddo
       call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
       if(k0/=0)stop 'Err 1215 grid overflow in module_airseaflux. See fort.xxx error files'
      endif                   !debug>

! Interpoler:
!      const1=1./meteo_londlt
!      const2=1./meteo_latdlt
      do j=0,jmax+1
      do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
!       deci=1.+( lon_t(i,j)*rad2deg-meteo_lonstr)*const1
!       i1=int(deci)
!       rapi=deci-i1
!       i1=i1-meteozoom_istr+1                              ! conservation la parallelisation

! Indice decimale j (latitude) dans grille aladin:
!       decj=1.+( lat_t(i,j)*rad2deg-meteo_latstr)*const2
!       j1=int(decj)
!       rapj=decj-j1
!       j1=j1-meteozoom_jstr+1                              ! conservation la parallelisation

!       deci=ij2meteo_i(i,j) ; i1=int(deci) ; rapi=deci-i1
!       decj=ij2meteo_j(i,j) ; j1=int(decj) ; rapj=decj-j1

        i1=int( ij2meteo_i(i,j) - (meteozoom_istr-1) ) !21-03-25
        j1=int( ij2meteo_j(i,j) - (meteozoom_jstr-1) )
        rapi=ij2meteo_i(i,j)-int(ij2meteo_i(i,j))
        rapj=ij2meteo_j(i,j)-int(ij2meteo_j(i,j))

       xy_t(i,j,1)= (1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )         &
                   +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)         &
                   +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )         &
                   +    rapi *    rapj *meteo_var(i1+1,j1+1)

      enddo
      enddo

      end subroutine airseaflux_interpolation

!..............................................................................

      subroutine airseaflux_read_only(vstart_,ttx_)
      implicit none
      integer vstart_,ttx_,ncid_
#ifdef synopsis
       subroutinetitle='airseaflux_read_only'
       subroutinedescription= &
       'Reads the meteo variable in the netcdf file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0)write(6,'(a,a5,a,a,i4)')'airseaflux_read_only ' &
                                            ,meteovarname(1),' '     &
                                            ,trim(texte80(ttx_)),vstart_

      status=nf_open(trim(texte80(ttx_)),nf_nowrite,ncid_);if(status/=0)stop 'erreur nf_open airseaflux_read_only'
                   status=nf_inq_varid(ncid_,meteovarname(1),var_id)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(2),var_id)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(3),var_id)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(4),var_id)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(5),var_id)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(6),var_id) !08-11-16
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(7),var_id) !17-09-21
      if(status/=0) then !------>
       write(6,'(8(a,1x))')trim(meteovarname(1)) &
                          ,trim(meteovarname(2)) &
                          ,trim(meteovarname(3)) &
                          ,trim(meteovarname(4)) &
                          ,trim(meteovarname(5)) &
                          ,trim(meteovarname(6)) &
                          ,trim(meteovarname(7)) &
                          ,'not found in the netcdf file' &
                          ,trim(texte80(ttx_))
       stop 'Erreur nf_inq_varid meteovarname 807'
      endif

      status=nf_inq_var(ncid_,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(status.ne.0)stop 'echec nf_inq_var 822'

      if(var_dims<=3) then !pmxpxm> !28-11-16
       varstart(3)=vstart_        ; varcount(3)=1          ! time
       varstart(2)=meteozoom_jstr ; varcount(2)=meteo_jmax ! lat
       varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax ! lon
      else                 !pmxpmx> !28-11-16
       varstart(4)=vstart_        ; varcount(4)=1          ! time
       varstart(3)=1              ; varcount(3)=1          ! height
       varstart(2)=meteozoom_jstr ; varcount(2)=meteo_jmax ! lat
       varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax ! lon
      endif                !pmxpmx>


!.................
! Debug: !11-05-16
      do k=1,var_dims
       if(varstart(k)<1) then !pmxpmx>
       write(6,'(a,i0,a,a)')                         & !11-05-16
        ' Err 1046 varstart(k)<1 dim ',k             &
       ,' meteovarname(1)=',trim(meteovarname(1))
       stop ' STOP subroutine airseaflux_read_only'
       endif                  !pmxpmx>
      enddo
!.................


      if(var_type==nf_real) then  !rrrrr>
       status=nf_get_vara_real(ncid_,var_id        &
                                    ,varstart(1:var_dims) &
                                    ,varcount(1:var_dims) &
             ,meteo_var(1:meteo_imax,1:meteo_jmax))
       if(status/=0)stop 'Erreur nf_get_vara_real airseaflux_read_only'
      endif                       !rrrrr>
      if(var_type==nf_short) then !iiiii>
       status=nf_get_vara_int(ncid_,var_id        &
                                   ,varstart(1:var_dims) &
                                   ,varcount(1:var_dims) &
          ,meteo_short(1:meteo_imax,1:meteo_jmax))
       if(status/=0) then !----->
        write(6,'(a,a)')'netcdf file=',trim(texte80(ttx_))
        write(6,'(5a)')'meteovarname ',meteovarname(:)
        write(6,*)'meteo_imax,meteo_jmax ',meteo_imax,meteo_jmax
        write(6,*)'ubound(meteo_short) ',ubound(meteo_short)
        write(6,*)'vardims=',var_dims
        write(6,*)'vstart_ ',vstart_
        write(6,*)'varstart ',varstart(1:var_dims)
        write(6,*)'varcount ',varcount(1:var_dims)
        stop 'Erreur nf_get_vara_int airseaflux_read_only'
       endif              !----->
       status=nf_get_att_real(ncid_,var_id,'scale_factor',var_scalefactor) !24-03-15
       if(status/=0) &
       stop 'error get scale_factor 843'
       status=nf_get_att_real(ncid_,var_id,'add_offset',var_addoffset) !24-03-15
       if(status/=0)stop 'error get add_offset 843'
       do j=1,meteo_jmax ; do i=1,meteo_imax
        meteo_var(i,j)=meteo_short(i,j)*var_scalefactor+var_addoffset
       enddo ; enddo
      endif                      !iiiii>


      status=nf_close(ncid_)                           ;if(status/=0)stop 'erreur nf_close airseaflux_read_only'


      end subroutine airseaflux_read_only

!..............................................................................

      subroutine airseaflux_get_time_from_binrecfile(loop_,t_,vstart1_,vstart2_,time1_,time2_)
      implicit none
      integer loop_,t_,vstart1_,vstart2_
      double precision time1_,time2_
#ifdef synopsis
       subroutinetitle='airseaflux_get_time_from_binrecfile'
       subroutinedescription= &
          'Gets next reading time and corresponding meteo file in the' &
       //' binary file. Increment the corresponding record number.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

 1100 continue
      nc=airseafile_nextrec(loop_)

!---------------------------------------------------------------------!05-05-18
      if(par%rank==0) then !0000000000000000000>
      open(unit=4,file=airseabinreclist(loop_)                     &
                 ,access='direct',recl=540,form='unformatted')
      read(4,rec=max(nc-1,1))texte80(2),vstart2_                   &
                                       ,airseafile_prvtime(loop_)  & ! Echeance avant
                                       ,time2_
      read(4,rec=nc         )texte80(1),vstart1_                   &
                                       ,airseafile_nextime(loop_)  & ! Echeance apres
                                       ,time1_                     &
                                       ,flag_cumul                 &
                                       ,scalefct                   &
                                       ,flag_refresh_interp !26-01-22
      close(4)
      endif                !0000000000000000000>
! par%rank=0 envoie les valeurs aux autres par%rank
      call mpi_bcast(vstart1_  ,1,mpi_integer,0,par%comm2d,ierr)
      call mpi_bcast(vstart2_  ,1,mpi_integer,0,par%comm2d,ierr)
      call mpi_bcast(flag_cumul,1,mpi_integer,0,par%comm2d,ierr)
      call mpi_bcast(time1_    ,1,mpi_double_precision,0,par%comm2d,ierr)
      call mpi_bcast(time2_    ,1,mpi_double_precision,0,par%comm2d,ierr)
      call mpi_bcast(scalefct  ,1,mpi_double_precision,0,par%comm2d,ierr)
      call mpi_bcast(flag_refresh_interp,1,mpi_integer,0,par%comm2d,ierr) !26-01-22

      x0=airseafile_prvtime(loop_)
      call mpi_bcast(x0,1,mpi_double_precision,0,par%comm2d,ierr)
      airseafile_prvtime(loop_)=x0

      x0=airseafile_nextime(loop_)
      call mpi_bcast(x0,1,mpi_double_precision,0,par%comm2d,ierr)
      airseafile_nextime(loop_)=x0

      texte250=texte80(1)
      k=len(texte250)
      call mpi_bcast(texte250,k,mpi_char,0,par%comm2d,ierr)
      texte80(1)=texte250

      texte250=texte80(2)
      k=len(texte250)
      call mpi_bcast(texte250,k,mpi_char,0,par%comm2d,ierr)
      texte80(2)=texte250
!---------------------------------------------------------------------!05-05-18

! Verifier la continuite des champs de la liste et remettre a jour la
! chaine d'interpolation si discontinuite detectee:
!     call airseaflux_checkfileconsistecy !05-05-16
! Ligne ci-dessus remplacee par call airseaflux_refresh_interp_proc
! Note: pourquoi conditionnE A loop_==1: parce que la premiere variable,SSR, est en principe cumulee. 
! En plus d'etre premiere, elle arrive donc en avance d'1/2 echeance sur les champs instantannES qui 
! beneficieront de la mise A jour de la regle d'interpolation
      if(loop_==1.and. & 
         flag_refresh_interp==1)call airseaflux_refresh_interp_proc !26-01-22

! Move forward du prochain numero de record:
      airseafile_nextrec(loop_)=airseafile_nextrec(loop_)+1


! En dehors du cas particulier de l'etat initial, si l'echeance ecmwf persiste A Etre avant 
! elapsedtime_now alors essayer l'echeance suivante: !31-03-16
      if(initial_main_status==1) then !pmxpmx> 
         if(elapsedtime_now>airseafile_nextime(loop_)) goto 1100 !31-03-16
      endif                           !pmxpmx>

! Verifications de l'etat initial:
      if(iteration3d==0) then !0000000000000000000>
        if(t_==0.and.airseafile_nextime(loop_)>elapsedtime_now) then  !------->
         write(6,*)'loop_,airseafile_nextime(loop_),elapsedtime_now'  &
                   ,loop_,airseafile_nextime(loop_),elapsedtime_now  
         stop 'erreur1 airseaflux_get_time_from_binrecfile'
        endif                                                         !------->
        if(t_==2.and.airseafile_nextime(loop_)<elapsedtime_now) then  !------->
         stop 'erreur2 airseaflux_get_time_from_binrecfile'
        endif                                                         !------->
        if(t_==2.and.airseafile_nextime(loop_)==elapsedtime_now) then !------->
         stop 'erreur3 airseaflux_get_time_from_binrecfile'
        endif                                                         !------->
      endif                   !0000000000000000000>

      end subroutine airseaflux_get_time_from_binrecfile

!..............................................................................

      subroutine airseaflux_fbk(ichoix)
      implicit none
      integer  ichoix
! Note sur les precipitations. Le model_ meteo donne generalement      !19/03/03
! les precipitations sous forme de hauteurs cumulees,
! (par exemple x milimetres de pluie). Dans le model_, PRECIPI_Z est homogene
! e une tendance, donc homogene e une hauteur par seconde.
#ifdef synopsis
       subroutinetitle='airseaflux_fbk'
       subroutinedescription= &
       'Bulk formulae driver for other than ECMWF cases'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCe
!     PARTIE I
!     INITIALISATION DES FORCAGES EN DEBUT DE SIMULATION:
!     DEBUT:
      if (ichoix==1) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(iairsea.ne.2.and.iairsea.ne.3)stop 'airseaflux_fbk.f erreur1'   !060104
      if(ifb.ne.1)    stop 'airseaflux_fbk.f erreur2'

      do k=1,nairsea                                                   !22/01/08
       i1=dateairsea(1,k)
       i2=dateairsea(2,k)
       i3=dateairsea(3,k)
       i4=dateairsea(4,k)
       i5=dateairsea(5,k)
       i6=dateairsea(6,k)
       call datetokount(i1,i2,i3,i4,i5,i6)    !03-11-10
       airseadt(k,2)=elapsedtime_out          !15-04-11
       airseadt(k,1)=airseainfo(k,1)*3600.    !15-04-11
      enddo

      i0=0                                                             !11/07/02
      i2=2                                                             !11/07/02

      k4=nairsea
! cas offline on ne lit que le flux solaire pour variables bio
      if(ioffline==2)k4=1                                            !16/06/06

! Soit:
! Lire et interpoller les fichiers meteo NETCDF de la grille native
      if(airseaoption==1) then !>>netcdf>>>                            !04-02-10
        call airseaflux_initial
        call read_interp_netcdf_airseafile(0) ! read previous fields
        call read_interp_netcdf_airseafile(2) ! read next fields
      endif                    !>>netcdf>>>                            !04-02-10

! Soit:
! Lire et interpoller les fichiers meteo BINAIRES de la grille native
      if(airseaoption==2)call read_interp_initial_file
! Soit:
! Lire les fichiers meteo deje interpolles sur la grille de symphonie:
      if(airseaoption==0)call read_initial_file

! interpolation temporelle entre 2 echeances:
      call airseaflux_bulkvar_now

! Calculs des flux par les formules bulk (cas special de l'etat initial)
      call bulk_formulae(1,0,1,relativewind)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     PARTIE I
!     INITIALISATION DES FORCAGES EN DEBUT DE SIMULATION:
!     FIN.
      endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!                           /    /    /

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!            PARTIE II:
!     EVOLUTION DES FORCAGES ATMOSPHERIQUES AU COURS DE LA SIMULATION
! DEBUT:
      if (ichoix==2) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      i0=0                                                             !11/07/02
      i2=2                                                             !11/07/02

      k4=nairsea
! cas offline on ne lit que le flux solaire pour variables bio
      if(ioffline==2)k4=1                                            !16/06/06

! Mise e jour des champs:
! Soit:
! Lire et interpoller les fichiers meteo NETCDF de la grille native:
      if(airseaoption==1) call read_interp_netcdf_airseafile(2) ! read next fields  !04-02-10

! Soit:
! Lire et interpoller les fichiers meteo BINAIRES de la grille native:
      if(airseaoption==2)call update_read_interp_file
! Soit:
! Les fichiers sont deje interpolles sur la grille symphonie:
      if(airseaoption==0)call update_read_file

! interpolation temporelle entre 2 echeances:
      call airseaflux_bulkvar_now

! Calculs des flux par les formules bulk
      call bulk_formulae(2,0,1,relativewind)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!            PARTIE II:
!     EVOLUTION DES FORCAGES ATMOSPHERIQUES AU COURS DE LA SIMULATION
!     FIN.
      endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end subroutine airseaflux_fbk

!______________________________________________________________________

      subroutine read_interp_initial_file
      implicit none
      integer dim_iloc,dim_jloc
#ifdef synopsis
       subroutinetitle='read_interp_initial_file'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Reset dimensions des tableaux de lecture des fichiers:
      meteo_imax=1
      meteo_jmax=1
      call allocate_forcages(1,2,meteo_imax,meteo_jmax,0) ! arg1=allouer arg2=meteo arg3,4,5=dimensions

! Lire et interpoller le fichier meteo e l'etat initial:

!......................................................................
! A l'etat initial il faut commencer par lire les listes ascii et faire
! les listes binaires record:

      if(par%rank==0) then !0000000000000> !13-10-10

      if(par%rank==0)write(6,*)'calcul liste binaire fichiers meteo:'
       do k1=1,k4 ! debut de boucle sur K1


        nc=1
        nc1=nc
        open(unit=3,file=airseafile(k1))
        write(texte30,'(a,i1)')trim(tmpdirname)//'listemeteo',k1             !05-11-09
        open(unit=4,file=texte30                                        &
                   ,access='direct'                                     &
                   ,recl=2*80*8                                         &
                   ,form='unformatted')
  131    read(3,'(a)',end=132)texte80(1)
         read(3,'(a)'        )texte80(2)
         write(4,rec=nc)texte80(1),texte80(2)
         nc1=max0(nc1,nc)
         nc=nc+1
         goto 131

  132   close(4)


        close(3)
       enddo     ! fin de boucle sur K1
       if(par%rank==0)write(6,*)'ok.'

      endif                !0000000000000>

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)  !13-10-10
#endif

!......................................................................


      do k1=1,k4    ! debut de boucle sur K1

       write(texte30,'(a,i1)')trim(tmpdirname)//'listemeteo',k1
        open(unit=4,file=texte30                                        &
                   ,access='direct'                                     &
                   ,recl=2*80*8                                         &
                   ,form='unformatted')

       do k3=0,2,2  ! debut de boucle sur K3

! Remplacees par: le:                                                  !22/01/08
        x2=(elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)
        x1=(elapsedtime_bef-airseadt(k1,2))/airseadt(k1,1)
        j2=int(x2)
        j1=int(x1)
        if(x2.lt.0.)j2=j2-1
        if(x1.lt.0.)j1=j1-1
        nc=int(1.+x2)+k3/2
        if(k3.eq.2.and.(j2-j1).eq.1)nc=nc-1
        nc=max0(nc,ncmin_airsea)

        if(nc.le.0)then !---- debug ----->                             !06/05/04
         if(par%rank==0)write(6,*)'Erreur l 795'
         if(par%rank==0)write(6,*)'je suis bloque dans aiseaflux_fbk choix 1'
         if(par%rank==0)write(6,*)'e la variable ne',k1
         if(par%rank==0)write(6,*)'car nc=',nc,' est <= 0 ce qui signifie que je suis'
         if(par%rank==0)write(6,*)'en avance sur la premiere echeance disponible.'
         if(par%rank==0)write(6,*)'pour continuer je fixe arbitrairement nc e 1'
         ncmin_airsea=1                                                !06/05/04
         nc=max0(nc,ncmin_airsea)                                      !06/05/04
         !pause!11-02-20
        endif           !---- debug ----->

        if(par%rank==0)write(6,*)'airseaflux_fbk choix 1 variable ',k1                 &
                 ,'echeance:',k3,' nc=',nc                              &
                 ,' ncdeci=',1.+x2                                     !02/02/08
        read(4,rec=nc)texte80(1),texte80(2)

! LECTURE DU FICHIER header (entete)
        if(par%rank==0)write(6,'(a,a60)')'sur le point de lire:',texte80(2)
        open(unit=3,file=texte80(2))
         read(3,*)
         read(3,*)dim_iloc       & ! METEO_IMAX
                 ,dim_jloc       & ! METEO_JMAX
                 ,meteo_latmin                                          &
                 ,meteo_lonmin                                          &
                 ,meteo_latmax                                          &
                 ,meteo_lonmax                                          &
                 ,x1,x2
        close(3)
        if(par%rank==0)write(6,*)'ok!'

!......
! Verifier dimensions:
         if(dim_iloc.gt.meteo_imax)then
         meteo_imax=dim_iloc
         call allocate_forcages(2,2,0,0,0) ! arg1=desallouer arg2=meteo
         call allocate_forcages(1,2,meteo_imax,meteo_jmax,0) ! arg1=allouer arg2=meteo arg3,4,5=dimensions
         endif
         if(dim_jloc.gt.meteo_jmax)then
         meteo_jmax=dim_jloc
         call allocate_forcages(2,2,0,0,0) ! arg1=desallouer arg2=meteo
         call allocate_forcages(1,2,meteo_imax,meteo_jmax,0) ! arg1=allouer arg2=meteo arg3,4,5=dimensions
         endif
!......

! LECTURE DU FICHIER BINAIRE RECORD de la VARIABLE
         if(par%rank==0)write(6,'(a,a60)')'sur le point de lire:',texte80(1)
         if(par%rank==0)write(6,*)'avec taille record=',4*meteo_imax*meteo_jmax
         open(unit=3,file=texte80(1)                                    &
                    ,access='direct'                                    &
                    ,recl=4*meteo_imax*meteo_jmax                       &
                    ,form='unformatted')
         read(3,rec=1)((meteo_var(i,j),i=1,meteo_imax),j=1,meteo_jmax)
         close(3)
         if(par%rank==0)write(6,*)'ok!'

            if(k1.eq.1)then !...........>
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

               ssr_w(i,j,k3)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1)

              enddo
              enddo
            endif           !...........>

            if(k1.eq.2)then !...........>
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

              snsf_w(i,j,k3)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1)

              enddo
              enddo
            endif           !...........>


            if(k1.eq.3)then !...........>
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

                 xy_t(i,j,1)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1)

              enddo
              enddo
            endif           !...........>

            if(k1.eq.4)then !...........>
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

                 xy_t(i,j,2)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1)

               uwind_t(i,j,k3)=                                         &
                               gridrotcos_t(i,j)*xy_t(i,j,1)                  &
                              -gridrotsin_t(i,j)*xy_t(i,j,2)

               vwind_t(i,j,k3)=                                         &
                               gridrotsin_t(i,j)*xy_t(i,j,1)                  &
                              +gridrotcos_t(i,j)*xy_t(i,j,2)

              enddo
              enddo
            endif           !...........>

            if(k1.eq.5)then !...........>
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

              pss_w(i,j,k3)=max(small1,                                 &
                             (1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1) )

              enddo
              enddo
            pss_mean(k3)=101300. ! pression de reference constante     !26/03/07
            endif           !...........>

            if(k1.eq.6)then !...........>
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

             teta2_t(i,j,k3)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1)

              enddo
              enddo
            endif           !...........>

            if(k1.eq.7)then !...........>
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

                q2_t(i,j,k3)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1)

              enddo
              enddo
            endif            !...........>

            if(k1.eq.8)then !...........>
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

           precipi_w(i,j,k3)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1)

              enddo
              enddo
            endif            !...........>

       enddo        ! fin de boucle sur K3

      close(4)

      enddo         ! fin de boucle sur K1

! Pour verif:
!     DO J=5,15
!     DO I=5,15
!     WRITE(66,'(16(D12.6,1X))')
!    &         SSR_Z(I,J,0),
!    &         SSR_Z(I,J,2),
!    &        SNSF_Z(I,J,0),
!    &        SNSF_Z(I,J,2),
!    &         UWIND_Z(I,J,0),
!    &         UWIND_Z(I,J,2),
!    &         VWIND_Z(I,J,0),
!    &         VWIND_Z(I,J,2),
!    &        PSS_Z(I,J,0),
!    &        PSS_Z(I,J,2),
!    &       TETA2_Z(I,J,0),
!    &       TETA2_Z(I,J,2),
!    &          Q2_Z(I,J,0),
!    &          Q2_Z(I,J,2),
!    &     PRECIPI_Z(I,J,0),
!    &     PRECIPI_Z(I,J,2)
!     ENDDO
!     ENDDO
!     stop 'verif1'

      end subroutine read_interp_initial_file

!______________________________________________________________________________

      subroutine read_initial_file
      implicit none
#ifdef synopsis
       subroutinetitle='read_initial_file'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Lire des fichiers meteo deje interpolles sur la grille de symphonie:

      do 1100 k1=1,k4                                                  !16/06/06
       do 2000 k3=0,2,2

! Commentees le                                                        !22/01/08
!cc     KOUNTMOD=AIRSEADT(K1,1) ! echeance
!cc     K2      =AIRSEADT(K1,2) ! repere NC=1
!cc     NC=INT( 1.+REAL(KOUNT-K2)/REAL(KOUNTMOD)+ K3/2 )
!cc     IF(MOD(KOUNT-K2,KOUNTMOD).EQ.0.AND.K3.EQ.2)NC=NC-1             !19/03/02
!cc     NC=MAX0(NC,NCMIN_AIRSEA)
! Remplacees par: le:                                                  !22/01/08
        x2=(elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)
        x1=(elapsedtime_bef-airseadt(k1,2))/airseadt(k1,1)
        j2=int(x2)
        j1=int(x1)
        if(x2.lt.0.)j2=j2-1
        if(x1.lt.0.)j1=j1-1
        nc=int(1.+x2)+k3/2
        if(k3.eq.2.and.(j2-j1).eq.1)nc=nc-1
        nc=max0(nc,ncmin_airsea)

        if(nc.le.0)then !---- debug ----->                             !06/05/04
         if(par%rank==0)write(6,*)'Erreur l 1129'
         if(par%rank==0)write(6,*)'je suis bloque dans aiseaflux_fbk choix 1'
         if(par%rank==0)write(6,*)'e la variable ne',k1
         if(par%rank==0)write(6,*)'car nc=',nc,' est <= 0 ce qui signifie que je suis'
         if(par%rank==0)write(6,*)'en avance sur la premiere echeance disponible.'
         if(par%rank==0)write(6,*)'pour continuer je fixe arbitrairement nc e 1'
         ncmin_airsea=1                                                !06/05/04
         nc=max0(nc,ncmin_airsea)                                      !06/05/04
         !pause!11-02-20
        endif           !---- debug ----->
        if(par%rank==0)write(6,*)'airseaflux_fbk choix 1 variable ',k1                 &
                 ,'echeance:',k3,' nc=',nc                              &
                 ,' ncdeci=',1.+x2                                     !02/02/08

        if(par%rank==0)write(6,*)'sur le point de lire:'                              !07/03/05
        if(par%rank==0)write(6,*)airseafile(k1)

        lrec=(imax+2)*(jmax+2)*4
        open(unit=3,file=airseafile(k1),access='direct',recl=lrec       &
         ,form='unformatted')
        read(3,rec=nc)anyvar2d
        close(3)

        if(par%rank==0)write(6,*)'lecture ok'                                         !07/03/05

! Commentees le                                                        !22/01/08
!       X1= AIRSEAINFO(K1,2)*REAL(KOUNTMOD)*DTI_FW
! Remplacees par: le:                                                  !22/01/08
        x1= airseainfo(k1,2)*airseadt(k1,1)                          & !30-04-11
       +(1.-airseainfo(k1,2))*1.

!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.1)then
              do j=0,jmax+1
              do i=0,imax+1
               ssr_w(i,j,k3)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.2)then
              do j=0,jmax+1
              do i=0,imax+1
               snsf_w(i,j,k3)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.3)then
              do j=0,jmax+1
              do i=0,imax+1
               xy_t(i,j,1)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.4)then
              do j=0,jmax+1
              do i=0,imax+1
               xy_t(i,j,2)=anyvar2d(i,j)/x1
              enddo
              enddo
              do j=1,jmax+1
              do i=1,imax+1
               uwind_t(i,j,k3)=                                         &
                               gridrotcos_t(i,j)*xy_t(i,j,1)                  &
                              -gridrotsin_t(i,j)*xy_t(i,j,2)

               vwind_t(i,j,k3)=                                         &
                               gridrotsin_t(i,j)*xy_t(i,j,1)                  &
                              +gridrotcos_t(i,j)*xy_t(i,j,2)
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.5)then
!cc           SUM1=0.
!cc           SUM2=0.
              do j=0,jmax+1
              do i=0,imax+1
               pss_w(i,j,k3)=max(anyvar2d(i,j)/x1,small1)              !29/01/03
!cc            SUM1=SUM1+MORPHO_Z(I,J,NR)
!cc            SUM2=SUM2+MORPHO_Z(I,J,NR)*PSS_Z(I,J,K3)
              enddo
              enddo
!cc         PSS_MEAN(K3)=SUM2/MAX(SUM1,1.D0) ! pression atms moyenne   !01/11/04
            pss_mean(k3)=101300. ! pression de reference constante     !26/03/07
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.6)then
              do j=0,jmax+1
              do i=0,imax+1
               teta2_t(i,j,k3)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.7)then
              do j=0,jmax+1
              do i=0,imax+1
               q2_t(i,j,k3)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.8)then                                            !19/03/03
              do j=0,jmax+1
              do i=0,imax+1
               precipi_w(i,j,k3)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif


 2000 continue
 1100 continue


      end subroutine read_initial_file

!______________________________________________________________________________

      subroutine update_read_file
      implicit none
#ifdef synopsis
       subroutinetitle='update_read_file'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! BOUCLE 1000: remise e jour des echeances:
! DEBUT:
      do 1000 k1=1,k4                                                  !16/06/06
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! Commentees le                                                        !22/01/08
!cccc KOUNTMOD=AIRSEADT(K1,1) ! echeance
!cccc K2      =AIRSEADT(K1,2) ! repere NC=1
!cccc IF( MOD(KOUNT-K2,KOUNTMOD).EQ.0 ) THEN
! Remplacees par: le:                                                  !22/01/08
      x2=(elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)
      x1=(elapsedtime_bef-airseadt(k1,2))/airseadt(k1,1)
      j2=int(x2)
      j1=int(x1)
      if(x2.lt.0.)j2=j2-1
      if(x1.lt.0.)j1=j1-1
      if(j2-j1.eq.1)  then                            !*******************>
! NC est le record au temps present ! NC+1 record pour echeance suivante:
      nc=max0(int(1.+x2),ncmin_airsea-1)

! Commentees le                                                        !22/01/08
! NC   enregistrement present
!cc     NC=MAX0(1+(KOUNT-K2)/KOUNTMOD
!cc  &         ,NCMIN_AIRSEA-1)                                        !06/05/04

        if(nc+1.le.0)then !---- debug ----->                           !06/05/04
        if(par%rank==0)write(6,*)'Erreur l 1281'
        if(par%rank==0)write(6,*)'je suis bloque dans aiseaflux_fbk choix 2'
        if(par%rank==0)write(6,*)'e la variable ne',k1
        if(par%rank==0)write(6,*)'car nc+1=',nc+1,' est <= 0 qui signifie que je suis'
        if(par%rank==0)write(6,*)'en avance sur la premiere echeance disponible.'
        if(par%rank==0)write(6,*)'pour continuer je fixe arbitrairement nc+1 e 1'
        ncmin_airsea=1
!cccccccNC=MAX0(1+(KOUNT-K2)/KOUNTMOD
!cccc&         ,NCMIN_AIRSEA-1)
        nc=max0(int(1.+x2),ncmin_airsea-1)                             !22/01/08
        !pause!11-02-20
        endif             !---- debug ----->
        if(par%rank==0)write(6,*)'airseaflux_fbk choix 2 variable ',k1                 &
                 ,'echeance:',i2,' nc+1=',nc+1

        lrec=(imax+2)*(jmax+2)*4
        open(unit=3,file=airseafile(k1),access='direct',recl=lrec       &
         ,form='unformatted')
! NC+1 enregistrement pour dans la prochaine echeance
        read(3,rec=nc+1)anyvar2d
        close(3)

! Commentees le                                                        !22/01/08
!       X1= AIRSEAINFO(K1,2)*REAL(KOUNTMOD)*DTI_FW
! Remplacees par: le:                                                  !22/01/08
        x1= airseainfo(k1,2)*airseadt(k1,1)                      &  !30-04-11
       +(1.-airseainfo(k1,2))*1.


!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.1)then
              do j=0,jmax+1
              do i=0,imax+1
               ssr_w(i,j,i0)=ssr_w(i,j,i2)
               ssr_w(i,j,i2)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.2)then
              do j=0,jmax+1
              do i=0,imax+1
               snsf_w(i,j,i0)=snsf_w(i,j,i2)
               snsf_w(i,j,i2)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.3)then
              do j=0,jmax+1
              do i=0,imax+1
               xy_t(i,j,1)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.4)then
              do j=0,jmax+1
              do i=0,imax+1
               xy_t(i,j,2)=anyvar2d(i,j)/x1
              enddo
              enddo
              do j=1,jmax+1
              do i=1,imax+1
               uwind_t(i,j,i0)=uwind_t(i,j,i2)
               vwind_t(i,j,i0)=vwind_t(i,j,i2)

               uwind_t(i,j,i2)=                                         &
                               gridrotcos_t(i,j)*xy_t(i,j,1)                  &
                              -gridrotsin_t(i,j)*xy_t(i,j,2)

               vwind_t(i,j,i2)=                                         &
                               gridrotsin_t(i,j)*xy_t(i,j,1)                  &
                              +gridrotcos_t(i,j)*xy_t(i,j,2)
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.5)then
!ccc          SUM1=0.
!ccc          SUM2=0.
              pss_mean(i0)=pss_mean(i2)
              do j=0,jmax+1
              do i=0,imax+1
               pss_w(i,j,i0)=pss_w(i,j,i2)
               pss_w(i,j,i2)=max(anyvar2d(i,j)/x1,small1)              !29/01/03
!ccc           SUM1=SUM1+MORPHO_Z(I,J,NR)
!ccc           SUM2=SUM2+MORPHO_Z(I,J,NR)*PSS_Z(I,J,I2)
              enddo
              enddo
            endif
!ccc        PSS_MEAN(I2)=SUM2/MAX(SUM1,1.D0) ! pression atms moyenne   !01/11/04
            pss_mean(i2)=101300. ! pression de reference constante     !26/03/07
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.6)then
              do j=0,jmax+1
              do i=0,imax+1
               teta2_t(i,j,i0)=teta2_t(i,j,i2)
               teta2_t(i,j,i2)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.7)then
              do j=0,jmax+1
              do i=0,imax+1
               q2_t(i,j,i0)=q2_t(i,j,i2)
               q2_t(i,j,i2)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !
            if(k1.eq.8)then                                            !19/03/03
              do j=0,jmax+1
              do i=0,imax+1
               precipi_w(i,j,i0)=precipi_w(i,j,i2)
               precipi_w(i,j,i2)=anyvar2d(i,j)/x1
              enddo
              enddo
            endif
!    .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   !

      endif                                             !*******************>

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! BOUCLE 1000: remise e jour des echeances:
! FIN.
 1000 continue
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end subroutine update_read_file

!______________________________________________________________________________

      subroutine update_read_interp_file
      implicit none
      integer dim_iloc,dim_jloc
#ifdef synopsis
       subroutinetitle='update_read_interp_file'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! BOUCLE 1000: remise e jour des echeances:
! DEBUT:
      do 1000 k1=1,k4                                                  !16/06/06
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      x2=(elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)
      x1=(elapsedtime_bef-airseadt(k1,2))/airseadt(k1,1)
      j2=int(x2)
      j1=int(x1)
      if(x2.lt.0.)j2=j2-1
      if(x1.lt.0.)j1=j1-1
      if(j2-j1.eq.1)  then                            !*******************>

! NC est le record au temps present ! NC+1 record pour echeance suivante:
      nc=max0(int(1.+x2),ncmin_airsea-1)

        if(nc+1.le.0)then !---- debug ----->                           !06/05/04
        if(par%rank==0)write(6,*)'Erreur 1437'
        if(par%rank==0)write(6,*)'je suis bloque dans aiseaflux_fbk choix 2'
        if(par%rank==0)write(6,*)'e la variable ne',k1
        if(par%rank==0)write(6,*)'car nc+1=',nc+1,' est <= 0 qui signifie que je suis'
        if(par%rank==0)write(6,*)'en avance sur la premiere echeance disponible.'
        if(par%rank==0)write(6,*)'pour continuer je fixe arbitrairement nc+1 e 1'
        ncmin_airsea=1
        nc=max0(int(1.+x2),ncmin_airsea-1)                             !22/01/08
        !pause!11-02-20
        endif             !---- debug ----->
        if(par%rank==0)write(6,*)'airseaflux_fbk choix 2 variable ',k1                 &
                 ,'echeance:',i2,' nc+1=',nc+1

       write(texte30,'(a,i1)')trim(tmpdirname)//'listemeteo',k1
        open(unit=4,file=texte30                                        &
                   ,access='direct'                                     &
                   ,recl=2*80*8                                         &
                   ,form='unformatted')
        read(4,rec=nc+1)texte80(1),texte80(2)
        close(4)
!       if(par%rank==0)write(6,'(a)')texte80(2)

        open(unit=3,file=texte80(2))
         read(3,*)
         read(3,*)dim_iloc       & ! METEO_IMAX
                 ,dim_jloc       & ! METEO_JMAX
                 ,meteo_latmin                                          &
                 ,meteo_lonmin                                          &
                 ,meteo_latmax                                          &
                 ,meteo_lonmax                                          &
                 ,x1,x2
        close(3)

!......
! Verifier dimensions:
         if(dim_iloc.gt.meteo_imax)then
         meteo_imax=dim_iloc
         call allocate_forcages(2,2,0,0,0) ! arg1=desallouer arg2=meteo
         call allocate_forcages(1,2,meteo_imax,meteo_jmax,0) ! arg1=allouer arg2=meteo arg3,4,5=dimensions
         endif
         if(dim_jloc.gt.meteo_jmax)then
         meteo_jmax=dim_jloc
         call allocate_forcages(2,2,0,0,0) ! arg1=desallouer arg2=meteo
         call allocate_forcages(1,2,meteo_imax,meteo_jmax,0) ! arg1=allouer arg2=meteo arg3,4,5=dimensions
         endif
!......

         open(unit=3,file=texte80(1)                                    &
                    ,access='direct'                                    &
                    ,recl=4*meteo_imax*meteo_jmax                       &
                    ,form='unformatted')
         read(3,rec=1)((meteo_var(i,j),i=1,meteo_imax),j=1,meteo_jmax)
         close(3)

            if(k1.eq.1)then !...........>
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

               ssr_w(i,j,i0)=ssr_w(i,j,i2)
               ssr_w(i,j,i2)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1)

              enddo
              enddo
            endif           !...........>

            if(k1.eq.2)then !...........>
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

              snsf_w(i,j,i0)=snsf_w(i,j,i2)
              snsf_w(i,j,i2)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1)

              enddo
              enddo
            endif           !...........>


            if(k1.eq.3)then !...........>
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

                 xy_t(i,j,1)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1)

              enddo
              enddo
            endif           !...........>

            if(k1.eq.4)then !...........>
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

                 xy_t(i,j,2)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1)

               uwind_t(i,j,i0)=uwind_t(i,j,i2)
               uwind_t(i,j,i2)=                                         &
                               gridrotcos_t(i,j)*xy_t(i,j,1)                  &
                              -gridrotsin_t(i,j)*xy_t(i,j,2)

               vwind_t(i,j,i0)=vwind_t(i,j,i2)
               vwind_t(i,j,i2)=                                         &
                               gridrotsin_t(i,j)*xy_t(i,j,1)                  &
                              +gridrotcos_t(i,j)*xy_t(i,j,2)

              enddo
              enddo
            endif           !...........>

            if(k1.eq.5)then !...........>
              pss_mean(i0)=pss_mean(i2)
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

              pss_w(i,j,i0)=pss_w(i,j,i2)
              pss_w(i,j,i2)=max(small1,                                 &
                             (1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1) )

              enddo
              enddo
            pss_mean(i2)=101300. ! pression de reference constante     !26/03/07
            endif           !...........>

            if(k1.eq.6)then !...........>
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

             teta2_t(i,j,i0)=teta2_t(i,j,i2)
             teta2_t(i,j,i2)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1)

              enddo
              enddo
            endif           !...........>

            if(k1.eq.7)then !...........>
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

                q2_t(i,j,i0)=q2_t(i,j,i2)
                q2_t(i,j,i2)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1)

              enddo
              enddo
            endif            !...........>

            if(k1.eq.8)then !...........>
              do j=0,jmax+1
              do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
              deci=1.+( lon_t(i,j)*180./pi-meteo_lonmin)                &
                            /(meteo_lonmax-meteo_lonmin)                &
                            *(meteo_imax-1.)
              i1=min0(max0(int(deci),1),meteo_imax-1)
              rapi=deci-i1 ! =0 si DECI=I1 =1 si DECI=I1+1
! Indice decimale j (latitude) dans grille aladin:
              decj=1.+( lat_t(i,j)*180./pi-meteo_latmin)                &
                            /(meteo_latmax-meteo_latmin)                &
                            *(meteo_jmax-1.)
              j1=min0(max0(int(decj),1),meteo_jmax-1)
              rapj=decj-j1 !=0 si DECJ=J1 =1 si DECJ=J1+1

           precipi_w(i,j,i0)=precipi_w(i,j,i2)
           precipi_w(i,j,i2)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                            +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                            +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                            +    rapi *    rapj *meteo_var(i1+1,j1+1)

              enddo
              enddo
            endif            !...........>

      endif                                             !*******************>

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! BOUCLE 1000: remise e jour des echeances:
! FIN.
 1000 continue
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end subroutine update_read_interp_file

!----------------------------------------------------------------------------

      subroutine airseaflux_initial
      implicit none
#ifdef synopsis
       subroutinetitle='airseaflux_initial'
       subroutinedescription= &
       'Driver for subroutines involved in initialsteps'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Idendifier le type de fichier netcdf e l'aide du nom de la liste des fichiers meteo
!     call airseaflux_model_flag

! Initialiser les identifiants (type integer) des variables meteo
      if(par%rank==0)write(6,*)&
      'airseaflux_initial call airseaflux_varidentifier'
      call airseaflux_varidentifier

! A l'etat initial, faire la liste binaire des fichiers meteo
      if(par%rank==0)write(6,*)&
      'airseaflux_initial call airseaflux_binary_file_list'
      call airseaflux_binary_file_list

! Devra t'on calculer l'albedo ?
      if(par%rank==0)write(6,*)&
      'airseaflux_initial call airseaflux_inquire_albedo'
      call airseaflux_inquire_albedo('initial') ! ialbedo

! Devra t'on deduire l'humidite specifique du dew point?
      if(par%rank==0)write(6,*)&
      'airseaflux_initial call airseaflux_sphum_or_dewp'
      call airseaflux_sphum_or_dewp('initial')

! Le flux infra rouge est t'il le bilan atmosphere/ocean (flag_net_ir=1) ou seulement le
! flux atmospherique (flag_net_ir=0):
      call airseaflux_inquire_longwave('initial') !28-11-16

! Reduire l'extraction des donnees a une zone englobant notre domaine numerique
! et definir les bouchetrous
      if(par%rank==0)write(6,*)&
      'airseaflux_initial call airseaflux_extract_zone'
      call airseaflux_extract_zone('initial')

      end subroutine airseaflux_initial

!----------------------------------------------------------------------------

      subroutine read_interp_netcdf_airseafile(t_)
      implicit none
      double precision scalecum_,time_,timearray_(1)
      integer dim_iloc,dim_jloc,dim_kloc,t_,loop_,ncid_,time1forecast_
#ifdef synopsis
       subroutinetitle='read_interp_netcdf_airseafile'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0)then
       write(6,*)'-------------------------------------------'
       write(6,*)'read_interp_netcdf_airseafile t_=',t_
      endif

! . . .
      call time_for_a_new_file(t_,3)!?  3 est le numero du vent dans notebook_airseaflux
      if(decision==1) then                            !*******************>    !15-04-11
! ouvrir le fichier netcdf:
      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(par%rank==0) then
       write(6,*)
       write(6,'(a,a)')' FILE:  ',trim(texte80(1))
      endif
      if(status.ne.0)stop 'echec fichier 2'
      time1forecast_=nint(airseadt(3,1)/3600.) ! time first forecast in hours

! lire la composante W->E du vent:
      ksecu=0
      if(flag_meteodata=='slimga') then
      ksecu=1
      varstart(4)=nc1            ; varcount(4)=1          ! time
      varstart(3)=1              ; varcount(3)=1          ! height_4
      stop 'cas slim e verifier pour passsage de varstart(4)'
      endif
      if(flag_meteodata=='ecmwf'.or.              &
         flag_meteodata=='arome') then
      ksecu=1
      varstart(3)=nc1            ; varcount(3)=1          ! time
      if(par%rank==0)write(6,*)'WE wind varstart     ',varstart(3)
      endif
      if(ksecu==0)stop 'meteo data file not understood 2'

      varstart(2)=meteozoom_jstr ; varcount(2)=meteo_jmax ! lat
      varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax ! lon
                   status=nf_inq_varid(ncid1,'u-wind',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'U10_GDS0_SFC',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'U10M',var_id)
      if(status/=0)stop 'echec id Vent WE'
      status=nf_get_vara_real(ncid1,var_id,varstart(1:3),varcount(1:3),meteo_var)     ! 25-01-11 ajout (1:3)
      if(status.ne.0)stop 'echec lecture u-wind'

! Bouchage des trous:
      call appli_bouchetrou_meteo

! Interpoler la composante W->E du vent:
      i0=0 ; i2=2
      const1=1./meteo_londlt                              !26-08-10
      const2=1./meteo_latdlt                              !26-08-10
      do j=0,jmax+1
      do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
      deci=1.+( lon_t(i,j)*rad2deg-meteo_lonstr)*const1
      i1=int(deci)
      rapi=deci-i1
      i1=i1-meteozoom_istr+1                              ! conservation la parallelisation

! Indice decimale j (latitude) dans grille aladin:
      decj=1.+( lat_t(i,j)*rad2deg-meteo_latstr)*const2
      j1=int(decj)
      rapj=decj-j1
      j1=j1-meteozoom_jstr+1                              ! conservation la parallelisation
!     i1=min(max(1,i1),meteo_imax)
!     j1=min(max(1,j1),meteo_jmax)

      xy_t(i,j,1)=  (1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                   +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                   +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                   +    rapi *    rapj *meteo_var(i1+1,j1+1)

      enddo
      enddo

! . . .

! lire la composante S->N du vent:
      ksecu=0
      if(flag_meteodata=='slimga') then
      ksecu=1
      varstart(4)=nc1            ; varcount(4)=1          ! time
      varstart(3)=1              ; varcount(3)=1          ! height_4
      endif
      if(flag_meteodata=='ecmwf'.or.              &
         flag_meteodata=='arome') then
      ksecu=1
      varstart(3)=nc1            ; varcount(3)=1          ! time
      if(par%rank==0)write(6,*)'SN wind varstart     ',varstart(3)
      endif
      if(ksecu==0)stop 'meteo data file not understood 3'
      varstart(2)=meteozoom_jstr ; varcount(2)=meteo_jmax ! lat
      varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax ! lon
                   status=nf_inq_varid(ncid1,'v-wind',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'V10_GDS0_SFC',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'V10M',var_id)
      if(status/=0)stop 'echec id Vent SN'
      status=nf_get_vara_real(ncid1,var_id,varstart(1:3),varcount(1:3),meteo_var)
      if(status.ne.0)stop 'echec lecture v-wind'

! Bouchage des trous:
      call appli_bouchetrou_meteo

! Interpoler la composante S->N du vent:
      i0=0 ; i2=2
      const1=1./meteo_londlt                              !26-08-10
      const2=1./meteo_latdlt                              !26-08-10
      do j=0,jmax+1
      do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
      deci=1.+( lon_t(i,j)*rad2deg-meteo_lonstr)*const1
      i1=int(deci)
      rapi=deci-i1
      i1=i1-meteozoom_istr+1                              ! conservation la parallelisation

! Indice decimale j (latitude) dans grille aladin:
      decj=1.+( lat_t(i,j)*rad2deg-meteo_latstr)*const2
      j1=int(decj)
      rapj=decj-j1
      j1=j1-meteozoom_jstr+1                              ! conservation la parallelisation
!     i1=min(max(1,i1),meteo_imax)
!     j1=min(max(1,j1),meteo_jmax)

      xy_t(i,j,2)=  (1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                   +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                   +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                   +    rapi *    rapj *meteo_var(i1+1,j1+1)

      uwind_t(i,j,i0)=uwind_t(i,j,i2)
      uwind_t(i,j,i2)=                                       &
                      gridrotcos_t(i,j)*xy_t(i,j,1)          &
                     -gridrotsin_t(i,j)*xy_t(i,j,2)

      vwind_t(i,j,i0)=vwind_t(i,j,i2)
      vwind_t(i,j,i2)=                                         &
                      gridrotsin_t(i,j)*xy_t(i,j,1)             &
                     +gridrotcos_t(i,j)*xy_t(i,j,2)

      enddo
      enddo


! Fermer le fichier netcdf:
      status=nf_close(ncid1)
      endif                                           !*******************>

! . . .
! lire le rayonnement courtes longueurs d'ondes descendant
      call time_for_a_new_file(t_,1)!?  1 est le numero du flux solaire dans notebook_airseaflux
      if(decision==1) then                            !*******************>    !15-04-11
! ouvrir le fichier netcdf:
      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(par%rank==0)then
       write(6,*)
       write(6,'(a,a)')' FILE:  ',trim(texte80(1))
      endif
! on ouvre aussi le fichier precedent pour le cas particulier ou celui-ci est necessaire pour decumuler
      status=nf_open(trim(texte80(2)),nf_nowrite,ncid2)
      if(status.ne.0)stop 'echec fichier 3'

      time1forecast_=nint(airseadt(1,1)/3600.) ! time first forecast in hours
      if(flag_meteotime==1) then !ttttttttttttttttttttttt>

! Lire le temps (en heures depuis le debut de la previ) de l'echeance:
                   status=nf_inq_varid(ncid1,'time',var_id)
      if(status/=0)stop 'echec id time'
      varstart(1)=nc1  ; varcount(1)=1          ! time
!     status=nf_get_vara_double(ncid1,var_id,varstart(1:1),varcount(1:1),time_)
      status=nf_get_vara_double(ncid1,var_id,varstart(1:1),varcount(1:1)&
                               ,timearray_(1))!26-03-11
      time_=timearray_(1)                  !26-03-11
      if(status/=0)stop 'echec lecture time'

      else                       !ttttttttttttttttttttttt>

       time_=time1forecast_

      endif                      !ttttttttttttttttttttttt>

! Traiter le cumul ou le non cumul:                           !25-08-10
                   status=nf_inq_varid(ncid1,'down-short-flux',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'SSRD',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'SO_GDS0_SFC',var_id)
      if(status/=0)stop 'echec id down-short-flux'
      texte90=''
      scalecum_=1.

      status=nf_get_att_text(ncid1,var_id,'units',texte90)
      if(status/=0)stop 'echec attribut units pour ssr'

      flag_cumul=0
      loop_=0
      if(texte90=='W m**-2 s'.or.texte90=='J m**-2') then !------>
! les unites indiquent que la quantite est cumulee
      flag_cumul=1
      scalecum_=1./airseadt(1,1) ! 30-04-11
      loop_=1
      if(time_==time1forecast_)loop_=0     ! sauf dans le cas de la premiere echeance de chaque fichier (k0=0)
      endif                         !------>

      do loop1=loop_,0,-1  ! 2 passages (loop1=1 puis loop1=0) si variables cumulee

      ksecu=0
      if(flag_meteodata=='slimga') then
      ksecu=1
      varstart(4)=nc1-loop1      ; varcount(4)=1          ! time
      varstart(3)=1              ; varcount(3)=1          ! height_4
      endif

      if(flag_meteodata=='ecmwf'.or.              &
         flag_meteodata=='arome') then

         if(time_==time1forecast_.and.nc1/=1)     &
         stop 'time_==time1forecast_.and.nc1/=1'

         ksecu=1
         varstart(3)=nc1-loop1      ; varcount(3)=1          ! time
         ncid_=ncid1

! Cas particulier: Premiere echeance d'un fichier qui n'est pas la premiere echeance de previ
! il faut decumuler par rapport e la derniere echeance du jour (fichier) precedent
! varstart(3)=dernier echeance du precedent fichier soit nint(airseainfo(1,3))
         if(flag_cumul==1) then !cumulcumulcumul>
          if(time_/=time1forecast_.and.nc1==1) then !cpcpcpcpcp>
           ksecu=1
           varstart(3)=max_meteo_time_counter  ; varcount(3)=1          ! time
           ncid_=ncid2
          endif                                     !cpcpcpcpcp>
         endif                  !cumulcumulcumul>

         if(par%rank==0)then
           write(6,*)'SSR varstart & loop1 ',varstart(3),loop1
           if(ncid_==ncid2)write(6,'(a,a)')'Previous file:',texte80(2)
         endif

      endif                              !11111111111>

      if(ksecu==0)stop 'meteo data file not understood'
      varstart(2)=meteozoom_jstr ; varcount(2)=meteo_jmax ! lat
      varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax ! lon
      status=nf_get_vara_real(ncid_,var_id,varstart(1:3),varcount(1:3),meteo_var)
      if(status/=0)stop 'echec lecture down-short-flux'

      if(scalecum_/=1.)meteo_var=meteo_var*scalecum_
      if(loop_==1.and.loop1==1)meteo_cum=meteo_var
      if(loop_==1.and.loop1==0)meteo_var=meteo_var-meteo_cum

      enddo                  ! end of loop on "loop1"

! Bouchage des trous:
      call appli_bouchetrou_meteo

! Interpoler le rayonnement courtes longueurs d'ondes descendant
      i0=0 ; i2=2
      const1=1./meteo_londlt                              !26-08-10
      const2=1./meteo_latdlt                              !26-08-10
      do j=0,jmax+1
      do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
      deci=1.+( lon_t(i,j)*rad2deg-meteo_lonstr)*const1
      i1=int(deci)
      rapi=deci-i1
      i1=i1-meteozoom_istr+1                              ! conservation la parallelisation

! Indice decimale j (latitude) dans grille aladin:
      decj=1.+( lat_t(i,j)*rad2deg-meteo_latstr)*const2
      j1=int(decj)
      rapj=decj-j1
      j1=j1-meteozoom_jstr+1                              ! conservation la parallelisation
!     i1=min(max(1,i1),meteo_imax)
!     j1=min(max(1,j1),meteo_jmax)

      ssr_w(i,j,i0)=ssr_w(i,j,i2)

      ssr_w(i,j,i2)= (1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )         &
                    +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)         &
                    +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )         &
                    +    rapi *    rapj *meteo_var(i1+1,j1+1)

      enddo
      enddo

! Fermer le fichier netcdf:
      status=nf_close(ncid1)
      status=nf_close(ncid2)
      endif                                           !*******************>

! . . .
! lire le rayonnement grandes longueurs d'ondes descendant
      call time_for_a_new_file(t_,2)!?  2 est le numero du flux thermique dans notebook_airseaflux
      if(decision==1) then                            !*******************>    !15-04-11
! ouvrir le fichier netcdf:
      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(par%rank==0) then
       write(6,*)
       write(6,'(a,a)')' FILE:  ',trim(texte80(1))
      endif
! on ouvre aussi le fichier precedent pour le cas particulier ou celui-ci est necessaire pour decumuler
      status=nf_open(trim(texte80(2)),nf_nowrite,ncid2)
      if(status.ne.0)stop 'echec fichier 3'

      time1forecast_=nint(airseadt(2,1)/3600.) ! time first forecast in hours
      if(flag_meteotime==1) then !ttttttttttttttttttttttt>

! Lire le temps (en heures depuis le debut de la previ) de l'echeance:
                   status=nf_inq_varid(ncid1,'time',var_id)
      if(status/=0)stop 'echec id time'
      varstart(1)=nc1  ; varcount(1)=1          ! time
!     status=nf_get_vara_double(ncid1,var_id,varstart(1:1),varcount(1:1),time_)
      status=nf_get_vara_double(ncid1,var_id,varstart(1:1),varcount(1:1),timearray_(1))!26-03-11
      time_=timearray_(1)              !26-03-11
      if(status/=0)stop 'echec lecture time'
      if(status.ne.0)stop 'echec fichier 4'

      else                      !ttttttttttttttttttttttt>

      time_=time1forecast_

      endif                      !ttttttttttttttttttttttt>

!> Traiter le cumul ou le non cumul:                           !25-08-10
                   status=nf_inq_varid(ncid1,'down-long-flux',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'STRD',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'TH_GDS0_SFC',var_id)
      if(status/=0)stop 'echec id down-long-flux'
      texte90=''
      scalecum_=1.

      status=nf_get_att_text(ncid1,var_id,'units',texte90)
      if(status/=0)stop 'echec attribut units pour snsf'

      flag_cumul=0
      loop_=0
      if(texte90=='W m**-2 s'.or.texte90=='J m**-2') then !------>
! les unites non indiquent que la quantite est cumulee
      flag_cumul=1
      scalecum_=1./airseadt(2,1)   !30-04-11
      loop_=1
      if(time_==time1forecast_)loop_=0     ! sauf dans le cas de la premiere echeance de chaque fichier (k0=0)
      endif                         !------>
!>

      do loop1=loop_,0,-1  ! 2 passages (loop1=1 puis loop1=0) si variables cumulee

      ksecu=0
      if(flag_meteodata=='slimga') then
      ksecu=1
      varstart(4)=nc1-loop1      ; varcount(4)=1          ! time
      varstart(3)=1              ; varcount(3)=1          ! height_4
      endif
      if(flag_meteodata=='ecmwf'.or.              &
         flag_meteodata=='arome') then

         if(time_==time1forecast_.and.nc1/=1)     &
         stop 'time_==time1forecast_.and.nc1/=1'

         ksecu=1
         varstart(3)=nc1-loop1      ; varcount(3)=1          ! time
         ncid_=ncid1

! Cas particulier: Premiere echeance d'un fichier qui n'est pas la premiere echeance de previ
! il faut decumuler par rapport e la derniere echeance du jour (fichier) precedent
! varstart(3)=dernier echeance du precedent fichier soit nint(airseainfo(1,3))
         if(flag_cumul==1) then          !cumulcumulcumul>
          if(time_/=time1forecast_.and.nc1==1) then !cpcpcpcpcp>
           ksecu=1
           varstart(3)=max_meteo_time_counter  ; varcount(3)=1          ! time
           ncid_=ncid2
          endif                                     !cpcpcpcpcp>
         endif                           !cumulcumulcumul>
         if(par%rank==0)then
           write(6,*)'IRD varstart & loop1 ',varstart(3),loop1
           if(ncid_==ncid2)write(6,'(a,a)')'Previous file:',texte80(2)
         endif

      endif                              !11111111111>

      if(ksecu==0)stop 'meteo data file not understood'
      varstart(2)=meteozoom_jstr ; varcount(2)=meteo_jmax ! lat
      varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax ! lon
      status=nf_get_vara_real(ncid_,var_id,varstart(1:3),varcount(1:3),meteo_var)
      if(status/=0)stop 'echec lecture down-long-flux'

      if(scalecum_/=1.)meteo_var=meteo_var*scalecum_
      if(loop_==1.and.loop1==1)meteo_cum=meteo_var
      if(loop_==1.and.loop1==0)meteo_var=meteo_var-meteo_cum

      enddo                  ! end of loop on "loop1"

! Bouchage des trous:
      call appli_bouchetrou_meteo

! Interpoler le rayonnement grandes longueurs d'ondes descendant
      i0=0 ; i2=2
      const1=1./meteo_londlt                              !26-08-10
      const2=1./meteo_latdlt                              !26-08-10
      do j=0,jmax+1
      do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
      deci=1.+( lon_t(i,j)*rad2deg-meteo_lonstr)*const1
      i1=int(deci)
      rapi=deci-i1
      i1=i1-meteozoom_istr+1                              ! conservation la parallelisation

! Indice decimale j (latitude) dans grille aladin:
      decj=1.+( lat_t(i,j)*rad2deg-meteo_latstr)*const2
      j1=int(decj)
      rapj=decj-j1
      j1=j1-meteozoom_jstr+1                              ! conservation la parallelisation
!     i1=min(max(1,i1),meteo_imax)
!     j1=min(max(1,j1),meteo_jmax)

      snsf_w(i,j,i0)=snsf_w(i,j,i2)

      snsf_w(i,j,i2)= (1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )         &
                     +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)         &
                     +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )         &
                     +    rapi *    rapj *meteo_var(i1+1,j1+1)

      enddo
      enddo

! Fermer le fichier netcdf:
      status=nf_close(ncid1)
      status=nf_close(ncid2)        !14-04-11
      endif                                           !*******************>

! . . . .
! lire la pression atmospherique:
      call time_for_a_new_file(t_,5)!?  5 est le numero de la pression dans notebook_airseaflux
      if(decision==1) then                            !*******************>    !15-04-11
! ouvrir le fichier netcdf:
      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(par%rank==0)then
       write(6,*)
       write(6,'(a,a)')' FILE:  ',trim(texte80(1))
      endif
      if(status.ne.0)stop 'echec fichier 5'

      varstart(3)=nc1            ; varcount(3)=1          ! time
      varstart(2)=meteozoom_jstr ; varcount(2)=meteo_jmax ! lat
      varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax ! lon
                   status=nf_inq_varid(ncid1,'pressure',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'SP',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'LNSP_GDS0_HYBL',var_id)
      if(status/=0)stop 'echec id pressure'

      if(par%rank==0)write(6,*)'Pressure varstart    ',varstart(3)

      status=nf_get_vara_real(ncid1,var_id,varstart(1:3),varcount(1:3),meteo_var)
      if(status.ne.0)stop 'echec lecture pressure'

! Bouchage des trous:
      call appli_bouchetrou_meteo

      status=nf_get_att_text(ncid1,var_id,'long_name',texte60)
      if(status/=0)stop 'Erreur long_name surface pressure'
      if(index(texte60,'Logarithm of surface pressure')/=0) then !lnlnlnlnln>
        do j=1,meteo_jmax
        do i=1,meteo_imax
          meteo_var(i,j)=exp(meteo_var(i,j))
        enddo
        enddo
      endif                                                      !lnlnlnlnln>

! Interpoler la pression
      i0=0 ; i2=2
      const1=1./meteo_londlt                              !26-08-10
      const2=1./meteo_latdlt                              !26-08-10
      do j=0,jmax+1
      do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
      deci=1.+( lon_t(i,j)*rad2deg-meteo_lonstr)*const1
      i1=int(deci)
      rapi=deci-i1
      i1=i1-meteozoom_istr+1                              ! conservation la parallelisation

! Indice decimale j (latitude) dans grille aladin:
      decj=1.+( lat_t(i,j)*rad2deg-meteo_latstr)*const2
      j1=int(decj)
      rapj=decj-j1
      j1=j1-meteozoom_jstr+1                              ! conservation la parallelisation
!     i1=min(max(1,i1),meteo_imax)
!     j1=min(max(1,j1),meteo_jmax)

      pss_w(i,j,i0)=pss_w(i,j,i2)
      pss_w(i,j,i2)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                   +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                   +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                   +    rapi *    rapj *meteo_var(i1+1,j1+1)

      enddo
      enddo
      pss_mean(i0)=pss_mean(i2)
      pss_mean(i2)=101300. ! pression de reference constante     !26/03/07

! Fermer le fichier netcdf:
      status=nf_close(ncid1)
      endif                                           !*******************>
! . . .

! lire l'humidite specifique ou dewpoint
      call time_for_a_new_file(t_,7)!?  7 est le numero de hum. spec. dans notebook_airseaflux
      if(decision==1) then                            !*******************>    !15-04-11
! ouvrir le fichier netcdf:
      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(par%rank==0)then
       write(6,*)
       write(6,'(a,a)')' FILE:  ',trim(texte80(1))
      endif
      if(status.ne.0)stop 'echec fichier 6'

      ksecu=0
      if(flag_meteodata=='slimga') then
      ksecu=1
      varstart(4)=nc1            ; varcount(4)=1          ! time
      varstart(3)=1              ; varcount(3)=1          ! height_4
      endif
      if(flag_meteodata=='ecmwf'.or.              &
         flag_meteodata=='arome') then
      ksecu=1
      varstart(3)=nc1            ; varcount(3)=1          ! time
      endif
      if(ksecu==0)stop 'meteo data file not understood'
      varstart(2)=meteozoom_jstr ; varcount(2)=meteo_jmax ! lat
      varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax ! lon

                   status=nf_inq_varid(ncid1,'specific_humidity',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'Q2_GDS0_SFC',var_id)
      if(status==0)flag_humspec_dewpoint='humspeci'
      if(status/=0)then !>>>>>>>      !19-04-12
       status=nf_inq_varid(ncid1,'D2M',var_id)
       if(status==0)flag_humspec_dewpoint='dewpoint'
      endif             !>>>>>>>

      if(par%rank==0)write(6,*)'Q2 or DWP varstart   ',varstart(3)

      if(status/=0)stop 'echec id specific_humidity'
      status=nf_get_vara_real(ncid1,var_id,varstart(1:3),varcount(1:3),meteo_var)
      if(status/=0)stop 'echec lecture specific_humidity ou dew point'

! Bouchage des trous:
      call appli_bouchetrou_meteo

! Interpoler l'humidite specifique ou dewpoint
      i0=0 ; i2=2
      const1=1./meteo_londlt                              !26-08-10
      const2=1./meteo_latdlt                              !26-08-10
      do j=0,jmax+1
      do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
      deci=1.+( lon_t(i,j)*rad2deg-meteo_lonstr)*const1
      i1=int(deci)
      rapi=deci-i1
      i1=i1-meteozoom_istr+1                              ! conservation la parallelisation

! Indice decimale j (latitude) dans grille aladin:
      decj=1.+( lat_t(i,j)*rad2deg-meteo_latstr)*const2
      j1=int(decj)
      rapj=decj-j1
      j1=j1-meteozoom_jstr+1                              ! conservation la parallelisation
!     i1=min(max(1,i1),meteo_imax)
!     j1=min(max(1,j1),meteo_jmax)

      q2_t(i,j,i0)=q2_t(i,j,i2)
      q2_t(i,j,i2)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                  +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                  +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                  +    rapi *    rapj *meteo_var(i1+1,j1+1)

      enddo
      enddo

      if(flag_humspec_dewpoint=='dewpoint')then ! on calcule l'humidite specifique                   !25-01-11 changement formules
!     const5=0.622*6.11*exp(5417.7530/273.16)        !  pression de vapeur saturante fction du pt de rosee
      do j=0,jmax+1
      do i=0,imax+1

      x1=6.1078*10.**((q2_t(i,j,i2)-273.15)*7.5/(q2_t(i,j,i2)-273.15+237.3))  ! pression vapeur saturante fct(point de rosee): unite mb
                                                                              ! verif tableau9 page 120 Queney
      x2=x1*100./pss_w(i,j,i2)      ! e/p avec  passe de mb en Pascal
      q2_t(i,j,i2)=0.622*x2/(1.-0.378*x2)                                     ! expression 52 p115 Queney q=f(e/p)
!     x1 = 6.11 * exp(5417.7530 * ((1/273.16) - (1/q2_t(i,j,i2))))   !  pression de vapeur saturante fction du pt de rosee
!     q2_t(i,j,i2)=           & ! humidite spec
!           0.622*x1          & ! 0.622 * Pression Vapeur Saturante
!           /pss_w(i,j,i2)      !  Pression Surface
! Equivalent e:
!      q2_t(i,j,i2)=const5*exp(-(5417.7530/q2_t(i,j,i2)))/pss_w(i,j,i2)

      enddo
      enddo
      endif

! Fermer le fichier netcdf:
      status=nf_close(ncid1)
      endif                                           !*******************>
! . . .

! lire la temperature e 2 m (egalement niveau 2 de la grille meteo):
      call time_for_a_new_file(t_,6)!?  6 est le numero de T pot dans notebook_airseaflux
      if(decision==1) then                            !*******************>    !15-04-11
! ouvrir le fichier netcdf:
      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(par%rank==0)then
       write(6,*)
       write(6,'(a,a)')' FILE:  ',trim(texte80(1))
      endif
      if(status.ne.0)stop 'echec fichier'

      ksecu=0
      if(flag_meteodata=='slimga') then
      ksecu=1
      varstart(4)=nc1            ; varcount(4)=1          ! time
      varstart(3)=2              ; varcount(3)=1          ! height_4
      endif
      if(flag_meteodata=='ecmwf'.or.              &
         flag_meteodata=='arome') then
      ksecu=1
      varstart(3)=nc1            ; varcount(3)=1          ! time
      if(par%rank==0)write(6,*)'Air Temp varstart    ',varstart(3)
      endif
      if(ksecu==0)stop 'meteo data file not understood'
      varstart(2)=meteozoom_jstr ; varcount(2)=meteo_jmax ! lat
      varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax ! lon
                   status=nf_inq_varid(ncid1,'temperature',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'T2_GDS0_SFC',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'T2M',var_id)
      if(status/=0)stop 'echec id temperature'
      status=nf_get_vara_real(ncid1,var_id,varstart(1:3),varcount(1:3),meteo_var)
      if(status/=0)stop 'echec lecture temperature'

! Bouchage des trous:
      call appli_bouchetrou_meteo

! Interpoler la temperature e 2 m
      i0=0 ; i2=2
      const3=2.*grav*rhoair ! correction de pression hydrostatique e 2m
      const1=1./meteo_londlt                              !26-08-10
      const2=1./meteo_latdlt                              !26-08-10
      do j=0,jmax+1
      do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
      deci=1.+( lon_t(i,j)*rad2deg-meteo_lonstr)*const1
      i1=int(deci)
      rapi=deci-i1
      i1=i1-meteozoom_istr+1                              ! conservation la parallelisation

! Indice decimale j (latitude) dans grille aladin:
      decj=1.+( lat_t(i,j)*rad2deg-meteo_latstr)*const2
      j1=int(decj)
      rapj=decj-j1
      j1=j1-meteozoom_jstr+1                              ! conservation la parallelisation
!     i1=min(max(1,i1),meteo_imax)
!     j1=min(max(1,j1),meteo_jmax)

      teta2_t(i,j,i0)=teta2_t(i,j,i2)
      teta2_t(i,j,i2)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                     +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                     +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                     +    rapi *    rapj *meteo_var(i1+1,j1+1)

! Transformer la temperature en temperature potentielle:
        teta2_t(i,j,i2)= &
        teta2_t(i,j,i2)*(1.e5/(pss_w(i,j,i2)-const3))**.286

      enddo
      enddo

! Fermer le fichier netcdf:
      status=nf_close(ncid1)
      endif                                           !*******************>

! lire les precipiations
      call time_for_a_new_file(t_,8)!?  8 est le numero de precipi dans notebook_airseaflux
      if(decision==1) then                            !*******************>    !15-04-11
! ouvrir le fichier netcdf:
      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(par%rank==0)then
       write(6,*)
       write(6,'(a,a)')' FILE:  ',trim(texte80(1))
      endif
      if(status.ne.0)stop 'echec fichier'
! on ouvre aussi le fichier precedent pour le cas particulier ou celui-ci est necessaire pour decumuler
      status=nf_open(trim(texte80(2)),nf_nowrite,ncid2)
      if(status.ne.0)stop 'echec fichier 3'

      time1forecast_=nint(airseadt(8,1)/3600.) ! time first forecast in hours
      if(flag_meteotime==1) then !ttttttttttttttttttttttt>
! Lire le temps (en heures depuis le debut de la previ) de l'echeance:
                   status=nf_inq_varid(ncid1,'time',var_id)
      if(status/=0)stop 'echec id time'
      varstart(1)=nc1  ; varcount(1)=1          ! time
      status=nf_get_vara_double(ncid1,var_id,varstart(1:1),varcount(1:1),timearray_(1))!26-03-11
      time_=timearray_(1)!26-03-11
      if(status/=0)stop 'echec lecture time'

      else                       !ttttttttttttttttttttttt>

       time_=time1forecast_

      endif                      !ttttttttttttttttttttttt>

!> Traiter le cumul ou le non cumul:                           !25-08-10
                   status=nf_inq_varid(ncid1,'precipitation',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'RR_GDS0_SFC',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'TP',var_id)
      if(status/=0)stop 'echec id precipitation'
      texte90=''
      scalecum_=1.

      status=nf_get_att_text(ncid1,var_id,'units',texte90)
      if(status/=0)stop 'echec attribut precipitation'

      flag_cumul=0
      loop_=0
      if(texte90=='kg/m2') then !------->
! les unites nous indiquent qu'il y a eu cumul
      flag_cumul=1
      scalecum_=1./1000./airseadt(8,1) ! Pas de temps du cumul et conversion en m !30-04-11
      loop_=1
      if(nc1==1)loop_=0                   ! sauf dans le cas de la premiere echeance de chaque fichier (k0=0)
      stop 'j''ai un doute sur la maniere dont precipi a ete cumulee'
      endif                     !------->

      if(texte90=='m')     then !------->
! les unites nous indiquent qu'il y a eu cumul
      flag_cumul=1
      scalecum_=1./airseadt(8,1)       ! Pas de temps du cumul !30-04-11
      loop_=1
      if(time_==time1forecast_)loop_=0     ! sauf dans le cas de la premiere echeance de chaque fichier (k0=0)
      endif

      if(texte90=='mm')     then !------->
! les unites nous indiquent qu'il y a eu cumul
! Conversion 'mm' en 'm'
      flag_cumul=1
      scalecum_=0.001/airseadt(8,1) ! Pas de temps du cumul + conversion "mm" en "m"
      loop_=1
      if(time_==time1forecast_)loop_=0     ! sauf dans le cas de la premiere echeance de chaque fichier (k0=0)
      endif

      if(texte90=='mm h**-1')     then !------->
       flag_cumul=0
       loop_=0
       scalecum_=1.
      endif                            !------->
!>

      do loop1=loop_,0,-1

      ksecu=0
      if(flag_meteodata=='slimga') then
      ksecu=1
      varstart(4)=nc1            ; varcount(4)=1          ! time
      varstart(3)=1              ; varcount(3)=1          ! height_4
      endif
      if(flag_meteodata=='ecmwf'.or.              &
         flag_meteodata=='arome') then

      if(time_==time1forecast_.and.nc1/=1)    &
      stop 'time_==time1forecast_.and.nc1/=1'

         ksecu=1
         varstart(3)=nc1-loop1      ; varcount(3)=1          ! time
         ncid_=ncid1

! Cas particulier: Premiere echeance d'un fichier qui n'est pas la premiere echeance de previ
! il faut decumuler par rapport e la derniere echeance du jour (fichier) precedent
! varstart(3)=dernier echeance du precedent fichier soit nint(airseainfo(1,3))
         if(flag_cumul==1) then !cumulcumulcumul>
          if(time_/=time1forecast_.and.nc1==1) then !cpcpcpcpcp>
           ksecu=1
           varstart(3)=max_meteo_time_counter  ; varcount(3)=1          ! time
           ncid_=ncid2
          endif                                     !cpcpcpcpcp>
         endif                  !cumulcumulcumul>

         if(par%rank==0)then
           write(6,*)'Rain varstart loop1  ',varstart(3),loop1
           if(ncid_==ncid2)write(6,'(a,a)')'Previous file:',texte80(2)
         endif

      endif                              !11111111111>

      if(ksecu==0)stop 'meteo data file not understood'
      varstart(2)=meteozoom_jstr ; varcount(2)=meteo_jmax ! lat
      varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax ! lon
      status=nf_get_vara_real(ncid_,var_id,varstart(1:3),varcount(1:3),meteo_var)
      if(status/=0)stop 'echec lecture precipitation'

      if(scalecum_/=1.)meteo_var=meteo_var*scalecum_
      if(loop_==1.and.loop1==1)meteo_cum=meteo_var
      if(loop_==1.and.loop1==0)meteo_var=meteo_var-meteo_cum

      enddo                  ! end of loop on "loop1"

! Bouchage des trous:
      call appli_bouchetrou_meteo

      if(texte90=='mm h**-1')     then !------->
! conversion mm/h en m/s:
        x1=1./1000./3600.
        meteo_var(:,:)=meteo_var(:,:)*x1
      endif                            !------->

! Interpoler les precipitations:
      i0=0 ; i2=2
      const1=1./meteo_londlt                              !26-08-10
      const2=1./meteo_latdlt                              !26-08-10
      do j=0,jmax+1
      do i=0,imax+1

! Indice decimale i (longitude) dans grille aladin:
      deci=1.+( lon_t(i,j)*rad2deg-meteo_lonstr)*const1
      i1=int(deci)
      rapi=deci-i1
      i1=i1-meteozoom_istr+1                              ! conservation la parallelisation

! Indice decimale j (latitude) dans grille aladin:
      decj=1.+( lat_t(i,j)*rad2deg-meteo_latstr)*const2
      j1=int(decj)
      rapj=decj-j1
      j1=j1-meteozoom_jstr+1                              ! conservation la parallelisation
!     i1=min(max(1,i1),meteo_imax)
!     j1=min(max(1,j1),meteo_jmax)


      precipi_w(i,j,i0)=precipi_w(i,j,i2)

      precipi_w(i,j,i2)=(1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )   &
                       +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)   &
                       +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )   &
                       +    rapi *    rapj *meteo_var(i1+1,j1+1)


      enddo
      enddo

! Fermer le fichier netcdf:
      status=nf_close(ncid1)
      status=nf_close(ncid2)      !14-04-11
      endif                                           !*******************>

      if(par%rank==0)  &
       write(6,*)'-------------------------------------------'

      end subroutine read_interp_netcdf_airseafile

!......................................................................................

      subroutine fichier_bouchetrou_meteo
      implicit none
      integer :: loop_, loop_file_
#ifdef synopsis
       subroutinetitle='fichier_bouchetrou_meteo'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      flag_stop=0  !05-05-18
      loop_file_=0 !13-09-22


! J'ai commente les lignes suivantes pour passer de ECMWF a MESO-NH dans
! une meme liste....
! Cette routine indentifie les points de la grille meteo en terre dans la grille meteo
! et en mer sur la grille du modele oceanique. L'interpolation remplacera ces points
! par les points de la grille meteo en mer les plus proches
!     if(flag_meteodata=='ecmwf') then !ffffffffffffffff>
!       texte30=airseabinreclist(1)
!       nc=1
!     else                             !ffffffffffffffff>
!       texte30=trim(tmpdirname)//'listemeteo1'
!     endif                            !ffffffffffffffff>
!     open(unit=4,file=texte30             &
!                ,access='direct'          &
!                ,recl=540                 &
!                ,form='unformatted')
!     read(4,rec=nc)texte80(1),i
!     close(4)


! loop_file_=1 veut dire : chercher d'abord dans le fichier de la variable !13-09-22
! loop_file_=2 veut dire : le LSM n'a as ete trouvE dans le fichier de variable => on cherche dans un eventuel grid file
      if(loop_file_==0)status=nf_open(trim(texte80(1)),nf_nowrite,ncid1) 
 3740 if(loop_file_==1)status=nf_open(trim(meteo_grid_file),nf_nowrite,ncid1) 
      if(status/=0)stop 'echec fichier routine fichier_bouchetrou_meteo'

      varstart(4)=1              ; varcount(4)=1          ! time
      varstart(3)=1              ; varcount(3)=1          ! height_4
      varstart(2)=meteozoom_jstr ; varcount(2)=meteo_jmax ! lat
      varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax ! lon
                     status=nf_inq_varid(ncid1,'land_sea_mask',var_id)
      if(status.ne.0)status=nf_inq_varid(ncid1,'lsm',var_id)
      if(status.ne.0)status=nf_inq_varid(ncid1,'LSM',var_id)
      if(status.ne.0)status=nf_inq_varid(ncid1,'LSM_GDS0_SFC',var_id)
      if(status.ne.0)status=nf_inq_varid(ncid1,'masqtermer',var_id)
      if(status.ne.0) then !pmxpmx>
         if(flag_lsm==0)  then !ooo>
            if(loop_file_==0) then !essayer-un-autre-fichier> !13-09-22
              status=nf_close(ncid1)
              loop_file_=1 
              goto 3740            ! revenir au point de depart
            endif                  !essayer-un-autre-fichier> !13-09-22
            write(6,'(4a)')'Cant not find land sea mask in file ',trim(texte80(1)) &
            ,' and in file ',trim(meteo_grid_file)
            stop 'echec id land_sea_mask'
         endif                 !ooo>
         if(flag_lsm==1)  goto 3115
      endif                !pmxpmx>

      status=nf_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(status.ne.0)stop 'echec nf_inq_var land_sea_mask'

      if(var_type==nf_real) then !rrrrrr>
       status=nf_get_vara_real(ncid1,var_id               &
                                    ,varstart(1:var_dims) &
                                    ,varcount(1:var_dims) &
                    ,meteo_var(1:meteo_imax,1:meteo_jmax))
       if(status.ne.0)stop 'echec lecture land_sea_mask'
      endif                      !rrrrrr>

      if(var_type==nf_short) then !iiiii>
       status=nf_get_vara_int(ncid1,var_id               &
                                   ,varstart(1:var_dims) &
                                   ,varcount(1:var_dims) &
                 ,meteo_short(1:meteo_imax,1:meteo_jmax))
       if(status.ne.0)stop 'echec lecture land_sea_mask'
       status=nf_get_att_real(ncid1,var_id,'scale_factor',var_scalefactor)
       if(status/=0) &
       stop 'error get scale_factor 2976'
       status=nf_get_att_real(ncid1,var_id,'add_offset',var_addoffset)
       if(status/=0)stop 'error get add_offset 2976'
       do j=1,meteo_jmax ; do i=1,meteo_imax
        meteo_var(i,j)=meteo_short(i,j)*var_scalefactor+var_addoffset
       enddo ; enddo

      endif                       !iiiii>

 3115 continue

      status=nf_close(ncid1)

      if(flag_lsm==1) return

! If SURFEX LSM convention then reverse LSM : (1=land, 0=sea)
      if(meteolandconvention==0)meteo_var(:,:)=1-meteo_var(:,:) !29-04-16

#ifdef bidon
!...........................................
! Enlarge the continental signal of the LSM
      allocate(meteo_var2(meteo_imax,meteo_jmax))
      do j=1,meteo_jmax
      jm1=max(j-1,1) ; jp1=min(j+1,meteo_jmax)
       do i=1,meteo_imax
        meteo_var2(i,j)=meteo_var(i,j)+meteo_var(i,jm1)+meteo_var(i,jp1)
       enddo
      enddo
      do i=1,meteo_imax
      im1=max(i-1,1) ; ip1=min(i+1,meteo_imax)
       do j=1,meteo_jmax
        meteo_var2(i,j)=meteo_var2(i,j)+meteo_var(im1,j)+meteo_var(ip1,j)
       enddo
      enddo

      do j=1,meteo_jmax
      do i=1,meteo_imax
       meteo_var(i,j)=meteo_var2(i,j)
      enddo
      enddo
      deallocate(meteo_var2)

#endif

!...........................................
! Reperage des trous
!     const1=1./meteo_londlt                              !26-08-10
!     const2=1./meteo_latdlt                              !26-08-10

      ksecu=0 !05-06-17
      do j=0,jmax+1 ; do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then !----->

       i0=int( ij2meteo_i(i,j) - (meteozoom_istr-1) ) !21-03-15
       j0=int( ij2meteo_j(i,j) - (meteozoom_jstr-1) )

        do j1=j0,j0+1 ; do i1=i0,i0+1
          if(i1<1)         ksecu=1 !05-06-17
          if(i1>meteo_imax)ksecu=1  
          if(j1<1)         ksecu=1
          if(j1>meteo_jmax)ksecu=1
          if(ksecu==0) then !pmxpmx> !05-06-17
              if(meteo_var(i1  ,j1  )>0.)meteo_var(i1  ,j1  )=-1.   !18-09-12
          else              !pmxpmx>
           write(10+par%rank,*)'------'
           write(10+par%rank,*)'par%rank,i,j loc',par%rank,i,j
           write(10+par%rank,*)'         i,j glb',i+par%timax(1),j+par%tjmax(1)
           write(10+par%rank,*)'meteo i,j,   loc',i1,j1
           write(10+par%rank,*)'meteo ijmax  loc',meteo_imax,meteo_jmax
           write(10+par%rank,*)'meteo i,j,   glb',i1+meteozoom_istr-1,j1+meteozoom_jstr-1
           write(10+par%rank,*)'meteo ijmax  glb',meteo_imax_full,meteo_jmax_full
           goto 3474
          endif             !pmxpmx>
        enddo ; enddo

       endif                      !----->
      enddo ; enddo

 3474 k10=0
      call mpi_allreduce(ksecu,k10,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k10/=0)stop 'Err 3501 meteo grid overflow' !05-06-17
              

!...........................................
! Trouver les boucheurs:
!     open(unit=4, &
!          file=trim(tmpdirname)//'bouchetrou_meteo_'//dom_c//'.out')   !#MPI

      do 3645 loop_=0,1 !14-06-18 
      if(loop_==1) then !pmx> 
        if(allocated(meteoplugs))deallocate(meteoplugs) !26-01-22
            allocate(meteoplugs(misspointmax,4))        !14-06-18
      endif             !pmx>
      misspointmax=0    !07-05-16
! Premier passage (loop_=0) evaluation d'une nbre de trous
! Allocation du tableau A la dimension misspointmax
! Second passage pour renseigner le tableau de relation (i,j) trou vers (i,j) boucheur

      do 1862 j=1,meteo_jmax
      do 1862 i=1,meteo_imax

       if(meteo_var(i,j)==-1.)then !eeeeeeeeeeeeeee>
        ksecu=0
        i10=1
        dist1=1e10
 1864   continue

        do 1863 k1=0,1

         j0=j-i10
         j2=j+i10
         j3=k1*(j2-j0)+(1-k1)

         i0=i-i10+k1
         i2=i+i10-k1
         i3=(1-k1)*(i2-i0)+k1

         do 1861 j1=j0,j2,j3
         do 1861 i1=i0,i2,i3

         if(i1<1         ) flag_stop=1 !05-05-18
         if(i1>meteo_imax) flag_stop=1 !05-05-18
         if(j1<1         ) flag_stop=1 !05-05-18
         if(j1>meteo_jmax) flag_stop=1 !05-05-18
         if(flag_stop==1) then !pmx> !05-05-18
          write(10+par%rank,*)'Meteo grid overflow in finding a land replacement value'
          write(10+par%rank,*)'par%rank=',par%rank
          write(10+par%rank,*)'Meteo grid relative to extraction zone i,j:',i,j
          write(10+par%rank,*)'is looking for fill value in        i1,j1:',i1,j1
          write(10+par%rank,*)'Meteo grid global i ,j :                   ', i+meteozoom_istr-1,j +meteozoom_jstr-1
          write(10+par%rank,*)'Meteo grid global i1,j1:                   ',i1+meteozoom_istr-1,j1+meteozoom_jstr-1
          write(10+par%rank,*)'meteo local imin imax:',1,meteo_imax
          write(10+par%rank,*)'meteo local jmin jmax:',1,meteo_jmax
          write(10+par%rank,*)'Extraction zoom coordinates in the meteo global file:'
          write(10+par%rank,*)'Min Max i:',meteozoom_istr,meteozoom_iend
          write(10+par%rank,*)'Min Max j:',meteozoom_jstr,meteozoom_jend
          write(10+par%rank,*)
          if(i1<1)         write(10+par%rank,*)'PROBLEM: local i1<1'
          if(i1>meteo_imax)write(10+par%rank,*)'PROBLEM: local i1>meteo_imax'
          if(j1<1)         write(10+par%rank,*)'PROBLEM: local j1<1'
          if(j1>meteo_jmax)write(10+par%rank,*)'PROBLEM: local j1>meteo_jmax'
          write(10+par%rank,*)'SOLUTION: increase the value of meteo_enlarged_bounds' &
                             ,' in module_airseaflux.F90'
          write(10+par%rank,*)
          goto 3678
         endif                 !pmx>
! si flag_stop=1 c'est qu'il n'y a pas de bouchage
! suffisament proche du trou et que l'algo cherche au dela de la zone d'extraction,
! ce qu'on ne permet pas pour avoir une parallelisation parfaite.
! Ce qu'on peu faire: on peut augmenter la taille de la
! zone d'extraction en augmentant meteo_enlarged_bounds. Mais cela revele surtout
! que les masques oceano et meteo sont trop differents et qu'il faut sans doute
! ameliorer le masque oceano

          if(meteo_var(i1,j1)==0.)then !%%%%%%%%%%%%%>
           dist2=sqrt(real(i-i1)**2+real(j-j1)**2)
           ksecu=1
             if(dist2.lt.dist1) then                    !>>>>>>>>>>>>>
              i4=i1
              j4=j1
              dist1=dist2
             endif                                      !>>>>>>>>>>>>>
          endif                         !%%%%%%%%%%%%%>
 1861    continue

 1863   continue
        i10=i10+1
        if(ksecu.eq.0)goto 1864
!       write(4,'(4(i4,1x))')i,j  ,i4,j4 ; misspointmax=misspointmax+1 !07-05-16
        misspointmax=misspointmax+1 !14-06-18
        if(loop_==1) then !m°v°m>   !14-06-18
           meteoplugs(misspointmax,1)=i
           meteoplugs(misspointmax,2)=j
           meteoplugs(misspointmax,3)=i4
           meteoplugs(misspointmax,4)=j4
        endif             !m°v°m>

       endif                           !eeeeeeeeeeeeeee>
 1862 continue

 3645 continue ! boucle loop_

!     close(4)

  3678   call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0)stop &
      'Err 3680 grid overflow in module_airseaflux. See fort.xxx files' !09-05-18

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif

      end subroutine fichier_bouchetrou_meteo

!.............................................................

      subroutine appli_bouchetrou_meteo
      implicit none
#ifdef synopsis
       subroutinetitle='appli_bouchetrou_meteo'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!     if(par%rank==0)write(6,*)'bouchage des trous de la grille meteo'

!     open(unit=3, &
!          file=trim(tmpdirname)//'bouchetrou_meteo_'//dom_c//'.out')   !#MPI
      do k=1,misspointmax !07-05-16
!      read(3,*,end=1760)i,j,i4,j4
!      meteo_var(i,j)=meteo_var(i4,j4)
       meteo_var(meteoplugs(k,1),meteoplugs(k,2))=  & !14-06-18
       meteo_var(meteoplugs(k,3),meteoplugs(k,4))
      enddo
 1760 continue
!     close(3)

      return
      end subroutine appli_bouchetrou_meteo

!.............................................................

      subroutine time_for_a_new_file(t_,k_)
      implicit none
      integer k_,t_
#ifdef synopsis
       subroutinetitle='time_for_a_new_file'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      decision=0

! Si les proc ne sont pas synchro en raison du subcycling on s'abstient afin
! de garantir la continuite horizontale des champs interpoles:
      if(subcycle_synchro==0)return

      x2=(elapsedtime_now-airseadt(k_,2))/airseadt(k_,1)  ! 14-04-11
      x1=(elapsedtime_bef-airseadt(k_,2))/airseadt(k_,1)  ! 14-04-11

      if(iteration3d==kount0.or.int(x2)/=int(x1)) then !>>>>>>>>>>>>>>>>>>>>>>>>>    !19-04-11

      decision=1

      x2=x2+0.5*t_                                              ! numero d'echeance absolu et en decimal
      nc=1+int(x2)                                                 ! numero du fichier (multi-echeance)

        if(nc.le.0)then !---- debug ----->                           !06/05/04
        if(par%rank==0)write(6,*)'Erreur l 2790'
        if(par%rank==0)write(6,*)'je suis bloque dans aiseaflux_fbk choix 2'
        if(par%rank==0)write(6,*)'car nc=',nc,' est <= 0 qui signifie que je suis'
        if(par%rank==0)write(6,*)'en avance sur la premiere echeance disponible.'
        if(par%rank==0)write(6,*)'pour continuer je fixe arbitrairement nc+1 e 1'
        ncmin_airsea=1
        nc=ncmin_airsea
        !pause!11-02-20
        endif            !---- debug ----->


        write(texte30,'(a,i0)')trim(tmpdirname)//'listemeteo',k_
        open(unit=4,file=texte30            &
                   ,access='direct'         &
                   ,recl=540                &
                   ,form='unformatted')
        read(4,rec=nc         )texte80(1),nc1
        read(4,rec=max(nc-1,1))texte80(2),i
        close(4)
!       if(par%rank==0)write(6,'(a,a,a,a)')'2 fichiers consecutifs:',trim(texte80(1)),' ',trim(texte80(2))

      else                                         !>>>>>>>>>>>>>>>>>>>>>>>>>

       decision=0                                     !14-04-11

      endif                                        !>>>>>>>>>>>>>>>>>>>>>>>>>

      end subroutine time_for_a_new_file

!....................................................................

      subroutine airseaflux_bulkvar_now
      implicit none
#ifdef synopsis
       subroutinetitle='airseaflux_bulkvar_now'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! interpolation temporelle entre 2 echeances
! DEBUT:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      i0=0  ! previous fields
      i2=2  ! next fields
            ! now=1

       k1=1                                                             !22/01/08
       x2=      (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)    &   !15-04-11
          -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
       x0=1.-x2
       do j=1,jmax
       do i=1,imax
        ssr_w(i,j,1)=x0*ssr_w(i,j,i0)+x2*ssr_w(i,j,i2)
       enddo
       enddo

      if(ioffline==2)return                                          !16/06/06

       k1=2                                                             !22/01/08
       x2=      (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)    &   !15-04-11
          -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
       x0=1.-x2
       do j=1,jmax
       do i=1,imax
        snsf_w(i,j,1)=x0*snsf_w(i,j,i0)+x2*snsf_w(i,j,i2)
       enddo
       enddo

       k1=3                                                             !22/01/08
       x2=      (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)    &   !15-04-11
          -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
       x0=1.-x2
!      x0=rampe*x0 ; x2=rampe* x2                                 !11-03-10
       do j=1,jmax
       do i=1,imax
        uwind_t(i,j,1)=x0*uwind_t(i,j,i0)+x2*uwind_t(i,j,i2)             !11-03-10
        vwind_t(i,j,1)=x0*vwind_t(i,j,i0)+x2*vwind_t(i,j,i2)
       enddo
       enddo

       k1=5                                                             !22/01/08
       x2=      (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)    &   !15-04-11
          -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
       x0=1.-x2
       do j=1,jmax
       do i=1,imax
        pss_w(i,j,1)=x0*pss_w(i,j,i0)+x2*pss_w(i,j,i2)
       enddo
       enddo
       pss_mean(1)=x0*pss_mean(i0)+ x2*pss_mean(i2)               !01/11/04

       k1=6                                                             !22/01/08
       x2=      (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)    &   !15-04-11
          -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
       x0=1.-x2
       do j=1,jmax
       do i=1,imax
        teta2_t(i,j,1)=x0*teta2_t(i,j,i0)+x2*teta2_t(i,j,i2)
       enddo
       enddo

       k1=7                                                             !22/01/08
       x2=      (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)    &   !15-04-11
          -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
       x0=1.-x2
       do j=1,jmax
       do i=1,imax
        q2_t(i,j,1)=x0*q2_t(i,j,i0)+x2*q2_t(i,j,i2)
       enddo
       enddo

       k1=8                                                             !22/01/08
       x2=      (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1)    &   !15-04-11
          -int( (elapsedtime_now-airseadt(k1,2))/airseadt(k1,1) )
       x0=1.-x2
       do j=1,jmax
       do i=1,imax
        precipi_w(i,j,1)=x0*precipi_w(i,j,i0)+x2*precipi_w(i,j,i2)
       enddo
       enddo

      if(irelaxsst==1)call airseaflux_relaxed_sst

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! interpolation temporelle entre 2 echeances
! FIN.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      end subroutine airseaflux_bulkvar_now

!...........................................................................
      subroutine airseaflux_glorys(case_)
      implicit none
      integer case_
#ifdef synopsis
       subroutinetitle='airseaflux_glorys'
       subroutinedescription='Principal driver for meteo "GLORYS" case'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(.not.allocated(varcode))allocate(varcode(nairsea)) !08-03-13

      if(case_==1) then !- Initial phase -> !15-03-13
       call airseaflux_initial
       call airseaflux_glorys_interp(0,case_)   ! read netcdf & grid interpolation
       call airseaflux_glorys_interp(2,case_)   ! read netcdf & grid interpolation
      endif            !- Initial phase ->

      if(case_==2)call airseaflux_glorys_interp(2,case_)   ! read netcdf & grid interpolation

      call airseaflux_flux_moveforward ! time linear interpolation

      end subroutine airseaflux_glorys

!...........................................................................

      subroutine airseaflux_model_flag
      implicit none
#ifdef synopsis
       subroutinetitle='airseaflux_model_flag'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Idendifier le type de fichier netcdf e l'aide du nom de la liste des fichiers meteo
        ksecu=0
        k1=1     ! 25-01-11

        if(index(airseafile(k1),'slimga')/=0) then
        ksecu=1
        flag_meteodata='slimga'
        endif

        if(index(airseafile(k1),'arome')/=0) then
        ksecu=1
        flag_meteodata='arome'
        endif

        if(index(airseafile(k1),'glorys')/=0) then !07-02-14
        ksecu=1
        flag_meteodata='glorys'
        endif


        if(ksecu==0) then
        if(par%rank==0)write(6,*)k1,airseafile(k1),'airseafile'
        stop 'The meteo data file is not understood 1'
        endif

      end subroutine airseaflux_model_flag

!...........................................................................

      subroutine airseaflux_binary_file_list
      implicit none
      double precision time_,timeforecast_,time1_,time2_
!     double precision valr8_(8)
      integer ncid_,flag_                                 &
             ,year_,month_,day_,hour_,minute_,second_     &
             ,loop1_,linemax_,loop2_,loop3_,vartype_,loop_var_
      integer , dimension(:) , allocatable :: tabi4_,nelt
      double precision , dimension(:) , allocatable :: tabr8_
#ifdef synopsis
       subroutinetitle='airseaflux_binary_file_list'
       subroutinedescription= &
          'Reads the text lists in notebook_airseaflux and creates' &
       //' the binary record direct access lists'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      cpu_seconds=MPI_Wtime ( )
      if(     flag_meteodata/='ecmwf'              &
         .and.flag_meteodata/='arome'              & !08-11-16
         .and.flag_meteodata/='glorys')            &
      stop 'Do not know how to make the lists for this flag_meteodata'

! A l'etat initial, faire la liste binaire des fichiers meteo
! Reset zero (important car mpi_sum plus tard)
      airseafile_nextrec(:)=0
      airseafile_nextime(:)=0.

!......... !26-01-22
! name of the record direct access binary file list:
      do loop_var_=1,nairsea !----> !26-01-22
        k2=index(airseafile(loop_var_),'/',back=.true.) ; k3=index(airseafile(loop_var_),' ')
        airseabinreclist(loop_var_)=trim(tmpdirname)//''//airseafile(loop_var_)(k2+1:k3-1)//'.binrec'
      enddo ! loop_var_      !----> !26-01-22

! Verifier que la liste n'existe pas dEjA:
       if(par%rank==0) then !-rank0 checks existing list->
        open(unit=3,file=airseabinreclist(1) &
            ,access='direct'                      &
            ,recl=540                             &
            ,form='unformatted'                   &
            ,status='new'                         &
            ,iostat=k0                            &
            )
        close(3)
       endif                !-rank0 checks existing list->
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)  !23-06-13
#endif
! Envoyer k0 A tous les autres proc:
      call mpi_bcast(k0,1,mpi_integer,0,par%comm2d,ierr) !19-04-16

! Si iostat=k0 est nul cela signifie que la liste binaire n'existe
! pas et qu'il faut la creer. Dans le cas contraire cela est interprEtE comme le
! fait qu l'utilisateur a intentionnellement placE des listes binaires (d'un precedent run)
! dans le repertoire tmp dans le but de ne pas recalculer les listes binaires. Le modEle
! prend alors l'initiative de considerer que les listes existent dEjA et prend alors 
! l'aiguillage 2453 qui lui permet de se positionner dans la liste binaire existante.
! Sinon il continue comme d'hab.
       if(k0/=0) then !existing old list in tmp directory >
        do loop_var_=1,nairsea         
         call airseaflux_oldlist(loop_var_)
        enddo
        goto 4505
       endif          !existing old list in tmp directory >

! flag_cumul et scalefct: !26-01-22
      allocate(flag_cumul_all(nairsea)) ; flag_cumul_all=0
      allocate(scalefct_all(nairsea))   ; scalefct_all=1.
      do loop_var_=1,nairsea !----> !26-01-22
       open(unit=3,file=airseafile(loop_var_))
            read(3,'(a)')texte80(1)
       close(3)
       call airseaflux_varname(loop_var_) ! donne meteovarname le nom de la variable du fichier netcdf
       status=nf_open(trim(texte80(1)),nf_share,ncid_) 
          call airseaflux_inquire_var(loop_var_,ncid_) ! returns flag_cumul scalefct
       status=nf_close(ncid_)
       flag_cumul_all(loop_var_)=flag_cumul !26-01-22
         scalefct_all(loop_var_)=scalefct   !26-01-22
      enddo ! loop_var_      !----> !26-01-22
!......... !26-01-22

! Boucle sur les variables meteo: une variable = un fichier "liste" a traiter
!     do loop1_=1,nairsea
      do loop1_=1,1 !Ne pas perdre du temps puisque fichier unique !26-01-22

! Get meteovarname, the variable to be found in the meteo netcdf file
      call airseaflux_varname(loop1_) ! donne meteovarname le nom de la variable du fichier netcdf

         open(unit=3,file=airseafile(loop1_))
! Combien y a t'il de lignes dans le fichier?
         linemax_=0
 1834    read(3,'(a)',end=1832)texte250 ; linemax_=linemax_+1 ; goto 1834
 1832    rewind(3)
         if(linemax_==0) then         !>>>>>
          write(6,'(a,a,a)')'File ',trim(airseafile(loop1_)),' is empty'
          stop ' stop airseaflux_binary_file_list 3284'
         endif                        !>>>>>

! Bien que la liste binaire ne soit faite qu'A partir d'un seul champs, on veut tout de meme connaitre le nom du fichier contenant les autres
! champs au cas oU chaque champ aurait son propre fichier. Pour cela on lit (on ne faire que lire le nom) les listes des autres champs:
         do loop_var_=2,nairsea ; open(unit=10+loop_var_,file=airseafile(loop_var_)) ; enddo !14-09-22

! Chaque proc va traiter une fraction du fichier. On commence par derouler les lignes !10-08-14
! qui ne concernent pas le proc
         do loop3_=1,int(real(par%rank  )/real(nbdom)*linemax_)
          read(3,*)                                             ! lire ces lignes pour rien
          do loop_var_=2,nairsea ; read(10+loop_var_,*) ; enddo ! lire ces lignes pour rien !14-09-22
         enddo

         write(texte60,'(a,i0)')trim(tmpdirname)//'tmpfile',par%rank
         open(unit=4,file=texte60,status='REPLACE')
          do loop3_=int(real(par%rank  )/real(nbdom)*linemax_)+1    &
                   ,int(real(par%rank+1)/real(nbdom)*linemax_)

! meteo file name given by the list:
          read(3,'(a)',end=332)texte80(1)
          do loop_var_=2,nairsea ; read(10+loop_var_,'(a)')texte80(loop_var_) ; enddo !14-09-22
! open meteo file:
!         status=nf_open(trim(texte80(1)),nf_nowrite,ncid_)
          status=nf_open(trim(texte80(1)),nf_share,ncid_) !23-03-15 sur conseil de Cyril
          if(status/=0) &
          stop 'error nf_open in subroutine airseaflux_binary_file_list'

! Verifier la coherence de la liste des fichiers meteo ! 25-04-16
!         call airseaflux_checkfileconsistecy(ncid_)

! Get the scalefactor and the flag_cumul of the variable:
!         call airseaflux_inquire_var(loop1_,ncid_) ! returns flag_cumul scalefct

! Get max_meteo_time_counter:
          max_meteo_time_counter=1 ; flag_meteotime=1
                       status=nf_inq_dimid(ncid_,'time_counter',dim_t_id)
          if(status/=0)status=nf_inq_dimid(ncid_,'TIME_COUNTER',dim_t_id)
          if(status/=0)status=nf_inq_dimid(ncid_,'time',dim_t_id)
          if(status/=0)status=nf_inq_dimid(ncid_,'TIME',dim_t_id)
          if(status/=0)status=nf_inq_dimid(ncid_,'Time',dim_t_id) !17-09-21
          if(status/=0)flag_meteotime=0
          if(status==0) then !>>>
                       status=nf_inq_dimlen(ncid_,dim_t_id,max_meteo_time_counter)
                    if(status/=0)stop 'Err 4267 nf_inq_dimlen'
          endif              !>>>
! Get time units:
                      status=nf_inq_varid(ncid_,'time',var_id)
         if(status/=0)status=nf_inq_varid(ncid_,'time_counter',var_id)
         if(status/=0)stop 'erreur nf_inq_varid time airsea'
         texte60=''
         status=nf_get_att_text(ncid_,var_id,'units',texte60)
         if(status/=0)stop 'erreur nf_get_att_text airsea'
         k=index(texte60,' since ')
         read(texte60(k+7 :k+10),*)year_

! Avant v250:
!        read(texte60(k+12:k+13),*)month_
!        read(texte60(k+15:k+16),*)day_
!        read(texte60(k+18:k+19),*)hour_
         k2=index(texte60,':')
         read(texte60(k2-2:k2-1),*)hour_
!        read(texte60(k+21:k+22),*)minute_ ! note: les secondes ne sont pas toujours presentes dans fichiers ecmwf

! Apres v250:
         k1=index(texte60,'-')             !31-03-19
         k2=index(texte60,'-',back=.true.) !31-03-19
         read(texte60(k1+1:k2-1),*)month_  !31-03-19
         read(texte60(k2+1:k2+2),*)day_    !31-03-19
         k2=index(texte60,':')             !31-03-19
         read(texte60(k2-2:k2-1),*)hour_   !31-03-19
         read(texte60(k2+1:k2+2),*)minute_ !31-03-19! note: les secondes ne sont pas toujours presentes dans fichiers ecmwf


! Elapsed time (seconds) since the reference date of the units:
         call datetokount(year_,month_,day_,hour_,minute_,0) ! returns elapsedtime_out

! Get time of each field contained in the opened meteo netcdf file:
!        status=nf_get_var_double(ncid_,var_id,valr8_(1:8))
!        status=nf_get_var_int(ncid_,var_id,valr8_(1:8))
!        if(status/=0)stop 'erreur nf_get_vara_double airsea'

         status=nf_inq_vartype(ncid_,var_id,vartype_) !14-02-16
         if(status/=0)stop 'Err 3707 nf_inq_vartype'
         do loop2_=1,max_meteo_time_counter

          status=-9999
          if(vartype_==nf_double) then
            status=nf_get_vara_double(ncid_,var_id,loop2_,1,x3)
            if(status/=0)stop 'err 3716A nf_get x3'
          endif
          if(vartype_==nf_real)then !14-02-16
           status=nf_get_vara_real(ncid_,var_id,loop2_,1,x3_r4)
           if(status/=0)stop 'err 3716B nf_get x3_r4'
           x3=x3_r4
          endif
          if(vartype_==nf_int)then !14-02-16
           status=nf_get_vara_int(ncid_,var_id,loop2_,1,k0)
           if(status/=0)stop 'err 3716B nf_get k0'
           x3=real(k0)
          endif
          if(status==-9999) then !>>>>>
            write(6,'(a,a)')'File=',trim(texte80(1))
            write(6,*)'vartype_  ',vartype_
            write(6,*)'nf_double ',nf_double
            write(6,*)'nf_real   ',nf_real
            write(6,*)'nf_int    ',nf_int
            stop 'err -9999 nf_get x3'
          endif                  !>>>>>

! Rustine pour contourner le bug netcdf qui se traduit par x3=0 (de maniere aleatoire)
!         x3=3.*loop2_
          x2=-999.
          if(index(texte60,'days')/=0)x2=86400.
          if(index(texte60,'hours')/=0)x2=3600.
          if(index(texte60,'seconds')/=0)x2=1.
          if(x2==-999.)stop 'airsea unites de temps incomprises'
          timeforecast_=x3*x2
!         timeforecast_=valr8_(1)*x2
!         timeforecast_=valr8_(loop2_)*x2
          time_=elapsedtime_out+timeforecast_
! write in the ascii temporary file:
!         write(4,'(a)')trim(texte80(1))  ! Nom fichier netcdf
          do loop_var_=1,nairsea !14-09-22
          write(4,'(a)')trim(texte80(loop_var_))  ! Nom fichier netcdf
          enddo
          write(4,*)loop2_                ! Numero d'echeance dans le fichier netdf
          write(4,*)time_                 ! temps en seconde du champ
          write(4,*)timeforecast_         ! temps en seconde depuis le debut du run ecmwf
          do loop_var_=1,nairsea              !26-01-22
          write(4,*)flag_cumul_all(loop_var_) ! Cumul ou pas cumul ?
          write(4,*)scalefct_all(loop_var_)   ! facteur d'echelle aditionnel
          enddo

         enddo ! loop2_

! Close meteo netcdf file:
          status=nf_close(ncid_)

          enddo ! loop3_
 332     close(4)
         close(3)
         do loop_var_=2,nairsea ; close(10+loop_var_) ; enddo !14-09-22

! C'est le proc zero qui a la mission d'assembler tous les fichiers tmp en un seul fichier liste binaire
! La barriere suivante permet de s'assurer que tous les fichiers individuels ont bien ete faits au moment
! d'attaquer la concatenation
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)
#endif
       if(par%rank==0) then !00000000000>
       do loop_var_=1,nairsea

        time2_=-999.
        open(unit=3,file=airseabinreclist(loop_var_) &
            ,access='direct'                      &
            ,recl=540                             &
            ,form='unformatted')

        nc=1
        do loop2_=0,nbdom-1 
         write(texte60,'(a,i0)')trim(tmpdirname)//'tmpfile',loop2_
         open(unit=4,file=texte60)

 3366     time1_=time2_
!         texte80(2)=texte80(1)             ! operation utile pour detection de listes desordonnees
          texte90=texte80(1)                ! operation utile pour detection de listes desordonnees
          do k=1,nairsea                    !14-09-22
          read(4,'(a)',end=3367)texte80(k)  ! Nom fichier netcdf
          enddo
          read(4,*)k10                      ! Numero d'echeance dans le fichier netdf
          read(4,*)time2_                   ! temps en seconde du champ
          read(4,*)timeforecast_            ! temps en seconde depuis le debut du run ecmwf
          do k=1,nairsea              !26-01-22
          read(4,*)flag_cumul_all(k)  ! Cumul ou pas cumul ?
          read(4,*)scalefct_all(k)    ! facteur d'echelle aditionnel
          enddo
          flag_cumul=flag_cumul_all(loop_var_)
          scalefct=scalefct_all(loop_var_)
          if(time1_==-999.)time1_=time2_
          if(flag_cumul==1.or.flag_meteo_average==1) then !>>>
           time_=0.5*(time1_+time2_)
          else                                            !>>>
           time_=time2_
          endif                                           !>>>
! Detection de listes desordonnees: !27-02-15
               if(time2_<time1_) then   !dbdbdbdb>
                write(6,'(a,a,a)')'File ',trim(texte80(1)) &
                ,' does not respect the chronological order'
                write(6,*)'Corresponding time2_',time2_
                write(6,'(a,a)')'Previous file ',trim(texte90) ! trim(texte80(2))
                write(6,*)'Corresponding time1_',time1_
                write(6,*)'loop2_=',loop2_
                stop ' Stop module_airseaflux'
               endif                    !dbdbdbdb>

          if(time_<=elapsedtime_now) then !>>>> ce test garantit entre autres que model3d ne relancera pas une
           airseafile_nextrec(loop_var_)=nc  !     procedure de remise a jour des champs a la premiere iteration
           airseafile_nextime(loop_var_)=time_
          endif                           !>>>>

          if(nc/=1) then !pmxpmx>
           write(3,rec=nc)texte80(loop_var_) & ! Nom fichier netcdf !14-09-22
                         ,k10                & ! Numero d'echeance dans le fichier netdf
                         ,time_              & ! temps en seconde dans le repere de S26
                         ,timeforecast_      & ! temps en seconde depuis le debut du run ecmwf
                         ,flag_cumul         & ! Cumul ou pas cumul ?
                         ,scalefct           & ! facteur d'echelle aditionnel
                         ,flag_refresh_interp  

          else           !pmxpmx>
           write(3,rec=nc,iostat=i10)        &
                          texte80(loop_var_) & ! Nom fichier netcdf !14-09-22
                         ,k10                & ! Numero d'echeance dans le fichier netdf
                         ,time_              & ! temps en seconde dans le repere de S26
                         ,timeforecast_      & ! temps en seconde depuis le debut du run ecmwf
                         ,flag_cumul         & ! Cumul ou pas cumul ?
                         ,scalefct           & ! facteur d'echelle aditionnel
                         ,datesim(1:6,1)       ! date repere nc=1
           if(i10/=0)stop 'Err 3889'
          endif          !pmxpmx>
          nc_meteo_max=nc !la derniere valeur de nc sera le max !26-01-22

! Une liste ascii pour virification:
!         if(loop1_==1)write(66,'(a)')trim(texte80(1))
!         if(loop1_==1)write(66,*)    k10
!         if(loop1_==1)write(66,*)time_,time_/86400.

          nc=nc+1


          goto 3366
 3367     close(4)

        enddo ! loop2_=0,nbdom-1 

! Finalement on utilise rec=1 pour archiver date "repere" pour elapsedtime_now=0
!       write(unit=3,rec=1)datesim(1:6,1)

        close(3)
        if(airseafile_nextrec(loop_var_)==0) then !>>>>
         write(6,*)'No meteo file for the period of simulation'
         stop 'Err 3940 airseaflux_binary_file_list'
        endif                                  !>>>>

       enddo ! loop_var_=1,nairsea
       endif                !00000000000>
       call mpi_bcast(nc_meteo_max,1,mpi_integer,0,par%comm2d,ierr)
 2453 continue
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)
#endif

      enddo ! loop1_ sur variables

      if(meteo_grid_file=='none') & ! Ce check n'est pertinent que si LSM est dans le fichier de variables !13-09-22
      call airseaflux_check_lsm     ! DeplacE le 01-07-22 avant "4505 continue" car on ne peut 
                                    ! pas verifier le masque si on lit les listes binaires du repertoire binary
 4505 continue !26-01-22
! Share with all other ranks:
      allocate(nelt(1))
      nelt = ubound(airseafile_nextime)    &
            -lbound(airseafile_nextime)+1
      call mpi_bcast(airseafile_nextime    & !19-04-16
             ,nelt(1)                      &
             ,mpi_double_precision         &
             ,0,par%comm2d,ierr)

      nelt = ubound(airseafile_nextrec)    &
            -lbound(airseafile_nextrec)+1
      call mpi_bcast(airseafile_nextrec    &
             ,nelt(1)                      &
             ,mpi_integer                  &
             ,0,par%comm2d,ierr)
      deallocate(nelt)
      if(allocated(flag_cumul_all))deallocate(flag_cumul_all) !30-06-22
      if(allocated(scalefct_all))  deallocate(scalefct_all)   !30-06-22

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)  !13-10-10
#endif
      cpu_seconds=MPI_Wtime ( ) - cpu_seconds
      if(par%rank==0) &
      write(6,*)'airseaflux_binary_file_listcpu_seconds=',cpu_seconds
      end subroutine airseaflux_binary_file_list

!...........................................................................

      subroutine airseaflux_extract_zone(txtcase_)
      implicit none
      character(len=*)txtcase_
      integer(kind=1)loop_
#ifdef synopsis
       subroutinetitle='airseaflux_extract_zone'
       subroutinedescription=                                      &
          'Determines the bounds of the reading zone in the meteo' &
       //' netcdf files'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(txtcase_=='initial') then !-- Initial Case -->

      if(    flag_meteodata=='ecmwf'                           &
         .or.flag_meteodata=='arome'                           &
         .or.flag_meteodata=='glorys') then !ffffffffffffffff>

        texte30=airseabinreclist(1)
!       nc=1
        nc=airseafile_nextrec(1) !26-01-22

      else                                  !ffffffffffffffff>

      x2=(elapsedtime_now-airseadt(1,2))/airseadt(1,1)
      nc=1+int(x2)

        if(nc.le.0)then !---- debug ----->                           !06/05/04
        if(par%rank==0)write(6,*)'Erreur l 3150'
        if(par%rank==0)write(6,*)'je suis bloque dans aiseaflux_fbk choix 2'
        if(par%rank==0)write(6,*)'car nc=',nc,' est <= 0 qui signifie que je suis'
        if(par%rank==0)write(6,*)'en avance sur la premiere echeance disponible.'
        if(par%rank==0)write(6,*)'pour continuer je fixe arbitrairement nc+1 e 1'
        ncmin_airsea=1
        nc=ncmin_airsea
        !pause!11-02-20
        endif            !---- debug ----->

        k1=1 ; write(texte30,'(a,i0)')trim(tmpdirname)//'listemeteo',k1

      endif                                 !ffffffffffffffff>

        open(unit=4,file=texte30                                        &
                   ,access='direct'          &
                   ,recl=540                 &
                   ,form='unformatted')
        read(4,rec=nc)texte80(1),i
        close(4)

      endif                        !-- Initial Case -->

      if(txtcase_=='current') then !-- Iterative Case -->
       if(par%rank==0)write(6,'(a,a)')                   &
        'Refreshing interpolation procedure' &
       ,' of meteo files regarding meteo and S grids relations'
      endif                        !-- Iterative Case -->

        if(par%rank==0)write(6,'(a,a)')                    &
                'Fichier lu par airseaflux_extract_zone: ' &
               ,trim(texte80(1))

      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(status.ne.0)stop 'echec fichier 1'

                   status=nf_inq_dimid(ncid1,'lon',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'longitude',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'g0_lon_1',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'x',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'west_east',dim_x_id) !17-09-21
      if(status/=0)stop 'Erreur dim_x_id variable meteo longitude'

                   status=nf_inq_dimid(ncid1,'g0_lat_0',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'latitude',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'lat',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'y',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'south_north',dim_y_id) !17-09-21
      if(status/=0)stop 'Erreur dim_y_id variable meteo latitude'


      status=nf_inq_dimlen(ncid1,dim_x_id,meteo_imax)
      status=nf_inq_dimlen(ncid1,dim_y_id,meteo_jmax)

! garder en memoire pour verifier la coherence des fichiers meteo mis dans la
! liste des fichiers de forcage meteo:
      meteo_imax_full=meteo_imax !25-04-16
      meteo_jmax_full=meteo_jmax !25-04-16

      if(txtcase_=='current') then !-- Iterative Case -->
       if(par%rank==0)write(6,'(a,a,a,2i6)')                   &
        'Refreshing interpolation procedure' &
       ,' of meteo files regarding meteo and S grids relations' &
       ,' meteo_imax_full,meteo_jmax_full=',meteo_imax_full,meteo_jmax_full
      endif                        !-- Iterative Case -->

! on part sur l'hypothese que l'on connait deje meteo_kmax
      meteo_kmax=2

! une premiere allocation pour pouvoir charger les tableaux lon lat
      call allocate_forcages(1,2,meteo_imax,meteo_jmax,meteo_kmax) ! arg1=allouer arg2=meteo arg3,4,5=dimensions

! lire les longitudes:
                   status=nf_inq_varid(ncid1,'lon',lon_id)
      if(status/=0)status=nf_inq_varid(ncid1,'longitude',lon_id)
      if(status/=0)status=nf_inq_varid(ncid1,'g0_lon_1',lon_id)
      if(status/=0)status=nf_inq_varid(ncid1,'nav_lon',lon_id)
      if(status/=0)stop 'Erreur lon_id variable meteo longitude'

      status=nf_inq_var(ncid1,lon_id,texte30,var_type,var_dims,tabdim,i4)
      if(status/=0)stop 'Erreur nf_inq_var meteo'

      if(var_type/=nf_float.and.var_type/=nf_double) &
      stop 'Erreur type sur longitude meteo'

      if(var_dims==1) then
       i1=meteo_imax ; j1=1 ; varstart(1)=1 ; varcount(1)=meteo_imax
      endif
      if(var_dims==2) then
       i1=meteo_imax ; j1=meteo_jmax
       varstart(1)=1 ; varcount(1)=meteo_imax ; varstart(2)=1 ; varcount(2)=meteo_jmax
      endif

      if(var_type==nf_double) then  !dddddddd>
       status=nf_get_vara_double(ncid1,lon_id,varstart(1:var_dims)   &
                                              ,varcount(1:var_dims)   &
                                             ,meteo_lon(1:i1,1:j1))
      endif                      !dddddddd>
      if(var_type==nf_real) then    !rrrrrrrr>
       allocate(meteo_lonlat_r4(meteo_imax,meteo_jmax))
       status=nf_get_vara_real(ncid1,lon_id,varstart(1:var_dims)      &
                                            ,varcount(1:var_dims)      &
                                           ,meteo_lonlat_r4(1:i1,1:j1))
       meteo_lon(1:i1,1:j1)=meteo_lonlat_r4(1:i1,1:j1)
       deallocate(meteo_lonlat_r4)
      endif                      !rrrrrrrr>
      if(status/=0)stop 'Erreur lecture longitude meteo'

! lire les latitudes:
                   status=nf_inq_varid(ncid1,'lat',lat_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude',lat_id)
      if(status/=0)status=nf_inq_varid(ncid1,'g0_lat_0',lat_id)
      if(status/=0)status=nf_inq_varid(ncid1,'nav_lat',lat_id)
      if(status/=0)stop 'Erreur lat_id variable meteo latitude'

      status=nf_inq_var(ncid1,lat_id,texte30,var_type,var_dims,tabdim,i4)

      if(var_type/=nf_float.and.var_type/=nf_double) &
      stop 'Erreur type sur latitude meteo'

      if(var_dims==1) then
       i1=1 ; j1=meteo_jmax ; varstart(1)=1 ; varcount(1)=meteo_jmax
      endif
      if(var_dims==2) then
       i1=meteo_imax ; j1=meteo_jmax
       varstart(1)=1 ; varcount(1)=meteo_imax ; varstart(2)=1 ; varcount(2)=meteo_jmax
      endif

      if(var_type==nf_double) then  !dddddddd>
       status=nf_get_vara_double(ncid1,lat_id,varstart(1:var_dims)   &
                                              ,varcount(1:var_dims)   &
                                             ,meteo_lat(1:i1,1:j1))
      endif                      !dddddddd>
      if(var_type==nf_real) then    !rrrrrrrr>
       allocate(meteo_lonlat_r4(meteo_imax,meteo_jmax))
       status=nf_get_vara_real(ncid1,lat_id,varstart(1:var_dims)        &
                                            ,varcount(1:var_dims)        &
                                           ,meteo_lonlat_r4(1:i1,1:j1))
       meteo_lat(1:i1,1:j1)=meteo_lonlat_r4(1:i1,1:j1)
       deallocate(meteo_lonlat_r4)
      endif                      !rrrrrrrr>
      if(status/=0)stop 'Erreur lecture latitude meteo'


      if(var_dims==1) then !1111111111111111111111>
       do j=1,meteo_jmax
       do i=1,meteo_imax
        meteo_lon(i,j)=meteo_lon(i,1)
        meteo_lat(i,j)=meteo_lat(1,j)
       enddo
       enddo
       meteo_lonstr=meteo_lon(1,1)
       meteo_latstr=meteo_lat(1,1)
       meteo_londlt=(meteo_lon(meteo_imax,1)-meteo_lon(1,1))/real(meteo_imax-1)
       meteo_latdlt=(meteo_lat(1,meteo_jmax)-meteo_lat(1,1))/real(meteo_jmax-1)

!      write(10+par%rank,*)meteo_lon(meteo_imax,1),meteo_lon(1,1),meteo_imax
!      write(10+par%rank,*)meteo_londlt
!      write(10+par%rank,*)meteo_lat(1,meteo_jmax),meteo_lat(1,1),meteo_jmax
!      write(10+par%rank,*)meteo_latdlt

      endif                !1111111111111111111111>

! A partir des min & max des (lon,lat) des 2 grilles, reduire la zone d'extraction:
       meteozoom_istr=999999 ; meteozoom_iend=-999999 ! first guess
       meteozoom_jstr=999999 ; meteozoom_jend=-999999 ! first guess

! Si grid overflow non autorisE alors verifier les limites des domaines
      if(flag_gridoverflow==0) then !pmwpmw>
       ksecu=0 !05-06-17
       if(minval(meteo_lon)>lonmin.or.     &
          minval(meteo_lat)>latmin.or.     &
          maxval(meteo_lon)<lonmax.or.     &
          maxval(meteo_lat)<latmax) then !>>>>>>
        write(6,*)'------------------------------------'   !08-01-15
        write(6,*)'warning: the ocean grid is not 100% ' & !08-01-15
                 ,'included in the meteo field grid.'    &
                 ,' Details in fort.xxx files'     

        write(10+par%rank,*)'------------------------------------'   !20-01-17
        write(10+par%rank,*)'WARNING the ocean grid is not 100% ' &
                           ,'included in the meteo field grid.'
        if(minval(meteo_lon)>lonmin) &
        write(10+par%rank,*)'meteolonmin > lonmin:',minval(meteo_lon),lonmin

        if(maxval(meteo_lon)<lonmax) &
        write(10+par%rank,*)'meteolonmax < lonmax:',maxval(meteo_lon),lonmax

        if(minval(meteo_lat)>latmin) &
        write(10+par%rank,*)'meteolatmin > latmin:',minval(meteo_lat),latmin

        if(maxval(meteo_lat)<latmax) &
        write(10+par%rank,*)'meteolatmax < latmax:',maxval(meteo_lat),latmax

        write(10+par%rank,*)'------------------------------------'   !08-01-15
        write(10+par%rank,*)'Ocean model grid par%rank:',par%rank
        write(10+par%rank,*)'lonmin lonmax ',lonmin,lonmax
        write(10+par%rank,*)'latmin latmax ',latmin,latmax
        write(10+par%rank,*)
        write(10+par%rank,*)'Meteo model:'
        write(10+par%rank,*)'lon min max ',minval(meteo_lon),maxval(meteo_lon)
        write(10+par%rank,*)'lat min max ',minval(meteo_lat),maxval(meteo_lat)
        write(10+par%rank,*)'------------------------------------'   !08-01-15

        ksecu=1
       endif
       k10=0
       call mpi_allreduce(ksecu,k10,1,mpi_integer,mpi_sum,par%comm2d,ierr)
       if(k10/=0)stop &
      'Consider accepting meteogrid overflow using flag_gridoverflow=1' !05-06-17

      endif                         !pmwpmw>

       ksecu=0
       do j=1,meteo_jmax
       do i=1,meteo_imax

       if(meteo_lon(i,j)>=lonmin.and.meteo_lon(i,j)<=lonmax.and.   &
          meteo_lat(i,j)>=latmin.and.meteo_lat(i,j)<=latmax)then

         meteozoom_istr=min(meteozoom_istr,i)
         meteozoom_jstr=min(meteozoom_jstr,j)
         meteozoom_iend=max(meteozoom_iend,i)
         meteozoom_jend=max(meteozoom_jend,j)
         ksecu=1

       endif

       enddo
       enddo

! Si ksecu=0 c'est que le proc est si petit, qu'aucun point meteo ne s'est trouve e l'interieur
! ou alors que le modele meteo est plus petit que la grille symphonie
! On refait le test differement. On cherche les points les plus proches des limites min max du domaine
       if(ksecu==0)then !000000000000000>
        do loop_=1,4 !10-07-18
         dist1=1.d20
         if(loop_==1)then ; x1=latmin*deg2rad ; x2=lonmin*deg2rad ; endif
         if(loop_==2)then ; x1=latmax*deg2rad ; x2=lonmin*deg2rad ; endif
         if(loop_==3)then ; x1=latmax*deg2rad ; x2=lonmax*deg2rad ; endif
         if(loop_==4)then ; x1=latmin*deg2rad ; x2=lonmax*deg2rad ; endif
         do j=1,meteo_jmax
         do i=1,meteo_imax

!         dist2=rayonterre*                                       &
!         acos( sin(meteo_lat(i,j)*deg2rad)*sin(lat_t(i1,j1))     &
!              +cos(meteo_lat(i,j)*deg2rad)*cos(lat_t(i1,j1))     &
!              *cos(lon_t(i1,j1)-meteo_lon(i,j)*deg2rad))
          dist2=rayonterre*acos( min(max(sin(meteo_lat(i,j)*deg2rad)*sin(x1)+cos(meteo_lat(i,j)*deg2rad)*cos(x1)*cos(x2-meteo_lon(i,j)*deg2rad),-1.),1.) )

          if(dist2<dist1)then !--->
           dist1=dist2 ; i2=i ; j2=j
          endif               !--->

         enddo
         enddo
         meteozoom_istr=min(i2,meteozoom_istr) ; meteozoom_iend=max(i2,meteozoom_iend)
         meteozoom_jstr=min(j2,meteozoom_jstr) ; meteozoom_jend=max(j2,meteozoom_jend)
       enddo  ! loop_
       endif            !000000000000000>



! on elargit un peu plus pour boucher les trous sans etre restreint par la taille
! reduite de la zone d'extraction, et puis aussi pour rattraper l'erreur liee au fait
! que plusieurs grilles meteo (point u, v, t) peuvent etre presentes.
       meteozoom_istr=max(meteozoom_istr-meteo_enlarged_bounds,1) !05-05-18
       meteozoom_iend=min(meteozoom_iend+meteo_enlarged_bounds,meteo_imax)
       meteozoom_jstr=max(meteozoom_jstr-meteo_enlarged_bounds,1)
       meteozoom_jend=min(meteozoom_jend+meteo_enlarged_bounds,meteo_jmax)

       meteo_imax=meteozoom_iend-meteozoom_istr+1
       meteo_jmax=meteozoom_jend-meteozoom_jstr+1

!      write(*,*)'meteozoom_istr=',meteozoom_istr
!      write(*,*)'meteozoom_jstr=',meteozoom_jstr
!      write(*,*)'meteozoom_iend=',meteozoom_iend
!      write(*,*)'meteozoom_jend=',meteozoom_jend
!      write(*,*)'meteo_imax=',meteo_imax
!      write(*,*)'meteo_jmax=',meteo_jmax

! Refaire l'allocation dynamique en fonction des dimensions reduites:
      call allocate_forcages(2,2,0,0,0) ! arg1=desallouer arg2=meteo
      call allocate_forcages(1,2,meteo_imax,meteo_jmax,meteo_kmax) ! arg1=allouer arg2=meteo arg3,4,5=dimensions

! Fermer le fichier netcdf:
      status=nf_close(ncid1)

      if(flag_meteodata=='glorys') then !---------->

        call airseaflux_sgrid_to_glorysgrid_driver ! glorys case
        call airseaflux_bouchetrou_driver          ! glorys case

      else                              !---------->

        call airseaflux_sgrid_to_meteogrid_perform('t',txtcase_)

! Identifier les points de la grille meteo e la fois en terre (sur grille meteo)
! et en mer (sur grille symphonie):
        call fichier_bouchetrou_meteo

        if(txtcase_=='current') then !-- Iterative Case -->
          if(par%rank==0)write(6,'(a,a)')        &
           'Refreshing interpolation procedure:' &
          ,' successful'
        endif                        !-- Iterative Case -->

      endif                             !---------->

      end subroutine airseaflux_extract_zone

!----------------------------------------------------------------------------

      subroutine airseaflux_sgrid_to_glorysgrid_driver
      implicit none
      integer :: flag_gridt_=0,flag_gridu_=0,flag_gridv_=0
      integer loop_
#ifdef synopsis
       subroutinetitle='airseaflux_sgrid_to_glorysgrid_driver'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


! Attention, GLORYS archive sur la grille C. Il y a donc 3 correspondances
! e etablir. On regarde les fichiers les un apres les autres, et progressivement
! les 3 flag_grid passent sur la valeur 1

      do loop_=1,nairsea
        open(unit=4,file=airseabinreclist(loop_) &
                   ,access='direct'              &
                   ,recl=540                     &
                   ,form='unformatted')
        read(4,rec=1)texte80(1),i
        close(4)


        if(loop_==ustrs_id) then !uuuuu>

         if(flag_gridu_==0) then
          flag_gridu_=1
          call airseaflux_sgrid_to_meteogrid_perform('u','glorys_u')
         endif

        endif                    !uuuuu>

        if(loop_==vstrs_id) then !vvvvv>

         if(flag_gridv_==0) then
          flag_gridv_=1
          call airseaflux_sgrid_to_meteogrid_perform('v','glorys_v')
         endif

        endif                    !vvvvv>

        if(loop_/=vstrs_id.and.loop_/=ustrs_id) then !ttttt>

         if(flag_gridt_==0) then
          flag_gridt_=1
          call airseaflux_sgrid_to_meteogrid_perform('t','glorys_t')
         endif

        endif                                        !ttttt>


        if(flag_gridt_==1.and.flag_gridu_==1.and.flag_gridv_==1)goto 2770

      enddo ! fin de boucle sur loop_

 2770 continue

      c_grid_shift=0
      if(lon_1_1_gridu> lon_1_1_gridt)c_grid_shift=-0.5
      if(lon_1_1_gridu< lon_1_1_gridt)c_grid_shift= 0.5

      if(lat_1_1_gridv> lat_1_1_gridt.and.c_grid_shift>=0.) &
      stop 'erreur1 c_grid_shift airseaflux_sgrid_to_glorysgrid_driver'
      if(lat_1_1_gridv< lat_1_1_gridt.and.c_grid_shift<=0.) &
      stop 'erreur2 c_grid_shift airseaflux_sgrid_to_glorysgrid_driver'
      if(lat_1_1_gridv==lat_1_1_gridt.and.c_grid_shift/=0.) &
      stop 'erreur3 c_grid_shift airseaflux_sgrid_to_glorysgrid_driver'

!     write(6,*)'c_grid_shift=',c_grid_shift
!     write(6,*)'lon_1_1_gridu lon_1_1_gridt=',lon_1_1_gridu,lon_1_1_gridt
!     write(6,*)'lat_1_1_gridv,lat_1_1_gridt=',lat_1_1_gridv,lat_1_1_gridt
!     stop'airseaflux_sgrid_to_glorysgrid_driver'

      end subroutine airseaflux_sgrid_to_glorysgrid_driver

!----------------------------------------------------------------------------

      subroutine airseaflux_sgrid_to_meteogrid_perform(text_,txtcase_)
      implicit none
      double precision deci_(0:1,0:1),decj_(0:1,0:1)    &
                        ,dy_(0:1,0:1),  dx_(0:1,0:1)    &
                        ,dlon_di_,dlon_dj_,dlon_dm_
      character*1 text_
      character(len=*)txtcase_
#ifdef synopsis
       subroutinetitle='airseaflux_sgrid_to_meteogrid_perform'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(status/=0)stop &
      'echec fichier airseaflux_sgrid_to_meteogrid_perform'

! lire les longitudes:
                   status=nf_inq_varid(ncid1,'lon',lon_id)
      if(status/=0)status=nf_inq_varid(ncid1,'longitude',lon_id)
      if(status/=0)status=nf_inq_varid(ncid1,'g0_lon_1',lon_id)
      if(status/=0)status=nf_inq_varid(ncid1,'nav_lon',lon_id)
      if(status/=0)stop 'Erreur lon_id variable meteo longitude'
      status=nf_inq_var(ncid1,lon_id,texte30,var_type,var_dims,tabdim,i4)
      if(status/=0)stop 'Erreur nf_inq_var meteo'

      if(var_type/=nf_float.and.var_type/=nf_double) &
      stop 'Erreur type sur longitude meteo'

! Recharger la longitude et la latitude depuis la zone reduite:
! Lire Longitude:

      if(var_dims==1) then
       i1=meteo_imax ; j1=1 ; varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax
      endif
      if(var_dims==2) then
       i1=meteo_imax ; j1=meteo_jmax
       varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax
       varstart(2)=meteozoom_jstr ; varcount(2)=meteo_jmax
      endif
      if(var_type==nf_double) then  !dddddddd>
       status=nf_get_vara_double(ncid1,lon_id,varstart(1:var_dims)   &
                                             ,varcount(1:var_dims)   &
                                             ,meteo_lon(1:i1,1:j1))
      endif                      !dddddddd>
      if(var_type==nf_real) then    !rrrrrrrr>
       allocate(meteo_lonlat_r4(meteo_imax,meteo_jmax))
       status=nf_get_vara_real(ncid1,lon_id,varstart(1:var_dims)      &
                                           ,varcount(1:var_dims)      &
                                           ,meteo_lonlat_r4(1:i1,1:j1))
       meteo_lon(1:i1,1:j1)=meteo_lonlat_r4(1:i1,1:j1)
       deallocate(meteo_lonlat_r4)
      endif                      !rrrrrrrr>
      if(status/=0)stop 'Erreur lecture longitude meteo'

      if(var_dims==1) then !ddddd> !19-03-15
       do j=1,meteo_jmax ; do i=1,meteo_imax !23-03-15
        meteo_lon(i,j)=meteo_lon(i,1)
       enddo             ; enddo
      endif                !ddddd>

! Lire Latitude:
                   status=nf_inq_varid(ncid1,'lat',lat_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude',lat_id)
      if(status/=0)status=nf_inq_varid(ncid1,'g0_lat_0',lat_id)
      if(status/=0)status=nf_inq_varid(ncid1,'nav_lat',lat_id)
      if(status/=0)stop 'Erreur lat_id variable meteo latitude'

      status=nf_inq_var(ncid1,lat_id,texte30,var_type,var_dims,tabdim,i4)

      if(var_type/=nf_float.and.var_type/=nf_double) &
      stop 'Erreur type sur latitude meteo'
      if(var_dims==1) then
       i1=1 ; j1=meteo_jmax ; varstart(1)=meteozoom_jstr ; varcount(1)=meteo_jmax
      endif
      if(var_dims==2) then
       i1=meteo_imax ; j1=meteo_jmax
       varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax
       varstart(2)=meteozoom_jstr ; varcount(2)=meteo_jmax
      endif

      if(var_type==nf_double) then  !dddddddd>
       status=nf_get_vara_double(ncid1,lat_id,varstart(1:var_dims)   &
                                             ,varcount(1:var_dims)   &
                                             ,meteo_lat(1:i1,1:j1))
      endif                      !dddddddd>
      if(var_type==nf_real) then    !rrrrrrrr>
       allocate(meteo_lonlat_r4(meteo_imax,meteo_jmax))
       status=nf_get_vara_real(ncid1,lat_id,varstart(1:var_dims)        &
                                           ,varcount(1:var_dims)        &
                                           ,meteo_lonlat_r4(1:i1,1:j1))
       meteo_lat(1:i1,1:j1)=meteo_lonlat_r4(1:i1,1:j1)
       deallocate(meteo_lonlat_r4)
      endif                      !rrrrrrrr>
      if(status/=0)stop 'Erreur lecture latitude meteo'

      if(var_dims==1) then !ddddd> !19-03-15
       do j=1,meteo_jmax ; do i=1,meteo_imax !23-03-15
        meteo_lat(i,j)=meteo_lat(1,j)
       enddo             ; enddo
      endif                !ddddd>

! Fermer le fichier netcdf:
      status=nf_close(ncid1)

      if(text_=='t') then
       lon_1_1_gridt=meteo_lon(1,1) ; lat_1_1_gridt=meteo_lat(1,1)
       lon_2_1_gridt=meteo_lon(2,1) ; lat_1_2_gridt=meteo_lat(1,2)
      endif
      if(text_=='u') then
       lon_1_1_gridu=meteo_lon(1,1) ; lat_1_1_gridu=meteo_lat(1,1)
      endif
      if(text_=='v') then
       lon_1_1_gridv=meteo_lon(1,1) ; lat_1_1_gridv=meteo_lat(1,1)
      endif

      if(text_=='u'.or.text_=='v')return

      if(.not.allocated(ij2meteo_i))   allocate(ij2meteo_i   (0:imax+1,0:jmax+1))
      if(.not.allocated(ij2meteo_j))   allocate(ij2meteo_j   (0:imax+1,0:jmax+1))
      if(.not.allocated(ij2meteo_teta))allocate(ij2meteo_teta(0:imax+1,0:jmax+1))


      do j=0,jmax+1                                                     !18-10-09
      do i=0,imax+1

      x2=real(meteo_imax/2)
      x3=real(meteo_jmax/2)
      deci=x2
      decj=x3

! ETAPE 1: trouver les coordonnees dans la grille ORCA:

! First guess: centre du domaine:
      k10=0
 1456 continue

! Principe suppose une relation lineaire entre lat lon et indice de grille.
!      dlon/di*Di+dlon/dj*Dj=Dlon
!      dlat/di*Di+dlat/dj*Dj=Dlat
! On cherche Di et Dj correspondant e Dlon=lon(i,j)-lonmeteo(i0,j0)
!                                et e Dlat=lat(i,j)-latmeteo(i0,j0)

      i1=int(deci)
      j1=int(decj)

      do j2=0,1
      do i2=0,1
      i0=i1+i2
      j0=j1+j2

      dlon_di_=meteo_lon(i0+1,j0  )-meteo_lon(i0-1,j0  )
      dlon_dj_=meteo_lon(i0  ,j0+1)-meteo_lon(i0  ,j0-1)
      dlon_dm_=rad2deg*lon_t(i  ,j)-meteo_lon(i0  ,j0)

      if(dlon_di_<-180.)dlon_di_=dlon_di_+360.
      if(dlon_di_> 180.)dlon_di_=dlon_di_-360.
      if(dlon_dj_<-180.)dlon_dj_=dlon_dj_+360.
      if(dlon_dj_> 180.)dlon_dj_=dlon_dj_-360.
      if(dlon_dm_<-180.)dlon_dm_=dlon_dm_+360.
      if(dlon_dm_> 180.)dlon_dm_=dlon_dm_-360.

! Determinant principal:
      x1=( dlon_di_                                             &
          *(meteo_lat(i0  ,j0+1)-meteo_lat(i0  ,j0-1))          &
          -(meteo_lat(i0+1,j0  )-meteo_lat(i0-1,j0  ))          &
          *dlon_dj_ )*0.25

      deci_(i2,j2)=min(max(                                  &
       i0+( dlon_dm_                                            &
          *(meteo_lat(i0  ,j0+1)-meteo_lat(i0,j0-1))            &
          -(rad2deg*lat_t(i,j)  -meteo_lat(i0,j0))              &
          *dlon_dj_)/x1*0.5    &
                   ,2.0001d0),meteo_imax-1.0001d0)

      decj_(i2,j2)=min(max(                                  &
       j0+( dlon_di_                                            &
          *(rad2deg*lat_t(i,j)  -meteo_lat(i0  ,j0))            &
          -(meteo_lat(i0+1,j0  )-meteo_lat(i0-1,j0))            &
          *dlon_dm_   )/x1*0.5    &
                   ,2.0001d0),meteo_jmax-1.0001d0)

      enddo
      enddo

      rapi=deci-i1
      rapj=decj-j1

      deci=(1.-rapi)*(1.-rapj)*deci_(0,0)   &
          +(1.-rapi)*    rapj *deci_(0,1)   &
          +    rapi *    rapj *deci_(1,1)   &
          +    rapi *(1.-rapj)*deci_(1,0)
      decj=(1.-rapi)*(1.-rapj)*decj_(0,0)   &
          +(1.-rapi)*    rapj *decj_(0,1)   &
          +    rapi *    rapj *decj_(1,1)   &
          +    rapi *(1.-rapj)*decj_(1,0)

! Si le point vise est different du first guess refaire le calcul
! avec un first guess donne par le dernier point vise:
      if(sqrt( (deci-x2)**2+(decj-x3)**2 ).gt.0.001)then
       x2=deci
       x3=decj
       k10=k10+1
       if(k10>20)then !!!!!!>
!       write(6,*)'par%rank=',par%rank
!       write(6,*)'(i,j)=   ',i,j
!       write(6,*)'deci decj',deci,decj
!       write(6,'(a,a)')'text_=',text_
        stop 'sgrid_to_glorysgrid_perform ne converge pas'
! Si le calcul ne converge pas car la grille de l'meteo est
! discontinue on peut essayer l'algo suivant: !21-06-15
         sum1=small1 ; sum2=0. ; sum3=0.
         do j3=1,meteo_jmax
! dist0 est une distance de reference mesuree entre (1,j3) et (2,j3) de l'meteo
          call lonlat2distance(meteo_lon(1,j3)*deg2rad,meteo_lat(1,j3)*deg2rad &
                              ,meteo_lon(2,j3)*deg2rad,meteo_lat(2,j3)*deg2rad,dist0)
         do i3=1,meteo_imax
          call lonlat2distance(lon_t(i,j),lat_t(i,j),meteo_lon(i3,j3)*deg2rad,meteo_lat(i3,j3)*deg2rad,dist1)
          sum1=sum1+exp(-dist1/dist0)
          sum2=sum2+exp(-dist1/dist0)*i3
          sum3=sum3+exp(-dist1/dist0)*j3
         enddo
         enddo
         deci=sum2/sum1 ; decj=sum3/sum1
         goto 4367
       endif          !!!!!!>
       goto 1456
      endif

! ETAPE 2: CALCULER L'ANGLE D'ORIENTATION LOCALE DE LA GRILLE ORCA

 4367 continue
      i1=int(deci) ; j1=int(decj)
      rapi=deci-i1 ; rapj=decj-j1
      do j2=0,1
      do i2=0,1
       i0=i1+i2
       j0=j1+j2

       dlon_di_=meteo_lon(i0+1,j0  )-meteo_lon(i0-1,j0  )
       if(dlon_di_<-180.)dlon_di_=dlon_di_+360.
       if(dlon_di_> 180.)dlon_di_=dlon_di_-360.

       dy_(i2,j2)= meteo_lat(i0+1,j0)-meteo_lat(i0-1,j0)
       dx_(i2,j2)=dlon_di_           &
                            *cos(meteo_lat(i0,j0)*deg2rad)
      enddo
      enddo
      x1=(1.-rapi)*(1.-rapj)*dx_(0,0)   &
        +(1.-rapi)*    rapj *dx_(0,1)   &
        +    rapi *    rapj *dx_(1,1)   &
        +    rapi *(1.-rapj)*dx_(1,0)
      y1=(1.-rapi)*(1.-rapj)*dy_(0,0)   &
        +(1.-rapi)*    rapj *dy_(0,1)   &
        +    rapi *    rapj *dy_(1,1)   &
        +    rapi *(1.-rapj)*dy_(1,0)

      ij2meteo_i(i,j)=deci+(meteozoom_istr-1) !21-03-15
      ij2meteo_j(i,j)=decj+(meteozoom_jstr-1)
      ij2meteo_teta(i,j)=atan2(y1,x1)

      enddo
      enddo

      call airseaflux_ij2meteo_mpi('za') !za=e0+z1 !23-04-15


! Dans le cas particulier de l'option "glorys" les fonctions de
! changement de repere sont basees sur des indices locaux:
      if(flag_meteodata=='glorys') then !glorys>
       ij2meteo_i(:,:)=ij2meteo_i(:,:)-(meteozoom_istr-1) !21-03-15
       ij2meteo_j(:,:)=ij2meteo_j(:,:)-(meteozoom_jstr-1)
      endif                             !glorys>

      end subroutine airseaflux_sgrid_to_meteogrid_perform

!..........................................................

      subroutine airseaflux_glorys_interp(t_,case_)
      implicit none
      integer loop_,t_,vstart1_,vstart2_,flag_rotation_,filval_,case_
      double precision time1_,time2_
      real i_shift_,j_shift_,filvalr4_
      character*1 text_
#ifdef synopsis
       subroutinetitle='airseaflux_glorys_interp'
       subroutinedescription= &
      'Reads meteo fields and interpolates on S grid'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      flag_rotation_=0

      do loop_=1,nairsea

      call airseaflux_decision(loop_,case_) ! donne feu vert pour lecture fichier !09-02-14

      if(decision==1) then                            !*******************>

      call airseaflux_get_time_from_binrecfile(loop_,t_,vstart1_,vstart2_,time1_,time2_)

! ouvrir le fichier netcdf:
      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(par%rank==0)write(6,'(a,a)')' FILE:  ',trim(texte80(1))
      if(status.ne.0)stop 'echec fichier 2'

      i_shift_=0. ; j_shift_=0.

      if(loop_==ustrs_id) then
        if(i_shift_==0) then ; text_='t' ; else ; text_='u' ; endif
        varcode(loop_)=1 ; i_shift_=c_grid_shift ; flag_rotation_=1
                   status=nf_inq_varid(ncid1,'sozotaux',var_id)
      endif
      if(loop_==vstrs_id) then
        if(i_shift_==0) then ; text_='t' ; else ; text_='v' ; endif
        varcode(loop_)=2 ; j_shift_=c_grid_shift
                   status=nf_inq_varid(ncid1,'sometauy',var_id)
      endif
      if(loop_==ssr_id) then
        text_='t' ; varcode(loop_)=3 ; ialbedo=0
                   status=nf_inq_varid(ncid1,'soceshwf',var_id)
      endif
      if(loop_==slhf_id) then
        text_='t' ; varcode(loop_)=4
                   status=nf_inq_varid(ncid1,'socelatf',var_id)
      endif
      if(loop_==netir_id) then
        text_='t' ; varcode(loop_)=5
                   status=nf_inq_varid(ncid1,'socelowf',var_id)
      endif
      if(loop_==sshf_id) then
        text_='t' ; varcode(loop_)=6
                   status=nf_inq_varid(ncid1,'socesenf',var_id)
      endif
      if(loop_==rain_id) then
        text_='t' ; varcode(loop_)=7
                   status=nf_inq_varid(ncid1,'sowaprec',var_id)
      endif

      if(status/=0) then
       write(6,'(a,a)')'Pour variable:',trim(airsea_standard_name(loop_))
       stop 'erreur varname airseaflux_glorys_interp'
      endif

      status=nf_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)
!     if(var_type/=nf_short) then !------->
!      write(6,*)'loop_=',loop_
!      stop &
!     'cas var_type pas prevu airseaflux_glorys_interp'
!     endif                       !------->

      varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax
      varstart(2)=meteozoom_jstr ; varcount(2)=meteo_jmax
      varstart(3)=vstart1_       ; varcount(3)=1
      varstart(4)=1              ; varcount(4)=1

      ksecu=0
      if(var_type==nf_short) then !ssssssss>
      ksecu=1

                   status=nf_get_att_int(ncid1,var_id,'_Fillvalue',filval_)
      if(status/=0)status=nf_get_att_int(ncid1,var_id,'_FillValue',filval_)
      if(status/=0)stop 'erreur _FillValue airseaflux_glorys_interp'
      status=nf_get_vara_int(ncid1,var_id,varstart(1:var_dims)        &
                                         ,varcount(1:var_dims)        &
                       ,meteo_short(1:meteo_imax,1:meteo_jmax))
      if(status/=0)stop &
      'erreur nf_get_vara_int airseaflux_glorys_interp'

      status=nf_get_att_real(ncid1,var_id,'scale_factor',var_scalefactor)
      if(status/=0) &
      stop 'error get scale_factor airseaflux_glorys_interp'
      status=nf_get_att_real(ncid1,var_id,'add_offset',var_addoffset)
      if(status/=0)stop 'error get add_offset irseaflux_flux_interp'

      texte90='' ; status=nf_get_att_text(ncid1,var_id,'units',texte90)
      if(status/=0)stop 'echec att units airseaflux_glorys_interp'

       filvalr4_=-9999.
       do j=1,meteo_jmax ; do i=1,meteo_imax
        meteo_var(i,j)=meteo_short(i,j)*var_scalefactor+var_addoffset
        if(meteo_short(i,j)==filval_)meteo_var(i,j)=filvalr4_
       enddo ; enddo
! Des que meteo_var est obtenu reset var_scalefactor=1. var_addoffset=0.
! pour etre compatible avec calculs plus tard
       var_scalefactor=1. ; var_addoffset=0.
      endif                       !ssssssss>
      if(var_type==nf_real) then  !rrrrrrrr>
      ksecu=1
                   status=nf_get_att_real(ncid1,var_id,'_Fillvalue',filvalr4_)
      if(status/=0)status=nf_get_att_real(ncid1,var_id,'_FillValue',filvalr4_)
      if(status/=0)stop 'erreur _FillValue airseaflux_glorys_interp'
      status=nf_get_vara_real(ncid1,var_id,varstart(1:var_dims)        &
                                          ,varcount(1:var_dims)        &
                          ,meteo_var(1:meteo_imax,1:meteo_jmax))
! A priori pas de var_scalefactor var_addoffset dans archivage real....:
      var_scalefactor=1. ; var_addoffset=0.
      endif                       !rrrrrrrr>
      if(ksecu==0)stop ' Erreur type 3905 module_airseaflux'

      texte90='' ; status=nf_get_att_text(ncid1,var_id,'units',texte90)
      if(status/=0)stop 'echec att units airseaflux_glorys_interp'
      if(varcode(loop_)==7.and.texte90=='kg m-2 s-1') &
      var_scalefactor=var_scalefactor*0.001

! Fermer le fichier netcdf
      status=nf_close(ncid1)

      call airseaflux_bouchetrou_apply(text_) ! glorys case
      call airseaflux_glorys_debug(filvalr4_) !12-03-13 ! BUG GLORYS MASQUE VARIABLE !!!!!????

      if(varcode(loop_)==1) then !----->
       do j=0,jmax+1 ; do i=0,imax+1
        deci=ij2meteo_i(i,j)+i_shift_ ; i1=int(deci) ; rapi=deci-i1
        decj=ij2meteo_j(i,j)+j_shift_ ; j1=int(decj) ; rapj=decj-j1
        wstress_u(i,j,0)=wstress_u(i,j,2)
        wstress_u(i,j,2)=((1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )  &
                         +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)  &
                         +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )  &
                         +    rapi *    rapj *meteo_var(i1+1,j1+1)) &
                         *var_scalefactor+var_addoffset
       enddo ; enddo
      endif             !----->

      if(varcode(loop_)==2) then !----->
       do j=0,jmax+1 ; do i=0,imax+1
        deci=ij2meteo_i(i,j)+i_shift_ ; i1=int(deci) ; rapi=deci-i1
        decj=ij2meteo_j(i,j)+j_shift_ ; j1=int(decj) ; rapj=decj-j1
        wstress_v(i,j,0)=wstress_v(i,j,2)
        wstress_v(i,j,2)=((1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )  &
                         +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)  &
                         +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )  &
                         +    rapi *    rapj *meteo_var(i1+1,j1+1)) &
                         *var_scalefactor+var_addoffset
       enddo ; enddo
      endif             !----->

      if(varcode(loop_)==3) then !----->
       do j=0,jmax+1 ; do i=0,imax+1
        deci=ij2meteo_i(i,j)+i_shift_ ; i1=int(deci) ; rapi=deci-i1
        decj=ij2meteo_j(i,j)+j_shift_ ; j1=int(decj) ; rapj=decj-j1
            ssr_w(i,j,0)=ssr_w(i,j,2)
            ssr_w(i,j,2)=((1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )  &
                         +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)  &
                         +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )  &
                         +    rapi *    rapj *meteo_var(i1+1,j1+1)) &
                         *var_scalefactor+var_addoffset
       enddo ; enddo
#ifdef parallele
      call get_type_echange('z0','ssr_w_z0_2',ssr_w,lbound(ssr_w),ubound(ssr_w),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(ssr_w,i2,mpi_neighbor_list(loop3)) !31-07-14
      enddo
      call loc_wait()
#endif
      if(flag_ssr24avr==1)call airseaflux_ssr24avr !30-07-14
      endif             !----->

      if(varcode(loop_)==4) then !----->
       do j=0,jmax+1 ; do i=0,imax+1
        deci=ij2meteo_i(i,j)+i_shift_ ; i1=int(deci) ; rapi=deci-i1
        decj=ij2meteo_j(i,j)+j_shift_ ; j1=int(decj) ; rapj=decj-j1
           slhf_w(i,j,0)=slhf_w(i,j,2)
           slhf_w(i,j,2)=((1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )  &
                         +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)  &
                         +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )  &
                         +    rapi *    rapj *meteo_var(i1+1,j1+1)) &
                         *var_scalefactor+var_addoffset
       enddo ; enddo
#ifdef parallele
      call get_type_echange('z0','slhf_w_z0_2',slhf_w,lbound(slhf_w),ubound(slhf_w),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(slhf_w,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()
#endif
      endif             !----->

      if(varcode(loop_)==5) then !----->
       do j=0,jmax+1 ; do i=0,imax+1
        deci=ij2meteo_i(i,j)+i_shift_ ; i1=int(deci) ; rapi=deci-i1
        decj=ij2meteo_j(i,j)+j_shift_ ; j1=int(decj) ; rapj=decj-j1
           snsf_w(i,j,0)=snsf_w(i,j,2)
           snsf_w(i,j,2)=((1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )  &
                         +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)  &
                         +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )  &
                         +    rapi *    rapj *meteo_var(i1+1,j1+1)) &
                         *var_scalefactor+var_addoffset
       enddo ; enddo
#ifdef parallele
      call get_type_echange('z0','snsf_w_z0_2',snsf_w,lbound(snsf_w),ubound(snsf_w),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(snsf_w,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()
#endif
      endif             !----->

      if(varcode(loop_)==6) then !----->
       do j=0,jmax+1 ; do i=0,imax+1
        deci=ij2meteo_i(i,j)+i_shift_ ; i1=int(deci) ; rapi=deci-i1
        decj=ij2meteo_j(i,j)+j_shift_ ; j1=int(decj) ; rapj=decj-j1
           sshf_w(i,j,0)=sshf_w(i,j,2)
           sshf_w(i,j,2)=((1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )  &
                         +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)  &
                         +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )  &
                         +    rapi *    rapj *meteo_var(i1+1,j1+1)) &
                         *var_scalefactor+var_addoffset
       enddo ; enddo
#ifdef parallele
      call get_type_echange('z0','sshf_w_z0_2',sshf_w,lbound(sshf_w),ubound(sshf_w),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(sshf_w,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()
#endif
      endif             !----->

      if(varcode(loop_)==7) then !----->
       do j=0,jmax+1 ; do i=0,imax+1
        deci=ij2meteo_i(i,j)+i_shift_ ; i1=int(deci) ; rapi=deci-i1
        decj=ij2meteo_j(i,j)+j_shift_ ; j1=int(decj) ; rapj=decj-j1
        precipi_w(i,j,0)=precipi_w(i,j,2)
        precipi_w(i,j,2)=((1.-rapi)*(1.-rapj)*meteo_var(i1  ,j1  )  &
                         +(1.-rapi)*    rapj *meteo_var(i1  ,j1+1)  &
                         +    rapi *(1.-rapj)*meteo_var(i1+1,j1  )  &
                         +    rapi *    rapj *meteo_var(i1+1,j1+1)) &
                         *var_scalefactor+var_addoffset
       enddo ; enddo
#ifdef parallele
      call get_type_echange('z0','precipi_w_z0_2',precipi_w,lbound(precipi_w),ubound(precipi_w),2,i2)
      do loop3=1,subcycle_exchange
        call echange_voisin(precipi_w,i2,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()
#endif
      endif             !----->


      endif                                           !*******************>    !

      enddo ! find de boucle sur loop_

! Une fois sortie de la boucle loop_ c'est le flag_rotation=1 qui nous
! dit de faire la rotation:
      if(flag_rotation_==1) then !++++++++++++++++++++++++++++>
! Passer des composantes alongaxis de la grille nemo
! aux direction zonales et meridiennes:
       do j=0,jmax+1
       do i=0,imax+1
        xy_t(i,j,1)=cos(ij2meteo_teta(i,j))*wstress_u(i,j,2) &
                   -sin(ij2meteo_teta(i,j))*wstress_v(i,j,2)

        xy_t(i,j,2)=sin(ij2meteo_teta(i,j))*wstress_u(i,j,2) &
                   +cos(ij2meteo_teta(i,j))*wstress_v(i,j,2)
       enddo
       enddo
! Passer des composantes zonales et meridiennes aux composantes
! alongaxis de la grille S:
       do j=0,jmax+1
       do i=0,imax+1
        xy_t(i,j,3)=gridrotcos_t(i,j)*xy_t(i,j,1)   &
                   -gridrotsin_t(i,j)*xy_t(i,j,2)

        xy_t(i,j,4)=gridrotsin_t(i,j)*xy_t(i,j,1)  &
                   +gridrotcos_t(i,j)*xy_t(i,j,2)
       enddo
       enddo
#ifdef parallele
      call get_type_echange('za','xy_t_za_3',xy_t,lbound(xy_t),ubound(xy_t),3,i3) !04-10-14
      call get_type_echange('za','xy_t_za_4',xy_t,lbound(xy_t),ubound(xy_t),4,i4)
      do loop3=1,subcycle_exchange
        call echange_voisin(xy_t,i3,mpi_neighbor_list(loop3)) !31-07-14
        call echange_voisin(xy_t,i4,mpi_neighbor_list(loop3))
      enddo
      call loc_wait()
#endif


! Interpolation point "u":
       do j=1,jmax
       do i=1,imax+1
        wstress_u(i,j,2)=0.5*(xy_t(i,j,3)+xy_t(i-1,j,3))
       enddo
       enddo
! Interpolation point "v":
       do j=1,jmax+1
       do i=1,imax
        wstress_v(i,j,2)=0.5*(xy_t(i,j,4)+xy_t(i,j-1,4))
       enddo
       enddo
      endif                      !++++++++++++++++++++++++++++>


      end subroutine airseaflux_glorys_interp

!..........................................................

      subroutine airseaflux_bouchetrou_driver
      implicit none
      integer :: flag_gridt_=0,flag_gridu_=0,flag_gridv_=0
      integer loop_
#ifdef synopsis
       subroutinetitle='airseaflux_bouchetrou_driver'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Attention, GLORYS archive sur la grille C. Il y a donc 3 correspondances
! e etablir. On regarde les fichiers les un apres les autres, et progressivement
! les 3 flag_grid passent sur la valeur 1

      do loop_=1,nairsea
        open(unit=4,file=airseabinreclist(loop_) &
                   ,access='direct'              &
                   ,recl=540                     &
                   ,form='unformatted')
        read(4,rec=1)texte80(1),i
        close(4)

        if(flag_gridt_==0) then !-------->
        if(index(texte80(1),'gridt')/=0.or. &
           index(texte80(1),'gridT')/=0.or. & !13-03-14
           index(texte80(1),'GRIDT')/=0) then !ooooo>
           flag_gridt_=1
           call airseaflux_bouchetrou_perform('t')
        endif                                 !ooooo>
        endif                   !-------->

        if(flag_gridu_==0) then !-------->
        if(index(texte80(1),'gridu')/=0.or. &
           index(texte80(1),'gridU')/=0.or. & !13-03-14
           index(texte80(1),'GRIDU')/=0) then !ooooo>
           flag_gridu_=1
           call airseaflux_bouchetrou_perform('u')
        endif                                 !ooooo>
        endif                   !-------->

        if(flag_gridv_==0) then !-------->
        if(index(texte80(1),'gridv')/=0.or. &
           index(texte80(1),'gridV')/=0.or. & !13-03-14
           index(texte80(1),'GRIDV')/=0) then !ooooo>
           flag_gridv_=1
           call airseaflux_bouchetrou_perform('v')
        endif                                 !ooooo>
        endif                   !-------->

        if(flag_gridt_==1.and.flag_gridu_==1.and.flag_gridv_==1)goto 2770

      enddo ! fin de boucle sur loop_

      if(flag_gridt_+flag_gridu_+flag_gridv_==0) then !00000>
      write(6,'(a,a)')'Apparement le nom du fichier meteo ne contient ' &
      ,'pas la suite de caracteres attendue: gridt GRIDT gridu etc...'
      stop ' STOP airseaflux_bouchetrou_driver'
      endif                                           !00000>

 2770 continue

      end subroutine airseaflux_bouchetrou_driver

!.........................................................................

      subroutine airseaflux_bouchetrou_perform(text_)
      implicit none
      character*1 text_
      integer filval_
      real filvalr4_
#ifdef synopsis
       subroutinetitle='airseaflux_bouchetrou_perform'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! ATTENTION ON NE DOIT PASSER PAR CET ROUTINE QUE POUR LES CHAMPS GLORYS
! Contrairement e un vrai modele meteo qui fournit un masque terre mer
! ici c'est la grille nemo qui fait valeur de masque meteo. On regarde
! si la valeur de la variable meteo correspond au "missing value"
      if(flag_meteodata/='glorys')stop &
      'erreur airseaflux_bouchetrou_perform'

      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(status/=0)stop 'echec fichier routine fichier_bouchetrou_meteo'

      if(text_=='u') then
                    status=nf_inq_varid(ncid1,'sozotaux',var_id)
       if(status/=0)stop 'erreur sozotaux airseaflux_bouchetrou_perform'
      endif
      if(text_=='v') then
                    status=nf_inq_varid(ncid1,'sometauy',var_id)
       if(status/=0)stop 'erreur sometauy airseaflux_bouchetrou_perform'
      endif
      if(text_=='t') then
!                   status=nf_inq_varid(ncid1,'soceshwf',var_id)
!                   status=nf_inq_varid(ncid1,'soceshwf',var_id)
!                   status=nf_inq_varid(ncid1,'socelatf',var_id)
                    status=nf_inq_varid(ncid1,'socelowf',var_id)
!                   status=nf_inq_varid(ncid1,'socesenf',var_id)
       if(status/=0)stop 'erreur soceshwf airseaflux_bouchetrou_perform'
      endif

      varstart(1)=meteozoom_istr ; varcount(1)=meteo_imax
      varstart(2)=meteozoom_jstr ; varcount(2)=meteo_jmax
      varstart(3:4)=1            ; varcount(3:4)=1

      status=nf_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)

      ksecu=0
      if(var_type==nf_short) then !ssssssss>
       ksecu=1
                    status=nf_get_att_int(ncid1,var_id,'_Fillvalue',filval_)
       if(status/=0)status=nf_get_att_int(ncid1,var_id,'_FillValue',filval_)
       if(status/=0) &
       stop 'erreur _FillValue airseaflux_bouchetrou_perform'

       status=nf_get_vara_int(ncid1,var_id,varstart(1:var_dims)        &
                                          ,varcount(1:var_dims)        &
                        ,meteo_short(1:meteo_imax,1:meteo_jmax))
       filvalr4_=-9999.
       do j=1,meteo_jmax ; do i=1,meteo_imax
        meteo_var(i,j)=meteo_short(i,j)
        if(meteo_short(i,j)==filval_)meteo_var(i,j)=filvalr4_
       enddo ; enddo
      endif                       !ssssssss>
      if(var_type==nf_real)  then !rrrrrrrr>
       ksecu=1
                    status=nf_get_att_real(ncid1,var_id,'_Fillvalue',filvalr4_)
       if(status/=0)status=nf_get_att_real(ncid1,var_id,'_FillValue',filvalr4_)
       if(status/=0) &
       stop 'erreur nf_get_att_real _Fillvalue filvalr4_'
       status=nf_get_vara_real(ncid1,var_id,varstart(1:var_dims)       &
                                           ,varcount(1:var_dims)       &
                           ,meteo_var(1:meteo_imax,1:meteo_jmax))
      endif                       !rrrrrrrr>

      if(ksecu==0)stop 'type inconnu 4175 module_airseaflux'
      if(status/=0)stop &
      'erreur nf_get_vara_int airseaflux_bouchetrou_perform'

      status=nf_close(ncid1)

      x0=0. ; y0=0.
      if(text_=='u')x0=c_grid_shift
      if(text_=='v')y0=c_grid_shift
      do j=0,jmax+1
      do i=0,imax+1
      if(mask_t(i,j,kmax+1)==1) then !----->

      deci=ij2meteo_i(i,j)+x0 ; i1=int(deci) ; rapi=deci-i1
      decj=ij2meteo_j(i,j)+y0 ; j1=int(decj) ; rapj=decj-j1

      if(meteo_var(i1  ,j1  )==filvalr4_)meteo_var(i1  ,j1  )=-filvalr4_
      if(meteo_var(i1  ,j1+1)==filvalr4_)meteo_var(i1  ,j1+1)=-filvalr4_
      if(meteo_var(i1+1,j1+1)==filvalr4_)meteo_var(i1+1,j1+1)=-filvalr4_
      if(meteo_var(i1+1,j1  )==filvalr4_)meteo_var(i1+1,j1  )=-filvalr4_

      endif                          !----->
      enddo
      enddo

      if(text_=='t') &
      open(unit=4,file=trim(tmpdirname)//'bouchetrou_meteo_gridt'//dom_c//'.out')
      if(text_=='u') &
      open(unit=4,file=trim(tmpdirname)//'bouchetrou_meteo_gridu'//dom_c//'.out')
      if(text_=='v') &
      open(unit=4,file=trim(tmpdirname)//'bouchetrou_meteo_gridv'//dom_c//'.out')

      do 1862 j=1,meteo_jmax
      do 1862 i=1,meteo_imax

       if(meteo_var(i,j)==-filvalr4_)then !eeeeeeeeeeeeeee>

        ksecu=0
        i10=1
        dist1=1e10
 1864   continue

        do 1863 k1=0,1

         j0=max(j-i10,1)
         j2=min(j+i10,meteo_jmax)
         j3=k1*(j2-j0)+(1-k1)

         i0=max(i-i10+k1,1)
         i2=min(i+i10-k1,meteo_imax)
         i3=(1-k1)*(i2-i0)+k1

         do 1861 j1=j0,j2,j3
         do 1861 i1=i0,i2,i3

         if(i1<1         )stop 'faute meteo fp5'  !01-07-10
         if(i1>meteo_imax)stop 'faute meteo fp6'
         if(j1<1         )stop 'faute meteo fp7'
         if(j1>meteo_jmax)stop 'faute meteo fp8'
! si le modele plante sur fautes 5 6 7 ou 8 c'est qu'il n'y a pas de bouchage
! suffisament proche du trou et que l'algo cherche au dele de la zone d'extraction,
! ce qu'on ne permet pas pour avoir une parallelisation parfaite.
! Ce qu'on peu faire: on peut augmenter la taille de la
! zone d'extraction en augmentant i0 (repere 1498). Mais cela revele surtout
! que les masques oceano et meteo sont trop differents et qu'il faut sans doute
! ameliorer le masque oceano

          if(     meteo_var(i1,j1)/=-filvalr4_                      &
             .and.meteo_var(i1,j1)/= filvalr4_)then !%%%%%%%%%%%%%>
           dist2=sqrt(real(i-i1)**2+real(j-j1)**2)
           ksecu=1
             if(dist2.lt.dist1) then                    !>>>>>>>>>>>>>
              i4=i1
              j4=j1
              dist1=dist2
             endif                                      !>>>>>>>>>>>>>
          endif                                      !%%%%%%%%%%%%%>
 1861    continue

 1863   continue
        i10=i10+1
        if(ksecu.eq.0)goto 1864
                 write(4,'(4(i4,1x))')i,j  ,i4,j4

       endif                           !eeeeeeeeeeeeeee>
 1862 continue

      close(4)

      end subroutine airseaflux_bouchetrou_perform

!.............................................................

      subroutine airseaflux_bouchetrou_apply(text_)
      implicit none
      character*1 text_
#ifdef synopsis
       subroutinetitle='airseaflux_bouchetrou_apply'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(text_=='t') &
      open(unit=3,file=trim(tmpdirname)//'bouchetrou_meteo_gridt'//dom_c//'.out')
      if(text_=='u') &
      open(unit=3,file=trim(tmpdirname)//'bouchetrou_meteo_gridu'//dom_c//'.out')
      if(text_=='v') &
      open(unit=3,file=trim(tmpdirname)//'bouchetrou_meteo_gridv'//dom_c//'.out')

      do j1=1,meteo_jmax
      do i1=1,meteo_imax
       read(3,*,end=1760)i,j,i4,j4
!      meteo_short(i,j)=meteo_short(i4,j4)
       meteo_var(i,j)=meteo_var(i4,j4)
      enddo
      enddo
 1760 continue

      close(3)

      end subroutine airseaflux_bouchetrou_apply

!.............................................................

      subroutine airseaflux_flux_moveforward ! time linear interpolation
      implicit none
      integer loop_
#ifdef synopsis
       subroutinetitle='airseaflux_flux_moveforward'
       subroutinedescription='Linear time interpolation of meteo fluxes'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      do loop_=1,nairsea

      x2=(elapsedtime_now          -airseafile_prvtime(loop_))      &
        /(airseafile_nextime(loop_)-airseafile_prvtime(loop_))
! Bornage car l'attente de la synchro subcycling peut faire sortir (legerement) de l'intervale [O:1]
      x2=min(max(x2,0.d00),1.d00) !09-02-14
      x0=1.-x2

      if(varcode(loop_)==1) then
       do j=1,jmax ; do i=1,imax+1
        wstress_u(i,j,1)=x0*wstress_u(i,j,0)+x2*wstress_u(i,j,2)
       enddo ; enddo
      endif

      if(varcode(loop_)==2) then
       do j=1,jmax+1 ; do i=1,imax
        wstress_v(i,j,1)=x0*wstress_v(i,j,0)+x2*wstress_v(i,j,2)
       enddo ; enddo
      endif

      if(varcode(loop_)==3) then
       do j=1,jmax ; do i=1,imax
        ssr_w(i,j,1)=x0*ssr_w(i,j,0)+x2*ssr_w(i,j,2)
       enddo ; enddo
      endif

      if(varcode(loop_)==4) then
       do j=1,jmax ; do i=1,imax
        slhf_w(i,j,1)=x0*slhf_w(i,j,0)+x2*slhf_w(i,j,2)
       enddo ; enddo
      endif

      if(varcode(loop_)==5) then
       do j=1,jmax ; do i=1,imax
        snsf_w(i,j,1)=x0*snsf_w(i,j,0)+x2*snsf_w(i,j,2)
       enddo ; enddo
      endif

      if(varcode(loop_)==6) then
       do j=1,jmax ; do i=1,imax
        sshf_w(i,j,1)=x0*sshf_w(i,j,0)+x2*sshf_w(i,j,2)
       enddo ; enddo
      endif

      if(varcode(loop_)==7) then
       do j=1,jmax ; do i=1,imax
        precipi_w(i,j,1)=x0*precipi_w(i,j,0)+x2*precipi_w(i,j,2)
       enddo ; enddo
      endif

      enddo ! fin de boucle sur loop_

      if(flag_1dv==1)call airseaflux_1dv_glorys !10-3-16 LEo

! rebuild the total stress magnitude for the tKe equation !21-05-10
      do j=1,jmax ; do i=1,imax
       wstress_w(i,j)=sqrt( (0.5*( wstress_u(i  ,j  ,1)          &
                                  +wstress_u(i+1,j  ,1)))**2     &
                           +(0.5*( wstress_v(i  ,j  ,1)          &
                                  +wstress_v(i  ,j+1,1)))**2 )
      enddo ; enddo

      end subroutine airseaflux_flux_moveforward ! time linear interpolation

!.............................................................

      subroutine airseaflux_glorys_debug(filval_) ! BUG GLORYS MASQUE VARIABLE !!!!!????
      implicit none
!     integer filval_
      real filval_
#ifdef synopsis
       subroutinetitle='airseaflux_glorys_debug'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Le masque filval de GLORYS est apparement bugue puisque fluctuant !!??
! Ceci prend en defaut la regle de bouchage des trous de notre algo
! On applique une roue de secours pour eviter les aberations.
! En pratique on met a zero les valeurs masquees
      do j=1,meteo_jmax
      do i=1,meteo_imax
!      if(meteo_short(i,j)==filval_) then
!        meteo_short(i,j)=0
!        write(6,*)'pb en ',par%rank,i,j,filval_
!        stop 'coucou'
!      endif
       if(meteo_var(i,j)==filval_)meteo_var(i,j)=0.
      enddo
      enddo

      end subroutine airseaflux_glorys_debug ! BUG GLORYS MASQUE VARIABLE !!!!!????

!.............................................................

      subroutine airseaflux_readnotebook
      implicit none
      integer ncid_ &
             ,vartype_ !17-09-21
      character*456 txt19_
      namelist/notebook_airseaflux/iairsea,flag_meteodata,airseaoption &!30-07-14
                                  ,irelaxsst,flag_meteo_average        &
                                  ,texte90,airseafile,flag_ssr24avr    &
                                  ,flag_abl,meteo_t0                   &
                                  ,flag_lsm,flag_p0m_filter            &
                                  ,relativewind,flag_abl2              &!24-11-15
                                  ,meteolandconvention                 &!29-04-16
                                  ,bulk_scheme                         &!02-05-16
                                  ,flag_gridoverflow                   &
                                  ,flag_meteo_land_plug                &!13-07-18
                                  ,flag_meteo_land_plug_wind           &!13-07-18
                                  ,flag_wstressbulk                    &!25-01-18
                                  ,meteo_cumul_modulo                  &!28-10-118
                                  ,meteo_grid_file                      
#ifdef synopsis
       subroutinetitle='airseaflux_readnotebook'
       subroutinedescription='Reads notebook_airseaflux'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Initialiation de quelques constantes necessaires A plusieurs
! fonctionnalites du modele:
      call initial_bulk_moon(0,0) ! Arguments sans importance !19-05-16

      iairsea=0 ; irelaxsst=0 ; ialbedo=0 ; nairsea=0 ; airseafile='s'
      flag_meteodata='ecmwf' ; flag_meteo_average=0 ; texte90='s'
      airseaoption=1 ; flag_ssr24avr=0 ; meteo_t0='file' ; flag_lsm=0
      flag_p0m_filter=0 ; relativewind=0. ; flag_abl2=0 ; flag_gridoverflow=0
      meteolandconvention=1                                             !29-04-16
      bulk_scheme=bulk_core                                             !02-05-16
      flag_wstressbulk=1                                                !25-01-18
      flag_meteo_land_plug=1                                            !13-07-18
      flag_meteo_land_plug_wind=-99 !27-12-21
       open(100,file=nomfichier(7))
       read(100,nml=notebook_airseaflux)
       close(100)
! Si aucune valeur pour flag_meteo_land_plug_wind n'est donnee dans
! notebook_airseaflux, par defaut on prend la valeur de flag_meteo_land_plug:
      if(flag_meteo_land_plug_wind==-99) &
         flag_meteo_land_plug_wind=flag_meteo_land_plug !27-12-21


      if(meteo_t0=='user'.and.meteo_cumul_modulo==0) & !28-10-18
      stop 'Err 5732 meteo_t0=="user".and.meteo_cumul_modulo==0'

      if(iairsea==0)return

      if(flag_abl2==1.and.flag_abl==0) stop &
      'airseaflux flag_abl2==1.and.flag_abl==0 inconsistent choice'
!     if(flag_abl2==1.and.relativewind/=1.) stop &
!     'airseaflux flag_abl2==1.and.relativewind/=1. inconsistent choice'
      

      if(flag_meteodata=='ecmwf')  then !eee>
       nairsea=8 ; iairsea=2
       if(flag_wstressbulk==0)nairsea=10 ! 2 champs de plus si wstress pas bulk !25-01-18
      endif                             !eee>
      if(flag_meteodata=='arome')  then !eee> !08-11-16
       nairsea=8 ; iairsea=2
       if(flag_wstressbulk==0)nairsea=10 ! 2 champs de plus si wstress pas bulk !25-01-18
      endif                             !eee>
      if(flag_meteodata=='glorys') then !ggg>
       nairsea=7 ; iairsea=1
      endif                             !ggg>
      if(nairsea==0)stop 'flag_meteodata not recognised'

      do k=1,nairsea
       airseafile(k)=trim(texte90)//trim(airseafile(k))
       if(par%rank==0)write(6,'(a,i0,a,a)')'airseafile(',k,')=',trim(airseafile(k))
      enddo

! Detection des dates min et max pour allocation de certains tableaux
      do loop1=1,nairsea
       open(unit=3,file=airseafile(loop1))

       do loop2=1,2
         if(loop2==1) then ! Premier fichier:
          read(3,'(a)')texte90
         endif
         if(loop2==2) then ! Dernier fichier
 100      read(3,'(a)',end=200)texte90 ; goto 100
         endif
 200     continue

         status=nf_open(trim(texte90),nf_nowrite,ncid_)
         if(status/=0) then
          write(6,'(a,a)')' Echec ouverture fichier ',trim(texte90)
          stop ' error nf_open in airseaflux_readnotebook'
         endif
                      status=nf_inq_dimid(ncid_,'time',var_id)
         if(status/=0)status=nf_inq_dimid(ncid_,'time_counter',var_id)
         if(status/=0)status=nf_inq_dimid(ncid_,'Time',var_id) !17-09-21
         if(status/=0)stop 'erreur nf_inq_dimid time airsea'

         status=nf_inq_dimlen(ncid_,var_id,max_time_counter) 
         if(status/=0)stop 'erreur nf_inq_dimlen max_time_counter'

                      status=nf_inq_varid(ncid_,'time',var_id)
         if(status/=0)status=nf_inq_varid(ncid_,'time_counter',var_id)
         if(status/=0)status=nf_inq_varid(ncid_,'Times',var_id) !17-09-21
         if(status/=0)stop 'erreur nf_inq_varid time airsea'

         status=nf_inq_vartype(ncid_,var_id,vartype_) !17-09-21
         if(status/=0)stop 'erreur nf_inq_vartype time'

         texte60=''
             status=nf_get_att_text(ncid_,var_id,'units',texte60)
          if(status/=0)stop 'erreur nf_get_att_text airsea units'

         k=index(texte60,' since ')
         read(texte60(k+7:k+10),*)k0
         airsea_year_min=min(airsea_year_min,k0)
         airsea_year_max=max(airsea_year_max,k0)

         status=nf_close(ncid_)

        enddo !loop2

       close(3)
      enddo !loop1

      end subroutine airseaflux_readnotebook

!.............................................................

      subroutine airseaflux_allocate
      implicit none
#ifdef synopsis
       subroutinetitle='airseaflux_allocate'
       subroutinedescription='Air-sea fluxes arrays allocation'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Dimensions des tableaux sshf et slhf
      if(ifb==0) then
          k1=0 ; k2=2   ! Cas flux precalcules
      endif
      if(ifb==1) then
       if(flag_abl==0) then
          k1=1 ; k2=1  ! Cas formule bulk
       else
          k1=0 ; k2=2  ! Cas formule bulk + modele de couche limite atmospherique
       endif
      endif
      if(allocated(sshf_w))    deallocate(sshf_w)    ; allocate(sshf_w    (0:imax+1,0:jmax+1,k1:k2)  ) ; sshf_w=0
      if(allocated(slhf_w))    deallocate(slhf_w)    ; allocate(slhf_w    (0:imax+1,0:jmax+1,k1:k2)  ) ; slhf_w=0

      k3=k1 ; k4=k2
      if(ifb==1.and.flag_wstressbulk==0) then !m0v0m> !25-01-18
          k3=0 ; k4=2   ! Cas wstress donnE par le modele meteo et non pas calculE par les formules bulk
      endif                                   !m0v0m>
      if(allocated(wstress_u)) deallocate(wstress_u) ; allocate(wstress_u (0:imax+1,0:jmax+1,k3:k4)  ) ; wstress_u=0 !24-11-15
      if(allocated(wstress_v)) deallocate(wstress_v) ; allocate(wstress_v (0:imax+1,0:jmax+1,k3:k4)  ) ; wstress_v=0 !24-11-15

      if(ifb==1) then !11111>
      if(allocated(q10_t))     deallocate(q10_t)      ; allocate(q10_t     (0:imax+1,0:jmax+1)    ) ; q10_t=0
      if(allocated(teta10_t))  deallocate(teta10_t)   ; allocate(teta10_t  (0:imax+1,0:jmax+1)    ) ; teta10_t=0
      if(allocated(teta2_t))   deallocate(teta2_t)    ; allocate(teta2_t   (0:imax+1,0:jmax+1,0:2)) ; teta2_t=0
      if(allocated(q2_t))      deallocate(q2_t)       ; allocate(q2_t      (0:imax+1,0:jmax+1,0:2)) ; q2_t=0
      if(allocated(tetastar_t))deallocate(tetastar_t) ; allocate(tetastar_t(0:imax+1,0:jmax+1,0:2)) ; tetastar_t=0
      if(allocated(ustar_t))   deallocate(ustar_t)    ; allocate(ustar_t   (0:imax+1,0:jmax+1,0:2)) ; ustar_t=0
      if(allocated(qstar_t))   deallocate(qstar_t)    ; allocate(qstar_t   (0:imax+1,0:jmax+1,0:2)) ; qstar_t=0
      endif           !11111>

      if(allocated(airseafile_prvtime))deallocate(airseafile_prvtime); allocate(airseafile_prvtime(dim_airsea)) ;airseafile_prvtime=0.
      if(allocated(airseafile_nextime))deallocate(airseafile_nextime); allocate(airseafile_nextime(dim_airsea)) ;airseafile_nextime=0.
      if(allocated(airseafile_prvtrec))deallocate(airseafile_prvtrec); allocate(airseafile_prvtrec(dim_airsea)) ;airseafile_prvtrec=0
      if(allocated(airseafile_nextrec))deallocate(airseafile_nextrec); allocate(airseafile_nextrec(dim_airsea)) ;airseafile_nextrec=0
      if(allocated(airseabinreclist))  deallocate(airseabinreclist)  ; allocate(airseabinreclist  (dim_airsea)) ;airseabinreclist='reset'
      if(allocated(kz_abl_w))  deallocate(kz_abl_w) ; allocate(kz_abl_w   (  imax  ,  jmax  ,0:2  )  ) ; kz_abl_w=0.

      if(flag_abl==1) then !ooo>

!     if(allocated(uwindabl_t))  deallocate(uwindabl_t)  ; allocate(uwindabl_t  ( 0:imax+1, 0:jmax+1,1:kmax_abl,1:2)) ; uwindabl_t=0.
!     if(allocated(vwindabl_t))  deallocate(vwindabl_t)  ; allocate(vwindabl_t  ( 0:imax+1, 0:jmax+1,1:kmax_abl,1:2)) ; vwindabl_t=0.
      if(allocated(uwindabl_t))  deallocate(uwindabl_t)  ; allocate(uwindabl_t  ( 0:imax+1, 0:jmax+1,1         ,1:2)) ; uwindabl_t=0.
      if(allocated(vwindabl_t))  deallocate(vwindabl_t)  ; allocate(vwindabl_t  ( 0:imax+1, 0:jmax+1,1         ,1:2)) ; vwindabl_t=0.
      if(allocated(wwindabl_w))  deallocate(wwindabl_w)  ; allocate(wwindabl_w  ( 0:imax+1, 0:jmax+1,1:kmax_abl+1  )) ; wwindabl_w=0.
      if(allocated(ablheight_t)) deallocate(ablheight_t) ; allocate(ablheight_t ( 0:imax+1, 0:jmax+1,1:2)) ; ablheight_t=0.

      if(allocated(zw_abl))deallocate(zw_abl);allocate(zw_abl(kmax_abl));zw_abl=0.
      do k=2,kmax_abl ; zw_abl(k)=zw_abl(k-1)+50. ; enddo

      if(allocated(teta0_t))     deallocate(teta0_t)     ; allocate(teta0_t     ( 0:imax+1, 0:jmax+1,1:2)) ; teta0_t=0.
      if(allocated(teta2delta_t))deallocate(teta2delta_t); allocate(teta2delta_t( 1:imax  , 1:jmax  ,kmax_abl))   ; teta2delta_t=0.
      if(allocated(q2delta_t))   deallocate(q2delta_t)   ; allocate(q2delta_t   ( 1:imax  , 1:jmax  ,kmax_abl))   ; q2delta_t=0.
      if(allocated(dt_abl_max))  deallocate(dt_abl_max)  ; allocate(dt_abl_max  (0:2))                     ; dt_abl_max=0.
        
      if(flag_abl2==1) then !222>
      if(allocated(uwinddelta_t))deallocate(uwinddelta_t) ; allocate(uwinddelta_t( 1:imax  , 1:jmax  ,kmax_abl))   ; uwinddelta_t=0.
      if(allocated(vwinddelta_t))deallocate(vwinddelta_t) ; allocate(vwinddelta_t( 1:imax  , 1:jmax  ,kmax_abl))   ; vwinddelta_t=0.
      endif                 !222>

      endif                !ooo>

      end subroutine airseaflux_allocate

!.............................................................

      subroutine airseaflux_varidentifier
      implicit none
#ifdef synopsis
       subroutinetitle='airseaflux_varidentifier'
       subroutinedescription='Defines the meteo variable identifiers'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! reset meteovarname
        meteovarname='undefined' !26-01-22

       do k=1,nairsea

        ksecu=0 !19-01-17

        k1=len_trim(airseafile(k))
!       write(6,'(a)')airseafile(k)
        if(airseafile(k)(k1-3:k1)=='_ssr')   then ; ssr_id=k   ; ksecu=1 ; endif
        if(airseafile(k)(k1-2:k1)=='_ir')    then ; ir_id=k    ; ksecu=1 ; endif
        if(airseafile(k)(k1-4:k1)=='_rain')  then ; rain_id=k  ; ksecu=1 ; endif
        if(airseafile(k)(k1-3:k1)=='_t2m')   then ; t2m_id=k   ; ksecu=1 ; endif
        if(airseafile(k)(k1-4:k1)=='_dp2m')  then ; dp2m_id=k  ; ksecu=1 ; endif
        if(airseafile(k)(k1-4:k1)=='_u10m')  then ; u10m_id=k  ; ksecu=1 ; endif
        if(airseafile(k)(k1-4:k1)=='_v10m')  then ; v10m_id=k  ; ksecu=1 ; endif
        if(airseafile(k)(k1-3:k1)=='_p0m')   then ; p0m_id=k   ; ksecu=1 ; endif
        if(airseafile(k)(k1-5:k1)=='_ustrs') then ; ustrs_id=k ; ksecu=1 ; endif
        if(airseafile(k)(k1-5:k1)=='_vstrs') then ; vstrs_id=k ; ksecu=1 ; endif
        if(airseafile(k)(k1-4:k1)=='_slhf')  then ; slhf_id=k  ; ksecu=1 ; endif
        if(airseafile(k)(k1-5:k1)=='_netir') then ; netir_id=k ; ksecu=1 ; endif
        if(airseafile(k)(k1-4:k1)=='_sshf')  then ; sshf_id=k  ; ksecu=1 ; endif

         if(par%rank==0.and.ksecu==0) then !>>>>> !19-03-15
          write(6,'(a,a)')'Invalid list name: ',trim(airseafile(k))
          write(6,'(a)')'Last characters do not respect the coding rule'
          stop ' subroutine airseaflux_varidentifier'
         endif                             !>>>>>

       enddo ! k

       u100m_id=k   ; k=k+1
       v100m_id=k   ; k=k+1
         t0m_id=k   ; k=k+1
         abl_id=k

!     write(6,*)ssr_id,ir_id,rain_id,t2m_id,dp2m_id,u10m_id,v10m_id,p0m_id
!     write(6,*)ssr_id,ustrs_id,vstrs_id,slhf_id,netir_id,sshf_id,rain_id

!     stop 'sof'

      end subroutine airseaflux_varidentifier

!.............................................................

      subroutine airseaflux_relaxed_sst
      implicit none
#ifdef synopsis
       subroutinetitle='airseaflux_relaxed_sst'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      x1=rho*cp    &
        /(8.*86400.) ! echelle de temps de rappel (en secondes)
!     x0=1.-rap_obc                                                    !19-03-10
      x0=1.-timeweightobc(trc_id)                                      !11-07-14
      x2=1.-x0                                                         !19-03-10
       do j=1,jmax
       do i=1,imax
        heatrelax_w(i,j,1)=( x0*temobc_t(i,j,kmax,0)  &
                            +x2*temobc_t(i,j,kmax,2)  &
                                  -tem_t(i,j,kmax,0)) &
                                *x1*dz_t(i,j,kmax,0)

       enddo
       enddo

      end subroutine airseaflux_relaxed_sst

!.............................................................

      subroutine airseaflux_decision(loop_,case_)
      implicit none
      integer loop_,case_
#ifdef synopsis
       subroutinetitle='airseaflux_decision'
       subroutinedescription= &
          'Detects if the time-interpolation between to consecutive'   &
       //' meteo fields will pass the "after" field and thus requires' &
       //' to read a new meteo file. In the former case decision=1'    &
       //' (otherwise decision=0)'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! decision=1 indique qu'il est temps de lire un nouveau fichier. decision=0 sinon.

      decision=0

! Si les proc ne sont pas synchro en raison du subcycling on s'abstient afin
! de garantir la continuite horizontale des champs interpoles:
      if(subcycle_synchro==0)return

! Jalon temps depasse = lire un nouveau fichier
      if(elapsedtime_now>airseafile_nextime(loop_))decision=1

! A l'etat initial lire imperativement
      if(case_==1)decision=1 ! Etat initial

      end subroutine airseaflux_decision

!.............................................................

      subroutine airseaflux_inquire_albedo(txtcase_)
      implicit none
      integer ncid_,varid_
      character(len=*)txtcase_
      
#ifdef synopsis
       subroutinetitle='airseaflux_inquire_albedo'
       subroutinedescription= &
         'Determines if albebo should be applied on solar flux'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Include albedo ?
! Var id = ssr_id
!     if(txtcase_=='initial') then   !--initial case-->

       nc=airseafile_nextrec(ssr_id) ! Get netcdf fil name
       if(par%rank==0) then
        write(6,'(a)')'Subroutine airseaflux_inquire_albedo:'
        write(6,'(a,a)')'Open file ',trim(airseabinreclist(ssr_id))
        write(6,'(a,i0)')'ssr_id=',ssr_id
        write(6,'(a,i0)')'nc=',nc
       endif
       open(unit=4,file=airseabinreclist(ssr_id)                     &
                  ,access='direct',recl=540,form='unformatted')
        read(4,rec=nc         )texte80(1)
       close(4)

!     endif                         !--initial case-->

      if(txtcase_=='current') then  !--iterative case-->
                if(par%rank==0)write(6,'(a,a)')                   &
                             'Refreshing interpolation procedure' &
                            ,' of meteo files regarding ALBEDO'
      endif                         !--iterative case-->

      if(par%rank==0)write(6,'(a,a)')' Pour savoir albedo ouverture ' &
                                     ,trim(texte80(1))

      status=nf_open(trim(texte80(1)),nf_nowrite,ncid_)
      if(status/=0)stop ' stop nf_open airseaflux_inquire_albedo'


      call airseaflux_varname(ssr_id) ! get meteovarname

                   status=nf_inq_varid(ncid_,meteovarname(1),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(2),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(3),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(4),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(5),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(6),varid_) !08-11-16
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(7),varid_) !17-09-21
      if(status/=0) then !---->
       write(6,'(8(a,1x))')trim(meteovarname(1)) &
                          ,trim(meteovarname(2)) &
                          ,trim(meteovarname(3)) &
                          ,trim(meteovarname(4)) &
                          ,trim(meteovarname(5)) &
                          ,trim(meteovarname(6)) &
                          ,trim(meteovarname(7)) &
                          ,'not found in the netcdf file' &
                          ,trim(texte80(1))
       stop ' Stop airseaflux_inquire_albedo nf_inq_varid'
      endif              !---->

      texte90=''   ! reset
                   status=nf_get_att_text(ncid_,varid_,'long_name',texte90)
      if(status/=0)status=nf_get_att_text(ncid_,varid_,'description',texte90) !17-09-21
      if(status/=0) &
      stop ' Stop nf_get_att_text airseaflux_inquire_albedo'

      if(par%rank==0)write(6,'(a,a)')' long_name= ',trim(texte90)

! Albedo pris en compte si la chaine de caractere selectionnee est
! presente, aurtrement di si "index" est /=0:
       if(index(texte90,'downward')/=0.or.             &
          index(texte90,'DOWNWARD')/=0.or.             & !17-09-21
          index(texte90,'Downwelling')/=0) then

                ialbedo=1

       else
                ialbedo=0

       endif
       if(index(texte90,'Bilan du rayonnement')/=0)ialbedo=0 !08-11-16

       if(par%rank==0) then
         write(6,*)'ialbedo=',ialbedo
!        stop 'koko'
       endif

       status=nf_close(ncid_)
       if(status/=0)stop ' stop nf_close airseaflux_inquire_albedo'


      end subroutine airseaflux_inquire_albedo

!..............................................................................

      subroutine airseaflux_sphum_or_dewp(txtcase_)
      implicit none
      integer ncid_,varid_
      character(len=*)txtcase_
#ifdef synopsis
       subroutinetitle='airseaflux_sphum_or_dewp'
       subroutinedescription= &
          'Determines if the specific humidity or/and dew point are' &
       //' contained in the meteo netcdf files.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif
! Specific humidity or dew point ?

      if(flag_meteodata=='glorys')return !20-04-14

!     if(txtcase_=='initial') then !--Initial Case-->

! Var id = dp2m_id
       if(dp2m_id==0)stop ' airseaflux_sphum_or_dewp dp2m_id non defini'
       nc=airseafile_nextrec(dp2m_id) ! Get netcdf fil name
       open(unit=4,file=airseabinreclist(dp2m_id)                     &
                  ,access='direct',recl=540,form='unformatted')
       read(4,rec=nc         )texte80(1)
       close(4)

!     endif                        !--Initial Case-->

      if(txtcase_=='current') then !--Iterative Case-->
        if(par%rank==0)write(6,'(a,a)')                    &
                      'Refreshing interpolation procedure' &
                    ,' of meteo files regarding sphum_or_dewp'
      endif                        !--Iterative Case-->

      if(par%rank==0)write(6,'(a,a)')'airseaflux_sphum_or_dewp open ',trim(texte80(1))

      status=nf_open(trim(texte80(1)),nf_nowrite,ncid_)
      if(status/=0)stop ' stop nf_open airseaflux_sphum_or_dewp'

      call airseaflux_varname(dp2m_id) ! Get meteovarname

                   status=nf_inq_varid(ncid_,meteovarname(1),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(2),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(3),varid_) !28-04-16
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(4),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(5),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(6),varid_) !17-09-21
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(7),varid_) !13-09-22
      if(status/=0) then
       write(6,'(8a)')'meteovarname=',meteovarname(:)
       stop ' Stop meteovarname airseaflux_sphum_or_dewp not found'
      endif

      texte90=''   ! reset
                   status=nf_get_att_text(ncid_,varid_,'long_name',texte90)
      if(status/=0)status=nf_get_att_text(ncid_,varid_,'description',texte90)
      if(status/=0)stop &
      ' Stop nf_get_att_text airseaflux_sphum_or_dewp'

!     if(index(texte90,'dewpoint')/=0) then

!          flag_humspec_dewpoint='dewpoint'

!     else

!          flag_humspec_dewpoint='humspeci'

!     endif

      flag_humspec_dewpoint='humspeci'            ! Par defaut

      if(index(texte90,'dewpoint')/=0)          &
      flag_humspec_dewpoint='dewpoint'

      if(index(texte90,'relative humidity')/=0) &
      flag_humspec_dewpoint='humrelat'
      if(index(texte90,'Humidite relative')/=0) &
      flag_humspec_dewpoint='humrelat'

      if(index(texte90,'QV at 2 M')/=0)flag_humspec_dewpoint='humspeci' !17-09-21 

      if(par%rank==0)write(6,'(a,a)')'flag_humspec_dewpoint=',flag_humspec_dewpoint

      status=nf_close(ncid_)
      if(status/=0)stop ' stop nf_close airseaflux_sphum_or_dewp'

      end subroutine airseaflux_sphum_or_dewp

!..............................................................................
      subroutine airseaflux_ssr24avr
      implicit none
      if(flag_ssr24avr==0)return
#ifdef synopsis
       subroutinetitle='airseaflux_ssr24avr'
       subroutinedescription='Computes a 24h average of the solar flux'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! 24h average of the solar flux
! Shortcoming warning: the average centre is 12h before the current time

! NOTE: at this stage ssr_w(:,:,2) contains the next 'still NOT averaged' ssr field

! Step 1: allocate the averaging array (if not already done)
      if(.not.allocated(ssr24prv_w)) then !------------->
! The dimension depends on the field sampling (1h --> dim=24 , 3h --> dim=8)
        dimssr24prv=0
        if(nint(( airseafile_nextime(ssr_id)                   &
                 -airseafile_prvtime(ssr_id))/3600.)==3)dimssr24prv=8
        if(nint(( airseafile_nextime(ssr_id)                   &
                 -airseafile_prvtime(ssr_id))/3600.)==1)dimssr24prv=24
        if(dimssr24prv==0) &
        stop 'airseaflux_ssr24prv unrecognised periodicity'
!        dimssr24prv=24 ! (15/11/21: comme les fréquences des forcages ssr sont
! hétérogènes tri-horaire pour 2011, horaire pour 2012 et +, on prend 24)

!       lb3=lbound(ssr_w) ; ub3=ubound(ssr_w)
!       allocate(ssr24prv_w(lb3(1):ub3(1),lb3(2):ub3(2),dimssr24prv))
        allocate(ssr24prv_w(0:imax+1,0:jmax+1,dimssr24prv)) ; ssr24prv_w=0. !04-10-14

      endif                               !------------->

! Move forward the 24hours previous 'not averaged' fields:
      do k=1,dimssr24prv-1
       do j=0,jmax+1 ; do i=0,imax+1
        ssr24prv_w(i,j,k)=ssr24prv_w(i,j,k+1)
       enddo ; enddo
      enddo
      do j=0,jmax+1 ; do i=0,imax+1 !DeplacEe le !11-03-16 
       ssr24prv_w(i,j,dimssr24prv)=ssr_w(i,j,2)
      enddo ; enddo

! Compute the 24h average stored into ssr_w(2):
! Note that at this stage ssr_w(2) already contains the most recent
! 'not averaged' field ssr24prv_w(:,:,dimssr24prv) consequently excluded from the
! following sum:
      do k=1,dimssr24prv-1
       do j=0,jmax+1 ; do i=0,imax+1
        ssr_w(i,j,2)=ssr_w(i,j,2)+ssr24prv_w(i,j,k)
       enddo ; enddo
      enddo
      x1=1./real(dimssr24prv)
      do j=0,jmax+1 ; do i=0,imax+1
       ssr_w(i,j,2)=ssr_w(i,j,2)*x1
      enddo ; enddo

! At this stage ssr_w(2) contains the 24h averaged "next" solar flux

      end subroutine airseaflux_ssr24avr
!..............................................................................
#ifdef checkmpi
      subroutine airseaflux_checkmpi_teta2
      use module_principal
      use module_parallele                                        !#MPI !16-09-09
      implicit none
      integer loop_
#ifdef synopsis
       subroutinetitle='airseaflux_checkmpi_teta2'
       subroutinedescription= &
       'Checks the mpi continuity of the teta2_t field'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Pour verifier la conservation de la parallelisation:
      do j=1,jmax
      do i=1,imax
       xy_t(i,j,1)=teta2_t(i,j,2)
      enddo
      enddo
      call obc_h_xyt1

      ksecu=0
      j1=1 ; j2=jmax

      if(par%tvoisin(sud)/=mpi_proc_null) then !sssssssssss>
      do i=1,imax
      if(mask_t(i,j1,kmax)==1) then
       if(xy_t(i,j1,1)/=teta2_t(i,j1,2))then
        ksecu=1
        write(10+par%rank,*)'mpi non conservé en par%rank i j =',par%rank,i,j1
        write(10+par%rank,*)'voisin sud',par%tvoisin(sud)
        k0=par%tvoisin(sud) ; j0=par%gtjmax(k0,2)-par%gtjmax(k0,1) ! j0=jmax du voisin
        write(10+par%rank,*)'coordonnées ',i,j0-1   ! jmax-1
        write(10+par%rank,*) &
         'teta2',xy_t(i,j1,1),teta2_t(i,j1,2),xy_t(i,j1,1)-teta2_t(i,j1,2)
        write(10+par%rank,*)'mask_t & h_w:',mask_t(i,j1,kmax),h_w(i,j1)
        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      endif
      enddo
      endif                        !sssssssssss>


      if(par%tvoisin(nord)/=mpi_proc_null) then !nnnnnnnnnnn>
      do i=1,imax
      if(mask_t(i,j2,kmax)==1) then
       if(xy_t(i,j2,1)/=teta2_t(i,j2,2))then
        ksecu=1
        write(10+par%rank,*)'mpi non conservé en par%rank i j =' &
         ,par%rank,i,j2
        write(10+par%rank,*)'          voisin en par%rank i j =' &
         ,par%tvoisin(nord),i,2
        write(10+par%rank,*) &
         'teta2',xy_t(i,j2,1),teta2_t(i,j2,2),xy_t(i,j2,1)-teta2_t(i,j2,2)
        write(10+par%rank,*)'mask_t & h_w:',mask_t(i,j2,kmax),h_w(i,j2)
        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      endif
      enddo
      endif                         !nnnnnnnnnnn>


      i1=1 ; i2=imax

      if(par%tvoisin(ouest)/=mpi_proc_null) then !ooooooooooo>
      do j=1,jmax
      if(mask_t(i1,j,kmax)==1) then
       if(xy_t(i1,j,1)/=teta2_t(i1,j,2)) then
        ksecu=1
        write(10+par%rank,*)'mpi non conservé en par%rank i j =' &
          ,par%rank,i1,j
        write(10+par%rank,*)'voisin ouest=',par%tvoisin(ouest)
        k0=par%tvoisin(ouest) ; i0=par%gtimax(k0,2)-par%gtimax(k0,1) ! i0=imax du voisin
        write(10+par%rank,*)'coordonnées ',i0-1,j ! imax-1,j
        write(10+par%rank,*) &
          'teta2',xy_t(i1,j,1),teta2_t(i1,j,2),xy_t(i1,j,1)-teta2_t(i1,j,2)
        write(10+par%rank,*)'mask_t & h_w:',mask_t(i1,j,kmax),h_w(i1,j)
        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      endif
      enddo
      endif                          !ooooooooooo>

      if(par%tvoisin(est)/=mpi_proc_null) then   !eeeeeeeeeee>
      do j=1,jmax
      if(mask_t(i2,j,kmax)==1) then
       if(xy_t(i2,j,1)/=teta2_t(i2,j,2)) then
        ksecu=1
        write(10+par%rank,*)'mpi non conservé en par%rank i j =' &
          ,par%rank,i2,j
        write(10+par%rank,*)'          voisin en par%rank i j =' &
          ,par%tvoisin(est),2,j
        write(10+par%rank,*) &
          'teta2',xy_t(i2,j,1),teta2_t(i2,j,2),xy_t(i2,j,1)-teta2_t(i2,j,2)
        write(10+par%rank,*)'mask_t & h_w:',mask_t(i2,j,kmax),h_w(i2,j)
        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      endif
      enddo
      endif                          !eeeeeeeeeee>

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
      if(ksecu==1)then
      write(6,*)'Erreurs mpi notifiees dans fichiers fort locaux'
        stop 'STOP dans airseaflux_checkmpi_teta2'
      endif

      end subroutine airseaflux_checkmpi_teta2
#endif
!..............................................................................
#ifdef checkmpi
      subroutine airseaflux_checkmpi_q2
      use module_principal
      use module_parallele                                        !#MPI !16-09-09
      implicit none
      integer loop_
#ifdef synopsis
       subroutinetitle='airseaflux_checkmpi_q2'
       subroutinedescription= &
       'Checks the mpi continuity of the q2_t field'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Pour verifier la conservation de la parallelisation:
      do j=1,jmax
      do i=1,imax
       xy_t(i,j,1)=q2_t(i,j,2)
      enddo
      enddo
      call obc_h_xyt1

      ksecu=0
      j1=1 ; j2=jmax

      if(par%tvoisin(sud)/=mpi_proc_null) then !sssssssssss>
      do i=1,imax
      if(mask_t(i,j1,kmax)==1) then
       if(xy_t(i,j1,1)/=q2_t(i,j1,2))then
        ksecu=1
        write(10+par%rank,*)'mpi non conservé en par%rank i j =',par%rank,i,j1
        write(10+par%rank,*)'voisin sud',par%tvoisin(sud)
        k0=par%tvoisin(sud) ; j0=par%gtjmax(k0,2)-par%gtjmax(k0,1) ! j0=jmax du voisin
        write(10+par%rank,*)'coordonnées ',i,j0-1   ! jmax-1
        write(10+par%rank,*) &
         'q2',xy_t(i,j1,1),q2_t(i,j1,2),xy_t(i,j1,1)-q2_t(i,j1,2)
        write(10+par%rank,*)'mask_t & h_w:',mask_t(i,j1,kmax),h_w(i,j1)
        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      endif
      enddo
      endif                        !sssssssssss>


      if(par%tvoisin(nord)/=mpi_proc_null) then !nnnnnnnnnnn>
      do i=1,imax
      if(mask_t(i,j2,kmax)==1) then
       if(xy_t(i,j2,1)/=q2_t(i,j2,2))then
        ksecu=1
        write(10+par%rank,*)'mpi non conservé en par%rank i j =' &
         ,par%rank,i,j2
        write(10+par%rank,*)'          voisin en par%rank i j =' &
         ,par%tvoisin(nord),i,2
        write(10+par%rank,*) &
         'q2',xy_t(i,j2,1),q2_t(i,j2,2),xy_t(i,j2,1)-q2_t(i,j2,2)
        write(10+par%rank,*)'mask_t & h_w:',mask_t(i,j2,kmax),h_w(i,j2)
        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      endif
      enddo
      endif                         !nnnnnnnnnnn>


      i1=1 ; i2=imax

      if(par%tvoisin(ouest)/=mpi_proc_null) then !ooooooooooo>
      do j=1,jmax
      if(mask_t(i1,j,kmax)==1) then
       if(xy_t(i1,j,1)/=q2_t(i1,j,2)) then
        ksecu=1
        write(10+par%rank,*)'mpi non conservé en par%rank i j =' &
          ,par%rank,i1,j
        write(10+par%rank,*)'voisin ouest=',par%tvoisin(ouest)
        k0=par%tvoisin(ouest) ; i0=par%gtimax(k0,2)-par%gtimax(k0,1) ! i0=imax du voisin
        write(10+par%rank,*)'coordonnées ',i0-1,j ! imax-1,j
        write(10+par%rank,*) &
          'q2',xy_t(i1,j,1),q2_t(i1,j,2),xy_t(i1,j,1)-q2_t(i1,j,2)
        write(10+par%rank,*)'mask_t & h_w:',mask_t(i1,j,kmax),h_w(i1,j)
        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      endif
      enddo
      endif                          !ooooooooooo>

      if(par%tvoisin(est)/=mpi_proc_null) then   !eeeeeeeeeee>
      do j=1,jmax
      if(mask_t(i2,j,kmax)==1) then
       if(xy_t(i2,j,1)/=q2_t(i2,j,2)) then
        ksecu=1
        write(10+par%rank,*)'mpi non conservé en par%rank i j =' &
          ,par%rank,i2,j
        write(10+par%rank,*)'          voisin en par%rank i j =' &
          ,par%tvoisin(est),2,j
        write(10+par%rank,*) &
          'q2',xy_t(i2,j,1),q2_t(i2,j,2),xy_t(i2,j,1)-q2_t(i2,j,2)
        write(10+par%rank,*)'mask_t & h_w:',mask_t(i2,j,kmax),h_w(i2,j)
        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      endif
      enddo
      endif                          !eeeeeeeeeee>

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
      if(ksecu==1)then
      write(6,*)'Erreurs mpi notifiees dans fichiers fort locaux'
        stop 'STOP dans airseaflux_checkmpi_q2'
      endif

      end subroutine airseaflux_checkmpi_q2
#endif
!..............................................................................
#ifdef checkmpi
      subroutine airseaflux_checkmpi_uwind
      use module_principal
      use module_parallele                                        !#MPI !16-09-09
      implicit none
      integer loop_
#ifdef synopsis
       subroutinetitle='airseaflux_checkmpi_uwind'
       subroutinedescription= &
       'Checks the mpi continuity of the uwind_t field'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Pour verifier la conservation de la parallelisation:
      do j=1,jmax
      do i=1,imax
       xy_t(i,j,1)=uwind_t(i,j,2)
      enddo
      enddo
      call obc_h_xyt1

      ksecu=0
      j1=1 ; j2=jmax

      if(par%tvoisin(sud)/=mpi_proc_null) then !sssssssssss>
      do i=1,imax
      if(mask_t(i,j1,kmax)==1) then
       if(xy_t(i,j1,1)/=uwind_t(i,j1,2))then
        ksecu=1
        write(10+par%rank,*)'mpi non conservé en par%rank i j =',par%rank,i,j1
        write(10+par%rank,*)'voisin sud',par%tvoisin(sud)
        k0=par%tvoisin(sud) ; j0=par%gtjmax(k0,2)-par%gtjmax(k0,1) ! j0=jmax du voisin
        write(10+par%rank,*)'coordonnées ',i,j0-1   ! jmax-1
        write(10+par%rank,*) &
         'uwind',xy_t(i,j1,1),uwind_t(i,j1,2),xy_t(i,j1,1)-uwind_t(i,j1,2)
        write(10+par%rank,*)'mask_t & h_w:',mask_t(i,j1,kmax),h_w(i,j1)
        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      endif
      enddo
      endif                        !sssssssssss>


      if(par%tvoisin(nord)/=mpi_proc_null) then !nnnnnnnnnnn>
      do i=1,imax
      if(mask_t(i,j2,kmax)==1) then
       if(xy_t(i,j2,1)/=uwind_t(i,j2,2))then
        ksecu=1
        write(10+par%rank,*)'mpi non conservé en par%rank i j =' &
         ,par%rank,i,j2
        write(10+par%rank,*)'          voisin en par%rank i j =' &
         ,par%tvoisin(nord),i,2
        write(10+par%rank,*) &
         'uwind',xy_t(i,j2,1),uwind_t(i,j2,2),xy_t(i,j2,1)-uwind_t(i,j2,2)
        write(10+par%rank,*)'mask_t & h_w:',mask_t(i,j2,kmax),h_w(i,j2)
        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      endif
      enddo
      endif                         !nnnnnnnnnnn>


      i1=1 ; i2=imax

      if(par%tvoisin(ouest)/=mpi_proc_null) then !ooooooooooo>
      do j=1,jmax
      if(mask_t(i1,j,kmax)==1) then
       if(xy_t(i1,j,1)/=uwind_t(i1,j,2)) then
        ksecu=1
        write(10+par%rank,*)'mpi non conservé en par%rank i j =' &
          ,par%rank,i1,j
        write(10+par%rank,*)'voisin ouest=',par%tvoisin(ouest)
        k0=par%tvoisin(ouest) ; i0=par%gtimax(k0,2)-par%gtimax(k0,1) ! i0=imax du voisin
        write(10+par%rank,*)'coordonnées ',i0-1,j ! imax-1,j
        write(10+par%rank,*) &
          'uwind',xy_t(i1,j,1),uwind_t(i1,j,2),xy_t(i1,j,1)-uwind_t(i1,j,2)
        write(10+par%rank,*)'mask_t & h_w:',mask_t(i1,j,kmax),h_w(i1,j)
        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      endif
      enddo
      endif                          !ooooooooooo>

      if(par%tvoisin(est)/=mpi_proc_null) then   !eeeeeeeeeee>
      do j=1,jmax
      if(mask_t(i2,j,kmax)==1) then
       if(xy_t(i2,j,1)/=uwind_t(i2,j,2)) then
        ksecu=1
        write(10+par%rank,*)'mpi non conservé en par%rank i j =' &
          ,par%rank,i2,j
        write(10+par%rank,*)'          voisin en par%rank i j =' &
          ,par%tvoisin(est),2,j
        write(10+par%rank,*) &
          'uwind',xy_t(i2,j,1),uwind_t(i2,j,2),xy_t(i2,j,1)-uwind_t(i2,j,2)
        write(10+par%rank,*)'mask_t & h_w:',mask_t(i2,j,kmax),h_w(i2,j)
        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      endif
      enddo
      endif                          !eeeeeeeeeee>

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
      if(ksecu==1)then
      write(6,*)'Erreurs mpi notifiees dans fichiers fort locaux'
        stop 'STOP dans airseaflux_checkmpi_uwind'
      endif

      end subroutine airseaflux_checkmpi_uwind
#endif

!..........................................................................

      subroutine airseaflux_ij2meteo_mpi(txt_) !23-04-15
      use module_principal
      use module_parallele
      implicit none
      integer loop_,varid1_,varid2_
      character*2 txt_
#ifdef synopsis
       subroutinetitle='airseaflux_ij2meteo_mpi'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

#ifdef parallele
       write(texte30,'(a,a)')'ij2meteo_i_',trim(txt_)
       call get_type_echange(txt_,trim(texte30),ij2meteo_i        &
                                        ,lbound(ij2meteo_i)       &
                                        ,ubound(ij2meteo_i),varid1_)


       write(texte30,'(a,a)')'ij2meteo_j_',trim(txt_)
       call get_type_echange(txt_,trim(texte30),ij2meteo_j        &
                                        ,lbound(ij2meteo_j)       &
                                        ,ubound(ij2meteo_j),varid2_)

      ! Echanges
      do loop_=1,subcycle_exchange
        call echange_voisin(ij2meteo_i,varid1_,mpi_neighbor_list(loop_))
        call echange_voisin(ij2meteo_j,varid2_,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine airseaflux_ij2meteo_mpi

!.................................................................

      subroutine airseaflux_2dnoiseremover(var_id_) !26-06-15
      implicit none
      integer var_id_
#ifdef synopsis
       subroutinetitle='airseaflux_2dnoiseremover'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! 5 points filter:
      allocate(meteo_var2(meteo_imax,meteo_jmax))
! Algo 1:
!     do j=2,meteo_jmax-1
!     do i=2,meteo_imax-1
!       meteo_var2(i,j)=0.5*    meteo_var(i  ,j  )    &
!                      +0.125*( meteo_var(i+1,j  )    &
!                              +meteo_var(i-1,j  )    &
!                              +meteo_var(i  ,j-1)    &
!                              +meteo_var(i  ,j+1))
!     enddo
!     enddo
! Algo 2:
      do j=1,meteo_jmax
      do i=2,meteo_imax-1
        meteo_var2(i,j)=0.5*   meteo_var(i  ,j  )    &
                       +0.25*( meteo_var(i+1,j  )    &
                              +meteo_var(i-1,j  ) )
      enddo
      enddo
      do j=1,meteo_jmax
      do i=2,meteo_imax-1
       meteo_var(i,j)=meteo_var2(i,j)
      enddo
      enddo

      do j=2,meteo_jmax-1
      do i=1,meteo_imax
        meteo_var2(i,j)=0.5*   meteo_var(i  ,j  )    &
                       +0.25*( meteo_var(i  ,j-1)    &
                              +meteo_var(i  ,j+1))
      enddo
      enddo
      do j=2,meteo_jmax-1
      do i=1,meteo_imax
       meteo_var(i,j)=meteo_var2(i,j)
      enddo
      enddo


      deallocate(meteo_var2)

      end subroutine airseaflux_2dnoiseremover

!.................................................................

      subroutine airseaflux_oldlist(loop1_)
      implicit none
      integer unit_,loop1_,flag_
      double precision time_

! Lire une liste binaire deja existante replacee intentionnellement
! dans le repertoire tmp) et se positionner entre les 2 records encadrant
! le temps initial elapsedtime_now

       flag_=0
       unit_=s_unit(7)
       open(unit=unit_,file=airseabinreclist(loop1_) &
               ,access='direct'                      &
               ,recl=540                             &
               ,form='unformatted'                   &
               ,status='old'                         &
               ,iostat=k0)
       if(k0/=0)stop 'Err 5998'

        i0=0
        nc=0

         do while (i0==0)

          nc=nc+1

         if(nc/=1) then !>>>>
          read(unit_,rec=nc,iostat=i0)       &
                            texte80(1)       & ! Nom fichier netcdf
                           ,k10              & ! Numero d'echeance dans le fichier netdf
                           ,time_
         else           !>>>>
          read(unit_,rec=nc,iostat=i0)     &
                          texte80(1)       & ! Nom fichier netcdf
                         ,k10              & ! Numero d'echeance dans le fichier netdf
                         ,time_            & ! temps en seconde dans le repere de S26
                         ,x1               & ! temps en seconde depuis le debut du run ecmwf
                         ,i                & ! Cumul ou pas cumul ?
                         ,x2               & ! facteur d'echelle aditionnel
                         ,i1,i2,i3,i4,i5,i6  ! date repere nc=1

          if(i1/=datesim(1,1).or. &
             i2/=datesim(2,1).or. &
             i3/=datesim(3,1).or. &
             i4/=datesim(4,1).or. &
             i5/=datesim(5,1).or. &
             i6/=datesim(6,1)) then !bugbug>
             write(6,*)'Impossible de repartir de la liste binaire' &
                      ,'car date repere nc=1:'                      &
                      ,i1,i2,i3,i4,i5,i6                            &
                      ,' ne correspond pas A date'                  &
                      ,' depart indiquEe dans notebook_time:'       &
                      ,datesim(:,1)
             stop 'Err 6057'

          endif                     !bugbug>

         endif          !>>>>

         if(i0==0) then !pmx>
          if(time_<=elapsedtime_now) then !>>>> 
           airseafile_nextrec(loop1_)=nc  
           airseafile_nextime(loop1_)=time_
           flag_=1
          endif                           !>>>>
         endif          !pmx>


         enddo

      close(unit_)



      if(flag_==0) STOP 'Pas de fichiers METEO autour du temps initial'

      if(par%rank==0) then
        write(6,*)'subroutine airseaflux_oldlist:'
        write(6,*)'   loop1_=',loop1_
        write(6,*)'   airseafile_nextime(loop1_)',airseafile_nextime(loop1_)
        write(6,*)'   airseafile_nextrec(loop1_)',airseafile_nextrec(loop1_)
      endif

      end subroutine airseaflux_oldlist

!.................................................................

      subroutine airseaflux_checkfileconsistecy
      implicit none
      integer ncid_

! Cette routine sert A detecter des inhomogeneites de fichiers dans
! la liste. Pour le moment on verifie simplement que les dimensions
! ne changent pas en cours de liste
!     if(par%rank==0)write(6,'(a,a)') &
!                   'airseaflux_checkfileconsistecy Open ' &
!                  ,trim(texte80(1))

      status=nf_open(trim(texte80(1)),nf_nowrite,ncid_)
      if(status/=0)stop 'Err 6223 airseaflux_checkfileconsistecy'

      k10=0 ! OK si k10 reste 0

                   status=nf_inq_dimid(ncid_,'lon',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'longitude',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'g0_lon_1',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'x',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'west_east',dim_x_id) !17-09-21
      if(status/=0)stop 'Err dim_x_id variable meteo longitude'

                   status=nf_inq_dimid(ncid_,'g0_lat_0',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'latitude',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'lat',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'y',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'south_north',dim_y_id) !17-09-21
      if(status/=0)stop 'Err dim_y_id variable meteo latitude'

!................................
! Verifier dimension 1:

      status=nf_inq_dimlen(ncid_,dim_x_id,i0)
      if(status/=0) &
      stop 'Err 6156 nf_inq_dimlen airseaflux_checkfileconsistecy'      

! Cas particulier du premier fichier qui donne sa valeur initiale A meteo_imax_full
      if(meteo_imax_full==0)meteo_imax_full=i0

      if(i0/=meteo_imax_full) then !pmxpmx>
        if(par%rank==0) &
        write(6,'(3a,2(i0,1x))')'Inhomogeneous dimension in meteo file ' &
        ,trim(texte80(1)),' different from previous dim1 val:' &
        ,meteo_imax_full,i0
        k10=1 ! Disconstinuite dans la liste
!      stop 'Err 6159 airseaflux_checkfileconsistecy'
      endif                        !pmxpmx>

!................................
! Verifier dimension 2:

      status=nf_inq_dimlen(ncid_,dim_y_id,i0)
      if(status/=0) &
      stop 'Err 6157 nf_inq_dimlen airseaflux_checkfileconsistecy'      

! Cas particulier du premier fichier qui donne sa valeur initiale A meteo_jmax_full
      if(meteo_jmax_full==0)meteo_jmax_full=i0

      if(i0/=meteo_jmax_full) then !pmxpmx>
        if(par%rank==0) &
       write(6,'(3a,2(i0,1x))')'Inhomogeneous dimension in meteo file ' &
       ,trim(texte80(1)),' different from previous dim2 val:' &
       ,meteo_jmax_full,i0
        k10=1 ! Disconstinuite dans la liste
!     stop 'Err 6160 airseaflux_checkfileconsistecy'
      endif                        !pmxpmx>

      status=nf_close(ncid_) 
      if(status/=0)stop 'Err 6278 airseaflux_checkfileconsistecy'

! Si discontinuite dans la liste alors mise a jour de la chaine
! d'interpolation:
      if(k10==1)call airseaflux_refresh_interp_proc !05-05-16

      end subroutine airseaflux_checkfileconsistecy

!.................................................................

      subroutine airseaflux_refresh_interp_proc !05-05-16
      implicit none

! On passe par cette subroutine car on vient de detecter un changement
! dans l'entete du fichier meteo qui entraine la mise a jour de la
! chaine d'interpolation
      if(par%rank==0)write(6,'(a,a)')                    &
       ' Discontinuity detected in the meteo file list ' &
      ,' - Updating interpolation procedure '

! Devra t'on calculer l'albedo ?
      call airseaflux_inquire_albedo('current') ! ialbedo

! Devra t'on deduire l'humidite specifique du dew point?
      call airseaflux_sphum_or_dewp('current')

! Le flux infra rouge est t'il le bilan atmosphere/ocean (flag_net_ir=1) ou seulement le
! flux atmospherique (flag_net_ir=0):
      call airseaflux_inquire_longwave('current') !28-11-16

! Reduire l'extraction des donnees a une zone englobant notre domaine
! numerique et definir les bouchetrous et les relations de passage d'une grille a
! l'autre (ecmwf vs s26 etc...)
      call allocate_forcages(2,2,0,0,0) ! deallocate des tableaux meteo
                                        ! avant remise a jour de la
                                        ! chaine d'interpolation
      call airseaflux_extract_zone('current')

      end subroutine airseaflux_refresh_interp_proc

!.................................................................

      subroutine airseaflux_inquire_longwave(txtcase_) !28-11-16
      implicit none
      integer ncid_,varid_,ir_or_netir_
      character(len=*)txtcase_
      
#ifdef synopsis
       subroutinetitle='airseaflux_inquire_longwave'
       subroutinedescription= &
         'Determines if albebo should be applied on solar flux'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(flag_meteodata=='glorys')return !19-01-17

! Include albedo ?
! Var id = ir_id
!     if(txtcase_=='initial') then   !--initial case-->

       ir_or_netir_=0
       if(ir_id/=0)   ir_or_netir_=ir_id    !19-01-17
       if(netir_id/=0)ir_or_netir_=netir_id !19-01-17
       if(par%rank==0.and.ir_or_netir_==0) &
       stop 'Err ir_id=netir_id=0 airseaflux_inquire_longwave'

       nc=airseafile_nextrec(ir_or_netir_) ! Get netcdf fil name
       if(par%rank==0) then
        write(6,'(a)')'Subroutine airseaflux_inquire_longwave:'
        write(6,'(a,a)')'Open file ',trim(airseabinreclist(ir_or_netir_))
        write(6,'(a,i0)')'ir_or_netir_=',ir_or_netir_
        write(6,'(a,i0)')'nc=',nc
       endif
       open(unit=4,file=airseabinreclist(ir_or_netir_)           &
                  ,access='direct',recl=540,form='unformatted')
        read(4,rec=nc         )texte80(1)
       close(4)

!     endif                         !--initial case-->

      if(txtcase_=='current') then  !--iterative case-->
                if(par%rank==0)write(6,'(a,a)')                   &
                             'Refreshing interpolation procedure' &
                            ,' of meteo files regarding longwave flux'
      endif                         !--iterative case-->

      if(par%rank==0)write(6,'(a,a)')' Pour savoir net ou descendant: ' &
                                     ,trim(texte80(1))

      status=nf_open(trim(texte80(1)),nf_nowrite,ncid_)
      if(status/=0)stop ' stop nf_open airseaflux_inquire_longwave'


      call airseaflux_varname(ir_or_netir_) ! get meteovarname

                   status=nf_inq_varid(ncid_,meteovarname(1),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(2),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(3),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(4),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(5),varid_)
      if(status/=0)status=nf_inq_varid(ncid_,meteovarname(6),varid_) !08-11-16
      if(status/=0) then !---->
       write(6,'(7(a,1x))')trim(meteovarname(1)) &
                          ,trim(meteovarname(2)) &
                          ,trim(meteovarname(3)) &
                          ,trim(meteovarname(4)) &
                          ,trim(meteovarname(5)) &
                          ,trim(meteovarname(6)) &
                          ,'not found in the netcdf file' &
                          ,trim(texte80(1))
       stop ' Stop airseaflux_inquire_longwave nf_inq_varid'
      endif              !---->

      texte90=''   ! reset
                   status=nf_get_att_text(ncid_,varid_,'long_name',texte90)
      if(status/=0)status=nf_get_att_text(ncid_,varid_,'description',texte90) !17-09-21
      if(status/=0) &
      stop ' Stop nf_get_att_text airseaflux_inquire_longwave'

      if(par%rank==0)write(6,'(a,a)')' long_name= ',trim(texte90)

! Albedo pris en compte si la chaine de caractere selectionnee est
! presente, aurtrement di si "index" est /=0:
       if(index(texte90,'Bilan du rayonnement')/=0) then
         flag_net_ir=1 ! Le flux IR est net (IR atmos + IR ocean) donc ne pas calculer IR montant avec SST
       else
         flag_net_ir=0 ! Le flux IR est seulement la composante atmospherique donc calculer IR montant avec SST
       endif
       if(index(texte90,'DOWNWARD')/=0)flag_net_ir=0 !flux descendant seulement !17-09-21

       if(par%rank==0) then
!        write(6,'(a,a)')'texte90',trim(texte90)
         write(6,*)'flag_net_ir=',flag_net_ir
!        stop 'koko' !15-02-17
       endif

       status=nf_close(ncid_)
       if(status/=0)stop ' stop nf_close airseaflux_inquire_longwave'


      end subroutine airseaflux_inquire_longwave

!..............................................................................

      subroutine airseaflux_check_lsm
      implicit none
      double precision airseafile_nextime_,time1_
      integer vstart1_,nc_str_,nc_end_,loop_var_,loop_rank_
      integer,dimension(:),allocatable :: tmp_refresh

! flag_lsm=1 = land-sea mask is not needed, 0 otherwise (run stops if lsm not found)
      if(flag_lsm==1)return

! Les domaines mpi se partagent le travail, de nc=nc_str_ A nc=nc_end_:
      nc_str_=  floor(real(nc_meteo_max)*real(par%rank  )/real(nbdom))
      nc_end_=ceiling(real(nc_meteo_max)*real(par%rank+1)/real(nbdom))
      nc_str_=min(max(nc_str_,1),nc_meteo_max)
      nc_end_=min(max(nc_end_,1),nc_meteo_max)
      if(par%rank==0)nc_str_=1
      if(par%rank==nbdom-1)nc_end_=nc_meteo_max

! tmp_refresh est un tableau representat flag_refresh_interp pour chaque valeur de nc
      allocate(tmp_refresh(nc_str_:nc_end_)) ; tmp_refresh=-99

!     write(10+par%rank,*)'nc_meteo_max',nc_meteo_max
!     write(10+par%rank,*)real(nc_meteo_max)*real(par%rank)/real(nbdom)
!     write(10+par%rank,*)real(nc_meteo_max)*real(par%rank+1)/real(nbdom)
!     write(10+par%rank,*)'nc_str_,nc_end_=',nc_str_,nc_end_

      open(unit=4,file=airseabinreclist(1)                     &
                 ,access='direct',recl=540,form='unformatted')

      do nc=nc_str_,nc_end_

      tmp_refresh(nc)=0 ! reset, par defaut 0

      read(4,rec=nc         )texte80(1),vstart1_                   &
                                       ,airseafile_nextime_        & !  Ne pas perdre la valeur de airseafile_nextime(:) en cours
                                       ,time1_                     &
                                       ,flag_cumul                 &
                                       ,scalefct
                                 ! Le 7eme element, non lu, est soit datesim(1:6,1) si nc=1 soit tmp_refresh(nc) sinon.


      status=nf_open(trim(texte80(1)),nf_nowrite,ncid1)
      if(status/=0)stop 'echec fichier routine fichier_bouchetrou_meteo'


! Lire les dimensions puis verifier si elles changent d'un fichier a l'autre
                   status=nf_inq_dimid(ncid1,'lon',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'longitude',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'g0_lon_1',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'x',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'west_east',dim_x_id) !17-09-21
      if(status/=0)stop 'Err 154 dim_x_id variable meteo longitude'

                   status=nf_inq_dimid(ncid1,'g0_lat_0',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'latitude',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'lat',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'y',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'south_north',dim_y_id) !17-09-21
      if(status/=0)stop 'Err 155 dim_y_id variable meteo latitude'

!................................
! Verifier dimension 1:
      i1=i2
      status=nf_inq_dimlen(ncid1,dim_x_id,i2)
      if(status/=0)stop 'Err 0156 nf_inq_dimlen i2'      
      if(nc/=nc_str_.and.i1/=i2) then !>>>
       tmp_refresh(nc)=1    ! grille differente, refaire l'interpolation
      endif                           !>>>

!................................
! Verifier dimension 2:
      j1=j2
      status=nf_inq_dimlen(ncid1,dim_y_id,j2)
      if(status/=0)stop 'Err 0157 nf_inq_dimlen j2'      
      if(nc/=nc_str_.and.j1/=j2) then !>>>
       tmp_refresh(nc)=1 ! grille differente, refaire l'interpolation
      endif                           !>>>

      if(tmp_refresh(nc)==1) then !pmx>
! Si tmp_refresh(nc)=1, c'est que les tests au dessus ont trouvE que les dimensions ont changees.
! Il faut donc desallouer les tableaux:
       if(allocated(meteo_var)) deallocate(meteo_var)
       if(allocated(meteo_var2))deallocate(meteo_var2)
      endif                       !pmx>

! Meme si les dimensions sont les memes, maintenant on verifie que le LSM ne change
! pas d'un fichier a l'autre
      if(.not.allocated(meteo_var)) then
               allocate(meteo_var (i2,j2)) ; meteo_var=1. 
      endif
      if(.not.allocated(meteo_var2)) then
               allocate(meteo_var2(i2,j2)) ; meteo_var2=1.  
      endif
      if(nc>nc_str_)meteo_var2(1:i2,1:j2)=meteo_var(1:i2,1:j2)

      varstart(4)=1              ; varcount(4)=1          ! time
      varstart(3)=1              ; varcount(3)=1          ! height_4
      varstart(2)=1              ; varcount(2)=j2 ! lat
      varstart(1)=1              ; varcount(1)=i2 ! lon
                     status=nf_inq_varid(ncid1,'land_sea_mask',var_id)
      if(status.ne.0)status=nf_inq_varid(ncid1,'lsm',var_id)
      if(status.ne.0)status=nf_inq_varid(ncid1,'LSM',var_id)
      if(status.ne.0)status=nf_inq_varid(ncid1,'LSM_GDS0_SFC',var_id)
      if(status.ne.0) then !pmxpmx>
            write(6,'(a,a)')'Cant not find land sea mask in file ',trim(texte80(1)) !28-04-16
            stop 'echec id land_sea_mask'
      endif                !pmxpmx>

      status=nf_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(status.ne.0)stop 'echec nf_inq_var land_sea_mask'

      if(var_type==nf_real) then !rrrrrr>
       status=nf_get_vara_real(ncid1,var_id               &
                                    ,varstart(1:var_dims) &
                                    ,varcount(1:var_dims) &
                                    ,meteo_var(1:i2,1:j2))
       if(status.ne.0)stop 'echec lecture land_sea_mask'
      endif                      !rrrrrr>

      if(var_type==nf_short) then !iiiii>
       status=nf_get_vara_int(ncid1,var_id               &
                                   ,varstart(1:var_dims) &
                                   ,varcount(1:var_dims) &
                 ,meteo_short(1:i2,1:j2))
       if(status.ne.0)stop 'echec lecture land_sea_mask'
       status=nf_get_att_real(ncid1,var_id,'scale_factor',var_scalefactor)
       if(status/=0) &
       stop 'error get scale_factor 2976'
       status=nf_get_att_real(ncid1,var_id,'add_offset',var_addoffset)
       if(status/=0)stop 'error get add_offset 2976'
       do j=1,j2 ; do i=1,i2
        meteo_var(i,j)=meteo_short(i,j)*var_scalefactor+var_addoffset
       enddo ; enddo
      endif                       !iiiii>

!     if(meteolandconvention==0)meteo_var(:,:)=1-meteo_var(:,:) 
! Convention: 0=mer, 1=terre
! Le LSM est reel. Valeurs entre 0 et 1 existent. On considere "terre" (trous à boucher) si LSM/=0, autrement dit LSM=0.1 sera considere en terre
! Pour savoir si le masque evolue (au sens oU notre bouchage l'entend) on arrondit donc d'abord les valeurs strictement positives A 1
!     meteo_var=ceiling(meteo_var)

      sum1=0.
      if(nc>nc_str_) then !ooo>
        do j=1,j2 ; do i=1,i2
         sum1=sum1+abs(meteo_var(i,j)-meteo_var2(i,j))
!        if(meteo_var(i,j)/=meteo_var2(i,j)) then
!         write(100+par%rank,*)nc &
!         ,i,j,meteo_var(i,j),meteo_var2(i,j),'trouvemi'
!        endif
        enddo ; enddo
       if(sum1/=0.)tmp_refresh(nc)=1
!      write(200+par%rank,*)nc,tmp_refresh(nc)

      endif               !ooo>

 2022 continue

      status=nf_close(ncid1)
      enddo ! nc=nc_str_,nc_end_
      close(4)

      if(allocated(meteo_var)) deallocate(meteo_var)
      if(allocated(meteo_var2))deallocate(meteo_var2)

! L'un apres l'autre, chaque domaine ecrit son flag_refresh_interp dans les listes binaires
      do loop_rank_=0,nbdom-1 
      if(loop_rank_==par%rank) then !--rank-->

       do loop_var_=1,nairsea
        open(unit=4,file=airseabinreclist(loop_var_)                     &
                  ,access='direct',recl=540,form='unformatted')
         do nc=max(nc_str_+1,2),nc_end_ ! attention en nc=1, le 7 element est datesim(1:6,1)
! note: la boucle sur nc commence en nc_str_+1 car en nc=nc_str_ tmp_refresh(nc) n'existe pas
           read(4,rec=nc)texte80(1)          & !1
                        ,vstart1_            & !2
                        ,airseafile_nextime_ & !3 Ne pas perdre la valeur de airseafile_nextime(:) en cours
                        ,time1_              & !4
                        ,flag_cumul          & !5
                        ,scalefct              !6
          write(4,rec=nc)texte80(1)          & !1
                        ,vstart1_            & !2
                        ,airseafile_nextime_ & !3 Ne pas perdre la valeur de airseafile_nextime(:) en cours
                        ,time1_              & !4
                        ,flag_cumul          & !5
                        ,scalefct            & !6
                        ,tmp_refresh(nc)       !7
         enddo ! nc
        close(4)
       enddo ! loop_var_=1,nairsea

      endif                         !--rank-->
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
      enddo  !loop_rank_=0,nbdom-1
      deallocate(tmp_refresh)


! LES LIGNES EN SUIVANT SERVENT UNIQUEMENT A FAIRE DES TESTS DE VALIDATION 
#ifdef bidon
      if(par%rank==0) then
       do loop_var_=1,nairsea
       write(10+par%rank+loop_var_,'(a)')trim(airseabinreclist(loop_var_))
        open(unit=4,file=airseabinreclist(loop_var_)                     &
                  ,access='direct',recl=540,form='unformatted')
         do nc=1,nc_meteo_max ! attention en nc=1, le 7 element est datesim(1:6,1)
! note: la boucle sur nc commence en nc_str_+1 car en nc=nc_str_ tmp_refresh(nc) n'existe pas
           read(4,rec=nc)texte80(1)          & !1
                        ,vstart1_            & !2
                        ,airseafile_nextime_ & !3 Ne pas perdre la valeur de airseafile_nextime(:) en cours
                        ,time1_              & !4
                        ,flag_cumul          & !5
                        ,scalefct            & !6
                        ,flag_refresh_interp

          write(10+par%rank+loop_var_,'(a,1x,i0,1x,i0)') &
                         trim(texte80(1))    & !1
                        ,vstart1_            & !2
                        ,flag_refresh_interp
         enddo ! nc
        enddo ! loop_var_
       endif
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
#endif
      stop 'beber'
#endif

      end subroutine airseaflux_check_lsm

!..............................................................................

      end module module_airseaflux
