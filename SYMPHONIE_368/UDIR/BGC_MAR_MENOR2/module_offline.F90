      module module_offline
!______________________________________________________________________
! SYMPHONIE Ocean Model
! release 367 - last update: 07-03-23
!______________________________________________________________________
!______________________________________________________________________
! Version date      Description des modifications
!         05/03/01: mise en service
!         08/01/03: correction pour eviter inconvenient des parties
!                   entieres sur nombre negatifŋ8029
!         27/01/03: VELOBC remplace VHZOBC
!         13/07/05: suite au bug trouvé par katell qui a montré qu'on ne
!                   pouvait pas faire marcher les bouees avec les options
!                   offline et forward on remplit tous les tableaux de vitesse
!                   quand l'option forward est signalée
!         10/05/06: fonctions compatibles avec double precision
!         19/09/06: 1 Possibilite de filtrage du courant de petit echelle quoique
!                   ce sous programme, filtre_uvprime, soit commenté par defaut.
!                   2 MOdification pour conservation du toit plus rigoureuse
!         05/03/07: debug du point precedent (supprime la multiplication par
!                   morpho incompatible avec les sources des fleuves)
!         03/04/07: ajout d'une securité pour que l'archivage ne puisse pas
!                   être active en même temps qu'une imbrication (utilise la
!                   même memoire vive)
!         21/04/07: passage à coordonnee curviligne
!         02/10/07: choix=-1 rien que pour connaitre ioffline suffisement tot
!                   dans la phase d'initialisation
!         13-06-09  - entrees/sorties au format netcdf
!                   - des boucles plus etendues pour calculs zta
!         15-06-09  Pour qu'xscan fonctionne on a revu la taille des boucles
!                   pour u et v
!         22-06-09  ecrire dxdy_c
!                   DZ_x(i,j)=0.5*(DZ_z(i-1,j)+DZ_z(i,j))
!                   DZ_y(i,j)=0.5*(DZ_z(i,j-1)+DZ_z(i,j))
!         23-06-09  Masquer T et S en i=0,i=imax+1,j=0,j=jmax+1
! 2009.3  02-10-09  Amenagements divers pour passer aux nouvelles variables
!         22-10-09  - dans cette version toutes les variables sont archivées de
!                   0 à imax+1 et de 0 à jmax+1
!                   - on archive le facteur d'echelle verticale à chaque echeance
!         23-10-09  suite du point precedent (modif dans la partie "lecture")
!         08-11-09  les dimensions des fichiers netcdf sont declarees dans
!                   module_principal
!         16-11-09  introduction du temps ecoulé dans les fichiers netcdf
! 2010.2  22-12-09  - ecrire dans tmp
!                   - graph_out_mod renommé module_graph
! 2010.6  02-02-10  renomme lon_t lat_t
!         05-02-10  blindage division par sumarchive
!         11-02-10  notebook_offline = nomfichier(20)
! 2010.7  16-02-10  vitesse forward = veldxdz(2)
!         05-03-10  'axis' remplacé par 'content'
! 2010.8  03-05-10  Nouveau schema forward
!         09-05-10  seul le proc 0 ecrit messages
! 2010.9  29-05-10  suppression offline_filtre_uvprime
!         11-06-10  Verifier que la date donnée dans notebook_offline n'est pas
!                   anterieure au debut du run
! 2010.10 16-06-10  Suite point precedent + modifs sur calcul omega et
!                   facteurs d'echelle
!         24-06-10  les noms des fichiers commence par 's'
! 2010.11 15-07-10  "status" declaré dans module_principal
!         18-07-10  lecture removetide
! 2010.12 27-09-10  suppression arg "key" dans appel à routine omega
! 2010.13 21-10-10  A la demande des utilisateurs qui veulent pouvoir exploiter
!                   facilement la vitesse (plutôt que le transport) en mode
!                   "lecture offline" pour faire des diags, une interpolation
!                   interpolation des vitesses est ajouté au cas case_=2
!         01-11-11  mpi: visualiser la frontiere des sous-domaines
!         03-11-10  des arguments passés dans date_to_kount
! 2010.14 05-11-05  ecriture du temps en double precision
!         11-11-10  modif ecriture du temps (utilisation de fonction nint)
!         18-11-10  ajout de la bathy+mask dans le fichier de grille
!         21-11-10  ecrire upwindriver_t dans le fichier de grille
! 2010.18 21-03-11  Suppression d'un attribut mal interprete par FERRET
!         27-03-11  Mode lecture adapté au fichier unique
! 2010.19 30-03-11  - Possibilité de mettre à jour les listes de fichiers en
!                   cours de simulation (applications operationnelles)
!                   - Gain cpu: inutile de calculer T et S quand le modele
!                   biogeochimique n'est pas activé.
! 2010.20 17-04-11  Calculs sur la base d'un temps en secondes
!         19-04-11  k2dfin renommé iteration2d_max_now
! 2010.21 20-04-11  Ecrire une liste unique de fichier "offline"
! 2010.22 05-05-11  Ecrire le numero du rank dans variable mpi_t
! 2010.23 12-05-11  Ecriture variable mpi_t: differencier par%rank=0 du masque
!                   continental
!         17-05-11  Oter la pression atmospherique de la ssh archivée dans le fichier
!                   netcdf
!         02-06-11  model_name porte le nom du modele
! 2010.24 23-11-11  stockage dzofl: blindage de la condition sous le fond
! 2010.25 08-02-12  kz_w renomme kh_w
!         26-02-12  allocation dynamique
!         17-03-12  commencer a prevoir d'archiver les lon lat en double precision
!         03-04-12  routine pour attributs generaux des fichiers netcdf
!         06-04-12  modifs pour norme netcdf comodo
!         21-06-12  Passage aux bibliotheques pnetcdf
!         27-06-12  attribut axis "T" pour le temps
! S25     07-09-12  debug lecture netcdf
! S26     26-09-12  Ajout (si demande dans notebook_offline) du courant eastward et northward
!         04-11-12  lecture des champs nemo
!         03-02-13  correction d'un bug introduit à la modif precedente
!         06-02-13  ecrire lon lat des points "f"
!         02-04-13  Debug runoff nemo
!         07-04-13  Debug runoff nemo
!         23-06-13  mangrove
!         03-07-13  - ssh-ib devient ssh_ib (pas de signe "-" svp!!!!)
!                   - Lecture des fichiers dans routine offline_read_file. Pour economiser
!                   calculs la nouvelle routine ne lit que la prochaine echeance et non
!                   pas les 2 echeances "avant" et "apres".
!         04-07-13  Seul par%rank=0 ecrit liste binrec
!         26-07-13  Ecriture de i et j dans le fichier grid.nc
!         27-07-13  Ecriture sqrt_dxdy masquee pour graphiques....
!         01-11-13  Format "short", archivage temporel irregulier
!         10-11-13  ecrire par%subcycle et dti_fw dans fichier grille
!         04-12-13  boucle sur nriver associee a un test de "localisation mpi"
!         29-12-13  De nouveaux champs dans le fichier tmp/grid.nc (angle de axes,...)
!         07-01-14  ajout dx_t masquee et dy_t masquee
!         28-03-14  ajouts pour nemo offline
!         02-04-14  time real.....
!         19-04-14  river_inout renomme rivertrc_inout
!         29-04-14  debug mpi
!         30-04-14  Si e3u(e3v) n'est pas present dans le fichier nemo, on utilise la valeur
!                   initiale
!         01-05-14  Ecrire la grille en double precision pour relecture precise par la
!                   procedure offline
!         02-05-14  Un type special pour la SSH dans l'ecriture du fichier offline
!         09-05-14  nc_==1 devient nc_<=1
!         17-05-14  modif appel EOS
!         10-06-14  ecriture gridrotcos_t gridrotsin_t dans fichier grille
!         19-06-14  eqs_state1 devient eos_author
!         26-06-14  modif ecritude CFL 2D
!         27-06-14  ecrire le mask en format short
!         09-07-14  notebook_offline distingue le nom de la directory des fichiers offline
!                   et le nom de la liste contenant les fichiers offline
!         11-07-14  - rap_obc devient timeweightobc
!                   - debug offline_nemo_obcext
!         13-07-14  notebook_offline partitionne en 2 zones.
!         15-07-14  nouveaux echanges
!         16-07-14  l'experience montre que s'appuyer sur iteration3d=0 est un
!                   test faillible pour offline_read car la bio peut faire des appels
!                   qui lui sont propres, d'ou offline_init_status
!         20-07-14  debug signe obc nemo offline obc nord & est
!         21-07-14  utiliser directory_offline pour acceder aux fichiers offline
!         23-07-14  - update dz deplace avant vel=veldxdz/dy/dz pour coherence entre
!                   veldxdz et nbre de courant du schema d'advection bio
!                   - Seuillage ssh_int pour eviter decouvrement du point ssh source river
!         25-07-14  Alarme/arret dans le cas ou la vitesse w n'est pas trouvee dans les fichiers
!                   nemo pour le cas nemo-offline
!                   Ajout des flux de surface moyens dans fichiers netcdf
!         01-08-14  debug test sur ofl_surflux
!         27-09-14  ajout offline_kz_convect_adjust qui augmente le Kz de nemo-mercator
!                   lorsque le gradient vertical de densite potentielle est nul
!         20-10-14  EOS: supprime possibilite de passer par densite potentielle zref local
!         22-10-14  affichage de i.j dans grid.nc
!         01-11-14  modif du critere de correction du kz de nemo en mode offline
!         02-11-14  blindage fonction sqrt
!         03-11-14  ofl_type et ofl_sshtype definis dans notebook_offline
!                   ecriture mask_u et mask_v dans grid.nc
!         05-11-14  blinder log10
!         28-11-14  ne pas ecrire sigma_w dans grid.nc s'il n'est pas declare
!         17-12-14  sigma_w remplacé par dsig_t
!         18-12-14  modif boucle pour convenance graphique
!         20-12-14  utiliser masque de surface pour archive salofl temofl
!         08-01-15  Pour convenances de comparaison de la ssh offline reconstituee
!                   avec celle des fichiers, on prefere ajouter le BI a sshobc
!         16-01-15  texte60 commente
!         02-02-15  modif long name gridrotcos et gridrotsin
!         06-02-15  sshofl_w borne decouvrement
!         05-03-15  ajout ofl_tke pour archive de la tke dans fichiers "offline"
!         16-03-15  call obc_int_anyvar2d parallelise anyvar2d pour esthetique graphique
!         19-03-15  debug reset tke_ofl
!         20-03-15  modif sur test _Fillvalue
!         16-04-15  La subroutine offline_read_file devient un driver de plusieurs
!                   sous-subroutines
!         18-04-15  Plus de message a l'ecran pour debugage occasionnel
!         23-04-15  Pour ne pas coincer la lecture du fichier restart on prevoit
!                   une possibilite de desalloction ordonnee depuis subroutine dyn_restart
!         25-05-15  calcul seconds_cpu_
!         26-06-15  call offline_write_oasisgrids
!         29-06-15  suite du point precedent (correction du contenu du fichier
!                   de grille)
!         05-07-15  Deux nouveaux champs sont stokes dans les fichiers offline.
!                   Le flux d'eau a la surface de l'ocean (dans la foulee de 
!                   l'apport d'eau douce dans la simulation directe) et d'autre part
!                   on archive la ssh instantanee delimitant l'intervale d'integration
!                   du transport pour pouvoir faire une verification rigoureuse de
!                   la variation de ssh en mode offline. La ssh est masquee pour
!                   eviter l'abaissement en amont des points sources des rivieres. 
!         17-10-15  ajout i_index et j_index pour points u et v dans le fichier grid.nc
!         07-11-15  ajout dans grid.nc de smoothed_h_minus_h
!         18-11-15  Ecrire H-Hnemo dans fichier de grid.nc
!         20-12-15  Bornes depth_u depth_v
!         22-12-15  test mpi_proc_null
!         03-03-16  boucles individuelles sur fluxbar_u fluxbar_v
!         24-03-16  ajout mask_w dans grid.nc
!         04-04-16  Debug mask mangrove
!         08-04-16  Ajout h_f
!         15-04-16  ti0tide devient un tableau 1D
!         21-06-16  debug lecture du fichier de vitesse verticale de
!                   mercator
!         01-07-16  Ajout offline_another_model_divergence_2 la version 2 de la surface
!                  libre du cas nemo bio offline. Plus chere mais plus efficace pour reduire 
!                  l'ecart entre la ssh diagnostiquee et la ssh de nemo
!         05-01-17 Bathy nemo dans grid.nc
!         18-01-17 Ajout dz_u et dz_v dans fichier de grille
!         19-01-17 debug point precedent
!         21-01-17 fichier netcdf modif attribut longname de la densite
!         07-02-17 message d'aide au debugage
!         09-02-17 ajout subroutine offline_oldlist
!         28-02-17 ajout d'une extension au nom des fichiers netcdf
!         07-03-17 reset texte80(11)
!         21-03-17 debug du nom 'netcdf' de la variable v_sn
!         25-04-17 ajout angle de grille dans grid.nc
!         05-05-17 seul rank=0 lit fichier puis broadcast
!         15-06-17 Couche "convective" on applique un fort Kz pour la simu bio offline 
!         21-07-17 modification sur attribut netcdf
!         09-11-17 ajout subroutine maximum_bottom_stress
!         17-11-17 lecture maxbotstress
!         21-01-18 pouvoir changer la periodicite de l'ecriture apres un restart
!         25-03-18 procedure stop pour occigen
!         26-03-18 adaptation A la grille verticale fusionnee
!         08-05-18 modification condition d'allocation du tableau oflshort
!         09-05-18 procedure stop pour occigen
!         22-05-18 modele offline demarre sur fichier "lonlatfile"
!         27-05-18 nf_clobber+ NF_64BIT_DATA
!         08-06-18 si flag_offline_binary=1 alors fichiers offline au format binaire
!         18-06-18 subroutine offlinegbin
!         05-07-18 nf_clobber+ NF_64BIT_DATA permet de creer des fichiers pour de grandes grilles
!         21-07-18 retour A NF_64BIT_OFFSET
!         23-10-18 ajout sig1dpom_t dans grid.nc
! v245    08-02-19 ne plus appliquer l'algo d'augmentation de difofl entre ksl et kmax
! v246    09-02-19 ajout kslofl_t
! v251    10-04-19 ajout sources sous marine
! v252    23-04-19 ajout reseau de canaux
! v253    06-05-19 hz_w(i,j,2)=h_w(i,j)+ssh_int_w(i,j,2) deplace apres echanges canaux de ssh_int_w
!         07-05-19 Pour ne pas diviser par dz=quasizero  on moyenne separement les couches fusionnees
!         08-05-19 - deplacemet du reset de ssh_int pour iteration3d=0
!                  - flux3d nuls dans couches ecrases, correction de la
!                    moyenne proportionnelle au flux3d ebauche
!                  - Diag de verification de la divergence 
! v255   25-05-19  Tous les coeurs font une partie de la liste. Dans le
!                  cas symphonie (et pas nemo) faire qu'une liste et la
!                  recopier sur les autres champs
!        26-05-19  Inverser l'ordre des fichiers pour pouvoir calculer
!                  des retrotrajectoires en mode offline
!        06-06-19  suite du point 25-05-19
! v257   05-07-19  zone1_mask dans grid.nc
! v259   09-09-19  ajout bioofl_t
!        11-09-19  ajout hsedofl_t
! v260   05-10-19  - ecrire tmp/grid.nc mEme en mode offline 2
!                  - Couches ecrasEes mode offline. Desormais le courant
!                  peut y etre non nul (pour ne pas y pieger des drifters) mais par
!                  contre on n'ajuste toujours pas le transport dans ces
!                  couches car peut Etre risquE du point de vue de la precision
! v261   22-10-19  Amelioration du nom des fichiers (comme dans graph_out.F90)
! v269   06-12-19  ajout kmergedr4_u et v
! v275   18-02-20  ajout z0_w dans fichier de grille
! v276   13-03-20  limiter la taille du fichier de grille en evitant les
!                  variables 3D inutiles
! v279   25-04-20  ajout d'un fillval
! v284   22-05-20  ajout upwzone0_t dans grid.nc
! v285   07-06-20  call mpi_bcast(texte30,len(texte30),mpi_character,0,par%comm2d,ierr) !07-06-20
! v287   17-07-20  utiliser tmpdirname
! v289   24-09-20  zone_mask renommE zone1_mask
! v290   20-10-20  detider uah et vah
! v294   17-12-20  ssh_w remplace ssh_int_w
! v295   01-02-21  ajoute de tauw et maxbotstress dans les fichiers offline
! v303   14-07-21  ajout de dz_u, dz_v, vel_u(0), vel_v(0) pour pouvoir calculer velbot en offline
! v309   27-09-21  correction attribut netcdf unites
! v310   13-10-21  modifs sur champ netcdf CFL2D
!        22-10-21  ajout kmerged_u , kmerged_v
! v328   23-02-22  appels aux EOS comme dans presgrad
! v351   15-08-22  convertir les fichiers offline binaires en fichiers netcdf
! v352   18-08-22  ecrire dx_t etc... en double precision
! v362   06-01-23  mises A jours Claude pour anr popnco
! v367   07-03-23  augmenter taille dde boucle lecture temobc, salobc
!............................................................................
!  _________                    .__                  .__              ! m°v°m
! /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____        !
! \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \       !
! /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/       !
!/_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >      !
!        \/\/          \/|__|        \/            \/        \/       !
!......................................................................
      use module_principal ; use module_parallele ; use module_s ; use module_webcanals
      use module_my_outputs ; use module_biology
      use sedimento , only : dzs,ksmi,ksma
      use pnetcdf
      implicit none
      real :: scalarval(1)
      character*200 filename_
      character*60 txt_units
      character*33 time_units
      character*15 :: date_shortname  ='000000000000000' &
                     ,date_shortname_b='000000000000000'
      integer :: recordlength_=250 ! 204+2*8+6*2
      integer :: ogcmtimecounter_,time_,time2_,loop_ &
              ,binaryfilelist_len,flag_ncfile=0
      integer kbegin_,kend_,kstep_,type_
      integer :: i1val(1)
      integer(kind=MPI_OFFSET_KIND) start(4)
      integer(kind=MPI_OFFSET_KIND) edge(4)
      integer(kind=MPI_OFFSET_KIND) :: start_mpikind,count_mpikind=1
      integer var_dims,var_type,tabdim(5) ! Que represente tabdim?
                                          ! Soit 1,2,3, l'ordre dans lequel
                                          ! sont listees les dimensions du fichiers
                                          ! netcdf. Tabdim indique (via le numero
                                          ! d'ordre) les dimensions dont dépend la
                                          ! variable
      integer,dimension(:,:,:),allocatable :: oflshort
      integer,dimension(:),allocatable :: binaryfile_inout,imaxbf,jmaxbf
      real,dimension(:,:),allocatable :: binaryfile_tab2d_r4
      real,dimension(:,:,:),allocatable :: binaryfile_tab3d_r4
      integer(kind=2),dimension(:,:,:),allocatable :: binaryfile_tab3d_i2
      character,dimension(:),allocatable :: binaryfile_name*200
      double precision dte_nemo_small,dte_nemo_big
      integer :: loop_nemo_ssh,flag_newalgo=0,offline_gridbin_status=0

! https://docs.google.com/document/d/1qCWQrn6EmHbOk1mlVIslo71GlCpBNx43LA9_wWrP1p8/edit

contains

!..............................................................................

      subroutine offline_inout(case_)
      implicit none
      integer case_
#ifdef synopsis
       subroutinetitle='offline_inout'
       subroutinedescription= &
       'Driver of the in/out "offline" subroutines'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! DRIVER:

!     if(par%rank==0)write(6,*)'BIDOUILLE offline_inout'
!     RETURN

      if(case_==2) then !-read-read->
! Physic is not computed, load the "offline" files
       call offline_read(case_) ; return
      endif             !-read-read->

      if(case_==1) then !-write-write->
! Physic is computed, "offline" files are created:
       call offline_write(case_) ; return
      endif             !-write-write->

      if(case_==0.or.case_==-1) then !iiiiiii>
! Initialization step: read notebook, allocate,...:
       call offline_initial(case_) ; return
      endif                          !iiiiiii>

      end subroutine offline_inout

!..............................................................................

      subroutine offline_initial(case_)
      implicit none
      integer case_,nc_
      double precision period_
#ifdef synopsis
       subroutinetitle='offline_initial'
       subroutinedescription='Reads notebook_offline'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! lecture de notebook_offline

      if(par%rank==0)write(6,'(a,a)')'Lecture de ',trim(nomfichier(20)) !11-02-10

      open(unit=3,file=nomfichier(20)) ! lecture de notebook_offline
! Passer par dessus les lignes de la zone 'namelist' pour lire en mode 'fichier' !13-07-14
! Note: la zone 'namelist' a ete lue par la routine set_parameters.F90
 1964 read(3,'(a)')texte90
      if(texte90(1:57)/=                                           &
      'Periodicity (hours) ! until yyyy / mm / dd / hh / mm / ss')goto 1964

      if(case_==-1.and.ioffline/=0)call offline_allocate('a') !26-02-12

      if(case_==-1)  goto 1000 ! On a IOFFLINE, sortie rapide.       !02/10/07
      if(ioffline==0)goto 1000 ! sortie rapide si procedure offline desactivee

      if(ioffline==1.and.iarchive==1) then
       if(par%rank==0)then
        write(6,*)'incompatibilité entre notebook_offline'
        write(6,*)'et notebook_obcforcing. on ne peut à la'
        write(6,*)'fois archiver des champs pour une future simulation'
        write(6,*)'physique imbriquée et une future simulation offline.'
        write(6,*)'faire le choix puis relancer.'
       endif
       stop ' dans offline_inout.f'
      endif

! Lire les dates figurant dans la zone 'fichier' du notebook_offline et construire
! le liste format fichier binaire:
      if(par%rank==0)open(unit=4,file=trim(tmpdirname)//'ofl_list'  &
                                ,access='direct'      &
                                ,recl=16              &
                                ,form='unformatted')

      nc_=0
  212 read(3,*,end=211)period_,i1,i2,i3,i4,i5,i6
      call datetokount(i1,i2,i3,i4,i5,i6)

      if(elapsedtime_out>elapsedtime_now) then !-------->
       nc_=nc_+1
       ofl_rec_max=nc_
       if(par%rank==0)write(4,rec=nc_)period_*3600.,elapsedtime_out
      endif                                    !-------->

      if(nc_<=1) then !1111111>                     ! 09-05-14
         ofl_period_now= period_*3600.              ! periodicite en cours
         ofl_period_next=period_*3600.              ! periodicite a suivre
         ofl_writime=elapsedtime_now+ofl_period_now ! Ecriture du prochain fichier necdf
         ofl_rec_now=1                              ! record de la periodicite en cours
         ofl_nextrec_time=elapsedtime_out           ! temps du changement de periodicite
         if(nc_==1.and.par%rank==0) then !........>
          write(6,*)'ofl_rec_now     =',ofl_rec_now
          write(6,*)'ofl_period_now  =',ofl_period_now
          write(6,*)'ofl_period_next =',ofl_period_next
          write(6,*)'ofl_writime     =',ofl_writime
          write(6,*)'ofl_nextrec_time=',ofl_nextrec_time
          write(6,*)'nc_             =',nc_
         endif                           !........>
      endif          !1111111>

      goto 212

  211 continue

      if(par%rank==0)close(4)
  213 continue
 1000 continue
      close(3)

! en mode lecture des fichiers le centre de gravité temporel est décallé
! d'une demi période en arrière:
      if(ioffline.eq.2) then !***>
!     if(index(nomfichier(20),'another_model')==0) then !-->           !04-11-12
      if(flag_nemoffline==0) then                       !-->
         offlinedt(2)=offlinedt(2)-offlinedt(1)/2.
      endif                                             !-->
      endif                  !***>

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)
#endif

      if(ioffline==0)return
      if(case_==-1)return                                           !02/10/07

      if(par%rank==0) then !#mpi-->>-->                       !09-05-10
      open(unit=3,file=trim(tmpdirname)//'messages',position='append')
      write(3,*)'-----------------------------------------------------'
      write(3,*)'subroutine offline_inout:'
      write(3,*)
      write(3,*)'procedure offline activée'
      write(3,*)'date du premier archivage en s:',offlinedt(2)
      write(3,*)'periode d`archivage en s:',offlinedt(1)
      if(ioffline.eq.1)                                                 &
      write(3,*)'phase préalable de création des fichiers'
      if(ioffline.eq.2)write(3,*)                                       &
      'la dynamique n`est plus calculée mais lue dans des fichiers'
      close(3)
      endif                !#mpi-->>-->                       !09-05-10
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !09-05-10
#endif

! Faire la liste format binaire acces direct:
      if(ioffline==2)call offline_readthelists(0) !30-03-11

      if(ioffline.eq.1) then !11111111111111111>

! Si phase ecriture, ecrire le fichier de grille:
      if(mpi_hole_plugging=='none')call offline_write_grid(2) ! arg 2 fichier offline/grid.nc

      if(par%rank==0) then !0000000000>                                !20-04-11
       write(6,*)'si le model_ plante ici c''est que le FICHIER'
       write(6,*)trim(offlinefile(1)),' n''existe pas'
       open(unit=4,file=trim(offlinefile(1)))                          !04-11-12
       write(4,'(a)')'files created by the offline mode:'
       close(4)
       write(6,*)'ok'
      endif                !0000000000>                                !20-04-11
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)
#endif
      endif                  !1111111111111111>

      end subroutine offline_initial

!..............................................................................

      subroutine offline_write(case_)
      implicit none
      integer case_,iostatus_
#ifdef synopsis
       subroutinetitle='offline_write'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Mise à jour de la periodicite?:
!     if(ofl_rec_now<ofl_rec_max) then !)))))))))>
       if(elapsedtime_now>=ofl_nextrec_time) then !>>>>>>>>>>>
  437   ofl_rec_now=ofl_rec_now+1
        open(unit=4,file=trim(tmpdirname)//'ofl_list'  &
                   ,access='direct'      &
                   ,recl=16              &
                   ,form='unformatted')
        read(4,rec=ofl_rec_now,IOSTAT=iostatus_)ofl_period_next,ofl_nextrec_time
        if(     iostatus_/=0  &      ! iostatus=0  signifie pas d'erreur
           .and.iostatus_/=36 &      ! iostatus=36 signifie fin du fichier (on garde les valeurs precedentes de ofl_period_next et ofl_nextrec_time)
          ) stop 'err 513 ofl_list' ! iostatus/=0 et /=36 est une erreur d'un autre type qui impose l'arret du modele

        close(4)
        if(ofl_nextrec_time<=elapsedtime_now)goto 437 ! Si prochaine echeance anterieure au temps present voir echeance suivante !21-01-18

        if(ofl_period_now<=0.and.ofl_period_next>0)then !spspsp>
! Cas particuier: la phase de spin-up (pas d'ecriture) vient d'etre depassee:
         ofl_period_now=ofl_period_next
         ofl_writime=elapsedtime_now+ofl_period_now
        endif                                           !spspsp>

       endif                                      !>>>>>>>>>>>
!     endif                            !)))))))))>

      if(ofl_period_now<=0)return ! on ne fait rien tant que la phase de spin-up n'est pas depassee


      if(elapsedtime_now>=ofl_writime) then !1111111111111111111111111>

      call s_cpu('before_offline',1)

!.....calcul de la date                                                    !27/01/03
      call elapsedtimetodate(elapsedtime_now-0.5*ofl_period_now  &
                             ,i5,i6,i7,i3,i2,i1)

      write(texte80(1),                            &
       '(1x,i2,1x,a9,1x,i4,1x,a,i2,a1,i2,a1,i2)')  &
       i7,month(i6),i5,'h:m:s ',i3,':',i2,':',i1

!.....puis afficher la date à l'écran:                                     !27/01/03
      if(par%rank==0)write(6,*)'--------------------------------'
      if(par%rank==0)write(6,*)'fichier offline date:'
      if(par%rank==0)write(6,'(a33)')texte80(1)(1:33)

!.....écrire année année mois dans TEXTE90:                                !10/07/03
      i0=i7+100*i6+10000*i5
      write(texte90(1:8),'(i8)')i0

!.....écrire heure seconde minute dans TEXTE90:
      i0=i1+100*i2+10000*i3
      i0=i0+1000000 ! on ajoute cette constante pour forcer            !11/07/03
                    ! l'ecriture des caracteres "0" qd I0=0
      write(texte90(9:15),'(i7)')i0
!.....puis on efface le "1" avec un caractere de séparation, "_"
!     écrire ":" dans TEXTE90:
      write(texte90(9:9),'(a1)')'_'                                    !28/07/04

      if(par%rank==0) then !m[°u°]m> !05-05-17
       k=s_unit(7) 
       open(unit=k,file='output_file_extension')
        texte30='' !22-10-19
        read(k,*,end=471)texte30
  471  close(k)
      endif                !m[°u°]m> 
      call mpi_bcast(texte30,len(texte30),mpi_character,0,par%comm2d,ierr) !07-06-20

! Nom des fichiers de variables
!     k=index(offlinefile(1),'/',back=.true.)                                    !03-02-13
      texte250=trim(directory_offline)//texte90(1:15) &
                                      //trim(texte30) &           !28-02-17
                                      //'.symphonie.nc'//char(0)  !21-07-14!22-10-19
      if(flag_offline_binary==1)date_shortname=texte90(1:15)
      if(par%rank==0) then !00000000000> !20-04-11
       open(unit=4,file=trim(offlinefile(1)) &    !04-11-12
                  ,position='append')
       k=index(texte250,".nc")
       write(4,'(a)')texte250(1:k+2)
       close(4)
      endif                !00000000000> !20-04-11
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes
#endif

               call offline_endofaverage !18-07-10
               call offline_write_var  ! en mode "time unlimited"


! RESET GENERAL:
                                 temofl_t(:,:,:,1)=0.
                                 salofl_t(:,:,:,1)=0.
               if(ale_selected==1)dzofl_t(:,:,:,1)=0.
                                 velofl_u(:,:,:,1)=0.
                                 velofl_v(:,:,:,1)=0.
                         fluxbarsum_ofl_u(:,:)=0.
                         fluxbarsum_ofl_v(:,:)=0.
                                 sshofl_w(:,:,1)=0.
                                 dfvofl_w(:,:,:,-1)=0.
           if(flag_ksloffline==1)kslofl_t(:,:,1)=0. !09-02-19
                                 w0mofl_w(:,:,1)=0.
                   if(ofl_tke==1)tkeofl_w(:,:,:,-1)=0.
                   if(ofl_bio==1) then !bbb>
                          bioofl_t=0. !09-09-19
                          hsedofl_t=0. !11-09-19
                   endif               !bbb>
      if(flag_groundwater==1)w_keq1_ofl_w(:,:,1)=0. !10-04-19

       if(ofl_surflux==1) then !---> !01-08-14
        slhf_aver_w(:,:)=0.
        sshf_aver_w(:,:)=0.
        snsf_aver_w(:,:)=0.
        ssr_aver_w(:,:)=0.
        precipi_aver_w(:,:)=0.
        wstress_aver_u(:,:)=0.
        wstress_aver_v(:,:)=0.
       endif                   !--->
      if(flag_maxbotstress==1) then ! m[°v°]m > !09-11-17
         maxbotstress_w=0.
         stresswave_w=0.    ! 01-02-21
         stressc_w=0.
      endif                         ! m[°v°]m >
       
       sumarchive=small1
       ofl_period_now=ofl_period_next
       ofl_writime=ofl_writime+ofl_period_now

       call s_cpu('offline_written',1)

      endif        ! 1111111111111111111111111111111111111111111111111>

!...............................................................................
! L'INTEGRATION TEMPORELLE DES VARIABLES COMMENCE UNE ECHEANCE AVANT NC=1
! DEBUT:
      if(elapsedtime_now>=ofl_writime-ofl_period_now) then      !4444444>  !17-04-11
!...............................................................................

      sumarchive=sumarchive+dti_now     !17-04-11
      do k=1,kmax
       do j=1,jmax ; do i=1,imax
         temofl_t(i,j,k,1)=temofl_t(i,j,k,1)+tem_t(i,j,k,1)*dti_now !17-04-11
         salofl_t(i,j,k,1)=salofl_t(i,j,k,1)+sal_t(i,j,k,1)*dti_now !17-04-11
        enddo ; enddo
       do j=1,jmax ; do i=1,imax+1
        velofl_u(i,j,k,1)=velofl_u(i,j,k,1)+veldydz_u(i,j,k,1)*dti_now !17-04-11
       enddo ; enddo
       do j=1,jmax+1 ; do i=1,imax
        velofl_v(i,j,k,1)=velofl_v(i,j,k,1)+veldxdz_v(i,j,k,1)*dti_now !17-04-11
       enddo ; enddo
      enddo

      const1=-1./grav/rho     !17-05-11
      do 28 j=0,jmax+1
      do 28 i=0,imax+1
! dans cette version prendre zta_ext_z(i,j,2) car c'est la seule à
! etre calculee sur les c.l. de type "Z0". Dans symphonie2009,
! on peut prendre zta_int
       sshofl_w(i,j,1)=sshofl_w(i,j,1)              &
!      +dti_now*(     ssh_int_w(i,j,1)              &
       +dti_now*(         ssh_w(i,j,1)              & !17-12-20
                 -const1*(pss_w(i,j,1)-pss_mean(1)) & !17-05-11 ! OTER LE B.I.
                     )
       w0mofl_w(i,j,1)=w0mofl_w(i,j,1)+dti_now*omega_w(i,j,kmax+1,1) !04-07-15
   28 continue
      if(flag_groundwater==1) then !oooo> !11-04-19
       do j=1,jmax ; do i=1,imax
        w_keq1_ofl_w(i,j,1)=w_keq1_ofl_w(i,j,1)+dti_now*omega_w(i,j,1,1) !11-04-19
       enddo         ; enddo
      endif                        !oooo> !11-04-19

      if(flag_ksloffline==1) then !m°v°m>
       do j=1,jmax   ; do i=1,imax  
        kslofl_t(i,j,1)=kslofl_t(i,j,1)+dti_now*ksl_t(i,j)            !09-02-19
       enddo         ; enddo
      endif                       !m°v°m>

      do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
       dfvofl_w(i,j,k,-1)=dfvofl_w(i,j,k,-1)+kh_w(i,j,k)*dti_now !17-04-11
      enddo ; enddo ; enddo
! Commentees le !08-02-19
!     do j=1,jmax ; do i=1,imax
!      do k=kmax,ksl_t(i,j),-1 !Couche "convective" on applique un fort Kz pour la simu bio offline !15-06-17
!       dfvofl_w(i,j,k,-1)=dfvofl_w(i,j,k,-1)+10.*dti_now 
!      enddo
!      do k=ksl_t(i,j)-1,kmin_w(i,j)+1,-1 ! couche standard
!       dfvofl_w(i,j,k,-1)=dfvofl_w(i,j,k,-1)+kh_w(i,j,k)*dti_now !17-04-11
!      enddo
! Note en k=kmax+1 et k=kmin_w kz=0
!     enddo ; enddo

      do j=1,jmax ; do i=1,imax+1
       fluxbarsum_ofl_u(i,j)=fluxbarsum_ofl_u(i,j)+fluxbar_sumt_u(i,j,1)*dti_now !17-04-11
      enddo ; enddo
      do j=1,jmax+1 ; do i=1,imax
       fluxbarsum_ofl_v(i,j)=fluxbarsum_ofl_v(i,j)+fluxbar_sumt_v(i,j,1)*dti_now !17-04-11
      enddo ; enddo

      if(ale_selected==1) then !-ale-ale-ale-ale->
       do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
        dzofl_t(i,j,k,1)=dzofl_t(i,j,k,1)+dz_t(i,j,k,1)*dti_now !17-04-11
       enddo ; enddo ; enddo
      endif                    !-ale-ale-ale-ale->
      if(ofl_tke==1) then      !tke-tke-tke-> !05-03-15
       do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
        tkeofl_w(i,j,k,-1)=tkeofl_w(i,j,k,-1)+tken_w(i,j,k)*dti_now !17-04-11
       enddo ; enddo ; enddo
      endif                    !tke-tke-tke->

      if(ofl_bio==1) then      !bio-bio-bio-> !09-09-19
       if(ubound(bioofl_t,4)/=vbmax) then 
         deallocate(bioofl_t)
         allocate(bioofl_t(0:imax+1,0:jmax+1,kmax+1,vbmax))
         bioofl_t=0
       endif
       do vb=1,vbmax ; do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
        bioofl_t(i,j,k,vb)=bioofl_t(i,j,k,vb)+bio_t(i,j,k,vb)*dti_now 
       enddo ; enddo ; enddo ; enddo
       ! epaisseur de couche ! 11-09-19
       do j=1,jmax ; do i=1,imax ; do k=ksmi(i,j),ksma(i,j)
         hsedofl_t(i,j)=hsedofl_t(i,j) + dzs(k,i,j)*dti_now
       enddo ; enddo ; enddo
      endif                    !bio-bio-bio-> !09-09-19

      if(ofl_surflux==1) then !---> !01-08-14
       do j=0,jmax+1
       do i=0,imax+1
           slhf_aver_w(i,j)=slhf_aver_w(i,j)+dti_now*slhf_w(i,j,1)
           sshf_aver_w(i,j)=sshf_aver_w(i,j)+dti_now*sshf_w(i,j,1)
           snsf_aver_w(i,j)=snsf_aver_w(i,j)+dti_now*snsf_w(i,j,1)
            ssr_aver_w(i,j)=ssr_aver_w(i,j)+dti_now*ssr_w(i,j,1)
        precipi_aver_w(i,j)=precipi_aver_w(i,j)+dti_now*precipi_w(i,j,1)
       enddo
       enddo
       do j=1,jmax
       do i=1,imax+1
        wstress_aver_u(i,j)=wstress_aver_u(i,j)+dti_now*wstress_u(i,j,1)
       enddo
       enddo
       do j=1,jmax+1
       do i=1,imax
        wstress_aver_v(i,j)=wstress_aver_v(i,j)+dti_now*wstress_v(i,j,1)
       enddo
       enddo
      endif                   !--->

!Cumul du stress de fond max
      if(flag_maxbotstress==1)call dragcoef_maxbotstress !09-11-17


!...............................................................................
! L'INTEGRATION TEMPORELLE DES VARIABLES COMMENCE UNE ECHEANCE AVANT NC=1
! FIN.
      endif                                                   !4444444>
!...............................................................................

      end subroutine offline_write

!...............................................................................

      subroutine offline_read(case_)
      implicit none
      integer case_
#ifdef synopsis
       subroutinetitle='offline_read'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!     if(iteration3d==0) then !-initial-case->
      if(offline_init_status==0) then !-initial-case-> !16-07-14


! Lecture champ "0"
       nc=ofl_rec_now
       call offline_read_file
       ofl_rec_now=ofl_rec_now+1


! Lecture champ "2"
       ofl_readtime_prev=ofl_readtime_next
       ofl_period_prev=  ofl_period_next
       nc=ofl_rec_now
       call offline_read_file
       ofl_rec_now=ofl_rec_now+1

       offline_init_status=1 ! initial phase is done 

      else                    !-iterative-phase->

       if(elapsedtime_now>ofl_readtime_next) then !>>>>>>>>>>

        ofl_readtime_prev=ofl_readtime_next
        ofl_period_prev  =ofl_period_next
        nc=ofl_rec_now
        call offline_read_file
        ofl_rec_now=ofl_rec_now+1

       endif                                      !>>>>>>>>>>


      endif                   !-iterative-phase->


!......................................................................
! interpolation temporelle de la physique entre 2 archivages:
! debut:
!......................................................................
! Pourcentage des 2 écheances:
       if(ofl_period_prev<ofl_period_next) then !>>>>>>
        x1=ofl_readtime_prev+ofl_period_prev
        rap=1.-max( (x1-elapsedtime_now  )                   &
                   /(x1-ofl_readtime_prev) , 0.d00 )
       else                                     !>>>>>>
        x1=ofl_period_prev+2.*ofl_readtime_prev-ofl_readtime_next
        rap=1.-min( (ofl_readtime_next-elapsedtime_now)      &
                   /(ofl_readtime_next-x1             ),1.d00)
       endif                                    !>>>>>>

       x2=rap      ! weighted time interpolation "next" field
       x0=1.-x2    ! weighted time interpolation "previous" field

!      rap_obc=rap          ! bidouille pour que les interpolations temporelles fonctionnent
       timeweightobc(:)=rap ! bidouille pour que les interpolations temporelles fonctionnent !11-07-14
                            ! dans graph_out pour visualiser ssh-sshobc.....


! interpolation sshobc si cas nemo ou cas symphonie "anciens" fichers
      if(flag_nemoffline==1.or.flag_newalgo==0) then !------->
       do j=0,jmax+1 ; do i=0,imax+1 ! j=1,jmax !25-09-14 extension des boucles pour confort graphique
         sshobc_w(i,j,1)=(1.-timeweightobc(1))*sshobc_w(i,j,0) &
                            +timeweightobc(1) *sshobc_w(i,j,2)
!       if(i+par%timax(1)==280.and.j+par%tjmax(1)==90)write(90,*)elapsedtime_now/86400.,sshobc_w(i,j,1)
       enddo ; enddo
      endif                                          !------->

      if(iteration3d==0) then !****************>            !17-04-11

       if(flag_newalgo==1) then !11111> !08-05-19
        do j=0,jmax+1 ; do i=0,imax+1 
         ssh_int_w(i,j,:)=                                            &
         0.5*((1.-timeweightobc(1))*(sshobc_w(i,j,0)+sshobc_w(i,j,1)) &
                 +timeweightobc(1) *(sshobc_w(i,j,2)+sshobc_w(i,j,1)))
        enddo         ; enddo
       else                     !11111>
        do j=0,jmax+1 ; do i=0,imax+1 
         ssh_int_w(i,j,:)=sshobc_w(i,j,1)
        enddo         ; enddo
       endif                    !11111>
       do j=0,jmax+1 ! 1,jmax !29-04-14
       do i=0,imax+1 ! 1,imax
         hz_w(i,j,0)=h_w(i,j)+ssh_int_w(i,j,0)
         hz_w(i,j,1)=h_w(i,j)+ssh_int_w(i,j,1)
         hz_w(i,j,2)=h_w(i,j)+ssh_int_w(i,j,2)
         do k=1,kmax
          dz_t(i,j,k,0)=hz_w(i,j,0)*dsig_t(i,j,k)
          dz_t(i,j,k,1)=hz_w(i,j,1)*dsig_t(i,j,k)
          dz_t(i,j,k,2)=hz_w(i,j,2)*dsig_t(i,j,k)
         enddo
       enddo
       enddo
      endif                   !****************>

! Update T,S:
       do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1 !07-03-23
          tem_t(i,j,k,1)=(1.-timeweightobc(1))*temobc_t(i,j,k,0)    &
                            +timeweightobc(1) *temobc_t(i,j,k,2)
          sal_t(i,j,k,1)=(1.-timeweightobc(1))*salobc_t(i,j,k,0)    &
                            +timeweightobc(1) *salobc_t(i,j,k,2)
       enddo ; enddo ; enddo

!.......
! Update densite ?
       if(drifter_onoff==1.and.flag_buo_w==4) then !---> !07-03-23
! Si drifters avec vitesse de flottabilite calculee en ligne alors il faut calculer la densite de l'eau !07-03-23
       if(eos_comprs==0) then !000000>
! Compute rhp the potential density at z=0m and the a1, a2, ....,a7 coefficients of the Jackett et al 2006 (JAOT) EOS
         if(eos_author==0)call equation_of_state_linear(now)
         if(eos_author==3)call equation_of_state_rhp_a1_a7(now) ! rhp=densite potentielle (zref=0)
       endif                  !000000>
       if(eos_comprs==1) then !111111>
         if(eos_author==3)call equation_of_state_full_a1_a7(now) ! rhp=densite potentielle (zref=0)
       endif                  !111111>
       rho_t(0:imax+1,0:jmax+1,1:kmax)=rhp_t(0:imax+1,0:jmax+1,1:kmax)+rho !26-01-23
      endif                                       !--->

! Update U,V:
      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax+1

      veldydz_u(i,j,k, 1)=(1.-timeweightobc(1))*velobc_u(i,j,k,0)  & !16-02-10
                             +timeweightobc(1) *velobc_u(i,j,k,2)
      veldxdz_v(i,j,k, 1)=(1.-timeweightobc(1))*velobc_v(i,j,k,0)  & !16-02-10
                              +timeweightobc(1)*velobc_v(i,j,k,2)

      enddo ; enddo ; enddo

! Update Kz:
      sum1=0.
      do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
           kh_w(i,j,k)=(1.-timeweightobc(1))*dfvofl_w(i,j,k,-1)  &
                          +timeweightobc(1) *dfvofl_w(i,j,k,0 )
      enddo ; enddo ; enddo
! Update tken
      if(ofl_tke==1) then !tke-tke-> !05-03-15
       do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
           tken_w(i,j,k)=(1.-timeweightobc(1))*tkeofl_w(i,j,k,-1)   &
                            +timeweightobc(1) *tkeofl_w(i,j,k,0 )
       enddo ; enddo ; enddo
      endif               !tke-tke-> !05-03-15

! Nemo rivers:
      if(flag_nemoffline==1) then  !------->
       do kr=1,nriver
        riverflux(kr,1)=(1.-timeweightobc(1))*riverflux(kr,0)  &
                           +timeweightobc(1) *riverflux(kr,2) !02-04-13
       enddo
      endif                        !------->

! water flux through at the ocean surface:
      if(flag_nemoffline==1.or.flag_newalgo==1) then !wwwww>
       do j=1,jmax ; do i=1,imax
        omega_w(i,j,kmaxp1,1)=(1.-timeweightobc(1))*w0mofl_w(i,j,-1)  &
                                 +timeweightobc(1) *w0mofl_w(i,j,0 )
       enddo ; enddo
      endif                                          !wwwww>
! groundwater flux:
      if(flag_groundwater==1) then !m°v°m> !11-04-19
       do j=1,jmax ; do i=1,imax
        omega_w(i,j,1,1)=(1.-timeweightobc(1))*w_keq1_ofl_w(i,j,-1)  &
                            +timeweightobc(1) *w_keq1_ofl_w(i,j,0 )
       enddo ; enddo
      endif                        !m°v°m> !11-04-19

! Convective surface layer depth
      if(flag_ksloffline==1) then !m°v°m> !09-02-19
! https://docs.google.com/document/d/17xbkkh_KwMwHof8Q0iDcKnOcS2ZIt_ztgkW0EzswARo/edit
       do j=1,jmax ; do i=1,imax
        ksl_t(i,j)=nint( (1.-timeweightobc(1))*kslofl_t(i,j,-1)  &
                            +timeweightobc(1) *kslofl_t(i,j,0 )) !09-02-19
        do k=ksl_t(i,j),kmax !Couche "convective" on applique un fort Kz pour la simu bio offline 
         kh_w(i,j,k)=10.
        enddo
       enddo ; enddo
      endif                       !m°v°m> 

      if(flag_maxbotstress==1) then ! m[°v°]m > !17-11-17
       if(.not.allocated(maxbotstress_w))allocate(maxbotstress_w(imax,jmax))
       do j=1,jmax ; do i=1,imax
        maxbotstress_w(i,j)=(1.-timeweightobc(1))*maxbotstress_bef_w(i,j) &
                               +timeweightobc(1) *maxbotstress_aft_w(i,j)
       enddo ; enddo
      endif                         ! m[°v°]m > !17-11-17

! Calculer les composantes du transport:
!     do j=1,jmax+1
!     do i=1,imax+1
!      fluxbar_u(i,j,1)=0.
!      fluxbar_v(i,j,1)=0.
!     enddo
!     enddo

      fluxbar_u(:,:,1)=0. !03-03-16
      fluxbar_v(:,:,1)=0.

      do k=1,kmax
      do j=1,jmax !03-03-16
      do i=1,imax+1
       fluxbar_u(i,j,1)=fluxbar_u(i,j,1)+veldydz_u(i,j,k,1)  
      enddo
      enddo
      enddo

      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax !03-03-16
       fluxbar_v(i,j,1)=fluxbar_v(i,j,1)+veldxdz_v(i,j,k,1)
      enddo
      enddo
      enddo

! Elevation de la surface:
      do j=0,jmax+1 ! 1,jmax !25-09-14 Extension des boucles
      do i=0,imax+1 ! 1,imax !         pour confort graphique
! Bascule exceptionnellement placée ici pour ne pas avoir à
! modifier le call zthickness dans model_offline...
        ssh_int_w(i,j,0)=ssh_int_w(i,j,1)
        ssh_int_w(i,j,1)=ssh_int_w(i,j,2)
      enddo
      enddo

      do j=1,jmax
      do i=1,imax

       ! MOD Thom 
       !ssh_int_w(i,j,2)=max(wetdry_cst3-h_w(i,j) ,           & ! seuil zones intertidales !02-03-24
        ssh_int_w(i,j,2)=max(0.001*wetdry_cst3-h_w(i,j) ,           & ! seuil zones intertidales !23-07-14
       (ssh_int_w(i,j,0)-dti_lp*(                                   &
              ( fluxbar_u(i+1,j  ,1)-fluxbar_u(i,j,1)               &
               +fluxbar_v(i  ,j+1,1)-fluxbar_v(i,j,1) )/dxdy_t(i,j) &
               +omega_w(i,j,kmaxp1,1)                               &
               -omega_w(i,j,     1,1)                               & !11-04-19
                                ) )*mask_t(i,j,kmaxp1))  
! La multiplication par le masque sert a eviter l'abaissement de la ssh en amont
! des points source des rivieres

      enddo
      enddo
   
      if(nbcanal>0) then !m°v°m>
       call webcanals_gt2_to_gt1_sshint !23-04-19
       call webcanals_gt1_to_gt2_sshint !23-04-19
      endif              !m°v°m>

      do j=1,jmax
      do i=1,imax
! deplacE apres call webcanals_gt2_to_gt1_sshint le 06-05-19
        hz_w(i,j,2)=h_w(i,j)+ssh_int_w(i,j,2) !16-06-10
      enddo
      enddo

! VERIFICATION DU CALCUL DANS LE CAS SYMPHONIE
!     if(flag_newalgo==1) then !----->
!      do j=1,jmax ; do i=1,imax
!      i=30-par%timax(1)  ; j=3-par%tjmax(1)
!      if(i>0.and.i<imax+1.and.j>0.and.j<jmax+1)                 &
!       write(10+par%rank,*)                                     &
!       (1.-timeweightobc(1))*(sshobc_w(i,j,1)-sshobc_w(i,j,0))/ofl_period_prev &
!          +timeweightobc(1) *(sshobc_w(i,j,2)-sshobc_w(i,j,1))/ofl_period_next &
!      ,(ssh_int_w(i,j,2)-ssh_int_w(i,j,0))/dti_lp
!      enddo       ; enddo
!      write(10+par%rank,*)elapsedtime_now,ssh_int_w(i,j-1,2),ssh_int_w(i,j,2)
!     endif                    !----->

! CAS NEMO: AJUSTER LA DIVERGENCE DU TRANSPORT
      if(flag_nemoffline==1) then
!        call offline_another_model_divergence
         call offline_another_model_divergence_2 !01-07-16
      endif

!     do j=1,jmax ; do i=1,imax
!       hz_w(i,j,2)=h_w(i,j)+ssh_int_w(i,j,2) !01-07-16
!     enddo ; enddo

#ifdef bidon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Les lignes suivantes servent A verifier l'apport de la correction du
! toit libre
      sum1=0. ; sum2=0. ; sum3=0.
      sum4=0. ; sum5=0.
      do j=1,jmax
      do i=1,imax
       sum1=sum1+mask_t(i,j,kmax)
       sum2=sum2+mask_t(i,j,kmax)*(ssh_int_w(i,j,2)-sshobc_w(i,j,1))
       sum3=sum3+mask_t(i,j,kmax)*(ssh_int_w(i,j,2)-sshobc_w(i,j,1))**2
       sum4=sum4+mask_t(i,j,kmax)* ssh_int_w(i,j,2)
       sum5=sum5+mask_t(i,j,kmax)*  sshobc_w(i,j,1)
      enddo
      enddo
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum3,sum3glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum4,sum4glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum5,sum5glb,1,mpi_double_precision,  & !#MPI
                         mpi_sum,par%comm2d,ierr)
      if(par%rank==0) then                                         
        write(66,'(6(1x,e14.7))')elapsedtime_now/86400.            &
       ,sum2glb/sum1glb                                            &
       ,sqrt(sum3glb)/sum1glb                                      &
       ,sqrt( sum3glb/sum1glb - (sum2glb/sum1glb)**2 )             &
       ,sum4glb/sum1glb                                            &
       ,sum5glb/sum1glb             
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif

! Pour confort graphique une obc sur ssh_int sur couronne extra-peripherique
      call obc_ssh_int_clamped !25-09-14 ! "clamped" ssh_int=sshobc
!     call obc_ssh_int_0grd    !25-09-14 ! zero gradient sur ssh_int

! Update dz  ! place avant calcul vel=veldxdz/dy/dz le !23-07-14
!            ! pour maintenir coherence entre veldxdz et le nombre courant schema d'adv bio
      do k=1,kmax
      do j=0,jmax+1                    !21-10-10
      do i=0,imax+1                    !21-10-10
       dz_t(i,j,k,0)=dz_t(i,j,k,1)
       dz_t(i,j,k,1)=dz_t(i,j,k,2)
      enddo
      enddo
      enddo

! Placee ici cette operation permet de visualiser le courant en prenant compte de
! l'ajustement de la suface libre "nemo offline"
! D'autre part, placer le calcul de vel apres la mise a jour de dz(1), conserve
! la coherence du nombre de courant calcule dans le schema d'advection des traceurs
! bio, car la relation vel=veldxdz/dy/dz utilisera les memes valeurs pour dz.
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
! Parce que velbot "offline" a besoin de dz_u(dnow)
       dz_u(i,j,k,dnow)=0.5*(dz_t(i  ,j,k,1)      &
                            +dz_t(i-1,j,k,1))     !14-07-21
      enddo ; enddo ; enddo
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
           vel_u(i,j,k,1)=                                     &
       veldydz_u(i,j,k,1)/dy_u(i,j)/dz_u(i,j,k,dnow)
      enddo ; enddo ; enddo
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
! Parce que velbot "offline" a besoin de vel_u(0)
           vel_u(i,j,k,0)=vel_u(i,j,k,1) !14-07-21
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
! Parce que velbot "offline" a besoin de dz_v(dnow)
       dz_v(i,j,k,dnow)=0.5*(dz_t(i,j  ,k,1)      &
                            +dz_t(i,j-1,k,1))     !14-07-21
      enddo ; enddo ; enddo
      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
           vel_v(i,j,k,1)=                                     &
       veldxdz_v(i,j,k,1)/dx_v(i,j)/dz_v(i,j,k,dnow)
      enddo ; enddo ; enddo
      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
! Parce que velbot "offline" a besoin de vel_v(0)
           vel_v(i,j,k,0)=vel_v(i,j,k,1) !14-07-21
      enddo ; enddo ; enddo

      do k=1,kmax
      do j=1,jmax
      do i=1,imax
       dz_t(i,j,k,2)=hz_w(i,j,2)*dsig_t(i,j,k)
      enddo
      enddo
      enddo

      call obc_dz                      !21-10-10

! compute omega:   !16-06-10
      call omega   !27-09-10

!     i=191-par%timax(1)
!     j=200-par%tjmax(1)
!     if(i>1.and.i<imax.and.j>1.and.j<jmax) then
!       do k=kmax+1,1,-1
!         write(100,*)k,omega_w(i,j,k,1) &
!                            ,(1.-timeweightobc(1))*wobc_w(i,j,k,-1)  &
!                                +timeweightobc(1) *wobc_w(i,j,k,0 )
!       enddo
!     endif
!......................................................................
! Deduire Omega de la vitesse.
! fin.
!......................................................................




!......................................................................
! Cas forward uniquement sinon stop:
! Début:
      if(itimebio.eq.0)stop 'cas leapfrog pour bio non autorise!'
!......................................................................


      end subroutine offline_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine offline_netcdf_defvar
      implicit none
      include 'netcdf.inc'
#ifdef synopsis
       subroutinetitle='offline_netcdf_defvar'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(texte80(7)(1:4)=='real') then !!!!!!!!!!!!!!!>

      status=nf_def_var(ncid,texte80(1),nf_real,k0,vardim,varid(k10))
      if(status/=0)call offline_netcdf_error(1,1)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'units',len_trim(texte80(2)),texte80(2))
      if(status/=0)call offline_netcdf_error(1,2)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'long_name',len_trim(texte80(3)),texte80(3))
      if(status/=0)call offline_netcdf_error(1,3)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'standard_name',len_trim(texte80(4)),texte80(4))
      if(status/=0)call offline_netcdf_error(1,4)

      status=nf_put_att_text(ncid,varid(k10)                            &
!           ,'axis',len_trim(texte80(5)),texte80(5))
            ,'content',len_trim(texte80(5)),texte80(5))                !05-03-10
      if(status/=0)call offline_netcdf_error(1,5)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'associate',len_trim(texte80(6)),texte80(6))
      if(status/=0)call offline_netcdf_error(1,6)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'coordinates',len_trim(texte80(6)),texte80(6))
      if(status/=0)call offline_netcdf_error(1,6)

      status=nf_put_att_real(ncid,varid(k10)                            &
            ,'missing_value',nf_real,1,filval)
      if(status/=0)call offline_netcdf_error(1,10)

      status=nf_put_att_real(ncid,varid(k10)                            &
            ,'_FillValue',nf_real,1,filval)
      if(status/=0)call offline_netcdf_error(1,11)

      if(cgridshift(1)/=0.) then !ccccccccccccccccc>
      status=nf_put_att_real(ncid,varid(k10)                            &
            ,'C_grid_i_axis_shift',nf_real,1,cgridshift(1))
      if(status/=0)call netcdf_error(1,13)
      endif                      !ccccccccccccccccc>

      if(cgridshift(2)/=0.) then !ccccccccccccccccc>
      status=nf_put_att_real(ncid,varid(k10)                            &
            ,'C_grid_j_axis_shift',nf_real,1,cgridshift(2))
      if(status/=0)call netcdf_error(1,14)
      endif                      !ccccccccccccccccc>

      if(cgridshift(3)/=0.) then !ccccccccccccccccc>
      if(index(texte80(6),"depth")/=0) then !->->->
       status=nf_put_att_real(ncid,varid(k10)                            &
             ,'C_grid_k_axis_shift',nf_real,1,cgridshift(3))
       if(status/=0)call netcdf_error(1,15)
      endif                                 !->->->
      endif                      !ccccccccccccccccc>

      return
      endif                            !!!!!!!!!!!!!!!>


      if(texte80(7)(1:6)=='double') then !**************>

        if(texte80(1)(1:4)=='time') then  !§§§§§§§§§§§§§§§>

      status=nf_def_var(ncid,texte80(1),nf_double,k0,vardim,varid(k10)) !05-11-10
      if(status/=0)call offline_netcdf_error(1,1)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'units',len_trim(texte80(2)),texte80(2))
      if(status/=0)call offline_netcdf_error(1,2)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'time_origin',len_trim(texte80(8)),texte80(8))
      if(status/=0)call offline_netcdf_error(1,8)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'calendar',len_trim(texte80(9)),texte80(9))
      if(status/=0)call offline_netcdf_error(1,9)

!     status=nf_put_att_text(ncid,varid(k10)                            &!21-03-11
!           ,'axis',len_trim(texte80(5)),texte80(5))
!     if(status/=0)call offline_netcdf_error(1,5)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'associate',len_trim(texte80(6)),texte80(6))
      if(status/=0)call offline_netcdf_error(1,6)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'coordinates',len_trim(texte80(6)),texte80(6))
      if(status/=0)call offline_netcdf_error(1,6)

        return
        endif                             !§§§§§§§§§§§§§§§>

      return
      endif                              !**************>

      if(texte80(7)(1:7)=='integer') then !!!!!!!!!!!!!!!>

      status=nf_def_var(ncid,texte80(1),nf_int,k0,vardim,varid(k10))
      if(status/=0)call offline_netcdf_error(1,1)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'units',len_trim(texte80(2)),texte80(2))
      if(status/=0)call offline_netcdf_error(1,2)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'long_name',len_trim(texte80(3)),texte80(3))
      if(status/=0)call offline_netcdf_error(1,3)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'standard_name',len_trim(texte80(4)),texte80(4))
      if(status/=0)call offline_netcdf_error(1,4)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'axis',len_trim(texte80(5)),texte80(5))
      if(status/=0)call offline_netcdf_error(1,5)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'associate',len_trim(texte80(6)),texte80(6))
      if(status/=0)call offline_netcdf_error(1,6)

      status=nf_put_att_text(ncid,varid(k10)                            &
            ,'coordinates',len_trim(texte80(6)),texte80(6))
      if(status/=0)call offline_netcdf_error(1,6)

      if(cgridshift(1)/=0.) then !ccccccccccccccccc>
      status=nf_put_att_real(ncid,varid(k10)                            &
            ,'C_grid_i_axis_shift',nf_real,1,cgridshift(1))
      if(status/=0)call netcdf_error(1,13)
      endif                      !ccccccccccccccccc>

      if(cgridshift(2)/=0.) then !ccccccccccccccccc>
      status=nf_put_att_real(ncid,varid(k10)                            &
            ,'C_grid_j_axis_shift',nf_real,1,cgridshift(2))
      if(status/=0)call netcdf_error(1,14)
      endif                      !ccccccccccccccccc>

      if(cgridshift(3)/=0.) then !ccccccccccccccccc>
      if(index(texte80(6),"depth")/=0) then !->->->
       status=nf_put_att_real(ncid,varid(k10)                            &
             ,'C_grid_k_axis_shift',nf_real,1,cgridshift(3))
       if(status/=0)call netcdf_error(1,15)
      endif                                 !->->->
      endif                      !ccccccccccccccccc>

      return
      endif                            !!!!!!!!!!!!!!!>

      end subroutine offline_netcdf_defvar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine offline_netcdf_error(ki1,ki2)
      implicit none
      integer ki1,ki2
#ifdef synopsis
       subroutinetitle='offline_netcdf_error'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! diag netcdf errors in defining netcdf variables: begin
      if(ki1.eq.1) then

        if(par%rank==0)write(6,'(a,a)')'oups! Failed define netcdf :',TRIM(TEXTE80(1))
        if(ki2.eq.1)write(*,*)'in creating: varid'
        if(ki2.eq.2)write(*,*)'in creating: units'
        if(ki2.eq.3)write(*,*)'in creating: long_name'
        if(ki2.eq.4)write(*,*)'in creating: standard name'
        if(ki2.eq.5)write(*,*)'in creating: axis'
        if(ki2.eq.6)write(*,*)'in creating: associate or coordinates'
        if(ki2.eq.8)write(*,*)'in creating: time origin'
        if(ki2.eq.9)write(*,*)'in creating: calendar'
        if(ki2.eq.10)write(*,*)'in creating: missing_value'
        if(ki2.eq.11)write(*,*)'in creating: _FillValue'
        if(ki2.eq.13)write(*,*)'in creating: C_grid_Oi_axis_shift'
        if(ki2.eq.14)write(*,*)'in creating: C_grid_Oj_axis_shift'
        if(ki2.eq.15)write(*,*)'in creating: C_grid_Ok_axis_shift'


      stop ' stop in offline_inout.f'
      endif
! diag netcdf errors in defining netcdf variables: end



      end subroutine offline_netcdf_error

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine offline_endofaverage        !18-07-10
      implicit none
#ifdef synopsis
       subroutinetitle='offline_endofaverage'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


      do j=0,jmax+1   ! Bornes de boucles en adequation avec
      do i=0,imax+1   ! (plus loin) division de sshobc par sumarchive
       xy_t(i,j,1)=0.
      enddo
      enddo
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,1)=0.
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,1)=0.
      enddo
      enddo

      if(removetide==1) then !-rmtide-rmtide-rmtide->

      do 999 ktide=1,kmaxtide
! signal: F=F1*cos(freq*t)+F2*sin(freq*t)
! moyenne du signal sur 48h=F1*( sin(freq*t1)-sin(freq*t2) )/freq
!                          -F2*( cos(freq*t1)-cos(freq*t2) )/freq avec t1-t2=48h
!                          puis on divise par 48h

! Temps present:
      time1=frqtide(ktide)*                            &
        (elapsedtime_now-ti0tide(ktide))               &  !15-04-16
           +v0tide(ktide)+utide(ktide,1)
! Temps echeance precedente:
      time2=frqtide(ktide)*                            &
        (elapsedtime_now-ti0tide(ktide)-sumarchive)    &  !15-04-16
           +v0tide(ktide)+utide(ktide,1)

      x2=ftide(ktide,1)*(cos(time1)-cos(time2))/frqtide(ktide)        &
!                                              /(sumarchive*dti_fw)
                                               / sumarchive           !17-04-11
      x1=ftide(ktide,1)*(sin(time1)-sin(time2))/frqtide(ktide)        &
!                                              /(sumarchive*dti_fw)
                                               / sumarchive           !17-04-11

      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,1)=xy_u(i,j,1)+veltidecosout_u(i,j,ktide)*x1        &
                              -veltidesinout_u(i,j,ktide)*x2
      enddo
      enddo

      do j=1,jmax+1
      do i=1,imax
        xy_v(i,j,1)=xy_v(i,j,1)+veltidecosout_v(i,j,ktide)*x1        &
                               -veltidesinout_v(i,j,ktide)*x2
      enddo
      enddo

      do j=0,jmax+1 !18-12-14
      do i=0,imax+1
        xy_t(i,j,1)=xy_t(i,j,1)+sshtidecosout_w(i,j,ktide)*x1         &
                               -sshtidesinout_w(i,j,ktide)*x2
      enddo
      enddo

  999 continue

      endif                  !-rmtide-rmtide-rmtide->


      do j=0,jmax+1
      do i=0,imax+1
!      sshofl_w(i,j,1)=sshofl_w(i,j,1)/sumarchive     &
!                       -mask_t(i,j,kmax)*xy_t(i,j,1)  ! remove tide

       sshofl_w(i,j,1)=max(sshofl_w(i,j,1)/sumarchive     &
                            -mask_t(i,j,kmax)*xy_t(i,j,1) & ! remove tide
                          ,wetdry_cst3-h_w(i,j))            !06-02-15

       w0mofl_w(i,j,1)=w0mofl_w(i,j,1)/sumarchive    !04-07-15 

      enddo
      enddo
      if(flag_groundwater==1) then !pmx> !11-04-19
       do j=1,jmax ; do i=1,imax
         w_keq1_ofl_w(i,j,1)=w_keq1_ofl_w(i,j,1)/sumarchive
       enddo       ; enddo
      endif                        !pmx> !11-04-19
      if(flag_ksloffline==1) then !m°v°m>
       do j=1,jmax   ; do i=1,imax  
        kslofl_t(i,j,1)=kslofl_t(i,j,1)/sumarchive    !09-02-19
       enddo         ; enddo
      endif                       !m°v°m>
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
       temofl_t(i,j,k,1)=temofl_t(i,j,k,1)/sumarchive
       salofl_t(i,j,k,1)=salofl_t(i,j,k,1)/sumarchive
      enddo
      enddo
      enddo
      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1
       velofl_u(i,j,k,1)=velofl_u(i,j,k,1)/sumarchive  &
                          -mask_u(i,j,k)*xy_u(i,j,1)  ! remove tide
      enddo
      enddo
      enddo
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
       velofl_v(i,j,k,1)=velofl_v(i,j,k,1)/sumarchive  &
                          -mask_v(i,j,k)*xy_v(i,j,1)  ! remove tide
      enddo
      enddo
      enddo

      if(removetide==1) then !-rmtide-rmtide-rmtide->
       do j=1,jmax ; do i=1,imax+1 !20-10-20
         fluxbarsum_ofl_u(i,j)                            &
        =fluxbarsum_ofl_u(i,j)                            &
                  -mask_u(i,j,kmax)*xy_u(i,j,1)           & ! remove tide
              *(iteration2d_max_now*dy_u(i,j)*sumarchive) &
         *(h_u(i,j)+0.5*(sshofl_w(i,j,1)+sshofl_w(i-1,j,1)))
       enddo      ; enddo
       do j=1,jmax+1 ; do i=1,imax !20-10-20
         fluxbarsum_ofl_v(i,j)                            &
        =fluxbarsum_ofl_v(i,j)                            &
                  -mask_v(i,j,kmax)*xy_v(i,j,1)           & ! remove tide
              *(iteration2d_max_now*dx_v(i,j)*sumarchive) &
         *(h_v(i,j)+0.5*(sshofl_w(i,j,1)+sshofl_w(i,j-1,1)))
       enddo        ; enddo
      endif                  !-rmtide-rmtide-rmtide->

      if(ale_selected==1) then !-ale-ale-ale->
       do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
        dzofl_t(i,j,k,1)=mask_t(i,j,k)*dzofl_t(i,j,k,1)/sumarchive & !23-11-11
                    +(1.-mask_t(i,j,k))*1.e-6                        !23-11-11
       enddo ; enddo ; enddo
      endif                    !-ale-ale-ale->

      if(ofl_surflux==1) then !---> !01-08-14
      do j=0,jmax+1
      do i=0,imax+1
       slhf_aver_w(i,j)= slhf_aver_w(i,j)/sumarchive
       sshf_aver_w(i,j)= sshf_aver_w(i,j)/sumarchive
       snsf_aver_w(i,j)= snsf_aver_w(i,j)/sumarchive
       ssr_aver_w(i,j)=  ssr_aver_w(i,j)/sumarchive
       precipi_aver_w(i,j)=  precipi_aver_w(i,j)/sumarchive
       wstress_aver_u(i,j)= wstress_aver_u(i,j)/sumarchive
       wstress_aver_v(i,j)= wstress_aver_v(i,j)/sumarchive
      enddo
      enddo
      endif                   !--->

      if(ofl_bio==1) then      !bio-bio-bio-> !09-09-19
       do vb=1,vbmax ; do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
        bioofl_t(i,j,k,vb)=bioofl_t(i,j,k,vb)/sumarchive
       enddo ; enddo ; enddo ; enddo
       ! epaisseur de sediment ! 11_09_19
       do j=1,jmax ; do i=1,imax
        hsedofl_t(i,j)=hsedofl_t(i,j)/sumarchive
       enddo ; enddo
      endif                    !bio-bio-bio-> !09-09-19

      if(flag_maxbotstress==1) then ! m[°v°]m > !09-11-17
         maxbotstress_w=maxbotstress_w/sumarchive
         stresswave_w=stresswave_w/sumarchive    ! 01-02-21
         stressc_w=stressc_w/sumarchive
      endif                         ! m[°v°]m >

      end subroutine offline_endofaverage

!...................................................................................

      subroutine offline_readthelists(choice_   )      !30-03-11
      implicit none
      include 'netcdf.inc'
      integer choice_,nc_,loop_,loopmax_,ncid_ &
             ,year_,month_,day_,hour_,minute_,second_
      integer :: max_time_counter_,flag_,varnum_,flagbinreclist_ &
                ,linemax_,loop2_,loop1_
      double precision :: time1_=0.,time2_=0.
      double precision,dimension(:),allocatable :: time_counter_
      real,dimension(:),allocatable :: time_counter4_    !02-04-14
#ifdef synopsis
       subroutinetitle='offline_readthelists'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)
#endif

      flag_stop=0 !25-03-18

! listes binaires - binary lists

      cpu_seconds=MPI_Wtime ( )

      ofl_rec_now=-999  ; ofl_nextrec_time=-999.

      texte80(1)=trim(tmpdirname)//'liste_offline_u.binrec'
      texte80(2)=trim(tmpdirname)//'liste_offline_v.binrec'
      texte80(3)=trim(tmpdirname)//'liste_offline_t.binrec'
      texte80(4)=trim(tmpdirname)//'liste_offline_s.binrec'
      texte80(5)=trim(tmpdirname)//'liste_offline_ssh.binrec'
      texte80(6)=trim(tmpdirname)//'liste_offline_kz.binrec'
      texte80(7)=trim(tmpdirname)//'liste_offline_w.binrec'  !04-07-15

! Verifier que les listes ne sont pas deja presentes:
       if(par%rank==0) then !-rank0 checks existing list->
        open(unit=3,file=texte80(1)   &
            ,access='direct'          &
            ,recl=540                 &
            ,form='unformatted'       &
            ,status='new'             &
            ,iostat=flagbinreclist_   &
            )
        close(3)
       endif                !-rank0 checks existing list->
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)  !23-06-13
#endif
! Envoyer flagbinreclist_ A tous les autres proc:
      call mpi_bcast(flagbinreclist_,1,mpi_integer,0,par%comm2d,ierr) !19-04-16

      if(flag_nemoffline==0) then !------------->
! Cas standard symphonie
       loopmax_=1
       varnum_=7 !04-07-15
       if(ofl_tke==1) then !>>>>
        varnum_=varnum_+1
        texte80(varnum_)=trim(tmpdirname)//'liste_offline_tke.binrec' !05-03-15
       endif               !>>>>
      else                        !------------->
! Cas nemoffline:
       loopmax_=7                               !11-03-15
       varnum_=7
      endif                       !------------->

      do loop_=1,loopmax_

! Si flagbinreclist_ est nul cela signifie que la liste binaire n'existe
! pas et qu'il faut la creer. Dans le cas contraire cela est interprEtE comme le
! fait qu l'utilisateur a intentionnellement placE des listes binaires (d'un precedent run)
! dans le repertoire tmp dans le but de ne pas recalculer les listes binaires. Le modEle
! prend alors l'initiative de considerer que les listes existent dEjA et prend alors 
! l'aiguillage 2453 qui lui permet de se positionner dans la liste binaire existante.
! Sinon il continue comme d'hab.
       if(flagbinreclist_/=0) then !existing old list in tmp directory >
        call offline_oldlist(loop_)
        goto 2453
       endif          !existing old list in tmp directory >


! Faire la liste format binaire acces direct:
        nc_   =1
        open(unit=3,file=trim(offlinefile(min(loop_,loopmax_))))

! Combien y a t'il de lignes dans le fichier?
         linemax_=0
 2834    read(3,'(a)',end=2832)texte250 ; linemax_=linemax_+1 ; goto 2834
 2832    rewind(3)
         if(linemax_==0) then !>>>>> !23-04-14
          write(6,'(a,a,a)')'File ',trim(obcfile(1,1)),' is empty'
          stop ' stop ogcm_lire_les_listes erreur 2205'
         endif                !>>>>>

! Chaque proc va traiter une fraction du fichier. On commence par !25-05-19
! derouler les lignes qui ne concernent pas le proc
         do loop2_=1,int(real(par%rank  )/real(nbdom)*linemax_)
          read(3,*) ! lire ces lignes pour rien
         enddo
         write(texte30,'(a,i0)')trim(tmpdirname)//'tmpfile',par%rank
         open(unit=7,file=texte30,status='REPLACE')

        do loop2_=int(real(par%rank  )/real(nbdom)*linemax_)+1    &
                 ,int(real(par%rank+1)/real(nbdom)*linemax_)

!        if(par%rank==0)write(6,*)int(real(par%rank)/real(nbdom)*linemax_)+1,loop2_,int(real(par%rank+1)/real(nbdom)*linemax_)

 1131    read(3,'(a)',end=132)filename_
         if(filename_=='files created by the offline mode:')goto 1131



! Open the netcdf ogcm file and get the time_counter dimension then store in the binary file as much as time_counter:
! so that each line of the binary file corresponds to the ogcm output periodicity:
        !!print *, par%rank," K1=",K1," ",trim(texte250)
        max_time_counter_=1
        status=nf_open(trim(filename_),nf_nowrite,ncid_)
            if(status/=0) then !>>>>>>>>>>>
               write(6,'(i4,a,a)')par%rank,' ouverture ',trim(filename_) ; flag_stop=1 ; goto 2731
            endif              !>>>>>>>>>>>
                         status=nf_inq_dimid(ncid_,'time_counter',dim_t_id)
            if(status/=0)status=nf_inq_dimid(ncid_,'TIME_COUNTER',dim_t_id)
            if(status/=0)status=nf_inq_dimid(ncid_,'time',dim_t_id)
            if(status/=0)status=nf_inq_dimid(ncid_,'TIME',dim_t_id)
            if(status==0)status=nf_inq_dimlen(ncid_,dim_t_id,max_time_counter_)
! float time_counter(time_counter) ;
! time_counter:calendar = "gregorian" ;
! time_counter:standard_name = "time" ;
! time_counter:long_name = "Time axis" ;
! time_counter:axis = "T" ;
! time_counter:units = "seconds since 2006-10-11 00:00:00" ;
! time_counter:time_origin = "2006-OCT-11 00:00:00" ;

                         status=nf_inq_varid(ncid_,'time',var_id)
            if(status/=0)status=nf_inq_varid(ncid_,'time_counter',var_id)
            if(status/=0)stop 'erreur nf_inq_varid time module_offline'

            txt_units=''
            status=nf_get_att_text(ncid_,var_id,'units',txt_units);if(status/=0)stop 'erreur nf_get_att_text module_offline'
            k=index(txt_units,'since')
            read(txt_units(k+6:k+9),*)year_
            read(txt_units(k+11:k+12),*)month_
            read(txt_units(k+14:k+15),*)day_
            read(txt_units(k+17:k+18),*)hour_
            read(txt_units(k+20:k+21),*)minute_
            read(txt_units(k+23:k+24),*)second_
            call datetokount(year_,month_,day_,hour_,minute_,second_) ! donne elapsedtime_out

            status=nf_inq_var(ncid_,var_id,texte30,var_type,var_dims,tabdim,i4)

            if(.not.allocated(time_counter_)) then
              allocate(time_counter_(max_time_counter_))
              if(var_type==5)allocate(time_counter4_(max_time_counter_))
            endif
            if(ubound(time_counter_,dim=1)<max_time_counter_) then !...>
             deallocate(time_counter_)
               allocate(time_counter_(max_time_counter_))
               if(var_type==5) then !rrrrrrrrr>
                  deallocate(time_counter4_)
                    allocate(time_counter4_(max_time_counter_))
               endif                !rrrrrrrrr>
            endif                                                  !...>

               if(var_type==6)status=nf_get_var_double(ncid_,var_id,time_counter_(1:max_time_counter_))
               if(var_type==5) then !rrrrrrrrr>
                 status=nf_get_var_real(ncid_,var_id,time_counter4_(1:max_time_counter_))
                 time_counter_(1:max_time_counter_)=time_counter4_(1:max_time_counter_)
               endif                !rrrrrrrrr>
            if(status/=0)stop 'stop nf_get_var_double time_counter_'

            do ogcmtimecounter_=1,max_time_counter_
!            status=nf_get_vara_double(ncid_,var_id,ogcmtimecounter_,1,x1); &
!               if(status/=0)stop 'erreur1 nf_get_vara_double module_offline'
             x1=time_counter_(ogcmtimecounter_)
             x2=-999.
             if(index(txt_units,'days')/=0)x2=86400.
             if(index(txt_units,'hours')/=0)x2=3600.
             if(index(txt_units,'seconds')/=0)x2=1.
             if(x2==-999.)stop 'offline unites time incorrectes'
             time1_=time2_
             time2_=elapsedtime_out+x1*x2

            status=nf_inq_varid(ncid_,'cumulativetime',var_id)
            if(status==0) then !xxxxxxx>
               status=nf_get_vara_double(ncid_,var_id,ogcmtimecounter_,1,x3); &
                    if(status/=0)stop 'erreur2 nf_get_vara_double module_offline'
               status=nf_get_att_text(ncid_,var_id,'units',txt_units); &
                    if(status/=0)stop 'erreur nf_get_att_text module_offline'
             if(index(txt_units,'hours')==0) &
             stop 'Unites cumulativetime incorrectes'
             ofl_period_next=x3*3600.
             ofl_readtime_next=time2_-0.5*ofl_period_next
            endif              !xxxxxxx>

!           if(index(nomfichier(20),'another_model')/=0) then !nnnn>
            if(flag_nemoffline==1) then                       !nnnn>
             ofl_period_next=time2_-time1_ ! time1_ n'est pas connu au premier coup, donc
             ofl_readtime_next=time2_      ! ofl_period_next defini a posteriori quand nc_=2
! Verifier que les fichiers sont bien dans l'ordre chronologique:
             if(nc_>1.and.time1_>=time2_) then !bugbug>
              write(6,*)'Erreur chronologie des fichiers ogcm' &
                       ,' details dans fichiers fort'
              write(10+par%rank,*)'nc_   ',nc_
              write(10+par%rank,*)'time1_',time1_
              write(10+par%rank,*)'time2_',time2_
              write(10+par%rank,'(a,a)')'File: ',trim(filename_)
              stop 'STOP in subroutine offline_readthelists'
             endif                             !bugbug>
            endif                                             !nnnn>

!              write(4,rec=nc_)filename_           &
!                             ,ogcmtimecounter_    &
!                             ,ofl_readtime_next   &
!                             ,ofl_period_next
!                    if(nc_==2) then !222222>
!                       read(4,rec=1)filename_,i1,x1,x2
!                      write(4,rec=1)filename_,i1,x1,ofl_period_next,datesim(1:6,1) !09-02-17
!                    endif           !222222>
               write(7,'(a)')trim(filename_)
               write(7,*)nc_,' nc_'
               write(7,*)ogcmtimecounter_,' ogcmtimecounter_'
               write(7,*)ofl_readtime_next,' ofl_readtime_next'
               write(7,*)ofl_period_next,' ofl_period_next'

             ofl_rec_max=nc_
!            if(ofl_readtime_next<=elapsedtime_now) then !----->
!              ofl_rec_now=nc_
!              ofl_period_now=ofl_period_next
!            endif                                       !----->

             nc_=nc_+1
            enddo ! fin de boucle sur ogcmtimecounter_

         status=nf_close(ncid_)

        enddo ! loop2

  132   close(7) ; close(3)

!     write(6,*)'C ',par%rank,loop_
!      if(.not.allocated(time_counter_)) then !--->
!       write(6,'(a,a,a)')'LISTE '   &
!                        ,trim(offlinefile(min(loop_,loopmax_))) &
!                        ,' VIDE!!!!'
!       stop ' STOP subroutine offline_readthelists'
!      else                                   !--->
!       deallocate(time_counter_)
!      endif                                  !--->
       if(allocated(time_counter_)) deallocate(time_counter_) !06-06-19
       if(allocated(time_counter4_))deallocate(time_counter4_)

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)
#endif
! Le rank 0 assemble les listes intermediaires
      if(par%rank==0) then !00000000000000000000000000>
       open(unit=4,file=texte80(loop_) &       !22-12-09
                  ,access='direct'     &
                  ,recl=recordlength_  &
                  ,form='unformatted')
        nc_=1
        do loop1_=0,nbdom-1
         write(texte30,'(a,i0)')trim(tmpdirname)//'tmpfile',loop1_
         open(unit=7,file=texte30)
 1692         read(7,'(a)',end=1677)filename_
              read(7,*)i0
              read(7,*)ogcmtimecounter_
              read(7,*)ofl_readtime_next
              read(7,*)ofl_period_next

             if(ofl_readtime_next<=elapsedtime_now) then !----->
               ofl_rec_now=nc_
               ofl_period_now=ofl_period_next
             endif                                       !----->

               write(4,rec=nc_)filename_           &
                              ,ogcmtimecounter_    &
                              ,ofl_readtime_next   &
                              ,ofl_period_next
                     if(nc_==2) then !222222>
                        read(4,rec=1)filename_,i1,x1,x2
                       write(4,rec=1)filename_,i1,x1,ofl_period_next,datesim(1:6,1) !09-02-17
                     endif           !222222>
               ofl_rec_max=nc_
               nc_=nc_+1
               goto 1692
 1677   close(7)
        enddo !loop1_
       close(4)

! Dans le cas SYMPHONIE seule la premiere liste a ete faite car toutes
! les autres sont identiques, le rank0 se charge de faire les restantes
! en recopiant la liste 1
      if(flag_nemoffline==0) then !------------->
         do i=2,varnum_
           open(unit=4,file=texte80(1),access='direct',recl=recordlength_,form='unformatted')
           open(unit=7,file=texte80(i),access='direct',recl=recordlength_,form='unformatted')
           j=1
           read(4,rec=j)filename_,ogcmtimecounter_,ofl_readtime_next,ofl_period_next
          write(7,rec=j)filename_,ogcmtimecounter_,ofl_readtime_next,ofl_period_next,datesim(1:6,1)
           do j=2,ofl_rec_max
             read(4,rec=j)filename_,ogcmtimecounter_,ofl_readtime_next,ofl_period_next
            write(7,rec=j)filename_,ogcmtimecounter_,ofl_readtime_next,ofl_period_next
           enddo       
           close(4) ; close(7)
         enddo ! boucle k de 2 A varnum_
      endif                       !------------->
      

! LIGNES TEST (A COMMENTER) POUR VERIFIER LE CONTENU DES LISTES
!      open(unit=4,file=texte80(loop_) &       !22-12-09
!                 ,access='direct'     &
!                 ,recl=recordlength_  &
!                 ,form='unformatted')
!       
!       write(texte30,'(a,i0)')'veriflistebinaire',loop_
!       open(unit=7,file=texte30)
!       do k=1,nc_-1
!              read(4,rec=k)filename_           &
!                          ,ogcmtimecounter_    &
!                          ,ofl_readtime_next   &
!                          ,ofl_period_next
!              write(7,'(a)')trim(filename_)
!              write(7,*)ogcmtimecounter_
!              write(7,*)ofl_readtime_next
!              write(7,*)ofl_period_next
!       enddo
!      close(7)
!      close(4)

      endif                !00000000000000000000000000>
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)
#endif


 2453 continue

      enddo ! fin de boucle sur loop_ sur les listes ascii
 2731 call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0) stop 'ERREUR OFFLINE LISTES BINAIRES' !09-05-18
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)
#endif

!     stop 'KIKO'

! Le proc 0 va envoyer ses valeurs aux autres proc.
      call mpi_bcast(ofl_rec_max      ,1,mpi_integer         ,0,par%comm2d,ierr)
      call mpi_bcast(ofl_readtime_next,1,mpi_double_precision,0,par%comm2d,ierr)
      call mpi_bcast(ofl_rec_now      ,1,mpi_integer         ,0,par%comm2d,ierr)
      call mpi_bcast(ofl_period_now   ,1,mpi_double_precision,0,par%comm2d,ierr)
      call mpi_bcast(ofl_period_next  ,1,mpi_double_precision,0,par%comm2d,ierr)

! WARNING:
      flag_stop=0
      if(ofl_rec_now<0.or.ofl_rec_now==ofl_rec_max) then !-warning->
       flag_stop=1
       if(par%rank==0) then !wwww>
        write(10+par%rank,*)'Pas de fichier offline dispo etat initial'
        write(10+par%rank,*)'ofl_rec_now=',ofl_rec_now
        write(10+par%rank,*)'ofl_rec_max=',ofl_rec_max
       endif                !www>
      endif                                              !-warning->
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0) stop 'ERREUR OFFLINE LISTES BINAIRES see fortxxx files'

      cpu_seconds=MPI_Wtime ( ) - cpu_seconds
      if(par%rank==0)write(6,*)'cpu_seconds offline_readthelists =',cpu_seconds

      if(ofl_reversedtime==1)call offline_reversethelists(varnum_) !26-05-19

      end subroutine offline_readthelists

!...................................................................................

      subroutine offline_allocate(txt_)
      implicit none
      character*1 txt_
#ifdef synopsis
       subroutinetitle='offline_allocate'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(txt_=='a')  then !>>>>
         if(ioffline==2) then
          allocate(w0mofl_w(0:imax+1,0:jmax+1,-1:0)) 
         else
          allocate(w0mofl_w(0:imax+1,0:jmax+1,1)) 
         endif
         w0mofl_w=0.
      endif               !>>>>
      if(txt_=='d')deallocate(w0mofl_w)

      if(txt_=='a')  then !>>>>
         if(ioffline==2) then
!         allocate(w_keq1_ofl_w(imax,jmax,-1:0)) 
         else
          if(flag_groundwater==1)allocate(w_keq1_ofl_w(imax,jmax,1))  !10-04-19
         endif
       if(allocated(w_keq1_ofl_w))w_keq1_ofl_w=0.
      endif               !>>>>
      if(txt_=='d') then  !--->
        if(allocated(w_keq1_ofl_w))deallocate(w_keq1_ofl_w) !10-04-19
      endif               !--->


      if(flag_ksloffline==1) then !m°v°m>
      if(txt_=='a')  then !>>>>
         if(ioffline==2) then
          allocate(kslofl_t(1:imax  ,1:jmax  ,-1:0))  !09-02-19
         else
          allocate(kslofl_t(1:imax  ,1:jmax  ,1)) !09-02-19
         endif
         kslofl_t=real(kmax) !09-02-19
      endif               !>>>>
      if(txt_=='d')deallocate(kslofl_t) !09-02-19
      endif                       !m°v°m>

      if(txt_=='a')  allocate(temofl_t        (0:imax+1,0:jmax+1,kmax,1)) ; if(txt_=='a') temofl_t=0
      if(txt_=='d')deallocate(temofl_t)

      if(txt_=='a')  allocate(salofl_t        (0:imax+1,0:jmax+1,kmax,1)) ; if(txt_=='a') salofl_t=0
      if(txt_=='d')deallocate(salofl_t)

      if(txt_=='a')  allocate(velofl_u        (0:imax+1,0:jmax+1,kmax,1)) ; if(txt_=='a') velofl_u=0
      if(txt_=='d')deallocate(velofl_u)

      if(txt_=='a')  allocate(velofl_v        (0:imax+1,0:jmax+1,kmax,1)) ; if(txt_=='a') velofl_v=0
      if(txt_=='d')deallocate(velofl_v)

      if(txt_=='a')  allocate(dfvofl_w        (0:imax+1,0:jmax+1,kmax+1,-1:0)) ; if(txt_=='a') dfvofl_w=0
      if(txt_=='d')deallocate(dfvofl_w)

      if(txt_=='a')  allocate(sshofl_w        (0:imax+1,0:jmax+1,1))      ; if(txt_=='a') sshofl_w=0
      if(txt_=='d')deallocate(sshofl_w)

      if(txt_=='a')  allocate(velbarofl_u     (0:imax+1,0:jmax+1,1))      ; if(txt_=='a') velbarofl_u=0
      if(txt_=='d')deallocate(velbarofl_u)

      if(txt_=='a')  allocate(velbarofl_v     (0:imax+1,0:jmax+1,1))      ; if(txt_=='a') velbarofl_v=0
      if(txt_=='d')deallocate(velbarofl_v)

      if(txt_=='a')  allocate(fluxbarsum_ofl_u(imax+1,jmax+1))            ; if(txt_=='a') fluxbarsum_ofl_u=0
      if(txt_=='d')deallocate(fluxbarsum_ofl_u)

      if(txt_=='a')  allocate(fluxbarsum_ofl_v(imax+1,jmax+1))            ; if(txt_=='a') fluxbarsum_ofl_v=0
      if(txt_=='d')deallocate(fluxbarsum_ofl_v)

      if(ale_selected==1) then !---->
       if(txt_=='a')  allocate(dzofl_t(0:imax+1,0:jmax+1,kmax,1)) ; if(txt_=='a') dzofl_t=0
       if(txt_=='d')deallocate(dzofl_t)

      endif                    !---->
      if(ofl_surflux==1)then   ! ***
       if(txt_=='a')  allocate(slhf_aver_w   (0:imax+1,0:jmax+1)) ; if(txt_=='a') slhf_aver_w=0
       if(txt_=='d')deallocate(slhf_aver_w)

       if(txt_=='a')  allocate(sshf_aver_w   (0:imax+1,0:jmax+1)) ; if(txt_=='a') sshf_aver_w=0
       if(txt_=='d')deallocate(sshf_aver_w)

       if(txt_=='a')  allocate(snsf_aver_w   (0:imax+1,0:jmax+1)) ; if(txt_=='a') snsf_aver_w=0
       if(txt_=='d')deallocate(snsf_aver_w)

       if(txt_=='a')  allocate(ssr_aver_w    (0:imax+1,0:jmax+1)) ; if(txt_=='a') ssr_aver_w=0
       if(txt_=='d')deallocate(ssr_aver_w)

       if(txt_=='a')  allocate(precipi_aver_w(0:imax+1,0:jmax+1)) ; if(txt_=='a') precipi_aver_w=0
       if(txt_=='d')deallocate(precipi_aver_w)

       if(txt_=='a')  allocate(wstress_aver_u(0:imax+1,0:jmax+1)) ; if(txt_=='a') wstress_aver_u=0
       if(txt_=='d')deallocate(wstress_aver_u)

       if(txt_=='a')  allocate(wstress_aver_v(0:imax+1,0:jmax+1)) ; if(txt_=='a') wstress_aver_v=0
       if(txt_=='d')deallocate(wstress_aver_v)

      endif                    ! ***

      if(ofl_tke==1) then !05-03-15
       if(txt_=='a')  allocate(tkeofl_w(0:imax+1,0:jmax+1,kmax+1,-1:0)) ; if(txt_=='a') tkeofl_w=0
       if(txt_=='d')deallocate(tkeofl_w)
      endif               !19-03-15
      if(ofl_bio==1) then !09-09-19 et 11-09-19
       if(txt_=='a')  allocate(hsedofl_t(0:imax+1,0:jmax+1)) ; if(txt_=='a') hsedofl_t=0
       if(txt_=='a')  allocate(bioofl_t(0:imax+1,0:jmax+1,kmax+1,vbmax)) ; if(txt_=='a') bioofl_t=0
       if(txt_=='d')deallocate(bioofl_t)
       if(txt_=='d')deallocate(hsedofl_t)
      endif              

      if(txt_=='a')return
      if(txt_=='d')return
      stop 'offline_allocate Err on txt_'

      end subroutine offline_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine offline_write_grid(case_)
      implicit none
      double precision dlon_di_,dlon_dj_
      integer looproc_,looprocmax_,case_
      character txt100_*100,txt_type_*6
#ifdef synopsis
       subroutinetitle='offline_write_grid'
       subroutinedescription='Produces a netcdf grid file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!     if(flag_offline_binary==1) return ! pas de fichier grille si format binaire !08-06-18

! On n'ecrase pas le fichier de grille deja existant quand ioffline=2 !20-12-14
!     if(ioffline==2)return !20-12-14
      if(ioffline==2.and.case_==2)return !05-10-19

      offline_gridbin_status=0

! ECRIRE LE FICHIER DE GRILLE:


      if(flag_offline_binary==0)txt_type_='double'
      if(flag_offline_binary==1)txt_type_='real'

      do loop_netcdf=flag_offline_binary,1 ! do loop_netcdf=0,1 si netcdf, do loop_netcdf=1,1 si binaire

      count_netcdfvar=0

      filval=-9999.
      texte80(3)='none'
      texte80(4)='none'
      texte80(11)='none' !07-03-17

! Nom du fichier:
! Choix entre tmp/grid.nc et offlinefile(1)(1:k)//'grid.nc'
      if(case_==2) then !------->
       txt100_=trim(directory_offline)//'grid.nc' !21-07-14
      endif             !------->
      if(trim(offlinefile(1))=='s'.or.case_==1)txt100_=trim(tmpdirname)//'grid.nc'

      if(flag_offline_binary==0) then !netcdf-case->

      if(loop_netcdf==0) then  !§§§§§§§>
!      status=nfmpi_create(par%comm2d,txt100_ ,nf_clobber, MPI_INFO_NULL,ncid)
       status=nfmpi_create(par%comm2d,txt100_ ,nf_clobber + NF_64BIT_OFFSET, MPI_INFO_NULL,ncid) !10-09-12
      else                     !§§§§§§§>
!      status=nfmpi_open(par%comm2d,txt100_ ,nf_write, MPI_INFO_NULL,ncid)
       status=nfmpi_open(par%comm2d,txt100_ ,nf_write + NF_64BIT_OFFSET, MPI_INFO_NULL,ncid) !10-09-12
      endif                    !§§§§§§§>
      if(status/=0) then
        write(6,'(a,a)')'Error nfmpi_create ',trim(txt100_)
        stop ' stop module_offline'
      endif

      if(loop_netcdf==0)call netcdf_dim

! Reset VARDIM
      vardim(1)=i_t_dim ; vardim(2)=j_t_dim
      vardim(3)=k_t_dim ; vardim(4)=time_dim

! write variables dimensions:
      !call graph_out_trueaxis_nfmpi

      endif                           !netcdf-case->

      if(loop_netcdf==1) then !--------->
       if(txt_type_=='real')anyvar2d(0:imax+1,0:jmax+1)    =lon_t(0:imax+1,0:jmax+1)*rad2deg
       if(txt_type_=='double')anyv3d(0:imax+1,0:jmax+1,1,1)=lon_t(0:imax+1,0:jmax+1)*rad2deg
      endif                  !--------->
      texte80(1)='longitude_t'
      texte80(2)='degrees_east'                              ! units
      texte80(3)='longitude' ; texte80(4)='longitude'       ! long_name ; standard_name
      texte80(5)='YX' ; texte80(7)=txt_type_
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !08-06-18

      if(loop_netcdf==1) then !--------->
       if(txt_type_=='real')anyvar2d(0:imax+1,0:jmax+1)    =lon_u(0:imax+1,0:jmax+1)*rad2deg
       if(txt_type_=='double')anyv3d(0:imax+1,0:jmax+1,1,1)=lon_u(0:imax+1,0:jmax+1)*rad2deg
      endif                  !--------->
      texte80(1)='longitude_u'
      texte80(2)='degrees_east'                              ! units
      texte80(3:4)='longitude_at_u_location'                 ! long_name ; standard_name
      texte80(5)='YX' ;  texte80(7)=txt_type_
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') !08-06-18

      if(loop_netcdf==1) then !--------->
       if(txt_type_=='real')anyvar2d(0:imax+1,0:jmax+1)    =lon_v(0:imax+1,0:jmax+1)*rad2deg
       if(txt_type_=='double')anyv3d(0:imax+1,0:jmax+1,1,1)=lon_v(0:imax+1,0:jmax+1)*rad2deg
      endif                  !--------->
      texte80(1)='longitude_v'
      texte80(2)='degrees_east'                              ! units
      texte80(3:4)='longitude_at_v_location'                 ! long_name ; standard_name
      texte80(5)='YX' ; texte80(7)=txt_type_
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offlinegbin('_v') !08-06-18

      if(loop_netcdf==1) then !--------->
       if(txt_type_=='real')anyvar2d(0:imax+1,0:jmax+1)    =lon_f(0:imax+1,0:jmax+1)*rad2deg
       if(txt_type_=='double')anyv3d(0:imax+1,0:jmax+1,1,1)=lon_f(0:imax+1,0:jmax+1)*rad2deg
      endif                  !--------->
      texte80(1)='longitude_f'
      texte80(2)='degrees_east'                              ! units
      texte80(3:4)='longitude_at_f_location'                 ! long_name ; standard_name
      texte80(5)='YX' ; texte80(7)=txt_type_
      if(flag_offline_binary==0)call netcdf_main('_f')
      if(flag_offline_binary==1)call offlinegbin('_f') !08-06-18

      if(loop_netcdf==1) then !--------->
       if(txt_type_=='real')anyvar2d(0:imax+1,0:jmax+1)    =lat_t(0:imax+1,0:jmax+1)*rad2deg
       if(txt_type_=='double')anyv3d(0:imax+1,0:jmax+1,1,1)=lat_t(0:imax+1,0:jmax+1)*rad2deg
      endif                  !--------->
      texte80(1)='latitude_t'
      texte80(2)='degrees_north'                              ! units
      texte80(3)='latitude' ; texte80(4)='latitude'       ! latg_name ; standard_name
      texte80(5)='YX' ; texte80(7)=txt_type_
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !08-06-18


      if(loop_netcdf==1) then !--------->
       if(txt_type_=='real')anyvar2d(0:imax+1,0:jmax+1)    =lat_u(0:imax+1,0:jmax+1)*rad2deg
       if(txt_type_=='double')anyv3d(0:imax+1,0:jmax+1,1,1)=lat_u(0:imax+1,0:jmax+1)*rad2deg
      endif                  !--------->
      texte80(1)='latitude_u'
      texte80(2)='degrees_north'                              ! units
      texte80(3:4)='latitude_at_u_location'                   ! latg_name ; standard_name
      texte80(5)='YX' ; texte80(7)=txt_type_
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') !08-06-18


      if(loop_netcdf==1) then !--------->
       if(txt_type_=='real')anyvar2d(0:imax+1,0:jmax+1)    =lat_v(0:imax+1,0:jmax+1)*rad2deg
       if(txt_type_=='double')anyv3d(0:imax+1,0:jmax+1,1,1)=lat_v(0:imax+1,0:jmax+1)*rad2deg
      endif                  !--------->
      texte80(1)='latitude_v'
      texte80(2)='degrees_north'                              ! units
      texte80(3:4)='latitude_at_v_location'                   ! latg_name ; standard_name
      texte80(5)='YX' ; texte80(7)=txt_type_
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offlinegbin('_v') !08-06-18

      if(loop_netcdf==1) then !--------->
       if(txt_type_=='real')anyvar2d(0:imax+1,0:jmax+1)    =lat_f(0:imax+1,0:jmax+1)*rad2deg
       if(txt_type_=='double')anyv3d(0:imax+1,0:jmax+1,1,1)=lat_f(0:imax+1,0:jmax+1)*rad2deg
      endif                  !--------->
      texte80(1)='latitude_f'
      texte80(2)='degrees_north'                              ! units
      texte80(3:4)='latitude_at_f_location'                   ! latg_name ; standard_name
      texte80(5)='YX' ; texte80(7)=txt_type_
      if(flag_offline_binary==0)call netcdf_main('_f')
      if(flag_offline_binary==1)call offlinegbin('_f') !08-06-18

      if(loop_netcdf==1) then !--------->  !22-10-14
      do j=0,jmax+1 ; do i=0,imax+1
       anyv3d(i,j,1,1)=(i*1000.+j)/1000.
      enddo         ; enddo
      endif                  !--------->
      texte80(1)='i.j'
      texte80(2)='i.j'                               ! units
      texte80(3:4)='local_indexes'                   ! latg_name ; standard_name
      texte80(5)='YX' ; texte80(7)=txt_type_
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !08-06-18

      if(loop_netcdf==1) then !--------->
       anyvar2d(0:imax+1,0:jmax+1)=coriolis_t(0:imax+1,0:jmax+1)
      endif                  !--------->
      texte80(1)='coriolis_t'
      texte80(2)='s-1'
      texte80(3:4)='coriolis_parameter_at_t_location'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !08-06-18

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; do i=0,imax+1
        anyvar2d(i,j)=mask_t(i,j,kmax)*z0_w(i,j)+(1-mask_t(i,j,kmax))*filval !18-02-20
       enddo         ; enddo
      endif                  !--------->
      texte80(1)='z0_w'
      texte80(2)='m'
      texte80(3:4)='Bottom_roughness_length'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') 

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then
            anyvar2d(i,j)=h_w(i,j) ! ecrire en real
          anyv3d(i,j,1,1)=h_w(i,j) ! ecrire en double
       else
            anyvar2d(i,j)=-9999.
          anyv3d(i,j,1,1)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='hm_w' ; texte80(2)='m'                     ! variable ; units
      texte80(3:4)='model_sea_floor_depth_below_geoid'
!     texte80(5)='YX' ; texte80(7)='real'
      texte80(5)='YX' ; texte80(7)=txt_type_ !01-05-14
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offlinegbin('_w') !18-06-18

      if(loop_netcdf==1) then !---------> !08-04-16
       do j=0,jmax+1 ; do i=0,imax+1
            anyvar2d(i,j)=h_f(i,j) ! ecrire en real
          anyv3d(i,j,1,1)=h_f(i,j) ! ecrire en double
       enddo         ; enddo
      endif                  !--------->
      texte80(1)='h_f' ; texte80(2)='m'                     ! variable ; units
      texte80(3:4)='model_sea_floor_depth_below_geoid_at_f_point'
!     texte80(5)='YX' ; texte80(7)='real'
      texte80(5)='YX' ; texte80(7)=txt_type_
      if(flag_offline_binary==0)call netcdf_main('_f')
      if(flag_offline_binary==1)call offlinegbin('_f') !18-06-18

      if(loop_netcdf==1) then !---------> !07-11-15
       write(texte30,'(a,i0)')trim(tmpdirname)//'err_h',par%rank !17-07-20
       open(unit=3,file=trim(texte30)) !07-11-15
        xy_t(:,:,0)=0. !05-10-19
        do j=0,jmax+1 ; do i=0,imax+1
!        read(3,*)xy_t(i,j,0)
         read(3,*,end=2070)xy_t(i,j,0) !05-10-19
        enddo         ; enddo
!      close(3)
 2070  close(3)
       do j=0,jmax+1 ; do i=0,imax+1
        anyvar2d(i,j)=    xy_t(i,j,0)*mask_t(i,j,kmax)-(1.-mask_t(i,j,kmax))*9999.
          anyv3d(i,j,1,1)=xy_t(i,j,0)*mask_t(i,j,kmax)-(1.-mask_t(i,j,kmax))*9999.
       enddo ; enddo
      endif                  !--------->
      texte80(1)='h_err' ; texte80(2)='m'                     ! variable ; units
      texte80(3:4)='smoothed_h_minus_h' ; texte80(11)='mask_obc_z1' !14-01-17
      texte80(5)='YX' ; texte80(7)=txt_type_ !01-05-14
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offlinegbin('_w') !18-06-18

      if(mergebathy_filename/='nomerge') then !pmpmpmpmpm> !18-11-15
       if(loop_netcdf==1) then !---------> 
!       write(texte30,'(a,i0)')'tmp/h-hnemo-',par%rank
        write(texte30,'(a,i0)')'tmp/hnemo-',par%rank !05-01-17
        open(unit=3,file=trim(texte30)) !18-11-15
         do j=1,jmax ; do i=1,imax
          read(3,*)xy_t(i,j,0)
         enddo         ; enddo
        close(3)
        do j=0,jmax+1 ; do i=0,imax+1
         anyvar2d(i,j)=    xy_t(i,j,0)*mask_t(i,j,kmax)-(1.-mask_t(i,j,kmax))*9999.
           anyv3d(i,j,1,1)=xy_t(i,j,0)*mask_t(i,j,kmax)-(1.-mask_t(i,j,kmax))*9999.
        enddo ; enddo
       endif                  !--------->
       texte80(1)='hogcm' ; texte80(2)='m'                     ! variable ; units
       texte80(3:4)='OGCMbathy'
       texte80(5)='YX' ; texte80(7)=txt_type_ !01-05-14
       if(flag_offline_binary==0)call netcdf_main('_w')
       if(flag_offline_binary==1)call offlinegbin('_w') !18-06-18

       if(loop_netcdf==1) then !---------> 
        do j=0,jmax+1 ; do i=0,imax+1
         xy_t(i,j,0)=h_w(i,j)-xy_t(i,j,0)
         anyvar2d(i,j)=    xy_t(i,j,0)*mask_t(i,j,kmax)-(1.-mask_t(i,j,kmax))*9999.
           anyv3d(i,j,1,1)=xy_t(i,j,0)*mask_t(i,j,kmax)-(1.-mask_t(i,j,kmax))*9999.
        enddo ; enddo
       endif                  !--------->
       texte80(1)='h_minus_hogcm' ; texte80(2)='m'                     ! variable ; units
       texte80(3:4)='Sbathy_minus_OGCMbathy'
       texte80(5)='YX' ; texte80(7)=txt_type_ !01-05-14
       if(flag_offline_binary==0)call netcdf_main('_w')
       if(flag_offline_binary==1)call offlinegbin('_w') !18-06-18
      endif                                   !pmpmpmpmpm>

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
            anyvar2d(i,j)=h_w(i,j) ! ecrire en real
          anyv3d(i,j,1,1)=h_w(i,j) ! ecrire en double
      enddo
      enddo
      endif                  !--------->
      texte80(1)='h_w' ; texte80(2)='m'                     ! variable ; units
      texte80(3:4)='sea_floor_depth'
!     texte80(5)='YX' ; texte80(7)='real'
      texte80(5)='YX' ; texte80(7)=txt_type_ !01-05-14
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offlinegbin('_w') !18-06-18


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       anyvar2d(i,j)=dy_u(max0(i,1),max0(j,1))
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dy_u'   ; texte80(2)='m'                    ! variable ; units
      texte80(3:4)='cell_y_size_at_u_location'
      texte80(7)='real'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') !18-06-18

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; do i=0,imax+1
        anyvar2d(i,j)=dy_v(max0(i,1),max0(j,1)) !10-03-16
       enddo ; enddo
      endif                  !--------->
      texte80(1)='dy_v'   ; texte80(2)='m'                    ! variable ; units
      texte80(3:4)='cell_y_size_at_v_location'
      texte80(7)='real'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offlinegbin('_v') !18-06-18

      s_typevar=0 ; call donne_le_type(dy_t(1,1)) !18-08-22 
      if(s_typevar/=4.and.s_typevar/=8)stop 'type dy_t pas prevu'
      if(loop_netcdf==1) then !--------->
       if(s_typevar==4)anyvar2d(0:imax+1,0:jmax+1)    =dy_t(0:imax+1,0:jmax+1) 
       if(s_typevar==8)  anyv3d(0:imax+1,0:jmax+1,1,1)=dy_t(0:imax+1,0:jmax+1)
      endif                  !--------->
      if(s_typevar==4)texte80(7)='real'
      if(s_typevar==8)texte80(7)='double'
      texte80(1)='dy_t'   ; texte80(2)='m'                    ! variable ; units
      texte80(3:4)='cell_y_size_at_t_location'
      texte80(5)='YX'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; do i=0,imax+1
        anyvar2d(i,j)=dx_u(max0(i,1),max0(j,1))
       enddo ; enddo
      endif                  !--------->
      texte80(1)='dx_u'   ; texte80(2)='m'                    ! variable ; units
      texte80(3:4)='cell_x_size_at_u_location'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') !18-06-18

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       anyvar2d(i,j)=dx_v(max0(i,1),max0(j,1))
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dx_v'   ; texte80(2)='m'                    ! variable ; units
      texte80(3:4)='cell_x_size_at_v_location'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offlinegbin('_v') !18-06-18

      s_typevar=0 ; call donne_le_type(dx_t(1,1)) !18-08-22 
      if(s_typevar/=4.and.s_typevar/=8)stop 'type dx_t pas prevu'
      if(loop_netcdf==1) then !--------->
       if(s_typevar==4)anyvar2d(0:imax+1,0:jmax+1)    =dx_t(0:imax+1,0:jmax+1) 
       if(s_typevar==8)  anyv3d(0:imax+1,0:jmax+1,1,1)=dx_t(0:imax+1,0:jmax+1)
      endif                  !--------->
      if(s_typevar==4)texte80(7)='real'
      if(s_typevar==8)texte80(7)='double'
      texte80(1)='dx_t'   ; texte80(2)='m'                    ! variable ; units
      texte80(3:4)='cell_x_size_at_t_location'
      texte80(5)='YX'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18


      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=dxdy_t(i,j)
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dxdy_t'   ; texte80(2)='m2'                 ! variable ; units
      texte80(3:4)='cell_area'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

      if(loop_netcdf==1) then !---------> !27-07-13
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then
          anyvar2d(i,j)=sqrt(dxdy_t(i,j))
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='sqrt_dxdy'   ; texte80(2)='m'                 ! variable ; units
      texte80(3:4)='mesh_size'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

      if(loop_netcdf==1) then !---------> !25-12-13
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then
          anyvar2d(i,j)=dx_t(i,j)/(dy_t(i,j)+small1)
       else
          anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dx_over_dy'   ; texte80(2)='none'                 ! variable ; units
      texte80(3:4)='dx_over_dy'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

!.....
      if(loop_netcdf==1) then !--------->
        do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=i+par%timax(1)      !26-07-13
        enddo ; enddo
      endif                   !--------->
      texte80(1)='i_index_t'   ; texte80(2)='none'                 ! variable ; units
      texte80(3:4)='i_grid_index_t'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

!.....
      if(loop_netcdf==1) then !--------->
        do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=i+par%timax(1)-0.5
        enddo ; enddo
      endif                   !--------->
      texte80(1)='i_index_u'   ; texte80(2)='none'                 ! variable ; units
      texte80(3:4)='i_grid_index_u'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') !18-06-18

!.....
      if(loop_netcdf==1) then !--------->
        do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=i+par%timax(1)
        enddo ; enddo
      endif                   !--------->
      texte80(1)='i_index_v'   ; texte80(2)='none'                 ! variable ; units
      texte80(3:4)='i_grid_index_v'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offlinegbin('_v') !18-06-18

!.....
      if(loop_netcdf==1) then !--------->
        do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=j+par%tjmax(1)      !26-07-13
        enddo ; enddo
      endif                  !--------->
      texte80(1)='j_index_t'   ; texte80(2)='none'                 ! variable ; units
      texte80(3:4)='j_grid_index_t'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

!.....
      if(loop_netcdf==1) then !--------->
        do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=j+par%tjmax(1)      !26-07-13
        enddo ; enddo
      endif                  !--------->
      texte80(1)='j_index_u'   ; texte80(2)='none'                 ! variable ; units
      texte80(3:4)='j_grid_index_u'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') !18-06-18

!.....
      if(loop_netcdf==1) then !--------->
        do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=j+par%tjmax(1)-0.5
        enddo ; enddo
      endif                  !--------->
      texte80(1)='j_index_v'   ; texte80(2)='none'                 ! variable ; units
      texte80(3:4)='j_grid_index_v'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offlinegbin('_v') !18-06-18


      if(loop_netcdf==1) then !--------->                 !29-12-13

! Calcul des angles des axes de la grille:
      do j=1,jmax
      do i=1,imax

! Angle1
       dlon_di_=lon_t(i+1,j)-lon_t(i-1,j)
       if(dlon_di_<-pi)dlon_di_=dlon_di_+2.*pi
       if(dlon_di_> pi)dlon_di_=dlon_di_-2.*pi
       anyvar2d(i,j)=atan2(lat_t(i+1,j)-lat_t(i-1,j)  &
                        ,dlon_di_*cos(lat_t(i,j)))

! Angle2-angle1
       dlon_dj_=lon_t(i,j+1)-lon_t(i,j-1)
       if(dlon_dj_<-pi)dlon_dj_=dlon_dj_+2.*pi
       if(dlon_dj_> pi)dlon_dj_=dlon_dj_-2.*pi
       anyvar2d(i,j)=-anyvar2d(i,j)                  &
                  +atan2(lat_t(i,j+1)-lat_t(i,j-1)   &
                        ,dlon_dj_*cos(lat_t(i,j)))

! Conversion degres:
       anyvar2d(i,j)=anyvar2d(i,j)*rad2deg
       if(anyvar2d(i,j)<-180.)anyvar2d(i,j)=anyvar2d(i,j)+360.

! Masque terre/mer
       if(mask_t(i,j,kmaxp1)==0)anyvar2d(i,j)=filval

      enddo
      enddo
! CL
      do j=1,jmax
       anyvar2d(imax+1,j)=anyvar2d(imax,j)
       anyvar2d(0     ,j)=anyvar2d(1   ,j)
      enddo
      do i=0,imax+1
       anyvar2d(i,jmax+1)=anyvar2d(i,jmax)
       anyvar2d(i,0     )=anyvar2d(i,1   )
      enddo

      endif                  !--------->
      texte80(1)='angle_axis'   ; texte80(2)='degrees'          ! variable ; units
      texte80(3:4)='angle_axis'
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

      var_validmin=-2  ; var_validmax=99999
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
!      anyvar3dint(i,j,1)=par%rank                     !12-05-11
!      if(mask_t(i,j,kmax+1)==0)anyvar3dint(i,j,1)=-1  !12-05-11
       anyvar2d(i,j)=par%rank                     !12-05-11
       if(mask_t(i,j,kmax+1)==0)anyvar2d(i,j)=-1  !12-05-11
      enddo
      enddo
      endif                  !--------->
      texte80(1)='mpi_grid'   ; texte80(2)='none'         ! variable ; units
      texte80(3)='mpi subdomain number'                   ! long_name
      texte80(4)='mpi_subdomain_number'                   !  standard_name
      texte80(5)='YX' ; texte80(7)='short'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       anyvar2d(i,j)=par%subcycle                     !10-11-13
       if(mask_t(i,j,kmax+1)==0)anyvar2d(i,j)=filval
!      if(i==1   )anyvar2d(i,j)=filval
!      if(i==imax)anyvar2d(i,j)=filval
       if(j==1   )anyvar2d(i,j)=filval
       if(j==jmax)anyvar2d(i,j)=filval
      enddo
      enddo
      endif                  !--------->
      texte80(1)='par%subcycle' ; texte80(2)='count_per_principalcount' ! variable ; units
      texte80(3:4)='subcycling_frequency'                   ! long_name
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       anyvar2d(i,j)=dti_fw                    !10-11-13
       if(mask_t(i,j,kmax+1)==0)anyvar2d(i,j)=filval
!      if(i==1   )anyvar2d(i,j)=filval
!      if(i==imax)anyvar2d(i,j)=filval
       if(j==1   )anyvar2d(i,j)=filval
       if(j==jmax)anyvar2d(i,j)=filval
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dti_fw' ; texte80(2)='s' ! variable ; units
      texte80(3:4)='internal_mode_time_step'       ! long_name
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
! La valeur maximum theorique locale du pas de temps dte_lp (respectivement dte_fw) est l'expression en suivant
! en utilisant cfl_reduce=2 pour un schema FB (c.a.d remplacer 0.5*cfl_reduce par 1)
!        dte_lp_max=    cfl_reduce*sqrt(2.)*min(dx_t(i,j),dy_t(i,j))/(2.*sqrt(max(grav*(h_w(i,j)+cfl_sshmax),small1))+cfl_umax) 
!        dte_fw_max=0.5*cfl_reduce*sqrt(2.)*min(dx_t(i,j),dy_t(i,j))/(2.*sqrt(max(grav*(h_w(i,j)+cfl_sshmax),small1))+cfl_umax) 
       anyvar2d(i,j)= & !13-10-21
          mask_t(i,j,kmax)*sqrt(2.)*min(dx_t(i,j),dy_t(i,j))/(2.*sqrt(max(grav*(h_w(i,j)+cfl_sshmax),small1))+cfl_umax) &
      +(1-mask_t(i,j,kmax))*filval

      enddo
      enddo
      endif                  !--------->
      texte80(1)='CFL2D' ; texte80(2)='s'       ! variable ; units
      texte80(3:4)='stability_limit_for_external_time_step_(dte_fw)' ! long_name !13-10-21
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then
        anyvar2d(i,j)=upwindriver_t(i,j)
       else
        anyvar2d(i,j)=-9999.
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='upwindriver_t'   ; texte80(2)='none'         ! variable ; units
      texte80(3)='river mount area with upwind scheme'           ! long_name
      texte80(4)='river_mount_area_with_upwind scheme'           ! long_name
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

! write Land Sea Mask & River Input Number !01-11-14
      filval=-32767.
      var_validmin=-1. ; var_validmax=nriver
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       anyvar2d(i,j)=-mask_t(i,j,kmax)  &
                  +(1-mask_t(i,j,kmax))*filval
      enddo
      enddo
      do kr=1,nriver
       if(rivertrc_inout(kr)==1) then !----->
        anyvar2d(iriver(kr,1),jriver(kr,1))=real(kr)
       endif                          !----->
      enddo
      endif                  !--------->
      texte80(1)='LSMRIN' ; texte80(2)='integer'     !03-07-13
      texte80(3:4)='Sea=-1_Land=0_RiverInputNumber>0'
      texte80(5)='TYX'  ; texte80(7)='short'
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offlinegbin('_w') !18-06-18

!.....

! mask_t
      var_validmin=0. ; var_validmax=1.
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
!     do k=1,kmax
      do j=0,jmax+1 ; do i=0,imax+1
!      anyvar3d(i,j,k)=mask_t(i,j,k)
       anyvar2d(i,j)=mask_t(i,j,kmax)
      enddo         ; enddo
!     enddo
      endif                  !--------->
      texte80(1)='mask_t'   ; texte80(2)='none'         ! variable ; units
      texte80(3:4)='Land_Sea_Mask_T_nodes'
!     texte80(5)='ZYX' ; texte80(7)='short'             !27-06-14
      texte80(5)='YX'  ; texte80(7)='short'             !27-06-14
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

      var_validmin=-2. ; var_validmax=50.
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=mask_f(i,j,kmax)
      enddo
      enddo
      endif                  !--------->
      texte80(1)='mask_f'   ; texte80(2)='none'              ! variable ; units
      texte80(3:4)='Land_Sea_Mask_F_nodes'
!     texte80(5)='YX' ; texte80(7)='real'
      texte80(5)='YX' ; texte80(7)='short'                     !27-06-14
      if(flag_offline_binary==0)call netcdf_main('_f')
      if(flag_offline_binary==1)call offlinegbin('_f') !18-06-18

      if(allocated(zone1_mask)) then !m°v°m> !05-07-19
       var_validmin=-2. ; var_validmax=50.
       var_addoffset=0. ; var_scalefactor=1.
       if(loop_netcdf==1) then !--------->
       do j=0,jmax+1
       do i=0,imax+1
           anyvar2d(i,j)=zone1_mask(i,j)
       enddo
       enddo
       endif                  !--------->
       texte80(1)='zone1_mask'   ; texte80(2)='none'              !  variable ; units
       texte80(3:4)=texte80(1)
!      texte80(5)='YX' ; texte80(7)='real'
       texte80(5)='YX' ; texte80(7)='short'                     !27-06-14
       if(flag_offline_binary==0)call netcdf_main('_t')
       if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18
      endif                         !m°v°m> !05-07-19
!...
      if(allocated(zone2_mask)) then !m°v°m> 
       var_validmin=-2. ; var_validmax=50.
       var_addoffset=0. ; var_scalefactor=1.
       if(loop_netcdf==1) then !--------->
       do j=0,jmax+1
       do i=0,imax+1
           anyvar2d(i,j)=zone2_mask(i,j)
       enddo
       enddo
       endif                  !--------->
       texte80(1)='zone2_mask'   ; texte80(2)='none'              !  variable ; units
       texte80(3:4)=texte80(1)
!      texte80(5)='YX' ; texte80(7)='real'
       texte80(5)='YX' ; texte80(7)='short'                     !27-06-14
       if(flag_offline_binary==0)call netcdf_main('_t')
       if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18
      endif                         !m°v°m>
!...
      if(allocated(zone3_mask)) then !m°v°m> 
       var_validmin=-2. ; var_validmax=50.
       var_addoffset=0. ; var_scalefactor=1.
       if(loop_netcdf==1) then !--------->
       do j=0,jmax+1
       do i=0,imax+1
           anyvar2d(i,j)=zone3_mask(i,j)
       enddo
       enddo
       endif                  !--------->
       texte80(1)='zone3_mask'   ; texte80(2)='none'              !  variable ; units
       texte80(3:4)=texte80(1)
!      texte80(5)='YX' ; texte80(7)='real'
       texte80(5)='YX' ; texte80(7)='short'                     !27-06-14
       if(flag_offline_binary==0)call netcdf_main('_t')
       if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18
      endif                         !m°v°m>
!...
      if(allocated(zone4_mask)) then !m°v°m> 
       var_validmin=-2. ; var_validmax=50.
       var_addoffset=0. ; var_scalefactor=1.
       if(loop_netcdf==1) then !--------->
       do j=0,jmax+1
       do i=0,imax+1
           anyvar2d(i,j)=zone4_mask(i,j)
       enddo
       enddo
       endif                  !--------->
       texte80(1)='zone4_mask'   ; texte80(2)='none'              !  variable ; units
       texte80(3:4)=texte80(1)
!      texte80(5)='YX' ; texte80(7)='real'
       texte80(5)='YX' ; texte80(7)='short'                     !27-06-14
       if(flag_offline_binary==0)call netcdf_main('_t')
       if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18
      endif                         !m°v°m>
!...
      if(allocated(zone5_mask)) then !m°v°m> 
       var_validmin=-2. ; var_validmax=50.
       var_addoffset=0. ; var_scalefactor=1.
       if(loop_netcdf==1) then !--------->
       do j=0,jmax+1
       do i=0,imax+1
           anyvar2d(i,j)=zone5_mask(i,j)
       enddo
       enddo
       endif                  !--------->
       texte80(1)='zone5_mask'   ; texte80(2)='none'              !  variable ; units
       texte80(3:4)=texte80(1)
!      texte80(5)='YX' ; texte80(7)='real'
       texte80(5)='YX' ; texte80(7)='short'                     !27-06-14
       if(flag_offline_binary==0)call netcdf_main('_t')
       if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18
      endif                         !m°v°m>
!...
      if(allocated(zone6_mask)) then !m°v°m> 
       var_validmin=-2. ; var_validmax=50.
       var_addoffset=0. ; var_scalefactor=1.
       if(loop_netcdf==1) then !--------->
       do j=0,jmax+1
       do i=0,imax+1
           anyvar2d(i,j)=zone6_mask(i,j)
       enddo
       enddo
       endif                  !--------->
       texte80(1)='zone6_mask'   ; texte80(2)='none'              !  variable ; units
       texte80(3:4)=texte80(1)
!      texte80(5)='YX' ; texte80(7)='real'
       texte80(5)='YX' ; texte80(7)='short'                     !27-06-14
       if(flag_offline_binary==0)call netcdf_main('_t')
       if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18
      endif                         !m°v°m>


      var_validmin=1.  ; var_validmax=real(kmax+1) !03-12-14
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=kmin_w(i,j)
      enddo
      enddo
      endif                  !--------->
      texte80(1)='kmin_w'   ; texte80(2)='none'
      texte80(3)='first_level_from_bottom_w_location'
      texte80(4)=texte80(3)
      texte80(5)='YX' ; texte80(7)='short'
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offlinegbin('_w') !18-06-18

      var_validmin=1.  ; var_validmax=real(kmax+1) !03-12-14
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=kmerged_t(i,j)
      enddo ; enddo
      endif                  !--------->
      texte80(1)='kmerged_t'   ; texte80(2)='none'
      texte80(3)=texte80(1)
      texte80(4)=texte80(3)
      texte80(5)='YX' ; texte80(7)='short'
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offlinegbin('_w') !18-06-18

! mask_u
      var_validmin=1.  ; var_validmax=real(kmax+1) !03-12-14
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
         do j=0,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=mask_u(i,j,kmax)
        enddo ; enddo
      endif                  !--------->
      texte80(1)='mask_u'   ; texte80(2)='none'
      texte80(3)='Land_Sea_Mask_U_nodes'
      texte80(4)=texte80(3)
      texte80(5)='YX' ; texte80(7)='short'
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') !18-06-18

! kmin_u
      var_validmin=1.  ; var_validmax=real(kmax+1) !03-12-14
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=kmin_u(i,j)
      enddo
      enddo
      endif                  !--------->
      texte80(1)='kmin_u'   ; texte80(2)='none'
      texte80(3)='first_level_from_bottom_u_location'
      texte80(4)=texte80(3)
      texte80(5)='YX' ; texte80(7)='short'
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') !18-06-18

! kmerged_u !22-10-21
      var_validmin=1.  ; var_validmax=real(kmax+1) !03-12-14
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1 ; do i=1,imax+1
          anyvar2d(i,j)=kmerged_u(i,j)
      enddo ; enddo
      endif                  !--------->
      texte80(1)='kmerged_u'   ; texte80(2)='none'
      texte80(3:4)=texte80(1) ; texte80(5)='YX' ; texte80(7)='short'
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') !18-06-18

! kmergedr4_u !06-12-19
      var_validmin=1.  ; var_validmax=real(kmax+1) 
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
        do j=0,jmax+1 ; do i=1,imax+1
          anyvar2d(i,j)=kmergedr4_u(i,j)
        enddo ; enddo
      endif                  !--------->
      texte80(1)='kmergedr4_u'   ; texte80(2)='none'
      texte80(3)=texte80(1) ; texte80(4)=texte80(3)
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') 

! dsigmerged_u !02-11-21
      var_validmin=1.  ; var_validmax=real(kmax+1) 
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
        do j=0,jmax+1 ; do i=1,imax+1
          anyvar2d(i,j)=dsigmerged_u(i,j)
        enddo ; enddo
      endif                  !--------->
      texte80(1)='dsigmerged_u'   ; texte80(2)='none'
      texte80(3)=texte80(1) ; texte80(4)=texte80(3)
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') 

! pgfratio_u !03-11-21
      if(loop_netcdf==1) then !--------->
        do j=1,jmax   ; do i=1,imax+1
          anyvar2d(i,j)=pgfratio_u(i,j)
        enddo ; enddo
      endif                  !--------->
      texte80(1)='pgfratio_u'   ; texte80(2)='none'
      texte80(3)=texte80(1) ; texte80(4)=texte80(3)
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') 

! mask_v
      var_validmin=1.  ; var_validmax=real(kmax+1) !03-12-14
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; do i=0,imax+1
        anyvar2d(i,j)=mask_v(i,j,kmax)
       enddo         ; enddo
      endif                  !--------->
      texte80(1)='mask_v'   ; texte80(2)='none'
      texte80(3)='Land_Sea_Mask_V_nodes'
      texte80(4)=texte80(3)
      texte80(5)='YX' ; texte80(7)='short'
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offlinegbin('_v') !18-06-18

! kmin_v
      var_validmin=1.  ; var_validmax=real(kmax+1) !03-12-14
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
          anyvar2d(i,j)=kmin_v(i,j)
      enddo
      enddo
      endif                  !--------->
      texte80(1)='kmin_v'   ; texte80(2)='none'
      texte80(3)='first_level_from_bottom_v_location'
      texte80(4)=texte80(3)
      texte80(5)='YX' ; texte80(7)='short'
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offlinegbin('_v') !18-06-18

! kmerged_v
      var_validmin=1.  ; var_validmax=real(kmax+1) !03-12-14
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
      do j=1,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=kmerged_v(i,j)
      enddo ; enddo
      endif                  !--------->
      texte80(1)='kmerged_v'   ; texte80(2)='none'
      texte80(3:4)=texte80(1) ; texte80(5)='YX' ; texte80(7)='short'
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') !18-06-18

! kmergedr4_v !06-12-19
      var_validmin=1.  ; var_validmax=real(kmax+1) 
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
        do j=1,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=kmergedr4_v(i,j)
        enddo ; enddo
      endif                  !--------->
      texte80(1)='kmergedr4_v'   ; texte80(2)='none'
      texte80(3)=texte80(1) ; texte80(4)=texte80(3)
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offlinegbin('_v') 

! dsigmerged_v !02-11-21
      var_validmin=1.  ; var_validmax=real(kmax+1) 
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
        do j=1,jmax+1 ; do i=0,imax+1
          anyvar2d(i,j)=dsigmerged_v(i,j)
        enddo ; enddo
      endif                  !--------->
      texte80(1)='dsigmerged_v'   ; texte80(2)='none'
      texte80(3)=texte80(1) ; texte80(4)=texte80(3)
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offlinegbin('_v') 

! pgfratio_v !03-11-21
      if(loop_netcdf==1) then !--------->
        do j=1,jmax+1 ; do i=1,imax  
          anyvar2d(i,j)=pgfratio_v(i,j)
        enddo ; enddo
      endif                  !--------->
      texte80(1)='pgfratio_v'   ; texte80(2)='none'
      texte80(3)=texte80(1) ; texte80(4)=texte80(3)
      texte80(5)='YX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offlinegbin('_v') 

      if(coef_diss_mangrove>0.) then !mmmmmmmmmmmmmmmmmmmmmmmmmmm> !23-06-13
       filval=-9999.
       if(loop_netcdf==1) then !--------->  Mangrove Sea Mask
        anyvar2d(0,:)=0. ; anyvar2d(imax+1,:)=0. ; anyvar2d(:,0)=0. ; anyvar2d(:,jmax+1)=0.
        do j=1,jmax ! bornes coherentes avec dimensions mask_mangrove_t
        do i=1,imax
          if(mask_t(i,j,kmax+1)==1) then
           anyvar2d(i,j)=mask_mangrove_t(i,j)
          else
           anyvar2d(i,j)=-9999. !04-04-16
          endif
        enddo
        enddo
       endif                  !--------->
       texte80(1)='mask_mangrove_t'   ; texte80(2)='none'   ! variable ; units
       texte80(3)='mangrove mask'                           ! long_name
       texte80(4)='mangrove_mask'                           ! standard_name
       texte80(5)='YX' ; texte80(7)='real'
       if(flag_offline_binary==0)call netcdf_main('_t')
       if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18
      endif                          !mmmmmmmmmmmmmmmmmmmmmmmmmmm> !23-06-13

!     s_typevar=0 ; call donne_le_type(gridrotcos_t(1,1)) ; if(s_typevar/=4.and.s_typevar/=8) stop 'type gridrotcos_t pas prevu' !18-08-22
      if(loop_netcdf==1) then !--------->
                       anyvar2d(0:imax+1,0:jmax+1)    =gridrotcos_t(0:imax+1,0:jmax+1) 
!      if(s_typevar==4)anyvar2d(0:imax+1,0:jmax+1)    =gridrotcos_t(0:imax+1,0:jmax+1) 
!      if(s_typevar==8)  anyv3d(0:imax+1,0:jmax+1,1,1)=gridrotcos_t(0:imax+1,0:jmax+1)
      endif                  !--------->
      texte80(7)='real'
!     if(s_typevar==4)texte80(7)='real'
!     if(s_typevar==8)texte80(7)='double'
      texte80(1)='gridrotcos_t'   ; texte80(2)='none'              ! variable ; units
      texte80(3)='grid_rotation_cosinus_term'  
      texte80(4)='u_Oi=gridrotcos_t*u_eastward-gridrotsin_t*v_northward_v  v_Oj=gridrotsin_t*u_eastward+gridrotcos_t*v_northward'
!'
      texte80(5)='YX' 
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

!     s_typevar=0 ; call donne_le_type(gridrotsin_t(1,1)) ; if(s_typevar/=4.and.s_typevar/=8) stop 'type gridrotsin_t pas prevu' !18-08-22
      if(loop_netcdf==1) then !--------->
                       anyvar2d(0:imax+1,0:jmax+1)    =gridrotsin_t(0:imax+1,0:jmax+1) 
!      if(s_typevar==4)anyvar2d(0:imax+1,0:jmax+1)    =gridrotsin_t(0:imax+1,0:jmax+1) 
!      if(s_typevar==8)  anyv3d(0:imax+1,0:jmax+1,1,1)=gridrotsin_t(0:imax+1,0:jmax+1)
      endif                  !--------->
      texte80(7)='real'
!     if(s_typevar==4)texte80(7)='real'
!     if(s_typevar==8)texte80(7)='double'
      texte80(1)='gridrotsin_t'   ; texte80(2)='none'            ! variable ; units
      texte80(5)='YX'
      texte80(3)='grid_rotation_sinus_term' 
      texte80(4)='u_Oi=gridrotcos_t*u_eastward-gridrotsin_t*v_northward_v  v_Oj=gridrotsin_t*u_eastward+gridrotcos_t*v_northward'
!'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       anyvar2d(i,j)=atan2(gridrotsin_t(i,j),gridrotcos_t(i,j))*rad2deg
      enddo
      enddo
      endif                  !--------->
      texte80(1)='grid_angle'   ; texte80(2)='°'            ! variable ; units
      texte80(5)='YX' ; texte80(7)='real'
      texte80(3)='Angle_from_Oi_axis_to_OE_axis' ; &
         texte80(4)='https://docs.google.com/presentation/d/1FmAXNCdY_vL5KCUvkY0AxOQ6u145xQc1q-3cfXL1shU/edit#slide=id.p'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18


      if(loop_netcdf==1) then !--------->
      do k=1,kmax+1
      do j=0,jmax+1
      do i=0,imax+1
          anyvar3d(i,j,k)=depth_w(i,j,k)
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='depth_w'
      texte80(2)='m'                             ! units
      texte80(3:4)='depth_at_w_location'    ! long_name
      texte80(5)='ZYX' ; texte80(7)='real'
!     texte80(5)='ZYX' ; texte80(7)=txt_type_ !01-05-14
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offlinegbin('_w') !18-06-18


! REPERE SUPPRESSION V3D !13-03-20
!#ifdef bidon
      x1_r4=minval(upwzone0_t) !22-05-20
      call mpi_allreduce(x1_r4,x2_r4,1,mpi_real,mpi_min,par%comm2d ,ierr)
! Le test if(x2_r4<1.) en suivant pour ne pas ecrire inutilement un tableau inutilise
      if(x2_r4<1.) then !pmx> 
       if(loop_netcdf==1) then !--------->
       do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
           anyvar3d(i,j,k)=upwzone0_t(i,j,k)
       enddo ; enddo ; enddo
       endif                   !--------->
       texte80(1)='upwzone0_t'
       texte80(2)='none'                        ! units
       texte80(3:4)='upwind_zone0_3D'
       texte80(5)='ZYX' ; texte80(7)='real'
       if(flag_offline_binary==0)call netcdf_main('_t')
       if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18
      endif             !pmx> 

      if(loop_netcdf==1) then !--------->
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
          anyvar3d(i,j,k)=depth_t(i,j,k)
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='depth_t'
      texte80(2)='m'                        ! units
      texte80(3:4)='depth_at_t_location'
      texte80(5)='ZYX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

      filval=-9999. !25-04-20
      if(loop_netcdf==1) then !--------->
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
          anyvar3d(i,j,k)=max(k,kmin_w(i,j))*mask_t(i,j,kmax) &
                                  +filval*(1-mask_t(i,j,kmax))
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='klevel'
      texte80(2)='integer'                        ! units
      texte80(3:4)='vertical_level_number_t_location'
      texte80(5)='ZYX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

      if(loop_netcdf==1) then !--------->
      do k=1,kmax
      do j=0,jmax+1
      do i=1,imax+1 !20-12-15
          anyvar3d(i,j,k)=depth_u(i,j,k)
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='depth_u'
      texte80(2)='m'                        ! units
      texte80(3:4)='depth_at_u_location'
      texte80(5)='ZYX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') !18-06-18

      if(loop_netcdf==1) then !--------->
      do k=1,kmax
      do j=1,jmax+1 !20-12-15
      do i=0,imax+1
          anyvar3d(i,j,k)=depth_v(i,j,k)
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='depth_v'
      texte80(2)='m'                        ! units
      texte80(3:4)='depth_at_v_location'    ! long_name
      texte80(5)='ZYX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offlinegbin('_v') !18-06-18

      if(allocated(sig1dpom_t)) then !ooooo> !23-10-18
      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
! ON BIDOUILLE POUR APPARENCE FERRET
        x0=-sig1dpom_t(i,j,kmax)
        do k=1,kmax
         anyv1d(k,1)=sig1dpom_t(i,j,k)+x0 ! pour avoir sig(kmax)=0
        enddo
        x0=1./abs(anyv1d(1,1))
        do k=1,kmax
         anyv1d(k,1)=anyv1d(k,1)*x0 ! pour avoir sig(1)=-1
        enddo
        do k=1,kmax 
          anyvar3d(i,j,k) =anyv1d(k,1) ! sig1dpom_t(i,j,k) ! ecrire en real
           anyv3d(i,j,k,1)=anyv1d(k,1) ! sig1dpom_t(i,j,k) ! ecrire en double
        enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='sig1dpom_t'
      texte80(2)='?'                             ! units
      texte80(3:4)='s_coordinate_at_w_location'  ! long_name
      texte80(5)='ZYX' ; texte80(7)=txt_type_ !01-05-14
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18
      endif                       !ooooo> !28-11-14

      if(allocated(dsig_t)) then !ooooo> !28-11-14
      if(loop_netcdf==1) then !--------->
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
          anyvar3d(i,j,k) =dsig_t(i,j,k) ! ecrire en real
           anyv3d(i,j,k,1)=dsig_t(i,j,k) ! ecrire en double
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dsig_t'
      texte80(2)='none'                             ! units
      texte80(3:4)='delta_sigma_at_t_location'  ! long_name
      texte80(5)='ZYX' ; texte80(7)=txt_type_ !01-05-14
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18
      endif                       !ooooo> !28-11-14
      filval=-9999.
      if(loop_netcdf==1) then !--------->
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then
          anyvar3d(i,j,k)=dz_t(i,j,k,1)
       else
          anyvar3d(i,j,k)=filval
       endif
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dz_t'   ; texte80(2)='m'                            ! variable ; units
      texte80(3:4)='cell_thickness'                ! long_name
      texte80(5)='ZYX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offlinegbin('_t') !18-06-18

      filval=-9999.
      if(loop_netcdf==1) then !--------->
      do k=1,kmax
      do j=1,jmax                    !19-01-17
      do i=1,imax+1
       anyvar3d(i,j,k)=dz_u(i,j,k,1) !18-01-17
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dz_u'   ; texte80(2)='m'                            ! variable ; units
      texte80(3:4)='cell_thickness_u'                ! long_name
      texte80(5)='ZYX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') !18-06-18

      filval=-9999.
      if(loop_netcdf==1) then !--------->
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax                    !18-01-17
       anyvar3d(i,j,k)=dz_v(i,j,k,1) !18-01-17
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='dz_v'   ; texte80(2)='m'                            ! variable ; units
      texte80(3:4)='cell_thickness_v'                ! long_name
      texte80(5)='ZYX' ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offlinegbin('_v') !18-06-18

#ifdef bidon
! mask_w !24-03-16
      var_validmin=0. ; var_validmax=1.
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
       do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
        anyvar3d(i,j,k)=mask_t(i,j,k)
       enddo ; enddo ; enddo
       anyvar3d(:,:,kmax+1)=anyvar3d(:,:,kmax)
      endif                  !--------->
      texte80(1)='mask_w'   ; texte80(2)='none'        ! variable ; units
      texte80(3)='omega land sea binary mask'          ! long_name
      texte80(4)='omega_land_sea_binary_mask'          ! standard_name
      texte80(5)='ZYX' ; texte80(7)='short'            
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offlinegbin('_w') !18-06-18
#endif

#ifdef bidon
! mask_u
      var_validmin=0. ; var_validmax=1.
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
       anyvar3d(i,j,k)=mask_u(i,j,k)
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='mask_u'   ; texte80(2)='none'              ! variable ; units
      texte80(3)='Oi vel land sea binary mask'
      texte80(4)='Oi_vel_land_sea_binary_mask'
      texte80(5)='ZYX' ; texte80(7)='short'                  !27-06-14
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offlinegbin('_u') !18-06-18

! mask_v
      var_validmin=0. ; var_validmax=1.
      var_addoffset=0. ; var_scalefactor=1.
      if(loop_netcdf==1) then !--------->
      do k=1,kmax
      do j=0,jmax+1
      do i=0,imax+1
       anyvar3d(i,j,k)=mask_v(i,j,k)
      enddo
      enddo
      enddo
      endif                  !--------->
      texte80(1)='mask_v'   ; texte80(2)='none'              ! variable ; units
      texte80(3)='Oj vel land sea binary mask'
      texte80(4)='Oj_vel_land_sea_binary_mask'
      texte80(5)='ZYX' ; texte80(7)='short'                  !27-06-14
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offlinegbin('_v') !18-06-18
#endif

! REPERE SUPPRESSION V3D
!#endif


      if(loop_netcdf==0) then !>>>>>>>>>>>>>>>>>>>>

! Definition of variables: done.
      call netcdf_general_attributes(ncid) !03-04-12
      status=nfmpi_enddef(ncid)
      endif                  !>>>>>>>>>>>>>>>>>>>

      status=nfmpi_close(ncid)


!      endif                         !procprocprocproc>
!#ifdef parallele
!     call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!      call barriere(looproc_   ,3)
!#endif
!      enddo ! fin de boucle sur looproc_

      enddo  ! fin de boucle sur loop_netcdf

      end subroutine offline_write_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine offline_write_var    !29-12-11
      implicit none
      double precision tem_,tem2_,tem3_,tem4_,sal_,sal2_
      integer looproc_,looprocmax_
!      include 'netcdf.inc'
#ifdef synopsis
       subroutinetitle='offline_write_var'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      filval=-9999.

      do loop_netcdf=flag_offline_binary,1 ! do loop_netcdf=0,1 si netcdf, do loop_netcdf=1,1 si binaire

      count_netcdfvar=0

      if(flag_offline_binary==0) then !netcdf-case->

       if(loop_netcdf==0) then  !§§§§§§§>
        status=nfmpi_create(par%comm2d,texte250,nf_clobber+ NF_64BIT_OFFSET,MPI_INFO_NULL,ncid) !27-05-18
       else                     !§§§§§§§>
        status=nfmpi_open(par%comm2d,texte250,nf_write+ NF_64BIT_OFFSET,MPI_INFO_NULL,ncid)
       endif                    !§§§§§§§>
       if(status/=0) then
        write(6,'(a,a)')'Error nfmpi_create ',trim(texte250)
        stop ' stop module_offline'
       endif
       if(loop_netcdf==0)call netcdf_dim

! Reset VARDIM
      vardim(1)=i_t_dim ; vardim(2)=j_t_dim
      vardim(3)=k_t_dim ; vardim(4)=time_dim

! time:
      k0=1 ; vardim(1)=time_dim                                         ! 1D ; dim1
      texte80(1)='time'                                                 ! variable
      call kount_to_date(0)                          !16-11-09  time origin corresponds to kount=0
!     write(texte80(2),'(a14,i4,a1,a3,4(a1,i2))')                    & !units
!     'seconds since ',i5,'-',month(i6)(1:3),'-',i7,' ',i3,':',i2,':',i1

      write(texte80(2),'(a14,i4,a1,i2,4(a1,i2))')                    & !units
      'seconds since ',i5,'-',i6,'-',i7,' ',i3,':',i2,':',i1
      if(texte80(2)(20:20)==' ')texte80(2)(20:20)='0'
      if(texte80(2)(23:23)==' ')texte80(2)(23:23)='0'
      if(texte80(2)(26:26)==' ')texte80(2)(26:26)='0'
      if(texte80(2)(29:29)==' ')texte80(2)(29:29)='0'
      if(texte80(2)(32:32)==' ')texte80(2)(32:32)='0'

      texte80(8)=texte80(2)(14:33)                                      ! time_origin : kount=0
      texte80(9)='gregorian'                                            ! calendar
      texte80(3:4)='time'                                               ! long_name
      texte80(5)='T'  ; texte80(6)='time'                               ! axis ; associate
      texte80(7)='double'
      if(flag_offline_binary==0)call netcdf_main('_t')

! length of the time average
      k0=1 ; vardim(1)=time_dim                                         ! 1D ; dim1
      texte80(1)='cumulativetime'
      texte80(2)='hours'
      texte80(3:4)='length_of_the_time_average'
      texte80(5)='T'  ; texte80(6)='time'                               ! axis ; associate
      texte80(7)='double'
      if(flag_offline_binary==0)call netcdf_main('_t')

! write variables dimensions:
      call graph_out_trueaxis_nfmpi

      endif                           !netcdf-case->

      texte80(10)='none'

      if(flag_maxbotstress==1) then ! m[°v°]m > !09-11-17
       filval=-32767.
       if(loop_netcdf==1) then !--------->
        anyvar2d(0     ,:     )=filval
        anyvar2d(imax+1,:     )=filval
        anyvar2d(:     ,jmax+1)=filval
        anyvar2d(:     ,0     )=filval
        do j=1,jmax   ; do i=1,imax  
         anyvar2d(i,j)=maxbotstress_w(i,j)*mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
        enddo         ; enddo
       endif                  !--------->
       texte80(1)='maxbotstress_w' ; texte80(2)='?'
       texte80(3)='maxbotstress_w'
       texte80(4)='maxbotstress_w'
       texte80(5)='TYX'  ; texte80(7)='real'
       if(flag_offline_binary==0)call netcdf_main('_w')
       if(flag_offline_binary==1)call offline_bin('_w') !08-06-18

       if(loop_netcdf==1) then !--------->  !01-02-21
        anyvar2d(0     ,:     )=filval
        anyvar2d(imax+1,:     )=filval
        anyvar2d(:     ,jmax+1)=filval
        anyvar2d(:     ,0     )=filval
        do j=1,jmax   ; do i=1,imax
         anyvar2d(i,j)=stresswave_w(i,j)*mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
        enddo         ; enddo
       endif                  !--------->
       texte80(1)='tauw' ; texte80(2)='?'
       texte80(3)='tauw'
       texte80(4)='tauw'
       texte80(5)='TYX'  ; texte80(7)='real'
       if(flag_offline_binary==0)call netcdf_main('_w')
       if(flag_offline_binary==1)call offline_bin('_w') 

       if(loop_netcdf==1) then !--------->
        anyvar2d(0     ,:     )=filval
        anyvar2d(imax+1,:     )=filval
        anyvar2d(:     ,jmax+1)=filval
        anyvar2d(:     ,0     )=filval
        do j=1,jmax   ; do i=1,imax
         anyvar2d(i,j)=stressc_w(i,j)*mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
        enddo         ; enddo
       endif                  !--------->
       texte80(1)='tauc' ; texte80(2)='?'
       texte80(3)='tauc'
       texte80(4)='tauc'
       texte80(5)='TYX'  ; texte80(7)='real'
       if(flag_offline_binary==0)call netcdf_main('_w')
       if(flag_offline_binary==1)call offline_bin('_w') !  01-02-21 fin

      endif                         ! m[°v°]m > !09-11-17

      if(ofl_surflux==1) then !-cl-cl-cl-cl-cl-cl-cl->
! Flux de surface moyen

      filval=-32767.

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(  mask_t(i,j,kmax)==1) then
          anyvar2d(i,j)=slhf_aver_w(i,j)
       else
          anyvar2d(i,j)=filval
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='slhf' ; texte80(2)='W/m**2'
      texte80(3)='Latent_Heat_Flux'
      texte80(4)='Latent_Heat_Flux'
      texte80(5)='TYX'  ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offline_bin('_w') !08-06-18

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(  mask_t(i,j,kmax)==1) then
          anyvar2d(i,j)=sshf_aver_w(i,j)
       else
          anyvar2d(i,j)=filval
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='sshf' ; texte80(2)='W/m**2'
      texte80(3)='Sensible_Heat_Flux'
      texte80(4)='Sensible_Heat_Flux'
      texte80(5)='TYX'  ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offline_bin('_w') !08-06-18

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(  mask_t(i,j,kmax)==1) then
          anyvar2d(i,j)=snsf_aver_w(i,j)
       else
          anyvar2d(i,j)=filval
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='snsf' ; texte80(2)='W/m**2'
      texte80(3)='Longwave_Radiative_Flux'
      texte80(4)='Longwave_Radiative_Flux'
      texte80(5)='TYX'  ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offline_bin('_w') !08-06-18

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(  mask_t(i,j,kmax)==1) then
          anyvar2d(i,j)=ssr_aver_w(i,j)
       else
          anyvar2d(i,j)=filval
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='ssr' ; texte80(2)='W/m**2'
      texte80(3)='Solar_Radiative_Flux'
      texte80(4)='Solar_Radiative_Flux'
      texte80(5)='TYX'  ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offline_bin('_w') !08-06-18

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(  mask_t(i,j,kmax)==1) then
          anyvar2d(i,j)=precipi_aver_w(i,j)
       else
          anyvar2d(i,j)=filval
       endif
      enddo
      enddo
      endif                  !--------->
      texte80(1)='precip' ; texte80(2)='m/s' !27-09-21
      texte80(3)='Precipitation'
      texte80(4)='Precipitation'
      texte80(5)='TYX'  ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offline_bin('_w') !08-06-18

      if(loop_netcdf==1) then !--------->
        do j=0,jmax+1
        do i=0,imax+1
          anyvar2d(i,j)=wstress_aver_u(i,j)*mask_u(i,j,kmax) &
                                 +filval*(1-mask_u(i,j,kmax))
        enddo
        enddo
        call obc_int_anyvar2d('u1') !16-03-15
      endif                  !--------->
      texte80(1)='wstressu' ; texte80(2)='N/m**2'
      texte80(3)='X_Windstress'
      texte80(4)='X_Windstress'
      texte80(5)='TYX'  ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offline_bin('_u') !08-06-18

      if(loop_netcdf==1) then !--------->
        do j=0,jmax+1
        do i=0,imax+1
          anyvar2d(i,j)=wstress_aver_v(i,j)*mask_v(i,j,kmax) &
                                 +filval*(1-mask_v(i,j,kmax))
        enddo
        enddo
        call obc_int_anyvar2d('v1') !16-03-15
      endif                  !--------->
      texte80(1)='wstressv' ; texte80(2)='N/m**2'
      texte80(3)='Y_Windstress'
      texte80(4)='Y_Windstress'
      texte80(5)='TYX'  ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offline_bin('_v') !08-06-18

      endif                   !-cl-cl-cl-cl-cl-cl-cl->

!.....
! write time-averaged sea surface elevation (inverse barometer effect removed)
      if(.not.allocated(sshofl_w))stop 'Error: sshofl not allocated'
      filval=-32767.
      if(ofl_sshtype=='short') then !oooooooo> !02-05-14
       var_validmin=-3. ; var_validmax=3.
       var_addoffset=0.5*(var_validmin+var_validmax)
       var_scalefactor=abs(var_validmax-var_addoffset)/(abs(filval)-1.)
      else                       !oooooooo>
       var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
      endif                      !oooooooo>
      inv_scalefactor=1./var_scalefactor

      if(loop_netcdf==1) then !--------->
      do j=0,jmax+1
      do i=0,imax+1
       if(mask_t(i,j,kmax+1)==1) then
        anyvar2d(i,j)=inv_scalefactor*(-var_addoffset        &
        +min(max(sshofl_w(i,j,1),var_validmin),var_validmax))
       else
        anyvar2d(i,j)=filval
       endif
      enddo
      enddo
! Attention que le masque vient d'effacer la valeur de zta à la source des
! fleuves, or cette info est importante pour la bonne conversion de vhz en vel
! et (surtout) reciproquement. On reintroduit donc ces valeurs:
      do kr=1,nriver
      if(rivertrc_inout(kr)==1) then !-mpi-river-mpi-river-> !04-12-13
       anyvar2d(iriver(kr,1),jriver(kr,1))=           &
        inv_scalefactor*(-var_addoffset               &
        +min(max(                                     &
               sshofl_w(iriver(kr,1),jriver(kr,1),1), &
               var_validmin),                  &
               var_validmax))
      endif                       !-mpi-river-mpi-river-> !04-12-13
      enddo
      endif                  !--------->
!     texte80(1)='ssh-ib' ; texte80(2)='m'     !06-04-12
      texte80(1)='ssh_ib' ; texte80(2)='m'     !03-07-13
      texte80(3)='Sea surface height minus inverse barometer'
      texte80(4)='sea_surface_height_above_geoid'
!     texte80(5)='TYX'  ; texte80(7)=ofl_type
      texte80(5)='TYX'  ; texte80(7)=ofl_sshtype !02-05-14
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offline_bin('_w') !08-06-18
!.....

!.....
! write instantaneous ssh
      filval=-9999.
      var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
      inv_scalefactor=1./var_scalefactor
      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; do i=0,imax+1
        anyvar2d(i,j)=ssh_int_w(i,j,1)
       enddo ; enddo
      endif                   !--------->
      texte80(1)='ssh_inst' ; texte80(2)='m' 
      texte80(3:4)='instantaneous_ssh'
      texte80(5)='TYX'  ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offline_bin('_w') !08-06-18
!.....

!.....
! write time-averaged surface vertical velocity (for water budget) !04-07-15
      if(.not.allocated(w0mofl_w))stop 'Error: w0mofl not allocated'
      filval=-9999.
      var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
      inv_scalefactor=1./var_scalefactor

      if(loop_netcdf==1) then !--------->
       do j=0,jmax+1 ; do i=0,imax+1
        anyvar2d(i,j)=w0mofl_w(i,j,1)*mask_t(i,j,kmax+1)
       enddo ; enddo
      endif                   !--------->
      texte80(1)='w0m' ; texte80(2)='m/s' 
      texte80(3:4)='surface_vertical_velocity'
      texte80(5)='TYX'  ; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offline_bin('_w') !08-06-18
!.....

!.....
      if(flag_groundwater==1) then !GRWGRWGWR> 11-04-19
! write time-averaged bottom vertical velocity
       if(.not.allocated(w_keq1_ofl_w)) &
            stop 'Error: w_keq1_ofl not allocated'
       filval=-9999.
       var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
       inv_scalefactor=1./var_scalefactor

       if(loop_netcdf==1) then !--------->
        do j=1,jmax ; do i=1,imax
         anyvar2d(i,j)=w_keq1_ofl_w(i,j,1)*mask_t(i,j,kmax+1)
        enddo ; enddo
       endif                   !--------->
       texte80(1)='wbottom' ; texte80(2)='m/s'
       texte80(3:4)='bottom_vertical_velocity'
       texte80(5)='TYX'  ; texte80(7)='real' ; texte80(11)='mask_obc_z1'
       if(flag_offline_binary==0)call netcdf_main('_w')
       if(flag_offline_binary==1)call offline_bin('_w') !08-06-18
      endif                        !GRWGRWGWR> 11-04-19
!.....

!.....
! write time-averaged convective_surface_layer_vertical_index !09-02-19
      if(flag_ksloffline==1) then !pmx>
      if(.not.allocated(kslofl_t))stop 'Error: kslofl not allocated'
      filval=-9999.
      var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
      inv_scalefactor=1./var_scalefactor

      if(loop_netcdf==1) then !--------->
       do j=1,jmax   ; do i=1,imax  
        anyvar2d(i,j)=kslofl_t(i,j,1)*mask_t(i,j,kmax+1)
       enddo ; enddo
      endif                   !--------->
      texte80(1)='ksl' ; texte80(2)='???' 
      texte80(3:4)='convective_surface_layer_vertical_index'
      texte80(5)='TYX'  ; texte80(7)='real'
      texte80(11)='mask_obc_z1'
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offline_bin('_w') 
      endif                       !pmx>
!.....

!.....
! write temperature:
      filval=-32767.
      if(ofl_type=='short') then !oooooooo>
       var_validmin=-2. ; var_validmax=50.
       var_addoffset=0.5*(var_validmin+var_validmax)
       var_scalefactor=abs(var_validmax-var_addoffset)/(abs(filval)-1.)
      else                       !oooooooo>
       var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
      endif                      !oooooooo>
      inv_scalefactor=1./var_scalefactor

      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax
        do i=1,imax

         anyvar3d(i,j,k)=inv_scalefactor*(-var_addoffset+        &
         min(max(temofl_t(i,j,k,1),var_validmin),var_validmax)   &
                          )*mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval !20-12-14
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='tem'   ; texte80(2)='degrees_Celsius'         ! variable ; units
      texte80(3:4)='sea_water_potential_temperature'
      texte80(5)='TZYX' ; texte80(7)=ofl_type
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offline_bin('_t') !08-06-18

!.....
! write salinity:
      filval=-32767.
      if(ofl_type=='short') then !oooooooo>
       var_validmin=0. ; var_validmax=41.
       var_addoffset=0.5*(var_validmin+var_validmax)
       var_scalefactor=abs(var_validmax-var_addoffset)/(abs(filval)-1.)
      else                       !oooooooo>
       var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
      endif                      !oooooooo>
      inv_scalefactor=1./var_scalefactor
      if(loop_netcdf==1) then !=======>
        do k=1,kmax
        do j=1,jmax
        do i=1,imax

         anyvar3d(i,j,k)=inv_scalefactor*(-var_addoffset+        &
         min(max(salofl_t(i,j,k,1),var_validmin),var_validmax)   &
                          )*mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval !20-12-14
        enddo
        enddo
        enddo
      endif                  !=======>
      texte80(1)='sal'   ; texte80(2)='1e-3'              ! variable ; units
      texte80(3:4)='sea_water_salinity'
      texte80(5)='TZYX' ; texte80(7)=ofl_type

      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offline_bin('_t') !08-06-18

!.....
! write bio_t !09-09-19
      if(ofl_bio==1) then !biobiobiobiobio>
       do vb=1,vbmax

        filval=-32767.
        if(ofl_type=='short') then !oooooooo>
         var_validmin=0. ; var_validmax=100.
         var_addoffset=0.5*(var_validmin+var_validmax)
         var_scalefactor=abs(var_validmax-var_addoffset)/(abs(filval)-1.)
        else                       !oooooooo>
         var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
        endif                      !oooooooo>
        inv_scalefactor=1./var_scalefactor
        if(loop_netcdf==1) then !=======>
         do k=1,kmax
         do j=1,jmax
         do i=1,imax
          anyvar3d(i,j,k)=inv_scalefactor*(-var_addoffset+        &
          min(max(bioofl_t(i,j,k,vb),var_validmin),var_validmax)   &
                           )*mask_t(i,j,kmax)+(1-mask_t(i,j,kmax))*filval
         enddo
         enddo
         enddo
        endif                  !=======>
        write(texte80(1),'(a3,i0)')'bio',vb
        texte80(2)='???'              ! variable ; units
        texte80(3:4)=texte80(1)
        texte80(5)='TZYX' ; texte80(7)=ofl_type

        if(flag_offline_binary==0)call netcdf_main('_t')
        if(flag_offline_binary==1)call offline_bin('_t') !08-06-18

       enddo ! vb loop

       ! write epaisseur de couche sedim  ! 11-09-19
        filval=-9999.
        var_validmin=0. ; var_validmax=30. ; var_scalefactor=1. ; var_addoffset=0.
        inv_scalefactor=1./var_scalefactor

        if(loop_netcdf==1) then !--------->
         do j=1,jmax   ; do i=1,imax  
          anyvar2d(i,j)=hsedofl_t(i,j)*mask_t(i,j,kmax+1) &
                        +filval*(1-mask_t(i,j,kmax+1))
         enddo ; enddo
        endif                   !--------->
        texte80(1)='h_sed' ; texte80(2)='m' 
        texte80(3:4)='sediment_layer_thickness'
        texte80(5)='TYX'  ; texte80(7)='real'
        texte80(11)='mask_obc_z1'
        if(flag_offline_binary==0)call netcdf_main('_w')
        if(flag_offline_binary==1)call offline_bin('_w') 

      endif               !biobiobiobiobio>
!.....
      if(ofl_rhp==1) then !>>>>>>>>>>>>>>>>>
      filval=-32767.
      if(ofl_type=='short') then !oooooooo>
       var_validmin=0. ; var_validmax=50.
       var_addoffset=0.5*(var_validmin+var_validmax)
       var_scalefactor=abs(var_validmax-var_addoffset)/(abs(filval)-1.)
      else                       !oooooooo>
       var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
      endif                      !oooooooo>
      inv_scalefactor=1./var_scalefactor
      if(loop_netcdf==1) then !=======>
       if(eos_author==3) then
        if(eos_tkezref>=0.) then
! Note que si eos_tkezref<0 la profondeur de reference sera 0
          call equation_of_state_potzref_jmfwg('ofl',1) !20-10-14
        else
          call equation_of_state_potloc_jmfwg('ofl',1) !23-02-22
        endif
       else
        stop 'module_offline cas eos pas prevu'
       endif
       do k=1,kmax
       do j=1,jmax
       do i=1,imax
       if(  mask_t(i,j,kmax)==1) then !20-12-14

         anyvar3d(i,j,k)=inv_scalefactor*(-var_addoffset+      &
         min(max(anyv3d(i,j,k,1)+rho-1000.,var_validmin),var_validmax))

        else
         anyvar3d(i,j,k)=filval
        endif
       enddo
       enddo
       enddo
      endif                  !=======>
      texte80(1)='rhp'   ; texte80(2)='kg m**-3'              ! variable ; units
      write(texte30,'(i0)')nint(max(eos_tkezref,zero))                  !20-10-14
      texte80(3:4)='sea_water_potential_density_reference_level_' &     !21-01-17
                  //trim(texte30)//'_meters'
      texte80(5)='TZYX' ; texte80(7)=ofl_type
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offline_bin('_t') !08-06-18
      endif               !>>>>>>>>>>>>>>>>>

!.....
! write dz
      if(ale_selected==1) then !aaaaaaaaaaaaaaaaaaaaaaaaaa>

      filval=-32767.
      if(ofl_type=='short') then !oooooooo>
       var_validmin=-3. ; var_validmax=4. ! Attention pour echelle log10
       var_addoffset=0.5*(var_validmin+var_validmax)
       var_scalefactor=abs(var_validmax-var_addoffset)/(abs(filval)-1.)
      else                       !oooooooo>
       var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
      endif                      !oooooooo>
      inv_scalefactor=1./var_scalefactor
      if(loop_netcdf==1) then !=======>
        if(ofl_type=='short') then !oooooooo>
         do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
          anyvar3d(i,j,k)=log10(dzofl_t(i,j,k,1))
          anyvar3d(i,j,k)=inv_scalefactor*(-var_addoffset+ &
          min(max(anyvar3D(i,j,k),var_validmin),var_validmax))
         enddo ; enddo ; enddo
        else                       !oooooooo>
         do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
          anyvar3d(i,j,k)=dzofl_t(i,j,k,1)
         enddo ; enddo ; enddo
        endif                      !oooooooo>
      endif                  !=======>
! ATTENTION si dz change un jour de nom: la condition netcdf_obc ne doit pas s'appliquer
      if(ofl_type=='short') then !oooooooo>
        texte80(1)='log10dz' ; texte80(2)='m log10'
        texte80(3:4)='cell_thickness_log10'
        texte80(10)='dz=10**(log10dz*scale_factor+add_offset)'
      else                       !oooooooo>
        texte80(1)='dz' ; texte80(2)='m'
        texte80(3:4)='cell_thickness'
      endif                      !oooooooo>
      texte80(5)='TZYX' ; texte80(7)=ofl_type
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offline_bin('_t') !08-06-18

      endif                    !aaaaaaaaaaaaaaaaaaaaaaaaaa>

!.....
! write Kz
      filval=-32767.
      if(ofl_type=='short') then !oooooooo>
       var_validmin=-7. ; var_validmax=2. ! Attention pour echelle log10
       var_addoffset=0.5*(var_validmin+var_validmax)
       var_scalefactor=abs(var_validmax-var_addoffset)/(abs(filval)-1.)
      else                       !oooooooo>
       var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
      endif                      !oooooooo>
      inv_scalefactor=1./var_scalefactor

      if(loop_netcdf==1) then !=======>
! note: division par sumarchive car cette variable n'est pas passe (pour
! convenance) par routine end_of_average
       x1=1./sumarchive
        if(ofl_type=='short') then !oooooooo>
         do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
          if(mask_t(i,j,k)==1) then  !******>
           anyvar3d(i,j,k)=log10(max(dfvofl_w(i,j,k,-1)*x1,small1)) !05-11-14
           anyvar3d(i,j,k)=inv_scalefactor*(-var_addoffset+ &
           min(max(anyvar3D(i,j,k),var_validmin),var_validmax))
          else                       !******>
           anyvar3d(i,j,k)=filval
          endif                      !******>
         enddo ; enddo ; enddo
        else                       !oooooooo>
         do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
          if(mask_t(i,j,k)==1) then  !******>
           anyvar3d(i,j,k)=dfvofl_w(i,j,k,-1)*x1
          else                       !******>
           anyvar3d(i,j,k)=filval
          endif                      !******>
         enddo ; enddo ; enddo
        endif                      !oooooooo>
      endif                  !=======>
      if(ofl_type=='short') then !oooooooo>
       texte80(1)='log10difv' ; texte80(2)='m2 s-1 log10'
       texte80(3:4)='ocean_vertical_tracer_diffusivity_log10'
       texte80(10)='Kz=10**(log10difv*scale_factor+add_offset)'
      else                       !oooooooo>
       texte80(1)='difv'    ; texte80(2)='m2 s-1'
       texte80(3:4)='ocean_vertical_tracer_diffusivity'
      endif                      !oooooooo>
      texte80(5)='TZYX' ; texte80(7)=ofl_type
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offline_bin('_w') !08-06-18
      texte80(10)='none'

!.....
! write tken
      if(ofl_tke==1) then !>>>>>>>>>>>>>>>>>>>>>>>> !05-03-15
      filval=-32767.
      if(ofl_type=='short') then !oooooooo>
       var_validmin=-7. ; var_validmax=2. ! Attention pour echelle log10
       var_addoffset=0.5*(var_validmin+var_validmax)
       var_scalefactor=abs(var_validmax-var_addoffset)/(abs(filval)-1.)
      else                       !oooooooo>
       var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
      endif                      !oooooooo>
      inv_scalefactor=1./var_scalefactor

      if(loop_netcdf==1) then !=======>
! note: division par sumarchive car cette variable n'est pas passe (pour
! convenance) par routine end_of_average
       x1=1./sumarchive
        if(ofl_type=='short') then !oooooooo>
         do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
          if(mask_t(i,j,k)==1) then  !******>
           anyvar3d(i,j,k)=log10(max(tkeofl_w(i,j,k,-1)*x1,small1)) !05-11-14
           anyvar3d(i,j,k)=inv_scalefactor*(-var_addoffset+ &
           min(max(anyvar3D(i,j,k),var_validmin),var_validmax))
          else                       !******>
           anyvar3d(i,j,k)=filval
          endif                      !******>
         enddo ; enddo ; enddo
        else                       !oooooooo>
         do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
          if(mask_t(i,j,k)==1) then  !******>
           anyvar3d(i,j,k)=tkeofl_w(i,j,k,-1)*x1
          else                       !******>
           anyvar3d(i,j,k)=filval
          endif                      !******>
         enddo ; enddo ; enddo
        endif                      !oooooooo>
      endif                  !=======>
      if(ofl_type=='short') then !oooooooo>
       texte80(1)='log10tke' ; texte80(2)='m2 s-2 log10'
       texte80(3:4)='turbulent_kinetic_energy_log10'
       texte80(10)='tke=10**(log10tke*scale_factor+add_offset)'
      else                       !oooooooo>
       texte80(1)='tke'    ; texte80(2)='m2 s-2'
       texte80(3:4)='turbulent_kinetic_energy'
      endif                      !oooooooo>
      texte80(5)='TZYX' ; texte80(7)=ofl_type
      if(flag_offline_binary==0)call netcdf_main('_w')
      if(flag_offline_binary==1)call offline_bin('_w') !08-06-18
      texte80(10)='none'
      endif               !>>>>>>>>>>>>>>>>>>>>>>>> !05-03-15

!.....
! write U
      filval=-32767. ! Pour calculer var_scalefactor puis sera mis à zéro
      if(ofl_type=='short') then !oooooooo>
       var_validmin=-3. ; var_validmax=3.
       var_addoffset=0.5*(var_validmin+var_validmax)
       var_scalefactor=abs(var_validmax-var_addoffset)/(abs(filval)-1.)
      else                       !oooooooo>
       var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
      endif                      !oooooooo>
      inv_scalefactor=1./var_scalefactor
      filval=0.
      if(loop_netcdf==1) then !=======>

       if(ale_selected==1) then !sssss>
        do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
         anyv3d(i,j,k,1)=0.5*(dzofl_t(i,j,k,1)+dzofl_t(i-1,j,k,1))
        enddo ; enddo ; enddo
       else                     !sssss>
!       if(fgrid_or_wgrid==wgrid_case) then !wwwww>
         do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
          anyv3d(i,j,k,1)=0.5*(                                        &
          dsig_t(i  ,j,k)*(h_w(i  ,j)+sshofl_w(i  ,j,1)) &
         +dsig_t(i-1,j,k)*(h_w(i-1,j)+sshofl_w(i-1,j,1)))
!         (sigma_w(i  ,j,k+1)-sigma_w(i  ,j,k))*(h_w(i  ,j)+sshofl_w(i  ,j,1)) &
!        +(sigma_w(i-1,j,k+1)-sigma_w(i-1,j,k))*(h_w(i-1,j)+sshofl_w(i-1,j,1)))
         enddo ; enddo ; enddo
!       else                                !wwwww>
!        do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
!         anyv3d(i,j,k,1)=                                   &
!         (h_u(i,j)+0.5*(sshofl_w(i-1,j,1)+sshofl_w(i,j,1))) &
!         *0.5*(dsig_f(i,j,k)+dsig_f(i,j+1,k))
!        enddo ; enddo ; enddo
!       endif                               !wwwww>
       endif                    !sssss>

! Pour ne pas diviser par dz=quasizero  on moyenne separement les
! couches fusionnees: !07-05-19
        do j=1,jmax ; do i=1,imax+1
!           sum1=0. ; sum2=0.
!           do k=1,kmerged_u(i,j)
!            sum1=sum1+anyv3d(i,j,k,1)
!            sum2=sum2+velofl_u(i,j,k,1)
!           enddo
!           x1=sum2/max(sum1*dy_u(i,j),small1)
!           do k=1,kmerged_u(i,j)
!            anyvar3d(i,j,k)=x1
!           enddo
!           do k=kmerged_u(i,j)+1,kmax 
            do k=1,kmax !07-03-24 
             anyvar3d(i,j,k)=                                         &
             velofl_u(i,j,k,1)                                        &
                /dy_u(i,j)                                            &
              /anyv3d(i,j,k,1)
            enddo
        enddo ; enddo

        do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
          anyvar3d(i,j,k)=inv_scalefactor*(-var_addoffset+ &
          min(max(anyvar3D(i,j,k),var_validmin),var_validmax))
        enddo ; enddo ; enddo

      endif                  !=======>
      texte80(1)='u'    ; texte80(2)='m s-1'                          ! variable ; units
      texte80(3:4)='sea_water_x_velocity_at_u_location'
      texte80(5)='TZYX' ; texte80(7)=ofl_type
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offline_bin('_u') !08-06-18

!.....
! write V
      filval=-32767. ! Pour calculer var_scalefactor puis sera mis à zéro
      if(ofl_type=='short') then !oooooooo>
       var_validmin=-3. ; var_validmax=3.
       var_addoffset=0.5*(var_validmin+var_validmax)
       var_scalefactor=abs(var_validmax-var_addoffset)/(abs(filval)-1.)
      else                       !oooooooo>
       var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
      endif                      !oooooooo>
      inv_scalefactor=1./var_scalefactor
      filval=0.

      if(loop_netcdf==1) then !=======>

       if(ale_selected==1) then !sssss>
        do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
         anyv3d(i,j,k,1)=0.5*(dzofl_t(i,j,k,1)+dzofl_t(i,j-1,k,1))
        enddo ; enddo ; enddo
       else                     !sssss>
!       if(fgrid_or_wgrid==wgrid_case) then !wwwww>
         do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
          anyv3d(i,j,k,1)=0.5*(                          &
          dsig_t(i,j-1,k)*(h_w(i,j-1)+sshofl_w(i,j-1,1)) &
         +dsig_t(i,j  ,k)*(h_w(i,j  )+sshofl_w(i,j  ,1)) )
         enddo ; enddo ; enddo
!       else                                !wwwww>
!        do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
!         anyv3d(i,j,k,1)=                                   &
!         (h_v(i,j)+0.5*(sshofl_w(i,j-1,1)+sshofl_w(i,j,1))) &
!         *0.5*(dsig_f(i,j,k)+dsig_f(i+1,j,k))
!        enddo ; enddo ; enddo
!       endif                               !wwwww>
       endif                    !sssss>

! Pour ne pas diviser par dz=quasizero  on moyenne separement les
! couches fusionnees: !07-05-19
       do j=1,jmax+1 ; do i=1,imax
!           sum1=0. ; sum2=0.
!           do k=1,kmerged_v(i,j)
!            sum1=sum1+anyv3d(i,j,k,1)
!            sum2=sum2+velofl_v(i,j,k,1)
!           enddo
!           x1=sum2/max(sum1*dx_v(i,j),small1)
!           do k=1,kmerged_v(i,j)
!            anyvar3d(i,j,k)=x1
!           enddo
!           do k=kmerged_v(i,j)+1,kmax 
            do k=1,kmax !07-03-24 
             anyvar3d(i,j,k)=                        &
             velofl_v(i,j,k,1)                       &
                /dx_v(i,j)                           &
              /anyv3d(i,j,k,1)
            enddo
       enddo ; enddo

       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
          anyvar3d(i,j,k)=inv_scalefactor*(-var_addoffset+ &
          min(max(anyvar3D(i,j,k),var_validmin),var_validmax))
       enddo ; enddo ; enddo

      endif                  !=======>
      texte80(1)='v'    ; texte80(2)='m s-1'                          ! variable ; units
      texte80(3:4)='sea_water_y_velocity_at_v_location'
      texte80(5)='TZYX' ; texte80(7)=ofl_type
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offline_bin('_v') !08-06-18

      if(ofl_rotation==1) then !>>>>>>>>>>>>>>>>>>>>>>>>
!.....
! write U rotated:
      filval=-32767. ! Pour calculer var_scalefactor puis sera mis à zéro
      if(ofl_type=='short') then !oooooooo>
       var_validmin=-3. ; var_validmax=3.
       var_addoffset=0.5*(var_validmin+var_validmax)
       var_scalefactor=abs(var_validmax-var_addoffset)/(abs(filval)-1.)
      else                       !oooooooo>
       var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
      endif                      !oooooooo>
      inv_scalefactor=1./var_scalefactor
      filval=0.

      if(loop_netcdf==1) then !=======>

       if(ale_selected==1) then !sssss>
        do k=1,kmax ; do j=1,jmax ; do i=1,imax
         anyv3d(i,j,k,1)=dzofl_t(i,j,k,1)
        enddo ; enddo ; enddo
       else                     !sssss>
        if(fgrid_or_wgrid==wgrid_case) then !wwwww>
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
          anyv3d(i,j,k,1)=dsig_t(i,j,k)*(h_w(i,j)+sshofl_w(i,j,1))
         enddo ; enddo ; enddo
        else                                !wwwww>
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
          anyv3d(i,j,k,1)=0.25*( dsig_f(i,j  ,k)+dsig_f(i+1,j  ,k)  &
                                +dsig_f(i,j+1,k)+dsig_f(i+1,j+1,k)) &
                              *(h_w(i,j)+sshofl_w(i,j,1))
         enddo ; enddo ; enddo
        endif                               !wwwww>
       endif                    !sssss>

       do k=1,kmax ; do j=1,jmax ; do i=1,imax

             anyvar3d(i,j,k)=0.5*(                        &
            ( velofl_u(i  ,j  ,k,1)                       &
             +velofl_u(i+1,j  ,k,1))*gridrotcos_t(i,j)    &
           +( velofl_v(i  ,j  ,k,1)                       &
             +velofl_v(i  ,j+1,k,1))*gridrotsin_t(i,j)    &
              )/dy_t(i,j)                                 &
             /anyv3d(i,j,k,1)

          anyvar3d(i,j,k)=inv_scalefactor*(-var_addoffset+ &
          min(max(anyvar3D(i,j,k),var_validmin),var_validmax))

       enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='u_we'    ; texte80(2)='m s-1'                     ! variable ; units
      texte80(3:4)='eastward_current_at_t_location'
      texte80(5)='TZYX' ; texte80(7)=ofl_type
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offline_bin('_t') !08-06-18

!.....
! write V rotated:
      filval=-32767. ! Pour calculer var_scalefactor puis sera mis à zéro
      if(ofl_type=='short') then !oooooooo>
       var_validmin=-3. ; var_validmax=3.
       var_addoffset=0.5*(var_validmin+var_validmax)
       var_scalefactor=abs(var_validmax-var_addoffset)/(abs(filval)-1.)
      else                       !oooooooo>
       var_validmin=-1.e10 ; var_validmax=1.e10 ; var_scalefactor=1. ; var_addoffset=0.
      endif                      !oooooooo>
      inv_scalefactor=1./var_scalefactor
      filval=0.

      if(loop_netcdf==1) then !=======>

       if(ale_selected==1) then !sssss>
        do k=1,kmax ; do j=1,jmax ; do i=1,imax
         anyv3d(i,j,k,1)=dzofl_t(i,j,k,1)
        enddo ; enddo ; enddo
       else                     !sssss>
        if(fgrid_or_wgrid==wgrid_case) then !wwwww>
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
          anyv3d(i,j,k,1)=dsig_t(i,j,k)*(h_w(i,j)+sshofl_w(i,j,1))
         enddo ; enddo ; enddo
        else                                !wwwww>
         do k=1,kmax ; do j=1,jmax ; do i=1,imax
          anyv3d(i,j,k,1)=0.25*( dsig_f(i,j  ,k)+dsig_f(i+1,j  ,k)  &
                                +dsig_f(i,j+1,k)+dsig_f(i+1,j+1,k)) &
                              *(h_w(i,j)+sshofl_w(i,j,1))
         enddo ; enddo ; enddo
        endif                               !wwwww>
       endif                    !sssss>

       do k=1,kmax ; do j=1,jmax ; do i=1,imax
             anyvar3d(i,j,k)=0.5*(                        &

           -( velofl_u(i  ,j  ,k,1)                       &
             +velofl_u(i+1,j  ,k,1))*gridrotsin_t(i,j)    &
           +( velofl_v(i  ,j  ,k,1)                       &
             +velofl_v(i  ,j+1,k,1))*gridrotcos_t(i,j)    &

              )/dx_t(i,j)                                 &
             /anyv3d(i,j,k,1)

          anyvar3d(i,j,k)=inv_scalefactor*(-var_addoffset+ &
          min(max(anyvar3D(i,j,k),var_validmin),var_validmax))
       enddo ; enddo ; enddo
      endif                  !=======>
      texte80(1)='v_sn'    ; texte80(2)='m s-1'   !21-03-17
      texte80(3:4)='northward_current_at_t_location'
      texte80(5)='TZYX' ; texte80(7)=ofl_type
      if(flag_offline_binary==0)call netcdf_main('_t')
      if(flag_offline_binary==1)call offline_bin('_t') !08-06-18

      endif                    !>>>>>>>>>>>>>>>>>>>>>>>>

      if(loop_netcdf==1) then !=======>
! note: division par sumarchive car cette variable n'est pas passe (pour
! convenance) par routine end_of_average
      do j=1,jmax
       do i=1,imax+1
       anyvar2d(i,j)=fluxbarsum_ofl_u(i,j)    &
           /(iteration2d_max_now*dy_u(i,j)*sumarchive)
       enddo
      enddo
      endif                  !=======>
      texte80(1)='uah'    ; texte80(2)='m2 s-1'                         ! variable ; units
      texte80(3:4)='sea_water_x_transport_at_u_location'
      texte80(5)='TYX'; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_u')
      if(flag_offline_binary==1)call offline_bin('_u') !08-06-18

      if(loop_netcdf==1) then !=======>
! note: division par sumarchive car cette variable n'est pas passe (pour
! convenance) par routine end_of_average
      do j=1,jmax+1
       do i=1,imax
       anyvar2d(i,j)=fluxbarsum_ofl_v(i,j)    &
           /(iteration2d_max_now*dx_v(i,j)*sumarchive)
       enddo
      enddo
      endif                  !=======>
      texte80(1)='vah'    ; texte80(2)='m2 s-1'                         ! variable ; units
      texte80(3:4)='sea_water_y_transport_at_v_location'
      texte80(5)='TYX'; texte80(7)='real'
      if(flag_offline_binary==0)call netcdf_main('_v')
      if(flag_offline_binary==1)call offline_bin('_v') !08-06-18

      if(flag_offline_binary==0) then !netcdf-case->
       if(loop_netcdf==0) then !>>>>>>>>>>>>>>>>>>>>
        call netcdf_general_attributes(ncid)
! Definition of variables: done.
        status=nfmpi_enddef(ncid)
       endif                  !>>>>>>>>>>>>>>>>>>>
       status=nfmpi_close(ncid)
      endif                           !netcdf-case->


      enddo  ! fin de boucle sur loop_netcdf

      end subroutine offline_write_var

!.............................................................

      subroutine offline_another_model_divergence
      implicit none
      double precision implicit_fric_,explicit_fric_
#ifdef synopsis
       subroutinetitle='offline_another_model_divergence'
       subroutinedescription= &
       'Adjust depth-averaged currents in order to restore nemo ssh'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Calcul des perturbations de vitesse moyenne de la pseudo surface libre
! pour equilibrage de la divergence

      if(cfl_reduce>0.9) then   !21-06-16
      write(6,*)'cfl_reduce=',cfl_reduce
      stop 'ERR offline_another_model_divergence cfl_reduce>0.9'
      endif

      x1=dte_lp*2.6d-3*1.d-3    ! dt*cd*moduleU
      implicit_fric_=1./(1.+x1) ! 1/(1+dt*cd*moduleU)
      explicit_fric_=1.-x1      ! 1-dt*cd*moduleU

! Equation pour u: du/dt=-gDssh'/Dx avec ssh' ecart entre la ssh deduite du calcul
! et de la ssh de nemo
      do j=1,jmax   ! Attention la parallelisation impose des boucles
      do i=1,imax+1 ! plus grandes pour la partie moveforward que la boucle
       velbar_u(i,j,1)=velbar_u(i,j,2) ! du calcul suivant
      enddo
      enddo
      do j=2,jmax-1
      do i=2,imax
       velbar_u(i,j,2)=(                                             &
!      velbar_u(i,j,1)*explicit_fric_                                &
       velbar_u(i,j,1)                                               &
                      -dte_lp*grav*(                                 &
                                    +ssh_int_w(i  ,j  ,1)            &
                                    -ssh_int_w(i-1,j  ,1)            &
                                     -sshobc_w(i  ,j  ,1)            &
                                     +sshobc_w(i-1,j  ,1))/dx_u(i,j) &
       )*mask_u(i,j,kmax+1)*implicit_fric_
      enddo
      enddo
! Equation pour v: dv/dt=-gDssh'/Dy avec ssh' ecart entre la ssh deduite du calcul
! et de la ssh de nemo
      do j=1,jmax+1 ! Attention la parallelisation impose des boucles
      do i=1,imax   ! plus grandes pour la partie moveforward que la boucle
       velbar_v(i,j,1)=velbar_v(i,j,2) ! du calcul suivant
      enddo
      enddo
      do j=2,jmax
      do i=2,imax-1
       velbar_v(i,j,2)=(                                             &
!      velbar_v(i,j,1)*explicit_fric_                                &
       velbar_v(i,j,1)                                               &
                      -dte_lp*grav*(                                 &
                                    +ssh_int_w(i  ,j  ,1)            &
                                    -ssh_int_w(i  ,j-1,1)            &
                                     -sshobc_w(i  ,j  ,1)            &
                                     +sshobc_w(i  ,j-1,1))/dy_v(i,j) &
       )*mask_v(i,j,kmax+1)*implicit_fric_
      enddo
      enddo

!#ifdef parallele
!      ub3=ubound(velbar_u) ; lb3=lbound(velbar_u)
!      call echange('x ',velbar_u,lb3,ub3,2) ! 2 pour echange 3eme arg = 2 ! C.L. i=imax+1 i=1 j=jmax j=1
!      ub3=ubound(velbar_v) ; lb3=lbound(velbar_v)
!      call echange('y ',velbar_v,lb3,ub3,2) ! 2 pour echange 3eme arg = 2 !C.L. i=imax+1 i=1 j=jmax j=1
!#endif
! Nouveaux echanges:
!     call offline_nemo_obcext2  ! C.L sur velbar(2) !25-09-14
      call obc_ext_velbar_mpi(2) ! echange x sur velbar_u(:,:,2) et y sur velbar_v(:,:,2) !15-07-14

! Ajouter a l'equation de la divergence la contribution de la divergence de l'anomalie
! de courant. Attention, comme le pas de temps des equations precedentes est plus petit
! (CFL) que le pas de temps adopte pour la divergence, on normalise le courant par le
! rapport des pas de temps. C'est ce courant normalise qui sera ajoute au courant NEMO

      x0=0.5*2.*dte_lp/dti_lp ! redimensionnement du courant pour etirer la
                              ! divergence "petit pas de temps" a une divergence
                              ! "grand pas de temps" pour conservation des traceurs
      do j=1,jmax+1
      do i=1,imax+1
       xy_u(i,j,0)=x0*dy_u(i,j)*h_u(i,j)*(velbar_u(i,j,1)+velbar_u(i,j,2))
       xy_v(i,j,0)=x0*dx_v(i,j)*h_v(i,j)*(velbar_v(i,j,1)+velbar_v(i,j,2))
      enddo
      enddo

! Dssh/Dt=-div(courant):
      call offline_nemo_obcext ! C.L sur xy_u(0) et xy_v(0) !29-04-14
!     sum1=small1
!     sum2=0.
!     sum3=0.
      do j=1,jmax
      do i=1,imax

       ssh_int_w(i,j,2)=                                      &
       ssh_int_w(i,j,2)                                       &
           -dti_lp                                            &
                  *( xy_u(i+1,j  ,0)-xy_u(i,j,0)              &
                    +xy_v(i  ,j+1,0)-xy_v(i,j,0) )/dxdy_t(i,j)

! commente car deplace apres l'appel A la subroutine
       hz_w(i,j,2)=h_w(i,j)+ssh_int_w(i,j,2)

!      sum1=sum1+mask_t(i,j,kmax)
!      sum2=sum2+mask_t(i,j,kmax)*(ssh_int_w(i,j,2)-sshobc_w(i,j,1))
!      sum3=sum3+mask_t(i,j,kmax)*(ssh_int_w(i,j,2)-sshobc_w(i,j,1))**2


      enddo
      enddo


!     i=imax/2 ; j=jmax/2
!     if(mask_t(i,j,kmax)==1) &
!     write(10+par%rank,*)real(elapsedtime_now)/86400.,ssh_int_w(i,j,2),sshobc_w(i,j,1)

!     call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,  & !#MPI
!                        mpi_sum,par%comm2d,ierr)
!     call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,  & !#MPI
!                        mpi_sum,par%comm2d,ierr)
!     call mpi_allreduce(sum3,sum3glb,1,mpi_double_precision,  & !#MPI
!                        mpi_sum,par%comm2d,ierr)

!     if(par%rank==0)                                 &
!       write(66,*)elapsedtime_now/86400.             &
!      ,sum2glb/sum1glb                               &
!      ,sqrt( sum3glb/sum1glb - (sum2glb/sum1glb)**2 )

!       write(10+par%rank,*)elapsedtime_now/86400.    &
!      ,sum2/sum1                               &
!      ,sqrt( sum3/sum1 - (sum2/sum1)**2 )


#ifdef bidon
! ON NE PASSE PLUS PAR CES LIGNES
      if(flag_nemoffline==1) then                       !--nemo-case->
! Cas où les fleuves sont introduits à la surface. La divergence du transport
! n'est pas équilibrée par une variation de la surface mais par un flux via
! omega, celui ci introduira les C.L. pour les traceurs
       do kr=1,nriver
        if(riverdir(kr)==0) then !00000>
         i=iriver(kr,1) ; j=jriver(kr,1)
          if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) then !mpimpimpi> !04-12-13
           ssh_int_w(i,j,2)=                       & !02-04-13
           ssh_int_w(i,j,2)-dti_lp*riverflux(kr,1)
                hz_w(i,j,2)=h_w(i,j)+ssh_int_w(i,j,2)
          endif                                              !mpimpimpi>
        endif                    !00000>
       enddo
      endif                                             !--nemo-case->
#ifdef checkmpi
      call offline_check_mpi_conservation !13-01-11
#endif
#endif


! AJUSTER LE COURANT U
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,1)=small1
      enddo
      enddo
      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,1)=xy_u(i,j,1)+mask_u(i,j,k)*abs(veldydz_u(i,j,k,1))
      enddo
      enddo
      enddo
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,2)=xy_u(i,j,0)     &
                  /xy_u(i,j,1)
      enddo
      enddo
      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1
            veldydz_u(i,j,k,1)=     &
            veldydz_u(i,j,k,1)+     &
               xy_u(i,j,2)          &
              *mask_u(i,j,k)        &
       *abs(veldydz_u(i,j,k,1))
      enddo
      enddo
      enddo

! AJUSTER LE COURANT V
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,1)=small1
      enddo
      enddo
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,1)=xy_v(i,j,1)+mask_v(i,j,k)*abs(veldxdz_v(i,j,k,1))
      enddo
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,2)=xy_v(i,j,0)     &
                  /xy_v(i,j,1)
      enddo
      enddo
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
            veldxdz_v(i,j,k,1)=     &
            veldxdz_v(i,j,k,1)+     &
                xy_v(i,j,2)         &
              *mask_v(i,j,k)        &
       *abs(veldxdz_v(i,j,k,1))
      enddo
      enddo
      enddo

#ifdef bidon
! Verification:
      if(iteration3d==10) then
      do j=1,jmax+1
      do i=1,imax+1
       xy_u(i,j,1)=0.
       xy_v(i,j,1)=0.
       do k=1,kmax
        xy_u(i,j,1)=xy_u(i,j,1)+veldydz_u(i,j,k,1)
        xy_v(i,j,1)=xy_v(i,j,1)+veldxdz_v(i,j,k,1)
       enddo
      enddo
      enddo
      do j=1,jmax
      do i=1,imax
      if(mask_t(i,j,kmaxp1)==1) then
       write(10+par%rank,*)   &
       (ssh_int_w(i,j,2)-ssh_int_w(i,j,0))/dti_lp &
      +(xy_u(i+1,j,1)-xy_u(i,j,1)+xy_v(i,j+1,1)-xy_v(i,j,1))/dxdy_t(i,j) &
       ,-omega_w(i,j,kmaxp1,1)
      endif
      enddo
      enddo
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)
#endif
      stop 'bibi2'
      endif
#endif

      end subroutine offline_another_model_divergence

!------------------------------------------------------------------------------------

      subroutine offline_check_mpi_conservation
      implicit none
#ifdef synopsis
       subroutinetitle='offline_check_mpi_conservation'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Pour verifier la conservation de la parallelisation:
      j1=1 ; j2=jmax
      do i=1,imax
       xy_t(i,j1,1)=ssh_int_w(i,j1,2)
       xy_t(i,j2,1)=ssh_int_w(i,j2,2)
      enddo
      i1=1 ; i2=imax
      do j=1,jmax
       xy_t(i1,j,1)=ssh_int_w(i1,j,2)
       xy_t(i2,j,1)=ssh_int_w(i2,j,2)
      enddo
!#ifdef parallele
!      lb3=lbound(ssh_int_w) ; ub3=ubound(ssh_int_w)
!      call echange('z0',ssh_int_w,lb3,ub3,2) ! 2 pour echange 3eme arg = 2 !C.L. i=imax+1 i=0 j=jmax+1 j=0
!#endif
! Nouveaux echanges !15-07-14
      call obc_ssh_int_mpi(2) ! echange z0 sur ssh_int_w(:,:,2)
      ksecu=0
      j1=1 ; j2=jmax

      if(par%tvoisin(sud)/=mpi_proc_null) then !sssssssssss>
      do i=1,imax
       if(xy_t(i,j1,1)/=ssh_int_w(i,j1,2))then
        ksecu=1
        write(10+par%rank,*)'mpi non conservé en i j =',i,j1,'du proc',par%rank, &
         'voisin sud',par%tvoisin(sud),'coordonnées',i,jmax-1, &
         'ssh',xy_t(i,j1,1),ssh_int_w(i,j1,2),xy_t(i,j1,1)-ssh_int_w(i,j1,2)
       endif
      enddo
      endif                        !sssssssssss>


      if(par%tvoisin(nord)/=mpi_proc_null) then !nnnnnnnnnnn>
      do i=1,imax
       if(xy_t(i,j2,1)/=ssh_int_w(i,j2,2))then
        ksecu=1
        write(10+par%rank,*)'mpi non conservé en i j =',i,j2,'du proc',par%rank, &
         'voisin nord',par%tvoisin(nord),'coordonnées',i,2, &
         'ssh',xy_t(i,j2,1),ssh_int_w(i,j2,2),xy_t(i,j2,1)-ssh_int_w(i,j2,2)
       endif
      enddo
      endif                         !nnnnnnnnnnn>


      i1=1 ; i2=imax

      if(par%tvoisin(ouest)/=mpi_proc_null) then !ooooooooooo>
      do j=1,jmax
       if(xy_t(i1,j,1)/=ssh_int_w(i1,j,2)) then
        ksecu=1
        write(10+par%rank,*)'mpi non conservé en i j =',i1,j,'du proc',par%rank, &
          'voisin ouest',par%tvoisin(ouest),'coordonnées',imax-1,j, &
          'ssh',xy_t(i1,j,1),ssh_int_w(i1,j,2),xy_t(i1,j,1)-ssh_int_w(i1,j,2)
       endif
      enddo
      endif                          !ooooooooooo>

      if(par%tvoisin(est)/=mpi_proc_null) then   !eeeeeeeeeee>
      do j=1,jmax
       if(xy_t(i2,j,1)/=ssh_int_w(i2,j,2)) then
        ksecu=1
        write(*,*)'mpi non conservé en i j =',i2,j,'du proc',par%rank, &
          'voisin Est',par%tvoisin(est),'coordonnées',2,j, &
          'ssh',xy_t(i2,j,1),ssh_int_w(i2,j,2),xy_t(i2,j,1)-ssh_int_w(i2,j,2)
       endif
      enddo
      endif                          !eeeeeeeeeee>

      if(ksecu==1)then
      write(*,*)'Erreurs mpi notifiees dans fichiers fort locaux'
        stop 'STOP dans offline_check_mpi_conservation'
      endif

      end subroutine offline_check_mpi_conservation

!............................................................................................

      subroutine offline_nemo_obcext
      implicit none
#ifdef synopsis
       subroutinetitle='offline_nemo_obcext'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! C.L. = inversion de l'equation de continuité suivante:
      if(par%tvoisin(nord ) == mpi_proc_null) then !----->
       j=jmax
       x0=dte_lp/dti_lp
       do i=1,imax !11-07-14
          xy_v(i,j+1,0)= x0*(                                   & !20-07-14
!                            ssh_ext_w(i,j,1)                   &
!                            ssh_int_w(i,j,2)                   &
                        0.5*(ssh_int_w(i,j,2)+ssh_int_w(i,j,1)) &
                             -sshobc_w(i,j,1))                  &
                        *sqrt(grav*h_w(i,j))                    &
         *dy_v(i,j+1)                                           &
       *mask_v(i,j+1,kmaxp1)
       enddo
      endif                                        !----->

      if(par%tvoisin(sud  ) == mpi_proc_null) then !----->
       j=1
       x0=dte_lp/dti_lp
       do i=1,imax !11-07-14
          xy_v(i,j,0)=-x0*(                                   &
!                          ssh_ext_w(i,j,1)                   &
!                          ssh_int_w(i,j,2)                   &
                      0.5*(ssh_int_w(i,j,2)+ssh_int_w(i,j,1)) &
                           -sshobc_w(i,j,1))                  &
                      *sqrt(grav*h_w(i,j))                    &
         *dy_v(i,j)                                           &
       *mask_v(i,j,kmaxp1)
       enddo
      endif                                        !----->

      if(par%tvoisin(est  ) == mpi_proc_null) then !----->
       i=imax
       x0=dte_lp/dti_lp
       do j=1,jmax
          xy_u(i+1,j,0)= x0*(                                   & !20-07-14
!                            ssh_ext_w(i,j,1)                   &
!                            ssh_int_w(i,j,2)                   &
                        0.5*(ssh_int_w(i,j,2)+ssh_int_w(i,j,1)) &
                             -sshobc_w(i,j,1))                  &
                        *sqrt(grav*h_w(i,j))                    &
         *dy_u(i+1,j)                                           &
       *mask_u(i+1,j,kmaxp1)
       enddo
      endif                                        !----->

      if(par%tvoisin(ouest) == mpi_proc_null) then !----->
       i=1
       x0=dte_lp/dti_lp
       do j=1,jmax
          xy_u(i,j,0)=-x0*(                                   &
!                          ssh_ext_w(i,j,1)                   &
!                          ssh_int_w(i,j,2)                   &
                      0.5*(ssh_int_w(i,j,2)+ssh_int_w(i,j,1)) &
                           -sshobc_w(i,j,1))                  &
                      *sqrt(grav*h_w(i,j))                    &
         *dy_u(i,j)                                           &
       *mask_u(i,j,kmaxp1)
       enddo
      endif                                        !----->


      end subroutine offline_nemo_obcext

!..........................................................................

      subroutine offline_nemo_obcext2 !25-09-14
      implicit none
#ifdef synopsis
       subroutinetitle='offline_nemo_obcext2'
       subroutinedescription= &
        'Radiative boundary conditions for the 2D momentum' &
       ' equations used to restore the nemo ssh'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       x0=grav*dte_lp/86400. ! 24h restoring time scale

! C.L. = inversion de l'equation de continuité suivante:
      if(par%tvoisin(nord ) == mpi_proc_null) then !----->
       j=jmax
       do i=1,imax !11-07-14
           velbar_v(i,j+1,2)=(velbar_v(i,j+1,1)+x0*(              &
                        0.5*(ssh_int_w(i,j,2)+ssh_int_w(i,j,1))   &
                             -sshobc_w(i,j,1) )/dy_v(i,j+1)       &
                                            )*mask_v(i,j+1,kmaxp1)
       enddo
      endif                                        !----->

      if(par%tvoisin(sud  ) == mpi_proc_null) then !----->
       j=1
       do i=1,imax !11-07-14
             velbar_v(i,j,2)=(velbar_v(i,j,1)-x0*(              &
                        0.5*(ssh_int_w(i,j,2)+ssh_int_w(i,j,1)) &
                             -sshobc_w(i,j,1) )/dy_v(i,j)       &
                                            )*mask_v(i,j,kmaxp1)

       enddo
      endif                                        !----->

      if(par%tvoisin(est  ) == mpi_proc_null) then !----->
       i=imax
       do j=1,jmax
           velbar_u(i+1,j,2)=(velbar_u(i+1,j,1)+x0*(              &
                        0.5*(ssh_int_w(i,j,2)+ssh_int_w(i,j,1))   &
                             -sshobc_w(i,j,1) )/dx_u(i+1,j)       &
                                            )*mask_u(i+1,j,kmaxp1)
       enddo
      endif                                        !----->

      if(par%tvoisin(ouest) == mpi_proc_null) then !----->
       i=1
       do j=1,jmax
             velbar_u(i,j,2)=(velbar_u(i,j,1)-x0*(                &
                        0.5*(ssh_int_w(i,j,2)+ssh_int_w(i,j,1))   &
                             -sshobc_w(i,j,1) )/dx_u(i,j)       &
                                            )*mask_u(i,j,kmaxp1)

       enddo

      endif                                        !----->


      end subroutine offline_nemo_obcext2

!..........................................................................

      subroutine offline_kz_convect_adjust !27-09-14
#ifdef synopsis
       subroutinetitle='offline_kz_convect_adjust'
       subroutinedescription= &
        'In case of arbiratry limitation on kz in NEMO-MERCATOR,'   &
       ' restores a big value on kz in case of neutral or instable' &
       ' vertical gradient of potential density'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Computes sea water potential density from temobc(2) and salobc(2)
        if(eos_tkezref>=0.) then
         call equation_of_state_potzref_jmfwg('obc',2) !23-02-22
        else
         call equation_of_state_potloc_jmfwg('obc',2)
        endif

      do k=1,kmax ; do j=1,jmax ; do i=1,imax
       rhp_t(i,j,k)=anyv3d(i,j,k,1)
      enddo       ; enddo       ; enddo

! Case method 1 or 2: Applies a large value (here 10) to kz if drhp/dz>=0.
! with x1 an error bar on drhp/dz
!     x1=-0.00001 ! error bar on drhp/dz
      x1=0.

! Case method 2 or 3: applies a large value (here 10) to kz if kz>=x2
      x2=0.0099

      do k=2,kmax
      do j=1,jmax
      do i=1,imax

! Methode 1 uniquement basee sur signe de drhp/dz
!      if( (  rhp_t(i,j,k)  -rhp_t(i,j,k-1))      &
!         /(depth_t(i,j,k)-depth_t(i,j,k-1))>=x1)dfvofl_w(i,j,k,0)=10.

! Methode 2 basee sur signe de drhp/dz et sur seuil Kz nemo:
       if(dfvofl_w(i,j,k,0)>=x2) then !-------->
        if( (  rhp_t(i,j,k)  -rhp_t(i,j,k-1))      &
           /(depth_t(i,j,k)-depth_t(i,j,k-1))>=x1)dfvofl_w(i,j,k,0)=10.
       endif                          !-------->

! Methode 3 basee uniquement sur seuil Kz nemo:
!      if(dfvofl_w(i,j,k,0)>=x2)dfvofl_w(i,j,k,0)=10.

      enddo
      enddo
      enddo

      end subroutine offline_kz_convect_adjust

!..........................................................................

      subroutine offline_read_file  !16-04-15
      implicit none
#ifdef synopsis
       subroutinetitle='offline_read_file'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Get ssh (egalement ici qu'on lit la prochaine periode et prochain rdv)
      call offline_read_file_ssh ! et aussi ofl_readtime_next et ofl_period_next

! Get temperature & salinity
      call offline_read_file_ts

! Ici aiguillage entre 2 cas de modèles:
      if(flag_nemoffline==0) then                       !---Aiguillage--->
           call offline_read_file_dz
      else                                              !---Aiguillage--->           !04-11-12
! Cas modele NEMO: (fichier de grille. Inconvenient dz_u et dz_v de l'etat initial)
           call offline_read_file_e1e2e3
      endif                                             !---Aiguillage--->           !04-11-12

! Get vertical mixing coef:
      call offline_read_file_kz

! Get surface vertical velocity:
      if(flag_nemoffline==1) then
         call offline_read_file_w_nemo
      else
         call offline_read_file_w
      endif

! Get convective surface layer vertical index !09-02-19
      if(flag_ksloffline==1)call offline_read_file_ksl !09-02-19

      if(flag_maxbotstress==1)call offline_read_file_botstress !17-11-17

! Get velocity
      call offline_read_file_vel

! Reversed time mode: change fluxes et velocity sign !26-05-19
      if(ofl_reversedtime==1) then !m°v°m>
            velobc_u(:,:,:,2)=-velobc_u(:,:,:,2)
            velobc_v(:,:,:,2)=-velobc_v(:,:,:,2)
            w0mofl_w(:,:,0)  =-w0mofl_w(:,:,0)
        w_keq1_ofl_w(:,:,0)  =-w_keq1_ofl_w(:,:,0)
       if(flag_nemoffline==1)riverflux(:,2)=-riverflux(:,2)
      endif                        !m°v°m>

      end subroutine offline_read_file

!------------------------------------------------------------------------

      subroutine offline_read_file_ssh
      implicit none
      integer nc_

      nc_=nc
! Phase initial du cas symphonie: lire 2 champs instantanes consecutifs
      if(flag_nemoffline==0.and.offline_init_status==0)nc_=nc-1

! Si le code plante ici c'est que l'on demarre trop pres de l'etat initial
! Avoir en effet en tete que la lecture de ssh instantannee necessite 1 champs de plus
      if(nc_<1)stop 'Err nc_<1 offline_read_file_ssh'

! Get ssh:
 3794 continue
      if(par%rank==0)write(6,*)'Get SSH with nc=',nc_
      open(unit=3,file=trim(tmpdirname)//'liste_offline_ssh.binrec'           &
                 ,access='direct',recl=recordlength_,form='unformatted')
      read(3,rec=nc_)filename_,ogcmtimecounter_,ofl_readtime_next,ofl_period_next
      close(3)
      if(par%rank==0) then  !.....
       write(6,'(a,a)')'offline fichier: ',trim(filename_)
       write(6,*)'ogcmtimecounter_=',ogcmtimecounter_
       write(6,*)'ofl_readtime_next in days=',ofl_readtime_next/86400.
       write(6,*)'ofl_period_next in hours=',ofl_period_next/3600.
      endif                 !.....

      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
      if(status/=0)stop ' stop offline_inout.f erreur 1'

      start(3)=ogcmtimecounter_  ; edge(3)=1        ! time
      start(2)=1+par%tjmax(1)    ; edge(2)=jmax+2   ! j
      start(1)=1+par%timax(1)    ; edge(1)=imax+2   ! i

      flag_newalgo=0 ! par defaut algo compatible avec anciens fichiers
                   status=nfmpi_inq_varid(ncid1,'ssh_inst',var_id)    !05-07-15
      if(status==0)flag_newalgo=1 ! Si =1 algo nouveaux fichiers active
      if(status/=0)status=nfmpi_inq_varid(ncid1,'ssh',var_id)   
      if(status/=0)status=nfmpi_inq_varid(ncid1,'ssh-bi',var_id)      !17-05-11
      if(status/=0)status=nfmpi_inq_varid(ncid1,'ssh-ib',var_id)      !06-04-12
      if(status/=0)status=nfmpi_inq_varid(ncid1,'sossheig',var_id)    !05-11-12
      if(status/=0)status=nfmpi_inq_varid(ncid1,'ssh_ib',var_id)      !03-07-13
      if(status/=0)stop ' stop offline_inout.f var_id ssh'

      status=nfmpi_get_att_real(ncid1,var_id,'scale_factor',scalarval(1))
      if(status/=0)scalarval(1)=1.
      var_scalefactor=scalarval(1)
                   status=nfmpi_get_att_real(ncid1,var_id,'add_offset',scalarval(1))
      if(status/=0)status=nfmpi_get_att_real(ncid1,var_id,'offset',scalarval(1))
      if(status/=0)scalarval(1)=0.
      var_addoffset=scalarval(1)
      status=nfmpi_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)

      if(var_type==nf_real)  then !rrrrrrrrr
      status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:3),edge(1:3) &!07-09-12
                                          ,anyvar2d(0:imax+1,0:jmax+1))
      endif                       !rrrrrrrrr

      if(var_type==nf_short) then !sssssssss
!       if(iteration3d_restart==0) then !08-05-18
         if(.not.allocated(oflshort))allocate(oflshort(0:imax+1,0:jmax+1,kmax+1))
!       endif
       status=nfmpi_get_vara_int_all(ncid1,var_id,start(1:3),edge(1:3) &
                                        ,oflshort(0:imax+1,0:jmax+1,1))
       do j=0,jmax+1 ; do i=0,imax+1
        anyvar2d(i,j)=oflshort(i,j,1)*var_scalefactor+var_addoffset
       enddo         ; enddo
      endif                       !sssssssss
      if(status/=0)stop ' stop offline_inout.f get ssh'

      status=nfmpi_close(ncid1)

      if(flag_nemoffline==1) then !-- Cas nemo offline -->
       do j=0,jmax+1 ; do i=0,imax+1
        sshobc_w(i,j,0)=sshobc_w(i,j,2)
        sshobc_w(i,j,2)=max(anyvar2d(i,j)*mask_t(i,j,kmax) &   !20-12-14
                                ,wetdry_cst2-h_w(i,j))
       enddo ; enddo
      else                        !-- Cas Symphonie -->
       do j=0,jmax+1 ; do i=0,imax+1
        sshobc_w(i,j,0)=sshobc_w(i,j,1)
        sshobc_w(i,j,1)=sshobc_w(i,j,2)
        sshobc_w(i,j,2)=max(anyvar2d(i,j)*mask_t(i,j,kmax) &   !20-12-14
                                ,wetdry_cst2-h_w(i,j))
        sshobc_w(i,j,0)=0.
        sshobc_w(i,j,1)=0.
        sshobc_w(i,j,2)=0.
       enddo ; enddo
      endif                       !--               -->

      nc_=nc_+1
      if(nc_==nc)goto 3794 ! Phase initial symphonie lecture 2 echeances

      end subroutine offline_read_file_ssh

!.........................................................................................

      subroutine offline_read_file_ts
      implicit none

! Get temperature
!     if(imodelbio==1) then !(provided that the biogeochemical model is running)!30-03-11

      if(par%rank==0)write(6,*)'Get temperature with nc=',nc
      open(unit=3,file=trim(tmpdirname)//'liste_offline_t.binrec'           &
                 ,access='direct',recl=recordlength_,form='unformatted')
      read(3,rec=nc)filename_,ogcmtimecounter_
      close(3)
      if(par%rank==0)write(6,'(a,a)')'offline fichier:',trim(filename_)
      if(par%rank==0)write(6,*)'ogcmtimecounter_=',ogcmtimecounter_

       start(4)=ogcmtimecounter_ ; edge(4)=1        ! time
       start(3)=1                ; edge(3)=kmax     ! k
       start(2)=1+par%tjmax(1)   ; edge(2)=jmax+2   ! j
       start(1)=1+par%timax(1)   ; edge(1)=imax+2   ! i

      if(vert_axis_conv_direc=='up') then  !05-11-12
       kbegin_=1 ; kend_=kmax ; kstep_=1
      endif
      if(vert_axis_conv_direc=='dw') then !05-11-12
       kend_=1 ; kbegin_=kmax ; kstep_=-1
      endif

      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
      if(status/=0)stop ' stop offline_inout.f erreur 2'

                    status=nfmpi_inq_varid(ncid1,'tem',var_id)
       if(status/=0)status=nfmpi_inq_varid(ncid1,'votemper',var_id)
       if(status/=0)stop ' stop offline_inout.f var_id tem'

       status=nfmpi_get_att_real(ncid1,var_id,'scale_factor',scalarval(1))
       if(status/=0)scalarval(1)=1.
       var_scalefactor=scalarval(1)
                    status=nfmpi_get_att_real(ncid1,var_id,'add_offset',scalarval(1))
       if(status/=0)status=nfmpi_get_att_real(ncid1,var_id,'offset',scalarval(1))
       if(status/=0)scalarval(1)=0.
       var_addoffset=scalarval(1)
       status=nfmpi_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)

      if(var_type==nf_real)  then !rrrrrrrrr
       status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:4),edge(1:4) &
                      ,anyvar3d(0:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
      endif                       !rrrrrrrrr

      if(var_type==nf_short) then !sssssssss
!       if(iteration3d_restart==0) then !08-05-18
         if(.not.allocated(oflshort))allocate(oflshort(0:imax+1,0:jmax+1,kmax+1))
!       endif
        status=nfmpi_get_vara_int_all(ncid1,var_id,start(1:4),edge(1:4) &
                      ,oflshort(0:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
       do k=kbegin_,kend_,kstep_  ; do j=0,jmax+1 ; do i=0,imax+1
        anyvar3d(i,j,k)=oflshort(i,j,k)*var_scalefactor+var_addoffset
       enddo         ; enddo         ; enddo
      endif                       !sssssssss


       if(status/=0)stop ' stop offline_inout.f get tem'

       do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1 !07-03-23
        temobc_t(i,j,k,0)=temobc_t(i,j,k,2)
        temobc_t(i,j,k,2)=anyvar3d(i,j,k)
       enddo ; enddo ; enddo

      status=nfmpi_close(ncid1)


! Get salinity:
!     if(imodelbio==1) then !(provided that the biogeochemical model is running)!30-03-11
      if(par%rank==0)write(6,*)'Get salinity with nc=',nc
      open(unit=3,file=trim(tmpdirname)//'liste_offline_s.binrec'           &
                 ,access='direct',recl=recordlength_,form='unformatted')
      read(3,rec=nc)filename_,ogcmtimecounter_
      close(3)
      if(par%rank==0)write(6,'(a,a)')'offline fichier:',trim(filename_)
      if(par%rank==0)write(6,*)'ogcmtimecounter_=',ogcmtimecounter_

       start(4)=ogcmtimecounter_ ; edge(4)=1        ! time
       start(3)=1                ; edge(3)=kmax     ! k
       start(2)=1+par%tjmax(1)   ; edge(2)=jmax+2   ! j
       start(1)=1+par%timax(1)   ; edge(1)=imax+2   ! i

      if(vert_axis_conv_direc=='up') then  !05-11-12
       kbegin_=1 ; kend_=kmax ; kstep_=1
      endif
      if(vert_axis_conv_direc=='dw') then !05-11-12
       kend_=1 ; kbegin_=kmax ; kstep_=-1
      endif

      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
      if(status/=0)stop ' stop offline_inout.f erreur 3'

                    status=nfmpi_inq_varid(ncid1,'sal',var_id)
       if(status/=0)status=nfmpi_inq_varid(ncid1,'vosaline',var_id)
       if(status/=0)stop ' stop offline_inout.f var_id sal'

       status=nfmpi_get_att_real(ncid1,var_id,'scale_factor',scalarval(1))
       if(status/=0)scalarval(1)=1.
       var_scalefactor=scalarval(1)
                    status=nfmpi_get_att_real(ncid1,var_id,'add_offset',scalarval(1))
       if(status/=0)status=nfmpi_get_att_real(ncid1,var_id,'offset',scalarval(1))
       if(status/=0)scalarval(1)=0.
       var_addoffset=scalarval(1)
       status=nfmpi_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)

      if(var_type==nf_real)  then !rrrrrrrrr
       status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:4),edge(1:4) &
                      ,anyvar3d(0:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
      endif                       !rrrrrrrrr
      if(var_type==nf_short) then !sssssssss
!       if(iteration3d_restart==0) then !08-05-18
         if(.not.allocated(oflshort))allocate(oflshort(0:imax+1,0:jmax+1,kmax+1))
!       endif
        status=nfmpi_get_vara_int_all(ncid1,var_id,start(1:4),edge(1:4) &
                      ,oflshort(0:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
       do k=kbegin_,kend_,kstep_  ; do j=0,jmax+1 ; do i=0,imax+1
        anyvar3d(i,j,k)=oflshort(i,j,k)*var_scalefactor+var_addoffset
       enddo         ; enddo         ; enddo
      endif                       !sssssssss
      if(status/=0)stop ' stop offline_inout.f get sal'

       do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1 !07-03-23
        salobc_t(i,j,k,0)=salobc_t(i,j,k,2)
        salobc_t(i,j,k,2)=anyvar3d(i,j,k)
       enddo ; enddo ; enddo

      status=nfmpi_close(ncid1)

      end subroutine offline_read_file_ts

!.........................................................................................

      subroutine offline_read_file_dz
      implicit none

! Cas modele S:
      open(unit=3,file=trim(tmpdirname)//'liste_offline_t.binrec'           &
                 ,access='direct',recl=recordlength_,form='unformatted')
      read(3,rec=nc)filename_,ogcmtimecounter_
      close(3)

      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
      if(status/=0)stop ' stop offline_inout.f erreur 4'

! Get dz_t:
       start(4)=ogcmtimecounter_ ; edge(4)=1        ! time
       start(3)=1                ; edge(3)=kmax     ! k
       start(2)=1+par%tjmax(1)   ; edge(2)=jmax+2   ! j
       start(1)=1+par%timax(1)   ; edge(1)=imax+2   ! i

      if(vert_axis_conv_direc=='up') then  !05-11-12
       kbegin_=1 ; kend_=kmax ; kstep_=1
      endif
      if(vert_axis_conv_direc=='dw') then !05-11-12
       kend_=1 ; kbegin_=kmax ; kstep_=-1
      endif

       status=nfmpi_inq_varid(ncid1,'dz',var_id)
       if(status==0) then !-dz_available->
       status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:4),edge(1:4) &
                     ,anyvar3d(0:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
       if(status/=0)stop ' stop offline_inout.f erreur get dz'
       status=nfmpi_close(ncid1)
       endif              !-dz_available->

       if(status/=0) then !-dz_not-available-use-dsig_t->

       status=nfmpi_close(ncid1)

       status=nfmpi_open(par%comm2d                             &
!                        ,trim(directory_offline)//'grid.nc'  & !21-07-14
                         ,trim(lonlatfile)  & !22-05-18
                         ,nf_nowrite,MPI_INFO_NULL,ncid1)
       if(status/=0)stop ' stop offline_inout.f erreur 2609'

       status=nfmpi_inq_varid(ncid1,'dsig_t',var_id)
       if(status/=0)stop 'error nfmpi_inq_varid dsig_t'
       status=nfmpi_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)
       if(var_type/=nf_double)stop 'error var_type dsig_t'

       start(3)=1 ; edge(3)=kmax
       status=nfmpi_get_vara_double_all(ncid1,var_id            &
                                             ,start(1:var_dims) &
                                              ,edge(1:var_dims) &
                     ,anyv3d(0:imax+1,0:jmax+1,1:kmax,1))
       if(status/=0)stop ' error nfmpi_get_vara_double_all dsig_t'
        do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
         anyvar3d(i,j,k)=anyv3d(i,j,k,1)*(h_w(i,j)+sshobc_w(i,j,2))
        enddo ; enddo ; enddo

       status=nfmpi_close(ncid1)

       endif              !-dz_not-available-use-dsig_t->

      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax+1
! Pour NE PAS INTERFERER AVEC le dz des calculs, on passe par un tableau anyv3d
        anyv3d(i,j,k,1)=0.5*(anyvar3d(i-1,j  ,k)+anyvar3d(i,j,k))  ! dz_u
        anyv3d(i,j,k,2)=0.5*(anyvar3d(i  ,j-1,k)+anyvar3d(i,j,k))  ! dz_v
      enddo
      enddo
      enddo

      end subroutine offline_read_file_dz

!.........................................................................................

      subroutine offline_read_file_e1e2e3
      implicit none

      status=nfmpi_open(par%comm2d,initialgridfile_txt,nf_nowrite,MPI_INFO_NULL,ncid1)
      if(status/=0)stop 'erreur ouverture fichier initialgridfile_txt'
       if(hori_axis_conv_start=='v') &
       stop 'Cas hori_axis_conv_start==0 non prevu'

! Get e3u:
      status=nfmpi_inq_varid(ncid1,'e3u',var_id)
      if(status==0)then !uuuuuuuuuuuuuuuuuuuuuuuu>

       if(hori_axis_conv_start=='t'.and.hori_axis_conv_end=='v') then !--->
        start(4)=1              ; edge(4)=1        ! time
        start(3)=1              ; edge(3)=kmax     ! k
        start(2)=1+par%tjmax(1) ; edge(2)=jmax+2   ! j
        start(1)=1+par%timax(1) ; edge(1)=imax+1   ! i
        i1=1 ; i2=imax+1
       else                                                           !--->
        stop 'cas non prevu 5011'
       endif

       if(vert_axis_conv_direc=='up') then  !05-11-12
        kbegin_=1 ; kend_=kmax ; kstep_=1
       endif
       if(vert_axis_conv_direc=='dw') then !05-11-12
        kend_=1 ; kbegin_=kmax ; kstep_=-1
       endif

       status=nfmpi_inq_vartype(ncid1,var_id,type_)
       if(type_==6) then !>>>
        status=                                                     &
        nfmpi_get_vara_double_all(ncid1,var_id,start(1:4),edge(1:4) &
                ,anyv3d(i1:i2,0:jmax+1,kbegin_:kend_:kstep_,1))
       else              !>>>
        stop 'offline_inout cas e3u non double non prevu'
       endif             !>>>
      if(status/=0)stop 'offline_inout erreur lecture e3u'

      else              !uuuuuuuuuuuuuuuuuuuuuuuu>

       stop 'e3u non trouve dans le fichier initialgridfile_txt'

       do k=1,kmax ; do j=0,jmax+1 ; do i=1,imax+1
         anyv3d(i,j,k,1)=dz_u(i,j,k,1) ! valeur initiale !30-04-14
       enddo ; enddo ; enddo


! Ces lignes ont été commentées à cause de la retro-action instable
! de la surface libre sur e3u et e3v
!      do k=1,kmax ; do j=0,jmax+1 ; do i=1,imax+1
!        anyv3d(i,j,k,1)=0.5*(dz_t(i,j,k,1)+dz_t(i-1,j,k,1)) !=e3u
!      enddo ; enddo ; enddo
! Partial step e3u deduit du partial step de e3t
!      do j=0,jmax+1 ; do i=1,imax+1
!        anyv3d(i,j,kmin_u(i,j),1)=min(dz_t(i  ,j,kmin_u(i,j),1)      &
!                                     ,dz_t(i-1,j,kmin_u(i,j),1))
!      enddo ; enddo

!      do k=1,kmax ; do j=0,jmax+1 ; do i=1,imax+1
!       write(50+par%rank,*)anyv3d(i,j,k,1),dz_u(i,j,k,1)
!      enddo ; enddo ; enddo


      endif             !uuuuuuuuuuuuuuuuuuuuuuuu>


! Get e3v:
      status=nfmpi_inq_varid(ncid1,'e3v',var_id)
      if(status==0)then !vvvvvvvvvvvvvvvvvvvvvvvvv>

       if(hori_axis_conv_start=='t'.and.hori_axis_conv_end=='v') then !--->
        start(4)=1              ; edge(4)=1        ! time
        start(3)=1              ; edge(3)=kmax     ! k
        start(2)=1+par%tjmax(1) ; edge(2)=jmax+1   ! j
        start(1)=1+par%timax(1) ; edge(1)=imax+2   ! i
        j1=1 ; j2=jmax+1
       else                                                           !--->
        stop 'cas non prevu 5012'
       endif

       if(vert_axis_conv_direc=='up') then  !05-11-12
        kbegin_=1 ; kend_=kmax ; kstep_=1
       endif
       if(vert_axis_conv_direc=='dw') then !05-11-12
        kend_=1 ; kbegin_=kmax ; kstep_=-1
       endif

       status=nfmpi_inq_vartype(ncid1,var_id,type_)
       if(type_==6) then !>>>
        status=                                                     &
        nfmpi_get_vara_double_all(ncid1,var_id,start(1:4),edge(1:4) &
                ,anyv3d(0:imax+1,j1:j2,kbegin_:kend_:kstep_,2))
       else              !>>>
        stop 'offline_inout cas e3v non double non prevu'
       endif             !>>>
      if(status/=0)stop 'offline_inout erreur lecture e3u'

      else              !vvvvvvvvvvvvvvvvvvvvvvvvvvv>

       stop 'e3v non trouve dans le fichier initialgridfile_txt'

        do k=1,kmax ; do j=1,jmax+1 ; do i=0,imax+1
         anyv3d(i,j,k,2)=dz_v(i,j,k,1) ! valeur initiale !30-04-14
       enddo ; enddo ; enddo
!      do k=1,kmax ; do j=1,jmax+1 ; do i=0,imax+1
!        anyv3d(i,j,k,2)=0.5*(dz_t(i,j,k,1)+dz_t(i,j-1,k,1)) !=e3v
!      enddo ; enddo ; enddo
! Partial step e3v deduit du partial step de e3t
!      do j=1,jmax+1 ; do i=0,imax+1
!        anyv3d(i,j,kmin_v(i,j),2)=min(dz_t(i,j  ,kmin_v(i,j),1)      &
!                                     ,dz_t(i,j-1,kmin_v(i,j),1))
!      enddo ; enddo

!       do k=1,kmax ; do j=1,jmax+1 ; do i=0,imax+1
!        write(60+par%rank,*)anyv3d(i,j,k,2),dz_v(i,j,k,1)
!       enddo ; enddo ; enddo

!      if(par%rank==0) then
!       k=1 ; j=1 ; i=0
!       write(77,*)'module_offline',dz_t(i,j,k,1),dz_t(i,j-1,k,1)
!      endif

!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)
!#endif
!     stop 'coucou'


      endif             !vvvvvvvvvvvvvvvvvvvvvvvvvvv>


! Get e2u:
      if(hori_axis_conv_start=='t'.and.hori_axis_conv_end=='v') then !--->
       start(3)=1              ; edge(3)=1        ! time
       start(2)=1+par%tjmax(1) ; edge(2)=jmax+2   ! j
       start(1)=1+par%timax(1) ; edge(1)=imax+1   ! i
       i1=1 ; i2=imax+1
      else                                                           !--->
       stop 'cas non prevu 3011'
      endif                                                          !--->
      status=nfmpi_inq_varid(ncid1,'e2u',var_id)
      if(status/=0)stop ' stop offline_inout.f erreur var_id e2u'
      status=nfmpi_inq_vartype(ncid1,var_id,type_)
      if(type_==6) then !>>>
       status=                                                     &
       nfmpi_get_vara_double_all(ncid1,var_id,start(1:3),edge(1:3) &
                                ,dy_u(i1:i2,0:jmax+1))
!                               ,xy_u(i1:i2,0:jmax+1,1))
      else              !>>>
       stop 'offline_inout cas e2u non double non prevu'
      endif             !>>>
      if(status/=0)stop 'offline_inout erreur lecture e2u'

! Get e1v:
      if(hori_axis_conv_start=='t'.and.hori_axis_conv_end=='v') then !--->
       start(3)=1              ; edge(3)=1        ! time
       start(2)=1+par%tjmax(1) ; edge(2)=jmax+1   ! j
       start(1)=1+par%timax(1) ; edge(1)=imax+2   ! i
       j1=1 ; j2=jmax+1
      else                                                           !--->
       stop 'cas non prevu 3044'
      endif                                                          !--->
      status=nfmpi_inq_varid(ncid1,'e1v',var_id)
      if(status/=0)stop ' stop offline_inout.f erreur var_id e1v'
      status=nfmpi_inq_vartype(ncid1,var_id,type_)
      if(type_==6) then !>>>
       status=                                                     &
       nfmpi_get_vara_double_all(ncid1,var_id,start(1:3),edge(1:3) &
                                ,dx_v(0:imax+1,j1:j2))
!                               ,xy_v(0:imax+1,j1:j2,1))
      else              !>>>
       stop 'offline_inout cas e1v non double non prevu'
      endif             !>>>
      if(status/=0)stop 'offline_inout erreur lecture e1v'

      status=nfmpi_close(ncid1)

      end subroutine offline_read_file_e1e2e3

!.........................................................................................

      subroutine offline_read_file_kz
      implicit none
      integer :: logscale_=0

! Get vertical mixing coef:
      if(par%rank==0)write(6,*)'Get Kz with nc=',nc
      open(unit=3,file=trim(tmpdirname)//'liste_offline_kz.binrec'           &
                 ,access='direct',recl=recordlength_,form='unformatted')
      read(3,rec=nc)filename_,ogcmtimecounter_
      close(3)
      if(par%rank==0)write(6,'(a,a)')'offline fichier kz: ',trim(filename_)
      if(par%rank==0)write(6,*)'ogcmtimecounter_=',ogcmtimecounter_

       start(4)=ogcmtimecounter_ ; edge(4)=1        ! time
       start(3)=1                ; edge(3)=kmax+1   ! k
       start(2)=1+par%tjmax(1)   ; edge(2)=jmax+2   ! j
       start(1)=1+par%timax(1)   ; edge(1)=imax+2   ! i

! CAS NEMO: autant de point t que de point Kz:
      if(vert_axis_conv_end/=vert_axis_conv_start)edge(3)=kmax

      if(vert_axis_conv_direc=='up') then  !05-11-12
       kbegin_=1 ; kend_=kmax+1 ; kstep_=1
      endif
      if(vert_axis_conv_direc=='dw') then !05-11-12
       kbegin_=kmax+1 ; kstep_=-1 ; kend_=kbegin_-edge(3)+1 ! kend_=2 si cas NEMO
      endif

      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
      if(status/=0)stop ' stop offline_inout.f erreur 5'
                   status=nfmpi_inq_varid(ncid1,'difv',var_id)
      if(status/=0)status=nfmpi_inq_varid(ncid1,'votkeavt',var_id)
      if(status/=0) then !logloglog?
                        status=nfmpi_inq_varid(ncid1,'log10difv',var_id)
           if(status/=0)status=nfmpi_inq_varid(ncid1,'vokzln10',var_id)
                    if(status==0)logscale_=1
      endif              !logloglog?
      if(status/=0)stop ' stop offline_inout.f erreur var_id difv'

       status=nfmpi_get_att_real(ncid1,var_id,'scale_factor',scalarval(1))
       if(status/=0)scalarval(1)=1.
       var_scalefactor=scalarval(1)
                    status=nfmpi_get_att_real(ncid1,var_id,'add_offset',scalarval(1))
       if(status/=0)status=nfmpi_get_att_real(ncid1,var_id,'offset',scalarval(1))
       if(status/=0)scalarval(1)=0.
       var_addoffset=scalarval(1)
       status=nfmpi_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)

      if(var_type==nf_real) then  !rrrrrrrrrr
                   status=nfmpi_get_att_real(ncid1,var_id,'_FillValue',scalarval(1))
      if(status/=0)status=nfmpi_get_att_real(ncid1,var_id,'_Fillvalue',scalarval(1)) !20-03-15
      if(status/=0)stop ' module_offline err 3969 _FillValue not found'
       filval=scalarval(1)
       status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:4),edge(1:4) &
                      ,anyvar3d(0:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
       if(logscale_==1) then !llllllll> 13-03-15
        do k=kbegin_,kend_ ,kstep_ ; do j=0,jmax+1 ; do i=0,imax+1
         if(anyvar3d(i,j,k)==filval) then !--->
            anyvar3d(i,j,k)=1.e-6
         else                             !--->
            anyvar3d(i,j,k)=10**anyvar3d(i,j,k)
         endif                            !--->
        enddo         ; enddo         ; enddo
       else                  !llllllll>
        do k=kbegin_,kend_ ,kstep_ ; do j=0,jmax+1 ; do i=0,imax+1
         if(anyvar3d(i,j,k)==filval)anyvar3d(i,j,k)=1.e-6
        enddo         ; enddo         ; enddo
       endif                 !llllllll>
      endif                       !rrrrrrrrrr
      if(var_type==nf_short) then !ssssssssss
                   status=nfmpi_get_att_int(ncid1,var_id,'_FillValue',i1val(1))
      if(status/=0)stop ' module_offline err 3980 _FillValue not found'
      filvalshort=i1val(1)
        status=nfmpi_get_vara_int_all(ncid1,var_id,start(1:4),edge(1:4) &
                      ,oflshort(0:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
       do k=kbegin_,kend_ ,kstep_ ; do j=0,jmax+1 ; do i=0,imax+1
        anyvar3d(i,j,k)=oflshort(i,j,k)*var_scalefactor+var_addoffset
                     if(oflshort(i,j,k)==filvalshort)anyvar3d(i,j,k)=0.
       enddo         ; enddo         ; enddo
       if(logscale_==1) then !llllllll> 13-03-15
        do k=kbegin_,kend_ ,kstep_ ; do j=0,jmax+1 ; do i=0,imax+1
            anyvar3d(i,j,k)=10**anyvar3d(i,j,k)
         if(oflshort(i,j,k)==filvalshort)anyvar3d(i,j,k)=0.
        enddo         ; enddo         ; enddo
       endif                 !llllllll>
      endif                       !ssssssssss

      if(status/=0)stop ' stop offline_inout.f erreur get difv'
      status=nfmpi_close(ncid1)

      do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
       dfvofl_w(i,j,k,-1)=dfvofl_w(i,j,k,0)
       dfvofl_w(i,j,k,0) =anyvar3d(i,j,k)
      enddo ; enddo ; enddo

! CAS NEMO: autant de point t que de point Kz:
      if(flag_nemoffline==1) then !111>     ! cas nemo offline
        if(vert_axis_conv_end/=vert_axis_conv_start) then !222>
                       if(vert_axis_conv_direc=='dw') then  !333>
                                     dfvofl_w(:,:,1,0)=dfvofl_w(:,:,2,0)
                                                      endif !333>
                                                    endif  !222>
      endif                       !111>


      if(ofl_tke==1) then !tke-tke-tke-> !05-03-15
! Get tke: !05-03-15
      if(par%rank==0)write(6,*)'Get tke with nc=',nc
      open(unit=3,file=trim(tmpdirname)//'liste_offline_tke.binrec'           &
                 ,access='direct',recl=recordlength_,form='unformatted')
      read(3,rec=nc)filename_,ogcmtimecounter_
      close(3)
      if(par%rank==0)write(6,'(a,a)')'offline fichier tke: ',trim(filename_)
      if(par%rank==0)write(6,*)'ogcmtimecounter_=',ogcmtimecounter_

       start(4)=ogcmtimecounter_ ; edge(4)=1        ! time
       start(3)=1                ; edge(3)=kmax+1   ! k
       start(2)=1+par%tjmax(1)   ; edge(2)=jmax+2   ! j
       start(1)=1+par%timax(1)   ; edge(1)=imax+2   ! i

! CAS NEMO: autant de point t que de point Kz:
      if(vert_axis_conv_end/=vert_axis_conv_start)edge(3)=kmax

      if(vert_axis_conv_direc=='up') then  !05-11-12
       kbegin_=1 ; kend_=kmax+1 ; kstep_=1
      endif
      if(vert_axis_conv_direc=='dw') then !05-11-12
       kbegin_=kmax+1 ; kstep_=-1 ; kend_=kbegin_-edge(3)+1 ! kend_=2 si cas NEMO
      endif

      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
      if(status/=0)stop ' stop offline_inout.f erreur 5'
                   status=nfmpi_inq_varid(ncid1,'tke',var_id)
      if(status/=0) then !logloglog?
                        status=nfmpi_inq_varid(ncid1,'log10tke',var_id)
                    if(status==0)logscale_=1
      endif              !logloglog?
      if(status/=0)stop ' stop offline_inout.f erreur var_id tke'

       status=nfmpi_get_att_real(ncid1,var_id,'scale_factor',scalarval(1))
       if(status/=0)scalarval(1)=1.
       var_scalefactor=scalarval(1)
                    status=nfmpi_get_att_real(ncid1,var_id,'add_offset',scalarval(1))
       if(status/=0)status=nfmpi_get_att_real(ncid1,var_id,'offset',scalarval(1))
       if(status/=0)scalarval(1)=0.
       var_addoffset=scalarval(1)
       status=nfmpi_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)

      if(var_type==nf_real) then  !rrrrrrrrrr
       status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:4),edge(1:4) &
                      ,anyvar3d(0:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
      endif                       !rrrrrrrrrr
      if(var_type==nf_short) then !ssssssssss
        status=nfmpi_get_vara_int_all(ncid1,var_id,start(1:4),edge(1:4) &
                      ,oflshort(0:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
       do k=kbegin_,kend_ ,kstep_ ; do j=0,jmax+1 ; do i=0,imax+1
        anyvar3d(i,j,k)=oflshort(i,j,k)*var_scalefactor+var_addoffset
        anyvar3d(i,j,k)=10**anyvar3d(i,j,k)
       enddo         ; enddo         ; enddo
       if(logscale_==0)stop 'offline Kz must be en log10 scale'
      endif                       !ssssssssss

      if(status/=0)stop ' stop offline_inout.f erreur get difv'
      status=nfmpi_close(ncid1)

      do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
       tkeofl_w(i,j,k,-1)=tkeofl_w(i,j,k,0)
       tkeofl_w(i,j,k,0) =anyvar3d(i,j,k)
      enddo ; enddo ; enddo

! CAS NEMO: autant de point t que de point Kz:
      if(flag_nemoffline==1) then !111>   ! cas nemo offline
        if(vert_axis_conv_end/=vert_axis_conv_start) then !222>
                       if(vert_axis_conv_direc=='dw') then  !333>
                                     tkeofl_w(:,:,1,0)=tkeofl_w(:,:,2,0)
                                                      endif !333>
                                                    endif !222>
      endif                       !111>
      endif               !tke-tke-tke-> !05-03-15

      end subroutine offline_read_file_kz

!.........................................................................................

      subroutine offline_read_file_w_nemo
      implicit none

       if(.not.allocated(w0mofl_w)) then
        allocate(w0mofl_w(0:imax+1,0:jmax+1,-1:0))
        w0mofl_w=0. ! attention on ne met à zéro qu'au premier passage !!!!
       endif

! wobc_w n'est pas necessaire. Je l'utilise A l'occasion pour verifier
! que la vitesse verticale diagnostiquee par S26 est bien la mEme que
! celle archivee dans NEMO. Decommenter les lignes correspondantes pour
! realiser ce test
!      if(.not.allocated(wobc_w)) then
!       allocate(wobc_w(0:imax+1,0:jmax+1,1:kmax+1,-1:0))
!       wobc_w=0.
!      endif

! ATTENTION: la grille NEMO compte 75 niveaux, que l'on considere U,V ou
! W. Les fichiers ne tiennent donc pas compte de ce qu'il n'y a en
! principe pas le meme nombre de niveaux pour W et les autres variables.
! Autrement dit, quand on numerote A l'envers (cas 'dw') on part de la
! surface numerotee kmax+1 dans la convention S26 pour descendre au
! niveau k=2 (et non pas k=1) car seuls kmax niveaux ont ete archives
! dans le fichier NEMO
       if(vert_axis_conv_direc=='up') then  !04-07-15
        kbegin_=1 ; kend_=kmax ; kstep_=1
       endif
       if(vert_axis_conv_direc=='dw') then   !04-07-15
        kend_=2 ; kbegin_=kmax+1 ; kstep_=-1 !21-06-16
       endif

      if(par%rank==0)write(6,*)'Get W with nc=',nc
      open(unit=3,file=trim(tmpdirname)//'liste_offline_w.binrec'           & !11-03-15
                 ,access='direct',recl=recordlength_,form='unformatted')
      read(3,rec=nc)filename_,ogcmtimecounter_
      close(3)
      if(par%rank==0)write(6,'(a,a)')'offline fichier W: ',trim(filename_)
      if(par%rank==0)write(6,*)'ogcmtimecounter_=',ogcmtimecounter_

! filename_ est le même que celui de la variable Kz
      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
      status=nfmpi_inq_varid(ncid1,'vovecrtz',var_id)
      if(status==0)then !000000000>

       status=nfmpi_get_att_real(ncid1,var_id,'scale_factor',scalarval(1))
       if(status/=0)scalarval(1)=1.
       var_scalefactor=scalarval(1)
                    status=nfmpi_get_att_real(ncid1,var_id,'add_offset',scalarval(1))
       if(status/=0)status=nfmpi_get_att_real(ncid1,var_id,'offset',scalarval(1))
       if(status/=0)scalarval(1)=0.
       var_addoffset=scalarval(1)
       status=nfmpi_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)

       if(var_type==nf_real) then  !rrrrrrrrrr
        status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:4),edge(1:4) &
                      ,anyvar3d(0:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
       endif                       !rrrrrrrrrr
       if(var_type==nf_short) then !ssssssssss
        status=nfmpi_get_vara_int_all(ncid1,var_id,start(1:4),edge(1:4) &
                      ,oflshort(0:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
        do k=kbegin_,kend_ ,kstep_ ; do j=0,jmax+1 ; do i=0,imax+1
         anyvar3d(i,j,k)=oflshort(i,j,k)*var_scalefactor+var_addoffset
        enddo         ; enddo         ; enddo
       endif                       !ssssssssss

       if(status/=0)stop ' stop offline_inout.f erreur get w'

       do j=1,jmax ; do i=1,imax
        w0mofl_w(i,j,-1)=w0mofl_w(i,j,0)
        w0mofl_w(i,j,0)= anyvar3d(i,j,kmaxp1)*mask_t(i,j,kmaxp1)
       enddo ; enddo
! wobc_w n'est pas necessaire. Je l'utilise A l'occasion pour verifier
! que la vitesse verticale diagnostiquee par S26 est bien la mEme que
! celle archivee dans NEMO. Decommenter les lignes correspondantes pour
! realiser ce test
!      do k=2,kmax+1 ; do j=1,jmax ; do i=1,imax
!       wobc_w(i,j,k,-1)=wobc_w(i,j,k,0)
!       wobc_w(i,j,k,0)= anyvar3d(i,j,k)*mask_t(i,j,kmaxp1)
!      enddo ; enddo ; enddo

      else              !000000000>
       stop ' Vertical Velocity not found in Nemo file' !25-07-14
      endif             !000000000>
      status=nfmpi_close(ncid1) ! attention le close doit bien etre ici...

      end subroutine offline_read_file_w_nemo

!------------------------------------------------------------------------

      subroutine offline_read_file_w
      implicit none

! Get surface vertical velocity
       if(.not.allocated(w0mofl_w)) then
        allocate(w0mofl_w(0:imax+1,0:jmax+1,-1:0))
        w0mofl_w=0. ! attention on ne met à zéro qu'au premier passage !!!!
       endif

      if(par%rank==0)write(6,*)'Get omega(kmax) with nc=',nc
      open(unit=3,file=trim(tmpdirname)//'liste_offline_w.binrec'     &
                 ,access='direct',recl=recordlength_,form='unformatted')
      read(3,rec=nc)filename_,ogcmtimecounter_,ofl_readtime_next,ofl_period_next
      close(3)
      if(par%rank==0) then  !.....
       write(6,'(a,a)')'offline fichier: ',trim(filename_)
       write(6,*)'ogcmtimecounter_=',ogcmtimecounter_
       write(6,*)'ofl_readtime_next in days=',ofl_readtime_next/86400.
       write(6,*)'ofl_period_next in hours=',ofl_period_next/3600.
      endif                 !.....

      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
      if(status/=0)stop ' stop offline_inout.f erreur 1'

      start(3)=ogcmtimecounter_  ; edge(3)=1        ! time
      start(2)=1+par%tjmax(1)    ; edge(2)=jmax+2   ! j
      start(1)=1+par%timax(1)    ; edge(1)=imax+2   ! i

                   status=nfmpi_inq_varid(ncid1,'w0m',var_id)
      if(status/=0)stop ' stop offline_inout.f var_id w0m'

      status=nfmpi_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(status/=0)stop ' Err 4454 nfmpi_inq_var'
      if(var_type==nf_real)  then !rrrrrrrrr
      status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:3),edge(1:3) &
                                          ,anyvar2d(0:imax+1,0:jmax+1))
      else                        !rrrrrrrrr
       stop 'Err w0m type/=nf_real'
      endif                       !rrrrrrrrr

! Avant de fermer detecter si sources sous marine presentes dans le fichier
         status=nfmpi_inq_varid(ncid1,'wbottom',var_id)
      if(status/=0) then !grw> !11-04-19
           flag_groundwater=0
      else               !grw> !11-04-19
           flag_groundwater=1
      endif              !grw> !11-04-19


      status=nfmpi_close(ncid1)

      do j=0,jmax+1
      do i=0,imax+1
       w0mofl_w(i,j,-1)=w0mofl_w(i,j,0)
       w0mofl_w(i,j,0) =anyvar2d(i,j)
      enddo
      enddo

! Get bottom vertical velocity
      if(flag_groundwater==1) then !GRWGRWGRW> 11-04-19

       if(.not.allocated(w_keq1_ofl_w)) then
        allocate(w_keq1_ofl_w(0:imax+1,0:jmax+1,-1:0))
        w_keq1_ofl_w=0. ! attention on ne met à zéro qu'au premier passage !!!!
       endif

      if(par%rank==0)write(6,*)'Get omega(k=1) with nc=',nc
      open(unit=3,file=trim(tmpdirname)//'liste_offline_w.binrec'     &
                 ,access='direct',recl=recordlength_,form='unformatted')
      read(3,rec=nc)filename_,ogcmtimecounter_,ofl_readtime_next,ofl_period_next
      close(3)
      if(par%rank==0) then  !.....
       write(6,'(a,a)')'offline fichier: ',trim(filename_)
       write(6,*)'ogcmtimecounter_=',ogcmtimecounter_
       write(6,*)'ofl_readtime_next in days=',ofl_readtime_next/86400.
       write(6,*)'ofl_period_next in hours=',ofl_period_next/3600.
      endif                 !.....

      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
      if(status/=0)stop ' stop offline_inout.f erreur 1'

      start(3)=ogcmtimecounter_  ; edge(3)=1        ! time
      start(2)=1+par%tjmax(1)    ; edge(2)=jmax+2   ! j
      start(1)=1+par%timax(1)    ; edge(1)=imax+2   ! i

                   status=nfmpi_inq_varid(ncid1,'wbottom',var_id)
      if(status/=0)stop ' stop offline_inout.f var_id wbottom'

      status=nfmpi_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(status/=0)stop ' Err 4454 nfmpi_inq_var'
      if(var_type==nf_real)  then !rrrrrrrrr
      status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:3),edge(1:3) &
                                          ,anyvar2d(0:imax+1,0:jmax+1))
      else                        !rrrrrrrrr
       stop 'Err wbottom type/=nf_real'
      endif                       !rrrrrrrrr

      status=nfmpi_close(ncid1)

      do j=1,jmax ; do i=1,imax
       w_keq1_ofl_w(i,j,-1)=w_keq1_ofl_w(i,j,0)
       w_keq1_ofl_w(i,j,0) =anyvar2d(i,j)
      enddo ; enddo

      endif                        !GRWGRWGRW> 11-04-19

      end subroutine offline_read_file_w

!.........................................................................................

      subroutine offline_read_file_ksl !09-02-19
      implicit none

! Get convective surface layer vertical index
       if(.not.allocated(kslofl_t)) then
        allocate(kslofl_t(0:imax+1,0:jmax+1,-1:0))
        kslofl_t=real(kmax)
       endif

      if(par%rank==0)write(6,*)'Get ksl with nc=',nc
      open(unit=3,file=trim(tmpdirname)//'liste_offline_ksl.binrec'     &
                 ,access='direct',recl=recordlength_,form='unformatted')
      read(3,rec=nc)filename_,ogcmtimecounter_,ofl_readtime_next,ofl_period_next
      close(3)
      if(par%rank==0) then  !.....
       write(6,'(a,a)')'offline fichier: ',trim(filename_)
       write(6,*)'ogcmtimecounter_=',ogcmtimecounter_
       write(6,*)'ofl_readtime_next in days=',ofl_readtime_next/86400.
       write(6,*)'ofl_period_next in hours=',ofl_period_next/3600.
      endif                 !.....

      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
      if(status/=0)stop ' stop offline_inout.f erreur 1'

      start(3)=ogcmtimecounter_  ; edge(3)=1        ! time
      start(2)=1+par%tjmax(1)    ; edge(2)=jmax+2   ! j
      start(1)=1+par%timax(1)    ; edge(1)=imax+2   ! i

                   status=nfmpi_inq_varid(ncid1,'ksl',var_id)
      if(status/=0)stop ' stop offline_inout.f var_id ksl'

      status=nfmpi_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(status/=0)stop ' Err 14454 nfmpi_inq_var'
      if(var_type==nf_real)  then !rrrrrrrrr
      status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:3),edge(1:3) &
                                          ,anyvar2d(0:imax+1,0:jmax+1))
      else                        !rrrrrrrrr
       stop 'Err ksl type/=nf_real'
      endif                       !rrrrrrrrr

      status=nfmpi_close(ncid1)

      do j=0,jmax+1
      do i=0,imax+1
       kslofl_t(i,j,-1)=kslofl_t(i,j,0)
       kslofl_t(i,j,0) =anyvar2d(i,j)
      enddo
      enddo

      end subroutine offline_read_file_ksl

!.........................................................................................

      subroutine offline_read_file_vel
      implicit none

! https://docs.google.com/document/d/1qCWQrn6EmHbOk1mlVIslo71GlCpBNx43LA9_wWrP1p8/edit

! Get u velocity: (and apply scale factors)
      if(par%rank==0)write(6,*)'Get u velocity with nc=',nc
      open(unit=3,file=trim(tmpdirname)//'liste_offline_u.binrec'           &
                 ,access='direct',recl=recordlength_,form='unformatted')
      read(3,rec=nc)filename_,ogcmtimecounter_
      close(3)
      if(par%rank==0)write(6,'(a,a)')'offline fichier u: ',trim(filename_)
      if(par%rank==0)write(6,*)'ogcmtimecounter_=',ogcmtimecounter_

       start(4)=ogcmtimecounter_ ; edge(4)=1        ! time
       start(3)=1                ; edge(3)=kmax     ! k
       start(2)=1+par%tjmax(1)   ; edge(2)=jmax+2   ! j
       start(1)=1+par%timax(1)   ; edge(1)=imax+1   ! i
      if(vert_axis_conv_direc=='up') then  !05-11-12
       kbegin_=1 ; kend_=kmax ; kstep_=1
      endif
      if(vert_axis_conv_direc=='dw') then !05-11-12
       kbegin_=kmax ; kend_=1 ; kstep_=-1
      endif

      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
      if(status/=0)stop ' stop offline_inout.f erreur 6'
                   status=nfmpi_inq_varid(ncid1,'u',var_id)
      if(status/=0)status=nfmpi_inq_varid(ncid1,'vozocrtx',var_id)
      if(status/=0)stop ' stop offline_inout.f erreur var_id u'

       status=nfmpi_get_att_real(ncid1,var_id,'scale_factor',scalarval(1))
       if(status/=0)scalarval(1)=1.
       var_scalefactor=scalarval(1)
                    status=nfmpi_get_att_real(ncid1,var_id,'add_offset',scalarval(1))
       if(status/=0)status=nfmpi_get_att_real(ncid1,var_id,'offset',scalarval(1))
       if(status/=0)scalarval(1)=0.
       var_addoffset=scalarval(1)
       status=nfmpi_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)

      if(var_type==nf_real) then  !rrrrrrrrrr
      status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:4),edge(1:4) &
                     ,anyvar3d(1:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
      endif                       !rrrrrrrrrr
      if(var_type==nf_short) then !ssssssssss
        status=nfmpi_get_vara_int_all(ncid1,var_id,start(1:4),edge(1:4) &
                      ,oflshort(1:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
        do k=kbegin_,kend_ ,kstep_ ; do j=0,jmax+1 ; do i=1,imax+1
         anyvar3d(i,j,k)=oflshort(i,j,k)*var_scalefactor+var_addoffset
        enddo         ; enddo         ; enddo
      endif                       !ssssssssss

      if(status/=0)stop ' stop offline_inout.f erreur get u'
      status=nfmpi_close(ncid1)

! On n'applique pas le masque pour S en raison des points fleuves mais on
! le fait pour nemo car on ne presume pas de son fillvalue (eventuellement /=0)
!     if(index(nomfichier(20),'another_model')/=0) then ! cas nemo offline
      if(flag_nemoffline==1) then                       ! cas nemo offline
       do k=1,kmax ; do j=1,jmax ; do i=1,imax+1
        anyvar3d(i,j,k)=anyvar3d(i,j,k)*mask_u(i,j,k)
       enddo       ; enddo       ; enddo
      endif                                             ! cas nemo offline

      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1

           velobc_u(i,j,k,0)=velobc_u(i,j,k,2)
           velobc_u(i,j,k,2)=                                          &
           anyvar3d(i,j,k)                                             &
              *dy_u(i,j)                                               &
            *anyv3d(i,j,k,1)

      enddo
      enddo
      enddo

! Get v velocity: (and apply scale factors)
      open(unit=3,file=trim(tmpdirname)//'liste_offline_v.binrec'           &
                 ,access='direct',recl=recordlength_,form='unformatted')
      read(3,rec=nc)filename_,ogcmtimecounter_
      close(3)
      if(par%rank==0)write(6,'(a,a)')'offline fichier v: ',trim(filename_)
      if(par%rank==0)write(6,*)'ogcmtimecounter_=',ogcmtimecounter_

      start(4)=ogcmtimecounter_ ; edge(4)=1        ! time
      start(3)=1                ; edge(3)=kmax     ! k
      start(2)=1+par%tjmax(1)   ; edge(2)=jmax+1   ! j
      start(1)=1+par%timax(1)   ; edge(1)=imax+2   ! i

      if(vert_axis_conv_direc=='up') then  !05-11-12
       kbegin_=1 ; kend_=kmax ; kstep_=1
      endif
      if(vert_axis_conv_direc=='dw') then !05-11-12
       kbegin_=kmax ; kend_=1 ; kstep_=-1
      endif

      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
      if(status/=0)stop ' stop offline_inout.f erreur 7'

                   status=nfmpi_inq_varid(ncid1,'v',var_id)
      if(status/=0)status=nfmpi_inq_varid(ncid1,'vomecrty',var_id)
      if(status/=0)stop ' stop offline_inout.f erreur var_id v'

       status=nfmpi_get_att_real(ncid1,var_id,'scale_factor',scalarval(1))
       if(status/=0)scalarval(1)=1.
       var_scalefactor=scalarval(1)
                    status=nfmpi_get_att_real(ncid1,var_id,'add_offset',scalarval(1))
       if(status/=0)status=nfmpi_get_att_real(ncid1,var_id,'offset',scalarval(1))
       if(status/=0)scalarval(1)=0.
       var_addoffset=scalarval(1)
       status=nfmpi_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)

      if(var_type==nf_real) then  !rrrrrrrrrr
      status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:4),edge(1:4) &
                     ,anyvar3d(0:imax+1,1:jmax+1,kbegin_:kend_:kstep_))
      endif                       !rrrrrrrrrr
      if(var_type==nf_short) then !ssssssssss
      status=nfmpi_get_vara_int_all(ncid1,var_id,start(1:4),edge(1:4) &
                    ,oflshort(0:imax+1,1:jmax+1,kbegin_:kend_:kstep_))
        do k=kbegin_,kend_ ,kstep_ ; do j=1,jmax+1 ; do i=0,imax+1
         anyvar3d(i,j,k)=oflshort(i,j,k)*var_scalefactor+var_addoffset
        enddo         ; enddo         ; enddo
      endif                       !ssssssssss

      if(status/=0)stop ' stop offline_inout.f erreur get v'
      status=nfmpi_close(ncid1)

! On n'applique pas le masque pour S en raison des points fleuves mais on
! le fait pour nemo car on ne presume pas de son fillvalue (eventuellement /=0)
!     if(index(nomfichier(20),'another_model')/=0) then ! cas nemo offline
      if(flag_nemoffline==1) then                       ! cas nemo offline
       do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax
!       if(mask_v(i,j,k)==1.and.abs(anyvar3d(i,j,k))>9.)stop 'bug v'
!        if(par%rank==0.and.i==6) then
!         if(j==69.or.j==70) then
!          write(6,*)'V=',k,mask_v(i,j,k),anyvar3d(i,j,k)
!         endif
!        endif
        anyvar3d(i,j,k)=anyvar3d(i,j,k)*mask_v(i,j,k)
       enddo       ; enddo       ; enddo
      endif                                             ! cas nemo offline

      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
          velobc_v(i,j,k,0)=velobc_v(i,j,k,2)
          velobc_v(i,j,k,2)=                                          &
          anyvar3d(i,j,k)                                             &
             *dx_v(i,j)                                               &
            *anyv3d(i,j,k,2)
      enddo
      enddo
      enddo

!     if(par%rank==4) then
!      i=imax/2 ; j=jmax/2
!      do k=kmax,1,-1
!       write(6,*)k,mask_v(i,j,k),velobc_v(i,j,k,2)/anyv3d(i,j,k,1)/dx_v(i,j),anyvar3d(i,j,k)
!      enddo
!     stop 'Verifier la grille verticale pour v nemo'
!     endif

! CAS NEMO: lire les flux des rivieres à l'etat initial pour les localiser:
!     if(index(nomfichier(20),'another_model')/=0   &
      if(flag_nemoffline==1.and.iteration3d==0) then !---runoff---runoff--->

! Get runoff:
      open(unit=3,file=trim(tmpdirname)//'liste_offline_ssh.binrec'           &
                 ,access='direct',recl=recordlength_,form='unformatted')
      read(3,rec=nc)filename_,ogcmtimecounter_
      close(3)
      if(par%rank==0)write(6,'(a,a)')'offline fichier: ',trim(filename_)
      if(par%rank==0)write(6,*)'ogcmtimecounter_=',ogcmtimecounter_

      start(3)=ogcmtimecounter_ ; edge(3)=1        ! time
      start(2)=1+par%tjmax(1)   ; edge(2)=jmax+2   ! j
      start(1)=1+par%timax(1)   ; edge(1)=imax+2   ! i

      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
      if(status/=0)stop ' stop offline_inout.f erreur 8'

                   status=nfmpi_inq_varid(ncid1,'sorunoff',var_id)
!     if(status/=0)stop ' stop offline_inout.f var_id sorunoff'
      if(status/=0) then !>>>>>>>>
         status=nfmpi_close(ncid1)
         goto 3374
      endif              !>>>>>>>>

      status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:3),edge(1:3) &!07-09-12
                                          ,anyvar2d(0:imax+1,0:jmax+1))
      if(status/=0)stop ' stop offline_inout.f get ssh'
      status=nfmpi_close(ncid1)

! 2 Solutions pour les rivieres de NEMO. Soit le notebook_river ignore
! l'existence des runoff nemo (nriver=0) soit le notebook_river connait
! l'existence des rivieres (nriver/=0)


      if(nriver==0) then !z z z z z z z z z>
       do j=0,jmax+1
       do i=0,imax+1
        if(anyvar2d(i,j)/=0.)then !rrrrrrrr>
         nriver=nriver+1
         iriver(nriver,1)=i
         jriver(nriver,1)=j
         riverdir(nriver)=0
         rivertrc_inout(nriver)=1
! Attention aux unites, des m/s pour etre homogene a omega
         riverflux(nriver,0)=riverflux(nriver,2)
         riverflux(nriver,2)=-abs(anyvar2d(i,j))*0.001 ! EN M/S (=omega) !02-04-13
        endif                     !rrrrrrrr>
       enddo
       enddo
      else               !z z z z z z z z z> !07-04-13

       do kr=1,nriver
        if(rivertrc_inout(kr)==1) then !******>
          i=iriver(kr,1) ; j=jriver(kr,1)
          riverflux(kr,0)=riverflux(kr,2)
          riverflux(kr,2)=-abs(anyvar2d(i,j))*0.001 ! EN M/S (=omega) !02-04-13
          if(anyvar2d(i,j)==0.) then !---->
           write(6,*)'Stop no nemo runoff at location given' &
                    ,' in notebook_rivers'
          endif                      !---->
        endif                       !******>
       enddo

      endif              !z z z z z z z z z>

 3374 continue

      endif                                             !---runoff---runoff--->


! CAS NEMO: Appliquer une correction sur le Kz de nemo-mercator en situation de convection
      if(     flag_nemoffline==1 &
         .and.flag_kz_enhanced==1)call offline_kz_convect_adjust !26-09-14
      if(flag_nemoffline==1)return

! Get transport:
      if(par%rank==0)write(6,*)'Get transport with nc=',nc
      open(unit=3,file=trim(tmpdirname)//'liste_offline_u.binrec'           &
                 ,access='direct',recl=recordlength_,form='unformatted')
      read(3,rec=nc)filename_,ogcmtimecounter_
      close(3)
      if(par%rank==0)write(6,'(a,a)')'offline fichier transport: ',trim(filename_)
      if(par%rank==0)write(6,*)'ogcmtimecounter_=',ogcmtimecounter_

      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
      if(status/=0)stop ' stop offline_inout.f erreur 6'

! U:
      status=nfmpi_inq_varid(ncid1,'uah',var_id)
      if(status/=0)stop ' stop offline_inout.f erreur var_id uah'
      start(3)=ogcmtimecounter_ ; edge(3)=1        ! time
      start(2)=1+par%tjmax(1)   ; edge(2)=jmax+2   ! j
      start(1)=1+par%timax(1)   ; edge(1)=imax+1   ! i
      status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:3),edge(1:3) &
                                          ,anyvar2d(1:imax+1,0:jmax+1))
      if(status/=0)stop ' stop offline_inout.f erreur get uah'

      do j=0,jmax+1 ; do i=1,imax+1
       xy_u(i,j,1)=anyvar2d(i,j)*dy_u(i,j)
      enddo ; enddo

! V:
      status=nfmpi_inq_varid(ncid1,'vah',var_id)
      if(status/=0)stop ' stop offline_inout.f erreur var_id vah'
      start(3)=ogcmtimecounter_ ; edge(3)=1        ! time
      start(2)=1+par%tjmax(1)   ; edge(2)=jmax+1   ! j
      start(1)=1+par%timax(1)   ; edge(1)=imax+2   ! i
      status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:3),edge(1:3) &
                                          ,anyvar2d(0:imax+1,1:jmax+1))

      do j=1,jmax+1 ; do i=0,imax+1
       xy_v(i,j,1)=anyvar2d(i,j)*dx_v(i,j)
      enddo ; enddo

      status=nfmpi_close(ncid1)

! Assurer que l'integrale de flux3d est bien egal au flux2d
      do j=0,jmax+1 ; do i=1,imax+1
! Avant 05-10-19 le courant est annulE dans les couches ecrasEes. Par la suite
! on prefere simplement ne pas corriger le flux barotrope dans ces couches (pour ne pas risquer une imprecision) et en mEme temps garder
! un courant non nul pour ne pas creer une zone de courant nulle dans laquelle viendraient se pieger des drifters (qui seraient alors ! immobiles)
! Avant 05-10-19:
! Pour eviter pb de precision on force le flux des couches ecrasees A zero !08-05-19
!      do k=1,kmin_u(i,j)-1
!       velobc_u(i,j,k,2)=0.
!      enddo
! On ajuste le transport sur les couches restantes
! Ne connaissant pas dz_u a priori, on pondere par la valeur absolue 08-05-19
! du courant afin d'eviter d'attribuer une correction de flux trop forte
! A la couche kmin qui peut etre extrement petite
! A partir du 05-10-19: corriger de kmin A kmax un flux malgres tout calculE de 1 A kmax
       sum1=0.  ; sum2=0.
       do k=kmin_u(i,j),kmax
        sum1=sum1+(abs(velobc_u(i,j,k,2))+small1)
       enddo
       do k=1,kmax
        sum2=sum2+velobc_u(i,j,k,2)
       enddo
       x1=(xy_u(i,j,1)-sum2)/sum1
       do k=kmin_u(i,j),kmax
        velobc_u(i,j,k,2)=velobc_u(i,j,k,2)+x1*(abs(velobc_u(i,j,k,2))+small1)
       enddo
        velobc_u(i,j,:,:)=0.
      enddo ; enddo

      do j=1,jmax+1 ; do i=0,imax+1
! Pour eviter pb de precision on force le flux des couches ecrasees A zero
!      do k=1,kmin_v(i,j)-1
!       velobc_v(i,j,k,2)=0.
!      enddo
! On ajuste le transport sur les couches restantes
! Ne connaissant pas dz_u a priori, on pondere par la valeur absolue
! du courant afin d'eviter d'attribuer une correction de flux trop forte
! A la couche kmin qui peut etre extrement petite
       sum1=0.  ; sum2=0.
       do k=kmin_v(i,j),kmax
        sum1=sum1+(abs(velobc_v(i,j,k,2))+small1)
       enddo
       do k=1,kmax
        sum2=sum2+velobc_v(i,j,k,2)
       enddo
       x1=(xy_v(i,j,1)-sum2)/sum1
       do k=kmin_v(i,j),kmax
        velobc_v(i,j,k,2)=velobc_v(i,j,k,2)+x1*(abs(velobc_v(i,j,k,2))+small1)
       enddo
        velobc_v(i,j,:,:)=0.
      enddo ; enddo

#ifdef bidon
! Pour test seulement (bidonner sinon): verifier que la divergence des flux volumiques equilibre la variation de ssh_inst !08-05-19
      do j=1,jmax ; do i=1,imax

      x1=0.
      if(allocated(w0mofl_w))x1=-w0mofl_w(i,j,0)
      if(allocated(w_keq1_ofl_w))x1=x1+w_keq1_ofl_w(i,j,0)

       sum1=0. ; sum2=0. ; sum3=0. ; sum4=0.
       do k=1,kmax
        sum1=sum1+velobc_u(i  ,j  ,k,2)
        sum2=sum2+velobc_u(i+1,j  ,k,2)
        sum3=sum3+velobc_v(i  ,j  ,k,2)
        sum4=sum4+velobc_v(i  ,j+1,k,2)
       enddo
       

!     if(i+par%timax(1)==14.and.j+par%tjmax(1)==38) &
      write(10+par%rank,*)'flux div, dssh/dt',          &
      real(x1-( xy_u(i+1,j,1)-xy_u(i,j,1)               &
               +xy_v(i,j+1,1)-xy_v(i,j,1))/dxdy_t(i,j)) &
     ,real(x1-(      sum2-sum1+sum4-sum3 )/dxdy_t(i,j)) &
        ,real( (sshobc_w(i,j,2)-sshobc_w(i,j,1))/(ofl_readtime_next-ofl_readtime_prev) )
      enddo ; enddo
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)
#endif
      stop 'coco90'
#endif
      end subroutine offline_read_file_vel

!------------------------------------------------------------------------

      subroutine offline_read_file_botstress
      implicit none

! Get surface vertical velocity
       if(.not.allocated(maxbotstress_bef_w)) then !>>>
        allocate(maxbotstress_bef_w(imax,jmax)) ; maxbotstress_bef_w=0.
       endif                                       !>>>
       if(.not.allocated(maxbotstress_aft_w)) then !>>>
        allocate(maxbotstress_aft_w(imax,jmax)) ; maxbotstress_aft_w=0.
       endif                                       !>>>

      if(par%rank==0)write(6,*)'Get maxbotstress_bef_w with nc=',nc
! NOTE: je n'ai pas cree de liste particuliere pour maxbotstress je lis dans
! celle pour w....
      open(unit=3,file=trim(tmpdirname)//'liste_offline_w.binrec'     &
                 ,access='direct',recl=recordlength_,form='unformatted')
      read(3,rec=nc)filename_,ogcmtimecounter_,ofl_readtime_next,ofl_period_next
      close(3)
      if(par%rank==0) then  !.....
       write(6,'(a,a)')'offline fichier: ',trim(filename_)
       write(6,*)'ogcmtimecounter_=',ogcmtimecounter_
       write(6,*)'ofl_readtime_next in days=',ofl_readtime_next/86400.
       write(6,*)'ofl_period_next in hours=',ofl_period_next/3600.
      endif                 !.....

      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
      if(status/=0)stop ' stop offline_inout.f erreur 1'

      start(3)=ogcmtimecounter_  ; edge(3)=1        ! time
      start(2)=1+par%tjmax(1)    ; edge(2)=jmax+2   ! j
      start(1)=1+par%timax(1)    ; edge(1)=imax+2   ! i

                   status=nfmpi_inq_varid(ncid1,'maxbotstress_w',var_id)
      if(status/=0)stop ' stop offline_inout.f var_id maxbotstress'

      status=nfmpi_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)
      if(status/=0)stop ' Err 4454 nfmpi_inq_var'
      if(var_type==nf_real)  then !rrrrrrrrr
      status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:3),edge(1:3) &
                                          ,anyvar2d(0:imax+1,0:jmax+1))
      else                        !rrrrrrrrr
       stop 'Err maxbotstress type/=nf_real'
      endif                       !rrrrrrrrr

      status=nfmpi_close(ncid1)

      do j=1,jmax ; do i=1,imax
       maxbotstress_bef_w(i,j)=maxbotstress_aft_w(i,j)
       maxbotstress_aft_w(i,j)=          anyvar2d(i,j)
      enddo       ; enddo

      end subroutine offline_read_file_botstress


! ---------------------------------------------------------------------------------------------------------


!      subroutine offline_read_file_tenfonw
!      implicit none
!      integer :: logscale_=0
!
!! Get vertical mixing coef:
!      if(par%rank==0)write(6,*)'Get Kz with nc=',nc
!      open(unit=3,file=trim(tmpdirname)//'liste_offline_kz.binrec'           &
!                 ,access='direct',recl=recordlength_,form='unformatted')
!      read(3,rec=nc)filename_,ogcmtimecounter_
!      close(3)
!      if(par%rank==0)write(6,'(a,a)')'offline fichier kz: ',trim(filename_)
!      if(par%rank==0)write(6,*)'ogcmtimecounter_=',ogcmtimecounter_
!
!       start(4)=ogcmtimecounter_ ; edge(4)=1        ! time
!       start(3)=1                ; edge(3)=kmax+1   ! k
!       start(2)=1+par%tjmax(1)   ; edge(2)=jmax+2   ! j
!       start(1)=1+par%timax(1)   ; edge(1)=imax+2   ! i
!
!! CAS NEMO: autant de point t que de point Kz:
!      if(vert_axis_conv_end/=vert_axis_conv_start)edge(3)=kmax
!
!      if(vert_axis_conv_direc=='up') then  !05-11-12
!       kbegin_=1 ; kend_=kmax+1 ; kstep_=1
!      endif
!      if(vert_axis_conv_direc=='dw') then !05-11-12
!       kbegin_=kmax+1 ; kstep_=-1 ; kend_=kbegin_-edge(3)+1 ! kend_=2 si cas NEMO
!      endif
!
!      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
!      if(status/=0)stop ' stop offline_inout.f erreur 5'
!                   status=nfmpi_inq_varid(ncid1,'difv',var_id)
!      if(status/=0)status=nfmpi_inq_varid(ncid1,'votkeavt',var_id)
!      if(status/=0) then !logloglog?
!                        status=nfmpi_inq_varid(ncid1,'log10difv',var_id)
!           if(status/=0)status=nfmpi_inq_varid(ncid1,'vokzln10',var_id)
!                    if(status==0)logscale_=1
!      endif              !logloglog?
!      if(status/=0)stop ' stop offline_inout.f erreur var_id difv'
!
!       status=nfmpi_get_att_real(ncid1,var_id,'scale_factor',scalarval(1))
!       if(status/=0)scalarval(1)=1.
!       var_scalefactor=scalarval(1)
!                    status=nfmpi_get_att_real(ncid1,var_id,'add_offset',scalarval(1))
!       if(status/=0)status=nfmpi_get_att_real(ncid1,var_id,'offset',scalarval(1))
!       if(status/=0)scalarval(1)=0.
!       var_addoffset=scalarval(1)
!       status=nfmpi_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)
!
!      if(var_type==nf_real) then  !rrrrrrrrrr
!                   status=nfmpi_get_att_real(ncid1,var_id,'_FillValue',scalarval(1))
!      if(status/=0)status=nfmpi_get_att_real(ncid1,var_id,'_Fillvalue',scalarval(1)) !20-03-15
!      if(status/=0)stop ' module_offline err 3969 _FillValue not found'
!       filval=scalarval(1)
!       status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:4),edge(1:4) &
!                      ,anyvar3d(0:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
!       if(logscale_==1) then !llllllll> 13-03-15
!        do k=kbegin_,kend_ ,kstep_ ; do j=0,jmax+1 ; do i=0,imax+1
!         if(anyvar3d(i,j,k)==filval) then !--->
!            anyvar3d(i,j,k)=1.e-6
!         else                             !--->
!            anyvar3d(i,j,k)=10**anyvar3d(i,j,k)
!         endif                            !--->
!        enddo         ; enddo         ; enddo
!       else                  !llllllll>
!        do k=kbegin_,kend_ ,kstep_ ; do j=0,jmax+1 ; do i=0,imax+1
!         if(anyvar3d(i,j,k)==filval)anyvar3d(i,j,k)=1.e-6
!        enddo         ; enddo         ; enddo
!       endif                 !llllllll>
!      endif                       !rrrrrrrrrr
!      if(var_type==nf_short) then !ssssssssss
!                   status=nfmpi_get_att_int(ncid1,var_id,'_FillValue',i1val(1))
!      if(status/=0)stop ' module_offline err 3980 _FillValue not found'
!      filvalshort=i1val(1)
!        status=nfmpi_get_vara_int_all(ncid1,var_id,start(1:4),edge(1:4) &
!                      ,oflshort(0:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
!       do k=kbegin_,kend_ ,kstep_ ; do j=0,jmax+1 ; do i=0,imax+1
!        anyvar3d(i,j,k)=oflshort(i,j,k)*var_scalefactor+var_addoffset
!                     if(oflshort(i,j,k)==filvalshort)anyvar3d(i,j,k)=0.
!       enddo         ; enddo         ; enddo
!       if(logscale_==1) then !llllllll> 13-03-15
!        do k=kbegin_,kend_ ,kstep_ ; do j=0,jmax+1 ; do i=0,imax+1
!            anyvar3d(i,j,k)=10**anyvar3d(i,j,k)
!         if(oflshort(i,j,k)==filvalshort)anyvar3d(i,j,k)=0.
!        enddo         ; enddo         ; enddo
!       endif                 !llllllll>
!      endif                       !ssssssssss
!
!      if(status/=0)stop ' stop offline_inout.f erreur get difv'
!      status=nfmpi_close(ncid1)
!
!      do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
!       dfvofl_w(i,j,k,-1)=dfvofl_w(i,j,k,0)
!       dfvofl_w(i,j,k,0) =anyvar3d(i,j,k)
!      enddo ; enddo ; enddo
!
!! CAS NEMO: autant de point t que de point Kz:
!      if(flag_nemoffline==1) then !111>     ! cas nemo offline
!        if(vert_axis_conv_end/=vert_axis_conv_start) then !222>
!                       if(vert_axis_conv_direc=='dw') then  !333>
!                                     dfvofl_w(:,:,1,0)=dfvofl_w(:,:,2,0)
!                                                      endif !333>
!                                                    endif  !222>
!      endif                       !111>
!
!
!      if(ofl_tke==1) then !tke-tke-tke-> !05-03-15
!! Get tke: !05-03-15
!      if(par%rank==0)write(6,*)'Get tke with nc=',nc
!      open(unit=3,file=trim(tmpdirname)//'liste_offline_tke.binrec'           &
!                 ,access='direct',recl=recordlength_,form='unformatted')
!      read(3,rec=nc)filename_,ogcmtimecounter_
!      close(3)
!      if(par%rank==0)write(6,'(a,a)')'offline fichier tke: ',trim(filename_)
!      if(par%rank==0)write(6,*)'ogcmtimecounter_=',ogcmtimecounter_
!
!       start(4)=ogcmtimecounter_ ; edge(4)=1        ! time
!       start(3)=1                ; edge(3)=kmax+1   ! k
!       start(2)=1+par%tjmax(1)   ; edge(2)=jmax+2   ! j
!       start(1)=1+par%timax(1)   ; edge(1)=imax+2   ! i
!
!! CAS NEMO: autant de point t que de point Kz:
!      if(vert_axis_conv_end/=vert_axis_conv_start)edge(3)=kmax
!
!      if(vert_axis_conv_direc=='up') then  !05-11-12
!       kbegin_=1 ; kend_=kmax+1 ; kstep_=1
!      endif
!      if(vert_axis_conv_direc=='dw') then !05-11-12
!       kbegin_=kmax+1 ; kstep_=-1 ; kend_=kbegin_-edge(3)+1 ! kend_=2 si cas NEMO
!      endif
!
!      status=nfmpi_open(par%comm2d,filename_,nf_nowrite, MPI_INFO_NULL,ncid1)
!      if(status/=0)stop ' stop offline_inout.f erreur 5'
!                   status=nfmpi_inq_varid(ncid1,'tke',var_id)
!      if(status/=0) then !logloglog?
!                        status=nfmpi_inq_varid(ncid1,'log10tke',var_id)
!                    if(status==0)logscale_=1
!      endif              !logloglog?
!      if(status/=0)stop ' stop offline_inout.f erreur var_id tke'
!
!       status=nfmpi_get_att_real(ncid1,var_id,'scale_factor',scalarval(1))
!       if(status/=0)scalarval(1)=1.
!       var_scalefactor=scalarval(1)
!                    status=nfmpi_get_att_real(ncid1,var_id,'add_offset',scalarval(1))
!       if(status/=0)status=nfmpi_get_att_real(ncid1,var_id,'offset',scalarval(1))
!       if(status/=0)scalarval(1)=0.
!       var_addoffset=scalarval(1)
!       status=nfmpi_inq_var(ncid1,var_id,texte30,var_type,var_dims,tabdim,i4)
!
!      if(var_type==nf_real) then  !rrrrrrrrrr
!       status=nfmpi_get_vara_real_all(ncid1,var_id,start(1:4),edge(1:4) &
!                      ,anyvar3d(0:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
!      endif                       !rrrrrrrrrr
!      if(var_type==nf_short) then !ssssssssss
!        status=nfmpi_get_vara_int_all(ncid1,var_id,start(1:4),edge(1:4) &
!                      ,oflshort(0:imax+1,0:jmax+1,kbegin_:kend_:kstep_))
!       do k=kbegin_,kend_ ,kstep_ ; do j=0,jmax+1 ; do i=0,imax+1
!        anyvar3d(i,j,k)=oflshort(i,j,k)*var_scalefactor+var_addoffset
!        anyvar3d(i,j,k)=10**anyvar3d(i,j,k)
!       enddo         ; enddo         ; enddo
!       if(logscale_==0)stop 'offline Kz must be en log10 scale'
!      endif                       !ssssssssss
!
!      if(status/=0)stop ' stop offline_inout.f erreur get difv'
!      status=nfmpi_close(ncid1)
!
!      do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
!       tkeofl_w(i,j,k,-1)=tkeofl_w(i,j,k,0)
!       tkeofl_w(i,j,k,0) =anyvar3d(i,j,k)
!      enddo ; enddo ; enddo
!
!! CAS NEMO: autant de point t que de point Kz:
!      if(flag_nemoffline==1) then !111>   ! cas nemo offline
!        if(vert_axis_conv_end/=vert_axis_conv_start) then !222>
!                       if(vert_axis_conv_direc=='dw') then  !333>
!                                     tkeofl_w(:,:,1,0)=tkeofl_w(:,:,2,0)
!                                                      endif !333>
!                                                    endif !222>
!      endif                       !111>
!      endif               !tke-tke-tke-> !05-03-15
!
!      end subroutine offline_read_file_tenfonw




















































!.........................................................................................

      subroutine offline_write_oasisgrids !26-06-15
      implicit none
      integer ncid_
      integer(kind=MPI_OFFSET_KIND) :: iglbp2,jglbp2
      integer(kind=MPI_OFFSET_KIND) start(2)
      integer(kind=MPI_OFFSET_KIND) edge(2)
#ifdef synopsis
       subroutinetitle='offline_write_oasisgrids'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      iglbp2=iglb+2
      jglbp2=jglb+2

! Fichier lon,lat

      status=nfmpi_create(par%comm2d,trim(tmpdirname)//'oasis_grid.nc' &
      ,nf_clobber + NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)

      texte30='x_'//cl_grid_name 
      status=nfmpi_def_dim(ncid_,trim(texte30),iglbp2,i_t_dim)

      texte30='y_'//cl_grid_name 
      status=nfmpi_def_dim(ncid_,trim(texte30),jglbp2,j_t_dim)

      vardim(1)=i_t_dim 
      vardim(2)=j_t_dim 

      texte30=cl_grid_name//'.lat'
      status=nfmpi_def_var(ncid_,trim(texte30),nf_double,2,vardim(1:2),var_id)
      texte30=cl_grid_name//'.lon'
      status=nfmpi_def_var(ncid_,trim(texte30),nf_double,2,vardim(1:2),var_id)

      status=nfmpi_enddef(ncid_)

      istr=1                     ; jstr=1
      iend=imax-2                ; jend=jmax-2
      if(par%timax(1)==0)istr=0               !04-02-13
      if(imax+par%timax(1)==iglb)iend=imax+1

      if(par%tjmax(1)==0)jstr=0
      if(jmax+par%tjmax(1)==jglb)jend=jmax+1

      start(1)=1+par%timax(1)+istr ; start(2)=1+par%tjmax(1)+jstr
       edge(1)=iend-istr+1         ;  edge(2)=jend-jstr+1

      texte30=cl_grid_name//'.lat'
      status=nfmpi_inq_varid(ncid_,trim(texte30),var_id)
      anyv3d(istr:iend,jstr:jend,1,1)=lat_t(istr:iend,jstr:jend)*rad2deg
      status=nfmpi_put_vara_double_all(ncid_,var_id              &
                               ,start(1:2),edge(1:2)             &
                               ,anyv3d(istr:iend,jstr:jend,1,1))

      texte30=cl_grid_name//'.lon' 
      status=nfmpi_inq_varid(ncid_,trim(texte30),var_id)
      anyv3d(istr:iend,jstr:jend,1,1)=lon_t(istr:iend,jstr:jend)*rad2deg!29-06-15
      status=nfmpi_put_vara_double_all(ncid_,var_id              &
                               ,start(1:2),edge(1:2)             &
                               ,anyv3d(istr:iend,jstr:jend,1,1))

      status=nfmpi_close(ncid_)


! Fichier mask

      status=nfmpi_create(par%comm2d,trim(tmpdirname)//'oasis_mask.nc' &
      ,nf_clobber + NF_64BIT_OFFSET, MPI_INFO_NULL,ncid_)

      texte30='x_'//cl_grid_name 
      status=nfmpi_def_dim(ncid_,trim(texte30),iglbp2,i_t_dim)

      texte30='y_'//cl_grid_name 
      status=nfmpi_def_dim(ncid_,trim(texte30),jglbp2,j_t_dim)

      vardim(1)=i_t_dim 
      vardim(2)=j_t_dim 

      texte30=cl_grid_name//'.msk'
      status=nfmpi_def_var(ncid_,trim(texte30),nf_int,2,vardim(1:2),var_id)

      status=nfmpi_enddef(ncid_)

      istr=1                     ; jstr=1
      iend=imax-2                ; jend=jmax-2
      if(par%timax(1)==0)istr=0               !04-02-13
      if(imax+par%timax(1)==iglb)iend=imax+1

      if(par%tjmax(1)==0)jstr=0
      if(jmax+par%tjmax(1)==jglb)jend=jmax+1

      start(1)=1+par%timax(1)+istr ; start(2)=1+par%tjmax(1)+jstr
       edge(1)=iend-istr+1         ;  edge(2)=jend-jstr+1

      texte30=cl_grid_name//'.msk'
      status=nfmpi_inq_varid(ncid_,trim(texte30),var_id)
      anyv3dint(istr:iend,jstr:jend,1)=1-mask_t(istr:iend,jstr:jend,kmax)
      status=nfmpi_put_vara_int_all(ncid_,var_id              &
                               ,start(1:2),edge(1:2)             &
                               ,anyv3dint(istr:iend,jstr:jend,1))


      status=nfmpi_close(ncid_)

      end subroutine offline_write_oasisgrids

!.........................................................................................

      subroutine offline_another_model_divergence_2 !01-07-16
      implicit none
      double precision implicit_fric_,explicit_fric_,reduced_grav_,dt_,rap2_
      integer loop_,loopmax_
#ifdef synopsis
       subroutinetitle='offline_another_model_divergence'
       subroutinedescription= &
       'Adjust depth-averaged currents in order to restore nemo ssh'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Calcul des perturbations de vitesse moyenne de la pseudo surface libre
! pour equilibrage de la divergence
! Details dans:
! https://docs.google.com/document/d/1gh3vxqZ5ma-YzOOunahy0eorT2KCz1og2Yun_OBxH0g/edit 

      if(cfl_reduce>1.91) then   !21-06-16
       write(6,*)'cfl_reduce=',cfl_reduce
       stop 'ERR offline_another_model_divergence_2 cfl_reduce>1.9'
      endif
      if(cfl_reduce<1.) then     !21-06-16
! Cette subroutine a une stabilite comparable A celle du scheme FB, donc
! utiliser si possible un grand pas de temps pour economiser de la cpu
       write(6,*)'cfl_reduce=',cfl_reduce
       stop 'WARNING offline_another_model_divergence_2 cfl_reduce < 1'
      endif

! Pour eviter l'emergence du mode numerique on preserve la continuite
! temporelle du forcage ssh_int. Le pourcentage timeweightobc(1) vient d'etre calculE
! avant l'appel A cette routine pour le temps elapsedtime_now
! On refait ce calcul pour le temps elapsedtime_now+dti_fw et on obtient
! le pourcentage rap2_ de l'iteration suivante. On interpole lineairement
! entre les 2
       if(ofl_period_prev<ofl_period_next) then !>>>>>>
        x1=ofl_readtime_prev+ofl_period_prev
        rap2_=1.-max( (x1-elapsedtime_now-dti_fw  )                   &
                     /(x1-ofl_readtime_prev) , 0.d00 )
       else                                     !>>>>>>
        x1=ofl_period_prev+2.*ofl_readtime_prev-ofl_readtime_next
        rap2_=1.-min( (ofl_readtime_next-elapsedtime_now-dti_fw)      &
                     /(ofl_readtime_next-x1             ),1.d00)
       endif                                    !>>>>>>
  

! Parametre de friction:
      x1=dte_lp*2.6d-3*1.d-3    ! dt*cd*moduleU
!     x1=0.
      implicit_fric_=1./(1.+x1) ! 1/(1+dt*cd*moduleU)
      explicit_fric_=1.-x1      ! 1-dt*cd*moduleU


! Tableaux pour faire la moyenne
      xy_u(:,:,0)=0.
      xy_v(:,:,0)=0.

     

      loopmax_=2*nint( dti_lp / dte_lp)
!     loopmax_=10
      dt_=dti_lp/loopmax_
!     reduced_grav_=0.1*grav*(dte_fw/dt_)**2
      reduced_grav_=    grav*(dte_fw/dt_)**2

      

! Les loopmax_ iterations avec un pas de temps unitiare de dte_lp representent un
! temps egal A dti_lp qui est la periode d'integration voulue sur cette
! sequence....

!     sshobc_w=0.
!     velbar_u=0.
!     velbar_v=0.
!     do i=1,imax
!     do j=1,jmax
!      ssh_int_w(i,j,2)=0.001*cos(pi*real(i))*cos(pi*real(j))
!     enddo
!     enddo
!     call graph_out
!     stop 'gigi'

      do iteration2d=1,loopmax_


       rap=rap2_*real(iteration2d)/real(loopmax_)    &
            +(1.-real(iteration2d)/real(loopmax_))*timeweightobc(1)

! interpolation sshobc si cas nemo ou cas symphonie "anciens" fichers
       do j=0,jmax+1 ; do i=0,imax+1 ! j=1,jmax !25-09-14 extension des boucles pour confort graphique
         sshobc_w(i,j,1)=(1.-rap)*sshobc_w(i,j,0) &
                            +rap *sshobc_w(i,j,2)

!       if(i+par%timax(1)==280.and.j+par%tjmax(1)==90) &
!       write(80,*)(elapsedtime_now+dti_fw*real(iteration2d)/real(loopmax_))/86400.   &
!                  ,real(rap),sshobc_w(i,j,1)

       enddo ; enddo

       



      do j=2,jmax-1
      do i=2,imax
       velbar_u(i,j,2)=(                                             &
       velbar_u(i,j,1)                                               &
                -dt_*reduced_grav_*( ssh_int_w(i  ,j  ,2)            &
                                    -ssh_int_w(i-1,j  ,2)            &
                                     -sshobc_w(i  ,j  ,1)            &
                                     +sshobc_w(i-1,j  ,1))/dx_u(i,j) &


        +0.02*(velbar_u(i+1,j  ,1)  &
              +velbar_u(i-1,j  ,1)  &
              +velbar_u(i  ,j+1,1)  &
              +velbar_u(i  ,j-1,1)  &
           -4.*velbar_u(i  ,j  ,1)) &

       )*mask_u(i,j,kmax+1)*implicit_fric_
      enddo
      enddo

      do j=2,jmax
      do i=2,imax-1
       velbar_v(i,j,2)=(                                             &
       velbar_v(i,j,1)                                               &
                -dt_*reduced_grav_*( ssh_int_w(i  ,j  ,2)            &
                                    -ssh_int_w(i  ,j-1,2)            &
                                     -sshobc_w(i  ,j  ,1)            &
                                     +sshobc_w(i  ,j-1,1))/dy_v(i,j) &

        +0.02*(velbar_v(i+1,j  ,1)  &
              +velbar_v(i-1,j  ,1)  &
              +velbar_v(i  ,j+1,1)  &
              +velbar_v(i  ,j-1,1)  &
           -4.*velbar_v(i  ,j  ,1)) &

       )*mask_v(i,j,kmax+1)*implicit_fric_
      enddo
      enddo

      call obc_ext_velbar_mpi(2) ! echange x sur velbar_u(:,:,2) et y sur velbar_v(:,:,2) !15-07-14

      do j=1,jmax
      do i=1,imax

       ssh_int_w(i,j,2)=                                      &
       ssh_int_w(i,j,2)                                       &
           -dt_                                            &
                  *( & !ooo>
           dy_u(i+1,j)*h_u(i+1,j)*velbar_u(i+1,j,2)    &
          -dy_u(i  ,j)*h_u(i  ,j)*velbar_u(i  ,j,2)    &
          +dx_v(i,j+1)*h_v(i,j+1)*velbar_v(i,j+1,2)    &
          -dx_v(i,j  )*h_v(i,j  )*velbar_v(i,j  ,2)    &
                   ) & !ooo>
                   /dxdy_t(i,j)

      enddo
      enddo

        do j=1,jmax ; do i=1,imax+1
         xy_u(i,j,0)=xy_u(i,j,0)+velbar_u(i,j,2)
        enddo ; enddo
        do j=1,jmax+1 ; do i=1,imax
         xy_v(i,j,0)=xy_v(i,j,0)+velbar_v(i,j,2)
        enddo ; enddo

!     call graph_out

! Update des vitesses
      do j=1,jmax ; do i=1,imax+1
       velbar_u(i,j,1)=velbar_u(i,j,2)
      enddo ; enddo
      do j=1,jmax+1 ; do i=1,imax
       velbar_v(i,j,1)=velbar_v(i,j,2)
      enddo ; enddo

      enddo ! fin de boucle sur loop_


! Apres la boucle iterative, transformer les cumuls des courants en des
! transports moyens:
        x0=1./real(loopmax_)
        do j=1,jmax ; do i=1,imax+1
         xy_u(i,j,0)=xy_u(i,j,0)*dy_u(i,j)*h_u(i,j)*x0
        enddo ; enddo
        do j=1,jmax+1 ; do i=1,imax
         xy_v(i,j,0)=xy_v(i,j,0)*dx_v(i,j)*h_v(i,j)*x0
        enddo ; enddo


! Deplacer apres l'appel de la routine
      do j=1,jmax
      do i=1,imax
       hz_w(i,j,2)=h_w(i,j)+ssh_int_w(i,j,2)
      enddo
      enddo

! AJUSTER LE COURANT U
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,1)=small1
      enddo
      enddo
      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,1)=xy_u(i,j,1)+mask_u(i,j,k)*abs(veldydz_u(i,j,k,1))
      enddo
      enddo
      enddo
      do j=1,jmax
      do i=1,imax+1
       xy_u(i,j,2)=xy_u(i,j,0)     &
                  /xy_u(i,j,1)
      enddo
      enddo
      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1
            veldydz_u(i,j,k,1)=     &
            veldydz_u(i,j,k,1)+     &
               xy_u(i,j,2)          &
              *mask_u(i,j,k)        &
       *abs(veldydz_u(i,j,k,1))
      enddo
      enddo
      enddo

! AJUSTER LE COURANT V
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,1)=small1
      enddo
      enddo
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,1)=xy_v(i,j,1)+mask_v(i,j,k)*abs(veldxdz_v(i,j,k,1))
      enddo
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
       xy_v(i,j,2)=xy_v(i,j,0)     &
                  /xy_v(i,j,1)
      enddo
      enddo
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax
            veldxdz_v(i,j,k,1)=     &
            veldxdz_v(i,j,k,1)+     &
                xy_v(i,j,2)         &
              *mask_v(i,j,k)        &
       *abs(veldxdz_v(i,j,k,1))
      enddo
      enddo
      enddo

#ifdef bidon
! Verification:
      if(iteration3d==2) then
      do j=1,jmax+1
      do i=1,imax+1
       xy_u(i,j,1)=0.
       xy_v(i,j,1)=0.
       do k=1,kmax
        xy_u(i,j,1)=xy_u(i,j,1)+veldydz_u(i,j,k,1)
        xy_v(i,j,1)=xy_v(i,j,1)+veldxdz_v(i,j,k,1)
       enddo
      enddo
      enddo
      do j=1,jmax
      do i=1,imax
      if(mask_t(i,j,kmaxp1)==1) then
       write(10+par%rank,*)   &
       (ssh_int_w(i,j,2)-ssh_int_w(i,j,0))/dti_lp &
      +(xy_u(i+1,j,1)-xy_u(i,j,1)+xy_v(i,j+1,1)-xy_v(i,j,1))/dxdy_t(i,j) &
       ,-omega_w(i,j,kmaxp1,1)
      endif
      enddo
      enddo
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)
#endif
      stop 'bibi2'
      endif
#endif

!     call graph_out
!     stop 'bibi4'
!     if(iteration3d==2)stop 'bibi3'

      end subroutine offline_another_model_divergence_2

!.........................................................................................

      subroutine offline_oldlist(loop_) !09-02-17
      implicit none
      integer loop_,flag_,unit_

! Si la liste binaire est presente dans le tmp alors se positionner
! dans celle-ci pour demarrer la simulation:

        if(par%rank==0)write(6,'(a,a)')'In offline_oldlist open ',trim(texte80(loop_))

        flag_=0
        unit_=s_unit(7)
        open(unit=unit_,file=trim(texte80(loop_)) &
                       ,access='direct'     &
                       ,recl=recordlength_  &
                       ,form='unformatted'  &
                       ,status='old'        &
                       ,iostat=k0)
         if(k0/=0)stop 'Err 2424'
         nc=0

         do while (flag_==0)

          nc=nc+1

!              ofl_rec_now=     ofl_rec_next
          ofl_readtime_prev=ofl_readtime_next 
            ofl_period_now=  ofl_period_next


!        ofl_rec_next=nc
         if(nc/=1) then !pmxpmx>
          read(unit_,rec=nc,iostat=i0)filename_,i1          &
                           ,ofl_readtime_next &
                             ,ofl_period_next
         else           !pmxpmx>
! Verifier la date repere de nc=1:
          read(unit_,rec=nc,iostat=i0)filename_,i1          &
                           ,ofl_readtime_next &
                             ,ofl_period_next &
                             ,i1,i2,i3,i4,i5,i6

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
             stop 'Err 5252'
          endif                     !bugbug>

         endif          !pmxpmx>

         if(ofl_readtime_next<=elapsedtime_now) then
!            write(6,*)'bobi ofl_readtime_next,elapsedtime_now,nc',ofl_readtime_next,elapsedtime_now,nc
             ofl_rec_now=nc
         endif

         if(i0/=0)stop 'Err 5228'

         ofl_rec_max=nc

! On cherche les 2 records dont les temps encadrent le temps initial elapsedtime_now
         if(nc>1) then !pmx>
          if(     elapsedtime_now>=ofl_readtime_prev &
             .and.elapsedtime_now< ofl_readtime_next)flag_=1
         endif         !pmx>

         enddo

! Poursuivre la lecture jusqu'au bout du fichier pour determiner ofl_rec_max
         do while (i0==0)
          nc=nc+1
          read(unit_,rec=nc,iostat=i0)filename_,i1,x1,x2
          ofl_rec_max=nc-1
!         write(6,*)'nc i0 ofl_rec_max',nc,i0,ofl_rec_max
         enddo

         close(unit_)

        if(flag_==0) STOP 'Pas de fichiers OGCM autour du temps initial'

        if(par%rank==0) then !ooo>
         write(6,*)'k1=                           ',k1
         write(6,*)'elapsedtime_now=              ',elapsedtime_now
         write(6,*)'ofl_readtime_prev',ofl_readtime_prev
         write(6,*)'ofl_readtime_next',ofl_readtime_next
         write(6,*)'ofl_rec_now     ',ofl_rec_now
!        write(6,*)'ofl_rec_next     ',ofl_rec_next
        endif                !ooo>

      end subroutine offline_oldlist

!......................................................................................................

      subroutine offline_bin(txt2) !08-06-18
      implicit none
      character*2 txt2
      integer(kind=2),dimension(:,:,:),allocatable :: &
       anyshort3d       
      if(.not.allocated(anyshort3d)) allocate (anyshort3d (-1:imax+2,-1:jmax+2,0:kmax+1)) ; anyshort3d=0

! Si fichier offline binaire n'existe pas (offline_status=0), alors le creer avec le nom de la version du modele en "header"
! Noter que offline_binary_status est dans le fichier restart et donc que les fichiers ne sont pas ecrasE si le modele repart
! du restart
       if(offline_binary_status==0) then !>>> !08-06-18
        write(offline_binrec_name,'(a,a,i0,a)')trim(directory_offline),'offline',par%rank,'.bin'
!'
!       offline_binrec_name=trim(directory_offline)//'offline_'//dom_c//'.binrec'
        open(unit=8,file=trim(offline_binrec_name),access='sequential',form='unformatted')
!       write(8)model_name ! character*60
        write(8)iglb,jglb,imax,jmax,kmax
        write(time_units,'(a14,i4,a1,i2,4(a1,i2))')                    & !units
         'seconds since ',datesim(1,1),'-',datesim(2,1),'-',datesim(3,1) &
        ,' ',datesim(4,1),':',datesim(5,1),':',datesim(6,1)
        if(time_units(20:20)==' ')time_units(20:20)='0'
        if(time_units(23:23)==' ')time_units(23:23)='0'
        if(time_units(26:26)==' ')time_units(26:26)='0'
        if(time_units(29:29)==' ')time_units(29:29)='0'
        if(time_units(32:32)==' ')time_units(32:32)='0'
        write(8)time_units
!       close(8)
        offline_binary_status=1 
       endif                           !>>>

      ksecu=0
!     open(unit=8,file=trim(offline_binrec_name),access='sequential',form='unformatted',position='append')
        write(8)date_shortname
        write(8)elapsedtime_now ! real*8        temps ecoule en secondes
        write(8)ofl_period_now
        write(8)texte80(1)      ! character*200 nom de la variable
        write(8)texte80(2)      ! character*200 unites
        write(8)texte80(3)      ! character*200 long name
        write(8)texte80(4)      ! character*200 standard name
        write(8)texte80(5)      ! character*200 TZYX ou TYX
        write(8)texte80(7)      ! character*200 real, short, double
        write(8)txt2            ! character*2   _t _u _v _w
        write(8)filval
        write(8)var_scalefactor
        write(8)var_addoffset
        write(8)par%timax(1)    ! integer       i local --> i global
        write(8)par%tjmax(1)    ! integer       j local --> j global
        if(texte80(5)=='TYX') then  !2d2d>
             write(8)lbound(anyvar2d) ! 2 integer
             write(8)ubound(anyvar2d) ! 2 intger
                    write(8)anyvar2d  ! real 2d
             ksecu=1
        endif                       !2d2d
        if(texte80(5)=='TZYX') then !3d3d>
             write(8)lbound(anyvar3d) ! 3 integer
             write(8)ubound(anyvar3d) ! 3 intger
      if(trim(texte80(7))=='short') then
                    anyshort3d=nint(anyvar3d)
                    write(8)anyshort3d  !integer*2
          else
                    write(8)anyvar3d  ! real 3d
          endif
          
             ksecu=1
        endif                       !3d3d
!     close(8)
      if(ksecu==0) then ! debug >
       write(6,'(a,a)')'texte80(5)=',trim(texte80(5))
       stop ' unrecognised texte80(5)'
      endif             ! debug >

      end subroutine offline_bin

!......................................................................................................

      subroutine offlinegbin(txt2) !18-06-18
      implicit none
      character*2 txt2

! Fichier de grille binaire

        if(offline_gridbin_status==0) then !reset>
         open(unit=3 &
             ,file=trim(directory_offline)//'grid_'//dom_c//'.bin',access='sequential',form='unformatted')
             write(3)iglb,jglb,imax,jmax,kmax
             close(3)
           offline_gridbin_status=1
        endif                              !reset>

        open(unit=3 &
           ,file=trim(directory_offline)//'grid_'//dom_c//'.bin',access='sequential',form='unformatted',position='append')
        ksecu=0
        write(3)texte80(1)      ! character*200 nom de la variable
        write(3)texte80(2)      ! character*200 unites
        write(3)texte80(3)      ! character*200 long name
        write(3)texte80(4)      ! character*200 standard name
        write(3)texte80(5)      ! character*200 TZYX ou TYX
        write(3)txt2            ! character*2   _t _u _v _w
        write(3)par%timax(1)    ! integer       i local --> i global
        write(3)par%tjmax(1)    ! integer       j local --> j global
        write(3)elapsedtime_now ! real*8        temps ecoule en secondes
        write(3)year_now,month_now,day_now,hour_now,minute_now,second_now ! 6 integer pour la date
        if(texte80(5)=='TYX'.or.    &
           texte80(5)=='YX' ) then  !2d2d>
             write(3)lbound(anyvar2d) ! 2 integer
             write(3)ubound(anyvar2d) ! 2 intger
                    write(3)anyvar2d  ! real 2d
             ksecu=1
        endif                       !2d2d
        if(texte80(5)=='TZYX'.or.   &
           texte80(5)=='ZYX' ) then !3d3d>
             write(3)lbound(anyvar3d) ! 3 integer
             write(3)ubound(anyvar3d) ! 3 intger
                    write(3)anyvar3d  ! real 3d
             ksecu=1
        endif                       !3d3d
        if(ksecu==0) then ! debug >
         write(6,'(a,a)')'texte80(5)=',trim(texte80(5))
         stop ' unrecognised texte80(5)'
        endif             ! debug >
      close(3)

      end subroutine offlinegbin

!...............................................................................

      subroutine offline_bin2netcdf !15-08-22
      implicit none
      integer :: loop_ &
                ,iglbbf_,jglbbf_,kmaxbf_  & ! iglb,jglb,imax,jmax,kmax de la simulation qui a produit les fichiers
                ,error_,flag_reset_=0
      integer,dimension(:),allocatable :: imaxbf,jmaxbf
      character*2 txt2
!     character*4 an_texte

!............................................
! combien de fichiers binaires dans la liste?
      binaryfilelist_len=0
      open(unit=3,file='binary2netcdf_list')
 7124 read(3,*,end=7125)texte90
      binaryfilelist_len=binaryfilelist_len+1
      goto 7124
 7125 close(3)

      
!............................................
! le tableau binaryfile_inout dit si un fichier binaire emarge (=1)
! ou n'emarge pas (=0) sur le rank.
      allocate(binaryfile_inout(binaryfilelist_len)) ; binaryfile_inout=0
      allocate(binaryfile_name(binaryfilelist_len)) ; binaryfile_name='none' 
!'
      allocate(imaxbf(binaryfilelist_len)) ; imaxbf=0
      allocate(jmaxbf(binaryfilelist_len)) ; jmaxbf=0

      open(unit=3,file='binary2netcdf_list')
       do loop_=1,binaryfilelist_len
        read(3,'(a)')binaryfile_name(loop_)
       enddo
      close(3)

      do loop_=1,binaryfilelist_len
        open(unit=3,file=trim(binaryfile_name(loop_)) &
                   ,access='sequential'               &
                   ,form='unformatted')
         read(3,iostat=error_)iglbbf_,jglbbf_,imaxbf(loop_),jmaxbf(loop_),kmaxbf_  ! iglb,jglb,imax,jmax,kmax de la simulation qui a produit les fichiers
         if(error_/=0)stop 'STOP ERR binary2netcdf'
         read(3)time_units

         read(3)date_shortname
         read(3)elapsedtime_now ! real*8        temps ecoule en secondes
         read(3)ofl_period_now
         read(3)texte80(1)      ! character*200 nom de la variable
         read(3)texte80(2)      ! character*200 unites
         read(3)texte80(3)      ! character*200 long name
         read(3)texte80(4)      ! character*200 standard name
         read(3)texte80(5)      ! character*200 TZYX ou TYX
         read(3)texte80(7)      ! character*200 real, short, double
         read(3)txt2            ! character*2   _t _u _v _w
         read(3)filval
         read(3)var_scalefactor
         read(3)var_addoffset
         read(3)i0 ! par%timax(1) de la simulation qui a produit les fichiers
         read(3)j0 ! par%tjmax(1) de la simulation qui a produit les fichiers

! juste pour faire du tri
!        an_texte=date_shortname(1:4)

         if(texte80(5)=='TYX') then  !2d2d>
          read(3)i1,j1 ! lbound(anyvar2d) ! 2 integer
          read(3)i2,j2 ! ubound(anyvar2d) ! 2 intger
         endif                       !2d2d
         if(texte80(5)=='TZYX') then !3d3d>
          read(3)i1,j1,k1 ! lbound(anyvar3d) ! 3 integer
          read(3)i2,j2,k2 ! ubound(anyvar3d) ! 3 intger
         endif                       !3d3d

! Debut et fin des indices globaux des tableaux binaires 
         i1=i1+i0 ! debut d'indice i global des tableaux du fichier binaire
         i2=i2+i0 ! fin   d'indice i global des tableaux du fichier binaire
         j1=j1+j0 ! debut d'indice j global des tableaux du fichier binaire
         j2=j2+j0 ! fin   d'indice j global des tableaux du fichier binaire
! Les memes mais avec les indices locaux du convertisseur:
         i1=i1-par%timax(1)
         i2=i2-par%timax(1)
         j1=j1-par%tjmax(1)
         j2=j2-par%tjmax(1)
        close(3)
        binaryfile_inout(loop_)=0
        do i3=i1,i2,i2-i1
        do j3=j1,j2,j2-j1
         if(i3>=0.and.i3<=imax+1.and.j3>=0.and.j3<=jmax+1) then
           binaryfile_inout(loop_)=1 ! Ce fichier binaire emarge sur le rank du convertisseur
         endif
        enddo
        enddo
      enddo ! loop_


! Securite PNETCF: attention parce qu'avec pnetcdf tous les rank doivent
! faire les memes actions, il n'est pas possible qu'un rank soit sans aucun
! fichier binaire A traiter. Si A l'issue de la precedente etape aucun
! fichier n'est identifiE comme emargeant sur le rank (c.a.d. maxval(binaryfile_inout)=0),
! alors on force le premier fichier binaire A Etre traitE en tant que fichier emargeant
       if(maxval(binaryfile_inout)==0)binaryfile_inout(1)=1

!......................................................
! Ouvrir, une fois pour toutes, tous les fichiers binaires emargeant sur le rank
      do loop_=1,binaryfilelist_len
       if(binaryfile_inout(loop_)==1) then
        open(unit=10+loop_                            &
                   ,file=trim(binaryfile_name(loop_)) &
                   ,access='sequential'               &
                   ,form='unformatted')
        read(10+loop_,iostat=error_)iglbbf_,jglbbf_,imaxbf(loop_),jmaxbf(loop_),kmaxbf_  ! iglb,jglb,imax,jmax,kmax de la simulation qui a produit les fichiers
        if(error_/=0)stop 'STOP ERROR 1 binary2netcdf'
        read(10+loop_)time_units
       endif
      enddo ! loop_

!......................................................
! lire les fichiers concernEs:
 7217 continue
      flag_reset_=0 ! sera =1 apres lecture du premier fichier
      do loop_=1,binaryfilelist_len ! boucle sur les fichiers binaires
       if(binaryfile_inout(loop_)==1) then !pmx>

         read(10+loop_,iostat=error_)date_shortname
         if(error_/=0)goto 7219
         read(10+loop_)elapsedtime_now ! real*8        temps ecoule en secondes
         read(10+loop_)ofl_period_now

!....
! Detecter le moment de passer A un nouveau fichier netcdf:
         if(date_shortname(1:15)/=date_shortname_b(1:15))flag_ncfile=0
         if(flag_ncfile==0) then !ooo>
           call offline_bin2netcdf_step1 ! nf_create etc....
           flag_ncfile=1
         endif                   !ooo>
         date_shortname_b(1:15)=date_shortname(1:15)
!....

         read(10+loop_)texte80(1)      ! character*200 nom de la variable
         read(10+loop_)texte80(2)      ! character*200 unites
         read(10+loop_)texte80(3)      ! character*200 long name
         read(10+loop_)texte80(4)      ! character*200 standard name
         read(10+loop_)texte80(5)      ! character*200 TZYX ou TYX
         read(10+loop_)texte80(7)      ! character*200 real, short, double
         read(10+loop_)txt2            ! character*2   _t _u _v _w
         read(10+loop_)filval
         read(10+loop_)var_scalefactor
         read(10+loop_)var_addoffset
         read(10+loop_)i0 ! par%timax(1) de la simulation qui a produit les fichiers
         read(10+loop_)j0 ! par%tjmax(1) de la simulation qui a produit les fichiers

         if(texte80(5)=='TYX') then  !2d2d>
          if(flag_reset_==0)anyvar2d=filval ! reset des zones "proc inutiles"
          read(10+loop_)i1,j1 ! lbound(anyvar2d) ! 2 integer
          read(10+loop_)i2,j2 ! ubound(anyvar2d) ! 2 intger
          allocate(binaryfile_tab2d_r4(i1:i2,j1:j2))
          read(10+loop_,iostat=error_)binaryfile_tab2d_r4 ; if(error_/=0)stop 'STOP ERROR 2 binary2netcdf'
          i1=2 ; i2=imaxbf(loop_)-1 ; j1=2 ; j2=jmaxbf(loop_)-1
          do j3=j1,j2 ; do i3=i1,i2
           i=i3+i0-par%timax(1) ; j=j3+j0-par%tjmax(1)
           if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1)anyvar2d(i,j)=binaryfile_tab2d_r4(i3,j3)
          enddo      ; enddo
          flag_reset_=1
         endif                       !2d2d

         if(texte80(5)=='TZYX') then !3d3d>
          if(flag_reset_==0)anyvar3d=filval ! reset des zones "proc inutiles"
          read(10+loop_)i1,j1,k1 ! lbound(anyvar3d) ! 3 integer
          read(10+loop_)i2,j2,k2 ! ubound(anyvar3d) ! 3 intger
          if(trim(texte80(7))=='real')then
          allocate(binaryfile_tab3d_r4(i1:i2,j1:j2,k1:k2))
          read(10+loop_,iostat=error_)binaryfile_tab3d_r4 ; if(error_/=0)stop 'STOP ERROR 3 binary2netcdf'
          else if (trim(texte80(7))=='short')then
          allocate(binaryfile_tab3d_i2(i1:i2,j1:j2,k1:k2))
          read(10+loop_,iostat=error_)binaryfile_tab3d_i2 ; if(error_/=0)stop 'STOP ERROR 3 binary2netcdf'
          else
           stop 'erreur format sur texte80(7)'
          endif

          i1=2 ; i2=imaxbf(loop_)-1 ; j1=2 ; j2=jmaxbf(loop_)-1
          if(              i0==0      )i1=0
          if(imaxbf(loop_)+i0==iglbbf_)i2=imaxbf(loop_)+1
          if(              j0==0      )j1=0
          if(jmaxbf(loop_)+j0==jglbbf_)j2=jmaxbf(loop_)+1

          do k=k1,k2 ; do j3=j1,j2 ; do i3=i1,i2
           i=i3+i0-par%timax(1) ; j=j3+j0-par%tjmax(1)
          if(trim(texte80(7))=='real')then
           if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1)anyvar3d(i,j,k)=binaryfile_tab3d_r4(i3,j3,k)
          else if(trim(texte80(7))=='short')then
           if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1)anyvar3d(i,j,k)=binaryfile_tab3d_i2(i3,j3,k)
          endif
          enddo ; enddo      ; enddo
          flag_reset_=1
         endif                       !3d3d

         if(texte80(5)=='TYX') deallocate(binaryfile_tab2d_r4)
         if(texte80(5)=='TZYX'.and.trim(texte80(7))=='real')deallocate(binaryfile_tab3d_r4)
         if(texte80(5)=='TZYX'.and.trim(texte80(7))=='short')deallocate(binaryfile_tab3d_i2)

       endif                               !pmx>
      enddo ! loop_ boucle sur les fichiers binaires

      call offline_bin2netcdf_step2(txt2) ! entete et ecriture du champ

      goto 7217 ! Lire un nouveau champ

 7219 continue

      status=nfmpi_close(ncid)
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)
#endif
      stop 'Binary files converted into netcdf files !!!'

      end subroutine offline_bin2netcdf

!...............................................................................

      subroutine offline_bin2netcdf_step1 !15-08-22
      implicit none

!       status=nf_redef(ncid)

      status=nfmpi_close(ncid)
      status=nfmpi_create(par%comm2d &
                         ,date_shortname//'.nc' &
                         ,nf_clobber+ NF_64BIT_OFFSET,MPI_INFO_NULL,ncid) !27-05-18
      if(status/=0)stop 'STOP offline_bin2netcdf nfmpi_create'
      count_netcdfvar=0
      loop_netcdf=0
      call netcdf_dim
      call netcdf_general_attributes(ncid)

! time:
      k0=1 ; vardim(1)=time_dim                                         ! 1D ; dim1
      texte80(1)='time'                                                 ! variable
      texte80(2)=time_units
      texte80(8)=time_units(14:33)                                      ! time_origin : kount=0
      texte80(9)='gregorian'                                            ! calendar
      texte80(3:4)='time'                                               ! long_name
      texte80(5)='T'  ; texte80(6)='time'                               ! axis ; associate
      texte80(7)='double'

      loop_netcdf=0 ! faire entete
      call netcdf_main('_t')
         status=nfmpi_enddef(ncid)
      if(status/=0)stop 'Err offline_bin2netcdf_step1 nfmpi_enddef'
         status=nfmpi_close(ncid)
      if(status/=0)stop 'Err offline_bin2netcdf_step1 nfmpi_close'

      status=nfmpi_open(par%comm2d   &
                       ,date_shortname//'.nc' &
                       ,nf_write+ NF_64BIT_OFFSET,MPI_INFO_NULL,ncid)
      if(status/=0)stop 'Err offline_bin2netcdf_step1 nfmpi_open'

      loop_netcdf=1 ! ecrire le champ
      count_netcdfvar=count_netcdfvar-1 ! (parce que netcdf_main rajoutera 1...)
      call netcdf_main('_t')

!     status=nfmpi_close(ncid)

! length of the time average
      status=nfmpi_redef(ncid)
      if(status/=0)stop 'Err offline_bin2netcdf_step1 nfmpi_redef'
      k0=1 ; vardim(1)=time_dim                                         ! 1D ; dim1
      texte80(1)='cumulativetime'
      texte80(2)='hours'
      texte80(3:4)='length_of_the_time_average'
      texte80(5)='T'  ; texte80(6)='time'                               ! axis ; associate
      texte80(7)='double'

      loop_netcdf=0 ! faire entete
      call netcdf_main('_t')
      status=nfmpi_enddef(ncid)

      loop_netcdf=1 ! ecrire le champ
      count_netcdfvar=count_netcdfvar-1 ! (parce que netcdf_main rajoutera 1...)
      call netcdf_main('_t')


      end subroutine offline_bin2netcdf_step1

!...............................................................................

      subroutine offline_bin2netcdf_step2(txt2) !15-08-22
      implicit none
      character*2 txt2

      status=nfmpi_redef(ncid)
      if(status/=0)stop 'Err offline_bin2netcdf_step2 nfmpi_redef'

      loop_netcdf=0 ! faire entete
      call netcdf_main(txt2)

      status=nfmpi_enddef(ncid)
      if(status/=0)stop 'Err offline_bin2netcdf_step2 nfmpi_enddef'

      loop_netcdf=1 ! ecrire le champ
      count_netcdfvar=count_netcdfvar-1 ! (parce que netcdf_main rajoutera 1...)
      call netcdf_main(txt2)

      end subroutine offline_bin2netcdf_step2

!...............................................................................

      subroutine offline_reversethelists(varnum_) !26-05-19
      implicit none
      integer varnum_

! Details dans
! https://docs.google.com/document/d/1jjmXtsN-rAfKzEdnYFlz3p__6ieAxbZNuvwKfu5yGck/edit

      if(par%rank==0) then !0000000000000>

      do loop1=1,varnum_

! Etape 1 faire une copie
       open(unit=4,file=texte80(loop1) &       !22-12-09
                  ,access='direct'     &
                  ,recl=recordlength_  &
                  ,form='unformatted')
       open(unit=3,file='tmp/tmpreversedlist' &
                  ,access='direct'     &
                  ,recl=recordlength_  &
                  ,form='unformatted')

        do nc=1,ofl_rec_max
           read(4,rec=nc)filename_,ogcmtimecounter_,ofl_readtime_next,ofl_period_next
           if(nc==1)          x0=ofl_readtime_next 
           if(nc==ofl_rec_max)x1=ofl_readtime_next 
          write(3,rec=nc)filename_,ogcmtimecounter_,ofl_readtime_next,ofl_period_next
        enddo ! nc

! Etape 2 changer l'ordre des noms et apporter le shift temporel expliquE dans
! https://docs.google.com/document/d/1jjmXtsN-rAfKzEdnYFlz3p__6ieAxbZNuvwKfu5yGck/edit
! le shift est elapsedtime_end-(x0+x1)+ofl_period_next
        do nc=1,ofl_rec_max
           read(4,rec=nc)filename_,ogcmtimecounter_,ofl_readtime_next,ofl_period_next
           ofl_readtime_next=ofl_readtime_next+elapsedtime_end-(x0+x1)+ofl_period_next
           read(3,rec=ofl_rec_max-nc+1)filename_
          write(4,rec=nc)filename_,ogcmtimecounter_,ofl_readtime_next,ofl_period_next
          if(ofl_readtime_next<=elapsedtime_now) then !----->
             ofl_rec_now=nc 
             ofl_period_now=ofl_period_next
          endif                                       !----->
        enddo ! nc

        close(3) ; close(4)
      enddo  ! loop1

! Si la reorganisation du fichier a reussi, au lancement de la
! simulation, le fichier physique associe correspond A la date de fin et
! reciproquement. TEST pour le verifier (commenter sinon):
!     do loop1=1,varnum_
!      open(unit=4,file=texte80(loop1) &       !22-12-09
!                 ,access='direct'     &
!                 ,recl=recordlength_  &
!                 ,form='unformatted')
!      write(100+loop1,'(a)')trim(texte80(loop1))
!      do nc=1,ofl_rec_max
!      read(4,rec=nc)filename_,ogcmtimecounter_,ofl_readtime_next,ofl_period_next
!      write(100+loop1,'(a)')trim(filename_)
!      write(100+loop1,*)ofl_readtime_next
!      enddo
!      close(4)
!     enddo


      endif                !0000000000000>

! Le proc 0 va envoyer ses valeurs aux autres proc.
      call mpi_bcast(ofl_rec_max      ,1,mpi_integer ,0,par%comm2d,ierr)
      call mpi_bcast(ofl_readtime_next,1,mpi_double_precision,0,par%comm2d,ierr)
      call mpi_bcast(ofl_rec_now      ,1,mpi_integer ,0,par%comm2d,ierr)
      call mpi_bcast(ofl_period_now   ,1,mpi_double_precision,0,par%comm2d,ierr)
      call mpi_bcast(ofl_period_next  ,1,mpi_double_precision,0,par%comm2d,ierr)

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)
#endif
!     stop ' offline_reversethelists'
      end subroutine offline_reversethelists

!...............................................................................
      end module module_offline
