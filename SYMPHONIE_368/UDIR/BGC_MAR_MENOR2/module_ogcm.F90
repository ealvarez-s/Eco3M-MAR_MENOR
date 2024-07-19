      module module_ogcm
!______________________________________________________________________
! SYMPHONIE ocean model
! release 255 - last update: 18-05-19
!______________________________________________________________________

! Les notes du developpeur:
! https://docs.google.com/document/d/1VMmzC7jXV7sycAxoFCUcOBxVbydbyXD3imgM0wAnAFQ/edit

      use module_principal ; use module_parallele ; use module_forcages 
      use module_s
      implicit none
      include 'netcdf.inc'

      double precision, dimension(:,:,:), allocatable :: ogcmvar_r8
      integer         , dimension(:,:,:), allocatable :: ogcmvar_int    !21-04-14
      integer         , dimension(:,:)  , allocatable :: ogcmgridnum !16-05-19 si plusieurs grille nmo successives (psy2, psy4 etc...)
      character(len=250) :: ogcm_grid_file_name='none'
      integer   :: ogcm_time_counter,ogcm_var_ndims     &
                  ,ogcmzoom_istr , ogcmzoom_iend    &
                  ,ogcmzoom_jstr , ogcmzoom_jend    &
                  ,ogcmzoom_istrt, ogcmzoom_iendt   &
                  ,ogcmzoom_jstrt, ogcmzoom_jendt   &
                  ,ogcmzoom_istru, ogcmzoom_iendu   &
                  ,ogcmzoom_jstru, ogcmzoom_jendu   &
                  ,ogcmzoom_istrv, ogcmzoom_iendv   &
                  ,ogcmzoom_jstrv, ogcmzoom_jendv   &
                  ,ogcm_enlarged_bounds=20          &
                  ,ogcmgridcode=1                   & !16-05-19
                  ,flag_updatebouchetrou                    !17-05-19
!     integer :: vel_id=2,trc_id=1,ssh_id=3,var_num=3
      character*60 txt_units
       real,dimension(:,:,:,:),allocatable ::    &
              correc_t
       real,dimension(:,:,:,:),allocatable ::    &
              correc_s
       real,dimension(:,:,:),allocatable ::    &
              correc_zeta
       real time_init_correc,time_end1_correc,time_end2_correc,time_end3_correc &
         ,time_end4_correc,time_end5_correc

! Generalized interpolation of ogcm outputs. Compatible with nemo and Smodel.

!...............................................................................
! Version date:     Description des modifications:
!         08/12/08: Mise en service
!         10/12/08: Ajout d'une routine produisant un melange de la
!                    bathy de depart et de la bathy de l'ogcm
!         05/03/09: adaptation ay cas psi2v3
! 2009.2  04-09-09  Fermer les fichiers netcdf
!         07-09-09  - ogcm_rec_max remplace nc1
!                    - Parallelisation
!         08-09-09  mise e jour de la routine par_gather
! 2009.3  01-10-09  utilisation des nouveaux facteurs d'echelle verticale
!         05-10-09  ajout d'un "ifdef parallele"
!         08-10-09  compatibilite avec f95: adequation des longueurs de
!                   chaines de charactere passees en argument de allocate_global
!         18-10-09  Modification du first guess pour parallelisation parfaite
!         22-10-09  lire la missing_value des fichiers
!         04-11-09  ecrire les fichiers temporaire dans repertoire tmp
!         06-11-09  chronologie des routines revues pour garantir parallelisation
!                   parfaite: "call ogcm_grille_sous_le_fond" doit se faire avant
!                   le bouchage des trous
!         12-11-09  generalisation de l'interpolation: fonctionne pour nemo et
!                   Smodel
!         13-11-09  amelioration de la routine grille_sous_le_fond. On ne
!                   pas un point masque.
! 2010.2  17-12-09  ajout d'un close(3) oublie
!         18-12-09  remplacer angle0 par un calcul deduit de lat lon
! 2010.3  09-01-10  centrage de la rotation au point z
!         13-01-10  prevoir nom de ssh cas Smodel
!         14-01-10  prevoir zones decouvertes
! 2010.5  29-01-10  debug bouche trou ssh
! 2010.7  22-02-10  adaptation de la routine aux fichiers "medoc"
!         23-02-10  Suite du point precedent
! 2010.8  08-05-10  de nouvelles possibilites de noms pour T et S
!         03-06-10  nouvelle possibilite de noms pour cas nemo
! 2010.10 16-06-10  nouveau cas (modele Somot)
! 2010.11 16-07-10  modif affichage ecran
!         29-07-10  amelioration de l'extrapolation en cas de grille "en dessous"
!                   de la grille ogcm.
! 2010.12 20-09-10  Possibilite de calcul en simple precision
!         27-09-10  "grille sous le fond" forfait pas de 100 e 10
!         29-09-10  adaptation aux derniers fichiers offline
! 2010.14 19-11-10  le bouchage des trous du cas "nemo_z" est differencie du
!                   cas "sympa": il est maintenant 3d
!         22-11-10  suite du point precedent
! 2010.18 08-02-11  Debug routine get_bathy_ogcm
! 2010.19 13-04-11  Version adaptee e la lecture des champs MERCATOR PSY4V1R3
! 2010.20 17-04-11  Calculs sur la base d'un temps en secondes
! 2010.22 28-04-11  Adapte au fait que fichiers offline Smodel n'ont pas
!                   de scale factor ni de addoffset
! 2010.24 02-09-11  Le cas nemo_z (psi4) trouve les latitudes longitudes dans un
!                   fichier de grille
!         19-11-11  Adaptation du cas particulier oe la cellule de surface est
!                   completement au dessus de la cellule de surface de l'ogcm, suite au cas
!                   de bug trouve par Claude et Thomas. On cesse de translater les
!                   niveaux (pour avoir z=0 en kmaxp1) et a la place on force z2(max_z)
!                   a etre au dessus de depth_w(kmaxp1)
!         23-11-11  - Priorite a _FillValue plutot qu'a missingvalue pour u et v
!                   - Suite du point 19-11-11:
!                      a- une routine donnant la profondeur ogcm_z dans le cas "sympa"
!                      b- l'option de recallage de ogcm_z sur 0m en surface est supprimee
!                      c- passage du ncid en argument de la routine get_dim
! 2010.25 08-02-12  lecture 'ssh_bi'
!         16-02-12  lecture 'h_w'
!         17-02-12  Compte des coordonnees ALE la routine grille_sous_lefond doit
!                   etre appelee a chaque fois dans le cas d'une imbrication S-S
!         29-02-12  nouveaux fichiers mercators...
!         04-04-12  lon lat de l'ogcm en real ou en double
!         06-04-12  ssh-bi devient ssh-ib
!         21-06-12  10 fichiers dans le notebook_obcforcing cas nemo_z
! S26     19-09-12  Affichage a l'ecran seulement si par%rank==0
!         08-12-12  debug ncid1 en ncid_ dans routine interp_ts
!         14-02-13  blindage rotation symphonie (voir dlon_di_ etc...)
!         18-02-13  gestion Fill Value
!         24-02-13  amelioration (et non debugage!) de l'interpolation de T,S en surface
!         25-02-13  suite du point precedent
!         23-06-13  Seul le proc 0 ecrit les listes binaires
!         03-07-13  ssh-ib devient ssh_ib (pas de signe "-" svp!!)
!         29-11-13  simplifier notebook_obcfocing (les fichiers netcdf renseignent
!                   eux meme les aspects temporels)
!                   zoom d'extraction adapte a la taille du domaine mpi local
!         04-12-13  modif format ascii du fichier bathy_nested
!         05-12-13  la routine ogcm_get_bathy_ogcm es appelee depuis initial_mask_bathy
!                   voir egalement notebook_bathy
!         06-12-13  traitement des bancs decouvrants
!         17-03-14  Cas nemo-z, tous les procs participent a la construction des
!                   listes binaires
!         03-04-14  supprime call mpifinalize qui empeche l'arret du code si erreur
!         17-04-14  debug cas sympa
!         21-04-14  modif dans procedure merge bathy
!         23-04-14   - stop si liste ascii vide
!                    - on envisage egalement la possibilite de double-masquage dans nemo
!         23-06-14  On reduit de max_x max_y de un pour eviter depassement de "fichiers"
!                   En effet (et helas) on applique un dimensionnement unique (pour
!                   u v t) qui conduisait dans le cas des fichiers S a sortir des bornes
!                   du fichier S pour les vitesses
!         02-07-14  Passage aux nouveaux echanges
!         17-07-14  - ajout de ogcm_enlarged_bounds
!                   - Coherence de ce parametre avec la limite d'extension du chenillard de
!                     bouchage de trous.
!                   - Reecriture (plus claire) de la boucle de lecture du fichier de bouchage
!                     algo de bouchage vertical revu
!         25-08-14  - messages ecran
!                   - cas sympa adapte a boucle do var_count_=1,var_num
!         05-10-14  lignes bizarres commentees....
!         20-10-14  tmpfile ouvert avec status 'replace'
!         25-11-14  - extrapolation T,S cas nemo_z. Hybride extrapolation horizontale/verticale
!                   - correction bug trouve par Malek
!         27-11-14  attention a ne pas systematiser le cas particulier "simed"
!         29-12-14  interpolation du masque de l'ogcm
!         28-01-15  un nouvel algo d'extrapolation de T et S cas nemo
!         29-01-15  - suite du point precedent. La comparaison sur le cas bobshelf semble plaider
!                   en faveur de l'algo initial 100% extrapolation horizontale
!                   - interpolation de la grille OGCM. Prevoir le cas ou le masque n'est pas
!                   dans le fichier grille de nemo
!         01-02-15  - ajout d'un call obc_h(0)
!                   - diviser par max(mask_t,small1) pris en defaut en i=1 etc...
!         27-02-15  Detection de listes ogcm desordonnees
!         03-03-15  Detection de listes ogcm desordonnees SUITE
!         19-05-15  un format de plus pour variable time....
!         21-06-15  si le calcul de correspondance de grilles ne converge pas car
!                   la grille ogcm est discontinue, un algo alternatif est possible
!         30-06-15  correction d'un bug "dormant" dans construction des listes
!         12-10-15  ajout du cas 1dv
!         15-10-15  Possibilite d'interpoler un champ basse frequence pour T et S
!                   afin de constituer un etat de reference pour le PGF
!         31-10-15  suite du point precedent. Plus stricte dans les zones de bancs decouvrants
!         05-11-15  methode d'extrapolation T,S nemo sous le fond: algo 3 permet un meilleur
!                   recollement sous le fond mais cette option est toutefois commentee car elle 
!                   semble avoir la propriete d'accentuer les gradients horizontaux pres du fond
!         11-11-15  ajout flag_refstate
!         13-11-15  rejet du premier nivau au dessus du fond de mercator
!         14-11-15  Algo 3 modifie
!         18-11-15  Ecrire H-Hnemo dans fichier grid.nc
!         27-11-15  continuite cas 1DV
!         28-11-15  ajout d'un debugueur de listes ogcm
!         29-11-15  ajout ogcm_1dv_geos
!         01-12-15  plus de verifications des listes
!         03-12-15  cas ou l'ogcm est en double precision
!         15-12-15  Messages ecran
!         13-02-16  accepter les listes grilles de plus d'un seul fichier
!         15-02-16  accepter une nouvelle orthographe por fillvalue
!         18-04-16  Possibilite de ne pas refaire les listes binaires si elles
!                   sont dEjA dans le repertoire tmp
!         19-04-16  Seul le rank0 teste l'existence des listes binaires dans tmp
!         28-04-16  Augmentation de la taille des reccords des listes binaires
!         20-06-16  cas particulier flag_1dv=1 et mergebathy
!         25-11-16  debug cas lon lat "1d" 
!         05-01-17  bathy nemo dans grid.nc
!         16-02-17  trim(texte..
!         07-03-17  modif sur algo bouchage vertical et avertisseur si
!                   boucle d'elargissement de recherche horizontale des trous
!                   est insuffisante
!         14-04-17  ajout initialisation salinite dans reservoir
!         16-10-17  initialisation dans la couche fusionnee
!         15-11-17  - recherche les correspondances de grille dans mask=1 uniquement
!                   - etat de reference produit par une grille s26 fusionnee
!         24-03-18  procedure d'arret adaptee A occigen
!         05-05-18  Seul par%rank=0 lit le fichier puis envoie les valeurs aux autres par%rank par mpi_bcast
!         09-05-18  procedure d'arret adaptee A occigen
!         16-05-18  - procedure d'arret adaptee A occigen (SUITE)
!                   - ajout noms de variables indeso
!         20-05-18  detection listes desordonnees du cas "sympa"
!         22-05-18  Aide A la detection d'incoherences de grille
!         06-06-18  call ogcm_grille_sous_le_fond du cas sympa est de nouveau calculE uniquement
!                   au premier passage (cas initial) 
!         25-09-18  adaptation A la norme copernicus
!         17-11-18  ajout call ogcm_stop('blabla',1 ou 2)
! 245     23-01-19  ajout merge avec bathy ogcm cas sympa
! 255     15-05-19  pour moins de cpu changement de l'ordre des operations pour
!                   le passage des conventions nemo-symphonie pour la variable z
!         16-05-19  Enchainer des fichiers NEMO-MERCATOR differents dans une meme simulation
!         18-05-19  suite point precedent
!...............................................................................
!    _________                    .__                  .__             ! m°v°m 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................


contains


      subroutine ogcm_initial
      implicit none
      integer loop_
      character*1 txt_
#ifdef synopsis
       subroutinetitle='ogcm_initial'
       subroutinedescription= &
          'OGCM (often MERCATOR) interpolation. Initial step: reads' &
       //' notebook_obcforcing, builds binary file lists, determines' &
       //' the bounds of the extracted area'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif
! Attention cette routine n'interpole pas les champs ogcm!
! Ce qu'elle fait
! 1. lecture des listes du notebook_obcforcing et production des listes binaires acces direct
! 2. calcul coordonnees du zoom d'extraction adaptees aux dimensions locales du domaine mpi
! 3. fichiers correspondance des coordonnees s et ogcm
       if (.not. allocated(correc_t)) allocate (correc_t(imax,jmax,kmax,2))
       if (.not. allocated(correc_s)) allocate (correc_s(imax,jmax,kmax,2))
       if (.not. allocated(correc_zeta)) allocate (correc_zeta(imax,jmax,2))
       correc_t=0. ; correc_s=0. ; correc_zeta=0.

      if(par%rank==0)write(6,*)'Routine ogcm_initial'

!....................
! Faire les listes binaires acces direct:
#ifdef parallele
      cpu_seconds=MPI_Wtime ( )
#endif
      if(.not.allocated(ogcm_readtime_next)) then
               allocate(ogcm_readtime_next(0:var_num)) ; ogcm_readtime_next=0.
      endif
      if(.not.allocated(ogcm_readtime_prev)) then
               allocate(ogcm_readtime_prev(0:var_num)) ; ogcm_readtime_prev=0.
      endif
      if(.not.allocated(ogcm_period_prev))  then
               allocate(ogcm_period_prev  (0:var_num)) ; ogcm_period_prev=0.
      endif
      if(.not.allocated(ogcm_period_next))  then
               allocate(ogcm_period_next  (0:var_num)) ; ogcm_period_next=0.
      endif
      if(.not.allocated(ogcm_rec_next))     then
               allocate(ogcm_rec_next     (0:var_num)) ; ogcm_rec_next=0 
      endif
      if(.not.allocated(ogcm_rec_prev))     then
               allocate(ogcm_rec_prev     (0:var_num)) ; ogcm_rec_prev=0 
      endif
      if(.not.allocated(ogcmgridnum))       then
               allocate(ogcmgridnum       (3,2)) ; ogcmgridnum=1 !16-05-19
      endif
      call ogcm_lire_les_listes
#ifdef parallele
      cpu_seconds=MPI_Wtime ( ) - cpu_seconds
#endif

!....................
! ZOOMS D'EXTRACTION
! FICHIERS HR_TO_LR
      if(par%rank==0)write(6,*)'Calcul coordonnees zoom d''extraction'
      call ogcm_grid_to_grid_pathway(1,1) ! indice 1=type grille (t,u,v) indice 2= rec number
      call ogcm_grid_to_grid_pathway(2,1)
      call ogcm_grid_to_grid_pathway(3,1)

      end subroutine ogcm_initial

!...............................................................................

      subroutine ogcm_grid_to_grid_pathway(loop_,recnum_) !15-05-19
      implicit none
      character*1 txt_
      integer loop_, recnum_

      if(loop_==1)txt_='t' ; if(loop_==2)txt_='u' ; if(loop_==3)txt_='v'

       call ogcm_get_gridfilename(txt_,recnum_) !16-05-19

! Ouvrir le fichier des parametres de grille:
       status=nf_open(trim(ogcm_grid_file_name),nf_nowrite,ncid1)
       if(status/=0)stop ' echec nf_open ogcm_grid_file_name t'

! Lecture des valeurs des dimensions maxi pour allocation memoire
       call ogcm_get_dim(ncid1,txt_)

       call ogcm_allocate(1,1,max_x,max_y,max_z) ! arg1=allouer arg2=ogcm arg3,4,5=dimensions

! Lire les longitudes & latitudes:
       ogcmzoom_istr=1 ; ogcmzoom_jstr=1
       call ogcm_get_lonlat(ncid1,txt_)

! coordonnees du zoom:
       call ogcm_zoom_coordinate(txt_,0)

! Fichier hr_to_lr:
       call ogcm_hr_to_lr_coord(loop_)

! Fermer le fichier des parametres de grille:
       status=nf_close(ncid1)
       call ogcm_allocate(2,1,0,0,0) ! arg1=desallouer arg2=ogcm

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)  !23-06-13
#endif

      end subroutine ogcm_grid_to_grid_pathway

!...............................................................................

      subroutine ogcm_interp(var_id_,nc_,t_)
      implicit none
      integer var_id_,nc_,t_
      real*4,dimension(:),allocatable ::                              &
             ogcm_lat_1d, ogcm_lon_1d
#ifdef synopsis
       subroutinetitle='ogcm_interp'
       subroutinedescription= &
          'OGCM (often MERCATOR) interpolation general driver.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      nc=nc_

      if(par%rank==0) then
        write(6,*)
        write(6,*)'interpolation champs ogcm:'
        write(6,*)
      endif

      flag_stop=0
      if(nc_   <1.or.nc_   >ogcm_rec_max)then !--------->                  !07-09-09
       if(par%rank==0) then
        write(6,*)
        write(6,*)'Probleme:'
        write(6,*)'Il n''y a pas de fichiers ogcm a la date demandee.'
        write(6,*)'Verifier l''adequation entre la date de la simu et'
        write(6,*)'les parametrages du notebook du forcage par l''ogcm'
        write(6,*)'et les listes de fichiers ogcm.'
        write(6,*)'ncmin nc_    ncmax: 1',nc_   ,ogcm_rec_max
        write(6,*)'var_id_=',var_id_
        write(6,*)'t_=',t_
        flag_stop=1 !16-05-18
       endif
      endif                               !--------->                  !07-09-09
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0) stop 'Error 321 in subroutine ogcm_interp' !09-05-18

! Ouvrir le fichier des parametres de grille:
!     call ogcm_get_gridfilename('t') ! donne ogcm_grid_file_name (character)
!     status=nf_open(trim(ogcm_grid_file_name),nf_nowrite,ncid1)
!     if(status/=0)stop ' echec nf_open ogcm_grid_file_name t'


! Detection de changement de grille pour remise A jour des parametres d'interpolation:
! https://docs.google.com/document/d/1VMmzC7jXV7sycAxoFCUcOBxVbydbyXD3imgM0wAnAFQ/edit
      call ogcm_gridchangedetection(var_id_,nc_)

! Dimensions (max):
      call ogcm_zoom_coordinate('m',1)            ! 'm' pour dim max pour allocation unique
      if(var_id_==ssh_id) then !>>>>
        call ogcm_allocate(1,1,max_x,max_y,1)     ! arg1=allouer arg2=ogcm arg3,4,5=dimensions
      else                     !>>>>
        call ogcm_allocate(1,1,max_x,max_y,max_z) ! arg1=allouer arg2=ogcm arg3,4,5=dimensions
      endif                    !>>>>

! Dimensions grille 't':
      call ogcm_zoom_coordinate('t',1)

! Fermer le fichier des parametres de grille:
!     status=nf_close(ncid1)

! Interpolation Traceurs et SSH:
      if(var_id_==trc_id) then !trctrctrctrctrc>

!---- TEMPERATURE START -------------------------------
        if(par%rank==0)write(6,*)'interpolation de la temperature'

! Ouvrir le fichier contenant t:
        call ogcm_get_varfilename(var_id_,'t',nc_) ! donne texte250
        call ogcm_interp_ts(1,t_,nc_)


!---- SALINITY START -------------------------------
        if(par%rank==0)write(6,*)'interpolation de la salinite'

! Ouvrir le fichier contenant s:
        call ogcm_get_varfilename(var_id_,'s',nc_) ! donne texte250
        call ogcm_interp_ts(2,t_,nc_)

        call set_rivers_reservoir_sal(t_) !Salinity in rivers reservoirs !14-04-17

#ifdef bidonref
! Reference state (if any)
        if(t_==0) then !refref> (initial state only)
         if(flag_refstate==1) then !---> !11-11-15
          texte250=obcfile(9,1)
          call ogcm_interp_ts(3,t_,nc_) ! first guess temref_z
          call ogcm_interp_ts(4,t_,nc_) ! first guess salref_z

! Interpolation of temref_z salref_z on s grid  gives temref_t salref_t
          id_tem=0 ; id_sal=1
          call quick_tsref_interp(0) ! return anyv3d(i,j,k,id_tem) anyv3d(i,j,k,id_sal)
          do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
           temref_t(i,j,k)=anyv3d(i,j,k,id_tem)
           salref_t(i,j,k)=anyv3d(i,j,k,id_sal)
          enddo       ; enddo              ; enddo

         endif                     !--->
        endif          !refref
#endif

      endif                    !trctrctrctrctrc>


      if(var_id_==ssh_id) then !sshsshsshsshssh>
!---- SSH START -------------------------------
        if(par%rank==0)  &
        write(6,*)'interpolation de l''elevation de la surface'

! Ouvrir le fichier contenant ssh:
        call ogcm_get_varfilename(var_id_,'e',nc_) ! donne texte250
        call ogcm_interp_ssh(t_,nc_)

      endif                    !sshsshsshsshssh>


      if(var_id_==vel_id) then !velvelvelvelvel>
!---- U VELOCITY START -------------------------------

! Interpolation des courants:
      if(par%rank==0)write(6,*)
      if(par%rank==0)write(6,*)'interpolation de u'

! Ouvrir le fichier des parametres de grille:
      call ogcm_get_gridfilename('u',1) ! donne ogcm_grid_file_name (character)
      status=nf_open(trim(ogcm_grid_file_name),nf_nowrite,ncid1)
      if(status/=0)stop ' echec nf_open ogcm_grid_file_name u'

! Lecture les valeurs des dimensions
      call ogcm_zoom_coordinate('u',1)

! Fermer le fichier des parametres de grille:
      status=nf_close(ncid1)

! Ouvrir le fichier contenant u:
      call ogcm_get_varfilename(var_id_,'u',nc_) ! donne texte250
      if(par%rank==0)write(6,'(a,a)')'fichier variable ',trim(texte250)
      status=nf_open(trim(texte250),nf_nowrite,ncid1)
      if(status.ne.0)stop 'echec fichier ogcm pour u'

      call ogcm_interp_uv(1,t_,nc_)

! Fermer le fichier contenant u:
      status=nf_close(ncid1)                                            !04-09-09

!---- V VELOCITY START -------------------------------
      if(par%rank==0)write(6,*)
      if(par%rank==0)write(6,*)'interpolation de v'

! Ouvrir le fichier des parametres de grille:
      call ogcm_get_gridfilename('v',1) ! donne ogcm_grid_file_name (character)
      status=nf_open(trim(ogcm_grid_file_name),nf_nowrite,ncid1)
      if(status/=0)stop ' echec nf_open ogcm_grid_file_name v'

! Lecture les valeurs des dimensions
      call ogcm_zoom_coordinate('v',1)

! Fermer le fichier des parametres de grille:
      status=nf_close(ncid1)

! Ouvrir le fichier contenant v:
      call ogcm_get_varfilename(var_id_,'v',nc_) ! donne texte250
      status=nf_open(trim(texte250),nf_nowrite,ncid1)
      if(status.ne.0)stop 'echec fichier ogcm pour v'

      call ogcm_interp_uv(2,t_,nc_)

      if(par%rank==0)write(6,*)
      if(par%rank==0)write(6,*)'courant moyen'
      call ogcm_get_vmeanobc(t_)

! Fermer le fichier contenant v:
      status=nf_close(ncid1)                                            !04-09-09

      endif                    !velvelvelvelvel>

      call ogcm_allocate(2,1,0,0,0) ! arg1=desallouer arg2=ogcm

      end subroutine ogcm_interp
!------------------------------------------------------------------------------

      subroutine ogcm_1dv_geos(t_)
      implicit none
      integer t_
#ifdef synopsis
       subroutinetitle='ogcm_1dv_geos'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Modele 1DV: tableau velobc chargE avec le courant geostrophique !29-11-15
      id_tem=0 ; id_sal=1 ; id_rhf=2
      do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
       anyv3d(i,j,k,id_tem)=temobc_t(i,j,k,t_)
       anyv3d(i,j,k,id_sal)=salobc_t(i,j,k,t_)
      enddo       ; enddo         ; enddo
! Attention ne pas calculer EOS en i=j=1 et en imax jmax car les lon lat de
! ces points sont fausses du fait de la cyclicite mpi du domaine 1DV
! On ne calcule de gradients qu'entre i=2 et 3 et j=2 et 3
      call equation_of_state_full_anyv(2,3,2,3)



! ATTENTION LE POINT (i,j)=(3,3) EST COHERENT AVEC ogcm_1dv_0grd
      i=3 ; j=3 ; k=kmax+1


! Geostrophic current at depth_t(k) from depth_w(k)

! Barotropic geostrophic current:
      velobc_u(i,j,k,t_)=-grav*(sshobc_w(i,j,t_)-sshobc_w(i,j-1,t_))/( coriolis_t(i,j)*dy_t(i,j) )
      velobc_v(i,j,k,t_)=+grav*(sshobc_w(i,j,t_)-sshobc_w(i-1,j,t_))/( coriolis_t(i,j)*dx_t(i,j) )

! Geostrophic current at depth_w(k) ("thermal wind relation")
      x0=grav/(rho*coriolis_t(i,j)*dy_t(i,j))
      x1=grav/(rho*coriolis_t(i,j)*dx_t(i,j))
      do k=kmax,1,-1

       velobc_u(i,j,k,t_)=velobc_u(i,j,k+1,t_)  &
       -x0*dz_t(i,j,k,1)*(anyv3d(i,j,k,id_rhf)-anyv3d(i,j-1,k,id_rhf))

       velobc_v(i,j,k,t_)=velobc_v(i,j,k+1,t_)  &
       +x1*dz_t(i,j,k,1)*(anyv3d(i,j,k,id_rhf)-anyv3d(i-1,j,k,id_rhf))

      enddo

! Geostrophic current at depth_t(k) from depth_w(k)
      do k=1,kmax
       velobc_u(i,j,k,t_)=0.5*(velobc_u(i,j,k,t_)+velobc_u(i,j,k+1,t_))
       velobc_v(i,j,k,t_)=0.5*(velobc_v(i,j,k,t_)+velobc_v(i,j,k+1,t_))
      enddo
     
      end subroutine ogcm_1dv_geos

!------------------------------------------------------------------------------
      subroutine ogcm_1dv_0grd(t_) !29-11-15
      implicit none
      integer t_
#ifdef synopsis
       subroutinetitle='ogcm_1dv_0grd'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      do k=2,var_num ! debuger
       if(ogcm_readtime_next(k)/=ogcm_readtime_next(k-1)) then
        write(6,*)' Err OGCM variables not synchronous'
        stop 'ogcm_1dv_0grd'
       endif 
      enddo

! ATTENTION LE POINT (i,j)=(3,3) EST COHERENT AVEC ogcm_1dv_geos

! Tracers:
        ub4=ubound(temobc_t) ; lb4=lbound(temobc_t)
        do k=lb4(3),ub4(3) ; do j=lb4(2),ub4(2) ; do i=lb4(1),ub4(1)      
         temobc_t(i,j,k,t_)=temobc_t(3,3,k,t_)
         salobc_t(i,j,k,t_)=salobc_t(3,3,k,t_)
        enddo ; enddo ; enddo

#ifdef bidonref
        if(allocated(temref_t)) then
         ub3=ubound(temref_t) ; lb3=lbound(temref_t)
         do k=lb3(3),ub3(3) ; do j=lb3(2),ub3(2) ; do i=lb3(1),ub3(1)      
          temref_t(i,j,k)=temref_t(3,3,k) !27-11-15
          salref_t(i,j,k)=salref_t(3,3,k)
         enddo ; enddo ; enddo
        endif
#endif

! Velocities:
        ub4=ubound(velobc_u) ; lb4=lbound(velobc_u)
        do k=lb4(3),ub4(3) ; do j=lb4(2),ub4(2) ; do i=lb4(1),ub4(1)      
         velobc_u(i,j,k,t_)=velobc_u(3,3,k,t_)
        enddo ; enddo ; enddo

        ub4=ubound(velobc_v) ; lb4=lbound(velobc_v)
        do k=lb4(3),ub4(3) ; do j=lb4(2),ub4(2) ; do i=lb4(1),ub4(1)      
         velobc_v(i,j,k,t_)=velobc_v(3,3,k,t_)
        enddo ; enddo ; enddo

        ub3=ubound(velbarobc_u) ; lb3=lbound(velbarobc_u)
        do j=lb3(2),ub3(2) ; do i=lb3(1),ub3(1)      
         velbarobc_u(i,j,t_)=velbarobc_u(3,3,t_)
        enddo ; enddo 

        ub3=ubound(velbarobc_v) ; lb3=lbound(velbarobc_v)
        do j=lb3(2),ub3(2) ; do i=lb3(1),ub3(1)      
         velbarobc_v(i,j,t_)=velbarobc_v(3,3,t_)
        enddo ; enddo 

! ssh:
        ub3=ubound(sshobc_w) ; lb3=lbound(sshobc_w)
        do j=lb3(2),ub3(2) ; do i=lb3(1),ub3(1)      
         sshobc_w(i,j,t_)=sshobc_w(3,3,t_)
        enddo ; enddo 


      end subroutine ogcm_1dv_0grd

!------------------------------------------------------------------------------

      subroutine ogcm_find_zmin_zmax_Smodel
      use module_principal
      use module_parallele !#MPI
      implicit none
      double precision zmaxglob,zminglob !#MPI
#ifdef synopsis
       subroutinetitle='ogcm_find_zmin_zmax_Smodel'
       subroutinedescription= &
      'Finds min and max values of depth_w the depth of grid points'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Pour etre sur que la grille OGCM encadre bien verticalement la grille Smodel
! on repositionne (de maniere conservative) les niveaux extremes de l'ogcm:

      zmax=-1.d10
      zmin= 1.d10
      do j=1,jmax
      do i=1,imax
       if(  mask_t(i,j,kmaxp1).eq.1) then
        do k=kmin_w(i,j),kmaxp1
         zmax=max(zmax, depth_w(i,j,k))
         zmin=min(zmin, depth_w(i,j,k))
        enddo
       endif
      enddo
      enddo
      zmax=zmax+0.001
      zmin=zmin-0.001


#ifdef parallele
      call mpi_allreduce(zmin,zminglob,1,mpi_double_precision,          &
           mpi_min,par%comm2d ,ierr)
      zmin=zminglob
      call mpi_allreduce(zmax,zmaxglob,1,mpi_double_precision,          &
           mpi_max,par%comm2d ,ierr)
      zmax=zmaxglob
#endif

      end subroutine ogcm_find_zmin_zmax_Smodel

!------------------------------------------------------------------------------

      subroutine ogcm_adapt_zmin_zmax_ogcm
      use module_principal
      use module_forcages
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='ogcm_adapt_zmin_zmax_ogcm'
       subroutinedescription= &
          'For convenience the surface level of the OGCM is shifted'   &
       //' upward if lower than zmax the highest level of the S model.'&
       //' The OGCM variable of the surface level is thus possibly'    &
       //' ajusted using a linear extrapolation'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0)write(6,*)'repositionne la surface de l''ogcm'
! Le niveau le plus haut de l'ogcm sera au-dessus de celui de modele
! La variable corespondante de l'ogcm est recalculee.

       do j=1,max_y
       do i=1,max_x

        if(zmax.gt.ogcm_z(i,j,max_z)) then !----->

          ogcm_var3d(i,j,max_z)=                                        &
          ogcm_var3d(i,j,max_z)                                         &
       +( ogcm_var3d(i,j,max_z)-ogcm_var3d(i,j,max_z-1) )               &
       /(     ogcm_z(i,j,max_z)-    ogcm_z(i,j,max_z-1) )               &
       *(                zmax      -ogcm_z(i,j,max_z)   )

         ogcm_z(i,j,max_z)=zmax

        endif                              !----->

!       if(zmin.lt.ogcm_z(i,j,1    )) then !----->

!        ogcm_var3d(i,j,1)=
!    &   ogcm_var3d(i,j,1)
!    & +(ogcm_var3d(i,j,1)-ogcm_var3d(i,j,2) )
!    & /(    ogcm_z(i,j,1)    -ogcm_z(i,j,2) )
!    & *(             zmin    +ogcm_z(i,j,1) )

!        ogcm_z(i,j,1    )=zmin

!       endif                              !----->

       enddo
       enddo

      return
      end subroutine ogcm_adapt_zmin_zmax_ogcm

!------------------------------------------------------------------------------

      subroutine ogcm_fichier_bouchetrou
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
#ifdef synopsis
       subroutinetitle='ogcm_fichier_bouchetrou'
       subroutinedescription= &
          'Finds nearest replacing sea grid points for OGCM land grid' &
       //' points concerned by the interpolation process. Stores'      &
       //' coordinates in "bouchetrou" file.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      open(unit=4,file=trim(tmpdirname)//'bouchetrou_ssh'//dom_c//'.out')   !#MPI

      k=max_z
      do 1862 j=1,max_y
      do 1862 i=1,max_x

       if(ogcm_var3d(i,j,k).lt.-1.e9)then !eeeeeeeeeeeeeee>
        ksecu=0
        i10=1
        dist1=1e10
 1864   continue

        do 1863 k1=0,1

         j0=max0(j-i10,1)
         j2=min0(j+i10,max_y)
         j3=k1*(j2-j0)+(1-k1)

         i0=max0(i-i10+k1,1)
         i2=min0(i+i10-k1,max_x)
         i3=(1-k1)*(i2-i0)+k1

         do 1861 j1=j0,j2,j3
         do 1861 i1=i0,i2,i3
          if(abs(ogcm_var3d(i1,j1,k)).lt.1.e9)then !%%%%%%%%%%%%%>
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

      return
      end subroutine ogcm_fichier_bouchetrou

!------------------------------------------------------------------------------
#ifdef bidon
      subroutine ogcm_appli_bouchetrou_uv(ichoix)
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
      integer ichoix
#ifdef synopsis
       subroutinetitle='ogcm_appli_bouchetrou_uv'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


      if(par%rank==0)write(6,*)'bouchage des trous'

      if(ichoix.eq.1)then
       open(unit=3,file='bouchetrou_u'//dom_c//'.out')     !#MPI
       open(unit=4,file='bouchesurf_u'//dom_c//'.out')     !#MPI
      endif
      if(ichoix.eq.2)then
       open(unit=3,file='bouchetrou_v'//dom_c//'.out')     !#MPI
       open(unit=4,file='bouchesurf_v'//dom_c//'.out')     !#MPI
      endif

      do j1=1,max_y
      do i1=1,max_x
       read(4,*,end=1761)i,j
       ogcm_var3d(i,j,1)=0.
      enddo
      enddo
 1761 continue
      close(4)

      do k1=1,max_z
      do j1=1,max_y
      do i1=1,max_x
       read(3,*,end=1760)i,j,k2,i4,j4,k3
       ogcm_var3d(i,j,k2)=ogcm_var3d(i4,j4,k3)
      enddo
      enddo
      enddo
 1760 continue
      close(3)

      return
      end subroutine ogcm_appli_bouchetrou_uv
#endif
!------------------------------------------------------------------------------

      subroutine ogcm_appli_bouchetrou
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
#ifdef synopsis
       subroutinetitle='ogcm_appli_bouchetrou'
       subroutinedescription= &
          'Reads "bouchetrou" file and replaces land values'    &
       //' concerned by the interpolation process with nearest' &
       //' OGCM sea values. 3D arrays.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0)write(6,*)'bouchage des trous'

!     open(unit=3,file='bouchetrou_t'//dom_c//'.out')     !#MPI
      open(unit=3,file=trim(tmpdirname)//'bouchetrou_ssh'//dom_c//'.out')   !#MPI
      do j1=1,max_y
      do i1=1,max_x
       read(3,*,end=1760)i,j,i4,j4
       do k=1,max_z
        ogcm_var3d(i,j,k)=ogcm_var3d(i4,j4,k)
            ogcm_z(i,j,k)=    ogcm_z(i4,j4,k)
       enddo
      enddo
      enddo
!     ENDDO
 1760 continue
      close(3)

      return
      end subroutine ogcm_appli_bouchetrou

!------------------------------------------------------------------------------

      subroutine ogcm_appli_bouchetrou_ssh
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
#ifdef synopsis
       subroutinetitle='ogcm_appli_bouchetrou_ssh'
       subroutinedescription= &
          'Reads "bouchetrou" file and replaces land values'    &
       //' concerned by the interpolation process with nearest' &
       //' OGCM sea values. 2D arrays.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0)write(6,*)'bouchage des trous'

!     open(unit=3,file='bouchetrou_ssh.out')
      open(unit=3,file=trim(tmpdirname)//'bouchetrou_ssh'//dom_c//'.out')     !#MPI
      do j1=1,max_y
      do i1=1,max_x
       read(3,*,end=1760)i,j,i4,j4
       ogcm_var2d(i,j)=ogcm_var2d(i4,j4)
      enddo
      enddo
 1760 continue
      close(3)

      return
      end subroutine ogcm_appli_bouchetrou_ssh
!------------------------------------------------------------------------------

      subroutine ogcm_interp_ssh(t_   ,nc_   )
      use module_principal
      use module_parallele
      use module_forcages
      implicit none
      include 'netcdf.inc'
      integer t_   ,nc_
#ifdef synopsis
       subroutinetitle='ogcm_interp_ssh'
       subroutinedescription= &
          'Reads SSH in the OGCM file and interpolates on the S grid.' &
       //' Stores the result in sshobc_w'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      k=index(texte250," ")
      if(par%rank==0)write(6,'(a,a)')'fichier variable ',trim(texte250)
! lire l'identifiant du fichier:
!     status=nf_open(texte250(1:k-1),nf_nowrite,ncid1)
      status=nf_open(trim(texte250),nf_nowrite,ncid1)
      if(status.ne.0)stop 'echec ouverture fichier'

      if(par%rank==0)write(6,*)'About to load SSH'
                    status=nf_inq_varid(ncid1,'sossheig',var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'SOSSHEIG',var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'surf_el',var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'ssh',var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'SSH',var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'ssh_bi',var_id)      !08-02-12
       if(status/=0)status=nf_inq_varid(ncid1,'ssh-bi',var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'ssh-ib',var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'ssh_ib',var_id)   !03-07-13
       if(status/=0)status=nf_inq_varid(ncid1,'zos',var_id)   !25-09-18
       if(status/=0)stop 'interp_ssh erreur sur nom ssh_ib'

!-- FACTEUR & BIAIS ->
       status=nf_get_att_real(ncid1,var_id,'scale_factor',var_scalefactor)
       if(status/=0)var_scalefactor=1.
       status=nf_get_att_real(ncid1,var_id,'add_offset',var_addoffset)
       if(status/=0)var_addoffset=0.
       if(par%rank==0)write(6,*)var_scalefactor,var_addoffset

!-- LECTURE SSH
      status=nf_inq_varndims(ncid1,var_id,ogcm_var_ndims)
      varstart(1:2)=1 ; varcount(1)=max_x ; varcount(2)=max_y
      varstart(3)=ogcm_time_counter ; varcount(3)=1

      varstart(1)=ogcmzoom_istr     ; varcount(1)=max_x
      varstart(2)=ogcmzoom_jstr     ; varcount(2)=max_y
      varstart(3)=ogcm_time_counter ; varcount(3)=1

      if(par%rank==0)write(6,*)'SSH =ogcm_time_counter=',ogcm_time_counter

      status=nf_inq_vartype(ncid1,var_id,k0)
      if(par%rank==0)write(6,*)"SSH nf_inq_vartype k0=",k0
      ksecu=0
      if(k0==5) then !........>
       ksecu=1
       status=nf_get_vara_real(ncid1,var_id                &
                              ,varstart(1:ogcm_var_ndims)  &
                              ,varcount(1:ogcm_var_ndims)  &
                            ,ogcm_var2d(1:max_x,1:max_y))
      endif          !........>
      if(k0==3) then !........>
       ksecu=1
       status=nf_get_vara_int(ncid1,var_id                &
                             ,varstart(1:ogcm_var_ndims)  &
                             ,varcount(1:ogcm_var_ndims)  &
                              ,short2d(1:max_x,1:max_y))
       do j=1,max_y
       do i=1,max_x
        ogcm_var2d(i,j)=var_scalefactor*short2d(i,j)+var_addoffset
!       if(j==jmax/2.and.short2d(i,j))write(66,*)i,ogcm_var2d(i,j),short2d(i,j)
       enddo
       enddo
      endif          !........>
      if(status/=0)stop 'error nf_get_vara for ogcm ssh'
      if(ksecu==0) stop 'interp_ssh type non identifie'

! remplir les trous:
      call ogcm_appli_bouchetrou_ssh
      if(par%rank==0)write(6,*)'calcul interpolation'

      open(unit=3, &
           file=trim(tmpdirname)//'hr_to_lr_coord_1_'//dom_c//'.out')

! Interpolation:
      do j=0,jmax+1
      do i=0,imax+1

!     read(3,*)i0,j0,deci,decj,x0,rapi,rapj
      read(3,*)i0,j0,i1,j1,x0,rapi,rapj
      if(i0.ne.i)stop 'erreur1 sur i'
      if(j0.ne.j)stop 'erreur1 sur j'

      if(  mask_t(i,j,kmaxp1).eq.1) then !>>>>>>>>>>>>

!      i1=int(deci)
!      rapi=deci-i1
!      j1=int(decj)
!      rapj=decj-j1

       sshobc_w(i,j,t_   )=                                             &
              (1.-rapi)*(1.-rapj)*ogcm_var2d(i1  ,j1  )                 &
             +(1.-rapi)*    rapj *ogcm_var2d(i1  ,j1+1)                 &
             +    rapi *(1.-rapj)*ogcm_var2d(i1+1,j1  )                 &
             +    rapi *    rapj *ogcm_var2d(i1+1,j1+1)

!       sshobc_w(i,j,t_)=max(real(sshobc_w(i,j,t_),kind=8)),1.d-3*wetdry_cst3-h_w(i,j))!08-05-14

      endif                         !>>>>>>>>>>>>>>
      enddo   ! fin J
      enddo   ! Fin I

      close(3)

! Fermer le fichier contenant t, s et ssh:
      status=nf_close(ncid1)                                            !04-09-09


!#ifdef parallele
!      lb3=lbound(sshobc_w) ; ub3=ubound(sshobc_w)
!      call echange('z0',sshobc_w,lb3,ub3,t_)
!#endif
!mpi: echange "z0" sur sshobc_w(:,:,t_):
      call obc_sshobc_mpi(t_) !mpi: echange "z0" sur sshobc_w(:,:,t_) !02-07-14

      end subroutine ogcm_interp_ssh

!------------------------------------------------------------------------------

      subroutine ogcm_transformer_1dv_en_3d
      use module_principal
      use module_forcages
      implicit none
#ifdef synopsis
       subroutinetitle='ogcm_transformer_1dv_en_3d'
       subroutinedescription= &
          'Transforms 1DV depth array (NEMO case) into a 3D array'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Charger tout le tableau 3d:
      do k=1,max_z
      do j=1,max_y
      do i=1,max_x
       ogcm_z(i,j,k)=ogcm_z(1,1,k)
      enddo
      enddo
      enddo

      end subroutine ogcm_transformer_1dv_en_3d

!------------------------------------------------------------------------------

      subroutine ogcm_changer_le_signe_de_z
      use module_principal
      use module_forcages
      implicit none
#ifdef synopsis
       subroutinetitle='ogcm_changer_le_signe_de_z'
       subroutinedescription= &
          'Apply a negative sign to depth levels under the' &
       //'  sea surface (in NEMO case only)'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      do k=1,max_z
      do j=1,max_y
      do i=1,max_x
       ogcm_z(i,j,k)=-ogcm_z(i,j,k)
      enddo
      enddo
      enddo

      end subroutine ogcm_changer_le_signe_de_z

!------------------------------------------------------------------------------

      subroutine ogcm_renverser_axe_vertical_z
      use module_principal
      use module_forcages
      implicit none
#ifdef synopsis
       subroutinetitle='ogcm_changer_le_signe_de_z'
       subroutinedescription= &
          'Vertical indexes of z array are re-ordered respecting an' &
       //' increasing upward convention. Practically from 1 at the' &
       //' the bottom to max_z at the surface. NEMO case only.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      do k=1,max_z/2

       do j=1,max_y
       do i=1,max_x

         x1=ogcm_z(i,j,k)
            ogcm_z(i,j,k)=ogcm_z(i,j,max_z+1-k)
                          ogcm_z(i,j,max_z+1-k)=x1

       enddo
       enddo

      enddo

      end subroutine ogcm_renverser_axe_vertical_z

!------------------------------------------------------------------------------

      subroutine ogcm_renverser_axe_vertical_var
      use module_principal
      use module_forcages
      implicit none
#ifdef synopsis
       subroutinetitle='ogcm_renverser_axe_vertical_var'
       subroutinedescription= &
          'Vertical indexes of the variable array are re-ordered'   &
       //' respecting an increasing upward convention. Practically' &
       //' from 1 at the bottom to max_z at the surface. NEMO case' &
       //' only.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      do k=1,max_z/2

       do j=1,max_y
       do i=1,max_x

         x1=ogcm_var3d(i,j,k)
            ogcm_var3d(i,j,k)=ogcm_var3d(i,j,max_z+1-k)
                              ogcm_var3d(i,j,max_z+1-k)=x1

       enddo
       enddo

      enddo

      end subroutine ogcm_renverser_axe_vertical_var

!------------------------------------------------------------------------------
#ifdef bidon
      subroutine ogcm_top_level_at_zero_meter
      use module_principal
      use module_forcages
      implicit none

      do k=1,max_z
      do j=1,max_y
      do i=1,max_x
       ogcm_z(i,j,k)=ogcm_z(i,j,k)-ogcm_z(i,j,max_z)
      enddo
      enddo
      enddo

      end subroutine ogcm_top_level_at_zero_meter
#endif
!------------------------------------------------------------------------------

      subroutine ogcm_bottom_boundary_condition
      use module_principal
      use module_forcages
      implicit none
#ifdef synopsis
       subroutinetitle='ogcm_bottom_boundary_condition'
       subroutinedescription=                                          &
          'Values under the bottom (identified as such with the filval'&
       //' value) are replaced by the first sea value above the sea'   &
       //' bottom'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Avec cette nouvelle routine generique, les points "sous le fond" de
! l'ogcm, y compris dans le cas mercator_ingv, sont quasi-colles, de
! sorte que l'extrapolation sous le fond suivra la meme procedure, que
! l'ogcm soit en coordonnee z ou en coordonnee sigma:
      do j=1,max_y
      do i=1,max_x
       do k=max_z-1,1,-1
        if(ogcm_var3d(i,j,k)==filval) then

             ogcm_z(i,j,k)=    ogcm_z(i,j,k+1)
         ogcm_var3d(i,j,k)=ogcm_var3d(i,j,k+1)

        endif
       enddo
      enddo
      enddo

      end subroutine ogcm_bottom_boundary_condition

!------------------------------------------------------------------------------

      subroutine ogcm_interp_ts(ichoix,t_   ,nc_   )
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
      include 'netcdf.inc'
      double precision zup_,zdw_,vup_,vdw_,dzsurf_
      integer ichoix,t_,nc_,ncid_,valmin_,valmax_
#ifdef synopsis
       subroutinetitle='ogcm_interp_ts'
       subroutinedescription= &
          'Reads T or S field in the netcdf OGCM file and interpolates' &
       //' on the S grid and stores in temobc_t ans salobc_t'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      k=index(texte250," ")
      if(par%rank==0)write(6,'(a,a)')'fichier variable:',trim(texte250)
!     write(6,'(a,a)')'fichier variable:',texte250(1:k-1)
! lire l'identifiant du fichier:
!     status=nf_open(texte250(1:k-1),nf_nowrite,ncid_)
      status=nf_open(trim(texte250),nf_nowrite,ncid_)
      if(status/=0) then !15-12-15
        write(6,'(a,i0,a,a)')'par%rank=',par%rank   &
                            ,' Err nf_open ',texte250(1:k-1)
        stop 'Error nf_open file1 in ogcm_interp_ts'
      endif

      if(ichoix==1) then !1111111111111111111>

      if(obc_ogcm_type(1:6)=='nemo_z'.or.   &
         obc_ogcm_type(1:4)=='ncom')then !nenenenenenenenenenenenenenenenennenenenenen>

      texte30='deptht'
      status=nf_inq_varid(ncid_,trim(texte30),var_id)
      if(status/=0)then                 ! 03-06-10
      texte30='gdept'                   ! 03-06-10
      status=nf_inq_varid(ncid_,trim(texte30),var_id)     ! 03-06-10
      endif
      if(status/=0)then
      texte30='depth'                 ! ncom
      status=nf_inq_varid(ncid_,trim(texte30),var_id)     ! 03-06-10
      endif                                               ! 03-06-10
      if(status/=0)stop 'interp_ts erreur sur nom variable profondeur'

      status=nf_inq_vartype(ncid_,var_id,k0)
      !print *,"ogcm_z type k0=",k0,trim(texte30)
      if((k0/=5).and.(k0/=6))  &
       stop 'interp_ts n''a pas prevu ogcm_z non real'

! si cartesien z est un tableau 1DV:
      status=nf_get_var_real(ncid_,var_id,ogcm_z(1:1,1:1,1:max_z))

! passage conventions nemo-symphonie:
! changer le signe !15-05-19
       ogcm_z(1,1,1:max_z)=-ogcm_z(1,1,1:max_z)
! axe vertical pointant vers le haut !15-05-19
       do k=1,max_z/2
         x1=ogcm_z(1,1,k)
            ogcm_z(1,1,k)=ogcm_z(1,1,max_z+1-k)
                          ogcm_z(1,1,max_z+1-k)=x1
       enddo

! Le tableau de profondeur 1DV est transforme en un tableau 3D:
      call ogcm_transformer_1dv_en_3d

! On recalle les surfaces des 2 model_s l'une par rapport e l'autre:
!     call top_level_at_zero_meter ! commenter le 23-11-11

      endif                                !nenenenenenenenenenenenenenenenenenenenenen>

      if(obc_ogcm_type(1:5)=='sympa') then
       call ogcm_get_z_sympa(nc_   ,'t')
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)  !23-06-13
#endif
! comme routine precedente ecrase max_x max_y on rappelle le calcul de max_x max_y
       call ogcm_zoom_coordinate('t',1)
      endif

! On recalle les surfaces des 2 model_s l'une par rapport e l'autre:
!     call top_level_at_zero_meter ! commenter le 23-11-11

      endif              !1111111111111111111>


      if(ichoix==1.and.par%rank==0) &
      write(6,*)'About to read ogcm temperature'
      if(ichoix==2.and.par%rank==0) &
      write(6,*)'About to read ogcm salinity'
      if(par%rank==0)write(6,*)'T S ogcm_time_counter=',ogcm_time_counter
      if(ichoix==3.and.par%rank==0) &
      write(6,*)'About to read ogcm temperature for reference state'
      if(ichoix==4.and.par%rank==0) &
      write(6,*)'About to read ogcm salinity for reference state'

!-- IDENTIFIANT T ou S ->
!     if(ichoix==1.and.obc_ogcm_type(1:6)=='nemo_z') then
      if(obc_ogcm_type(1:6)=='nemo_z') then !nemonemo>

        if(ichoix==1.or.ichoix==3) then
                      status=nf_inq_varid(ncid_,'votemper',var_id)
         if(status/=0)status=nf_inq_varid(ncid_,'temperature',var_id) !03-12-15
         if(status/=0)status=nf_inq_varid(ncid_,'thetao',var_id) !25-09-18
         if(status/=0)stop 'interp_ts erreur sur nom variable T'
        endif
        if(ichoix==2.or.ichoix==4) then
                      status=nf_inq_varid(ncid_,'vosaline',var_id)
         if(status/=0)status=nf_inq_varid(ncid_,'salinity',var_id) !03-12-15
         if(status/=0)status=nf_inq_varid(ncid_,'so',var_id) !25-09-18
         if(status/=0)stop 'interp_ts erreur sur nom variable S'
        endif

      endif                                 !nemonemo>

      if(ichoix==1.and.obc_ogcm_type(1:4)=='ncom') then
        texte30='water_temp'
        status=nf_inq_varid(ncid_,trim(texte30),var_id)
        if(status/=0)stop 'interp_ts erreur sur nom variable water_temp'
      endif
      if(ichoix==2.and.obc_ogcm_type(1:4)=='ncom') then
       texte30='salinity'
       status=nf_inq_varid(ncid_,trim(texte30),var_id)
       if(status/=0)stop 'interp_ts erreur sur nom ncom salinity'
      endif

      if(obc_ogcm_type(1:5)=='sympa') then !symphonie> !15-11-17
       if(ichoix==1.or.ichoix==3) then
        texte30='t'
        status=nf_inq_varid(ncid_,trim(texte30),var_id)
        if(status/=0)status=nf_inq_varid(ncid_,'T',var_id)
        if(status/=0)status=nf_inq_varid(ncid_,'tem',var_id)             !07-05-10
        if(status/=0)status=nf_inq_varid(ncid_,'temp',var_id)            !07-05-10
        if(status/=0)status=nf_inq_varid(ncid_,'temperature',var_id)     !07-05-10
        if(status/=0)stop 'interp_ts err sur nom variable temperature'
       endif
       if(ichoix==2.or.ichoix==4) then
        texte30='s'
        status=nf_inq_varid(ncid_,trim(texte30),var_id)
        if(status/=0)status=nf_inq_varid(ncid_,'S',var_id)
        if(status/=0)status=nf_inq_varid(ncid_,'sal',var_id)             !07-05-10
        if(status/=0)status=nf_inq_varid(ncid_,'sali',var_id)            !07-05-10
        if(status/=0)status=nf_inq_varid(ncid_,'salinity',var_id)        !07-05-10
        if(status/=0)stop 'interp_ts erreur sur nom variable salinity'
       endif
      endif                                !symphonie> !15-11-17

!-- FACTEUR & BIAIS ->
       status=nf_get_att_real(ncid_,var_id,'scale_factor',var_scalefactor)
       if(status/=0)var_scalefactor=1.
       status=nf_get_att_real(ncid_,var_id,'add_offset',var_addoffset)
       if(status/=0)var_addoffset=0.
       if(par%rank==0)write(6,*)var_scalefactor,var_addoffset
       status=nf_get_att_text(ncid_,var_id,'units',txt_units)
       if(status/=0)txt_units='none'
       if(txt_units(1:1)=='K')var_addoffset=var_addoffset-273.15

!-- FILL VALUE ->
                   status=nf_get_att_real(ncid_,var_id,'_FillValue',filval)   !16-06-10
      if(status/=0)status=nf_get_att_real(ncid_,var_id,'_Fillvalue',filval)   !15-02-16
      if(status/=0)status=nf_get_att_real(ncid_,var_id,'missing_value',filval)

!-- LECTURE T ou S -->
      status=nf_inq_varndims(ncid_,var_id,ogcm_var_ndims)
      if(status/=0) then ! debug
       write(10+par%rank,*)'ichoix=',ichoix
       write(10+par%rank,'(a,a)')'FILE: ',trim(texte250)
       write(10+par%rank,*)'ncid_',ncid_
       write(10+par%rank,*)'var_id',var_id
       write(10+par%rank,*)'ogcm_var_ndims',ogcm_var_ndims
       stop 'Err 1287 module_ogcm nf_inq_varndims see fort file'
      endif              ! debug
      varstart(1)=ogcmzoom_istr     ; varcount(1)=max_x
      varstart(2)=ogcmzoom_jstr     ; varcount(2)=max_y
      varstart(3)=1                 ; varcount(3)=max_z
      varstart(4)=ogcm_time_counter ; varcount(4)=1

      status=nf_inq_vartype(ncid_,var_id,k0)
      if(status/=0)stop 'Err 1293 nf_inq_vartype module_ogcm'
      ksecu=0
      if(k0==nf_real.or.k0==nf_double)then !555555555> !03-12-15
       ksecu=1
       if(k0==nf_real)                                              & !03-12-15
       status=nf_get_vara_real(ncid_,var_id                         &
                              ,varstart(1:ogcm_var_ndims)           &
                              ,varcount(1:ogcm_var_ndims)           &
                              ,ogcm_var3d(1:max_x,1:max_y,1:max_z))

        if(status/=0) then ! debug
           if(varstart(1)-1+varcount(1)>max_x) &
           write(10+par%rank,*)'depassement de memoire dim x'
           if(varstart(2)-1+varcount(2)>max_y) &
           write(10+par%rank,*)'depassement de memoire dim y'
           if(varstart(3)-1+varcount(3)>max_z) &
           write(10+par%rank,*)'depassement de memoire dim z'
           write(10+par%rank,*)'ogcm_var_ndims',ogcm_var_ndims 
           write(10+par%rank,*)'varstart(1:ogcm_var_ndims)',varstart(1:ogcm_var_ndims)
           write(10+par%rank,*)'varcount(1:ogcm_var_ndims)',varcount(1:ogcm_var_ndims)
           write(10+par%rank,*)'max_x,max_y,max_z',max_x,max_y,max_z
           write(10+par%rank,*)'ichoix',ichoix
           write(10+par%rank,*)'t_',t_
           write(10+par%rank,'(a,a)')'FILE: ',trim(texte250)
           stop 'Err 1302 nf_get_vara_real see fort files'
        endif              ! debug

       if(k0==nf_double) then !dddd> !03-12-15
       if(allocated(ogcmvar_r8))deallocate(ogcmvar_r8)
          allocate (ogcmvar_r8(1:max_x,1:max_y,1:max_z))
          status=nf_get_vara_double(ncid_,var_id            &
                              ,varstart(1:ogcm_var_ndims)   &
                              ,varcount(1:ogcm_var_ndims)   &
                    ,ogcmvar_r8(1:max_x,1:max_y,1:max_z))
                     ogcm_var3d(1:max_x,1:max_y,1:max_z)    &
                    =ogcmvar_r8(1:max_x,1:max_y,1:max_z)
          deallocate(ogcmvar_r8)
       endif                  !dddd>

! Cas bizarre des fichiers simed possedant un double masque:
!     do k=1,max_z ; do j=1,max_y ; do i=1,max_x
!      if(ogcm_var3d(i,j,k)==0.)ogcm_var3d(i,j,k)=filval
!     enddo        ; enddo        ; enddo

      if(var_scalefactor/=1.or.var_addoffset/=0.) then !--->
       do k=1,max_z
       do j=1,max_y
       do i=1,max_x
        if(ogcm_var3d(i,j,k)/=filval) &                               !23-04-14
           ogcm_var3d(i,j,k)=var_scalefactor*ogcm_var3d(i,j,k)+var_addoffset
       enddo
       enddo
       enddo
      endif                                            !--->
      endif         !555555555>
      if(k0==3)then !333333333>
       status=nf_get_att_int(ncid_,var_id,'valid_min',valmin_)
       status=nf_get_att_int(ncid_,var_id,'valid_max',valmax_)
       ksecu=1
       status=nf_get_vara_int(ncid_,var_id                         &
                             ,varstart(1:ogcm_var_ndims)           &
                             ,varcount(1:ogcm_var_ndims)           &
                              ,short3d(1:max_x,1:max_y,1:max_z))
! Cas bizarre des fichiers simed possedant un double masque:
!     do k=1,max_z ; do j=1,max_y ; do i=1,max_x
!      if(short3d(i,j,k)==0)short3d(i,j,k)=nint(filval)
!     enddo        ; enddo        ; enddo

       do k=1,max_z
       do j=1,max_y
       do i=1,max_x
        if(short3d(i,j,k)==nint(filval)) then                      !23-04-14
         ogcm_var3d(i,j,k)=filval
        else
         ogcm_var3d(i,j,k)=var_scalefactor*short3d(i,j,k)+var_addoffset
        endif
       enddo
       enddo
       enddo
      endif         !333333333>
      if(status/=0)stop 'error nf_get_vara ogcm T or S'
      if(ksecu==0)stop 'interp_ts type variable non identifie'

      if(par%rank==0) then
       if(ichoix==1)write(6,*)'Ogcm temperature is loaded'
       if(ichoix==3)write(6,*)'Ogcm reference temperature is loaded'
       if(ichoix==2)write(6,*)'Ogcm salinity is loaded'
       if(ichoix==4)write(6,*)'Ogcm reference salinity is loaded'
      endif

      if(obc_ogcm_type(1:6)=='nemo_z'.or.obc_ogcm_type(1:4)=='ncom') then !nenenen>
! Changer le sens de l'axe vertical
        call ogcm_renverser_axe_vertical_var
      endif                                 !nenenen>

! Niveaux sous le fond "colles" au niveau le plus profond:               !19-11-10
      if(obc_ogcm_type(1:5)=='sympa')call ogcm_bottom_boundary_condition

! Changer la valeur de masque
      do k=1,max_z
      do j=1,max_y
      do i=1,max_x
!      if(ogcm_var3d(i,j,k)==filval.or.ogcm_var3d(i,j,k)==0.) ogcm_var3d(i,j,k)=1.e11 !23-04-14
       if(ogcm_var3d(i,j,k)==filval) ogcm_var3d(i,j,k)=1.e11 !23-04-14
      enddo
      enddo
      enddo

! Interpolation. Dans le cas initial on fait 2 passages. Le passage 1 sert
! e detecter les trous et e les boucher.
! Principe de l'interpolation: sur l'horizontale on applique un schema
! lineaire dans chaque direction. Sur la verticale, on suppose que le profil
! est lineaire entre 2 valeurs consecutives et que la valeur au point de grille
! de Smodel correspond e la moyenne de l'ogcm (d'oe l'integrale) sur la
! maille de Smodel. Ceci permet de mieux restituer la pression de l'ogcm.

! CommentE le 06-06-18
!     if(ichoix==1.and. &
!        obc_ogcm_type(1:5)=='sympa')call ogcm_grille_sous_le_fond   !17-02-12

      flag=0
      if(ichoix==1) then !ooo>
         if(t_==0.and.iteration3d==kount0)flag=1
         if(flag_updatebouchetrou==1)           flag=1 !17-05-19
      endif              !ooo>
      if(flag==1) then !m°v°m> !17-05-19
           if(obc_ogcm_type(1:5)=='sympa') then  !sssssss>
            call ogcm_grille_sous_le_fond !06-06-18
            call ogcm_repere_les_trous
            call ogcm_fichier_bouchetrou
           endif                                 !sssssss>
           if(obc_ogcm_type(1:6)=='nemo_z'.or.  &
              obc_ogcm_type(1:4)=='ncom')   then !nnnnnnn>
            call ogcm_grille_sous_le_fond
            call ogcm_repere_les_trous_ts
            call ogcm_fichier_bouchetrou_ts
           endif                                 !nnnnnnn>
      endif            !m°v°m> !17-05-19


! remplir les trous:
       if(obc_ogcm_type(1:5)=='sympa') call ogcm_appli_bouchetrou
       if(obc_ogcm_type(1:6)=='nemo_z'                         &
      .or.obc_ogcm_type(1:4)=='ncom') call ogcm_appli_bouchetrou_ts    !19-11-10

      if(par%rank==0)write(6,*)'calcul interpolation'

      if(ichoix==1.and.obc_ogcm_type(1:5)=='sympa') call ogcm_interp_core1(ichoix,t_,nc_)
      if(ichoix==2.and.obc_ogcm_type(1:5)=='sympa') call ogcm_interp_core1(ichoix,t_,nc_)
      if(ichoix==1.and.obc_ogcm_type(1:6)=='nemo_z')call ogcm_interp_core3(ichoix,t_,nc_)
      if(ichoix==2.and.obc_ogcm_type(1:6)=='nemo_z')call ogcm_interp_core3(ichoix,t_,nc_)
      if(ichoix==3)call ogcm_interp_core2(ichoix,t_,nc_)
      if(ichoix==4)call ogcm_interp_core2(ichoix,t_,nc_)

! Fermer le fichier contenant t, s et ssh:
      status=nf_close(ncid_)                                            !04-09-09

! mpi: echanges z0 sur temobc_t(:,:,:,t_) et salobc_t(:,:,:,t_)
      if(ichoix==1)call obc_scal_temobc(t_)
      if(ichoix==2)call obc_scal_salobc(t_)
      if(ichoix==3)call obc_scal_temref('z0')
      if(ichoix==4)call obc_scal_salref('z0')

      end subroutine ogcm_interp_ts

!------------------------------------------------------------------------------

      subroutine ogcm_interp_uv(ichoix,t_   ,nc_   )
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
      include 'netcdf.inc'
      double precision zup_,zdw_,vup_,vdw_,dlon_di_
      integer ichoix,t_   ,nc_
#ifdef synopsis
       subroutinetitle='ogcm_interp_uv'
       subroutinedescription= &
          'Reads U and V in the netcdf OGCM file and interpolates on' &
       //' the S grid and stores in velobc_u and velobc_v.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(obc_ogcm_type(1:6)=='nemo_z'.or.obc_ogcm_type(1:4)=='ncom')then !nenenenenen>

      if(ichoix==1)texte30='depthu'
      if(ichoix==2)texte30='depthv'

      status=nf_inq_varid(ncid1,trim(texte30),var_id)
      if(status/=0)then                 ! 03-06-10
      texte30='gdept'                   ! 03-06-10
      status=nf_inq_varid(ncid1,trim(texte30),var_id)     ! 03-06-10
      endif                                               ! 03-06-10
      if(status/=0)then                 ! 03-06-10
      texte30='depth'                   ! 03-06-10
      status=nf_inq_varid(ncid1,trim(texte30),var_id)     ! 03-06-10
      endif
      if(status/=0)status=nf_inq_varid(ncid1,'deptht',var_id)
      if(status/=0)stop 'interp_uv erreur sur nom variable profondeur'

      status=nf_inq_vartype(ncid1,var_id,k0)     ! 03-06-10
      !print *,"nf_inq_vartype k0=",k0,trim(texte30)
      if((k0/=5).and.(k0/=6))   &
       stop 'interp_uv n''a pas prevu ogcm_z non real'

! si cartesien z est un tableau 1DV:
      status=nf_get_var_real(ncid1,var_id,ogcm_z(1:1,1:1,1:max_z))

! passage conventions nemo-symphonie:
! changer le signe !15-05-19
       ogcm_z(1,1,1:max_z)=-ogcm_z(1,1,1:max_z)
! axe vertical pointant vers le haut !15-05-19
       do k=1,max_z/2
         x1=ogcm_z(1,1,k)
            ogcm_z(1,1,k)=ogcm_z(1,1,max_z+1-k)
                          ogcm_z(1,1,max_z+1-k)=x1
       enddo

! Le tableau de profondeur 1DV est transforme en un tableau 3D:
      call ogcm_transformer_1dv_en_3d

      endif                                !nenenenenen>

      if(obc_ogcm_type(1:5)=='sympa')then  !sssssssssss>

      if(ichoix==1) then !1111>
                   call ogcm_get_z_sympa(nc_,'u')
! Comme la routine precedente remet max_x et max_y aux valeurs des points 't'
! on recalcule max_x max_y des points 'u'
                   call ogcm_zoom_coordinate('u',1)
      endif              !1111>
      if(ichoix==2) then !2222>
                   call ogcm_get_z_sympa(nc_,'v')
! Comme la routine precedente remet max_x et max_y aux valeurs des points 't'
! on recalcule max_x max_y des points 'v'
                   call ogcm_zoom_coordinate('v',1)
      endif              !2222>

      endif                                !sssssssssss>

! On recalle les surfaces des 2 model_s l'une par rapport e l'autre:
!     call top_level_at_zero_meter

! ICHOIX=1 : Interpolation de U
! ICHOIX=2 : Interpolation de V

      if(ichoix==1.and.par%rank==0)write(6,*)'About to read ogcm u'
      if(ichoix==2.and.par%rank==0)write(6,*)'About to read ogcm v'
!     if(ichoix==1)call ogcm_get_dim(ncid1,'u')
!     if(ichoix==2)call ogcm_get_dim(ncid1,'v')

!--- IDENTIFIANTS ->
      if(ichoix==1.and.obc_ogcm_type(1:6)=='nemo_z') then
       texte30='vozocrtx'
       status=nf_inq_varid(ncid1,trim(texte30),var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'u',var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'U',var_id)  !16-05-18
       if(status/=0)status=nf_inq_varid(ncid1,'uo',var_id) !25-09-18
       if(status/=0)stop 'interp_uv erreur sur nom variable courant 1'
      endif
      if(ichoix==2.and.obc_ogcm_type(1:6)=='nemo_z') then
       texte30='vomecrty'
       status=nf_inq_varid(ncid1,trim(texte30),var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'v',var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'V',var_id) !16-05-18
       if(status/=0)status=nf_inq_varid(ncid1,'vo',var_id) !25-09-18
       if(status/=0)stop 'interp_uv erreur sur nom variable courant 2'
      endif
      if(ichoix==1.and.obc_ogcm_type(1:4)=='ncom') then
       texte30='water_u'
       status=nf_inq_varid(ncid1,trim(texte30),var_id)
       if(status/=0)stop 'interp_uv erreur sur nom variable courant 3'
      endif
      if(ichoix==2.and.obc_ogcm_type(1:4)=='ncom') then
       texte30='water_v'
       status=nf_inq_varid(ncid1,trim(texte30),var_id)
       if(status/=0)stop 'interp_uv erreur sur nom variable courant 4'
      endif

      if(ichoix==1.and.obc_ogcm_type(1:5)=='sympa') then
       texte30='u'
       status=nf_inq_varid(ncid1,trim(texte30),var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'U',var_id)
       if(status/=0)stop 'interp_uv erreur sur nom variable courant 5'
      endif
      if(ichoix==2.and.obc_ogcm_type(1:5)=='sympa') then
       texte30='v'
       status=nf_inq_varid(ncid1,trim(texte30),var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'V',var_id)
       if(status/=0)stop 'interp_uv erreur sur nom variable courant 6'
      endif

!--- FACTEUR et BIAIS ->
       status=nf_get_att_real(ncid1,var_id,'scale_factor',var_scalefactor)
!      if(status/=0)stop 'interp_ogcm erreur var_scalefactor ssh'
       if(status/=0)var_scalefactor=1.
       status=nf_get_att_real(ncid1,var_id,'add_offset',var_addoffset)
!      if(status/=0)stop 'interp_ogcm erreur var_addoffset ssh'
       if(status/=0)var_addoffset=0.
       if(par%rank==0)write(6,*)var_scalefactor,var_addoffset

!--- FILL VALUE ->
                   status=nf_get_att_real(ncid1,var_id,'_FillValue',filval)  !23-11-11
      if(status/=0)status=nf_get_att_real(ncid1,var_id,'_Fillvalue',filval)  !15-02-16
      if(status/=0)status=nf_get_att_real(ncid1,var_id,'missing_value',filval)
!     if(status/=0)stop 'interp_ogcm erreur missing_value'
      if(status/=0)filval=0.

!--- LECTURE U ou V ->
      status=nf_inq_varndims(ncid1,var_id,ogcm_var_ndims)
      varstart(1)=ogcmzoom_istr     ; varcount(1)=max_x
      varstart(2)=ogcmzoom_jstr     ; varcount(2)=max_y
      varstart(3)=1                 ; varcount(3)=max_z
      varstart(4)=ogcm_time_counter ; varcount(4)=1

      status=nf_inq_vartype(ncid1,var_id,k0)
      ksecu=0
      if(k0==5) then !........>
       ksecu=1
       status=nf_get_vara_real(ncid1,var_id                         &
                              ,varstart(1:ogcm_var_ndims)           &
                              ,varcount(1:ogcm_var_ndims)           &
                              ,ogcm_var3d(1:max_x,1:max_y,1:max_z))
      endif          !........>
      if(k0==3) then !........>
       ksecu=1
       status=nf_get_vara_int(ncid1,var_id                         &
                             ,varstart(1:ogcm_var_ndims)           &
                             ,varcount(1:ogcm_var_ndims)           &
                             ,short3d(1:max_x,1:max_y,1:max_z))

       do k=1,max_z
       do j=1,max_y
       do i=1,max_x
        if(short3d(i,j,k)==nint(filval)) then
          ogcm_var3d(i,j,k)=filval
        else
          ogcm_var3d(i,j,k)=var_scalefactor*short3d(i,j,k)+var_addoffset
        endif
       enddo
       enddo
       enddo
      endif          !........>
      if(ksecu==0)stop 'interp_uv type non identifie'
      if(status/=0)stop 'error nf_get_vara ogcm u or v'

      if(ichoix==1.and.par%rank==0)write(6,*)'Ogcm u is loaded'
      if(ichoix==2.and.par%rank==0)write(6,*)'Ogcm v is loaded'


      if(obc_ogcm_type(1:6)=='nemo_z'.or.obc_ogcm_type(1:4)=='ncom') then !nenenen>
! Changer le sens de l'axe vertical
        call ogcm_renverser_axe_vertical_var
! Niveaux sous le fond "colles" au niveau le plus profond:
        call ogcm_bottom_boundary_condition
      endif                                 !nenenen>

! Changer la valeur de masque
      do k=1,max_z
      do j=1,max_y
      do i=1,max_x
       if(ogcm_var3d(i,j,k)==filval)ogcm_var3d(i,j,k)=0.              !22-10-09
      enddo
      enddo
      enddo

      loop1=1

! remplir les trous:
      if(loop1.eq.1.and.par%rank==0)write(6,*)'calcul interpolation'

      if(ichoix==1)                                                     &
        open(unit=3,file=trim(tmpdirname)//'hr_to_lr_coord_2_'//dom_c//'.out')
      if(ichoix==2)                                                     &
        open(unit=3,file=trim(tmpdirname)//'hr_to_lr_coord_3_'//dom_c//'.out')

! Interpolation:
      do j=0,jmax+1
      do i=0,imax+1

!     read(3,*)i0,j0,deci,decj,x0,rapi,rapj
      read(3,*)i0,j0,i1,j1,x0,rapi,rapj
      if(i0.ne.i)stop 'erreur1 sur i'
      if(j0.ne.j)stop 'erreur1 sur j'

      if(  mask_t(i,j,kmaxp1).eq.1) then !>>>>>>>>>>>>

!      i1=int(deci)
!      j1=int(decj)

       do k=kmax,kmin_w(i,j),-1

        sum1=0.
        sum2=0.

        do k1=max_z-1,1,-1

! Z1 & Z2  profondeurs dans la grille grossiere en K1 & K1+1
        z1= (1.-rapi)*(1.-rapj)*ogcm_z(i1  ,j1  ,k1)                    &
           +(1.-rapi)*    rapj *ogcm_z(i1  ,j1+1,k1)                    &
           +    rapi *(1.-rapj)*ogcm_z(i1+1,j1  ,k1)                    &
           +    rapi *    rapj *ogcm_z(i1+1,j1+1,k1)
        z2= (1.-rapi)*(1.-rapj)*ogcm_z(i1  ,j1  ,k1+1)                  &
           +(1.-rapi)*    rapj *ogcm_z(i1  ,j1+1,k1+1)                  &
           +    rapi *(1.-rapj)*ogcm_z(i1+1,j1  ,k1+1)                  &
           +    rapi *    rapj *ogcm_z(i1+1,j1+1,k1+1)

         if(k1==max_z-1)z2=max(z2,depth_w(i,j,kmaxp1))          !19-11-11

! Pour s'assurer que les niveaux de surface des 2 codes coincident bien,
! c'est z-ssh que l'on considere et pas z:
!        zup_   =min( depth_w(i,j,k+1)-depth_w(i,j,kmaxp1) ,z2)
!        zdw_   =max( depth_w(i,j,k  )-depth_w(i,j,kmaxp1) ,z1)
         zup_   =min( depth_w(i,j,k+1)                     ,z2) !19-11-11
         zdw_   =max( depth_w(i,j,k  )                     ,z1) !19-11-11

         if(zup_   -zdw_   .gt.zero) then !-_-_-_-_-_-_>

           x2=(1.-rapi)*(1.-rapj)*ogcm_var3d(i1  ,j1  ,k1+1)            &
             +(1.-rapi)*    rapj *ogcm_var3d(i1  ,j1+1,k1+1)            &
             +    rapi *(1.-rapj)*ogcm_var3d(i1+1,j1  ,k1+1)            &
             +    rapi *    rapj *ogcm_var3d(i1+1,j1+1,k1+1)
           x1=(1.-rapi)*(1.-rapj)*ogcm_var3d(i1  ,j1  ,k1  )            &
             +(1.-rapi)*    rapj *ogcm_var3d(i1  ,j1+1,k1  )            &
             +    rapi *(1.-rapj)*ogcm_var3d(i1+1,j1  ,k1  )            &
             +    rapi *    rapj *ogcm_var3d(i1+1,j1+1,k1  )

           rapk=(zdw_   -z1)/(z2-z1)

           vdw_   =(1.-rapk)*x1+rapk*x2

           rapk=(zup_   -z1)/(z2-z1)

           vup_   =(1.-rapk)*x1+rapk*x2

          sum1=sum1+(zup_   -zdw_   )
          sum2=sum2+(zup_   -zdw_   )*0.5*(vup_   +vdw_   )

         endif          !-_-_-_-_-_-_>

        enddo ! fin K1

        anyv3d(i,j,k,ichoix)=sum2/max(sum1,small2)

       enddo   ! fin K

! Dans la couche fusionnee: !16-10-17
        anyv3d(i,j,1:kmin_w(i,j)-1,ichoix)=anyv3d(i,j,kmin_w(i,j),ichoix)

      endif                         !>>>>>>>>>>>>>>
      enddo   ! fin J
      enddo   ! Fin I
      close(3)

      if(ichoix.eq.1) return


! Rotation
      open(unit=3, &
      file=trim(tmpdirname)//'hr_to_lr_coord_1_'//dom_c//'.out')

      do j=0,jmax+1
      do i=0,imax+1

!     read(3,*)i0,j0,deci,decj,x0
      read(3,*)i0,j0,i1,j1,x0
      if(i0.ne.i)stop 'erreur1 sur i'
      if(j0.ne.j)stop 'erreur1 sur j'

! x0 rotation modele externe vers Smodel                             !18-12-09
       dlon_di_=lon_t(i+1,j)-lon_t(i-1,j)
       if(dlon_di_<-pi)dlon_di_=dlon_di_+2.*pi
       if(dlon_di_> pi)dlon_di_=dlon_di_-2.*pi
      x0=x0*deg2rad                          & ! rotation modele externe
        -atan2(  lat_t(i+1,j)-lat_t(i-1,j)   & ! rotation Smodel     !09-01-10
                ,dlon_di_*cos(lat_t(i,j))  )
!              ,(lon_t(i+1,j)-lon_t(i-1,j))*cos(lat_t(i,j))  )

      do k=1,kmaxp1-1
       anyv3d(i,j,k,3)=cos(x0)*anyv3d(i,j,k,1)                          &
                      -sin(x0)*anyv3d(i,j,k,2)
       anyv3d(i,j,k,4)=sin(x0)*anyv3d(i,j,k,1)                          &
                      +cos(x0)*anyv3d(i,j,k,2)
      enddo

      enddo
      enddo

      close(3)                                                          !17-12-09

! Grille C:
      do k=1,kmax
       do j=1,jmax+1
       do i=1,imax+1
        velobc_u(i,j,k,t_   )=0.5*(anyv3d(i-1,j  ,k,3)+anyv3d(i,j,k,3)) &
       *  mask_u(i,j,k)
        velobc_v(i,j,k,t_   )=0.5*(anyv3d(i  ,j-1,k,4)+anyv3d(i,j,k,4)) &
       *  mask_v(i,j,k)
       enddo
       enddo
      enddo

      end subroutine ogcm_interp_uv

!------------------------------------------------------------------------------

      subroutine ogcm_get_bathy_ogcm
      use module_principal ; use module_parallele ; use module_forcages
      use module_global
      implicit none
      double precision zup_,zdw_,vup_,vdw_
      integer :: ncid_,flag_i_,flag_j_  &
                ,flag_bathy_=0          &  !21-04-14
                ,flag_mask_=0              !21-04-14
#ifdef synopsis
       subroutinetitle='ogcm_get_bathy_ogcm'
       subroutinedescription= &
         ' Interpolates the OGCM sea-land mask and bathymetry and'   &
       //' merges with the S mask and bathymetry within a lateral'   &
       //' layer boundary in order to improve the continuity of the' &
       //' S grid with the OGCM grid'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif
 
      if(par%rank==0)write(6,'(a,a,a)')'subroutine ogcm_get_bathy_ogcm' &
      ,' obc_ogcm_type=',trim(obc_ogcm_type)

      if(par%rank==0)write(6,'(a,a)')'open ',trim(mergebathy_filename)
      status=nf_open(trim(mergebathy_filename),nf_nowrite,ncid_)
      if(status/=0)stop ' Erreur nf_open mergebathy_filename'

                   status=nf_inq_dimid(ncid_,'x',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'ni_t',dim_x_id) !23-01-19
      if(status/=0)stop ' Erreur nf_inq_dimid x get_bathy_ogcm'
                   status=nf_inq_dimid(ncid_,'y',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'nj_t',dim_y_id) !23-01-19
      if(status/=0)stop ' Erreur nf_inq_dimid y get_bathy_ogcm'
                   status=nf_inq_dimid(ncid_,'z',dim_z_id)
      if(status/=0)status=nf_inq_dimid(ncid_,'nk_t',dim_z_id) !23-01-19
      if(status/=0)status=nf_inq_dimid(ncid_,'deptht',dim_z_id) !21-04-14
      if(status/=0)stop ' Erreur nf_inq_dimid z get_bathy_ogcm'

      status=nf_inq_dimlen(ncid_,dim_x_id,max_x)
      if(status/=0)stop ' Erreur nf_inq_dimlen max_x'
      status=nf_inq_dimlen(ncid_,dim_y_id,max_y)
      if(status/=0)stop ' Erreur nf_inq_dimlen max_y'
      status=nf_inq_dimlen(ncid_,dim_z_id,max_z)
      if(status/=0)stop ' Erreur nf_inq_dimlen max_z'

      call ogcm_allocate(1,1,max_x,max_y,max_z) ! arg1=allouer arg2=ogcm arg3,4,5=dimensions
      ogcmzoom_istr=1 ; ogcmzoom_jstr=1
      call ogcm_get_lonlat(ncid_,'t')
      call ogcm_hr_to_lr_coord(1)

      flag=0
      if(obc_ogcm_type(1:5)=='sympa') flag=1
      if(obc_ogcm_type(1:6)=='nemo_z')flag=2
      if(flag==0)stop 'Undefined obc_ogcm_type in ogcm_get_bathy_ogcm'

!ssssssssssssss>
      if(obc_ogcm_type(1:5)=='sympa')  then !sympa-sympa-sympa>   !23-01-19

        flag_mask_=0 !le 23-01-19 on commence par ignorer l'aspect land/sea mask. Seule la bathy est considErEe
! LIRE BATHY OGCM sympa
       varstart(1)=ogcmzoom_istr     ; varcount(1)=max_x
       varstart(2)=ogcmzoom_jstr     ; varcount(2)=max_y
       allocate(ogcmvar_r8(max_x,max_y,1))
       status=nf_inq_varid(ncid_,'h_w',var_id) 
       if(status/=0)stop 'Err 1892 nf_inq_varid'
       status=nf_get_vara_double(ncid_,var_id,varstart(1:2),varcount(1:2),ogcmvar_r8(1:max_x,1:max_y,1))
       if(status/=0)stop 'Err 1895 nf_get_vara_double'
       open(unit=3, &
       file=trim(tmpdirname)//'hr_to_lr_coord_1_'//dom_c//'.out')
        do j=0,jmax+1 ; do i=0,imax+1
        read(3,*)i0,j0,i1,j1,x0,rapi,rapj
        if(i0.ne.i)stop 'erreur1 sur i'
        if(j0.ne.j)stop 'erreur1 sur j'
! Interpolation bathy:
        xy_t(i,j,1)=                                               &
               (1.-rapi)*(1.-rapj)*ogcmvar_r8(i1  ,j1  ,1)         &
              +(1.-rapi)*    rapj *ogcmvar_r8(i1  ,j1+1,1)         &
              +    rapi *(1.-rapj)*ogcmvar_r8(i1+1,j1  ,1)         &
              +    rapi *    rapj *ogcmvar_r8(i1+1,j1+1,1)

        enddo ; enddo   
        flag_bathy_=1
       close(3)
       deallocate(ogcmvar_r8)

      endif                                 !sympa-sympa-sympa>   !23-01-19
!ssssssssssssss>

!nnnnnnnnnnnnnn>
      if(obc_ogcm_type(1:6)=='nemo_z') then !nemo-nemo-nemo>   !23-01-19

! LIRE LAND-SEA MASK (tmask)
      status=nf_inq_varid(ncid_,'tmask',var_id)
      if(status==0) then !>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !29-01-15
!     if(status/=0)stop ' erreur nf_inq_varid tmaskutil'

      status=nf_inq_varndims(ncid_,var_id,ogcm_var_ndims)
      if(status/=0)stop ' erreur 1639 nf_inq_varndims tmask'
      if(ogcm_var_ndims<3) &
      stop ' cas ogcm_var_ndims<3 pas prevu pour tmask'

      varstart(1)=ogcmzoom_istr     ; varcount(1)=max_x
      varstart(2)=ogcmzoom_jstr     ; varcount(2)=max_y
      varstart(3)=1                 ; varcount(3)=max_z
      varstart(4)=1                 ; varcount(4)=1

      status=nf_inq_vartype(ncid_,var_id,k0)
      if(status/=0)stop ' erreur 1647 nf_inq_vartype tmask'
      if(k0/=nf_byte)stop ' erreur 1648 type tmaskutil non prevu'

      allocate(ogcmvar_int(max_x,max_y,varcount(3)))
      status=nf_get_vara_int(ncid_,var_id                       &
                          ,varstart(1:ogcm_var_ndims)           &
                          ,varcount(1:ogcm_var_ndims)           &
                          ,ogcmvar_int(1:max_x,1:max_y,1:varcount(3)))
      if(status/=0)stop ' erreur 1655 nf_get_vara_int ogcmvar_int'

        flag_mask_=1

      else               !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        flag_mask_=0

      endif              !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


! LIRE BATHY ou E3T:
      varstart(1)=ogcmzoom_istr     ; varcount(1)=max_x
      varstart(2)=ogcmzoom_jstr     ; varcount(2)=max_y
      varstart(3)=1 !  varcount(3) depend du nom de la vaiable a suivre...
      varstart(4)=1                 ; varcount(4)=1
                   status=nf_inq_varid(ncid_,'hdepw',var_id)
      if(status==0)varcount(3)=1
      if(status/=0)then !>>>>>!25-11-14
             status=nf_inq_varid(ncid_,'e3t',var_id)
             if(status==0)varcount(3)=max_z
             if(status/=0)stop ' error nf_inq_varid e3t'
      endif             !>>>>>!25-11-14

      status=nf_inq_varndims(ncid_,var_id,ogcm_var_ndims)
      if(status/=0)stop ' erreur nf_inq_varndims'

      status=nf_inq_vartype(ncid_,var_id,k0)
      if(status/=0)stop ' erreur nf_inq_vartype k0'
      if(k0/=6)stop ' Routine get_bathy double precision svp!'

      if(k0==6) then !dpdpdpdpdpdp>
          allocate(ogcmvar_r8(max_x,max_y,varcount(3)))
       status=nf_get_vara_double(ncid_,var_id                          &
                              ,varstart(1:ogcm_var_ndims)              &
                              ,varcount(1:ogcm_var_ndims)              &
                            ,ogcmvar_r8(1:max_x,1:max_y,1:varcount(3)))
        if(status/=0)stop ' erreur nf_get_vara_double ogcmvar_r8'

! si varcount(3)=max_z on a affaire a e3t. La bathy se deduit de son integrale verticale:
       if(varcount(3)==max_z) then !iii>
        do j=1,max_y ; do i=1,max_y
         ogcmvar_r8(i,j,1)=ogcmvar_r8(i,j,1)*ogcmvar_int(i,j,1)
        enddo        ; enddo
        do k=2,max_z
         do j=1,max_y ; do i=1,max_x
          ogcmvar_r8(i,j,1)=ogcmvar_r8(i,j,1) &
                           +ogcmvar_r8(i,j,k)*ogcmvar_int(i,j,k)
         enddo        ; enddo
        enddo
       endif                       !iii>

       open(unit=3, &
       file=trim(tmpdirname)//'hr_to_lr_coord_1_'//dom_c//'.out')

        do j=0,jmax+1
        do i=0,imax+1
        read(3,*)i0,j0,i1,j1,x0,rapi,rapj
        if(i0.ne.i)stop 'erreur1 sur i'
        if(j0.ne.j)stop 'erreur1 sur j'

! Interpolation bathy:
        xy_t(i,j,1)=                                               &
               (1.-rapi)*(1.-rapj)*ogcmvar_r8(i1  ,j1  ,1)         &
              +(1.-rapi)*    rapj *ogcmvar_r8(i1  ,j1+1,1)         &
              +    rapi *(1.-rapj)*ogcmvar_r8(i1+1,j1  ,1)         &
              +    rapi *    rapj *ogcmvar_r8(i1+1,j1+1,1)

! Interpolation mask:
        if(flag_mask_==1) then !ooo>
        xy_t(i,j,2)=                                              &
               (1.-rapi)*(1.-rapj)*ogcmvar_int(i1  ,j1  ,1)       &
              +(1.-rapi)*    rapj *ogcmvar_int(i1  ,j1+1,1)       &
              +    rapi *(1.-rapj)*ogcmvar_int(i1+1,j1  ,1)       &
              +    rapi *    rapj *ogcmvar_int(i1+1,j1+1,1)
        endif                  !ooo>

        enddo   ! fin J
        enddo   ! Fin I
        flag_bathy_=1

       close(3)

      if(allocated(ogcmvar_r8)) deallocate(ogcmvar_r8)
      if(allocated(ogcmvar_int))deallocate(ogcmvar_int) !29-01-15
      endif          !dpdpdpdpdpdp>

      endif                                 !nemo-nemo-nemo>
!nnnnnnnnnnnnnn>

! Reprise du tronc commun sympa & nemo_z

      status=nf_close(ncid_)
      call ogcm_allocate(2,1,0,0,0)             ! arg1=desallouer arg2=ogcm

      flag_i_=1 ; flag_j_=1               !05-02-13
      if(jperiodicboundary)flag_j_=0
      if(iperiodicboundary)flag_i_=0


      do j=1,jmax
      do i=1,imax

        dist=1.e10

        if(flag_i_==1) then !>>>>>>>>>>>>>
         do j1=1,jglb                                             !10-06-09

         i1=1
           if(glob_mask  (i1,j1)==1)                       & !01-02-15
           dist=min(dist                                   &
           ,sqrt(flag_i_*real(i+par%timax(1)-i1)**2+       &
                 flag_j_*real(j+par%tjmax(1)-j1)**2))
!          /max(un*glob_mask  (i1,j1),small1))

         i1=iglb
           if(glob_mask  (i1,j1)==1)                       & !01-02-15
           dist=min(dist                                   &
           ,sqrt(flag_i_*real(i+par%timax(1)-i1)**2+       &
                 flag_j_*real(j+par%tjmax(1)-j1)**2))
!          /max(un*glob_mask  (i1,j1),small1))


         enddo ! J1
        endif               !>>>>>>>>>>>>>

        if(flag_j_==1) then !>>>>>>>>>>>>>
         do i1=1,iglb                                              !10-06-09

         j1=1
           if(glob_mask  (i1,j1)==1)                     & !01-02-15
           dist=min(dist                                 &
           ,sqrt(        real(i+par%timax(1)-i1)**2+     &
                         real(j+par%tjmax(1)-j1)**2))
!          /max(un*glob_mask  (i1,j1),small1))

         j1=jglb                                                   !04-09-09
           if(glob_mask  (i1,j1)==1)                     & !01-02-15
           dist=min(dist                                 &
           ,sqrt(        real(i+par%timax(1)-i1)**2+     &
                         real(j+par%tjmax(1)-j1)**2))
!          /max(un*glob_mask  (i1,j1),small1))

         enddo ! I1
        endif               !>>>>>>>>>>>>>

       rap=min(max((real(mergebathy_sponge)-dist)/real(mergebathy_sponge),zero),un)
       x1=rap*flag_bathy_
       x2=rap*flag_mask_

       h_w(i,j)=(1.-x1)*h_w(i,j)+x1*xy_t(i,j,1)
       mask_t(i,j,kmaxp1)=nint( (1.-x2)*mask_t(i,j,kmaxp1)+x2*xy_t(i,j,2) ) !21-04-14

! Decommenter lignes suivantes pour tester (interpolation sur toute la grille)
! que l'interpolation des champs ogcm marche bien....
!      if(i<imax/2)h_w(i,j)=6666.
!      h_w(i,j)=xy_t(i,j,1)
!      mask_t(i,j,kmaxp1)=nint( xy_t(i,j,2) )

      enddo ! i
      enddo ! j

      if(flag_1dv==1) then !pmxpmx>
           h_w(:,:)=xy_t(imax/2,jmax/2,1) !20-06-16
      endif                !pmxpmx>

! OBC:
      call obc_h(0) !01-02-15

! Masque 3D:
      lb3=lbound(mask_t)
      do k=lb3(3),kmax
       mask_t(:,:,k)=mask_t(:,:,kmax+1) !29-12-14
      enddo

!     write(texte30,'(a,i0)')'tmp/h-hnemo-',par%rank
      write(texte30,'(a,i0)')'tmp/hnemo-',par%rank
      open(unit=3,file=trim(texte30)) !18-11-15
       do j=1,jmax ; do i=1,imax
!       write(3,*)h_w(i,j)-xy_t(i,j,1)
        write(3,*)         xy_t(i,j,1) !05-01-17
       enddo         ; enddo
      close(3)

      end subroutine ogcm_get_bathy_ogcm

!------------------------------------------------------------------------------
#ifdef bidon
      subroutine ogcm_visu_matlab_scalar_horiz(k_i1,k_i2,k_j1,k_j2)
      use module_principal
      implicit none
      integer k_i1,k_i2,k_j1,k_j2
#ifdef synopsis
       subroutinetitle='ogcm_visu_matlab_scalar_horiz'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      open(unit=16,file='plan_horiz.dat')
         do i=k_i1,k_i2
         do j=k_j1,k_j2
!         WRITE(16,*)STREAMF_R(I,J,0)
!         WRITE(16,*)XY_Z(I,J,1)
          write(16,*)xy_t(i,j,2)
         enddo
         enddo
      close(16)

      open(unit=16,file='plan_horiz.m')
      write(16,*)'fid=fopen(''plan_horiz.dat'');'
      write(16,*)                                                       &
       'a=fscanf(fid,''%g %g'',[',k_j2-k_j1+1,k_i2-k_i1+1,']);'
      write(16,*)'pcolor(a)'
      write(16,*)'shading interp;'
!     hold on
!     WRITE(16,*)"contourf(a,20,'LineStyle','none')"

      close(16)

      return
      end subroutine ogcm_visu_matlab_scalar_horiz

!------------------------------------------------------------------------------

      subroutine ogcm_visu_matlab_vector_horiz(k_i1,k_i2,k_j1,k_j2)
      use module_principal
      implicit none
      integer k_i1,k_i2,k_j1,k_j2
#ifdef synopsis
       subroutinetitle='ogcm_visu_matlab_vector_horiz'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


      open(unit=16,file='plan_vector_u.dat')
      open(unit=17,file='plan_vector_v.dat')
         do i=k_i1,k_i2
         do j=k_j1,k_j2
          write(16,*)0.5*(xy_u(i,j,1)+xy_u(i+1,j  ,1))
          write(17,*)0.5*(xy_v(i,j,1)+xy_v(i  ,j+1,1))
         enddo
         enddo
      close(16)
      close(17)

      open(unit=16,file='plan_vec_hor.m')
      write(16,*)'fid=fopen(''plan_vector_u.dat'');'
      write(16,*)                                                       &
       'u=fscanf(fid,''%g %g'',[',k_j2-k_j1+1,k_i2-k_i1+1,']);'
      write(16,*)'fid=fopen(''plan_vector_v.dat'');'
      write(16,*)                                                       &
       'v=fscanf(fid,''%g %g'',[',k_j2-k_j1+1,k_i2-k_i1+1,']);'
      write(16,*)'quiver(u,v,2)'
      close(16)


      return
      end subroutine ogcm_visu_matlab_vector_horiz
#endif
!------------------------------------------------------------------------------

      subroutine ogcm_get_vmeanobc(t_   )
      use module_principal
      implicit none
      integer t_
#ifdef synopsis
       subroutinetitle='ogcm_get_vmeanobc'
       subroutinedescription= &
          'Computes the depth-average of the OGCM 3D current and' &
       //' the result in velbarobc_u and velbarobc_v'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! velbar_u
      k=1
      do j=1,jmax
      do i=1,imax+1
       velbarobc_u(i,j,t_)=velobc_u(i,j,k,t_)*dz_u(i,j,k,1)
      enddo
      enddo
      do k=2,kmax
      do j=1,jmax
      do i=1,imax+1
       velbarobc_u(i,j,t_)=    &
       velbarobc_u(i,j,t_)+velobc_u(i,j,k,t_)*dz_u(i,j,k,1)
      enddo
      enddo
      enddo
      do j=1,jmax
      do i=1,imax+1
       velbarobc_u(i,j,t_)=velbarobc_u(i,j,t_)/hz_u(i,j,1)
      enddo
      enddo

! velbar_v
      k=1
      do j=1,jmax+1
      do i=1,imax
       velbarobc_v(i,j,t_)=velobc_v(i,j,k,t_)*dz_v(i,j,k,1)
      enddo
      enddo
      do k=2,kmax
      do j=1,jmax+1
      do i=1,imax
       velbarobc_v(i,j,t_)=    &
       velbarobc_v(i,j,t_)+velobc_v(i,j,k,t_)*dz_v(i,j,k,1)
      enddo
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
       velbarobc_v(i,j,t_)=velbarobc_v(i,j,t_)/hz_v(i,j,1)
      enddo
      enddo

      end subroutine ogcm_get_vmeanobc

!------------------------------------------------------------------------------

      subroutine ogcm_lire_les_listes
      use module_principal
      use module_parallele !#mpi
      implicit none
      include 'netcdf.inc'
      double precision :: tm1_=0.,tm2_=0.
      integer ncid_,loop_,year_,month_,day_,hour_,minute_,second_ &
      ,linemax_,loop2_,loop1_,var_count_,unit_
      integer,dimension(1) :: nelt
      integer :: ncmaxprev_ !16-05-19

#ifdef synopsis
       subroutinetitle='ogcm_lire_les_listes'
       subroutinedescription= &
         'Reads the file lists given in notebook_obcforcing and'&
      //' and creates the binary file lists'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      var_count_=1 ! reset defaut pour tous modeles

! A l'etat initial on lit les listes ascii et on fabrique les listes binaires
! acces record:

      if(par%rank==0)write(6,*)'Routine ogcm_lire_les_listes'

      flag_stop=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cas ogcm = nemo_z
      if(obc_ogcm_type(1:6)=='nemo_z') then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ogcm_rec_max=0

      liste_ogcm_t     =trim(tmpdirname)//''//'binrec'//'_liste_ogcm_t'
      liste_ogcm_s     =trim(tmpdirname)//''//'binrec'//'_liste_ogcm_s'
      liste_ogcm_u     =trim(tmpdirname)//''//'binrec'//'_liste_ogcm_u'
      liste_ogcm_v     =trim(tmpdirname)//''//'binrec'//'_liste_ogcm_v'
      liste_ogcm_ssh   =trim(tmpdirname)//''//'binrec'//'_liste_ogcm_ssh'
      liste_ogcm_grid_t=trim(tmpdirname)//''//'binrec'//'_liste_ogcm_grid_t' !21-06-12
      liste_ogcm_grid_u=trim(tmpdirname)//''//'binrec'//'_liste_ogcm_grid_u' !21-06-12
      liste_ogcm_grid_v=trim(tmpdirname)//''//'binrec'//'_liste_ogcm_grid_v' !21-06-12

      do k1=1,8 !10 ! debut de boucle K1

        if(k1<=2)          var_count_=trc_id
        if(k1>=3.and.k1<=4)var_count_=vel_id
        if(k1==5)          var_count_=ssh_id
        if(k1>=6)          var_count_=0

        if(k1==1) texte30=trim(liste_ogcm_t)
        if(k1==2) texte30=trim(liste_ogcm_s)
        if(k1==3) texte30=trim(liste_ogcm_u)
        if(k1==4) texte30=trim(liste_ogcm_v)
        if(k1==5) texte30=trim(liste_ogcm_ssh)
        if(k1==6) texte30=trim(liste_ogcm_grid_t)
        if(k1==7) texte30=trim(liste_ogcm_grid_u)
        if(k1==8) texte30=trim(liste_ogcm_grid_v)

        if(par%rank==0) then !--rank0 checks existing list-->
         unit_=s_unit(7)
         open(unit=unit_,file=texte30        &
                        ,access='direct'     &
                        ,recl=541            & !28-04-16
                        ,form='unformatted'  &
                        ,status='new'        &
                        ,iostat=k0)
        close(unit_)
       endif                !--rank0 checks existing list-->
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
        call ogcm_oldlist(var_count_) 
        goto 2453
       endif          !existing old list in tmp directory >

       if(par%rank==0)write(6,'(a,a)')'Faire la liste ',trim(texte30)



        if(k1<6) then !----->
         k10=ogcm_rec_prev(var_count_) ! sert e verifier que ogcm_rec_prev(var_count_) ne depend pas de la variable...
         ogcm_rec_prev(var_count_)=-999
        endif         !----->

        open(unit=3,file=obcfile(k1,1))

! Combien y a t'il de lignes dans le fichier?
         linemax_=0
 1834    read(3,'(a)',end=1832)texte250 ; linemax_=linemax_+1 ; goto 1834
 1832    rewind(3)
         if(linemax_==0) then         !>>>>> debug !23-04-14
          write(6,'(a,a,a)')'File ',trim(obcfile(k1,1)),' is empty'
          call ogcm_stop('ogcm_lire_les_listes erreur 1928',1)
         endif                        !>>>>>
!        if(k1>5)linemax_=1  !13-02-16 ! commentE le !16-05-19 pour possibilite plusieurs grilles differentes en cours de simulation

! Chaque proc va traiter une fraction du fichier. On commence par derouler les lignes
! qui ne concernent pas le proc
         do loop2_=1,int(real(par%rank  )/real(nbdom)*linemax_)
          read(3,*) ! lire ces lignes pour rien
         enddo

!        write(6,*)'DECOUPAGE ',par%rank   &
!                              ,int(real(par%rank  )/real(nbdom)*linemax_)+1    &
!                              ,int(real(par%rank+1)/real(nbdom)*linemax_)


        write(texte60,'(a,i0)')trim(tmpdirname)//'tmpfile',par%rank
         open(unit=4,file=texte60,status='REPLACE') !20-10-14

         do loop2_=int(real(par%rank  )/real(nbdom)*linemax_)+1    &
                  ,int(real(par%rank+1)/real(nbdom)*linemax_)

         read(3,'(a)',end=332)texte250
!Si le nom se termine par "blabla.nc , i" avec i un numero "code" de grille (si plusieurs grilles successives) alors tronquer le nom etc..!16-05-19
         ogcmgridcode=1 ! valeur par defaut si pas de code numero de grille acollE au nom du fichier
         if(index(texte250,',')/=0) then !>>>
           j1=index(texte250,',') ; j2=len(trim(texte250))
           texte250=trim(texte250(1:j1-1))
           write(txtformat,'(a,i0,a,i0,a)')'(',j1,'x,i',j2-j1,')'
           backspace 3
           read(3,trim(txtformat))ogcmgridcode
         endif                           !>>>

! Open the netcdf ogcm file and get the time_counter dimension then store in the binary file as much as time_counter:
! so that each line of the binary file corresponds to the ogcm output periodicity:
            max_time_counter=1
            status=nf_open(trim(texte250),nf_nowrite,ncid_)
            if(status/=0) then !ccccc>
             write(6,'(a,a)')'Echec lecture ',trim(texte250) !03-04-14
             call ogcm_stop('nf_open in subroutine ogcm_lire_les_listes',1)
            endif              !ccccc>

! Fichiers de variables seulement (fichier de grille non concernes)
         if(k1<6) then !>>>>>>>

                         status=nf_inq_dimid(ncid_,'time_counter',dim_t_id)
            if(status/=0)status=nf_inq_dimid(ncid_,'TIME_COUNTER',dim_t_id)
            if(status/=0)status=nf_inq_dimid(ncid_,'time',dim_t_id)
            if(status/=0)status=nf_inq_dimid(ncid_,'TIME',dim_t_id)
            if(status/=0)call ogcm_stop('2234 nf_inq_dimid time',1)


            if(status==0)status=nf_inq_dimlen(ncid_,dim_t_id,max_time_counter)
            if(status/=0)call ogcm_stop('Err 2238 nf_inq_dimlen',1)

                         status=nf_inq_varid(ncid_,'time',var_id)
            if(status/=0)status=nf_inq_varid(ncid_,'time_counter',var_id)
            if(status/=0)call ogcm_stop('erreur nf_inq_varid time',1)

            txt_units=''
            status=nf_get_att_text(ncid_,var_id,'units',txt_units)
            if(status/=0)call ogcm_stop('nf_get_att_text units ogcm',1)

            k=index(txt_units,'since')
            read(txt_units(k+6:k+9),*)year_
            read(txt_units(k+11:k+12),*)month_
            read(txt_units(k+14:k+15),*)day_
            read(txt_units(k+17:k+18),*)hour_
            read(txt_units(k+20:k+21),*)minute_
            read(txt_units(k+23:k+24),*)second_
            call datetokount(year_,month_,day_,hour_,minute_,second_) ! donne elapsedtime_out

            status=nf_inq_vartype(ncid_,var_id,var_nftype) ! deplace avant endif le 30-06-15
            if(status/=0)call ogcm_stop('2258 nf_inq_vartype',1)

         endif         !>>>>>>>

            do loop_=1,max_time_counter

! Fichiers de variables seulement (fichier de grille non concernes)
            if(k1<6) then !>>>>>>>

            if(var_nftype==nf_double) then !----->
             status=nf_get_vara_double(ncid_,var_id,loop_,1,x1)
             if(status/=0)call ogcm_stop('2267 nf_get_vara_double',1)
            else                           !----->
              if(var_nftype==nf_real) then !......>
               status=nf_get_vara_real(ncid_,var_id,loop_,1,x1_r4)
               if(status/=0)call ogcm_stop('2271',1)
               x1=x1_r4
              else                         !......>

               if(var_nftype==nf_int.or.  &
                  var_nftype==nf_short) then   !gggggg> !19-05-15
                 status=nf_get_vara_int(ncid_,var_id,loop_,1,k10)
                 if(status/=0)call ogcm_stop('2277',1)
                 x1=k10
               else                            !gggggg>
                   write(6,'(a,a)')'Pb var_nftype fichier ',trim(texte250)
                   write(6,*)'var_nftype=',var_nftype
                   write(6,*)'nf_double= ',nf_double
                   write(6,*)'nf_real=   ',nf_real
                   write(6,*)'nf_int=    ',nf_int
                   write(6,*)'nf_short   ',nf_short
                   call ogcm_stop('2414',1)
               endif                           !gggggg>

              endif                        !......>
            endif                          !----->

             if(status/=0)call ogcm_stop('2420',1)
             x2=-999.
             if(index(txt_units,'days')/=0)x2=86400.
             if(index(txt_units,'hours')/=0)x2=3600.
             if(index(txt_units,'seconds')/=0)x2=1.
             if(x2==-999.)call ogcm_stop('2425',1)
             ogcm_readtime_next(var_count_)=elapsedtime_out+x1*x2

            endif         !>>>>>>>

             write(4,'(a)')trim(texte250)
             write(4,*)loop_,ogcm_readtime_next(var_count_),ogcmgridcode !16-05-19

            enddo ! fin de boucle sur loop_

            status=nf_close(ncid_)
         enddo ! fin de boucle sur loop2_

  332   close(4)
        close(3)



! C'est le proc zero qui a la mission d'assembler tous les fichiers tmp en un seul fichier liste binaire
! La barriere suivante permet de s'assurer que tous les fichiers individuels ont bien ete fait au moment
! d'attaquer la concatenation
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)  !23-06-13
#endif
      if(par%rank==0) then !00000000000>

        if(k1==1) texte30=trim(liste_ogcm_t)
        if(k1==2) texte30=trim(liste_ogcm_s)
        if(k1==3) texte30=trim(liste_ogcm_u)
        if(k1==4) texte30=trim(liste_ogcm_v)
        if(k1==5) texte30=trim(liste_ogcm_ssh)
        if(k1==6) texte30=trim(liste_ogcm_grid_t)
        if(k1==7) texte30=trim(liste_ogcm_grid_u)
        if(k1==8) texte30=trim(liste_ogcm_grid_v)
        unit_=s_unit(7)
        open(unit=unit_,file=texte30,access='direct',recl=541 & !270   & !28-04-16 !16-05-19
                                                    ,form='unformatted')

        nc=1
        do loop1_=0,nbdom-1
         write(texte60,'(a,i0)')trim(tmpdirname)//'tmpfile',loop1_
         open(unit=4,file=texte60)
 2016          tm1_=tm2_
               read(4,'(a)',end=2002)texte250
               read(4,*,end=2002)loop_,tm2_,ogcmgridcode !16-05-19
               if(nc==1)tm1_=tm2_ !03-03-15

               if(k1<6) then !>>>>

! Detection de listes desordonnees: !27-02-15
               if(tm2_<=tm1_.and.nc/=1) then  !dbdbdbdb> !03-03-15!01-12-15
                write(10+par%rank,'(a,a,a)')'File ',trim(texte250) &
                ,' does not respect the chronological order'
                write(10+par%rank,'(a,a)')'List:',trim(texte30)
                write(10+par%rank,*)'tm1_,tm2_',tm1_,tm2_
                write(10+par%rank,*)'nc',nc
                write(10+par%rank,'(a,a)')'File list',trim(texte60)
                call ogcm_stop('2478',1)
                goto 2453
               endif                          !dbdbdbdb>

! Cas la variable time correspond au moment de l'ecriture du champs
               if(ogcm_time_shift==1) then
                  ogcm_readtime_next(var_count_)=0.5*(tm1_+tm2_)
                    ogcm_period_next(var_count_)=     tm2_-tm1_
                    if(nc==2)ogcm_period_prev(var_count_) &
                            =ogcm_period_next(var_count_)
               endif

! Cas la variable time est centree sur le baricentre de l'echantillonage
               if(ogcm_time_shift==0) then
                 ogcm_readtime_next(var_count_)=tm2_
                  if(nc==2) then ! cas particulier >
                     ogcm_period_next(var_count_)=tm2_-tm1_
                     ogcm_period_prev(var_count_)=tm2_-tm1_
                  else           ! cas general >
! Attention ogcm_period_prev(var_count_) ne fait pas la bascule il faut
! donc faire ogcm_period_next(var_count_)=2.*(tm2_-tm1_)-ogcm_period_next(var_count_)
! et non pas ogcm_period_next(var_count_)=2.*(tm2_-tm1_)-ogcm_period_prev(var_count_)
                     ogcm_period_next(var_count_)=2.*(tm2_-tm1_)-ogcm_period_next(var_count_)!17-11-18
                  endif
               endif

                write(unit_,rec=nc)texte250,loop_             &
                              ,ogcm_readtime_next(var_count_) &
                              ,ogcm_period_next(var_count_)   &
                              ,ogcmgridcode !16-05-19
                ogcm_rec_max=max(ogcm_rec_max,nc)

                if(nc/=1.and.ogcm_period_next(var_count_)<=0.) then !28-11-15 ! debuger
                 write(10+par%rank,'(a,a)')'tmpfile=',trim(texte60)
                 write(10+par%rank,'(a,a)')'file=',trim(texte250)
                 write(10+par%rank,*)'loop_,tm2_',loop_,tm2_
                 write(10+par%rank,*)'nc',nc
                 write(10+par%rank,*)'var_count_',var_count_
                 write(10+par%rank,*)'ogcm_period_next',ogcm_period_next(var_count_)
                 write(10+par%rank,*)'ogcm_period_prev',ogcm_period_prev(var_count_)
                 write(10+par%rank,*)'tm1_ tm2_',tm1_,tm2_
                 write(10+par%rank,*)'tm2_-tm1_',tm2_-tm1_
                 write(10+par%rank,*)'ogcm_time_shift=',ogcm_time_shift
                 write(10+par%rank,*)'Inconsistent Ogcm file list'
                 if(ogcm_time_shift==0) then !ooo>
                 write(10+par%rank,*)'One field is possibly missing'
                 write(10+par%rank,*)'missing. In that case repair the'
                 write(10+par%rank,*)'defective file, or (but it is not'
                 write(10+par%rank,*)'rigourous) set ogcm_time_shift=1'
                 write(10+par%rank,*)'in notebook_obcforcing'
                 endif                       !ooo>
                 call ogcm_stop('ogcm_period_next(var_count_)<=0',1) !17-11-18
                endif !28-11-15 ! debuger
      


               else          !>>>>
                write(unit_,rec=nc)texte250,loop_,-999.,-999.,ogcmgridcode !16-05-19
!               goto 2002 ! commentE le 16-05-19 pour possibilitE plusieurs grilles successives
               endif         !>>>>

! TRES IMPORTANT: ogcm_readtime_next est calcule pour tous les nc (du debut a la fin)
! et ecrit dans le fichier binaire alors que ogcm_readtime_prev n'est calcule qu'a
! condition de preceder elapsedtime_now, de sorte que la derniere valeur sera celle
! precedent le temps present, autrement dit l'echeance "0" pour les tableaux OBC

               if(nc==2) then !nc=2 nc=2 nc=2>
                if(ogcm_time_shift==1) then
                 if(1.5*tm1_-0.5*tm2_<=elapsedtime_now) then
                  ogcm_rec_prev(var_count_)=nc-1
                  ogcm_readtime_prev(var_count_)=1.5*tm1_-0.5*tm2_
                 endif
                endif
               endif          !nc=2 nc=2 nc=2>

               if(ogcm_readtime_next(var_count_)<=elapsedtime_now) then !----->
                     ogcm_rec_prev(var_count_)=nc
                ogcm_readtime_prev(var_count_)=ogcm_readtime_next(var_count_)
                  ogcm_period_prev(var_count_)=  ogcm_period_next(var_count_)
               endif                                                    !----->


! Refaire proprement nc=1 en se servant de nc=2
               if(nc==2) then !nc=2 nc=2 nc=2>
                 read(unit_,rec=1)texte250,i0,x0,x2,ogcmgridcode !16-05-19
                 if(ogcm_time_shift==1)x0=1.5*tm1_-0.5*tm2_
                write(unit_,rec=1,iostat=i10)texte250,i0,x0      &
                                 ,ogcm_period_prev(var_count_)   &
                                 ,datesim(1:6,1) &! date repere nc=1
                                 ,ogcmgridcode !16-05-19
                   if(i10/=0)call ogcm_stop('2553',1)
               endif          !nc=2 nc=2 nc=2>
               nc=nc+1
               goto 2016
 2002    close(4)
         enddo ! loop1_

! Finalement on utilise rec=1 pour archiver date "repere" pour elapsedtime_now=0
!       write(unit=unit_,rec=1)datesim(1:6,1)

        close(unit_)


      endif                !00000000000>

2453 continue ! point de RDV si listes binaires dEjA existante
     call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
     if(k0/=0)stop 'Err 2453 ogcm file list'
     

#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)  !23-06-13
#endif

! A ce stade nc-1 est la valeur maximale de record de la liste constituee
! on memorise cette valeur dans ncmaxprev_ et on verifie que les listes des
! variables ont bien toute la meme longueur:
      if(k1>1.and.k1<6) then
       if(nc-1/=ncmaxprev_) then !>>>>
        write(6,'(a,a)')trim(texte30) &
        ,' length is different from the others' 
        write(6,*)'line number of these 2 lists:',ncmaxprev_,nc-1
        call ogcm_stop('2586',1)
       endif                     !>>>>
      endif
      ncmaxprev_=nc-1

      enddo ! fin de boucle K1

#ifdef parallele
      nelt = ubound(ogcm_readtime_next)    &
            -lbound(ogcm_readtime_next)+1
      call mpi_bcast(ogcm_readtime_next    &
             ,nelt(1)                      &
             ,mpi_double_precision         &
             ,0,par%comm2d,ierr)

      nelt = ubound(ogcm_readtime_prev)   &
            -lbound(ogcm_readtime_prev)+1
      call mpi_bcast(ogcm_readtime_prev    &
             ,nelt(1)                      &
             ,mpi_double_precision         &
             ,0,par%comm2d,ierr)

      nelt = ubound(ogcm_period_next)   &
            -lbound(ogcm_period_next)+1
      call mpi_bcast(ogcm_period_next    &
             ,nelt(1)                    &
             ,mpi_double_precision       &
             ,0,par%comm2d,ierr)

      nelt = ubound(ogcm_period_prev)   &
            -lbound(ogcm_period_prev)+1
      call mpi_bcast(ogcm_period_prev    &
             ,nelt(1)                    &
             ,mpi_double_precision       &
             ,0,par%comm2d,ierr)

      nelt = ubound(ogcm_rec_prev)   &
            -lbound(ogcm_rec_prev)+1
      call mpi_bcast(ogcm_rec_prev    &
             ,nelt(1)                 &
             ,mpi_integer             &
             ,0,par%comm2d,ierr)
#endif

      do var_count_=1,var_num
       ogcm_rec_next(var_count_)=ogcm_rec_prev(var_count_)+1
      enddo

!     do k=1,var_num
!     write(6,*)'rank var ogcm_readtime_next ',par%rank,k,ogcm_readtime_next(k)/3600.
!     write(6,*)'rank var ogcm_readtime_prev ',par%rank,k,ogcm_readtime_prev(k)/3600.
!     write(6,*)'rank var ogcm_period_next   ',par%rank,k,ogcm_period_next(k)/3600.
!     write(6,*)'rank var ogcm_period_prev   ',par%rank,k,ogcm_period_prev(k)/3600.
!     write(6,*)'rank var ogcm_rec_prev      ',par%rank,k,ogcm_rec_prev(k)
!     write(6,*)'rank var ogcm_rec_next      ',par%rank,k,ogcm_rec_next(k)
!     enddo
!     stop'kokot'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cas ogcm = nemo_z
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)  !23-06-13
#endif
      call ogcm_stop('Stop or not?',2)
      return
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cas ogcm = sympa
      if(obc_ogcm_type(1:5)=='sympa') then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        var_count_=1                                                   !25-08-14
        nc=1
        ogcm_rec_max=nc                                                !07-09-09
        open(unit=3,file=obcfile(1,1))

! Combien y a t'il de lignes dans le fichier?
         linemax_=0
 2834    read(3,'(a)',end=2832)texte250 ; linemax_=linemax_+1 ; goto 2834
 2832    rewind(3)
         if(linemax_==0) then !>>>>> !23-04-14
          write(6,'(a,a,a)')'File ',trim(obcfile(1,1)),' is empty'
          stop ' stop ogcm_lire_les_listes erreur 2205'
         endif                !>>>>>

! Chaque proc va traiter une fraction du fichier. On commence par derouler les lignes
! qui ne concernent pas le proc
         do loop2_=1,int(real(par%rank  )/real(nbdom)*linemax_)
          read(3,*) ! lire ces lignes pour rien
         enddo

        write(texte60,'(a,i0)')trim(tmpdirname)//'tmpfile',par%rank
        liste_ogcm_uvts=trim(tmpdirname)//''//'binrec'//'_liste_ogcm_uvts'        !13-10-10
        liste_ogcm_grid=trim(tmpdirname)//''//'binrec'//'_liste_ogcm_grid'        !13-10-10

        open(unit=4,file=texte60,status='REPLACE') !20-10-14

        do loop2_=int(real(par%rank  )/real(nbdom)*linemax_)+1    &
                 ,int(real(par%rank+1)/real(nbdom)*linemax_)


         read(3,'(a)',end=232)texte250

! Open the netcdf ogcm file and get the time_counter dimension then store in the binary file as much as time_counter:
! so that each line of the binary file corresponds to the ogcm output periodicity:
            max_time_counter=1
            status=nf_open(trim(texte250),nf_nowrite,ncid_)
            if(status/=0) then
               write(6,'(a,a)')'Echec ouverture ',trim(texte250) !25-08-14
               stop 'error nf_open in subroutine ogcm_lire_les_listes'
            endif
                         status=nf_inq_dimid(ncid_,'time_counter',dim_t_id)
            if(status/=0)status=nf_inq_dimid(ncid_,'TIME_COUNTER',dim_t_id)
            if(status/=0)status=nf_inq_dimid(ncid_,'time',dim_t_id)
            if(status/=0)status=nf_inq_dimid(ncid_,'TIME',dim_t_id)
            if(status==0)status=nf_inq_dimlen(ncid_,dim_t_id,max_time_counter)

                         status=nf_inq_varid(ncid_,'time',var_id)
!           if(status/=0)stop 'erreur nf_inq_varid time module_ogcm'
            if(status/=0) then !>>>
             write(6,'(a,a)')'err nf_inq_varid time in file: ',trim(texte250) ; flag_stop=1 ; goto 2731 !24-03-18
            endif              !>>>
            txt_units=''
            status=nf_get_att_text(ncid_,var_id,'units',txt_units)
            if(status/=0)stop 'erreur nf_get_att_text units ogcm'
            k=index(txt_units,'since')
            read(txt_units(k+6:k+9),*)year_
            read(txt_units(k+11:k+12),*)month_
            read(txt_units(k+14:k+15),*)day_
            read(txt_units(k+17:k+18),*)hour_
            read(txt_units(k+20:k+21),*)minute_
            read(txt_units(k+23:k+24),*)second_
            call datetokount(year_,month_,day_,hour_,minute_,second_) ! donne elapsedtime_out

            do loop_=1,max_time_counter

             status=nf_inq_varid(ncid_,'time',var_id)
             status=nf_get_vara_double(ncid_,var_id,loop_,1,x1)
             if(status/=0)stop ' stop module_ogcm nf_get_vara_double x1'
             txt_units=''
             status=nf_get_att_text(ncid_,var_id,'units',txt_units)
             if(status/=0)stop 'erreur nf_get_att_text units ogcm'
             x2=-999.
             if(index(txt_units,'days')/=0)x2=86400.
             if(index(txt_units,'hours')/=0)x2=3600.
             if(index(txt_units,'seconds')/=0)x2=1.
             if(x2==-999.)stop 'offline unites time incorrectes'
             ogcm_readtime_next(var_count_)=elapsedtime_out+x1*x2

             status=nf_inq_varid(ncid_,'cumulativetime',var_id)
             if(status/=0)stop ' stop1 module_ogcm cumulativetime'
             status=nf_get_vara_double(ncid_,var_id,loop_,1,x1)
             if(status/=0)stop ' stop2 module_ogcm cumulativetime'
             txt_units=''
             status=nf_get_att_text(ncid_,var_id,'units',txt_units)
             if(status/=0)stop 'erreur nf_get_att_text units ogcm'
             x2=-999.
             if(index(txt_units,'days')/=0)x2=86400.
             if(index(txt_units,'hours')/=0)x2=3600.
             if(index(txt_units,'seconds')/=0)x2=1.
             if(x2==-999.)stop 'offline unites time incorrectes'
             ogcm_period_next(var_count_)=x1*x2
             ogcm_readtime_next(var_count_)=ogcm_readtime_next(var_count_)-0.5*ogcm_period_next(var_count_) !23-04-14
!            write(6,*)'ogcm_period_next(var_count_) ',ogcm_period_next(var_count_)
!            write(6,*)'ogcm_readtime_next(var_count_)',ogcm_readtime_next(var_count_)

!            write(4,rec=nc)texte250,loop_,ogcm_readtime_next(var_count_),ogcm_period_next(var_count_)

             write(4,'(a)')trim(texte250)
             write(4,*)loop_
             write(4,*)ogcm_readtime_next(var_count_),ogcm_period_next(var_count_)

             if(ogcm_readtime_next(var_count_)<=elapsedtime_now) then          !----->
              ogcm_rec_prev(var_count_)=nc
              ogcm_period_prev(var_count_)=ogcm_period_next(var_count_)
             endif                                                 !----->

            enddo
            status=nf_close(ncid_)

         enddo ! fin de boucle sur loop2_
  2731   call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
         if(k0/=0)stop 'ERREUR OGCM LISTES BINAIRES' !09-05-18


  232   close(4)
        close(3)

! C'est le proc zero qui a la mission d'assembler tous les fichiers tmp en un seul fichier liste binaire
! La barriere suivante permet de s'assurer que tous les fichiers individuels ont bien ete fait au moment
! d'attaquer la concatenation
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)  !23-06-13
#endif

        if(par%rank==0) then !00000000000000000000>

         open(unit=3,file=liste_ogcm_uvts & !13-10-10
                   ,access='direct'       &
                   ,recl=541              & !28-04-16
                   ,form='unformatted')

        nc=1
        do loop1_=0,nbdom-1
         write(texte60,'(a,i0)')trim(tmpdirname)//'tmpfile',loop1_
         open(unit=4,file=texte60)
 2017    tm1_=tm2_                           !20-05-18
         read(4,'(a)',end=2003)texte250
         read(4,*,end=2003)loop_
         read(4,*,end=2003)ogcm_readtime_next(var_count_),ogcm_period_next(var_count_)
! Detection de listes desordonnees:                !20-05-18
               tm2_=ogcm_readtime_next(var_count_) !20-05-18
               if(nc==1)tm1_=tm2_                  !20-05-18
               if(tm2_<=tm1_.and.nc/=1) then  !m°v°m>
                write(6,'(a,a,a)')'File ',trim(texte250) &
                ,' does not respect the chronological order'
                write(6,'(a,a)')'List:',trim(texte60)
                write(6,*)'tm1_,tm2_',tm1_,tm2_
!               stop ' Stop module_ogcm'
                flag_stop=1
                goto 2820
               endif                          !m°v°m>

          write(3,rec=nc)texte250,loop_,ogcm_readtime_next(var_count_) &
                                       ,ogcm_period_next(var_count_)   &
                                       ,ogcmgridcode !16-05-19

          if(ogcm_readtime_next(var_count_)<=elapsedtime_now) then !----->
               ogcm_rec_prev(var_count_)=nc
          ogcm_readtime_prev(var_count_)=ogcm_readtime_next(var_count_)
            ogcm_period_prev(var_count_)=  ogcm_period_next(var_count_)
          endif                                                    !----->

          ogcm_rec_max=max0(ogcm_rec_max,nc)
          nc=nc+1
         goto 2017
 2003    close(4)
        enddo ! loop1_

        close(3)

      endif                !00000000000000000000>
 2820 call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0)stop 'Err 2820 ogcm file list'


#ifdef parallele
      nelt = ubound(ogcm_readtime_next)   &
            -lbound(ogcm_readtime_next)+1
      call mpi_bcast(ogcm_readtime_next    &
             ,nelt(1)                      &
             ,mpi_double_precision         &
             ,0,par%comm2d,ierr)

      nelt = ubound(ogcm_readtime_prev)   &
            -lbound(ogcm_readtime_prev)+1
      call mpi_bcast(ogcm_readtime_prev    &
             ,nelt(1)                      &
             ,mpi_double_precision         &
             ,0,par%comm2d,ierr)

      nelt = ubound(ogcm_period_next)   &
            -lbound(ogcm_period_next)+1
      call mpi_bcast(ogcm_period_next    &
             ,nelt(1)                    &
             ,mpi_double_precision       &
             ,0,par%comm2d,ierr)

      nelt = ubound(ogcm_period_prev)   &
            -lbound(ogcm_period_prev)+1
      call mpi_bcast(ogcm_period_prev    &
             ,nelt(1)                    &
             ,mpi_double_precision       &
             ,0,par%comm2d,ierr)

      nelt = ubound(ogcm_rec_prev)   &
            -lbound(ogcm_rec_prev)+1
      call mpi_bcast(ogcm_rec_prev    &
             ,nelt(1)                 &
             ,mpi_integer             &
             ,0,par%comm2d,ierr)
#endif

      do var_count_=1,var_num
       ogcm_rec_next(var_count_)=ogcm_rec_prev(var_count_)+1
      enddo

! Dans le cas sympa il y a homogeneite temporelle de toutes les variables donc
! valeurs identiques quelque soit var_count_
      do var_count_=2,var_num !25-08-14
       ogcm_readtime_next(var_count_)=ogcm_readtime_next(1)
       ogcm_readtime_prev(var_count_)=ogcm_readtime_prev(1)
       ogcm_period_next(var_count_)=  ogcm_period_next(1)
       ogcm_period_prev(var_count_)=  ogcm_period_prev(1)
       ogcm_rec_prev(var_count_)=     ogcm_rec_prev(1)
       ogcm_rec_next(var_count_)=     ogcm_rec_next(1)
      enddo

      if(par%rank==0) then !00000000000000000000>

       open(unit=3,file=obcfile(2,1))
       open(unit=4,file=liste_ogcm_grid           & !13-10-10
                  ,access='direct'                &
                  ,recl=541                       & !28-04-16
                  ,form='unformatted')
       read(3,'(a)')texte250
       write(4,rec=1)texte250
       close(4)
       close(3)

      endif                !00000000000000000000>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cas ogcm = sympa
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)  !23-06-13
#endif
      return
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cas ogcm = ncom
      if(obc_ogcm_type(1:4)=='ncom') then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     do k1=1,3 ! debut de boucle K1
      do k1=1,4 ! debut de boucle K1
        nc=1
        ogcm_rec_max=nc                                                    !07-09-09
        open(unit=3,file=obcfile(k1,1))

        liste_ogcm_t=trim(tmpdirname)//''//'binrec'//'_liste_ogcm_t' !13-10-10
        liste_ogcm_s=trim(tmpdirname)//''//'binrec'//'_liste_ogcm_s' !13-10-10
        liste_ogcm_u=trim(tmpdirname)//''//'binrec'//'_liste_ogcm_u' !13-10-10
        liste_ogcm_v=trim(tmpdirname)//''//'binrec'//'_liste_ogcm_v' !13-10-10
        if(k1==1)texte30=trim(liste_ogcm_t)         !13-10-10
        if(k1==2)texte30=trim(liste_ogcm_s)         !13-10-10
        if(k1==3)texte30=trim(liste_ogcm_u)         !13-10-10
        if(k1==4)texte30=trim(liste_ogcm_v)         !13-10-10

        if(par%rank==0)open(unit=4,file=texte30  &
                           ,access='direct'      &
                           ,recl=541             & !28-04-16
                           ,form='unformatted')
  131    read(3,'(a)',end=132)texte250
         if(par%rank==0)write(4,rec=nc)texte250
         ogcm_rec_max=max0(ogcm_rec_max,nc)
         nc=nc+1
         goto 131

  132   if(par%rank==0)close(4)
        close(3)

      enddo ! fin de boucle K1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cas ogcm = ncom
#ifdef parallele
      call mpi_barrier(par%comm2d,k_out)  !23-06-13
#endif
      return
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return
 2176 stop 'stop error opening binary ogcm file list'
      end subroutine ogcm_lire_les_listes

!------------------------------------------------------------------------------

      subroutine ogcm_hr_to_lr_coord(k_)
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
      double precision deci_(0:1,0:1),decj_(0:1,0:1)    &
                        ,dy_(0:1,0:1),  dx_(0:1,0:1)    &
                        ,dlon_di_,dlon_dj_,dlon_dm_
      integer k_,loop1_   ,loop2_
#ifdef synopsis
       subroutinetitle='ogcm_hr_to_lr_coord'
       subroutinedescription= &
          'For all horizontal locations (i,j) on the S grid,'      &
       //' finds the corresponding grid indexes on the OGCM grid.' &
       //' Stores the result in "hr_to_lr_coord_" files.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(k_==1)open(unit=3,file=trim(tmpdirname)//'hr_to_lr_coord_1_'//dom_c//'.out')
      if(k_==2)open(unit=3,file=trim(tmpdirname)//'hr_to_lr_coord_2_'//dom_c//'.out')
      if(k_==3)open(unit=3,file=trim(tmpdirname)//'hr_to_lr_coord_3_'//dom_c//'.out')

      do j=0,jmax+1                                                     !18-10-09
      do i=0,imax+1

      if(mask_t(i,j,kmax)==1) then !pmx> !15-11-17

      x2=real(max_x/2)
      x3=real(max_y/2)
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

      dlon_di_= ogcm_lon(i0+1,j0  )- ogcm_lon(i0-1,j0  )
      dlon_dj_= ogcm_lon(i0  ,j0+1)- ogcm_lon(i0  ,j0-1)
      dlon_dm_=rad2deg*lon_t(i  ,j)- ogcm_lon(i0  ,j0)

      if(dlon_di_<-180.)dlon_di_=dlon_di_+360.
      if(dlon_di_> 180.)dlon_di_=dlon_di_-360.
      if(dlon_dj_<-180.)dlon_dj_=dlon_dj_+360.
      if(dlon_dj_> 180.)dlon_dj_=dlon_dj_-360.
      if(dlon_dm_<-180.)dlon_dm_=dlon_dm_+360.
      if(dlon_dm_> 180.)dlon_dm_=dlon_dm_-360.

! Determinant principal:
      x1=( dlon_di_                                             &
          *( ogcm_lat(i0  ,j0+1)- ogcm_lat(i0  ,j0-1))          &
          -( ogcm_lat(i0+1,j0  )- ogcm_lat(i0-1,j0  ))          &
          *dlon_dj_ )*0.25

       if(x1==0) then !debug> !15-11-17
        write(10+par%rank,*)x1,i,j,i0,j0,dlon_di_,ogcm_lat(i0  ,j0+1)-ogcm_lat(i0  ,j0-1),ogcm_lat(i0+1,j0  )-ogcm_lat(i0-1,j0  ),dlon_dj_
        stop 'Err 2932, x1==0 in ogcm_hr_to_lr_coord . See fort files'
       endif          !debug>

      deci_(i2,j2)=min(max(                                  &
       i0+( dlon_dm_                                            &
          *( ogcm_lat(i0  ,j0+1)- ogcm_lat(i0,j0-1))            &
          -(rad2deg*lat_t(i,j)  - ogcm_lat(i0,j0))              &
          *dlon_dj_)/x1*0.5    &
                   ,2.0001d0),max_x-1.0001d0)

      decj_(i2,j2)=min(max(                                  &
       j0+( dlon_di_                                            &
          *(rad2deg*lat_t(i,j)  - ogcm_lat(i0  ,j0))            &
          -( ogcm_lat(i0+1,j0  )- ogcm_lat(i0-1,j0))            &
          *dlon_dm_   )/x1*0.5    &
                   ,2.0001d0),max_y-1.0001d0)

!      if(par%rank==0)write(6,*)deci_(i2,j2),decj_(i2,j2)
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
         if(par%rank==0)write(6,*)'par%rank=',par%rank
         if(par%rank==0)write(6,*)'(i,j)=   ',i,j
         if(par%rank==0)write(6,*)'deci decj',deci,decj
         stop 'hr_to_lr ne converge pas dans le forfait'
! Si le calcul ne converge pas car la grille de l'ogcm est
! discontinue on peut essayer l'algo suivant: !21-06-15
         sum1=small1 ; sum2=0. ; sum3=0.
         do j3=1,max_y
! dist0 est une distance de reference mesuree entre (1,j3) et (2,j3) de l'ogcm
          call lonlat2distance(ogcm_lon(1,j3)*deg2rad,ogcm_lat(1,j3)*deg2rad &
                              ,ogcm_lon(2,j3)*deg2rad,ogcm_lat(2,j3)*deg2rad,dist0)
         do i3=1,max_x
          call lonlat2distance(lon_t(i,j),lat_t(i,j),ogcm_lon(i3,j3)*deg2rad,ogcm_lat(i3,j3)*deg2rad,dist1)
          sum1=sum1+exp(-dist1/dist0)
          sum2=sum2+exp(-dist1/dist0)*i3
          sum3=sum3+exp(-dist1/dist0)*j3
         enddo
         enddo
         deci=sum2/sum1 ; decj=sum3/sum1
         goto 2908
       endif          !!!!!!>
       goto 1456
      endif

! ETAPE 2: CALCULER L'ANGLE D'ORIENTATION LOCALE DE LA GRILLE ORCA
 2908 continue
      i1=int(deci) ; j1=int(decj)
      rapi=deci-i1 ; rapj=decj-j1
      do j2=0,1
      do i2=0,1
       i0=i1+i2
       j0=j1+j2

       dlon_di_= ogcm_lon(i0+1,j0  )- ogcm_lon(i0-1,j0  )
       if(dlon_di_<-180.)dlon_di_=dlon_di_+360.
       if(dlon_di_> 180.)dlon_di_=dlon_di_-360.

       dy_(i2,j2)=  ogcm_lat(i0+1,j0)- ogcm_lat(i0-1,j0)
       dx_(i2,j2)=dlon_di_           &
                            *cos( ogcm_lat(i0,j0)*deg2rad)
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

!     ij2meteo_i(i,j)=deci
!     ij2meteo_j(i,j)=decj
!     ij2meteo_teta(i,j)=atan2(y1,x1)*rad2deg

!     write(66,*)i,j,ij2meteo_i(i,j),ij2meteo_j(i,j),ij2meteo_teta(i,j)
!     write(66,*)lon_t(i,j)*rad2deg,lat_t(i,j)*rad2deg
!     write(67,*) ogcm_lon(nint(deci),nint(decj)) &
!               , ogcm_lat(nint(deci),nint(decj))

!     write(3,'(2(i4,1x),3(f9.4,1x))')i,j,deci,decj,atan2(y1,x1)*rad2deg
! Pour la conservation de mpi on calcule les relations de passages dans
! la grande grille ogcm commune a tous les procs. En prevision du
! d'extraction on retranche ici le shift (ogcmzoom_istr,ogcmzoom_,str)
!     write(3,'(2(i4,1x),3(f9.4,1x))')i,j                  &
!                                    ,deci-ogcmzoom_istr+1 &
!                                    ,decj-ogcmzoom_jstr+1 &
!                                    ,atan2(y1,x1)*rad2deg
!     write(3,'(2(i5,1x),5(e14.7,1x))')i,j                 &
!                                    ,deci-ogcmzoom_istr+1 &
!                                    ,decj-ogcmzoom_jstr+1 &
!                                    ,atan2(y1,x1)*rad2deg &
!                                    ,deci-int(deci)      &
!                                    ,decj-int(decj)

      write(3,'(4(i5,1x),3(e14.7,1x))')i,j                      &
                                     ,int(deci)-ogcmzoom_istr+1 &
                                     ,int(decj)-ogcmzoom_jstr+1 &
                                     ,atan2(y1,x1)*rad2deg      &
                                     ,deci-int(deci)            &
                                     ,decj-int(decj)

      else                         !pmx> !15-11-17

      write(3,'(4(i5,1x),3(e14.7,1x))')i,j,1,1,0.,0.001,0.001 ! valeurs bidons dans le masque continental

      endif                        !pmx> !15-11-17


      enddo
      enddo

      close(3)

      end subroutine ogcm_hr_to_lr_coord

!------------------------------------------------------------------------------
#ifdef bidon
      subroutine ogcm_zero_divergence(t_   )
      use module_principal
      use module_parallele !#MPI
      use module_global
      implicit none
      integer t_
#ifdef synopsis
       subroutinetitle='ogcm_zero_divergence'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      call allocate_global('a','glob_u_x            ',iglb,jglb,0)     !08-10-09
      call allocate_global('a','glob_v_y            ',iglb,jglb,0)
      call allocate_global('a','glob_h_u            ',iglb,jglb,0)
      call allocate_global('a','glob_h_v            ',iglb,jglb,0)
      call allocate_global('a','glob_sf0_r          ',iglb,jglb,0)
      call allocate_global('a','glob_sf1_r          ',iglb,jglb,0)
      call allocate_global('a','glob_mask           ',iglb,jglb,0)

!_________________________________
! PREMIERE ESTIMATION DU TRANSPORT:
      do j=1,jmax+1
      do i=1,imax+1
       i1=i+par%timax(1)
       j1=j+par%tjmax(1)
        glob_h_u(i1,j1)=  mask_u(i,j,kmaxp1)*h_u(i,j)
        glob_u_x(i1,j1)=0.
        glob_h_v(i1,j1)=  mask_v(i,j,kmaxp1)*h_v(i,j)
        glob_v_y(i1,j1)=0.
       glob_sf0_r(i1,j1)=0.
       glob_sf1_r(i1,j1)=0.

      enddo
      enddo

! Attention dans allocate on dimensionne glob_mask   depuis i=j=0 donc
! la boucle doit partir de zero si on veut que gather fonctionne...
      do j=1,jmax+1
      do i=1,imax+1
       i1=i+par%timax(1)
       j1=j+par%tjmax(1)
       glob_mask  (i1,j1)=  mask_f(i,j,kmaxp1)
      enddo
      enddo

      do k=1,kmaxp1-1
      do j=1,jmax+1
      do i=1,imax+1
       i1=i+par%timax(1)
       j1=j+par%tjmax(1)
! Composante du transport*DY axe Oi:
       glob_u_x(i1,j1)=                                                 &
       glob_u_x(i1,j1)                                                  &
                +    dy_u(i,j)*dz_u(i,j,k,1)*velobc_u(i,j,k,t_   )      &
                                            *  mask_u(i,j,k)
! Composante du transport*DX axe Oj:
       glob_v_y(i1,j1)=                                                 &
       glob_v_y(i1,j1)                                                  &
                +    dx_v(i,j)*dz_v(i,j,k,1)*velobc_v(i,j,k,t_   )      &
                                            *  mask_v(i,j,k)
      enddo
      enddo
      enddo



! Tous les sous-domaines "voient" le grand domaine:
      ind(1)=par%timax(1)+1 ; ind(2)=par%timax(1)+imax+1
      ind(3)=par%tjmax(1)+1 ; ind(4)=par%tjmax(1)+jmax
      call  par_gatherall_2d(glob_h_u,                                  &
                      lbound(glob_h_u),                                 &
                      ubound(glob_h_u),                                 &
                             ind,par%nbdom)                  !08-09-09

      ind(1)=par%timax(1)+1 ; ind(2)=par%timax(1)+imax
      ind(3)=par%tjmax(1)+1 ; ind(4)=par%tjmax(1)+jmax+1
      call par_gatherall_2d(glob_h_v,                                   &
                     lbound(glob_h_v),                                  &
                     ubound(glob_h_v),                                  &
                            ind,par%nbdom)                  !08-09-09

      ind(1)=par%timax(1)+1 ; ind(2)=par%timax(1)+imax+1
      ind(3)=par%tjmax(1)+1 ; ind(4)=par%tjmax(1)+jmax
      call par_gatherall_2d(glob_u_x,                                   &
                     lbound(glob_u_x),                                  &
                     ubound(glob_u_x),                                  &
                            ind,par%nbdom)                  !08-09-09

      ind(1)=par%timax(1)+1 ; ind(2)=par%timax(1)+imax
      ind(3)=par%tjmax(1)+1 ; ind(4)=par%tjmax(1)+jmax+1
      call par_gatherall_2d(glob_v_y,                                   &
                     lbound(glob_v_y),                                  &
                     ubound(glob_v_y),                                  &
                            ind,par%nbdom)                  !08-09-09

      ind(1)=par%timax(1)+1 ; ind(2)=par%timax(1)+imax+1
      ind(3)=par%tjmax(1)+1 ; ind(4)=par%tjmax(1)+jmax+1
      call par_gatherall_2d(glob_mask  ,                                &
                     lbound(glob_mask  ),                               &
                     ubound(glob_mask  ),                               &
                            ind,par%nbdom)                  !08-09-09



!______________________________________________________________
! ASSURER QUE LE BILAN NET SUR LES FRONTIERES OUVERTES EST NUL:
! Debut:
!______________________________________________________________

      som0=0.
      som2=0.
      do i=1,iglb
       som0=som0+glob_h_v(i,jglb+1)                                     &
                +glob_h_v(i,1     )
       som2=som2+glob_v_y(i,jglb+1)                                     &
                -glob_v_y(i,1     )
      enddo
      do j=1,jglb
       som0=som0+glob_h_u(iglb+1,j)                                     &
                +glob_h_u(1     ,j)
       som2=som2+glob_u_x(iglb+1,j)                                     &
                -glob_u_x(1     ,j)
      enddo

      if(par%rank==0)write(6,*)'bilan obc avant:',som2
      lagrange_ssh=-som2/max(som0,small2)


      som2=0.
      do i=1,iglb

       glob_v_y(i,jglb+1)=glob_v_y(i,jglb+1)                            &
            +lagrange_ssh*glob_h_v(i,jglb+1)

       glob_v_y(i,1)=glob_v_y(i,1)                                      &
       -lagrange_ssh*glob_h_v(i,1)

       som2=som2+glob_v_y(i,jglb+1)                                     &
                -glob_v_y(i,1     )

      enddo

      do j=1,jglb

         glob_u_x(iglb+1,j)=glob_u_x(iglb+1,j)                          &
              +lagrange_ssh*glob_h_u(iglb+1,j)

         glob_u_x(1,j)=glob_u_x(1,j)                                    &
         -lagrange_ssh*glob_h_u(1,j)

       som2=som2+glob_u_x(iglb+1,j)                                     &
                -glob_u_x(1     ,j)

      enddo
      if(par%rank==0)write(6,*)'bilan obc apres:',som2

! Reset des tableaux de fonction de courant:
      do j=1,jglb+1
      do i=1,iglb+1
       glob_sf0_r(i,j)=0.
       glob_sf1_r(i,j)=0.
      enddo
      enddo

! Deduire du courant la fonction de courant sur les bords:
      i=1
      do j=1,jglb
        glob_sf0_r(i,j+1)=glob_sf0_r(i,j)-glob_u_x(i,j)
        glob_sf1_r(i,j+1)=glob_sf0_r(i,j+1)
      enddo

      j=jglb+1
      do i=1,iglb
        glob_sf0_r(i+1,j)=glob_sf0_r(i,j)+glob_v_y(i,j)
        glob_sf1_r(i+1,j)=glob_sf0_r(i+1,j)
      enddo

      i=iglb+1
      do j=jglb,1,-1
        glob_sf0_r(i,j)=glob_sf0_r(i,j+1)+glob_u_x(i,j)
        glob_sf1_r(i,j)=glob_sf0_r(i,j)
      enddo

      j=1
      do i=iglb,1,-1
        glob_sf0_r(i,j)=glob_sf0_r(i+1,j)-glob_v_y(i,j)
        glob_sf1_r(i,j)=glob_sf0_r(i,j)
      enddo
      if(par%rank==0)write(6,*)'fonction de courant au bout d''un tour:'               &
       ,glob_sf1_r(1,1)

!______________________________________________________________
! ASSURER QUE LE BILAN NET SUR LES FRONTIERES OUVERTES EST NUL:
! FIN.
!______________________________________________________________


!____________________________________________________________
! PREMIER AJUSTEMENT SANS DISTINCTION DES POINTS TERRE ET MER
      if(par%rank==0)write(6,*)'passage 1: sans distinction du masque'

      const1=1.99/4.
      do loop1=1,3000

      do j=2,jglb
      do i=2,iglb

       glob_sf1_r(i,j)=glob_sf0_r(i,j)+const1*(                         &

               glob_sf0_r(i+1,j  )                                      &
              +glob_sf1_r(i-1,j  )                                      &
              +glob_sf0_r(i  ,j+1)                                      &
              +glob_sf1_r(i  ,j-1)                                      &
           -4.*glob_sf0_r(i  ,j  )                                      &

        + glob_u_x(i,j)-glob_u_x(i,j-1)                                 &
        - glob_v_y(i,j)+glob_v_y(i-1,j)                                 &
                                                )

      enddo
      enddo

      sum1=0.
      do j=2,jglb
      do i=2,iglb
      sum1=sum1+abs(glob_sf1_r(i,j)-glob_sf0_r(i,j))
                    glob_sf0_r(i,j)=glob_sf1_r(i,j)
      enddo
      enddo
      if(mod(loop1,500).eq.0.and.par%rank==0)write(6,*)sum1

      enddo ! fin de boucle LOOP1

!____________________________________________________________
! ASSURER QUE LA FONCTION DE COURANT EST CONSTANTE EN TERRE:
      if(par%rank==0)write(6,*)'passage 2: points en terre uniquement'

      const1=1.99/4.
      do loop1=1,3000

      do j=2,jglb
      do i=2,iglb

       if(glob_mask  (i,j).eq.0) then !>>>>>>>>>>>

       glob_sf1_r(i,j)=glob_sf0_r(i,j)+const1*(                         &

              (glob_sf0_r(i+1,j  )-glob_sf0_r(i  ,j  ))                 &
          *(1-glob_mask  (i+1,j  ))                                     &

             +(glob_sf1_r(i-1,j  )-glob_sf0_r(i  ,j  ))                 &
          *(1-glob_mask  (i-1,j  ))                                     &

             +(glob_sf0_r(i  ,j+1)-glob_sf0_r(i  ,j  ))                 &
          *(1-glob_mask  (i  ,j+1))                                     &

             +(glob_sf1_r(i  ,j-1)-glob_sf0_r(i  ,j  ))                 &
          *(1-glob_mask  (i  ,j-1))                                     &
                                                )

       endif                          !>>>>>>>>>>>

      enddo
      enddo

      sum1=0.
      do j=2,jglb
      do i=2,iglb
      if(glob_mask  (i,j).eq.0) then !>>>>>>>>>>>
       sum1=sum1+abs(glob_sf1_r(i,j)-glob_sf0_r(i,j))
       glob_sf0_r(i,j)=glob_sf1_r(i,j)
      endif                          !>>>>>>>>>>>
      enddo
      enddo
      if(mod(loop1,500).eq.0.and.par%rank==0)write(6,*)sum1

      enddo ! fin de boucle LOOP1

!____________________________________________________________
! DEUXIEME AJUSTEMENT: POINTS DE MER UNIQUEMENT AVEC PRISE EN
! COMPTE DE LA CONDITION AUX LIMITES SOLIDES
      if(par%rank==0)write(6,*)'passage 3: points en mer uniquement'
      const1=1.99/4.
      do loop1=1,3000

      do j=2,jglb
      do i=2,iglb

       if(glob_mask  (i,j).eq.1) then !>>>>>>>>>>>

       glob_sf1_r(i,j)=glob_sf0_r(i,j)+const1*(                         &

               glob_sf0_r(i+1,j  )                                      &
              +glob_sf1_r(i-1,j  )                                      &
              +glob_sf0_r(i  ,j+1)                                      &
              +glob_sf1_r(i  ,j-1)                                      &
           -4.*glob_sf0_r(i  ,j  )                                      &

        + glob_u_x(i,j)-glob_u_x(i,j-1)                                 &
        - glob_v_y(i,j)+glob_v_y(i-1,j)                                 &
                                                )

       endif                          !>>>>>>>>>>>

      enddo
      enddo

      sum1=0.
      do j=2,jglb
      do i=2,iglb

      if(glob_mask  (i,j).eq.1) then !>>>>>>>>>>>

       sum1=sum1+abs(glob_sf1_r(i,j)-glob_sf0_r(i,j))
       glob_sf0_r(i,j)=glob_sf1_r(i,j)

      endif                          !>>>>>>>>>>>

      enddo
      enddo
      if(mod(loop1,500).eq.0.and.par%rank==0)write(6,*)sum1

      enddo ! fin de boucle LOOP1


!____________________________________________________________
! DEDUIRE LE TRANSPORT DE LA FONCTION DE COURANT:
      do j=1,jmax
      do i=1,imax+1
       i1=i+par%timax(1)
       j1=j+par%tjmax(1)
       xy_u(i,j,2)= glob_sf1_r(i1,j1  )                                 &
                   -glob_sf1_r(i1,j1+1)
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
       i1=i+par%timax(1)
       j1=j+par%tjmax(1)
       xy_v(i,j,2)=-glob_sf1_r(i1  ,j1)                                 &
                   +glob_sf1_r(i1+1,j1)
      enddo
      enddo

      call allocate_global('d','glob_u_x            ',iglb,jglb,0)     !08-10-09
      call allocate_global('d','glob_v_y            ',iglb,jglb,0)
      call allocate_global('d','glob_h_u            ',iglb,jglb,0)
      call allocate_global('d','glob_h_v            ',iglb,jglb,0)
      call allocate_global('d','glob_sf0_r          ',iglb,jglb,0)
      call allocate_global('d','glob_sf1_r          ',iglb,jglb,0)
      call allocate_global('d','glob_mask           ',iglb,jglb,0)

!____________________________________________________________
! INTRODUIRE LE NOUVEAU TRANSPORT DANS LE COURANT 3D
! La correction depend de z: en pratique elle est proportionnelle
! e l'intensite du courant
      do j=1,jmax
      do i=1,imax+1
         sum1=0.
         sum2=0.
         do k=kmin_u(i,j),kmaxp1-1
          sum1=sum1+abs(velobc_u(i,j,k,t_   ))                          &
                         *    dy_u(i,j)*dz_u(i,j,k,1)
          sum2=sum2+    dy_u(i,j)*dz_u(i,j,k,1)*velobc_u(i,j,k,t_   )
         enddo
         x1=(xy_u(i,j,2)-sum2)/max(sum1,small2)
         sum3=0.
         do k=kmin_u(i,j),kmaxp1-1
          velobc_u(i,j,k,t_   )=velobc_u(i,j,k,t_   )                   &
                        +x1*abs(velobc_u(i,j,k,t_   ))
          sum3=sum3+    dy_u(i,j)*dz_u(i,j,k,1)*velobc_u(i,j,k,t_   )
         enddo
      enddo
      enddo

      do j=1,jmax+1
      do i=1,imax
         sum1=0.
         sum2=0.
         do k=kmin_v(i,j),kmaxp1-1
          sum1=sum1+abs(velobc_v(i,j,k,t_   ))                          &
                         *    dx_v(i,j)*dz_v(i,j,k,1)
          sum2=sum2+    dx_v(i,j)*dz_v(i,j,k,1)*velobc_v(i,j,k,t_   )
         enddo
         x1=(xy_v(i,j,2)-sum2)/max(sum1,small2)

         sum3=0.
         do k=kmin_v(i,j),kmaxp1-1
          velobc_v(i,j,k,t_   )=velobc_v(i,j,k,t_   )                   &
                        +x1*abs(velobc_v(i,j,k,t_   ))
          sum3=sum3+    dx_v(i,j)*dz_v(i,j,k,1)*velobc_v(i,j,k,t_   )
         enddo
      enddo
      enddo


      call mpi_barrier(par%comm2d,k_out)      ! synchro processes


      end subroutine ogcm_zero_divergence
#endif
!------------------------------------------------------------------------------
      subroutine ogcm_repere_les_trous
      use module_principal
      use module_parallele
      use module_forcages
      implicit none
#ifdef synopsis
       subroutinetitle='ogcm_repere_les_trous'
       subroutinedescription= &
          'Identification of the OGCM land grid points that will' &
       //' be involved in the interpolation process.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0)write(6,*)'reperage des trous'

      open(unit=3, &
           file=trim(tmpdirname)//'hr_to_lr_coord_1_'//dom_c//'.out')
      k1=max_z
      do j=0,jmax+1
      do i=0,imax+1

!     read(3,*)i0,j0,deci,decj,x0,rapi,rapj
      read(3,*)i0,j0,i1,j1,x0,rapi,rapj
      if(i.ne.i0)stop 'erreur2 sur i'
      if(j.ne.j0)stop 'erreur2 sur j'

      if(  mask_t(i,j,kmaxp1).eq.1) then !>>>>>>>>>>>>

!      i1=int(deci)
!      j1=int(decj)

          if(ogcm_var3d(i1  ,j1  ,k1)> 1.e10)                           &
             ogcm_var3d(i1  ,j1  ,k1)=-1.e10
          if(ogcm_var3d(i1+1,j1  ,k1)> 1.e10)                           &
             ogcm_var3d(i1+1,j1  ,k1)=-1.e10
          if(ogcm_var3d(i1  ,j1+1,k1)> 1.e10)                           &
             ogcm_var3d(i1  ,j1+1,k1)=-1.e10
          if(ogcm_var3d(i1+1,j1+1,k1)> 1.e10)                           &
             ogcm_var3d(i1+1,j1+1,k1)=-1.e10

      endif                          !>>>>>>>>>>>>

      enddo
      enddo
      close(3)

      end subroutine ogcm_repere_les_trous

!------------------------------------------------------------------------------
      subroutine ogcm_grille_sous_le_fond
      use module_principal
      use module_parallele
      use module_forcages
      implicit none
#ifdef synopsis
       subroutinetitle='ogcm_grille_sous_le_fond'
       subroutinedescription=                            &
          'For all S grid points under the OGCM bottom,' &
       //' identifies the nearest OGCM grid point above' &
       //' the bottom.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! cette routine repere les points de Smodel qui sont sous le fond
! de l'ogcm. On cherche dans l'ogcm le point le plus proche dont
! la profondeur est en dessous du fond de Smodel. Ce point
! va donne un profil de reference pour extrapolation sous le fond

      if(par%rank==0)write(6,*)'reperage cas sous le fond'

      open(unit=3 &
      ,file=trim(tmpdirname)//'hr_to_lr_coord_1_'//dom_c//'.out')
      open(unit=4,file=trim(tmpdirname)//'sous_le_fond_'//dom_c//'.out')
      k1=max_z
      do j=0,jmax+1
      do i=0,imax+1

!     read(3,*)i0,j0,deci,decj,x0,rapi,rapj
      read(3,*)i0,j0,i1,j1,x0,rapi,rapj

      if(i.ne.i0)stop 'erreur2 sur i'
      if(j.ne.j0)stop 'erreur2 sur j'

      if(  mask_t(i,j,kmaxp1).eq.1) then !>>>>>>>>>>>>

!      i1=int(deci)
!      rapi=deci-i1
!      j1=int(decj)
!      rapj=decj-j1

         do loop1=0,10                                                  !27-09-10
          do j2=max0(j1-loop1,1),min0(j1+loop1,max_y)
          do i2=max0(i1-loop1,1),min0(i1+loop1,max_x)

           if(abs(ogcm_var3d(i2,j2,max_z))<1.e10)then !---->            !13-11-09

!           if(ogcm_z(i2,j2,1)<depth_t(i,j,kmin_w(i,j))) then !***>
            if(ogcm_z(i2,j2,1)<depth_w(i,j,kmin_w(i,j))) then !***> !23-11-11
             write(4,'(5(i4,1x))')i,j,i2,j2,loop1
             goto 1818
            endif                                             !***<

           endif                                      !----<

          enddo
          enddo
         enddo
            loop1=-1
            write(4,'(5(i4,1x))')i,j,i1,j1,loop1

 1818    continue

      endif                          !>>>>>>>>>>>>

      enddo
      enddo
      close(3)
      close(4)

      end subroutine ogcm_grille_sous_le_fond

!---------------------------------------------------------------------

      subroutine ogcm_get_dim(ncid_   ,txt_   )
      use module_principal
      use module_parallele
      implicit none
      include 'netcdf.inc'
      character*1 txt_
      integer ncid_
#ifdef synopsis
       subroutinetitle='ogcm_get_dim'
       subroutinedescription= &
       'Reads the dimensions in the OGCM netcdf file header'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

#ifdef bidon
! bidonnE le 18-05-19
      if(iteration3d/=0) then !oooooooooooooooooooooooooooooooooo>

         if(txt_=='m') then !mmmmmm>
          max_x=    ogcmzoom_iendt-ogcmzoom_istrt+1
          max_y=    ogcmzoom_jendt-ogcmzoom_jstrt+1
          max_x=max(ogcmzoom_iendu-ogcmzoom_istru+1,max_x)
          max_y=max(ogcmzoom_jendu-ogcmzoom_jstru+1,max_y)
          max_x=max(ogcmzoom_iendv-ogcmzoom_istrv+1,max_x)
          max_y=max(ogcmzoom_jendv-ogcmzoom_jstrv+1,max_y)
          return
         endif              !mmmmmm>
         if(txt_=='t') then !tttttt>
          max_x=ogcmzoom_iendt-ogcmzoom_istrt+1
          max_y=ogcmzoom_jendt-ogcmzoom_jstrt+1
          return
         endif              !tttttt>
         if(txt_=='u') then !uuuuuu>
          max_x=ogcmzoom_iendu-ogcmzoom_istru+1
          max_y=ogcmzoom_jendu-ogcmzoom_jstru+1
          return
         endif              !uuuuuu>
         if(txt_=='v') then !vvvvvv>
          max_x=ogcmzoom_iendv-ogcmzoom_istrv+1
          max_y=ogcmzoom_jendv-ogcmzoom_jstrv+1
          return
         endif              !vvvvvv>

         stop 'erreur argument ogcm_get_dim'
      endif                   !oooooooooooooooooooooooooooooooooo>
#endif

! Lecture des valeurs des dimensions

      if(obc_ogcm_type(1:5)=='sympa') then  !sssssssss>

      if(txt_   =='t') then !ttttttttttttttttttttttttttttttttttttttt>

                    status=nf_inq_dimid(ncid_   ,'x_zhl',dim_x_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'x_ZHL',dim_x_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'imax_t',dim_x_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'ni_t',dim_x_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'ni',dim_x_id)
       if(status/=0)stop 'get_dim erreur dim_x_id t'

                    status=nf_inq_dimid(ncid_   ,'y_zhl',dim_y_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'y_ZHL',dim_y_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'jmax_t',dim_y_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nj_t',dim_y_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nj',dim_y_id)
       if(status/=0)stop 'get_dim erreur dim_y_id t'

                    status=nf_inq_dimid(ncid_   ,'z_zhl',dim_z_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'z_ZHL',dim_z_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'kmax_t',dim_z_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nk_t',dim_z_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nk',dim_z_id)
       if(status/=0)stop 'get_dim erreur dim_z_id t'

      status=nf_inq_dimlen(ncid_,dim_x_id,max_x)
      status=nf_inq_dimlen(ncid_,dim_y_id,max_y)
      status=nf_inq_dimlen(ncid_,dim_z_id,max_z)
      endif                 !ttttttttttttttttttttttttttttttttttttttt>


      if(txt_   =='u') then !uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu>

                    status=nf_inq_dimid(ncid_   ,'x_xhl',dim_x_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'x_XHL',dim_x_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'imax_u',dim_x_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'ni_u',dim_x_id)
       if(status/=0)stop 'interp_ogcm erreur dim_x_id u'

                    status=nf_inq_dimid(ncid_   ,'y_xhl',dim_y_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'y_XHL',dim_y_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'jmax_u',dim_y_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nj_u',dim_y_id)
       if(status/=0)stop 'interp_ogcm erreur dim_y_id u'

                    status=nf_inq_dimid(ncid_   ,'z_xhl',dim_z_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'z_XHL',dim_z_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'kmax_u',dim_z_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nk_u',dim_z_id)
       if(status/=0)stop 'interp_ogcm erreur dim_z_id u'

      status=nf_inq_dimlen(ncid_   ,dim_x_id,max_x)
      status=nf_inq_dimlen(ncid_   ,dim_y_id,max_y)
      status=nf_inq_dimlen(ncid_   ,dim_z_id,max_z)
      endif                 !uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu>

      if(txt_   =='v') then !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv>

                    status=nf_inq_dimid(ncid_   ,'x_yhl',dim_x_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'x_YHL',dim_x_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'imax_v',dim_x_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'ni_v',dim_x_id)
       if(status/=0)stop 'interp_ogcm erreur dim_x_id v'
                    status=nf_inq_dimid(ncid_   ,'y_yhl',dim_y_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'y_YHL',dim_y_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'jmax_v',dim_y_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nj_v',dim_y_id)
       if(status/=0)stop 'interp_ogcm erreur dim_y_id v'
                    status=nf_inq_dimid(ncid_   ,'z_yhl',dim_z_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'z_YHL',dim_z_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'kmax_v',dim_z_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nk_v',dim_z_id)
       if(status/=0)stop 'interp_ogcm erreur dim_z_id v'

      status=nf_inq_dimlen(ncid_   ,dim_x_id,max_x)
      status=nf_inq_dimlen(ncid_   ,dim_y_id,max_y)
      status=nf_inq_dimlen(ncid_   ,dim_z_id,max_z)
      endif                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv>

      if(txt_   =='m') then !mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm>
! Trouver les dimensions maximum pour allocater une fois pour toutes
                    status=nf_inq_dimid(ncid_   ,'imax_t',dim_x_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'ni_t',dim_x_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'ni',dim_x_id)
       if(status/=0)stop 'get_dim erreur dim_x_id'
       status=nf_inq_dimlen(ncid_   ,dim_x_id,max_x)

                    status=nf_inq_dimid(ncid_   ,'imax_u',dim_x_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'ni_u',dim_x_id)
       if(status/=0)stop 'get_dim erreur dim_x_id'
       status=nf_inq_dimlen(ncid_   ,dim_x_id,i)
       max_x=max(max_x,i)

                    status=nf_inq_dimid(ncid_   ,'imax_v',dim_x_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'ni_v',dim_x_id)
       if(status/=0)stop 'get_dim erreur dim_x_id'
       status=nf_inq_dimlen(ncid_   ,dim_x_id,i)
       max_x=max(max_x,i)

                    status=nf_inq_dimid(ncid_   ,'jmax_t',dim_y_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nj_t',dim_y_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nj',dim_y_id)
       if(status/=0)stop 'get_dim erreur dim_y_id'
       status=nf_inq_dimlen(ncid_   ,dim_y_id,max_y)

                    status=nf_inq_dimid(ncid_   ,'jmax_u',dim_y_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nj_u',dim_y_id)
       if(status/=0)stop 'get_dim erreur dim_y_id'
       status=nf_inq_dimlen(ncid_   ,dim_y_id,i)
       max_y=max(max_y,i)

                    status=nf_inq_dimid(ncid_   ,'jmax_v',dim_y_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nj_v',dim_y_id)
       if(status/=0)stop 'get_dim erreur dim_y_id'
       status=nf_inq_dimlen(ncid_   ,dim_y_id,i)
       max_y=max(max_y,i)

                    status=nf_inq_dimid(ncid_   ,'kmax_t',dim_z_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nk_t',dim_z_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nk',dim_z_id)
       if(status/=0)stop 'interp_ogcm erreur dim_z_id'
       status=nf_inq_dimlen(ncid_   ,dim_z_id,max_z)

                    status=nf_inq_dimid(ncid_   ,'kmax_u',dim_z_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nk_u',dim_z_id)
       if(status/=0)stop 'interp_ogcm erreur dim_z_id'
       status=nf_inq_dimlen(ncid_   ,dim_z_id,i)
       max_z=max(max_z,i)

                    status=nf_inq_dimid(ncid_   ,'kmax_v',dim_z_id)
       if(status/=0)status=nf_inq_dimid(ncid_   ,'nk_v',dim_z_id)
       if(status/=0)stop 'interp_ogcm erreur dim_z_id'
       status=nf_inq_dimlen(ncid_   ,dim_z_id,i)
       max_z=max(max_z,i)
       if(par%rank==0)write(6,*)'Dimensions max=',max_x,max_y,max_z
      endif                 !mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm>

      endif                                 !sssssssss>

      if(obc_ogcm_type(1:6)=='nemo_z') then !nenenenen>

      if(txt_   =='t') then !ttttttttttttttttttttttttttttttttttttttt>
                     status=nf_inq_dimid(ncid_   ,'x',dim_x_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'longitude',dim_x_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'lon',dim_x_id) !16-06-18
        if(status/=0)stop 'erreur lecture dim_x_id 1'
                     status=nf_inq_dimid(ncid_   ,'y',dim_y_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'latitude',dim_y_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'lat',dim_y_id) !16-05-18
        if(status/=0)stop 'erreur lecture dim_y_id 1'
                     status=nf_inq_dimid(ncid_   ,'deptht',dim_z_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'z',dim_z_id)       !03-06-10
        if(status/=0)status=nf_inq_dimid(ncid_   ,'depth',dim_z_id)
        if(status/=0)stop 'erreur lecture dim_z_id 1470'              !03-06-10
        status=nf_inq_dimlen(ncid_   ,dim_x_id,max_x)
        status=nf_inq_dimlen(ncid_   ,dim_y_id,max_y)
        status=nf_inq_dimlen(ncid_   ,dim_z_id,max_z)
      endif                 !ttttttttttttttttttttttttttttttttttttttt>
      if(txt_   =='u') then !uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu>
                     status=nf_inq_dimid(ncid_   ,'x',dim_x_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'longitude',dim_x_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'lon',dim_x_id) !16-05-18
        if(status/=0)stop 'erreur lecture dim_x_id 2'
                     status=nf_inq_dimid(ncid_   ,'y',dim_y_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'latitude',dim_y_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'lat',dim_y_id) !16-05-18
        if(status/=0)stop 'erreur lecture dim_y_id 2'
                     status=nf_inq_dimid(ncid_   ,'depthu',dim_z_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'z',dim_z_id)         !03-06-10
        if(status/=0)status=nf_inq_dimid(ncid_   ,'deptht',dim_z_id)    !29-02-12
        if(status/=0)status=nf_inq_dimid(ncid_   ,'depth',dim_z_id)
        if(status/=0)stop 'erreur lecture dim_z_id 1471'                !03-06-10
        status=nf_inq_dimlen(ncid_   ,dim_x_id,max_x)
        status=nf_inq_dimlen(ncid_   ,dim_y_id,max_y)
        status=nf_inq_dimlen(ncid_   ,dim_z_id,max_z)
!     if(par%rank==0)write(6,*)'Routine get dim pour nemo u max_x...',max_x,max_y,max_z
      endif                 !uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu>
      if(txt_   =='v') then !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv>
                     status=nf_inq_dimid(ncid_   ,'x',dim_x_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'longitude',dim_x_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'lon',dim_x_id) !16-05-18
        if(status/=0)stop 'erreur lecture dim_x_id 3'
                     status=nf_inq_dimid(ncid_   ,'y',dim_y_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'latitude',dim_y_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'lat',dim_y_id) !16-05-18
        if(status/=0)stop 'erreur lecture dim_y_id 3'
                     status=nf_inq_dimid(ncid_   ,'depthv',dim_z_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'z',dim_z_id)         !03-06-10
        if(status/=0)status=nf_inq_dimid(ncid_   ,'deptht',dim_z_id)    !29-02-12
        if(status/=0)status=nf_inq_dimid(ncid_   ,'depth',dim_z_id)
        if(status/=0)stop 'erreur lecture dim_z_id 1472'                !03-06-10
        status=nf_inq_dimlen(ncid_   ,dim_x_id,max_x)
        status=nf_inq_dimlen(ncid_   ,dim_y_id,max_y)
        status=nf_inq_dimlen(ncid_   ,dim_z_id,max_z)
      endif                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv>
      if(txt_   =='m') then !mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm>
                     status=nf_inq_dimid(ncid_   ,'x',dim_x_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'longitude',dim_x_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'lon',dim_x_id) !16-05-18
        if(status/=0)stop 'erreur lecture dim_x_id 4'
                     status=nf_inq_dimid(ncid_   ,'y',dim_y_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'latitude',dim_y_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'lat',dim_y_id) !16-05-18
        if(status/=0)stop 'erreur lecture dim_y_id 4'
                     status=nf_inq_dimid(ncid_   ,'z',dim_z_id)             !03-06-10
        if(status/=0)status=nf_inq_dimid(ncid_   ,'deptht',dim_z_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'depthu',dim_z_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'depthv',dim_z_id)
        if(status/=0)status=nf_inq_dimid(ncid_   ,'depth',dim_z_id)
        if(status/=0)stop 'erreur lecture dim_z_id'
        status=nf_inq_dimlen(ncid_   ,dim_x_id,max_x)
        status=nf_inq_dimlen(ncid_   ,dim_y_id,max_y)
        status=nf_inq_dimlen(ncid_   ,dim_z_id,max_z)
      endif                 !mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm>

      endif                                 !nenenenen>

      if(obc_ogcm_type(1:4)=='ncom') then !ncomncom
                    status=nf_inq_dimid(ncid_   ,'lon',dim_x_id)
       if(status/=0)stop 'erreur lecture dim_x_id 5'
                    status=nf_inq_dimid(ncid_   ,'lat',dim_y_id)
       if(status/=0)stop 'erreur lecture dim_y_id'
                    status=nf_inq_dimid(ncid_   ,'depth',dim_z_id)             !03-06-10
       if(status/=0)stop 'erreur lecture dim_z_id'
      status=nf_inq_dimlen(ncid_   ,dim_x_id,max_x)
      status=nf_inq_dimlen(ncid_   ,dim_y_id,max_y)
      status=nf_inq_dimlen(ncid_   ,dim_z_id,max_z)
      endif                                 !ncomncom>

      end subroutine ogcm_get_dim

!------------------------------------------------------------------------------

      subroutine ogcm_repere_les_trous_ts                                 !19-11-10
      use module_principal
      use module_parallele
      use module_forcages
      implicit none
#ifdef synopsis
       subroutinetitle='ogcm_repere_les_trous_ts'
       subroutinedescription= &
          'Identification of the OGCM land grid points that will' &
       //' be involved in the interpolation process.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


! Ces lignes servent a rejeter le premier point au dessus du fond: !13-11-15
      do j=1,max_y
      do i=1,max_x
       k=1
       do while (ogcm_var3d(i,j,k)>1.e10.and.k<max_z)
        k=k+1
       enddo
! on sort de la boucle while avec le premier niveau au dessus du fond
! on le rejete:
       ogcm_var3d(i,j,k)=1.e11
      enddo
      enddo

      if(par%rank==0)write(6,*)'reperage des trous 3d'

      flag_stop=0
      open(unit=3 &
      ,file=trim(tmpdirname)//'hr_to_lr_coord_1_'//dom_c//'.out')
      do j=0,jmax+1
      do i=0,imax+1

      read(3,*)i0,j0,i1,j1,x0

!---- Section detection des bugs --------->
      if(i/=i0.or.j/=j0) then
       write(10+par%rank,*) &
       'Err 4004 i/=i0.or.j/=j0 ogcm_repere_les_trous_ts'
       flag_stop=1
      endif
      if(i1<1.or.i1+1>max_x.or.j1<1.or.j1+1>max_y) then
       write(10+par%rank,*) &
        'Err 4009 i1<1.or.i1+1>max_x.or.j1<1.or.j1+1>max_y' &
       ,' in ogcm_repere_les_trous_ts'
       write(10+par%rank,*)'max_x,max_y',max_x,max_y
       write(10+par%rank,*)'ogcm i1,j1',i1,j1
       write(10+par%rank,*)'symphonie i,j glob',i+par%timax(1),j+par%tjmax(1)
       write(10+par%rank,*)'OGCM grid extension is possibly too small'
       flag_stop=1
      endif
      if(flag_stop==1) goto 4015
!---- Section detection des bugs --------->

      if(  mask_t(i,j,kmaxp1).eq.1) then !>>>>>>>>>>>>

        do k1=1,max_z

          if(ogcm_var3d(i1  ,j1  ,k1)> 1.e10)  &
             ogcm_var3d(i1  ,j1  ,k1)=         &
            -ogcm_var3d(i1  ,j1  ,k1)
          if(ogcm_var3d(i1+1,j1  ,k1)> 1.e10)  &
             ogcm_var3d(i1+1,j1  ,k1)=         &
            -ogcm_var3d(i1+1,j1  ,k1)           
          if(ogcm_var3d(i1  ,j1+1,k1)> 1.e10)  &
             ogcm_var3d(i1  ,j1+1,k1)=         &
            -ogcm_var3d(i1  ,j1+1,k1)
          if(ogcm_var3d(i1+1,j1+1,k1)> 1.e10)  &
             ogcm_var3d(i1+1,j1+1,k1)=         &
            -ogcm_var3d(i1+1,j1+1,k1) 


         enddo

      endif                          !>>>>>>>>>>>>

      enddo
      enddo
 4015 close(3)
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0) stop 'Err 4046 in module_ogcm. See fort.xx files' !22-05-18

      end subroutine ogcm_repere_les_trous_ts

!------------------------------------------------------------------------------

      subroutine ogcm_fichier_bouchetrou_ts                         !19-11-10
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
#ifdef synopsis
       subroutinetitle='ogcm_fichier_bouchetrou_ts'
       subroutinedescription= &
          'Finds nearest replacing sea grid points for OGCM land grid' &
       //' points concerned by the interpolation process. Stores'      &
       //' coordinates in "bouchetrou" file.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0)write(6,*)'reperage des bouche trous 3d'

      open(unit=4,file=trim(tmpdirname)//'bouchetrou_ts'//dom_c//'.out')   !#MPI
      open(unit=3,file=trim(tmpdirname)//'bouchetrou_ssh'//dom_c//'.out')   !#MPI

      do 1862 j=1,max_y
      do 1862 i=1,max_x

! Attention la boucle k doit etre au centre
      i10=1
      do 1862 k=max_z,1,-1

       if(ogcm_var3d(i,j,k).lt.-1.e9)then !eeeeeeeeeeeeeee>
        ksecu=0
        dist1=1e10
 1864   continue

        do 1863 k1=0,1

         j0=max0(j-i10,1)
         j2=min0(j+i10,max_y)
         j3=k1*(j2-j0)+(1-k1)

         i0=max0(i-i10+k1,1)
         i2=min0(i+i10-k1,max_x)
         i3=(1-k1)*(i2-i0)+k1

         do 1861 j1=j0,j2,j3
         do 1861 i1=i0,i2,i3
          if(abs(ogcm_var3d(i1,j1,k)).lt.1.e9)then !%%%%%%%%%%%%%>
           dist2=sqrt(real(i-i1)**2+real(j-j1)**2)
           ksecu=1
             if(dist2.lt.dist1) then                    !>>>>>>>>>>>>>
              i4=i1
              j4=j1
              k4=k
              dist1=dist2
             endif                                      !>>>>>>>>>>>>>
          endif                                      !%%%%%%%%%%%%%>
 1861    continue ! fin de boucles i1,j1

 1863   continue  ! fin de boucle k1

! Si pas de bouche trou horizontal a proximite alors bouche trou vertical:
        if(i10 > ogcm_enlarged_bounds-2) then !---------> !17-07-14
! Le "-2" de ogcm_enlarged_bounds-2 est intuitif et doit servir a ce que le
! point de bouchage ne depende pas du nombre de proc mpi

       if(k==max_z) then  !07-03-17
           write(10+par%rank,*)par%rank,i,j
        stop 'No substitue on the surface, enlarge ogcm_enlarged_bounds'
       endif

! A condition de ne pas se retrouver ici en k=max_z, on est sur qu'il
! existera une valeur (bouchee ou non bouchee, peu importe) en
! ogcm_var3d(i,j,max_z) car le bouchage commence par le haut et 
! renseigne au fur et A mesure ogcm_var3d(i,j,k). Par consequent
! quand on se retrouve ici en k/=max_z on sait que i,j,k+1 est un
! boucheur possible. Je crois donc que le test if(abs(ogcm_var3d(i,j,max_z))<1.e9)
! est non seulement inutile mais penalisant. Je propose donc de le
! commenter:

         ksecu=1
!        if(abs(ogcm_var3d(i,j,max_z))<1.e9)then !fffffffff>
! Pourquoi k4=k+1? Parce que la lecture du fichier commence par le haut
! et renseignera au fur et a mesure ogcm_var3d des niveaux superieurs...
          i4=i ; j4=j ; k4=k+1    ! claude !17-07-14
          goto 2632               ! claude !17-07-14
!        endif                                   !fffffffff>

        endif                                 !---------->

        if(ksecu.eq.0)then !kkkkkkk>
           i10=i10+1
           goto 1864
        endif              !kkkkkkk>

 2632   continue
                    write(4,'(6(i4,1x))')i,j,k,i4,j4,k4
        if(k==max_z)write(3,'(4(i4,1x))')i,j,i4,j4

       endif                           !eeeeeeeeeeeeeee>

 1862 continue ! fin de boucle k
      close(4)
      close(3)

      end subroutine ogcm_fichier_bouchetrou_ts

!------------------------------------------------------------------------------

      subroutine ogcm_appli_bouchetrou_ts                          !19-11-10
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
      integer loop_
#ifdef synopsis
       subroutinetitle='ogcm_appli_bouchetrou_ts'
       subroutinedescription= &
          'Reads "bouchetrou" file and replaces land values'    &
       //' concerned by the interpolation process with nearest' &
       //' OGCM sea values. 3D arrays.'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0)write(6,*)'bouchage des trous'

      open(unit=3,file=trim(tmpdirname)//'bouchetrou_ts'//dom_c//'.out')   !#MPI

      do loop_=1,max_z*max_y*max_x                        !17-07-14
       read(3,*,end=1760)i,j,k,i4,j4,k4
         if(k==max_z) then !>>>>>>>                   ! Cet algo est possible car le fichier
              ogcm_var3d(i,j,k)=ogcm_var3d(i4,j4,k4)  ! bouchetrou est rempli depuis k=max_z
         else              !>>>>>>>

! Algo 1: melange d'extrapolation horizontale et verticale
!         x1=max(min( (ogcm_z(i,j,k+1)-ogcm_z(i,j,k))/200.,un),zero) ! Si DZ>200m 100% extrapolation horizontale

!         x1=1.                                                !29-01-15

!             ogcm_var3d(i,j,k)                              & !25-11-14
!            =ogcm_var3d(i,j,k+1)*(1.-x1)                    & !extrapolation verticale
!            +ogcm_var3d(i4,j4,k4 )  *x1                       !extrapolation horizontale

! Algo 2: Le gradient vertical en (i,j) est donné par le gradient       !28-01-15
! vertical en (i4,j4): on espere ainsi limiter le risque de creer une
! un gradient vertical de densite incoherent
!        if(k4==max_z)stop 'ogcm_appli_bouchetrou_ts k4+1>max_z'

!        ogcm_var3d(i ,j ,k )=ogcm_var3d(i,j,k+1)              &
!                   -(ogcm_z(i ,j ,k+1 )    -ogcm_z(i ,j ,k )) &
!                   /(ogcm_z(i4,j4,k4+1)    -ogcm_z(i4,j4,k4)) &
!               *(ogcm_var3d(i4,j4,k4+1)-ogcm_var3d(i4,j4,k4))

! Algo 3: le trou est bouche avec la valeur horizontale la plus proche !05-11-15
! mais on assure le recollement vertical avec le premier niveau au
! dessus du fond dont la position est donnee par k2:
         if(k4==k.and.k/=max_z) then !---> !14-11-15

! rap=0 si Dz=0, rap=1 si Dz>100m
          rap=min((ogcm_z(i,j,k+1)-ogcm_z(i,j,k))/100.,1.) 

          ogcm_var3d(i,j,k)=              &
               rap*  ogcm_var3d(i4,j4,k4) & ! bouchage simple
          +(1.-rap)*(ogcm_var3d(i4,j4,k4)+ogcm_var3d(i,j,k+1)-ogcm_var3d(i4,j4,k+1))
          
         else                     !--->

          ogcm_var3d(i,j,k)=ogcm_var3d(i4,j4,k4)  

         endif          !--->

! Algo 4: pour amoindrir les discontinuites du bouchage on moyenne les
! points non masques a distance similaire
!       if(k==k4) then !ppppppp>
!        i0=max( abs(i4-i) , abs(j4-j) )
!        sum1=0.  ; sum2=0.
!        do j1=j-i0,j+i0,2*i0
!        do i1=i-i0,i+i0
!         if(i1>=1.and.i1<=max_x.and.j1>=1.and.j1<=max_y) then
!          if(abs(ogcm_var3d(i1,j1,k4))<1.e10) then
!           sum1=sum1+1.
!           sum2=sum2+ogcm_var3d(i1,j1,k4)
!          endif          
!         endif          
!        enddo
!        enddo
!        do j1=j-i0+1,j+i0-1
!        do i1=i-i0  ,i+i0,2*i0
!         if(i1>=1.and.i1<=max_x.and.j1>=1.and.j1<=max_y) then
!          if(abs(ogcm_var3d(i1,j1,k4))<1.e10) then
!           sum1=sum1+1.
!           sum2=sum2+ogcm_var3d(i1,j1,k4)
!          endif          
!         endif          
!        enddo
!        enddo
!        ogcm_var3d(i,j,k)=sum2/sum1
!       else           !ppppppp>
!        ogcm_var3d(i,j,k)=ogcm_var3d(i4,j4,k4)  
!       endif          !ppppppp>
         

         endif             !>>>>>>>
      enddo
 1760 continue
      close(3)

      end subroutine ogcm_appli_bouchetrou_ts

!------------------------------------------------------------------------------

      subroutine ogcm_get_z_sympa(nc_   ,txt_   ) !23-11-11
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
      include 'netcdf.inc'
      character txt_   *1
      integer nc_,ncid_,statusprev_
#ifdef synopsis
       subroutinetitle='ogcm_get_z_sympa'
       subroutinedescription= &
       'Computes the depth of the levels of the ogcm in S case'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      ksecu=0

! Compute the depth of the levels of the ogcm in the "sympa" case
! using dz and the bathymetry of the "sympa" grid

! Ouvrir le fichier de grille
        open(unit=4,file=liste_ogcm_grid  &  !13-10-10
                   ,access='direct'       &
                   ,recl=541              &  !28-04-16
                   ,form='unformatted')
        read(4,rec=1)texte250
        close(4)
        status=nf_open(trim(texte250),nf_nowrite,ncid_   )
        if(status/=0)stop 'erreur ouverture fichier grille sympa'

!       call ogcm_get_dim(ncid_   ,'t') ! get max_x max_y max_z
        call ogcm_zoom_coordinate('t',1)

                     status=nf_inq_varid(ncid_   ,'h_t',var_id)
        if(status/=0)status=nf_inq_varid(ncid_   ,'h_w',var_id) !16-02-12
        if(status/=0)stop 'erreur identifiant h_t'
        status=nf_inq_varndims(ncid_,var_id,ogcm_var_ndims)
        if(status/=0)stop ' erreur 3442'

! Lire la bathy de l'ogcm sympa
        varstart(1)=ogcmzoom_istrt ; varcount(1)=max_x
        varstart(2)=ogcmzoom_jstrt ; varcount(2)=max_y
        status=nf_get_vara_real(ncid_,var_id                 &
                               ,varstart(1:ogcm_var_ndims)   &
                               ,varcount(1:ogcm_var_ndims)   &
                               ,ogcm_var2d(1:max_x,1:max_y))

        if(status/=0)stop 'erreur lecture h_t ou h_w'
! Fermer le fichier de grille:
        status=nf_close(ncid_   ) ! fermer le fichier de grille

! Ouvrir le fichier de variable:
        open(unit=4,file=liste_ogcm_uvts      &  !13-10-10
                   ,access='direct'           &
                   ,recl=541                  &  !28-04-16
                   ,form='unformatted')
        read(4,rec=nc_   )texte250
        close(4)
        status=nf_open(trim(texte250),nf_nowrite,ncid_   ) ! re-ouvrir le fichier de variable
        if(status/=0)stop 'erreur ouverture fichier sympa'

! Lire dz (au point t)
      status=nf_inq_varid(ncid_   ,'dz',var_id)

      statusprev_=status
      if(status==0)then !00000000>
      write(6,*)'Ne pas deduire dephu d''une demi somme'              &
               ,' de deptht pour plusieurs raisons: '                 &
               ,' 1 la grille peut etre basee sur points vorticite '  &
               ,' 2 les bornes istrt et istru ne permettent pas de '  &
               ,'savoir qui precede qui'
      stop ' stop donc 3563'
       status=nf_inq_varndims(ncid_,var_id,ogcm_var_ndims)
       if(status/=0)stop ' erreur 3470'
       varstart(1)=ogcmzoom_istrt ; varcount(1)=max_x
       varstart(2)=ogcmzoom_jstrt ; varcount(2)=max_y
       varstart(3)=1              ; varcount(3)=max_z
       status=nf_get_vara_real(ncid_,var_id                 &
                              ,varstart(1:ogcm_var_ndims)   &
                              ,varcount(1:ogcm_var_ndims)   &
                              ,ogcm_var3d(1:max_x,1:max_y,1:max_z))
       if(status/=0)stop ' erreur 3479'
! Passer de dz e z: z=-h_t+integrale(dz)
       do j=1,max_y
       do i=1,max_x
        sum1=-ogcm_var2d(i,j)
        do k=1,max_z
         sum1=sum1+ogcm_var3d(i,j,k)
         ogcm_z(i,j,k)=sum1-0.5*ogcm_var3d(i,j,k)
        enddo
       enddo
       enddo
      endif             !00000000>

      status=nf_close(ncid_   ) ! fermer le fichier de variable

! Si status=0 cela signifie que dz n'est pas contenu dans le fichier de variable.
! A la place, on lit z dans le fichier (unique) de grille. C'est pas top car les
! niveaux peuvent bouger beaucoup (grille ALE) mais c'est mieux que rien.
      if(statusprev_/=0)then !-------->

! Ouvrir le fichier de grille:
        open(unit=4,file=liste_ogcm_grid    &  !13-10-10
                   ,access='direct'         &
                   ,recl=541                &  !28-04-16
                   ,form='unformatted')
        read(4,rec=1)texte250
        close(4)
        if(par%rank==0)write(6,'(a)')trim(texte250)

        status=nf_open(trim(texte250),nf_nowrite,ncid_   )
        if(status/=0)stop 'erreur 879'

        if(txt_=='t')status=nf_inq_varid(ncid_,'depth_t',var_id)
        if(txt_=='u')status=nf_inq_varid(ncid_,'depth_u',var_id)
        if(txt_=='v')status=nf_inq_varid(ncid_,'depth_v',var_id)
        if(status/=0)stop 'erreur 881'
        status=nf_inq_varndims(ncid_,var_id,ogcm_var_ndims)
        if(status/=0)stop ' erreur 3516'

        if(txt_=='t') then !tttttt>
         varstart(1)=ogcmzoom_istrt ; varstart(2)=ogcmzoom_jstrt
        endif              !tttttt>
        if(txt_=='u') then !uuuuuu>
         varstart(1)=ogcmzoom_istru ; varstart(2)=ogcmzoom_jstru
        endif              !uuuuuu>
        if(txt_=='v') then !vvvvvv>
         varstart(1)=ogcmzoom_istrv ; varstart(2)=ogcmzoom_jstrv
        endif              !vvvvvv>
        varcount(1)=max_x ; varcount(2)=max_y
        varstart(3)=1     ; varcount(3)=max_z

        status=nf_inq_vartype(ncid_,var_id,k0)

        if(k0==nf_real) then !------->
         status=nf_get_vara_real(ncid_,var_id                 &
                                ,varstart(1:ogcm_var_ndims)   &
                                ,varcount(1:ogcm_var_ndims)   &
                                ,ogcm_z(1:max_x,1:max_y,1:max_z))
        else                 !------->
         stop ' E3502 ogcm_z must be real'    !23-06-14
        endif                !------->
        if(status/=0) then !>>>>>>
         write(6,*)'E3501 par%rank=',par%rank                   & !23-06-14
                  ,'max_x max_y max_z ',max_x,max_y,max_z       &
                  ,'varstart(1:3) ',varstart(1:3)               &
                  ,'varcount(1:3) ',varcount(1:3)               &
                  ,'ogcm_var_ndims=',ogcm_var_ndims             &
                  ,'vartype nf_real ',k0,nf_real                &
                  ,'vstar+vcount-1  ',varstart(1)+varcount(1)-1 &
                                     ,varstart(2)+varcount(2)-1 &
                  ,'txt_=',txt_
         stop ' erreur 3501'
        endif              !>>>>>>
! Fermer le fichier de grille:
        status=nf_close(ncid_   )
        if(status/=0)stop ' erreur 3504'

      endif                  !------->

      if(txt_   =='t')ksecu=1 ! temoin de passage
      if(txt_   =='u')ksecu=1 ! temoin de passage
      if(txt_   =='v')ksecu=1 ! temoin de passage
! Verifier le temoin de passage avant de sortir:
      if(ksecu==0)stop 'get_ogcm_z_sympa bad argument'

      end subroutine ogcm_get_z_sympa
!..........................................................................
      subroutine ogcm_error_diag(t_or_s_)
      use module_principal
      use module_parallele
      use module_forcages
      implicit none
      integer k_,t_or_s_
#ifdef synopsis
       subroutinetitle='ogcm_error_diag'
       subroutinedescription='Gives error diagnostic if any'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(par%rank==0)write(6,*)'//////////////////////////////////////////////'
      if(t_or_s_==1.and.par%rank==0)write(6,*)'Pas de valeurs pour interpoler T'
      if(t_or_s_==2.and.par%rank==0)write(6,*)'Pas de valeurs pour interpoler S'
      if(par%rank==0)write(6,*)'sum1=',sum1
      if(par%rank==0)write(6,*)'sum2=',sum2
      if(par%rank==0)write(6,*)'par%rank i j k',par%rank,i,j,k
      if(par%rank==0)write(6,'(a,a)')'dom_c=',dom_c
      if(par%rank==0)write(6,*)'depth_w(k:k+1)',depth_w(i,j,k:k+1)
      if(par%rank==0)write(6,*)'ogcm_z1d min max=',ogcm_z1d(1),ogcm_z1d(max_z)
      if(par%rank==0)write(6,'(a,5i4)')'Indices de correspondances sous le fond=',i0,j0,i2,j2,k10
      if(par%rank==0)write(6,*)'Bathy ogcm en i2 j2:',ogcm_var2d(i2,j2)
      if(par%rank==0)write(6,*)'ogcm_z(i2,j2,1:max_z)=------>'
      do k_=max_z,1,-1
      if(par%rank==0)write(6,*)k_,ogcm_z(i2,j2,k_)
      enddo


      stop 'STOP dans error_diag'

      end subroutine ogcm_error_diag

!.................................................................................

      subroutine ogcm_allocate(ki1_,ki2_,dim1_,dim2_,dim3_)
      use module_forcages
      implicit none
      integer ki1_,ki2_,dim1_,dim2_,dim3_
#ifdef synopsis
       subroutinetitle='ogcm_allocate'
       subroutinedescription='Arrays allocation for OGCM interpolation'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


       if(ki2_==1) then !>>>>>>>>>>>>>>>>>>>>>>>

        if(ki1_==1) then !111111111111111>

         if(allocated( ogcm_var3d ))deallocate(ogcm_var3d ) !16-05-19
         if(allocated( short3d    ))deallocate(short3d    )
         if(allocated( ogcm_var2d ))deallocate(ogcm_var2d )
         if(allocated( ogcm_lat   ))deallocate(ogcm_lat   )
         if(allocated( ogcm_lon   ))deallocate(ogcm_lon   )
         if(allocated( short2d    ))deallocate(short2d    )
         if(allocated( ogcm_z     ))deallocate(ogcm_z     )
         if(allocated( ogcm_var1d ))deallocate(ogcm_var1d )
         if(allocated( ogcm_z1d   ))deallocate(ogcm_z1d   )

         allocate(ogcm_var3d    (dim1_,dim2_,dim3_)) ; ogcm_var3d=0 
         allocate(short3d       (dim1_,dim2_,dim3_)) ; short3d=0
         allocate(ogcm_var2d    (dim1_,dim2_      )) ; ogcm_var2d=0
         allocate(ogcm_lat      (dim1_,dim2_      )) ; ogcm_lat=0
         allocate(ogcm_lon      (dim1_,dim2_      )) ; ogcm_lon=0
         allocate(short2d       (dim1_,dim2_      )) ; short2d=0
         allocate(ogcm_z        (dim1_,dim2_,dim3_)) ; ogcm_z=0
         allocate(ogcm_var1d    (            dim3_)) ; ogcm_var1d=0
         allocate(ogcm_z1d      (            dim3_)) ; ogcm_z1d=0

        endif            !111111111111111>

        if(ki1_==2) then !222222222222222>

         if(.not.allocated( ogcm_var3d ))stop 'ERREUR 3606 A'
         if(.not.allocated( short3d    ))stop 'ERREUR 3606 B'
         if(.not.allocated( ogcm_var2d ))stop 'ERREUR 3606 C'
         if(.not.allocated( ogcm_lat   ))stop 'ERREUR 3606 D'
         if(.not.allocated( ogcm_lon   ))stop 'ERREUR 3606 E'
         if(.not.allocated( short2d    ))stop 'ERREUR 3606 F'
         if(.not.allocated( ogcm_z     ))stop 'ERREUR 3606 G'
         if(.not.allocated( ogcm_var1d ))stop 'ERREUR 3606 I'
         if(.not.allocated( ogcm_z1d   ))stop 'ERREUR 3606 J'

         deallocate(ogcm_var3d)
         deallocate(ogcm_var2d)
         deallocate(ogcm_lat)
         deallocate(ogcm_lon)
         deallocate(ogcm_z  )
         deallocate(ogcm_var1d)
         deallocate(ogcm_z1d)
         deallocate(short2d)
         deallocate(short3d)

        endif            !222222222222222>

       endif            !>>>>>>>>>>>>>>>>>>>>>>>

      end subroutine ogcm_allocate

!.........................................................................

      subroutine ogcm_get_lonlat(ncid_,gridlocation_)       !04-04-12
      use module_principal
      use module_forcages
      implicit none
      include 'netcdf.inc'
      integer type_,ncid_
      character*1 gridlocation_
      real*4,dimension(:,:),allocatable ::                             &
             ogcm2d_r4_
      real*4,dimension(:),allocatable ::                               &
             ogcm1d_r4_
      real*8,dimension(:),allocatable ::                               &
             ogcm1d_r8_
#ifdef synopsis
       subroutinetitle='ogcm_get_lonlat'
       subroutinedescription='Reads lon lat in OGCM netcdf file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Load ogcm longitude:
       if(gridlocation_=='t') then !tttttttttt>
                    status=nf_inq_varid(ncid_,'glamt',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lon_zhl',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lon_ZHL',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lon_t',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'longitude_t',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'nav_lon',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'longitude',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lon',var_id) !16-05-18
       endif                       !tttttttttt>
       if(gridlocation_=='u') then !uuuuuuuuuu>
                    status=nf_inq_varid(ncid_,'glamu',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lon_xhl',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lon_XHL',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lon_u',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'longitude_u',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'nav_lon',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'longitude',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lon',var_id) !16-05-18
       endif                       !uuuuuuuuuu>
       if(gridlocation_=='v') then !vvvvvvvvvv>
                    status=nf_inq_varid(ncid_,'glamv',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lon_yhl',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lon_YHL',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lon_v',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'longitude_v',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'nav_lon',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'longitude',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lon',var_id) !16-05-18
       endif                       !vvvvvvvvvv>
       if(status/=0)stop 'interp_ogcm erreur var_id longitude'

       status=nf_inq_vartype(ncid_,var_id,type_)
       if(type_/=nf_float.and.type_/=nf_double)     &
       stop 'ogcm longitude type not accepted'

       status=nf_inq_varndims(ncid_,var_id,ogcm_var_ndims)
       if(status/=0)stop ' erreur nf_inq_varndim ogcm_var_ndims'
       varstart(1)=ogcmzoom_istr ; varcount(1)=max_x
       varstart(2)=ogcmzoom_jstr ; varcount(2)=max_y
       varstart(3)=1             ; varcount(3)=1

       if (ogcm_var_ndims==1) then

         if(type_==nf_float) then  !----------->
          allocate(ogcm1d_r4_(max_x))
          status=nf_get_vara_real(ncid_,var_id          &
                                       ,varstart(1:1)   & !25-11-16
                                       ,varcount(1:1)   &
                                       ,ogcm1d_r4_(1:max_x))
          do j=1,max_y ; ogcm_lon(:,j)=ogcm1d_r4_(:) ; enddo
          if(status/=0)stop 'error nf_get_vara_real ogcm2d_r4_ lon'
          deallocate(ogcm1d_r4_)
         endif                     !----------->                     
         if(type_==nf_double) then !----------->
          allocate(ogcm1d_r8_(max_x))
          status=nf_get_vara_double(ncid_,var_id                     &
                                         ,varstart(1:1) &
                                         ,varcount(1:1) &
                                         ,ogcm1d_r8_(1:max_x))
          do j=1,max_y ; ogcm_lon(:,j)=ogcm1d_r8_(:) ; enddo
          if(status/=0)stop 'error nf_get_vara_real ogcm_lon'
          deallocate(ogcm1d_r8_)
         endif                     !----------->

       else

         if(type_==nf_float) then  !----------->
          allocate(ogcm2d_r4_(max_x,max_y))
          status=nf_get_vara_real(ncid_,var_id                       &
                                       ,varstart(1:ogcm_var_ndims)   &
                                       ,varcount(1:ogcm_var_ndims)   &
                                       ,ogcm2d_r4_(1:max_x,1:max_y))
          ogcm_lon(1:max_x,1:max_y)=ogcm2d_r4_(1:max_x,1:max_y)
          if(status/=0)stop 'error nf_get_vara_real ogcm2d_r4_ lon'
          deallocate(ogcm2d_r4_)
         endif                     !----------->
         if(type_==nf_double) then !----------->
          status=nf_get_vara_double(ncid_,var_id                     &
                                         ,varstart(1:ogcm_var_ndims) &
                                         ,varcount(1:ogcm_var_ndims) &
                                         ,ogcm_lon(1:max_x,1:max_y))
          if(status/=0)stop 'error nf_get_vara_real ogcm_lon'
         endif                     !----------->

       endif                     !----------->
       if(status/=0)stop 'can not get ogcm longitude'


! Load ogcm latitude:
       if(gridlocation_=='t') then !tttttttttt>
                    status=nf_inq_varid(ncid_,'gphit',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lat_zhl',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lat_ZHL',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lat_t',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'latitude_t',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'nav_lat',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'latitude',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lat',var_id) !16-05-18
       endif                       !tttttttttt>
       if(gridlocation_=='u') then !uuuuuuuuuu>
                    status=nf_inq_varid(ncid_,'gphiu',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lat_xhl',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lat_XHL',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lat_u',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'latitude_u',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'nav_lat',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'latitude',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lat',var_id) !16-05-18
       endif                       !uuuuuuuuuu>
       if(gridlocation_=='v') then !vvvvvvvvvv>
                    status=nf_inq_varid(ncid_,'gphiv',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lat_yhl',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lat_YHL',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lat_v',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'latitude_v',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'nav_lat',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'latitude',var_id)
       if(status/=0)status=nf_inq_varid(ncid_,'lat',var_id) !16-05-18
       endif                       !vvvvvvvvvv>
       if(status/=0)stop 'interp_ogcm erreur var_id latitude'

       status=nf_inq_vartype(ncid_,var_id,type_)
       if(type_/=nf_float.and.type_/=nf_double)     &
       stop 'ogcm latitude type not accepted'

       status=nf_inq_varndims(ncid_,var_id,ogcm_var_ndims)
       if(status/=0)stop ' erreur nf_inq_varndim ogcm_var_ndims'
       varstart(1)=ogcmzoom_istr ; varcount(1)=max_x
       varstart(2)=ogcmzoom_jstr ; varcount(2)=max_y
       varstart(3)=1             ; varcount(3)=1

       if (ogcm_var_ndims==1) then

         if(type_==nf_float) then  !----------->
          allocate(ogcm1d_r4_(max_y))
          status=nf_get_vara_real(ncid_,var_id         &
                                       ,varstart(2:2)  &
                                       ,varcount(2:2)  &
                                       ,ogcm1d_r4_(1:max_y))
          do i=1,max_x ; ogcm_lat(i,:)=ogcm1d_r4_(:) ; enddo
          if(status/=0)stop 'error nf_get_var_real ogcm2d_r4_ lat'
          deallocate(ogcm1d_r4_)
         endif                     !----------->
         if(type_==nf_double) then !----------->
          allocate(ogcm1d_r8_(max_y))
           status=nf_get_vara_double(ncid_,var_id         &
                                          ,varstart(2:2)  &
                                          ,varcount(2:2)  &
                                          ,ogcm1d_r8_(1:max_y))
          do i=1,max_x ; ogcm_lat(i,:)=ogcm1d_r8_(:) ; enddo
          if(status/=0)stop 'error nf_get_var_real ogcm_lat'
          deallocate(ogcm1d_r4_)
         endif 

       else

         if(type_==nf_float) then  !----------->
          allocate(ogcm2d_r4_(max_x,max_y))
          status=nf_get_vara_real(ncid_,var_id                      &
                                       ,varstart(1:ogcm_var_ndims)  &
                                       ,varcount(1:ogcm_var_ndims)  &
                                       ,ogcm2d_r4_(1:max_x,1:max_y))
          ogcm_lat(1:max_x,1:max_y)=ogcm2d_r4_(1:max_x,1:max_y)
          if(status/=0)stop 'error nf_get_var_real ogcm2d_r4_ lat'
          deallocate(ogcm2d_r4_)
         endif                     !----------->
         if(type_==nf_double) then !----------->
           status=nf_get_vara_double(ncid_,var_id                      &
                                          ,varstart(1:ogcm_var_ndims)  &
                                          ,varcount(1:ogcm_var_ndims)  &
                                          ,ogcm_lat(1:max_x,1:max_y))
          if(status/=0)stop 'error nf_get_var_real ogcm_lat'
         endif                     !----------->

       endif                     !----------->
       if(status/=0)stop 'can not get ogcm latitude'

      end subroutine ogcm_get_lonlat

!.........................................................................

      subroutine ogcm_zoom_coordinate(txt_,case_)
      implicit none
      character*1 txt_
      integer case_
#ifdef synopsis
       subroutinetitle='ogcm_zoom_coordinate'
       subroutinedescription='Gets the bounds of the reading zone'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(case_==1) then !ooooooooooooooooooooooo>

       if(txt_=='t') then !----->
        ogcmzoom_istr=ogcmzoom_istrt
        ogcmzoom_iend=ogcmzoom_iendt
        ogcmzoom_jstr=ogcmzoom_jstrt
        ogcmzoom_jend=ogcmzoom_jendt
       endif              !----->
       if(txt_=='u') then !----->
        ogcmzoom_istr=ogcmzoom_istru
        ogcmzoom_iend=ogcmzoom_iendu
        ogcmzoom_jstr=ogcmzoom_jstru
        ogcmzoom_jend=ogcmzoom_jendu
       endif              !----->
       if(txt_=='v') then !----->
        ogcmzoom_istr=ogcmzoom_istrv
        ogcmzoom_iend=ogcmzoom_iendv
        ogcmzoom_jstr=ogcmzoom_jstrv
        ogcmzoom_jend=ogcmzoom_jendv
       endif              !----->

       if(txt_/='m') then !----->
        max_x=ogcmzoom_iend-ogcmzoom_istr+1
        max_y=ogcmzoom_jend-ogcmzoom_jstr+1
       else               !------>
        max_x=    ogcmzoom_iendt-ogcmzoom_istrt+1
        max_y=    ogcmzoom_jendt-ogcmzoom_jstrt+1
        max_x=max(ogcmzoom_iendu-ogcmzoom_istru+1,max_x)
        max_y=max(ogcmzoom_jendu-ogcmzoom_jstru+1,max_y)
        max_x=max(ogcmzoom_iendv-ogcmzoom_istrv+1,max_x)
        max_y=max(ogcmzoom_jendv-ogcmzoom_jstrv+1,max_y)
       endif

      endif                   !ooooooooooooooooooooooo>

      if(case_==0) then       !iiiiiiiiiiiiiiiiiiiiiii>

! A partir des min & max des (lon,lat) des 2 grilles, reduire la zone d'extraction:
       ogcmzoom_istr=999999 ; ogcmzoom_iend=-999999 ! first guess
       ogcmzoom_jstr=999999 ; ogcmzoom_jend=-999999 ! first guess

       ksecu=0
       do j=1,max_y
       do i=1,max_x

       if(ogcm_lon(i,j)>=lonmin.and.ogcm_lon(i,j)<=lonmax.and.   &
          ogcm_lat(i,j)>=latmin.and.ogcm_lat(i,j)<=latmax)then

         ogcmzoom_istr=min(ogcmzoom_istr,i)
         ogcmzoom_jstr=min(ogcmzoom_jstr,j)
         ogcmzoom_iend=max(ogcmzoom_iend,i)
         ogcmzoom_jend=max(ogcmzoom_jend,j)
         ksecu=1

       endif

       enddo
       enddo

! Si ksecu=0 c'est que le proc est si petit, qu'aucun point ogcm ne s'est trouve e l'interieur
! On refait le test differement. On cherche le point le plus proche.
       if(ksecu==0)then !000000000000000>
        dist1=1.d20 ; i1=imax/2 ; j1=jmax/2
        do j=1,max_y
        do i=1,max_x

         dist2=rayonterre*                                         &
         acos( sin(ogcm_lat(i,j)*deg2rad)*sin(lat_t(i1,j1))     &
              +cos(ogcm_lat(i,j)*deg2rad)*cos(lat_t(i1,j1))     &
              *cos(lon_t(i1,j1)-ogcm_lon(i,j)*deg2rad))

         if(dist2<dist1)then !--->
          dist1=dist2 ; i2=i ; j2=j
         endif               !--->

        enddo
        enddo
        ogcmzoom_istr=i2 ; ogcmzoom_iend=i2
        ogcmzoom_jstr=j2 ; ogcmzoom_jend=j2
       endif            !000000000000000>


! on elargit un peu plus pour boucher les trous sans etre restreint par la taille
! reduite de la zone d'extraction, et puis aussi pour rattraper l'erreur liee au fait
! que plusieurs grilles ogcm (point u, v, t) peuvent etre presentes.
! elargissement a 10 lignes / 10 colonnes supplementaires - repere 1233 :
!      ogcm_enlarged_bounds est defini a 10 points en tete de module_ogcm

       ogcmzoom_istr=max(ogcmzoom_istr-ogcm_enlarged_bounds,1      )
       ogcmzoom_iend=min(ogcmzoom_iend+ogcm_enlarged_bounds,max_x-1) !23-06-14
       ogcmzoom_jstr=max(ogcmzoom_jstr-ogcm_enlarged_bounds,1      )
       ogcmzoom_jend=min(ogcmzoom_jend+ogcm_enlarged_bounds,max_y-1) !23-06-14

       if(txt_=='t') then !----->
        ogcmzoom_istrt=ogcmzoom_istr
        ogcmzoom_iendt=ogcmzoom_iend
        ogcmzoom_jstrt=ogcmzoom_jstr
        ogcmzoom_jendt=ogcmzoom_jend
       endif              !----->
       if(txt_=='u') then !----->
        ogcmzoom_istru=ogcmzoom_istr
        ogcmzoom_iendu=ogcmzoom_iend
        ogcmzoom_jstru=ogcmzoom_jstr
        ogcmzoom_jendu=ogcmzoom_jend
       endif              !----->
       if(txt_=='v') then !----->
        ogcmzoom_istrv=ogcmzoom_istr
        ogcmzoom_iendv=ogcmzoom_iend
        ogcmzoom_jstrv=ogcmzoom_jstr
        ogcmzoom_jendv=ogcmzoom_jend
       endif              !----->


!     max_x=ogcmzoom_iend-ogcmzoom_istr+1
!     max_y=ogcmzoom_jend-ogcmzoom_jstr+1

      endif                   !iiiiiiiiiiiiiiiiiiiiiii>


!     write(6,*)'zoom i ',par%rank,ogcmzoom_istr,ogcmzoom_iend
!     write(6,*)'zoom j ',par%rank,ogcmzoom_jstr,ogcmzoom_jend
!     write(6,*)'APRES',par%rank,max_x,max_y

      end subroutine ogcm_zoom_coordinate

!.............................................................................

      subroutine ogcm_get_gridfilename(txt_,recnum_) ! donne ogcm_grid_file_name (character)
      implicit none
      character*1 txt_
      integer recnum_ !16-05-19
#ifdef synopsis
       subroutinetitle='ogcm_get_gridfilename'
       subroutinedescription='Gets the name of the OGCM grid file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Obtenir le nom du fichier de grille:
      texte90='none'

       if(obc_ogcm_type(1:5)=='sympa') texte90=liste_ogcm_grid

      if(txt_=='t') then !ttttttttttt>
       if(obc_ogcm_type(1:4)=='ncom')  texte90=liste_ogcm_t
       if(obc_ogcm_type(1:6)=='nemo_z')texte90=liste_ogcm_grid_t
      endif              !ttttttttttt>
      if(txt_=='u') then !uuuuuuuuuuu>
       if(obc_ogcm_type(1:4)=='ncom')  texte90=liste_ogcm_u
       if(obc_ogcm_type(1:6)=='nemo_z')texte90=liste_ogcm_grid_u
      endif              !uuuuuuuuuuu>
      if(txt_=='v') then !vvvvvvvvvvv>
       if(obc_ogcm_type(1:4)=='ncom')  texte90=liste_ogcm_v
       if(obc_ogcm_type(1:6)=='nemo_z')texte90=liste_ogcm_grid_v
      endif              !vvvvvvvvvvv>
      if(texte90=='none')stop ' nom grd invalide dans module_ogcm'

      open(unit=4,file=trim(texte90) &
                 ,access='direct'    &
                 ,recl=541           & !28-04-16
                 ,form='unformatted')
      read(4,rec=recnum_)ogcm_grid_file_name 
      close(4)

      if(par%rank==0)write(6,'(a,a)')'Grid file=',trim(ogcm_grid_file_name)

      end subroutine ogcm_get_gridfilename

!.............................................................................

      subroutine ogcm_get_varfilename(var_id_,txt_,nc_) ! donne texte250
      implicit none
      character*1 txt_
      integer nc_,var_id_
#ifdef synopsis
       subroutinetitle='ogcm_get_varfilename'
       subroutinedescription='Gets the name of the variable netcdf file'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Obtenir le nom du fichier de variable:
      texte90='none'

      if(obc_ogcm_type(1:5)=='sympa') texte90=liste_ogcm_uvts

      if(txt_=='e') then !eeeeeeeeeee> ! 'e' pour elevation de la surface
       if(obc_ogcm_type(1:4)=='ncom')  texte90=liste_ogcm_ssh
       if(obc_ogcm_type(1:6)=='nemo_z')texte90=liste_ogcm_ssh
      endif              !eeeeeeeeeee>
      if(txt_=='s') then !sssssssssss>
       if(obc_ogcm_type(1:4)=='ncom')  texte90=liste_ogcm_s
       if(obc_ogcm_type(1:6)=='nemo_z')texte90=liste_ogcm_s
      endif              !sssssssssss>
      if(txt_=='t') then !ttttttttttt>
       if(obc_ogcm_type(1:4)=='ncom')  texte90=liste_ogcm_t
       if(obc_ogcm_type(1:6)=='nemo_z')texte90=liste_ogcm_t
      endif              !ttttttttttt>
      if(txt_=='u') then !uuuuuuuuuuu>
       if(obc_ogcm_type(1:4)=='ncom')  texte90=liste_ogcm_u
       if(obc_ogcm_type(1:6)=='nemo_z')texte90=liste_ogcm_u
      endif              !uuuuuuuuuuu>
      if(txt_=='v') then !vvvvvvvvvvv>
       if(obc_ogcm_type(1:4)=='ncom')  texte90=liste_ogcm_v
       if(obc_ogcm_type(1:6)=='nemo_z')texte90=liste_ogcm_v
      endif              !vvvvvvvvvvv>
      if(texte90=='none')stop ' nom var invalide dans module_ogcm'

!-----------------------------------------------------------------------------------------!
! Seul par%rank=0 lit le fichier puis envoie les valeurs aux autres par%rank par mpi_bcast
      if(par%rank==0) then !--[m0v0m]--> !05-05-18
       open(unit=4,file=trim(texte90) &
                  ,access='direct'    &
                  ,recl=541           & !28-04-16
                  ,form='unformatted')
       read(4,rec=nc_   )texte250,ogcm_time_counter           &
                                 ,ogcm_readtime_next(var_id_) &
                                   ,ogcm_period_next(var_id_) &
                                   ,ogcmgridcode !16-05-19
       close(4)
       write(6,'(a,a)')'Variable file=',trim(texte250)
      endif                !--[m0v0m]--> !05-05-18

      k=len(texte250)
      call mpi_bcast(texte250,k,mpi_char,0,par%comm2d,ierr)

      call mpi_bcast(ogcm_time_counter,1,mpi_integer,0,par%comm2d,ierr)

      x0=ogcm_readtime_next(var_id_)
      call mpi_bcast(x0,1,mpi_double_precision,0,par%comm2d,ierr)
      ogcm_readtime_next(var_id_)=x0

      x0=ogcm_period_next(var_id_)
      call mpi_bcast(x0,1,mpi_double_precision,0,par%comm2d,ierr)
      ogcm_period_next(var_id_)=x0


!-----------------------------------------------------------------------------------------!

      end subroutine ogcm_get_varfilename

!------------------------------------------------------------------------------

      subroutine ogcm_interp_core1(ichoix,t_,nc_)
      implicit none
      double precision zup_,zdw_,vup_,vdw_,dzsurf_
      integer ichoix,t_,nc_,ncid_,valmin_,valmax_
#ifdef synopsis
       subroutinetitle='ogcm_interp_core1'
       subroutinedescription= &
          'Reads T or S field in the netcdf OGCM file and interpolates' &
       //' on the S grid and stores in temobc_t ans salobc_t'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif




      open(unit=3 &
             ,file=trim(tmpdirname)//'hr_to_lr_coord_1_'//dom_c//'.out')
      open(unit=4,file=trim(tmpdirname)//'sous_le_fond_'//dom_c//'.out')

! Interpolation:
      do j=0,jmax+1
      do i=0,imax+1

!     read(3,*)i0,j0,deci,decj,x0,rapi,rapj
      read(3,*)i0,j0,i1,j1,x0,rapi,rapj
      if(i0.ne.i)stop 'erreur1 sur i'
      if(j0.ne.j)stop 'erreur1 sur j'

      if(  mask_t(i,j,kmaxp1).eq.1) then !>>>>>>>>>>>>

      read(4,*)i0,j0,i2,j2,k10
! charger le profil de reference pour extrapolation sous le fond de l'ogcm:
      if(i0.ne.i)stop 'erreur2 sur i'
      if(j0.ne.j)stop 'erreur2 sur j'

      if(k10>=0)then !exexexexexex>
        do k1=1,max_z
           ogcm_z1d(k1)=    ogcm_z(i2,j2,k1)
         ogcm_var1d(k1)=ogcm_var3d(i2,j2,k1)
        enddo
      endif         !exexexexexex>

! Si pas de point assez profond dans l'ogcm alors extrapolation e gradient nul:
      if(k10==-1)then !trop_profond>
        do k1=1,max_z
           ogcm_z1d(k1)=-real(max_z-k1)/real(max_z-1)*hmax
         ogcm_var1d(k1)=0.
        enddo
      endif           !trop_profond>

!      i1=int(deci)
!      rapi=deci-i1
!      j1=int(decj)
!      rapj=decj-j1

       k0=kmaxp1
       ksecu=0
       const2=0. ! constante d'ajustement du profil de reference et de la grille grossiere
       do k=kmaxp1-1,kmin_w(i,j),-1

        sum1=0.
        sum2=0.

        do k1=max_z-1,1,-1

! Z1 & Z2  profondeurs dans la grille grossiere en K1 & K1+1
        z1= (1.-rapi)*(1.-rapj)*ogcm_z(i1  ,j1  ,k1)                    &
           +(1.-rapi)*    rapj *ogcm_z(i1  ,j1+1,k1)                    &
           +    rapi *(1.-rapj)*ogcm_z(i1+1,j1  ,k1)                    &
           +    rapi *    rapj *ogcm_z(i1+1,j1+1,k1)
        z2= (1.-rapi)*(1.-rapj)*ogcm_z(i1  ,j1  ,k1+1)                  &
           +(1.-rapi)*    rapj *ogcm_z(i1  ,j1+1,k1+1)                  &
           +    rapi *(1.-rapj)*ogcm_z(i1+1,j1  ,k1+1)                  &
           +    rapi *    rapj *ogcm_z(i1+1,j1+1,k1+1)
         if(k1==max_z-1) then !ssss> ! pour mettre les surfaces des 2 modeles au meme niveau:
           dzsurf_=z2-depth_t(i,j,kmax  ) ! 25-02-13
         endif                !ssss>

! Pour s'assurer que les niveaux de surface des 2 codes coincident bien,
! c'est z-ssh que l'on considere et pas z:
         zup_   =min( depth_w(i,j,k+1)+dzsurf_*real(k+1-kmin_w(i,j))/real(kmaxp1-kmin_w(i,j)),z2) !25-02-13
         zdw_   =max( depth_w(i,j,k  )+dzsurf_*real(k  -kmin_w(i,j))/real(kmaxp1-kmin_w(i,j)),z1) !25-02-13

         if(zup_   -zdw_   .gt.zero) then !-_-_-_-_-_-_>

           x2=(1.-rapi)*(1.-rapj)*ogcm_var3d(i1  ,j1  ,k1+1)            &
             +(1.-rapi)*    rapj *ogcm_var3d(i1  ,j1+1,k1+1)            &
             +    rapi *(1.-rapj)*ogcm_var3d(i1+1,j1  ,k1+1)            &
             +    rapi *    rapj *ogcm_var3d(i1+1,j1+1,k1+1)
           x1=(1.-rapi)*(1.-rapj)*ogcm_var3d(i1  ,j1  ,k1  )            &
             +(1.-rapi)*    rapj *ogcm_var3d(i1  ,j1+1,k1  )            &
             +    rapi *(1.-rapj)*ogcm_var3d(i1+1,j1  ,k1  )            &
             +    rapi *    rapj *ogcm_var3d(i1+1,j1+1,k1  )


           rapk=(zdw_   -z1)/(z2-z1)

           vdw_   =(1.-rapk)*x1+rapk*x2

           rapk=(zup_   -z1)/(z2-z1)

           vup_   =(1.-rapk)*x1+rapk*x2

          sum1=sum1+(zup_   -zdw_   )
          sum2=sum2+(zup_   -zdw_   )*0.5*(vup_   +vdw_   )

!         if(i+par%timax(1)==383.and.j+par%tjmax(1)==411.and.ichoix==4) then
!          write(64,*)'-----------------'
!          write(64,*)'k=',k
!          write(64,*)'Zs26',depth_w(i,j,k),depth_w(i,j,k+1)
!          write(64,*)'ZNem',z1,z2
!          write(64,*)'vup_,vdw_',vup_,vdw_
!          write(64,*)'i1 j1 k1',i1,j1,k1
!          write(64,*)'ogcm_var3d(i1  ,j1  ,k1  )',ogcm_var3d(i1  ,j1  ,k1  )
!          write(64,*)'ogcm_var3d(i1+1,j1  ,k1  )',ogcm_var3d(i1+1,j1  ,k1  )
!          write(64,*)'ogcm_var3d(i1  ,j1+1,k1  )',ogcm_var3d(i1  ,j1+1,k1  )
!          write(64,*)'ogcm_var3d(i1+1,j1+1,k1  )',ogcm_var3d(i1+1,j1+1,k1  )
!          write(64,*)'dom_c',dom_c
!          write(64,*)'par%rank',par%rank
!         endif

         endif          !-_-_-_-_-_-_>

        enddo ! fin K1

       if(sum1>0.001) then !sssssssssss>
! Cas normal oe la colonne de la grille grossiere englobe la grille fine locale:

        k0=k
        if(ichoix==1)temobc_t(i,j,k,t_)=sum2/max(sum1,small2)
        if(ichoix==2)salobc_t(i,j,k,t_)=sum2/max(sum1,small2)
#ifdef bidonref
        if(ichoix==3)temref_t(i,j,k   )=sum2/max(sum1,small2)
        if(ichoix==4)salref_t(i,j,k   )=sum2/max(sum1,small2)
#endif

       else                !sssssssssss>

! Si k0>=kmaxp1 cela signifie que nous sommes dans une zone decouvertes. Dans cette situation on !14-01-10
! utilise la valeur de surface de l'ogcm:
        if(k0>=kmaxp1) then !kkkkkkkkkkkkk>

           x2=(1.-rapi)*(1.-rapj)*ogcm_var3d(i1  ,j1  ,max_z)           &
             +(1.-rapi)*    rapj *ogcm_var3d(i1  ,j1+1,max_z)           &
             +    rapi *(1.-rapj)*ogcm_var3d(i1+1,j1  ,max_z)           &
             +    rapi *    rapj *ogcm_var3d(i1+1,j1+1,max_z)

         if(ichoix==1)temobc_t(i,j,k,t_)=x2
         if(ichoix==2)salobc_t(i,j,k,t_)=x2
#ifdef bidonref
         if(ichoix==3)temref_t(i,j,k)=x2
         if(ichoix==4)salref_t(i,j,k)=x2
#endif

        else            !kkkkkkkkkkkkk>


      if(ksecu==0) then !00000000000>          !29-07-10
        sum1=0.
        sum2=0.
        do k1=max_z-1,1,-1 ! debut boucle K1
         z1=ogcm_z1d(k1)
         z2=ogcm_z1d(k1+1)
         if(k1==max_z-1)z2=max(z2,depth_w(i,j,kmaxp1))          !23-11-11
!        zup_   =min(depth_w(i,j,k0+1)-depth_w(i,j,kmaxp1),z2)
!        zdw_   =max(depth_w(i,j,k0  )-depth_w(i,j,kmaxp1),z1)
         zup_   =min(depth_w(i,j,k0+1)                    ,z2)  !23-11-11
         zdw_   =max(depth_w(i,j,k0  )                    ,z1)
         if(zup_   -zdw_   .gt.zero) then !-_-_-_-_-_-_>
           x2=ogcm_var1d(k1+1)
           x1=ogcm_var1d(k1  )
           rapk=(zdw_   -z1)/(z2-z1)
           vdw_   =(1.-rapk)*x1+rapk*x2
           rapk=(zup_   -z1)/(z2-z1)
           vup_   =(1.-rapk)*x1+rapk*x2
           sum1=sum1+(zup_   -zdw_   )
           sum2=sum2+(zup_   -zdw_   )*0.5*(vup_   +vdw_   )
         endif                            !-_-_-_-_-_-_>
        enddo               ! fin K1
        if(ichoix==1)const2=temobc_t(i,j,k0,t_)-sum2/max(sum1,small2)
        if(ichoix==2)const2=salobc_t(i,j,k0,t_)-sum2/max(sum1,small2)
#ifdef bidonref
        if(ichoix==3)const2=temref_t(i,j,k0)   -sum2/max(sum1,small2)
        if(ichoix==4)const2=salref_t(i,j,k0)   -sum2/max(sum1,small2)
#endif
        ksecu=1
      endif             !00000000000>          !29-07-10

        sum1=0.
        sum2=0.
        do k1=max_z-1,1,-1 ! debut boucle K1
         z1=ogcm_z1d(k1)
         z2=ogcm_z1d(k1+1)
         if(k1==max_z-1)z2=max(z2,depth_w(i,j,kmaxp1))          !23-11-11
!        zup_   =min(depth_w(i,j,k+1)-depth_w(i,j,kmaxp1),z2)
!        zdw_   =max(depth_w(i,j,k  )-depth_w(i,j,kmaxp1),z1)
         zup_   =min(depth_w(i,j,k+1)                    ,z2)   !23-11-11
         zdw_   =max(depth_w(i,j,k  )                    ,z1)
         if(zup_   -zdw_   .gt.zero) then !-_-_-_-_-_-_>
           x2=ogcm_var1d(k1+1)
           x1=ogcm_var1d(k1  )
           rapk=(zdw_   -z1)/(z2-z1)
           vdw_   =(1.-rapk)*x1+rapk*x2
           rapk=(zup_   -z1)/(z2-z1)
           vup_   =(1.-rapk)*x1+rapk*x2
           sum1=sum1+(zup_   -zdw_   )
           sum2=sum2+(zup_   -zdw_   )*0.5*(vup_   +vdw_   )
         endif                            !-_-_-_-_-_-_>
        enddo               ! fin K1
        if(sum1==0)call ogcm_error_diag(ichoix)
        if(ichoix==1)temobc_t(i,j,k,t_)=const2+sum2/max(sum1,small2) !29-07-10
        if(ichoix==2)salobc_t(i,j,k,t_)=const2+sum2/max(sum1,small2)
#ifdef bidonref
        if(ichoix==3)temref_t(i,j,k   )=const2+sum2/max(sum1,small2) !29-07-10
        if(ichoix==4)salref_t(i,j,k   )=const2+sum2/max(sum1,small2)
#endif

        endif           !kkkkkkkkkkkkk>

       endif               !sssssssssss>

       enddo   ! fin K

      endif                         !>>>>>>>>>>>>>>
      enddo   ! fin J
      enddo   ! Fin I
      close(3)
      close(4)

      end subroutine ogcm_interp_core1

!................................................................

      subroutine ogcm_interp_core2(ichoix,t_,nc_)
      implicit none
      integer ichoix,t_,nc_
#ifdef synopsis
       subroutinetitle='ogcm_interp_core2'
       subroutinedescription= &
          'Reads T or S field in the netcdf OGCM file and interpolates' &
       //' on the S grid and stores in temref_z ans salref_z'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      open(unit=3 &
             ,file=trim(tmpdirname)//'hr_to_lr_coord_1_'//dom_c//'.out')

! Interpolation de l'etat de reference sur tableaux "geo z" temref_z et salref_z
      do j=0,jmax+1
      do i=0,imax+1

      read(3,*)i0,j0,i1,j1,x5,rapi,rapj
      if(i0.ne.i)stop 'erreur1 sur i'
      if(j0.ne.j)stop 'erreur1 sur j'

      if(mask_t(i,j,kmax)==1) then !>>>>>>>>>>>>

       k1=max_z
       do k=kmax,kmin_w(i,j),-1

        z1= (1.-rapi)*(1.-rapj)*ogcm_z(i1  ,j1  ,k1)       &
           +(1.-rapi)*    rapj *ogcm_z(i1  ,j1+1,k1)       &
           +    rapi *(1.-rapj)*ogcm_z(i1+1,j1  ,k1)       &
           +    rapi *    rapj *ogcm_z(i1+1,j1+1,k1)

        do while ( z1>zref_z(k) .and. k1>1)

         z2=z1
         k1=k1-1
         z1= (1.-rapi)*(1.-rapj)*ogcm_z(i1  ,j1  ,k1)      &
            +(1.-rapi)*    rapj *ogcm_z(i1  ,j1+1,k1)      &
            +    rapi *(1.-rapj)*ogcm_z(i1+1,j1  ,k1)      &
            +    rapi *    rapj *ogcm_z(i1+1,j1+1,k1)

        enddo
        if(k1==max_z) then !ooo>
           x0=(1.-rapi)*(1.-rapj)*ogcm_var3d(i1  ,j1  ,k1  )            &
             +(1.-rapi)*    rapj *ogcm_var3d(i1  ,j1+1,k1  )            &
             +    rapi *(1.-rapj)*ogcm_var3d(i1+1,j1  ,k1  )            &
             +    rapi *    rapj *ogcm_var3d(i1+1,j1+1,k1  )
        else               !ooo>

           x2=(1.-rapi)*(1.-rapj)*ogcm_var3d(i1  ,j1  ,k1+1)            &
             +(1.-rapi)*    rapj *ogcm_var3d(i1  ,j1+1,k1+1)            &
             +    rapi *(1.-rapj)*ogcm_var3d(i1+1,j1  ,k1+1)            &
             +    rapi *    rapj *ogcm_var3d(i1+1,j1+1,k1+1)
           x1=(1.-rapi)*(1.-rapj)*ogcm_var3d(i1  ,j1  ,k1  )            &
             +(1.-rapi)*    rapj *ogcm_var3d(i1  ,j1+1,k1  )            &
             +    rapi *(1.-rapj)*ogcm_var3d(i1+1,j1  ,k1  )            &
             +    rapi *    rapj *ogcm_var3d(i1+1,j1+1,k1  )
           
! Cas particulier zref_z est en dessous de k1 ou cas z2=z1 qui peut
! etre associe aux couches fusionnees (si ogcm=s26 couches fusionees)
           if(z2==z1) then !pmx>
            rap=0.
! par contre cette situation devant entrainer k1=1 on active le
! debugage dans le cas contraire:
              if(k1/=1) then !debug> !15-11-17
               write(10+par%rank,*)'par%rank,i,j,k pour temref',par%rank,i,j,k
               write(10+par%rank,*)'i1,j1,k1',i1,j1,k1
               write(10+par%rank,*)'z2       =',z2
               write(10+par%rank,*)'z1       =',z1
               write(10+par%rank,*)'zref_z(k)=',zref_z(k)
               stop 'Err 5128 ogcm_interp_core2 z2-z1=0 see fort files'
              endif          !debug>
           else            !pmx>
! CAS STANDARD:
            rap=max((zref_z(k)-z1)/(z2-z1),0.) ! cas standard
           endif           !pmx>

           x0=rap*x2+(1.-rap)*x1

        endif              !ooo>

#ifdef bidonref
        if(ichoix==3)temref_z(i,j,k)=x0
        if(ichoix==4)salref_z(i,j,k)=x0
#endif

       enddo   ! fin K

      endif                         !>>>>>>>>>>>>>>
      enddo   ! fin J
      enddo   ! Fin I

      close(3)

      end subroutine ogcm_interp_core2

!................................................................

      subroutine ogcm_interp_core3(ichoix,t_,nc_)
      implicit none
      integer ichoix,t_,nc_
#ifdef synopsis
       subroutinetitle='ogcm_interp_core2'
       subroutinedescription= &
          'Reads T or S field in the netcdf OGCM file and interpolates' &
       //' on the S grid and stores in temref_z ans salref_z'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Le principe de cette interpolation: egaler les deux integrales
! verticales (mercator, s26) pour avoir la meme quantite de sel ou de
! chaleur cumulee dans les deux colonnes (mercator, s26). C'est un
! principe de conservation de bilan.

      open(unit=3 &
             ,file=trim(tmpdirname)//'hr_to_lr_coord_1_'//dom_c//'.out')

! Interpolation de l'etat de reference sur tableaux "geo z" temref_z et salref_z
      do j=0,jmax+1
      do i=0,imax+1

      read(3,*)i0,j0,i1,j1,x5,rapi,rapj
      if(i0.ne.i)stop 'erreur1 sur i'
      if(j0.ne.j)stop 'erreur1 sur j'

      if(mask_t(i,j,kmax)==1) then !>>>>>>>>>>>>

        k1=max_z                                          ! niveau de surface de mercator
        z2=(1.-rapi)*(1.-rapj)*    ogcm_z(i1  ,j1  ,k1) & ! z correspondant
          +(1.-rapi)*    rapj *    ogcm_z(i1  ,j1+1,k1) &
          +    rapi *(1.-rapj)*    ogcm_z(i1+1,j1  ,k1) &
          +    rapi *    rapj *    ogcm_z(i1+1,j1+1,k1)
        x2=(1.-rapi)*(1.-rapj)*ogcm_var3d(i1  ,j1  ,k1) & ! Valeur correspondate
          +(1.-rapi)*    rapj *ogcm_var3d(i1  ,j1+1,k1) &
          +    rapi *(1.-rapj)*ogcm_var3d(i1+1,j1  ,k1) &
          +    rapi *    rapj *ogcm_var3d(i1+1,j1+1,k1)

        if(depth_w(i,j,1)>z2) then !sssssss>
! Cas particulier zone seche:
! Si le fond de s26 est au dessus de la surface de mercator. 

         if(ichoix==1)temobc_t(i,j,:,t_)=x2
         if(ichoix==2)salobc_t(i,j,:,t_)=x2
        
        else                       !sssssss>

! Cas standard:
! On ajuste la surface du champ mercator a celle du modele:
       z2=depth_w(i,j,kmax+1)
       sum1=0.  ; sum3=0.  

        k1=max_z-1                                        ! niveau mercattor

        z1=(1.-rapi)*(1.-rapj)*    ogcm_z(i1  ,j1  ,k1) & ! z correspondant
          +(1.-rapi)*    rapj *    ogcm_z(i1  ,j1+1,k1) & ! au point (i,j) de s26
          +    rapi *(1.-rapj)*    ogcm_z(i1+1,j1  ,k1) &
          +    rapi *    rapj *    ogcm_z(i1+1,j1+1,k1)
        z1=min(z1,z2) ! On a ajuste z2 a la ssh de s26 on assure donc que z1 reste sous z2

        x1=(1.-rapi)*(1.-rapj)*ogcm_var3d(i1  ,j1  ,k1) & ! Valeur correspondante
          +(1.-rapi)*    rapj *ogcm_var3d(i1  ,j1+1,k1) &
          +    rapi *(1.-rapj)*ogcm_var3d(i1+1,j1  ,k1) &
          +    rapi *    rapj *ogcm_var3d(i1+1,j1+1,k1)

       do k=kmax,kmin_w(i,j),-1

        do while ( z1>depth_w(i,j,k).and. k1>1) !wwwwwwww>

         sum1=sum1+(z2-z1)*0.5*(x1+x2) ! integrale verticale (trapeze) colonne mercator

         z2=z1
         x2=x1
         k1=k1-1
         z1=(1.-rapi)*(1.-rapj)*    ogcm_z(i1  ,j1  ,k1) & ! z
           +(1.-rapi)*    rapj *    ogcm_z(i1  ,j1+1,k1) &
           +    rapi *(1.-rapj)*    ogcm_z(i1+1,j1  ,k1) &
           +    rapi *    rapj *    ogcm_z(i1+1,j1+1,k1)
         z1=min(z1,z2) ! On a ajuste z2 a la ssh de s26 on assure donc que z1 reste sous z2

         x1=(1.-rapi)*(1.-rapj)*ogcm_var3d(i1  ,j1  ,k1) & ! val(z)
           +(1.-rapi)*    rapj *ogcm_var3d(i1  ,j1+1,k1) &
           +    rapi *(1.-rapj)*ogcm_var3d(i1+1,j1  ,k1) &
           +    rapi *    rapj *ogcm_var3d(i1+1,j1+1,k1)

        enddo                                   !wwwwwwww>

! On sort de la boucle "while" avec z1<depth_w(i,j,k)<z2
! Le reste d'integrale de mercator va de z2 ou la valeur est x2 a depth_w a laquelle il
! faut donner une valeur (x1). Si le z(1) de  mercator reste au dessus du
! niveau i,j,k alors max(rap,0) prolonge mercator d'un gradient nul
        rap=max((depth_w(i,j,k)-z1)/(z2-z1),0.)
! sum11 est le reste  de integrale:
        sum11=(z2-depth_w(i,j,k))*0.5*(x2+rap*x2+(1.-rap)*x1)

! x0 valeur au point i,j,k de s26:
        x0=(sum1+sum11-sum3)/dz_t(i,j,k,1)

! Poursuite integrale sur colonne s26
        sum3=sum3+x0*dz_t(i,j,k,1)

        if(ichoix==1)temobc_t(i,j,k,t_)=x0
        if(ichoix==2)salobc_t(i,j,k,t_)=x0

       enddo   ! fin k

! Couche fusionnee: !16-10-17
       if(ichoix==1)temobc_t(i,j,1:kmin_w(i,j)-1,t_)=temobc_t(i,j,kmin_w(i,j),t_)
       if(ichoix==2)salobc_t(i,j,1:kmin_w(i,j)-1,t_)=salobc_t(i,j,kmin_w(i,j),t_)

       endif                      !sssssss>

      endif                         !>>>>>>>>>>>>>>
      enddo   ! fin j
      enddo   ! fin i

      close(3)

      end subroutine ogcm_interp_core3

!...................................................................

      subroutine ogcm_oldlist(var_count_) !18-04-16
      implicit none
      integer var_count_,flag_,unit_

! Si la liste binaire est presente dans le tmp alors se positionner
! dans celle-ci pour demarrer la simulation:

        if(var_count_==0)return

        flag_=0
        unit_=s_unit(7)
        open(unit=unit_,file=texte30        &
                       ,access='direct'     &
                       ,recl=541            & !28-04-16
                       ,form='unformatted'  &
                       ,status='old'        &
                       ,iostat=k0)
         if(k0/=0)stop 'Err 2424'
         nc=0

         do while (flag_==0)

          nc=nc+1

               ogcm_rec_prev(var_count_)=     ogcm_rec_next(var_count_)
          ogcm_readtime_prev(var_count_)=ogcm_readtime_next(var_count_) 
            ogcm_period_prev(var_count_)=  ogcm_period_next(var_count_)


         ogcm_rec_next(var_count_)=nc
         if(nc/=1) then !pmxpmx>
          read(unit_,rec=nc,iostat=i0)texte250,i1          &
                           ,ogcm_readtime_next(var_count_) &
                             ,ogcm_period_next(var_count_)
         else           !pmxpmx>
! Verifier la date repere de nc=1:
          read(unit_,rec=nc,iostat=i0)texte250,i1          &
                           ,ogcm_readtime_next(var_count_) &
                             ,ogcm_period_next(var_count_) &
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

         if(i0/=0)stop 'Err 5228'

         ogcm_rec_max=nc

! On cherche les 2 records dont les temps encadrent le temps initial elapsedtime_now
         if(nc>1) then !pmx>
          if(     elapsedtime_now>=ogcm_readtime_prev(var_count_) &
             .and.elapsedtime_now< ogcm_readtime_next(var_count_))flag_=1
         endif         !pmx>

         enddo

! Poursuivre la lecture jusqu'au bout du fichier pour determiner ogcm_rec_max
         do while (i0==0)
          nc=nc+1
          read(unit_,rec=nc,iostat=i0)texte250,i1,x1,x2
          ogcm_rec_max=nc-1
!         write(6,*)'nc i0 ogcm_rec_max',nc,i0,ogcm_rec_max
         enddo

         close(unit_)

        if(flag_==0) STOP 'Pas de fichiers OGCM autour du temps initial'

        if(par%rank==0) then !ooo>
         write(6,*)'k1=                           ',k1
         write(6,*)'var_count_=                   ',var_count_
         write(6,*)'elapsedtime_now=              ',elapsedtime_now
         write(6,*)'ogcm_readtime_prev(var_count_)',ogcm_readtime_prev(var_count_)
         write(6,*)'ogcm_readtime_next(var_count_)',ogcm_readtime_next(var_count_)
         write(6,*)'ogcm_rec_prev(var_count_)     ',ogcm_rec_prev(var_count_)
         write(6,*)'ogcm_rec_next(var_count_)     ',ogcm_rec_next(var_count_)
        endif                !ooo>

      end subroutine ogcm_oldlist

!...................................................................

      subroutine ogcm_stop(txt_,case_)
      implicit none
      character(len=*) :: txt_
      integer :: case_,flag_stop_glb_=0

      if(case_==1) then !case1>
       flag_stop=1
       write(10+par%rank,'(a)')'.......................'
       write(10+par%rank,'(a)')' Err module_ogcm '
       write(10+par%rank,'(a)')trim(txt_)
      endif             !case1>

      if(case_==2) then !case2>
       call mpi_allreduce(flag_stop,flag_stop_glb_,1,mpi_integer,mpi_sum,par%comm2d,ierr)
       if(flag_stop_glb_/=0) stop 'Err ogcm_stop see fort.xx files'
      endif             !case2>


      end subroutine ogcm_stop

!...................................................................

      subroutine ogcm_gridchangedetection(var_id_,nc_)  !17-05-19
      implicit none
      integer nc_,var_id_,recnum_

! Detection de changement de grille pour remise A jour des parametres d'interpolation.
! https://docs.google.com/document/d/1VMmzC7jXV7sycAxoFCUcOBxVbydbyXD3imgM0wAnAFQ/edit


! Rappel: les variables sont regroupees par identifiant: 
! T,S pour trc_id pour grille t
! SSH pour ssh_id pour grille t (frequence eventuellement differente de trc_id) 
! U,V pour vel_id pour grille u et v => si grille u change alors grille v change aussi 

       flag_updatebouchetrou=0 !17-05-19

! Obtenir le nom du fichier de variable:
      texte90='none'

      if(obc_ogcm_type(1:5)=='sympa') texte90=liste_ogcm_uvts
      if(var_id_==ssh_id) then !eeeeeeeeeee> ! 'e' pour elevation de la surface
       if(obc_ogcm_type(1:4)=='ncom')  texte90=liste_ogcm_ssh
       if(obc_ogcm_type(1:6)=='nemo_z')texte90=liste_ogcm_ssh
      endif                    !eeeeeeeeeee>
      if(var_id_==trc_id) then !ttttttttttt>
       if(obc_ogcm_type(1:4)=='ncom')  texte90=liste_ogcm_t
       if(obc_ogcm_type(1:6)=='nemo_z')texte90=liste_ogcm_t
      endif                    !ttttttttttt>
      if(var_id_==vel_id) then !uuuuuuuuuuu>
       if(obc_ogcm_type(1:4)=='ncom')  texte90=liste_ogcm_u
       if(obc_ogcm_type(1:6)=='nemo_z')texte90=liste_ogcm_u
      endif                    !uuuuuuuuuuu>
      if(texte90=='none')stop ' nom var invalide dans module_ogcm'

!-----------------------------------------------------------------------------------------!
! Seul par%rank=0 lit le fichier puis envoie les valeurs aux autres par%rank par mpi_bcast
      if(par%rank==0) then !--[m0v0m]--> !05-05-18
       open(unit=4,file=trim(texte90) &
                  ,access='direct'    &
                  ,recl=541           & !28-04-16
                  ,form='unformatted')
       read(4,rec=nc_   )texte250,ogcm_time_counter           &
                                 ,ogcm_readtime_next(var_id_) &
                                   ,ogcm_period_next(var_id_) &
                                   ,ogcmgridcode !16-05-19
       close(4)
      endif                !--[m0v0m]--> !05-05-18

! le tableau ogcmgridnum est le numero de record pour acceder au fichier
! de grille associee A la variable considerEe. Archivant 2 etats
! successifs il permet de savoir, variable par variable, si la grille associee
! est changee afin de declencher une remise A jour des parametres d'interpolation
! Details dans:
! https://docs.google.com/document/d/1VMmzC7jXV7sycAxoFCUcOBxVbydbyXD3imgM0wAnAFQ/edit
      call mpi_bcast(ogcmgridcode,1,mpi_integer,0,par%comm2d,ierr) !16-05-19
      ogcmgridnum(var_id_,1)=ogcmgridnum(var_id_,2)
      ogcmgridnum(var_id_,2)=ogcmgridcode 

! Si changement de grille refaire les parametres d'interpolation sauf
! pour salinite (qui va avec temperature), sachant que SSH n'est pas
! consideree comme allant avec T puisque frequence differente possible
! et sachant que grilles U et V sont liEes.
! https://docs.google.com/document/d/1VMmzC7jXV7sycAxoFCUcOBxVbydbyXD3imgM0wAnAFQ/edit
! Le driver  ogcm_grid_to_grid_pathway permet de relancer toute la chaine de remise A jour 
! des parametres d'interpolation
      if(ogcmgridnum(var_id_,2)/=ogcmgridnum(var_id_,1)) then !m°v°m>
        recnum_=ogcmgridnum(var_id_,2)
        if(par%rank==0)write(6,*) &
        'Detection changement de grille var_id=',var_id_

       if(var_id_==trc_id)call ogcm_grid_to_grid_pathway(1,recnum_)
       if(var_id_==ssh_id)call ogcm_grid_to_grid_pathway(1,recnum_)
       if(var_id_==vel_id)call ogcm_grid_to_grid_pathway(2,recnum_)
       if(var_id_==vel_id)call ogcm_grid_to_grid_pathway(3,recnum_)

! Dans le cas des traceurs il faudra egalement remettre A jour les trous
! et leur bouchage:
       if(var_id_==trc_id)flag_updatebouchetrou=1 !17-05-19

      endif                                                   !m°v°m>


!-----------------------------------------------------------------------------------------!

      end subroutine ogcm_gridchangedetection

!.....................................................................

      end module module_ogcm
