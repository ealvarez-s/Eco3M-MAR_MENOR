      subroutine initial_tide
!______________________________________________________________________
! SYMPHONIE ocean model
! release 310 last update: 03-11-21
!______________________________________________________________________

      use module_principal
      use module_parallele !#mpi
      use module_systeme
      implicit none
#ifdef synopsis
       subroutinetitle='initial_tide'
       subroutinedescription= &
       'Computes the initial state of the tidal fields'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!    _________                    .__                  .__             !m°v°m 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !
!...............................................................................
! version date      description des modifications
!         01/01/02: version initiale
!         07/01/02: bienvenue à TIDESTEP3
!         27/01/02: extension accrue des fichiers d'entree de maree
!         18/02/02: reset de tableaux
!         22/02/02: amenagement 3D PASSETIDE augmente de taille
!         26/02/02: PASSETIDE prend en compte la possibilité de commencer la
!                   sommation par la marée et de continuer par le thermohalin.
!         05/03/02: compatibilite avec f95
!         26/08/02: U1TIDE_X U2TIDE_X V1TIDE_Y V2TIDE_Y sont maintenant des
!                   composantes de courant moyen et non plus de transport
!         24/07/03: nouveaux formats de lecture des fichiers de marée de
!                   Florent, + lecture de fichiers de courant de marée
!         26/07/03: suite.....
!         21/08/03: ajout d'une possibilité de selectionner les forces de
!                   marée: ondes + tout potentiel
!                           ondes + potentiel astro
!                           ondes + potentiel charge
!                           ondes + aucun potentiel
!                   c.a.d. bienvenue à TIDEFORCES
!                   supression TIDESTEP1
!                   bienvenue à TIDEVEL (choix entre courant fichier et
!                                        courant deduit d'un model_ lineaire
!                                        spectral)
!         28/07/04: MTDE et NTDE remplacent MECO et NECO
!                   Bienvenue à ONOFF_TIDE
!                   Bienvenue à ANALYSETIDE_X
!                   Bienvenue à ANALYSETIDE_Y
!                   TIDESTEP2 TIDESTEP3 IN_OUT_TIDE ne sont plus lus dans
!                   notebook_tide mais forcés en dur dans la routine
!                   d'initialisation
!         22/02/06: debugage "latitude factor" de la classe "long period"
!         05/04/07: Cas d'une simulation imbriquee où la maree est deja modelisee
!                   dans le model_ maman. Dans ce cas, la maree aux limites ouvertes
!                   vient du model_ maman et pas de mog2d: seul le potentiel et
!                   l'analyse sont activés.
!         18/01/08: Bienvenue à ZTAMEATIDE (moyennes spatiales par harmonique)
!                           et ZTATOTIDE_Z (sla de maree totale)
!         12-05-09  Mieux calculer la moyenne de la maree
!         16-05-09  vers en version netcdf finalisee...
!         13-06-09  introduction de C.L. de type Z0 pour l'amplitude et la phase
!                   de la marée. Requis pour Z1TIDE et Z2TIDE puis mu_outputs.F
!         26-06-09  boucle sur imax,jmax remplacee par imax,jmax
! 2009.3  02-10-09  dte_lp remplace dte
!         09-10-09  - lecture des fichiers netcdf de Florent Lyard
!                   - parallelisation
!         08-11-09  optimobc_tide est renommé tide_analysis
! 2010.2  16-12-09  interpolation en ligne de la marée
!         22-12-09  initialisation parametres nodaux
!         27-12-09  suppression cas tidestep3 inutilisé
!         28-12-09  suppression tidevel
! 2010.3  08-01-10  - ne passe pas par calcul des parametres nodaux si pas de maree
!                   - le cas i2dh n'existe plus
! 2010.6  02-02-10  renomme lon_t, lat_t
! 2010.11 05-07-10  lecture adaptee aux fichiers d'analyse
!         13-07-10  lecture des fichiers de maree sirocco
! 2010.12 20-09-10  Possibilité de calcul en simple precision
! 2010.19 12-04-11  obc mpi sur tableaux de forcage de la maree
! 2010.20 16-04-11  Calculs sur la base d'un temps en secondes
! 2010.22 10-05-11  Ne passe pas par la routine tide_find_missing_value si proc
!                   totalement masqué
! 2010.24 15-10-11  "ifdef bidouille" commenté en attendant que Cyril ajoute une
!                    fonction echange de plus afin d'eviter l'erreur de parallelisation
!         25-10-11  sécurité sur boucle i10 de recherche de comblement de trou
!         27-10-11  debug du chenillard
!         20-11-11  - Possibilite de surdeterminer l'analyse avec une contrainte de rappel
!                   aux champs de maree utilises en forcage, via la constate weight_
!                   - Adapter la lecture des fichiers de reanalyse au calcul parallele
! 2010.25 26-02-12 supression mtde et ntde
!         28-03-12  echange compatible avec pgf Stelios
! S26     07-02-13  lecture des champs FES2012
!         15-02-13  Forcer le code avec une reanalyse S
!         25-02-13  suppresion de l'appel a tide_analysis(5)
!         12-09-13  modifs notebook_tide
!         11-11-13  champs FES2012 suite: correction faute sur lecture lon lat courant Oy
!         16-11-13  possibilite de lire LSA dans le fichier S reanalysis
!         26-03-14  debug interpolation
!         24-07-14  nouveaux echanges
!         11-11-14  le cas "grande grille de fes2012" est trop couteux en cpu on
!                   revient par conséquent à l'algo initial
!         12-11-14  pour conservation mpi
!         18-12-14  modif boucle pour convenance graphique
!         19-04-15  Une meilleure prise en compte de la periodicite de la grille fes2012
!         10-04-16  cas academiques: aucun forcage de maree mais analyse harmonique
!         09-09-16  sshtotide inutile
!         25-09-16  Possibilite d'ajouter "A la main" un offset sur la phase des ondes de maree
!         26-11-16  If SLA not provided (unfortunatly the case with FES fields) then use 
!                   R.D. Ray (Marine Geodesy 1998) approximation
!         11-03-17  Ne pas apliquer le potentiel de Ray 1998 sur des solutions re-analysee
!         12-03-17  Verifier la coherence des fichiers d'entree 
!         20-03-17  suppression sshmeatide
!         29-03-17  un fichier separE pour le LSA
!         02-05-17  debug cas particulier absence de courant dans fichier FES (ne pas
!                   empecher de charger les tableaux "out"
!         21-05-17  champ maree multipliE par le masque
!         11-06-17  - sur occigen passer par flag_stop pour arreter la simulation
!                   - lire filval
!         28-06-17  Modifs thom
!         15-04-18  verification de la coherence de notebook_tide
!         17-04-18  ajout subroutine tide_ssh_bias_correction 
! v260    12-10-19  amelioration procedure debugage
! v273    26-01-20  suppression test sur used_unused_dom (occigen se perd dans mpi)
!         28-01-20  ajout de nouveaux attributs nom de variable netcdf
! v284    20-05-20  Bouchage de la grille de maree conditionnE au fait que cela
!                   concerne uniquement des points symphonie en mer 
! v287    19-08-20  initialisation tideana_spinup
! v296    11-02-21  adaptation des noms netcdf aux fichiers FES2014b
! v310    03-11-21  correction ssh : pas si =filval
!...............................................................................

! La parallelisation doit prevoir le calcul des moyennes spatiales des harmoniques'

! REFERENCE:
!
! La recomposition du champ d'elevation de la surface et des courants de marée,
! ainsi que du potentiel astronomique, à partir des composantes harmoniques et
! des paramètres nodaux est détaillé dans l'Annexe de Pairaud et al, 2008:
!
! Pairaud I. L., Lyard F., Auclair F., Letellier T., Marsaleix P., 2008,
! Dynamics of the semi-diurnal and quarter-diurnal internal tides in the
! Bay of Biscay. Part 1: Barotropic tides,
! Continental Shelf Research,28, 1294-1315. doi:10.1016/j.csr.2008.03.004

! reset:
      tidestep2=0                                                      !28/07/04
      tidestep3=0                                                      !28/07/04
      in_out_tide=1                                                    !28/07/04
      tidepotential_w(:,:,:)=0.                            !28/07/04
!     sshtotide_w(:,:,:)=0.                                !18/01/08
      if(kmaxtide==0)return                                        !18/02/02

      analysetide_u(:,:,:)=0.                              !28/07/04
      analysetide_u(:,:,:)=0.                              !28/07/04
      analysetide_v(:,:,:)=0.                              !28/07/04
      analysetide_v(:,:,:)=0.                              !28/07/04
      analysetide_w(:,:,:)=0.
      analysetide_w(:,:,:)=0.
      potidecos_w(:,:,:)=0.
      potidesin_w(:,:,:)=0.
      sshtidecos_w(:,:,:)=0.
      sshtidesin_w(:,:,:)=0.
      veltidecos_u(:,:,:)=0.
      veltidesin_u(:,:,:)=0.
      veltidecos_v(:,:,:)=0.
      veltidesin_v(:,:,:)=0.

      sshtidecosout_w(:,:,:)=0.
      sshtidesinout_w(:,:,:)=0.
      veltidecosout_u(:,:,:)=0.
      veltidesinout_u(:,:,:)=0.
      veltidecosout_v(:,:,:)=0.
      veltidesinout_v(:,:,:)=0.

!.......
! initialiser les parametres nodaux et les kounts de rdv:              !08-01-10
      call tide_nodal_parameters                                       !22-12-09

!     rap=(       real(kount)-kounttide_prev_rdv)                       &
!        /(kounttide_next_rdv-kounttide_prev_rdv)
      rap=(   elapsedtime_now-tidenodal_prev_rdv)           & ! 16-04-11
         /(tidenodal_next_rdv-tidenodal_prev_rdv)

      do ktide=1,kmaxtide
       utide(ktide,1)=(1.-rap)*utide(ktide,0)+rap*utide(ktide,2)
       ftide(ktide,1)=(1.-rap)*ftide(ktide,0)+rap*ftide(ktide,2)
      enddo
!.......


! Analyse harmonique:
! définir le modulo d'echantillonage en nbre d'itération
!     tideana_modulo=max0( nint(tideana_delta*3600./dti_fw) , 1 )      !08-01-10
      tideana_modulo=max(tideana_delta*3600.,un)                       !16-04-11
!     tideana_nextime=elapsedtime_now+tideana_spinup*86400.            !06-11-13
      tideana_nextime=             0.+tideana_spinup*86400.            !19-08-20
! Note: bien que revenant au meme car tideana_nextime est dans le fichier
! restart (si on repart d'un restart on ne refait pas la partie du spinup deja
! faite au run prededent), il est tout de meme plus logique que tideana_nextime soit referencE
! par rapport au temps 0 (19-08-20)

!......................................................................!22/02/02
! PASSETIDE permet d'initier les sommations afin d'empiler les composantes
! harmoniques. 1er  indice: le numero de l'harmonique
!              2eme indice: permet de distinguer les variables zeta courant
!                           d'une part (K=1) du potentiel d'autre part (K=2).
      do k=1,2
      do ktide=2,kmaxtide
      passetide(ktide,k)=1.
      enddo
      enddo
      passetide(1,1)=0.
      passetide(1,2)=0.

!*******************************************************************************
! LECTURE DES FICHIERS D'ENTREE: AMPLITUDE ET PHASE DES ONDES          !18/01/02
! ET DU POTENTIEL DE CHARGE DE LA MAREE SUR ELLE MEME
! VOIR APEL PAGES 214,215 Attention Apel utilise des colatitudes et nous des
! latitudes
! DEBUT:
!*******************************************************************************


! Note chaque parametre fait l'objet de sa propre boucle pour eventuellement
! remettre A jour les regles d'interpolation

      if(tideforces<4) then ! obc obc >
!....................
! Lire le fichier ssh:
       do ktide=1,kmaxtide        !  m[0_0]m  >>>>

        if(tide_interpolation==0) then               !--------------->  !12-09-13
         call tide_previous_run_analysis('h0') ! case: interpolation preprocessed on s grid
        else                                         !--------------->
         call tide_netcdf_read_h(1)            ! case: perform the interpolation on s grid
        endif                                        !--------------->

       enddo ! boucle sur ktide   !  m[0_0]m  >>>>
      endif                 ! obc obc >

!....................
! Lire le fichier du potentiel de charge:
      if(tideforces==0.or.                       &
         tideforces==2.or.                       &
         tideforces==4) then !ppppppppppppp>

       do ktide=1,kmaxtide        !  (°L°)  >>>>
   
        if(index(nametide(ktide,5),'none')==0) then ! m[0v0]m >

         call tide_netcdf_read_h(2) ! potentiel de charge     !16-05-09

        else                                        ! m[0v0]m >

        ! le test if(tide_interpolation/=0) signifie qu'une solution deja
        ! interpolee est interpretee comme une re-analyse qui peut contenir des
        ! defauts (par ex dans les zones intertidales, les fleuves etc...) qu'on
        ! ne souhaite pas voir introduits dans le potentiel de charge. La methode
        ! de Ray 1998 n'est par consequent utilisee que pour les solutions de
        ! type FES
         if(tide_interpolation/=0) then !pmxpmx> !11-03-17
          ! If SLA not provided (unfortunatly the case with FES fields) then use Ray(1998)
          ! but an accurate SLA field is a much better option.
          ! R. D. Ray (1998) Ocean self¿attraction and loading in numerical tidal
          ! models, Marine Geodesy, 21:3, 181-192, http://dx.doi.org/10.1080/01490419809388134
           do j=0,jmax+1 ; do i=0,imax+1
            potidecos_w(i,j,ktide)=0.085*sshtidecos_w(i,j,ktide)
            potidesin_w(i,j,ktide)=0.085*sshtidesin_w(i,j,ktide)
           enddo ; enddo
         endif                          !pmxpmx> !11-03-17


        endif                                       ! m[0v0]m >

       enddo  ! ktide loop        !  (°L°)  >>>>


      endif                  !pppppppppppp>


!                         /   /   /

!...............................................................................
! Ajout du potentiel astronomique au potentiel de charge:
      if(tideforces==0.or.                       &
         tideforces==1.or.                       &
         tideforces==4) then !-ASTRO->

      do ktide=1,kmaxtide        !  )°o°(  >>>>

       x0=0.
       x1=0.
       x2=0.
       if(nutide(ktide).eq.0)x0=1.
       if(nutide(ktide).eq.1)x1=1.
       if(nutide(ktide).eq.2)x2=1.

!      const1=0.7*equitide(ktide)*ftide(ktide)
       const1=0.7*equitide(ktide)
       if(tideforces==2)const1=0.                                     !21/08/03
       if(tideforces==3)const1=0.                                     !21/08/03
       if(tideforces==5)const1=0.                                     !10-04-16

       do j=0,jmax+1
       do i=0,imax+1

! le potentiel composé de l'effet de la marée sur elle même et des effets des
! astres:
! Note: on passe par les coefs complexes pour effectuer l'addition

         potidecos_w(i,j,ktide)=potidecos_w(i,j,ktide)                & !24-12-09
         +const1*( x0*0.5*(1.-3.*sin( lat_t(i,j))**2)                 & !effets des astres!22/02/06
                     +x1*     sin(2.* lat_t(i,j))                     &
                     +x2*        cos( lat_t(i,j))**2 )*               &
                    cos(nutide(ktide)*lon_t(i,j))

         potidesin_w(i,j,ktide)=potidesin_w(i,j,ktide)                & !24-12-09
         -const1*( x0*0.5*(1.-3.*sin( lat_t(i,j))**2)                 & !effets des astres!22/02/06
                     +x1*     sin(2.* lat_t(i,j))                     &
                     +x2*        cos( lat_t(i,j))**2 )*               &
                    sin(nutide(ktide)*lon_t(i,j))


       enddo
       enddo

      enddo  ! ktide loop        !  )°o°(  >>>>

      endif               !-ASTRO->



!...............................................................................
! lecture courant
       do ktide=1,kmaxtide        !  @°v°@  >>>>

         if(obc2dtype==1.and.trim(nametide(ktide,2))=='none') then !ooo> !15-04-18
stop 'notebook_tide inconsistent choice: obc2dtype==1 and no uv file??? Hum....'
         endif                                                    !ooo>

!     if(tideforces<4.and.obc2dtype>0.) then !obc obc> !02-05-17
      if(tideforces<4.and.trim(nametide(ktide,2))/='none') then !obc obc> !02-05-17


! Lire le fichier des vitesses:
        if(tide_interpolation==0) then               !--------------->  !12-09-13
         call tide_previous_run_analysis('u0') ! case: interpolation preprocessed on s grid
         call tide_previous_run_analysis('v0') ! case: interpolation preprocessed on s grid
        else                                         !--------------->
         if (obc2dtype>0.) call tide_netcdf_read_uv  ! case: perform the interpolation on s grid !28-06-17
        endif                                        !--------------->
      endif                                  !obc obc>

        if(nametide(ktide,3)/='none') then                ! r r r r r r r r r> !12-09-13
! Si une simulation precedente a fourni une re-analyse de la maree
! celle ci servira eventuellement à corriger les sorties des effets de la maree.
         call tide_previous_run_analysis('h1')
         call tide_previous_run_analysis('u1')
         call tide_previous_run_analysis('v1')
        else                                              ! r r r r r r r r r>
! Sinon ce sont les tableaux de forcage qui serviront à corriger les sorties
! des effets de la marée:
         do j=0,jmax+1 ; do i=0,imax+1 
          sshtidecosout_w(i,j,ktide)=sshtidecos_w(i,j,ktide)
          sshtidesinout_w(i,j,ktide)=sshtidesin_w(i,j,ktide)
         enddo         ; enddo
         do j=1,jmax   ; do i=1,imax+1
          veltidecosout_u(i,j,ktide)=veltidecos_u(i,j,ktide)
          veltidesinout_u(i,j,ktide)=veltidesin_u(i,j,ktide)
         enddo         ; enddo
         do j=1,jmax+1 ; do i=1,imax
          veltidecosout_v(i,j,ktide)=veltidecos_v(i,j,ktide)
          veltidesinout_v(i,j,ktide)=veltidesin_v(i,j,ktide)
         enddo         ; enddo
        endif                                             ! r r r r r r r r r>

       enddo  ! ktide loop        !  @°v°@  >>>>

! Correction eventuelle de la SSH 
       do ktide=1,kmaxtide
         if(index(nametide(ktide,6),'none')==0)call tide_ssh_bias_correction('add') !17-04-18
       enddo

! Echange mpi:
      call initial_tide_mpi_obc    !#mpi

      end subroutine initial_tide
!__________________________________________________________________________________
      subroutine tide_netcdf_read_h(ki1)
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
      include 'netcdf.inc'
      integer ki1,file1orfile2_,status1_,status2_,londim_,latdim_  &
             ,tabdim_(4)
      real , dimension(:) , allocatable :: lonr4_,latr4_
      character*30 varname_amplitude_,varname_phase_,bidon_name_
#ifdef synopsis
       subroutinetitle='tide_netcdf_read_h'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(ki1==1) then !>>>>>>>>>>
! On commence par verifier que la variable est bien dans le premier
! fichier, sinon on regarde dans le deuxieme fichier.
        varname_amplitude_='Ha' ; varname_phase_='Hg'
!................
! Ces lignes servent a savoir si Ha et Hg sont present dans le
! nametide(ktide,1) ou le ficher nametide(ktide,2)
       file1orfile2_=-999 ; k0=0
  379  k0=k0+1
       status=nf_open(nametide(ktide,k0),nf_nowrite,ncid1)
!                      status1_=nf_inq_varid(ncid1,trim(varname_amplitude_),var_id)
                       status1_=nf_inq_varid(ncid1,'Ha',var_id)
        if(status1_/=0)status1_=nf_inq_varid(ncid1,'elevation_a',var_id)
                       status2_=nf_inq_varid(ncid1,'Hg',var_id)
        if(status2_/=0)status2_=nf_inq_varid(ncid1,'elevation_G',var_id)
        if(status1_==0.and.status2_==0)file1orfile2_=k0
        if(status1_/=status2_)stop 'Erreur varname tide_netcdf_read_h'
       status=nf_close(ncid1)
       if(file1orfile2_==-999.and.k0==1) goto 379
       if(file1orfile2_==-999.and.k0==2) then !!!!!!!!!!>
        write(6,'(a,a,a)')'Tidal variables: ',varname_amplitude_,varname_phase_
        write(6,*)'not found in the specified following files:'
        write(6,'(a)')trim(nametide(ktide,1))
        write(6,'(a)')trim(nametide(ktide,2))
        stop 'STOP1 dans tide_netcdf_read_h'
       endif                                  !!!!!!!!!!>
      endif           !>>>>>>>>>>

      if(ki1==2) then !>>>>>>>>>>
        varname_amplitude_='LSAa' ; varname_phase_='LSAg' 
        file1orfile2_=5 !29-03-17
      endif           !>>>>>>>>>>

! lire l'identifiant du fichier:
      if(par%rank==0) then !>>>>>>>
      write(6,*)'dans tide_netcdf_read_h je lis:'
      write(6,'(a)')nametide(ktide,file1orfile2_)
      endif                !>>>>>>>

      status=nf_open(nametide(ktide,file1orfile2_),nf_nowrite,ncid1)
      if(status/=0) then !ooo>
         write(6,'(a,a)') 'File ',trim(nametide(ktide,file1orfile2_))
         stop ' Err 482 fichier maree absent'
         status=nf_close(ncid1)
         return
      endif              !ooo>


! lire l'identifiant des dimensions:
                   status=nf_inq_dimid(ncid1,'x_ZHL',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'imax_t',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'x',dim_x_id)            !13-07-10
      if(status/=0)status=nf_inq_dimid(ncid1,'nx',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'longitude',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'ni_t',dim_x_id)            !13-07-10
      if(status/=0)stop 'err 494 nom pour dim_x_id'

                   status=nf_inq_dimid(ncid1,'y_ZHL',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'jmax_t',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'y',dim_y_id)            !13-07-10
      if(status/=0)status=nf_inq_dimid(ncid1,'ny',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'latitude',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'nj_t',dim_y_id)
      if(status/=0)stop 'err 502 nom pour dim_y_id'

! Lecture des valeurs des dimensions
      status=nf_inq_dimlen(ncid1,dim_x_id,tide_imax)
      status=nf_inq_dimlen(ncid1,dim_y_id,tide_jmax)

! Garder en memoire la taille totale du domaine avant redecoupage:
      tidefull_imax=tide_imax
      tidefull_jmax=tide_jmax

!     k0=0
!     status=nf_inq_attlen(ncid1,nf_global,"grid_type",k0)
!     texte90='unknown'
!     if(k0/=0)status=nf_get_att_text(ncid1,nf_global,"grid_type",texte90(1:k0))
! Si l'attribut grid_type est "regular" (cas des fichiers de maree produits par LEGOS)
! on peut faire l'hypothèse que longitude (respectivement latitude) ne depend
! que de x (respect. y) et ainsi eviter de charger un tableau trop grand pour la
! capacité en memoire vive du modèle. La subroutine suivante explore un petit morceau
! du fichier pour verifier cette hypothèse et definir la fenêtre d'extraction:
!     if(texte90(1:7)=='regular')call tide_netcdf_zoom
      call tide_netcdf_zoom(1)


      call allocate_forcages(1,7,tide_imax,tide_jmax,2) ! arg1=allouer arg2=tide !09-10-09

! lire longitude et latitude:
                   status=nf_inq_varid(ncid1,'lon_t',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'longitude_t',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'longitude',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lon_ZHL',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lon',var_id)              !13-07-10
      if(status/=0)stop 'err nom longitude'

      status=nf_inq_var(ncid1,var_id,bidon_name_,k0,londim_,tabdim_,i4)

      if(londim_==2) then !--------------->
      if(k0/=nf_double)stop 'err 539 type longitude'
       varstart(1)=tidezoom_istr ; varcount(1)=tide_imax
       varstart(2)=tidezoom_jstr ; varcount(2)=tide_jmax
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:2)     &
                                 ,varcount(1:2)     &
                                 ,multi2d_lon(1:tide_imax,1:tide_jmax))
      endif               !--------------->
      if(londim_==1) then !--------------->
       if(k0/=nf_double.and.k0/=nf_real)stop 'Err 549 lon type'
       varstart(1)=tidezoom_istr ; varcount(1)=tide_imax
       if(k0==nf_double)                            &
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:1)     &
                                 ,varcount(1:1)     &
                                 ,multi2d_lon(1:tide_imax,1))
       if(k0==nf_real) then !rrrrr>
          allocate(lonr4_(tide_imax))
           status=nf_get_vara_real(ncid1             &
                                  ,var_id            &
                                  ,varstart(1:1)     &
                                  ,varcount(1:1)     &
                                  ,lonr4_(1:tide_imax))
        multi2d_lon(1:tide_imax,1)=lonr4_(1:tide_imax)
        deallocate(lonr4_)
       endif                !rrrrr>
       
       do j=1,tide_jmax
       do i=1,tide_imax
        multi2d_lon(i,j)=multi2d_lon(i,1)
       enddo
       enddo
      endif               !--------------->
      if(londim_<1.or.londim_>2)stop 'erreur dim lon tide_netcdf_zoom'
      if(status/=0) then !------>
        write(6,*)'erreur lecture longitude 2 londim_=',londim_ &
                  ,'tidezoom_istr=',tidezoom_istr &
                  ,'tide_imax=',tide_imax         &
                  ,'tidefull_imax=',tidefull_imax &
                  ,'k0=',k0                       &
                  ,'rank',par%rank
        stop 'stop dans tide_netcdf_read_h'
      endif              !------>

       do j=1,tide_jmax
       do i=1,tide_imax
        if(multi2d_lon(i,j)>180.)multi2d_lon(i,j)=multi2d_lon(i,j)-360.
       enddo
       enddo

                   status=nf_inq_varid(ncid1,'lat_t',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude_t',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lat_ZHL',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lat',var_id)
      if(status/=0)stop 'erreur nom latitude'


      status=nf_inq_var(ncid1,var_id,bidon_name_,k0,latdim_,tabdim_,i4)

      if(latdim_==2) then !--------------->
      if(k0/=nf_double)stop 'err 597 type latitude'
       varstart(1)=tidezoom_istr ; varcount(1)=tide_imax
       varstart(2)=tidezoom_jstr ; varcount(2)=tide_jmax
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:2)     &
                                 ,varcount(1:2)     &
                                 ,multi2d_lat(1:tide_imax,1:tide_jmax))
      endif               !--------------->
      if(latdim_==1) then !--------------->
       if(k0/=nf_double.and.k0/=nf_real)stop 'Err 607 lat type'
       varstart(1)=tidezoom_jstr ; varcount(1)=tide_jmax
       if(k0==nf_double)                            &
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:1)     &
                                 ,varcount(1:1)     &
                                 ,multi2d_lat(1,1:tide_jmax))
       if(k0==nf_real) then !rrrrr>
          allocate(latr4_(tide_jmax))
           status=nf_get_vara_real(ncid1             &
                                  ,var_id            &
                                  ,varstart(1:1)     &
                                  ,varcount(1:1)     &
                                  ,latr4_(1:tide_jmax))
        multi2d_lat(1,1:tide_jmax)=latr4_(1:tide_jmax)
        deallocate(latr4_)
       endif                !rrrrr>
       do j=1,tide_jmax
       do i=1,tide_imax
        multi2d_lat(i,j)=multi2d_lat(1,j)
       enddo
       enddo
      endif               !--------------->
      if(latdim_<1.or.latdim_>2)stop 'erreur dim lat tide_netcdf_zoom'
      if(status/=0)stop 'erreur lecture latitude'

! correspondances pour passer de la grille du modele à la grille des données:
      if(ktide==1)call tide_hr_to_lr_grid(1)

! lire amplitude:
      if(ki1==1) then !--SSH--> !11-02-21
                    status=nf_inq_varid(ncid1,'Ha',var_id)          !11-02-21
       if(status/=0)status=nf_inq_varid(ncid1,'elevation_a',var_id) !11-02-21
      endif           !--SSH--> !11-02-21
      if(ki1==2) then !--LSA--> !11-02-21
                    status=nf_inq_varid(ncid1,'LSAa',var_id)          !11-02-21
       if(status/=0)status=nf_inq_varid(ncid1,'SAL_amplitude',var_id) !11-02-21
      endif           !--LSA--> !11-02-21
      if(status/=0) stop 'erreur nom amplitude'
      varstart(2)=tidezoom_jstr ; varcount(2)=tide_jmax
      varstart(1)=tidezoom_istr ; varcount(1)=tide_imax
      status=nf_inq_vartype(ncid1,var_id,var_nftype) !12-03-17
      if(var_nftype/=nf_real)stop 'Err 593 tide var type not recognized'
      status=nf_get_vara_real(ncid1,var_id,varstart(1:2)       &
                                          ,varcount(1:2)       &
                                          ,multi2d_var(:,:,1))
      if(status/=0)stop 'erreur lecture amplitude ssh ou potentiel'

! lire retard de phase:
      if(ki1==1) then !--SSH--> !11-02-21
                    status=nf_inq_varid(ncid1,'Hg',var_id)          !11-02-21
       if(status/=0)status=nf_inq_varid(ncid1,'elevation_G',var_id) !11-02-21
      endif           !--SSH--> !11-02-21
      if(ki1==2) then !--LSA--> !11-02-21
                    status=nf_inq_varid(ncid1,'LSAg',var_id)        !11-02-21
       if(status/=0)status=nf_inq_varid(ncid1,'SAL_phase',var_id)   !11-02-21
      endif           !--LSA--> !11-02-21
      if(status/=0)stop 'erreur nom phase'
      varstart(2)=tidezoom_jstr ; varcount(2)=tide_jmax
      varstart(1)=tidezoom_istr ; varcount(1)=tide_imax
      status=nf_inq_vartype(ncid1,var_id,var_nftype) !12-03-17
      if(var_nftype/=nf_real)stop 'Err 606 tide var type not recognized'
      status=nf_get_vara_real(ncid1,var_id,varstart(1:2)    &
                                          ,varcount(1:2)    &
                                          ,multi2d_var(:,:,2))

      if(status/=0)stop 'erreur lecture retard phase ssh ou potentiel'
      call initial_tide_add_phase_offset !25-09-16

! lire valeur par defaut quand la grille est masquée
                   status=nf_get_att_real(ncid1,var_id,'missing_value',filval)
      if(status/=0)status=nf_get_att_real(ncid1,var_id,'_FillValue',filval)
!     if(status/=0)stop 'err 665 valeur de masquage'
! Il est possible que cet attribut n'existe pas si la notion de land sea mask
! n'est pas justifiee pour le champ considere et dans ce cas filval prend la
! valeur par defaut -9999.
      if(status/=0)filval=-9999. !11-02-21

!$ Find  values before interpolation:
!     if(ktide==1.and.used_unused_dom==1)call tide_find_missing_value(1)!10-05-11
      if(ktide==1)call tide_find_missing_value(1)!26-01-20

! Transformer les amplitude/retard_de_phase en amplitudes complexes:
      x0=pi/180.
      do j=1,tide_jmax
      do i=1,tide_imax
       x1=multi2d_var(i,j,1)
       x2=multi2d_var(i,j,2)*x0
          multi2d_var(i,j,1)=x1*cos(x2)
          multi2d_var(i,j,2)=x1*sin(x2)
      enddo
      enddo

!$ fill missing values with nearest values:
      if(par%rank==0)write(6,*)'call tide_fill_value'
      call tide_fill_value(1)

!interpoler:
      if(par%rank==0)write(6,*)'interpolation debut'
      open(unit=3,file=trim(tmpdirname)//'tide_hr_to_lr_grid_1_'//dom_c//'.out')
      do j=0,jmax+1
      do i=0,imax+1
       read(3,*)i0,j0,deci,decj,x0
       if(i0/=i)stop 'erreur fichier tide_hr_to_lr_grid_1'
       if(j0/=j)stop 'erreur fichier tide_hr_to_lr_grid_1'

       i1=min0(max0(int(deci),2),tide_imax-2)
       j1=min0(max0(int(decj),2),tide_jmax-2)
       call tide_interp_poly3(1) ! oui ca coute cher un appel dans une boucle!
       xy_t(i,j,1)=x1
       call tide_interp_poly3(2) ! oui ca coute cher un appel dans une boucle!
       xy_t(i,j,2)=x1

      enddo
      enddo
      close(3)

! Transformer les amplitudes complexes en amplitude/retard_de_phase (d'où le signe -):
      if(ki1==1)then !>>>
       do j=0,jmax+1
       do i=0,imax+1
!       amptide_z(i,j,ktide)=  sqrt(xy_z(i,j,2)**2+xy_z(i,j,1)**2)
!       phitide_z(i,j,ktide)=-atan2(xy_z(i,j,2)  , xy_z(i,j,1)   )
        sshtidecos_w(i,j,ktide)=xy_t(i,j,1)
        sshtidesin_w(i,j,ktide)=xy_t(i,j,2)
       enddo
       enddo
      else           !>>>
       do j=0,jmax+1
       do i=0,imax+1
!       apotide_z(i,j,ktide)=  sqrt(xy_z(i,j,2)**2+xy_z(i,j,1)**2)
!       ppotide_z(i,j,ktide)=-atan2(xy_z(i,j,2)  , xy_z(i,j,1)   )
        potidecos_w(i,j,ktide)=xy_t(i,j,1)
        potidesin_w(i,j,ktide)=xy_t(i,j,2)
       enddo
       enddo
      endif          !>>>

!     write(*,*)'APRES',amptide_z(135,285,ktide) &
!                     ,-phitide_z(135,285,ktide)*180./pi

      status=nf_close(ncid1)

      call allocate_forcages(2,7,0,0,0) ! arg1=desallouer arg2=tide

!     write(6,*)par%rank,'fin tide_netcdf_read_h'

      end subroutine tide_netcdf_read_h

!__________________________________________________________________________________

      subroutine tide_netcdf_read_uv
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
      include 'netcdf.inc'
      double precision dlon_di_
      character*30 uname_amplitude_,uname_phase_    &
                  ,vname_amplitude_,vname_phase_,bidon_name_
      integer status1_,status2_,status3_,status4_,file1orfile2_ &
             ,londim_,latdim_,tabdim_(4)
#ifdef synopsis
       subroutinetitle='tide_netcdf_read_uv'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

        uname_amplitude_='Ua' ; uname_phase_='Ug'
        vname_amplitude_='Va' ; vname_phase_='Vg'

       file1orfile2_=-999 ; k0=0 ; k1=0
  379  k0=k0+1
       status=nf_open(nametide(ktide,k0),nf_nowrite,ncid1)
                       status1_=nf_inq_varid(ncid1,'Ua',var_id)
        if(status1_/=0)status1_=nf_inq_varid(ncid1,'eastward_a',var_id)
                       status2_=nf_inq_varid(ncid1,'Ug',var_id)
        if(status2_/=0)status2_=nf_inq_varid(ncid1,'eastward_G',var_id)
                       status3_=nf_inq_varid(ncid1,'Va',var_id)
        if(status3_/=0)status3_=nf_inq_varid(ncid1,'northward_a',var_id)
                       status4_=nf_inq_varid(ncid1,'Vg',var_id)
        if(status4_/=0)status4_=nf_inq_varid(ncid1,'northward_G',var_id)
        if(status1_==0.and.status2_==0.and.               &
           status3_==0.and.status4_==0)file1orfile2_=k0
       status=nf_close(ncid1)
       if(file1orfile2_==-999.and.k0==1) goto 379
       if(file1orfile2_==-999.and.k0==2) then !!!!!!!!!!>
        if(k1==1) then !1111111111>
        write(6,'(a,a,a,a,a)')'Tidal variables:',trim(uname_amplitude_) &
                                                ,trim(vname_amplitude_) &
                                                ,trim(uname_phase_)     &
                                                ,trim(vname_phase_)
        write(6,*)'not found in the specified following files:'
        write(6,'(a)')trim(nametide(ktide,1))
        write(6,'(a)')trim(nametide(ktide,2))
        stop 'STOP2 dans tide_netcdf_read_uv'
        endif          !1111111111>
! si noms FES ne marchent pas essayer noms S:
        uname_amplitude_='OiVa' ; uname_phase_='OiVg'
        vname_amplitude_='OjVa' ; vname_phase_='OjVg'
        k0=0 ; k1=k1+1 ; goto 379
       endif                                  !!!!!!!!!!>


! lire l'identifiant du fichier:
      if(par%rank==0) then !>>>>>
       write(6,*)'dans tide_netcdf_read_uv je lis:'
       write(6,'(a)')nametide(ktide,file1orfile2_)
      endif                !>>>>>

      status=nf_open(nametide(ktide,file1orfile2_),nf_nowrite,ncid1)
      if(status/=0)then
         stop ' Err 795 fichier maree absent'
         status=nf_close(ncid1)
         return
      endif

! Interpoler composante u: debut

! lire l'identifiant des dimensions:
                   status=nf_inq_dimid(ncid1,'x_XHL',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'imax_u',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'ni_u',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'x',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'nx',dim_x_id)
      if(status/=0)stop 'erreur nom pour dim_x_id'
                   status=nf_inq_dimid(ncid1,'y_XHL',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'jmax_u',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'nj_u',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'y',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'ny',dim_y_id)
      if(status/=0)stop 'erreur nom pour dim_y_id'

! Lecture des valeurs des dimensions
      status=nf_inq_dimlen(ncid1,dim_x_id,tide_imax)
      status=nf_inq_dimlen(ncid1,dim_y_id,tide_jmax)

! Garder en memoire la taille totale du domaine avant redecoupage:
      tidefull_imax=tide_imax
      tidefull_jmax=tide_jmax
      call tide_netcdf_zoom(2)

      call allocate_forcages(1,7,tide_imax,tide_jmax,2) ! arg1=allouer arg2=tide !09-10-09

! lire longitude et latitude vitesse Ox:
                   status=nf_inq_varid(ncid1,'lon_XHL',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lon_u',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'longitude_u',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lon',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'longitude',var_id)
      if(status/=0)stop 'erreur nom pour longitude'

      status=nf_inq_var(ncid1,var_id,bidon_name_,k0,londim_,tabdim_,i4)
      if(k0/=nf_double)stop 'err 833 type longitude'

      if(londim_==2) then !--------------->
       varstart(1)=tidezoom_istr ; varcount(1)=tide_imax
       varstart(2)=tidezoom_jstr ; varcount(2)=tide_jmax
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:2)     &
                                 ,varcount(1:2)     &
                                 ,multi2d_lon(1:tide_imax,1:tide_jmax))
      endif               !--------------->
      if(londim_==1) then !--------------->
       varstart(1)=tidezoom_istr ; varcount(1)=tide_imax
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:1)     &
                                 ,varcount(1:1)     &
                                 ,multi2d_lon(1:tide_imax,1))
       do j=1,tide_jmax
       do i=1,tide_imax
        multi2d_lon(i,j)=multi2d_lon(i,1)
       enddo
       enddo
      endif               !--------------->
      if(londim_<1.or.londim_>2)stop 'erreur dim lon tide_netcdf_zoom'
      if(status/=0)stop 'erreur lecture longitude 3'

                   status=nf_inq_varid(ncid1,'lat_XHL',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lat_u',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude_u',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lat',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude',var_id)
      if(status/=0)stop 'erreur nom pour latitude'

      status=nf_inq_var(ncid1,var_id,bidon_name_,k0,latdim_,tabdim_,i4)
      if(k0/=nf_double)stop 'erreur type latitude'

      if(latdim_==2) then !--------------->
       varstart(1)=tidezoom_istr ; varcount(1)=tide_imax
       varstart(2)=tidezoom_jstr ; varcount(2)=tide_jmax
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:2)     &
                                 ,varcount(1:2)     &
                                 ,multi2d_lat(1:tide_imax,1:tide_jmax))
      endif               !--------------->
      if(latdim_==1) then !--------------->
       varstart(1)=tidezoom_jstr ; varcount(1)=tide_jmax
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:1)     &
                                 ,varcount(1:1)     &
                                 ,multi2d_lat(1,1:tide_jmax))
       do j=1,tide_jmax
       do i=1,tide_imax
        multi2d_lat(i,j)=multi2d_lat(1,j)
       enddo
       enddo
      endif               !--------------->
      if(latdim_<1.or.latdim_>2)stop 'erreur dim lat tide_netcdf_zoom'
      if(status/=0)stop 'erreur lecture latitude'

! correspondances pour passer de la grille du modele à la grille des données:
      if(ktide==1)call tide_hr_to_lr_grid(2)

! lire l'amplitude du courant Ox:
                   status=nf_inq_varid(ncid1,'Ua',var_id)         !11-02-21
      if(status/=0)status=nf_inq_varid(ncid1,'OiVa',var_id)       !11-02-21
      if(status/=0)status=nf_inq_varid(ncid1,'eastward_a',var_id) !11-02-21
      if(status/=0)stop 'erreur nom amplitude courant de maree Oi'
      varstart(2)=tidezoom_jstr ; varcount(2)=tide_jmax
      varstart(1)=tidezoom_istr ; varcount(1)=tide_imax
      status=nf_inq_vartype(ncid1,var_id,var_nftype) !12-03-17
      if(var_nftype/=nf_real)stop 'Err 854 tide var type not recognized'
      status=nf_get_vara_real(ncid1,var_id,varstart(1:2)       &
                                          ,varcount(1:2)       &
                                          ,multi2d_var(:,:,1))
      if(status/=0)stop 'erreur lecture amplitude U maree'

! lire le retard de phase du courant Ox:
                   status=nf_inq_varid(ncid1,'Ug',var_id)         !11-02-21
      if(status/=0)status=nf_inq_varid(ncid1,'OiVg',var_id)       !11-02-21
      if(status/=0)status=nf_inq_varid(ncid1,'eastward_G',var_id) !11-02-21
      if(status/=0)stop 'erreur nom retard phase courant de maree Oi'
      varstart(2)=tidezoom_jstr ; varcount(2)=tide_jmax
      varstart(1)=tidezoom_istr ; varcount(1)=tide_imax
      status=nf_inq_vartype(ncid1,var_id,var_nftype) !12-03-17
      if(var_nftype/=nf_real)stop 'Err 866 tide var type not recognized'
      status=nf_get_vara_real(ncid1,var_id,varstart(1:2)       &
                                          ,varcount(1:2)       &
                                          ,multi2d_var(:,:,2))
      if(status/=0)stop 'erreur lecture phase U maree'
      call initial_tide_add_phase_offset !25-09-16

! lire la valeur par defaut quand la grille est masquée
                   status=nf_get_att_real(ncid1,var_id,'missing_value',filval)
      if(status/=0)status=nf_get_att_real(ncid1,var_id,'_FillValue',filval)
      if(status/=0)stop 'err 925 sur valeur de masquage'

!$ Find missing values before interpolation:
!     if(ktide==1.and.used_unused_dom==1)call tide_find_missing_value(2)!10-05-11
      if(ktide==1)call tide_find_missing_value(2)!26-01-20

! Transformer les amplitude/retard_de_phase en amplitudes complexes:
      x0=pi/180.
      do j=1,tide_jmax
      do i=1,tide_imax
       x1=multi2d_var(i,j,1)
       x2=multi2d_var(i,j,2)*x0
          multi2d_var(i,j,1)=x1*cos(x2)
          multi2d_var(i,j,2)=x1*sin(x2)
      enddo
      enddo

!$ fill missing values with nearest values:
      if(par%rank==0)write(6,*)'tide_fill_value'
      call tide_fill_value(2)

! interpolation:
      if(par%rank==0)write(6,*)'interpolation debut'
      open(unit=3,file=trim(tmpdirname)//'tide_hr_to_lr_grid_2_'//dom_c//'.out')
      do j=0,jmax+1
      do i=0,imax+1
       read(3,*)i0,j0,deci,decj,x0
       if(i0/=i)stop 'erreur fichier tide_hr_to_lr_grid_1'
       if(j0/=j)stop 'erreur fichier tide_hr_to_lr_grid_1'
       i1=min0(max0(int(deci),2),tide_imax-2)
       j1=min0(max0(int(decj),2),tide_jmax-2)
       call tide_interp_poly3(1) ! oui ca coute cher un appel dans une boucle!
       xy_t(i,j,1)=x1   ! u1*cos(wt)
       call tide_interp_poly3(2) ! oui ca coute cher un appel dans une boucle!
       xy_t(i,j,2)=x1   ! u2*sin(wt)
      enddo
      enddo
      close(3)

      call allocate_forcages(2,7,0,0,0) ! arg1=allouer arg2=tide !09-10-09

! Interpoler composante u: Fin

! Interpoler composante v: debut

! lire l'identifiant des dimensions:
                   status=nf_inq_dimid(ncid1,'x_YHL',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'imax_v',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'ni_v',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'x',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'nx',dim_x_id)
      if(status/=0)stop 'erreur nom pour dim_x_id'
                   status=nf_inq_dimid(ncid1,'y_YHL',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'jmax_v',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'nj_v',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'y',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'ny',dim_y_id)
      if(status/=0)stop 'erreur nom pour dim_y_id'

! Lecture des valeurs des dimensions
      status=nf_inq_dimlen(ncid1,dim_x_id,tide_imax)
      status=nf_inq_dimlen(ncid1,dim_y_id,tide_jmax)

! Garder en memoire la taille totale du domaine avant redecoupage:
      tidefull_imax=tide_imax
      tidefull_jmax=tide_jmax
      call tide_netcdf_zoom(3)

      call allocate_forcages(1,7,tide_imax,tide_jmax,2) ! arg1=allouer arg2=tide !09-10-09

! lire longitude et latitude vitesse Oy: !11-11-13
                   status=nf_inq_varid(ncid1,'lon_YHL',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lon_v',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'longitude_v',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'longitude',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lon',var_id)
      if(status/=0)stop 'erreur nom pour longitude'

      status=nf_inq_var(ncid1,var_id,bidon_name_,k0,londim_,tabdim_,i4)
      if(k0/=nf_double)stop 'err 1000 type longitude'

      if(londim_==2) then !--------------->
       varstart(1)=tidezoom_istr ; varcount(1)=tide_imax
       varstart(2)=tidezoom_jstr ; varcount(2)=tide_jmax
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:2)     &
                                 ,varcount(1:2)     &
                                 ,multi2d_lon(1:tide_imax,1:tide_jmax))
      endif               !--------------->
      if(londim_==1) then !--------------->
       varstart(1)=tidezoom_istr ; varcount(1)=tide_imax
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:1)     &
                                 ,varcount(1:1)     &
                                 ,multi2d_lon(1:tide_imax,1))
       do j=1,tide_jmax
       do i=1,tide_imax
        multi2d_lon(i,j)=multi2d_lon(i,1)
       enddo
       enddo
      endif               !--------------->
      if(londim_<1.or.londim_>2)stop 'erreur dim lon tide_netcdf_zoom'
      if(status/=0)stop 'erreur lecture longitude 3'

                   status=nf_inq_varid(ncid1,'lat_YHL',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lat_v',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude_v',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lat',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude',var_id)
      if(status/=0)stop 'erreur nom pour latitude'

      status=nf_inq_var(ncid1,var_id,bidon_name_,k0,latdim_,tabdim_,i4)
      if(k0/=nf_double)stop 'erreur type latitude'

      if(latdim_==2) then !--------------->
       varstart(1)=tidezoom_istr ; varcount(1)=tide_imax
       varstart(2)=tidezoom_jstr ; varcount(2)=tide_jmax
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:2)     &
                                 ,varcount(1:2)     &
                                 ,multi2d_lat(1:tide_imax,1:tide_jmax))
      endif               !--------------->
      if(latdim_==1) then !--------------->
       varstart(1)=tidezoom_jstr ; varcount(1)=tide_jmax
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:1)     &
                                 ,varcount(1:1)     &
                                 ,multi2d_lat(1,1:tide_jmax))
       do j=1,tide_jmax
       do i=1,tide_imax
        multi2d_lat(i,j)=multi2d_lat(1,j)
       enddo
       enddo
      endif               !--------------->
      if(latdim_<1.or.latdim_>2)stop 'erreur dim lat tide_netcdf_zoom'
      if(status/=0)stop 'erreur lecture latitude'

! correspondances pour passer de la grille du modele à la grille des données:
      if(ktide==1)call tide_hr_to_lr_grid(3)

! lire l'amplitude du courant Oy:
!       vname_amplitude_='Va' ; vname_phase_='Vg'
!       vname_amplitude_='OjVa' ; vname_phase_='OjVg'
                   status=nf_inq_varid(ncid1,'Va',var_id)          !11-02-21
      if(status/=0)status=nf_inq_varid(ncid1,'OjVa',var_id)        !11-02-21
      if(status/=0)status=nf_inq_varid(ncid1,'northward_a',var_id) !11-02-21
      if(status/=0)stop 'erreur nom amplitude V maree'
      varstart(2)=tidezoom_jstr ; varcount(2)=tide_jmax
      varstart(1)=tidezoom_istr ; varcount(1)=tide_imax
      status=nf_inq_vartype(ncid1,var_id,var_nftype) !12-03-17
      if(var_nftype/=nf_real)stop 'Err1021 tide var type not recognized'
      status=nf_get_vara_real(ncid1,var_id,varstart(1:2)       &
                                          ,varcount(1:2)       &
                                          ,multi2d_var(:,:,1))
      if(status/=0)stop 'erreur lecture amplitude V maree'

! lire le retard de phase du courant Oy:
                   status=nf_inq_varid(ncid1,'Vg',var_id)          !11-02-21
      if(status/=0)status=nf_inq_varid(ncid1,'OjVg',var_id)        !11-02-21
      if(status/=0)status=nf_inq_varid(ncid1,'northward_G',var_id) !11-02-21
      if(status/=0)stop 'erreur nom phase V maree'
      varstart(2)=tidezoom_jstr ; varcount(2)=tide_jmax
      varstart(1)=tidezoom_istr ; varcount(1)=tide_imax
      status=nf_inq_vartype(ncid1,var_id,var_nftype) !12-03-17
      if(var_nftype/=nf_real)stop 'Err1033 tide var type not recognized'
      status=nf_get_vara_real(ncid1,var_id,varstart(1:2)       &
                                          ,varcount(1:2)       &
                                          ,multi2d_var(:,:,2))
      call initial_tide_add_phase_offset !25-09-16

! lire la valeur par defaut quand la grille est masquée
                   status=nf_get_att_real(ncid1,var_id,'missing_value',filval)
      if(status/=0)status=nf_get_att_real(ncid1,var_id,'_FillValue',filval)
      if(status/=0)stop 'err 1091 sur valeur de masquage'

!$ Find missing values before interpolation:
!     if(ktide==1.and.used_unused_dom==1)call tide_find_missing_value(3)!10-05-11
      if(ktide==1)call tide_find_missing_value(3)!10-05-11

! Transformer les amplitude/retard_de_phase en amplitudes complexes:
      x0=pi/180.
      do j=1,tide_jmax
      do i=1,tide_imax
       x1=multi2d_var(i,j,1)
       x2=multi2d_var(i,j,2)*x0
          multi2d_var(i,j,1)=x1*cos(x2)
          multi2d_var(i,j,2)=x1*sin(x2)
      enddo
      enddo

!$ fill missing values with nearest values:
      if(par%rank==0)write(6,*)'tide_fill_value'
      call tide_fill_value(3)

! interpolation:
      if(par%rank==0)write(6,*)'interpolation debut'
      open(unit=3,file=trim(tmpdirname)//'tide_hr_to_lr_grid_3_'//dom_c//'.out')
      do j=0,jmax+1
      do i=0,imax+1
       read(3,*)i0,j0,deci,decj,x0
       if(i0/=i)stop 'erreur fichier tide_hr_to_lr_grid_1'
       if(j0/=j)stop 'erreur fichier tide_hr_to_lr_grid_1'
       i1=min0(max0(int(deci),2),tide_imax-2)
       j1=min0(max0(int(decj),2),tide_jmax-2)
       call tide_interp_poly3(1) ! oui ca coute cher un appel dans une boucle!
       xy_t(i,j,3)=x1  ! v1*cos(wt)
       call tide_interp_poly3(2) ! oui ca coute cher un appel dans une boucle!
       xy_t(i,j,4)=x1  ! v2*sin(wt)
      enddo
      enddo
      close(3)

! A ce stade les 2 composantes sont interpolées au centre de la grille C
! Il reste à: 1 appliquer la rotation pour aligner les composantes du courant
! selon les axes du modeles, 2 passer du centre de la maille aux facettes,
! 3 revenir des amplitudes complexes aux amplitudes/retard_de_phase

! Rotation
      open(unit=3,file=trim(tmpdirname)//'tide_hr_to_lr_grid_1_'//dom_c//'.out')
      do j=0,jmax+1
      do i=0,imax+1

      read(3,*)i0,j0,deci,decj,x0

      if(i0.ne.i)stop 'erreur1 sur i'
      if(j0.ne.j)stop 'erreur1 sur j'
! x0 rotation modele externe vers symphonie
      dlon_di_=lon_t(i+1,j)-lon_t(i-1,j)      !14-02-13
      if(dlon_di_<-pi)dlon_di_=dlon_di_+2.*pi
      if(dlon_di_> pi)dlon_di_=dlon_di_-2.*pi
      x0=x0*deg2rad                        & ! rotation modele externe
        -atan2(  lat_t(i+1,j)-lat_t(i,j)   & ! rotation symphonie
               ,dlon_di_*cos(lat_t(i,j))  )
!              ,(lon_t(i+1,j)-lon_t(i,j))*cos(lat_t(i,j))  )

       x1=cos(x0)*xy_t(i,j,1)-sin(x0)*xy_t(i,j,3)  !   U1=u1*cos(x0)-v1*sin(x0)  ! cos(wt)
       x3=sin(x0)*xy_t(i,j,1)+cos(x0)*xy_t(i,j,3)  !   V1=u1*sin(x0)+v1*cos(x0)  ! cos(wt)

! nouvelles composantes tournées:
       xy_t(i,j,1)=x1     ! U1 ! cos(wt)
       xy_t(i,j,3)=x3     ! V1 ! cos(wt)

       x2=cos(x0)*xy_t(i,j,2)-sin(x0)*xy_t(i,j,4)  ! U2=u2*cos(x0)-v2*sin(x0) ! sin(wt)
       x4=sin(x0)*xy_t(i,j,2)+cos(x0)*xy_t(i,j,4)  ! V2=u2*sin(x0)+v2*cos(x0) ! sin(wt)

! nouvelles composantes tournées:
       xy_t(i,j,2)=x2    ! U2 ! sin(wt)
       xy_t(i,j,4)=x4    ! V2 ! sin(wt)

      enddo
      enddo
      close(3)                                                          !17-12-09

! Grille C:
! Echange z0 sur xy_t(:,:,1:4)              !12-11-14
      call obc_h_xy_t_z0(1)
      call obc_h_xy_t_z0(2)
      call obc_h_xy_t_z0(3)
      call obc_h_xy_t_z0(4)
      do j=1,jmax+1
      do i=1,imax+1
       veltidecos_u(i,j,ktide)=0.5*(xy_t(i,j,1)+xy_t(i-1,j  ,1))*mask_u(i,j,kmax+1) ! U1_x ! cos(wt)
       veltidesin_u(i,j,ktide)=0.5*(xy_t(i,j,2)+xy_t(i-1,j  ,2))*mask_u(i,j,kmax+1) ! U2_x ! sin(wt)
       veltidecos_v(i,j,ktide)=0.5*(xy_t(i,j,3)+xy_t(i  ,j-1,3))*mask_v(i,j,kmax+1) ! V1_y ! cos(wt)
       veltidesin_v(i,j,ktide)=0.5*(xy_t(i,j,4)+xy_t(i  ,j-1,4))*mask_v(i,j,kmax+1) ! V2_y ! sin(wt)
      enddo
      enddo

      status=nf_close(ncid1)

      call allocate_forcages(2,7,0,0,0) ! desallouer

      end subroutine tide_netcdf_read_uv

!-------------------------------------------------------------------------------

      subroutine tide_find_missing_value(case_)
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
      real trou_
      integer case_

      trou_   =-filval
      if(filval==0)trou_   =-9999.

      if(par%rank==0)write(6,*)'tide_find_missing_value debut ',case_

      if(case_==1) then
      open(unit=3,file=trim(tmpdirname)//'tide_misval_h_'//dom_c//'.out',STATUS='REPLACE')
      open(unit=4,file=trim(tmpdirname)//'tide_hr_to_lr_grid_1_'//dom_c//'.out')
      endif
      if(case_==2)then
      open(unit=3,file=trim(tmpdirname)//'tide_misval_u_'//dom_c//'.out')
      open(unit=4,file=trim(tmpdirname)//'tide_hr_to_lr_grid_2_'//dom_c//'.out')
      endif
      if(case_==3)then
      open(unit=3,file=trim(tmpdirname)//'tide_misval_v_'//dom_c//'.out')
      open(unit=4,file=trim(tmpdirname)//'tide_hr_to_lr_grid_3_'//dom_c//'.out')
      endif

      do j=0,jmax+1
      do i=0,imax+1

       read(4,*)i0,j0,deci,decj,x0

! Il pouvait se produire qu'un proc inutile soit associE A une extraction dans
! la grille de maree elle aussi 100% inutile, mettant en echec le processus de
! recherche de bouchage de l'etape suivante. Depuis on conditionne la necessite
! de bouchage au fait de concerner seulement les points en mer avec le test
! suivant: !20-05-20
       if(mask_t(i,j,kmax)==1) then !pmx> !20-05-20

         if(i0/=i)stop 'erreur fichier tide_hr_to_lr_grid'
         if(j0/=j)stop 'erreur fichier tide_hr_to_lr_grid'
         i1=min0(max0(int(deci),2),tide_imax-2)
         j1=min0(max0(int(decj),2),tide_jmax-2)
         do j2=j1-1,j1+2
         do i2=i1-1,i1+2
          if(multi2d_var(i2,j2,1)==filval)multi2d_var(i2,j2,1)=trou_
         enddo
         enddo

       endif                        !pmx> !20-05-20

      enddo
      enddo
      close(4)


      do 1862 j=1,tide_jmax
      do 1862 i=1,tide_imax

       if(multi2d_var(i,j,1)==trou_   )then !§§§§§§§§§§§§§§§>
        ksecu=0
        i10=1
        dist1=1e10
 1864   continue

        do 1863 k1=0,1

         j0=j-i10
         j2=j+i10
         j0=min(max(j0,1),tide_jmax)      !27-10-11
         j2=min(max(j2,1),tide_jmax)      !27-10-11
         j3=k1*(j2-j0)+(1-k1)

         i0=i-i10+k1
         i2=i+i10-k1
         i0=min(max(i0,1),tide_imax)      !27-10-11
         i2=min(max(i2,1),tide_imax)      !27-10-11
         i3=(1-k1)*(i2-i0)+k1

         do 1861 j1=j0,j2,j3
         do 1861 i1=i0,i2,i3

! commentée le 25-10-11
!        if(i1<1        )stop 'faute tide 5'  !01-07-10
!        if(i1>tide_imax)stop 'faute tide 6'
!        if(j1<1        )stop 'faute tide 7'
!        if(j1>tide_jmax)stop 'faute tide 8'

! si le modele plante sur fautes 5 6 7 ou 8 c'est qu'il n'y a pas de bouchage
! suffisament proche du trou et que l'algo cherche au delà de la zone d'extraction,
! ce qu'on ne permet pas pour avoir une parallelisation parfaite.
! Ce qu'on peu faire: on peut augmenter la taille de la
! zone d'extraction en augmentant i0 (repere 1233). Mais cela revele surtout
! que les masques oceano et tide sont trop differents et qu'il faut sans doute
! ameliorer le masque oceano

          if(multi2d_var(i1,j1,1)/=filval.and.          &
             multi2d_var(i1,j1,1)/=trou_   )then !%%%%%%%%%%%%%>
           dist2=sqrt(real(i-i1)**2+real(j-j1)**2)
           ksecu=1
             if(dist2.lt.dist1) then                    !>>>>>>>>>>>>>
              i4=i1
              j4=j1
              dist1=dist2
             endif                                      !>>>>>>>>>>>>>
          endif                                !%%%%%%%%%%%%%>
 1861    continue

 1863   continue
        i10=i10+1
        if(i10>999)  then !m°v°m>

! Le point i,j de la grille de maree n'a pas de boucheur? 
! Produire un rapport d'erreur

          write(10+par%rank,*)'PROBLEME: forfait bouchage trou dEpassE' !12-10-19
          write(10+par%rank,'(a,a)')'dom_c:',trim(dom_c)
          write(10+par%rank,*)'par%rank=',par%rank
          write(10+par%rank,*)'Le point i,j= ',i,j &
                             ,' de la grille de marEe n''a pas trouvE' &
                             ,' de bouchage proche'
          write(10+par%rank,*)'tidezoom_istr,tidezoom_jstr',tidezoom_istr,tidezoom_jstr
          write(10+par%rank,*)'tide_imax,tide_jmax',tide_imax,tide_jmax
          write(10+par%rank,*)'i+tidezoom_istri-1,j+tidezoom_jstr-1',i+tidezoom_istr-1,j+tidezoom_jstr-1
          write(10+par%rank,*)'Lon,Lat de ce point',multi2d_lon(i,j),multi2d_lat(i,j)
          write(10+par%rank,*)'Dans subroutine tide_find_missing_value' &
         ,' case_=',case_

           write(10+par%rank,*)'Liste des points qui aurait pu etre choisis:'
           do j5=1,tide_jmax ; do i5=1,tide_imax
              if(multi2d_var(i5,j5,1)/=filval.and.   &
                 multi2d_var(i5,j5,1)/=trou_   )     &
            write(10+par%rank,*)i5,j5,multi2d_var(i5,j5,1)
           enddo            ; enddo
           write(10+par%rank,*)'Si cette liste est vide, cela signifie que' &
           ,' l extraction dans la grille du modele de maree encadrant le'  &
           ,' le sousdomaine mpi est entierement masquee.'

          flag_stop=1 ; goto 1250
        endif             !m°v°m>
        if(ksecu.eq.0)goto 1864
                 write(3,'(4(i4,1x))')i,j  ,i4,j4

       endif                           !§§§§§§§§§§§§§§§>
 1862 continue
 1250 close(3)

      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) ;
      if(k0/=0) stop &
      'STOP tide_find_missing_value i10>999 See fort.xxx files' !25-10-11

      end subroutine tide_find_missing_value

!-------------------------------------------------------------------------------

      subroutine tide_fill_value(ki1)
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
      integer ki1

      if(ki1==1) then
      open(unit=3,file=trim(tmpdirname)//'tide_misval_h_'//dom_c//'.out')
      endif
      if(ki1==2)then
      open(unit=3,file=trim(tmpdirname)//'tide_misval_u_'//dom_c//'.out')
     endif
      if(ki1==3)then
      open(unit=3,file=trim(tmpdirname)//'tide_misval_v_'//dom_c//'.out')
      endif

  931 read(3,*,end=928)i,j,i1,j1
      multi2d_var(i,j,1)=multi2d_var(i1,j1,1)
      multi2d_var(i,j,2)=multi2d_var(i1,j1,2)
      goto 931

  928 close(3)


      end subroutine tide_fill_value

!------------------------------------------------------------------------------

      subroutine tide_hr_to_lr_grid(ki1)
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
      integer ki1
      double precision deci_(0:1,0:1),decj_(0:1,0:1)    &
                        ,dy_(0:1,0:1),  dx_(0:1,0:1)    &
                        ,dlon_di_,dlon_dj_,dlon_dm_
#ifdef synopsis
       subroutinetitle='tide_find_missing_value'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(ki1==1)                                                        &
      open(unit=3,file=trim(tmpdirname)//'tide_hr_to_lr_grid_1_'//dom_c//'.out')
      if(ki1==2)                                                        &
      open(unit=3,file=trim(tmpdirname)//'tide_hr_to_lr_grid_2_'//dom_c//'.out')
      if(ki1==3)                                                        &
      open(unit=3,file=trim(tmpdirname)//'tide_hr_to_lr_grid_3_'//dom_c//'.out')

      if(par%rank==0)write(6,*)'tide_hr_to_lr_grid',ki1

      do j=0,jmax+1                                                     !18-10-09
      do i=0,imax+1

! First guess. Dans quel grand secteur se trouve le point?
! Le point de chaque début de ligne fait l'objet d'un pretraitement
! visant à acceler la convergence de l'algo general qui peut se
! trouver pris en defaut par une grille non periodique ou tres
! non lineaire. Ensuite le first guess des points suivants correpond
! au resultat precedent.
      if(i==0) then !0000000000000000000>
       dist1=1.d20
       do j1=3,tide_jmax-2,(tide_jmax-5)/5
       do i1=3,tide_imax-2,(tide_imax-5)/5
        call lonlat2distance(lon_t(i,j),lat_t(i,j)       &
                            ,multi2d_lon(i1,j1)*deg2rad  &
                            ,multi2d_lat(i1,j1)*deg2rad  &
                            ,dist2)
        if(dist2<dist1) then
         dist1=dist2 ; x2=i1 ; x3=j1
        endif
       enddo
       enddo
      deci=x2
      decj=x3
      endif         !0000000000000000000>

! ETAPE 1: trouver les coordonnées dans la grille ORCA:

! First guess: centre du domaine:
      k10=0
 1456 continue

! Principe suppose une relation lineaire entre lat lon et indice de grille.
!      dlon/di*Di+dlon/dj*Dj=Dlon
!      dlat/di*Di+dlat/dj*Dj=Dlat
! On cherche Di et Dj correspondant à Dlon=lon(i,j)-lonmeteo(i0,j0)
!                                et à Dlat=lat(i,j)-latmeteo(i0,j0)

! Algo special grille globale fes2012 dont la periodicite semble etre lon(1)=lon(tide_imax)
      if(deci>tide_imax)deci=deci-tide_imax+1
      if(deci<1        )deci=deci+tide_imax-1
      decj=min(max(decj,2.),tide_jmax-2.)
      i1=int(deci)
      j1=int(decj)

      do j2=0,1
      do i2=0,1

      i0=i1+i2
      j0=j1+j2
      im1=i0-1
      ip1=i0+1

! Algo special grille globale fes2012 dont la periodicite semble etre lon(1)=lon(tide_imax) !19-04-15
      if(i0 >tide_imax)i0= i0 -tide_imax+1
      if(im1>tide_imax)im1=im1-tide_imax+1
      if(ip1>tide_imax)ip1=ip1-tide_imax+1
      if(i0 <1        )i0= i0 +tide_imax-1
      if(im1<1        )im1=im1+tide_imax-1
      if(ip1<1        )ip1=ip1+tide_imax-1

      dlon_di_= multi2d_lon(ip1 ,j0  )- multi2d_lon(im1 ,j0  )
      dlon_dj_= multi2d_lon(i0  ,j0+1)- multi2d_lon(i0  ,j0-1)
      dlon_dm_=rad2deg*lon_t(i  ,j)   - multi2d_lon(i0  ,j0)

      if(dlon_di_<-180.)dlon_di_=dlon_di_+360.
      if(dlon_di_> 180.)dlon_di_=dlon_di_-360.
      if(dlon_dj_<-180.)dlon_dj_=dlon_dj_+360.
      if(dlon_dj_> 180.)dlon_dj_=dlon_dj_-360.
      if(dlon_dm_<-180.)dlon_dm_=dlon_dm_+360.
      if(dlon_dm_> 180.)dlon_dm_=dlon_dm_-360.

! Determinant principal:
      x1=( dlon_di_                                             &
          *( multi2d_lat(i0  ,j0+1)- multi2d_lat(i0  ,j0-1))    &
          -( multi2d_lat(ip1 ,j0  )- multi2d_lat(im1 ,j0  ))    &
          *dlon_dj_ )*0.25

      deci_(i2,j2)=                                             &
       i0+( dlon_dm_                                            &
          *( multi2d_lat(i0  ,j0+1)- multi2d_lat(i0,j0-1))      &
          -(rad2deg*lat_t(i,j)     - multi2d_lat(i0,j0))        &
          *dlon_dj_)/x1*0.5

      decj_(i2,j2)=                                             &
       j0+( dlon_di_                                            &
          *(rad2deg*lat_t(i,j)     - multi2d_lat(i0  ,j0))      &
          -( multi2d_lat(ip1 ,j0  )- multi2d_lat(im1 ,j0))      &
          *dlon_dm_   )/x1*0.5

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


! Ces lignes prevoit le fait que la grille externe peut
! etre periodique selon Oi (par exemple en cas de domaine circulaire)
!     if(deci>tide_imax-1.0001d0)deci=deci+4-tide_imax
!     if(deci<          2.0001d0)deci=deci+tide_imax-4
! Bornes pour eviter depassement de memoire:
!     deci=min(max(deci,2.0001d0),tide_imax-1.0001d0)

! Si le point visé est different du first guess refaire le calcul
! avec un first guess donné par le dernier point visé:
      if(sqrt( (deci-x2)**2+(decj-x3)**2 ).gt.0.001)then
       x2=deci
       x3=decj
       k10=k10+1
       if(k10>20)then !!!!!!>
         write(6,*)'PB convergence hr_to_lr tide' &
                  ,'par%rank i j deci decj='      &
                  ,par%rank,i,j,deci,decj         &
       ,'Indices globaux ',i+par%timax(1),j+par%tjmax(1) &
       ,'tide_imax tide_jmax',tide_imax,tide_jmax &
                  , 'ki1=',ki1
         stop 'hr_to_lr ne converge pas dans le forfait'
!        write(6,*)'ATTENTION BIDOUILLE DANS INITIAL_TIDE'
!        GOTO 1266 ! BIDOUILLE
       endif          !!!!!!>
       goto 1456
      endif
!1266 CONTINUE ! BIDOUILLE

! ETAPE 2: CALCULER L'ANGLE D'ORIENTATION LOCALE DE LA GRILLE FES2012
      if(deci>tide_imax)deci=deci-tide_imax+1
      if(deci<1        )deci=deci+tide_imax-1
      decj=min(max(decj,2.),tide_jmax-2.)

       i1=int(deci) ; j1=int(decj)
       rapi=deci-i1 ; rapj=decj-j1
       do j2=0,1
       do i2=0,1
        i0=i1+i2
        j0=j1+j2
! Algo special grille globale fes2012 dont la periodicite semble etre lon(1)=lon(tide_imax)
        if(i0 >tide_imax)i0= i0 -tide_imax+1
        if(im1>tide_imax)im1=im1-tide_imax+1
        if(ip1>tide_imax)ip1=ip1-tide_imax+1
        if(i0 <1        )i0= i0 +tide_imax-1
        if(im1<1        )im1=im1+tide_imax-1
        if(ip1<1        )ip1=ip1+tide_imax-1

        dlon_di_= multi2d_lon(ip1 ,j0  )- multi2d_lon(im1 ,j0  )
        if(dlon_di_<-180.)dlon_di_=dlon_di_+360.
        if(dlon_di_> 180.)dlon_di_=dlon_di_-360.

        dy_(i2,j2)=  multi2d_lat(ip1 ,j0)- multi2d_lat(im1 ,j0)
        dx_(i2,j2)=dlon_di_           &
                             *cos( multi2d_lat(i0,j0)*deg2rad)
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

       write(3,'(2(i4,1x),3(f9.4,1x))')i,j,deci,decj     &
             ,atan2(y1,x1)*rad2deg*tide_flagrotation  ! Si convention 'earthaxis' tide_flagrotation=0

      enddo
      enddo

      close(3)

      return
      end subroutine tide_hr_to_lr_grid

!------------------------------------------------------------------------------

      subroutine tide_interp_poly3(ki1)
      use module_principal
      use module_forcages
      implicit none
      double precision,dimension(4,4):: y_
      double precision,dimension(4)::   x_
      double precision,dimension(4)::   val_
      integer ichoix,ki1
#ifdef synopsis
       subroutinetitle='tide_interp_poly3'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

        x_(1)=i1-1.
        x_(2)=i1
        x_(3)=i1+1.
        x_(4)=i1+2.

        do j2=1,4
        j0=j1-2+j2


        y_(1,1)=multi2d_var(i1-1,j0,ki1)
        y_(2,1)=multi2d_var(i1  ,j0,ki1)
        y_(3,1)=multi2d_var(i1+1,j0,ki1)
        y_(4,1)=multi2d_var(i1+2,j0,ki1)

          do j3=2,4
           do i3=1,5-j3
            y_(i3,j3)=(y_(i3+1,j3-1)-y_(i3,j3-1))              &
                     /(x_(i3+j3-1)  -x_(i3))
           enddo
          enddo

! val_    est l'interpolation polynome p3 direction x
       val_(j2)=y_(1,4)
       do i3=2,4
        val_(j2)=val_(j2)*(deci-x_(i3))+y_(i3,5-i3)
       enddo

       enddo ! fin de boucle sur j2

! interpolation direction j:
        x_(1)=j1-1.
        x_(2)=j1
        x_(3)=j1+1.
        x_(4)=j1+2.
        y_(1,1)=val_(1)
        y_(2,1)=val_(2)
        y_(3,1)=val_(3)
        y_(4,1)=val_(4)
          do j3=2,4
           do i3=1,5-j3
            y_(i3,j3)=(y_(i3+1,j3-1)-y_(i3,j3-1))              &
                     /(x_(i3+j3-1)  -x_(i3))
           enddo
          enddo
        x1=y_(1,4)
        do i3=2,4
         x1=x1*(decj-x_(i3))+y_(i3,5-i3)
        enddo

      end subroutine tide_interp_poly3

!----------------------------------------------------------------------

      subroutine tide_netcdf_zoom(ichoix)
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
      include 'netcdf.inc'
      integer ichoix,tabdim_(4),londim_,latdim_,loop_
      character*40 bidon_name_
      real , dimension(:) , allocatable :: lonr4_, latr4_
#ifdef synopsis
       subroutinetitle='tide_netcdf_zoom'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Si ktide>1 ne pas recalculer le zoom mais recharger le zoom archive:
      if(ktide>1) then !ne pas recalculer les coordd u zoom pour ktide>1' !11-11-14
        if(ichoix==1) then !11-11-14
         tidezoom_iend=tidezoom_iend_t
         tidezoom_jend=tidezoom_jend_t
         tidezoom_istr=tidezoom_istr_t
         tidezoom_jstr=tidezoom_jstr_t
        endif
        if(ichoix==2) then !11-11-14
         tidezoom_iend=tidezoom_iend_u
         tidezoom_jend=tidezoom_jend_u
         tidezoom_istr=tidezoom_istr_u
         tidezoom_jstr=tidezoom_jstr_u
        endif
        if(ichoix==3) then !11-11-14
         tidezoom_iend=tidezoom_iend_v
         tidezoom_jend=tidezoom_jend_v
         tidezoom_istr=tidezoom_istr_v
         tidezoom_jstr=tidezoom_jstr_v
        endif
        tide_imax=tidezoom_iend-tidezoom_istr+1
        tide_jmax=tidezoom_jend-tidezoom_jstr+1

      if(par%rank==0) then
      write(6,*)'ktide ichoix ',ktide,ichoix
      write(6,*)'tidezoom_iend,tidezoom_istr ',tidezoom_iend,tidezoom_istr
      write(6,*)'tidezoom_jend,tidezoom_jstr ',tidezoom_jend,tidezoom_jstr
      write(6,*)'tide_imax,tide_jmax ',tide_imax,tide_jmax
      endif

      return ! sortir
      endif            !ne pas recalculer les coordd u zoom pour ktide>1'

! Si ktide=1: calculer les coordonnees du zoom

! Le fichier FES2012 étant énorme, on évite que tous les proc allouent en meme temps
! trop de memoire vive: (lignes finalement commentees le 11-11-14)
!     do loop_=0,nbdom-1
!     if(loop_==par%rank) then !--------------------------------------->

      call allocate_forcages(1,7,tide_imax,tide_jmax,2) ! arg1=allouer arg2=tide !09-10-09

! lire longitude et latitude:
      if(ichoix==1) then !111111111111111>
                   status=nf_inq_varid(ncid1,'lon_t',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'longitude_t',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'longitude',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lon_ZHL',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lon',var_id)              !13-07-10
      endif              !111111111111111>
      if(ichoix==2) then !222222222222222>
                   status=nf_inq_varid(ncid1,'lon_u',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'longitude_u',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lon_XHL',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lon',var_id)              !13-07-10
      if(status/=0)status=nf_inq_varid(ncid1,'longitude',var_id)
      endif              !222222222222222>
      if(ichoix==3) then !333333333333333>
                   status=nf_inq_varid(ncid1,'lon_v',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'longitude_v',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lon_YHL',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lon',var_id)              !13-07-10
      if(status/=0)status=nf_inq_varid(ncid1,'longitude',var_id)
      endif              !333333333333333>

      if(status/=0)stop 'err 1681 nom longitude'
      status=nf_inq_var(ncid1,var_id,bidon_name_,k0,londim_,tabdim_,i4)

      if(londim_==2) then !--------------->
      if(k0/=nf_double)stop 'err 1683 type longitude'
       varstart(1)=1 ; varcount(1)=tide_imax
       varstart(2)=1 ; varcount(2)=tide_jmax
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:2)     &
                                 ,varcount(1:2)     &
                                 ,multi2d_lon(1:tide_imax,1:tide_jmax))
      endif               !--------------->
      if(londim_==1) then !--------------->
      if(k0/=nf_double.and.k0/=nf_real)stop 'err 1695 type longitude'
       varstart(1)=1 ; varcount(1)=tide_imax
       if(k0==nf_double)                            &
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:1)     &
                                 ,varcount(1:1)     &
                                 ,multi2d_lon(1:tide_imax,1))
       if(k0==nf_real) then !rrrrr>
          allocate(lonr4_(tide_imax))
           status=nf_get_vara_real(ncid1             &
                                  ,var_id            &
                                  ,varstart(1:1)     &
                                  ,varcount(1:1)     &
                                  ,lonr4_(1:tide_imax))
        multi2d_lon(1:tide_imax,1)=lonr4_(1:tide_imax)
        deallocate(lonr4_)
       endif                !rrrrr>
       do j=1,tide_jmax
       do i=1,tide_imax
        if(multi2d_lon(i,1)>180.)multi2d_lon(i,1)=multi2d_lon(i,1)-360. !26-03-14
        multi2d_lon(i,j)=multi2d_lon(i,1)
       enddo
       enddo
      endif               !--------------->
      if(londim_<1.or.londim_>2)stop 'erreur dim lon tide_netcdf_zoom'

      if(status/=0)stop 'erreur lecture longitude 1'

      if(ichoix==1) then !111111111111111>
                   status=nf_inq_varid(ncid1,'lat_t',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude_t',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lat_ZHL',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lat',var_id)
      endif              !111111111111111>
      if(ichoix==2) then !222222222222222>
                   status=nf_inq_varid(ncid1,'lat_u',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude_u',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lat_XHL',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lat',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude',var_id)
      endif              !222222222222222>
      if(ichoix==3) then !333333333333333>
                   status=nf_inq_varid(ncid1,'lat_v',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude_v',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lat_YHL',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'lat',var_id)
      if(status/=0)status=nf_inq_varid(ncid1,'latitude',var_id)
      endif              !333333333333333>

      if(status/=0)stop 'erreur nom latitude'

      status=nf_inq_var(ncid1,var_id,bidon_name_,k0,latdim_,tabdim_,i4)

      if(latdim_==2) then !--------------->
      if(k0/=nf_double)stop 'err 1748 type latitude'
       varstart(1)=1 ; varcount(1)=tide_imax
       varstart(2)=1 ; varcount(2)=tide_jmax
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:2)     &
                                 ,varcount(1:2)     &
                                 ,multi2d_lat(1:tide_imax,1:tide_jmax))
      endif               !--------------->
      if(latdim_==1) then !--------------->
      if(k0/=nf_double.and.k0/=nf_real)stop 'err 1758 type latitude'
       varstart(1)=1 ; varcount(1)=tide_jmax
       if(k0==nf_double)                            &
       status=nf_get_vara_double( ncid1             &
                                 ,var_id            &
                                 ,varstart(1:1)     &
                                 ,varcount(1:1)     &
                                 ,multi2d_lat(1,1:tide_jmax))
       if(k0==nf_real) then !rrrrr>
          allocate(latr4_(tide_jmax))
           status=nf_get_vara_real(ncid1             &
                                  ,var_id            &
                                  ,varstart(1:1)     &
                                  ,varcount(1:1)     &
                                  ,latr4_(1:tide_jmax))
        multi2d_lat(1,1:tide_jmax)=latr4_(1:tide_jmax)
        deallocate(latr4_)
       endif                !rrrrr>
       do j=1,tide_jmax
       do i=1,tide_imax
        multi2d_lat(i,j)=multi2d_lat(1,j)
       enddo
       enddo
      endif               !--------------->
      if(latdim_<1.or.latdim_>2)stop 'erreur dim lat tide_netcdf_zoom'

! A partir des min & max des (lon,lat) des 2 grilles, reduire la zone d'extraction:
       tidezoom_istr=999999 ; tidezoom_iend=-999999 ! first guess
       tidezoom_jstr=999999 ; tidezoom_jend=-999999 ! first guess

       ksecu=0
!      do j=1,tide_jmax-1
!      do i=1,tide_imax-1
       do j=1,tide_jmax
       do i=1,tide_imax

!      if(par%rank==0)write(10,*)multi2d_lon(i,j),multi2d_lat(i,j)

       if(multi2d_lon(i,j)>=lonmin.and.multi2d_lon(i,j)<=lonmax.and.   &
          multi2d_lat(i,j)>=latmin.and.multi2d_lat(i,j)<=latmax)then

         tidezoom_istr=min(tidezoom_istr,i)
         tidezoom_jstr=min(tidezoom_jstr,j)
         tidezoom_iend=max(tidezoom_iend,i)
         tidezoom_jend=max(tidezoom_jend,j)
         ksecu=1

       endif

       enddo
       enddo

! Si ksecu=0 c'est que le proc est si petit, qu'aucun point de maree ne s'est trouvé à l'interieur
! On refait le test differement. On cherche le point le plus proche.
       if(ksecu==0)then !000000000000000>
        dist1=1.d20 ; i1=imax/2 ; j1=jmax/2
        do j=1,tide_jmax
        do i=1,tide_imax

         dist2=rayonterre*                                         &
         acos( sin(multi2d_lat(i,j)*deg2rad)*sin(lat_t(i1,j1))     &
              +cos(multi2d_lat(i,j)*deg2rad)*cos(lat_t(i1,j1))     &
              *cos(lon_t(i1,j1)-multi2d_lon(i,j)*deg2rad))

         if(dist2<dist1)then !--->
          dist1=dist2 ; i2=i ; j2=j
         endif               !--->

        enddo
        enddo
        tidezoom_istr=i2 ; tidezoom_iend=i2
        tidezoom_jstr=j2 ; tidezoom_jend=j2
       endif            !000000000000000>

       if(tidezoom_istr<1)stop 'Err 1771 tidezoom_istr<1'
       if(tidezoom_jstr<1)stop 'Err 1771 tidezoom_jstr<1'
       if(tidezoom_iend>tidefull_imax) &
       stop 'Err 1771 tidezoom_iend>tidefull_imax'
       if(tidezoom_jend>tidefull_jmax) &
       stop 'Err 1771 tidezoom_jend>tidefull_jmax'


! on elargit un peu plus pour boucher les trous sans être restreint par la taille
! reduite de la zone d'extraction:
       i0=10 ! elargissement à 10 lignes / 10 colonnes supplementaires - repere 1233
!      tidezoom_istr=max0(tidezoom_istr-i0,1)
!      tidezoom_iend=min0(tidezoom_iend+i0,tidefull_imax)
       tidezoom_istr=tidezoom_istr-i0
       tidezoom_iend=tidezoom_iend+i0
! Cet algo pour la direction i vise A prendre en compte une grille
! symphonie depassant de par et d'autre de longitude=0 et une grille FES avec 0<lon<360
       if(tidezoom_istr<1.or. &
          tidezoom_iend>tidefull_imax) then !ooo>
           tidezoom_istr=1 ; tidezoom_iend=tidefull_imax
       endif                                              !ooo>
       tidezoom_jstr=max0(tidezoom_jstr-i0,1)
       tidezoom_jend=min0(tidezoom_jend+i0,tidefull_jmax)

! En attendant que Cyril nous nous une nouvelle fonction d'echange (pour u et v)
! on force le code à passer par les lignes suivantes:
! 15-10-11
#ifdef bidouille
! Pour verifier la conservation de la parallelisation on ne peut pas passer
! par un redecoupage du domaine:
       tidezoom_istr=1 ; tidezoom_iend=tidefull_imax
       tidezoom_jstr=1 ; tidezoom_jend=tidefull_jmax
#endif

       tide_imax=tidezoom_iend-tidezoom_istr+1
       tide_jmax=tidezoom_jend-tidezoom_jstr+1


!      write(*,*)'ichoix=',ichoix
!      write(*,*)'tidezoom_istr=',tidezoom_istr
!      write(*,*)'tidezoom_jstr=',tidezoom_jstr
!      write(*,*)'tidezoom_iend=',tidezoom_iend
!      write(*,*)'tidezoom_jend=',tidezoom_jend
!      write(*,*)'tide_imax=',tide_imax
!      write(*,*)'tide_jmax=',tide_jmax

      call allocate_forcages(2,7,0,0,0)

!     write(6,*)'tidezoom_istr,tidezoom_iend',par%rank,tidezoom_istr,tidezoom_iend
!     write(6,*)'tidezoom_jstr,tidezoom_jend',par%rank,tidezoom_jstr,tidezoom_jend
!     write(6,*)'lonmin,lonmax',par%rank,lonmin,lonmax
!     write(6,*)'latmin,latmax',par%rank,latmin,latmax

!      endif                    !--------------------------------------->
!#ifdef parallele
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif
!      enddo ! fin de boucle pour loop_

! Memoriser le zoom pour le pas le recalculer pour les autres
! harmoniques
      if(ichoix==1) then !11-11-14
       tidezoom_iend_t=tidezoom_iend
       tidezoom_jend_t=tidezoom_jend
       tidezoom_istr_t=tidezoom_istr
       tidezoom_jstr_t=tidezoom_jstr
      endif
      if(ichoix==2) then !11-11-14
       tidezoom_iend_u=tidezoom_iend
       tidezoom_jend_u=tidezoom_jend
       tidezoom_istr_u=tidezoom_istr
       tidezoom_jstr_u=tidezoom_jstr
      endif
      if(ichoix==3) then !11-11-14
       tidezoom_iend_v=tidezoom_iend
       tidezoom_jend_v=tidezoom_jend
       tidezoom_istr_v=tidezoom_istr
       tidezoom_jstr_v=tidezoom_jstr
      endif


      end subroutine tide_netcdf_zoom

!__________________________________________________________________________________

      subroutine tide_previous_run_analysis(txt2_)
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
      include 'netcdf.inc'
      integer istr_,iend_,jstr_,jend_
      character txt2_*2 , txt_*1
      txt_(1:1)=txt2_(1:1)
#ifdef synopsis
       subroutinetitle='tide_previous_run_analysis'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(txt2_(2:2)=='0')texte250=nametide(ktide,1)
      if(txt2_(2:2)=='1')texte250=nametide(ktide,3)

! lire l'identifiant du fichier:
      if(par%rank==0) then !>>>>>
       write(6,*)'dans tide_previous_run_analysis je lis:'
       write(6,'(a)')trim(texte250)
       write(6,'(a,a)')'pour la variable ',txt_
      endif                !>>>>>

      status=nf_open(trim(texte250),nf_nowrite,ncid1)
      if(status/=0)then
         write(6,'(3a)')'File ',trim(texte250),' not found'
         stop ' 1930 fichier maree absent'
         status=nf_close(ncid1)
         return
      endif

! lire l'identifiant des dimensions:
      if(txt_   =='l'.or. &
         txt_   =='h') then ! h h h h h h >
                   status=nf_inq_dimid(ncid1,'x_ZHL',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'imax_t',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'ni_t',dim_x_id)
      endif                 ! h h h h h h >
      if(txt_   =='u') then ! u u u u u u >
                   status=nf_inq_dimid(ncid1,'x_XHL',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'imax_u',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'ni_u',dim_x_id)
      endif                 ! u u u u u u >
      if(txt_   =='v') then ! v v v v v v >
                   status=nf_inq_dimid(ncid1,'x_YHL',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'imax_v',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'ni_v',dim_x_id)
      endif                 ! v v v v v v >
      if(status/=0)stop 'erreur nom pour dim_x_id'

      if(txt_   =='l'.or. &
         txt_   =='h') then ! h h h h h h >
                   status=nf_inq_dimid(ncid1,'y_ZHL',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'jmax_t',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'nj_t',dim_y_id)
      endif                 ! h h h h h h >
      if(txt_   =='u') then ! u u u u u u >
                   status=nf_inq_dimid(ncid1,'y_XHL',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'jmax_u',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'nj_u',dim_y_id)
      endif                 ! u u u u u u >
      if(txt_   =='v') then ! v v v v v v >
                   status=nf_inq_dimid(ncid1,'y_YHL',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'jmax_v',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'nj_v',dim_y_id)
      endif                 ! v v v v v v >
      if(status/=0)stop 'erreur nom pour dim_y_id'

! Lecture des valeurs des dimensions
      status=nf_inq_dimlen(ncid1,dim_x_id,tide_imax)
      status=nf_inq_dimlen(ncid1,dim_y_id,tide_jmax)

      i0=0 ; j0=0
      if(txt_=='u')i0=1
      if(txt_=='v')j0=1
      if(tide_imax/=iglb+2-i0)stop 'Erreur dimension Oi fichier analyse'
      if(tide_jmax/=jglb+2-j0)stop 'Erreur dimension Oj fichier analyse'
      tide_imax=imax+2-i0 ; tide_jmax=jmax+2-j0

      if(txt_=='l'.or. &
         txt_=='h') then !t-t-t-t-t->
       istr_=0 ; iend_=imax+1 ; jstr_=0 ; jend_=jmax+1
      endif              !t-t-t-t-t->
      if(txt_=='u') then !u-u-u-u-u->
       istr_=1 ; iend_=imax+1 ; jstr_=0 ; jend_=jmax+1
      endif              !u-u-u-u-u->
      if(txt_=='v') then !v-v-v-v-v->
       istr_=0 ; iend_=imax+1 ; jstr_=1 ; jend_=jmax+1
      endif              !v-v-v-v-v->

      allocate(multi2d_var(istr_:iend_,jstr_:jend_,2))

! reset:
      varstart(:)=1 ; varcount(:)=1

! lire amplitude:
      if(txt_   =='h') then ! h h h h h h >
                    status=nf_inq_varid(ncid1,'Ha',var_id)
      endif                 ! h h h h h h >
      if(txt_   =='l') then ! h h h h h h >
                    status=nf_inq_varid(ncid1,'LSAa',var_id)
      endif                 ! h h h h h h >
      if(txt_   =='u') then ! u u u u u u >
                    status=nf_inq_varid(ncid1,'Ua',var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'OiVa',var_id)
      endif                 ! u u u u u u >
      if(txt_   =='v') then ! v v v v v v >
                    status=nf_inq_varid(ncid1,'Va',var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'OjVa',var_id)
      endif                 ! v v v v v v >
      if(status/=0)stop 'erreur nom amplitude analyse'

      status=nf_inq_vartype(ncid1,var_id,var_nftype) !12-03-17
      if(status/=0)stop 'Err 1929 nf_inq_vartype'
      if(var_nftype/=nf_real) then !pmxpmx>
        write(6,*)'Type of tide variable not recognized'
        stop 'Err 1930 var_nftype/=nf_real'
      endif                        !pmxpmx>

      varstart(1)=1+par%timax(1) ; varcount(1)=iend_-istr_+1
      varstart(2)=1+par%tjmax(1) ; varcount(2)=jend_-jstr_+1

      status=nf_get_vara_real(ncid1,var_id,varstart       &
                                          ,varcount       &
                                          ,multi2d_var(:,:,1))
      if(status/=0)stop 'erreur lecture amplitude analyse'

      status=nf_get_att_real(ncid1,var_id,'_FillValue',filval) !11-06-17
      if(status/=0)stop 'err filval amplitude analyse'

! lire retard de phase:
      if(txt_   =='h') then ! h h h h h h >
                    status=nf_inq_varid(ncid1,'Hg',var_id)
      endif                 ! h h h h h h >
      if(txt_   =='l') then ! h h h h h h >
                    status=nf_inq_varid(ncid1,'LSAg',var_id)
      endif                 ! h h h h h h >
      if(txt_   =='u') then ! u u u u u u >
                    status=nf_inq_varid(ncid1,'Ug',var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'OiVg',var_id)
      endif                 ! u u u u u u >
      if(txt_   =='v') then ! v v v v v v >
                    status=nf_inq_varid(ncid1,'Vg',var_id)
       if(status/=0)status=nf_inq_varid(ncid1,'OjVg',var_id)
      endif                 ! v v v v v v >
      if(status/=0)stop 'erreur nom retard de phase analyse'

      status=nf_inq_vartype(ncid1,var_id,var_nftype) !12-03-17
      if(status/=0)stop 'Err 1962 nf_inq_vartype'
      if(var_nftype/=nf_real) then !pmxpmx>
        write(6,*)'Type of tide variable not recognized'
        stop 'Err 1965 var_nftype/=nf_real'
      endif                        !pmxpmx>

      varstart(1)=1+par%timax(1) ; varcount(1)=iend_-istr_+1
      varstart(2)=1+par%tjmax(1) ; varcount(2)=jend_-jstr_+1
      status=nf_get_vara_real(ncid1,var_id,varstart    &
                                          ,varcount    &
                                          ,multi2d_var(:,:,2))
      if(status/=0)stop 'erreur lecture retard phase analyse'

! Transformer les amplitude/retard_de_phase en amplitudes complexes:
      flag_stop=0

! CAS TABLEAUX POUR FORCER AUX FRONTIERES:
      if(txt2_(2:2)=='0')then !s-s-s-s-s-s-s-s-s-s-s-s-s-s-s-s-s-s->

      if(txt_   =='h') then ! h h h h h h >
      do j=0,jmax+1
      do i=0,imax+1
       sshtidecos_w(i,j,ktide)=multi2d_var(i,j,1)*cos(multi2d_var(i,j,2)*deg2rad)*mask_t(i,j,kmax) !21-05-17
       sshtidesin_w(i,j,ktide)=multi2d_var(i,j,1)*sin(multi2d_var(i,j,2)*deg2rad)*mask_t(i,j,kmax)
       if(mask_t(i,j,kmax)==1) then !ooo>
        if(multi2d_var(i,j,1)==filval.or. &
           multi2d_var(i,j,2)==filval) then
            flag_stop=1
            write(10+par%rank,*)'Err 2013 Tide filval found on mask_t=1' !12-03-17
            write(10+par%rank,*)'ktide=',ktide
            write(10+par%rank,*)'par%rank',par%rank
            write(10+par%rank,*)'i,j,loc,imax,jmax',i,j,imax,jmax
            write(10+par%rank,*)'i,j,iglb,jglb',i+par%timax(1),j+par%tjmax(1),iglb,jglb
            write(10+par%rank,*)'multi2d_var(i,j,1:2)',multi2d_var(i,j,1:2)
            write(10+par%rank,*)'filval',filval
            write(10+par%rank,*)'ssh cos sin',sshtidecos_w(i,j,ktide),sshtidesin_w(i,j,ktide)
            write(10+par%rank,*)'Err 2013 Tide filval found on mask_h=1' !12-03-17
        endif
       endif                        !ooo>
      enddo
      enddo
      endif                 ! h h h h h h >
      if(txt_   =='l') then ! l l l l l l > !16-11-13
      do j=0,jmax+1
      do i=0,imax+1
      potidecos_w(i,j,ktide)=multi2d_var(i,j,1)*cos(multi2d_var(i,j,2)*deg2rad)*mask_t(i,j,kmax) !21-05-17
      potidesin_w(i,j,ktide)=multi2d_var(i,j,1)*sin(multi2d_var(i,j,2)*deg2rad)*mask_t(i,j,kmax) !21-05-17
       if(mask_t(i,j,kmax)==1) then !ooo>
        if(multi2d_var(i,j,1)==filval.or. &
           multi2d_var(i,j,2)==filval) then
           write(10+par%rank,*) 'Err 2027 Tide filval found on mask_t=1' !12-03-17
           flag_stop=1
        endif
       endif                        !ooo>
      enddo
      enddo
      endif                 ! l l l l l l >
      if(txt_   =='u') then ! u u u u u u >
      do j=0,jmax+1
      do i=1,imax+1
      veltidecos_u(i,j,ktide)=multi2d_var(i,j,1)*cos(multi2d_var(i,j,2)*deg2rad)*mask_u(i,j,kmax) !21-05-17
      veltidesin_u(i,j,ktide)=multi2d_var(i,j,1)*sin(multi2d_var(i,j,2)*deg2rad)*mask_u(i,j,kmax) !21-05-17
       if(mask_u(i,j,kmax)==1) then !ooo>
        if(multi2d_var(i,j,1)==filval.or. &
           multi2d_var(i,j,2)==filval) then !(°L°)>
           !if(j+par%tjmax(1)/=0.and.j+par%tjmax(1)/=jglb)then !pmx>
           if(j+par%tjmax(1)/=0.and.j+par%tjmax(1)/=jglb+1)then !thom> !28-06-17
            flag_stop=1
            write(10+par%rank,*)'Err 2040 Tide filval found on mask_u=1'
            write(10+par%rank,*)'par%rank',par%rank
            write(10+par%rank,*)'i,j,loc,imax,jmax',i,j,imax,jmax
            write(10+par%rank,*)'i,j,iglb,jglb',i+par%timax(1),j+par%tjmax(1),iglb,jglb
            write(10+par%rank,*)'multi2d_var(i,j,1:2)',multi2d_var(i,j,1:2)
            write(10+par%rank,*)'filval',filval
            write(10+par%rank,*)'Err 2040 Tide filval found on mask_u=1' !12-03-17
           endif                                              !pmx>
        endif                               !(°L°)>
       endif                        !ooo>
      enddo
      enddo
      endif                 ! u u u u u u >
      if(txt_   =='v') then ! v v v v v v >
      do j=1,jmax+1
      do i=0,imax+1
      veltidecos_v(i,j,ktide)=multi2d_var(i,j,1)*cos(multi2d_var(i,j,2)*deg2rad)*mask_v(i,j,kmax) !21-05-17
      veltidesin_v(i,j,ktide)=multi2d_var(i,j,1)*sin(multi2d_var(i,j,2)*deg2rad)*mask_v(i,j,kmax) !21-05-17
       if(mask_v(i,j,kmax)==1) then !ooo>
        if(multi2d_var(i,j,1)==filval.or. &
           multi2d_var(i,j,2)==filval) then !(O_O)>
           !if(i+par%timax(1)/=0.and.i+par%timax(1)/=iglb)then !pmx>
           if(i+par%timax(1)/=0.and.i+par%timax(1)/=iglb+1)then !thom> !28-06-17
            flag_stop=1
            write(10+par%rank,*)'Err 2053 Tide filval found on mask_v=1'
            write(10+par%rank,*)'par%rank',par%rank
            write(10+par%rank,*)'i,j,loc,imax,jmax',i,j,imax,jmax
            write(10+par%rank,*)'i,j,iglb,jglb',i+par%timax(1),j+par%tjmax(1),iglb,jglb
            write(10+par%rank,*)'multi2d_var(i,j,1:2)',multi2d_var(i,j,1:2)
            write(10+par%rank,*)'filval',filval
            write(10+par%rank,*)'Err 2053 Tide filval found on mask_v=1' !12-03-17
           endif                                              !pmx>
        endif                               !(O_O)>
       endif                        !ooo>
      enddo
      enddo
      endif                 ! v v v v v v >

      endif                   !s-s-s-s-s-s-s-s-s-s-s-s-s-s-s-s-s-s->

! CAS TABLEAUX POUR DETIDER LE SIGNAL
      if(txt2_(2:2)=='1')then !o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o->

      if(txt_   =='h') then ! h h h h h h >
      do j=0,jmax+1
      do i=0,imax+1
      sshtidecosout_w(i,j,ktide)=multi2d_var(i,j,1)*cos(multi2d_var(i,j,2)*deg2rad)*mask_t(i,j,kmax) !21-05-17
      sshtidesinout_w(i,j,ktide)=multi2d_var(i,j,1)*sin(multi2d_var(i,j,2)*deg2rad)*mask_t(i,j,kmax) !21-05-17
       if(mask_t(i,j,kmax)==1) then !ooo>
        if(multi2d_var(i,j,1)==filval.or. &
           multi2d_var(i,j,2)==filval) then
           write(10+par%rank,*)'Err 2069 Tide filval found on mask_t=1' !12-03-17
           flag_stop=1
        endif
       endif                        !ooo>
      enddo
      enddo
      endif                 ! h h h h h h >
      if(txt_   =='u') then ! u u u u u u >
      do j=0,jmax+1
      do i=1,imax+1
      veltidecosout_u(i,j,ktide)=multi2d_var(i,j,1)*cos(multi2d_var(i,j,2)*deg2rad)*mask_u(i,j,kmax) !21-05-17
      veltidesinout_u(i,j,ktide)=multi2d_var(i,j,1)*sin(multi2d_var(i,j,2)*deg2rad)*mask_u(i,j,kmax) !21-05-17
       if(mask_u(i,j,kmax)==1) then !ooo>
        if(multi2d_var(i,j,1)==filval.or. &
           multi2d_var(i,j,2)==filval) then
          !if(j+par%tjmax(1)/=0.and.j+par%tjmax(1)/=jglb)then
          if(j+par%tjmax(1)/=0.and.j+par%tjmax(1)/=jglb+1)then !thom> !28-06-17
           write(10+par%rank,*) 'Err 2082 Tide filval found on mask_u=1' !12-03-17
           flag_stop=1
          endif
        endif
       endif                        !ooo>
      enddo
      enddo
      endif                 ! u u u u u u >
      if(txt_   =='v') then ! v v v v v v >
      do j=1,jmax+1
      do i=0,imax+1
      veltidecosout_v(i,j,ktide)=multi2d_var(i,j,1)*cos(multi2d_var(i,j,2)*deg2rad)*mask_v(i,j,kmax) !21-05-17
      veltidesinout_v(i,j,ktide)=multi2d_var(i,j,1)*sin(multi2d_var(i,j,2)*deg2rad)*mask_v(i,j,kmax) !21-05-17
       if(mask_v(i,j,kmax)==1) then !ooo>
        if(multi2d_var(i,j,1)==filval.or. &
           multi2d_var(i,j,2)==filval) then
          !if(i+par%timax(1)/=0.and.i+par%timax(1)/=iglb)then !pmx>
          if(i+par%timax(1)/=0.and.i+par%timax(1)/=iglb+1)then !thom> !28-06-17
           write(10+par%rank,*) 'Err 2097 Tide filval found on mask_v=1' !12-03-17
           flag_stop=1
          endif
        endif
       endif                        !ooo>
      enddo
      enddo
      endif                 ! v v v v v v >

      endif                   !o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o->

!     call allocate_forcages(2,7,0,0,0) ! arg1=desallouer arg2=tide
      deallocate(multi2d_var)

      status=nf_close(ncid1)

      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0) then
       write(6,*)'par%rank flag_stop k0',par%rank,flag_stop,k0
       stop 'Err tide_previous_run_analysis see fort file' !11-06-17
      endif

      end subroutine tide_previous_run_analysis

!__________________________________________________________________________________

      subroutine initial_tide_mpi_obc
#ifdef parallele
      use module_principal
      use module_parallele
      implicit none
      integer loop_ ,potidecos_w_z0_id_ ,potidesin_w_z0_id_   &
                   ,sshtidecos_w_z0_id_,sshtidesin_w_z0_id_   &
                   ,veltidecos_u_u1_id_,veltidesin_u_u1_id_   &
                   ,veltidecos_v_v1_id_,veltidesin_v_v1_id_
#ifdef synopsis
       subroutinetitle='initial_tide_mpi_obc'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       call get_type_echange('z0','potidecos_w_z0'    &       !24-07-14
                                  ,potidecos_w        &
                           ,lbound(potidecos_w)       &
                           ,ubound(potidecos_w)       &
                                  ,potidecos_w_z0_id_)

       call get_type_echange('z0','potidesin_w_z0'    &
                                  ,potidesin_w        &
                           ,lbound(potidesin_w)       &
                           ,ubound(potidesin_w)       &
                                  ,potidesin_w_z0_id_)

       call get_type_echange('z0','sshtidecos_w_z0'    &
                                  ,sshtidecos_w        &
                           ,lbound(sshtidecos_w)       &
                           ,ubound(sshtidecos_w)       &
                                  ,sshtidecos_w_z0_id_)

       call get_type_echange('z0','sshtidesin_w_z0'    &
                                  ,sshtidesin_w        &
                           ,lbound(sshtidesin_w)       &
                           ,ubound(sshtidesin_w)       &
                                  ,sshtidesin_w_z0_id_)

       call get_type_echange('u1','veltidecos_u_u1'    &
                                  ,veltidecos_u        &
                           ,lbound(veltidecos_u)       &
                           ,ubound(veltidecos_u)       &
                                  ,veltidecos_u_u1_id_)

       call get_type_echange('u1','veltidesin_u_u1'    &
                                  ,veltidesin_u        &
                           ,lbound(veltidesin_u)       &
                           ,ubound(veltidesin_u)       &
                                  ,veltidesin_u_u1_id_)

       call get_type_echange('v1','veltidecos_v_v1'    &
                                  ,veltidecos_v        &
                           ,lbound(veltidecos_v)       &
                           ,ubound(veltidecos_v)       &
                                  ,veltidecos_v_v1_id_)

       call get_type_echange('v1','veltidesin_v_v1'    &
                                  ,veltidesin_v        &
                           ,lbound(veltidesin_v)       &
                           ,ubound(veltidesin_v)       &
                                  ,veltidesin_v_v1_id_)

      ! Echanges
      do loop_=1,subcycle_exchange

        call echange_voisin(potidecos_w,potidecos_w_z0_id_,mpi_neighbor_list(loop_))
        call echange_voisin(potidesin_w,potidesin_w_z0_id_,mpi_neighbor_list(loop_))

        call echange_voisin(sshtidecos_w,sshtidecos_w_z0_id_,mpi_neighbor_list(loop_))
        call echange_voisin(sshtidesin_w,sshtidesin_w_z0_id_,mpi_neighbor_list(loop_))

        call echange_voisin(veltidecos_u,veltidecos_u_u1_id_,mpi_neighbor_list(loop_))
        call echange_voisin(veltidesin_u,veltidesin_u_u1_id_,mpi_neighbor_list(loop_))

        call echange_voisin(veltidecos_v,veltidecos_v_v1_id_,mpi_neighbor_list(loop_))
        call echange_voisin(veltidesin_v,veltidesin_v_v1_id_,mpi_neighbor_list(loop_))

      enddo
      call loc_wait() ! ----> important sinon pas d'echanges
#endif

      end subroutine initial_tide_mpi_obc

!.................................................................
      subroutine initial_tide_add_phase_offset
      use module_principal
      use module_parallele

! Add phase off set
!     if(ktide==1)multi2d_var(:,:,2)=multi2d_var(:,:,2)-8.
!     if(ktide==2)multi2d_var(:,:,2)=multi2d_var(:,:,2)-8.
!     if(ktide==3)multi2d_var(:,:,2)=multi2d_var(:,:,2)-8.
!     if(ktide==4)multi2d_var(:,:,2)=multi2d_var(:,:,2)-8.
!     if(ktide==5)multi2d_var(:,:,2)=multi2d_var(:,:,2)-8.
!     if(ktide==6)multi2d_var(:,:,2)=multi2d_var(:,:,2)-8.
!     if(ktide==7)multi2d_var(:,:,2)=multi2d_var(:,:,2)-8.
!     if(ktide==8)multi2d_var(:,:,2)=multi2d_var(:,:,2)-8.
!     if(ktide==9)multi2d_var(:,:,2)=multi2d_var(:,:,2)-8.

      end subroutine initial_tide_add_phase_offset

!.................................................................

      subroutine tide_ssh_bias_correction(txt_) !17-04-18
      use module_principal
      use module_parallele !#MPI
      use module_forcages
      implicit none
      include 'netcdf.inc'
      integer istr_,iend_,jstr_,jend_
      character*3 txt_

      if(index(nametide(ktide,6),'none')/=0)return
      texte250=nametide(ktide,6)

! lire l'identifiant du fichier:
      if(par%rank==0) then !>>>>>
       write(6,*)'Ouverture de:'
       write(6,'(a)')trim(texte250)
      endif                !>>>>>

      status=nf_open(trim(texte250),nf_nowrite,ncid1)
      if(status/=0)then
         write(6,'(3a)')'File ',trim(texte250),' not found'
         stop ' 1964 fichier maree absent'
         status=nf_close(ncid1)
         return
      endif

! lire l'identifiant des dimensions:
                   status=nf_inq_dimid(ncid1,'x_ZHL',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'imax_t',dim_x_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'ni_t',dim_x_id)
      if(status/=0)stop 'erreur nom pour dim_x_id'

                   status=nf_inq_dimid(ncid1,'y_ZHL',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'jmax_t',dim_y_id)
      if(status/=0)status=nf_inq_dimid(ncid1,'nj_t',dim_y_id)
      if(status/=0)stop 'erreur nom pour dim_y_id'

! Lecture des valeurs des dimensions
      status=nf_inq_dimlen(ncid1,dim_x_id,i1)
      status=nf_inq_dimlen(ncid1,dim_y_id,j1)

      if(i1/=iglb+2)stop 'Erreur dimension Oi fichier ssh correction'
      if(j1/=jglb+2)stop 'Erreur dimension Oj fichier ssh correction'

      istr_=0 ; iend_=imax+1 ; jstr_=0 ; jend_=jmax+1
      allocate(multi2d_var(istr_:iend_,jstr_:jend_,2))

! reset:
      varstart(:)=1 ; varcount(:)=1

! lire correction amplitude cosinus:
         status=nf_inq_varid(ncid1,'Hcos-ObcHcos',var_id)
      if(status/=0)stop 'err nf_inq_varid Hcos-ObcHcos'

      status=nf_inq_vartype(ncid1,var_id,var_nftype) !12-03-17
      if(status/=0)stop 'Err nf_inq_vartype Hcos-ObcHcos'
      if(var_nftype/=nf_real) then !pmxpmx>
        write(6,*)'Type Hcos-ObcHcos not recognized'
        stop 'Err Hcos-ObcHcos var_nftype/=nf_real'
      endif                        !pmxpmx>

      varstart(1)=1+par%timax(1) ; varcount(1)=iend_-istr_+1
      varstart(2)=1+par%tjmax(1) ; varcount(2)=jend_-jstr_+1

      status=nf_get_vara_real(ncid1,var_id,varstart       &
                                          ,varcount       &
                                          ,multi2d_var(:,:,1))
      if(status/=0)stop 'Err nf_get_vara_real Hcos-ObcHcos'

      status=nf_get_att_real(ncid1,var_id,'_FillValue',filval) !11-06-17
      if(status/=0)stop 'err filval amplitude analyse'

! lire correction amplitude sinus
         status=nf_inq_varid(ncid1,'Hsin-ObcHsin',var_id)
      if(status/=0)stop 'err nf_inq_varid Hsin-ObcHsin'

      status=nf_inq_vartype(ncid1,var_id,var_nftype) !12-03-17
      if(status/=0)stop 'Err Hsin-ObcHsin nf_inq_vartype'
      if(var_nftype/=nf_real) then !pmxpmx>
        write(6,*)'Type Hsin-ObcHsin not recognized'
        stop 'Err Hsin-ObcHsin var_nftype/=nf_real'
      endif                        !pmxpmx>

      varstart(1)=1+par%timax(1) ; varcount(1)=iend_-istr_+1
      varstart(2)=1+par%tjmax(1) ; varcount(2)=jend_-jstr_+1
      status=nf_get_vara_real(ncid1,var_id,varstart    &
                                          ,varcount    &
                                          ,multi2d_var(:,:,2))
      if(status/=0)stop 'Err nf_get_vara_real Hsin-ObcHsin'
      status=nf_close(ncid1)

! CORRIGER SSH AVEC LES TERMES DE CORRECTION:
       k0=0
       if(txt_=='add')k0=1  ! Ajoute la correction en debut de simulation
       if(txt_=='rmv')k0=-1 ! Enleve en fin de simu pour ecrire le "vrai" forcage fes dans le fichier final
       if(k0==0)stop 'Err unrecognized arg in tide_ssh_bias_correction'
       do j=0,jmax+1 ; do i=0,imax+1
        if(multi2d_var(i,j,1)==filval)multi2d_var(i,j,1)=0. !03-11-21
        if(multi2d_var(i,j,2)==filval)multi2d_var(i,j,2)=0. !03-11-21
        sshtidecos_w(i,j,ktide)=sshtidecos_w(i,j,ktide)-k0*multi2d_var(i,j,1)
        sshtidesin_w(i,j,ktide)=sshtidesin_w(i,j,ktide)-k0*multi2d_var(i,j,2)
       enddo         ; enddo


      deallocate(multi2d_var)
      end subroutine tide_ssh_bias_correction
