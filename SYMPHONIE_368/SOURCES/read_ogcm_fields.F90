      subroutine read_ogcm_fields(case_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 358 - last update: 28-10-22
!______________________________________________________________________
!    _________                    .__                  .__             !(°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !m°v°m 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!______________________________________________________________________
      use module_principal
      use module_parallele !#MPI
      use module_ogcm
      implicit none
      integer case_,tm_,var_count_
#ifdef synopsis
       subroutinetitle='read_ogcm_fields'
       subroutinedescription=' Driver for ogcm forcing'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date      Description des modifications
!         16/06/01: le tableau OBCFILE devient bidimensionnel.
!         27/06/01: KI1 n'est plus passé en argument.
!         27/06/01: bug sur taille de LREC
!         25/10/01: les variables 3D vont maintenant de 1 à NR-1
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         14/03/02: amenagement pour model_ inverse. On lit le vrai courant
!                   et le courant moyen et non plus le courant tilde et le
!                   transport, une conversion est donc necessaire
!                   + debugage de l'initialisation (initial_with_obc initialisait
!                   eventuellement le model_ avec NC-1 artificiel entrainant
!                   une incompatibilité avec les obc
!         26/08/02: VBAROBC remplacé par VMEAOBC
!                   & VHZOBC remplacé par VELOBC
!         22/01/04: on precise que les fichiers sont FORM='UNFORMATTED'
!         28/03/04: redemarrage du model_ avec un fichier particulier (case_=0)
!         15/04/04: suite du point precedent. Attention a ce que l'echeance
!                   suivante ne soit pas celle du fichier de redemarrage.
!         16/04/04: multiplier par le masque pour éviter les "flags" de VIFOP
!         06/05/04: option pour commencer avant la 1ere echeance disponible
!         21/06/04: ajout d'une aide pour bien définir notebook_time dans le cas
!                   du fichier redemarrage d'read_ogcm_fields
!                   + ajout d'une contrainte de moyenne de sse = 0
!         28/01/05: ajustement convectif des champs de forcage
!         08/03/05: debugage suite au point precedent. Les variables K2 et
!                   KOUNTMOD sont abandonnées.
!         12/04/05: ajustement convectif: en + on restratifie tres legerement
!                   pour eviter d'entretenir de l'instabilité convective
!                   aux conditions aux limites
!         27/07/05: l'option moyenne sse=0 finalement commentée
!         30/09/05: modif pour r8 possible (ANYV3D remplace ANYVAR3D)
!         20/01/06: conservation du toit: creation de ztaobc_cum pour
!                   suivi de l'ogcm
!         21/04/06: Melange de l'ogcm par la marée dans la couche de fond
!         23/10/06: Un "close(3)" qui manquait, mais qui n'était pas un bug...
!         11/12/06: option case_.eq.-1 pour faire un fichier restart
!         13/04/07: possibilité d'initialiser le courant avec le courant
!                   geostrophique
!         17/04/07: Passage à coordonnees curvilignes (voir ajout dx_y et cie...)
!         21/04/07: Passage à coordonnees curvilignes (voir ajout dx_y et cie...)
!         30/11/07: affichages à l'ecran
!         23/01/08: renouvellement des echeances par modulo d'entier supprimé
!         08/12/08: ajout d'une filiere d'interpolation online
!         29/12/08: securite relative a IOBC_LR placée dans case_=1
!         12-05-09: Par defaut, on prend la ssh moyenne de l'ogcm
!         29-05-09: Parallelisation
!         31-05-09: Parallelisation: traitement de la moyenne spatiale de la ssh
! 2009.3  01-10-09: utilisation des nouveaux facteurs d'echelle verticale
!         02-10-09: dti_lp remplace dti
!         05-10-09: ajout d'un "ifdef parallele"
! 2010.2  27-12-09: obc_af renommé read_ogcm_fields
! 2010.9  06-06-10  tridiagonalsolver renommé tridiagonalsolver
! 2010.11 16-07-10  temobc & salobc renommés temobc & salobc
! 2010.13 03-11-10  des arguments passés dans date_to_kount
! 2010.14 22-11-10  prise en compte equation d'etat non lineaire
!         23-11-10  suppression calcul incompatible avec equation d'etat compressible
! 2010.20 15-04-11  le renouvellement des fichiers est décidé sur la base d'un temps
!                   en secondes et non plus sur la basede kount
! 2010.22 30-04-11  Ajuster velobc pour controler temobc*velobc avec la subroutine
!                   obc_ts_fluxes
! 2010.25 01-04-12  module_ogcm remplace interp_ogcm
! S.26    29-11-13  simplification notebook_obcforcing
!         09-02-14  Ne remettre a jour les champs externes que si subcycle_synchro=1
!         29-04-14  Modif des messages
!         17-05-14  Pour ajustement convectif appeler la meme EOS que le reste
!                   du modele
!         13-06-14  generalisation de la routine obc_ts_fluxes
!         11-08-14  gestion de la singularite a l'equateur dans la parametrisation
!                   de l'epaisseur de la couche turbulente
!         28-10-14  Melange convectif de l'ogcm base sur EOS potentielle a pression
!                   de reference locale
!         26-11-14  mineur: taille de boucle
!         27-10-15  amelioration de l'algo de melange de la maree qui complete l'interpolation 
!                   de mercator
!         29-11-15  ajout du cas 1DV
!         20-04-16  debug obc_ts_fluxes
!         20-03-17  suppression sshobc_cum
!         26-04-17  dans l'evaluation de l'epaisseur de la couche de
!                   melange d la maree il est plus robuste de remplacer cd par une constante
!         10-11-18  tidal_ogcm_mix: un evaluation plus precise du courant de maree total
!                   pour eviter la surestimation de l'ancien algo
! v253    02-05-19  Evitement division par dz
! v259    03-08-19  l'homogeneisation de T,S dans la couche limite de fond de la maree
!                   obtenue par une methode de moyenne sur la verticale simple
! v269    11-12-19  if(flag_ogcmtidemixing==1)call tidal_ogcm_mix
! v309    13-09-21  message pour aide au debugage
!         22-09-21  revenir A la maree maximale dans subroutine tidal_ogcm_mix(tm_)
! v327    19-02-22  if(flag_ogcm_instab==1) call ajust_convect(
!                   if(flag_ogcm_instab==2) call ajust_convect_mix(
!                   if(flag_ogcm_instab==3) call ajust_convect_min(
! v328    23-02-22  appels aux EOS comme dans presgrad
! v358    28-10-22  if(obctime_order==4) 
!...............................................................................


      var_count_=1  ! reset defaut tous modeles


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                        PARTIE I
!                INITIALISATION DU MODELE
!     appelée dans initial_main.F90
! DEBUT:
      if (case_==1) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! securité:
      if(iobc_lr.eq.1)stop 'je ne dois pas passé par read_ogcm_fields.f'         !29/12/08

! Cas interpolation des fichiers ogcm source:
      if(obc_option==2)call ogcm_initial

! note: si kount=k2 NC serait pile egal à 1
! dans un cas general l'etat initial ne coincide pas dans le temps
! avec un etat analysé donné mais se trouve quelque part entre
! deux états d'où la prise en compte de deux reccord: NC et NC+1

!...........................................................................
! ETAPE 1:
! la boucle 2000 permet de charger les 2 echeances qui encadrent l'iteration
! de l'etat initial
! ajout le:                                                            !23/01/08

      do  2000 tm_=0,2,2

       do var_count_=1,var_num    ! Boucle sur variables traceurs, vitesses, ssh
        if(tm_==0)nc=ogcm_rec_prev(var_count_)
        if(tm_==2)nc=ogcm_rec_next(var_count_)
        if(nc.le.0)then !---- debug ----->                             !06/05/04
         write(6,*)'-----------------'
         write(6,'(a,a,a,a)')' No OGCM field available at the'        & !29-04-14
          ,' departure time. The inconsistency between the departure' &
          ,' time in notebook_time and the forcing lists of'          &
          ,' notebook_obcforcing must be corrected.'
         write(6,*)'nc=',nc
         write(6,*)'var_count_=',var_count_
         write(6,*)'tm_=',tm_
         write(6,*)'ogcm_rec_prev(var_count_)=',ogcm_rec_prev(var_count_)
         write(6,*)'ogcm_rec_next(var_count_)=',ogcm_rec_next(var_count_)
!        write(6,*)'ogcm_readtime_prev(var_count_)=',ogcm_readtime_prev(var_count_) !13-09-21
!        call elapsedtime2date(ogcm_readtime_prev(var_count_),i1,i2,i3,i4,i5,i6)
!        write(6,*)'and corresponding date',i1,i2,i3,i4,i5,i6
         write(6,*)'ogcm_readtime_next(var_count_)=',ogcm_readtime_next(var_count_) !13-09-21
         call elapsedtime2date(ogcm_readtime_next(var_count_),i1,i2,i3,i4,i5,i6)
         write(6,*)'and corresponding date',i1,i2,i3,i4,i5,i6

         stop ' STOP read_ogcm_fields choix 1'
        endif           !---- debug ----->
! Cas interpolation des fichiers:
        if(obc_option==2)call ogcm_interp(var_count_,nc,tm_)
! Cas fichiers déjà interpolés sur la grille:
!       if(obc_option==0)call already_interpolated(nc,tm_)           !26-04-11

       enddo ! var_count_ loop end

! Calcul densite pour ajustement convectif:
      if(eos_author==0) then !eqeqeqeqeqeqeqeqeq>  !22-11-10

! Equation d'état linéaire:
      const1=-rho*alp_t
      const2= rho*alp_s
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax+1
      x2=(const1*(temobc_t(i,j,k,tm_)-t0)       & ! densite
         +const2*(salobc_t(i,j,k,tm_)-s0))
      anyv3d(i,j,k,1)=x2
      enddo
      enddo
      enddo

      else                   !eqeqeqeqeqeqeqeqeq>   !22-11-10
! Equation d'état non linéaire:

      if(eos_author==3) then
! charger anyv3d(i,j,k,1) avec la densite potentielle:
        if(eos_tkezref>=0.) then
          call equation_of_state_potzref_jmfwg('obc',tm_) !23-02-22
        else
          call equation_of_state_potloc_jmfwg('obc',tm_) !28-10-14
        endif
      else
       stop 'read_ogcm_fields eos_author/=3 pas encore prevu'
      endif

      endif                  !eqeqeqeqeqeqeqeqeq>  !22-11-10

      if(flag_ogcm_instab==1) call ajust_convect(tm_)     !19-02-22
      if(flag_ogcm_instab==2) call ajust_convect_mix(tm_) !19-02-22
      if(flag_ogcm_instab==3) call ajust_convect_min(tm_) !19-02-22

      if(flag_ogcmtidemixing==1)call tidal_ogcm_mix(tm_) !11-12-19

!..............................................................................
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax+1
       velobc_u(i,j,k,tm_)=velobc_u(i,j,k,tm_)*mask_u(i,j,k)    !16/04/04
       velobc_v(i,j,k,tm_)=velobc_v(i,j,k,tm_)*mask_v(i,j,k)    !16/04/04
      enddo
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax+1
       velbarobc_u(i,j,tm_)=velbarobc_u(i,j,tm_)*mask_u(i,j,kmax+1) !16/04/04
       velbarobc_v(i,j,tm_)=velbarobc_v(i,j,tm_)*mask_v(i,j,kmax+1) !16/04/04
      enddo
      enddo
!..............................................................................

      if(flag_1dv==1) then !1dv1dv1dv> !29-11-15
        call ogcm_1dv_geos(tm_) !1DV case: velobc=geostrophic current
        call ogcm_1dv_0grd(tm_) !1DV case: zero x-y-gradient fields
      endif                !1dv1dv1dv>

 2000 continue ! tm_ loop

      if(obctime_order==4) then !m°v°m> !28-10-22

        temobc_t(:,:,:,obctime_aft2)= temobc_t(:,:,:,obctime_aft) 
        salobc_t(:,:,:,obctime_aft2)= salobc_t(:,:,:,obctime_aft) 
        velobc_u(:,:,:,obctime_aft2)= velobc_u(:,:,:,obctime_aft) 
        velobc_v(:,:,:,obctime_aft2)= velobc_v(:,:,:,obctime_aft) 
          sshobc_w(:,:,obctime_aft2)=   sshobc_w(:,:,obctime_aft) 
       velbarobc_u(:,:,obctime_aft2)=velbarobc_u(:,:,obctime_aft) 
       velbarobc_v(:,:,obctime_aft2)=velbarobc_v(:,:,obctime_aft) 

        temobc_t(:,:,:,obctime_bef2)= temobc_t(:,:,:,obctime_bef) 
        salobc_t(:,:,:,obctime_bef2)= salobc_t(:,:,:,obctime_bef) 
        velobc_u(:,:,:,obctime_bef2)= velobc_u(:,:,:,obctime_bef) 
        velobc_v(:,:,:,obctime_bef2)= velobc_v(:,:,:,obctime_bef) 
          sshobc_w(:,:,obctime_bef2)=   sshobc_w(:,:,obctime_bef) 
       velbarobc_u(:,:,obctime_bef2)=velbarobc_u(:,:,obctime_bef) 
       velbarobc_v(:,:,obctime_bef2)=velbarobc_v(:,:,obctime_bef) 

      endif                     !m°v°m> !28-10-22

!     call obc_ts_fluxes(-500.,-0.05,-999.,10.,2) ! ATTENTION ARG EN REAL
!     call obc_ts_fluxes(-500.,-0.05,-999.,10.,2) ! ATTENTION ARG EN REAL
!     call obc_ts_fluxes(-500.,+0.1 ,-999.,365.,2) !13-06-14


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                        PARTIE I
!                INITIALISATION DU MODELE
! FIN.
      return
      endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!                       /   /   /

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                        PARTIE II
!                MISE A JOUR DE L'ANALYSE
! appelée dans model_3D.F(au debut)
! DEBUT:
      if (case_==2) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      tm_=2 !28-10-22
      if(obctime_order==4)tm_=obctime_aft2

      do var_count_=1,var_num ! boucle sur variables: 1=T&S 2=vel, 3=ssh

      if(elapsedtime_now>ogcm_readtime_next(var_count_)     &
                         .and.subcycle_synchro==1) then !*******************> !09-02-14

! Lignes deplacee dans module_ogcm.F90 routine ogcm_get_varfilename
       ogcm_readtime_prev(var_count_)=ogcm_readtime_next(var_count_)
         ogcm_period_prev(var_count_)=  ogcm_period_next(var_count_)

           ogcm_rec_next(var_count_)=ogcm_rec_next(var_count_)+1
        nc=ogcm_rec_next(var_count_)

      if(par%rank==0) then
       if(var_count_==trc_id)write(6,*)'read_ogcm_fields traceurs nc=',nc
       if(var_count_==vel_id)write(6,*)'read_ogcm_fields vitesses nc=',nc
       if(var_count_==ssh_id)write(6,*)'read_ogcm_fields ssh nc=',nc
      endif

!......................................................
! MISE A JOUR DES TABLEAUX.
! DEBUT:
!......................................................
      if(var_count_==ssh_id) then !ssssssss>

       if(obctime_order==4) then !ooo> !28-10-22
         sshobc_w(:,:,obctime_bef2)=sshobc_w(:,:,obctime_bef)
       endif                     !ooo> !28-10-22

        do j=1,jmax ; do i=1,imax
         sshobc_w(i,j,0)=sshobc_w(i,j,2)
        enddo ; enddo

       if(obctime_order==4) then !ooo> !28-10-22
         sshobc_w(:,:,obctime_aft)=sshobc_w(:,:,obctime_aft2)
       endif                     !ooo> !28-10-22

      endif                       !ssssssss>
      if(var_count_==trc_id) then !tttttttt>

       if(obctime_order==4) then !ooo> !28-10-22
         temobc_t(:,:,:,obctime_bef2)=temobc_t(:,:,:,obctime_bef)
         salobc_t(:,:,:,obctime_bef2)=salobc_t(:,:,:,obctime_bef)
       endif                     !ooo> !28-10-22

        do k=1,kmax ; do j=1,jmax ; do i=1,imax
         temobc_t(i,j,k,0)=temobc_t(i,j,k,2)
         salobc_t(i,j,k,0)=salobc_t(i,j,k,2)
        enddo ; enddo ; enddo

       if(obctime_order==4) then !ooo> !28-10-22
         temobc_t(:,:,:,obctime_aft)=temobc_t(:,:,:,obctime_aft2)
         salobc_t(:,:,:,obctime_aft)=salobc_t(:,:,:,obctime_aft2)
       endif                     !ooo> !28-10-22

      endif                       !tttttttt>
      if(var_count_==vel_id) then !vvvvvvvv>

       if(obctime_order==4) then !ooo> !28-10-22
            velobc_u(:,:,:,obctime_bef2)= velobc_u(:,:,:,obctime_bef)
            velobc_v(:,:,:,obctime_bef2)= velobc_v(:,:,:,obctime_bef)
           velbarobc_u(:,:,obctime_bef2)=velbarobc_u(:,:,obctime_bef)
           velbarobc_v(:,:,obctime_bef2)=velbarobc_v(:,:,obctime_bef)
       endif                     !ooo> !28-10-22

        do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax+1
         velobc_u(i,j,k,0)=velobc_u(i,j,k,2)
         velobc_v(i,j,k,0)=velobc_v(i,j,k,2)
        enddo ; enddo ; enddo
        do j=1,jmax+1 ; do i=1,imax+1
         velbarobc_u(i,j,0)=velbarobc_u(i,j,2)
         velbarobc_v(i,j,0)=velbarobc_v(i,j,2)
        enddo ; enddo

       if(obctime_order==4) then !ooo> !28-10-22
            velobc_u(:,:,:,obctime_aft)= velobc_u(:,:,:,obctime_aft2)
            velobc_v(:,:,:,obctime_aft)= velobc_v(:,:,:,obctime_aft2)
           velbarobc_u(:,:,obctime_aft)=velbarobc_u(:,:,obctime_aft2)
           velbarobc_v(:,:,obctime_aft)=velbarobc_v(:,:,obctime_aft2)
       endif                     !ooo> !28-10-22

      endif                       !vvvvvvvv>
!......................................................
! MISE A JOUR DES TABLEAUX.
! FIN.
!......................................................

        if(nc.le.0)then !---- debug ----->                           !06/05/04
         write(6,*)'-----------------'
         write(6,'(a,a,a,a)')' No OGCM field available at the'        & !29-04-14
          ,' current time. The inconsistency between '         &
          ,' notebook_time and the forcing lists of'          &
          ,' notebook_obcforcing must be corrected.'
         stop ' STOP read_ogcm_fields choix 1'
        endif             !---- debug ----->

! Cas interpolation des fichiers:
      if(obc_option.eq.2)call ogcm_interp(var_count_,nc,tm_) !28-10-22

!..................................................

!.................................................
      if(var_count_==trc_id) then !ttttttttttttttt>

! densite pour ajustement convectif:
      if(eos_author==0) then !eqeqeqeqeqeqeqeqeq>  !22-11-10

! Equation d'état linéaire:
      const1=-rho*alp_t
      const2= rho*alp_s
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax+1
      x2=(const1*(temobc_t(i,j,k,tm_)-t0)       & ! densite
         +const2*(salobc_t(i,j,k,tm_)-s0))
      anyv3d(i,j,k,1)=x2
      enddo
      enddo
      enddo

      else                   !eqeqeqeqeqeqeqeqeq>   !22-11-10
! Equation d'état non linéaire:

      if(eos_author==3) then
! charger anyv3d(i,j,k,1) avec la densite potentielle:
        if(eos_tkezref>=0.) then
          call equation_of_state_potzref_jmfwg('obc',tm_) !23-02-22
        else
          call equation_of_state_potloc_jmfwg('obc',tm_) !28-10-14
        endif
      else
       stop 'read_ogcm_fields eos_author/=3 pas encore prevu'
      endif

      endif                  !eqeqeqeqeqeqeqeqeq>  !22-11-10

      if(flag_ogcm_instab==1) call ajust_convect(tm_)     !19-02-22
      if(flag_ogcm_instab==2) call ajust_convect_mix(tm_) !19-02-22
      if(flag_ogcm_instab==3) call ajust_convect_min(tm_) !19-02-22

      if(flag_ogcmtidemixing==1)call tidal_ogcm_mix(tm_) !11-12-19

      endif                       !ttttttttttttttt>

!..................................................
      if(var_count_==vel_id) then !vvvvvvvvvvvvvvv>
       do k=1,kmax
       do j=1,jmax+1
       do i=1,imax+1
        velobc_u(i,j,k,tm_)=velobc_u(i,j,k,tm_)*mask_u(i,j,k)
        velobc_v(i,j,k,tm_)=velobc_v(i,j,k,tm_)*mask_v(i,j,k)
       enddo
       enddo
       enddo
       do j=1,jmax+1
       do i=1,imax+1
        velbarobc_u(i,j,tm_)=velbarobc_u(i,j,tm_)*mask_u(i,j,kmax+1)
        velbarobc_v(i,j,tm_)=velbarobc_v(i,j,tm_)*mask_v(i,j,kmax+1)
       enddo
       enddo
      endif                       !vvvvvvvvvvvvvvv>
!..............................................................................

!     call obc_ts_fluxes(-500.,-0.05,-999.,10.,2) ! ATTENTION ARG EN REAL
!     call obc_ts_fluxes(-500.,+0.1 ,-999.,365.,2) !13-06-14

       if(flag_1dv==1) then !1d1d1d> !29-11-15
! Le modele 1DV fait l'hypothese que toutes les variables ogcm sont au meme temps afin
! de pouvoir traiter dans le meme temps, l'interpolation de T,S,ssh, le courant, le courant
! geostrophique, puis de mettre tous les champs obtenus egaux A la valeur de i=2, j=2 pour
! obtenir des champs constants compatibles avec la continuite mpi 1DV. Toutes ces operations
! se font A la fin quand var_count_=var_num
        if(var_count_==var_num) then !vvvvv> 
         call ogcm_1dv_geos(tm_) !1DV case: velobc=geostrophic current
         call ogcm_1dv_0grd(tm_) !1DV case: zero x-y-gradient fields
        endif                        !vvvvv>
       endif                !1d1d1d>

      endif                                             !*******************>
      enddo  ! fin de boucle sur var_count_  1=traceurs, 2=vitesses, 3=ssh

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                        PARTIE II
!                MISE A JOUR DE L'ANALYSE
! FIN.
      return
      endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!------------------------------------------------------------------
      if(case_==0 )stop 'choix 0 dans read_ogcm_fields supprimé'
      if(case_==-1)stop 'choix -1 dans read_ogcm_fields supprimé'

      end subroutine read_ogcm_fields
!__________________________________________________________________






!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   !28/01/05
      subroutine ajust_convect(ki1)
      use module_principal
      implicit none
      integer ki1
#ifdef synopsis
       subroutinetitle='ajust_convect'
       subroutinedescription= &
          'Returns a well-mixed T S profile whenever a positive' &
       //' vertical gradient of potential sea water density is'  &
       //' encountered'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      const1=-rho*alp_t
      const2= rho*alp_s

! En entrée ANYV3D contient un profil de densite fixé à priori.
! En sortie ANYV3D contient le même profil mais ajusté convectivement.
! Les tableaux TOBC et SOBC(i,j,k,KI1) lui sont ses temperature
! et salinité associées

! traitement de l'instabilite convective
      do j=1,jmax !26-11-14
      do i=1,imax !25-11-14
      if(  mask_t(i,j,kmax).eq.1) then !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



!.............................................>
 50   continue
      k1=0
      do k=kmax,kmin_w(i,j)+1,-1

       k5=k
       do k2=k-1,kmin_w(i,j),-1
        if(anyv3d(i,j,k2,1).gt.anyv3d(i,j,k2+1,1))then !--------->
         k5=k2+1
         goto 60
        endif                                          !--------->
        if(k2.eq.kmin_w(i,j))k5=1   ! cas instable jusqu'au fond
       enddo


   60 if(k.gt.k5.and.                                                   &
         abs(anyv3d(i,j,k,1)-anyv3d(i,j,k5,1)).gt.1.e-5)then !*******>


       x3=0.
       x2=0.
       x1=0.
       x0=0.
       do k2=k,k5,-1
        x0=x0+dz_t(i,j,k2 ,1)
        x1=x1+dz_t(i,j,k2 ,1)*anyv3d(i,j,k2,1)
        x2=x2+dz_t(i,j,k2 ,1)*  temobc_t(i,j,k2,ki1)
        x3=x3+dz_t(i,j,k2 ,1)*  salobc_t(i,j,k2,ki1)
       enddo
       do k2=k,k5,-1
        anyv3d(i,j,k2,1)    =x1/x0
          temobc_t(i,j,k2,ki1)=x2/x0
          salobc_t(i,j,k2,ki1)=x3/x0
       enddo
       k1=1

      endif                                                  !*******>



      enddo


      if(k1.eq.1) goto 50    ! Test : il restait de l'instabilite à la precedente iteration


! Pour ne pas entretenir de l'instabilité convective aux frontières ouvertes  !12/04/05
! (avec un gros niveau de turbulence) on "re-stratifie" tres legerement la colonne
! quand celle ci est parfaitement mélangée
!     do k=kmin_w(i,j)+1,kmax

! Partie supprimée car l'inversion de la densité via cette equation d'etat simple !23-11-10
! n'est plus valide si equation d'etat complete (notamment parce que la constante
! rho est une densité potentielle dans ele premier cas et une densité dans le second)
! et ce d'autant plus qu'à l'etat initial la routine initial_state_eq n'a pas encore
! fourni les parametres optimaux

!     anyv3d(i,j,k,1)=min(anyv3d(i,j,k  ,1)                             &
!                        ,anyv3d(i,j,k-1,1)                             &
!                 -1.e-6*( depth_t(i,j,k)- depth_t(i,j,k-1)))

!     temobc_t(i,j,k,ki1)=t0+(anyv3d(i,j,k,1)                             &
!                    -const2*(salobc_t(i,j,k,ki1)-s0))/const1

! Eventuellement, remplacer par:                                        !23-11-10
!     sal_t(i,j,k,ki1)=min(sal_t(i,j,k  ,ki1)                        &
!                         ,sal_t(i,j,k-1,ki1)                        &
!                 -1.e-6*( depth_t(i,j,k)- depth_t(i,j,k-1)))

!     enddo

!.............................................>


      endif                          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      enddo
      enddo

! fin du traitement de l'instabilite convective
      end subroutine ajust_convect
!___________________________________________________________________________________
      subroutine tidal_ogcm_mix(tm_)
      use module_principal ; use module_parallele
      implicit none
      integer tm_
#ifdef synopsis
       subroutinetitle='tidal_ogcm_mix'
       subroutinedescription= &
          'Computes the thickness of the tide driven bottom turbulent' &
       //' layer and mixes within the former the ogcm T & S fields'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      if(kmaxtide.eq.0)return

! Suggestion pour plus tard: l'epaisseur de la couche
! d'Ekman pourrait être deduite de Csanady 1982:
! D= 0.1*ustar/f (page 12, paragraphe 1.7)
! avec ustar la vitesse de friction = (stress/rho)**0.5
! le stress de fond etant deduit du CD et du module des
! vitesse de maree de mog2d
! Le choix entre la moyenne et le max de courant de maree
! a ete tranche en faveur du maximum


!#ifdef bidon
!------------------------------------------
! Methode maree maximale: !22-09-21
      do j=1,jmax
      do i=1,imax

      x1=0.
      do ktide=1,kmaxtide

      x1=x1                                                             &
       +((veltidecos_u(i  ,j  ,ktide)+veltidecos_u(i+1,j  ,ktide))/2.)**2       &
       +((veltidesin_u(i  ,j  ,ktide)+veltidesin_u(i+1,j  ,ktide))/2.)**2       &
       +((veltidecos_v(i  ,j  ,ktide)+veltidecos_v(i  ,j+1,ktide))/2.)**2       &
       +((veltidesin_v(i  ,j  ,ktide)+veltidesin_v(i  ,j+1,ktide))/2.)**2

      enddo

! XY_Z(I,J,1)=D la hauteur
!     xy_t(i,j,1)=0.1*sqrt(cdb_t(i,j)*sqrt(x1))      &
!     xy_t(i,j,1)=0.1*sqrt( 2.5d-3   *sqrt(x1))      & !26-04-17
      xy_t(i,j,1)=0.1*sqrt( 1.d-3    *sqrt(x1))      & !10-11-18
             /max(abs(coriolis_t(i,j)),5.d-5)          !11-08-14
! Note: le choix de cd=1.d-3 et non pas cdb_t(:,:) vient de ce que !10-11-18
! le courant de maree "estimE" est 2D et que le cisaillement vertical
! du courant ne m'est pas connu. J'emploie donc un coef de frottement
! typique d'un modele de maree 2D.

      enddo
      enddo
!#endif

#ifdef bidon
!------------------------------------------
! Methode maree moyenne au temps considErE:
! Le principe de cet algo est de calculer le courant de maree !10-11-18
! A differents intervalle de temps (2h A priori) sur une periode
! de ... de maintenant A maintenant+... et de calculer le module
! moyen pour en deduire un effet de melange moyen A
! appliquer A l'OGCM
! anyvar3d(:,:,1:2) sont les composantes u et v du courant de maree instantanE total
! anyvar2d est la moyenne temporelle sur la fenetre de ... de ce dernier
! reset:

! note: pour avoir l'effet moyen du module il n'est pas necessaire
! d'integrer sur une periode entiere de maree, une demi-periode suffit.
! Donc pour une onde diurne la moyenne sur 12h suffit. k10 explore
! donc une fenetre de 12h
      anyvar2d=0. 

      do k10=1,6 ! boucle sur 12h avec un pas de temps de 2h

! time2 va de "maintenant" jusqu'A un jour plus tard
      time2=elapsedtime_now+(k10-1)*3600.*2.

! Courant de maree total au temps time2
      do ktide=1,kmaxtide

       time1=frqtide(ktide)*(time2-ti0tide(ktide))+v0tide(ktide)+utide(ktide,1) 
       const3=cos(time1)*ftide(ktide,1)
       const4=sin(time1)*ftide(ktide,1)                          

      do j=1,jmax ; do i=1,imax+1

       anyvar3d(i,j,1)=                         &
       anyvar3d(i,j,1)*passetide(ktide,1)       &
             +( veltidecos_u(i,j,ktide)*const3  &
               +veltidesin_u(i,j,ktide)*const4)

      enddo ; enddo

      do j=1,jmax+1 ; do i=1,imax

       anyvar3d(i,j,2)=                         &
       anyvar3d(i,j,2)*passetide(ktide,1)       &
             +( veltidecos_v(i,j,ktide)*const3  &
               +veltidesin_v(i,j,ktide)*const4)

      enddo ; enddo

      enddo ! fin de boucle sur ktide

! Moyenne temporelle du module du courant
      do j=1,jmax ; do i=1,imax
       anyvar2d(i,j)=                                     &
       anyvar2d(i,j)+                                     &
       sqrt( (0.5*(anyvar3d(i,j,1)+anyvar3d(i+1,j,1)))**2 &
            +(0.5*(anyvar3d(i,j,2)+anyvar3d(i,j+1,2)))**2 )

      enddo       ; enddo

      enddo ! fin de boucle temporelle sur k10
      anyvar2d=anyvar2d/6. ! division par le nbre d'iteration temporelle
      do j=1,jmax ; do i=1,imax
! XY_Z(I,J,1)=D la hauteur
!     xy_t(i,j,1)=0.1*sqrt(cdb_t(i,j)*sqrt(x1))      &
!     xy_t(i,j,1)=0.1*sqrt( 2.5d-3   *sqrt(x1))      & !26-04-17
      xy_t(i,j,1)=0.1*sqrt( 1.d-3*anyvar2d(i,j))     & !10-11-18 
             /max(abs(coriolis_t(i,j)),5.d-5)          !11-08-14
! Note: le choix de cd=1.d-3 et non pas cdb_t(:,:) vient de ce que !10-11-18
! le courant de maree "estimE" est 2D et que le cisaillement vertical
! du courant ne m'est pas connu. J'emploie donc un coef de frottement
! typique d'un modele de maree 2D.
      enddo       ; enddo
#endif


! MELANGER LA COUCHE DE FOND
      do j=1,jmax ; do i=1,imax
        k1=0
        do k=1,kmax
         if(h_w(i,j)+depth_t(i,j,k)<xy_t(i,j,1))k1=k
        enddo
        if(k1>1) then !>>>>>
          sum1=0.
          sum2=0.
          sum3=0.
          do k=1,k1
           sum1=sum1+dz_t(i,j,k,1)
           sum2=sum2+dz_t(i,j,k,1)*temobc_t(i,j,k,tm_)
           sum3=sum3+dz_t(i,j,k,1)*salobc_t(i,j,k,tm_)
          enddo
          sum1=max(sum1,small1)
          do k=1,k1
           temobc_t(i,j,k,tm_)=sum2/sum1
           salobc_t(i,j,k,tm_)=sum3/sum1
          enddo
        endif         !>>>>>

      enddo       ; enddo

! Les lignes suivantes ont ete "bidonnees" le 03-08-19 apres avoir constate qu'elles produisaient
! un melange insuffisant (les couches fusionnees persistaient A etre stratifiees) 
! occasionnant un plantage de la simulation de l'estuaire du komo
#ifdef bidon

! Calcul du coef de diffusivite verticale induit par la maree:
!c    X1=50.     ! epaisseur de la couche limite de fond
!     x3=5.      ! epaisseur de la transition
!     x2=1000./dti_lp ! valeur max de K

! coef1 --->
      do j=1,jmax ; do i=1,imax

       x1=xy_t(i,j,1)   ! epaisseur de la couche limite de fond
       x3=0.1*x1+small1 ! epaisseur de la transition !27-10-15
!      x3=0.0*x1+small1 ! epaisseur de la transition !27-10-15
       

       do k=1,kmax ! +1

!       x2=1000.                                 & ! Valeur max de K !27-10-15
!       *( depth_t(i,j,k  )- depth_t(i,j,k-1))   &
!            /dti_lp     
!       anyv3d(i,j,k,1)=x2*(1.-tanh((h_w(i,j)+ depth_w(i,j,k)-x1)/x3))/2.
! cette variante arrete vraiment le melange au dessus de D ! !10-11-18
!       anyv3d(i,j,k,1)=max(x2*(0.99-tanh((h_w(i,j)+ depth_w(i,j,k)-x1)/x3))/2.,0.) !10-11-18

!        tridia_in(i,j,k,1)=-max(1000.*(0.99-tanh((h_w(i,j)+depth_w(i,j,k)-x1)/x3))/2.,0.)*mask_t(i,j,k)
         tridia_in(i,j,k,1)=-max(10000.*dz_t(i,j,k,1)*(0.99-tanh((h_w(i,j)+depth_w(i,j,k)-x1)/x3))/2.,0.)*mask_t(i,j,k) !02-05-19

       enddo

      enddo ; enddo


      do k=1,kmax ; do j=1,jmax ; do i=1,imax

! coef1 --->
!     tridia_in(i,j,k,1)=                                               &
!       -anyv3d(i,j,k  ,1)*dti_lp/                                      &
!       ( depth_t(i,j,k  )- depth_t(i,j,k-1))                           &
!                                    *  mask_t(i,j,k)                  !28/02/03

! coef4 --->
      tridia_in(i,j,k,4)=temobc_t(i,j,k,tm_)   &
                            *dz_t(i,j,k,1)     &
                          *mask_t(i,j,k)

      enddo ; enddo ; enddo

! coef3 --->
      do k=1,kmax-1 ; do j=1,jmax ; do i=1,imax
         tridia_in(i,j,k  ,3) &
        =tridia_in(i,j,k+1,1)
      enddo       ; enddo       ; enddo

!......................................................
! LES CONDITIONS AUX LIMITES SUR Coef1 Coef3 et Coef4:
      do j=1,jmax
      do i=1,imax
       tridia_in(i,j,1   ,1)=0. !02-05-19
       tridia_in(i,j,kmax,3)=0.
      enddo
      enddo
!......................................................

! coef2 --->
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

      tridia_in(i,j,k,2)=dz_t(i,j,k,1)-( &
                                         tridia_in(i,j,k,1)  &
                                        +tridia_in(i,j,k,3)  &
                                       )*mask_t(i,j,k)

      enddo ; enddo ; enddo

      call tridiagonalsolver(1,1,imax,1,jmax,kmax)

      do k=1,kmax ; do j=1,jmax ; do i=1,imax
      temobc_t(i,j,k,tm_)=tridia_out(i,j,k)         &
                             *mask_t(i,j,k)
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=1,imax
! coef4 --->
      tridia_in(i,j,k,4)=salobc_t(i,j,k,tm_)   &
                            *dz_t(i,j,k,1)     &
                          *mask_t(i,j,k)
      enddo ; enddo ; enddo

      call tridiagonalsolver(1,1,imax,1,jmax,kmax)

      do k=1,kmax ; do j=1,jmax ; do i=1,imax
      salobc_t(i,j,k,tm_)=tridia_out(i,j,k)     &
                             *mask_t(i,j,k)
      enddo ; enddo ; enddo
#endif

      end subroutine tidal_ogcm_mix
!___________________________________________________________________________________
#ifdef bidon
      subroutine read_ogcm_fields_geo(ki3)
      use module_principal
      implicit none
      integer ki3
#ifdef synopsis
       subroutinetitle='read_ogcm_fields_geo'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! OPTION: LE COURANT EST DEDUIT de la ssh et de T et S via l'equilibre
! hydrostatique.

! Recapitulatiof de l'utilisation des tableaux temporaire:
! Anomalie de densité:               ANYV3D(I,J,K,3)
! Anomalie de pression:              ANYV3D(I,J,K,4)
! Gradient de pression direction Ox: ANYV3D(I,J,K,1)
! Gradient de pression direction Oy: ANYV3D(I,J,K,2)

! calcul de la pression:

      const1=-rho*alp_t
      const2= rho*alp_s
      const3=-const1*t0-const2*s0
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
      anyv3d(i,j,k,3)=const1*temobc_t(i,j,k,ki3)                          &
                     +const2*salobc_t(i,j,k,ki3)                          &
                     +const3
      enddo
      enddo
      enddo

! integrale par les rectangle:
! calcul de la pression:
      do 9 j=1,jmax
      do 9 i=1,imax
      xy_t(i,j,1)=0.
    9 continue
      do 10 k=kmax,1,-1
       do 10 j=1,jmax
       do 10 i=1,imax
       x2=grav*anyv3d(i,j,k,3)*dz_t(i,j,k,1)
       xy_t(i,j,1)=xy_t(i,j,1)+x2
!ccccccPRESSURE_Z(I,J,K)=XY_Z(I,J,1)-0.5*X2
       anyv3d(i,j,k,4)=xy_t(i,j,1)-0.5*x2
   10 continue


!     calcul du gradient en x
      const1=1./(   rho)
      const2= grav/(2.*rho)

      do 14 k=kmax,1,-1
      do 14 j=1,jmax
      do 14 i=2,imax

      anyv3d(i,j,k,1)=const1*(  anyv3d(i,j,k,4)  -anyv3d(i-1,j,k,4))    &
                     +const2*( depth_t(i,j,k) -  depth_t(i-1,j,k))      &
                            *(  anyv3d(i,j,k,3) + anyv3d(i-1,j,k,3))

   14 continue

!     calcul du gradient en y
      const1=1./(   dyb*rho)
      const2= grav/(2.*dyb*rho)

      do 19 k=kmax,1,-1
      do 19 j=2,jmax
      do 19 i=1,imax

      anyv3d(i,j,k,2)=const1*(  anyv3d(i,j,k,4)-  anyv3d(i,j-1,k,4))    &
                     +const2*( depth_t(i,j,k) -  depth_t(i,j-1,k))      &
                            *(  anyv3d(i,j,k,3) + anyv3d(i,j-1,k,3))

   19 continue

        do k=1,kmax

        do j=2,jmax-1
        do i=2,imax
           if (  mask_u(i,j,k).ne.0) then
              sum1=0.
              do i1=i-1,i
              do j1=j,j+1
                 sum1=sum1+                                             &
                    (-anyv3d(i1,j1,k,2)                                 &
                     -grav*(sshobc_w(i1,j1,ki3)-sshobc_w(i1,j1-1,ki3)))    &
                   /(0.5*(coriolis_t(i1,j1)   +coriolis_t(i1,j1-1)))    &
                             /dy_v(i1,j1)                               &
                         *  mask_v(i1,j1,k)
              enddo
              enddo
              velobc_u(i,j,k,ki3)=sum1/4.
            else
              velobc_u(i,j,k,ki3)=0.
            endif
         enddo
         enddo

! OBC:
        do j=2,jmax-1
              velobc_u(imax+1,j,k,ki3)=velobc_u(imax,j,k,ki3)           &
             *  mask_u(imax+1,j,k)

              velobc_u(1     ,j,k,ki3)=velobc_u(2   ,j,k,ki3)           &
             *  mask_u(1     ,j,k)
        enddo

        do i=1,imax+1
              velobc_u(i,jmax,k,ki3)=velobc_u(i,jmax-1,k,ki3)           &
             *  mask_u(i,jmax,k)

              velobc_u(i,1   ,k,ki3)=velobc_u(i,2     ,k,ki3)           &
             *  mask_u(i,1   ,k)
        enddo


         do j=2,jmax
         do i=2,imax-1
            if (  mask_v(i,j,k).ne.0) then
               sum1=0.
               do i1=i,i+1
               do j1=j-1,j
                  sum1=sum1+                                            &
                   (+anyv3d(i1,j1,k,1)                                  &
                    +grav*(sshobc_w(i1,j1,ki3)-sshobc_w(i1-1,j1,ki3)))     &
                    /(0.5*(coriolis_t(i1,j1)   +coriolis_t(i1-1,j1)))   &
                              /dx_u(i1,j1)                              &
                          *  mask_u(i1,j1,k)
               enddo
               enddo
              velobc_v(i,j,k,ki3)=sum1/4.
            else
              velobc_v(i,j,k,ki3)=0.
            endif
         enddo
         enddo

! OBC:
         do j=2,jmax
              velobc_v(imax,j,k,ki3)=velobc_v(imax-1,j,k,ki3)           &
             *  mask_v(imax,j,k)
              velobc_v(1   ,j,k,ki3)=velobc_v(2     ,j,k,ki3)           &
             *  mask_v(1   ,j,k)
         enddo
         do i=1,imax
              velobc_v(i,jmax+1,k,ki3)=velobc_v(i,jmax,k,ki3)           &
             *  mask_v(i,jmax+1,k)
              velobc_v(i,1     ,k,ki3)=velobc_v(i,2   ,k,ki3)           &
             *  mask_v(i,1     ,k)
         enddo

         enddo

! Courant moyen:
         do j=1,jmax+1
         do i=1,imax+1
          velbarobc_u(i,j,ki3)=0.
          velbarobc_v(i,j,ki3)=0.
         enddo
         enddo
         do j=1,jmax+1
         do i=1,imax+1
          velbarobc_u(i,j,ki3)=velbarobc_u(i,j,ki3)                     &
                                 +dz_u(i,j,ki3,1)/hz_u(i,j,1)           &
                             *velobc_u(i,j,k,ki3)
          velbarobc_v(i,j,ki3)=velbarobc_v(i,j,ki3)                     &
                                 +dz_v(i,j,ki3,1)/hz_v(i,j,1)           &
                             *velobc_v(i,j,k,ki3)
         enddo
         enddo

      end subroutine read_ogcm_fields_geo

!----------------------------------------------------------------------!

      subroutine already_interpolated(nc_loc,t_loc)
      stop 'Routine already_interpolated obsolete'
#ifdef bidon
      use module_principal
      implicit none
      integer nc_loc,t_loc
#ifdef synopsis
       subroutinetitle='already_interpolated'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! longueur des enregistrements des tableaux 2D:
      lrec=(imax+1)*(jmax+1)*4 ! 27/06/01
      write(6,*)'valeur decimale du rec=',                              &
!      ( 1.+     (kount-obcafdt(2))/      obcafdt(1)  ) + real(t_loc)/2.
       (1.+(elapsedtime_now-obcafdt(2))/obcafdt(1))+real(t_loc)/2.     !17-04-11
      write(6,*)'dans read_ogcm_fields nc_loc=',nc_loc
      write(6,*)'dans read_ogcm_fields lrec 2d =',lrec

! elevation de la surface:
      write(6,'(a)')obcfile(1,1)
      open(unit=3,file=obcfile(1,1),access='direct'                     &
          ,form='unformatted'                                          & !22/01/04
          ,recl=lrec)
      read(3,rec=nc_loc)((                                                  &
       sshobc_w(i,j,t_loc)                                             &
       ,i=1,imax+1),j=1,jmax+1)
      close(3)                                                         !23/10/06

! courant moyen composante X
      write(6,'(a)')obcfile(2,1)
      open(unit=3,file=obcfile(2,1),access='direct'                     &
          ,form='unformatted'                                          & !22/01/04
          ,recl=lrec)
      read(3,rec=nc_loc)((                                                  &
       velbarobc_u(i,j,t_loc)                                          &
       ,i=1,imax+1),j=1,jmax+1)
      close(3)

! courant moyen composante Y
      write(6,'(a)')obcfile(3,1)
      open(unit=3,file=obcfile(3,1),access='direct'                     &
          ,form='unformatted'                                          & !22/01/04
          ,recl=lrec)
      read(3,rec=nc_loc)((                                                  &
       velbarobc_v(i,j,t_loc)                                          &
       ,i=1,imax+1),j=1,jmax+1)
      close(3)

! longueur des enregistrements des tableaux 3D:
      lrec=(imax+1)*(jmax+1)*(kmax)*4                                  !25/10/01
      write(6,*)'dans read_ogcm_fields lrec 3d =',lrec

      write(6,'(a)')obcfile(4,1)
      open(unit=3,file=obcfile(4,1),access='direct'                     &
          ,form='unformatted'                                          & !22/01/04
          ,recl=lrec)
      read(3,rec=nc_loc)(((                                                 &
       temobc_t(i,j,k,t_loc)                                             &
       ,i=1,imax+1),j=1,jmax+1),k=1,kmax)
      close(3)


      write(6,'(a)')obcfile(5,1)
      open(unit=3,file=obcfile(5,1),access='direct'                     &
          ,form='unformatted'                                          & !22/01/04
          ,recl=lrec)
      read(3,rec=nc_loc)(((                                                 &
       salobc_t(i,j,k,t_loc)                                             &
       ,i=1,imax+1),j=1,jmax+1),k=1,kmax)
      close(3)

! courant 3D VRAI composante X
      write(6,'(a)')obcfile(6,1)
      open(unit=3,file=obcfile(6,1),access='direct'                     &
          ,form='unformatted'                                          & !22/01/04
          ,recl=lrec)
      read(3,rec=nc_loc)(((                                                 &
       velobc_u(i,j,k,t_loc)                                           &
       ,i=1,imax+1),j=1,jmax+1),k=1,kmax)
      close(3)

! courant 3D VRAI composante Y
      write(6,'(a)')obcfile(7,1)
      open(unit=3,file=obcfile(7,1),access='direct'                     &
          ,form='unformatted'                                          & !22/01/04
          ,recl=lrec)
      read(3,rec=nc_loc)(((                                                 &
       velobc_v(i,j,k,t_loc)                                           &
       ,i=1,imax+1),j=1,jmax+1),k=1,kmax)
      close(3)

#endif
      end subroutine already_interpolated
#endif
!----------------------------------------------------------------------!
#ifdef bidon
      subroutine obc_ts_fluxes(z_,dtem_,dsal_,time_,t_) !13-06-14
      use module_principal
      use module_parallele !#mpi
      use module_systeme
      implicit none
      double precision lagrange_flux_t_,lagrange_flux_m_   &
                      ,lagrange_flux_s_
      real z_,dtem_,dsal_,time_
      integer status_argument_,mski1_,mski2_,mskj1_   &
      ,mskj2_,t_,loop_
      double precision, dimension(:,:), allocatable :: matcoef_
#ifdef synopsis
       subroutinetitle='obc_ts_fluxes'
       subroutinedescription=  &
          'Ajusts the entering lateral boundary current in order' &
       //' to "control" (the big word!) the T and/or S budget'
 
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


! Cette routine modifie le champs de vitesse velobc pour
! modifier les flux de temperature et de sel au travers des
! frontieres ouvertes sans modifier le volume d'eau global
! Les arguments:
! z_:    epaisseur sur laquelle est modifie le courant
! dtem_: increment moyen de temperature sur le domaine total
! dsal_: increment moyen de salinite sur le domaine total
! time_: echelle de temps en jours associee a l'increment

! Voir détails:
! https://docs.google.com/document/d/13fgO20h4d5LJaN3Q4nEUqLOSVaTlySXjQmi9WHk8C4s/edit

      status_argument_=0
      if(z_>zero)stop 'obc_ts_fluxes z_ doit etre negatif'
      status_argument_=1

! Par défaut hypothèse proc non bordé par une frontière ouverte -> msk...=0
      mski1_=0 ; mski2_=0 ; mskj1_=0 ; mskj2_=0
! Si possibilite de frontière ouverte (voir tout de mEme valeur de mask_t) le msk concerné=1
       if(obcstatus(ieq1   )==1)mski1_=1 ! frontière ouest (i=1) ouverte
       if(obcstatus(ieqimax)==1)mski2_=1 ! frontière est (i=imax) ouverte
       if(obcstatus(jeq1   )==1)mskj1_=1 ! frontière sud (j=1) ouverte
       if(obcstatus(jeqjmax)==1)mskj2_=1 ! frontière nord (j=jmax) ouverte

! La pertubation du courant est de la forme: dz*dy*f(z)*sign(obc)*(alpha*T+beta)
! On satisfait 2 contraintes.
! La contrainte sur le flux d'eau est SOMME( u*dz*dy*sign(obc) )=0 où sign(OBC)=+-1 selon
! que l'axe de grille pointe en dehors (-1) ou en dedans (+1) du domaine
! La contrainte du flux de chaleur est SOMME(u*dz*dy*sign(obc)*T)=V*dtem_/time_
       sum1=0. ; sum2=0. ; sum3=0. ; sum4=0. ; sum5=0. ; sum6=0.
       do i=1,imax,imax-1
        if(i==1   )then !-->
         k0=mski1_      ! (1 ou 0)
        endif           !-->
        if(i==imax)then !..>
         k0=-mski2_     ! (-1 ou 0)
        endif           !..>
        do k=1,kmax
        do j=1,jmax

! f(z)*sign_obc
          anyv3d(i,j,k,0)=max(zero,1-depth_t(i,j,k)/z_)             &
          *sqrt( (velobc_u(i,j,k,t_)+velobc_u(i+1,j  ,k,t_))**2     &
                +(velobc_v(i,j,k,t_)+velobc_v(i  ,j+1,k,t_))**2)*k0


          x1=dy_t(i,j)*dz_t(i,j,k,now)*anyv3d(i,j,k,0)*mask_t(i,j,k) &
                                                      *mask_j_w(j)

          sum1=sum1+x1
          sum2=sum2+x1*temobc_t(i,j,k,t_)
          sum3=sum3+x1*temobc_t(i,j,k,t_)**2
          sum4=sum4+x1*salobc_t(i,j,k,t_)
          sum5=sum5+x1*salobc_t(i,j,k,t_)**2
          sum6=sum6+x1*temobc_t(i,j,k,t_)*salobc_t(i,j,k,t_)

        enddo ! j
        enddo ! k

       enddo ! i

       do j=1,jmax,jmax-1
       if(j==1   )then !----->
         k0=mskj1_     ! (1 ou 0)
       endif           !----->
       if(j==jmax)then !----->
         k0=-mskj2_    ! (-1 ou 0)
       endif           !----->
        do k=1,kmax
        do i=1,imax

! Comme le coin est en commun utiliser un autre 4eme argument (1) dans anyv3d pour ne
! pas ecraser l'autre:
          anyv3d(i,j,k,1)=max(zero,1-depth_t(i,j,k)/z_)             &
          *sqrt( (velobc_u(i,j,k,t_)+velobc_u(i+1,j  ,k,t_))**2     &
                +(velobc_v(i,j,k,t_)+velobc_v(i  ,j+1,k,t_))**2)*k0

          x1=dx_t(i,j)*dz_t(i,j,k,now)*anyv3d(i,j,k,1)*mask_t(i,j,k) &
                                                      *mask_i_w(i)

          sum1=sum1+x1
          sum2=sum2+x1*temobc_t(i,j,k,t_)
          sum3=sum3+x1*temobc_t(i,j,k,t_)**2
          sum4=sum4+x1*salobc_t(i,j,k,t_)
          sum5=sum5+x1*salobc_t(i,j,k,t_)**2
          sum6=sum6+x1*temobc_t(i,j,k,t_)*salobc_t(i,j,k,t_)

       enddo ! i
       enddo ! k
       enddo ! j

#ifdef parallele
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum3,sum3glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum4,sum4glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum5,sum5glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum6,sum6glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
#else
      sum1glb=sum1
      sum2glb=sum2
      sum3glb=sum3
      sum4glb=sum4
      sum5glb=sum5
      sum6glb=sum6
#endif

! x1=Determinant principal:
!     x1=sum3glb*sum1glb-sum2glb*sum2glb
!     lagrange_flux_t_= dtem_*grid_volumeglb/(time_*86400.)*sum1glb/x1
!     lagrange_flux_m_=-dtem_*grid_volumeglb/(time_*86400.)*sum2glb/x1
!     write(6,*)'SOLUTION A ',lagrange_flux_t_,lagrange_flux_m_

      if(dsal_==-999.0.or.dtem_==-999.0) then !nnnnn>
       nbequa=2 ; nbinco=2 ; nbsparse=4
      else                                    !nnnnn>
       nbequa=3 ; nbinco=3 ; nbsparse=9
      endif                                   !nnnnn>

      allocate(matcoef_(nbinco,nbinco))

      call allocate_systeme(1,        & ! Flag "allocate"
                            nbequa,   & ! Nombre d'equations (>= nb d'inconnues)
                            nbsparse, & ! Nombre d'elements non nuls dans la matrice creuse
                            nbinco  )   ! Nombre d'inconnues

      if(dsal_==-999.) then !>>>>>>>
       matcoef_(1,1)=sum2glb ; matcoef_(1,2)=sum1glb ; matrix_b(1)=0.
       matcoef_(2,1)=sum3glb ; matcoef_(2,2)=sum2glb ; matrix_b(2)=dtem_*grid_volumeglb/(time_*86400.)
      endif                 !>>>>>>>
      if(dtem_==-999.) then !>>>>>>>
       matcoef_(1,1)=sum4glb ; matcoef_(1,2)=sum1glb ; matrix_b(1)=0.
       matcoef_(2,1)=sum5glb ; matcoef_(2,2)=sum4glb ; matrix_b(2)=dsal_*grid_volumeglb/(time_*86400.)
      endif                 !>>>>>>>
      if(dsal_/=-999.and.  &
         dtem_/=-999.) then !>>>>>>>
       matcoef_(1,1)=sum2glb ; matcoef_(1,2)=sum4glb ; matcoef_(1,3)=sum1glb ; matrix_b(1)=0.
       matcoef_(2,1)=sum3glb ; matcoef_(2,2)=sum6glb ; matcoef_(2,3)=sum2glb ; matrix_b(2)=dtem_*grid_volumeglb/(time_*86400.)
       matcoef_(3,1)=sum6glb ; matcoef_(3,2)=sum5glb ; matcoef_(3,3)=sum4glb ; matrix_b(3)=dsal_*grid_volumeglb/(time_*86400.)
      endif                 !>>>>>>>


      ksparse=0
      do kequation=1,nbequa
      do kinconnue=1,nbinco
           ksparse=ksparse+1
           isparse(ksparse)=kequation
           jsparse(ksparse)=kinconnue
       sparsevalue(ksparse)=matcoef_(kequation,kinconnue)
      enddo
      enddo

      call leastsquaresolver(0,0,0)

      if(dsal_==-999.) then !>>>>>>>
       lagrange_flux_t_=tab7(1)
       lagrange_flux_m_=tab7(2)
       lagrange_flux_s_=0.
      endif                 !>>>>>>>
      if(dtem_==-999.) then !>>>>>>>
       lagrange_flux_s_=tab7(1)
       lagrange_flux_m_=tab7(2)
       lagrange_flux_t_=0.
      endif                 !>>>>>>>
      if(dsal_/=-999.and.  &
         dtem_/=-999.) then !>>>>>>>
       lagrange_flux_t_=tab7(1)
       lagrange_flux_s_=tab7(2)
       lagrange_flux_m_=tab7(3)
      endif                 !>>>>>>>

!     write(6,*)'SOLUTION B ',lagrange_flux_t_,lagrange_flux_m_

      call allocate_systeme(2,0,0,0)
      deallocate(matcoef_)

        sum1=0. ; sum2=0. ; sum3=0.
        do i=1,imax,imax-1
         if(i==1) then
          i1=1 ; i3=1 ; i2=i1+i3*nint(sponge_l)
         endif
         if(i==imax) then
          i1=imax+1 ; i3=-1 ; i2=i1+i3*nint(sponge_l)
         endif

         do j=1,jmax
         do i4=i1,i2,i3
          velbarobc_u(i4,j,t_)=0.
         enddo
         enddo
         do k=1,kmax
          do j=1,jmax
! x1 est la correction du courant au point t reportée sur les 2 points u contigus:
          x1=anyv3d(i,j,k,0)*(lagrange_flux_t_*temobc_t(i,j,k,t_)    &
                             +lagrange_flux_s_*salobc_t(i,j,k,t_)    &
                             +lagrange_flux_m_)                      &
            *mask_u(i1,j,k) ! masque sur i1=1 ou i1=imax+1 sinon A cause des eponges !20-04-16
                            ! on corrige aussi des points A cOtE de frontieres fermEes

            do i4=i1,i2,i3

             velobc_u(i4,j,k,t_)=                      &
             velobc_u(i4,j,k,t_)                       &
              +mask_u(i4,j,k)*x1*min((1.-i3*real(i4-i1-i3)/sponge_l),un)
! Note: la fonction de decroissance=1 en i=1 et i=2 pour obc

             velbarobc_u(i4,j,t_)=                     &
             velbarobc_u(i4,j,t_)+velobc_u(i4,j,k,t_)  &
                                     *dz_u(i4,j,k,now) !20-04-16

            enddo

! Pour verification
           sum1=sum1+x1*dy_t(i,j)*dz_t(i,j,k,1)*mask_t(i,j,k)*mask_j_w(j)
           sum2=sum2+x1*dy_t(i,j)*dz_t(i,j,k,1)*mask_t(i,j,k)*mask_j_w(j)*temobc_t(i,j,k,t_)
           sum3=sum3+x1*dy_t(i,j)*dz_t(i,j,k,1)*mask_t(i,j,k)*mask_j_w(j)*salobc_t(i,j,k,t_)

          enddo
         enddo
         do j=1,jmax
         do i4=i1,i2,i3
          velbarobc_u(i4,j,t_)=velbarobc_u(i4,j,t_)/hz_u(i4,j,1)
         enddo
         enddo

        enddo  ! loop i


        do j=1,jmax,jmax-1
         if(j==1) then
          j1=1 ; j3=1 ; j2=j1+j3*nint(sponge_l)
         endif
         if(j==jmax) then
          j1=jmax+1 ; j3=-1 ; j2=j1+j3*nint(sponge_l)
         endif

         do i=1,imax
         do j4=j1,j2,j3
          velbarobc_v(i,j4,t_)=0.
         enddo
         enddo

         do k=1,kmax
          do i=1,imax
! x1 est la correction du courant au point t reportée sur les 2 points v contigus:
          x1=anyv3d(i,j,k,1)*(lagrange_flux_t_*temobc_t(i,j,k,t_)    &
                             +lagrange_flux_s_*salobc_t(i,j,k,t_)    &
                             +lagrange_flux_m_)                      &
            *mask_v(i,j1,k) ! masque sur i1=1 ou i1=imax+1 sinon A cause des eponges !20-04-16
                            ! on corrige aussi des points A cOtE de frontieres fermEes

            do j4=j1,j2,j3

             velobc_v(i,j4,k,t_)=                      &
             velobc_v(i,j4,k,t_)                       &
              +mask_v(i,j4,k)*x1*min(1.-j3*real(j4-j1-j3)/sponge_l,un)

            velbarobc_v(i,j4,t_)=                     &
            velbarobc_v(i,j4,t_)+velobc_v(i,j4,k,t_)  &
                                    *dz_v(i,j4,k,now)

            enddo

! Pour verification
            sum1=sum1+x1*dx_t(i,j)*dz_t(i,j,k,1)*mask_t(i,j,k)*mask_i_w(i)
            sum2=sum2+x1*dx_t(i,j)*dz_t(i,j,k,1)*mask_t(i,j,k)*mask_i_w(i)*temobc_t(i,j,k,t_)
            sum3=sum3+x1*dx_t(i,j)*dz_t(i,j,k,1)*mask_t(i,j,k)*mask_i_w(i)*salobc_t(i,j,k,t_)
        enddo
        enddo
         do i=1,imax
         do j4=j1,j2,j3
          velbarobc_v(i,j4,t_)=velbarobc_v(i,j4,t_)/hz_v(i,j4,1)
         enddo
         enddo

        enddo ! loop j

! Pour verification
#ifdef parallele
      call mpi_allreduce(sum3,sum3glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
#else
      sum3glb=sum3
      sum2glb=sum2
      sum1glb=sum1
#endif
      if(par%rank==0) then
       write(6,*)'sum1glb=',sum1glb
       write(6,*)'sum2glb=',sum2glb,dtem_*grid_volumeglb/(time_*86400.)
       write(6,*)'sum3glb=',sum3glb,dsal_*grid_volumeglb/(time_*86400.)
      endif

      if(status_argument_==0)stop 'arg error in obc_ts_fluxes'
!#ifdef parallele
!     call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#!endif
!     stop 'SENEGAL'
      end subroutine obc_ts_fluxes
#endif
!.......................................................................

      subroutine ajust_convect_mix(tm_)
      use module_principal
      use module_parallele
      implicit none
      integer tm_,flag2_

! En entrée ANYV3D contient un profil de densite fixé à priori.
! En sortie ANYV3D contient le même profil mais ajusté convectivement.
! Les tableaux TOBC et SOBC(i,j,k,KI1) lui sont ses temperature
! et salinité associées

! traitement de l'instabilite convective
      do j=1,jmax ; do i=1,imax 
      if(mask_t(i,j,kmax)==1) then !>>>

       flag2_=0 ! ne sera pas remis A zero 
  47   flag=0   ! est remis A zero
       do k=kmax,kmin_w(i,j)+1,-1

        if(anyv3d(i,j,k,1)-anyv3d(i,j,k-1,1)>1.e-10)then !--------->
! Instabilite detectee
        k1=k-1 ! base (provisoire) de la couche instable

        flag2_=1 ! ne sera pas remis A zero
        flag=1
             sum0= &
               dz_t(i,j,k  ,1)                      &
              +dz_t(i,j,k-1,1) 
             sum1=                                  &
               dz_t(i,j,k  ,1)*anyv3d(i,j,k  ,1)    &
              +dz_t(i,j,k-1,1)*anyv3d(i,j,k-1,1)     

! OU est la base de la couche instable ?
! La densite melangEe sum1/sum0 est elle > ou < A la densite de la 
! couche en dessous ?
! Cette boucle en suivant pour repondre A ces questions:
         do k3=k-1,kmin_w(i,j)+1,-1
          if(sum1/sum0-anyv3d(i,j,k3-1,1)>1.e-10) then  !pmx>
! instable:
           k1=k3-1 ! mise A jour position de la base de la couche instable
           sum0=sum0+dz_t(i,j,k3-1,1)
           sum1=sum1+dz_t(i,j,k3-1,1)*anyv3d(i,j,k3-1,1)
          else                                      !pmx>
! stable:  
           goto 48  
          endif                                     !pmx>
         enddo

! Finaliser le melange entre le niveau k (sommet couche instable) et le niveau k1 (base)
   48    continue

!        sum1=0. ! facultatif car dEjA calculE
         sum2=0. 
         sum3=0.
         do k4=k1,k
!          sum1=sum1+dz_t(i,j,k4,1)  *anyv3d(i,j,k4,1) ! facultatif
           sum2=sum2+dz_t(i,j,k4,1)*temobc_t(i,j,k4,tm_)
           sum3=sum3+dz_t(i,j,k4,1)*salobc_t(i,j,k4,tm_)
         enddo
         do k4=k1,k
             anyv3d(i,j,k4,1)  =sum1/sum0
           temobc_t(i,j,k4,tm_)=sum2/sum0
           salobc_t(i,j,k4,tm_)=sum3/sum0
         enddo

        endif                                          !--------->

       enddo ! boucle k
       if(flag==1) goto 47

       if(flag2_==1) then !-renforcer-la -stabilite->
! flag2_=1 indique que l'algo d'ajustement A EtE appliquE.
! Parce que l'algo ne garantit pas strictement la neutralite
! en raison d'imprecisions d'origines diverses, on ajoute une
! surstratification arbitraire sur salobc et temobc. Cette correction
! est la plus petite possible, en tenant compte de ce que ces tableaux 
! sont en simple precision:
        do k=1,kmax
         salobc_t(i,j,k,tm_)=salobc_t(i,j,k,tm_)-0.001*depth_t(i,j,k)/hmax
         temobc_t(i,j,k,tm_)=temobc_t(i,j,k,tm_)+0.001*depth_t(i,j,k)/hmax
        enddo
       endif              !-renforcer-la -stabilite->
      

      endif                        !>>> ! test sur mask_t
      enddo ; enddo

! fin du traitement de l'instabilite convective
      end subroutine ajust_convect_mix

!.......................................................................

      subroutine ajust_convect_min(tm_)
      use module_principal
      use module_parallele
      implicit none
      integer tm_

! traitement simple de l'instabilite convective
      do k=kmax-1,1,-1
       do j=1,jmax ; do i=1,imax 

        if(anyv3d(i,j,k+1,1)>anyv3d(i,j,k,1))then !--------->
           anyv3d(i,j,k,1)    =anyv3d(i,j,k+1,1)
         temobc_t(i,j,k,tm_)=temobc_t(i,j,k+1,tm_) 
         salobc_t(i,j,k,tm_)=salobc_t(i,j,k+1,tm_) 
        endif                                     !--------->

       enddo       ; enddo
      enddo

      end subroutine ajust_convect_min

!.......................................................................
