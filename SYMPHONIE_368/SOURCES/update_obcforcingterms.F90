      subroutine update_obcforcingterms(initstatus_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 358 - last update: 28-10-22
!...............................................................................
!  _________                    .__                  .__              ! m°v°m
! /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____  
! \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \       !
! /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/       !
!/_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >      !
!        \/\/          \/|__|        \/            \/        \/       !
!...............................................................................

      use module_principal

      use module_principal
      use module_parallele
      implicit none
      integer initstatus_,flag1_,flag2_,var_count_,flag_tide_
#ifdef synopsis
       subroutinetitle='update_obcforcingterms'
       subroutinedescription= &
       'Linear time interpolation of the external forcing T,S,U,V,SSH' &
       //' terms and stores the results in velbarobc velobc temobc '   &
       //'salobc sshobc (last index = 1)'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!...............................................................................
! Version Date      Description des modifications
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         26/02/02: amenagement pour ajouter le courant de maree au courant
!                   thermohalin.
!                   + reamenagement des boucles pour economiser du temps cpu
!         26/08/02: VBAROBC remplace par VMEAOBC
!                  & VHZOBC remplace par VELOBC
!         10/09/02: cas initial different du cas standart
!         16/12/02: Ajout d'une securite sur l'extension de la zone de mise
!                   e jour de la salinite et de la temperature.
!                   Mise e jour des commentaires.
!         31/08/03: amenagement pour recevoir nouvelle filiere nest_inout.F
!                   et le cas oe la maree a ete calculee par un model_
!                   "symphonie maman"
!         01/11/04: Ajout du barometre inverse
!         19/07/05: utilisation des tableaux SPO_J1_X et cie pour calculer les
!                   boucles
!         14/03/06: destruction de lignes inutilisees depuis longtemps
!                   Maintenant il faut choisir en filiere obc_wv ou read_ogcm_fields car
!                   on ne cumule plus les 2.
!         28/03/06: RAP est garde en memoire grace a RAP_OBC
!         26/03/07: BI_ONOFF (barometre inverse) est initialise dans set_parameters.F
!                   On n'ajoute pas la maree si imbrication symphonie->symphonie
!         05/06/07: prevoir le cas oe "il n'y a pas" d'ogcm: pour que T et S ne soient
!                   pas nuls X0=1 et non 0.
!         18/06/07: mis e jour pas de temps du champs de reference.
!         23/01/08: Renouvellement des read_ogcm_fields par modulo d'entier supprime
! 2010.2  27-12-09  dtvar_obc renomme update_obcforcingterms
! 2010.7  04-03-10  prise en compte ssh d'equilibre de Stokes, et courant moyen de Stokes
! 2010.11 16-07-10  temobc & salobc renommes temobc & salobc
! 2010.20 15-04-11  Renouvellement des fichiers de forcage sur la base d'un temps
!                   en secondes
! 2010.23 17-05-11  - Ajout de commentaires
!                   - Prise en compte velstoke y compris pour la vitesse du mode interne
! S.26    29-11-13  interpolation temporelle de l'ogcm: possibilite de fichiers
!                   irregulierement archives dans le temps
!         09-02-14  Blindage mpi subcycling
!         03-11-14  Maree non definie en phase d'initialisation
!         28-11-15  modif calcul timeweightobc
!         19-09-16  ajout offset_sshobc
! v299    21-03-21  avec le calcul mpi, spo_i et spo_j sont inutiles
! v345    07-04-22  offset_sshobc deplacE dans module_ogcm le 07-04-22
! v358    28-10-22  obctime_order==4
!...............................................................................

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                MISE A JOUR DE L'ANALYSE
!                     AU COURS DU TEMPS
! Ce programme interpole au temps t les differentes echeances
! DEBUT:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      var_count_=1 ! reset defaut tous modeles

      x2=0.                                                            !31/08/03
      x0=1.                                                            !05/06/07

      if(iobc_ogcm==1)then !eeeeeeeeeeeeeeeeeeeeeeeeeeeeeee>           !31/08/03

! Pourcentage des 2 echeances: !25-11-13
! Le pourcentage n'est pas lineaire mais depend du poid de
! chaque echeance, ce dernier etant donne par la duree de la
! moyenne de chaque champs, sachant qu'avec S l'archivage peut
! etre resserre sur des periodes choisies

      do var_count_=1,var_num ! boucles sur variables traceurs, vitesses, ssh

! Cet algo est dangereux si l'utilisateur n'est pas rigoureux pour la construction
! des listes.
#ifdef bidon
       if(ogcm_period_prev(var_count_)<ogcm_period_next(var_count_)) then !>>>>>>              !29-11-13

        x1=ogcm_readtime_prev(var_count_)+ogcm_period_prev(var_count_)

        timeweightobc(var_count_)=                            &
             1.-max( (x1-elapsedtime_now  )                   &
                    /(x1-ogcm_readtime_prev(var_count_)),zero)

       else                                                             !>>>>>>

        x1=      ogcm_period_prev(var_count_)   &
           +2.*ogcm_readtime_prev(var_count_)   &
              -ogcm_readtime_next(var_count_)

        timeweightobc(var_count_)=                                     &
           1.-min( (ogcm_readtime_next(var_count_)-elapsedtime_now)    &
                  /(ogcm_readtime_next(var_count_)-x1             ),un)

       endif                                                           !>>>>>>

! Bornage car l'attente de la synchro peut faire sortir timeweightobc(var_count_) de l'intervale [0:1]
               timeweightobc(var_count_)=                &
       min(max(timeweightobc(var_count_),zero),un)             !09-02-14

#endif
! On prefere par consequent cet algo plus robuste mais moins "fin" que le precedent
       timeweightobc(var_count_)=(elapsedtime_now-ogcm_readtime_prev(var_count_)) & !28-11-15
                 /(ogcm_readtime_next(var_count_)-ogcm_readtime_prev(var_count_))


       dtobc(var_count_)=ogcm_readtime_next(var_count_)  &
                        -ogcm_readtime_prev(var_count_)

      enddo ! var_count_ fin de boucle sur variables traceurs, vitesses, ssh

      endif                !eeeeeeeeeeeeeeeeeeeeeeeeeeeeeee>           !31/08/03

      if(nest_onoff_in.eq.1)stop 'filiere nest_onoff_in supprimee'     !15-04-11


      flag_tide_=1                                                            !26/03/07
!     if(kmaxtide==0.or.nest_onoff_in.eq.1)flag_tide_=0
      if(kmaxtide==0)   flag_tide_=0
      if(initstatus_==0)flag_tide_=0 ! Maree non definie en phase d'initialisation!03-11-14


      const1=-bi_onoff/grav/rho                                        !01/11/04

!**********************************************************************!10/09/02
! HORS INITIALISATION:
! Debut:
      if(initstatus_.eq.0.or.initstatus_.eq.1) then
!**********************************************************************!10/09/02

! Concerning the Stokes velocity and ssh. The boundary condition makes the hypothesis
! of a stationnary eulerian current. The momentum equation is thus supposed to be
! -f(uhat+us)=-gd(ssh-sshs)/dx
! from which (and supposing no other source of external forcing) we assume that
! a staionnary solution is uhat=-us and ssh=sshs. As a consequence, the velbarobc
! velocity is set to -velbarstokes + other forcing velocities and the sshobc elevation
! is set to sshstokes + other forcing ssh...
! As far as the ogcm velocity is concerned, one has to distinguish 3 cases.
! Case 1 (wave_obc_type=1): the ogcm provides eulerian velocities and no wave effect
! Case 2: the ogcm provides lagrangian velocities (wave effect considered)
! Case 3: the ogcm provides eulerian velocities and the wave effect considered

! In Case 1, we assume that the ogcm provides eulerian velocity because the wave
! effect is not included in its physic. The ogcm provides uhat(T,S) where
! (T,S) means that T and S is the major forcing in the OGCM.
! In Case 2 and as the OBC concern the eulerian velocity, the ogcm velocities
! should be corrected in order to represent an eulerian velocity field. A solution
! (applied here) is to remove the Stokes fields of the present run from the OBC
! field. The ogcm provides uL=uhat(T,S)+uhat(wave)+uS. If we remove uS from uL we
! obtain uhat(T,S)+uhat(wave), that is a consistent forcing term for the eulerian
! velocity of the model.
! In case I, the ogcm provides uhat(T,S) only and thus we need to add uhat(wave)
! and the latter is assumed to be -uS. Case 1 and Case 2 both lead to substract
! the Stokes Velocity.
! Concerning the Sea Surface Elevation:
! In Case 1, the ogcm provides ssh(T,S) and thus we need to add ssh(wave)
! and the latter is assumed to be sshstokes.
! In Case 2, the ogcm provides ssh(T,S)+ssh(wave).
! Consequently Case 1 leads to add sshstokes but not Case 2.
! In case 3 both sshstokes and stokes currents are ignored by the OBC.

      flag1_=0  ; flag2_=0
      if(wave_obc_type==1) then !18-05-11
       flag1_=1  ; flag2_=1  ! Case 1 No wave effect in the ogcm eulerian current
      endif
      if(wave_obc_type==2) then !18-05-11
       flag1_=1  ; flag2_=0  ! Case 2 wave effect in the ogcm Lagrangian current
      endif
      if(wave_obc_type==3) then !18-05-11
       flag1_=0  ; flag2_=0  ! Case 3 wave effect in the ogcm Eulerian current
      endif



! Explication sur valeur loop1 dans boucle 100
! Si initstatus_=0 cas initial domaine entier tableaux eponges avec argument unique 5
! Si initstatus_=1 cas obc par obc argument tableaux eponge va de 1 e 4
!     do 100 loop1=initstatus_*1+(1-initstatus_)*5                 & !19/07/05
!                 ,initstatus_*4+(1-initstatus_)*5                   !19/07/05


! Avec mpi cela n'a plus de sens d'utiliser spo_i et spo_j et c'est meme
! potentiellement genant avec l'option rappel en fonction de la
! resolution et de grilles compliquees comme la bipolaire mondiale (voir
! application Ile de Java). PAr consequent la boucle sur loop1 est
! remplacee par loop1=5 (rappel sur toute la grille)
      loop1=5 !21-03-21


      if(obctime_order==2) then !-obctime_order=2-> !28-10-22

      x2=timeweightobc(vel_id) ; x0=1.-x2

!................................>
! COMPOSANTE X DU TRANSPORT:
      do j=1,jmax ; do i=1,imax+1 !28-10-22

      velbarobc_u(i,j,1)=x2*velbarobc_u(i,j,2)           &
                        +x0*velbarobc_u(i,j,0)           &
                +flag_tide_*velbarobc_u(i,j,1)           & ! tide
                        -velbarstokes_u(i,j,1)*flag1_   ! Stokes current !04-03-10 !18-05-11

      enddo ; enddo

!................................>
! COMPOSANTE Y DU TRANSPORT:
      do j=1,jmax+1 ; do i=1,imax !28-10-22

      velbarobc_v(i,j,1)=x2*velbarobc_v(i,j,2)  &
                        +x0*velbarobc_v(i,j,0)  &
                +flag_tide_*velbarobc_v(i,j,1)  & ! tide
                        -velbarstokes_v(i,j,1)*flag1_   ! Stokes current !04-03-10 !18-05-11

      enddo ; enddo


!................................>
! ELEVATION DE LA SURFACE
      x2=timeweightobc(ssh_id) ; x0=1.-x2
      do j=1,jmax ; do i=1,imax !28-10-22

      sshobc_w(i,j,1)=x2*sshobc_w(i,j,2)              &
                     +x0*sshobc_w(i,j,0)              &
            +flag_tide_* sshobc_w(i,j,1)              & !terme de maree             !26/02/02
                   +const1*(pss_w(i,j,1)-pss_mean(1)) & !01/11/04
                     +sshstokes_w(i,j)*flag2_           ! Stokes SSH !04-03-10 !18-05-11
! offset_sshobc deplacE dans module_ogcm le 07-04-22
!                    +offset_sshobc                     !19-09-16

      enddo ; enddo

!................................>
! COMPOSANTE X DU COURANT:
      x2=timeweightobc(vel_id) ; x0=1.-x2
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1 !28-10-22

      velobc_u(i,j,k,1)=x2*velobc_u(i,j,k,2)  &
                       +x0*velobc_u(i,j,k,0)  &
               +flag_tide_*velobc_u(i,j,k,1)  &   !terme de maree
                       -velstokes_u(i,j,k,1)*flag1_   ! Stokes current !04-03-10 !18-05-11

      enddo ; enddo ; enddo

!................................>
! COMPOSANTE Y DU COURANT:
      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax !28-10-22

      velobc_v(i,j,k,1)=x2*velobc_v(i,j,k,2)  &
                       +x0*velobc_v(i,j,k,0)  &
               +flag_tide_*velobc_v(i,j,k,1)  &  !terme de maree
                       -velstokes_v(i,j,k,1)*flag1_   ! Stokes current !04-03-10 !18-05-11

      enddo ; enddo ; enddo

!................................>
! TEMPERATURE ET SALINITE
      x2=timeweightobc(trc_id) ; x0=1.-x2
      do k=1,kmax ; do j=1,jmax ; do i=1,imax !28-10-22

      temobc_t(i,j,k,1)=x2*temobc_t(i,j,k,2)      &
                       +x0*temobc_t(i,j,k,0)

      salobc_t(i,j,k,1)=x2*salobc_t(i,j,k,2)      &
                       +x0*salobc_t(i,j,k,0)

      enddo ; enddo ; enddo

      endif                     !-obctime_order=2-> !28-10-22

      if(obctime_order==4) then !-obctime_order=4-> !28-10-22

      time_r4=timeweightobc(vel_id)
      x0_r4=- time_r4   *(time_r4-1)*(time_r4-2)*0.166666666666667  
      x1_r4= (time_r4+1)*(time_r4-1)*(time_r4-2)*0.5                
      x2_r4=-(time_r4+1)* time_r4   *(time_r4-2)*0.5                
      x3_r4= (time_r4+1)* time_r4   *(time_r4-1)*0.166666666666667  

!................................>
! COMPOSANTE X
      do j=1,jmax ; do i=1,imax+1 !28-10-22

       velbarobc_u(i,j,1)=                                &
                      x0_r4*velbarobc_u(i,j,obctime_bef2) &
                     +x1_r4*velbarobc_u(i,j,obctime_bef)  &
                     +x2_r4*velbarobc_u(i,j,obctime_aft)  &
                     +x3_r4*velbarobc_u(i,j,obctime_aft2) &
                +flag_tide_*velbarobc_u(i,j,1)           & ! tide
                        -velbarstokes_u(i,j,1)*flag1_   ! Stokes current !04-03-10 !18-05-11

      enddo ; enddo
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1 !28-10-22

       velobc_u(i,j,k,1)=                               &
                     x0_r4*velobc_u(i,j,k,obctime_bef2) &
                    +x1_r4*velobc_u(i,j,k,obctime_bef)  &
                    +x2_r4*velobc_u(i,j,k,obctime_aft)  &
                    +x3_r4*velobc_u(i,j,k,obctime_aft2) &
               +flag_tide_*velobc_u(i,j,k,1)  &   !terme de maree
                       -velstokes_u(i,j,k,1)*flag1_   ! Stokes current !04-03-10 !18-05-11

      enddo ; enddo ; enddo

!................................>
! COMPOSANTE Y
      do j=1,jmax+1 ; do i=1,imax !28-10-22

       velbarobc_v(i,j,1)=                                &
                      x0_r4*velbarobc_v(i,j,obctime_bef2) &
                     +x1_r4*velbarobc_v(i,j,obctime_bef)  &
                     +x2_r4*velbarobc_v(i,j,obctime_aft)  &
                     +x3_r4*velbarobc_v(i,j,obctime_aft2) &
                +flag_tide_*velbarobc_v(i,j,1)  & ! tide
                        -velbarstokes_v(i,j,1)*flag1_   ! Stokes current !04-03-10 !18-05-11

      enddo ; enddo
      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax !28-10-22

       velobc_v(i,j,k,1)=                               &
                     x0_r4*velobc_v(i,j,k,obctime_bef2) &
                    +x1_r4*velobc_v(i,j,k,obctime_bef)  &
                    +x2_r4*velobc_v(i,j,k,obctime_aft)  &
                    +x3_r4*velobc_v(i,j,k,obctime_aft2) &
               +flag_tide_*velobc_v(i,j,k,1)  &  !terme de maree
                       -velstokes_v(i,j,k,1)*flag1_   ! Stokes current !04-03-10 !18-05-11

      enddo ; enddo ; enddo


!................................>
! ELEVATION DE LA SURFACE
      time_r4=timeweightobc(ssh_id)
      x0_r4=- time_r4   *(time_r4-1)*(time_r4-2)*0.166666666666667  
      x1_r4= (time_r4+1)*(time_r4-1)*(time_r4-2)*0.5                
      x2_r4=-(time_r4+1)* time_r4   *(time_r4-2)*0.5                
      x3_r4= (time_r4+1)* time_r4   *(time_r4-1)*0.166666666666667  

      do j=1,jmax ; do i=1,imax !28-10-22

       sshobc_w(i,j,1)=                             &
                   x0_r4*sshobc_w(i,j,obctime_bef2) &
                  +x1_r4*sshobc_w(i,j,obctime_bef)  &
                  +x2_r4*sshobc_w(i,j,obctime_aft)  &
                  +x3_r4*sshobc_w(i,j,obctime_aft2) &
            +flag_tide_* sshobc_w(i,j,1)              & !terme de maree             !26/02/02
                   +const1*(pss_w(i,j,1)-pss_mean(1)) & !01/11/04
                     +sshstokes_w(i,j)*flag2_           ! Stokes SSH !04-03-10 !18-05-11
! offset_sshobc deplacE dans module_ogcm le 07-04-22
!                    +offset_sshobc                     !19-09-16

      enddo ; enddo

!................................>
! TEMPERATURE ET SALINITE
      time_r4=timeweightobc(trc_id)
      x0_r4=- time_r4   *(time_r4-1)*(time_r4-2)*0.166666666666667  
      x1_r4= (time_r4+1)*(time_r4-1)*(time_r4-2)*0.5                
      x2_r4=-(time_r4+1)* time_r4   *(time_r4-2)*0.5                
      x3_r4= (time_r4+1)* time_r4   *(time_r4-1)*0.166666666666667  

      do k=1,kmax ; do j=1,jmax ; do i=1,imax !28-10-22

       temobc_t(i,j,k,1)=                               &
                     x0_r4*temobc_t(i,j,k,obctime_bef2) &
                    +x1_r4*temobc_t(i,j,k,obctime_bef)  &
                    +x2_r4*temobc_t(i,j,k,obctime_aft)  &
                    +x3_r4*temobc_t(i,j,k,obctime_aft2)  

       salobc_t(i,j,k,1)=                               &
                     x0_r4*salobc_t(i,j,k,obctime_bef2) &
                    +x1_r4*salobc_t(i,j,k,obctime_bef)  &
                    +x2_r4*salobc_t(i,j,k,obctime_aft)  &
                    +x3_r4*salobc_t(i,j,k,obctime_aft2)  

      enddo ; enddo ; enddo

      endif                     !-obctime_order=4-> !28-10-22

! 100 continue ! fin de boucle sur loop1

!**********************************************************************!10/09/02
! HORS INITIALISATION:
! FIN
      return
      endif
!**********************************************************************!10/09/02


!                           /     /     /


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                MISE A JOUR DE L'ANALYSE
!                     AU COURS DU TEMPS
! Ce programme interpole au temps t les differentes echeances
! FIN.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end subroutine update_obcforcingterms
