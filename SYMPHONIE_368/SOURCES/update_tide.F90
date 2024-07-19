      subroutine update_tide
!______________________________________________________________________
! SYMPHONIE ocean model
! release 301 - last update: 21-03-21
!...............................................................................
!  _________                    .__                  .__              ! m°v°m
! /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____  
! \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \       !
! /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/       !
!/_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >      !
!        \/\/          \/|__|        \/            \/        \/       !
!...............................................................................
      use module_principal ; use module_subcycle
      implicit none
#ifdef synopsis
       subroutinetitle='update_tide'
       subroutinedescription='Updates tidal forcing terms'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date      Description des modifications
!         01/01/02: version initiale
!         15/01/02: introduction d'un choix entre un forcage sur toutes les
!                   frontieres et un forcage sur les frontieres seulement
!                   où les ondes sont incidentes.
!         21/01/02: corrections de bugs + nouveau cas IN_OUT_TIDE=3
!         24/01/02: introduction potentiel + effet de la maree sur elle-meme
!         22/02/02: amenagements 3D PASSETIDE augmente de dimension. Voir
!                   egalement initial_tide
!                   appels aux routines d'analyse
!         25/02/02: reorganisation des boucles pour gagner de la cpu
!                   prise en compte de la maree sur le courant 3D total
!         26/08/02: VBAROBC remplacé par VMEAOBC
!                   & VHZOBC remplacé par VELOBC
!         30/04/04: remettre la rampe en courant X
!         28/07/04: On sort immediatement de la routine si ONOFF_TIDE=0
!                   Suppression de lignes commentées
!                   Amenagement de la filiere analyse harmonique
!         28/12/05: debug rampe externe: dte_lp remplace dti_fw
!         10/05/06: fonctions compatibles avec double precision
!         26/03/07: rampe "commune" calculée dans chronos
!                   cos(time1) et sin(time1) passés dans const3 et const4
!                   pour economie de cpu
!                   On ne calcule pas les obc (mais le reste si) en cas
!                   d'imbrication symphonie -> symphonie
!         18/01/08: La solution totale de SLA est stockee. But: soustraction
!                   d'une moyenne de reference degagee de l'effet de maree dans
!                   la routine botprescont.F
!         23/12/08: On recalle correctement la mise à jour de la conditions aux
!                   limites compte tenu du schéma d'interpolation temporelle dans
!                   obc_ext.F. En pratique la C.L. est calculée pour le temps t+1
!         06/05/09: TI0TIDE
!         12-05-09: tendance de la moyenne pour contrainte dans obc_ext.F
! 2009.3  02-10-09: dte_lp_lp remplace dte_lp
!         08-11-09: optimobc_tide est renommé tide_analysis
! 2010.2  21-12-09: - parametres nodaux apparaissent explicitement dans les calculs
!                   contrairement à la version precedente etaient une fois pour
!                   toute figés dans les conditions initiales
!                   - modif dimensions parametres nodaux
! 2010.20 16-04-11  Calculs sur la base d'un temps en secondes
! 2010.25 26-02-12  suppression onoff_tide
! S26     06-11-13  - ajout tideanalysis_nextime (requis subcycling)
!                   - subcycling
!         15-04-16  ti0tide devient tableau 1DV
!         09-09-16  tableau sshtotide inutile
!         20-03-17  suprresion sshmeataide
! v301    21-03-21  avec le calcul mpi, spo_i et spo_j sont inutiles
!...............................................................................

      if(kmaxtide.eq.0)return                                        !28/07/04

! une fois par an lire les parametres nodaux suivants:
!     if(kount>kounttide_next_rdv)call tide_nodal_parameters           !22-12-09
      if(elapsedtime_now>tidenodal_next_rdv)call tide_nodal_parameters !16-04-11


! interpoler lineairement les parametres nodaux entre 2 echeances annuelles:
!     rap=(       real(kount)-kounttide_prev_rdv)                       &
!        /(kounttide_next_rdv-kounttide_prev_rdv)
      rap=(   elapsedtime_now-tidenodal_prev_rdv)           & ! 16-04-11
         /(tidenodal_next_rdv-tidenodal_prev_rdv)

      do ktide=1,kmaxtide
       utide(ktide,1)=(1.-rap)*utide(ktide,0)+rap*utide(ktide,2)
       ftide(ktide,1)=(1.-rap)*ftide(ktide,0)+rap*ftide(ktide,2)
      enddo


! Note: pour le moment on ne traite pas la composante 3D du courant car
! à priori cela est inutile. A voir.

! Provisoire: PENSER A ENLEVER CES LIGNES!
! pleine grille:
!     sshmeatide(0)=sshmeatide(1)                                      !16-05-09
      do 999 ktide=1,kmaxtide

! mise à jour du potentiel:                                            !24/01/02
! Attention la rampe sur le potentiel est appliqué à la fin
! de la boucle 999 (economise cpu) et sous certaines conditions
! (pas de rampe dans le cas d'une imbrication symphonie->symphonie)

! Temps t:  (pour potentiel)
      time1=frqtide(ktide)*(elapsedtime_now-ti0tide(ktide))           & !06/05/09
           +v0tide(ktide)+utide(ktide,1)                                !21-12-09
      const3=cos(time1)*ftide(ktide,1)                                  !24-12-09
      const4=sin(time1)*ftide(ktide,1)                                  !24-12-09
      do j=1,jmax
       do i=1,imax
        tidepotential_w(i,j,1)=                                           &
        tidepotential_w(i,j,1)*passetide(ktide,2)                         & !24-12-09
       +potidecos_w(i,j,ktide)*const3                                 &
       +potidesin_w(i,j,ktide)*const4
!       apotide_z(i,j,ktide)*cos(time1+ppotide_z(i,j,ktide))            &
!                           *ftide(ktide,1)                             !21-12-09
       enddo
      enddo

! forcage sur frontieres ondes incidentes:
      if(in_out_tide.eq.0) then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>
      write(6,*)                                                        &
       'forcage sur frontieres ondes incidentes n''est pas pret'
      write(6,*)'modifier notebook_tide en consequence'
      stop ' dans update_tide!'
! forcage sur frontieres ondes incidentes:
      endif                     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%>


! Forcages aux conditions aux limites (à condiions que ce ne soit pas
!                                      un model_ imbrique)
      if(nest_onoff_in.ne.1) then !**********************************> !26/03/07

! Temps t+1:  (pour C.L. barotrope)
!     time1=frqtide(ktide)*((kount+1-ti0tide)*dti_fw+0.5*dte_lp)      & !06/05/09
      time1=frqtide(ktide)*((elapsedtime_aft-ti0tide(ktide))+0.5*dte_lp) & !15-04-16
            +v0tide(ktide)+utide(ktide,1)                               !21-12-09
      const3=cos(time1)*rampe                                         & !23/12/08
                       *ftide(ktide,1)                                  !21-12-09
      const4=sin(time1)*rampe                                         &
                       *ftide(ktide,1)                                  !21-12-09

!     do 888 loop1=1,4 ! boucle sur 4 frontières
! Avec mpi cela n'a plus de sens d'utiliser spo_i et spo_j et c'est meme
! potentiellement genant avec l'option rappel en fonction de la
! resolution et de grilles compliquees comme la bipolaire mondiale (voir
! application Ile de Java). PAr consequent la boucle sur loop1 est
! remplacee par loop1=5 (rappel sur toute la grille)
      loop1=5 !21-03-21

!................................>
! COMPOSANTE X DU TRANSPORT:
      do j=spo_j1_u(loop1),spo_j2_u(loop1)
      do i=spo_i1_u(loop1),spo_i2_u(loop1)

       velbarobc_u(i,j,1)=                                              &
       velbarobc_u(i,j,1)*passetide(ktide,1)+                           &
         (veltidecos_u(i,j,ktide)*const3+ & ! verifier qu'à l'initialisation
          veltidesin_u(i,j,ktide)*const4) ! U1 et U2 sont bien * par MORPHO

      enddo
      enddo

!................................>
! COMPOSANTE Y DU TRANSPORT:
      do j=spo_j1_v(loop1),spo_j2_v(loop1)
      do i=spo_i1_v(loop1),spo_i2_v(loop1)

        velbarobc_v(i,j,1)=                                             &
        velbarobc_v(i,j,1)*passetide(ktide,1)+ & ! verifier qu'à l'initialisation
          (veltidecos_v(i,j,ktide)*const3+   & ! V1 et V2 sont bien
           veltidesin_v(i,j,ktide)*const4)   ! multiplies par MORPHO

      enddo
      enddo

!................................>
! ELEVATION DE LA SURFACE

      do j=spo_j1_t(loop1),spo_j2_t(loop1)
      do i=spo_i1_t(loop1),spo_i2_t(loop1)

        sshobc_w(i,j,1)=                                                &
        sshobc_w(i,j,1)*passetide(ktide,1)+                             &
       (sshtidecos_w(i,j,ktide)*const3+                                     &
        sshtidesin_w(i,j,ktide)*const4)

! La solution totale de SLA est stockee. Usage: soustraction d'une moyenne
! de reference degagee de l'effet de maree dans la routine botprescon.F
!        sshtotide_w(i,j,1)=sshobc_w(i,j,1)                            !18/01/08

      enddo
      enddo

! 888 continue

!................................>
! Les moyennes de SLA sur tout le domaine:                             !18/01/08
! Moyenne simple:
!           sshmeatide(1)=                                              &
!           sshmeatide(1)*passetide(ktide,1)+                           &
!      z1meatide(ktide,1)*const3+                                       &
!      z2meatide(ktide,1)*const4
! Moyenne ponderee par H:
!           sshmeatide(2)=                                              &
!           sshmeatide(2)*passetide(ktide,1)+                           &
!      z1meatide(ktide,2)*const3+                                       &
!      z2meatide(ktide,2)*const4

! forcage sur toutes les frontieres: factorisation cosinus et sinus
      endif                     !**************************************>

  999 continue
! ZTATIDE_CUM est la tendance du cumul pour terme de correction dans
! obc_ext.F =  d(ssh_cumul)/dt avec ssh_cumul=ssh_averaged*grid_area
!     sshtide_cum=grid_area*(sshmeatide(1)-sshmeatide(0))/dti_now  !16-04-11

      if(nest_onoff_in.ne.1)then !>>>>>>>
       if(rampe.lt.un)then        !------>
       do j=1,jmax
       do i=1,imax
        tidepotential_w(i,j,1)=tidepotential_w(i,j,1)*rampe
       enddo
       enddo
      endif                       !------>
      endif                      !>>>>>>>

!*******************************************************************************
! CAS TRIDIMENSIONNEL CALCUL DU COURANT TOTAL
! DEBUT:
      if(i2dh.eq.1) then
      if(nest_onoff_in.ne.1) then                                      !26/03/07
!*******************************************************************************

!     do 200 loop1=1,4
! Avec mpi cela n'a plus de sens d'utiliser spo_i et spo_j et c'est meme
! potentiellement genant avec l'option rappel en fonction de la
! resolution et de grilles compliquees comme la bipolaire mondiale (voir
! application Ile de Java). PAr consequent la boucle sur loop1 est
! remplacee par loop1=5 (rappel sur toute la grille)
      loop1=5 !21-03-21

      do 200 k=1,kmax

!................................>
! COMPOSANTE X DU COURANT:
      do j=spo_j1_u(loop1),spo_j2_u(loop1)
      do i=spo_i1_u(loop1),spo_i2_u(loop1)

        velobc_u(i,j,k,1)=velbarobc_u(i,j,1)

      enddo
      enddo

!................................>
! COMPOSANTE Y DU COURANT:
      do j=spo_j1_v(loop1),spo_j2_v(loop1)
      do i=spo_i1_v(loop1),spo_i2_v(loop1)

        velobc_v(i,j,k,1)=velbarobc_v(i,j,1)

      enddo
      enddo

  200 continue

!*******************************************************************************
! CAS TRIDIMENSIONNEL CALCUL DU COURANT TOTAL
! FIN.
      endif
      endif
!*******************************************************************************


!_______________________________________________________________________________
! APPEL AUX PROGRAMMES D'ANALYSE EN HARMONIQUES
! DEBUT:
      if(tideana_yesno==1) then                                        !28/07/04
!_______________________________________________________________________________

      if(subcycle_synchro==1) then                   !>>>>>>>           !06-11-13
       if(elapsedtime_now>tideana_nextime) then      !------>           !06-11-13

            call tide_analysis(4)                                       !16-04-11
            tideana_nextime=tideana_nextime+tideana_modulo

       endif                                         !------>           !06-11-13
      endif                                          !>>>>>>>           !06-11-13

!_______________________________________________________________________________
! APPEL AUX PROGRAMMES D'ANALYSE EN HARMONIQUES
! FIN.
      endif
!_______________________________________________________________________________


      end subroutine update_tide
