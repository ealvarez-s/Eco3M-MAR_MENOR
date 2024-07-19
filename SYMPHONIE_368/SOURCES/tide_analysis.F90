      subroutine tide_analysis(case_)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 339 - last update: 24-03-22
!______________________________________________________________________

      use module_principal
      use module_parallele !#mpi
      use module_systeme
      implicit none
      integer case_,nbite_tide,nbtot_tide,nbobc_tide
#ifdef synopsis
       subroutinetitle='tide_analysis'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!...............................................................................
! Version date      Description des modifications
!         02/01/02: version initiale
!         26/01/02: bienvenue à ICHOIX=0
!         01/03/02: debug ATAN2
!         05/03/02: passage à norme f95
!         26/08/02: VBAROBC remplacé par VMEAOBC
!                   prise en compte de ce que U1TIDE_X et cie sont maintenant
!                   des composantes de courant moyen et non plus transport
!         14/04/03: SMALL2 remplace SMALL1
!         01/10/03: on ne cree plus le graphique dans ce sous programme
!         09/07/04: amenagement pour ne pas calculer le produit de la matrice
!                   principale par sa transposée pour chaque point de grille
!                   (puisque celle ci est la même en tous points)
!         10/07/04: Le calcul de la matrice principale ne se fait qu'une fois.
!         28/07/04: On sort immediatement de la routine si ONOFF_TIDE=0
!                   Suppression des parties qui ne servaient plus...
!                   Le tableau ANALYSETIDE est mis à zéro dans initial_tide et non
!                   plus dans ICHOIX.EQ.5
!                   On sort immediatement de la routine si TIDEANA_YESNO=0
!                   On analyse également les composante du courant moyen de la
!                   marée
!         15/01/08: sortir l'analyse dans des fichiers semblables au forcage
!                   t-ugo2d
!         20/05/08: Ne pas ecraser time.dat quand on repart d'un fichier restart
!         10/05/09: Sortir l'analyse dans un fichier netcdf Amplitude, retard
!                   de Phase
! 2009.3  13-10-09: sorties netcdf: ajout de la grille 4
!         16-10-09: ecriture fichier pour version parallele
!         17-10-09: suite...
!         08-11-09: - ecriture fichier revue
!                   - optimobc_tide est renommé tide_analysis
!         16-12-09: plus de rigueur dans la logique des noms des fichiers
!                   ecrivant en netcdf
!         21-12-09: modif dimensions parametres nodaux
! 2010.2  24-12-09: la routine tide_netcdf_create_file de creation de fichier est
!                   dans une routine separee
!         27-12-09: - nettoyage des case_ inutilisés
!                   - modif dans section "Remove the nodal factors"
! 2010.8  19-05-10  En raison du changement de nom de la variable ssh du mode
!                   barotrope on choisit de faire porter l'analyse sur ssh_int
!                   au temps t=2 pour pouvoir fonctionner en mode nr<=3
! 2010.10 16-06-10  leastsquaresolver renommé leastsquaresolver
! 2010.11 27-06-10  - stockage info pour accelerer execution
!                   - En cas d'analyse sur une periode trop courte la proximité
!                   des frequences amène à un systeme instable (equations
!                   redondantes): on ajoute des equations supplementaires
!                   stabilisatrices grace au renforcement de la diagonale principale
!                   - time.dat devient dom_c//time.dat
!         01-07-10  modif sur time2 & time.dat ecrit dans tmp
!         18-07-10  archivage dans fichier time.dat remplacé par un calcul en
!                   ligne dans tideanalysismatrix
!         28-07-10  l'analyse porte sur ssh_int_w(i,j,2) car dans le cas d'une
!                   modelisation 2d, ssh_int_w(i,j,1) n'est pas determiné
! 2010.20 16-04-11  Calculs sur la base d'un temps en seconde
! 2010.24 19-11-11  L'analyse harmonique comporte une possibilité de rappel vers
!                   la solution de forcage
! 2010.25 26-02-12  Suppression onoff_tide
! S26     25-02-13  Suppression choix 5
!         11-11-14  Introduction d'un compteur de champs analyses decidant si 
!                   le nombre est suffisant pour calculer l'analyse
!         06-02-15  Calcul de l'analyse de la SSH sur 0:imax+1 0:jmax+1
!         15-04-16  ti0tide devient tableau 1DV
!         23-09-16  possibilite (si decommentE) d'analyser:
!                   ssh_int_w(i,j,1)+0.5*rhpzavr_w(i,j)*hz_w(i,j,1)/rho plutot que
!                   ssh_int_w pour mieux separer maree externe et maree interne
! v257    03-07-19  analyse harmonique ssh: multiplication wetmask 
! v276    03-04-20  analyse harmonique 3D
! v282    08-05-20  ne pas ajouter la densite 2D dans l'analyse 3d
!                   Ajout de l'analyse de la moyenne temporelle
! v296    24-02-21  signal maree interne dans la SSH
! v302    30-04-21  utilisation velbot_u et velbot_v
! v339    24-03-22  l'analyse harmonique porte sur rho_t la densitE utilisEe par le 
!                   gradient de pression
!...............................................................................
!    _________                    .__                  .__             !m°v°w
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................



      if(kmaxtidep1==0)return                                        !28/07/04
      if(tideana_yesno==0)return                                   !28/07/04

!.............................................................................
! CALCUL PROGRESSIF DU PRODUIT DE lA MATRICE PRINCIPALE PAR LA MATRICE EGALITE
      if(case_==4) then                             !oooooo>
       if(subcycle_synchro==1) then                  !>>>>>>>
        if(elapsedtime_now>tideana_nextime) then      !------> 

! Analyse en harmoniques du signal en elevation de la surface:
! pour eviter un stokage volumineux et gagner de la cpu le produit de la
! transposée de la matrice principale par la matrice egalite est réalisé
! progressivement au cours de la simulation.

      tideanalysis_count=tideanalysis_count+1 !11-11-14


      kinconnue=0
      do 426 ktide=1,kmaxtidep1

      time2=frqtide(ktide)*(elapsedtime_now-ti0tide(ktide))    &  !15-04-16
            +v0tide(ktide)+utide(ktide,1)

! Cas particulier de l'harmonique ktide=kmaxtide+1 qui correspond A la frequence nulle
! cette derniere n'etant pas compatible avec sin(time2) tout le temps nul car cela conduit la matrice
! principale A contenir une colonne et une ligne 100% nulle. On applique donc un protocole particulier:
! un coup sur 2 (cos(time2),sin(time2))=(1,0) et le reste du temps =(0,1)
! ETAPE 1: (voir ETAPE 2 dans le calcul du produit de la transposee)
       if(ktide==kmaxtide+1) then !m°v°m> !09-05-20
         if(modulo(tideanalysis_count,2)==0) then !>>>
          time2=0.
         else                                     !>>>
          time2=0.5*pi
         endif                                    !>>>
       endif                      !m°v°m>

      do 426 k1=1,2
       kinconnue=kinconnue+1

       x1=(     mod(k1,2) *cos(time2)                    &
           +(1.-mod(k1,2))*sin(time2) )*ftide(ktide,1)

! comme la routine est appelee avant external_mode on peut considerer
! que velbar(1) et ssh_int(2) sont la solution au t "now" d'où
! le calcul de time2 basé sur kount (et non pas kount+1 ou autre...)
      do j=1,jmax   ; do i=1,imax+1 !06-02-15
        analysetide_u(i,j,kinconnue)=analysetide_u(i,j,kinconnue)     & !28/07/04
         +x1*velbar_u(i,j,1)                                            !27-06-10
      enddo       ; enddo

      do j=1,jmax+1 ; do i=1,imax   !06-02-15
        analysetide_v(i,j,kinconnue)=analysetide_v(i,j,kinconnue)     & !28/07/04
         +x1*velbar_v(i,j,1)
      enddo       ; enddo


      if(kmax==1) then !pmxpmx> !15-09-16

! CAS simulation 2D:
       do j=0,jmax+1 ; do i=0,imax+1   !06-02-15
        analysetide_w(i,j,kinconnue)=analysetide_w(i,j,kinconnue)     &
        +x1*ssh_int_w(i,j,2)  ! ssh_int(i,j,1) n'existe pas en mode 2d !28-07-10
       enddo       ; enddo

      else             !pmxpmx> !15-09-16

! CAS simulation 3D:
       do j=0,jmax+1 ; do i=0,imax+1   !06-02-15

        analysetide_w(i,j,kinconnue)=analysetide_w(i,j,kinconnue)     &
        +x1*ssh_int_w(i,j,1)                                          &
           *wetmask_t(i,j) !03-07-19

       enddo       ; enddo

! Ici on analyse separement la part du signal de SSH qui est associEe Aux gradients de densite !24-02-21
! Explications dans:
! https://docs.google.com/document/d/1_mf-UA8AZqasnL-ifM5JWyst7YnqqU0KoVyoswwFcQc/edit
       if(rhp_zavr_xy==1) then !>>>
        x2=x1*(-0.5/rho)
        do j=2,jmax-1 ; do i=2,imax-1   !06-02-15
         analysetide2_w(i,j,kinconnue)=analysetide2_w(i,j,kinconnue)     &
                                        +x2*rhpzavr_w(i,j)*hz_w(i,j,1)
!       +x1*(                -0.5*rhpzavr_w(i,j)*hz_w(i,j,1)/rho)
        enddo       ; enddo
       endif                   !>>>

! Analyse du courant 3d et de la densite potentielle:
      if(flag_tide3d_analysis==1) then !m°v°m> 03-04-20
!      do k=1,kmax ; do j=1,jmax   ; do i=1,imax+1 
!       analysetide3d_u(i,j,k,kinconnue)=analysetide3d_u(i,j,k,kinconnue)+x1*vel_u(i,j,k,1)                  
!      enddo       ; enddo         ; enddo
!      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax 
!       analysetide3d_v(i,j,k,kinconnue)=analysetide3d_v(i,j,k,kinconnue)+x1*vel_v(i,j,k,1)
!      enddo       ; enddo       ; enddo
       do j=1,jmax   ; do i=1,imax+1 
        do k=1,kmerged_u(i,j)
         analysetide3d_u(i,j,k,kinconnue)=analysetide3d_u(i,j,k,kinconnue)+x1*velbot_u(i,j) !30-04-21
        enddo
        do k=kmerged_u(i,j)+1,kmax
         analysetide3d_u(i,j,k,kinconnue)=analysetide3d_u(i,j,k,kinconnue)+x1*vel_u(i,j,k,1)                  
        enddo
       enddo         ; enddo
       do j=1,jmax+1 ; do i=1,imax 
        do k=1,kmerged_v(i,j)
         analysetide3d_v(i,j,k,kinconnue)=analysetide3d_v(i,j,k,kinconnue)+x1*velbot_v(i,j) !30-04-21
        enddo
        do k=kmerged_v(i,j)+1,kmax
         analysetide3d_v(i,j,k,kinconnue)=analysetide3d_v(i,j,k,kinconnue)+x1*vel_v(i,j,k,1)
        enddo
       enddo       ; enddo
       do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1  
        analysetide3d_t(i,j,k,kinconnue)=analysetide3d_t(i,j,k,kinconnue)+x1*rho_t(i,j,k) !24-03-22
       enddo       ; enddo         ; enddo
      endif                            !m°v°m> 03-04-20

      endif            !pmxpmx> !15-09-16

  426 continue



! Calcul progressif du produit de la matrice principale par sa transposee: !18-07-10
      do kinconnue=1,2*kmaxtidep1
      do kequation=1,2*kmaxtidep1

      k1=int(real(kequation+1)/2.)
      k2=int(real(kinconnue+1)/2.)

      time1=frqtide(k1)*(elapsedtime_now-ti0tide(k1))            &  !16-04-11
            +v0tide(k1)+utide(k1,1)

      time2=frqtide(k2)*(elapsedtime_now-ti0tide(k2))            &  !16-04-11
            +v0tide(k2)+utide(k2,1)

! Cas particulier de l'harmonique ktide=kmaxtide+1 qui correspond A la frequence nulle
! cette derniere n'etant pas compatible avec sin(time2) tout le temps nul car cela conduit la matrice
! principale A contenir une colonne et une ligne 100% nulle. On applique donc un protocole particulier:
! un coup sur 2 (cos(time2),sin(time2))=(1,0) et le reste du temps =(0,1)
! ETAPE 1: (voir ETAPE 2 dans le calcul du produit de la transposee)
       if(k2==kmaxtide+1) then !m°v°m> !09-05-20
         if(modulo(tideanalysis_count,2)==0) then !>>>
          time2=0.
         else                                     !>>>
          time2=0.5*pi
         endif                                    !>>>
       endif                   !m°v°m>
       if(k1==kmaxtide+1) then !m°v°m> !09-05-20
         if(modulo(tideanalysis_count,2)==0) then !>>>
          time1=0.
         else                                     !>>>
          time1=0.5*pi
         endif                                    !>>>
       endif                   !m°v°m>

      tideanalysismatrix(kequation,kinconnue)=                         &
      tideanalysismatrix(kequation,kinconnue)+                         &
       (mod(kequation,2)*cos(time1)+(1.-mod(kequation,2))*sin(time1))*ftide(k1,1) &
      *(mod(kinconnue,2)*cos(time2)+(1.-mod(kinconnue,2))*sin(time2))*ftide(k2,1)


      enddo
      enddo

! Determiner le temps du prochain champs A analyser:
      tideana_nextime=tideana_nextime+tideana_modulo

      return
      endif                               !------> 
      endif                              !>>>>>>>
      endif                             !oooooo>
!...................................................................

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! DECOMPOSER LE SIGNAL EN HARMONIQUES PRODUIT TRANSPOSE MATRICE EGALITE
! REALISE AU COURS DE LA SIMULATION (ANALYSETIDE_Z)
! DEBUT:
      if(case_==6) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      call allocate_systeme(1,1,1,2*kmaxtidep1)            !21/01/08

      do kinconnue=1,2*kmaxtidep1                          !18-07-10
      do kequation=1,2*kmaxtidep1
       tab3(kequation,kinconnue)=tideanalysismatrix(kequation,kinconnue)
      enddo
      enddo

!     do kequation=1,2*kmaxtidep1
!      write(10+par%rank,*)(real(tab3(kequation,kinconnue)),kinconnue=1,2*kmaxtidep1)
!     enddo

      ksecu=0                                                          !09/07/04
      nbinco=2*kmaxtidep1

      if(par%rank==0)write(6,'(a,a,i0,a)')'Analyse harmonique pour SSH' &
      ,' basee sur ',tideanalysis_count,' echantillons'
      do 900 j3=0,jmax+1 !06-02-15
      do 900 i3=0,imax+1

      if(  mask_t(i3,j3,kmax+1).eq.1) then !££££££££££££££££££££££££>

            kinconnue=0
            do ktide=1,kmaxtidep1
            do k1=1,2
               kinconnue=kinconnue+1

                       tab4(kinconnue)=analysetide_w(i3,j3,kinconnue)

            enddo
            enddo

!.....................................................................
! Au premier coup (KSECU=1) on calcule le produit par la transposée
! Au deuxieme coup (KSECU.GE.2) on relie le resultat:
                           ksecu=ksecu+1                               !09/07/04
                           call leastsquaresolver(min0(ksecu,2),0,0)   !09/07/04
!.....................................................................

                           do 903 ktide=1,kmaxtidep1
                            k1=(ktide-1)*2+1
                            k2=(ktide-1)*2+2
                            sshtidecosout_w(i3,j3,ktide)=tab7(k1)
                            sshtidesinout_w(i3,j3,ktide)=tab7(k2)

  903                      continue


      endif                            !££££££££££££££££££££££££>
  900 continue


      if(kmax>1.and.rhp_zavr_xy==1) then !3DCASE> !24-02-21
      if(par%rank==0)write(6,'(a,a,i0,a)')'Analyse harmonique SSH3D' &
      ,' basee sur ',tideanalysis_count,' echantillons'
      do j3=2,jmax-1 ; do i3=2,imax-1
      if( mask_t(i3,j3,kmax+1).eq.1) then !-LSM-LSM-LSM>

            kinconnue=0
            do ktide=1,kmaxtidep1
            do k1=1,2
               kinconnue=kinconnue+1
                    tab4(kinconnue)=analysetide2_w(i3,j3,kinconnue)
            enddo
            enddo

                           call leastsquaresolver(min0(ksecu,2),0,0)   !09/07/04

                           do ktide=1,kmaxtidep1
                            k1=(ktide-1)*2+1
                            k2=(ktide-1)*2+2
                            ssh3dtidecosout_w(i3,j3,ktide)=tab7(k1)
                            ssh3dtidesinout_w(i3,j3,ktide)=tab7(k2)
                           enddo


      endif                               !-LSM-LSM-LSM>
      enddo ; enddo
!...........
! Finaliser: faire la continuite mpi 
       do ktide=1,kmaxtidep1
          xy_t(:,:,1:2)=0.
          do j=2,jmax-1 ; do i=2,imax-1
           xy_t(i,j,1)=ssh3dtidecosout_w(i,j,ktide) ; xy_t(i,j,2)=ssh3dtidesinout_w(i,j,ktide)
          enddo         ; enddo
          call obc_ext_xy_t('za',1) ; call obc_ext_xy_t('za',2)
          do j=0,jmax+1 ; do i=0,imax+1
           ssh3dtidecosout_w(i,j,ktide)=xy_t(i,j,1) ; ssh3dtidesinout_w(i,j,ktide)=xy_t(i,j,2)
          enddo         ; enddo
        enddo
!...........
      endif           !3DCASE>


      if(par%rank==0)write(6,'(a,a,i0,a)')'Analyse harmonique pour U' &
      ,' basee sur ',tideanalysis_count,' echantillons'
!     write(6,*)'Analyse harmonique pour U ',par%rank
      do j3=1,jmax
      do i3=1,imax+1
      if(  mask_u(i3,j3,kmax+1).eq.1) then !££££££££££££££££££££££££>

            kinconnue=0
            do ktide=1,kmaxtidep1
            do k1=1,2
               kinconnue=kinconnue+1

                       tab4(kinconnue)=analysetide_u(i3,j3,kinconnue) 

            enddo
            enddo
            call leastsquaresolver(2,0,0)

                           do ktide=1,kmaxtidep1
                             k1=(ktide-1)*2+1
                             k2=(ktide-1)*2+2
                             veltidecosout_u(i3,j3,ktide)=tab7(k1)
                             veltidesinout_u(i3,j3,ktide)=tab7(k2)
                           enddo


      endif                            !££££££££££££££££££££££££>

      enddo
      enddo

      if(par%rank==0)write(6,'(a,a,i0,a)')'Analyse harmonique pour V' &
      ,' basee sur ',tideanalysis_count,' echantillons'
!     write(6,*)'Analyse harmonique pour V ',par%rank
      do j3=1,jmax+1
      do i3=1,imax
      if(  mask_v(i3,j3,kmax+1).eq.1) then !££££££££££££££££££££££££>

            kinconnue=0
            do ktide=1,kmaxtidep1
            do k1=1,2
               kinconnue=kinconnue+1

                       tab4(kinconnue)=analysetide_v(i3,j3,kinconnue) 

            enddo
            enddo
            call leastsquaresolver(2,0,0)

                           do ktide=1,kmaxtidep1
                             k1=(ktide-1)*2+1
                             k2=(ktide-1)*2+2
                             veltidecosout_v(i3,j3,ktide)=tab7(k1)
                             veltidesinout_v(i3,j3,ktide)=tab7(k2)
                           enddo

      endif                            !££££££££££££££££££££££££>

      enddo
      enddo

! Analyse du courant 3d et de la densite potentielle:
      if(flag_tide3d_analysis==1) then !m°v°m> 03-04-20

      if(par%rank==0)write(6,'(a,a,i0,a)')'Analyse harmonique pour RHP' &
      ,' basee sur ',tideanalysis_count,' echantillons'
      do k3=1,kmax ; do j3=0,jmax+1 ; do i3=0,imax+1
      if(  mask_t(i3,j3,kmax+1).eq.1) then !££££££££££££££££££££££££>
            kinconnue=0
            do ktide=1,kmaxtidep1
            do k1=1,2
               kinconnue=kinconnue+1
                    tab4(kinconnue)=analysetide3d_t(i3,j3,k3,kinconnue)  
            enddo
            enddo
                           call leastsquaresolver(min0(ksecu,2),0,0)   !09/07/04
                           do ktide=1,kmaxtidep1
                            k1=(ktide-1)*2+1
                            k2=(ktide-1)*2+2
                            rhptidecosout_t(i3,j3,k3,ktide)=tab7(k1)
                            rhptidesinout_t(i3,j3,k3,ktide)=tab7(k2)
                           enddo

      endif                            !££££££££££££££££££££££££>
      enddo ; enddo ; enddo 

      if(par%rank==0)write(6,'(a,a,i0,a)')'Analyse harmonique vel_u' &
      ,' basee sur ',tideanalysis_count,' echantillons'
      do k3=1,kmax ; do j3=1,jmax ; do i3=1,imax+1
      if(  mask_u(i3,j3,kmax+1).eq.1) then !££££££££££££££££££££££££>
            kinconnue=0
            do ktide=1,kmaxtidep1
            do k1=1,2
               kinconnue=kinconnue+1
                    tab4(kinconnue)=analysetide3d_u(i3,j3,k3,kinconnue)  
            enddo
            enddo
            call leastsquaresolver(2,0,0)
                           do ktide=1,kmaxtidep1
                             k1=(ktide-1)*2+1
                             k2=(ktide-1)*2+2
                             vel3dtidecosout_u(i3,j3,k3,ktide)=tab7(k1)
                             vel3dtidesinout_u(i3,j3,k3,ktide)=tab7(k2)
                           enddo
      endif                            !££££££££££££££££££££££££>
      enddo ; enddo ; enddo

      if(par%rank==0)write(6,'(a,a,i0,a)')'Analyse harmonique vel_v' &
      ,' basee sur ',tideanalysis_count,' echantillons'
      do k3=1,kmax ; do j3=1,jmax+1 ; do i3=1,imax
      if(  mask_v(i3,j3,kmax+1).eq.1) then !££££££££££££££££££££££££>
            kinconnue=0
            do ktide=1,kmaxtidep1
            do k1=1,2
               kinconnue=kinconnue+1
                    tab4(kinconnue)=analysetide3d_v(i3,j3,k3,kinconnue)  
            enddo
            enddo
            call leastsquaresolver(2,0,0)
                           do ktide=1,kmaxtidep1
                             k1=(ktide-1)*2+1
                             k2=(ktide-1)*2+2
                             vel3dtidecosout_v(i3,j3,k3,ktide)=tab7(k1)
                             vel3dtidesinout_v(i3,j3,k3,ktide)=tab7(k2)
                           enddo
      endif                            !££££££££££££££££££££££££>
      enddo ; enddo ; enddo

      endif                            !m°v°m> 03-04-20

! LIBERER LA MEMOIRE DU SOLVER:
      call allocate_systeme(2,0,0,0) ! desallouer

!     write(6,*)'Avant tide_netcdf ',par%rank
! Avant de faire faire le chemin inverse de l'initialisation aux variables
! de l'analyse harmonique (afin d'être semblable aux fichiers de forcage) il
! peut être interessant d'utiliser telles quelles ces variables, pour par
! exemple enlever le signal de maree d'une solution globale....
      call tide_netcdf_create_file                                     !16-12-09

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! DECOMPOSER LE SIGNAL EN HARMONIQUES PRODUIT TRANSPOSE MATRICE EGALITE
! REALISE AU COURS DE LA SIMULATION (ANALYSETIDE_Z)
! FIN.
      return
      endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(case_.eq.2) then
      stop 'tide_analysis choix 2 supprimé!'
      endif
      if(case_.eq.3) then
      stop 'tide_analysis choix 3 supprimé!'
      endif
      if(case_.eq.0) then
      stop 'tide_analysis choix 0 supprimé!'                            !27-12-09
      endif

      end subroutine tide_analysis

!......................................................................
#ifdef bidon
      subroutine tide_analysis_ssh3d
      use module_principal ; use module_parallele
      implicit none

! Determiner la part du signal de SSH qui n'est pas liE A la marEe externe
! https://docs.google.com/document/d/1_mf-UA8AZqasnL-ifM5JWyst7YnqqU0KoVyoswwFcQc/edit

      k=1
      do j=1,jmax ; do i=1,imax
       xy_t(i,j,id_ssh3d)=dz_t(i,j,k,1)*anyv3d(i,j,k,id_prs)
      enddo       ; enddo
      do k=2,kmax
       do j=1,jmax ; do i=1,imax
        xy_t(i,j,id_ssh3d)=xy_t(i,j,id_ssh3d)+dz_t(i,j,k,1)*anyv3d(i,j,k,id_prs)
        if(par%rank==114.and.i==imax/2.and.j==jmax/2.and.mask_t(imax/2,jmax/2,kmax)==1) then
        write(66,*)k,anyv3d(i,j,k,id_prs),anyv3d(i,j,1,id_prs)
        write(67,*)k,xy_t(i,j,id_ssh3d)
        endif
       enddo       ; enddo
      enddo
      do j=1,jmax ; do i=1,imax
       xy_t(i,j,id_ssh3d)=-xy_t(i,j,id_ssh3d)/(rho*grav*hz_w(i,j,1))
       if(par%rank==114.and.i==imax/2.and.j==jmax/2.and.mask_t(imax/2,jmax/2,kmax)==1) then
       write(67,*)xy_t(i,j,id_ssh3d),anyv3d(i,j,1,id_prs)/rho/grav
       endif
      enddo       ; enddo
      if(par%rank==114)stop '??????'
      end subroutine tide_analysis_ssh3d
#endif
