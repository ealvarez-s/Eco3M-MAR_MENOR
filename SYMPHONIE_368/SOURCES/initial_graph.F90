      subroutine initial_graph
!______________________________________________________________________
! SYMPHONIE ocean model
! release 296 - last update: 12-02-21
!______________________________________________________________________
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='initial_graph'
       subroutinedescription='Reads notebook_graph'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!......................................................................
! Version Date      Description des modifications
!         27/07/03: mis en service
!         29/07/03: Appel à CALL MOYENNE_TEMPS
!         01/08/03: impression KLOW_Z, SPONGE_Z
!         02/08/03: impression VELOBC, ZTAOBC, VMEAOBC, TOBC, SOBC
!         09/08/03: impression RELAX_X RELAX_Y
!         01/10/03: impression harmonique 3D
!         02/10/03: impression harmonique 3D debug
!         10/10/03: impression streamf_r(0) et streamf_r(1)
!                  + debug i/o form='unformatted'
!         21/10/03: impression RELAX_R
!         16/03/04: impression du niveau vertical
!         12/01/05: impression ecart des variables calculées au variables
!                   forcantes
!         19/05/05: visualisation precipitations
!         03/10/05: modifs pour double precision compatible avec fichiers
!                   pvwave simple
!         13/03/06: impression densite de forcage, densite-densite de forcage
!         13/04/06: choix du repertoire des sorties graphiques
!         19/01/07: Compatible avec model_ 1DV
!         05/03/07: VARIABLESDESC construit et non plus lu dans notebook_graph
! 2010.6  02-02-10  renomme lon_t lat_t
!         11-02-10  notebook_graph = nomfichier(21)
! 2010.7  13-02-10  lecture d'un nouveau notebook_graph
! 2010.9  02-06-10  ne pas ecrire le fichier binaire
! S26     07-02-13  seul proc 0 ecrit à l'ecran
!         19-04-15  dimension varid lue dans notebook_graph
!         16-05-15  routine separee pour lire notebook_graph et vite connaitre
!                   dim_varid
!         07-09-17  modif chaine de caracteres dirgraph
! v296    12-02-21  messages dans fichiers fort.xxx
!...............................................................................
!    _________                    .__                  .__             !   (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................

      if(flag3d.eq.0)return                                               !19/01/07

      if(idate_output.eq.2)call date_output(4)                         !13-02-10

!..............................................................................
! Nom du repertoire des sorties graphiques
! Debut:
!..............................................................................
      if(dirgraph(1:1).eq.'*')then !>>>>>>>>>>>>>>
! Si ancien notebook alors par default c'est ../GRAPHIQUES/
       lname4=14
       dirgraph(1:lname4)='../graphiques/'
      else                         !>>>>>>>>>>>>>>
       do k=1,90
         if(dirgraph(k:k).eq.' ') then !-------->
          lname4=k-1
          goto 21
         endif                         !-------->
       enddo
   21  continue
       if(dirgraph(lname4:lname4).eq.'/')lname4=lname4-1
       lname4=lname4+1
       dirgraph(lname4:lname4)='/'
      endif                        !>>>>>>>>>>>>>>
      dirgraph=dirgraph(1:lname4) !07-09-17
!..............................................................................
! Nom du repertoire des sorties graphiques
! Fin.
!..............................................................................

!                            /   /   /

!..............................................................................
! Définir le nb de variables 2D (L1), variables 3D 1/2 niveaux (L2),
! variables 3D niveaux entiers (L3)
! Début:
!..............................................................................
! reset des compteurs:
      k1=0
      l1=0
      l2=0
      l3=0

! variables 2D scalaires:

      do 201 k=1,grh_nb(1)
      k1=k1+1

      if(grh_out_var(k1).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

         if(k.eq.1)then !------->
! 4 flux IR
           l1=l1+1
           l1=l1+1
           l1=l1+1
           l1=l1+1
! precipi                                                              !19/05/05
           l1=l1+1
         endif          !------->

         if(k.eq.2)then !------->
! amplitude onde maree
          do ktide=1,kmaxtide
          l1=l1+1
          enddo
         endif          !------->

         if(k.eq.3)then !------->
! phase onde maree
          do ktide=1,kmaxtide
          l1=l1+1
          enddo
         endif          !------->

         if(k.eq.4)then !------->
! amplitude onde maree
          do ktide=1,kmaxtide
          l1=l1+1
          enddo
         endif          !------->

         if(k.eq.5)then !------->
! phase onde maree
          do ktide=1,kmaxtide
          l1=l1+1
          enddo
         endif          !------->

         if(k.eq.6)then !------->
! elevation de la surface
          l1=l1+1
         endif          !------->

         if(k.eq.7)then !------->
! kmin_z
          l1=l1+1
         endif          !------->

         if(k.eq.8)then !------->
! sponge_z
          l1=l1+1
         endif          !------->

         if(k.eq.9)then !------->
! zta_obc_z
          l1=l1+1
         endif          !------->

         if(k.eq.10)then!------->
! relax_x                                                              !09/08/03
          l1=l1+1
! relax_y                                                              !09/08/03
          l1=l1+1
! relax_r                                                              !21/10/03
          l1=l1+1
         endif          !------->

         if(k.eq.11)then!------->
! streamf_r(0)                                                         !10/10/03
          l1=l1+1
! streamf_r(1)                                                         !10/10/03
          l1=l1+1
         endif          !------->

         if(k.eq.12)then!------->
! zta-ztaobc:                                                          !12/01/05
          l1=l1+1
         endif          !------->

         if(k.eq.13)then!------->
! bancs decouvrants                                                    !05/04/07
          l1=l1+1
         endif          !------->

      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

  201 continue

      l1_sca=l1

! variables 2D vecteurs:

      do 202 k=1,grh_nb(2)
      k1=k1+1

      if(grh_out_var(k1).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

! tension du vent
         if(k.eq.1)then !------->
           l1=l1+1
           l1=l1+1
         endif          !------->

! courant moyen
         if(k.eq.2)then !------->
           l1=l1+1
           l1=l1+1
         endif          !------->

! courant marée à phi=0
         if(k.eq.3)then !------->
          do ktide=1,kmaxtide
          l1=l1+1
          l1=l1+1
          enddo
         endif          !------->

! courant marée à phi=pi/2
         if(k.eq.4)then !------->
          do ktide=1,kmaxtide
          l1=l1+1
          l1=l1+1
          enddo
         endif          !------->

! courant moyen forcant
         if(k.eq.5)then !------->
           l1=l1+1
           l1=l1+1
         endif          !------->

! courant moyen - courant moyen forcant
         if(k.eq.6)then !------->                                      !12/01/05
           l1=l1+1
           l1=l1+1
         endif          !------->

      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

  202 continue

      l1_vec=l1-l1_sca

! variables 3D 1/2 niveaux scalaires:

      do 203 k=1,grh_nb(3)
      k1=k1+1

      if(grh_out_var(k1).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

         if(k.eq.1)then !------->
           l2=l2+1
         endif          !------->

         if(k.eq.2)then !------->
           l2=l2+1
         endif          !------->

         if(k.eq.3)then !------->
           l2=l2+1
         endif          !------->

         if(k.eq.4)then !------->
           do vb=1,vbmax
            l2=l2+1
           enddo
         endif          !------->

         if(k.eq.5)then !------->
! impression température forcante
           l2=l2+1
         endif          !------->

         if(k.eq.6)then !------->
! impression salinité forcante
           l2=l2+1
         endif          !------->

         if(k.eq.7)then !------->
!cimpression analyse harmonique 3D                                     !01/10/03
!c        DO K3=1,ANALYSE_MAX ! amplitudes                             !02/10/03
!c         L2=L2+1
!c        ENDDO
!c        DO K3=1,ANALYSE_MAX ! phases                                 !02/10/03
!c         L2=L2+1
!c        ENDDO
         endif          !------->

         if(k.eq.8)then !------->                                      !16/03/04
! impression niveau vertical
           l2=l2+1
         endif          !------->

         if(k.eq.9)then !------->                                      !12/01/05
! impression temperature - temperature de forcage
           l2=l2+1
         endif          !------->

         if(k.eq.10)then !------->                                     !12/01/05
! impression salinite - salinite de forcage
           l2=l2+1
         endif           !------->

         if(k.eq.11)then !------->                                     !13/03/06
! impression densite de forcage
           l2=l2+1
         endif           !------->

         if(k.eq.12)then !------->                                     !13/03/06
! impression densite - densite de forcage
           l2=l2+1
         endif           !------->

      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

  203 continue

      l2_sca=l2

! variables 3D 1/2 niveaux vecteurs:

      do 204 k=1,grh_nb(4)
      k1=k1+1

      if(grh_out_var(k1).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

! courant 3D:
         if(k.eq.1)then !------->
           l2=l2+1
           l2=l2+1
         endif          !------->

! courant 3D géostrophique:
         if(k.eq.2)then !------->
           l2=l2+1
           l2=l2+1
         endif          !------->

! courant de forcage velobc
         if(k.eq.3)then !------->
           l2=l2+1
           l2=l2+1
         endif          !------->

! courant - courant de forcage
         if(k.eq.4)then !------->                                      !12/01/05
           l2=l2+1
           l2=l2+1
         endif          !------->

      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

  204 continue

      l2_vec=l2-l2_sca

! variables 3D niveaux entiers scalaires:

      do 205 k=1,grh_nb(5)
      k1=k1+1

      if(grh_out_var(k1).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

         if(k.eq.1)then !------->
           l3=l3+1
         endif          !------->

         if(k.eq.2)then !------->
           l3=l3+1
         endif          !------->

         if(k.eq.3)then !------->
           l3=l3+1
         endif          !------->

         if(k.eq.4)then !------->
           l3=l3+1
         endif          !------->

         if(k.eq.5)then !------->
           l3=l3+1
         endif          !------->

      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

  205 continue

      l3_sca=l3

! variables 3D niveaux entiers vecteurs:
      do 206 k=1,grh_nb(6)
      k1=k1+1
      if(grh_out_var(k1).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>


      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>
  206 continue

      l3_vec=l3-l3_sca

      if(l1+l2+l3.gt.dimgrh_titre) then !pmx> !12-02-21
       write(10+par%rank,*)'l1+l2+l3 > dimgrh_titre'
       write(10+par%rank,*)'dimgrh_titre=',dimgrh_titre
       write(10+par%rank,*)'l1+l2+l3=',l1+l2+l3
       write(10+par%rank,*)'Increase dimgrh_titre in module_parameter'
       stop 'STOP initial_graph see fort.xxx files'
      endif                             !pmx>


!..............................................................................
! Définir le nb de variables 2D (L1), variables 3D 1/2 niveaux (L2),
! variables 3D niveaux entiers (L3)
! Fin.
!..............................................................................

!                            /   /   /

!..............................................................................
! Définir les titres des variables
! Début:
!..............................................................................
! reset des compteurs:
      k1=0
      k2=0

! titres variables 2D scalaires:

      do 101 k=1,grh_nb(1)
      k1=k1+1

      if(grh_out_var(k1).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

         if(k.eq.1)then !------->
           k2=k2+1
           grh_titre(k2)='flux radiatif ir'
           k2=k2+1
           grh_titre(k2)='flux radiatif solaire'
           k2=k2+1
           grh_titre(k2)='flux chaleur latente'
           k2=k2+1
           grh_titre(k2)='flux chaleur sensible'
           k2=k2+1                                                     !19/05/05
           grh_titre(k2)='precipitation m/jour'
         endif          !------->

         if(k.eq.2)then !------->
          do ktide=1,kmaxtide
          k2=k2+1
          write(grh_titre(k2),'(a,i3)')'amplitude onde',ktide
          enddo
         endif          !------->

         if(k.eq.3)then !------->
          do ktide=1,kmaxtide
          k2=k2+1
          write(grh_titre(k2),'(a,i3)')'phase onde',ktide
          enddo
         endif          !------->

         if(k.eq.4)then !------->
          do ktide=1,kmaxtide
          k2=k2+1
          write(grh_titre(k2),'(a,i3)')'amplitude potentiel',ktide
          enddo
         endif          !------->

         if(k.eq.5)then !------->
          do ktide=1,kmaxtide
          k2=k2+1
          write(grh_titre(k2),'(a,i3)')'phase potentiel',ktide
          enddo
         endif          !------->

         if(k.eq.6)then !------->
          k2=k2+1
          grh_titre(k2)='surface elevation'
         endif          !------->

         if(k.eq.7)then !------->
          k2=k2+1
          grh_titre(k2)='bottom level number'
         endif          !------->

         if(k.eq.8)then !------->
          k2=k2+1
          grh_titre(k2)='sponge layer in day'
         endif          !------->

         if(k.eq.9)then !------->
          k2=k2+1
          grh_titre(k2)='forcing elevation'
         endif          !------->

         if(k.eq.10)then!------->
          k2=k2+1
          grh_titre(k2)='relax x mpv'                                  !09/08/03
          k2=k2+1
          grh_titre(k2)='relax y mpv'                                  !09/08/03
          k2=k2+1
          grh_titre(k2)='relax r mpv'                                  !21/10/03
         endif          !------->

         if(k.eq.11)then!------->
          k2=k2+1
          grh_titre(k2)='fonction courant ebauche'                     !10/10/03
          k2=k2+1
          grh_titre(k2)='fonction courant analyse'                     !10/10/03
         endif          !------->

         if(k.eq.12)then!------->
          k2=k2+1
          grh_titre(k2)='elevation - elevation forcage'                !12/01/05
         endif          !------->

         if(k.eq.13)then!------->
          k2=k2+1
          grh_titre(k2)='bancs decouvrants'                            !05/04/07
         endif          !------->

      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

  101 continue

! titres variables 2D vecteurs:

      do 102 k=1,grh_nb(2)
      k1=k1+1

      if(grh_out_var(k1).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

         if(k.eq.1)then !------->
           k2=k2+1
           grh_titre(k2)='windstress oi'
           k2=k2+1
           grh_titre(k2)='windstress oj'
         endif          !------->

         if(k.eq.2)then !------->
           k2=k2+1
           grh_titre(k2)='courant moyen oi'
           k2=k2+1
           grh_titre(k2)='courant moyen oj'
         endif          !------->

         if(k.eq.3)then !------->
          do ktide=1,kmaxtide
          k2=k2+1
          write(grh_titre(k2),'(a,i3)')'cour oi maree phi=0 ',ktide
          k2=k2+1
          write(grh_titre(k2),'(a,i3)')'cour oj maree phi=0 ',ktide
          enddo
         endif          !------->

         if(k.eq.4)then !------->
          do ktide=1,kmaxtide
          k2=k2+1
          write(grh_titre(k2),'(a,i3)')'cour oi maree pi/2 ',ktide
          k2=k2+1
          write(grh_titre(k2),'(a,i3)')'cour oj maree pi/2 ',ktide
          enddo
         endif          !------->

         if(k.eq.5)then !------->
           k2=k2+1
           grh_titre(k2)='courant moyen forcant oi'
           k2=k2+1
           grh_titre(k2)='courant moyen forcant oj'
         endif          !------->

         if(k.eq.6)then !------->                                      !12/01/05
           k2=k2+1
           grh_titre(k2)='courant moyen - forcage oi'
           k2=k2+1
           grh_titre(k2)='courant moyen - forcage oj'
         endif          !------->

      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

  102 continue

! titres variables 3D 1/2 niveaux scalaires:

      do 103 k=1,grh_nb(3)
      k1=k1+1

      if(grh_out_var(k1).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

         if(k.eq.1)then !------->
           k2=k2+1
           grh_titre(k2)='temperature'
         endif          !------->

         if(k.eq.2)then !------->
           k2=k2+1
           grh_titre(k2)='salinite'
         endif          !------->

         if(k.eq.3)then !------->
           k2=k2+1
           grh_titre(k2)='densite'
         endif          !------->

         if(k.eq.4)then !------->
           do vb=1,vbmax
            k2=k2+1
            if(vb.le.9)then
             write(grh_titre(k2),'(a,i1)')'traceur',vb
            else
             write(grh_titre(k2),'(a,i2)')'traceur',vb
            endif
           enddo
         endif          !------->

         if(k.eq.5)then !------->
           k2=k2+1
           grh_titre(k2)='forcing temperature'
         endif          !------->

         if(k.eq.6)then !------->
           k2=k2+1
           grh_titre(k2)='forcing salinity'
         endif          !------->

         if(k.eq.7)then !------->                                      !01/10/03
!         DO K3=1,ANALYSE_MAX
!          K2=K2+1
!          WRITE(GRH_TITRE(K2),'(A,I3)')'Amplitude onde3D',K3
!         ENDDO
!         DO K3=1,ANALYSE_MAX
!          K2=K2+1
!          WRITE(GRH_TITRE(K2),'(A,I3)')'Phase onde3D',K3
!         ENDDO
         endif          !------->

         if(k.eq.8)then !------->                                      !16/03/04
           k2=k2+1
           grh_titre(k2)='vertical level of the stepped sigma grid'
         endif          !------->

         if(k.eq.9)then !------->                                      !12/01/05
           k2=k2+1
           grh_titre(k2)='temperature - forcage'
         endif          !------->

         if(k.eq.10)then !------->                                      !12/01/05
           k2=k2+1
           grh_titre(k2)='salinite - forcage'
         endif          !------->

         if(k.eq.11)then !------->                                      !13/03/06
           k2=k2+1
           grh_titre(k2)='forcing density'
         endif          !------->

         if(k.eq.12)then !------->                                      !13/03/06
           k2=k2+1
           grh_titre(k2)='densite - forcage'
         endif          !------->

      endif                        !§§§§Â§§§§§§§§§§§§§§§§§§§§§§§§§§§>

  103 continue

! titres variables 3D 1/2 niveaux vecteurs:

      do 104 k=1,grh_nb(4)
      k1=k1+1

      if(grh_out_var(k1).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

         if(k.eq.1)then !------->
           k2=k2+1
           grh_titre(k2)='courant 3d oi'
           k2=k2+1
           grh_titre(k2)='courant 3d oj'
         endif          !------->

         if(k.eq.2)then !------->
           k2=k2+1
           grh_titre(k2)='courant geostrophique oi'
           k2=k2+1
           grh_titre(k2)='courant geostrophique oj'
         endif          !------->

         if(k.eq.3)then !------->
           k2=k2+1
           grh_titre(k2)='courant forcant oi'
           k2=k2+1
           grh_titre(k2)='courant forcant oj'
         endif          !------->

         if(k.eq.4)then !------->                                      !12/01/05
           k2=k2+1
           grh_titre(k2)='courant - forcage oi'
           k2=k2+1
           grh_titre(k2)='courant - forcage oj'
         endif          !------->

      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

  104 continue

! titres variables 3D niveaux entiers scalaires:

      do 105 k=1,grh_nb(5)
      k1=k1+1

      if(grh_out_var(k1).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

         if(k.eq.1)then !------->
           k2=k2+1
           grh_titre(k2)='diffusivite verticale'
         endif          !------->

         if(k.eq.2)then !------->
           k2=k2+1
           grh_titre(k2)='energie cinetique turbulente'
         endif          !------->

         if(k.eq.3)then !------->
           k2=k2+1
           grh_titre(k2)='longueur de melange'
         endif          !------->

         if(k.eq.4)then !------->
           k2=k2+1
           grh_titre(k2)='longueur de dissipation'
         endif          !------->

         if(k.eq.5)then !------->
           k2=k2+1
           grh_titre(k2)='omega'
         endif          !------->

      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>

  105 continue

! titres variables 3D niveaux entiers vecteurs:
      do 106 k=1,grh_nb(6)
      k1=k1+1
      if(grh_out_var(k1).eq.1)then !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>


      endif                        !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§>
  106 continue

!..............................................................................
! Définir les titres des variables
! Fin.
!..............................................................................

!                            /   /   /

!..............................................................................
! Création du fichier 3Dgraph_init.out'
! Début:
!..............................................................................

! Supprimer le !02-06-10

!..............................................................................
! Création du fichier 3Dgraph_init.out'
! Fin.
!..............................................................................

      end subroutine initial_graph

!.........................................................................

      subroutine initial_graph_notebook !16-05-15
      use module_principal ; use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='initial_graph_notebook'
       subroutinedescription='Reads notebook_graph'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      do k=1,kmax_grh
      grh_out_var(k)=0
      enddo
      do k=1,dimgrh_titre
      grh_titre(k)='reset'
      enddo

      open(unit=3,file=nomfichier(21))      !lecture de notebook_graph !11-02-10
      read(3,*)

      read(3,*)texte30
      if(texte30(1:6)/='Output') then  !13-02-10
      write(*,*)'Your notebook_graph file does not correspond to this'
      write(*,*)'model version. A new notebook_graph file is available'
      write(*,*)'in the "configbox" directory'
      stop
      endif

      read(3,'(a)')dirgraph     !13/04/06
      read(3,*)dim_varid                 !19-04-15 
           if(dim_varid==0)dim_varid=100 !19-04-15
           allocate( varid(dim_varid)) ; varid=0
      read(3,*)idate_output     !13-02-10
      read(3,*)const1           !13-02-10
      read(3,*)                 !13-02-10

      graphperiod=1.d10
      if(idate_output.eq.1)graphperiod=const1*3600.*24              !15-04-11

      read(3,*)

      k=0
      do 38 loop1=1,6
      read(3,*)grh_nb(loop1)
       do 40 loop2=1,grh_nb(loop1)
       k=k+1
       if(k.gt.kmax_grh)stop 'initial_graph augmenter kmax_grh'
       read(3,*)grh_out_var(k)
   40  continue
       read(3,*)!ligne de separation

   38 continue
      close(3)

      end subroutine initial_graph_notebook

