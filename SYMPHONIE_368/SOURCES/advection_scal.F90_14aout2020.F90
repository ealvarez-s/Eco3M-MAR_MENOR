      subroutine advection_scal
!______________________________________________________________________
! SYMPHONIE ocean model
! release 287 - last update: 17-07-20
!______________________________________________________________________

      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='advection_scal'
       subroutinedescription='Driver for T & S advection subroutines'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!.......................................................................
!    _________                    .__                  .__             ! m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !
!.......................................................................
! Version date      Description des modifications
!         14/07/01: passage à sigma généralisée
!         28/07/01: passage à sigma généralisée + advection "densité"
!         23/08/01: des "return" ajoutés à la fin de chaque choix
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         27/05/02: la puissance de la fonction de stabilité n'est plus
!                   fixée à 2 mais est un bouton (CONST3) reglable.
!         29/07/02: on impose un schema upwind au point source des
!                   fleuves. Notons que ceci a entrainé une reorganisation
!                   generale des boucles.
!         14/12/04: nouveau calcul de la fonction de stabilité
!         04/01/05: derniere mouture.
!         06/01/05: derniere mouture +CL upwind
!         07/01/05: suite du point precedent.
!         31/01/05: une autre obc pour le temps t-1
!         14/06/05: modération empirique du limiteur de flux et différenciée
!                   selon la verticale et l'horizontale
!         19/06/05: suite: les coef de ponderation seront lus dans
!                   notebook_advection
!         12/11/05: advection upstream dans les canaux
!         04/05/06: amenagement pour compilation en double precision
!         22/05/06: dans les canaux advection champs pris au temp t-1
!         04/07/06: pour liberer de la memoire vive on emploie anyvar3D
!                   a la place de adhy...
!         16/10/06: Nouveau cas (ichoix=3) où la fonction de stabilite est
!                   calculee à partir du champ de densite de sorte que
!                   l'effet diffusif du schema est identique pour T et S
!                   pour une meilleure conservation des characteristiques des
!                   masses d'eau (en gros le fameux diagramme T/S).
!         20/12/06: Nouveau schema en ichoix=4
!                   Meme type de schema que ichoix=3 sauf que quand la fonction
!                   de stabilité = 1 le schema upstream est forward. En milieu
!                   estuarien la fonction de stabilité = 1
!                   ichoix=1 ou =2 stoppé!
!         09/05/07: Introduction des poids ADVMIX_X & ADVMIX_Y (maintient la
!                   consistence energetique avec l'ordre du terme correctif
!                   du gradient de pression
!         21/12/07: ajout cas 5 pour etre coherent, du point de la conservation
!                   de l'energie avec une nouvelle version du gradient de pression
!         01/02/09: dans choix 5 le schema upwind de la condition limite laterale
!                   est modifié: harmonisation avec le schema interieur.
!         10/11/08: version en accord avec les nouvelles conditions aux limites
!                   (voir egalement obc_ts.F)
!         02/04/09: Seuil min de detection estuaire passe de 20 à 25psu
! 2009.3  03-10-09  utilisation des nouveaux facteurs d'echelle verticale
!                   before remplace itimets
!         10-10-09  tem_c remplace thz_c
!         15-10-09  Le filtre temporel est calculé en même temps que l'advection
!                   application i1d
!         16-10-09  application de la nouvelle terminologie du time stepping
! 2010.3  15-01-10  subroutine advection_ts renommée advection_scal
! 2010.8  12-03-10  mixer le time filter 3L et 4L avec tfc1 et tfc2
!         17-05-10  ajout schema up3
!         19-05-10  Dans panache schema upwind
! 2010.9  31-05-10  ajouts de nouveaux schemas
!         01-06-10  gain cpu
!         10-06-10  schema d'advection vertical "3 ponderations"
!         13-06-10  appel à obc_river(3) donnant flux aux embouchures des fleuves
! 2010.11 20-07-10  option cross averaged par defaut dans vertical advection
!         11-09-10  si coordonnee ale choix du schema up3 de base
!                   sinon (c.a.d. grille sigma) choix du schema up3 modifié pour
!                   reduire le melange vertical diapycnal
! 2010.13 09-10-10  Mise à jour de la fonction de detection du panache
! 2010.14 09-11-10  upwindriver_t est un tableau qui permet d'appliquer un
!                   schema 100% upwind à proximité de l'embouchure des fleuves
!         16-11-10  limiteur incoherence hydrostatique sur calcul x6
! 2010.15 31-12-10  - Evolution de la grille sigma-step (transition continue des marches)
!                   en remplacant kmin_w par kmin_w
!                   - Extrapolation verticale simplifiee et regroupee dans advection_scal_vbc
! 2010.16 12-01-11  L'update de T et S est regroupe en une seule etape,
!                   (suppression de l'etape intermediaire apres
!                   l'advection). Le but est de simplifier le calcul
!                   du filtre temporel d'ordre eleve
! 2010.17 14-01-11  modif extension boucle calcul dans subroutine advection_scal_vbc
! 2010.18 10-02-11  amelioration calcul derivee seconde dans advection verticale up3
!         22-02-11  Retour en arriere par rapport à kmin
! 2010.25 18-01-12  Le schema d'advection est laissé au choix de l'utilisateur
! S.26    24-04-14  schema upg
!         10-05-14  upg suite
!         20-08-14  - c'est plus exact (a cause notamment de stokes) de se baser sur
!                   veldydz_u plutot que sur vel_u pour calculer certains parametres
!                   de l'advection (signe, nombre de courant, etc...)
!                   - Possibilite d'advection verticale implicite
!         21-08-14  Suite du point precedent. Un tableau special nombre de courant
!                   a ete cree. On l'utilise pour x5
!         22-08-14  - adaptation du schema upg a des mailles rectangulaires
!                   - aiguillage vers up2 dependant de wetmask
!         08-10-14  prise en compte upwindriver_t dans schema upm
!         13-12-14  adaptation du masquage a la grille sigstep
!         17-12-14  a condition de decommenter les lignes datees, possibilite de
!                   passage a upwind si [d(rhp)/dz]x*=cst > 0
!         19-12-14  Suite du point precedent: amelioration de l'indice de robustesse
!                   du schema d'advection verticale
!         22-12-14  indice robutesse adv verticale: par defaut premiere couche choix C2
!                   indice robutesse adv horizontale: suppression cas particulier kmin
!         09-04-15  evaporation/precipition = terme puit/source ssh --> il faut traiter
!                   le cas omega_w(kmax+1)/=0
!         12-04-15  suite du point precedent. Flux vertical advection de surface 
!                   100% explicite pour coherence avec flux de sel explicite dans
!                   vertmix_sal
!         09-07-15  Amelioration de la condition aux limites au fond de l'algo 
!                   de ponderation up2/c2 du schema d'advection verticale
!         10-07-15  expnum (notebook_advection) permet de deplacer les poids du
!                   schema hybride up3/up2
!         05-08-15  methode vst continue vs discontinue
!         07-11-15  la soustraction de temref etc.... est integree dans
!                   les calculs du schema upm. Le schema up3 continue d'avoir des
!                   routines de soustraction et de retablissement
!         09-11-15  ne pas passer dans certains routines si methode etatde reference non activee
!         11-11-15  ajout flag_refstate
!         14-11-15  advection_scal_rmv_ref appele apres advection_scal_diffactor
!         12-12-15  velef calculE A partir de grad(Tprime+Tref) et non pas grad(Tprime) seulement
!         05-01-16  advection pas de temps sEparEs
!         31-01-16  modif sur message ecran
!         03-02-16  velef multipliE par 2 
!         14-02-16  Ne pas calculer l'advection si modele 1DV
!         19-04-16  Fonctions de ponderation C2/UP2 & UP3/UP2 ElevEes au carrE
!         07-08-16  call couple_modes_currentnumber 
!         09-08-16  advection T,S direction Oi Oj separee + diverses mises A jour
!         15-08-16  suite du point precedent
!         19-08-16  - nbre de courant divisE par min(dz(k),dz-k-1))
!                   - coef d'hybridation adv verticale elevE au carrE
!         09-04-16  Premiere iteration: permettre grand nombre de sous_iterations advectives
!         14-11-16  Ajout schema quickest pour l'advection verticale de T et S
!         21-11-16  Subcycling schema quickest: modif sur timing du terme diffusif
!         22-11-16  Etape de finalisation quickest groupee apres les 3 directions
!         30-12-16  dti_w commentEs
!         13-01-17  anyv3d(i,j,k,id_dxdydz) contient dxdy_t(:,:)*dz_t(:,:,:,0)
!         16-01-17  Traceur en facteur de du/dx et dv/dy au temps before
!         26-01-17  Variable pivot de divergence = variable instantannee
!         28-01-17  schema quick, etat de reference
!         06-02-17  Le cumul des termes de divergence doit etre fait avant la
!                   ligne de l'advetion partielle
!         16-02-17  Pas d'advection si modele 1DV 
!         23-02-17  ofactors remplace ofactort dans advection sal
!         25-02-17  suppression if sur iadvec_ts
!         11-03-17  dti max: algo equivalent mais plus compacte pour ajustement 
!                   eventuel de wetmask
!         16-03-17  suppresion d'un appel inutile  obc_river(3..
!         19-03-17  routines de bilan exact de flux T S aux frontieres laterales ouvertes
!         08-04-17  s-z coordinate
!         25-04-17  s-z coordinate suite
!         29-04-17  Si alarme grand nombre de courant declenchee refaire le diagnostique
!                   de maniere plus detaillee pour aide au debugage
!         15-05-17  Arret d'urgence si loopmax<0
!         11-09-17  x2 remplacee par cst_
!         07-10-17  fusion deplacee avant restitution du defaut de divergence
!         17-10-17  suite du point precedent: etape annulee car prise en
!                   charge dans la matrice implicite de melange vertical
!         07-12-17  kundermin_t reduit le calcul de la moyennes aux seuls points utiles
!         05-05-18  Modif de l'aide au debugage si loopmaxts>50 
!         20-05-18  Ajout schema c2limited pour iadvec_ts==2
!         21-05-18  Ajout d'un test d'instabilite de pente
!         23-05-18  ajout lim c1 pour quickest
!         26-05-18  ajout cn_power
!         07-10-18  - diffusivite quick horizontale
!                   - T et S densite effective (schema upmx)
! v245    29-01-19 call my_outputs_obcsaltflux
! v247    25-02-19 advection verticale implicite en kmax+1
! v251    10-04-19 ajout sources sous marine
! v253    01-05-19 evitement de la division par dz
! v256    11-06-19 call my_outputs_zonesalttempflux
! v257    19-06-19 advection verticale implicite dans les zones decouvrantes
! v284    21-05-20 ajout flag_ofactor
! v287    17-07-20 utiliser tmpdirname
!..............................................................................

! Pas d'advection si modele 1DV !16-02-17
      if(flag_1dv==1) then !>>>>> !16-02-17
       tem_t(:,:,:,2)=tem_t(:,:,:,0)
       sal_t(:,:,:,2)=sal_t(:,:,:,0)
       return
      endif                !>>>>> !16-02-17

! Temperature & salinity advection diffusion
! The choice of the advection scheme is decided in notebook_advection !18-01-12
      call advection_scal_upm_sub_drv !25-02-17

      end subroutine advection_scal

!________________________________________________________________________________

      subroutine advection_scal_rmv_ref
      use module_principal
      implicit none

      if(flag_refstate==0) return !11-11-15

#ifdef bidonref
! remove reference state:
      do k=1,kmax ; do j=-1,jmax+2 ; do i=-1,imax+2
       tem_t(i,j,k,0)=tem_t(i,j,k,0)-temref_t(i,j,k)
       tem_t(i,j,k,1)=tem_t(i,j,k,1)-temref_t(i,j,k)
       sal_t(i,j,k,0)=sal_t(i,j,k,0)-salref_t(i,j,k)
       sal_t(i,j,k,1)=sal_t(i,j,k,1)-salref_t(i,j,k)
      enddo ; enddo ; enddo
#endif

      end subroutine advection_scal_rmv_ref

!________________________________________________________________________________

      subroutine advection_scal_rst_ref
      use module_principal
      implicit none

      if(flag_refstate==0) return !11-11-15

#ifdef bidonref
! restoreference state:
      do k=1,kmax ; do j=-1,jmax+2 ; do i=-1,imax+2
       tem_t(i,j,k,0)=tem_t(i,j,k,0)+temref_t(i,j,k)
       tem_t(i,j,k,1)=tem_t(i,j,k,1)+temref_t(i,j,k)
       sal_t(i,j,k,0)=sal_t(i,j,k,0)+salref_t(i,j,k)
       sal_t(i,j,k,1)=sal_t(i,j,k,1)+salref_t(i,j,k)
      enddo ; enddo ; enddo
#endif

      end subroutine advection_scal_rst_ref

!________________________________________________________________________________
      subroutine advection_scal_up3  ! 11-08-10
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='advection_scal_up3'
       subroutinedescription=                                     &
          'Computes temperature and salinity advection diffusion' &
       //' horizontal fluxes using a "3 points" advection scheme' &
       //' hybrid QUICK-UPWIND'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif
      stop 'je ne dois plus passer par advection_scal_up3'
#ifdef bidon

! Ce schema implique une C.L. de type "z2" pour depth_t
!     call obc_depth(2)

      call advection_scal_rmv_ref

      ksecu=1


      const1=0.5*dti_lp*flag3d                                            !15-10-09

      do j=1,jmax
      do i=1,imax

       anyv3d(i,j,1 ,3)=0.
       anyv3d(i,j,kmax+1,3)=0.
       anyv3d(i,j,1 ,6)=0.
       anyv3d(i,j,kmax+1,6)=0.

      enddo
      enddo

      if(cst_adv_ver>0) then !>>>>>>>>>>>>>>>>>>>>>>>>>>      !01-06-10
      call advection_scal_vbc ! Vertical Boundary Condition
      endif             !>>>>>>>>>>>>>>>>>>>>>>>>>>

      const2=cst_adv_hor/3.
      do k=1,kmax
      kp1=min0(k+1,kmax)
      km1=max0(k-1,1)

      do j=1,jmax
      do i=1,imax+1

      i1=i-( nint(sign(un,veldydz_u(i,j,k,1))) + 1 )/2
!     x5=min(max((min(sal_t(i,j,k,1),sal_t(i-1,j,k,1))-25.)/5.,zero),un)
!     x5=min(max(                                                    &
!     (sal_t(i1,j,k,0)-advsmin)*(advsmax-sal_t(i1,j,k,0))*advsstp    & !09-10-10
!                ,zero),un)

! Selection up3/up2 predefinie par la zone d'influence du fleuve
!     x5=0.5*(upwindriver_t(i,j)+upwindriver_t(i-1,j))                 !09-11-10
! Selection up3/up2 donnee par le nbre de courant:
!     x5=max(1.-abs(vel_u(i,j,k,1)*dti_lp/dx_u(i,j)),zero)
! Attention prendre veldydz plutot que vel_u a cause entre autres de stokes
!     x5=max(1.-abs(veldydz_u(i,j,k,1)*dti_lp              &
!                     /dxdy_u(i,j)                         &
!                       /dz_u(i,j,k,2)),zero) !20-08-14
! Note sur le nombre de courant. On divise par dxdy_u et dz_u mais en toute
! rigueur il faudrait prendre la plus grande de deux estimations. L'une etant
! divisee par dxdy_t et dz_t en (i,j) et l'autre en (i-1,j) puisque le flux
! du point u concerne les 2 cellules _t en (i,j) et (i-1,j). A voir comme future
! amelioration possible...
! Within dried areas (wetmask=0) the up2 scheme is 100% used
      x5=max(1.-abs(veldtodx_u(i,j,k,1)),zero) & ! 21-08-14
                    *wetmask_u(i,j)              ! 22-08-14

! Flux T:
!#ifdef bidon
      anyv3d(i,j,k,1)=-veldydz_u(i,j,k,1)*const1*(                     &

            x5*(  tem_t(i   ,j,k,1)+tem_t(i-1 ,j,k,1)                  &

       -const2*( mask_t(i1+1,j,kmax)                                   &!31-12-10
                *(tem_t(i1+1,j,k,0)                                    &
                 -tem_t(i1  ,j,k,0))                                   &

              -(  tem_t(i1  ,j,k,0)                                    &
                 -tem_t(i1-1,j,k,0)                                    &
               )*mask_t(i1-1,j,kmax) ) )                               &!31-12-10

       +(1.-x5)*2.*tem_t(i1,j,k,0)  )

! Flux S:
      anyv3d(i,j,k,4)=-veldydz_u(i,j,k,1)*const1*(                     &

            x5*(  sal_t(i   ,j,k,1)+sal_t(i-1 ,j,k,1)                  &

       -const2*( mask_t(i1+1,j,kmax)                                   &
                *(sal_t(i1+1,j,k,0)                                    &
                 -sal_t(i1  ,j,k,0))                                   &

              -(  sal_t(i1  ,j,k,0)                                    &
                 -sal_t(i1-1,j,k,0)                                    &
               )*mask_t(i1-1,j,kmax) ) )                               &

       +(1.-x5)*2.*sal_t(i1,j,k,0)  )
!#endif

#ifdef bidon
      anyv3d(i,j,k,1)=-veldydz_u(i,j,k,1)*dti_lp    &
      *( dz_t(i  ,j,k,1)*(tem_t(i  ,j,k,1))         &
        +dz_t(i-1,j,k,1)*(tem_t(i-1,j,k,1)))        &
      /( dz_t(i  ,j,k,1)                            &
        +dz_t(i-1,j,k,1)                 )            

      anyv3d(i,j,k,4)=-veldydz_u(i,j,k,1)*dti_lp    &
      *( dz_t(i  ,j,k,1)*(sal_t(i  ,j,k,1))         &
        +dz_t(i-1,j,k,1)*(sal_t(i-1,j,k,1)))        &
      /( dz_t(i  ,j,k,1)                            &
        +dz_t(i-1,j,k,1)                 )            
#endif


      enddo
      enddo

      do j=1,jmax+1
      do i=1,imax

      j1=j-( nint(sign(un,veldxdz_v(i,j,k,1))) + 1 )/2
!     x5=min(max((min(sal_t(i,j,k,1),sal_t(i,j-1,k,1))-25.)/5.,zero),un)!02/04/09
!     x5=min(max(                                                    &
!     (sal_t(i,j1,k,0)-advsmin)*(advsmax-sal_t(i,j1,k,0))*advsstp    & !09-10-10
!                ,zero),un)
! Selection up3/up2 predefinie par la zone d'influence du fleuve
!     x5=0.5*(upwindriver_t(i,j)+upwindriver_t(i,j-1))                 !09-11-10
! Selection up3/up2 donnee par le nbre de courant:
!     x5=max(1.-abs(vel_v(i,j,k,1)*dti_lp/dy_v(i,j)),zero)
! Attention prendre veldydz plutot que vel_u a cause entre autres de stokes
!     x5=max(1.-abs(veldxdz_v(i,j,k,1)*dti_lp      &
!                     /dxdy_v(i,j)                 &
!                       /dz_v(i,j,k,2)),zero) !20-08-14
! Within dried areas (wetmask=0) the up2 scheme is 100% used
      x5=max(1.-abs(veldtody_v(i,j,k,1)),zero)  & ! 21-08-14
                    *wetmask_v(i,j)               ! 22-08-14

!#ifdef bidon
! Flux T:
      anyv3d(i,j,k,2)=-veldxdz_v(i,j,k,1)*const1*(                     &

             x5*( tem_t(i,j   ,k,1)+tem_t(i,j-1 ,k,1)                  &


       -const2*( mask_t(i,j1+1,kmax)                                   &
                *(tem_t(i,j1+1,k,0)                                    &
                 -tem_t(i,j1  ,k,0))                                   &

              -(  tem_t(i,j1  ,k,0)                                    &
                 -tem_t(i,j1-1,k,0)                                    &
               )*mask_t(i,j1-1,kmax) ) )                               &

       +(1.-x5)*2.*tem_t(i,j1,k,0)  )

! Flux S:
      anyv3d(i,j,k,5)=-veldxdz_v(i,j,k,1)*const1*(                     &

             x5*( sal_t(i,j   ,k,1)+sal_t(i,j-1  ,k,1)                 &


       -const2*( mask_t(i,j1+1,kmax)                                   &
                *(sal_t(i,j1+1,k,0)                                    &
                 -sal_t(i,j1  ,k,0))                                   &

              -(  sal_t(i,j1  ,k,0)                                    &
                 -sal_t(i,j1-1,k,0)                                    &
               )*mask_t(i,j1-1,kmax) ) )                               &

       +(1.-x5)*2.*sal_t(i,j1,k,0)  )
!#endif

#ifdef bidon
      anyv3d(i,j,k,2)=-veldxdz_v(i,j,k,1)*dti_lp    &
      *( dz_t(i,j  ,k,1)*(tem_t(i,j  ,k,1))         &
        +dz_t(i,j-1,k,1)*(tem_t(i,j-1,k,1)))        &
      /( dz_t(i,j  ,k,1)                            &
        +dz_t(i,j-1,k,1)                 )            

      anyv3d(i,j,k,5)=-veldxdz_v(i,j,k,1)*dti_lp    &
      *( dz_t(i,j  ,k,1)*(sal_t(i,j  ,k,1))         &
        +dz_t(i,j-1,k,1)*(sal_t(i,j-1,k,1)))        &
      /( dz_t(i,j  ,k,1)                            &
        +dz_t(i,j-1,k,1)                 )            
#endif

      enddo
      enddo

      enddo

      call advection_scal_rst_ref
      call advection_scal_c2_ref

! VERTICAL ADVECTION
      call advection_scal_zaxis


#endif
      end subroutine advection_scal_up3
!________________________________________________________________________________

      subroutine advection_scal_c2_ref
      use module_principal
      use module_parallele
      implicit none

! 
      if(flag_refstate==0) return !11-11-15
#ifdef bidonref

       do k=1,kmax
       do j=1,jmax
       do i=1,imax+1

         anyv3d(i,j,k,1)=anyv3d(i,j,k,1)           &
                     -veldydz_u(i,j,k,1)*dti_lp    &
                      *temref_u(i,j,k)


         anyv3d(i,j,k,4)=anyv3d(i,j,k,4)           &
                     -veldydz_u(i,j,k,1)*dti_lp    &
                      *salref_u(i,j,k)

       enddo
       enddo
       enddo

       do k=1,kmax
       do j=1,jmax+1
       do i=1,imax

         anyv3d(i,j,k,2)=anyv3d(i,j,k,2)         &
                     -veldxdz_v(i,j,k,1)*dti_lp  &
                      *temref_v(i,j,k)

         anyv3d(i,j,k,5)=anyv3d(i,j,k,5)         &
                     -veldxdz_v(i,j,k,1)*dti_lp  &
                      *salref_v(i,j,k)

       enddo
       enddo
       enddo
#endif

      end subroutine advection_scal_c2_ref
   
!________________________________________________________________________________

      subroutine advection_scal_zaxis
      use module_principal
      use module_parallele
      implicit none
      double precision :: grdamp_=2. ! lower gradient amplification
#ifdef synopsis
       subroutinetitle='advection_scal_zaxis'
       subroutinedescription=                                          &
          'Temperature and Salinity vertical advection diffusion'      &
       //' fluxes using a hybrid c2-up2 scheme depending on the'       &
       //' dimensionless vertical velocity' 
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! anyv3d(:,:,:,0) is based on dimensionless velocity (or "current number")
      stop 'Ne plus passer par zaxis depuis 22-11-16'

!#ifdef bidon
      do j=1,jmax 
      do i=1,imax
! weighted c2/up2 scheme based on the dimensionless velocity 
! See also cst_adv_ver ranging [0:1] in notebook_advection:
!    set cst_adv_ver=1 to fully use the hybrid scheme. 
!    set cst_adv_ver=0 for a upwind scheme only.
! Within dried areas (wetmask=0) the up2 scheme is 100% used
!        anyv3d(i,j,k,0)=1.-upwindriver_t(i,j)+upwindriver_t(i,j)*(   & !08-10-14
!                       (1.-wetmask_t(i,j)+wetmask_t(i,j)*cst_adv_ver & !22-08-14
!                       *min(1.d0,abs( omega_w(i,j,k  ,1)*dti_lp      &
!                                   /( depth_t(i,j,k  )               &
!                                     -depth_t(i,j,k-1)))  ))  )

!        anyv3d(i,j,k,0)=1. ! 100%UP2
!        anyv3d(i,j,k,0)=0. ! 100%C2

! Cas general
      do k=kmin_w(i,j)+2,kmax                         !22-12-14

! Robustness index:
       anyv3d(i,j,k,0)=1.-(1.- &!ooooo>
! Robustness regarding convection: 
! https://docs.google.com/document/d/1WsIhfKI1GP269wglgIM6YUKXdVrONILKupKP9Gui9hI/edit
                  (1.-min(max(                                       &
              -1.+min((  rhp_t(i,j,k  )-  rhp_t(i,j,k-1))            &
                     /(depth_t(i,j,k  )-depth_t(i,j,k-1)),-small1)   &
         /min(grdamp_*(  rhp_t(i,j,k-1)-  rhp_t(i,j,k-2))            &
                     /(depth_t(i,j,k-1)-depth_t(i,j,k-2)),-small1)   &
                  ,zero),un))                                        &
! Robustness regarding rivers and dry zones:
                      *upwindriver_t(i,j)   & ! river zone influence
                      *wetmask_t(i,j)       & ! drying zone
                      *cst_adv_ver          & ! arbitrary user choice (notebook_advection)
! Robustness regarding vertical current number:
                      *max(1.-abs( omega_w(i,j,k  ,1)*dti_lp      &
!                               /( depth_t(i,j,k  )               &
!                                 -depth_t(i,j,k-1))) ,zero)      & !19-08-16
                                /min( dz_t(i,j,k  ,2)             &
                                     ,dz_t(i,j,k-1,2))) ,zero)    &  
                         )**2 !ooooo> !19-08-16


!        anyv3d(i,j,k,0)=0. ! 100%UP2
!        anyv3d(i,j,k,0)=1. ! 100%C2
       enddo

! Les cas particuliers:

! Cas particulier k=kmin_w(i,j)+1 qui ne connait pas la stratif entre
! kmin_w(i,j) et kmin_w(i,j)-1 et du coup on essaie quand meme de detecter
! la presence d'un bi-couche mais en utilisant la couche du dessus !09-07-15
      k=kmin_w(i,j)+1
! Robustness index:
      anyv3d(i,j,k,0)=1.-(1.- &!ooooo>
! Robustness regarding convection: 
                  (1.-min(max(                                       &
              -1.+min((  rhp_t(i,j,k  )-  rhp_t(i,j,k-1))            &
                     /(depth_t(i,j,k  )-depth_t(i,j,k-1)),-small1)   &
         /min(grdamp_*(  rhp_t(i,j,k+1)-  rhp_t(i,j,k  ))            &
                     /(depth_t(i,j,k+1)-depth_t(i,j,k  )),-small1)   &
                  ,zero),un))                                        &
! Robustness regarding rivers and dry zones:
                      *upwindriver_t(i,j)   & ! river zone influence
                      *wetmask_t(i,j)       & ! drying zone
                      *cst_adv_ver          & ! arbitrary user choice (notebook_advection)
! Robustness regarding vertical current number:
                      *max(1.-abs( omega_w(i,j,k  ,1)*dti_lp      &
!                               /( depth_t(i,j,k  )               &
!                                 -depth_t(i,j,k-1))) ,zero)      &
                                /min( dz_t(i,j,k  ,2)             &
                                     ,dz_t(i,j,k-1,2))) ,zero)    &  

                         )**2 !ooooo>

!     anyv3d(i,j,k,0)=0.

! Cas particuliers k=1:kmin_w(i,j) et k=kmax:
!      anyv3d(i,j,kmin_w(i,j),0)=1.
       anyv3d(i,j,1:kmin_w(i,j),0)=1.
       anyv3d(i,j,kmax+1       ,0)=1. ! A la surface schema 100% explicite pour
      enddo                           ! garantir la coherence du flux d'eau douce
      enddo                           ! qui, dans vertmix_sal, est calcule avec 
                                      ! avec sal_t(:,:,kmax,1). (Conservation bilan 
                                      ! de sel) !12-04-15

!#endif
      
! Bidouille pour renforcer le c2 (autrement dit rapprocher de 1) :
!     anyv3d(:,:,:,0)=1.-(1.-anyv3d(:,:,:,0))**2

! Explicit "c2" (centred 2 points) part of the vertical advection:
      do k=2,kmax   
      do j=1,jmax
      do i=1,imax

! Flux T:
        anyv3d(i,j,k,3)=-anyv3d(i,j,k,0)*omega_w(i,j,k,1)*dti_lp  &
! 4 points: 
!                     *( 0.5625*(tem_t(i,j,k  ,1)+tem_t(i,j,k-1,1))    &
!                       -0.0625*(tem_t(i,j,k+1,1)+tem_t(i,j,k-2,1)) )
! half half average:
!                      *0.5*(tem_t(i,j,k,1)+tem_t(i,j,k-1,1))
! Energy conserving average:
                  *( ( dz_t(i,j,k  ,1)*tem_t(i,j,k  ,1)        &
                      +dz_t(i,j,k-1,1)*tem_t(i,j,k-1,1))       &
                    /( dz_t(i,j,k  ,1) +dz_t(i,j,k-1,1)) )      
! "cross" average:
!                 *( ( dz_t(i,j,k-1,1)*tem_t(i,j,k  ,1)        &
!                     +dz_t(i,j,k  ,1)*tem_t(i,j,k-1,1))       &
!                   /( dz_t(i,j,k  ,1) +dz_t(i,j,k-1,1)) )      

! Flux S:
        anyv3d(i,j,k,6)=-anyv3d(i,j,k,0)*omega_w(i,j,k,1)*dti_lp &
! half half average:
!                      *0.5*(sal_t(i,j,k,1)+sal_t(i,j,k-1,1))
! Energy conserving average:
                  *( ( dz_t(i,j,k  ,1)*sal_t(i,j,k  ,1)        &
                      +dz_t(i,j,k-1,1)*sal_t(i,j,k-1,1))       &
                    /( dz_t(i,j,k  ,1) +dz_t(i,j,k-1,1)) )      
! "cross" average:
!                 *( ( dz_t(i,j,k-1,1)*sal_t(i,j,k  ,1)        &
!                     +dz_t(i,j,k  ,1)*sal_t(i,j,k-1,1))       &
!                   /( dz_t(i,j,k  ,1) +dz_t(i,j,k-1,1)) )     &

      enddo
      enddo
      enddo

!...............
! Conditions limites en k=1 et en k=kmax+1:
! A la surface on suppose une condition de gradient nul (ce qui signifie que !09-04-15
! la variation de T,S due aux flux de surface passera par le membre de droite
! dans vertmix.F90) sans modifier le protocole habituel. 
      k=kmax+1 !12-04-15
      do j=1,jmax
      do i=1,imax
        anyv3d(i,j,k,3)=-anyv3d(i,j,k,0)*omega_w(i,j,k,1)*dti_lp*tem_t(i,j,k-1,1)
        anyv3d(i,j,k,6)=-anyv3d(i,j,k,0)*omega_w(i,j,k,1)*dti_lp*sal_t(i,j,k-1,1)
      enddo
      enddo
! Flux nul A travers le niveau 1: !03-01-16
      do j=1,jmax
      do i=1,imax
        anyv3d(i,j,1,3)=0.
        anyv3d(i,j,1,6)=0.
      enddo
      enddo
!...............

       do k=1,kmax ; do j=1,jmax ; do i=1,imax 

         tem_t(i,j,k,2)=tem_t(i,j,k,2)+( &

                              anyv3d(i,j,k+1,3) &
                             -anyv3d(i,j,k  ,3) &
         
                                       )/dz_t(i,j,k,before)

         sal_t(i,j,k,2)=sal_t(i,j,k,2)+( &

                              anyv3d(i,j,k+1,6) &
                             -anyv3d(i,j,k  ,6) &
         
                                       )/dz_t(i,j,k,before)

       enddo ; enddo ; enddo


! Implicit "up2" (upwind 2 points) part of the vertical advection: !20-08-14
      do k=1,kmax-1 !kmax !09-04-15
      do j=1,jmax
      do i=1,imax

      tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
       (1.-anyv3d(i,j,k+1,0))*( omega_w(i,j,k+1,1)+abs(omega_w(i,j,k+1,1))) &
      +(1.-anyv3d(i,j,k  ,0))*(-omega_w(i,j,k  ,1)+abs(omega_w(i,j,k  ,1))) )

      tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
       (1.-anyv3d(i,j,k  ,0))*(-omega_w(i,j,k,1)-abs(omega_w(i,j,k,1)))

      tridia_in(i,j,k,3)=tridia_in(i,j,k,3)+0.5*dti_lp*                &
       (1.-anyv3d(i,j,k+1,0))*( omega_w(i,j,k+1,1)-abs(omega_w(i,j,k+1,1)))

      enddo
      enddo
      enddo

! A la surface on suppose une condition de gradient nul (ce qui signifie que !09-04-15
! la variation de T,S due aux flux de surface passera par le membre de droite
! dans vertmix.F90) sans modifier le protocole habituel.
! Pour obtenir les lignes suivantes, on part des lignes precedentes (cas "dans la couche")
! A1*T(k-1)+A2*T(k)+A3*T(k+1) et on remplace T(k+1) par T(k) pour avoir la condition de
! gradient nul et deduire les nouveaux coefficients: A1=A1 ; A2=A2+A3 ; A3=0. 
! Notons qu'on a force le schema a etre implicite, autrement dit anyv3d(i,j,k+1,0)=0
! et disparait du calcul. On obtient:
      k=kmax
      do j=1,jmax ; do i=1,imax

      tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
       (1.-anyv3d(i,j,k  ,0))*(-omega_w(i,j,k,1)-abs(omega_w(i,j,k,1)))

      tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
      +(1.-anyv3d(i,j,k+1,0))*2.*omega_w(i,j,k+1,1)                    &
      +(1.-anyv3d(i,j,k  ,0))*( -omega_w(i,j,k  ,1)+abs(omega_w(i,j,k  ,1))) )

      enddo ; enddo

      end subroutine advection_scal_zaxis

!________________________________________________________________________________

      subroutine advection_scal_vertical    !10-06-10
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='advection_scal_vertical'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!     stop 'advection_vetical'

      const1=dti_lp*flag3d

       do k=2,kmax
       do j=1,jmax
       do i=1,imax

       k1=k-( nint(sign(un,omega_w(i,j,k,1))) + 1 )/2
!      x5=min(max((min(sal_t(i,j,k,1),sal_t(i,j,k-1,1))-25.)/5.,zero),un)
!      x5=min(max(                                                    &
!      (sal_t(i,j,k1,0)-advsmin)*(advsmax-sal_t(i,j,k1,0))*advsstp    & !09-10-10
!                 ,zero),un)
       x5=upwindriver_t(i,j)                                          !09-11-10

! Flux T:
       anyv3d(i,j,k,3)=-omega_w(i,j,k,1)*const1*(              &


! half half average:
!           x5*0.5*(   tem_t(i,j,k   ,1)+tem_t(i,j,k-1,1) )    &

! Energy conserving average:
!               x5*( ( dz_t(i,j,k  ,1)*tem_t(i,j,k  ,1)        &
!                     +dz_t(i,j,k-1,1)*tem_t(i,j,k-1,1))       &
!                   /( dz_t(i,j,k  ,1) +dz_t(i,j,k-1,1)) )     &

! "cross" average:
                x5*( ( dz_t(i,j,k-1,1)*tem_t(i,j,k  ,1)        &
                      +dz_t(i,j,k  ,1)*tem_t(i,j,k-1,1))       &
                    /( dz_t(i,j,k  ,1) +dz_t(i,j,k-1,1)) )     &

           +(1.-x5)*tem_t(i,j,k1,0) )



! Flux S:
       anyv3d(i,j,k,6)=-omega_w(i,j,k,1)*const1*(             &

! half half average:
!          x5*0.5*(   sal_t(i,j,k   ,1)+sal_t(i,j,k-1,1) )    &

! Energy conserving average:
!               x5*( ( dz_t(i,j,k  ,1)*sal_t(i,j,k  ,1)        &
!                     +dz_t(i,j,k-1,1)*sal_t(i,j,k-1,1))       &
!                   /( dz_t(i,j,k  ,1) +dz_t(i,j,k-1,1)) )     &

! "cross" average:
                x5*( ( dz_t(i,j,k-1,1)*sal_t(i,j,k  ,1)        &
                      +dz_t(i,j,k  ,1)*sal_t(i,j,k-1,1))       &
                    /( dz_t(i,j,k  ,1) +dz_t(i,j,k-1,1)) )     &

          +(1.-x5)*sal_t(i,j,k1,0) )


      enddo
      enddo
      enddo

      end subroutine advection_scal_vertical

!________________________________________________________________________________

      subroutine advection_scal_vbc
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='advection_scal_vbc'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      do j=-1,jmax+2     !14-01-11
      do i=-1,imax+2

       tem_t(i,j,kmax+1,0)=tem_t(i,j,kmax,0)
       sal_t(i,j,kmax+1,0)=sal_t(i,j,kmax,0)

!      tem_t(i,j,kmin_w(i,j)-1,0)=tem_t(i,j,kmin_w(i,j)  ,0)
!      sal_t(i,j,kmin_w(i,j)-1,0)=sal_t(i,j,kmin_w(i,j)  ,0)
       tem_t(i,j,0,0)=tem_t(i,j,1,0)
       sal_t(i,j,0,0)=sal_t(i,j,1,0)

      enddo
      enddo

      end subroutine advection_scal_vbc

!________________________________________________________________________________
!#ifdef bidon
! CETTE VERSION S'UTILISE AVEC LES ROUTINES QUI SOUSTRAIENT ET RETABLISSENT
! L'ETAT DE REFERENCE

      subroutine advection_scal_upm  ! 24-04-14
      use module_principal
      use module_parallele
      implicit none
      double precision :: velef_,flagmsk_=1.,scalefac_
      integer :: t_=1 , flagconvec_=1
#ifdef synopsis
       subroutinetitle='advection_scal_upm'
       subroutinedescription=                                     &
          'Computes temperature and salinity advection diffusion' &
       //' horizontal fluxes using a "3 points" advection scheme' &
       //' hybrid QUICK-UPWIND'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      ksecu=1

      do j=1,jmax
      do i=1,imax

       anyv3d(i,j,1     ,3)=0.
       anyv3d(i,j,kmax+1,3)=0.
       anyv3d(i,j,1     ,6)=0.
       anyv3d(i,j,kmax+1,6)=0.

      enddo
      enddo

      if(cst_adv_ver>0)call advection_scal_vbc ! Vertical Boundary Condition

      call advection_scal_diffactor ! produit anyv3d(:,:,:,0) et anyv3d(:,:,:,7)
      call advection_scal_rmv_ref   !14-11-15

      const1=0.5*dti_lp*flag3d                                            !15-10-09
      const2=cst_adv_hor/3.
      x2=0.5*cst_adv_hor/3.

      do j=1,jmax
      do i=1,imax+1

       flagmsk_=wetmask_u(i,j)*0.5*(upwindriver_t(i,j)+upwindriver_t(i-1,j))

       k=kmin_u(i,j) ; flagconvec_=1
!      if( (  rhp_t(i,j,k)-  rhp_t(i-1,j,k))            & ! commente le 22-12-14
!         *(depth_t(i,j,k)-depth_t(i-1,j,k))>0.)flagconvec_=0
! flagconvec_=0 identifie, pour la couche du fond, les zones convectives 
! ou on applique un schema 100% upwind
! Note 1: seul le signe nous interesse, autrement dit remplacer la division par
! Dz par une multiplication ne change pas le test et evite les divisions par zero...
! Note 2: au premier coup flagconvec_ est eventuellement nul. L'operation
! flagconvec_=1 executee a la fin du calcul interieur a la boucle k supprime
! le flag convectif pour k>kmin

       scalefac_=0.5/dx_u(i,j)*dy_u(i,j)

      do k=1,kmax

! Flux T:
!     velef_=0.5*(anyv3d(i,j,k,0)+anyv3d(i,j+1,k,0))*(tem_t(i,j,k,1)-tem_t(i-1,j,k,1))/dx_u(i,j) & !22-08-14
!               *dy_u(i,j)*dz_u(i,j,k,1)*mask_u(i,j,kmax)
      velef_=scalefac_*(anyv3d(i,j,k,0)+anyv3d(i,j+1,k,0))    &
                      *( tem_t(i,j,k,1)- tem_t(i-1,j,k,1)     & 
#ifdef bidonref
                     +temref_t(i,j,k)-temref_t(i-1,j,k)       &
#endif
                                                         )    & !12-12-15
            *dz_u(i,j,k,1)                                    &
          *mask_u(i,j,kmax)

! Within dried areas (wetmask=0) the up2 scheme is 100% used
      x5=max(1.-abs(velef_*dti_lp                      &
                      /dxdy_u(i,j)                     &
                        /dz_u(i,j,k,2)),zero)          &
                                             **expnum  & !10-07-15
        *flagmsk_*flagconvec_                             !17-12-14
! Bidouille upw
!     x5=0.

      anyv3d(i,j,k,1)=    const1*(                                 &


! Schema QUICK en dehors influence des fleuves (selection par x5)
      x5*(                                                             &
! Partie advective QUICK:
!     -veldydz_u(i,j,k,1)*(tem_t(i,j,k,t_)+tem_t(i-1 ,j,k,t_))         &
      -veldydz_u(i,j,k,1)*                                             &
          (tem_t(i,j,k,t_)+tem_t(i-1 ,j,k,t_))                         &!    centre simple
!     2.*(dz_t(i,j,k,1)*tem_t(i,j,k,t_)+dz_t(i-1,j,k,1)*tem_t(i-1,j,k,t_)) &!centre nrj conserv.
!       /(dz_t(i,j,k,1)                +dz_t(i-1,j,k,1))                   &
! Partie diffusive QUICK:
       +x2*(                                                           &
             (velef_+abs(velef_))*((tem_t(i  ,j,k,0)-tem_t(i-1,j,k,0)) &
                                  -(tem_t(i-1,j,k,0)-tem_t(i-2,j,k,0)) &
                                                   *mask_t(i-2,j,kmax))&

            +(velef_-abs(velef_))*((tem_t(i+1,j,k,0)-tem_t(i,j,k,0))   &
                                  *mask_t(i+1,j,kmax)                  &
                                  -(tem_t(i,j,k,0)-tem_t(i-1,j,k,0)))  &

                                                                )  & !()x2
                                                                )  & !()x5
! Schema UPWIND dans les zones sous influence des fleuves (selection par x5)
      +(1.-x5)*(                                                   &
           -veldydz_u(i,j,k,1)* (tem_t(i,j,k,0)+tem_t(i-1,j,k,0))  &
       +abs(veldydz_u(i,j,k,1))*(tem_t(i,j,k,0)-tem_t(i-1,j,k,0))) &

                                )   !()const1

! Flux S:
!     velef_=0.5*(anyv3d(i,j,k,7)+anyv3d(i,j+1,k,7))*(sal_t(i,j,k,1)-sal_t(i-1,j,k,1))/dx_u(i,j) &
!               *dy_u(i,j)*dz_u(i,j,k,1)*mask_u(i,j,kmax)
      velef_=scalefac_*(anyv3d(i,j,k,7)+anyv3d(i,j+1,k,7))    &
                      *( sal_t(i,j,k,1)- sal_t(i-1,j,k,1)     & 
#ifdef bidonref
                     +salref_t(i,j,k)-salref_t(i-1,j,k)       
#endif
                                                         )    & !12-12-15
            *dz_u(i,j,k,1)                                    &
          *mask_u(i,j,kmax)

! Within dried areas (wetmask=0) the up2 scheme is 100% used
      x5=max(1.-abs(velef_*dti_lp                       &
                      /dxdy_u(i,j)                      &
                        /dz_u(i,j,k,2)),zero)           & !21-08-14
                                             **expnum   & !10-07-15
        *flagmsk_*flagconvec_                             !17-12-14
! Bidouille upw
!     x5=0.

      anyv3d(i,j,k,4)=   const1*(                                   &

! Schema QUICK en dehors influence des fleuves (selection par x5)
      x5*(                                                             &
! Partie advective QUICK:
!      -veldydz_u(i,j,k,1)*(sal_t(i,j,k,t_)+sal_t(i-1 ,j,k,t_))         &
      -veldydz_u(i,j,k,1)*                                             &
          (sal_t(i,j,k,t_)+sal_t(i-1,j,k,t_))                          &!    centre simple
!     2.*(dz_t(i,j,k,1)*sal_t(i,j,k,t_)+dz_t(i-1,j,k,1)*sal_t(i-1,j,k,t_)) &!centre nrj conserv.
!       /(dz_t(i,j,k,1)                +dz_t(i-1,j,k,1))                   &

! Partie diffusive QUICK:
       +x2*(                                                           &
             (velef_+abs(velef_))*((sal_t(i  ,j,k,0)-sal_t(i-1,j,k,0)) &
                                  -(sal_t(i-1,j,k,0)-sal_t(i-2,j,k,0)) &
                                                   *mask_t(i-2,j,kmax))&

            +(velef_-abs(velef_))*((sal_t(i+1,j,k,0)-sal_t(i,j,k,0))   &
                                  *mask_t(i+1,j,kmax)                  &
                                  -(sal_t(i,j,k,0)-sal_t(i-1,j,k,0)))  &

                                                                )  & !()x2
                                                                )  & !()x5
! Schema UPWIND dans les zones sous influence des fleuves (selection par x5)
      +(1.-x5)*(                                                   &
           -veldydz_u(i,j,k,1)* (sal_t(i,j,k,0)+sal_t(i-1,j,k,0))  &
       +abs(veldydz_u(i,j,k,1))*(sal_t(i,j,k,0)-sal_t(i-1,j,k,0))) &

                                )   !()const1

      flagconvec_=1
      enddo !k
      enddo !i
      enddo !j


      const1=0.5*dti_lp*flag3d                                            !15-10-09
      const2=cst_adv_hor/3.
      x2=0.5*cst_adv_hor/3.
      do j=1,jmax+1
      do i=1,imax

       flagmsk_=wetmask_v(i,j)*0.5*(upwindriver_t(i,j)+upwindriver_t(i,j-1))

       k=kmin_v(i,j) ; flagconvec_=1
!      if( (  rhp_t(i,j,k)-  rhp_t(i,j-1,k))            & ! commente le 22-12-14
!         *(depth_t(i,j,k)-depth_t(i,j-1,k))>0.)flagconvec_=0
! Note: au premier coup flagconvec_ est eventuellement nul. L'operation
! flagconvec_=1 executee a la fin du calcul interieur a la boucle k supprime
! le masquage convectif pour k>kmin

      scalefac_=0.5/dy_v(i,j)*dx_v(i,j)

      do k=1,kmax

! Flux T:
!     velef_=0.5*(anyv3d(i,j,k,0)+anyv3d(i+1,j,k,0))*(tem_t(i,j,k,1)-tem_t(i,j-1,k,1))/dy_v(i,j) &
!               *dx_v(i,j)*dz_v(i,j,k,1)*mask_v(i,j,kmax)
      velef_=scalefac_*(anyv3d(i,j,k,0)+anyv3d(i+1,j,k,0))   &
                      *( tem_t(i,j,k,1)- tem_t(i,j-1,k,1)    &
#ifdef bidonref
                     +temref_t(i,j,k)-temref_t(i,j-1,k)      &
#endif
                                                         )   &
                *dz_v(i,j,k,1)                               &
              *mask_v(i,j,kmax)

! Within dried areas (wetmask=0) the up2 scheme is 100% used
      x5=max(1.-abs(velef_*dti_lp                      &
                      /dxdy_v(i,j)                     &
                        /dz_v(i,j,k,2)),zero)          & !21-08-14
                                             **expnum  & !10-07-15
        *flagmsk_*flagconvec_                            !17-12-14
! Bidouille upw
!     x5=0.

      anyv3d(i,j,k,2)=   const1*(                                   &

! Schema QUICK en dehors influence des fleuves
      x5*(                                                             &
! Partie advective QUICK:
!     -veldxdz_v(i,j,k,1)*(tem_t(i,j,k,t_)+tem_t(i,j-1,k,t_))          &
      -veldxdz_v(i,j,k,1)*                                             &
          (tem_t(i,j,k,t_)+tem_t(i,j-1,k,t_))                          &!    centre simple
!     2.*(dz_t(i,j,k,1)*tem_t(i,j,k,t_)+dz_t(i,j-1,k,1)*tem_t(i,j-1,k,t_)) &!centre nrj conserv.
!       /(dz_t(i,j,k,1)                +dz_t(i,j-1,k,1))                   &


! Partie diffusive QUICK:
       +x2*(                                                           &
             (velef_+abs(velef_))*((tem_t(i,j  ,k,0)-tem_t(i,j-1,k,0)) &
                                  -(tem_t(i,j-1,k,0)-tem_t(i,j-2,k,0)) &
                                                   *mask_t(i,j-2,kmax))&

            +(velef_-abs(velef_))*((tem_t(i,j+1,k,0)-tem_t(i,j,k,0))   &
                                  *mask_t(i,j+1,kmax)                  &
                                  -(tem_t(i,j,k,0)-tem_t(i,j-1,k,0)))  &

                                                                   )  & !()x2
                                                                   )  & !()x5
! Schema UPWIND dans les zones sous influence des fleuves:
      +(1.-x5)*(                                                      &
          -veldxdz_v(i,j,k,1)* (tem_t(i,j,k,0)+tem_t(i,j-1,k,0))      &
      +abs(veldxdz_v(i,j,k,1))*(tem_t(i,j,k,0)-tem_t(i,j-1,k,0)))     &

                                ) !()const1

! Flux S:
!     velef_=0.5*(anyv3d(i,j,k,7)+anyv3d(i+1,j,k,7))*(sal_t(i,j,k,1)-sal_t(i,j-1,k,1))/dy_v(i,j) &
!              *dx_v(i,j)*dz_v(i,j,k,1)*mask_v(i,j,kmax)
      velef_=scalefac_*(anyv3d(i,j,k,7)+anyv3d(i+1,j,k,7))   &
                      *( sal_t(i,j,k,1)- sal_t(i,j-1,k,1)    &
#ifdef bidonref
                     +salref_t(i,j,k)-salref_t(i,j-1,k)      &
#endif
                                                         )   &
                *dz_v(i,j,k,1)                               &
              *mask_v(i,j,kmax)

! Within dried areas (wetmask=0) the up2 scheme is 100% used
      x5=max(1.-abs(velef_*dti_lp                      &
                      /dxdy_v(i,j)                     &
                        /dz_v(i,j,k,2)),zero)          & !21-08-14
                                             **expnum  & !10-07-15
        *flagmsk_*flagconvec_                             !17-12-14
! Bidouille upw
!     x5=0.

      anyv3d(i,j,k,5)=   const1*(                                   &

! Schema QUICK en dehors influence des fleuves
      x5*(                                                             &
! Partie advective QUICK:
!     -veldxdz_v(i,j,k,1)*(sal_t(i,j,k,t_)+sal_t(i,j-1,k,t_))          &
      -veldxdz_v(i,j,k,1)*                                             &
          (sal_t(i,j,k,t_)+sal_t(i,j-1,k,t_))                          &!    centre simple
!     2.*(dz_t(i,j,k,1)*sal_t(i,j,k,t_)+dz_t(i,j-1,k,1)*sal_t(i,j-1,k,t_)) &!centre nrj conserv.
!       /(dz_t(i,j,k,1)                +dz_t(i,j-1,k,1))                   &

! Partie diffusive QUICK:
       +x2*(                                                           &
             (velef_+abs(velef_))*((sal_t(i,j  ,k,0)-sal_t(i,j-1,k,0)) &
                                  -(sal_t(i,j-1,k,0)-sal_t(i,j-2,k,0)) &
                                                   *mask_t(i,j-2,kmax))&

            +(velef_-abs(velef_))*((sal_t(i,j+1,k,0)-sal_t(i,j,k,0))   &
                                  *mask_t(i,j+1,kmax)                  &
                                  -(sal_t(i,j,k,0)-sal_t(i,j-1,k,0)))  &

                                                                   )  & !()x2
                                                                   )  & !()x5
! Schema UPWIND dans les zones sous influence des fleuves:
      +(1.-x5)*(                                                      &
          -veldxdz_v(i,j,k,1)* (sal_t(i,j,k,0)+sal_t(i,j-1,k,0))      &
      +abs(veldxdz_v(i,j,k,1))*(sal_t(i,j,k,0)-sal_t(i,j-1,k,0)))     &

                                ) !()const1

!      if(i+par%timax(1)==30.and.j+par%tjmax(1)==3.and.k==kmax) then
!       write(6,*)'mask x5 veldxdz',mask_v(i,j,k),real(x5),real(veldxdz_v(i,j,k,1)),real(sal_t(i,j-1,k,0))
!       write(6,*)'mask x5 veldxdz',anyv3d(i,j,k,5)
!      endif

      flagconvec_=1
      enddo !k
      enddo !i
      enddo !j

! Methode vst discontinue
      if(vststep==1) then !>>>>>> !05-08-15
! Note: je passe par une boucle de k=1,kmin-1 plutot que par
! une boucle complete avec multiplication par mask car cette
! seconde option aurait pour effet d'annuler l'input des fleuves
        do j=1,jmax
        do i=1,imax+1
        do k=1,kmin_u(i,j)-1
         anyv3d(i,j,k,1)=0.
         anyv3d(i,j,k,4)=0.
        enddo
        enddo
        enddo
        do j=1,jmax+1
        do i=1,imax
        do k=1,kmin_v(i,j)-1
         anyv3d(i,j,k,2)=0.
         anyv3d(i,j,k,5)=0.
        enddo
        enddo
        enddo
      endif               !>>>>>>

      call advection_scal_rst_ref
      call advection_scal_c2_ref

! Vertical advection
      call advection_scal_zaxis

      end subroutine advection_scal_upm

!#endif
!________________________________________________________________________________

!________________________________________________________________________________

      subroutine advection_scal_diffactor ! 24-04-14
      use module_principal
      use module_parallele
      implicit none
      double precision velx_,vely_,grdx_,grdy_
#ifdef synopsis
       subroutinetitle='advection_scal_diffactor'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      stop 'Je ne dois plus passer par advection_scal_diffactor'
#ifdef bidon

! compute veldtodx_u and veldtody_v:
      call couple_modes_currentnumber !07-08-16
!     write(6,*)veldtodx_u(imax/2,jmax/2,kmax/2,1)
!     write(6,*)veldtody_v(imax/2,jmax/2,kmax/2,1)
!     stop 'naz'
! "+2" range OBC on veldtodx_u and veldtody_v
      call obc_int_cn_o4

      x1=cos(45.*deg2rad)
      x2=sin(45.*deg2rad)
      x3=sqrt(2.)
! Pourquoi la multiplication par x3=sqrt(2)? Parce que alpha est relatif a dT/di dT/dj
! et que le calcul suivant se base sur dT/did dT/djd ou id et jd sont les indices
! diagonaux et que la proportionalite de l'accroissement dans ces deux espaces est sqrt(2)
! Attention que le calcul suivant suppose dx=dy ce qui ne sera pas toujours le cas.
      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax+1

#ifdef bidon
! Pour T
       anyv3d(i,j,k,0)=x3*mask_f(i,j,k)*(               &

             0.5*( (vel_u(i,j,k,1)+vel_u(i,j-1,k,1))*x1     &   !(velx_
                  -(vel_v(i,j,k,1)+vel_v(i-1,j,k,1))*x2 )   &
                  *(tem_t(i,j-1,k,1)-tem_t(i-1,j,k,1))      &   ! *grdx_

            +0.5*( (vel_u(i,j,k,1)+vel_u(i,j-1,k,1))*x2     &   ! vely_
                  +(vel_v(i,j,k,1)+vel_v(i-1,j,k,1))*x1 )   &   !
                  *(tem_t(i,j,k,1)-tem_t(i-1,j-1,k,1))    ) &   ! grady_)

            / ( small1 +                                    &   ! /
                (tem_t(i,j-1,k,1)-tem_t(i-1,j,k,1))**2      &   !(gradx_**2
               +(tem_t(i,j,k,1)-tem_t(i-1,j-1,k,1))**2 )        ! grady-**2)

! Pour S
       anyv3d(i,j,k,7)=x3*mask_f(i,j,k)*(               &

             0.5*( (vel_u(i,j,k,1)+vel_u(i,j-1,k,1))*x1     &   !(velx_
                  -(vel_v(i,j,k,1)+vel_v(i-1,j,k,1))*x2 )   &
                  *(sal_t(i,j-1,k,1)-sal_t(i-1,j,k,1))      &   ! *grdx_

            +0.5*( (vel_u(i,j,k,1)+vel_u(i,j-1,k,1))*x2     &   ! vely_
                  +(vel_v(i,j,k,1)+vel_v(i-1,j,k,1))*x1 )   &   !
                  *(sal_t(i,j,k,1)-sal_t(i-1,j-1,k,1))    ) &   ! grady_)

            / ( small1 +                                    &   ! /
                (sal_t(i,j-1,k,1)-sal_t(i-1,j,k,1))**2      &   !(gradx_**2
               +(sal_t(i,j,k,1)-sal_t(i-1,j-1,k,1))**2 )        ! grady-**2)

#endif

! decommenter lignes suivantes pour verifier la coherence du masque
!      if(mask_f(i,j,k)==1) then
!      if(mask_t(i,j,k)  *mask_t(i-1,j,k)       &
!        *mask_t(i,j-1,k)*mask_t(i-1,j-1,k)==0) &
!        stop ' masques incoherents 2316'
!      endif

#ifdef bidon
! Pour T
       anyv3d(i,j,k,0)=                                            &

        ( vel_u(i,j  ,k,1)*(tem_t(i,j  ,k,1)-tem_t(i-1,j  ,k,1))   &
         +vel_u(i,j-1,k,1)*(tem_t(i,j-1,k,1)-tem_t(i-1,j-1,k,1))   & ! 2*U*dT/dx

         +vel_v(i  ,j,k,1)*(tem_t(i  ,j,k,1)-tem_t(i  ,j-1,k,1))   &
         +vel_v(i-1,j,k,1)*(tem_t(i-1,j,k,1)-tem_t(i-1,j-1,k,1)) ) & ! 2*V*dT/dy

       /( mask_u(i,j  ,k)*(tem_t(i,j  ,k,1)-tem_t(i-1,j  ,k,1))**2 & ! 2*dT/dx**2
         +mask_u(i,j-1,k)*(tem_t(i,j-1,k,1)-tem_t(i-1,j-1,k,1))**2 &

         +mask_v(i  ,j,k)*(tem_t(i  ,j,k,1)-tem_t(i  ,j-1,k,1))**2 & ! 2*dT/dy**2
         +mask_v(i-1,j,k)*(tem_t(i-1,j,k,1)-tem_t(i-1,j-1,k,1))**2+small1 )

! Pour S
       anyv3d(i,j,k,7)=                                            &

        ( vel_u(i,j  ,k,1)*(sal_t(i,j  ,k,1)-sal_t(i-1,j  ,k,1))   &
         +vel_u(i,j-1,k,1)*(sal_t(i,j-1,k,1)-sal_t(i-1,j-1,k,1))   & ! 2*U*dT/dx

         +vel_v(i  ,j,k,1)*(sal_t(i  ,j,k,1)-sal_t(i  ,j-1,k,1))   &
         +vel_v(i-1,j,k,1)*(sal_t(i-1,j,k,1)-sal_t(i-1,j-1,k,1)) ) & ! 2*V*dT/dy

       /( mask_u(i,j  ,k)*(sal_t(i,j  ,k,1)-sal_t(i-1,j  ,k,1))**2 & ! 2*dT/dx**2
         +mask_u(i,j-1,k)*(sal_t(i,j-1,k,1)-sal_t(i-1,j-1,k,1))**2 &

         +mask_v(i  ,j,k)*(sal_t(i  ,j,k,1)-sal_t(i  ,j-1,k,1))**2 & ! 2*dT/dy**2
         +mask_v(i-1,j,k)*(sal_t(i-1,j,k,1)-sal_t(i-1,j-1,k,1))**2+small1 )

       x10=anyv3d(i,j,k,0)
       x11=anyv3d(i,j,k,7)

#endif

!#ifdef bidon
!NOTES: id_tcn=0 et id_scn=7
! J'ai verifie que cet algo conduisait A velef=veldydz_u si bruit en Oi sur U ou T,S
! J'ai verifie que cet algo conduisait A velef=veldxdz_v si bruit en Oj sur V ou T,S
! MAIS J'AI aussi constate que la moitie de veldxdz si jamais il y a du bruit A la fois
! sur Oi et Oj et donc cela conduit plus tard A considerer le double de la vitesse
! reconstituee...
! Pour T
       anyv3d(i,j,k,id_tcn)=(1./dti_lp)*                               &

       ( veldtodx_u(i,j  ,k,1)*(tem_t(i,j  ,k,1)-tem_t(i-1,j  ,k,1))   &
        +veldtodx_u(i,j-1,k,1)*(tem_t(i,j-1,k,1)-tem_t(i-1,j-1,k,1))   &! 2*U*dT/dx

        +veldtody_v(i  ,j,k,1)*(tem_t(i  ,j,k,1)-tem_t(i  ,j-1,k,1))   &
        +veldtody_v(i-1,j,k,1)*(tem_t(i-1,j,k,1)-tem_t(i-1,j-1,k,1)) ) &! 2*V*dT/dy

       /( (mask_u(i,j  ,k)*(tem_t(i,j  ,k,1)-tem_t(i-1,j  ,k,1))/dx_u(i,j  ))**2 & !22-06-14
         +(mask_u(i,j-1,k)*(tem_t(i,j-1,k,1)-tem_t(i-1,j-1,k,1))/dx_u(i,j-1))**2 &
         +(mask_v(i  ,j,k)*(tem_t(i  ,j,k,1)-tem_t(i  ,j-1,k,1))/dy_v(i  ,j))**2 &
         +(mask_v(i-1,j,k)*(tem_t(i-1,j,k,1)-tem_t(i-1,j-1,k,1))/dy_v(i-1,j))**2 &
                        +small2 )

! Pour S
       anyv3d(i,j,k,id_scn)=(1./dti_lp)*                               &

       ( veldtodx_u(i,j  ,k,1)*(sal_t(i,j  ,k,1)-sal_t(i-1,j  ,k,1))   &
        +veldtodx_u(i,j-1,k,1)*(sal_t(i,j-1,k,1)-sal_t(i-1,j-1,k,1))   &! 2*U*dT/dx

        +veldtody_v(i  ,j,k,1)*(sal_t(i  ,j,k,1)-sal_t(i  ,j-1,k,1))   &
        +veldtody_v(i-1,j,k,1)*(sal_t(i-1,j,k,1)-sal_t(i-1,j-1,k,1)) ) &! 2*V*dT/dy

       /( (mask_u(i,j  ,k)*(sal_t(i,j  ,k,1)-sal_t(i-1,j  ,k,1))/dx_u(i,j  ))**2 &
         +(mask_u(i,j-1,k)*(sal_t(i,j-1,k,1)-sal_t(i-1,j-1,k,1))/dx_u(i,j-1))**2 &
         +(mask_v(i  ,j,k)*(sal_t(i  ,j,k,1)-sal_t(i  ,j-1,k,1))/dy_v(i  ,j))**2 &
         +(mask_v(i-1,j,k)*(sal_t(i-1,j,k,1)-sal_t(i-1,j-1,k,1))/dy_v(i-1,j))**2 &
                        +small2 )


!        if(anyv3d(i,j,k,id_scn)/=anyv3d(i,j,k,id_tcn)) &
!        stop 'anyv3d(i,j,k,id_scn)/=anyv3d(i,j,k,id_tcn)'
!#endif

#ifdef bidon
! Cette discretisation cherche A annuler la contribution d'un champ de vitesse
! buitE A la construction de velef car la composante bruitEe intervient autrement
! (A 100%) dans la diffusivitEe
! Temperature Current Number:
      anyv3d(i,j,k,id_tcn)=(0.5/dti_lp)* & ! Le 0.5 est pour contrer 2* devant dT/dx et dt/dy

       ( & !pmxpmx>

        (veldtodx_u(i,j  ,k,1)+veldtodx_u(i,j-1,k,1))*( & !pjp1>

         mask_u(i,j  ,kmax)*(tem_t(i,j  ,k,1)-tem_t(i-1,j  ,k,1))   &
        +mask_u(i,j-1,kmax)*(tem_t(i,j-1,k,1)-tem_t(i-1,j-1,k,1))   &! 2*U*2*dT/dx

                                                      ) & !pjp1>

       +(veldtody_v(i  ,j,k,1)+veldtody_v(i-1,j,k,1))*( & !pjp2>

        +mask_v(i  ,j,kmax)*(tem_t(i  ,j,k,1)-tem_t(i  ,j-1,k,1))   &
        +mask_v(i-1,j,kmax)*(tem_t(i-1,j,k,1)-tem_t(i-1,j-1,k,1))   &! 2*V*2*dT/dy

                                                      ) & !pjp2>

       ) & !pmxpmx>

       /( (mask_u(i,j  ,kmax)*(tem_t(i,j  ,k,1)-tem_t(i-1,j  ,k,1))/dx_u(i,j  ))**2 & 
         +(mask_u(i,j-1,kmax)*(tem_t(i,j-1,k,1)-tem_t(i-1,j-1,k,1))/dx_u(i,j-1))**2 &
         +(mask_v(i  ,j,kmax)*(tem_t(i  ,j,k,1)-tem_t(i  ,j-1,k,1))/dy_v(i  ,j))**2 &
         +(mask_v(i-1,j,kmax)*(tem_t(i-1,j,k,1)-tem_t(i-1,j-1,k,1))/dy_v(i-1,j))**2 &
                        +small2 )

! Pour S
! Salinity Current Number:
       anyv3d(i,j,k,id_scn)=(0.5/dti_lp)* & ! Le 0.5 est pour contrer 2* devant dS/dx et dt/dy

       ( & !pmxpmx>

        (veldtodx_u(i,j  ,k,1)+veldtodx_u(i,j-1,k,1))*( & !pjp1>

         mask_u(i,j  ,kmax)*(sal_t(i,j  ,k,1)-sal_t(i-1,j  ,k,1))   &
        +mask_u(i,j-1,kmax)*(sal_t(i,j-1,k,1)-sal_t(i-1,j-1,k,1))   &! 2*U*2*dS/dx

                                                      ) & !pjp1>

       +(veldtody_v(i  ,j,k,1)+veldtody_v(i-1,j,k,1))*( & !pjp2>

        +mask_v(i  ,j,kmax)*(sal_t(i  ,j,k,1)-sal_t(i  ,j-1,k,1))   &
        +mask_v(i-1,j,kmax)*(sal_t(i-1,j,k,1)-sal_t(i-1,j-1,k,1))   &! 2*V*2*dS/dy

                                                      ) & !pjp2>

       ) & !pmxpmx>

       /( (mask_u(i,j  ,kmax)*(sal_t(i,j  ,k,1)-sal_t(i-1,j  ,k,1))/dx_u(i,j  ))**2 & 
         +(mask_u(i,j-1,kmax)*(sal_t(i,j-1,k,1)-sal_t(i-1,j-1,k,1))/dx_u(i,j-1))**2 &
         +(mask_v(i  ,j,kmax)*(sal_t(i  ,j,k,1)-sal_t(i  ,j-1,k,1))/dy_v(i  ,j))**2 &
         +(mask_v(i-1,j,kmax)*(sal_t(i-1,j,k,1)-sal_t(i-1,j-1,k,1))/dy_v(i-1,j))**2 &
                        +small2 )
#endif

!      if(iteration3d==10) then
!       write(10+par%rank,*)x10,anyv3d(i,j,k,0)
!       write(20+par%rank,*)x11,anyv3d(i,j,k,7)
!      endif

      enddo
      enddo
      enddo
!      if(iteration3d==10) stop 'yep'
#endif
      end subroutine advection_scal_diffactor

!................................................................

      subroutine advection_scal_ofactor
      use module_principal
      use module_parallele
      implicit none

! Dans le cas ou le critere d'orthogonalite du courant n'est pas applique (flag_ofactor==0)
! alors fixer les tableaux anyv3d(:,:,:,id_ofactort) anyv3d(:,:,:,id_ofactors) a 1 et sortir:
      if(flag_ofactor==0) then !m°v°m> !21-05-20
         anyv3d(:,:,:,id_ofactort)=1.
         anyv3d(:,:,:,id_ofactors)=1.
       RETURN
      endif                    !m°v°m>

! Commentaires dans:
!https://docs.google.com/document/d/1-4NdFRKCug7u4-reURBVsHOeHMtf5p_8RlCYDLwAZjA/edit

      do k=1,kmax
       do j=2,jmax
       do i=2,imax
! u:
        xy_f(i,j,1)=0.5*( veldydz_u(i,j  ,k,1) &
                              /dz_u(i,j  ,k,1) &
                              /dy_u(i,j  )     &
                         +veldydz_u(i,j-1,k,1) &
                              /dz_u(i,j-1,k,1) &
                              /dy_u(i,j-1)     )
! v:
        xy_f(i,j,2)=0.5*( veldxdz_v(i,j  ,k,1) &
                              /dz_v(i,j  ,k,1) &
                              /dx_v(i,j  )     &
                         +veldxdz_v(i,j-1,k,1) &
                              /dz_v(i,j-1,k,1) &
                              /dx_v(i,j-1)     )
! dT/dx:
        xy_f(i,j,3)=0.5*(mask_u(i,j  ,kmax)*(tem_t(i,j  ,k,1)-tem_t(i-1,j  ,k,1))/dx_u(i,j  ) &
                        +mask_u(i,j-1,kmax)*(tem_t(i,j-1,k,1)-tem_t(i-1,j-1,k,1))/dx_u(i,j-1))

! dT/dy:
        xy_f(i,j,4)=0.5*(mask_v(i  ,j,kmax)*(tem_t(i  ,j,k,1)-tem_t(i  ,j-1,k,1))/dy_v(i  ,j) &
                        +mask_v(i-1,j,kmax)*(tem_t(i-1,j,k,1)-tem_t(i-1,j-1,k,1))/dy_v(i-1,j))

! dS/dx:
        xy_f(i,j,5)=0.5*(mask_u(i,j  ,kmax)*(sal_t(i,j  ,k,1)-sal_t(i-1,j  ,k,1))/dx_u(i,j  ) &
                        +mask_u(i,j-1,kmax)*(sal_t(i,j-1,k,1)-sal_t(i-1,j-1,k,1))/dx_u(i,j-1))

! dS/dy:
        xy_f(i,j,6)=0.5*(mask_v(i  ,j,kmax)*(sal_t(i  ,j,k,1)-sal_t(i  ,j-1,k,1))/dy_v(i  ,j) &
                        +mask_v(i-1,j,kmax)*(sal_t(i-1,j,k,1)-sal_t(i-1,j-1,k,1))/dy_v(i-1,j))

       enddo
       enddo
       do j=2,jmax
       do i=2,imax
         anyv3d(i,j,k,id_ofactort)=                             &
         abs( xy_f(i,j,1)*xy_f(i,j,3)+xy_f(i,j,2)*xy_f(i,j,4))  & !   abs( u.dT/dx+v.dT/dy)
       /(abs( xy_f(i,j,1)*xy_f(i,j,3)+xy_f(i,j,2)*xy_f(i,j,4))  & ! /(abs( u.dT/dx+v.dT/dy)
        +abs(-xy_f(i,j,2)*xy_f(i,j,3)+xy_f(i,j,1)*xy_f(i,j,4))  & !  +abs(-v.dT/dx+u.dT/dy))
        +small2)
         anyv3d(i,j,k,id_ofactors)=                             &
         abs( xy_f(i,j,1)*xy_f(i,j,5)+xy_f(i,j,2)*xy_f(i,j,6))  & !   abs( u.dS/dx+v.dS/dy)
       /(abs( xy_f(i,j,1)*xy_f(i,j,5)+xy_f(i,j,2)*xy_f(i,j,6))  & ! /(abs( u.dS/dx+v.dS/dy)
        +abs(-xy_f(i,j,2)*xy_f(i,j,5)+xy_f(i,j,1)*xy_f(i,j,6))  & !  +abs(-v.dS/dx+u.dS/dy))
        +small2)
       enddo
       enddo
      enddo ! k loop

!.....................
! Border conditions:
      if(obcstatus(ieq1)==1)   anyv3d(1     ,:     ,:,id_ofactort)=anyv3d(2   ,:   ,:,id_ofactort)
      if(obcstatus(jeq1)==1)   anyv3d(:     ,1     ,:,id_ofactort)=anyv3d(:   ,2   ,:,id_ofactort)
      if(obcstatus(ieqimax)==1)anyv3d(imax+1,:     ,:,id_ofactort)=anyv3d(imax,:   ,:,id_ofactort)
      if(obcstatus(jeqjmax)==1)anyv3d(:     ,jmax+1,:,id_ofactort)=anyv3d(:   ,jmax,:,id_ofactort)

      if(obcstatus(ieq1)==1)   anyv3d(1     ,:     ,:,id_ofactors)=anyv3d(2   ,:   ,:,id_ofactors)
      if(obcstatus(jeq1)==1)   anyv3d(:     ,1     ,:,id_ofactors)=anyv3d(:   ,2   ,:,id_ofactors)
      if(obcstatus(ieqimax)==1)anyv3d(imax+1,:     ,:,id_ofactors)=anyv3d(imax,:   ,:,id_ofactors)
      if(obcstatus(jeqjmax)==1)anyv3d(:     ,jmax+1,:,id_ofactors)=anyv3d(:   ,jmax,:,id_ofactors)

      call obc_int_anyv3d(id_ofactort,'r')
      call obc_int_anyv3d(id_ofactors,'r')
!.....................

!     do k=1,kmax
!     do j=0,jmax+1
!     do i=0,imax+1
!      i1=max(i,1) ; j1=max(j,1)
!      tem_t(i,j,k,1)=anyv3d(i1,j1,k,2)
!     enddo
!     enddo
!     enddo
!     call graph_out
!     stop 'hihan'

      end subroutine advection_scal_ofactor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                       END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine advection_scal_upcfl
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='advection_scal_upcfl'
       subroutinedescription= &
          'Computes the stability condition of the upwind advection' &
       //' scheme'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      do  k=1,kmax    
      do  j=1,jmax
      do  i=1,imax

      x1=sign(un,veldydz_u(i+1,j,k,1))
      x2=sign(un,veldydz_u(i,  j,k,1))
      x3=sign(un,veldxdz_v(i,j+1,k,1))
      x4=sign(un,veldxdz_v(i,j  ,k,1))
      x5=sign(un,omega_w(i,j,k+1,1))
      x6=sign(un,omega_w(i,j,k  ,1))

! Coef en facteur de bio_z(i,j,k):
      const0=(dz_t(i,j,k,0)+                                &
       (  veldydz_u(i+1,j  ,k  ,1)*dti_fw*(-1.-x1)                   &
         +veldydz_u(i  ,j  ,k  ,1)*dti_fw*( 1.-x2)                   &
         +veldxdz_v(i  ,j+1,k  ,1)*dti_fw*(-1.-x3)                   &
         +veldxdz_v(i  ,j  ,k  ,1)*dti_fw*( 1.-x4)  )/dxdy_t(i,j)    &
           +omega_w(i  ,j  ,k+1,1)*dti_fw*(-1.-x5)                   &
           +omega_w(i  ,j  ,k  ,1)*dti_fw*( 1.-x6)                   &
       )/dz_t(i,j,k,2)

! Coef en facteur de bio_t(i+1,j,k):
      const1= &
      veldydz_u(i+1,j  ,k  ,1)*dti_fw*(-1.+x1)/dxdy_t(i,j)/dz_t(i,j,k,2)

! Coef en facteur de bio_t(i-1,j,k):
      const2=  &
      veldydz_u(i  ,j  ,k  ,1)*dti_fw*( 1.+x2)/dxdy_t(i,j)/dz_t(i,j,k,2)

! Coef en facteur de bio_t(i,j+1,k):
      const3=  &
      veldxdz_v(i  ,j+1,k  ,1)*dti_fw*(-1.+x3)/dxdy_t(i,j)/dz_t(i,j,k,2)

! Coef en facteur de bio_t(i,j-1,k):
      const4=  &
      veldxdz_v(i  ,j  ,k  ,1)*dti_fw*( 1.+x4)/dxdy_t(i,j)/dz_t(i,j,k,2)

! Coef en facteur de bio_t(i,j,k+1):
      const5=  &
      omega_w(i  ,j  ,k+1,1)*dti_fw*(-1.+x5)/dz_t(i,j,k,2)

! Coef en facteur de bio_t(i,j,k-1):
      const6=  &
      omega_w(i  ,j  ,k  ,1)*dti_fw*( 1.+x6)/dz_t(i,j,k,2)

      if(i+par%timax(1)==186.and.&
         j+par%tjmax(1)==432.and.k==13)write(6,'(8(1x,f10.5))') &
       const0,const1,const2,const3,const4,const5,const6         &
      ,const0+const1+const2+const3+const4+const5+const6          

      enddo
      enddo
      enddo

      end subroutine advection_scal_upcfl

!.................................................................

      subroutine advection_scal_upm_sub_drv
      use module_principal
      use module_parallele
      implicit none

       call advection_scal_dz        ! computes the "true" dz
       call advection_scal_dtisub    ! computes loopmaxts et dti_lpsub
!      call advection_scal_diffactor ! computes anyv3d(:,:,:,id_tcn) anyv3d(:,:,:,id_scn)
       call advection_scal_ofactor   ! computes anyv3d(:,:,:,id_ofactor)
       call advection_scal_hybcoef   ! computes anyv3d(:,:,:,id_hybcoefu) anyv3d(:,:,:,id_hybcoefv)

! Ici anyv3d(:,:,:,4:5) est disponible
!      if(flag_ts_effectivedensity==0) then !m°v°m> 07-10-18
       if(flag_ts_quicklim==0) then !m°v°m> 07-10-18
        call advection_scal_upm_subt  
        call advection_scal_upm_subs  
       else                         !m°v°m> 07-10-18
!       call advection_scal_upmx_subt  
!       call advection_scal_upmx_subs  
        call advection_scal_uplim_subt  
        call advection_scal_uplim_subs  
       endif                        !m°v°m> 07-10-18

       if(cst_adv_ver>0)call advection_scal_vbc ! Vertical Boundary Condition
!      call advection_scal_zaxis
!      call advection_scal_z_laxwen   !12-11-16
       if(iadvec_ts==1) then !-vertical-advection->
          call advection_scal_z_quickest              !14-11-16
       else                  !-vertical-advection->
          call advection_scal_z_c2limited             !20-05-18
       endif                 !-vertical-advection->



      end subroutine advection_scal_upm_sub_drv

!.................................................................

      subroutine advection_scal_dz !22-08-16
      use module_principal ; use module_parallele
      implicit none

!     if(flag_merged_levels==0) then !pmxpmx>

! Standard case:
!      do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
!!!!!   anyv3d(i,j,k,id_dxdydz)=dxdy_t(i,j)*dz_t(i,j,k,0)
!       anyv3d(i,j,k,id_dxdydz)=dxdy_t(i,j)*min(dz_t(i,j,k,0),dz_t(i,j,k,2))
!      enddo ; enddo ; enddo

! Ici k=kmerged_t est la plus haute des 2 cellules fusionnEes:
       do j=0,jmax+1 ; do i=0,imax+1
        do k=kmerged_t(i,j),kmax
         anyv3d(i,j,k,id_dxdydz)=dxdy_t(i,j)*min(dz_t(i,j,k,0),dz_t(i,j,k,2))
        enddo
        do k=1,kmerged_t(i,j)-1
         anyv3d(i,j,k,id_dxdydz)=anyv3d(i,j,kmerged_t(i,j),id_dxdydz)
        enddo
       enddo ; enddo

!     else                           !pmxpmx>

! Merged levels case:
!      do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
!       anyv3d(i,j,k,id_dxdydz)=dxdy_t(i,j)*dz_t(i,j,min(k,int(truekmax_t(i,j))),0)   ! tkmax+0
!       anyv3d(i,j,k,id_dxdydz)=dxdy_t(i,j)*dz_t(i,j,min(k,int(truekmax_t(i,j))+1),0) ! tkmax+1 au choix...
!      enddo ; enddo ; enddo

!     endif                          !pmxpmx>

      end subroutine advection_scal_dz

!.................................................................

      subroutine advection_scal_dtisub
      use module_principal
      use module_parallele
      implicit none

! Cet algo suppose que la couche merged est 100% melangee avec les couches sous-jacente
! et donc que l'on a affaire A une couche dont l'epaisseur peut etre consideree comme
! z_w(kmerged+1)-z_w(1) et non pas dz(kmerged) (methode 2). Si on garde dz(kmerged) (methode 1) l'algo considere
! un volume plus petit pour plus de securite (mais evntuellement un dt inutilement trop petit)
! D'autre part les flux qui entrent dans ce volume sont ceux du niveau k=kmerged mais egalement
! ceux sous kmerged, kundermin_t etant le niveau du premier flux lateral non nul
      x1=0.
! Niveaux standard:
      do j=1,jmax ; do i=1,imax
        do k=kmerged_t(i,j)+1,kmax
         x1=max(x1,max(max(max(abs(veldydz_u(i  ,j  ,k,1))  &
                              ,abs(veldydz_u(i+1,j  ,k,1))) &
                              ,abs(veldxdz_v(i  ,j  ,k,1))) &
                              ,abs(veldxdz_v(i  ,j+1,k,1))) & 
                    *dti_lp                                 &
                    *wetmask_t(i,j)                         &
                       /anyv3d(i,j,k,id_dxdydz))
        enddo
      enddo       ; enddo
! Niveaux <= kmerged
      do j=1,jmax ; do i=1,imax
        sum1=0. ; sum2=0. ; sum3=0. ; sum4=0. ; sum5=0. ; sum6=0.
        do k=kundermin_t(i,j),kmerged_t(i,j) !07-12-17
         sum1=sum1+veldydz_u(i  ,j  ,k,1)
         sum2=sum2+veldydz_u(i+1,j  ,k,1)
         sum3=sum3+veldxdz_v(i  ,j  ,k,1)
         sum4=sum4+veldxdz_v(i  ,j+1,k,1)
         sum5=sum5     +dz_t(i  ,j  ,k,0)
         sum6=sum6     +dz_t(i  ,j  ,k,2)
        enddo
         x1=max(x1,max(max(max(abs(sum1)  &
                              ,abs(sum2)) &
                              ,abs(sum3)) &
                              ,abs(sum4)) & 
                    *dti_lp                                  &
                    *wetmask_t(i,j)                          &
!                      /anyv3d(i,j,kmerged_t(i,j),id_dxdydz) &          ! methode 1
                      /(dxdy_t(i,j)*min(sum5,sum6))          &          ! methode 2
               ) ! fermeture max
      enddo       ; enddo


      call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_max,par%comm2d,ierr)
      loopmaxts=int(x2+1) ! =min(int(x2)+1,4)
      invloopmaxts=1./real(loopmaxts)
!     write(6,*)x2,loopmaxts
!     stop 'fififif'
      dti_lpsub=dti_lp/loopmaxts

!-----------------------------------------------------------------------
      if(loopmaxts>50.or.loopmaxts<0) then !pmxpmx> !04-09-16!15-05-17

! On refait le exactement meme test mais on examine point par point le depassement !05-05-18
! Niveaux standard:
      do j=1,jmax ; do i=1,imax
        do k=kmerged_t(i,j)+1,kmax

         x2=       max(max(max(abs(veldydz_u(i  ,j  ,k,1))  &
                              ,abs(veldydz_u(i+1,j  ,k,1))) &
                              ,abs(veldxdz_v(i  ,j  ,k,1))) &
                              ,abs(veldxdz_v(i  ,j+1,k,1))) & 
                    *dti_lp                                 &
                    *wetmask_t(i,j)                         &
                       /anyv3d(i,j,k,id_dxdydz) 

             if(int(x2+1)>50) then !(0v0)>
               write(10+par%rank,*)'---------------'
               write(10+par%rank,*)'x2=           ',x2
               write(10+par%rank,*)'i,j,k loc     ',i,j,k
               write(10+par%rank,*)'i,j,k glb     ',i+par%timax(1),j+par%tjmax(1),k
               write(10+par%rank,*)'hz_w(i,j,0:2) ',real(hz_w(i,j,0:2))
               write(10+par%rank,*)'dz_t(i,j,k,1) ',dz_t(i,j,k,1)
               write(10+par%rank,*)'dz "anyv3d"   ',anyv3d(i,j,k,id_dxdydz)/dxdy_t(i,j)
               write(10+par%rank,*)'volume        ',anyv3d(i,j,k,id_dxdydz)
               write(10+par%rank,*)'kmin_w kmerge ',kmin_w(i,j),kmerged_t(i,j)
               write(10+par%rank,*)'wetmask_t(i,j)',wetmask_t(i,j)
               write(10+par%rank,*)'veldydz_u(i  ,j  ,k,1)',veldydz_u(i  ,j  ,k,1)
               write(10+par%rank,*)'    vel_u(i  ,j  ,k,1)',vel_u(i  ,j  ,k,1)
               write(10+par%rank,*)'veldydz_u(i+1,j  ,k,1)',veldydz_u(i+1,j  ,k,1)
               write(10+par%rank,*)'    vel_u(i+1,j  ,k,1)',vel_u(i+1,j  ,k,1)
               write(10+par%rank,*)'veldxdz_v(i  ,j  ,k,1)',veldxdz_v(i ,j  ,k,1)
               write(10+par%rank,*)'    vel_v(i  ,j  ,k,1)',vel_v(i  ,j  ,k,1)
               write(10+par%rank,*)'veldxdz_v(i  ,j+1,k,1)',veldxdz_v(i ,j+1,k,1)
               write(10+par%rank,*)'    vel_v(i  ,j+1,k,1)',vel_v(i  ,j+1,k,1)
             endif                 !(0v0)>

        enddo
      enddo       ; enddo
! Niveaux <= kmerged
      do j=1,jmax ; do i=1,imax
        sum1=0. ; sum2=0. ; sum3=0. ; sum4=0. ; sum5=0. ; sum6=0.
        do k=kundermin_t(i,j),kmerged_t(i,j) !07-12-17
         sum1=sum1+veldydz_u(i  ,j  ,k,1)
         sum2=sum2+veldydz_u(i+1,j  ,k,1)
         sum3=sum3+veldxdz_v(i  ,j  ,k,1)
         sum4=sum4+veldxdz_v(i  ,j+1,k,1)
         sum5=sum5     +dz_t(i  ,j  ,k,0)
         sum6=sum6     +dz_t(i  ,j  ,k,2)
        enddo
         x2=       max(max(max(abs(sum1)  &
                              ,abs(sum2)) &
                              ,abs(sum3)) &
                              ,abs(sum4)) & 
                    *dti_lp                                  &
                    *wetmask_t(i,j)                          &
                      /(dxdy_t(i,j)*min(sum5,sum6))                     ! methode 2

             if(int(x2+1)>50) then !(0v0)>
               write(10+par%rank,*)'---------------'
               write(10+par%rank,*)'x2=           ',x2
               write(10+par%rank,*)'i,j,k loc     ',i,j,k
               write(10+par%rank,*)'i,j,k glb     ',i+par%timax(1),j+par%tjmax(1),k
               write(10+par%rank,*)'hz_w(i,j,0:2) ',real(hz_w(i,j,0:2))
               write(10+par%rank,*)'dz_t(i,j,k,1) ',dz_t(i,j,k,1)
               write(10+par%rank,*)'dz "anyv3d"   ',anyv3d(i,j,k,id_dxdydz)/dxdy_t(i,j)
               write(10+par%rank,*)'volume        ',anyv3d(i,j,k,id_dxdydz)
               write(10+par%rank,*)'kmin_w kmerge ',kmin_w(i,j),kmerged_t(i,j)
               write(10+par%rank,*)'wetmask_t(i,j)',wetmask_t(i,j)
               write(10+par%rank,*)'veldydz_u(i  ,j  ,k,1)',veldydz_u(i  ,j  ,k,1)
               write(10+par%rank,*)'    vel_u(i  ,j  ,k,1)',vel_u(i  ,j  ,k,1)
               write(10+par%rank,*)'veldydz_u(i+1,j  ,k,1)',veldydz_u(i+1,j  ,k,1)
               write(10+par%rank,*)'    vel_u(i+1,j  ,k,1)',vel_u(i+1,j  ,k,1)
               write(10+par%rank,*)'veldxdz_v(i  ,j  ,k,1)',veldxdz_v(i ,j  ,k,1)
               write(10+par%rank,*)'    vel_v(i  ,j  ,k,1)',vel_v(i  ,j  ,k,1)
               write(10+par%rank,*)'veldxdz_v(i  ,j+1,k,1)',veldxdz_v(i ,j+1,k,1)
               write(10+par%rank,*)'    vel_v(i  ,j+1,k,1)',vel_v(i  ,j+1,k,1)
             endif                 !(0v0)>

      enddo       ; enddo

           call graph_out

           if(loopmaxts<0) then !>>>
            stop 'loopmaxts<0'
           else                 !>>>
            write(6,*)'loopmaxts ',loopmaxts,x2
            stop 'loopmaxts>50 see fort.XXX error files'
           endif                !>>>

!         endif                 !ooo>
      endif                !pmxpmx>
!-----------------------------------------------------------------------

      if(par%rank==0) then
        open(unit=3,file=trim(tmpdirname)//'dti_ts',position='append') !17-07-20
         write(3,*)real(elapsedtime_now/86400.),loopmaxts,real(x2),real(dti_lpsub)
        close(3)
      endif
      end subroutine advection_scal_dtisub

!..............................................................................

      subroutine advection_scal_upm_subt  ! 05-01-16
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      double precision :: velef_,flagmsk_=1.,scalefac_,cst_
      integer :: t_=1        &
                ,id_flx_=1   &  ! x flux identifier
                ,id_fly_=1   &  ! y flux identifier
                ,id_tbef_=3  &  ! temperature    "before" identifier
!               ,id_tdiv=4  &  ! T*divergence(velocity) !26-01-17
                ,loop_
#ifdef synopsis
       subroutinetitle='advection_scal_upm_subs'
       subroutinedescription=                                     &
          'Computes temperature and teminity advection diffusion' &
       //' horizontal fluxes using a "3 points" advection scheme' &
       //' hybrid QUICK-UPWIND'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Cas particulier du modele 1DV
!     if(flag_1dv==1) then !>>>>> !14-02-16
!      tem_t(:,:,:,2)=tem_t(:,:,:,0)
!      return
!     endif

#ifdef bidon
      sum1=0.
      sum2=0.
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
       sum1=sum1+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)
       sum2=sum2+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*tem_t(i,j,k,1)
      enddo
      enddo
      enddo
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      if(par%rank==0)write(6,*)sum2glb,sum2glb/sum1glb
#endif



      cst_=0.5*cst_adv_hor/3. !11-09-17


!#ifdef bidon
      do k=1,kmax ; do j=-1,jmax+2 ; do i=-1,imax+2
       anyv3d(i,j,k,id_tbef_)=tem_t(i,j,k,before)
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,id_tdiv)=0. ! reset cumul temperature*divergence(courant) !26-01-17
      enddo ; enddo ; enddo

      do loop_=1,loopmaxts !*********>

! Flux Oi
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax
      do i=1,imax+1

      anyv3d(i,j,k,id_flx_)=0.5*( & !prime>


! Schema QUICK en dehors influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      anyv3d(i,j,k,id_hybcoefu)*( &
! Partie advective QUICK:
      -veldydz_u(i,j,k,1) &
!       *(anyv3d(i,j,k,id_tnow_)+anyv3d(i-1,j,k,id_tnow_)) &
! Si le nombre de courant est faible il n'y a pas de risque A maintenir la
! partie centree sur tem_t(now) ce qui permet de maintenir un leap-frog centrE
! classique meme en subdivisant le pas de temps, ce qui permet A priori de 
! maintenir les proprietes de conservation associeEs
        *(    tem_t(i,j,k,now)+tem_t(i-1,j,k,now)      &
#ifdef bidonref
       +2.*temref_u(i,j,k)                             &
          -temref_t(i,j,k) -temref_t(i-1,j,k)          &
#endif
                                                 )     &

! Partie diffusive QUICK remise A jour aux sous-pas de temps advectif:
       +cst_*0.5*(anyv3d(i,j,k,id_ofactort)+anyv3d(i,j+1,k,id_ofactort)) &

! efficacite densitaire
!       *abs(-alp_t*(   tem_t(i,j,k,1)-tem_t(i-1,j,k,1))+alp_s*(   sal_t(i,j,k,1)-sal_t(i-1,j,k,1))) &
!       /   (+alp_t*abs(tem_t(i,j,k,1)-tem_t(i-1,j,k,1))+alp_s*abs(sal_t(i,j,k,1)-sal_t(i-1,j,k,1))+small1) &

                                                               *( & !m°v°m>
             (veldydz_u(i,j,k,1)+abs(veldydz_u(i,j,k,1)))*( &           !pmxpmx>

               (   anyv3d(i  ,j,k,id_tbef_)   -anyv3d(i-1,j,k,id_tbef_)  &

! correction anyv3d(i  ,j,k,id_tbef_) par rapport niveau z(i-1,j,k)
      +( anyv3d(i  ,j,kp1,id_tbef_)- anyv3d(i  ,j,km1,id_tbef_))  &
      /(depth_t(i  ,j,kp1)         -depth_t(i  ,j,km1)         )  &
      *(depth_t(i-1,j,k  )         -depth_t(i  ,j,k  )         )  &

               ) &

              -(   anyv3d(i-1,j,k,id_tbef_)-anyv3d(i-2,j,k,id_tbef_)  &

! correction anyv3d(i-2,j,k,id_tbef_) par rapport niveau z(i-1,j,k)
      -( anyv3d(i-2,j,kp1,id_tbef_)- anyv3d(i-2,j,km1,id_tbef_))  &
      /(depth_t(i-2,j,kp1)         -depth_t(i-2,j,km1)         )  &
      *(depth_t(i-1,j,k  )         -depth_t(i-2,j,k  )         )  &

               ) &
                                           *mask_t(i-2,j,kmax)      ) & !pmxpmx>

            +(veldydz_u(i,j,k,1)-abs(veldydz_u(i,j,k,1)))*( &           !pmxpmx>

               (   anyv3d(i+1,j,k,id_tbef_)-anyv3d(i,j,k,id_tbef_)    &

! correction anyv3d(i+1,j,k,id_tbef_) par rapport niveau z(i  ,j,k)
      +( anyv3d(i+1,j,kp1,id_tbef_)- anyv3d(i+1,j,km1,id_tbef_))  &
      /(depth_t(i+1,j,kp1)         -depth_t(i+1,j,km1)         )  &
      *(depth_t(i  ,j,k  )         -depth_t(i+1,j,k  )         )  &

               )   &
                  *mask_t(i+1,j,kmax)                                 &

              -(   anyv3d(i,j,k,id_tbef_)-anyv3d(i-1,j,k,id_tbef_)    &

! correction anyv3d(i-1,j,k,id_tbef_) par rapport niveau z(i  ,j,k)
      -( anyv3d(i-1,j,kp1,id_tbef_)- anyv3d(i-1,j,km1,id_tbef_))  &
      /(depth_t(i-1,j,kp1)         -depth_t(i-1,j,km1)         )  &
      *(depth_t(i  ,j,k  )         -depth_t(i-1,j,k  )         )  &

                                                                  ))  & !pmxpmx>

                                                                )  & !m°v°m>()cst_
                                                                )  & !()anyv3d(i,j,k,id_hybcoefu)
! Schema UPWIND dans les zones sous influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      +(1.-anyv3d(i,j,k,id_hybcoefu))*( & !pmx>
           -veldydz_u(i,j,k,1)* ( & !ppp>
               anyv3d(i,j,k,id_tbef_)+anyv3d(i-1,j,k,id_tbef_)   &
                                ) & !ppp>
       +abs(veldydz_u(i,j,k,1))*( & !mmm>
               anyv3d(i,j,k,id_tbef_)-anyv3d(i-1,j,k,id_tbef_)   &
                                ) & !mmm>
                                      ) & !pmx>

       ) !prime>

!     -veldydz_u(i,j,k,1)*dti_lpsub*temref_u(i,j,k)


      enddo !k
      enddo !i
      enddo !j

! Partial Oi Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_tdiv)= & !26-01-17
       anyv3d(i,j,k,id_tdiv)  &
      +anyv3d(i,j,k,id_tbef_)*dti_lpsub*(veldydz_u(i+1,j,k,1)           &
                                        -veldydz_u(i  ,j,k,1))/dxdy_t(i,j)

!       anyv3d(i,j,k,id_tbef_)                       &
        anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !ooo> !25-04-17

        anyv3d(i,j,k,id_tbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i+1,j,k,id_flx_)                &
             -anyv3d(i  ,j,k,id_flx_)                &

         +(veldydz_u(i+1,j,k,1)                      &
          -veldydz_u(i  ,j,k,1))*anyv3d(i,j,k,id_tbef_) & !26-01-17


         )/dxdy_t(i,j) ) & !ooo>

         +(1-mask_t(i,j,k))*tem_t(i,j,k,before)*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

! Note: le terme anyv3d(i,j,k,id_tdiv) ne fait pas l'objet d'un traitement similaire pour le masque car
! si mask_t=0 alors veldydz_u(i:i+1,...)=0 (question: et les fleuves ?)
      enddo
      enddo
      enddo

!     call my_outputs_zonesalttempflux('it',id_flx_)

      call vertmix_merged_levels_any(0,id_tbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_tbef_)
      call obc_int_anyv3d(id_tbef_,'zb')

! Flux Oj
      cst_=0.5*cst_adv_hor/3.
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax+1
      do i=1,imax

      anyv3d(i,j,k,id_fly_)=0.5*( & !prime>

! Schema QUICK en dehors influence des fleuves
      anyv3d(i,j,k,id_hybcoefv)*( &
! Partie advective QUICK:
      -veldxdz_v(i,j,k,1)  &
!       *(anyv3d(i,j,k,id_tnow_)+anyv3d(i,j-1,k,id_tnow_)) &
! Si le nombre de courant est faible il n'y a pas de risque A maintenir la
! partie centree sur tem_t(now) ce qui permet de maintenir un leap-frog centrE
! classique meme en subdivisant le pas de temps ce qui permet A priori de 
! maintenir les proprietes de conservation associeE
        *(    tem_t(i,j,k,now)+tem_t(i,j-1,k,now)  &
#ifdef bidonref
       +2.*temref_v(i,j,k)                         &
          -temref_t(i,j,k) -temref_t(i,j-1,k)      &
#endif
                                                 ) &

! Partie diffusive QUICK:
       +cst_*0.5*(anyv3d(i,j,k,id_ofactort)+anyv3d(i+1,j,k,id_ofactort)) &

! efficacite densitaire
!       *abs(-alp_t*(   tem_t(i,j,k,1)-tem_t(i,j-1,k,1))+alp_s*(   sal_t(i,j,k,1)-sal_t(i,j-1,k,1))) &
!       /   (+alp_t*abs(tem_t(i,j,k,1)-tem_t(i,j-1,k,1))+alp_s*abs(sal_t(i,j,k,1)-sal_t(i,j-1,k,1))+small1) &

                                                                  *(  & !m°v°m>

             (veldxdz_v(i,j,k,1)+abs(veldxdz_v(i,j,k,1)))*(  &       !ooo>

            (   anyv3d(i,j  ,k,id_tbef_)-anyv3d(i,j-1,k,id_tbef_)  &

! correction anyv3d(i,j  ,k,id_tbef_) par rapport niveau z(i,j-1,k)
      +( anyv3d(i,j  ,kp1,id_tbef_)- anyv3d(i,j  ,km1,id_tbef_))  &
      /(depth_t(i,j  ,kp1)         -depth_t(i,j  ,km1)         )  &
      *(depth_t(i,j-1,k  )         -depth_t(i,j  ,k  )         )  &

            ) &

           -(   anyv3d(i,j-1,k,id_tbef_)-anyv3d(i,j-2,k,id_tbef_)  &

! correction anyv3d(i,j-2,k,id_tbef_) par rapport niveau z(i,j-1,k)
      -( anyv3d(i,j-2,kp1,id_tbef_)- anyv3d(i,j-2,km1,id_tbef_))  &
      /(depth_t(i,j-2,kp1)         -depth_t(i,j-2,km1)         )  &
      *(depth_t(i,j-1,k  )         -depth_t(i,j-2,k  )         )  &

            ) &
                                        *mask_t(i,j-2,kmax)       )& !ooo>

            +(veldxdz_v(i,j,k,1)-abs(veldxdz_v(i,j,k,1)))*(  &       !ooo>

            (   anyv3d(i,j+1,k,id_tbef_)-anyv3d(i,j,k,id_tbef_)    &

! correction anyv3d(i,j+1,k,id_tbef_) par rapport niveau z(i,j  ,k)
      +( anyv3d(i,j+1,kp1,id_tbef_)- anyv3d(i,j+1,km1,id_tbef_))  &
      /(depth_t(i,j+1,kp1)         -depth_t(i,j+1,km1)         )  &
      *(depth_t(i,j  ,k  )         -depth_t(i,j+1,k  )         )  &

            )   &
               *mask_t(i,j+1,kmax)                                 &

           -(   anyv3d(i,j,k,id_tbef_)-anyv3d(i,j-1,k,id_tbef_)    &

! correction anyv3d(i,j-1,k,id_tbef_) par rapport niveau z(i,j  ,k)
      -( anyv3d(i,j-1,kp1,id_tbef_)- anyv3d(i,j-1,km1,id_tbef_))  &
      /(depth_t(i,j-1,kp1)         -depth_t(i,j-1,km1)         )  &
      *(depth_t(i,j  ,k  )         -depth_t(i,j-1,k  )         )  &

                                                               ))  & !ooo>

                                                                   )  & !m°v°m>()cst_
                                                                   )  & !()anyv3d(i,j,k,id_hybcoefv)
! Schema UPWIND dans les zones sous influence des fleuves:
      +(1.-anyv3d(i,j,k,id_hybcoefv))*( & !pmx>
          -veldxdz_v(i,j,k,1)* ( & !ppp>
              anyv3d(i,j,k,id_tbef_)+anyv3d(i,j-1,k,id_tbef_)   &
                               ) & !ppp>
      +abs(veldxdz_v(i,j,k,1))*( & !mmm>
              anyv3d(i,j,k,id_tbef_)-anyv3d(i,j-1,k,id_tbef_)   &
                               ) & !mmm>
                                      ) & !pmx>

       ) !prime>

!      -veldxdz_v(i,j,k,1)*dti_lpsub*temref_v(i,j,k)

      enddo !k
      enddo !i
      enddo !j

! Partial Horizontal Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_tdiv)= & !26-01-17
       anyv3d(i,j,k,id_tdiv)  &
      +anyv3d(i,j,k,id_tbef_)*dti_lpsub*(veldxdz_v(i,j+1,k,1)                 &
                                        -veldxdz_v(i,j  ,k,1))/dxdy_t(i,j)

        anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !ooo>
        anyv3d(i,j,k,id_tbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i,j+1,k,id_fly_)                &
             -anyv3d(i,j  ,k,id_fly_)                &

         +(veldxdz_v(i,j+1,k,1)                      &
          -veldxdz_v(i,j  ,k,1))*anyv3d(i,j,k,id_tbef_) & !26-01-17

         )/dxdy_t(i,j)  ) & !ooo>
 
         +(1-mask_t(i,j,k))*tem_t(i,j,k,before)*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

      enddo
      enddo
      enddo

!     call my_outputs_zonesalttempflux('jt',id_fly_) 

      call vertmix_merged_levels_any(0,id_tbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_tbef_)
      call obc_int_anyv3d(id_tbef_,'zb')

!#endif
      enddo               !*********>

      do k=1,kmax ; do j=1,jmax  ; do i=1,imax  

           tem_t(i,j,k,2)                              &
       =  anyv3d(i,j,k,id_tbef_)                    !  &

!     -wetmask_t(i,j)*tem_t(i,j,k,now)*dti_lp*(        &

!                                              veldydz_u(i+1,j,k,1)   &
!                                             -veldydz_u(i  ,j,k,1)   &
!                                             +veldxdz_v(i,j+1,k,1)   &
!                                             -veldxdz_v(i,j  ,k,1)   &

!                                             )/(dz_t(i,j,k,before)*dxdy_t(i,j))

      enddo ; enddo ; enddo

      end subroutine advection_scal_upm_subt


!..............................................................................

      subroutine advection_scal_upm_subs  ! 05-01-16
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      double precision :: velef_,flagmsk_=1.,scalefac_,cst_
      integer :: t_=1        &
                ,id_flx_=1   &  ! x flux identifier
                ,id_fly_=1   &  ! y flux identifier
                ,id_sbef_=3  &  ! salinity    "before" identifier
!               ,id_sdiv=5  &  ! S*divergence(velocity)
                ,loop_
#ifdef synopsis
       subroutinetitle='advection_scal_upm_subs'
       subroutinedescription=                                     &
          'Computes salperature and salinity advection diffusion' &
       //' horizontal fluxes using a "3 points" advection scheme' &
       //' hybrid QUICK-UPWIND'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Cas particulier du modele 1DV
!     if(flag_1dv==1) then !>>>>> !14-02-16
!      sal_t(:,:,:,2)=sal_t(:,:,:,0)
!      return
!     endif

      cst_=0.5*cst_adv_hor/3.


!#ifdef bidon
      do k=1,kmax ; do j=-1,jmax+2 ; do i=-1,imax+2
       anyv3d(i,j,k,id_sbef_)=sal_t(i,j,k,before) 
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,id_sdiv)=0. ! reset cumul salinity*divergence(courant)
      enddo ; enddo ; enddo


      do loop_=1,loopmaxts !*********>


! Flux Oi
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax
      do i=1,imax+1

      anyv3d(i,j,k,id_flx_)=0.5*( & !prime>


! Schema QUICK en dehors influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      anyv3d(i,j,k,id_hybcoefu)*( &
! Partie advective QUICK:
      -veldydz_u(i,j,k,1) &
!       *(anyv3d(i,j,k,id_tnow_)+anyv3d(i-1,j,k,id_tnow_)) &
! Si le nombre de courant est faible il n'y a pas de risque A maintenir la
! partie centree sur sal_t(now) ce qui permet de maintenir un leap-frog centrE
! classique meme en subdivisant le pas de salps, ce qui permet A priori de 
! maintenir les proprietes de conservation associeEs
        *(    sal_t(i,j,k,now)+sal_t(i-1,j,k,now)      &
#ifdef bidonref
       +2.*salref_u(i,j,k)                             &
          -salref_t(i,j,k) -salref_t(i-1,j,k)          &
#endif
                                                 )     &

! Partie diffusive QUICK remise A jour aux sous-pas de salps advectif:
       +cst_*0.5*(anyv3d(i,j,k,id_ofactors)+anyv3d(i,j+1,k,id_ofactors))  &

! Efficacite densitaire
!       *abs(-alp_t*(   tem_t(i,j,k,1)-tem_t(i-1,j,k,1))+alp_s*(   sal_t(i,j,k,1)-sal_t(i-1,j,k,1))) &
!       /   (+alp_t*abs(tem_t(i,j,k,1)-tem_t(i-1,j,k,1))+alp_s*abs(sal_t(i,j,k,1)-sal_t(i-1,j,k,1))+small1) &

                                                               *( & !m°v°m>
             (veldydz_u(i,j,k,1)+abs(veldydz_u(i,j,k,1)))*( &           !pmxpmx>

               (   anyv3d(i  ,j,k,id_sbef_)-anyv3d(i-1,j,k,id_sbef_)  &

! correction anyv3d(i  ,j,k,id_sbef_) par rapport niveau z(i-1,j,k)
      +( anyv3d(i  ,j,kp1,id_sbef_)- anyv3d(i  ,j,km1,id_sbef_))  &
      /(depth_t(i  ,j,kp1)         -depth_t(i  ,j,km1)         )  &
      *(depth_t(i-1,j,k  )         -depth_t(i  ,j,k  )         )  &

               ) &

              -(   anyv3d(i-1,j,k,id_sbef_)-anyv3d(i-2,j,k,id_sbef_)  &

! correction anyv3d(i-2,j,k,id_sbef_) par rapport niveau z(i-1,j,k)
      -( anyv3d(i-2,j,kp1,id_sbef_)- anyv3d(i-2,j,km1,id_sbef_))  &
      /(depth_t(i-2,j,kp1)         -depth_t(i-2,j,km1)         )  &
      *(depth_t(i-1,j,k  )         -depth_t(i-2,j,k  )         )  &

               ) &
                                           *mask_t(i-2,j,kmax)      ) & !pmxpmx>

            +(veldydz_u(i,j,k,1)-abs(veldydz_u(i,j,k,1)))*( &           !pmxpmx>

               (   anyv3d(i+1,j,k,id_sbef_)-anyv3d(i,j,k,id_sbef_)    &

! correction anyv3d(i+1,j,k,id_sbef_) par rapport niveau z(i  ,j,k)
      +( anyv3d(i+1,j,kp1,id_sbef_)- anyv3d(i+1,j,km1,id_sbef_))  &
      /(depth_t(i+1,j,kp1)         -depth_t(i+1,j,km1)         )  &
      *(depth_t(i  ,j,k  )         -depth_t(i+1,j,k  )         )  &

               )   &
                  *mask_t(i+1,j,kmax)                                 &

              -(   anyv3d(i,j,k,id_sbef_)-anyv3d(i-1,j,k,id_sbef_)    &

! correction anyv3d(i-1,j,k,id_sbef_) par rapport niveau z(i  ,j,k)
      -( anyv3d(i-1,j,kp1,id_sbef_)- anyv3d(i-1,j,km1,id_sbef_))  &
      /(depth_t(i-1,j,kp1)         -depth_t(i-1,j,km1)         )  &
      *(depth_t(i  ,j,k  )         -depth_t(i-1,j,k  )         )  &

                                                                  ))  & !pmxpmx>

                                                                )  & !m°v°m>()cst_
                                                                )  & !()anyv3d(i,j,k,id_hybcoefu)
! Schema UPWIND dans les zones sous influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      +(1.-anyv3d(i,j,k,id_hybcoefu))*( & !pmx>
           -veldydz_u(i,j,k,1)* ( & !ppp>
               anyv3d(i,j,k,id_sbef_)+anyv3d(i-1,j,k,id_sbef_)   &
                                ) & !ppp>
       +abs(veldydz_u(i,j,k,1))*( & !mmm>
               anyv3d(i,j,k,id_sbef_)-anyv3d(i-1,j,k,id_sbef_)   &
                                ) & !mmm>
                                      ) & !pmx>

       ) !prime>


      enddo !k
      enddo !i
      enddo !j

! Partial Oi Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax


! Cumul du pivot obligatoirement avant advection partielle 06-02-17
! Note que pour economiser un peu de cpu la multiplication par wetmask_t(i,j)/dz_t(i,j,k,before),
! commune a toutes les operations, sera appliquee une fois pour toute a la fin...
       anyv3d(i,j,k,id_sdiv)= &
       anyv3d(i,j,k,id_sdiv)  &
      +anyv3d(i,j,k,id_sbef_)*dti_lpsub*(veldydz_u(i+1,j,k,1)           &
                                        -veldydz_u(i  ,j,k,1))/dxdy_t(i,j)

        anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !ooo>
        anyv3d(i,j,k,id_sbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i+1,j,k,id_flx_)                &
             -anyv3d(i  ,j,k,id_flx_)                &

         +(veldydz_u(i+1,j,k,1)                      &
          -veldydz_u(i  ,j,k,1))*anyv3d(i,j,k,id_sbef_) & !25-01-17

         )/dxdy_t(i,j)  ) & !ooo>

         +(1-mask_t(i,j,k))*sal_t(i,j,k,before)*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

      enddo
      enddo
      enddo

!     call my_outputs_zonesalttempflux('is',id_flx_) !29-01-19

      call vertmix_merged_levels_any(0,id_sbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_sbef_)
      call obc_int_anyv3d(id_sbef_,'zb')

! Flux Oj
      cst_=0.5*cst_adv_hor/3.
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax+1
      do i=1,imax

      anyv3d(i,j,k,id_fly_)=0.5*( & !prime>

! Schema QUICK en dehors influence des fleuves
      anyv3d(i,j,k,id_hybcoefv)*( &
! Partie advective QUICK:
      -veldxdz_v(i,j,k,1)  &
!       *(anyv3d(i,j,k,id_tnow_)+anyv3d(i,j-1,k,id_tnow_)) &
! Si le nombre de courant est faible il n'y a pas de risque A maintenir la
! partie centree sur sal_t(now) ce qui permet de maintenir un leap-frog centrE
! classique meme en subdivisant le pas de salps ce qui permet A priori de 
! maintenir les proprietes de conservation associeE
        *(    sal_t(i,j,k,now)+sal_t(i,j-1,k,now)  &
#ifdef bidonref
       +2.*salref_v(i,j,k)                         &
          -salref_t(i,j,k) -salref_t(i,j-1,k)      &
#endif
                                                 ) &

! Partie diffusive QUICK:
       +cst_*0.5*(anyv3d(i,j,k,id_ofactors)+anyv3d(i+1,j,k,id_ofactors))  &

! Efficacite densitaire
!       *abs(-alp_t*(   tem_t(i,j,k,1)-tem_t(i,j-1,k,1))+alp_s*(   sal_t(i,j,k,1)-sal_t(i,j-1,k,1))) &
!       /   (+alp_t*abs(tem_t(i,j,k,1)-tem_t(i,j-1,k,1))+alp_s*abs(sal_t(i,j,k,1)-sal_t(i,j-1,k,1))+small1) &

                                                                  *( & !m°v°m>

             (veldxdz_v(i,j,k,1)+abs(veldxdz_v(i,j,k,1)))*(  &       !ooo>

            (   anyv3d(i,j  ,k,id_sbef_)-anyv3d(i,j-1,k,id_sbef_)  &

! correction anyv3d(i,j  ,k,id_sbef_) par rapport niveau z(i,j-1,k)
      +( anyv3d(i,j  ,kp1,id_sbef_)- anyv3d(i,j  ,km1,id_sbef_))  &
      /(depth_t(i,j  ,kp1)         -depth_t(i,j  ,km1)         )  &
      *(depth_t(i,j-1,k  )         -depth_t(i,j  ,k  )         )  &

            ) &

           -(   anyv3d(i,j-1,k,id_sbef_)-anyv3d(i,j-2,k,id_sbef_)  &

! correction anyv3d(i,j-2,k,id_sbef_) par rapport niveau z(i,j-1,k)
      -( anyv3d(i,j-2,kp1,id_sbef_)- anyv3d(i,j-2,km1,id_sbef_))  &
      /(depth_t(i,j-2,kp1)         -depth_t(i,j-2,km1)         )  &
      *(depth_t(i,j-1,k  )         -depth_t(i,j-2,k  )         )  &

            ) &
                                        *mask_t(i,j-2,kmax)       )& !ooo>

            +(veldxdz_v(i,j,k,1)-abs(veldxdz_v(i,j,k,1)))*(  &       !ooo>

            (   anyv3d(i,j+1,k,id_sbef_)-anyv3d(i,j,k,id_sbef_)    &

! correction anyv3d(i,j+1,k,id_sbef_) par rapport niveau z(i,j  ,k)
      +( anyv3d(i,j+1,kp1,id_sbef_)- anyv3d(i,j+1,km1,id_sbef_))  &
      /(depth_t(i,j+1,kp1)         -depth_t(i,j+1,km1)         )  &
      *(depth_t(i,j  ,k  )         -depth_t(i,j+1,k  )         )  &

            )   &
               *mask_t(i,j+1,kmax)                                 &

           -(   anyv3d(i,j,k,id_sbef_)-anyv3d(i,j-1,k,id_sbef_)    &

! correction anyv3d(i,j-1,k,id_sbef_) par rapport niveau z(i,j  ,k)
      -( anyv3d(i,j-1,kp1,id_sbef_)- anyv3d(i,j-1,km1,id_sbef_))  &
      /(depth_t(i,j-1,kp1)         -depth_t(i,j-1,km1)         )  &
      *(depth_t(i,j  ,k  )         -depth_t(i,j-1,k  )         )  &

                                                               ))  & !ooo>

                                                                   )  & !m°v°m>()cst_
                                                                   )  & !()anyv3d(i,j,k,id_hybcoefv)
! Schema UPWIND dans les zones sous influence des fleuves:
      +(1.-anyv3d(i,j,k,id_hybcoefv))*( & !pmx>
          -veldxdz_v(i,j,k,1)* ( & !ppp>
              anyv3d(i,j,k,id_sbef_)+anyv3d(i,j-1,k,id_sbef_)   &
                               ) & !ppp>
      +abs(veldxdz_v(i,j,k,1))*( & !mmm>
              anyv3d(i,j,k,id_sbef_)-anyv3d(i,j-1,k,id_sbef_)   &
                               ) & !mmm>
                                      ) & !pmx>

       ) !prime>

      enddo !k
      enddo !i
      enddo !j

! Partial Horizontal Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_sdiv)= &
       anyv3d(i,j,k,id_sdiv)  &
      +anyv3d(i,j,k,id_sbef_)*dti_lpsub*(veldxdz_v(i,j+1,k,1)   &
                                        -veldxdz_v(i,j  ,k,1))/dxdy_t(i,j)

        anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !ooo>
        anyv3d(i,j,k,id_sbef_)*dz_t(i,j,k,before)   & !01-05-19    

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i,j+1,k,id_fly_)                &
             -anyv3d(i,j  ,k,id_fly_)                &

         +(veldxdz_v(i,j+1,k,1)                      &
          -veldxdz_v(i,j  ,k,1))*anyv3d(i,j,k,id_sbef_) & !25-01-17

         )/dxdy_t(i,j)  ) & !ooo>

         +(1-mask_t(i,j,k))*sal_t(i,j,k,before)*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

      enddo
      enddo
      enddo

!     call my_outputs_zonesalttempflux('js',id_fly_) 
      call vertmix_merged_levels_any(0,id_sbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_sbef_)
      call obc_int_anyv3d(id_sbef_,'zb')

!#endif
      enddo               !*********>

      do k=1,kmax ; do j=1,jmax  ; do i=1,imax  

           sal_t(i,j,k,2)                              &
       =  anyv3d(i,j,k,id_sbef_)                    !  &

!     -wetmask_t(i,j)*sal_t(i,j,k,now)*dti_lp*(        &

!                                              veldydz_u(i+1,j,k,1)   &
!                                             -veldydz_u(i  ,j,k,1)   &
!                                             +veldxdz_v(i,j+1,k,1)   &
!                                             -veldxdz_v(i,j  ,k,1)   &

!                                             )/(dz_t(i,j,k,before)*dxdy_t(i,j))

      enddo ; enddo ; enddo

      end subroutine advection_scal_upm_subs

!..............................................................................
#ifdef bidon
      subroutine advection_scal_upmx_subt  ! 05-01-16
      use module_principal
      use module_parallele
      implicit none
      double precision :: velef_,flagmsk_=1.,scalefac_,cst_
      integer :: t_=1        &
                ,id_flx_=1   &  ! x flux identifier
                ,id_fly_=1   &  ! y flux identifier
                ,id_tbef_=3  &  ! temperature    "before" identifier
!               ,id_tdiv=4  &  ! T*divergence(velocity) !26-01-17
                ,loop_
#ifdef synopsis
       subroutinetitle='advection_scal_upmx_subs'
       subroutinedescription=                                     &
          'Computes temperature and teminity advection diffusion' &
       //' horizontal fluxes using a "3 points" advection scheme' &
       //' hybrid QUICK-UPWIND'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Cas particulier du modele 1DV
!     if(flag_1dv==1) then !>>>>> !14-02-16
!      tem_t(:,:,:,2)=tem_t(:,:,:,0)
!      return
!     endif

#ifdef bidon
      sum1=0.
      sum2=0.
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
       sum1=sum1+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)
       sum2=sum2+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*tem_t(i,j,k,1)
      enddo
      enddo
      enddo
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      if(par%rank==0)write(6,*)sum2glb,sum2glb/sum1glb
#endif



      cst_=0.5*cst_adv_hor/3. !11-09-17


!#ifdef bidon
      do k=1,kmax ; do j=-1,jmax+2 ; do i=-1,imax+2
       anyv3d(i,j,k,id_tbef_)=tem_t(i,j,k,before)
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,id_tdiv)=0. ! reset cumul temperature*divergence(courant) !26-01-17
      enddo ; enddo ; enddo

      do loop_=1,loopmaxts !*********>

! Flux Oi
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax
      do i=1,imax+1

      anyv3d(i,j,k,id_flx_)=0.5*( & !prime>


! Schema QUICK en dehors influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      anyv3d(i,j,k,id_hybcoefu)*( &
! Partie advective QUICK:
      -veldydz_u(i,j,k,1) &
!       *(anyv3d(i,j,k,id_tnow_)+anyv3d(i-1,j,k,id_tnow_)) &
! Si le nombre de courant est faible il n'y a pas de risque A maintenir la
! partie centree sur tem_t(now) ce qui permet de maintenir un leap-frog centrE
! classique meme en subdivisant le pas de temps, ce qui permet A priori de 
! maintenir les proprietes de conservation associeEs
        *(    tem_t(i,j,k,now)+tem_t(i-1,j,k,now)      &
#ifdef bidonref
       +2.*temref_u(i,j,k)                             &
          -temref_t(i,j,k) -temref_t(i-1,j,k)          &
#endif
                                                 )     &

! Partie diffusive QUICK remise A jour aux sous-pas de temps advectif:
       +cst_*0.5*(anyv3d(i,j,k,id_ofactort)+anyv3d(i,j+1,k,id_ofactort)) &

! efficacite densitaire
        *abs(-alp_t*(   tem_t(i,j,k,1)-tem_t(i-1,j,k,1))+alp_s*(   sal_t(i,j,k,1)-sal_t(i-1,j,k,1))) &
        /   (+alp_t*abs(tem_t(i,j,k,1)-tem_t(i-1,j,k,1))+alp_s*abs(sal_t(i,j,k,1)-sal_t(i-1,j,k,1))+small1) &

                                                               *( & !m°v°m>
             (veldydz_u(i,j,k,1)+abs(veldydz_u(i,j,k,1)))*( &           !pmxpmx>

               (   anyv3d(i  ,j,k,id_tbef_)   -anyv3d(i-1,j,k,id_tbef_)  &

! correction anyv3d(i  ,j,k,id_tbef_) par rapport niveau z(i-1,j,k)
      +( anyv3d(i  ,j,kp1,id_tbef_)- anyv3d(i  ,j,km1,id_tbef_))  &
      /(depth_t(i  ,j,kp1)         -depth_t(i  ,j,km1)         )  &
      *(depth_t(i-1,j,k  )         -depth_t(i  ,j,k  )         )  &

               ) &

              -(   anyv3d(i-1,j,k,id_tbef_)-anyv3d(i-2,j,k,id_tbef_)  &

! correction anyv3d(i-2,j,k,id_tbef_) par rapport niveau z(i-1,j,k)
      -( anyv3d(i-2,j,kp1,id_tbef_)- anyv3d(i-2,j,km1,id_tbef_))  &
      /(depth_t(i-2,j,kp1)         -depth_t(i-2,j,km1)         )  &
      *(depth_t(i-1,j,k  )         -depth_t(i-2,j,k  )         )  &

               ) &
                                           *mask_t(i-2,j,kmax)      ) & !pmxpmx>

            +(veldydz_u(i,j,k,1)-abs(veldydz_u(i,j,k,1)))*( &           !pmxpmx>

               (   anyv3d(i+1,j,k,id_tbef_)-anyv3d(i,j,k,id_tbef_)    &

! correction anyv3d(i+1,j,k,id_tbef_) par rapport niveau z(i  ,j,k)
      +( anyv3d(i+1,j,kp1,id_tbef_)- anyv3d(i+1,j,km1,id_tbef_))  &
      /(depth_t(i+1,j,kp1)         -depth_t(i+1,j,km1)         )  &
      *(depth_t(i  ,j,k  )         -depth_t(i+1,j,k  )         )  &

               )   &
                  *mask_t(i+1,j,kmax)                                 &

              -(   anyv3d(i,j,k,id_tbef_)-anyv3d(i-1,j,k,id_tbef_)    &

! correction anyv3d(i-1,j,k,id_tbef_) par rapport niveau z(i  ,j,k)
      -( anyv3d(i-1,j,kp1,id_tbef_)- anyv3d(i-1,j,km1,id_tbef_))  &
      /(depth_t(i-1,j,kp1)         -depth_t(i-1,j,km1)         )  &
      *(depth_t(i  ,j,k  )         -depth_t(i-1,j,k  )         )  &

                                                                  ))  & !pmxpmx>

                                                                )  & !m°v°m>()cst_
                                                                )  & !()anyv3d(i,j,k,id_hybcoefu)
! Schema UPWIND dans les zones sous influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      +(1.-anyv3d(i,j,k,id_hybcoefu))*( & !pmx>
           -veldydz_u(i,j,k,1)* ( & !ppp>
               anyv3d(i,j,k,id_tbef_)+anyv3d(i-1,j,k,id_tbef_)   &
                                ) & !ppp>
       +abs(veldydz_u(i,j,k,1))*( & !mmm>
               anyv3d(i,j,k,id_tbef_)-anyv3d(i-1,j,k,id_tbef_)   &
                                ) & !mmm>
                                      ) & !pmx>

       ) !prime>

!     -veldydz_u(i,j,k,1)*dti_lpsub*temref_u(i,j,k)

      enddo !k
      enddo !i
      enddo !j

! Partial Oi Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax


! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_tdiv)= & !26-01-17
       anyv3d(i,j,k,id_tdiv)  &
      +anyv3d(i,j,k,id_tbef_)*dti_lpsub*(veldydz_u(i+1,j,k,1)           &
                                        -veldydz_u(i  ,j,k,1))/dxdy_t(i,j)

!       anyv3d(i,j,k,id_tbef_)                       &
        anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !ooo> !25-04-17

        anyv3d(i,j,k,id_tbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i+1,j,k,id_flx_)                &
             -anyv3d(i  ,j,k,id_flx_)                &

         +(veldydz_u(i+1,j,k,1)                      &
          -veldydz_u(i  ,j,k,1))*anyv3d(i,j,k,id_tbef_) & !26-01-17


         )/dxdy_t(i,j) ) & !ooo>

         +(1-mask_t(i,j,k))*tem_t(i,j,k,before)*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

! Note: le terme anyv3d(i,j,k,id_tdiv) ne fait pas l'objet d'un traitement similaire pour le masque car
! si mask_t=0 alors veldydz_u(i:i+1,...)=0 (question: et les fleuves ?)

      enddo
      enddo
      enddo

      call vertmix_merged_levels_any(0,id_tbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_tbef_)
      call obc_int_anyv3d(id_tbef_,'zb')

! Flux Oj
      cst_=0.5*cst_adv_hor/3.
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax+1
      do i=1,imax

      anyv3d(i,j,k,id_fly_)=0.5*( & !prime>

! Schema QUICK en dehors influence des fleuves
      anyv3d(i,j,k,id_hybcoefv)*( &
! Partie advective QUICK:
      -veldxdz_v(i,j,k,1)  &
!       *(anyv3d(i,j,k,id_tnow_)+anyv3d(i,j-1,k,id_tnow_)) &
! Si le nombre de courant est faible il n'y a pas de risque A maintenir la
! partie centree sur tem_t(now) ce qui permet de maintenir un leap-frog centrE
! classique meme en subdivisant le pas de temps ce qui permet A priori de 
! maintenir les proprietes de conservation associeE
        *(    tem_t(i,j,k,now)+tem_t(i,j-1,k,now)  &
#ifdef bidonref
       +2.*temref_v(i,j,k)                         &
          -temref_t(i,j,k) -temref_t(i,j-1,k)      &
#endif
                                                 ) &

! Partie diffusive QUICK:
       +cst_*0.5*(anyv3d(i,j,k,id_ofactort)+anyv3d(i+1,j,k,id_ofactort)) &

! efficacite densitaire
        *abs(-alp_t*(   tem_t(i,j,k,1)-tem_t(i,j-1,k,1))+alp_s*(   sal_t(i,j,k,1)-sal_t(i,j-1,k,1))) &
        /   (+alp_t*abs(tem_t(i,j,k,1)-tem_t(i,j-1,k,1))+alp_s*abs(sal_t(i,j,k,1)-sal_t(i,j-1,k,1))+small1) &

                                                                  *(  & !m°v°m>

             (veldxdz_v(i,j,k,1)+abs(veldxdz_v(i,j,k,1)))*(  &       !ooo>

            (   anyv3d(i,j  ,k,id_tbef_)-anyv3d(i,j-1,k,id_tbef_)  &

! correction anyv3d(i,j  ,k,id_tbef_) par rapport niveau z(i,j-1,k)
      +( anyv3d(i,j  ,kp1,id_tbef_)- anyv3d(i,j  ,km1,id_tbef_))  &
      /(depth_t(i,j  ,kp1)         -depth_t(i,j  ,km1)         )  &
      *(depth_t(i,j-1,k  )         -depth_t(i,j  ,k  )         )  &

            ) &

           -(   anyv3d(i,j-1,k,id_tbef_)-anyv3d(i,j-2,k,id_tbef_)  &

! correction anyv3d(i,j-2,k,id_tbef_) par rapport niveau z(i,j-1,k)
      -( anyv3d(i,j-2,kp1,id_tbef_)- anyv3d(i,j-2,km1,id_tbef_))  &
      /(depth_t(i,j-2,kp1)         -depth_t(i,j-2,km1)         )  &
      *(depth_t(i,j-1,k  )         -depth_t(i,j-2,k  )         )  &

            ) &
                                        *mask_t(i,j-2,kmax)       )& !ooo>

            +(veldxdz_v(i,j,k,1)-abs(veldxdz_v(i,j,k,1)))*(  &       !ooo>

            (   anyv3d(i,j+1,k,id_tbef_)-anyv3d(i,j,k,id_tbef_)    &

! correction anyv3d(i,j+1,k,id_tbef_) par rapport niveau z(i,j  ,k)
      +( anyv3d(i,j+1,kp1,id_tbef_)- anyv3d(i,j+1,km1,id_tbef_))  &
      /(depth_t(i,j+1,kp1)         -depth_t(i,j+1,km1)         )  &
      *(depth_t(i,j  ,k  )         -depth_t(i,j+1,k  )         )  &

            )   &
               *mask_t(i,j+1,kmax)                                 &

           -(   anyv3d(i,j,k,id_tbef_)-anyv3d(i,j-1,k,id_tbef_)    &

! correction anyv3d(i,j-1,k,id_tbef_) par rapport niveau z(i,j  ,k)
      -( anyv3d(i,j-1,kp1,id_tbef_)- anyv3d(i,j-1,km1,id_tbef_))  &
      /(depth_t(i,j-1,kp1)         -depth_t(i,j-1,km1)         )  &
      *(depth_t(i,j  ,k  )         -depth_t(i,j-1,k  )         )  &

                                                               ))  & !ooo>

                                                                   )  & !m°v°m>()cst_
                                                                   )  & !()anyv3d(i,j,k,id_hybcoefv)
! Schema UPWIND dans les zones sous influence des fleuves:
      +(1.-anyv3d(i,j,k,id_hybcoefv))*( & !pmx>
          -veldxdz_v(i,j,k,1)* ( & !ppp>
              anyv3d(i,j,k,id_tbef_)+anyv3d(i,j-1,k,id_tbef_)   &
                               ) & !ppp>
      +abs(veldxdz_v(i,j,k,1))*( & !mmm>
              anyv3d(i,j,k,id_tbef_)-anyv3d(i,j-1,k,id_tbef_)   &
                               ) & !mmm>
                                      ) & !pmx>

       ) !prime>

!      -veldxdz_v(i,j,k,1)*dti_lpsub*temref_v(i,j,k)

      enddo !k
      enddo !i
      enddo !j

! Partial Horizontal Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_tdiv)= & !26-01-17
       anyv3d(i,j,k,id_tdiv)  &
      +anyv3d(i,j,k,id_tbef_)*dti_lpsub*(veldxdz_v(i,j+1,k,1)                 &
                                        -veldxdz_v(i,j  ,k,1))/dxdy_t(i,j)

        anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !ooo>
        anyv3d(i,j,k,id_tbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i,j+1,k,id_fly_)                &
             -anyv3d(i,j  ,k,id_fly_)                &

         +(veldxdz_v(i,j+1,k,1)                      &
          -veldxdz_v(i,j  ,k,1))*anyv3d(i,j,k,id_tbef_) & !26-01-17

         )/dxdy_t(i,j)  ) & !ooo>
 
         +(1-mask_t(i,j,k))*tem_t(i,j,k,before)*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

      enddo
      enddo
      enddo

      call vertmix_merged_levels_any(0,id_tbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_tbef_)
      call obc_int_anyv3d(id_tbef_,'zb')

!#endif
      enddo               !*********>

      do k=1,kmax ; do j=1,jmax  ; do i=1,imax  

           tem_t(i,j,k,2)                              &
       =  anyv3d(i,j,k,id_tbef_)                    !  &

!     -wetmask_t(i,j)*tem_t(i,j,k,now)*dti_lp*(        &

!                                              veldydz_u(i+1,j,k,1)   &
!                                             -veldydz_u(i  ,j,k,1)   &
!                                             +veldxdz_v(i,j+1,k,1)   &
!                                             -veldxdz_v(i,j  ,k,1)   &

!                                             )/(dz_t(i,j,k,before)*dxdy_t(i,j))

      enddo ; enddo ; enddo
      end subroutine advection_scal_upmx_subt


!..............................................................................

      subroutine advection_scal_upmx_subs  ! 05-01-16
      use module_principal
      use module_parallele
      implicit none
      double precision :: velef_,flagmsk_=1.,scalefac_,cst_
      integer :: t_=1        &
                ,id_flx_=1   &  ! x flux identifier
                ,id_fly_=1   &  ! y flux identifier
                ,id_sbef_=3  &  ! salinity    "before" identifier
!               ,id_sdiv=5  &  ! S*divergence(velocity)
                ,loop_
#ifdef synopsis
       subroutinetitle='advection_scal_upmx_subs'
       subroutinedescription=                                     &
          'Computes salperature and salinity advection diffusion' &
       //' horizontal fluxes using a "3 points" advection scheme' &
       //' hybrid QUICK-UPWIND'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Cas particulier du modele 1DV
!     if(flag_1dv==1) then !>>>>> !14-02-16
!      sal_t(:,:,:,2)=sal_t(:,:,:,0)
!      return
!     endif

      cst_=0.5*cst_adv_hor/3.


!#ifdef bidon
      do k=1,kmax ; do j=-1,jmax+2 ; do i=-1,imax+2
       anyv3d(i,j,k,id_sbef_)=sal_t(i,j,k,before) 
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,id_sdiv)=0. ! reset cumul salinity*divergence(courant)
      enddo ; enddo ; enddo


      do loop_=1,loopmaxts !*********>


! Flux Oi
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax
      do i=1,imax+1

      anyv3d(i,j,k,id_flx_)=0.5*( & !prime>


! Schema QUICK en dehors influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      anyv3d(i,j,k,id_hybcoefu)*( &
! Partie advective QUICK:
      -veldydz_u(i,j,k,1) &
!       *(anyv3d(i,j,k,id_tnow_)+anyv3d(i-1,j,k,id_tnow_)) &
! Si le nombre de courant est faible il n'y a pas de risque A maintenir la
! partie centree sur sal_t(now) ce qui permet de maintenir un leap-frog centrE
! classique meme en subdivisant le pas de salps, ce qui permet A priori de 
! maintenir les proprietes de conservation associeEs
        *(    sal_t(i,j,k,now)+sal_t(i-1,j,k,now)      &
#ifdef bidonref
       +2.*salref_u(i,j,k)                             &
          -salref_t(i,j,k) -salref_t(i-1,j,k)          &
#endif
                                                 )     &

! Partie diffusive QUICK remise A jour aux sous-pas de salps advectif:
       +cst_*0.5*(anyv3d(i,j,k,id_ofactors)+anyv3d(i,j+1,k,id_ofactors))  &

! Efficacite densitaire
        *abs(-alp_t*(   tem_t(i,j,k,1)-tem_t(i-1,j,k,1))+alp_s*(   sal_t(i,j,k,1)-sal_t(i-1,j,k,1))) &
        /   (+alp_t*abs(tem_t(i,j,k,1)-tem_t(i-1,j,k,1))+alp_s*abs(sal_t(i,j,k,1)-sal_t(i-1,j,k,1))+small1) &

                                                               *( & !m°v°m>
             (veldydz_u(i,j,k,1)+abs(veldydz_u(i,j,k,1)))*( &           !pmxpmx>

               (   anyv3d(i  ,j,k,id_sbef_)-anyv3d(i-1,j,k,id_sbef_)  &

! correction anyv3d(i  ,j,k,id_sbef_) par rapport niveau z(i-1,j,k)
      +( anyv3d(i  ,j,kp1,id_sbef_)- anyv3d(i  ,j,km1,id_sbef_))  &
      /(depth_t(i  ,j,kp1)         -depth_t(i  ,j,km1)         )  &
      *(depth_t(i-1,j,k  )         -depth_t(i  ,j,k  )         )  &

               ) &

              -(   anyv3d(i-1,j,k,id_sbef_)-anyv3d(i-2,j,k,id_sbef_)  &

! correction anyv3d(i-2,j,k,id_sbef_) par rapport niveau z(i-1,j,k)
      -( anyv3d(i-2,j,kp1,id_sbef_)- anyv3d(i-2,j,km1,id_sbef_))  &
      /(depth_t(i-2,j,kp1)         -depth_t(i-2,j,km1)         )  &
      *(depth_t(i-1,j,k  )         -depth_t(i-2,j,k  )         )  &

               ) &
                                           *mask_t(i-2,j,kmax)      ) & !pmxpmx>

            +(veldydz_u(i,j,k,1)-abs(veldydz_u(i,j,k,1)))*( &           !pmxpmx>

               (   anyv3d(i+1,j,k,id_sbef_)-anyv3d(i,j,k,id_sbef_)    &

! correction anyv3d(i+1,j,k,id_sbef_) par rapport niveau z(i  ,j,k)
      +( anyv3d(i+1,j,kp1,id_sbef_)- anyv3d(i+1,j,km1,id_sbef_))  &
      /(depth_t(i+1,j,kp1)         -depth_t(i+1,j,km1)         )  &
      *(depth_t(i  ,j,k  )         -depth_t(i+1,j,k  )         )  &

               )   &
                  *mask_t(i+1,j,kmax)                                 &

              -(   anyv3d(i,j,k,id_sbef_)-anyv3d(i-1,j,k,id_sbef_)    &

! correction anyv3d(i-1,j,k,id_sbef_) par rapport niveau z(i  ,j,k)
      -( anyv3d(i-1,j,kp1,id_sbef_)- anyv3d(i-1,j,km1,id_sbef_))  &
      /(depth_t(i-1,j,kp1)         -depth_t(i-1,j,km1)         )  &
      *(depth_t(i  ,j,k  )         -depth_t(i-1,j,k  )         )  &

                                                                  ))  & !pmxpmx>

                                                                )  & !m°v°m>()cst_
                                                                )  & !()anyv3d(i,j,k,id_hybcoefu)
! Schema UPWIND dans les zones sous influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      +(1.-anyv3d(i,j,k,id_hybcoefu))*( & !pmx>
           -veldydz_u(i,j,k,1)* ( & !ppp>
               anyv3d(i,j,k,id_sbef_)+anyv3d(i-1,j,k,id_sbef_)   &
                                ) & !ppp>
       +abs(veldydz_u(i,j,k,1))*( & !mmm>
               anyv3d(i,j,k,id_sbef_)-anyv3d(i-1,j,k,id_sbef_)   &
                                ) & !mmm>
                                      ) & !pmx>

       ) !prime>


      enddo !k
      enddo !i
      enddo !j

! Partial Oi Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax


! Cumul du pivot obligatoirement avant advection partielle 06-02-17
! Note que pour economiser un peu de cpu la multiplication par wetmask_t(i,j)/dz_t(i,j,k,before),
! commune a toutes les operations, sera appliquee une fois pour toute a la fin...
       anyv3d(i,j,k,id_sdiv)= &
       anyv3d(i,j,k,id_sdiv)  &
      +anyv3d(i,j,k,id_sbef_)*dti_lpsub*(veldydz_u(i+1,j,k,1)           &
                                        -veldydz_u(i  ,j,k,1))/dxdy_t(i,j)

        anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !ooo>
        anyv3d(i,j,k,id_sbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i+1,j,k,id_flx_)                &
             -anyv3d(i  ,j,k,id_flx_)                &

         +(veldydz_u(i+1,j,k,1)                      &
          -veldydz_u(i  ,j,k,1))*anyv3d(i,j,k,id_sbef_) & !25-01-17

         )/dxdy_t(i,j)  ) & !ooo>

         +(1-mask_t(i,j,k))*sal_t(i,j,k,before)*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

      enddo
      enddo
      enddo

      call vertmix_merged_levels_any(0,id_sbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_sbef_)
      call obc_int_anyv3d(id_sbef_,'zb')


! Flux Oj
      cst_=0.5*cst_adv_hor/3.
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax+1
      do i=1,imax

      anyv3d(i,j,k,id_fly_)=0.5*( & !prime>

! Schema QUICK en dehors influence des fleuves
      anyv3d(i,j,k,id_hybcoefv)*( &
! Partie advective QUICK:
      -veldxdz_v(i,j,k,1)  &
!       *(anyv3d(i,j,k,id_tnow_)+anyv3d(i,j-1,k,id_tnow_)) &
! Si le nombre de courant est faible il n'y a pas de risque A maintenir la
! partie centree sur sal_t(now) ce qui permet de maintenir un leap-frog centrE
! classique meme en subdivisant le pas de salps ce qui permet A priori de 
! maintenir les proprietes de conservation associeE
        *(    sal_t(i,j,k,now)+sal_t(i,j-1,k,now)  &
#ifdef bidonref
       +2.*salref_v(i,j,k)                         &
          -salref_t(i,j,k) -salref_t(i,j-1,k)      &
#endif
                                                 ) &

! Partie diffusive QUICK:
       +cst_*0.5*(anyv3d(i,j,k,id_ofactors)+anyv3d(i+1,j,k,id_ofactors))  &

! Efficacite densitaire
        *abs(-alp_t*(   tem_t(i,j,k,1)-tem_t(i,j-1,k,1))+alp_s*(   sal_t(i,j,k,1)-sal_t(i,j-1,k,1))) &
        /   (+alp_t*abs(tem_t(i,j,k,1)-tem_t(i,j-1,k,1))+alp_s*abs(sal_t(i,j,k,1)-sal_t(i,j-1,k,1))+small1) &

                                                                  *( & !m°v°m>

             (veldxdz_v(i,j,k,1)+abs(veldxdz_v(i,j,k,1)))*(  &       !ooo>

            (   anyv3d(i,j  ,k,id_sbef_)-anyv3d(i,j-1,k,id_sbef_)  &

! correction anyv3d(i,j  ,k,id_sbef_) par rapport niveau z(i,j-1,k)
      +( anyv3d(i,j  ,kp1,id_sbef_)- anyv3d(i,j  ,km1,id_sbef_))  &
      /(depth_t(i,j  ,kp1)         -depth_t(i,j  ,km1)         )  &
      *(depth_t(i,j-1,k  )         -depth_t(i,j  ,k  )         )  &

            ) &

           -(   anyv3d(i,j-1,k,id_sbef_)-anyv3d(i,j-2,k,id_sbef_)  &

! correction anyv3d(i,j-2,k,id_sbef_) par rapport niveau z(i,j-1,k)
      -( anyv3d(i,j-2,kp1,id_sbef_)- anyv3d(i,j-2,km1,id_sbef_))  &
      /(depth_t(i,j-2,kp1)         -depth_t(i,j-2,km1)         )  &
      *(depth_t(i,j-1,k  )         -depth_t(i,j-2,k  )         )  &

            ) &
                                        *mask_t(i,j-2,kmax)       )& !ooo>

            +(veldxdz_v(i,j,k,1)-abs(veldxdz_v(i,j,k,1)))*(  &       !ooo>

            (   anyv3d(i,j+1,k,id_sbef_)-anyv3d(i,j,k,id_sbef_)    &

! correction anyv3d(i,j+1,k,id_sbef_) par rapport niveau z(i,j  ,k)
      +( anyv3d(i,j+1,kp1,id_sbef_)- anyv3d(i,j+1,km1,id_sbef_))  &
      /(depth_t(i,j+1,kp1)         -depth_t(i,j+1,km1)         )  &
      *(depth_t(i,j  ,k  )         -depth_t(i,j+1,k  )         )  &

            )   &
               *mask_t(i,j+1,kmax)                                 &

           -(   anyv3d(i,j,k,id_sbef_)-anyv3d(i,j-1,k,id_sbef_)    &

! correction anyv3d(i,j-1,k,id_sbef_) par rapport niveau z(i,j  ,k)
      -( anyv3d(i,j-1,kp1,id_sbef_)- anyv3d(i,j-1,km1,id_sbef_))  &
      /(depth_t(i,j-1,kp1)         -depth_t(i,j-1,km1)         )  &
      *(depth_t(i,j  ,k  )         -depth_t(i,j-1,k  )         )  &

                                                               ))  & !ooo>

                                                                   )  & !m°v°m>()cst_
                                                                   )  & !()anyv3d(i,j,k,id_hybcoefv)
! Schema UPWIND dans les zones sous influence des fleuves:
      +(1.-anyv3d(i,j,k,id_hybcoefv))*( & !pmx>
          -veldxdz_v(i,j,k,1)* ( & !ppp>
              anyv3d(i,j,k,id_sbef_)+anyv3d(i,j-1,k,id_sbef_)   &
                               ) & !ppp>
      +abs(veldxdz_v(i,j,k,1))*( & !mmm>
              anyv3d(i,j,k,id_sbef_)-anyv3d(i,j-1,k,id_sbef_)   &
                               ) & !mmm>
                                      ) & !pmx>

       ) !prime>

      enddo !k
      enddo !i
      enddo !j

! Partial Horizontal Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_sdiv)= &
       anyv3d(i,j,k,id_sdiv)  &
      +anyv3d(i,j,k,id_sbef_)*dti_lpsub*(veldxdz_v(i,j+1,k,1)   &
                                        -veldxdz_v(i,j  ,k,1))/dxdy_t(i,j)

        anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !ooo>
        anyv3d(i,j,k,id_sbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i,j+1,k,id_fly_)                &
             -anyv3d(i,j  ,k,id_fly_)                &

         +(veldxdz_v(i,j+1,k,1)                      &
          -veldxdz_v(i,j  ,k,1))*anyv3d(i,j,k,id_sbef_) & !25-01-17

         )/dxdy_t(i,j)  ) & !ooo>

         +(1-mask_t(i,j,k))*sal_t(i,j,k,before)*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

      enddo
      enddo
      enddo

      call vertmix_merged_levels_any(0,id_sbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_sbef_)
      call obc_int_anyv3d(id_sbef_,'zb')

!#endif
      enddo               !*********>

      do k=1,kmax ; do j=1,jmax  ; do i=1,imax  

           sal_t(i,j,k,2)                              &
       =  anyv3d(i,j,k,id_sbef_)                    !  &

!     -wetmask_t(i,j)*sal_t(i,j,k,now)*dti_lp*(        &

!                                              veldydz_u(i+1,j,k,1)   &
!                                             -veldydz_u(i  ,j,k,1)   &
!                                             +veldxdz_v(i,j+1,k,1)   &
!                                             -veldxdz_v(i,j  ,k,1)   &

!                                             )/(dz_t(i,j,k,before)*dxdy_t(i,j))


      enddo ; enddo ; enddo


      end subroutine advection_scal_upmx_subs
#endif
!................................................................

      subroutine advection_scal_hybcoef
      use module_principal
      use module_parallele
      implicit none

      do k=1,kmax

       do j=0,jmax+1
       do i=0,imax+1
!       xy_t(i,j,1)=1./dxdy_t(i,j)     &
!                       /dz_t(i,j,k,2)
        xy_t(i,j,1)=1./anyv3d(i,j,k,id_dxdydz)
       enddo
       enddo

       do j=1,jmax
       do i=1,imax+1
        anyv3d(i,j,k,id_hybcoefu)=max(1.-(      & !pmx>
              veldydz_u(i,j,k,1)*dti_lp*max(xy_t(i,j,1),xy_t(i-1,j,1)) &
                                         )**2   & !pmx>
                                  ,0.)          &
       *wetmask_u(i,j)*0.5*(upwindriver_t(i,j)+upwindriver_t(i-1,j)) & 
                         *0.5*(upwzone0_t(i,j,k) +upwzone0_t(i-1,j,k)) !21-05-20
! Note: en situation standard 0.5*(upwzone0_t(i,j,k)+upwzone0_t(i-1,j,k))=1
!       et devient 0 si schema 100% upwind

       enddo
       enddo

       do j=1,jmax+1
       do i=1,imax
        anyv3d(i,j,k,id_hybcoefv)=max(1.-(      & !pmx>
              veldxdz_v(i,j,k,1)*dti_lp*max(xy_t(i,j,1),xy_t(i,j-1,1)) &
                                         )**2   & !pmx>
                                  ,0.)          &
       *wetmask_v(i,j)*0.5*(upwindriver_t(i,j)+upwindriver_t(i,j-1))   & 
                         *0.5*(upwzone0_t(i,j,k) +upwzone0_t(i,j-1,k)) !21-05-20

       enddo
       enddo

      enddo
!     anyv3d(:,:,:,id_hybcoefu)=1. ! 100% up3
!     anyv3d(:,:,:,id_hybcoefv)=1. ! 100% up3
!     anyv3d(:,:,:,id_hybcoefu)=0. ! 100% up2
!     anyv3d(:,:,:,id_hybcoefv)=0. ! 100% up2

! Flaguer les zones d'instabilite de stratification le long du niveau de fond
! Le gradient est instable si drho/dz>0 (ou, histoire d'economiser une division, drho*dz>0)
! Si l'instabilite est declaree, le flux du point u voisin "en dessous" est upwind.
! L'algo regarde donc le signe du produit dz*drho au point voisin "en dessous".
! Il faut donc detecter le signe de la pente pour trouver le voisin du dessous.
! Pour eviter echange mpi on suppose que la pente d(zb)/di a le meme signe aux deux points encadrant le point i
! On regarde le gradient "en dessous", soit entre i-2 et i-1 si pente d(zb)/di > 0
!                                        et entre i+1 et i   si pente d(zb)/di < 0
! on a choisit de considerer depth_w(:,:,1) (parce que immobile, mais deph_t(:,:,1) aurait pu faire)

!     return

! Assurer la continuite mpi de rhp(:,:,1) en donnant des valeurs en i=imax+1, imax+2 etc.....
       call obc_scal_bottom_rhp('zb')  ! echange z1 et z2

       do j=1,jmax ; do i=1,imax+1

        if(depth_w(i,j,1)-depth_w(i-1,j,1)>0) then !m°v°m>
      
! si dz/di>0 alors dz*drho>0 si drho>0
           if(mask_t(i-2,j,1)*(rhp_t(i-1,j,1)-rhp_t(i-2,j,1))>0.) then !ooo> !21-05-18
               do k=1,kmerged_u(i,j)
                         anyv3d(i,j,k,id_hybcoefu)=0.
               enddo
           endif                                                       !ooo>

        else                                       !m°v°m>

! si dz/di<0 alors dz*drho>0 si drho<0
           if(mask_t(i+1,j,1)*(rhp_t(i+1,j,1)-rhp_t(i  ,j,1))<0.) then !ooo>
               do k=1,kmerged_u(i,j)
                         anyv3d(i,j,k,id_hybcoefu)=0.
               enddo
           endif                                                       !ooo>

        endif                                      !m°v°m>

       enddo      ; enddo

       do j=1,jmax+1 ; do i=1,imax

        if(depth_w(i,j,1)-depth_w(i,j-1,1)>0) then !m°v°m>
      
! si dz/di>0 alors dz*drho>0 si drho>0
           if(mask_t(i,j-2,1)*(rhp_t(i,j-1,1)-rhp_t(i,j-2,1))>0.) then !ooo>
               do k=1,kmerged_v(i,j)
                         anyv3d(i,j,k,id_hybcoefv)=0.
               enddo
           endif                                                       !ooo>

        else                                       !m°v°m>

! si dz/di<0 alors dz*drho>0 si drho<0
           if(mask_t(i,j+1,1)*(rhp_t(i,j+1,1)-rhp_t(i  ,j,1))<0.) then !ooo>
               do k=1,kmerged_v(i,j)
                         anyv3d(i,j,k,id_hybcoefv)=0.
               enddo
           endif                                     !ooo>

        endif                                      !m°v°m>

       enddo      ; enddo

      end subroutine advection_scal_hybcoef

!................................................................
      subroutine advection_scal_z_laxwen !12-11-16
      use module_principal
      use module_parallele
      implicit none
      double precision :: grdamp_=2. ! lower gradient amplification
      integer :: id_cn_=0           &! current number "limited"
                ,id_we_=1           &! Explicit limited vertical velocity
                ,id_wi_=2           &! residual implicit vertical velocity
                ,loop_,loopmax_     & 
                ,looplimit_=1        ! nombre d'iteration maximum la procedure iterative
#ifdef synopsis
       subroutinetitle='advection_scal_z_laxwen'
       subroutinedescription='advection_scal_z_laxwen'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      stop 'Ne plus passer par laxwen depuis 22-11-16'

      xy_t(:,:,1)=small2 ! reset small2 evite que le nombre d'iteration soit nulle
      do k=2,kmax ; do j=1,jmax ; do i=1,imax
! x0=Dz/Dt
       x0=min(dz_t(i,j,k,0),dz_t(i,j,k-1,0))/dti_lp
! Omega bornee par + ou - N fois dz/dt ici 2 fois
       anyv3d(i,j,k,id_we_)=max(min(omega_w(i,j,k,1),looplimit_*x0),-looplimit_*x0)*wetmask_t(i,j)
! max over z of the decimal current number:
       xy_t(i,j,1)=max(xy_t(i,j,1),abs(anyv3d(i,j,k,id_we_))/x0)
      enddo       ; enddo       ; enddo
! Au fond et A la surface:
      do j=1,jmax ; do i=1,imax
       anyv3d(i,j,1     ,id_we_)=omega_w(i,j,1     ,1)
       anyv3d(i,j,kmax+1,id_we_)=omega_w(i,j,kmax+1,1) ! A la surface w explicite pour flux de surface
      enddo       ; enddo

!....................................................................
! Decommenter ces lignes pour avoir la valeur max sur tout le domaine
!     x1=0.
!     do j=1,jmax ; do i=1,imax
!     x1=max(x1,xy_t(i,j,1))
!     enddo ; enddo 
!     call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_max,par%comm2d,ierr)
!     if(par%rank==1)write(66,*)x2
!....................................................................

!     write(6,*)xy_t(imax/2,jmax/2,1)

! Vitesse verticale residuelle implicite:
      do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,id_wi_)=omega_w(i,j,k,1)-anyv3d(i,j,k,id_we_)
!      if(i==imax/2.and.j==jmax/2)write(6,*)k,real(omega_w(i,j,k,1)),real(anyv3d(i,j,k,id_we_)),real(anyv3d(i,j,k,id_wi_))
      enddo       ; enddo       ; enddo

      do j=1,jmax ; do i=1,imax ! Boucle i,j No1

! ceiling de xy_t(i,j,1) nombre de sous-iterations advectives
       loopmax_=min( ceiling(xy_t(i,j,1)) , looplimit_ )
       dti_lpsub=dti_lp/loopmax_

! Cas general
      do k=kmin_w(i,j)+2,kmax                         !22-12-14

! Current number with additional limititations
       anyv3d(i,j,k,id_cn_)=1.-(    &!ooooo>
! Robustness regarding convection: 
! https://docs.google.com/document/d/1WsIhfKI1GP269wglgIM6YUKXdVrONILKupKP9Gui9hI/edit
!                 (1.-min(max(                                       &
!             -1.+min((  rhp_t(i,j,k  )-  rhp_t(i,j,k-1))            &
!                    /(depth_t(i,j,k  )-depth_t(i,j,k-1)),-small1)   &
!        /min(grdamp_*(  rhp_t(i,j,k-1)-  rhp_t(i,j,k-2))            &
!                    /(depth_t(i,j,k-1)-depth_t(i,j,k-2)),-small1)   &
!                 ,0.),1.))                                        &
! Robustness regarding rivers and dry zones:
                       upwindriver_t(i,j)   & ! river zone influence
                      *wetmask_t(i,j)       & ! wet/dry zone
! Robustness regarding vertical current number:
        *max(1.-abs(anyv3d(i,j,k,id_we_)*dti_lpsub               &
                 /min(dz_t(i,j,k,0),dz_t(i,j,k-1,0))),0.)        &  
                              ) !ooooo>


!        anyv3d(i,j,k,id_cn_)=0. ! 100%C2
!        anyv3d(i,j,k,id_cn_)=1. ! 100%UP2
       enddo

! Les cas particuliers:

! Cas particulier k=kmin_w(i,j)+1 qui ne connait pas la stratif entre
! kmin_w(i,j) et kmin_w(i,j)-1 et du coup on essaie quand meme de detecter
! la presence d'un bi-couche mais en utilisant la couche du dessus !09-07-15
      k=kmin_w(i,j)+1
! Current number with additional limititations
      anyv3d(i,j,k,id_cn_)=1.-(    &!ooooo>
! Robustness regarding convection: 
!                 (1.-min(max(                                       &
!             -1.+min((  rhp_t(i,j,k  )-  rhp_t(i,j,k-1))            &
!                    /(depth_t(i,j,k  )-depth_t(i,j,k-1)),-small1)   &
!        /min(grdamp_*(  rhp_t(i,j,k+1)-  rhp_t(i,j,k  ))            &
!                    /(depth_t(i,j,k+1)-depth_t(i,j,k  )),-small1)   &
!                 ,0.),1.))                                        &
! Robustness regarding rivers and dry zones:
                       upwindriver_t(i,j)   & ! river zone influence
                      *wetmask_t(i,j)       & ! wet/dry zone
! Robustness regarding vertical current number:
        *max(1.-abs(anyv3d(i,j,k,id_we_)*dti_lpsub               &
                 /min(dz_t(i,j,k,0),dz_t(i,j,k-1,0))),0.)        &  
                              ) !ooooo>


! Cas particuliers k=1:kmin_w(i,j) et k=kmax:
       anyv3d(i,j,1:kmin_w(i,j),id_cn_)=0.
       anyv3d(i,j,kmax+1       ,id_cn_)=0. ! A la surface schema 100% explicite pour
                                           ! garantir la coherence du flux d'eau douce
                                           ! qui, dans vertmix_sal, est calcule avec 
                                           ! avec sal_t(:,:,kmax,1). (Conservation bilan 
                                           ! de sel) !12-04-15

! Flux centres frozen:
       do k=2,kmax

        anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)  &
                  *( ( dz_t(i,j,k  ,1)*tem_t(i,j,k  ,1)             &
                      +dz_t(i,j,k-1,1)*tem_t(i,j,k-1,1))            &
                    /( dz_t(i,j,k  ,1) +dz_t(i,j,k-1,1)) )      

        anyv1d(k,2)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)  &
                  *( ( dz_t(i,j,k  ,1)*sal_t(i,j,k  ,1)             &
                      +dz_t(i,j,k-1,1)*sal_t(i,j,k-1,1))            &
                    /( dz_t(i,j,k  ,1) +dz_t(i,j,k-1,1)) )      

       enddo
       anyv1d(1     ,1:4)=0.
       anyv1d(kmax+1,3:4)=0.
       k=kmax+1
       anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)*tem_t(i,j,k-1,1)
       anyv1d(k,2)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)*sal_t(i,j,k-1,1)

! Partie Iterative:
       k1=1
       do loop_=1,loopmax_ ! boucle iterative >>>

       do k=2,kmax !nexon>

! Flux upwind:
        anyv1d(k,3)=-0.5*anyv3d(i,j,k,id_cn_)*( & !pmxpmx>

        (anyv3d(i,j,k,id_we_)+abs(anyv3d(i,j,k,id_we_)))*tem_t(i,j,k-1,2)&
       +(anyv3d(i,j,k,id_we_)-abs(anyv3d(i,j,k,id_we_)))*tem_t(i,j,k  ,2)&

                                              )   !pmxpmx>

        anyv1d(k,4)=-0.5*anyv3d(i,j,k,id_cn_)*( & !pmxpmx>

        (anyv3d(i,j,k,id_we_)+abs(anyv3d(i,j,k,id_we_)))*sal_t(i,j,k-1,2)&
       +(anyv3d(i,j,k,id_we_)-abs(anyv3d(i,j,k,id_we_)))*sal_t(i,j,k  ,2)&

                                              )   !pmxpmx>

       enddo       !nexon>

       if(loop_==loopmax_)k1=1-loopmax_ ! Raccourcir etape finalisation
       do k=1,kmax

         tem_t(i,j,k,2)=tem_t(i,j,k,2)+dti_lpsub*( &

                              anyv1d(k+1,3)+anyv1d(k+1,1) &
                             -anyv1d(k  ,3)-anyv1d(k  ,1) &
         
            +k1*tem_t(i,j,k,now)*( anyv3d(i,j,k+1,id_we_)    &
                                  -anyv3d(i,j,k  ,id_we_))   &

                                                 )/dz_t(i,j,k,before)

         sal_t(i,j,k,2)=sal_t(i,j,k,2)+dti_lpsub*( & 

                              anyv1d(k+1,4)+anyv1d(k+1,2) &
                             -anyv1d(k  ,4)-anyv1d(k  ,2) &
         
            +k1*sal_t(i,j,k,now)*( anyv3d(i,j,k+1,id_we_)    &
                                  -anyv3d(i,j,k  ,id_we_))   &

                                                 )/dz_t(i,j,k,before)

       enddo


       enddo               ! boucle iterative >>>


      enddo ; enddo ! Boucle i,j No1
! 


!...............
! La vitesse residuelle implicite est omega_w moins omega_explicite 


! Implicit "up2" (upwind 2 points) part of the vertical advection: !20-08-14
      do k=1,kmax-1 !kmax !09-04-15
      do j=1,jmax
      do i=1,imax

!     tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
!      anyv3d(i,j,k+1,id_cn_)*( anyv3d(i,j,k+1,id_wi_)+abs(anyv3d(i,j,k+1,id_wi_))) &
!     +anyv3d(i,j,k  ,id_cn_)*(-anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

!     tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
!      anyv3d(i,j,k  ,id_cn_)*(-anyv3d(i,j,k,id_wi_)-abs(anyv3d(i,j,k,id_wi_)))

!     tridia_in(i,j,k,3)=tridia_in(i,j,k,3)+0.5*dti_lp*                &
!      anyv3d(i,j,k+1,id_cn_)*( anyv3d(i,j,k+1,id_wi_)-abs(anyv3d(i,j,k+1,id_wi_)))

      tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
       ( anyv3d(i,j,k+1,id_wi_)+abs(anyv3d(i,j,k+1,id_wi_))) &
      +(-anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

      tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
       (-anyv3d(i,j,k,id_wi_)-abs(anyv3d(i,j,k,id_wi_)))

      tridia_in(i,j,k,3)=tridia_in(i,j,k,3)+0.5*dti_lp*                &
       ( anyv3d(i,j,k+1,id_wi_)-abs(anyv3d(i,j,k+1,id_wi_)))

      enddo
      enddo
      enddo

! A la surface on suppose une condition de gradient nul (ce qui signifie que !09-04-15
! la variation de T,S due aux flux de surface passera par le membre de droite
! dans vertmix.F90) sans modifier le protocole habituel.
! Pour obtenir les lignes suivantes, on part des lignes precedentes (cas "dans la couche")
! A1*T(k-1)+A2*T(k)+A3*T(k+1) et on remplace T(k+1) par T(k) pour avoir la condition de
! gradient nul et deduire les nouveaux coefficients: A1=A1 ; A2=A2+A3 ; A3=0. 
! Notons qu'on a force le schema a etre implicite, autrement dit anyv3d(i,j,k+1,0)=0
! et disparait du calcul. On obtient:
      k=kmax
      do j=1,jmax ; do i=1,imax

!     tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
!      anyv3d(i,j,k  ,id_cn_)*(-anyv3d(i,j,k,id_wi_)-abs(anyv3d(i,j,k,id_wi_)))

!     tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
!     +anyv3d(i,j,k+1,id_cn_)*2.*anyv3d(i,j,k+1,id_wi_)                &
!     +anyv3d(i,j,k  ,id_cn_)*( -anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

      tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
       (-anyv3d(i,j,k,id_wi_)-abs(anyv3d(i,j,k,id_wi_)))

      tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
      +2.*anyv3d(i,j,k+1,id_wi_)                &
      +( -anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

      enddo ; enddo

      end subroutine advection_scal_z_laxwen

!................................................................

      subroutine advection_scal_z_c2limited
      use module_principal
      use module_parallele
      implicit none
!     double precision :: grdamp_=2. ! lower gradient amplification
      integer :: id_cn_=0           &! current number "limited"
                ,id_we_=1           &! Explicit limited vertical velocity
                ,id_wi_=2           &! residual implicit vertical velocity
!               ,id_tdiv=4         &! cumulated T*div(current) !26-01-17
!               ,id_sdiv=5         &! cumulated S*div(current)
                ,loop_,loopmax_     &
                ,looplimit_=1        ! nombre d'iteration maximum la procedure iterative
#ifdef synopsis
       subroutinetitle='advection_scal_z_quickest'
       subroutinedescription='advection_scal_z_quickest'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!     x1=0. ! Pas obligatoire seulement pour diag du maximum
      xy_t(:,:,1)=small2 ! reset small2 evite que le nombre d'iteration soit nulle
      do k=2,kmax ; do j=1,jmax ; do i=1,imax

!   x0=Dz/Dt
        x0=min(dz_t(i,j,k,0),dz_t(i,j,k-1,0))/dti_lp
!   Omega bornee par + ou - N fois dz/dt ici 2 fois
        anyv3d(i,j,k,id_we_)=max(min(omega_w(i,j,k,1),looplimit_*x0),-looplimit_*x0)*wetmask_t(i,j) &
                            *wetmask_wi_t(i,j) !19-06-19
!   max over z of the decimal current number:
        xy_t(i,j,1)=max(xy_t(i,j,1),abs(anyv3d(i,j,k,id_we_))/x0)
!   x1=max(x1,xy_t(i,j,1)) ! Pas obligatoire seulement pour diag du maximum

      enddo  ; enddo ; enddo

! Fond w explicite = 0
      do j=1,jmax ; do i=1,imax
        anyv3d(i,j,1,id_we_)=0.
      enddo       ; enddo

! Couche merged (dans laquelle w est 100% implicite, donc w explicite = 0)
      if(flag_merged_levels==1) then !pmxpmx> !18-02-17
       do j=1,jmax ; do i=1,imax
        do k=2,kmerged_t(i,j)
         anyv3d(i,j,k,id_we_)=0.
        enddo
       enddo       ; enddo
      endif                          !pmxpmx> !18-02-17

! Surface:
      do j=1,jmax ; do i=1,imax
!      anyv3d(i,j,kmax+1,id_we_)=omega_w(i,j,kmax+1,1) ! A la surface w explicite pour flux de surface
       anyv3d(i,j,kmax+1,id_we_)=0. ! 25-02-19
      enddo       ; enddo

!....................................................................
! Decommenter ces lignes pour avoir la valeur max sur tout le domaine !30-12-16
!     call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_max,par%comm2d,ierr)
!     if(par%rank==0) then !>>>>>
!       open(unit=3,file='tmp/dti_w',position='append')
!        write(3,*)real(elapsedtime_now/86400.),real(x2)
!       close(3)
!     endif                !>>>>> 
!....................................................................

! Vitesse verticale residuelle implicite:
      do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,id_wi_)=omega_w(i,j,k,1)-anyv3d(i,j,k,id_we_)
      enddo       ; enddo       ; enddo

! Dans la couche merged, profil de vitesse verticale implicite de type "gradient de flux constant"
!     if(flag_merged_levels==1) then !pmxpmx> !18-02-17
!      do j=1,jmax ; do i=1,imax
!       do k=2,kmerged_t(i,j)
!        anyv3d(i,j,k,id_wi_)=anyv3d(i,j,kmerged_t(i,j)+1,id_wi_)       &
!                          *(depth_w(i,j,k)            -depth_w(i,j,1)) &
!                          /(depth_w(i,j,kmerged_t(i,j)+1)-depth_w(i,j,1))  
!       enddo
!      enddo       ; enddo
!     endif                          !pmxpmx> !18-02-17

      do j=1,jmax ; do i=1,imax ! Boucle i,j No1

! ceiling de xy_t(i,j,1) nombre de sous-iterations advectives
       loopmax_=min( ceiling(xy_t(i,j,1)) , looplimit_ )
       dti_lpsub=dti_lp/loopmax_

! Limiteurs C1,C2 ou C3
! https://docs.google.com/document/d/1WsIhfKI1GP269wglgIM6YUKXdVrONILKupKP9Gui9hI/edit
      do k=2,kmax-1
!      if(rhp_t(i,j,k  )-rhp_t(i,j,k-1)>=0) then
       if(rhp_t(i,j,k  )-rhp_t(i,j,k-1)>0.01*(rhp_t(i,j,k+1)-rhp_t(i,j,k))) then !m°v°m>
! if(rhp_t(i,j,k)*(depth_t(i,j,k+1)-depth_t(i,j,k-1)-dti_fw*0.5*(omega_w(i,j,k 
!!! suggestion de test: if((rhp(k)-rhp(k-1))/(small1+abs(rhp(kmax)-rhp(kmin)))>-1e-n afin de remplacer le test
!!! >=0 par une valeur legerement negative anticipative confrontee au gradient normalisE par rhp(surf)-rhp(bot)
            anyv1d(k,1)=0.
       else                                                                      !m°v°m>
            anyv1d(k,1)=1.
       endif

! Proposition de fonction continue: anyv1d(k,1)=0. (100 %up) si (rhp(k)-rhp(k-1))/(rhp(k+1)-rhp(k))=0.
!                                   anyv1d(k,1)=1. (100 %c2) si (rhp(k)-rhp(k-1))/(rhp(k+1)-rhp(k))=0.1
! anyv1d(k,1)=min(max( (rhp_t(i,j,k)-rhp_t(i,j,k-1))/min(0.1*(rhp_t(i,j,k+1)-rhp_t(i,j,k)),small1),0.),1.)
! Autrement dit si rhp(k)-rhp(k-1) est nul alors 100% upwind
!                  rhp(k)-rhp(k-1) est 10 fois plus petite que 
!                  la pente du point au dessus, alors 100% c2 et
!                  commence la transition continue vers upwind



      enddo                                                                      !m°v°m>
!     anyv1d(kmax+1,1)=1. ! inutile si C1
!     anyv1d(kmax  ,1)=1. ! inutile si C1
      anyv1d(     1,1)=1.

      do k=2,kmax  ! do k=kmin_w(i,j)+1,kmax                         !22-12-14


! Current number with additional limititations
       anyv3d(i,j,k,id_cn_)=( &         !power>
                              1.-( &!ooooo>
! Robustness regarding convection: 
! https://docs.google.com/document/d/1WsIhfKI1GP269wglgIM6YUKXdVrONILKupKP9Gui9hI/edit
!                      anyv1d(k+1,1)*         & ! limiteur C3 cette ligne et les 2 lignes suivantes
!                      anyv1d(k  ,1)*         & ! limiteur C2 cette ligne et la ligne suivante
                       anyv1d(k-1,1)*         & ! limiteur C1 cette ligne seulement
! Robustness regarding rivers and dry zones:
                       upwindriver_t(i,j)   & ! river zone influence
                         *upwzone0_t(i,j,k) & ! additional requests for upwind. Note upwzone0_t=1 in standard conditions
                      *wetmask_t(i,j)       & ! wet/dry zone
! Robustness regarding vertical current number:
        *max(1.-abs(anyv3d(i,j,k,id_we_)*dti_lpsub               &
                 /min(dz_t(i,j,k,0),dz_t(i,j,k-1,0))),0.)        &  
                                 ) &!ooooo>                        
                            )**cn_power !power !26-05-18

!        anyv3d(i,j,k,id_cn_)=0. ! 100%C2
!        anyv3d(i,j,k,id_cn_)=1. ! 100%UP2

       enddo


! Cas particuliers k=1:kmin_w(i,j) et k=kmax:
!      anyv3d(i,j,1:kmin_w(i,j),id_cn_)=0.
       anyv3d(i,j,1            ,id_cn_)=0.
!      anyv3d(i,j,kmax+1       ,id_cn_)=0. ! A la surface schema 100% explicite pour
!                                          ! garantir la coherence du flux d'eau douce
!                                          ! qui, dans vertmix_sal, est calcule avec 
!                                          ! avec sal_t(:,:,kmax,1). (Conservation bilan 
!                                          ! de sel) !12-04-15
       anyv3d(i,j,kmax+1       ,id_cn_)=1. ! A la surface schema upwind mais noter
! surtout que le choix nc=0 ou 1 n'a pas d'importance si la vitesse verticale est
! 100% implicite en kmax+1
! Voir mes commentaires à propos du passage en implicite en kmax+1:
! https://docs.google.com/document/d/1fIzAG9mo_zvtTdVauu7yw-VkIqFRldqq4tF2b6ZOnzI/edit


! Flux centres frozen:
       do k=2,kmax

        anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)  &
                  *( ( dz_t(i,j,k  ,1)*tem_t(i,j,k  ,1)             &
                      +dz_t(i,j,k-1,1)*tem_t(i,j,k-1,1))            &
                    /( dz_t(i,j,k  ,1) +dz_t(i,j,k-1,1)) )      

        anyv1d(k,2)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)  &
                  *( ( dz_t(i,j,k  ,1)*sal_t(i,j,k  ,1)             &
                      +dz_t(i,j,k-1,1)*sal_t(i,j,k-1,1))            &
                    /( dz_t(i,j,k  ,1) +dz_t(i,j,k-1,1)) )      

       enddo
       anyv1d(1     ,1:4)=0.
       anyv1d(kmax+1,3:4)=0.
       k=kmax+1
       anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)*tem_t(i,j,k-1,1)
       anyv1d(k,2)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)*sal_t(i,j,k-1,1)

! Partie Iterative:
       do loop_=1,loopmax_ ! boucle iterative >>>


! A priori lignes suivantes inutiles si schema implicite en kmax+1
       tem_t(i,j,kmax+1,2)=tem_t(i,j,kmax,2)
       tem_t(i,j,0     ,2)=tem_t(i,j,1   ,2)
       sal_t(i,j,kmax+1,2)=sal_t(i,j,kmax,2)
       sal_t(i,j,0     ,2)=sal_t(i,j,1   ,2)
       do k=2,kmax !nexon>

        x1=0.5*(anyv3d(i,j,k,id_we_)+abs(anyv3d(i,j,k,id_we_)))
        x2=0.5*(anyv3d(i,j,k,id_we_)-abs(anyv3d(i,j,k,id_we_)))

! Flux upwind:
        anyv1d(k,3)=  &
          -anyv3d(i,j,k,id_cn_) *(x1*tem_t(i,j,k-1,2)+x2*tem_t(i,j,k  ,2)) !&
!      +(1-anyv3d(i,j,k,id_cn_))*( &
!              x1*(tem_t(i,j,k  ,2)-2.*tem_t(i,j,k-1,2)+tem_t(i,j,k-2,2))&!21-11-16
!             +x2*(tem_t(i,j,k+1,2)-2.*tem_t(i,j,k  ,2)+tem_t(i,j,k-1,2))&
!                                 )*(1+anyv3d(i,j,k,id_cn_))/6. ! QUICKEST
!!                                                         )/6. ! UBS-UPW


        anyv1d(k,4)=  &
          -anyv3d(i,j,k,id_cn_) *(x1*sal_t(i,j,k-1,2)+x2*sal_t(i,j,k  ,2)) !&
!      +(1-anyv3d(i,j,k,id_cn_))*( &
!              x1*(sal_t(i,j,k  ,2)-2.*sal_t(i,j,k-1,2)+sal_t(i,j,k-2,2))&!21-11-16
!             +x2*(sal_t(i,j,k+1,2)-2.*sal_t(i,j,k  ,2)+sal_t(i,j,k-1,2))&
!                                )*(1+anyv3d(i,j,k,id_cn_))/6. ! QUICKEST
!!                                                         )/6. ! UBS-UPW

       enddo       !nexon>

       do k=1,kmax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_tdiv)= & !26-01-17
       anyv3d(i,j,k,id_tdiv)+dti_lpsub*(anyv3d(i,j,k+1,id_we_)    &
                                        -anyv3d(i,j,k  ,id_we_))*tem_t(i,j,k,2) 

         tem_t(i,j,k,2)=tem_t(i,j,k,2)+dti_lpsub*wetmask_t(i,j)*( &

                              anyv1d(k+1,3)+anyv1d(k+1,1) &
                             -anyv1d(k  ,3)-anyv1d(k  ,1) &
         
                 +tem_t(i,j,k,2)*( anyv3d(i,j,k+1,id_we_)    & !26-01-17
                                  -anyv3d(i,j,k  ,id_we_))   &

                                                 )/dz_t(i,j,k,before)

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_sdiv)= & !25-01-17
       anyv3d(i,j,k,id_sdiv)+dti_lpsub*(anyv3d(i,j,k+1,id_we_)    &
                                        -anyv3d(i,j,k  ,id_we_))*sal_t(i,j,k,2) 

         sal_t(i,j,k,2)=sal_t(i,j,k,2)+dti_lpsub*wetmask_t(i,j)*( & 

                              anyv1d(k+1,4)+anyv1d(k+1,2) &
                             -anyv1d(k  ,4)-anyv1d(k  ,2) &
         
                 +sal_t(i,j,k,2)*( anyv3d(i,j,k+1,id_we_)    & !25-01-17
                                  -anyv3d(i,j,k  ,id_we_))   &

                                                 )/dz_t(i,j,k,before)

       enddo


       enddo               ! boucle iterative >>>


      enddo ; enddo ! Boucle i,j No1

!#ifdef bidon
#ifdef bidon
      do k=1,kmax ; do j=1,jmax  ; do i=1,imax !22-11-16  


           tem_t(i,j,k,2)                                &
          =tem_t(i,j,k,2)                                &

!     -wetmask_t(i,j)*tem_t(i,j,k,now)*dti_lp*(          &
!     -wetmask_t(i,j)*tem_t(i,j,k,before)*dti_lp*(       & !16-01-17

!              ( veldydz_u(i+1,j  ,k  ,1)                &
!               -veldydz_u(i  ,j  ,k  ,1)                &
!               +veldxdz_v(i  ,j+1,k  ,1)                &
!               -veldxdz_v(i  ,j  ,k  ,1) )/dxdy_t(i,j)  &
!                  +anyv3d(i  ,j  ,k+1,id_we_)           &
!                  -anyv3d(i  ,j  ,k  ,id_we_)           &

!                                             )/dz_t(i,j,k,before)
      -wetmask_t(i,j)*anyv3d(i,j,k,id_tdiv)/dz_t(i,j,k,before) !26-01-17

           sal_t(i,j,k,2)                                &
          =sal_t(i,j,k,2)                                &

!     -wetmask_t(i,j)*sal_t(i,j,k,now)*dti_lp*(          &
!     -wetmask_t(i,j)*sal_t(i,j,k,before)*dti_lp*(       & !16-01-17

!              ( veldydz_u(i+1,j  ,k  ,1)                &
!               -veldydz_u(i  ,j  ,k  ,1)                &
!               +veldxdz_v(i  ,j+1,k  ,1)                &
!               -veldxdz_v(i  ,j  ,k  ,1) )/dxdy_t(i,j)  &
!                  +anyv3d(i  ,j  ,k+1,id_we_)           &
!                  -anyv3d(i  ,j  ,k  ,id_we_)           &

!                                             )/dz_t(i,j,k,before)
      -wetmask_t(i,j)*anyv3d(i,j,k,id_sdiv)/dz_t(i,j,k,before) !25-01-17

      enddo ; enddo ; enddo
#endif
      check5=anyv3d(imax/2,jmax/2,kmax-1,id_tdiv)
      check6=anyv3d(imax/2,jmax/2,kmax-1,id_sdiv)
!#endif

!...............
! La vitesse residuelle implicite est omega_w moins omega_explicite 


! Implicit "up2" (upwind 2 points) part of the vertical advection: !20-08-14
      do k=2,kmax-1 !kmax !09-04-15 !10-04-19
      do j=1,jmax
      do i=1,imax

!     tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
!      anyv3d(i,j,k+1,id_cn_)*( anyv3d(i,j,k+1,id_wi_)+abs(anyv3d(i,j,k+1,id_wi_))) &
!     +anyv3d(i,j,k  ,id_cn_)*(-anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

!     tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
!      anyv3d(i,j,k  ,id_cn_)*(-anyv3d(i,j,k,id_wi_)-abs(anyv3d(i,j,k,id_wi_)))

!     tridia_in(i,j,k,3)=tridia_in(i,j,k,3)+0.5*dti_lp*                &
!      anyv3d(i,j,k+1,id_cn_)*( anyv3d(i,j,k+1,id_wi_)-abs(anyv3d(i,j,k+1,id_wi_)))

      tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
       ( anyv3d(i,j,k+1,id_wi_)+abs(anyv3d(i,j,k+1,id_wi_))) &
      +(-anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

      tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
       (-anyv3d(i,j,k,id_wi_)-abs(anyv3d(i,j,k,id_wi_)))

      tridia_in(i,j,k,3)=tridia_in(i,j,k,3)+0.5*dti_lp*                &
       ( anyv3d(i,j,k+1,id_wi_)-abs(anyv3d(i,j,k+1,id_wi_)))

      enddo
      enddo
      enddo

! A la surface on suppose une condition de gradient nul (ce qui signifie que !09-04-15
! la variation de T,S due aux flux de surface passera par le membre de droite
! dans vertmix.F90) sans modifier le protocole habituel.
! Pour obtenir les lignes suivantes, on part des lignes precedentes (cas "dans la couche")
! A1*T(k-1)+A2*T(k)+A3*T(k+1) et on remplace T(k+1) par T(k) pour avoir la condition de
! gradient nul et deduire les nouveaux coefficients: A1=A1 ; A2=A2+A3 ; A3=0. 
! Notons qu'on a force le schema a etre implicite, autrement dit anyv3d(i,j,k+1,0)=0
! et disparait du calcul. On obtient:
      k=kmax
      do j=1,jmax ; do i=1,imax

!     tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
!      anyv3d(i,j,k  ,id_cn_)*(-anyv3d(i,j,k,id_wi_)-abs(anyv3d(i,j,k,id_wi_)))

!     tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
!     +anyv3d(i,j,k+1,id_cn_)*2.*anyv3d(i,j,k+1,id_wi_)                &
!     +anyv3d(i,j,k  ,id_cn_)*( -anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

      tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
       (-anyv3d(i,j,k,id_wi_)-abs(anyv3d(i,j,k,id_wi_)))

      tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
      +2.*anyv3d(i,j,k+1,id_wi_)                &
      +( -anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

      enddo ; enddo

! Meme raisonnement au fond pour ce qui concerne les sources sous marine
! conduisant a A2=A2+A1 !10-04-19
      k=1
      do j=1,jmax ; do i=1,imax

      tridia_in(i,j,k,3)=tridia_in(i,j,k,3)+0.5*dti_lp*             &
       ( anyv3d(i,j,k+1,id_wi_)-abs(anyv3d(i,j,k+1,id_wi_)))

      tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(            &
       (   anyv3d(i,j,k+1,id_wi_)+abs(anyv3d(i,j,k+1,id_wi_)))      &
      +(-2*anyv3d(i,j,k  ,id_wi_)                            ) )

      enddo ; enddo

      end subroutine advection_scal_z_c2limited

!................................................................
#ifdef bidon
      subroutine advection_scal_surfavr(id_)
      use module_principal
      use module_parallele
      implicit none
      integer id_


! No merged levels if flag_merged_levels=0
      if(flag_merged_levels==0)return

! Moyenner anyv3d sur les couches de surface "collees" de int(truekmax_t))+1 A kmax

      do j=1,jmax
      do i=1,imax
       k1=int(truekmax_t(i,j))+1
       sum1=small2
       sum2=0.
       do k=k1,kmax
! On choisit 4eme arg=2 car a priori coherent avec critere de conservation des traceurs
        sum1=sum1+dz_t(i,j,k,2) 
        sum2=sum2+dz_t(i,j,k,2)*anyv3d(i,j,k,id_)
       enddo
       x1=sum2/sum1
       do k=k1,kmax
        anyv3d(i,j,k,id_)=x1
!       if(j==jmax/2)write(64,*)i,depth_t(i,j,k) ! A superposer au fichier tmp/z2dv_OE_
       enddo
      enddo
      enddo

!     stop 'test fichier'

      end subroutine advection_scal_surfavr
#endif
!................................................................

      subroutine advection_scal_z_quickest !12-11-16
      use module_principal
      use module_parallele
      implicit none
!     double precision :: grdamp_=2. ! lower gradient amplification
      integer :: id_cn_=0           &! current number "limited"
                ,id_we_=1           &! Explicit limited vertical velocity
                ,id_wi_=2           &! residual implicit vertical velocity
!               ,id_tdiv=4         &! cumulated T*div(current) !26-01-17
!               ,id_sdiv=5         &! cumulated S*div(current)
                ,loop_,loopmax_     &
                ,looplimit_=1        ! nombre d'iteration maximum la procedure iterative
#ifdef synopsis
       subroutinetitle='advection_scal_z_quickest'
       subroutinedescription='advection_scal_z_quickest'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


!     x1=0. ! Pas obligatoire seulement pour diag du maximum
      xy_t(:,:,1)=small2 ! reset small2 evite que le nombre d'iteration soit nulle
      do k=2,kmax ; do j=1,jmax ; do i=1,imax

!   x0=Dz/Dt
        x0=min(dz_t(i,j,k,0),dz_t(i,j,k-1,0))/dti_lp
!   Omega bornee par + ou - N fois dz/dt ici 2 fois
        anyv3d(i,j,k,id_we_)=max(min(omega_w(i,j,k,1),looplimit_*x0),-looplimit_*x0)*wetmask_t(i,j)  &
                            *wetmask_wi_t(i,j) !19-06-19

!   max over z of the decimal current number:
        xy_t(i,j,1)=max(xy_t(i,j,1),abs(anyv3d(i,j,k,id_we_))/x0)
!   x1=max(x1,xy_t(i,j,1)) ! Pas obligatoire seulement pour diag du maximum

      enddo  ; enddo ; enddo

! Fond w explicite = 0
      do j=1,jmax ; do i=1,imax
        anyv3d(i,j,1,id_we_)=0.
      enddo       ; enddo

! Couche merged (dans laquelle w est 100% implicite, donc w explicite = 0)
      if(flag_merged_levels==1) then !pmxpmx> !18-02-17
       do j=1,jmax ; do i=1,imax
        do k=2,kmerged_t(i,j)
         anyv3d(i,j,k,id_we_)=0.
        enddo
       enddo       ; enddo
      endif                          !pmxpmx> !18-02-17

! Surface:
      do j=1,jmax ; do i=1,imax
!      anyv3d(i,j,kmax+1,id_we_)=omega_w(i,j,kmax+1,1) ! A la surface w explicite pour flux de surface
       anyv3d(i,j,kmax+1,id_we_)=0. !25-02-19
      enddo       ; enddo

!....................................................................
! Decommenter ces lignes pour avoir la valeur max sur tout le domaine !30-12-16
!     call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_max,par%comm2d,ierr)
!     if(par%rank==0) then !>>>>>
!       open(unit=3,file='tmp/dti_w',position='append')
!        write(3,*)real(elapsedtime_now/86400.),real(x2)
!       close(3)
!     endif                !>>>>> 
!....................................................................

! Vitesse verticale residuelle implicite:
      do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,id_wi_)=omega_w(i,j,k,1)-anyv3d(i,j,k,id_we_)
      enddo       ; enddo       ; enddo

! Dans la couche merged, profil de vitesse verticale implicite de type "gradient de flux constant"
!     if(flag_merged_levels==1) then !pmxpmx> !18-02-17
!      do j=1,jmax ; do i=1,imax
!       do k=2,kmerged_t(i,j)
!        anyv3d(i,j,k,id_wi_)=anyv3d(i,j,kmerged_t(i,j)+1,id_wi_)       &
!                          *(depth_w(i,j,k)            -depth_w(i,j,1)) &
!                          /(depth_w(i,j,kmerged_t(i,j)+1)-depth_w(i,j,1))  
!       enddo
!      enddo       ; enddo
!     endif                          !pmxpmx> !18-02-17

      do j=1,jmax ; do i=1,imax ! Boucle i,j No1

! ceiling de xy_t(i,j,1) nombre de sous-iterations advectives
       loopmax_=min( ceiling(xy_t(i,j,1)) , looplimit_ )
       dti_lpsub=dti_lp/loopmax_

! Limiteurs C1,C2 ou C3
! https://docs.google.com/document/d/1WsIhfKI1GP269wglgIM6YUKXdVrONILKupKP9Gui9hI/edit
      do k=2,kmax-1
!      if(rhp_t(i,j,k  )-rhp_t(i,j,k-1)>=0) then
       if(rhp_t(i,j,k  )-rhp_t(i,j,k-1)>0.01*(rhp_t(i,j,k+1)-rhp_t(i,j,k))) then !m°v°m>
! if(rhp_t(i,j,k)*(depth_t(i,j,k+1)-depth_t(i,j,k-1)-dti_fw*0.5*(omega_w(i,j,k 
!!! suggestion de test: if((rhp(k)-rhp(k-1))/(small1+abs(rhp(kmax)-rhp(kmin)))>-1e-n afin de remplacer le test
!!! >=0 par une valeur legerement negative anticipative confrontee au gradient normalisE par rhp(surf)-rhp(bot)
            anyv1d(k,1)=0.
       else                                                                      !m°v°m>
            anyv1d(k,1)=1.
       endif

! Proposition de fonction continue: anyv1d(k,1)=0. (100 %up) si (rhp(k)-rhp(k-1))/(rhp(k+1)-rhp(k))=0.
!                                   anyv1d(k,1)=1. (100 %c2) si (rhp(k)-rhp(k-1))/(rhp(k+1)-rhp(k))=0.1
! anyv1d(k,1)=min(max( (rhp_t(i,j,k)-rhp_t(i,j,k-1))/min(0.1*(rhp_t(i,j,k+1)-rhp_t(i,j,k)),small1),0.),1.)
! Autrement dit si rhp(k)-rhp(k-1) est nul alors 100% upwind
!                  rhp(k)-rhp(k-1) est 10 fois plus petite que 
!                  la pente du point au dessus, alors 100% c2 et
!                  commence la transition continue vers upwind

      enddo                                                                      !m°v°m>
!     anyv1d(kmax+1,1)=1. ! inutile si C1
!     anyv1d(kmax  ,1)=1. ! inutile si C1
      anyv1d(     1,1)=1.

      do k=2,kmax  ! do k=kmin_w(i,j)+1,kmax                         !22-12-14

! Current number with additional limititations
       anyv3d(i,j,k,id_cn_)=1.-( &!ooooo>
! Robustness regarding convection: 
! https://docs.google.com/document/d/1WsIhfKI1GP269wglgIM6YUKXdVrONILKupKP9Gui9hI/edit
!                      anyv1d(k+1,1)*         & ! limiteur C3 cette ligne et les 2 lignes suivantes
!                      anyv1d(k  ,1)*         & ! limiteur C2 cette ligne et la ligne suivante
                       anyv1d(k-1,1)*         & ! limiteur C1 cette ligne seulement !23-05-18
! Robustness regarding rivers and dry zones:
                       upwindriver_t(i,j)   & ! river zone influence
                         *upwzone0_t(i,j,k) & ! additional requests for upwind. Note upwzone0_t=1 in standard conditions
                      *wetmask_t(i,j)       & ! wet/dry zone
! Robustness regarding vertical current number:
        *max(1.-abs(anyv3d(i,j,k,id_we_)*dti_lpsub               &
                 /min(dz_t(i,j,k,0),dz_t(i,j,k-1,0))),0.)        &  
                                ) !ooooo>

!        anyv3d(i,j,k,id_cn_)=0. ! 100%C2
!        anyv3d(i,j,k,id_cn_)=1. ! 100%UP2

       enddo


! Cas particuliers k=1:kmin_w(i,j) et k=kmax:
!      anyv3d(i,j,1:kmin_w(i,j),id_cn_)=0.
       anyv3d(i,j,1            ,id_cn_)=0.
! https://docs.google.com/document/d/1fIzAG9mo_zvtTdVauu7yw-VkIqFRldqq4tF2b6ZOnzI/edit
!      anyv3d(i,j,kmax+1       ,id_cn_)=0. ! A la surface schema 100% explicite pour
!                                          ! garantir la coherence du flux d'eau douce
!                                          ! qui, dans vertmix_sal, est calcule avec 
!                                          ! avec sal_t(:,:,kmax,1). (Conservation bilan 
!                                          ! de sel) !12-04-15
       anyv3d(i,j,kmax+1       ,id_cn_)=1. ! A la surface schema upwind mais noter
! surtout que le choix nc=0 ou 1 n'a pas d'importance si la vitesse verticale est
! 100% implicite en kmax+1
! Voir mes commentaires à propos du passage en implicite en kmax+1:
! https://docs.google.com/document/d/1fIzAG9mo_zvtTdVauu7yw-VkIqFRldqq4tF2b6ZOnzI/edit


! Flux centres frozen:
       do k=2,kmax

        anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)  &
                  *( ( dz_t(i,j,k  ,1)*tem_t(i,j,k  ,1)             &
                      +dz_t(i,j,k-1,1)*tem_t(i,j,k-1,1))            &
                    /( dz_t(i,j,k  ,1) +dz_t(i,j,k-1,1)) )      

        anyv1d(k,2)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)  &
                  *( ( dz_t(i,j,k  ,1)*sal_t(i,j,k  ,1)             &
                      +dz_t(i,j,k-1,1)*sal_t(i,j,k-1,1))            &
                    /( dz_t(i,j,k  ,1) +dz_t(i,j,k-1,1)) )      

       enddo
       anyv1d(1     ,1:4)=0.
       anyv1d(kmax+1,3:4)=0.
       k=kmax+1
       anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)*tem_t(i,j,k-1,1)
       anyv1d(k,2)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)*sal_t(i,j,k-1,1)

! Partie Iterative:
       do loop_=1,loopmax_ ! boucle iterative >>>

! A priori lignes suivantes inutiles si schema implicite en kmax+1
       tem_t(i,j,kmax+1,2)=tem_t(i,j,kmax,2)
       tem_t(i,j,0     ,2)=tem_t(i,j,1   ,2)
       sal_t(i,j,kmax+1,2)=sal_t(i,j,kmax,2)
       sal_t(i,j,0     ,2)=sal_t(i,j,1   ,2)
       do k=2,kmax !nexon>

        x1=0.5*(anyv3d(i,j,k,id_we_)+abs(anyv3d(i,j,k,id_we_)))
        x2=0.5*(anyv3d(i,j,k,id_we_)-abs(anyv3d(i,j,k,id_we_)))

! Flux upwind:
        anyv1d(k,3)=  &
          -anyv3d(i,j,k,id_cn_) *(x1*tem_t(i,j,k-1,2)+x2*tem_t(i,j,k  ,2))&
       +(1-anyv3d(i,j,k,id_cn_))*( &
               x1*(tem_t(i,j,k  ,2)-2.*tem_t(i,j,k-1,2)+tem_t(i,j,k-2,2))&!21-11-16
              +x2*(tem_t(i,j,k+1,2)-2.*tem_t(i,j,k  ,2)+tem_t(i,j,k-1,2))&
                                  )*(1+anyv3d(i,j,k,id_cn_))/6. ! QUICKEST
!                                                          )/6. ! UBS-UPW


!       anyv1d(k,4)=-anyv3d(i,j,k,id_cn_)*(x1*sal_t(i,j,k-1,2)+x2*sal_t(i,j,k  ,2))&
        anyv1d(k,4)=  &
          -anyv3d(i,j,k,id_cn_) *(x1*sal_t(i,j,k-1,2)+x2*sal_t(i,j,k  ,2))&
       +(1-anyv3d(i,j,k,id_cn_))*( &
               x1*(sal_t(i,j,k  ,2)-2.*sal_t(i,j,k-1,2)+sal_t(i,j,k-2,2))&!21-11-16
              +x2*(sal_t(i,j,k+1,2)-2.*sal_t(i,j,k  ,2)+sal_t(i,j,k-1,2))&
                                 )*(1+anyv3d(i,j,k,id_cn_))/6. ! QUICKEST
!                                                          )/6. ! UBS-UPW

       enddo       !nexon>

       do k=1,kmax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_tdiv)= & !26-01-17
       anyv3d(i,j,k,id_tdiv)+dti_lpsub*(anyv3d(i,j,k+1,id_we_)    &
                                        -anyv3d(i,j,k  ,id_we_))*tem_t(i,j,k,2) 

         tem_t(i,j,k,2)=tem_t(i,j,k,2)+dti_lpsub*wetmask_t(i,j)*( &

                              anyv1d(k+1,3)+anyv1d(k+1,1) &
                             -anyv1d(k  ,3)-anyv1d(k  ,1) &
         
                 +tem_t(i,j,k,2)*( anyv3d(i,j,k+1,id_we_)    & !26-01-17
                                  -anyv3d(i,j,k  ,id_we_))   &

                                                 )/dz_t(i,j,k,before)

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_sdiv)= & !25-01-17
       anyv3d(i,j,k,id_sdiv)+dti_lpsub*(anyv3d(i,j,k+1,id_we_)    &
                                        -anyv3d(i,j,k  ,id_we_))*sal_t(i,j,k,2) 

         sal_t(i,j,k,2)=sal_t(i,j,k,2)+dti_lpsub*wetmask_t(i,j)*( & 

                              anyv1d(k+1,4)+anyv1d(k+1,2) &
                             -anyv1d(k  ,4)-anyv1d(k  ,2) &
         
                 +sal_t(i,j,k,2)*( anyv3d(i,j,k+1,id_we_)    & !25-01-17
                                  -anyv3d(i,j,k  ,id_we_))   &

                                                 )/dz_t(i,j,k,before)

       enddo


       enddo               ! boucle iterative >>>


      enddo ; enddo ! Boucle i,j No1

#ifdef bidon
      do k=1,kmax ; do j=1,jmax  ; do i=1,imax !22-11-16  

           tem_t(i,j,k,2)                                &
          =tem_t(i,j,k,2)                                &

!     -wetmask_t(i,j)*tem_t(i,j,k,now)*dti_lp*(          &
!     -wetmask_t(i,j)*tem_t(i,j,k,before)*dti_lp*(       & !16-01-17

!              ( veldydz_u(i+1,j  ,k  ,1)                &
!               -veldydz_u(i  ,j  ,k  ,1)                &
!               +veldxdz_v(i  ,j+1,k  ,1)                &
!               -veldxdz_v(i  ,j  ,k  ,1) )/dxdy_t(i,j)  &
!                  +anyv3d(i  ,j  ,k+1,id_we_)           &
!                  -anyv3d(i  ,j  ,k  ,id_we_)           &

!                                             )/dz_t(i,j,k,before)
      -wetmask_t(i,j)*anyv3d(i,j,k,id_tdiv)/dz_t(i,j,k,before) !26-01-17

           sal_t(i,j,k,2)                                &
          =sal_t(i,j,k,2)                                &

!     -wetmask_t(i,j)*sal_t(i,j,k,now)*dti_lp*(          &
!     -wetmask_t(i,j)*sal_t(i,j,k,before)*dti_lp*(       & !16-01-17

!              ( veldydz_u(i+1,j  ,k  ,1)                &
!               -veldydz_u(i  ,j  ,k  ,1)                &
!               +veldxdz_v(i  ,j+1,k  ,1)                &
!               -veldxdz_v(i  ,j  ,k  ,1) )/dxdy_t(i,j)  &
!                  +anyv3d(i  ,j  ,k+1,id_we_)           &
!                  -anyv3d(i  ,j  ,k  ,id_we_)           &

!                                             )/dz_t(i,j,k,before)
      -wetmask_t(i,j)*anyv3d(i,j,k,id_sdiv)/dz_t(i,j,k,before) !25-01-17

      enddo ; enddo ; enddo
#endif
      check5=anyv3d(imax/2,jmax/2,kmax-1,id_tdiv)
      check6=anyv3d(imax/2,jmax/2,kmax-1,id_sdiv)

!     if(flag_merged_levels==1)call vertmix_merged_levels_tem(0)
!     if(flag_merged_levels==1)call vertmix_merged_levels_sal(0)

!...............
! La vitesse residuelle implicite est omega_w moins omega_explicite 

! Implicit "up2" (upwind 2 points) part of the vertical advection: !20-08-14
      do k=2,kmax-1 !kmax !09-04-15 !10-04-19
      do j=1,jmax
      do i=1,imax

!     tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
!      anyv3d(i,j,k+1,id_cn_)*( anyv3d(i,j,k+1,id_wi_)+abs(anyv3d(i,j,k+1,id_wi_))) &
!     +anyv3d(i,j,k  ,id_cn_)*(-anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

!     tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
!      anyv3d(i,j,k  ,id_cn_)*(-anyv3d(i,j,k,id_wi_)-abs(anyv3d(i,j,k,id_wi_)))

!     tridia_in(i,j,k,3)=tridia_in(i,j,k,3)+0.5*dti_lp*                &
!      anyv3d(i,j,k+1,id_cn_)*( anyv3d(i,j,k+1,id_wi_)-abs(anyv3d(i,j,k+1,id_wi_)))

      tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
       ( anyv3d(i,j,k+1,id_wi_)+abs(anyv3d(i,j,k+1,id_wi_))) &
      +(-anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

      tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
       (-anyv3d(i,j,k,id_wi_)-abs(anyv3d(i,j,k,id_wi_)))

      tridia_in(i,j,k,3)=tridia_in(i,j,k,3)+0.5*dti_lp*                &
       ( anyv3d(i,j,k+1,id_wi_)-abs(anyv3d(i,j,k+1,id_wi_)))

      enddo
      enddo
      enddo

! A la surface on suppose une condition de gradient nul (ce qui signifie que !09-04-15
! la variation de T,S due aux flux de surface passera par le membre de droite
! dans vertmix.F90) sans modifier le protocole habituel.
! Pour obtenir les lignes suivantes, on part des lignes precedentes (cas "dans la couche")
! A1*T(k-1)+A2*T(k)+A3*T(k+1) et on remplace T(k+1) par T(k) pour avoir la condition de
! gradient nul et deduire les nouveaux coefficients: A1=A1 ; A2=A2+A3 ; A3=0. 
! Notons qu'on a force le schema a etre implicite, autrement dit anyv3d(i,j,k+1,0)=0
! et disparait du calcul. On obtient:
      k=kmax
      do j=1,jmax ; do i=1,imax

!     tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
!      anyv3d(i,j,k  ,id_cn_)*(-anyv3d(i,j,k,id_wi_)-abs(anyv3d(i,j,k,id_wi_)))

!     tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
!     +anyv3d(i,j,k+1,id_cn_)*2.*anyv3d(i,j,k+1,id_wi_)                &
!     +anyv3d(i,j,k  ,id_cn_)*( -anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

      tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
       (-anyv3d(i,j,k,id_wi_)-abs(anyv3d(i,j,k,id_wi_)))

      tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
      +2.*anyv3d(i,j,k+1,id_wi_)                &
      +( -anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

      enddo ; enddo

! Meme raisonnement au fond pour ce qui concerne les sources sous marine
! conduisant a A2=A2+A1 !10-04-19
      k=1
      do j=1,jmax ; do i=1,imax

      tridia_in(i,j,k,3)=tridia_in(i,j,k,3)+0.5*dti_lp*             &
       ( anyv3d(i,j,k+1,id_wi_)-abs(anyv3d(i,j,k+1,id_wi_)))

      tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(            &
       (   anyv3d(i,j,k+1,id_wi_)+abs(anyv3d(i,j,k+1,id_wi_)))      &
      +(-2*anyv3d(i,j,k  ,id_wi_)                            ) )

      enddo ; enddo


      end subroutine advection_scal_z_quickest

!................................................................

      subroutine advection_scal_uplim_subt  ! 05-01-16
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      double precision :: velef_,flagmsk_=1.,scalefac_,cst_
      integer :: t_=1        &
                ,id_flx_=1   &  ! x flux identifier
                ,id_fly_=1   &  ! y flux identifier
                ,id_tbef_=3  &  ! temperature    "before" identifier
!               ,id_tdiv=4  &  ! T*divergence(velocity) !26-01-17
                ,loop_
#ifdef synopsis
       subroutinetitle='advection_scal_uplim_subs'
       subroutinedescription=                                     &
          'Computes temperature and teminity advection diffusion' &
       //' horizontal fluxes using a "3 points" advection scheme' &
       //' hybrid QUICK-UPWIND'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Cas particulier du modele 1DV
!     if(flag_1dv==1) then !>>>>> !14-02-16
!      tem_t(:,:,:,2)=tem_t(:,:,:,0)
!      return
!     endif

#ifdef bidon
      sum1=0.
      sum2=0.
      do k=1,kmax
      do j=1,jmax
      do i=1,imax
       sum1=sum1+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)
       sum2=sum2+dz_t(i,j,k,1)*dxdy_t(i,j)*mask_t(i,j,k)*mask_i_w(i)*mask_j_w(j)*tem_t(i,j,k,1)
      enddo
      enddo
      enddo
      call mpi_allreduce(sum1,sum1glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      call mpi_allreduce(sum2,sum2glb,1,mpi_double_precision,mpi_sum,par%comm2d,ierr)
      if(par%rank==0)write(6,*)sum2glb,sum2glb/sum1glb
#endif



      cst_=0.5*cst_adv_hor/3. !11-09-17


!#ifdef bidon
      do k=1,kmax ; do j=-1,jmax+2 ; do i=-1,imax+2
       anyv3d(i,j,k,id_tbef_)=tem_t(i,j,k,before)
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,id_tdiv)=0. ! reset cumul temperature*divergence(courant) !26-01-17
      enddo ; enddo ; enddo

      do loop_=1,loopmaxts !*********>

! Flux Oi
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax
      do i=1,imax+1

      anyv3d(i,j,k,id_flx_)=0.5*( & !prime>


! Schema QUICK en dehors influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      anyv3d(i,j,k,id_hybcoefu)*( &
! Partie advective QUICK:
      -veldydz_u(i,j,k,1) &
!       *(anyv3d(i,j,k,id_tnow_)+anyv3d(i-1,j,k,id_tnow_)) &
! Si le nombre de courant est faible il n'y a pas de risque A maintenir la
! partie centree sur tem_t(now) ce qui permet de maintenir un leap-frog centrE
! classique meme en subdivisant le pas de temps, ce qui permet A priori de 
! maintenir les proprietes de conservation associeEs
        *(    tem_t(i,j,k,now)+tem_t(i-1,j,k,now)      &
#ifdef bidonref
       +2.*temref_u(i,j,k)                             &
          -temref_t(i,j,k) -temref_t(i-1,j,k)          &
#endif
                                                 )     &

! Partie diffusive QUICK remise A jour aux sous-pas de temps advectif:
       +cst_*0.5*(anyv3d(i,j,k,id_ofactort)+anyv3d(i,j+1,k,id_ofactort)) &

! efficacite densitaire
!       *abs(-alp_t*(   tem_t(i,j,k,1)-tem_t(i-1,j,k,1))+alp_s*(   sal_t(i,j,k,1)-sal_t(i-1,j,k,1))) &
!       /   (+alp_t*abs(tem_t(i,j,k,1)-tem_t(i-1,j,k,1))+alp_s*abs(sal_t(i,j,k,1)-sal_t(i-1,j,k,1))+small1) &

                                                               *( & !m°v°m>
             (veldydz_u(i,j,k,1)+abs(veldydz_u(i,j,k,1)))*( &           !pmxpmx>

               (   anyv3d(i  ,j,k,id_tbef_)   -anyv3d(i-1,j,k,id_tbef_)  &

! correction anyv3d(i  ,j,k,id_tbef_) par rapport niveau z(i-1,j,k)
      +( anyv3d(i  ,j,kp1,id_tbef_)- anyv3d(i  ,j,km1,id_tbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i-1,j,k  )         -depth_t(i  ,j,k  )         )    &
      /(depth_t(i  ,j,kp1)         -depth_t(i  ,j,km1)         )    &
          )       )      & !minmax>

               ) &

              -(   anyv3d(i-1,j,k,id_tbef_)-anyv3d(i-2,j,k,id_tbef_)  &

! correction anyv3d(i-2,j,k,id_tbef_) par rapport niveau z(i-1,j,k)
      -( anyv3d(i-2,j,kp1,id_tbef_)- anyv3d(i-2,j,km1,id_tbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i-1,j,k  )         -depth_t(i-2,j,k  )         )    &
      /(depth_t(i-2,j,kp1)         -depth_t(i-2,j,km1)         )    &
          )       )      & !minmax>

               ) &
                                           *mask_t(i-2,j,kmax)      ) & !pmxpmx>

            +(veldydz_u(i,j,k,1)-abs(veldydz_u(i,j,k,1)))*( &           !pmxpmx>

               (   anyv3d(i+1,j,k,id_tbef_)-anyv3d(i,j,k,id_tbef_)    &

! correction anyv3d(i+1,j,k,id_tbef_) par rapport niveau z(i  ,j,k)
      +( anyv3d(i+1,j,kp1,id_tbef_)- anyv3d(i+1,j,km1,id_tbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i  ,j,k  )         -depth_t(i+1,j,k  )         )    &
      /(depth_t(i+1,j,kp1)         -depth_t(i+1,j,km1)         )    &
          )       )      & !minmax>

               )   &
                  *mask_t(i+1,j,kmax)                                 &

              -(   anyv3d(i,j,k,id_tbef_)-anyv3d(i-1,j,k,id_tbef_)    &

! correction anyv3d(i-1,j,k,id_tbef_) par rapport niveau z(i  ,j,k)
      -( anyv3d(i-1,j,kp1,id_tbef_)- anyv3d(i-1,j,km1,id_tbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i  ,j,k  )         -depth_t(i-1,j,k  )         )    &
      /(depth_t(i-1,j,kp1)         -depth_t(i-1,j,km1)         )    &
          )       )      & !minmax>

                                                                  ))  & !pmxpmx>

                                                                )  & !m°v°m>()cst_
                                                                )  & !()anyv3d(i,j,k,id_hybcoefu)
! Schema UPWIND dans les zones sous influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      +(1.-anyv3d(i,j,k,id_hybcoefu))*( & !pmx>
           -veldydz_u(i,j,k,1)* ( & !ppp>
               anyv3d(i,j,k,id_tbef_)+anyv3d(i-1,j,k,id_tbef_)   &
                                ) & !ppp>
       +abs(veldydz_u(i,j,k,1))*( & !mmm>
               anyv3d(i,j,k,id_tbef_)-anyv3d(i-1,j,k,id_tbef_)   &
                                ) & !mmm>
                                      ) & !pmx>

       ) !prime>

!     -veldydz_u(i,j,k,1)*dti_lpsub*temref_u(i,j,k)


      enddo !k
      enddo !i
      enddo !j

! Partial Oi Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_tdiv)= & !26-01-17
       anyv3d(i,j,k,id_tdiv)  &
      +anyv3d(i,j,k,id_tbef_)*dti_lpsub*(veldydz_u(i+1,j,k,1)           &
                                        -veldydz_u(i  ,j,k,1))/dxdy_t(i,j)

!       anyv3d(i,j,k,id_tbef_)                       &
        anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !ooo> !25-04-17

        anyv3d(i,j,k,id_tbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i+1,j,k,id_flx_)                &
             -anyv3d(i  ,j,k,id_flx_)                &

         +(veldydz_u(i+1,j,k,1)                      &
          -veldydz_u(i  ,j,k,1))*anyv3d(i,j,k,id_tbef_) & !26-01-17


         )/dxdy_t(i,j) ) & !ooo>

         +(1-mask_t(i,j,k))*tem_t(i,j,k,before)*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

! Note: le terme anyv3d(i,j,k,id_tdiv) ne fait pas l'objet d'un traitement similaire pour le masque car
! si mask_t=0 alors veldydz_u(i:i+1,...)=0 (question: et les fleuves ?)
      enddo
      enddo
      enddo

!     call my_outputs_zonesalttempflux('it',id_flx_)

      call vertmix_merged_levels_any(0,id_tbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_tbef_)
      call obc_int_anyv3d(id_tbef_,'zb')

! Flux Oj
      cst_=0.5*cst_adv_hor/3.
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax+1
      do i=1,imax

      anyv3d(i,j,k,id_fly_)=0.5*( & !prime>

! Schema QUICK en dehors influence des fleuves
      anyv3d(i,j,k,id_hybcoefv)*( &
! Partie advective QUICK:
      -veldxdz_v(i,j,k,1)  &
!       *(anyv3d(i,j,k,id_tnow_)+anyv3d(i,j-1,k,id_tnow_)) &
! Si le nombre de courant est faible il n'y a pas de risque A maintenir la
! partie centree sur tem_t(now) ce qui permet de maintenir un leap-frog centrE
! classique meme en subdivisant le pas de temps ce qui permet A priori de 
! maintenir les proprietes de conservation associeE
        *(    tem_t(i,j,k,now)+tem_t(i,j-1,k,now)  &
#ifdef bidonref
       +2.*temref_v(i,j,k)                         &
          -temref_t(i,j,k) -temref_t(i,j-1,k)      &
#endif
                                                 ) &

! Partie diffusive QUICK:
       +cst_*0.5*(anyv3d(i,j,k,id_ofactort)+anyv3d(i+1,j,k,id_ofactort)) &

! efficacite densitaire
!       *abs(-alp_t*(   tem_t(i,j,k,1)-tem_t(i,j-1,k,1))+alp_s*(   sal_t(i,j,k,1)-sal_t(i,j-1,k,1))) &
!       /   (+alp_t*abs(tem_t(i,j,k,1)-tem_t(i,j-1,k,1))+alp_s*abs(sal_t(i,j,k,1)-sal_t(i,j-1,k,1))+small1) &

                                                                  *(  & !m°v°m>

             (veldxdz_v(i,j,k,1)+abs(veldxdz_v(i,j,k,1)))*(  &       !ooo>

            (   anyv3d(i,j  ,k,id_tbef_)-anyv3d(i,j-1,k,id_tbef_)  &

! correction anyv3d(i,j  ,k,id_tbef_) par rapport niveau z(i,j-1,k)
      +( anyv3d(i,j  ,kp1,id_tbef_)- anyv3d(i,j  ,km1,id_tbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j-1,k  )         -depth_t(i,j  ,k  )         )    &
      /(depth_t(i,j  ,kp1)         -depth_t(i,j  ,km1)         )    &
          )       )      & !minmax>

            ) &

           -(   anyv3d(i,j-1,k,id_tbef_)-anyv3d(i,j-2,k,id_tbef_)  &

! correction anyv3d(i,j-2,k,id_tbef_) par rapport niveau z(i,j-1,k)
      -( anyv3d(i,j-2,kp1,id_tbef_)- anyv3d(i,j-2,km1,id_tbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j-1,k  )         -depth_t(i,j-2,k  )         )    &
      /(depth_t(i,j-2,kp1)         -depth_t(i,j-2,km1)         )    &
          )       )      & !minmax>

            ) &
                                        *mask_t(i,j-2,kmax)       )& !ooo>

            +(veldxdz_v(i,j,k,1)-abs(veldxdz_v(i,j,k,1)))*(  &       !ooo>

            (   anyv3d(i,j+1,k,id_tbef_)-anyv3d(i,j,k,id_tbef_)    &

! correction anyv3d(i,j+1,k,id_tbef_) par rapport niveau z(i,j  ,k)
      +( anyv3d(i,j+1,kp1,id_tbef_)- anyv3d(i,j+1,km1,id_tbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j  ,k  )         -depth_t(i,j+1,k  )         )    &
      /(depth_t(i,j+1,kp1)         -depth_t(i,j+1,km1)         )    &
          )       )      & !minmax>

            )   &
               *mask_t(i,j+1,kmax)                                 &

           -(   anyv3d(i,j,k,id_tbef_)-anyv3d(i,j-1,k,id_tbef_)    &

! correction anyv3d(i,j-1,k,id_tbef_) par rapport niveau z(i,j  ,k)
      -( anyv3d(i,j-1,kp1,id_tbef_)- anyv3d(i,j-1,km1,id_tbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j  ,k  )         -depth_t(i,j-1,k  )         )    &
      /(depth_t(i,j-1,kp1)         -depth_t(i,j-1,km1)         )    &
          )       )      & !minmax>

                                                               ))  & !ooo>

                                                                   )  & !m°v°m>()cst_
                                                                   )  & !()anyv3d(i,j,k,id_hybcoefv)
! Schema UPWIND dans les zones sous influence des fleuves:
      +(1.-anyv3d(i,j,k,id_hybcoefv))*( & !pmx>
          -veldxdz_v(i,j,k,1)* ( & !ppp>
              anyv3d(i,j,k,id_tbef_)+anyv3d(i,j-1,k,id_tbef_)   &
                               ) & !ppp>
      +abs(veldxdz_v(i,j,k,1))*( & !mmm>
              anyv3d(i,j,k,id_tbef_)-anyv3d(i,j-1,k,id_tbef_)   &
                               ) & !mmm>
                                      ) & !pmx>

       ) !prime>

!      -veldxdz_v(i,j,k,1)*dti_lpsub*temref_v(i,j,k)

      enddo !k
      enddo !i
      enddo !j

! Partial Horizontal Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_tdiv)= & !26-01-17
       anyv3d(i,j,k,id_tdiv)  &
      +anyv3d(i,j,k,id_tbef_)*dti_lpsub*(veldxdz_v(i,j+1,k,1)                 &
                                        -veldxdz_v(i,j  ,k,1))/dxdy_t(i,j)

        anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !ooo>
        anyv3d(i,j,k,id_tbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i,j+1,k,id_fly_)                &
             -anyv3d(i,j  ,k,id_fly_)                &

         +(veldxdz_v(i,j+1,k,1)                      &
          -veldxdz_v(i,j  ,k,1))*anyv3d(i,j,k,id_tbef_) & !26-01-17

         )/dxdy_t(i,j)  ) & !ooo>
 
         +(1-mask_t(i,j,k))*tem_t(i,j,k,before)*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

      enddo
      enddo
      enddo

!     call my_outputs_zonesalttempflux('jt',id_fly_) 

      call vertmix_merged_levels_any(0,id_tbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_tbef_)
      call obc_int_anyv3d(id_tbef_,'zb')

!#endif
      enddo               !*********>

      do k=1,kmax ; do j=1,jmax  ; do i=1,imax  

           tem_t(i,j,k,2)                              &
       =  anyv3d(i,j,k,id_tbef_)                    !  &

!     -wetmask_t(i,j)*tem_t(i,j,k,now)*dti_lp*(        &

!                                              veldydz_u(i+1,j,k,1)   &
!                                             -veldydz_u(i  ,j,k,1)   &
!                                             +veldxdz_v(i,j+1,k,1)   &
!                                             -veldxdz_v(i,j  ,k,1)   &

!                                             )/(dz_t(i,j,k,before)*dxdy_t(i,j))

      enddo ; enddo ; enddo

      end subroutine advection_scal_uplim_subt


!..............................................................................

      subroutine advection_scal_uplim_subs  ! 05-01-16
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      double precision :: velef_,flagmsk_=1.,scalefac_,cst_
      integer :: t_=1        &
                ,id_flx_=1   &  ! x flux identifier
                ,id_fly_=1   &  ! y flux identifier
                ,id_sbef_=3  &  ! salinity    "before" identifier
!               ,id_sdiv=5  &  ! S*divergence(velocity)
                ,loop_
#ifdef synopsis
       subroutinetitle='advection_scal_uplim_subs'
       subroutinedescription=                                     &
          'Computes salperature and salinity advection diffusion' &
       //' horizontal fluxes using a "3 points" advection scheme' &
       //' hybrid QUICK-UPWIND'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Cas particulier du modele 1DV
!     if(flag_1dv==1) then !>>>>> !14-02-16
!      sal_t(:,:,:,2)=sal_t(:,:,:,0)
!      return
!     endif

      cst_=0.5*cst_adv_hor/3.


!#ifdef bidon
      do k=1,kmax ; do j=-1,jmax+2 ; do i=-1,imax+2
       anyv3d(i,j,k,id_sbef_)=sal_t(i,j,k,before) 
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,id_sdiv)=0. ! reset cumul salinity*divergence(courant)
      enddo ; enddo ; enddo


      do loop_=1,loopmaxts !*********>


! Flux Oi
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax
      do i=1,imax+1

      anyv3d(i,j,k,id_flx_)=0.5*( & !prime>


! Schema QUICK en dehors influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      anyv3d(i,j,k,id_hybcoefu)*( &
! Partie advective QUICK:
      -veldydz_u(i,j,k,1) &
!       *(anyv3d(i,j,k,id_tnow_)+anyv3d(i-1,j,k,id_tnow_)) &
! Si le nombre de courant est faible il n'y a pas de risque A maintenir la
! partie centree sur sal_t(now) ce qui permet de maintenir un leap-frog centrE
! classique meme en subdivisant le pas de salps, ce qui permet A priori de 
! maintenir les proprietes de conservation associeEs
        *(    sal_t(i,j,k,now)+sal_t(i-1,j,k,now)      &
#ifdef bidonref
       +2.*salref_u(i,j,k)                             &
          -salref_t(i,j,k) -salref_t(i-1,j,k)          &
#endif
                                                 )     &

! Partie diffusive QUICK remise A jour aux sous-pas de salps advectif:
       +cst_*0.5*(anyv3d(i,j,k,id_ofactors)+anyv3d(i,j+1,k,id_ofactors))  &

! Efficacite densitaire
!       *abs(-alp_t*(   tem_t(i,j,k,1)-tem_t(i-1,j,k,1))+alp_s*(   sal_t(i,j,k,1)-sal_t(i-1,j,k,1))) &
!       /   (+alp_t*abs(tem_t(i,j,k,1)-tem_t(i-1,j,k,1))+alp_s*abs(sal_t(i,j,k,1)-sal_t(i-1,j,k,1))+small1) &

                                                               *( & !m°v°m>
             (veldydz_u(i,j,k,1)+abs(veldydz_u(i,j,k,1)))*( &           !pmxpmx>

               (   anyv3d(i  ,j,k,id_sbef_)-anyv3d(i-1,j,k,id_sbef_)  &

! correction anyv3d(i  ,j,k,id_sbef_) par rapport niveau z(i-1,j,k)
      +( anyv3d(i  ,j,kp1,id_sbef_)- anyv3d(i  ,j,km1,id_sbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i-1,j,k  )         -depth_t(i  ,j,k  )         )    &
      /(depth_t(i  ,j,kp1)         -depth_t(i  ,j,km1)         )    &
          )       )      & !minmax>

               ) &

              -(   anyv3d(i-1,j,k,id_sbef_)-anyv3d(i-2,j,k,id_sbef_)  &

! correction anyv3d(i-2,j,k,id_sbef_) par rapport niveau z(i-1,j,k)
      -( anyv3d(i-2,j,kp1,id_sbef_)- anyv3d(i-2,j,km1,id_sbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i-1,j,k  )         -depth_t(i-2,j,k  )         )    &
      /(depth_t(i-2,j,kp1)         -depth_t(i-2,j,km1)         )    &
          )       )      & !minmax>

               ) &
                                           *mask_t(i-2,j,kmax)      ) & !pmxpmx>

            +(veldydz_u(i,j,k,1)-abs(veldydz_u(i,j,k,1)))*( &           !pmxpmx>

               (   anyv3d(i+1,j,k,id_sbef_)-anyv3d(i,j,k,id_sbef_)    &

! correction anyv3d(i+1,j,k,id_sbef_) par rapport niveau z(i  ,j,k)
      +( anyv3d(i+1,j,kp1,id_sbef_)- anyv3d(i+1,j,km1,id_sbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i  ,j,k  )         -depth_t(i+1,j,k  )         )    &
      /(depth_t(i+1,j,kp1)         -depth_t(i+1,j,km1)         )    &
          )       )      & !minmax>

               )   &
                  *mask_t(i+1,j,kmax)                                 &

              -(   anyv3d(i,j,k,id_sbef_)-anyv3d(i-1,j,k,id_sbef_)    &

! correction anyv3d(i-1,j,k,id_sbef_) par rapport niveau z(i  ,j,k)
      -( anyv3d(i-1,j,kp1,id_sbef_)- anyv3d(i-1,j,km1,id_sbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i  ,j,k  )         -depth_t(i-1,j,k  )         )    &
      /(depth_t(i-1,j,kp1)         -depth_t(i-1,j,km1)         )    &
          )       )      & !minmax>

                                                                  ))  & !pmxpmx>

                                                                )  & !m°v°m>()cst_
                                                                )  & !()anyv3d(i,j,k,id_hybcoefu)
! Schema UPWIND dans les zones sous influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      +(1.-anyv3d(i,j,k,id_hybcoefu))*( & !pmx>
           -veldydz_u(i,j,k,1)* ( & !ppp>
               anyv3d(i,j,k,id_sbef_)+anyv3d(i-1,j,k,id_sbef_)   &
                                ) & !ppp>
       +abs(veldydz_u(i,j,k,1))*( & !mmm>
               anyv3d(i,j,k,id_sbef_)-anyv3d(i-1,j,k,id_sbef_)   &
                                ) & !mmm>
                                      ) & !pmx>

       ) !prime>


      enddo !k
      enddo !i
      enddo !j

! Partial Oi Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax


! Cumul du pivot obligatoirement avant advection partielle 06-02-17
! Note que pour economiser un peu de cpu la multiplication par wetmask_t(i,j)/dz_t(i,j,k,before),
! commune a toutes les operations, sera appliquee une fois pour toute a la fin...
       anyv3d(i,j,k,id_sdiv)= &
       anyv3d(i,j,k,id_sdiv)  &
      +anyv3d(i,j,k,id_sbef_)*dti_lpsub*(veldydz_u(i+1,j,k,1)           &
                                        -veldydz_u(i  ,j,k,1))/dxdy_t(i,j)

        anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !ooo>
        anyv3d(i,j,k,id_sbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i+1,j,k,id_flx_)                &
             -anyv3d(i  ,j,k,id_flx_)                &

         +(veldydz_u(i+1,j,k,1)                      &
          -veldydz_u(i  ,j,k,1))*anyv3d(i,j,k,id_sbef_) & !25-01-17

         )/dxdy_t(i,j)  ) & !ooo>

         +(1-mask_t(i,j,k))*sal_t(i,j,k,before)*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

      enddo
      enddo
      enddo

!     call my_outputs_zonesalttempflux('is',id_flx_) !29-01-19

      call vertmix_merged_levels_any(0,id_sbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_sbef_)
      call obc_int_anyv3d(id_sbef_,'zb')

! Flux Oj
      cst_=0.5*cst_adv_hor/3.
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax+1
      do i=1,imax

      anyv3d(i,j,k,id_fly_)=0.5*( & !prime>

! Schema QUICK en dehors influence des fleuves
      anyv3d(i,j,k,id_hybcoefv)*( &
! Partie advective QUICK:
      -veldxdz_v(i,j,k,1)  &
!       *(anyv3d(i,j,k,id_tnow_)+anyv3d(i,j-1,k,id_tnow_)) &
! Si le nombre de courant est faible il n'y a pas de risque A maintenir la
! partie centree sur sal_t(now) ce qui permet de maintenir un leap-frog centrE
! classique meme en subdivisant le pas de salps ce qui permet A priori de 
! maintenir les proprietes de conservation associeE
        *(    sal_t(i,j,k,now)+sal_t(i,j-1,k,now)  &
#ifdef bidonref
       +2.*salref_v(i,j,k)                         &
          -salref_t(i,j,k) -salref_t(i,j-1,k)      &
#endif
                                                 ) &

! Partie diffusive QUICK:
       +cst_*0.5*(anyv3d(i,j,k,id_ofactors)+anyv3d(i+1,j,k,id_ofactors))  &

! Efficacite densitaire
!       *abs(-alp_t*(   tem_t(i,j,k,1)-tem_t(i,j-1,k,1))+alp_s*(   sal_t(i,j,k,1)-sal_t(i,j-1,k,1))) &
!       /   (+alp_t*abs(tem_t(i,j,k,1)-tem_t(i,j-1,k,1))+alp_s*abs(sal_t(i,j,k,1)-sal_t(i,j-1,k,1))+small1) &

                                                                  *( & !m°v°m>

             (veldxdz_v(i,j,k,1)+abs(veldxdz_v(i,j,k,1)))*(  &       !ooo>

            (   anyv3d(i,j  ,k,id_sbef_)-anyv3d(i,j-1,k,id_sbef_)  &

! correction anyv3d(i,j  ,k,id_sbef_) par rapport niveau z(i,j-1,k)
      +( anyv3d(i,j  ,kp1,id_sbef_)- anyv3d(i,j  ,km1,id_sbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j-1,k  )         -depth_t(i,j  ,k  )         )    &
      /(depth_t(i,j  ,kp1)         -depth_t(i,j  ,km1)         )    &
          )       )      & !minmax>

            ) &

           -(   anyv3d(i,j-1,k,id_sbef_)-anyv3d(i,j-2,k,id_sbef_)  &

! correction anyv3d(i,j-2,k,id_sbef_) par rapport niveau z(i,j-1,k)
      -( anyv3d(i,j-2,kp1,id_sbef_)- anyv3d(i,j-2,km1,id_sbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j-1,k  )         -depth_t(i,j-2,k  )         )    &
      /(depth_t(i,j-2,kp1)         -depth_t(i,j-2,km1)         )    &
          )       )      & !minmax>

            ) &
                                        *mask_t(i,j-2,kmax)       )& !ooo>

            +(veldxdz_v(i,j,k,1)-abs(veldxdz_v(i,j,k,1)))*(  &       !ooo>

            (   anyv3d(i,j+1,k,id_sbef_)-anyv3d(i,j,k,id_sbef_)    &

! correction anyv3d(i,j+1,k,id_sbef_) par rapport niveau z(i,j  ,k)
      +( anyv3d(i,j+1,kp1,id_sbef_)- anyv3d(i,j+1,km1,id_sbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j  ,k  )         -depth_t(i,j+1,k  )         )    &
      /(depth_t(i,j+1,kp1)         -depth_t(i,j+1,km1)         )    &
          )       )      & !minmax>

            )   &
               *mask_t(i,j+1,kmax)                                 &

           -(   anyv3d(i,j,k,id_sbef_)-anyv3d(i,j-1,k,id_sbef_)    &

! correction anyv3d(i,j-1,k,id_sbef_) par rapport niveau z(i,j  ,k)
      -( anyv3d(i,j-1,kp1,id_sbef_)- anyv3d(i,j-1,km1,id_sbef_))    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j  ,k  )         -depth_t(i,j-1,k  )         )    &
      /(depth_t(i,j-1,kp1)         -depth_t(i,j-1,km1)         )    &
          )       )      & !minmax>

                                                               ))  & !ooo>

                                                                   )  & !m°v°m>()cst_
                                                                   )  & !()anyv3d(i,j,k,id_hybcoefv)
! Schema UPWIND dans les zones sous influence des fleuves:
      +(1.-anyv3d(i,j,k,id_hybcoefv))*( & !pmx>
          -veldxdz_v(i,j,k,1)* ( & !ppp>
              anyv3d(i,j,k,id_sbef_)+anyv3d(i,j-1,k,id_sbef_)   &
                               ) & !ppp>
      +abs(veldxdz_v(i,j,k,1))*( & !mmm>
              anyv3d(i,j,k,id_sbef_)-anyv3d(i,j-1,k,id_sbef_)   &
                               ) & !mmm>
                                      ) & !pmx>

       ) !prime>

      enddo !k
      enddo !i
      enddo !j

! Partial Horizontal Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_sdiv)= &
       anyv3d(i,j,k,id_sdiv)  &
      +anyv3d(i,j,k,id_sbef_)*dti_lpsub*(veldxdz_v(i,j+1,k,1)   &
                                        -veldxdz_v(i,j  ,k,1))/dxdy_t(i,j)

        anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !ooo>
        anyv3d(i,j,k,id_sbef_)*dz_t(i,j,k,before)   & !01-05-19    

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i,j+1,k,id_fly_)                &
             -anyv3d(i,j  ,k,id_fly_)                &

         +(veldxdz_v(i,j+1,k,1)                      &
          -veldxdz_v(i,j  ,k,1))*anyv3d(i,j,k,id_sbef_) & !25-01-17

         )/dxdy_t(i,j)  ) & !ooo>

         +(1-mask_t(i,j,k))*sal_t(i,j,k,before)*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

      enddo
      enddo
      enddo

!     call my_outputs_zonesalttempflux('js',id_fly_) 
      call vertmix_merged_levels_any(0,id_sbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_sbef_)
      call obc_int_anyv3d(id_sbef_,'zb')

!#endif
      enddo               !*********>

      do k=1,kmax ; do j=1,jmax  ; do i=1,imax  

           sal_t(i,j,k,2)                              &
       =  anyv3d(i,j,k,id_sbef_)                    !  &

!     -wetmask_t(i,j)*sal_t(i,j,k,now)*dti_lp*(        &

!                                              veldydz_u(i+1,j,k,1)   &
!                                             -veldydz_u(i  ,j,k,1)   &
!                                             +veldxdz_v(i,j+1,k,1)   &
!                                             -veldxdz_v(i,j  ,k,1)   &

!                                             )/(dz_t(i,j,k,before)*dxdy_t(i,j))

      enddo ; enddo ; enddo

      end subroutine advection_scal_uplim_subs
