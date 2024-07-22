










      subroutine advection_scal
!______________________________________________________________________
! SYMPHONIE ocean model
! release 365 - last update: 26-01-23
!______________________________________________________________________

      use module_principal
      use module_parallele
      implicit none

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
!         14-08-20 obc_int_anyv3d renomme obc_mpi_anyv3d
! v289    24-09-20 my_outputs_zonesalttempflux renommE my_outputs_zone1salttempflux
! v296    25-02-21 nouvelle subroutine advection_scal_ofactor !25-02-21
! v297    05-03-21 rivieres de surface et de fond integrees au schema d'advection
! v299    11-03-21 Ajout d'une variante plus proche de l'algo quickest "officiel" (iadvec_ts=2)
!                  La version precedente est maintenue (subroutines xxxx_odl) avec iadvec_ts=1 
! v300    22-03-21 Suite point precedent
! v301    29-04-21 x2_r4=x1_r4*ratio_negdif_ver !29-04-21
! v302    05-05-21 ajout neg dif sur la composante upwind
!         09-05-21 nouvelles fonctions temf_t et salf_t: filtre 3 points
! v303    31-07-21 conditions limites surf et fond des champs "negdif"
! v305    23-08-21 call advection_scal_vertical_laxwen      !23-08-21
! v309    26-08-21 subroutine advection_scal_lw_tem         !26-08-21
!         30-09-21 subroutine advection_scal_negdif_hori 
!                  subroutine advection_scal_negdif_vert
!         04-10-21 dif neg: distinguer directions horizontales et verticale
! v310    02-11-21 Appel A vertmix_merged_levels_any supprimE 
! v314    29-11-21 correction du test if(iadvec_ts_hor==iadvec_ts_hor_quickest2)
! v322    24-01-22 correction bug: ratio_negdif_hor remplace ratio_negdif_ver
! v325    12-02-22 ajout call my_outputs_3dflux_s
! v330    24-02-22 amelioration du limiteur
! v331    27-02-22 amelioration du limiteur avec int(wetmask...)
! v332    28-02-22 *upwindwetdry_t(i,j)  wet/dry zone 
! v334    03-03-22 limiteur C0: CL sur rhp_t
! v336    08-03-22 - Detection de presence de panache (si dS/di > 5 upwind) 
!                  - limiteur detectant inversion de gradient  
! v350    18-07-22 - possibilite de quickest-2 pour l'advection verticale 
!                  et suppression (temporaire?) du schema c2
!                  - mise A jour du filtre 5 point sur la verticale
!         20-07-22 if(quick_filter_points==
! v356    06-10-22 cas quick_filter_points=0
!         07-10-22 si temf_t et salf_t ont ete utilises par l'advection horizontale retablir zero
! v361    31-12-22 advection horizontale implicite
! v362    04-01-23 work_status_mpi_min_
! v365    26-01-23 modif echange mpi sur densite depuis qu'elle est definie sur 0:imax+1, 0:jmax+1 !26-01-23
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


      end subroutine advection_scal_rmv_ref

!________________________________________________________________________________

      subroutine advection_scal_rst_ref
      use module_principal
      implicit none

      if(flag_refstate==0) return !11-11-15


      end subroutine advection_scal_rst_ref

!________________________________________________________________________________
      subroutine advection_scal_up3  ! 11-08-10
      use module_principal
      use module_parallele
      implicit none
      stop 'je ne dois plus passer par advection_scal_up3'
      end subroutine advection_scal_up3
!________________________________________________________________________________

      subroutine advection_scal_c2_ref
      use module_principal
      use module_parallele
      implicit none

! 
      if(flag_refstate==0) return !11-11-15

      end subroutine advection_scal_c2_ref
   
!________________________________________________________________________________

      subroutine advection_scal_zaxis
      use module_principal
      use module_parallele
      implicit none
      double precision :: grdamp_=2. ! lower gradient amplification

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

      subroutine advection_scal_vbc
      use module_principal
      implicit none

! Completement reecrit le 05-03-21
       tem_t(:,:,kmax+1,2)=tem_t(:,:,kmax,2)
       tem_t(:,:,0,2)=tem_t(:,:,1,2)
! Cette condition S(0),S(kmax+1) etait non seulement inutile car advection_scal_vertical
! force le flux vertical de sel en k=1 et k=kmax+1 A etre nul (la vitesse verticale est 
! 100% explicite, le schema 100% upwind et anyv1d(kmax+1,4)=anyv1d(1,4)=0) mais en plus
! c'est nuisible A la continuite de la derivee seconde du quickest en k=2 et k=kmax
!      sal_t(:,:,kmax+1,2)=0.
!      sal_t(:,:,0,2)=0.
! On revient par consequent A la mEme condition que pour la temperature !29-11-21
       sal_t(:,:,kmax+1,2)=sal_t(:,:,kmax,2)                            !29-11-21
       sal_t(:,:,0,2)     =sal_t(:,:,1,2)                               !29-11-21

      do kr=1,nriver
       if(riverdir(kr)==0) then  !m°v°m> 
        if(rivertrc_inout(kr)==1)     then !rrrrrrrrrrrrrrrrrr>  
                 i=iriver(kr,1)
                 j=jriver(kr,1)
                 tem_t(i,j,kmax+1,2)=river_t(kr,1)
        endif                              !rrrrrrrrrrrrrrrrrr>
       endif                     !m°v°m> 
       if(riverdir(kr)==-1) then !m°v°m> 
        if(rivertrc_inout(kr)==1)     then !rrrrrrrrrrrrrrrrrr>  
                 i=iriver(kr,1)
                 j=jriver(kr,1)
                 tem_t(i,j,0,2)=river_t(kr,1) !29-11-21
        endif                              !rrrrrrrrrrrrrrrrrr>
       endif                     !m°v°m> 
      enddo ! kr loop

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
!     -veldydz_u(i,j,k,id_veltot)*(tem_t(i,j,k,t_)+tem_t(i-1 ,j,k,t_))         &
      -veldydz_u(i,j,k,id_veltot)*                                             &
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
           -veldydz_u(i,j,k,id_veltot)* (tem_t(i,j,k,0)+tem_t(i-1,j,k,0))  &
       +abs(veldydz_u(i,j,k,id_veltot))*(tem_t(i,j,k,0)-tem_t(i-1,j,k,0))) &

                                )   !()const1

! Flux S:
!     velef_=0.5*(anyv3d(i,j,k,7)+anyv3d(i,j+1,k,7))*(sal_t(i,j,k,1)-sal_t(i-1,j,k,1))/dx_u(i,j) &
!               *dy_u(i,j)*dz_u(i,j,k,1)*mask_u(i,j,kmax)
      velef_=scalefac_*(anyv3d(i,j,k,7)+anyv3d(i,j+1,k,7))    &
                      *( sal_t(i,j,k,1)- sal_t(i-1,j,k,1)     & 
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
!      -veldydz_u(i,j,k,id_veltot)*(sal_t(i,j,k,t_)+sal_t(i-1 ,j,k,t_))         &
      -veldydz_u(i,j,k,id_veltot)*                                             &
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
           -veldydz_u(i,j,k,id_veltot)* (sal_t(i,j,k,0)+sal_t(i-1,j,k,0))  &
       +abs(veldydz_u(i,j,k,id_veltot))*(sal_t(i,j,k,0)-sal_t(i-1,j,k,0))) &

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
!     -veldxdz_v(i,j,k,id_veltot)*(tem_t(i,j,k,t_)+tem_t(i,j-1,k,t_))          &
      -veldxdz_v(i,j,k,id_veltot)*                                             &
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
          -veldxdz_v(i,j,k,id_veltot)* (tem_t(i,j,k,0)+tem_t(i,j-1,k,0))      &
      +abs(veldxdz_v(i,j,k,id_veltot))*(tem_t(i,j,k,0)-tem_t(i,j-1,k,0)))     &

                                ) !()const1

! Flux S:
!     velef_=0.5*(anyv3d(i,j,k,7)+anyv3d(i+1,j,k,7))*(sal_t(i,j,k,1)-sal_t(i,j-1,k,1))/dy_v(i,j) &
!              *dx_v(i,j)*dz_v(i,j,k,1)*mask_v(i,j,kmax)
      velef_=scalefac_*(anyv3d(i,j,k,7)+anyv3d(i+1,j,k,7))   &
                      *( sal_t(i,j,k,1)- sal_t(i,j-1,k,1)    &
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
!     -veldxdz_v(i,j,k,id_veltot)*(sal_t(i,j,k,t_)+sal_t(i,j-1,k,t_))          &
      -veldxdz_v(i,j,k,id_veltot)*                                             &
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
          -veldxdz_v(i,j,k,id_veltot)* (sal_t(i,j,k,0)+sal_t(i,j-1,k,0))      &
      +abs(veldxdz_v(i,j,k,id_veltot))*(sal_t(i,j,k,0)-sal_t(i,j-1,k,0)))     &

                                ) !()const1

!      if(i+par%timax(1)==30.and.j+par%tjmax(1)==3.and.k==kmax) then
!       write(6,*)'mask x5 veldxdz',mask_v(i,j,k),real(x5),real(veldxdz_v(i,j,k,id_veltot)),real(sal_t(i,j-1,k,0))
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

      stop 'Je ne dois plus passer par advection_scal_diffactor'
      end subroutine advection_scal_diffactor

!................................................................

      subroutine advection_scal_ofactor !25-02-21
      use module_principal
      use module_parallele
      implicit none

      if(flag_negdif_hor==1)return !car les subroutines concernees n'utilisent pas anyv3d(:,:,:,id_ofactor

! Dans le cas ou le critere d'orthogonalite du courant n'est pas applique (flag_ofactor==0)
! alors fixer les tableaux anyv3d(:,:,:,id_ofactort) anyv3d(:,:,:,id_ofactors) a 1 et sortir:
      if(flag_ofactor==0) then !m°v°m> !30-09-21
         anyv3d(:,:,:,id_ofactort)=1.
         anyv3d(:,:,:,id_ofactors)=1.
       RETURN
      endif                    !m°v°m>

! ajouts corrections sigma, moyenne sur 2 temps, inversion du numerateur 
! pour favoriser la valeur 1 en situations ambigues !25-02-21

! Commentaires dans:
! https://docs.google.com/document/d/1-4NdFRKCug7u4-reURBVsHOeHMtf5p_8RlCYDLwAZjA/edit

      do k=1,kmax
         kp1=min(k+1,kmax)
         km1=max(k-1,1)
       do j=2,jmax
       do i=2,imax


         anyv3d(i,j,k,id_ofactort)=                             &
        1.-mask_f(i,j,k)*     &
        abs( & !ABS> abs(-v.dT/dx+u.dT/dy))
           -(veldxdz_v(i,j  ,k,id_veltot)*invdx_v(i,j)+veldxdz_v(i,j-1,k,id_veltot)*invdx_v(i,j-1)) &
      *( & !PMX>     
         (tem_t(i  ,j,k,1)+tem_t(i  ,j-1,k,1)+tem_t(i  ,j,k,0)+tem_t(i  ,j-1,k,0)) & !  tem_v(i  ,j) 
        -(tem_t(i-1,j,k,1)+tem_t(i-1,j-1,k,1)+tem_t(i-1,j,k,0)+tem_t(i-1,j-1,k,0)) & ! -tem_v(i-1,j) 
        -(depth_v(i,j,k)-depth_v(i-1,j,k))*( & !ooo>                                 ! -(z_v(i,j,k)-z_v(i-1,j,k))*[
         (tem_t(i  ,j,kp1,1)+tem_t(i  ,j-1,kp1,1)+tem_t(i  ,j,kp1,0)+tem_t(i  ,j-1,kp1,0)) & !  tem_v(i  ,j,kp1)
        +(tem_t(i-1,j,kp1,1)+tem_t(i-1,j-1,kp1,1)+tem_t(i-1,j,kp1,0)+tem_t(i-1,j-1,kp1,0)) & ! +tem_v(i-1,j,kp1)
        -(tem_t(i  ,j,km1,1)+tem_t(i  ,j-1,km1,1)+tem_t(i  ,j,km1,0)+tem_t(i  ,j-1,km1,0)) & ! -tem_v(i  ,j,km1)
        -(tem_t(i-1,j,km1,1)+tem_t(i-1,j-1,km1,1)+tem_t(i-1,j,km1,0)+tem_t(i-1,j-1,km1,0)) & ! -tem_v(i-1,j,km1)
                                           ) & !ooo>                                 ! ]
        /(depth_v(i,j,kp1)+depth_v(i-1,j,kp1)-depth_v(i,j,km1)-depth_v(i-1,j,km1)) & ! /(   z_v(i,j,kp1)+  z_v(i-1,j,kp1)  -z_v(i,j,km1)  -z_v(i-1,j,km1))
          ) & !PMX>
           *invdx_f(i,j)  &

           +(veldydz_u(i,j  ,k,id_veltot)*invdy_u(i,j)+veldydz_u(i,j-1,k,id_veltot)*invdy_u(i,j-1)) &
      *( & !PMX>     
         (tem_t(i,j  ,k,1)+tem_t(i-1,j  ,k,1)+tem_t(i,j  ,k,0)+tem_t(i-1,j  ,k,0)) & !  tem_u(i,j) 
        -(tem_t(i,j-1,k,1)+tem_t(i-1,j-1,k,1)+tem_t(i,j-1,k,0)+tem_t(i-1,j-1,k,0)) & !  tem_u(i,j-1) 
        -(depth_u(i,j,k)-depth_u(i,j-1,k))*( & !ooo>                                 ! -(z_u(i,j)-z_u(i,j-1)*[
         (tem_t(i,j  ,kp1,1)+tem_t(i-1,j  ,kp1,1)+tem_t(i,j  ,kp1,0)+tem_t(i-1,j  ,kp1,0)) & !  tem_u(i,j  ,kp1)
        +(tem_t(i,j-1,kp1,1)+tem_t(i-1,j-1,kp1,1)+tem_t(i,j-1,kp1,0)+tem_t(i-1,j-1,kp1,0)) & ! +tem_u(i,j-1,kp1)
        -(tem_t(i,j  ,km1,1)+tem_t(i-1,j  ,km1,1)+tem_t(i,j  ,km1,0)+tem_t(i-1,j  ,km1,0)) & ! -tem_u(i,j  ,km1)
        -(tem_t(i,j-1,km1,1)+tem_t(i-1,j-1,km1,1)+tem_t(i,j-1,km1,0)+tem_t(i-1,j-1,km1,0)) & ! -tem_u(i,j-1,km1)
                                           ) & !ooo>                                 ! ]
        /(depth_u(i,j,kp1)+depth_u(i,j-1,kp1)-depth_u(i,j,km1)-depth_u(i,j-1,km1)) & !   /(   z_u(i,j,kp1)  +z_u(i,j-1,kp1)  -z_u(i,j,km1)  -z_u(i,j-1,km1))
           ) & !PMX>
           *invdy_f(i,j) &
            ) & !ABS> abs(-v.dT/dx+u.dT/dy))
           /max( & ! PARENTHESE 1
        abs( & !ABS> abs(-v.dT/dx+u.dT/dy))
           -(veldxdz_v(i,j  ,k,id_veltot)*invdx_v(i,j)+veldxdz_v(i,j-1,k,id_veltot)*invdx_v(i,j-1)) &
      *( & !PMX>     
         (tem_t(i  ,j,k,1)+tem_t(i  ,j-1,k,1)+tem_t(i  ,j,k,0)+tem_t(i  ,j-1,k,0)) & !  tem_v(i  ,j) 
        -(tem_t(i-1,j,k,1)+tem_t(i-1,j-1,k,1)+tem_t(i-1,j,k,0)+tem_t(i-1,j-1,k,0)) & ! -tem_v(i-1,j) 
        -(depth_v(i,j,k)-depth_v(i-1,j,k))*( & !ooo>                                 ! -(z_v(i,j,k)-z_v(i-1,j,k))*[
         (tem_t(i  ,j,kp1,1)+tem_t(i  ,j-1,kp1,1)+tem_t(i  ,j,kp1,0)+tem_t(i  ,j-1,kp1,0)) & !  tem_v(i  ,j,kp1)
        +(tem_t(i-1,j,kp1,1)+tem_t(i-1,j-1,kp1,1)+tem_t(i-1,j,kp1,0)+tem_t(i-1,j-1,kp1,0)) & ! +tem_v(i-1,j,kp1)
        -(tem_t(i  ,j,km1,1)+tem_t(i  ,j-1,km1,1)+tem_t(i  ,j,km1,0)+tem_t(i  ,j-1,km1,0)) & ! -tem_v(i  ,j,km1)
        -(tem_t(i-1,j,km1,1)+tem_t(i-1,j-1,km1,1)+tem_t(i-1,j,km1,0)+tem_t(i-1,j-1,km1,0)) & ! -tem_v(i-1,j,km1)
                                           ) & !ooo>                                 ! ]
        /(depth_v(i,j,kp1)+depth_v(i-1,j,kp1)-depth_v(i,j,km1)-depth_v(i-1,j,km1)) & ! /(   z_v(i,j,kp1)+  z_v(i-1,j,kp1)  -z_v(i,j,km1)  -z_v(i-1,j,km1))
          ) & !PMX>
           *invdx_f(i,j)  &

           +(veldydz_u(i,j  ,k,id_veltot)*invdy_u(i,j)+veldydz_u(i,j-1,k,id_veltot)*invdy_u(i,j-1)) &
      *( & !PMX>     
         (tem_t(i,j  ,k,1)+tem_t(i-1,j  ,k,1)+tem_t(i,j  ,k,0)+tem_t(i-1,j  ,k,0)) & !  tem_u(i,j) 
        -(tem_t(i,j-1,k,1)+tem_t(i-1,j-1,k,1)+tem_t(i,j-1,k,0)+tem_t(i-1,j-1,k,0)) & !  tem_u(i,j-1) 
        -(depth_u(i,j,k)-depth_u(i,j-1,k))*( & !ooo>                                 ! -(z_u(i,j)-z_u(i,j-1)*[
         (tem_t(i,j  ,kp1,1)+tem_t(i-1,j  ,kp1,1)+tem_t(i,j  ,kp1,0)+tem_t(i-1,j  ,kp1,0)) & !  tem_u(i,j  ,kp1)
        +(tem_t(i,j-1,kp1,1)+tem_t(i-1,j-1,kp1,1)+tem_t(i,j-1,kp1,0)+tem_t(i-1,j-1,kp1,0)) & ! +tem_u(i,j-1,kp1)
        -(tem_t(i,j  ,km1,1)+tem_t(i-1,j  ,km1,1)+tem_t(i,j  ,km1,0)+tem_t(i-1,j  ,km1,0)) & ! -tem_u(i,j  ,km1)
        -(tem_t(i,j-1,km1,1)+tem_t(i-1,j-1,km1,1)+tem_t(i,j-1,km1,0)+tem_t(i-1,j-1,km1,0)) & ! -tem_u(i,j-1,km1)
                                           ) & !ooo>                                 ! ]
        /(depth_u(i,j,kp1)+depth_u(i,j-1,kp1)-depth_u(i,j,km1)-depth_u(i,j-1,km1)) & !   /(   z_u(i,j,kp1)  +z_u(i,j-1,kp1)  -z_u(i,j,km1)  -z_u(i,j-1,km1))
           ) & !PMX>
           *invdy_f(i,j) &
            ) & !ABS> abs(-v.dT/dx+u.dT/dy))

       +abs( & !ABS> abs( u.dT/dx+v.dT/dy)
            (veldydz_u(i,j  ,k,id_veltot)*invdy_u(i,j)+veldydz_u(i,j-1,k,id_veltot)*invdy_u(i,j-1)) &
      *( & !PMX>     
         (tem_t(i  ,j,k,1)+tem_t(i  ,j-1,k,1)+tem_t(i  ,j,k,0)+tem_t(i  ,j-1,k,0)) & !  tem_v(i  ,j) 
        -(tem_t(i-1,j,k,1)+tem_t(i-1,j-1,k,1)+tem_t(i-1,j,k,0)+tem_t(i-1,j-1,k,0)) & ! -tem_v(i-1,j) 
        -(depth_v(i,j,k)-depth_v(i-1,j,k))*( & !ooo>                                 ! -(z_v(i,j,k)-z_v(i-1,j,k))*[
         (tem_t(i  ,j,kp1,1)+tem_t(i  ,j-1,kp1,1)+tem_t(i  ,j,kp1,0)+tem_t(i  ,j-1,kp1,0)) & !  tem_v(i  ,j,kp1)
        +(tem_t(i-1,j,kp1,1)+tem_t(i-1,j-1,kp1,1)+tem_t(i-1,j,kp1,0)+tem_t(i-1,j-1,kp1,0)) & ! +tem_v(i-1,j,kp1)
        -(tem_t(i  ,j,km1,1)+tem_t(i  ,j-1,km1,1)+tem_t(i  ,j,km1,0)+tem_t(i  ,j-1,km1,0)) & ! -tem_v(i  ,j,km1)
        -(tem_t(i-1,j,km1,1)+tem_t(i-1,j-1,km1,1)+tem_t(i-1,j,km1,0)+tem_t(i-1,j-1,km1,0)) & ! -tem_v(i-1,j,km1)
                                           ) & !ooo>                                 ! ]
        /(depth_v(i,j,kp1)+depth_v(i-1,j,kp1)-depth_v(i,j,km1)-depth_v(i-1,j,km1)) & ! /(   z_v(i,j,kp1)+  z_v(i-1,j,kp1)  -z_v(i,j,km1)  -z_v(i-1,j,km1))
          ) & !PMX>
           *invdx_f(i,j)  &

           +(veldxdz_v(i,j  ,k,id_veltot)*invdx_v(i,j)+veldxdz_v(i,j-1,k,id_veltot)*invdx_v(i,j-1)) &
      *( & !PMX>     
         (tem_t(i,j  ,k,1)+tem_t(i-1,j  ,k,1)+tem_t(i,j  ,k,0)+tem_t(i-1,j  ,k,0)) & !  tem_u(i,j) 
        -(tem_t(i,j-1,k,1)+tem_t(i-1,j-1,k,1)+tem_t(i,j-1,k,0)+tem_t(i-1,j-1,k,0)) & !  tem_u(i,j-1) 
        -(depth_u(i,j,k)-depth_u(i,j-1,k))*( & !ooo>                                 ! -(z_u(i,j)-z_u(i,j-1)*[
         (tem_t(i,j  ,kp1,1)+tem_t(i-1,j  ,kp1,1)+tem_t(i,j  ,kp1,0)+tem_t(i-1,j  ,kp1,0)) & !  tem_u(i,j  ,kp1)
        +(tem_t(i,j-1,kp1,1)+tem_t(i-1,j-1,kp1,1)+tem_t(i,j-1,kp1,0)+tem_t(i-1,j-1,kp1,0)) & ! +tem_u(i,j-1,kp1)
        -(tem_t(i,j  ,km1,1)+tem_t(i-1,j  ,km1,1)+tem_t(i,j  ,km1,0)+tem_t(i-1,j  ,km1,0)) & ! -tem_u(i,j  ,km1)
        -(tem_t(i,j-1,km1,1)+tem_t(i-1,j-1,km1,1)+tem_t(i,j-1,km1,0)+tem_t(i-1,j-1,km1,0)) & ! -tem_u(i,j-1,km1)
                                           ) & !ooo>                                 ! ]
        /(depth_u(i,j,kp1)+depth_u(i,j-1,kp1)-depth_u(i,j,km1)-depth_u(i,j-1,km1)) & !   /(   z_u(i,j,kp1)  +z_u(i,j-1,kp1)  -z_u(i,j,km1)  -z_u(i,j-1,km1))
           ) & !PMX>
           *invdy_f(i,j) &
            ) & !ABS> abs( u.dT/dx+v.dT/dy)

            ,small2)   ! PARENTHESE 1

         anyv3d(i,j,k,id_ofactors)=                             &
        1.-mask_f(i,j,k)*     &
        abs( & !ABS> abs(-v.dS/dx+u.dS/dy))
           -(veldxdz_v(i,j  ,k,id_veltot)*invdx_v(i,j)+veldxdz_v(i,j-1,k,id_veltot)*invdx_v(i,j-1)) &
      *( & !PMX>     
         (sal_t(i  ,j,k,1)+sal_t(i  ,j-1,k,1)+sal_t(i  ,j,k,0)+sal_t(i  ,j-1,k,0)) & !  sal_v(i  ,j) 
        -(sal_t(i-1,j,k,1)+sal_t(i-1,j-1,k,1)+sal_t(i-1,j,k,0)+sal_t(i-1,j-1,k,0)) & ! -sal_v(i-1,j) 
        -(depth_v(i,j,k)-depth_v(i-1,j,k))*( & !ooo>                                 ! -(z_v(i,j,k)-z_v(i-1,j,k))*[
         (sal_t(i  ,j,kp1,1)+sal_t(i  ,j-1,kp1,1)+sal_t(i  ,j,kp1,0)+sal_t(i  ,j-1,kp1,0)) & !  sal_v(i  ,j,kp1)
        +(sal_t(i-1,j,kp1,1)+sal_t(i-1,j-1,kp1,1)+sal_t(i-1,j,kp1,0)+sal_t(i-1,j-1,kp1,0)) & ! +sal_v(i-1,j,kp1)
        -(sal_t(i  ,j,km1,1)+sal_t(i  ,j-1,km1,1)+sal_t(i  ,j,km1,0)+sal_t(i  ,j-1,km1,0)) & ! -sal_v(i  ,j,km1)
        -(sal_t(i-1,j,km1,1)+sal_t(i-1,j-1,km1,1)+sal_t(i-1,j,km1,0)+sal_t(i-1,j-1,km1,0)) & ! -sal_v(i-1,j,km1)
                                           ) & !ooo>                                 ! ]
        /(depth_v(i,j,kp1)+depth_v(i-1,j,kp1)-depth_v(i,j,km1)-depth_v(i-1,j,km1)) & ! /(   z_v(i,j,kp1)+  z_v(i-1,j,kp1)  -z_v(i,j,km1)  -z_v(i-1,j,km1))
          ) & !PMX>
           *invdx_f(i,j)  &

           +(veldydz_u(i,j  ,k,id_veltot)*invdy_u(i,j)+veldydz_u(i,j-1,k,id_veltot)*invdy_u(i,j-1)) &
      *( & !PMX>     
         (sal_t(i,j  ,k,1)+sal_t(i-1,j  ,k,1)+sal_t(i,j  ,k,0)+sal_t(i-1,j  ,k,0)) & !  sal_u(i,j) 
        -(sal_t(i,j-1,k,1)+sal_t(i-1,j-1,k,1)+sal_t(i,j-1,k,0)+sal_t(i-1,j-1,k,0)) & !  sal_u(i,j-1) 
        -(depth_u(i,j,k)-depth_u(i,j-1,k))*( & !ooo>                                 ! -(z_u(i,j)-z_u(i,j-1)*[
         (sal_t(i,j  ,kp1,1)+sal_t(i-1,j  ,kp1,1)+sal_t(i,j  ,kp1,0)+sal_t(i-1,j  ,kp1,0)) & !  sal_u(i,j  ,kp1)
        +(sal_t(i,j-1,kp1,1)+sal_t(i-1,j-1,kp1,1)+sal_t(i,j-1,kp1,0)+sal_t(i-1,j-1,kp1,0)) & ! +sal_u(i,j-1,kp1)
        -(sal_t(i,j  ,km1,1)+sal_t(i-1,j  ,km1,1)+sal_t(i,j  ,km1,0)+sal_t(i-1,j  ,km1,0)) & ! -sal_u(i,j  ,km1)
        -(sal_t(i,j-1,km1,1)+sal_t(i-1,j-1,km1,1)+sal_t(i,j-1,km1,0)+sal_t(i-1,j-1,km1,0)) & ! -sal_u(i,j-1,km1)
                                           ) & !ooo>                                 ! ]
        /(depth_u(i,j,kp1)+depth_u(i,j-1,kp1)-depth_u(i,j,km1)-depth_u(i,j-1,km1)) & !   /(   z_u(i,j,kp1)  +z_u(i,j-1,kp1)  -z_u(i,j,km1)  -z_u(i,j-1,km1))
           ) & !PMX>
           *invdy_f(i,j) &
            ) & !ABS> abs(-v.dS/dx+u.dS/dy))
           /max( & ! PARENTHESE 1
        abs( & !ABS> abs(-v.dS/dx+u.dS/dy))
           -(veldxdz_v(i,j  ,k,id_veltot)*invdx_v(i,j)+veldxdz_v(i,j-1,k,id_veltot)*invdx_v(i,j-1)) &
      *( & !PMX>     
         (sal_t(i  ,j,k,1)+sal_t(i  ,j-1,k,1)+sal_t(i  ,j,k,0)+sal_t(i  ,j-1,k,0)) & !  sal_v(i  ,j) 
        -(sal_t(i-1,j,k,1)+sal_t(i-1,j-1,k,1)+sal_t(i-1,j,k,0)+sal_t(i-1,j-1,k,0)) & ! -sal_v(i-1,j) 
        -(depth_v(i,j,k)-depth_v(i-1,j,k))*( & !ooo>                                 ! -(z_v(i,j,k)-z_v(i-1,j,k))*[
         (sal_t(i  ,j,kp1,1)+sal_t(i  ,j-1,kp1,1)+sal_t(i  ,j,kp1,0)+sal_t(i  ,j-1,kp1,0)) & !  sal_v(i  ,j,kp1)
        +(sal_t(i-1,j,kp1,1)+sal_t(i-1,j-1,kp1,1)+sal_t(i-1,j,kp1,0)+sal_t(i-1,j-1,kp1,0)) & ! +sal_v(i-1,j,kp1)
        -(sal_t(i  ,j,km1,1)+sal_t(i  ,j-1,km1,1)+sal_t(i  ,j,km1,0)+sal_t(i  ,j-1,km1,0)) & ! -sal_v(i  ,j,km1)
        -(sal_t(i-1,j,km1,1)+sal_t(i-1,j-1,km1,1)+sal_t(i-1,j,km1,0)+sal_t(i-1,j-1,km1,0)) & ! -sal_v(i-1,j,km1)
                                           ) & !ooo>                                 ! ]
        /(depth_v(i,j,kp1)+depth_v(i-1,j,kp1)-depth_v(i,j,km1)-depth_v(i-1,j,km1)) & ! /(   z_v(i,j,kp1)+  z_v(i-1,j,kp1)  -z_v(i,j,km1)  -z_v(i-1,j,km1))
          ) & !PMX>
           *invdx_f(i,j)  &

           +(veldydz_u(i,j  ,k,id_veltot)*invdy_u(i,j)+veldydz_u(i,j-1,k,id_veltot)*invdy_u(i,j-1)) &
      *( & !PMX>     
         (sal_t(i,j  ,k,1)+sal_t(i-1,j  ,k,1)+sal_t(i,j  ,k,0)+sal_t(i-1,j  ,k,0)) & !  sal_u(i,j) 
        -(sal_t(i,j-1,k,1)+sal_t(i-1,j-1,k,1)+sal_t(i,j-1,k,0)+sal_t(i-1,j-1,k,0)) & !  sal_u(i,j-1) 
        -(depth_u(i,j,k)-depth_u(i,j-1,k))*( & !ooo>                                 ! -(z_u(i,j)-z_u(i,j-1)*[
         (sal_t(i,j  ,kp1,1)+sal_t(i-1,j  ,kp1,1)+sal_t(i,j  ,kp1,0)+sal_t(i-1,j  ,kp1,0)) & !  sal_u(i,j  ,kp1)
        +(sal_t(i,j-1,kp1,1)+sal_t(i-1,j-1,kp1,1)+sal_t(i,j-1,kp1,0)+sal_t(i-1,j-1,kp1,0)) & ! +sal_u(i,j-1,kp1)
        -(sal_t(i,j  ,km1,1)+sal_t(i-1,j  ,km1,1)+sal_t(i,j  ,km1,0)+sal_t(i-1,j  ,km1,0)) & ! -sal_u(i,j  ,km1)
        -(sal_t(i,j-1,km1,1)+sal_t(i-1,j-1,km1,1)+sal_t(i,j-1,km1,0)+sal_t(i-1,j-1,km1,0)) & ! -sal_u(i,j-1,km1)
                                           ) & !ooo>                                 ! ]
        /(depth_u(i,j,kp1)+depth_u(i,j-1,kp1)-depth_u(i,j,km1)-depth_u(i,j-1,km1)) & !   /(   z_u(i,j,kp1)  +z_u(i,j-1,kp1)  -z_u(i,j,km1)  -z_u(i,j-1,km1))
           ) & !PMX>
           *invdy_f(i,j) &
            ) & !ABS> abs(-v.dS/dx+u.dS/dy))

       +abs( & !ABS> abs( u.dS/dx+v.dS/dy)
            (veldydz_u(i,j  ,k,id_veltot)*invdy_u(i,j)+veldydz_u(i,j-1,k,id_veltot)*invdy_u(i,j-1)) &
      *( & !PMX>     
         (sal_t(i  ,j,k,1)+sal_t(i  ,j-1,k,1)+sal_t(i  ,j,k,0)+sal_t(i  ,j-1,k,0)) & !  sal_v(i  ,j) 
        -(sal_t(i-1,j,k,1)+sal_t(i-1,j-1,k,1)+sal_t(i-1,j,k,0)+sal_t(i-1,j-1,k,0)) & ! -sal_v(i-1,j) 
        -(depth_v(i,j,k)-depth_v(i-1,j,k))*( & !ooo>                                 ! -(z_v(i,j,k)-z_v(i-1,j,k))*[
         (sal_t(i  ,j,kp1,1)+sal_t(i  ,j-1,kp1,1)+sal_t(i  ,j,kp1,0)+sal_t(i  ,j-1,kp1,0)) & !  sal_v(i  ,j,kp1)
        +(sal_t(i-1,j,kp1,1)+sal_t(i-1,j-1,kp1,1)+sal_t(i-1,j,kp1,0)+sal_t(i-1,j-1,kp1,0)) & ! +sal_v(i-1,j,kp1)
        -(sal_t(i  ,j,km1,1)+sal_t(i  ,j-1,km1,1)+sal_t(i  ,j,km1,0)+sal_t(i  ,j-1,km1,0)) & ! -sal_v(i  ,j,km1)
        -(sal_t(i-1,j,km1,1)+sal_t(i-1,j-1,km1,1)+sal_t(i-1,j,km1,0)+sal_t(i-1,j-1,km1,0)) & ! -sal_v(i-1,j,km1)
                                           ) & !ooo>                                 ! ]
        /(depth_v(i,j,kp1)+depth_v(i-1,j,kp1)-depth_v(i,j,km1)-depth_v(i-1,j,km1)) & ! /(   z_v(i,j,kp1)+  z_v(i-1,j,kp1)  -z_v(i,j,km1)  -z_v(i-1,j,km1))
          ) & !PMX>
           *invdx_f(i,j)  &

           +(veldxdz_v(i,j  ,k,id_veltot)*invdx_v(i,j)+veldxdz_v(i,j-1,k,id_veltot)*invdx_v(i,j-1)) &
      *( & !PMX>     
         (sal_t(i,j  ,k,1)+sal_t(i-1,j  ,k,1)+sal_t(i,j  ,k,0)+sal_t(i-1,j  ,k,0)) & !  sal_u(i,j) 
        -(sal_t(i,j-1,k,1)+sal_t(i-1,j-1,k,1)+sal_t(i,j-1,k,0)+sal_t(i-1,j-1,k,0)) & !  sal_u(i,j-1) 
        -(depth_u(i,j,k)-depth_u(i,j-1,k))*( & !ooo>                                 ! -(z_u(i,j)-z_u(i,j-1)*[
         (sal_t(i,j  ,kp1,1)+sal_t(i-1,j  ,kp1,1)+sal_t(i,j  ,kp1,0)+sal_t(i-1,j  ,kp1,0)) & !  sal_u(i,j  ,kp1)
        +(sal_t(i,j-1,kp1,1)+sal_t(i-1,j-1,kp1,1)+sal_t(i,j-1,kp1,0)+sal_t(i-1,j-1,kp1,0)) & ! +sal_u(i,j-1,kp1)
        -(sal_t(i,j  ,km1,1)+sal_t(i-1,j  ,km1,1)+sal_t(i,j  ,km1,0)+sal_t(i-1,j  ,km1,0)) & ! -sal_u(i,j  ,km1)
        -(sal_t(i,j-1,km1,1)+sal_t(i-1,j-1,km1,1)+sal_t(i,j-1,km1,0)+sal_t(i-1,j-1,km1,0)) & ! -sal_u(i,j-1,km1)
                                           ) & !ooo>                                 ! ]
        /(depth_u(i,j,kp1)+depth_u(i,j-1,kp1)-depth_u(i,j,km1)-depth_u(i,j-1,km1)) & !   /(   z_u(i,j,kp1)  +z_u(i,j-1,kp1)  -z_u(i,j,km1)  -z_u(i,j-1,km1))
           ) & !PMX>
           *invdy_f(i,j) &
            ) & !ABS> abs( u.dS/dx+v.dS/dy)

            ,small2)   ! PARENTHESE 1

       enddo
       enddo
      enddo

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

      call obc_mpi_anyv3d(0,id_ofactort,'r')
      call obc_mpi_anyv3d(0,id_ofactors,'r')
!.....................
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

      do  k=1,kmax    
      do  j=1,jmax
      do  i=1,imax

      x1=sign(un,veldydz_u(i+1,j,k,id_veltot))
      x2=sign(un,veldydz_u(i,  j,k,id_veltot))
      x3=sign(un,veldxdz_v(i,j+1,k,id_veltot))
      x4=sign(un,veldxdz_v(i,j  ,k,id_veltot))
      x5=sign(un,omega_w(i,j,k+1,1))
      x6=sign(un,omega_w(i,j,k  ,1))

! Coef en facteur de bio_z(i,j,k):
      const0=(dz_t(i,j,k,0)+                                &
       (  veldydz_u(i+1,j  ,k  ,id_veltot)*dti_fw*(-1.-x1)                   &
         +veldydz_u(i  ,j  ,k  ,id_veltot)*dti_fw*( 1.-x2)                   &
         +veldxdz_v(i  ,j+1,k  ,id_veltot)*dti_fw*(-1.-x3)                   &
         +veldxdz_v(i  ,j  ,k  ,id_veltot)*dti_fw*( 1.-x4)  )/dxdy_t(i,j)    &
           +omega_w(i  ,j  ,k+1,1)*dti_fw*(-1.-x5)                   &
           +omega_w(i  ,j  ,k  ,1)*dti_fw*( 1.-x6)                   &
       )/dz_t(i,j,k,2)

! Coef en facteur de bio_t(i+1,j,k):
      const1= &
      veldydz_u(i+1,j  ,k  ,id_veltot)*dti_fw*(-1.+x1)/dxdy_t(i,j)/dz_t(i,j,k,2)

! Coef en facteur de bio_t(i-1,j,k):
      const2=  &
      veldydz_u(i  ,j  ,k  ,id_veltot)*dti_fw*( 1.+x2)/dxdy_t(i,j)/dz_t(i,j,k,2)

! Coef en facteur de bio_t(i,j+1,k):
      const3=  &
      veldxdz_v(i  ,j+1,k  ,id_veltot)*dti_fw*(-1.+x3)/dxdy_t(i,j)/dz_t(i,j,k,2)

! Coef en facteur de bio_t(i,j-1,k):
      const4=  &
      veldxdz_v(i  ,j  ,k  ,id_veltot)*dti_fw*( 1.+x4)/dxdy_t(i,j)/dz_t(i,j,k,2)

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
      use module_principal ; use module_parallele ; use module_my_outputs
      implicit none


! Advection horizontale
       call advection_scal_dz        ! computes the "true" dz
       call advection_scal_dtisub    ! computes loopmaxts et dti_lpsub
       call advection_scal_ofactor   ! computes anyv3d(:,:,:,id_ofactor)

       call advection_scal_hybcoef   ! computes anyv3d(:,:,:,id_hybcoefu) anyv3d(:,:,:,id_hybcoefv)

       if(flag_negdif_hor==1) then !>>> 
         call advection_scal_negdif_hori   !-> temf_t & salf_t

         call advection_scal_tem_hornegdif !30-09-21  
         call advection_scal_sal_hornegdif !30-09-21  
 
!       if(loopmaxts>looplimit_hor) then !>>>
        if(sch_imp_ts_u_glb==1.or.sch_imp_ts_v_glb==1) then !ooo>
!        call advection_scal_tem_horimpli
!        call advection_scal_sal_horimpli
         call advection_scal_ts_horimpli
        endif                                               !ooo>

       else                       !>>>
         stop 'Le cas looplimit_hor reste A faire'
! note: si temf_t et salf_t ont ete utilises par l'advection verticale retablir zero: 07-10-22
         if(flag_negdif_ver==1) then ; temf_t=0. ; salf_t=0. ; endif
         call advection_scal_tem  
         call advection_scal_sal  
       endif                      !>>>

! Advection verticale
       call advection_scal_vbc ! Vertical Boundary Condition !05-03-21

       if(flag_negdif_ver==1) then ! Diffusion negative
          call advection_scal_negdif_vert !-> temf_t & salf_t
       else
! note: si temf_t et salf_t ont ete utilises par l'advection horizontale retablir zero: 07-10-22
        if(flag_negdif_hor==1) then ; temf_t=0. ; salf_t=0. ; endif
       endif 

       call advection_scal_vertical        !20-05-18

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
      do k=1,kmax
      do j=1,jmax ; do i=1,imax
!       do k=kmerged_t(i,j)+1,kmax
         x1=max(x1,max(max(max(abs(veldydz_u(i  ,j  ,k,id_veltot))  &
                              ,abs(veldydz_u(i+1,j  ,k,id_veltot))) &
                              ,abs(veldxdz_v(i  ,j  ,k,id_veltot))) &
                              ,abs(veldxdz_v(i  ,j+1,k,id_veltot))) & 
                    *dti_lp                                 &
                    *wetmask_t(i,j)                         &
                       /anyv3d(i,j,k,id_dxdydz))
!       enddo



      enddo       ; enddo
      enddo
! Niveaux <= kmerged
!     do j=1,jmax ; do i=1,imax
!       sum1=0. ; sum2=0. ; sum3=0. ; sum4=0. ; sum5=0. ; sum6=0.
!       do k=kundermin_t(i,j),kmerged_t(i,j) !07-12-17
!        sum1=sum1+veldydz_u(i  ,j  ,k,id_veltot)
!        sum2=sum2+veldydz_u(i+1,j  ,k,id_veltot)
!        sum3=sum3+veldxdz_v(i  ,j  ,k,id_veltot)
!        sum4=sum4+veldxdz_v(i  ,j+1,k,id_veltot)
!        sum5=sum5     +dz_t(i  ,j  ,k,0)
!        sum6=sum6     +dz_t(i  ,j  ,k,2)
!       enddo
!        x1=max(x1,max(max(max(abs(sum1)  &
!                             ,abs(sum2)) &
!                             ,abs(sum3)) &
!                             ,abs(sum4)) & 
!                   *dti_lp                                  &
!                   *wetmask_t(i,j)                          &
!                     /(dxdy_t(i,j)*min(sum5,sum6))          &          ! methode 2
!              ) ! fermeture max
!     enddo       ; enddo


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
      do k=1,kmax
      do j=1,jmax ; do i=1,imax
!       do k=kmerged_t(i,j)+1,kmax

         x2=       max(max(max(abs(veldydz_u(i  ,j  ,k,id_veltot))  &
                              ,abs(veldydz_u(i+1,j  ,k,id_veltot))) &
                              ,abs(veldxdz_v(i  ,j  ,k,id_veltot))) &
                              ,abs(veldxdz_v(i  ,j+1,k,id_veltot))) & 
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
               write(10+par%rank,*)'veldydz_u(i  ,j  ,k,id_veltot)',veldydz_u(i  ,j  ,k,id_veltot)
               write(10+par%rank,*)'    vel_u(i  ,j  ,k,1)',vel_u(i  ,j  ,k,1)
               write(10+par%rank,*)'veldydz_u(i+1,j  ,k,id_veltot)',veldydz_u(i+1,j  ,k,id_veltot)
               write(10+par%rank,*)'    vel_u(i+1,j  ,k,1)',vel_u(i+1,j  ,k,1)
               write(10+par%rank,*)'veldxdz_v(i  ,j  ,k,id_veltot)',veldxdz_v(i ,j  ,k,id_veltot)
               write(10+par%rank,*)'    vel_v(i  ,j  ,k,1)',vel_v(i  ,j  ,k,1)
               write(10+par%rank,*)'veldxdz_v(i  ,j+1,k,id_veltot)',veldxdz_v(i ,j+1,k,id_veltot)
               write(10+par%rank,*)'    vel_v(i  ,j+1,k,1)',vel_v(i  ,j+1,k,1)
             endif                 !(0v0)>

!       enddo
      enddo       ; enddo
      enddo

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

!     looplimit_hor=1
       if(loopmaxts>looplimit_hor) then !pmx>
         id_velexp=0 ; id_velimp=2
         call advection_scal_imp_exp
       else                             !pmx>
        id_velexp=id_veltot ! la vitesse explicite et la vitesse totale c'est la meme chose
       endif                            !pmx>

      end subroutine advection_scal_dtisub

!..............................................................................

      subroutine advection_scal_hybcoef !11-03-21
      use module_principal
      use module_parallele
      implicit none
      integer :: cn_power_=1

      cn_power_=1                      ! valeur standard pour le quickest
      if(iadvec_ts_hor==iadvec_ts_hor_quickest2)cn_power_=2 ! valeur qui affaiblit la composante UP2 du quickest quand cn petit

      do k=1,kmax

       do j=0,jmax+1
       do i=0,imax+1
        xy_t(i,j,1)=1./anyv3d(i,j,k,id_dxdydz)
       enddo
       enddo

       do j=1,jmax
       do i=1,imax+1
        anyv3d(i,j,k,id_hybcoefu)=max(1.-abs(            & !pmx>
              veldydz_u(i,j,k,id_veltot)*dti_lp*max(xy_t(i,j,1),xy_t(i-1,j,1)) &
                                            )**cn_power_ & !pmx> !11-03-21
                                  ,0.)          &

! Detection de presence de panache (si dS/di > 5 upwind) 
       *(1.-min(max( 0.4*abs(sal_t(i,j,k,0)-sal_t(i-1,j,k,0))-1. ,0.),1.)) & !08-03-22

       *0.125*(  upwindriver_t(i,j) +upwindriver_t(i-1,j)  ) & 
              *(upwindwetdry_t(i,j)+upwindwetdry_t(i-1,j)  ) & !28-02-22 
                  *(upwzone0_t(i,j,k)  +upwzone0_t(i-1,j,k)) !21-05-20

! Note: en situation standard 0.5*(upwzone0_t(i,j,k)+upwzone0_t(i-1,j,k))=1
!       et devient 0 si schema 100% upwind

       enddo
       enddo

       do j=1,jmax+1
       do i=1,imax
        anyv3d(i,j,k,id_hybcoefv)=max(1.-abs(            & !pmx>
              veldxdz_v(i,j,k,id_veltot)*dti_lp*max(xy_t(i,j,1),xy_t(i,j-1,1)) &
                                            )**cn_power_ & !pmx> !11-03-21
                                  ,0.)          &

! Detection de presence de panache (si dS/dj > 5 upwind) 
       *(1.-min(max( 0.4*abs(sal_t(i,j,k,0)-sal_t(i,j-1,k,0))-1. ,0.),1.)) & !08-03-22

!      *wetmask_v(i,j)*0.5*(upwindriver_t(i,j)+upwindriver_t(i,j-1)) & !27-02-22 
       *0.125*(  upwindriver_t(i,j) +upwindriver_t(i,j-1)  ) &  
              *(upwindwetdry_t(i,j)+upwindwetdry_t(i,j-1)  ) & !28-02-22 
                  *(upwzone0_t(i,j,k)  +upwzone0_t(i,j-1,k))  !21-05-20

       enddo
       enddo

      enddo
!     anyv3d(:,:,:,id_hybcoefu)=1. ! 100% up3
!     anyv3d(:,:,:,id_hybcoefv)=1. ! 100% up3
!     anyv3d(:,:,:,id_hybcoefu)=0. ! 100% up2
!     anyv3d(:,:,:,id_hybcoefv)=0. ! 100% up2

      if(flag_negdif_hor==0) then !-ancien-limiteur-du-cas-sans-negdif->
!---------------------------------------
! LIMITEUR 2D (DE BULLES D'EAUX DENSES)

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
!      call obc_scal_bottom_rhp('zb')  ! echange z1 et z2
       call obc_scal_bottom_rhp('z2')  !26-01-23

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

      endif                       !-ancien-limiteur-du-cas-sans-negdif->

      if(flag_negdif_hor==1) then !-nouveau-limiteur-du-cas-AVEC-negdif->
!............
! LIMITEUR 3D

      call obc_scal_rhp ! obc et mpi 'zb'


! LIMITEUR Flux Oi
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1

      if(veldydz_u(i,j,k,id_veltot)>0) then !m°v°m>
! PARTIE SI U POSITIF:

!     if(abs(rhp_t(i-1,j,k)-rhp_t(i-2,j,k))*mask_u(i-1,j,k)  &
!      <=abs(rhp_t(i  ,j,k)-rhp_t(i-1,j,k))*mask_u(i  ,j,k)  & !24-02-22
!       *(1.-anyv3d(i,j,k,id_hybcoefu)))                     &
!            anyv3d(i,j,k,id_hybcoefu)=0. ! schema upwind

! Similaire au test precedent mais en plus s'annule si gradient de signe inversE !08-03-22
      if(   (rhp_t(i-1,j,k)-rhp_t(i-2,j,k))*(rhp_t(i  ,j,k)-rhp_t(i-1,j,k))*mask_u(i-1,j,k)  &
       <=   (rhp_t(i  ,j,k)-rhp_t(i-1,j,k))*(rhp_t(i  ,j,k)-rhp_t(i-1,j,k))  & !24-02-22
        *(1.-anyv3d(i,j,k,id_hybcoefu)))                     &
             anyv3d(i,j,k,id_hybcoefu)=0. ! schema upwind

      else                          !m°v°m>
! PARTIE SI U NEGATIF:

!     if(abs(rhp_t(i+1,j,k)-rhp_t(i  ,j,k))*mask_u(i+1,j,k)  &
!      <=abs(rhp_t(i  ,j,k)-rhp_t(i-1,j,k))*mask_u(i  ,j,k)  & !24-02-22
!       *(1.-anyv3d(i,j,k,id_hybcoefu)))                     &
!            anyv3d(i,j,k,id_hybcoefu)=0. ! schema upwind

! Similaire au test precedent mais en plus s'annule si gradient de signe inversE !08-03-22
      if(   (rhp_t(i+1,j,k)-rhp_t(i  ,j,k))*(rhp_t(i  ,j,k)-rhp_t(i-1,j,k))*mask_u(i+1,j,k)  &
       <=   (rhp_t(i  ,j,k)-rhp_t(i-1,j,k))*(rhp_t(i  ,j,k)-rhp_t(i-1,j,k))                  & !24-02-22
        *(1.-anyv3d(i,j,k,id_hybcoefu)))                     &
             anyv3d(i,j,k,id_hybcoefu)=0. ! schema upwind

      endif                         !m°v°m>

      enddo ; enddo ; enddo 

! LIMITEUR Flux Oj
      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax

      if(veldxdz_v(i,j,k,id_veltot)>0) then !m°v°m>
! PARTIE SI V POSITIF:

!     if(abs(rhp_t(i,j-1,k)-rhp_t(i,j-2,k))*mask_v(i,j-1,k)  &
!      <=abs(rhp_t(i  ,j,k)-rhp_t(i,j-1,k))*mask_v(i,j  ,k)  & !24-02-22
!       *(1.-anyv3d(i,j,k,id_hybcoefv)))                     &
!            anyv3d(i,j,k,id_hybcoefv)=0. ! schema upwind

      if(   (rhp_t(i,j-1,k)-rhp_t(i,j-2,k))*(rhp_t(i  ,j,k)-rhp_t(i,j-1,k))*mask_v(i,j-1,k)  &
       <=   (rhp_t(i  ,j,k)-rhp_t(i,j-1,k))*(rhp_t(i  ,j,k)-rhp_t(i,j-1,k))  & !24-02-22
        *(1.-anyv3d(i,j,k,id_hybcoefv)))                     &
             anyv3d(i,j,k,id_hybcoefv)=0. ! schema upwind

      else                          !m°v°m>
! PARTIE SI V NEGATIF:

!     if(abs(rhp_t(i,j+1,k)-rhp_t(i  ,j,k))*mask_v(i,j+1,k)  &
!      <=abs(rhp_t(i  ,j,k)-rhp_t(i,j-1,k))*mask_v(i,j  ,k)  & !24-02-22
!       *(1.-anyv3d(i,j,k,id_hybcoefv)))                     &
!            anyv3d(i,j,k,id_hybcoefv)=0. ! schema upwind

      if(   (rhp_t(i,j+1,k)-rhp_t(i  ,j,k))*(rhp_t(i  ,j,k)-rhp_t(i,j-1,k))*mask_v(i,j+1,k)  &
       <=   (rhp_t(i  ,j,k)-rhp_t(i,j-1,k))*(rhp_t(i  ,j,k)-rhp_t(i,j-1,k))  & !24-02-22
        *(1.-anyv3d(i,j,k,id_hybcoefv)))                     &
             anyv3d(i,j,k,id_hybcoefv)=0. ! schema upwind

      endif                         !m°v°m>

      enddo ; enddo ; enddo 

      endif                       !-nouveau-limiteur-du-cas-AVEC-negdif->

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

      subroutine advection_scal_vertical !09-05-21
      use module_principal ; use module_parallele ; use module_my_outputs
      implicit none
!     double precision :: grdamp_=2. ! lower gradient amplification
      integer :: id_cn_=0           &! current number "limited"
                ,id_we_=1           &! Explicit limited vertical velocity
                ,id_wi_=2           &! residual implicit vertical velocity
!               ,id_tdiv=4          &! cumulated T*div(current) !26-01-17
!               ,id_sdiv=5          &! cumulated S*div(current)
                ,loop_,loopmax_     &
                ,looplimit_=1       & ! nombre d'iteration maximum la procedure iterative
                ,cnpower_=1         &
                ,cst1_=1            &
                ,cst2_=1  

!     if(flag_negdif_ver==0) then !m°v°m> !29-11-21 la coherence avec flag_negdif_ver est verifiEe dans set_parameters
!     if(iadvec_ts_ver==iadvec_ts_ver_c2) then !m°v°m> !29-11-21
!      cnpower_=2 ; cst1_=0 ; cst2_=0 !-- SCHEMA: (1-cn**2)*CE   +(cn**2)*UPWIND pour le cas iadvec_ts_ver_c2

      cnpower_=0 ! pour test debug                            !18-07-22
      if(iadvec_ts_ver==iadvec_ts_ver_quickest)  then !m°v°m> !18-07-22
       cnpower_=1 ; cst1_=1 ; cst2_=1 !-- SCHEMA: QUICKEST
      endif                                           !m°v°m>
      if(iadvec_ts_ver==iadvec_ts_ver_quickest2) then !m°v°m> !18-07-22
       cnpower_=2 ; cst1_=1 ; cst2_=0 !-- SCHEMA: QUICKEST-2: (1-cn**2)*QUICK+(cn**2)*UPWIND
      endif                                           !m°v°m>
      if(iadvec_ts_ver==iadvec_ts_ver_c2)  then       !m°v°m> !06-10-22
       cnpower_=1 ; cst1_=0 ; cst2_=0 !-- SCHEMA: c2lim
! Note: oui c'est normal d'avoir cnpower_=1 Details dans https://docs.google.com/document/d/1nccFYyPX3k_fsdX53wbLxvLv_nk0Pyt8l8kXTGXPPGs/edit?usp=sharing

      endif                                           !m°v°m> !06-10-22
      if(cnpower_==0)stop 'Err 2398 advection_scal_vertical' !test debug !18-07-22

!     x1=0. ! Pas obligatoire seulement pour diag du maximum
      xy_t(:,:,1)=small2 ! reset small2 evite que le nombre d'iteration soit nulle
      do k=2,kmax ; do j=1,jmax ; do i=1,imax

!   x0=Dz/Dt
        x0=min(dz_t(i,j,k,0),dz_t(i,j,k-1,0))/dti_lp
!   Omega bornee par + ou - N fois dz/dt ici 2 fois
        anyv3d(i,j,k,id_we_)=max(min(omega_w(i,j,k,1),looplimit_*x0),-looplimit_*x0) &
!                           *upwindriver_t(i,j) &
                            *wetmask_t(i,j)     &
                            *wetmask_wi_t(i,j)     !19-06-19
!   max over z of the decimal current number:
        xy_t(i,j,1)=max(xy_t(i,j,1),abs(anyv3d(i,j,k,id_we_))/x0)
!   x1=max(x1,xy_t(i,j,1)) ! Pas obligatoire seulement pour diag du maximum

      enddo  ; enddo ; enddo


! Fond w explicite
      do j=1,jmax ; do i=1,imax
!       anyv3d(i,j,1,id_we_)=0.
        anyv3d(i,j,1,id_we_)=omega_w(i,j,1,1) ! peut etre non nul si source sous marine !05-03-21
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
!      anyv3d(i,j,kmax+1,id_we_)=0. ! 25-02-19
! Schema explicite pour les flux et le reste (schema intertidale n°6) est implicite
       anyv3d(i,j,kmax+1,id_we_)=0.5*(omega_evaprec_w(i,j,0)+omega_evaprec_w(i,j,1)) !05-03-21
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


! Limiteur C0 !24-02-22
      anyv1d(:,1)=1. 
      rhp_t(i,j,0)     =rhp_t(i,j,1)    !03-03-22
      rhp_t(i,j,kmax+1)=rhp_t(i,j,kmax) !03-03-22

      do k=2,kmax

      if(anyv3d(i,j,k,id_we_)>0) then !m°v°m>
! PARTIE SI U POSITIF:

      if(abs(rhp_t(i,j,k-1)-rhp_t(i,j,k-2))  &
       <=abs(rhp_t(i,j,k  )-rhp_t(i,j,k-1))  & 
      *min(abs(anyv3d(i,j,k,id_we_)*dti_lpsub/min(dz_t(i,j,k,0),dz_t(i,j,k-1,0))),1.)) & !nc
             anyv1d(k,1)=0. ! schema upwind

      else                          !m°v°m>
! PARTIE SI U NEGATIF:

      if(abs(rhp_t(i,j,k+1)-rhp_t(i,j,k  ))  &
       <=abs(rhp_t(i,j,k  )-rhp_t(i,j,k-1))  & 
      *min(abs(anyv3d(i,j,k,id_we_)*dti_lpsub/min(dz_t(i,j,k,0),dz_t(i,j,k-1,0))),1.)) & !nc
             anyv1d(k,1)=0. ! schema upwind

      endif                         !m°v°m>

      enddo
!     anyv1d(:,1)=1. ! BIDOUILLE PATRICK ZERO LIMITEUR 
!     anyv1d(:,1)=0. ! BIDOUILLE PATRICK UPWIND

      do k=2,kmax  ! do k=kmin_w(i,j)+1,kmax                         !22-12-14


! Current number with additional limititations
       anyv3d(i,j,k,id_cn_)=( &         !power>
                              1.-( &!ooooo>
! Robustness regarding convection: 
! https://docs.google.com/document/d/1WsIhfKI1GP269wglgIM6YUKXdVrONILKupKP9Gui9hI/edit
!                      anyv1d(k+1,1)*         & ! limiteur C3 cette ligne et les 2 lignes suivantes
!                      anyv1d(k  ,1)*         & ! limiteur C2 cette ligne et la ligne suivante
!                      anyv1d(k-1,1)*         & ! limiteur C1 cette ligne seulement
                       anyv1d(k  ,1)*         & ! limiteur C0 cette ligne seulement
! Robustness regarding rivers and dry zones:
                       upwindriver_t(i,j)   & ! river zone influence
                         *upwzone0_t(i,j,k) & ! additional requests for upwind. Note upwzone0_t=1 in standard conditions
                     *upwindwetdry_t(i,j)   & ! wet/dry zone !28-02-22
! Robustness regarding vertical current number:
        *max(1.-abs(anyv3d(i,j,k,id_we_)*dti_lpsub               &
                 /min(dz_t(i,j,k,0),dz_t(i,j,k-1,0))),0.)        &  
                                 ) &!ooooo>                        
                            )**cnpower_  !09-05-21

!        anyv3d(i,j,k,id_cn_)=0. ! 100%C2
!        anyv3d(i,j,k,id_cn_)=1. ! 100%UP2

       enddo


! Cas particuliers k=1:kmin_w(i,j) et k=kmax:
!      anyv3d(i,j,1:kmin_w(i,j),id_cn_)=0.
!      anyv3d(i,j,1            ,id_cn_)=0.
!      anyv3d(i,j,kmax+1       ,id_cn_)=0. ! A la surface schema 100% explicite pour
!                                          ! garantir la coherence du flux d'eau douce
!                                          ! qui, dans vertmix_sal, est calcule avec 
!                                          ! avec sal_t(:,:,kmax,1). (Conservation bilan 
!                                          ! de sel) !12-04-15
       anyv3d(i,j,1            ,id_cn_)=1. !05-03-21
       anyv3d(i,j,kmax+1       ,id_cn_)=1. ! A la surface schema upwind 
! Noter que le choix nc=0 ou 1 n'aurait pas d'importance si la vitesse verticale est
! 100% implicite en kmax+1 . Ceci etant dit, noter que la vitesse
! verticale associee au flux atmospherique de surface est potentiellement explicite 
! (ca depend des versions, donc A verifier) et donc cela a de l'importance.
! Voir mes commentaires à propos du passage en implicite en kmax+1:
! https://docs.google.com/document/d/1fIzAG9mo_zvtTdVauu7yw-VkIqFRldqq4tF2b6ZOnzI/edit
! PAR CONTRE CELA A DE L'IMPORTANCE depuis la diffusion negative car celle-ci
! s'annule si nc=1, et par consequent elle n'interfere pas dans le flux de surface


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
!      anyv1d(1     ,1:4)=0. !05-03-21
!      anyv1d(kmax+1,3:4)=0. !05-03-21

! note: si anyv3d(i,j,kmax+1,id_cn_)=1 alors en fait anyv1d(kmax+1,1:2)=0.
       anyv1d(kmax+1,1:2)=0. !29-11-21
       anyv1d(1     ,1:2)=0.
!      k=kmax+1
!      anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)*tem_t(i,j,k-1,1)
!      anyv1d(k,2)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)*sal_t(i,j,k-1,1)

! Partie Iterative:
       do loop_=1,loopmax_ ! boucle iterative >>>


! Les lignes suivantes sont commentees car les conditions aux limites
! sont dans la subroutine advection_scal_vbc  :
! A priori lignes suivantes inutiles si schema implicite en kmax+1
!      tem_t(i,j,kmax+1,2)=tem_t(i,j,kmax,2)
!      tem_t(i,j,0     ,2)=tem_t(i,j,1   ,2)
!      sal_t(i,j,kmax+1,2)=sal_t(i,j,kmax,2)
!      sal_t(i,j,0     ,2)=sal_t(i,j,1   ,2)
! Si schema implicite pour la partie de omega associee au schema intertidate n° 6
! et schema explicite pour la partie flux alors les hypotheses sur T et S en kmax+1
! sont sal(kmax+1)=0 et tem(kmax+1)=tem(kmax) ou temobc ou triver ou autre A definir
!      tem_t(i,j,kmax+1,2)=tem_t(i,j,kmax,2)
!      tem_t(i,j,kmax+1,2)=0.5*(tem_t(i,j,kmax,0)+tem_t(i,j,kmax,1))
!      tem_t(i,j,kmax+1,2)=timeweightobc(trc_id)*temobc_t(i,j,kmax,2)+(1.-timeweightobc(trc_id))*temobc_t(i,j,kmax,0)
!      sal_t(i,j,kmax+1,2)=0.


! La boucle k de 1 1 kmax+1 pour inclure la condition limite de surface
! et fond en T et S de la subroutine advection_scal_vbc
       do k=2,kmax   !nexon> !05-03-21

! Flux upwind:
        anyv1d(k,3)=  &
          -anyv3d(i,j,k,id_cn_)*0.5*( & !ooo>
                anyv3d(i,j,k,id_we_) *(tem_t(i,j,k,2)+tem_t(i,j,k-1,2)) &
           -abs(anyv3d(i,j,k,id_we_))*( & !m°v°m>
                                       tem_t(i,j,k,2)-tem_t(i,j,k-1,2)  & !Dif
          -(1.-anyv3d(i,j,k,id_cn_))*(temf_t(i,j,k) -temf_t(i,j,k-1)  ) & !Neg Dif. !04-05-21
                                      ) & !m°v°m>
                                    ) & !ooo>

       +(1-anyv3d(i,j,k,id_cn_)) &
       *(1+anyv3d(i,j,k,id_cn_)*cst2_)*cst1_  & !09-05-21
                                *(  & !quick>
          0.5*(anyv3d(i,j,k,id_we_)+abs(anyv3d(i,j,k,id_we_)))*( & !ooo>
               tem_t(i,j,k  ,2)-2.* tem_t(i,j,k-1,2)+ tem_t(i,j,k-2,2) &
             -temf_t(i,j,k)    +2.*temf_t(i,j,k-1)  -temf_t(i,j,k-2)   &
                                                               ) & !ooo>

         +0.5*(anyv3d(i,j,k,id_we_)-abs(anyv3d(i,j,k,id_we_)))*( & !xxx>
               tem_t(i,j,k+1,2)-2.* tem_t(i,j,k  ,2)+ tem_t(i,j,k-1,2) &
             -temf_t(i,j,k+1)  +2.*temf_t(i,j,k  )  -temf_t(i,j,k-1)   &
                                                               ) & !xxx>
!                                )/6. !quick>
                                 )*quick_coef !quick> !18-07-22

        anyv1d(k,4)=  &
          -anyv3d(i,j,k,id_cn_)*0.5*( & !ooo>
                anyv3d(i,j,k,id_we_) *(sal_t(i,j,k,2)+sal_t(i,j,k-1,2)) &
           -abs(anyv3d(i,j,k,id_we_))*( & !m°v°m>
                                       sal_t(i,j,k,2)-sal_t(i,j,k-1,2)  & !Dif
          -(1.-anyv3d(i,j,k,id_cn_))*(salf_t(i,j,k) -salf_t(i,j,k-1)  ) & !Neg Dif. !04-05-21
                                      ) & !m°v°m>
                                    ) & !ooo>

       +(1-anyv3d(i,j,k,id_cn_)) &
       *(1+anyv3d(i,j,k,id_cn_)*cst2_)*cst1_  & !09-05-21
                                *(  & !quick>
          0.5*(anyv3d(i,j,k,id_we_)+abs(anyv3d(i,j,k,id_we_)))*( & !ooo>
               sal_t(i,j,k  ,2)-2.* sal_t(i,j,k-1,2)+ sal_t(i,j,k-2,2) &
             -salf_t(i,j,k)    +2.*salf_t(i,j,k-1)  -salf_t(i,j,k-2)   &
                                                               ) & !ooo>

         +0.5*(anyv3d(i,j,k,id_we_)-abs(anyv3d(i,j,k,id_we_)))*( & !xxx>
               sal_t(i,j,k+1,2)-2.* sal_t(i,j,k  ,2)+ sal_t(i,j,k-1,2) &
             -salf_t(i,j,k+1)  +2.*salf_t(i,j,k  )  -salf_t(i,j,k-1)   &
                                                               ) & !xxx>
!                                )/6. !quick>
                                 )*quick_coef !quick> !18-07-22

       enddo       !nexon>
! Attention que quelque soit le signe du courant, le flux explicite de sel est nul au fond et a la surface
! (sachant que l'implicite a lui la condition d/dz=0 de la zone intertidale schema 6)
       anyv1d(kmax+1,4)=0.
       anyv1d(1     ,4)=0.
       do k=1,kmax+1,kmax
        anyv1d(k,3)=  &
          -anyv3d(i,j,k,id_cn_)*0.5*( & !ooo>
                anyv3d(i,j,k,id_we_) *(tem_t(i,j,k,2)+tem_t(i,j,k-1,2)) &
           -abs(anyv3d(i,j,k,id_we_))*( & !m°v°m>
                                       tem_t(i,j,k,2)-tem_t(i,j,k-1,2)  & !Dif
!         -(1.-anyv3d(i,j,k,id_cn_))*(temf_t(i,j,k) -temf_t(i,j,k-1)  ) & !Neg Dif. !04-05-21
                                      ) & !m°v°m>
                                    )   !ooo>
       enddo


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

      end subroutine advection_scal_vertical

!................................................................

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

!     x1=0. ! Pas obligatoire seulement pour diag du maximum
      xy_t(:,:,1)=small2 ! reset small2 evite que le nombre d'iteration soit nulle
      do k=2,kmax ; do j=1,jmax ; do i=1,imax

!   x0=Dz/Dt
        x0=min(dz_t(i,j,k,0),dz_t(i,j,k-1,0))/dti_lp
!   Omega bornee par + ou - N fois dz/dt ici 2 fois
        anyv3d(i,j,k,id_we_)=max(min(omega_w(i,j,k,1),looplimit_*x0),-looplimit_*x0) &
!                           *upwindriver_t(i,j) &
                            *wetmask_t(i,j)     &
                            *wetmask_wi_t(i,j)     !19-06-19
!   max over z of the decimal current number:
        xy_t(i,j,1)=max(xy_t(i,j,1),abs(anyv3d(i,j,k,id_we_))/x0)
!   x1=max(x1,xy_t(i,j,1)) ! Pas obligatoire seulement pour diag du maximum

      enddo  ; enddo ; enddo

! Fond w explicite
      do j=1,jmax ; do i=1,imax
!       anyv3d(i,j,1,id_we_)=0.
        anyv3d(i,j,1,id_we_)=omega_w(i,j,1,1) ! peut etre non nul si source sous marine
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
!      anyv3d(i,j,kmax+1,id_we_)=0. ! 25-02-19
! Schema explicite pour les flux et le reste (schema intertidale n°6) est implicite
       anyv3d(i,j,kmax+1,id_we_)=0.5*(omega_evaprec_w(i,j,0)+omega_evaprec_w(i,j,1))
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
                            ) !**cn_power !power !26-05-18

!        anyv3d(i,j,k,id_cn_)=0. ! 100%C2
!        anyv3d(i,j,k,id_cn_)=1. ! 100%UP2

       enddo


! Cas particuliers k=1:kmin_w(i,j) et k=kmax:
!      anyv3d(i,j,1:kmin_w(i,j),id_cn_)=0.
!      anyv3d(i,j,1            ,id_cn_)=0.
!      anyv3d(i,j,kmax+1       ,id_cn_)=0. ! A la surface schema 100% explicite pour
!                                          ! garantir la coherence du flux d'eau douce
!                                          ! qui, dans vertmix_sal, est calcule avec 
!                                          ! avec sal_t(:,:,kmax,1). (Conservation bilan 
!                                          ! de sel) !12-04-15
       anyv3d(i,j,1            ,id_cn_)=1.
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
!      anyv1d(1     ,1:4)=0.
!      anyv1d(kmax+1,3:4)=0.

       anyv1d(1     ,1:2)=0.
       k=kmax+1
       anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)*tem_t(i,j,k-1,1)
       anyv1d(k,2)=-(1-anyv3d(i,j,k,id_cn_))*anyv3d(i,j,k,id_we_)*sal_t(i,j,k-1,1)

! Partie Iterative:
       do loop_=1,loopmax_ ! boucle iterative >>>


! Les lignes suivantes sont commentees car les conditions aux limites
! sont dans la subroutine advection_scal_vbc  :
! A priori lignes suivantes inutiles si schema implicite en kmax+1
!      tem_t(i,j,kmax+1,2)=tem_t(i,j,kmax,2)
!      tem_t(i,j,0     ,2)=tem_t(i,j,1   ,2)
!      sal_t(i,j,kmax+1,2)=sal_t(i,j,kmax,2)
!      sal_t(i,j,0     ,2)=sal_t(i,j,1   ,2)
! Si schema implicite pour la partie de omega associee au schema intertidate n° 6
! et schema explicite pour la partie flux alors les hypotheses sur T et S en kmax+1
! sont sal(kmax+1)=0 et tem(kmax+1)=tem(kmax) ou temobc ou triver ou autre A definir
!      tem_t(i,j,kmax+1,2)=tem_t(i,j,kmax,2)
!      tem_t(i,j,kmax+1,2)=0.5*(tem_t(i,j,kmax,0)+tem_t(i,j,kmax,1))
!      tem_t(i,j,kmax+1,2)=timeweightobc(trc_id)*temobc_t(i,j,kmax,2)+(1.-timeweightobc(trc_id))*temobc_t(i,j,kmax,0)
!      sal_t(i,j,kmax+1,2)=0.


! La boucle k de 1 1 kmax+1 pour inclure la condition limite de surface
! et fond en T et S de la subroutine advection_scal_vbc
       do k=2,kmax   !nexon>

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


        anyv1d(k,4)=  &
          -anyv3d(i,j,k,id_cn_) *(x1*sal_t(i,j,k-1,2)+x2*sal_t(i,j,k  ,2))&
       +(1-anyv3d(i,j,k,id_cn_))*( &
               x1*(sal_t(i,j,k  ,2)-2.*sal_t(i,j,k-1,2)+sal_t(i,j,k-2,2))&!21-11-16
              +x2*(sal_t(i,j,k+1,2)-2.*sal_t(i,j,k  ,2)+sal_t(i,j,k-1,2))&
                                 )*(1+anyv3d(i,j,k,id_cn_))/6. ! QUICKEST
!                                                          )/6. ! UBS-UPW


       enddo       !nexon>

! En k=1 et en k=kmax+1
       do k=1,kmax+1,kmax ! en k=1 et en k=kmax+1
        x1=0.5*(anyv3d(i,j,k,id_we_)+abs(anyv3d(i,j,k,id_we_)))
        x2=0.5*(anyv3d(i,j,k,id_we_)-abs(anyv3d(i,j,k,id_we_)))
        anyv1d(k,3)=  &
          -anyv3d(i,j,k,id_cn_) *(x1*tem_t(i,j,k-1,2)+x2*tem_t(i,j,k  ,2)) 
! Attention que quelque soit le signe du courant, le flux explicite de sel est nul au fond et a la surface
! (sachant que l'implicite a lui la condition d/dz=0 de la zone intertidale schema 6)
        anyv1d(k,4)=0.
       enddo


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

      end subroutine advection_scal_z_quickest

!................................................................

      subroutine advection_scal_tem  !11-03-21
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      double precision :: velef_,flagmsk_=1.,scalefac_,cst_,cst2_=2.,cst1_=1.
      integer :: t_=1        &
                ,id_flx_=1   &  ! x flux identifier
                ,id_fly_=1   &  ! y flux identifier
                ,id_tbef_=3  &  ! temperature    "before" identifier
!               ,id_tdiv=4   &  ! T*divergence(velocity) !26-01-17
                ,loop_

! https://docs.google.com/document/d/1QLszDYbzYsIavgI1MgexRbpfEtN6rZDSzRiqdkCsIKE/edit?usp=sharing

!     cst_=1./6. ! ou 1./8.
      cst_=quick_coef !18-07-22
      if(iadvec_ts_hor==iadvec_ts_hor_quickest2) then !29-11-21
! Quickest-2
!        cst_=1./8. ! ou 1./8.
         cst2_=1.
         cst1_=0.
      else
! Quickest
         cst2_=2.
         cst1_=1.
      endif

! A quoi sert cst2_ dans l'expression (cst2-anyv3d(i,j,k,id_hybcoefu)*cst1) ?
! Elle sert a choisir entre QUICK et QUICKEST
! vaut 2-anyv3d(i,j,k,id_hybcoefu)=1+nc pour QUICKEST et vaut 1 pour QUICK



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
      -veldydz_u(i,j,k,id_veltot) &
        *(    tem_t(i,j,k,now)+tem_t(i-1,j,k,now)      &
                                                 )     &
! Partie diffusive QUICK remise A jour aux sous-pas de temps advectif:
       +cst_*max(anyv3d(i,j,k,id_ofactort),anyv3d(i,j+1,k,id_ofactort)) &
       *(cst2_-anyv3d(i,j,k,id_hybcoefu)*cst1_) & !11-03-21 ! https://docs.google.com/document/d/1QLszDYbzYsIavgI1MgexRbpfEtN6rZDSzRiqdkCsIKE/edit?usp=sharing

                                                               *( & !m°v°m>
             (veldydz_u(i,j,k,id_veltot)+abs(veldydz_u(i,j,k,id_veltot)))*( &           !pmxpmx>

               (   anyv3d(i  ,j,k,id_tbef_)   -anyv3d(i-1,j,k,id_tbef_)  &
!                 -temf_t(i  ,j,k         )   +temf_t(i-1,j,k         )  &

! correction anyv3d(i  ,j,k,id_tbef_) par rapport niveau z(i-1,j,k)
      +( anyv3d(i  ,j,kp1,id_tbef_)- anyv3d(i  ,j,km1,id_tbef_)     &
!       -temf_t(i  ,j,kp1         )+ temf_t(i  ,j,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i-1,j,k  )         -depth_t(i  ,j,k  )         )    &
      /(depth_t(i  ,j,kp1)         -depth_t(i  ,j,km1)         )    &
          )       )      & !minmax>

               ) &

              -(   anyv3d(i-1,j,k,id_tbef_)-anyv3d(i-2,j,k,id_tbef_)  &
!                 -temf_t(i-1,j,k         )+temf_t(i-2,j,k         )  &

! correction anyv3d(i-2,j,k,id_tbef_) par rapport niveau z(i-1,j,k)
      -( anyv3d(i-2,j,kp1,id_tbef_)- anyv3d(i-2,j,km1,id_tbef_)     &
!       -temf_t(i-2,j,kp1         )+ temf_t(i-2,j,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i-1,j,k  )         -depth_t(i-2,j,k  )         )    &
      /(depth_t(i-2,j,kp1)         -depth_t(i-2,j,km1)         )    &
          )       )      & !minmax>

               ) &
                                           *mask_t(i-2,j,kmax)      ) & !pmxpmx>

            +(veldydz_u(i,j,k,id_veltot)-abs(veldydz_u(i,j,k,id_veltot)))*( &           !pmxpmx>

               (   anyv3d(i+1,j,k,id_tbef_)-anyv3d(i,j,k,id_tbef_)    &
!                 -temf_t(i+1,j,k         )+temf_t(i,j,k         )    &

! correction anyv3d(i+1,j,k,id_tbef_) par rapport niveau z(i  ,j,k)
      +( anyv3d(i+1,j,kp1,id_tbef_)- anyv3d(i+1,j,km1,id_tbef_)     &
!       -temf_t(i+1,j,kp1         )+ temf_t(i+1,j,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i  ,j,k  )         -depth_t(i+1,j,k  )         )    &
      /(depth_t(i+1,j,kp1)         -depth_t(i+1,j,km1)         )    &
          )       )      & !minmax>

               )   &
                  *mask_t(i+1,j,kmax)                                 &

              -(   anyv3d(i,j,k,id_tbef_)-anyv3d(i-1,j,k,id_tbef_)    &
!                 -temf_t(i,j,k         )+temf_t(i-1,j,k         )    &

! correction anyv3d(i-1,j,k,id_tbef_) par rapport niveau z(i  ,j,k)
      -( anyv3d(i-1,j,kp1,id_tbef_)- anyv3d(i-1,j,km1,id_tbef_)     &
!       -temf_t(i-1,j,kp1         )+ temf_t(i-1,j,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i  ,j,k  )         -depth_t(i-1,j,k  )         )    &
      /(depth_t(i-1,j,kp1)         -depth_t(i-1,j,km1)         )    &
          )       )      & !minmax>

                                                                  ))  & !pmxpmx>

                                                                )  & !m°v°m>()cst_
                                                                )  & !()anyv3d(i,j,k,id_hybcoefu)
! Schema UPWIND dans les zones sous influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      +(1.-anyv3d(i,j,k,id_hybcoefu))*( & !pmx>
           -veldydz_u(i,j,k,id_veltot)* ( & !ppp>
               anyv3d(i,j,k,id_tbef_)+anyv3d(i-1,j,k,id_tbef_)   &
                                ) & !ppp>
       +abs(veldydz_u(i,j,k,id_veltot))*( & !mmm>
               anyv3d(i,j,k,id_tbef_)-anyv3d(i-1,j,k,id_tbef_)   &
!            -(temf_t(i,j,k)         -temf_t(i-1,j,k))*anyv3d(i,j,k,id_hybcoefu) & !05-05-21
                                ) & !mmm>
                                      ) & !pmx>

       ) !prime>

!     -veldydz_u(i,j,k,id_veltot)*dti_lpsub*temref_u(i,j,k)


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
      +anyv3d(i,j,k,id_tbef_)*dti_lpsub*(veldydz_u(i+1,j,k,id_veltot)           &
                                        -veldydz_u(i  ,j,k,id_veltot))/dxdy_t(i,j)

!       anyv3d(i,j,k,id_tbef_)                       &
        anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !ooo> !25-04-17

        anyv3d(i,j,k,id_tbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i+1,j,k,id_flx_)                &
             -anyv3d(i  ,j,k,id_flx_)                &

         +(veldydz_u(i+1,j,k,id_veltot)                      &
          -veldydz_u(i  ,j,k,id_veltot))*anyv3d(i,j,k,id_tbef_) & !26-01-17


                               )/dxdy_t(i,j) ) & !ooo>
                                                /dz_t(i,j,k,before) &!02-11-21 !retablir division par dz si appel A vertmix_merged_levels_any supprimE 

         +(1-mask_t(i,j,k))*tem_t(i,j,k,before)!*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

! Note: le terme anyv3d(i,j,k,id_tdiv) ne fait pas l'objet d'un traitement similaire pour le masque car
! si mask_t=0 alors veldydz_u(i:i+1,...)=0 (question: et les fleuves ?)
      enddo
      enddo
      enddo


      !call vertmix_merged_levels_any(0,id_tbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_tbef_)
      call obc_mpi_anyv3d(1,id_tbef_,'zb')

! Flux Oj
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax+1
      do i=1,imax

      anyv3d(i,j,k,id_fly_)=0.5*( & !prime>

! Schema QUICK en dehors influence des fleuves
      anyv3d(i,j,k,id_hybcoefv)*( &
! Partie advective QUICK:
      -veldxdz_v(i,j,k,id_veltot)  &
        *(    tem_t(i,j,k,now)+tem_t(i,j-1,k,now)  &
                                                 ) &
! Partie diffusive QUICK:
       +cst_*max(anyv3d(i,j,k,id_ofactort),anyv3d(i+1,j,k,id_ofactort)) &
       *(cst2_-anyv3d(i,j,k,id_hybcoefv)*cst1_) & !11-03-21 ! https://docs.google.com/document/d/1QLszDYbzYsIavgI1MgexRbpfEtN6rZDSzRiqdkCsIKE/edit?usp=sharing

                                                                  *(  & !m°v°m>

             (veldxdz_v(i,j,k,id_veltot)+abs(veldxdz_v(i,j,k,id_veltot)))*(  &       !ooo>

            (   anyv3d(i,j  ,k,id_tbef_)-anyv3d(i,j-1,k,id_tbef_)  &
!              -temf_t(i,j  ,k         )+temf_t(i,j-1,k         )  &

! correction anyv3d(i,j  ,k,id_tbef_) par rapport niveau z(i,j-1,k)
      +( anyv3d(i,j  ,kp1,id_tbef_)- anyv3d(i,j  ,km1,id_tbef_)     &
!       -temf_t(i,j  ,kp1         )+ temf_t(i,j  ,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j-1,k  )         -depth_t(i,j  ,k  )         )    &
      /(depth_t(i,j  ,kp1)         -depth_t(i,j  ,km1)         )    &
          )       )      & !minmax>

            ) &

           -(   anyv3d(i,j-1,k,id_tbef_)-anyv3d(i,j-2,k,id_tbef_)  &
!              -temf_t(i,j-1,k         )+temf_t(i,j-2,k         )  &

! correction anyv3d(i,j-2,k,id_tbef_) par rapport niveau z(i,j-1,k)
      -( anyv3d(i,j-2,kp1,id_tbef_)- anyv3d(i,j-2,km1,id_tbef_)     &
!       -temf_t(i,j-2,kp1         )+ temf_t(i,j-2,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j-1,k  )         -depth_t(i,j-2,k  )         )    &
      /(depth_t(i,j-2,kp1)         -depth_t(i,j-2,km1)         )    &
          )       )      & !minmax>

            ) &
                                        *mask_t(i,j-2,kmax)       )& !ooo>

            +(veldxdz_v(i,j,k,id_veltot)-abs(veldxdz_v(i,j,k,id_veltot)))*(  &       !ooo>

            (   anyv3d(i,j+1,k,id_tbef_)-anyv3d(i,j,k,id_tbef_)    &
!              -temf_t(i,j+1,k         )+temf_t(i,j,k         )    &

! correction anyv3d(i,j+1,k,id_tbef_) par rapport niveau z(i,j  ,k)
      +( anyv3d(i,j+1,kp1,id_tbef_)- anyv3d(i,j+1,km1,id_tbef_)     &
!       -temf_t(i,j+1,kp1         )+ temf_t(i,j+1,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j  ,k  )         -depth_t(i,j+1,k  )         )    &
      /(depth_t(i,j+1,kp1)         -depth_t(i,j+1,km1)         )    &
          )       )      & !minmax>

            )   &
               *mask_t(i,j+1,kmax)                                 &

           -(   anyv3d(i,j,k,id_tbef_)-anyv3d(i,j-1,k,id_tbef_)    &
!              -temf_t(i,j,k         )+temf_t(i,j-1,k         )    &

! correction anyv3d(i,j-1,k,id_tbef_) par rapport niveau z(i,j  ,k)
      -( anyv3d(i,j-1,kp1,id_tbef_)- anyv3d(i,j-1,km1,id_tbef_)     &
!       -temf_t(i,j-1,kp1         )+ temf_t(i,j-1,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j  ,k  )         -depth_t(i,j-1,k  )         )    &
      /(depth_t(i,j-1,kp1)         -depth_t(i,j-1,km1)         )    &
          )       )      & !minmax>

                                                               ))  & !ooo>

                                                                   )  & !m°v°m>()cst_
                                                                   )  & !()anyv3d(i,j,k,id_hybcoefv)
! Schema UPWIND dans les zones sous influence des fleuves:
      +(1.-anyv3d(i,j,k,id_hybcoefv))*( & !pmx>
          -veldxdz_v(i,j,k,id_veltot)* ( & !ppp>
              anyv3d(i,j,k,id_tbef_)+anyv3d(i,j-1,k,id_tbef_)   &
                               ) & !ppp>
      +abs(veldxdz_v(i,j,k,id_veltot))*( & !mmm>
              anyv3d(i,j,k,id_tbef_)-anyv3d(i,j-1,k,id_tbef_)   &
!           -(temf_t(i,j,k)         -temf_t(i,j-1,k))*anyv3d(i,j,k,id_hybcoefv) & !05-05-21
                               ) & !mmm>
                                      ) & !pmx>

       ) !prime>

!      -veldxdz_v(i,j,k,id_veltot)*dti_lpsub*temref_v(i,j,k)

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
      +anyv3d(i,j,k,id_tbef_)*dti_lpsub*(veldxdz_v(i,j+1,k,id_veltot)                 &
                                        -veldxdz_v(i,j  ,k,id_veltot))/dxdy_t(i,j)

        anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !ooo>
        anyv3d(i,j,k,id_tbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i,j+1,k,id_fly_)                &
             -anyv3d(i,j  ,k,id_fly_)                &

         +(veldxdz_v(i,j+1,k,id_veltot)                      &
          -veldxdz_v(i,j  ,k,id_veltot))*anyv3d(i,j,k,id_tbef_) & !26-01-17

                              )/dxdy_t(i,j)  ) & !ooo>
                                                /dz_t(i,j,k,before) &!02-11-21 !retablir division par dz si appel A vertmix_merged_levels_any supprimE 
 
         +(1-mask_t(i,j,k))*tem_t(i,j,k,before)!*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

      enddo
      enddo
      enddo


      !call vertmix_merged_levels_any(0,id_tbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_tbef_)
      call obc_mpi_anyv3d(1,id_tbef_,'zb')

!#endif
      enddo               !*********>

      do k=1,kmax ; do j=1,jmax  ; do i=1,imax  

           tem_t(i,j,k,2)                              &
       =  anyv3d(i,j,k,id_tbef_)                    !  &

      enddo ; enddo ; enddo

      end subroutine advection_scal_tem 


!..............................................................................

      subroutine advection_scal_sal  !11-03-21
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      double precision :: velef_,flagmsk_=1.,scalefac_,cst_,cst2_=2.,cst1_=1.
      integer :: t_=1        &
                ,id_flx_=1   &  ! x flux identifier
                ,id_fly_=1   &  ! y flux identifier
                ,id_sbef_=3  &  ! temperature    "before" identifier
!               ,id_sdiv=4   &  ! T*divergence(velocity) !26-01-17
                ,loop_

! https://docs.google.com/document/d/1QLszDYbzYsIavgI1MgexRbpfEtN6rZDSzRiqdkCsIKE/edit?usp=sharing

!     cst_=1./6. ! ou 1./8.
      cst_=quick_coef !18-07-22
      if(iadvec_ts_hor==iadvec_ts_hor_quickest2) then
! Quickest-2
!        cst_=1./8. ! ou 1./8.
         cst2_=1.
         cst1_=0.
      else
! Quickest
         cst2_=2.
         cst1_=1.
      endif

! A quoi sert cst2_ dans l'expression (cst2-anyv3d(i,j,k,id_hybcoefu)*cst1) ?
! Elle sert a choisir entre QUICK et QUICKEST
! vaut 2-anyv3d(i,j,k,id_hybcoefu)=1+nc pour QUICKEST et vaut 1 pour QUICK



      do k=1,kmax ; do j=-1,jmax+2 ; do i=-1,imax+2
       anyv3d(i,j,k,id_sbef_)=sal_t(i,j,k,before)
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,id_sdiv)=0. ! reset cumul temperature*divergence(courant) !26-01-17
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
      -veldydz_u(i,j,k,id_veltot) &
        *(    sal_t(i,j,k,now)+sal_t(i-1,j,k,now)      &
                                                 )     &
! Partie diffusive QUICK remise A jour aux sous-pas de temps advectif:
       +cst_*max(anyv3d(i,j,k,id_ofactort),anyv3d(i,j+1,k,id_ofactort)) &
       *(cst2_-anyv3d(i,j,k,id_hybcoefu)*cst1_) & !11-03-21 ! https://docs.google.com/document/d/1QLszDYbzYsIavgI1MgexRbpfEtN6rZDSzRiqdkCsIKE/edit?usp=sharing

                                                               *( & !m°v°m>
             (veldydz_u(i,j,k,id_veltot)+abs(veldydz_u(i,j,k,id_veltot)))*( &           !pmxpmx>

               (   anyv3d(i  ,j,k,id_sbef_)   -anyv3d(i-1,j,k,id_sbef_)  &
!                 -salf_t(i  ,j,k         )   +salf_t(i-1,j,k         )  &

! correction anyv3d(i  ,j,k,id_sbef_) par rapport niveau z(i-1,j,k)
      +( anyv3d(i  ,j,kp1,id_sbef_)- anyv3d(i  ,j,km1,id_sbef_)     &
!       -salf_t(i  ,j,kp1         )+ salf_t(i  ,j,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i-1,j,k  )         -depth_t(i  ,j,k  )         )    &
      /(depth_t(i  ,j,kp1)         -depth_t(i  ,j,km1)         )    &
          )       )      & !minmax>

               ) &

              -(   anyv3d(i-1,j,k,id_sbef_)-anyv3d(i-2,j,k,id_sbef_)  &
!                 -salf_t(i-1,j,k         )+salf_t(i-2,j,k         )  &

! correction anyv3d(i-2,j,k,id_sbef_) par rapport niveau z(i-1,j,k)
      -( anyv3d(i-2,j,kp1,id_sbef_)- anyv3d(i-2,j,km1,id_sbef_)     &
!       -salf_t(i-2,j,kp1         )+ salf_t(i-2,j,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i-1,j,k  )         -depth_t(i-2,j,k  )         )    &
      /(depth_t(i-2,j,kp1)         -depth_t(i-2,j,km1)         )    &
          )       )      & !minmax>

               ) &
                                           *mask_t(i-2,j,kmax)      ) & !pmxpmx>

            +(veldydz_u(i,j,k,id_veltot)-abs(veldydz_u(i,j,k,id_veltot)))*( &           !pmxpmx>

               (   anyv3d(i+1,j,k,id_sbef_)-anyv3d(i,j,k,id_sbef_)    &
!                 -salf_t(i+1,j,k         )+salf_t(i,j,k         )    &

! correction anyv3d(i+1,j,k,id_sbef_) par rapport niveau z(i  ,j,k)
      +( anyv3d(i+1,j,kp1,id_sbef_)- anyv3d(i+1,j,km1,id_sbef_)     &
!       -salf_t(i+1,j,kp1         )+ salf_t(i+1,j,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i  ,j,k  )         -depth_t(i+1,j,k  )         )    &
      /(depth_t(i+1,j,kp1)         -depth_t(i+1,j,km1)         )    &
          )       )      & !minmax>

               )   &
                  *mask_t(i+1,j,kmax)                                 &

              -(   anyv3d(i,j,k,id_sbef_)-anyv3d(i-1,j,k,id_sbef_)    &
!                 -salf_t(i,j,k         )+salf_t(i-1,j,k         )    &

! correction anyv3d(i-1,j,k,id_sbef_) par rapport niveau z(i  ,j,k)
      -( anyv3d(i-1,j,kp1,id_sbef_)- anyv3d(i-1,j,km1,id_sbef_)     &
!       -salf_t(i-1,j,kp1         )+ salf_t(i-1,j,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i  ,j,k  )         -depth_t(i-1,j,k  )         )    &
      /(depth_t(i-1,j,kp1)         -depth_t(i-1,j,km1)         )    &
          )       )      & !minmax>

                                                                  ))  & !pmxpmx>

                                                                )  & !m°v°m>()cst_
                                                                )  & !()anyv3d(i,j,k,id_hybcoefu)
! Schema UPWIND dans les zones sous influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      +(1.-anyv3d(i,j,k,id_hybcoefu))*( & !pmx>
           -veldydz_u(i,j,k,id_veltot)* ( & !ppp>
               anyv3d(i,j,k,id_sbef_)+anyv3d(i-1,j,k,id_sbef_)   &
                                ) & !ppp>
       +abs(veldydz_u(i,j,k,id_veltot))*( & !mmm>
               anyv3d(i,j,k,id_sbef_)-anyv3d(i-1,j,k,id_sbef_)   &
!            -(salf_t(i,j,k)         -salf_t(i-1,j,k))*anyv3d(i,j,k,id_hybcoefu) & !05-05-21
                                ) & !mmm>
                                      ) & !pmx>

       ) !prime>

!     -veldydz_u(i,j,k,id_veltot)*dti_lpsub*temref_u(i,j,k)


      enddo !k
      enddo !i
      enddo !j

! Partial Oi Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_sdiv)= & !26-01-17
       anyv3d(i,j,k,id_sdiv)  &
      +anyv3d(i,j,k,id_sbef_)*dti_lpsub*(veldydz_u(i+1,j,k,id_veltot)           &
                                        -veldydz_u(i  ,j,k,id_veltot))/dxdy_t(i,j)

!       anyv3d(i,j,k,id_sbef_)                       &
        anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !ooo> !25-04-17

        anyv3d(i,j,k,id_sbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i+1,j,k,id_flx_)                &
             -anyv3d(i  ,j,k,id_flx_)                &

         +(veldydz_u(i+1,j,k,id_veltot)                      &
          -veldydz_u(i  ,j,k,id_veltot))*anyv3d(i,j,k,id_sbef_) & !26-01-17


                               )/dxdy_t(i,j) ) & !ooo>
                                                /dz_t(i,j,k,before) &!02-11-21 !retablir division par dz si appel A vertmix_merged_levels_any supprimE 

         +(1-mask_t(i,j,k))*sal_t(i,j,k,before)!*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

! Note: le terme anyv3d(i,j,k,id_sdiv) ne fait pas l'objet d'un traitement similaire pour le masque car
! si mask_t=0 alors veldydz_u(i:i+1,...)=0 (question: et les fleuves ?)
      enddo
      enddo
      enddo


      !call vertmix_merged_levels_any(0,id_sbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_sbef_)
      call obc_mpi_anyv3d(1,id_sbef_,'zb')

! Flux Oj
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax+1
      do i=1,imax

      anyv3d(i,j,k,id_fly_)=0.5*( & !prime>

! Schema QUICK en dehors influence des fleuves
      anyv3d(i,j,k,id_hybcoefv)*( &
! Partie advective QUICK:
      -veldxdz_v(i,j,k,id_veltot)  &
        *(    sal_t(i,j,k,now)+sal_t(i,j-1,k,now)  &
                                                 ) &
! Partie diffusive QUICK:
       +cst_*max(anyv3d(i,j,k,id_ofactort),anyv3d(i+1,j,k,id_ofactort)) &
       *(cst2_-anyv3d(i,j,k,id_hybcoefv)*cst1_) & !11-03-21 ! https://docs.google.com/document/d/1QLszDYbzYsIavgI1MgexRbpfEtN6rZDSzRiqdkCsIKE/edit?usp=sharing

                                                                  *(  & !m°v°m>

             (veldxdz_v(i,j,k,id_veltot)+abs(veldxdz_v(i,j,k,id_veltot)))*(  &       !ooo>

            (   anyv3d(i,j  ,k,id_sbef_)-anyv3d(i,j-1,k,id_sbef_)  &
!              -salf_t(i,j  ,k         )+salf_t(i,j-1,k         )  &

! correction anyv3d(i,j  ,k,id_sbef_) par rapport niveau z(i,j-1,k)
      +( anyv3d(i,j  ,kp1,id_sbef_)- anyv3d(i,j  ,km1,id_sbef_)     &
!       -salf_t(i,j  ,kp1         )+ salf_t(i,j  ,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j-1,k  )         -depth_t(i,j  ,k  )         )    &
      /(depth_t(i,j  ,kp1)         -depth_t(i,j  ,km1)         )    &
          )       )      & !minmax>

            ) &

           -(   anyv3d(i,j-1,k,id_sbef_)-anyv3d(i,j-2,k,id_sbef_)  &
!              -salf_t(i,j-1,k         )+salf_t(i,j-2,k         )  &

! correction anyv3d(i,j-2,k,id_sbef_) par rapport niveau z(i,j-1,k)
      -( anyv3d(i,j-2,kp1,id_sbef_)- anyv3d(i,j-2,km1,id_sbef_)     &
!       -salf_t(i,j-2,kp1         )+ salf_t(i,j-2,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j-1,k  )         -depth_t(i,j-2,k  )         )    &
      /(depth_t(i,j-2,kp1)         -depth_t(i,j-2,km1)         )    &
          )       )      & !minmax>

            ) &
                                        *mask_t(i,j-2,kmax)       )& !ooo>

            +(veldxdz_v(i,j,k,id_veltot)-abs(veldxdz_v(i,j,k,id_veltot)))*(  &       !ooo>

            (   anyv3d(i,j+1,k,id_sbef_)-anyv3d(i,j,k,id_sbef_)    &
!              -salf_t(i,j+1,k         )+salf_t(i,j,k         )    &

! correction anyv3d(i,j+1,k,id_sbef_) par rapport niveau z(i,j  ,k)
      +( anyv3d(i,j+1,kp1,id_sbef_)- anyv3d(i,j+1,km1,id_sbef_)     &
!       -salf_t(i,j+1,kp1         )+ salf_t(i,j+1,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j  ,k  )         -depth_t(i,j+1,k  )         )    &
      /(depth_t(i,j+1,kp1)         -depth_t(i,j+1,km1)         )    &
          )       )      & !minmax>

            )   &
               *mask_t(i,j+1,kmax)                                 &

           -(   anyv3d(i,j,k,id_sbef_)-anyv3d(i,j-1,k,id_sbef_)    &
!              -salf_t(i,j,k         )+salf_t(i,j-1,k         )    &

! correction anyv3d(i,j-1,k,id_sbef_) par rapport niveau z(i,j  ,k)
      -( anyv3d(i,j-1,kp1,id_sbef_)- anyv3d(i,j-1,km1,id_sbef_)     &
!       -salf_t(i,j-1,kp1         )+ salf_t(i,j-1,km1         )     &
                                                               )    &
      *min(0.5,max(-0.5, & !minmax>     !22-05-20
       (depth_t(i,j  ,k  )         -depth_t(i,j-1,k  )         )    &
      /(depth_t(i,j-1,kp1)         -depth_t(i,j-1,km1)         )    &
          )       )      & !minmax>

                                                               ))  & !ooo>

                                                                   )  & !m°v°m>()cst_
                                                                   )  & !()anyv3d(i,j,k,id_hybcoefv)
! Schema UPWIND dans les zones sous influence des fleuves:
      +(1.-anyv3d(i,j,k,id_hybcoefv))*( & !pmx>
          -veldxdz_v(i,j,k,id_veltot)* ( & !ppp>
              anyv3d(i,j,k,id_sbef_)+anyv3d(i,j-1,k,id_sbef_)   &
                               ) & !ppp>
      +abs(veldxdz_v(i,j,k,id_veltot))*( & !mmm>
              anyv3d(i,j,k,id_sbef_)-anyv3d(i,j-1,k,id_sbef_)   &
!           -(salf_t(i,j,k)         -salf_t(i,j-1,k))*anyv3d(i,j,k,id_hybcoefv) & !05-05-21
                               ) & !mmm>
                                      ) & !pmx>

       ) !prime>

!      -veldxdz_v(i,j,k,id_veltot)*dti_lpsub*temref_v(i,j,k)

      enddo !k
      enddo !i
      enddo !j

! Partial Horizontal Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_sdiv)= & !26-01-17
       anyv3d(i,j,k,id_sdiv)  &
      +anyv3d(i,j,k,id_sbef_)*dti_lpsub*(veldxdz_v(i,j+1,k,id_veltot)                 &
                                        -veldxdz_v(i,j  ,k,id_veltot))/dxdy_t(i,j)

        anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !ooo>
        anyv3d(i,j,k,id_sbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i,j+1,k,id_fly_)                &
             -anyv3d(i,j  ,k,id_fly_)                &

         +(veldxdz_v(i,j+1,k,id_veltot)                      &
          -veldxdz_v(i,j  ,k,id_veltot))*anyv3d(i,j,k,id_sbef_) & !26-01-17

                              )/dxdy_t(i,j)  ) & !ooo>
                                                /dz_t(i,j,k,before) &!02-11-21 !retablir division par dz si appel A vertmix_merged_levels_any supprimE 
 
         +(1-mask_t(i,j,k))*sal_t(i,j,k,before)!*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

      enddo
      enddo
      enddo


      !call vertmix_merged_levels_any(0,id_sbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_sbef_)
      call obc_mpi_anyv3d(1,id_sbef_,'zb')

!#endif
      enddo               !*********>

      do k=1,kmax ; do j=1,jmax  ; do i=1,imax  

           sal_t(i,j,k,2)                              &
       =  anyv3d(i,j,k,id_sbef_)                    !  &

      enddo ; enddo ; enddo

      end subroutine advection_scal_sal 

!..............................................................................

      subroutine advection_scal_negdif_hori !30-09-21
      use module_principal ; use module_parallele
      implicit none
      integer loop_

! Filtre horizontal 3 points 2 temps 2 directions en une passe
       do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1

        temf_t(i,j,k)=( & !oooo>

       ( (tem_t(i,j,k,1)+tem_t(i,j,k,0))*(1.-0.25*(mask_t(i+1,j,k)+mask_t(i-1,j,k)))+0.25*((tem_t(i+1,j,k,1)+tem_t(i+1,j,k,0))*mask_t(i+1,j,k)+(tem_t(i-1,j,k,1)+tem_t(i-1,j,k,0))*mask_t(i-1,j,k) ) ) &
      *(1.-0.25*(mask_t(i,j+1,k)+mask_t(i,j-1,k))) & 

           +0.25*( & !pmx(   
                   mask_t(i,j+1,k)   &
                   *( (tem_t(i,j+1,k,1)+tem_t(i,j+1,k,0))*(1.-0.25*(mask_t(i+1,j+1,k)+mask_t(i-1,j+1,k)))+0.25*((tem_t(i+1,j+1,k,1)+tem_t(i+1,j+1,k,0))*mask_t(i+1,j+1,k)+(tem_t(i-1,j+1,k,1)+tem_t(i-1,j+1,k,0))*mask_t(i-1,j+1,k) ) )  &
                  +mask_t(i,j-1,k)   &
                   *( (tem_t(i,j-1,k,1)+tem_t(i,j-1,k,0))*(1.-0.25*(mask_t(i+1,j-1,k)+mask_t(i-1,j-1,k)))+0.25*((tem_t(i+1,j-1,k,1)+tem_t(i+1,j-1,k,0))*mask_t(i+1,j-1,k)+(tem_t(i-1,j-1,k,1)+tem_t(i-1,j-1,k,0))*mask_t(i-1,j-1,k) ) )  &
                 ) & !pmx)
                      ) & !oooo>
                       *0.5*ratio_negdif_hor !*0.5 pour 1/2 somme des 2 temps !24-01-22

        salf_t(i,j,k)=( & !oooo>

       ( (sal_t(i,j,k,1)+sal_t(i,j,k,0))*(1.-0.25*(mask_t(i+1,j,k)+mask_t(i-1,j,k)))+0.25*((sal_t(i+1,j,k,1)+sal_t(i+1,j,k,0))*mask_t(i+1,j,k)+(sal_t(i-1,j,k,1)+sal_t(i-1,j,k,0))*mask_t(i-1,j,k) ) ) &
      *(1.-0.25*(mask_t(i,j+1,k)+mask_t(i,j-1,k))) & 

           +0.25*( & !pmx(   
                   mask_t(i,j+1,k)   &
                   *( (sal_t(i,j+1,k,1)+sal_t(i,j+1,k,0))*(1.-0.25*(mask_t(i+1,j+1,k)+mask_t(i-1,j+1,k)))+0.25*((sal_t(i+1,j+1,k,1)+sal_t(i+1,j+1,k,0))*mask_t(i+1,j+1,k)+(sal_t(i-1,j+1,k,1)+sal_t(i-1,j+1,k,0))*mask_t(i-1,j+1,k) ) )  &
                  +mask_t(i,j-1,k)   &
                   *( (sal_t(i,j-1,k,1)+sal_t(i,j-1,k,0))*(1.-0.25*(mask_t(i+1,j-1,k)+mask_t(i-1,j-1,k)))+0.25*((sal_t(i+1,j-1,k,1)+sal_t(i+1,j-1,k,0))*mask_t(i+1,j-1,k)+(sal_t(i-1,j-1,k,1)+sal_t(i-1,j-1,k,0))*mask_t(i-1,j-1,k) ) )  &
                 ) & !pmx)
                      ) & !oooo>
                       *0.5*ratio_negdif_hor !*0.5 pour 1/2 somme des 2 temps !24-01-22


       enddo ; enddo ; enddo

! Conditions limites aux frontieres ouvertes:
       if(obcstatus(ieq1)==1) then !--->
        temf_t(-1,:,:)=temf_t(0,:,:)
        salf_t(-1,:,:)=salf_t(0,:,:)
       endif                       !--->
       if(obcstatus(ieqimax)==1) then !--->
        temf_t(imax+2,:,:)=temf_t(imax+1,:,:)
        salf_t(imax+2,:,:)=salf_t(imax+1,:,:)
       endif                          !--->
       if(obcstatus(jeq1)==1) then !--->
        temf_t(:,-1,:)=temf_t(:,0,:)
        salf_t(:,-1,:)=salf_t(:,0,:)
       endif                       !--->
       if(obcstatus(jeqjmax)==1) then !--->
        temf_t(:,jmax+2,:)=temf_t(:,jmax+1,:)
        salf_t(:,jmax+2,:)=salf_t(:,jmax+1,:)
       endif                          !--->

! Conditions limites mpi:
      call get_type_echange('z2','temf_t_z',anyvar3d,lbound(temf_t),ubound(temf_t),k1)
      call get_type_echange('z2','salf_t_z',anyvar3d,lbound(salf_t),ubound(salf_t),k2)
! Exchanges:
      do loop_=1,subcycle_exchange
        call echange_voisin(temf_t,k1,mpi_neighbor_list(loop_))
        call echange_voisin(salf_t,k2,mpi_neighbor_list(loop_))
      enddo
      call loc_wait() ! ----> important sinon pas d'echanges

      end subroutine advection_scal_negdif_hori

!..............................................................................

      subroutine advection_scal_negdif_vert
      use module_principal ; use module_parallele
      implicit none
      integer loop_

      if(quick_filter_points==3) then !m°v°m> !20-07-22
! Filtre vertical 3 points et 2 temps !30-09-21

       do k=2,kmax-1 ; do j=-1,jmax+2 ; do i=-1,imax+2
! Note: je prends le temps 1 et non pas le temps 0 car dans le cas de
! l'antidiffusion, il me semble que c'est le temps "1" qui va eviter
! l'amplification du mode numerique. On garde quand meme une EMA coef 0.5

        temf_t(i,j,k)=ratio_negdif_ver*(                                    &
                               0.25*(tem_t(i,j,k  ,1)+tem_t(i,j,k  ,0)) &
                             +0.125*(tem_t(i,j,k+1,1)+tem_t(i,j,k+1,0)  &
                                    +tem_t(i,j,k-1,1)+tem_t(i,j,k-1,0)) &
                                   )
        salf_t(i,j,k)=ratio_negdif_ver*(                                    &
                               0.25*(sal_t(i,j,k  ,1)+sal_t(i,j,k  ,0)) &
                             +0.125*(sal_t(i,j,k+1,1)+sal_t(i,j,k+1,0)  &
                                    +sal_t(i,j,k-1,1)+sal_t(i,j,k-1,0)) &
                                   )
       enddo ; enddo ; enddo
! CL surface et fond:
       do k=kmax,kmax+1
        do j=-1,jmax+2 ; do i=-1,imax+2
! Note1 critique: ici par facilite on n'utilise tem_t "non filtrE dans le temps" contrairement
! A l'etape precedente. 
! Note2: on utilise tem_t(:,:,kmax+1,1) qui est a priori correctement defini par advection_scal_vbc
         temf_t(i,j,k)=temf_t(i,j,kmax-1)+(tem_t(i,j,k,1)-tem_t(i,j,kmax-1,1))*ratio_negdif_ver !31-07-21
         salf_t(i,j,k)=salf_t(i,j,kmax-1)+(sal_t(i,j,k,1)-sal_t(i,j,kmax-1,1))*ratio_negdif_ver !31-07-21
        enddo ; enddo
       enddo
       do k=0,1
        do j=-1,jmax+2 ; do i=-1,imax+2
         temf_t(i,j,k)=temf_t(i,j,2)+(tem_t(i,j,k,1)-tem_t(i,j,2,1))*ratio_negdif_ver !31-07-21
         salf_t(i,j,k)=salf_t(i,j,2)+(sal_t(i,j,k,1)-sal_t(i,j,2,1))*ratio_negdif_ver !31-07-21
        enddo ; enddo
       enddo

      return
      endif                           !m°v°m> !20-07-22

      if(quick_filter_points==5) then !m°v°m> !20-07-22
! Filtre 5 points:

       x0=1./16.
       do k=3,kmax-2 ; do j=-1,jmax+2 ; do i=-1,imax+2
! Note: je prends le temps 1 et non pas le temps 0 car dans le cas de
! l'antidiffusion, il me semble que c'est le temps "1" qui va eviter
! l'amplification du mode numerique. On garde quand meme une EMA coef 0.5

        temf_t(i,j,k)=    &
                      0.5 & ! parce que moyenne sur temps 0 et 1 !18-07-22
                     *ratio_negdif_ver*(    (tem_t(i,j,k  ,1)+tem_t(i,j,k  ,0)) &
                                    +x0*(  -(tem_t(i,j,k-2,1)+tem_t(i,j,k-2,0)) &
                                         +4*(tem_t(i,j,k-1,1)+tem_t(i,j,k-1,0)) &
                                         -6*(tem_t(i,j,k  ,1)+tem_t(i,j,k  ,0)) &
                                         +4*(tem_t(i,j,k+1,1)+tem_t(i,j,k+1,0)) &
                                           -(tem_t(i,j,k+2,1)+tem_t(i,j,k+2,0)) &
                                        ))

        salf_t(i,j,k)=    &
                      0.5 & ! parce que moyenne sur temps 0 et 1 !18-07-22
                     *ratio_negdif_ver*(    (sal_t(i,j,k  ,1)+sal_t(i,j,k  ,0)) &
                                    +x0*(  -(sal_t(i,j,k-2,1)+sal_t(i,j,k-2,0)) &
                                         +4*(sal_t(i,j,k-1,1)+sal_t(i,j,k-1,0)) &
                                         -6*(sal_t(i,j,k  ,1)+sal_t(i,j,k  ,0)) &
                                         +4*(sal_t(i,j,k+1,1)+sal_t(i,j,k+1,0)) &
                                           -(sal_t(i,j,k+2,1)+sal_t(i,j,k+2,0)) &
                                        ))
       enddo ; enddo ; enddo
! CL surface et fond:
       do k=kmax-1,kmax+1
        do j=-1,jmax+2 ; do i=-1,imax+2
         temf_t(i,j,k)=temf_t(i,j,kmax-2)+(tem_t(i,j,k,1)-tem_t(i,j,kmax-2,1))*ratio_negdif_ver !18-07-22
         salf_t(i,j,k)=salf_t(i,j,kmax-2)+(sal_t(i,j,k,1)-sal_t(i,j,kmax-2,1))*ratio_negdif_ver !18-07-22
        enddo ; enddo
       enddo
       do k=0,2
        do j=-1,jmax+2 ; do i=-1,imax+2
         temf_t(i,j,k)=temf_t(i,j,3)+(tem_t(i,j,k,1)-tem_t(i,j,3,1))*ratio_negdif_ver !18-07-22
         salf_t(i,j,k)=salf_t(i,j,3)+(sal_t(i,j,k,1)-sal_t(i,j,3,1))*ratio_negdif_ver !18-07-22
        enddo ; enddo
       enddo

      return
      endif                           !m°v°m> !20-07-22

! Pas de diffusion du tout si ratio_negdif_ver=1.
      if(quick_filter_points==0) then !m°v°m>  !06-10-22

       if(ratio_negdif_ver<1.)stop 'Bizarre ratio_negdif_ver<1'       
       if(iadvec_ts_ver/=iadvec_ts_ver_c2) &
       stop 'iadvec_ts_ver/=iadvec_ts_ver_c2'

       do k=1,kmax   ; do j=-1,jmax+2 ; do i=-1,imax+2
        temf_t(i,j,k)= &
         tem_t(i,j,k,2)+ & ! voir note 1
             (ratio_negdif_ver-1.)*(                                    &
                               0.25*(tem_t(i,j,k  ,1)+tem_t(i,j,k  ,0)) &
                             +0.125*(tem_t(i,j,k+1,1)+tem_t(i,j,k+1,0)  &
                                    +tem_t(i,j,k-1,1)+tem_t(i,j,k-1,0)) &
                                   )
        salf_t(i,j,k)= &
         sal_t(i,j,k,2)+ & ! voir note 1
             (ratio_negdif_ver-1.)*(                                    &
                               0.25*(sal_t(i,j,k  ,1)+sal_t(i,j,k  ,0)) &
                             +0.125*(sal_t(i,j,k+1,1)+sal_t(i,j,k+1,0)  &
                                    +sal_t(i,j,k-1,1)+sal_t(i,j,k-1,0)) &
                                   )
! note 1: on prend le temps 2 car la routine est appelee juste avant
! l'advection verticale et que le temps 2 est le "temps present" de
! cette etape. Si ratio_negdif_ver=1. la partie diffusion negative est
! nulle et tem_t(i,j,k,2)-temf_t(i,j,k) sera strictement nul dans
! l'etape de l'advection verticale
       enddo ; enddo ; enddo
      return
      endif                           !m°v°m>

      stop 'Err 4128 quick_filter_points undefined'

      end subroutine advection_scal_negdif_vert

!..............................................................................
!................................................................
!..............................................................................
!..............................................................................

      subroutine advection_scal_tem_hornegdif  !30-09-21
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      double precision :: velef_,flagmsk_=1.,scalefac_,cst_,cst2_=2.,cst1_=1.
      integer :: t_=1        &
                ,id_flx_=1   &  ! x flux identifier
                ,id_fly_=1   &  ! y flux identifier
                ,id_tbef_=3  &  ! temperature    "before" identifier
!               ,id_tdiv=4   &  ! T*divergence(velocity) !26-01-17
                ,loop_

! https://docs.google.com/document/d/1QLszDYbzYsIavgI1MgexRbpfEtN6rZDSzRiqdkCsIKE/edit?usp=sharing

!     cst_=1./6. ! ou 1./8.
      cst_=quick_coef !18-07-22
      if(iadvec_ts_hor==iadvec_ts_hor_quickest2) then !29-11-21
! Quickest-2
!        cst_=1./8. ! ou 1./8.
         cst2_=1.
         cst1_=0.
      else
! Quickest
         cst2_=2.
         cst1_=1.
      endif

! A quoi sert cst2_ dans l'expression (cst2-anyv3d(i,j,k,id_hybcoefu)*cst1) ?
! Elle sert a choisir entre QUICK et QUICKEST
! vaut 2-anyv3d(i,j,k,id_hybcoefu)=1+nc pour QUICKEST et vaut 1 pour QUICK



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
      -veldydz_u(i,j,k,id_velexp) &
        *(    tem_t(i,j,k,now)+tem_t(i-1,j,k,now)      &
                                                 )     &
! Partie diffusive QUICK remise A jour aux sous-pas de temps advectif:
       +cst_ &!*max(anyv3d(i,j,k,id_ofactort),anyv3d(i,j+1,k,id_ofactort)) &
       *(cst2_-anyv3d(i,j,k,id_hybcoefu)*cst1_) & !11-03-21 ! https://docs.google.com/document/d/1QLszDYbzYsIavgI1MgexRbpfEtN6rZDSzRiqdkCsIKE/edit?usp=sharing

                                                               *( & !m°v°m>
             (veldydz_u(i,j,k,id_velexp)+abs(veldydz_u(i,j,k,id_velexp)))*( &           !pmxpmx>

               (   anyv3d(i  ,j,k,id_tbef_)   -anyv3d(i-1,j,k,id_tbef_)  &
                  -temf_t(i  ,j,k         )   +temf_t(i-1,j,k         )  &

! correction anyv3d(i  ,j,k,id_tbef_) par rapport niveau z(i-1,j,k)
!     +( anyv3d(i  ,j,kp1,id_tbef_)- anyv3d(i  ,j,km1,id_tbef_)     &
!       -temf_t(i  ,j,kp1         )+ temf_t(i  ,j,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i-1,j,k  )         -depth_t(i  ,j,k  )         )    &
!     /(depth_t(i  ,j,kp1)         -depth_t(i  ,j,km1)         )    &
!         )       )      & !minmax>

               ) &

              -(   anyv3d(i-1,j,k,id_tbef_)-anyv3d(i-2,j,k,id_tbef_)  &
                  -temf_t(i-1,j,k         )+temf_t(i-2,j,k         )  &

! correction anyv3d(i-2,j,k,id_tbef_) par rapport niveau z(i-1,j,k)
!     -( anyv3d(i-2,j,kp1,id_tbef_)- anyv3d(i-2,j,km1,id_tbef_)     &
!       -temf_t(i-2,j,kp1         )+ temf_t(i-2,j,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i-1,j,k  )         -depth_t(i-2,j,k  )         )    &
!     /(depth_t(i-2,j,kp1)         -depth_t(i-2,j,km1)         )    &
!         )       )      & !minmax>

               ) &
                                           *mask_t(i-2,j,kmax)      ) & !pmxpmx>

            +(veldydz_u(i,j,k,id_velexp)-abs(veldydz_u(i,j,k,id_velexp)))*( &           !pmxpmx>

               (   anyv3d(i+1,j,k,id_tbef_)-anyv3d(i,j,k,id_tbef_)    &
                  -temf_t(i+1,j,k         )+temf_t(i,j,k         )    &

! correction anyv3d(i+1,j,k,id_tbef_) par rapport niveau z(i  ,j,k)
!     +( anyv3d(i+1,j,kp1,id_tbef_)- anyv3d(i+1,j,km1,id_tbef_)     &
!       -temf_t(i+1,j,kp1         )+ temf_t(i+1,j,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i  ,j,k  )         -depth_t(i+1,j,k  )         )    &
!     /(depth_t(i+1,j,kp1)         -depth_t(i+1,j,km1)         )    &
!         )       )      & !minmax>

               )   &
                  *mask_t(i+1,j,kmax)                                 &

              -(   anyv3d(i,j,k,id_tbef_)-anyv3d(i-1,j,k,id_tbef_)    &
                  -temf_t(i,j,k         )+temf_t(i-1,j,k         )    &

! correction anyv3d(i-1,j,k,id_tbef_) par rapport niveau z(i  ,j,k)
!     -( anyv3d(i-1,j,kp1,id_tbef_)- anyv3d(i-1,j,km1,id_tbef_)     &
!       -temf_t(i-1,j,kp1         )+ temf_t(i-1,j,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i  ,j,k  )         -depth_t(i-1,j,k  )         )    &
!     /(depth_t(i-1,j,kp1)         -depth_t(i-1,j,km1)         )    &
!         )       )      & !minmax>

                                                                  ))  & !pmxpmx>

                                                                )  & !m°v°m>()cst_
                                                                )  & !()anyv3d(i,j,k,id_hybcoefu)
! Schema UPWIND dans les zones sous influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      +(1.-anyv3d(i,j,k,id_hybcoefu))*( & !pmx>
           -veldydz_u(i,j,k,id_velexp)* ( & !ppp>
               anyv3d(i,j,k,id_tbef_)+anyv3d(i-1,j,k,id_tbef_)   &
                                ) & !ppp>
       +abs(veldydz_u(i,j,k,id_velexp))*( & !mmm>
               anyv3d(i,j,k,id_tbef_)-anyv3d(i-1,j,k,id_tbef_)   &
             -(temf_t(i,j,k)         -temf_t(i-1,j,k))*anyv3d(i,j,k,id_hybcoefu) & !05-05-21
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
       anyv3d(i,j,k,id_tdiv)= & !26-01-17
       anyv3d(i,j,k,id_tdiv)  &
      +anyv3d(i,j,k,id_tbef_)*dti_lpsub*(veldydz_u(i+1,j,k,id_velexp)           &
                                        -veldydz_u(i  ,j,k,id_velexp))/dxdy_t(i,j)

!       anyv3d(i,j,k,id_tbef_)                       &
        anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !ooo> !25-04-17

        anyv3d(i,j,k,id_tbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i+1,j,k,id_flx_)                &
             -anyv3d(i  ,j,k,id_flx_)                &

         +(veldydz_u(i+1,j,k,id_velexp)                      &
          -veldydz_u(i  ,j,k,id_velexp))*anyv3d(i,j,k,id_tbef_) & !26-01-17


                               )/dxdy_t(i,j) ) & !ooo>
                                                /dz_t(i,j,k,before) &!02-11-21 !retablir division par dz si appel A vertmix_merged_levels_any supprimE 

         +(1-mask_t(i,j,k))*tem_t(i,j,k,before)!*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

! Note: le terme anyv3d(i,j,k,id_tdiv) ne fait pas l'objet d'un traitement similaire pour le masque car
! si mask_t=0 alors veldydz_u(i:i+1,...)=0 (question: et les fleuves ?)
      enddo
      enddo
      enddo


      !call vertmix_merged_levels_any(0,id_tbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_tbef_)
      call obc_mpi_anyv3d(1,id_tbef_,'zb')

! Flux Oj
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax+1
      do i=1,imax

      anyv3d(i,j,k,id_fly_)=0.5*( & !prime>

! Schema QUICK en dehors influence des fleuves
      anyv3d(i,j,k,id_hybcoefv)*( &
! Partie advective QUICK:
      -veldxdz_v(i,j,k,id_velexp)  &
        *(    tem_t(i,j,k,now)+tem_t(i,j-1,k,now)  &
                                                 ) &
! Partie diffusive QUICK:
       +cst_ &!*max(anyv3d(i,j,k,id_ofactort),anyv3d(i+1,j,k,id_ofactort)) &
       *(cst2_-anyv3d(i,j,k,id_hybcoefv)*cst1_) & !11-03-21 ! https://docs.google.com/document/d/1QLszDYbzYsIavgI1MgexRbpfEtN6rZDSzRiqdkCsIKE/edit?usp=sharing

                                                                  *(  & !m°v°m>

             (veldxdz_v(i,j,k,id_velexp)+abs(veldxdz_v(i,j,k,id_velexp)))*(  &       !ooo>

            (   anyv3d(i,j  ,k,id_tbef_)-anyv3d(i,j-1,k,id_tbef_)  &
               -temf_t(i,j  ,k         )+temf_t(i,j-1,k         )  &

! correction anyv3d(i,j  ,k,id_tbef_) par rapport niveau z(i,j-1,k)
!     +( anyv3d(i,j  ,kp1,id_tbef_)- anyv3d(i,j  ,km1,id_tbef_)     &
!       -temf_t(i,j  ,kp1         )+ temf_t(i,j  ,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i,j-1,k  )         -depth_t(i,j  ,k  )         )    &
!     /(depth_t(i,j  ,kp1)         -depth_t(i,j  ,km1)         )    &
!         )       )      & !minmax>

            ) &

           -(   anyv3d(i,j-1,k,id_tbef_)-anyv3d(i,j-2,k,id_tbef_)  &
               -temf_t(i,j-1,k         )+temf_t(i,j-2,k         )  &

! correction anyv3d(i,j-2,k,id_tbef_) par rapport niveau z(i,j-1,k)
!     -( anyv3d(i,j-2,kp1,id_tbef_)- anyv3d(i,j-2,km1,id_tbef_)     &
!       -temf_t(i,j-2,kp1         )+ temf_t(i,j-2,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i,j-1,k  )         -depth_t(i,j-2,k  )         )    &
!     /(depth_t(i,j-2,kp1)         -depth_t(i,j-2,km1)         )    &
!         )       )      & !minmax>

            ) &
                                        *mask_t(i,j-2,kmax)       )& !ooo>

            +(veldxdz_v(i,j,k,id_velexp)-abs(veldxdz_v(i,j,k,id_velexp)))*(  &       !ooo>

            (   anyv3d(i,j+1,k,id_tbef_)-anyv3d(i,j,k,id_tbef_)    &
               -temf_t(i,j+1,k         )+temf_t(i,j,k         )    &

! correction anyv3d(i,j+1,k,id_tbef_) par rapport niveau z(i,j  ,k)
!     +( anyv3d(i,j+1,kp1,id_tbef_)- anyv3d(i,j+1,km1,id_tbef_)     &
!       -temf_t(i,j+1,kp1         )+ temf_t(i,j+1,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i,j  ,k  )         -depth_t(i,j+1,k  )         )    &
!     /(depth_t(i,j+1,kp1)         -depth_t(i,j+1,km1)         )    &
!         )       )      & !minmax>

            )   &
               *mask_t(i,j+1,kmax)                                 &

           -(   anyv3d(i,j,k,id_tbef_)-anyv3d(i,j-1,k,id_tbef_)    &
               -temf_t(i,j,k         )+temf_t(i,j-1,k         )    &

! correction anyv3d(i,j-1,k,id_tbef_) par rapport niveau z(i,j  ,k)
!     -( anyv3d(i,j-1,kp1,id_tbef_)- anyv3d(i,j-1,km1,id_tbef_)     &
!       -temf_t(i,j-1,kp1         )+ temf_t(i,j-1,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i,j  ,k  )         -depth_t(i,j-1,k  )         )    &
!     /(depth_t(i,j-1,kp1)         -depth_t(i,j-1,km1)         )    &
!         )       )      & !minmax>

                                                               ))  & !ooo>

                                                                   )  & !m°v°m>()cst_
                                                                   )  & !()anyv3d(i,j,k,id_hybcoefv)
! Schema UPWIND dans les zones sous influence des fleuves:
      +(1.-anyv3d(i,j,k,id_hybcoefv))*( & !pmx>
          -veldxdz_v(i,j,k,id_velexp)* ( & !ppp>
              anyv3d(i,j,k,id_tbef_)+anyv3d(i,j-1,k,id_tbef_)   &
                               ) & !ppp>
      +abs(veldxdz_v(i,j,k,id_velexp))*( & !mmm>
              anyv3d(i,j,k,id_tbef_)-anyv3d(i,j-1,k,id_tbef_)   &
            -(temf_t(i,j,k)         -temf_t(i,j-1,k))*anyv3d(i,j,k,id_hybcoefv) & !05-05-21
                               ) & !mmm>
                                      ) & !pmx>

       ) !prime>

!      -veldxdz_v(i,j,k,id_velexp)*dti_lpsub*temref_v(i,j,k)

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
      +anyv3d(i,j,k,id_tbef_)*dti_lpsub*(veldxdz_v(i,j+1,k,id_velexp)                 &
                                        -veldxdz_v(i,j  ,k,id_velexp))/dxdy_t(i,j)

        anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !ooo>
        anyv3d(i,j,k,id_tbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i,j+1,k,id_fly_)                &
             -anyv3d(i,j  ,k,id_fly_)                &

         +(veldxdz_v(i,j+1,k,id_velexp)                      &
          -veldxdz_v(i,j  ,k,id_velexp))*anyv3d(i,j,k,id_tbef_) & !26-01-17

                              )/dxdy_t(i,j)  ) & !ooo>
                                                /dz_t(i,j,k,before) &!02-11-21 !retablir division par dz si appel A vertmix_merged_levels_any supprimE 
 
         +(1-mask_t(i,j,k))*tem_t(i,j,k,before)!*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

      enddo
      enddo
      enddo


      !call vertmix_merged_levels_any(0,id_tbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_tbef_)
      call obc_mpi_anyv3d(1,id_tbef_,'zb')

!#endif
      enddo               !*********>

      do k=1,kmax ; do j=1,jmax  ; do i=1,imax  

           tem_t(i,j,k,2)                              &
       =  anyv3d(i,j,k,id_tbef_)                    !  &

      enddo ; enddo ; enddo

      end subroutine advection_scal_tem_hornegdif 


!..............................................................................

      subroutine advection_scal_sal_hornegdif  !30-09-21
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      double precision :: velef_,flagmsk_=1.,scalefac_,cst_,cst2_=2.,cst1_=1.
      integer :: t_=1        &
                ,id_flx_=1   &  ! x flux identifier
                ,id_fly_=1   &  ! y flux identifier
                ,id_sbef_=3  &  ! temperature    "before" identifier
!               ,id_sdiv=5   &  ! S*divergence(velocity) !26-01-17
                ,loop_

! https://docs.google.com/document/d/1QLszDYbzYsIavgI1MgexRbpfEtN6rZDSzRiqdkCsIKE/edit?usp=sharing

!     cst_=1./6. ! ou 1./8.
      cst_=quick_coef !18-07-22
      if(iadvec_ts_hor==iadvec_ts_hor_quickest2) then
! Quickest-2
!        cst_=1./8. ! ou 1./8.
         cst2_=1.
         cst1_=0.
      else
! Quickest
         cst2_=2.
         cst1_=1.
      endif

! A quoi sert cst2_ dans l'expression (cst2-anyv3d(i,j,k,id_hybcoefu)*cst1) ?
! Elle sert a choisir entre QUICK et QUICKEST
! vaut 2-anyv3d(i,j,k,id_hybcoefu)=1+nc pour QUICKEST et vaut 1 pour QUICK



      do k=1,kmax ; do j=-1,jmax+2 ; do i=-1,imax+2
       anyv3d(i,j,k,id_sbef_)=sal_t(i,j,k,before)
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,id_sdiv)=0. ! reset cumul temperature*divergence(courant) !26-01-17
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
      -veldydz_u(i,j,k,id_velexp) &
        *(    sal_t(i,j,k,now)+sal_t(i-1,j,k,now)      &
                                                 )     &
! Partie diffusive QUICK remise A jour aux sous-pas de temps advectif:
       +cst_ &!*max(anyv3d(i,j,k,id_ofactort),anyv3d(i,j+1,k,id_ofactort)) &
       *(cst2_-anyv3d(i,j,k,id_hybcoefu)*cst1_) & !11-03-21 ! https://docs.google.com/document/d/1QLszDYbzYsIavgI1MgexRbpfEtN6rZDSzRiqdkCsIKE/edit?usp=sharing

                                                               *( & !m°v°m>
             (veldydz_u(i,j,k,id_velexp)+abs(veldydz_u(i,j,k,id_velexp)))*( &           !pmxpmx>

               (   anyv3d(i  ,j,k,id_sbef_)   -anyv3d(i-1,j,k,id_sbef_)  &
                  -salf_t(i  ,j,k         )   +salf_t(i-1,j,k         )  &

! correction anyv3d(i  ,j,k,id_sbef_) par rapport niveau z(i-1,j,k)
!     +( anyv3d(i  ,j,kp1,id_sbef_)- anyv3d(i  ,j,km1,id_sbef_)     &
!       -salf_t(i  ,j,kp1         )+ salf_t(i  ,j,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i-1,j,k  )         -depth_t(i  ,j,k  )         )    &
!     /(depth_t(i  ,j,kp1)         -depth_t(i  ,j,km1)         )    &
!         )       )      & !minmax>

               ) &

              -(   anyv3d(i-1,j,k,id_sbef_)-anyv3d(i-2,j,k,id_sbef_)  &
                  -salf_t(i-1,j,k         )+salf_t(i-2,j,k         )  &

! correction anyv3d(i-2,j,k,id_sbef_) par rapport niveau z(i-1,j,k)
!     -( anyv3d(i-2,j,kp1,id_sbef_)- anyv3d(i-2,j,km1,id_sbef_)     &
!       -salf_t(i-2,j,kp1         )+ salf_t(i-2,j,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i-1,j,k  )         -depth_t(i-2,j,k  )         )    &
!     /(depth_t(i-2,j,kp1)         -depth_t(i-2,j,km1)         )    &
!         )       )      & !minmax>

               ) &
                                           *mask_t(i-2,j,kmax)      ) & !pmxpmx>

            +(veldydz_u(i,j,k,id_velexp)-abs(veldydz_u(i,j,k,id_velexp)))*( &           !pmxpmx>

               (   anyv3d(i+1,j,k,id_sbef_)-anyv3d(i,j,k,id_sbef_)    &
                  -salf_t(i+1,j,k         )+salf_t(i,j,k         )    &

! correction anyv3d(i+1,j,k,id_sbef_) par rapport niveau z(i  ,j,k)
!     +( anyv3d(i+1,j,kp1,id_sbef_)- anyv3d(i+1,j,km1,id_sbef_)     &
!       -salf_t(i+1,j,kp1         )+ salf_t(i+1,j,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i  ,j,k  )         -depth_t(i+1,j,k  )         )    &
!     /(depth_t(i+1,j,kp1)         -depth_t(i+1,j,km1)         )    &
!         )       )      & !minmax>

               )   &
                  *mask_t(i+1,j,kmax)                                 &

              -(   anyv3d(i,j,k,id_sbef_)-anyv3d(i-1,j,k,id_sbef_)    &
                  -salf_t(i,j,k         )+salf_t(i-1,j,k         )    &

! correction anyv3d(i-1,j,k,id_sbef_) par rapport niveau z(i  ,j,k)
!     -( anyv3d(i-1,j,kp1,id_sbef_)- anyv3d(i-1,j,km1,id_sbef_)     &
!       -salf_t(i-1,j,kp1         )+ salf_t(i-1,j,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i  ,j,k  )         -depth_t(i-1,j,k  )         )    &
!     /(depth_t(i-1,j,kp1)         -depth_t(i-1,j,km1)         )    &
!         )       )      & !minmax>

                                                                  ))  & !pmxpmx>

                                                                )  & !m°v°m>()cst_
                                                                )  & !()anyv3d(i,j,k,id_hybcoefu)
! Schema UPWIND dans les zones sous influence des fleuves (selection par anyv3d(i,j,k,id_hybcoefu))
      +(1.-anyv3d(i,j,k,id_hybcoefu))*( & !pmx>
           -veldydz_u(i,j,k,id_velexp)* ( & !ppp>
               anyv3d(i,j,k,id_sbef_)+anyv3d(i-1,j,k,id_sbef_)   &
                                ) & !ppp>
       +abs(veldydz_u(i,j,k,id_velexp))*( & !mmm>
               anyv3d(i,j,k,id_sbef_)-anyv3d(i-1,j,k,id_sbef_)   &
             -(salf_t(i,j,k)         -salf_t(i-1,j,k))*anyv3d(i,j,k,id_hybcoefu) & !05-05-21
                                ) & !mmm>
                                      ) & !pmx>

       ) !prime>

!     -veldydz_u(i,j,k,id_velexp)*dti_lpsub*temref_u(i,j,k)


      enddo !k
      enddo !i
      enddo !j

! Partial Oi Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_sdiv)= & !26-01-17
       anyv3d(i,j,k,id_sdiv)  &
      +anyv3d(i,j,k,id_sbef_)*dti_lpsub*(veldydz_u(i+1,j,k,id_velexp)           &
                                        -veldydz_u(i  ,j,k,id_velexp))/dxdy_t(i,j)

!       anyv3d(i,j,k,id_sbef_)                       &
        anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !ooo> !25-04-17

        anyv3d(i,j,k,id_sbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i+1,j,k,id_flx_)                &
             -anyv3d(i  ,j,k,id_flx_)                &

         +(veldydz_u(i+1,j,k,id_velexp)                      &
          -veldydz_u(i  ,j,k,id_velexp))*anyv3d(i,j,k,id_sbef_) & !26-01-17


                               )/dxdy_t(i,j) ) & !ooo>
                                                /dz_t(i,j,k,before) &!02-11-21 !retablir division par dz si appel A vertmix_merged_levels_any supprimE 

         +(1-mask_t(i,j,k))*sal_t(i,j,k,before)!*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

! Note: le terme anyv3d(i,j,k,id_sdiv) ne fait pas l'objet d'un traitement similaire pour le masque car
! si mask_t=0 alors veldydz_u(i:i+1,...)=0 (question: et les fleuves ?)
      enddo
      enddo
      enddo


      !call vertmix_merged_levels_any(0,id_sbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_sbef_)
      call obc_mpi_anyv3d(1,id_sbef_,'zb')

! Flux Oj
      do k=1,kmax
      km1=max(k-1,1) ; kp1=min(k+1,kmax)
      do j=1,jmax+1
      do i=1,imax

      anyv3d(i,j,k,id_fly_)=0.5*( & !prime>

! Schema QUICK en dehors influence des fleuves
      anyv3d(i,j,k,id_hybcoefv)*( &
! Partie advective QUICK:
      -veldxdz_v(i,j,k,id_velexp)  &
        *(    sal_t(i,j,k,now)+sal_t(i,j-1,k,now)  &
                                                 ) &
! Partie diffusive QUICK:
       +cst_ &!*max(anyv3d(i,j,k,id_ofactort),anyv3d(i+1,j,k,id_ofactort)) &
       *(cst2_-anyv3d(i,j,k,id_hybcoefv)*cst1_) & !11-03-21 ! https://docs.google.com/document/d/1QLszDYbzYsIavgI1MgexRbpfEtN6rZDSzRiqdkCsIKE/edit?usp=sharing

                                                                  *(  & !m°v°m>

             (veldxdz_v(i,j,k,id_velexp)+abs(veldxdz_v(i,j,k,id_velexp)))*(  &       !ooo>

            (   anyv3d(i,j  ,k,id_sbef_)-anyv3d(i,j-1,k,id_sbef_)  &
               -salf_t(i,j  ,k         )+salf_t(i,j-1,k         )  &

! correction anyv3d(i,j  ,k,id_sbef_) par rapport niveau z(i,j-1,k)
!     +( anyv3d(i,j  ,kp1,id_sbef_)- anyv3d(i,j  ,km1,id_sbef_)     &
!       -salf_t(i,j  ,kp1         )+ salf_t(i,j  ,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i,j-1,k  )         -depth_t(i,j  ,k  )         )    &
!     /(depth_t(i,j  ,kp1)         -depth_t(i,j  ,km1)         )    &
!         )       )      & !minmax>

            ) &

           -(   anyv3d(i,j-1,k,id_sbef_)-anyv3d(i,j-2,k,id_sbef_)  &
               -salf_t(i,j-1,k         )+salf_t(i,j-2,k         )  &

! correction anyv3d(i,j-2,k,id_sbef_) par rapport niveau z(i,j-1,k)
!     -( anyv3d(i,j-2,kp1,id_sbef_)- anyv3d(i,j-2,km1,id_sbef_)     &
!       -salf_t(i,j-2,kp1         )+ salf_t(i,j-2,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i,j-1,k  )         -depth_t(i,j-2,k  )         )    &
!     /(depth_t(i,j-2,kp1)         -depth_t(i,j-2,km1)         )    &
!         )       )      & !minmax>

            ) &
                                        *mask_t(i,j-2,kmax)       )& !ooo>

            +(veldxdz_v(i,j,k,id_velexp)-abs(veldxdz_v(i,j,k,id_velexp)))*(  &       !ooo>

            (   anyv3d(i,j+1,k,id_sbef_)-anyv3d(i,j,k,id_sbef_)    &
               -salf_t(i,j+1,k         )+salf_t(i,j,k         )    &

! correction anyv3d(i,j+1,k,id_sbef_) par rapport niveau z(i,j  ,k)
!     +( anyv3d(i,j+1,kp1,id_sbef_)- anyv3d(i,j+1,km1,id_sbef_)     &
!       -salf_t(i,j+1,kp1         )+ salf_t(i,j+1,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i,j  ,k  )         -depth_t(i,j+1,k  )         )    &
!     /(depth_t(i,j+1,kp1)         -depth_t(i,j+1,km1)         )    &
!         )       )      & !minmax>

            )   &
               *mask_t(i,j+1,kmax)                                 &

           -(   anyv3d(i,j,k,id_sbef_)-anyv3d(i,j-1,k,id_sbef_)    &
               -salf_t(i,j,k         )+salf_t(i,j-1,k         )    &

! correction anyv3d(i,j-1,k,id_sbef_) par rapport niveau z(i,j  ,k)
!     -( anyv3d(i,j-1,kp1,id_sbef_)- anyv3d(i,j-1,km1,id_sbef_)     &
!       -salf_t(i,j-1,kp1         )+ salf_t(i,j-1,km1         )     &
!                                                              )    &
!     *min(0.5,max(-0.5, & !minmax>     !22-05-20
!      (depth_t(i,j  ,k  )         -depth_t(i,j-1,k  )         )    &
!     /(depth_t(i,j-1,kp1)         -depth_t(i,j-1,km1)         )    &
!         )       )      & !minmax>

                                                               ))  & !ooo>

                                                                   )  & !m°v°m>()cst_
                                                                   )  & !()anyv3d(i,j,k,id_hybcoefv)
! Schema UPWIND dans les zones sous influence des fleuves:
      +(1.-anyv3d(i,j,k,id_hybcoefv))*( & !pmx>
          -veldxdz_v(i,j,k,id_velexp)* ( & !ppp>
              anyv3d(i,j,k,id_sbef_)+anyv3d(i,j-1,k,id_sbef_)   &
                               ) & !ppp>
      +abs(veldxdz_v(i,j,k,id_velexp))*( & !mmm>
              anyv3d(i,j,k,id_sbef_)-anyv3d(i,j-1,k,id_sbef_)   &
            -(salf_t(i,j,k)         -salf_t(i,j-1,k))*anyv3d(i,j,k,id_hybcoefv) & !05-05-21
                               ) & !mmm>
                                      ) & !pmx>

       ) !prime>

!      -veldxdz_v(i,j,k,id_velexp)*dti_lpsub*temref_v(i,j,k)

      enddo !k
      enddo !i
      enddo !j

! Partial Horizontal Advection
      do k=1,kmax
      do j=1,jmax
      do i=1,imax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_sdiv)= & !26-01-17
       anyv3d(i,j,k,id_sdiv)  &
      +anyv3d(i,j,k,id_sbef_)*dti_lpsub*(veldxdz_v(i,j+1,k,id_velexp)                 &
                                        -veldxdz_v(i,j  ,k,id_velexp))/dxdy_t(i,j)

        anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !ooo>
        anyv3d(i,j,k,id_sbef_)*dz_t(i,j,k,before)    & !01-05-19

       +dti_lpsub*wetmask_t(i,j)*(                   &

              anyv3d(i,j+1,k,id_fly_)                &
             -anyv3d(i,j  ,k,id_fly_)                &

         +(veldxdz_v(i,j+1,k,id_velexp)                      &
          -veldxdz_v(i,j  ,k,id_velexp))*anyv3d(i,j,k,id_sbef_) & !26-01-17

                              )/dxdy_t(i,j)  ) & !ooo>
                                                /dz_t(i,j,k,before) &!02-11-21 !retablir division par dz si appel A vertmix_merged_levels_any supprimE 
 
         +(1-mask_t(i,j,k))*sal_t(i,j,k,before)!*dz_t(i,j,k,before) ! Cette ligne pour ne pas perdre la valeur du champs dans le masque

      enddo
      enddo
      enddo


      !call vertmix_merged_levels_any(0,id_sbef_) !07-12-17

! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_sbef_)
      call obc_mpi_anyv3d(1,id_sbef_,'zb')

!#endif
      enddo               !*********>

      do k=1,kmax ; do j=1,jmax  ; do i=1,imax  

           sal_t(i,j,k,2)                              &
       =  anyv3d(i,j,k,id_sbef_)                    !  &

      enddo ; enddo ; enddo

      end subroutine advection_scal_sal_hornegdif 

!..............................................................................

      subroutine advection_scal_tem_horimpli
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      integer :: t_=1        &
                ,id_tbef_=3  &  ! temperature    "before" identifier
                ,loop_       &
                ,flag_mpi_=0 &
                ,work_status_=0     & !0=rien A signaler, 1=travail demandE       , 2=travail fait
                ,send_status_=0     & !0=rien A signaler, 1=envoi mpi demandE     , 2=envoi mpi rEalisE
                ,recv_status_=0     & !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
                ,status_voisin_upst & ! 0= pas d'echange avec voisin amont, 1= echange
                ,status_voisin_down   ! 0= pas d'echange avec voisin aval , 1= echange




  
      if(iperiodicboundary==.true.)stop 'conditions periodiques A faire'
      if(jperiodicboundary==.true.)stop 'conditions periodiques A faire'

! flag_mpi_=0: je ne fais rien
! flag_mpi_=1: je fais le calcul
! flag_mpi_=2: j'ai fait le calcul, j'envoie la valeur au voisin
! flag_mpi_=3: j'ai tout fait, je m'arrete

      do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
!      anyv3d(i,j,k,id_tbef_)=tem_t(i,j,k,before) ! pour schema 100% implicite
       anyv3d(i,j,k,id_tbef_)=tem_t(i,j,k,2)      ! pour schema hybride explicite/implicite
      enddo ; enddo ; enddo

! Schema 100% implicite decommenter les lignes suivantes, sinon le cumul de anyv3d(i,j,k,id_tdiv) se poursuit
!     do k=1,kmax ; do j=1,jmax ; do i=1,imax 
!      anyv3d(i,j,k,id_tdiv)=0. ! reset cumul temperature*divergence(courant) !26-01-17
!     enddo ; enddo ; enddo

!...........................
! Advection axe Oi
!...........................

! Passage boucle i croissante pour la partie positive de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0       ! 0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0       ! 0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0       ! 0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement positive (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(maxval(veldydz_u(1,1:jmax,1:kmax,id_velimp))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_u en i=imax-1 du present rank est vel_u en i=1 du voisin aval
      if(maxval(veldydz_u(imax-1,1:jmax,1:kmax,id_velimp))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(veldydz_u(1:imax,1:jmax,1:kmax,id_velimp))<=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                          !>>>  

      do loop_=1,nbdom_imax !-boucle-mpi->

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(est)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(ouest)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

         call advection_scal_impobc(1,id_tbef_,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

!      anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       *anyv3d(i,j,k,id_tbef_)                              &
!       +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_lp &
!       *anyv3d(i-1,j,k,id_tbef_)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_lp &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !LSM>

         anyv3d(i,j,k,id_tbef_) &

        +( & !haut>

         (veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_lp &
        *(anyv3d(i-1,j,k,id_tbef_)-anyv3d(i,j,k,id_tbef_))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(1,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

       anyv3d(i,j,k,id_tdiv)= & 
       anyv3d(i,j,k,id_tdiv)  &
      +anyv3d(i,j,k,id_tbef_)*dti_lp*0.5*( &
                                          veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)) &
                                         -veldydz_u(i  ,j,k,id_velimp)-abs(veldydz_u(i  ,j,k,id_velimp)) &
                                         )/dxdy_t(i,j)
      enddo ; enddo ; enddo

! Passage boucle i decroissante pour la partie negative de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement negative (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(minval(veldydz_u(imax+1,1:jmax,1:kmax,id_velimp))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_u en i=3 du present rank est vel_u en i=imax+1 du voisin aval
      if(minval(veldydz_u(3,1:jmax,1:kmax,id_velimp))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(veldydz_u(2:imax+1,1:jmax,1:kmax,id_velimp))>=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                            !>>>  

      do loop_=1,nbdom_imax !-boucle-mpi->

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(ouest)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(est)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

      call advection_scal_impobc(2,id_tbef_,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=1,jmax ; do i=imax,1,-1

!      anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       *anyv3d(i,j,k,id_tbef_)                              &
!       -(veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)))*dti_lp &
!       *anyv3d(i+1,j,k,id_tbef_)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       -(veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)))*dti_lp &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !LSM>

         anyv3d(i,j,k,id_tbef_) &

        +( & !haut>

         (-veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)))*dti_lp &
        *(anyv3d(i+1,j,k,id_tbef_)-anyv3d(i,j,k,id_tbef_))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +(-veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(2,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

       anyv3d(i,j,k,id_tdiv)= & 
       anyv3d(i,j,k,id_tdiv)  &
      +anyv3d(i,j,k,id_tbef_)*dti_lp*0.5*( &
                                          veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)) &
                                         -veldydz_u(i  ,j,k,id_velimp)+abs(veldydz_u(i  ,j,k,id_velimp)) &
                                         )/dxdy_t(i,j)
      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans advection_scal_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_tbef_)
!     call obc_mpi_anyv3d(1,id_tbef_,'zb')

!...........................
! Advection axe Oj
!...........................

! Passage boucle j croissante pour la partie positive de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement positive (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(maxval(veldxdz_v(1:imax,1,1:kmax,id_velimp))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                 status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_v en j=jmax-1 du present rank est vel_v en j=1 du voisin aval
      if(maxval(veldxdz_v(1:imax,jmax-1,1:kmax,id_velimp))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(veldxdz_v(1:imax,1:jmax,1:kmax,id_velimp))<=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                          !>>>  

      do loop_=1,nbdom_jmax !-boucle-mpi->

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(nord)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(sud)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

         call advection_scal_impobc(3,id_tbef_,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

!      anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       *anyv3d(i,j,k,id_tbef_)                              &
!       +(veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_lp &
!       *anyv3d(i,j-1,k,id_tbef_)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       +(veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_lp &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !LSM>

         anyv3d(i,j,k,id_tbef_) &

        +( & !haut>

         ( veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_lp &
        *(anyv3d(i,j-1,k,id_tbef_)-anyv3d(i,j,k,id_tbef_))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +( veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(3,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1


      enddo                 !-boucle-mpi->

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

       anyv3d(i,j,k,id_tdiv)= & 
       anyv3d(i,j,k,id_tdiv)  &
      +anyv3d(i,j,k,id_tbef_)*dti_lp*0.5*( &
                                          veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)) &
                                         -veldxdz_v(i,j  ,k,id_velimp)-abs(veldxdz_v(i,j  ,k,id_velimp)) &
                                         )/dxdy_t(i,j)
      enddo ; enddo ; enddo

! Passage boucle j decroissante pour la partie negative de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement negative (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(minval(veldxdz_v(1:imax,jmax+1,1:kmax,id_velimp))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_v en j=3 du present rank est vel_v en j=jmax+1 du voisin aval
      if(minval(veldxdz_v(1:imax,3,1:kmax,id_velimp))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                 status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(veldxdz_v(1:imax,2:jmax+1,1:kmax,id_velimp))>=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                            !>>>  

      do loop_=1,nbdom_jmax !-boucle-mpi->

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(sud)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(nord)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

      call advection_scal_impobc(4,id_tbef_,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=jmax,1,-1 ; do i=1,imax

!      anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       *anyv3d(i,j,k,id_tbef_)                              &
!       -(veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_lp &
!       *anyv3d(i,j+1,k,id_tbef_)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       -(veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_lp &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       anyv3d(i,j,k,id_tbef_)=mask_t(i,j,k)*( & !LSM>

         anyv3d(i,j,k,id_tbef_) &

        +( & !haut>

         (-veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_lp &
        *(anyv3d(i,j+1,k,id_tbef_)-anyv3d(i,j,k,id_tbef_))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +(-veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(4,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

       anyv3d(i,j,k,id_tdiv)= & 
       anyv3d(i,j,k,id_tdiv)  &
      +anyv3d(i,j,k,id_tbef_)*dti_lp*0.5*( &
                                          veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)) &
                                         -veldxdz_v(i,j  ,k,id_velimp)+abs(veldxdz_v(i,j  ,k,id_velimp)) &
                                         )/dxdy_t(i,j)
      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans advection_scal_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_tbef_)
!     call obc_mpi_anyv3d(1,id_tbef_,'zb')

! Tout A la fin (apres advection Oj)
      do k=1,kmax ; do j=1,jmax  ; do i=1,imax  
           tem_t(i,j,k,2)                              &
         =anyv3d(i,j,k,id_tbef_)                    !  &
      enddo ; enddo ; enddo


      end subroutine advection_scal_tem_horimpli

!..............................................................................

      subroutine advection_scal_sal_horimpli
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      integer :: t_=1        &
                ,id_sbef_=3  &  ! salinity "before" identifier
                ,loop_       &
                ,flag_mpi_=0 &
                ,work_status_=0     & !0=rien A signaler, 1=travail demandE       , 2=travail fait
                ,send_status_=0     & !0=rien A signaler, 1=envoi mpi demandE     , 2=envoi mpi rEalisE
                ,recv_status_=0     & !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
                ,status_voisin_upst & ! 0= pas d'echange avec voisin amont, 1= echange
                ,status_voisin_down   ! 0= pas d'echange avec voisin aval , 1= echange




! flag_mpi_=0: je ne fais rien
! flag_mpi_=1: je fais le calcul
! flag_mpi_=2: j'ai fait le calcul, j'envoie la valeur au voisin
! flag_mpi_=3: j'ai tout fait, je m'arrete

      do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+1
!      anyv3d(i,j,k,id_sbef_)=sal_t(i,j,k,before) ! pour schema 100% implicite
       anyv3d(i,j,k,id_sbef_)=sal_t(i,j,k,2)      ! pour schema hybride explicite/implicite
      enddo ; enddo ; enddo

! Schema 100% implicite decommenter les lignes suivantes, sinon le cumul de anyv3d(i,j,k,id_sdiv) se poursuit
!     do k=1,kmax ; do j=1,jmax ; do i=1,imax 
!      anyv3d(i,j,k,id_sdiv)=0. ! reset cumul temperature*divergence(courant) !26-01-17
!     enddo ; enddo ; enddo

!...........................
! Advection axe Oi
!...........................

! Passage boucle i croissante pour la partie positive de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0       ! 0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0       ! 0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0       ! 0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement positive (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(maxval(veldydz_u(1,1:jmax,1:kmax,id_velimp))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_u en i=imax-1 du present rank est vel_u en i=1 du voisin aval
      if(maxval(veldydz_u(imax-1,1:jmax,1:kmax,id_velimp))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(veldydz_u(1:imax,1:jmax,1:kmax,id_velimp))<=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                          !>>>  

      do loop_=1,nbdom_imax !-boucle-mpi->
!     write(10+par%rank,*)'wrs',loop_,work_status_,recv_status_,send_status_

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(est)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(ouest)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

         call advection_scal_impobc(1,id_sbef_,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

!      anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       *anyv3d(i,j,k,id_sbef_)                              &
!       +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_lp &
!       *anyv3d(i-1,j,k,id_sbef_)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_lp &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*sal_t(i,j,k,before)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !LSM>

         anyv3d(i,j,k,id_sbef_) &

        +( & !haut>

         (veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_lp &
        *(anyv3d(i-1,j,k,id_sbef_)-anyv3d(i,j,k,id_sbef_))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*sal_t(i,j,k,before)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(1,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

       anyv3d(i,j,k,id_sdiv)= & 
       anyv3d(i,j,k,id_sdiv)  &
      +anyv3d(i,j,k,id_sbef_)*dti_lp*0.5*( &
                                          veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)) &
                                         -veldydz_u(i  ,j,k,id_velimp)-abs(veldydz_u(i  ,j,k,id_velimp)) &
                                         )/dxdy_t(i,j)
      enddo ; enddo ; enddo

! Passage boucle i decroissante pour la partie negative de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement negative (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(minval(veldydz_u(imax+1,1:jmax,1:kmax,id_velimp))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_u en i=3 du present rank est vel_u en i=imax+1 du voisin aval
      if(minval(veldydz_u(3,1:jmax,1:kmax,id_velimp))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(veldydz_u(2:imax+1,1:jmax,1:kmax,id_velimp))>=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                            !>>>  

      do loop_=1,nbdom_imax !-boucle-mpi->

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(ouest)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(est)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

      call advection_scal_impobc(2,id_sbef_,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=1,jmax ; do i=imax,1,-1

!      anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       *anyv3d(i,j,k,id_sbef_)                              &
!       -(veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)))*dti_lp &
!       *anyv3d(i+1,j,k,id_sbef_)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       -(veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)))*dti_lp &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*sal_t(i,j,k,before)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !LSM>

         anyv3d(i,j,k,id_sbef_) &

        +( & !haut>

         (-veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)))*dti_lp &
        *(anyv3d(i+1,j,k,id_sbef_)-anyv3d(i,j,k,id_sbef_))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +(-veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*sal_t(i,j,k,before)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(2,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

       anyv3d(i,j,k,id_sdiv)= & 
       anyv3d(i,j,k,id_sdiv)  &
      +anyv3d(i,j,k,id_sbef_)*dti_lp*0.5*( &
                                          veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)) &
                                         -veldydz_u(i  ,j,k,id_velimp)+abs(veldydz_u(i  ,j,k,id_velimp)) &
                                         )/dxdy_t(i,j)
      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans advection_scal_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_sbef_)
!     call obc_mpi_anyv3d(1,id_sbef_,'zb')

!...........................
! Advection axe Oj
!...........................

! Passage boucle j croissante pour la partie positive de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement positive (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(maxval(veldxdz_v(1:imax,1,1:kmax,id_velimp))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                 status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_v en j=jmax-1 du present rank est vel_v en j=1 du voisin aval
      if(maxval(veldxdz_v(1:imax,jmax-1,1:kmax,id_velimp))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(veldxdz_v(1:imax,1:jmax,1:kmax,id_velimp))<=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                          !>>>  

      do loop_=1,nbdom_jmax !-boucle-mpi->

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(nord)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(sud)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

         call advection_scal_impobc(3,id_sbef_,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

!      anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       *anyv3d(i,j,k,id_sbef_)                              &
!       +(veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_lp &
!       *anyv3d(i,j-1,k,id_sbef_)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       +(veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_lp &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*sal_t(i,j,k,before)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !LSM>

         anyv3d(i,j,k,id_sbef_) &

        +( & !haut>

         ( veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_lp &
        *(anyv3d(i,j-1,k,id_sbef_)-anyv3d(i,j,k,id_sbef_))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +( veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*sal_t(i,j,k,before)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(3,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1


      enddo                 !-boucle-mpi->

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

       anyv3d(i,j,k,id_sdiv)= & 
       anyv3d(i,j,k,id_sdiv)  &
      +anyv3d(i,j,k,id_sbef_)*dti_lp*0.5*( &
                                          veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)) &
                                         -veldxdz_v(i,j  ,k,id_velimp)-abs(veldxdz_v(i,j  ,k,id_velimp)) &
                                         )/dxdy_t(i,j)
      enddo ; enddo ; enddo

! Passage boucle j decroissante pour la partie negative de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement negative (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(minval(veldxdz_v(1:imax,jmax+1,1:kmax,id_velimp))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_v en j=3 du present rank est vel_v en j=jmax+1 du voisin aval
      if(minval(veldxdz_v(1:imax,3,1:kmax,id_velimp))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                 status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(veldxdz_v(1:imax,2:jmax+1,1:kmax,id_velimp))>=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                            !>>>  

      do loop_=1,nbdom_jmax !-boucle-mpi->

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(sud)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(nord)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

      call advection_scal_impobc(4,id_sbef_,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=jmax,1,-1 ; do i=1,imax

!      anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       *anyv3d(i,j,k,id_sbef_)                              &
!       -(veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_lp &
!       *anyv3d(i,j+1,k,id_sbef_)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       -(veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_lp &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*sal_t(i,j,k,before)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       anyv3d(i,j,k,id_sbef_)=mask_t(i,j,k)*( & !LSM>

         anyv3d(i,j,k,id_sbef_) &

        +( & !haut>

         (-veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_lp &
        *(anyv3d(i,j+1,k,id_sbef_)-anyv3d(i,j,k,id_sbef_))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +(-veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*sal_t(i,j,k,before)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(4,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

       anyv3d(i,j,k,id_sdiv)= & 
       anyv3d(i,j,k,id_sdiv)  &
      +anyv3d(i,j,k,id_sbef_)*dti_lp*0.5*( &
                                          veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)) &
                                         -veldxdz_v(i,j  ,k,id_velimp)+abs(veldxdz_v(i,j  ,k,id_velimp)) &
                                         )/dxdy_t(i,j)
      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans advection_scal_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur anyv3d(i,j,k,id_sbef_)
!     call obc_mpi_anyv3d(1,id_sbef_,'zb')

! Tout A la fin (apres advection Oj)
      do k=1,kmax ; do j=1,jmax  ; do i=1,imax  
           sal_t(i,j,k,2)                              &
         =anyv3d(i,j,k,id_sbef_)                    !  &
      enddo ; enddo ; enddo


      end subroutine advection_scal_sal_horimpli

!..............................................................................

!....................................................
!....................................................
! VERSION POUR tableaux tem_t et sal_t
!31-12-22
      subroutine advection_scal_impobc(case_,idvar_,send_status_,recv_status_)
      use module_principal
      use module_parallele
      implicit none
      integer,parameter :: tagouest_   =5000,  tagest_   =5010   &
                            ,tagsud_   =6000, tagnord_   =6010
      integer,parameter :: tagsudouest_   =15000, tagsudest_   =15010  &
                         ,tagnordouest_   =16000, tagnordest_   =16010 &
                         ,nexchgmax_   =2
!     integer nexchg_   ,warning_   ,var_
      integer nexchg_,loop_,idvar_,case_,send_status_,recv_status_
      integer,dimension(nexchgmax_   ) :: tabreq_
      integer,dimension(mpi_status_size,nexchgmax_   ) :: tstatus_


      if(.not.allocated(adv_send_est)) then !>>>
               allocate(adv_send_est  (2*kmax*jmax)) ; adv_send_est=0.
               allocate(adv_send_ouest(2*kmax*jmax)) ; adv_send_ouest=0.
               allocate(adv_recv_est  (2*kmax*jmax)) ; adv_recv_est=0.
               allocate(adv_recv_ouest(2*kmax*jmax)) ; adv_recv_ouest=0.
      endif                                    !>>>
      if(.not.allocated(adv_send_nord)) then !>>>
               allocate(adv_send_nord  (2*kmax*imax)) ; adv_send_nord=0.
               allocate(adv_send_sud   (2*kmax*imax)) ; adv_send_sud=0.
               allocate(adv_recv_nord  (2*kmax*imax)) ; adv_recv_nord=0.
               allocate(adv_recv_sud   (2*kmax*imax)) ; adv_recv_sud=0.
      endif                                    !>>>


      if(case_==1) then !---echange-Oi-boucle-croissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere Est:
      if (par%tvoisin(est) /= mpi_proc_null) then !pmx>

!      write(10+par%rank,*)'ENVOI send_status_=',send_status_

       if(send_status_==1) then !>>>          
        k0=0
        i=imax-2
        do k=1,kmax ; do j=1,jmax
         k0=k0+1
         adv_send_est(k0)=tem_t(i,j,k,2)
         k0=k0+1
         adv_send_est(k0)=sal_t(i,j,k,2)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_est(1)     & ! envoyé au proc Est
                ,size(adv_send_est(1:k0)) &
                ,mpi_double                    &
                ,par%tvoisin(est)            &
                ,tagest_                     &
                ,par%comm2d                  &
                ,tabreq_   (nexchg_   )      &
                ,ierr)

      endif                                       !pmx>

! Frontiere Ouest:
      if (par%tvoisin(ouest) /= mpi_proc_null) then !m°v°m>

!      write(10+par%rank,*)'RECOI recv_status_=',recv_status_

      if(recv_status_==1) then !>>>          
       k0=2*kmax*jmax
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_ouest(1)    &
               ,size(adv_recv_ouest(1:k0)) &
               ,mpi_double                         &
               ,par%tvoisin(ouest)                           &
               ,tagest_                                      &
               ,par%comm2d                                   &
               ,tabreq_   (nexchg_   )                       &
               ,ierr)

      endif                                         !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(ouest) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       i=0
       k0=0
       do k=1,kmax ; do j=1,jmax
        k0=k0+1
        tem_t(i,j,k,2)=adv_recv_ouest(k0)
        k0=k0+1
        sal_t(i,j,k,2)=adv_recv_ouest(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                         !m°v°m>


      endif             !---echange-Oi-boucle-croissante--->

      if(case_==2) then !---echange-Oi-boucle-decroissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere ouest:
      if (par%tvoisin(ouest) /= mpi_proc_null) then !pmx>

       if(send_status_==1) then !>>>          
        k0=0
        i=3
        do k=1,kmax ; do j=1,jmax
         k0=k0+1
         adv_send_ouest(k0)=tem_t(i,j,k,2)
         k0=k0+1
         adv_send_ouest(k0)=sal_t(i,j,k,2)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_ouest(1)     & ! envoyé au proc Est
                ,size(adv_send_ouest(1:k0)) &
                ,mpi_double                 &
                ,par%tvoisin(ouest)         &
                ,tagouest_                  &
                ,par%comm2d                 &
                ,tabreq_   (nexchg_   )     &
                ,ierr)

      endif                                       !pmx>

! Frontiere Est:
      if (par%tvoisin(est) /= mpi_proc_null) then !m°v°m>

      if(recv_status_==1) then !>>>          
       k0=2*kmax*jmax
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_est(1)     &
               ,size(adv_recv_est(1:k0)) &
               ,mpi_double               &
               ,par%tvoisin(est)         &
               ,tagouest_                &
               ,par%comm2d               &
               ,tabreq_   (nexchg_   )   &
               ,ierr)

      endif                                       !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(est) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       i=imax+1
       k0=0
       do k=1,kmax ; do j=1,jmax
        k0=k0+1
        tem_t(i,j,k,2)=adv_recv_est(k0)
        k0=k0+1
        sal_t(i,j,k,2)=adv_recv_est(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                       !m°v°m>


      endif             !---echange-Oi-boucle-decroissante--->

      if(case_==3) then !---echange-Oj-boucle-croissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere Nord:
      if (par%tvoisin(nord) /= mpi_proc_null) then !pmx>

!      write(10+par%rank,*)'ENVOI send_status_=',send_status_

       if(send_status_==1) then !>>>          
        k0=0
        j=jmax-2
        do k=1,kmax ; do i=1,imax
         k0=k0+1
         adv_send_nord(k0)=tem_t(i,j,k,2)
         k0=k0+1
         adv_send_nord(k0)=sal_t(i,j,k,2)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_nord(1)     & ! envoyé au proc Est
                ,size(adv_send_nord(1:k0)) &
                ,mpi_double                    &
                ,par%tvoisin(nord)            &
                ,tagnord_                     &
                ,par%comm2d                  &
                ,tabreq_   (nexchg_   )      &
                ,ierr)

      endif                                       !pmx>

! Frontiere Sud:
      if (par%tvoisin(sud) /= mpi_proc_null) then !m°v°m>

!      write(10+par%rank,*)'RECOI recv_status_=',recv_status_

      if(recv_status_==1) then !>>>          
       k0=2*kmax*imax
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_sud(1)    &
               ,size(adv_recv_sud(1:k0)) &
               ,mpi_double                         &
               ,par%tvoisin(sud)                           &
               ,tagnord_                                      &
               ,par%comm2d                                   &
               ,tabreq_   (nexchg_   )                       &
               ,ierr)

      endif                                         !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(sud) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       j=0
       k0=0
       do k=1,kmax ; do i=1,imax
        k0=k0+1
        tem_t(i,j,k,2)=adv_recv_sud(k0)
        k0=k0+1
        sal_t(i,j,k,2)=adv_recv_sud(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                         !m°v°m>


      endif             !---echange-Oj-boucle-croissante--->

      if(case_==4) then !---echange-Oj-boucle-decroissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere Sud:
      if (par%tvoisin(sud) /= mpi_proc_null) then !pmx>

       if(send_status_==1) then !>>>          
        k0=0
        j=3
        do k=1,kmax ; do i=1,imax
         k0=k0+1
         adv_send_sud(k0)=tem_t(i,j,k,2)
         k0=k0+1
         adv_send_sud(k0)=sal_t(i,j,k,2)
        enddo       ; enddo
        send_status_=2 ! pour memoriser le fait que l'envoi est fait
       else                     !>>>          
        k0=1
       endif                    !>>>          


      nexchg_   =nexchg_   +1
      call mpi_issend(adv_send_sud(1)     & ! envoyé au proc Est
                ,size(adv_send_sud(1:k0)) &
                ,mpi_double                 &
                ,par%tvoisin(sud)         &
                ,tagsud_                    &
                ,par%comm2d                 &
                ,tabreq_   (nexchg_   )     &
                ,ierr)

      endif                                       !pmx>

! Frontiere Nord:
      if (par%tvoisin(nord) /= mpi_proc_null) then !m°v°m>

      if(recv_status_==1) then !>>>          
       k0=2*kmax*imax
      else                     !>>>
       k0=1 
      endif                    !>>>

      nexchg_   =nexchg_   +1
      call mpi_irecv(adv_recv_nord(1)     &
               ,size(adv_recv_nord(1:k0)) &
               ,mpi_double               &
               ,par%tvoisin(nord)         &
               ,tagsud_                  &
               ,par%comm2d               &
               ,tabreq_   (nexchg_   )   &
               ,ierr)

      endif                                       !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)


      if (par%tvoisin(nord) /= mpi_proc_null) then !m°v°m>
      if(recv_status_==1) then !>>>          
       j=jmax+1
       k0=0
       do k=1,kmax ; do i=1,imax
        k0=k0+1
        tem_t(i,j,k,2)=adv_recv_nord(k0)
        k0=k0+1
        sal_t(i,j,k,2)=adv_recv_nord(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                       !m°v°m>


      endif             !---echange-Oj-boucle-decroissante--->

      end subroutine advection_scal_impobc

!....................................................
!31-12-22
      subroutine advection_scal_flagobc(case_,flag_)
      use module_principal
      use module_parallele
      implicit none
      integer,parameter :: tagouest_   =5000,  tagest_   =5010   &
                            ,tagsud_   =6000, tagnord_   =6010
      integer,parameter :: tagsudouest_   =15000, tagsudest_   =15010  &
                         ,tagnordouest_   =16000, tagnordest_   =16010 &
                         ,nexchgmax_   =2
!     integer nexchg_   ,warning_   ,var_
      integer nexchg_,loop_,case_,flag_
      integer,dimension(nexchgmax_   ) :: tabreq_
      integer,dimension(mpi_status_size,nexchgmax_   ) :: tstatus_

      if(case_==1) then !---echange-Oi-boucle-croissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere Est:
      if (par%tvoisin(est) /= mpi_proc_null) then !pmx>

      nexchg_   =nexchg_   +1
      call mpi_issend(flag_                  & ! envoyé au proc Est
                ,1                           &
                ,mpi_integer                 &
                ,par%tvoisin(est)            &
                ,tagest_                     &
                ,par%comm2d                  &
                ,tabreq_   (nexchg_   )      &
                ,ierr)

      endif                                       !pmx>

! Frontiere Ouest:
      if (par%tvoisin(ouest) /= mpi_proc_null) then !m°v°m>

      nexchg_   =nexchg_   +1
      call mpi_irecv(flag_                 &
               ,1                          &
               ,mpi_integer                &
               ,par%tvoisin(ouest)         &
               ,tagest_                    &
               ,par%comm2d                 &
               ,tabreq_   (nexchg_   )     &
               ,ierr)

      endif                                         !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)

      endif             !---echange-Oi-boucle-croissante--->

      if(case_==2) then !---echange-Oi-boucle-decroissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere ouest:
      if (par%tvoisin(ouest) /= mpi_proc_null) then !pmx>

      nexchg_   =nexchg_   +1
      call mpi_issend(flag_                 & ! envoyé au proc Est
                ,1                          &
                ,mpi_integer                &
                ,par%tvoisin(ouest)         &
                ,tagouest_                  &
                ,par%comm2d                 &
                ,tabreq_   (nexchg_   )     &
                ,ierr)

      endif                                       !pmx>

! Frontiere Est:
      if (par%tvoisin(est) /= mpi_proc_null) then !m°v°m>

      nexchg_   =nexchg_   +1
      call mpi_irecv(flag_               &
               ,1                        &
               ,mpi_integer              &
               ,par%tvoisin(est)         &
               ,tagouest_                &
               ,par%comm2d               &
               ,tabreq_   (nexchg_   )   &
               ,ierr)

      endif                                       !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)

      endif             !---echange-Oi-boucle-decroissante--->

      if(case_==3) then !---echange-Oj-boucle-croissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere Nord:
      if (par%tvoisin(nord) /= mpi_proc_null) then !pmx>

      nexchg_   =nexchg_   +1
      call mpi_issend(flag_                  & ! envoyé au proc Nord
                ,1                           &
                ,mpi_integer                 &
                ,par%tvoisin(nord)            &
                ,tagnord_                     &
                ,par%comm2d                  &
                ,tabreq_   (nexchg_   )      &
                ,ierr)

      endif                                       !pmx>

! Frontiere Sud:
      if (par%tvoisin(sud) /= mpi_proc_null) then !m°v°m>

      nexchg_   =nexchg_   +1
      call mpi_irecv(flag_                 &
               ,1                          &
               ,mpi_integer                &
               ,par%tvoisin(sud)         &
               ,tagnord_                    &
               ,par%comm2d                 &
               ,tabreq_   (nexchg_   )     &
               ,ierr)

      endif                                         !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)

      endif             !---echange-Oj-boucle-croissante--->

      if(case_==4) then !---echange-Oj-boucle-decroissante--->
!.............................................................................................
! Echanger les informations des zones d'echanges:
      nexchg_   =0

! Frontiere Sud:
      if (par%tvoisin(sud) /= mpi_proc_null) then !pmx>

      nexchg_   =nexchg_   +1
      call mpi_issend(flag_                 & ! envoyé au proc Nord
                ,1                          &
                ,mpi_integer                &
                ,par%tvoisin(sud)         &
                ,tagsud_                  &
                ,par%comm2d                 &
                ,tabreq_   (nexchg_   )     &
                ,ierr)

      endif                                       !pmx>

! Frontiere Nord:
      if (par%tvoisin(nord) /= mpi_proc_null) then !m°v°m>

      nexchg_   =nexchg_   +1
      call mpi_irecv(flag_               &
               ,1                        &
               ,mpi_integer              &
               ,par%tvoisin(nord)         &
               ,tagsud_                &
               ,par%comm2d               &
               ,tabreq_   (nexchg_   )   &
               ,ierr)

      endif                                       !m°v°m>

      if(nexchg_   >nexchgmax_   ) &
      stop 'adv_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_  (:,1:nexchg_   ) &
                      ,ierr)

      endif             !---echange-Oj-boucle-decroissante--->

      end subroutine advection_scal_flagobc

!..............................................................
!31-12-22
      subroutine advection_scal_imp_exp
      use module_principal ; use module_parallele
      implicit none

      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1

       veldydz_u(i,j,k,id_velexp)=max(min(veldydz_u(i,j,k,id_veltot) &

        ,+looplimit_hor*min(     anyv3d(i  ,j,k,id_dxdydz)             &
                         /max(wetmask_t(i  ,j),1.e-10)                 &
                                ,anyv3d(i-1,j,k,id_dxdydz)             &
                         /max(wetmask_t(i-1,j),1.e-10)                 &
                           )/dti_lp)                                   &
        ,-looplimit_hor*min(     anyv3d(i  ,j,k,id_dxdydz)             &
                         /max(wetmask_t(i  ,j),1.e-10)                 &
                                ,anyv3d(i-1,j,k,id_dxdydz)             &
                         /max(wetmask_t(i-1,j),1.e-10)                 &
                           )/dti_lp)  


       veldydz_u(i,j,k,id_velimp)=veldydz_u(i,j,k,id_veltot)-veldydz_u(i,j,k,id_velexp)

!      if(veldydz_u(i,j,k,id_velimp)/=0.)write(10+par%rank,*)i,j,k,real(veldydz_u(i,j,k,id_velimp)),' veldydz_u'

      enddo ; enddo ; enddo

!...

      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax

       veldxdz_v(i,j,k,id_velexp)=max(min(veldxdz_v(i,j,k,id_veltot) &

        ,+looplimit_hor*min(     anyv3d(i,j  ,k,id_dxdydz)             &
                         /max(wetmask_t(i,j  ),1.e-10)                 &
                                ,anyv3d(i,j-1,k,id_dxdydz)             &
                         /max(wetmask_t(i,j-1),1.e-10)                 &
                           )/dti_lp)    &
        ,-looplimit_hor*min(     anyv3d(i,j  ,k,id_dxdydz)             &
                         /max(wetmask_t(i,j  ),1.e-10)                 &
                                ,anyv3d(i,j-1,k,id_dxdydz)             &
                         /max(wetmask_t(i,j-1),1.e-10)                 &
                           )/dti_lp)  


       veldxdz_v(i,j,k,id_velimp)=veldxdz_v(i,j,k,id_veltot)-veldxdz_v(i,j,k,id_velexp)

!      if(veldxdz_v(i,j,k,id_velimp)/=0.)write(10+par%rank,*)i,j,k,real(veldxdz_v(i,j,k,id_velimp)),' veldxdz_v'
      enddo ; enddo ; enddo

      sch_imp_ts_u_loc=0
      if(maxval(abs(veldydz_u(:,:,:,id_velimp)))>0.)sch_imp_ts_u_loc=1
      call mpi_allreduce(sch_imp_ts_u_loc,sch_imp_ts_u_glb,1,mpi_integer,mpi_max,par%comm2d,ierr)

      sch_imp_ts_v_loc=0
      if(maxval(abs(veldxdz_v(:,:,:,id_velimp)))>0.)sch_imp_ts_v_loc=1
      call mpi_allreduce(sch_imp_ts_v_loc,sch_imp_ts_v_glb,1,mpi_integer,mpi_max,par%comm2d,ierr)

!     if(sch_imp_ts_u_loc/=0.or.sch_imp_ts_v_loc/=0) then
!      write(10+par%rank,*)'sch_imp_ts_u_loc',sch_imp_ts_u_loc
!      write(10+par%rank,*)'sch_imp_ts_u_glb',sch_imp_ts_u_glb
!      write(10+par%rank,*)'sch_imp_ts_v_loc',sch_imp_ts_v_loc
!      write(10+par%rank,*)'sch_imp_ts_v_glb',sch_imp_ts_v_glb
!     endif

      if(mod(iteration3d,100)==0) then !-statistiques->

       k0=0 ; k1=0
       do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
        k0=k0+mask_u(i,j,kmax)
        if(veldydz_u(i,j,k,id_velimp)/=0.)k1=k1+1
       enddo       ; enddo       ; enddo
       call mpi_allreduce(k0,k2,1,mpi_integer,mpi_sum,par%comm2d,ierr)
       call mpi_allreduce(k1,k3,1,mpi_integer,mpi_sum,par%comm2d,ierr)
       x1_r4=real(k3)/real(k2)

       k0=0 ; k1=0
       do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
        k0=k0+mask_v(i,j,kmax)
        if(veldxdz_v(i,j,k,id_velimp)/=0.)k1=k1+1
       enddo       ; enddo       ; enddo
       call mpi_allreduce(k0,k2,1,mpi_integer,mpi_sum,par%comm2d,ierr)
       call mpi_allreduce(k1,k3,1,mpi_integer,mpi_sum,par%comm2d,ierr)
       x2_r4=real(k3)/real(k2)

       if(par%rank==0) then !-rk0->
        open(unit=3,file=trim(tmpdirname)//'dti_stat_ts' &
                   ,position='append')
          write(3,*)real(elapsedtime_now/86400.) &
                   ,real(x1_r4*100.) & ! % de points concernes par adv implicite Oi
                   ,real(x2_r4*100.)   ! % de points concernes par adv implicite Oj
       close(3)
       endif                !-rk0->

      endif                            !-statistiques->

      end subroutine advection_scal_imp_exp

!..............................................................................
!31-12-22
      subroutine advection_scal_ts_horimpli
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      integer :: t_=1        &
                ,loop_       &
                ,flag_mpi_=0 &
                ,work_status_=0      & !0=rien A signaler, 1=travail demandE       , 2=travail fait
                ,send_status_=0      & !0=rien A signaler, 1=envoi mpi demandE     , 2=envoi mpi rEalisE
                ,recv_status_=0      & !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
                ,status_voisin_upst  & ! 0= pas d'echange avec voisin amont, 1= echange
                ,status_voisin_down  & ! 0= pas d'echange avec voisin aval , 1= echange
                ,work_status_mpi_min_  !04-01-23

      if(iperiodicboundary==.true.)stop 'conditions periodiques A faire'
      if(jperiodicboundary==.true.)stop 'conditions periodiques A faire'

! flag_mpi_=0: je ne fais rien
! flag_mpi_=1: je fais le calcul
! flag_mpi_=2: j'ai fait le calcul, j'envoie la valeur au voisin
! flag_mpi_=3: j'ai tout fait, je m'arrete

! Schema 100% implicite decommenter les lignes suivantes, sinon le cumul de anyv3d(i,j,k,id_tdiv) se poursuit
!     do k=1,kmax ; do j=1,jmax ; do i=1,imax 
!      anyv3d(i,j,k,id_tdiv)=0. ! reset cumul temperature*divergence(courant) !26-01-17
!     enddo ; enddo ; enddo

!...........................
! Advection axe Oi
!...........................

! Passage boucle i croissante pour la partie positive de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0       ! 0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0       ! 0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0       ! 0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement positive (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(maxval(veldydz_u(1,1:jmax,1:kmax,id_velimp))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_u en i=imax-1 du present rank est vel_u en i=1 du voisin aval
      if(maxval(veldydz_u(imax-1,1:jmax,1:kmax,id_velimp))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(veldydz_u(1:imax,1:jmax,1:kmax,id_velimp))<=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                          !>>>  

      do loop_=1,nbdom_imax !-boucle-mpi->

       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 10 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(est)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(ouest)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

         call advection_scal_impobc(1,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

!      tem_t(i,j,k,2)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       *tem_t(i,j,k,2)                              &
!       +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_lp &
!       *tem_t(i-1,j,k,2)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_lp &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       tem_t(i,j,k,2)=mask_t(i,j,k)*( & !LSM>

         tem_t(i,j,k,2) &

        +( & !haut>

         (veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_lp &
        *(tem_t(i-1,j,k,2)-tem_t(i,j,k,2))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    

       sal_t(i,j,k,2)=mask_t(i,j,k)*( & !LSM>

         sal_t(i,j,k,2) &

        +( & !haut>

         (veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_lp &
        *(sal_t(i-1,j,k,2)-sal_t(i,j,k,2))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +(veldydz_u(i,j,k,id_velimp)+abs(veldydz_u(i,j,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*sal_t(i,j,k,before)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(1,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

!      call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
!      if(work_status_mpi_min_==2)goto 10 ! Des que tous les travaux sont accomplis, faire l'etape suivante

      enddo                 !-boucle-mpi->
   10 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

       anyv3d(i,j,k,id_tdiv)= & 
       anyv3d(i,j,k,id_tdiv)  &
      +tem_t(i,j,k,2)*dti_lp*0.5*( &
                                          veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)) &
                                         -veldydz_u(i  ,j,k,id_velimp)-abs(veldydz_u(i  ,j,k,id_velimp)) &
                                         )/dxdy_t(i,j)

       anyv3d(i,j,k,id_sdiv)= & 
       anyv3d(i,j,k,id_sdiv)  &
      +sal_t(i,j,k,2)*dti_lp*0.5*( &
                                          veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)) &
                                         -veldydz_u(i  ,j,k,id_velimp)-abs(veldydz_u(i  ,j,k,id_velimp)) &
                                         )/dxdy_t(i,j)

      enddo ; enddo ; enddo

! Passage boucle i decroissante pour la partie negative de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement negative (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(minval(veldydz_u(imax+1,1:jmax,1:kmax,id_velimp))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_u en i=3 du present rank est vel_u en i=imax+1 du voisin aval
      if(minval(veldydz_u(3,1:jmax,1:kmax,id_velimp))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(veldydz_u(2:imax+1,1:jmax,1:kmax,id_velimp))>=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                            !>>>  

      do loop_=1,nbdom_imax !-boucle-mpi->

       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 11 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(ouest)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(est)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

      call advection_scal_impobc(2,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=1,jmax ; do i=imax,1,-1

!      tem_t(i,j,k,2)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       *tem_t(i,j,k,2)                              &
!       -(veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)))*dti_lp &
!       *tem_t(i+1,j,k,2)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       -(veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)))*dti_lp &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       tem_t(i,j,k,2)=mask_t(i,j,k)*( & !LSM>

         tem_t(i,j,k,2) &

        +( & !haut>

         (-veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)))*dti_lp &
        *(tem_t(i+1,j,k,2)-tem_t(i,j,k,2))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +(-veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    

       sal_t(i,j,k,2)=mask_t(i,j,k)*( & !LSM>

         sal_t(i,j,k,2) &

        +( & !haut>

         (-veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)))*dti_lp &
        *(sal_t(i+1,j,k,2)-sal_t(i,j,k,2))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +(-veldydz_u(i+1,j,k,id_velimp)+abs(veldydz_u(i+1,j,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*sal_t(i,j,k,before)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(2,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

!      call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
!      if(work_status_mpi_min_==2)goto 11 ! Des que tous les travaux sont accomplis, faire l'etape suivante
      enddo                 !-boucle-mpi->
   11 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

       anyv3d(i,j,k,id_tdiv)= & 
       anyv3d(i,j,k,id_tdiv)  &
      +tem_t(i,j,k,2)*dti_lp*0.5*( &
                                          veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)) &
                                         -veldydz_u(i  ,j,k,id_velimp)+abs(veldydz_u(i  ,j,k,id_velimp)) &
                                         )/dxdy_t(i,j)

       anyv3d(i,j,k,id_sdiv)= & 
       anyv3d(i,j,k,id_sdiv)  &
      +sal_t(i,j,k,2)*dti_lp*0.5*( &
                                          veldydz_u(i+1,j,k,id_velimp)-abs(veldydz_u(i+1,j,k,id_velimp)) &
                                         -veldydz_u(i  ,j,k,id_velimp)+abs(veldydz_u(i  ,j,k,id_velimp)) &
                                         )/dxdy_t(i,j)
      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans advection_scal_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur tem_t(i,j,k,2)
!     call obc_mpi_tem_t(1,2,'zb')

!...........................
! Advection axe Oj
!...........................

! Passage boucle j croissante pour la partie positive de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement positive (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(maxval(veldxdz_v(1:imax,1,1:kmax,id_velimp))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                 status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_v en j=jmax-1 du present rank est vel_v en j=1 du voisin aval
      if(maxval(veldxdz_v(1:imax,jmax-1,1:kmax,id_velimp))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(veldxdz_v(1:imax,1:jmax,1:kmax,id_velimp))<=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                          !>>>  

      do loop_=1,nbdom_jmax !-boucle-mpi->

       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 20 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(nord)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(sud)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

         call advection_scal_impobc(3,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

!      tem_t(i,j,k,2)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       *tem_t(i,j,k,2)                              &
!       +(veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_lp &
!       *tem_t(i,j-1,k,2)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       +(veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_lp &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       tem_t(i,j,k,2)=mask_t(i,j,k)*( & !LSM>

         tem_t(i,j,k,2) &

        +( & !haut>

         ( veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_lp &
        *(tem_t(i,j-1,k,2)-tem_t(i,j,k,2))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +( veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    

       sal_t(i,j,k,2)=mask_t(i,j,k)*( & !LSM>

         sal_t(i,j,k,2) &

        +( & !haut>

         ( veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_lp &
        *(sal_t(i,j-1,k,2)-sal_t(i,j,k,2))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +( veldxdz_v(i,j,k,id_velimp)+abs(veldxdz_v(i,j,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*sal_t(i,j,k,before)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(3,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1


!      call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
!      if(work_status_mpi_min_==2)goto 20 ! Des que tous les travaux sont accomplis, faire l'etape suivante
      enddo                 !-boucle-mpi->
   20 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

       anyv3d(i,j,k,id_tdiv)= & 
       anyv3d(i,j,k,id_tdiv)  &
      +tem_t(i,j,k,2)*dti_lp*0.5*( &
                                          veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)) &
                                         -veldxdz_v(i,j  ,k,id_velimp)-abs(veldxdz_v(i,j  ,k,id_velimp)) &
                                         )/dxdy_t(i,j)

       anyv3d(i,j,k,id_sdiv)= & 
       anyv3d(i,j,k,id_sdiv)  &
      +sal_t(i,j,k,2)*dti_lp*0.5*( &
                                          veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)) &
                                         -veldxdz_v(i,j  ,k,id_velimp)-abs(veldxdz_v(i,j  ,k,id_velimp)) &
                                         )/dxdy_t(i,j)
      enddo ; enddo ; enddo

! Passage boucle j decroissante pour la partie negative de la vitesse
! reset:
      status_voisin_upst=1 ! Par defaut on suppose echange avec voisin amont
      status_voisin_down=1 ! Par defaut on suppose echange avec voisin aval
      work_status_=0   !0=rien A signaler, 1=travail demandE,   2=travail fait
      send_status_=0   !0=rien A signaler, 1=envoi mpi demandE, 2=envoi mpi rEalisE
      recv_status_=0   !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
      flag_mpi_=0

!...............................
! Connaitre le status des voisins amont et aval (echange ou pas echange)
! 1- Le proc voisin amont (selon sign du courant donc) n'existe pas
! 2- Aucune des vitesses de la premiere rangee est strictement negative (le proc ne depend donc pas du proc voisin amont)
! Voisin amont:
      if(minval(veldxdz_v(1:imax,jmax+1,1:kmax,id_velimp))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: vel_v en j=3 du present rank est vel_v en j=jmax+1 du voisin aval
      if(minval(veldxdz_v(1:imax,3,1:kmax,id_velimp))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                 status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(veldxdz_v(1:imax,2:jmax+1,1:kmax,id_velimp))>=0.) then !>>> ! raison 1
        work_status_=2 ; send_status_=2 ; recv_status_=2
      endif                                                            !>>>  

      do loop_=1,nbdom_jmax !-boucle-mpi->

       call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
       if(work_status_mpi_min_==2)goto 21 ! Des que tous les travaux sont accomplis, faire l'etape suivante

       flag_mpi_=0
!      if(work_status_==2.and.send_status_==0.and.par%tvoisin(sud)/=mpi_proc_null)send_status_=1
!      if(work_status_==1.and.par%tvoisin(nord)/=mpi_proc_null)recv_status_=1
       if(work_status_==2.and.send_status_==0.and.status_voisin_down==1)send_status_=1
       if(work_status_==1.and.                    status_voisin_upst==1)recv_status_=1

      call advection_scal_impobc(4,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=jmax,1,-1 ; do i=1,imax

!      tem_t(i,j,k,2)=mask_t(i,j,k)*( & !haut>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       *tem_t(i,j,k,2)                              &
!       -(veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_lp &
!       *tem_t(i,j+1,k,2)                            &
!                                           ) & !haut>
!       /( & !bas>
!        (2.*dxdy_t(i,j)*dz_t(i,j,k,before))                 &
!       -(veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_lp &
!        ) & !bas>
!          +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    
! Noter que le calcul suivant (equivalent) est prEfErE pour garantir que
! le champ soit STRICTEMENT inchangE si la vitesse est nulle (pour la
! continuite mpi)
       tem_t(i,j,k,2)=mask_t(i,j,k)*( & !LSM>

         tem_t(i,j,k,2) &

        +( & !haut>

         (-veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_lp &
        *(tem_t(i,j+1,k,2)-tem_t(i,j,k,2))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +(-veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*tem_t(i,j,k,before)    


       sal_t(i,j,k,2)=mask_t(i,j,k)*( & !LSM>

         sal_t(i,j,k,2) &

        +( & !haut>

         (-veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_lp &
        *(sal_t(i,j+1,k,2)-sal_t(i,j,k,2))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_t(i,j)*dz_t(i,j,k,before)) &

        +(-veldxdz_v(i,j+1,k,id_velimp)+abs(veldxdz_v(i,j+1,k,id_velimp)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_t(i,j,k))*sal_t(i,j,k,before)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(4,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

!      call mpi_allreduce(work_status_,work_status_mpi_min_,1,mpi_integer,mpi_min,par%comm2d ,ierr)
!      if(work_status_mpi_min_==2)goto 21 ! Des que tous les travaux sont accomplis, faire l'etape suivante
      enddo                 !-boucle-mpi->
   21 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=1,jmax ; do i=1,imax

       anyv3d(i,j,k,id_tdiv)= & 
       anyv3d(i,j,k,id_tdiv)  &
      +tem_t(i,j,k,2)*dti_lp*0.5*( &
                                          veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)) &
                                         -veldxdz_v(i,j  ,k,id_velimp)+abs(veldxdz_v(i,j  ,k,id_velimp)) &
                                         )/dxdy_t(i,j)

       anyv3d(i,j,k,id_sdiv)= & 
       anyv3d(i,j,k,id_sdiv)  &
      +sal_t(i,j,k,2)*dti_lp*0.5*( &
                                          veldxdz_v(i,j+1,k,id_velimp)-abs(veldxdz_v(i,j+1,k,id_velimp)) &
                                         -veldxdz_v(i,j  ,k,id_velimp)+abs(veldxdz_v(i,j  ,k,id_velimp)) &
                                         )/dxdy_t(i,j)
      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans advection_scal_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur tem_t(i,j,k,2)
!     call obc_mpi_tem_t(1,2,'zb')

      end subroutine advection_scal_ts_horimpli

!..............................................................................
