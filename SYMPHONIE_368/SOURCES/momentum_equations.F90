      subroutine momentum_equations(texte_   )
!______________________________________________________________________
! SYMPHONIE ocean model
! release 367 - last update: 21-02-23
!______________________________________________________________________

      use module_principal ; use module_parallele ; use module_mangrove
      implicit none
      character texte_*20
      integer t_
#ifdef synopsis
       subroutinetitle='momentum_equations'
       subroutinedescription='Computes the 3D momentum equations'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!...............................................................................
! Version date      Description des modifs
!         19/07/01: passage à la coordonnée sigma généralisée
!         28/07/01: changement de nom pour le courant moyen (VMEAN)
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         26/08/02: VHZOBC remplacé par VELOBC
!         10/02/03: les vitesses doivent être multipliées par leur mask
!                   + divers arrangements pour les besoins de la grille hybride
!                   X1 X2 X9 supprimés
!         02/11/03: bilan energetique des termes turbulents
!         30/03/06: appel (au choix) à TRIDIAGSYST_DBP
!         22/07/06: reecriture totale: tous les processus sont regroupes
!         07/08/06: Modes externes parallelles
!         07/11/06: Debug OBC (boucle 1,imax & 1,jmax) devient 1:imax+1,1:jmax+1
!         19/04/07: Passage à coordonnee curviligne
!         18/02/08: Avant seule l'advection horizontale de u et v etait hybride
!                   centré upwind (qui est la même chose que centré diffusif).
!                   La methode est generalisée à l'advection verticale. Le terme
!                   diffusif vertical est au temps t+1 alors que pour l'horizontale
!                   il est à t-1
!         15/03/08: correction eponge velobc
!         02/04/08: Pour eviter que la modif sur difver_x et difver_y (qui ajoute
!                   l'effet diffusif de l'advection) n'ait un effet sur le calcul
!                   de la turbulence, on utilise un tableau temporaire pour calculer
!                   le melange vertical.
!         03/04/08: retour en arriere vis de la diffusion ajoutee à l'advection
!                   verticale depuis l'adoption de l'option 4 pour la diffusion
!                   horizontale qui n'est plus du même type que celle ajoutée le
!                   18/02/08. Cette diffusion est supprimée, seule subsiste donc la
!                   diffusion verticale de la fermeture turbulente.
!         11/11/08: Division par facteur d'echelle verticale au temps t+1 (et non
!                   plus au temps t)
!         28/03/09: WET/DRY option
!         27/04/09: WET/DRY option: une option plus robuste suite au retour
!                   d'experience de Gabriel Jorda
!         06/05/09: Retour d'experience de Gabriel Jorda. min(h_z(i),h_z(i-1) plutôt
!                   que h_x
! 2009.3  30-09-09  difver_x et difver_y supprimés et remplacés par des 1/2 sommes
!                   de kz_z
!                   Suppresion de omega_x et omega_y
!         04-10-09  uv_upd.F devient momentum_equations.F
! 2010.3  08-01-10  finalisation masquage bancs decouvrant
!         14-01-10  suppression constante inutilisée
! 2010.7  15-02-10  La variable d'etat devient vel_u et vel_v
!         02-03-10  Ajout forces de Stokes. Moyenne verticale du terme de rotation
!                   pour couplage 3d2d desormais inclus dans dif3d2d
! 2010.8  10-03-10  4 levels time filter
!         12-03-10  mixer le 3L et le 4L avec tfc1 et tfc2
!         02-04-10  parametrisation ardhuin suite: frottement sur le fond
! 2010.9  06-06-10  tridiagonalsolver renommé tridiagonalsolver
! 2010.11 30-06-10  suppresion de la soustraction de velavr dans le
!                   terme de tension de fond
! 2010.16 12-01-11  Filtre temporel ordre élevé
! 2010.23 18-05-11  Rappel OGCM implicite
! 2010.25 08-02-12  kz_w renomme km_w
!         30-03-12  Distribution vertical du stress de surface sur une epaisseur
!                   dependante de la hauteur des vagues
! S.26    04-09-12  gradient de pression donné ssh_int(1)
!         27-01-13  advection qdm o4
!         10-06-13  debug terme de rappel (*dz)
!         23-06-13  ajout mangrove
!         13-08-14  Amelioration wetdrying scheme: au lieu d'annuler l'integralite
!                   du courant mode interne (qui entraine la chute de la tke et 
!                   donc du kz) wetmask est applique en amont, dans l'equation 
!                   des moments, aux seuls termes sources. Les termes turbulents
!                   restent actifs, ce qui permet le developpement d'un profil
!                   cisaille maintenant un haut niveau de turlence, comme attendu.
!         05-11-14  retour en arriere aout 2014
!         12-11-14  diffusivitee donnee par gradient de vitesse (retour en arriere)
!         01-12-14  remplacer dsig par dz/hz
!         10-07-15  rhpzavr est une moyenne verticale ponderee de rhp dependant de x,y
!         28-10-15  modif sur diffusion upwind
!         05-01-16  Advection pas de temps sEparEs
!         19-01-16  Amenagements pour cas test baroclinic jet:
!                    1- Possibilite d'une friction de fond lineaire
!                    2- Ajout momentum_zonal_restoring 
!         20-01-16  baroclinic jet: la moyenne zonale esr calculEe sur velobc-vel
!         28-01-16  nouveau time-stepping hybride LF-FB
!         09-02-16  filtre FB et LP: nouveaux noms pour les coef tfilterfb et tfilterlf
!         14-02-16  - Suppression de lignes inutiles
!                   - Pas d'advection en 1DV
!         10-03-16  Seuillage argument de la fonction exp
!         06-08-16  advections selon Oi et Oj separees
!         15-08-16  dz(:,:,:,0) (plutot que dz(:,:,:,2)) semble plus coherent avec le dz de l'advection partielle
!         19-08-16  nbre de courant divisE par min(dz(k),dz(k-1)
!         22-08-16  dz(:,:,:,0) plutot que dz(:,:,:,2) suite....
!         28-09-16  distribution verticale pour wstresb
!         20-11-16  Mises A jour du coef "smago like" de la viscosite bi-harmonique
!         25-11-16  Suite du point precedent (methode 3 corrigee) et C.L. continentales pour bi-harmonique
!         19-12-16  distribution lineraire de la QDM apportee par le vent et les vagues en surface
!         23-12-16  momentum_input_depth est la profondeur de penetration
!                   des inputs de QDM, defini dans notebook_visco
!         16-01-17  vitesses en facteur de du/dx et dv/dy au temps before
!         01-02-17  advection verticale pour u et v
!         06-02-17  cumul du pivot doit etre avant ligne advection partielle
!         15-02-17  advection si flag_1dv=0
!         24-02-17  Mise A jour de la distribution verticale de l'injection de QDM
!         06-04-17  Re_introduction cst_adv_vel dans bilaplacien
!         08-04-17  s-z coordinate
!         25-04-17  s-z coordinate suite
!         03-05-17  pas de melange si kmerge=1
!         14-05-17  biharm_c2*0.5 biharm_c1 pour un coef viscositE basEe sur une
!                   proportion du gradient de vitesse et de la vitesse
!         15-05-17  arret d'urgence si loopmax_ out of range
!         07-10-17  fusion avant restitution du terme de divergence
!         17-10-17  Distribution homogene du stress de fond sur les couches fusionnees:
!                   Chaque couche recoit une fraction (dz/(ztop-zfond) du stress total:
!         24-10-17  suite du point precedent
!         10-11-17  dztimea et dztimeb permettent d'adapter la chronologie des variables
!                   au cas sans time-splitting
!         07-12-17  ajout kundermin
!         09-06-18  viscosite biharmonique hybride: constante & gradient de vitesse
!         29-08-18  if(flag_timesplitting==1) then !m°v°m> 
!         03-09-18  Pas d'advection si flag_timesplitting_adve_uv=1 et iteration=0 
!         18-01-19  ajout subroutine momentum_replay_ncumax
! v246    21-02-19  Modif du dimensionnement du coef de melange dans la  couche fusionnEe
! v247    01-03-19 !call vertmix_merged_levels_uv(2) !01-03-19
! v253    01-05-19  Evitement de la division par dz
! v258    19-07-19  ne pas employer des dz dans le masque pour les flux diffusif
! v269    04-12-19  On ne moyenne plus kmin et kmerged (seule l'advection isolement
!                   est moyennee). Le stress de fond est appliquE seulement entre
!                   -h+dz(kmerged) autrement dit la couche kmerged peut ne pas etre
!                   freinEe si la couche kmin est tres epaisse
! v270    12-12-19  suite point precedent (debug sur calcul velbot)
! v280    28-04-20  temps du facteur d'echelle dz devant presgrad
!                   temps de l'advection verticale cetree
! v287    17-07-20  utiliser tmpdirname
!         14-08-20  obc_int_anyv3d renomme obc_mpi_anyv3d
! v292    19-11-20  subroutine momemtum_logprofile !19-11-20
! v310    17-10-21  possibilite (si on enleve les commentaires) d'utiliser upwindriver
!         20-10-21  ajout subroutines momentum_exp_adv_u_v2 , subroutine momentum_exp_adv_v_v2 
!         04-11-21  VQS: - remplacer la methode de la couche moyenne par un limiteur de flux
!                        - call vertmix_merged_levels_uv(2) de nouveau en service mais pour 
!                          un melange partiel proportionnel au chevauchement de la couche
!                          virtuelle de fond et la couche kmerged
!        08-11-21   faire dti_u et dti_v pour nouveau schema d'advection
! v321   17-01-22   nouvelles subroutines momentum_exp_adv_u, momentum_exp_adv_v
!                   ,momentum_uv_minmax(t_) 
! v335   05-03-22   Advection upwind en zones intertidales
!        07-03-22   A la surface vitesse verticale 100% implicite et vel(kmax+1)=0
! v344   05-04-22   nu doit tendre vers 1 quand wetmask tend vers 0 (pour basculer sur advection upwind).
! v345   06-04-22   *wetmask
! v348   11-06-22   - condition no slip ou free slip
!                   - macro rugosite
! v349   28-06-22   retablissement du point 17-01-22 c.a.d. viscosite donnEe par vitesse advectante
!        29-06-22   if(coastal_viscosity>0.) call dragcoef_no_slip_condition_3d 
! v361   31-12-22   advection horizontale implicite
!        01-01-23   subroutine momentum_bihar_limiter pour limiter la diffusion biharmonique par un critere CFL
! v362   04-01-23   work_status_mpi_min_ 
! v363   15-01-23   debug vitesse implicit
! v367   21-02-23   algo fric_t favorisant la continuite (ecrit mais commente car pas vraiment concluant)
!...............................................................................
!    _________                    .__                  .__             !m°v°m 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !

      if(texte_   (1:19)=='after external mode') then !AAAAAAAAAAAAAAAAAA>

      call momentum_bihar_limiter !01-01-23

      if(flag_adve3d==0) then !-prime->
         call momentum_exp_adv_u_prime !17-01-22
         call momentum_exp_adv_v_prime !17-01-22
      endif                   !-prime->
      if(flag_adve3d==1) stop 'UPWIND 3D A FAIRE' !29-06-22
      if(flag_adve3d==2) then !-modul->
         call momentum_exp_adv_u_modul !17-01-22
         call momentum_exp_adv_v_modul !17-01-22
      endif                   !-modul->

!-----------------------------------------------------------------------------------------
! Vertical distribution of the surface stress:   !30-03-12 !24-02-17
! Echelle de distance:
      if(iwve==1) then !with-waves->
       do j=1,jmax ; do i=1,imax
        xy_t(i,j,1)=1./(0.64*hsw_wave_t(i,j,1)) ! indexee sur hauteur vagues deferlante
!       xy_t(i,j,1)=1./(      hs_wave_t(i,j,1)) ! indexee sur hauteur significative des vagues
!       xy_t(i,j,1)=1./momentum_input_depth     ! empirique
       enddo       ; enddo
      else             !no-waves->
       do j=1,jmax ; do i=1,imax
        xy_t(i,j,1)=1./momentum_input_depth !23-12-16
       enddo       ; enddo
      endif            !no-waves->
! Forme du profil
! Profil "gradient d'une fonction lineaire" (avec normalisation de l'integrale)
! Debut:
!     Forme du profil: max(1+z/Hs,0)
      do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,0)=max(0.,1.+(depth_w(i,j,k)-depth_w(i,j,kmax+1))*xy_t(i,j,1))   ! Linear profil
      enddo         ; enddo       ; enddo
!     Normalisation: 1/(Fsurf-Fbot)
      x0=0.5/rho ! anticipation de la demi-somme "grille C" + division par rho
      do j=1,jmax ; do i=1,imax
       xy_t(i,j,0)=x0/(anyv3d(i,j,kmax+1,0)-anyv3d(i,j,1,0))
      enddo       ; enddo
!     Differentiation verticale: dF/dz 
      do k=1,kmax ; do j=1,jmax ; do i=1,imax
       anyv3d(i,j,k,0)=xy_t(i,j,0)*(anyv3d(i,j,k+1,0)-anyv3d(i,j,k,0))
      enddo       ; enddo       ; enddo
! Fin.
! Les 3 profils de Uchiyama et al 2010, Eq 53:
! Debut:
!     do j=1,jmax ; do i=1,imax
!      xy_t(i,j,0)=small2 ! Pour la normalisation de l'integrale
!     enddo       ; enddo
!     do k=1,kmax ; do j=1,jmax ; do i=1,imax
!      anyv3d(i,j,k,0)= &
!     1.-tanh(xy_t(i,j,1)*(depth_w(i,j,kmax+1)-depth_w(i,j,k)))**4  ! Type I   profil (Uchiyama et al 2010, Eq 53)
!     1.-tanh(xy_t(i,j,1)*(depth_w(i,j,kmax+1)-depth_w(i,j,k)))**2  ! Type II  profil (Uchiyama et al 2010, Eq 53)
!        cosh(xy_t(i,j,1)*(depth_w(i,j,k)+h_w(i,j)))                ! Type III profil (Uchiyama et al 2010, Eq 53)
!      xy_t(i,j,0)=xy_t(i,j,0)+anyv3d(i,j,k,0) ! Pour la normalisation de l'integrale
!     enddo       ; enddo       ; enddo
!     x0=0.5/rho ! anticipation de la demi-somme "grille C" + division par rho
!     do j=1,jmax ; do i=1,imax
!      xy_t(i,j,0)=x0/xy_t(i,j,0)
!     enddo       ; enddo
!     do k=1,kmax ; do j=1,jmax ; do i=1,imax
!      anyv3d(i,j,k,0)=xy_t(i,j,0)*anyv3d(i,j,k,0)
!     enddo       ; enddo       ; enddo
! Fin.
      check0=anyv3d(imax/2,jmax/2,kmax/2,0)
! Verification de la normalisation de l'integrale 
!     i=imax/2 ; j=jmax/2
!     sum2=0.
!     do k=1,kmax
!       sum2=sum2+anyv3d(i,j,k,0)/x0
!       write(66,*)depth_t(i,j,k),anyv3d(i,j,k,0)/x0,sum2 ! sum2=1?
!     enddo
!-----------------------------------------------------------------------------------------

! Vertical distribution of the wave bottom injection of momentum
      if(iwve==1)call momentum_wstresbprofile ! bloque anyv3d(:,:,:,1)

!...............................
! PRECALCULS:
! 1 Divers 2D phase 1:
      if(flag_linearfric==1) then !fricfricfric>
         fric_t=coef_linearfric   !19-01-16
      else                        !fricfricfric>
       do j=1,jmax ; do i=1,imax
!       fric_t(i,j)=cdb_t(i,j)*sqrt(0.5*(                 &
!                  vel_u(i+1,j  ,kmin_u(i+1,j  ),1)**2    &
!                 +vel_u(i  ,j  ,kmin_u(i  ,j  ),1)**2    &
!                 +vel_v(i  ,j+1,kmin_v(i  ,j+1),1)**2    &
!                 +vel_v(i  ,j  ,kmin_v(i  ,j  ),1)**2 ))

! Cet algo favorise les extremas:
        fric_t(i,j)=cdb_t(i,j)*sqrt(0.5*(  &
                   velbot_u(i+1,j  )**2    &
                  +velbot_u(i  ,j  )**2    &
                  +velbot_v(i  ,j+1)**2    &
                  +velbot_v(i  ,j  )**2 ))

! Cet algo favorise la continuite: !21-02-23
!       fric_t(i,j)=cdb_t(i,j)*sqrt(  &
!                       (0.5*(velbot_u(i+1,j  )+velbot_u(i  ,j  )))**2 &
!                      +(0.5*(velbot_v(i  ,j+1)+velbot_v(i  ,j  )))**2 &
!                                  )
       enddo ; enddo
      endif                       !fricfricfric>

! ETAPE 2: precalcul des flux advectifs horizontaux et de coriolis:

!$ u component of the velocity:

      id_gradssh=1
      if(flag_timesplitting==1) then !m°v°m> !29-08-18
       if(iteration3d==0) then !pjppjp>
        gradssh_u(:,:,0)=gradssh_u(:,:,1)
        gradssh_v(:,:,0)=gradssh_v(:,:,1)
       endif                   !pjppjp>
       do j=2,jmax-1
       do i=2,imax
! gradssh_u(i,j,0) gradssh entre t-1 et t
! gradssh_u(i,j,1) gradssh entre t et t+1
       xy_u(i,j,id_gradssh)=(gradssh_u(i,j,0)+gradssh_u(i,j,1))/hz_u(i,j,1)/dti_lp
! Noter que Somme verticale de dz_u(dnow)/hz_u(1)=1 et donc que l'integrale verticale
! de xy_u(i,j,1) fois le pas de temps est gradssh_u(i,j,0)+gradssh_u(i,j,1)
       enddo
       enddo
      else                           !m°v°m>
       do j=2,jmax-1 ; do i=2,imax
        xy_u(i,j,id_gradssh)=0.
       enddo ; enddo
      endif                          !m°v°m>

      do j=2,jmax-1
      do i=2,imax
! omega_w(:,:,kmax)/=0 impose une C.L. vel_u(i,j,kmax+1,1)
       vel_u(i,j,kmax+1,1)=vel_u(i,j,kmax,1)
      enddo
      enddo

      const2=-0.5*dti_lp                                                !30-09-09

! ANNULER DES TERMES:
!      presgrad_u=0.
!     presgradb_u=0.
!     xy_u(:,:,id_gradssh)=0.
!     tfc1=0. ; tfc2=0. 
!     km_w=0.
!     fric_t=0.
!     dz_u(:,:,:,dafter) =dz_u(:,:,:,dnow)
!     dz_u(:,:,:,dbefore)=dz_u(:,:,:,dnow)

! Selecteur LEAP-FROG ou FORWARD-BACKWARD
      if(timestep_type==timestep_leapfrog)t_=1
      if(timestep_type==timestep_forwbckw)t_=0
!     write(6,*)'t_=',t_
!     stop 'mimi'

      id_coriolis1=5
      if(check0/=anyv3d(imax/2,jmax/2,kmax/2,0)) &
      stop 'Err momentum u anyv3d corrupted'
      if(check7/=anyv3d(imax/2,jmax/2,kmax-1,id_udiv)) then 
       write(10+par%rank,*)'Err momentum check7/=anyv3d' 
       write(10+par%rank,*)'check7=',check7
       write(10+par%rank,*)'anyv3d(imax/2,jmax/2,kmax-1,id_udiv)=',anyv3d(imax/2,jmax/2,kmax-1,id_udiv)
       write(10+par%rank,*)'id_udiv=',id_udiv
       write(10+par%rank,*)'imax/2,jmax/2,kmax-1 glb',imax/2+par%timax(1),jmax/2+par%tjmax(1),kmax-1
       write(10+par%rank,*)'imax/2,jmax/2,kmax-1 loc',imax/2,jmax/2,kmax-1
       stop 'Err momentum check7/=anyv3d see fort.xxx files'
      endif

      do k=1,kmax   ! debut boucle 2 K

      do j=2,jmax
      do i=1,imax


! Rotation (coriolis + grid stretching) term: (1/2)*fstar*(v+VS)
!$ Rotation factor fstar = H*dxdy*f+H*(v*de2/di-u*de1/dj )
!$ including 0.25* for the 4 points average
      xy_t(i,j,id_coriolis1)=                                          &
                    0.25*dz_t(i,j,k,now)*(coriolis_t(i,j)*dxdy_t(i,j)  &
                                          +0.5*(                       &
       (vel_v(i,j,k,1)+vel_v(i  ,j+1,k,1))*(dy_u(i+1,j  )-dy_u(i,j))   &
      -(vel_u(i,j,k,1)+vel_u(i+1,j  ,k,1))*(dx_v(i  ,j+1)-dx_v(i,j)))) &

                     *( vel_v(i,j,k,1)+      vel_v(i,j+1,k,1)          &
                 +velstokes_v(i,j,k,1)+velstokes_v(i,j+1,k,1))

      enddo
      enddo

! ANNULER DES TERMES:
!     xy_t(:,:,id_coriolis1)=0.

      do j=2,jmax-1 ! debut boucle 2 J
      do i=2,imax   ! debut boucle 2 I

! Au membre de droite du systeme tridiagonal:
      tridia_in(i,j,k,4)=( &  !pppp>

       dz_u(i,j,k,dbefore)*( vel_u(i,j,k,2)                   & !05-01-16

                        +tfilterfb*( vel_u(i,j,k,before2)         &
                                    -vel_u(i,j,k,before ) ))      &
                                 
      +dz_u(i,j,k,dnow  )*tfilterlf*( vel_u(i,j,k,before2)         &
                                     -vel_u(i,j,k,before ))        &

      +dz_u(i,j,k,dafter)*( (tfilterlf-tfilterfb)*vel_u(i,j,k,now   )  &
                                     -tfilterlf *vel_u(i,j,k,before)) & 

         -anyv3d(i,j,k,id_udiv)     & !01-05-19

!      +dti_lp*(                                                      &
       +dti_lp*( & !xxxx>                                !13-08-14


! Coriolis
      +(xy_t(i,j,id_coriolis1)+xy_t(i-1,j,id_coriolis1))/dxdy_u(i,j)  &

! Surface stress with a depth distribution function:
       +wstress_u(i,j,1)*(anyv3d(i,j,k,0)+anyv3d(i-1,j,k,0))      &    !30-03-12

       +dz_u(i,j,k,dnow  )*( & !pmxpmx>

! 2d pressure gradient
         +xy_u(i,j,id_gradssh)                             &  
! Stokes forces:                                           & !
              +stokesforces_u(i,j,k)/dx_u(i,j)             & !
! Rappel vers ogcm:                                        & !
              +(sponge_u(i,j,1)*(velobc_u(i,j,k,1)         & !
!                                  -vel_u(i,j,k,now  )     & !
                                                       ))  &
                          ) & !pmxpmx>

! Baroclinic pressure gradient
          -0.5*( dz_u(i,j,k,1)*presgrad_u(i,j,k,1)            & !28-04-20
                +dz_u(i,j,k,0)*presgrad_u(i,j,k,0))/dx_u(i,j) & !10-11-17


                             ) & !xxxx>
           )*mask_u(i,j,k) !pppp>

! Turbulence:
! coef1 --->
      tridia_in(i,j,k,1)=                                             & !15-02-10
         mask_u(i,j,k)                                                &
                    *const2*(   km_w(i,j,k  )+   km_w(i-1,j,k  ))/    & !15-02-10
                            (depth_u(i,j,k  )-depth_u(i  ,j,k-1))

! coef3 --->
      tridia_in(i,j,k,3)=                                             &
         mask_u(i,j,k)                                                &
                    *const2*(   km_w(i,j,k+1)+   km_w(i-1,j,k+1))/    & !15-02-10
                            (depth_u(i,j,k+1)-depth_u(i  ,j,k  ))



      enddo ! fin boucle 2 I
      enddo ! fin boucle 2 J
      enddo ! fin boucle 2 K

! Augmenter le melange entre la couche kmin et la couche kmerged
      x0=const2*100.
      do j=2,jmax-1 ; do i=2,imax   
!      if(kmerged_u(i,j)>1) then ! (°L°) !>
!       k=kmerged_u(i,j) ! Note: si kmerged=kmin alors pas de fusion....
!                        ! La boucle en suivant annule tridia_in(i,j,kmin_u(i,j),1)
!       do k=kmin_u(i,j)+1,kmerged_u(i,j)
!       do k=2,kmerged_u(i,j)
        do k=2,kmin_u(i,j) !04-12-19 ! on ne moyenne plus kmin et kmerged
!        tridia_in(i,j,k  ,1)=x0/(depth_u(i,j,k)-depth_u(i,j,k-1))*mask_u(i,j,kmax)
! La precedente ligne pouvant conduire A des quantitEs non calculables
! par rapport A la precision de la machine dans l'etape du solveur, on
! adopte un autre dimensionnement du melange dans la couche fusionnee   
         tridia_in(i,j,k,1)=-100.*mask_u(i,j,kmax) !21-02-19
         tridia_in(i,j,k-1,3)=tridia_in(i,j,k,1)
        enddo
!      endif                     ! (°L°) !>
      enddo         ; enddo

!......................................................
! Bottom and surface boundary conditions for coef1 & coef3:
      do j=2,jmax-1
      do i=2,imax
       tridia_in(i,j,kmax       ,3)=0.  !surface
!      tridia_in(i,j,kmin_u(i,j),1)=0.  !bottom
       tridia_in(i,j,1          ,1)=0.  !bottom
      enddo
      enddo
!......................................................

! coef2 --->
      if(flag_merged_levels==0.and.flag_z0_macro==1) then !m°v°m> !11-06-22

! CAS PARTICULIER de la rugosite de macro echelle:
! https://docs.google.com/document/d/1QkYuJuw6InVor7g8NKgxRPjSJPrS6TLwC4kxOggG3c4/edit?usp=sharing
! coef2 --->

      tfb0=1-tfilterfb
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax

       tridia_in(i,j,k,2)=dz_u(i,j,k,dafter)*tfb0   & !15-02-10!05-11-10

         +mask_u(i,j,k)*(  & !mmmm>
                         -tridia_in(i,j,k,1)-tridia_in(i,j,k,3)    &
                         +dti_lp*sponge_u(i,j,1)*dz_u(i,j,k,dafter) &

! friction sur plusieurs niveaux
        +dti_lp*0.5*(fric_t(i,j)+fric_t(i-1,j))*(sigma_fric_wu(i,j,k+1)-sigma_fric_wu(i,j,k)) &

                        )    !mmmm>

      enddo ; enddo ; enddo

      else                                                !m°v°m>

! coef2 --->
! CAS GENERAL:
      tfb0=1-tfilterfb
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
       tridia_in(i,j,k,2)=dz_u(i,j,k,dafter)*tfb0   & !15-02-10!05-11-10
         +mask_u(i,j,k)*(  & !mmmm>
                         -tridia_in(i,j,k,1)-tridia_in(i,j,k,3)    &
                         +dti_lp*sponge_u(i,j,1)*dz_u(i,j,k,dafter) &
                        )    !mmmm>
      enddo ; enddo ; enddo
!......................................................
! Bottom boundary conditions for coef2:
      const7=dti_lp/rho
      const8=dti_lp*0.5
      do j=2,jmax-1 ; do i=2,imax
! Distribution homogene du stress de fond sur les couches fusionnees:
! Chaque couche recoit une fraction (dz/(ztop-zfond) du stress total:
       k1=kmerged_u(i,j) !17-10-17
       do k=1,k1
! v269: le stress est appliquE entre le fond et -h_u+dz_u(kmerged) !04-12-19
        tridia_in(i,j,k,2)=                              &
        tridia_in(i,j,k,2)                               &
       +const8*(fric_t(i,j)+fric_t(i-1,j))               &
         *max(min(0.5*(depth_w(i-1,j,k+1)+depth_w(i,j,k+1)),-h_u(i,j)+dz_u(i,j,k1,dnow)) & !04-12-19
                 -0.5*(depth_w(i-1,j,k  )+depth_w(i,j,k  )),0.)                          &
       /dz_u(i,j,k1,dnow)                                  
! Note: l'expression ci-dessus dans max(min( est equivalent A dz_u(k) tant que la facette superieure
! de la couche k est en dessous de -h+dz(kmerged) et sinon
! alors l'epaisseur retenue est zup-zdown avec zup=-h+dz(kmerged) et zdown est z de la facette inferieure
! la couche k
       enddo
      enddo ; enddo

      endif                                               !m°v°m>

!......................................................
! Add wave bottom momentum injection in u equation
      if(iwve==1)call momentum_add_wstresb_u

! Completer la matrice avec les coef de l'advection verticale implicite:
! Dans le cas flag_timesplitting_adve_uv=1:
! A la prmiere iteration les champs u,v,omega,dz  necessaires A l'advection n'ont
! pas ete calculEs, on ne calcule donc pas l'advection 
      flag=1
      if(flag_1dv==1)flag=0                                      ! pas d'advection si modele 1D
      if(flag_timesplitting_adve_uv==1.and.iteration3d==0)flag=0 ! pas d'advection si flag_timesplitting_adve_uv=1 et iteration=0 !03-09-18
      if(flag==1)call momentum_imp_adv_u !01-02-17!15-02-17


! Resoudre le systeme:
 163  call tridiagonalsolver(2,2,imax,2,jmax-1,kmax)                     !30/03/06

      do k=1,kmax
      do j=2,jmax-1
      do i=2,imax
       vel_u(i,j,k,2)=tridia_out(i,j,k)*mask_u(i,j,k) &
                                    *wetmask_u(i,j) !06-04-22
      enddo
      enddo
      enddo

!     if(mod(iteration3d,10)==0) then !>>>
!      i=201-par%timax(1)
!      j=216-par%tjmax(1)
!      if(i>2.and.i<imax.and.j>1.and.j<jmax) then
!        write(10+par%rank,*)'---',iteration3d,kmin_u(i,j),kmergedr4_u(i,j),kmerged_u(i,j)-1
!        do k=kmax,kmin_u(i,j),-1
!         write(10+par%rank,*)  &
!                       real(depth_u(i,j,k)) &
!             ,real(vel_u(i,j,k,2)) &
!             ,real(km_w(i,j,k  )+km_w(i-1,j,k  )) &
!             ,real(km_w(i,j,k  )) &
!             ,real(km_w(i-1,j,k  ))
!        enddo
!      endif
!      j=jmax/2
!      k=kmax/2
!      write(6,*)'rap',vel_u(i,j,k,2)
!      stop 'darta'       
!     endif                      !>>>

! Verifier la divergence du courant (colonne 1 = colonne 3 /= colonne 2)
!     i=imax/2
!     j=jmax/2
!     k=kmax/2
!     write(500+par%rank,*)                                     &
!           -omega_w(i    ,j    ,k+kpt,1  )                     &
!           +omega_w(i    ,j    ,k-kmt,1  )                     &
!       -( veldydz_u(i+ipt,j    ,k    ,now)                     &
!         -veldydz_u(i-imt,j    ,k    ,now)                     &
!         +veldxdz_v(i    ,j+jpt,k    ,now)                     &
!         -veldxdz_v(i    ,j-jmt,k    ,now)    )/dxdy_t(i,j)    &
!       ,-(     dz_t(i    ,j    ,k    ,after )                  &
!              -dz_t(i    ,j    ,k    ,before) )/dti_lp         &
!       ,-(     dz_t(i    ,j    ,k    ,dafter )                 &
!              -dz_t(i    ,j    ,k    ,dbefore) )/dti_lp    

! On commente cette ligne car on a augmente le melange en amont
! afin que le frottement sur le fond implicite tienne compte de
! ce que le cd est fusionne

!$ v component of the velocity:
      id_gradssh=1
      if(flag_timesplitting==1) then !m°v°m> !29-08-18
       do j=2,jmax
       do i=2,imax-1
! gradssh_u(i,j,0) gradssh entre t-1 et t
! gradssh_u(i,j,1) gradssh entre t et t+1
       xy_v(i,j,id_gradssh)=(gradssh_v(i,j,0)+gradssh_v(i,j,1))/hz_v(i,j,1)/dti_lp
! Noter que Somme verticale de dz_u(dnow)/hz_u(1)=1 et donc que l'integrale verticale
! de xy_u(i,j,1) fois le pas de temps est gradssh_u(i,j,0)+gradssh_u(i,j,1)
       enddo
       enddo
      else                           !m°v°m> !29-08-18
       do j=2,jmax ; do i=2,imax-1
        xy_v(i,j,id_gradssh)=0.
       enddo ; enddo
      endif                          !m°v°m> !29-08-18

      do j=2,jmax
      do i=2,imax-1
! omega_w(:,:,kmax)/=0 impose une C.L. vel_v(i,j,kmax+1,1)
       vel_v(i,j,kmax+1,1)=vel_v(i,j,kmax,1)
      enddo
      enddo

! ANNULER DES TERMES:
!      presgrad_v=0.
!     presgradb_v=0.
!      xy_v(:,:,id_gradssh)=0.
!     tfc1=0. ; tfc2=0. 
!     write(6,*)'tfc1,tfc2',tfc1,tfc2
!     stop 'coco'
!     km_w=0.
!     fric_t=0.
!     dz_v(:,:,:,dafter) =dz_v(:,:,:,dnow)
!     dz_v(:,:,:,dbefore)=dz_v(:,:,:,dnow)

      const2=-0.5*dti_lp                                                !30-09-09

      id_coriolis2=6

      if(check0/=anyv3d(imax/2,jmax/2,kmax/2,0)) &
      stop 'Err momentum v anyv3d corrupted'
      if(check8/=anyv3d(imax/2,jmax/2,kmax-1,id_vdiv)) &
      stop 'Err momentum v check8/=anyv3d'


      do k=1,kmax   ! debut boucle 2 K

      do j=1,jmax
      do i=2,imax


! Rotation (coriolis + grid stretching) term: -(1/2)*fstar*(u+US)
!$ Rotation factor fstar = H*dxdy*f+H*(v*de2/di-u*de1/dj )
!$ including 0.25* for the 4 points average
      xy_t(i,j,id_coriolis2)=                                          &
                   -0.25*dz_t(i,j,k,now)*(coriolis_t(i,j)*dxdy_t(i,j)  &
                                          +0.5*(                       &
       (vel_v(i,j,k,1)+vel_v(i  ,j+1,k,1))*(dy_u(i+1,j  )-dy_u(i,j))   &
      -(vel_u(i,j,k,1)+vel_u(i+1,j  ,k,1))*(dx_v(i  ,j+1)-dx_v(i,j)))) &

                     *( vel_u(i,j,k,1)+      vel_u(i+1,j,k,1)          &
                 +velstokes_u(i,j,k,1)+velstokes_u(i+1,j,k,1))

      enddo
      enddo

! ANNULER DES TERMES:
!     xy_t(:,:,id_coriolis2)=0.

      do j=2,jmax   ! debut boucle 2 J
      do i=2,imax-1 ! debut boucle 2 I

! Au membre de droite du systeme tridiagonal:
      tridia_in(i,j,k,4)=(                                            &

       dz_v(i,j,k,dbefore)*( vel_v(i,j,k,2)                   & !05-01-16

                        +tfilterfb*( vel_v(i,j,k,before2)         &
                                    -vel_v(i,j,k,before ) ))      &
                                 
      +dz_v(i,j,k,dnow  )*tfilterlf*( vel_v(i,j,k,before2)         &
                                    -vel_v(i,j,k,before ))        &

      +dz_v(i,j,k,dafter)*( (tfilterlf-tfilterfb)*vel_v(i,j,k,now   )  &
                                     -tfilterlf *vel_v(i,j,k,before)) & 


      -anyv3d(i,j,k,id_vdiv)      & !01-05-19

       +dti_lp*(                                       & !13-08-14


! Coriolis:
        +(xy_t(i,j,id_coriolis2)+xy_t(i,j-1,id_coriolis2))/dxdy_v(i,j)  &

! Surface stress with a depth distribution function:
       +wstress_v(i,j,1)*(anyv3d(i,j,k,0)+anyv3d(i,j-1,k,0))      &

       +dz_v(i,j,k,dnow  )*( & !pmxpmx>

! 2d pressure gradient
         +xy_v(i,j,id_gradssh)                                     &  

! Stokes forces:                                                   & !
                         +stokesforces_v(i,j,k)/dy_v(i,j)          & !
! Rappel vers ogcm:                                                & !
                         +(sponge_v(i,j,1)*(velobc_v(i,j,k,1)      & !
!                                             -vel_v(i,j,k,now  )  & !
                                                               ))  &
                          ) & !pmxpmx>

! Baroclinic pressure gradient:                                    & !
          -0.5*( dz_v(i,j,k,1)*presgrad_v(i,j,k,1)            & !28-04-20
                +dz_v(i,j,k,0)*presgrad_v(i,j,k,0))/dy_v(i,j) &

                               ))*mask_v(i,j,k)

! Turbulence:
! coef1 --->
      tridia_in(i,j,k,1)=                                             &
         mask_v(i,j,k)                                                &
                    *const2*(   km_w(i,j,k  )+   km_w(i,j-1,k  ))/    & !15-02-10
                            (depth_v(i,j,k  )-depth_v(i,j  ,k-1))

! coef3 --->
      tridia_in(i,j,k,3)=                                             &
         mask_v(i,j,k)                                                &
                    *const2*(   km_w(i,j,k+1)+   km_w(i,j-1,k+1))/    & !15-02-10
                            (depth_v(i,j,k+1)-depth_v(i,j  ,k  ))

      enddo ! fin boucle 2 I
      enddo ! fin boucle 2 J
      enddo ! fin boucle 2 K

! Augmenter le melange entre la couche kmin et la couche kmerged
      x0=const2*100.
      do j=2,jmax ; do i=2,imax-1 ! debut boucle 2 I
!      if(kmerged_v(i,j)>1) then ! (oVo) !> !03-05-17
!       k=kmerged_v(i,j)
!       do k=kmin_v(i,j)+1,kmerged_v(i,j)
!       do k=2,kmerged_v(i,j)
        do k=2,kmin_v(i,j) !04-12-19 on ne moyenne pas kmin et kmerged
!        tridia_in(i,j,k  ,1)=x0/(depth_v(i,j,k)-depth_v(i,j,k-1))*mask_v(i,j,kmax)
! La precedente ligne pouvant conduire A des quantitEs non calculables
! par rapport A la precision de la machine dans l'etape du solveur, on
! adopte un autre dimensionnement du melange dans la couche fusionnee
         tridia_in(i,j,k,1)=-100.*mask_v(i,j,kmax) !21-02-19
         tridia_in(i,j,k-1,3)=tridia_in(i,j,k,1)
        enddo
!      endif                     ! (oVo) !>
      enddo       ; enddo

!......................................................
! Bottom and surface boundary conditions for coef1 & coef3:
      do j=2,jmax
      do i=2,imax-1
       tridia_in(i,j,kmax       ,3)=0.
!      tridia_in(i,j,kmin_v(i,j),1)=0.
       tridia_in(i,j,1          ,1)=0.
      enddo
      enddo
!......................................................

! coef2 --->
      if(flag_merged_levels==0.and.flag_z0_macro==1) then !m°v°m> !11-06-22
! CAS PARTICULIER de la rugosite de macro echelle:
! https://docs.google.com/document/d/1QkYuJuw6InVor7g8NKgxRPjSJPrS6TLwC4kxOggG3c4/edit?usp=sharing
! coef2 --->
      tfb0=1-tfilterfb
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
       tridia_in(i,j,k,2)= dz_v(i,j,k,dafter)*tfb0  & 
         +mask_v(i,j,k)*(  & !mmmm>
                         -tridia_in(i,j,k,1)-tridia_in(i,j,k,3)     &
                         +dti_lp*sponge_v(i,j,1)*dz_v(i,j,k,dafter) &

! friction sur plusieurs niveaux:
       +dti_lp*0.5*(fric_t(i,j)+fric_t(i,j-1))*(sigma_fric_wv(i,j,k+1)-sigma_fric_wv(i,j,k)) &

                        )    !mmmm>
      enddo ; enddo ; enddo

      else                                               !m°v°m>

! CAS GENERAL
! coef2 --->
      tfb0=1-tfilterfb
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
       tridia_in(i,j,k,2)= dz_v(i,j,k,dafter)*tfb0  & 
         +mask_v(i,j,k)*(  & !mmmm>
                         -tridia_in(i,j,k,1)-tridia_in(i,j,k,3)    &
                         +dti_lp*sponge_v(i,j,1)*dz_v(i,j,k,dafter) &
                        )    !mmmm>
      enddo ; enddo ; enddo
! Bottom boundary conditions for coef2:
      const7=dti_lp/rho
      const8=dti_lp*0.5
      do j=2,jmax ; do i=2,imax-1
! Distribution homogene du stress de fond sur les couches fusionnees:
       k1=kmerged_v(i,j)
       do k=1,k1
! v269: le stress est appliquE entre le fond et -h_u+dz_u(kmerged) !04-12-19
        tridia_in(i,j,k,2)=                              &
        tridia_in(i,j,k,2)                               &
       +const8*(fric_t(i,j)+fric_t(i,j-1))               &
         *max(min(0.5*(depth_w(i,j-1,k+1)+depth_w(i,j,k+1)),-h_v(i,j)+dz_v(i,j,k1,dnow)) & !04-12-19
                 -0.5*(depth_w(i,j-1,k  )+depth_w(i,j,k  )),0.)                          &
       /dz_v(i,j,k1,dnow)      

       enddo
      enddo ; enddo

      endif                                              !m°v°m>

!......................................................
! Add wave bottom momentum injection in v equation
      if(iwve==1)call momentum_add_wstresb_v

! Completer la matrice avec les coef de l'advection verticale implicite:
! A la prmiere iteration les champs u,v,omega,dz  necessaires A l'advection n'ont
! pas ete calculEs, on ne calcule donc pas l'advection !03-09-18
      flag=1
      if(flag_1dv==1)flag=0                                      ! pas d'advection si modele 1D
      if(flag_timesplitting_adve_uv==1.and.iteration3d==0)flag=0 ! pas d'advection si flag_timesplitting_adve_uv=1 et iteration=0
      if(flag==1)call momentum_imp_adv_v !01-02-17!15-02-17

! Resoudre le systeme:
      call tridiagonalsolver(3,2,imax-1,2,jmax,kmax)                     !30/03/06

      do k=1,kmax
      do j=2,jmax
      do i=2,imax-1
       vel_v(i,j,k,2)=tridia_out(i,j,k)*mask_v(i,j,k) &
                                    *wetmask_v(i,j)     !06-04-22
      enddo
      enddo
      enddo

      if(coef_diss_mangrove>0.)call mangrove_3d_friction !23-06-13
      if(coastal_viscosity>0.) call dragcoef_no_slip_condition_3d !29-06-22

! Melange des couches fusionnEes:
      call vertmix_merged_levels_uv(2) !01-03-19!04-11-21

! Appliquer un profil log dans la couche fusionnee:
!     if(flag_merged_levels==1)call momemtum_logprofile !19-11-20

!     call momentum_uv_minmax(2) !16-01-22

      endif                                           !AAAAAAAAAAAAAAA>

      if(texte_   (1:20)=='before external mode') then !BBBBBBBBBBBBBBB>

      stop 'QUE FAIS JE ICI ?'

      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1

        veldydz_u(i,j,k,1)=(vel_u(i,j,k,1)                            &
                     +velstokes_u(i,j,k,1)                            &
                           )*dy_u(i,j)*dz_u(i,j,k,dnow)

      enddo
      enddo
      enddo


      do k=1,kmax
      do j=1,jmax+1
      do i=1,imax

        veldxdz_v(i,j,k,1)=(vel_v(i,j,k,1)                            &
                     +velstokes_v(i,j,k,1)                            &
                           )*dx_v(i,j)*dz_v(i,j,k,dnow)

      enddo
      enddo
      enddo

      endif                                           !BBBBBBBBBBBBBBB>

      end subroutine momentum_equations

!.......................................................................

      subroutine momentum_zonal_restoring !19-01-16
      use module_principal ; use module_parallele
      integer :: moduliter_=1
!     stop 'momentum_zonal_restoring'

      if(relax_ts<=0.)return
      if(mod(iteration3d,moduliter_)/=0) return

!...................................................................
! Rappel moyenne zonale cas testjet barocline comodo !17-10-13!15-01-16
!...................................................................
! Attention les dimensions de ces tableaux sont communes avec celles du
! rappel de T
!
      if(.not.allocated(iaveraged_in ))allocate(iaveraged_in (jglb,kmax))
      if(.not.allocated(iaveraged_out))allocate(iaveraged_out(jglb,kmax))

! u component:
!     velobc_u=0.5 ; vel_u=-0.5 ! pour verif u
      x2=timeweightobc(vel_id) ; x0=1.-x2
      iaveraged_in=0.
      do k=1,kmax
       do j=2,jmax-1 
        j1=j+par%tjmax(1)
        do i=2,imax
        iaveraged_in(j1,k)=iaveraged_in(j1,k)+mask_i_u(i)*mask_j_u(j)*(&!ooo>!20-01-16
                           x2*velobc_u(i,j,k,2)                        &
                          +x0*velobc_u(i,j,k,0)                        &
                                -vel_u(i,j,k,0)                        &
                                                                      ) !ooo>

        enddo
       enddo
      enddo
      call mpi_allreduce( iaveraged_in         &
                         ,iaveraged_out        &
                         ,jglb*kmax            &
                         ,mpi_double_precision & 
                         ,mpi_sum              &
                         ,par%comm2d ,ierr)


      x0=dti_lp*moduliter_*relax_ts/real(iglb-2) 
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax
       vel_u(i,j,k,2)=vel_u(i,j,k,2)+iaveraged_out(j+par%tjmax(1),k)*x0
!     write(30+par%rank,*)iaveraged_out(j+par%tjmax(1),k)/real(iglb-2) ! pour verif u doit donner 1
      enddo ; enddo ; enddo
!     stop 'verif u'

! v component:
!     velobc_v=0.5 ; vel_v=-0.5 ! pour verif v
      x2=timeweightobc(vel_id) ; x0=1.-x2
      iaveraged_in=0.
      do k=1,kmax
       do j=2,jmax
        j1=j+par%tjmax(1)
        do i=2,imax-1
        iaveraged_in(j1,k)=iaveraged_in(j1,k)+mask_i_v(i)*mask_j_v(j)*(&!ooo>
                           x2*velobc_v(i,j,k,2)                        &
                          +x0*velobc_v(i,j,k,0)                        &
                                -vel_v(i,j,k,0)                        &
                                                                      ) !ooo>
        enddo
       enddo
      enddo
      call mpi_allreduce( iaveraged_in         &
                         ,iaveraged_out        &
                         ,jglb*kmax            &
                         ,mpi_double_precision & 
                         ,mpi_sum              &
                         ,par%comm2d ,ierr)

      x0=dti_lp*moduliter_*relax_ts/real(iglb-2)
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
       vel_v(i,j,k,2)=vel_v(i,j,k,2)+iaveraged_out(j+par%tjmax(1),k)*x0
!     write(40+par%rank,*)iaveraged_out(j+par%tjmax(1),k)/real(iglb-2) ! pour verif v doit donner 1
      enddo ; enddo ; enddo
!     stop 'verif v'

!...................................................................
! Rappel moyenne zonale cas testjet barocline comodo !17-10-13
!...................................................................

      end subroutine momentum_zonal_restoring !19-01-16

!....................................................................
#ifdef bidon
      subroutine momentum_simplified
      use module_principal 
      implicit none

!     if(iteration3d==10) then
! Termes annulEs dans momentum_equation:
           cdb_t=0.
           vel_u=0.
!          vel_v=0.
      coriolis_t=0.
            km_w=0.
       wetmask_u=1.
       wetmask_v=1.
       wetmask_t=1.
       gradssh_u=0.
       gradssh_v=0.
      presgrad_u=0.
      presgrad_v=0.
      velobc_u=0.
      velobc_v=0.
      sponge_u=0.
      sponge_v=0.
      tfilterfb=0.
      tfilterlf=0.

      sum0=0.
      sum2=0.
      sum3=0.
      sum4=0.
      do k=1,kmax
      do j=1,jmax
      do i=1,imax+1
       sum0=sum0+vel_u(i,j,k,0)*dz_u(i,j,k,dbefore)*dxdy_u(i,j)*mask_u(i,j,k)
       sum2=sum2+vel_u(i,j,k,2)*dz_u(i,j,k,dafter)*dxdy_u(i,j)*mask_u(i,j,k)
       sum3=sum3+vel_v(i,j,k,0)*dz_v(i,j,k,dbefore)*dxdy_v(i,j)*mask_v(i,j,k)
       sum4=sum4+vel_v(i,j,k,2)*dz_v(i,j,k,dafter)*dxdy_v(i,j)*mask_v(i,j,k)
      enddo
      enddo
      enddo
      write(6,*)'sum0,sum2=',sum0,sum2
      write(6,*)'sum3,sum4=',sum3,sum4

!     endif
      end subroutine momentum_simplified
#endif
!....................................................................
      subroutine momentum_wstresbprofile !28-09-16
#ifdef stokes
      use module_principal ; use module_parallele
      implicit none

      if(id_wb/=1)stop 'id_wb/=1'

! Vertical distribution of the wave bottom injection of momentum
      do j=1,jmax ; do i=1,imax
       xy_t(i,j,1)=2.*pi/max(ubw(i,j)*t_wave_t(i,j,1),2.*pi) ! inverse de la demi-excursion
                                                             ! Michaud et al, OS, 2012, Eq27
!      xy_t(i,j,1)=1.
      enddo       ; enddo
      x0=0.5*dti_lp/rho
      do k=1,kmax+1 ; do j=1,jmax ; do i=1,imax
        anyv3d(i,j,k,id_wb)=x0*exp(-(depth_w(i,j,k)+h_w(i,j))*xy_t(i,j,1))
      enddo         ; enddo       ; enddo
      check1=anyv3d(imax/2,jmax/2,kmax/2,id_wb)
#endif
      end subroutine momentum_wstresbprofile
!....................................................................
      subroutine momentum_add_wstresb_u !28-09-16
#ifdef stokes
      use module_principal ; use module_parallele
      implicit none

! Add wave bottom momentum injection in u equation
      if(check1/=anyv3d(imax/2,jmax/2,kmax/2,id_wb)) &
      stop 'Err momentum_add_wstresb_u anyv3d id_wb corrupted'

      do k=1,kmax ; do j=2,jmax-1  ; do i=2,imax

         tridia_in(i,j,k,4)= &
         tridia_in(i,j,k,4)  &
        +wstresb_u(i,j)      &
        *wetmask_u(i,j)*( anyv3d(i,j,k  ,id_wb)+anyv3d(i-1,j,k  ,id_wb) &
                         -anyv3d(i,j,k+1,id_wb)-anyv3d(i-1,j,k+1,id_wb))
! Note anyv3d(:,:,k,0)-anyv3d(:,:,k+1,0) > 0 c'est donc bien dans ce
! sens que se fait la diffence

      enddo         ; enddo       ; enddo

#endif
      end subroutine momentum_add_wstresb_u
!....................................................................
      subroutine momentum_add_wstresb_v !28-09-16
#ifdef stokes
      use module_principal ; use module_parallele
      implicit none

      if(check1/=anyv3d(imax/2,jmax/2,kmax/2,id_wb)) &
      stop 'Err momentum_add_wstresb_v anyv3d id_wb corrupted'

! Add wave bottom momentum injection in v equation
!     sum1=0.
      do k=1,kmax ; do j=2,jmax  ; do i=2,imax-1

         tridia_in(i,j,k,4)= &
         tridia_in(i,j,k,4)  &
        +wstresb_v(i,j)      &
        *wetmask_v(i,j)*( anyv3d(i,j,k  ,id_wb)+anyv3d(i,j-1,k  ,id_wb) &
                         -anyv3d(i,j,k+1,id_wb)-anyv3d(i,j-1,k+1,id_wb))
! Note anyv3d(:,:,k,0)-anyv3d(:,:,k+1,0) > 0 c'est donc bien dans ce
! sens que se fait la diffence
!       if(i==imax/2.and.j==jmax/2.and.mask_v(i,j,k)==1) then
!       sum1=sum1+( anyv3d(i,j,k  ,id_wb)+anyv3d(i,j-1,k  ,id_wb) &
!                  -anyv3d(i,j,k+1,id_wb)-anyv3d(i,j-1,k+1,id_wb))
!       write(10+par%rank,*)depth_v(i,j,k)    &
!                      ,( anyv3d(i,j,k  ,id_wb)+anyv3d(i,j-1,k  ,id_wb) &
!                        -anyv3d(i,j,k+1,id_wb)-anyv3d(i,j-1,k+1,id_wb)),sum1*rho/dti_lp
!       endif

      enddo         ; enddo       ; enddo

#endif
      end subroutine momentum_add_wstresb_v

!.......................................................................


      subroutine momentum_imp_adv_u !01-02-17
      use module_principal ; use module_parallele
      implicit none
      integer :: id_wi_=6       

      if(check3/=anyv3d(imax/2,jmax/2,kmax/2,id_wi_)) then
       write(6,*)'check3/=anyv3d(id_wi_)',check3,anyv3d(imax/2,jmax/2,kmax/2,id_wi_),par%rank
       stop 'momentum_imp_adv_u avy3d corrupted'
      endif

! Implicit "up2" (upwind 2 points) part of the vertical advection: !20-08-14
      do k=1,kmax-1 !kmax !09-04-15
      do j=2,jmax-1 ; do i=2,imax  

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


      k=kmax
      do j=2,jmax-1 ; do i=2,imax  

      tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
       (-anyv3d(i,j,k,id_wi_)-abs(anyv3d(i,j,k,id_wi_)))

! Condition de gradient nul
! Pour obtenir les lignes suivantes, on part des lignes precedentes (cas "dans la couche")
! A1*T(k-1)+A2*T(k)+A3*T(k+1) et on remplace T(k+1) par T(k) pour avoir la condition de
! gradient nul et deduire les nouveaux coefficients: A1=A1 ; A2=A2+A3 ; A3=0. 
!     tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
!     +2.*anyv3d(i,j,k+1,id_wi_)                &
!     +( -anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

!......
! Condition courant nul: !07-03-22
      tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(     &
       ( anyv3d(i,j,k+1,id_wi_)+abs(anyv3d(i,j,k+1,id_wi_))) &
      +(-anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

!     tridia_in(i,j,k,3)=tridia_in(i,j,k,3)+0. inutile donc commentee
!......

      enddo ; enddo

      end subroutine momentum_imp_adv_u 

!...................................................................

      subroutine momentum_imp_adv_v !01-02-17
      use module_principal ; use module_parallele
      implicit none
      integer :: id_wi_=7       

! Implicit "up2" (upwind 2 points) part of the vertical advection: !20-08-14
      if(check4/=anyv3d(imax/2,jmax/2,kmax/2,id_wi_)) &
      stop 'momentum_imp_adv_v avy3d corrupted'

      do k=1,kmax-1 !kmax !09-04-15
      do j=2,jmax ; do i=2,imax-1  

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


      k=kmax
      do j=2,jmax ; do i=2,imax-1  

      tridia_in(i,j,k,1)=tridia_in(i,j,k,1)+0.5*dti_lp*                &
       (-anyv3d(i,j,k,id_wi_)-abs(anyv3d(i,j,k,id_wi_)))

! Condition de gradient nul
! Pour obtenir les lignes suivantes, on part des lignes precedentes (cas "dans la couche")
! A1*T(k-1)+A2*T(k)+A3*T(k+1) et on remplace T(k+1) par T(k) pour avoir la condition de
! gradient nul et deduire les nouveaux coefficients: A1=A1 ; A2=A2+A3 ; A3=0. 
!     tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(               &
!     +2.*anyv3d(i,j,k+1,id_wi_)                &
!     +( -anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

!......
! Condition courant nul: !07-03-22
      tridia_in(i,j,k,2)=tridia_in(i,j,k,2)+0.5*dti_lp*(     &
       ( anyv3d(i,j,k+1,id_wi_)+abs(anyv3d(i,j,k+1,id_wi_))) &
      +(-anyv3d(i,j,k  ,id_wi_)+abs(anyv3d(i,j,k  ,id_wi_))) )

!     tridia_in(i,j,k,3)=tridia_in(i,j,k,3)+0. inutile donc commentee
!......

      enddo ; enddo

      end subroutine momentum_imp_adv_v 

!...................................................................

      subroutine momentum_crushed_cells
      use module_principal ; use module_parallele
      implicit none

      stop 'momentum_crushed_cells'

! Les couches ecrasees (epaisseur nulle) prennent la valeur de la premiere couche non nulle
      do j=2,jmax-1
      do i=2,imax
       do k=1,kmin_u(i,j)-1
        vel_u(i,j,k,2)=vel_u(i,j,kmin_u(i,j),2)
       enddo
      enddo
      enddo
      do j=2,jmax
      do i=2,imax-1
       do k=1,kmin_v(i,j)-1
        vel_v(i,j,k,2)=vel_v(i,j,kmin_v(i,j),2)
       enddo
      enddo
      enddo

      end subroutine momentum_crushed_cells

!...................................................................

      subroutine momentum_replay_ncumax(id_veldydz_u_,id_veldxdz_v_)
      use module_principal ; use module_parallele
      implicit none
      integer :: id_veldydz_u_,id_veldxdz_v_

! Precedement x2 a ete identifie comme le nombre de courant max. A quel
! point de grille correspond cette valeur? On rejoue l'algo et on 
! cherche le point correspondant A x2

      do k=1,kmax !04-11-21
      do j=2,jmax-1 ; do i=2,imax
!       do k=kmerged_u(i,j)+1,kmax

         x1=       max(max(max(abs(anyv3d(i  ,j  ,k,id_veldydz_u_))  &
                              ,abs(anyv3d(i-1,j  ,k,id_veldydz_u_))) &
                              ,abs(anyv3d(i  ,j  ,k,id_veldxdz_v_))) &
                              ,abs(anyv3d(i  ,j+1,k,id_veldxdz_v_))) & 
                    *dti_lp                                          &
                    *wetmask_u(i,j)                                  &
                      /(dxdy_u(i,j)*min(dz_u(i,j,k,dbefore),dz_u(i,j,k,dafter))) 

          if(x1==x2) then !>>>
           write(10+par%rank,*)'momentum_replay_ncumax regular cells'
           write(10+par%rank,*)'x1=',x1
           write(10+par%rank,*)'i,j,k loc',i,j,k
           write(10+par%rank,*)'i,j glb',i+par%timax(1),j+par%tjmax(1)
           write(10+par%rank,*)'kmin_u(i,j)',kmin_u(i,j)
           write(10+par%rank,*)'kmerged_u(i,j)',kmerged_u(i,j)
           write(10+par%rank,*)'wetmask_u(i,j)',wetmask_u(i,j)
           write(10+par%rank,*)'anyv3d(i  ,j  ,k,id_veldydz_u_)',anyv3d(i  ,j  ,k,id_veldydz_u_)
           write(10+par%rank,*)'anyv3d(i-1,j  ,k,id_veldydz_u_)',anyv3d(i-1,j  ,k,id_veldydz_u_)
           write(10+par%rank,*)'anyv3d(i  ,j  ,k,id_veldxdz_v_)',anyv3d(i  ,j  ,k,id_veldxdz_v_)
           write(10+par%rank,*)'anyv3d(i  ,j+1,k,id_veldxdz_v_)',anyv3d(i  ,j+1,k,id_veldxdz_v_)
           write(10+par%rank,*)'veldydz_u(i-1:i+1,j,k,1)',real(veldydz_u(i-1:i+1,j,k,1))
           write(10+par%rank,*)'veldxdz_v(i-1:i,j  ,k,1)',real(veldxdz_v(i-1:i,j,k,1))
           write(10+par%rank,*)'veldxdz_v(i-1:i,j+1,k,1)',real(veldxdz_v(i-1:i,j+1,k,1))
           write(10+par%rank,*)'vel_u(i-1:i+1,j,k,1)',real(vel_u(i-1:i+1,j,k,1))
           write(10+par%rank,*)'vel_v(i-1:i,j  ,k,1)',real(vel_v(i-1:i,j,k,1))
           write(10+par%rank,*)'vel_v(i-1:i,j+1,k,1)',real(vel_v(i-1:i,j+1,k,1))
           write(10+par%rank,*)'dz_u(i-1:i+1,j,k,0:2)',real(dz_u(i-1:i+1,j,k,0:2))
           write(10+par%rank,*)'dz_v(i-1:i,j  ,k,0:2)',real(dz_v(i-1:i,j,k,0:2))
           write(10+par%rank,*)'dz_v(i-1:i,j+1,k,0:2)',real(dz_v(i-1:i,j+1,k,0:2))
! Note:
!      anyv3d(i,j,k,id_veldydz_u_)=0.5*(veldydz_u(i,j,k,1)+veldydz_u(i+1,j,k,1)) 
!      anyv3d(i,j,k,id_veldxdz_v_)=0.5*(veldxdz_v(i-1,j,k,1)+veldxdz_v(i,j,k,1))
          endif           !>>>

      enddo       ; enddo
      enddo
#ifdef bidon
! Niveaux <= kmerged
      do j=2,jmax-1 ; do i=2,imax
        sum1=0. ; sum2=0. ; sum3=0. ; sum4=0. ; sum5=0. ; sum6=0.
        do k=kundermin_u(i,j),kmerged_u(i,j) !07-12-17
         sum1=sum1+anyv3d(i  ,j  ,k,id_veldydz_u_)
         sum2=sum2+anyv3d(i-1,j  ,k,id_veldydz_u_)
         sum3=sum3+anyv3d(i  ,j  ,k,id_veldxdz_v_)
         sum4=sum4+anyv3d(i  ,j+1,k,id_veldxdz_v_)
         sum5=sum5     +dz_u(i,j,k,dbefore)
         sum6=sum6     +dz_u(i,j,k,dafter)
        enddo
         x1=       max(max(max(abs(sum1)  &
                              ,abs(sum2)) &
                              ,abs(sum3)) &
                              ,abs(sum4)) & 
                    *dti_lp                                  &
                    *wetmask_u(i,j)                          &
                      /(dxdy_u(i,j)*min(sum5,sum6))            ! methode 2
                 ! fermeture max
          if(x1==x2) then !>>>
           write(10+par%rank,*)'momentum_replay_ncumax merged cells'
           write(10+par%rank,*)'x1=',x1
           write(10+par%rank,*)'k,kundermin_u(i,j),kmerged_u(i,j)',k,kundermin_u(i,j),kmerged_u(i,j)
           write(10+par%rank,*)'sum1',sum1
           write(10+par%rank,*)'sum2',sum2
           write(10+par%rank,*)'sum3',sum3
           write(10+par%rank,*)'sum4',sum4
          endif           !>>>
      enddo       ; enddo
#endif

      end subroutine momentum_replay_ncumax
!..........................................................
      subroutine momentum_replay_ncvmax(id_veldydz_u_,id_veldxdz_v_)
      use module_principal ; use module_parallele
      implicit none
      integer :: id_veldydz_u_,id_veldxdz_v_


! Niveaux standard:
      do k=1,kmax !04-11-21
      do j=2,jmax ; do i=2,imax-1
!       do k=kmerged_v(i,j)+1,kmax
         x1=       max(max(max(abs(anyv3d(i  ,j  ,k,id_veldydz_u_))  &
                              ,abs(anyv3d(i+1,j  ,k,id_veldydz_u_))) &
                              ,abs(anyv3d(i  ,j  ,k,id_veldxdz_v_))) &
                              ,abs(anyv3d(i  ,j-1,k,id_veldxdz_v_))) & 
                    *dti_lp                                          &
                    *wetmask_v(i,j)                                  &
                      /(dxdy_v(i,j)*min(dz_v(i,j,k,dbefore),dz_v(i,j,k,dafter))) 

          if(x1==x2) then !>>>
           write(10+par%rank,*)'momentum_replay_ncvmax regular cells'
           write(10+par%rank,*)'x1=',x1
           write(10+par%rank,*)'k,kmerged_v(i,j)+1',k,kmerged_v(i,j)+1
           write(10+par%rank,*)'i,j loc',i,j
           write(10+par%rank,*)'i,j glb',i+par%timax(1),j+par%tjmax(1)
           write(10+par%rank,*)'anyv3d(i  ,j ,k,id_veldydz_u_)',anyv3d(i  ,j  ,k,id_veldydz_u_)
           write(10+par%rank,*)'anyv3d(i+1,j ,k,id_veldydz_u_)',anyv3d(i+1,j  ,k,id_veldydz_u_)
           write(10+par%rank,*)'anyv3d(i  ,j ,k,id_veldxdz_v_)',anyv3d(i  ,j  ,k,id_veldxdz_v_)
           write(10+par%rank,*)'anyv3d(i ,j-1,k,id_veldxdz_v_)',anyv3d(i  ,j-1,k,id_veldxdz_v_)
          endif           !>>>


      enddo       ; enddo
      enddo
#ifdef bidon
! Niveaux <= kmerged
      do j=2,jmax ; do i=2,imax-1
        sum1=0. ; sum2=0. ; sum3=0. ; sum4=0. ; sum5=0. ; sum6=0.
        do k=kundermin_v(i,j),kmerged_v(i,j) !07-12-17
         sum1=sum1+anyv3d(i  ,j  ,k,id_veldydz_u_)
         sum2=sum2+anyv3d(i+1,j  ,k,id_veldydz_u_)
         sum3=sum3+anyv3d(i  ,j  ,k,id_veldxdz_v_)
         sum4=sum4+anyv3d(i  ,j-1,k,id_veldxdz_v_)
         sum5=sum5  +dz_v(i,j,k,dbefore)
         sum6=sum6  +dz_v(i,j,k,dafter)
        enddo
         x1=       max(max(max(abs(sum1)  &
                              ,abs(sum2)) &
                              ,abs(sum3)) &
                              ,abs(sum4)) & 
                    *dti_lp                                  &
                    *wetmask_v(i,j)                          &
                      /(dxdy_v(i,j)*min(sum5,sum6))           ! methode 2

          if(x1==x2) then !>>>
           write(10+par%rank,*)'momentum_replay_ncvmax merged cells'
           write(10+par%rank,*)'x1=',x1
           write(10+par%rank,*)'k,kundermin_v(i,j),kmerged_v(i,j)',k,kundermin_v(i,j),kmerged_v(i,j)
           write(10+par%rank,*)'sum1',sum1
           write(10+par%rank,*)'sum2',sum2
           write(10+par%rank,*)'sum3',sum3
           write(10+par%rank,*)'sum4',sum4
          endif           !>>>

      enddo       ; enddo
#endif

      end subroutine momentum_replay_ncvmax

!..........................................................

      subroutine momentum_bottom_current !06-12-19
      use module_principal ; use module_parallele
      implicit none

! Courant de fond

!............................
! Cas des couches sigma standard
! velbot=vel(k=1)

!............................
! Cas des couches fusionneEs:
! La vitesse de fond est la vitesse moyennee entre z=-h et z=-h+dz(kmerged)
! Details dans:
! https://docs.google.com/presentation/d/1aUGTYvHSOEZIDCHzs-EBPKCYGJ9cXM1ntRLuaw9Ra0w/edit?folder=0BxfDfpz8eP5VU2RKRGZjU01IMkU#slide=id.g6c10187b84_0_1

! La partie entiere de kmergedr4 = kmerge (sauf cas particuliers avec kmerged-kmin>1)
! La partie decimale de kmergedr4 indique la fraction de la couche kmerge qui intersecte la couche de fond
! definie comme la couche entre le fond et une hauteur egale A dz(kmerged) cette derniere etant
! continue sur l'horizontale

      if(flag_merged_levels==0.and.flag_z0_macro==0) then !SIGMA>
       do j=1,jmax ; do i=1,imax+1
        velbot_u(i,j)=0.5*(vel_u(i,j,1,0)+vel_u(i,j,1,1)) !12-12-19
       enddo         ; enddo
       do j=1,jmax+1 ; do i=1,imax
        velbot_v(i,j)=0.5*(vel_v(i,j,1,0)+vel_v(i,j,1,1))
       enddo         ; enddo
      endif                                               !SIGMA>

      if(flag_merged_levels==0.and.flag_z0_macro==1) then !SIGMA-et-Z0MACRO> !11-06-22
! https://docs.google.com/document/d/1QkYuJuw6InVor7g8NKgxRPjSJPrS6TLwC4kxOggG3c4/edit?usp=sharing
       k0=0
       do k=1,kmax
        do j=1,jmax ; do i=1,imax+1
         velbot_u(i,j)=                    &
         velbot_u(i,j)*k0                  &
      +0.5*(vel_u(i,j,k,0)+vel_u(i,j,k,1)) &
                *( sigma_fric_wu(i,j,k+1)  &
                  -sigma_fric_wu(i,j,k  ))
        enddo         ; enddo
        do j=1,jmax+1 ; do i=1,imax
         velbot_v(i,j)=                    &
         velbot_v(i,j)*k0                  &
      +0.5*(vel_v(i,j,k,0)+vel_v(i,j,k,1)) &
                *( sigma_fric_wv(i,j,k+1)  &
                  -sigma_fric_wv(i,j,k  ))
        enddo         ; enddo
        k0=1
       enddo ! boucle sur K
      endif                                               !SIGMA-et-Z0MACRO>

      if(flag_merged_levels==1) then !MERGED-LAYERS>
! velbot_u:
      do j=1,jmax ; do i=1,imax+1
       sum1=0. ; sum2=0.
       do k=kmin_u(i,j),int(kmergedr4_u(i,j))-1
        sum1=sum1+dz_u(i,j,k,dnow)
        sum2=sum2+dz_u(i,j,k,dnow)*0.5*(vel_u(i,j,k,0)+vel_u(i,j,k,1))
       enddo
       k=int(kmergedr4_u(i,j))
       sum1=sum1+(kmergedr4_u(i,j)-int(kmergedr4_u(i,j)))*dz_u(i,j,k,dnow)
       sum2=sum2+(kmergedr4_u(i,j)-int(kmergedr4_u(i,j)))*dz_u(i,j,k,dnow)*0.5*(vel_u(i,j,k,0)+vel_u(i,j,k,1))
       velbot_u(i,j)=sum2/sum1
      enddo         ; enddo

! velbot_v:
      do j=1,jmax+1 ; do i=1,imax
       sum1=0. ; sum2=0.
       do k=kmin_v(i,j),int(kmergedr4_v(i,j))-1
        sum1=sum1+dz_v(i,j,k,dnow)
        sum2=sum2+dz_v(i,j,k,dnow)*0.5*(vel_v(i,j,k,0)+vel_v(i,j,k,1))
       enddo
       k=int(kmergedr4_v(i,j))
       sum1=sum1+(kmergedr4_v(i,j)-int(kmergedr4_v(i,j)))*dz_v(i,j,k,dnow)
       sum2=sum2+(kmergedr4_v(i,j)-int(kmergedr4_v(i,j)))*dz_v(i,j,k,dnow)*0.5*(vel_v(i,j,k,0)+vel_v(i,j,k,1))
       velbot_v(i,j)=sum2/sum1
      enddo         ; enddo

      endif                          !MERGED-LAYERS>

      end subroutine momentum_bottom_current
!.................................................................
      subroutine momemtum_logprofile !19-11-20
      use module_principal ; use module_parallele
      implicit none

! Appliquer un profil log a l'interieur de la "couche de fond" comprise entre kmin et kmergedr4, 
! kmergedr4 etant decimal pour ne couvrir que la portion de la cellule mordant sur la couche de fond
! qui est aussi la couche correspondant A velbot

! vel_u:
      do j=2,jmax-1 ; do i=2,imax
       sum2=0. ; sum3=0.
       do k=kmin_u(i,j),int(kmergedr4_u(i,j))-1
        sum2=sum2+dz_u(i,j,k,dafter)*vel_u(i,j,k,2)
        sum3=sum3+dz_u(i,j,k,dafter)*log((h_u(i,j)+depth_u(i,j,k))/(0.5*(z0_w(i,j)+z0_w(i-1,j))))
       enddo
       k=int(kmergedr4_u(i,j))
       sum2=sum2+(kmergedr4_u(i,j)-int(kmergedr4_u(i,j)))*dz_u(i,j,k,dafter)*vel_u(i,j,k,2)
       sum3=sum3+(kmergedr4_u(i,j)-int(kmergedr4_u(i,j)))*dz_u(i,j,k,dafter)*log((h_u(i,j)+depth_u(i,j,k))/(0.5*(z0_w(i,j)+z0_w(i-1,j))))

! Explications sur une alternative possible basee sur la definition de la profondeur suivante:
! depth_u(i,j,k)+((kmergedr4_u(i,j)-int(kmergedr4_u(i,j)))-0.5)*dz_u(i,j,k,dafter)) :
! si kmergedr4_u(i,j)-int(kmergedr4_u(i,j))=0 alors la profondeur est depth_u(k)-0.5dz_u(k) soit la facette inferieure de la cellule k
! si kmergedr4_u(i,j)-int(kmergedr4_u(i,j))=1 alors la profondeur est depth_u(k)+0.5dz_u(k) soit la facette superieure de la cellule k

! x0 est l'amplitude pour u=x0*log((h+z)/z0)
       x0=sum2/sum3
       do k=kmin_u(i,j),int(kmergedr4_u(i,j))-1
        vel_u(i,j,k,2)=x0*log((h_u(i,j)+depth_u(i,j,k))/(0.5*(z0_w(i,j)+z0_w(i-1,j))))
       enddo
       k=int(kmergedr4_u(i,j))
                   vel_u(i,j,k,2)= &
            (kmergedr4_u(i,j)-int(kmergedr4_u(i,j))) *x0*log((h_u(i,j)+depth_u(i,j,k))/(0.5*(z0_w(i,j)+z0_w(i-1,j)))) &
       +(1.-(kmergedr4_u(i,j)-int(kmergedr4_u(i,j))))*vel_u(i,j,k,2)
      enddo         ; enddo

! vel_v:
      do j=2,jmax   ; do i=2,imax-1
       sum2=0. ; sum3=0.
       do k=kmin_v(i,j),int(kmergedr4_v(i,j))-1
        sum2=sum2+dz_v(i,j,k,dafter)*vel_v(i,j,k,2)
        sum3=sum3+dz_v(i,j,k,dafter)*log((h_v(i,j)+depth_v(i,j,k))/(0.5*(z0_w(i,j)+z0_w(i,j-1))))
       enddo
       k=int(kmergedr4_v(i,j))
       sum2=sum2+(kmergedr4_v(i,j)-int(kmergedr4_v(i,j)))*dz_v(i,j,k,dafter)*vel_v(i,j,k,2)
       sum3=sum3+(kmergedr4_v(i,j)-int(kmergedr4_v(i,j)))*dz_v(i,j,k,dafter)*log((h_v(i,j)+depth_v(i,j,k))/(0.5*(z0_w(i,j)+z0_w(i,j-1))))
! x0 est l'amplitude pour u=x0*log((h+z)/z0)
       x0=sum2/sum3
       do k=kmin_v(i,j),int(kmergedr4_v(i,j))-1
         vel_v(i,j,k,2)=x0*log((h_v(i,j)+depth_v(i,j,k))/(0.5*(z0_w(i,j)+z0_w(i,j-1))))
       enddo
       k=int(kmergedr4_v(i,j))
                   vel_v(i,j,k,2)= &
            (kmergedr4_v(i,j)-int(kmergedr4_v(i,j))) *x0*log((h_v(i,j)+depth_v(i,j,k))/(0.5*(z0_w(i,j)+z0_w(i,j-1)))) &
       +(1.-(kmergedr4_v(i,j)-int(kmergedr4_v(i,j))))*vel_v(i,j,k,2)
      enddo         ; enddo

      end subroutine momemtum_logprofile

!.......................................................................

      subroutine momentum_uv_minmax(t_) !16-01-22
      use module_principal ; use module_parallele
      implicit none
      integer t_

! Partie U

       x1=9999. ; x2=-x1 ; x3=9999. ; x4=-x1 
       do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax

!      if(mask_u(i,j,k)==1) then

        if(vel_u(i,j,k,t_)<x1) then
         x1=vel_u(i,j,k,t_) ; i1=i ; j1=j ; k1=k
        endif
        if(vel_u(i,j,k,t_)>x2) then
         x2=vel_u(i,j,k,t_) ; i2=i ; j2=j ; k2=k
        endif

!      endif

       enddo ; enddo ; enddo


      call mpi_allreduce(x1,x0,1,mpi_double_precision,       & !03/04/09
           mpi_min,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MIN U=',x0
      if(x1==x0) then !>>>
           write(6,*)'MIN U val',x0,vel_u(i1,j1,k1,t_)
           write(6,*)'MIN U rank i j k',par%rank,i1,j1,k1
           write(6,*)'MIN U glob i j k',i1+par%timax(1),j1+par%tjmax(1),k1
           write(6,*)'MIN U dz_u',dz_u(i1,j1,k1,t_)
           write(6,*)'MIN U h_u',h_u(i1,j1)
           write(6,*)'MIN U kmin',kmin_u(i1,j1)
           write(6,*)'MIN U kmerged_v',kmerged_u(i1,j1)
      endif           !>>>

      call mpi_allreduce(x2,x0,1,mpi_double_precision,       & !03/04/09
           mpi_max,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MAX U=',x0
      if(x2==x0) then !>>>
           write(6,*)'MAX U val',x0,vel_u(i2,j2,k2,t_)
           write(6,*)'MAX U rank i j k',par%rank,i2,j2,k2
           write(6,*)'MAX U glob i j k',i2+par%timax(1),j2+par%tjmax(1),k2
           write(6,*)'MAX U dz_u',dz_u(i2,j2,k2,t_)
           write(6,*)'MAX U h_u',h_u(i2,j2)
           write(6,*)'MAX U kmin',kmin_u(i2,j2)
           write(6,*)'MAX U kmerged_v',kmerged_u(i2,j2)
      endif           !>>>


! Partie V

       x1=9999. ; x2=-x1 ; x3=9999. ; x4=-x1 
       do k=1,kmax ; do j=2,jmax ; do i=2,imax-1
!      if(mask_v(i,j,k)==1) then

        if(vel_v(i,j,k,t_)<x1) then
         x1=vel_v(i,j,k,t_) ; i1=i ; j1=j ; k1=k
        endif
        if(vel_v(i,j,k,t_)>x2) then
         x2=vel_v(i,j,k,t_) ; i2=i ; j2=j ; k2=k
        endif

!      endif
       enddo ; enddo ; enddo


      call mpi_allreduce(x1,x0,1,mpi_double_precision,       & !03/04/09
           mpi_min,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MIN V=',x0
      if(x1==x0) then !>>>
           write(6,*)'MIN V val',x0,vel_v(i1,j1,k1,t_)
           write(6,*)'MIN V rank i j k',par%rank,i1,j1,k1
           write(6,*)'MIN V glob i j k',i1+par%timax(1),j1+par%tjmax(1),k1
           write(6,*)'MIN V dz_v',dz_v(i1,j1,k1,t_)
           write(6,*)'MIN V h_v',h_v(i1,j1)
           write(6,*)'MIN V kmin',kmin_v(i1,j1)
           write(6,*)'MIN V kmerged_v',kmerged_v(i1,j1)
      endif           !>>>

      call mpi_allreduce(x2,x0,1,mpi_double_precision,       & !03/04/09
           mpi_max,par%comm2d ,ierr)
!     if(par%rank==0)write(6,*)'MAX V=',x0
      if(x2==x0) then !>>>
           write(6,*)'MAX V val',x0,vel_v(i2,j2,k2,t_)
           write(6,*)'MAX V rank i j k',par%rank,i2,j2,k2
           write(6,*)'MAX V glob i j k',i2+par%timax(1),j2+par%tjmax(1),k2
           write(6,*)'MAX V dz_v',dz_v(i2,j2,k2,t_)
           write(6,*)'MAX V h_v',h_v(i2,j2)
           write(6,*)'MAX V kmin',kmin_v(i2,j2)
           write(6,*)'MAX V kmerged_v',kmerged_v(i2,j2)
      endif           !>>>


      end subroutine momentum_uv_minmax

!.......................................................................

      subroutine momentum_exp_adv_u_prime !20-10-21
      use module_principal ; use module_parallele
      implicit none
      real*4 :: diflxfactor_
      double precision dti_lpsub_
! Note: anyv3d(:,:,:,0) pris pour distribution verticale du stress de surface
!       anyv3d(:,:,:,3) pris pour nombre de courant vertical
      integer :: id_cnu_=1       &   ! current number relative to u velocity identifier
                ,id_cnv_=2       &   ! current number relative to v velocity identifier
                ,id_veldydz_u_=3 &   ! veldydz_u identifier
                ,id_veldxdz_v_=4 &   ! veldxdz_v identifier
                ,id_ubef_=0      &   ! u velocity identifier
                ,id_ubef2_=9     &   ! u velocity identifier before update
!               ,id_udiv=5       &   ! A definir
                ,looplimit_=1    &   ! nombre d'iteration maximum pour l'advection verticale
                ,id_we_=1        &   ! id_cnu_ n'est plus occupe quand id_we_  entre en service
                ,id_wi_=6        &   
                ,id_cnw_=2       &   ! id_cnv_ n'est plus occupe quand id_cnw_ entre en service
                ,loop_,loopmax_  &   ! u est entirerement terminEe quand commence l'etape v et donc
                ,loopmax_hor_u_      

      diflxfactor_=cst_adv_vel/16. 
! note: dans le mode externe on divise par 32 car l'advection y est traitee dans les 2 directions simultanement

! Cas particulier du modele 1DV
      if(flag_1dv==1) then !>>>>> !14-02-16
       vel_u(:,:,:,2)=vel_u(:,:,:,0)
       return
      endif                !>>>>>

!..............................................
! vel_u:
      do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+2
       anyv3d(i,j,k,id_ubef_)=vel_u(i,j,k,before)
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax  
       anyv3d(i,j,k,id_udiv)=0.
      enddo ; enddo ; enddo

!.............................................
! flux facettes
      do k=1,kmax ; do j=2,jmax-1 ; do i=1,imax
       anyv3d(i,j,k,id_veldydz_u_)=0.5*(veldydz_u(i,j,k,1)+veldydz_u(i+1,j,k,1)) 
      enddo       ; enddo         ; enddo

      do k=1,kmax ; do j=2,jmax ; do i=2,imax
       anyv3d(i,j,k,id_veldxdz_v_)=0.5*(veldxdz_v(i-1,j,k,1)+veldxdz_v(i,j,k,1))
      enddo       ; enddo         ; enddo

!..............................................
! Determiner le nombre de sous-iterations d'advection
! Cet algo suppose que la couche merged est 100% melangee avec les couches sous-jacente
! et donc que l'on a affaire A une couche dont l'epaisseur peut etre consideree comme
! z_w(kmerged+1)-z_w(1) et non pas dz(kmerged) (methode 2). Si on garde dz(kmerged) (methode 1) l'algo considere
! un volume plus petit pour plus de securite (mais evntuellement un dt inutilement trop petit)
! D'autre part les flux qui entrent dans ce volume sont ceux du niveau k=kmerged mais egalement
! ceux sous kmerged, kundermin_u etant le niveau du premier flux lateral non nul
      x1=0.
! Niveaux standard:
      do k=1,kmax !08-11-21
      do j=2,jmax-1 ; do i=2,imax
!       do k=kmerged_u(i,j)+1,kmax
         x1=max(x1,max(max(max(abs(anyv3d(i  ,j  ,k,id_veldydz_u_))  &
                              ,abs(anyv3d(i-1,j  ,k,id_veldydz_u_))) &
                              ,abs(anyv3d(i  ,j  ,k,id_veldxdz_v_))) &
                              ,abs(anyv3d(i  ,j+1,k,id_veldxdz_v_))) & 
                    *dti_lp                                          &
                    *wetmask_u(i,j)                                  &
                      /(dxdy_u(i,j)*min(dz_u(i,j,k,dbefore),dz_u(i,j,k,dafter))))
      enddo       ; enddo
      enddo

      call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_max,par%comm2d,ierr)
!     loopmax_=min(int(x2)+1,100) ! grande valeur expres pour ne pas borner....
      loopmax_=    int(x2)+1      
      dti_lpsub_=dti_lp/loopmax_
      loopmax_hor_u_=loopmax_ ! faire une copie de loopmax_ pour plus tard 

      if(par%rank==0)then
       open(unit=3,file=trim(tmpdirname)//'dti_u',position='append') !17-07-20
        write(3,*)real(elapsedtime_now/86400.),loopmax_,real(x2),real(dti_lpsub_)
       close(3)
      endif

      if(loopmax_>50.or.loopmax_<0) then !>>>
        call momentum_replay_ncumax(id_veldydz_u_,id_veldxdz_v_) !18-01-19
        call graph_out
        stop 'loopmax_ u advection out of range' !15-05-17
      endif                              !>>>

!.............................................
! flux explictes limitEs:
      if(loopmax_>looplimit_hor) then !-limiter-les-flux-explicites->

      do k=1,kmax ; do j=2,jmax-1 ; do i=1,imax

       anyv3d(i,j,k,id_veldydz_u_)=max(min( &

       0.5*(veldydz_u(i,j,k,1)+veldydz_u(i+1,j,k,1))  &

        ,+looplimit_hor*min(dxdy_u(i  ,j)*min(dz_u(i  ,j,k,dbefore),dz_u(i  ,j,k,dafter))/max(wetmask_u(i  ,j),1.e-10) &
                           ,dxdy_u(i+1,j)*min(dz_u(i+1,j,k,dbefore),dz_u(i+1,j,k,dafter))/max(wetmask_u(i+1,j),1.e-10) &
                           )/dti_lp)                                   &

        ,-looplimit_hor*min(dxdy_u(i  ,j)*min(dz_u(i  ,j,k,dbefore),dz_u(i  ,j,k,dafter))/max(wetmask_u(i  ,j),1.e-10) &
                           ,dxdy_u(i+1,j)*min(dz_u(i+1,j,k,dbefore),dz_u(i+1,j,k,dafter))/max(wetmask_u(i+1,j),1.e-10) &
                           )/dti_lp)                                    

      enddo       ; enddo         ; enddo

      do k=1,kmax ; do j=2,jmax ; do i=2,imax

       anyv3d(i,j,k,id_veldxdz_v_)=max(min( &

       0.5*(veldxdz_v(i-1,j,k,1)+veldxdz_v(i,j,k,1)) &

        ,+looplimit_hor*min(dxdy_u(i,j  )*min(dz_u(i,j  ,k,dbefore),dz_u(i,j  ,k,dafter))/max(wetmask_u(i,j  ),1.e-10) &
                           ,dxdy_u(i,j-1)*min(dz_u(i,j-1,k,dbefore),dz_u(i,j-1,k,dafter))/max(wetmask_u(i,j-1),1.e-10) &
                           )/dti_lp)                                   &

        ,+looplimit_hor*min(dxdy_u(i,j  )*min(dz_u(i,j  ,k,dbefore),dz_u(i,j  ,k,dafter))/max(wetmask_u(i,j  ),1.e-10) &
                           ,dxdy_u(i,j-1)*min(dz_u(i,j-1,k,dbefore),dz_u(i,j-1,k,dafter))/max(wetmask_u(i,j-1),1.e-10) &
                           )/dti_lp)                                    

      enddo       ; enddo         ; enddo


      endif                           !-limiter-les-flux-explicites->

!..............................................
! Calculer les nombres de courant correspondant
! NOTE: je considere dti_lp parce que la partie centree de l'advection
! est "figee" sur vel_u et non par sur anyv3d(id_ubef)
      do k=1,kmax
      do j=2,jmax-1 ; do i=1,imax

! Methode 1 lineaire:
!                anyv3d(i,j,k,id_cnu_)=                 &
!            abs(anyv3d(i,j,k,id_veldydz_u_))*dti_lp    &
!               /min(dxdy_u(i+1,j  )*dz_u(i+1,j,k,dbefore)    &
!                   ,dxdy_u(i  ,j  )*dz_u(i  ,j,k,dbefore))

! Methode 2 maintenue centree jusqu'a 0.9 (ou autre):
!                anyv3d(i,j,k,id_cnu_)=min(1.,max(0.,   &
!           (abs(anyv3d(i,j,k,id_veldydz_u_))*dti_lp    &
!               /min(dxdy_u(i+1,j  )*dz_u(i+1,j,k,dbefore)    & !22-08-16
!                   ,dxdy_u(i  ,j  )*dz_u(i  ,j,k,dbefore))   &
!               -0.9)*10.))

! Methode 3 : (j'ai ajoute abs pour pouvoir enlever **2)
!                anyv3d(i,j,k,id_cnu_)=min(1.,          & !25-11-16
!           (abs(anyv3d(i,j,k,id_veldydz_u_))*dti_lp           &
!               /min(dxdy_u(i+1,j  )*dz_u(i+1,j,k,dbefore)     &
!                   ,dxdy_u(i  ,j  )*dz_u(i  ,j,k,dbefore))) ) & 
!                *wetmask_u(i+1,j  )                           & !05-03-22   
!                *wetmask_u(i  ,j  )

! Methode 4: (nu-1)*upwzone + 1 , =nu si upwzone=1 et =1 si upwzone=0 !05-04-22
      anyv3d(i,j,k,id_cnu_)=min(1., & !min>

      (      & !(nu-1)
       (abs(anyv3d(i,j,k,id_veldydz_u_))*dti_lp          &
           /min(dxdy_u(i+1,j  )*dz_u(i+1,j,k,dbefore)    &
               ,dxdy_u(i  ,j  )*dz_u(i  ,j,k,dbefore)))  &
        -1.) & !(nu-1)
            *wetmask_u(i+1,j)                            &
            *wetmask_u(i  ,j)                            &
         *(1.-sponge_u(i+1,j,0))                         &
         *(1.-sponge_u(i  ,j,0))                         &

         +1.                                               &

                               )      !min>   

      enddo ; enddo
      enddo

      do k=1,kmax
      do j=2,jmax ; do i=2,imax


! Methode 1 lineaire:
!                anyv3d(i,j,k,id_cnv_)=                 &
!            abs(anyv3d(i,j,k,id_veldxdz_v_))*dti_lp    &
!               /min(dxdy_u(i,j-1)*dz_u(i,j-1,k,dbefore)      &
!                   ,dxdy_u(i,j  )*dz_u(i,j  ,k,dbefore))

! Methode 2 maintenue centree jusqu'a 0.9 (ou autre):
!                anyv3d(i,j,k,id_cnv_)=min(1.,max(0.,   &
!           (abs(anyv3d(i,j,k,id_veldxdz_v_))*dti_lp    &
!               /min(dxdy_u(i,j-1)*dz_u(i,j-1,k,dbefore)      &
!                   ,dxdy_u(i,j  )*dz_u(i,j  ,k,dbefore))     &
!               -0.9)*10.))

! Methode 3 : (j'ai ajoute abs pour pouvoir enlever **2)
!                anyv3d(i,j,k,id_cnv_)=min(1.,               &
!           (abs(anyv3d(i,j,k,id_veldxdz_v_))*dti_lp         &
!               /min(dxdy_u(i,j-1)*dz_u(i,j-1,k,dbefore)     &
!                   ,dxdy_u(i,j  )*dz_u(i,j  ,k,dbefore))) ) &    
!                *wetmask_u(i,j-1)                           & !05-03-22  
!                *wetmask_u(i,j  )


! Methode 4 : (nv-1)*upwzone + 1 , =nv si upwzone=1 et =1 si upwzone=0 !05-04-22
      anyv3d(i,j,k,id_cnv_)=min(1., & !min>

      (      & !(nu-1)
       (abs(anyv3d(i,j,k,id_veldxdz_v_))*dti_lp         &
           /min(dxdy_u(i,j-1)*dz_u(i,j-1,k,dbefore)     &
               ,dxdy_u(i,j  )*dz_u(i,j  ,k,dbefore)))   &    
        -1.) & !(nu-1)
            *wetmask_u(i,j-1)                           & 
            *wetmask_u(i,j  )                           &
         *(1.-sponge_u(i,j-1,0))                        & 
         *(1.-sponge_u(i,j  ,0))                        &

         +1.                                            &

                               )      !min>   

      enddo ; enddo
      enddo

!     if(mod(iteration3d,100)==0)write(6,*)'bidouille advection u'
!                anyv3d(:,:,:,id_cnu_)=1.
!                anyv3d(:,:,:,id_cnv_)=1.


!     i=imax/2 ; j=jmax/2 ; k=kmax/2
!     write(6,*)'--------------------------------------------'
!     write(6,*)'anyv3d(i,j,k,id_cnu_)',anyv3d(i,j,k,id_cnu_)
!     write(6,*)'anyv3d(i,j,k,id_cnv_)',anyv3d(i,j,k,id_cnv_)
!     write(6,*)'--------------------------------------------'
!     write(6,*)'loopmax_',loopmax_
!     write(6,*)'--------------------------------------------'
!     write(6,*)'anyv3d(i,j,k,id_veldydz_u_)',anyv3d(i,j,k,id_veldydz_u_)
!     write(6,*)'anyv3d(i,j,k,id_veldxdz_v_)',anyv3d(i,j,k,id_veldxdz_v_)
!     write(6,*)'--------------------------------------------'

      do loop_=1,loopmax_ !*********>

      do k=1,kmax   

! Flux direction Oi
      do j=2,jmax-1 ; do i=1,imax

      xflux_t(i,j)=                                                    &
      
              (1.- anyv3d(i,j,k,id_cnu_))*( &     !cccccc> ! adv scheme selector

! Ox Advective flux at t point:
                     -anyv3d(i,j,k,id_veldydz_u_)                      &
                   *( (vel_u(i  ,j,k,now)+  vel_u(i+1,j,k,now))*0.5625 &
                     -(vel_u(i-1,j,k,now)+  vel_u(i+2,j,k,now))*0.0625)&
! Ox Diffusive flux:
! Details sur le codage du bi-harmonique dans S:
! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit
!.......................................................................................
! Coefficient constant (dans l'espace des indices):
!     +cst_adv_vel*(dxdy_u(i,j)*dz_u(i,j,k,dnow)+dxdy_u(i+1,j)*dz_u(i+1,j,k,dnow))/(32*dti_lp) &
! Coefficient "smago like": !20-11-16
!     +(dy_u(i,j)*dz_u(i,j,k,dnow)+dy_u(i+1,j)*dz_u(i+1,j,k,dnow))    &
!     *abs(anyv3d(i+1,j,k,id_ubef_)-anyv3d(i,j,k,id_ubef_))*cst_adv_vel/64. &
! shema 3 hybride :
!     +(dy_u(i,j)*dz_u(i,j,k,dnow)+dy_u(i+1,j)*dz_u(i+1,j,k,dnow))*0.5*diflxfactor_*( &
!!    +(dy_t(i,j)*dz_t(i,j,k,dnow)                               )    *diflxfactor_*( & !voir note1 !19-07-19
!!        biharm_c2*0.5*abs(anyv3d(i+1,j,k,id_ubef_)-anyv3d(i,j,k,id_ubef_))  & 
!!       +biharm_c1*dx_t(i,j)*inv_dti_lp                                  ) & !09-06-18
!     +dy_t(i,j)*dz_t(i,j,k,dnow)*diflxfactor_*0.5*abs(anyv3d(i+1,j,k,id_ubef_)-anyv3d(i,j,k,id_ubef_)) & !28-06-22 
!     +diflxfactor_*0.5*abs(veldydz_u(i,j,k,1)-veldydz_u(i+1,j,k,1)) & !17-01-22 
      +diflxfactor_*0.5*abs(veldydz_u(i,j,k,id_bihar_lim)-veldydz_u(i+1,j,k,id_bihar_lim)) & !01-01-23
!.......................................................................................
       *( & !ooooo>
         (anyv3d(i+1,j,k,id_ubef_)-anyv3d(i  ,j,k,id_ubef_))*3.        &
        -(anyv3d(i+2,j,k,id_ubef_)-anyv3d(i-1,j,k,id_ubef_))           &
         *mask_u(i+2,j,kmax)      *mask_u(i-1,j,kmax)                  &
        ) & !ooooo>
                                          ) &     !cccccc>

                                   +anyv3d(i,j,k,id_cnu_)*( & !upupup>

        -0.5*( & !xxxx>
                    anyv3d(i,j,k,id_veldydz_u_)                      &
              *(    anyv3d(i,j,k,id_ubef_)+anyv3d(i+1,j,k,id_ubef_)) &
               +abs(anyv3d(i,j,k,id_veldydz_u_))                     & 
                  *(anyv3d(i,j,k,id_ubef_)-anyv3d(i+1,j,k,id_ubef_)) &
             ) & !xxxx>
                                                          )   !upupup>

      enddo ; enddo
! note1: dz_t plutot que demi_somme dz_u car dans masque dz_u peut etre n'importe quoi

! Advection partielle direction Oi:
      do j=2,jmax-1 ; do i=2,imax  


! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_udiv)=     &
       anyv3d(i,j,k,id_udiv)      &  
      +anyv3d(i,j,k,id_ubef_)*dti_lpsub_*( anyv3d(i  ,j,k,id_veldydz_u_) &
                                          -anyv3d(i-1,j,k,id_veldydz_u_))/dxdy_u(i,j)

! Couches fusionnees: memoriser champs avant update pour plus tard ne moyenner que la contrib advective
      !anyv3d(i,j,k,id_ubef2_)=anyv3d(i,j,k,id_ubef_) !04-12-19

      anyv3d(i,j,k,id_ubef_)=mask_u(i,j,k)*( &  !ooo> !25-04-17
      anyv3d(i,j,k,id_ubef_)*dz_u(i,j,k,dbefore)          & !01-05-19

       +dti_lpsub_*wetmask_u(i,j)*(                       & !--->

         (                                               & !pmx>
                         xflux_t(i  ,j)                  &
                        -xflux_t(i-1,j)                  &

      +anyv3d(i,j,k,id_ubef_)*( & 
                          anyv3d(i  ,j,k,id_veldydz_u_)  &
                         -anyv3d(i-1,j,k,id_veldydz_u_)  &
                        ) &
                                           )/dxdy_u(i,j) & !pmx>

                                 )                        & !--->

                                            ) & !ooo>      
                                             /dz_u(i,j,k,dbefore) & !04-11-21
      +(1-mask_u(i,j,k))*vel_u(i,j,k,before)!*dz_u(i,j,k,dbefore)

      enddo ; enddo

      enddo ! boucle k

! Sous le fond
     !call vertmix_merged_levels_u(dbefore,id_ubef_,id_ubef2_) !04-12-19

! ICI CONDITION MPI U1 U2 anyv3d(i,j,k,id_ubef_)
      call obc_mpi_anyv3d(1,id_ubef_,'u1')
      call obc_mpi_anyv3d(1,id_ubef_,'u2')

      do k=1,kmax   

! Flux direction Oj
      do j=2,jmax ; do i=2,imax

      yflux_f(i,j)=                                                    &

              (1.- anyv3d(i,j,k,id_cnv_))*( &     !cccccc> ! adv scheme selector

! Oy Advective flux at p point:
                                      -anyv3d(i,j,k,id_veldxdz_v_)     &
                *(  (vel_u(i  ,j-1,k,now)+  vel_u(i,j  ,k,now))*0.5625 &
                   -(vel_u(i  ,j-2,k,now)+  vel_u(i,j+1,k,now))*0.0625)&
! Oy Diffusive flux:
! Details sur le codage du bi-harmonique dans S:
! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit
!........................................................................................
! Coefficient constant (dans l'espace des indices):
!     +cst_adv_vel*(dxdy_u(i,j)*dz_u(i,j,k,dnow)+dxdy_u(i,j-1)*dz_u(i,j-1,k,dnow))/(32*dti_lp) &
! Coefficient "smago like"
!     +(dy_u(i,j)*dz_u(i,j,k,dnow)+dy_u(i,j-1)*dz_u(i,j-1,k,dnow))    &
!     *abs(anyv3d(i,j,k,id_ubef_)-anyv3d(i,j-1,k,id_ubef_))*cst_adv_vel/64. &
! Schema 3 hybride:
!     +(dy_u(i,j)*dz_u(i,j,k,dnow)+dy_u(i,j-1)*dz_u(i,j-1,k,dnow))*0.5*diflxfactor_*( &
!!     +mask_f(i,j,k)*(dy_u(i,j)*dz_u(i,j,k,dnow)+dy_u(i,j-1)*dz_u(i,j-1,k,dnow))*0.5*diflxfactor_*( & !note2 !19-07-19
!!         biharm_c2*0.5*abs(anyv3d(i,j,k,id_ubef_)-anyv3d(i,j-1,k,id_ubef_))   &
!!        +biharm_c1*dy_f(i,j)*inv_dti_lp                                 ) &
!      +mask_f(i,j,k)*diflxfactor_*0.5*abs(veldxdz_v(i-1,j,k,1)-veldxdz_v(i,j,k,1)) & !17-01-22
       +mask_f(i,j,k)*diflxfactor_*0.5*abs(veldxdz_v(i-1,j,k,id_bihar_lim)-veldxdz_v(i,j,k,id_bihar_lim)) & !01-01-23
!!!                  +diflxfactor_*0.5*abs(veldydz_u(i,j-1,k,1)-veldydz_u(i,j,k,1)) & !11-06-22
!      +0.5*(dy_u(i,j  )*dz_u(i,j  ,k,dnow)*mask_u(i,j  ,kmax)  & !28-06-22
!           +dy_u(i,j-1)*dz_u(i,j-1,k,dnow)*mask_u(i,j-1,kmax)) &
!      *diflxfactor_*0.5*abs(anyv3d(i,j,k,id_ubef_)-anyv3d(i,j-1,k,id_ubef_))   &
!      *(1+flag_free_slip*(mask_f(i,j,k)-1)) & ! free slip si flag_free_slip=1        !11-06-22 
!........................................................................................
       *( & !ooooo> 
         (anyv3d(i  ,j  ,k,id_ubef_)-anyv3d(i  ,j-1,k,id_ubef_))*3. &
        -(anyv3d(i  ,j+1,k,id_ubef_)-anyv3d(i  ,j-2,k,id_ubef_))    &
         *mask_u(i  ,j+1,kmax)      *mask_u(i  ,j-2,kmax)           & !25-11-16
        ) & !ooooo>
                                          ) &     !cccccc>

                                   +anyv3d(i,j,k,id_cnv_)*( & !upupup>

        -0.5*( & !yyyy>
                anyv3d(i,j,k,id_veldxdz_v_)                      &
!!!!!!!        0.5*(veldxdz_v(i-1,j,k,1)+veldxdz_v(i,j,k,1))     &
              *(anyv3d(i,j-1,k,id_ubef_)+anyv3d(i,j,k,id_ubef_)) &
           +abs(anyv3d(i,j,k,id_veldxdz_v_))                     &
              *(anyv3d(i,j-1,k,id_ubef_)-anyv3d(i,j,k,id_ubef_)) &
             ) & !yyyy>
                                                          )   !upupup>

      enddo ; enddo

! note2: le flux advectif yflux_f(i,j) est nul si mask_t(i-1,j)=mask_t(i,j)=0 (car alors veldxdz_v(i-1:i,j)=0)
! flux diffusif nul (si *mask_f(i,j,k)) car dz_u peut etre n'importe quoi dans le masque et compromettre CFL
! Par consequent cela correspond A une condition "free-slip" A la cOte.


! Advection partielle direction Oj:
      do j=2,jmax-1 ; do i=2,imax  

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_udiv)=     &
       anyv3d(i,j,k,id_udiv)      &  
      +anyv3d(i,j,k,id_ubef_)*dti_lpsub_*( anyv3d(i,j+1,k,id_veldxdz_v_) &
                                          -anyv3d(i,j  ,k,id_veldxdz_v_))/dxdy_u(i,j)

! Couches fusionnees: memoriser champs avant update pour plus tard ne moyenner que la contrib advective
      !anyv3d(i,j,k,id_ubef2_)=anyv3d(i,j,k,id_ubef_) !04-12-19

      anyv3d(i,j,k,id_ubef_)=mask_u(i,j,k)*( & !ooo>
      anyv3d(i,j,k,id_ubef_)*dz_u(i,j,k,dbefore) & !01-05-19

       +dti_lpsub_*wetmask_u(i,j)*(           & !--->

          (                                  & !pmx>
           +yflux_f(i,j+1)-yflux_f(i  ,j)    &

      +anyv3d(i,j,k,id_ubef_)*( & 
          anyv3d(i,j+1,k,id_veldxdz_v_)      &
         -anyv3d(i,j  ,k,id_veldxdz_v_)      &
                        ) &
           )/dxdy_u(i,j)                     & !pmx>

                   )                          & !--->
                                          )  & !ooo>
                                             /dz_u(i,j,k,dbefore) & !04-11-21
      +(1-mask_u(i,j,k))*vel_u(i,j,k,before)!*dz_u(i,j,k,dbefore)

      enddo ; enddo

      enddo ! k loop

! Sous le fond
     !call vertmix_merged_levels_u(dbefore,id_ubef_,id_ubef2_) !04-12-19

! ICI CONDITION MPI U1 U2 anyv3d(i,j,k,id_ubef_)
      call obc_mpi_anyv3d(1,id_ubef_,'u1')
      call obc_mpi_anyv3d(1,id_ubef_,'u2')

      enddo               !*********> ! loopmax_


      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax  
        vel_u(i,j,k,2)=anyv3d(i,j,k,id_ubef_) 
      enddo ; enddo ; enddo

! Partie implicite de l'advection horizontale:
      if(loopmax_hor_u_>looplimit_hor)call momentum_adv_u_impli(id_veldydz_u_,id_veldxdz_v_)

! ---- ADVECTION VERTICALE ----

! Advection verticale explicite:

!     x1=0. ! Pas obligatoire seulement pour diag du maximum
      xy_u(:,:,1)=small2 ! reset small2 evite que le nombre d'iteration soit nulle
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax  
! x0=Dz/Dt
       x0=min(dz_u(i,j,k,dbefore),dz_u(i,j,k-1,dbefore))/dti_lp
! Omega bornee par + ou - N fois dz/dt ici 2 fois
       anyv3d(i,j,k,id_we_)=                                         &
            max(min(0.5*(omega_w(i,j,k,1)+omega_w(i-1,j,k,1))        &
               ,looplimit_*x0),-looplimit_*x0)*wetmask_u(i,j)        &
                                           *(1.-sponge_u(i,j,0)) !05-04-22
! max over z of the decimal current number:
       xy_u(i,j,1)=max(xy_u(i,j,1),abs(anyv3d(i,j,k,id_we_))/x0)
!      x1=max(x1,xy_u(i,j,1)) ! Pas obligatoire seulement pour diag du maximum

      enddo       ; enddo       ; enddo

! Fond et couche merged (dans laquelle w est 100% implicite, donc w explitite = 0)
      do j=2,jmax-1 ; do i=2,imax  
       do k=1,kmerged_u(i,j)
        anyv3d(i,j,k,id_we_)=0.
       enddo
      enddo       ; enddo

! Surface:
      do j=2,jmax-1 ; do i=2,imax  
!      anyv3d(i,j,kmax+1,id_we_)=0.5*(omega_w(i,j,kmax+1,1)+omega_w(i-1,j,kmax+1,1)) ! A la surface w explicite pour flux de surface
! Depuis !07-03-22 vitesse de surface 100% implicite combinee A condition vel(kmax+1)=0 dans triadia_in(:,:,kmax,3)
       anyv3d(i,j,kmax+1,id_we_)=0. !07-03-22
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
      do k=1,kmax+1 ; do j=2,jmax-1 ; do i=2,imax  
       anyv3d(i,j,k,id_wi_)=0.5*(omega_w(i,j,k,1)+omega_w(i-1,j,k,1))-anyv3d(i,j,k,id_we_)
      enddo       ; enddo       ; enddo


! Dans la couche merged, profil de vitesse verticale implicite de type "gradient de flux constant"
      do j=2,jmax-1 ; do i=2,imax  

       do k=2,kmerged_u(i,j)

! Noter qu'on estime depth_w_u (depth_w au point u) comme etant depth_u(k)-0.5*dz_u(k)
         anyv3d(i,j,k,id_wi_)=anyv3d(i,j,kmerged_u(i,j)+1,id_wi_)         &
      *(depth_u(i,j,k            )-0.5*dz_u(i,j,k,dnow)+h_u(i,j)) &
      /(depth_u(i,j,kmerged_u(i,j)+1)-0.5*dz_u(i,j,kmerged_u(i,j),dnow)+h_u(i,j))

       enddo

      enddo       ; enddo
      check3=anyv3d(imax/2,jmax/2,kmax/2,id_wi_)

      if(flag_timesplitting==1) then !AVEC time-splitting> !28-04-20
       i0=0 ; x1_r4=0.5
      else                           !SANS time-splitting>
       i0=1 ; x1_r4=0.25 
      endif                          !SANS time-splitting>

      do j=2,jmax-1 ; do i=2,imax  ! Boucle i,j No1

! ceiling de xy_u(i,j,1) nombre de sous-iterations advectives
       loopmax_=min( ceiling(xy_u(i,j,1)) , looplimit_ )
       dti_lpsub=dti_lp/loopmax_

!      do k=kmerged_u(i,j)+1,kmax                         !22-12-14
       do k=2,kmax                         !22-12-14

!      if(i+par%timax(1)==39.and.j+par%tjmax(1)==10) &
!      write(6,*)'w',k,anyv3d(i,j,k,id_we_),anyv3d(i,j,k,id_wi_)

! Current number with additional limititations
        anyv3d(i,j,k,id_cnw_)=1.-( &!ooooo>
! Robustness regarding rivers and dry zones:
                       wetmask_u(i,j)       & ! wet/dry zone
!      *0.5*(upwindriver_t(i,j)+upwindriver_t(i-1,j)) & ! upwind  zone !17-10-21
! Robustness regarding vertical current number:
        *max(1.-abs(anyv3d(i,j,k,id_we_)*dti_lpsub               &
                 /min(dz_u(i,j,k,dbefore),dz_u(i,j,k-1,dbefore))),0.)        &
                                ) !ooooo>

!        anyv3d(i,j,k,id_cnw_)=0. ! 100%C2
!        anyv3d(i,j,k,id_cnw_)=1. ! 100%UP2

       enddo

! Cas particuliers k=1:kmerged_u(i,j) et k=kmax:
!      anyv3d(i,j,1:kmerged_u(i,j),id_cnw_)=0.
       anyv3d(i,j,1            ,id_cnw_)=0.
       anyv3d(i,j,kmax+1       ,id_cnw_)=0. ! A la surface schema 100% explicite (voir ce qui est
                                            ! fait pour le sel

! Flux centres frozen:
       do k=2,kmax
        anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cnw_))*anyv3d(i,j,k,id_we_)  &
                  *x1_r4*(     vel_u(i,j,k  ,1)+vel_u(i,j,k-1,1)     &
                          +i0*(vel_u(i,j,k  ,0)+vel_u(i,j,k-1,0)) ) !28-04-20
       enddo
       anyv1d(1     ,1)=0.  
       anyv1d(1     ,3)=0.  
       anyv1d(kmax+1,3)=0. 
       k=kmax+1
       anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cnw_))*anyv3d(i,j,k,id_we_)*vel_u(i,j,k-1,1)

! Partie Iterative:
       do loop_=1,loopmax_ ! boucle iterative >>>

       vel_u(i,j,kmax+1,2)=vel_u(i,j,kmax,2)
       vel_u(i,j,0     ,2)=vel_u(i,j,1   ,2)

       do k=2,kmax !nexon>

        x1=0.5*(anyv3d(i,j,k,id_we_)+abs(anyv3d(i,j,k,id_we_)))
        x2=0.5*(anyv3d(i,j,k,id_we_)-abs(anyv3d(i,j,k,id_we_)))

! Flux upwind:
        anyv1d(k,3)=  &
          -anyv3d(i,j,k,id_cnw_) *(x1*vel_u(i,j,k-1,2)+x2*vel_u(i,j,k  ,2))&
       +(1-anyv3d(i,j,k,id_cnw_))*( &
               x1*(vel_u(i,j,k  ,2)-2.*vel_u(i,j,k-1,2)+vel_u(i,j,k-2,2))&!21-11-16
              +x2*(vel_u(i,j,k+1,2)-2.*vel_u(i,j,k  ,2)+vel_u(i,j,k-1,2))&
                                  )*(1+anyv3d(i,j,k,id_cnw_))/6. ! QUICKEST
!                                                           )/6. ! UBS-UPW


       enddo       !nexon>

       do k=1,kmax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_udiv)= & !26-01-17
       anyv3d(i,j,k,id_udiv)+dti_lpsub*(anyv3d(i,j,k+1,id_we_)    &
                                       -anyv3d(i,j,k  ,id_we_))*vel_u(i,j,k,2) 

         vel_u(i,j,k,2)=vel_u(i,j,k,2)+dti_lpsub*wetmask_u(i,j)*( &

                              anyv1d(k+1,3)+anyv1d(k+1,1) &
                             -anyv1d(k  ,3)-anyv1d(k  ,1) &

                 +vel_u(i,j,k,2)*( anyv3d(i,j,k+1,id_we_)    & !26-01-17
                                  -anyv3d(i,j,k  ,id_we_))   &

                                                 )/dz_u(i,j,k,dbefore)


       enddo

       enddo               ! boucle iterative >>>

      enddo       ; enddo          ! Boucle i,j No1

#ifdef bidon
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax  

           vel_u(i,j,k,2)                                &
          =vel_u(i,j,k,2)                                &

      -wetmask_u(i,j)*anyv3d(i,j,k,id_udiv)/dz_u(i,j,k,dbefore) !26-01-17

      enddo         ; enddo       ; enddo ! Boucle i,j No2
#endif
      check7=anyv3d(imax/2,jmax/2,kmax-1,id_udiv)

! Partie implicite de l'advection horizontale:
!     if(loopmax_hor_u_>looplimit_hor)call momentum_adv_u_impli(id_veldydz_u_,id_veldxdz_v_)

      end subroutine momentum_exp_adv_u_prime

!...................................................................

      subroutine momentum_exp_adv_v_prime !20-10-21
      use module_principal ; use module_parallele
      implicit none
      real*4 :: diflxfactor_
      double precision dti_lpsub_
! Note: anyv3d(:,:,:,0) pris pour distribution verticale du stress de surface
!       anyv3d(:,:,:,3) pris pour nombre de courant vertical
      integer :: id_cnu_=1       &   ! current number relative to u velocity identifier
                ,id_cnv_=2       &   ! current number relative to v velocity identifier
                ,id_veldydz_u_=3 &   ! veldydz_u identifier
                ,id_veldxdz_v_=4 &   ! veldxdz_v identifier
                ,id_vbef_=0      &   ! u velocity identifier
                ,id_vbef2_=9     &   ! u velocity identifier before update
!               ,id_vdiv=5       &   ! A definir
                ,looplimit_=1    &   ! nombre d'iteration maximum pour l'advection verticale
                ,id_we_=1        &   ! id_cnu_ n'est plus occupe quand id_we_  entre en service
                ,id_wi_=7        &   
                ,id_cnw_=2       &   ! id_cnv_ n'est plus occupe quand id_cnw_ entre en service
                ,loopmax_hor_v_  &   
                ,loop_,loopmax_      ! u est entirerement terminEe quand commence l'etape v et donc
                                     ! l'espace anyv3d(id_vbef_) est libre quand commence l'etape v

      diflxfactor_=cst_adv_vel/16.
! note: dans le mode externe on divise par 32 car l'advection y est traitee dans les 2 directions simultanement

! Cas particulier du modele 1DV
      if(flag_1dv==1) then !>>>>> !14-02-16
       vel_v(:,:,:,2)=vel_v(:,:,:,0)
       return
      endif                !>>>>>

!..............................................
! vel_v:
      do k=1,kmax ; do j=0,jmax+2 ; do i=0,imax+1
       anyv3d(i,j,k,id_vbef_)=vel_v(i,j,k,before)
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1  
       anyv3d(i,j,k,id_vdiv)=0.
      enddo ; enddo ; enddo
!..............................................
! flux facettes:
      do k=1,kmax ; do j=2,jmax ; do i=2,imax
       anyv3d(i,j,k,id_veldydz_u_)=0.5*(veldydz_u(i,j-1,k,1)+veldydz_u(i,j,k,1))      
      enddo       ; enddo       ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=2,imax-1
       anyv3d(i,j,k,id_veldxdz_v_)=0.5*(veldxdz_v(i,j,k,1)+veldxdz_v(i,j+1,k,1))        
      enddo       ; enddo       ; enddo

!..............................................
! Determiner le nombre de sous-iterations d'advection
! Cet algo suppose que la couche merged est 100% melangee avec les couches sous-jacente
! et donc que l'on a affaire A une couche dont l'epaisseur peut etre consideree comme
! z_w(kmerged+1)-z_w(1) et non pas dz(kmerged) (methode 2). Si on garde dz(kmerged) (methode 1) l'algo considere
! un volume plus petit pour plus de securite (mais evntuellement un dt inutilement trop petit)
! D'autre part les flux qui entrent dans ce volume sont ceux du niveau k=kmerged mais egalement
! ceux sous kmerged, kundermin_v etant le niveau du premier flux lateral non nul
      x1=0.
! Niveaux standard:
      do k=1,kmax !08-11-21
      do j=2,jmax ; do i=2,imax-1
!       do k=kmerged_v(i,j)+1,kmax
         x1=max(x1,max(max(max(abs(anyv3d(i  ,j  ,k,id_veldydz_u_))  &
                              ,abs(anyv3d(i+1,j  ,k,id_veldydz_u_))) &
                              ,abs(anyv3d(i  ,j  ,k,id_veldxdz_v_))) &
                              ,abs(anyv3d(i  ,j-1,k,id_veldxdz_v_))) & 
                    *dti_lp                                          &
                    *wetmask_v(i,j)                                  &
                      /(dxdy_v(i,j)*min(dz_v(i,j,k,dbefore),dz_v(i,j,k,dafter))))
      enddo       ; enddo
      enddo

      call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_max,par%comm2d,ierr)
!     loopmax_=min(int(x2)+1,100) ! grande valeur expres pour ne pas borner....
      loopmax_=    int(x2)+1      
      dti_lpsub_=dti_lp/loopmax_
      loopmax_hor_v_=loopmax_ ! faire une copie de loopmax_ pour plus tard 

      if(par%rank==0)then
       open(unit=3,file=trim(tmpdirname)//'dti_v',position='append') !17-07-20
        write(3,*)real(elapsedtime_now/86400.),loopmax_,real(x2),real(dti_lpsub_)
       close(3)
      endif

      if(loopmax_>50.or.loopmax_<0) then !>>>
        call momentum_replay_ncvmax(id_veldydz_u_,id_veldxdz_v_) !18-01-19
        call graph_out
        stop 'loopmax_ v advection out of range' !15-05-17
      endif                              !>>>

!.............................................
! flux explictes limitEs:
      if(loopmax_>looplimit_hor) then !-limiter-les-flux-explicites->

      do k=1,kmax ; do j=2,jmax ; do i=2,imax

       anyv3d(i,j,k,id_veldydz_u_)=max(min( &

       0.5*(veldydz_u(i,j,k,1)+veldydz_u(i,j-1,k,1))  &

        ,+looplimit_hor*min(dxdy_v(i  ,j)*min(dz_v(i  ,j,k,dbefore),dz_v(i  ,j,k,dafter))/max(wetmask_v(i  ,j),1.e-10) &
                           ,dxdy_v(i-1,j)*min(dz_v(i-1,j,k,dbefore),dz_v(i-1,j,k,dafter))/max(wetmask_v(i-1,j),1.e-10) &
                           )/dti_lp)                                   &

        ,-looplimit_hor*min(dxdy_v(i  ,j)*min(dz_v(i  ,j,k,dbefore),dz_v(i  ,j,k,dafter))/max(wetmask_v(i  ,j),1.e-10) &
                           ,dxdy_v(i-1,j)*min(dz_v(i-1,j,k,dbefore),dz_v(i-1,j,k,dafter))/max(wetmask_v(i-1,j),1.e-10) &
                           )/dti_lp)                                    

      enddo       ; enddo         ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=2,imax-1
       anyv3d(i,j,k,id_veldxdz_v_)=max(min( &

       0.5*(veldxdz_v(i,j+1,k,1)+veldxdz_v(i,j,k,1)) &

        ,+looplimit_hor*min(dxdy_v(i,j  )*min(dz_v(i,j  ,k,dbefore),dz_v(i,j  ,k,dafter))/max(wetmask_v(i,j  ),1.e-10) &
                           ,dxdy_v(i,j+1)*min(dz_v(i,j+1,k,dbefore),dz_v(i,j+1,k,dafter))/max(wetmask_v(i,j+1),1.e-10) &
                           )/dti_lp)                                   &

        ,-looplimit_hor*min(dxdy_v(i,j  )*min(dz_v(i,j  ,k,dbefore),dz_v(i,j  ,k,dafter))/max(wetmask_v(i,j  ),1.e-10) &
                           ,dxdy_v(i,j+1)*min(dz_v(i,j+1,k,dbefore),dz_v(i,j+1,k,dafter))/max(wetmask_v(i,j+1),1.e-10) &
                           )/dti_lp)                                    
      enddo       ; enddo         ; enddo

      endif                           !-limiter-les-flux-explicites->

!..............................................
! Calculer les nombres de courant correspondant
! NOTE: je considere dti_lp parce que la partie centree de l'advection
! est "figee" sur vel_v et non par sur anyv3d(id_ubef)
      do k=1,kmax
      do j=2,jmax ; do i=2,imax

! Methode 3 
!      anyv3d(i,j,k,id_cnu_)=min(1.,                &
!         (abs(anyv3d(i,j,k,id_veldydz_u_))*dti_lp  &
!      /min(dxdy_v(i-1,j)*dz_v(i-1,j,k,dbefore)     &
!          ,dxdy_v(  i,j)*dz_v(i  ,j,k,dbefore))) ) &   
!       *wetmask_v(i-1,j)                           & !05-03-22
!       *wetmask_v(i  ,j)

! Methode 4: (nu-1)*upwzone + 1 , =nu si upwzone=1 et =1 si upwzone=0 !05-04-22
       anyv3d(i,j,k,id_cnu_)=min(1., & !min>

      (      & !(nu-1)
          (abs(anyv3d(i,j,k,id_veldydz_u_))*dti_lp  &
       /min(dxdy_v(i-1,j)*dz_v(i-1,j,k,dbefore)     &
           ,dxdy_v(  i,j)*dz_v(i  ,j,k,dbefore)))   &   
        -1.) & !(nu-1)
            *wetmask_v(i-1,j)                       & 
            *wetmask_v(i  ,j)                       &
         *(1.-sponge_v(i-1,j,0))                    & 
         *(1.-sponge_v(i  ,j,0))                    &

         +1.                                        &

                                )      !min>   

      enddo ; enddo
      enddo

      do k=1,kmax
      do j=1,jmax ; do i=2,imax-1

! Methode 3 
!      anyv3d(i,j,k,id_cnv_)=min(1.,                &
!         (abs(anyv3d(i,j,k,id_veldxdz_v_))*dti_lp  &
!      /min(dxdy_v(i,j+1)*dz_v(i,j+1,k,dbefore)     &
!          ,dxdy_v(i,j  )*dz_v(i,j  ,k,dbefore))) ) &    
!       *wetmask_v(i,j+1)                           & !05-03-22
!       *wetmask_v(i,j  )

! Methode 4 : (nv-1)*upwzone + 1 , =nv si upwzone=1 et =1 si upwzone=0 !05-04-22
       anyv3d(i,j,k,id_cnv_)=min(1., & !min>

      (      & !(nv-1)
          (abs(anyv3d(i,j,k,id_veldxdz_v_))*dti_lp  &
       /min(dxdy_v(i,j+1)*dz_v(i,j+1,k,dbefore)     &
           ,dxdy_v(i,j  )*dz_v(i,j  ,k,dbefore)))   &    
        -1.) & !(nv-1)
            *wetmask_v(i,j+1)                       &
            *wetmask_v(i,j  )                       &
         *(1.-sponge_v(i,j+1,0))                    &
         *(1.-sponge_v(i,j  ,0))                    &

         +1.                                        &

                                )      !min>   

      enddo       ; enddo
      enddo


      do loop_=1,loopmax_ !*********>

      do k=1,kmax   

! Flux direction Oi
      do j=2,jmax ; do i=2,imax

      xflux_f(i,j)=                                                   &

               (1.-anyv3d(i,j,k,id_cnu_))*( &     !cccccc> ! adv scheme selector

! Ox Advective flux at p point:
                     -anyv3d(i,j,k,id_veldydz_u_)                      &
                 *(  (vel_v(i-1,j  ,k,now)+vel_v(i  ,j,k,now))*0.5625  &
                    -(vel_v(i-2,j  ,k,now)+vel_v(i+1,j,k,now))*0.0625) &
! Ox Diffusive flux at p point:
! Details sur le codage du bi-harmonique dans S:
! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit
!.........................................................................................
! Coef constant:
!     +cst_adv_vel*(dxdy_v(i,j)*dz_v(i,j,k,dnow)+dxdy_v(i-1,j)*dz_v(i-1,j,k,dnow))/(32*dti_lp) &
! Coef smago like:
!     +(dx_v(i,j)*dz_v(i,j,k,dnow)+dx_v(i-1,j)*dz_v(i-1,j,k,dnow))       &
!    *abs(anyv3d(i  ,j,k,id_vbef_ )-anyv3d(i-1,j,k,id_vbef_ ))*cst_adv_vel/64. &
! schema 3 hybride
!     +(dx_v(i,j)*dz_v(i,j,k,dnow)+dx_v(i-1,j)*dz_v(i-1,j,k,dnow))*0.5*diflxfactor_*( &
!!    +mask_f(i,j,k)*(dx_v(i,j)*dz_v(i,j,k,dnow)+dx_v(i-1,j)*dz_v(i-1,j,k,dnow))*0.5*diflxfactor_*( & !note3 !19-07-19
!!      biharm_c2*0.5*abs(anyv3d(i  ,j,k,id_vbef_ )-anyv3d(i-1,j,k,id_vbef_ ))  &
!!     +biharm_c1*dx_f(i,j)*inv_dti_lp                                    ) &
!     +mask_f(i,j,k)*diflxfactor_*0.5*abs(veldydz_u(i,j-1,k,1)-veldydz_u(i,j,k,1)) & !17-01-22
      +mask_f(i,j,k)*diflxfactor_*0.5*abs(veldydz_u(i,j-1,k,id_bihar_lim)-veldydz_u(i,j,k,id_bihar_lim)) & !01-01-23
!!!!                +diflxfactor_*0.5*abs(veldxdz_v(i-1,j,k,1)-veldxdz_v(i,j,k,1)) & !11-06-22
!      +0.5*(dx_v(i  ,j)*dz_v(i  ,j,k,dnow)*mask_v(i  ,j,kmax)  & !28-06-22
!           +dx_v(i-1,j)*dz_v(i-1,j,k,dnow)*mask_v(i-1,j,kmax)) &
!      *diflxfactor_*0.5*abs(anyv3d(i,j,k,id_vbef_)-anyv3d(i-1,j,k,id_vbef_))  &
!      *(1+flag_free_slip*(mask_f(i,j,k)-1)) & ! free slip si flag_free_slip=1       !11-06-22 
!.........................................................................................
       *( & !ooooo>
         (anyv3d(i  ,j,k,id_vbef_ )-anyv3d(i-1,j,k,id_vbef_ ))*3.      &
        -(anyv3d(i+1,j,k,id_vbef_ )-anyv3d(i-2,j,k,id_vbef_ ))         &
         *mask_v(i+1,j,kmax)       *mask_v(i-2,j,kmax)                 &
        ) & !ooooo>

                                          ) &     !cccccc>

                                   +anyv3d(i,j,k,id_cnu_)*( & !upupup>

         -0.5*( & !xxxx>
                 anyv3d(i,j,k,id_veldydz_u_)                       &
               *(anyv3d(i-1,j,k,id_vbef_)+anyv3d(i,j,k,id_vbef_))  &
            +abs(anyv3d(i,j,k,id_veldydz_u_))                      &
               *(anyv3d(i-1,j,k,id_vbef_)-anyv3d(i,j,k,id_vbef_))  &
                         ) & !xxxx>
                                                          )   !upupup>


      enddo ; enddo
! note3: flux diffusif nul (si *mask_f) car dz_v peut etre n'importe quoi dans le masque

! Advection partielle direction Oi:
      do j=2,jmax ; do i=2,imax-1

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_vdiv)=     &
       anyv3d(i,j,k,id_vdiv)      &  
      +anyv3d(i,j,k,id_vbef_)*dti_lpsub_*( anyv3d(i+1,j,k,id_veldydz_u_) &
                                          -anyv3d(i  ,j,k,id_veldydz_u_))/dxdy_v(i,j)

! Couches fusionnees: memoriser champs avant update pour plus tard ne moyenner que la contrib advective
      !anyv3d(i,j,k,id_vbef2_)=anyv3d(i,j,k,id_vbef_) !04-12-19

      anyv3d(i,j,k,id_vbef_)=mask_v(i,j,k)*( & !ooo>
      anyv3d(i,j,k,id_vbef_)*dz_v(i,j,k,dbefore)          & !01-05-19

       +dti_lpsub_*wetmask_v(i,j)*( & !--->

             (                                            & !pmx>
                         xflux_f(i+1,j)                   &
                        -xflux_f(i  ,j)                   &

      +anyv3d(i,j,k,id_vbef_)*( anyv3d(i+1,j,k,id_veldydz_u_)   & !16-01-17
                               -anyv3d(i  ,j,k,id_veldydz_u_) ) &

                                           )/dxdy_v(i,j)  & !pmx>

               ) & !--->
                                           ) & !ooo>
                                              /dz_v(i,j,k,dbefore) & !04-11-21
       +(1-mask_v(i,j,k))*vel_v(i,j,k,before)!*dz_v(i,j,k,dbefore)

      enddo ; enddo

      enddo ! boucle k

! Sous le fond
     !call vertmix_merged_levels_v(dbefore,id_vbef_,id_vbef2_) !04-12-19

! ICI CONDITION MPI v1 v2 anyv3d(i,j,k,id_vbef_)
      call obc_mpi_anyv3d(1,id_vbef_,'v1')
      call obc_mpi_anyv3d(1,id_vbef_,'v2')

      do k=1,kmax   

! Flux direction Oj
      do j=1,jmax ; do i=2,imax-1

      yflux_t(i,j)=                                                  &

               (1.-anyv3d(i,j,k,id_cnv_))*( &     !cccccc> ! adv scheme selector

! Oy Advective flux at t point:
                     -anyv3d(i,j,k,id_veldxdz_v_)                      &
                  *( (vel_v(i,j  ,k,now)+vel_v(i,j+1,k,now))*0.5625    &
                    -(vel_v(i,j-1,k,now)+vel_v(i,j+2,k,now))*0.0625)   &
! Oy Diffusive flux at t point:
! Details sur le codage du bi-harmonique dans S:
! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit
!.........................................................................................
! Coef constant:
!     +cst_adv_vel*(dxdy_v(i,j)*dz_v(i,j,k,dnow)+dxdy_v(i,j+1)*dz_v(i,j+1,k,dnow))/(32*dti_lp) &
! Coef smago like:
!     +(dx_v(i,j)*dz_v(i,j,k,dnow)+dx_v(i,j+1)*dz_v(i,j+1,k,dnow))       &
!    *abs(anyv3d(i,j+1,k,id_vbef_ )-anyv3d(i,j  ,k,id_vbef_ ))*cst_adv_vel/64. &
! schema 3 hybride
!     +(dx_v(i,j)*dz_v(i,j,k,dnow)+dx_v(i,j+1)*dz_v(i,j+1,k,dnow))*0.5*diflxfactor_*( &
!!    +(dx_t(i,j)*dz_t(i,j,k,dnow)                               )    *diflxfactor_*( & !note4 !19-07-19
!!      biharm_c2*0.5*abs(anyv3d(i,j+1,k,id_vbef_ )-anyv3d(i,j  ,k,id_vbef_ ))  &
!!     +biharm_c1*dy_t(i,j)*inv_dti_lp                                    ) &
!!!   +diflxfactor_*0.5*abs(veldxdz_v(i,j,k,1)-veldxdz_v(i,j+1,k,1)) & !17-01-22
!     +dx_t(i,j)*dz_t(i,j,k,dnow)*diflxfactor_*0.5*abs(anyv3d(i,j+1,k,id_vbef_)-anyv3d(i,j,k,id_vbef_)) & !28-06-22
!     +diflxfactor_*0.5*abs(veldxdz_v(i,j,k,1)-veldxdz_v(i,j+1,k,1)) & !17-01-22
      +diflxfactor_*0.5*abs(veldxdz_v(i,j,k,id_bihar_lim)-veldxdz_v(i,j+1,k,id_bihar_lim)) & !01-01-23
!.........................................................................................
       *( & !ooooo>
         (anyv3d(i,j+1,k,id_vbef_ )-anyv3d(i,j  ,k,id_vbef_ ))*3.      &
        -(anyv3d(i,j+2,k,id_vbef_ )-anyv3d(i,j-1,k,id_vbef_ ))         &
         *mask_v(i,j+2,kmax)       *mask_v(i,j-1,kmax)                 &
        ) & !ooooo>

                                          ) &     !cccccc>

                                   +anyv3d(i,j,k,id_cnv_)*( & !upupup>

       -0.5*( & !yyyy>
                 anyv3d(i,j,k,id_veldxdz_v_)                      &
               *(anyv3d(i,j,k,id_vbef_)+anyv3d(i,j+1,k,id_vbef_)) &
            +abs(anyv3d(i,j,k,id_veldxdz_v_))                     &
               *(anyv3d(i,j,k,id_vbef_)-anyv3d(i,j+1,k,id_vbef_)) &
                         ) & !yyyy>
                                                            )   !upupup>

      enddo ; enddo
! note4: dz_t plutot que demi somme dz_v car dz_v peut etre n'importe quoi dans le masque

! Advection partielle direction Oj:
      do j=2,jmax ; do i=2,imax-1

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_vdiv)=     &
       anyv3d(i,j,k,id_vdiv)      &  
      +anyv3d(i,j,k,id_vbef_)*dti_lpsub_*( anyv3d(i,j  ,k,id_veldxdz_v_) &
                                          -anyv3d(i,j-1,k,id_veldxdz_v_))/dxdy_v(i,j)

! Couches fusionnees: memoriser champs avant update pour plus tard ne moyenner que la contrib advective
      !anyv3d(i,j,k,id_vbef2_)=anyv3d(i,j,k,id_vbef_) !04-12-19

      anyv3d(i,j,k,id_vbef_)=mask_v(i,j,k)*( & !ooo>
      anyv3d(i,j,k,id_vbef_)*dz_v(i,j,k,dbefore)      & !01-05-19

       +dti_lpsub_*wetmask_v(i,j)*(                   & !--->

             (                                           & !pmx>
                        yflux_t(i,j  )                   &
                       -yflux_t(i,j-1)                   &

      +anyv3d(i,j,k,id_vbef_)*(anyv3d(i,j  ,k,id_veldxdz_v_)   & !16-01-17
                              -anyv3d(i,j-1,k,id_veldxdz_v_) ) &

                                          )/dxdy_v(i,j)  & !pmx>

                                 ) & !--->
                                           ) & !ooo>
                                              /dz_v(i,j,k,dbefore) & !04-11-21
       +(1-mask_v(i,j,k))*vel_v(i,j,k,before)!*dz_v(i,j,k,dbefore)

      enddo ; enddo

      enddo ! k loop

! Sous le fond
     !call vertmix_merged_levels_v(dbefore,id_vbef_,id_vbef2_) !04-12-19

! ICI CONDITION MPI U1 U2 anyv3d(i,j,k,id_vbef_)
      call obc_mpi_anyv3d(1,id_vbef_,'v1')
      call obc_mpi_anyv3d(1,id_vbef_,'v2')


      enddo               !*********> ! loopmax_


      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1  
        vel_v(i,j,k,2)=anyv3d(i,j,k,id_vbef_) 
      enddo ; enddo ; enddo

! Partie implicite de l'advection horizontale:
      if(loopmax_hor_v_>looplimit_hor)call momentum_adv_v_impli(id_veldydz_u_,id_veldxdz_v_)

! ---- ADVECTION VERTICALE ----

! Advection verticale explicite:

!     x1=0. ! Pas obligatoire seulement pour diag du maximum
      xy_v(:,:,1)=small2 ! reset small2 evite que le nombre d'iteration soit nulle
      do k=2,kmax ; do j=2,jmax ; do i=2,imax-1  
! x0=Dz/Dt
       x0=min(dz_v(i,j,k,dbefore),dz_v(i,j,k-1,dbefore))/dti_lp
! Omega bornee par + ou - N fois dz/dt ici 2 fois
       anyv3d(i,j,k,id_we_)=                                    &
       max(min(0.5*(omega_w(i,j,k,1)+omega_w(i,j-1,k,1))        &
              ,looplimit_*x0),-looplimit_*x0)*wetmask_v(i,j)    &
                                          *(1.-sponge_v(i,j,0)) !05-04-22
! max over z of the decimal current number:
       xy_v(i,j,1)=max(xy_v(i,j,1),abs(anyv3d(i,j,k,id_we_))/x0)
!      x1=max(x1,xy_v(i,j,1)) ! Pas obligatoire seulement pour diag du maximum

      enddo       ; enddo       ; enddo

! Fond et couche merged (dans laquelle w est 100% implicite, donc w explitite = 0)
      do j=2,jmax ; do i=2,imax-1  
       do k=1,kmerged_v(i,j)
        anyv3d(i,j,k,id_we_)=0.
       enddo
      enddo       ; enddo

! Surface:
      do j=2,jmax ; do i=2,imax-1  
!      anyv3d(i,j,kmax+1,id_we_)=0.5*(omega_w(i,j,kmax+1,1)+omega_w(i,j-1,kmax+1,1)) ! A la surface w explicite pour flux de surface
! Depuis !07-03-22 vitesse de surface 100% implicite combinee A condition vel(kmax+1)=0 dans triadia_in(:,:,kmax,3)
       anyv3d(i,j,kmax+1,id_we_)=0. !07-03-22
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
      do k=1,kmax+1 ; do j=2,jmax ; do i=2,imax-1  
       anyv3d(i,j,k,id_wi_)=0.5*(omega_w(i,j,k,1)+omega_w(i,j-1,k,1))-anyv3d(i,j,k,id_we_)
      enddo       ; enddo       ; enddo

! Dans la couche merged, profil de vitesse verticale implicite de type "gradient de flux constant"
      do j=2,jmax ; do i=2,imax-1

       do k=2,kmerged_v(i,j)

! Noter qu'on estime depth_w_u (depth_w au point u) comme etant depth_u(k)-0.5*dz_u(k)
         anyv3d(i,j,k,id_wi_)=       anyv3d(i,j,kmerged_v(i,j)+1,id_wi_)    &
      *(depth_v(i,j,k            )-0.5*dz_v(i,j,k,dnow)+h_v(i,j)) &
      /(depth_v(i,j,kmerged_v(i,j)+1)-0.5*dz_v(i,j,kmerged_v(i,j),dnow)+h_v(i,j))

       enddo

      enddo       ; enddo
      check4=anyv3d(imax/2,jmax/2,kmax/2,id_wi_)

      if(flag_timesplitting==1) then !AVEC time-splitting> !28-04-20
       i0=0 ; x1_r4=0.5
      else                           !SANS time-splitting>
       i0=1 ; x1_r4=0.25
      endif                          !SANS time-splitting>

      do j=2,jmax ; do i=2,imax-1  ! Boucle i,j No1

! ceiling de xy_v(i,j,1) nombre de sous-iterations advectives
       loopmax_=min( ceiling(xy_v(i,j,1)) , looplimit_ )
       dti_lpsub=dti_lp/loopmax_

!      do k=kmerged_v(i,j)+1,kmax                         !22-12-14
       do k=2,kmax                         !22-12-14

! Current number with additional limititations
        anyv3d(i,j,k,id_cnw_)=1.-( &!ooooo>
! Robustness regarding rivers and dry zones:
                       wetmask_v(i,j)       & ! wet/dry zone
!      *0.5*(upwindriver_t(i,j)+upwindriver_t(i,j-1)) & ! upwind  zone
! Robustness regarding vertical current number:
        *max(1.-abs(anyv3d(i,j,k,id_we_)*dti_lpsub               &
                 /min(dz_v(i,j,k,dbefore),dz_v(i,j,k-1,dbefore))),0.)        &
                                ) !ooooo>

!        anyv3d(i,j,k,id_cnw_)=0. ! 100%C2
!        anyv3d(i,j,k,id_cnw_)=1. ! 100%UP2

       enddo

! Cas particuliers k=1:kmerged_v(i,j) et k=kmax:
!      anyv3d(i,j,1:kmerged_v(i,j),id_cnw_)=0.
       anyv3d(i,j,1            ,id_cnw_)=0.
       anyv3d(i,j,kmax+1       ,id_cnw_)=0. ! A la surface schema 100% explicite (voir ce qui est
                                            ! fait pour le sel

! Flux centres frozen:
       do k=2,kmax
        anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cnw_))*anyv3d(i,j,k,id_we_)  &
                  *x1_r4*(     vel_v(i,j,k  ,1)+vel_v(i,j,k-1,1)     &  !28-04-20
                          +i0*(vel_v(i,j,k  ,0)+vel_v(i,j,k-1,0)) )
       enddo

       anyv1d(1     ,1)=0.  
       anyv1d(1     ,3)=0.  
       anyv1d(kmax+1,3)=0. 
       k=kmax+1
       anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cnw_))*anyv3d(i,j,k,id_we_)*vel_v(i,j,k-1,1)

! Partie Iterative:
       do loop_=1,loopmax_ ! boucle iterative >>>

       vel_v(i,j,kmax+1,2)=vel_v(i,j,kmax,2)
       vel_v(i,j,0     ,2)=vel_v(i,j,1   ,2)

       do k=2,kmax !nexon>

        x1=0.5*(anyv3d(i,j,k,id_we_)+abs(anyv3d(i,j,k,id_we_)))
        x2=0.5*(anyv3d(i,j,k,id_we_)-abs(anyv3d(i,j,k,id_we_)))

! Flux upwind:
        anyv1d(k,3)=  &
          -anyv3d(i,j,k,id_cnw_) *(x1*vel_v(i,j,k-1,2)+x2*vel_v(i,j,k  ,2))&
       +(1-anyv3d(i,j,k,id_cnw_))*( &
               x1*(vel_v(i,j,k  ,2)-2.*vel_v(i,j,k-1,2)+vel_v(i,j,k-2,2))&!21-11-16
              +x2*(vel_v(i,j,k+1,2)-2.*vel_v(i,j,k  ,2)+vel_v(i,j,k-1,2))&
                                  )*(1+anyv3d(i,j,k,id_cnw_))/6. ! QUICKEST
!                                                           )/6. ! UBS-UPW


       enddo       !nexon>

       do k=1,kmax


! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_vdiv)= & !26-01-17
       anyv3d(i,j,k,id_vdiv)+dti_lpsub*(anyv3d(i,j,k+1,id_we_)    &
                                        -anyv3d(i,j,k  ,id_we_))*vel_v(i,j,k,2) 

         vel_v(i,j,k,2)=vel_v(i,j,k,2)+dti_lpsub*wetmask_v(i,j)*( &

                              anyv1d(k+1,3)+anyv1d(k+1,1) &
                             -anyv1d(k  ,3)-anyv1d(k  ,1) &

                 +vel_v(i,j,k,2)*( anyv3d(i,j,k+1,id_we_)    & !26-01-17
                                  -anyv3d(i,j,k  ,id_we_))   &

                                                 )/dz_v(i,j,k,dbefore)

       enddo

       enddo               ! boucle iterative >>>

      enddo       ; enddo          ! Boucle i,j No1

#ifdef bidon
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1  

           vel_v(i,j,k,2)                                &
          =vel_v(i,j,k,2)                                &

      -wetmask_v(i,j)*anyv3d(i,j,k,id_vdiv)/dz_v(i,j,k,dbefore) !26-01-17

      enddo         ; enddo      ; enddo  ! Boucle i,j No2
#endif
      check8=anyv3d(imax/2,jmax/2,kmax-1,id_vdiv)

! Partie implicite de l'advection horizontale:
!     if(loopmax_hor_v_>looplimit_hor)call momentum_adv_v_impli(id_veldydz_u_,id_veldxdz_v_)

      end subroutine momentum_exp_adv_v_prime

!...................................................................

      subroutine momentum_exp_adv_u_modul !20-10-21
      use module_principal ; use module_parallele
      implicit none
      real*4 :: diflxfactor_
      double precision dti_lpsub_
! Note: anyv3d(:,:,:,0) pris pour distribution verticale du stress de surface
!       anyv3d(:,:,:,3) pris pour nombre de courant vertical
      integer :: id_cnu_=1       &   ! current number relative to u velocity identifier
                ,id_cnv_=2       &   ! current number relative to v velocity identifier
                ,id_veldydz_u_=3 &   ! veldydz_u identifier
                ,id_veldxdz_v_=4 &   ! veldxdz_v identifier
                ,id_ubef_=0      &   ! u velocity identifier
                ,id_ubef2_=9     &   ! u velocity identifier before update
!               ,id_udiv=5       &   ! A definir
                ,looplimit_=1    &   ! nombre d'iteration maximum pour l'advection verticale
                ,id_we_=1        &   ! id_cnu_ n'est plus occupe quand id_we_  entre en service
                ,id_wi_=6        &   
                ,id_cnw_=2       &   ! id_cnv_ n'est plus occupe quand id_cnw_ entre en service
                ,loopmax_hor_u_  &   
                ,loop_,loopmax_      ! u est entirerement terminEe quand commence l'etape v et donc
                                     ! l'espace anyv3d(id_ubef_) est libre quand commence l'etape v

      diflxfactor_=cst_adv_vel/16. 
! note: dans le mode externe on divise par 32 car l'advection y est traitee dans les 2 directions simultanement

! Cas particulier du modele 1DV
      if(flag_1dv==1) then !>>>>> !14-02-16
       vel_u(:,:,:,2)=vel_u(:,:,:,0)
       return
      endif                !>>>>>

!..............................................
! vel_u:
      do k=1,kmax ; do j=0,jmax+1 ; do i=0,imax+2
       anyv3d(i,j,k,id_ubef_)=vel_u(i,j,k,before)
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax  
       anyv3d(i,j,k,id_udiv)=0.
      enddo ; enddo ; enddo

!.............................................
! flux facettes
      do k=1,kmax ; do j=2,jmax-1 ; do i=1,imax
       anyv3d(i,j,k,id_veldydz_u_)=0.5*(veldydz_u(i,j,k,1)+veldydz_u(i+1,j,k,1)) 
      enddo       ; enddo         ; enddo

      do k=1,kmax ; do j=2,jmax ; do i=2,imax
       anyv3d(i,j,k,id_veldxdz_v_)=0.5*(veldxdz_v(i-1,j,k,1)+veldxdz_v(i,j,k,1))
      enddo       ; enddo         ; enddo

!..............................................
! Determiner le nombre de sous-iterations d'advection
! Cet algo suppose que la couche merged est 100% melangee avec les couches sous-jacente
! et donc que l'on a affaire A une couche dont l'epaisseur peut etre consideree comme
! z_w(kmerged+1)-z_w(1) et non pas dz(kmerged) (methode 2). Si on garde dz(kmerged) (methode 1) l'algo considere
! un volume plus petit pour plus de securite (mais evntuellement un dt inutilement trop petit)
! D'autre part les flux qui entrent dans ce volume sont ceux du niveau k=kmerged mais egalement
! ceux sous kmerged, kundermin_u etant le niveau du premier flux lateral non nul
      x1=0.
! Niveaux standard:
      do k=1,kmax !08-11-21
      do j=2,jmax-1 ; do i=2,imax
!       do k=kmerged_u(i,j)+1,kmax
         x1=max(x1,max(max(max(abs(anyv3d(i  ,j  ,k,id_veldydz_u_))  &
                              ,abs(anyv3d(i-1,j  ,k,id_veldydz_u_))) &
                              ,abs(anyv3d(i  ,j  ,k,id_veldxdz_v_))) &
                              ,abs(anyv3d(i  ,j+1,k,id_veldxdz_v_))) & 
                    *dti_lp                                          &
                    *wetmask_u(i,j)                                  &
                      /(dxdy_u(i,j)*min(dz_u(i,j,k,dbefore),dz_u(i,j,k,dafter))))
      enddo       ; enddo
      enddo

      call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_max,par%comm2d,ierr)
!     loopmax_=min(int(x2)+1,100) ! grande valeur expres pour ne pas borner....
      loopmax_=    int(x2)+1      
      dti_lpsub_=dti_lp/loopmax_
      loopmax_hor_u_=loopmax_ ! faire une copie de loopmax_ pour plus tard 

      if(par%rank==0)then
       open(unit=3,file=trim(tmpdirname)//'dti_u',position='append') !17-07-20
        write(3,*)real(elapsedtime_now/86400.),loopmax_,real(x2),real(dti_lpsub_)
       close(3)
      endif

      if(loopmax_>50.or.loopmax_<0) then !>>>
        call momentum_replay_ncumax(id_veldydz_u_,id_veldxdz_v_) !18-01-19
        call graph_out
        stop 'loopmax_ u advection out of range' !15-05-17
      endif                              !>>>

!.............................................
! flux explictes limitEs:
      if(loopmax_>looplimit_hor) then !-limiter-les-flux-explicites->

      do k=1,kmax ; do j=2,jmax-1 ; do i=1,imax
       anyv3d(i,j,k,id_veldydz_u_)=max(min( &

       0.5*(veldydz_u(i,j,k,1)+veldydz_u(i+1,j,k,1))  &

        ,+looplimit_hor*min(dxdy_u(i  ,j)*min(dz_u(i  ,j,k,dbefore),dz_u(i  ,j,k,dafter))/max(wetmask_u(i  ,j),1.e-10) &
                           ,dxdy_u(i+1,j)*min(dz_u(i+1,j,k,dbefore),dz_u(i+1,j,k,dafter))/max(wetmask_u(i+1,j),1.e-10) &
                           )/dti_lp)                                   &

        ,-looplimit_hor*min(dxdy_u(i  ,j)*min(dz_u(i  ,j,k,dbefore),dz_u(i  ,j,k,dafter))/max(wetmask_u(i  ,j),1.e-10) &
                           ,dxdy_u(i+1,j)*min(dz_u(i+1,j,k,dbefore),dz_u(i+1,j,k,dafter))/max(wetmask_u(i+1,j),1.e-10) &
                           )/dti_lp)                                    
      enddo       ; enddo         ; enddo

      do k=1,kmax ; do j=2,jmax ; do i=2,imax

       anyv3d(i,j,k,id_veldxdz_v_)=max(min( &

       0.5*(veldxdz_v(i-1,j,k,1)+veldxdz_v(i,j,k,1)) &

        ,+looplimit_hor*min(dxdy_u(i,j  )*min(dz_u(i,j  ,k,dbefore),dz_u(i,j  ,k,dafter))/max(wetmask_u(i,j  ),1.e-10) &
                           ,dxdy_u(i,j-1)*min(dz_u(i,j-1,k,dbefore),dz_u(i,j-1,k,dafter))/max(wetmask_u(i,j-1),1.e-10) &
                           )/dti_lp)                                   &

        ,+looplimit_hor*min(dxdy_u(i,j  )*min(dz_u(i,j  ,k,dbefore),dz_u(i,j  ,k,dafter))/max(wetmask_u(i,j  ),1.e-10) &
                           ,dxdy_u(i,j-1)*min(dz_u(i,j-1,k,dbefore),dz_u(i,j-1,k,dafter))/max(wetmask_u(i,j-1),1.e-10) &
                           )/dti_lp)                                    

      enddo       ; enddo         ; enddo

      endif                           !-limiter-les-flux-explicites->

!..............................................
! Calculer les nombres de courant correspondant
! NOTE: je considere dti_lp parce que la partie centree de l'advection
! est "figee" sur vel_u et non par sur anyv3d(id_ubef)
      do k=1,kmax
      do j=2,jmax-1 ; do i=1,imax

! Methode 1 lineaire:
!                anyv3d(i,j,k,id_cnu_)=                 &
!            abs(anyv3d(i,j,k,id_veldydz_u_))*dti_lp    &
!               /min(dxdy_u(i+1,j  )*dz_u(i+1,j,k,dbefore)    &
!                   ,dxdy_u(i  ,j  )*dz_u(i  ,j,k,dbefore))

! Methode 2 maintenue centree jusqu'a 0.9 (ou autre):
!                anyv3d(i,j,k,id_cnu_)=min(1.,max(0.,   &
!           (abs(anyv3d(i,j,k,id_veldydz_u_))*dti_lp    &
!               /min(dxdy_u(i+1,j  )*dz_u(i+1,j,k,dbefore)    & !22-08-16
!                   ,dxdy_u(i  ,j  )*dz_u(i  ,j,k,dbefore))   &
!               -0.9)*10.))

! Methode 3 : (j'ai ajoute abs pour pouvoir enlever **2)
!                anyv3d(i,j,k,id_cnu_)=min(1.,                 & !25-11-16
!           (abs(anyv3d(i,j,k,id_veldydz_u_))*dti_lp           &
!               /min(dxdy_u(i+1,j  )*dz_u(i+1,j,k,dbefore)     &
!                   ,dxdy_u(i  ,j  )*dz_u(i  ,j,k,dbefore))) ) & 
!                *wetmask_u(i+1,j  )                           & !05-03-22     
!                *wetmask_u(i  ,j  )

! Methode 4: (nu-1)*upwzone + 1 , =nu si upwzone=1 et =1 si upwzone=0 !05-04-22
      anyv3d(i,j,k,id_cnu_)=min(1., & !min>

      (      & !(nu-1)
       (abs(anyv3d(i,j,k,id_veldydz_u_))*dti_lp          &
           /min(dxdy_u(i+1,j  )*dz_u(i+1,j,k,dbefore)    &
               ,dxdy_u(i  ,j  )*dz_u(i  ,j,k,dbefore)))  &
        -1.) & !(nu-1)
            *wetmask_u(i+1,j)                            &
            *wetmask_u(i  ,j)                            &
         *(1.-sponge_u(i+1,j,0))                         &
         *(1.-sponge_u(i  ,j,0))                         &

         +1.                                               &

                               )      !min>   

      enddo ; enddo
      enddo

      do k=1,kmax
      do j=2,jmax ; do i=2,imax


! Methode 1 lineaire:
!                anyv3d(i,j,k,id_cnv_)=                 &
!            abs(anyv3d(i,j,k,id_veldxdz_v_))*dti_lp    &
!               /min(dxdy_u(i,j-1)*dz_u(i,j-1,k,dbefore)      &
!                   ,dxdy_u(i,j  )*dz_u(i,j  ,k,dbefore))

! Methode 2 maintenue centree jusqu'a 0.9 (ou autre):
!                anyv3d(i,j,k,id_cnv_)=min(1.,max(0.,   &
!           (abs(anyv3d(i,j,k,id_veldxdz_v_))*dti_lp    &
!               /min(dxdy_u(i,j-1)*dz_u(i,j-1,k,dbefore)      &
!                   ,dxdy_u(i,j  )*dz_u(i,j  ,k,dbefore))     &
!               -0.9)*10.))

! Methode 3 : (j'ai ajoute abs pour pouvoir enlever **2)
!                anyv3d(i,j,k,id_cnv_)=min(1.,               &
!           (abs(anyv3d(i,j,k,id_veldxdz_v_))*dti_lp         &
!               /min(dxdy_u(i,j-1)*dz_u(i,j-1,k,dbefore)     &
!                   ,dxdy_u(i,j  )*dz_u(i,j  ,k,dbefore))) ) &   
!                *wetmask_u(i,j-1)                           & !05-03-22     
!                *wetmask_u(i,j  )

! Methode 4 : (nv-1)*upwzone + 1 , =nv si upwzone=1 et =1 si upwzone=0 !05-04-22
      anyv3d(i,j,k,id_cnv_)=min(1., & !min>

      (      & !(nu-1)
       (abs(anyv3d(i,j,k,id_veldxdz_v_))*dti_lp         &
           /min(dxdy_u(i,j-1)*dz_u(i,j-1,k,dbefore)     &
               ,dxdy_u(i,j  )*dz_u(i,j  ,k,dbefore)))   &    
        -1.) & !(nu-1)
            *wetmask_u(i,j-1)                           & 
            *wetmask_u(i,j  )                           &
         *(1.-sponge_u(i,j-1,0))                        & 
         *(1.-sponge_u(i,j  ,0))                        &

         +1.                                            &

                               )      !min>   


      enddo ; enddo
      enddo

!     if(mod(iteration3d,100)==0)write(6,*)'bidouille advection u'
!                anyv3d(:,:,:,id_cnu_)=1.
!                anyv3d(:,:,:,id_cnv_)=1.


!     i=imax/2 ; j=jmax/2 ; k=kmax/2
!     write(6,*)'--------------------------------------------'
!     write(6,*)'anyv3d(i,j,k,id_cnu_)',anyv3d(i,j,k,id_cnu_)
!     write(6,*)'anyv3d(i,j,k,id_cnv_)',anyv3d(i,j,k,id_cnv_)
!     write(6,*)'--------------------------------------------'
!     write(6,*)'loopmax_',loopmax_
!     write(6,*)'--------------------------------------------'
!     write(6,*)'anyv3d(i,j,k,id_veldydz_u_)',anyv3d(i,j,k,id_veldydz_u_)
!     write(6,*)'anyv3d(i,j,k,id_veldxdz_v_)',anyv3d(i,j,k,id_veldxdz_v_)
!     write(6,*)'--------------------------------------------'

      do loop_=1,loopmax_ !*********>

      do k=1,kmax   

! Flux direction Oi
      do j=2,jmax-1 ; do i=1,imax

      xflux_t(i,j)=                                                    &
      
              (1.- anyv3d(i,j,k,id_cnu_))*( &     !cccccc> ! adv scheme selector

! Ox Advective flux at t point:
                     -anyv3d(i,j,k,id_veldydz_u_)                      &
                   *( (vel_u(i  ,j,k,now)+  vel_u(i+1,j,k,now))*0.5625 &
                     -(vel_u(i-1,j,k,now)+  vel_u(i+2,j,k,now))*0.0625)&
! Ox Diffusive flux:
! Details sur le codage du bi-harmonique dans S:
! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit
!.......................................................................................
! Coefficient constant (dans l'espace des indices):
!     +cst_adv_vel*(dxdy_u(i,j)*dz_u(i,j,k,dnow)+dxdy_u(i+1,j)*dz_u(i+1,j,k,dnow))/(32*dti_lp) &
! Coefficient "smago like": !20-11-16
!     +(dy_u(i,j)*dz_u(i,j,k,dnow)+dy_u(i+1,j)*dz_u(i+1,j,k,dnow))    &
!     *abs(anyv3d(i+1,j,k,id_ubef_)-anyv3d(i,j,k,id_ubef_))*cst_adv_vel/64. &
! shema 3 hybride :
!     +(dy_u(i,j)*dz_u(i,j,k,dnow)+dy_u(i+1,j)*dz_u(i+1,j,k,dnow))*0.5*diflxfactor_*( &
!!    +(dy_t(i,j)*dz_t(i,j,k,dnow)                               )    *diflxfactor_*( & !voir note1 !19-07-19
!!        biharm_c2*0.5*abs(anyv3d(i+1,j,k,id_ubef_)-anyv3d(i,j,k,id_ubef_))  & 
!!       +biharm_c1*dx_t(i,j)*inv_dti_lp                                  ) & !09-06-18
!     +diflxfactor_*0.5*(abs(veldydz_u(i,j,k,1))+abs(veldydz_u(i+1,j,k,1))) & !17-01-22 
      +diflxfactor_*0.5*(abs(veldydz_u(i,j,k,id_bihar_lim))+abs(veldydz_u(i+1,j,k,id_bihar_lim))) & !01-01-23
!.......................................................................................
       *( & !ooooo>
         (anyv3d(i+1,j,k,id_ubef_)-anyv3d(i  ,j,k,id_ubef_))*3.        &
        -(anyv3d(i+2,j,k,id_ubef_)-anyv3d(i-1,j,k,id_ubef_))           &
         *mask_u(i+2,j,kmax)      *mask_u(i-1,j,kmax)                  &
        ) & !ooooo>
                                          ) &     !cccccc>

                                   +anyv3d(i,j,k,id_cnu_)*( & !upupup>

        -0.5*( & !xxxx>
                    anyv3d(i,j,k,id_veldydz_u_)                      &
              *(    anyv3d(i,j,k,id_ubef_)+anyv3d(i+1,j,k,id_ubef_)) &
               +abs(anyv3d(i,j,k,id_veldydz_u_))                     & 
                  *(anyv3d(i,j,k,id_ubef_)-anyv3d(i+1,j,k,id_ubef_)) &
             ) & !xxxx>
                                                          )   !upupup>

      enddo ; enddo
! note1: dz_t plutot que demi_somme dz_u car dans masque dz_u peut etre n'importe quoi

! Advection partielle direction Oi:
      do j=2,jmax-1 ; do i=2,imax  


! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_udiv)=     &
       anyv3d(i,j,k,id_udiv)      &  
      +anyv3d(i,j,k,id_ubef_)*dti_lpsub_*( anyv3d(i  ,j,k,id_veldydz_u_) &
                                          -anyv3d(i-1,j,k,id_veldydz_u_))/dxdy_u(i,j)

! Couches fusionnees: memoriser champs avant update pour plus tard ne moyenner que la contrib advective
      !anyv3d(i,j,k,id_ubef2_)=anyv3d(i,j,k,id_ubef_) !04-12-19

      anyv3d(i,j,k,id_ubef_)=mask_u(i,j,k)*( &  !ooo> !25-04-17
      anyv3d(i,j,k,id_ubef_)*dz_u(i,j,k,dbefore)          & !01-05-19

       +dti_lpsub_*wetmask_u(i,j)*(                       & !--->

         (                                               & !pmx>
                         xflux_t(i  ,j)                  &
                        -xflux_t(i-1,j)                  &

      +anyv3d(i,j,k,id_ubef_)*( & 
                          anyv3d(i  ,j,k,id_veldydz_u_)  &
                         -anyv3d(i-1,j,k,id_veldydz_u_)  &
                        ) &
                                           )/dxdy_u(i,j) & !pmx>

                                 )                        & !--->

                                            ) & !ooo>      
                                             /dz_u(i,j,k,dbefore) & !04-11-21
      +(1-mask_u(i,j,k))*vel_u(i,j,k,before)!*dz_u(i,j,k,dbefore)

      enddo ; enddo

      enddo ! boucle k

! Sous le fond
     !call vertmix_merged_levels_u(dbefore,id_ubef_,id_ubef2_) !04-12-19

! ICI CONDITION MPI U1 U2 anyv3d(i,j,k,id_ubef_)
      call obc_mpi_anyv3d(1,id_ubef_,'u1')
      call obc_mpi_anyv3d(1,id_ubef_,'u2')

      do k=1,kmax   

! Flux direction Oj
      do j=2,jmax ; do i=2,imax

      yflux_f(i,j)=                                                    &

              (1.- anyv3d(i,j,k,id_cnv_))*( &     !cccccc> ! adv scheme selector

! Oy Advective flux at p point:
                                      -anyv3d(i,j,k,id_veldxdz_v_)     &
                *(  (vel_u(i  ,j-1,k,now)+  vel_u(i,j  ,k,now))*0.5625 &
                   -(vel_u(i  ,j-2,k,now)+  vel_u(i,j+1,k,now))*0.0625)&
! Oy Diffusive flux:
! Details sur le codage du bi-harmonique dans S:
! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit
!........................................................................................
! Coefficient constant (dans l'espace des indices):
!     +cst_adv_vel*(dxdy_u(i,j)*dz_u(i,j,k,dnow)+dxdy_u(i,j-1)*dz_u(i,j-1,k,dnow))/(32*dti_lp) &
! Coefficient "smago like"
!     +(dy_u(i,j)*dz_u(i,j,k,dnow)+dy_u(i,j-1)*dz_u(i,j-1,k,dnow))    &
!     *abs(anyv3d(i,j,k,id_ubef_)-anyv3d(i,j-1,k,id_ubef_))*cst_adv_vel/64. &
! Schema 3 hybride:
!     +(dy_u(i,j)*dz_u(i,j,k,dnow)+dy_u(i,j-1)*dz_u(i,j-1,k,dnow))*0.5*diflxfactor_*( &
!!     +mask_f(i,j,k)*(dy_u(i,j)*dz_u(i,j,k,dnow)+dy_u(i,j-1)*dz_u(i,j-1,k,dnow))*0.5*diflxfactor_*( & !note2 !19-07-19
!!         biharm_c2*0.5*abs(anyv3d(i,j,k,id_ubef_)-anyv3d(i,j-1,k,id_ubef_))   &
!!        +biharm_c1*dy_f(i,j)*inv_dti_lp                                 ) &
!      +mask_f(i,j,k)*diflxfactor_*0.5*(abs(veldxdz_v(i-1,j,k,1))+abs(veldxdz_v(i,j,k,1))) & !17-01-22
! note: mask_f ne sert A rien si veldxdz_v(i-1:i,....
!     +diflxfactor_*0.5*(abs(veldxdz_v(i-1,j,k,1))+abs(veldxdz_v(i,j,k,1))) & !17-01-22
      +diflxfactor_*0.5*(abs(veldxdz_v(i-1,j,k,id_bihar_lim))+abs(veldxdz_v(i,j,k,id_bihar_lim))) & !01-01-23
!........................................................................................
       *( & !ooooo> 
         (anyv3d(i  ,j  ,k,id_ubef_)-anyv3d(i  ,j-1,k,id_ubef_))*3. &
        -(anyv3d(i  ,j+1,k,id_ubef_)-anyv3d(i  ,j-2,k,id_ubef_))    &
         *mask_u(i  ,j+1,kmax)      *mask_u(i  ,j-2,kmax)           & !25-11-16
        ) & !ooooo>
                                          ) &     !cccccc>

                                   +anyv3d(i,j,k,id_cnv_)*( & !upupup>

        -0.5*( & !yyyy>
                anyv3d(i,j,k,id_veldxdz_v_)                      &
!!!!!!!        0.5*(veldxdz_v(i-1,j,k,1)+veldxdz_v(i,j,k,1))     &
              *(anyv3d(i,j-1,k,id_ubef_)+anyv3d(i,j,k,id_ubef_)) &
           +abs(anyv3d(i,j,k,id_veldxdz_v_))                     &
              *(anyv3d(i,j-1,k,id_ubef_)-anyv3d(i,j,k,id_ubef_)) &
             ) & !yyyy>
                                                          )   !upupup>

      enddo ; enddo

! note2: le flux advectif yflux_f(i,j) est nul si mask_t(i-1,j)=mask_t(i,j)=0 (car alors veldxdz_v(i-1:i,j)=0)
! flux diffusif nul (si *mask_f(i,j,k)) car dz_u peut etre n'importe quoi dans le masque et compromettre CFL
! Par consequent cela correspond A une condition "free-slip" A la cOte.


! Advection partielle direction Oj:
      do j=2,jmax-1 ; do i=2,imax  

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_udiv)=     &
       anyv3d(i,j,k,id_udiv)      &  
      +anyv3d(i,j,k,id_ubef_)*dti_lpsub_*( anyv3d(i,j+1,k,id_veldxdz_v_) &
                                          -anyv3d(i,j  ,k,id_veldxdz_v_))/dxdy_u(i,j)

! Couches fusionnees: memoriser champs avant update pour plus tard ne moyenner que la contrib advective
      !anyv3d(i,j,k,id_ubef2_)=anyv3d(i,j,k,id_ubef_) !04-12-19

      anyv3d(i,j,k,id_ubef_)=mask_u(i,j,k)*( & !ooo>
      anyv3d(i,j,k,id_ubef_)*dz_u(i,j,k,dbefore) & !01-05-19

       +dti_lpsub_*wetmask_u(i,j)*(           & !--->

          (                                  & !pmx>
           +yflux_f(i,j+1)-yflux_f(i  ,j)    &

      +anyv3d(i,j,k,id_ubef_)*( & 
          anyv3d(i,j+1,k,id_veldxdz_v_)      &
         -anyv3d(i,j  ,k,id_veldxdz_v_)      &
                        ) &
           )/dxdy_u(i,j)                     & !pmx>

                   )                          & !--->
                                          )  & !ooo>
                                             /dz_u(i,j,k,dbefore) & !04-11-21
      +(1-mask_u(i,j,k))*vel_u(i,j,k,before)!*dz_u(i,j,k,dbefore)

      enddo ; enddo

      enddo ! k loop

! Sous le fond
     !call vertmix_merged_levels_u(dbefore,id_ubef_,id_ubef2_) !04-12-19

! ICI CONDITION MPI U1 U2 anyv3d(i,j,k,id_ubef_)
      call obc_mpi_anyv3d(1,id_ubef_,'u1')
      call obc_mpi_anyv3d(1,id_ubef_,'u2')

      enddo               !*********> ! loopmax_


      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax  
        vel_u(i,j,k,2)=anyv3d(i,j,k,id_ubef_) 
      enddo ; enddo ; enddo

! Partie implicite de l'advection horizontale:
      if(loopmax_hor_u_>looplimit_hor)call momentum_adv_u_impli(id_veldydz_u_,id_veldxdz_v_)

! ---- ADVECTION VERTICALE ----

! Advection verticale explicite:

!     x1=0. ! Pas obligatoire seulement pour diag du maximum
      xy_u(:,:,1)=small2 ! reset small2 evite que le nombre d'iteration soit nulle
      do k=2,kmax ; do j=2,jmax-1 ; do i=2,imax  
! x0=Dz/Dt
       x0=min(dz_u(i,j,k,dbefore),dz_u(i,j,k-1,dbefore))/dti_lp
! Omega bornee par + ou - N fois dz/dt ici 2 fois
       anyv3d(i,j,k,id_we_)=                                         &
            max(min(0.5*(omega_w(i,j,k,1)+omega_w(i-1,j,k,1))        &
               ,looplimit_*x0),-looplimit_*x0)*wetmask_u(i,j)        &
                                           *(1.-sponge_u(i,j,0)) !05-04-22
! max over z of the decimal current number:
       xy_u(i,j,1)=max(xy_u(i,j,1),abs(anyv3d(i,j,k,id_we_))/x0)
!      x1=max(x1,xy_u(i,j,1)) ! Pas obligatoire seulement pour diag du maximum

      enddo       ; enddo       ; enddo

! Fond et couche merged (dans laquelle w est 100% implicite, donc w explitite = 0)
      do j=2,jmax-1 ; do i=2,imax  
       do k=1,kmerged_u(i,j)
        anyv3d(i,j,k,id_we_)=0.
       enddo
      enddo       ; enddo

! Surface:
      do j=2,jmax-1 ; do i=2,imax  
!      anyv3d(i,j,kmax+1,id_we_)=0.5*(omega_w(i,j,kmax+1,1)+omega_w(i-1,j,kmax+1,1)) ! A la surface w explicite pour flux de surface
! Depuis !07-03-22 vitesse de surface 100% implicite combinee A condition vel(kmax+1)=0 dans triadia_in(:,:,kmax,3)
       anyv3d(i,j,kmax+1,id_we_)=0. !07-03-22
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
      do k=1,kmax+1 ; do j=2,jmax-1 ; do i=2,imax  
       anyv3d(i,j,k,id_wi_)=0.5*(omega_w(i,j,k,1)+omega_w(i-1,j,k,1))-anyv3d(i,j,k,id_we_)
      enddo       ; enddo       ; enddo


! Dans la couche merged, profil de vitesse verticale implicite de type "gradient de flux constant"
      do j=2,jmax-1 ; do i=2,imax  

       do k=2,kmerged_u(i,j)

! Noter qu'on estime depth_w_u (depth_w au point u) comme etant depth_u(k)-0.5*dz_u(k)
         anyv3d(i,j,k,id_wi_)=anyv3d(i,j,kmerged_u(i,j)+1,id_wi_)         &
      *(depth_u(i,j,k            )-0.5*dz_u(i,j,k,dnow)+h_u(i,j)) &
      /(depth_u(i,j,kmerged_u(i,j)+1)-0.5*dz_u(i,j,kmerged_u(i,j),dnow)+h_u(i,j))

       enddo

      enddo       ; enddo
      check3=anyv3d(imax/2,jmax/2,kmax/2,id_wi_)

      if(flag_timesplitting==1) then !AVEC time-splitting> !28-04-20
       i0=0 ; x1_r4=0.5
      else                           !SANS time-splitting>
       i0=1 ; x1_r4=0.25 
      endif                          !SANS time-splitting>

      do j=2,jmax-1 ; do i=2,imax  ! Boucle i,j No1

! ceiling de xy_u(i,j,1) nombre de sous-iterations advectives
       loopmax_=min( ceiling(xy_u(i,j,1)) , looplimit_ )
       dti_lpsub=dti_lp/loopmax_

!      do k=kmerged_u(i,j)+1,kmax                         !22-12-14
       do k=2,kmax                         !22-12-14

!      if(i+par%timax(1)==39.and.j+par%tjmax(1)==10) &
!      write(6,*)'w',k,anyv3d(i,j,k,id_we_),anyv3d(i,j,k,id_wi_)

! Current number with additional limititations
        anyv3d(i,j,k,id_cnw_)=1.-( &!ooooo>
! Robustness regarding rivers and dry zones:
                       wetmask_u(i,j)       & ! wet/dry zone
!      *0.5*(upwindriver_t(i,j)+upwindriver_t(i-1,j)) & ! upwind  zone !17-10-21
! Robustness regarding vertical current number:
        *max(1.-abs(anyv3d(i,j,k,id_we_)*dti_lpsub               &
                 /min(dz_u(i,j,k,dbefore),dz_u(i,j,k-1,dbefore))),0.)        &
                                ) !ooooo>

!        anyv3d(i,j,k,id_cnw_)=0. ! 100%C2
!        anyv3d(i,j,k,id_cnw_)=1. ! 100%UP2

       enddo

! Cas particuliers k=1:kmerged_u(i,j) et k=kmax:
!      anyv3d(i,j,1:kmerged_u(i,j),id_cnw_)=0.
       anyv3d(i,j,1            ,id_cnw_)=0.
       anyv3d(i,j,kmax+1       ,id_cnw_)=0. ! A la surface schema 100% explicite (voir ce qui est
                                            ! fait pour le sel

! Flux centres frozen:
       do k=2,kmax
        anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cnw_))*anyv3d(i,j,k,id_we_)  &
                  *x1_r4*(     vel_u(i,j,k  ,1)+vel_u(i,j,k-1,1)     &
                          +i0*(vel_u(i,j,k  ,0)+vel_u(i,j,k-1,0)) ) !28-04-20
       enddo
       anyv1d(1     ,1)=0.  
       anyv1d(1     ,3)=0.  
       anyv1d(kmax+1,3)=0. 
       k=kmax+1
       anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cnw_))*anyv3d(i,j,k,id_we_)*vel_u(i,j,k-1,1)

! Partie Iterative:
       do loop_=1,loopmax_ ! boucle iterative >>>

       vel_u(i,j,kmax+1,2)=vel_u(i,j,kmax,2)
       vel_u(i,j,0     ,2)=vel_u(i,j,1   ,2)

       do k=2,kmax !nexon>

        x1=0.5*(anyv3d(i,j,k,id_we_)+abs(anyv3d(i,j,k,id_we_)))
        x2=0.5*(anyv3d(i,j,k,id_we_)-abs(anyv3d(i,j,k,id_we_)))

! Flux upwind:
        anyv1d(k,3)=  &
          -anyv3d(i,j,k,id_cnw_) *(x1*vel_u(i,j,k-1,2)+x2*vel_u(i,j,k  ,2))&
       +(1-anyv3d(i,j,k,id_cnw_))*( &
               x1*(vel_u(i,j,k  ,2)-2.*vel_u(i,j,k-1,2)+vel_u(i,j,k-2,2))&!21-11-16
              +x2*(vel_u(i,j,k+1,2)-2.*vel_u(i,j,k  ,2)+vel_u(i,j,k-1,2))&
                                  )*(1+anyv3d(i,j,k,id_cnw_))/6. ! QUICKEST
!                                                           )/6. ! UBS-UPW


       enddo       !nexon>

       do k=1,kmax

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_udiv)= & !26-01-17
       anyv3d(i,j,k,id_udiv)+dti_lpsub*(anyv3d(i,j,k+1,id_we_)    &
                                       -anyv3d(i,j,k  ,id_we_))*vel_u(i,j,k,2) 

         vel_u(i,j,k,2)=vel_u(i,j,k,2)+dti_lpsub*wetmask_u(i,j)*( &

                              anyv1d(k+1,3)+anyv1d(k+1,1) &
                             -anyv1d(k  ,3)-anyv1d(k  ,1) &

                 +vel_u(i,j,k,2)*( anyv3d(i,j,k+1,id_we_)    & !26-01-17
                                  -anyv3d(i,j,k  ,id_we_))   &

                                                 )/dz_u(i,j,k,dbefore)


       enddo

       enddo               ! boucle iterative >>>

      enddo       ; enddo          ! Boucle i,j No1

#ifdef bidon
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax  

           vel_u(i,j,k,2)                                &
          =vel_u(i,j,k,2)                                &

      -wetmask_u(i,j)*anyv3d(i,j,k,id_udiv)/dz_u(i,j,k,dbefore) !26-01-17

      enddo         ; enddo       ; enddo ! Boucle i,j No2
#endif
      check7=anyv3d(imax/2,jmax/2,kmax-1,id_udiv)

! Partie implicite de l'advection horizontale:
!     if(loopmax_hor_u_>looplimit_hor)call momentum_adv_u_impli(id_veldydz_u_,id_veldxdz_v_)

      end subroutine momentum_exp_adv_u_modul

!...................................................................

      subroutine momentum_exp_adv_v_modul !20-10-21
      use module_principal ; use module_parallele
      implicit none
      real*4 :: diflxfactor_
      double precision dti_lpsub_
! Note: anyv3d(:,:,:,0) pris pour distribution verticale du stress de surface
!       anyv3d(:,:,:,3) pris pour nombre de courant vertical
      integer :: id_cnu_=1       &   ! current number relative to u velocity identifier
                ,id_cnv_=2       &   ! current number relative to v velocity identifier
                ,id_veldydz_u_=3 &   ! veldydz_u identifier
                ,id_veldxdz_v_=4 &   ! veldxdz_v identifier
                ,id_vbef_=0      &   ! u velocity identifier
                ,id_vbef2_=9     &   ! u velocity identifier before update
!               ,id_vdiv=5       &   ! A definir
                ,looplimit_=1    &   ! nombre d'iteration maximum pour l'advection verticale
                ,id_we_=1        &   ! id_cnu_ n'est plus occupe quand id_we_  entre en service
                ,id_wi_=7        &   
                ,id_cnw_=2       &   ! id_cnv_ n'est plus occupe quand id_cnw_ entre en service
                ,loopmax_hor_v_  &   
                ,loop_,loopmax_      ! u est entirerement terminEe quand commence l'etape v et donc
                                     ! l'espace anyv3d(id_vbef_) est libre quand commence l'etape v

      diflxfactor_=cst_adv_vel/16.
! note: dans le mode externe on divise par 32 car l'advection y est traitee dans les 2 directions simultanement

! Cas particulier du modele 1DV
      if(flag_1dv==1) then !>>>>> !14-02-16
       vel_v(:,:,:,2)=vel_v(:,:,:,0)
       return
      endif                !>>>>>

!..............................................
! vel_v:
      do k=1,kmax ; do j=0,jmax+2 ; do i=0,imax+1
       anyv3d(i,j,k,id_vbef_)=vel_v(i,j,k,before)
      enddo ; enddo ; enddo

      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1  
       anyv3d(i,j,k,id_vdiv)=0.
      enddo ; enddo ; enddo
!..............................................
! flux facettes:
      do k=1,kmax ; do j=2,jmax ; do i=2,imax
       anyv3d(i,j,k,id_veldydz_u_)=0.5*(veldydz_u(i,j-1,k,1)+veldydz_u(i,j,k,1))      
      enddo       ; enddo       ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=2,imax-1
       anyv3d(i,j,k,id_veldxdz_v_)=0.5*(veldxdz_v(i,j,k,1)+veldxdz_v(i,j+1,k,1))        
      enddo       ; enddo       ; enddo

!..............................................
! Determiner le nombre de sous-iterations d'advection
! Cet algo suppose que la couche merged est 100% melangee avec les couches sous-jacente
! et donc que l'on a affaire A une couche dont l'epaisseur peut etre consideree comme
! z_w(kmerged+1)-z_w(1) et non pas dz(kmerged) (methode 2). Si on garde dz(kmerged) (methode 1) l'algo considere
! un volume plus petit pour plus de securite (mais evntuellement un dt inutilement trop petit)
! D'autre part les flux qui entrent dans ce volume sont ceux du niveau k=kmerged mais egalement
! ceux sous kmerged, kundermin_v etant le niveau du premier flux lateral non nul
      x1=0.
! Niveaux standard:
      do k=1,kmax !08-11-21
      do j=2,jmax ; do i=2,imax-1
!       do k=kmerged_v(i,j)+1,kmax
         x1=max(x1,max(max(max(abs(anyv3d(i  ,j  ,k,id_veldydz_u_))  &
                              ,abs(anyv3d(i+1,j  ,k,id_veldydz_u_))) &
                              ,abs(anyv3d(i  ,j  ,k,id_veldxdz_v_))) &
                              ,abs(anyv3d(i  ,j-1,k,id_veldxdz_v_))) & 
                    *dti_lp                                          &
                    *wetmask_v(i,j)                                  &
                      /(dxdy_v(i,j)*min(dz_v(i,j,k,dbefore),dz_v(i,j,k,dafter))))
      enddo       ; enddo
      enddo

      call mpi_allreduce(x1,x2,1,mpi_double_precision,mpi_max,par%comm2d,ierr)
!     loopmax_=min(int(x2)+1,100) ! grande valeur expres pour ne pas borner....
      loopmax_=    int(x2)+1      
      dti_lpsub_=dti_lp/loopmax_
      loopmax_hor_v_=loopmax_ ! faire une copie de loopmax_ pour plus tard 

      if(par%rank==0)then
       open(unit=3,file=trim(tmpdirname)//'dti_v',position='append') !17-07-20
        write(3,*)real(elapsedtime_now/86400.),loopmax_,real(x2),real(dti_lpsub_)
       close(3)
      endif

      if(loopmax_>50.or.loopmax_<0) then !>>>
        call momentum_replay_ncvmax(id_veldydz_u_,id_veldxdz_v_) !18-01-19
        call graph_out
        stop 'loopmax_ v advection out of range' !15-05-17
      endif                              !>>>

!.............................................
! flux explictes limitEs:
      if(loopmax_>looplimit_hor) then !-limiter-les-flux-explicites->

      do k=1,kmax ; do j=2,jmax ; do i=2,imax

       anyv3d(i,j,k,id_veldydz_u_)=max(min( &

       0.5*(veldydz_u(i,j,k,1)+veldydz_u(i,j-1,k,1))  &

        ,+looplimit_hor*min(dxdy_v(i  ,j)*min(dz_v(i  ,j,k,dbefore),dz_v(i  ,j,k,dafter))/max(wetmask_v(i  ,j),1.e-10) &
                           ,dxdy_v(i-1,j)*min(dz_v(i-1,j,k,dbefore),dz_v(i-1,j,k,dafter))/max(wetmask_v(i-1,j),1.e-10) &
                           )/dti_lp)                                   &

        ,-looplimit_hor*min(dxdy_v(i  ,j)*min(dz_v(i  ,j,k,dbefore),dz_v(i  ,j,k,dafter))/max(wetmask_v(i  ,j),1.e-10) &
                           ,dxdy_v(i-1,j)*min(dz_v(i-1,j,k,dbefore),dz_v(i-1,j,k,dafter))/max(wetmask_v(i-1,j),1.e-10) &
                           )/dti_lp)                                    

      enddo       ; enddo         ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=2,imax-1

       anyv3d(i,j,k,id_veldxdz_v_)=max(min( &

       0.5*(veldxdz_v(i,j+1,k,1)+veldxdz_v(i,j,k,1)) &

        ,+looplimit_hor*min(dxdy_v(i,j  )*min(dz_v(i,j  ,k,dbefore),dz_v(i,j  ,k,dafter))/max(wetmask_v(i,j  ),1.e-10) &
                           ,dxdy_v(i,j+1)*min(dz_v(i,j+1,k,dbefore),dz_v(i,j+1,k,dafter))/max(wetmask_v(i,j+1),1.e-10) &
                           )/dti_lp)                                   &

        ,-looplimit_hor*min(dxdy_v(i,j  )*min(dz_v(i,j  ,k,dbefore),dz_v(i,j  ,k,dafter))/max(wetmask_v(i,j  ),1.e-10) &
                           ,dxdy_v(i,j+1)*min(dz_v(i,j+1,k,dbefore),dz_v(i,j+1,k,dafter))/max(wetmask_v(i,j+1),1.e-10) &
                           )/dti_lp)                                    
      enddo       ; enddo         ; enddo

      endif                           !-limiter-les-flux-explicites->

!..............................................
! Calculer les nombres de courant correspondant
! NOTE: je considere dti_lp parce que la partie centree de l'advection
! est "figee" sur vel_v et non par sur anyv3d(id_ubef)
      do k=1,kmax
      do j=2,jmax ; do i=2,imax

! Methode 3 
!      anyv3d(i,j,k,id_cnu_)=min(1.,                &
!         (abs(anyv3d(i,j,k,id_veldydz_u_))*dti_lp  &
!      /min(dxdy_v(i-1,j)*dz_v(i-1,j,k,dbefore)     &
!          ,dxdy_v(  i,j)*dz_v(i  ,j,k,dbefore))) ) &   
!       *wetmask_v(i-1,j)                           & !05-03-22
!       *wetmask_v(i  ,j)

! Methode 4: (nu-1)*upwzone + 1 , =nu si upwzone=1 et =1 si upwzone=0 !05-04-22
       anyv3d(i,j,k,id_cnu_)=min(1., & !min>

      (      & !(nu-1)
          (abs(anyv3d(i,j,k,id_veldydz_u_))*dti_lp  &
       /min(dxdy_v(i-1,j)*dz_v(i-1,j,k,dbefore)     &
           ,dxdy_v(  i,j)*dz_v(i  ,j,k,dbefore)))   &   
        -1.) & !(nu-1)
            *wetmask_v(i-1,j)                       & 
            *wetmask_v(i  ,j)                       &
         *(1.-sponge_v(i-1,j,0))                    & 
         *(1.-sponge_v(i  ,j,0))                    &

         +1.                                        &

                                )      !min>   

      enddo ; enddo
      enddo

      do k=1,kmax
      do j=1,jmax ; do i=2,imax-1

! Methode 3 
!      anyv3d(i,j,k,id_cnv_)=min(1.,                &
!         (abs(anyv3d(i,j,k,id_veldxdz_v_))*dti_lp  &
!      /min(dxdy_v(i,j+1)*dz_v(i,j+1,k,dbefore)     &
!          ,dxdy_v(i,j  )*dz_v(i,j  ,k,dbefore))) ) &   
!       *wetmask_v(i,j+1)                           & !05-03-22
!       *wetmask_v(i,j  )

! Methode 4 : (nv-1)*upwzone + 1 , =nv si upwzone=1 et =1 si upwzone=0 !05-04-22
       anyv3d(i,j,k,id_cnv_)=min(1., & !min>

      (      & !(nv-1)
          (abs(anyv3d(i,j,k,id_veldxdz_v_))*dti_lp  &
       /min(dxdy_v(i,j+1)*dz_v(i,j+1,k,dbefore)     &
           ,dxdy_v(i,j  )*dz_v(i,j  ,k,dbefore)))   &    
        -1.) & !(nv-1)
            *wetmask_v(i,j+1)                       &
            *wetmask_v(i,j  )                       &
         *(1.-sponge_v(i,j+1,0))                    &
         *(1.-sponge_v(i,j  ,0))                    &

         +1.                                        &

                                )      !min>   

      enddo       ; enddo
      enddo


      do loop_=1,loopmax_ !*********>

      do k=1,kmax   

! Flux direction Oi
      do j=2,jmax ; do i=2,imax

      xflux_f(i,j)=                                                   &

               (1.-anyv3d(i,j,k,id_cnu_))*( &     !cccccc> ! adv scheme selector

! Ox Advective flux at p point:
                     -anyv3d(i,j,k,id_veldydz_u_)                      &
                 *(  (vel_v(i-1,j  ,k,now)+vel_v(i  ,j,k,now))*0.5625  &
                    -(vel_v(i-2,j  ,k,now)+vel_v(i+1,j,k,now))*0.0625) &
! Ox Diffusive flux at p point:
! Details sur le codage du bi-harmonique dans S:
! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit
!.........................................................................................
! Coef constant:
!     +cst_adv_vel*(dxdy_v(i,j)*dz_v(i,j,k,dnow)+dxdy_v(i-1,j)*dz_v(i-1,j,k,dnow))/(32*dti_lp) &
! Coef smago like:
!     +(dx_v(i,j)*dz_v(i,j,k,dnow)+dx_v(i-1,j)*dz_v(i-1,j,k,dnow))       &
!    *abs(anyv3d(i  ,j,k,id_vbef_ )-anyv3d(i-1,j,k,id_vbef_ ))*cst_adv_vel/64. &
! schema 3 hybride
!     +(dx_v(i,j)*dz_v(i,j,k,dnow)+dx_v(i-1,j)*dz_v(i-1,j,k,dnow))*0.5*diflxfactor_*( &
!!    +mask_f(i,j,k)*(dx_v(i,j)*dz_v(i,j,k,dnow)+dx_v(i-1,j)*dz_v(i-1,j,k,dnow))*0.5*diflxfactor_*( & !note3 !19-07-19
!!      biharm_c2*0.5*abs(anyv3d(i  ,j,k,id_vbef_ )-anyv3d(i-1,j,k,id_vbef_ ))  &
!!     +biharm_c1*dx_f(i,j)*inv_dti_lp                                    ) &
!     +mask_f(i,j,k)*diflxfactor_*0.5*(abs(veldydz_u(i,j-1,k,1))+abs(veldydz_u(i,j,k,1))) & !17-01-22
! note: mask_f ne sert A rien si veldydz_u(j-1:j....
!     +diflxfactor_*0.5*(abs(veldydz_u(i,j-1,k,1))+abs(veldydz_u(i,j,k,1))) & !17-01-22
      +diflxfactor_*0.5*(abs(veldydz_u(i,j-1,k,id_bihar_lim))+abs(veldydz_u(i,j,k,id_bihar_lim))) & !01-01-23
!.........................................................................................
       *( & !ooooo>
         (anyv3d(i  ,j,k,id_vbef_ )-anyv3d(i-1,j,k,id_vbef_ ))*3.      &
        -(anyv3d(i+1,j,k,id_vbef_ )-anyv3d(i-2,j,k,id_vbef_ ))         &
         *mask_v(i+1,j,kmax)       *mask_v(i-2,j,kmax)                 &
        ) & !ooooo>

                                          ) &     !cccccc>

                                   +anyv3d(i,j,k,id_cnu_)*( & !upupup>

         -0.5*( & !xxxx>
                 anyv3d(i,j,k,id_veldydz_u_)                       &
               *(anyv3d(i-1,j,k,id_vbef_)+anyv3d(i,j,k,id_vbef_))  &
            +abs(anyv3d(i,j,k,id_veldydz_u_))                      &
               *(anyv3d(i-1,j,k,id_vbef_)-anyv3d(i,j,k,id_vbef_))  &
                         ) & !xxxx>
                                                          )   !upupup>


      enddo ; enddo
! note3: flux diffusif nul (si *mask_f) car dz_v peut etre n'importe quoi dans le masque

! Advection partielle direction Oi:
      do j=2,jmax ; do i=2,imax-1

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_vdiv)=     &
       anyv3d(i,j,k,id_vdiv)      &  
      +anyv3d(i,j,k,id_vbef_)*dti_lpsub_*( anyv3d(i+1,j,k,id_veldydz_u_) &
                                          -anyv3d(i  ,j,k,id_veldydz_u_))/dxdy_v(i,j)

! Couches fusionnees: memoriser champs avant update pour plus tard ne moyenner que la contrib advective
      !anyv3d(i,j,k,id_vbef2_)=anyv3d(i,j,k,id_vbef_) !04-12-19

      anyv3d(i,j,k,id_vbef_)=mask_v(i,j,k)*( & !ooo>
      anyv3d(i,j,k,id_vbef_)*dz_v(i,j,k,dbefore)          & !01-05-19

       +dti_lpsub_*wetmask_v(i,j)*( & !--->

             (                                            & !pmx>
                         xflux_f(i+1,j)                   &
                        -xflux_f(i  ,j)                   &

      +anyv3d(i,j,k,id_vbef_)*( anyv3d(i+1,j,k,id_veldydz_u_)   & !16-01-17
                               -anyv3d(i  ,j,k,id_veldydz_u_) ) &

                                           )/dxdy_v(i,j)  & !pmx>

               ) & !--->
                                           ) & !ooo>
                                              /dz_v(i,j,k,dbefore) & !04-11-21
       +(1-mask_v(i,j,k))*vel_v(i,j,k,before)!*dz_v(i,j,k,dbefore)

      enddo ; enddo

      enddo ! boucle k

! Sous le fond
     !call vertmix_merged_levels_v(dbefore,id_vbef_,id_vbef2_) !04-12-19

! ICI CONDITION MPI v1 v2 anyv3d(i,j,k,id_vbef_)
      call obc_mpi_anyv3d(1,id_vbef_,'v1')
      call obc_mpi_anyv3d(1,id_vbef_,'v2')

      do k=1,kmax   

! Flux direction Oj
      do j=1,jmax ; do i=2,imax-1

      yflux_t(i,j)=                                                  &

               (1.-anyv3d(i,j,k,id_cnv_))*( &     !cccccc> ! adv scheme selector

! Oy Advective flux at t point:
                     -anyv3d(i,j,k,id_veldxdz_v_)                      &
                  *( (vel_v(i,j  ,k,now)+vel_v(i,j+1,k,now))*0.5625    &
                    -(vel_v(i,j-1,k,now)+vel_v(i,j+2,k,now))*0.0625)   &
! Oy Diffusive flux at t point:
! Details sur le codage du bi-harmonique dans S:
! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit
!.........................................................................................
! Coef constant:
!     +cst_adv_vel*(dxdy_v(i,j)*dz_v(i,j,k,dnow)+dxdy_v(i,j+1)*dz_v(i,j+1,k,dnow))/(32*dti_lp) &
! Coef smago like:
!     +(dx_v(i,j)*dz_v(i,j,k,dnow)+dx_v(i,j+1)*dz_v(i,j+1,k,dnow))       &
!    *abs(anyv3d(i,j+1,k,id_vbef_ )-anyv3d(i,j  ,k,id_vbef_ ))*cst_adv_vel/64. &
! schema 3 hybride
!     +(dx_v(i,j)*dz_v(i,j,k,dnow)+dx_v(i,j+1)*dz_v(i,j+1,k,dnow))*0.5*diflxfactor_*( &
!!    +(dx_t(i,j)*dz_t(i,j,k,dnow)                               )    *diflxfactor_*( & !note4 !19-07-19
!!      biharm_c2*0.5*abs(anyv3d(i,j+1,k,id_vbef_ )-anyv3d(i,j  ,k,id_vbef_ ))  &
!!     +biharm_c1*dy_t(i,j)*inv_dti_lp                                    ) &
!     +diflxfactor_*0.5*(abs(veldxdz_v(i,j,k,1))+abs(veldxdz_v(i,j+1,k,1))) & !17-01-22
      +diflxfactor_*0.5*(abs(veldxdz_v(i,j,k,id_bihar_lim))+abs(veldxdz_v(i,j+1,k,id_bihar_lim))) & !01-01-23
!.........................................................................................
       *( & !ooooo>
         (anyv3d(i,j+1,k,id_vbef_ )-anyv3d(i,j  ,k,id_vbef_ ))*3.      &
        -(anyv3d(i,j+2,k,id_vbef_ )-anyv3d(i,j-1,k,id_vbef_ ))         &
         *mask_v(i,j+2,kmax)       *mask_v(i,j-1,kmax)                 &
        ) & !ooooo>
                                          ) &     !cccccc>

                                   +anyv3d(i,j,k,id_cnv_)*( & !upupup>

       -0.5*( & !yyyy>
                 anyv3d(i,j,k,id_veldxdz_v_)                      &
               *(anyv3d(i,j,k,id_vbef_)+anyv3d(i,j+1,k,id_vbef_)) &
            +abs(anyv3d(i,j,k,id_veldxdz_v_))                     &
               *(anyv3d(i,j,k,id_vbef_)-anyv3d(i,j+1,k,id_vbef_)) &
                         ) & !yyyy>
                                                            )   !upupup>

      enddo ; enddo
! note4: dz_t plutot que demi somme dz_v car dz_v peut etre n'importe quoi dans le masque

! Advection partielle direction Oj:
      do j=2,jmax ; do i=2,imax-1

! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_vdiv)=     &
       anyv3d(i,j,k,id_vdiv)      &  
      +anyv3d(i,j,k,id_vbef_)*dti_lpsub_*( anyv3d(i,j  ,k,id_veldxdz_v_) &
                                          -anyv3d(i,j-1,k,id_veldxdz_v_))/dxdy_v(i,j)

! Couches fusionnees: memoriser champs avant update pour plus tard ne moyenner que la contrib advective
      !anyv3d(i,j,k,id_vbef2_)=anyv3d(i,j,k,id_vbef_) !04-12-19

      anyv3d(i,j,k,id_vbef_)=mask_v(i,j,k)*( & !ooo>
      anyv3d(i,j,k,id_vbef_)*dz_v(i,j,k,dbefore)      & !01-05-19

       +dti_lpsub_*wetmask_v(i,j)*(                   & !--->

             (                                           & !pmx>
                        yflux_t(i,j  )                   &
                       -yflux_t(i,j-1)                   &

      +anyv3d(i,j,k,id_vbef_)*(anyv3d(i,j  ,k,id_veldxdz_v_)   & !16-01-17
                              -anyv3d(i,j-1,k,id_veldxdz_v_) ) &

                                          )/dxdy_v(i,j)  & !pmx>

                                 ) & !--->
                                           ) & !ooo>
                                              /dz_v(i,j,k,dbefore) & !04-11-21
       +(1-mask_v(i,j,k))*vel_v(i,j,k,before)!*dz_v(i,j,k,dbefore)

      enddo ; enddo

      enddo ! k loop

! Sous le fond
     !call vertmix_merged_levels_v(dbefore,id_vbef_,id_vbef2_) !04-12-19

! ICI CONDITION MPI U1 U2 anyv3d(i,j,k,id_vbef_)
      call obc_mpi_anyv3d(1,id_vbef_,'v1')
      call obc_mpi_anyv3d(1,id_vbef_,'v2')


      enddo               !*********> ! loopmax_


      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1  
        vel_v(i,j,k,2)=anyv3d(i,j,k,id_vbef_) 
      enddo ; enddo ; enddo

! Partie implicite de l'advection horizontale:
      if(loopmax_hor_v_>looplimit_hor)call momentum_adv_v_impli(id_veldydz_u_,id_veldxdz_v_)

! ---- ADVECTION VERTICALE ----

! Advection verticale explicite:

!     x1=0. ! Pas obligatoire seulement pour diag du maximum
      xy_v(:,:,1)=small2 ! reset small2 evite que le nombre d'iteration soit nulle
      do k=2,kmax ; do j=2,jmax ; do i=2,imax-1  
! x0=Dz/Dt
       x0=min(dz_v(i,j,k,dbefore),dz_v(i,j,k-1,dbefore))/dti_lp
! Omega bornee par + ou - N fois dz/dt ici 2 fois
       anyv3d(i,j,k,id_we_)=                                    &
       max(min(0.5*(omega_w(i,j,k,1)+omega_w(i,j-1,k,1))        &
              ,looplimit_*x0),-looplimit_*x0)*wetmask_v(i,j)    &
                                          *(1.-sponge_v(i,j,0)) !05-04-22
! max over z of the decimal current number:
       xy_v(i,j,1)=max(xy_v(i,j,1),abs(anyv3d(i,j,k,id_we_))/x0)
!      x1=max(x1,xy_v(i,j,1)) ! Pas obligatoire seulement pour diag du maximum

      enddo       ; enddo       ; enddo

! Fond et couche merged (dans laquelle w est 100% implicite, donc w explitite = 0)
      do j=2,jmax ; do i=2,imax-1  
       do k=1,kmerged_v(i,j)
        anyv3d(i,j,k,id_we_)=0.
       enddo
      enddo       ; enddo

! Surface:
      do j=2,jmax ; do i=2,imax-1  
!      anyv3d(i,j,kmax+1,id_we_)=0.5*(omega_w(i,j,kmax+1,1)+omega_w(i,j-1,kmax+1,1)) ! A la surface w explicite pour flux de surface
! Depuis !07-03-22 vitesse de surface 100% implicite combinee A condition vel(kmax+1)=0 dans triadia_in(:,:,kmax,3)
       anyv3d(i,j,kmax+1,id_we_)=0. !07-03-22
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
      do k=1,kmax+1 ; do j=2,jmax ; do i=2,imax-1  
       anyv3d(i,j,k,id_wi_)=0.5*(omega_w(i,j,k,1)+omega_w(i,j-1,k,1))-anyv3d(i,j,k,id_we_)
      enddo       ; enddo       ; enddo

! Dans la couche merged, profil de vitesse verticale implicite de type "gradient de flux constant"
      do j=2,jmax ; do i=2,imax-1

       do k=2,kmerged_v(i,j)

! Noter qu'on estime depth_w_u (depth_w au point u) comme etant depth_u(k)-0.5*dz_u(k)
         anyv3d(i,j,k,id_wi_)=       anyv3d(i,j,kmerged_v(i,j)+1,id_wi_)    &
      *(depth_v(i,j,k            )-0.5*dz_v(i,j,k,dnow)+h_v(i,j)) &
      /(depth_v(i,j,kmerged_v(i,j)+1)-0.5*dz_v(i,j,kmerged_v(i,j),dnow)+h_v(i,j))

       enddo

      enddo       ; enddo
      check4=anyv3d(imax/2,jmax/2,kmax/2,id_wi_)

      if(flag_timesplitting==1) then !AVEC time-splitting> !28-04-20
       i0=0 ; x1_r4=0.5
      else                           !SANS time-splitting>
       i0=1 ; x1_r4=0.25
      endif                          !SANS time-splitting>

      do j=2,jmax ; do i=2,imax-1  ! Boucle i,j No1

! ceiling de xy_v(i,j,1) nombre de sous-iterations advectives
       loopmax_=min( ceiling(xy_v(i,j,1)) , looplimit_ )
       dti_lpsub=dti_lp/loopmax_

!      do k=kmerged_v(i,j)+1,kmax                         !22-12-14
       do k=2,kmax                         !22-12-14

! Current number with additional limititations
        anyv3d(i,j,k,id_cnw_)=1.-( &!ooooo>
! Robustness regarding rivers and dry zones:
                       wetmask_v(i,j)       & ! wet/dry zone
!      *0.5*(upwindriver_t(i,j)+upwindriver_t(i,j-1)) & ! upwind  zone
! Robustness regarding vertical current number:
        *max(1.-abs(anyv3d(i,j,k,id_we_)*dti_lpsub               &
                 /min(dz_v(i,j,k,dbefore),dz_v(i,j,k-1,dbefore))),0.)        &
                                ) !ooooo>

!        anyv3d(i,j,k,id_cnw_)=0. ! 100%C2
!        anyv3d(i,j,k,id_cnw_)=1. ! 100%UP2

       enddo

! Cas particuliers k=1:kmerged_v(i,j) et k=kmax:
!      anyv3d(i,j,1:kmerged_v(i,j),id_cnw_)=0.
       anyv3d(i,j,1            ,id_cnw_)=0.
       anyv3d(i,j,kmax+1       ,id_cnw_)=0. ! A la surface schema 100% explicite (voir ce qui est
                                            ! fait pour le sel

! Flux centres frozen:
       do k=2,kmax
        anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cnw_))*anyv3d(i,j,k,id_we_)  &
                  *x1_r4*(     vel_v(i,j,k  ,1)+vel_v(i,j,k-1,1)     &  !28-04-20
                          +i0*(vel_v(i,j,k  ,0)+vel_v(i,j,k-1,0)) )
       enddo

       anyv1d(1     ,1)=0.  
       anyv1d(1     ,3)=0.  
       anyv1d(kmax+1,3)=0. 
       k=kmax+1
       anyv1d(k,1)=-(1-anyv3d(i,j,k,id_cnw_))*anyv3d(i,j,k,id_we_)*vel_v(i,j,k-1,1)

! Partie Iterative:
       do loop_=1,loopmax_ ! boucle iterative >>>

       vel_v(i,j,kmax+1,2)=vel_v(i,j,kmax,2)
       vel_v(i,j,0     ,2)=vel_v(i,j,1   ,2)

       do k=2,kmax !nexon>

        x1=0.5*(anyv3d(i,j,k,id_we_)+abs(anyv3d(i,j,k,id_we_)))
        x2=0.5*(anyv3d(i,j,k,id_we_)-abs(anyv3d(i,j,k,id_we_)))

! Flux upwind:
        anyv1d(k,3)=  &
          -anyv3d(i,j,k,id_cnw_) *(x1*vel_v(i,j,k-1,2)+x2*vel_v(i,j,k  ,2))&
       +(1-anyv3d(i,j,k,id_cnw_))*( &
               x1*(vel_v(i,j,k  ,2)-2.*vel_v(i,j,k-1,2)+vel_v(i,j,k-2,2))&!21-11-16
              +x2*(vel_v(i,j,k+1,2)-2.*vel_v(i,j,k  ,2)+vel_v(i,j,k-1,2))&
                                  )*(1+anyv3d(i,j,k,id_cnw_))/6. ! QUICKEST
!                                                           )/6. ! UBS-UPW


       enddo       !nexon>

       do k=1,kmax


! Cumul du pivot obligatoirement avant advection partielle 06-02-17
       anyv3d(i,j,k,id_vdiv)= & !26-01-17
       anyv3d(i,j,k,id_vdiv)+dti_lpsub*(anyv3d(i,j,k+1,id_we_)    &
                                        -anyv3d(i,j,k  ,id_we_))*vel_v(i,j,k,2) 

         vel_v(i,j,k,2)=vel_v(i,j,k,2)+dti_lpsub*wetmask_v(i,j)*( &

                              anyv1d(k+1,3)+anyv1d(k+1,1) &
                             -anyv1d(k  ,3)-anyv1d(k  ,1) &

                 +vel_v(i,j,k,2)*( anyv3d(i,j,k+1,id_we_)    & !26-01-17
                                  -anyv3d(i,j,k  ,id_we_))   &

                                                 )/dz_v(i,j,k,dbefore)

       enddo

       enddo               ! boucle iterative >>>

      enddo       ; enddo          ! Boucle i,j No1

#ifdef bidon
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1  

           vel_v(i,j,k,2)                                &
          =vel_v(i,j,k,2)                                &

      -wetmask_v(i,j)*anyv3d(i,j,k,id_vdiv)/dz_v(i,j,k,dbefore) !26-01-17

      enddo         ; enddo      ; enddo  ! Boucle i,j No2
#endif
      check8=anyv3d(imax/2,jmax/2,kmax-1,id_vdiv)

! Partie implicite de l'advection horizontale:
!     if(loopmax_hor_v_>looplimit_hor)call momentum_adv_v_impli(id_veldydz_u_,id_veldxdz_v_)

      end subroutine momentum_exp_adv_v_modul

!...................................................................

      subroutine momentum_adv_u_impli(id_veldydz_exp_u_,id_veldxdz_exp_v_) !31-12-22
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      integer :: t_=1        &
                ,loop_       &
                ,flag_mpi_=0 &
                ,work_status_=0     & !0=rien A signaler, 1=travail demandE       , 2=travail fait
                ,send_status_=0     & !0=rien A signaler, 1=envoi mpi demandE     , 2=envoi mpi rEalisE
                ,recv_status_=0     & !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
                ,status_voisin_upst & ! 0= pas d'echange avec voisin amont, 1= echange
                ,status_voisin_down & ! 0= pas d'echange avec voisin aval , 1= echange
                ,work_status_mpi_min_ !04-01-23
      integer id_veldydz_imp_u_  &
             ,id_veldxdz_imp_v_  &
             ,id_veldydz_exp_u_  &
             ,id_veldxdz_exp_v_

!.............................................
! flux facettes residuels implicites pour equation de u

! On reprend le meme emplacement memoire desormais vacant:
      id_veldydz_imp_u_=id_veldydz_exp_u_  
      id_veldxdz_imp_v_=id_veldxdz_exp_v_  

      do k=1,kmax ; do j=2,jmax-1 ; do i=1,imax
       anyv3d(i,j,k,id_veldydz_imp_u_)= &
      -anyv3d(i,j,k,id_veldydz_exp_u_)+0.5*(veldydz_u(i,j,k,1)+veldydz_u(i+1,j,k,1)) 
      enddo       ; enddo         ; enddo

      do k=1,kmax ; do j=2,jmax ; do i=2,imax
       anyv3d(i,j,k,id_veldxdz_imp_v_)= &
      -anyv3d(i,j,k,id_veldxdz_exp_v_)+0.5*(veldxdz_v(i-1,j,k,1)+veldxdz_v(i,j,k,1))
      enddo       ; enddo         ; enddo

      if(iperiodicboundary==.true.)stop 'conditions periodiques A faire'
      if(jperiodicboundary==.true.)stop 'conditions periodiques A faire'

! flag_mpi_=0: je ne fais rien
! flag_mpi_=1: je fais le calcul
! flag_mpi_=2: j'ai fait le calcul, j'envoie la valeur au voisin
! flag_mpi_=3: j'ai tout fait, je m'arrete

! Schema 100% implicite decommenter les lignes suivantes, sinon le cumul de anyv3d(i,j,k,id_udiv) se poursuit
!     do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax 
!      anyv3d(i,j,k,id_udiv)=0. ! reset cumul temperature*divergence(courant) !26-01-17
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
!    noter: debut de boucle en i=2 utilisant U centree sur point T en i-1=1
! Voisin amont:
      if(maxval(anyv3d(1,2:jmax-1,1:kmax,id_veldydz_imp_u_))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: U au point T en i=imax-1 du present rank est U au point T en i=1 du voisin aval
      if(maxval(anyv3d(imax-1,2:jmax-1,1:kmax,id_veldydz_imp_u_))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(anyv3d(1:imax-1,2:jmax-1,1:kmax,id_veldydz_imp_u_))<=0.) then !>>> ! raison 1
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

         call momentum_adv_u_impobc(1,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax

       vel_u(i,j,k,2)=mask_u(i,j,k)*( & !LSM>

         vel_u(i,j,k,2) &

        +( & !haut>

         (anyv3d(i-1,j,k,id_veldydz_imp_u_)+abs(anyv3d(i-1,j,k,id_veldydz_imp_u_)))*dti_lp &
        *(vel_u(i-1,j,k,2)-vel_u(i,j,k,2))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_u(i,j)*dz_u(i,j,k,dbefore)) &

        +(anyv3d(i-1,j,k,id_veldydz_imp_u_)+abs(anyv3d(i-1,j,k,id_veldydz_imp_u_)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_u(i,j,k))*vel_u(i,j,k,2)    


      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(1,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->
   10 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax

       anyv3d(i,j,k,id_udiv)= & 
       anyv3d(i,j,k,id_udiv)  &
      +vel_u(i,j,k,2)*dti_lp*0.5*( &
                                          anyv3d(i  ,j,k,id_veldydz_imp_u_)+abs(anyv3d(i  ,j,k,id_veldydz_imp_u_)) &
                                         -anyv3d(i-1,j,k,id_veldydz_imp_u_)-abs(anyv3d(i-1,j,k,id_veldydz_imp_u_)) &
                                         )/dxdy_u(i,j)

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
!    noter: la boucle finit en i=imax utilisant une vitesse U centree sur point T en i=imax 
! Voisin amont:
      if(minval(anyv3d(imax,2:jmax-1,1:kmax,id_veldydz_imp_u_))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: U centree sur point T en i=2 du present rank est U centree sur point T en i=imax du voisin aval
      if(minval(anyv3d(2,2:jmax-1,1:kmax,id_veldydz_imp_u_))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(anyv3d(2:imax,2:jmax-1,1:kmax,id_veldydz_imp_u_))>=0.) then !>>> ! raison 1
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

      call momentum_adv_u_impobc(2,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=2,jmax-1 ; do i=imax,2,-1

       vel_u(i,j,k,2)=mask_u(i,j,k)*( & !LSM>

         vel_u(i,j,k,2) &

        +( & !haut>

         (-anyv3d(i  ,j,k,id_veldydz_imp_u_)+abs(anyv3d(i  ,j,k,id_veldydz_imp_u_)))*dti_lp &
        *(vel_u(i+1,j,k,2)-vel_u(i,j,k,2))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_u(i,j)*dz_u(i,j,k,dbefore)) &

        +(-anyv3d(i  ,j,k,id_veldydz_imp_u_)+abs(anyv3d(i  ,j,k,id_veldydz_imp_u_)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_u(i,j,k))*vel_u(i,j,k,2)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(2,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->
   11 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax

       anyv3d(i,j,k,id_udiv)= & 
       anyv3d(i,j,k,id_udiv)  &
      +vel_u(i,j,k,2)*dti_lp*0.5*( &
                                          anyv3d(i  ,j,k,id_veldydz_imp_u_)-abs(anyv3d(i  ,j,k,id_veldydz_imp_u_)) &
                                         -anyv3d(i-1,j,k,id_veldydz_imp_u_)+abs(anyv3d(i-1,j,k,id_veldydz_imp_u_)) &
                                         )/dxdy_u(i,j)

      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans momentum_adv_u_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur vel_u(i,j,k,2)
!     call obc_mpi_vel_u(1,2,'zb')

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
!    noter: la boucle commence en j=2 avec vitesse centree "F" en j=2
! Voisin amont:
      if(maxval(anyv3d(2:imax,2,1:kmax,id_veldxdz_imp_v_))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                   status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: V centree "F" en j=jmax du present rank est V centree "F" en j=2 du voisin aval
      if(maxval(anyv3d(2:imax,jmax,1:kmax,id_veldxdz_imp_v_))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(anyv3d(2:imax,2:jmax-1,1:kmax,id_veldxdz_imp_v_))<=0.) then !>>> ! raison 1
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

         call momentum_adv_u_impobc(3,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax

       vel_u(i,j,k,2)=mask_u(i,j,k)*( & !LSM>

         vel_u(i,j,k,2) &

        +( & !haut>

         ( anyv3d(i,j,k,id_veldxdz_imp_v_)+abs(anyv3d(i,j,k,id_veldxdz_imp_v_)))*dti_lp &
        *(vel_u(i,j-1,k,2)-vel_u(i,j,k,2))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_u(i,j)*dz_u(i,j,k,dbefore)) &

        +( anyv3d(i,j,k,id_veldxdz_imp_v_)+abs(anyv3d(i,j,k,id_veldxdz_imp_v_)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_u(i,j,k))*vel_u(i,j,k,2)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(3,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1


      enddo                 !-boucle-mpi->
   20 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax

       anyv3d(i,j,k,id_udiv)= & 
       anyv3d(i,j,k,id_udiv)  &
      +vel_u(i,j,k,2)*dti_lp*0.5*( &
                                          anyv3d(i,j+1,k,id_veldxdz_imp_v_)+abs(anyv3d(i,j+1,k,id_veldxdz_imp_v_)) &
                                         -anyv3d(i,j  ,k,id_veldxdz_imp_v_)-abs(anyv3d(i,j  ,k,id_veldxdz_imp_v_)) &
                                         )/dxdy_u(i,j)

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
!    noter: la boucle finit en j=jmax-1 en utilisant V centree "F" en j=jmax
! Voisin amont:
      if(minval(anyv3d(2:imax,jmax,1:kmax,id_veldxdz_imp_v_))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: V centree "F" en j=2 du present rank est V centree "F" en j=jmax du voisin aval
      if(minval(anyv3d(2:imax,2,1:kmax,id_veldxdz_imp_v_))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                   status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(anyv3d(2:imax,3:jmax,1:kmax,id_veldxdz_imp_v_))>=0.) then !>>> ! raison 1
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

      call momentum_adv_u_impobc(4,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=jmax-1,2,-1 ; do i=2,imax

       vel_u(i,j,k,2)=mask_u(i,j,k)*( & !LSM>

         vel_u(i,j,k,2) &

        +( & !haut>

         (-anyv3d(i,j+1,k,id_veldxdz_imp_v_)+abs(anyv3d(i,j+1,k,id_veldxdz_imp_v_)))*dti_lp &
        *(vel_u(i,j+1,k,2)-vel_u(i,j,k,2))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_u(i,j)*dz_u(i,j,k,dbefore)) &

        +(-anyv3d(i,j+1,k,id_veldxdz_imp_v_)+abs(anyv3d(i,j+1,k,id_veldxdz_imp_v_)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_u(i,j,k))*vel_u(i,j,k,2)    


      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(4,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->
   21 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=2,jmax-1 ; do i=2,imax

       anyv3d(i,j,k,id_udiv)= & 
       anyv3d(i,j,k,id_udiv)  &
      +vel_u(i,j,k,2)*dti_lp*0.5*( &
                                          anyv3d(i,j+1,k,id_veldxdz_imp_v_)-abs(anyv3d(i,j+1,k,id_veldxdz_imp_v_)) &
                                         -anyv3d(i,j  ,k,id_veldxdz_imp_v_)+abs(anyv3d(i,j  ,k,id_veldxdz_imp_v_)) &
                                         )/dxdy_u(i,j)

      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans momentum_adv_u_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur vel_u(i,j,k,2)
!     call obc_mpi_vel_u(1,2,'zb')

      end subroutine momentum_adv_u_impli

!..............................................................................

      subroutine momentum_adv_u_impobc(case_,idvar_,send_status_,recv_status_) !31-12-22
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
        i=imax-1
        do k=1,kmax ; do j=2,jmax-1
         k0=k0+1
         adv_send_est(k0)=vel_u(i,j,k,2)
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
!      k0=2*kmax*jmax
!      k0=  (kmax-1)*(jmax-2) ! puisque do k=2,kmax ; do j=2,jmax-1
       k0=   kmax   *(jmax-2) ! puisque do k=1,kmax ; do j=2,jmax-1
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
       i=1
       k0=0
       do k=1,kmax ; do j=2,jmax-1
        k0=k0+1
        vel_u(i,j,k,2)=adv_recv_ouest(k0)
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
        do k=1,kmax ; do j=2,jmax-1
         k0=k0+1
         adv_send_ouest(k0)=vel_u(i,j,k,2)
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
!      k0=2*kmax*jmax
!      k0=  (kmax-1)*(jmax-2) ! puisque do k=2,kmax ; do j=2,jmax-1
       k0=   kmax   *(jmax-2) ! puisque do k=1,kmax ; do j=2,jmax-1
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
       do k=1,kmax ; do j=2,jmax-1
        k0=k0+1
        vel_u(i,j,k,2)=adv_recv_est(k0)
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
        j=jmax-1
        do k=1,kmax ; do i=2,imax
         k0=k0+1
         adv_send_nord(k0)=vel_u(i,j,k,2)
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
!      k0=2*kmax*imax
       k0=   kmax   *(imax-1) ! puisque do k=1,kmax ; do i=2,imax
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
       j=1
       k0=0
       do k=1,kmax ; do i=2,imax
        k0=k0+1
        vel_u(i,j,k,2)=adv_recv_sud(k0)
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
        j=2
        do k=1,kmax ; do i=2,imax
         k0=k0+1
         adv_send_sud(k0)=vel_u(i,j,k,2)
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
!      k0=2*kmax*imax
       k0=   kmax   *(imax-1) ! puisque do k=1,kmax ; do i=2,imax
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
       j=jmax
       k0=0
       do k=1,kmax ; do i=2,imax
        k0=k0+1
        vel_u(i,j,k,2)=adv_recv_nord(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                       !m°v°m>


      endif             !---echange-Oj-boucle-decroissante--->

      end subroutine momentum_adv_u_impobc

!....................................................

      subroutine momentum_adv_v_impli(id_veldydz_exp_u_,id_veldxdz_exp_v_) !31-12-22
      use module_principal
      use module_parallele
      use module_my_outputs
      implicit none
      integer :: t_=1        &
                ,loop_       &
                ,flag_mpi_=0 &
                ,work_status_=0     & !0=rien A signaler, 1=travail demandE       , 2=travail fait
                ,send_status_=0     & !0=rien A signaler, 1=envoi mpi demandE     , 2=envoi mpi rEalisE
                ,recv_status_=0     & !0=rien A signaler, 1=reception mpi demandEe, 2=reception mpi rEalisEe
                ,status_voisin_upst & ! 0= pas d'echange avec voisin amont, 1= echange
                ,status_voisin_down & ! 0= pas d'echange avec voisin aval , 1= echange
                ,work_status_mpi_min_
      integer id_veldydz_imp_u_  &
             ,id_veldxdz_imp_v_  &
             ,id_veldydz_exp_u_  &
             ,id_veldxdz_exp_v_

      if(iperiodicboundary==.true.)stop 'conditions periodiques A faire'
      if(jperiodicboundary==.true.)stop 'conditions periodiques A faire'

!.............................................
! flux facettes residuels implicites pour equation de v

! On reprend le meme emplacement memoire desormais vacant:
      id_veldydz_imp_u_=id_veldydz_exp_u_  
      id_veldxdz_imp_v_=id_veldxdz_exp_v_  

      do k=1,kmax ; do j=2,jmax ; do i=2,imax
       anyv3d(i,j,k,id_veldydz_imp_u_)= &
      -anyv3d(i,j,k,id_veldydz_exp_u_)+0.5*(veldydz_u(i,j-1,k,1)+veldydz_u(i,j,k,1))      
      enddo       ; enddo       ; enddo

      do k=1,kmax ; do j=1,jmax ; do i=2,imax-1
       anyv3d(i,j,k,id_veldxdz_imp_v_)= &
      -anyv3d(i,j,k,id_veldxdz_exp_v_)+0.5*(veldxdz_v(i,j,k,1)+veldxdz_v(i,j+1,k,1))        
      enddo       ; enddo       ; enddo

! flag_mpi_=0: je ne fais rien
! flag_mpi_=1: je fais le calcul
! flag_mpi_=2: j'ai fait le calcul, j'envoie la valeur au voisin
! flag_mpi_=3: j'ai tout fait, je m'arrete

! Schema 100% implicite decommenter les lignes suivantes, sinon le cumul de anyv3d(i,j,k,id_vdiv) se poursuit
!     do k=1,kmax ; do j=2,jmax ; do i=2,imax-1 
!      anyv3d(i,j,k,id_vdiv)=0. ! reset cumul temperature*divergence(courant) !26-01-17
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
!    noter: debut de boucle en i=2 utilisant U centree sur point F en i=2
! Voisin amont:
      if(maxval(anyv3d(2,2:jmax,1:kmax,id_veldydz_imp_u_))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: U au point F en i=imax du present rank est U au point F en i=2 du voisin aval
      if(maxval(anyv3d(imax,2:jmax,1:kmax,id_veldydz_imp_u_))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(anyv3d(2:imax-1,2:jmax,1:kmax,id_veldydz_imp_u_))<=0.) then !>>> ! raison 1
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

         call momentum_adv_v_impobc(1,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1

       vel_v(i,j,k,2)=mask_v(i,j,k)*( & !LSM>

         vel_v(i,j,k,2) &

        +( & !haut>

         (anyv3d(i  ,j,k,id_veldydz_imp_u_)+abs(anyv3d(i  ,j,k,id_veldydz_imp_u_)))*dti_lp &
        *(vel_v(i-1,j,k,2)-vel_v(i,j,k,2))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_v(i,j)*dz_v(i,j,k,dbefore)) &

        +(anyv3d(i  ,j,k,id_veldydz_imp_u_)+abs(anyv3d(i  ,j,k,id_veldydz_imp_u_)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_v(i,j,k))*vel_v(i,j,k,2)    


      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(1,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->
   10 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1

       anyv3d(i,j,k,id_vdiv)= & 
       anyv3d(i,j,k,id_vdiv)  &
      +vel_v(i,j,k,2)*dti_lp*0.5*( &
                                          anyv3d(i+1,j,k,id_veldydz_imp_u_)+abs(anyv3d(i+1,j,k,id_veldydz_imp_u_)) &
                                         -anyv3d(i  ,j,k,id_veldydz_imp_u_)-abs(anyv3d(i  ,j,k,id_veldydz_imp_u_)) &
                                         )/dxdy_v(i,j)

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
!    noter: la boucle finit en i=imax-1 utilisant une vitesse U centree sur point F en i=imax 
! Voisin amont:
      if(minval(anyv3d(imax,2:jmax,1:kmax,id_veldydz_imp_u_))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(est)==mpi_proc_null)                      status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: U centree sur point F en i=2 du present rank est U centree sur point F en i=imax du voisin aval
      if(minval(anyv3d(2,2:jmax,1:kmax,id_veldydz_imp_u_))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(ouest)==mpi_proc_null)               status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(anyv3d(3:imax,2:jmax,1:kmax,id_veldydz_imp_u_))>=0.) then !>>> ! raison 1
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

      call momentum_adv_v_impobc(2,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=2,jmax ; do i=imax-1,2,-1

       vel_v(i,j,k,2)=mask_v(i,j,k)*( & !LSM>

         vel_v(i,j,k,2) &

        +( & !haut>

         (-anyv3d(i+1,j,k,id_veldydz_imp_u_)+abs(anyv3d(i+1,j,k,id_veldydz_imp_u_)))*dti_lp &
        *(vel_v(i+1,j,k,2)-vel_v(i,j,k,2))                   &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_v(i,j)*dz_v(i,j,k,dbefore)) &

        +(-anyv3d(i+1,j,k,id_veldydz_imp_u_)+abs(anyv3d(i+1,j,k,id_veldydz_imp_u_)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_v(i,j,k))*vel_v(i,j,k,2)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(2,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->
   11 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1

       anyv3d(i,j,k,id_vdiv)= & 
       anyv3d(i,j,k,id_vdiv)  &
      +vel_v(i,j,k,2)*dti_lp*0.5*( &
                                          anyv3d(i+1,j,k,id_veldydz_imp_u_)-abs(anyv3d(i+1,j,k,id_veldydz_imp_u_)) &
                                         -anyv3d(i  ,j,k,id_veldydz_imp_u_)+abs(anyv3d(i  ,j,k,id_veldydz_imp_u_)) &
                                         )/dxdy_v(i,j)

      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans momentum_adv_v_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur vel_v(i,j,k,2)
!     call obc_mpi_vel_v(1,2,'zb')

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
!    noter: la boucle commence en j=2 avec vitesse centree "T" en j=1
! Voisin amont:
      if(maxval(anyv3d(2:imax-1,1,1:kmax,id_veldxdz_imp_v_))<=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                   status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: V centree "T" en j=jmax-1 du present rank est V centree "T" en j=1 du voisin aval
      if(maxval(anyv3d(2:imax-1,jmax-1,1:kmax,id_veldxdz_imp_v_))<=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant mEme d'avoir commencE dans le cas des vitesses positives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(maxval(anyv3d(2:imax-1,1:jmax-1,1:kmax,id_veldxdz_imp_v_))<=0.) then !>>> ! raison 1
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

         call momentum_adv_v_impobc(3,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1

       vel_v(i,j,k,2)=mask_v(i,j,k)*( & !LSM>

         vel_v(i,j,k,2) &

        +( & !haut>

         ( anyv3d(i,j-1,k,id_veldxdz_imp_v_)+abs(anyv3d(i,j-1,k,id_veldxdz_imp_v_)))*dti_lp &
        *(vel_v(i,j-1,k,2)-vel_v(i,j,k,2))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_v(i,j)*dz_v(i,j,k,dbefore)) &

        +( anyv3d(i,j-1,k,id_veldxdz_imp_v_)+abs(anyv3d(i,j-1,k,id_veldxdz_imp_v_)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_v(i,j,k))*vel_v(i,j,k,2)    

      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(3,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1


      enddo                 !-boucle-mpi->
   20 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1

       anyv3d(i,j,k,id_vdiv)= & 
       anyv3d(i,j,k,id_vdiv)  &
      +vel_v(i,j,k,2)*dti_lp*0.5*( &
                                          anyv3d(i,j  ,k,id_veldxdz_imp_v_)+abs(anyv3d(i,j  ,k,id_veldxdz_imp_v_)) &
                                         -anyv3d(i,j-1,k,id_veldxdz_imp_v_)-abs(anyv3d(i,j-1,k,id_veldxdz_imp_v_)) &
                                         )/dxdy_v(i,j)

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
!    noter: la boucle finit en j=jmax en utilisant V centree "T" en j=jmax
! Voisin amont:
      if(minval(anyv3d(2:imax-1,jmax,1:kmax,id_veldxdz_imp_v_))>=0.)status_voisin_upst=0 !raison 2
      if(par%tvoisin(nord)==mpi_proc_null)                     status_voisin_upst=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier

! Voisin aval:
! Note: V centree "T" en j=2 du present rank est V centree "T" en j=jmax du voisin aval
      if(minval(anyv3d(2:imax-1,2,1:kmax,id_veldxdz_imp_v_))>=0.)status_voisin_down=0 !raison 2
      if(par%tvoisin(sud)==mpi_proc_null)                   status_voisin_down=0 !raison 1
! Note: attention A l'ordre des raisons. La raison 1 est plus forte donc elle s'applique en dernier
!...............................

! Si il n'y a pas de voisin amont, le travail commence dEs la premiere iteration:
      if(status_voisin_upst==0)work_status_=1

! Les raisons pour decreter que le travail est fini avant d'avoir commencE dans le cas des vitesses negatives
! 1- Aucune vitesse positive sur toute la boucle 3D
      if(minval(anyv3d(2:imax-1,2:jmax,1:kmax,id_veldxdz_imp_v_))>=0.) then !>>> ! raison 1
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

      call momentum_adv_v_impobc(4,2,send_status_,recv_status_)

      if(work_status_==1) then !-work->
      do k=1,kmax ; do j=jmax,2,-1 ; do i=2,imax-1

       vel_v(i,j,k,2)=mask_v(i,j,k)*( & !LSM>

         vel_v(i,j,k,2) &

        +( & !haut>

         (-anyv3d(i,j  ,k,id_veldxdz_imp_v_)+abs(anyv3d(i,j  ,k,id_veldxdz_imp_v_)))*dti_lp &
        *(vel_v(i,j+1,k,2)-vel_v(i,j,k,2))                    &

         ) & !haut>      
        /( & !bas>

         (2.*dxdy_v(i,j)*dz_v(i,j,k,dbefore)) &

        +(-anyv3d(i,j  ,k,id_veldxdz_imp_v_)+abs(anyv3d(i,j  ,k,id_veldxdz_imp_v_)))*dti_lp &

         ) & !bas>

                                            ) & !LSM>
                                             +(1-mask_v(i,j,k))*vel_v(i,j,k,2)    


      enddo       ; enddo       ; enddo
      flag_mpi_=1    ! pour faire une demande de travail au voisin
      work_status_=2 ! travail fait
      endif                 !-work->

      call advection_scal_flagobc(4,flag_mpi_)
! si advection_scal_flagobc renvoie flag_mpi_=1 alors que work_status_=0 c'est que work_status_ doit passer =1
      if(work_status_==0.and.flag_mpi_==1)work_status_=1

      enddo                 !-boucle-mpi->
   21 continue

! Cumul du pivot obligatoirement apres advection partielle 
      do k=1,kmax ; do j=2,jmax ; do i=2,imax-1

       anyv3d(i,j,k,id_vdiv)= & 
       anyv3d(i,j,k,id_vdiv)  &
      +vel_v(i,j,k,2)*dti_lp*0.5*( &
                                          anyv3d(i,j  ,k,id_veldxdz_imp_v_)-abs(anyv3d(i,j  ,k,id_veldxdz_imp_v_)) &
                                         -anyv3d(i,j-1,k,id_veldxdz_imp_v_)+abs(anyv3d(i,j-1,k,id_veldxdz_imp_v_)) &
                                         )/dxdy_v(i,j)

      enddo ; enddo ; enddo

! Pas d'echange car les echanges ont eu lieu dans momentum_adv_v_impobc
! ICI CONDITION MPI ZB=Z1Z2 sur vel_v(i,j,k,2)
!     call obc_mpi_vel_v(1,2,'zb')

      end subroutine momentum_adv_v_impli

!..............................................................................

      subroutine momentum_adv_v_impobc(case_,idvar_,send_status_,recv_status_) !31-12-22
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
        i=imax-1
        do k=1,kmax ; do j=2,jmax
         k0=k0+1
         adv_send_est(k0)=vel_v(i,j,k,2)
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
!      k0=2*kmax*jmax
       k0=   kmax   *(jmax-1) ! puisque do k=1,kmax ; do j=2,jmax
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
       i=1
       k0=0
       do k=1,kmax ; do j=2,jmax
        k0=k0+1
        vel_v(i,j,k,2)=adv_recv_ouest(k0)
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
        i=2
        do k=1,kmax ; do j=2,jmax
         k0=k0+1
         adv_send_ouest(k0)=vel_v(i,j,k,2)
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
!      k0=2*kmax*jmax
!      k0=  (kmax-1)*(jmax-2) ! puisque do k=2,kmax ; do j=2,jmax-1
       k0=   kmax   *(jmax-1) ! puisque do k=1,kmax ; do j=2,jmax
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
       i=imax
       k0=0
       do k=1,kmax ; do j=2,jmax
        k0=k0+1
        vel_v(i,j,k,2)=adv_recv_est(k0)
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
        j=jmax-1
        do k=1,kmax ; do i=2,imax-1
         k0=k0+1
         adv_send_nord(k0)=vel_v(i,j,k,2)
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
!      k0=2*kmax*imax
       k0=   kmax   *(imax-2) ! puisque do k=1,kmax ; do i=2,imax-1
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
       j=1
       k0=0
       do k=1,kmax ; do i=2,imax-1
        k0=k0+1
        vel_v(i,j,k,2)=adv_recv_sud(k0)
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
        do k=1,kmax ; do i=2,imax-1
         k0=k0+1
         adv_send_sud(k0)=vel_v(i,j,k,2)
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
!      k0=2*kmax*imax
       k0=   kmax   *(imax-2) ! puisque do k=1,kmax ; do i=2,imax-1
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
       do k=1,kmax ; do i=2,imax-1
        k0=k0+1
        vel_v(i,j,k,2)=adv_recv_nord(k0)
       enddo       ; enddo
       recv_status_=2 ! pour memoriser le fait que la reception est faite
      endif                    !>>>          
      endif                                       !m°v°m>


      endif             !---echange-Oj-boucle-decroissante--->

      end subroutine momentum_adv_v_impobc

!....................................................
#ifdef bidon
      subroutine momentum_bihar_lim_u !31-12-22
      use module_principal
      implicit none

! Calculer la valeur abolue (bien noter la valeur absolue) des
! composantes des flux en la limitant par un nombre de courant de 1,
! pour etre plus tard utilisee dans le coef de viscosite biharmonique

      stop 'pour utiliser cette routine il faut...'
! remplacer les lignes abs(vel( blablabla ) par les tabeaux suivant
! Bien noter que la valeur absolue est deja introduite dans le tableau
! REVERIFIER les limites des boucles svp

! Limiter les flux de la diffusion bi-harmonique en imposant la limite  d'un de courant unitaire 
! sous jacente A diflxfactor_=cst_adv_vel/16. Voir egalement:
! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit?usp=sharing
      do k=1,kmax ; do j=2,jmax-1 ; do i=1,imax+1 ! u u
       veldydz_u(i,j,k,id_abs_u_lim)=min(abs(veldydz_u(i,j,k,1)),dy_u(i,j)*dz_u(i,j,k,1)*dx_u(i,j)*inv_dti_fw)
      enddo       ; enddo         ; enddo
      do k=1,kmax ; do j=2,jmax ; do i=1,imax ! u v
       veldxdz_v(i,j,k,id_abs_v_lim)=min(abs(veldxdz_v(i,j,k,1)),dy_v(i,j)*dz_v(i,j,k,1)*dx_v(i,j)*inv_dti_fw)
      enddo       ; enddo         ; enddo

      end subroutine momentum_bihar_lim_u
#endif
!........................................
#ifdef bidon

      subroutine momentum_bihar_lim_v !31-12-22
      use module_principal
      implicit none

! Calculer la valeur abolue (bien noter la valeur absolue) des
! composantes des flux en la limitant par un nombre de courant de 1,
! pour etre plus tard utilisee dans le coef de viscosite biharmonique

      stop 'pour utiliser cette routine il faut...'
! remplacer les lignes abs(vel( blablabla ) par les tabeaux suivant
! Bien noter que la valeur absolue est deja introduite dans le tableau
! REVERIFIER les limites des boucles svp

! Limiter les flux de la diffusion bi-harmonique en imposant la limite  d'un de courant unitaire 
! sous jacente A diflxfactor_=cst_adv_vel/16. Voir egalement:
! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit?usp=sharing
      do k=1,kmax ; do j=1,jmax ; do i=2,imax ! v u
       veldydz_u(i,j,k,id_abs_u_lim)=min(abs(veldydz_u(i,j,k,1)),dy_u(i,j)*dz_u(i,j,k,1)*dx_u(i,j)*inv_dti_fw)
      enddo       ; enddo         ; enddo
      do k=1,kmax ; do j=1,jmax+1 ; do i=2,imax-1 ! v v
       veldxdz_v(i,j,k,id_abs_v_lim)=min(abs(veldxdz_v(i,j,k,1)),dy_v(i,j)*dz_v(i,j,k,1)*dx_v(i,j)*inv_dti_fw)
      enddo       ; enddo         ; enddo

      end subroutine momentum_bihar_lim_v

#endif
!........................................
      subroutine momentum_bihar_limiter !01-01-23
      use module_principal
      use module_parallele
      implicit none

! Limiter les flux de la diffusion bi-harmonique en imposant la limite  d'un de courant unitaire 
! sous jacente A diflxfactor_=cst_adv_vel/16. Voir egalement:
! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit?usp=sharing

! Note: un nombre de courant unitaire (avec dt =  dti_fw comme ici) correspond 
! A un nombre de courant de 2 avec dt=dti_lp. Le schema d'advection utilisant dti_lp, la limite suivante 
! est donc atteinte alors qu'on a deja 2 sous pas de temps advectifs

! Les boucles sont faites pour convenir pour les 2 equations des moments (celle pour u et celle pour v)
      do k=1,kmax ; do j=1,jmax ; do i=1,imax+1 ! u u
               veldydz_u(i,j,k,id_bihar_lim)= &
       max(min(veldydz_u(i,j,k,1),+dy_u(i,j)*dz_u(i,j,k,1)*dx_u(i,j)*inv_dti_fw) &
                                 ,-dy_u(i,j)*dz_u(i,j,k,1)*dx_u(i,j)*inv_dti_fw)
      enddo       ; enddo         ; enddo
      do k=1,kmax ; do j=1,jmax+1 ; do i=1,imax ! u v
               veldxdz_v(i,j,k,id_bihar_lim)= &
       max(min(veldxdz_v(i,j,k,1),+dy_v(i,j)*dz_v(i,j,k,1)*dx_v(i,j)*inv_dti_fw) &
                                 ,-dy_v(i,j)*dz_v(i,j,k,1)*dx_v(i,j)*inv_dti_fw)
      enddo       ; enddo         ; enddo

      end subroutine momentum_bihar_limiter !31-12-22
