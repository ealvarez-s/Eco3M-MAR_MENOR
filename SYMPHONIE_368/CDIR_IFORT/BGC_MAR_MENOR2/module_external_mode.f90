










      module module_external_mode
      use module_principal ; use module_parallele ; use module_mangrove ; use module_s
      use module_q ; use module_webcanals
      implicit none
      real*4 :: fbfc1=0.02, fbfc3=0.003  , fbfc2
!     real*4 :: fbfc1=0.01, fbfc3=0.0015 , fbfc2
!     real*4 :: fbfc1=0.  , fbfc3=0.     , fbfc2
      integer :: iteration2d_prev_=0
!     fbfc2=1.-fbfc1

!     integer noise_mode
!     double precision, dimension(:,:) , allocatable :: velbar0_u_,velbar0_v_
!______________________________________________________________________
! SYMPHONIE ocean model
! release 367 - last update: 21-02-23
!______________________________________________________________________
!...............................................................................
! Version date     Description des modifications
!         26/12/02: mise en service
!         04/03/03: passage � la coordonnee hybride
!         25/03/03: suppression de CALL ANIMATION
!         22/02/04: ajout gradient de pression atmospherique de surface
!         05/07/06: viscosite horizontale calculee sur vitesse prises au temps t-1
!                   filtre d'asselin remplace par diffusion temporelle decentree
!                   en arriere
!         22/07/06: grosses modifs: le mode interne n'utilise plus la moyenne
!                   de certains processus barotrope. Donc disparition de trans_x
!                   et de trans_y et de extavera, ce qui permet de regrouper des
!                   termes figes initialiement sur 2 tableau xy(1) et xy(2) en un
!                   tableau xy(1)
!         07/08/06: Modes externes parallelles
!         10/08/06: La deduction de velbar � partir de fluxbar ayant �t� diff�r�e
!                   pour que la diffusion horizontale soit calcul�e � partir
!                   des courants au temps t-1, il faut deplacer le calcul de
!                   fric_z pour conserver un calcul au temps t du module de la
!                   friction
!         15/09/06: Coef de diffusion horizontal non constant, qui implique
!                   par ailleurs la prise en compte d'un terme de couplage du
!                   mode externe par le mode interne: DIF3D2D
!         20/09/06: debug du point precedent, les constantes devant �tre divisee
!                   par 2
!         03/10/06: 3 options de diffusion horizontale possibles:
!                   op1: Coef Diffusivit� constante
!                   op2: Coef diffusivit� proportionnelle �: abs(u+u)*dx/4
!                        qui combin� � un schema d'advection centr� tend
!                        vers le schema d'advection upwind.
!                   op3: comme op2 mais (abs(u)+abs(u))*dx/2 mieux adapt�
!                        au bruit 2dx
!         09/11/06: Mise en place d'un filtre de m�lange des 2 modes externes
!         17/04/07: Passage � coordonnees curvilignes (voir ajout dx_y et cie...)
!         31/05/07: bienvenue � ZTAREFSAVE_I VBRREFSAVE_I ZTAREFSAVE_J
!                   VBRREFSAVE_J
!         22/02/08: En attendant un posttraitement adapt�, la version
!                   par defaut supprime le barometre inverse. Pour le
!                   reintroduire, agir directement dans le code (dans la
!                   presente routine). Agir egalement dans set_paremeters.F
!         02/04/08: option 4 adopt�e pour melange horizontal
!         09/03/09: regroupement des obc pour clarifier le travail de parallelisation
!         28/04/09: Bancs decouvrants revus
!         09/05/09: suppression de verrou sur dxb different de dyb
!         04-06-09  Conditions aux limites sur i=0,=imax+1 j=0,=jmax+1 regroupee
!                   avec C.L. du mode externe
!         06-06-09  Boucle pour calcul zta_z demarrees � i=0 et j=0 pour
!                   parallelisation
! 2009.2  01-09-09  externe_fast renomme externe
!         16-09-09  Synchro des proc � la fin de chaque iteration
! 2009.3  02-10-09  dte_lp remplace dte
!         03-10-09  externe .F devient external_mode.F
!         05-10-09  ajout d'un "ifdef parallele"
!         14-10-09  re-introduction d'une parenth�se betement supprim�e lors de la
!                   construction de la version parallelisee
!         16-10-09  barometre inverse activ�
! 2010.2  25-12-09  potentide_z renomm� tidepotential_z
!         26-12-09  nettoyage de constantes inutilisees
! 2010.3  12-01-10  - utilisation wetdry_cts2
!                   - termes baroclines masqu�s par wetmask
! 2010.3  14-01-10  cst_adv_vel n'est plus multipli� par 0.25 (attention � reduire
!                   la valeur en entree de notebook_advection
! 2010.7  22-02-10  La variable d'etat devient velbar
!         26-02-10  - Ajout velbarstokes
!                   - Plus de rigueur dans la terminologie des extensions des
!                     tableaux utilis�s pour coriolis
!         03-03-10  Stokes Vortex Force
!         04-03-10  Regroupement des pressions d'equilibre dans presgrad.F
! 2010.8  11-03-10  4 levels time filter
!         12-03-10  mixer le filtre 3L et 4L avec tfc1 tfc2
!         17-03-10  - retour filtre 2L
!                   - etiquette de compilation pour stokes
!         18-03-10  Appliquer le wetmask sur les vitesses de stokes
!         22-03-10  - re-ecriture pour diminuer temps de calculs.
!                   - ajout de hssh_w = ssh_ext_w + h_w
!         24-03-10  Pour garantir le flux de masse aux rivieres la C.L.
!                   des rivieres est directement appliquee � fluxbar
!         02-04-10  parametrisation ardhuin suite: frottement sur le fond
!         03-05-10  - Pour gain cpu, dif3d2d est n�glig� (� l'essai)
!                   - Pour gain cpu, force de vortex fig� sur calcul complet issu
!                     du mode interne (routine stokesvortexforce.F90). (� l'essai)
! 2010.12 18-09-10  Le message d'alerte sur iwve est pass� dans routine set_parameters
! 2010.16 12-01-11  Filtre temporel orde �lev�
! 2010.17 13-01-11  Un protocole de verification de la conservation de mpi
! 2010.18 28-01-11  Possibilit� de tester le vrai filtre d'Assslin
!         15-02-11  Asselin s'appelle avant "Move Forward"
! 2010.19 13-04-11  L'option de complation checkmpi permet de tester la parallelisation
!                   en remplacement de l'ancienne option de compilation bidouille
! 2010.20 19-04-11  - k2dfin renomme iteration2d_max_now
!                   - suppression lignes mode 1 "double"
! 2010.24 01-09-11  routine check mpi
!         15-12-11  Amelioration de la routine de verification de mpi
! 2010.25 23-02-12  wstresb remplace C.L. "rappel vers 1.5*us"
!         02-03-12  reorganisation de boucles pour accelerer le calcul
!         28-03-12  echange compatible avec pgf Stelios
! S26.1   04-09-12  Cas sans mode externe
!         27-01-13  advection qdm o4
!         01-03-13  calcul ssh_avr_w moyenne de la ssh sur cycle barotrope
!         10-04-13  Ajouter l'increment de force de vortex qui decoule de l'increment
!                   sur le courant moyen
!         15-05-13  flag_externalmode_lf_or_fb d�clar� dans module_principal
!         23-06-13  ajout mangrove
!         30-11-13  blindage wetdrying
!         14-03-14  La modif precedente est annulee pour etre reportee dans z_thickness.F90
!         24-04-14  Un meilleur affichage des messages des erreurs checkmpi
!         30-06-14  nouveaux echanges
!         11-11-14  l'introduction de fluxbar(:,:,0) pour completer fluxbar(:,:,1)
!                   renforce les proprietes de conservation de la QDM
!         13-11-14  nouveaux noms de routines suite a re-organisation obc_ext
!         12-01-15  boucles sur wetmask aggrandies pour coherence avec calcul
!                   du flux dans le cas de la prise en compte du forcage par les vagues
!         09-04-15  evaporation/precipition = terme puit/source ssh
!         10-07-15  rhpzavr est une moyenne verticale ponderee dependant de x,y
!         18-08-15  - ajout d'une condition aux limites pour le flux de QDM
!                   - boucles individuelles pour xflux_t et cie
!         27-11-15  check mpi cas 1DV
!         22-12-15  test mpi_proc_null
!         19-01-16  Cas test baroclinic jet:
!                   1- Possibilite d'une friction de fond lineaire
!                   2- zonal restoring
!         20-01-16  Cas test baroclinic jet: moyenne zonale calculEe sur velbarobc-velbar
!                   velbot3d2d calculE avec velbar_u et non plus velavr
!         28-01-16  Nouveau filtre hybride 2dt-4dt remplace suppression du mode
!                   numerique par fourrier
!         30-01-16  flux advectif sur 2 temps
!         28-03-16  correction d'un bug introduit dans la version precedente
!         31-03-16  suppression call mpi_barrier(par%comm2d,k_out) 
!         09-04-16  ajout d'une routine pour forcage academique par une onde de surface
!         13-05-16  - Module friction moyennE sur 2 temps
!                   - wetmask_u et wetmask_v basEs sur min de ssh(0) et ssh(1) 
!         10-08-16  nouvel algo zone intertidale
!         22-09-16  call stokesforces_vortex_external commentE pour cause d'instabilitE
!         31-10-16  nouvel algo zone intertidale: d'autres variantes
!         11-03-17  Zone intertidale: schema 4
!         06-04-17  diflx est divisE par 8 pour etre homogene avec le mode interne
!         08-04-17  s-z coordinate
!         30-04-17  filtrage du mode numerique dans velbot
!         14-05-17  - biharm_c1 biharm_c2 pour un coef viscositE basEe sur une
!                   proportion du gradient de vitesse et de la vitesse
!                   - biharm_2dfactor permet d'amplifier la viscositE du mode externe
!         15-05-17  Stop si ssh_int_w(imax/2,jmax/2,2) is NAN
!         25-05-17  appel A obc_ext_qdm commentE
!         04-06-17  call s_stop_by_interrupteur_file
!         07-06-17  filtre temporel multipliE par wetmask
!         08-06-17  - aiguillage velbot si kmax=1 ou kmax/=1
!                   - filtre temporel integrE dans le calcul global
!                   - bornage H banc decouvrant
!         18-10-17  merge avec version m0v0m
!         15-04-18  sponge_u et sponge_v etaient annulees par erreur. 
!         05-05-18  Arret urgence occigen
!         09-06-18  viscosite biharmonique hybride: constante & gradient de vitesse
!         23-08-18  suppression flag_nh2d
!         06-09-18  appel external_adve2d contenant 2 schemas d'advection
!         10-01-19  stop compatible avec occigen
! v246    12-02-19  call graph_out avant stop
! v247    24-02-19  wetdry schema 5
! v251    09-04-19  sources d'eau douce sous marine
! v258    19-07-19  Le flux diffusif est nul (si *mask_f(i,j,kmax)) car dans le masque 
!                   la hauteur de la colonne peut etre n'importe quoi et compromettre la CFL
! v269    06-12-19  Ajout velbot_u et velbot_v
! v272    31-12-19  flag_stop
! v294    18-12-20  mask_f(i,j,kmax) 
! v310    19-10-21  if(flag_adve2d==2) then !-c4-biharm-version2-> 
! v311    14-11-21  Cas d'une onde entrant par j=jglb 
! v321    17-01-22  petite modification pour renforcer la robustesse du
!                   schema d'advection
! v346    26-04-22  Empecher ssh de passer sous le fond 
! v347    17-05-22  if(flag_bathy_update==1)call update_bathy 
! v348    11-06-22  flux "no slip" 
! v349   28-06-22   retablissement du point 17-01-22 c.a.d. viscosite donnEe par vitesse advectante
!        29-06-22   if(coastal_viscosity>0.) call dragcoef_no_slip_condition_2d 
! v367   21-02-23   algo fric_t favorisant la continuite (ecrit mais commente car pas vraiment concluant)
!...............................................................................
!    _________                    .__                  .__             ! (�o�) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! m�v�m 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................


contains

!...............................................................................

      subroutine external_mode_driver
      implicit none

      k=kmax !18-12-20

!     flag_externalmode_lf_or_fb='lf' ! Leap-Frog Scheme
      flag_externalmode_lf_or_fb='fb' ! Forward-Backward Scheme

       if(flag_externalmode_lf_or_fb=='fb') then
                    call external_mode_fb
         return
       endif

       if(flag_externalmode_lf_or_fb=='lf') then
                    call external_mode_lf
         return
       endif

      stop 'erreur sur flag_externalmode_lf_or_fb'
      end subroutine external_mode_driver

!...............................................................................

      subroutine external_mode_fb
      implicit none
!     double precision diflxfactor_,grav_over_rho_
      double precision grav_over_rho_

      fbfc2=1.-fbfc1 !08-06-17

      gradssh_u(:,:,0)=gradssh_u(:,:,1)
      gradssh_u(:,:,1)=0.
      gradssh_v(:,:,0)=gradssh_v(:,:,1)
      gradssh_v(:,:,1)=0.

! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit#
!     diflxfactor_=biharm_2dfactor/32.                 !09-06-18
! Note on divise par 32 alors que dans le mode interne on divise par 16 parce que le mode externe
! traite les 2 directions A la fois ce qui divise par 2 la CFL de la viscosite biharmonique

! Note sur diflxfactor_: il permet d'oter le mode numerique dans une direction en
! une iteration si velbar_u=dx/dte_lp si cst_vel_adv=1, d'oter un bruit
! sur les 2 directions si velbar_u=dx/dte_lp et cst_vel_adv=0.5 
! Si on estime que velbar_u est toujours tres inferieur A dx/dte_lp on
! peut amplifier la viscosite avec biharm_2dfactor dont la valeur par
! defaut est 1

      grav_over_rho_=grav/rho

      if(iteration2d_max_now==1) then !---------------->
! No time_splitting case:
       stop 'cas external_mode_no_time_splitting pas teste en FB'
       call external_mode_no_time_splitting ! 04-09-12
       goto 117
      endif                           !---------------->


!$ Compute the external mode equations
      call external_mode_surfacewave

!..........
!$ reset:
      extmodtype='single'

       do j=1,jmax ; do i=1,imax+1
        fluxbar_sumt_u(i,j,0)=fluxbar_sumt_u(i,j,1)
        fluxbar_sumt_u(i,j,1)=0.
       enddo ; enddo
       do j=1,jmax+1 ; do i=1,imax
        fluxbar_sumt_v(i,j,0)=fluxbar_sumt_v(i,j,1)
        fluxbar_sumt_v(i,j,1)=0.
       enddo ; enddo

!..........
!$ Frozen part of the velocity at the firt level above the bottom:
       if(kmax/=1) then ! m[�v�]m -->
! CAS 3D !08-06-17
        do j=1,jmax ; do i=1,imax+1
         velbot3d2d_u(i,j)=                  &
             velbot_u(i,j)                   & !06-12-19
!          0.5*(vel_u(i,j,kmin_u(i,j),0)     & !30-04-17
!              +vel_u(i,j,kmin_u(i,j),1))    &
                       -velbar_u(i,j,1) !20-01-16 
        enddo ; enddo
        do j=1,jmax+1 ; do i=1,imax
         velbot3d2d_v(i,j)=                  &
             velbot_v(i,j)                   & !06-12-19
!          0.5*(vel_v(i,j,kmin_v(i,j),0)     & !30-04-17
!              +vel_v(i,j,kmin_v(i,j),1))    &
                       -velbar_v(i,j,1) 
        enddo ; enddo

       else             ! m[�v�]m -->
! CAS 2D !08-06-17
        velbot3d2d_u=0.  ; velbot3d2d_v=0.
       endif            ! m[�v�]m -->


!-------------------------------------------------------------------------------
!$ Begin the iterative process of the external (2 D) mode:
!     do 9 iteration2d=1,iteration2d_max_now                            !19-04-11
      do 9 iteration2d=iteration2d_begin,iteration2d_max_now                            !19-04-11
!     if(mod(iteration2d,50)==0.and.par%rank==0) &
!       write(6,*)'iteration2d',iteration2d      &
!                ,'temps minutes',real(iteration2d*dte_fw/60.)
!-------------------------------------------------------------------------------

! Lignes suivantes placEes avant un possible calcul de l'advection2d de reference (cas schema anomalie d'advection) !06-09-18
      do j=1,jmax ; do i=1,imax+1
! fluxbar(:,:,0) whom divergence balances ssh_w variation on 2 time steps:
        fluxbar_u(i,j,0)=0.5*(fluxbar_u(i,j,0)+fluxbar_u(i,j,1))
      enddo       ; enddo
      do j=1,jmax+1 ; do i=1,imax
! fluxbar(:,:,0) whom divergence balances ssh_w variation on 2 time steps:
        fluxbar_v(i,j,0)=0.5*(fluxbar_v(i,j,0)+fluxbar_v(i,j,1))
      enddo       ; enddo

!.......................
!$ Frozen forcing terms: (sur 2 iterations pour filtrage du mode numerique)
      if(iteration2d==1) then !ooo>
              if(flag_timesplitting_adve_uv==1) then !m�v�m> !si advection barotrope de reference Ajouter la divergence
                  call external_adve2d                       ! la divergence des flux advecif de l'iteration 1
                  call external_add_adve_ref                 ! au terme adve3d2d (avant procedure filtrage frozenterms)
              endif                                  !m�v�m> 
                        call external_frozen_terms_1 ! Apres possible mise a jour de adve3d2d !06-09-18 !  
      endif                   !ooo>
      if(iteration2d==2)call external_frozen_terms_2 !12-04-13

! If wave/current effect add to the vortex force the increment related to the delta of velbar !12-04-13

! CommentE pour cause d'instabilitE le 22-09-16
!     if(iwve==1) then !--->
!      if(iteration2d>=2)call stokesforces_vortex_external   ! forcement apres call external_frozen_terms_2
!     endif            !--->

!     noise_mode=-noise_mode

!.....

!$ Coriolis terms & advective fluxes of the external mode momentum equations:
      id_ffreq=0
      do j=1,jmax
      do i=1,imax
!$ fstar = Rotation Factor = H*dxdy*f+H*(v*de2/di-u*de1/dj ):
      xy_t(i,j,id_ffreq)=0.25*( & !>>>>> !28-03-16
                  (h_w(i,j)+ssh_w(i,j,1))*(coriolis_t(i,j)*dxdy_t(i,j) &
                       +0.5*( (velbar_v(i  ,j,1)+velbar_v(i  ,j+1,1))  &
                             *(    dy_u(i+1,j  )-    dy_u(i  ,j    ))  &
                             -(velbar_u(i  ,j,1)+velbar_u(i+1,j  ,1))  &
                             *(    dx_v(i  ,j+1)-    dx_v(i  ,j    ))))&
                              )   !>>>>>


      enddo
      enddo

      id_coriolis1=1
      id_coriolis2=2
      do j=1,jmax
      do i=1,imax
!$ u momentum equation:
!$ Coriolis: -H*dxdy*( f + (v*de2/di-u*de1/dj) ) * (v+Vs)
      xy_t(i,j,id_coriolis1)=-xy_t(i,j,id_ffreq)*(velbar_v(i,j,1)+velbar_v(i  ,j+1,1))
      enddo
      enddo

      do j=1,jmax
      do i=1,imax
!$ v momentum equation:
!$ Coriolis: +H*dxdy*( f + (v*de2/di-u*de1/dj) )*(u+Us)
      xy_t(i,j,id_coriolis2)=xy_t(i,j,id_ffreq)*(velbar_u(i,j,1)+velbar_u(i+1,j  ,1))
      enddo
      enddo

! Bottom stress amplitude at t grid nodes:
      if(flag_linearfric==1) then !fricfricfric>
         fric_t=0.5*coef_linearfric !19-01-16
      else                        !fricfricfric>
!       call dragcoef
        do j=1,jmax ; do i=1,imax
! Note: ici (contrairement au cas dans momentum_equation) fric_t est d'avance multiplie par 0.5 pour anticiper la demi-somme
!        fric_t(i,j)=0.5*cdb_t(i,j)*sqrt(0.5*(                     &
!                +(velbot3d2d_u(i+1,j  )+velbar_u(i+1,j  ,1))**2   &
!                +(velbot3d2d_u(i  ,j  )+velbar_u(i  ,j  ,1))**2   &
!                +(velbot3d2d_v(i  ,j+1)+velbar_v(i  ,j+1,1))**2   &
!                +(velbot3d2d_v(i  ,j  )+velbar_v(i  ,j  ,1))**2)) & !18/04/01
!                +small1                                             !27-11-14
! Note: ici (contrairement au cas dans momentum_equation) fric_t est d'avance multiplie par 0.5 pour anticiper la demi-somme
! Cet algo favorise les extremas 
         fric_t(i,j)=0.5*cdb_t(i,j)*sqrt(0.25*( &        !oooo>      !13-05-16
                  (velbot3d2d_u(i+1,j  )+velbar_u(i+1,j  ,1))**2   &
                 +(velbot3d2d_u(i  ,j  )+velbar_u(i  ,j  ,1))**2   &
                 +(velbot3d2d_v(i  ,j+1)+velbar_v(i  ,j+1,1))**2   &
                 +(velbot3d2d_v(i  ,j  )+velbar_v(i  ,j  ,1))**2   &
                 +(velbot3d2d_u(i+1,j  )+velbar_u(i+1,j  ,0))**2   &
                 +(velbot3d2d_u(i  ,j  )+velbar_u(i  ,j  ,0))**2   &
                 +(velbot3d2d_v(i  ,j+1)+velbar_v(i  ,j+1,0))**2   &
                 +(velbot3d2d_v(i  ,j  )+velbar_v(i  ,j  ,0))**2   &
                                                  )  ) & !oooo>
                 +small1                                             !27-11-14
! Note: ici (contrairement au cas dans momentum_equation) fric_t est d'avance multiplie par 0.5 pour anticiper la demi-somme
! Cet algo favorise la continuite: !21-02-23
!       fric_t(i,j)=0.5*cdb_t(i,j)*sqrt(  &
!               (0.25*(                                                &
!                       velbot3d2d_u(i+1,j  )+velbar_u(i+1,j  ,1)      &
!                      +velbot3d2d_u(i  ,j  )+velbar_u(i  ,j  ,1)      &
!                      +velbot3d2d_u(i+1,j  )+velbar_u(i+1,j  ,0)      &
!                      +velbot3d2d_u(i  ,j  )+velbar_u(i  ,j  ,0)      &
!                                                                ))**2 &
!              +(0.25*(                                                &
!                       velbot3d2d_v(i  ,j+1)+velbar_v(i  ,j+1,1)      &
!                      +velbot3d2d_v(i  ,j  )+velbar_v(i  ,j  ,1)      &
!                      +velbot3d2d_v(i  ,j+1)+velbar_v(i  ,j+1,0)      &
!                      +velbot3d2d_v(i  ,j  )+velbar_v(i  ,j  ,0)      &
!                                                                ))**2 &
!                                  )                                   &
!                +small1                                             !27-11-14

        enddo ; enddo
      endif                       !fricfricfric>

      call external_adve2d !xflux, yflux ! 06-09-18 

!     call obc_ext_qdmflux !18-08-15

! toolbox to check conservation properties : !11-11-14
! call external_check_advection_conservation (included in the present file)

!     do j=1,jmax
!     do i=1,imax
!      ssh_ext_w(i,j,1)=0.5*(ssh_w(i,j,1)+ssh_w(i,j,0))
!      ssh_avr_w(i,j,1)=ssh_avr_w(i,j,1)+0.5*(ssh_w(i,j,1)+ssh_w(i,j,0))
!     enddo
!     enddo

!---------------------------------------------------------------------------------------
!     frozenterm3d_u=0.
!     frozenterm2d_u=0.
!     frozenterm2d_v=0.
!     xy_t(:,:,1:2)=0. ! Coriolis
!     xflux_t=0.
!     yflux_f=0.
!     yflux_t=0.
!     xflux_f=0.
!     sponge_u=0.
!     sponge_v=0.
!     wetmask_u=1.
!     fric_t=0.

!$ Begin the Momentum equations:
!     do j=0,jmax+1
!     do i=0,imax+1
!      wetmask_t(i,j)=max(min(                                         &
!           (h_w(i,j)+min(ssh_w(i,j,1),ssh_w(i,j,0))                   &
!                         )/wetdry_cst2-1.                             &
!                             ,1.),0.)
!     enddo
!     enddo


      do j=1,jmax   !12-01-15
      do i=1,imax+1
!      wetmask_u(i,j)=max(min(                                          &
!           (h_u(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i-1,j,1)))/wetdry_cst2-1. &
!                             ,un),zero)                                ! Wetdry attenuation !28/04/09

!      wetmask_u(i,j)=max(min(                                          &
!           (h_u(i,j)+0.5*(                                             &
!         min(ssh_w(i,j,1)+ssh_w(i-1,j,1),ssh_w(i,j,0)+ssh_w(i-1,j,0))  & !13-05-16
!                         ))/wetdry_cst2-1.                             &
!                             ,un),zero)                                ! Wetdry attenuation !28/04/09

!      wetmask_u(i,j)=0.5*(wetmask_t(i,j)+wetmask_t(i-1,j))

! https://docs.google.com/presentation/d/1c361LU5mr_Q6IGKD2ipFJ5qnSvVNygkVajIz1ggZSyo/edit#slide=id.p
! Schema 1:
!      wetmask_u(i,j)=max(min(                                        & !10-08-16
!           (h_u(i,j)+0.5*(min(ssh_w(i,j,1)+ssh_w(i-1,j,1),ssh_w(i,j,0)+ssh_w(i-1,j,0)))) & 
!          /(max(h_u(i,j)-h_w(i,j),h_u(i,j)-h_w(i-1,j))+wetdry_cst2)  &
!                             ,1.),0.)                   
! Schema 2:
!      wetmask_u(i,j)=max(min(                                                                   & !31-10-16
!      (h_u(i,j)-wetdry_cst2+0.5*(min(ssh_w(i,j,1)+ssh_w(i-1,j,1),ssh_w(i,j,0)+ssh_w(i-1,j,0)))) & 
!     /(max(h_u(i,j)-h_w(i,j),h_u(i,j)-h_w(i-1,j))+wetdry_cst2)  &
!                             ,1.),0.)                   
! Schema 3:
!      wetmask_u(i,j)=max(min(                                                                   & !31-10-16
!      (h_u(i,j)-0.1*wetdry_cst2+0.5*(min(ssh_w(i,j,1)+ssh_w(i-1,j,1),ssh_w(i,j,0)+ssh_w(i-1,j,0)))) & 
!     /(max(h_u(i,j)-h_w(i,j),h_u(i,j)-h_w(i-1,j))+wetdry_cst2)  &
!                             ,1.),0.)                   

! https://docs.google.com/presentation/d/1c361LU5mr_Q6IGKD2ipFJ5qnSvVNygkVajIz1ggZSyo/edit#slide=id.p
! Schema 4:
!      wetmask_u(i,j)=max(min(                                          & !11-03-17
!        ( min(h_w(i,j)                   ,h_w(i-1,j)                 ) &
!     +0.5*min(ssh_w(i,j,0)+ssh_w(i-1,j,0),ssh_w(i,j,1)+ssh_w(i-1,j,1)) &
!        )/wetdry_cst2                                                  &
!                             ,1.),0.)                   

! https://docs.google.com/presentation/d/1c361LU5mr_Q6IGKD2ipFJ5qnSvVNygkVajIz1ggZSyo/edit#slide=id.p
! Schema 5:
       wetmask_u(i,j)=max(min(                                          & !24-02-19
         ( min(h_w(i,j)                   ,h_w(i-1,j)                 ) &
      +0.5*max(ssh_w(i,j,0)+ssh_w(i,j,1),ssh_w(i-1,j,0)+ssh_w(i-1,j,1)) &
         )/wetdry_cst2                                                  &
      -0.1 & ! Cette constante pour atteindre 0 avant que la hauteur soit nulle: Dz/wetdry_cst2-0.5=0 pour Dz=0.5*wetdry_cst2
                              ,1.),0.)                   
      enddo
      enddo

!..............................
! Hydrostatic surface pressure gradient force (PGF)
      do j=2,jmax-1
      do i=2,imax

       pgf_u(i,j,0)=pgf_u(i,j,1)

       pgf_u(i,j,1)=                                            & ! Half 2D Pressure Gradient Force U component
        0.5*(h_u(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i-1,j,1)))*(      & ! 0.5*(h+ssh)*[

        grav_over_rho_*( ssh_w(i  ,j,1)*(rho+rhpzavr_w(i  ,j))  & ! g/rh0*d(ssh+rhobar)/di
                        -ssh_w(i-1,j,1)*(rho+rhpzavr_w(i-1,j))) &

                                                         )/dx_u(i,j)   ! ]/dx 
      enddo
      enddo
      do j=2,jmax
      do i=2,imax-1

       pgf_v(i,j,0)=pgf_v(i,j,1)

       pgf_v(i,j,1)=                                            & ! Half 2D Pressure Gradient Force U component
        0.5*(h_v(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i,j-1,1)))*(      & ! 0.5*(h+ssh)*[

        grav_over_rho_*( ssh_w(i  ,j,1)*(rho+rhpzavr_w(i  ,j))  & ! g/rh0*d(ssh+rhobar)/dj
                        -ssh_w(i,j-1,1)*(rho+rhpzavr_w(i,j-1))) &

                                                         )/dy_v(i,j)   ! ]/dy
      enddo
      enddo
!     if(flag_nh2d==1) then !m[0v0]m>
!        call q_add_nhpgf2d  ! Add non-hydrostatic PGF part
!     endif                 !m[0v0]m>

      if(iteration3d==0.and.iteration2d==0) then ! [�v�] >
            pgf_u(:,:,0)=pgf_u(:,:,1)
            pgf_v(:,:,0)=pgf_v(:,:,1)
      endif                                      ! [�v�] >
!..............................

! ANNULER DES TERMES:
!     pgf_u=0.


      do j=2,jmax-1
      do i=2,imax

!$ u component:
         velbar_u(i,j,2)=(                                 &

         velbar_u(i,j,0)*(h_u(i,j)+0.5*(ssh_w(i,j,-1)+ssh_w(i-1,j,-1))) &

! - Flux diffusif temporel "now"
       -(                                                   &
                            fbfc1*( velbar_u(i,j,1 )        & !FB
                                   -velbar_u(i,j,-1))       & !FB
                        +fbfc3*(    velbar_u(i,j,0 )        & ! LP
                                   -velbar_u(i,j,-1))       & ! LP

           )*(h_u(i,j)+0.5*(ssh_w(i,j, 0)+ssh_w(i-1,j, 0))) &

! + Flux diffusif temporel "after"
       +(                                                  &
                           +fbfc1*(                        &
                                   -velbar_u(i,j,0 ))      & !FB !08-06-17
                        +fbfc3*(    velbar_u(i,j,1 )       & ! LP
                                   -velbar_u(i,j,0 ))      & ! LP

           )*(h_u(i,j)+0.5*(ssh_w(i,j, 1)+ssh_w(i-1,j, 1)))            &


        -dte_lp*(           & ! open "(" 1


         +frozenterm3d_u(i,j,1)  & ! le facteur H/dx appliquE en amont dans routine frozen

         +pgf_u(i,j,0)+pgf_u(i,j,1)                                   & ! ssh+rhbar contrib to Pressure Gradient Force

         +( &             !pjppjppjp>

             xy_t(i,j,id_coriolis1)+xy_t(i-1,j,id_coriolis1)          &! Coriolis

            -xflux_t(i,j)+xflux_t(i-1,j)-yflux_f(i,j+1)+yflux_f(i,j)  &! Advection diffusion

          )*invdxdy_u(i,j) & !pjppjppjp>

                 +(fric_t(i,j)+fric_t(i-1,j))*velbot3d2d_u(i,j)        &! Bottom str. (frozen part)
                                           -frozenterm2d_u(i,j,1)      &! Surf+Bot stresses

                 )             & ! close ")" 1

      )/( (fric_t(i,j)+fric_t(i-1,j))*dte_lp                           &!Implicit terms
!        +(h_u(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i-1,j,1)))*wetmask_u(i,j)  &
        +((h_u(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i-1,j,1))-wetdry_cst2)*wetmask_u(i,j)+wetdry_cst2) & ! (H-Hmin)*wetmask+Hmin !08-06-17
         *fbfc2                                                        &!08-06-17
         *(1.+dte_lp*sponge_u(i,j,2)) )                                &

          *mask_u(i,j,kmax)                                            &!Land/sea mask
       *wetmask_u(i,j)


       gradssh_u(i,j,1)=gradssh_u(i,j,1)-dte_lp*pgf_u(i,j,1)*wetmask_u(i,j)

      enddo
      enddo


!$ v component:
      do j=1,jmax+1 !12-01-15
      do i=1,imax
!      wetmask_v(i,j)=max(min(                                          &
!           (h_v(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i,j-1,1)))/wetdry_cst2-1. &
!                             ,un),zero)                                ! Wetdry attenuation !28/04/09

!      wetmask_v(i,j)=max(min(                                         &
!           (h_v(i,j)+0.5*(                                            &
!         min(ssh_w(i,j,1)+ssh_w(i,j-1,1),ssh_w(i,j,0)+ssh_w(i,j-1,0)) & !13-05-16
!                         ))/wetdry_cst2-1.                            &
!                             ,un),zero)                                ! Wetdry attenuation !28/04/09

!      wetmask_v(i,j)=0.5*(wetmask_t(i,j)+wetmask_t(i,j-1))

! https://docs.google.com/presentation/d/1c361LU5mr_Q6IGKD2ipFJ5qnSvVNygkVajIz1ggZSyo/edit#slide=id.p
! Schema 1:
!      wetmask_v(i,j)=max(min(                                        & !10-08-16
!           (h_v(i,j)+0.5*(min(ssh_w(i,j,1)+ssh_w(i,j-1,1),ssh_w(i,j,0)+ssh_w(i,j-1,0)))) & 
!          /(max(h_v(i,j)-h_w(i,j),h_v(i,j)-h_w(i,j-1))+wetdry_cst2)  &
!                             ,1.),0.)                   
! Schema 2:
!      wetmask_v(i,j)=max(min(                                        & !31-10-16
!      (h_v(i,j)-wetdry_cst2+0.5*(min(ssh_w(i,j,1)+ssh_w(i,j-1,1),ssh_w(i,j,0)+ssh_w(i,j-1,0)))) & 
!      /(max(h_v(i,j)-h_w(i,j),h_v(i,j)-h_w(i,j-1))+wetdry_cst2)  &
!                             ,1.),0.)                   

! Schema 3:
!      wetmask_v(i,j)=max(min(                                        & !31-10-16
!      (h_v(i,j)-0.1*wetdry_cst2+0.5*(min(ssh_w(i,j,1)+ssh_w(i,j-1,1),ssh_w(i,j,0)+ssh_w(i,j-1,0)))) & 
!      /(max(h_v(i,j)-h_w(i,j),h_v(i,j)-h_w(i,j-1))+wetdry_cst2)  &
!                             ,1.),0.)                   

! Schema 4:
!      wetmask_v(i,j)=max(min(                                          & !11-03-17
!        ( min(h_w(i,j)                   ,h_w(i,j-1)                 ) &
!     +0.5*min(ssh_w(i,j,0)+ssh_w(i,j-1,0),ssh_w(i,j,1)+ssh_w(i,j-1,1)) &
!        )/wetdry_cst2                                                  &
!                             ,1.),0.)                   

! https://docs.google.com/presentation/d/1c361LU5mr_Q6IGKD2ipFJ5qnSvVNygkVajIz1ggZSyo/edit#slide=id.p
! Schema 5:
       wetmask_v(i,j)=max(min(                                          & !24-02-19
         ( min(h_w(i,j)                   ,h_w(i,j-1)                 ) &
      +0.5*max(ssh_w(i,j,0)+ssh_w(i,j,1),ssh_w(i,j-1,0)+ssh_w(i,j-1,1)) &
         )/wetdry_cst2                                                  &
      -0.1 & ! Cette constante pour atteindre 0 avant que la hauteur soit nulle: Dz/wetdry_cst2-0.5=0 pour Dz=0.5*wetdry_cst2
                              ,1.),0.)                   

      enddo
      enddo

! ANNULER DES TERMES:
!     pgf_v=0.

      do j=2,jmax
      do i=2,imax-1

         velbar_v(i,j,2)=(                                             &


         velbar_v(i,j,0)*(h_v(i,j)+0.5*(ssh_w(i,j,-1)+ssh_w(i,j-1,-1))) &

! - Flux diffusif temporel "now"
       -(                                                   &
                            fbfc1*( velbar_v(i,j,1 )        & !FB
                                   -velbar_v(i,j,-1))       & !FB
                        +fbfc3*(    velbar_v(i,j,0 )        & ! LP
                                   -velbar_v(i,j,-1))       & ! LP

           )*(h_v(i,j)+0.5*(ssh_w(i,j, 0)+ssh_w(i,j-1, 0))) &

! + Flux diffusif temporel "after"
       +(                                                  &
                           +fbfc1*(                        &
                                   -velbar_v(i,j,0 ))      & !FB !08-06-17
                        +fbfc3*(    velbar_v(i,j,1 )       & ! LP
                                   -velbar_v(i,j,0 ))      & ! LP

           )*(h_v(i,j)+0.5*(ssh_w(i,j, 1)+ssh_w(i,j-1, 1)))            &

        -dte_lp*(                                                      &

         +frozenterm3d_v(i,j,1)  & ! le facteur H/dx appliquE en amont dans routine frozen
         +pgf_v(i,j,0)+pgf_v(i,j,1)                                   & ! ssh+rhbar contrib to Pressure Gradient Force

         +( &             !pmpmpmpmpm>  

             xy_t(i,j,id_coriolis2)+xy_t(i,j-1,id_coriolis2)           &!Coriolis

            -yflux_t(i,j)+yflux_t(i,j-1)-xflux_f(i+1,j)+xflux_f(i,j)   &!Advection Diffusion

          )*invdxdy_v(i,j) & !pmpmpmpmpm>

            +(fric_t(i,j)+fric_t(i  ,j-1))*velbot3d2d_v(i,j)           &!Bottom str. (Frozen part)
                                        -frozenterm2d_v(i,j,1)         &!Surf+Bot stresses

                 )                                                     &

      )/(  (fric_t(i,j)+fric_t(i,j-1))*dte_lp                          &
!         +(h_v(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i,j-1,1)))*wetmask_v(i,j) &
         +((h_v(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i,j-1,1))-wetdry_cst2)*wetmask_v(i,j)+wetdry_cst2 ) & ! (H-Hmin)*wetmask+Hmin !08-06-17
          *fbfc2                                                       &
          *(1.+dte_lp*sponge_v(i,j,2)))                                &

          *mask_v(i,j,kmax)                                            &!Sea/land mask
       *wetmask_v(i,j)

       gradssh_v(i,j,1)=gradssh_v(i,j,1)-dte_lp*pgf_v(i,j,1)*wetmask_v(i,j)

!---------------------------------------------------------------------------------------
!$ End of the Momentum equations:
      enddo
      enddo
!---------------------------------------------------------------------------------------

      if(coef_diss_mangrove>0.)call mangrove_2d_friction !23-06-13
      if(coastal_viscosity>0.) call dragcoef_no_slip_condition_2d !29-06-22


! Zonal restoring (baroclinic jet case)
      if(relaxtype_ts==40)call external_zonal_restoring !19-01-16

!     if(nbcanal>0)call webcanals_mpi_2 ! (1) ! canaux envoie velbar_u(2) et ssh_w(1) aux connexions en mer
      if(nbcanal>0)call webcanals_gt1_to_gt2_velbar

!$ (Physical) Open boundary conditions:
!     call obc_ext('openb')
      call obc_ext_fb !13-11-14

!     call external_kelvinwave(0) ! Generateur d'ondes de Kelvin:
! Si place apres call obc_ext_fb faire mpi apres diffusion sinon l'inverse
!     if(flag_nh2d==1)call q_wavebreaking

!$ Barotropic fluxes:
       do j=1,jmax
       do i=1,imax+1
! cette bascule est necessaire A cet emplacement pour la routine NH (sachant que fluxbar(0) a un double emploi)
        fluxbar_u(i,j,0)=fluxbar_u(i,j,1) !11-11-14
       fluxbar_u(i,j,1)=velbar_u(i,j,2)                        &
           *(h_u(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i-1,j  ,1)))     &
           *dy_u(i,j)
      enddo
      enddo

      do j=1,jmax+1
      do i=1,imax
! cette bascule est necessaire A cet emplacement pour la routine NH (sachant que fluxbar(0) a un double emploi)
        fluxbar_v(i,j,0)=fluxbar_v(i,j,1) !11-11-14
       fluxbar_v(i,j,1)=velbar_v(i,j,2)                      &
           *(h_v(i,j)+0.5*(ssh_w(i,j,1)+ssh_w(i,j-1,1)))     &
           *dx_v(i,j)

      enddo
      enddo

!     call obc_ext('river')
      call obc_ext_river    !13-11-14
!     if(nbcanal>0)call webcanals_mpi_2(1) ! canaux envoie fluxbar aux connexions en mer

!****************************************************************

! La condition "radiative" d'onde incidente est calculee dans obc_ext

!     if(iteration2d==1)call q_poisson_ondes_diag ! donne u_amplitude_,periode etc...
!     stop 'coco'
!     if(par%rank==0) then !>>>
! Cas onde incidente i croissant
!      x1=cos(0.-2.*pi*iteration2d*dte_fw/period_-0.5*pi)
!      i=1
!      do j=1,jmax
!         velbar_u(i,j,2)=u_amplitude_*x1
!        fluxbar_u(i,j,1)=velbar_u(i,j,2)  &
!                            *dy_u(i,j)    &
!                            *(h_u(i,j)+ssh_amplitude_*x1)
!      enddo
!     endif                !>>>
!     stop 'cjbpt'
!****************************************************************

!.....
!$ Sea Surface Height:
      do j=1,jmax
      do i=1,imax

!       if(i+par%timax(1)==1.and.j==jmax/2)write(6,*)'ssh_w(i,j,2)',ssh_w(i,j,2)

! attention modifier ce schema implique de changer aussi l'inversion de
! cette equation dans obc_ext_fb
         ssh_w(i,j,2)=                                        &
      max(-h_w(i,j),                                          & !26-04-22
         ssh_w(i,j,1)                                         &
              -dte_fw*(                         & !dtdtdt>
                        ( fluxbar_u(i+1,j  ,1)                &
                         -fluxbar_u(i  ,j  ,1)                &
                         +fluxbar_v(i  ,j+1,1)                &
                         -fluxbar_v(i  ,j  ,1))*invdxdy_t(i,j)   &

                           +omega_w(i,j,kmax+1,1)             & !09-04-15
                           -omega_w(i,j,     1,1)             & !09-04-19
                      )                         & !dtdtdt>      
                       *mask_t(i,j,kmax)) !18-04-19 
! 18-04-19: note sur *mask_t(i,j,kmax)
! Seule la partie entre parenth�se est multipliEe par le masque donc A priori
! le flux d'entree des rivieres n'est pas empeche. La condition limite "river" sur la
! ssh "masquEe" n'est A priori pas empechee non plus puisque demeure alors ssh_w(i,j,2)=ssh_w(i,j,1)

!        +( &
!           (ssh_w(i+1,j,1)-ssh_w(i,j,1))*mask_u(i+1,j,kmax) &
!          +(ssh_w(i-1,j,1)-ssh_w(i,j,1))*mask_u(i  ,j,kmax) &
!          +(ssh_w(i,j+1,1)-ssh_w(i,j,1))*mask_v(i,j+1,kmax) &
!          +(ssh_w(i,j-1,1)-ssh_w(i,j,1))*mask_v(i,j  ,kmax) &
!         )*0.1 

!!     if(isnan(ssh_w(i,j,2)))flag_stop=1 !05-05-18

      enddo
      enddo

       if(nbcanal>0) then
           call webcanals_gt2_to_gt1_ssh
           call webcanals_gt1_to_gt2_ssh
       endif

       call obc_ssh ! this is not a physical boundary condition

!#ifdef 1
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif

!13-04-11

!$ Move forward:
      do j=1,jmax
      do i=1,imax+1
        fluxbar_sumt_u(i,j,1)=fluxbar_sumt_u(i,j,1)+fluxbar_u(i,j,1)
      enddo
      enddo

      do j=0,jmax+1
      do i=0,imax+2
        velbar_u(i,j,-1)=velbar_u(i,j, 0)
        velbar_u(i,j, 0)=velbar_u(i,j, 1)
        velbar_u(i,j, 1)=velbar_u(i,j, 2)
      enddo
      enddo

      do j=1,jmax+1
      do i=1,imax
        fluxbar_sumt_v(i,j,1)=fluxbar_sumt_v(i,j,1)+fluxbar_v(i,j,1)
      enddo
      enddo

      do j=0,jmax+2
      do i=0,imax+1
        velbar_v(i,j,-1)=velbar_v(i,j, 0)
        velbar_v(i,j, 0)=velbar_v(i,j, 1)
        velbar_v(i,j, 1)=velbar_v(i,j, 2)
      enddo
      enddo

      do j=0,jmax+1                !26-11-14
      do i=0,imax+1
        ssh_w(i,j,-1)=ssh_w(i,j,0)
        ssh_w(i,j,0 )=ssh_w(i,j,1)
        ssh_w(i,j,1 )=ssh_w(i,j,2)
      enddo
      enddo

!***************************************************
!     MISE A JOUR DES TABLEAUX AU TEMPS T+1 - FIN
!***************************************************


!-------------------------------------------------------------------------------
!$ End of the iterative process of the external (2 D) mode:
      if(mod(iteration2d,nh2d_graph_period)==0) then
!          call q_u3d_from_u2d
           call graph_out
      endif
!     if(mod(iteration2d,100)==0)call graph_out
!#ifdef 1
!      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
!#endif

    9 continue

!-------------------------------------------------------------------------------

      if(flag_bathy_update==1)call update_bathy !17-05-22

  117 continue

! Valeur de l'elevation de la surface dans le mode interne temps t et t+1
!     x1=0.001*wetdry_cst3
      do j=0,jmax+1                                                   !06-06-09
      do i=0,imax+1
       ssh_int_w(i,j,2)=ssh_w(i,j,1)
      enddo
      enddo

!     if(isnan(ssh_int_w(imax/2,jmax/2,2))) then ! (�v�) ??>
!      call s_stop_by_interrupteur_file('Err isnan in module_external') !04-06-17
!     endif                                      ! (�v�) ??>

      flag_stop=0
      if(isnan(ssh_int_w(imax/2,jmax/2,2)))flag_stop=1
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0) then  !>>> !12-02-19
       call graph_out ; 
       stop 'Err isnan(ssh_int_w(imax/2,jmax/2,2))'
      endif           !>>>

      end subroutine external_mode_fb

!------------------------------------------------------------------------------------

      subroutine external_mode_lf
      implicit none
      double precision diflxfactor_
      stop 'subroutine external_mode_lf discarded'
! Entres autres raisons:
! rhpzavr_w?

      end subroutine external_mode_lf

!------------------------------------------------------------------------------------

      subroutine external_mode_checkmpi
      use module_principal
      use module_parallele                                        !#MPI !16-09-09
      implicit none
!     integer idi_ssh_z0,loop_
       if(flag_1dv==1) then
         call check_mpi_1dv_external !27-11-15
         return
       endif

! Pour verifier la conservation de la parallelisation:
      j1=1 ; j2=jmax
      do i=1,imax
       xy_t(i,j1,1)=ssh_w(i,j1,2)
       xy_t(i,j2,1)=ssh_w(i,j2,2)
      enddo
      i1=1 ; i2=imax
      do j=1,jmax
       xy_t(i1,j,1)=ssh_w(i1,j,2)
       xy_t(i2,j,1)=ssh_w(i2,j,2)
      enddo

      call obc_ssh_checkmpi !30-06-14

      ksecu=0
      j1=1 ; j2=jmax

      if(par%tvoisin(sud)/=mpi_proc_null) then !sssssssssss>
      do i=1,imax
       if(xy_t(i,j1,1)/=ssh_w(i,j1,2))then
        ksecu=1
        write(10+par%rank,*)'mpi non conserv� en par%rank i j =',par%rank,i,j1
        k0=par%tvoisin(sud) ; j0=par%gtjmax(k0,2)-par%gtjmax(k0,1) ! j0=jmax du voisin
        write(10+par%rank,*)'                   voisin sud i j=',par%tvoisin(sud),i,j0-1 !jmax-1
        write(10+par%rank,*)'ssh',xy_t(i,j1,1),ssh_w(i,j1,2),xy_t(i,j1,1)-ssh_w(i,j1,2)

        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'i j global =',i+par%timax(1),j1+par%tjmax(1)
        write(10+par%rank,*)'Mask=',mask_t(i,j1,kmax)
        j=j1
        write(10+par%rank,*)'fluxbar_u(i+1,j,1)=',fluxbar_u(i+1,j,1)
        write(10+par%rank,*)'fluxbar_u(i  ,j,1)=',fluxbar_u(i,j,1)
        write(10+par%rank,*)'fluxbar_v(i,j+1,1)=',fluxbar_v(i,j+1,1)
        write(10+par%rank,*)'fluxbar_v(i,j  ,1)=',fluxbar_v(i,j,1)
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      enddo
      endif                        !sssssssssss>


      if(par%tvoisin(nord)/=mpi_proc_null) then !nnnnnnnnnnn>
      do i=1,imax
       if(xy_t(i,j2,1)/=ssh_w(i,j2,2))then
        ksecu=1
        write(10+par%rank,*)'mpi non conserv� en par%rank i j =' &
         ,par%rank,i,j2
        write(10+par%rank,*)'          voisin en par%rank i j =' &
         ,par%tvoisin(nord),i,2
        write(10+par%rank,*)'ssh',xy_t(i,j2,1),ssh_w(i,j2,2),xy_t(i,j2,1)-ssh_w(i,j2,2)
        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'i j global =',i+par%timax(1),j2+par%tjmax(1)
        write(10+par%rank,*)'Mask=',mask_t(i,j2,kmax)
        j=j2
        write(10+par%rank,*)'fluxbar_u(i+1,j,1)=',fluxbar_u(i+1,j,1)
        write(10+par%rank,*)'fluxbar_u(i  ,j,1)=',fluxbar_u(i,j,1)
        write(10+par%rank,*)'fluxbar_v(i,j+1,1)=',fluxbar_v(i,j+1,1)
        write(10+par%rank,*)'fluxbar_v(i,j  ,1)=',fluxbar_v(i,j,1)
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      enddo
      endif                         !nnnnnnnnnnn>


      i1=1 ; i2=imax

      if(par%tvoisin(ouest)/=mpi_proc_null) then !ooooooooooo>
      do j=1,jmax
       if(xy_t(i1,j,1)/=ssh_w(i1,j,2)) then
        ksecu=1
        write(10+par%rank,*)'mpi non conserv� en par%rank i j =' &
          ,par%rank,i1,j
        k0=par%tvoisin(ouest) ; i0=par%gtimax(k0,2)-par%gtimax(k0,1) ! i0=imax du voisin
        write(10+par%rank,*)'                voisin ouest i j =',par%tvoisin(ouest),i0-1,j ! imax-1,j
        write(10+par%rank,*)'ssh',xy_t(i1,j,1),ssh_w(i1,j,2),xy_t(i1,j,1)-ssh_w(i1,j,2)
        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'i j global =',i1+par%timax(1),j+par%tjmax(1)
        write(10+par%rank,*)'Mask=',mask_t(i1,j,kmax)
        i=i1
        write(10+par%rank,*)'fluxbar_u(i+1,j,1)=',fluxbar_u(i+1,j,1)
        write(10+par%rank,*)'fluxbar_u(i  ,j,1)=',fluxbar_u(i,j,1)
        write(10+par%rank,*)'fluxbar_v(i,j+1,1)=',fluxbar_v(i,j+1,1)
        write(10+par%rank,*)'fluxbar_v(i,j  ,1)=',fluxbar_v(i,j,1)
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      enddo
      endif                          !ooooooooooo>

      if(par%tvoisin(est)/=mpi_proc_null) then   !eeeeeeeeeee>
      do j=1,jmax
       if(xy_t(i2,j,1)/=ssh_w(i2,j,2)) then
        ksecu=1
        write(10+par%rank,*)'mpi non conserv� en par%rank i j =' &
          ,par%rank,i2,j
        write(10+par%rank,*)'          voisin en par%rank i j =' &
          ,par%tvoisin(est),2,j
        write(10+par%rank,*)'ssh',xy_t(i2,j,1),ssh_w(i2,j,2),xy_t(i2,j,1)-ssh_w(i2,j,2)
        write(10+par%rank,*)'iteration3d & 2d=',iteration3d,iteration2d!23-04-14
        write(10+par%rank,*)'i j global =',i2+par%timax(1),j+par%tjmax(1)
        write(10+par%rank,*)'Mask=',mask_t(i2,j,kmax)
        i=i2
        write(10+par%rank,*)'fluxbar_u(i+1,j,1)=',fluxbar_u(i+1,j,1)
        write(10+par%rank,*)'fluxbar_u(i  ,j,1)=',fluxbar_u(i,j,1)
        write(10+par%rank,*)'fluxbar_v(i,j+1,1)=',fluxbar_v(i,j+1,1)
        write(10+par%rank,*)'fluxbar_v(i,j  ,1)=',fluxbar_v(i,j,1)
        write(10+par%rank,*)'-----------------------------------------'!23-04-14
       endif
      enddo
      endif                          !eeeeeeeeeee>

      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
      flag_stop=0
      if(ksecu==1)then
      write(6,*)'Erreurs mpi notifiees dans fichiers fort locaux'
!       stop 'STOP dans external_mode_checkmpi'
      flag_stop=1
      endif
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr)
      if(k0/=0)stop 'STOP dans external_mode_checkmpi' !31-12-19

      end subroutine external_mode_checkmpi

!..........................................................................

      subroutine external_mode_no_time_splitting ! 04-09-12
      use module_principal
      implicit none

       do j=1,jmax
       do i=1,imax+1
        fluxbar_u(i,j,1)=0.
       enddo
       enddo
       do k=1,kmax
       do j=1,jmax
       do i=1,imax+1
        fluxbar_u(i,j,1)=fluxbar_u(i,j,1)+veldydz_u(i,j,k,1)
       enddo
       enddo
       enddo

       do j=1,jmax+1
       do i=1,imax
        fluxbar_v(i,j,1)=0.
       enddo
       enddo
       do k=1,kmax
       do j=1,jmax+1
       do i=1,imax
        fluxbar_v(i,j,1)=fluxbar_v(i,j,1)+veldxdz_v(i,j,k,1)
       enddo
       enddo
       enddo

      end subroutine external_mode_no_time_splitting ! 04-09-12

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

      subroutine external_kelvinwave(case_)
      use module_principal
      use module_parallele
      implicit none
      double precision c0_,radius_,time_,x_,y_,omega_ &
                      ,kx_,ssh_,ssh0_,period_
      integer case_


! Generateur d'ondes de Kelvin:
        period_=0.1*43200. ! 12h
        ssh0_=0.01

        if(par%rank==0.and.iteration2d==1)write(6,*) &
        'bidouille appel a kelvinwave'


!-----------------------------------------------
! CAS FORCAGE A LA FRONTIERE i=1

      if(case_==0) then !0000000000000>

      if(par%timax(1)==0) then !------>
        i=1  ! OBC i=1
!       time_=elapsedtime_now+ iteration2d   *0.5*dte_lp    ! time now
        time_=elapsedtime_now+(iteration2d+1)*0.5*dte_lp    ! time after
        omega_=2.*pi/period_

         do j=2,jmax-1

          c0_=sqrt(grav*h_u(i,j))
          radius_=c0_/coriolis_t(i,j)
          kx_=omega_/c0_
          y_=real(j-1)*dyb
          x_=real(i-1)*dxb

! Solution analytique pour ssh:
          ssh_=ssh0_*exp(-coriolis_t(i,j)*y_/c0_)  &
                    *sin(kx_*x_-omega_*time_)

          velbar_u(i,j,after)=grav/c0_*ssh_

         enddo
      endif                    !------>

      endif             !0000000000000>
!-----------------------------------------------



!-----------------------------------------------
! CAS SOLUTION ANALYTIQUE POUR ARCHIVAGE DANS FICHIER NETCDF
      if(case_==1) then !11111111111111>

!        time_=elapsedtime_now-0.5*dte_lp
         time_=elapsedtime_now
         omega_=2.*pi/period_

         write(6,*)'time_=',time_

         do j1=1,jglb
         do i1=1,iglb
         i=i1-par%timax(1)
         j=j1-par%tjmax(1)
         if(i>=1.and.i<=imax+1.and.j>=0.and.j<=jmax+1)then !mpimpimpi>

          c0_=sqrt(grav*h_u(i,j))
          radius_=c0_/coriolis_t(i,j)
          kx_=omega_/c0_
          y_=real(j1-2)*dyb ! y=0 pour le premier point en mer
          x_=real(i1-1)*dxb ! x=0 sur le point de C.L. pour u

! Solution analytique pour ssh:
          ssh_=ssh0_*exp(-coriolis_t(i,j)*y_/c0_)  &
                    *sin(kx_*x_-omega_*time_)

          velbarobc_u(i,j,1)=grav/c0_*ssh_
             sshobc_w(i,j,1)=ssh_

! c'est velbar(1) qui doit correspondre � velbarobc...
!         if(j==10)write(10+par%rank,*)i1,velbarobc_u(i,j,1),velbar_u(i,j,1),vel_u(i,j,kmax,1)
!         if(j==10)write(10+par%rank,*)i1,velbarobc_u(i,j,1),velbar_u(i,j,0:1)

         endif                                             !mpimpimpi>
         enddo
         enddo


      endif             !11111111111111>
!-----------------------------------------------

      end subroutine external_kelvinwave

!...............................................................................

      subroutine external_frozen_terms_1 !12-04-13
      implicit none

! Iteration2d=1 ==> frozenterms(1) = 0.5*(frozenterms(t-1)+frozenterms(t))
! Iteration2d>1 ==> frozenterms(1) = frozenterms(t)

!.......................
!$ Frozen forcing terms:
      do j=2,jmax-1
      do i=2,imax

!$ u momentum equations:
      frozenterm3d_u(i,j,0)=  &

          adve3d2d_u(i,j)*wetmask_u(i,j)     & !  -momentum advection !03-09-18

      +( &                       !pmxpmx>

                -sponge_u(i,j,2)*velbarobc_u(i,j,1)*dx_u(i,j) & 

            +( pres3d2d_u(i,j)                                & !  -pressure: density, air, tide, stokes
!!!!          +adve3d2d_u(i,j)                                & !  -momentum advection
      -stokesforces3d2d_u(i,j) ) *wetmask_u(i,j)              & !  -Stokes Vortex

       )*hz_u(i,j,1)/dx_u(i,j)   !pmxpmx>

!     if(i==imax/2.and.j==jmax/2)write(10+par%rank,*)adve3d2d_u(i,j)

!........

      frozenterm3d_u(i,j,1)=0.5*( frozenterm3d_u(i,j,1)        &
                                 +frozenterm3d_u(i,j,0))


      frozenterm2d_u(i,j,0)=                                  &
                               ( wstress_u(i,j,1)             &
                                                  )/rho
      frozenterm2d_u(i,j,1)=0.5*( frozenterm2d_u(i,j,1)        &
                                 +frozenterm2d_u(i,j,0))

      enddo
      enddo

      do j=2,jmax
      do i=2,imax-1

!$ v momentum equations:
      frozenterm3d_v(i,j,0)=  &

          adve3d2d_v(i,j)*wetmask_v(i,j)     & !  -momentum advection !03-09-18

      +( &                       !pmxpmx>

                -sponge_v(i,j,2)*velbarobc_v(i,j,1)*dy_v(i,j) & 

            +( pres3d2d_v(i,j)                                & !  -pressure: density, air, tide, stokes
!!!!          +adve3d2d_v(i,j)                                & !  -momentum advection
      -stokesforces3d2d_v(i,j) ) *wetmask_v(i,j)              & !  -Stokes Vortex

       )*hz_v(i,j,1)/dy_v(i,j)   !pmxpmx>

!........

      frozenterm3d_v(i,j,1)=0.5*( frozenterm3d_v(i,j,1)       &
                                 +frozenterm3d_v(i,j,0))


      frozenterm2d_v(i,j,0)=                                  &
                               ( wstress_v(i,j,1)             &
                                                  )/rho
      frozenterm2d_v(i,j,1)=0.5*( frozenterm2d_v(i,j,1)       &
                                 +frozenterm2d_v(i,j,0))

      enddo
      enddo
!...........

! Very first iteration (quasi initial state)
      if(iteration3d==0)call external_frozen_terms_2

      end subroutine external_frozen_terms_1

!...............................................................................

      subroutine external_frozen_terms_2 !12-04-13
      implicit none

! Iteration2d=1 ==> frozenterms(1) = 0.5*(frozenterms(t-1)+frozenterms(t))
! Iteration2d>1 ==> frozenterms(1) = frozenterms(t)

!.......................
!$ Frozen forcing terms:
      do j=2,jmax-1
      do i=2,imax
       frozenterm3d_u(i,j,1)=frozenterm3d_u(i,j,0)
       frozenterm2d_u(i,j,1)=frozenterm2d_u(i,j,0)
      enddo
      enddo

      do j=2,jmax
      do i=2,imax-1
       frozenterm3d_v(i,j,1)=frozenterm3d_v(i,j,0)
       frozenterm2d_v(i,j,1)=frozenterm2d_v(i,j,0)
      enddo
      enddo
!...........


      end subroutine external_frozen_terms_2

!...............................................................................
!...............................................................................

      subroutine external_zonal_restoring
      implicit none
      integer moduliter_

      if(relax_ts<=0.)return

!     stop 'external_zonal_restoring'
!...................................................................
! Rappel moyenne zonale cas testjet barocline comodo !17-10-13!15-01-16
!...................................................................
! Attention les dimensions de ces tableaux sont communes avec celles du
! rappel de T
!
      moduliter_=iteration2d_max_now
      if(mod(iteration2d,moduliter_)/=0) return
      
      if(.not.allocated(iaverag1d_in ))allocate(iaverag1d_in (jglb,1))
      if(.not.allocated(iaverag1d_out))allocate(iaverag1d_out(jglb,1))

! u component:
!     velbarobc_u=0.5 ; velbar_u=-0.5 ! Pour verif u2d 
      x2=timeweightobc(vel_id) ; x0=1.-x2
      iaverag1d_in=0. 
      do j=2,jmax-1 
       j1=j+par%tjmax(1)
       do i=2,imax
        iaverag1d_in(j1,1)=iaverag1d_in(j1,1)+mask_i_u(i)*mask_j_u(j)*(&!ooo>!20-01-16
                           x2*velbarobc_u(i,j,2)                       &
                          +x0*velbarobc_u(i,j,0)                       &
                                -velbar_u(i,j,0)                       &
                                                                      ) !ooo>

       enddo
      enddo
      call mpi_allreduce( iaverag1d_in         &
                         ,iaverag1d_out        &
                         ,jglb                 &
                         ,mpi_double_precision & 
                         ,mpi_sum              &
                         ,par%comm2d ,ierr)


      x0=dte_fw*moduliter_*relax_ts/real(iglb-2)

      do j=2,jmax-1 ; do i=2,imax
      velbar_u(i,j,2)=velbar_u(i,j,2)+iaverag1d_out(j+par%tjmax(1),1)*x0
!     write(10+par%rank,*),iaverag1d_out(j+par%tjmax(1),1)/real(iglb-2) ! Pour verif u2d. Doit donner 1
      enddo ; enddo
!     stop 'verif u2d'

! v component:
!     velbarobc_v=0.5 ; velbar_v=-0.5 ! Pour verif v2d 
      x2=timeweightobc(vel_id) ; x0=1.-x2
      iaverag1d_in=0.
      do j=2,jmax
       j1=j+par%tjmax(1)
       do i=2,imax-1
        iaverag1d_in(j1,1)=iaverag1d_in(j1,1)+mask_i_v(i)*mask_j_v(j)*(&!ooo>
                           x2*velbarobc_v(i,j,2)                       &
                          +x0*velbarobc_v(i,j,0)                       &
                                -velbar_v(i,j,0)                       &
                                                                      ) !ooo>
       enddo
      enddo
      call mpi_allreduce( iaverag1d_in         &
                         ,iaverag1d_out        &
                         ,jglb                 &
                         ,mpi_double_precision & 
                         ,mpi_sum              &
                         ,par%comm2d ,ierr)

      x0=dte_fw*moduliter_*relax_ts/real(iglb-2)
      do j=2,jmax ; do i=2,imax-1
      velbar_v(i,j,2)=velbar_v(i,j,2)+iaverag1d_out(j+par%tjmax(1),1)*x0
!     write(20+par%rank,*),iaverag1d_out(j+par%tjmax(1),1)/real(iglb-2) ! Pour verif v2d. Doit donner 1
      enddo ; enddo
!     stop 'verif v2d'

!...................................................................
! Rappel moyenne zonale cas testjet barocline comodo !17-10-13
!...................................................................

      end subroutine external_zonal_restoring


!...............................................................................

      subroutine external_mode_surfacewave !09-04-16
      implicit none






      end subroutine external_mode_surfacewave

!...............................................................................

      subroutine external_adve2d
      implicit none
      double precision diflxfactor_



      if(flag_adve2d==0) then !-c4-biharm-version0-> !19-10-21!17-01-22
!//////////////////////////////////////////////////////////////////////////////////////////
! Advection centrEe O4 + Diffusion biharmonique

! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit#
      diflxfactor_=biharm_2dfactor/32.                 !09-06-18
! Note on divise par 32 alors que dans le mode interne on divise par 16 parce que le mode externe
! traite les 2 directions A la fois ce qui divise par 2 la CFL de la viscosite biharmonique

! Note sur diflxfactor_: il permet d'oter le mode numerique dans une direction en
! une iteration si velbar_u=dx/dte_lp si cst_vel_adv=1, d'oter un bruit
! sur les 2 directions si velbar_u=dx/dte_lp et cst_vel_adv=0.5 
! Si on estime que velbar_u est toujours tres inferieur A dx/dte_lp on
! peut amplifier la viscosite avec biharm_2dfactor dont la valeur par
! defaut est 1

!............
! u equation:
      do j=2,jmax-1
      do i=1,imax

      xflux_t(i,j)=                                               &
! Ox Advective flux at t point:
                -0.125*(fluxbar_u(i,j,0) +fluxbar_u(i+1,j,0))     & !11-11-14
!                      *(velbar_u(i,j,now)+velbar_u(i+1,j,now))   &
                 *( (velbar_u(i  ,j,now)+velbar_u(i+1,j,now)         &
                    +velbar_u(i  ,j,0  )+velbar_u(i+1,j,0  ))*1.125  & !30-01-16
                   -(velbar_u(i-1,j,now)+velbar_u(i+2,j,now)         &
                    +velbar_u(i-1,j,0  )+velbar_u(i+2,j,0  ))*0.125) &
! Ox Diffusive flux:
!......................................................................
!     +diflxfactor_*abs(dy_t(i,j)*(h_w(i,j)+ssh_w(i,j,1)))*( & !m�v�m>
!          biharm_e1*dx_t(i,j)*inv_dte_lp                        & ! biviscosite constante 09-06-18
!         +biharm_e2*0.5*abs(velbar_u(i,j,0)-velbar_u(i+1,j,0))  & ! biviscosite gradient de vitesse
!                                                          ) & !m�v�m>
      +diflxfactor_*0.5*abs(fluxbar_u(i,j,0)-fluxbar_u(i+1,j,0)) & !17-01-22
! Note: le mode externe traite les 2 directions simultanement ce qui rend 2 fois
! plus severe le critere CFL et explique que la viscosite est divisee par 2 par rapport
! a celle du mode interne.
!......................................................................
      *( & !ooo>
!         (velbar_u(i+1,j  ,before )-velbar_u(i  ,j  ,before ))*3. &
!        -(velbar_u(i+2,j  ,before )-velbar_u(i-1,j  ,before ))    &
!           *mask_u(i+2,j  ,kmax)     *mask_u(i-1,j  ,kmax)        & !22-07-17
          (velbar_u(i+1,j  ,before)-velbar_u(i  ,j  ,before)      &
          +velbar_u(i+1,j  ,now   )-velbar_u(i  ,j  ,now   ))*1.5 &
         -(velbar_u(i+2,j  ,before)-velbar_u(i-1,j  ,before)      &
          +velbar_u(i+2,j  ,now   )-velbar_u(i-1,j  ,now   ))*0.5 &
            *mask_u(i+2,j  ,kmax)    *mask_u(i-1,j  ,kmax)        & !22-07-17
       )   !ooo>
      enddo
      enddo

!............
! u equation:
      do j=2,jmax
      do i=2,imax
   
      yflux_f(i,j)=                                                &
! Oy Advective flux at p point:
                -0.125*(fluxbar_v(i-1,j  ,0) +fluxbar_v(i,j,0))    &
!                      *(velbar_u(i  ,j-1,now)+velbar_u(i,j,now))  &
               *(  (velbar_u(i  ,j-1,now)+velbar_u(i,j  ,now)         &
                   +velbar_u(i  ,j-1,0  )+velbar_u(i,j  ,0  ))*1.125  &
                  -(velbar_u(i  ,j-2,now)+velbar_u(i,j+1,now)         &
                   +velbar_u(i  ,j-2,0  )+velbar_u(i,j+1,0  ))*0.125) &
! Oy Diffusive flux:
!......................................................................
!       +diflxfactor_*mask_f(i,j,kmax) & !voir notes !19-07-19!18-12-20
!                  *abs(dx_f(i,j)                                      &
!                      *(h_f(i,j)+0.25*(ssh_w(i  ,j  ,1)               &! schema cst   (et sa difinition de +diflxfactor_)
!                                      +ssh_w(i-1,j  ,1)               &
!                                      +ssh_w(i  ,j-1,1)               &
!                                      +ssh_w(i-1,j-1,1))))*( & !m�v�m>
!            biharm_e1*dy_f(i,j)*inv_dte_lp                       & ! biviscosite constante 09-06-18
!           +biharm_e2*0.5*abs(velbar_u(i,j-1,0)-velbar_u(i,j,0)) & ! biviscosite gradient de vitesse
!                                                           ) & !m�v�m>
        +diflxfactor_*mask_f(i,j,kmax)*0.5*abs(fluxbar_v(i-1,j  ,0)-fluxbar_v(i,j,0)) & !17-01-22
!                        +diflxfactor_*0.5*abs(fluxbar_u(i  ,j-1,0)-fluxbar_u(i,j,0)) & !11-06-22
!      *(1+flag_free_slip*(mask_f(i,j,kmax)-1)) & ! free slip si flag_free_slip=1      !11-06-22 
!......................................................................
      *( & !ooo>
!        (velbar_u(i  ,j  ,before )-velbar_u(i  ,j-1,before ))*3. &
!       -(velbar_u(i  ,j+1,before )-velbar_u(i  ,j-2,before ))    &
!          *mask_u(i  ,j+1,kmax   )  *mask_u(i  ,j-2,kmax)        & !22-07-17
         (velbar_u(i  ,j  ,before )-velbar_u(i  ,j-1,before )      &
         +velbar_u(i  ,j  ,now    )-velbar_u(i  ,j-1,now    ))*1.5 &
        -(velbar_u(i  ,j+1,before )-velbar_u(i  ,j-2,before )      &
         +velbar_u(i  ,j+1,now    )-velbar_u(i  ,j-2,now    ))*0.5 &
           *mask_u(i  ,j+1,kmax   )  *mask_u(i  ,j-2,kmax)        & !22-07-17
       )   !ooo>
      enddo
      enddo
! Notes: le flux advectif de yflux_f(i,j) est nul si mask_t(i-1,j)=mask_t(i,j)=0 (car alors fluxbar_v(i-1:i,j)=0)
! Le flux diffusif est nul (si *mask_f(i,j,kmax)) car dans le masque la hauteur de la colonne peut etre
! n'importe quoi et compromettre la CFL
! Par consequent cela correspond A une condition "free-slip" A la cOte. 


!     if(iteration2d==4.and.mod(iteration3d,10)==0) then
!      i=imax/2 ; j=jmax/2
!      if(mask_t(i,j,kmax)==1) &
!      write(10+par%rank,*) &
!     ,abs(dy_t(i,j)*(h_w(i,j)+ssh_w(i,j,1)))*( & !m�v�m>
!          biharm_e1*dx_t(i,j)*inv_dte_lp                        & ! biviscosite constante 09-06-18
!         +biharm_e2*0.5*abs(velbar_u(i,j,0)-velbar_u(i+1,j,0))  & ! biviscosite gradient de vitesse
!                                                          ) & !m�v�m>
!     ,diflxfactor_*0.5*abs(fluxbar_u(i,j,0) +fluxbar_u(i+1,j,0)) &
!     ,' trouvemi'
!     endif

!............
! v equation:
      do j=2,jmax
      do i=2,imax

! boucle stricte: do i=2,imax | do j=2,jmax
      xflux_f(i,j)=                                                &
! Ox Advective flux at p point:
                -0.125*(fluxbar_u(i  ,j-1,0) +fluxbar_u(i,j,0))    &
!                      *(velbar_v(i-1,j  ,now)+velbar_v(i,j,now))  &
                 *( (velbar_v(i-1,j  ,now)+velbar_v(i  ,j,now)         &
                    +velbar_v(i-1,j  ,0  )+velbar_v(i  ,j,0  ))*1.125  &
                   -(velbar_v(i-2,j  ,now)+velbar_v(i+1,j,now)         &
                    +velbar_v(i-2,j  ,0  )+velbar_v(i+1,j,0  ))*0.125) &
! Ox Diffusive flux at p point:
!..........................................................................
!       +diflxfactor_*mask_f(i,j,kmax) & !voir notes 19-07-19 !18-12-20
!                  *abs(dy_f(i,j)                                      &
!                      *(h_f(i,j)+0.25*(ssh_w(i  ,j  ,1)               &! schema cst   (et sa difinition de +diflxfactor_)
!                                      +ssh_w(i-1,j  ,1)               &
!                                      +ssh_w(i  ,j-1,1)               &
!                                      +ssh_w(i-1,j-1,1))))*( & !m�v�m>
!            biharm_e1*dx_f(i,j)*inv_dte_lp                       & ! biviscosite constante 09-06-18
!           +biharm_e2*0.5*abs(velbar_v(i,j,0)-velbar_v(i-1,j,0)) & ! biviscosite gradient de vitesse
!                                                           ) & !m�v�m>
        +diflxfactor_*mask_f(i,j,kmax)*0.5*abs(fluxbar_u(i  ,j-1,0)-fluxbar_u(i,j,0)) & !17-01-22
!                        +diflxfactor_*0.5*abs(fluxbar_v(i-1,j  ,0)-fluxbar_v(i,j,0)) & !11-06-22
!      *(1+flag_free_slip*(mask_f(i,j,kmax)-1)) & ! free slip si flag_free_slip=1       !11-06-22 
!..........................................................................
      *( & !ooo>
!        (velbar_v(i  ,j  ,before )-velbar_v(i-1,j  ,before ))*3. &
!       -(velbar_v(i+1,j  ,before )-velbar_v(i-2,j  ,before ))    &
!          *mask_v(i+1,j  ,kmax)     *mask_v(i-2,j  ,kmax)        & !22-07-17
         (velbar_v(i  ,j  ,before )-velbar_v(i-1,j  ,before )      &
         +velbar_v(i  ,j  ,now    )-velbar_v(i-1,j  ,now    ))*1.5 &
        -(velbar_v(i+1,j  ,before )-velbar_v(i-2,j  ,before )      &
         +velbar_v(i+1,j  ,now    )-velbar_v(i-2,j  ,now    ))*0.5 &
           *mask_v(i+1,j  ,kmax)     *mask_v(i-2,j  ,kmax)        & !22-07-17
       )   !ooo>
      enddo
      enddo
! Note le flux advectif de xflux_f(i,j) est nul si mask_t(i,j-1)=mask_t(i,j)=0 (car alors fluxbar_u(i,j-1:j)=0)
! Le flux diffusif est nul (si *mask_f(i,j,kmax)) car dans le masque la hauteur de la colonne peut etre
! n'importe quoi et compromettre la CFL
! Par consequent cela correspond A une condition "free-slip" A la cOte. 

!............
! v equation:
      do j=1,jmax
      do i=2,imax-1
      yflux_t(i,j)=                                                &
! Oy Advective flux at t point:
                -0.125*(fluxbar_v(i,j,0) +fluxbar_v(i,j+1,0))      &
!                      *(velbar_v(i,j,now)+velbar_v(i,j+1,now))    &
                 *( (velbar_v(i,j  ,now)+velbar_v(i,j+1,now)         &
                    +velbar_v(i,j  ,0  )+velbar_v(i,j+1,0  ))*1.125  &
                   -(velbar_v(i,j-1,now)+velbar_v(i,j+2,now)         &
                    +velbar_v(i,j-1,0  )+velbar_v(i,j+2,0  ))*0.125) &
! Oy Diffusive flux at t point:
!.......................................................................
!     +diflxfactor_*abs(dx_t(i,j)*(h_w(i,j)+ssh_w(i,j,1)))*( & !m�v�m> 
!             biharm_e1*dy_t(i,j)*inv_dte_lp                       & ! biviscosite constante 09-06-18
!            +biharm_e2*0.5*abs(velbar_v(i,j,0)-velbar_v(i,j+1,0)) & ! biviscosite gradient de vitesse 
!                                                          ) & !m�v�m>
      +diflxfactor_*0.5*abs(fluxbar_v(i,j,0)-fluxbar_v(i,j+1,0)) & !17-01-22
!.......................................................................
      *( & !ooo>
!        (velbar_v(i  ,j+1,before )-velbar_v(i  ,j  ,before ))*3. &
!       -(velbar_v(i  ,j+2,before )-velbar_v(i  ,j-1,before ))    &
!          *mask_v(i  ,j+2,kmax)     *mask_v(i  ,j-1,kmax)        & !22-07-17
         (velbar_v(i  ,j+1,before )-velbar_v(i  ,j  ,before )      &
         +velbar_v(i  ,j+1,now    )-velbar_v(i  ,j  ,now    ))*1.5 &
        -(velbar_v(i  ,j+2,before )-velbar_v(i  ,j-1,before )      &
         +velbar_v(i  ,j+2,now    )-velbar_v(i  ,j-1,now    ))*0.5 &
           *mask_v(i  ,j+2,kmax)     *mask_v(i  ,j-1,kmax)        & !22-07-17
       )   !ooo>

      enddo
      enddo

      return
      endif                   !-c4-biharm-version0->

      if(flag_adve2d==1) then !-up2->
!//////////////////////////////////////////////////////////////////////////////////////////
! Schema UPWIND

!............
! u equation:
      do j=2,jmax-1
      do i=1,imax
       xflux_t(i,j)=  &
              -0.25*( & !m�v�m>
                        (fluxbar_u(i,j,0)+fluxbar_u(i+1,j,0))     & 
                        *(velbar_u(i,j,0)+ velbar_u(i+1,j,0))     &
                    +abs(fluxbar_u(i,j,0)+fluxbar_u(i+1,j,0))     & 
                        *(velbar_u(i,j,0)- velbar_u(i+1,j,0))     &
                    )   !m�v�m>
      enddo
      enddo

!............
! u equation:
      do j=2,jmax
      do i=2,imax
       yflux_f(i,j)=   &
              -0.25*(  & !m�v�m>
                        (fluxbar_v(i-1,j  ,0)+fluxbar_v(i,j,0))  &
                        *(velbar_u(i  ,j-1,0)+ velbar_u(i,j,0))  &
                    +abs(fluxbar_v(i-1,j  ,0)+fluxbar_v(i,j,0))  &
                        *(velbar_u(i  ,j-1,0)- velbar_u(i,j,0))  &
                    )   !m�v�m>
      enddo
      enddo

!............
! v equation:
      do j=2,jmax
      do i=2,imax
       xflux_f(i,j)=   &
              -0.25*(  & !m�v�m>
                        (fluxbar_u(i  ,j-1,0)+fluxbar_u(i,j,0))  &
                        *(velbar_v(i-1,j  ,0)+ velbar_v(i,j,0))  &
                    +abs(fluxbar_u(i  ,j-1,0)+fluxbar_u(i,j,0))  &
                        *(velbar_v(i-1,j  ,0)- velbar_v(i,j,0))  &
                     )   !m�v�m>
      enddo
      enddo

!............
! v equation:
      do j=1,jmax
      do i=2,imax-1
       yflux_t(i,j)=   &
              -0.25*(  & !m�v�m>
                       (fluxbar_v(i,j,0)+fluxbar_v(i,j+1,0))   &
                       *(velbar_v(i,j,0)+ velbar_v(i,j+1,0))   &
                   +abs(fluxbar_v(i,j,0)+fluxbar_v(i,j+1,0))   &
                       *(velbar_v(i,j,0)- velbar_v(i,j+1,0))   &
                    )   !m�v�m>
      enddo
      enddo

      endif                   !-up2->


      if(flag_adve2d==2) then !-c4-biharm-version2-> !19-10-21!17-01-22
!//////////////////////////////////////////////////////////////////////////////////////////
! Advection centrEe O4 + Diffusion biharmonique

! https://docs.google.com/document/d/1VtN0Tp0FfiqZyr_L01QYtsZa47CmkMYnv6nizqnZhXc/edit#
      diflxfactor_=biharm_2dfactor/32.                 !09-06-18
! Note on divise par 32 alors que dans le mode interne on divise par 16 parce que le mode externe
! traite les 2 directions A la fois ce qui divise par 2 la CFL de la viscosite biharmonique

! Note sur diflxfactor_: il permet d'oter le mode numerique dans une direction en
! une iteration si velbar_u=dx/dte_lp si cst_vel_adv=1, d'oter un bruit
! sur les 2 directions si velbar_u=dx/dte_lp et cst_vel_adv=0.5 
! Si on estime que velbar_u est toujours tres inferieur A dx/dte_lp on
! peut amplifier la viscosite avec biharm_2dfactor dont la valeur par
! defaut est 1

!............
! u equation:
      do j=2,jmax-1
      do i=1,imax

      xflux_t(i,j)=                                               &
! Ox Advective flux at t point:
                -0.125*(fluxbar_u(i,j,0) +fluxbar_u(i+1,j,0))     & !11-11-14
!                      *(velbar_u(i,j,now)+velbar_u(i+1,j,now))   &
                 *( (velbar_u(i  ,j,now)+velbar_u(i+1,j,now)         &
                    +velbar_u(i  ,j,0  )+velbar_u(i+1,j,0  ))*1.125  & !30-01-16
                   -(velbar_u(i-1,j,now)+velbar_u(i+2,j,now)         &
                    +velbar_u(i-1,j,0  )+velbar_u(i+2,j,0  ))*0.125) &
! Ox Diffusive flux:
!......................................................................
!     +diflxfactor_*abs(dy_t(i,j)*(h_w(i,j)+ssh_w(i,j,1)))*( & !m�v�m>
!          biharm_e1*dx_t(i,j)*inv_dte_lp                        & ! biviscosite constante 09-06-18
!         +biharm_e2*0.5*abs(velbar_u(i,j,0)-velbar_u(i+1,j,0))  & ! biviscosite gradient de vitesse
!                                                          ) & !m�v�m>
      +diflxfactor_*0.5*(abs(fluxbar_u(i,j,0))+abs(fluxbar_u(i+1,j,0))) & !17-01-22
! Note: le mode externe traite les 2 directions simultanement ce qui rend 2 fois
! plus severe le critere CFL et explique que la viscosite est divisee par 2 par rapport
! a celle du mode interne.
!......................................................................
      *( & !ooo>
!         (velbar_u(i+1,j  ,before )-velbar_u(i  ,j  ,before ))*3. &
!        -(velbar_u(i+2,j  ,before )-velbar_u(i-1,j  ,before ))    &
!           *mask_u(i+2,j  ,kmax)     *mask_u(i-1,j  ,kmax)        & !22-07-17
          (velbar_u(i+1,j  ,before)-velbar_u(i  ,j  ,before)      &
          +velbar_u(i+1,j  ,now   )-velbar_u(i  ,j  ,now   ))*1.5 &
         -(velbar_u(i+2,j  ,before)-velbar_u(i-1,j  ,before)      &
          +velbar_u(i+2,j  ,now   )-velbar_u(i-1,j  ,now   ))*0.5 &
            *mask_u(i+2,j  ,kmax)    *mask_u(i-1,j  ,kmax)        & !22-07-17
       )   !ooo>
      enddo
      enddo

!............
! u equation:
      do j=2,jmax
      do i=2,imax
   
      yflux_f(i,j)=                                                &
! Oy Advective flux at p point:
                -0.125*(fluxbar_v(i-1,j  ,0) +fluxbar_v(i,j,0))    &
!                      *(velbar_u(i  ,j-1,now)+velbar_u(i,j,now))  &
               *(  (velbar_u(i  ,j-1,now)+velbar_u(i,j  ,now)         &
                   +velbar_u(i  ,j-1,0  )+velbar_u(i,j  ,0  ))*1.125  &
                  -(velbar_u(i  ,j-2,now)+velbar_u(i,j+1,now)         &
                   +velbar_u(i  ,j-2,0  )+velbar_u(i,j+1,0  ))*0.125) &
! Oy Diffusive flux:
!......................................................................
!       +diflxfactor_*mask_f(i,j,kmax) & !voir notes !19-07-19!18-12-20
!                  *abs(dx_f(i,j)                                      &
!                      *(h_f(i,j)+0.25*(ssh_w(i  ,j  ,1)               &! schema cst   (et sa difinition de +diflxfactor_)
!                                      +ssh_w(i-1,j  ,1)               &
!                                      +ssh_w(i  ,j-1,1)               &
!                                      +ssh_w(i-1,j-1,1))))*( & !m�v�m>
!            biharm_e1*dy_f(i,j)*inv_dte_lp                       & ! biviscosite constante 09-06-18
!           +biharm_e2*0.5*abs(velbar_u(i,j-1,0)-velbar_u(i,j,0)) & ! biviscosite gradient de vitesse
!                                                           ) & !m�v�m>
        +diflxfactor_*0.5*(abs(fluxbar_v(i-1,j,0))+abs(fluxbar_v(i,j,0))) & !17-01-22
!......................................................................
      *( & !ooo>
!        (velbar_u(i  ,j  ,before )-velbar_u(i  ,j-1,before ))*3. &
!       -(velbar_u(i  ,j+1,before )-velbar_u(i  ,j-2,before ))    &
!          *mask_u(i  ,j+1,kmax   )  *mask_u(i  ,j-2,kmax)        & !22-07-17
         (velbar_u(i  ,j  ,before )-velbar_u(i  ,j-1,before )      &
         +velbar_u(i  ,j  ,now    )-velbar_u(i  ,j-1,now    ))*1.5 &
        -(velbar_u(i  ,j+1,before )-velbar_u(i  ,j-2,before )      &
         +velbar_u(i  ,j+1,now    )-velbar_u(i  ,j-2,now    ))*0.5 &
           *mask_u(i  ,j+1,kmax   )  *mask_u(i  ,j-2,kmax)        & !22-07-17
       )   !ooo>
      enddo
      enddo
! Notes: le flux advectif de yflux_f(i,j) est nul si mask_t(i-1,j)=mask_t(i,j)=0 (car alors fluxbar_v(i-1:i,j)=0)
! Le flux diffusif est nul (si *mask_f(i,j,kmax)) car dans le masque la hauteur de la colonne peut etre
! n'importe quoi et compromettre la CFL
! Par consequent cela correspond A une condition "free-slip" A la cOte. 


!     if(iteration2d==4.and.mod(iteration3d,10)==0) then
!      i=imax/2 ; j=jmax/2
!      if(mask_t(i,j,kmax)==1) &
!      write(10+par%rank,*) &
!     ,abs(dy_t(i,j)*(h_w(i,j)+ssh_w(i,j,1)))*( & !m�v�m>
!          biharm_e1*dx_t(i,j)*inv_dte_lp                        & ! biviscosite constante 09-06-18
!         +biharm_e2*0.5*abs(velbar_u(i,j,0)-velbar_u(i+1,j,0))  & ! biviscosite gradient de vitesse
!                                                          ) & !m�v�m>
!     ,diflxfactor_*0.5*abs(fluxbar_u(i,j,0) +fluxbar_u(i+1,j,0)) &
!     ,' trouvemi'
!     endif

!............
! v equation:
      do j=2,jmax
      do i=2,imax

! boucle stricte: do i=2,imax | do j=2,jmax
      xflux_f(i,j)=                                                &
! Ox Advective flux at p point:
                -0.125*(fluxbar_u(i  ,j-1,0) +fluxbar_u(i,j,0))    &
!                      *(velbar_v(i-1,j  ,now)+velbar_v(i,j,now))  &
                 *( (velbar_v(i-1,j  ,now)+velbar_v(i  ,j,now)         &
                    +velbar_v(i-1,j  ,0  )+velbar_v(i  ,j,0  ))*1.125  &
                   -(velbar_v(i-2,j  ,now)+velbar_v(i+1,j,now)         &
                    +velbar_v(i-2,j  ,0  )+velbar_v(i+1,j,0  ))*0.125) &
! Ox Diffusive flux at p point:
!..........................................................................
!       +diflxfactor_*mask_f(i,j,kmax) & !voir notes 19-07-19 !18-12-20
!                  *abs(dy_f(i,j)                                      &
!                      *(h_f(i,j)+0.25*(ssh_w(i  ,j  ,1)               &! schema cst   (et sa difinition de +diflxfactor_)
!                                      +ssh_w(i-1,j  ,1)               &
!                                      +ssh_w(i  ,j-1,1)               &
!                                      +ssh_w(i-1,j-1,1))))*( & !m�v�m>
!            biharm_e1*dx_f(i,j)*inv_dte_lp                       & ! biviscosite constante 09-06-18
!           +biharm_e2*0.5*abs(velbar_v(i,j,0)-velbar_v(i-1,j,0)) & ! biviscosite gradient de vitesse
!                                                           ) & !m�v�m>
        +diflxfactor_*0.5*(abs(fluxbar_u(i,j-1,0))+abs(fluxbar_u(i,j,0))) & !17-01-22
!..........................................................................
      *( & !ooo>
!        (velbar_v(i  ,j  ,before )-velbar_v(i-1,j  ,before ))*3. &
!       -(velbar_v(i+1,j  ,before )-velbar_v(i-2,j  ,before ))    &
!          *mask_v(i+1,j  ,kmax)     *mask_v(i-2,j  ,kmax)        & !22-07-17
         (velbar_v(i  ,j  ,before )-velbar_v(i-1,j  ,before )      &
         +velbar_v(i  ,j  ,now    )-velbar_v(i-1,j  ,now    ))*1.5 &
        -(velbar_v(i+1,j  ,before )-velbar_v(i-2,j  ,before )      &
         +velbar_v(i+1,j  ,now    )-velbar_v(i-2,j  ,now    ))*0.5 &
           *mask_v(i+1,j  ,kmax)     *mask_v(i-2,j  ,kmax)        & !22-07-17
       )   !ooo>
      enddo
      enddo
! Note le flux advectif de xflux_f(i,j) est nul si mask_t(i,j-1)=mask_t(i,j)=0 (car alors fluxbar_u(i,j-1:j)=0)
! Le flux diffusif est nul (si *mask_f(i,j,kmax)) car dans le masque la hauteur de la colonne peut etre
! n'importe quoi et compromettre la CFL
! Par consequent cela correspond A une condition "free-slip" A la cOte. 

!............
! v equation:
      do j=1,jmax
      do i=2,imax-1
      yflux_t(i,j)=                                                &
! Oy Advective flux at t point:
                -0.125*(fluxbar_v(i,j,0) +fluxbar_v(i,j+1,0))      &
!                      *(velbar_v(i,j,now)+velbar_v(i,j+1,now))    &
                 *( (velbar_v(i,j  ,now)+velbar_v(i,j+1,now)         &
                    +velbar_v(i,j  ,0  )+velbar_v(i,j+1,0  ))*1.125  &
                   -(velbar_v(i,j-1,now)+velbar_v(i,j+2,now)         &
                    +velbar_v(i,j-1,0  )+velbar_v(i,j+2,0  ))*0.125) &
! Oy Diffusive flux at t point:
!.......................................................................
!     +diflxfactor_*abs(dx_t(i,j)*(h_w(i,j)+ssh_w(i,j,1)))*( & !m�v�m> 
!             biharm_e1*dy_t(i,j)*inv_dte_lp                       & ! biviscosite constante 09-06-18
!            +biharm_e2*0.5*abs(velbar_v(i,j,0)-velbar_v(i,j+1,0)) & ! biviscosite gradient de vitesse 
!                                                          ) & !m�v�m>
      +diflxfactor_*0.5*(abs(fluxbar_v(i,j,0))+abs(fluxbar_v(i,j+1,0))) & !17-01-22
!.......................................................................
      *( & !ooo>
!        (velbar_v(i  ,j+1,before )-velbar_v(i  ,j  ,before ))*3. &
!       -(velbar_v(i  ,j+2,before )-velbar_v(i  ,j-1,before ))    &
!          *mask_v(i  ,j+2,kmax)     *mask_v(i  ,j-1,kmax)        & !22-07-17
         (velbar_v(i  ,j+1,before )-velbar_v(i  ,j  ,before )      &
         +velbar_v(i  ,j+1,now    )-velbar_v(i  ,j  ,now    ))*1.5 &
        -(velbar_v(i  ,j+2,before )-velbar_v(i  ,j-1,before )      &
         +velbar_v(i  ,j+2,now    )-velbar_v(i  ,j-1,now    ))*0.5 &
           *mask_v(i  ,j+2,kmax)     *mask_v(i  ,j-1,kmax)        & !22-07-17
       )   !ooo>

      enddo
      enddo

      return
      endif                   !-c4-biharm-version2->

      if(flag_adve2d==-1) then !-zero->
!//////////////////////////////////////////////////////////////////////////////////////////
! RIEN
       xflux_f=0.  ; xflux_t=0.  ; yflux_f=0.  ; yflux_t=0.  ; return
      endif                    !-zero->

      end subroutine external_adve2d

!.....................................

      subroutine external_add_adve_ref                 ! au terme adve3d2d (avant procedure filtrage frozenterms)
      implicit none

! Ajoute la reference d'advection barotrope (celle en iteration2d=1)
      do j=2,jmax-1 ; do i=2,imax
       adve3d2d_u(i,j)= &
       adve3d2d_u(i,j)  &
       +invdxdy_u(i,j)*(xflux_t(i,j)-xflux_t(i-1,j)+yflux_f(i,j+1)-yflux_f(i,j))
!      if(i==imax/2.and.j==jmax/2)write(20+par%rank,*)adve3d2d_u(i,j),invdxdy_u(i,j),real(1./dxdy_u(i,j))
       if(i==imax/2.and.j==jmax/2)write(20+par%rank,*)invdxdy_v(i,j),real(1./dxdy_v(i,j)),invdxdy_t(i,j),real(1./dxdy_t(i,j))
      enddo         ; enddo

      do j=2,jmax ; do i=2,imax-1
       adve3d2d_v(i,j)= &
       adve3d2d_v(i,j)  &
       +invdxdy_v(i,j)*(yflux_t(i,j)-yflux_t(i,j-1)+xflux_f(i+1,j)-xflux_f(i,j))
      enddo       ; enddo

      end subroutine external_add_adve_ref                 ! au terme adve3d2d (avant procedure filtrage frozenterms)

!.....................................

      end module module_external_mode
