










      module module_drifter
      use module_principal
      use module_parallele
      use pnetcdf
      implicit none
      integer :: drifter_dim2=9       &  !15-05-15
                ,rungekutta_loop=1    &
                ,rungekutta_order=2   &
                ,tag_gt1gt2_canal=10  &
                ,tag_gt2gt1_canal=11  &
                ,flag_stay_in_water=1 &  !21-09-20
                ,nbobuffercanalmax=1  &
                ,new_kbumax           &
                ,nexchgmax_webc=0     &
                ,id_csf               & ! identifiant pour: Corey Shape Factor
                ,id_dequi             & ! identifiant pour: Equivalent particle diameter
                ,id_density           & ! identifiant pour: particule density
                ,loop_buo_w           &
                ,loop_buo_max=1
!     double precision :: deci1,decj1,deck1
      double precision :: drifter_l1_tmp,drifter_l2_tmp,drifter_l3_tmp &
                         ,deci_drf_min,decj_drf_min,deci_drf_max,decj_drf_max &
                         ,udrifter(4)=0.,vdrifter(4)=0.,wdrifter(4)=0.  &
                         ,wrandom(2)=0.  
      real*4  :: random_variability=0.,deck_drf_min,deck_drf_max
!     integer,dimension(mpi_status_size) :: status_gt1gt2_canal  &  
!                                          ,status_gt2gt1_canal     
      integer,dimension(:),allocatable :: tabreq, kbumax_glb1,kbumax_glb2
      integer,dimension(:,:),allocatable :: tstatus
!______________________________________________________________________
! SYMPHONIE ocean model
! release 368 - last update: 30-03-23
!______________________________________________________________________

!______________________________________________________________________
! version date      Description des modifications
!         22/03/02: bienvenue à ICHOIX
!         27/03/02: debugage sigma generalisee
!         02/04/02: debug (mineur) date
!                   ecriture d'un fichier "fond de carte" pour visu gnu
!         07/10/02: ajout d'une date de "lacher" de drifter
!                   modif du format du fichier notebook_drifter
!         08/10/02: debug test sur kount
!         09/10/02: debug cas KBOMAX=0
!         14/10/02: debug cas notebook_drifter inexistant
!         04/05/06: amenagement pour compiler en double precision
!         11/03/09: un fichier par drifter. Ecriture avec format
!         27/04/09: notebook_drifter passe dans NOMFICHIER(16)
! 2010.5  29-01-10: nouvelles coordonnées & nouvelle terminologie
! 2010.6  02-02-10: renomme lon_t lat_t
! 2010.13 21-10-10  latlon_to_ij(1) devient latlon_to_ij('glb')
!         03-11-10  des arguments passés dans date_to_kount
! 2010.20 15-04-11  Calculs sur la base d'un temps exprimés en secondes
! 2010.25 29-01-12  La routine drifter remplace la routine bouees.F90
! 2010.25 02-02-12  Ajout echange mpi des coin
!         25-02-12  allocation dynamique
!         21-03-12  ajout de clefs de compilation "parallele"
! S.26    02-01-14  seul par%rank==0 ecrit a l'ecran
!         03-03-15  - ajout drifter_dim2, 2eme dimensiou du tableau drifter_l
!                   - drifter_l(:,5) contient la date de depart des drifter
!         15-05-15  x1,x2,x3 deviennent posx_,posy_,posz_
!         30-11-15  modifs drifter_initial
! v267    18-11-19  ajout derive de stokes
! v268    23-11-19  - calcul du deplacement vertical revu pour eviter les singularites de la couche fusionnee
!                   - division par dx, dy remplacees par multiplications par invdx, invdy
! v287    28-08-20  modif Claude
! v289    21-09-20  - elapsedtime_bef remplacE par elapsedtime_now-dt_drf
!                   - Permettre ou non (selon flag_stay_in_water) de deposer les particules sur le fond
! v290    30-10-20  correction d'un bug sur le mouvement aleatoire
! v292    15-11-20  ajout d'un message pour completer le stop 4
!         16-11-20  les tableaux de la buffer zone se reallouent en ligne si memoire 
!                   insuffisante en cours de simulation
! v345    06-04-22  plus de precision sur le format d'ecriture des fichiers individuels
! v353    04-09-22  1- Les echanges sont fait en une etape et non plus dans une boucle sur 
!                      le 2eme indice de drifter_l(:,:), voir details du gain de calcul dans:
!                      https://docs.google.com/document/d/13gPeQ15viDKH0fkoXPlK27Z-8huNH_z_hLNgtrTu3D4/edit?usp=sharing
!                   2- Les drifters peuvent transiter par les webcanaux
! v354    12-09-22  if(maxval(drifter(:,4))>9999999.)flag_stop=1
!         16-09-22  Ajout RK4 (et quelques modifs sur RK1 et RK2)
! v355    24-09-22  - schema random turbulent w en 2 temps "implicite"
!                   - schema random turbulent: possibilite qu'il soit indexE sur KH
!         01-10-22  Vitesse verticale en m/s (et non plus en indice/s) pour plus de precision 
!                   en maillage vertical inhommogene (typiquement la coord VQS)
! v357    24-10-22  Ecrire en netcdf
! v359    06-12-22  Ecrire la vitesse de chute dans le fichier netcdf
! v360    11-12-22  plus de messages d'erreurs
! v365    26-01-23  utilisation rho_t
! v367    28-02-23  - ajout de i,j,k dans fichier sortie netcdf
!                   - vitesse de chute en ligne DOI: 10.1021/acs.est.8b06794
!         07-03-23  - vitesse de chute en ligne suite
! v368    30-03-23  dans le cas RK2 passer en predictor/corrector quand la vitesse du temps intermediaire est nulle
!         20-04-23  traiter le cas k1<1 dans la renumerotation
!______________________________________________________________________

contains

!...............................................................................

      subroutine drifter_update
!     use module_principal
!     use module_parallele
      implicit none
      double precision zdrifter_now,zdrifter_after
      real*4 reynolds_,cd_,waterdensity_


      if(modulo(iteration3d,dt_drf_over_dti_fw)/=0) return !24-11-19

! Si la particule doit rester dans la colonne d'eau on applique un seuil (ici deck_drf_min)
! coherent sur la position verticale, sinon ce seuil est zero pour autoriser drifter_l(kbu,3)
! a devenir plus petit que 1 pour ensuite ne plus deplacer la particule.
!     if(flag_stay_in_water==1) then !>>>>
!      deck_drf_min=0.51
!     else                           !>>>>
!      deck_drf_min=0. !21-09-20
!     endif                          !>>>>

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! MISE A JOUR DE LA POSITION DES BOUEES                                !22/03/02
! DEBUT:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!.................................................
! Definition d'une echelle de vitesse verticale turbulente archive dans anyvar3d
      if(drifter_random_w/=0)call drifter_turbulence_w !24-09-22

      do 29 kbu=1,kbumax !21-09-20

! Si, parce que flag_stay_in_water=0, par effet de la sedimentation, deck=drifter_l(kbu,3)<1, alors
! on considere que la particule s'est deposee sur le fond et reste immobile. Dans cette eventualitE
! le calcul est conditionnE au fait que drifter_l(kbu,3)>=1.
      if(drifter_l(kbu,3)>0.5) then !-flag_stay_in_water-> !21-09-20 !01-10-22

! Note: on calule le nombre aleatoire dans boucle 29 (un nombre par drifter) mais pas dans la boucle du runge-kuta
      if(drifter_random_w==1) then !-tke-scheme->
! random_variability nombre aleatoire muktipliant anyvar3d. 
         call random_number(random_variability)
! La normalisation en suivant permet que l'energie cinetique deduite de
! la vitesse verticale aleatoire est bien tken_w/3 en moyenne. Voir la
! routine fortran test_random en fin du present fichier.
         random_variability=2.*sqrt(3.)*(random_variability-0.5)
      endif                        !-tke-scheme->
      if(drifter_random_w==2) then !-kh-scheme->
         call random_number(random_variability)
         random_variability=4.8*(random_variability-0.5)
      endif                        !-kh-scheme->

      do 19 rungekutta_loop=1,rungekutta_order ! rungekutta_order=2

      if(elapsedtime_now>=drifter_l(kbu,5)) then !%%%%%%%%%%%%%%%%%%%%%%> !03-03-15

      if(rungekutta_loop==1) then !pmx>

       deci=drifter_l(kbu,1)-par%timax(1)
       decj=drifter_l(kbu,2)-par%tjmax(1)
       deck=drifter_l(kbu,3)

      else                        !pmx>

       deci=drifter_l1_tmp-par%timax(1)
       decj=drifter_l2_tmp-par%tjmax(1)
       deck=drifter_l3_tmp
! Attention que la position temporaire de la particule est en dehors du sousdomaine mpi (auquel cas c'est la position du depart qui est reprise)
! On note toutefois que les bornes 0.501, imax+0.499 sont beaucoup plus tolerantes que deci_drf_min etc... qui ont pour
! consequence qu'il n'y a au depart du calcul aucun drifter en deca de i=1.5 ou au dela de i=imax-0.5, ce qui signifie qu'il y a
! tres peu de chance (si le pas de temps dt_drf est raisonable) de se trouver contraint par le critere suivant.
!         if(deci<=1.or.deci>=imax.or.decj<=1.or.decj>=jmax) then
          if(deci<=0.501.or.deci>=real(imax)+0.499   &
         .or.decj<=0.501.or.decj>=real(jmax)+0.499) then !-securite-mpi->
            deci=drifter_l(kbu,1)-par%timax(1)
            decj=drifter_l(kbu,2)-par%tjmax(1)
            deck=drifter_l(kbu,3)
          endif                                          !-securite-mpi->

      endif                       !pmx>

      i=int(deci)
      j=int(decj)
      k=int(deck)
      rapi=deci-real(i)
      rapj=decj-real(j)
      rapk=deck-real(k)

!     if(i<1.or.i>imax.or.j<1.or.j>jmax) &
!     if(deci<1.or.deci>real(imax).or.decj<1.or.decj>real(jmax)) then !>>>
       if(deci<=0.501.or.deci>=real(imax)+0.499   &
      .or.decj<=0.501.or.decj>=real(jmax)+0.499) then !-securite-mpi->
         flag_stop=1
         write(10+par%rank,*)'................'
         write(10+par%rank,*)'kbu',kbu   
         write(10+par%rank,*)'drifter_l(1:4)',real(drifter_l(kbu,1:4)) !11-12-22
         write(10+par%rank,*)'deci,decj glb',deci+par%timax(1),decj+par%tjmax(1)
         write(10+par%rank,*)'deci,decj loc',deci,decj
         write(10+par%rank,*)'imax,jmax',imax,jmax
         write(10+par%rank,*)'iglb,jglb',iglb,jglb
         write(10+par%rank,*)'rungekutta_loop',rungekutta_loop
         write(10+par%rank,*)'iteration3d',iteration3d
         write(10+par%rank,*)'par%rank',par%rank
!     stop 'indices en dehors des limites dans drifter_update'
      endif                                          !-securite-mpi->

!            /  /  /


! Recapitulation:
! DECI DECJ DECK sont les 3 indices (decimaux) caracterisant
! la position de la drifter sur la grille. Attention au fait
! que la grille C comporte en fait plusieurs grille: points
! de courant et points de masse qui se distinguent sur l'horizontal
! et points de courant horizontal et vertical qui se distinguent
! sur la vertical. Le triplet d'indice repere la drifter dans
! la grille des points de temperature c'est à dire points _Z et
! niveaux verticaux intermediares. Ceci explique les "glissements"
! d'indice pour aller chercher une composante u du courant (DECI+0.5)
! une composante v (DECJ+0.5) une composante omega (DECK+0.5).

      if(flag_buo_w==4) then !-VITESSE PROPRE CAS 4-> !07-03-23
! Cas oU la vitesse de chute se deduit de:
! Effects of Particle Properties on the Settling and Rise Velocities of
! Microplastics in Freshwater under Laboratory Conditions
!      DOI: 10.1021/acs.est.8b06794
!      Environ. Sci. Technol. 2019, 53, 1958¿1966

! Pour le moment on ne remet A jour la vitesse de chute qu'A la premiere boucle
! mais cela peut etre modifiE A l'avenir...
      if(rungekutta_loop==1) then !pmx>

      k1=max0(k  ,1)
      k2=min0(k+1,kmax)

! densite de l'eau au niveau de la particule:
      waterdensity_=                                       & 
       ( (1.-rapi)*(1.-rapj)*(1.-rapk)*mask_t(i  ,j  ,kmax)*rho_t(i  ,j  ,k1)   &
        +(1.-rapi)*(1.-rapj)*    rapk *mask_t(i  ,j  ,kmax)*rho_t(i  ,j  ,k2)   &
        +(1.-rapi)*    rapj *(1.-rapk)*mask_t(i  ,j+1,kmax)*rho_t(i  ,j+1,k1)   &
        +(1.-rapi)*    rapj *    rapk *mask_t(i  ,j+1,kmax)*rho_t(i  ,j+1,k2)   &
        +    rapi *(1.-rapj)*(1.-rapk)*mask_t(i+1,j  ,kmax)*rho_t(i+1,j  ,k1)   &
        +    rapi *(1.-rapj)*    rapk *mask_t(i+1,j  ,kmax)*rho_t(i+1,j  ,k2)   &
        +    rapi *    rapj *(1.-rapk)*mask_t(i+1,j+1,kmax)*rho_t(i+1,j+1,k1)   &
        +    rapi *    rapj *    rapk *mask_t(i+1,j+1,kmax)*rho_t(i+1,j+1,k2) ) &
      /max( &
         (1.-rapi)*(1.-rapj)*(1.-rapk)*mask_t(i  ,j  ,kmax)   &
        +(1.-rapi)*(1.-rapj)*    rapk *mask_t(i  ,j  ,kmax)   &
        +(1.-rapi)*    rapj *(1.-rapk)*mask_t(i  ,j+1,kmax)   &
        +(1.-rapi)*    rapj *    rapk *mask_t(i  ,j+1,kmax)   &
        +    rapi *(1.-rapj)*(1.-rapk)*mask_t(i+1,j  ,kmax)   &
        +    rapi *(1.-rapj)*    rapk *mask_t(i+1,j  ,kmax)   &
        +    rapi *    rapj *(1.-rapk)*mask_t(i+1,j+1,kmax)   &
        +    rapi *    rapj *    rapk *mask_t(i+1,j+1,kmax),1.e-10)  

       if(waterdensity_==0.) then
         write(10+par%rank,*)'Err: waterdensity_=0'
         write(10+par%rank,*)'i,j,k',i,j,k
         write(10+par%rank,*)'rho_t(i  ,j  ,k1)',rho_t(i  ,j  ,k1)
         write(10+par%rank,*)'rho_t(i+1,j+1,k2)',rho_t(i+1,j+1,k2)
         flag_stop=1 
       endif

       do loop_buo_w=1,loop_buo_max ! Note: boucle iterative mais loop_buo_max peut etre egal A 1

       reynolds_= & ! EQ9 : abs(ws)*dequi/kinvis
                abs(drifter_l(kbu,id_wdrifter)) & ! vitesse de flottabilite de la particule
                   *drifter_l(kbu,id_dequi)     & ! Equivalent  particle diameter    
                   /1.e-6                         ! Kinematic viscosity (m**2/s)

       cd_=4.7/sqrt(max(reynolds_,1.e-10))+sqrt(drifter_l(kbu,id_csf)) ! EQ11: CD=4.7/sqrt(reynolds)+sqrt(csf)

! EQ14: ws=sqrt((4./3.)*(dequi/cd)*(abs(rhos-rhow)/rhow)*grav) 
!   *sign(1.,rhow-rhos) ! signe de Ws = signe du delta de densite

      drifter_l(kbu,id_wdrifter) & ! EQ14: vitesse de flottabilite de la particule
      =sqrt( &!-SQRT->
            4./3.*(drifter_l(kbu,id_dequi)/cd_) &
                 *(abs(drifter_l(kbu,id_density)-waterdensity_)      &
                                                /waterdensity_)*grav &
           ) &!-SQRT->
           *sign(1.,waterdensity_-drifter_l(kbu,id_density)) ! signe de Ws = signe du delta de densite

!       write(10+par%rank,*)iteration3d,kbu,' iteration3d,kbu'
!       write(10+par%rank,*)'i,j,k',i,j,k
!       write(10+par%rank,*)'rapi,rapj,rapk',rapi,rapj,rapk
!       write(10+par%rank,*)'k1,k2',k1,k2
!       write(10+par%rank,*)'rho_t(i  ,j  ,k1)',rho_t(i  ,j  ,k1)
!       write(10+par%rank,*)reynolds_,cd_,waterdensity_,' Re,Cd,Rhow'
!       write(10+par%rank,*)drifter_l(kbu,id_csf),' CSF'
!       write(10+par%rank,*)drifter_l(kbu,id_dequi),' Dequi'
!       write(10+par%rank,*)drifter_l(kbu,id_density),' rhoParticule'
!       write(10+par%rank,*)drifter_l(kbu,id_wdrifter),' WS'

       enddo ! loop_buo_w

      endif                       !pmx>

      endif                  !-VITESSE PROPRE CAS 4-> !07-03-23


!            /  /  /

!.......................................................................!
! interpolation de la vitesse U au point (DECI,DECJ,DECK)
! (attention au glissement d'indice sur i)
! Debut:
!.......................................................................!

      i=int(deci+0.5)
      rapi=deci+0.5-real(i)


! à ce stade K est compris entre 0 et NR-1 donc attention à ne pas
! chercher des niveaux de courant qui n'existent pas:
      k1=max0(k  ,1)
      k2=min0(k+1,kmax)

! x1*dx = premiere composante de la vitesse sur le point decimal (DECI,DECJ,DECK)
      udrifter(rungekutta_loop)=                              & !16-09-22 
         (1.-rapi)*(1.-rapj)*(1.-rapk)*(vel_u(i  ,j  ,k1,1)   &
                                 +velstokes_u(i  ,j  ,k1,1)   &
                                    )*invdx_u(i  ,j  )        & !23-11-19
        +(1.-rapi)*(1.-rapj)*    rapk *(vel_u(i  ,j  ,k2,1)   &
                                 +velstokes_u(i  ,j  ,k2,1)   &
                                    )*invdx_u(i  ,j  )        &
        +(1.-rapi)*    rapj *(1.-rapk)*(vel_u(i  ,j+1,k1,1)   &
                                 +velstokes_u(i  ,j+1,k1,1)   &
                                    )*invdx_u(i  ,j+1)        &
        +(1.-rapi)*    rapj *    rapk *(vel_u(i  ,j+1,k2,1)   &
                                 +velstokes_u(i  ,j+1,k2,1)   &
                                    )*invdx_u(i  ,j+1)        &
        +    rapi *(1.-rapj)*(1.-rapk)*(vel_u(i+1,j  ,k1,1)   &
                                 +velstokes_u(i+1,j  ,k1,1)   &
                                    )*invdx_u(i+1,j  )        &
        +    rapi *(1.-rapj)*    rapk *(vel_u(i+1,j  ,k2,1)   &
                                 +velstokes_u(i+1,j  ,k2,1)   &
                                    )*invdx_u(i+1,j  )        &
        +    rapi *    rapj *(1.-rapk)*(vel_u(i+1,j+1,k1,1)   &
                                 +velstokes_u(i+1,j+1,k1,1)   &
                                    )*invdx_u(i+1,j+1)        &
        +    rapi *    rapj *    rapk *(vel_u(i+1,j+1,k2,1)   &
                                 +velstokes_u(i+1,j+1,k2,1)   &
                                    )*invdx_u(i+1,j+1)

!.......................................................................!
! interpolation de la vitesse U au point (DECI,DECJ,DECK)
! (attention au glissement d'indice sur i)
! Fin.
!.......................................................................!


!            /  /  /

!.......................................................................!
! interpolation de la vitesse V au point (DECI,DECJ,DECK)
! (attention au glissement d'indice sur j)
! Debut:
!.......................................................................!

      i=int(deci)
      j=int(decj+0.5)
      rapi=deci    -real(i)
      rapj=decj+0.5-real(j)

      vdrifter(rungekutta_loop)=                            & !16-09-22 
         (1.-rapi)*(1.-rapj)*(1.-rapk)*(vel_v(i  ,j  ,k1,1) &
                                 +velstokes_v(i  ,j  ,k1,1) &
                                    )*invdy_v(i  ,j  )      & !23-11-19
        +(1.-rapi)*(1.-rapj)*    rapk *(vel_v(i  ,j  ,k2,1) &
                                 +velstokes_v(i  ,j  ,k2,1) &
                                    )*invdy_v(i  ,j  )      &
        +(1.-rapi)*    rapj *(1.-rapk)*(vel_v(i  ,j+1,k1,1) &
                                 +velstokes_v(i  ,j+1,k1,1) &
                                    )*invdy_v(i  ,j+1)      &
        +(1.-rapi)*    rapj *    rapk *(vel_v(i  ,j+1,k2,1) &
                                 +velstokes_v(i  ,j+1,k2,1) &
                                    )*invdy_v(i  ,j+1)      &
        +    rapi *(1.-rapj)*(1.-rapk)*(vel_v(i+1,j  ,k1,1) &
                                 +velstokes_v(i+1,j  ,k1,1) &
                                    )*invdy_v(i+1,j  )      &
        +    rapi *(1.-rapj)*    rapk *(vel_v(i+1,j  ,k2,1) &
                                 +velstokes_v(i+1,j  ,k2,1) &
                                    )*invdy_v(i+1,j  )      &
        +    rapi *    rapj *(1.-rapk)*(vel_v(i+1,j+1,k1,1) &
                                 +velstokes_v(i+1,j+1,k1,1) &
                                    )*invdy_v(i+1,j+1)      &
        +    rapi *    rapj *    rapk *(vel_v(i+1,j+1,k2,1) &
                                 +velstokes_v(i+1,j+1,k2,1) &
                                    )*invdy_v(i+1,j+1)

!.......................................................................!
! interpolation de la vitesse V au point (DECI,DECJ,DECK)
! (attention au glissement d'indice sur j)
! Fin.
!.......................................................................!


!            /  /  /

!.......................................................................!
! interpolation de la vitesse omega au point (DECI,DECJ,DECK)
! (attention au glissement d'indice sur k)
! Debut:
!.......................................................................!

      j=int(decj)
      k=int(deck+0.5)
      rapj=decj    -real(j)
      rapk=deck+0.5-real(k)

      k=min0(max0(k,1),kmax)

! x3*dz = troisieme composante de la vitesse sur le point decimal (DECI,DECJ,DECK)
! A partir de la version 268 la division par dz n'est pas appliquee individuellement pour eviter ponctuellement une division par zero dans !23-11-19
! une couche fusionnee, mais globalement. La division globale evite les singularites et ce d'autant plus que deck est bornE de sorte qu'une 
! particule n'est en principe jamais A 100% dans la couche fusionnee.
! On note au passage que le calcul ne contient plus qu'une seule division et pas 8 comme dans l'ancien schema

! Position z de la particule au temps now !01-10-22
      if(rungekutta_loop==1) then !m°v°m>
        zdrifter_now=                                             & 
         (1.-rapi)*(1.-rapj)*(1.-rapk)*depth_w(i  ,j  ,k  )       & 
        +(1.-rapi)*(1.-rapj)*    rapk *depth_w(i  ,j  ,k+1)       &
        +(1.-rapi)*    rapj *(1.-rapk)*depth_w(i  ,j+1,k  )       &
        +(1.-rapi)*    rapj *    rapk *depth_w(i  ,j+1,k+1)       &
        +    rapi *(1.-rapj)*(1.-rapk)*depth_w(i+1,j  ,k  )       &
        +    rapi *(1.-rapj)*    rapk *depth_w(i+1,j  ,k+1)       &
        +    rapi *    rapj *(1.-rapk)*depth_w(i+1,j+1,k  )       &
        +    rapi *    rapj *    rapk *depth_w(i+1,j+1,k+1)        
      endif                       !m°v°m>

! Note: attention wdrifter est en m/s et pas en indice/s comme u et v parce que le maillage vertical 
! est trop irregulier pour que l'approche en indice/s soit assez precise, surtout avec la coordonnee VQS.
! Linconvenient est de devoir convertir z en deck (plus de calculs donc). 
      wdrifter(rungekutta_loop)=                                    & !01-10-22
         (1.-rapi)*(1.-rapj)*(1.-rapk)*omega_w(i  ,j  ,k  ,1)       & !courant

        +(1.-rapi)*(1.-rapj)*    rapk *omega_w(i  ,j  ,k+1,1)       &

        +(1.-rapi)*    rapj *(1.-rapk)*omega_w(i  ,j+1,k  ,1)       &

        +(1.-rapi)*    rapj *    rapk *omega_w(i  ,j+1,k+1,1)       &

        +    rapi *(1.-rapj)*(1.-rapk)*omega_w(i+1,j  ,k  ,1)       &

        +    rapi *(1.-rapj)*    rapk *omega_w(i+1,j  ,k+1,1)       &

        +    rapi *    rapj *(1.-rapk)*omega_w(i+1,j+1,k  ,1)       &

        +    rapi *    rapj *    rapk *omega_w(i+1,j+1,k+1,1)          

! ligne suivante commentee le !30-03-23
!       +min(flag_buo_w,1)*drifter_l(kbu,id_wdrifter) !vitesse de flottabilite de la particule !28-02-23

      if(drifter_random_w/=0.and.rungekutta_loop<=2)          & 
         wrandom(rungekutta_loop)=                            & !01-10-22 
         random_variability*( & !ooo>
         (1.-rapi)*(1.-rapj)*(1.-rapk)*anyvar3d(i  ,j  ,k  )  & !vitesse turbulente aleatoire
        +(1.-rapi)*(1.-rapj)*    rapk *anyvar3d(i  ,j  ,k+1)  &
        +(1.-rapi)*    rapj *(1.-rapk)*anyvar3d(i  ,j+1,k  )  &
        +(1.-rapi)*    rapj *    rapk *anyvar3d(i  ,j+1,k+1)  &
        +    rapi *(1.-rapj)*(1.-rapk)*anyvar3d(i+1,j  ,k  )  &
        +    rapi *(1.-rapj)*    rapk *anyvar3d(i+1,j  ,k+1)  &
        +    rapi *    rapj *(1.-rapk)*anyvar3d(i+1,j+1,k  )  &
        +    rapi *    rapj *    rapk *anyvar3d(i+1,j+1,k+1)  &
                            )   !ooo>


! Note: le time stepping du mouvement aleatoire n'est pas un runge-kutta
! mais un schema "implicite" en 2 temps pour eviter d'etre ejecter des
! couches limites. C'est wrandom(2) qui fait un final le mouvement.
! Necessitant 2 boucles de temps il ne peut pas etre "implicite" en cas
! de RK1 (il est alors explicite en avant et peu performant)
!.............!16-09-22
! RK1
      if(rungekutta_order==1) then !-rk1->

! RK1 ou RK4: ajouter la vitesse de chute pour former la vitesse verticale totale maintenant: !30-03-23
        wdrifter(rungekutta_loop)= &
        wdrifter(rungekutta_loop)  &
        +min(flag_buo_w,1)*drifter_l(kbu,id_wdrifter) !vitesse de flottabilite de la particule !30-03-23

        zdrifter_after=zdrifter_now+dt_drf*wdrifter(1) &
                                   +dt_drf*wrandom(1)   
      endif                        !-rk1->
!.............!16-09-22
! RK2
      if(rungekutta_order==2) then !-rk2->
! Attention avec RK2 wdrifter ne contient que la vitesse omega , la vitesse de flottabilite etant par consequent ajoutee en suivant !30-03-23
       if(rungekutta_loop==1) then !-rk2-step1->
        zdrifter_after=zdrifter_now+dt_drf*wdrifter(1) &
                                   +dt_drf*wrandom(1)  &
                                   +dt_drf*min(flag_buo_w,1)*drifter_l(kbu,id_wdrifter) !vitesse de flottabilite de la particule !30-03-23
       endif                       !-rk2-step1->
       if(rungekutta_loop==2) then !-rk2-step2->
! Attention ici wdrifter(2) ne contient que la vitesse omega (et pas la vitesse de flottabilite)
! Si wdrifter(2)=0 alors appliquer predictor/corrector pour la partie wdrifter (ce qui revient A ne pas prendre en compte wdrifter dans le deplacement)
! Le but est d'eviter l'echouage autre que l'echouage provoquE par la vitesse de flottabilite, les particules "neutres" ne s'echouant jamais.
       if(wdrifter(2)/=0.) then !ooo> !30-03-23
        zdrifter_after=zdrifter_now+dt_drf*0.5*(wdrifter(1)+wdrifter(2)) &
                                   +dt_drf*wrandom(2) & ! schema "implicite" pour ne pas etre ejecter des couches limites
                                   +dt_drf*min(flag_buo_w,1)*drifter_l(kbu,id_wdrifter) !vitesse de flottabilite de la particule !30-03-23
       else                     !ooo> !30-03-23
        zdrifter_after=zdrifter_now &
                                   +dt_drf*wrandom(2) & ! schema "implicite" pour ne pas etre ejecter des couches limites
                                   +dt_drf*min(flag_buo_w,1)*drifter_l(kbu,id_wdrifter) !vitesse de flottabilite de la particule !30-03-23
       endif                    !ooo> !30-03-23
       endif                       !-rk2-step2->
      endif                        !-rk2->
!.............!16-09-22
! RK4
      if(rungekutta_order==4) then !-rk4->

! RK1 ou RK4: ajouter la vitesse de chute pour former la vitesse verticale totale maintenant: !30-03-23
        wdrifter(rungekutta_loop)= &
        wdrifter(rungekutta_loop)  &
        +min(flag_buo_w,1)*drifter_l(kbu,id_wdrifter) !vitesse de flottabilite de la particule !30-03-23

       if(rungekutta_loop==1) then !-rk4-step1->
        zdrifter_after=zdrifter_now+dt_drf*0.5*wdrifter(1) &
                                   +dt_drf*wrandom(1)
       endif                       !-rk4-step1->
       if(rungekutta_loop==2) then !-rk4-step2->
        zdrifter_after=zdrifter_now+dt_drf*0.5*wdrifter(2)  
       endif                       !-rk4-step2->
       if(rungekutta_loop==3) then !-rk4-step3->
        zdrifter_after=zdrifter_now+dt_drf*wdrifter(3)  
       endif                       !-rk4-step3->
       if(rungekutta_loop==4) then !-rk4-step4->
        zdrifter_after=zdrifter_now+dt_drf*(wdrifter(1)+2.*wdrifter(2)+2.*wdrifter(3)+wdrifter(4))/6. &
                                   +dt_drf*wrandom(2) ! schema "implicite" pour ne pas etre ejecter des couches limites
       endif                       !-rk4-step4->
      endif                        !-rk4->

       do while ( zdrifter_after                            &
                 > (1.-rapi)*(1.-rapj)*depth_w(i  ,j  ,k+1) &
                  +(1.-rapi)*    rapj *depth_w(i  ,j+1,k+1) &
                  +    rapi *(1.-rapj)*depth_w(i+1,j  ,k+1) &
                  +    rapi *    rapj *depth_w(i+1,j+1,k+1) &
           .and.k+1<kmax+1)
        k=k+1
       enddo
       do while ( zdrifter_after                            &
                 < (1.-rapi)*(1.-rapj)*depth_w(i  ,j  ,k  ) & 
                  +(1.-rapi)*    rapj *depth_w(i  ,j+1,k  ) &
                  +    rapi *(1.-rapj)*depth_w(i+1,j  ,k  ) &
                  +    rapi *    rapj *depth_w(i+1,j+1,k  ) &
           .and.k>1 )
        k=k-1
       enddo

! Pour la conversion de z en deck:
       rap=min(max(    (zdrifter_after                   &
           -( (1.-rapi)*(1.-rapj)*depth_w(i  ,j  ,k  )   & 
             +(1.-rapi)*    rapj *depth_w(i  ,j+1,k  )   &
             +    rapi *(1.-rapj)*depth_w(i+1,j  ,k  )   &
             +    rapi *    rapj *depth_w(i+1,j+1,k  ))) &
          /( &
            ( (1.-rapi)*(1.-rapj)*depth_w(i  ,j  ,k+1)   &
             +(1.-rapi)*    rapj *depth_w(i  ,j+1,k+1)   &
             +    rapi *(1.-rapj)*depth_w(i+1,j  ,k+1)   &
             +    rapi *    rapj *depth_w(i+1,j+1,k+1))  &
           -( (1.-rapi)*(1.-rapj)*depth_w(i  ,j  ,k  )   & 
             +(1.-rapi)*    rapj *depth_w(i  ,j+1,k  )   &
             +    rapi *(1.-rapj)*depth_w(i+1,j  ,k  )   &
             +    rapi *    rapj *depth_w(i+1,j+1,k  ))  &
           ),0.),1.)


!.......................................................................!
! interpolation de la vitesse omega au point (DECI,DECJ,DECK)
! (attention au glissement d'indice sur k)
! Fin.
!.......................................................................!


!.............!16-09-22
! RK1
      if(rungekutta_order==1) then !-rk1->
        drifter_l(kbu,1)=        drifter_l(kbu,1)+udrifter(1)*dt_drf
        drifter_l(kbu,2)=        drifter_l(kbu,2)+vdrifter(1)*dt_drf
        drifter_l(kbu,3)=min(max(rap*(k+1)+(1.-rap)*k-0.5,deck_drf_min),deck_drf_max)
      endif                        !-rk1->

!.............!16-09-22
! RK2
      if(rungekutta_order==2) then !-rk2->
       if(rungekutta_loop==1) then !-rk2-step1->
        drifter_l1_tmp=        drifter_l(kbu,1)+udrifter(1)*dt_drf
        drifter_l2_tmp=        drifter_l(kbu,2)+vdrifter(1)*dt_drf
        drifter_l3_tmp=min(max(rap*(k+1)+(1.-rap)*k-0.5,deck_drf_min),deck_drf_max)
       endif                       !-rk2-step1->
       if(rungekutta_loop==2) then !-rk2-step2->
! Le but du predictor/corrector est d'eviter l'echouage (autre que l'echouage provoquE par la vitesse de flottabilite, les particules
! "neutres" ne s'echouant jamais). Donc si le point d'arrivee au temps 1 est dans le masque (u=0,v=0) on bascule sur le predictor/corrector 
! a la place de RK2. La vitesse "predictor" etant nulle, cela equivaut A ne pas deplacer la particule !30-03-23
        if(udrifter(2)/=0.) & !predictor/corrector anti echouage si u(2)=0 !30-03-23
        drifter_l(kbu,1)=        drifter_l(kbu,1)+0.5*(udrifter(1)+udrifter(2))*dt_drf
        if(vdrifter(2)/=0.) & !predictor/corrector anti echouage si v(2)=0 !30-03-23
        drifter_l(kbu,2)=        drifter_l(kbu,2)+0.5*(vdrifter(1)+vdrifter(2))*dt_drf
        drifter_l(kbu,3)=min(max(rap*(k+1)+(1.-rap)*k-0.5,deck_drf_min),deck_drf_max)
       endif                       !-rk2-step2->
      endif                        !-rk2->

!.............!16-09-22
! RK4
      if(rungekutta_order==4) then !-rk4->
       if(rungekutta_loop==1) then !-rk4-step1->
        drifter_l1_tmp=        drifter_l(kbu,1)+0.5*udrifter(1)*dt_drf
        drifter_l2_tmp=        drifter_l(kbu,2)+0.5*vdrifter(1)*dt_drf
        drifter_l3_tmp=min(max(rap*(k+1)+(1.-rap)*k-0.5,deck_drf_min),deck_drf_max)
       endif                       !-rk4-step1->
       if(rungekutta_loop==2) then !-rk4-step2->
        drifter_l1_tmp=        drifter_l(kbu,1)+0.5*udrifter(2)*dt_drf
        drifter_l2_tmp=        drifter_l(kbu,2)+0.5*vdrifter(2)*dt_drf
        drifter_l3_tmp=min(max(rap*(k+1)+(1.-rap)*k-0.5,deck_drf_min),deck_drf_max)
       endif                       !-rk4-step2->
       if(rungekutta_loop==3) then !-rk4-step3->
        drifter_l1_tmp=        drifter_l(kbu,1)+udrifter(3)*dt_drf
        drifter_l2_tmp=        drifter_l(kbu,2)+vdrifter(3)*dt_drf
        drifter_l3_tmp=min(max(rap*(k+1)+(1.-rap)*k-0.5,deck_drf_min),deck_drf_max)
       endif                       !-rk4-step3->
       if(rungekutta_loop==4) then !-rk4-step4->
        drifter_l(kbu,1)=        drifter_l(kbu,1)+(udrifter(1)+2.*udrifter(2)+2.*udrifter(3)+udrifter(4))*dt_drf/6.
        drifter_l(kbu,2)=        drifter_l(kbu,2)+(vdrifter(1)+2.*vdrifter(2)+2.*vdrifter(3)+vdrifter(4))*dt_drf/6.
        drifter_l(kbu,3)=min(max(rap*(k+1)+(1.-rap)*k-0.5,deck_drf_min),deck_drf_max)
       endif                       !-rk4-step4->
      endif                        !-rk4->

      endif                                   !%%%%%%%%%%%%%%%%%%%%%%> !03-03-15
   19 continue

      endif                         !-flag_stay_in_water->  !21-09-20
   29 continue ! boucle sur kbu                             !21-09-20

      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) 
      if(k0/=0)stop 'Erreur drifters. Voir fort.xxx'

      call drifter_exchange_driver

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! MISE A JOUR DE LA POSITION DES BOUEES                                !22/03/02
! FIN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end subroutine drifter_update

!................................................................................

      subroutine drifter_initial
!     use module_principal
!     use module_parallele
      implicit none
      integer, dimension(:) , allocatable :: id_fieldrifter
      real   , dimension(:) , allocatable :: val_fieldrifter
      double precision posx_,posy_,posz_ !15-05-15
      real vertical_buoyancy_velocity
      integer kbuglb_,kglbmax_,in_or_out_,k_
      character*60 txt_spy_

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! INITIALISATION DES BOUEES                                            !22/03/02
! DEBUT:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


! Reset:
      txt_spy_   ='drifter_initial'
      kbuglb_   =0
      kbu=0
      kbumax=0               ! Nb de drifter dans le proc (inconnu à ce stade)
      in_or_out_   =0        ! 0=out 1=in
! bornes min max pour decider d'un echange mpi:
      deci_drf_min=1.5      ! 1. dans precedente version
      decj_drf_min=1.5      ! 1. dans precedente version
      deci_drf_max=imax-0.5 ! imax dans precedente version
      decj_drf_max=jmax-0.5 ! jmax dans precedente version 

      if(par%rank==0)write(6,*)'lecture de ',nomfichier(16) !02-01-14
      open(unit=3,file=nomfichier(16)) ! lecture de notebook_drifter    !27/04/09
      read(3,*)drifter_onoff
  400 continue

      if(drifter_onoff==0) then
        close(3) !26-01-23
        return
      endif
!     if(drifter_onoff==1) then !>>> !26-01-23
! note 1: si vitesse verticale propre, l'evaluation de cette derniere peut necessiter la connaissance de la densite vraie rho_t
! note 2: rho_t est aussi possiblement allouE dans set_parameters si analyse harmonique 3D activee
!       if(.not.allocated(rho_t)) then !pmx>
!         allocate(rho_t(0:imax+1,0:jmax+1,kmax)) !26-01-23
!       endif                          !pmx>
!      endif                    !>>>



      read(3,'(a)')drifter_initial_mode
      read(3,*)drifter_out_sampling
      if(drifter_out_sampling>0.) &
      drifter_out_sampling=max(drifter_out_sampling,dti_fw) !24-11-19
      read(3,*)drifter_output_files    !24-11-19

!........................ 
! allocation des tableaux (une fois connu drifter_output_files)
      nbomax=2400000 ; nbobuffermax=24000 ; nbobuffermax_c=2400  !fixer les dimensions
      call drifter_allocate(0) !16-11-20
!........................ 

      read(3,*)rungekutta_order        !19-11-19
      read(3,*)dt_drf_over_dti_fw      !24-11-19
        dt_drf=dt_drf_over_dti_fw*dti_fw

      read(3,*)flag_buo_w,id_wdrifter  !22-11-19
! Cas ou plusieurs emplacement son utilises pour calculer Ws !28-02-23
! note: meme si flag_buo_w>1 l'emplacement id_wdrifter servira a stocker Ws
      if(flag_buo_w>1) then !m°v°m> !28-02-23
       backspace 3
       allocate(id_fieldrifter(flag_buo_w-1))
       allocate(val_fieldrifter(flag_buo_w-1))
       read(3,*)flag_buo_w                   &
               ,id_wdrifter                  & ! indice de stockage de Ws
               ,id_fieldrifter(1:flag_buo_w-1) ! entrees supplementaires
      endif                 !m°v°m> !28-02-23

      if(flag_buo_w==4) then !m°v°m> !07-03-23
! Cas ou la vitesse de chute se deduit de:
! Effects of Particle Properties on the Settling and Rise Velocities of
! Microplastics in Freshwater under Laboratory Conditions
!      DOI: 10.1021/acs.est.8b06794
!      Environ. Sci. Technol. 2019, 53, 1958¿1966
       id_csf    =id_fieldrifter(1) ! EQ1: Corey Shape Factor
       id_dequi  =id_fieldrifter(2) ! EQ4: Equivalent particle diameter
       id_density=id_fieldrifter(3) ! particule density
      endif                  !m°v°m>
      

      read(3,*)drifter_random_w        !22-11-19
      read(3,*)flag_stay_in_water      !21-09-20
      if(flag_stay_in_water==1) then
       deck_drf_min=0.51
      else
       deck_drf_min=0.
      endif
      deck_drf_max=kmax+0.49

      if(rungekutta_order==1.and.drifter_random_w/=0) then !-warning->
       write(10+par%rank,*) &
       'Without a high-order scheme, random motion has the disadvantage' &
      ,' of clustering drifters at the edges of turbulent vertical' &
      ,' layers. If you are really sure of your choice, replace' &
      ,' rungekutta_order (currently 1 in notebook_drifter) by -1 ' &
      ,' However, I advise you to choose RK2'
       stop 'rungekutta_order==1.and.drifter_random_w/=0 see fort.xxx'
      endif                                                !-warning->
      if(rungekutta_order==-1)rungekutta_order=1

      if(drifter_initial_mode(1:7)=='fortran') then !------->
           call drifter_initial_fortran
           goto 488
      endif                                         !------->

      read(3,*)                                                        !07/10/02

  369 continue

! Initial location of drifters:
!     read(3,*,end=363)posx_,posy_,posz_,k0
!     read(3,*,end=363)posx_,posy_,posz_,k0,x10 ! bidouille patrick oU x10 est vitesse de la particule
      read(3,*,end=363)posx_,posy_,posz_,k_  
      if(flag_buo_w<=1) then !m°v°m> !28-02-23
        read(3,*)vertical_buoyancy_velocity       ! vertical_buoyancy_velocity est vitesse de la particule
      else                   !m°v°m> !28-02-23
       do k=1,flag_buo_w-1
! Note: si flag_buo_w=4 alors lire les parametres dans cet ordre svp:
!      id_csf     ! val_fieldrifter(1)=Corey Shape Factor
!      id_dequi   ! val_fieldrifter(2)=Equivalent particle diameter
!      id_density ! val_fieldrifter(3)=particule density
        read(3,*)val_fieldrifter(k)
       enddo
      endif                  !m°v°m> !28-02-23


! Initial time of drifters:
      read(3,*)i1,i2,i3,i4,i5,i6 !An,Mois,Jour,Heure,Minute,Seconde !03-03-15
      call datetokount(i1,i2,i3,i4,i5,i6) !An,Mois,Jour,Heure,Minute,Seconde

      kbuglb_   =kbuglb_   +1

      if(k_==0)then !§§§§§§§§> ! position en lat. lon. prof(<0)
         latit1=posx_*pi/180.
         longi1=posy_*pi/180.
         call latlon_to_ij('loc')                                     !21-10-10
      else          !§§§§§§§§> ! position i j k
! conversion des indices globaux en indices locaux
       deci=posx_-par%timax(1) ;  decj=posy_-par%tjmax(1)
      endif         !§§§§§§§§> ! position i j k

      if(deci>deci_drf_min.and.deci<deci_drf_max.and.  & ! test: la drifter est elle dans le domaine?
         decj>decj_drf_min.and.decj<decj_drf_max) then   !>>>>>>>>>>>>>>>>> !30-11-15

         in_or_out_   =1  ! 1=in
! Si oui, incrementer le compteur de drifter:
         kbu=kbu+1
         kbumax=kbumax+1
         if(kbu>nbomax)call drifter_error_message(1)

! Saisir la position en indices globaux dans tableau drifter_l:
         drifter_l(kbu,1)=deci+par%timax(1)
         drifter_l(kbu,2)=decj+par%tjmax(1)

        if(k_==0)then !-drifter_l(kbu,3)->
! Convertir la profondeur (m) en indice k:
         i=int(deci)
         j=int(decj)
         rapi=deci-real(i)
         rapj=decj-real(j)
         do k=1,kmax-1
          z1=(1.-rapi)*(1.-rapj)*depth_t(i  ,j  ,k)     &
            +(1.-rapi)*    rapj *depth_t(i  ,j+1,k)     &
            +    rapi *(1.-rapj)*depth_t(i+1,j  ,k)     &
            +    rapi *    rapj *depth_t(i+1,j+1,k)
          z2=(1.-rapi)*(1.-rapj)*depth_t(i  ,j  ,k+1)   &
            +(1.-rapi)*    rapj *depth_t(i  ,j+1,k+1)   &
            +    rapi *(1.-rapj)*depth_t(i+1,j  ,k+1)   &
            +    rapi *    rapj *depth_t(i+1,j+1,k+1)
            if(posz_>=z1.and.posz_<=z2) then  !>>>>>>
              drifter_l(kbu,3)=k+(posz_-z1)/(z2-z1)
            endif                       !>>>>>>
                 if(k==1.and.posz_<=z1)drifter_l(kbu,3)=1
            if(k==kmax-1.and.posz_>=z2)drifter_l(kbu,3)=kmax
         enddo
        else          !-drifter_l(kbu,3)->
         drifter_l(kbu,3)=posz_
        endif         !-drifter_l(kbu,3)->

         drifter_l(kbu,4)=real( kbuglb_    )

         drifter_l(kbu,5)=elapsedtime_out !03-03-15

! note: meme si flag_buo_w>1 l'emplacement id_wdrifter servira a stocker Ws
         if(flag_buo_w==1)drifter_l(kbu,id_wdrifter)=vertical_buoyancy_velocity ! bidouille patrick vitesse de la particule 

         if(flag_buo_w>1) then !>>> !28-02-23
          do k=1,flag_buo_w-1
           drifter_l(kbu,id_fieldrifter(k))=val_fieldrifter(k)
          enddo
         endif                 !>>>

      endif                                 !>>>>>>>>>>>>>>>>>

      goto 369

 363  close(3)
      kbumax_glb=kbuglb_ !24-10-22

!     kbumax_glb=kbumax_glb+nbdom !+nbdom car les proc vides ecriront dans une zone "poubelle"

!     if(par%rank==154) then
!      drifter_l(kbu,2)=real(jmax)+0.6+par%tjmax(1)
!      do kbu=1,kbumax
!       write(154,*)drifter_l(kbu,1:3),' trouvemoi'
!      enddo
!     endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! INITIALISATION DES BOUEES                                            !22/03/02
! FIN.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  488 continue

! Empecher les drifters d'Etre totalement dans la couche fusionnee:  23-11-19
      do kbu=1,kbumax
         i=int(drifter_l(kbu,1))-par%timax(1)
         j=int(drifter_l(kbu,2))-par%tjmax(1)
         if(i>=0.and.i<=imax+1.and.j>=0.and.j<=jmax+1) &
               drifter_l(kbu,3)=min(             & !ooo>
           max(drifter_l(kbu,3),real(min(kmin_w(i,j),kmin_w(i+1,j),kmin_w(i,j+1),kmin_w(i+1,j+1)))) & !23-11-19
                                    ,real(kmax))   !ooo>
      enddo

! Stopper le modele si la numerotation dans drifter_l(:,4) s'approche trop de la limite de !12-09-22
! precision necessaire A distinguer n et n+1,  a priori 16 777 216, qu'on va anticiper en se 
! limitant A 9 999 999, et recommander d'utiliser en plus drifter_l(:,5) pour prolonger au delA la numerotation
      flag_stop=0
      if(maxval(drifter_l(:,4))>9999999.)flag_stop=1 !12-09-22
      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) 
      if(k0/=0) then !debug>
        write(10+par%rank,*) &
                   'WARNING: The numbering in drifter_l(:,4) exceeds' &
      ,' 9 999 999 and is approaching the limit of precision allowed' &
      ,' with real4 (i.e. 16 777 216). ' &
      ,' Eventually complete with drifter_l(:,5) to extend the' &
      ,' numbering. Cancel this stop to continue.'  
        stop ' Warning in drifter_initial, see fort.xxx files'
      endif          !debug>

! mpi:
      call drifter_exchange_driver

      if(drifter_output_files==2)call drifter_netcdf_file(0) !24-10-22

      end subroutine drifter_initial

!..................................................................................

      subroutine drifter_error_message(choice_   )
!     use module_principal
      implicit none
      integer choice_

      if(choice_   ==1) then !111111111>
       write(6,*)'vous demandez plus de drifter qu''il en est déclaré'
       write(6,*)'dans parameter, actuellement dimensionné à ',nbomax
       write(6,*)'mettre nbomax à au moins ',kbu
       stop ' STOP dans subroutine drifter'
      endif                  !111111111>

      if(choice_   ==2) then !222222222>
       write(6,*)'nbobuffermax trop petit. Modifiez parameter!'
       write(6,*)'nbobuffermax dans parameter=',nbobuffermax
       write(6,*)'Choisir max de:',nd_send_est,nd_send_ouest    &
                                  ,nd_send_nord,nd_send_sud,nd_send_out
       stop ' STOP dans subroutine drifter'
      endif                  !222222222>

      if(choice_   ==3) then !333333333>
       write(6,*)'vous demandez plus de drifter qu''il en est déclaré'
       write(6,*)'dans parameter, actuellement dimensionné à ',nbomax
       write(6,*)'mettre nbomax à au moins ',new_kbumax
       stop ' STOP dans subroutine drifter_error_message cas 3'
      endif                  !333333333>

      if(choice_   ==4) then !444444444>
       write(6,*)'nbobuffermax_c trop petit. Modifiez parameter!' !15-11-20
       write(6,*)'nbobuffermax_c dans parameter=',nbobuffermax_c
       write(6,*)'Choisir au moins:',max(nd_send_sudouest,nd_send_sudest,nd_send_nordouest,nd_send_nordest)
       stop ' STOP dans subroutine drifter_error_message cas 4'
      endif                  !444444444>

      end subroutine drifter_error_message

!..................................................................................

      subroutine drifter_initial_fortran
      use module_principal
      use module_parallele
      implicit none
      integer loop_
      character*60 txt_spy_

      txt_spy_   ='drifter_initial_fortran'

      kbu=0
      do j=2,jmax-1
      do i=2,imax-1
       if(mask_t(i,j,kmaxp1)==1) then
!      if(mod(i+par%timax(1),3)==0.and.mod(j+par%tjmax(1),3)==0) then
        kbu=kbu+1
!      endif

       endif
      enddo
      enddo
      call barriere(iteration3d,4,txt_spy_   )                  !27-05-11
      if(kbu>nbomax)call drifter_error_message(1)

      kbu=0
      do j=2,jmax-1
      do i=2,imax-1
       if(mask_t(i,j,kmaxp1)==1) then
        do loop_=-50,50
         kbu=kbu+1
!        call random_number(x1_r4) ; x1_r4=0.99*(x1_r4-0.5)
         drifter_l(kbu,1)=real(i)+par%timax(1) !+x1_r4
!        call random_number(x1_r4) ; x1_r4=0.99*(x1_r4-0.5)
         drifter_l(kbu,2)=real(j)+par%tjmax(1) !+x1_r4
!        drifter_l(kbu,3)=real(kmax)
!        drifter_l(kbu,3)=real(25)
         call random_number(x1_r4) ; x1_r4=0.99*(x1_r4-0.5)
         drifter_l(kbu,3)=real(kmin_w(i,j)+2)+x1_r4
!        drifter_l(kbu,5)=-0.0 
        enddo !loop
!       if(i==imax/2.and.j==jmax/2) then
!        kbu=kbu+1
!        drifter_l(kbu,1)=real(i)+par%timax(1)
!        drifter_l(kbu,2)=real(j)+par%tjmax(1)
!        drifter_l(kbu,3)=real(kmax)
!        drifter_l(kbu,6)=0. ! -0.01
!       endif
       endif
      enddo
      enddo
      kbumax=kbu

! Numero d'identifiant drifter_l(kbu,4)
! On numerote en commencant par le rank 0 par ordre croissant de kbu
! on poursuit avec rank 1 A partir de la derniere valeur et ainsi de
! suite....
      k10=0 
      do loop1=0,nbdom-1
      if(par%rank==loop1) then !rank>
       do kbu=1,kbumax
        k10=k10+1
        drifter_l(kbu,4)=k10
       enddo
       if(kbumax>0) then   !message>
        write(6,*)'par%rank,kbumax,drifter_l(kbumax,4)',par%rank,kbumax,drifter_l(kbumax,4) 
       endif               !message>
      endif                    !rank>
! envoyer aux autres rank la derniere valeur de la numerotation (max de k10)
      k0=k10
      call mpi_allreduce(k0,k10,1,mpi_integer,mpi_max,par%comm2d ,ierr)

      enddo ! loop1

      kbumax_glb=k10 !24-10-22
!     kbumax_glb=kbumax_glb+nbdom !+nbdom car les proc vides ecriront dans une zone "poubelle"

      end subroutine drifter_initial_fortran

!..................................................................................
      subroutine drifter_who_is_out
!     use module_principal
!     use module_parallele
      implicit none
      integer,parameter :: tagouest_   =5000,  tagest_   =5010   &
                            ,tagsud_   =6000, tagnord_   =6010
      integer,parameter :: tagsudouest_   =15000, tagsudest_   =15010  &
                         ,tagnordouest_   =16000, tagnordest_   =16010 &
                         ,nexchgmax_   =16
      integer nexchg_   ,warning_
      integer,dimension(nexchgmax_   ) :: tabreq_
      integer,dimension(mpi_status_size,nexchgmax_   ) :: tstatus_
      character*60 txt_spy_

! Cette routine recense les numeros des drifters qui doivent migrer vers un
! domaine voisin:
      txt_spy_   ='drifter_who_is_out'

! Etape 1: identifier les drifters sortant
! et produire nd_send_nord nd_send_sud nd_send_est nd_send_ouest
 1012 continue  ! point d'envoi du goto 1012

      warning_   =0

      nd_send_sudouest=0
      nd_send_sudest=0
      nd_send_nordouest=0
      nd_send_nordest=0

      nd_send_nord=0
      nd_send_sud=0
      nd_send_est=0
      nd_send_ouest=0
      nd_send_out=0

      do kbu=1,kbumax  ! kbu loop

! Indices locaux:
      deci=drifter_l(kbu,1)-par%timax(1)
      decj=drifter_l(kbu,2)-par%tjmax(1)

! la zone d'echange Coin Nord-Est
      if(deci>deci_drf_max.and.decj>decj_drf_max) then !nenenene>

      if (par%tvoisin(nordest) /= mpi_proc_null) then !-mpn-mpn->

       nd_send_nordest=nd_send_nordest+1
        if(nd_send_nordest<=nbobuffermax_c) then !>>>>>>
         drifter_send_order_nordest(nd_send_nordest)=kbu
         drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
! note: je fais -abs(drifter...) plutot que -drifter... pour parer au fait que le calcul peut etre refait suite au goto 1012 (et retablirait un signe positif errone)
        endif                                   !>>>>>>

      else                                        !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out<=nbobuffermax) then !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
        endif                             !>>>>>>

      endif                                       !-mpn-mpn->

      goto 1968
      endif                                        !nenenene>

! la zone d'echange Coin Nord-Ouest
      if(deci<deci_drf_min.and.decj>decj_drf_max) then !nonononono>

      if (par%tvoisin(nordouest) /= mpi_proc_null) then !-mpn-mpn->

       nd_send_nordouest=nd_send_nordouest+1
        if(nd_send_nordouest<=nbobuffermax_c) then !>>>>>>
         drifter_send_order_nordouest(nd_send_nordouest)=kbu
         drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
        endif                                   !>>>>>>

      else                                              !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out<=nbobuffermax) then !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
        endif                             !>>>>>>

      endif                                             !-mpn-mpn->

      goto 1968
      endif                               !nonononono>

! la zone d'echange Coin Sud-Ouest
      if(deci<deci_drf_min.and.decj<decj_drf_min) then !sososososo>

      if (par%tvoisin(sudouest) /= mpi_proc_null) then !-mpn-mpn->

       nd_send_sudouest=nd_send_sudouest+1
        if(nd_send_sudouest<=nbobuffermax_c) then !>>>>>>
         drifter_send_order_sudouest(nd_send_sudouest)=kbu
         drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
        endif                                   !>>>>>>

      else                                              !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out<=nbobuffermax) then !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
        endif                             !>>>>>>

      endif                                             !-mpn-mpn->

      goto 1968
      endif                       !sososososo>

! la zone d'echange Coin Sud-est
      if(deci>deci_drf_max.and.decj<decj_drf_min) then !sesesesese>

      if (par%tvoisin(sudest) /= mpi_proc_null) then !-mpn-mpn->

       nd_send_sudest=nd_send_sudest+1
        if(nd_send_sudest<=nbobuffermax_c) then !>>>>>>
         drifter_send_order_sudest(nd_send_sudest)=kbu
         drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
        endif                                   !>>>>>>

      else                                              !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out<=nbobuffermax) then !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
        endif                             !>>>>>>

      endif                                             !-mpn-mpn->

      goto 1968
      endif                       !sesesesese>

! la zone d'echange Est:
      if(deci>deci_drf_max) then !eeeeeee>
      if (par%tvoisin(est) /= mpi_proc_null) then !-mpn-mpn->

        nd_send_est=nd_send_est+1
        if(nd_send_est<=nbobuffermax) then !>>>>>>
         drifter_send_order_est(nd_send_est)=kbu
         drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
        endif                             !>>>>>>

      else                                        !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out<=nbobuffermax) then !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
        endif                             !>>>>>>

      endif                                       !-mpn-mpn->
      goto 1968
      endif                    !eeeeeee>


! la zone d'echange Ouest:
      if(deci<deci_drf_min)        then !ooooooo>
      if (par%tvoisin(ouest) /= mpi_proc_null) then !-mpn-mpn->

        nd_send_ouest=nd_send_ouest+1
        if(nd_send_ouest<=nbobuffermax) then !>>>>>>>
          drifter_send_order_ouest(nd_send_ouest)=kbu
          drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
        endif                               !>>>>>>>

      else                                        !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out<=nbobuffermax) then !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
        endif                             !>>>>>>

      endif                                       !-mpn-mpn->
      goto 1968
      endif                    !ooooooo>


! la zone d'echange Nord:
      if(decj>decj_drf_max) then !nnnnnnnn>
      if (par%tvoisin(nord) /= mpi_proc_null) then !-mpn-mpn->

        nd_send_nord=nd_send_nord+1
        if(nd_send_nord<=nbobuffermax) then !>>>>>>>>
          drifter_send_order_nord(nd_send_nord)=kbu
          drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
        endif                              !>>>>>>>>

      else                                        !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out<=nbobuffermax) then !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
        endif                             !>>>>>>

      endif                                       !-mpn-mpn->
      goto 1968
      endif                    !nnnnnnnn>

! la zone d'echange Sud:
      if(decj<decj_drf_min)         then !ssssssss>
      if (par%tvoisin(sud) /= mpi_proc_null) then !-mpn-mpn->

        nd_send_sud=nd_send_sud+1
        if(nd_send_sud<=nbobuffermax) then !>>>>>>>>>
          drifter_send_order_sud(nd_send_sud)=kbu
          drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
        endif                             !>>>>>>>>>

      else                                        !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out<=nbobuffermax) then !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-abs(drifter_l(kbu,4))
        endif                             !>>>>>>

      endif                                       !-mpn-mpn->
      goto 1968
      endif                    !ssssssss>

 1968 continue    ! SORTIE

      enddo       ! kbu loop

!.............................................................................................
! Etape 2:
! Les domaines indiquent à leurs voisins le nombre de drifters s'appretant à traverser les
! frontieres. L'objectif est qu'en routine "obc" on n'echange que le nombre vrai de drifters
! migrant d'un domaine à l'autre et non pas l'integralité des tableaux "buffer" surdimensionné
! à nbobuffermax
      nexchg_   =0
! Coin Nord-Est:
      if (par%tvoisin(nordest) /= mpi_proc_null) then
!     if (par%tvoisin(nord) /= mpi_proc_null) then
!     if (par%tvoisin(est)  /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nd_send_nordest,1,mpi_integer,par%tvoisin(nordest)  &
              ,tagnordest_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nd_recv_nordest,1,mpi_integer,par%tvoisin(nordest)  &
              ,tagsudouest_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
!     endif
      endif
! Coin Nord-Ouest:
      if (par%tvoisin(nordouest) /= mpi_proc_null) then
!     if (par%tvoisin(nord)  /= mpi_proc_null) then
!     if (par%tvoisin(ouest) /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nd_send_nordouest,1,mpi_integer,par%tvoisin(nordouest)  &
              ,tagnordouest_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nd_recv_nordouest,1,mpi_integer,par%tvoisin(nordouest)  &
              ,tagsudest_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
      endif
! Coin Sud-Ouest:
      if (par%tvoisin(sudouest) /= mpi_proc_null) then
!     if (par%tvoisin(sud  ) /= mpi_proc_null) then
!     if (par%tvoisin(ouest) /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nd_send_sudouest,1,mpi_integer,par%tvoisin(sudouest)  &
              ,tagsudouest_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nd_recv_sudouest,1,mpi_integer,par%tvoisin(sudouest)  &
              ,tagnordest_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
      endif
! Coin Sud-Est:
      if (par%tvoisin(sudest) /= mpi_proc_null) then
!     if (par%tvoisin(sud) /= mpi_proc_null) then
!     if (par%tvoisin(est) /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nd_send_sudest,1,mpi_integer,par%tvoisin(sudest)  &
              ,tagsudest_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nd_recv_sudest,1,mpi_integer,par%tvoisin(sudest)  &
              ,tagnordouest_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
      endif

! Frontiere Est:
      if (par%tvoisin(est) /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nd_send_est,1,mpi_integer,par%tvoisin(est)  &
              ,tagest_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nd_recv_est,1,mpi_integer,par%tvoisin(est)  &
              ,tagouest_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
      endif
! Frontiere Ouest:
      if (par%tvoisin(ouest) /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nd_send_ouest,1,mpi_integer,par%tvoisin(ouest)  &
              ,tagouest_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nd_recv_ouest,1,mpi_integer,par%tvoisin(ouest)  &
              ,tagest_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
      endif
! Frontiere Nord:
      if (par%tvoisin(nord) /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nd_send_nord,1,mpi_integer,par%tvoisin(nord)  &
              ,tagnord_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nd_recv_nord,1,mpi_integer,par%tvoisin(nord)  &
              ,tagsud_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
      endif
! Frontiere sud:
      if (par%tvoisin(sud) /= mpi_proc_null) then
       nexchg_   =nexchg_   +1
       call mpi_issend(nd_send_sud,1,mpi_integer,par%tvoisin(sud)  &
              ,tagsud_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
       nexchg_   =nexchg_   +1
       call mpi_irecv(nd_recv_sud,1,mpi_integer,par%tvoisin(sud)  &
              ,tagnord_   ,par%comm2d,tabreq_   (nexchg_   ),ierr)
      endif

      if(nexchg_   >nexchgmax_   ) &
      stop 'drifter_who_is_out nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_   (:,1:nexchg_   ) &
                      ,ierr)

! A l'issu de ce calcul on verifie que:
! ( par%rank , nd_send_est ) = ( par%voisin(ouest) , nc_recv_ouest )
! etc....
! Ces parametres vont permettre de definir la "size" des tableaux à echanger dans l'etape suivante



!-!-!-!-!-!-!-!-!-!-!-!>
! VERIFICATION DE LA TAILLE MEMOIRE DES TABLEAUX ET REALLOCATION SI NECESSAIRE
      if(nbobuffermax<max(nd_send_est,nd_send_ouest,nd_send_nord,nd_send_sud,nd_send_out &
                         ,nd_recv_est,nd_recv_ouest,nd_recv_nord,nd_recv_sud)) then !m°v°m>

      warning_=1
! Memoire insuffisante detectee pour la buffer zone. Reallouer les tableaux avec la dimension necessaire
      nbobuffermax=max(nd_send_est,nd_send_ouest,nd_send_nord,nd_send_sud,nd_send_out &
                      ,nd_recv_est,nd_recv_ouest,nd_recv_nord,nd_recv_sud)
       call drifter_allocate(-1) ! desallouer les tableaux buffer
       call drifter_allocate(1)  ! reallouer  avec la nouvelle valeur de nbobuffermax
       write(6,*)'rank ',par%rank &
                ,' realloue buffer zone avec nbobuffermax=',nbobuffermax
      endif                   !m°v°m>

       if(nbobuffermax_c<max(nd_send_sudouest,nd_send_sudest,nd_send_nordouest,nd_send_nordest &
                            ,nd_recv_sudouest,nd_recv_sudest,nd_recv_nordouest,nd_recv_nordest)) then !m°O°m>

       warning_=2
! Memoire insuffisante detectee pour les coins de la buffer zone. Reallouer les tableaux avec la dimension necessaire
       nbobuffermax_c=max(nd_send_sudouest,nd_send_sudest,nd_send_nordouest,nd_send_nordest &
                         ,nd_recv_sudouest,nd_recv_sudest,nd_recv_nordouest,nd_recv_nordest)
       call drifter_allocate(-2) ! desallouer les tableaux coins buffer
       call drifter_allocate(2)  ! reallouer avec la nouvelle valeur de nbobuffermax_c
       write(6,*)'rank ',par%rank &
                ,' realloue coin buffer zone avec nbobuffermax_c=',nbobuffermax_c
      endif                   !m°O°m>

! Partager le warning avec tous les ranks puis aviser
      call mpi_allreduce(warning_,k0,1,mpi_integer,mpi_max,par%comm2d,ierr)
      warning_=k0
      if(warning_/=0)goto 1012 ! REFAIRE LE CALCUL AVEC LES NOUVELLES DIMENSIONS
!-!-!-!-!-!-!-!-!-!-!-!>

!     call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
      end subroutine drifter_who_is_out
!..................................................................................
      subroutine drifter_obc(var_   )
!     use module_principal
!     use module_parallele
      implicit none
      integer,parameter :: tagouest_   =5000,  tagest_   =5010   &
                            ,tagsud_   =6000, tagnord_   =6010
      integer,parameter :: tagsudouest_   =15000, tagsudest_   =15010  &
                         ,tagnordouest_   =16000, tagnordest_   =16010 &
                         ,nexchgmax_   =16
!     integer nexchg_   ,warning_   ,var_
      integer nexchg_   ,var_,loop_
      integer,dimension(nexchgmax_   ) :: tabreq_
      integer,dimension(mpi_status_size,nexchgmax_   ) :: tstatus_
      character*60 txt_spy_

      txt_spy_   ='drifter_obc'

!.............................................................................................
! Etape 1:
! Charger les tableaux d'echanges avec les drifters positionnes sur les bords des sous domaines:
      do loop_=1,drifter_dim2 ! loop_

      do k1=1,nd_send_nordest
       k2=drifter_send_order_nordest(k1)
!      k3=k1+(loop_-1)*nd_send_nordest
       k3=loop_+(k1-1)*drifter_dim2
       drifter_send_nordest(k3)=drifter_l(k2,loop_)
      enddo
      do k1=1,nd_send_nordouest
       k2=drifter_send_order_nordouest(k1)
!      k3=k1+(loop_-1)*nd_send_nordouest
       k3=loop_+(k1-1)*drifter_dim2
       drifter_send_nordouest(k3)=drifter_l(k2,loop_)
      enddo
      do k1=1,nd_send_sudouest
       k2=drifter_send_order_sudouest(k1)
!      k3=k1+(loop_-1)*nd_send_sudouest
       k3=loop_+(k1-1)*drifter_dim2
       drifter_send_sudouest(k3)=drifter_l(k2,loop_)
      enddo
      do k1=1,nd_send_sudest
       k2=drifter_send_order_sudest(k1)
!      k3=k1+(loop_-1)*nd_send_sudest
       k3=loop_+(k1-1)*drifter_dim2
       drifter_send_sudest(k3)=drifter_l(k2,loop_)
      enddo

      do k1=1,nd_send_est
       k2=drifter_send_order_est(k1)
!      k3=k1+(loop_-1)*nd_send_est
       k3=loop_+(k1-1)*drifter_dim2
       drifter_send_est(k3)=drifter_l(k2,loop_)
      enddo
      do k1=1,nd_send_ouest
       k2=drifter_send_order_ouest(k1)
!      k3=k1+(loop_-1)*nd_send_ouest
       k3=loop_+(k1-1)*drifter_dim2
       drifter_send_ouest(k3)=drifter_l(k2,loop_)
      enddo
      do k1=1,nd_send_nord
       k2=drifter_send_order_nord(k1)
!      k3=k1+(loop_-1)*nd_send_nord
       k3=loop_+(k1-1)*drifter_dim2
       drifter_send_nord(k3)=drifter_l(k2,loop_)
      enddo
      do k1=1,nd_send_sud
       k2=drifter_send_order_sud(k1)
!      k3=k1+(loop_-1)*nd_send_sud
       k3=loop_+(k1-1)*drifter_dim2
       drifter_send_sud(k3)=drifter_l(k2,loop_)
      enddo

      enddo ! loop_ 

!.............................................................................................
! Etape 2:
! Echanger les informations des zones d'echanges:

      nexchg_   =0

! Coin Nord-Est:
      if (par%tvoisin(nordest) /= mpi_proc_null) then
!     if (par%tvoisin(nord) /= mpi_proc_null) then
!     if (par%tvoisin(est)  /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nd_send_nordest*drifter_dim2)
      call mpi_issend(drifter_send_nordest(1)    & ! envoyé au proc Est
                ,size(drifter_send_nordest(1:k)) &
                ,mpi_real                        &
                ,par%tvoisin(nordest)            &
                ,tagnordest_                     &
                ,par%comm2d                      &
                ,tabreq_   (nexchg_   )          &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nd_recv_nordest*drifter_dim2)
      call mpi_irecv(drifter_recv_nordest(1)     & ! tableau recu
               ,size(drifter_recv_nordest(1:k))  &
               ,mpi_real                         &
               ,par%tvoisin(nordest)             &
               ,tagsudouest_                     &
               ,par%comm2d                       &
               ,tabreq_   (nexchg_   )           &
               ,ierr)
      endif

! Coin Nord-Ouest:
      if (par%tvoisin(nordouest) /= mpi_proc_null) then
!     if (par%tvoisin(nord ) /= mpi_proc_null) then
!     if (par%tvoisin(ouest) /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nd_send_nordouest*drifter_dim2)
      call mpi_issend(drifter_send_nordouest(1)    & ! envoyé au proc Est
                ,size(drifter_send_nordouest(1:k)) &
                ,mpi_real                          &
                ,par%tvoisin(nordouest)            &
                ,tagnordouest_                     &
                ,par%comm2d                        &
                ,tabreq_   (nexchg_   )            &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nd_recv_nordouest*drifter_dim2)
      call mpi_irecv(drifter_recv_nordouest(1)     & ! tableau recu
               ,size(drifter_recv_nordouest(1:k))  &
               ,mpi_real                           &
               ,par%tvoisin(nordouest)             &
               ,tagsudest_                         &
               ,par%comm2d                         &
               ,tabreq_   (nexchg_   )             &
               ,ierr)
      endif

! Coin Sud-Ouest:
      if (par%tvoisin(sudouest) /= mpi_proc_null) then
!     if (par%tvoisin(sud  ) /= mpi_proc_null) then
!     if (par%tvoisin(ouest) /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nd_send_sudouest*drifter_dim2)
      call mpi_issend(drifter_send_sudouest(1)    & ! envoyé au proc Est
                ,size(drifter_send_sudouest(1:k)) &
                ,mpi_real                         &
                ,par%tvoisin(sudouest)            &
                ,tagsudouest_                     &
                ,par%comm2d                       &
                ,tabreq_   (nexchg_   )           &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nd_recv_sudouest*drifter_dim2)
      call mpi_irecv(drifter_recv_sudouest(1)     & ! tableau recu
               ,size(drifter_recv_sudouest(1:k))  &
               ,mpi_real                          &
               ,par%tvoisin(sudouest)             &
               ,tagnordest_                       &
               ,par%comm2d                        &
               ,tabreq_   (nexchg_   )            &
               ,ierr)
      endif

! Coin Sud-Est:
      if (par%tvoisin(sudest) /= mpi_proc_null) then
!     if (par%tvoisin(sud) /= mpi_proc_null) then
!     if (par%tvoisin(est) /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nd_send_sudest*drifter_dim2)
      call mpi_issend(drifter_send_sudest(1)    & ! envoyé au proc Est
                ,size(drifter_send_sudest(1:k)) &
                ,mpi_real                       &
                ,par%tvoisin(sudest)            &
                ,tagsudest_                     &
                ,par%comm2d                     &
                ,tabreq_   (nexchg_   )         &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nd_recv_sudest*drifter_dim2)
      call mpi_irecv(drifter_recv_sudest(1)     & ! tableau recu
               ,size(drifter_recv_sudest(1:k))  &
               ,mpi_real                        &
               ,par%tvoisin(sudest)             &
               ,tagnordouest_                   &
               ,par%comm2d                      &
               ,tabreq_   (nexchg_   )          &
               ,ierr)
      endif


! Frontiere Est:
      if (par%tvoisin(est) /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nd_send_est*drifter_dim2)
!     call mpi_issend(drifter_send_est(1:k)  & ! envoyé au proc Est
      call mpi_issend(drifter_send_est(1)    & ! envoyé au proc Est
                ,size(drifter_send_est(1:k)) &
                ,mpi_real                       &
                ,par%tvoisin(est)                           &
                ,tagest_                                    &
                ,par%comm2d                                 &
                ,tabreq_   (nexchg_   )                     &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nd_recv_est*drifter_dim2)
!     call mpi_irecv(drifter_recv_est(1:k)   & ! tableau recu du proc Ouest
      call mpi_irecv(drifter_recv_est(1)     & ! tableau recu du proc Ouest
               ,size(drifter_recv_est(1:k))  &
               ,mpi_real                        &
               ,par%tvoisin(est)                            &
               ,tagouest_                                   &
               ,par%comm2d                                  &
               ,tabreq_   (nexchg_   )                      &
               ,ierr)
      endif

! Frontiere Ouest:
      if (par%tvoisin(ouest) /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nd_send_ouest*drifter_dim2)
!     call mpi_issend(drifter_send_ouest(1:k)  &
      call mpi_issend(drifter_send_ouest(1)    &
                ,size(drifter_send_ouest(1:k)) &
                ,mpi_real                         &
                ,par%tvoisin(ouest)                           &
                ,tagouest_                                    &
                ,par%comm2d                                   &
                ,tabreq_   (nexchg_   )                       &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nd_recv_ouest*drifter_dim2)
!     call mpi_irecv(drifter_recv_ouest(1:k)  &
      call mpi_irecv(drifter_recv_ouest(1)    &
               ,size(drifter_recv_ouest(1:k)) &
               ,mpi_real                         &
               ,par%tvoisin(ouest)                           &
               ,tagest_                                      &
               ,par%comm2d                                   &
               ,tabreq_   (nexchg_   )                       &
               ,ierr)
      endif

! Frontiere nord:
      if (par%tvoisin(nord) /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nd_send_nord*drifter_dim2)
!     call mpi_issend(drifter_send_nord(1:k)  &
      call mpi_issend(drifter_send_nord(1)    &
                ,size(drifter_send_nord(1:k)) &
                ,mpi_real                        &
                ,par%tvoisin(nord)                           &
                ,tagnord_                                    &
                ,par%comm2d                                  &
                ,tabreq_   (nexchg_   )                      &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nd_recv_nord*drifter_dim2)
!     call mpi_irecv(drifter_recv_nord(1:k)  &
      call mpi_irecv(drifter_recv_nord(1)    &
               ,size(drifter_recv_nord(1:k)) &
               ,mpi_real                        &
               ,par%tvoisin(nord)                           &
               ,tagsud_                                     &
               ,par%comm2d                                  &
               ,tabreq_   (nexchg_   )                      &
               ,ierr)
      endif

! Frontiere sud:
      if (par%tvoisin(sud) /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nd_send_sud*drifter_dim2)
      call mpi_issend(drifter_send_sud(1)     &
                ,size(drifter_send_sud(1:k))  &
                ,mpi_real                        &
                ,par%tvoisin(sud)                            &
                ,tagsud_                                     &
                ,par%comm2d                                  &
                ,tabreq_   (nexchg_   )                      &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nd_recv_sud*drifter_dim2)
      call mpi_irecv(drifter_recv_sud(1)     &
               ,size(drifter_recv_sud(1:k))  &
               ,mpi_real                        &
               ,par%tvoisin(sud)                            &
               ,tagnord_                                    &
               ,par%comm2d                                  &
               ,tabreq_   (nexchg_   )                      &
               ,ierr)
      endif

      if(nexchg_   >nexchgmax_   ) &
      stop 'drifter_obc nexchg_   >nexchgmax_   '

      ! La salle d'attente
      call mpi_waitall(nexchg_                     &
                      ,tabreq_   (1:nexchg_   )    &
                      ,tstatus_   (:,1:nexchg_   ) &
                      ,ierr)


      call barriere(iteration3d,4,txt_spy_   )                  !27-05-11

      end subroutine drifter_obc
!....................................................................

      subroutine drifter_exchange_driver
!     use module_principal
      implicit none
      integer loop_

      call drifter_who_is_out ! recense les drifter à sortir

      call drifter_obc(0) !04-09-22

      if(nbcanal>0)call drifter_webcanals !04-09-22

      call drifter_update_numbering

! archivage des trajectoires dans des fichiers
      if(drifter_out_sampling>0.) then !m°v°m>

      if(  int( elapsedtime_now        /drifter_out_sampling)       &
         /=int((elapsedtime_now-dt_drf)/drifter_out_sampling)) then !ooo>

        if(drifter_output_files<2) then !->
         call drifter_ascii_file
        else                            !->
         call drifter_netcdf_file(1) !24-10-22
        endif                           !->

      endif                                                         !ooo>

      endif                            !m°v°m>

      end subroutine drifter_exchange_driver
!.................................................................
      subroutine drifter_gnuplot_output
!     use module_principal
!     use module_parallele
      implicit none

      if(mod(iteration3d,50)==0) then
      do k=0,7
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
      if(k==par%rank) then

      if(par%rank==0) then
       open(unit=10,file='toto0')
       open(unit=11,file='toto1')
       open(unit=12,file='toto2')
       open(unit=13,file='toto3')
       open(unit=14,file='toto4')
       open(unit=15,file='toto5')
       open(unit=16,file='toto6')
       open(unit=17,file='toto7')
       write(10,*)-10,-10
       write(11,*)-10,-10
       write(12,*)-10,-10
       write(13,*)-10,-10
       write(14,*)-10,-10
       write(15,*)-10,-10
       write(15,*)-10,-10
       write(16,*)-10,-10
       write(17,*)-10,-10
      else
       open(unit=10,file='toto0',position='append')
       open(unit=11,file='toto1',position='append')
       open(unit=12,file='toto2',position='append')
       open(unit=13,file='toto3',position='append')
       open(unit=14,file='toto4',position='append')
       open(unit=15,file='toto5',position='append')
       open(unit=16,file='toto6',position='append')
       open(unit=17,file='toto7',position='append')
      endif

      do kbu=1,kbumax
      k1=int(drifter_l(kbu,4)/1000000)
      if(k1<0)stop 'fou1'
      if(k1>7)stop 'fou2'
      write(10+k1,*)drifter_l(kbu,1),drifter_l(kbu,2)
      enddo

      close(10)
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      close(17)
      endif
      call mpi_barrier(par%comm2d,k_out)      ! synchro processes      !16-09-09
      enddo
      endif

      end subroutine drifter_gnuplot_output
!..................................................................................

      subroutine drifter_ascii_file
      use module_principal
      use module_parallele
      implicit none

      if(kbumax==0)return

      if(drifter_output_files==1) then !rank-files>
        write(texte30,'(a,i0)')'tmp/drifters_rank',par%rank
        open(unit=3,file=trim(texte30),position='append')
      endif                            !rank-files>
      do 145 kbu=1,kbumax

! Deduire lon, lat, z des indices i,j,k:
      deci=drifter_l(kbu,1)-par%timax(1)
      decj=drifter_l(kbu,2)-par%tjmax(1)
      deck=drifter_l(kbu,3)

      i=int(deci)
      j=int(decj)
      k=int(deck)

      rapi=deci-real(i)
      rapj=decj-real(j)
      rapk=deck-real(k)

      k1=max0(k  ,1)
      k2=min0(k+1,kmax)

      x1=(1.-rapi)*(1.-rapj)*lon_t(i  ,j  )   &
        +(1.-rapi)*    rapj *lon_t(i  ,j+1)   &
        +    rapi *(1.-rapj)*lon_t(i+1,j  )   &
        +    rapi *    rapj *lon_t(i+1,j+1)

      x2=(1.-rapi)*(1.-rapj)*lat_t(i  ,j  )   &
        +(1.-rapi)*    rapj *lat_t(i  ,j+1)   &
        +    rapi *(1.-rapj)*lat_t(i+1,j  )   &
        +    rapi *    rapj *lat_t(i+1,j+1)

      x3=(1.-rapi)*(1.-rapj)*(1.-rapk)*depth_t(i  ,j  ,k1)   &
        +(1.-rapi)*(1.-rapj)*    rapk *depth_t(i  ,j  ,k2)   &
        +(1.-rapi)*    rapj *(1.-rapk)*depth_t(i  ,j+1,k1)   &
        +(1.-rapi)*    rapj *    rapk *depth_t(i  ,j+1,k2)   &
        +    rapi *(1.-rapj)*(1.-rapk)*depth_t(i+1,j  ,k1)   &
        +    rapi *(1.-rapj)*    rapk *depth_t(i+1,j  ,k2)   &
        +    rapi *    rapj *(1.-rapk)*depth_t(i+1,j+1,k1)   &
        +    rapi *    rapj *    rapk *depth_t(i+1,j+1,k2)


      if(drifter_output_files==0) then !individual-files>
        write(texte30,'(i8)')10000000+int(drifter_l(kbu,4))
        open(unit=3,file=trim(tmpdirname)//'drifter'//texte30(2:8),position='append')
        write(3,'(e12.5,6(1x,e14.7))')        & !11/03/09 !06-04-22
          elapsedtime_now/86400.              & ! temps en jours depuis le temps ref !15-04-11
         ,drifter_l(kbu,1)                    & ! position i
         ,drifter_l(kbu,2)                    & ! position j
         ,drifter_l(kbu,3)                    & ! position k
         ,x1*rad2deg                          & ! position longitude
         ,x2*rad2deg                          & ! position latitude
         ,x3                                    ! profondeur
        close(3)
      endif                            !individual-files>

      if(drifter_output_files==1) then !rank-files>
!       write(3,'(                     &
!                i4,1x,5(i2,1x)        & 
!               ,f8.2,1x               & 
!               ,i10                   & 
!               ,3(1x,f13.6),1x        &
!               ,f12.1,1x              & 
!               ,f9.6,1x               & 
!                )')                   &
!       year_now,month_now,day_now     &
!      ,hour_now,minute_now,second_now & ! date
!      ,elapsedtime_now/86400.         & ! temps EcoulE en jours
!      ,nint(drifter_l(kbu,4),kind=8)  & ! IDentite du drifter
!      ,x1*rad2deg,x2*rad2deg,x3       & ! longitude, latitude, profondeur
!      ,drifter_l(kbu,5)               & ! temps A la position initiale du drifter
!      ,drifter_l(kbu,id_wdrifter)       ! vitesse individuelle de flottabilitE
        write(texte30,'(a2,i0,a1)')' I',nint(drifter_l(kbu,4),kind=8),'D'
!'
        write(3,'(e12.5,3(1x,e14.7),a,1x)')   & !11/03/09 !06-04-22
          elapsedtime_now/86400.              & ! temps en jours depuis le temps ref !15-04-11
         ,drifter_l(kbu,1)                    & ! position i
         ,drifter_l(kbu,2)                    & ! position j
         ,drifter_l(kbu,3)                    & ! position k
         ,trim(texte30)
      endif                            !rank-files>

  145 continue

      if(drifter_output_files==1)close(3) !rank-files>

      end subroutine drifter_ascii_file

!.........................................................................................

      subroutine drifter_allocate(case_)
!     use module_principal
      implicit none
      integer case_

       if(case_==0) then !---> ! etat initial !16-11-20
        if(drifter_output_files/=2) then !ooo>
! sorties en fichiers ascii:
         allocate(drifter_l(nbomax,drifter_dim2));drifter_l=0. !03-03-15
        else                             !ooo>
! sorties en fichiers netcdf "facon Ariane": l'emplacement drifter_l(:,0) servira A stocker temporairement lon,lat,z
         allocate(drifter_l(nbomax,0:drifter_dim2));drifter_l=0.  !24-10-22
         allocate(kbumax_glb1(0:nbdom-1)) ; kbumax_glb1=0
         allocate(kbumax_glb2(0:nbdom-1)) ; kbumax_glb2=0
        endif                            !ooo>
       endif             !--->

       if(case_==0.or. &      !    etat initial !16-11-20
          case_==1) then !>>> ! ou reallocation en cours de run de la buffer zone 

        allocate(drifter_send_nord           (nbobuffermax*drifter_dim2)      )
        allocate(drifter_recv_nord           (nbobuffermax*drifter_dim2)      )
        allocate(drifter_send_sud            (nbobuffermax*drifter_dim2)      )
        allocate(drifter_recv_sud            (nbobuffermax*drifter_dim2)      )
        allocate(drifter_send_est            (nbobuffermax*drifter_dim2)      )
        allocate(drifter_recv_est            (nbobuffermax*drifter_dim2)      )
        allocate(drifter_send_ouest          (nbobuffermax*drifter_dim2)      )
        allocate(drifter_recv_ouest          (nbobuffermax*drifter_dim2)      )

        allocate(drifter_send_order_nord     (nbobuffermax)      )
        allocate(drifter_send_order_sud      (nbobuffermax)      )
        allocate(drifter_send_order_est      (nbobuffermax)      )
        allocate(drifter_send_order_ouest    (nbobuffermax)      )
        allocate(drifter_send_order_out      (nbobuffermax)      )

       endif             !>>> 

       if(case_==-1) then !-1-1-> ! desallouer la buffer zone !16-11-20
        deallocate(drifter_send_nord       )
        deallocate(drifter_recv_nord       )
        deallocate(drifter_send_sud        )
        deallocate(drifter_recv_sud        )
        deallocate(drifter_send_est        )
        deallocate(drifter_recv_est        )
        deallocate(drifter_send_ouest      )
        deallocate(drifter_recv_ouest      )

        deallocate(drifter_send_order_nord )
        deallocate(drifter_send_order_sud  )
        deallocate(drifter_send_order_est  )
        deallocate(drifter_send_order_ouest)
        deallocate(drifter_send_order_out  )
       endif              !-1-1-> ! desallouer la buffer zone


       if(case_==0.or. &       !    etat initial !16-11-20
          case_==2) then !ooo> ! ou reallocation en cours de run des coins la buffer zone 

        allocate(drifter_send_nordest        (nbobuffermax_c*drifter_dim2)    )
        allocate(drifter_send_nordouest      (nbobuffermax_c*drifter_dim2)    )
        allocate(drifter_send_sudest         (nbobuffermax_c*drifter_dim2)    )
        allocate(drifter_send_sudouest       (nbobuffermax_c*drifter_dim2)    )
        allocate(drifter_recv_nordest        (nbobuffermax_c*drifter_dim2)    )
        allocate(drifter_recv_nordouest      (nbobuffermax_c*drifter_dim2)    )
        allocate(drifter_recv_sudest         (nbobuffermax_c*drifter_dim2)    )
        allocate(drifter_recv_sudouest       (nbobuffermax_c*drifter_dim2)    )

        allocate(drifter_send_order_nordest  (nbobuffermax_c) )
        allocate(drifter_send_order_nordouest(nbobuffermax_c) )
        allocate(drifter_send_order_sudest   (nbobuffermax_c) )
        allocate(drifter_send_order_sudouest (nbobuffermax_c) )

       endif             !ooo> 
       if(case_==-2) then !-2-2-> ! desallouer les coins de la buffer zone !16-11-20

        deallocate(drifter_send_nordest        )
        deallocate(drifter_send_nordouest      )
        deallocate(drifter_send_sudest         )
        deallocate(drifter_send_sudouest       )
        deallocate(drifter_recv_nordest        )
        deallocate(drifter_recv_nordouest      )
        deallocate(drifter_recv_sudest         )
        deallocate(drifter_recv_sudouest       )

        deallocate(drifter_send_order_nordest  )
        deallocate(drifter_send_order_nordouest)
        deallocate(drifter_send_order_sudest   )
        deallocate(drifter_send_order_sudouest )

       endif              !-2-2-> ! desallouer les coins de la buffer zone !16-11-20

      end subroutine drifter_allocate

!...............................................

!......................................................................

      subroutine drifter_webcanals !04-09-22
      implicit none
      integer :: loop_,icnl_,point_,warning_=0,nexchg_=0 & !,nexchgmax_webc  &
       ,canalrank_send_,canalrank_recv_,loop_gridtype,gridtype_s,gridtype_r &
       ,tag_,new_kbumax_,rk_

! NB: la convention des drifters est differentes de celles des traceurs, a savoir que le point qui envoie est le point le plus extreme (donc sender=2)
      sender=2  
      receiver=1 
      gridtype1=1
      gridtype2=2
      warning_=0

!.....................
! RESET, ETAT INITIAL:
      if(.not.allocated(tabreq)) then !-INITIAL->
! Combien de connections multi-rank?:
      nexchg_=0
      do loop_gridtype=1,2
       if(loop_gridtype==1) then !>>>
        gridtype_s=gridtype1
        gridtype_r=gridtype2
       else                      !>>>
        gridtype_s=gridtype2
        gridtype_r=gridtype1
       endif                     !>>>
       do icnl_=1,nbcanal  
       do point_=1,2
        canalrank_recv_=canalrank(icnl_,point_,gridtype_r,receiver)
        canalrank_send_=canalrank(icnl_,point_,gridtype_s,sender)
! Pour etablir une connection mpi: 1 il faut que le rank des 2 connecteurs consideres soient positifs 
!                                  2 il faut que le rank des 2 connecteurs soient differents
!                                  3 il faut que l'un des 2 connecteurs soit par%rank
        if(canalrank_recv_>=0.and.canalrank_send_>=0) then                !test de connection 1>
         if(canalrank_recv_/=canalrank_send_) then                        !test de connection 2>
          if(canalrank_recv_==par%rank.or.canalrank_send_==par%rank) then !test de connection 3>
                   nexchg_=nexchg_+1
          endif                                                           !test de connection 3>
         endif                                                            !test de connection 2>
        endif                                                             !test de connection 1>
       enddo ! point_
       enddo ! icnl_
       enddo ! loop_gridtype
       nexchgmax_webc=nexchg_
       allocate(tabreq(nexchgmax_webc)) ; tabreq=0 ! chaque rank echange avec potententiellement tous les autres ranks sauf avec lui meme
       allocate(tstatus(mpi_status_size,nexchgmax_webc)) ; tstatus=-999
       allocate(nd_send_canal(0:nbdom-1)) ; nd_send_canal=-999 ! chaque rank echange avec potententiellement tous les autres ranks sauf avec lui meme
       allocate(nd_recv_canal(0:nbdom-1)) ; nd_recv_canal=-999 ! chaque rank echange avec potententiellement tous les autres ranks sauf avec lui meme
      endif                           !-INITIAL->

 73   continue ! Point de renvoi vers un nouveau depart suite A une augmentation de la taille des tableaux buffer
      if(warning_==2) then !-memoire-insuffisante-detectee->
       deallocate(drifter_send_order_canal)
       deallocate(drifter_send_canal)
       deallocate(drifter_recv_canal)
      endif                !-memoire-insuffisante-detectee->
      if(.not.allocated(drifter_send_order_canal)) then !RESET>
       allocate(drifter_send_order_canal(nbobuffercanalmax,0:nbdom-1)) ; drifter_send_order_canal=0
       allocate(drifter_send_canal(drifter_dim2*nbobuffercanalmax,0:nbdom-1)) ; drifter_send_canal=0
       allocate(drifter_recv_canal(drifter_dim2*nbobuffercanalmax,0:nbdom-1)) ; drifter_send_canal=0
       warning_=0
      endif                                             !RESET>

!.....................

! les echanges se font point par point de connexion (et non pas tous les
! points de connexion A la fois) car l'echange de tous les points A la fois est complique
! par le fait que envoyeurs et recepteurs ne sont pas symetriques (je
! veux dire que par exemple le ! rank1 peut envoyer au rank 10 seulement, alors que le rank 10 peut
! recevoir de rank1 et de rank5 etc...)

! Le sender est le canal: gridtype=gridtype_s
! le receiver est la mer: gridtype=gridtype_r


! Dans le cas des drifters, on a une etape supplementaire qui consiste A
! chercher combien de drifters sont A exporter et A envoyer ce nombre au
! recepteur afin que celui-ci sache combien d'element il doit recevoir

      nd_send_canal=0        ! tableau(par rank destinataire)
      nd_recv_canal=0
      nexchg_   =0

      do loop_gridtype=1,2

      if(loop_gridtype==1) then !>>>
       gridtype_s=gridtype1
       gridtype_r=gridtype2
      else                      !>>>
       gridtype_s=gridtype2
       gridtype_r=gridtype1
      endif                     !>>>

      do icnl_=1,nbcanal  

      do point_=1,2

      canalrank_recv_=canalrank(icnl_,point_,gridtype_r,receiver)
      canalrank_send_=canalrank(icnl_,point_,gridtype_s,sender)

      if(canalrank_recv_==par%rank.or.canalrank_send_==par%rank) then !test de connection 3>
      if(canalrank_recv_>=0.and.canalrank_send_>=0)              then !test de connection 1>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148
        isend=i_canalcoord(icnl_,point_,gridtype_s,sender)   ! -par%timax(1)
        jsend=j_canalcoord(icnl_,point_,gridtype_s,sender)   ! -par%tjmax(1)
        irecv=i_canalcoord(icnl_,point_,gridtype_r,receiver) ! -par%timax(1)
        jrecv=j_canalcoord(icnl_,point_,gridtype_r,receiver) ! -par%tjmax(1)


        do kbu=1,kbumax 

! test de presence d'un drifter sur un point de connexion:
        if(nint(drifter_l(kbu,1))==isend.and. &
           nint(drifter_l(kbu,2))==jsend) then !m°v°m>

         if(canalrank_recv_/=canalrank_send_) then !test de connection 2>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI

           nd_send_canal(canalrank_recv_) &
          =nd_send_canal(canalrank_recv_)+1


         endif                                     !test de connection 2> 

        endif !m°v°m> test de presence d'un drifter sur un point de connexion
        enddo ! boucle sur kbu

! Regle globale multi-rank de numerotation de tag pour un envoi de RK1 a RK2: tag=RK2+nbdom*RK1
      tag_=canalrank_recv_+canalrank_send_*nbdom

      if(canalrank_recv_/=canalrank_send_) then !--mpi-->
      if(par%rank==canalrank_send_)  then         !-send->

            nexchg_=nexchg_+1
            call mpi_issend(nd_send_canal(canalrank_recv_)          &
                           ,1                                       &
                           ,mpi_integer                             &
                           ,canalrank_recv_                         &
                           ,tag_                                    &
                           ,par%comm2d                              &
                           ,tabreq(nexchg_)                         &
                           ,ierr)

      else  if (par%rank==canalrank_recv_) then !-receive->

            nexchg_=nexchg_+1
            call mpi_irecv(nd_recv_canal(canalrank_send_)            &
                          ,1                                         &
                          ,mpi_integer                               &
                          ,canalrank_send_                           &
                          ,tag_                                      &
                          ,par%comm2d                                &
                          ,tabreq(nexchg_)                           &
                          ,ierr)

      endif                                       !-receive->
      endif                                     !--mpi-->

      endif                                                           !test de connection 1>
      endif                                                           !test de connection 3>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  
      enddo ! loop_gridtype=1,2

      if(nexchg_   >nexchgmax_webc   ) then !>>>
       write(10+par%rank,*)'nexchg_,nexchgmax_webc',nexchg_,nexchgmax_webc 
       stop 'drifter_webcanals nexchg_   >nexchgmax_webc   '
      endif                             !>>>

      ! La salle d'attente
      call mpi_waitall(nexchg_  &
             ,tabreq(1:nexchg_) &
          ,tstatus(:,1:nexchg_) &
          ,ierr)


!-!-!-!-!-!-!-!-!-!-!-!>
! VERIFICATION DE LA TAILLE MEMOIRE DES TABLEAUX ET REALLOCATION SI NECESSAIRE
! Il y a 2 choses A verifier: 
! 1- la dimension requise pour les tableaux "send" 
! 2- la dimension requise pour les tableaux "recv" (on la connait suite A l'echange qui vient d'avoir lieu)
      if(nbobuffercanalmax<max(maxval(nd_send_canal),maxval(nd_recv_canal))) then !m°v°m>
! Memoire insuffisante detectee pour la buffer zone. Reallouer les tableaux avec la dimension necessaire
        nbobuffercanalmax=max(maxval(nd_send_canal),maxval(nd_recv_canal))
        warning_=2
      endif                   !m°v°m>
! Partager le warning avec tous les ranks puis aviser
      call mpi_allreduce(warning_,k0,1,mpi_integer,mpi_max,par%comm2d,ierr)
      warning_=k0
      if(warning_/=0)goto 73 ! REFAIRE LE CALCUL AVEC LES NOUVELLES AVEC NOUVELLES DIMENSIONS
!-!-!-!-!-!-!-!-!-!-!-!>


!.....................................................................
! Translater la position i,j des drifters migrant
! Numeroter les drifters migrant avec drifter_send_order_canal
! Marquer les drifters migrant en changeant le signe de drifter_l(:,4)
      nd_send_canal=0

      do loop_gridtype=1,2
      if(loop_gridtype==1) then !>>>
       gridtype_s=gridtype1
       gridtype_r=gridtype2
      else                      !>>>
       gridtype_s=gridtype2
       gridtype_r=gridtype1
      endif                     !>>>

      do icnl_=1,nbcanal  
      do point_=1,2

      canalrank_recv_=canalrank(icnl_,point_,gridtype_r,receiver)
      canalrank_send_=canalrank(icnl_,point_,gridtype_s,sender)

      if(canalrank_recv_==par%rank.or.canalrank_send_==par%rank) then !test de connection 3>
      if(canalrank_recv_>=0.and.canalrank_send_>=0)              then !test de connection 1>

! https://docs.google.com/presentation/d/1O5jjBXzC1JnDyF5TkFTlLLf-C27vM7lmsWmqpz1zurg/edit#slide=id.g568c99ce7e_0_148
        isend=i_canalcoord(icnl_,point_,gridtype_s,sender)   ! -par%timax(1)
        jsend=j_canalcoord(icnl_,point_,gridtype_s,sender)   ! -par%tjmax(1)
        irecv=i_canalcoord(icnl_,point_,gridtype_r,receiver) ! -par%timax(1)
        jrecv=j_canalcoord(icnl_,point_,gridtype_r,receiver) ! -par%tjmax(1)

        do kbu=1,kbumax 

! test que le drifter n'a pas ete deja shiftE (necessaire en raison de la precision des calculs real4 et de l'arrondi 
! nint qui peuvent (rarement mais bon...) remettre la particule sur isend,jsend de loop_gridtype=2
        if(drifter_l(kbu,4)>0) then !pas deja shifte>

! test de presence d'un drifter sur un point de connexion:
        if(nint(drifter_l(kbu,1))==isend.and. &
           nint(drifter_l(kbu,2))==jsend) then !m°v°m>

! On translate la particule avant l'echange (il me semble plus simple de le faire avant qu'apres)
! d'autant qu'il peut ne pas y avoir d'echange si canalrank_recv_==canalrank_send_
          drifter_l(kbu,1)=drifter_l(kbu,1)+real(irecv-isend)
          drifter_l(kbu,2)=drifter_l(kbu,2)+real(jrecv-jsend)

         if(canalrank_recv_/=canalrank_send_) then !test de connection 2>
! sender et receiver ont un RANK DIFFERENT = ECHANGES MPI


           nd_send_canal(canalrank_recv_) &
          =nd_send_canal(canalrank_recv_)+1

         drifter_send_order_canal(nd_send_canal(canalrank_recv_),canalrank_recv_)=kbu
         drifter_l(kbu,4)=-abs(drifter_l(kbu,4))

         endif                                     !test de connection 2> 

        endif !m°v°m> test de presence d'un drifter sur un point de connexion
        endif                        !pas deja shifte>

        enddo ! boucle sur kbu

      endif                                                           !test de connection 1>
      endif                                                           !test de connection 3>
      enddo ! point_=1,2
      enddo ! icnl_=1,nbcanal  
      enddo ! loop_gridtype=1,2
!.....................................................................


! CHARGER LE TABLEAU BUFFER D'ENVOI drifter_send_canal avec drifter_l
      do k2=0,nbdom-1
      do k1=1,nd_send_canal(k2)
       kbu=drifter_send_order_canal(k1,k2)
       do loop_=1,drifter_dim2 ! loop_
        k3=loop_+(k1-1)*drifter_dim2
        drifter_send_canal(k3,k2)=drifter_l(kbu,loop_)
       enddo
      enddo
      enddo


      nexchg_=0
      do loop_gridtype=1,2
       if(loop_gridtype==1) then !>>>
        gridtype_s=gridtype1
        gridtype_r=gridtype2
       else                      !>>>
        gridtype_s=gridtype2
        gridtype_r=gridtype1
       endif                     !>>>
       do icnl_=1,nbcanal  
       do point_=1,2
        canalrank_recv_=canalrank(icnl_,point_,gridtype_r,receiver)
        canalrank_send_=canalrank(icnl_,point_,gridtype_s,sender)

        if(canalrank_recv_>=0.and.canalrank_send_>=0) then                !test de connection 1>
         if(canalrank_recv_/=canalrank_send_) then                        !test de connection 2>

! Regle globale multi-rank de numerotation de tag pour un envoi de RK1 a RK2: tag=RK2+nbdom*RK1
      tag_=canalrank_recv_+canalrank_send_*nbdom

      if(par%rank==canalrank_send_)  then         !-send->

            nexchg_=nexchg_+1
            k=max(drifter_dim2*nd_send_canal(canalrank_recv_),1)
            call mpi_issend(drifter_send_canal(1,canalrank_recv_)   &
                           ,k                                       &
                           ,mpi_real                                &
                           ,canalrank_recv_                         &
                           ,tag_                                    &
                           ,par%comm2d                              &
                           ,tabreq(nexchg_)                         &
                           ,ierr)

      else  if (par%rank==canalrank_recv_) then !-receive->

            nexchg_=nexchg_+1
            k=max(drifter_dim2*nd_recv_canal(canalrank_send_),1)
            call mpi_irecv(drifter_recv_canal(1,canalrank_send_)   &
                          ,k                                       &
                          ,mpi_real                                &
                          ,canalrank_send_                         &
                          ,tag_                                    &
                          ,par%comm2d                              &
                          ,tabreq(nexchg_)                         &
                          ,ierr)


      endif                                       !-receive->

         endif                                                            !test de connection 2>
        endif                                                             !test de connection 1>
       enddo ! point_
       enddo ! icnl_
       enddo ! loop_gridtype


      if(nexchg_   >nexchgmax_webc   ) &
      stop 'drifter_webcanals nexchg_   >nexchgmax_webc   '

      ! La salle d'attente
      call mpi_waitall(nexchg_  &
             ,tabreq(1:nexchg_) &
          ,tstatus(:,1:nexchg_) &
          ,ierr)

      end subroutine drifter_webcanals

!......................................................................

      subroutine drifter_update_numbering !04-09-22
      implicit none
      integer loop_

      flag_stop=0

! RENUMEROTER LES DRIFTER
      new_kbumax=kbumax                                                &
        +nd_recv_est+nd_recv_ouest+nd_recv_nord+nd_recv_sud            &
        +nd_recv_nordest+nd_recv_nordouest                             &
        +nd_recv_sudest +nd_recv_sudouest                              &
      -(nd_send_est+nd_send_ouest+nd_send_nord+nd_send_sud+nd_send_out &
       +nd_send_nordest+nd_send_nordouest                              &
       +nd_send_sudest +nd_send_sudouest  )

      if(nbcanal>0) then !pmx>
       new_kbumax= &
       new_kbumax+sum(nd_recv_canal)-sum(nd_send_canal)
      endif              !pmx>
      if(kbumax>nbomax)call drifter_error_message(3)

! Partie 1: enlever les "sortants" de la numerotation

! L'algo facile mais consommateur de calculs
!     do k=1,kbumax
!      if(drifter_l(k,4)>0) then !-n'est pas "sortant"->
!       kbu=kbu+1 ! incrementer le compteur seulement si le drifter n'est pas "sortant"
!       drifter_l(kbu,1:drifter_dim2)=drifter_l(k,1:drifter_dim2)
!      endif                     !-n'est pas "sortant"->
!     enddo

      if(new_kbumax>0) then !ooo>

       k2=0
       k1=kbumax
       do kbu=1,kbumax
        if(drifter_l(kbu,4)<0) then !-"sortant"->
         k2=k2+1 ! k2 compte les "sortants"
! on cherche un drifter au bout de la liste pour prendre la place du "sortant", mais ce drifter ne doit pas etre lui meme un "sortant", d'ou la recherche selective a suivre avec do while:
         do while (drifter_l(k1,4)<0)
          k1=k1-1
         enddo
         if(k1>0) then !ooo> !20-04-23
          drifter_l(kbu,1:drifter_dim2)=drifter_l(k1,1:drifter_dim2)  
         endif         !ooo> !20-04-23
         k1=k1-1
        endif                       !-"sortant"->
       enddo ! kbu
! Reset de la numerotation de la suite (partie 2):
       kbu=kbumax-k2  ! (k2=nombre de "sortants")

      else                  !ooo>

! Reset de la numerotation de la suite (partie 2):
       kbu=0

      endif                 !ooo>

! Partie 2: ajouter les "entrants" a la suite
      if(nbcanal>0) then !canaux>
! entrant par les canaux:
       do k2=0,nbdom-1 ! canaux situes dans le rank=k2
       do k1=1,nd_recv_canal(k2)
        kbu=kbu+1

        do loop_=1,drifter_dim2 ! loop_
         drifter_l(kbu,loop_)=drifter_recv_canal(loop_+(k1-1)*drifter_dim2,k2)
        enddo
       enddo
       enddo
      endif              !canaux>
! entrant par le nord:
      do k1=1,nd_recv_nord
       kbu=kbu+1
       do loop_=1,drifter_dim2 ! loop_
        drifter_l(kbu,loop_)=drifter_recv_nord(loop_+(k1-1)*drifter_dim2)
       enddo
      enddo
! entrant par le sud:
      do k1=1,nd_recv_sud
       kbu=kbu+1
       do loop_=1,drifter_dim2 ! loop_
        drifter_l(kbu,loop_)=drifter_recv_sud(loop_+(k1-1)*drifter_dim2)
       enddo
      enddo
! entrant par l'est:
      do k1=1,nd_recv_est
       kbu=kbu+1
       do loop_=1,drifter_dim2 ! loop_
        drifter_l(kbu,loop_)=drifter_recv_est(loop_+(k1-1)*drifter_dim2)
       enddo
      enddo
! entrant par l'ouest:
      do k1=1,nd_recv_ouest
       kbu=kbu+1
       do loop_=1,drifter_dim2 ! loop_
        drifter_l(kbu,loop_)=drifter_recv_ouest(loop_+(k1-1)*drifter_dim2)
       enddo
      enddo
! entrant par le nordest
      do k1=1,nd_recv_nordest
       kbu=kbu+1
       do loop_=1,drifter_dim2 ! loop_
        drifter_l(kbu,loop_)=drifter_recv_nordest(loop_+(k1-1)*drifter_dim2)
       enddo
      enddo
! entrant par le nordouest
      do k1=1,nd_recv_nordouest
       kbu=kbu+1
       do loop_=1,drifter_dim2 ! loop_
        drifter_l(kbu,loop_)=drifter_recv_nordouest(loop_+(k1-1)*drifter_dim2)
       enddo
      enddo
! entrant par le sudest
      do k1=1,nd_recv_sudest
       kbu=kbu+1
       do loop_=1,drifter_dim2 ! loop_
        drifter_l(kbu,loop_)=drifter_recv_sudest(loop_+(k1-1)*drifter_dim2)
       enddo
      enddo
! entrant par le sudouest
      do k1=1,nd_recv_sudouest
       kbu=kbu+1
       do loop_=1,drifter_dim2 ! loop_
        drifter_l(kbu,loop_)=drifter_recv_sudouest(loop_+(k1-1)*drifter_dim2)
       enddo
      enddo

      if(kbu/=new_kbumax) then
        write(10+par%rank,*)'BIZARRE kbu/=new_kbumax'
        write(10+par%rank,*)'kbu               ',kbu
        write(10+par%rank,*)'new_kbumax        ',new_kbumax
        write(10+par%rank,*)'kbumax            ',kbumax
        write(10+par%rank,*)'iteration3d       ',iteration3d
        write(10+par%rank,*)'sum(nd_recv_canal)',sum(nd_recv_canal)
        write(10+par%rank,*)'sum(nd_send_canal)',sum(nd_send_canal)
        flag_stop=1
        STOP 'BIZARRE kbu/=new_kbumax'
      endif

! UPDATE KBUMAX
      kbumax=new_kbumax

! RETABLIR LE SIGNE POSITIF DU NUMERO DES PARTICULES:
      do kbu=1,kbumax
       drifter_l(kbu,4)=abs(drifter_l(kbu,4))
      enddo

      call mpi_allreduce(flag_stop,k0,1,mpi_integer,mpi_sum,par%comm2d,ierr) 
      if(k0/=0)stop 'ERR drifter_update_numbering BIZARRE kbu/=new_kbumax see fort.xxx files'
!'

      end subroutine drifter_update_numbering

!...............................................

      subroutine drifter_tester_rk !22-07-22
      use module_principal ; use module_parallele
      implicit none
      double precision :: dti_ref_=180. , dist_ref_=50.

! Donner un champ de vitesse
      vel_u=0.
      vel_v=0.
      omega_w=0.
      i0=iglb/2
      j0=jglb/2
! definir un champ de vitesse academique au point "t"
      do j=0,jmax+1
      do i=0,imax+1
       i1=i+par%timax(1)
       j1=j+par%tjmax(1)
       x1=sqrt(real(i1-i0)**2+real(j1-j0)**2)
       x2=real(i1-i0)
       x3=real(j1-j0)
       x0=atan2(x3,x2)

!      x5=(x1/50.)*dxb/dti_fw
       x5=(x1/dist_ref_)*dxb/dti_ref_

       xy_t(i,j,1)=-x5*sin(x0)
       xy_t(i,j,2)=+x5*cos(x0)
      enddo
      enddo
! passer sur points u et v
      do j=0,jmax+1
      do i=1,imax+1
       vel_u(i,j,:,:)=0.5*(xy_t(i-1,j,1)+xy_t(i,j,1)) 
      enddo
      enddo
      do j=1,jmax+1
      do i=0,imax+1
       vel_v(i,j,:,:)=0.5*(xy_t(i,j-1,2)+xy_t(i,j,2)) 
      enddo
      enddo

      call graph_out

! Nombre d'iteration pour faire un tour complet = ! 2.*pi*dti_ref_*dist_ref_/dti_fw
      write(6,*)'tour complet en ' &
        ,2.*pi*dti_ref_*dist_ref_/dti_fw,' kounts'


      do iteration3d=0,iteration3d_max
      elapsedtime_now=elapsedtime_now+dti_fw
       call drifter_update
      enddo
!     call graph_out
      stop 'COCO4'

      end subroutine drifter_tester_rk

!...............................................

      subroutine drifter_turbulence_w !24-09-22
      use module_principal ; use module_parallele
      implicit none
!     real*4  :: random_variability=0.
      integer :: loop_ &
                ,loopmax_=2
!               ,loopmax_=1

      if(drifter_random_w==1) then !-tke-scheme->
! Echelle de vitesse verticale indexee sur TKE:
       do j=1,jmax ; do i=1,imax
        do k=kmin_w(i,j)+1,kmax
! Pour un meilleur fonctionnement elle est bornee par le nombre de courant
! (si elle le depasse beaucoup, les particules ont du mal A s'approcher des bords des couches
! limites et ne developpent pas le nuage homogene attendu)
           anyvar3d(i,j,k)=min(sqrt(2*tken_w(i,j,k)/3.),(depth_t(i,j,k)-depth_t(i,j,k-1))/dt_drf)
        enddo
       enddo       ; enddo
      endif                        !-tke-scheme->
      if(drifter_random_w==2) then !-KH-scheme->
! Echelle de vitesse verticale indexee sur coef de melange turbulent
       do j=1,jmax ; do i=1,imax
        do k=kmin_w(i,j)+1,kmax
! Pour un meilleur fonctionnement elle est bornee par le nombre de courant
           anyvar3d(i,j,k)=min(sqrt(kh_w(i,j,k)/dt_drf),(depth_t(i,j,k)-depth_t(i,j,k-1))/dt_drf)
        enddo
       enddo       ; enddo
      endif                        !-KH-scheme->

!     write(6,*)'bidouille patrick kh'
!      do j=1,jmax ; do i=1,imax
!       do k=kmin_w(i,j)+1,kmax
!          anyvar3d(i,j,k)=min(0.01*max(min(-depth_w(i,j,k),depth_w(i,j,k)+10.),0.), &
!          anyvar3d(i,j,k)=min(0.01*max(min(depth_w(i,j,k)-depth_w(i,j,1),10.-depth_w(i,j,k)+depth_w(i,j,1)),0.), &
!                              (depth_t(i,j,k)-depth_t(i,j,k-1))/dt_drf)
!          anyvar3d(i,j,k)=max(min(depth_w(i,j,k)-depth_w(i,j,1),10.-depth_w(i,j,k)+depth_w(i,j,1)),0.)
!          if(i==imax/2.and.j==jmax/2)write(10+par%rank,*)k,depth_w(i,j,k),anyvar3d(i,j,k)
!       enddo
!      enddo       ; enddo

! Les choix pour la tubulence dans la couche VQS n'etant pas definitifs
! et ne garantissant par consequent pas que les champs soient toujours
! definis sous kmin_w, on utilise le calcul (prudent) suivant:
      do j=1,jmax ; do i=1,imax
       do k=2,kmin_w(i,j)
        anyvar3d(i,j,k)=0.
       enddo
      enddo       ; enddo

! c.l. en k=1 et k=kmax+1
      do j=1,jmax ; do i=1,imax !30-10-20
       anyvar3d(i,j,1)=0.  ; anyvar3d(i,j,kmax+1)=0.
      enddo       ; enddo

! continuite mpi de anyvar3d
      call obc_int_mpi_anyvar3d('z1')

      end subroutine drifter_turbulence_w !24-09-22

!...............................................

      subroutine drifter_netcdf_file(init_or_iter_) !24-10-22
      implicit none
      integer init_or_iter_,varid_,loop_,kbumax_
      integer(kind=MPI_OFFSET_KIND) :: kbumax_mpikind,len_,start(2),edge(2)
      double precision :: valr8(1)

      kbumax_glb1(par%rank)=kbumax
      call mpi_allreduce(kbumax_glb1,kbumax_glb2,nbdom,mpi_integer,mpi_sum,par%comm2d,ierr)

! Commenter cette ligne si on veut que kbumax_glb garde toujours la valeur obtenue A l'initialisation des drifters
!     kbumax_glb=sum(kbumax_glb2)

      kbumax_mpikind=kbumax_glb

      if(init_or_iter_==0) then !--->
       texte250='tmp/drifter_initial.nc'
      else                      !--->
       write(texte30,'(i4,5(i2))')year_now,month_now,day_now,hour_now,minute_now,second_now
       if(texte30(5:5)  ==' ')texte30(5:5)='0'
       if(texte30(7:7)  ==' ')texte30(7:7)='0'
       if(texte30(9:9)  ==' ')texte30(9:9)='0'
       if(texte30(11:11)==' ')texte30(11:11)='0'
       if(texte30(13:13)==' ')texte30(13:13)='0'
       texte250='tmp/drifter_'//trim(texte30)//'.nc'
      endif                     !--->

       status=nfmpi_create(par%comm2d,trim(texte250),nf_clobber+ NF_64BIT_OFFSET, MPI_INFO_NULL,ncid)
       if(status/=0)stop 'err drifters nfmpi_create'

       status=nfmpi_def_dim(ncid,'ntraj',kbumax_mpikind,i_t_dim)
       if(status/=0)stop 'Err drifter nfmpi_def_dim kbumax'

       status=nfmpi_def_dim(ncid,'time',nfmpi_unlimited,time_dim)
       if(status/=0)stop 'Err drifter nfmpi_def_dim time'

! Les parametres "initiaux" qui ne changent pas au cours du temps
      if(init_or_iter_==0) then !-initial->

       vardim(1)=i_t_dim ; vardim(2)=time_dim

       status=nfmpi_def_var(ncid,'init_i_index',nf_real,1,vardim,varid_) !28-02-23
       texte30='initial_i_index' ; len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)
       if(status/=0)stop 'Err drifter nfmpi_put_att_text long_name i'

       status=nfmpi_def_var(ncid,'init_j_index',nf_real,1,vardim,varid_) !28-02-23
       texte30='initial_j_index' ; len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)
       if(status/=0)stop 'Err drifter nfmpi_put_att_text long_name j'

       status=nfmpi_def_var(ncid,'init_k_index',nf_real,1,vardim,varid_) !28-02-23
       texte30='initial_k_index' ; len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)
       if(status/=0)stop 'Err drifter nfmpi_put_att_text long_name k'

       status=nfmpi_def_var(ncid,'init_lon',nf_real,1,vardim,varid_)
       texte30='initial longitude' ; len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)
       if(status/=0)stop 'Err drifter nfmpi_put_att_text long_name lon'

       status=nfmpi_def_var(ncid,'init_lat',nf_real,1,vardim,varid_)
       texte30='initial latitude' ; len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)
       if(status/=0)stop 'Err drifter nfmpi_put_att_text long_name lat'

       status=nfmpi_def_var(ncid,'init_z',nf_real,1,vardim,varid_)
       texte30='initial depth' ; len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)
       if(status/=0)stop 'Err drifter nfmpi_put_att_text long_name z'

       status=nfmpi_def_var(ncid,'init_identifier',nf_real,1,vardim,varid_)
       texte30='initial drifter ID' ; len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)
       if(status/=0) &
       stop 'Err drifter nfmpi_put_att_text long_name init_id'

       status=nfmpi_def_var(ncid,'init_time',nf_real,1,vardim,varid_)

       write(texte30,'(a30,i4,5(a1,i2))')       &
         'initial time in seconds since '       &
        ,datesim(1,1),'-',datesim(2,1),'-',datesim(3,1),' ' &
        ,datesim(4,1),':',datesim(5,1),':',datesim(6,1)
       if(texte30(36:36)==' ')texte30(36:36)='0'          
       if(texte30(39:39)==' ')texte30(39:39)='0'       
       if(texte30(42:42)==' ')texte30(42:42)='0'
       if(texte30(45:45)==' ')texte30(45:45)='0'
       if(texte30(48:48)==' ')texte30(48:48)='0'
       len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)

! Les parametres archivEs dans drifter_l(:,6:drifter_dim2)
       do loop_=6,drifter_dim2
        write(texte30,'(a,i0)')'init_param',loop_ ; len_=len_trim(texte30)
        status=nfmpi_def_var(ncid,trim(texte30),nf_real,1,vardim,varid_)
        status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)
       enddo

      endif                     !-initial->

! Les parametres qui evoluent au cours du temps:
      if(init_or_iter_/=0) then !-en-cours-de-simulation->

       vardim(1)=time_dim

       status=nfmpi_def_var(ncid,'time',nf_double,1,vardim,varid_)
       write(texte30,'(a30,i4,5(a1,i2))')       &
         'current time in seconds since '       &
        ,datesim(1,1),'-',datesim(2,1),'-',datesim(3,1),' ' &
        ,datesim(4,1),':',datesim(5,1),':',datesim(6,1)
       if(texte30(36:36)==' ')texte30(36:36)='0'          
       if(texte30(39:39)==' ')texte30(39:39)='0'       
       if(texte30(42:42)==' ')texte30(42:42)='0'
       if(texte30(45:45)==' ')texte30(45:45)='0'
       if(texte30(48:48)==' ')texte30(48:48)='0'
       len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)
       if(status/=0)stop 'Err drifter nfmpi_put_att_text long_name time'

       vardim(1)=i_t_dim ; vardim(2)=time_dim

       status=nfmpi_def_var(ncid,'i_index',nf_real,2,vardim,varid_)
       texte30='i_index' ; len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)

       status=nfmpi_def_var(ncid,'j_index',nf_real,2,vardim,varid_)
       texte30='j_index' ; len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)

       status=nfmpi_def_var(ncid,'k_index',nf_real,2,vardim,varid_)
       texte30='k_index' ; len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)

       status=nfmpi_def_var(ncid,'longitude',nf_real,2,vardim,varid_)
       texte30='longitude' ; len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)

       status=nfmpi_def_var(ncid,'latitude',nf_real,2,vardim,varid_)
       texte30='latitude' ; len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)

       status=nfmpi_def_var(ncid,'depth',nf_real,2,vardim,varid_)
       texte30='depth' ; len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)

       status=nfmpi_def_var(ncid,'identifier',nf_real,2,vardim,varid_)
       texte30='drifter ID' ; len_=len_trim(texte30)
       status=nfmpi_put_att_text(ncid,varid_,'long_name',len_,texte30)

      endif                     !-en-cours-de-simulation->

! fin de la partie "entete"
       status=nfmpi_enddef(ncid)

! Ecriture des parametres initiaux
      if(init_or_iter_==0) then !-initial->

       start(1)=1+sum(kbumax_glb2(0:par%rank-1))
       edge(1)=kbumax
       kbumax_=max(kbumax,1)
       if(edge(1)==0)start(1)=1

       status=nfmpi_inq_varid(ncid,'init_i_index',varid_) !28-02-23
       if(status/=0)stop 'drifter nfmpi_inq_varid init_i'
       status=nfmpi_put_vara_real_all(ncid,varid_,start(1:1),edge(1:1),drifter_l(1:kbumax_,1))
       if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all init_lon'

       status=nfmpi_inq_varid(ncid,'init_j_index',varid_) !28-02-23
       if(status/=0)stop 'drifter nfmpi_inq_varid init_j'
       status=nfmpi_put_vara_real_all(ncid,varid_,start(1:1),edge(1:1),drifter_l(1:kbumax_,2))
       if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all init_lon'

       status=nfmpi_inq_varid(ncid,'init_k_index',varid_) !28-02-23
       if(status/=0)stop 'drifter nfmpi_inq_varid init_k'
       status=nfmpi_put_vara_real_all(ncid,varid_,start(1:1),edge(1:1),drifter_l(1:kbumax_,3))
       if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all init_lon'

       call drifter_ij_to_lon
       status=nfmpi_inq_varid(ncid,'init_lon',varid_)
       if(status/=0)stop 'drifter nfmpi_inq_varid init_lon'
       status=nfmpi_put_vara_real_all(ncid,varid_,start(1:1),edge(1:1),drifter_l(1:kbumax_,0))
       if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all init_lon'

       call drifter_ij_to_lat
       status=nfmpi_inq_varid(ncid,'init_lat',varid_)
       if(status/=0)stop 'drifter nfmpi_inq_varid init_lat'
       status=nfmpi_put_vara_real_all(ncid,varid_,start(1:1),edge(1:1),drifter_l(1:kbumax_,0))
       if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all init_lat'

       call drifter_ijk_to_z
       status=nfmpi_inq_varid(ncid,'init_z',varid_)
       if(status/=0)stop 'drifter nfmpi_inq_varid init_z'
       status=nfmpi_put_vara_real_all(ncid,varid_,start(1:1),edge(1:1),drifter_l(1:kbumax_,0))
       if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all init_z'

       status=nfmpi_inq_varid(ncid,'init_identifier',varid_)
       if(status/=0)stop 'drifter nfmpi_inq_varid init_id'
       status=nfmpi_put_vara_real_all(ncid,varid_,start(1:1),edge(1:1),drifter_l(1:kbumax_,4))
       if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all init_id'

       status=nfmpi_inq_varid(ncid,'init_time',varid_)
       if(status/=0)stop 'drifter nfmpi_inq_varid init_time'
       status=nfmpi_put_vara_real_all(ncid,varid_,start(1:1),edge(1:1),drifter_l(1:kbumax_,5))
       if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all init_time'

       do loop_=6,drifter_dim2 !06-12-22
        write(texte30,'(a,i0)')'init_param',loop_ ; len_=len_trim(texte30)
        status=nfmpi_inq_varid(ncid,trim(texte30),varid_)
        if(status/=0)stop 'drifter nfmpi_inq_varid 2676'
        status=nfmpi_put_vara_real_all(ncid,varid_,start(1:1),edge(1:1),drifter_l(1:kbumax_,loop_))
        if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all 2678'
       enddo

      endif                     !-initial->

! Ecriture des parametres qui evoluent dans le temsp
      if(init_or_iter_/=0) then !-evoluent-en-cours-de-simulation->

       start(1)=1
        edge(1)=1
        if(par%rank==0)edge(1)=0

! Ecrire le temps
       valr8(1)=elapsedtime_now
       status=nfmpi_inq_varid(ncid,'time',varid_)
       status=nfmpi_put_vara_double_all(ncid,varid_,start(1:1),edge(1:1),valr8)

       start(1)=1+sum(kbumax_glb2(0:par%rank-1))
       edge(1)=kbumax
       kbumax_=max(kbumax,1)
       if(edge(1)==0)start(1)=1

       start(2)=1 ; edge(2)=1

       status=nfmpi_inq_varid(ncid,'i_index',varid_) !28-02-23
       if(status/=0)stop 'drifter nfmpi_inq_varid i'
       status=nfmpi_put_vara_real_all(ncid,varid_,start(1:2),edge(1:2),drifter_l(1:kbumax_,1))
       if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all i'

       status=nfmpi_inq_varid(ncid,'j_index',varid_) !28-02-23
       if(status/=0)stop 'drifter nfmpi_inq_varid j'
       status=nfmpi_put_vara_real_all(ncid,varid_,start(1:2),edge(1:2),drifter_l(1:kbumax_,2))
       if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all j'

       status=nfmpi_inq_varid(ncid,'k_index',varid_) !28-02-23
       if(status/=0)stop 'drifter nfmpi_inq_varid k'
       status=nfmpi_put_vara_real_all(ncid,varid_,start(1:2),edge(1:2),drifter_l(1:kbumax_,3))
       if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all j'

       call drifter_ij_to_lon
       status=nfmpi_inq_varid(ncid,'longitude',varid_)
       if(status/=0)stop 'drifter nfmpi_inq_varid lon'
       status=nfmpi_put_vara_real_all(ncid,varid_,start(1:2),edge(1:2),drifter_l(1:kbumax_,0))
       if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all lon'

       call drifter_ij_to_lat
       status=nfmpi_inq_varid(ncid,'latitude',varid_)
       if(status/=0)stop 'drifter nfmpi_inq_varid lat'
       status=nfmpi_put_vara_real_all(ncid,varid_,start(1:2),edge(1:2),drifter_l(1:kbumax_,0))
       if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all lat'

       call drifter_ijk_to_z
       status=nfmpi_inq_varid(ncid,'depth',varid_)
       if(status/=0)stop 'drifter nfmpi_inq_varid z'
       status=nfmpi_put_vara_real_all(ncid,varid_,start(1:2),edge(1:2),drifter_l(1:kbumax_,0))
       if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all z'

       status=nfmpi_inq_varid(ncid,'identifier',varid_)
       if(status/=0)stop 'drifter nfmpi_inq_varid id'
       status=nfmpi_put_vara_real_all(ncid,varid_,start(1:2),edge(1:2),drifter_l(1:kbumax_,4))
       if(status/=0)stop 'Err drifter nfmpi_put_vara_real_all id'

      endif                     !-evoluent-en-cours-de-simulation->
      
       status=nfmpi_close(ncid)
       if(status/=0)stop 'err drifters nfmpi_close'

      end subroutine drifter_netcdf_file

!...............................................

      subroutine drifter_ij_to_lon
      implicit none

      lb2=lbound(drifter_l)
      if(lb2(2)>0) stop 'drifter_ij_to_lon lb2(2)>0'
      
! Deduire la longitude des indices:
      do kbu=1,kbumax
       deci=drifter_l(kbu,1)-par%timax(1)
       decj=drifter_l(kbu,2)-par%tjmax(1)
       i=int(deci)
       j=int(decj)
       rapi=deci-real(i)
       rapj=decj-real(j)
       
       x1=0. ; x2=0. ; x3=0.
       if(lon_t(i  ,j+1)-lon_t(i  ,j  )> pi)x1=-2*pi
       if(lon_t(i  ,j+1)-lon_t(i  ,j  )<-pi)x1= 2*pi
       if(lon_t(i+1,j  )-lon_t(i  ,j  )> pi)x2=-2*pi
       if(lon_t(i+1,j  )-lon_t(i  ,j  )<-pi)x2= 2*pi
       if(lon_t(i+1,j+1)-lon_t(i  ,j  )> pi)x3=-2*pi
       if(lon_t(i+1,j+1)-lon_t(i  ,j  )<-pi)x3= 2*pi

       drifter_l(kbu,0)=rad2deg*( & !>>>
         (1.-rapi)*(1.-rapj)* lon_t(i  ,j  )       &
        +(1.-rapi)*    rapj *(lon_t(i  ,j+1)+x1)   &
        +    rapi *(1.-rapj)*(lon_t(i+1,j  )+x2)   &
        +    rapi *    rapj *(lon_t(i+1,j+1)+x3)   &
                                )   !>>>
      enddo ! kbu

      end subroutine drifter_ij_to_lon

!...............................................

      subroutine drifter_ij_to_lat
      implicit none

! Deduire la latitude des indices:
      do kbu=1,kbumax
       deci=drifter_l(kbu,1)-par%timax(1)
       decj=drifter_l(kbu,2)-par%tjmax(1)
       i=int(deci)
       j=int(decj)
       rapi=deci-real(i)
       rapj=decj-real(j)
       
       drifter_l(kbu,0)=rad2deg*( & !>>>
         (1.-rapi)*(1.-rapj)*lat_t(i  ,j  )       &
        +(1.-rapi)*    rapj *lat_t(i  ,j+1)       &
        +    rapi *(1.-rapj)*lat_t(i+1,j  )       &
        +    rapi *    rapj *lat_t(i+1,j+1)       &
                                )   !>>>
      enddo ! kbu

      end subroutine drifter_ij_to_lat

!...............................................

      subroutine drifter_ijk_to_z
      implicit none

! Deduire z des indices
      do kbu=1,kbumax

       deci=drifter_l(kbu,1)-par%timax(1)
       decj=drifter_l(kbu,2)-par%tjmax(1)
       deck=drifter_l(kbu,3)
       i=int(deci)
       j=int(decj)
       k=int(deck)
       rapi=deci-real(i)
       rapj=decj-real(j)
       rapk=deck-real(k)
       k1=max0(k  ,1)
       k2=min0(k+1,kmax)

       drifter_l(kbu,0)=                                     &
         (1.-rapi)*(1.-rapj)*(1.-rapk)*depth_t(i  ,j  ,k1)   &
        +(1.-rapi)*(1.-rapj)*    rapk *depth_t(i  ,j  ,k2)   &
        +(1.-rapi)*    rapj *(1.-rapk)*depth_t(i  ,j+1,k1)   &
        +(1.-rapi)*    rapj *    rapk *depth_t(i  ,j+1,k2)   &
        +    rapi *(1.-rapj)*(1.-rapk)*depth_t(i+1,j  ,k1)   &
        +    rapi *(1.-rapj)*    rapk *depth_t(i+1,j  ,k2)   &
        +    rapi *    rapj *(1.-rapk)*depth_t(i+1,j+1,k1)   &
        +    rapi *    rapj *    rapk *depth_t(i+1,j+1,k2)

      enddo

      end subroutine drifter_ijk_to_z

!...............................................
      end module module_drifter


! Une routine pour tester la normalisatio du nombre aleatoire
