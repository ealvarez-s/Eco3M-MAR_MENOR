      module module_drifter
      implicit none
      integer :: drifter_dim2=9     &  !15-05-15
                ,rungekutta_loop=1  &
                ,rungekutta_order=2
!     double precision :: deci1,decj1,deck1
      double precision :: drifter_l1,drifter_l2,drifter_l3
!______________________________________________________________________
! SYMPHONIE ocean model
! release 268 - last update: 23-11-19
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
!______________________________________________________________________

contains

!...............................................................................

      subroutine drifter_update
      use module_principal
      use module_parallele
      implicit none


      if(modulo(iteration3d,dt_drf_over_dti_fw)/=0) return !24-11-19
      dt_drf=dt_drf_over_dti_fw*dti_fw


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! MISE A JOUR DE LA POSITION DES BOUEES                                !22/03/02
! DEBUT:
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!.................................................
! Definition d'une echelle de vitesse verticale turbulente
      if(drifter_random_w==1) then !m°v°m> !22-11-19

         do k=2,kmax ; do j=1,jmax ; do i=1,imax
          anyvar3d(i,j,k)=sqrt(2*tken_w(i,j,k)/3.)
         enddo       ; enddo       ; enddo
! c.l. en k=1 et k=kmax+1
         do j=1,jmax ; do i=1,jmax
          anyvar3d(i,j,1)=0.  ; anyvar3d(i,j,kmax+1)=0.
         enddo       ; enddo

      else                         !m°v°m> !22-11-19

          anyvar3d=0. ! pas de vitesse turbulente aleatoire

      endif                        !m°v°m> !22-11-19
!..........

      do 19 kbu=1,kbumax

      if(drifter_random_w==1) then !m°v°m> !22-11-19
! x3_r4 nombre aleatoire muktipliant anyvar3d. 
         call random_number(x3_r4)
! La normalisation en suivant permet que l'energie cinetique deduite de
! la vitesse verticale aleatoire est bien tken_w/3 en moyenne. Voir la
! routine fortran test_random en fin du present fichier.
         x3_r4=2.*sqrt(3.)*(x3_r4-0.5)
      else                         !m°v°m> !22-11-19
         x3_r4=0.
      endif                        !m°v°m> !22-11-19

! BIDOUILLE PATRICK
!      if(nint(drifter_l(kbu,4))==1)rungekutta_order=1
!      if(nint(drifter_l(kbu,4))==2)rungekutta_order=2

      do 19 rungekutta_loop=1,rungekutta_order ! rungekutta_order=2

      if(elapsedtime_now>drifter_l(kbu,5)) then !%%%%%%%%%%%%%%%%%%%%%%> !03-03-15

!     write(10+par%Rank,*)kbu,flag_buo_w,id_wdrifter,drifter_l(kbu,id_wdrifter)

      if(rungekutta_loop==1) then !m°v°m> !18-11-19
! temps 1: vitesses interpolees au point de depart de la particule
       deci=drifter_l(kbu,1)-par%timax(1)
       decj=drifter_l(kbu,2)-par%tjmax(1)
       deck=drifter_l(kbu,3)
      endif                       !m°v°m> !18-11-19
      if(rungekutta_loop==2) then !m°v°m> !18-11-19
! temps 2: vitesses interpolees a mi-chemin entre depart et arrivee (approchee a la premiere iteration du RK2)
       deci=0.5*(drifter_l(kbu,1)+drifter_l1)-par%timax(1)
       decj=0.5*(drifter_l(kbu,2)+drifter_l2)-par%tjmax(1)
       deck=0.5*(drifter_l(kbu,3)+drifter_l3)
!         Attention qu'a l'etape 1 du RK2 la particule est peut etre sortie du
!         sousdomaine mpi (auquel cas c'est la position du depart qui est reprise)
          if(deci<=1.or.deci>=imax.or.decj<=1.or.decj>=jmax) then
            deci=drifter_l(kbu,1)-par%timax(1)
            decj=drifter_l(kbu,2)-par%tjmax(1)
            deck=drifter_l(kbu,3)
          endif
      endif                       !m°v°m> !18-11-19

      i=int(deci)
      j=int(decj)
      k=int(deck)
      rapi=deci-real(i)
      rapj=decj-real(j)
      rapk=deck-real(k)

      if(i<1.or.i>imax.or.j<1.or.j>jmax) &
      stop 'indices en dehors des limites dans drifter_update'

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

      x1=(1.-rapi)*(1.-rapj)*(1.-rapk)*(vel_u(i  ,j  ,k1,1)   &
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

! x2*dy = deuxieme composante de la vitesse sur le point decimal (DECI,DECJ,DECK)

      x2=(1.-rapi)*(1.-rapj)*(1.-rapk)*(vel_v(i  ,j  ,k1,1) &
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

      x3=( & !ooo>
         (1.-rapi)*(1.-rapj)*(1.-rapk)*(omega_w(i  ,j  ,k  ,1)       & !courant
                             +flag_buo_w*drifter_l(kbu,id_wdrifter)  & !vitesse de flottabilite de la particule
                                      +anyvar3d(i  ,j  ,k  )*x3_r4 ) & !vitesse turbulente aleatoire

        +(1.-rapi)*(1.-rapj)*    rapk *(omega_w(i  ,j  ,k+1,1)       &
                             +flag_buo_w*drifter_l(kbu,id_wdrifter)  &
                                      +anyvar3d(i  ,j  ,k+1)*x3_r4 ) &

        +(1.-rapi)*    rapj *(1.-rapk)*(omega_w(i  ,j+1,k  ,1)       &
                             +flag_buo_w*drifter_l(kbu,id_wdrifter)  &
                                      +anyvar3d(i  ,j+1,k  )*x3_r4 ) &

        +(1.-rapi)*    rapj *    rapk *(omega_w(i  ,j+1,k+1,1)       &
                             +flag_buo_w*drifter_l(kbu,id_wdrifter)  &
                                      +anyvar3d(i  ,j+1,k+1)*x3_r4 ) &

        +    rapi *(1.-rapj)*(1.-rapk)*(omega_w(i+1,j  ,k  ,1)       &
                             +flag_buo_w*drifter_l(kbu,id_wdrifter)  &
                                      +anyvar3d(i+1,j  ,k  )*x3_r4 ) &

        +    rapi *(1.-rapj)*    rapk *(omega_w(i+1,j  ,k+1,1)       &
                             +flag_buo_w*drifter_l(kbu,id_wdrifter)  &
                                      +anyvar3d(i+1,j  ,k+1)*x3_r4 ) &

        +    rapi *    rapj *(1.-rapk)*(omega_w(i+1,j+1,k  ,1)       &
                             +flag_buo_w*drifter_l(kbu,id_wdrifter)  &
                                      +anyvar3d(i+1,j+1,k  )*x3_r4 ) &

        +    rapi *    rapj *    rapk *(omega_w(i+1,j+1,k+1,1)       &
                             +flag_buo_w*drifter_l(kbu,id_wdrifter)  &
                                      +anyvar3d(i+1,j+1,k+1)*x3_r4 ) &
         ) & !ooo>
        /( & !ppp> !23-11-19
         (1.-rapi)*(1.-rapj)*(1.-rapk)*(depth_t(i  ,j  ,k  )         &
                                       -depth_t(i  ,j  ,k-1))        &

        +(1.-rapi)*(1.-rapj)*    rapk *(depth_t(i  ,j  ,k+1)         &
                                       -depth_t(i  ,j  ,k  ))        &

        +(1.-rapi)*    rapj *(1.-rapk)*(depth_t(i  ,j+1,k  )         &
                                       -depth_t(i  ,j+1,k-1))        &

        +(1.-rapi)*    rapj *    rapk *(depth_t(i  ,j+1,k+1)         &
                                       -depth_t(i  ,j+1,k  ))        &

        +    rapi *(1.-rapj)*(1.-rapk)*(depth_t(i+1,j  ,k  )         &
                                       -depth_t(i+1,j  ,k-1))        &

        +    rapi *(1.-rapj)*    rapk *(depth_t(i+1,j  ,k+1)         &
                                       -depth_t(i+1,j  ,k  ))        &

        +    rapi *    rapj *(1.-rapk)*(depth_t(i+1,j+1,k  )         &
                                       -depth_t(i+1,j+1,k-1))        &

        +    rapi *    rapj *    rapk *(depth_t(i+1,j+1,k+1)         &
                                       -depth_t(i+1,j+1,k  ))        &
         )   !ppp>

!.......................................................................!
! interpolation de la vitesse omega au point (DECI,DECJ,DECK)
! (attention au glissement d'indice sur k)
! Fin.
!.......................................................................!


      if(rungekutta_loop==2.or.rungekutta_order==1) then !m°v°m> !18-11-19
! temps 2 du RK2 ou temps 1 du RK1 (position d'arrivee finale)
       drifter_l(kbu,1)=        drifter_l(kbu,1)+x1*dt_drf
       drifter_l(kbu,2)=        drifter_l(kbu,2)+x2*dt_drf
!      drifter_l(kbu,3)=min(max(drifter_l(kbu,3)+x3*dt_drf,1.),real(kmax))
               drifter_l(kbu,3)=min(             & !ooo>
           max(drifter_l(kbu,3)+x3*dt_drf,real(min(kmin_w(i,j),kmin_w(i+1,j),kmin_w(i,j+1),kmin_w(i+1,j+1)))) & !23-11-19
                                    ,real(kmax))   !ooo>
!      drifter_l(kbu,3)=kmax ! bidouille patrick pour k=kmax
      else                                               !m°v°m> !18-11-19
! temps 1 du RK2 (position d'arrivee first guess)
       drifter_l1=        drifter_l(kbu,1)+x1*dt_drf
       drifter_l2=        drifter_l(kbu,2)+x2*dt_drf
!      drifter_l3=min(max(drifter_l(kbu,3)+x3*dt_drf,1.),real(kmax))
           drifter_l3=min(            & !ooo>
       max(drifter_l(kbu,3)+x3*dt_drf,real(min(kmin_w(i,j),kmin_w(i+1,j),kmin_w(i,j+1),kmin_w(i+1,j+1)))) &
                          ,real(kmax))  !ooo>
!      drifter_l3=kmax ! bidouille patrick pour k=kmax
      endif                                              !m°v°m> !18-11-19

!     if(rungekutta_loop==1.and.kbu==1) then
!           write(10+par%Rank,*)kbu,x3_r4,anyvar3d(i,j,k), 'tonton'
!           write(10+par%Rank,*)anyvar3d(i,j,k),i,j,k,tken_w(i,j,k),imax,jmax
!     endif

      endif                                   !%%%%%%%%%%%%%%%%%%%%%%> !03-03-15
   19 continue

#ifdef parallele
      call drifter_exchange_driver
#endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! MISE A JOUR DE LA POSITION DES BOUEES                                !22/03/02
! FIN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      end subroutine drifter_update

!................................................................................

      subroutine drifter_initial
      use module_principal
      use module_parallele
      implicit none
      double precision posx_,posy_,posz_ !15-05-15
      real vertical_buoyancy_velocity
!     integer loop_   ,kglbmax_   ,in_or_out_
      integer kbuglb_   ,kglbmax_   ,in_or_out_
      character*60 txt_spy_
#ifdef synopsis
       subroutinetitle='drifter_initial'
       subroutinedescription='Reads notebook_drifter'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

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

      if(par%rank==0)write(6,*)'lecture de ',nomfichier(16) !02-01-14
      open(unit=3,file=nomfichier(16)) ! lecture de notebook_drifter    !27/04/09
      read(3,*)drifter_onoff
  400 continue

      if(drifter_onoff==0)return

      call drifter_allocate !25-02-12

      read(3,'(a)')drifter_initial_mode
      read(3,*)drifter_out_sampling
      drifter_out_sampling=max(drifter_out_sampling,dti_fw) !24-11-19
      read(3,*)drifter_output_files    !24-11-19
      read(3,*)rungekutta_order        !19-11-19
      read(3,*)dt_drf_over_dti_fw      !24-11-19
      read(3,*)flag_buo_w,id_wdrifter  !22-11-19
      read(3,*)drifter_random_w        !22-11-19

      if(drifter_initial_mode(1:7)=='fortran') then !------->
           call drifter_initial_fortran
           goto 488
      endif                                         !------->

      read(3,*)                                                        !07/10/02

  369 continue

! Initial location of drifters:
!     read(3,*,end=363)posx_,posy_,posz_,k0
!     read(3,*,end=363)posx_,posy_,posz_,k0,x10 ! bidouille patrick oU x10 est vitesse de la particule
      read(3,*,end=363)posx_,posy_,posz_,k0  
      read(3,*)vertical_buoyancy_velocity         ! bidouille patrick oU vertical_buoyancy_velocity est vitesse de la particule
! Initial time of drifters:
      read(3,*)i1,i2,i3,i4,i5,i6 !An,Mois,Jour,Heure,Minute,Seconde !03-03-15
      call datetokount(i1,i2,i3,i4,i5,i6) !An,Mois,Jour,Heure,Minute,Seconde

      kbuglb_   =kbuglb_   +1

      if(k0==0)then !§§§§§§§§> ! position en lat. lon. prof(<0)

         latit1=posx_*pi/180.
         longi1=posy_*pi/180.
         call latlon_to_ij('loc')                                     !21-10-10

         if(deci>2.and.deci<real(imax).and.  & ! test: la drifter est elle dans le domaine?
            decj>2.and.decj<real(jmax)) then   !>>>>>>>>>>>>>>>>> !30-11-15

         in_or_out_   =1  ! 1=in
! Si oui, incrementer le compteur de drifter:
         kbu=kbu+1
         kbumax=kbumax+1
!        kbuglb_   =kbuglb_   +1
         if(kbu>nbomax)call drifter_error_message(1)

! Saisir la position en indices globaux dans tableau drifter_l:
         drifter_l(kbu,1)=deci+par%timax(1)
         drifter_l(kbu,2)=decj+par%tjmax(1)
         drifter_l(kbu,4)=real( kbuglb_    )
         drifter_l(kbu,5)=elapsedtime_out !03-03-15
         if(flag_buo_w==1)drifter_l(kbu,id_wdrifter)=vertical_buoyancy_velocity ! bidouille patrick vitesse de la particule 

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
               if(k==1.and.posz_<=z1)drifter_l(kbu,3)=un
            if(k==kmax-1.and.posz_>=z2)drifter_l(kbu,3)=kmax
         enddo

         endif                                 !>>>>>>>>>>>>>>>>>

      else           !§§§§§§§§> ! position i j k

! conversion des indices globaux en indices locaux
         deci=posx_-par%timax(1) ;  decj=posy_-par%tjmax(1)

! Tester si cette drifter est dans le proc:
         if(deci>2.and.deci<real(imax).and.  & ! test: la drifter est elle dans le domaine?
            decj>2.and.decj<real(jmax)) then   !\\\\\\\\\\\> !30-11-15

         in_or_out_   =1  ! 1=in
! Si oui, incrementer le compteur de drifter:
         kbu=kbu+1
         kbumax=kbumax+1
!        kbuglb_   =kbuglb_   +1
         if(kbu>nbomax)call drifter_error_message(1)
! Saisir la position en indices globaux dans tableau drifter_l:
         drifter_l(kbu,1)=posx_
         drifter_l(kbu,2)=posy_
         drifter_l(kbu,3)=posz_                                         !29-01-10
         drifter_l(kbu,4)=real( kbuglb_    )
         drifter_l(kbu,5)=elapsedtime_out !03-03-15
         if(flag_buo_w==1)drifter_l(kbu,id_wdrifter)=vertical_buoyancy_velocity ! bidouille patrick oU x10 est vitesse de la particule 

         endif                                   !\\\\\\\\\\>

      endif          !§§§§§§§§>


      goto 369

 363  close(3)

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

! mpi:
#ifdef parallele
      call drifter_exchange_driver
#endif

      end subroutine drifter_initial

!..................................................................................

      subroutine drifter_error_message(choice_   )
      use module_principal
      implicit none
      integer choice_   ,new_kbumax_
#ifdef synopsis
       subroutinetitle='drifter_error_message'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

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
       new_kbumax_   =                                             &
       kbumax+ nd_recv_est+nd_recv_ouest+nd_recv_nord+nd_recv_sud  &
             -(nd_send_est+nd_send_ouest+nd_send_nord+nd_send_sud)
       write(6,*)'vous demandez plus de drifter qu''il en est déclaré'
       write(6,*)'dans parameter, actuellement dimensionné à ',nbomax
       write(6,*)'mettre nbomax à au moins ',new_kbumax_
       stop ' STOP dans subroutine drifter_error_message cas 3'
      endif                  !333333333>

      if(choice_   ==4) then !444444444>
       stop 'choix 4'
      endif                  !444444444>

      end subroutine drifter_error_message

!..................................................................................

      subroutine drifter_initial_fortran
      use module_principal
      use module_parallele
      implicit none
      character*60 txt_spy_
#ifdef synopsis
       subroutinetitle='drifter_initial_fortran'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      txt_spy_   ='drifter_initial_fortran'

      kbu=0
      do j=1,jmax
      do i=1,imax
       if(mask_t(i,j,kmaxp1)==1) then

       if(mod(i+par%timax(1),3)==0.and.mod(j+par%tjmax(1),3)==0) then
        kbu=kbu+1
       endif

       endif
      enddo
      enddo
#ifdef parallele
      call barriere(iteration3d,4,txt_spy_   )                  !27-05-11
#endif
      if(kbu>nbomax)call drifter_error_message(1)

      kbu=0
      do j=1,jmax
      do i=1,imax
       if(mask_t(i,j,kmaxp1)==1) then
       if(mod(i+par%timax(1),3)==0.and.mod(j+par%tjmax(1),3)==0) then
        kbu=kbu+1
        drifter_l(kbu,1)=real(i)+par%timax(1)
        drifter_l(kbu,2)=real(j)+par%tjmax(1)
        drifter_l(kbu,3)=real(kmax)
!       drifter_l(kbu,4)=(real(i)-1)*jmax+real(j)+par%rank*imax*jmax
        drifter_l(kbu,4)=(real(i)-1)*jmax+real(j)+(par%rank  )*1000000
        drifter_l(kbu,5)=0. !03-03-15
       endif
       endif
      enddo
      enddo
      kbumax=kbu

      end subroutine drifter_initial_fortran

!..................................................................................
#ifdef parallele
      subroutine drifter_who_is_out
      use module_principal
      use module_parallele
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
#ifdef synopsis
       subroutinetitle='drifter_who_is_out'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

! Cette routine recense les numeros des drifters qui doivent migrer vers un
! domaine voisin:
      txt_spy_   ='drifter_who_is_out'

! Etape 1: identifier les drifters sortant
! et produire nd_send_nord nd_send_sud nd_send_est nd_send_ouest
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
      if(deci>real(imax).and.decj>real(jmax)) then !nenenene>

      if (par%tvoisin(nordest) /= mpi_proc_null) then !-mpn-mpn->

       nd_send_nordest=nd_send_nordest+1
        if(nd_send_nordest>nbobuffermax_c) then !>>>>>>
          warning_   =2
        else                                    !>>>>>>
         drifter_send_order_nordest(nd_send_nordest)=kbu
         drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                                   !>>>>>>

      else                                        !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out>nbobuffermax) then !>>>>>>
          warning_   =1
        else                              !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                             !>>>>>>

      endif                                       !-mpn-mpn->

      goto 1968
      endif                                        !nenenene>

! la zone d'echange Coin Nord-Ouest
      if(deci<1.and.decj>real(jmax)) then !nonononono>

      if (par%tvoisin(nordouest) /= mpi_proc_null) then !-mpn-mpn->

       nd_send_nordouest=nd_send_nordouest+1
        if(nd_send_nordouest>nbobuffermax_c) then !>>>>>>
          warning_   =2
        else                                    !>>>>>>
         drifter_send_order_nordouest(nd_send_nordouest)=kbu
         drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                                   !>>>>>>

      else                                              !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out>nbobuffermax) then !>>>>>>
          warning_   =1
        else                              !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                             !>>>>>>

      endif                                             !-mpn-mpn->

      goto 1968
      endif                               !nonononono>

! la zone d'echange Coin Sud-Ouest
      if(deci<1.and.decj<1.) then !sososososo>

      if (par%tvoisin(sudouest) /= mpi_proc_null) then !-mpn-mpn->

       nd_send_sudouest=nd_send_sudouest+1
        if(nd_send_sudouest>nbobuffermax_c) then !>>>>>>
          warning_   =2
        else                                    !>>>>>>
         drifter_send_order_sudouest(nd_send_sudouest)=kbu
         drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                                   !>>>>>>

      else                                              !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out>nbobuffermax) then !>>>>>>
          warning_   =1
        else                              !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                             !>>>>>>

      endif                                             !-mpn-mpn->

      goto 1968
      endif                       !sososososo>

! la zone d'echange Coin Sud-est
      if(deci>real(imax).and.decj<1.) then !sesesesese>

      if (par%tvoisin(sudest) /= mpi_proc_null) then !-mpn-mpn->

       nd_send_sudest=nd_send_sudest+1
        if(nd_send_sudest>nbobuffermax_c) then !>>>>>>
          warning_   =2
        else                                    !>>>>>>
         drifter_send_order_sudest(nd_send_sudest)=kbu
         drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                                   !>>>>>>

      else                                              !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out>nbobuffermax) then !>>>>>>
          warning_   =1
        else                              !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                             !>>>>>>

      endif                                             !-mpn-mpn->

      goto 1968
      endif                       !sesesesese>

! la zone d'echange Est:
      if(deci>real(imax)) then !eeeeeee>
      if (par%tvoisin(est) /= mpi_proc_null) then !-mpn-mpn->

        nd_send_est=nd_send_est+1
        if(nd_send_est>nbobuffermax) then !>>>>>>
          warning_   =1
        else                              !>>>>>>
         drifter_send_order_est(nd_send_est)=kbu
         drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                             !>>>>>>

      else                                        !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out>nbobuffermax) then !>>>>>>
          warning_   =1
        else                              !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                             !>>>>>>

      endif                                       !-mpn-mpn->
      goto 1968
      endif                    !eeeeeee>


! la zone d'echange Ouest:
      if(deci<1.)        then !ooooooo>
      if (par%tvoisin(ouest) /= mpi_proc_null) then !-mpn-mpn->

        nd_send_ouest=nd_send_ouest+1
        if(nd_send_ouest>nbobuffermax) then !>>>>>>>
          warning_   =1
        else                                !>>>>>>>
          drifter_send_order_ouest(nd_send_ouest)=kbu
          drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                               !>>>>>>>

      else                                        !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out>nbobuffermax) then !>>>>>>
          warning_   =1
        else                              !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                             !>>>>>>

      endif                                       !-mpn-mpn->
      goto 1968
      endif                    !ooooooo>


! la zone d'echange Nord:
      if(decj>real(jmax)) then !nnnnnnnn>
      if (par%tvoisin(nord) /= mpi_proc_null) then !-mpn-mpn->

        nd_send_nord=nd_send_nord+1
        if(nd_send_nord>nbobuffermax) then !>>>>>>>>
          warning_   =1
        else                               !>>>>>>>>
          drifter_send_order_nord(nd_send_nord)=kbu
          drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                              !>>>>>>>>

      else                                        !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out>nbobuffermax) then !>>>>>>
          warning_   =1
        else                              !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                             !>>>>>>

      endif                                       !-mpn-mpn->
      goto 1968
      endif                    !nnnnnnnn>

! la zone d'echange Sud:
      if(decj<1.)         then !ssssssss>
      if (par%tvoisin(sud) /= mpi_proc_null) then !-mpn-mpn->

        nd_send_sud=nd_send_sud+1
        if(nd_send_sud>nbobuffermax) then !>>>>>>>>>
          warning_   =1
        else                              !>>>>>>>>>
          drifter_send_order_sud(nd_send_sud)=kbu
          drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                             !>>>>>>>>>

      else                                        !-mpn-mpn->

        nd_send_out=nd_send_out+1
        if(nd_send_out>nbobuffermax) then !>>>>>>
          warning_   =1
        else                              !>>>>>>
         drifter_send_order_out(nd_send_out)=kbu
         drifter_l(kbu,4)=-drifter_l(kbu,4)
        endif                             !>>>>>>

      endif                                       !-mpn-mpn->
      goto 1968
      endif                    !ssssssss>

 1968 continue    ! SORTIE

      enddo       ! kbu loop

      if(warning_   ==1)call drifter_error_message(2)
      if(warning_   ==2)call drifter_error_message(4)

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

!     do kbu=1,kbumax
!       write(60+par%rank,*)kbu,drifter_l(kbu,4)
!     enddo

!.............................................................................................

      end subroutine drifter_who_is_out
#endif
!..................................................................................
#ifdef parallele
      subroutine drifter_obc(var_   )
      use module_principal
      use module_parallele
      implicit none
      integer,parameter :: tagouest_   =5000,  tagest_   =5010   &
                            ,tagsud_   =6000, tagnord_   =6010
      integer,parameter :: tagsudouest_   =15000, tagsudest_   =15010  &
                         ,tagnordouest_   =16000, tagnordest_   =16010 &
                         ,nexchgmax_   =16
!     integer nexchg_   ,warning_   ,var_
      integer nexchg_   ,var_
      integer,dimension(nexchgmax_   ) :: tabreq_
      integer,dimension(mpi_status_size,nexchgmax_   ) :: tstatus_
      character*60 txt_spy_
#ifdef synopsis
       subroutinetitle='drifter_obc'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      txt_spy_   ='drifter_obc'

!.............................................................................................
! Etape 1:
! Charger les tableaux d'echanges avec les drifters positionnes sur les bords des sous domaines:
      do k1=1,nd_send_nordest
       k2=drifter_send_order_nordest(k1)
       drifter_send_nordest(k1)=drifter_l(k2,var_   )
      enddo
      do k1=1,nd_send_nordouest
       k2=drifter_send_order_nordouest(k1)
       drifter_send_nordouest(k1)=drifter_l(k2,var_   )
      enddo
      do k1=1,nd_send_sudouest
       k2=drifter_send_order_sudouest(k1)
       drifter_send_sudouest(k1)=drifter_l(k2,var_   )
      enddo
      do k1=1,nd_send_sudest
       k2=drifter_send_order_sudest(k1)
       drifter_send_sudest(k1)=drifter_l(k2,var_   )
      enddo

      do k1=1,nd_send_est
       k2=drifter_send_order_est(k1)
       drifter_send_est(k1)=drifter_l(k2,var_   )
      enddo
      do k1=1,nd_send_ouest
       k2=drifter_send_order_ouest(k1)
       drifter_send_ouest(k1)=drifter_l(k2,var_   )
      enddo
      do k1=1,nd_send_nord
       k2=drifter_send_order_nord(k1)
       drifter_send_nord(k1)=drifter_l(k2,var_   )
      enddo
      do k1=1,nd_send_sud
       k2=drifter_send_order_sud(k1)
       drifter_send_sud(k1)=drifter_l(k2,var_   )
      enddo

!.............................................................................................
! Etape 2:
! Echanger les informations des zones d'echanges:

      nexchg_   =0

! Coin Nord-Est:
      if (par%tvoisin(nordest) /= mpi_proc_null) then
!     if (par%tvoisin(nord) /= mpi_proc_null) then
!     if (par%tvoisin(est)  /= mpi_proc_null) then
      nexchg_   =nexchg_   +1
      k=max(1,nd_send_nordest)
      call mpi_issend(drifter_send_nordest(1)    & ! envoyé au proc Est
                ,size(drifter_send_nordest(1:k)) &
                ,mpi_real                        &
                ,par%tvoisin(nordest)            &
                ,tagnordest_                     &
                ,par%comm2d                      &
                ,tabreq_   (nexchg_   )          &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nd_recv_nordest)
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
      k=max(1,nd_send_nordouest)
      call mpi_issend(drifter_send_nordouest(1)    & ! envoyé au proc Est
                ,size(drifter_send_nordouest(1:k)) &
                ,mpi_real                          &
                ,par%tvoisin(nordouest)            &
                ,tagnordouest_                     &
                ,par%comm2d                        &
                ,tabreq_   (nexchg_   )            &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nd_recv_nordouest)
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
      k=max(1,nd_send_sudouest)
      call mpi_issend(drifter_send_sudouest(1)    & ! envoyé au proc Est
                ,size(drifter_send_sudouest(1:k)) &
                ,mpi_real                         &
                ,par%tvoisin(sudouest)            &
                ,tagsudouest_                     &
                ,par%comm2d                       &
                ,tabreq_   (nexchg_   )           &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nd_recv_sudouest)
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
      k=max(1,nd_send_sudest)
      call mpi_issend(drifter_send_sudest(1)    & ! envoyé au proc Est
                ,size(drifter_send_sudest(1:k)) &
                ,mpi_real                       &
                ,par%tvoisin(sudest)            &
                ,tagsudest_                     &
                ,par%comm2d                     &
                ,tabreq_   (nexchg_   )         &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nd_recv_sudest)
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
      k=max(1,nd_send_est)
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
      k=max(1,nd_recv_est)
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
      k=max(1,nd_send_ouest)
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
      k=max(1,nd_recv_ouest)
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
      k=max(1,nd_send_nord)
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
      k=max(1,nd_recv_nord)
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
      k=max(1,nd_send_sud)
      call mpi_issend(drifter_send_sud(1)     &
                ,size(drifter_send_sud(1:k))  &
                ,mpi_real                        &
                ,par%tvoisin(sud)                            &
                ,tagsud_                                     &
                ,par%comm2d                                  &
                ,tabreq_   (nexchg_   )                      &
                ,ierr)

      nexchg_   =nexchg_   +1
      k=max(1,nd_recv_sud)
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


! Renumeroter en comblant les vides laissés par les drifters qui quittent le domaine:
#ifdef parallele
      call barriere(iteration3d,4,txt_spy_   )                  !27-05-11
#endif
      call drifter_incoming(var_   )

      end subroutine drifter_obc
#endif
!....................................................................

      subroutine drifter_incoming(var_   )
      use module_principal
      use module_parallele
      implicit none
      integer var_   ,new_kbumax_
#ifdef synopsis
       subroutinetitle='drifter_incoming'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      new_kbumax_   =kbumax                                            &
        +nd_recv_est+nd_recv_ouest+nd_recv_nord+nd_recv_sud            &
        +nd_recv_nordest+nd_recv_nordouest                             &
        +nd_recv_sudest +nd_recv_sudouest                              &
      -(nd_send_est+nd_send_ouest+nd_send_nord+nd_send_sud+nd_send_out &
       +nd_send_nordest+nd_send_nordouest                              &
       +nd_send_sudest +nd_send_sudouest  )

      if(new_kbumax_   >nbomax)call drifter_error_message(3)

! Les drifters entrants sont ajoutés dans les emplacements
! laissés vacants par les drifters sortants

      i1=nd_send_est ; i2=i1+nd_send_ouest ; i3=i2+nd_send_nord
      i4=i3+nd_send_sud
      i5=i4+nd_send_nordest
      i6=i5+nd_send_nordouest
      i7=i6+nd_send_sudouest
      i8=i7+nd_send_sudest
      i9=i8+nd_send_out

      j1=nd_recv_est ; j2=j1+nd_recv_ouest ; j3=j2+nd_recv_nord
      j4=j3+nd_recv_sud
      j5=j4+nd_recv_nordest
      j6=j5+nd_recv_nordouest
      j7=j6+nd_recv_sudouest
      j8=j7+nd_recv_sudest

      k=0
      k2=0
      k0=kbumax
      k1=kbumax

 10   k=k+1

      if(k>i9.and.k>j8)goto 20

! Donne le numero kbu d'un emplacement libre:
      if(         k<=i1)kbu=drifter_send_order_est(k)
      if(k>i1.and.k<=i2)kbu=drifter_send_order_ouest(k-i1)
      if(k>i2.and.k<=i3)kbu=drifter_send_order_nord(k-i2)
      if(k>i3.and.k<=i4)kbu=drifter_send_order_sud(k-i3)

      if(k>i4.and.k<=i5)kbu=drifter_send_order_nordest(k-i4)
      if(k>i5.and.k<=i6)kbu=drifter_send_order_nordouest(k-i5)
      if(k>i6.and.k<=i7)kbu=drifter_send_order_sudouest(k-i6)
      if(k>i7.and.k<=i8)kbu=drifter_send_order_sudest(k-i7)

      if(k>i8.and.k<=i9)kbu=drifter_send_order_out(k-i8)
      if(k>i9) then  !>>>>>>     ! Quand tous les espaces vacants sont comblés
       k0=k0+1                   ! et qu'il reste des entrants (c.a.d. k<=j4)
       kbu=k0                    ! on les ajoute en bout de liste. Le nbre
      endif          !>>>>>>     ! total de drifter augmente car k0>kbumax

! Les derniers emplacements ne doivent pas être modifiés pour que l'algo
! du cas k>j4 fonctionne normalement. Par consequent je demande un
! autre kbu si kbu>new_kbumax_   . Par ailleurs, les drifters vacants
! déjà en dehors des limites sont tres bien où ils sont: qu'il y restent.
      if(kbu>new_kbumax_   )goto 10
      k2=k2+1

!     if(par%rank==1)write(*,*)'emplacement libre:',kbu

! Attribut la valeur d'un drifter entrant:
      if(          k2<=j1)drifter_l(kbu,var_   )=drifter_recv_est(k2)
      if(k2>j1.and.k2<=j2)drifter_l(kbu,var_   )=drifter_recv_ouest(k2-j1)
      if(k2>j2.and.k2<=j3)drifter_l(kbu,var_   )=drifter_recv_nord(k2-j2)
      if(k2>j3.and.k2<=j4)drifter_l(kbu,var_   )=drifter_recv_sud(k2-j3)
      if(k2>j4.and.k2<=j5)drifter_l(kbu,var_   )=drifter_recv_nordest(k2-j4)
      if(k2>j5.and.k2<=j6)drifter_l(kbu,var_   )=drifter_recv_nordouest(k2-j5)
      if(k2>j6.and.k2<=j7)drifter_l(kbu,var_   )=drifter_recv_sudouest(k2-j6)
      if(k2>j7.and.k2<=j8)drifter_l(kbu,var_   )=drifter_recv_sudest(k2-j7)
      if(k2>j8) then  !)))))))>

        do while (drifter_l(k1,4)<0)
         k1=k1-1
        enddo
        if(var_   ==4) then !------>
           drifter_l(k1,var_   )=-drifter_l(k1,var_   )
        endif               !------>
        drifter_l(kbu,var_   )=drifter_l(k1,var_   )
        k1=k1-1
      endif          !))))))>

      goto 10


   20 continue

      end subroutine drifter_incoming

!.................................................................

      subroutine drifter_new_kbumax
      use module_principal
      use module_parallele
      implicit none
      integer kbumax_b_glb_   ,kbumax_a_glb_   ,kbu_exchg_glb_
      character*60 txt_spy_
      txt_spy_   ='drifter_new_kbumax'
#ifdef synopsis
       subroutinetitle='drifter_new_kbumax'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif


#ifdef parallele
! Nbre global de drifter avant mise à jour ("b"efore)
!     call mpi_allreduce(kbumax,kbumax_b_glb_   ,1,mpi_integer,mpi_sum,par%comm2d,ierr)
!     call mpi_allreduce(nd_send_est+nd_send_ouest+nd_send_nord+nd_send_sud &
!             ,kbu_exchg_glb_   ,1,mpi_integer,mpi_sum,par%comm2d,ierr)
#endif

! Mise à jour du nombre total de drifter par subdomaine
      kbumax=kbumax                                                    &
        +nd_recv_est+nd_recv_ouest+nd_recv_nord+nd_recv_sud            &
        +nd_recv_nordest+nd_recv_nordouest                             &
        +nd_recv_sudest +nd_recv_sudouest                              &
      -(nd_send_est+nd_send_ouest+nd_send_nord+nd_send_sud+nd_send_out &
       +nd_send_nordest+nd_send_nordouest                              &
       +nd_send_sudest +nd_send_sudouest  )

      if(kbumax>nbomax)call drifter_error_message(1)

! Retablir le signe positif de la numerotation:
      do kbu=1,kbumax
       drifter_l(kbu,4)=abs(drifter_l(kbu,4))
      enddo

      end subroutine drifter_new_kbumax

!.................................................................
#ifdef parallele
      subroutine drifter_exchange_driver
      use module_principal
      implicit none
      integer loop_
#ifdef synopsis
       subroutinetitle='drifter_exchange_driver'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

      call drifter_who_is_out ! recense les drifter à sortir

      call drifter_obc(1)     ! changement de subdomaine variable 1
      call drifter_obc(2)     ! changement de subdomaine variable 2
      call drifter_obc(3)     ! changement de subdomaine variable 3

      do loop_=5,drifter_dim2  ! changement de subdomaine variable 5 a drifter_dim2
       call drifter_obc(loop_) !03-03-15
      enddo


      call drifter_obc(4) ! DANGER! La variable 4 s'echange en dernier!

      call drifter_new_kbumax

!     call drifter_gnuplot_output

      if(drifter_out_sampling>0.)call drifter_ascii_file

      end subroutine drifter_exchange_driver
#endif
!.................................................................
#ifdef parallele
      subroutine drifter_gnuplot_output
      use module_principal
      use module_parallele
      implicit none
#ifdef synopsis
       subroutinetitle='drifter_gnuplot_output'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

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
#endif
!..................................................................................

      subroutine drifter_ascii_file
      use module_principal
      use module_parallele
      implicit none

      if(kbumax==0)return

! archivage des trajectoires dans fichiers individuels ascii
      if(  int(elapsedtime_now/drifter_out_sampling)           &
         /=int(elapsedtime_bef/drifter_out_sampling) ) then  !-------->

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
        write(3,'(e12.5,6(1x,e13.6))')  & !11/03/09
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
        write(3,'(                     &
                 i4,1x,5(i2,1x)        & 
                ,f6.2,1x               & 
                ,i10                   & 
                ,3(1x,f13.6),1x        &
                ,f13.6,1x              & 
                ,f9.6,1x               & 
                 )')                   &
        year_now,month_now,day_now     &
       ,hour_now,minute_now,second_now & ! date
       ,elapsedtime_now/86400.         & ! temps EcoulE en jours
       ,nint(drifter_l(kbu,4),kind=8)  & ! IDentite du drifter
       ,x1*rad2deg,x2*rad2deg,x3       & ! longitude, latitude, profondeur
       ,drifter_l(kbu,5)               & ! temps A la position initiale du drifter
       ,drifter_l(kbu,id_wdrifter)       ! vitesse individuelle de flottabilitE
      endif                            !rank-files>

  145 continue

      if(drifter_output_files==1)close(3) !rank-files>

      endif                                !--------->

      end subroutine drifter_ascii_file

!.........................................................................................

      subroutine drifter_allocate
      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='drifter_allocate'
       subroutinedescription=' '
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

       allocate(drifter_l                   (nbomax,drifter_dim2));drifter_l=0. !03-03-15
       allocate(drifter_send_nord           (nbobuffermax)      )
       allocate(drifter_recv_nord           (nbobuffermax)      )
       allocate(drifter_send_sud            (nbobuffermax)      )
       allocate(drifter_recv_sud            (nbobuffermax)      )
       allocate(drifter_send_est            (nbobuffermax)      )
       allocate(drifter_recv_est            (nbobuffermax)      )
       allocate(drifter_send_ouest          (nbobuffermax)      )
       allocate(drifter_recv_ouest          (nbobuffermax)      )
       allocate(drifter_send_nordest        (nbobuffermax_c)    )
       allocate(drifter_send_nordouest      (nbobuffermax_c)    )
       allocate(drifter_send_sudest         (nbobuffermax_c)    )
       allocate(drifter_send_sudouest       (nbobuffermax_c)    )
       allocate(drifter_recv_nordest        (nbobuffermax_c)    )
       allocate(drifter_recv_nordouest      (nbobuffermax_c)    )
       allocate(drifter_recv_sudest         (nbobuffermax_c)    )
       allocate(drifter_recv_sudouest       (nbobuffermax_c)    )
       allocate(drifter_send_order_nord     (nbobuffermax)      )
       allocate(drifter_send_order_sud      (nbobuffermax)      )
       allocate(drifter_send_order_est      (nbobuffermax)      )
       allocate(drifter_send_order_ouest    (nbobuffermax)      )
       allocate(drifter_send_order_out      (nbobuffermax)      )
       allocate(drifter_send_order_nordest  (nbobuffermax_c)    )
       allocate(drifter_send_order_nordouest(nbobuffermax_c)    )
       allocate(drifter_send_order_sudest   (nbobuffermax_c)    )
       allocate(drifter_send_order_sudouest (nbobuffermax_c)    )

      end subroutine drifter_allocate

!...............................................
#ifdef bidon
      subroutine drifter_isobathe
      use module_principal
      use module_parallele
      implicit none

       do kbu=1,kbumax
        deci=drifter_l(kbu,1)-par%timax(1)
        decj=drifter_l(kbu,2)-par%tjmax(1)
        i=int(deci)
        j=int(decj)
        rapi=deci-real(i)
        rapj=decj-real(j)
        drifter_l(kbu,6)=(1.-rapi)*(1.-rapj)*h_w(i  ,j  )   &
                        +(1.-rapi)*    rapj *h_w(i  ,j+1)   &
                        +    rapi *(1.-rapj)*h_w(i+1,j  )   &
                        +    rapi *    rapj *h_w(i+1,j+1)
       enddo
      end subroutine drifter_isobathe
#endif

!...............................................

      end module module_drifter


! Une routine pour tester la normalisatio du nombre aleatoire
#ifdef bidon
        program test_random
        implicit none
        real :: x3_r4
        double precision :: valmean , sum1=0. , tken=5. , w0 , w , sum2=0. &
         ,sum3=0.,sum4=0.
        integer :: loop

! Ce programme constitue une serie aleatoire de vitesse verticale w
! Son energie cinetique est 0.5*(w**2).  On verifie que la moyenne de 
! cette energie cinetique est bien tken/3 comme attendu pour peu que le
! nombre aleatoire x3_r4 produit par la fonction fortran random_number 
! est ensuite transforme en 2.*sqrt(3.)*(x3_r4-0.5)
        w0=sqrt(2.*tken/3.)
        sum1=0.
        sum4=0.
        do loop=1,10000000
         call random_number(x3_r4)
         x3_r4=2.*sqrt(3.)*(x3_r4-0.5)
         w=w0*x3_r4
         sum1=sum1+1.
         sum4=sum4+0.5*w**2
        enddo
        write(6,*)'moyenne 0.5*(w**2)',sum4/sum1
        write(6,*)'tken/3.           ',tken/3.

        end
#endif
