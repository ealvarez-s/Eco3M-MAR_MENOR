










      subroutine wetdry_mask_airseafluxes
!______________________________________________________________________
! SYMPHONIE model
! release 364 - last update: 17-01-23
!______________________________________________________________________

      use module_principal
      implicit none

!______________________________________________________________________
! version   date     Description des modifications
! 2010.3    05-01-10 mis en service
!           12-01-10 wetmask_t pour T et S
! 2010.8    18-03-10 Extension des boucles pour courant de Stokes aux
!                    frontieres
! S26       31-10-15 ajout de precipi_w
!           12-01-16 masque traceur: pour proprietes de conservation du
!                    schema d'advection sous-itErE le masque des traceurs
!                    est soit 0 soit 1 mais pas entre les deux
!           16-07-16 modif de la fonction de decroissance de wetmask_t
!           15-08-16 min(hz_w(0),hz_w(2)) pour calcul wetmask
!           20-09-16 wetmask multipliE par mask afin que multiplier par
!                    wetmask_u ou wetmask_v wetmask englobe aussi le masque fixe
!           29-04-17 wetmask_t: minimum de 3 echeances Leap-Frog
!           19-01-19 dans le cas sans mode splitting, wetmask_u et _v
!                    sont calculEs avec le meme algo que wetmask_u et _v
!                    du mode externe
! v247      25-02-19 modif dans le cadre du schema 6 du wetdrying barocline
! v247      26-02-19 pour un calcul correct de dti_ts max il est important
!                    de multiplier wetmask_t par la masque terre/mer sinon
!                    le point "traceur" en aval d'un point source de riviere
!                    contraint inutilement le pas de temps des traceurs
!           28-02-19 amenagements pour cas NH
! v257      19-06-19 wetmask_wi_t permet de passer A un schema d'advection verticale implicite
!           20-06-19 Attention au role du land/sea mask dans le calcul de wetmask_wi_t
! v261      17-10-19 wetmask_u et v du cas nosplitting multipliEs par mask_u et v
! v327      19-02-22 Ne pas empecher les precipitations sur la zone intertidale
! v332      28-02-22 upwindwetdry_t(i,j) =1 tant que h+ssh>2*wetdry_cst1 
! v343      03-04-22  - il manquait le = dans <= de if(min(min(hz_w(i,j,2),hz_w(i,j,0)),hz_w(i,j,1))<=wetdry_cst3)
!                     - ajout calcul wetmask_t pour cas nonsplitting
! v364      17-01-23 boucles augmentees pour schema d'advection horizontale implicite
!______________________________________________________________________
!    _________                    .__                  .__             !m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !
!..............................................................................

!$ This routines computes the wetdry mask for the 3D grid and then
!$ cancels the airsea fluxes above dried aeras.

!$ Wetdry mask for baroclinic velocities:
      do j=0,jmax+1   !18-03-10
      do i=0,imax+1
!      wetmask_t(i,j)=max(min( hz_w(i,j,2)/wetdry_cst1-1. , un ),zero)
!      wetmask_t(i,j)=max(min( hz_w(i,j,2)/wetdry_cst1    , 1.),0.) !16-07-16
!      wetmask_t(i,j)=max(min( min(hz_w(i,j,2),hz_w(i,j,0))/wetdry_cst1    , 1.),0.) !15-08-16
!      wetmask_t(i,j)=max(min( min(min(hz_w(i,j,2),hz_w(i,j,0)),hz_w(i,j,1)) /wetdry_cst1    , 1.),0.) !29-04-17

       wetmask_t(i,j)=max(min( min(min(hz_w(i,j,2),hz_w(i,j,0)),hz_w(i,j,1)) /wetdry_cst1 -0.1  , 1.),0.) !25-02-19

! upwindwetdry_t(i,j) =1 tant que h+ssh>2*wetdry_cst1 !28-02-22
! upwindwetdry_t(i,j) =0 tant que h+ssh<  wetdry_cst1 ! => adv T,S 100% schema upwind
       upwindwetdry_t(i,j)=max(min( min(hz_w(i,j,2),hz_w(i,j,0),hz_w(i,j,1)) /(2*wetdry_cst1) -0.5  , 1.),0.) 

      enddo
      enddo

! Pour le cas SANS mode splitting voir subroutine wetdry_mask_nosplitting_uv !19-01-19
! pour le calcul de wetmask_u et wetmask_v
! et https://docs.google.com/presentation/d/1c361LU5mr_Q6IGKD2ipFJ5qnSvVNygkVajIz1ggZSyo/edit#slide=id.g4c9ed15b6c_0_1
      if(flag_timesplitting==1)  then  !-time-splitting-case-> !19-01-19

!$ Wetdry mask for u:
       do j=1,jmax    !18-03-10
       do i=1,imax+1
        wetmask_u(i,j)=min(wetmask_t(i-1,j),wetmask_t(i,j))*mask_u(i,j,kmax) !20-09-16
       enddo
       enddo
!$ Wetdry mask for v:
       do j=1,jmax+1  !18-03-10
       do i=1,imax
        wetmask_v(i,j)=min(wetmask_t(i,j-1),wetmask_t(i,j))*mask_v(i,j,kmax) !20-09-16
       enddo
       enddo

!     else                             !-time-splitting-case-> !19-01-19
   
! Lignes commentees le 28-10-19 car cela changait la valeur de wetmask_u et wetmask_v dans la
! conditions aux limites ouvertes pour q
! On passe par ici si pas de mode splitting et pour etape traceurs et turbulence car l'etape !28-02-19
! wetmask_u wetmask_v des vitesses baroclines se produit dans la subroutine
! wetdry_mask_nosplitting_uv. On applique wetmask_u=1 pour ne pas reduire les
! flux veldxdz veldydz de l'etape d'advection des traceurs comme c'est la cas actuellement
! du cas hydrostatique (en periode d'essai)
!       wetmask_u=1. !28-02-19
!       wetmask_v=1.
     
      endif                            !-time-splitting-case-> !19-01-19

! En zone de decouvrement (dEs le debut du processus) l'advection verticale de T,S est implicite
      do j=1,jmax ; do i=1,imax !19-06-19
! la manip capilotractee sur mask a pour but de contrer le fait que wetmask=0 dans le masque continental !20-06-19
! qui conduirait (si on ne le prenant pas en consideration) A avoir wetmask_wi_t=0 en tout point contigu
! de la cote (entre autres les canaux des rivieres)
       wetmask_wi_t(i,j)=(int( wetmask_u(i  ,j  ))*mask_u(i  ,j  ,kmax)+1-mask_u(i  ,j  ,kmax))  &
                        *(int( wetmask_u(i+1,j  ))*mask_u(i+1,j  ,kmax)+1-mask_u(i+1,j  ,kmax))  &
                        *(int( wetmask_v(i  ,j  ))*mask_v(i  ,j  ,kmax)+1-mask_v(i  ,j  ,kmax))  &
                        *(int( wetmask_v(i  ,j+1))*mask_v(i  ,j+1,kmax)+1-mask_v(i  ,j+1,kmax))   
      enddo         ; enddo

!$ cancels the airsea fluxes in dried aeras:
      do j=1,jmax
      do i=1,imax
           ssr_w(i,j,1)=    ssr_w(i,j,1)*wetmask_t(i,j)
          snsf_w(i,j,1)=   snsf_w(i,j,1)*wetmask_t(i,j)
          slhf_w(i,j,1)=   slhf_w(i,j,1)*wetmask_t(i,j)
          sshf_w(i,j,1)=   sshf_w(i,j,1)*wetmask_t(i,j)
! CommentEe le 19-02-22 pour ne pas empecher les precipitations sur la zone intertidale
!      precipi_w(i,j,1)=precipi_w(i,j,1)*wetmask_t(i,j) !31-10-15
       wstress_w(i,j)=  wstress_w(i,j)  *wetmask_t(i,j)
      enddo
      enddo

! Lignes commentees le 28-02-19 car a priori inutiles
      do j=2,jmax-1
      do i=2,imax
       wstress_u(i,j,1)=wstress_u(i,j,1)*wetmask_u(i,j)
      enddo
      enddo

      do j=2,jmax
      do i=2,imax-1
       wstress_v(i,j,1)=wstress_v(i,j,1)*wetmask_v(i,j)
      enddo
      enddo

!$ Wetdry mask for T & S:
!     do j=1,jmax   ; do i=1,imax
      do j=0,jmax+1 ; do i=0,imax+1 !boucles augmentees pour schema d'advection horizontale implicite !17-01-23
        if(min(min(hz_w(i,j,2),hz_w(i,j,0)),hz_w(i,j,1))<=wetdry_cst3) then !03-04-22
          wetmask_t(i,j)=0.
        else
          wetmask_t(i,j)=1.*mask_t(i,j,kmax) !26-02-19
        endif
      enddo ; enddo

      end subroutine wetdry_mask_airseafluxes

!.......................................................................

      subroutine wetdry_mask_nosplitting_uv !19-01-19
      use module_principal
      use module_parallele
      implicit none

! Dans le cas sans mode splitting, wetmask_u et _v sont calculEs avec le meme 
! algo que wetmask_u et _v du mode externe
! Algo "schema 4". Details dans:
! https://docs.google.com/presentation/d/1c361LU5mr_Q6IGKD2ipFJ5qnSvVNygkVajIz1ggZSyo/edit#slide=id.p

!$ Wetdry mask for u:
!      do j=1,jmax ; do i=1,imax+1
!        wetmask_u(i,j)=max(min( & !19-01-19
!        ( min(h_w(i,j),h_w(i-1,j)) +0.5*min(ssh_int_w(i,j,0)+ssh_int_w(i-1,j,0)  &
!                                           ,ssh_int_w(i,j,1)+ssh_int_w(i-1,j,1)) &
!        )/wetdry_cst1                                                            &
!                             ,1.),0.)
!      enddo ; enddo
!$ Wetdry mask for v:
!      do j=1,jmax+1 ; do i=1,imax
!        wetmask_v(i,j)=max(min( & !19-01-19
!        ( min(h_w(i,j),h_w(i,j-1)) +0.5*min(ssh_int_w(i,j,0)+ssh_int_w(i,j-1,0)  &
!                                           ,ssh_int_w(i,j,1)+ssh_int_w(i,j-1,1)) &
!        )/wetdry_cst1                                                            &
!                             ,1.),0.)
!      enddo ; enddo

! https://docs.google.com/presentation/d/1c361LU5mr_Q6IGKD2ipFJ5qnSvVNygkVajIz1ggZSyo/edit#slide=id.p
! Schema 5:

!$ Wetdry mask for u:
      do j=1,jmax ; do i=1,imax+1

       wetmask_u(i,j)=max(min(                                          & !28-02-19
         ( min(h_w(i,j)                   ,h_w(i-1,j)                 ) &
      +0.5*max(ssh_int_w(i,j,0)+ssh_int_w(i,j,1),ssh_int_w(i-1,j,0)+ssh_int_w(i-1,j,1)) &
         )/wetdry_cst1                                                  &
      -0.1 & ! Cette constante pour atteindre 0 avant que la hauteur soit nulle: Dz/wetdry_cst1-0.5=0 pour Dz=0.5*wetdry_cst1
                              ,1.),0.) &
         *mask_u(i,j,kmax) !17-10-19

      enddo ; enddo

!$ Wetdry mask for v:
      do j=1,jmax+1 ; do i=1,imax

       wetmask_v(i,j)=max(min(                                          & !28-02-19
         ( min(h_w(i,j)                   ,h_w(i,j-1)                 ) &
      +0.5*max(ssh_int_w(i,j,0)+ssh_int_w(i,j,1),ssh_int_w(i,j-1,0)+ssh_int_w(i,j-1,1)) &
         )/wetdry_cst1                                                  &
      -0.1 & ! Cette constante pour atteindre 0 avant que la hauteur soit nulle: Dz/wetdry_cst1-0.5=0 pour Dz=0.5*wetdry_cst1
                              ,1.),0.) &
         *mask_v(i,j,kmax) !17-10-19

      enddo ; enddo

! A priori wetmask_t(i,j) ne sert A rien si T et S ne sont pas calculEes !03-04-22
! sauf A la rigueur dans graph_out.F90 pour masquer des champs....
!$ Wetdry mask for T & S:
      do j=1,jmax
      do i=1,imax

        if(min(min(hz_w(i,j,2),hz_w(i,j,0)),hz_w(i,j,1))<=wetdry_cst3) then !29-04-17
          wetmask_t(i,j)=0.
        else
          wetmask_t(i,j)=1.*mask_t(i,j,kmax) !26-02-19
        endif

      enddo
      enddo

      end subroutine wetdry_mask_nosplitting_uv !19-01-19
