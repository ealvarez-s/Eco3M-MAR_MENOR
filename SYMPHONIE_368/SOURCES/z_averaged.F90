      subroutine z_averaged(i1_,i2_,j1_,j2_,ki2)
!______________________________________________________________________
! SYMPHONIE ocean model
! release 285 - last update: 05-06-20
!______________________________________________________________________

      use module_principal
      implicit none
      integer i1_,i2_,j1_,j2_,ki2                          !07/06/07
#ifdef synopsis
       subroutinetitle='z_averaged'
       subroutinedescription= &
       'Computes the depth-averaged current and stores the result in' &
       //' velavr_u and velavr_v' 
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!.......................................................................
! Version date      Description des modifications
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         11/08/06: appel à vm_to_vmhz supprimé
!         07/06/07: * l'argument KI2 designe l'iteration consideree
!                   et la sommation porte sur vhz, une boucle 2D
!                   etant rajoutee pour revenir au courant moyen
!                   * ajout du choix de bornes min et max pour le deroulement
!                   des boucles i et j
!         12/06/07: utilisation des tableaux hzdx_y et hzdy_x
! 2009.3  20-10-09  Avec la grille generique, on reprend la demi-somme des
!                   niveaux qui avait été remise en cause sur les marches
!                   d'escalier de la version sigma/step. Cela signifie que
!                   la premiere cellule de courant "sous le fond" n'a pas une
!                   epaisseur nulle. Cela implique que vmean_ n'a plus tout
!                   à fait le même sens dans ce cas de figure. Pour une grille
!                   "classique" c'est la moyenne du courant, pour une grille
!                   hybride c'est la valeur "barotrope" des cellules non masquées.
!                   Dans ce cas la valeur vmean est le courant moyen amplifiée au
!                   prorata du ration H/H0 où H est l'epaisseur totale et
!                   H0 (zeroleveldepth) l'epaisseur où le courant est non nul.
! 2010.7  28-02-10  la moyenne est calculée à partir de vel car veldxdz
!                   peut contenir le courant de stokes
! 2010.8  21-05-10  - Un calcul plus rapide pour zeroleveldepth
!                   - Le calcul du courant moyen de stokes est identique à celui
!                     du courant moyen
!                   - les entiers passés en arguments ne sont plus utilises
!         26-05-10  correction sur calcul courant moyen de stokes
! S.26    09-04-13  zeroleveldepth est desormais calcule plus tot, c.a.d. dans
!                   z_levels afin de pouvoir etre utilise dans d'autres routines
!                   Cela permet notamment de calculer velbarstokes dans le module_wave
!                   et non plus dans la presente routine
!         07-04-14  boucle stricte
!         05-08-15  methode vst continue
!         21-04-17  coordonnee s-z
! v285    05-06-20  suppression zeroleveldepth_u et v
!...............................................................................
!    _________                    .__                  .__                     ! 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____               !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \              ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/              !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >             !
!          \/\/          \/|__|        \/            \/        \/              !
!...............................................................................

      if(flag3d.eq.1) then !111111111111111111>

       sum1=0

      k0=0
      do k=1,kmax

! si method vst discontinue k1=k sinon (sigma ou vst continue) k1=kmax !05-08-15
       k1=vststep*k+(1-vststep)*kmax 
       do j=1,jmax
       do i=1,imax+1
          velavr_u(i,j,1)=velavr_u(i,j,1)*k0               &
                            +vel_u(i,j,k,1)                &
                             *dz_u(i,j,k,1)                &
                           *mask_u(i,j,k ) 
!                          *mask_u(i,j,k1) 
       enddo
       enddo
       do j=1,jmax+1
       do i=1,imax
          velavr_v(i,j,1)=velavr_v(i,j,1)*k0               &
                            +vel_v(i,j,k,1)                &
                             *dz_v(i,j,k,1)                &
                           *mask_v(i,j,k ) 
!                          *mask_v(i,j,k1) 
       enddo
       enddo
       k0=1
      enddo ! fin de boucle sur k

      do j=1,jmax !07-04-14
      do i=1,imax+1
         velavr_u(i,j,1)=velavr_u(i,j,1)               & !28-02-10
                            /hz_u(i,j,1) !05-06-20
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
         velavr_v(i,j,1)=velavr_v(i,j,1)               & !28-02-10
                             /hz_v(i,j,1)
      enddo
      enddo

#ifdef bidon
! Calcul du courant moyenné sur la verticale:
      k=1
      do j=1,jmax
      do i=1,imax+1
          velavr_u(i,j,1)=   vel_u(i,j,k,1)                &
                             *dz_u(i,j,k,1)
      enddo
      enddo
      do j=1,jmax
      do i=1,imax+1
      do k=2,kmax
          velavr_u(i,j,1)=velavr_u(i,j,1)                  &
                            +vel_u(i,j,k,1)                &
                             *dz_u(i,j,k,1)
      enddo
      enddo
      enddo
      k=1
      do j=1,jmax+1
      do i=1,imax
          velavr_v(i,j,1)=   vel_v(i,j,k,1)                &
                             *dz_v(i,j,k,1)
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
      do k=2,kmax
          velavr_v(i,j,1)=velavr_v(i,j,1)                  &
                            +vel_v(i,j,k,1)                &
                             *dz_v(i,j,k,1)
      enddo
      enddo
      enddo
      do j=1,jmax !07-04-14
      do i=1,imax+1
         velavr_u(i,j,1)=velavr_u(i,j,1)               & !28-02-10
                            /hz_u(i,j,1)
      enddo
      enddo
      do j=1,jmax+1
      do i=1,imax
         velavr_v(i,j,1)=velavr_v(i,j,1)               & !28-02-10
                            /hz_v(i,j,1)
      enddo
      enddo
#endif



      return
      else             !1111111111000000000>


! note: pour I1D=0 passer par z_averaged permet d'annuler VEL(NR)
! ce qui est indispensable pour un calcul 1D correct de coriolis

      do j=1,jmax+1
        do i=1,imax+1
         velavr_u(i,j,1)=0.
         velavr_v(i,j,1)=0.
       enddo
      enddo

      return
      endif             !000000000000000000>

      end subroutine z_averaged
