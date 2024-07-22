










      subroutine adve3dto2d ! version centree
!______________________________________________________________________
! SYMPHONIE ocean model
! release 321 - last update: 16-01-22
!______________________________________________________________________
      use module_principal
      use module_parallele
      implicit none

! computes the frozen advective terms used in the external mode

!..............................................................................
! Version date      Description des modifications
!         01/11/01: inversion de l'ordre des boucles: J passe avant I
!         11/08/06: Le tableau vmeanhz_x et vmeanhz_y occupaient un
!                   emplacement memoire inutile. Des tableaux temporaires
!                   2D les remplacent desormais.
!         20/04/07: Passage à la coordonnée curviligne
!         23/04/07: suite du point precedent
! 2009.3  30-09-09  utilisation des nouveaux facteurs d'echelle verticale
! 2010.7  28-02-10  ajout velstokes et velbarstokes
! 2010.8  03-05-10  tableaux advx renommés xflux etc..., division par hzdy_u etc...
! S.26    06-12-13  application boucles strictes
!         10-11-14  possibilite schema upwind. Par defaut schema centre
!         28-11-14  remplacement hzdx par son produit hz*dx
!         05-08-15  methode vst continue
!         03-09-19  si flag_timesplitting_adve_uv==1 adve3d2d devient la moyenne verticale de l'advection 3D
! v321    16-01-22  stop option flag_timesplitting_adve_uv=1
!...............................................................................
!    _________                    .__                  .__             ! (°o°) 
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       ! m°v°m 
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      ! 
!...............................................................................


! ATTENTION flag_timesplitting_adve_uv=0 EST LE SCHEMA PAR DEFAUT

      if(flag_timesplitting_adve_uv==1) then !m°v°m>
!/////////////////////////////////////////////////////////////////////
! ADVECTION TOTALE AVANT MODE EXTERNE

      stop & !16-01-22
      ' NE PAS PASSER ICI (flag_timesplitting_adve_uv=1 dans adve3D2D)'

      endif                                  !m°v°m>

!/////////////////////////////////////////////////////////////////////
! ADVECTION PARTIELLE AVANT MODE EXTERNE (termes de couplage uniquement)
      if(flag_timesplitting_adve_uv==0) then !m°v°m>

! ATTENTION SI ON REVIENT A L'ANCIEN SCHEMA IL FAUT TENIR QUAND MEME COMPTE DE CE QUE
! FROZEN EST DIFFERENT PAR rapport A la prise en compte des facteurs d'echelle dans la nouvelle
! version

! La divergence du flux prime est 1/dxdy*( d(dydzu'u')/di +d(dxdzv'u') )
! Les lignes suivantes ne calcule que la partie entre parenthses soit d(dydzu'u')/di +d(dxdzv'u')
! Il faudra appliquer le facteur 1/dxdy sachant que dans l'ancienne version la division par dx
! se faisait dans le calcul de frozenterme en consequence de quoi on ne divisait que par dy A
! la fin de cette routine. Dans l'ancienne version frozenterm faisait la multiplication par hz
! de sorte que par anticipation, adve3d2d devait etre une moyenne et non une integrale en consequence
! de quoi on anticipait A la fin de cette routine par une division par H. Avec la nouvelle version
! frozen n'applique pas les operations susmentionnees et donc: 
! 1) il faut diviser par dxdy A la fin et 2) on ne divise pas par HZ A la fin

      do j=2,jmax-1
      do i=2,imax
       adve3d2d_u(i,j)=0.
      enddo
      enddo
      do j=2,jmax
      do i=2,imax-1
       adve3d2d_v(i,j)=0.
      enddo
      enddo

      do 14 k=1,kmax

      k1=vststep*k+(1-vststep)*kmax  !05-08-15

      do j=1,jmax
      do i=1,imax+1

       xy_u(i,j,1)=vel_u(i,j,k,1)-velavr_u(i,j,1)*mask_u(i,j,k1)

       xy_u(i,j,2)=(vel_u(i,j,k,1) +velstokes_u(i,j,k,1)              & !28-02-10
               -(velavr_u(i,j,1)+velbarstokes_u(i,j,1))*mask_u(i,j,k1) &
                   )*dz_u(i,j,k,now)                                  &
                    *dy_u(i,j)
      enddo
      enddo

      do j=1,jmax+1
      do i=1,imax

       xy_v(i,j,1)=vel_v(i,j,k,1)-velavr_v(i,j,1)*mask_v(i,j,k1)

       xy_v(i,j,2)=(vel_v(i,j,k,1) +velstokes_v(i,j,k,1)              & !28-02-10
               -(velavr_v(i,j,1)+velbarstokes_v(i,j,1))*mask_v(i,j,k1) &
                   )*dz_v(i,j,k,now)                                  &
                    *dx_v(i,j)

      enddo
      enddo

      do j=2,jmax-1
      do i=2,imax

! interne force externe via advection pour equation U
       adve3d2d_u(i,j)=                                                &
       adve3d2d_u(i,j)+(                                               &

                   (xy_u(i  ,j  ,2)+xy_u(i+1,j  ,2))                   &
                  *(xy_u(i  ,j  ,1)+xy_u(i+1,j  ,1))                   &

                 - (xy_u(i-1,j  ,2)+xy_u(i  ,j  ,2))                   &
                  *(xy_u(i-1,j  ,1)+xy_u(i  ,j  ,1))                   &

                 + (xy_v(i-1,j+1,2)+xy_v(i  ,j+1,2))                   &
                  *(xy_u(i  ,j  ,1)+xy_u(i  ,j+1,1))                   &

                 - (xy_v(i-1,j  ,2)+xy_v(i  ,j  ,2))                   &
                  *(xy_u(i  ,j-1,1)+xy_u(i  ,j  ,1))                   &

                                                      )*mask_u(i,j,k1)

      enddo
      enddo

      do j=2,jmax
      do i=2,imax-1

! interne force externe via advection pour equation V
       adve3d2d_v(i,j)=                                                &
       adve3d2d_v(i,j)+(                                               &

                   (xy_u(i+1,j  ,2)+xy_u(i+1,j-1,2))                   &
                  *(xy_v(i+1,j  ,1)+xy_v(i  ,j  ,1))                   &

                 - (xy_u(i  ,j  ,2)+xy_u(i  ,j-1,2))                   &
                  *(xy_v(i  ,j  ,1)+xy_v(i-1,j  ,1))                   &

                 + (xy_v(i  ,j+1,2)+xy_v(i  ,j  ,2))                   &
                  *(xy_v(i  ,j+1,1)+xy_v(i  ,j  ,1))                   &

                 - (xy_v(i  ,j  ,2)+xy_v(i  ,j-1,2))                   &
                  *(xy_v(i  ,j  ,1)+xy_v(i  ,j-1,1))                   &

                                                      )*mask_v(i,j,k1)

      enddo
      enddo

   14 continue

! Reste à appliquer la multiplication par 0.25 et anticiper
! sur le fait que dans le mode externe ce terme sera multiplié
! par (H+sse)*DY_u (respectivement (H+sse)*DX_v)
! ATTENTION SI ON REVIENT A L'ANCIEN SCHEMA IL FAUT TENIR QUAND MEME COMPTE DE CE QUE
! FROZEN EST DIFFERENT PAR rapport A la prise en compte des facteurs d'echelle dans la nouvelle
! version
      do j=2,jmax-1
      do i=2,imax
       adve3d2d_u(i,j)=0.25*adve3d2d_u(i,j)*invdxdy_u(i,j) !06-09-18
!      adve3d2d_u(i,j)=0.25*adve3d2d_u(i,j)/(hz_u(i,j,1)*dy_u(i,j)) !28-11-14
      enddo
      enddo
      do j=2,jmax
      do i=2,imax-1
       adve3d2d_v(i,j)=0.25*adve3d2d_v(i,j)*invdxdy_v(i,j) !06-09-18
!      adve3d2d_v(i,j)=0.25*adve3d2d_v(i,j)/(hz_v(i,j,1)*dx_v(i,j)) !28-11-14
      enddo
      enddo


      return
      endif                                  !m°v°m>

      end subroutine adve3dto2d
