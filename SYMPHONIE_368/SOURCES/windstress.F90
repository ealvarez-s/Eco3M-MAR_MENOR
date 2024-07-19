      subroutine windstress
!______________________________________________________________________
! SYMPHONIE ocean model
! release 290 - last update: 23-10-20
!______________________________________________________________________

      use module_principal
      implicit none
#ifdef synopsis
       subroutinetitle='windstress'
       subroutinedescription= &
       'Computes a "simple" wind stress when realistic meteolorolical' &
       //' fields are not available'
       call main_synopsis(subroutinetitle,subroutinedescription)
#endif

!..............................................................................
! modifs 01/11/01: inversion de l'ordre des boucles: J passe avant I
!        21/03/02: par defaut le vent est nul. Ca evite les erreurs
!                  d'etourderie comme par ex lancer model_streamf en
!                  pensant (à tort) que le vent est nul.
!        08-10-14  possibilite de construire des flux a partir des formules bulk
!        25-01-15  Attention par defaut wstress_ =0 pour prévoir le cas des
!                  simulation forcee par les vagues seule (sans modele meteo)
!                  car l'effet taw-two est un delta qui s'ajoute a wstress_...
!        05-08-16  ajout d'un commentaire
!        23-11-15  3 arguments passent dans bulk_formulae
!        24-08-16  deplacement de l'etiquette bidon pour permettre calcul wstress_w
!..............................................................................

!...........................................................................
!        25-01-15  Attention par defaut wstress_ =0 pour prévoir le cas des
!                  simulation forcee par les vagues seule (sans modele meteo)
!                  car l'effet taw-two est un delta qui s'ajoute a wstress_...
!                  wstress=rho*(ustar**2)
!                  wstress_u(:,:,1)=0. ; wstress_v(:,:,1)=0. ! reset necessaire a module_wave
! v290   23-10-20  La contrainte enoncee ci-dessus disparait avec la multiplication par 
!                  iairsea, desormais, de wstress_u dans module_wave, ce
!                  qui evite ce reset toutes les iterations.
!...............................................................................
!    _________                    .__                  .__             ! m°v°m
!   /   _____/__.__. _____ ______ |  |__   ____   ____ |__| ____       !
!   \_____  <   |  |/     \\____ \|  |  \ /  _ \ /    \|  |/ __ \      ! 
!   /        \___  |  Y Y  \  |_> >   Y  (  <_> )   |  \  \  ___/      !
!  /_______  / ____|__|_|  /   __/|___|  /\____/|___|  /__|\___  >     !
!          \/\/          \/|__|        \/            \/        \/      !
!......................................................................! m[°o°]m

      RETURN !23-10-20

#ifdef bidon
      if(iteration3d==0) then !0000000000>

! Reset, declaration au demarrage....
      if(.not.allocated(ssr_w)) then
              allocate (ssr_w     (0:imax+1,0:jmax+1,0:2) ) ; ssr_w=0
!     else
!      write(6,*)lbound(ssr_w) 
!      write(6,*)ubound(ssr_w) 
!      stop 'coucou'
      endif

      if(allocated(snsf_w))deallocate(snsf_w)
      allocate (snsf_w(0:imax+1,0:jmax+1,0:2) ) ; snsf_w=0

      if(allocated(precipi_w))deallocate(precipi_w)
              allocate (precipi_w (0:imax+1,0:jmax+1,0:2) ) ; precipi_w=0

      if(allocated(pss_w))deallocate(pss_w)
              allocate (pss_w     (0:imax+1,0:jmax+1,0:2) ) ; pss_w=0.

      if(allocated(uwind_t))deallocate(uwind_t)
              allocate (uwind_t   (0:imax+1,0:jmax+1,0:2) ) ; uwind_t=0

      if(allocated(vwind_t))deallocate(vwind_t)
              allocate (vwind_t   (0:imax+1,0:jmax+1,0:2) ) ; vwind_t=0

      if(allocated(q2_t))deallocate(q2_t)
            allocate (q2_t      (0:imax+1,0:jmax+1,0:2) ) ; q2_t=0

      if(allocated(q10_t))deallocate(q10_t)
            allocate (q10_t      (0:imax+1,0:jmax+1) ) ; q10_t=0

      if(allocated(teta2_t))deallocate(teta2_t)
              allocate(teta2_t   (0:imax+1,0:jmax+1,0:2) ) ; teta2_t=0

      if(allocated(teta10_t))deallocate(teta10_t)
              allocate(teta10_t   (0:imax+1,0:jmax+1) ) ; teta10_t=0

      if(allocated(pss_mean))deallocate(pss_mean)
              allocate(pss_mean  (0:2) )                   ; pss_mean=0

      if(allocated(tetastar_t))deallocate(tetastar_t)
      allocate(tetastar_t(0:imax+1,0:jmax+1,0:2)) ; tetastar_t=0

      if(allocated(ustar_t))deallocate(ustar_t)
      allocate(ustar_t(0:imax+1,0:jmax+1,0:2)) ; ustar_t=0

      if(allocated(qstar_t))deallocate(qstar_t)
      allocate(qstar_t(0:imax+1,0:jmax+1,0:2)) ; qstar_t=0

! Valeurs academiques
      do j=0,jmax+1
      do i=0,imax+1
       uwind_t(i,j,0:2)=10.                      !m/s
       vwind_t(i,j,0:2)=0.                       !m/s
       pss_w  (i,j,0:2)=101300.                  ! Pascals
       teta2_t(i,j,0:2)=273.15+20.               !  Kelvins
       q2_t   (i,j,0:2)=0.001                    ! specific humidity
      enddo
      enddo

      call bulk_formulae(1,0,1,relativewind) !23-11-15

       i=imax/2
       j=jmax/2
       write(6,*)'wstress_u=            ',wstress_u(i,j,1)
       write(6,*)'wstress_v=            ',wstress_v(i,j,1)
       write(6,*)'Flux chaleur latente= ',slhf_w(i,j,1)      
       write(6,*)'Flux chaleur sensible=',sshf_w(i,j,1)

!      stop 'coco'

      else                  !0000000000>
    
! Phase ietrative
       call bulk_formulae(2,0,1,relativewind) !23-11-15

      endif                 !0000000000>

#endif

#ifdef bidon
!___________________________________________________
! 0.1 est un ordre de grandeur classique pour une tension de vent
! (en gros correspond à un vent de 10m/s)
! Une parametriation simple de la tension de vent
! en fonction de la vitesse du vent à 10 mètres:
! module du vent:
!     X3=10.*AMIN1(1.,REAL(KOUNT-KOUNT0)*DTI_FW/(17.*3600.))
      x3=0.                                                            !21/03/02
!     X3=10.
! direction du vent (convention météo: Du Nord vers le Sud   = 0°)
! direction du vent (convention météo: De l'Est vers l'Ouest =90°)
      x2=150.*pi/180.
!     X2=360.*PI/180.*REAL(KOUNT-KOUNT0)/(24.*3600.)

! drag coef surface
      x1=(0.8+0.065*x3)*1.0e-3
      do j=0,jmax+1
      do i=0,imax+1
      xy_t(i,j,1)=-x1*rhoair*x3*x3*sin(x2)    ! Ouest vers Est
      xy_t(i,j,2)=-x1*rhoair*x3*x3*cos(x2)    ! Sud vers Nord
      enddo
      enddo
              do j=1,jmax+1
              do i=1,imax+1

               wstress_u(i,j,1)=                                        &
                     0.5*cos(angle0)*(xy_t(i,j,1)+xy_t(i-1,j,1))        &
                    -0.5*sin(angle0)*(xy_t(i,j,2)+xy_t(i-1,j,2))

               wstress_v(i,j,1)=                                        &
                     0.5*sin(angle0)*(xy_t(i,j,1)+xy_t(i,j-1,1))        &
                    +0.5*cos(angle0)*(xy_t(i,j,2)+xy_t(i,j-1,2))
              enddo
              enddo
!             WRITE(6,*)'coucou',WSTRESS_X(1,1,1),WSTRESS_Y(1,1,1)
!___________________________________________________


!..............................................................
! CAS ACADEMIQUE D'UN CYCLONE
! LE VENT DERIVE D'UNE FONCTION DE COURANT:
!     X1=REAL(MECO)/2.
!     X2=REAL(NECO)*REAL(KOUNT)/100.
!     DO I=1,MECO+1
!     DO J=1,NECO+1
!     DIST=SQRT( (REAL(I)-X1)**2 + (REAL(J)-X2)**2 )
!     DIST0=20.
!     XY_R(I,J,1)=1.E4*AMAX1(EXP(- (DIST/DIST0)**2 )-0.01,0.)
!     ENDDO
!     ENDDO

!     DO I=1,MECO+1
!     DO J=1,NECO
!     WSTRESS_X(I,J)=(XY_R(I,J,1)-XY_R(I,J+1,1))/DYB
!     ENDDO
!     ENDDO

!     DO I=1,MECO
!     DO J=1,NECO+1
!     WSTRESS_Y(I,J)=(XY_R(I+1,J,1)-XY_R(I,J,1))/DXB
!     ENDDO
!     ENDDO
!..............................................................
#endif

#ifdef bidon
      x1=0.2 ! OE wstress
      x2=0.  ! SN wstress
      do j=1,jmax   ; do i=1,imax+1
       wstress_u(i,j,1)=0.5*(gridrotcos_t(i,j)+gridrotcos_t(i-1,j))*x1       &
                       -0.5*(gridrotsin_t(i,j)+gridrotsin_t(i-1,j))*x2
      enddo ; enddo
      do j=1,jmax+1 ; do i=1,imax
       wstress_v(i,j,1)=0.5*(gridrotsin_t(i,j)+gridrotsin_t(i,j-1))*x1       &
                       +0.5*(gridrotcos_t(i,j)+gridrotcos_t(i,j-1))*x2
      enddo ; enddo
#endif

! Module de la tension de vent:

      do j=1,jmax
      do i=1,imax

      wstress_w(i,j)=sqrt(                                              &
       ((wstress_u(i,j,1)+wstress_u(i+1,j,1))/2.)**2+                   &
       ((wstress_v(i,j,1)+wstress_v(i,j+1,1))/2.)**2)

      enddo
      enddo

!     CALL GRAPH_OUT
!     WRITE(6,*)'DANS SBR WINSTRESS'
!     WRITE(6,*)'-----------------------'

      end subroutine windstress
